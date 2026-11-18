#!/usr/bin/env python3
"""Prepare identical paired inputs and run frozen T1K baselines, without truth."""
import argparse
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

QUALITY_BYTES = bytes(range(33, 127))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def records(stream, mate):
    while True:
        header = stream.readline()
        if not header:
            return
        sequence, plus, quality = (stream.readline().rstrip(b'\r\n') for _ in range(3))
        if not header.startswith(b'@') or not plus.startswith(b'+') or not sequence or len(sequence) != len(quality):
            raise ValueError('Malformed four-line FASTQ')
        if sequence.translate(None, b'ACGTRYSWKMBDHVNacgtryswkmbdhvn') or quality.translate(None, QUALITY_BYTES):
            raise ValueError('Invalid sequence/quality alphabet')
        fields = header[1:].split()
        if not fields:
            raise ValueError('FASTQ read name missing')
        name = fields[0]
        if name.endswith((b'/1', b'/2')):
            if int(name[-1:]) != mate:
                raise ValueError('FASTQ mate suffix disagrees with input slot')
            name = name[:-2]
        if len(fields) > 1 and re.match(rb'[12]:[YN]:', fields[1]):
            if int(fields[1][:1]) != mate:
                raise ValueError('CASAVA mate label disagrees with input slot')
        yield name, len(sequence)


def verify_pairs(r1, r2):
    count, minimum, maximum = 0, 10**9, 0
    with gzip.open(r1, 'rb') as first, gzip.open(r2, 'rb') as second:
        for a, b in itertools.zip_longest(records(first, 1), records(second, 2)):
            if a is None or b is None or a[0] != b[0]:
                raise ValueError('FASTQ mate names or counts differ')
            count += 1
            minimum = min(minimum, a[1], b[1])
            maximum = max(maximum, a[1], b[1])
    if not count:
        raise ValueError('No read pairs')
    return {'pairs': count, 'minimum_read_length': minimum, 'maximum_read_length': maximum}


def prepare_lai(manifest, index, output):
    samples = json.loads(manifest.read_text())['samples']
    if not 0 <= index < len(samples):
        raise IndexError('LAI sample index outside manifest')
    row = samples[index]
    if row['sample_key'] != 'lai.' + row['donor'] or not re.fullmatch(r'(HG|NA)[0-9]+', row['donor']):
        raise ValueError('Invalid LAI sample identity')
    record = {'sample_key': row['sample_key'], 'donor': row['donor'], 'assay': 'targeted_capture',
              'status': 'preparing', 'source_manifest_sha256': sha(manifest),
              'driver_sha256': sha(__file__), 'started': time.time()}
    destination = output / row['sample_key'] / 'PAIRS.json'
    for mate in (1, 2):
        source = row['reads'][str(mate)]
        path = Path(source['path'])
        if path.stat().st_size != source['bytes'] or sha(path) != source['sha256']:
            raise ValueError('Read source differs from verified download receipt')
        record['r' + str(mate)] = str(path)
        record['r' + str(mate) + '_sha256'] = source['sha256']
    record.update(verify_pairs(record['r1'], record['r2']))
    record.update(status='paired_fastq_verified', finished=time.time())
    save(destination, record)
    return destination


def validate_runtime(config):
    artifacts = {item['path']: item['sha256'] for item in config['artifacts']}
    if len(artifacts) != len(config['artifacts']):
        raise ValueError('Duplicate pinned baseline artifact')
    required = {config[key] for key in ('run_t1k', 'original_reference', 'genomic_reference', 'aliases')}
    required.update(str(Path(config['decoder_dir']) / (name + '.py')) for name in ('decode_t1k', 'anchor_t1k_coarse'))
    if not required <= artifacts.keys():
        raise ValueError('Baseline runtime uses an unpinned artifact')
    for item in config['artifacts']:
        if sha(item['path']) != item['sha256']:
            raise ValueError('Pinned baseline artifact changed: ' + Path(item['path']).name)


def extend_labels(calls, aliases):
    """Keep original full labels/suffixes alongside legacy numeric projections."""
    for gene in calls.values():
        for resolution in gene['resolutions'].values():
            for slot in resolution:
                full = set()
                for token in slot.get('raw', '').split(','):
                    candidate = token if token in aliases else token.removeprefix('HLA-')
                    if candidate in aliases:
                        full.update(aliases[candidate].get('genomic_labels', []))
                        full.update(aliases[candidate].get('cds_labels', []))
                    elif '*' in candidate:
                        full.add(candidate)
                slot['source_full_labels_including_suffixes'] = sorted(full)
    return calls


def run_t1k(pairs_file, runtime_file, output, threads):
    pairs = json.loads(pairs_file.read_text())
    runtime = json.loads(runtime_file.read_text())
    if pairs['status'] != 'paired_fastq_verified' or not re.fullmatch(r'(lai|wgs)\.(HG|NA)[0-9]+', pairs['sample_key']):
        raise ValueError('Unverified input or unsafe sample key')
    if pairs['sample_key'] != ('wgs.' if pairs['assay'] == 'wgs' else 'lai.') + pairs['donor'] or pairs['assay'] not in ('wgs', 'targeted_capture'):
        raise ValueError('Paired input identity or assay mismatch')
    validate_runtime(runtime)
    for mate in (1, 2):
        if sha(pairs['r' + str(mate)]) != pairs['r' + str(mate) + '_sha256']:
            raise ValueError('Prepared read hash changed')
    output.mkdir(parents=True, exist_ok=True)
    completion = output / 'COMPLETE.json'
    if completion.exists():
        old = json.loads(completion.read_text())
        expected = {'t1k_original/PREDICTIONS.json', 't1k_original/calls_genotype.tsv',
                    't1k_genomic_ipd/PREDICTIONS.json', 't1k_genomic_ipd/calls_genotype.tsv',
                    'T1K_GENOMIC_COARSE_ANCHORED.json'}
        if (old.get('status') == 'complete' and old.get('sample_key') == pairs['sample_key']
                and old.get('pairs_sha256') == sha(pairs_file)
                and old.get('runtime_sha256') == sha(runtime_file)
                and old.get('driver_sha256') == sha(__file__)
                and set(old.get('output_sha256', {})) == expected):
            if all(sha(output / name) == checksum for name, checksum in old['output_sha256'].items()):
                return old
        raise ValueError('Existing completion differs from current pinned inputs')
    sys.path.insert(0, runtime['decoder_dir'])
    from decode_t1k import decode
    from anchor_t1k_coarse import anchor
    record = {'status': 'running', 'sample_key': pairs['sample_key'], 'donor': pairs['donor'],
              'assay': pairs['assay'], 'pairs_sha256': sha(pairs_file), 'runtime_sha256': sha(runtime_file),
              'driver_sha256': sha(__file__), 'started': time.time(), 'commands': [],
              'scope': 'Baseline inference only; published truth not loaded, no accuracy scoring.',
              'database_version': 'IPD-IMGT/HLA 3.65.0', 'truth_database_version': '3.59.0 (not loaded)',
              'nomenclature_reconciliation_required_before_scoring': True}
    save(output / 'RUN.json', record)
    all_calls = {}
    try:
        for method, reference, alias_path in [('t1k_original', runtime['original_reference'], None),
                                               ('t1k_genomic_ipd', runtime['genomic_reference'], runtime['aliases'])]:
            directory = output / method
            directory.mkdir(exist_ok=True)
            prefix = directory / 'calls'
            preset = 'hla-wgs' if pairs['assay'] == 'wgs' else 'hla'
            command = [runtime['run_t1k'], '-1', pairs['r1'], '-2', pairs['r2'], '--preset', preset,
                       '-f', reference, '-t', str(threads), '--alleleDigitUnits', '4', '--alleleDelimiter', ':',
                       '-o', str(prefix)]
            record['commands'].append(command)
            save(output / 'RUN.json', record)
            with (directory / 'run.log').open('w') as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            table = directory / 'calls_genotype.tsv'
            if not table.is_file() or not table.stat().st_size:
                raise ValueError('Missing T1K genotype table')
            aliases = json.loads(Path(alias_path).read_text()) if alias_path else {}
            calls = extend_labels(decode(table.read_text(), aliases), aliases)
            save(directory / 'PREDICTIONS.json', {'method': method, 'sample_key': pairs['sample_key'],
                'calls': calls, 'raw_table_sha256': sha(table), 'reference_sha256': sha(reference),
                'note': 'Numeric projections are legacy decoder output; full source labels/suffixes are retained.'})
            all_calls[method] = calls
        anchored, decisions = anchor(all_calls['t1k_original'], all_calls['t1k_genomic_ipd'])
        save(output / 'T1K_GENOMIC_COARSE_ANCHORED.json', {'method': 't1k-genomic-ipd-coarse-anchor',
             'sample_key': pairs['sample_key'], 'calls': anchored, 'decisions': decisions,
             'scope': 'Truth-free composition of frozen original and genomic-IPD predictions.'})
        files = [p for p in output.rglob('*') if p.is_file() and p.name in ('PREDICTIONS.json', 'calls_genotype.tsv', 'T1K_GENOMIC_COARSE_ANCHORED.json')]
        record.update(status='complete', finished=time.time(),
                      output_sha256={str(p.relative_to(output)): sha(p) for p in files})
        save(completion, record)
        save(output / 'RUN.json', record)
        return record
    except BaseException as exc:
        record.update(status='failed', error_type=type(exc).__name__, finished=time.time())
        save(output / 'RUN.json', record)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare-lai')
    prep.add_argument('--manifest', type=Path, required=True)
    prep.add_argument('--index', type=int, required=True)
    prep.add_argument('--output', type=Path, required=True)
    run = sub.add_parser('run-t1k')
    run.add_argument('--pairs', type=Path, required=True)
    run.add_argument('--runtime', type=Path, required=True)
    run.add_argument('--output', type=Path, required=True)
    run.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if args.action == 'prepare-lai':
        result = prepare_lai(args.manifest, args.index, args.output)
        print(json.dumps({'status': 'paired_fastq_verified', 'manifest': str(result)}))
    else:
        result = run_t1k(args.pairs, args.runtime, args.output, args.threads)
        print(json.dumps({'status': result['status'], 'sample_key': result['sample_key']}))
