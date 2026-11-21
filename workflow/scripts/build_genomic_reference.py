#!/usr/bin/env python3
"""Build the pinned IPD genomic reference (genomic.fa) from the public IPD-IMGT/HLA hla.dat.

Method (the one that made our pinned file, sha256 9202c3d3...):
  1. baseline = T1K t1k-build.pl on hla.dat, DNA mode: truncated entries for every allele
     ("<prefix>_dna_seq.fa");
  2. genomic  = T1K ParseDatFile.pl hla.dat --mode genome --gene HLA: full genomic entries for the alleles
     that have a full genomic sequence in IPD;
  3. every baseline entry keeps its name and position (the order of the pinned file is given by --order, see below); where a genomic entry of the same name exists
     (and its coding coordinates are valid) it replaces the baseline entry; all other baseline
     entries stay unchanged, including entries with inherited coordinate defects.
This file is a plain re-implementation of that combination step (Python only).

Order of the records: T1K's ParseDatFile.pl appends the "rescued partial" alleles by iterating a Perl hash,
so the order of its DNA-mode output differs from run to run (names, headers and sequences are identical).
The pinned file has one fixed order; it is stored as the allele name list resources/ipd/genomic_order.txt.gz
(IPD allele names only) and given with --order, which makes the output byte-identical to the pinned file.

  python3 build_genomic_reference.py --baseline hla_dna_seq.fa --genomic genomic_mode.fa --out genomic.fa
"""
import argparse
import json
from pathlib import Path


def records(text, validate_coordinates=True):
    result = {}
    name, header, parts = None, None, []

    def save():
        if name is None:
            return
        if name in result:
            raise ValueError('Duplicate allele: ' + name)
        sequence = ''.join(parts).upper()
        fields = header.split()
        count = int(fields[1])
        coordinates = list(map(int, fields[2:]))
        if not sequence or count < 1 or len(coordinates) != 2 * count:
            raise ValueError('Incomplete reference annotation: ' + name)
        previous = -1
        for start, end in zip(coordinates[::2], coordinates[1::2]):
            if validate_coordinates and not previous < start <= end < len(sequence):
                raise ValueError('Invalid reference coordinates: ' + name)
            previous = end
        result[name] = (header, sequence)
    for line in text.splitlines():
        if line.startswith('>'):
            save()
            header = line[1:]
            name, parts = header.split()[0], []
        else:
            parts.append(line.strip())
    save()
    if not result:
        raise ValueError('Empty reference')
    return result


def invalid_annotations(data):
    invalid = {}
    for name, (header, sequence) in data.items():
        try:
            records('>' + header + '\n' + sequence + '\n')
        except ValueError as error:
            invalid[name] = str(error)
    return invalid


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--baseline', type=Path, required=True)
    ap.add_argument('--genomic', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--order', type=Path, help='text file (optionally .gz) with one allele name per line: output order')
    ap.add_argument('--summary', type=Path)
    a = ap.parse_args()
    original = records(a.baseline.read_text(), validate_coordinates=False)
    genomic = records(a.genomic.read_text(), validate_coordinates=False)
    genomic_invalid = invalid_annotations(genomic)
    genomic = {n: r for n, r in genomic.items() if n not in genomic_invalid}
    combined = dict(original)
    replaced = []
    for name, record in genomic.items():
        if name in combined:
            combined[name] = record
            replaced.append(name)
    if not replaced or set(combined) != set(original):
        raise SystemExit('candidate-preserving genome replacement failed')
    if a.order:
        import gzip
        op = gzip.open if str(a.order).endswith('.gz') else open
        with op(a.order, 'rt') as fo:
            names = [x.strip() for x in fo if x.strip()]
        if sorted(names) != sorted(combined) or len(set(names)) != len(names):
            raise SystemExit('--order is not a permutation of the allele names')
        combined = {n: combined[n] for n in names}
    with a.out.open('w') as f:
        for header, sequence in combined.values():
            f.write('>' + header + '\n' + sequence + '\n')
    s = dict(candidates=len(original), genomic_replacements=len(replaced),
             retained_original_candidates=len(original) - len(replaced),
             unused_genomic_candidates=len(set(genomic) - set(original)), rejected_genomic=len(genomic_invalid))
    if a.summary:
        a.summary.write_text(json.dumps(s, indent=2) + '\n')
    print(json.dumps(s))


if __name__ == '__main__':
    main()
