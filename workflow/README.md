# DogoHLA Portable Typing Workflow

A modular, containerized CWL v1.2 workflow for high-accuracy HLA typing, custom pangenome panel construction, and mosaic candidate discovery.

## 1. Overview

The DogoHLA workflow packages the complete HLA typing workflow developed for diverse population cohorts:
- **Part 1 (`extract-genes`):** Extracts full genomic and coding sequences for 8 key HLA genes (`HLA-A`, `HLA-B`, `HLA-C`, `HLA-DPA1`, `HLA-DPB1`, `HLA-DQA1`, `HLA-DQB1`, `HLA-DRB1`) from phased diploid assemblies using Immuannot v3 annotation models and 2,000 bp flanking sequence rules.
- **Part 2 (`build-panel`):** Augments the baseline IPD-IMGT/HLA reference with population-specific novel HLA candidate alleles (deduplicated, non-overlapping with IPD, aliased as `HLA-<gene>*PG<sha12>`).
- **Part 3 (`repull-reads`):** Recruits MHC short reads from a CRAM/BAM: the regional recruitment (chr6:28.5-33.5 Mb plus all `chr6_*` and `HLA-*` contigs, with mate rescue) plus every whole-genome read pair that carries at least 3 bait 31-mers.
- **Part 4 (`type-t1k`):** Performs linear-reference HLA typing using T1K across baseline (Arm A) and augmented (Arm B) reference panels.
- **Part 5 (`patch-builder`):** Discovers unannotated and recombinant alleles from the T1K candidate reads, applying the frozen R4n abstain rule to prevent over-prediction (arms C4R4n and D4R4n).
- **Part 6 (`score`):** Evaluates typing accuracy against the truth sequences of the cohort, per copy and arm (A, B, C4, C4R4n, D4, D4R4n).

Parts 3 to 6 run the frozen reference code vendored in `dogohla/frozen/`. The scripts in `scripts/` only call it. Every vendored file is checked against `dogohla/frozen/MANIFEST.sha256` at each start; a changed file stops the run.

---

## 2. Requirements & Dependencies

### Hardware
- **CPU:** 8+ cores recommended (16 cores optimal for whole-genome read extraction and alignment).
- **Memory:** 64 GB for `repull-reads` from a 30x CRAM (our test runs used all 64 GB they were given); 16 GB is enough for the other parts.
- **Storage:** ~50 GB scratch space per sample during CRAM/FASTQ processing.

### Software
- **Container Engine:** Singularity / Apptainer (v3.8+) or Docker / Podman.
- **Workflow Runner:** `cwltool` (v3.1+) or `toil-cwl-runner`.
- **Container Image:** `ghcr.io/asian-hla-consortium/dogohla-workflow:latest` (or built via `containers/Apptainer.def` / `containers/Dockerfile`).

### Input files you need
| Input | What it is |
|---|---|
| CRAM or BAM of the sample, **with its index** | The index must sit next to the file: `sample.cram.crai` (create with `samtools index sample.cram`), or `sample.bam.bai`. The workflow stops at the start if it is missing. |
| GRCh38 reference FASTA **with its `.fai`** | `GRCh38_full_analysis_set_plus_decoy_hla.fa` and `GRCh38_full_analysis_set_plus_decoy_hla.fa.fai` (`samtools faidx`). `scripts/fetch_resources.sh` makes both. |
| `genomic.fa`, `bait_kmers.txt` | IPD-IMGT/HLA 3.65.0 genomic alleles in T1K format and the 31-mer bait set; `bait_kmers.txt` is built by `scripts/fetch_resources.sh`; `genomic.fa` is supplied by the consortium. |
| `confounded_v1/` | Confounded windows, shipped in `resources/confounded_v1`. |
| Cohort resource folder (`builder_inputs`) | The cohort files the builder and the scorer read (`SOURCE_RECORDS.json`, `FOLDS.json`, `mosaic_classes.tsv`, `panels_atlas.py`, `INPUTS.sha256.tsv`, `msa/`). Supplied by the consortium; it is not downloaded by `fetch_resources.sh`. The donor must be one of the cohort donors in `FOLDS.json`. |
| Panel `reference.fa` and `aliases.json` (optional) | The arm-B reference (IPD plus cohort panel) and its alias file. Without them only arm A and C4R4n run. |

---

## 3. Quick Start (Single Sample Typing)

The commands below are run from this `workflow/` folder.

### Step 0: Get the container image
`cwltool --singularity` looks for the image as a file named after the `dockerPull` name in the folder `CWL_SINGULARITY_CACHE`. If it is not there, cwltool pulls it itself, with default temporary folders and without a progress message (about 10 minutes). Pull it once yourself, with the temporary and cache folders on an ordinary filesystem (your home directory), not on an encrypted or ACL-managed one:
<!-- cmd: pull -->
```bash
export SINGULARITY_TMPDIR=$HOME/sif-tmp SINGULARITY_CACHEDIR=$HOME/sif-cache CWL_SINGULARITY_CACHE=$HOME/sif-images
mkdir -p $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $CWL_SINGULARITY_CACHE
singularity pull --name $CWL_SINGULARITY_CACHE/ghcr.io_asian-hla-consortium_dogohla-workflow:latest.sif docker://ghcr.io/asian-hla-consortium/dogohla-workflow:latest
```
Keep `CWL_SINGULARITY_CACHE` set in the shell where you run `cwltool`.

### Step 1: Fetch Reference Resources
```bash
singularity exec $CWL_SINGULARITY_CACHE/ghcr.io_asian-hla-consortium_dogohla-workflow:latest.sif bash scripts/fetch_resources.sh resources
```
This downloads or builds, and checks the sha256 of each result:
- `resources/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa` and its `.fai` (GRCh38 with decoys and HLA contigs; 1000 Genomes).
- `resources/ipd/genomic.fa`: IPD-IMGT/HLA v3.65.0 genomic reference in T1K format. This file is checked, not built: place the pinned `genomic.fa` (sha256 `9202c3d3...`, supplied by the consortium) in `resources/ipd/` before running the script. The public IPD data alone gives the same alleles with shorter sequences, so it does not reproduce the pinned file.
- `resources/bait/bait_kmers.txt`: precomputed 31-mer bait set (888,049 k-mers).
- `resources/immuannot/Data-ipd3590-kir2130`: Immuannot reference data.

### Step 2: Prepare the input files
`type_inputs.yml` (the CRAM index `.crai` and the reference `.fai` are found next to their files):
<!-- file: type_inputs.yml -->
```yaml
donor: "SAMPLE01"
cram:
  class: File
  path: "/path/to/sample.cram"
reference_genome:
  class: File
  path: "/path/to/resources/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
ipd_genomic:
  class: File
  path: "/path/to/resources/ipd/genomic.fa"
bait_kmers:
  class: File
  path: "/path/to/resources/bait/bait_kmers.txt"
panel_reference:
  class: File
  path: "/path/to/panel/reference.fa"
panel_aliases:
  class: File
  path: "/path/to/panel/aliases.json"
builder_inputs:
  class: Directory
  path: "/path/to/cohort_resource"
confounded_dir:
  class: Directory
  path: "resources/confounded_v1"
threads: 8
```
Leave out `panel_reference` and `panel_aliases` to run arm A and C4R4n only.

`evaluate_inputs.yml` (scores the outputs of the typing run; the file names are the ones the typing run writes to `out/type/`):
<!-- file: evaluate_inputs.yml -->
```yaml
donor: "SAMPLE01"
calls_a:
  class: File
  path: "out/type/armA_genotype.tsv"
calls_b:
  class: File
  path: "out/type/armB_genotype.tsv"
builder_dir_c:
  class: Directory
  path: "out/type/builder_C4R4n"
builder_dir_d:
  class: Directory
  path: "out/type/builder_D4R4n"
panel_reference:
  class: File
  path: "/path/to/panel/reference.fa"
ipd_genomic:
  class: File
  path: "/path/to/resources/ipd/genomic.fa"
builder_inputs:
  class: Directory
  path: "/path/to/cohort_resource"
```
Without a panel leave out `calls_b`, `builder_dir_d` and `panel_reference`.

### Step 3: Run the workflow
cwltool writes its temporary files to `$TMPDIR`; point it to a folder with about 100 GB free.
<!-- cmd: run -->
```bash
cwltool --singularity --outdir out/type type.cwl type_inputs.yml
cwltool --singularity --outdir out/evaluate evaluate.cwl evaluate_inputs.yml
```
The first command takes about 5.5 hours (see below), the second about 10 minutes.

---

## 4. Running on a Slurm Cluster

For distributed execution across multiple samples, use `toil-cwl-runner`:
```bash
toil-cwl-runner \
  --batchSystem slurm \
  --singularity \
  --clean always \
  --logFile workflow.log \
  type.cwl type_inputs.yml
```

### Estimated Runtime
Measured on our test runs (8 CPUs):
- **`extract-genes` with annotation:** about 4 minutes per assembly.
- **`repull-reads`:** about 5 hours per 30x whole-genome CRAM (most of it reading the whole file).
- **`type-t1k` + `patch-builder`:** about 40 minutes per sample (both reference arms); **`score`:** about 10 minutes.
- **Total per sample from CRAM:** about 6 hours, almost all of it read fetching.

---

## 5. Workflow Outputs

`type.cwl` writes, in `out/type/`:
- `<donor>_r1.fq.gz`, `<donor>_r2.fq.gz`, `<donor>_PAIRS.json`: the repulled read set and its counts.
- `armA_genotype.tsv`: T1K calls on the IPD reference (arm A). `armB_genotype.tsv`: T1K calls on the panel reference (arm B), only with a panel.
- `<donor>_C4R4n_final_calls.tsv`, `<donor>_D4R4n_final_calls.tsv`: final calls after the builder and the R4n abstain rule (`D4R4n` only with a panel); columns `gene`, `slot1`, `slot2`.
- `builder_C4R4n/`, `builder_D4R4n/`: the builder folders, per gene `<gene>/proposals.json`, `<gene>/additions.fa` (novel allele sequences), `<gene>/OVERRIDES.json` (calls of arm C4/D4), `<gene>/OVERRIDES_R4n.json` (calls of arm C4R4n/D4R4n), further abstain variants, and `COMPLETE.json`.

`evaluate.cwl` writes, in `out/evaluate/`:
- `score_rows.tsv`: one row per truth copy and arm (A, B, C4, C4R4n, D4, D4R4n), all columns of the frozen scorer.
- `score_avail.tsv`: which arm tables were found.
