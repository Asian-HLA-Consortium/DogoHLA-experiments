# DogoHLA Portable Typing Workflow

A modular, containerized CWL v1.2 workflow for high-accuracy HLA typing, custom pangenome panel construction, and mosaic candidate discovery.

## 1. Overview

The DogoHLA workflow packages the complete HLA typing workflow developed for diverse population cohorts:
- **Part 1 (`extract-genes`):** Extracts full genomic and coding sequences for 8 key HLA genes (`HLA-A`, `HLA-B`, `HLA-C`, `HLA-DPA1`, `HLA-DPB1`, `HLA-DQA1`, `HLA-DQB1`, `HLA-DRB1`) from phased diploid assemblies using Immuannot v3 annotation models and 2,000 bp flanking sequence rules.
- **Part 2 (`build-panel`):** Augments the baseline IPD-IMGT/HLA reference with population-specific novel HLA candidate alleles (deduplicated, non-overlapping with IPD, aliased as `HLA-<gene>*PG<sha12>`).
- **Part 3 (`repull-reads`):** Recruits MHC short reads from CRAM/BAM or WGS FASTQ using rule-2 regional pulling plus whole-genome candidate pair extraction with 31-mer bait filtering.
- **Part 4 (`type-t1k`):** Performs linear-reference HLA typing using T1K across baseline (Arm A) and augmented (Arm B) reference panels.
- **Part 5 (`patch-builder`):** Discovers unannotated and recombinant alleles using read pileups against pass-1 alleles, applying the frozen R4n abstain rule to prevent over-prediction.
- **Part 6 (`score`):** Evaluates typing accuracy against assembly-derived truth sequences across 2-field, 3-field, and 4-field resolutions.

All modules use frozen reference algorithms with SHA-256 integrity verification (`dogohla/frozen/MANIFEST.sha256`) checked at execution start.

---

## 2. Requirements & Dependencies

### Hardware
- **CPU:** 8+ cores recommended (16 cores optimal for whole-genome read extraction and alignment).
- **Memory:** 64 GB for `repull-reads` from a 30x CRAM (our test runs used all 64 GB they were given); 16 GB is enough for the other parts.
- **Storage:** ~50 GB scratch space per sample during CRAM/FASTQ processing.

### Software
- **Container Engine:** Singularity / Apptainer (v3.8+) or Docker / Podman.
- **Workflow Runner:** `cwltool` (v3.1+) or `toil-cwl-runner` for Slurm clusters.
- **Container Image:** `ghcr.io/asian-hla-consortium/dogohla-workflow:latest` (or built via `containers/Apptainer.def` / `containers/Dockerfile`).

---

## 3. Quick Start (Single Sample Typing)

### Step 1: Fetch Reference Resources
```bash
bash scripts/fetch_resources.sh
```
This downloads and verifies:
- Pinned IPD-IMGT/HLA v3.65.0 genomic reference (`genomic.fa`).
- Precomputed 31-mer bait dictionary (`bait_kmers.txt`, 888,049 k-mers).
- Immuannot reference models (`Data-ipd3590-kir2130`).

### Step 2: Prepare Input YAML
Create `inputs.yml`:
```yaml
sample_id: "SAMPLE01"
cram_or_bam:
  class: File
  path: "/path/to/sample.cram"
reference_genome:
  class: File
  path: "/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa"
ipd_genomic:
  class: File
  path: "/path/to/resources/ipd/genomic.fa"
bait_kmers:
  class: File
  path: "/path/to/resources/bait/bait_kmers.txt"
confounded_dir:
  class: Directory
  path: "resources/confounded_v1"
threads: 8
```

### Step 3: Run Workflow with CWLTool
```bash
cwltool --singularity type.cwl inputs.yml
```

---

## 4. Running on a Slurm Cluster

For distributed execution across multiple samples, use `toil-cwl-runner`:
```bash
toil-cwl-runner \
  --batchSystem slurm \
  --singularity \
  --clean always \
  --logFile workflow.log \
  type.cwl inputs.yml
```

### Estimated Runtime
Measured on our test runs (8 CPUs):
- **`extract-genes` with annotation:** about 4 minutes per assembly.
- **`repull-reads`:** about 5 hours per 30x whole-genome CRAM (most of it reading the whole file).
- **`type-t1k` + `patch-builder` + `score`:** about 13 minutes per sample (both reference arms).
- **Total per sample from CRAM:** about 5.5 hours, almost all of it read fetching.

---

## 5. Workflow Outputs

The output directory contains:
- `calls_genotype.tsv`: Base linear-reference typing calls from T1K.
- `final_calls.tsv`: Post-builder typing calls incorporating R4n overrides.
- `<gene>/additions.fa`: Predicted novel allele sequences for unrepresented haplotypes.
- `<gene>/proposals.json`: Candidate patch proposals and evidence metrics.
- `COMPLETE.json`: Run receipt with input file digests and execution timestamps.
