# Running DogoHLA Workflow on NIG Supercomputer (Japan)

This document provides instructions for collaborator Hassan Sibroe Daanaa to execute the DogoHLA CWL v1.2 workflow on the National Institute of Genetics (NIG) supercomputer in Mishima, Japan.

---

## 1. Environment & Prerequisites at NIG

The workflow is portable and contains zero hardcoded cluster paths. It can be run either via Singularity/Apptainer or inside a local micromamba environment.

### Option A: Singularity / Apptainer (Recommended)
NIG compute nodes support Singularity containers.
```bash
# Pull published workflow container from GHCR
singularity pull docker://ghcr.io/asian-hla-consortium/dogohla-workflow:latest
```

### Option B: Native Environment via Micromamba
If running uncontainerized (`--no-container`), create the environment using:
```bash
micromamba create -n dogohla -c conda-forge -c bioconda \
  python=3.12 samtools=1.21 minimap2=2.28 bwa jellyfish toil \
  pysam edlib numpy scipy pandas
```

---

## 2. Resource Preparation on NIG Storage

Download the public reference resources using the staged helper script:
```bash
bash scripts/fetch_resources.sh
```
This fetches:
- `genomic.fa`: IPD-IMGT/HLA 3.65.0 genomic alleles.
- `bait_kmers.txt`: MHC 31-mer bait dictionary (`sha256: 50b3b7324a5f5bca...`).
- `Data-ipd3590-kir2130`: Immuannot reference models.

Place your local Japanese WGS CRAMs/BAMs and phased assemblies in your project directory (e.g. `/home/hdaanaa/hla/data/`). Note that Japanese donor data must not leave NIG.

---

## 3. Slurm Submission Script (`run_dogohla_nig.sbatch`)

Create the submission script tailored to NIG's Slurm configuration:

```bash
#!/bin/bash
#SBATCH --job-name=dogohla_wf
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/dogohla_%j.out
#SBATCH --error=logs/dogohla_%j.err

set -euo pipefail

# Environment setup
export TMPDIR=${SLURM_TMPDIR:-/tmp}
export SINGULARITY_CACHEDIR=$HOME/.singularity/cache

SAMPLE_ID="SAMPLE_JPT01"
WORKFLOW_DIR="$HOME/dogohla/workflow"
INPUTS="$WORKFLOW_DIR/examples/${SAMPLE_ID}_inputs.yml"
OUTDIR="$HOME/dogohla/results/${SAMPLE_ID}"
JOBSTORE="$HOME/dogohla/jobstore-${SAMPLE_ID}"

mkdir -p "$OUTDIR" logs

echo "Starting DogoHLA workflow for $SAMPLE_ID on $(hostname) at $(date)"

# Execute with Toil CWL runner using Singularity
toil-cwl-runner \
  --batchSystem single_machine \
  --maxCores 16 \
  --maxMemory 60G \
  --jobStore "$JOBSTORE" \
  --workDir "$TMPDIR" \
  --singularity \
  --clean onSuccess \
  --logLevel INFO \
  --logFile "$OUTDIR/toil.log" \
  --outdir "$OUTDIR" \
  "$WORKFLOW_DIR/type.cwl" "$INPUTS"

echo "Completed DogoHLA workflow at $(date)"
```

---

## 4. Input Configuration YAML

Example `SAMPLE_JPT01_inputs.yml`:
```yaml
sample_id: "SAMPLE_JPT01"
cram_or_bam:
  class: File
  path: "/path/to/nig/data/SAMPLE_JPT01.cram"
reference_genome:
  class: File
  path: "/path/to/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
ipd_genomic:
  class: File
  path: "/path/to/resources/ipd/genomic.fa"
bait_kmers:
  class: File
  path: "/path/to/resources/bait/bait_kmers.txt"
confounded_dir:
  class: Directory
  path: "/path/to/workflow/resources/confounded_v1"
threads: 16
```

---

## 5. Verification & Integrity

The workflow automatically validates that all vendored frozen modules match their registered SHA-256 hashes at runtime before proceeding. If any file fails verification, execution halts immediately with an integrity assertion error.
