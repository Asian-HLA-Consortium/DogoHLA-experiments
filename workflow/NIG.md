# Running DogoHLA Workflow on NIG Supercomputer (Japan)

This document provides instructions for collaborator Hassan Sibroe Daanaa to execute the DogoHLA CWL v1.2 workflow on the National Institute of Genetics (NIG) supercomputer in Mishima, Japan.

---

## 1. Environment & Prerequisites at NIG

The workflow is portable and contains zero hardcoded cluster paths. It can be run either via Singularity/Apptainer or inside a local micromamba environment.

### Option A: Singularity / Apptainer (Recommended)
NIG compute nodes support Singularity containers.
```bash
# Pull the published workflow container from GHCR into the folder where cwltool looks for it.
# Temporary and cache folders must be on an ordinary filesystem (for example your home directory),
# not on an encrypted or ACL-managed one; otherwise `singularity exec` can fail with
# "bad superblock for squashfs".
export SINGULARITY_TMPDIR=$HOME/sif-tmp SINGULARITY_CACHEDIR=$HOME/sif-cache CWL_SINGULARITY_CACHE=$HOME/sif-images
mkdir -p $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $CWL_SINGULARITY_CACHE
singularity pull --name $CWL_SINGULARITY_CACHE/ghcr.io_asian-hla-consortium_dogohla-workflow:latest.sif docker://ghcr.io/asian-hla-consortium/dogohla-workflow:latest
```
`cwltool --singularity` finds the image by this file name in `CWL_SINGULARITY_CACHE`. Without it, cwltool
pulls the image itself (about 10 minutes, no message) with default temporary folders. Use the `APPTAINER_`
variables if you run Apptainer.

### Option B: Native Environment via Micromamba
If running uncontainerized (`--no-container`), create the environment using:
```bash
micromamba create -n dogohla -c conda-forge -c bioconda \
  python=3.12 samtools=1.21 minimap2=2.28 bwa jellyfish toil \
  pysam edlib numpy scipy pandas
```

---

## 2. Resource Preparation on NIG Storage

Download and build the public reference resources with the helper script, inside the image:
```bash
singularity exec $CWL_SINGULARITY_CACHE/ghcr.io_asian-hla-consortium_dogohla-workflow:latest.sif bash scripts/fetch_resources.sh resources
```
This stages, each checked by sha256:
- `resources/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa` and its `.fai`: GRCh38 reference with decoys and HLA contigs.
- `resources/ipd/genomic.fa`: IPD-IMGT/HLA 3.65.0 genomic alleles (checked, not built: place the pinned file supplied by the consortium there first).
- `resources/bait/bait_kmers.txt`: MHC 31-mer bait dictionary (`sha256: 50b3b7324a5f5bca...`).
- `resources/immuannot/Data-ipd3590-kir2130`: Immuannot reference data.

Each CRAM or BAM needs its index next to it (`samtools index sample.cram` gives `sample.cram.crai`); the
workflow stops at the start if it is missing. The builder also needs the cohort resource folder
(`builder_inputs`) that the consortium supplies.

Place your local Japanese WGS CRAMs/BAMs and phased assemblies in your project directory (e.g. `/path/to/data/`). Note that Japanese donor data must not leave NIG.

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
export TMPDIR=${SLURM_TMPDIR:-$HOME/tmp}
mkdir -p "$TMPDIR"
export SINGULARITY_CACHEDIR=$HOME/sif-cache SINGULARITY_TMPDIR=$HOME/sif-tmp CWL_SINGULARITY_CACHE=$HOME/sif-images

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

# Score the sample afterwards with evaluate.cwl (see the README)

echo "Completed DogoHLA workflow at $(date)"
```

---

## 4. Input Configuration YAML

Example `SAMPLE_JPT01_inputs.yml` (the CRAM index `.crai` and the reference `.fai` are found next to their files):
```yaml
donor: "SAMPLE_JPT01"
cram:
  class: File
  path: "/path/to/nig/data/SAMPLE_JPT01.cram"
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
  path: "/path/to/workflow/resources/confounded_v1"
threads: 16
```
Leave out `panel_reference` and `panel_aliases` to run arm A only.

---

## 5. Verification & Integrity

The workflow automatically validates that all vendored frozen modules match their registered SHA-256 hashes at runtime before proceeding. If any file fails verification, execution halts immediately with an integrity assertion error.
