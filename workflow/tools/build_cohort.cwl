cwlVersion: v1.2
class: CommandLineTool
label: Build the cohort resource folder for new-person typing and scoring
doc: |
  Builds, from the explicit panel (reference.fa + aliases.json) and, optionally, the part-1 records of the
  user's own phased assemblies, the folder that the frozen builder and scorer read as `builder_inputs`:
  SOURCE_RECORDS.json, FOLDS.json, mosaic_classes.tsv, panels_atlas.py, msa/ and INPUTS.sha256.tsv.
  Panel sequences enter as records of the donor "PANEL" (the panel must not contain the typed sample or its
  relatives). Truth records are used only by the scorer. Donors without a family table entry are singletons.
  Alignment: MAFFT 7.526 `--auto` on the unique sequences, then the vendored atlas classifier (dogohla/atlas).

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 16000

hints:
  DockerRequirement:
    dockerPull: ghcr.io/asian-hla-consortium/dogohla-workflow:latest

inputs:
  script:
    type: File
    default:
      class: File
      location: ../scripts/build_cohort.py
    inputBinding:
      position: 1
  package:
    type: Directory
    default:
      class: Directory
      location: ../dogohla
    inputBinding:
      prefix: --package-dir
      position: 2
  donors:
    type: string[]
    inputBinding:
      prefix: --donors
      position: 3
    doc: Sample ids that will be typed (and scored when truth records are given)
  ipd_genomic:
    type: File
    inputBinding:
      prefix: --ipd-genomic
      position: 4
  panel_reference:
    type: File?
    inputBinding:
      prefix: --panel-reference
      position: 5
  panel_aliases:
    type: File?
    inputBinding:
      prefix: --panel-aliases
      position: 6
  truth_records:
    type: File[]?
    inputBinding:
      prefix: --truth-records
      position: 7
    doc: Part-1 record JSON files (extract-genes) of the user's assemblies
  families:
    type: File?
    inputBinding:
      prefix: --families
      position: 8
    doc: Optional table donor<TAB>family id (same id = relatives); without it every donor is a singleton
  truth_resource:
    type: string?
    inputBinding:
      prefix: --truth-resource
      position: 9
    doc: Label for the `cohort` column of the score rows (default USER)
  mafft_dir:
    type: Directory?
    inputBinding:
      prefix: --mafft-prefix
      position: 11
    doc: Folder that contains bin/mafft (MAFFT 7.526); not needed when the image already has mafft
  threads:
    type: int
    default: 8
    inputBinding:
      prefix: --threads
      position: 10

baseCommand: [python3]

arguments:
  - prefix: --out-dir
    valueFrom: cohort
    position: 20

outputs:
  cohort:
    type: Directory
    outputBinding:
      glob: cohort
