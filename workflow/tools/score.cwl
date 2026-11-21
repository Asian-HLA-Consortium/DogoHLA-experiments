cwlVersion: v1.2
class: CommandLineTool
label: Score the typing arms of one donor with the frozen scorer
doc: |
  Runs the frozen all-gene scorer on one donor and writes one row per truth copy and arm
  (all columns of the frozen scorer) for the arms A, B (T1K tables), C4, C4R4n (builder on the IPD
  run) and D4, D4R4n (builder on the panel run). Arms B, D4 and D4R4n need the panel inputs.
  The script is a thin wrapper around the vendored frozen code in dogohla/frozen.

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 2
    ramMin: 8000

hints:
  DockerRequirement:
    dockerPull: ghcr.io/asian-hla-consortium/dogohla-workflow:latest

inputs:
  script:
    type: File
    default:
      class: File
      location: ../scripts/score.py
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
  donor:
    type: string
    inputBinding:
      prefix: --donor
      position: 3
  calls_a:
    type: File
    inputBinding:
      prefix: --calls-a
      position: 4
    doc: T1K genotype table on the IPD reference (arm A)
  calls_b:
    type: File?
    inputBinding:
      prefix: --calls-b
      position: 5
    doc: T1K genotype table on the panel reference (arm B)
  builder_dir_c:
    type: Directory
    inputBinding:
      prefix: --builder-dir-c
      position: 6
  builder_dir_d:
    type: Directory?
    inputBinding:
      prefix: --builder-dir-d
      position: 7
  panel_reference:
    type: File?
    inputBinding:
      prefix: --panel-reference
      position: 8
  ipd_genomic:
    type: File
    inputBinding:
      prefix: --ipd-genomic
      position: 9
  builder_inputs:
    type: Directory
    inputBinding:
      prefix: --builder-inputs
      position: 10
  params:
    type: File
    default:
      class: File
      location: ../resources/builder/PARAMS.json
    inputBinding:
      prefix: --params
      position: 11

baseCommand: [python3]

arguments:
  - prefix: --out-rows
    valueFrom: score_rows.tsv
    position: 20
  - prefix: --out-avail
    valueFrom: score_avail.tsv
    position: 21

outputs:
  rows_tsv:
    type: File
    outputBinding:
      glob: score_rows.tsv
  avail_tsv:
    type: File
    outputBinding:
      glob: score_avail.tsv
