cwlVersion: v1.2
class: CommandLineTool
label: Score HLA calls against truth sequences
doc: |
  Scores final genotype calls against truth sequences from extract-genes.
  Calculates gene-body recovery of IPD-absent copies, 4-field named accuracy,
  and family-clustered bootstrap confidence intervals.

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 2
    ramMin: 4000

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
  calls:
    type: File[]
    inputBinding:
      prefix: --calls
      position: 2
  truth_records:
    type: File[]
    inputBinding:
      prefix: --truth-records
      position: 3
  family_table:
    type: File?
    inputBinding:
      prefix: --family-table
      position: 4
  arm:
    type: string
    default: "final"
    inputBinding:
      prefix: --arm
      position: 5

baseCommand: [python3]

arguments:
  - prefix: --out-rows
    valueFrom: score_rows.tsv
  - prefix: --out-summary
    valueFrom: score_summary.json

outputs:
  rows_tsv:
    type: File
    outputBinding:
      glob: score_rows.tsv
  summary_json:
    type: File
    outputBinding:
      glob: score_summary.json
