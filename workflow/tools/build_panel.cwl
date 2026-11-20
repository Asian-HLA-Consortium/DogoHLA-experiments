cwlVersion: v1.2
class: CommandLineTool
label: Build deduplicated composite T1K panel reference
doc: |
  Combines pinned IPD genomic alleles with cohort HLA sequences from extract-genes,
  filters excluded samples/relatives, removes IPD-identical alleles and duplicates,
  and synthesises arm B reference.fa with aliases.json.

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
      location: ../scripts/build_panel.py
  ipd_genomic:
    type: File
    doc: IPD-IMGT/HLA 3.65.0 genomic.fa
  records:
    type: File[]
    doc: Gene record JSON files from extract-genes
  exclude_samples:
    type: string?
    doc: Comma-separated list of sample IDs to exclude

baseCommand: [python3]

arguments:
  - $(inputs.script.path)
  - prefix: --ipd-genomic
    valueFrom: $(inputs.ipd_genomic.path)
  - prefix: --records
    valueFrom: $(inputs.records)
  - prefix: --exclude-samples
    valueFrom: $(inputs.exclude_samples || "")
  - prefix: --out-reference
    valueFrom: reference.fa
  - prefix: --out-aliases
    valueFrom: aliases.json
  - prefix: --out-summary
    valueFrom: panel_summary.json

outputs:
  reference:
    type: File
    outputBinding:
      glob: reference.fa
  aliases:
    type: File
    outputBinding:
      glob: aliases.json
  summary:
    type: File
    outputBinding:
      glob: panel_summary.json
