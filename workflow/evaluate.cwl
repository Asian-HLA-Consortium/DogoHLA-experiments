cwlVersion: v1.2
class: Workflow
label: DogoHLA Evaluation Workflow (score one donor, all arms)
doc: |
  Scores the outputs of type.cwl for one donor with the frozen scorer: one row per truth copy and arm,
  for the arms A, B, C4, C4R4n, D4, D4R4n (B, D4 and D4R4n only when the panel inputs are given).
  The truth of the donor and the family groups (FOLDS.json) come from the cohort resource folder,
  the same folder that type.cwl gets as builder_inputs.

requirements:
  InlineJavascriptRequirement: {}

inputs:
  donor:
    type: string
  calls_a:
    type: File
    doc: armA_calls_tsv of type.cwl
  calls_b:
    type: File?
    doc: armB_calls_tsv of type.cwl
  builder_dir_c:
    type: Directory
    doc: builder_dir_c of type.cwl
  builder_dir_d:
    type: Directory?
    doc: builder_dir_d of type.cwl
  panel_reference:
    type: File?
  ipd_genomic:
    type: File
  builder_inputs:
    type: Directory

steps:
  score:
    run: tools/score.cwl
    in:
      donor: donor
      calls_a: calls_a
      calls_b: calls_b
      builder_dir_c: builder_dir_c
      builder_dir_d: builder_dir_d
      panel_reference: panel_reference
      ipd_genomic: ipd_genomic
      builder_inputs: builder_inputs
    out: [rows_tsv, avail_tsv]

outputs:
  score_rows:
    type: File
    outputSource: score/rows_tsv
  score_avail:
    type: File
    outputSource: score/avail_tsv
