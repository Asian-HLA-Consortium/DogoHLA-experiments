cwlVersion: v1.2
class: Workflow
label: DogoHLA Evaluation Workflow (score over samples)
doc: |
  Scores final genotype calls against truth sequences across multiple samples,
  computing gene-body recovery of IPD-absent copies and named accuracy.

requirements:
  InlineJavascriptRequirement: {}

inputs:
  calls:
    type: File[]
    doc: Final calls TSVs across samples
  truth_records:
    type: File[]
    doc: Truth records JSONs from extract-genes
  family_table:
    type: File?
    doc: Optional FOLDS.json for family component clustering
  arm:
    type: string
    default: "D4R4n"

steps:
  score:
    run: tools/score.cwl
    in:
      calls: calls
      truth_records: truth_records
      family_table: family_table
      arm: arm
    out: [rows_tsv, summary_json]

outputs:
  score_rows:
    type: File
    outputSource: score/rows_tsv
  score_summary:
    type: File
    outputSource: score/summary_json
