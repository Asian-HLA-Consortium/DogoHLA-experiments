cwlVersion: v1.2
class: CommandLineTool
label: Run T1K short-read HLA typing
doc: |
  Runs T1K against the specified linear reference (IPD genomic or augmented panel).
  Produces <out_prefix>_genotype.tsv and the candidate reads <out_prefix>_candidate_1.fq / _2.fq.

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 16000

hints:
  DockerRequirement:
    dockerPull: ghcr.io/asian-hla-consortium/dogohla-workflow:latest

inputs:
  reads_r1:
    type: File
    inputBinding:
      prefix: "-1"
      position: 1
  reads_r2:
    type: File
    inputBinding:
      prefix: "-2"
      position: 2
  reference:
    type: File
    inputBinding:
      prefix: -f
      position: 3
  threads:
    type: int
    default: 8
    inputBinding:
      prefix: -t
      position: 4
  out_prefix:
    type: string
    default: calls
    inputBinding:
      prefix: -o
      position: 5
    doc: Prefix of the T1K output files (use a different one per arm, e.g. armA / armB)

baseCommand: [run-t1k]

arguments:
  - --preset
  - hla-wgs
  - --alleleDigitUnits
  - "4"
  - --alleleDelimiter
  - ":"

outputs:
  genotype_tsv:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_genotype.tsv
  allele_tsv:
    type: File?
    outputBinding:
      glob: $(inputs.out_prefix)_allele.tsv
  candidate_1:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_candidate_1.fq
  candidate_2:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_candidate_2.fq
