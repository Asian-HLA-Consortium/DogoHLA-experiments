cwlVersion: v1.2
class: CommandLineTool
label: Run T1K short-read HLA typing
doc: |
  Runs T1K against the specified linear reference (IPD genomic or augmented panel).
  Produces calls_genotype.tsv and candidate read alignments.

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

baseCommand: [run-t1k]

arguments:
  - --preset
  - hla-wgs
  - --alleleDigitUnits
  - "4"
  - --alleleDelimiter
  - ":"
  - -o
  - calls

outputs:
  genotype_tsv:
    type: File
    outputBinding:
      glob: calls_genotype.tsv
  allele_tsv:
    type: File?
    outputBinding:
      glob: calls_allele.tsv
  candidate_1:
    type: File?
    outputBinding:
      glob: calls_candidate_1.fq
  candidate_2:
    type: File?
    outputBinding:
      glob: calls_candidate_2.fq
