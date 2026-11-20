cwlVersion: v1.2
class: CommandLineTool
label: Extract 8 HLA gene sequences and exon annotations from assembly
doc: |
  Annotates HLA genes using Immuannot and extracts the 8 HLA genes with 2,000 bp
  flanks, parsing core sequences and CDS intervals. Output matches SOURCE_RECORDS format.

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 8000

hints:
  DockerRequirement:
    dockerPull: ghcr.io/asian-hla-consortium/dogohla-workflow:latest

inputs:
  script:
    type: File
    default:
      class: File
      location: ../scripts/extract_genes.py
  assembly:
    type: File
    doc: Phased assembly FASTA
  gtf:
    type: File
    doc: Immuannot GTF annotation
  sample:
    type: string
    doc: Sample identifier
  haplotype:
    type: string
    default: "1"
    doc: Haplotype identifier
  threads:
    type: int
    default: 4

baseCommand: [python3]

arguments:
  - $(inputs.script.path)
  - --assembly
  - $(inputs.assembly.path)
  - --gtf
  - $(inputs.gtf.path)
  - --sample
  - $(inputs.sample)
  - --haplotype
  - $(inputs.haplotype)
  - --out-fasta
  - $(inputs.sample)_$(inputs.haplotype)_genes.fa
  - --out-records
  - $(inputs.sample)_$(inputs.haplotype)_records.json

outputs:
  fasta:
    type: File
    outputBinding:
      glob: $(inputs.sample)_$(inputs.haplotype)_genes.fa
  records:
    type: File
    outputBinding:
      glob: $(inputs.sample)_$(inputs.haplotype)_records.json
