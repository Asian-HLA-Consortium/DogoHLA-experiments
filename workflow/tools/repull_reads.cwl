cwlVersion: v1.2
class: CommandLineTool
label: Repull HLA reads via MHC regional pull and bait 31-mer matching
doc: |
  Rule-2 read set generation: MHC regional pull UNION genome-wide pairs
  with >= 3 bait 31-mers. Accepts WGS CRAM or paired FASTQ.

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
      location: ../scripts/repull_reads.py
    inputBinding:
      position: 1
  donor:
    type: string
    inputBinding:
      prefix: --donor
      position: 2
  bait_kmers:
    type: File
    inputBinding:
      prefix: --bait
      position: 3
  threads:
    type: int
    default: 8
    inputBinding:
      prefix: --threads
      position: 4
  cram:
    type: File?
    inputBinding:
      prefix: --cram
      position: 10
  reference_genome:
    type: File?
    inputBinding:
      prefix: --reference-genome
      position: 11
  ipd_genomic:
    type: File?
    inputBinding:
      prefix: --ipd-genomic
      position: 12
  regional_r1:
    type: File?
    inputBinding:
      prefix: --regional-r1
      position: 13
  regional_r2:
    type: File?
    inputBinding:
      prefix: --regional-r2
      position: 14
  candidate_r1:
    type: File?
    inputBinding:
      prefix: --candidate-r1
      position: 15
  candidate_r2:
    type: File?
    inputBinding:
      prefix: --candidate-r2
      position: 16

baseCommand: [python3]

arguments:
  - prefix: --out-r1
    valueFrom: $(inputs.donor)_r1.fq.gz
    position: 20
  - prefix: --out-r2
    valueFrom: $(inputs.donor)_r2.fq.gz
    position: 21
  - prefix: --out-pairs
    valueFrom: $(inputs.donor)_PAIRS.json
    position: 22

outputs:
  repulled_r1:
    type: File
    outputBinding:
      glob: $(inputs.donor)_r1.fq.gz
  repulled_r2:
    type: File
    outputBinding:
      glob: $(inputs.donor)_r2.fq.gz
  pairs_json:
    type: File
    outputBinding:
      glob: $(inputs.donor)_PAIRS.json
