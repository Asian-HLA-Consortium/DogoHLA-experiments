cwlVersion: v1.2
class: Workflow
label: DogoHLA part 1 for one assembly haplotype (annotate-immuannot + extract-genes)
doc: |
  Immuannot annotation of one phased assembly (or MHC span) FASTA, then extraction of the 8 genes.
  Used by scripts/run_panel.sh, which runs many of these as separate processes.

requirements:
  InlineJavascriptRequirement: {}

inputs:
  assembly:
    type: File
  immuannot_ref:
    type: Directory
  sample:
    type: string
  haplotype:
    type: string

steps:
  annotate:
    run: tools/annotate_immuannot.cwl
    in:
      assembly: assembly
      immuannot_ref: immuannot_ref
      sample: sample
      haplotype: haplotype
    out: [gtf]
  extract:
    run: tools/extract_genes.cwl
    in:
      assembly: assembly
      gtf: annotate/gtf
      sample: sample
      haplotype: haplotype
    out: [fasta, records]

outputs:
  gtf:
    type: File
    outputSource: annotate/gtf
  records:
    type: File
    outputSource: extract/records
