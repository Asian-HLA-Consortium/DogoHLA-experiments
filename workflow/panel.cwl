cwlVersion: v1.2
class: Workflow
label: DogoHLA Panel Workflow (annotate-immuannot + extract-genes + build-panel)
doc: |
  Constructs an augmented cohort reference panel from haplotype-resolved assemblies.
  Runs Immuannot on each assembly FASTA to generate HLA annotations GTF,
  scatters extract-genes over cohort assemblies and GTFs,
  and compiles deduplicated arm-B reference.fa and aliases.json.

requirements:
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  assemblies:
    type: File[]
    doc: Phased assembly FASTAs or MHC span FASTAs
  immuannot_ref:
    type: Directory
    doc: Immuannot reference bundle (Data-ipd3590-kir2130)
  samples:
    type: string[]
    doc: Sample IDs
  haplotypes:
    type: string[]
    doc: Haplotype identifiers (1 or 2)
  ipd_genomic:
    type: File
    doc: IPD-IMGT/HLA 3.65.0 genomic.fa
  exclude_samples:
    type: string?
    doc: Optional comma-separated sample IDs to exclude

steps:
  annotate:
    run: tools/annotate_immuannot.cwl
    scatter: [assembly, sample, haplotype]
    scatterMethod: dotproduct
    in:
      assembly: assemblies
      immuannot_ref: immuannot_ref
      sample: samples
      haplotype: haplotypes
    out: [gtf]

  extract:
    run: tools/extract_genes.cwl
    scatter: [assembly, gtf, sample, haplotype]
    scatterMethod: dotproduct
    in:
      assembly: assemblies
      gtf: annotate/gtf
      sample: samples
      haplotype: haplotypes
    out: [fasta, records]

  build:
    run: tools/build_panel.cwl
    in:
      ipd_genomic: ipd_genomic
      records: extract/records
      exclude_samples: exclude_samples
    out: [reference, aliases, summary]

outputs:
  panel_reference:
    type: File
    outputSource: build/reference
  aliases:
    type: File
    outputSource: build/aliases
  summary:
    type: File
    outputSource: build/summary
  extracted_records:
    type: File[]
    outputSource: extract/records
  annotations:
    type: File[]
    outputSource: annotate/gtf
