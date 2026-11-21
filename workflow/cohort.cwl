cwlVersion: v1.2
class: Workflow
label: DogoHLA cohort resources for new-person scoring (annotate + extract-genes + build-cohort)
doc: |
  Builds the `builder_inputs` folder for type.cwl and evaluate.cwl from the user's own phased assemblies
  (the truth) and the explicit panel. Per assembly haplotype: Immuannot annotation (part 1), extraction of
  the 8 genes (part 1), then build-cohort. For typing only, without assemblies, run tools/build_cohort.cwl
  directly (leave out `truth_records`).

requirements:
  ScatterFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  assemblies:
    type: File[]
    doc: Phased assembly FASTAs (one per haplotype) or MHC span FASTAs
  samples:
    type: string[]
    doc: Sample id of each assembly (same order as assemblies)
  haplotypes:
    type: string[]
    doc: Haplotype ("1" or "2") of each assembly
  immuannot_ref:
    type: Directory
  ipd_genomic:
    type: File
  panel_reference:
    type: File?
  panel_aliases:
    type: File?
  donors:
    type: string[]
    doc: Sample ids to be typed and scored
  families:
    type: File?
  mafft_dir:
    type: Directory?
  truth_resource:
    type: string?
  threads:
    type: int
    default: 8

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
  cohort:
    run: tools/build_cohort.cwl
    in:
      donors: donors
      ipd_genomic: ipd_genomic
      panel_reference: panel_reference
      panel_aliases: panel_aliases
      truth_records: extract/records
      families: families
      mafft_dir: mafft_dir
      truth_resource: truth_resource
      threads: threads
    out: [cohort]

outputs:
  cohort:
    type: Directory
    outputSource: cohort/cohort
  truth_records:
    type: File[]
    outputSource: extract/records
