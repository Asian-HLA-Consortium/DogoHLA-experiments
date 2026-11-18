cwlVersion: v1.2
class: Workflow
label: DogoHLA Typing Workflow (repull-reads + type-t1k + patch-builder per sample)
doc: |
  Types a sample using the best DogoHLA configuration:
  1. Repulls missing HLA reads (Rule-2).
  2. Runs T1K against linear reference (IPD alone for Arm A, panel for Arm B).
  3. Runs patch-builder with R4n abstain rule to reconstruct missing alleles.

requirements:
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  donor:
    type: string
  cram:
    type: File?
  reference_genome:
    type: File?
  regional_r1:
    type: File?
  regional_r2:
    type: File?
  candidate_r1:
    type: File?
  candidate_r2:
    type: File?
  bait_kmers:
    type: File
  ipd_genomic:
    type: File
  panel_reference:
    type: File?
  confounded_dir:
    type: Directory
  threads:
    type: int
    default: 8

steps:
  repull:
    run: tools/repull_reads.cwl
    in:
      donor: donor
      cram: cram
      reference_genome: reference_genome
      ipd_genomic: ipd_genomic
      regional_r1: regional_r1
      regional_r2: regional_r2
      candidate_r1: candidate_r1
      candidate_r2: candidate_r2
      bait_kmers: bait_kmers
      threads: threads
    out: [repulled_r1, repulled_r2, pairs_json]

  type_armA:
    run: tools/type_t1k.cwl
    in:
      reads_r1: repull/repulled_r1
      reads_r2: repull/repulled_r2
      reference: ipd_genomic
      threads: threads
    out: [genotype_tsv, candidate_1, candidate_2]

  builder_armA:
    run: tools/patch_builder.cwl
    in:
      donor: donor
      pass1_calls: type_armA/genotype_tsv
      reference: ipd_genomic
      reads_r1: repull/repulled_r1
      reads_r2: repull/repulled_r2
      confounded_dir: confounded_dir
      arm:
        valueFrom: "C4R4n"
      threads: threads
    out: [final_calls_tsv, builder_dir]

  type_armB:
    run: tools/type_t1k.cwl
    when: $(inputs.reference != null)
    in:
      reads_r1: repull/repulled_r1
      reads_r2: repull/repulled_r2
      reference: panel_reference
      threads: threads
    out: [genotype_tsv, candidate_1, candidate_2]

  builder_armB:
    run: tools/patch_builder.cwl
    when: $(inputs.reference != null)
    in:
      donor: donor
      pass1_calls: type_armB/genotype_tsv
      reference: panel_reference
      reads_r1: repull/repulled_r1
      reads_r2: repull/repulled_r2
      confounded_dir: confounded_dir
      arm:
        valueFrom: "D4R4n"
      threads: threads
    out: [final_calls_tsv, builder_dir]

outputs:
  repulled_r1:
    type: File
    outputSource: repull/repulled_r1
  repulled_r2:
    type: File
    outputSource: repull/repulled_r2
  pairs_json:
    type: File
    outputSource: repull/pairs_json
  armA_calls_tsv:
    type: File
    outputSource: type_armA/genotype_tsv
  c4r4n_final_calls_tsv:
    type: File
    outputSource: builder_armA/final_calls_tsv
  armB_calls_tsv:
    type: File?
    outputSource: type_armB/genotype_tsv
  d4r4n_final_calls_tsv:
    type: File?
    outputSource: builder_armB/final_calls_tsv
