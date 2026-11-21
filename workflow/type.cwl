cwlVersion: v1.2
class: Workflow
label: DogoHLA Typing Workflow (repull-reads + type-t1k + patch-builder per sample)
doc: |
  Types one sample with the best DogoHLA configuration:
  1. Repulls the HLA reads (regional MHC recruitment with mate rescue + genome-wide bait pairs) from the CRAM.
  2. Runs T1K on the IPD reference (arm A) and, when a panel reference is given, on the panel (arm B).
  3. Runs the frozen patch-builder on each T1K run and applies the R4n abstain rule (C4R4n, D4R4n).
  The CRAM index (.crai) must sit next to the CRAM and the reference FASTA needs its .fai next to it.
  Outputs: both T1K tables, the builder folders of both arms, the final calls of both arms, the repulled reads.

requirements:
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  donor:
    type: string
  cram:
    type: File?
    secondaryFiles:
      - pattern: .crai
        required: $(self.basename.endsWith('.cram'))
      - pattern: .bai
        required: $(self.basename.endsWith('.bam'))
  reference_genome:
    type: File?
    secondaryFiles:
      - pattern: .fai
        required: true
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
    doc: Arm-B T1K reference (IPD plus cohort panel). Without it arm B and D4R4n are skipped.
  panel_aliases:
    type: File?
    doc: aliases.json that belongs to panel_reference (needed with it)
  builder_inputs:
    type: Directory
    doc: Cohort resource folder for the builder (SOURCE_RECORDS.json, FOLDS.json, msa/, ...)
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
      out_prefix:
        valueFrom: "armA"
    out: [genotype_tsv, candidate_1, candidate_2]

  builder_armA:
    run: tools/patch_builder.cwl
    in:
      donor: donor
      pass1_calls: type_armA/genotype_tsv
      pass1_candidate_r1: type_armA/candidate_1
      pass1_candidate_r2: type_armA/candidate_2
      reference: ipd_genomic
      ipd_genomic: ipd_genomic
      builder_inputs: builder_inputs
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
      out_prefix:
        valueFrom: "armB"
    out: [genotype_tsv, candidate_1, candidate_2]

  builder_armB:
    run: tools/patch_builder.cwl
    when: $(inputs.reference != null)
    in:
      donor: donor
      pass1_calls: type_armB/genotype_tsv
      pass1_candidate_r1: type_armB/candidate_1
      pass1_candidate_r2: type_armB/candidate_2
      reference: panel_reference
      aliases: panel_aliases
      ipd_genomic: ipd_genomic
      builder_inputs: builder_inputs
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
  armB_calls_tsv:
    type: File?
    outputSource: type_armB/genotype_tsv
  c4r4n_final_calls_tsv:
    type: File
    outputSource: builder_armA/final_calls_tsv
  d4r4n_final_calls_tsv:
    type: File?
    outputSource: builder_armB/final_calls_tsv
  builder_dir_c:
    type: Directory
    outputSource: builder_armA/builder_dir
  builder_dir_d:
    type: Directory?
    outputSource: builder_armB/builder_dir
