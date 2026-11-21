cwlVersion: v1.2
class: CommandLineTool
label: Run the frozen patch-builder (one arm)
doc: |
  Runs the frozen all-gene patch builder on the pass-1 T1K run of one arm (arm C = IPD reference,
  arm D = IPD plus cohort panel) and applies its overrides to the pass-1 calls. The script is a thin
  wrapper: the builder itself is the vendored frozen code in dogohla/frozen, checked against
  MANIFEST.sha256 at every start.
  `arm` is C4, C4R4n, D4 or D4R4n; the final calls file holds the calls of that arm. The builder folder
  `builder_<arm>` holds the overrides of every variant (OVERRIDES.json = C4/D4, OVERRIDES_R4n.json = C4R4n/D4R4n).

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
      location: ../scripts/patch_builder.py
    inputBinding:
      position: 1
  package:
    type: Directory
    default:
      class: Directory
      location: ../dogohla
    inputBinding:
      prefix: --package-dir
      position: 2
    doc: The workflow package folder (verify.py, wrapper_env.py, frozen/)
  donor:
    type: string
    inputBinding:
      prefix: --donor
      position: 3
  arm:
    type: string
    default: "C4R4n"
    inputBinding:
      prefix: --arm
      position: 4
  pass1_calls:
    type: File
    inputBinding:
      prefix: --pass1-calls
      position: 5
  pass1_candidate_r1:
    type: File
    inputBinding:
      prefix: --pass1-candidate-r1
      position: 6
  pass1_candidate_r2:
    type: File
    inputBinding:
      prefix: --pass1-candidate-r2
      position: 7
  reference:
    type: File
    inputBinding:
      prefix: --reference
      position: 8
    doc: T1K reference FASTA of this arm (IPD genomic.fa for C, the panel reference.fa for D)
  aliases:
    type: File?
    inputBinding:
      prefix: --aliases
      position: 9
    doc: aliases.json of the panel (arm D)
  ipd_genomic:
    type: File
    inputBinding:
      prefix: --ipd-genomic
      position: 10
  builder_inputs:
    type: Directory
    inputBinding:
      prefix: --builder-inputs
      position: 11
    doc: Cohort resource folder (SOURCE_RECORDS.json, FOLDS.json, mosaic_classes.tsv, panels_atlas.py, msa/, INPUTS.sha256.tsv)
  params:
    type: File
    default:
      class: File
      location: ../resources/builder/PARAMS.json
    inputBinding:
      prefix: --params
      position: 12
  confounded_dir:
    type: Directory
    inputBinding:
      prefix: --confounded-dir
      position: 13
  threads:
    type: int
    default: 8
    inputBinding:
      prefix: --threads
      position: 14

baseCommand: [python3]

arguments:
  - prefix: --out-dir
    valueFrom: builder_$(inputs.arm)
    position: 20
  - prefix: --out-final-calls
    valueFrom: $(inputs.donor)_$(inputs.arm)_final_calls.tsv
    position: 21

outputs:
  final_calls_tsv:
    type: File
    outputBinding:
      glob: $(inputs.donor)_$(inputs.arm)_final_calls.tsv
  builder_dir:
    type: Directory
    outputBinding:
      glob: builder_$(inputs.arm)
