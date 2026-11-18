cwlVersion: v1.2
class: CommandLineTool
label: Run patch-builder with R4n abstain rule
doc: |
  Discovers novel mosaic patches from mapped candidate reads, applies the R4n
  abstain rule to drop unresolved novel proposals, and produces final calls.

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 16000

hints:
  DockerRequirement:
    dockerPull: dogohla/workflow:1.0.0

inputs:
  script:
    type: File
    default:
      class: File
      location: ../scripts/patch_builder.py
  donor:
    type: string
  pass1_calls:
    type: File
  reference:
    type: File
  reads_r1:
    type: File
  reads_r2:
    type: File
  confounded_dir:
    type: Directory
  arm:
    type: string
    default: "C4R4n"
  threads:
    type: int
    default: 8

baseCommand: [python3]

arguments:
  - $(inputs.script.path)
  - prefix: --donor
    valueFrom: $(inputs.donor)
  - prefix: --pass1-calls
    valueFrom: $(inputs.pass1_calls.path)
  - prefix: --reference
    valueFrom: $(inputs.reference.path)
  - prefix: --reads-r1
    valueFrom: $(inputs.reads_r1.path)
  - prefix: --reads-r2
    valueFrom: $(inputs.reads_r2.path)
  - prefix: --confounded-dir
    valueFrom: $(inputs.confounded_dir.path)
  - prefix: --arm
    valueFrom: $(inputs.arm)
  - prefix: --out-dir
    valueFrom: builder_out
  - prefix: --out-final-calls
    valueFrom: $(inputs.donor)_$(inputs.arm)_final_calls.tsv
  - prefix: --threads
    valueFrom: $(inputs.threads)

outputs:
  final_calls_tsv:
    type: File
    outputBinding:
      glob: $(inputs.donor)_$(inputs.arm)_final_calls.tsv
  builder_dir:
    type: Directory
    outputBinding:
      glob: builder_out
