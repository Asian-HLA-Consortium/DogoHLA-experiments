cwlVersion: v1.2
class: CommandLineTool
label: Immuannot HLA Gene Annotation on Phased Assembly FASTA
doc: |
  Annotates HLA genes using frozen Immuannot (commit 31362c370123f6706a04ffdda6a3bf0c108328e0)
  and pinned reference bundle (Data-ipd3590-kir2130). Produces a compressed GTF (<sample>.<haplotype>.gtf.gz).

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
      location: ../scripts/annotate_immuannot.py
  assembly:
    type: File
    doc: Phased assembly FASTA or MHC span FASTA
  immuannot_ref:
    type: Directory
    doc: Immuannot reference data directory (e.g. Data-ipd3590-kir2130)
  sample:
    type: string
    doc: Sample / donor identifier
  haplotype:
    type: string
    default: "1"
    doc: Haplotype identifier (1 or 2)
  threads:
    type: int
    default: 4
    doc: Number of minimap2 compute threads

baseCommand: [python3]

arguments:
  - $(inputs.script.path)
  - --assembly
  - $(inputs.assembly.path)
  - --immuannot-ref
  - $(inputs.immuannot_ref.path)
  - --sample
  - $(inputs.sample)
  - --haplotype
  - $(inputs.haplotype)
  - --threads
  - $(inputs.threads)

outputs:
  gtf:
    type: File
    outputBinding:
      glob: "*.gtf.gz"
