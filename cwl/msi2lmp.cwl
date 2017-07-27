cwlVersion: v1.0
class: CommandLineTool
label: Convert Materials Studio car/mdf files + forcefield into LAMMPS data file

# requirements:
#   InitialWorkDirRequirement:
#     listing:
#       - entryname: result.txt
#         entry: |
#             outdir=$(runtime.outdir)

baseCommand: /home/oliver/academic/research/droplet/gen_droplet/scripts/msi2lmp.exe
inputs:
    material_name: 
        type: string
        inputBinding:
            position: 1
    class_number: # Use 1
        type: int
        inputBinding:
            position: 2
            prefix: -class
    verbosity: # Use 3
        type: int
        inputBinding:
            position: 3
            prefix: -print
    forcefield: 
        type: string
        inputBinding:
            position: 4
            prefix: -frc
    ignore_errors: # Use true
        type: boolean
        inputBinding:
            position: 5
            prefix: -i

outputs: []
# outputs:
#     out_file: 
#         type: File
#         outputBinding:
#             glob: result.txt
