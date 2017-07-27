# Generate a mica substrate of the desired size in angstroms

cwlVersion: v1.0
class: Workflow
inputs:
    nx: int
    ny: int
    nz: int
    output_name: string
    car_file: File
    mdf_file: File
    ff_file: File # forcefield

outputs:
    classout:
        type: File
        outputSource: lammps_data/$(inputs.output_name)

steps:
    msi2lmp:
        run: msi2lmp.cwl
    extract_coords:
        run: extract_coords.cwl
    duplicate:
        run: duplicate.cwl



      
