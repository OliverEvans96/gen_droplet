#!/bin/bash
# Generate mica - tie together individual steps

#  1. Start with ${name}ural data from STRUCTURE.car and STRUCTURE.mdf (must have same name)  and forcefield data from FF.frc
#  2. Generate a LAMMPS input file from the ${name}ure files using:
#      ./msi2lmp.exe STRUCTURE -class 1 -print 3 -frc FF -i
#          -msi2lmp can be found in this folder or in your LAMMPS directory
#          -both ${name}ure and forcefield files can be provided without suffixes
#          -class 1 for 12/6 LJ, class 2 for 9/6
#          -print 3 for verbose output
#          -i ignores errors and treats them as warnings and proceeds to generate a file
#          -> This generates STRUCTURE.data (or STRUCTURE.lammps05 if using msi2lmp as provided by INTERFACE)
#  3. Create a simplified ${name}ure file by running:
#      ./car2coords.awk STRUCTURE.car COORDS.car
#  4. Duplicate unit cell as many times as desired by running:
#      lattice.py NX NY NZ COORDS.car STRUCTURE.data DATAFILE.data
#          -NX, NY, NZ are the number of times that the cell will be repeated in the x y and z directions.
#  5. That's it! (hopefully) DATAFILE.data should be readable by LAMMPS using read_data.

nx=10
ny=20
nz=1

name="mica_${nx}x${ny}"

BASE='/home/oliver/academic/research/droplet/gen_droplet'
BASE='..'

STRUCTURE="${BASE}/materials/mica15_single_layer"
COORDS="${BASE}/tmp/mica15_coords"
FF="${BASE}/forcefields/cvff_interface_v1_5.frc"
DATAFILE="${BASE}/lammps_data/${name}.data"

echo
echo "Step 0: Copy ${name}ure files here."
cp "${STRUCTURE}.car" "${name}.car"
cp "${STRUCTURE}.mdf" "${name}.mdf"


echo
echo "Step 1: Create initial LAMMPS data file"
echo ../scripts/msi2lmp_gcc32.exe "${name}" -class 1 -print 3 -frc "$FF"
../scripts/msi2lmp_gcc32.exe "${name}" -class 1 -print 3 -frc "$FF"

echo "Step 1.5: Delete copies of structures."
echo rm ${name}.car ${name}.mdf
rm ${name}.car ${name}.mdf

echo
echo "Step 2: Simplify coordinates for duplication"
echo ../scripts/car2coords.awk "${STRUCTURE}.car" "${COORDS}.car"
../scripts/car2coords.awk "${STRUCTURE}.car" "${COORDS}.car"

echo
echo "Step 3: Duplicate and generate LAMMPS file"
echo python ../scripts/lattice.py "$nx" "$ny" "$nz" "${COORDS}.car" "${name}.lammps05" "${DATAFILE}"
python ../scripts/lattice.py "$nx" "$ny" "$nz" "${COORDS}.car" "${name}.lammps05" "${DATAFILE}"

echo
echo "Delete temporary files"
echo rm -f "${name}.lammps05" "${COORDS}.car"
rm -f "struct.lammps05" "${COORDS}.car"

echo
echo "Generated '${DATAFILE}'"

