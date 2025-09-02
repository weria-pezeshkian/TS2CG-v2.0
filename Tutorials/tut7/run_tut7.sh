#!/bin/bash

#===============================================================================================
files="files"  
#===============================================================================================

# Step 1: Create a point-folder using PCG, only requires .str in the working directory.

echo "Running TS2CG PCG with initial .str file..."

TS2CG PCG -str input.str -function analytical_shape -WPointDir

if [ $? -ne 0 ]; then
    echo "Error: TS2CG PCG (point-folder creation) failed."
    exit 1
fi
# ===============================================================================================

# Step 2: Modify the point-folder using DOP, requires domain_input.txt which sets c0 for the lipids.

echo "Modifying point-folder using TS2CG DOP with domain_input.txt..."

TS2CG DOP -s domain_input.txt -ni input_DOP.str -k 10

if [ $? -ne 0 ]; then
    echo "Error: TS2CG DOP (point-folder modification) failed."
    exit 1
fi
# ===============================================================================================

# Step 3: Build the bilayer based on the modified point-folder and lipid library.

echo "Building bilayer using TS2CG PCG with modified .str file and Martini3.LIB..."

TS2CG PCG -dts point -str input_DOP.str -defout system -LLIB $files/Martini3.LIB  

if [ $? -ne 0 ]; then
    echo "Error: TS2CG PCG (bilayer building) failed."
    exit 1
fi

# ===============================================================================================
# Edit the topology file
# ===============================================================================================

> topol.top
printf "#include \"$files/martini3/martini_v3.0.4.itp\"\n" >> topol.top
printf "#include \"$files/martini3/martini_v3.0_phospholipids.itp\"\n" >> topol.top
printf "#include \"files/martini3/martini_v3.0_CDLs.itp\"\n" >> topol.top

if [ -f system.top ]; then
    cat system.top >> topol.top
else
    echo "Warning: system.top does not exist. Skipping append."
fi

# ===============================================================================================
# Prompt the user to execute the GROMACS part
# ===============================================================================================

read -p "Do you want to execute the GROMACS simulation part? (yes/no): " user_input

if [[ "$user_input" == "yes" ]]; then
     if which gmx &> /dev/null; then
         gmx_exec=gmx
     elif which gmx_d &> /dev/null; then
         gmx_exec=gmx_d
     elif which gmx_mpi &> /dev/null; then
         gmx_exec=gmx_mpi
     else
         echo "GROMACS executable not found. Please install it and try again."
         exit 1
     fi

     if [ -d "output" ]; then
     mv output "#output"
     fi

     mkdir -p output

     set -o errexit
     set -o nounset

    # ===============================================================================================
    # Vacuum energy minimization
    # ===============================================================================================

    ## soft-core
    $gmx_exec grompp -f $files/mdp/1DFourierShape/em_1.mdp -c system.gro -r system.gro -p topol.top -o output/em_1.tpr
    $gmx_exec mdrun -s output/em_1.tpr -v -deffnm output/em_1

    ## regular
    $gmx_exec grompp -f $files/mdp/1DFourierShape/em_2.mdp -c output/em_1.gro -r output/em_1.gro -p topol.top -o output/em_2.tpr
    $gmx_exec mdrun -s output/em_2.tpr -v -deffnm output/em_2

    #===============================================================================================
    # Vacuum equilibration
    #===============================================================================================

    $gmx_exec grompp -f $files/mdp/1DFourierShape/eq_v.mdp -c output/em_2.gro -p topol.top -r output/em_2.gro -o output/eq_v.tpr 
    $gmx_exec mdrun -v -deffnm output/eq_v

else
    echo "Skipping the GROMACS part."
fi



