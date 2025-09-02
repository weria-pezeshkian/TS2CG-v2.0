#!/bin/bash

# Tutorial 8_2: Place proteins based on curvature
# This script generates a vesicle with proteins using TS2CG and runs the outputs with GROMACS.

#===============================================================================================
# Function to display help text
#===============================================================================================

function show_help {
    echo "Usage: script.sh [-h]"
    echo "No options will only execute the TS2CG part of the Tutorial."
    echo "Options:"
    echo "  -h         Display this help message and exit."
    exit
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h) 
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
    shift
done

#===============================================================================================
files="files"  
#===============================================================================================
# TS2CG Commands
#===============================================================================================

echo "Running TS2CG steps..."

# Generate point file with bilayer thickness
TS2CG PLM -TSfile Small_Curved.tsi -bilayerThickness 3.8 -rescalefactor 4 4 4 -o point_SmallCurved || {
    echo "Error: TS2CG PLM failed."
    exit 1
}

# Add proteins
TS2CG INU -p point_SmallCurved -n 10 -t 1 -r 5 -c -0.5 -k 100 -o point_SmallCurved_new || {
    echo "Error: TS2CG INU failed."
    exit 1
}

# Build the system
TS2CG PCG -str input_SmallCurved.str -Bondlength 0.2 -LLIB "$files/Martini3.LIB" -dts ./point_SmallCurved_new -incdirtype Local -defout system_SmallCurved || {
    echo "Error: TS2CG PCG failed."
    exit 1
}

#===============================================================================================
# Edit the topology file
#===============================================================================================

echo "Editing topology file..."

topol_file="topol.top"
> "$topol_file"

printf "#include \"$files/martini3/martini_v3.0.4.itp\"\n" >> "$topol_file"
printf "#include \"$files/martini3/protein1.itp\"\n" >> "$topol_file"
printf "#include \"$files/martini3/protein2.itp\"\n" >> "$topol_file"
printf "#include \"$files/martini3/martini_v3.0_phospholipids.itp\"\n" >> "$topol_file"

if [ -f system_SmallCurved.top ]; then
    cat system_SmallCurved.top >> "$topol_file"
else
    echo "Warning: system.top does not exist. Skipping append."
fi

#===============================================================================================
# Prompt the user to execute the GROMACS part
#===============================================================================================

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

     if [ -d "output_2" ]; then
     mv output_2 "#output_2"
     fi

     mkdir -p output_2

     set -o errexit
     set -o nounset

    #===============================================================================================
    # Vacuum energy minimization
    #===============================================================================================

## regular
    $gmx_exec grompp -f $files/mdp/vesicle/em_2.mdp -c system_SmallCurved.gro -r system_SmallCurved.gro -p topol.top -o output_2/em_2.tpr
    $gmx_exec mdrun -s output_2/em_2.tpr -v -deffnm output_2/em_2

    #===============================================================================================
    # Vacuum equilibration
    #===============================================================================================
 
 $gmx_exec grompp -f $files/mdp/vesicle/eq_v.mdp -c output_2/em_2.gro -p topol.top -r output_2/em_2.gro -o output_2/eq_v.tpr  -maxwarn 1 
    $gmx_exec mdrun -v -s output_2/eq_v.tpr -deffnm output_2/eq_v 

else
    echo "Skipping the GROMACS part."
fi
