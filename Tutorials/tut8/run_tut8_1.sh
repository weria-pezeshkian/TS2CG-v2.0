#!/bin/bash

# Tutorial 8_1: Adding more Proteins to the Membrane
# This script generates a POPC vesicle containing 8 proteins (five copies of protein1 and three copies of protein2)
# and runs the TS2CG outputs using GROMACS.

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

TS2CG PLM -TSfile Sphere.tsi -bilayerThickness 3.8 -rescalefactor 4 4 4 -o point_vesicle  || {
    echo "Error: TS2CG PLM failed."
    exit 1
}

# Add proteins

TS2CG INU -p point_vesicle -n 3 -t 1 -r 5 -o point_vesicle_new || {
    echo "Error: TS2CG INU for protein type 1 failed."
    exit 1
}

TS2CG INU -p point_vesicle_new -n 2 -t 2 -r 5 -o point_vesicle_new2 || {
    echo "Error: TS2CG INU for protein type 2 failed."
    exit 1
}

# Build the system

TS2CG PCG -str input_vesicle.str -Bondlength 0.2 -LLIB "$files/Martini3.LIB" -dts ./point_vesicle_new2 -incdirtype Local -defout system_vesicle || {

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

if [ -f system_vesicle.top ]; then
    cat system_vesicle.top >> "$topol_file"
else
    echo "Warning: system_vesicle.top does not exist. Skipping append."
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

     if [ -d "output_1" ]; then
     mv output_1 "#output_1"
     fi

     mkdir -p output_1

     set -o errexit
     set -o nounset

    #===============================================================================================
    # Vacuum energy minimization
    #===============================================================================================
  
## regular
    $gmx_exec grompp -f $files/mdp/vesicle/em_2.mdp -c system_vesicle.gro -r system_vesicle.gro -p topol.top -o output_1/em_2.tpr
    $gmx_exec mdrun -s output_1/em_2.tpr -v -deffnm output_1/em_2

    #===============================================================================================
    # Vacuum equilibration
    #===============================================================================================
  
  $gmx_exec grompp -f $files/mdp/vesicle/eq_v.mdp -c output_1/em_2.gro -p topol.top -r output_1/em_2.gro -o output_1/eq_v.tpr  -maxwarn 1 
    $gmx_exec mdrun -v -s output_1/eq_v.tpr -deffnm output_1/eq_v 

else
    echo "Skipping the GROMACS part."
fi
