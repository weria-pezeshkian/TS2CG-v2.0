#!/bin/bash

## Tutorial 10: Simulating a Membrane 
#===============================================================================================
## Function to display help text
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
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
    shift
done

#===============================================================================================
files="files"  
#===============================================================================================
# Use SOL to solvate the system
#===============================================================================================

TS2CG SOL -in eq_v.gro -tem ./files/water.gro -o SOL.gro -Rcutoff 0.32 || { echo "Error: TS2CG SOL command failed"; exit 1; }

#===============================================================================================
## Edit the topology file
#===============================================================================================

> topol.top

printf "#include \"$files/martini3/martini_v3.0.4.itp\"\n" >> topol.top
printf "#include \"$files/martini3/martini_v3.0_phospholipids.itp\"\n" >> topol.top
printf "#include \"$files/martini3/martini_v3.0_solvents.itp\"\n" >> topol.top

if [ -f system.top ]; then
    cat system.top >> topol.top
else
    echo "Warning: system.top does not exist. Skipping append."
fi


topol_file="topol.top"
info_file="info.txt"

# Check if info.txt exists
if [ ! -f "$info_file" ]; then
    echo "$info_file does not exist."
    exit 1
fi

# Create a temporary file to store updated topol_file content
temp_file=$(mktemp)

# Remove matching lines from topol_file and save to temp_file
while IFS= read -r line; do
    grep -vF -- "$line" "$topol_file" >> "$temp_file"
done < "$info_file"

# Overwrite topol_file with updated content from temp_file
mv "$temp_file" "$topol_file"

# Append content of info.txt to topol.top
cat "$info_file" >> "$topol_file"
echo "Appended contents of $info_file to $topol_file."

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
    # Energy minimization
    #===============================================================================================

    ## regular
    $gmx_exec grompp -f $files/mdp/1DFourierShape/em_2.mdp -c SOL.gro -r SOL.gro -p topol.top -o output_1/em2_sol.tpr -maxwarn 1
    $gmx_exec mdrun -s output_1/em2_sol.tpr -v -deffnm output_1/em2_sol

    #===============================================================================================
    # Equilibration
    #===============================================================================================

    $gmx_exec grompp -f $files/mdp/1DFourierShape/eq_1.mdp -c output_1/em2_sol.gro -r output_1/em2_sol.gro -p topol.top -o output_1/eq_1.tpr
    $gmx_exec mdrun -v -s output_1/eq_1.tpr -deffnm output_1/eq_1

    $gmx_exec grompp -f $files/mdp/1DFourierShape/eq_2.mdp -c output_1/eq_1.gro -r output_1/eq_1.gro -p topol.top -o output_1/eq_2.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_1/eq_2.tpr -deffnm output_1/eq_2

    $gmx_exec grompp -f $files/mdp/1DFourierShape/eq_3.mdp -c output_1/eq_2.gro -r output_1/eq_2.gro -p topol.top -o output_1/eq_3.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_1/eq_3.tpr -deffnm output_1/eq_3

    $gmx_exec grompp -f $files/mdp/1DFourierShape/eq_4.mdp -c output_1/eq_3.gro -r output_1/eq_3.gro -p topol.top -o output_1/eq_4.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_1/eq_4.tpr -deffnm output_1/eq_4

    $gmx_exec grompp -f $files/mdp/1DFourierShape/eq_5.mdp -c output_1/eq_4.gro -r output_1/eq_4.gro -p topol.top -o output_1/eq_5.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_1/eq_5.tpr -deffnm output_1/eq_5

    #===============================================================================================
    # Production run
    #===============================================================================================

    $gmx_exec grompp -f $files/mdp/1DFourierShape/md.mdp -c output_1/eq_5.gro -p topol.top  -o output_1/md.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_1/md.tpr -deffnm output_1/md

else
    echo "Skipping the GROMACS part."
fi
