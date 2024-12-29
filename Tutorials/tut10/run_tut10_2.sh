#!/bin/bash

## Tutorial 10: Simulating a Membrane with Wall
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

> topol_2.top

printf "#include \"$files/martini3/martini_v3.0.4.itp\"\n" >> topol_2.top
printf "#include \"Wall.itp\"\n" >> topol_2.top
printf "#include \"$files/martini3/martini_v3.0_phospholipids.itp\"\n" >> topol_2.top
printf "#include \"$files/martini3/martini_v3.0_solvents.itp\"\n" >> topol_2.top


if [ -f system_2.top ]; then
    cat system_2.top >> topol_2.top
else
    echo "Warning: system.top does not exist. Skipping append."
fi


topol_file="topol_2.top"
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

# Append content of info.txt to topol_2.top
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

    # Check if the output_2 folder exists
    if [ -d "output_2" ]; then
    # Backup the existing output_2 folder by prefixing with #
       mv output_2 "#output_2"  # Rename the existing folder
       echo "Backup of existing output_2 folder created as #output_2"
    fi

    mkdir -p output_2

    set -o errexit
    set -o nounset

    #===============================================================================================
    # Energy minimization
    #===============================================================================================

    ## regular
    $gmx_exec grompp -f $files/mdp/Wall/em_2.mdp -c SOL.gro -r SOL.gro -p topol_2.top -o output_2/em2_sol.tpr -maxwarn 1
    $gmx_exec mdrun -s output_2/em2_sol.tpr -v -deffnm output_2/em2_sol

    #===============================================================================================
    # Making an index file
    #===============================================================================================

    # Input file (structure file, e.g., .gro or .pdb)
    structure_file="em2_sol.gro"
    # Output index file
    output_index="index.ndx"

    # Run GROMACS to create an index file
    echo "Running GROMACS make_ndx to create and name the Lipids group..."

    $gmx_exec make_ndx -f output_2/$structure_file -o output_2/$output_index << EOF

2|3
q
EOF
     echo "Index file created: $output_index"

    #===============================================================================================
    # Equilibration
    #===============================================================================================

    $gmx_exec grompp -f $files/mdp/Wall/eq_1.mdp -c output_2/em2_sol.gro -r output_2/em2_sol.gro -p topol_2.top -n output_2/index.ndx  -o output_2/eq_1.tpr
    $gmx_exec mdrun -v -s output_2/eq_1.tpr -deffnm output_2/eq_1

    $gmx_exec grompp -f $files/mdp/Wall/eq_2.mdp -c output_2/eq_1.gro -r output_2/eq_1.gro -p topol_2.top -n output_2/index.ndx  -o output_2/eq_2.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_2/eq_2.tpr -deffnm output_2/eq_2

    $gmx_exec grompp -f $files/mdp/Wall/eq_3.mdp -c output_2/eq_2.gro -r output_2/eq_2.gro -p topol_2.top -n output_2/index.ndx  -o output_2/eq_3.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_2/eq_3.tpr -deffnm output_2/eq_3

    $gmx_exec grompp -f $files/mdp/Wall/eq_4.mdp -c output_2/eq_3.gro -r output_2/eq_3.gro -p topol_2.top -n output_2/index.ndx  -o output_2/eq_4.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_2/eq_4.tpr -deffnm output_2/eq_4

    $gmx_exec grompp -f $files/mdp/Wall/eq_5.mdp -c output_2/eq_4.gro -r output_2/eq_4.gro -p topol_2.top -n output_2/index.ndx  -o output_2/eq_5.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_2/eq_5.tpr -deffnm output_2/eq_5

    #===============================================================================================
    # Production run
    #===============================================================================================

    $gmx_exec grompp -f $files/mdp/Wall/md.mdp -c output_2/eq_5.gro -p topol_2.top  -n output_2/index.ndx  -o output_2/md.tpr -maxwarn 1
    $gmx_exec mdrun -v -s output_2/md.tpr -deffnm output_2/md

else
    echo "Skipping the GROMACS part."
fi
