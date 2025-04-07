#include <iostream>
#include <iomanip>  // For std::setw and std::left
#include "SimDef.h"
#include "help.h"

help::help(std::string version, std::string exe)
{
    int size = exe.size();

    if (exe.at(size - 1) == 'M' && exe.at(size - 2) == 'L' && exe.at(size - 3) == 'P')
    {
        std::cout << "==========================================================================================================\n\n";
        std::cout << "-- Pointillism\n";
        std::cout << "-- Version:  " << SoftWareVersion << "\n";
        std::cout << "-- Groningen Biomolecular Sciences and Biotechnology Institute and Zernike Institute for Advanced Materials,\n-- University of Groningen, Groningen, Netherlands\n";
        std::cout << "-- For more information contact Weria Pezeshkian: w.pezeshkian@rug.nl and weria.pezeshkian@gmail.com\n";
        std::cout << "-- citation: Pezeshkian, W., König, M., Wassenaar, T.A. et al. Backmapping triangulated surfaces to coarse-grained membrane models. Nat Commun 11, 2296 (2020).\n";
        std::cout << "-- citation: F. Schuhmann, J. A. Stevens, N. Rahmani, I. Lindahl, C. M. Brown, C. Brasnett, D. Anastasiou, A. Bravo Vidal, B. Geiger, S. J. Marrink and W. Pezeshkian. TS2CG as a membrane builder. (2025).\n\n";
        std::cout << "==========================================================================================================\n\n";
        std::cout << "-- Pointillism (PLM), reads a triangulated surface input file e.g., typical outputs of dynamically \n-- triangulated surfaces simulations, and generates two sets of points that represent two smooth surfaces \n-- (upper and lower monolayers of a bilayer).\n";
        std::cout << "-- Simplest triangulated surfaces file that can be read by this script should be as: \n\n";
        std::cout << "     Box_X   Box_Y   Box_Z\n";
        std::cout << "     No_of_vertex\n";
        std::cout << "        0   Ver1_X   Ver1_Y   Ver1_Z\n";
        std::cout << "        1   Ver2_X   Ver2_Y   Ver2_Z\n";
        std::cout << "      ...    ...      ...      ...\n";
        std::cout << "      ...    ...      ...      ...\n";
        std::cout << "        i   Veri_X   Veri_Y   Veri_Z\n";
        std::cout << "      ...  ... ... ...\n";
        std::cout << "      No_of_triangles\n";
        std::cout << "        0   Ver1_id   Ver2_id   Ver3_id\n";
        std::cout << "        1   Ver1_id   Ver2_id   Ver3_id\n\n";
        std::cout << "-- Note: the number of the output points is always larger or equal to the number of the vertices \n-- in the input triangulated surface. With option -Mashno you can tune how many points you want (No_of_vertex*4^Mashno).\n";
        std::cout << "-- option -r check allows you to get some information about the input triangulated surface.\n";
        std::cout << "-- option -o allows you to specify a unique folder for the output files.\n";
        std::cout << "-- option -smooth allows you to handle highly rough surfaces; However, try to avoid this option, since it has not been tested well yet.\n";
        std::cout << "-- option -monolayer allows you to generate a monolayer instead of a bilayer. It should be either 1 or -1, which can invert the shape of the monolayer.\n\n";
        std::cout << "==========================================================================================================\n\n";
        std::cout << "-- Pointillism: generating double layer points from a triangulated surface\n\n";
        std::cout << "==========================================================================================================\n\n";
        std::cout << "-------------------------------------------------------------------------------\n";
        std::cout << std::left << std::setw(20) << "  option"
                  << std::setw(15) << "type"
                  << std::setw(20) << "default"
                  << "description\n";
        std::cout << "-------------------------------------------------------------------------------\n";

        // Print the options with proper alignment
        std::cout << std::left << std::setw(20) << Def_rescalefactor
                  << std::setw(15) << "rx ry rz"
                  << std::setw(20) << "(1 1 1)"
                  << "rescaling factor (Note: there is no guarantee to get proper surface if rx, ry, rz are not equal).\n";

        std::cout << std::left << std::setw(20) << Def_BilayerThickness
                  << std::setw(15) << "double"
                  << std::setw(20) << "3.8"
                  << "bilayer thickness\n";

        std::cout << std::left << std::setw(20) << Def_Monolayer
                  << std::setw(15) << "int"
                  << std::setw(20) << "0"
                  << "to generate monolayer instead (1/-1)\n";

        std::cout << std::left << std::setw(20) << Def_TaskName
                  << std::setw(15) << "string"
                  << std::setw(20) << "PLM"
                  << "function(PLM/in_out/check/add_pbc)\n";

        std::cout << std::left << std::setw(20) << Def_SmoothingFlag
                  << std::setw(15) << "------"
                  << std::setw(20) << "no"
                  << "might be necessary for rough surfaces\n";

        std::cout << std::left << std::setw(20) << Def_OutputFolderName
                  << std::setw(15) << "string"
                  << std::setw(20) << "point"
                  << "name of the output folder\n";

        std::cout << std::left << std::setw(20) << Def_resizebox
                  << std::setw(15) << "------"
                  << std::setw(20) << "no"
                  << "find a better box for the system\n";

        std::cout << std::left << std::setw(20) << Def_TSfile
                  << std::setw(15) << "string"
                  << std::setw(20) << "TS.tsi"
                  << "TS file name (three file format types, *.q/*.tsi/*.dat)\n";

        std::cout << std::left << std::setw(20) << Def_DegreeOfMeshing
                  << std::setw(15) << "int"
                  << std::setw(20) << "1"
                  << "number of Mosaicing, your point number grows as 4^Mashno\n";
        
        std::cout << std::left << std::setw(20) << Def_PrintLessPutput
                  << std::setw(15) << "bool"
                  << std::setw(20) << "false"
                  << "print less outputs\n";
        
        std::cout << std::left << std::setw(20) << Def_AlgType
                  << std::setw(15) << "string"
                  << std::setw(20) << "Type1"
                  << "algorithm type for Mosaicing (Type1 and Type2); no difference has been reported yet\n\n";

        std::cout << "-- Note: using Mashno [1-4], unless you know what you are doing.\n\n";
        std::cout << "-- basic example: PLM  -TSfile Traj1.tsi -bilayerThickness 4  -rescalefactor 3 3 3\n";
    }
}

help::~help()
{
    // Destructor implementation (if needed)
}
