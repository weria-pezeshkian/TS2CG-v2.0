#include <iostream>
#include <iomanip>  // For std::setw and std::left
#include "help.h"
#include "Def.h"

help::help(std::string exe)
{
    int size = exe.size();
    {
        std::cout << "\n";
        std::cout << "==========================================================================================================\n\n";
        std::cout << "-- " << SoftWareName << "\n";
        std::cout << "-- Version:  " << SoftWareVersion << "\n";
        std::cout << "-- Niels Bohr International Academy, \n Niels Bohr Institute, \n University of Copenhagen, Copenhagen, Denmark\n";
        std::cout << "-- For more information contact Weria Pezeshkian: w.pezeshkian@nbi.ku.dk and weria.pezeshkian@gmail.com\n";
        std::cout << "-- citation: Pezeshkian, W., König, M., Wassenaar, T.A. et al. Backmapping triangulated surfaces to coarse-grained membrane models. Nat Commun 11, 2296 (2020).\n\n";
        std::cout << "-- citation: F. Schuhmann, J. A. Stevens, N. Rahmani, I. Lindahl, C. M. Brown, C. Brasnett, D. Anastasiou, A. Bravo Vidal, B. Geiger, S. J. Marrink and W. Pezeshkian. TS2CG as a membrane builder. (2025).\n\n";
        std::cout << "==========================================================================================================\n\n";
        std::cout << "-- With option -Bondlength, you can change the initial bond guess. Large Bondlength may generate an unstable structure.\n";
        std::cout << "-- With option -renorm, the molar ratio of the lipid will be renormalized.\n";
        std::cout << "==========================================================================================================\n\n";
        std::cout << "------------ This script converts Pointillism outputs to a CG model -------------------\n";
        std::cout << "-------------------------------------------------------------------------------\n";
        std::cout << std::left << std::setw(20) << "  option"
                  << std::setw(15) << "type"
                  << std::setw(20) << "default"
                  << "description\n";
        std::cout << "-------------------------------------------------------------------------------\n";

        // Convert macro-defined STRExt to std::string for concatenation
        std::string inputFileDefault = "input." + std::string(STRExt);

        // Print the options with proper alignment
        std::cout << std::left << std::setw(20) << G_POINT_FOLDER
                  << std::setw(15) << "string"
                  << std::setw(20) << "point"
                  << "dts folder address\n";

        std::cout << std::left << std::setw(20) << G_STR_FILE_TAG
                  << std::setw(15) << "string"
                  << std::setw(20) << inputFileDefault
                  << "input file\n";

        std::cout << std::left << std::setw(20) << G_DEFAULT_TAG
                  << std::setw(15) << "string"
                  << std::setw(20) << "output"
                  << "output files prefix\n";

        std::cout << std::left << std::setw(20) << G_BOND_LENGTH
                  << std::setw(15) << "double"
                  << std::setw(20) << "0.1"
                  << "initial bond guess\n";

        std::cout << std::left << std::setw(20) << LipidLibraryFileName
                  << std::setw(15) << "string"
                  << std::setw(20) << "no"
                  << "CG lipid library file name\n";

        std::cout << std::left << std::setw(20) << G_renormalized_lipid_ratio
                  << std::setw(15) << "------"
                  << std::setw(20) << "no"
                  << "renormalized the lipid molar ratio\n";

        std::cout << std::left << std::setw(20) << "-iter"
                  << std::setw(15) << "double"
                  << std::setw(20) << "4"
                  << "the number of point selection is iter*number of the point\n";

        std::cout << std::left << std::setw(20) << INC_DIR_TYPE
                  << std::setw(15) << "string"
                  << std::setw(20) << "Global"
                  << "the type of protein direction data (Local/Global)\n";

        std::cout << std::left << std::setw(20) << "-Wall"
                  << std::setw(15) << "------"
                  << std::setw(20) << "off"
                  << "a flag to create a wall around the membrane\n";

        std::cout << std::left << std::setw(20) << FunctionType
                  << std::setw(15) << "string"
                  << std::setw(20) << "backmap"
                  << "backmap/analytical_shape\n";

        std::cout << std::left << std::setw(20) << "-WallBName"
                  << std::setw(15) << "string"
                  << std::setw(20) << "WL"
                  << "Name of the Wall beads\n";

        std::cout << std::left << std::setw(20) << "-WPointDir"
                  << std::setw(15) << "bool"
                  << std::setw(20) << "false"
                  << "Just write the folder\n";

        std::cout << std::left << std::setw(20) << G_SKIP_LIPID_PLACEMENT
                  << std::setw(15) << "bool"
                  << std::setw(20) << "false"
                  << "skipping lipid placement step\n";

        std::cout << std::left << std::setw(20) << G_RADIUS_CUT_OFF
                  << std::setw(15) << "double"
                  << std::setw(20) << "0.5"
                  << "distance between proteins and lipids\n";

        std::cout << std::left << std::setw(20) << G_KEEP_POINTS_CLOSE_TO_PROTEINS
                  << std::setw(15) << "bool"
                  << std::setw(20) << "false"
                  << "keep points close to proteins\n";
        
        std::cout << std::left << std::setw(20) << G_PRINT_LESS_OUTPUTS
                  << std::setw(15) << "bool"
                  << std::setw(20) << "false"
                  << "print less outputs\n";
        std::cout << "=========================================================================== \n";
        std::cout << "basic example:  "<<ExecutableName<<" "<<G_POINT_FOLDER<<"  point "<<G_STR_FILE_TAG<<" input.str \n";
    }
}

help::~help()
{
    // Destructor implementation (if needed)
}

