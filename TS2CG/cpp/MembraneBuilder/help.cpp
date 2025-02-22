#include <iostream>
#include "help.h"
#include "Def.h"

help::help(std::string exe)
{
    int size= exe.size();
 {
     

    std::cout<<"\n";
     std::cout<<"=========================================================================================================="<<"\n\n";
     std::cout<<"-- "<<SoftWareName<<"\n";
     std::cout<<"-- Version:  "<<SoftWareVersion<<"\n";
     std::cout<<"-- Niels Bohr International Academy, \n Niels Bohr Institute, \n University of Copenhagen, Copenhagen, Denmark"<<"\n";
     std::cout<<"-- For more information contact Weria Pezeshkian: w.pezeshkian@nbi.ku.dk and weria.pezeshkian@gmail.com"<<"\n";
     std::cout<<"-- citation: Pezeshkian, W., König, M., Wassenaar, T.A. et al. Backmapping triangulated surfaces to coarse-grained membrane models. Nat Commun 11, 2296 (2020)."<<"\n\n";
     std::cout<<"=========================================================================================================="<<"\n\n";
     std::cout<<"-- With option -Bondlength, you can chnage the initial bond guess. Large Bondlength may generate an unstable structure ";
     std::cout<<"-- With  option -renorm  the molar ratio of the lipid will be renormalized  "<<"\n";
     std::cout<<"-- To get higher denisty, you may increase -Mashno value or reduce -ap value in PLM command.   "<<"\n\n";
     std::cout<<"=========================================================================================================="<<"\n\n";
     std::cout<<"------------ This script convert Pointillism outputs to a CG model -------------------"<<"\n";
     std::cout<<"-------------------------------------------------------------------------------"<<"\n";
     std::cout<<"  option            type        default            description "<<"\n";
     std::cout<<"-------------------------------------------------------------------------------"<<"\n";
     std::cout<<G_POINT_FOLDER<<"           string       point                dts folder address "<<"\n";
     std::cout<<G_STRExt<<"             string       input."<<STRExt<<"            input file "<<"\n";
     std::cout<<G_DEFAULT_TAG<<"        string       output               output files prefix "<<"\n";
     std::cout<<G_BOND_LENGTH<<"       double       0.1                  initial bond guess;  "<<"\n";
     std::cout<<LipidLibraryFileName<<"            string       no                   CG lipid library file name;  "<<"\n";
     std::cout<<G_renormalized_lipid_ratio<<"           ------       no                   renormalized the lipid molar ratio  "<<"\n";
     std::cout<<" -iter             double       4                    the number of point selection is iter*number of the point  "<<"\n";
     std::cout<<INC_DIR_TYPE<<"     string       Global               the type of protein direction data (Local/Global)  "<<"\n";
     std::cout<<" -Wall             ------       off                  a flag to create a wall around the membrane  "<<"\n";
     std::cout<<FunctionType<<"         string       backmap              backmap/analytical_shape  "<<"\n";
     std::cout<<" -WallBName        string       WL                   Name of the Wall beads  "<<"\n\n";
     std::cout<<" -WPointDir        bool         false                Just write the folder  "<<"\n\n";
     std::cout<<G_SKIP_LIPID_PLACEMENT<<"    bool         false                skipping lipid placement step  "<<"\n\n";
     std::cout<<G_RADIUS_CUT_OFF<<"    double         0.5              distance between proteins and lipids  "<<"\n\n";
     std::cout<<G_KEEP_POINTS_CLOSE_TO_PROTEINS<<"    bool         false            keep points close to proteins  "<<"\n\n";

     
//analytical_shape
     
     
    
     std::cout<< "basic example: PCG -dts point -str input.str -seed 39234  -Bondlength 0.15 "<<"\n\n";




   }

    

}

help::~help()
{
    
}
