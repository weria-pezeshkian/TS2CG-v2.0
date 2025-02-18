#include <iostream>
#include "SimDef.h"
#include "help.h"
help::help(std::string version,std::string exe)
{
    int size= exe.size();

    if(exe.at(size-1)=='M' && exe.at(size-2)=='L' && exe.at(size-3)=='P')
    {
        std::cout<<"=========================================================================================================="<<"\n\n";
        std::cout<<"-- Pointillism  "<<"\n";
        std::cout<<"-- Version:  "<<SoftWareVersion<<"\n";
        std::cout<<"-- Groningen Biomolecular Sciences and Biotechnology Institute and Zernike Institute for Advanced Materials,\n-- University of Groningen, Groningen, Netherlands"<<"\n";
        std::cout<<"-- For more information contact Weria Pezeshkian: w.pezeshkian@rug.nl and weria.pezeshkian@gmail.com"<<"\n";
        std::cout<<"-- citation: Pezeshkian, W., König, M., Wassenaar, T.A. et al. Backmapping triangulated surfaces to coarse-grained membrane models. Nat Commun 11, 2296 (2020)."<<"\n\n";
        std::cout<<"=========================================================================================================="<<"\n\n";
        std::cout<<"-- Pointillism (PLM), reads a triangulated surface input file e.g., typical outputs of dynamically \n-- triangulated surfaces simulations, and generates two sets of points that represent a two smooth surfaces \n-- (upper and lower monolayers of a bilayer). "<<"\n";
        std::cout<<"-- Simplest triangulated surfaces file that can be read by this script should be as: "<<"\n\n";
        std::cout<<"     Box_X   Box_Y   Box_Z "<<"\n";
        std::cout<<"     No_of_vertex  "<<"\n";
        std::cout<<"        0   Ver1_X   Ver1_Y   Ver1_Z"<<"\n";
        std::cout<<"        1   Ver2_X   Ver2_Y   Ver2_Z"<<"\n";
        std::cout<<"      ...    ...      ...      ..."<<"\n";
        std::cout<<"      ...    ...      ...      ..."<<"\n";
        std::cout<<"        i   Veri_X   Veri_Y   Veri_Z"<<"\n";
        std::cout<<"      ...  ... ... ..."<<"\n";
        std::cout<<"      No_of_triangles "<<"\n";
        std::cout<<"        0   Ver1_id   Ver2_id   Ver3_id"<<"\n";
        std::cout<<"        1   Ver1_id   Ver2_id   Ver3_id"<<"\n\n";
        std::cout<<"-- Note: the number of the output points is always larger or equal the number of the vertices \n-- in the input triangulated surface.  With option -Mashno you can tune how many points you want (No_of_vertex*4^Mashno) "<<"\n";
        std::cout<<"-- option -r check allows you can get some information about the input triangulated surface   "<<"\n";
        std::cout<<"-- option -o allows you to specify a unique folder for the output files    "<<"\n";
        std::cout<<"-- option -smooth allows you to handle highly rough surfaces; However, try to avoid this option, since has it not been tested well yet.    "<<"\n";
        std::cout<<"-- option -monolayer allows you to generate a monolayer instead of bilayer. It should be either 1 or -1, which can invert the shape of the monolayer    "<<"\n\n";
        std::cout<<"=========================================================================================================="<<"\n\n";
        std::cout<<"-- Pointillism: generating double layer points from a closed triangulated surface   "<<"\n\n";
        std::cout<<"=========================================================================================================="<<"\n\n";
        std::cout<<"-------------------------------------------------------------------------------"<<"\n";
        std::cout<<"  option                type        default            description "<<"\n";
        std::cout<<"-------------------------------------------------------------------------------"<<"\n";
        std::cout<<"  -rescalefactor        rx ry rz      (1 1 1)           rescaling factor (Note: there is no guarantee to get proper surface if rx,ry,rz are not equal. "<<"\n";
        std::cout<<"  -bilayerThickness     double         3.8              bilayer thickness "<<"\n";
        std::cout<<"  -monolayer            int             0               to generate monolayer instead (1/-1)  "<<"\n";
        std::cout<<"  -r                    string         PLM              function(PLM/in_out/check/add_pbc) "<<"\n";
        std::cout<<"  -smooth               ------         no               might be necessary for rough surfaces "<<"\n";
        std::cout<<"  -o                    string         point            name of the output folder "<<"\n";
        std::cout<<"  -resizebox            ------         no               find a better box for the system  "<<"\n";
        std::cout<<"  -TSfile               string         TS.tsi           TS file name (three file format type, *.q/*.tsi/*.dat)"<<"\n";
        std::cout<<"  -Mashno               int             1               number of Mosaicing, your point number grows as 4^Mashno  "<<"\n";
        std::cout<<"  -AlgType              string         Type1            algorithm type for Mosaicing (Type1 and Type2); no difference has been reported yet  "<<"\n\n";
        std::cout<<"-- Note: using Mashno [1-4], unless you know what you are doing.    "<<"\n\n";
        std::cout<<"-- basic example: PLM  -TSfile Traj1.tsi -bilayerThickness 4  -rescalefactor 3 3 3 \n";



    }


}

help::~help()
{

}
