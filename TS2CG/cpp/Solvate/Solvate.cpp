


#include "Argument.h"
#include "Solvate.h"
#include "Nfunction.h"
#include "GroFile.h"
#include "GenerateUnitCells.h"
#include "Solvate.h"

Solvate::Solvate(Argument *pArg)
{
    Nfunction f; // several function to be used
    
// get the input parameters from the argument class
    std::string ingrofilename  = pArg->GetIn_GroFileName();
    std::string outgrofilename = pArg->GetOut_GroFileName();
    std::string temfilename = pArg->GetTemplateFileName();
    std::string NegName = pArg->GetNegativeIonName();
    std::string PosName = pArg->GetPositiveIonName();
    double cutoff = pArg->GetRCutOff();
    int seed = pArg->GetSeed();
    std::vector<int> ion = pArg->GetIon();
    double db = pArg->GetDB();
    double usize  = pArg->GetUCELLSize();

//======
    
    if(f.FileExist(ingrofilename)==false)
    {
        std::cout<<"---> error: file "<<ingrofilename<< " do not exist. \n";
        exit(0);
    }
    if(f.FileExist(temfilename)==false)
    {
        std::cout<<"file "<<temfilename<< " do not exist. Error: template file should be given  \n";
        exit(0);
    }

        GroFile InGro = GroFile(ingrofilename); // read the gro file
        std::vector<bead*> Sysbead = InGro.GetpAllBeads(); // get all the beads in the system gro file in the vector=
        Vec3D *FBox = InGro.GetBox(); // get the box info
        Bring2Box(Sysbead,FBox);  // removing box crossing of the beads. Grofile could have it

    //-- generate unit cells to check the distance between the created solvent beads and the system beads.
        GenerateUnitCells UCELL(Sysbead,pArg,FBox, cutoff, usize);
        
        //==  read the template water gro file that will be used to put water beads
        GroFile TemGro = GroFile(temfilename);
        std::vector<bead*> Wbead = TemGro.GetpAllBeads();
        Vec3D *WBox = TemGro.GetBox();  // box of the template water beads box. much smaller then FBox
        Bring2Box(Wbead,WBox); // removing box crossing of the water beads. Grofile could have it
    //-- info to see how many copy of the water box is needed. some will go out of the box size, so it should be counted for after
        int nBox_X = int((*FBox)(0)/(*WBox)(0))+1;
        int nBox_Y = int((*FBox)(1)/(*WBox)(1))+1;
        int nBox_Z = int((*FBox)(2)/(*WBox)(2))+1;
        
    //-- a vector to store all the generated water beads
        std::vector<bead> FullWaterBead;
        for (int i=0;i<nBox_X;i++)
        for (int j=0;j<nBox_Y;j++)
        for (int k=0;k<nBox_Z;k++)
        {
            for (std::vector<bead *>::iterator it = Wbead.begin() ; it != Wbead.end(); ++it)
            {
                double x=(*it)->GetXPos()+((*WBox)(0))*double(i)+db;
                double y=(*it)->GetYPos()+((*WBox)(1))*double(j)+db;
                double z=(*it)->GetZPos()+((*WBox)(2))*double(k)+db;

                Vec3D Pos(x,y,z);
                if(UCELL.anythingaround (Pos)!=true)// remove beads that overlaps with system beads
                if(x>0 && y>0 && z>0 && x<(*FBox)(0) && y<(*FBox)(1) && z<(*FBox)(2))// remove beads that are not inside the box
                {

                        bead TB = *(*it);
                
                        TB.UpdatePos(FBox,x,y,z);
                        FullWaterBead.push_back(TB);
                }
            }
        }///
//== FullWaterBead is now being filled with water beads; note beads that are crossing the box is removed and also the one which overlaps with the system beads
    //create a function for ion placement, we may choose different placement
    
    // adding ions
    std::vector<bead> PreBeads = InGro.GetAllBeads();
    // generate ions and return ions and water beads.
    std::vector<bead> pWaterIonBeads = AddIons(FullWaterBead, ion.at(0),ion.at(1), PosName, NegName, seed);
    // add the water beads to the main beads container for gro productions
    for (std::vector<bead>::iterator it = pWaterIonBeads.begin() ; it != pWaterIonBeads.end(); ++it)
            PreBeads.push_back(*it);
    
    // write the final file
    TemGro.RenewBeads(PreBeads);
    TemGro.UpdateBox(*FBox);
    TemGro.WriteGroFile(outgrofilename);
}
Solvate::~Solvate()
{
    
}
// a function to bring all the beads outside of the box, inside of the box
void Solvate::Bring2Box(std::vector<bead*>& Sysbead, Vec3D* pBox) {
    
    for (std::vector<bead*>::iterator it = Sysbead.begin(); it != Sysbead.end(); ++it) {
        Vec3D beadPosition((*it)->GetXPos(), (*it)->GetYPos(), (*it)->GetZPos());

        for (int i = 0; i < 3; ++i) {
            int nx = static_cast<int>(beadPosition(i) / (*pBox)(i));

            if (beadPosition(i) >= 0)
                beadPosition(i) -= nx * (*pBox)(i);
            else
                beadPosition(i) += (nx + 1) * (*pBox)(i);
        }

        (*it)->UpdatePos(beadPosition(0), beadPosition(1), beadPosition(2));
    }
}
std::vector<bead> Solvate::AddIons(std::vector<bead>& FullWaterBead, int Nposion, int Nnegion, const std::string& pName, const std::string& nName, int seed) {
    std::vector<bead> outBeads;  // Final output bead vector containing water and ions
    const int numer_of_water = FullWaterBead.size();
    const int numer_of_total_ions = Nposion + Nnegion;

    if (numer_of_total_ions > numer_of_water) {
        std::cout << "---> error: total number of requested ions is larger than the total generated water beads\n";
        std::cout << "   ---> total requested ions " << numer_of_total_ions << "\n";
        std::cout << "   ---> total requested water beads " << numer_of_water << "\n";
        exit(0);
    }

    // shuffle the water beads to select ions without being localaized
    std::srand(seed);
    std::random_shuffle(FullWaterBead.begin(), FullWaterBead.end());

    int nn = 0;
    int np = 0;
    int id = 1;

    for (std::vector<bead>::iterator it = FullWaterBead.begin(); it != FullWaterBead.end(); ++it) {
        if (id > numer_of_water - numer_of_total_ions && id <= numer_of_water - numer_of_total_ions + Nposion) {
            it->UpdateBeadName(pName);
            it->UpdateResName("ION");
            ++np;
        } else if (id > numer_of_water - numer_of_total_ions) {
            it->UpdateBeadName(nName);
            it->UpdateResName("ION");
            ++nn;
        }
        ++id;
    }

    // just to check that the number of requested is equal to the generated one
    std::cout << "---> created ions " << np << " positive  " << nn << " negative ions\n";

    // Store all the beads in a new container to maintain the order.
    outBeads = FullWaterBead;

    // Report some info about numbers
    std::ofstream info("info.txt");
    if (info.is_open()) {
        info << "W    " << numer_of_water - numer_of_total_ions << "\n";
        if (Nposion != 0)
            info << pName << "    " << Nposion << "\n";
        if (Nnegion != 0)
            info << nName << "    " << Nnegion << "\n";
        info.close();
    }

    // Display info
    std::cout << "-------------------------------------------------------\n";
    std::cout << "-------------- generated ion and solvent --------------\n";
    std::cout << "W    " << numer_of_water - numer_of_total_ions << "\n";
    if (Nposion != 0)
        std::cout << pName << "    " << Nposion << "\n";
    if (Nnegion != 0)
        std::cout << nName << "    " << Nnegion << "\n";
    std::cout << "-------------------------------------------------------\n";

    return outBeads;
}

