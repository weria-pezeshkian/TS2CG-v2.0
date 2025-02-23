  #if !defined(AFX_BackMap_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_BackMap_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "BackMap.h"
#include "GroFile.h"
#include "GenerateUnitCells.h"
#include "Def.h"
#include "PDBFile.h"
#include "PointBasedBlueprint.h"
/*
 1) read the point
 2) exclude the point base of exclsuion data
 3) place proteins or generate them
 4) exclude point based on proteins
 5) place the lipids
 
 
 further extentions
 2.) angle and theta is not used yet (not needed), you can rotate it in the gro file.
 3.) pattern based protein insertion
 
 
 
 */
BackMap::BackMap(Argument *pArgu)
{
    m_monolayer = false;  // this is false
    m_Warning=0;
    srand (pArgu->GetSeed());
    std::cout<<"\n";
    std::cout<<"███████████████████████████████████████████████████████████████  \n";
    std::cout<<"                      "<< SoftWareName <<"              "<<"\n";
    std::cout<<"                         Version:  "<<SoftWareVersion<<"  \n";
    std::cout<<"███████████████████████████████████████████████████████████████  \n";
    Nfunction f;      // In this class there are some useful function and we can use it.

    //====== getting data points to create cg membrane
    std::cout<<"---> attempting to obtain point data \n";
    PointBasedBlueprint SurfDataPoint(pArgu);
    std::vector<point*>  pPointUp = SurfDataPoint.m_pPointUp;
    std::vector<point*>  pPointDown = SurfDataPoint.m_pPointDown;
    std::vector<inclusion*>  pInc = SurfDataPoint.m_pInc;
    std::vector<exclusion*>  pExc = SurfDataPoint.m_pExc;
    Vec3D *pBox = SurfDataPoint.m_pBox;
    m_pBox = pBox;
    Wall *pWall = SurfDataPoint.m_pWall;
    m_monolayer = SurfDataPoint.m_monolayer; // if pPointDown is empty, m_monolayer is false,
    //if m_monolayer is not true, we check the PCG option and see if monolayer is defined
    if(m_monolayer==false && pArgu->GetMonolayer()==true)
    {
        pPointDown.clear();
        m_monolayer=true;
        std::cout<<"---> note: while we have points to build both monolayers, since monolayer is defined in PCG commond, we create only one layer \n";
    }
    std::cout<<"---> point data has been obtained \n";

    //======== OutPut file name declaration and finding input file names ========================================
    std::string gname = pArgu->GetGeneralOutputFilename();   // get the generic name for the outputs
    m_FinalOutputGroFileName =gname+".gro";                     // create the output gro file name based on the generic name
    m_FinalTopologyFileName=gname+".top";
    
    //======
    std::cout<<"---> attempting to generate molecule type \n";
    GenerateMolType  MOLTYPE(pArgu);   // using the str file, the included gro file in the str and the lib file, different mol types will be generated. proteins and lipids are treated as mols.
    m_map_MolName2MoleculesType = MOLTYPE.GetMolType();   // a map containing all the mol types: (name, MolType)
    std::cout<<"---> molecule types have been generated \n";
  
    //== we should exclude points and get rid of exclusion. This is done by making the area of the point zero.
    //== this could be made more efficient but is not needed as exclusion should not inlcude to many points
    ExcludePointsUsingExclusion(pExc, pPointUp, pPointDown);
    
    //==== now we need to place the proteins
    // first we read str file to find protein info
    std::string strfilename = pArgu->GetStructureFileName();    // str file name
    if(FindProteinList(strfilename)==false) // this data will be stored in m_map_IncID2ProteinLists map(tsi_protein_id,ProteinList)
        std::exit(0);
    //== before going further, it would be better to check if the info about
    //proteins in str file also exist in the PLM output and also in the mol type
    if(CheckProteinInfo (m_map_IncID2ProteinLists, m_map_MolName2MoleculesType, pInc)==false) //
        std::exit(0);
    

    std::vector<inclusion> RandomInc;
    if(pInc.size()!=0)
    {
        //== I do not know why do we need this
        for ( std::map<int,ProteinList>::iterator it = m_map_IncID2ProteinLists.begin(); it != m_map_IncID2ProteinLists.end(); it++ )
        {
            int number = 0;
            int id = (it->second).ID;
            for ( std::vector<inclusion*>::iterator it1 = pInc.begin(); it1 != pInc.end(); it1++ )
            {
                if(id==(*it1)->GetTypeID())
                    number++;
            }
            (it->second).created = number;
        }
    }
    else  // if the inclusion file is empty, we generate proteins according to the data in the str file
    {
        RandomInc = CreateRandomInclusion(pPointUp, pBox);
        for (std::vector<inclusion>::iterator it2 = RandomInc.begin() ; it2 != RandomInc.end(); ++it2)
            pInc.push_back(&(*it2));
    }
    
//== better to place the proteins and make then remove the excluded points  then generate the domains. Note, we should tell the user some domain could be removed due to protein and exclsuion
    //== however, if a error happen due to user, the error will be revealed later about domain. We can add a check function later
    m_ResID = 1; // we set the resid of the first mol to 2, the number will be updated when a mol is generated
    m_InclusionDirectionType = pArgu->GetInclusionDirectionType(); //Note: in the normal condition PLM always write global so only applicable if you want to change the point folder manually
    if(m_InclusionDirectionType=="Local")
    {
        std::cout<<"The inclusion direction type is: "<<m_InclusionDirectionType<<"\n";
        std::cout<<" Note: in the normal condition PLM always write global so only applicable if you want to change the point folder manually\n";
    }
    //== Placing the inclusions;
    if(PlaceProteins(pPointUp,pInc)==false) // this creates all the protein beads and put them in m_FinalBeads
        std::exit(0);
    std::cout<<"---> proteins are placed, now we remove points that are close to the proteins \n";
    
    //=== removing points closeby the proteins
    if(!pArgu->m_KEEP_POINTS_CLOSE_TO_PROTEINS)
    {
        double RCutOff = pArgu->GetRCutOff();     /// That will be counted as a cutoff for the protein_lipid distance
        std::vector<bead*> tempropbeads;
    
        for ( std::vector<bead>::iterator it = m_FinalBeads.begin(); it != m_FinalBeads.end(); it++ )
            tempropbeads.push_back(&(*it));
        
        bool RMpoint = RemovePointsCloseToBeadList(pPointUp, pPointDown, tempropbeads, RCutOff, pBox);
        
        if(RMpoint==false){
            std::exit(0);
        }
    }
    
    std::vector<Domain*> pAllDomain;
    bool Renormalizedlipidratio = pArgu->GetRenorm();
    m_Iter = pArgu->GetIter();  // how many iteration should be made to make sure enough lipid is placed.
    GenDomains GENDOMAIN(strfilename,pPointUp,pPointDown,Renormalizedlipidratio);  // this somehow reads the lipids
    pAllDomain = GENDOMAIN.GetDomains();

if(pArgu->Skip_LipidPlacement() == false){
    std::cout<<"---> generating domains using the input files \n";



 //   bool Renormalizedlipidratio = pArgu->GetRenorm();
  //  m_Iter = pArgu->GetIter();  // how many iteration should be made to make sure enough lipid is placed.
   // GenDomains GENDOMAIN(strfilename,pPointUp,pPointDown,Renormalizedlipidratio);  // this somehow reads the lipids
   // pAllDomain = GENDOMAIN.GetDomains();
    //
    std::cout<<"---> now,  the domain info will be used to place lipids \n";

    std::cout<<"expected time:  ";
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ ) // all the domains
    std::cout<<"█";
    
    std::cout<<"\n";
    std::cout<<"remaining time: ";
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ ) // all the domains
    {
        if((*it)->GetDomainPoint().size()!=0)// this is only valid if
            GenLipidsForADomain(*it);
        std::cout<<"█";
    }
    std::cout<<"\n";

    std::string sms = InfoDomain(pAllDomain);
    std::cout<<sms;
}//    if(pArgu->Skip_LipidPlacement() == false){
else{
    std::cout<<"--> Note: We have skipped lipid placement as requested. \n";
}
    //=============== write the wall info
    std::cout<<"---> attempting to make the wall beads \n";
    std::vector<bead> WB = pWall->GetWallBead();
    if((pWall->GetWallPoint()).size()>0 && pWall->GetState()==true)
    {
        PDBFile pdb;
        std::string pdbfile = "Wall.pdb";
        pdb.WritePDBFile(pdbfile, pWall->GetWallPoint());
    }
    else if(pWall->GetState()==true)
    {
        std::cout<<"Note ----> No wall.pdb file will be generated since the total created wall beads are zero \n";
    }
    for (std::vector<bead>::iterator it = WB.begin() ; it != WB.end(); ++it)
    {
        m_FinalBeads.push_back((*it));
    }
    std::cout<<"---> attempting to write the final gro file \n";
    WriteFinalGroFile(pBox);
    std::cout<<"---> attempting to write the final topology file \n";
    bool gentop = GenTopologyFile(pAllDomain,(WB.size()));
    if(m_Warning==0)
    {
        std::cout<<" ██████████████████████████████████████████████████████████████  \n";
        std::cout<<" █████████  Seems everything went well. Well done! ████████████  \n";
        std::cout<<" ██████████████████████████████████████████████████████████████  \n";
    }
    else
    {
        std::cout<<" █████████  outputs have been generated, but there were "<< m_Warning<<" warnings in the process █████████████  \n";
    }
    Welldone();
    
    


}
BackMap::~BackMap()
{
    
}
void BackMap::GenLipid(MolType moltype, int listid, Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2)
{
     Tensor2 LG = TransferMatLG(Normal, t1, t2);
      std::vector<bead> vbeads = moltype.Beads;
            for ( std::vector<bead>::iterator it = vbeads.begin(); it != vbeads.end(); it++ )
            {
                Vec3D BPos((*it).GetXPos(),(*it).GetYPos(),(*it).GetZPos());
                Vec3D vX = LG*BPos+ Pos;
                int beadid = (m_FinalBeads.size()+1);
                bead TemB(beadid, (*it).GetBeadName(), (*it).GetBeadType(), (*it).GetResName(), m_ResID, vX(0), vX(1),vX(2));
                TemB.BringBeadInBox(m_pBox);
                m_FinalBeads.push_back(TemB);
            }
    
    m_ResID++;
    
    return;
}
void BackMap::GenProtein(MolType moltype, int listid, Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2)
{
        Tensor2 LG = TransferMatLG(Normal, t1, t2);
        Tensor2 GL = LG.Transpose(LG);

        //===== to fit to the protein diretion
        Vec3D LocalDir;
        if(m_InclusionDirectionType=="Global") //Note: in the normal condition PLM always write global so only applicable if you want to change the point folder manually
        LocalDir = GL*Dir;
        else if(m_InclusionDirectionType=="Local")
        LocalDir = Dir;
        double C=LocalDir(0);
        double S= LocalDir(1);
        Tensor2 Rot=Rz(C,S);
        Vec3D DH= Normal*((m_map_IncID2ProteinLists.at(listid)).Z0);
        double phi = (m_map_IncID2ProteinLists.at(listid)).Phi;
        double theta = (m_map_IncID2ProteinLists.at(listid)).Theta;

        //
        std::vector<bead> vbeads = moltype.Beads;
        for ( std::vector<bead>::iterator it = vbeads.begin(); it != vbeads.end(); it++ )
        {
            Vec3D BPos((*it).GetXPos(),(*it).GetYPos(),(*it).GetZPos());
            Vec3D vX = LG*(Rot*BPos)+ Pos+DH;
            int beadid = (m_FinalBeads.size()+1);
            bead TemB(beadid, (*it).GetBeadName(), (*it).GetBeadType(), (*it).GetResName(), m_ResID, vX(0), vX(1),vX(2));
            TemB.BringBeadInBox(m_pBox);
            m_FinalBeads.push_back(TemB);
        }
        m_ResID++;
    return;
}
void BackMap::WriteFinalGroFile(Vec3D *pBox)
{
    FILE *fgro;
    fgro = fopen(m_FinalOutputGroFileName.c_str(), "w");
    /// resid  res name   noatom   x   y   z
    const char* Title=" System ";
    int Size=m_FinalBeads.size();
    
    fprintf(fgro,  "%s\n",Title);
    fprintf(fgro, "%5d\n",Size);
    int i=0;
    for (std::vector<bead>::iterator it = m_FinalBeads.begin() ; it != m_FinalBeads.end(); ++it)
    {
        i++;
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos();
        std::string rname  = (*it).GetResName();
        std::string bname = (*it).GetBeadName();
        
        const char* A1=rname.c_str();
        const char* A2=bname.c_str();
        int resid=(*it).GetResid();
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid%100000,A1,A2,i%100000,x,y,z );
        
    }

    fprintf(fgro,  "%10.5f%10.5f%10.5f\n",(*pBox)(0),(*pBox)(1),(*pBox)(2) );
    fclose(fgro);
    
    return;
}
Tensor2 BackMap::Rz(double cos, double sin)
{
    Tensor2  R('O');
    R(0,0) = cos;
    R(0,1) = -sin;
    R(1,0) = sin;
    R(1,1) = cos;
    R(2,2) = 1;

    return R;
}
//=== we can make this function better
std::vector<inclusion> BackMap::CreateRandomInclusion(std::vector<point*> &pPointUp, Vec3D *pBox)
{
    
    std::vector<inclusion>  RandomInc;
    std::vector<ExcludedVolumeBeads>  ExcludeBeads;

    int totinc = 0;    // total number of inclusions of all types to be created
    int totcreated = 0;
    //==== reading input file to find number of requested proteins
    if(m_map_IncID2ProteinLists.size()!=0)
    {
        std::cout<<"---> According to the data and area we generate  "<<m_map_IncID2ProteinLists.size()<<" proteins types \n";
    }
//===== find the total number demanded of each proteins and how many proteins is assked to create
    for ( std::map<int,ProteinList>::iterator it = m_map_IncID2ProteinLists.begin(); it != m_map_IncID2ProteinLists.end(); it++ )
    {
        (it->second).created = 0;
        double ratio = (it->second).Ratio;
        std::string type = (it->second).ProteinName;
//======================== error message =================================
        if (m_map_MolName2MoleculesType.count(type) == 0)
        {
            std::cout<<"---> error: mol type "<<type<<" does not exist in the map \n";
            std::cout<<"    ----------------------<error>----------------------------- \n";
            std::cout<<"         below list exist in the lib and gro \n";
            for ( std::map<std::string , MolType>::iterator it = m_map_MolName2MoleculesType.begin(); it != m_map_MolName2MoleculesType.end(); it++ )
            {
                std::cout<<"                      -) "<<(*it).first<<"\n";
            }
            std::cout<<"    ---------------------------------------------------------- \n";
            exit(0);
        }
//=========================================================
        double membrane_total_area = 0 ;
        for ( std::vector<point*>::iterator it = pPointUp.begin(); it != pPointUp.end(); it++ )
            membrane_total_area+= (*it)->GetArea();
        
        double area = (m_map_MolName2MoleculesType.at(type)).molarea;
        int neededno = membrane_total_area/area*ratio;   // A/a_pro=number of possible proteins. multipled by the ratio give the demanded number
        (it->second).Maxno = neededno;                     // updated the max no of inc
        totinc+=neededno;                                   // too see how many inc is needed.
        
        std::cout<<"---> We are attempting to generate  "<<neededno<<" proteins of type "<<type<<"  \n" ;

    }
//==== end of evalutaing data in the str file to find number of proteins

    int id=0;
    int s=0;
    //=== checking and finding random position for proteins without overlapping. overlapping between proteins not with lipids
    while (totcreated<totinc && s<(pPointUp.size())) // to select a point randomly
    {
        s++;
        bool accept = true;
        //=== randomly selecting a point
        point* temPoint = pPointUp.at((rand()%pPointUp.size()));
        int pointid = temPoint->GetID();  // id of the randomly selected point
        int tid=0;                          // type id of a radon protein type that will be selected later
        double Rpro = 0;                    // radii of the chosen protein type

        ProteinList* pRandomProList;
        MolType* pRandomProtype;
        //== we find a protein randomly
        while (true)
        {
            std::map<int,ProteinList>::iterator it = m_map_IncID2ProteinLists.begin();
            std::advance(it, std::rand() % m_map_IncID2ProteinLists.size());
            
            pRandomProList = &(it->second);
            if(pRandomProList->created<pRandomProList->Maxno)
            {
                tid = pRandomProList->ID;
                //== this take the protein in the str and find a protype from the lib or gro
                pRandomProtype = &(m_map_MolName2MoleculesType.at(pRandomProList->ProteinName));// find the pro type asscoaited with the name
                double area = pRandomProtype->molarea;
                Rpro=sqrt(area/acos(-1));
                break;
                
            }
        }

//==== if there is an overlap between this protein and others that have been placed
        Vec3D XR1 =temPoint->GetPos();
        for ( std::vector<ExcludedVolumeBeads>::iterator it = ExcludeBeads.begin(); it != ExcludeBeads.end(); it++ )
        {
            double dist = it->R + Rpro; // R1 (it->R); R1+ Rpro is the minimum distance between two proteins
            Vec3D dR =it->X-XR1;   // dx between the chosen random point (a point to place protein) and the other selelected points
            if(dR.dot(dR,dR)<dist*dist)  // R1+R will be the minimum distance between two proteins
                accept = false;
            
            if(accept == true ) // check if in PBC the problems happens
            {
                for(int i=0;i<3;i++)
                    if(fabs((*pBox)(i))-fabs(dR(i))<fabs(dR(i)))
                        dR(i) = fabs((*pBox)(i))-fabs(dR(i));
                if(dR.dot(dR,dR)<dist*dist)  // R1+R will be the minimum distance between two proteins
                    accept = false;
            }
        }
        
// create the inclusion in the position of the selected point if everything was fine
        if(accept==true)
        {
            double d1= double(rand()%1000)/1000;
            double d2= double(rand()%1000)/1000;
            //
            Vec3D D(d1,d2,0);
            D=D*(1/(D.norm()));

            Vec3D N = temPoint->GetNormal();
            Vec3D T1 =   temPoint->GetP1();
            Vec3D T2 =   temPoint->GetP2();
            Tensor2 LG = TransferMatLG(N,T1,T2);
            D=LG*D;
            //
            id++;
            inclusion inc(id, tid,pointid,D);
            RandomInc.push_back(inc);
            totcreated++;
            pRandomProList->created = pRandomProList->created +1;
            
            
            ExcludedVolumeBeads Ex;
            Ex.X =temPoint->GetPos();
            Ex.R = Rpro;
            ExcludeBeads.push_back(Ex);
        }
    }
    
    return RandomInc;
}
bool BackMap::FindProteinList(std::string filename)
{
    //===  this function was updated in Dec, 2023.
    std::cout<<"---> reading protein information from the str file "<<"\n";
    Nfunction f;
    std::ifstream strfile;
    strfile.open(filename.c_str());
    if(strfile.good()==false)
    {
        std::cout<<"---> error (3R22-D): while opening the str file with  name: "<<filename<<"  ; check if the file exist "<<"\n";
        exit(0);
    }
    //== Read the file until reaching the protein section
    bool proteinflag=false;
    std::string str;
    while (true)
    {
        std::getline (strfile,str);
        if(strfile.eof()){
            break;
        }
        str = f.trim(str);
        str.erase(std::remove(str.begin(), str.end(), '['), str.end());
        str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
        str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
        if(str=="ProteinList")
        {
            proteinflag = true;
            break;
        }
    }
    //== now read the protein information
    if(proteinflag==true)
    {
      while (true)
      {
        std::getline (strfile,str);
        std::vector<std::string> Line = f.split(str);
        if(Line.size()!=0)
        {
            if(Line[0]=="End" || Line[0]=="END")
            {
                break;
            }
            else if(Line.size()>=6)
            {
                ProteinList CP;
                CP.ProteinName = Line.at(0);
                CP.ID = f.String_to_Int(Line.at(1));
                CP.Ratio = f.String_to_Double(Line.at(2));
                CP.Phi=f.String_to_Double(Line.at(3));
                CP.Theta=f.String_to_Double(Line.at(4));
                CP.Z0=f.String_to_Double(Line.at(5));
                m_map_IncID2ProteinLists.insert(std::pair<int,ProteinList>(CP.ID, CP));
            }
            else
            {
            std::cout<<"---> error (3R22-E): protein information in the str file is incomplete "<<"\n";
            std::cout<<"---> other possible causs of this error. 1) the protein section in str file does not end with word End "<<"\n";
            exit(0);
            }
        }
      }//while (true)
    }
    else
    {
    std::cout<<"---> warning: the str file does not contain any information about proteins "<<"\n";
        m_Warning++;
    }
    strfile.close();

    
    if(m_map_IncID2ProteinLists.size()!=0)
    {
      std::cout<<"\n";
      std::cout<<"          |------------------------------------------------------------|\n";
      std::cout<<"          |     Protein List and ID have been read from the input file |\n";
      std::cout<<"          |------------------------------------------------------------| \n";
      for ( std::map<int,ProteinList>::iterator it = m_map_IncID2ProteinLists.begin(); it != m_map_IncID2ProteinLists.end(); it++ )
      {
          std::cout <<"          |        inclusion with id "<< it->first  <<" is mapped to ------> "<< (it->second).ProteinName<< std::endl ;
      }
      std::cout<<"          |------------------------------------------------------------| \n\n";
    }
    return true;
}
Tensor2  BackMap::TransferMatLG(Vec3D Normal, Vec3D t1, Vec3D t2)
{
    Tensor2  GL(t1,t2,Normal);
    Tensor2 LG=GL.Transpose(GL);
    return LG;
    
}

//=== since 2024
// a function that use up the exclsuions by making the area of that specific points zero.
void BackMap::ExcludePointsUsingExclusion(std::vector<exclusion*> &pExc, std::vector<point*> &m_pPointUp, std::vector<point*> &m_pPointDown)
{
    
    if(pExc.size()!=0)
    {
        std::cout<<" Note: we are excluding points based on exclusion, If it is slow, contact the developer \n";
        
        for ( std::vector<exclusion*>::iterator it = pExc.begin(); it != pExc.end(); it++ )
        {
            int pointid=(*it)->GetPointID();
            if(pointid<0 || pointid>m_pPointUp.size()){
                std::cout<<"---> error: id = PCG23456: please report to the developer with the error id name \n";
                exit(-1);
            }
            
            point *Up_p1=m_pPointUp.at(pointid);
            Vec3D Pos = Up_p1->GetPos();
            Vec3D N = Up_p1->GetNormal();
            double R = (*it)->GetRadius();
            if(R!=0)
            Up_p1->UpdateArea(0);

            for ( std::vector<point*>::iterator it1 = m_pPointUp.begin(); it1 != m_pPointUp.end(); it1++ )
            {
                Vec3D Pos1 = (*it1)->GetPos();
                Vec3D DP = Pos1-Pos;
                double dist = DP.norm();
                Vec3D UnitDP =DP*(1/dist);
                double sinT = fabs((UnitDP*N).norm());
                double cosT = fabs(N.dot(UnitDP,N));

                if(dist*sinT<=R && dist*cosT<6)
                {
                    (*it1)->UpdateArea(0);

                }
            }
            for ( std::vector<point*>::iterator it1 = m_pPointDown.begin(); it1 != m_pPointDown.end(); it1++ )
            {
                Vec3D Pos1 = (*it1)->GetPos();
                Vec3D DP = Pos1-Pos;
                double dist = DP.norm();
                Vec3D UnitDP =DP*(1/dist);
                double sinT = fabs((UnitDP*N).norm());
                double cosT = fabs(N.dot(UnitDP,N));

                if(dist*sinT<=R && dist*cosT< 6)
                {
                    (*it1)->UpdateArea(0);
                    
                }
            }
        }
    }
}
//=== check if all the incs have been mapped to a proteinlist and the protein list exist in the moltype
bool BackMap::CheckProteinInfo (std::map<int , ProteinList>& plist, std::map<std::string , MolType>& moltype, std::vector<inclusion*> & incs)
{
    
    //== first find all the inc type id and store them
    std::vector<int> inc_typeid;
    for ( std::vector<inclusion*>::iterator it = incs.begin(); it != incs.end(); it++ )
    {
        int id = (*it)->GetTypeID();
        if(count(inc_typeid.begin(), inc_typeid.end(), id)<=0)
            inc_typeid.push_back(id);
    }
    //== check if for each id in incs, there is one in proteinlist
    for ( std::vector<int>::iterator it = inc_typeid.begin(); it != inc_typeid.end(); it++ )
    {
        
        if(plist.find(*it)!=plist.end())
        {
            // now we get the coresponding protein name and check if in the mol def file (either lib or gro files) such a protein exit
            std::string pname = ((plist.find(*it))->second).ProteinName;
            if(moltype.find(pname)!=moltype.end())
            {
                // do we need to check anything more?? for now no
            }
            else
            {
                std::cout<<"error 2522024_a--> there is an inclusion with type id "<<*it<<" and protein name "<<pname<<" however the stracture does not exist, perhaps include a gro file of it \n";
                return false;
            }
        }
        else
        {
            std::cout<<"error 2522024_b--> there is an inclusion with type id "<<*it<<" however it is not mapped to a protein in the str file \n";
            return false;
        }
            
    }
    
    std::cout<<"---> protein information in the str file, lib/gro files and inclusion info matches well \n";
    return true;
}

bool BackMap::PlaceProteins(std::vector<point*> &pPointUp, std::vector<inclusion*>  &pInc)
{
    for ( std::map<int,ProteinList>::iterator it1 = m_map_IncID2ProteinLists.begin(); it1 != m_map_IncID2ProteinLists.end(); it1++ )
    {
        int plistid=it1->first;
    for ( std::vector<inclusion*>::iterator it = pInc.begin(); it != pInc.end(); it++ )
    {
        int id = (*it)->GetTypeID();
        if(plistid==id)
        {
            std::string ptype=(m_map_IncID2ProteinLists.at(id)).ProteinName;
            int pointid = (*it)->GetPointID();
            Vec3D  Dir =  (*it)->GetDirection();
            point *Up_p1=pPointUp.at(pointid);
            Vec3D N = Up_p1->GetNormal();
            Vec3D Pos = Up_p1->GetPos();
            Vec3D T1 =   Up_p1->GetP1();
            Vec3D T2 =   Up_p1->GetP2();

            if (m_map_MolName2MoleculesType.count(ptype) == 0)
            {
                std::cout << "---> error: molecule name " <<ptype<<" does not exist in the attached gro files \n";
                return false;
            }
            
            GenProtein(m_map_MolName2MoleculesType.at(ptype), id, Pos, N, Dir, T1,T2);
        }
    }
    plistid++;
    }
    
    return true;
}
bool BackMap::RemovePointsCloseToBeadList(std::vector<point*> &pPointUp, std::vector<point*> &pPointDown, std::vector<bead*> vpbeads, double RCutOff, Vec3D* pBox)
{
    
    GenerateUnitCells GCNT(vpbeads, pBox,RCutOff,1.0);
    GCNT.Generate();
    // Here, we try to remove the points that are covered by the proteins. We do it by setting the A=0
  
    for ( std::vector<point*>::iterator it = pPointUp.begin(); it != pPointUp.end(); it++ )
    {
        bool rem = false;
        Vec3D Pos1 = (*it)->GetPos();
        Vec3D N =   (*it)->GetNormal();
        Vec3D Pos2 = Pos1 - N*1.5;
        rem = GCNT.anythingaround(Pos1);  //121945
        if(rem==false)
            rem = GCNT.anythingaround(Pos2);
        if(rem==true)
            (*it)->UpdateArea(0);

    }
    if(m_monolayer == false)
    for ( std::vector<point*>::iterator it = pPointDown.begin(); it != pPointDown.end(); it++ )
    {
        bool rem = false;
        Vec3D Pos1 = (*it)->GetPos();
        Vec3D N =   (*it)->GetNormal();
        Vec3D Pos2 = Pos1 - N*1.5;
        rem = GCNT.anythingaround(Pos1);  //121945
        if(rem==false)
            rem = GCNT.anythingaround(Pos2);
        if(rem==true)
            (*it)->UpdateArea(0);
    }
    
    return true;
}
bool BackMap::GenLipidsForADomain(Domain *pdomain)
{

    std::vector<DomainLipid*> pdomainlipids = pdomain->GetpDomainLipids();
    std::vector<point*>  dpoint = pdomain->GetDomainPoint();
   // int CreateTotalLipids = pdomain->GetDomainTotalLipid();

int iteration = 0;
int NoMadeTotalLipid = 0;
int NoMadeLipid = 0 ;  // temperory because our strcture cannot increase it

for (std::vector<DomainLipid*>::iterator it = pdomainlipids.begin(); it != pdomainlipids.end(); it++ )
{ //=============
    while(true)
    {
        iteration++;
        if(iteration>m_Iter*(dpoint.size()))
        {
            std::cout<<"---> Warning: with "<< m_Iter <<" iterations, we could not place the expected number of the lipids \n";
            std::cout<<" if you are unhappy, increase the number of the iteration with option -iter, or regenerate the points \n";
            m_Warning++;
            break;
        }
        int pointid = rand()%(dpoint.size());   // selecting a random point
        point* Ran_point = dpoint.at(pointid);  // select a random point
        double area = Ran_point->GetArea();

        double prob=area/((*it)->Ap);   // Ap is the area per lipid for that specific lipid
        double rn = double(rand()%(1000000))/1000000.0;

        if(prob>rn && NoMadeLipid < (*it)->MaxNo )
        {
            NoMadeTotalLipid++;
            NoMadeLipid++;
            (*it)->no_created=(*it)->no_created+1;
            {//================================ Create a single lipid at the point position ======================================
                Vec3D  Dir(0,0,0);
                Vec3D N =    Ran_point->GetNormal();
                Vec3D T1 =   Ran_point->GetP1();
                Vec3D T2 =   Ran_point->GetP2();
                Vec3D Pos =  Ran_point->GetPos();
                std::string ltype = (*it)->Name;
                if (m_map_MolName2MoleculesType.count(ltype) == 0)
                {
                    std::cout << " \n---> error: molecule name " <<ltype<<" does not exist in the lib files \n";
                    exit(0);
                }
                GenLipid(m_map_MolName2MoleculesType.at(ltype), 0, Pos, N, Dir, T1, T2);
                Ran_point->UpdateArea(0);
            }//============================================================================================================
        }
        if((*it)->MaxNo==NoMadeLipid )
        {
            NoMadeLipid = 0;
            break;
        }
    } //while(true)
}
    
return true;
}
//=============== make topology file
bool BackMap::GenTopologyFile(std::vector<Domain*> pdomains, int WBead_no)
{
    //==========================================================================================================
    //==========================================================================================================
    std::ofstream Topgro;
    Topgro.open(m_FinalTopologyFileName.c_str());
    Topgro<<" ;This file was generated by TS2CG membrane builder script i.e., PCG \n";
    Topgro<<" [ system ] \n";
    Topgro<<" Expect a large membrane \n";
    Topgro<<" [ molecules ] \n";
    
    for ( std::map<int,ProteinList>::iterator it = m_map_IncID2ProteinLists.begin(); it != m_map_IncID2ProteinLists.end(); it++ )
        Topgro<<(it->second).ProteinName<<"   "<<(it->second).created<<"\n";

    int layer = 0;
    
    
    for ( std::vector<Domain*>::iterator it = pdomains.begin(); it != pdomains.end(); it++ )
    {
        layer++;
        std::vector<DomainLipid> DL = (*it)->GetDomainLipids();
        
        if(layer%2!=0 || m_monolayer == false)
        {
            Topgro <<"; domain "<<(*it)->GetDomainID() <<" \n";
            if(layer%2!=0 && m_monolayer == false)
            Topgro  <<" ;  in the upper monolayer \n";
            else if(layer%2==0 && m_monolayer == false)
            Topgro <<" ;  in the lower monolayer \n";
            for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
            {
                if((*it2).no_created == 0)
                    Topgro <<"   ;  ";
                Topgro <<"     "<<(*it2).Name<<"  "<<(*it2).no_created<<"     "<<std::endl ;
                
            }
        }

    }
    if(WBead_no!=0)
    Topgro<<"Wall    "<<WBead_no<<"\n";
    
 
    return true;
}
void BackMap::Welldone()
{
//std::cout << "\n ████████████████ Well Done ██████████████████\n";
}

std::string BackMap::InfoDomain(std::vector<Domain*> pAllDomain)
{
    Nfunction f;
    int layer = 0;
    std::string message;
    message+= "        |--------------------------------------------------| \n";
    message+= "        |        Information on the Generated Lipids       |\n";
    message+= "        |--------------------------------------------------| \n";
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {

        layer++;
        std::vector<DomainLipid> DL = (*it)->GetDomainLipids();
        
        if(layer%2!=0 || m_monolayer == false)// this means that if it is monolayer, we do not print the inner monolayer info
        {
                message+="        | -> For domain with ID "+f.Int_to_String((*it)->GetDomainID())+" \n";
            for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
            {
                if((*it2).MaxNo==(*it2).no_created)
                {
                message+="        |   -> PCG created "+f.Int_to_String((*it2).no_created)+"  molecules of "+(*it2).Name+" \n" ;
                }
                else
                {
                message+="        |      -> warning "+f.Int_to_String((*it2).MaxNo)+"  molecules of "+(*it2).Name+" was planned BUT we succeed to create  "+f.Int_to_String((*it2).no_created)+"  \n" ;
                m_Warning++;
                }
            }
        }
        if(layer%2!=0 && m_monolayer == false)
        {
                message+="        |     In the upper monolayer \n";
                message+="        |     \n";

        }
        else if(layer%2==0 && m_monolayer == false)
        {
                message+="        |     In the lower monolayer \n";
                message+="        |     \n";
        }
        if(layer%2!=0 && m_monolayer == true)
                message+="        |     In the  monolayer \n";
    }
    message+= "        |-------------------------------------------------- \n";
    return message;
}


#endif



