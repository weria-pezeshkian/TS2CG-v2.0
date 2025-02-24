

#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdlib>
#include "Edit_configuration.h"
#include "WriteFiles.h"
#include "help.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "VMDOutput.h"
#include "Surface_Mosaicing.h"
#include "In_OR_Out.h"
#include "Traj_XXX.h"


/*
 This has been updated in Sept 2023. This software is now both backmapping and also input generator.
 */
Edit_configuration::Edit_configuration( std::vector <std::string> Arguments) :
                    m_LessOutPut (false),           // print as you wish
                    m_smooth (false),               // Smoothing flag
                    m_Folder("point"),               // Default output folder name
                    m_TaskName("PLM"),             // Default task name
                    m_monolayer(0),                // Monolayer flag (0 for bilayer)
                    m_FindnewBox(false),          // Flag to find the smallest possible box
                    m_MosAlType("Type1"),          // Default algorithm type
                    m_BilayerThickness(3.8),       // Default bilayer thickness
                    m_Iteration(-1),               // Default degree of meshing (disabled)
                    m_MeshFileName("TS.q"),        // Default mesh file name
                    m_AP(0.62)                    // Default area per lipid
{
    // Initialize the Variables to their default values
    InitializeVariables();
    // updates variables based on the command line arguments
    UpdateVariables(Arguments);
    if(!ValidateVariable()){
        std::cerr << "---> error: bad inputs.\n";
        exit(0);
    }
    
    
    m_pCurvatureCalculations = &m_CurvatureCalculations;

//// do the jobs
  if(m_TaskName == "in_out")
  {
      In_OR_Out in_out(m_MeshFileName);
  }
  else if(m_TaskName=="check")
  {
        check(m_MeshFileName);
  }
  else if(m_TaskName=="minimize")
  {
        Minimize(m_MeshFileName);
  }
  else if(m_TaskName=="vertexinfo")
  {
          // a new feature in TS2CG2.0 and above. Finding geometric info on a vertex
          VertexInfo(m_MeshFileName);
  }
  else if (m_TaskName=="PLM")
  {
        double H=m_BilayerThickness/2.0;
        const int dir_err = system(("mkdir -p "+m_Folder).c_str());
        if (-1 == dir_err)
        {
            std::cout<<"error--> creating directory  "<<m_Folder<<"\n";
            exit(1);
        }
        if(!m_LessOutPut){
          const int dir_err2 = system(("mkdir -p "+m_Folder+"visualization_data").c_str());
          if (-1 == dir_err2)
          {
            std::cout<<"error--> creating directory  visualization_data"<<"\n";
            exit(1);
          }
        }
                BackMapOneLayer(1 , m_MeshFileName, H);
                if(m_monolayer==0)
                BackMapOneLayer(-1 , m_MeshFileName, H);
  }
  else {
      std::cout<<" error--> unrecognized Task \n";
  }
}
Edit_configuration::~Edit_configuration()
{
    
}
void Edit_configuration::Rescaling(Vec3D zoom , MESH *pMesh)
{
    Vec3D L = pMesh->m_Box;
    (*(pMesh->m_pBox))(0) = zoom(0)*L(0);
    (*(pMesh->m_pBox))(1) = zoom(1)*L(1);
    (*(pMesh->m_pBox))(2) = zoom(2)*L(2);
    for (std::vector<vertex *>::iterator it = (pMesh->m_pActiveV).begin() ; it != (pMesh->m_pActiveV).end(); ++it)
    {
        double x = (*it)->GetVXPos();
        double y = (*it)->GetVYPos();
        double z = (*it)->GetVZPos();
        (*it)->UpdateVXPos(zoom(0)*x);
        (*it)->UpdateVYPos(zoom(1)*y);
        (*it)->UpdateVZPos(zoom(2)*z);
        
    }
    return;
}
bool Edit_configuration::FileExist (const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}
void  Edit_configuration::UpdateGeometry(MESH *pmesh)
{
    Curvature CurvatureCalculations;
    for (std::vector<triangle *>::iterator it = (pmesh->m_pActiveT).begin() ; it != (pmesh->m_pActiveT).end(); ++it)
    {
    (*it)->UpdateNormal_Area(m_pBox);
    }
    
    //===== Prepare links:  normal vector and shape operator
    for (std::vector<links *>::iterator it = (pmesh->m_pHL).begin() ; it != (pmesh->m_pHL).end(); ++it)
    {
            (*it)->UpdateNormal();
            (*it)->UpdateShapeOperator(m_pBox);
    }
    //======= Prepare vertex:  area and normal vector and curvature of surface vertices not the edge one
    for (std::vector<vertex *>::iterator it = (pmesh->m_pSurfV).begin() ; it != (pmesh->m_pSurfV).end(); ++it)
        CurvatureCalculations.SurfVertexCurvature(*it);
    //====== edge links should be updated
    for (std::vector<links *>::iterator it = (pmesh->m_pEdgeL).begin() ; it != (pmesh->m_pEdgeL).end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);
    for (std::vector<vertex *>::iterator it = (pmesh->m_pEdgeV).begin() ; it != (pmesh->m_pEdgeV).end(); ++it)
        CurvatureCalculations.EdgeVertexCurvature(*it);
}
bool Edit_configuration::check(std::string file){

    // generating the mesh
    std::string domyfile;
    CreateMashBluePrint BluePrint;
    MeshBluePrint meshblueprint;
    meshblueprint = BluePrint.MashBluePrintFromInput_Top(domyfile,m_MeshFileName);
    m_Mesh.GenerateMesh(meshblueprint);
    m_pMesh = &m_Mesh;
    m_pBox=m_pMesh->m_pBox;
    
    UpdateGeometry(m_pMesh);

        double Tarea  = 0;
        double totalgaussianC = 0;
        double pi = acos(-1);
        for (std::vector<vertex *>::iterator it = (m_pMesh->m_pSurfV).begin() ; it != (m_pMesh->m_pSurfV).end(); ++it)
        {
            
            Vec3D X((*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
            Vec3D normal = (*it)->GetNormalVector();
            Tarea+= (*it)->GetArea();
            std::vector <double> C = (*it)->GetCurvature();
            totalgaussianC+= (C.at(0))*(C.at(1))*((*it)->GetArea());
        }
        int nv = (m_pMesh->m_pSurfV).size();
        int nt = (m_pMesh->m_pActiveT).size();
        int nl = (m_pMesh->m_pHL).size();
        std::cout<<" ************************************************************************************** \n";
        std::cout<<" Total area of the surface is "<<Tarea<<" nm^2 you can use -rescalefactor to increase it \n";
        //std::cout<<" total gaussian curvature /2PI is: "<<totalgaussianC/(2*pi)<<" \n";
       // std::cout<<" total mean curvature /2PI is: "<<totalmeanC/(2*pi)<<" \n";
        std::cout<<" Euler characteristic: "<<nv+nt-nl<<" \n";
        std::cout<<" topology genus "<<(2-(nv+nt-nl))/2<<" \n";


    return true;
}
void Edit_configuration::VertexInfo(std::string file){
    // generating the mesh
    std::string domyfile;
    CreateMashBluePrint BluePrint;
    MeshBluePrint meshblueprint;
    meshblueprint = BluePrint.MashBluePrintFromInput_Top(domyfile,m_MeshFileName);
    m_Mesh.GenerateMesh(meshblueprint);
    m_pMesh = &m_Mesh;
    m_pBox=m_pMesh->m_pBox;
    
        UpdateGeometry(m_pMesh);
        int index;
        std::cout<<" Please provide the index of the vertex (note, it will start from zero): \n";
        std::cin>>index;
        
        if(index<0 || index>(m_pMesh->m_pActiveV).size())
        {
            std::cout<<" ---> error:the index number is invalid \n";
            exit(0);
        }
        vertex *pv = (m_pMesh->m_pActiveV).at(index);
        

        double x = pv->GetVXPos();
        double y = pv->GetVYPos();
        double z = pv->GetVZPos();
        double A = pv->GetArea();
        Tensor2  L2GT = pv->GetL2GTransferMatrix();
        Tensor2  G2LT = pv->GetG2LTransferMatrix();
        Vec3D normal  = pv->GetNormalVector();
        std::vector <double> C = pv->GetCurvature();
        
        std::cout<<" vertex position is "<<x<<"  "<<y<<"  "<<z<<" \n ";
        std::cout<<" vertex area is "<<A<<" \n ";
        std::cout<<" vertex normal is "<<normal(0)<<"  "<<normal(1)<<"  "<<normal(2)<<" \n ";

        Vec3D t1(1,0,0);
        Vec3D t2(0,1,0);
        Vec3D n(0,0,1);
        
        Vec3D T1 = L2GT*t1;
        Vec3D T2 = L2GT*t2;
        Vec3D N = L2GT*n;
        std::cout<<" first direction "<<T1(0)<<"  "<<T1(1)<<"  "<<T1(2)<<" \n ";
        std::cout<<" second direction "<<T2(0)<<"  "<<T2(1)<<"  "<<T2(2)<<" \n ";
        std::cout<<" normal "<<N(0)<<"  "<<N(1)<<"  "<<N(2)<<" \n ";
}


/// ====== Updates in 2025
void Edit_configuration::UpdateVariables(const std::vector<std::string>& Arguments) {
    Nfunction f;  // Utility object for pre-made functions
    std::vector<std::string> v_error;  // To store error messages

    // Check if help is called
    for (size_t i = 1; i < Arguments.size(); ++i) {
        if (Arguments[i] == Def_HelpCall) {
            help helpmessage(SoftWareVersion, Arguments[0]);
            exit(0);  // Exit without generating a log file
        }
    }

    // Open log file
    std::ofstream log("plm.log");
    if (!log.is_open()) {
        std::cerr << "---> error: Could not open log file.\n";
        return;
    }

    // Log the command used to generate outputs
    log << "-- Outputs have been generated by the following command:\n";
    for (size_t i = 0; i < Arguments.size(); ++i) {
        log << Arguments[i] << "  ";
    }
    log << "\n";

    // Process arguments
    for (size_t i = 1; i < Arguments.size(); i += 2) {
        try {
            if (Arguments[i] == Def_TSfile) {
                m_MeshFileName = Arguments[i + 1];  // TS file name
            } else if (Arguments[i] == Def_BilayerThickness) {
                m_BilayerThickness = f.String_to_Double(Arguments.at(i + 1));
            } else if (Arguments[i] == Def_AlgType) {
                m_MosAlType = Arguments.at(i + 1);  // Algorithm type
            } else if (Arguments[i] == Def_rescalefactor) {
                m_Zoom(0) = f.String_to_Double(Arguments.at(i + 1));
                m_Zoom(1) = f.String_to_Double(Arguments.at(i + 2));
                m_Zoom(2) = f.String_to_Double(Arguments.at(i + 3));
                i += 2;  // Skip next two arguments
            } else if (Arguments[i] == Def_AreaPerLipid) {
                m_AP = f.String_to_Double(Arguments[i + 1]);
            } else if (Arguments[i] == Def_DegreeOfMeshing) {
                m_Iteration = f.String_to_Int(Arguments.at(i + 1));
            } else if (Arguments[i] == Def_TaskName) {
                m_TaskName = Arguments[i + 1];
            } else if (Arguments[i] == Def_OutputFolderName) {
                m_Folder = Arguments[i + 1];
            } else if (Arguments[i] == Def_SmoothingFlag) {
                m_smooth = true;
                --i;  // No additional argument for this flag
            } else if (Arguments[i] == Def_resizebox) {
                m_FindnewBox = true;
                --i;  // No additional argument for this flag
            } else if (Arguments[i] == Def_Monolayer) {
                m_monolayer = f.String_to_Int(Arguments[i + 1]);
                m_BilayerThickness = 0;
            } else if (Arguments[i] == Def_PrintLessPutput) {
                    m_LessOutPut = true;
                    --i;  // No additional argument for this flag
            } else {
                std::string error = "---> error: Unrecognized argument < " + Arguments[i] + " >";
                v_error.push_back(error);
                std::cout << error << "\n";
                break;
            }
        } catch (const std::exception& e) {
            std::string error = "---> error processing argument " + Arguments[i] + ": " + e.what();
            v_error.push_back(error);
            std::cout << error << "\n";
            break;
        }
    }

    // Log errors and exit if any
    if (!v_error.empty()) {
        for (size_t i = 0; i < v_error.size(); ++i) {
            log << v_error[i] << "\n";
        }
        log.close();
        exit(0);
    }

    // Validate TS file existence
    if (!Nfunction::FileExist(m_MeshFileName)) {
        std::string error = "---> error: TS file " + m_MeshFileName + " does not exist in the folder.";
        std::cout << error << "\n";
        log << error << "\n";
        log.close();
        exit(0);
    }

    log.close();
}
bool Edit_configuration::ValidateVariable(){
    
    
    if (m_Iteration < 1 ) {
        std::cerr << "---> error: Iteration must be larger then 0.\n";
        return false;
    }
    if (m_Zoom(0) <= 0 || m_Zoom(1) <= 0 || m_Zoom(2) <= 0 ) {
        std::cerr << "---> error: rescaling should be larger than zero \n";
        return false;
    }
    if (!(m_monolayer == 0 || m_monolayer == 1 || m_monolayer == -1)) {
        std::cerr << "---> error: monolayer should be either 0, or 1, it is set to = "<<m_monolayer<<" \n";
        return false;
    }
    if (m_AP <= 0 ) {
        std::cerr << "---> error: area per molecules should be larger then zero. \n";
        return false;
    }
    if (m_BilayerThickness < 0 ) {
        std::cerr << "---> error: bilayer thickness should be positive. \n";
        return false;
    }
    else if(m_BilayerThickness == 0){
        std::cerr << "---> warnning: bilayer thickness is set to zero. \n";
    }
    else if(m_BilayerThickness > 8){
        std::cerr << "---> warnning: bilayer thickness is unrealistic. Do you know what are you doing? \n";
    }
    

    return true;
}
void Edit_configuration::InitializeVariables() {
    // Initialize variables with default values
    m_Zoom(0) = 1;                  // Default zoom factor in x-direction
    m_Zoom(1) = 1;                  // Default zoom factor in y-direction
    m_Zoom(2) = 1;                  // Default zoom factor in z-direction
}
double  Edit_configuration::PPBCM_Cluster(double Lx, std::vector <double> X)
{
    double cm=0;
    
    
    std::vector <double> X1,X2;
    for (std::vector<double>::iterator it = X.begin() ; it != X.end(); ++it)
    {
        if((*it)>Lx/2)
        {
            X2.push_back((*it));
        }
        else
        {
            X1.push_back((*it));
        }
        
    }
    double cm1=0;
    double cm2=0;
    for (std::vector<double>::iterator it = X1.begin() ; it != X1.end(); ++it)
    {
        cm1+=(*it)/X1.size();
        
    }
    for (std::vector<double>::iterator it = X2.begin() ; it != X2.end(); ++it)
    {
        cm2+=(*it)/X2.size();
        
    }
    
        bool boxout=false;
    
    if(fabs(cm2-cm1)>Lx/2)
        boxout=true;
    if(boxout==false)
    {
        cm=(cm2*X2.size()+cm1*X1.size())/X.size();
    }
    else
    {
        cm=(cm2*X2.size()+(cm1+Lx)*X1.size())/X.size();
    }

    return cm;
    
}
void Edit_configuration::UpdateBoxSize(MESH* pmesh)
{
    // Finding the center of geometry
    Vec3D CM(0,0,0);
    int totalvNo = (pmesh->m_pActiveV).size();
    for (std::vector<vertex *>::iterator it1 = (pmesh->m_pActiveV).begin() ; it1 != (pmesh->m_pActiveV).end(); ++it1)
    {
        Vec3D X ((*it1)->GetVXPos(),(*it1)->GetVYPos(),(*it1)->GetVZPos());
        CM = CM+X*(1.0/double(totalvNo));
    }
    // Finding the max distance from the CMG
    Vec3D MaxB(0,0,0);
    for (std::vector<vertex *>::iterator it = (pmesh->m_pActiveV).begin() ; it != (pmesh->m_pActiveV).end(); ++it)
    {

            double dx = (*it)->GetVXPos()-CM(0);
            double dy = (*it)->GetVYPos()-CM(1);
            double dz = (*it)->GetVZPos()-CM(2);
        
            if(MaxB(0)<fabs(dx))
                MaxB(0)=fabs(dx);
            if(MaxB(1)<fabs(dy))
                MaxB(1)=fabs(dy);
            if(MaxB(2)<fabs(dz))
                MaxB(2)=fabs(dz);

    }
    double H = m_BilayerThickness/2;
    Vec3D DB ((6+H)/m_Zoom(0),(6+H)/m_Zoom(1),(6+H)/m_Zoom(2));
    (*m_pBox) = MaxB*2+DB;
    for (std::vector<vertex *>::iterator it1 = (pmesh->m_pActiveV).begin() ; it1 != (pmesh->m_pActiveV).end(); ++it1)
    {
        Vec3D X ((*it1)->GetVXPos(),(*it1)->GetVYPos(),(*it1)->GetVZPos());
        (*it1)->UpdateVXPos((*it1)->GetVXPos()-CM(0)+((*m_pBox)(0))*0.5);
        (*it1)->UpdateVYPos((*it1)->GetVYPos()-CM(1)+((*m_pBox)(1))*0.5);
        (*it1)->UpdateVZPos((*it1)->GetVZPos()-CM(2)+((*m_pBox)(2))*0.5);

    }
    
}
//=== the backmapping function
//==================================
void Edit_configuration::BackMapOneLayer(int layer , std::string file, double H)
{
    
//----> generating the mesh
    CreateMashBluePrint BluePrint;
    MeshBluePrint meshblueprint;
    meshblueprint = BluePrint.MashBluePrintFromInput_Top(file,file);
    MESH Mesh;
    Mesh.GenerateMesh(meshblueprint);
    MESH *pMesh =&Mesh;
    // Mesh is ready
    m_pBox=pMesh->m_pBox;
    Rescaling(m_Zoom,pMesh);
    UpdateGeometry(pMesh);
    
    if(m_FindnewBox==true)
        UpdateBoxSize(pMesh);
    
    // this means that we need to find the optiomal iteration value
    double totalvNo = double((pMesh->m_pActiveV).size());

//----> if the number of the iteration is not specified, we find an optimal one
    int Iteration = 0;
    if(m_Iteration==-1 )
    {
            double Tarea  = 0;
            for (std::vector<vertex *>::iterator it = (pMesh->m_pActiveV).begin() ; it != (pMesh->m_pActiveV).end(); ++it)
            Tarea+= (*it)->GetArea();
            double requiredNo = Tarea/(m_AP*totalvNo);
            Iteration = int (log(requiredNo)/log(4))+1;
            if(Iteration<5)
            std::cout<<"We will increase the number of the available points by 4^"<<Iteration<<"\n";
            else
            std::cout<<"The requested rescaling requires  4^"<<Iteration<<" points, this may take some time to finish \n";
    }
    else
        Iteration = m_Iteration;
//----> Moving each vertex in the direction of the normal vector
    for (std::vector<vertex *>::iterator it = (pMesh->m_pActiveV).begin() ; it != (pMesh->m_pActiveV).end(); ++it)
    {
        double x = (*it)->GetVXPos();
        double y = (*it)->GetVYPos();
        double z = (*it)->GetVZPos();
        Vec3D normal = (*it)->GetNormalVector();
        x=x+layer*H*(normal(0));
        y=y+layer*H*(normal(1));
        z=z+layer*H*(normal(2));
        
        (*it)->UpdateVXPos(x);
        (*it)->UpdateVYPos(y);
        (*it)->UpdateVZPos(z);
    }
    UpdateGeometry(pMesh);
//-----------> increasing the number of points, i.e., vertices
    std::vector <Surface_Mosaicing> Vmos;
    for (int j=0;j<Iteration;j++)
    {
        Surface_Mosaicing  MOS(m_MosAlType,m_smooth);
        Vmos.push_back(MOS);
    }

    for (int j=0;j<Iteration;j++)
    {
        std::cout<<" Iteration number "<<j+1<<" total is "<<Iteration<<"\n";

        // Here will cause error when
        (Vmos.at(j)).PerformMosaicing(pMesh);
        pMesh = (Vmos.at(j)).m_pMesh;
    }
  
    m_pBox = pMesh->m_pBox;

    double Lx=(*m_pBox)(0);
    double Ly=(*m_pBox)(1);
    double Lz=(*m_pBox)(2);
   
    if(!m_LessOutPut) {
        Vec3D BoxSides=(*m_pBox);
        std::string filename;
        
        if(layer==1)
        filename=m_Folder+"visualization_data/Upper";
        if(layer==-1)
        filename=m_Folder+"visualization_data/Lower";
        VMDOutput GRO(BoxSides, pMesh->m_pActiveV , pMesh->m_pActiveL, filename);
        GRO.WriteGro();
        GRO.WriteGro2();
        WriteFiles vtu(m_pBox);
        std::string fi;
        if(layer==1)
        fi=m_Folder+"visualization_data/Upper.vtu";
        if(layer==-1)
        fi=m_Folder+"visualization_data/Lower.vtu";
        vtu.Writevtu(pMesh->m_pActiveV,pMesh->m_pActiveT,pMesh->m_pActiveL,fi);
        Traj_XXX TSI(m_pBox);
        TSI.WriteTSI(0,"extended.tsi",pMesh->m_pActiveV,pMesh->m_pActiveT,pMesh->m_pInclusion,pMesh->m_pExclusion);

    }
    //=============
    std::string     UFUpper = m_Folder+"/OuterBM.dat";
    std::string     UFInner = m_Folder+"/InnerBM.dat";

    FILE *BMFile1;
    if(layer==1)
    BMFile1 = fopen(UFUpper.c_str(), "w");
    FILE *BMFile2;
    if(layer==-1)
    BMFile2 = fopen(UFInner.c_str(), "w");

    const char* Cbox="Box";
  
    if(layer==1)
    fprintf(BMFile1,  "%s%12.3f%12.3f%12.3f\n",Cbox,Lx,Ly,Lz);
    if(layer==2)
    fprintf(BMFile2,  "%s%12.3f%12.3f%12.3f\n",Cbox,Lx,Ly,Lz);
    
    const char* STR1="< Point NoPoints";
    const char* STR2=">";
    
    int NoPoints=(pMesh->m_pActiveV).size();
    if(layer==1)
    fprintf(BMFile1,  "%s%10d%s\n",STR1,NoPoints,STR2);
    if(layer==-1)
    fprintf(BMFile2,  "%s%10d%s\n",STR1,NoPoints,STR2);
    const char* Cont="< id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2 vtype >";
    if(layer==1)
    fprintf(BMFile1,  "%s\n",Cont);
    if(layer==-1)
    fprintf(BMFile2,  "%s\n",Cont);
    
    {
        const char* lay="< Outer >";
        if(layer==1)
        fprintf(BMFile1,  "%s\n",lay);
    }
    {
        const char* lay="< Inner >";
        if(layer==-1)
        fprintf(BMFile2,  "%s\n",lay);
    }
    int i=0;
    double dr=1;
    if (m_monolayer==-1)
    dr=-1;
   for (std::vector<vertex *>::iterator it = (pMesh->m_pActiveV).begin() ; it != (pMesh->m_pActiveV).end(); ++it)
    {
        double area=(*it)->GetArea();
        Vec3D normal=(*it)->GetNormalVector();
        normal = normal*(layer)*(dr);
        //    inclusion* inc = (*it)->GetInclusion();
        Tensor2  L2G = (*it)->GetL2GTransferMatrix();
        Vec3D LD1(1,0,0);
        Vec3D GD1 = L2G*LD1;
        Vec3D LD2(0,1,0);
        Vec3D GD2 = L2G*LD2;
        double x=(*it)->GetVXPos();
        double y=(*it)->GetVYPos();
        double z=(*it)->GetVZPos();
        int domain = (*it)->GetDomainID();
        int vtype = (*it)->m_VertexType;
        double c1,c2;
        if((*it)->m_VertexType==0)
        {
            std::vector <double> curvature = (*it)->GetCurvature();
            c1 = curvature[0];
            c2 = curvature[1];
        }
        else
        {
            c1 = (*it)->m_Normal_Curvature;
            c2 = 0;
        }
        
        if(layer==1)
        {

            fprintf(BMFile1,"%10d%5d%10.3f%10.3f%10.3f%10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%10d\n",i,domain,area,x,y,z,normal(0),normal(1),normal(2),GD1(0),GD1(1),GD1(2),GD2(0),GD2(1),GD2(2),c1,c2,vtype);
        }
        if(layer==-1)
        {
            fprintf(BMFile2,"%10d%5d%10.3f%10.3f%10.3f%10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%10d\n",i,domain,area,x,y,z,normal(0),normal(1),normal(2),GD1(0),GD1(1),GD1(2),GD2(0),GD2(1),GD2(2),c1,c2,vtype);
        }
        
        i++;
        
    }
if(layer==1)
{
    
    
   if ((pMesh->m_pInclusion).size()!=0){
    FILE *IncFile;
    IncFile = fopen((m_Folder+"/IncData.dat").c_str(), "w");
    
    
    
    const char* CHAR1 ="< Inclusion NoInc   ";
    NoPoints=((pMesh->m_pInclusion)).size();
    const char* CHAR2 ="   >";
    
    fprintf(IncFile,  "%s%5d%s\n",CHAR1,NoPoints,CHAR2);
    
    const char* CHAR3 ="< id typeid pointid lx ly lz  >";
    fprintf(IncFile,  "%s\n",CHAR3);
    
    
    
    i=0;
      
    for (std::vector<inclusion *>::iterator it = (pMesh->m_pInclusion).begin() ; it != (pMesh->m_pInclusion).end(); ++it)
    {
        
        int intypeid = ((*it)->GetInclusionTypeID());
        vertex* ver = (*it)->Getvertex();
        Tensor2  L2G = ver->GetL2GTransferMatrix();
        Vec3D LD = (*it)->GetLDirection();
        Vec3D GD = L2G*LD;
        int verid=ver->GetVID();
        fprintf(IncFile,  "%12d%12d%12d%8.3f%8.3f%8.3f\n",i,intypeid,verid,GD(0),GD(1),GD(2));
        i++;
    }
   }
    //==== We write Exclusion data
    
    if((pMesh->m_pInclusion).size()!=0)
    {
    FILE *ExcFile;
    ExcFile = fopen((m_Folder+"/ExcData.dat").c_str(), "w");
    
    
    
    const char* CHAR1 ="< Exclusion NoExc   ";
    NoPoints=(pMesh->m_pExclusion).size();
    const char* CHAR2 ="   >";
    
    fprintf(ExcFile,  "%s%5d%s\n",CHAR1,NoPoints,CHAR2);
    
    const char* CHAR3 ="< id  pointid r >";
    fprintf(ExcFile,  "%s\n",CHAR3);
    
    
    
    i=0;
    
    for (std::vector<exclusion *>::iterator it = (pMesh->m_pExclusion).begin() ; it != (pMesh->m_pExclusion).end(); ++it)
    {
        
        vertex* ver = (*it)->Getvertex();
        int verid=ver->GetVID();
        double R=(*it)->GetRadius();
        fprintf(ExcFile,  "%5d%10d%8.3f\n",i,verid,R);
        i++;
    }
    }
}

}
void Edit_configuration::Minimize(std::string file){
    std::cout<<" error---> this function has been removed \n";
    exit(0);
}
