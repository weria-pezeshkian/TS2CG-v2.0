

#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdlib>
#include "Curvature.h"
#include "In_OR_Out.h"
#include "WriteFiles.h"
#include "help.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "VMDOutput.h"
#include "Surface_Mosaicing.h"


/*
 this class is takes a TS file, ask for a point and checks if the point is inside the TS file.
 */
In_OR_Out::In_OR_Out(std::string filename)
{
    
    std::cout<<" Some fun, is a point inside or outside the TS \n";


    if(FileExist(filename) == false)
    {
        std::cout<<"Error: The TS file does not exist \n";
    }
    else
    {

        Initialize(filename);   // read this TS file
        m_pBox=m_pMesh->m_pBox;

        /// make the box size good; Find box size function create one
        ///
        /// 2) create a function to put the triangles inside a CNT cells
        ///
        /// 3) One line from the point to the box edge should give you what you need (only checking a few boxes )
        ///
        ///
        
        while(true)
        {
            std::cout<<" Enter 3 numbers for X, Y and Z coordinate of your point \n";
            double x,y,z;
            std::cin>>x>>y>>z;
            int inside = 0;
            for (std::vector<triangle *>::iterator it = (m_pMesh->m_pActiveT).begin() ; it != (m_pMesh->m_pActiveT).end(); ++it)
            {
                (*it)->UpdateNormal_Area(m_pBox);
                Vec3D N = (*it)->GetAreaVector();
                vertex *v = (*it)->GetV1();
                Vec3D X0(v->GetVXPos(),v->GetVYPos(),v->GetVZPos());
                Vec3D X(x,y,z);
                double value = N.dot((X-X0),N);
                if(value<0)
                {
                    inside--;
                }
                else if(value>0)
                {
                    inside++;
                }
                else
                {
                    std::cout<<" Point is on the surface \n";
                    std::cout<<0<<"\n";
                    break;
                }



            }
            if(inside>0)
            {
                std::cout<<" Point is outside \n";
                std::cout<<inside<<"\n";
            }
            else if(inside<0)
            {
                std::cout<<" Point is inside \n";
                std::cout<<inside<<"\n";
            }
            else
            {
                std::cout<<" unexpected \n";
            }
            
            
            std::cout<<" ****************** New try ************** \n";

        }


        

    }
}

In_OR_Out::~In_OR_Out()
{
    
}
void  In_OR_Out::UpdateGeometry(MESH *pmesh)
{
    Curvature CurvatureCalculations;
    for (std::vector<triangle *>::iterator it = (pmesh->m_pActiveT).begin() ; it != (pmesh->m_pActiveT).end(); ++it)
    (*it)->UpdateNormal_Area(m_pBox);

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
bool In_OR_Out::FileExist (const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}
void In_OR_Out::Initialize(std::string file)
{
    // generating the mesh
    std::string domyfile;
    CreateMashBluePrint BluePrint;
    MeshBluePrint meshblueprint;
    meshblueprint = BluePrint.MashBluePrintFromInput_Top(domyfile,file);
    std::vector<double> x;
    //m_Mesh.GenerateMesh(meshblueprint,0,0,x);
    m_pMesh = &m_Mesh;
    
    
    
}
