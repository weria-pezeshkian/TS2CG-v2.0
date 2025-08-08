#if !defined(AFX_Edit_configuration_H_8F4B21B8_D13C_9321_TT23_124095086224__INCLUDED_)
#define AFX_Edit_configuration_H_8F4B21B8_D13C_9321_TT23_124095086224__INCLUDED_
/*
 *  Copyright Weria Pezeshkian (weria.pezeshkian@gmail.com), last update 2023.
 */
#include "SimDef.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "Curvature.h"
#include "exclusion.h"
#include "MESH.h"
#include "CreateMashBluePrint.h"

class Edit_configuration
{
public:
    
	Edit_configuration( std::vector <std::string> arg);
	 ~Edit_configuration();

private:

    Vec3D *m_pBox;
    std::vector<inclusion*>   m_pInc;
    std::vector<exclusion*>   m_pExc;

    int  m_Iteration;  // how many time we should increase the number of triangles each time, *4T, therefore N_T = N_T*4^m_Iteration
    Vec3D m_Zoom ;
    std::string m_Folder ;  // Folder to write output in
    double m_AP ;
    bool m_FindnewBox;
    bool m_smooth;
    int m_monolayer;
    bool m_LessOutPut;
    void Rescaling(Vec3D zoom , MESH *pMesh);   // rescale the position and the box based on a vector
    void UpdateGeometry(MESH *pmesh);  // updates curvature, area etc of each triangle, vertex etc
    std::string m_MosAlType;
    void BackMapOneLayer(int layer , std::string file, double);
    bool check(std::string file);     // a function to check how the ts file looklike and do nothing
    void VertexInfo(std::string file);     // gives info about a vertex 

    void Minimize(std::string file);   // may not work well, needs optimization
    
    double  PPBCM_Cluster(double , std::vector <double>);
    bool FileExist (const std::string& name); // check if a file exist


    
//======= since 2023
    // Initialize Variables default values
    void InitializeVariables();
    // updates variables based on the command line arguments
    void UpdateVariables(const std::vector<std::string>& Arguments);
    void UpdateBoxSize(MESH* pmesh, double edge_distance);

    bool ValidateVariable();
private:
    std::string m_MeshFileName;
    std::string m_TaskName;
    Curvature m_CurvatureCalculations;
    Curvature * m_pCurvatureCalculations;
    bool m_Health;
    double m_BilayerThickness;
    double m_BoxDist;
    MESH          m_Mesh;
    MESH          *m_pMesh;

    
};


#endif
