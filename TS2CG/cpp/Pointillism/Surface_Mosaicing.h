#if !defined(AFX_Surface_Mosaicing_H_7B4B21B8_D13C_9321_BF23_124095086234__INCLUDED_)
#define AFX_Surface_Mosaicing_H_7B4B21B8_D13C_9321_BF23_124095086234__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "links.h"
#include "triangle.h"
#include "Vec3D.h"
#include "inclusion.h"
#include "MESH.h"

class Surface_Mosaicing
{
public:
    Surface_Mosaicing(std::string altype, bool smooth);
    Surface_Mosaicing();
	 ~Surface_Mosaicing();




public:

//    inline std::vector <links *> GetLinks()             {return m_pFL;}


    
public:
    void PerformMosaicing(MESH * pMesh);
    void RoughnessOfALink(links *l, double *linklength, double *midpointdistance);
private:

Vec3D m_Box;
Vec3D *m_pBox;
    std::string m_AlgorithmType;
    
    
private:
    void UpdateGeometry(MESH *pmesh);
    void MosaicOneRound(MESH * pMesh);
    void BestEstimateOfMidPointPossition(links *, double *x, double *y,double *z);
    Tensor2 NormalCoord(Vec3D N);

private:
    std::vector<inclusion* > m_Inc;

    std::vector<vertex* > m_pFV;
    std::vector<triangle* > m_pFT;
    std::vector<links* >  m_pFL;
    bool m_smooth;
public:
    MESH *m_pMesh;
    MESH m_Mesh;

    
    // since 2023
private:
    void  GenerateMidVForAllLinks(std::vector<links *> vlink);

//---
    bool m_Mash_IS_Smooth;
    
};


#endif
