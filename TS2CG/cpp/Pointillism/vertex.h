#if !defined(AFX_vertex_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_vertex_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_
#include "SimDef.h"
#include "Vec3D.h"
#include "Tensor2.h"
#include "inclusion.h"
/*******************
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Vertex object: can be defined by and integer and three double numbers for the positions of the vertices
 Also there are a few lists pointing to the nighbouring links and trinagles and also the inclusion if the vertex own one.
 *******************/
class links;
class triangle;
class vertex
{
public:
	vertex(int id, double x, double y, double z);
	vertex(int id);
    vertex();

	 ~vertex();

// A set of function to get vertex variables (since all of them are private)
	    inline int GetVID()                                 {return m_ID;}
        inline double GetVXPos()                            {return m_X;}
        inline double GetVYPos()                            {return m_Y;}
        inline double GetVZPos()                            {return m_Z;}
        inline double GetKappa()                            {return m_kappa;}
        inline double GetKappaG()                           {return m_kappaG;}
        inline double GetArea()                             {return m_Area;}
        inline Tensor2  GetL2GTransferMatrix()              {return m_T_Local_2_Global;}
        inline Tensor2  GetG2LTransferMatrix()              {return m_T_Global_2_Local;}
        inline Vec3D GetNormalVector()                      {return m_Normal;}
        inline std::vector <double> GetCurvature()          {return m_Curvature;}// surface curvature
        inline double GetEnergy()                           {return m_Energy;}
        inline std::vector <links *> GetVLinkList()             {return m_VLinkList;}
        inline std::vector <triangle *> GetVTraingleList()         {return m_VTraingleList;}
        inline std::vector <vertex *> GetVNeighbourVertex()     {return m_VNeighbourVertex;}
        inline inclusion* GetInclusion()                    {return m_pInclusion;}
        inline bool VertexOwnInclusion()                    {return m_OwnInclusion;}
        inline Vec3D *GetBox()                              {return m_pBox;}
        inline int GetGroup()                               {return m_Group;}
        inline std::string GetGroupName()                               {return m_GroupName;}
    inline int GetDomainID()                               {return m_DomainID;}
    inline bool  GetIsFullDomain()                               {return m_IsFullDomain;}
    
public:
  // A set of functions to update vertex variables
  void UpdateVXPos(double x);    // a function for update the x position of a vertex
  void UpdateVYPos(double y);   // a function for update the y position of a vertex
  void UpdateVZPos(double z);   // a function for updates the z position of a vertex
  void NOPBCUpdatePos(Vec3D z);
  void UpdateGroupName(std::string z); // A vertex can only have one group name, is different from group id
  void UpdateKappa(double z1,double z2);
  void UpdateBox(Vec3D *z);
  void UpdateCurvature(double,double); // vis
  void UpdateEnergy(double); 
  void UpdateNormal_Area(Vec3D,double); // vis
  void UpdateL2GTransferMatrix(Tensor2 v);
  void UpdateG2LTransferMatrix(Tensor2 v);
  void UpdateInclusion(inclusion* );
  void UpdateOwnInclusion(bool );
  void UpdateVID(int); // this should never be called, only for temporarily vis functions
  void AddtoLinkList(links* z);
  void AddtoTraingleList(triangle * z);
  void AddtoNeighbourVertex(vertex* z);
  void RemoveFromLinkList(links* z);
  void RemoveFromTraingleList(triangle * z);
  void RemoveFromNeighbourVertex(vertex* z);
  void UpdateGroup(int z);
    void UpdateIsFullDomain(bool z);
    void UpdateDomainID(int );                    // update domain ID, only in version 1.1 and above

private:

    int m_ID;         // ID of the vertex, a unique number
    double m_X;       // X coordinate
    double m_Y;       // Y coordinate
    double m_Z;       // Z coordinate
    double m_kappa;   // bending rigidity
    double m_kappaG;  // gaussian  rigidity
    std::vector <triangle *> m_VTraingleList;  // A list holding all the nighbouring triangles
    std::vector <links *> m_VLinkList;      // A list holding all the nighbouring edges (linkss)
    std::vector <vertex *> m_VNeighbourVertex; // A list holding all the nighbouring vertexes
    inclusion *m_pInclusion;                    // pointer to an inclusion that the vertex hold (could be empty)
    bool m_OwnInclusion;                        // to check if the vertex own any inclusion
    double m_Area;                              // area of the vertex
    int m_Group;            // Id of a group that the vertex belong too
    private:
    Vec3D m_Normal;
    Vec3D *m_pBox;
    std::vector<double> m_Curvature;
    double m_Energy;
    Tensor2  m_T_Local_2_Global;         //  Local to global transformation matrix
    Tensor2  m_T_Global_2_Local;        //  global to local transformation matrix
    std::string m_GroupName;
    
    int m_DomainID;
    bool m_IsFullDomain;
    
    
public:
    // lets have them public for now
    //================================================
    // development of Aug 2023; for membranes with hole
    //a*kg^2+b*k^2    note: any spontaneous one will come from inclusions
    
    // we might need to define n and t and b as well, for now lets ignore it
    //================================================
    double m_Geodesic_Curvature;          // Edge Vertex Curvature
    double m_Normal_Curvature;          // Edge Vertex Curvature
    double m_KGC;                          // model parameters for Geodesic Curvature
    double m_KNC;                           // model parameters for normal Curvature
    double m_VLength;                       // length of the vertex
    double m_Lambda;                   // line tension
    int m_VertexType;                   // 0 surface vertex; 1 edge vertex;
    links * m_pEdgeLink;
    links * m_pPrecedingEdgeLink;// preceding link at the edge
    
    
    
};


#endif
