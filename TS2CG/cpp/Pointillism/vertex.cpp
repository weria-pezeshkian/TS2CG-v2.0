 #if !defined(AFX_vertex_CPP_7F4A21C7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_vertex_CPP_7F4A21C7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include "vertex.h"
#include "links.h"
#include "triangle.h"
#include "Nfunction.h"
vertex::vertex()
{
}
vertex::vertex(int id, double x, double y, double z)
{

m_X = x;
m_Y = y;
m_Z = z;
m_ID = id;
m_kappa=1.0;
    m_Group = 0;
    m_OwnInclusion = false;
    m_kappaG = 0;
    m_GroupName = "system";
    
    
    m_Geodesic_Curvature = 0;
    m_Normal_Curvature = 0;
    m_KGC = 0;
    m_KNC = 0;
    m_VertexType = 0;
    m_VLength = 0;                       // length of the vertex
    m_Lambda = 0;                   // line tension
    
    
}
vertex::vertex(int id)
{

m_ID=id;
m_X=0;
m_Y=0;
m_Z=0;
m_kappa=1.0;
    m_Group = 0;
    m_OwnInclusion = false;
    m_kappaG = 0;
    m_GroupName = "system";
    
    
    m_Geodesic_Curvature = 0;
    m_Normal_Curvature = 0;
    m_KGC = 0;
    m_KNC = 0;
    m_VertexType = 0;
    m_VLength = 0;                       // length of the vertex
    m_Lambda = 0;                   // line tension

}
vertex::~vertex()
{
    
}
void vertex::UpdateGroupName(std::string x)
{
    m_GroupName=x;
}
void vertex::UpdateBox(Vec3D *x)
{
m_pBox=x;
}
void vertex::UpdateVID(int i)
{
    m_ID=i;
}
void vertex::UpdateVXPos(double x)
{
	if(x>=(*m_pBox)(0))
	{
		m_X=x-(*m_pBox)(0);
	}
	else if(x<0)
	{
		m_X=x+(*m_pBox)(0);
	}
	else
	{
		m_X=x;
	}
    
}
void vertex::UpdateIsFullDomain(bool x)
{
    m_IsFullDomain=x;
}
void vertex::UpdateDomainID(int x)
{
    m_DomainID = x;
}
void vertex::UpdateVYPos(double x)
{
	if(x>=(*m_pBox)(1))
	{
		m_Y=x-(*m_pBox)(1);
	}
	else if(x<0)
	{
		m_Y=x+(*m_pBox)(1);
	}
	else
	{
		m_Y=x;
	}
}
void vertex::UpdateVZPos(double x)
{
	if(x>=(*m_pBox)(2))
	{
		m_Z=x-(*m_pBox)(2);
	}
	else if(x<0)
	{
		m_Z=x+(*m_pBox)(2);
	}
	else
	{
		m_Z=x;
	}
}
void vertex::AddtoLinkList(links* z)
{
m_VLinkList.push_back(z);
}
void vertex::UpdateKappa(double z1, double z2)
{
    m_kappa=z1;
    m_kappaG=z2;

}
void vertex::RemoveFromLinkList(links* z)
{

m_VLinkList.erase(std::remove(m_VLinkList.begin(), m_VLinkList.end(), z), m_VLinkList.end());
}
void vertex::AddtoTraingleList(triangle* z)
{
m_VTraingleList.push_back(z);
}
void vertex::RemoveFromTraingleList(triangle* z)
{
m_VTraingleList.erase(std::remove(m_VTraingleList.begin(), m_VTraingleList.end(), z), m_VTraingleList.end());
}
void vertex::AddtoNeighbourVertex(vertex* z)
{
m_VNeighbourVertex.push_back(z);
}
void vertex::RemoveFromNeighbourVertex(vertex* z)
{
m_VNeighbourVertex.erase(std::remove(m_VNeighbourVertex.begin(), m_VNeighbourVertex.end(), z), m_VNeighbourVertex.end());
}
void vertex::UpdateInclusion(inclusion * inc)
{
m_pInclusion=inc;
}
void vertex::UpdateOwnInclusion(bool own)
{
m_OwnInclusion=own;
}
void vertex::UpdateGroup(int z)
{
    m_Group=z;
}
void vertex::UpdateEnergy(double z)
{
m_Energy=z;
}
//////////// Functions related to Curvature energy and normal
void vertex::UpdateCurvature(double x,double y)
{
    m_Curvature.clear();
    m_Curvature.push_back(x);
    m_Curvature.push_back(y);
}
void vertex::UpdateNormal_Area(Vec3D v,double a)
{
    m_Normal=v;
    m_Area=a;
}
void vertex::UpdateL2GTransferMatrix(Tensor2 v)
{
    m_T_Local_2_Global = v;
}
void vertex::UpdateG2LTransferMatrix(Tensor2 v)
{
    m_T_Global_2_Local = v;
}
void vertex::NOPBCUpdatePos(Vec3D X)
{
    m_X= X(0);
    m_Y= X(1);
    m_Z= X(2);
}
#endif



