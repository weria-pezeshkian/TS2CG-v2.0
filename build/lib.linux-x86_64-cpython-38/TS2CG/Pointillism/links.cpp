

#include "links.h"

links::links(int id, vertex *v1, vertex *v2, triangle *t1)
{
    m_IntEnergy = 0;
    m_T1=t1;
    m_V1=v1;
    m_V2=v2;
    m_SimTimeStep=-1;
    m_ID=id;
    m_LinkSide = 0;
    m_mirorflag=false;
    m_Show=true;
    m_EdgeSize = 0;
    m_LinkType = 0;
}
links::links(int id)
{
    m_IntEnergy = 0;
    m_ID=id;
    m_LinkSide = 0;
    m_mirorflag=false;
    m_Show=true;
    m_SimTimeStep=-1;
    m_EdgeSize = 0;
    m_LinkType = 0;
}

links::~links()
{
    
}
void links::UpdateTriangle(triangle *v)
{
    m_T1=v;
}
void links::UpdateV0(vertex *v0)
{
    m_V0=v0;
    
}
void links::UpdateV(vertex *v1,vertex *v2,vertex *v3)
{
    m_V3=v3;
    m_V2=v2;
    m_V1=v1;
}
void links::UpdateSimTimeStep(int v)
{
    m_SimTimeStep=v;
}
void links::UpdateV3(vertex *v3)
{
    m_V3=v3;

}

void links::UpdateMirrorLink(links* v)
{
    m_mirorlink=v;
}
void links::UpdateNeighborLink1(links* v)
{
    m_neighborlink1=v;
}
void links::UpdateNeighborLink2(links* v)
{
    m_neighborlink2=v;
}
void links::UpdateVisualize(bool v)
{
    m_Show=v;
}
void links::UpdateMirrorFlag(bool v)
{
    m_mirorflag=v;
}
void links::UpdateNormal()
{

   if(this->GetMirrorFlag()==true)
   {
       Vec3D v2=(m_mirorlink->GetTriangle())->GetNormalVector();
       Vec3D v1=m_T1->GetNormalVector();
       m_Normal=v1+v2;
       double norm=m_Normal.norm();
       m_Normal=m_Normal*(1.0/norm);
       if(norm==0)
       {
           std::cout<<"error 2022----> one of the normals has zero size; normal link cannot be defined  \n";
           exit(0);
       }
       m_mirorlink->PutNormal(m_Normal);
    }
    else
    {
        // this is an edge link
        std::cout<<" link type "<<m_LinkType<<" \n";
        std::cout<<"error ----> normal vector for edge links has not been defined   \n";
        exit(0);
    }

}
void links::PutNormal(Vec3D n)
{
    m_Normal=n;
    
}
void links::PutShapeOperator(Vec3D Be,double He)
{
    m_Be=Be;
    m_He=He;
    
}
void links::UpdateIntEnergy(double en)
{
    m_IntEnergy = en;
}
void links::UpdateShapeOperator(Vec3D *pBox)
{
    UpdateEdgeVector(pBox);
   if(this->GetMirrorFlag()==true)
   {
       Vec3D Re = m_EdgeVector;
       Re=Re*(1.0/m_EdgeSize);
       Vec3D Be=m_Normal*Re;

    //====== Finding the size of the Be vector to make it normaized; this should not be needed
    // just have it so in case.
       double size=Be.norm();
       if(size!=0)
       {
           size=1.0/size;
       }
       else
       {
           std::cout<<" error 7634---> this should not happen \n";
           exit(0);
       }
       Be=Be*size;
       Vec3D Nf1=(m_mirorlink->GetTriangle())->GetNormalVector();
       Vec3D Nf2=m_T1->GetNormalVector();

       
//=== this is different from the orginal paper; it is faster
//==========
       double sign=Re.dot(Nf1*Nf2,Re);
       double tangle=Re.dot(Nf1,Nf2);
       double He=0;

       if(tangle<1)
       {
           if(sign>0)
           {
               He=-m_EdgeSize*sqrt(0.5*(1.0-tangle));
           }
           else if(sign<0)
           {
            He=m_EdgeSize*sqrt(0.5*(1.0-tangle));  //He=2*cos(m_Dihedral/2.0)*renorm;
           }
           else
           {
            He=0;
           }
       }
       else if(tangle>=1 && tangle<1.01)
       {
           // in case some numerical probelm happens; 1.01 is too large however,
           He=0;
   
       }
       else if(tangle>1.01)
       {
           std::cout<<"error--->: somthing wrong with this link \n";
           exit(0);
       }
	
       m_Be=Be;
       m_He=He;
       m_mirorlink->PutShapeOperator(m_Be,m_He);


  }
  else
  {

  }


}
void links::UpdateEdgeVector(Vec3D *pBox)
{
        double x1=m_V1->GetVXPos();
        double y1=m_V1->GetVYPos();
        double z1=m_V1->GetVZPos();
        double x2=m_V2->GetVXPos();
        double y2=m_V2->GetVYPos();
        double z2=m_V2->GetVZPos();
        
        double dx1=x2-x1;
        if(fabs(dx1)>(*pBox)(0)/2.0)
        {
            if(dx1<0)
            dx1=(*pBox)(0)+dx1;
            else if(dx1>0)
            dx1=dx1-(*pBox)(0);
        }
        double dy1=y2-y1;
        if(fabs(dy1)>(*pBox)(1)/2.0)
        {
            if(dy1<0)
            dy1=(*pBox)(1)+dy1;
            else if(dy1>0)
            dy1=dy1-(*pBox)(1);
        }
        double dz1=z2-z1;
        if(fabs(dz1)>(*pBox)(2)/2.0)
        {
            if(dz1<0)
            dz1=(*pBox)(2)+dz1;
            else if(dz1>0)
            dz1=dz1-(*pBox)(2);
        }
      
        Vec3D Re(dx1,dy1,dz1);
        m_EdgeVector = Re;
        m_EdgeSize=(m_EdgeVector.norm());
    
        if(this->GetMirrorFlag()==true)
        {
            m_mirorlink->PutEdgeVector(m_EdgeVector*(-1),m_EdgeSize);
        }

}
void links::PutEdgeVector(Vec3D v, double l)
{
    m_EdgeSize   = l;
    m_EdgeVector = v;
}
void links::Flip()
{
    
   if(this->GetMirrorFlag()==true)
   {
    
    triangle *T2 = m_mirorlink->GetTriangle();
    vertex  *V4 = m_mirorlink->GetV3();
    vertex  *v1 = m_V1;
    vertex  *v2 = m_V2;
    vertex  *v3 = m_V3;
    vertex  *v4 = V4;
    links *l1=this->GetNeighborLink1();
    links *l2=this->GetNeighborLink2();
    links *l3=m_mirorlink->GetNeighborLink1();
    links *l4=m_mirorlink->GetNeighborLink2();

       m_V1->RemoveFromNeighbourVertex(m_V2);
       m_V2->RemoveFromNeighbourVertex(m_V1);
       V4->AddtoNeighbourVertex(m_V3);
       m_V3->AddtoNeighbourVertex(V4);

       m_V1->RemoveFromLinkList(this);
       m_V2->RemoveFromLinkList(m_mirorlink);
       V4->AddtoLinkList(this);
       m_V3->AddtoLinkList(m_mirorlink);
       m_V1->RemoveFromTraingleList(T2);
       m_V2->RemoveFromTraingleList(m_T1);
       m_V3->AddtoTraingleList(T2);
       V4->AddtoTraingleList(m_T1);

       this->UpdateNeighborLink1(l2);
       this->UpdateNeighborLink2(l3);
       m_mirorlink->UpdateNeighborLink1(l4);
       m_mirorlink->UpdateNeighborLink2(l1);
       l1->UpdateV3(V4);
       l2->UpdateV3(V4);
       l3->UpdateV3(m_V3);
       l4->UpdateV3(m_V3);
       l1->UpdateNeighborLink1(m_mirorlink);
       l1->UpdateNeighborLink2(l4);
       l2->UpdateNeighborLink1(l3);
       l2->UpdateNeighborLink2(this);
    
    l3->UpdateNeighborLink1(this);
    l3->UpdateNeighborLink2(l2);
    
    l4->UpdateNeighborLink1(l1);
    l4->UpdateNeighborLink2(m_mirorlink);
    

    l1->UpdateTriangle(T2);
    l3->UpdateTriangle(m_T1);
    l4->UpdateTriangle(T2);
    l2->UpdateTriangle(m_T1);
 
  
    m_V1=v4;
    m_V2=v3;
    m_V3=v1;
    V4=v2;
    m_mirorlink->UpdateV(m_V2,m_V1,V4);



    
    int id2=T2->GetTriID();
    int id1=m_T1->GetTriID();
    triangle tm1(id1,m_V1,m_V2,m_V3);
    triangle tm2(id2,m_V2,m_V1,V4);
      *m_T1 = tm1;
   
      *(m_mirorlink->GetTriangle()) = tm2;
    
   
   }
    else
    {
        std::cout<<"error---> a link without a mirror, possibly an edge link, is asked to be flipped, such an action is not possible \n";
        exit(0);
    }
    
}
// this function are from old time; should be removed soon
void links::ReadLinkFromFile(std::ifstream *inputfile,std::vector <vertex *> pv, std::vector <links *> pL, std::vector <triangle *> pT)
{
    int id,n,v1,v2,v3,l1,l2;
    std::string word;
    bool readok=true;
    (*inputfile)>>word>>id;
    if(id!=m_ID || word!="Link")
        readok=false;
    
    (*inputfile)>>word>>v1>>v2>>v3;
    if(id!=m_ID || word!="Vertices")
        readok=false;
    m_V1=pv.at(v1);
    m_V2=pv.at(v2);
    m_V3=pv.at(v3);
    
    (*inputfile)>>word>>n;
    if(id!=m_ID || word!="Triangle")
        readok=false;
    m_T1=pT.at(n);

    (*inputfile)>>word>>l1>>l2;;
    if(id!=m_ID || word!="Neighborlinks")
        readok=false;
    
    m_neighborlink1=pL.at(l1);
    m_neighborlink2=pL.at(l2);
    
    (*inputfile)>>word>>n;
    if(id!=m_ID || word!="Visualization")
        readok=false;
    m_Show = n;
    
    (*inputfile)>>word>>n;
    if(id!=m_ID || word!="MirorLink")
        readok=false;
    m_mirorlink = pL.at(n);
    if(m_mirorflag == false)
    {
        m_LinkSide = 1;
        m_mirorlink->UpdateMirrorFlag(true);
    }
    else
       m_LinkSide = 2;
       m_mirorflag = true;
    
    
    
    if(readok!=true)
    {
        std::cout<<"Error: wrong link id in reading the file: perhpas the file is broken \n";
    }
}
void links::WriteLinkToFile(std::ofstream *output)
{

    (*output)<<"Link "<<m_ID<<"\n";
    (*output)<<"Vertices "<<m_V1->GetVID()<<" "<<m_V2->GetVID()<<" "<<m_V3->GetVID()<<"\n";
    (*output)<<"Triangle "<<m_T1->GetTriID()<<"\n";
    (*output)<<"Neighborlinks  "<<m_neighborlink1->GetID()<<" "<<m_neighborlink2->GetID()<<" \n";
    (*output)<<"Visualization "<<m_Show<<"\n";
    (*output)<<"MirorLink "<<m_mirorlink->GetID()<<"\n";

    
    
}




