


#include "Surface_Mosaicing.h"
#include "Nfunction.h"
#include "VMDOutput.h"
#include "WriteFiles.h"
#include "Curvature.h"

Surface_Mosaicing::Surface_Mosaicing(std::string altype, bool smooth)
{
        m_AlgorithmType = altype;
        m_smooth = smooth;

}Surface_Mosaicing::Surface_Mosaicing()
{
    m_AlgorithmType = "Type1";
    m_smooth = false;
}
void Surface_Mosaicing::PerformMosaicing(MESH * pMesh)
{
    m_pBox = pMesh->m_pBox;
    m_Mesh.m_Box = *m_pBox;
    m_Mesh.m_pBox = &(m_Mesh.m_Box);
    MosaicOneRound(pMesh);
    UpdateGeometry(m_pMesh);
}
Surface_Mosaicing::~Surface_Mosaicing()
{
    
}
void  Surface_Mosaicing::GenerateMidVForAllLinks(std::vector<links *> vlink)
{
    //=========================================
        //===== for each link, generate a mid vertex
    //=========================================
        int inisize = (m_Mesh.m_Vertex).size();
        int id=inisize;
        for (std::vector<links *>::iterator it = vlink.begin() ; it != vlink.end(); ++it)
        {
                double x,y,z;
                BestEstimateOfMidPointPossition((*it), &x, &y,&z);
                if(isnan(x))
                {
                    std::cout<<"error---> estimate of the mid point is bad "<<x<<"  "<<y<<"  "<<z<<"\n";
                    exit(0);
                }
                vertex v(id,x,y,z);
                v.UpdateBox(m_pBox);
            //==== for version 1.1 and above
                int dom1 = ((*it)->GetV1())->GetDomainID();
                int dom2 = ((*it)->GetV2())->GetDomainID();
                if((*it)->GetV1()->m_VertexType == 1 && (*it)->GetV2()->m_VertexType == 1) {
                    v.m_VertexType = 1;
                }
                else{
                    v.m_VertexType = 0;
                }
                double domain = 0;
                if (dom1==dom2)
                {
                    domain = dom1;
                }
                else
                {
                    v.UpdateIsFullDomain(false);
                    bool dtype1 =  ((*it)->GetV1())->GetIsFullDomain();
                    bool dtype2 =  ((*it)->GetV2())->GetIsFullDomain();
                    if(dtype2==false)
                        domain = dom1;
                    else if(dtype1==false)
                        domain = dom2;
                    else
                        domain = dom1;
                }
                v.UpdateDomainID(domain);
            (m_Mesh.m_Vertex).push_back(v);
            id++;
        }
    // note, we should first find all the link mid point and then copy them
    id=inisize;
    for (std::vector<links *>::iterator it = vlink.begin() ; it != vlink.end(); ++it)
    {

        (*it)->UpdateV0(&((m_Mesh.m_Vertex)[id]));
        if((*it)->GetMirrorFlag()==true)
        {
        ((*it)->GetMirrorLink())->UpdateV0(&((m_Mesh.m_Vertex)[id]));
        }
        id++;
    }
}
void Surface_Mosaicing::MosaicOneRound(MESH * pMesh)
{    ///
//---> First we copy the old vertices into the new vertices only the position, box and incs
    m_Mesh.m_Inclusion = pMesh->m_Inclusion;
    m_Mesh.m_Exclusion = pMesh->m_Exclusion;

    for (std::vector<vertex *>::iterator it = (pMesh->m_pActiveV).begin() ; it != (pMesh->m_pActiveV).end(); ++it)
    {
        double x=(*it)->GetVXPos();
        double y=(*it)->GetVYPos();
        double z=(*it)->GetVZPos();
        int id=(*it)->GetVID();
        vertex v(id,x,y,z);
        if((*it)->VertexOwnInclusion()==true)
        {
            v.UpdateOwnInclusion(true);
            int incid = ((*it)->GetInclusion())->GetID();
            v.UpdateInclusion(&((m_Mesh.m_Inclusion)[incid]));
        }
        v.UpdateDomainID((*it)->GetDomainID());
        v.UpdateBox(m_pBox);
        (m_Mesh.m_Vertex).push_back(v);
    }

//------> finding the mid point
    std::vector<links *> vlink = pMesh->m_pHL;
    std::vector<links *> elink = pMesh->m_pEdgeL;
    vlink.insert(vlink.end(), elink.begin(), elink.end());
    GenerateMidVForAllLinks(vlink);
//-------> Now we have all the vertices
//-------> we can give the inclusions a vertex
    for (std::vector<inclusion>::iterator it = (m_Mesh.m_Inclusion).begin() ; it != (m_Mesh.m_Inclusion).end(); ++it)
    {
        int vid = (it->Getvertex())->GetVID();
        if(vid>=(m_Mesh.m_Vertex).size())
            std::cout<<"error 3234---> should not happen \n";
        it->Updatevertex(&((m_Mesh.m_Vertex)[vid]));
        (m_Mesh.m_pInclusion).push_back(&(*it));

    }
    for (std::vector<exclusion>::iterator it = (m_Mesh.m_Exclusion).begin() ; it != (m_Mesh.m_Exclusion).end(); ++it)
    {
        int vid = (it->Getvertex())->GetVID();
        if(vid>=(m_Mesh.m_Vertex).size())
            std::cout<<"error 3234---> should not happen \n";
        it->Updatevertex(&((m_Mesh.m_Vertex)[vid]));
        (m_Mesh.m_pExclusion).push_back(&(*it));

    }

//----> generate new triangles
    
        int tid=0;
        for (std::vector<links *>::iterator it = (pMesh->m_pActiveL).begin() ; it != (pMesh->m_pActiveL).end(); ++it)
        {
            triangle *t1=(*it)->GetTriangle();
            if(t1->GetGotMashed()==false)
            {
                links * l1 = (*it)->GetNeighborLink1();
                links * l2 = (*it)->GetNeighborLink2();
                vertex *V1=&((m_Mesh.m_Vertex)[(((*it)->GetV1())->GetVID())]);
                vertex *V2=&((m_Mesh.m_Vertex)[(((*it)->GetV2())->GetVID())]);
                vertex *V3=&((m_Mesh.m_Vertex)[(((*it)->GetV3())->GetVID())]);
                vertex *VM0=(*it)->GetV0();
                vertex *VM1=l1->GetV0();
                vertex *VM2=l2->GetV0();
                triangle T1(tid,V1,VM0,VM2);
                tid++;
                triangle T2(tid,VM0,V2,VM1);
                tid++;
                triangle T3(tid,VM0,VM1,VM2);
                tid++;
                triangle T4(tid,VM2,VM1,V3);
                tid++;
                (m_Mesh.m_Triangle).push_back(T1);
                (m_Mesh.m_Triangle).push_back(T2);
                (m_Mesh.m_Triangle).push_back(T3);
                (m_Mesh.m_Triangle).push_back(T4);
                t1->UpdateGotMashed(true);
            }
        }
    //===
    for (std::vector<vertex >::iterator it = (m_Mesh.m_Vertex).begin() ; it != (m_Mesh.m_Vertex).end(); ++it)
    {
        (m_Mesh.m_pActiveV).push_back(&(*it));
    }
    for (std::vector<triangle >::iterator it = (m_Mesh.m_Triangle).begin() ; it != (m_Mesh.m_Triangle).end(); ++it)
    {
        (m_Mesh.m_pActiveT).push_back(&(*it));
    }
    
    
    int li=-1;
    
    for (std::vector<triangle*>::iterator it = (m_Mesh.m_pActiveT).begin() ; it != (m_Mesh.m_pActiveT).end(); ++it)
    {
        ((*it)->GetV1())->AddtoTraingleList((*it));
        ((*it)->GetV1())->AddtoNeighbourVertex(((*it)->GetV2()));
        ((*it)->GetV2())->AddtoTraingleList((*it));
        ((*it)->GetV2())->AddtoNeighbourVertex(((*it)->GetV3()));
        ((*it)->GetV3())->AddtoTraingleList((*it));
        ((*it)->GetV3())->AddtoNeighbourVertex(((*it)->GetV1()));
        
        /// create links
        li++;
        int id1=li;
        li++;
        int id2=li;
        li++;
        int id3=li;
        
        links l1(id1,(*it)->GetV1(),(*it)->GetV2(),(*it));
        l1.UpdateV3((*it)->GetV3());
        
        links l2(id2,(*it)->GetV2(),(*it)->GetV3(),(*it));
        l2.UpdateV3((*it)->GetV1());
        
        links l3(id3,(*it)->GetV3(),(*it)->GetV1(),(*it));
        l3.UpdateV3((*it)->GetV2());
        (m_Mesh.m_Links).push_back(l1);
        (m_Mesh.m_Links).push_back(l2);
        (m_Mesh.m_Links).push_back(l3);
        
    }
    li=-1;
    for (std::vector<triangle*>::iterator it = (m_Mesh.m_pActiveT).begin() ; it != (m_Mesh.m_pActiveT).end(); ++it)
    {
        li++;
        int id1=li;
        li++;
        int id2=li;
        li++;
        int id3=li;
        links * l1=&((m_Mesh.m_Links).at(id1));
        links * l2=&((m_Mesh.m_Links).at(id2));
        links * l3=&((m_Mesh.m_Links).at(id3));
        l1->UpdateNeighborLink1(l2);
        l1->UpdateNeighborLink2(l3);
        l2->UpdateNeighborLink1(l3);
        l2->UpdateNeighborLink2(l1);
        l3->UpdateNeighborLink1(l1);
        l3->UpdateNeighborLink2(l2);
        
        
        ((*it)->GetV1())->AddtoLinkList(l1);
        ((*it)->GetV2())->AddtoLinkList(l2);
        ((*it)->GetV3())->AddtoLinkList(l3);
        
    }
    for (std::vector<links>::iterator it = (m_Mesh.m_Links).begin() ; it != (m_Mesh.m_Links).end(); ++it)
    {
        (m_Mesh.m_pActiveL).push_back(&(*it));
    }
    for (std::vector<links>::iterator it = (m_Mesh.m_Links).begin() ; it != (m_Mesh.m_Links).end(); ++it)
    {
        bool foundM=false;
        if((it)->GetMirrorFlag()==true)
        {
            (m_Mesh.m_pMHL).push_back(it->GetMirrorLink());
            (m_Mesh.m_pHL).push_back(&(*it));
            foundM = true;
        }
        else
        {
            vertex *v1=it->GetV1();
            vertex *v2=it->GetV2();
            
            std::vector<links*>  lList = v2->GetVLinkList();
            for (std::vector<links*>::iterator it2 = lList.begin() ; it2 != lList.end(); ++it2)
            {
                if(((*it2)->GetV2())->GetVID()==v1->GetVID())
                {
                    it->UpdateMirrorLink((*it2));
                    (*it2)->UpdateMirrorLink(&(*it));
                    it->UpdateMirrorFlag(true);
                    (*it2)->UpdateMirrorFlag(true);
                    foundM = true;
                    break;
                }
            }
        }
        if(foundM == false)
        {
            (m_Mesh.m_pEdgeL).push_back(&(*it));
            it->m_LinkType = 1;
        }
    }
    //==== Getting the edge vertex from link
    for (std::vector<links*>::iterator it = (m_Mesh.m_pEdgeL).begin() ; it != (m_Mesh.m_pEdgeL).end(); ++it)
    {
        (m_Mesh.m_pEdgeV).push_back((*it)->GetV1());
        ((*it)->GetV1())->m_pEdgeLink = *it;
        ((*it)->GetV2())->m_pPrecedingEdgeLink = *it;
        ((*it)->GetV2())->AddtoNeighbourVertex((*it)->GetV1());
        ((*it)->GetV1())->m_VertexType = 1;
    }
    for (std::vector<vertex*>::iterator it = (m_Mesh.m_pActiveV).begin() ; it != (m_Mesh.m_pActiveV).end(); ++it)
    {
        bool isedge = false;
        for (std::vector<vertex*>::iterator it2 = (m_Mesh.m_pEdgeV).begin() ; it2 != (m_Mesh.m_pEdgeV).end(); ++it2)
        {
            if((*it2)->GetVID()==(*it)->GetVID())
                isedge = true;
        }
        if(isedge==false)
        (m_Mesh.m_pSurfV).push_back(*it);
    }
    m_pMesh = &m_Mesh;
}

void Surface_Mosaicing::BestEstimateOfMidPointPossition(links *l, double *X, double *Y,double *Z)
{
    double x=0;
    double y=0;
    double z=0;
    vertex * pv1=l->GetV1();
    vertex * pv2=l->GetV2();

    Vec3D *pBox=pv1->GetBox();
    double x1=pv1->GetVXPos();
    double y1=pv1->GetVYPos();
    double z1=pv1->GetVZPos();
               
    double x2=pv2->GetVXPos();
    double y2=pv2->GetVYPos();
    double z2=pv2->GetVZPos();

    double xmid=(x1+x2)/2.0;
    double ymid=(y1+y2)/2.0;
    double zmid=(z1+z2)/2.0;


        if(fabs(x1-x2)>(*pBox)(0)/2)
                xmid=xmid+(*pBox)(0)/2;
        if(fabs(y1-y2)>(*pBox)(1)/2)
                ymid=ymid+(*pBox)(1)/2;
        if(fabs(z1-z2)>(*pBox)(2)/2)
            zmid=zmid+(*pBox)(2)/2;

     Vec3D geodesic_dir(x2-x1,y2-y1,z2-z1);
     for (int i=0;i<3;i++)
     {
         if(fabs(geodesic_dir(i))>(*pBox)(i)/2)
         {
             if(geodesic_dir(i)<0)
             geodesic_dir(i)=geodesic_dir(i)+(*pBox)(i);
             else if(geodesic_dir(i)>0)
             geodesic_dir(i)=geodesic_dir(i)-(*pBox)(i);
         }
     }

    double Linklenght= geodesic_dir.norm();
    geodesic_dir = geodesic_dir*(1/geodesic_dir.norm());

     Vec3D Lo_geoV1=(pv1->GetG2LTransferMatrix())*geodesic_dir;
     Lo_geoV1(2)=0;
     Lo_geoV1=Lo_geoV1*(1/(Lo_geoV1.norm()));
     Vec3D Glo_geoV1=(pv1->GetL2GTransferMatrix())*Lo_geoV1;
    
    Vec3D Lo_geoV2=(pv2->GetG2LTransferMatrix())*geodesic_dir;
    Lo_geoV2(2)=0;
    Lo_geoV2=Lo_geoV2*(1/(Lo_geoV2.norm()));
    Vec3D Glo_geoV2=(pv2->GetL2GTransferMatrix())*Lo_geoV2;

    
    Tensor2 Hous = NormalCoord(geodesic_dir);

    Vec3D t_2=Hous*Glo_geoV2;
    Vec3D t_1=Hous*Glo_geoV1;


    t_2=t_2*(1/t_2(2));
    t_1=t_1*(1/t_1(2));

    t_2=t_2*(Linklenght/2.0);
    t_1=t_1*(Linklenght/2.0);
    
    /// test case to see if inclduign curvature make it better
    Vec3D Dr(0,0,0);

    if(m_AlgorithmType == "Type2")
    {


    Vec3D N1=pv1->GetNormalVector();
    Vec3D N2=pv2->GetNormalVector();
    std::vector <double> C1=pv1->GetCurvature();
    std::vector <double> C2=pv2->GetCurvature();

        double Cos1=Lo_geoV1(0);
        double Sin1=Lo_geoV1(1);
        
        double Cos2=Lo_geoV2(0);
        double Sin2=Lo_geoV2(1);

        double Curve1=C1.at(0)*Cos1*Cos1+C1.at(1)*Sin1*Sin1;
        double Curve2=C2.at(0)*Cos2*Cos2+C2.at(1)*Sin2*Sin2;
        
        
        double D2X_1=Curve1*(t_1.dot(t_1,t_1))*(2*N1(2)*t_1(0)/Linklenght-N1(0));
        double D2Y_1=Curve1*(t_1.dot(t_1,t_1))*(2*N1(2)*t_1(1)/Linklenght-N1(1));
        double D2X_2=Curve2*(t_2.dot(t_2,t_2))*(2*N2(2)*t_2(0)/Linklenght-N2(0));
        double D2Y_2=Curve2*(t_2.dot(t_2,t_2))*(2*N2(2)*t_2(1)/Linklenght-N2(1));
    double X_0=(D2X_1+D2X_2+5*(t_1(0)-t_2(0)))/16.0;
    double Y_0=(D2Y_1+D2Y_2+5*(t_1(1)-t_2(1)))/16.0;
  
        
        Dr(0)=X_0;
        Dr(1)=Y_0;
    }
  else if(m_AlgorithmType == "Type1")
  {
 
      Dr(0)=(t_1(0)-t_2(0))/4;
      Dr(1)=(t_1(1)-t_2(1))/4;

  }
  else
    {
        std::cout<<"Error: 12344! \n";
    }
    // For highly rough surfaces
    {
        double drsize=Dr.norm();
       // std::cout<<drsize<<"   "<<Linklenght <<"\n";
     if(m_smooth==true)
     {
         if(drsize>0.5*Linklenght)
         {
             
             Dr = Dr*(0.2*Linklenght/drsize);
         }
     }
    else
    {
        if(drsize>0.5*Linklenght)
        {
            // no much doing
           /* pv1->UpdateOwnInclusion(true);
            pv1->UpdateInclusion(m_Inc.at(0));
            pv2->UpdateOwnInclusion(true);
            pv2->UpdateInclusion(m_Inc.at(0));*/
          
          std::cout<<"warning: the surfaces is very rough, you may use option -smooth \n";
        }
    }
        
        
    }
    Vec3D GDr=(Hous.Transpose(Hous))*Dr;
    
    x=xmid+GDr(0);
    y=ymid+GDr(1);
    z=zmid+GDr(2);
    *X=x;
    *Y=y;
    *Z=z;
    
}
void  Surface_Mosaicing::UpdateGeometry(MESH *pmesh)
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
// This is for minimazation
void Surface_Mosaicing::RoughnessOfALink(links *l, double *linklength, double *midpointdistance)
{
    double x=0;
    double y=0;
    double z=0;
    vertex * pv1=l->GetV1();
    vertex * pv2=l->GetV2();
    Vec3D *pBox=pv1->GetBox();
    double x1=pv1->GetVXPos();
    double y1=pv1->GetVYPos();
    double z1=pv1->GetVZPos();
    
    double x2=pv2->GetVXPos();
    double y2=pv2->GetVYPos();
    double z2=pv2->GetVZPos();
    
    double xmid=(x1+x2)/2.0;
    double ymid=(y1+y2)/2.0;
    double zmid=(z1+z2)/2.0;
    
    
    if(fabs(x1-x2)>(*pBox)(0)/2)
        xmid=xmid+(*pBox)(0)/2;
    if(fabs(y1-y2)>(*pBox)(1)/2)
        ymid=ymid+(*pBox)(1)/2;
    if(fabs(z1-z2)>(*pBox)(2)/2)
        zmid=zmid+(*pBox)(2)/2;
    
    Vec3D geodesic_dir(x2-x1,y2-y1,z2-z1);
    
    for (int i=0;i<3;i++)
    {
        if(fabs(geodesic_dir(i))>(*pBox)(i)/2)
        {
            if(geodesic_dir(i)<0)
                geodesic_dir(i)=geodesic_dir(i)+(*pBox)(i);
            else if(geodesic_dir(i)>0)
                geodesic_dir(i)=geodesic_dir(i)-(*pBox)(i);
        }
    }
    
    double Linklenght= geodesic_dir.norm();
    geodesic_dir = geodesic_dir*(1/geodesic_dir.norm());
    
    Vec3D Lo_geoV1=(pv1->GetG2LTransferMatrix())*geodesic_dir;
    Lo_geoV1(2)=0;
    Lo_geoV1=Lo_geoV1*(1/(Lo_geoV1.norm()));
    Vec3D Glo_geoV1=(pv1->GetL2GTransferMatrix())*Lo_geoV1;
    
    Vec3D Lo_geoV2=(pv2->GetG2LTransferMatrix())*geodesic_dir;
    Lo_geoV2(2)=0;
    Lo_geoV2=Lo_geoV2*(1/(Lo_geoV2.norm()));
    Vec3D Glo_geoV2=(pv2->GetL2GTransferMatrix())*Lo_geoV2;
    
    
    
    
    Tensor2 Hous;
    Vec3D Zk;
    Zk(2)=1.0;
    double SignT=1;
    
    if((1+geodesic_dir(2))>(1-geodesic_dir(2)))
    {
        Zk=Zk+geodesic_dir;
        SignT=-1;
        
    }
    else if((1+geodesic_dir(2))<=(1-geodesic_dir(2)))
    {
        Zk=Zk-geodesic_dir;
    }
    
    Zk=Zk*(1.0/Zk.norm());
    
    
    Tensor2 I('I');
    Tensor2 W=Hous.makeTen(Zk);
    
    Hous=(I-(W*(W.Transpose(W)))*2)*SignT;
    
    
    
    
    Vec3D t_2=Hous*Glo_geoV2;
    Vec3D t_1=Hous*Glo_geoV1;
    
    t_2=t_2*(1/t_2(2));
    t_1=t_1*(1/t_1(2));
    
    
    
    t_2=t_2*(Linklenght/2.0);
    t_1=t_1*(Linklenght/2.0);
    
    /// test case to see if inclduign curvature make it better
    Vec3D Dr(0,0,0);
    
    if(m_AlgorithmType == "Type2")
    {
        
        
        Vec3D N1=pv1->GetNormalVector();
        Vec3D N2=pv2->GetNormalVector();
        std::vector <double> C1=pv1->GetCurvature();
        std::vector <double> C2=pv2->GetCurvature();
        
        double Cos1=Lo_geoV1(0);
        double Sin1=Lo_geoV1(1);
        
        double Cos2=Lo_geoV2(0);
        double Sin2=Lo_geoV2(1);
        
        double Curve1=C1.at(0)*Cos1*Cos1+C1.at(1)*Sin1*Sin1;
        double Curve2=C2.at(0)*Cos2*Cos2+C2.at(1)*Sin2*Sin2;
        
        
        double D2X_1=Curve1*(t_1.dot(t_1,t_1))*(2*N1(2)*t_1(0)/Linklenght-N1(0));
        double D2Y_1=Curve1*(t_1.dot(t_1,t_1))*(2*N1(2)*t_1(1)/Linklenght-N1(1));
        double D2X_2=Curve2*(t_2.dot(t_2,t_2))*(2*N2(2)*t_2(0)/Linklenght-N2(0));
        double D2Y_2=Curve2*(t_2.dot(t_2,t_2))*(2*N2(2)*t_2(1)/Linklenght-N2(1));
        
        
        
        
        
        
        double X_0=(D2X_1+D2X_2+5*(t_1(0)-t_2(0)))/16.0;
        double Y_0=(D2Y_1+D2Y_2+5*(t_1(1)-t_2(1)))/16.0;
        
        
        Dr(0)=X_0;
        Dr(1)=Y_0;
    }
    else if(m_AlgorithmType == "Type1")
    {
        
        Dr(0)=(t_1(0)-t_2(0))/4;
        Dr(1)=(t_1(1)-t_2(1))/4;
        
    }
    else
    {
        std::cout<<"Error: 12344! \n";
    }

        double drsize=Dr.norm();
        *linklength = Linklenght;
        *midpointdistance = drsize;

    
}
Tensor2 Surface_Mosaicing::NormalCoord(Vec3D N)
{

    Tensor2 Hous;
    if(N(2)==-1)
    {
        N(1)=0.00000001;
        N=N*(1/N.norm());
        
    }
        Vec3D Zk;
        Zk(2)=1.0;
        Zk=Zk+N;
        Zk=Zk*(1.0/Zk.norm());
        
        Tensor2 I('I');
        Tensor2 W=Hous.makeTen(Zk);
        Hous=(I-W*2)*(-1);
    return Hous;
    
}

