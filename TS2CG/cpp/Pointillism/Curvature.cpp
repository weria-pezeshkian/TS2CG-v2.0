


#include "Curvature.h"
#include "Tensor2.h"
Curvature::Curvature()
{

}
Curvature::~Curvature()
{

}
void Curvature::SurfVertexCurvature(vertex * pvertex)
{
    m_pVertex=pvertex;
    std::vector<triangle *> Ntr=m_pVertex->GetVTraingleList();
    double Area=0.0;
    Vec3D Normal;
    for (std::vector<triangle *>::iterator it = Ntr.begin() ; it != Ntr.end(); ++it)
    {
        Vec3D v=(*it)->GetNormalVector();
        double a=(*it)->GetArea();
        v=v*a;
        Normal=Normal+v;
        Area+=a;
    }
    Area=Area/3.0;
    if(Area==0)
    {
        std::cout<<Ntr.size()<<"\n";
        std::string sms=" error----> vertex has a zero area \n";
        std::cout<<sms<<"\n";
        exit(0);
    }
    else if(Area<0)
    {
        std::string sms=" error----> vertex has a negetive area \n";
        std::cout<<sms<<"\n";
        exit(0);
    }
    double no=Normal.norm();
    no=1.0/no;
    Normal=Normal*no;
    m_pVertex->UpdateNormal_Area(Normal,Area);
        ///=======
    //=== Shape Operator
    //========
    Tensor2  SV;
    Tensor2 IT('I');

    Tensor2 P=IT-IT.makeTen(Normal);

    Tensor2 Pt=P.Transpose(P);
    std::vector<links *> NLinks=m_pVertex->GetVLinkList();

    // what if the edge vertex has one triangle?
    for (std::vector<links *>::iterator it = NLinks.begin() ; it != NLinks.end(); ++it)
    {
       if((*it)->GetMirrorFlag()==true)
       {
           Vec3D ve=(*it)->GetNormal();
           double we=ve.dot(Normal,ve);
           Vec3D Be=(*it)->GetBe();
           double he=(*it)->GetHe();
           Vec3D Se = P*Be;

           // ff should be 1 but just for sake of numerical errors
           double ff=Se.norm();
           if(ff==0)
           {
               std::cout<<"-----> Error: projection is zero error"<<"\n";
               exit(0);
           }
           else
           {
               Se=Se*(1.0/ff);
           }

           Tensor2 Q=P.makeTen(Se);

           SV=SV+(Q)*(we*he);
       }//if((*it)->GetMirrorFlag()==true)
    }


    /*
   // method in the paper
    Tensor2 SV2;
    for (std::vector<links *>::iterator it = NLinks.begin() ; it != NLinks.end(); ++it)
    {
       if((*it)->GetMirrorFlag()==true)
       {
           Vec3D ve=(*it)->GetNormal();
           double we=ve.dot(Normal,ve);
           Vec3D Be=(*it)->GetBe();
           double he=(*it)->GetHe();
           Tensor2 Se;
           Se=(Se.makeTen(Be))*he;
           SV2=SV2+(Pt*(Se*P))*we;
       }//if((*it)->GetMirrorFlag()==true)
    }
    */

    // end paper


    ///=============
    //==== Find Curvature and local frame
    //=============

    Tensor2 Hous = Householder(Normal);
    Tensor2 LSV;// Local SV
    LSV=(Hous.Transpose(Hous))*(SV*Hous);    // LSV is the curvature matrix in the local frame, the it is a 2*2 minor matrix since all the 3 component are zero.

    double b=LSV(0,0)+LSV(1,1);
    double c=LSV(0,0)*LSV(1,1)-LSV(1,0)*LSV(0,1);


    double delta=b*b-4*c;
    double c1,c2;
    if(delta>0.0)
    {
    delta=sqrt(delta);
    c1=b+delta;
    c1=0.5*c1;      // c1 always will be larger then c2
    c2=b-delta;
    c2=0.5*c2;


    }
    else if (fabs(delta)<0.0001)
    {

    c1=0.5*b;
    c2=c1;

    }
    else
    {
        c1=0;
        c2=0;
        std::cout<<"WARNING: faild to find curvature on vertex "<<pvertex->GetVID()<<"  because delta is "<<delta<<"  c1 and c2 are set to 100 \n";
        std::cout<<" if you face this too much, you should stop the job and .... \n";
    }

    // if(true) in general we do not need this, only if we have directional inclsuions
    {
    Tensor2 EigenvMat('O');

    double p=LSV(0,0);
    double q=LSV(0,1);


    double size=sqrt(q*q+(c1-p)*(c1-p));                   // The Eigenvectors can be calculated using this equation LSV*R=c1*R


        if(size == 0.0){
            size = 1;
            q = 1;
        }

    if(size == 0.0){
        size = 1;
        q = 1;
    }

    EigenvMat(0,0)=q/size;                                  // only one of them needs to be calculated, one is normal vector and the other is perpendicular to first one
    EigenvMat(1,0)=(c1-p)/size;
    EigenvMat(0,1)=-EigenvMat(1,0);
    EigenvMat(1,1)=EigenvMat(0,0);
    EigenvMat(2,2)=1;
        // std::cout<<(Hous*Normal)(0)<<"   "<<(Hous*Normal)(1)<<"   "<<(Hous*Normal)(2)<<"   \n";

 // Tensor2 TransferMatGL=EigenvMat*Hous;   /// This matrix transfers vectors from Global coordinate to local coordinate//
 // Tensor2 TransferMatLG = TransferMatGL.Transpose(TransferMatGL);//// This matrix transfers vectors from local coordinate to global coordinate


        ///  this is correct, We can check by applying transpose(E)*t1 = (1,0,0)

     Tensor2 TransferMatLG=Hous*EigenvMat;   /// This matrix transfers vectors from local coordinate to global coordinate

       /* std::cout<<"===============\n";
        std::cout<<TransferMatLG(0,0)<<" "<<TransferMatLG(0,1)<<" "<<TransferMatLG(0,2)<<" \n";
        std::cout<<TransferMatLG(1,0)<<" "<<TransferMatLG(1,1)<<" "<<TransferMatLG(1,2)<<" \n";
        std::cout<<TransferMatLG(2,0)<<" "<<TransferMatLG(2,1)<<" "<<TransferMatLG(2,2)<<" \n";
        */

     Tensor2 TransferMatGL=TransferMatLG.Transpose(TransferMatLG);   /// This matrix transfers vectors from Global coordinate to local coordinate

        m_pVertex->UpdateL2GTransferMatrix(TransferMatLG);
        m_pVertex->UpdateG2LTransferMatrix(TransferMatGL);
    }
    c1=c1/Area;
    c2=c2/Area;
    m_pVertex->UpdateCurvature(c1,c2);

}
void Curvature::EdgeVertexCurvature(vertex * pvertex)
{
    // first we obtain the vertex area and normal. Area is not important here as the vertex is an edge vertex
    m_pVertex=pvertex;
    std::vector<triangle *> Ntr=m_pVertex->GetVTraingleList();

    double Area=0.0;
    Vec3D Normal;
    for (std::vector<triangle *>::iterator it = Ntr.begin() ; it != Ntr.end(); ++it)
    {
        double a=(*it)->GetArea();
        Vec3D v=(*it)->GetNormalVector();
        v=v*a;
        Normal=Normal+v;
        Area+=a;
    }
    Area=Area/3.0;
    // just in case; this should never happen unless the inputs are wrong
    if(Area<=0)
    {
        std::cout<<Ntr.size()<<"\n";
        std::string sms=" error----> bad area for vertex \n";
        std::cout<<sms<<"\n";
        exit(0);
    }
    double no=Normal.norm();
    no=1.0/no;
    Normal=Normal*no;
    m_pVertex->UpdateNormal_Area(Normal,Area);



    // the shape of the system                          //         v
                                                        //     l1 / \ l2

    links* link1 = m_pVertex->m_pPrecedingEdgeLink;
    links* link2 = m_pVertex->m_pEdgeLink;

    double l1 = link1->m_EdgeSize;
    double l2 = link2->m_EdgeSize;
    double vlenght = 0.5*(l1+l2);
    m_pVertex->m_VLength = vlenght;
    Vec3D L1=link1->m_EdgeVector;
    Vec3D L2=link2->m_EdgeVector;
    L1 = L1*(1/l1);
    L2 = L2*(1/l2);
    Vec3D Norm = (L1-L2)*(1.0/vlenght);    // dT/ds = -Norm ; the size is curvature
    Vec3D Tv = Norm*Normal;  // T at each vertex
    Tv = Tv*(1/Tv.norm());
    Vec3D P = Normal*Tv;
    double cn = Norm.dot(Norm,Normal);
    double cg = Norm.dot(Norm,P);
    m_pVertex->m_Geodesic_Curvature = cg;
    m_pVertex->m_Normal_Curvature =cn;
    Tensor2 TransferMatGL(Tv,P,Normal);     // P1,P2,N is for other surfaces
    Tensor2 TransferMatLG=TransferMatGL.Transpose(TransferMatGL);
    m_pVertex->UpdateL2GTransferMatrix(TransferMatLG);
    m_pVertex->UpdateG2LTransferMatrix(TransferMatGL);




#if TEST_MODE == Enabled

       Vec3D p1(1,0,0);
       Vec3D p2(0,1,0);
        Vec3D p3(0,0,1);

       p1=TransferMatLG*p1;
       p2=TransferMatLG*p2;
        p3=TransferMatLG*p3;
    std::cout<<"=========== \n ";
    std::cout<<p1(0)<<"  "<<p1(1)<<"  "<<p1(2)<<" \n ";
    std::cout<<Tv(0)<<"  "<<Tv(1)<<"  "<<Tv(2)<<" \n ";
    std::cout<<p2(0)<<"  "<<p2(1)<<"  "<<p2(2)<<" \n ";
    std::cout<<P(0)<<"  "<<P(1)<<"  "<<P(2)<<" \n ";
    std::cout<<p3(0)<<"  "<<p3(1)<<"  "<<p3(2)<<" \n ";
    std::cout<<Normal(0)<<"  "<<Normal(1)<<"  "<<Normal(2)<<" \n ";


#endif




}
Tensor2 Curvature::Householder(Vec3D N)
{

    Tensor2 Hous;
    Vec3D Zk;
    Zk(2)=1.0;
    Zk=Zk+N;
    Zk=Zk*(1.0/Zk.norm());

    Tensor2 I('I');
    Tensor2 W=Hous.makeTen(Zk);
    Hous=(I-W*2)*(-1);

    return Hous;




    /*
    Zk(0) = 0;Zk(1) = 0;
    Zk(2)=1.0;
    double SignT=1;
    if((1+N(2))>(1-N(2)))   // 1+/-Normal(2) is always larger then 1
    {
     Zk=Zk+N;
    SignT=-1;
     }
    else if((1+N(2))<=(1-N(2)))
     Zk=Zk-N;
    Zk=Zk*(1.0/Zk.norm());
    W=Hous.makeTen(Zk);
    Hous=(I-W*2)*(SignT);


    std::cout<<"====== old ===========\n";
    std::cout<<Hous(0,0)<<"   "<<Hous(0,1)<<"   "<<Hous(0,2)<<"\n";
    std::cout<<Hous(1,0)<<"   "<<Hous(1,1)<<"   "<<Hous(1,2)<<"\n";
    std::cout<<Hous(2,0)<<"   "<<Hous(2,1)<<"   "<<Hous(2,2)<<"\n";*/
   // double SignT=1;

}
/// normal vector update
Vec3D Curvature::Calculate_Vertex_Normal(vertex *pvertex)
{
    // first we obtain the vertex area and normal.
    std::vector<triangle *> Ntr=pvertex->GetVTraingleList();

    double Area=0.0;
    Vec3D Normal;
    for (std::vector<triangle *>::iterator it = Ntr.begin() ; it != Ntr.end(); ++it)
    {
        double a=(*it)->GetArea();
        Vec3D v=(*it)->GetNormalVector();
        v=v*a;
        Normal=Normal+v;
        Area+=a;
    }
    Area=Area/3.0;
    // just in case; this should never happen unless the inputs are wrong
    if(Area==0)
    {
        std::cout<<Ntr.size()<<"\n";
        std::string sms=" error----> vertex has a zero area \n";
        std::cout<<sms<<"\n";
        exit(0);
    }
    else if(Area<0)
    {
        std::string sms=" error----> vertex has a negetive area \n";
        std::cout<<sms<<"\n";
        exit(0);
    }
    double no=Normal.norm();
    no=1.0/no;
    Normal=Normal*no;
    pvertex->UpdateNormal_Area(Normal,Area);


    return Normal;
}






