#include "MESH.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
MESH class for quick access
 */
MESH::MESH()
{


}
MESH::~MESH()
{

}
void MESH::GenerateMesh(MeshBluePrint meshblueprint)
{
    m_Box = meshblueprint.simbox;
    m_pBox = &m_Box;



    // Making vertices
    for (std::vector<Vertex_Map>::iterator it = (meshblueprint.bvertex).begin() ; it != (meshblueprint.bvertex).end(); ++it)
    {
            vertex v(it->id,(it->x+0.00000001*(rand() % 100)),it->y,it->z);
           // vertex v(it->id,it->x,it->y,it->z);
            v.UpdateBox(m_pBox);
            v.UpdateGroup(it->domain);
            v.UpdateDomainID(it->domain);
            v.UpdateKappa(0,0);
            v.m_Lambda = 0;
            v.m_KGC = 0;
            v.m_KNC = 0;
            m_Vertex.push_back(v);
    }
//===== Make exclution [since June, 2023]

       for (std::vector<vertex>::iterator it = m_Vertex.begin() ; it != m_Vertex.end(); ++it)
                   m_pActiveV.push_back(&(*it));

    // Making triangles
    int t=0;
    int temid = 0;
    for (std::vector<Triangle_Map>::iterator it = (meshblueprint.btriangle).begin() ; it != (meshblueprint.btriangle).end(); ++it)
    {
        bool pr=true;

        int vid1 = ((meshblueprint.btriangle).at(t)).v1;
        int vid2 = ((meshblueprint.btriangle).at(t)).v2;
        int vid3 = ((meshblueprint.btriangle).at(t)).v3;


        triangle T(temid,&(m_Vertex.at(it->v1)),&(m_Vertex.at(it->v2)),&(m_Vertex.at(it->v3)));
        m_Triangle.push_back(T);
            temid++;
        t++;
    }
    //make inclusions
    for (std::vector<Inclusion_Map>::iterator it = (meshblueprint.binclusion).begin() ; it != (meshblueprint.binclusion).end(); ++it)
    {
        inclusion Tinc(it->id);
        if(m_Vertex.size()<it->vid+1)
        {
        std::cout<<"----> Error: Inclusion vertex id is out of range "<<std::endl;
            exit(0);
        }
        Tinc.Updatevertex(&(m_Vertex.at(it->vid)));
        Tinc.UpdateInclusionTypeID(it->tid);
        Vec3D D(it->x,it->y,0);
        Tinc.UpdateLocalDirection(D);
        m_Inclusion.push_back(Tinc);
        (m_Vertex.at(it->vid)).UpdateOwnInclusion(true);
    }
    //make exclusion
    for (std::vector<Exclusion_Map>::iterator it = (meshblueprint.bexclusion).begin() ; it != (meshblueprint.bexclusion).end(); ++it)
    {
        exclusion Texc(it->id);
        if(m_Vertex.size()<it->vid+1)
        {
        std::cout<<"----> error: exclusion vertex id is out of range "<<std::endl;
            exit(0);
        }
        Texc.Updatevertex(&(m_Vertex.at(it->vid)));
        Texc.UpdateRadius(it->R);
        m_Exclusion.push_back(Texc);
    }
    for (std::vector<exclusion>::iterator it = m_Exclusion.begin() ; it != m_Exclusion.end(); ++it)
        m_pExclusion.push_back(&(*it));

    for (std::vector<inclusion>::iterator it = m_Inclusion.begin() ; it != m_Inclusion.end(); ++it)
        m_pInclusion.push_back(&(*it));

    t=0;
    for (std::vector<triangle>::iterator it = m_Triangle.begin() ; it != m_Triangle.end(); ++it)
        m_pActiveT.push_back(&(*it));

// end Make exclution


    int li=-1;

    for (std::vector<triangle*>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it)
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
        m_Links.push_back(l1);
        m_Links.push_back(l2);
        m_Links.push_back(l3);

    }
    li=-1;
    for (std::vector<triangle*>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it)
    {
        li++;
        int id1=li;
        li++;
        int id2=li;
        li++;
        int id3=li;
        links * l1=&(m_Links.at(id1));
        links * l2=&(m_Links.at(id2));
        links * l3=&(m_Links.at(id3));
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
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        bool foundM=false;
        if((it)->GetMirrorFlag()==true)
        {
            m_pMHL.push_back(it->GetMirrorLink());
            m_pHL.push_back(&(*it));
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
            m_pEdgeL.push_back(&(*it));
            it->m_LinkType = 1;
        }

    }
    int edgelink=0;
    for (std::vector<links*>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it)
        edgelink++;

    if(edgelink!=0)
    {
        std::cout<<"----> Note: "<<edgelink<<" links at the edge \n";
        std::cout<<"----> Note (Warning): the system is not closed! \n";

    }
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {

        m_pActiveL.push_back(&(*it));
    }
    //==== Getting the edge vertex from link
    for (std::vector<links*>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it)
    {
        m_pEdgeV.push_back((*it)->GetV1());
        ((*it)->GetV1())->m_pEdgeLink = *it;
        ((*it)->GetV2())->m_pPrecedingEdgeLink = *it;
        ((*it)->GetV2())->AddtoNeighbourVertex((*it)->GetV1());
        ((*it)->GetV1())->m_VertexType = 1;
    }
    for (std::vector<vertex*>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it)
    {
        bool isedge = false;
        for (std::vector<vertex*>::iterator it2 = m_pEdgeV.begin() ; it2 != m_pEdgeV.end(); ++it2)
        {
            if((*it2)->GetVID()==(*it)->GetVID())
                isedge = true;
        }
        if(isedge==false)
        m_pSurfV.push_back(*it);
    }
    for (std::vector<inclusion*>::iterator it = m_pInclusion.begin() ; it != m_pInclusion.end(); ++it)
        ((*it)->Getvertex())->UpdateInclusion((*it));

    // =======
    // ==== info of the mesh




    std::cout<<"---> active vertex "<<m_pActiveV.size()<<" surf vertex "<<m_pSurfV.size()<<"  edge vertex "<<m_pEdgeV.size()<<"";
    std::cout<<" inclusions "<<m_pInclusion.size()<<"  ";
    std::cout<<" total links "<<m_Links.size()<<"  ";
    std::cout<<" triangles "<<m_pActiveT.size()<<"  \n";





    //WritevtuFiles VTU(pState);
    //std::string file="ini_Mesh.vtu";
    //VTU.Writevtu(m_pAllV,m_pAllT,m_pAllLinks,file);

//==== checking the mesh; this can go to another function. for now we keep it here.
    int no_repeated_link = 0;
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        for (std::vector<links>::iterator it1 = it+1 ; it1 != m_Links.end(); ++it1)
        {
            if(it->GetV1()==it1->GetV1() && it->GetV2()==it1->GetV2())
            {
                no_repeated_link++;

            }

        }
    }
    if(no_repeated_link!=0)
    {
        std::cout<<" error---> approximatly  "<<no_repeated_link/3<<" triangles was found to be inconsisent in their orientation \n";
        exit(0);
    }
}
//===========================================================
// Note, the converted blue print will not have the exclusions
//
MeshBluePrint MESH::Convert_Mesh_2_BluePrint(MESH *mesh)
{
    MeshBluePrint BluePrint;
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    Vec3D simbox;

    // vertex member of the blue print
    std::vector<vertex*> pV = mesh->m_pActiveV;
    for (std::vector<vertex *>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        Vertex_Map tvm;
        tvm.x = (*it)->GetVXPos();
        tvm.y = (*it)->GetVYPos();
        tvm.z = (*it)->GetVZPos();
        tvm.id = (*it)->GetVID();
        tvm.domain = (*it)->GetGroup();
        bvertex.push_back(tvm);
    }
    // triangle map member of the blue print
    std::vector<triangle*> pT = mesh->m_pActiveT;
    for (std::vector<triangle *>::iterator it = pT.begin() ; it != pT.end(); ++it)
    {
        Triangle_Map ttm;
        ttm.v1 = ((*it)->GetV1())->GetVID();
        ttm.v2 = ((*it)->GetV2())->GetVID();
        ttm.v3 = ((*it)->GetV3())->GetVID();
        ttm.id = (*it)->GetTriID();
        btriangle.push_back(ttm);
    }

    // inclusion map member of the blue print
    std::vector<inclusion*> pInc = mesh->m_pInclusion;
    for (std::vector<inclusion *>::iterator it = pInc.begin() ; it != pInc.end(); ++it)
    {
        Inclusion_Map tim;
        tim.x = ((*it)->GetLDirection())(0);
        tim.y = ((*it)->GetLDirection())(1);
        tim.vid = ((*it)->Getvertex())->GetVID();
        tim.tid = ((*it)->GetInclusionTypeID());
        tim.id = ((*it)->GetID());
        binclusion.push_back(tim);

    }
    // Add other map into the mesh map
    BluePrint.bvertex = bvertex;
    BluePrint.btriangle = btriangle;
    BluePrint.binclusion = binclusion;
    BluePrint.simbox = *(mesh->m_pBox);

    return BluePrint;
}

