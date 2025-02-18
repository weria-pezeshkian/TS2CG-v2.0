

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include "CreateMashBluePrint.h"
#include "Nfunction.h"
CreateMashBluePrint::CreateMashBluePrint()
{
}
MeshBluePrint CreateMashBluePrint::MashBluePrintFromInput_Top(std::string inputfilename, std::string topfilename)
{
// initializing the box size
    m_Box(0) = 1;
    m_Box(1) = 1;
    m_Box(2) = 1;

    m_InputFileName = inputfilename;
    m_TopologyFileName = topfilename;

    //==== We read the inclusions type from the input file, we create inclusion types
    ReadTopology(topfilename);

    m_MeshBluePrint.bvertex = m_VertexMap;
    m_MeshBluePrint.btriangle = m_TriangleMap;
    m_MeshBluePrint.binclusion = m_InclusionMap;
    m_MeshBluePrint.simbox = m_Box;
    m_MeshBluePrint.bexclusion = m_ExclusionMap;
    return m_MeshBluePrint;
}
CreateMashBluePrint::~CreateMashBluePrint()
{

}
void CreateMashBluePrint::ReadTopology(std::string file)
{

        std::string ext = file.substr(file.find_last_of(".") + 1);
        if(ext=="tsi")
        {
            Read_TSIFile(file);

        }
        else if(ext=="top")
        {
            Read_Mult_QFile(file);
            GenerateIncFromInputfile();
        }


}
void CreateMashBluePrint::Read_TSIFile(std::string tsifile)
{
    Nfunction f;
    std::ifstream tsi;
    tsi.open(tsifile.c_str());
    std::string str;
    int nver,ntr,ninc,nexc,id;
    while (true)
    {
        tsi>>str;
        if(tsi.eof())
            break;
        if(str=="version")
        {
            getline(tsi,str);
        }
        else if(str=="box")
        {
            getline(tsi,str);
            std::vector<std::string> S = f.split(str);
            if(S.size()<3)
            {
                std::cout<<"---> Error, information of the box is not sufficent in the tsi file \n";
                exit(0);
            }
            else
            {
                m_Box(0) = f.String_to_Double(S.at(0));
                m_Box(1) = f.String_to_Double(S.at(1));
                m_Box(2) = f.String_to_Double(S.at(2));
            }
#if DEBUG_MODE == Enabled
            std::cout<<"----> box was read  "<<std::endl;
            std::cout<<m_Box(0)<<"  "<<m_Box(1)<<"  "<<m_Box(2)<<"  "<<str<<std::endl;
#endif
        }
        else if(str=="vertex")
        {
            tsi>>nver;
            getline(tsi,str);

            for (int i=0;i<nver;i++)
            {
                getline(tsi,str);
                std::vector<std::string> S = f.split(str);

                if(S.size()<4)
                {
                    std::cout<<"error ---> information of the vertex "<<i<<" is not sufficent in the tsi file \n";
                    exit(0);
                }
                else
                {
                    Vertex_Map v;
                    v.id=i;
                    v.include = true;
                    v.x=f.String_to_Double(S[1]);
                    v.y=f.String_to_Double(S[2]);
                    v.z=f.String_to_Double(S[3]);
                    if(S.size()>4)
                        v.domain=f.String_to_Int(S.at(4));
                    else
                        v.domain=0;
                    m_VertexMap.push_back(v);
                }
            }

        }
        else if(str=="triangle")
        {
            tsi>>ntr;
            getline(tsi,str);
            for (int i=0;i<ntr;i++)
            {
                getline(tsi,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<4)
                {
                    std::cout<<"error ---> information of the triangles  "<<i<<" is not sufficent in the tsi file \n";
                    exit(0);
                }
                else
                {
                    Triangle_Map t;
                    t.id=f.String_to_Int(S.at(0));
                    t.v1=f.String_to_Int(S.at(1));
                    t.v2=f.String_to_Int(S.at(2));
                    t.v3=f.String_to_Int(S.at(3));
                    m_TriangleMap.push_back(t);

                }
            }
        }
        else if(str=="inclusion")
        {
            tsi>>ninc;
            getline(tsi,str);
            for (int i=0;i<ninc;i++)
            {
                getline(tsi,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<5)
                {
                    std::cout<<"error ---> information of the inclusion "<<i<<" is not sufficent in the tsi file \n";
                    exit(0);
                }
                else
                {
                    Inclusion_Map inc;
                    inc.id=f.String_to_Int(S.at(0));
                    inc.tid=f.String_to_Int(S.at(1));
                    inc.vid=f.String_to_Int(S.at(2));
                    double x=f.String_to_Double(S.at(3));
                    double y=f.String_to_Double(S.at(4));
                    double norm = sqrt(x*x+y*y);
                    x=x/norm;
                    y=y/norm;
                    inc.x = x; inc.y = y;
                    m_InclusionMap.push_back(inc);
                }
            }
        }
        else if(str=="exclusion")
        {
            tsi>>nexc;
            getline(tsi,str);
            for (int i=0;i<nexc;i++)
            {
                getline(tsi,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<3)
                {
                    std::cout<<"error ---> information of the exclusion at line "<<i<<" is not sufficent in the tsi file \n";
                    exit(0);
                }
                Exclusion_Map tem;
                tem.id = f.String_to_Int(S[0]);
                tem.vid =f.String_to_Int(S[1]);
                tem.R = f.String_to_Double(S[2]);
                m_ExclusionMap.push_back(tem);
            }

        }
        else
        {
            std::cout<<"error ---> "<<str<<" is unidentified key word for tsi file \n";
            exit(0);
        }
    }
}
void CreateMashBluePrint::Read_Mult_QFile(std::string topfile)
{
    Nfunction f;
    //== read the top file and store all the q files with the group name.
    std::vector<std::string> qfiles;
    std::vector<int> groupid;
    std::string str;
    int id;
    std::ifstream top;
    top.open(topfile.c_str());
    while (true)
    {
        top>>str>>id;
        if(top.eof())
            break;
        qfiles.push_back(str);
        groupid.push_back(id);
#if TEST_MODE == Enabled
        std::cout<<"----> q file in the top file: "<<str<<std::endl;
#endif
        getline(top,str);
    }
    top.close();
    // read each q file
    double Lx,Ly,Lz;
    int vid = 0;
    int tid = 0;
    for (int fi=0;fi<qfiles.size();fi++)
    {
        std::ifstream Qs;
        Qs.open((qfiles.at(fi)).c_str());

        // first line is the box size and it should only contain 3 numbers;
        getline(Qs,str);
        std::vector<std::string> b = f.split(str);
        if(b.size()>3)
        {
            std::cout<<"---> Error: box information in the file "<<qfiles.at(fi)<<" is not correct "<<std::endl;
            exit(0);
        }
        // The final box size will be the largest box in all the q files
        if(m_Box(0)<f.String_to_Double(b[0]))
            m_Box(0)=f.String_to_Double(b[0]);
        if(m_Box(1)<f.String_to_Double(b[1]))
            m_Box(1)=f.String_to_Double(b[1]);
        if(m_Box(2)<f.String_to_Double(b[2]))
            m_Box(2)=f.String_to_Double(b[2]);
#if DEBUG_MODE == Enabled
        std::cout<<"----> box was read  "<<std::endl;
        std::cout<<m_Box(0)<<"  "<<m_Box(1)<<"  "<<m_Box(2)<<"  "<<str<<std::endl;
#endif
        // reading the number of the vertices in this file
        getline(Qs,str);
        b.clear();
        b = f.split(str);
        if(b.size()>1)
        {
            std::cout<<"----> Error: number of vertices in the file "<<qfiles.at(fi)<<" is not correct "<<std::endl;
            exit(0);
        }
        int NV = f.String_to_Int(b[0]);
        for (int i=0;i<NV;i++)
        {
            getline(Qs,str);
            b.clear();
            b = f.split(str);
            if(b.size()>5 || b.size()<4)
            {
                std::cout<<"----> Error: Line "<<i+2<<", info of a vertex in the file "<<qfiles.at(fi)<<" is not correct.  "<<std::endl;
                exit(0);
            }
            Vertex_Map v;
            v.id=vid;
            v.include =true;
            v.x=f.String_to_Double(b[1]);
            v.y=f.String_to_Double(b[2]);
            v.z=f.String_to_Double(b[3]);
            if(b.size()==5)
            v.domain=f.String_to_Int(b[4]);
            else
            v.domain=0;
            m_VertexMap.push_back(v);
            vid++;
        }
#if DEBUG_MODE == Enabled
        std::cout<<"----> vertex section was read  "<<std::endl;
#endif
        getline(Qs,str);
        b.clear();
        b = f.split(str);
        if(b.size()>1)
        {
            std::cout<<"----> Error: number of triangle in the file "<<qfiles.at(fi)<<" is not correct "<<str<<std::endl;
            exit(0);
        }
        int nt=f.String_to_Int(b[0]);
        for (int i=0;i<nt;i++)
        {
            getline(Qs,str);
            b.clear();
            b = f.split(str);
            if(b.size()>5 || b.size()<4)
            {
                std::cout<<"----> Error: Line "<<i+2<<", info of a triangle in the file "<<qfiles.at(fi)<<" is not correct.  "<<std::endl;
                exit(0);
            }
            Triangle_Map t;
            t.id=tid;
            t.v1=f.String_to_Int(b[1]);
            t.v2=f.String_to_Int(b[2]);
            t.v3=f.String_to_Int(b[3]);
            m_TriangleMap.push_back(t);
            tid++;
        }
        Qs.close();
    }
    std::cout<<"triangle is read "<<"\n";

}
void CreateMashBluePrint::GenerateIncFromInputfile()
{

}
void CreateMashBluePrint::WriteCreateMashBluePrintLog()
{

}
