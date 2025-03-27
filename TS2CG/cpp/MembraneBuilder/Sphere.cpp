  #if !defined(AFX_Sphere_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_Sphere_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "Sphere.h"
#include "GroFile.h"
#include "GenerateUnitCells.h"
#include "Def.h"
#include "PDBFile.h"
#include "GenDomains.h"
Sphere::Sphere(Argument *pArgu)
{
    m_WallBox.push_back(0);
    m_WallBox.push_back(1);
    m_WallBox.push_back(0);
    m_WallBox.push_back(1);
    m_DL = 0.0;
    m_Box(0) = 15;m_Box(1) = 15;m_Box(2) = 15;
    m_R = 5;
    m_pBox = &m_Box;
    m_Thickness = 4;
    m_Density = 2.5;
    m_WallDensityup = 1;
    m_WallDensityin = 1;
    std::string ifilename = pArgu->GetStructureFileName();
    Initialize(ifilename);

    //=============== We need the H wall distance, and the only info needed about the wall here
    Wall CWall = pArgu->GetWall();
    double Hwall = CWall.GetH();

    //*********

    m_Point1 = CalculateArea_MakePoints(1, 1/m_Density,m_Thickness/2,m_DL,false);
    m_Point2 = CalculateArea_MakePoints(-1, 1/m_Density,m_Thickness/2,m_DL,false);

    
    m_WallPoint1 = CalculateArea_MakePoints(1, 1/m_WallDensityup,m_Thickness/2,Hwall,true);
    m_WallPoint2 = CalculateArea_MakePoints(-1, 1/m_WallDensityin,m_Thickness/2, Hwall,true);

}
Sphere::~Sphere()
{
    
}
std::vector<point> Sphere::CalculateArea_MakePoints(int layer, double APL,double H, double DL, bool wall)
{
    std::vector<point> Cpoints;
    double pi = acos(-1);
    double TotalArea = 0;
    std::vector <double> Curv;
    
    if(layer==1)
    {
    Curv.push_back(1.0/(m_R+H));
    Curv.push_back(1.0/(m_R+H));
        TotalArea = 4*pi*(m_R+H)*(m_R+H);
    }
    else if(layer==-1)
    {
        Curv.push_back(1.0/(m_R-H));
        Curv.push_back(1.0/(m_R-H));
        TotalArea = 4*pi*(m_R-H)*(m_R-H);

    }
    else
    {
        std::cout<<"error 308 \n";
    }

    int Npoints = TotalArea/APL;
    APL =  TotalArea/double(Npoints);
    int beadid = 0;
    double DT = FindDeltaTheta(Npoints);
    double T=0;
    double R=m_R+double(layer)*H;
    for (int i=0;i<=PI/DT;i++)
    {
        T+=DT;
        int M=2*PI*sin(T)/DT;
        
        
        for (int j=0;j<M;j++)
        {
            double Phi=DT/sin(T)*j;
            double x=R*cos(Phi)*sin(T);
            double y=R*sin(Phi)*sin(T);
            double z=R*cos(T);
            
            Vec3D Pos(x,y,z);
            Vec3D N=Pos*(double(layer)/Pos.norm());
            Vec3D BoxC(m_Box(0)/2, m_Box(1)/2, m_Box(2)/2 );
            Pos=BoxC+Pos+N*(DL);
            
            Vec3D P1(-N(1),N(0),0);
            Vec3D P2 = N*P1;
            
            
            point p(beadid, APL, Pos, N, P1, P2 , Curv);
            Cpoints.push_back(p);
            beadid++;
            

        }
        
    }
   // std::cout<<Npoints<<"  "<<Cpoints.size()<<"\n";
    for (std::vector<point >::iterator it = Cpoints.begin() ; it != Cpoints.end(); ++it)
    {
            it->UpdateArea(TotalArea/double(Cpoints.size()));
    }
    return Cpoints;

}
double Sphere::FindDeltaTheta(int N)
{
    
    double deltaTheta=PI/2;
    int TotPoint=0;
    int s=2;
    
    while(true)
    {
        int M=PI/deltaTheta;
        TotPoint =0;
        for (int i=0;i<M;i++)
            TotPoint+=2*PI*sin(i*deltaTheta)/deltaTheta;
        
        if(TotPoint>N)
            break;
        else
        {
            s++;
            deltaTheta=PI/s;
        }
        
    }

    int Q = PI/deltaTheta;
    deltaTheta = PI/double(Q);
    
    return deltaTheta;
}
void Sphere::Initialize(std::string filename)
{
    
    bool OK=true;
    Nfunction f;
    std::ifstream file;
    file.open(filename.c_str());
    bool flag = false;
    std::string str;

    while (true)
    {
        std::getline (file,str);
        if(file.eof())
            break;
        
        std::vector<std::string> Line = f.split(str);
        if(Line.size()!=0 && (Line.at(0)).at(0)!=';')
        {
            if((Line.at(0)).at(0)=='[' && flag==false)
            {
                str = f.trim(str);
                str.erase(std::remove(str.begin(), str.end(), '['), str.end());
                str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
                str.erase(std::remove(str.begin(), str.end(), ' '), str.end());

                if(str=="ShapeData")
                    flag = true;
            }
            else if((Line.at(0))=="End" && flag==true)
            {
                flag=false;
            }
            else if(flag==true)
            {
                if(Line.at(0)=="Box")
                {
                    if(Line.size()<4)
                        std::cout<<" Error: Box information in the str file is not correct \n";
                    else
                    {
                        m_Box(0) = f.String_to_Double(Line.at(1));
                        m_Box(1) = f.String_to_Double(Line.at(2));
                        m_Box(2) = f.String_to_Double(Line.at(3));

                    }
                }
                else if(Line.at(0)=="ShapeType")
                {
                }
                else if(Line.at(0)=="Density")
                {
                    if(Line.size()<2)
                    {
                        std::cout<<" Error: Density information in the str file is not correct \n";
                    }
                    else
                    {
                        m_Density = f.String_to_Double(Line.at(1));
                    }
                }
                else if(Line.at(0)=="WallDensity")
                {
                    if(Line.size()<2)
                    {
                        std::cout<<" Error: Density information in the str file is not correct \n";
                    }
                    else
                    {
                        m_WallDensityup = f.String_to_Double(Line.at(1));
                        m_WallDensityin = f.String_to_Double(Line.at(2));

                    }
                }
                else if(Line.at(0)=="Thickness")
                {
                    if(Line.size()<2)
                        std::cout<<" Error: Thickness information in the str file is not correct \n";
                    else
                    {
                        m_Thickness = f.String_to_Double(Line.at(1));
                        
                    }
                }
                else if(Line.at(0)=="DL")
                {
                    if(Line.size()<2)
                        std::cout<<" Error: DL information in the str file is not correct \n";
                    else
                    {
                        m_DL = f.String_to_Double(Line.at(1));
                        
                    }
                }
                else if(Line.at(0)=="Radius")
                {
                    if(Line.size()<2)
                        std::cout<<" Error: Radius information in the str file is not correct \n";
                    else
                    {
                        m_R = f.String_to_Double(Line.at(1));
                        
                    }
                }

            }


        }
            
            
        }
        
    
    file.close();

    
}






#endif



