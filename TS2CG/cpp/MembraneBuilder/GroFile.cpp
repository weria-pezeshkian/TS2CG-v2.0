

#include <stdio.h>
#include <algorithm>
#include "GroFile.h"
#include "Nfunction.h"

GroFile::GroFile(std::string gmxfilename) {
    m_GroFileName = gmxfilename;
    ReadGroFile(m_GroFileName);

}


GroFile::~GroFile() {
    
}

void GroFile::AddBead(bead b) {

    m_AllBeads.push_back(b);
}


#if GroRead2024 == Enabled

void GroFile::ReadGroFile(std::string file) {
    
    if (Nfunction::SubstringFromRight(file,'.') != "gro"){
        file = file + ".gro";
    }
    
    std::ifstream FGRO(file.c_str());
    if (!FGRO.is_open()) {
        std::cout << "---> Error: Could not open the file: " << file << std::endl;
        exit(0);  // Return on error instead of throwing an exception
    }
    getline (FGRO,m_Title);
    m_Title.erase(std::remove(m_Title.begin(), m_Title.end(), ' '), m_Title.end());
    
    std::string str,line;
    int NoBeads;
    FGRO>>NoBeads;
    getline (FGRO,str);

    double x,y,z,v1,v2,v3;
    int resid, secondInt;
    std::string firstStr, secondStr;
        
    for (int i=0; i<NoBeads; i++) {
        
        getline (FGRO,line);
        resid = std::stoi(line.substr(0, 5));  // First integer (column 0-4)
        firstStr = line.substr(5, 5);             // First string (column 5-9)
        secondStr = line.substr(10, 5);           // Second string (column 10-14)
        secondInt = std::stoi(line.substr(15, 5)); // Second integer (column 15-19)
        firstStr.erase(std::remove(firstStr.begin(), firstStr.end(), ' '), firstStr.end());
        secondStr.erase(std::remove(secondStr.begin(), secondStr.end(), ' '), secondStr.end());
        std::istringstream iss(line.substr(20));
        iss >> x >> y >> z >> v1 >> v2 >> v3;
        
        bead make_Bead(i, secondStr, secondStr, firstStr, resid, x, y, z);
        m_AllBeads.push_back(make_Bead);
        
    }
    float Lx,Ly,Lz;
    FGRO>>Lx>>Ly>>Lz;
    FGRO.close();
    
    m_Box(0)=Lx; m_Box(1)=Ly; m_Box(2)=Lz;
    m_pBox = &m_Box;
    double xcm =0;
    double ycm =0;
    double zcm =0;
    
    double no_beads = double(m_AllBeads.size());
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it) {
        (*it).UpdateBox(m_pBox);
        xcm += (*it).GetXPos()/no_beads;
        ycm += (*it).GetYPos()/no_beads;
        zcm += (*it).GetZPos()/no_beads;
    }
    
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it) {
        (*it).UpdateBox(m_pBox);
        (*it).UpdateXPos((*it).GetXPos()-xcm);
        (*it).UpdateYPos((*it).GetYPos()-ycm);
        (*it).UpdateZPos((*it).GetZPos()-zcm);
    }
    
    
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it) {
        m_pAllBeads.push_back(&(*it));
    }
    
    return;
}
#elif
void GroFile::ReadGroFile(std::string file)
{
    Nfunction f;
    char str[1000];
    if(file.size()<4)
    {
        file=file+".gro";
    }
    else if(file.at(file.size()-1)=='o' && file.at(file.size()-2)=='r' && file.at(file.size()-3)=='g')
    {
        
    }
    else
    {
        file=file+".gro";
    }
    
    
    std::ifstream FGRO;
    FGRO.open(file.c_str());
    std::string str1;
    getline (FGRO,str1);
    getline (FGRO,str1);
    FILE *fgro;
    fgro = fopen(file.c_str(), "r");
    
    if (fgro == NULL){
        printf(" Error: Could not open file %s",file.c_str());
    }
    bool check = fgets(str, 1000, fgro);
    m_Title = str;
    m_Title.pop_back();
    check = fgets(str, 1000, fgro);
    
    int NoBeads = atoi(str);
    
    // Reading the atoms from gro file
    
    float x,y,z,v1,v2,v3;
    char *A1;
    char *A2;
    
    char a[200];
    char b[200];
    int resid,beadid;
    std::string beadtype = "MDBeads";
    
    for (int i=0; i<NoBeads; i++) //NoBeads
    {
        
        getline (FGRO,str1);
        std::vector <std::string> l=f.split(str1);
        
        int readafile = fscanf(fgro, "%d%s%s%d%f%f%f",&resid,a,b,&beadid,&x,&y,&z);
        check = fgets(str, 1000, fgro);
        
        if(l.size()==7 || l.size()==10)
        {
            x=atof((l.at(4)).c_str());
            y=atof((l.at(5)).c_str());
            z=atof((l.at(6)).c_str());
        }
        else if(l.size()==6 || l.size()==9)
        {
            x=atof((l.at(3)).c_str());
            y=atof((l.at(4)).c_str());
            z=atof((l.at(5)).c_str());
        }
        else if(l.size()==5 || l.size()==8)
        {
            x=atof((l.at(2)).c_str());
            y=atof((l.at(3)).c_str());
            z=atof((l.at(4)).c_str());
        }
        else
        {
            std::cout<<"Warning: Perhaps error, something wrong with"<<file <<"file \n";
        }
        std::string bt = b;
        std::string beadname;
        if(bt.size()>0)
            beadname.push_back(bt.at(0));
        if(bt.size()>1)
            beadname.push_back(bt.at(1));
        
        std::string resname = a;
        bead Be(i, beadname, beadtype, resname, resid, x, y, z);
        m_AllBeads.push_back(Be);
        
        //  std::cout<<resid<<"  "<<resname<<"  "<<beadname<<" "<<beadid<<" "<<x<<"  "<<y<<"  "<<z<<"  \n";
        
        
    }
    FGRO.close();
    float Lx,Ly,Lz;
    int readafile = fscanf(fgro, "%f%f%f",&Lx,&Ly,&Lz);
    // std::cout<<Lx<<" "<<Ly<<" "<<Lz<<" \n";
    
    fclose(fgro);
    
    m_Box(0)=Lx; m_Box(1)=Ly; m_Box(2)=Lz;
    m_pBox = &m_Box;
    double xcm =0;
    double ycm =0;
    double zcm =0;
    
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        (*it).UpdateBox(m_pBox);
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos();
        
        xcm+=x/double(m_AllBeads.size());
        ycm+=y/double(m_AllBeads.size());
        zcm+=z/double(m_AllBeads.size());
        
        
        
    }
    
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        (*it).UpdateBox(m_pBox);
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos();
        (*it).UpdateXPos(x-xcm);
        (*it).UpdateYPos(y-ycm);
        (*it).UpdateZPos(z-zcm);
        
        
    }
    
    
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        m_pAllBeads.push_back(&(*it));
        
    }
    
}
#endif

/*void GroFile::ReadGroFile(std::string file)
{
    Nfunction f;
    if(file.size()<4)
    {
        file=file+".gro";
    }
    else if(file.at(file.size()-1)=='o' && file.at(file.size()-2)=='r' && file.at(file.size()-3)=='g')
    {
        
    }
    else
    {
        file=file+".gro";
    }
    
    std::string str;
    std::ifstream FGRO;
    FGRO.open(file.c_str());
    std::string title;
    getline (FGRO,title);
    m_Title = title;
    int atomno;
    FGRO>>atomno;
    getline (FGRO,str);
    std::string resname,aname;
    for (int i=0; i<atomno; i++) //NoBeads
    {
        char *cr  = new char [5];
        FGRO.read(cr,5);
        int resid = atoi(cr);
        FGRO.read(cr,5);
        resname = cr;
        resname.erase(remove_if(resname.begin(), resname.end(), isspace), resname.end());
        FGRO.read(cr,5);
        aname = cr;
        aname.erase(remove_if(aname.begin(), aname.end(), isspace), aname.end());
        FGRO.read(cr,5);
        int atomno = atoi(cr);
        char *Xr  = new char [8];
        FGRO.read(Xr,8);
        double X=atof(Xr);
        FGRO.read(Xr,8);
        double Y=atof(Xr);
        FGRO.read(Xr,8);
        double Z=atof(Xr);
        getline (FGRO,title);
        
        bead Be(i, aname, aname, resname, resid, X, Y, Z);
        m_AllBeads.push_back(Be);
        
    }
    float Lx,Ly,Lz;
    char *Lr  = new char [10];
    FGRO.read(Lr,10);
    Lx=atof(Lr);
    FGRO.read(Lr,10);
    Ly=atof(Lr);
    FGRO.read(Lr,10);
    Lz=atof(Lr);
    FGRO.close();
    
    m_Box(0)=Lx; m_Box(1)=Ly; m_Box(2)=Lz;
    m_pBox = &m_Box;
    double xcm =0;
    double ycm =0;
    double zcm =0;

    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        (*it).UpdateBox(m_pBox);
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos();
        
        xcm+=x/double(m_AllBeads.size());
        ycm+=y/double(m_AllBeads.size());
        zcm+=z/double(m_AllBeads.size());


        
    }
    

    
    
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        m_pAllBeads.push_back(&(*it));
  
    }
    
}
 */

void GroFile::WriteGroFile(std::string file) {

    if (Nfunction::SubstringFromRight(file,'.') != "gro"){
        file = file + ".gro";
    }
    
    FILE *fgro;
    fgro = fopen(file.c_str(), "w");
    
    
    /// resid  res name   noatom   x   y   z
    const char* Title="dmc gmx file handler";
    int Size=m_AllBeads.size();
    
    fprintf(fgro,  "%s\n",Title);
    fprintf(fgro, "%5d\n",Size);
    int i=0;
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        i++;
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos();
        int resid=(*it).GetResid();
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid,((*it).GetResName()).c_str(),((*it).GetBeadName()).c_str(),i,x,y,z );

    }
    

    fprintf(fgro,  "%10.5f%10.5f%10.5f\n",m_Box(0),m_Box(1),m_Box(2) );
    fclose(fgro);
    
    return;
}


