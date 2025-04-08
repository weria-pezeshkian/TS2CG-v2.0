#include <iostream>
#include "GroFile.h"
#include "Nfunction.h"
// a class that has functions to read and write gro files
GroFile::GroFile(std::string gmxfilename) {
    m_GroFileName = gmxfilename;
    ReadGroFile(m_GroFileName);
}

GroFile::~GroFile() {
}

void GroFile::RenewBeads(std::vector<bead> vB) {
    m_pAllBeads.clear();
    m_AllBeads.clear();
    m_AllBeads = vB;

    for (std::vector<bead>::iterator it = m_AllBeads.begin(); it != m_AllBeads.end(); ++it)
        m_pAllBeads.push_back(&(*it));
}

void GroFile::AddBead(bead b) {
    m_AllBeads.push_back(b);
}

void GroFile::UpdateBox(Vec3D box) {
    m_Box = box;
}

void GroFile::ReadGroFile(std::string file) {
    Nfunction f;
    char str[1000];
    if (file.size() < 4) {
        file = file + ".gro";
    } else if (file.at(file.size() - 1) == 'o' && file.at(file.size() - 2) == 'r' && file.at(file.size() - 3) == 'g') {
    } else {
        file = file + ".gro";
    }

    std::ifstream FGRO;
    FGRO.open(file.c_str());
    std::string str1;
    getline(FGRO, str1);
    getline(FGRO, str1);
    FILE *fgro;
    fgro = fopen(file.c_str(), "r");

    if (fgro == NULL) {
        printf(" Error: Could not open file %s", file.c_str());
    }

    bool check = fgets(str, 1000, fgro);
    m_Title = str;
    m_Title.pop_back();
    check = fgets(str, 1000, fgro);

    int NoBeads = atoi(str);

    float x, y, z, v1, v2, v3;
    char *A1;
    char *A2;

    char a[200];
    char b[200];
    int resid, beadid;
    std::string beadtype = "MDBeads";

    for (int i = 0; i < NoBeads; i++) {
        getline(FGRO, str1);
        std::vector<std::string> l = f.split(str1);

        int readafile = fscanf(fgro, "%d%s%s%d%f%f%f", &resid, a, b, &beadid, &x, &y, &z);
        check = fgets(str, 1000, fgro);

        if (l.size() == 7 || l.size() == 10) {
            x = atof((l.at(4)).c_str());
            y = atof((l.at(5)).c_str());
            z = atof((l.at(6)).c_str());
        } else if (l.size() == 6 || l.size() == 9) {
            x = atof((l.at(3)).c_str());
            y = atof((l.at(4)).c_str());
            z = atof((l.at(5)).c_str());
        } else if (l.size() == 5 || l.size() == 8) {
            x = atof((l.at(2)).c_str());
            y = atof((l.at(3)).c_str());
            z = atof((l.at(4)).c_str());
        } else {
            std::cout << "Warning: Perhaps error, something wrong with " << file << " file \n";
        }

        std::string bt = b;
        std::string beadname;
        if (bt.size() > 0)
            beadname.push_back(bt.at(0));
        if (bt.size() > 1)
            beadname.push_back(bt.at(1));

        std::string resname = a;
        bead Be(i, beadname, beadtype, resname, resid, x, y, z);
        m_AllBeads.push_back(Be);
    }

    FGRO.close();
    float Lx, Ly, Lz;
    int readafile = fscanf(fgro, "%f%f%f", &Lx, &Ly, &Lz);

    fclose(fgro);

    m_Box(0) = Lx;
    m_Box(1) = Ly;
    m_Box(2) = Lz;
    m_pBox = &m_Box;

    for (std::vector<bead>::iterator it = m_AllBeads.begin(); it != m_AllBeads.end(); ++it) {
        (*it).UpdateBox(m_pBox);
    }

    for (std::vector<bead>::iterator it = m_AllBeads.begin(); it != m_AllBeads.end(); ++it) {
        m_pAllBeads.push_back(&(*it));
    }
}

void GroFile::WriteGroFile(std::string file) {
    if (file.size() < 4) {
        file = file + ".gro";
    } else if (file.at(file.size() - 1) == 'o' && file.at(file.size() - 2) == 'r' && file.at(file.size() - 3) == 'g') {
    } else {
        file = file + ".gro";
    }

    FILE *fgro;
    fgro = fopen(file.c_str(), "w");

    if (!fgro) {
        std::cout << "Error opening file: " << file << std::endl;
        return;
    }

    const char *Title = "dmc gmx file handler";
    int Size = m_AllBeads.size();

    fprintf(fgro, "%s\n", Title);
    fprintf(fgro, "%5d\n", Size);

    int i = 0;
    for (std::vector<bead>::iterator it = m_AllBeads.begin(); it != m_AllBeads.end(); ++it) {
        i++;
        double x = (*it).GetXPos();
        double y = (*it).GetYPos();
        double z = (*it).GetZPos();
        i = i % 100000;
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n", i, ((*it).GetResName()).c_str(), ((*it).GetBeadName()).c_str(), i, x, y, z);
    }

    fprintf(fgro, "%10.5f%10.5f%10.5f\n", m_Box(0), m_Box(1), m_Box(2));
    fclose(fgro);
}
