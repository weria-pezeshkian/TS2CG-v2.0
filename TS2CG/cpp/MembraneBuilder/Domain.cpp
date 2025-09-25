  #if !defined(AFX_Domain_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_Domain_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "Domain.h"
#include "GenerateUnitCells.h"
#include "Def.h"

/*
 

 
 
 
 */

Domain::Domain(int domaintypeid, std::vector<point*>  point1)
{
    m_point = point1;
    m_DomainTypeID = domaintypeid;
    m_DomainTotalLipid = 0;
}
Domain::~Domain()
{
    
}
void Domain::AddADomainLipid(std::string name, double Ap, double Ratio)
{
    DomainLipid DL;
    DL.Name = name;
    DL.Ap=Ap;
    DL.Ratio = Ratio;
    m_AllDomainLipids.push_back(DL);
}
void Domain::Configure(bool renorm)
{
  //== this function 1) obtain the total ratio and possibly make it 1. 2) obtains max lipids of each one in this domain

    
    // obtain the total ratio
    double totalratio=0;
    for ( std::vector<DomainLipid>::iterator it = m_AllDomainLipids.begin(); it != m_AllDomainLipids.end(); it++ )
        totalratio+=(*it).Ratio ;

    if (renorm) {
        std::cout << "--> You requested renormalization of the lipid ratio. "
                  << "The new values are:\n";

        if (totalratio != 0.0) {
            for (auto &lipid : m_AllDomainLipids) {
                lipid.Ratio /= totalratio;
                std::cout << "--> " << lipid.Name << "   " << lipid.Ratio << '\n';
            }
            totalratio = 1.0;
        } else {
            std::cerr << "Error: total ratio is zero, cannot renormalize.\n";
        }
    }
    else if (std::abs(totalratio - 1.0) > 1e-6) {
        std::cerr << "Warning: the total lipid percentage for domain "
                  << m_DomainTypeID << " is " << totalratio
                  << " (not 1). Make sure this is intentional, "
                  << "or use the -renorm option.\n";
    }
    //== obtaining the total area of the domain points
    
    double Tarea = 0; // total area of the domain point;
    for ( std::vector<point*>::iterator it = m_point.begin(); it != m_point.end(); it++ )
    {
        Tarea+= (*it)->GetArea();
    }
    Tarea = Tarea*totalratio;
    double avAP = 0;
    // update the max no lipid of each lipid type in the domain
    //**************
    //************  Ratio_new_lipid = ratio/(ratio*a1+ratio*a2+......)
    //******************
    for ( std::vector<DomainLipid>::iterator it = m_AllDomainLipids.begin(); it != m_AllDomainLipids.end(); it++ )
    {
        avAP+= ((*it).Ap)*((*it).Ratio);
    }
    for ( std::vector<DomainLipid>::iterator it = m_AllDomainLipids.begin(); it != m_AllDomainLipids.end(); it++ )
    {
        if(avAP!=0)
        (*it).MaxNo= int(((*it).Ratio)*Tarea/avAP);
        else
        (*it).MaxNo=0;
        (*it).no_created = 0;
    }
    
    m_DomainTotalLipid = 0;
    for ( std::vector<DomainLipid>::iterator it = m_AllDomainLipids.begin(); it != m_AllDomainLipids.end(); it++ )
    { m_DomainTotalLipid+= (*it).MaxNo;
       // std::cout<<"max domain lipid "<<m_DomainTotalLipid<<"  "<<(*it).MaxNo<<"\n";
    }
    
    if(m_DomainTotalLipid>m_point.size())
    {
        std::cout<<" Error: Not enough point for the domain to place lipid: "<<m_DomainTotalLipid<<"  lipid should be created  "<<m_point.size()<<" points is available \n";
        std::exit(0);
    }
    for ( std::vector<DomainLipid>::iterator it = m_AllDomainLipids.begin(); it != m_AllDomainLipids.end(); it++ )
    {
        m_pAllDomainLipids.push_back(&(*it));
    }
}







#endif



