#if !defined(AFX_Domain_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_Domain_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <math.h>
#include <list>
#include <map>
#include <iomanip>
#include <valarray>
#include "LMatrix.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "point.h"
#include "Data_Structure.h"
/*
-------------------------------------------------------------------------------
 Class: Domain

 Purpose:
   Represents a lipid domain, which consists of:
     - A unique domain type ID
     - A set of spatial points belonging to the domain
     - A list of lipid species in the domain
     - Bookkeeping information (total lipid count, pointers for fast access)

 Responsibilities:
   - Store the set of points that make up the domain.
   - Configure lipid distribution within the domain by renormalizing ratios
   - Provide accessors for domain data (points, lipids, IDs).

 Key Methods:
   - AddADomainLipid(): Add a lipid species with a given area and ratio.
   - Configure():       Renormalize lipid ratios (if requested) and
                        calculate how many lipids of each type can fit.

 Notes:
   - Does not use C++11 features to maintain compatibility with older compilers.
   - Inline getters provide fast read-only access to members.
-------------------------------------------------------------------------------
*/

class Domain
{
public:
    
	Domain(int domaintypeid, std::vector<point*>  point);
	virtual ~Domain();
    
    
    
    // Inline getters
    inline  std::vector<point*>            GetDomainPoint()                const  {return m_point;}
    inline  std::vector<DomainLipid>       GetDomainLipids()               const  {return m_AllDomainLipids;}
    inline  std::vector<DomainLipid*>      GetpDomainLipids()              const  {return m_pAllDomainLipids;}
    inline  int                            GetDomainID()                   const  {return m_DomainTypeID;}
    inline  int                            GetDomainTotalLipid()           const  {return m_DomainTotalLipid;}


    // Public methods
public:
    void AddADomainLipid(std::string name, double Ap, double Ratio);
    void Configure(bool);
    
private:
    int m_DomainTypeID;
    std::vector<point*>  m_point;
    std::vector<DomainLipid>  m_AllDomainLipids;
    std::vector<DomainLipid*>  m_pAllDomainLipids;
    int m_DomainTotalLipid;
};

#endif // AFX_Domain_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_
