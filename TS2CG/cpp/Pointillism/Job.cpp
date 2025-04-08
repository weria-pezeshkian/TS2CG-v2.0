#include "Job.h"
#include "Nfunction.h"
#include "Edit_configuration.h"

Job::Job(std::vector <std::string> argument)
{
    /*
     *  Copyright Weria Pezeshkian (weria.pezeshkian@gmail.com), updated 2023.
     * This is a class to call different functions: at the moment, only one exist
     */
     
        std::string exacutable=argument.at(0);
        if(TrimNameFromPath(exacutable) == Binary_Name)
        {
            Edit_configuration B(argument);
        }
        else
        {
            std::cout<<" ---> error: unrecognized executable binary name: "<<argument.at(0)<<"\n";
        }
}
Job::~Job()
{
    
}
std::string Job::TrimNameFromPath(std::string &s)
{
    std::string bname;
    
    for (int i=s.size()-1;i>=0;i--)
    {
        if(s[i]=='/')
            break;
        else
            bname.insert(bname.begin(), 1, s[i]);
    }
    return bname;
}
