
#include "Argument.h"
#include "help.h"
#include "Nfunction.h"

Argument::Argument(std::vector <std::string> argument)
         :  m_Argument(argument),
            m_BondL(0.1),
            m_Seed(9474),
            m_Health(true),
            m_Monolayer(false),
            m_WPointDir(false),
            m_DTSFolder("point"),
            m_LipidLibrary("no"),
            m_GeneralOutputFilename("output"),
            m_InclusionDirectionType("Global"),
            m_Function("backmap"),
            m_RCutOff(0.5),
            m_Iter(5),
            m_Renorm(false),
            m_SkipLipids(false),
            m_KEEP_POINTS_CLOSE_TO_PROTEINS(false),
            m_PRINT_LESS_OUTPUT(false)
{


    m_SoftWareVersion = SoftWareVersion;
    m_ArgCon = 1;
    std::string Arg1;
    Nfunction f;

        for (long i=1;i<m_Argument.size();i=i+2)
        {
            Arg1 = m_Argument[i];
            if(Arg1 == G_POINT_FOLDER) //  "-dts"
            {
                if(m_Argument.size()>i+1){
                    m_DTSFolder = m_Argument[i+1];
                }
                else{
                    std::cout<<"---> error: 332 \n";
                    exit(0);
                }
            }
            else if (Arg1==FunctionType) // (Arg1=="-function")
            {
                if(m_Argument.size()>i+1){
                    m_Function = m_Argument.at(i+1);
                }
                if(m_Function=="1dsin")
                {
                    m_1DSinState.Lx = Nfunction::String_to_Double(m_Argument.at(i+2));
                    m_1DSinState.Ly = Nfunction::String_to_Double(m_Argument.at(i+3));
                    m_1DSinState.Lz = Nfunction::String_to_Double(m_Argument.at(i+4));
                    m_1DSinState.Omega = Nfunction::String_to_Double(m_Argument.at(i+5));
                    m_1DSinState.A = Nfunction::String_to_Double(m_Argument.at(i+6));
                    m_1DSinState.H = Nfunction::String_to_Double(m_Argument.at(i+7));
                    
                    m_1DSinState.APL = Nfunction::String_to_Double(m_Argument.at(i+8));
                    m_1DSinState.APW = Nfunction::String_to_Double(m_Argument.at(i+9));
                    i=i+8;

                }
            }
            else if(Arg1==G_DEFAULT_TAG)
            {
                if(m_Argument.size()>i+1){
                    m_GeneralOutputFilename = m_Argument.at(i+1);
                }
                else{
                    std::cout<<"---> error: 332 \n";
                    exit(0);
                }
            }
            else if(Arg1 == G_SKIP_LIPID_PLACEMENT)
            {
                
                m_SkipLipids = true;
                i=i-1;

            }
            else if(Arg1 == G_KEEP_POINTS_CLOSE_TO_PROTEINS)
            {
                
                m_KEEP_POINTS_CLOSE_TO_PROTEINS = true;
                i=i-1;

            }
            else if(Arg1 == G_PRINT_LESS_OUTPUTS)
            {
                m_PRINT_LESS_OUTPUT = true;
                i=i-1;

            }
            else if(Arg1 == G_HELPEx)
            {
                help helpmessage(m_Argument.at(0));
                exit(0);
            }
            else if(Arg1 == G_renormalized_lipid_ratio)
            {
                i=i-1;
                m_Renorm = true;

            }
            else if(Arg1=="-monolayer")
            {
                i=i-1;
                m_Monolayer = true;

            }
            else if(Arg1 == INC_DIR_TYPE)
            {
                m_InclusionDirectionType = m_Argument.at(i+1);
                if(m_InclusionDirectionType!="Local" || m_InclusionDirectionType!="Global")
                {
                    std::cout<<"Error: The inclusion direction type is unknown \n";
                }
                
            }
            else if(Arg1=="-iter")
            {
               m_Iter = f.String_to_Double(m_Argument.at(i+1));
                
            }
            else if(Arg1=="-Wall")
            {
                m_Wall.UpdateState(true);
                i=i-1;
            }
            else if(Arg1=="-WallDen")
            {
                m_Wall.UpdateDen(f.String_to_Double(m_Argument.at(i+1)));
                if(f.String_to_Double(m_Argument.at(i+1))>1)
                {
                    std::cout<<" Warning: the density of the wall beads is larger than 1: this has no effect and will act as 1 \n";
                }
            }
            else if(Arg1=="-WallBName")
            {
                m_Wall.UpdateBeadName(m_Argument.at(i+1));
            }
            else if(Arg1=="-WPointDir")
            {
                m_WPointDir = true;
                i=i-1;
            }
            else if(Arg1=="-WallH")
            {
                m_Wall.UpdateH(f.String_to_Double(m_Argument.at(i+1)));
            }
            else if(Arg1=="-WallBin")
            {
                m_Wall.UpdateCellSize(f.String_to_Double(m_Argument.at(i+1)));
            }
            else if(Arg1=="-WallUniform")
            {
                m_Wall.UpdateUniform(true);
                i=i-1;
            }
            else if(Arg1==LipidLibraryFileName)// "-LLIB"
            {
                m_LipidLibrary = m_Argument.at(i+1);
                
                if(m_LipidLibrary.substr(m_LipidLibrary.find_last_of(".") + 1) != LIBExt)
                {

                    m_LipidLibrary = m_LipidLibrary + "." +LIBExt;
                }
                if (f.FileExist (m_LipidLibrary)!=true)
                {
                    
                    std::cout<<"---> error: lipid library file, with name "<<m_LipidLibrary<<" does not exist \n";
                    m_Health = false;
                }
            }
            else if(Arg1 == G_STR_FILE_TAG)
            {
                m_StrFileName = m_Argument.at(i+1);
                if(m_StrFileName.substr(m_StrFileName.find_last_of(".") + 1) != STRExt)
                {
                    m_StrFileName = m_StrFileName + "." + STRExt;
                }
            }
            else if(Arg1 == GET_SEED) // "-seed"
            {
                m_Seed = f.String_to_Double(m_Argument.at(i+1));
            }
            else if(Arg1 == G_RADIUS_CUT_OFF) {
                m_RCutOff = f.String_to_Double(m_Argument.at(i+1));
            }
            else if(Arg1 == G_BOND_LENGTH) {
                m_BondL = f.String_to_Double(m_Argument.at(i+1));
            }
            else
            {
                std::cout << "---> error processing argument:"<<Arg1<<". \n";
                m_ArgCon=0;
                m_Health = false;
            }
        }
        
        /// checking if the defined files exists.
        if(m_Health == true)
        {
            if (f.FileExist (m_StrFileName)!=true)
            {
                std::cout<<"---> error, str file, with name . "<<m_StrFileName<<" . does not exist \n";
                m_Health = false;
            }
            if (f.FileExist (m_LipidLibrary)!=true)
            {
                std::cout<<" Note (warning): lipid library file, is not provided, some lipids may not exist  \n";
            }
        }
    
    std::ofstream log("pcg.log");  // Open the file directly in the constructor
    if (log.is_open())  // Ensure the file opened successfully
    {
        for (long i = 0; i < m_Argument.size(); i++)
        {
            log << m_Argument[i] << "  ";  // Use array-style indexing
        }
        log.close();  // Close the file after writing
    }
    else
    {
        std::cout << "Error: Could not open log file for writing.\n";
    }
    
    if(!m_Health){
        std::cout << "---> error: there was some errors in the commandline \n";
        std::cout<<"\n"<<"*** For more information and tips execute ./PCG -h ***"<<"\n";
        exit(0);
    }
    if(!ValidateVariables()){
        std::cout << "---> error: bad input data. \n";
        std::cout<<"\n"<<"*** For more information and tips execute ./PCG -h ***"<<"\n";
        exit(0);
    }
 
}
Argument::~Argument() {
   
}
bool Argument::ValidateVariables(){
    
    if(m_BondL <= 0 ){
        std::cout << "---> error: bond length should be positive and larger then zero.\n";
        return false;
    }
    if(m_Iter <= 0 ){
        std::cout << "---> error: iterations should be positive and larger then zero.\n";
        return false;
    }
    if(m_RCutOff <= 0 ){
        std::cout << "---> error: cutoff distance should be positive and larger then zero.\n";
        return false;
    }
    if(m_RCutOff <= 0 ){
        std::cout << "---> error: cutoff distance should be positive and larger then zero.\n";
        return false;
    }
    return true;
}
