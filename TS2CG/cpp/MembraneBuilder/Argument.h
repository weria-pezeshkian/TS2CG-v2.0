#if !defined(AFX_ARGUMENT_H_7F4A21B8_C13C_11D3_BF23_124095086234__INCLUDED_)
#define AFX_ARGUMENT_H_7F4A21B8_C13C_11D3_BF23_124095086234__INCLUDED_
#include "Def.h"
#include "Wall.h"
/**
 * @struct Shape_1DSin
 * @brief Structure representing the parameters for a 1D sine wave shape.
 *
 * This structure stores parameters used for generating 1D sine wave patterns,
 * which include dimensions, wave frequency (Omega), amplitude, and additional
 * properties such as APL (Area Per Lipid) and APW (Area Per Wave).
 */
struct Shape_1DSin {
    double Lx;   ///< Length in the x-dimension
    double Ly;   ///< Length in the y-dimension
    double Lz;   ///< Length in the z-dimension
    int Omega;   ///< Frequency of the wave
    double A;    ///< Amplitude of the sine wave
    double H;    ///< Height of the sine wave
    double APL;  ///< Area per lipid
    double APW;  ///< wall bead
};

/**
 * @class Argument
 * @brief Class to handle and process command-line arguments for the PCG program.
 *
 * The Argument class is responsible for parsing and validating command-line arguments,
 * setting configuration options, and providing access to these settings. It also checks
 * file existence for important inputs like structure files and lipid libraries.
 */
class Argument
{
public:
    /**
     * @brief Constructor that initializes the Argument object with command-line inputs.
     *
     * @param arg Vector of strings representing the command-line arguments.
     */
    Argument(std::vector<std::string> arg);

    /**
     * @brief Destructor for the Argument class.
     *
     * Cleans up any resources used, including closing log files or freeing allocated memory.
     */
    ~Argument();

    // Getters for various program settings and arguments
    inline const std::vector<std::string> GetArgumentString() const { return m_Argument; }
    inline const int GetArgCon() const { return m_ArgCon; }
    inline const std::string GetDTSFolder() const { return m_DTSFolder; }
    inline const std::string GetLipidLibrary() const { return m_LipidLibrary; }
    inline const std::string GetStructureFileName() const { return m_StrFileName; }
    inline const std::string GetInclusionDirectionType() const { return m_InclusionDirectionType; }
    inline const std::string GetGeneralOutputFilename() const { return m_GeneralOutputFilename; }
    inline const bool GetHealth() const { return m_Health; }
    inline const std::string GetSoftwareVersion() const { return m_SoftWareVersion; }
    inline const std::string GetFunction() const { return m_Function; }
    inline const int GetSeed() const { return m_Seed; }
    inline const double GetBond_length() const { return m_BondL; }
    inline const bool GetRenorm() const { return m_Renorm; }
    inline const double GetIter() const { return m_Iter; }
    inline const Wall GetWall() const { return m_Wall; }
    inline const double GetRCutOff() const { return m_RCutOff; }
    inline Shape_1DSin Get1DSinState() const { return m_1DSinState; }
    inline bool GetMonolayer() const { return m_Monolayer; }
    inline bool Skip_LipidPlacement() const { return m_SkipLipids; }

    bool m_WPointDir; ///< Flag for wall point direction, public to allow direct modification

private:
    std::vector<std::string> m_Argument; ///< Vector storing command-line arguments
    std::string m_DTSFolder;             ///< DTS folder name
    std::string m_InclusionDirectionType; ///< Inclusion direction type (Local/Global)
    std::string m_StrFileName;           ///< Structure file name
    std::string m_LipidLibrary;          ///< Lipid library file name
    std::string m_GeneralOutputFilename; ///< General output filename
    std::string m_SoftWareVersion;       ///< Software version string
    std::string m_Function;              ///< Function type (e.g., "backmap", "1dsin")
    
    int m_ArgCon;                        ///< Argument condition flag (for error checking)
    bool m_Health;                       ///< Health flag to check if inputs are valid
    int m_Seed;                          ///< Seed for random number generation
    double m_BondL;                      ///< Bond length
    bool m_Renorm;                       ///< Renormalization flag
    bool m_Monolayer;                    ///< Monolayer flag
    double m_Iter;                       ///< Number of iterations for the algorithm
    double m_RCutOff;                    ///< Cutoff distance for interactions
    bool m_SkipLipids;                    ///if true, do not place any lipid,

    Wall m_Wall;                         ///< Wall object storing wall-related data and settings
    Shape_1DSin m_1DSinState;            ///< Shape configuration for the 1D sine wave
};

#endif // !defined(AFX_ARGUMENT_H_7F4A21B8_C13C_11D3_BF23_124095086234__INCLUDED_)
