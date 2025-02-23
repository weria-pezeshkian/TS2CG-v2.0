#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <list>
#include <map>
#include <iomanip>
#include <valarray>
#include <stdio.h>

#define PI 3.14159265359
#define S60 0.8660254037844
#define SQ3 1.73205080757
#define Precision               8
#define Enabled                 1
#define Disabled                2
//==== version and file extenions 
#define SoftWareVersion                 "2.0"
#define ExecutableName                  "PCG"
#define SoftWareName                    "CG Membrane builder"
#define LIBExt                          "LIB"
#define STRExt                          "str"


//===== command line key tages
#define G_STR_FILE_TAG                  "str"
#define FunctionType                    "-function"
#define G_HELPEx                        "-h"
#define LipidLibraryFileName            "-LLIB"
#define GET_SEED                        "-seed"
#define INC_DIR_TYPE                    "-incdirtype"
#define G_SKIP_LIPID_PLACEMENT            "-skip_lipids"
#define G_POINT_FOLDER                  "-dts"
#define G_DEFAULT_TAG                   "-defout"
#define G_BOND_LENGTH                   "-Bondlength"
#define G_RADIUS_CUT_OFF                   "-Rcutoff"
#define G_renormalized_lipid_ratio          "-renorm"
#define G_KEEP_POINTS_CLOSE_TO_PROTEINS          "-keep"
#define G_PRINT_LESS_OUTPUTS                    "-less"



#define TESTMODE Disabled
#define EigenLib Disabled
#define NoLib Enabled
#define GroRead2024 Enabled
