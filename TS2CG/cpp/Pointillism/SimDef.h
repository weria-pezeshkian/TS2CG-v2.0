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


/// The keyword used for the commondline options, can be changed
#define Def_TSfile     "-TSfile"   // the name of the meshfile name
#define Def_HelpCall     "-h"   // help message
#define Def_BilayerThickness   "-bilayerThickness"
#define Def_AlgType             "-AlgType"
#define Def_rescalefactor       "-rescalefactor"
#define Def_AreaPerLipid        "-ap"
#define Def_DegreeOfMeshing        "-Mashno"
#define Def_TaskName                "-r"
#define Def_OutputFolderName        "-o"
#define Def_SmoothingFlag        "-smooth"
#define Def_resizebox           "-resizebox"
#define Def_resizeboxdist           "-b_dist"
#define Def_Monolayer           "-monolayer"
#define Def_PrintLessPutput           "-less"


#define KBT 1
#define PI 3.14159265359
#define S60 0.8660254037844
#define SQ3 1.73205080757
#define SoftWareVersion   "version 2.0"
#define Precision       8
#define Enabled   1
#define Disabled  2
#define Binary_Name     "PLM"
// commond line options name


#define G_tsiPrecision     "18.10"



#define Libarmadillo Disabled
#define EigenLib Enabled
#define NoLib Enabled
#define TEST_MODE  Disabled
#define RNGTYPE  UNIFROMTYPE0
#define BACKMAP  Enabled
