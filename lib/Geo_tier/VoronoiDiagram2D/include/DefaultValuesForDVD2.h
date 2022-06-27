#ifndef _DEFAULT_VALUES_FOR_DVD2_
#define _DEFAULT_VALUES_FOR_DVD2_

#include "rg_Const.h"
#include "SimpleTypesForDVD2.h"


namespace V {
namespace GeometryTier {

//const PolynomialSolverType POLYNOMIAL_SOLVER = STURMS_SEQ;
//const PolynomialSolverType POLYNOMIAL_SOLVER = JENKINS_TRAUB;
const PolynomialSolverType POLYNOMIAL_SOLVER = PolynomialSolverType::LABSRC;

/*******************           USED TOLERANCE IN DVD                ******************/
const double TOLERANCE_OF_EDGE_FLIP_ROOT                       = resNeg4;
const double TOLERANCE_OF_GAVRILOVA_INEQUALITY                 = resNeg4;
const double TOLERANCE_OF_SHADOW_FLIP                          = resNeg4;
const double TOLERANCE_OF_DISK_INTERSECTION                    = resNeg6;
const double TOLERANCE_OF_EQUAL_RADIUS                         = resNeg4;
const double TOLERANCE_OF_COLLISION_PREDICTION                 = resNeg6;
const double TOLERANCE_OF_EQUAL_VELOCITY_VEC                   = resNeg6;
const double TOLERANCE_OF_TOP_ELEMENT_WITH_IDENTICAL_KEY_IN_PQ = resNeg4;

const int    MAX_ITERATION_OF_FALSE_POSITION                   = 100;
const double TOLERANCE_OF_FALSE_POSITION                       = resNeg12;
/*************************************************************************************/

/*******************           USED TOLERANCE IN DVDSIM             ******************/
const double TOLERANCE_OF_ZERO_VELOCITY = resNeg6;
//const double TOLERANCE_OF_DISK_INTERSECTION = resNeg6;
/*************************************************************************************/


/******************* USED CHARACTER FOR PARSING EVENTSEQ IN DVDSIM  ******************/
const char EDGE_FLIP        = 'E';
const char SHADOW_EDGE_FLIP = 'S';
const char COLLISION        = 'C';
const char VELOCITY_CHANGE  = 'V';
const char TELEPORT_BTW_CONTAINERS = 'T';
const char BOUNCE_SHADOW_TELPORT   = 'B';
/*************************************************************************************/








const string DVD_2D_VERSION = "VDRC_DVD_OF_DISKS_IN_2D_1.0";

const double INITIAL_EVENT_SEQUENCE_FILE_VERSION = 1.0; // initial version

const double EVENT_SEQUENCE_FILE_FORMAT_VERSION = 1.1; // delete shadow flip, add collision pos., chage to .evtsq





}
}

#endif