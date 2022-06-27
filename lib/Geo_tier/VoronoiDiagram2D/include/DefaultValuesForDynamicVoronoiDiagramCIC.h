#ifndef _DEFAULT_VALUES_FOR_DYNAMIC_VORONOI_DIAGRAM_CIC_
#define _DEFAULT_VALUES_FOR_DYNAMIC_VORONOI_DIAGRAM_CIC_

#include "rg_Const.h"
#include "SimpleTypesForDynamicVoronoiDiagramCIC.h"

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

}
}

#endif