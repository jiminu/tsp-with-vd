#ifndef _CONST_FOR_SHORTEST_PATH_FINDER_
#define _CONST_FOR_SHORTEST_PATH_FINDER_

#include "rg_Const.h"
#include <string>
using namespace std;

enum SolverForBarrierProblemType {TPT_METHOD, EXPAND_PATH_ALONG_BARRIER_METHOD};

const int INPUT_MAX_NUM_OF_CANDIDATE_PATHS_BY_TP_TREE = 5;
const SolverForBarrierProblemType solverForBarrierProblemType = EXPAND_PATH_ALONG_BARRIER_METHOD;

enum ShortestPathFinder2DAlgorithmStep
{
    SPF_TIME_TOTAL,
    SPF_TIME_CONSTRUCT_VD_FAMILY,
    SPF_TIME_FILTER_QT_FACE_BY_ELLIPSE,
    SPF_TIME_FIND_PATH_BY_BFS,
    SPF_TIME_FIND_PATH_BY_DFS,
    SPF_TIME_FIND_PATH_BY_ALL_PATH_SEARCH,
    SPF_TIME_SIZE_PHASE_ONE,
};


#endif