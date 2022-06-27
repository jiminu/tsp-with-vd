#ifndef _CONSTFORVORONOIDIAGRAM3D_H
#define _CONSTFORVORONOIDIAGRAM3D_H

#include "rg_Const.h"
#include "rg_Point3D.h"

namespace V {

namespace GeometryTier {


typedef rg_Point3D rg_Vector3D;

const rg_INT  NUM_DEFINING_CELLS_OF_ONE_VERTEX = 4;
const rg_INT  NUM_INCIDENT_EDGES_OF_ONE_VERTEX = 4;
const rg_INT  NUM_INCIDENT_FACES_OF_ONE_EDGE   = 3;

const rg_INT NUM_CELLS_TO_DEFINE_EDGE = 3;
const rg_INT NUM_FACES_TO_SHARE_EDGE  = 3;
const rg_INT NUM_EDGES_TO_BE_PROPAGATED_BY_NEW_VERTEX = 3;

const rg_INT NUM_INITIAL_EDGES   = 4;
const rg_INT NUM_INITIAL_FACES   = 6;
const rg_INT NUM_PARTIALEDGE_PER_EDGE = 3;

// RATIO_OF_BOXFORVD_TO_BOXFORVERTEX defines the size of bounding box for voronoi diagram 
// to represent voronoi vertices in the infinite region. The bounding box for voronoi diagram 
// is RATIO_OF_BOXFORVD_TO_BOXFORVERTEX times as big as that for voronoi vertices.
const rg_REAL RATIO_OF_BOXFORVD_TO_BOXFORVERTEX = 2.0;


// The following constant defines the type of voronoi edge to be obtained three spheres.
const rg_FLAG LINEAR_TYPE_EDGE    = 0;
const rg_FLAG NONLINEAR_TYPE_EDGE = 1;

const rg_FLAG NON_EDGE                     = 0;
const rg_FLAG PARABOLIC_OR_HYPERBOLIC_EDGE = 1;
const rg_FLAG CIRCULAR_OR_ELLIPTIC_EDGE    = 2;
const rg_FLAG LINEAR_EDGE                  = 3;
const rg_FLAG EDGE_BY_VERTEX_DEGENERACY    = 4;
const rg_FLAG EDGE_WITHOUT_VERTICES        = 5;

const rg_FLAG NONELLIPTIC_EDGE = 1;
const rg_FLAG ELLIPTIC_EDGE    = 2;


const rg_FLAG NOT_FOUND                 = 0;
const rg_FLAG FRONT_FEASIBLE_REGION     = 1;
const rg_FLAG FRONTREAR_FEASIBLE_REGION = 2;
const rg_FLAG REAR_FEASIBLE_REGION      = 3;


const rg_FLAG NORMAL_STATE             = 1;
const rg_FLAG NOT_FOUND_INITIAL_VERTEX = 2;
const rg_FLAG VERTEX_DEGREE_VIOLATION  = 3;

//const rg_REAL SYMBOLIC_PERTURBATION = H_RADIUS*0.1;

//const int numOfEvalPtsOnCrv = 5;
//const int numOfGridPts = 11;

enum PositionOnEdge {START, END, MID, ALL};

const rg_FLAG FULLY_TANGIBLE     = 2;
const rg_FLAG PARTIALLY_TANGIBLE = 1;
const rg_FLAG NON_TANGIBLE       = 0;

const rg_INT  ERROR_INFINITE_LOOP = 100;

const rg_INT  NUM_POINTS_ON_EDGE = 50;


typedef struct {
            char*    m_version;
            const char*    m_pdbFilename;
            const char*    m_date;
            char*    m_author;
            char*    m_organization;
            char*    m_operatingSystem;
            rg_INT   m_computingAlgorithm;  //  1: edge_tracing,  2: region_expansion
            char*    m_OS;
            char*    m_CPU;
            char*    m_RAM;
            rg_REAL  m_computingTime;
        } HeaderOfVDSFile;


} // namespace GeometryTier

} // namespace V


#endif


