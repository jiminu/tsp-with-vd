#include "BetaFace.h"

#include "BetaVertex.h"
#include "BetaEdge.h"
#include "BetaCell.h"

#include "rg_RelativeOp.h"

#include "FunctionsForVoronoiDiagram3D.h"

#include "VDEdge.h"
using namespace V::GeometryTier;



BetaFace::BetaFace()
: m_leftCell(rg_NULL), m_rightCell(rg_NULL),
  m_visited(rg_FALSE)
{
    for ( int i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_edge[i]        = rg_NULL;
        m_orientation[i] = rg_TRUE;
    }

    m_vEdge = rg_NULL;
}



BetaFace::BetaFace(const rg_INT& ID)
: TopologicalEntity(ID),
  m_leftCell(rg_NULL), m_rightCell(rg_NULL),
  m_visited(rg_FALSE)
{
    for ( int i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_edge[i]        = rg_NULL;
        m_orientation[i] = rg_TRUE;
    }

    m_vEdge = rg_NULL;
}



BetaFace::BetaFace(BetaEdge** edge)
: m_leftCell(rg_NULL), m_rightCell(rg_NULL),
  m_visited(rg_FALSE)
{
    for ( int i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_edge[i]        = edge[i];
        m_orientation[i] = rg_TRUE;
    }

    m_vEdge = rg_NULL;
}



BetaFace::BetaFace(BetaCell* leftCell, BetaCell* rightCell)
: m_leftCell(leftCell), m_rightCell(rightCell),
  m_visited(rg_FALSE)
{
    for ( int i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_edge[i]        = rg_NULL;
        m_orientation[i] = rg_TRUE;
    }

    m_vEdge = rg_NULL;
}



BetaFace::BetaFace(const rg_INT& ID, BetaCell* leftCell, BetaCell* rightCell)
: TopologicalEntity(ID),
  m_leftCell(leftCell), m_rightCell(rightCell),
  m_visited(rg_FALSE)
{
    for ( int i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_edge[i]        = rg_NULL;
        m_orientation[i] = rg_TRUE;
    }

    m_vEdge = rg_NULL;
}



BetaFace::BetaFace(const BetaFace& betaFace)
: TopologicalEntity(betaFace.m_ID),
  m_leftCell(betaFace.m_leftCell), m_rightCell(betaFace.m_rightCell),
  m_betaSpan(betaFace.m_betaSpan),
  m_visited(betaFace.m_visited)
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_edge[i]        = betaFace.m_edge[i];
        m_orientation[i] = betaFace.m_orientation[i];
    }

    m_vEdge = betaFace.m_vEdge;
}



BetaFace::~BetaFace()
{
    m_leftCell    = rg_NULL;
    m_rightCell   = rg_NULL;
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_edge[i] = rg_NULL;
    }

    //if ( m_vEdge != rg_NULL ) {
    //    m_vEdge->disconnectBetaFace(this);
    //}
}



rg_BOOL BetaFace::isVirtual() const
{
    rg_BOOL isVirtualFace = rg_FALSE;
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        BetaVertex* vertex = rg_NULL;
        if ( m_orientation[i] )  {
            vertex = m_edge[i]->getEndVertex();
        }
        else  {
            vertex = m_edge[i]->getStartVertex();
        }
        
        if ( vertex->isVirtual() ) {
            isVirtualFace = rg_TRUE;
            break;
        }
    }

    return isVirtualFace;
}



BetaCell*  BetaFace::getLeftCell() const
{
    return m_leftCell;
}



BetaCell*  BetaFace::getRightCell() const
{
    return m_rightCell;
}



BetaCell*  BetaFace::getOppositeCell(BetaCell* incidentCell) const
{
    BetaCell* oppositeCell = rg_NULL;
    if ( m_leftCell == incidentCell ) {
        oppositeCell = m_rightCell;
    }
    else if ( m_rightCell == incidentCell ) {
        oppositeCell = m_leftCell;
    }

    return oppositeCell;
}



BetaEdge** BetaFace::getEdges()
{
    return m_edge;
}



BetaEdge*  BetaFace::getEdge(const rg_INT& i) const
{
    if ( i >= 0 && i<EIWDS_NUM_EDGE_ON_FACE ) {
        return m_edge[i];
    }
    else {
        return rg_NULL;
    }
}



rg_INT BetaFace::getPosOfEdge(BetaEdge* edge) const
{
    for ( rg_INT i_edge=0; i_edge<EIWDS_NUM_EDGE_ON_FACE; i_edge++ )  {
        if ( edge == m_edge[i_edge] )  
            return i_edge;
    }

    return rg_UNKNOWN;
}



rg_BOOL* BetaFace::getEdgeOrientations()
{
    return m_orientation;
}



rg_BOOL BetaFace::getEdgeOrientation(const rg_INT& i) const
{
    if ( i >= 0 && i<EIWDS_NUM_EDGE_ON_FACE ) {
        return m_orientation[i];
    }
    else {
        return rg_UNKNOWN;
    }
}



rg_BOOL BetaFace::getEdgeOrientation(const BetaEdge* const edge) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ )  {
        if ( edge == m_edge[i] )  
            return m_orientation[i];
    }

    return rg_UNKNOWN;
}



rg_BetaSpan  BetaFace::getBetaSpan() const
{
    return m_betaSpan;
}



rg_REAL BetaFace::getMinValueForValidSimplex() const
{
    return m_betaSpan.getLowerValueOfInterval(1);
}



rg_INT      BetaFace::getDepth() const
{
    rg_INT depth = rg_UNKNOWN;

    rg_dList<BetaVertex*> vertexOnThisFace;
    searchVerticesInWholeWorld( vertexOnThisFace );

    vertexOnThisFace.reset4Loop();
    while ( vertexOnThisFace.setNext4Loop() ) {
        BetaVertex* currVtx = vertexOnThisFace.getEntity();
        if ( depth < currVtx->getDepth() ) {
            depth = currVtx->getDepth();
        }
    }

    return depth;
}



void BetaFace::setLeftCell(BetaCell*  leftCell)
{
    m_leftCell = leftCell;
}



void BetaFace::setRightCell(BetaCell* rightCell)
{
    m_rightCell = rightCell;
}



void BetaFace::setEdges(BetaEdge** edges)
{
    for ( int i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_edge[i] = edges[i];
    }
}



void BetaFace::setEdge(const rg_INT& i, BetaEdge* edge)
{
    if ( i >= 0 && i<EIWDS_NUM_EDGE_ON_FACE ) {
        m_edge[i] = edge;
    }
}



void BetaFace::setEdge(const rg_INT& i, const rg_BOOL& edgeOrientation, BetaEdge* edge)
{
    if ( i >= 0 && i<EIWDS_NUM_EDGE_ON_FACE ) {
        m_edge[i]        = edge;
        m_orientation[i] = edgeOrientation;
    }
}



void BetaFace::setEdgeOrientations(rg_BOOL* edgeOrientation)
{
    for ( int i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_orientation[i] = edgeOrientation[i];
    }
}



void BetaFace::setEdgeOrientation(const rg_INT& i, const rg_BOOL& edgeOrientation)
{
    if ( i >= 0 && i<EIWDS_NUM_EDGE_ON_FACE ) {
        m_orientation[i] = edgeOrientation;
    }
}



void BetaFace::addEdge(BetaEdge* edge)
{
	for( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ )  {
		if( m_edge[i] == rg_NULL )  {
			m_edge[i] = edge;
			break;
		}
	}
}



void BetaFace::setBetaSpan(const rg_BetaSpan& betaSpan)
{
    m_betaSpan = betaSpan;
}

void BetaFace::shiftBetaSpan(const rg_REAL& delta)
{
	m_betaSpan.shiftBetaSpan(delta);
}


BetaFace& BetaFace::operator=(const BetaFace& betaFace)
{
    if ( this == &betaFace ) {
        return *this;
    }

    m_ID        = betaFace.m_ID;
    m_leftCell  = betaFace.m_leftCell;
    m_rightCell = betaFace.m_rightCell;
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        m_edge[i]        = betaFace.m_edge[i];
        m_orientation[i] = betaFace.m_orientation[i];
    }
    m_betaSpan  = betaFace.m_betaSpan;

    m_visited   = betaFace.m_visited;

    m_vEdge = betaFace.m_vEdge;

    return *this;
}





rg_BOOL BetaFace::isIsolatedFace() const
{
    if ( m_leftCell == rg_NULL && m_rightCell == rg_NULL ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_BOOL BetaFace::isOnConvexHull() const
{
    if ( m_leftCell == rg_NULL || m_rightCell == rg_NULL ) {
        return rg_FALSE;
    }
    else {
        if (    (!m_leftCell->isVirtual() &&  m_rightCell->isVirtual())
             || ( m_leftCell->isVirtual() && !m_rightCell->isVirtual()) ) {
            return rg_TRUE;
        }
        else {
            return rg_FALSE;
        }
    }
}



rg_BOOL BetaFace::isSingularEllipticFace(const rg_REAL& beta) const
{
    if ( isIsolatedFace() ) {
        return rg_FALSE;
    }
    else {
        if ( this->getBoundingState(beta) == SINGULAR_SIMPLEX ) {
            if (    m_leftCell->getBoundingState(beta) == INTERIOR_SIMPLEX 
                && m_rightCell->getBoundingState(beta) == INTERIOR_SIMPLEX ) {
                return rg_TRUE;
            }
        }
    }

    return rg_FALSE;
}


rg_BOOL   BetaFace::isIn2AdjacencyAnomaly() const
{
    if ( isIsolatedFace() ) {
        return rg_FALSE;
    }
    else {
        rg_INT    multiAdjacency = m_leftCell->isMultiplicity();
        BetaCell* neighbor       = m_leftCell->getNeighborCellWithMultiplicity();

        if ( multiAdjacency == 2 && neighbor == m_rightCell) {
            return rg_TRUE;
        }
        else {
            return rg_FALSE;
        }    
    }
}



rg_BOOL   BetaFace::isIn3AdjacencyAnomaly() const
{
    if ( isIsolatedFace() ) {
        return rg_FALSE;
    }
    else {
        rg_INT    multiAdjacency = m_leftCell->isMultiplicity();
        BetaCell* neighbor       = m_leftCell->getNeighborCellWithMultiplicity();

        if ( multiAdjacency == 3 && neighbor == m_rightCell) {
            return rg_TRUE;
        }
        else {
            return rg_FALSE;
        }    
    }
}



rg_BOOL   BetaFace::isIn4AdjacencyAnomaly() const
{
    if ( isIsolatedFace() ) {
        return rg_FALSE;
    }
    else {
        rg_INT    multiAdjacency = m_leftCell->isMultiplicity();
        BetaCell* neighbor       = m_leftCell->getNeighborCellWithMultiplicity();

        if ( multiAdjacency == 4 && neighbor == m_rightCell) {
            return rg_TRUE;
        }
        else {
            return rg_FALSE;
        }    
    }
}



rg_BOOL BetaFace::isThere(const BetaVertex* const vertex) const
{
    rg_BOOL isThereVertex = rg_FALSE;
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        if ( m_orientation[i] )  {
            if ( vertex == m_edge[i]->getEndVertex() ) {
                isThereVertex =  rg_TRUE;
                break;
            }
        }
        else  {
            if ( vertex == m_edge[i]->getStartVertex() ) {
                isThereVertex =  rg_TRUE;
                break;
            }
        }
    }

    return isThereVertex;

//     if ( vertex == m_edge[0]->getStartVertex() )
//         return rg_TRUE;
// 
//     if ( vertex == m_edge[0]->getEndVertex() )
//         return rg_TRUE;
// 
//     if ( m_orientation[1] == rg_TRUE )  {
//         if ( vertex == m_edge[1]->getEndVertex() )
//             return rg_TRUE;
//         else 
//             return rg_FALSE;
//     }
//     else  {
//         if ( vertex == m_edge[1]->getStartVertex() )
//             return rg_TRUE;
//         else 
//             return rg_FALSE;
//     }
}



rg_BOOL BetaFace::isThere(const BetaEdge* const   edge)   const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ )  {
        if ( edge == m_edge[i] )  {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_BOOL BetaFace::isAdjacent(const BetaFace* const   face)   const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ )  {
        if ( face->isThere( m_edge[i] ) ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}




BetaEdge* BetaFace::findSharingEdgeWith(const BetaFace* const adjacentFace) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ )  {
        if ( adjacentFace->isThere( m_edge[i] ) ) {
            return m_edge[i];
        }
    }

    return rg_NULL;
}



BetaEdge* BetaFace::findEdge(const BetaVertex* const vertex1, const BetaVertex* const vertex2) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ )  {
        if ( m_edge[i] == rg_NULL ) {
            continue;
        }
        else {
            if (    m_edge[i]->getStartVertex() == vertex1 
                 && m_edge[i]->getEndVertex()   == vertex2   ) {
                return m_edge[i];
            }
            else if (    m_edge[i]->getStartVertex() == vertex2 
                      && m_edge[i]->getEndVertex()   == vertex1   ) {
                return m_edge[i];
            }
            else {
            }
        }
    }

    return rg_NULL;
}



void BetaFace::findEdgeIncidentToVertex(const BetaVertex* const vertex, BetaEdge** edge) const
{
    rg_INT i_edge = 0;
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ )  {
        if ( vertex == m_edge[i]->getStartVertex() || vertex == m_edge[i]->getEndVertex() ) {
            edge[i_edge] = m_edge[i];
            i_edge++;
        }
    }
}


    
BetaEdge*   BetaFace::findMateEdgeOfVertex(const BetaVertex* const vertex) const
{
    BetaEdge* mateEdge = rg_NULL;

    if ( isThere(vertex) ) {
        for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ )  {
            if ( vertex != m_edge[i]->getStartVertex() && vertex != m_edge[i]->getEndVertex() ) {
                mateEdge = m_edge[i];
                break;
            }
        }
    }

    return mateEdge;
}



BetaVertex* BetaFace::findMateVertexOfEdge(const BetaEdge* const edge) const
{
    rg_INT posOfNextEdge = rg_NULL;
    if ( m_edge[0] == edge ) {
        posOfNextEdge = 1;
    }
    else if ( m_edge[1] == edge ) {
        posOfNextEdge = 2;
    }
    else if ( m_edge[2] == edge ) {
        posOfNextEdge = 0;
    }
    else {
        return rg_NULL;
    }

    if ( m_orientation[posOfNextEdge] ) {
        return m_edge[posOfNextEdge]->getEndVertex();
    }
    else {
        return m_edge[posOfNextEdge]->getStartVertex();
    }
}


    
BetaFace* BetaFace::findFaceIncidentToEdgeInCell(const BetaEdge* const edge, const BetaCell* const cell) const
{
    BetaFace*   faceIncidentToEdgeInCell = rg_NULL;

    BetaVertex* mateVertex = findMateVertexOfEdge(edge);
    if ( cell == m_leftCell || cell == m_rightCell ) {
        faceIncidentToEdgeInCell = cell->getMateFace( mateVertex );
    }

    return faceIncidentToEdgeInCell;
}



void BetaFace::computeBetaSpan()
{
    BetaVertex* vertex[3];
    searchVerticesInIntraWorld( vertex );

    Sphere  tangentSphere[2];
    rg_INT  numTangentSphere = computeMinTangentSpheresTo3SpheresOutside_GLOBAL(
                                  vertex[0]->getBall(), vertex[1]->getBall(), vertex[2]->getBall(), tangentSphere );
                               
    rg_REAL minOfCurrSimplex = tangentSphere[0].getRadius();
    rg_REAL maxOfCurrSimplex = tangentSphere[1].getRadius();


    if ( isIsolatedFace() ) {    
        m_betaSpan.setNumBetaInterVal(3);
        m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
        m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX, minOfCurrSimplex, maxOfCurrSimplex);
        m_betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX, maxOfCurrSimplex, REAL_PLUS_INFINITE);
    }
    else {
        if ( numTangentSphere == 1 ) {
            computeBetaSpanOfNonEllipticFace(tangentSphere[0]);
        }
        else if ( numTangentSphere == 2 ) {
            computeBetaSpanOfEllipticFace(tangentSphere[0], tangentSphere[1]);
        }
        else {
            // do nothing.
        }

//         rg_REAL minOfHigherSimplex = REAL_PLUS_INFINITE;
//         rg_REAL maxOfHigherSimplex = REAL_MINUS_INFINITE;
// 
//         rg_dList<BetaCell*> finiteCellList;
//         searchFiniteCellsInWholeWorld( finiteCellList );
//         finiteCellList.reset4Loop();
//         while ( finiteCellList.setNext4Loop() ) {
//             rg_REAL radius = finiteCellList.getEntity()->getMinTangentSphere().getRadius();
// 
//             if ( radius < minOfHigherSimplex ) {
//                 minOfHigherSimplex = radius;
//             }
// 
//             if ( radius > maxOfHigherSimplex ) {
//                 maxOfHigherSimplex = radius;
//             }
//         }
// 
// 
//         rg_FLAG isEllipticFaceAttached                  = rg_UNKNOWN;
//         rg_BOOL isMaxOfCurrSimGreaterThanMaxOfHigherSim = rg_FALSE;
//         //  The attachedness for a q-face corresponding to circular or elliptic v-edge
//         //  is checked as the following block.
//         //  if ( numTangentSphere == 2 ) { ... }
//         //  The block also checks for whether the center of min. tangent sphere, tangentSphere[0],
//         //  of q-face is placed on the corresponding v-edge or not.
// 
// 
//         //  In the corresponding primary structure(sphere Voronoi diagram),
//         //  three balls define the elliptic or circular Voronoi edge.
//         if ( numTangentSphere == 2 ) {
//             // When three balls define an elliptic Voronoi edge, 
//             // tangentSphere[0] and tangentSphere[1] are respectively the minimum and maximum spheres tangent to three balls.
//             rg_Point3D startVoronoiVtx   = m_leftCell->getMinTangentSphere().getCenter();
//             rg_Point3D endVoronoiVtx     = m_rightCell->getMinTangentSphere().getCenter();
// 
//             Plane centerPlane;
//             centerPlane.definePlaneByThreePoints( vertex[0]->getBall().getCenter(), 
//                                                   vertex[1]->getBall().getCenter(), 
//                                                   vertex[2]->getBall().getCenter()  );
//             rg_Point3D vectorMaxToMin    = tangentSphere[0].getCenter() - tangentSphere[1].getCenter();
//             rg_Point3D normalOfEdgePlane = vectorMaxToMin.crossProduct( centerPlane.getNormal() );
//             normalOfEdgePlane.normalize();
// 
//             rg_Point3D axisPt = (tangentSphere[0].getCenter() + tangentSphere[1].getCenter())*0.5;
//             Plane      edgePlane(normalOfEdgePlane, axisPt);
// 
//             rg_REAL    angleFromStartToEnd = edgePlane.computeAngleInCCW(startVoronoiVtx, axisPt, endVoronoiVtx);
//             rg_REAL    angleFromStartToMax = edgePlane.computeAngleInCCW(startVoronoiVtx, axisPt, tangentSphere[1].getCenter());
//             rg_REAL    angleFromStartToMin = edgePlane.computeAngleInCCW(startVoronoiVtx, axisPt, tangentSphere[0].getCenter());
// 
//             //  check for whether the center of tangentSphere[1] is placed on the v-edge or not.
//             if ( angleFromStartToEnd > angleFromStartToMax ) {
//                 //maxOfHigherSimplex = maxOfCurrSimplex;
//                 isMaxOfCurrSimGreaterThanMaxOfHigherSim = rg_TRUE;
//             }
// 
//             //  check for whether the center of tangentSphere[0] is placed on the v-edge or not.
//             if ( angleFromStartToEnd > angleFromStartToMin ) {
//                 isEllipticFaceAttached = rg_FALSE;
//             }
//             else {
//                 isEllipticFaceAttached = rg_TRUE;
//             }
//         }
// 
// 
// 
//         if ( isOnConvexHull() ) {
//             if ( isAttached(isEllipticFaceAttached, tangentSphere[0]) ) {
//                 m_betaSpan.setNumBetaInterVal(2);
//                 m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
//                 m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  REAL_PLUS_INFINITE);
//             }
//             else {
//                 m_betaSpan.setNumBetaInterVal(3);
//                 m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
//                 m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
//                 m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  REAL_PLUS_INFINITE);
//             }
//         }
//         else {
//             if ( isAttached(isEllipticFaceAttached, tangentSphere[0]) ) {
//                 if ( isMaxOfCurrSimGreaterThanMaxOfHigherSim ) {
//                     m_betaSpan.setNumBetaInterVal(4);
//                     m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
//                     m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
//                     m_betaSpan.setBetaInterval(2, SINGULAR_SIMPLEX,   maxOfHigherSimplex,  maxOfCurrSimplex);
//                     m_betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxOfCurrSimplex,    REAL_PLUS_INFINITE);
//                 }
//                 else {
//                     m_betaSpan.setNumBetaInterVal(3);
//                     m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
//                     m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
//                     m_betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);
//                 }
//             }
//             else {
//                 if ( isMaxOfCurrSimGreaterThanMaxOfHigherSim ) {
//                     m_betaSpan.setNumBetaInterVal(5);
//                     m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
//                     m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
//                     m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
//                     m_betaSpan.setBetaInterval(3, SINGULAR_SIMPLEX,   maxOfHigherSimplex,  maxOfCurrSimplex);
//                     m_betaSpan.setBetaInterval(4, INTERIOR_SIMPLEX,   maxOfCurrSimplex,    REAL_PLUS_INFINITE);
//                 }
//                 else {
//                     m_betaSpan.setNumBetaInterVal(4);
//                     m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
//                     m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
//                     m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
//                     m_betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);
//                 }
//             }
//         }
    }
}



void BetaFace::computeBetaSpanOfNonEllipticFace(const Sphere& minTangentSphere)
{
    rg_REAL minOfCurrSimplex = minTangentSphere.getRadius();

    rg_REAL minOfHigherSimplex = REAL_PLUS_INFINITE;
    rg_REAL maxOfHigherSimplex = REAL_MINUS_INFINITE;

    BetaCell* incidentCell[2] = { m_leftCell, m_rightCell };
    for ( rg_INT i=0; i<2; i++ ) {
        if ( incidentCell[i]->isVirtual() ) {
            continue;
        }

        rg_REAL radius = incidentCell[i]->getMinTangentSphere().getRadius();

        if ( radius < minOfHigherSimplex ) {
            minOfHigherSimplex = radius;
        }

        if ( radius > maxOfHigherSimplex ) {
            maxOfHigherSimplex = radius;
        }
    }

    //rg_dList<BetaCell*> finiteCellList;
    //searchFiniteCellsInWholeWorld( finiteCellList );
    //finiteCellList.reset4Loop();
    //while ( finiteCellList.setNext4Loop() ) {
    //    rg_REAL radius = finiteCellList.getEntity()->getMinTangentSphere().getRadius();

    //    if ( radius < minOfHigherSimplex ) {
    //        minOfHigherSimplex = radius;
    //    }

    //    if ( radius > maxOfHigherSimplex ) {
    //        maxOfHigherSimplex = radius;
    //    }
    //}


    if ( isOnConvexHull() ) {
        if ( isNonEllipticFaceAttached() ) {
            m_betaSpan.setNumBetaInterVal(2);
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
            m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  REAL_PLUS_INFINITE);
        }
        else {
            m_betaSpan.setNumBetaInterVal(3);
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
            m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  REAL_PLUS_INFINITE);
        }
    }
    else {
        if ( isNonEllipticFaceAttached() ) {
            m_betaSpan.setNumBetaInterVal(3);
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
            m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);
        }
        else {
            m_betaSpan.setNumBetaInterVal(4);
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
            m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
            m_betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);
        }
    }
}



void BetaFace::computeBetaSpanOfEllipticFace(const Sphere& minTangentSphere, const Sphere& maxTangentSphere)
{
    rg_REAL minOfCurrSimplex   = minTangentSphere.getRadius();
    rg_REAL maxOfCurrSimplex   = maxTangentSphere.getRadius();


    rg_REAL minOfHigherSimplex = REAL_PLUS_INFINITE;
    rg_REAL maxOfHigherSimplex = REAL_MINUS_INFINITE;

    BetaCell* incidentCell[2] = { m_leftCell, m_rightCell };
    for ( rg_INT i=0; i<2; i++ ) {
        if ( incidentCell[i]->isVirtual() ) {
            continue;
        }

        rg_REAL radius = incidentCell[i]->getMinTangentSphere().getRadius();

        if ( radius < minOfHigherSimplex ) {
            minOfHigherSimplex = radius;
        }

        if ( radius > maxOfHigherSimplex ) {
            maxOfHigherSimplex = radius;
        }
    }

    //rg_dList<BetaCell*> finiteCellList;
    //searchFiniteCellsInWholeWorld( finiteCellList );
    //finiteCellList.reset4Loop();
    //while ( finiteCellList.setNext4Loop() ) {
    //    rg_REAL radius = finiteCellList.getEntity()->getMinTangentSphere().getRadius();

    //    if ( radius < minOfHigherSimplex ) {
    //        minOfHigherSimplex = radius;
    //    }

    //    if ( radius > maxOfHigherSimplex ) {
    //        maxOfHigherSimplex = radius;
    //    }
    //}



    rg_FLAG isEllipticFaceAttached                  = rg_UNKNOWN;
    rg_BOOL isMaxOfCurrSimGreaterThanMaxOfHigherSim = rg_FALSE;

    // When three balls define an elliptic Voronoi edge, 
    // tangentSphere[0] and tangentSphere[1] are respectively the minimum and maximum spheres tangent to three balls.
    rg_Point3D startVoronoiVtx   = m_leftCell->getMinTangentSphere().getCenter();
    rg_Point3D endVoronoiVtx     = m_rightCell->getMinTangentSphere().getCenter();

    BetaVertex* vertex[3];
    searchVerticesInIntraWorld( vertex );

    Plane centerPlane;
    centerPlane.definePlaneByThreePoints( vertex[0]->getBall().getCenter(), 
                                          vertex[1]->getBall().getCenter(), 
                                          vertex[2]->getBall().getCenter()  );

    rg_Point3D vectorMaxToMin    = minTangentSphere.getCenter() - maxTangentSphere.getCenter();
    rg_Point3D normalOfEdgePlane = vectorMaxToMin.crossProduct( centerPlane.getNormal() );
    normalOfEdgePlane.normalize();

    rg_Point3D axisPt = (minTangentSphere.getCenter() + maxTangentSphere.getCenter())*0.5;
    Plane      edgePlane(normalOfEdgePlane, axisPt);

    rg_REAL    angleFromStartToEnd = edgePlane.computeAngleInCCW(startVoronoiVtx, axisPt, endVoronoiVtx);
    rg_REAL    angleFromStartToMax = edgePlane.computeAngleInCCW(startVoronoiVtx, axisPt, maxTangentSphere.getCenter());
    rg_REAL    angleFromStartToMin = edgePlane.computeAngleInCCW(startVoronoiVtx, axisPt, minTangentSphere.getCenter());

    //  check for whether the center of maxTangentSphere is placed on the v-edge or not.
    if ( angleFromStartToEnd > angleFromStartToMax ) {
        isMaxOfCurrSimGreaterThanMaxOfHigherSim = rg_TRUE;
    }

    //  check for whether the center of minTangentSphere is placed on the v-edge or not.
    if ( angleFromStartToEnd > angleFromStartToMin ) {
        isEllipticFaceAttached = rg_FALSE;
    }
    else {
        isEllipticFaceAttached = rg_TRUE;
    }


    if ( isEllipticFaceAttached ) {
        if ( isMaxOfCurrSimGreaterThanMaxOfHigherSim ) {
            m_betaSpan.setNumBetaInterVal(4);
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
            m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, SINGULAR_SIMPLEX,   maxOfHigherSimplex,  maxOfCurrSimplex);
            m_betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxOfCurrSimplex,    REAL_PLUS_INFINITE);
        }
        else {
            m_betaSpan.setNumBetaInterVal(3);
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
            m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);
        }
    }
    else {
        if ( isMaxOfCurrSimGreaterThanMaxOfHigherSim ) {
            m_betaSpan.setNumBetaInterVal(5);
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
            m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
            m_betaSpan.setBetaInterval(3, SINGULAR_SIMPLEX,   maxOfHigherSimplex,  maxOfCurrSimplex);
            m_betaSpan.setBetaInterval(4, INTERIOR_SIMPLEX,   maxOfCurrSimplex,    REAL_PLUS_INFINITE);
        }
        else {
            m_betaSpan.setNumBetaInterVal(4);
            m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
            m_betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
            m_betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
            m_betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);
        }
    }
}



rg_BOOL BetaFace::isAttached(const rg_FLAG& isEllipticFaceAttached, const Sphere& minTangentSphere) const
{
    if ( isEllipticFaceAttached == rg_UNKNOWN ) {
        //  This attachedness check is for a q-face corresponding to linear, parabolic, or hyperbolic v-edge.
        rg_BOOL isAttachedOrNot = rg_FALSE;

        BetaVertex* vertex[3];
        searchVerticesInIntraWorld( vertex );

        Plane centerPlane;
        centerPlane.definePlaneByThreePoints( vertex[0]->getBall().getCenter(), 
                                              vertex[1]->getBall().getCenter(), 
                                              vertex[2]->getBall().getCenter()  );

        if ( m_leftCell->isVirtual() ) {
            rg_Point3D endVoronoiVtx     = m_rightCell->getMinTangentSphere().getCenter();

			rg_REAL    dist = centerPlane.distanceFromPoint( endVoronoiVtx );
            if ( rg_NEG( dist ) ) {
                isAttachedOrNot = rg_TRUE;
            }
        }
        else if ( m_rightCell->isVirtual() ) {
            rg_Point3D startVoronoiVtx   = m_leftCell->getMinTangentSphere().getCenter();

			rg_REAL    dist = centerPlane.distanceFromPoint( startVoronoiVtx );
            if ( rg_POS( dist ) ) {
                isAttachedOrNot = rg_TRUE;
            }
        }
        else {
            rg_Point3D startVoronoiVtx   = m_leftCell->getMinTangentSphere().getCenter();
            rg_Point3D endVoronoiVtx     = m_rightCell->getMinTangentSphere().getCenter();

			rg_REAL    distToStart = centerPlane.distanceFromPoint( startVoronoiVtx );
			rg_REAL    distToEnd   = centerPlane.distanceFromPoint( endVoronoiVtx );
            
            if ( rg_NEG( distToStart ) && rg_POS( distToEnd ) ) {
                isAttachedOrNot = rg_FALSE;
            }
            else {
                isAttachedOrNot = rg_TRUE;
            }

        }
//         Sphere ballInVicinity[2] = { m_leftCell->getMateVertex(this)->getBall(), 
//                                      m_rightCell->getMateVertex(this)->getBall() };
// 
//         rg_BOOL isAttachedOrNot = rg_FALSE;
//         for ( rg_INT i=0; i<2; i++ ) {
//             if ( minTangentSphere.isThereIntersectionWith( ballInVicinity[i]) ) {
//                 isAttachedOrNot = rg_TRUE;
//                 break;
//             }
//         }
// 
//     
//         if ( !isAttachedOrNot ) {
//             rg_INT multiplicity[2] = { m_leftCell->isMultiplicity(), m_rightCell->isMultiplicity() };
// 
//             if ( multiplicity[0] == 3 && multiplicity[1] == 3 ) {
//                 if ( m_leftCell != m_rightCell->getNeighborCellWithMultiplicity() ) {
//                     rg_REAL signedVolumeOfCellInVicinity[2] = { m_leftCell->computeSignedVolume(), 
//                                                                 m_rightCell->computeSignedVolume() };
// 
//                     if ( signedVolumeOfCellInVicinity[0]<0.0 || signedVolumeOfCellInVicinity[1]<0.0) {
//                         isAttachedOrNot = rg_TRUE;
//                     }
//                 }
//             }
//             else if ( multiplicity[0] == 3 && multiplicity[1] != 3 ) {
//                 if ( m_leftCell->computeSignedVolume() < 0.0 ) {
//                     isAttachedOrNot = rg_TRUE;
//                 }
//             }
//             else if ( multiplicity[0] != 3 && multiplicity[1] == 3 ) {
//                 if ( m_rightCell->computeSignedVolume() < 0.0 ) {
//                     isAttachedOrNot = rg_TRUE;
//                 }
//             }
//             else {
//                 // do nothing.
//             }
//         }

        return isAttachedOrNot;
    }
    else if ( isEllipticFaceAttached == rg_TRUE ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_BOOL BetaFace::isNonEllipticFaceAttached() const
{
    rg_BOOL isAttached = rg_FALSE;

    BetaVertex* vertex[3];
    searchVerticesInIntraWorld( vertex );

    Plane centerPlane;
    centerPlane.definePlaneByThreePoints( vertex[0]->getBall().getCenter(), 
                                          vertex[1]->getBall().getCenter(), 
                                          vertex[2]->getBall().getCenter()  );

    if ( m_leftCell->isVirtual() ) {
        rg_Point3D endVoronoiVtx     = m_rightCell->getMinTangentSphere().getCenter();

		rg_REAL    dist = centerPlane.distanceFromPoint( endVoronoiVtx );
        if ( rg_NEG( dist ) ) {
            isAttached = rg_TRUE;
        }
    }
    else if ( m_rightCell->isVirtual() ) {
        rg_Point3D startVoronoiVtx   = m_leftCell->getMinTangentSphere().getCenter();

		rg_REAL    dist = centerPlane.distanceFromPoint( startVoronoiVtx );
        if ( rg_POS( dist ) ) {
            isAttached = rg_TRUE;
        }
    }
    else {
        rg_Point3D startVoronoiVtx   = m_leftCell->getMinTangentSphere().getCenter();
        rg_Point3D endVoronoiVtx     = m_rightCell->getMinTangentSphere().getCenter();

		rg_REAL    distToStart = centerPlane.distanceFromPoint( startVoronoiVtx );
		rg_REAL    distToEnd   = centerPlane.distanceFromPoint( endVoronoiVtx );
        
        if ( rg_NEG( distToStart ) && rg_POS( distToEnd ) ) {
            isAttached = rg_FALSE;
        }
        else {
            isAttached = rg_TRUE;
        }
    }

    return isAttached;
}




///////////////////////////////////////////////////////////////////////////
//  Quasi-operator
//  
//  1. Q(f, Y | W): Query primitives for intra-world.
void BetaFace::searchVerticesInIntraWorld( BetaVertex** vertexArray, BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        if ( m_orientation[i] )  {
            vertexArray[i] = m_edge[i]->getStartVertex();
        }
        else  {
            vertexArray[i] = m_edge[i]->getEndVertex();
        }
    }
}



void BetaFace::searchVerticesInIntraWorld( rg_dList<BetaVertex*>& vertexList, BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        if ( m_orientation[i] )  {
            vertexList.add( m_edge[i]->getEndVertex() );
        }
        else  {
            vertexList.add( m_edge[i]->getStartVertex() );
        }
    }
}



void BetaFace::searchEdgesInIntraWorld(    rg_dList<BetaEdge*>&   edgeList,   BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        edgeList.add( m_edge[i] );
    }
}



void BetaFace::searchFacesInIntraWorld(    rg_dList<BetaFace*>&   faceList,   BetaFace* world ) const
{
    BetaCell* incidentCell[2] = {m_leftCell, m_rightCell};
    for ( rg_INT i_cell=0; i_cell<2; i_cell++ ) {
	    BetaFace** faceOnCell = incidentCell[i_cell]->getFaces();
        for ( rg_INT i_face = 0; i_face<EIWDS_NUM_FACE_ON_CELL; i_face++ ) {
            if ( this != faceOnCell[i_face] ) {
                faceList.add( faceOnCell[i_face] );
            }
        }
    }
}



void BetaFace::searchCellsInIntraWorld(    rg_dList<BetaCell*>&   cellList,   BetaFace* world ) const
{
    cellList.add( m_leftCell );
    cellList.add( m_rightCell );
}



void BetaFace::searchFiniteVerticesInIntraWorld( rg_dList<BetaVertex*>& finiteVertexList, BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        BetaVertex* vertex = rg_NULL;
        if ( m_orientation[i] )  {
            vertex = m_edge[i]->getEndVertex();
        }
        else  {
            vertex = m_edge[i]->getStartVertex();
        }


        if ( !vertex->isVirtual() ) {
            finiteVertexList.add( vertex );
        }
    }
}



void BetaFace::searchFiniteEdgesInIntraWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        if ( !m_edge[i]->isVirtual() ) {
            finiteEdgeList.add( m_edge[i] );
        }
    }
}



void BetaFace::searchFiniteFacesInIntraWorld( rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world ) const
{
    BetaCell* incidentCell[2] = {m_leftCell, m_rightCell};
    for ( rg_INT i_cell=0; i_cell<2; i_cell++ ) {
        if ( incidentCell[i_cell]->isVirtual() )  {
            continue;
        }

	    BetaFace** faceOnCell = incidentCell[i_cell]->getFaces();
        for ( rg_INT i_face = 0; i_face<EIWDS_NUM_FACE_ON_CELL; i_face++ ) {
            if ( this != faceOnCell[i_face] ) {
                finiteFaceList.add( faceOnCell[i_face] );
            }
        }
    }
}


void BetaFace::searchFiniteCellsInIntraWorld( rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world ) const
{
    if ( !m_leftCell->isVirtual() ) {
        finiteCellList.add( m_leftCell );
    }

    if ( !m_rightCell->isVirtual() ) {
        finiteCellList.add( m_rightCell );
    }
}

//  2. Query primitives for inter-world.
// void BetaFace::inquireSmallWorldsAtFace( rg_dList<BetaFace*>& smallWorldList, BetaFace* world )
// {
// }




//  3. Q(f, Y): Query primitives for whole-world.
void BetaFace::searchVerticesInWholeWorld( rg_dList<BetaVertex*>& vertexList ) const
{
    searchVerticesInIntraWorld( vertexList );
}



void BetaFace::searchEdgesInWholeWorld(    rg_dList<BetaEdge*>&   edgeList   ) const
{
    searchEdgesInIntraWorld( edgeList );
}



void BetaFace::searchFacesInWholeWorld(    rg_dList<BetaFace*>&   faceList   ) const
{
    searchFacesInIntraWorld( faceList );
}



void BetaFace::searchCellsInWholeWorld(    rg_dList<BetaCell*>&   cellList   ) const
{
    searchCellsInIntraWorld( cellList );
}



void BetaFace::searchFiniteVerticesInWholeWorld( rg_dList<BetaVertex*>& finiteVertexList ) const
{
    searchFiniteVerticesInIntraWorld(finiteVertexList);
}



void BetaFace::searchFiniteEdgesInWholeWorld(    rg_dList<BetaEdge*>&   finiteEdgeList   ) const
{
    searchFiniteEdgesInIntraWorld(finiteEdgeList);
}



void BetaFace::searchFiniteFacesInWholeWorld(    rg_dList<BetaFace*>&   finiteFaceList   ) const
{
    searchFiniteFacesInIntraWorld(finiteFaceList);
}



void BetaFace::searchFiniteCellsInWholeWorld(    rg_dList<BetaCell*>&   finiteCellList   ) const
{
    searchFiniteCellsInIntraWorld(finiteCellList);
}




//
//  END of Quasi-operator
///////////////////////////////////////////////////////////////////////////};



///////////////////////////////////////////////////////////////////////////
//  Beta-operator
//  
void BetaFace::searchFacesInBetaComplex( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const
{
    rg_dList<BetaFace*> finiteFaceList;
    searchFiniteFacesInWholeWorld( finiteFaceList );

    BetaFace* currFace = rg_NULL;
    finiteFaceList.reset4Loop();
    while ( finiteFaceList.setNext4Loop() ) {
        currFace = finiteFaceList.getEntity();

        if ( currFace->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            continue;
        }
        else {
            betaFaceList.add( currFace );
        }
    }
}



void BetaFace::searchCellsInBetaComplex( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList ) const
{
    rg_dList<BetaCell*> finiteCellList;
    searchFiniteCellsInWholeWorld( finiteCellList );

    BetaCell* currCell = rg_NULL;
    finiteCellList.reset4Loop();
    while ( finiteCellList.setNext4Loop() ) {
        currCell = finiteCellList.getEntity();

        if ( currCell->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            continue;
        }
        else {
            betaCellList.add( currCell );
        }
    }
}




void BetaFace::searchFacesInBetaShape(   const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const
{
    rg_dList<BetaFace*> finiteFaceList;
    searchFiniteFacesInWholeWorld( finiteFaceList );

    BetaFace* currFace = rg_NULL;
    finiteFaceList.reset4Loop();
    while ( finiteFaceList.setNext4Loop() ) {
        currFace = finiteFaceList.getEntity();

        rg_INT state = currFace->getBoundingState(beta);
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaFaceList.add( currFace );
        }
    }
}



//
//  END of Quasi-operator
///////////////////////////////////////////////////////////////////////////};

void BetaFace::connectVEdge(VDEdge* v_edge)
{
    m_vEdge = v_edge;
}



void BetaFace::disconnectVEdge(VDEdge* v_edge)
{
    if ( m_vEdge == v_edge ) {
        m_vEdge = rg_NULL;
    }
}


