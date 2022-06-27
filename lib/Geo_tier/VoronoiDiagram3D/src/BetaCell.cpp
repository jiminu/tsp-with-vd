#include "BetaCell.h"

#include "BetaVertex.h"
#include "BetaEdge.h"
#include "BetaFace.h"
#include "VDVertex.h"
using namespace V::GeometryTier;



BetaCell::BetaCell()
: m_visited(rg_FALSE)
{
    rg_INT i=0;
    for ( i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
        m_face[i] = 0;
    }

    for ( i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
        m_vertex[i] = 0;
    }

    m_vVertex = rg_NULL;
}



BetaCell::BetaCell(const rg_INT& ID)
: TopologicalEntity(ID),
  m_visited(rg_FALSE)
{
    rg_INT i=0;
    for ( i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
        m_face[i] = 0;
    }

    for ( i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
        m_vertex[i] = 0;
    }

    m_vVertex = rg_NULL;
}



BetaCell::BetaCell(const BetaCell& betaCell)
: TopologicalEntity(betaCell.m_ID),
  m_betaSpan(betaCell.m_betaSpan),
  m_minTangentSphere(betaCell.m_minTangentSphere), m_visited(betaCell.m_visited)
{
    rg_INT i=0;
    for ( i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
        m_face[i]   = betaCell.m_face[i];
    }

    for ( i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
        m_vertex[i] = betaCell.m_vertex[i];
    }

    m_vVertex = betaCell.m_vVertex;
}



BetaCell::~BetaCell()
{
    rg_INT i=0;
    for ( i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
        m_face[i] = rg_NULL;
    }

    for ( i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
        m_vertex[i] = rg_NULL;
    }

    //if ( m_vVertex != rg_NULL ) {
    //    m_vVertex->disconnectBetaCell(this);
    //}
}




rg_BOOL BetaCell::isVirtual() const
{
    rg_BOOL isVirtualCell = rg_FALSE;
    for ( rg_INT i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ )  {
        if ( m_vertex[i]->isVirtual() ) {
            isVirtualCell = rg_TRUE;
            break;
        }
    }

    return isVirtualCell;
}



BetaFace** BetaCell::getFaces() 
{
    return m_face;
}



BetaFace*    BetaCell::getFace(const rg_INT& i) const
{
    if ( i>=0 && i<EIWDS_NUM_FACE_ON_CELL ) {
        return m_face[i];
    }
    else {
        return rg_NULL;
    }
}



rg_INT BetaCell::getFacePos(const BetaFace* const face) const
{
    rg_INT facePos = -1;
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
        if ( face == m_face[i] ) {
            facePos = i;
            break;
        }
    }

    return facePos;
}



BetaCell* BetaCell::getNeighborCell(const rg_INT& i) const
{
    BetaCell* neighbor = rg_NULL;

    if ( i>=0 && i<EIWDS_NUM_FACE_ON_CELL ) {
        if ( this == m_face[i]->getLeftCell() ) {
            neighbor = m_face[i]->getRightCell(); 
        }
        else {
            neighbor = m_face[i]->getLeftCell(); 
        }
    }

    return neighbor;
}



BetaVertex** BetaCell::getVertices() 
{
    return m_vertex;
}



BetaVertex*  BetaCell::getVertex(const rg_INT& i) const
{
    if ( i>=0 && i<EIWDS_NUM_VERTEX_ON_CELL ) {
        return m_vertex[i];
    }
    else {
        return rg_NULL;
    }
}


    
rg_INT BetaCell::getVertexPos(const BetaVertex* const vertex) const
{
    rg_INT vertexPos = -1;
    for ( rg_INT i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
        if ( vertex == m_vertex[i] ) {
            vertexPos = i;
            break;
        }
    }

    return vertexPos;
}



rg_BetaSpan BetaCell::getBetaSpan() const
{
    return m_betaSpan;
}



rg_REAL BetaCell::getMinValueForValidSimplex() const
{
    return m_betaSpan.getLowerValueOfInterval(1);
}



//Sphere BetaCell::getMinTangentSphere() const
//{
//    return m_minTangentSphere;
//}

rg_INT BetaCell::getDepth() const
{
    rg_INT depth = rg_UNKNOWN;

    for ( rg_INT i=0; i<4; i++ ) {
        if ( depth < m_vertex[i]->getDepth() ) {
            depth = m_vertex[i]->getDepth();
        }
    }

    return depth;
}



void BetaCell::setFaces(BetaFace** faces)
{
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
        m_face[i] = faces[i];
    }
}



void BetaCell::setFace(const rg_INT& i, BetaFace* face)
{
    if ( i>=0 && i<EIWDS_NUM_FACE_ON_CELL ) {
        m_face[i] = face;
    }
}



void BetaCell::setVertices(BetaVertex** vertices)
{
    for ( rg_INT i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
        m_vertex[i] = vertices[i];
    }
}



void BetaCell::setVertex(const rg_INT& i, BetaVertex* vertex)
{
    if ( i>=0 && i<EIWDS_NUM_VERTEX_ON_CELL ) {
        m_vertex[i] = vertex;
    }
}



void BetaCell::setFaceAndItsMateVertex(const rg_INT& i, BetaFace* face, BetaVertex* vertex)
{
    if ( i>=0 && i<EIWDS_NUM_VERTEX_ON_CELL ) {
        m_vertex[i] = vertex;
        m_face[i]   = face;
    }
}



void BetaCell::addFaceAndItsMateVertex(BetaFace* face, BetaVertex* vertex)
{
    for ( rg_INT i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
        if ( m_vertex[i] == rg_NULL ) {
            m_vertex[i] = vertex;
            m_face[i]   = face;
        }
    }
}



void BetaCell::setBetaSpan(const rg_BetaSpan& betaSpan)
{
    m_betaSpan = betaSpan;
}



void BetaCell::setMinTangentSphere(const Sphere& sphere)
{
    m_minTangentSphere = sphere;
}

void BetaCell::shiftBetaSpan(const rg_REAL& delta)
{
	m_betaSpan.shiftBetaSpan(delta);
}



BetaCell& BetaCell::operator=(const BetaCell& betaCell)
{
    if ( this == &betaCell ) {
        return *this;
    }

    m_ID = betaCell.m_ID;

    rg_INT i=0;
    for ( i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
        m_face[i]   = betaCell.m_face[i];
    }
    for ( i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
        m_vertex[i] = betaCell.m_vertex[i];
    }
    m_betaSpan = betaCell.m_betaSpan;

    m_minTangentSphere = betaCell.m_minTangentSphere;
    m_visited          = betaCell.m_visited;

    m_vVertex = betaCell.m_vVertex;

    return *this;
}



rg_BOOL BetaCell::isThere(BetaVertex* vertex) const
{
    for (rg_INT i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ )  {
        if ( vertex == m_vertex[i] ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_BOOL BetaCell::isThere(BetaEdge* edge) const
{
    for (rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL-1; i++ )  {
        if ( m_face[i]->isThere( edge ) ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_BOOL BetaCell::isThere(BetaFace* face) const
{
    for (rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        if ( m_face[i] == face ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_BOOL BetaCell::isAdjacent(const BetaCell* const  cell) const
{
    for (rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        if ( this == m_face[i]->getLeftCell() ) {
            if ( cell == m_face[i]->getRightCell() ) {
                return rg_TRUE;
            }
        }
        else {
            if ( cell == m_face[i]->getLeftCell() ) {
                return rg_TRUE;
            }
        }
    }

    return rg_FALSE;
}



rg_BOOL BetaCell::isIntersectedAtEdgeWith(BetaFace* face) const
{
    BetaEdge** edgesOnFace = face->getEdges();
    for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++) {
        if ( isThere(edgesOnFace[i]) ) {
            return rg_TRUE;
        }
    }
    return rg_FALSE;
}



rg_BOOL BetaCell::isIntersectedAtEdgeWith(BetaCell* cell) const
{
    rg_dList<BetaEdge*> edgeList;
    searchEdgesInIntraWorld( edgeList );

    BetaEdge* currEdge = rg_NULL;
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        currEdge = edgeList.getEntity();

        if ( isThere( currEdge ) ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}




// rg_BOOL BetaCell::isOnConvexHull() const
// {
//     return rg_FALSE;
// }



rg_INT BetaCell::isMultiplicity() const
{
    BetaCell* neighbor[EIWDS_NUM_FACE_ON_CELL];
    getNeighborCells(neighbor);

    rg_INT multiplicity = 1;
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
        for ( rg_INT j=i+1; j<EIWDS_NUM_FACE_ON_CELL; j++ ) {
            if ( neighbor[i] == neighbor[j] ) {
                multiplicity++;
            }
        }

        if ( multiplicity != 1 ) {
            break;
        }
    }

    return multiplicity;
}



BetaCell* BetaCell::getNeighborCellWithMultiplicity() const
{
    BetaCell* neighbor[EIWDS_NUM_FACE_ON_CELL];
    getNeighborCells(neighbor);

    BetaCell* multipleNeighbor = rg_NULL;
    rg_INT multiplicity = 0;
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
        for ( rg_INT j=i+1; j<EIWDS_NUM_FACE_ON_CELL; j++ ) {
            if ( neighbor[i] == neighbor[j] ) {
                multiplicity++;
                multipleNeighbor = neighbor[i];
            }
        }

        if ( multiplicity != 0 ) {
            multiplicity++;
            break;
        }
    }

    return multipleNeighbor;
}




BetaFace* BetaCell::findSharingFaceWith(const BetaCell* const neighborCell) const
{
    for (rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        if ( this == m_face[i]->getLeftCell() ) {
            if ( neighborCell == m_face[i]->getRightCell() ) {
                return m_face[i];
            }
        }
        else {
            if ( neighborCell == m_face[i]->getLeftCell() ) {
                return m_face[i];
            }
        }
    }

    return rg_NULL;
}



BetaEdge* BetaCell::findEdge(const BetaVertex* const vertex1, const BetaVertex* const vertex2) const 
{
    BetaEdge* edgeToBeFound = rg_NULL;
    for (rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL-1; i++ )  {
        edgeToBeFound = m_face[i]->findEdge( vertex1, vertex2 );

        if ( edgeToBeFound != rg_NULL ) {
            break;
        }
    }

    return edgeToBeFound;
}


void BetaCell::findEdgesIncidentToVertex(const BetaVertex* const vertex, BetaEdge** edge) const
{
    rg_INT i_edge = 0;
    for ( rg_INT i_vtx=0; i_vtx<EIWDS_NUM_VERTEX_ON_CELL; i_vtx++ ) {
        if ( vertex == m_vertex[i_vtx] ) {
            continue;
        }

	    BetaEdge** edgesOnFace = m_face[i_vtx]->getEdges();
        for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ )  {
            if ( edgesOnFace[i]->isVisited() ) {
                continue;
            }

            if ( vertex == edgesOnFace[i]->getStartVertex() || vertex == edgesOnFace[i]->getEndVertex() ) {
                edgesOnFace[i]->isVisited(rg_TRUE);
                edge[i_edge] = edgesOnFace[i];
                i_edge++;
            }
        }        
    }

    for ( i_edge=0; i_edge<EIWDS_NUM_EDGE_ON_FACE; i_edge++ )  {
        edge[i_edge]->isVisited( rg_FALSE );
    }
}


void BetaCell::findFacesIncidentToVertex(const BetaVertex* const vertex, BetaFace** face) const
{
    rg_INT i_face = 0;
    for ( rg_INT i_vtx=0; i_vtx<EIWDS_NUM_VERTEX_ON_CELL; i_vtx++ ) {
        if ( vertex != m_vertex[i_vtx] ) {
            face[i_face] = m_face[i_vtx];
            i_face++;
        }
    }
}



void BetaCell::findFacesIncidentToEdge(const BetaEdge* const edge, BetaFace** face) const
{
    // thisCell이 givenEdge로 bounding되어있지 않으면 null값 return.
    if ( isThere((BetaEdge*)edge) == rg_FALSE ) {
        face[0] = rg_NULL;
        face[1] = rg_NULL;
        return;
    }

    BetaVertex* vertexOnEdge[2] = { edge->getStartVertex(), edge->getEndVertex() };
    rg_BOOL          isOnEdge[4]     = { rg_FALSE, rg_FALSE, rg_FALSE, rg_FALSE };

    for ( rg_INT i=0; i<2; i++ )  {
        for ( rg_INT i_vtx=0; i_vtx<EIWDS_NUM_VERTEX_ON_CELL; i_vtx++ ) {
            if ( vertexOnEdge[i] == m_vertex[i_vtx] ) {
                isOnEdge[i_vtx] = rg_TRUE;
                break;
            }
        }
    }

    rg_INT i_face = 0;
    for ( rg_INT i_vtx=0; i_vtx<EIWDS_NUM_VERTEX_ON_CELL; i_vtx++ ) {
        if ( !isOnEdge[i_vtx] ) {
            face[i_face++] = m_face[i_vtx];
        }
    }



    // face[]가 givenEdge기준으로 CCW 순서로 나오도록 한다.
    // face[0] 에서 givenEdge의 방향은
    //              thisCell이 face[0]의 LeftCell일 경우:  (-) 방향이어야 하고,
    //              thisCell이 face[0]의 RightCell일 경우: (+) 방향이어야 한다.
    // 위 조건을 만족하지 않을 경우 face[0]과 face[1]을 swapping해준다.

    if ( this == face[0]->getLeftCell() ) {
        if ( face[0]->getEdgeOrientation( edge ) == rg_TRUE ) {
            BetaFace* tempFace  = face[0];
            face[0]             = face[1];
            face[1]             = tempFace;
        }
    }
    else { // if ( this == face[0]->getRightCell() ) {
        if ( face[0]->getEdgeOrientation( edge ) == rg_FALSE ) {
            BetaFace* tempFace  = face[0];
            face[0]             = face[1];
            face[1]             = tempFace;
        }
    }
}



rg_REAL BetaCell::computeAngleBtwFacesIncidnetToEdge( const BetaEdge*   const edge ) const
{
    rg_REAL angle = -rg_REAL_INFINITY;
    if ( isThere((BetaEdge*)edge) == rg_FALSE ) {
        angle = -rg_REAL_INFINITY;
    }
    else {
        rg_Point3D  startPoint = edge->getStartVertex()->getBall().getCenter();
        rg_Point3D  endPoint   = edge->getEndVertex()->getBall().getCenter();
        rg_Point3D  normal     = endPoint - startPoint;
        normal.normalize();
        Plane       basePlane(normal, startPoint);


        BetaFace* incidentFace[2] = {rg_NULL, rg_NULL};
        findFacesIncidentToEdge(edge, incidentFace);

        BetaVertex* mateVertex[2]          = { incidentFace[0]->findMateVertexOfEdge(edge), incidentFace[1]->findMateVertexOfEdge(edge) };
        rg_Point3D  matePoint[2]           = { mateVertex[0]->getBall().getCenter(),        mateVertex[1]->getBall().getCenter() };
        rg_Point3D  projectedMatePoint[2]  = { basePlane.projectPointOnPlane(matePoint[0]), basePlane.projectPointOnPlane(matePoint[1]) };

        angle = basePlane.computeAngleInCCW(projectedMatePoint[0], startPoint, projectedMatePoint[1]);

        rg_REAL volume = computeSignedVolume();
        if ( rg_LT( volume, 0.0 ) ) {
            rg_REAL betaValueIntoInterior = m_minTangentSphere.getRadius();

            rg_BOOL isCommonFace[2] = {rg_FALSE, rg_FALSE};
            for ( rg_INT j=0; j<2; j++ ) {
                Interval_REAL regularInterval;
                incidentFace[j]->getpBetaSpan()->findBetaIntervalOfGivenState(REGULAR_SIMPLEX, regularInterval);

                if ( betaValueIntoInterior == regularInterval.getLowerValue() ) {
                    isCommonFace[j] = rg_FALSE;
                }
                else if ( betaValueIntoInterior == regularInterval.getUpperValue() ) {
                    isCommonFace[j] = rg_TRUE;
                }
            }

            if ( isCommonFace[0] != isCommonFace[1] ) {
                angle = -(ANGLE_360_DEGREE - angle);
            }
        }
    }

    return angle;
}



void BetaCell::getNeighborCells(BetaCell** neighborCells) const
{
    for (rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        if ( this == m_face[i]->getLeftCell() ) {
            neighborCells[i] = m_face[i]->getRightCell(); 
        }
        else {
            neighborCells[i] = m_face[i]->getLeftCell();
        }
    }
}



BetaCell* BetaCell::getNeighborCell(const BetaFace* const face) const
{
    for (rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        if ( face == m_face[i] ) {
            if ( this == m_face[i]->getLeftCell() ) {
                return m_face[i]->getRightCell(); 
            }
            else {
                return m_face[i]->getLeftCell(); 
            }
        }
    }

    return rg_NULL;
}


BetaVertex* BetaCell::getMateVertex(const BetaFace* const face) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        if ( face == m_face[i] ) {
            return m_vertex[i];
        }
    }

    return rg_NULL;
}



BetaFace* BetaCell::getMateFace(const BetaVertex* const vertex) const
{
    BetaFace* mateFace = rg_NULL;
    for ( rg_INT i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ )  {
        if ( vertex == m_vertex[i] ) {
            mateFace =  m_face[i]; 
            break;
        }
    }

    return mateFace;
}



BetaCell* BetaCell::getMateCell(const rg_INT& i_vtx) const
{
    if ( i_vtx>=0 && i_vtx<EIWDS_NUM_VERTEX_ON_CELL ) {
        if ( this == m_face[i_vtx]->getLeftCell() ) {
            return m_face[i_vtx]->getRightCell(); 
        }
        else {
            return m_face[i_vtx]->getLeftCell(); 
        }
    }

    return rg_NULL;
}



BetaCell* BetaCell::getMateCell(const BetaVertex* const vertex) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ )  {
        if ( vertex == m_vertex[i] ) {
            if ( this == m_face[i]->getLeftCell() ) {
                return m_face[i]->getRightCell(); 
            }
            else {
                return m_face[i]->getLeftCell(); 
            }
        }
    }

    return rg_NULL;
}

void BetaCell::updateRadiusOfMinimumTangentSphere(const rg_REAL& deltaOfBallRadii)
{
	rg_REAL originalRadius = m_minTangentSphere.getRadius();
	m_minTangentSphere.setRadius(originalRadius + deltaOfBallRadii);
}

void BetaCell::computeBetaSpan()
{
    rg_REAL minOfCurrSimplex = m_minTangentSphere.getRadius();

    m_betaSpan.setNumBetaInterVal(2);
    m_betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
    m_betaSpan.setBetaInterval(1, INTERIOR_SIMPLEX, minOfCurrSimplex, REAL_PLUS_INFINITE);
}

rg_REAL BetaCell::computeSignedVolume_old() const
{
    rg_Point3D pt[EIWDS_NUM_VERTEX_ON_CELL] = { m_vertex[0]->getBall().getCenter(),
                                                m_vertex[1]->getBall().getCenter(),
                                                m_vertex[2]->getBall().getCenter(),
                                                m_vertex[3]->getBall().getCenter() };

	rg_REAL lenSide[3]     = { pt[1].distance(pt[2]), pt[2].distance(pt[3]), pt[3].distance(pt[1]) };	
	rg_REAL halfPerimeter  = (lenSide[0] + lenSide[1] + lenSide[2])*0.5;
	rg_REAL squaredArea    = halfPerimeter*(halfPerimeter-lenSide[0])*(halfPerimeter-lenSide[1])*(halfPerimeter-lenSide[2]);
	rg_REAL areaOfBaseTriangle = sqrt( squaredArea );

	Plane   basePlane;
	basePlane.definePlaneByThreePoints(pt[1], pt[2], pt[3]);

	rg_REAL heightOfTetrahedron = basePlane.distanceFromPoint( pt[0] );

	rg_REAL signedVolume = (areaOfBaseTriangle*heightOfTetrahedron)/3.0;

    return signedVolume;
}

rg_REAL BetaCell::computeSignedVolume() const
{
	// comment about orientation...
//     rg_Point3D pt[EIWDS_NUM_VERTEX_ON_CELL] = { m_vertex[0]->getBall().getCenter(),
//                                                 m_vertex[1]->getBall().getCenter(),
//                                                 m_vertex[2]->getBall().getCenter(),
//                                                 m_vertex[3]->getBall().getCenter() };
// 
// 
// 	rg_REAL x0 = pt[1].getX();
// 	rg_REAL y0 = pt[1].getY();
// 	rg_REAL z0 = pt[1].getZ();
// 
// 	rg_REAL x1 = pt[2].getX();
// 	rg_REAL y1 = pt[2].getY();
// 	rg_REAL z1 = pt[2].getZ();
// 
// 	rg_REAL x2 = pt[3].getX();
// 	rg_REAL y2 = pt[3].getY();
// 	rg_REAL z2 = pt[3].getZ();
// 
// 	rg_REAL x3 = pt[0].getX();
// 	rg_REAL y3 = pt[0].getY();
// 	rg_REAL z3 = pt[0].getZ();
// 
// 	// Optimized by Maple 12
// 	// Before optimization: 23 additions 48 multiplications
// 	// Afeter optimization: 17 additions 16 multiplications 7 assignments
// 
// 	rg_REAL t1, t2, t3, t4, t5, t6;
// 	t6 =  z0-z2;
// 	t5 =  z1-z0; 
// 	t4 =  z1-z3; 
// 	t3 =  z2-z1; 
// 	t2 =  z2-z3; 
// 	t1 = -z3+z0; 
// 	rg_REAL signedVolume 
// 		= ((-t5*y2-t6*y1-t3*y0)*x3+(t5*y3+t1*y1-t4*y0)*x2+(t6*y3-t1*y2+t2*y0)*x1+(t3*y3+t4*y2-t2*y1)*x0)/6.0;
// 	
// 
//     return signedVolume;

	// Positive orientation in QT:
	// Seeing from m_vertex[0], m_vertex[1], m_vertex[2] and m_vertex[3] are ordered in CCW.
	

	// Positive orientation in formula for computing the signed volume of a tetrahedron:
	// Seeing from a vertex, the three remaining vertices should be ordered in CW.
	// Hence, we reorder the vertices as follows
	// so that the formula could produce the positive volume for a QT cell in positive orientation.

	// i.e.
	// Positive orientation in QT     : v[0]->(v[1], v[2], v[3]) in CCW
	// Positive orientation in Formula: v[3]->(v[0], v[1], v[2]) in  CW
	
	rg_Point3D pt[EIWDS_NUM_VERTEX_ON_CELL] = { m_vertex[3]->getBall().getCenter(),
                                                m_vertex[0]->getBall().getCenter(),
                                                m_vertex[1]->getBall().getCenter(),
                                                m_vertex[2]->getBall().getCenter() };


	rg_REAL x0 = pt[0].getX();
	rg_REAL y0 = pt[0].getY();
	rg_REAL z0 = pt[0].getZ();

	rg_REAL x1 = pt[1].getX();
	rg_REAL y1 = pt[1].getY();
	rg_REAL z1 = pt[1].getZ();

	rg_REAL x2 = pt[2].getX();
	rg_REAL y2 = pt[2].getY();
	rg_REAL z2 = pt[2].getZ();

	rg_REAL x3 = pt[3].getX();
	rg_REAL y3 = pt[3].getY();
	rg_REAL z3 = pt[3].getZ();

	// Optimized by Maple 12
	// Before optimization: 23 additions 48 multiplications
	// Afeter optimization: 17 additions 16 multiplications 7 assignments

	rg_REAL t1, t2, t3, t4, t5, t6;
	t6 =  z0-z2;
	t5 =  z1-z0; 
	t4 =  z1-z3; 
	t3 =  z2-z1; 
	t2 =  z2-z3; 
	t1 = -z3+z0; 
	rg_REAL signedVolume 
		= ((-t5*y2-t6*y1-t3*y0)*x3+(t5*y3+t1*y1-t4*y0)*x2+(t6*y3-t1*y2+t2*y0)*x1+(t3*y3+t4*y2-t2*y1)*x0)/6.0;
	

    return signedVolume;
}



///////////////////////////////////////////////////////////////////////////
//  Quasi-operator
//  
//  1. Q(c, Y | W): Query primitives for intra-world.
void BetaCell::searchVerticesInIntraWorld(rg_dList<BetaVertex*>& vertexList, BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ )  {
        vertexList.add( m_vertex[i] );
    }
}



void BetaCell::searchEdgesInIntraWorld(   rg_dList<BetaEdge*>&   edgeList,   BetaFace* world ) const
{
    for ( rg_INT i_face=0; i_face<(EIWDS_NUM_FACE_ON_CELL-1); i_face++ )  {
        BetaEdge** edgeOnFace = m_face[i_face]->getEdges();

        for ( rg_INT i_edge=0; i_edge<EIWDS_NUM_EDGE_ON_FACE; i_edge++ ) {
            if ( edgeOnFace[i_edge]->isVisited() == rg_FALSE ) {
                edgeOnFace[i_edge]->isVisited(rg_TRUE);
                edgeList.add( edgeOnFace[i_edge] );
            }
        }
    }

    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        edgeList.getEntity()->isVisited( rg_FALSE );
    }
}



void BetaCell::searchFacesInIntraWorld(   rg_dList<BetaFace*>&   faceList,   BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        faceList.add( m_face[i] );
    }
}



void BetaCell::searchCellsInIntraWorld(   rg_dList<BetaCell*>&   cellList,   BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_VERTEX_ON_CELL; i++ )  {
        if ( this == m_face[i]->getLeftCell() ) {
            cellList.add( m_face[i]->getRightCell() );
        }
        else {
            cellList.add( m_face[i]->getLeftCell() );
        }
    }
}



void BetaCell::searchFiniteVerticesInIntraWorld( rg_dList<BetaVertex*>& finiteVertexList, BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        if ( !m_vertex[i]->isVirtual() ) {
            finiteVertexList.add( m_vertex[i] );
        }
    }
}



void BetaCell::searchFiniteEdgesInIntraWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world ) const
{
    rg_dList<BetaEdge*> edgeList;
    searchEdgesInIntraWorld( edgeList );

    BetaEdge* currEdge = rg_NULL;
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        currEdge = edgeList.getEntity();

        if ( !currEdge->isVirtual() ) {
            finiteEdgeList.add( currEdge );
        }
    }
}



void BetaCell::searchFiniteFacesInIntraWorld(    rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world ) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        if ( !m_face[i]->isVirtual() ) {
            finiteFaceList.add( m_face[i] );
        }
    }
}



void BetaCell::searchFiniteCellsInIntraWorld( rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world) const
{
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
        BetaCell* neighbor = rg_NULL;
        if ( this == m_face[i]->getLeftCell() ) {
            neighbor = m_face[i]->getRightCell();
        }
        else {
            neighbor = m_face[i]->getLeftCell();
        }

        if ( !neighbor->isVirtual() ) {
            finiteCellList.add( neighbor );
        }
    }
}

//  2. Query primitives for inter-world.
// void BetaCell::inquireSmallWorldsAtCell(  rg_dList<BetaFace*>& smallWorldList, BetaFace* world )
// {
// 
// }




//  3. Q(c, Y): Query primitives for whole-world.
void BetaCell::searchVerticesInWholeWorld(rg_dList<BetaVertex*>& vertexList ) const
{
    searchVerticesInIntraWorld(vertexList);
}



void BetaCell::searchEdgesInWholeWorld(   rg_dList<BetaEdge*>&   edgeList   ) const
{
    searchEdgesInIntraWorld(edgeList);
}



void BetaCell::searchFacesInWholeWorld(   rg_dList<BetaFace*>&   faceList   ) const
{
    searchFacesInIntraWorld( faceList );
}



void BetaCell::searchCellsInWholeWorld(   rg_dList<BetaCell*>&   cellList   ) const
{
    searchCellsInIntraWorld( cellList );
}



void BetaCell::searchFiniteVerticesInWholeWorld( rg_dList<BetaVertex*>& finiteVertexList ) const
{
    searchFiniteVerticesInIntraWorld( finiteVertexList );
}



void BetaCell::searchFiniteEdgesInWholeWorld(    rg_dList<BetaEdge*>&   finiteEdgeList   ) const
{
    searchFiniteEdgesInIntraWorld( finiteEdgeList );
}



void BetaCell::searchFiniteFacesInWholeWorld(    rg_dList<BetaFace*>&   finiteFaceList   ) const
{
    searchFiniteFacesInIntraWorld( finiteFaceList );
}




void BetaCell::searchFiniteCellsInWholeWorld( rg_dList<BetaCell*>&   finiteCellList   ) const
{
    searchFiniteCellsInIntraWorld( finiteCellList );
}


//
//  END of Quasi-operator
///////////////////////////////////////////////////////////////////////////};




///////////////////////////////////////////////////////////////////////////
//  Beta-operator
//  
void BetaCell::searchCellsInBetaComplex( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList ) const
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



void BetaCell::searchEdgesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaEdge*>& betaEdgeList ) const
{
    rg_dList<BetaEdge*> finiteEdgeList;
    searchFiniteEdgesInWholeWorld( finiteEdgeList );

    BetaEdge* currEdge = rg_NULL;
    finiteEdgeList.reset4Loop();
    while ( finiteEdgeList.setNext4Loop() ) {
        currEdge = finiteEdgeList.getEntity();

        if ( currEdge->getBoundingState(beta) == boundingState ) {
            betaEdgeList.add( currEdge );
        }
    }
}



void BetaCell::searchFacesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaFace*>& betaFaceList ) const
{
    rg_dList<BetaFace*> finiteFaceList;
    searchFiniteFacesInWholeWorld( finiteFaceList );

    BetaFace* currFace = rg_NULL;
    finiteFaceList.reset4Loop();
    while ( finiteFaceList.setNext4Loop() ) {
        currFace = finiteFaceList.getEntity();

        if ( currFace->getBoundingState(beta) == boundingState ) {
            betaFaceList.add( currFace );
        }
    }
}



void BetaCell::searchCellsInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaCell*>& betaCellList ) const
{
    rg_dList<BetaCell*> finiteCellList;
    searchFiniteCellsInWholeWorld( finiteCellList );

    BetaCell* currCell = rg_NULL;
    finiteCellList.reset4Loop();
    while ( finiteCellList.setNext4Loop() ) {
        currCell = finiteCellList.getEntity();

        if ( currCell->getBoundingState(beta) == boundingState ) {
            betaCellList.add( currCell );
        }
    }
}




//
//  END of Quasi-operator
///////////////////////////////////////////////////////////////////////////};


// function to connect with vertex in VD
void BetaCell::connectVVertex(VDVertex* v_vtx)
{
    m_vVertex = v_vtx;
}

void BetaCell::disconnectVVertex(VDVertex* v_vtx)
{
    if ( m_vVertex == v_vtx ) {
        m_vVertex = rg_NULL;
    }
}
