#include "VDEdge.h"

#include "VDVertex.h"
#include "VDPartialEdge.h"
#include "VDLoop.h"
#include "VDFace.h"
#include "VDCell.h"
#include "rg_BallGenerator.h"
#include "rg_RBzCurve3D.h"
#include "Gate.h"
#include "BetaFace.h"
using namespace V::GeometryTier;


#include <float.h>

#include <set>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
VDEdge::VDEdge()
: m_startVertex(       rg_NULL ),
  m_endVertex(         rg_NULL ),
  m_partEdge(          rg_NULL ),
  m_edgeEquation(      rg_NULL ),
  m_numOfPointsOnEdge( 0 ),
  m_pointsOnEdge(      rg_NULL ),
  m_minRadius( 0 ),
  m_isOnInfinity( rg_FALSE ),
  m_ProbeTangibility( rg_UNKNOWN),
  m_edgeType( rg_UNKNOWN )
{
    m_visited  = rg_FALSE;
    m_betaFace = rg_NULL;
}

VDEdge::VDEdge( const rg_INT&        ID, 
                      VDVertex*      startVertex, 
                      VDVertex*      endVertex )
: TopologicalEntity(   ID ),
  m_startVertex(       startVertex ),
  m_endVertex(         endVertex ),
  m_partEdge(          rg_NULL ),
  m_edgeEquation(      rg_NULL ),
  m_numOfPointsOnEdge( 0 ),
  m_pointsOnEdge(      rg_NULL ),
  m_minRadius( 0 ),
  m_isOnInfinity( rg_FALSE ),
  m_ProbeTangibility( rg_UNKNOWN),
  m_edgeType( rg_UNKNOWN )
{
    m_visited  = rg_FALSE;
    m_betaFace = rg_NULL;
}


VDEdge::VDEdge(       VDVertex*      startVertex, 
                      VDVertex*      endVertex, 
                      VDPartialEdge* partialEdge, 
                      rg_Curve*      edgeEquation, 
                const rg_INT&        numPtsOnEdge, 
                      rg_Point3D*    ptsOnEdge )
: m_startVertex(       startVertex ),
  m_endVertex(         endVertex ),
  m_partEdge(          partialEdge ),
  m_edgeEquation(      edgeEquation ),
  m_numOfPointsOnEdge( numPtsOnEdge ),
  m_minRadius( 0 ),
  m_isOnInfinity( rg_FALSE ),
  m_ProbeTangibility( rg_UNKNOWN),
  m_edgeType( rg_UNKNOWN )
{
    if ( m_numOfPointsOnEdge == 0 )
    {
        m_pointsOnEdge = rg_NULL;
    }
    else 
    {
        //  Modified by Youngsong Cho 2005. 5. 24
        //  m_pointsOnEdge = new rg_Point3D[ m_numOfPointsOnEdge ];
        //  for (rg_INT i=0; i<m_numOfPointsOnEdge; i++)
        //    m_pointsOnEdge[i] = ptsOnEdge[i];
        m_pointsOnEdge = ptsOnEdge;
    }

    m_visited  = rg_FALSE;
    m_betaFace = rg_NULL;
}


VDEdge::VDEdge( const rg_INT&        ID,      
                      VDVertex*      startVertex, 
                      VDVertex*      endVertex, 
                      VDPartialEdge* partialEdge, 
                      rg_Curve*      edgeEquation, 
                const rg_INT&        numPtsOnEdge, 
                      rg_Point3D*    ptsOnEdge )
: TopologicalEntity(   ID ),
  m_startVertex(       startVertex ),
  m_endVertex(         endVertex ),
  m_partEdge(          partialEdge ),
  m_edgeEquation(      edgeEquation ),
  m_numOfPointsOnEdge( numPtsOnEdge ),
  m_minRadius( 0 ),
  m_isOnInfinity( rg_FALSE ),
  m_ProbeTangibility( rg_UNKNOWN),
  m_edgeType( rg_UNKNOWN )
{
    if ( m_numOfPointsOnEdge == 0 )
    {
        m_pointsOnEdge = rg_NULL;
    }
    else 
    {
        //  Modified by Youngsong Cho 2005. 5. 24
        //  m_pointsOnEdge = new rg_Point3D[ m_numOfPointsOnEdge ];
        //  for (rg_INT i=0; i<m_numOfPointsOnEdge; i++)
        //      m_pointsOnEdge[i] = ptsOnEdge[i];
        m_pointsOnEdge = ptsOnEdge;
    }

    m_visited  = rg_FALSE;
    m_betaFace = rg_NULL;
}


VDEdge::VDEdge( const TopologicalEntity& aTopoEntity,    
                      VDVertex*          startVertex, 
                      VDVertex*          endVertex, 
                      VDPartialEdge*     partialEdge, 
                      rg_Curve*          edgeEquation, 
                const rg_INT&            numPtsOnEdge, 
                      rg_Point3D*        ptsOnEdge )
: TopologicalEntity(   aTopoEntity ),
  m_startVertex(       startVertex ),
  m_endVertex(         endVertex ),
  m_partEdge(          partialEdge ),
  m_edgeEquation(      edgeEquation ),
  m_numOfPointsOnEdge( numPtsOnEdge ),
  m_minRadius( 0 ),
  m_isOnInfinity( rg_FALSE ),
  m_ProbeTangibility( rg_UNKNOWN),
  m_edgeType( rg_UNKNOWN )
{
    if ( m_numOfPointsOnEdge == 0 )
    {
        m_pointsOnEdge = rg_NULL;
    }
    else 
    {
        //  Modified by Youngsong Cho 2005. 5. 24
        //  m_pointsOnEdge = new rg_Point3D[ m_numOfPointsOnEdge ];
        //  for (rg_INT i=0; i<m_numOfPointsOnEdge; i++)
        //      m_pointsOnEdge[i] = ptsOnEdge[i];
        m_pointsOnEdge = ptsOnEdge;
    }

    m_visited  = rg_FALSE;
    m_betaFace = rg_NULL;
}


VDEdge::VDEdge( const VDEdge& anEdge )
: TopologicalEntity(   anEdge ),
  m_startVertex(       anEdge.m_startVertex ),
  m_endVertex(         anEdge.m_endVertex ),
  m_partEdge(          anEdge.m_partEdge ),
  m_numOfPointsOnEdge( anEdge.m_numOfPointsOnEdge ),
  m_minRadius(		   anEdge.m_minRadius ),
  m_isOnInfinity(      anEdge.m_isOnInfinity ),
  m_ProbeTangibility(  anEdge.m_ProbeTangibility),
  m_edgeType( anEdge.m_edgeType )
{
    m_minMaxTangentSphereByBalls[0] = anEdge.m_minMaxTangentSphereByBalls[0];
    m_minMaxTangentSphereByBalls[1] = anEdge.m_minMaxTangentSphereByBalls[1];
    m_edgePlane                     = anEdge.m_edgePlane;
    m_pointForAngleDist             = anEdge.m_pointForAngleDist;
    

    //  edgeEquation은 rg_Curve의 포인터이고, 
    //  다음은 객체의 복사가 아닌 단순히 포인터에 대한 복사만 이루어짐.
    //  rg_Curve에 대한 virtual function의 정의가 필요할 것으로 생각됨.
    //  임시로 Sphere set voronoi diagram 에 대한 edge equation 객체의 복사를 수행하게 함. 
    // 2004-4-6 by Youngsong Cho
    if ( anEdge.m_edgeEquation != rg_NULL )
    {
        if ( m_edgeEquation != rg_NULL )
            delete m_edgeEquation;

        m_edgeEquation = new rg_RBzCurve3D( *((rg_RBzCurve3D*) anEdge.m_edgeEquation) );
    }
    else
    {
        m_edgeEquation = rg_NULL;
    }

    

    if ( m_numOfPointsOnEdge == 0 )
    {
        m_pointsOnEdge = rg_NULL;
    }
    else 
    {
        m_pointsOnEdge = new rg_Point3D[ m_numOfPointsOnEdge ];
        for (rg_INT i=0; i<m_numOfPointsOnEdge; i++)
            m_pointsOnEdge[i] = anEdge.m_pointsOnEdge[i];
    }

    m_visited  = anEdge.m_visited;
    m_betaFace = anEdge.m_betaFace;
}


VDEdge::~VDEdge()
{
    if ( m_edgeEquation != rg_NULL ) {
        delete m_edgeEquation;
    }

    if ( m_pointsOnEdge != rg_NULL ) {
        delete [] m_pointsOnEdge;
    }

    //if ( m_betaFace != rg_NULL ) {
    //    m_betaFace->disconnectVEdge(this);
    //}
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
rg_BOOL VDEdge::isVirtual() const
{
    rg_dList<VDCell*> cellList;
    inquireIncidentCells( cellList );

    rg_BOOL isVirtualEdge = rg_FALSE;
    cellList.reset4Loop();
    while ( cellList.setNext4Loop() ) {
        VDCell* currCell = cellList.getEntity();

        if ( currCell->getGenerator() == rg_NULL ) {
            isVirtualEdge = rg_TRUE;
            break;
        }
    }

    return isVirtualEdge;
}


VDVertex* VDEdge::getStartVertex() const
{
    return m_startVertex;
}

VDVertex* VDEdge::getEndVertex() const
{
    return m_endVertex;
}

VDVertex* VDEdge::getOppositeVertex(const VDVertex* const vertex) const
{
    if ( vertex == m_startVertex ) {
        return m_endVertex;
    }
    else if ( vertex == m_endVertex ) {
        return m_startVertex;
    }
    else {
        return rg_NULL;
    }
}


VDCell* VDEdge::getStartVCell() const
{
    rg_dList<VDCell*> cellsDefiningThisEdge;
    inquireIntoOnlyEdgeSharingCellsInCCW( cellsDefiningThisEdge );


    set<VDCell*> setOfCellsDefiningThisEdge;
    cellsDefiningThisEdge.reset4Loop();
    while ( cellsDefiningThisEdge.setNext4Loop() ) {
        setOfCellsDefiningThisEdge.insert( cellsDefiningThisEdge.getEntity() );
    }


    rg_dList<VDCell*> cellsIncidnetToStartVtx;
    m_startVertex->inquireIncidentCells( cellsIncidnetToStartVtx );

    VDCell* startCell = rg_NULL;

    cellsIncidnetToStartVtx.reset4Loop();
    while ( cellsIncidnetToStartVtx.setNext4Loop() ) {
        VDCell* currCell = cellsIncidnetToStartVtx.getEntity();

        if ( setOfCellsDefiningThisEdge.find( currCell ) == setOfCellsDefiningThisEdge.end() ) {
            startCell = currCell;
            break;
        }
    }

    return startCell;
}



VDCell* VDEdge::getEndVCell() const
{
    rg_dList<VDCell*> cellsDefiningThisEdge;
    inquireIntoOnlyEdgeSharingCellsInCCW( cellsDefiningThisEdge );


    set<VDCell*> setOfCellsDefiningThisEdge;
    cellsDefiningThisEdge.reset4Loop();
    while ( cellsDefiningThisEdge.setNext4Loop() ) {
        setOfCellsDefiningThisEdge.insert( cellsDefiningThisEdge.getEntity() );
    }


    rg_dList<VDCell*> cellsIncidnetToEndVtx;
    m_endVertex->inquireIncidentCells( cellsIncidnetToEndVtx );

    VDCell* endCell = rg_NULL;

    cellsIncidnetToEndVtx.reset4Loop();
    while ( cellsIncidnetToEndVtx.setNext4Loop() ) {
        VDCell* currCell = cellsIncidnetToEndVtx.getEntity();

        if ( setOfCellsDefiningThisEdge.find( currCell ) == setOfCellsDefiningThisEdge.end() ) {
            endCell = currCell;
            break;
        }
    }

    return endCell;
}



VDPartialEdge* VDEdge::getPartialEdge() const
{
    return m_partEdge;
}


rg_Curve* VDEdge::getEdgeEquation() const
{
    return m_edgeEquation;
}


rg_INT VDEdge::getNumberOfPointsOnEdge() const
{
    return m_numOfPointsOnEdge;
}

rg_Point3D* VDEdge::getPointsOnEdge() const
{
    return m_pointsOnEdge;
}

rg_Point3D VDEdge::getPointOnEdge(const rg_INT& i) const
{
    return m_pointsOnEdge[i];
}

rg_Point3D VDEdge::getPointOnEdgeByEquation(const rg_REAL& t) const
{
    rg_RBzCurve3D* curve = (rg_RBzCurve3D*) m_edgeEquation;

    if ( rg_GE( t, 1.0) )
        return curve->evaluatePt( 1.0 );
    else if ( rg_LE( t, 0.0 ) )
        return curve->evaluatePt( 0.0 );
    else
        return curve->evaluatePt(t);
}

rg_REAL	VDEdge::getMinRadius() const
{
	return m_minRadius;
}

rg_REAL VDEdge::getMaxRadius() const
{
	if(m_maxPosition == START)
	{
		return m_startVertex->getRadiusOfTangentSphere();
	}
	else
	{
		if( m_endVertex->isOnInfinity() )
			return DBL_MAX;
		else
			return m_endVertex->getRadiusOfTangentSphere();
	}
}

PositionOnEdge VDEdge::getMinPosition() const
{
	return m_minPosition;
}

PositionOnEdge VDEdge::getMaxPosition() const
{
	return m_maxPosition;
}



rg_FLAG VDEdge::getEdgeType() const
{
	return m_edgeType;
}



Sphere VDEdge::getMinMaxTangentSphereByBalls(const rg_INT& i) const
{
    return m_minMaxTangentSphereByBalls[i];
}



Plane  VDEdge::getEdgePlane() const
{
    return m_edgePlane;
}



rg_Point3D VDEdge::getPointForAngleDistance() const
{
    return m_pointForAngleDist;
}



rg_FLAG VDEdge::isOnInfinity() const
{
    return m_isOnInfinity;
}

rg_FLAG VDEdge::isBoundedEdge() const
{
    if ( m_startVertex == rg_NULL && m_endVertex == rg_NULL )
        return rg_TRUE;

    if ( m_endVertex->isOnInfinity() && !m_startVertex->isOnInfinity() )
        return rg_FALSE;
    else if ( !m_endVertex->isOnInfinity() && m_startVertex->isOnInfinity() )
        return rg_FALSE;
    else
        return rg_TRUE;
}



rg_FLAG VDEdge::isBounded() const
{
    if ( m_startVertex->isOnInfinity()==rg_TRUE || m_endVertex->isOnInfinity()==rg_TRUE )
        return rg_FALSE;
    else 
        return rg_TRUE;
}



rg_FLAG VDEdge::isTangible() const
{
    return m_ProbeTangibility;
}

void    VDEdge::isTangible(const rg_FLAG& tangibility)
{
    m_ProbeTangibility = tangibility;
}

rg_FLAG	VDEdge::isMonotone() const
{
	if( m_minPosition == MID )
		return rg_FALSE;
	else
		return rg_TRUE;
}

    
//rg_BOOL VDEdge::isIncidentTo(VDFace* v_face) const
//{
//    rg_dList<VDEdge*> boundingEdges;
//    v_face->inquireBoundingEdges( boundingEdges );
//
//    rg_BOOL isFaceIncidentToThisEdge = rg_FALSE;
//    boundingEdges.reset4Loop();
//    while ( boundingEdges.setNext4Loop() ) {
//        VDEdge* currEdge = boundingEdges.getEntity();
//
//        if ( currEdge == this ) {
//            isFaceIncidentToThisEdge = rg_TRUE;
//            break;
//        }
//    }
//        
//    return isFaceIncidentToThisEdge;
//}
//

///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void VDEdge::setStartVertex( VDVertex* startVertex )
{
    m_startVertex = startVertex;
}

void VDEdge::setEndVertex( VDVertex* endVertex )
{
    m_endVertex = endVertex;
}


void VDEdge::setPartialEdge( VDPartialEdge* partialEdge )
{
    m_partEdge = partialEdge;
}

void VDEdge::setMinRadius(const rg_REAL& radius)
{
	m_minRadius = radius;
}

void VDEdge::setMinPosition(const PositionOnEdge& position)
{
	m_minPosition = position;
}

void VDEdge::setMaxPosition(const PositionOnEdge& position)
{
	m_maxPosition = position;
}

void VDEdge::setEdgeEquation( rg_Curve* edgeEquation )
{
    //  edgeEquation은 rg_Curve의 포인터이고, 
    //  다음은 객체의 복사가 아닌 단순히 포인터에 대한 복사만 이루어짐.
    //  rg_Curve에 대한 virtual function의 정의가 필요할 것으로 생각됨.
    //  임시로 Sphere set voronoi diagram 에 대한 edge equation 객체의 복사를 수행하게 함. 
    // 2004-4-6 by Youngsong Cho
    if ( edgeEquation != rg_NULL )
    {
        if ( m_edgeEquation != rg_NULL )
            delete m_edgeEquation;

        //  Modified by Youngsong Cho 2005. 5. 24
        //  m_edgeEquation = new rg_RBzCurve3D( *((rg_RBzCurve3D*) edgeEquation) );
        m_edgeEquation = edgeEquation;
    }
    else
    {
        m_edgeEquation = rg_NULL;
    }
    
}

void VDEdge::setPointsOnEdge( const rg_INT&     numPtsOnEdge,
                                    rg_Point3D* ptsOnEdge )
{
    if ( m_pointsOnEdge != rg_NULL )
        delete [] m_pointsOnEdge;

    m_numOfPointsOnEdge = numPtsOnEdge;

    if ( m_numOfPointsOnEdge == 0 )
    {
        m_pointsOnEdge = rg_NULL;
    }
    else 
    {
        //  Modified by Youngsong Cho 2005. 5. 24
        //  m_pointsOnEdge = new rg_Point3D[ m_numOfPointsOnEdge ];
        //  for (rg_INT i=0; i<m_numOfPointsOnEdge; i++)
        //      m_pointsOnEdge[i] = ptsOnEdge[i];
        m_pointsOnEdge = ptsOnEdge;
    }
}

void VDEdge::setInfinity()
{
	m_isOnInfinity = rg_TRUE;
}


void VDEdge::setEdge(       VDVertex*      startVertex, 
                            VDVertex*      endVertex, 
                            VDPartialEdge* partialEdge, 
                            rg_Curve*      edgeEquation, 
                      const rg_INT&        numPtsOnEdge, 
                            rg_Point3D*    ptsOnEdge )
{
    m_startVertex  = startVertex;
    m_endVertex    = endVertex;

    m_partEdge     = partialEdge;


    //  edgeEquation은 rg_Curve의 포인터이고, 
    //  다음은 객체의 복사가 아닌 단순히 포인터에 대한 복사만 이루어짐.
    //  rg_Curve에 대한 virtual function의 정의가 필요할 것으로 생각됨.
    //  임시로 Sphere set voronoi diagram 에 대한 edge equation 객체의 복사를 수행하게 함. 
    // 2004-4-6 by Youngsong Cho
    if ( edgeEquation != rg_NULL )
    {
        if ( m_edgeEquation != rg_NULL )
            delete m_edgeEquation;

        //  Modified by Youngsong Cho 2005. 5. 24
        //  m_edgeEquation = new rg_RBzCurve3D( *((rg_RBzCurve3D*) edgeEquation) );
        m_edgeEquation = edgeEquation;
    }
    else
    {
        m_edgeEquation = rg_NULL;
    }
    


    if ( m_pointsOnEdge != rg_NULL )
        delete [] m_pointsOnEdge;

    m_numOfPointsOnEdge = numPtsOnEdge;

    if ( m_numOfPointsOnEdge == 0 )
    {
        m_pointsOnEdge = rg_NULL;
    }
    else 
    {
        //  Modified by Youngsong Cho 2005. 5. 24
        //  m_pointsOnEdge = new rg_Point3D[ m_numOfPointsOnEdge ];
        //  for (rg_INT i=0; i<m_numOfPointsOnEdge; i++)
        //      m_pointsOnEdge[i] = ptsOnEdge[i];
        m_pointsOnEdge = ptsOnEdge;
    }
}

void VDEdge::setEdge( const rg_INT&        ID,      
                            VDVertex*      startVertex, 
                            VDVertex*      endVertex, 
                            VDPartialEdge* partialEdge, 
                            rg_Curve*      edgeEquation, 
                      const rg_INT&        numPtsOnEdge, 
                            rg_Point3D*    ptsOnEdge )
{
    m_ID           = ID;

    m_startVertex  = startVertex;
    m_endVertex    = endVertex;

    m_partEdge     = partialEdge;

    //  edgeEquation은 rg_Curve의 포인터이고, 
    //  다음은 객체의 복사가 아닌 단순히 포인터에 대한 복사만 이루어짐.
    //  rg_Curve에 대한 virtual function의 정의가 필요할 것으로 생각됨.
    //  임시로 Sphere set voronoi diagram 에 대한 edge equation 객체의 복사를 수행하게 함. 
    // 2004-4-6 by Youngsong Cho
    if ( edgeEquation != rg_NULL )
    {
        if ( m_edgeEquation != rg_NULL )
            delete m_edgeEquation;

        //  Modified by Youngsong Cho 2005. 5. 24
        //  m_edgeEquation = new rg_RBzCurve3D( *((rg_RBzCurve3D*) edgeEquation) );
        m_edgeEquation = edgeEquation;
    }
    else
    {
        m_edgeEquation = rg_NULL;
    }

    

    if ( m_pointsOnEdge != rg_NULL )
        delete [] m_pointsOnEdge;

    m_numOfPointsOnEdge = numPtsOnEdge;

    if ( m_numOfPointsOnEdge == 0 )
    {
        m_pointsOnEdge = rg_NULL;
    }
    else 
    {
        //  Modified by Youngsong Cho 2005. 5. 24
        //  m_pointsOnEdge = new rg_Point3D[ m_numOfPointsOnEdge ];
        //  for (rg_INT i=0; i<m_numOfPointsOnEdge; i++)
        //      m_pointsOnEdge[i] = ptsOnEdge[i];
        m_pointsOnEdge = ptsOnEdge;
    }
}



void VDEdge::setEdgeType(const rg_FLAG& edgeType)
{
	m_edgeType = edgeType;
}


void VDEdge::setMinMaxTangentSphereByBalls(Sphere* minMaxSphere)
{
    m_minMaxTangentSphereByBalls[0] = minMaxSphere[0];
    m_minMaxTangentSphereByBalls[1] = minMaxSphere[1];
}



void VDEdge::setMinMaxTangentSphereByBalls(const Sphere& sphere1, const Sphere& sphere2)
{
    m_minMaxTangentSphereByBalls[0] = sphere1;
    m_minMaxTangentSphereByBalls[1] = sphere2;
}



void VDEdge::setEdgePlane(const Plane& edgePlane)
{
    m_edgePlane = edgePlane;
}
    


void VDEdge::setPointForAngleDistance(const rg_Point3D& point)
{
    m_pointForAngleDist = point;
}
///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
VDEdge& VDEdge::operator =( const VDEdge& anEdge )
{
    if ( this == &anEdge )
        return *this;

    m_ID           = anEdge.m_ID;

    m_startVertex  = anEdge.m_startVertex;
    m_endVertex    = anEdge.m_endVertex;

    m_partEdge     = anEdge.m_partEdge;



    //  edgeEquation은 rg_Curve의 포인터이고, 
    //  다음은 객체의 복사가 아닌 단순히 포인터에 대한 복사만 이루어짐.
    //  rg_Curve에 대한 virtual function의 정의가 필요할 것으로 생각됨.
    //  임시로 Sphere set voronoi diagram 에 대한 edge equation 객체의 복사를 수행하게 함. 
    // 2004-4-6 by Youngsong Cho
    if ( anEdge.m_edgeEquation != rg_NULL )
    {
        if ( m_edgeEquation != rg_NULL )
            delete m_edgeEquation;

        m_edgeEquation = new rg_RBzCurve3D( *((rg_RBzCurve3D*) anEdge.m_edgeEquation) );
    }
    else
    {
        m_edgeEquation = rg_NULL;
    }


    if ( m_pointsOnEdge != rg_NULL )
        delete [] m_pointsOnEdge;

    m_numOfPointsOnEdge = anEdge.m_numOfPointsOnEdge;

    if ( m_numOfPointsOnEdge == 0 )
    {
        m_pointsOnEdge = rg_NULL;
    }
    else 
    {
        m_pointsOnEdge = new rg_Point3D[ m_numOfPointsOnEdge ];
        for (rg_INT i=0; i<m_numOfPointsOnEdge; i++)
            m_pointsOnEdge[i] = anEdge.m_pointsOnEdge[i];
    }

    m_isOnInfinity     = anEdge.m_isOnInfinity;
    m_ProbeTangibility = anEdge.m_ProbeTangibility;

	m_edgeType = anEdge.m_edgeType;
    m_minMaxTangentSphereByBalls[0] = anEdge.m_minMaxTangentSphereByBalls[0];
    m_minMaxTangentSphereByBalls[1] = anEdge.m_minMaxTangentSphereByBalls[1];
    m_edgePlane                     = anEdge.m_edgePlane;
    m_pointForAngleDist             = anEdge.m_pointForAngleDist;

    m_visited  = anEdge.m_visited;
    m_betaFace = anEdge.m_betaFace;

    return *this;
}

void VDEdge::fileOutVDEdge( ofstream& fout )
{
    fout << "e" << m_ID << "\t";


    fout << "v" << m_startVertex->getID() << "\t";
    fout << "v" << m_endVertex->getID() << "\t";
    fout << "pe" << m_partEdge->getID() << "\t";
    fout << "pe" << m_partEdge->getNextPartialEdgeInRadialCycle()->getID() << "\t";
    fout << "pe" << m_partEdge->getNextPartialEdgeInRadialCycle()->getNextPartialEdgeInRadialCycle()->getID() << "\t";

    if ( m_endVertex->isOnInfinity() || m_startVertex->isOnInfinity() )
    {
        fout << "(un-bounded)" << endl;
        return;
    }
    else
        fout << "(bounded)" << "\t";

    rg_RBzCurve3D* curve = (rg_RBzCurve3D*) m_edgeEquation; 
    rg_Point3D ctrlPt[3] = { curve->getCtrlPt(0), curve->getCtrlPt(1), curve->getCtrlPt(2) };
    
    fout << "Curve" << "\t";
    fout << ctrlPt[0].getX() << "\t";
    fout << ctrlPt[0].getY() << "\t";
    fout << ctrlPt[0].getZ() << "\t";
    fout << ctrlPt[1].getX() << "\t";
    fout << ctrlPt[1].getY() << "\t";
    fout << ctrlPt[1].getZ() << "\t";
    fout << ctrlPt[2].getX() << "\t";
    fout << ctrlPt[2].getY() << "\t";
    fout << ctrlPt[2].getZ() << "\t";

    fout << curve->getWeight(0) << "\t";
    fout << curve->getWeight(1) << "\t";
    fout << curve->getWeight(2) << endl;
}

///////////////////////////////////////////////////////////////////////////////
//
//  topological operators..
void VDEdge::inquireIntoOnlyEdgeSharingCellsInCCW( VDCell** onlyEdgeSharingCellsInCCW ) const
{
    rg_INT i=0;
    VDPartialEdge* currPartEdge = m_partEdge;

    do 
    {
        if ( currPartEdge->isRightOrientationInLoop() )
            onlyEdgeSharingCellsInCCW[i] = currPartEdge->getLoop()->getFace()->getRightCell();
        else
            onlyEdgeSharingCellsInCCW[i] = currPartEdge->getLoop()->getFace()->getLeftCell();

        i++;
        currPartEdge = currPartEdge->getNextPartialEdgeInRadialCycle();

    } while ( currPartEdge != m_partEdge );
}

void VDEdge::inquireIntoOnlyEdgeSharingCellsInCCW( rg_dList<VDCell*>& onlyEdgeSharingCellsInCCW ) const
{
    VDPartialEdge* currPartEdge = m_partEdge;

    do 
    {
        if ( currPartEdge->isRightOrientationInLoop() )
            onlyEdgeSharingCellsInCCW.addTail( currPartEdge->getLoop()->getFace()->getRightCell() );
        else
            onlyEdgeSharingCellsInCCW.addTail( currPartEdge->getLoop()->getFace()->getLeftCell() );

        currPartEdge = currPartEdge->getNextPartialEdgeInRadialCycle();

    } while ( currPartEdge != m_partEdge );
}

rg_dList<VDCell*>* VDEdge::inquireIntoOnlyEdgeSharingCellsInCCW()
{
    rg_dList<VDCell*>* onlyEdgeSharingCellsInCCW = new rg_dList<VDCell*>;

    VDPartialEdge* currPartEdge = m_partEdge;

    do 
    {
        if ( currPartEdge->isRightOrientationInLoop() )
            onlyEdgeSharingCellsInCCW->addTail( currPartEdge->getLoop()->getFace()->getRightCell() );
        else
            onlyEdgeSharingCellsInCCW->addTail( currPartEdge->getLoop()->getFace()->getLeftCell() );

        currPartEdge = currPartEdge->getNextPartialEdgeInRadialCycle();

    } while ( currPartEdge != m_partEdge );

    return onlyEdgeSharingCellsInCCW;
}

void VDEdge::inquireIncidentEdgesAtStartVertex( rg_dList<VDEdge*>& incidentEdgesAtStartVertex )
{
    for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
    {
        if ( this != m_startVertex->getIncidentEdge( i ) )
        {
            incidentEdgesAtStartVertex.addTail( m_startVertex->getIncidentEdge( i ) );
        }
    }
}

void VDEdge::inquireIncidentEdgesAtEndVertex( rg_dList<VDEdge*>& incidentEdgesAtEndVertex )
{
    for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
    {
        if ( this != m_endVertex->getIncidentEdge( i ) )
        {
            incidentEdgesAtEndVertex.addTail( m_endVertex->getIncidentEdge( i ) );
        }
    }
}




void  VDEdge::inquireBoundingVertices(rg_dList<VDVertex*>& vertexList)
{
    vertexList.add( m_startVertex );
    vertexList.add( m_endVertex );
}



void VDEdge::inquireIncidentEdges( rg_dList<VDEdge*>& edgeList )
{
    rg_INT i=0;
	for ( i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)  {
        if ( this != m_startVertex->getIncidentEdge( i ) )
            edgeList.addTail( m_startVertex->getIncidentEdge( i ) );
    }

    for ( i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)  {
        if ( this != m_endVertex->getIncidentEdge( i ) )  
            edgeList.addTail( m_endVertex->getIncidentEdge( i ) );
    }
}




void  VDEdge::inquireIncidentFaces(rg_dList<VDFace*>& faceList) const
{
    VDPartialEdge* currPrEdge = m_partEdge;

    do  {
        if ( currPrEdge->getLoop() != rg_NULL ) {
            faceList.add( currPrEdge->getLoop()->getFace() );
        }

        currPrEdge = currPrEdge->getNextPartialEdgeInRadialCycle();
    } while ( m_partEdge != currPrEdge );
}



void  VDEdge::inquireIncidentCells(rg_dList<VDCell*>& cellList) const
{
    VDPartialEdge* currPrEdge = m_partEdge;

    do {
        if ( currPrEdge->isRightOrientationInLoop() ) {
            cellList.addTail( currPrEdge->getLoop()->getFace()->getRightCell() );
        }
        else {
            cellList.addTail( currPrEdge->getLoop()->getFace()->getLeftCell() );
        }

        currPrEdge = currPrEdge->getNextPartialEdgeInRadialCycle();
    } while ( currPrEdge != rg_NULL && currPrEdge != m_partEdge );

    if ( cellList.getSize() < 3 )  {
        cellList.removeAll();

        cellList.addTail( m_partEdge->getLoop()->getFace()->getRightCell() );
        cellList.addTail( m_partEdge->getLoop()->getFace()->getLeftCell() );
    }

    //if ( m_partEdge->getNextPartialEdgeInRadialCycle() != rg_NULL ) {
    //    if ( m_partEdge->isRightOrientationInLoop() ) {
    //        cellList.addTail( m_partEdge->getLoop()->getFace()->getLeftCell() );
    //        cellList.addTail( m_partEdge->getLoop()->getFace()->getRightCell() );
    //    }
    //    else {
    //        cellList.addTail( m_partEdge->getLoop()->getFace()->getRightCell() );
    //        cellList.addTail( m_partEdge->getLoop()->getFace()->getLeftCell() );
    //    }

    //    VDPartialEdge* nextPrEdge = m_partEdge->getNextPartialEdgeInRadialCycle();
    //    if ( nextPrEdge->isRightOrientationInLoop() ) {
    //        cellList.addTail( m_partEdge->getLoop()->getFace()->getRightCell() );
    //    }
    //    else {
    //        cellList.addTail( m_partEdge->getLoop()->getFace()->getLeftCell() );
    //    }
    //}
    //else {
    //    cellList.addTail( m_partEdge->getLoop()->getFace()->getRightCell() );
    //    cellList.addTail( m_partEdge->getLoop()->getFace()->getLeftCell() );
    //}
}


rg_BOOL VDEdge::isBoundingCell( VDCell* currCell ) const
{
    rg_dList<VDCell*> cellList;
    inquireIncidentCells( cellList );

    if ( cellList.isInList( currCell ) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}


    
rg_BOOL VDEdge::isIncidentTo(   VDFace* currFace ) const
{
    rg_dList<VDFace*> faceList;
    inquireIncidentFaces( faceList );

    if ( faceList.isInList( currFace ) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}


    
VDCell* VDEdge::getMateCellOfStartVertex() const
{

    VDCell* cellsToDefineEndVtx[4] = { rg_NULL, rg_NULL, rg_NULL, rg_NULL };
    m_endVertex->inquireAllCellsToDefineThisVertex( cellsToDefineEndVtx );

    rg_dList<VDCell*> cellToDefineEdge;
    inquireIncidentCells( cellToDefineEdge );


    VDCell* mateCellOfStartVtx = rg_NULL;

    for ( rg_INT i=0; i<4; i++ ) {
        if ( !cellToDefineEdge.isInList( cellsToDefineEndVtx[i] ) ) {
            mateCellOfStartVtx = cellsToDefineEndVtx[i];
            break;
        }
    }

    return mateCellOfStartVtx;
}
    


VDCell* VDEdge::getMateCellOfEndVertex() const
{
    VDCell* cellsToDefineStartVtx[4] = { rg_NULL, rg_NULL, rg_NULL, rg_NULL };
    m_startVertex->inquireAllCellsToDefineThisVertex( cellsToDefineStartVtx );

    rg_dList<VDCell*> cellToDefineEdge;
    inquireIncidentCells( cellToDefineEdge );


    VDCell* mateCellOfEndVtx = rg_NULL;
    for ( rg_INT i=0; i<4; i++ ) {
        if ( !cellToDefineEdge.isInList( cellsToDefineStartVtx[i] ) ) {
            mateCellOfEndVtx = cellsToDefineStartVtx[i];
            break;
        }
    }

    return mateCellOfEndVtx;
}


//rg_INT VDEdge::searchIncidentCells(rg_dList<VDCell*>& cellList) const
//{
//    set <VDCell*> incidentCells;
//    
//    VDPartialEdge* currPartEdge = m_partEdge;
//    do {
//        VDLoop* currLoop = currPartEdge->getLoop();
//
//        if ( currLoop != rg_NULL ) {
//            VDFace* currFace = currLoop->getFace();
//
//        }
//        if ( currPartEdge->isRightOrientationInLoop() )
//            cellList.addTail( currPartEdge->getLoop()->getFace()->getRightCell() );
//        else
//            cellList.addTail( currPartEdge->getLoop()->getFace()->getLeftCell() );
//
//        currPartEdge = currPartEdge->getNextPartialEdgeInRadialCycle();
//
//    } while ( currPartEdge != m_partEdge );
//
//}



rg_INT VDEdge::getNumCellsToDefineThisEdge() const
{
    set<VDCell*> cellsDefineThisEdge;

    VDPartialEdge* currPrEdge  = m_partEdge;
    do {
        VDLoop* currLoop = currPrEdge->getLoop();
        if ( currLoop != rg_NULL ) {
            VDFace* currFace = currLoop->getFace();

            if ( currFace->getRightCell() != rg_NULL ) {
                cellsDefineThisEdge.insert( currFace->getRightCell() );
            }
            if ( currFace->getLeftCell() != rg_NULL ) {
                cellsDefineThisEdge.insert( currFace->getLeftCell() );
            }
        }

        currPrEdge = currPrEdge->getNextPartialEdgeInRadialCycle();
    } while ( currPrEdge != m_partEdge );

    return cellsDefineThisEdge.size();
}


    
rg_BOOL VDEdge::findPositionOfBlockingSphericalProbeFrom(VDVertex* startingVtx, 
                                                    const rg_REAL& radiusOfBlockingSphericalProbe,
                                                    Sphere& position ) const
{
    if ( m_betaFace == rg_NULL ) return rg_FALSE;
    rg_INT state = m_betaFace->getBoundingState( radiusOfBlockingSphericalProbe );
    if ( state==EXTERIOR_SIMPLEX || state==INTERIOR_SIMPLEX ) {
        return rg_FALSE;
    }


    VDCell* cellsToDefineEdge[3];
    inquireIntoOnlyEdgeSharingCellsInCCW( cellsToDefineEdge );
    Sphere gateBall[3] = { cellsToDefineEdge[0]->getGenerator()->getBall(), 
                           cellsToDefineEdge[1]->getGenerator()->getBall(), 
                           cellsToDefineEdge[2]->getGenerator()->getBall() };

    Sphere trimmingSphere[2];
	rg_INT numTrimmingSphere = computeSphereWithGivenRadiusAndTangentTo3Spheres(radiusOfBlockingSphericalProbe, gateBall, trimmingSphere);

    if ( numTrimmingSphere == 1 ) {
        position = trimmingSphere[0];
    }
    else {
        //  construct edge-plane
        Gate gateOfEdge( cellsToDefineEdge[0]->getGenerator(), cellsToDefineEdge[1]->getGenerator(), cellsToDefineEdge[2]->getGenerator(), 
                         getStartVCell()->getGenerator(), m_startVertex );
        gateOfEdge.obtainGeometricProperties();

        Plane      edgePlane( gateOfEdge.getNormalOfEdgePlane(), gateOfEdge.getMinTangentSphere().getCenter() );
    
    
        //  find point for new vertex to trim the edge.
        rg_Point3D startPoint = m_startVertex->getPoint();
        rg_Point3D axisPoint  = gateOfEdge.getAxisPoint(); 


        rg_REAL    angleToTrimmedSphere[2];
        for ( rg_INT i=0; i<numTrimmingSphere; i++ ) {
            angleToTrimmedSphere[i] = edgePlane.computeAngleInCCW(startPoint, axisPoint, trimmingSphere[i].getCenter() );
        }

        if ( startingVtx == m_startVertex ) {
            if ( angleToTrimmedSphere[0] < angleToTrimmedSphere[1] ) {
                position = trimmingSphere[0];
            }
            else {
                position = trimmingSphere[1];
            }
        }
        else {
            if ( angleToTrimmedSphere[1] < angleToTrimmedSphere[0] ) {
                position = trimmingSphere[1];
            }
            else {
                position = trimmingSphere[0];
            }
        }
    }

    return rg_TRUE;
}



void  VDEdge::connectBetaFace(BetaFace* b_face)
{
    m_betaFace = b_face;
}



void  VDEdge::disconnectBetaFace(BetaFace* b_face)
{
    if ( m_betaFace == b_face ) {
        m_betaFace = rg_NULL;
    }
}


