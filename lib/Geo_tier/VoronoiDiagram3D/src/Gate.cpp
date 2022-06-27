#include "Gate.h"
#include "FunctionsForVoronoiDiagram3D.h"
#include <float.h>
#include "rg_BallGenerator.h"
#include "rg_RelativeOp.h"
#include "VDVertex.h"
#include "rg_TMatrix3D.h"
#include "FunctionsForVoronoiDiagram3D.h"
using namespace V::GeometryTier;

#include "Plane.h"

///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
Gate::Gate()
: m_gateBallCCW1(   rg_NULL ), 
  m_gateBallCCW2(   rg_NULL ), 
  m_gateBallCCW3(   rg_NULL ), 
  m_startBall(      rg_NULL ),
  m_startVertex(    rg_NULL ), 
  m_nodePosInStack( rg_NULL ),
  m_edgeType(       rg_UNKNOWN),
  m_endBall(        rg_NULL),
  m_minAngularDistance(rg_MAX_REAL),
  m_maxValidAngle(-rg_MAX_REAL),
  m_numVertexBy4Balls(0),
  m_distanceBetVsAndPc(0.0)
{
    m_tolerance = resNeg10;
    m_numExcludedBalls = 0;
    m_excludedBall     = rg_NULL;

}



Gate::Gate( BallGenerator*   gateBallCCW1, 
            BallGenerator*   gateBallCCW2, 
            BallGenerator*   gateBallCCW3, 
            BallGenerator*   startBall,
            VDVertex* startVertex  )
: m_gateBallCCW1(   gateBallCCW1 ), 
  m_gateBallCCW2(   gateBallCCW2 ), 
  m_gateBallCCW3(   gateBallCCW3 ), 
  m_startBall(      startBall ),
  m_startVertex(    startVertex ), 
  m_nodePosInStack( rg_NULL ),
  m_edgeType(       rg_UNKNOWN),
  m_endBall(        rg_NULL),
  m_minAngularDistance(rg_MAX_REAL),
  m_maxValidAngle(-rg_MAX_REAL),
  m_numVertexBy4Balls(0),
  m_distanceBetVsAndPc(0.0)
{
    m_tolerance = resNeg10;
    m_numExcludedBalls = 0;
    m_excludedBall     = rg_NULL;
}



Gate::Gate( BallGenerator*         gateBallCCW1, 
            BallGenerator*         gateBallCCW2, 
            BallGenerator*         gateBallCCW3, 
            BallGenerator*         startBall,
            VDVertex*       startVertex, 
            rg_dNode<Gate>* nodePos )
: m_gateBallCCW1(   gateBallCCW1 ), 
  m_gateBallCCW2(   gateBallCCW2 ), 
  m_gateBallCCW3(   gateBallCCW3 ), 
  m_startBall(      startBall ),
  m_startVertex(    startVertex ), 
  m_nodePosInStack( nodePos ),
  m_edgeType(       rg_UNKNOWN),
  m_endBall(        rg_NULL),
  m_minAngularDistance(rg_MAX_REAL),
  m_maxValidAngle(-rg_MAX_REAL),
  m_numVertexBy4Balls(0),
  m_distanceBetVsAndPc(0.0)
{
    m_tolerance = resNeg10;
    m_numExcludedBalls = 0;
    m_excludedBall     = rg_NULL;
}



Gate::Gate( const Gate& aGate )
: m_gateBallCCW1(   aGate.m_gateBallCCW1 ), 
  m_gateBallCCW2(   aGate.m_gateBallCCW2 ), 
  m_gateBallCCW3(   aGate.m_gateBallCCW3 ), 
  m_startBall(      aGate.m_startBall ),
  m_startVertex(    aGate.m_startVertex ), 
  m_nodePosInStack( aGate.m_nodePosInStack ),
  m_edgeType(       aGate.m_edgeType),
  m_endBall(        aGate.m_endBall),
  m_emptySphereAtEndVertex( aGate.m_emptySphereAtEndVertex),
  m_minAngularDistance(aGate.m_minAngularDistance),
  m_maxValidAngle(aGate.m_maxValidAngle),
  m_normalOfCenterPlane(aGate.m_normalOfCenterPlane),
  m_normalOfEdgePlane(aGate.m_normalOfEdgePlane),
  m_axisPoint(aGate.m_axisPoint),
  m_vecAxisToStart(aGate.m_vecAxisToStart),
  m_numVertexBy4Balls(aGate.m_numVertexBy4Balls),
  m_centerPlane(aGate.m_centerPlane),
  m_distanceBetVsAndPc(aGate.m_distanceBetVsAndPc)
{
	//  2006. 8. 17  By Youngsong Cho
	m_minTangentSphere[0] = aGate.m_minTangentSphere[0];
	m_minTangentSphere[1] = aGate.m_minTangentSphere[1];

	m_endBallList.duplicateList(aGate.m_endBallList);
    m_tolerance = resNeg10;
    m_numExcludedBalls = aGate.m_numExcludedBalls;
    m_excludedBall     = new BallGenerator*[ m_numExcludedBalls ];
    for ( rg_INT i=0; i<m_numExcludedBalls; i++ )
        m_excludedBall[i] = aGate.m_excludedBall[i];
}



Gate::~Gate()
{
    if ( m_excludedBall != rg_NULL )
        delete [] m_excludedBall;
}

///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
BallGenerator* Gate::getGateBallCCW1() const
{
    return m_gateBallCCW1;
}

BallGenerator* Gate::getGateBallCCW2() const
{
    return m_gateBallCCW2;
}

BallGenerator* Gate::getGateBallCCW3() const
{
    return m_gateBallCCW3;
}

BallGenerator* Gate::getStartBall() const
{
    return m_startBall;
}


VDVertex* Gate::getStartVertex() const
{
    return m_startVertex;
}


rg_dNode<Gate>* Gate::getNodePosInStack() const
{
    return m_nodePosInStack;
}



rg_FLAG Gate::getEdgeType() const
{
	return m_edgeType;
}



BallGenerator* Gate::getEndBall() const
{
	return m_endBall;
}



Sphere  Gate::getEmptyTangentSphereAtEndVertex() const
{
	return m_emptySphereAtEndVertex;
}





rg_FLAG Gate::isTheSameGate( const Gate&     theOtherGate, 
                                   VDVertex* endVertexOfTheOtherGate  ) const
{
    if ( m_startVertex != endVertexOfTheOtherGate )
        return rg_FALSE;

    BallGenerator* thisGateConfiguration[3]   = { m_gateBallCCW1, m_gateBallCCW2, m_gateBallCCW3 };
    BallGenerator* targetGateConfiguration[3] = { theOtherGate.m_gateBallCCW1, theOtherGate.m_gateBallCCW2, theOtherGate.m_gateBallCCW3 };

    sort3PointersInAscendingPowers( thisGateConfiguration );
    sort3PointersInAscendingPowers( targetGateConfiguration );

    for (rg_INT i=0; i<3; i++)
    {
        if ( thisGateConfiguration[i] != targetGateConfiguration[i] )
            return rg_FALSE;
    }

    return rg_TRUE;
}

rg_FLAG Gate::isTheSameGate( BallGenerator*   gateBallCCW1, 
                             BallGenerator*   gateBallCCW2, 
                             BallGenerator*   gateBallCCW3, 
                             VDVertex* startVertex  ) const
{
    if ( m_startVertex != startVertex )
        return rg_FALSE;

    BallGenerator* thisGateConfiguration[3]   = { m_gateBallCCW1, m_gateBallCCW2, m_gateBallCCW3 };
    BallGenerator* targetGateConfiguration[3] = { gateBallCCW1, gateBallCCW2, gateBallCCW3 };

    sort3PointersInAscendingPowers( thisGateConfiguration );
    sort3PointersInAscendingPowers( targetGateConfiguration );

    for (rg_INT i=0; i<3; i++)
    {
        if ( thisGateConfiguration[i] != targetGateConfiguration[i] )
            return rg_FALSE;
    }

    return rg_TRUE;
}



//rg_dList<BallGenerator*>* Gate::getEndBallList() 
rg_dList< pair<BallGenerator*, Sphere> >* Gate::getEndBallList() 
{
    return &m_endBallList;
}



BallGenerator* Gate::getSingleEndBall() const
{
//    return m_endBallList.getFirstEntity();
    return m_endBallList.getFirstEntity().first;
}




//void Gate::addEndBall(BallGenerator* endBall)
//{
//    m_endBallList.add( endBall );
//}



void Gate::emptyEndBallList()
{
    m_endBallList.removeAll();
}



rg_INT Gate::getNumBallsTangentToSphere() const
{
    return m_endBallList.getSize() + 3;
}



BallGenerator** Gate::getBallsTangentToSphere()
{
    rg_INT          numBallsTangentToSphere = m_endBallList.getSize() + 3;
    BallGenerator** ballsTangentToSphere    = new BallGenerator*[ numBallsTangentToSphere ];
    ballsTangentToSphere[0] = m_gateBallCCW1;
    ballsTangentToSphere[1] = m_gateBallCCW3;
    ballsTangentToSphere[2] = m_gateBallCCW2;

    rg_INT i=3;
    m_endBallList.reset4Loop();
    while ( m_endBallList.setNext4Loop() )  {
        ballsTangentToSphere[i++] = m_endBallList.getEntity().first;
    }

    return ballsTangentToSphere;
}



rg_Point3D* Gate::getNormalizedTangentPointsBetweenBallsAndSphere()
{

    rg_INT      numBallsTangentToSphere = m_endBallList.getSize() + 3;
    rg_Point3D* ptsOnTangentSphere      = new rg_Point3D[ numBallsTangentToSphere ];
    rg_Point3D  centerOfEmptyTangentSphere = m_emptySphereAtEndVertex.getCenter();

    BallGenerator* gateBall[3] = {m_gateBallCCW1, m_gateBallCCW3, m_gateBallCCW2};
    rg_INT i = 0;
    for ( i=0; i<3; i++ ){
        //ptsOnTangentSphere[i] = gateBall[i]->getCenter() - m_emptySphereAtEndVertex.getCenter();
        ptsOnTangentSphere[i] = gateBall[i]->getCenter() - centerOfEmptyTangentSphere;
        ptsOnTangentSphere[i].normalize();

        //ptsOnTangentSphere[i] = ptsOnTangentSphere[i]*m_emptySphereAtEndVertex.getRadius();
    }

    pair<BallGenerator*, Sphere>* currEndBall = rg_NULL;
    m_endBallList.reset4Loop();
    while ( m_endBallList.setNext4Loop() )  {
        currEndBall = m_endBallList.getpEntity();

        //ptsOnTangentSphere[i] = currEndBall->first->getCenter() - currEndBall->second.getCenter();
        ptsOnTangentSphere[i] = currEndBall->first->getCenter() - centerOfEmptyTangentSphere;
        ptsOnTangentSphere[i].normalize();

        //ptsOnTangentSphere[i] = ptsOnTangentSphere[i]*(currEndBall->second.getRadius());

        i++;
    }

    return ptsOnTangentSphere;

    /*    
    rg_dList<BallGenerator*> ballsTangentToSphere;
    getBallsTangentToSphere( ballsTangentToSphere );

    rg_Point3D  centerOfEmptyTangentSphere = m_emptySphereAtEndVertex.getCenter();

    rg_INT      numBallsTangentToSphere = m_endBallList.getSize() + 3;
    rg_Point3D* ptsOnTangentSphere      = new rg_Point3D[ numBallsTangentToSphere ];

    rg_INT i = 0;
    BallGenerator* currBall = rg_NULL;
    ballsTangentToSphere.reset4Loop();
    while ( ballsTangentToSphere.setNext4Loop() )  {
        currBall = ballsTangentToSphere.getEntity();

        ptsOnTangentSphere[i] = currBall->getCenter() - centerOfEmptyTangentSphere;
        ptsOnTangentSphere[i].normalize();

        i++;
    }

    return ptsOnTangentSphere;
    */
}



void Gate::getBallsTangentToSphere(rg_dList<BallGenerator*>& ballsTangentToSphere)
{
    ballsTangentToSphere.addTail(m_gateBallCCW1);
    ballsTangentToSphere.addTail(m_gateBallCCW3);
    ballsTangentToSphere.addTail(m_gateBallCCW2);

    m_endBallList.reset4Loop();
    while ( m_endBallList.setNext4Loop() )  {
        ballsTangentToSphere.addTail( m_endBallList.getEntity().first );
    }
}



rg_FLAG Gate::getNumVerticesBy4Balls() const
{
    return m_numVertexBy4Balls;
}



rg_INT Gate::getNumExcludedBalls() const
{
    return m_numExcludedBalls;
}



BallGenerator** Gate::getExcludedBalls() const
{
    return m_excludedBall;
}



void Gate::setExcludedBalls(const rg_INT& numExcludedBalls )
{
    if ( m_excludedBall != rg_NULL )
        delete [] m_excludedBall;

    m_numExcludedBalls = numExcludedBalls;
    m_excludedBall = new BallGenerator*[ m_numExcludedBalls ];
    for ( rg_INT i=0; i<m_numExcludedBalls; i++ )
        m_excludedBall[i] = rg_NULL;
}



void Gate::setExcludedBalls( const rg_INT& numExcludedBalls, BallGenerator** excludedBalls)
{
    m_numExcludedBalls = numExcludedBalls;
    m_excludedBall = new BallGenerator*[ m_numExcludedBalls ];
    for ( rg_INT i=0; i<m_numExcludedBalls; i++ )
        m_excludedBall[i] = excludedBalls[i];
}



void Gate::setOneExcludedBall(const rg_INT& i, BallGenerator* oneExcludedBall)
{
    if ( m_excludedBall != rg_NULL )
        m_excludedBall[i] = oneExcludedBall;
}



void Gate::setCheckOfAllExcludedBalls(const rg_FLAG& check)
{
    for ( rg_INT i=0; i<m_numExcludedBalls; i++ )
        m_excludedBall[i]->m_checkForBucket = check;
}



///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void Gate::setGateBallCCW1( BallGenerator* gateBallCCW1 )
{
    m_gateBallCCW1 = gateBallCCW1;
}

void Gate::setGateBallCCW2( BallGenerator* gateBallCCW2 )
{
    m_gateBallCCW2 = gateBallCCW2;
}

void Gate::setGateBallCCW3( BallGenerator* gateBallCCW3 )
{
    m_gateBallCCW3 = gateBallCCW3;
}

void Gate::setStartBall( BallGenerator* startBall )
{
    m_startBall = startBall;
}


void Gate::setStartVertex( VDVertex* startVertex )
{
    m_startVertex = startVertex;
}

void Gate::setNodePosInStack( rg_dNode<Gate>* nodePos )
{
    m_nodePosInStack = nodePos;
}

void Gate::setGate( BallGenerator*         gateBallCCW1, 
                    BallGenerator*         gateBallCCW2, 
                    BallGenerator*         gateBallCCW3, 
                    BallGenerator*         startBall, 
                    VDVertex*       startVertex  )
{
    m_gateBallCCW1   = gateBallCCW1;
    m_gateBallCCW2   = gateBallCCW2;
    m_gateBallCCW3   = gateBallCCW3;

    m_startBall      = startBall;
    m_startVertex    = startVertex;
}


void Gate::setGate( BallGenerator*         gateBallCCW1, 
                    BallGenerator*         gateBallCCW2, 
                    BallGenerator*         gateBallCCW3, 
                    BallGenerator*         startBall, 
                    VDVertex*       startVertex, 
                    rg_dNode<Gate>* nodePos )
{
    m_gateBallCCW1   = gateBallCCW1;
    m_gateBallCCW2   = gateBallCCW2;
    m_gateBallCCW3   = gateBallCCW3;

    m_startBall      = startBall;
    m_startVertex    = startVertex;
    
    m_nodePosInStack = nodePos;
}



void Gate::setEdgeType(const rg_FLAG& edgeType)
{
	m_edgeType = edgeType;
}



void Gate::setEndBall(BallGenerator* endBall)
{
	m_endBall = endBall;
}



void Gate::setEmptyTangentSphereAtEndVertex(const Sphere& emptySphere)
{
	m_emptySphereAtEndVertex = emptySphere;
}




///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
Gate& Gate::operator =( const Gate& aGate )
{
    if ( this == &aGate )
        return *this;

    m_gateBallCCW1   = aGate.m_gateBallCCW1;
    m_gateBallCCW2   = aGate.m_gateBallCCW2;
    m_gateBallCCW3   = aGate.m_gateBallCCW3;

    m_startBall      = aGate.m_startBall;
    m_startVertex    = aGate.m_startVertex;

    m_nodePosInStack = aGate.m_nodePosInStack;

	m_edgeType       = aGate.m_edgeType;
	m_endBall        = aGate.m_endBall;
	m_emptySphereAtEndVertex = aGate.m_emptySphereAtEndVertex;
    m_numVertexBy4Balls = aGate.m_numVertexBy4Balls;

	//  2006. 8. 17  By Youngsong Cho
	m_emptySphereAtEndVertex = aGate.m_emptySphereAtEndVertex;
	m_minAngularDistance     = aGate.m_minAngularDistance;
	m_maxValidAngle          = aGate.m_maxValidAngle;
	m_normalOfCenterPlane    = aGate.m_normalOfCenterPlane;
	m_normalOfEdgePlane      = aGate.m_normalOfEdgePlane;
	m_axisPoint              = aGate.m_axisPoint;
	m_vecAxisToStart         = aGate.m_vecAxisToStart;

	m_minTangentSphere[0] = aGate.m_minTangentSphere[0];
	m_minTangentSphere[1] = aGate.m_minTangentSphere[1];

	m_endBallList.duplicateList(aGate.m_endBallList);
    m_tolerance = resNeg10;


    m_numExcludedBalls = aGate.m_numExcludedBalls;
    m_excludedBall     = new BallGenerator*[ m_numExcludedBalls ];
    for ( rg_INT i=0; i<m_numExcludedBalls; i++ )
        m_excludedBall[i] = aGate.m_excludedBall[i];
    

    m_centerPlane        = aGate.m_centerPlane;
    m_distanceBetVsAndPc = aGate.m_distanceBetVsAndPc;


    return *this;
}



//  2006. 8. 17  By Youngsong Cho
rg_INT  Gate::getNumOfEndBall() const
{
    return m_endBallList.getSize();
}


rg_Point3D Gate::getNormalOfCenterPlane() const
{
    return m_normalOfCenterPlane;
}



rg_REAL Gate::getCoeffOfCenterPlane() const
{
    return -(m_normalOfCenterPlane.innerProduct( m_minTangentSphere[0].getCenter() ));
}



Sphere Gate::getMinTangentSphere() const
{
    return m_minTangentSphere[0];
}



Sphere Gate::getMinTangentSphere(const rg_INT& i) const
{
    return m_minTangentSphere[i];
}



rg_Point3D Gate::getAxisPoint() const
{
    return m_axisPoint;
}


Plane Gate::getCenterPlane() const
{
    return m_centerPlane;
}



Plane   Gate::getEdgePlane() const
{
    Plane edgePlane(m_normalOfEdgePlane, m_axisPoint);
    return edgePlane;
}


void Gate::obtainGeometricProperties()
{
    Sphere gateBall[3] = { m_gateBallCCW1->getBall(), m_gateBallCCW2->getBall(), m_gateBallCCW3->getBall() }; 
    rg_Point3D centerOfGB[3] = { gateBall[0].getCenter(), gateBall[1].getCenter(), gateBall[2].getCenter() };
    rg_REAL    radiusOfGB[3] = { gateBall[0].getRadius(), gateBall[1].getRadius(), gateBall[2].getRadius() };



    ////////////////////////////////////////////////////
    //
    //  for test
    Sphere gateBallNew[3] = { m_gateBallCCW1->getBall(), m_gateBallCCW2->getBall(), m_gateBallCCW3->getBall() }; 

	rg_INT i=0;
	for ( i=0; i<2; i++ )  {
        if ( gateBallNew[i] < gateBallNew[i+1] )  {
            Sphere tempSphere = gateBallNew[i];
            gateBallNew[i]   = gateBallNew[i+1];
            gateBallNew[i+1] = tempSphere;
        }
    }
    if ( gateBallNew[0] < gateBallNew[1] )  {
        Sphere tempSphere = gateBallNew[0];
        gateBallNew[0] = gateBallNew[1];
        gateBallNew[1] = tempSphere;
    }
    
    rg_Point3D centerOfGBNew[3] = { gateBallNew[0].getCenter(), gateBallNew[1].getCenter(), gateBallNew[2].getCenter() };
    rg_REAL    radiusOfGBNew[3] = { gateBallNew[0].getRadius(), gateBallNew[1].getRadius(), gateBallNew[2].getRadius() };
    
    for ( i=0; i<3; i++ )  {
        radiusOfGBNew[i] = radiusOfGBNew[i] - radiusOfGBNew[2];
    }


    rg_REAL    halfVertexAngleOfCone = acos( (radiusOfGBNew[0]-radiusOfGBNew[1])/centerOfGBNew[0].distance(centerOfGBNew[1]) );

    rg_Point3D vec12New = centerOfGBNew[1] - centerOfGBNew[0];
    vec12New.normalize();  
    Plane      basePlane( vec12New, centerOfGBNew[0] );

    rg_REAL    distBetBasePlaneAndC3 = basePlane.distanceFromPoint(centerOfGBNew[2]);
    rg_Point3D projectedC3           = basePlane.projectPointOnPlane(centerOfGBNew[2]);

    rg_REAL distBetC1AndProjectedC3                    = centerOfGBNew[0].distance( projectedC3 );
    rg_REAL distBetprojectedC3AndSlantinglyProjectedC3 = distBetBasePlaneAndC3*tan( halfVertexAngleOfCone );
    rg_REAL distBetC1AndSlantinglyProjectedC3          = distBetC1AndProjectedC3 + distBetprojectedC3AndSlantinglyProjectedC3;
    rg_REAL distBetC1AndSlantinglyProjectedVertex      = radiusOfGBNew[0]*cos(halfVertexAngleOfCone);


    if ( distBetC1AndSlantinglyProjectedC3 >= distBetC1AndSlantinglyProjectedVertex )  {
    //  this edge is non-elliptic.
        m_centerPlane.definePlaneByThreePoints(centerOfGB[0], centerOfGB[1], centerOfGB[2]);
        m_distanceBetVsAndPc = m_centerPlane.distanceFromPoint( m_startVertex->getPoint() );
    }
    //
    ////////////////////////////////////////////////////




    //  If the gate balls define a linear edge,
    //  we use Euclidean distance between start and end vertices as the angular distance.
    if ( rg_EQ(radiusOfGB[0], radiusOfGB[1]) && rg_EQ(radiusOfGB[0], radiusOfGB[2]) && rg_EQ(radiusOfGB[1], radiusOfGB[2]) )
    {
	    rg_Point3D vec21 = centerOfGB[1] - centerOfGB[0];
	    rg_Point3D vec31 = centerOfGB[2] - centerOfGB[0];
	    m_normalOfCenterPlane = vec21.crossProduct(vec31);
        m_normalOfCenterPlane.normalize();
        
        rg_Point3D startVertex    = m_startVertex->getPoint();
        rg_REAL    distBetPcAndVs =   m_normalOfCenterPlane.innerProduct( startVertex )
                                    - m_normalOfCenterPlane.innerProduct( centerOfGB[0] );

        rg_Point3D center = startVertex - (distBetPcAndVs*m_normalOfCenterPlane);
        rg_REAL    radius = gateBall[0].distance( center );    
        m_minTangentSphere[0].setSphere( center, radius );

        m_axisPoint = centerOfGB[0];
        rg_Point3D vecPaCt = center - m_axisPoint;
        vecPaCt.normalize();
        m_normalOfEdgePlane = vecPaCt.crossProduct( m_normalOfCenterPlane );
        m_normalOfEdgePlane.normalize();

        m_edgeType = LINEAR_EDGE; 

        return;
    }

    ///////////////////////////////////////////////////////////////////////////
    //  1. find the normal of center plane.
    //    1.1 test if the configuration of centers of gate balls make a triangle.
    //        if a, b, and c (a<b<c) is the length of edges of a triangle,
    //        then c < a + b is always held.
    //    1.2 If the configuration of centers make a triangle,
    //        the gate balls may define a linear, parabolic, hyperbolic, or elliptic edge.
    //        Otherwise, the gate balls define a circular edge.

    rg_REAL length[3] = { centerOfGB[0].distance( centerOfGB[1] ),
                          centerOfGB[0].distance( centerOfGB[2] ),
                          centerOfGB[1].distance( centerOfGB[2] ) };

    //  normal of center plane for a linear, parabolic, hyperbolic, or elliptic edge.
    if ( rg_GT(length[0]+length[1], length[2]) && rg_GT(length[0]+length[2], length[1]) && rg_GT(length[1]+length[2], length[0]) )  {
	    rg_Point3D vec21 = centerOfGB[1] - centerOfGB[0];
	    rg_Point3D vec31 = centerOfGB[2] - centerOfGB[0];

	    m_normalOfCenterPlane = vec21.crossProduct(vec31);
        m_normalOfCenterPlane.normalize();
    }
    //  normal of center plane for a circular edge.
    else  {
        Sphere startBall = m_startBall->getBall();

	    rg_Point3D vecBg1Bg2 = centerOfGB[1]         - centerOfGB[0];
	    rg_Point3D vecBg1Bs  = startBall.getCenter() - centerOfGB[0];

	    m_normalOfCenterPlane = vecBg1Bg2.crossProduct(vecBg1Bs);
        m_normalOfCenterPlane.normalize();

        rg_Point3D vecBsBg1 = centerOfGB[0] - startBall.getCenter();
        rg_Point3D vecBsBg2 = centerOfGB[1] - startBall.getCenter();
        rg_Point3D vecBsBg3 = centerOfGB[2] - startBall.getCenter();

        vecBsBg1.normalize();
        vecBsBg2.normalize();
        vecBsBg3.normalize();
        rg_REAL angleCg1VsCg2 = angleCCWOfTwoUnitVectors( m_normalOfCenterPlane, 
                                                          vecBsBg1, vecBsBg2 );
        rg_REAL angleCg1VsCg3 = angleCCWOfTwoUnitVectors( m_normalOfCenterPlane, 
                                                          vecBsBg1, vecBsBg3 );

        if ( angleCg1VsCg2 > angleCg1VsCg3  )  
            m_normalOfCenterPlane = (-m_normalOfCenterPlane);
    }

    
    ///////////////////////////////////////////////////////////////////////////
    //
    //  find the intersection between three gate balls and center pland
    //    for the center of tangent circle, axis point, and normal of edge plane.
    //

    //  project gate balls into center plane.
	rg_TMatrix3D trMatrix;
	rg_TMatrix3D invMatrix;

	trMatrix.translate(-centerOfGB[0]);

    rg_FLAG bReverse = rg_EQ( m_normalOfCenterPlane.innerProduct( rg_Point3D(0.0, 0.0, 1.0) ), -1);
    if( !bReverse )
	    trMatrix.rotate(m_normalOfCenterPlane, rg_Point3D(0.0, 0.0, 1.0));
    else
        trMatrix.rotateY(rg_PI);

    if( !bReverse )
	    invMatrix.rotate(rg_Point3D(0.0, 0.0, 1.0), m_normalOfCenterPlane);
    else
        invMatrix.rotateY(-rg_PI);
	invMatrix.translate( centerOfGB[0] );


    rg_Point3D  trGateCenter[3] = { trMatrix*centerOfGB[0], trMatrix*centerOfGB[1], trMatrix*centerOfGB[2] };
	rg_Circle2D transformedCircle1( trGateCenter[0].getX(), trGateCenter[0].getY(), radiusOfGB[0] );
	rg_Circle2D transformedCircle2( trGateCenter[1].getX(), trGateCenter[1].getY(), radiusOfGB[1] );
	rg_Circle2D transformedCircle3( trGateCenter[2].getX(), trGateCenter[2].getY(), radiusOfGB[2] );

	rg_INT numOfCC = 0;
	rg_Circle2D tangentCircle[2];
    numOfCC = computeCircleTangentTo3CirclesOutside_GLOBAL( transformedCircle1, transformedCircle2, 
                                                            transformedCircle3, tangentCircle);   

    rg_Point3D centerOfTangentCircle;

    // The gate balls define a parabolic, or hyperbolic edge.
    if(numOfCC == 1)  {      
        ///////////////////////////////////////////////////////////////////////////
        //  2. find the center of tangent circle.

        centerOfTangentCircle = invMatrix*tangentCircle[0].getCenterPt();

        m_minTangentSphere[0].setSphere( centerOfTangentCircle, tangentCircle[0].getRadius() );

        
        ///////////////////////////////////////////////////////////////////////////
        //  3. find the axis point.

        rg_Point3D startVertex = m_startVertex->getPoint();
        if ( centerOfTangentCircle != startVertex )  {
            rg_REAL distanceVsAndCP = m_normalOfCenterPlane.innerProduct( startVertex )
                                      - m_normalOfCenterPlane.innerProduct( centerOfTangentCircle );
            m_axisPoint = startVertex - (distanceVsAndCP*m_normalOfCenterPlane);
        }
        else  {
            Sphere anotherTangentSphere[2];
            computeSphereWithGivenRadiusAndTangentTo3Spheres( 
                                                2.*m_minTangentSphere[0].getRadius(),
                                                gateBall, anotherTangentSphere );

            rg_REAL distanceViAndCP = m_normalOfCenterPlane.innerProduct( anotherTangentSphere[0].getCenter() )
                                      - m_normalOfCenterPlane.innerProduct( centerOfTangentCircle );

            m_axisPoint = anotherTangentSphere[0].getCenter() - (distanceViAndCP*m_normalOfCenterPlane);           
        }


        ///////////////////////////////////////////////////////////////////////////
        //  4. find the normal of edge plane.

        rg_Point3D vecPaCt = centerOfTangentCircle - m_axisPoint;
        vecPaCt.normalize();

        m_normalOfEdgePlane = vecPaCt.crossProduct( m_normalOfCenterPlane );
        m_normalOfEdgePlane.normalize();

		m_vecAxisToStart = startVertex - m_axisPoint;
		m_vecAxisToStart.normalize();    

		rg_Point3D vecCtPa = m_axisPoint - m_minTangentSphere[0].getCenter();
		vecCtPa.normalize();

		m_maxValidAngle = angleCCWOfTwoUnitVectors( m_normalOfEdgePlane, m_vecAxisToStart, vecCtPa );


        m_edgeType = PARABOLIC_OR_HYPERBOLIC_EDGE;	
	}
	else if(numOfCC == 2) // The gate defines an elliptic edge.
	{

    
    
        ///////////////////////////////////////////////////////////////////////////
        //
        //  2. find the center of tangent circle.
        //
        rg_Point2D vecGc1Tc1 = transformedCircle1.getCenterPt() - tangentCircle[0].getCenterPt();
        rg_Point2D vecGc2Tc1 = transformedCircle2.getCenterPt() - tangentCircle[0].getCenterPt();
        rg_Point2D vecGc3Tc1 = transformedCircle3.getCenterPt() - tangentCircle[0].getCenterPt();

        rg_REAL angleGc1Tc1Gc2 = angleCCW( vecGc1Tc1, vecGc2Tc1 );
        rg_REAL angleGc1Tc1Gc3 = angleCCW( vecGc1Tc1, vecGc3Tc1 );

        rg_Point3D theOtherCenter;
        if ( rg_LT( angleGc1Tc1Gc2, angleGc1Tc1Gc3 ) )  {
            centerOfTangentCircle = invMatrix*tangentCircle[0].getCenterPt();
            theOtherCenter        = invMatrix*tangentCircle[1].getCenterPt();

            m_minTangentSphere[0].setSphere( centerOfTangentCircle, tangentCircle[0].getRadius() );
            m_minTangentSphere[1].setSphere( theOtherCenter, tangentCircle[1].getRadius() );

        }
        else {
            centerOfTangentCircle = invMatrix*tangentCircle[1].getCenterPt();
            theOtherCenter        = invMatrix*tangentCircle[0].getCenterPt();

            m_minTangentSphere[0].setSphere( centerOfTangentCircle, tangentCircle[1].getRadius() );
            m_minTangentSphere[1].setSphere( theOtherCenter, tangentCircle[0].getRadius() );
        }


        ///////////////////////////////////////////////////////////////////////////
        //
        //  3. find the axis point.
        //
        m_axisPoint = ( centerOfTangentCircle + theOtherCenter ) / 2.0;

        ///////////////////////////////////////////////////////////////////////////
        //
        //  4. find the normal of edge plane.
        //
        rg_Point3D vecPaCt = centerOfTangentCircle - m_axisPoint;

        m_normalOfEdgePlane = vecPaCt.crossProduct( m_normalOfCenterPlane );
        m_normalOfEdgePlane.normalize();
		
		m_vecAxisToStart = m_startVertex->getPoint() - m_axisPoint;
		m_vecAxisToStart.normalize();    


        m_edgeType = CIRCULAR_OR_ELLIPTIC_EDGE;	
	}
    else 
    {
        m_edgeType = NON_EDGE;	
    }    
}



void Gate::determineEdgeType()
{
    Sphere     gateBall[3]   = { m_gateBallCCW1->getBall(), m_gateBallCCW2->getBall(), m_gateBallCCW3->getBall() }; 
    rg_Point3D centerOfGB[3] = { gateBall[0].getCenter(), gateBall[1].getCenter(), gateBall[2].getCenter() };
    rg_REAL    radiusOfGB[3] = { gateBall[0].getRadius(), gateBall[1].getRadius(), gateBall[2].getRadius() };



    Sphere gateBallSorted[3] = { gateBall[0], gateBall[1], gateBall[2] }; 

	rg_INT i=0;
    for ( i=0; i<2; i++ )  {
        if ( gateBallSorted[i] < gateBallSorted[i+1] )  {
            Sphere tempSphere   = gateBallSorted[i];
            gateBallSorted[i]   = gateBallSorted[i+1];
            gateBallSorted[i+1] = tempSphere;
        }
    }
    if ( gateBallSorted[0] < gateBallSorted[1] )  {
        Sphere tempSphere = gateBallSorted[0];
        gateBallSorted[0] = gateBallSorted[1];
        gateBallSorted[1] = tempSphere;
    }
    
    rg_Point3D centerOfSortedGB[3] = { gateBallSorted[0].getCenter(), gateBallSorted[1].getCenter(), gateBallSorted[2].getCenter() };
    rg_REAL    radiusOfSortedGB[3] = { gateBallSorted[0].getRadius(), gateBallSorted[1].getRadius(), gateBallSorted[2].getRadius() };
    for ( i=0; i<3; i++ )  
        radiusOfSortedGB[i] = radiusOfSortedGB[i] - radiusOfSortedGB[2];


    rg_REAL    diffBetR1AndR2        = radiusOfSortedGB[0]-radiusOfSortedGB[1];
    rg_REAL    distBetC1AndC2        = centerOfSortedGB[0].distance(centerOfSortedGB[1]);
    rg_REAL    halfVertexAngleOfCone = asin( diffBetR1AndR2/distBetC1AndC2 );

    rg_Point3D vec12 = centerOfSortedGB[1] - centerOfSortedGB[0];
    vec12.normalize();  
    Plane      basePlane( vec12, centerOfSortedGB[0] );

    rg_FLAG isNonElliptic = rg_TRUE;

    rg_REAL    distBetBasePlaneAndC3 = basePlane.distanceFromPoint(centerOfSortedGB[2]);

    if ( distBetBasePlaneAndC3 < 0.0 )
        isNonElliptic = rg_TRUE;
    else {   
        rg_Point3D projectedC3 = basePlane.projectPointOnPlane(centerOfSortedGB[2]);

        rg_REAL distBetC1AndProjectedC3                    = centerOfSortedGB[0].distance( projectedC3 );
        rg_REAL distBetprojectedC3AndSlantinglyProjectedC3 = distBetBasePlaneAndC3*tan( halfVertexAngleOfCone );
        rg_REAL distBetC1AndSlantinglyProjectedC3          = distBetC1AndProjectedC3 + distBetprojectedC3AndSlantinglyProjectedC3;
        rg_REAL distBetC1AndSlantinglyProjectedVertex      = radiusOfSortedGB[0]*cos(halfVertexAngleOfCone);

        if ( distBetC1AndSlantinglyProjectedC3 >= distBetC1AndSlantinglyProjectedVertex )  
            isNonElliptic = rg_TRUE;
        else
            isNonElliptic = rg_FALSE;
    }



    if ( isNonElliptic )  {
    //  this edge is non-elliptic.
        m_centerPlane.definePlaneByThreePoints(centerOfGB[0], centerOfGB[1], centerOfGB[2]);
        m_distanceBetVsAndPc = m_centerPlane.distanceFromPoint( m_startVertex->getPoint() );
        m_normalOfCenterPlane = m_centerPlane.getNormal();
        
        if ( rg_EQ(radiusOfSortedGB[0], radiusOfSortedGB[2]) )
            m_edgeType = LINEAR_EDGE; 
        else
            m_edgeType = PARABOLIC_OR_HYPERBOLIC_EDGE;	
    }
    else {
    //  this edge is elliptic.
        ///////////////////////////////////////////////////////////////////////////
        //  1. find the normal of center plane.
        //    1.1 test if the configuration of centers of gate balls make a triangle.
        //        if a, b, and c (a<b<c) is the length of edges of a triangle,
        //        then c < a + b is always held.
        //    1.2 If the configuration of centers make a triangle,
        //        the gate balls may define a linear, parabolic, hyperbolic, or elliptic edge.
        //        Otherwise, the gate balls define a circular edge.

        rg_REAL length[3] = { centerOfGB[0].distance( centerOfGB[1] ),
                              centerOfGB[0].distance( centerOfGB[2] ),
                              centerOfGB[1].distance( centerOfGB[2] ) };
        
        //  normal of center plane for an elliptic edge.
        if ( rg_GT(length[0]+length[1], length[2]) && rg_GT(length[0]+length[2], length[1]) && rg_GT(length[1]+length[2], length[0]) )  {
	        rg_Point3D vec21 = centerOfGB[1] - centerOfGB[0];
	        rg_Point3D vec31 = centerOfGB[2] - centerOfGB[0];

	        m_normalOfCenterPlane = vec21.crossProduct(vec31);
            m_normalOfCenterPlane.normalize();
        }
        //  normal of center plane for a circular edge.
        else  {
            Sphere startBall = m_startBall->getBall();

	        rg_Point3D vecBg1Bg2 = centerOfGB[1]         - centerOfGB[0];
	        rg_Point3D vecBg1Bs  = startBall.getCenter() - centerOfGB[0];

	        m_normalOfCenterPlane = vecBg1Bg2.crossProduct(vecBg1Bs);
            m_normalOfCenterPlane.normalize();

            rg_Point3D vecBsBg1 = centerOfGB[0] - startBall.getCenter();
            rg_Point3D vecBsBg2 = centerOfGB[1] - startBall.getCenter();
            rg_Point3D vecBsBg3 = centerOfGB[2] - startBall.getCenter();

            vecBsBg1.normalize();
            vecBsBg2.normalize();
            vecBsBg3.normalize();
            rg_REAL angleCg1VsCg2 = angleCCWOfTwoUnitVectors( m_normalOfCenterPlane, 
                                                              vecBsBg1, vecBsBg2 );
            rg_REAL angleCg1VsCg3 = angleCCWOfTwoUnitVectors( m_normalOfCenterPlane, 
                                                              vecBsBg1, vecBsBg3 );

            if ( angleCg1VsCg2 > angleCg1VsCg3  )  
                m_normalOfCenterPlane = (-m_normalOfCenterPlane);
        }


        ///////////////////////////////////////////////////////////////////////////
        //
        //  find the intersection between three gate balls and center pland
        //    for the center of tangent circle, axis point, and normal of edge plane.
        //
        //  project gate balls into center plane.
	    rg_TMatrix3D trMatrix;
	    rg_TMatrix3D invMatrix;

	    trMatrix.translate(-centerOfGB[0]);

        rg_FLAG bReverse = rg_EQ( m_normalOfCenterPlane.innerProduct( rg_Point3D(0.0, 0.0, 1.0) ), -1);
        if( !bReverse )
	        trMatrix.rotate(m_normalOfCenterPlane, rg_Point3D(0.0, 0.0, 1.0));
        else
            trMatrix.rotateY(rg_PI);

        if( !bReverse )
	        invMatrix.rotate(rg_Point3D(0.0, 0.0, 1.0), m_normalOfCenterPlane);
        else
            invMatrix.rotateY(-rg_PI);
	    invMatrix.translate( centerOfGB[0] );


        rg_Point3D  trGateCenter[3] = { trMatrix*centerOfGB[0], trMatrix*centerOfGB[1], trMatrix*centerOfGB[2] };
	    rg_Circle2D transformedCircle1( trGateCenter[0].getX(), trGateCenter[0].getY(), radiusOfGB[0] );
	    rg_Circle2D transformedCircle2( trGateCenter[1].getX(), trGateCenter[1].getY(), radiusOfGB[1] );
	    rg_Circle2D transformedCircle3( trGateCenter[2].getX(), trGateCenter[2].getY(), radiusOfGB[2] );

	    rg_INT numOfCC = 0;
	    rg_Circle2D tangentCircle[2];
        numOfCC = computeCircleTangentTo3CirclesOutside_GLOBAL( transformedCircle1, transformedCircle2, 
                                                                transformedCircle3, tangentCircle);   

        rg_Point3D centerOfTangentCircle;


        ///////////////////////////////////////////////////////////////////////////
        //
        //  2. find the center of tangent circle.
        //
        rg_Point2D vecGc1Tc1 = transformedCircle1.getCenterPt() - tangentCircle[0].getCenterPt();
        rg_Point2D vecGc2Tc1 = transformedCircle2.getCenterPt() - tangentCircle[0].getCenterPt();
        rg_Point2D vecGc3Tc1 = transformedCircle3.getCenterPt() - tangentCircle[0].getCenterPt();

        rg_REAL angleGc1Tc1Gc2 = angleCCW( vecGc1Tc1, vecGc2Tc1 );
        rg_REAL angleGc1Tc1Gc3 = angleCCW( vecGc1Tc1, vecGc3Tc1 );

        rg_Point3D theOtherCenter;
        if ( rg_LT( angleGc1Tc1Gc2, angleGc1Tc1Gc3 ) )  {
            centerOfTangentCircle = invMatrix*tangentCircle[0].getCenterPt();
            theOtherCenter        = invMatrix*tangentCircle[1].getCenterPt();

            m_minTangentSphere[0].setSphere( centerOfTangentCircle, tangentCircle[0].getRadius() );
            m_minTangentSphere[1].setSphere( theOtherCenter, tangentCircle[1].getRadius() );

        }
        else {
            centerOfTangentCircle = invMatrix*tangentCircle[1].getCenterPt();
            theOtherCenter        = invMatrix*tangentCircle[0].getCenterPt();

            m_minTangentSphere[0].setSphere( centerOfTangentCircle, tangentCircle[1].getRadius() );
            m_minTangentSphere[1].setSphere( theOtherCenter, tangentCircle[0].getRadius() );
        }


        ///////////////////////////////////////////////////////////////////////////
        //
        //  3. find the axis point.
        //
        m_axisPoint = ( centerOfTangentCircle + theOtherCenter ) / 2.0;

        ///////////////////////////////////////////////////////////////////////////
        //
        //  4. find the normal of edge plane.
        //
        rg_Point3D vecPaCt = centerOfTangentCircle - m_axisPoint;

        m_normalOfEdgePlane = vecPaCt.crossProduct( m_normalOfCenterPlane );
        m_normalOfEdgePlane.normalize();
		
		m_vecAxisToStart = m_startVertex->getPoint() - m_axisPoint;
		m_vecAxisToStart.normalize();    


        m_edgeType = CIRCULAR_OR_ELLIPTIC_EDGE;	
    }
}



rg_REAL Gate::angleCCWOfTwoUnitVectors(const rg_Point3D& normal, 
                                       const rg_Point3D& vector1, 
                                       const rg_Point3D& vector2)
{
	rg_REAL angle = vector1.innerProduct(vector2);

	if ( angle > 1.0 ) 
		angle = 0.0;
	else if ( angle < -1.0 )
		angle = 2.0;
    else 
		angle = 1.0 - angle;


    rg_Point3D cross = vector1.crossProduct( vector2 );
    cross.normalize();

    rg_REAL orientation = normal.innerProduct( cross );
	if( rg_POS( orientation ) )  //less than PI (cross product)
		return angle;
	else
		return (4.0 - angle);

}



rg_REAL Gate::angleCCW(const rg_Point2D& vector1, const rg_Point2D& vector2)
{
	rg_Point2D vec1 = vector1.getUnitVector();
	rg_Point2D vec2 = vector2.getUnitVector();
	rg_Point2D pt2pt( vec1 - vec2 );
	rg_REAL length = pt2pt.magnitude();

	rg_REAL cosine = (2. - length*length)/(2.);

	if( vec1 * vec2 > 0.0 )  //less than PI (cross product)
	{
		return acos( cosine );
	}
	else
	{
		return 2.*rg_PI - acos( cosine );
	}
}



void Gate::findEmptyTangentSphereWithMinAngleForParabolicEdge(rg_dList<BallGenerator*>& candidateBalls)
{
    Sphere gateBall[3] = { m_gateBallCCW1->getBall(), m_gateBallCCW2->getBall(), m_gateBallCCW3->getBall() }; 

    BallGenerator* currBall = rg_NULL;
    candidateBalls.reset4Loop();
	while( candidateBalls.setNext4Loop() )	{
	    currBall = candidateBalls.getEntity();

        Sphere  emptyTangentSphere[2];
        rg_INT  numTS = computeSphereTangentTo4SpheresOutside( gateBall[0], gateBall[1], gateBall[2], 
                                                               currBall->getBall(), emptyTangentSphere);


        for ( rg_INT i=0; i<numTS; i++ )  
        {
            if ( currBall == m_startBall )  {
                if ( numTS == 1 )
                    continue;

                if ( emptyTangentSphere[i].getCenter() == m_startVertex->getPoint() )
                    continue;
            }

            //  Modified by Youngsong Cho, 2006. 10. 14
            ///*
            rg_REAL coffCenterPlane    = m_normalOfCenterPlane.innerProduct(m_minTangentSphere[0].getCenter());
            rg_REAL distanceBetVsAndPc = m_normalOfCenterPlane.innerProduct(m_startVertex->getPoint()) - coffCenterPlane;
            rg_REAL distanceBetVeAndPc = m_normalOfCenterPlane.innerProduct(emptyTangentSphere[i].getCenter()) - coffCenterPlane;
            rg_REAL angularDistance    = distanceBetVeAndPc - distanceBetVsAndPc;

            if ( angularDistance < 0.0 )
                continue;

            if ( m_endBallList.getSize() > 0 )  {            
                // v_i:    emptyTangentSphere[i].getCenter()
                // v_star: current candidate end vertex
                rg_REAL distanceBetViAndVstar = m_emptySphereAtEndVertex.getCenter().distance( emptyTangentSphere[i].getCenter() );
            
                if ( rg_ZERO(distanceBetViAndVstar, m_tolerance) )  {
                    m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                    if ( emptyTangentSphere[i].getRadius() > m_emptySphereAtEndVertex.getRadius())
                        m_emptySphereAtEndVertex = emptyTangentSphere[i];
                }
                else  { 
                    if ( angularDistance < m_minAngularDistance ) {
                        m_minAngularDistance = angularDistance;
                        m_endBallList.removeAll();

                        m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                        m_emptySphereAtEndVertex = emptyTangentSphere[i];
                        m_numVertexBy4Balls      = numTS;
                    }
                }
            }
            else  {
                m_minAngularDistance = angularDistance;

                m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                m_emptySphereAtEndVertex = emptyTangentSphere[i];
                m_numVertexBy4Balls      = numTS;
            }
        }
    }
}



void Gate::findEmptyTangentSphereWithMinAngleForLinearEdge(rg_dList<BallGenerator*>& candidateBalls)
{
    Sphere gateBall[3] = { m_gateBallCCW1->getBall(), m_gateBallCCW2->getBall(), m_gateBallCCW3->getBall() }; 

    BallGenerator* currBall = rg_NULL;
    candidateBalls.reset4Loop();
	while( candidateBalls.setNext4Loop() )	{
	    currBall = candidateBalls.getEntity();

        Sphere  emptyTangentSphere[2];
        rg_INT  numTS = computeSphereTangentTo4SpheresOutside( gateBall[0], gateBall[1], gateBall[2], 
                                                               currBall->getBall(), emptyTangentSphere);


        for ( rg_INT i=0; i<numTS; i++ )  
        {
            if ( currBall == m_startBall )  {
                if ( numTS == 1 )
                    continue;

                if ( emptyTangentSphere[i].getCenter() == m_startVertex->getPoint() )
                    continue;
            }

            //  Modified by Youngsong Cho, 2006. 10. 14
            ///*
            rg_REAL coffCenterPlane    = m_normalOfCenterPlane.innerProduct(m_minTangentSphere[0].getCenter());
            rg_REAL distanceBetVsAndPc = m_normalOfCenterPlane.innerProduct(m_startVertex->getPoint()) - coffCenterPlane;
            rg_REAL distanceBetVeAndPc = m_normalOfCenterPlane.innerProduct(emptyTangentSphere[i].getCenter()) - coffCenterPlane;
            rg_REAL angularDistance    = distanceBetVeAndPc - distanceBetVsAndPc;

            if ( angularDistance < 0.0 )
                continue;

            if ( m_endBallList.getSize() > 0 )  {            
                // v_i:    emptyTangentSphere[i].getCenter()
                // v_star: current candidate end vertex
                rg_REAL distanceBetViAndVstar = m_emptySphereAtEndVertex.getCenter().distance( emptyTangentSphere[i].getCenter() );
            
                if ( rg_ZERO(distanceBetViAndVstar, m_tolerance) )  {
                    m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                    if ( emptyTangentSphere[i].getRadius() > m_emptySphereAtEndVertex.getRadius())
                        m_emptySphereAtEndVertex = emptyTangentSphere[i];
                }
                else  { 
                    if ( angularDistance < m_minAngularDistance ) {
                        m_minAngularDistance = angularDistance;
                        m_endBallList.removeAll();

                        m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                        m_emptySphereAtEndVertex = emptyTangentSphere[i];
                        m_numVertexBy4Balls      = numTS;
                    }
                }
            }
            else  {
                m_minAngularDistance = angularDistance;

                m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                m_emptySphereAtEndVertex = emptyTangentSphere[i];
                m_numVertexBy4Balls      = numTS;
            }
        }

        /*
        for ( rg_INT i=0; i<numTS; i++ )  {
            if ( currBall == m_startBall )  {
                if ( numTS == 1 )
                    continue;

                if ( emptyTangentSphere[i].getCenter() == m_startVertex->getPoint() )
                    continue;
            }

            rg_Point3D vecAxisToCurrEnd = emptyTangentSphere[i].getCenter() - m_startVertex->getPoint();
            rg_REAL    orientation = m_normalOfCenterPlane.innerProduct( vecAxisToCurrEnd.getUnitVector() );
            if ( !rg_POS( orientation ) )
                continue;

            rg_REAL    angularDistance = vecAxisToCurrEnd.magnitude();
            if ( rg_EQ(angularDistance, m_minAngularDistance, m_tolerance) )  {
                //m_endBallList.add( currBall );
                m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                if ( emptyTangentSphere[i].getRadius() > m_emptySphereAtEndVertex.getRadius())
                    m_emptySphereAtEndVertex = emptyTangentSphere[i];
            }
            else if ( rg_LT(angularDistance, m_minAngularDistance, m_tolerance) ) {
                m_minAngularDistance = angularDistance;
                m_endBallList.removeAll();
                //m_endBallList.add( currBall );
                m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );

                m_emptySphereAtEndVertex = emptyTangentSphere[i];
                m_numVertexBy4Balls      = numTS;
            }
            else {
            }
        }
        */
    }
}



void Gate::findEmptyTangentSphereWithMinAngleForEllipticEdge(rg_dList<BallGenerator*>& candidateBalls)
{
    Sphere gateBall[3] = { m_gateBallCCW1->getBall(), m_gateBallCCW2->getBall(), m_gateBallCCW3->getBall() }; 

    BallGenerator* currBall = rg_NULL;
    candidateBalls.reset4Loop();
	while( candidateBalls.setNext4Loop() )	{
	    currBall = candidateBalls.getEntity();

        Sphere  emptyTangentSphere[2];
        rg_INT  numTS = computeSphereTangentTo4SpheresOutside( gateBall[0], gateBall[1], gateBall[2], 
                                                               currBall->getBall(), emptyTangentSphere);


        for ( rg_INT i=0; i<numTS; i++ )  {
            
            if ( currBall == m_startBall && emptyTangentSphere[i].getCenter() == m_startVertex->getPoint() )
                continue;

            rg_Point3D vecAxisToCurrEnd = emptyTangentSphere[i].getCenter() - m_axisPoint;
            rg_REAL    angularDistance  = angleCCWOfTwoUnitVectors( m_normalOfEdgePlane, 
                                                                    m_vecAxisToStart, vecAxisToCurrEnd.getUnitVector() );

            if ( m_endBallList.getSize() > 0 )  {            
                // v_i:    emptyTangentSphere[i].getCenter()
                // v_star: current candidate end vertex
                rg_REAL distanceBetViAndVstar = m_emptySphereAtEndVertex.getCenter().distance( emptyTangentSphere[i].getCenter() );
            
                if ( rg_ZERO(distanceBetViAndVstar, m_tolerance) )  {
                    m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                    if ( emptyTangentSphere[i].getRadius() > m_emptySphereAtEndVertex.getRadius())
                        m_emptySphereAtEndVertex = emptyTangentSphere[i];
                }
                else  { 
                    if ( angularDistance < m_minAngularDistance ) {
                        m_minAngularDistance = angularDistance;
                        m_endBallList.removeAll();

                        m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                        m_emptySphereAtEndVertex = emptyTangentSphere[i];
                        m_numVertexBy4Balls      = numTS;
                    }
                }
            }
            else  {
                m_minAngularDistance = angularDistance;

                m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );
                m_emptySphereAtEndVertex = emptyTangentSphere[i];
                m_numVertexBy4Balls      = numTS;
            }

            /*
            if ( rg_EQ(angularDistance, m_minAngularDistance, m_tolerance) )  {
                //m_endBallList.add( currBall );
                m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );

                if ( emptyTangentSphere[i].getRadius() > m_emptySphereAtEndVertex.getRadius())
                    m_emptySphereAtEndVertex = emptyTangentSphere[i];
            }
            else if ( rg_LT(angularDistance, m_minAngularDistance, m_tolerance) ) {
                m_minAngularDistance = angularDistance;
                m_endBallList.removeAll();
                //m_endBallList.add( currBall );
                m_endBallList.add( pair<BallGenerator*, Sphere>(currBall, emptyTangentSphere[i]) );

                m_emptySphereAtEndVertex = emptyTangentSphere[i];
                m_numVertexBy4Balls      = numTS;
            }
            else {
            } 
            */
        }
    }
}


    
 bool Gate::doesDefineEllipticCurve(const Sphere& ball1, const Sphere& ball2, const Sphere& ball3)
{
    Sphere      gateBall[3]   = { ball1, ball2, ball3 }; 
    rg_Point3D  centerOfGB[3] = { gateBall[0].getCenter(), gateBall[1].getCenter(), gateBall[2].getCenter() };
    rg_REAL     radiusOfGB[3] = { gateBall[0].getRadius(), gateBall[1].getRadius(), gateBall[2].getRadius() };



    rg_Point3D normalOfCenterPlane;
	rg_Point3D vec21 = centerOfGB[1] - centerOfGB[0];
	rg_Point3D vec31 = centerOfGB[2] - centerOfGB[0];

	normalOfCenterPlane = vec21.crossProduct(vec31);
    normalOfCenterPlane.normalize();

    if ( rg_ZERO(normalOfCenterPlane.magnitude() ) ) {
        return true;
    }


	rg_TMatrix3D trMatrix;
	rg_TMatrix3D invMatrix;

	trMatrix.translate(-centerOfGB[0]);

    rg_FLAG bReverse = rg_EQ( normalOfCenterPlane.innerProduct( rg_Point3D(0.0, 0.0, 1.0) ), -1);
    if( !bReverse )
	    trMatrix.rotate(normalOfCenterPlane, rg_Point3D(0.0, 0.0, 1.0));
    else
        trMatrix.rotateY(rg_PI);

    if( !bReverse )
	    invMatrix.rotate(rg_Point3D(0.0, 0.0, 1.0), normalOfCenterPlane);
    else
        invMatrix.rotateY(-rg_PI);
	invMatrix.translate( centerOfGB[0] );


    rg_Point3D  trGateCenter[3] = { trMatrix*centerOfGB[0], trMatrix*centerOfGB[1], trMatrix*centerOfGB[2] };
	rg_Circle2D transformedCircle1( trGateCenter[0].getX(), trGateCenter[0].getY(), radiusOfGB[0] );
	rg_Circle2D transformedCircle2( trGateCenter[1].getX(), trGateCenter[1].getY(), radiusOfGB[1] );
	rg_Circle2D transformedCircle3( trGateCenter[2].getX(), trGateCenter[2].getY(), radiusOfGB[2] );

	rg_INT numOfCC = 0;
	rg_Circle2D tangentCircle[2];
    numOfCC = computeCircleTangentTo3CirclesOutside_GLOBAL( transformedCircle1, transformedCircle2, 
                                                            transformedCircle3, tangentCircle);   


    if ( numOfCC == 2 ) {
        return true;
        return false;
    }
    else if ( numOfCC == 1 ) {
        return false;
    }
    else {
        return false;
    }
}
