#include "VDFace.h"

#include "VDCell.h"
#include "VDLoop.h"
#include "VDPartialEdge.h"
#include "VDEdge.h"
#include "VDVertex.h"

#include "rg_BallGenerator.h"

#include "rg_RelativeOp.h"
#include "BetaEdge.h"
using namespace V::GeometryTier;


#include <float.h>
///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
VDFace::VDFace()
: m_rightCell(rg_NULL), m_leftCell(rg_NULL), m_faceEquation(rg_NULL), m_faceMesh(rg_NULL),
  m_isOnInfinity( rg_FALSE ),
  m_ProbeTangibility( rg_UNKNOWN), m_boundedness(rg_TRUE), m_triangularMesh(rg_NULL)
{
    m_visited  = rg_FALSE;
    m_betaEdge = rg_NULL;
}

VDFace::VDFace(const rg_INT& ID)
: TopologicalEntity(ID), m_rightCell(rg_NULL), m_leftCell(rg_NULL), m_faceEquation(rg_NULL), m_faceMesh(rg_NULL),
  m_isOnInfinity( rg_FALSE ),
  m_ProbeTangibility( rg_UNKNOWN), m_boundedness(rg_TRUE), m_triangularMesh(rg_NULL)
{
    m_visited  = rg_FALSE;
    m_betaEdge = rg_NULL;
}

VDFace::VDFace(const rg_INT& ID, VDCell* rightCell, VDCell* leftCell)
: TopologicalEntity(ID), m_rightCell(rightCell), m_leftCell(leftCell), m_faceEquation(rg_NULL), m_faceMesh(rg_NULL),
  m_isOnInfinity( rg_FALSE ),
  m_ProbeTangibility( rg_UNKNOWN), m_boundedness(rg_TRUE), m_triangularMesh(rg_NULL)
{
    m_visited  = rg_FALSE;
    m_betaEdge = rg_NULL;
}

VDFace::VDFace(const TopologicalEntity& aTopoEntity, VDCell* rightCell, VDCell* leftCell)
: TopologicalEntity(aTopoEntity), m_rightCell(rightCell), m_leftCell(leftCell), m_faceEquation(rg_NULL), m_faceMesh(rg_NULL),
  m_isOnInfinity( rg_FALSE ),
  m_ProbeTangibility( rg_UNKNOWN), m_boundedness(rg_TRUE), m_triangularMesh(rg_NULL)
{
    m_visited  = rg_FALSE;
    m_betaEdge = rg_NULL;
}

VDFace::VDFace(const VDFace& aFace)
: TopologicalEntity(aFace), m_rightCell(aFace.m_rightCell), m_leftCell(aFace.m_leftCell),
  m_isOnInfinity( aFace.m_isOnInfinity ),
  m_ProbeTangibility( aFace.m_ProbeTangibility ), m_boundedness(aFace.m_boundedness)
{
	m_loops.duplicateList(aFace.m_loops);

    //  faceEquation의 복사가 이루어져야 한다고 봄.
    m_faceEquation = aFace.m_faceEquation;

    if ( aFace.m_faceMesh != rg_NULL )
    {
        m_faceMesh = new rg_dList<rg_Point3D>;
        m_faceMesh->duplicateList(*aFace.m_faceMesh);
    }
    else
    {
        m_faceMesh = rg_NULL;
    }

    if ( aFace.m_triangularMesh != rg_NULL )
    {
        m_triangularMesh = new rg_dList<TriangleOnVDFace>;
        m_triangularMesh->duplicateList(*aFace.m_triangularMesh);
    }
    else
    {
        m_triangularMesh = rg_NULL;
    }

    m_visited  = aFace.m_visited;
    m_betaEdge = aFace.m_betaEdge;
}

VDFace::~VDFace()
{
    if ( m_faceMesh != rg_NULL )
        delete m_faceMesh;

    if ( m_triangularMesh != rg_NULL )
        delete m_triangularMesh;


    //if ( m_betaEdge != rg_NULL ) {
    //    m_betaEdge->disconnectVFace(this);
    //}
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
VDCell* VDFace::getRightCell() const
{
    return m_rightCell;
}

VDCell* VDFace::getLeftCell() const
{
    return m_leftCell;
}


VDLoop* VDFace::getOuterLoop() const
{
    if ( m_loops.getSize() > 0 && m_loops.getAt(0)->isOuterLoop() )
    {
        return m_loops.getAt(0);
    }

    VDLoop* currLoop = rg_NULL;

	m_loops.reset4Loop();
	while( m_loops.setNext4Loop() )
	{
	    currLoop = m_loops.getEntity();

        if ( currLoop->isOuterLoop() )  
            return currLoop;
    }

    return rg_NULL;
    
}

VDLoop* VDFace::getLoop(const rg_INT& i) const
{
    if ( i < m_loops.getSize() ) 
        return m_loops.getAt(i);
    else  
        return rg_NULL;
}

rg_INT  VDFace::getNumOfLoops() const
{
    return m_loops.getSize();
}


rg_dList<VDLoop*>* VDFace::getLoops()
{
    return &m_loops;
}


rg_Surface* VDFace::getFaceEquation() const
{
    return m_faceEquation;
}

rg_dList<rg_Point3D>* VDFace::getFaceMesh()
{
    return m_faceMesh;
}

rg_dList<TriangleOnVDFace>* VDFace::getTriangularMesh()
{
    return m_triangularMesh;
}

rg_FLAG VDFace::isTangible() const
{
    return m_ProbeTangibility;
}


    
rg_BOOL VDFace::isVirtual() const
{
    if ( m_leftCell->getGenerator() == rg_NULL ) {
        return rg_TRUE;
    }
    else if ( m_rightCell->getGenerator() == rg_NULL ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_BOOL VDFace::isUnbounded() const
{
    rg_BOOL isUnboundedFace = rg_FALSE;

    rg_dList<VDEdge*> edgeList;
    getOuterLoop()->searchBoundingEdges( edgeList );

    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        VDEdge* currEdge = edgeList.getEntity();

        if ( currEdge->isVirtual() ) {
            isUnboundedFace = rg_TRUE;
            break;
        }
    }

    return isUnboundedFace;
}



rg_FLAG VDFace::isOnInfinity() const
{
    return m_isOnInfinity;
}

void VDFace::setInfinity()
{
	m_isOnInfinity = rg_TRUE;
}

void    VDFace::isTangible(const rg_FLAG& tangibility)
{
    m_ProbeTangibility = tangibility;
}

rg_FLAG VDFace::isBounded() const
{
    return m_boundedness;
}

void VDFace::isBounded( const rg_FLAG& boundedness)
{
    m_boundedness = boundedness;
}




rg_FLAG VDFace::isThereInnerLoop() 
{
    VDLoop* currLoop = rg_NULL;
    m_loops.reset4Loop();
    while ( m_loops.setNext4Loop() )  {
        currLoop = m_loops.getEntity();

        if ( !currLoop->isOuterLoop() )
            return rg_TRUE;
    }

    return rg_FALSE;
}



///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void VDFace::setRightCell(VDCell* rightCell)
{
    m_rightCell = rightCell;
}

void VDFace::setLeftCell(VDCell* leftCell)
{
    m_leftCell = leftCell;
}


void VDFace::connectFaceWithRightCell(VDCell* rightCell)
{
    m_rightCell = rightCell;
    rightCell->addBoundingFace( this );
}

void VDFace::connectFaceWithLeftCell(VDCell* leftCell)
{
    m_leftCell = leftCell;
    leftCell->addBoundingFace( this );
}

void VDFace::connectFaceWithRightAndLeftCell(VDCell* rightCell, VDCell* leftCell)
{
    m_rightCell = rightCell;
    m_leftCell  = leftCell;

    rightCell->addBoundingFace( this );
    leftCell->addBoundingFace( this );
}



void VDFace::setLoops(const rg_dList<VDLoop*>& loops)
{
	m_loops.duplicateList(loops);
}

void VDFace::addLoop(VDLoop* loop)
{
    m_loops.addTail( loop );
}



void VDFace::setFaceEquation(rg_Surface* faceEquation)
{
    m_faceEquation = faceEquation;
}


///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
VDFace& VDFace::operator=(const VDFace& aFace)
{
    if ( this == &aFace )
        return *this;

    m_ID        = aFace.m_ID;
    m_rightCell = aFace.m_rightCell;
    m_leftCell  = aFace.m_leftCell;

	m_loops.duplicateList(aFace.m_loops);

    m_faceEquation = aFace.m_faceEquation;

    if ( aFace.m_faceMesh != rg_NULL )
    {
        m_faceMesh = new rg_dList<rg_Point3D>;
        m_faceMesh->duplicateList(*aFace.m_faceMesh);
    }
    else
    {
        m_faceMesh = rg_NULL;
    }


    if ( aFace.m_triangularMesh != rg_NULL )
    {
        m_triangularMesh = new rg_dList<TriangleOnVDFace>;
        m_triangularMesh->duplicateList(*aFace.m_triangularMesh);
    }
    else
    {
        m_triangularMesh = rg_NULL;
    }

    m_isOnInfinity     = aFace.m_isOnInfinity;

    m_ProbeTangibility = aFace.m_ProbeTangibility;
    m_boundedness      = aFace.m_boundedness;


    m_visited  = aFace.m_visited;
    m_betaEdge = aFace.m_betaEdge;

    return *this;
}



/*
rg_REAL VDFace::computeArea()
{
	if( isOnInfinity() )
	    return rg_MAX_REAL;

	if(m_faceMesh == NULL)
	{
		makeMeshForVoronoiFace();

		//unbounded face case
		if(m_faceMesh == NULL)
			return rg_MAX_REAL;
	}

	rg_Point3D vec1, vec2;
	rg_Point3D pt1, pt2, pt3;

	rg_REAL totalArea = 0.;
	rg_REAL triArea = 0.;
	m_faceMesh->reset4Loop();
	while(m_faceMesh->setNext4Loop())
	{
		//first threes are normals
		m_faceMesh->setNext4Loop();
		m_faceMesh->setNext4Loop();
		m_faceMesh->setNext4Loop();

		pt1 = m_faceMesh->getEntity();
		m_faceMesh->setNext4Loop();
		pt2 = m_faceMesh->getEntity();
		m_faceMesh->setNext4Loop();
		pt3 = m_faceMesh->getEntity();
		
		vec1 = pt2 - pt1;
		vec2 = pt3 - pt2;
		triArea =  vec1.crossProduct(vec2).magnitude(); 
		
		totalArea += triArea;
	}

    return totalArea/2.;
}

///////////////////////////////////////////////////////////////////////////////
//
//  topological operators..


///////////////////////////////////////////////////////////////////////////////
//
//  geometric operators..
rg_FLAG VDFace::isThereIntersectionWithSectionalYZPlane(const rg_REAL& xOfPlane)
{
    VDLoop* currLoop = rg_NULL;

	m_loops.reset4Loop();
	while(m_loops.setNext4Loop())
	{
		currLoop = m_loops.getEntity();

        rg_FLAG        existIntersection = rg_FALSE;
        VDPartialEdge* currPartEdge      = currLoop->getPartialEdge();
        do 
        {
            VDEdge* currEdge = currPartEdge->getOriginalEdge();

            if ( !currEdge->isBoundedEdge() )
            {
                currPartEdge = currPartEdge->getNextPartialEdgeInLoop();
                continue;
            }

            rg_REAL largeX = currEdge->getEndVertex()->getPoint().getX();
            rg_REAL smallX = currEdge->getStartVertex()->getPoint().getX();
            if ( rg_GT( smallX, largeX ) )
            {
                rg_REAL tempX( smallX );

                smallX = largeX;
                largeX = tempX;
            }

            if ( rg_GE(xOfPlane, smallX) && rg_LE(xOfPlane, largeX) )
            {
                existIntersection = rg_TRUE;
                break;
            }

            currPartEdge = currPartEdge->getNextPartialEdgeInLoop();

        } while ( currPartEdge != currLoop->getPartialEdge() );


        if ( existIntersection )
        {
            return rg_TRUE;
        }
    }

    return rg_FALSE;

}

void VDFace::makeMeshForVoronoiFace()
{
    VDLoop* outerLoop = getOuterLoop();
    if ( outerLoop == rg_NULL )
        return;
    
    if ( outerLoop->isClosedLoop() )
    {
        Sphere rightSphere = ((BallGeneratorCore*)m_rightCell->getGenerator())->getBall();
        Sphere leftSphere  = ((BallGeneratorCore*)m_leftCell->getGenerator())->getBall();

        rg_dList< OrientedCurve > curvesForVoronoiEdge;
        
        OrientedCurve currCurve;

        VDPartialEdge* startPartEdge = outerLoop->getPartialEdge();
        VDPartialEdge* currPartEdge  = startPartEdge;

		//disconnected case
		if(currPartEdge->getOriginalEdge()->getStartVertex() == NULL)
		{
            currCurve.m_first  = (rg_RBzCurve3D*)currPartEdge->getOriginalEdge()->getEdgeEquation();
            currCurve.m_second = -1;

            curvesForVoronoiEdge.addTail( currCurve );
		}
		else
		{
            int count = 0;
            do
            {
                rg_RBzCurve3D* curve = (rg_RBzCurve3D*)currPartEdge->getOriginalEdge()->getEdgeEquation();

                rg_Point3D vMid=( curve->getCtrlPt(0) + curve->getCtrlPt(2) ) *0.5;
                rg_Point3D vGap=vMid-curve->getCtrlPt(1);


                currCurve.m_first  = curve;
                currCurve.m_second = currPartEdge->isRightOrientationInLoop();

                curvesForVoronoiEdge.addTail( currCurve );

                currPartEdge = currPartEdge->getNextPartialEdgeInLoop();

                count++;
            } while (currPartEdge != startPartEdge);
        }

        m_faceMesh = triangulateVoronoiFace(rightSphere, leftSphere, curvesForVoronoiEdge);
    }
}

void VDFace::makeTriangularMeshForVoronoiFace()
{
    VDLoop* outerLoop = getOuterLoop();

    if ( outerLoop->isClosedLoop() )
    {
        Sphere rightSphere = ((BallGeneratorCore*)m_rightCell->getGenerator())->getBall();
        Sphere leftSphere  = ((BallGeneratorCore*)m_leftCell->getGenerator())->getBall();

        rg_dList< OrientedCurve > curvesForVoronoiEdge;
        
        OrientedCurve currCurve;

        VDPartialEdge* startPartEdge = outerLoop->getPartialEdge();
        VDPartialEdge* currPartEdge  = startPartEdge;

        int count = 0;
        do
        {
            rg_RBzCurve3D* curve = (rg_RBzCurve3D*)currPartEdge->getOriginalEdge()->getEdgeEquation();

            rg_Point3D vMid=( curve->getCtrlPt(0) + curve->getCtrlPt(2) ) *0.5;
            rg_Point3D vGap=vMid-curve->getCtrlPt(1);


            currCurve.m_first  = curve;
            currCurve.m_second = currPartEdge->isRightOrientationInLoop();

            curvesForVoronoiEdge.addTail( currCurve );

            currPartEdge = currPartEdge->getNextPartialEdgeInLoop();

            count++;
        } while (currPartEdge != startPartEdge);


        rg_dList<rg_Point3D>* triangles = triangulateVoronoiFace(rightSphere, leftSphere, curvesForVoronoiEdge, 5, 5);

        if ( triangles == rg_NULL )
        {
            m_triangularMesh = rg_NULL;
            return;
        }

        rg_Point3D rightBallCenter = ((BallGeneratorCore*)m_rightCell->getGenerator())->getCenter();
        rg_REAL    rightBallRadius = ((BallGeneratorCore*)m_rightCell->getGenerator())->getRadius();

        m_triangularMesh = new rg_dList<TriangleOnVDFace>;

        rg_Point3D VerticesAndNormalsTriangle[6];
        rg_INT     i = 0;

        triangles->reset4Loop();
	    while( triangles->setNext4Loop() )  
        {
            VerticesAndNormalsTriangle[i] = triangles->getEntity();
            i++;

            if ( i == 6 )
            {
                i=0;

                rg_REAL distanceToBall = 0.0;
                TriangleOnVDFace triangleOnVDF;
                for ( rg_INT j=0; j<3; j++)
                {
                    triangleOnVDF.setNormalOnVertex(j, VerticesAndNormalsTriangle[j]);
                    triangleOnVDF.setVertex(j, VerticesAndNormalsTriangle[3+j]);

                    rg_REAL distance = rightBallCenter.distance(VerticesAndNormalsTriangle[3+j]) - rightBallRadius;
                    triangleOnVDF.setDistance(j, distance);
                    distanceToBall += distance;
                }

                triangleOnVDF.setDistanceToGenerator( distanceToBall/3.0 );
                m_triangularMesh->addTail( triangleOnVDF );
            }
        }
    }
}

rg_REAL VDFace::getMaxDistanceBetweenTrianglesAndGenerator()
{
    rg_REAL maxDistance = -DBL_MAX;
    TriangleOnVDFace* triangle = rg_NULL;
    m_triangularMesh->reset4Loop();
	while( m_triangularMesh->setNext4Loop() )  
    {
        triangle = m_triangularMesh->getpEntity();

        rg_REAL distance = triangle->getDistanceToGenerator();

        if ( rg_GT(distance, maxDistance) )
            maxDistance = distance;
    }

    return maxDistance;
}



rg_dList<rg_Point3D>* VDFace::triangulateVoronoiFace(const Sphere& sphere1, const Sphere& sphere2, CrvList& bzCrvList, 
                                                    rg_INT numPointsOnCurve, 
                                                    rg_INT numGridPoints)
{
	Sphere s1, s2;
	//assume that s1.r < s2.r
	if( sphere1.getRadius() < sphere2.getRadius() )
	{
		s1 = sphere1;
		s2 = sphere2;
	}
	else
	{
		s1 = sphere2;
		s2 = sphere1;
	}

	rg_Point3D center  = (s1.getCenter()+s2.getCenter())/2;
	Vector3D   baseVec = (s1.getCenter() - s2.getCenter()).getUnitVector();
	

	//translate & rotate!
	rg_TMatrix3D trMatrix;
	trMatrix.translate(-center);
    const rg_FLAG bReverse=rg_EQ(baseVec.innerProduct(Vector3D(0,0,1)), -1);
    if( !bReverse )
	    trMatrix.rotate(baseVec,Vector3D(0,0,1));
    else
        trMatrix.rotateY(rg_PI);

	rg_TMatrix3D bkTrMatrix;
    if( !bReverse )
	    bkTrMatrix.rotate(Vector3D(0,0,1),baseVec);
    else
        bkTrMatrix.rotateY(-rg_PI);
//        bkTrMatrix.reflectToPoint(Vector3D(0,0,0));
	bkTrMatrix.translate(center);

	Pt3DVec polygon3;

	rg_REAL t = 0.0;
	rg_INT index = 0;
    OrientedCurve currCurve;
	rg_RBzCurve3D tmpCrv;
	rg_Point3D ctrlPt[3];

//    int curveID = 0;

	bzCrvList.reset4Loop();
	while(bzCrvList.setNext4Loop())
	{
        currCurve = bzCrvList.getEntity();
        tmpCrv    = *(currCurve.m_first);

		ctrlPt[0] = tmpCrv.getCtrlPt(0);
		ctrlPt[1] = tmpCrv.getCtrlPt(1);
		ctrlPt[2] = tmpCrv.getCtrlPt(2);

        if ( currCurve.m_second )
        {
		    tmpCrv.setCtrlPt(0, trMatrix*ctrlPt[0]);
		    tmpCrv.setCtrlPt(1, trMatrix*ctrlPt[1]);
		    tmpCrv.setCtrlPt(2, trMatrix*ctrlPt[2]);
        }
        else
        {
		    tmpCrv.setCtrlPt(0, trMatrix*ctrlPt[2]);
		    tmpCrv.setCtrlPt(1, trMatrix*ctrlPt[1]);
		    tmpCrv.setCtrlPt(2, trMatrix*ctrlPt[0]);
        }

		const rg_REAL minLengthForPtsOnCurve = 0.5;
		rg_REAL stepLength = tmpCrv.evaluatePt(0).distance(tmpCrv.evaluatePt(1./numPointsOnCurve));
		if( stepLength < minLengthForPtsOnCurve )
		{
			int numModifiedPtsOnCrv = ceil(numPointsOnCurve*stepLength/minLengthForPtsOnCurve); 
			for(int i=0; i<numModifiedPtsOnCrv; i++)
			{
				t = (rg_REAL)i/numModifiedPtsOnCrv;
				rg_Point3D ptOnCurve = tmpCrv.evaluatePt(t);
				polygon3.push_back( PointIndexPair(getPoint3d( ptOnCurve ), ++index) );
			}
		}
		else
		{
			for(int i=0; i<numPointsOnCurve; i++)
			{
				t = (rg_REAL)i/numPointsOnCurve;
				rg_Point3D ptOnCurve = tmpCrv.evaluatePt(t);
				polygon3.push_back( PointIndexPair(getPoint3d( ptOnCurve ), ++index) );
			}
		}
	}

	//input data code

	Pt3DVec::iterator iter;
    Polygon_CGAL polygon2;
	rg_REAL maxMagnitude = -1.;
//	rg_REAL minMagnitude = rg_MAX_REAL;
	rg_REAL x = 0;
	rg_REAL y = 0;
	rg_REAL magnitude = 0;

	for(iter = polygon3.begin(); iter != polygon3.end(); iter++)
	{
		x = iter->first.x();
		y = iter->first.y();
		polygon2.push_back( Point2d(x, y) );
		magnitude = x*x + y*y;
		if(magnitude > maxMagnitude)
			maxMagnitude = magnitude;
//		if(magnitude < minMagnitude)
//			minMagnitude = magnitude;
	}
	maxMagnitude = sqrt(maxMagnitude);

	Triangulation tr;

	rg_REAL z = 0.0;
	rg_REAL f = (center-s1.getCenter()).magnitude();
	rg_REAL r = s2.getRadius() - s1.getRadius();

	int maxLevel = 5; 
	rg_REAL theta = rg_PI/3.;
	rg_REAL cosine = 0;
	rg_REAL sine = 0;
	z = zValue(0,0,f,r);
	tr.insert( Point3d(0,0,z) );
	for(int k = 1; k <= maxLevel; k++)
	{
		for(int j = 0; j< 6*k; j++)
		{
			cosine = cos(theta/k*j);
			sine = sin(theta/k*j);

			x = cosine * maxMagnitude*1.1/maxLevel*k;
			y = sine   * maxMagnitude*1.1/maxLevel*k;

			z = zValue(x,y,f,r);
			tr.insert( Point3d(x,y,z) );
		}
	}


	int numOfPolyVtx = polygon3.size();
	Vertex_handle* vHandle = new Vertex_handle [numOfPolyVtx];
	for(int j = 0; j < numOfPolyVtx; j++)
	{
		//tr.insert(polygon3[j].first, polygon3[j+1].first); //constrained edge
		vHandle[j] = tr.insert(polygon3[j].first);
	}

	for(j=0; j < numOfPolyVtx-1; j++)
	{
		tr.insert(vHandle[j], vHandle[j+1]);
	}
	tr.insert(vHandle[numOfPolyVtx-1],vHandle[0]);
	delete [] vHandle;

	//normal computation
	Point3d tPt;
	Vertex_iterator vit = tr.vertices_begin();
	while( vit != tr.vertices_end() )
	{
		tPt = Point3d(vit->point().x(),
					  vit->point().y(),
					  zValue(vit->point().x(),vit->point().y(),f,r)
					 );
		vit->normal = getNormal(tPt,f,r);
		vit++;
	}

	rg_TMatrix3D rotMat;
    if( !bReverse )
	    rotMat.rotate(Vector3D(0,0,1),baseVec);
    else
        rotMat.rotateY(-rg_PI);

	rg_dList<rg_Point3D>* triangleList = new rg_dList<rg_Point3D>;


	rg_Point3D pt0, pt1, pt2;
	Vector3D vec1, vec2;
	Face_iterator it = tr.faces_begin(), beyond = tr.faces_end();
	while( it != beyond )
	{
		//Determine if triangle is inside or outside of the polygon.
		if( polygon2.bounded_side( centerOfMass(it) ) == CGAL::ON_BOUNDED_SIDE )
		{
			pt0 = getRGPoint3D(it->vertex(0)->point());
			pt0.setZ( zValue(pt0.getX(), pt0.getY(), f, r) );
			pt0 = bkTrMatrix*pt0;

			pt1 = getRGPoint3D(it->vertex(1)->point());
			pt1.setZ( zValue(pt1.getX(), pt1.getY(), f, r) );
			pt1 = bkTrMatrix*pt1;

			pt2 = getRGPoint3D(it->vertex(2)->point());
			pt2.setZ( zValue(pt2.getX(), pt2.getY(), f, r) );
			pt2 = bkTrMatrix*pt2;

//			vec1 = pt1 - pt0;
//			vec2 = pt2 - pt1;

//			Vector3D normal = (vec1*vec2).getUnitVector();
			//triangleList->add(normal);

			triangleList->add(rotMat*getRGPoint3D(it->vertex(0)->normal));
			triangleList->add(rotMat*getRGPoint3D(it->vertex(1)->normal));
			triangleList->add(rotMat*getRGPoint3D(it->vertex(2)->normal));

			triangleList->add(pt0);
			triangleList->add(pt1);
			triangleList->add(pt2);
		}
		it++;
	}


	return triangleList;
}

													
Point2d VDFace::centerOfMass(Face_iterator& it)
{
	Vector2d sumv(0.0, 0.0);
	
	sumv = Vector2d( it->vertex(0)->point().x() + it->vertex(1)->point().x() + it->vertex(2)->point().x(),
		             it->vertex(0)->point().y() + it->vertex(1)->point().y() + it->vertex(2)->point().y() );

	return Point2d(0.,0.) + sumv/3.;
}

Point3d VDFace::getPoint3d(const rg_Point3D& pt)
{
	return Point3d(pt.getX(), pt.getY(), pt.getZ());
}

rg_Point3D VDFace::getRGPoint3D(const Point3d& pt)
{
//	return rg_Point3D(pt.x().to_double(), pt.y().to_double(), pt.z().to_double());
	return rg_Point3D(pt.x(), pt.y(), pt.z());
}

rg_REAL VDFace::zValue(const leda_real& leda_x, const leda_real& leda_y, const rg_REAL& f, const rg_REAL& r)
{
    rg_REAL x = leda_x.to_double(), y = leda_y.to_double();
	return(2.0/(-4.0*r*r+16.0*f*f)*std::sqrt((-r*r+4.0*f*f)*(4.0*x*x+4.0*y*y+4.0*f*f-r*r))*r);
}

rg_REAL VDFace::zValue(const rg_REAL& x, const rg_REAL& y, const rg_REAL& f, const rg_REAL& r)
{
    return(2.0/(-4.0*r*r+16.0*f*f)*std::sqrt((-r*r+4.0*f*f)*(4.0*x*x+4.0*y*y+4.0*f*f-r*r))*r);
}


//Point3d getNormal(const rg_REAL& x,const rg_REAL& y,const rg_REAL& z,const rg_REAL& f,const rg_REAL& r)
Point3d VDFace::getNormal(const Point3d& pt, const rg_REAL& f,const rg_REAL& r)
{

	double x = pt.x();
	double y = pt.y();
	double z = pt.z();

	double dx = 1/(sqrt(x*x+y*y+z*z-2.0*z*f+f*f))*x-1/(sqrt(x*x+y*y+z*z+2.0*z*f+f*f))*x;
	double dy = 1/(sqrt(x*x+y*y+z*z-2.0*z*f+f*f))*y-1/(sqrt(x*x+y*y+z*z+2.0*z*f+f*f))*y;
	double dz = 1/(sqrt(x*x+y*y+z*z-2.0*z*f+f*f))*(2.0*z-2.0*f)/2.0-1/
				(sqrt(x*x+y*y+z*z+2.0*z*f+f*f))*(2.0*z+2.0*f)/2.0;

	double magnitude = sqrt(dx*dx+dy*dy+dz*dz);
	
	return Point3d(dx/magnitude, dy/magnitude, dz/magnitude);
}
*/

///////////////////////////////////////////////////////////////////////////////
//
//  topological quires..
//

void  VDFace::inquireBoundingVertices(rg_dList<VDVertex*>& vertexList) const
{
    VDLoop* currLoop = rg_NULL;

    m_loops.reset4Loop();
    while ( m_loops.setNext4Loop() )  {
        currLoop = m_loops.getEntity();

        VDPartialEdge* startPrEdge = currLoop->getPartialEdge();
        VDPartialEdge* currPrEdge  = startPrEdge;

        do  {
            if ( currPrEdge->isRightOrientationInLoop() )
                vertexList.add( currPrEdge->getOriginalEdge()->getEndVertex() );
            else
                vertexList.add( currPrEdge->getOriginalEdge()->getStartVertex() );

            currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
        } while ( currPrEdge != startPrEdge );
    }
}



void  VDFace::inquireBoundingEdges(rg_dList<VDEdge*>& edgeList) const
{
    VDLoop* currLoop = rg_NULL;

    m_loops.reset4Loop();
    while ( m_loops.setNext4Loop() )  {
        currLoop = m_loops.getEntity();

        VDPartialEdge* startPrEdge = currLoop->getPartialEdge();
        VDPartialEdge* currPrEdge  = startPrEdge;
        do  {
            edgeList.add( currPrEdge->getOriginalEdge() );

            currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
        } while ( currPrEdge != startPrEdge );
    }
}



void  VDFace::inquireIncidentFaces(rg_dList<VDFace*>& faceList)
{
    VDLoop* currLoop = rg_NULL;

    m_loops.reset4Loop();
    while ( m_loops.setNext4Loop() )  {
        currLoop = m_loops.getEntity();

        VDPartialEdge* startPrEdge = currLoop->getPartialEdge();
        VDPartialEdge* currPrEdge  = startPrEdge;
        do  {
            VDPartialEdge* nextPrEdgeInRadial     = currPrEdge->getNextPartialEdgeInRadialCycle();
            VDPartialEdge* nextNextPrEdgeInRadial = nextPrEdgeInRadial->getNextPartialEdgeInRadialCycle();

            faceList.addWithoutSame( nextPrEdgeInRadial->getLoop()->getFace() );
            faceList.addWithoutSame( nextNextPrEdgeInRadial->getLoop()->getFace() );

            currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
        } while ( currPrEdge != startPrEdge );
    }
}



void  VDFace::inquireIncidentCells(rg_dList<VDCell*>& cellList)
{
    cellList.add( m_rightCell );
    cellList.add( m_leftCell );
}



rg_BOOL VDFace::isIncidentTo(   VDCell* currCell ) const
{
    if ( m_leftCell == currCell ) {
        return rg_TRUE;
    }
    else if ( m_rightCell == currCell ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



void VDFace::connectBetaEdge(BetaEdge* b_edge)
{
    m_betaEdge = b_edge;
}



void VDFace::disconnectBetaEdge(BetaEdge* b_edge)
{
    if (  m_betaEdge == b_edge ) {
        m_betaEdge = rg_NULL;
    }
}



