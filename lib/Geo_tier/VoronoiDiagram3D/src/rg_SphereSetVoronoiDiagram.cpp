#include "rg_SphereSetVoronoiDiagram.h"


#include "VIDIC.h"
#include "Gate.h"

#include "FunctionsForVoronoiDiagram3D.h"
#include "BucketCellIndex.h"

#include "rg_TMatrix3D.h"
#include "rg_RBzCurve3D.h"
#include "rg_RevolvedSurface.h"
#include "rg_RelativeOp.h"

#include "Triangulation3D.h"

#include "PathInvestigator.h"

#include "BetaVertex.h"
#include "BetaEdge.h"
#include "BetaFace.h"
#include "BetaCell.h"

#include "Plane.h"

#include <float.h>
#include <time.h>
#include <math.h>
#include <stdlib.h> 

#define _CRTDBG_MAP_ALLOC  
#include <iostream>
#include <random>   
#include <cstdlib>
#include <cstdio>  
//#include <crtdbg.h>

using namespace V::GeometryTier;



///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
rg_SphereSetVoronoiDiagram::rg_SphereSetVoronoiDiagram()
: m_minPointOfBoundingBoxForBalls(DBL_MAX, DBL_MAX, DBL_MAX), m_maxPointOfBoundingBoxForBalls(-DBL_MAX, -DBL_MAX, -DBL_MAX),
  m_minPointOfBoundingBoxForVD(DBL_MAX, DBL_MAX, DBL_MAX), m_maxPointOfBoundingBoxForVD(-DBL_MAX, -DBL_MAX, -DBL_MAX),
  m_minRadiusOfBalls(DBL_MAX), m_maxRadiusOfBalls(-DBL_MAX), m_AvgRadiusOfBalls(0.0), 
  m_hasProbeTangibility(rg_FALSE), m_probeRadius(0.0) 
{
    m_status = NORMAL_STATE;

    m_sizeOfBucketElement = 7.0;

	rg_INT i=0;
    for ( i=0; i<NUM_TIMES_FOR_ANALYSIS; i++)
        computingTime[i] = 0;


    timeToBeExtracted[0] = 0;
    timeToBeExtracted[1] = 0;



    /*
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //
    m_bStartCheck = rg_FALSE;
    for ( i=0; i<2; i++ )  {
        m_numCandidateBallsOfElliptic[i] = 0;
        for ( rg_INT j=0; j<2; j++ )
            for ( rg_INT k=0; k<3; k++ )
                m_numCandidateBallsOfNonElliptic[i][j][k] = 0;
    }

    for ( i=0; i<2; i++ )  {
        m_numIntersectingBallsWithFrontFR[i] = 0;

        for ( rg_INT j=0; j<3; j++ )  {
            m_avgRadiusOfSphereForRf[i][j] = 0.0;
            
            m_numEdges[i][j] = 0;
            m_numEndBallsInPosition[i][j] = 0;
        }
    }

    for ( i=0; i<2; i++ )  {
        m_numBallsForUnboundedEdge[i] = 0;
        for ( rg_INT j=0; j<4; j++ )  {
            m_numSearchedBucketElements[i][j] = 0;
        }
    }

    m_frequencyOfCandidateBall = rg_NULL;
    //
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //*/

    ON_DEBUG = rg_TRUE;
    ON_DEBUG = rg_FALSE;

    if ( ON_DEBUG )  {
        foutEVDS.open("tracingEVDS.log");
        foutEVDS.precision(18);
    }
}



rg_SphereSetVoronoiDiagram::rg_SphereSetVoronoiDiagram(const rg_SphereSetVoronoiDiagram& VD)
{
    duplicate(VD);
}



rg_SphereSetVoronoiDiagram::~rg_SphereSetVoronoiDiagram()
{
    /*
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //
    if ( m_frequencyOfCandidateBall != rg_NULL )
        delete [] m_frequencyOfCandidateBall;
    //
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //*/
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
rg_dList< BallGenerator >* rg_SphereSetVoronoiDiagram::getGlobalGeneratorList() 
{
    return &m_GlobalGeneratorList;
}


rg_INT rg_SphereSetVoronoiDiagram::getNumOfBalls() const
{
    return m_GlobalGeneratorList.getSize();
}

rg_Point3D rg_SphereSetVoronoiDiagram::getMinPointOfBoundingBoxForBalls() const
{
    return m_minPointOfBoundingBoxForBalls;
}

rg_Point3D rg_SphereSetVoronoiDiagram::getMaxPointOfBoundingBoxForBalls() const
{
    return m_maxPointOfBoundingBoxForBalls;
}

rg_Point3D rg_SphereSetVoronoiDiagram::getMinPointOfBoundingBoxForVoronoiDiagram() const
{
    return m_minPointOfBoundingBoxForVD;
}

rg_Point3D rg_SphereSetVoronoiDiagram::getMaxPointOfBoundingBoxForVoronoiDiagram() const
{
    return m_maxPointOfBoundingBoxForVD;
}

rg_REAL rg_SphereSetVoronoiDiagram::getAverageRadiusOfBalls() const
{
    return m_AvgRadiusOfBalls;
}

rg_FLAG rg_SphereSetVoronoiDiagram::isProbeTangibilityCalculated() const
{
    return m_hasProbeTangibility;
}

rg_REAL rg_SphereSetVoronoiDiagram::getProbeRadius() const
{
    return m_probeRadius;
}

///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void rg_SphereSetVoronoiDiagram::setAvgRadiusOfBalls(const rg_REAL& avgRadius)
{
    m_AvgRadiusOfBalls = avgRadius;
}

void rg_SphereSetVoronoiDiagram::setMinPointOfBoundingBoxForBalls(const rg_Point3D& minPt)
{
    m_minPointOfBoundingBoxForBalls = minPt;
}

void rg_SphereSetVoronoiDiagram::setMaxPointOfBoundingBoxForBalls(const rg_Point3D& maxPt)
{
    m_maxPointOfBoundingBoxForBalls = maxPt;
}

void   rg_SphereSetVoronoiDiagram::setBucket(const rg_REAL& loadFactor)
{
	//m_bucket.constructBucket( loadFactor,
    //                          m_minPointOfBoundingBoxForBalls,
    //                          m_maxPointOfBoundingBoxForBalls,
    //                          m_GlobalGeneratorList );

    m_bucket.constructBucketBySizeOfElement( loadFactor,
                                             m_minPointOfBoundingBoxForBalls,
                                             m_maxPointOfBoundingBoxForBalls,
                                             m_GlobalGeneratorList );

}



BallGenerator* rg_SphereSetVoronoiDiagram::addBallGenerator(const BallGenerator& aBall)
{
    Sphere  ball = aBall.getBall();
    rg_FLAG isContained = rg_FALSE;

    BallGenerator* currBall = rg_NULL;
	m_GlobalGeneratorList.reset4Loop();
	while( m_GlobalGeneratorList.setNext4Loop() )  
    {
		currBall = m_GlobalGeneratorList.getpEntity();

        if ( currBall->isContainedIn( ball ) )
        {
            isContained = rg_TRUE;
            break;
        }
    }

    if ( isContained == rg_FALSE )
    {
        BallGenerator* ptrBall = m_GlobalGeneratorList.addTail( aBall );

        VDCell* ptrCell = m_GlobalCellList.addTail( VDCell( m_GlobalCellList.getSize() ) );

        connectCellAndGenerator( ptrCell, ptrBall );

        return ptrBall;
    }
    else
    {
        return rg_NULL;
    }
}



rg_BOOL rg_SphereSetVoronoiDiagram::addSphere(const Sphere& sphere, void* property, const rg_INT& IDFromInput)
{
    rg_REAL x      = sphere.getCenter().getX();
    rg_REAL y      = sphere.getCenter().getY();
    rg_REAL z      = sphere.getCenter().getZ();
    rg_REAL radius = sphere.getRadius();


	if( rg_GT( x + radius, m_maxPointOfBoundingBoxForBalls.getX() ) )
		m_maxPointOfBoundingBoxForBalls.setX( x + radius );
	if( rg_GT( y + radius, m_maxPointOfBoundingBoxForBalls.getY() ) )
		m_maxPointOfBoundingBoxForBalls.setY( y + radius );
	if( rg_GT( z + radius, m_maxPointOfBoundingBoxForBalls.getZ() ) )
		m_maxPointOfBoundingBoxForBalls.setZ( z + radius );

	if( rg_LT( x - radius, m_minPointOfBoundingBoxForBalls.getX() ) )
		m_minPointOfBoundingBoxForBalls.setX( x - radius );
	if( rg_LT( y - radius, m_minPointOfBoundingBoxForBalls.getY() ) )
		m_minPointOfBoundingBoxForBalls.setY( y - radius );
	if( rg_LT( z - radius, m_minPointOfBoundingBoxForBalls.getZ() ) )
		m_minPointOfBoundingBoxForBalls.setZ( z - radius );


    if ( rg_GT(radius, m_maxRadiusOfBalls ) )
        m_maxRadiusOfBalls = radius;
	if ( rg_LT(radius, m_minRadiusOfBalls ) )
        m_minRadiusOfBalls = radius;

    rg_REAL sumOfRadius = m_AvgRadiusOfBalls*m_GlobalGeneratorList.getSize();

    BallGenerator ball( m_GlobalGeneratorList.getSize(), sphere, property, IDFromInput);    
    
    BallGenerator* ptrBall = addBallGenerator( ball );


    if ( ptrBall != rg_NULL ) {
        m_AvgRadiusOfBalls = (sumOfRadius+radius)/m_GlobalGeneratorList.getSize();
        
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_BOOL rg_SphereSetVoronoiDiagram::addSphereWithoutItsValidityCheck(const Sphere& sphere, void* property, const rg_INT& IDFromInput)
{
    rg_REAL x = sphere.getCenter().getX();
    rg_REAL y = sphere.getCenter().getY();
    rg_REAL z = sphere.getCenter().getZ();
    rg_REAL radius = sphere.getRadius();


    if (rg_GT(x + radius, m_maxPointOfBoundingBoxForBalls.getX()))
        m_maxPointOfBoundingBoxForBalls.setX(x + radius);
    if (rg_GT(y + radius, m_maxPointOfBoundingBoxForBalls.getY()))
        m_maxPointOfBoundingBoxForBalls.setY(y + radius);
    if (rg_GT(z + radius, m_maxPointOfBoundingBoxForBalls.getZ()))
        m_maxPointOfBoundingBoxForBalls.setZ(z + radius);

    if (rg_LT(x - radius, m_minPointOfBoundingBoxForBalls.getX()))
        m_minPointOfBoundingBoxForBalls.setX(x - radius);
    if (rg_LT(y - radius, m_minPointOfBoundingBoxForBalls.getY()))
        m_minPointOfBoundingBoxForBalls.setY(y - radius);
    if (rg_LT(z - radius, m_minPointOfBoundingBoxForBalls.getZ()))
        m_minPointOfBoundingBoxForBalls.setZ(z - radius);


    if (rg_GT(radius, m_maxRadiusOfBalls))
        m_maxRadiusOfBalls = radius;
    if (rg_LT(radius, m_minRadiusOfBalls))
        m_minRadiusOfBalls = radius;

    rg_REAL sumOfRadius = m_AvgRadiusOfBalls*m_GlobalGeneratorList.getSize();

    BallGenerator ball(m_GlobalGeneratorList.getSize(), sphere, property, IDFromInput);

    ////////////////////////////////////////////////////////////////////////////
    // Youngsong 2017.12.05 
    BallGenerator* ptrBall = m_GlobalGeneratorList.addTail(ball);
    VDCell*        ptrCell = m_GlobalCellList.addTail(VDCell(m_GlobalCellList.getSize()));
    connectCellAndGenerator(ptrCell, ptrBall);
    //
    ////////////////////////////////////////////////////////////////////////////


    if (ptrBall != rg_NULL) {
        m_AvgRadiusOfBalls = (sumOfRadius + radius) / m_GlobalGeneratorList.getSize();

        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_INT rg_SphereSetVoronoiDiagram::getNumOfVoronoiCells()
{
    rg_dList< VDCell >* globalCellList = getGlobalCellList();

    rg_INT numVDCell = 0;
    VDCell* currCell = rg_NULL;
    globalCellList->reset4Loop();
    while ( globalCellList->setNext4Loop() )  {
        currCell = globalCellList->getpEntity();
        if ( currCell->getGenerator() == rg_NULL )
            continue;

        numVDCell++;
    }

    return numVDCell;
}



rg_INT rg_SphereSetVoronoiDiagram::getNumOfUnboundedVoronoiCells()
{
    /*
    rg_dList< VDCell >* globalCellList = getGlobalCellList();

    rg_INT numVDCell = 0;
    VDCell* currCell = rg_NULL;
    globalCellList->reset4Loop();
    while ( globalCellList->setNext4Loop() )  {
        currCell = globalCellList->getpEntity();
        if ( currCell->getGenerator() == rg_NULL )
            continue;

        if ( currCell->isBounded() == rg_TRUE )
            continue;

        numVDCell++;
    }
    */

    rg_dList< VDCell >* globalCellList = getGlobalCellList();
    VDCell* infinityCell     = globalCellList->getLastpEntity();
    rg_INT  numUnboundedCell = infinityCell->getNumOfBoundingFaces();

    return numUnboundedCell;
}



rg_INT rg_SphereSetVoronoiDiagram::getNumOfVoronoiFaces()
{
    rg_dList< VDFace >* globalFaceList = getGlobalFaceList();

    rg_INT numVDFace = 0;
    VDFace* currFace = rg_NULL;
    globalFaceList->reset4Loop();
    while ( globalFaceList->setNext4Loop() )  {
        currFace = globalFaceList->getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;

        numVDFace++;
    }

    return numVDFace;
}



rg_INT rg_SphereSetVoronoiDiagram::getNumOfVoronoiEdges()
{
    rg_dList< VDEdge >* globalEdgeList = getGlobalEdgeList();

    rg_INT numVDEdge = 0;
    VDEdge* currEdge = rg_NULL;
    globalEdgeList->reset4Loop();
    while ( globalEdgeList->setNext4Loop() )  {
        currEdge = globalEdgeList->getpEntity();

        if ( currEdge->isOnInfinity() == rg_TRUE  )
            continue;

        numVDEdge++;
    }

    return numVDEdge;
}



rg_INT rg_SphereSetVoronoiDiagram::getNumOfVoronoiVertices()
{
    rg_dList< VDVertex >* globalVertexList = getGlobalVerticeList();

    rg_INT numVDVertex = 0;
    VDVertex* currVertex = rg_NULL;
    globalVertexList->reset4Loop();
    while ( globalVertexList->setNext4Loop() )  {
        currVertex = globalVertexList->getpEntity();

        if ( currVertex->isOnInfinity() == rg_TRUE )
            continue;

        numVDVertex++;
    }

    return numVDVertex;
}



void rg_SphereSetVoronoiDiagram::duplicate(const rg_SphereSetVoronoiDiagram& origin)
{
    duplicateTopology( origin );

    duplicateGenerators( origin );
}

    

void rg_SphereSetVoronoiDiagram::duplicateSubsetByComplement(const rg_SphereSetVoronoiDiagram& origin, const rg_dList<VDEdge*>& edgesInSubset)
{
    //  map< origin_entity*, this_entity*>
    map<VDCell*, VDCell*>               cellMap;
    map<VDFace*, VDFace*>               faceMap;
    map<VDLoop*, VDLoop*>               loopMap;
    map<VDPartialEdge*, VDPartialEdge*> prEdgeMap;
    map<VDEdge*, VDEdge*>               edgeMap;
    map<VDVertex*, VDVertex*>           vertexMap;

    // make binary search tree with origin.entity as key
    makeEntityMapsForSubset( edgesInSubset, cellMap, faceMap, loopMap, prEdgeMap, edgeMap, vertexMap);

    setTopologyOfVerticesInSubset(     vertexMap, edgeMap );
    setTopologyOfEdgesInSubset(        prEdgeMap, edgeMap,   vertexMap);
    setTopologyOfPartialEdgesInSubset( loopMap,   prEdgeMap, edgeMap);
    setTopologyOfLoopsInSubset(        faceMap,   loopMap,   prEdgeMap);
    setTopologyOfFacesInSubset(        cellMap,   faceMap,   loopMap);
    setTopologyOfCellsInSubset(        cellMap,   faceMap);
}


    
void rg_SphereSetVoronoiDiagram::clean()
{
    VoronoiDiagram3D::clean();

    m_GlobalGeneratorList.removeAll();
}



void rg_SphereSetVoronoiDiagram::duplicateGenerators(const rg_SphereSetVoronoiDiagram& origin)
{
    //  The origin.VDCell_ID is identical to this.VDCell_ID.
    //  map< origin.VDCell_ID, this.BallGenerator*>
    map<rg_INT, BallGenerator*> ballMap;

    //  make binary search tree for ball.
    origin.m_GlobalGeneratorList.reset4Loop();
    while ( origin.m_GlobalGeneratorList.setNext4Loop() ) {
        BallGenerator* currOriBall  = origin.m_GlobalGeneratorList.getpEntity();
        BallGenerator* currThisBall = m_GlobalGeneratorList.add( *currOriBall );

        ballMap.insert( make_pair(currOriBall->getCell()->getID(), currThisBall) );
    }


    map<rg_INT, BallGenerator*>::const_iterator i_ball;
    map<rg_INT, BallGenerator*>::const_iterator i_empty = ballMap.end();

    m_GlobalCellList.reset4Loop();
    while ( m_GlobalCellList.setNext4Loop() ) {
        VDCell* currCell = m_GlobalCellList.getpEntity();

        i_ball = ballMap.find( currCell->getID() );
        if ( i_ball == i_empty ) {
            currCell->setGenerator( rg_NULL );
        }
        else {
            BallGenerator* currBall = i_ball->second;

            currCell->setGenerator( currBall );
            currBall->setCell( currCell );
        }
    }
}



void rg_SphereSetVoronoiDiagram::perturbGenerators(map<BallGenerator*, rg_Point3D>& geneCenterMap)
{
    VoronoiDiagram3D::cleanForRestart();

    random_device   rd;
    mt19937         gen(rd());
    uniform_real_distribution<> distr(0.0, resNeg6);

    map<BallGenerator*, rg_Point3D>::iterator i_gen;
    for (i_gen = geneCenterMap.begin(); i_gen != geneCenterMap.end(); i_gen++) {
        BallGenerator* currGen = i_gen->first;
        rg_Point3D     center = i_gen->second;

        rg_REAL x = center.getX() + distr(gen);
        rg_REAL y = center.getY() + distr(gen);
        rg_REAL z = center.getZ() + distr(gen);

        currGen->setCenter(rg_Point3D(x, y, z));

    }


//    rg_REAL randomMax = RAND_MAX;
//    srand((unsigned)time(NULL));

 //   map<BallGenerator*, rg_Point3D>::iterator i_gen;
 //   for ( i_gen = geneCenterMap.begin(); i_gen!=geneCenterMap.end(); i_gen++ ) {
 //       BallGenerator* currGen = i_gen->first;
 //       rg_Point3D     center  = i_gen->second;

	//	rg_REAL x = center.getX() + rand()/randomMax*resNeg6;
	//	rg_REAL y = center.getY() + rand()/randomMax*resNeg6;
	//	rg_REAL z = center.getZ() + rand()/randomMax*resNeg6;

 //       currGen->setCenter( rg_Point3D(x, y, z) );

 //   }
}



void rg_SphereSetVoronoiDiagram::makeEntityMapsForSubset( 
                                     //const rg_REAL& betaValue,
                                     const rg_dList<VDEdge*>&             edgesOnComponent, 
                                     map<VDCell*, VDCell*>&               cellMap,
                                     map<VDFace*, VDFace*>&               faceMap,
                                     map<VDLoop*, VDLoop*>&               loopMap,
                                     map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                     map<VDEdge*, VDEdge*>&               edgeMap,
                                     map<VDVertex*, VDVertex*>&           vertexMap)
{
    rg_INT i=0;

    edgesOnComponent.reset4Loop();
    while ( edgesOnComponent.setNext4Loop() ) {
        VDEdge* currEdge = edgesOnComponent.getEntity();

        //  for vertex
        VDVertex* currVtx[2] = { currEdge->getStartVertex(), currEdge->getEndVertex() };
        for ( i=0; i<2; i++ ) {
            if ( vertexMap.find( currVtx[i] ) != vertexMap.end() ) continue;
            
            VDVertex* com_vtx = createVertex( VDVertex( currVtx[i]->getID() ) );
            vertexMap.insert( make_pair( currVtx[i], com_vtx ) );
        }

        //  for edge
        VDEdge* com_edge = createEdge( VDEdge( currEdge->getID(), rg_NULL, rg_NULL  ) );
        edgeMap.insert( make_pair( currEdge, com_edge ) );


        //  for partial edge
        VDPartialEdge* currPrEdge[3] = { currEdge->getPartialEdge(), 
                                         currPrEdge[0]->getNextPartialEdgeInRadialCycle(), 
                                         currPrEdge[1]->getNextPartialEdgeInRadialCycle() }; 
        for ( i=0; i<3; i++ ) {
            if ( prEdgeMap.find( currPrEdge[i] ) !=  prEdgeMap.end() ) continue;

            VDPartialEdge* com_prEdge = createPrEdge( VDPartialEdge( currPrEdge[i]->getID(), rg_NULL ) );
            prEdgeMap.insert( make_pair( currPrEdge[i], com_prEdge ) );
        }


        //  for loop
        VDLoop* currLoop[3] = { currPrEdge[0]->getLoop(), currPrEdge[1]->getLoop(), currPrEdge[2]->getLoop() };
        for ( i=0; i<3; i++ ) {
            if ( currLoop[i] == rg_NULL ) continue;
            if ( loopMap.find( currLoop[i] ) !=  loopMap.end() ) continue;

            VDLoop* com_loop = createLoop( VDLoop( currLoop[i]->getID(), rg_NULL ) );
            loopMap.insert( make_pair( currLoop[i], com_loop ) );
        }


        //  for face
        VDFace* currFace[3] = { rg_NULL, rg_NULL, rg_NULL };
        for ( i=0; i<3; i++ ) {
            if ( currLoop[i] != rg_NULL ) {
                currFace[i] = currLoop[i]->getFace();
            }

            if ( currFace[i] == rg_NULL ) continue;
            if ( faceMap.find( currFace[i] ) !=  faceMap.end() ) continue;

            VDFace* com_face = createFace( VDFace( currFace[i]->getID() ) );
            faceMap.insert( make_pair( currFace[i], com_face ) );
        }


        //  for cell
        for ( i=0; i<3; i++ ) {
            if ( currFace[i] == rg_NULL ) continue;

            VDCell* currCell[2] = { currFace[i]->getLeftCell(), currFace[i]->getRightCell() };
            for ( rg_INT j=0; j<2; j++ ) {
                if ( currCell[j] == rg_NULL ) continue;
                if ( cellMap.find( currCell[j] ) !=  cellMap.end() ) continue;

                VDCell* com_cell = createCell( VDCell( currCell[j]->getID() ) );
                cellMap.insert( make_pair( currCell[j], com_cell ) );
            }
        }
    }



    //  for ball generator
    map<VDCell*, VDCell*>::iterator i_cell = cellMap.begin();
    for ( ; i_cell!=cellMap.end(); i_cell++ ) {
        VDCell* currCell = i_cell->first;
        VDCell* com_cell = i_cell->second;

        BallGenerator* currBall = currCell->getGenerator();
        if ( currBall != rg_NULL ) {
            BallGenerator* com_ball = m_GlobalGeneratorList.addTail( *currBall );

            com_cell->setGenerator( com_ball );
            com_ball->setCell( com_cell );
        }
    }
}



void rg_SphereSetVoronoiDiagram::setTopologyOfVerticesInSubset( 
                                          const map<VDVertex*, VDVertex*>&           vertexMap,
                                          const map<VDEdge*, VDEdge*>&               edgeMap )
{
    map<VDEdge*, VDEdge*>::const_iterator EDGE_NOT_FOUND = edgeMap.end();

    map<VDVertex*, VDVertex*>::const_iterator i_vtx = vertexMap.begin();
    for ( ; i_vtx != vertexMap.end(); i_vtx++ ) {
        VDVertex* currOriVtx  = i_vtx->first;
        VDVertex* currThisVtx = i_vtx->second;

        for ( rg_INT i=0; i<4; i++ ) {
            VDEdge*   incidentEdge = rg_NULL;
            map<VDEdge*, VDEdge*>::const_iterator i_edge = edgeMap.find( currOriVtx->getIncidentEdge(i) );
            if ( i_edge != EDGE_NOT_FOUND ) {
                incidentEdge = i_edge->second;
            }
            
            currThisVtx->setIncidentEdge(i, incidentEdge );
        }

        currThisVtx->isOnInfinity(             currOriVtx->isOnInfinity() );
        currThisVtx->setPoint(                 currOriVtx->getPoint() );
        currThisVtx->setRadiusOfTangentSphere( currOriVtx->getRadiusOfTangentSphere() );

        currThisVtx->connectBetaCell( currOriVtx->getBetaCell() );
        currThisVtx->connectQCell(    currOriVtx->getQCell() );
    }
}



void rg_SphereSetVoronoiDiagram::setTopologyOfEdgesInSubset(    
                                          const map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                          const map<VDEdge*, VDEdge*>&               edgeMap,
                                          const map<VDVertex*, VDVertex*>&           vertexMap)
{
    map<VDEdge*, VDEdge*>::const_iterator i_edge = edgeMap.begin();
    for ( ; i_edge != edgeMap.end(); i_edge++ ) {
        VDEdge* currOriEdge  = i_edge->first;
        VDEdge* currThisEdge = i_edge->second;

        VDVertex*      startVertex = vertexMap.find( currOriEdge->getStartVertex() )->second;
        VDVertex*      endVertex   = vertexMap.find( currOriEdge->getEndVertex()   )->second;
        VDPartialEdge* prEdge      = prEdgeMap.find( currOriEdge->getPartialEdge() )->second;

        currThisEdge->setEdgeType(    currOriEdge->getEdgeType() );
        currThisEdge->setStartVertex( startVertex );
        currThisEdge->setEndVertex(   endVertex );
        currThisEdge->setPartialEdge( prEdge );
        if ( currThisEdge->isOnInfinity() ) {
            currThisEdge->setInfinity();
        }

        currThisEdge->connectBetaFace( currOriEdge->getBetaFace() );
    }
}



void rg_SphereSetVoronoiDiagram::setTopologyOfPartialEdgesInSubset( 
                                              const map<VDLoop*, VDLoop*>&               loopMap,
                                              const map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                              const map<VDEdge*, VDEdge*>&               edgeMap)
{
    map<VDPartialEdge*, VDPartialEdge*>::const_iterator i_pre;
    map<VDPartialEdge*, VDPartialEdge*>::const_iterator PREDGE_NOT_FOUND = prEdgeMap.end();
    map<VDLoop*, VDLoop*>::const_iterator               LOOP_NOT_FOUND   = loopMap.end();


    map<VDPartialEdge*, VDPartialEdge*>::const_iterator i_prEdge;
    for ( i_prEdge = prEdgeMap.begin(); i_prEdge != prEdgeMap.end(); i_prEdge++ ) {
        VDPartialEdge* currOriPrEdge  = i_prEdge->first;
        VDPartialEdge* currThisPrEdge = i_prEdge->second;

        VDLoop* currLoop = rg_NULL;
        map<VDLoop*, VDLoop*>::const_iterator i_loop = loopMap.find( currOriPrEdge->getLoop() );
        if ( i_loop != LOOP_NOT_FOUND ) {
            currLoop = i_loop->second;
        }

        VDEdge* currEdge = edgeMap.find( currOriPrEdge->getOriginalEdge() )->second;

        VDPartialEdge* nextPrEdgeInRCycle = prEdgeMap.find( currOriPrEdge->getNextPartialEdgeInRadialCycle() )->second;

        VDPartialEdge* nextPrEdgeInLoop   = rg_NULL;
        i_pre = prEdgeMap.find( currOriPrEdge->getNextPartialEdgeInLoop() );
        if ( i_pre != PREDGE_NOT_FOUND ) {
            nextPrEdgeInLoop = i_pre->second;
        }

        VDPartialEdge* prevPrEdgeInLoop   = rg_NULL;
        i_pre = prEdgeMap.find( currOriPrEdge->getPreviousPartialEdgeInLoop() );
        if ( i_pre != PREDGE_NOT_FOUND ) {
            prevPrEdgeInLoop = i_pre->second;
        }

        currThisPrEdge->setPartialEdge( currLoop, currEdge, 
                                        nextPrEdgeInRCycle, nextPrEdgeInLoop, prevPrEdgeInLoop, 
                                        currOriPrEdge->isRightOrientationInLoop() );
    }
}



void rg_SphereSetVoronoiDiagram::setTopologyOfLoopsInSubset(    
                                          const map<VDFace*, VDFace*>&               faceMap,
                                          const map<VDLoop*, VDLoop*>&               loopMap,
                                          const map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap)
{
    map<VDPartialEdge*, VDPartialEdge*>::const_iterator i_prEdge;
    map<VDPartialEdge*, VDPartialEdge*>::const_iterator PREDGE_NOT_FOUND = prEdgeMap.end();


    map<VDLoop*, VDPartialEdge*> firstPrEdgeMap;
    for ( i_prEdge=prEdgeMap.begin(); i_prEdge!=prEdgeMap.end(); i_prEdge++ ) {
        VDPartialEdge* currThisPrEdge = i_prEdge->second;
        VDLoop*        currThisLoop   = currThisPrEdge->getLoop();

        if ( firstPrEdgeMap.find( currThisLoop ) == firstPrEdgeMap.end() ) {
            firstPrEdgeMap.insert( make_pair( currThisLoop, currThisPrEdge ) );
        }
    }


    map<VDLoop*, VDLoop*>::const_iterator i_loop = loopMap.begin();
    for ( ; i_loop != loopMap.end(); i_loop++ ) {
        VDLoop* currOriLoop  = i_loop->first;
        VDLoop* currThisLoop = i_loop->second;

        VDFace*        currFace   = faceMap.find( currOriLoop->getFace() )->second;
        VDPartialEdge* currPrEdge = rg_NULL;
        i_prEdge = prEdgeMap.find( currOriLoop->getPartialEdge() );
        if ( i_prEdge != PREDGE_NOT_FOUND ) {
            currPrEdge = i_prEdge->second;
        }
        else {
            currPrEdge = firstPrEdgeMap.find( currThisLoop )->second;
            //rg_dList<VDPartialEdge*> boundingPartialEdges;
            //currOriLoop->collectBoundingPartialEdgesInCCW( boundingPartialEdges );

            //boundingPartialEdges.reset4Loop();
            //while ( boundingPartialEdges.setNext4Loop() ) {
            //    VDPartialEdge* prEdge = boundingPartialEdges.getEntity();

            //    i_prEdge = prEdgeMap.find( prEdge );
            //    if ( i_prEdge != end_prEdge ) {
            //        currPrEdge = i_prEdge->second;
            //        break;
            //    }
            //}
        }

        currThisLoop->setLoop( currFace, currPrEdge, currOriLoop->isOuterLoop() );
    }
}



void rg_SphereSetVoronoiDiagram::setTopologyOfFacesInSubset(    
                                          const map<VDCell*, VDCell*>&               cellMap,
                                          const map<VDFace*, VDFace*>&               faceMap,
                                          const map<VDLoop*, VDLoop*>&               loopMap)
{
    map<VDLoop*, VDLoop*>::const_iterator LOOP_NOT_FOUND = loopMap.end();

    map<VDFace*, VDFace*>::const_iterator i_face = faceMap.begin();
    for ( ; i_face != faceMap.end(); i_face++ ) {
        VDFace* currOriFace  = i_face->first;
        VDFace* currThisFace = i_face->second;

        VDCell* leftCell  = cellMap.find( currOriFace->getLeftCell() )->second;
        VDCell* rightCell = cellMap.find( currOriFace->getRightCell() )->second;
        currThisFace->setLeftCell(  leftCell );
        currThisFace->setRightCell( rightCell );

        if ( currOriFace->isOnInfinity() ) currThisFace->setInfinity();
        currThisFace->isBounded(   currOriFace->isBounded() );


        rg_dList<VDLoop*>* originLoops = currOriFace->getLoops();
        originLoops->reset4Loop();
        while ( originLoops->setNext4Loop() ) {
            VDLoop* currOriLoop = originLoops->getEntity();

            map<VDLoop*, VDLoop*>::const_iterator i_loop = loopMap.find( currOriLoop );
            if ( i_loop != LOOP_NOT_FOUND ) {
                VDLoop* currThisLoop = i_loop->second;
                currThisFace->addLoop( currThisLoop );
            }
        }

        currThisFace->connectBetaEdge( currOriFace->getBetaEdge() );
    }
}



void rg_SphereSetVoronoiDiagram::setTopologyOfCellsInSubset(    
                                          const map<VDCell*, VDCell*>&               cellMap,
                                          const map<VDFace*, VDFace*>&               faceMap)
{
    map<VDCell*, VDCell*>::const_iterator i_cell = cellMap.begin();
    for ( ; i_cell != cellMap.end(); i_cell++ ) {
        VDCell* currOriCell  = i_cell->first;
        VDCell* currThisCell = i_cell->second;

        currThisCell->connectBetaVertex( currOriCell->getBetaVertex() );
        currThisCell->connectQVertex(    currOriCell->getQVertex() );

    }

    map<VDFace*, VDFace*>::const_iterator i_face = faceMap.begin();
    for ( ; i_face != faceMap.end(); i_face++ ) {
        VDFace* currThisFace = i_face->second;

        VDCell* leftCell  = currThisFace->getLeftCell();
        VDCell* rightCell = currThisFace->getRightCell();
    
        leftCell->addBoundingFace( currThisFace );
        rightCell->addBoundingFace( currThisFace );
    }
}




///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..


///////////////////////////////////////////////////////////////////////////////
//
//  constructing functions of sphere set voronoi diagram
rg_INT rg_SphereSetVoronoiDiagram::constructVoronoiDiagramOfSpheresByEdgeTracing()
{
    rg_INT statusOfVDS;

    map<BallGenerator*, rg_Point3D> geneCenterMap;
    m_GlobalGeneratorList.reset4Loop();
    while ( m_GlobalGeneratorList.setNext4Loop() ) {
        BallGenerator* currGenerator = m_GlobalGeneratorList.getpEntity();
        rg_Point3D center = currGenerator->getCenter();

        geneCenterMap.insert( make_pair(currGenerator, center) );
    }
    perturbGenerators(geneCenterMap);

    rg_INT NumRestart = 10;
    rg_INT count      = 0;
    do {
        statusOfVDS = constructSphereVoronoiDiagramRobustly();

        if ( statusOfVDS == NORMAL_STATE ) {
            break;
        }
        else if ( statusOfVDS == VERTEX_DEGREE_VIOLATION ) {
            count++;
            perturbGenerators( geneCenterMap );
        }
    } while ( count < NumRestart );

    if ( statusOfVDS == NORMAL_STATE ) {
        setGeometry();
    }

    // test perturbation
    //ofstream fout;
    //fout.open("test.txt", ios_base::app);
    //fout << count << endl;
    //fout.close();

    //statusOfVDS = constructSphereVoronoiDiagramByAcceleratedEdgeTracing();

    //start = clock();
    //statusOfVDS = constructSphereSetVoronoiDiagramByEdgeTracing();	
    //finish = clock();
    //computingTime[16] = finish - start;


    if ( ON_DEBUG ) {
        foutEVDS.close();
    }


    return statusOfVDS;
}



///////////////////////////////////////////////////////////////////////////////
//
//  Robust Edge-tracing Algorithm
//      by Youngsong Cho (2006. 7. 20)
//

rg_INT rg_SphereSetVoronoiDiagram::constructSphereVoronoiDiagramRobustly()
{
    //_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

    setBucket( m_sizeOfBucketElement );

    //_CrtDumpMemoryLeaks();

    //  construct BIG world.
    constructBigWorldRobustly();

    //  construct SMALL worlds.
    switch ( m_status )  {
        case NORMAL_STATE:
            constructSmallWorldRobustly();
            break;

        case NOT_FOUND_INITIAL_VERTEX:
            makeVirtualTwoBigBrothers();
            constructSmallWorldRobustly();
            break;

        default:
            break;
    }

    //  we have to deal with the disconnected big worlds.

    if ( m_status != VERTEX_DEGREE_VIOLATION ) 
        updateGlobalTopology();

    return m_status;
}



void rg_SphereSetVoronoiDiagram::constructBigWorldRobustly()
{
    rg_INT initialSizeOfVIDIC = 101;

    VIDIC            theVIDIC( initialSizeOfVIDIC );
    rg_dList< Gate > theGateStack;

    if ( locateInitialVoronoiVertexInBigWorldRobustly( theVIDIC, theGateStack ) == rg_FALSE )  {
        m_status = NOT_FOUND_INITIAL_VERTEX;
        return;
    }

    constructThisWorldRobustly( rg_NULL, theVIDIC, theGateStack );

    if ( ON_DEBUG )  {
        foutEVDS << endl << endl;
        foutEVDS << "complete the construction of big world!" << endl << endl;
    }
}



rg_FLAG rg_SphereSetVoronoiDiagram::locateInitialVoronoiVertexInBigWorldRobustly(VIDIC& theVIDIC, rg_dList<Gate>& theGateStack )
{
    rg_FLAG bSuccess = rg_FALSE;

    VIDIC            tempVIDIC( 31 );
    rg_dList< Gate > tempGateStack;

    //  1. set the configuration of virtual initial vertex
    makeBallConfigurationForVirtualVertex(tempVIDIC, tempGateStack );


    Sphere initialTangentSphere;
    Gate   initialGates[NUM_DEFINING_CELLS_OF_ONE_VERTEX];


    while ( tempGateStack.getSize() > 0 )
    {
        Gate currGate = tempGateStack.popFront();

        if ( currGate.getStartBall()->getID() == -4 
             && currGate.getGateBallCCW1()->getID() == 4 
             && currGate.getGateBallCCW2()->getID() == 8 
             && currGate.getGateBallCCW3()->getID() == 6 )
             int aaa = 1;


        findEndVertexRobustly( currGate );

        if ( currGate.getNumOfEndBall() == 1 )  {
            BallGenerator* endBall = currGate.getSingleEndBall();
            Sphere  tangentSphereByEndBall( currGate.getEmptyTangentSphereAtEndVertex() );

            if (    currGate.getGateBallCCW1()->getCell()->getID() > -2 && currGate.getGateBallCCW2()->getCell()->getID() > -2
                 && currGate.getGateBallCCW3()->getCell()->getID() > -2 && endBall->getCell()->getID() > -2 )
            {
                BallGenerator* gateBall[3] = { currGate.getGateBallCCW1(), currGate.getGateBallCCW2(), currGate.getGateBallCCW3()};

                initialTangentSphere = tangentSphereByEndBall;
                initialGates[0].setGate( gateBall[2], gateBall[1], gateBall[0], endBall, rg_NULL );
                initialGates[1].setGate( endBall, gateBall[0], gateBall[1], gateBall[2], rg_NULL );
                initialGates[2].setGate( endBall, gateBall[1], gateBall[2], gateBall[0], rg_NULL );
                initialGates[3].setGate( endBall, gateBall[2], gateBall[0], gateBall[1], rg_NULL );

                bSuccess = rg_TRUE;


                if ( ON_DEBUG )
                {
                    foutEVDS << "*** FIND INITIAL VERTEX" << endl;
                    foutEVDS << currGate.getStartBall()->getID()   << "\t"
                             << gateBall[0]->getID()           << "\t"
                             << gateBall[1]->getID()           << "\t"
                             << gateBall[2]->getID()           << "\t"
                             << endBall->getID()               << "\t";

                    rg_FLAG edgeType = currGate.getEdgeType();
                    if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
                        foutEVDS << "Para" << "\t";
                    else if ( edgeType == LINEAR_EDGE )
                        foutEVDS << "Line" << "\t";
                    else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
                        foutEVDS << "Elli" << "\t";
                    else 
                        foutEVDS << "Non" << "\t";

                    foutEVDS << endl;

                }
            
                break;
            }
        }

        switch ( currGate.getNumOfEndBall() )  {
            case 0:
                defineUnboundedEdge(rg_NULL, tempVIDIC, tempGateStack, currGate);
                break;
            case 1:
                defineBoundedEdge(rg_NULL, tempVIDIC, tempGateStack, currGate);
                break;
            default:
                solveDegeneracyAtEndVertexOfCurrentEdge(rg_NULL, tempVIDIC, tempGateStack, currGate);
                break;
        }
    }

    //  3. remove virtual initial vertex and its topology

	rg_INT i=0;
	for ( i=0; i<4; i++)  {
        m_GlobalCellList.popFront();
        m_GlobalGeneratorList.popFront();
    }

    VDCell* currCell = rg_NULL;
	m_GlobalCellList.reset4Loop();
	while( m_GlobalCellList.setNext4Loop() )
	{
		currCell = m_GlobalCellList.getpEntity();
        currCell->getBoungindFaces()->removeAll();
    }
    m_GlobalFaceList.removeAll();
    m_GlobalLoopList.removeAll();
    m_GlobalEdgeList.removeAll();
    m_GlobalPartialedgeList.removeAll();
    m_GlobalVertexList.removeAll();


    if ( bSuccess )  {    
        m_status = NORMAL_STATE;
        //  4. set valid initial vertex.
        //  insert initial voronoi vertex into GVL(Global Vertex List).
        VDVertex* ptrVertex = makeVoronoiVertex( rg_FALSE, initialTangentSphere );

        //  make and insert 4 gates of initial voronoi vertex into Gate-Stack.
        Gate* ptrGate = rg_NULL;        
        for ( i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        {
            ptrGate = theGateStack.addTail( initialGates[i] );

            ptrGate->setStartVertex( ptrVertex );
            ptrGate->setNodePosInStack( theGateStack.getTail() );

            ptrVertex->setGate( i, ptrGate );
        }   

        //  insert vertex configuration of initial voronoi vertex into VIDIC.
        theVIDIC.insertVIDICItem( VIDICItem( initialGates[0].getGateBallCCW1(),
                                             initialGates[0].getGateBallCCW2(),
                                             initialGates[0].getGateBallCCW3(),
                                             initialGates[0].getStartBall(),
                                             ptrVertex) 
                                 );
    }

    return bSuccess;
}



void rg_SphereSetVoronoiDiagram::makeBallConfigurationForVirtualVertex(VIDIC& theVIDIC, rg_dList<Gate>& theGateStack )
{
    //  1. set the configuration of virtual initial vertex
    rg_Point3D virtualVertexCoord(1.0, 1.0, 1.0);
    rg_REAL    virtualVertexRadius  = 0.73205080756887719;    

    rg_Point3D centerOfVirtualBall[4];
    centerOfVirtualBall[0].setPoint(2.0, 0.0, 0.0);
    centerOfVirtualBall[1].setPoint(0.0, 2.0, 0.0);
    centerOfVirtualBall[2].setPoint(0.0, 0.0, 2.0);
    centerOfVirtualBall[3].setPoint(0.0, 0.0, 0.0);

    rg_Point3D sizeOfBoundingBox = m_maxPointOfBoundingBoxForBalls - m_minPointOfBoundingBoxForBalls;
    rg_REAL    maxSize = sizeOfBoundingBox.getX();
    if ( maxSize < sizeOfBoundingBox.getY() )
        maxSize = sizeOfBoundingBox.getY();
    if ( maxSize < sizeOfBoundingBox.getZ() )
        maxSize = sizeOfBoundingBox.getZ();
    rg_Point3D translation(-maxSize, -maxSize, -maxSize);

    translation = m_minPointOfBoundingBoxForBalls + translation;
    
    virtualVertexCoord = virtualVertexCoord + translation;


    BallGenerator* ptrVirtualBall[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
    VDCell*        ptrVirtualCell[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
    rg_INT i=0;
	for (i=0; i<4; i++)  {
        centerOfVirtualBall[i] = centerOfVirtualBall[i] + translation;

        ptrVirtualBall[i] = m_GlobalGeneratorList.addHead( BallGenerator(-2-i, centerOfVirtualBall[i], 1.0) );
        ptrVirtualCell[i] = m_GlobalCellList.addHead( VDCell( -2-i ) );
        connectCellAndGenerator( ptrVirtualCell[i], ptrVirtualBall[i] );
    }


    defineInitialVoronoiVertex( ptrVirtualBall, Sphere(virtualVertexCoord, virtualVertexRadius), theVIDIC, theGateStack);

}


rg_FLAG rg_SphereSetVoronoiDiagram::constructThisWorldRobustly(
    VDLoop* thisWorld, 
    VIDIC& theVIDIC, 
    rg_dList<Gate>& theGateStack) 
{
    while ( theGateStack.getSize() > 0 ) {

        Gate currGate = theGateStack.popFront();

        findEndVertexRobustly( currGate );

        switch ( currGate.getNumOfEndBall() )  {
            case 0:
                defineUnboundedEdge(thisWorld, theVIDIC, theGateStack, currGate);
                break;
            case 1:
                defineBoundedEdge(thisWorld, theVIDIC, theGateStack, currGate);
                break;
            default:
                solveDegeneracyAtEndVertexOfCurrentEdge(thisWorld, theVIDIC, theGateStack, currGate);
                break;
        }

        if ( m_status != NORMAL_STATE ) { 
            return rg_FALSE;
        }
    }

    return rg_TRUE;
}



void rg_SphereSetVoronoiDiagram::findEndVertexRobustly(Gate& currGate)
{
    currGate.obtainGeometricProperties();
    //currGate.determineEdgeType();

    switch ( currGate.getEdgeType() )  {
        case PARABOLIC_OR_HYPERBOLIC_EDGE:
            findEndVertexForParabolicEdgeRobustly(currGate);
            break;
        case LINEAR_EDGE:
            findEndVertexForLinearEdgeRobustly(currGate);
            break;
        case CIRCULAR_OR_ELLIPTIC_EDGE:
            findEndVertexForEllipticEdgeRobustly(currGate);
            break;
        default:
            break;
    }
}




void rg_SphereSetVoronoiDiagram::findEndVertexForParabolicEdgeRobustly(Gate& currGate)
{
    //  This gate defines a Voronoi edge.
    BallGenerator* gateBall[3] = { currGate.getGateBallCCW1(), currGate.getGateBallCCW2(), currGate.getGateBallCCW3() };
    Sphere         ball[3]     = { gateBall[0]->getBall(),     gateBall[1]->getBall(),     gateBall[2]->getBall() };

    ///////////////////////////////////////////////////////////////////////////
    //
    rg_dList<BallGenerator*> checkedBallSet;

    rg_INT i=0;
	for (i=0; i<3; i++)  {
        gateBall[i]->m_checkForBucket = rg_TRUE;
        checkedBallSet.addTail( gateBall[i] );
    }
    currGate.setCheckOfAllExcludedBalls(rg_TRUE);

    
    rg_REAL gatePlane[4];
    computePlaneTangentTo3SphereFromCCWOutside( ball, gatePlane );
    rg_Point3D normalOfGatePlane(gatePlane[0], gatePlane[1], gatePlane[2]);


    //////////////////////////////////////////////////////////////////////////
    //
    //  For Front Feasible region of parabolic and hyperbolic edge.
    //

    //  construct sphere enclosing R1

    rg_Point3D ptOnGatePlane[3];
    for ( i=0; i<3; i++)  {
        ptOnGatePlane[i] = ball[i].getCenter() + ( ball[i].getRadius()*normalOfGatePlane );    
    }
    rg_Point3D center = computeCenterOfApproximatingSphereForR1( 
                            ptOnGatePlane, currGate.getNormalOfCenterPlane(), currGate.getCoeffOfCenterPlane() );
    
    Sphere approxFrontFeasibleRegion(center, center.distance(ptOnGatePlane[0]) );


    rg_dList<BallGenerator*> ballSetForR1;
    m_bucket.getCellsByMask( ballSetForR1, approxFrontFeasibleRegion );
    
    currGate.findEmptyTangentSphereWithMinAngleForParabolicEdge(ballSetForR1);

    checkedBallSet.mergeTail( ballSetForR1 );

    if ( currGate.getNumOfEndBall() != 0 )  {
        Sphere emptyTangentSphere = currGate.getEmptyTangentSphereAtEndVertex();

        // if candidate tangent sphere is intersected with R1.
        if ( emptyTangentSphere.isContainedIn( approxFrontFeasibleRegion ) == rg_FALSE )  {
            //  A CANDIDATE end ball is found.
            //  Therefore, we mask bucket for TS_R1 (tangent sphere from R1)

            rg_dList<BallGenerator*> ballSetForTS_C;
            m_bucket.getCellsByMask( ballSetForTS_C, emptyTangentSphere );
            
            currGate.findEmptyTangentSphereWithMinAngleForParabolicEdge(ballSetForTS_C);

            checkedBallSet.mergeTail( ballSetForTS_C );
        }  
        else  {
        }
    }
    else  {
        
        //////////////////////////////////////////////////////////////////////////
        //
        //  For Rear Feasible region of parabolic and hyperbolic edge.
        //
        rg_dList< NavigatingBucketCell > queueOfCells;

        rg_Point3D unitLength;

        initializePropagationInRearFeasibleRegionOfNonEllipticEdge( queueOfCells, unitLength, ball, gatePlane );

        rg_INT  normalSign = classifyNormalSign(gatePlane[0], gatePlane[1], gatePlane[2]);

        rg_dList< BallGenerator* > allCandidateBalls;


        NavigatingBucketCell* currBCell = rg_NULL;  
        rg_dNode< NavigatingBucketCell >* currNode = queueOfCells.getHead();

        do 
	    {
            currBCell = currNode->getpEntity();
            rg_dList< BallGenerator* > candidateBalls;

            if ( currBCell->m_status == ON_PLANE )   
                m_bucket.getBalls( candidateBalls, currBCell->m_index, gatePlane );
            else
                m_bucket.getBalls( candidateBalls, currBCell->m_index );


            if ( candidateBalls.getSize() > 0 )  {

                currGate.findEmptyTangentSphereWithMinAngleForParabolicEdge(candidateBalls);

                allCandidateBalls.mergeTail( candidateBalls );
                
                if ( currGate.getNumOfEndBall() != 0 ){
                    rg_dList<BallGenerator*> ballSetForR2;
                    m_bucket.getCellsByMask( ballSetForR2, currGate.getEmptyTangentSphereAtEndVertex() );
                    currGate.findEmptyTangentSphereWithMinAngleForParabolicEdge(ballSetForR2);

                    checkedBallSet.mergeTail( ballSetForR2 );
                    break;
                }
            }
       
            if ( currBCell->m_status == ON_PLANE )   
                propagateBucketCellsOnPlane( *currBCell, unitLength, queueOfCells );
            else
                propagateBucketCellsInPlusPlane( *currBCell, normalSign, queueOfCells );

            currNode = currNode->getNext();

        } while( currNode != queueOfCells.getHead() );

        checkedBallSet.mergeTail( allCandidateBalls );

        queueOfCells.reset4Loop();
	    while( queueOfCells.setNext4Loop() )  {
	        currBCell = queueOfCells.getpEntity();
            m_bucket.setMark(currBCell->m_index, rg_FALSE);
        }
    }


    currGate.setCheckOfAllExcludedBalls(rg_FALSE);
    checkedBallSet.reset4Loop();
	while( checkedBallSet.setNext4Loop() )  {
	    checkedBallSet.getEntity()->m_checkForBucket = rg_FALSE;
    }
}



void rg_SphereSetVoronoiDiagram::findEndVertexForLinearEdgeRobustly(Gate& currGate)
{
    //  This gate defines a Voronoi edge.
    BallGenerator* gateBall[3] = { currGate.getGateBallCCW1(), currGate.getGateBallCCW2(), currGate.getGateBallCCW3() };
    Sphere         ball[3]     = { gateBall[0]->getBall(),     gateBall[1]->getBall(),     gateBall[2]->getBall() };

    ///////////////////////////////////////////////////////////////////////////
    //
    rg_dList<BallGenerator*> checkedBallSet;

    rg_INT i=0;
	for ( i=0; i<3; i++)  {
        gateBall[i]->m_checkForBucket = rg_TRUE;
        checkedBallSet.addTail( gateBall[i] );
    }
    currGate.setCheckOfAllExcludedBalls(rg_TRUE);

    //////////////////////////////////////////////////////////////////////////
    //
    //  For Front Feasible Region
    //

    //  construct sphere enclosing R1
    rg_Point3D normalOfCenterPlane = currGate.getNormalOfCenterPlane();
    rg_Point3D ptOnGatePlane = ball[0].getCenter() + ( ball[0].getRadius()*normalOfCenterPlane );
    rg_Point3D center        = currGate.getMinTangentSphere().getCenter();
    Sphere approxFrontFeasibleRegion(center, center.distance(ptOnGatePlane) );

    /*
    Plane centerPlane = currGate.getCenterPlane();
    rg_Point3D coordStartVertex = currGate.getStartVertex()->getPoint();
    rg_Point3D ptOnGatePlane = ball[0].getCenter() + ( ball[0].getRadius()*centerPlane.getNormal() );
    rg_Point3D center        = centerPlane.projectPointOnPlane(coordStartVertex);
    Sphere approxFrontFeasibleRegion(center, center.distance(ptOnGatePlane) );
    */

    //rg_REAL gatePlane[4] = { normalOfCenterPlane.getX(), normalOfCenterPlane.getY(), normalOfCenterPlane.getZ(), 
    //                         currGate.getCoeffOfCenterPlane() };
    rg_REAL gatePlane[4] = { normalOfCenterPlane.getX(), normalOfCenterPlane.getY(), normalOfCenterPlane.getZ(), 
                             -normalOfCenterPlane.innerProduct(ptOnGatePlane) };


    rg_dList<BallGenerator*> ballSetForR1;
    m_bucket.getCellsByMask( ballSetForR1, approxFrontFeasibleRegion );
    
    currGate.findEmptyTangentSphereWithMinAngleForLinearEdge(ballSetForR1);

    checkedBallSet.mergeTail( ballSetForR1 );


    if ( currGate.getNumOfEndBall() != 0 )  {

        Sphere emptyTangentSphere = currGate.getEmptyTangentSphereAtEndVertex();
        
        if ( emptyTangentSphere.isContainedIn( approxFrontFeasibleRegion ) == rg_FALSE ) {

            rg_dList<BallGenerator*> ballSetForTS_C;
            m_bucket.getCellsByMask( ballSetForTS_C, emptyTangentSphere );

            currGate.findEmptyTangentSphereWithMinAngleForLinearEdge(ballSetForTS_C);

            checkedBallSet.mergeTail( ballSetForTS_C );

        }  
        else  {
        }
    }
    else  {
    
        //////////////////////////////////////////////////////////////////////////
        //
        //  For Rear Feasible Region
        //

        rg_dList< NavigatingBucketCell > queueOfCells;

        rg_Point3D unitLength;

        initializePropagationInRearFeasibleRegionOfNonEllipticEdge( queueOfCells, unitLength, ball, gatePlane );

        rg_INT  normalSign = classifyNormalSign(gatePlane[0], gatePlane[1], gatePlane[2]);



        NavigatingBucketCell* currBCell = rg_NULL;
    
        rg_dNode< NavigatingBucketCell >* currNode = queueOfCells.getHead();
        //currNode = currNode->getNext();
        do 
	    {
            currBCell = currNode->getpEntity();
            rg_dList< BallGenerator* > candidateBalls;

            if ( currBCell->m_status == ON_PLANE )   
                m_bucket.getBalls( candidateBalls, currBCell->m_index, gatePlane );
            else  
                m_bucket.getBalls( candidateBalls, currBCell->m_index );

            if ( candidateBalls.getSize() > 0 )  {

                currGate.findEmptyTangentSphereWithMinAngleForLinearEdge(candidateBalls);

                checkedBallSet.mergeTail( candidateBalls );               

                if ( currGate.getNumOfEndBall() != 0 )  {
                    rg_dList<BallGenerator*> ballSetForR2;
                    m_bucket.getCellsByMask( ballSetForR2, currGate.getEmptyTangentSphereAtEndVertex() );
                    currGate.findEmptyTangentSphereWithMinAngleForLinearEdge(ballSetForR2);

                    checkedBallSet.mergeTail( ballSetForR2 );
                    break;
                }
            }


            if ( currBCell->m_status == ON_PLANE )   
                propagateBucketCellsOnPlane( *currBCell, unitLength, queueOfCells );
            else  
                propagateBucketCellsInPlusPlane( *currBCell, normalSign, queueOfCells );

            currNode = currNode->getNext();
        
        } while( currNode != queueOfCells.getHead() );


        queueOfCells.reset4Loop();
	    while( queueOfCells.setNext4Loop() )  {
	        currBCell = queueOfCells.getpEntity();
            m_bucket.setMark(currBCell->m_index, rg_FALSE);
        }
    }


    currGate.setCheckOfAllExcludedBalls(rg_FALSE);
    checkedBallSet.reset4Loop();
	while( checkedBallSet.setNext4Loop() )  {
	    checkedBallSet.getEntity()->m_checkForBucket = rg_FALSE;
    }
}



void rg_SphereSetVoronoiDiagram::findEndVertexForEllipticEdgeRobustly(Gate& currGate)
{
    //  This gate defines a Voronoi edge.
    BallGenerator* gateBall[3] = { currGate.getGateBallCCW1(), currGate.getGateBallCCW2(), currGate.getGateBallCCW3() };
    Sphere         tangentSphereInCP[2] = { currGate.getMinTangentSphere(0), currGate.getMinTangentSphere(1) };

    //  define box filter
    rg_Point3D vecTS1ToTS2 = tangentSphereInCP[1].getCenter() - tangentSphereInCP[0].getCenter();
    vecTS1ToTS2.normalize();

    rg_Point3D farPointOfTS1 = tangentSphereInCP[0].getCenter() - (tangentSphereInCP[0].getRadius()*vecTS1ToTS2);
    rg_Point3D farPointOfTS2 = tangentSphereInCP[1].getCenter() + (tangentSphereInCP[1].getRadius()*vecTS1ToTS2);

    rg_Point3D center = (farPointOfTS1 + farPointOfTS2)/2.0;
    rg_REAL    radius = center.distance( farPointOfTS1 );
    Sphere boxfilter(center, radius);    


    currGate.setCheckOfAllExcludedBalls(rg_TRUE);
    gateBall[0]->m_checkForBucket = rg_TRUE;
    gateBall[1]->m_checkForBucket = rg_TRUE;
    gateBall[2]->m_checkForBucket = rg_TRUE;


    rg_dList<BallGenerator*> ballSetForBoxfilter;
    m_bucket.getCellsByMask( ballSetForBoxfilter, boxfilter );

    currGate.findEmptyTangentSphereWithMinAngleForEllipticEdge(ballSetForBoxfilter);


    gateBall[0]->m_checkForBucket = rg_FALSE;
    gateBall[1]->m_checkForBucket = rg_FALSE;
    gateBall[2]->m_checkForBucket = rg_FALSE;
    currGate.setCheckOfAllExcludedBalls(rg_FALSE);


    ballSetForBoxfilter.reset4Loop();
	while( ballSetForBoxfilter.setNext4Loop() )
	    ballSetForBoxfilter.getEntity()->m_checkForBucket = rg_FALSE;

}




void rg_SphereSetVoronoiDiagram::defineUnboundedEdge(VDLoop* thisWorld, 
                                                                    VIDIC& theVIDIC, 
                                                                    rg_dList<Gate>& theGateStack, 
                                                                    const Gate& currGate)
{
    VDVertex* ptrEndVertex = makeVoronoiVertex( rg_TRUE, Sphere( DBL_MAX, DBL_MAX, DBL_MAX, 0.0) );
 
    if ( currGate.getStartVertex()->getIncidentEdge( 2 ) != rg_NULL )  {
        currGate.getStartVertex()->removeAllGates();
    }
    
    BallGenerator* gateBall[3] = { currGate.getGateBallCCW1(), currGate.getGateBallCCW2(), currGate.getGateBallCCW3() };

    updateLocalTopologyByEdge(thisWorld, currGate, ptrEndVertex);

    if ( ON_DEBUG )
    {
        foutEVDS << currGate.getStartBall()->getID()    << "\t"
                 << currGate.getGateBallCCW1()->getID() << "\t"
                 << currGate.getGateBallCCW2()->getID() << "\t"
                 << currGate.getGateBallCCW3()->getID() << "\t"
                 << "\t";                           

        rg_FLAG edgeType = currGate.getEdgeType();
        if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
            foutEVDS << "Para" << "\t";
        else if ( edgeType == LINEAR_EDGE )
            foutEVDS << "Line" << "\t";
        else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
            foutEVDS << "Elli" << "\t";
        else 
            foutEVDS << "Non" << "\t";

        foutEVDS << "IV"                               << "\t"
                 << currGate.getStartVertex()->getID() << "\t"
                 << ptrEndVertex->getID()              << endl;
    }

}



void rg_SphereSetVoronoiDiagram::defineBoundedEdge(VDLoop* thisWorld, 
                                                                  VIDIC&          theVIDIC, 
                                                                  rg_dList<Gate>& theGateStack, 
                                                                  const Gate&     currGate)
{
    VDVertex*      ptrEndVertex  = rg_NULL;
    
    BallGenerator* gateBall[3]   = { currGate.getGateBallCCW1(), currGate.getGateBallCCW2(), currGate.getGateBallCCW3() };
    BallGenerator* endBall       = currGate.getSingleEndBall();
    Sphere         tangentSphere = currGate.getEmptyTangentSphereAtEndVertex();


    VIDICItem* existingVertexConfig = theVIDIC.findVIDICItem( gateBall[0], gateBall[1], gateBall[2], endBall,
                                                              tangentSphere.getCenter() );

    if ( existingVertexConfig != rg_NULL )  {
        ptrEndVertex = existingVertexConfig->getVertex();


        /*
        if ( ON_DEBUG )
        {
            foutEVDS << currGate.getStartBall()->getID()   << "\t"
                     << gateBall[0]->getID()           << "\t"
                     << gateBall[1]->getID()           << "\t"
                     << gateBall[2]->getID()           << "\t"
                     << endBall->getID()               << "\t";

            rg_FLAG edgeType = currGate.getEdgeType();
            if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
                foutEVDS << "Para" << "\t";
            else if ( edgeType == LINEAR_EDGE )
                foutEVDS << "Line" << "\t";
            else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
                foutEVDS << "Elli" << "\t";
            else 
                foutEVDS << "Non" << "\t";

            foutEVDS << "PV"                           << "\t"
                     << currGate.getStartVertex()->getID() << "\t"
                     << ptrEndVertex->getID()          << "\t"
                     << ptrEndVertex->getPoint().getX() << "\t"
                     << ptrEndVertex->getPoint().getY() << "\t"
                     << ptrEndVertex->getPoint().getZ() << "\t"
                     << ptrEndVertex->getRadiusOfTangentSphere() << endl;
        }
        */

        //  for debug. by Youngsong Cho 2005.03.09
        if ( ptrEndVertex->getIncidentEdge( 3 ) != rg_NULL )  {
            /*
            if ( ON_DEBUG )
            {
                foutEVDS << " ********  S T O P **************" << endl;
                VDEdge** edge = ptrEndVertex->getAllIncidentEdges();
                rg_INT j=0;
				for ( j=0; j<4; j++ )  {
                    foutEVDS << "\t" 
                             << "e" << edge[j]->getID() << "\t"
                             << "v" << edge[j]->getStartVertex()->getID() << "\t"
                             << "v" << edge[j]->getEndVertex()->getID() << "\t";
                    VDCell* edgeCell[3]  = {rg_NULL, rg_NULL, rg_NULL};
                    VDCell* startCell[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
                    VDCell* endCell[4]   = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
                    edge[j]->inquireIntoOnlyEdgeSharingCellsInCCW( edgeCell );
                    if (edge[j]->getStartVertex()->getIncidentEdge( 3 ) != rg_NULL)
                        edge[j]->getStartVertex()->inquireAllCellsToDefineThisVertex( startCell );
                    if (edge[j]->getEndVertex()->getIncidentEdge( 3 ) != rg_NULL)
                        edge[j]->getEndVertex()->inquireAllCellsToDefineThisVertex( endCell );

                    VDCell* cell[5];
                    rg_INT k=0;
					for ( k=0; k<3; k++ )
                        cell[k+1] = edgeCell[k];
                    for ( k=0; k<4; k++ )  {
                        if ( startCell[k] != edgeCell[0] && startCell[k] != edgeCell[1] && startCell[k] != edgeCell[2] )  {
                            cell[0] = startCell[k];
                            break;
                        }
                    }
                    for ( k=0; k<4; k++ )  {
                        if ( endCell[k] != edgeCell[0] && endCell[k] != edgeCell[1] && endCell[k] != edgeCell[2] )  {
                            cell[4] = endCell[k];
                            break;
                        }
                    }
                    for ( k=0; k<5; k++ )  {
                        if ( cell[k] == rg_NULL )  {
                            foutEVDS << "?" << "\t";
                        }
                        else {
                            foutEVDS << cell[k]->getGenerator()->getID() << "\t";
                        }
                    }
                    foutEVDS << edge[j]->getStartVertex()->getPoint().getX() << "\t"
                             << edge[j]->getStartVertex()->getPoint().getY() << "\t"
                             << edge[j]->getStartVertex()->getPoint().getZ() << "\t"
                             << edge[j]->getStartVertex()->getRadiusOfTangentSphere() << "\t";
                    foutEVDS << edge[j]->getEndVertex()->getPoint().getX() << "\t"
                             << edge[j]->getEndVertex()->getPoint().getY() << "\t"
                             << edge[j]->getEndVertex()->getPoint().getZ() << "\t"
                             << edge[j]->getEndVertex()->getRadiusOfTangentSphere() << endl;

                    foutEVDS <<  endl;
                }
            }
            */
            m_status = VERTEX_DEGREE_VIOLATION;
            return;
        }

        //  remove redundant gate in Gate-Stack
        rg_INT i=0;
		for (i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)  {
            Gate* gateOfEndVertex = ptrEndVertex->getGate( i );

            if ( gateOfEndVertex == rg_NULL )  continue;

            if ( gateOfEndVertex->isTheSameGate(currGate, ptrEndVertex) )  {
                theGateStack.kill( gateOfEndVertex->getNodePosInStack() );
                ptrEndVertex->setGate( i, rg_NULL );
                break;
            }
        }           
        
        if ( ptrEndVertex->getIncidentEdge( 2 ) != rg_NULL )  {
            //theVIDIC.removeVIDICItem( *existingVertexConfig );
            ptrEndVertex->removeAllGates();
        }       
    }
    else  {
        ptrEndVertex = makeVoronoiVertex( rg_FALSE, tangentSphere );

        /*
        if ( ON_DEBUG )
        {
            foutEVDS << currGate.getStartBall()->getID()   << "\t"
                     << gateBall[0]->getID()           << "\t"
                     << gateBall[1]->getID()           << "\t"
                     << gateBall[2]->getID()           << "\t"
                     << endBall->getID()               << "\t";

            rg_FLAG edgeType = currGate.getEdgeType();
            if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
                foutEVDS << "Para" << "\t";
            else if ( edgeType == LINEAR_EDGE )
                foutEVDS << "Line" << "\t";
            else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
                foutEVDS << "Elli" << "\t";
            else 
                foutEVDS << "Non" << "\t";

            foutEVDS << "NV"                           << "\t"
                     << currGate.getStartVertex()->getID() << "\t"
                     << ptrEndVertex->getID()          << "\t"
                     << ptrEndVertex->getPoint().getX() << "\t"
                     << ptrEndVertex->getPoint().getY() << "\t"
                     << ptrEndVertex->getPoint().getZ() << "\t"
                     << ptrEndVertex->getRadiusOfTangentSphere() << endl;
        }
        */

        //  make and insert 3 gates of new voronoi vertex into Gate-Stack.
        Gate currentGates[ 3 ] 
                 = { Gate( endBall, gateBall[0], gateBall[1], gateBall[2], ptrEndVertex ),
                     Gate( endBall, gateBall[1], gateBall[2], gateBall[0], ptrEndVertex ),
                     Gate( endBall, gateBall[2], gateBall[0], gateBall[1], ptrEndVertex ) };

        rg_INT i=0;
		for (i=0; i<3; i++)  {
            Gate* ptrGate = theGateStack.addTail( currentGates[i] );
            ptrGate->setNodePosInStack( theGateStack.getTail() );

            ptrEndVertex->setGate( i, ptrGate );
        }   
        ptrEndVertex->setGate( 3, rg_NULL );

        //  insert vertex configuration of initial voronoi vertex into VIDIC.
        theVIDIC.insertVIDICItem( gateBall[0], gateBall[1], gateBall[2], endBall, ptrEndVertex, currGate.getNumVerticesBy4Balls() );
    }
    
    
    if ( currGate.getStartVertex()->getIncidentEdge( 2 ) != rg_NULL )  {
        currGate.getStartVertex()->removeAllGates();
    }
    
    updateLocalTopologyByEdge(thisWorld, currGate, ptrEndVertex);
}



void rg_SphereSetVoronoiDiagram::solveDegeneracyAtEndVertexOfCurrentEdge(VDLoop* thisWorld, 
                                                                         VIDIC& theVIDIC, 
                                                                         rg_dList<Gate>& theGateStack, 
                                                                         Gate& currGate)
{
    rg_INT          numBallsTangentToSphere = currGate.getNumBallsTangentToSphere();
    BallGenerator** ballsTangentToSphere    = currGate.getBallsTangentToSphere();
    rg_Point3D*     ptsOnNormalizedSphere   = currGate.getNormalizedTangentPointsBetweenBallsAndSphere();
    
    

    Triangulation3D triangulationOfPtsOnSphere;
    triangulationOfPtsOnSphere.tetrahedralizePointsOnSphere(GIVEN_FACE_ON_CONVEXHULL, 
                                       numBallsTangentToSphere, ptsOnNormalizedSphere );
    
    Sphere     tangentSphereAtEndVertex = currGate.getEmptyTangentSphereAtEndVertex();
    
    rg_INT     numTetrahedra             = triangulationOfPtsOnSphere.getNumOfTetrahedra();
    VDVertex** vertexMappedByTetrahedron = new VDVertex*[ numTetrahedra ];
    rg_INT i=0;
	for ( i=0; i<numTetrahedra; i++ )  {
        vertexMappedByTetrahedron[i] = makeVoronoiVertex( rg_FALSE, tangentSphereAtEndVertex );
    }

    if ( ON_DEBUG )  {
        foutEVDS << "*** Degeneracy: ";
        for ( i=0; i<numBallsTangentToSphere; i++)
            foutEVDS << ballsTangentToSphere[i]->getID() << "\t";
        foutEVDS << endl;
    }

    //  determine the end vertex of curr gate.
    T3DFace faceForFirstEndVertex;
    triangulationOfPtsOnSphere.locateFaceByVertexID(faceForFirstEndVertex, 0, 1, 2);
    
    T3DTetrahedron* currTetra = faceForFirstEndVertex.getTetrahedron();
    rg_INT          currID    = currTetra->getID();
    rg_INT mateVertexPos=0;
	for ( mateVertexPos=0; mateVertexPos<4; mateVertexPos++ )  {
        T3DTetrahedron* neighbor = currTetra->getNeighbor(mateVertexPos);
               
        T3DVertex* vertexOnFace[3] = {rg_NULL, rg_NULL, rg_NULL};
        currTetra->getVerticesOnFace(mateVertexPos, vertexOnFace);
        BallGenerator* ballToDefineFace[3] = { ballsTangentToSphere[ vertexOnFace[0]->getID() ],
                                               ballsTangentToSphere[ vertexOnFace[1]->getID() ],
                                               ballsTangentToSphere[ vertexOnFace[2]->getID() ] };

        T3DVertex*     mateVertex = currTetra->getVertex(mateVertexPos);
        BallGenerator* startBall  = ballsTangentToSphere[ mateVertex->getID() ];

        if ( mateVertexPos == faceForFirstEndVertex.getMateVertexPos() )  {
            //  make Voronoi edge between this tetra and its neighbor.  
            ballToDefineFace[0] = currGate.getGateBallCCW1();
            ballToDefineFace[1] = currGate.getGateBallCCW2();
            ballToDefineFace[2] = currGate.getGateBallCCW3();
            //updateLocalTopologyByEdge(thisWorld, currGate.getEdgeType(), currGate.getStartVertex(), 
            //                          vertexMappedByTetrahedron[currTetra->getID()], ballToDefineFace );
            updateLocalTopologyByEdge(thisWorld, currGate, vertexMappedByTetrahedron[currTetra->getID()]);

            if ( ON_DEBUG )
            {
                BallGenerator* endBall = ballsTangentToSphere[ mateVertex->getID() ];
                
                VDVertex* ptrEndVertex = vertexMappedByTetrahedron[currTetra->getID()];
                
                foutEVDS << currGate.getStartBall()->getID()   << "\t"
                         << ballToDefineFace[0]->getID()           << "\t"
                         << ballToDefineFace[1]->getID()           << "\t"
                         << ballToDefineFace[2]->getID()           << "\t"
                         << endBall->getID()               << "\t";

                rg_FLAG edgeType = currGate.getEdgeType();
                if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
                    foutEVDS << "Para" << "\t";
                else if ( edgeType == LINEAR_EDGE )
                    foutEVDS << "Line" << "\t";
                else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
                    foutEVDS << "Elli" << "\t";
                else 
                    foutEVDS << "Non" << "\t";

                foutEVDS << "DV"                           << "\t"
                         << currGate.getStartVertex()->getID() << "\t"
                         << ptrEndVertex->getID()          << "\t"
                         << ptrEndVertex->getPoint().getX() << "\t"
                         << ptrEndVertex->getPoint().getY() << "\t"
                         << ptrEndVertex->getPoint().getZ() << "\t"
                         << ptrEndVertex->getRadiusOfTangentSphere() << endl;
            }
        }
        else {
            if ( neighbor == rg_NULL ) {
                //  add gate
                Gate* ptrGate = theGateStack.addHead( Gate( ballToDefineFace[0], ballToDefineFace[1], ballToDefineFace[2], 
                                                            startBall, vertexMappedByTetrahedron[currID] ) );
                ptrGate->setNodePosInStack( theGateStack.getTail() );
                vertexMappedByTetrahedron[currID]->setGate( ptrGate );

                rg_INT numExcludeBalls = numBallsTangentToSphere - 3;
                BallGenerator** excludedBalls = new BallGenerator*[numExcludeBalls];

                rg_INT j=0;
                for ( rg_INT i=0; i<numBallsTangentToSphere; i++)  {
                    if (    ballsTangentToSphere[i] != ballToDefineFace[0]
                         && ballsTangentToSphere[i] != ballToDefineFace[1]
                         && ballsTangentToSphere[i] != ballToDefineFace[2] )
                         excludedBalls[j++] = ballsTangentToSphere[i];
                }
                ptrGate->setExcludedBalls( numExcludeBalls, excludedBalls );
                delete [] excludedBalls;
            }
        }
    }
    currTetra->setCheck(rg_TRUE);

    
    rg_dList<T3DTetrahedron>* tetrahedra = triangulationOfPtsOnSphere.getTetrahedra();
    tetrahedra->reset4Loop();
    while ( tetrahedra->setNext4Loop() )  {
        currTetra = tetrahedra->getpEntity();
        currID    = currTetra->getID();

        if ( currTetra->isChecked() )
            continue;

        for ( mateVertexPos=0; mateVertexPos<4; mateVertexPos++ )  {
            T3DTetrahedron* neighbor = currTetra->getNeighbor(mateVertexPos);
            
            
            T3DVertex* vertexOnFace[3] = {rg_NULL, rg_NULL, rg_NULL};
            currTetra->getVerticesOnFace(mateVertexPos, vertexOnFace);
            BallGenerator* ballToDefineFace[3] = { ballsTangentToSphere[ vertexOnFace[0]->getID() ],
                                                   ballsTangentToSphere[ vertexOnFace[1]->getID() ],
                                                   ballsTangentToSphere[ vertexOnFace[2]->getID() ] };

            T3DVertex*     mateVertex = currTetra->getVertex(mateVertexPos);
            BallGenerator* startBall  = ballsTangentToSphere[ mateVertex->getID() ];

            if ( neighbor == rg_NULL ) {
                //  add gate
                Gate* ptrGate = theGateStack.addHead( Gate( ballToDefineFace[0], ballToDefineFace[1], ballToDefineFace[2], 
                                                            startBall, vertexMappedByTetrahedron[currID] ) );
                ptrGate->setNodePosInStack( theGateStack.getTail() );
                vertexMappedByTetrahedron[currID]->setGate( ptrGate );
                
                rg_INT numExcludeBalls = numBallsTangentToSphere - 3;
                BallGenerator** excludedBalls = new BallGenerator*[numExcludeBalls];

                rg_INT j=0;
                for ( rg_INT i=0; i<numBallsTangentToSphere; i++)  {
                    if (    ballsTangentToSphere[i] != ballToDefineFace[0]
                         && ballsTangentToSphere[i] != ballToDefineFace[1]
                         && ballsTangentToSphere[i] != ballToDefineFace[2] )
                         excludedBalls[j++] = ballsTangentToSphere[i];
                }
                ptrGate->setExcludedBalls( numExcludeBalls, excludedBalls );
                delete [] excludedBalls;
            }
            else  {
                if ( neighbor->isChecked() == rg_FALSE )
                    continue;
                
                //  make Voronoi edge between this tetra and its neighbor.  
                updateLocalTopologyByEdge(thisWorld, EDGE_BY_VERTEX_DEGENERACY, vertexMappedByTetrahedron[currID], 
                                          vertexMappedByTetrahedron[neighbor->getID()], ballToDefineFace );

                if ( ON_DEBUG )
                {
                    BallGenerator* endBall   = ballsTangentToSphere[ mateVertex->getID() ];
                    BallGenerator* startBall = ballsTangentToSphere[ neighbor->getMateVertex(currTetra)->getID() ];
                
                    VDVertex* ptrEndVertex = vertexMappedByTetrahedron[currTetra->getID()];
                
                    foutEVDS << startBall->getID()   << "\t"
                             << ballToDefineFace[0]->getID()           << "\t"
                             << ballToDefineFace[1]->getID()           << "\t"
                             << ballToDefineFace[2]->getID()           << "\t"
                             << endBall->getID()               << "\t";

                    rg_FLAG edgeType = currGate.getEdgeType();
                    if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
                        foutEVDS << "Para" << "\t";
                    else if ( edgeType == LINEAR_EDGE )
                        foutEVDS << "Line" << "\t";
                    else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
                        foutEVDS << "Elli" << "\t";
                    else 
                        foutEVDS << "Non" << "\t";

                    foutEVDS << "DV"                           << "\t"
                             << currGate.getStartVertex()->getID() << "\t"
                             << ptrEndVertex->getID()          << "\t"
                             << ptrEndVertex->getPoint().getX() << "\t"
                             << ptrEndVertex->getPoint().getY() << "\t"
                             << ptrEndVertex->getPoint().getZ() << "\t"
                             << ptrEndVertex->getRadiusOfTangentSphere() << endl;
                }
            }
        }
        
        currTetra->setCheck(rg_TRUE);
    }

    delete [] ballsTangentToSphere;
    delete [] ptsOnNormalizedSphere;
    delete [] vertexMappedByTetrahedron;
}



void rg_SphereSetVoronoiDiagram::updateLocalTopologyByEdge( VDLoop* thisWorld, const rg_FLAG& edgeType, 
                                                            VDVertex* startVertex, VDVertex* endVertex, 
                                                            BallGenerator** ballToDefineEdge )
{
    ///////////////////////////////////////////////////////////////
    //
    //  Create Voronoi edge.
    //
    VDEdge* ptrEdge = m_GlobalEdgeList.addTail( VDEdge( m_GlobalEdgeList.getSize(), startVertex, endVertex ) );
    ptrEdge->setEdgeType( edgeType );
    startVertex->setIncidentEdge( ptrEdge );
    endVertex->setIncidentEdge( ptrEdge );



    ///////////////////////////////////////////////////////////////
    //
    //  Create three partial edges.
    //
    VDPartialEdge* ptrPartialEdge[3] = {rg_NULL, rg_NULL, rg_NULL};
    rg_INT i=0;
	for ( i=0; i<3; i++ ) {
        ptrPartialEdge[i] = m_GlobalPartialedgeList.addTail( VDPartialEdge( m_GlobalPartialedgeList.getSize(), ptrEdge ) );

        ptrPartialEdge[i]->setNextPartialEdgeInLoop( ptrPartialEdge[i] );
        ptrPartialEdge[i]->setPreviousPartialEdgeInLoop( ptrPartialEdge[i] );
    }
    //  connect partial edges in Radial Cycle.
    ptrPartialEdge[0]->setNextPartialEdgeInRadialCycle( ptrPartialEdge[1] );
    ptrPartialEdge[1]->setNextPartialEdgeInRadialCycle( ptrPartialEdge[2] );
    ptrPartialEdge[2]->setNextPartialEdgeInRadialCycle( ptrPartialEdge[0] );

    ptrEdge->setPartialEdge( ptrPartialEdge[0] );



    ///////////////////////////////////////////////////////////////
    //
    //  Create Voronoi face and loop.
    //
    VDCell* gateCell[4] = { ballToDefineEdge[0]->getCell(), ballToDefineEdge[1]->getCell(), 
                            ballToDefineEdge[2]->getCell(), ballToDefineEdge[0]->getCell() };

    //  partial edge[0] <- gate ball 0 & gate ball 1
    //  partial edge[1] <- gate ball 1 & gate ball 2
    //  partial edge[2] <- gate ball 2 & gate ball 0

    for ( i=0; i<3; i++ )  {
        VDFace* ptrFace = gateCell[i]->findFaceToShareWith( gateCell[i+1] );
        VDLoop* ptrLoop = rg_NULL;

        if ( ptrFace == rg_NULL )  {
            ptrFace = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), gateCell[i+1], gateCell[i] ) );
            ptrLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrFace, rg_TRUE ) );

            //  define connection between VDCell and VDFace.
            gateCell[i]->addBoundingFace( ptrFace );
            gateCell[i+1]->addBoundingFace( ptrFace );

            //  define connection between VDFace and VDLoop.
            ptrFace->addLoop( ptrLoop );

            //  define connection between VDLoop and VDPartialEdge.
            ptrPartialEdge[i]->isRightOrientationInLoop(rg_TRUE);
        }
        else  {
            if ( thisWorld != rg_NULL && ptrFace == thisWorld->getFace() )  
                ptrLoop = thisWorld;
            else  
                ptrLoop = ptrFace->getLoop( 0 );

            //  define connection between VDLoop and VDPartialEdge.
            if (    ptrFace->getLeftCell() == gateCell[i] 
                 && ptrFace->getRightCell() == gateCell[i+1] )
                ptrPartialEdge[i]->isRightOrientationInLoop(rg_TRUE);
            else
                ptrPartialEdge[i]->isRightOrientationInLoop(rg_FALSE);
        }

        ptrLoop->addPartialEdgeToLoop( ptrPartialEdge[i] );
    }

}



void rg_SphereSetVoronoiDiagram::updateLocalTopologyByEdge(VDLoop* thisWorld, const Gate& currGate, VDVertex* endVertex)
{
    ///////////////////////////////////////////////////////////////
    //
    //  Create Voronoi edge.
    //
    VDVertex* startVertex = currGate.getStartVertex();
    VDEdge*   ptrEdge     = m_GlobalEdgeList.addTail( VDEdge( m_GlobalEdgeList.getSize(), startVertex, endVertex ) );
    ptrEdge->setEdgeType( currGate.getEdgeType() );
    ptrEdge->setMinMaxTangentSphereByBalls( currGate.getMinTangentSphere(0), currGate.getMinTangentSphere(1) );
    ptrEdge->setEdgePlane( currGate.getEdgePlane() );
    ptrEdge->setPointForAngleDistance( currGate.getAxisPoint() );

    startVertex->setIncidentEdge( ptrEdge );
    endVertex->setIncidentEdge( ptrEdge );


    BallGenerator* ballToDefineEdge[3]   = { currGate.getGateBallCCW1(), currGate.getGateBallCCW2(), currGate.getGateBallCCW3() };

    ///////////////////////////////////////////////////////////////
    //
    //  Create three partial edges.
    //
    VDPartialEdge* ptrPartialEdge[3] = {rg_NULL, rg_NULL, rg_NULL};
    rg_INT i=0;
	for ( i=0; i<3; i++ ) {
        ptrPartialEdge[i] = m_GlobalPartialedgeList.addTail( VDPartialEdge( m_GlobalPartialedgeList.getSize(), ptrEdge ) );

        ptrPartialEdge[i]->setNextPartialEdgeInLoop( ptrPartialEdge[i] );
        ptrPartialEdge[i]->setPreviousPartialEdgeInLoop( ptrPartialEdge[i] );
    }
    //  connect partial edges in Radial Cycle.
    ptrPartialEdge[0]->setNextPartialEdgeInRadialCycle( ptrPartialEdge[1] );
    ptrPartialEdge[1]->setNextPartialEdgeInRadialCycle( ptrPartialEdge[2] );
    ptrPartialEdge[2]->setNextPartialEdgeInRadialCycle( ptrPartialEdge[0] );

    ptrEdge->setPartialEdge( ptrPartialEdge[0] );



    ///////////////////////////////////////////////////////////////
    //
    //  Create Voronoi face and loop.
    //
    VDCell* gateCell[4] = { ballToDefineEdge[0]->getCell(), ballToDefineEdge[1]->getCell(), 
                            ballToDefineEdge[2]->getCell(), ballToDefineEdge[0]->getCell() };

    //  partial edge[0] <- gate ball 0 & gate ball 1
    //  partial edge[1] <- gate ball 1 & gate ball 2
    //  partial edge[2] <- gate ball 2 & gate ball 0

    for ( i=0; i<3; i++ )  {
        VDFace* ptrFace = gateCell[i]->findFaceToShareWith( gateCell[i+1] );
        VDLoop* ptrLoop = rg_NULL;

        if ( ptrFace == rg_NULL )  {
            ptrFace = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), gateCell[i+1], gateCell[i] ) );
            ptrLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrFace, rg_TRUE ) );

            //  define connection between VDCell and VDFace.
            gateCell[i]->addBoundingFace( ptrFace );
            gateCell[i+1]->addBoundingFace( ptrFace );

            //  define connection between VDFace and VDLoop.
            ptrFace->addLoop( ptrLoop );

            //  define connection between VDLoop and VDPartialEdge.
            ptrPartialEdge[i]->isRightOrientationInLoop(rg_TRUE);
        }
        else  {
            if ( thisWorld != rg_NULL && ptrFace == thisWorld->getFace() )  
                ptrLoop = thisWorld;
            else  
                ptrLoop = ptrFace->getLoop( 0 );

            //  define connection between VDLoop and VDPartialEdge.
            if (    ptrFace->getLeftCell() == gateCell[i] 
                 && ptrFace->getRightCell() == gateCell[i+1] )
                ptrPartialEdge[i]->isRightOrientationInLoop(rg_TRUE);
            else
                ptrPartialEdge[i]->isRightOrientationInLoop(rg_FALSE);
        }

        ptrLoop->addPartialEdgeToLoop( ptrPartialEdge[i] );
    }
}


void rg_SphereSetVoronoiDiagram::constructSmallWorldRobustly()
{
    rg_dList<BallGenerator*> isolatedBalls;
    collectIsolatedBalls( isolatedBalls );

    rg_dList<BallGenerator*> isolatedBallsInDisconnectedBigWorlds;

    BallGenerator* currBall = rg_NULL;
    while ( isolatedBalls.getSize() > 0 )  { 
        currBall = isolatedBalls.popFront();

        if ( currBall->getCell()->getNumOfBoundingFaces() != 0 )
            continue;

        //  1. find big brothers.        
        BallGenerator* bigBrother[2]     = { rg_NULL, rg_NULL };
        VDFace*        faceByBigBrothers = findBigBrothers( currBall, bigBrother );

        if ( faceByBigBrothers == rg_NULL )  {
            isolatedBallsInDisconnectedBigWorlds.add( currBall );
            continue;
        }


        //  2. find initial Voronoi vertex in the small world.
        rg_INT           initialSizeOfVIDIC = 11;
        VIDIC            theVIDIC( initialSizeOfVIDIC );
        rg_dList< Gate > theGateStack;
        VDVertex*        pVertex = rg_NULL;
        rg_FLAG          bValidSmallWorld = findInitialVertexInSmallWorld(currBall, bigBrother, pVertex, 
                                                                          theVIDIC, theGateStack );


        if ( bValidSmallWorld )  {
            if ( pVertex != rg_NULL )  {
                VDLoop* smallWorld = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), 
                                                                       faceByBigBrothers, rg_FALSE ) );
                faceByBigBrothers->addLoop( smallWorld );
                
                constructThisWorldRobustly( smallWorld, theVIDIC, theGateStack );
            }
            else  {
                constructSmallWorldWithoutVertices(currBall, bigBrother, faceByBigBrothers );
            }
        }
        else  {
                isolatedBalls.addTail( currBall );
        }

        if ( m_status == VERTEX_DEGREE_VIOLATION )
            return;
    }

    //  we have to deal with the disconnected big worlds.
}



void rg_SphereSetVoronoiDiagram::updateGlobalTopology()
{
    defineOrderOfPartialEdgesInAllLoops();
    setTopologyOnInfinity();

    rg_dList< VDCell >* globalCellList = getGlobalCellList();
    VDCell* infinityCell               = globalCellList->getLastpEntity();
    
    rg_dList<VDFace*>* boundingFaces = infinityCell->getBoungindFaces();
    VDFace* currFace = rg_NULL;
    boundingFaces->reset4Loop();
    while ( boundingFaces->setNext4Loop() )  {
        currFace = boundingFaces->getEntity();

        if ( currFace->getRightCell() == infinityCell ) 
            currFace->getLeftCell()->isBounded( rg_FALSE );
        else
            currFace->getRightCell()->isBounded( rg_FALSE );
    }



    VDEdge* currEdge = rg_NULL;
    m_GlobalEdgeList.reset4Loop();
    while ( m_GlobalEdgeList.setNext4Loop() )  {
        currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isOnInfinity() == rg_TRUE )  {
            rg_dList<VDFace*> faceList;
            currEdge->inquireIncidentFaces(faceList);

            faceList.reset4Loop();
            while ( faceList.setNext4Loop() )  {
                faceList.getEntity()->isBounded(rg_FALSE);
            }
        }
    }
}

//
//  End of Robust Edge-tracing Algorithm
//
///////////////////////////////////////////////////////////////////////////////

rg_INT rg_SphereSetVoronoiDiagram::constructSphereSetVoronoiDiagramByEdgeTracing()
{
   m_status = rg_TRUE;
    rg_INT initialSizeOfVIDIC = 101;

    VIDIC            theVIDIC( initialSizeOfVIDIC );
    rg_dList< Gate > theGateStack;


    //  0. find initial voronoi vertex.
    findInitialVertexForEdgeTracing( theVIDIC, theGateStack );


    while ( theGateStack.getSize() > 0 )
    {
        // 1.0 pop a gate;
        Gate currGate = theGateStack.popFront();


        findEndVertexBasedOnAngularDistance( currGate ); 


        // 1.2 EXPLORE GATE
        VDVertex* ptrEndVertex = exploreEdgeToBeMadeByCurrentGate( currGate, theVIDIC, theGateStack );

        if ( m_status != rg_TRUE )  {
            break;
        }

        setLocalTopologyByEndVertex( currGate, ptrEndVertex ); 

    }

    if ( m_status == rg_TRUE )  {
        setTopology();

        return m_status;
    }
    else 
    {
        return m_status;
    }
}



void rg_SphereSetVoronoiDiagram::shrinkBallsByMinimumRaidus()
{
    BallGenerator* currBall = rg_NULL;
	m_GlobalGeneratorList.reset4Loop();

	while( m_GlobalGeneratorList.setNext4Loop() )  {
		currBall = m_GlobalGeneratorList.getpEntity();

        rg_REAL radius = currBall->getRadius();
        currBall->setRadius( radius - m_minRadiusOfBalls );
    }
}



void rg_SphereSetVoronoiDiagram::recoverBallsByMinimumRaidus()
{
    BallGenerator* currBall = rg_NULL;
	m_GlobalGeneratorList.reset4Loop();

	while( m_GlobalGeneratorList.setNext4Loop() )  {
		currBall = m_GlobalGeneratorList.getpEntity();

        rg_REAL radius = currBall->getRadius();
        currBall->setRadius( radius + m_minRadiusOfBalls );
    }

    VDVertex* currVertex = rg_NULL;
	m_GlobalVertexList.reset4Loop();

	while( m_GlobalVertexList.setNext4Loop() )  {
		currVertex = m_GlobalVertexList.getpEntity();

        rg_REAL radius = currVertex->getRadiusOfTangentSphere();
        currVertex->setRadiusOfTangentSphere( radius - m_minRadiusOfBalls );
    }
}




void rg_SphereSetVoronoiDiagram::findInitialVertexForEdgeTracing( VIDIC&          theVIDIC,
                                                                  rg_dList<Gate>& theGateStack )
{

    //  1. set the configuration of virtual initial vertex
    rg_Point3D virtualVertexCoord(1.0, 1.0, 1.0);
    rg_REAL    virtualVertexRadius  = 0.73205080756887719;    

    rg_Point3D centerOfVirtualBall[4];
    centerOfVirtualBall[0].setPoint(2.0, 0.0, 0.0);
    centerOfVirtualBall[1].setPoint(0.0, 2.0, 0.0);
    centerOfVirtualBall[2].setPoint(0.0, 0.0, 2.0);
    centerOfVirtualBall[3].setPoint(0.0, 0.0, 0.0);

    rg_Point3D oldMinPointOfBoundingBoxForBalls = m_minPointOfBoundingBoxForBalls;

    rg_Point3D translation(-5.0, -5.0, -5.0);
    translation = m_minPointOfBoundingBoxForBalls + translation;
    m_minPointOfBoundingBoxForBalls = translation;
    
    virtualVertexCoord = virtualVertexCoord + translation;


    BallGenerator* ptrVirtualBall[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
    VDCell*        ptrVirtualCell[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
    rg_INT i=0;
	for (i=0; i<4; i++)  {
        centerOfVirtualBall[i] = centerOfVirtualBall[i] + translation;

        ptrVirtualBall[i] = m_GlobalGeneratorList.addHead( BallGenerator(centerOfVirtualBall[i], 1.0) );
        ptrVirtualCell[i] = m_GlobalCellList.addHead( VDCell( -2-i ) );
        connectCellAndGenerator( ptrVirtualCell[i], ptrVirtualBall[i] );
    }

    //  2. find valid initial vertex
    setBucket( m_sizeOfBucketElement );
    VIDIC            tempVIDIC( 31 );
    rg_dList< Gate > tempGateStack;

    defineInitialVoronoiVertex( ptrVirtualBall, Sphere(virtualVertexCoord, virtualVertexRadius), tempVIDIC, tempGateStack);

    Sphere initialTangentSphere;
    Gate   initialGates[NUM_DEFINING_CELLS_OF_ONE_VERTEX];

    while ( tempGateStack.getSize() > 0 )
    {
        // 1.0 pop a gate;
        Gate currGate = tempGateStack.popBack();


        // 1.1 GET THE END VERTEX
        findEndVertexBasedOnAngularDistanceInAcceleration( currGate ); 
        BallGenerator* endBall = currGate.getEndBall();
        Sphere  tangentSphereByEndBall( currGate.getEmptyTangentSphereAtEndVertex() );

        if ( endBall != rg_NULL )
        {
            if (    currGate.getGateBallCCW1()->getCell()->getID() > -2 && currGate.getGateBallCCW2()->getCell()->getID() > -2
                 && currGate.getGateBallCCW3()->getCell()->getID() > -2 && endBall->getID() > -2 )
            {
                BallGenerator* gateBall[3] = { currGate.getGateBallCCW1(), currGate.getGateBallCCW2(), currGate.getGateBallCCW3()};

                initialTangentSphere = tangentSphereByEndBall;
                initialGates[0].setGate( gateBall[2], gateBall[1], gateBall[0], endBall, rg_NULL );
                initialGates[1].setGate( endBall, gateBall[0], gateBall[1], gateBall[2], rg_NULL );
                initialGates[2].setGate( endBall, gateBall[1], gateBall[2], gateBall[0], rg_NULL );
                initialGates[3].setGate( endBall, gateBall[2], gateBall[0], gateBall[1], rg_NULL );

                break;
            }
        }

        // 1.2 EXPLORE GATE
        //VDVertex* ptrEndVertex = exploreEdgeToBeMadeByCurrentGate( currGate, endBall, 
        //                                                           tangentSphereByEndBall, 
        //                                                           tempVIDIC, tempGateStack, edgeType );
        VDVertex* ptrEndVertex = exploreEdgeToBeMadeByCurrentGate( currGate, tempVIDIC, tempGateStack );

        makeVoronoiEdge( currGate.getStartVertex(), ptrEndVertex );
        //setLocalTopologyByEndVertex( currentGate, ptrEndVertex ); 
    }
    m_bucket.removeBucketTable();

    //  3. remove virtual initial vertex and its topology
    m_minPointOfBoundingBoxForBalls = oldMinPointOfBoundingBoxForBalls;

    for ( i=0; i<4; i++)  {
        m_GlobalCellList.popFront();
        m_GlobalGeneratorList.popFront();
    }

    VDCell* currCell = rg_NULL;
	m_GlobalCellList.reset4Loop();
	while( m_GlobalCellList.setNext4Loop() )
	{
		currCell = m_GlobalCellList.getpEntity();
        currCell->getBoungindFaces()->removeAll();
    }
    m_GlobalEdgeList.removeAll();
    m_GlobalVertexList.removeAll();


    //  4. set valid initial vertex.
    //  insert initial voronoi vertex into GVL(Global Vertex List).
    VDVertex* ptrVertex = makeVoronoiVertex( rg_FALSE, initialTangentSphere );

    //  make and insert 4 gates of initial voronoi vertex into Gate-Stack.
    Gate* ptrGate = rg_NULL;        
    for ( i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
    {
        ptrGate = theGateStack.addTail( initialGates[i] );

        ptrGate->setStartVertex( ptrVertex );
        ptrGate->setNodePosInStack( theGateStack.getTail() );

        ptrVertex->setGate( i, ptrGate );
    }   

    //  insert vertex configuration of initial voronoi vertex into VIDIC.
    theVIDIC.insertVIDICItem( VIDICItem( initialGates[0].getGateBallCCW1(),
                                         initialGates[0].getGateBallCCW2(),
                                         initialGates[0].getGateBallCCW3(),
                                         initialGates[0].getStartBall(),
                                         ptrVertex) 
                             );
}





void rg_SphereSetVoronoiDiagram::defineInitialVoronoiVertex(       BallGenerator** surroundingBalls,
                                                             const Sphere&         tangentSphere,
                                                                   VIDIC&          theVIDIC,
                                                                   rg_dList<Gate>& theGateStack )        
{
    //  insert initial voronoi vertex into GVL(Global Vertex List).
    VDVertex* ptrVertex = makeVoronoiVertex( rg_FALSE, tangentSphere );


    //  make and insert 4 gates of initial voronoi vertex into Gate-Stack.
    Gate currentGates[NUM_DEFINING_CELLS_OF_ONE_VERTEX] 
             = { Gate( surroundingBalls[0], surroundingBalls[1], surroundingBalls[2], surroundingBalls[3], ptrVertex ),
                 Gate( surroundingBalls[0], surroundingBalls[2], surroundingBalls[3], surroundingBalls[1], ptrVertex ),
                 Gate( surroundingBalls[0], surroundingBalls[3], surroundingBalls[1], surroundingBalls[2], ptrVertex ),
                 Gate( surroundingBalls[2], surroundingBalls[1], surroundingBalls[3], surroundingBalls[0], ptrVertex ) };

    Gate* initialGate[4];
    rg_INT i=0;
	for (i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
    {
        Gate* ptrGate = theGateStack.addTail( currentGates[i] );
        ptrGate->setNodePosInStack( theGateStack.getTail() );

        ptrVertex->setGate( i, ptrGate );

        initialGate[i] = ptrGate;
    }   

    //  insert vertex configuration of initial voronoi vertex into VIDIC.
    theVIDIC.insertVIDICItem( VIDICItem(surroundingBalls, ptrVertex) );
}

VDVertex* rg_SphereSetVoronoiDiagram::makeVoronoiVertex( const rg_FLAG& onInfinity, 
                                                        const Sphere&  tangentSphere)
{
    VDVertex  vertex( m_GlobalVertexList.getSize(), 
                      onInfinity, 
                      tangentSphere.getCenter(), 
                      tangentSphere.getRadius() );

    return m_GlobalVertexList.addTail( vertex );
}




rg_FLAG rg_SphereSetVoronoiDiagram::findEndVertexBasedOnAngularDistance( Gate& gate )
{
    rg_Point3D normalOfCenterPlane;
    rg_Point3D normalOfEdgePlane;
    rg_Point3D axisPoint;
    Sphere     tangentSphereInCP[2];

    rg_FLAG    edgeType = computeEdgePlaneForAngleDistanceAccordingToEdgeType( 
                                 gate, normalOfCenterPlane, normalOfEdgePlane, axisPoint, tangentSphereInCP);
    gate.setEdgeType( edgeType );


    //  This gate defines a Voronoi edge.
    BallGenerator* ptrGateBall[3]  = { gate.getGateBallCCW1(), gate.getGateBallCCW2(), gate.getGateBallCCW3() };
    Sphere         gateBall[3]     = { ptrGateBall[0]->getBall(), ptrGateBall[1]->getBall(), ptrGateBall[2]->getBall() };

    BallGenerator* ptrStartBall = gate.getStartBall();
    rg_Point3D     startVertex  = gate.getStartVertex()->getPoint();
    rg_Point3D     vecPaVs      = startVertex - axisPoint;

    rg_Point3D     normalStartPlane = normalOfEdgePlane.crossProduct( vecPaVs );
    rg_REAL        coffStartPlane   = -normalStartPlane.innerProduct( startVertex );

    ///////////////////////////////////////////////////////////////////////////
    //
    //  candidateEndBall.first   : pointer to endBall
    //  candidateEndBall.second  : tangent sphere by endBall and gateBalls
    pair<BallGenerator*, Sphere> candidateEndBall;
    candidateEndBall.first   = rg_NULL;
    rg_REAL minAngularDistance = DBL_MAX;

    Sphere candidateTangentSpheres[2];
    rg_INT numOfTangentSpheres = 0;



    if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
    {
        rg_Point3D     vecCtPa  = axisPoint - tangentSphereInCP[0].getCenter();
        rg_REAL        maxValidAngle = angleCCW( normalOfEdgePlane, vecPaVs, vecCtPa);
        
        BallGenerator* currBall = rg_NULL;
        m_GlobalGeneratorList.reset4Loop();
	    while( m_GlobalGeneratorList.setNext4Loop() )
	    {
	        currBall = m_GlobalGeneratorList.getpEntity();

            //  exclude three gate balls to find end vertex.
            if ( currBall == ptrGateBall[0] || currBall == ptrGateBall[1] || currBall == ptrGateBall[2] )
                continue;

            numOfTangentSpheres = computeSphereTangentTo4SpheresOutside(
                                         gateBall[0], gateBall[1], gateBall[2], currBall->getBall(),
                                         candidateTangentSpheres);

            rg_INT i=0;
			for ( i=0; i<numOfTangentSpheres; i++ )  
            {
                if ( currBall == ptrStartBall )  {
                    if ( numOfTangentSpheres == 1 )
                        continue;
                    else
                        if ( candidateTangentSpheres[i].getCenter() == startVertex )
                            continue;
                }        
        
                rg_Point3D vecViPa = candidateTangentSpheres[i].getCenter() - axisPoint;
                rg_REAL    angularDistance = angleCCW( normalOfEdgePlane, vecPaVs, vecViPa);

                rg_REAL distance = normalStartPlane.innerProduct( candidateTangentSpheres[i].getCenter() )
                                   + coffStartPlane;

                if ( angularDistance < maxValidAngle && angularDistance < minAngularDistance)  {
                    candidateEndBall.first  = currBall;
                    candidateEndBall.second = candidateTangentSpheres[i];
                    minAngularDistance        = angularDistance;                
                }
            }
        }       
    }
    else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
    {
        BallGenerator* currBall = rg_NULL;
        m_GlobalGeneratorList.reset4Loop();
	    while( m_GlobalGeneratorList.setNext4Loop() )
	    {
	        currBall = m_GlobalGeneratorList.getpEntity();

            //  exclude three gate balls to find end vertex.
            if ( currBall == ptrGateBall[0] || currBall == ptrGateBall[1] || currBall == ptrGateBall[2] )
                continue;

            numOfTangentSpheres = computeSphereTangentTo4SpheresOutside(
                                         gateBall[0], gateBall[1], gateBall[2], currBall->getBall(),
                                         candidateTangentSpheres);

            rg_INT i=0;
			for ( i=0; i<numOfTangentSpheres; i++ )  
            {
                if ( currBall == ptrStartBall )  {
                    if ( numOfTangentSpheres == 1 )
                        continue;
                    else
                        if ( candidateTangentSpheres[i].getCenter() == startVertex )
                            continue;
                }        
        
                rg_Point3D vecViPa = candidateTangentSpheres[i].getCenter() - axisPoint;
                rg_REAL    angularDistance = angleCCW( normalOfEdgePlane, vecPaVs, vecViPa);

                if ( angularDistance < minAngularDistance )  {
                    candidateEndBall.first  = currBall;
                    candidateEndBall.second = candidateTangentSpheres[i];
                    minAngularDistance        = angularDistance;                
                }
            }
        }        
    }
    else if ( edgeType == LINEAR_EDGE )
    {
        BallGenerator* currBall = rg_NULL;
        m_GlobalGeneratorList.reset4Loop();
	    while( m_GlobalGeneratorList.setNext4Loop() )
	    {
	        currBall = m_GlobalGeneratorList.getpEntity();

            //  exclude three gate balls to find end vertex.
            if ( currBall == ptrGateBall[0] || currBall == ptrGateBall[1] || currBall == ptrGateBall[2] )
                continue;

            numOfTangentSpheres = computeSphereTangentTo4SpheresOutside(
                                         gateBall[0], gateBall[1], gateBall[2], currBall->getBall(),
                                         candidateTangentSpheres);

            rg_INT i=0;
			for ( i=0; i<numOfTangentSpheres; i++ )  
            {
                if ( currBall == ptrStartBall )  {
                    if ( numOfTangentSpheres == 1 )
                        continue;
                    else
                        if ( candidateTangentSpheres[i].getCenter() == startVertex )
                            continue;
                }        
        
                rg_Point3D vecViVs = candidateTangentSpheres[i].getCenter() - startVertex;
                
                rg_REAL    orientation = normalOfCenterPlane.innerProduct( vecViVs.getUnitVector() );
                if ( !rg_POS( orientation ) )
                    continue;

                rg_REAL    angularDistance = vecViVs.magnitude();

                if ( angularDistance < minAngularDistance )  {
                    candidateEndBall.first  = currBall;
                    candidateEndBall.second = candidateTangentSpheres[i];
                    minAngularDistance        = angularDistance;                
                }
            }
        }
    }
    else
    {}

    if ( candidateEndBall.first == rg_NULL )  {
        return rg_FALSE;
    }
    else  {
        gate.setEndBall(candidateEndBall.first);
        gate.setEmptyTangentSphereAtEndVertex(candidateEndBall.second);

        return rg_TRUE;
    }
}




rg_FLAG rg_SphereSetVoronoiDiagram::findEndVertexBasedOnAngularDistanceInAcceleration( Gate& gate )
{
    ///////////////////////////////////////////////////////////////////////////
    //
    //  prepare finding end vertex by computing edge plane.
    //
    start_in = clock();

    rg_FLAG    result = rg_FALSE;
    rg_Point3D normalOfCenterPlane;
    rg_Point3D normalOfEdgePlane;
    rg_Point3D axisPoint;
    Sphere     tangentSphereInCP[2];
    rg_FLAG    edgeType = computeEdgePlaneForAngleDistanceAccordingToEdgeType( 
                                 gate, normalOfCenterPlane, normalOfEdgePlane, axisPoint, tangentSphereInCP);
    gate.setEdgeType( edgeType );

    //  computingTime[2]:  preparing end vertex finding
    finish_in = clock();
    computingTime[2] += finish_in - start_in;


    ///////////////////////////////////////////////////////////////////////////
    //
    //  find end vertex.
    //    
    start_in = clock();
    if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )  {
        result =  findEndVertexOfParabolicAndHyperbolicEdge(gate, 
                                                            normalOfCenterPlane, 
                                                            normalOfEdgePlane,
                                                            axisPoint,
                                                            tangentSphereInCP[0]); 
        
        //  computingTime[3]:  finding end vertex of parabolic or hyperbolic edge.
        finish_in = clock();
        computingTime[3] += finish_in - start_in;
    }
    else if ( edgeType == LINEAR_EDGE )  {
        result = findEndVertexOfLinearEdgeInAcceleration( gate, 
                                                          normalOfCenterPlane, 
                                                          tangentSphereInCP[0]);
        
        //  computingTime[7]:  finding end vertex of linear edge.
        finish_in = clock();
        computingTime[7] += finish_in - start_in;
    }
    else  { // if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE ) 
        result = findEndVertexOfCircularAndEllipticEdgeInAcceleration(gate, 
                                                                      normalOfCenterPlane, 
                                                                      normalOfEdgePlane,
                                                                      axisPoint,
                                                                      tangentSphereInCP);

        //  computingTime[11]:  finding end vertex of circular or elliptic edge
        finish_in = clock();
        computingTime[11] += finish_in - start_in;
    }

    return result;
}





rg_FLAG rg_SphereSetVoronoiDiagram::findEndVertexOfParabolicAndHyperbolicEdge( 
                                                                Gate&       gate, 
                                                          const rg_Point3D& normalOfCenterPlane, 
                                                          const rg_Point3D& normalOfEdgePlane,
                                                          const rg_Point3D& axisPoint,
                                                          const Sphere&     minSphere)
{
    //  This gate defines a Voronoi edge.
    BallGenerator* gateBall[3]  = { gate.getGateBallCCW1(), gate.getGateBallCCW2(), gate.getGateBallCCW3() };
    Sphere         ball[3]      = { gateBall[0]->getBall(), gateBall[1]->getBall(), gateBall[2]->getBall() };
    BallGenerator* ptrStartBall = gate.getStartBall();
    
    rg_Point3D     startVertex  = gate.getStartVertex()->getPoint();
    rg_Point3D     vecPaVs      = startVertex - axisPoint;
    rg_Point3D     vecCtPa      = axisPoint - minSphere.getCenter();

    vecPaVs.normalize();    
    vecCtPa.normalize();
    rg_REAL maxValidAngle = angleCCWOfTwoUnitVectors( normalOfEdgePlane, vecPaVs, vecCtPa );



    ///////////////////////////////////////////////////////////////////////////
    //
    //  candidateEndBall.first   : pointer to endBall
    //  candidateEndBall.second  : tangent sphere by endBall and gateBalls
    pair<BallGenerator*, Sphere> candidateEndBall;
    candidateEndBall.first   = rg_NULL;
    rg_REAL minAngularDistance = DBL_MAX;


    rg_dList<BallGenerator*> checkedBallSet;

    rg_INT i=0;
	for ( i=0; i<3; i++)  {
        gateBall[i]->m_checkForBucket = rg_TRUE;
        checkedBallSet.addTail( gateBall[i] );
    }
    
    rg_REAL gatePlane[4];
    computePlaneTangentTo3SphereFromCCWOutside( ball, gatePlane );
    rg_Point3D normalOfGatePlane(gatePlane[0], gatePlane[1], gatePlane[2]);


    //////////////////////////////////////////////////////////////////////////
    //
    //  For Front Feasible region of parabolic and hyperbolic edge.
    //

    rg_INT numBucketElements = 0;

    //  construct sphere enclosing R1
    Sphere sphereForR1;

    rg_Point3D ptOnGatePlane[3];
    for ( i=0; i<3; i++)  {
        ptOnGatePlane[i] = ball[i].getCenter() + ( ball[i].getRadius()*normalOfGatePlane );    
    }

    rg_REAL    coeffCenterPlane = -(normalOfCenterPlane.innerProduct( minSphere.getCenter() ));
    rg_Point3D center           = computeCenterOfApproximatingSphereForR1( ptOnGatePlane,
                                                                           normalOfCenterPlane, coeffCenterPlane );
    sphereForR1.setCenter( center );
    sphereForR1.setRadius( center.distance( ptOnGatePlane[0] ) );


    rg_dList<BallGenerator*> ballSetForR1;
    numBucketElements = m_bucket.getCellsByMask( ballSetForR1, sphereForR1 );

    m_isFrontFeasibleRegion = rg_TRUE;
    findEmptyTangentSphereWithMinAngleInParabolicAndHyperbolicEdge( 
                                    candidateEndBall, minAngularDistance, ballSetForR1,
                                    ball, ptrStartBall, startVertex,  normalOfEdgePlane,
					                vecPaVs, axisPoint, maxValidAngle);
    m_isFrontFeasibleRegion = rg_FALSE;


    checkedBallSet.mergeTail( ballSetForR1 );


    rg_INT whereIsEndVertexFound = NOT_FOUND;


    if ( candidateEndBall.first != rg_NULL )  {
    
        // if candidate tangent sphere is intersected with R1.
        if ( candidateEndBall.second.isContainedIn( sphereForR1 ) == rg_FALSE )  {
            //  A CANDIDATE end ball is found.
            //  Therefore, we mask bucket for TS_R1 (tangent sphere from R1)

            rg_dList<BallGenerator*> ballSetForTS_C;
            numBucketElements += m_bucket.getCellsByMask( ballSetForTS_C, candidateEndBall.second );

            findEmptyTangentSphereWithMinAngleInParabolicAndHyperbolicEdge( 
                                            candidateEndBall, minAngularDistance, ballSetForTS_C,
                                            ball, ptrStartBall, startVertex,  normalOfEdgePlane,
					                        vecPaVs, axisPoint, maxValidAngle);

            checkedBallSet.mergeTail( ballSetForTS_C );

            whereIsEndVertexFound = FRONTREAR_FEASIBLE_REGION;
        }  
        else  {
            whereIsEndVertexFound = FRONT_FEASIBLE_REGION;
        }
    }
    else  {
        
        //////////////////////////////////////////////////////////////////////////
        //
        //  For Rear Feasible region of parabolic and hyperbolic edge.
        //
        rg_dList< NavigatingBucketCell > queueOfCells;

        rg_Point3D unitLength;

        initializePropagationInRearFeasibleRegionOfNonEllipticEdge( queueOfCells, unitLength, ball, gatePlane );

        rg_INT  normalSign = classifyNormalSign(gatePlane[0], gatePlane[1], gatePlane[2]);

        rg_dList< BallGenerator* > allCandidateBalls;


        NavigatingBucketCell* currBCell = rg_NULL;  
        rg_dNode< NavigatingBucketCell >* currNode = queueOfCells.getHead();

        do 
	    {
            currBCell = currNode->getpEntity();
            rg_dList< BallGenerator* > candidateBalls;

            if ( currBCell->m_status == ON_PLANE )   
                m_bucket.getBalls( candidateBalls, currBCell->m_index, gatePlane );
            else
                m_bucket.getBalls( candidateBalls, currBCell->m_index );


            if ( candidateBalls.getSize() > 0 )  {

                findEmptyTangentSphereWithMinAngleInParabolicAndHyperbolicEdge( 
                                                candidateEndBall, minAngularDistance, candidateBalls,
                                                ball, ptrStartBall, startVertex,  normalOfEdgePlane,
					                            vecPaVs, axisPoint, maxValidAngle);

                allCandidateBalls.mergeTail( candidateBalls );
                if ( candidateEndBall.first != rg_NULL )
                    break;
            }
       
            if ( currBCell->m_status == ON_PLANE )   
                propagateBucketCellsOnPlane( *currBCell, unitLength, queueOfCells );
            else
                propagateBucketCellsInPlusPlane( *currBCell, normalSign, queueOfCells );

            currNode = currNode->getNext();

        } while( currNode != queueOfCells.getHead() );

        checkedBallSet.mergeTail( allCandidateBalls );


        if ( candidateEndBall.first != rg_NULL )  {

            rg_dList<BallGenerator*> ballSetForR2;
            numBucketElements += m_bucket.getCellsByMask( ballSetForR2, candidateEndBall.second );

            findEmptyTangentSphereWithMinAngleInParabolicAndHyperbolicEdge( 
                                            candidateEndBall, minAngularDistance, ballSetForR2,
                                            ball, ptrStartBall, startVertex,  normalOfEdgePlane,
					                        vecPaVs, axisPoint, maxValidAngle);

            checkedBallSet.mergeTail( ballSetForR2 );
        }


        numBucketElements += queueOfCells.getSize();
        queueOfCells.reset4Loop();
	    while( queueOfCells.setNext4Loop() )  {
	        currBCell = queueOfCells.getpEntity();
            m_bucket.setMark(currBCell->m_index, rg_FALSE);
        }

        whereIsEndVertexFound = REAR_FEASIBLE_REGION;
    }

    
    if ( candidateEndBall.first == rg_NULL )  {
        /*
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //
        m_numSearchedBucketElements[0][3] += numBucketElements;

        m_numBallsForUnboundedEdge[0] += checkedBallSet.getSize();
        //
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //*/

        checkedBallSet.reset4Loop();
	    while( checkedBallSet.setNext4Loop() )  {
	        checkedBallSet.getEntity()->m_checkForBucket = rg_FALSE;
        }
    
        return rg_FALSE;
    }
    else  {
        gate.setEndBall( candidateEndBall.first );
        gate.setEmptyTangentSphereAtEndVertex( candidateEndBall.second );

        rg_INT numCheckedBalls = checkedBallSet.getSize();
        checkedBallSet.reset4Loop();
	    while( checkedBallSet.setNext4Loop() )  {
	        checkedBallSet.getEntity()->m_checkForBucket = rg_FALSE;
        }


        finish_in = clock();
        if ( whereIsEndVertexFound == FRONT_FEASIBLE_REGION )  {
            //  computingTime[4]:    finding end vertex of parabolic or hyperbolic edge in front feasible region
            computingTime[4] += finish_in - start_in;
        }
        else if ( whereIsEndVertexFound == FRONTREAR_FEASIBLE_REGION ) {
            //  computingTime[5]:    finding end vertex of parabolic or hyperbolic edge from front feasible region
            computingTime[5] += finish_in - start_in;
        }
        else if ( whereIsEndVertexFound == REAR_FEASIBLE_REGION ) {
            //  computingTime[6]:    finding end vertex of parabolic or hyperbolic edge in rear feasible region
            computingTime[6] += finish_in - start_in;
        }
        else {
        }


        /*
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //
        if ( m_bStartCheck )  {

            m_frequencyOfCandidateBall[numCheckedBalls]++;
            if ( whereIsEndVertexFound == FRONT_FEASIBLE_REGION )  {
                m_numCandidateBallsOfNonElliptic[1][0][0] += numCheckedBalls;

                m_avgRadiusOfSphereForRf[0][0] += sphereForR1.getRadius();
                m_numEdges[0][0]++;


                m_numSearchedBucketElements[0][0] += numBucketElements;
            }
            else if ( whereIsEndVertexFound == FRONTREAR_FEASIBLE_REGION ) {
                m_numCandidateBallsOfNonElliptic[1][0][1] += numCheckedBalls;

                m_avgRadiusOfSphereForRf[0][1] += sphereForR1.getRadius();
                m_numEdges[0][1]++;

                m_numSearchedBucketElements[0][1] += numBucketElements;
            }
            else if ( whereIsEndVertexFound == REAR_FEASIBLE_REGION ) {
                m_numCandidateBallsOfNonElliptic[1][0][2] += numCheckedBalls;

                m_avgRadiusOfSphereForRf[0][2] += sphereForR1.getRadius();
                m_numEdges[0][2]++;

                m_numSearchedBucketElements[0][2] += numBucketElements;
            }
            else {
            }


            rg_INT numCandidateEndBalls[3] = { 0, 0 , 0};
            BallGenerator* currBall = rg_NULL;
            m_feasibleBalls.reset4Loop();
            while ( m_feasibleBalls.setNext4Loop() )  {
                currBall = m_feasibleBalls.getEntity();

                currBall->m_checkForBucket = rg_TRUE;

                rg_REAL    distToGatePlane = normalOfGatePlane.innerProduct(currBall->getCenter()) + gatePlane[3];
                if ( rg_ABS(distToGatePlane) > currBall->getRadius() )  {
                    if ( distToGatePlane > 0 )
                        numCandidateEndBalls[2]++;
                    else
                        numCandidateEndBalls[0]++;
                }
                else {
                    numCandidateEndBalls[1]++;
                }
            }

            for ( rg_INT i=0; i<3; i++)  
                gateBall[i]->m_checkForBucket = rg_TRUE;

            m_GlobalGeneratorList.reset4Loop();
            while ( m_GlobalGeneratorList.setNext4Loop() )  {
                currBall = m_GlobalGeneratorList.getpEntity();

                if ( currBall->m_checkForBucket )
                    continue;

                rg_REAL    distToGatePlane = normalOfGatePlane.innerProduct(currBall->getCenter()) + gatePlane[3];
                if ( distToGatePlane + currBall->getRadius() > 0 )
                    numCandidateEndBalls[2]++;
            }

            m_numCandidateBallsOfNonElliptic[0][0][0] += numCandidateEndBalls[0];
            m_numCandidateBallsOfNonElliptic[0][0][1] += numCandidateEndBalls[1];
            m_numCandidateBallsOfNonElliptic[0][0][2] += numCandidateEndBalls[2];


            for ( i=0; i<3; i++)  
                gateBall[i]->m_checkForBucket = rg_FALSE;

            m_feasibleBalls.reset4Loop();
            while ( m_feasibleBalls.setNext4Loop() )  
                m_feasibleBalls.getEntity()->m_checkForBucket = rg_FALSE;


            Sphere endBall = candidateEndBall.first->getBall();
            rg_REAL    distBetEndAndGatePlane = normalOfGatePlane.innerProduct(endBall.getCenter()) + gatePlane[3];
            if ( rg_ABS(distBetEndAndGatePlane) > endBall.getRadius() )  {
                if ( distBetEndAndGatePlane > 0 )
                    m_numEndBallsInPosition[0][2]++;
                else  {
                    m_numEndBallsInPosition[0][0]++;
                    m_numIntersectingBallsWithFrontFR[0] += numCandidateEndBalls[0];
                }
            }
            else {
                m_numEndBallsInPosition[0][1]++;
                m_numIntersectingBallsWithFrontFR[0] += numCandidateEndBalls[1];
            }



            clock_t finish_extract = clock();
            timeToBeExtracted[0] += finish_extract - finish_in;
        }
        m_feasibleBalls.removeAll();
        //
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //*/

        return rg_TRUE;
    }
}




void rg_SphereSetVoronoiDiagram::initializePropagationInRearFeasibleRegionOfNonEllipticEdge(
                                               rg_dList< NavigatingBucketCell >& queueOfCells,
                                               rg_Point3D&                   unitLength,
                                               Sphere*                       gateBall,
                                               rg_REAL*                      gatePlane )
{
    rg_Point3D normalOfGatePlane(gatePlane[0], gatePlane[1], gatePlane[2]);
    rg_Point3D ptOnGatePlane[3];
    rg_INT i=0;
	for (i=0; i<3; i++)  {
        ptOnGatePlane[i] = gateBall[i].getCenter() + ( gateBall[i].getRadius()*normalOfGatePlane );    
    }
    rg_Point3D initPt = (ptOnGatePlane[0] + ptOnGatePlane[1] + ptOnGatePlane[2])/3.0;

    BucketCellIndex initIndex;
    
    m_bucket.getJustifiedBucketIndex(initPt, initIndex.m_i, initIndex.m_j, initIndex.m_k);
    rg_Point3D ptOfInitIndex( m_bucket.getMinPointOfBoundingBox() );
    rg_Point3D increment( initIndex.m_i * m_bucket.getXLengthOfCell(),
                          initIndex.m_j * m_bucket.getYLengthOfCell(),
                          initIndex.m_k * m_bucket.getZLengthOfCell() );
    ptOfInitIndex += increment;

    rg_REAL  distBetInitPtAndPg = normalOfGatePlane.innerProduct(ptOfInitIndex) + gatePlane[3];   
    

    rg_REAL  xUnitLength = gatePlane[0] * m_bucket.getXLengthOfCell();
    rg_REAL  yUnitLength = gatePlane[1] * m_bucket.getYLengthOfCell();
    rg_REAL  zUnitLength = gatePlane[2] * m_bucket.getZLengthOfCell();
    unitLength.setPoint( xUnitLength, yUnitLength, zUnitLength );

    //rg_INT faceMask[4] = {0, 0, 0, 0};
    m_bucket.setMark(initIndex, rg_TRUE);
    queueOfCells.add( NavigatingBucketCell( initIndex, ON_PLANE, distBetInitPtAndPg ) );

    /*
    BucketCellIndex neighbors[6];
    neighbors[0].setCellIndex(initIndex.m_i-1, initIndex.m_j,   initIndex.m_k);
    neighbors[1].setCellIndex(initIndex.m_i,   initIndex.m_j-1, initIndex.m_k);
    neighbors[2].setCellIndex(initIndex.m_i,   initIndex.m_j,   initIndex.m_k-1);
    neighbors[3].setCellIndex(initIndex.m_i+1, initIndex.m_j,   initIndex.m_k);
    neighbors[4].setCellIndex(initIndex.m_i,   initIndex.m_j+1, initIndex.m_k);
    neighbors[5].setCellIndex(initIndex.m_i,   initIndex.m_j,   initIndex.m_k+1);

    rg_REAL distViOfNeighbors[6];
    distViOfNeighbors[0] = distBetInitPtAndPg - unitLength.getX();
    distViOfNeighbors[1] = distBetInitPtAndPg - unitLength.getY();
    distViOfNeighbors[2] = distBetInitPtAndPg - unitLength.getZ();
    distViOfNeighbors[3] = distBetInitPtAndPg + unitLength.getX();
    distViOfNeighbors[4] = distBetInitPtAndPg + unitLength.getY();
    distViOfNeighbors[5] = distBetInitPtAndPg + unitLength.getZ();

    rg_REAL distBetPAndV[8];
    distBetPAndV[0]= distBetInitPtAndPg;
    distBetPAndV[1]= distBetInitPtAndPg + unitLength.getX();
    distBetPAndV[2]= distBetInitPtAndPg + unitLength.getX() + unitLength.getY();
    distBetPAndV[3]= distBetInitPtAndPg                     + unitLength.getY();
    distBetPAndV[4]= distBetPAndV[0] + unitLength.getZ();
    distBetPAndV[5]= distBetPAndV[1] + unitLength.getZ();
    distBetPAndV[6]= distBetPAndV[2] + unitLength.getZ();
    distBetPAndV[7]= distBetPAndV[3] + unitLength.getZ();


    //  construct masks of initial cell and its vertices.
    rg_INT cellMask = 0;
    rg_INT vertexMask[8];
    for (i=0; i<8; i++)  {
        if ( rg_ZERO( distBetPAndV[i] ) )
            vertexMask[i] = BUCKETVERTEX_MASK[i][2];
        else if ( rg_NEG( distBetPAndV[i] ) )
            vertexMask[i] = BUCKETVERTEX_MASK[i][1];
        else 
            vertexMask[i] = BUCKETVERTEX_MASK[i][0];

        cellMask += vertexMask[i];
    }


    //  decide if each neighbor should be navigated.
    for ( i=0; i<6; i++ )  {

        if ( m_bucket.isVisited(neighbors[i]) )
            continue;

        rg_INT faceStatus = cellMask&BUCKETCELL_FACE_FILTER[i];

        // f_ijkl > 0, neighbor[i] IN_PLUS_PLANE
        if ( faceStatus == BUCKETCELL_FACE_STATUS[i][0] )  {
            m_bucket.setMark(neighbors[i], rg_TRUE);
            queueOfCells.add( NavigatingBucketCell(neighbors[i], IN_PLUS_PLANE) );
        }
        // f_ijkl < 0, neighbor[i] IN_MINUS_PLANE -> No propagated
        else if ( faceStatus == BUCKETCELL_FACE_STATUS[i][1] )  {
            continue;
        }
        // neighbor[i] ON_PLANE -> the neighbor[i] intersects with the gate plane.
        else  {
            rg_INT faceMask[4];
            for (rg_INT j=0; j<4; j++)
                faceMask[j] = vertexMask[ BUCKETCELL_FACE[i][j] ];

            m_bucket.setMark(neighbors[i], rg_TRUE);
            queueOfCells.add( NavigatingBucketCell( neighbors[i], ON_PLANE, 
                                                    distViOfNeighbors[i], 
                                                    PREV_BUCKETCELL[i],
                                                    faceMask) );
        }
    }
    */
}




rg_Point3D rg_SphereSetVoronoiDiagram::computeCenterOfApproximatingSphereForR1(
                                                                          rg_Point3D* points, 
                                                                    const rg_Point3D& normalOfCenterPlane,
                                                                    const rg_REAL&    coeffCenterPlane )
{
    rg_REAL a1 = points[0].getX() - points[1].getX();
    rg_REAL a2 = points[0].getX() - points[2].getX();
    rg_REAL a3 = normalOfCenterPlane.getX();

    rg_REAL b1 = points[0].getY() - points[1].getY();
    rg_REAL b2 = points[0].getY() - points[2].getY();
    rg_REAL b3 = normalOfCenterPlane.getY();

    rg_REAL c1 = points[0].getZ() - points[1].getZ();
    rg_REAL c2 = points[0].getZ() - points[2].getZ();
    rg_REAL c3 = normalOfCenterPlane.getZ();

    rg_REAL d1 = ( points[1].innerProduct(points[1]) - points[0].innerProduct(points[0]) )/2.0;
    rg_REAL d2 = ( points[2].innerProduct(points[2]) - points[0].innerProduct(points[0]) )/2.0;
    rg_REAL d3 = coeffCenterPlane;

    rg_REAL denominator = (a2*b3*c1-a2*b1*c3-a3*b2*c1+b1*a3*c2-b3*a1*c2+a1*b2*c3);

    rg_REAL x = -(-b2*d3*c1+b2*c3*d1+c2*b1*d3-c2*b3*d1+d2*b3*c1-d2*b1*c3)/denominator;
    rg_REAL y =  (-a2*d3*c1+a2*c3*d1-a3*c2*d1-c3*a1*d2+d3*a1*c2+a3*d2*c1)/denominator;
    rg_REAL z = -(b1*a3*d2-a2*b1*d3+a2*b3*d1-b3*a1*d2-a3*b2*d1+a1*b2*d3)/denominator;

    return rg_Point3D(x,y,z);
}




rg_FLAG rg_SphereSetVoronoiDiagram::findEndVertexOfCircularAndEllipticEdgeInAcceleration( 
                                                                    Gate&       gate, 
                                                              const rg_Point3D& normalOfCenterPlane, 
                                                              const rg_Point3D& normalOfEdgePlane,
                                                              const rg_Point3D& axisPoint,
                                                              Sphere*           tangentSphereInCP)
{

    //  This gate defines a Voronoi edge.
    BallGenerator* gateBall[3]  = { gate.getGateBallCCW1(), 
                                    gate.getGateBallCCW2(), 
                                    gate.getGateBallCCW3() };

    Sphere ball[3] = { gateBall[0]->getBall(), gateBall[1]->getBall(), gateBall[2]->getBall() };
    BallGenerator* ptrStartBall = gate.getStartBall();
    
    rg_Point3D startVertex  = gate.getStartVertex()->getPoint();
    rg_Point3D vecPaVs      = startVertex - axisPoint;
    vecPaVs.normalize();


    ///////////////////////////////////////////////////////////////////////////
    //
    //  candidateEndBall.first   : pointer to endBall
    //  candidateEndBall.second  : tangent sphere by endBall and gateBalls
    pair<BallGenerator*, Sphere> candidateEndBall;
    candidateEndBall.first   = rg_NULL;
    rg_REAL minAngularDistance = DBL_MAX;



    //  define box filter
    rg_Point3D vecTS1ToTS2 = tangentSphereInCP[1].getCenter() - tangentSphereInCP[0].getCenter();
    vecTS1ToTS2.normalize();

    rg_Point3D farPointOfTS1 = tangentSphereInCP[0].getCenter() - (tangentSphereInCP[0].getRadius()*vecTS1ToTS2);
    rg_Point3D farPointOfTS2 = tangentSphereInCP[1].getCenter() + (tangentSphereInCP[1].getRadius()*vecTS1ToTS2);

    rg_Point3D center = (farPointOfTS1 + farPointOfTS2)/2.0;
    rg_REAL    radius = center.distance( farPointOfTS1 );
    Sphere boxfilter(center, radius);    


    /*
    //  define excluding plane.
    rg_REAL radiusOfThirdSphere = (tangentSphereInCP[0].getRadius() + tangentSphereInCP[1].getRadius())/2.0;
    Sphere  thirdSphere[2];
    
    computeSphereWithGivenRadiusAndTangentTo3Spheres(radiusOfThirdSphere, ball, thirdSphere);

    rg_Point3D vec1 = thirdSphere[0].getCenter()       - tangentSphereInCP[0].getCenter();
    rg_Point3D vec2 = tangentSphereInCP[1].getCenter() - tangentSphereInCP[0].getCenter();
    rg_Point3D directionTP = vec1.crossProduct( vec2 ).getUnitVector();

    rg_REAL excludingPlane[2][4];
    Sphere  tangentBall[2][3];
    if ( rg_GE( normalOfEdgePlane.innerProduct( directionTP ), 0.0 ) )  {
         tangentBall[0][0] = tangentSphereInCP[0]; 
         tangentBall[0][1] = thirdSphere[0]; 
         tangentBall[0][2] = tangentSphereInCP[1]; 
         tangentBall[1][0] = tangentSphereInCP[0]; 
         tangentBall[1][1] = thirdSphere[1]; 
         tangentBall[1][2] = tangentSphereInCP[1]; 
    }
    else  {
         tangentBall[0][0] = tangentSphereInCP[0]; 
         tangentBall[0][1] = thirdSphere[1]; 
         tangentBall[0][2] = tangentSphereInCP[1]; 
         tangentBall[1][0] = tangentSphereInCP[0]; 
         tangentBall[1][1] = thirdSphere[0]; 
         tangentBall[1][2] = tangentSphereInCP[1]; 
    }
       
    computePlaneTangentTo3SphereFromCCWOutside(tangentBall[0], excludingPlane[0]);
    computePlaneTangentTo3SphereFromCCWOutside(tangentBall[1], excludingPlane[1]);  

    rg_Point3D normalEP1(excludingPlane[0][0], excludingPlane[0][1], excludingPlane[0][2]);
    rg_Point3D normalEP2(excludingPlane[1][0], excludingPlane[1][1], excludingPlane[2][2]);
    rg_REAL distToEP[2];
    */


    gateBall[0]->m_checkForBucket = rg_TRUE;
    gateBall[1]->m_checkForBucket = rg_TRUE;
    gateBall[2]->m_checkForBucket = rg_TRUE;

    rg_dList<BallGenerator*> ballSetForBoxfilter;
    m_bucket.getCellsByMask( ballSetForBoxfilter, boxfilter );

    findEmptyTangentSphereWithMinAngleInEllipticAndCircularEdge( 
                        candidateEndBall, minAngularDistance, ballSetForBoxfilter,
                        ball, ptrStartBall, startVertex, normalOfEdgePlane,
					    vecPaVs, axisPoint);


    gateBall[0]->m_checkForBucket = rg_FALSE;
    gateBall[1]->m_checkForBucket = rg_FALSE;
    gateBall[2]->m_checkForBucket = rg_FALSE;


    if ( candidateEndBall.first == rg_NULL )  {
        ballSetForBoxfilter.reset4Loop();
	    while( ballSetForBoxfilter.setNext4Loop() )
	        ballSetForBoxfilter.getEntity()->m_checkForBucket = rg_FALSE;

        return rg_FALSE;
    }
    else  {
        gate.setEndBall( candidateEndBall.first );
        gate.setEmptyTangentSphereAtEndVertex( candidateEndBall.second );

        /*
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //
        m_numCandidateBallsOfElliptic[1] += ballSetForBoxfilter.getSize();
        //
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //*/

        ballSetForBoxfilter.reset4Loop();
	    while( ballSetForBoxfilter.setNext4Loop() )
	        ballSetForBoxfilter.getEntity()->m_checkForBucket = rg_FALSE;

        return rg_TRUE;
    }
}




rg_FLAG rg_SphereSetVoronoiDiagram::findEndVertexOfLinearEdgeInAcceleration(       Gate&       gate, 
                                                                             const rg_Point3D& normalOfCenterPlane, 
                                                                             const Sphere&     minSphere)
{
    //  This gate defines a Voronoi edge.
    BallGenerator* gateBall[3]  = { gate.getGateBallCCW1(), gate.getGateBallCCW2(), gate.getGateBallCCW3() };  
    Sphere         ball[3]      = { gateBall[0]->getBall(), gateBall[1]->getBall(), gateBall[2]->getBall() };
    BallGenerator* ptrStartBall = gate.getStartBall();
    rg_Point3D     startVertex  = gate.getStartVertex()->getPoint();

    ///////////////////////////////////////////////////////////////////////////
    //
    //  candidateEndBall.first   : pointer to endBall
    //  candidateEndBall.second  : tangent sphere by endBall and gateBalls
    pair<BallGenerator*, Sphere> candidateEndBall;
    candidateEndBall.first   = rg_NULL;
    rg_REAL minAngularDistance = DBL_MAX;


    //////////////////////////////////////////////////////////////////////////
    //
    //  For Front Feasible Region
    //
    rg_INT numBucketElements = 0;
    //  construct sphere enclosing R1
    Sphere sphereForR1;

    rg_Point3D ptOnGatePlane = ball[0].getCenter() + ( ball[0].getRadius()*normalOfCenterPlane );

    sphereForR1.setCenter( minSphere.getCenter() );   
    sphereForR1.setRadius( minSphere.getCenter().distance( ptOnGatePlane ) );

    //  SHRINK
    //sphereForR1.setRadius( minSphere.getCenter().distance( ptToIntersectTangentPlaneOfGate ) + m_minRadiusOfBalls );
    //sphereForR1.setRadius( minSphere.getCenter().distance( ptToIntersectTangentPlaneOfGate ) + m_maxRadiusOfBalls );

    rg_REAL gatePlane[4] = { normalOfCenterPlane.getX(), normalOfCenterPlane.getY(), normalOfCenterPlane.getZ(), 
                             -normalOfCenterPlane.innerProduct(ptOnGatePlane) };


    rg_dList<BallGenerator*> checkedBallSet;
    rg_INT i=0;
	for (i=0; i<3; i++)  {
        gateBall[i]->m_checkForBucket = rg_TRUE;
        checkedBallSet.add(gateBall[i]);
    }
 
    rg_dList<BallGenerator*> ballSetForR1;
    numBucketElements = m_bucket.getCellsByMask( ballSetForR1, sphereForR1 );
    
    m_isFrontFeasibleRegion = rg_TRUE;
    findEmptyTangentSphereWithMinAngleInLinearEdge( candidateEndBall, minAngularDistance, ballSetForR1,
                                                    ball, ptrStartBall, startVertex, normalOfCenterPlane );
    m_isFrontFeasibleRegion = rg_FALSE;

    checkedBallSet.mergeTail( ballSetForR1 );


    
    rg_INT whereIsEndVertexFound = NOT_FOUND;

    

    if (  candidateEndBall.first != rg_NULL )
    {
        if ( candidateEndBall.second.isContainedIn( sphereForR1 ) == rg_FALSE ) {

            rg_dList<BallGenerator*> ballSetForTS_C;
            numBucketElements += m_bucket.getCellsByMask( ballSetForTS_C, candidateEndBall.second );

            findEmptyTangentSphereWithMinAngleInLinearEdge( candidateEndBall, minAngularDistance, ballSetForTS_C,
                                                            ball, ptrStartBall, startVertex, normalOfCenterPlane );

            checkedBallSet.mergeTail( ballSetForTS_C );

            whereIsEndVertexFound = FRONTREAR_FEASIBLE_REGION;
        }  
        else  {
            whereIsEndVertexFound = FRONT_FEASIBLE_REGION;
        }
    }
    else  {
    
        //////////////////////////////////////////////////////////////////////////
        //
        //  For Rear Feasible Region
        //

        rg_dList< NavigatingBucketCell > queueOfCells;

        rg_Point3D unitLength;

        initializePropagationInRearFeasibleRegionOfNonEllipticEdge( queueOfCells, unitLength, ball, gatePlane );

        rg_INT  normalSign = classifyNormalSign(gatePlane[0], gatePlane[1], gatePlane[2]);



        NavigatingBucketCell* currBCell = rg_NULL;
    
        rg_dNode< NavigatingBucketCell >* currNode = queueOfCells.getHead();
        //currNode = currNode->getNext();
        do 
	    {
            currBCell = currNode->getpEntity();
            rg_dList< BallGenerator* > candidateBalls;

            if ( currBCell->m_status == ON_PLANE )   
                m_bucket.getBalls( candidateBalls, currBCell->m_index, gatePlane );
            else  
                m_bucket.getBalls( candidateBalls, currBCell->m_index );

            if ( candidateBalls.getSize() > 0 )  {

                findEmptyTangentSphereWithMinAngleInLinearEdge( candidateEndBall, minAngularDistance, candidateBalls,
                                                                ball, ptrStartBall, startVertex, normalOfCenterPlane );


                checkedBallSet.mergeTail( candidateBalls );               
                if ( candidateEndBall.first != rg_NULL )
                    break;
            }


            if ( currBCell->m_status == ON_PLANE )   
                propagateBucketCellsOnPlane( *currBCell, unitLength, queueOfCells );
            else  
                propagateBucketCellsInPlusPlane( *currBCell, normalSign, queueOfCells );

            currNode = currNode->getNext();
        
        } while( currNode != queueOfCells.getHead() );


        if ( candidateEndBall.first != rg_NULL )  {

            rg_dList<BallGenerator*> ballSetForR2;
            numBucketElements += m_bucket.getCellsByMask( ballSetForR2, candidateEndBall.second );

            findEmptyTangentSphereWithMinAngleInLinearEdge( candidateEndBall, minAngularDistance, ballSetForR2,
                                                            ball, ptrStartBall, startVertex, normalOfCenterPlane );

            checkedBallSet.mergeTail( ballSetForR2 );
        }

        numBucketElements += queueOfCells.getSize();
        queueOfCells.reset4Loop();
	    while( queueOfCells.setNext4Loop() )  {
	        currBCell = queueOfCells.getpEntity();
            m_bucket.setMark(currBCell->m_index, rg_FALSE);
        }

        whereIsEndVertexFound = REAR_FEASIBLE_REGION;
    }

    

    if ( candidateEndBall.first == rg_NULL )  {

        /*
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //
        m_numSearchedBucketElements[1][3] += numBucketElements;

        m_numBallsForUnboundedEdge[1] += checkedBallSet.getSize();
        //
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //*/

        checkedBallSet.reset4Loop();
	    while( checkedBallSet.setNext4Loop() )  {
	        checkedBallSet.getEntity()->m_checkForBucket = rg_FALSE;
        }

        return rg_FALSE;
    }
    else  {
        gate.setEndBall(candidateEndBall.first);
        gate.setEmptyTangentSphereAtEndVertex( candidateEndBall.second );

        rg_INT numCheckedBalls = checkedBallSet.getSize();
        checkedBallSet.reset4Loop();
	    while( checkedBallSet.setNext4Loop() )  {
	        checkedBallSet.getEntity()->m_checkForBucket = rg_FALSE;
        }


        finish_in = clock();
        if ( whereIsEndVertexFound == FRONT_FEASIBLE_REGION )  {
            //  computingTime[8]:    finding end vertex of linear edge in front feasible region
            computingTime[8] += finish_in - start_in;
        }
        else if ( whereIsEndVertexFound == FRONTREAR_FEASIBLE_REGION ) {
            //  computingTime[9]:    finding end vertex of linear edge from front feasible region
            computingTime[9] += finish_in - start_in;
        }
        else if ( whereIsEndVertexFound == REAR_FEASIBLE_REGION ) {
            //  computingTime[10]:    finding end vertex of linear edge in rear feasible region
            computingTime[10] += finish_in - start_in;
        }
        else {
        }


        /*
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //
        if ( m_bStartCheck )  {

            m_frequencyOfCandidateBall[numCheckedBalls]++;
            if ( whereIsEndVertexFound == FRONT_FEASIBLE_REGION )  {
                m_numCandidateBallsOfNonElliptic[1][1][0] += numCheckedBalls;

                m_avgRadiusOfSphereForRf[1][0] += sphereForR1.getRadius();
                m_numEdges[1][0]++;

                m_numSearchedBucketElements[1][0] += numBucketElements;
            }
            else if ( whereIsEndVertexFound == FRONTREAR_FEASIBLE_REGION ) {
                m_numCandidateBallsOfNonElliptic[1][1][1] += numCheckedBalls;

                m_avgRadiusOfSphereForRf[1][1] += sphereForR1.getRadius();
                m_numEdges[1][1]++;

                m_numSearchedBucketElements[1][1] += numBucketElements;
            }
            else if ( whereIsEndVertexFound == REAR_FEASIBLE_REGION ) {
                m_numCandidateBallsOfNonElliptic[1][1][2] += numCheckedBalls;

                m_avgRadiusOfSphereForRf[1][2] += sphereForR1.getRadius();
                m_numEdges[1][2]++;

                m_numSearchedBucketElements[1][2] += numBucketElements;
            }
            else {
            }


            rg_INT numCandidateEndBalls[3] = {0, 0, 0};
            BallGenerator* currBall = rg_NULL;
            m_feasibleBalls.reset4Loop();
            while ( m_feasibleBalls.setNext4Loop() )  {
                currBall = m_feasibleBalls.getEntity();

                currBall->m_checkForBucket = rg_TRUE;

                rg_REAL    distToGatePlane = normalOfCenterPlane.innerProduct(currBall->getCenter()) + gatePlane[3];
                if ( rg_ABS(distToGatePlane) > currBall->getRadius() )  {
                    if ( distToGatePlane > 0 )
                        numCandidateEndBalls[2]++;
                    else
                        numCandidateEndBalls[0]++;
                }
                else {
                    numCandidateEndBalls[1]++;
                }
            }

            for ( rg_INT i=0; i<3; i++)  
                gateBall[i]->m_checkForBucket = rg_TRUE;

            m_GlobalGeneratorList.reset4Loop();
            while ( m_GlobalGeneratorList.setNext4Loop() )  {
                currBall = m_GlobalGeneratorList.getpEntity();

                if ( currBall->m_checkForBucket )
                    continue;

                rg_REAL    distToGatePlane = normalOfCenterPlane.innerProduct(currBall->getCenter()) + gatePlane[3];
                if ( distToGatePlane + currBall->getRadius() > 0 )
                    numCandidateEndBalls[2]++;
            }

            m_numCandidateBallsOfNonElliptic[0][1][0] += numCandidateEndBalls[0];
            m_numCandidateBallsOfNonElliptic[0][1][1] += numCandidateEndBalls[1];
            m_numCandidateBallsOfNonElliptic[0][1][2] += numCandidateEndBalls[2];

            
            for ( i=0; i<3; i++)  
                gateBall[i]->m_checkForBucket = rg_FALSE;

            m_feasibleBalls.reset4Loop();
            while ( m_feasibleBalls.setNext4Loop() )  
                m_feasibleBalls.getEntity()->m_checkForBucket = rg_FALSE;


            Sphere endBall = candidateEndBall.first->getBall();
            rg_REAL    distBetEndAndGatePlane = normalOfCenterPlane.innerProduct(endBall.getCenter()) + gatePlane[3];
            if ( rg_ABS(distBetEndAndGatePlane) > endBall.getRadius() )  {
                if ( distBetEndAndGatePlane > 0 )
                    m_numEndBallsInPosition[1][2]++;
                else  {
                    m_numEndBallsInPosition[1][0]++;
                    m_numIntersectingBallsWithFrontFR[1] += numCandidateEndBalls[0];
                }
            }
            else {
                m_numEndBallsInPosition[1][1]++;
                m_numIntersectingBallsWithFrontFR[1] += numCandidateEndBalls[1];
            }
            

            clock_t finish_extract = clock();
            timeToBeExtracted[1] += finish_extract - finish_in;
        }
        m_feasibleBalls.removeAll();
        //
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //*/


        return rg_TRUE;
    }
}




void    rg_SphereSetVoronoiDiagram::propagateBucketCells(
                                         const BucketCellIndex& cell, 
                                         const rg_REAL&         distance,
                                         const rg_Point3D&      unitLength, 
                                         rg_dList< rg_Triple< rg_INT, BucketCellIndex, rg_REAL > >* cellList)
{
    BucketCellIndex cellToBePropagated[6];

    cellToBePropagated[0].setCellIndex(cell.m_i-1, cell.m_j, cell.m_k);
    cellToBePropagated[1].setCellIndex(cell.m_i, cell.m_j-1, cell.m_k);
    cellToBePropagated[2].setCellIndex(cell.m_i, cell.m_j, cell.m_k-1);
    cellToBePropagated[3].setCellIndex(cell.m_i+1, cell.m_j, cell.m_k);
    cellToBePropagated[4].setCellIndex(cell.m_i, cell.m_j+1, cell.m_k);
    cellToBePropagated[5].setCellIndex(cell.m_i, cell.m_j, cell.m_k+1);

    rg_REAL neighborBucket[6];
    neighborBucket[0] = distance - unitLength.getX();
    neighborBucket[1] = distance - unitLength.getY();
    neighborBucket[2] = distance - unitLength.getZ();
    neighborBucket[3] = distance + unitLength.getX();
    neighborBucket[4] = distance + unitLength.getY();
    neighborBucket[5] = distance + unitLength.getZ();

    rg_INT i=0;
	for ( i=0; i<6; i++)  {

        rg_FLAG posCell = m_bucket.getPositionOfCellWRTPlane(cellToBePropagated[i],
                                                             neighborBucket[i],
                                                             unitLength);

        if ( posCell == ON_PLANE || posCell == IN_PLUS_PLANE )  {
            cellList->addTail( rg_Triple< rg_INT, BucketCellIndex, rg_REAL >(posCell, cellToBePropagated[i], neighborBucket[i]) );

            m_bucket.setMark(cellToBePropagated[i], rg_TRUE);
        }
    }
}



void rg_SphereSetVoronoiDiagram::propagateBucketCellsOnPlane( 
                                                 const NavigatingBucketCell&       cell, 
                                                 const rg_Point3D&                 unitLength, 
                                                 rg_dList< NavigatingBucketCell >& queueOfCells)
{      
    
    BucketCellIndex neighbors[6];
    neighbors[0].setCellIndex(cell.m_index.m_i-1, cell.m_index.m_j,   cell.m_index.m_k);
    neighbors[1].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j-1, cell.m_index.m_k);
    neighbors[2].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k-1);
    neighbors[3].setCellIndex(cell.m_index.m_i+1, cell.m_index.m_j,   cell.m_index.m_k);
    neighbors[4].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j+1, cell.m_index.m_k);
    neighbors[5].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k+1);

    rg_REAL distViOfNeighbors[6];
    distViOfNeighbors[0] = cell.m_distBetPandVi - unitLength.getX();
    distViOfNeighbors[1] = cell.m_distBetPandVi - unitLength.getY();
    distViOfNeighbors[2] = cell.m_distBetPandVi - unitLength.getZ();
    distViOfNeighbors[3] = cell.m_distBetPandVi + unitLength.getX();
    distViOfNeighbors[4] = cell.m_distBetPandVi + unitLength.getY();
    distViOfNeighbors[5] = cell.m_distBetPandVi + unitLength.getZ();

    rg_REAL distBetPAndV[8];
    //cell.computeVisFromPrevCell(distBetPAndV, unitLength);
    //distBetPAndV[0]= cell.m_distBetPandVi;
    //distBetPAndV[1]= cell.m_distBetPandVi + unitLength.getX();
    //distBetPAndV[2]= cell.m_distBetPandVi + unitLength.getX() + unitLength.getY();
    //distBetPAndV[3]= cell.m_distBetPandVi                     + unitLength.getY();
    //distBetPAndV[4]= distBetPAndV[0] + unitLength.getZ();
    //distBetPAndV[5]= distBetPAndV[1] + unitLength.getZ();
    //distBetPAndV[6]= distBetPAndV[2] + unitLength.getZ();
    //distBetPAndV[7]= distBetPAndV[3] + unitLength.getZ();
   
    distBetPAndV[0]= cell.m_distBetPandVi;
    distBetPAndV[1]= distViOfNeighbors[3];
    distBetPAndV[2]= distViOfNeighbors[3] + unitLength.getY();
    distBetPAndV[3]= distViOfNeighbors[4];
    distBetPAndV[4]= distViOfNeighbors[5];
    distBetPAndV[5]= distBetPAndV[1] + unitLength.getZ();
    distBetPAndV[6]= distBetPAndV[2] + unitLength.getZ();
    distBetPAndV[7]= distBetPAndV[3] + unitLength.getZ();
    


    //  construct masks of initial cell and its vertices.
    rg_INT cellMask = 0;
    rg_INT vertexMask[8];
    rg_INT i=0;
	for ( i=0; i<8; i++)  {
        if ( rg_ZERO( distBetPAndV[i] ) )
            vertexMask[i] = BUCKETVERTEX_MASK[i][2];
        else if ( rg_NEG( distBetPAndV[i] ) )
            vertexMask[i] = BUCKETVERTEX_MASK[i][1];
        else 
            vertexMask[i] = BUCKETVERTEX_MASK[i][0];

        cellMask += vertexMask[i];
    }


    //  decide if each neighbor should be navigated.
    for ( i=0; i<6; i++ )  {

        if ( m_bucket.isVisited(neighbors[i]) )
            continue;

        rg_INT faceStatus = cellMask&BUCKETCELL_FACE_FILTER[i];

        // f_ijkl > 0, neighbor[i] IN_PLUS_PLANE
        if ( faceStatus == BUCKETCELL_FACE_STATUS[i][0] )  {
            m_bucket.setMark(neighbors[i], rg_TRUE);
            queueOfCells.add( NavigatingBucketCell(neighbors[i], IN_PLUS_PLANE) );
        }
        // f_ijkl < 0, neighbor[i] IN_MINUS_PLANE -> No propagated
        else if ( faceStatus == BUCKETCELL_FACE_STATUS[i][1] )  {
            continue;
        }
        // neighbor[i] ON_PLANE -> the neighbor[i] intersects with the gate plane.
        else  {
            //rg_INT faceMask[4];
            //for (rg_INT j=0; j<4; j++)
            //    faceMask[j] = vertexMask[ BUCKETCELL_FACE[i][j] ];

            m_bucket.setMark(neighbors[i], rg_TRUE);
            queueOfCells.add( NavigatingBucketCell( neighbors[i], ON_PLANE, distViOfNeighbors[i]) );
            //queueOfCells.add( NavigatingBucketCell( neighbors[i], ON_PLANE, 
            //                                        distViOfNeighbors[i], 
            //                                        PREV_BUCKETCELL[i],
            //                                        faceMask) );
        }
    }    
}



void rg_SphereSetVoronoiDiagram::propagateBucketCellsInPlusPlane( 
                                                 const NavigatingBucketCell&       cell, 
                                                 const rg_INT&                     normalSign,
                                                 rg_dList< NavigatingBucketCell >& queueOfCells)
{
    BucketCellIndex neighbors[3];

    // plane: ax + by + cz + d =0
    switch ( normalSign )  {
    case XGE_YGE_ZGE:  // a>=0, b>=0, c>=0
        neighbors[0].setCellIndex(cell.m_index.m_i+1, cell.m_index.m_j,   cell.m_index.m_k);
        neighbors[1].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j+1, cell.m_index.m_k);
        neighbors[2].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k+1);
        break;
    
    case XGE_YGE_ZLT:  // a>=0, b>=0, c<0
        neighbors[0].setCellIndex(cell.m_index.m_i+1, cell.m_index.m_j,   cell.m_index.m_k);
        neighbors[1].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j+1, cell.m_index.m_k);
        neighbors[2].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k-1);
        break;
    
    case XGE_YLT_ZGE:  // a>=0, b<0, c>=0
        neighbors[0].setCellIndex(cell.m_index.m_i+1, cell.m_index.m_j,   cell.m_index.m_k);
        neighbors[1].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j-1, cell.m_index.m_k);
        neighbors[2].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k+1);
        break;
    
    case XGE_YLT_ZLT:  // a>=0, b<0, c<0
        neighbors[0].setCellIndex(cell.m_index.m_i+1, cell.m_index.m_j,   cell.m_index.m_k);
        neighbors[1].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j-1, cell.m_index.m_k);
        neighbors[2].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k-1);
        break;
    
    case XLT_YGE_ZGE:  // a<0, b>=0, c>=0
        neighbors[0].setCellIndex(cell.m_index.m_i-1, cell.m_index.m_j,   cell.m_index.m_k);
        neighbors[1].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j+1, cell.m_index.m_k);
        neighbors[2].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k+1);
        break;
    
    case XLT_YGE_ZLT:  // a<0, b>=0, c<0
        neighbors[0].setCellIndex(cell.m_index.m_i-1, cell.m_index.m_j,   cell.m_index.m_k);
        neighbors[1].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j+1, cell.m_index.m_k);
        neighbors[2].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k-1);
        break;
    
    case XLT_YLT_ZGE:  // a<0, b<0, c>=0
        neighbors[0].setCellIndex(cell.m_index.m_i-1, cell.m_index.m_j,   cell.m_index.m_k);
        neighbors[1].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j-1, cell.m_index.m_k);
        neighbors[2].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k+1);
        break;

    case XLT_YLT_ZLT:  // a<0, b<0, c<0
        neighbors[0].setCellIndex(cell.m_index.m_i-1, cell.m_index.m_j,   cell.m_index.m_k);
        neighbors[1].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j-1, cell.m_index.m_k);
        neighbors[2].setCellIndex(cell.m_index.m_i,   cell.m_index.m_j,   cell.m_index.m_k-1);
        break;

    default:
        break;
    }

    rg_INT i=0;
	for ( i=0; i<3; i++ )  {

        if ( !m_bucket.setMarkAfterCheck(neighbors[i], rg_TRUE) )
            continue;

        queueOfCells.add( NavigatingBucketCell(neighbors[i], IN_PLUS_PLANE) );
    }
}



void rg_SphereSetVoronoiDiagram::computeTangentPlaneOfGate(const Gate& gate, rg_REAL* tangentPlane)
{
	int i=0;
    BallGenerator* cells[3] = {gate.getGateBallCCW1(), gate.getGateBallCCW2(), gate.getGateBallCCW3()};

	Sphere generator[3];
	for(i=0; i<3; i++)
	{
		generator[i] = cells[i]->getBall();
	}

	//make the index of the smallest generator  0.
	rg_REAL radius0 = generator[0].getRadius();
	rg_REAL radius1 = generator[1].getRadius();
	rg_REAL radius2 = generator[2].getRadius();

	if( radius1 <= radius2 && radius1 <= radius0 )
	{
		Sphere temp = generator[0];
		generator[0] = generator[1];
		generator[1] = generator[2];
		generator[2] = temp;
	}
	else if( radius2 <= radius1 && radius2 <= radius0 )
	{
		Sphere temp = generator[0];
		generator[0] = generator[2];
		generator[2] = generator[1];
		generator[1] = temp;
	}
	

	//normal vector of the plane passing through three centers of generators
	//toward infinity vertex
	rg_Point3D vec1 = generator[1].getCenter() - generator[0].getCenter();
	rg_Point3D vec2 = generator[2].getCenter() - generator[1].getCenter();
	rg_Point3D normalVec = vec1.crossProduct(vec2).getUnitVector();

	rg_REAL x1 = generator[0].getCenter().getX();
	rg_REAL y1 = generator[0].getCenter().getY();
	rg_REAL z1 = generator[0].getCenter().getZ();
	rg_REAL r1 = generator[0].getRadius();

	rg_REAL x2 = generator[1].getCenter().getX();
	rg_REAL y2 = generator[1].getCenter().getY();
	rg_REAL z2 = generator[1].getCenter().getZ();
	rg_REAL r2 = generator[1].getRadius();

	rg_REAL x3 = generator[2].getCenter().getX();
	rg_REAL y3 = generator[2].getCenter().getY();
	rg_REAL z3 = generator[2].getCenter().getZ();
	rg_REAL r3 = generator[2].getRadius();

	bool bPerturb_x = false;
	bool bPerturb_y = false;
	bool bPerturb_z = false;
	rg_REAL den = (y1*x3*z2+x1*y2*z3+x2*y3*z1-y3*x1*z2-x3*y2*z1-x2*y1*z3);
	if(rg_ZERO(den))
	{
		bPerturb_x = true;

		x1 = x1 + 1;
		x2 = x2 + 1;
		x3 = x3 + 1;
		den = (y1*x3*z2+x1*y2*z3+x2*y3*z1-y3*x1*z2-x3*y2*z1-x2*y1*z3);
	}
	if(rg_ZERO(den))
	{
		bPerturb_x = false;
		bPerturb_y = true;
		
		x1 = x1 - 1;
		x2 = x2 - 1;
		x3 = x3 - 1;

		y1 = y1 + 1;
		y2 = y2 + 1;
		y3 = y3 + 1;
		den = (y1*x3*z2+x1*y2*z3+x2*y3*z1-y3*x1*z2-x3*y2*z1-x2*y1*z3);
	}
	if(rg_ZERO(den))
	{
		bPerturb_y = false;
		bPerturb_z = true;

		y1 = y1 - 1;
		y2 = y2 - 1;
		y3 = y3 - 1;

		z1 = z1 + 1;
		z2 = z2 + 1;
		z3 = z3 + 1;
		den = (y1*x3*z2+x1*y2*z3+x2*y3*z1-y3*x1*z2-x3*y2*z1-x2*y1*z3);
	}

	rg_REAL a1 = -(-y3*z2+y2*z3-y2*z1+y1*z2+y3*z1-y1*z3);
	rg_REAL a2 = -(r2*y3*z1-r2*y1*z3+y2*z3*r1-y2*r3*z1-z2*y3*r1+z2*y1*r3);
	rg_REAL b1 = x2*z3-x2*z1-x3*z2+x1*z2-x1*z3+x3*z1;
	rg_REAL b2 = x2*z3*r1-x2*r3*z1+x3*r2*z1-x3*z2*r1-z3*x1*r2+r3*x1*z2;
	rg_REAL c1 = -(-x3*y2-x2*y1+x2*y3+y1*x3+x1*y2-y3*x1);
	rg_REAL c2 = -(x2*y3*r1-x2*y1*r3-y3*x1*r2-x3*y2*r1+x1*y2*r3+y1*x3*r2);

	rg_REAL d = 1/(c1*c1+b1*b1+a1*a1)*
			(-2.0*c1*c2-2.0*b1*b2-2.0*a1*a2+2.0*
			 sqrt(2.0*c1*c2*b1*b2+2.0*c1*c2*a1*a2+
				 2.0*b1*b2*a1*a2-c1*c1*b2*b2-c1*c1*a2*a2+
				 c1*c1*den*den-b1*b1*c2*c2-b1*b1*a2*a2+
				 b1*b1*den*den-a1*a1*c2*c2-a1*a1*b2*b2+
				 a1*a1*den*den))/2.0;

	rg_REAL a = (a1*d+a2)/den;
	rg_REAL b = (b1*d+b2)/den;
	rg_REAL c = (c1*d+c2)/den;

    // vector (a,b,c)\B0\A1 \BC\BC \B1\B8\C0\C7 \C1߽\C9\C0\B8\B7\CE \C1\A4\C0ǵǴ\C2 vector\BF\CD \B9\E6\C7\E2\C0\CC \C0\CFġ\C7ϴ\C2
    // \C1\F6\B8\A6 \B0˻\E7\C7Ͽ\A9\BE\DF \C7Ѵ\D9.
	if( rg_Point3D(a,b,c).innerProduct(normalVec) < 0 )
	{

		d = 1/(c1*c1+b1*b1+a1*a1)*(-2.0*c1*c2-2.0*b1*b2-2.0*a1*a2-2.0*
			sqrt(2.0*c1*c2*b1*b2+2.0*c1*c2*a1*a2+2.0*b1*b2*a1*a2-
				c1*c1*b2*b2-c1*c1*a2*a2+c1*c1*den*den-
				b1*b1*c2*c2-b1*b1*a2*a2+b1*b1*den*den-
				a1*a1*c2*c2-a1*a1*b2*b2+a1*a1*den*den))/2.0;

		a = (a1*d+a2)/den;
		b = (b1*d+b2)/den;
		c = (c1*d+c2)/den;
	}

	if(bPerturb_x)
	{
		d = d + a;
	}
	else if(bPerturb_y)
	{
		d = d + b;
	}
	else if(bPerturb_z)
	{
		d = d + c;
	}
    tangentPlane[0] = a;
    tangentPlane[1] = b;
    tangentPlane[2] = c;
    tangentPlane[3] = d;
}



rg_FLAG rg_SphereSetVoronoiDiagram::computeEdgePlaneForAngleDistanceAccordingToEdgeType(
                                                 const Gate&       gate, 
                                                       rg_Point3D& normalOfCenterPlane,
                                                       rg_Point3D& normalOfEdgePlane,
                                                       rg_Point3D& axisPoint,
                                                       Sphere*     tangentSphereOnCP)
{
    Sphere gateBall[3] = { gate.getGateBallCCW1()->getBall(), 
                           gate.getGateBallCCW2()->getBall(), 
                           gate.getGateBallCCW3()->getBall() }; 

    rg_Point3D centerOfGB[3] = { gateBall[0].getCenter(), gateBall[1].getCenter(), gateBall[2].getCenter() };
    rg_REAL    radiusOfGB[3] = { gateBall[0].getRadius(), gateBall[1].getRadius(), gateBall[2].getRadius() };

    //  If the gate balls define a linear edge,
    //  we use Euclidean distance between start and end vertices as the angular distance.
    if ( rg_EQ(radiusOfGB[0], radiusOfGB[1]) && rg_EQ(radiusOfGB[0], radiusOfGB[2]) && rg_EQ(radiusOfGB[1], radiusOfGB[2]) )
    {
	    rg_Point3D vec21 = centerOfGB[1] - centerOfGB[0];
	    rg_Point3D vec31 = centerOfGB[2] - centerOfGB[0];
	    normalOfCenterPlane = vec21.crossProduct(vec31);
        normalOfCenterPlane.normalize();
        
        rg_Point3D startVertex    = gate.getStartVertex()->getPoint();
        rg_REAL    distBetPcAndVs =   normalOfCenterPlane.innerProduct( startVertex )
                                    - normalOfCenterPlane.innerProduct( centerOfGB[0] );

        rg_Point3D center = startVertex - (distBetPcAndVs*normalOfCenterPlane);
        rg_REAL    radius = gateBall[0].distance( center );    
        tangentSphereOnCP[0].setSphere( center, radius );

        return LINEAR_EDGE; 
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

	    normalOfCenterPlane = vec21.crossProduct(vec31);
        normalOfCenterPlane.normalize();
    }
    //  normal of center plane for a circular edge.
    else  {
        Sphere startBall = gate.getGateBallCCW1()->getBall();

	    rg_Point3D vecBg1Bg2 = centerOfGB[1]         - centerOfGB[0];
	    rg_Point3D vecBg1Bs  = startBall.getCenter() - centerOfGB[0];

	    normalOfCenterPlane = vecBg1Bg2.crossProduct(vecBg1Bs);
        normalOfCenterPlane.normalize();

        rg_Point3D vecBsBg1 = centerOfGB[0] - startBall.getCenter();
        rg_Point3D vecBsBg2 = centerOfGB[1] - startBall.getCenter();
        rg_Point3D vecBsBg3 = centerOfGB[2] - startBall.getCenter();

        vecBsBg1.normalize();
        vecBsBg2.normalize();
        vecBsBg3.normalize();
        rg_REAL angleCg1VsCg2 = angleCCWOfTwoUnitVectors( normalOfCenterPlane, 
                                                          vecBsBg1, vecBsBg2 );
        rg_REAL angleCg1VsCg3 = angleCCWOfTwoUnitVectors( normalOfCenterPlane, 
                                                          vecBsBg1, vecBsBg3 );

        if ( angleCg1VsCg2 > angleCg1VsCg3  )  
            normalOfCenterPlane = (-normalOfCenterPlane);
    }

    /*
    rg_REAL centerPlane[4];
    centerPlane[0] = normalOfCenterPlane.getX();
    centerPlane[1] = normalOfCenterPlane.getY();
    centerPlane[2] = normalOfCenterPlane.getZ();
    centerPlane[3] = -normalOfCenterPlane.innerProduct( centerOfGB[0] );

	rg_INT numOfTSOnPc = 0;
    numOfTSOnPc = computeSphereWithItsCenterOnPlaneAndTangentTo3SpheresOutside( 
                         centerPlane, gateBall[0], gateBall[1], gateBall[2], tangentSphereOnCP);

    rg_REAL distToPc = normalOfCenterPlane.innerProduct(tangentSphereOnCP[0].getCenter()) + centerPlane[3];
    rg_REAL distToB1 = tangentSphereOnCP[0].distance( gateBall[0] );
    rg_REAL distToB2 = tangentSphereOnCP[0].distance( gateBall[1] );
    rg_REAL distToB3 = tangentSphereOnCP[0].distance( gateBall[2] );

    //  project gate balls into center plane.
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
    numOfCC = computeCircleTangentTo3CirclesOutside( transformedCircle1, transformedCircle2, 
                                                     transformedCircle3, tangentCircle);   
    rg_Point3D center[2];
    for(rg_INT i=0; i<2; i++)
        center[i] = invMatrix*tangentCircle[i].getCenterPt();

    if ( rg_TRUE )
        int aaa = 111;

    // The gate balls define a parabolic, or hyperbolic edge.
    if(numOfTSOnPc == 1) 
	{      
        ///////////////////////////////////////////////////////////////////////////
        //  2. find the center of tangent circle.
        //
        //  tangentSphereOnCP[0]
        
        ///////////////////////////////////////////////////////////////////////////
        //  3. find the axis point.

        rg_Point3D startVertex = gate.getStartVertex()->getPoint();
        if ( tangentSphereOnCP[0].getCenter() != startVertex )  {

            rg_REAL distBetVsAndPc = normalOfCenterPlane.innerProduct( startVertex ) + centerPlane[3];
            axisPoint = startVertex - (distBetVsAndPc*normalOfCenterPlane);
        }
        else  {
            
            numCenterMinTSEQVs++;
            Sphere anotherTangentSphere[2];
            rg_INT numTS = computeSphereWithGivenRadiusAndTangentTo3Spheres( 
                                                2.*tangentSphereOnCP[0].getRadius(),
                                                gateBall, anotherTangentSphere );

            rg_REAL distBetVsAndPc = normalOfCenterPlane.innerProduct( anotherTangentSphere[0].getCenter() ) + centerPlane[3];
            axisPoint = anotherTangentSphere[0].getCenter() - (distBetVsAndPc*normalOfCenterPlane);
        }


        ///////////////////////////////////////////////////////////////////////////
        //  4. find the normal of edge plane.

        rg_Point3D vecPaCt = tangentSphereOnCP[0].getCenter() - axisPoint;
        vecPaCt.normalize();

        normalOfEdgePlane = vecPaCt.crossProduct( normalOfCenterPlane );
        normalOfEdgePlane.normalize();


        return PARABOLIC_OR_HYPERBOLIC_EDGE;	
	}
	else if(numOfTSOnPc == 2) // The gate defines an elliptic edge.
	{
        ///////////////////////////////////////////////////////////////////////////
        //
        //  2. find the center of tangent circle.
        //
        rg_Vector3D vecTS1Bc1 = centerOfGB[0] - tangentSphereOnCP[0].getCenter();
        rg_Vector3D vecTS1Bc2 = centerOfGB[1] - tangentSphereOnCP[0].getCenter();
        rg_Vector3D vecTS1Bc3 = centerOfGB[2] - tangentSphereOnCP[0].getCenter();
        
        rg_REAL angleGc1Tc1Gc2 = angleCCW( normalOfCenterPlane, vecTS1Bc1, vecTS1Bc2 );
        rg_REAL angleGc1Tc1Gc3 = angleCCW( normalOfCenterPlane, vecTS1Bc1, vecTS1Bc3 );

        if ( rg_GT( angleGc1Tc1Gc2, angleGc1Tc1Gc3 ) )  {
            Sphere tempSphere = tangentSphereOnCP[0];
            tangentSphereOnCP[0] = tangentSphereOnCP[1];
            tangentSphereOnCP[1] = tempSphere;
        }
        
        ///////////////////////////////////////////////////////////////////////////
        //
        //  3. find the axis point.
        //
        axisPoint = ( tangentSphereOnCP[0].getCenter() + tangentSphereOnCP[1].getCenter() ) / 2.0;

        ///////////////////////////////////////////////////////////////////////////
        //
        //  4. find the normal of edge plane.
        //
        rg_Point3D vecPaCt = tangentSphereOnCP[0].getCenter() - axisPoint;

        normalOfEdgePlane = vecPaCt.crossProduct( normalOfCenterPlane );
        normalOfEdgePlane.normalize();


        return CIRCULAR_OR_ELLIPTIC_EDGE;	
	}
    else 
    {
        return NON_EDGE;	
    }    
    */


    
    ///////////////////////////////////////////////////////////////////////////
    //
    //  find the intersection between three gate balls and center pland
    //    for the center of tangent circle, axis point, and normal of edge plane.
    //

    //  project gate balls into center plane.
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
    numOfCC = computeCircleTangentTo3CirclesOutside( transformedCircle1, transformedCircle2, 
                                                     transformedCircle3, tangentCircle);   

    rg_Point3D centerOfTangentCircle;

    // The gate balls define a parabolic, or hyperbolic edge.
    if(numOfCC == 1)  {      
        ///////////////////////////////////////////////////////////////////////////
        //  2. find the center of tangent circle.

        centerOfTangentCircle = invMatrix*tangentCircle[0].getCenterPt();

        tangentSphereOnCP[0].setSphere( centerOfTangentCircle, tangentCircle[0].getRadius() );

        
        ///////////////////////////////////////////////////////////////////////////
        //  3. find the axis point.

        rg_Point3D startVertex = gate.getStartVertex()->getPoint();
        if ( centerOfTangentCircle != startVertex )  {
            rg_REAL distanceVsAndCP = normalOfCenterPlane.innerProduct( startVertex )
                                      - normalOfCenterPlane.innerProduct( centerOfTangentCircle );
            axisPoint = startVertex - (distanceVsAndCP*normalOfCenterPlane);
        }
        else  {
            Sphere anotherTangentSphere[2];
            rg_INT numTS = computeSphereWithGivenRadiusAndTangentTo3Spheres( 
                                                2.*tangentSphereOnCP[0].getRadius(),
                                                gateBall, anotherTangentSphere );

            rg_REAL distanceViAndCP = normalOfCenterPlane.innerProduct( anotherTangentSphere[0].getCenter() )
                                      - normalOfCenterPlane.innerProduct( centerOfTangentCircle );

            axisPoint = anotherTangentSphere[0].getCenter() - (distanceViAndCP*normalOfCenterPlane);
            
            //numCenterMinTSEQVs++;
            //BallGenerator* ptrGateBall[3] = { (BallGenerator*) gate.getGateBallCCW1()->getGenerator(), 
            //                                  (BallGenerator*) gate.getGateBallCCW2()->getGenerator(), 
            //                                  (BallGenerator*) gate.getGateBallCCW3()->getGenerator() }; 
            //BallGenerator* ptrStartBall = (BallGenerator*) gate.getStartBall()->getGenerator();
            //
            //rg_Point3D ptOnEdge;
            //rg_dNode<BallGenerator>* startNode = m_GlobalGeneratorList.getFirstpNode();
            //rg_dNode<BallGenerator>* currNode  = startNode;
            //Sphere tangentSphere[2];
            //rg_INT numTangentSphere;
            //do  {
            //    BallGenerator* currBall = currNode->getpEntity();
            //
            //    if (   currBall == ptrGateBall[0] || currBall == ptrGateBall[1] 
            //        || currBall == ptrGateBall[2] || currBall == ptrStartBall  )
            //        continue;
            //
            //    numTangentSphere = computeSphereTangentTo4SpheresOutside(
            //                                gateBall[0], gateBall[1], gateBall[2], currBall->getBall(),
            //                                tangentSphere);
            //    if ( numTangentSphere != 0 ) {
            //        ptOnEdge = tangentSphere[0].getCenter();
            //        break;
            //    }
            //
            //    currNode = currNode->getNext();
            //} while ( currNode != startNode );
            //
            //rg_REAL distanceViAndCP = normalOfCenterPlane.innerProduct( ptOnEdge )
            //                          - normalOfCenterPlane.innerProduct( centerOfTangentCircle );
            //
            //axisPoint = ptOnEdge - (distanceViAndCP*normalOfCenterPlane);
            
        }


        ///////////////////////////////////////////////////////////////////////////
        //  4. find the normal of edge plane.

        rg_Point3D vecPaCt = centerOfTangentCircle - axisPoint;
        vecPaCt.normalize();

        normalOfEdgePlane = vecPaCt.crossProduct( normalOfCenterPlane );
        normalOfEdgePlane.normalize();


        return PARABOLIC_OR_HYPERBOLIC_EDGE;	
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

            tangentSphereOnCP[0].setSphere( centerOfTangentCircle, tangentCircle[0].getRadius() );
            tangentSphereOnCP[1].setSphere( theOtherCenter, tangentCircle[1].getRadius() );

        }
        else {
            centerOfTangentCircle = invMatrix*tangentCircle[1].getCenterPt();
            theOtherCenter        = invMatrix*tangentCircle[0].getCenterPt();

            tangentSphereOnCP[0].setSphere( centerOfTangentCircle, tangentCircle[1].getRadius() );
            tangentSphereOnCP[1].setSphere( theOtherCenter, tangentCircle[0].getRadius() );
        }


        ///////////////////////////////////////////////////////////////////////////
        //
        //  3. find the axis point.
        //
        axisPoint = ( centerOfTangentCircle + theOtherCenter ) / 2.0;

        ///////////////////////////////////////////////////////////////////////////
        //
        //  4. find the normal of edge plane.
        //
        rg_Point3D vecPaCt = centerOfTangentCircle - axisPoint;

        normalOfEdgePlane = vecPaCt.crossProduct( normalOfCenterPlane );
        normalOfEdgePlane.normalize();


        return CIRCULAR_OR_ELLIPTIC_EDGE;	
	}
    else 
    {
        return NON_EDGE;	
    }    
    
    
}



rg_REAL rg_SphereSetVoronoiDiagram::angleCCW(const rg_Point2D& vector1, const rg_Point2D& vector2)
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

rg_REAL rg_SphereSetVoronoiDiagram::angleCCW(const rg_Point3D& normal, const rg_Point3D& vector1, const rg_Point3D& vector2)
{
    rg_Point3D vec1 = vector1.getUnitVector();
	rg_Point3D vec2 = vector2.getUnitVector();
	rg_REAL angle = vec1.innerProduct(vec2);

	if ( angle > 1.0 ) 
		angle = 0.0;
	else if ( angle < -1.0 )
		angle = 2.0;
    else 
		angle = 1.0 - angle;

    
    //rg_REAL angle = vec1.angle( vec2 );

    rg_Point3D cross = vec1.crossProduct( vec2 );
    cross.normalize();

    rg_REAL orientation = normal.innerProduct( cross );
	if( rg_POS( orientation ) )  //less than PI (cross product)
		return angle;
	else
		return (4.0 - angle);
    /*
    rg_Point3D vec1 = vector1.getUnitVector();
	rg_Point3D vec2 = vector2.getUnitVector();

	rg_REAL cosineTheta = 0.0;
	rg_REAL angle       = 0.0;

	cosineTheta = vec1.innerProduct(vec2);
	if ( cosineTheta > 1.0 ) 
		angle = 0.0;
	else if ( cosineTheta < -1.0 )
		angle = rg_PI;
	else
		angle = acos(cosineTheta);

    
    //rg_REAL angle = vec1.angle( vec2 );

    rg_Point3D cross = vec1.crossProduct( vec2 );
    cross.normalize();

    rg_REAL orientation = normal.innerProduct( cross );
	if( rg_POS( orientation ) )  //less than PI (cross product)
		return angle;
	else
		return 2.*rg_PI - angle;

    */
}



//rg_REAL rg_SphereSetVoronoiDiagram::angleCCWBetweenNormalizedTwoVectors(const rg_Point3D& normal, 
rg_REAL rg_SphereSetVoronoiDiagram::angleCCWOfTwoUnitVectors(const rg_Point3D& normal, 
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

    /*
	rg_REAL cosineTheta = vector1.innerProduct(vector2);
	rg_REAL angle       = 0.0;

	if ( cosineTheta > 1.0 ) 
		angle = 0.0;
	else if ( cosineTheta < -1.0 )
		angle = rg_PI;
	else
		angle = acos(cosineTheta);

    
    //rg_REAL angle = vec1.angle( vec2 );

    rg_Point3D cross = vector1.crossProduct( vector2 );
    cross.normalize();

    rg_REAL orientation = normal.innerProduct( cross );
	if( rg_POS( orientation ) )  //less than PI (cross product)
		return angle;
	else
		return 2.*rg_PI - angle;
    */
}




VDVertex* rg_SphereSetVoronoiDiagram::exploreEdgeToBeMadeByCurrentGate( const Gate           &gate,
                                                                              VIDIC          &theVIDIC,
                                                                              rg_dList<Gate> &theGateStack)
{
    VDVertex* ptrEndVertex = rg_NULL;
    
    BallGenerator* gateBall[3]   = { gate.getGateBallCCW1(), gate.getGateBallCCW2(), gate.getGateBallCCW3() };
    BallGenerator* endBall       = gate.getEndBall();
    Sphere         tangentSphere = gate.getEmptyTangentSphereAtEndVertex();
    rg_FLAG        edgeType      = gate.getEdgeType();


    if ( endBall != rg_NULL )  {
        //  Current gate has a end vertex. That means that the edge of the current gate is bounded.

        //  find the end vertex configuration in VIDIC.
        VIDICItem* existingVertexConfig = theVIDIC.findVIDICItem( gateBall[0], gateBall[1], gateBall[2], endBall,
                                                                  tangentSphere.getCenter() );

        if ( existingVertexConfig != rg_NULL )  {
            //  end vertex exist in VIDIC.
            ptrEndVertex = existingVertexConfig->getVertex();

            
            if ( ON_DEBUG )
            {
                foutEVDS << gate.getStartBall()->getID()   << "\t"
                         << gateBall[0]->getID()           << "\t"
                         << gateBall[1]->getID()           << "\t"
                         << gateBall[2]->getID()           << "\t"
                         << endBall->getID()               << "\t";

                if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
                    foutEVDS << "Para" << "\t";
                else if ( edgeType == LINEAR_EDGE )
                    foutEVDS << "Line" << "\t";
                else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
                    foutEVDS << "Elli" << "\t";
                else 
                    foutEVDS << "Non" << "\t";

                foutEVDS << "PV"                           << "\t"
                         << gate.getStartVertex()->getID() << "\t"
                         << ptrEndVertex->getID()          << "\t"
                         << ptrEndVertex->getPoint().getX() << "\t"
                         << ptrEndVertex->getPoint().getY() << "\t"
                         << ptrEndVertex->getPoint().getZ() << "\t"
                         << ptrEndVertex->getRadiusOfTangentSphere() << endl;
            }
            

            //  for debug. by Youngsong Cho 2005.03.09
            if ( ptrEndVertex->getIncidentEdge( 3 ) != rg_NULL )  {
                m_status = ERROR_INFINITE_LOOP;
                return rg_NULL;
            }

            //  remove redundant gate in Gate-Stack
            rg_INT i=0;
			for (i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)  {
                Gate* gateOfEndVertex = ptrEndVertex->getGate( i );

                if ( gateOfEndVertex == rg_NULL )  continue;

                if ( gateOfEndVertex->isTheSameGate(gate, ptrEndVertex) )
                {
                    theGateStack.kill( gateOfEndVertex->getNodePosInStack() );
                    ptrEndVertex->setGate( i, rg_NULL );
                    break;
                }
            }           
            
            //  if end vertex has 4 incident edges by edge to be made by current gate,
            //  remove the end vertex configuration in VIDIC.
            if ( ptrEndVertex->getIncidentEdge( 2 ) != rg_NULL )  {
                //theVIDIC.removeVIDICItem( *existingVertexConfig );
                ptrEndVertex->removeAllGates();
            }            
        }
        else  {
            ptrEndVertex = makeVoronoiVertex( rg_FALSE, tangentSphere );
    
            //  make and insert 3 gates of new voronoi vertex into Gate-Stack.
            Gate currentGates[ 3 ] 
                     = { Gate( endBall, gateBall[0], gateBall[1], gateBall[2], ptrEndVertex ),
                         Gate( endBall, gateBall[1], gateBall[2], gateBall[0], ptrEndVertex ),
                         Gate( endBall, gateBall[2], gateBall[0], gateBall[1], ptrEndVertex ) };

            rg_INT i=0;
			for (i=0; i<3; i++)  {
                Gate* ptrGate = theGateStack.addTail( currentGates[i] );
                ptrGate->setNodePosInStack( theGateStack.getTail() );

                ptrEndVertex->setGate( i, ptrGate );
            }   
            ptrEndVertex->setGate( 3, rg_NULL );

            //  insert vertex configuration of initial voronoi vertex into VIDIC.
            theVIDIC.insertVIDICItem( gateBall[0], gateBall[1], gateBall[2], endBall, ptrEndVertex );

            if ( ON_DEBUG )
            {
                foutEVDS << gate.getStartBall()->getID()   << "\t"
                         << gateBall[0]->getID()           << "\t"
                         << gateBall[1]->getID()           << "\t"
                         << gateBall[2]->getID()           << "\t"
                         << endBall->getID()               << "\t";

                if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
                    foutEVDS << "Para" << "\t";
                else if ( edgeType == LINEAR_EDGE )
                    foutEVDS << "Line" << "\t";
                else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
                    foutEVDS << "Elli" << "\t";
                else 
                    foutEVDS << "Non" << "\t";

                foutEVDS << "NV"                           << "\t"
                         << gate.getStartVertex()->getID() << "\t"
                         << ptrEndVertex->getID()          << "\t"
                         << ptrEndVertex->getPoint().getX() << "\t"
                         << ptrEndVertex->getPoint().getY() << "\t"
                         << ptrEndVertex->getPoint().getZ() << "\t"
                         << ptrEndVertex->getRadiusOfTangentSphere() << endl;
            }
        }
    }
    else  {  // => (endBall == rg_NULL)
        //  Current gate does not have a end vertex. That means that the edge of the current gate is unbounded.
        ptrEndVertex = makeVoronoiVertex( rg_TRUE, Sphere( DBL_MAX, DBL_MAX, DBL_MAX, 0.0) );

        if ( ON_DEBUG )
        {
            foutEVDS << gate.getStartBall()->getID()   << "\t"
                     << gateBall[0]->getID()           << "\t"
                     << gateBall[1]->getID()           << "\t"
                     << gateBall[2]->getID()           << "\t"
                     << "\t";                           

            if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
                foutEVDS << "Para" << "\t";
            else if ( edgeType == LINEAR_EDGE )
                foutEVDS << "Line" << "\t";
            else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
                foutEVDS << "Elli" << "\t";
            else 
                foutEVDS << "Non" << "\t";

            foutEVDS << "IV"                           << "\t"
                     << gate.getStartVertex()->getID() << "\t"
                     << ptrEndVertex->getID()          << endl;
        }
    }


    
    //  if start vertex of current gate has 4 incident edges by edge to make by current gate,
    //  remove the end vertex configuration in VIDIC.
    if ( gate.getStartVertex()->getIncidentEdge( 2 ) != rg_NULL )  {
        //theVIDIC.removeVIDICItem( gateBall[0],
        //                          gateBall[1],
        //                          gateBall[2],
        //                          gate.getStartBall(),
        //                          gate.getStartVertex() );

        gate.getStartVertex()->removeAllGates();
    }
    

    return ptrEndVertex;
}




VDEdge* rg_SphereSetVoronoiDiagram::setLocalTopologyByEndVertex( const Gate &gate, VDVertex* endVertex )
{
    VDEdge*         ptrEdge         = rg_NULL;
    VDPartialEdge** ptrPartialEdges = rg_NULL;


    ptrEdge         = makeVoronoiEdge( gate.getStartVertex(), endVertex );
    ptrEdge->setEdgeType( gate.getEdgeType() );

    ptrPartialEdges = make3VoronoiPartialEdges( ptrEdge );

    makeVoronoiFaceAndItsOuterLoop( gate, ptrPartialEdges );

    delete [] ptrPartialEdges;

    return ptrEdge;
}





void rg_SphereSetVoronoiDiagram::setTopology()
{
    defineOrderOfPartialEdgesInAllLoops();
    setTopologyOnInfinity();

    rg_dList< VDCell >* globalCellList = getGlobalCellList();
    VDCell* infinityCell               = globalCellList->getLastpEntity();
    
    rg_dList<VDFace*>* boundingFaces = infinityCell->getBoungindFaces();
    VDFace* currFace = rg_NULL;
    boundingFaces->reset4Loop();
    while ( boundingFaces->setNext4Loop() )  {
        currFace = boundingFaces->getEntity();

        if ( currFace->getRightCell() == infinityCell ) 
            currFace->getLeftCell()->isBounded( rg_FALSE );
        else
            currFace->getRightCell()->isBounded( rg_FALSE );
    }
}




void rg_SphereSetVoronoiDiagram::defineOrderOfPartialEdgesInAllLoops()
{
    VDLoop* currLoop = rg_NULL;

    m_GlobalLoopList.reset4Loop();
    while( m_GlobalLoopList.setNext4Loop() )
    {
        currLoop = m_GlobalLoopList.getpEntity();

        if ( ON_DEBUG )
        {
            foutEVDS << "Loop: " << currLoop->getID() << endl;

            if ( currLoop->getPartialEdge()->getOriginalEdge()->getEdgeType() != EDGE_WITHOUT_VERTICES)  {
                VDPartialEdge* startPartEdge = currLoop->getPartialEdge();
                VDPartialEdge* currPartEdge  = startPartEdge;

                do 
                {
                    if ( currPartEdge->isRightOrientationInLoop() == rg_TRUE )
                    {
                        foutEVDS << "\t" << "v" << currPartEdge->getOriginalEdge()->getStartVertex()->getID();
                        if ( currPartEdge->getOriginalEdge()->getStartVertex()->isOnInfinity() )
                            foutEVDS << " (inf)";
                        foutEVDS << "->" << "v" << currPartEdge->getOriginalEdge()->getEndVertex()->getID();
                        if ( currPartEdge->getOriginalEdge()->getEndVertex()->isOnInfinity() )
                            foutEVDS << " (inf)";

                        foutEVDS << "\t(+)\t" << currPartEdge->getID() << endl;
                    }
                    else //if ( currPartEdge->isRightOrientationInLoop() == rg_FALSE )
                    {
                        foutEVDS << "\t" << "v" << currPartEdge->getOriginalEdge()->getEndVertex()->getID();
                        if ( currPartEdge->getOriginalEdge()->getEndVertex()->isOnInfinity() )
                            foutEVDS << " (inf)";

                        foutEVDS << "->" << "v" << currPartEdge->getOriginalEdge()->getStartVertex()->getID();
                        if ( currPartEdge->getOriginalEdge()->getStartVertex()->isOnInfinity() )
                            foutEVDS << " (inf)";

                        foutEVDS << "\t(-)" << currPartEdge->getID() << endl;
                    }

            
                    currPartEdge = currPartEdge->getNextPartialEdgeInLoop();
                } while ( currPartEdge != startPartEdge );
            }
        }

        if ( currLoop->getPartialEdge()->getOriginalEdge()->getEdgeType() == EDGE_WITHOUT_VERTICES)  
            continue;

        rg_dList<VDPartialEdge*> multipleOuterLoop;
        //if ( currLoop->reorderOfPartialEdgesInLoop( multipleOuterLoop ) == rg_FALSE )
        if ( currLoop->constructLoopCycle( multipleOuterLoop ) == rg_FALSE )
        {
            if ( multipleOuterLoop.getSize() > 0 )  {
                //makeAdditiveOuterLoopOfFace( currLoop->getFace(), multipleOuterLoop );
                makeVDFaceOfAdditiveOuterLoop( currLoop->getFace(), multipleOuterLoop );
            }
        }

        //  for DEBUG
        if ( ON_DEBUG )
        {
            foutEVDS << "Loop: " << currLoop->getID() << endl;

            if ( currLoop->getPartialEdge()->getOriginalEdge()->getEdgeType() != EDGE_WITHOUT_VERTICES)  {
                VDPartialEdge* startPartEdge = currLoop->getPartialEdge();
                VDPartialEdge* currPartEdge  = startPartEdge;

                do 
                {
                    if ( currPartEdge->isRightOrientationInLoop() == rg_TRUE )
                    {
                        foutEVDS << "\t" << "v" << currPartEdge->getOriginalEdge()->getStartVertex()->getID();
                        if ( currPartEdge->getOriginalEdge()->getStartVertex()->isOnInfinity() )
                            foutEVDS << " (inf)";
                        foutEVDS << "->" << "v" << currPartEdge->getOriginalEdge()->getEndVertex()->getID();
                        if ( currPartEdge->getOriginalEdge()->getEndVertex()->isOnInfinity() )
                            foutEVDS << " (inf)";

                        foutEVDS << "\t(+)\t" << currPartEdge->getID() << endl;
                    }
                    else //if ( currPartEdge->isRightOrientationInLoop() == rg_FALSE )
                    {
                        foutEVDS << "\t" << "v" << currPartEdge->getOriginalEdge()->getEndVertex()->getID();
                        if ( currPartEdge->getOriginalEdge()->getEndVertex()->isOnInfinity() )
                            foutEVDS << " (inf)";

                        foutEVDS << "->" << "v" << currPartEdge->getOriginalEdge()->getStartVertex()->getID();
                        if ( currPartEdge->getOriginalEdge()->getStartVertex()->isOnInfinity() )
                            foutEVDS << " (inf)";

                        foutEVDS << "\t(-)" << currPartEdge->getID() << endl;
                    }

            
                    currPartEdge = currPartEdge->getNextPartialEdgeInLoop();
                } while ( currPartEdge != startPartEdge );
            }
        }
    }
}

void rg_SphereSetVoronoiDiagram::makeAdditiveOuterLoopOfFace(
                                     VDFace*                   ptrFace, 
                                     rg_dList<VDPartialEdge*>& multipleOuterLoop )
{
    VDLoop* ptrLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrFace, rg_TRUE ) );

    ptrFace->addLoop( ptrLoop );

    VDPartialEdge* currPrEdge = rg_NULL;

    multipleOuterLoop.reset4Loop();
    while( multipleOuterLoop.setNext4Loop() )
    {
        currPrEdge = multipleOuterLoop.getEntity();

        currPrEdge->setNextPartialEdgeInLoop( currPrEdge );
        currPrEdge->setPreviousPartialEdgeInLoop( currPrEdge );

        ptrLoop->addPartialEdgeToLoop( currPrEdge );
    }    
}


void rg_SphereSetVoronoiDiagram::makeVDFaceOfAdditiveOuterLoop(VDFace* ptrFace, rg_dList<VDPartialEdge*>& multipleOuterLoop)
{
    VDCell* rightCell = ptrFace->getRightCell();
    VDCell* leftCell  = ptrFace->getLeftCell();
    VDFace* ptrNewFace = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), rightCell, leftCell) );
    VDLoop* ptrNewLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrNewFace, rg_TRUE ) );

    rightCell->addBoundingFace( ptrNewFace );
    leftCell->addBoundingFace( ptrNewFace );

    ptrNewFace->addLoop( ptrNewLoop );

    VDPartialEdge* currPrEdge = rg_NULL;

    multipleOuterLoop.reset4Loop();
    while( multipleOuterLoop.setNext4Loop() )
    {
        currPrEdge = multipleOuterLoop.getEntity();

        currPrEdge->setNextPartialEdgeInLoop( currPrEdge );
        currPrEdge->setPreviousPartialEdgeInLoop( currPrEdge );

        ptrNewLoop->addPartialEdgeToLoop( currPrEdge );
    }    
}


void rg_SphereSetVoronoiDiagram::setTopologyOnInfinity()
{
    VDCell* virtualCell = m_GlobalCellList.addTail( VDCell( -1 ) );


    //  find all outer loops of unbounded faces
    rg_dList< VDLoop* > outerLoopsOfUnboundedFaces;

    VDLoop* currLoop = rg_NULL;
    VDFace* currFace = rg_NULL;
    m_GlobalFaceList.reset4Loop();
    while( m_GlobalFaceList.setNext4Loop() )
    {
        currFace = m_GlobalFaceList.getpEntity();

        rg_dList<VDLoop*>* loops = currFace->getLoops();
        loops->reset4Loop();
        while( loops->setNext4Loop() )
        {
            currLoop = loops->getEntity();
            if ( !currLoop->isOuterLoop() )
                continue;

            if ( !currLoop->isClosedLoop() )
                outerLoopsOfUnboundedFaces.add( currLoop );
        }
    }

    
    //  outer loop -> make one edge on infinity.
    //  an edge    -> make three partial edges.
    outerLoopsOfUnboundedFaces.reset4Loop();
    while( outerLoopsOfUnboundedFaces.setNext4Loop() )
    {
        currLoop = outerLoopsOfUnboundedFaces.getEntity();

        VDPartialEdge* firstPrEdge = currLoop->getPartialEdge();
        VDPartialEdge* currPrEdge  = firstPrEdge;

        do
        {
            if ( currPrEdge->getOriginalEdge()->isBoundedEdge() == rg_FALSE )
                break;

            currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
        } while ( currPrEdge != firstPrEdge );



        //  define the Voronoi edge which connects two infinite vertex in a unbounded face.
        VDVertex* startVertex = rg_NULL;
        VDVertex* endVertex   = rg_NULL;
        VDEdge* ptrEdge = rg_NULL;
        if ( currPrEdge->getNextPartialEdgeInLoop()->getOriginalEdge()->isBoundedEdge() == rg_FALSE )
        {
            if ( currPrEdge->isRightOrientationInLoop() )
                startVertex = currPrEdge->getOriginalEdge()->getEndVertex();
            else
                endVertex = currPrEdge->getOriginalEdge()->getEndVertex();

            if ( currPrEdge->getNextPartialEdgeInLoop()->isRightOrientationInLoop() )
                startVertex   = currPrEdge->getNextPartialEdgeInLoop()->getOriginalEdge()->getEndVertex();
            else
                endVertex   = currPrEdge->getNextPartialEdgeInLoop()->getOriginalEdge()->getEndVertex();
        }
        else if ( currPrEdge->getPreviousPartialEdgeInLoop()->getOriginalEdge()->isBoundedEdge() == rg_FALSE )
        {
            if ( currPrEdge->getPreviousPartialEdgeInLoop()->isRightOrientationInLoop() )
                startVertex   = currPrEdge->getPreviousPartialEdgeInLoop()->getOriginalEdge()->getEndVertex();
            else
                endVertex   = currPrEdge->getPreviousPartialEdgeInLoop()->getOriginalEdge()->getEndVertex();

            if ( currPrEdge->isRightOrientationInLoop() )
                startVertex = currPrEdge->getOriginalEdge()->getEndVertex();
            else
                endVertex = currPrEdge->getOriginalEdge()->getEndVertex();
        }
        else
        {}

        VDEdge edge( m_GlobalEdgeList.getSize(), startVertex, endVertex );
        ptrEdge = m_GlobalEdgeList.addTail( edge );
        ptrEdge->setInfinity();

        startVertex->setIncidentEdge( ptrEdge );
        endVertex->setIncidentEdge( ptrEdge );

        //  define three partial edges.
        VDPartialEdge* ptrPartialEdges[3] = {rg_NULL, rg_NULL, rg_NULL};

        rg_INT i=0;
		for ( i=0; i<3; i++ )
        {
            ptrPartialEdges[i] = m_GlobalPartialedgeList.addTail( VDPartialEdge( m_GlobalPartialedgeList.getSize(), ptrEdge ) );

            ptrPartialEdges[i]->setNextPartialEdgeInLoop( ptrPartialEdges[i] );
            ptrPartialEdges[i]->setPreviousPartialEdgeInLoop( ptrPartialEdges[i] );
        }
        ptrEdge->setPartialEdge( ptrPartialEdges[0] );

        ptrPartialEdges[0]->setNextPartialEdgeInRadialCycle( ptrPartialEdges[1] );
        ptrPartialEdges[1]->setNextPartialEdgeInRadialCycle( ptrPartialEdges[2] );
        ptrPartialEdges[2]->setNextPartialEdgeInRadialCycle( ptrPartialEdges[0] );



        //  define face on infinity and its outer loop.
        VDCell* ballsForEdge[4] = { currLoop->getFace()->getLeftCell(), 
                                    currLoop->getFace()->getRightCell(), 
                                    virtualCell,  
                                    currLoop->getFace()->getLeftCell() };

        //  partial edges[0] <- ballsForEdge 0 & ballsForEdge 1
        //  partial edges[1] <- ballsForEdge 1 & ballsForEdge 2
        //  partial edges[2] <- ballsForEdge 2 & ballsForEdge 0
        VDFace* faceIncidentToInfiniteEdge[3] = { currLoop->getFace(), 
                                                  ballsForEdge[1]->findFaceToShareWith( ballsForEdge[2] ), 
                                                  ballsForEdge[2]->findFaceToShareWith( ballsForEdge[3] )};
        for ( i=0; i<3; i++ )
        {
            VDFace* ptrFace = faceIncidentToInfiniteEdge[i];
            //VDFace* ptrFace = rg_NULL;
            //if ( ballsForEdge[i]->getID() != -1 )
            //    ptrFace = ballsForEdge[i]->findFaceToShareWith( ballsForEdge[i+1] );
            //else
            //    ptrFace = ballsForEdge[i+1]->findFaceToShareWith( ballsForEdge[i] );

            VDLoop* ptrLoop = rg_NULL;

            if ( ptrFace == rg_NULL )
            {
                ptrFace = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), ballsForEdge[i+1], ballsForEdge[i] ) );
                ptrLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrFace, rg_TRUE ) );
                
                ptrFace->setInfinity();


                //  define connection between VDCell and VDFace.
                ballsForEdge[i]->addBoundingFace( ptrFace );
                ballsForEdge[i+1]->addBoundingFace( ptrFace );

                //  define connection between VDFace and VDLoop.
                ptrFace->addLoop( ptrLoop );

                //  define connection between VDLoop and VDPartialEdge.
                //ptrPartialEdges[i]->setLoop(ptrLoop);
                ptrPartialEdges[i]->isRightOrientationInLoop(rg_TRUE);
            }
            else
            {
                ptrLoop = ptrFace->getLoop( 0 );

                //  define connection between VDLoop and VDPartialEdge.
                //ptrPartialEdges[i]->setLoop(ptrLoop);

                if (    ptrFace->getLeftCell()  == ballsForEdge[i] 
                     && ptrFace->getRightCell() == ballsForEdge[i+1] )
                    ptrPartialEdges[i]->isRightOrientationInLoop(rg_TRUE);
                else
                    ptrPartialEdges[i]->isRightOrientationInLoop(rg_FALSE);
            }

            ptrLoop->addPartialEdgeOnInfinityToLoop( ptrPartialEdges[i] );
        }
    }
}

void rg_SphereSetVoronoiDiagram::setGeometry()
{
    calculateBoundingBoxForVoronoiVertices();

    setEndVertexOfInfiniteEdge();

    //setGeometryOfVoronoiEdges();
}

void rg_SphereSetVoronoiDiagram::calculateBoundingBoxForVoronoiVertices()
{
    rg_REAL x;     
    rg_REAL y;     
    rg_REAL z;     

    rg_REAL min_X = m_minPointOfBoundingBoxForBalls.getX();
    rg_REAL max_X = m_maxPointOfBoundingBoxForBalls.getX();
    rg_REAL min_Y = m_minPointOfBoundingBoxForBalls.getY();
    rg_REAL max_Y = m_maxPointOfBoundingBoxForBalls.getY();
    rg_REAL min_Z = m_minPointOfBoundingBoxForBalls.getZ();
    rg_REAL max_Z = m_maxPointOfBoundingBoxForBalls.getZ();

    VDVertex* currVertex = rg_NULL;

    m_GlobalVertexList.reset4Loop();
    while( m_GlobalVertexList.setNext4Loop() )
    {
        currVertex = m_GlobalVertexList.getpEntity();
    
        if ( currVertex->isOnInfinity() )
            continue;

        x = currVertex->getPoint().getX();
        y = currVertex->getPoint().getY();
        z = currVertex->getPoint().getZ();

        if ( x < min_X )
            min_X = x;
        else if ( x > max_X )  
            max_X = x;
        else {}

        if ( y < min_Y )
            min_Y = y;
        else if ( y > max_Y )  
            max_Y = y;       
        else {}

        if ( z < min_Z )
            min_Z = z;
        else if ( z > max_Z )  
            max_Z = z;
        else {}
    }

    rg_REAL halfLengthOfDiagonal =   (max_X - min_X)*(max_X - min_X)
                                   + (max_Y - min_Y)*(max_Y - min_Y)
                                   + (max_Z - min_Z)*(max_Z - min_Z);

    halfLengthOfDiagonal = sqrt( halfLengthOfDiagonal ) / 2.0;

    m_minPointOfBoundingBoxForVD.setPoint(min_X-halfLengthOfDiagonal, min_Y-halfLengthOfDiagonal, min_Z-halfLengthOfDiagonal);
    m_maxPointOfBoundingBoxForVD.setPoint(max_X+halfLengthOfDiagonal, max_Y+halfLengthOfDiagonal, max_Z+halfLengthOfDiagonal);
}




void rg_SphereSetVoronoiDiagram::setPropertyOfVoronoiEdges()
{
    VDEdge* currEdge = rg_NULL;

	m_GlobalEdgeList.reset4Loop();
	while( m_GlobalEdgeList.setNext4Loop() )  
    {
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isOnInfinity() )
            continue;

		VDCell* cellsInCCWtoDefineThisEdge[NUM_CELLS_TO_DEFINE_EDGE]; 
		currEdge->inquireIntoOnlyEdgeSharingCellsInCCW( cellsInCCWtoDefineThisEdge );

		Sphere s1 = ((BallGenerator*)cellsInCCWtoDefineThisEdge[0]->getGenerator())->getBall();
		Sphere s2 = ((BallGenerator*)cellsInCCWtoDefineThisEdge[1]->getGenerator())->getBall();
		Sphere s3 = ((BallGenerator*)cellsInCCWtoDefineThisEdge[2]->getGenerator())->getBall();

        //rg_REAL* planeBySphereCenters = evaluatePlanePassingThroughSphereCenters(s1, s2, s3);
        Plane    planeBySphereCenters(s1.getCenter(), s2.getCenter(), s3.getCenter() );
        rg_Point3D startPt = currEdge->getStartVertex()->getPoint();
        rg_Point3D endPt   = currEdge->getEndVertex()->getPoint();

        if ( currEdge->isBoundedEdge() )
		{
		    ///////////////////////////////////////////////////////////////////////////
		    //  define the properties of tangent spheres on bounded edge.
			//rg_REAL posOfStartVertexOnPlane =   planeBySphereCenters[0]*startPt.getX()
			//								  + planeBySphereCenters[1]*startPt.getY()
			//								  + planeBySphereCenters[2]*startPt.getZ()
			//								  + planeBySphereCenters[3];
			//rg_REAL posOfEndVertexOnPlane   =   planeBySphereCenters[0]*endPt.getX()
			//								  + planeBySphereCenters[1]*endPt.getY()
			//								  + planeBySphereCenters[2]*endPt.getZ()
			//								  + planeBySphereCenters[3];
            rg_REAL posOfStartVertexOnPlane =   planeBySphereCenters.distanceFromPoint( startPt );
			rg_REAL posOfEndVertexOnPlane   =   planeBySphereCenters.distanceFromPoint( endPt );

			if ( rg_NEG( posOfStartVertexOnPlane*posOfEndVertexOnPlane ) )
			{
				Sphere smallestTangentSphere;
				//rg_Point3D normalOfPlaneBy3SphereCenters( planeBySphereCenters[0],
				//										  planeBySphereCenters[1],
				//										  planeBySphereCenters[2]);
                rg_Point3D normalOfPlaneBy3SphereCenters = planeBySphereCenters.getNormal();

				smallestTangentSphere = evaluateSmallestTangentSphereFor3Spheres( normalOfPlaneBy3SphereCenters,
																				  s1, s2, s3);

    			currEdge->setMinPosition(MID);            
    			currEdge->setMinRadius( smallestTangentSphere.getRadius() );

				if( currEdge->getStartVertex()->getRadiusOfTangentSphere() > 
					currEdge->getEndVertex()->getRadiusOfTangentSphere() )
				{
					currEdge->setMaxPosition(START);
				}
				else
				{
					currEdge->setMaxPosition(END);
				}

			}
			else
			{
				if( currEdge->getStartVertex()->getRadiusOfTangentSphere() > 
					currEdge->getEndVertex()->getRadiusOfTangentSphere() )
				{
					currEdge->setMaxPosition(START);
					currEdge->setMinPosition(END);
					currEdge->setMinRadius(currEdge->getEndVertex()->getRadiusOfTangentSphere());
				}
				else
				{
					currEdge->setMaxPosition(END);
					currEdge->setMinPosition(START);
					currEdge->setMinRadius(currEdge->getStartVertex()->getRadiusOfTangentSphere());
				}
			}
		}
		else
		{      
            ///////////////////////////////////////////////////////////////////////////
			//  define the properties of tangent spheres on unbounded edge.

			//rg_REAL posOfStartVertexOnPlane =   planeBySphereCenters[0]*startPt.getX()
			//								  + planeBySphereCenters[1]*startPt.getY()
			//								  + planeBySphereCenters[2]*startPt.getZ()
			//								  + planeBySphereCenters[3];
            rg_REAL posOfStartVertexOnPlane =   planeBySphereCenters.distanceFromPoint( startPt );

			if ( rg_NEG(posOfStartVertexOnPlane) )
			{
				Sphere smallestTangentSphere;
				//rg_Point3D normalOfPlaneBy3SphereCenters( planeBySphereCenters[0],
				//										  planeBySphereCenters[1],
				//										  planeBySphereCenters[2]);
                rg_Point3D normalOfPlaneBy3SphereCenters = planeBySphereCenters.getNormal();

				smallestTangentSphere = evaluateSmallestTangentSphereFor3Spheres( normalOfPlaneBy3SphereCenters,
																				  s1, s2, s3);

    			currEdge->setMaxPosition(END);
    			currEdge->setMinPosition(MID);            
    			currEdge->setMinRadius( smallestTangentSphere.getRadius() );

			}
			else
			{
    			currEdge->setMaxPosition(END);
    			currEdge->setMinPosition(START);            
    			currEdge->setMinRadius(currEdge->getStartVertex()->getRadiusOfTangentSphere());
			}
        
        }
	}
}



void rg_SphereSetVoronoiDiagram::setGeometryOfVoronoiEdges()
{
    VDEdge* currEdge = rg_NULL;

	m_GlobalEdgeList.reset4Loop();
	while( m_GlobalEdgeList.setNext4Loop() )  
    {
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isOnInfinity() )
            continue;

		VDCell* cellsInCCWtoDefineThisEdge[NUM_CELLS_TO_DEFINE_EDGE]; 
		currEdge->inquireIntoOnlyEdgeSharingCellsInCCW( cellsInCCWtoDefineThisEdge );

        Sphere gateBall[3] = { ((BallGenerator*)cellsInCCWtoDefineThisEdge[0]->getGenerator())->getBall(),
		                       ((BallGenerator*)cellsInCCWtoDefineThisEdge[1]->getGenerator())->getBall(),
		                       ((BallGenerator*)cellsInCCWtoDefineThisEdge[2]->getGenerator())->getBall() };

    	rg_RBzCurve3D* curve = rg_NULL; 
        if ( currEdge->getStartVertex() != rg_NULL )  {    
            
            rg_Point3D startPt = currEdge->getStartVertex()->getPoint();
            rg_Point3D endPt   = currEdge->getEndVertex()->getPoint();
            
            evaluateVoronoiEdge(curve, startPt, endPt, gateBall[0], gateBall[1], gateBall[2]);

            if ( currEdge->getEdgeType() == PARABOLIC_OR_HYPERBOLIC_EDGE )  {
				//  for partially tangible edges.
    		    rg_dList <rg_Point3D> listOfPointsOnCurve; 
				int i=0;
				for(i=0; i<=NUM_POINTS_ON_EDGE; i++)  { 
					listOfPointsOnCurve.add(curve->evaluatePt( (rg_REAL)i/NUM_POINTS_ON_EDGE ));
				}

				rg_INT      numOfPointOnCurve = listOfPointsOnCurve.getSize();
				rg_Point3D* pointsOnCurve     = listOfPointsOnCurve.getArray();


				currEdge->setEdgeEquation( curve );
				currEdge->setPointsOnEdge( numOfPointOnCurve, pointsOnCurve );
            }
            else if ( currEdge->getEdgeType() == LINEAR_EDGE )  {
				rg_INT      numOfPointOnCurve = 2;
				rg_Point3D* pointsOnCurve     = new rg_Point3D[numOfPointOnCurve];

				pointsOnCurve[0] = startPt;
				pointsOnCurve[1] = endPt;
                /*
				//  for partially tangible edges.
    		    rg_dList <rg_Point3D> listOfPointsOnCurve; 
				for(int i=0; i<=NUM_POINTS_ON_EDGE; i++)  { 
					listOfPointsOnCurve.add(curve->evaluatePt( (rg_REAL)i/NUM_POINTS_ON_EDGE ));
				}

				rg_INT      numOfPointOnCurve = listOfPointsOnCurve.getSize();
				rg_Point3D* pointsOnCurve     = listOfPointsOnCurve.getArray();
                */

				currEdge->setEdgeEquation( curve );
				currEdge->setPointsOnEdge( numOfPointOnCurve, pointsOnCurve );
            }
            else if ( currEdge->getEdgeType() == CIRCULAR_OR_ELLIPTIC_EDGE )  {
				//  for partially tangible edges.
                rg_INT numPointsOnEllipticEdge = NUM_POINTS_ON_EDGE*2;
    		    rg_dList <rg_Point3D> listOfPointsOnCurve; 
				int i=0;
				for(i=0; i<=numPointsOnEllipticEdge; i++)  { 
					listOfPointsOnCurve.add(curve->evaluatePt( (rg_REAL)i/numPointsOnEllipticEdge ));
				}

				rg_INT      numOfPointOnCurve = listOfPointsOnCurve.getSize();
				rg_Point3D* pointsOnCurve     = listOfPointsOnCurve.getArray();


				currEdge->setEdgeEquation( curve );
				currEdge->setPointsOnEdge( numOfPointOnCurve, pointsOnCurve );
            }
            else {
            }
        }
        else {
            Sphere s0, s1;
			Sphere* pSpheres = NULL;
			rg_INT numOfSpheres = 0;

			numOfSpheres = makeCircumcircleOnCenterPlane(gateBall[0], gateBall[1], gateBall[2],s0,s1);
			s0.setRadius(s0.getRadius()*0.9);
			evaluate_Voronoi_Vertex(&gateBall[0],&gateBall[1],&gateBall[2],&s0,pSpheres);
			if(pSpheres==NULL)
			{
				s0.setRadius(0);
				evaluate_Voronoi_Vertex(&gateBall[0],&gateBall[1],&gateBall[2],&s0,pSpheres);
			}

			rg_Point3D tVecStart = getTangentVectorOfEdgeAtVertex(s0.getCenter(),gateBall[0], gateBall[1], gateBall[2]);
			rg_Point3D tVecEnd = getTangentVectorOfEdgeAtVertex(pSpheres[0].getCenter(),gateBall[0], gateBall[1], gateBall[2]);
			
			rg_RBzCurve3D rbzCrv(2);
			makeConicByProjection(s0.getCenter(),tVecStart,pSpheres[0].getCenter(),tVecEnd,pSpheres[1].getCenter(),rbzCrv);

			curve = new rg_RBzCurve3D(rbzCrv);

	        rg_dList <rg_Point3D> listOfPointsOnCurve; 
            int i=0;
			for(i=0; i<20; i++)  { 
		        listOfPointsOnCurve.add(curve->evaluatePt( (float)i/20 ));
	        }
			curve->setWeight(1, -curve->getWeight(1));
			for(i=20; i>=0; i--)
			{
		        listOfPointsOnCurve.add(curve->evaluatePt( (float)i/20 ));
	        }
			curve->setWeight(1, -curve->getWeight(1));


	        rg_INT      numOfPointOnCurve = listOfPointsOnCurve.getSize();
            rg_Point3D* pointsOnCurve     = listOfPointsOnCurve.getArray();

            currEdge->setEdgeEquation( curve );
            currEdge->setPointsOnEdge( numOfPointOnCurve, pointsOnCurve );
        }
	}
}

void rg_SphereSetVoronoiDiagram::setEndVertexOfInfiniteEdge()
{
    VDEdge* currEdge = rg_NULL;

	m_GlobalEdgeList.reset4Loop();
	while( m_GlobalEdgeList.setNext4Loop() )  
    {
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isBoundedEdge() || currEdge->isOnInfinity() )  
            continue;

        VDCell* cellsInCCWSharingEdge[3]; 
        currEdge->inquireIntoOnlyEdgeSharingCellsInCCW( cellsInCCWSharingEdge );

        if ( currEdge->getStartVertex()->isOnInfinity() )
        {
            VDCell* tempCell = cellsInCCWSharingEdge[0];
            cellsInCCWSharingEdge[0] = cellsInCCWSharingEdge[2];
            cellsInCCWSharingEdge[2] = tempCell;
        }

        Sphere balls[3];
        int i=0;
		for ( i=0; i<3; i++)
            balls[i] = ((BallGenerator*)cellsInCCWSharingEdge[i]->getGenerator())->getBall();

        //rg_REAL* planeBySphereCenters = evaluatePlanePassingThroughSphereCenters(balls[0], balls[1], balls[2]);
        Plane planeBySphereCenters(balls[0].getCenter(), balls[1].getCenter(), balls[2].getCenter());


        sort3SpheresByRadiusInDescendingPowers( balls );
		
        rg_Point3D translateBall2( balls[2].getCenter() );
        rg_REAL    smallestRadius = balls[2].getRadius();
        // initilizing variables for calculating Vertex
        // 3 spheres are shrinked and translated by minimum sphere
        for ( i=0; i<3; i++)
        {
            balls[i].setCenter( balls[i].getCenter() - translateBall2 );
            balls[i].setRadius( balls[i].getRadius() - smallestRadius );
        }
        rg_Point3D minPoint( m_minPointOfBoundingBoxForVD - translateBall2 );
        rg_Point3D maxPoint( m_maxPointOfBoundingBoxForVD - translateBall2 );

        rg_dList<rg_Point3D> candidatePointList;

        evaluateEndVertexOfInfiniteEdge(minPoint, maxPoint, balls, candidatePointList);

        rg_Point3D point;
        candidatePointList.reset4Loop();
	    while( candidatePointList.setNext4Loop() )  
        {
		    point = candidatePointList.getEntity();

            point = point + translateBall2;

			//rg_REAL posOfPointOnPlane =   planeBySphereCenters[0]*point.getX()
			//						      + planeBySphereCenters[1]*point.getY()
            //                            + planeBySphereCenters[2]*point.getZ()
            //                            + planeBySphereCenters[3];
            rg_REAL posOfPointOnPlane =   planeBySphereCenters.distanceFromPoint( point );

			if ( rg_POS( posOfPointOnPlane ) )
			{
                if ( currEdge->getEndVertex()->isOnInfinity() )
                    currEdge->getEndVertex()->setPoint(point);
                else 
                    currEdge->getStartVertex()->setPoint(point);
                break;
            }
        }
    }
}

void rg_SphereSetVoronoiDiagram::defineEndVertexOfInfiniteEdge(VDEdge*           infiniteEdge, 
                                                               Sphere*           ballsInCCW,
                                                               const rg_Point3D& minPtOfBox,
                                                               const rg_Point3D& maxPtOfBox)
{
    if ( infiniteEdge->isBoundedEdge() )  
        return;

    //rg_REAL* planeBySphereCenters = evaluatePlanePassingThroughSphereCenters(ballsInCCW[0], ballsInCCW[1], ballsInCCW[2]);
    Plane planeBySphereCenters(ballsInCCW[0].getCenter(), ballsInCCW[1].getCenter(), ballsInCCW[2].getCenter());


    sort3SpheresByRadiusInDescendingPowers( ballsInCCW );
	
    rg_Point3D translateBall2( ballsInCCW[2].getCenter() );
    rg_REAL    smallestRadius = ballsInCCW[2].getRadius();
    // initilizing variables for calculating Vertex
    // 3 spheres are shrinked and translated by minimum sphere
    int i=0;
	for (i=0; i<3; i++)
    {
        ballsInCCW[i].setCenter( ballsInCCW[i].getCenter() - translateBall2 );
        ballsInCCW[i].setRadius( ballsInCCW[i].getRadius() - smallestRadius );
    }
    rg_Point3D minPoint( minPtOfBox - translateBall2 );
    rg_Point3D maxPoint( maxPtOfBox - translateBall2 );

    rg_dList<rg_Point3D> candidatePointList;

    evaluateEndVertexOfInfiniteEdge(minPoint, maxPoint, ballsInCCW, candidatePointList);

    rg_Point3D point;
    candidatePointList.reset4Loop();
	while( candidatePointList.setNext4Loop() )  
    {
		point = candidatePointList.getEntity();

        point = point + translateBall2;

		//rg_REAL posOfPointOnPlane =   planeBySphereCenters[0]*point.getX()
		//							+ planeBySphereCenters[1]*point.getY()
  //                                  + planeBySphereCenters[2]*point.getZ()
  //                                  + planeBySphereCenters[3];
        rg_REAL posOfPointOnPlane =   planeBySphereCenters.distanceFromPoint( point );

		if ( rg_POS( posOfPointOnPlane ) )
		{
            infiniteEdge->getEndVertex()->setPoint(point);
            break;
        }
    }
}

void rg_SphereSetVoronoiDiagram::defindGeometryOfVoronoiEdge(VDEdge* currEdge, Sphere* ballsInCCW)
{
	// define bounded edge equation.
    rg_RBzCurve3D* curve = rg_NULL; 

    rg_Point3D startPt = currEdge->getStartVertex()->getPoint();
    rg_Point3D endPt   = currEdge->getEndVertex()->getPoint();


    if ( LINEAR_TYPE_EDGE == evaluateVoronoiEdge(curve, startPt, endPt, ballsInCCW[0], ballsInCCW[1], ballsInCCW[2]) )
    {
		rg_INT      numOfPointOnCurve = 2;
		rg_Point3D* pointsOnCurve     = new rg_Point3D[numOfPointOnCurve];

		pointsOnCurve[0] = startPt;
		pointsOnCurve[1] = endPt;

		currEdge->setEdgeEquation( curve );
		currEdge->setPointsOnEdge( numOfPointOnCurve, pointsOnCurve );      
    }
    else
    {
		rg_dList <rg_Point3D> listOfPointsOnCurve; 
		int i=0;
		for(i=0; i<=NUM_POINTS_ON_EDGE; i++)  { 
			listOfPointsOnCurve.add(curve->evaluatePt( (rg_REAL)i/NUM_POINTS_ON_EDGE ));
		}

		rg_INT      numOfPointOnCurve = listOfPointsOnCurve.getSize();
		rg_Point3D* pointsOnCurve     = listOfPointsOnCurve.getArray();

		currEdge->setEdgeEquation( curve );
		currEdge->setPointsOnEdge( numOfPointOnCurve, pointsOnCurve );
    }
}

void rg_SphereSetVoronoiDiagram::evaluateEndVertexOfInfiniteEdge(const rg_Point3D& minPoint,
                                                                 const rg_Point3D& maxPoint,
                                                                 Sphere* balls,
                                                                 rg_dList<rg_Point3D>& candidatePointList)
{
    //  (x - x1)^2 + (y - y1)^2 + (z - z1)^2 = (r + r1)^2
    //  (x - x2)^2 + (y - y2)^2 + (z - z2)^2 = (r + r2)^2
    //  x^2        + y^2        + z^2        = r^2
    //   x = xL or xU, y = yL or yU, z = zL or zU

    rg_REAL x1 = balls[0].getCenter().getX();
    rg_REAL y1 = balls[0].getCenter().getY();
    rg_REAL z1 = balls[0].getCenter().getZ();
    rg_REAL r1 = balls[0].getRadius();
    
    rg_REAL x2 = balls[1].getCenter().getX();
    rg_REAL y2 = balls[1].getCenter().getY();
    rg_REAL z2 = balls[1].getCenter().getZ();
    rg_REAL r2 = balls[1].getRadius();

    rg_REAL e1 = (x1*x1 + y1*y1 + z1*z1 - r1*r1)/2.;
    rg_REAL e2 = (x2*x2 + y2*y2 + z2*z2 - r2*r2)/2.;

    rg_REAL a_xy = x1*y2 - x2*y1;
    rg_REAL a_yz = y1*z2 - y2*z1;
    rg_REAL a_xz = x1*z2 - x2*z1;

    rg_REAL b_xr = x2*r1 - x1*r2;
    rg_REAL b_yr = y2*r1 - y1*r2;
    rg_REAL b_zr = z2*r1 - z1*r2;

    //////////////////////////////////////////////////////////
    rg_REAL xL = minPoint.getX();        
    rg_REAL yL = minPoint.getY();        
    rg_REAL zL = minPoint.getZ();   
    
    rg_REAL xU = maxPoint.getX();        
    rg_REAL yU = maxPoint.getY();        
    rg_REAL zU = maxPoint.getZ();       
    
    rg_REAL c_xe = 0.0;
    rg_REAL c_ye = 0.0;
    rg_REAL c_ze = 0.0;

    rg_REAL x = 0.0;
    rg_REAL y = 0.0;
    rg_REAL z = 0.0;

    //  xL  //////////////////////////////////////////////////
    rg_REAL xL_roots[2];
    
    c_ye = y2*(x1*xL - e1) - y1*(x2*xL - e2);
    c_ze = z2*(x1*xL - e1) - z1*(x2*xL - e2);

    if ( evaluateInfiniteVertexOnOneBoundingFace(xL_roots, a_yz, b_yr, b_zr, c_ye, c_ze, xL) )
    {
        int i=0;
		for ( i=0; i<2; i++ )
        {
            if ( xL_roots[i] == -DBL_MAX )
                continue;

            y = -(b_zr*xL_roots[i] + c_ze)/a_yz;
            z =  (b_yr*xL_roots[i] + c_ye)/a_yz;

            if ( rg_GE(y, yL) && rg_LE(y, yU) && rg_GE(z, zL) && rg_LE(z, zU) )
                candidatePointList.addTail( rg_Point3D( xL, y, z) );
        }
    }

    //  xU  //////////////////////////////////////////////////
    rg_REAL xU_roots[2];
    
    c_ye = y2*(x1*xU - e1) - y1*(x2*xU - e2);
    c_ze = z2*(x1*xU - e1) - z1*(x2*xU - e2);

    if ( evaluateInfiniteVertexOnOneBoundingFace(xU_roots, a_yz, b_yr, b_zr, c_ye, c_ze, xU) )
    {
		int i=0;
        for ( i=0; i<2; i++ )
        {
            if ( xU_roots[i] == -DBL_MAX )
                continue;

            y = -(b_zr*xU_roots[i] + c_ze)/a_yz;
            z =  (b_yr*xU_roots[i] + c_ye)/a_yz;

            if ( rg_GE(y, yL) && rg_LE(y, yU) && rg_GE(z, zL) && rg_LE(z, zU) )
                candidatePointList.addTail( rg_Point3D( xU, y, z) );
        }
    }
    //  yL  //////////////////////////////////////////////////
    rg_REAL yL_roots[2];

    c_xe = x2*(y1*yL - e1) - x1*(y2*yL - e2);
    c_ze = z2*(y1*yL - e1) - z1*(y2*yL - e2);

    if ( evaluateInfiniteVertexOnOneBoundingFace(yL_roots, a_xz, b_xr, b_zr, c_xe, c_ze, yL) )
    {
		int i=0;
        for ( i=0; i<2; i++ )
        {
            if ( yL_roots[i] == -DBL_MAX )
                continue;

            x = -(b_zr*yL_roots[i] + c_ze)/a_xz;
            z =  (b_xr*yL_roots[i] + c_xe)/a_xz;

            if ( rg_GE(x, xL) && rg_LE(x, xU) && rg_GE(z, zL) && rg_LE(z, zU) )
                candidatePointList.addTail( rg_Point3D( x, yL, z) );
        }
    }    
    
    //  yU  //////////////////////////////////////////////////
    rg_REAL yU_roots[2];

    c_xe = x2*(y1*yU - e1) - x1*(y2*yU - e2);
    c_ze = z2*(y1*yU - e1) - z1*(y2*yU - e2);

    if ( evaluateInfiniteVertexOnOneBoundingFace(yU_roots, a_xz, b_xr, b_zr, c_xe, c_ze, yU) )
    {
		int i=0;
        for ( i=0; i<2; i++ )
        {
            if ( yU_roots[i] == -DBL_MAX )
                continue;

            x = -(b_zr*yU_roots[i] + c_ze)/a_xz;
            z =  (b_xr*yU_roots[i] + c_xe)/a_xz;

            if ( rg_GE(x, xL) && rg_LE(x, xU) && rg_GE(z, zL) && rg_LE(z, zU) )
                candidatePointList.addTail( rg_Point3D( x, yU, z) );
        }
    }    
    
    //  zL  //////////////////////////////////////////////////
    rg_REAL zL_roots[2];

    c_xe = x2*(z1*zL - e1) - x1*(z2*zL - e2);
    c_ye = y2*(z1*zL - e1) - y1*(z2*zL - e2);

    if ( evaluateInfiniteVertexOnOneBoundingFace(zL_roots, a_xy, b_xr, b_yr, c_xe, c_ye, zL) )
    {
		int i=0;
        for ( i=0; i<2; i++ )
        {
            if ( zL_roots[i] == -DBL_MAX )
                continue;

            x = -(b_yr*zL_roots[i] + c_ye)/a_xy;
            y =  (b_xr*zL_roots[i] + c_xe)/a_xy;

            if ( rg_GE(x, xL) && rg_LE(x, xU) && rg_GE(y, yL) && rg_LE(y, yU) )
                candidatePointList.addTail( rg_Point3D( x, y, zL) );
        }
    }    
    
    //  zL  //////////////////////////////////////////////////
    rg_REAL zU_roots[2];

    c_xe = x2*(z1*zU - e1) - x1*(z2*zU - e2);
    c_ye = y2*(z1*zU - e1) - y1*(z2*zU - e2);


    if ( evaluateInfiniteVertexOnOneBoundingFace(zU_roots, a_xy, b_xr, b_yr, c_xe, c_ye, zU) )
    {
		int i=0;
        for ( i=0; i<2; i++ )
        {
            if ( zU_roots[i] == -DBL_MAX )
                continue;

            x = -(b_yr*zU_roots[i] + c_ye)/a_xy;
            y =  (b_xr*zU_roots[i] + c_xe)/a_xy;

            if ( rg_GE(x, xL) && rg_LE(x, xU) && rg_GE(y, yL) && rg_LE(y, yU) )
                candidatePointList.addTail( rg_Point3D( x, y, zU) );
        }
    }    
    
}

rg_INT rg_SphereSetVoronoiDiagram::evaluateInfiniteVertexOnOneBoundingFace(rg_REAL* roots,
                                                                           const rg_REAL& a_ij, 
                                                                           const rg_REAL& b_ir, const rg_REAL& b_jr, 
                                                                           const rg_REAL& c_ie, const rg_REAL& c_je,
                                                                           const rg_REAL& boundingValue)
{
    rg_REAL a = b_ir*b_ir + b_jr*b_jr - a_ij*a_ij;
    rg_REAL b = 2*(b_ir*c_ie + b_jr*c_je);
    rg_REAL c = a_ij*a_ij*boundingValue*boundingValue + c_ie*c_ie + c_je*c_je;

    rg_REAL discriminant = b*b - 4.*a*c;
    if ( a == 0 || rg_NEG( discriminant ) )
    {
        return 0;
    }
    else
    {
        int numValidRoot = 2;
        roots[0] = ( -b + sqrt( discriminant ) )/ (2.*a);
        roots[1] = ( -b - sqrt( discriminant ) )/ (2.*a);

		int i=0;
		for ( i=0; i<2; i++ )
        {
            if ( rg_NEG( roots[i] ) )
            {
                roots[i] = -DBL_MAX;
                numValidRoot--;
            }
        }

        return numValidRoot;
    }
}

void rg_SphereSetVoronoiDiagram::setAttribute()
{
    setBoundednessOfFacesAndCells();
}

void rg_SphereSetVoronoiDiagram::setBoundednessOfFacesAndCells()
{
    VDEdge* currEdge = rg_NULL;

	m_GlobalEdgeList.reset4Loop();
	while( m_GlobalEdgeList.setNext4Loop() )  
    {
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isBoundedEdge() )  
            continue;

        VDPartialEdge* startPartEdge = currEdge->getPartialEdge();
        VDPartialEdge* currPartEdge  = startPartEdge;
        
        do 
        {
            VDFace* currFace = currPartEdge->getLoop()->getFace();

            currFace->isBounded(rg_FALSE);

            currFace->getLeftCell()->isBounded(rg_FALSE);
            currFace->getRightCell()->isBounded(rg_FALSE);

            currPartEdge = currPartEdge->getNextPartialEdgeInRadialCycle();
        } while ( currPartEdge != startPartEdge );
    }
}

VDEdge*  rg_SphereSetVoronoiDiagram::makeVoronoiEdge( VDVertex* startVertex, VDVertex* endVertex )
{
    //  insert edge by start and end vertices of gate into GEL(Global Edge List).
    VDEdge edge( m_GlobalEdgeList.getSize(), startVertex, endVertex );

    VDEdge* ptrEdge = m_GlobalEdgeList.addTail( edge );

    //  connect edge to vertex as incident edge of vertex.
    startVertex->setIncidentEdge( ptrEdge );
    endVertex->setIncidentEdge( ptrEdge );

    return ptrEdge;
}


VDPartialEdge** rg_SphereSetVoronoiDiagram::make3VoronoiPartialEdges( VDEdge* theEdge )
{    
    VDPartialEdge** ptrPartialEdges = new VDPartialEdge*[3];

    //  insert 3 partial edges of the edge into GPL(Global Partial edge List).

    for ( rg_INT i=0; i<3; i++ )
    {
        ptrPartialEdges[i] = m_GlobalPartialedgeList.addTail( VDPartialEdge( m_GlobalPartialedgeList.getSize(), theEdge ) );

        ptrPartialEdges[i]->setNextPartialEdgeInLoop( ptrPartialEdges[i] );
        ptrPartialEdges[i]->setPreviousPartialEdgeInLoop( ptrPartialEdges[i] );
    }

    //  connect one partial edge to the edge.
    theEdge->setPartialEdge( ptrPartialEdges[0] );

    //  connect partial edges in Radial Cycle.
    ptrPartialEdges[0]->setNextPartialEdgeInRadialCycle( ptrPartialEdges[1] );
    ptrPartialEdges[1]->setNextPartialEdgeInRadialCycle( ptrPartialEdges[2] );
    ptrPartialEdges[2]->setNextPartialEdgeInRadialCycle( ptrPartialEdges[0] );

    return ptrPartialEdges;
}



void rg_SphereSetVoronoiDiagram::makeVoronoiFaceAndItsOuterLoop( BallGenerator** ballToDefineEdge , VDPartialEdge** partialEdges )
{
    VDCell* gateCell[4] = { ballToDefineEdge[0]->getCell(), ballToDefineEdge[1]->getCell(), 
                            ballToDefineEdge[2]->getCell(), ballToDefineEdge[0]->getCell() };

    
    //  partial edges[0] <- gate ball 0 & gate ball 1
    //  partial edges[1] <- gate ball 1 & gate ball 2
    //  partial edges[2] <- gate ball 2 & gate ball 0
    for ( rg_INT i=0; i<3; i++ )
    {
        VDFace* ptrFace = gateCell[i]->findFaceToShareWith( gateCell[i+1] );
        VDLoop* ptrLoop = rg_NULL;

        if ( ptrFace == rg_NULL )
        {
            ptrFace = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), gateCell[i+1], gateCell[i] ) );
            ptrLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrFace, rg_TRUE ) );

            //  define connection between VDCell and VDFace.
            gateCell[i]->addBoundingFace( ptrFace );
            gateCell[i+1]->addBoundingFace( ptrFace );

            //  define connection between VDFace and VDLoop.
            ptrFace->addLoop( ptrLoop );

            //  define connection between VDLoop and VDPartialEdge.
            //partialEdges[i]->setLoop(ptrLoop);
            partialEdges[i]->isRightOrientationInLoop(rg_TRUE);
        }
        else
        {
            ptrLoop = ptrFace->getLoop( 0 );

            //  define connection between VDLoop and VDPartialEdge.
            //partialEdges[i]->setLoop(ptrLoop);

            if (    ptrFace->getLeftCell() == gateCell[i] 
                 && ptrFace->getRightCell() == gateCell[i+1] )
                partialEdges[i]->isRightOrientationInLoop(rg_TRUE);
            else
                partialEdges[i]->isRightOrientationInLoop(rg_FALSE);
        }

        ptrLoop->addPartialEdgeToLoop( partialEdges[i] );

    }
}



void rg_SphereSetVoronoiDiagram::makeVoronoiFaceAndItsOuterLoop( const Gate &gate, VDPartialEdge** partialEdges )
{
    VDCell* gateCell[4] = { gate.getGateBallCCW1()->getCell(), 
                            gate.getGateBallCCW2()->getCell(), 
                            gate.getGateBallCCW3()->getCell(), 
                            gate.getGateBallCCW1()->getCell() };

    //  partial edges[0] <- gate ball 0 & gate ball 1
    //  partial edges[1] <- gate ball 1 & gate ball 2
    //  partial edges[2] <- gate ball 2 & gate ball 0

    for ( rg_INT i=0; i<3; i++ )
    {
        VDFace* ptrFace = gateCell[i]->findFaceToShareWith( gateCell[i+1] );
        VDLoop* ptrLoop = rg_NULL;

        if ( ptrFace == rg_NULL )
        {
            ptrFace = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), gateCell[i+1], gateCell[i] ) );
            ptrLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrFace, rg_TRUE ) );

            //  define connection between VDCell and VDFace.
            gateCell[i]->addBoundingFace( ptrFace );
            gateCell[i+1]->addBoundingFace( ptrFace );

            //  define connection between VDFace and VDLoop.
            ptrFace->addLoop( ptrLoop );

            //  define connection between VDLoop and VDPartialEdge.
            //partialEdges[i]->setLoop(ptrLoop);
            partialEdges[i]->isRightOrientationInLoop(rg_TRUE);
        }
        else
        {
            ptrLoop = ptrFace->getLoop( 0 );

            //  define connection between VDLoop and VDPartialEdge.
            //partialEdges[i]->setLoop(ptrLoop);

            if (    ptrFace->getLeftCell() == gateCell[i] 
                 && ptrFace->getRightCell() == gateCell[i+1] )
                partialEdges[i]->isRightOrientationInLoop(rg_TRUE);
            else
                partialEdges[i]->isRightOrientationInLoop(rg_FALSE);
        }

        ptrLoop->addPartialEdgeToLoop( partialEdges[i] );

    }
}



rg_INT rg_SphereSetVoronoiDiagram::constructSphereVoronoiDiagramByAcceleratedEdgeTracing()
{
    //perturbBallGenerators();


    m_status = rg_TRUE;

    /*
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //
    rg_INT numBalls = m_GlobalGeneratorList.getSize();
    m_frequencyOfCandidateBall = new rg_INT[ numBalls ];
    for ( rg_INT i=0; i<numBalls; i++ )
        m_frequencyOfCandidateBall[i] = 0;
    //
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //*/


    start = clock();


    setBucket( m_sizeOfBucketElement );

    //  construct BIG world.
    rg_FLAG bConstruction = constructBigWorld();


    start_in = clock();

    
    //  construct SMALL worlds.
    if ( m_status )  {
        if ( bConstruction == rg_TRUE )
            constructSmallWorlds();
        else  {
            makeVirtualTwoBigBrothers();
            constructSmallWorlds();
        }
    }


    //  computingTime[15]: construct small worlds.
    finish_in = clock();
    computingTime[15] = finish_in - start_in;




    //m_bucket.removeBucketTable();
    

    finish = clock();
    computingTime[16] = finish - start;


    return m_status;
}



rg_FLAG rg_SphereSetVoronoiDiagram::constructBigWorld()
{
    ///////////////////////////////////////////////////////////////////////////
    //
    //  find initial voronoi vertex.
    //
    start_in = clock();

    rg_INT initialSizeOfVIDIC = 101;

    VIDIC            theVIDIC( initialSizeOfVIDIC );
    rg_dList< Gate > theGateStack;
  
    if ( findInitialVertexInBigWorld( theVIDIC, theGateStack ) == rg_FALSE )  
        return rg_FALSE;

    //  computingTime[0]:  initializing EVD(S)
    finish_in = clock();
    computingTime[0] = finish_in - start_in;
    m_bStartCheck = rg_TRUE;

    ///////////////////////////////////////////////////////////////////////////
    //
    //  trace Voronoi edges in BIG world.
    //
    while ( theGateStack.getSize() > 0 ) {

        ///////////////////////////////////////////////////////////////////////////
        //
        start_in2 = clock();

        Gate currGate = theGateStack.popFront();

        findEndVertexBasedOnAngularDistanceInAcceleration( currGate ); 

        //  computingTime[1]:  finding end vertex
        finish_in = clock();
        computingTime[1] += finish_in - start_in2;



        ///////////////////////////////////////////////////////////////////////////
        //
        start_in = clock();

        VDVertex* ptrEndVertex = exploreEdgeToBeMadeByCurrentGate( currGate, theVIDIC, theGateStack );

        if ( m_status != rg_TRUE )
            break;

        //  computingTime[12]: exploring edge
        finish_in = clock();
        computingTime[12] += finish_in - start_in;


        
        ///////////////////////////////////////////////////////////////////////////
        //
        start_in = clock();

        setLocalTopologyByEndVertex( currGate, ptrEndVertex ); 

        //  computingTime[13]: setting local topology
        finish_in = clock();
        computingTime[13] += finish_in - start_in;
    }



    ///////////////////////////////////////////////////////////////////////////
    //
    //  set global topology of BIG world.
    //
    start_in = clock();

    if ( m_status == rg_TRUE ) {
        setTopology();
    }

    //  computingTime[13]: setting global topology
    finish_in = clock();
    computingTime[14] += finish_in - start_in;
    m_bStartCheck = rg_FALSE;


    return rg_TRUE;

}



rg_FLAG rg_SphereSetVoronoiDiagram::findInitialVertexInBigWorld( VIDIC& theVIDIC, rg_dList<Gate>& theGateStack )
{
    rg_FLAG bSuccess = rg_FALSE;


    //  1. set the configuration of virtual initial vertex
    rg_Point3D virtualVertexCoord(1.0, 1.0, 1.0);
    rg_REAL    virtualVertexRadius  = 0.73205080756887719;    

    rg_Point3D centerOfVirtualBall[4];
    centerOfVirtualBall[0].setPoint(2.0, 0.0, 0.0);
    centerOfVirtualBall[1].setPoint(0.0, 2.0, 0.0);
    centerOfVirtualBall[2].setPoint(0.0, 0.0, 2.0);
    centerOfVirtualBall[3].setPoint(0.0, 0.0, 0.0);

    rg_Point3D oldMinPointOfBoundingBoxForBalls = m_minPointOfBoundingBoxForBalls;
    rg_Point3D sizeOfBoundingBox = m_maxPointOfBoundingBoxForBalls - m_minPointOfBoundingBoxForBalls;
    rg_REAL    maxSize = sizeOfBoundingBox.getX();
    if ( maxSize < sizeOfBoundingBox.getY() )
        maxSize = sizeOfBoundingBox.getY();
    if ( maxSize < sizeOfBoundingBox.getZ() )
        maxSize = sizeOfBoundingBox.getZ();
    rg_Point3D translation(-maxSize, -maxSize, -maxSize);

    translation = m_minPointOfBoundingBoxForBalls + translation;
    m_minPointOfBoundingBoxForBalls = translation;
    
    virtualVertexCoord = virtualVertexCoord + translation;


    BallGenerator* ptrVirtualBall[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
    VDCell*        ptrVirtualCell[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
    rg_INT i=0;
	for (i=0; i<4; i++)  {
        centerOfVirtualBall[i] = centerOfVirtualBall[i] + translation;

        ptrVirtualBall[i] = m_GlobalGeneratorList.addHead( BallGenerator(centerOfVirtualBall[i], 1.0) );
        ptrVirtualCell[i] = m_GlobalCellList.addHead( VDCell( -2-i ) );
        connectCellAndGenerator( ptrVirtualCell[i], ptrVirtualBall[i] );
    }

    //  2. find valid initial vertex
    VIDIC            tempVIDIC( 31 );
    rg_dList< Gate > tempGateStack;

    defineInitialVoronoiVertex( ptrVirtualBall, Sphere(virtualVertexCoord, virtualVertexRadius), tempVIDIC, tempGateStack);

    Sphere initialTangentSphere;
    Gate   initialGates[NUM_DEFINING_CELLS_OF_ONE_VERTEX];

    //setBucket(m_sizeOfBucketElement);

    while ( tempGateStack.getSize() > 0 )
    {
        // 1.0 pop a gate;
        Gate currGate = tempGateStack.popBack();

        if ( currGate.getStartBall()->getID() == -1 
             && currGate.getGateBallCCW1()->getID() == 3 
             && currGate.getGateBallCCW2()->getID() == 8 
             && currGate.getGateBallCCW3()->getID() == 547 )
             int aaa = 1;


        // 1.1 GET THE END VERTEX
        findEndVertexBasedOnAngularDistanceInAcceleration( currGate ); 
        BallGenerator* endBall = currGate.getEndBall();
        Sphere  tangentSphereByEndBall( currGate.getEmptyTangentSphereAtEndVertex() );

        if ( endBall != rg_NULL )
        {
            if (    currGate.getGateBallCCW1()->getCell()->getID() > -2 && currGate.getGateBallCCW2()->getCell()->getID() > -2
                 && currGate.getGateBallCCW3()->getCell()->getID() > -2 && endBall->getCell()->getID() > -2 )
            {
                BallGenerator* gateBall[3] = { currGate.getGateBallCCW1(), currGate.getGateBallCCW2(), currGate.getGateBallCCW3()};

                initialTangentSphere = tangentSphereByEndBall;
                initialGates[0].setGate( gateBall[2], gateBall[1], gateBall[0], endBall, rg_NULL );
                initialGates[1].setGate( endBall, gateBall[0], gateBall[1], gateBall[2], rg_NULL );
                initialGates[2].setGate( endBall, gateBall[1], gateBall[2], gateBall[0], rg_NULL );
                initialGates[3].setGate( endBall, gateBall[2], gateBall[0], gateBall[1], rg_NULL );

                bSuccess = rg_TRUE;


                if ( ON_DEBUG )
                {
                    foutEVDS << "*** FIND INITIAL VERTEX" << endl;
                    foutEVDS << currGate.getStartBall()->getID()   << "\t"
                             << gateBall[0]->getID()           << "\t"
                             << gateBall[1]->getID()           << "\t"
                             << gateBall[2]->getID()           << "\t"
                             << endBall->getID()               << "\t";

                    rg_FLAG edgeType = currGate.getEdgeType();
                    if ( edgeType == PARABOLIC_OR_HYPERBOLIC_EDGE )
                        foutEVDS << "Para" << "\t";
                    else if ( edgeType == LINEAR_EDGE )
                        foutEVDS << "Line" << "\t";
                    else if ( edgeType == CIRCULAR_OR_ELLIPTIC_EDGE )
                        foutEVDS << "Elli" << "\t";
                    else 
                        foutEVDS << "Non" << "\t";

                    foutEVDS << endl;

                }
                
                break;
            }
        }

        // 1.2 EXPLORE GATE
        VDVertex* ptrEndVertex = exploreEdgeToBeMadeByCurrentGate( currGate, tempVIDIC, tempGateStack );

        makeVoronoiEdge( currGate.getStartVertex(), ptrEndVertex );
    }
    //m_bucket.removeBucketTable();

    //  3. remove virtual initial vertex and its topology
    m_minPointOfBoundingBoxForBalls = oldMinPointOfBoundingBoxForBalls;

    for ( i=0; i<4; i++)  {
        m_GlobalCellList.popFront();
        m_GlobalGeneratorList.popFront();
    }

    VDCell* currCell = rg_NULL;
	m_GlobalCellList.reset4Loop();
	while( m_GlobalCellList.setNext4Loop() )
	{
		currCell = m_GlobalCellList.getpEntity();
        currCell->getBoungindFaces()->removeAll();
    }
    m_GlobalEdgeList.removeAll();
    m_GlobalVertexList.removeAll();


    if ( bSuccess )  {    
        //  4. set valid initial vertex.
        //  insert initial voronoi vertex into GVL(Global Vertex List).
        VDVertex* ptrVertex = makeVoronoiVertex( rg_FALSE, initialTangentSphere );

        //  make and insert 4 gates of initial voronoi vertex into Gate-Stack.
        Gate* ptrGate = rg_NULL;        
        for ( i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        {
            ptrGate = theGateStack.addTail( initialGates[i] );

            ptrGate->setStartVertex( ptrVertex );
            ptrGate->setNodePosInStack( theGateStack.getTail() );

            ptrVertex->setGate( i, ptrGate );
        }   

        //  insert vertex configuration of initial voronoi vertex into VIDIC.
        theVIDIC.insertVIDICItem( VIDICItem( initialGates[0].getGateBallCCW1(),
                                             initialGates[0].getGateBallCCW2(),
                                             initialGates[0].getGateBallCCW3(),
                                             initialGates[0].getStartBall(),
                                             ptrVertex) 
                                 );
    }

    return bSuccess;
}



void rg_SphereSetVoronoiDiagram::findEmptyTangentSphereWithMinAngleInParabolicAndHyperbolicEdge(
                                                                  pair<BallGenerator*, Sphere>& endBall,
                                                                  rg_REAL&                      minAngularDistance,
                                                                  rg_dList<BallGenerator*>&     candidateBalls,
                                                                  Sphere*                       gateBall,
                                                                  BallGenerator*                startBall,   
                                                                  const rg_Point3D&             coordStartVertex,
					                                              const rg_Point3D&             normalOfEdgePlane,
					                                              const rg_Vector3D&            vecPaVs,
					                                              const rg_Point3D&             axisPoint,
                                                                  const rg_REAL&                maxValidAngle)
{
    rg_INT  numTS = 0;
    Sphere  tangentSphere[2];

    BallGenerator* currBall = rg_NULL;

    candidateBalls.reset4Loop();
	while( candidateBalls.setNext4Loop() )	{
	    currBall = candidateBalls.getEntity();

        numTS = computeSphereTangentTo4SpheresOutside( gateBall[0], gateBall[1], gateBall[2], 
                                                       currBall->getBall(), tangentSphere);


        rg_FLAG bChecked = rg_FALSE;
        for ( rg_INT i=0; i<numTS; i++ )  
        {
            if ( currBall == startBall && tangentSphere[i].getCenter() == coordStartVertex )
                continue;
    
            rg_Point3D vecViPa         = tangentSphere[i].getCenter() - axisPoint;
            rg_REAL    angularDistance = angleCCWOfTwoUnitVectors( normalOfEdgePlane, 
                                                                   vecPaVs, vecViPa.getUnitVector() );

            if ( angularDistance < maxValidAngle ) {
                
                /*
                //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
                //
                if ( m_bStartCheck )  {
                    if ( m_isFrontFeasibleRegion )  {
                        if ( bChecked == rg_FALSE )  {
                            m_feasibleBalls.add( currBall );
                            bChecked = rg_TRUE;
                        }
                    }
                }
                //
                //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
                //*/

                if ( angularDistance < minAngularDistance) {
                    endBall.first    = currBall;
                    endBall.second   = tangentSphere[i];
                    minAngularDistance = angularDistance;                
                }
            }
        }
    }
}



void rg_SphereSetVoronoiDiagram::findEmptyTangentSphereWithMinAngleInLinearEdge(
                                                                  pair<BallGenerator*, Sphere>& endBall,
                                                                  rg_REAL&                      minAngularDistance,
                                                                  rg_dList<BallGenerator*>&     candidateBalls,
                                                                  Sphere*                       gateBall,
                                                                  BallGenerator*                startBall,   
                                                                  const rg_Point3D&             coordStartVertex,
					                                              const rg_Point3D&             normalOfCenterPlane)
{
    rg_INT  numTS = 0;
    Sphere  tangentSphere[2];

    BallGenerator* currBall = rg_NULL;

    candidateBalls.reset4Loop();
	while( candidateBalls.setNext4Loop() )	{
	    currBall = candidateBalls.getEntity();

        numTS = computeSphereTangentTo4SpheresOutside( gateBall[0], gateBall[1], gateBall[2], 
                                                       currBall->getBall(), tangentSphere);



        rg_FLAG bChecked = rg_FALSE;
        for ( rg_INT i=0; i<numTS; i++ )  
        {
            if ( currBall == startBall && tangentSphere[i].getCenter() == coordStartVertex )
                continue;
    
            rg_Point3D vecViVs     = tangentSphere[i].getCenter() - coordStartVertex;           
            rg_REAL    orientation = normalOfCenterPlane.innerProduct( vecViVs.getUnitVector() );
            if ( !rg_POS( orientation ) )
                continue;

            /*
            //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
            if ( m_bStartCheck )  {
                if ( m_isFrontFeasibleRegion )  {
                    if ( bChecked == rg_FALSE )  {
                        m_feasibleBalls.add( currBall );
                        bChecked = rg_TRUE;
                    }
                }
            }
            //
            //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
            //*/

            rg_REAL    angularDistance = vecViVs.magnitude();
            if ( angularDistance < minAngularDistance) {
                endBall.first    = currBall;
                endBall.second   = tangentSphere[i];
                minAngularDistance = angularDistance;                
            }
        }
    }
}



void rg_SphereSetVoronoiDiagram::findEmptyTangentSphereWithMinAngleInEllipticAndCircularEdge(
                                                                  pair<BallGenerator*, Sphere>& endBall,
                                                                  rg_REAL&                      minAngularDistance,
                                                                  rg_dList<BallGenerator*>&     candidateBalls,
                                                                  Sphere*                       gateBall,
                                                                  BallGenerator*                startBall,   
                                                                  const rg_Point3D&             coordStartVertex,
					                                              const rg_Point3D&             normalOfEdgePlane,
					                                              const rg_Vector3D&            vecPaVs,
					                                              const rg_Point3D&             axisPoint)
{
    rg_INT  numTS = 0;
    Sphere  tangentSphere[2];

    BallGenerator* currBall = rg_NULL;

    candidateBalls.reset4Loop();
	while( candidateBalls.setNext4Loop() )	{
	    currBall = candidateBalls.getEntity();

        numTS = computeSphereTangentTo4SpheresOutside( gateBall[0], gateBall[1], gateBall[2], 
                                                       currBall->getBall(), tangentSphere);


        /*
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //
        if ( numTS > 0)
            m_numCandidateBallsOfElliptic[0]++;
        //
        //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
        //*/

        for ( rg_INT i=0; i<numTS; i++ )  
        {
            if ( currBall == startBall && tangentSphere[i].getCenter() == coordStartVertex )
                continue;
    
            rg_Point3D vecViPa         = tangentSphere[i].getCenter() - axisPoint;
            rg_REAL    angularDistance = angleCCWOfTwoUnitVectors( normalOfEdgePlane, 
                                                                   vecPaVs, vecViPa.getUnitVector() );

            if ( angularDistance < minAngularDistance) {
                endBall.first    = currBall;
                endBall.second   = tangentSphere[i];
                minAngularDistance = angularDistance;                
            }
        }
    }
}



void rg_SphereSetVoronoiDiagram::makeVirtualTwoBigBrothers()
{
    //setBucket(m_sizeOfBucketElement);


    rg_Point3D sizeOfBox   = m_maxPointOfBoundingBoxForBalls - m_minPointOfBoundingBoxForBalls;
    rg_Point3D centerOfBox = ( m_minPointOfBoundingBoxForBalls + m_maxPointOfBoundingBoxForBalls ) / 2.0;

    rg_REAL maxLength = -1.0;
    if ( maxLength < sizeOfBox.getX() )
        maxLength = sizeOfBox.getX();
    if ( maxLength < sizeOfBox.getY() )
        maxLength = sizeOfBox.getY();
    if ( maxLength < sizeOfBox.getZ() )
        maxLength = sizeOfBox.getZ();

    Sphere bigBrother[2] = { Sphere(centerOfBox.getX(), centerOfBox.getY(), -maxLength*20., maxLength*10.),
                             Sphere(centerOfBox.getX(), centerOfBox.getY(), maxLength*20., maxLength*10.) };
        
    m_minPointOfBoundingBoxForBalls.setPoint(centerOfBox.getX()-maxLength*10., centerOfBox.getY()-maxLength*10., -maxLength*20.);
    m_maxPointOfBoundingBoxForBalls.setPoint(centerOfBox.getX()+maxLength*10., centerOfBox.getY()+maxLength*10., maxLength*20.);
    BallGenerator* ptrBall[2] = {rg_NULL, rg_NULL};
    VDCell*        ptrCell[2] = {rg_NULL, rg_NULL};
    for (rg_INT i=0; i<2; i++)  {
        ptrBall[i] = m_GlobalGeneratorList.addHead( BallGenerator(bigBrother[i].getCenter(), bigBrother[i].getRadius()) );
        ptrCell[i] = m_GlobalCellList.addHead( VDCell( -2-i ) );
        connectCellAndGenerator( ptrCell[i], ptrBall[i] );
    }

    VDFace* ptrFace = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), ptrCell[1], ptrCell[0] ) );
    VDLoop* ptrLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrFace, rg_TRUE ) );

    //  define connection between VDCell and VDFace.
    ptrCell[0]->addBoundingFace( ptrFace );
    ptrCell[1]->addBoundingFace( ptrFace );

}



void rg_SphereSetVoronoiDiagram::constructSmallWorlds()
{
    rg_dList<BallGenerator*> isolatedBalls;
    collectIsolatedBalls( isolatedBalls );

    rg_dList<BallGenerator*> isolatedBallsInDisconnectedBigWorlds;

    BallGenerator* currBall = rg_NULL;
    while ( isolatedBalls.getSize() > 0 )  { 
        currBall = isolatedBalls.popFront();

        if ( currBall->getCell()->getNumOfBoundingFaces() != 0 )
            continue;

        //  1. find big brothers.        
        BallGenerator* bigBrother[2]     = { rg_NULL, rg_NULL };
        VDFace*        faceByBigBrothers = findBigBrothers( currBall, bigBrother );

        if ( faceByBigBrothers == rg_NULL )  {
            isolatedBallsInDisconnectedBigWorlds.add( currBall );
            continue;
        }

        //  2. find initial Voronoi vertex in the small world.
        rg_INT           initialSizeOfVIDIC = 11;
        VIDIC            theVIDIC( initialSizeOfVIDIC );
        rg_dList< Gate > theGateStack;
        VDVertex*        pVertex = rg_NULL;
        rg_FLAG          bValidSmallWorld = findInitialVertexInSmallWorld(currBall, bigBrother, pVertex, 
                                                                          theVIDIC, theGateStack );

        if ( bValidSmallWorld )  {
            if ( pVertex != rg_NULL )  
                traceEdgesInSmallWorld( faceByBigBrothers, theVIDIC, theGateStack );
            else  
                constructSmallWorldWithoutVertices(currBall, bigBrother, faceByBigBrothers );                
                //setTopologyBySingleSmallBall(currBall, bigBrother, faceByBigBrothers );
        }
        else  {
                isolatedBalls.addTail( currBall );
        }
    }

    //  we have to deal with the disconnected big worlds.
}



void rg_SphereSetVoronoiDiagram::collectIsolatedBalls(rg_dList<BallGenerator*>& isolatedBalls)
{
    rg_dList<BallGenerator*> isoBalls;

    BallGenerator* currBall = rg_NULL;
    m_GlobalGeneratorList.reset4Loop();
    while ( m_GlobalGeneratorList.setNext4Loop() )  {
        currBall = m_GlobalGeneratorList.getpEntity();

        if ( currBall->getCell()->getNumOfBoundingFaces() == 0 )
            isoBalls.add( currBall );
    }

    rg_INT numBalls       = isoBalls.getSize();
    BallGenerator** balls = isoBalls.getArray();

    qsort( (void *) balls, numBalls, sizeof(BallGenerator*), compareBallsByRadiusInDescendingOrder);

    for (rg_INT i=0; i<numBalls; i++)
        isolatedBalls.add( balls[i] );
}



VDFace* rg_SphereSetVoronoiDiagram::findBigBrothers( BallGenerator* smallBrother, BallGenerator** bigBrother)
{
    rg_REAL        minDistance = DBL_MAX;
    BallGenerator* closestBall = rg_NULL;
    Sphere         smallBall   = smallBrother->getBall();

    BallGenerator* currBall = rg_NULL;
    m_GlobalGeneratorList.reset4Loop();
    while ( m_GlobalGeneratorList.setNext4Loop() )  {
        currBall = m_GlobalGeneratorList.getpEntity();

        if ( currBall->getCell()->getNumOfBoundingFaces() == 0 )
            continue;

        rg_REAL distance = smallBall.distanceToSphere( currBall->getBall() );
        if ( distance < minDistance )  {
            minDistance = distance;
            closestBall = currBall;
        }
    }

    bigBrother[0] = closestBall;
    
    Sphere  bigBall = bigBrother[0]->getBall();
    VDCell* bigCell = bigBrother[0]->getCell();

    rg_dList<VDFace*>* boundingFacesOfBigCell = bigCell->getBoungindFaces();
    BallGenerator*     candidateBall          = rg_NULL;

    minDistance = DBL_MAX;
    closestBall = rg_NULL;
    VDFace* faceByBigBrothers = rg_NULL;
    VDFace* currFace          = rg_NULL;
    boundingFacesOfBigCell->reset4Loop();
    while ( boundingFacesOfBigCell->setNext4Loop() )  {
        currFace = boundingFacesOfBigCell->getEntity();
    
        if ( currFace->isOnInfinity() )
            continue;

        if ( bigCell == currFace->getRightCell() )
            candidateBall = (BallGenerator*) currFace->getLeftCell()->getGenerator();
        else
            candidateBall = (BallGenerator*) currFace->getRightCell()->getGenerator();

        Sphere enclosingSphere( bigBall );
        enclosingSphere.makeContainingSphereWith(candidateBall->getBall());
        if ( smallBall.isContainedIn(enclosingSphere) == rg_FALSE )
            continue;

        Sphere tangentSphere[2];
        rg_INT numTS = makeCircumcircleOnCenterPlane( smallBall, bigBall, 
                                                      candidateBall->getBall(), 
                                                      tangentSphere[0], tangentSphere[1] );

        if ( numTS == 2 )  {
            rg_REAL distance = smallBall.distanceToSphere( candidateBall->getBall() );
            if ( distance < minDistance )  {
                minDistance       = distance;
                closestBall       = candidateBall;
                faceByBigBrothers = currFace;
            }
        }
    }    
    
    if ( faceByBigBrothers != rg_NULL )  {

        bigBrother[1] = closestBall;

        //  bigBrother[0] must be left ball of the face by two big brothers.
        if (bigCell == faceByBigBrothers->getRightCell() )  {
            candidateBall = bigBrother[0];
            bigBrother[0] = bigBrother[1];
            bigBrother[1] = candidateBall;
        }
    }

    return faceByBigBrothers;
}




rg_FLAG rg_SphereSetVoronoiDiagram::findInitialVertexInSmallWorld(BallGenerator*  smallBrother, 
                                                                  BallGenerator** bigBrother,
                                                                  VDVertex*&      ptrVertex,
                                                                  VIDIC&          theVIDIC, 
                                                                  rg_dList<Gate>& theGateStack )
{
    BallGenerator* gateBall[3] = { bigBrother[0], smallBrother, bigBrother[1] };   
    Sphere         ball[3]     = { gateBall[0]->getBall(), gateBall[1]->getBall(), gateBall[2]->getBall() };

    Sphere enclosingSphereOfBigBrothers( bigBrother[0]->getBall() );
    enclosingSphereOfBigBrothers.makeContainingSphereWith( bigBrother[1]->getBall() );

    rg_INT i=0;
	for ( i=0; i<3; i++)  {
        gateBall[i]->m_checkForBucket = rg_TRUE;
    }

    rg_dList<BallGenerator*> candidateBalls;
    m_bucket.getBallsIntersectingSphere( enclosingSphereOfBigBrothers, candidateBalls );

    Sphere tangentSphere[2];
    rg_INT numTS = 0;

    Sphere emptyTangentSphere[2];
    rg_INT numEmptyTS = 0;

    rg_FLAG bIntersected    = rg_FALSE;
    BallGenerator* endBall  = rg_NULL;
    BallGenerator* currBall = rg_NULL;
    candidateBalls.reset4Loop();
	while( candidateBalls.setNext4Loop() )  {
	    currBall = candidateBalls.getEntity();

        if ( currBall->getCell()->getNumOfBoundingFaces() != 0 )
            continue;


        numTS = computeSphereTangentTo4SpheresOutside(
                    ball[0], ball[1], ball[2], currBall->getBall(), tangentSphere);

        if ( numTS > 0 )  {
            currBall->m_checkForBucket = rg_TRUE;
            for ( i=0; i<numTS; i++)  {
                rg_dList<BallGenerator*> intersectingBalls;
                m_bucket.getBallsIntersectingSphere( tangentSphere[i], intersectingBalls );

                if ( intersectingBalls.getSize() == 0 )
                    emptyTangentSphere[numEmptyTS++] = tangentSphere[i];
                else
                    bIntersected = rg_TRUE;
            }            
            currBall->m_checkForBucket = rg_FALSE;

            if ( numEmptyTS > 0 )  {            
                endBall = currBall;
                break;
            }
        }
    }

    for ( i=0; i<3; i++)  {
        gateBall[i]->m_checkForBucket = rg_FALSE;
    }

    
    if ( endBall != rg_NULL )  {
        //  1. make center plane.
        rg_Vector3D vecStoB0 = ball[0].getCenter() - emptyTangentSphere[0].getCenter();
        rg_Vector3D vecStoB1 = ball[1].getCenter() - emptyTangentSphere[0].getCenter();

        rg_Vector3D normalCenterPlane = vecStoB0.crossProduct( vecStoB1 );
        normalCenterPlane.normalize();
        rg_REAL     dOfCenterPlane = -normalCenterPlane.innerProduct(ball[0].getCenter());

        //  2. construct tangent point on end ball.
        rg_Vector3D vecStoEndBall = endBall->getCenter() - emptyTangentSphere[0].getCenter();
        vecStoEndBall.normalize();
        rg_Point3D  tangentPtOnEndBall = emptyTangentSphere[0].getCenter() + (emptyTangentSphere[0].getRadius()*vecStoEndBall);

        rg_REAL distance = normalCenterPlane.innerProduct(tangentPtOnEndBall) + dOfCenterPlane;


        //  3. make gates.
        Gate   initialGates[NUM_DEFINING_CELLS_OF_ONE_VERTEX];

        if ( distance > 0 )  {
            initialGates[0].setGate( gateBall[2], gateBall[1], gateBall[0], endBall, rg_NULL );
            initialGates[1].setGate( endBall, gateBall[0], gateBall[1], gateBall[2], rg_NULL );
            initialGates[2].setGate( endBall, gateBall[1], gateBall[2], gateBall[0], rg_NULL );
            initialGates[3].setGate( endBall, gateBall[2], gateBall[0], gateBall[1], rg_NULL );
        }
        else {
            initialGates[0].setGate( gateBall[0], gateBall[1], gateBall[2], endBall, rg_NULL );
            initialGates[1].setGate( endBall, gateBall[1], gateBall[0], gateBall[2], rg_NULL );
            initialGates[2].setGate( endBall, gateBall[2], gateBall[1], gateBall[0], rg_NULL );
            initialGates[3].setGate( endBall, gateBall[0], gateBall[2], gateBall[1], rg_NULL );
        }

        ptrVertex = makeVoronoiVertex( rg_FALSE, emptyTangentSphere[0] );


        //  make and insert 4 gates of initial voronoi vertex into Gate-Stack.
        Gate* ptrGate = rg_NULL;        
        for ( i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        {
            ptrGate = theGateStack.addTail( initialGates[i] );

            ptrGate->setStartVertex( ptrVertex );
            ptrGate->setNodePosInStack( theGateStack.getTail() );

            ptrVertex->setGate( i, ptrGate );
        }   

        //  insert vertex configuration of initial voronoi vertex into VIDIC.
        theVIDIC.insertVIDICItem( VIDICItem( initialGates[0].getGateBallCCW1(),
                                             initialGates[0].getGateBallCCW2(),
                                             initialGates[0].getGateBallCCW3(),
                                             initialGates[0].getStartBall(),
                                             ptrVertex) 
                                 );
        return rg_TRUE;
    }
    else  {
        if ( bIntersected )
            return rg_FALSE;
        else
            return rg_TRUE;
    }
}




void rg_SphereSetVoronoiDiagram::constructSmallWorldWithoutVertices(BallGenerator*  smallBrother, 
                                                                    BallGenerator** bigBrother,
                                                                    VDFace*         faceByBigBrothers)
{
    rg_dList< VDFace* > smallerWorlds;
    VDFace* firstWorld = faceByBigBrothers;

    while ( smallBrother->getCell()->getNumOfBoundingFaces() == 0 )  {    
        VDFace* currWorld = rg_NULL;
        smallerWorlds.add(firstWorld);
        while ( smallerWorlds.getSize() > 0 )  {
            currWorld = smallerWorlds.popFront();

            BallGenerator* currBrother[2] = { (BallGenerator*) currWorld->getLeftCell()->getGenerator(), 
                                              (BallGenerator*) currWorld->getRightCell()->getGenerator() };

            if ( bigBrother[0] != currBrother[0] && bigBrother[1] != currBrother[1] )
                continue;

            VDFace* facesOfSmallWorld[2] = {rg_NULL, rg_NULL};
            setTopologyByTwoBigBallsInSmallWorldWithoutVertices(currBrother, currWorld, facesOfSmallWorld);

            if ( facesOfSmallWorld[0] != rg_NULL )  {       
                smallerWorlds.add( facesOfSmallWorld[0] );
                smallerWorlds.add( facesOfSmallWorld[1] );
            }
        }
    }
}




void rg_SphereSetVoronoiDiagram::setTopologyByTwoBigBallsInSmallWorldWithoutVertices(
                                                              BallGenerator** bigBrother,
                                                              VDFace*         faceByBigBrothers,
                                                              VDFace**        facesBySmallAndTwoBigBrothers)
{
    Sphere bigBall[2]    = { bigBrother[0]->getBall(),    bigBrother[1]->getBall() };   
    Sphere enclosingSphere(bigBall[0]);
    enclosingSphere.makeContainingSphereWith(bigBall[1]);

    bigBrother[0]->m_checkForBucket = rg_TRUE;
    bigBrother[1]->m_checkForBucket = rg_TRUE;

    rg_dList<BallGenerator*> ballsBetweenTwoBrothers;
    m_bucket.getBallsIntersectingSphere( enclosingSphere, ballsBetweenTwoBrothers );

    bigBrother[0]->m_checkForBucket = rg_FALSE;
    bigBrother[1]->m_checkForBucket = rg_FALSE;

    if ( ballsBetweenTwoBrothers.getSize() == 0 )
        return;


    BallGenerator* smallBrother        = rg_NULL;
    BallGenerator* largestSmallBrother = rg_NULL;
    Sphere largestTangentSphere[2];
    rg_INT numTS = 0;
    while ( ballsBetweenTwoBrothers.getSize() > 0 )  {
        smallBrother = ballsBetweenTwoBrothers.popFront();

        if ( smallBrother->getCell()->getNumOfBoundingFaces() != 0 )
            continue;

        numTS = makeCircumcircleOnCenterPlane( bigBall[0], bigBall[1], smallBrother->getBall(), 
                                               largestTangentSphere[0], largestTangentSphere[1] );
        if ( numTS == 2 )  {
            largestSmallBrother = smallBrother;
            break;
        }
    }

    while ( ballsBetweenTwoBrothers.getSize() > 0 )  {
        smallBrother = ballsBetweenTwoBrothers.popFront();
        Sphere smallBall(smallBrother->getBall());

        if ( smallBrother->getCell()->getNumOfBoundingFaces() != 0 )
            continue;
        
        if ( smallBall.isThereIntersectionWith( largestTangentSphere[0] ) )  {   
            numTS = makeCircumcircleOnCenterPlane( bigBall[0], bigBall[1], smallBall, 
                                           largestTangentSphere[0], largestTangentSphere[1] );

            if ( numTS == 2 )  
                largestSmallBrother = smallBrother;
        }
    }

    if ( largestSmallBrother != rg_NULL )  {   
        setTopologyBySingleSmallBall( largestSmallBrother, bigBrother, faceByBigBrothers );

        rg_dNode<VDFace>* lastNode = m_GlobalFaceList.getTail();

        facesBySmallAndTwoBigBrothers[0] = lastNode->getPrev()->getpEntity();
        facesBySmallAndTwoBigBrothers[1] = lastNode->getpEntity();
    }
}



void rg_SphereSetVoronoiDiagram::setTopologyBySingleSmallBall(BallGenerator*  smallBrother, 
                                                              BallGenerator** bigBrother,
                                                              VDFace*         faceByBigBrothers)
{
    ///////////////////////////////////////////////////////////////
    //
    //  1. make an edge without start and end vertices.
    VDEdge* ptrEdge = m_GlobalEdgeList.addTail( VDEdge( m_GlobalEdgeList.getSize(), rg_NULL, rg_NULL ) );
    ptrEdge->setEdgeType( EDGE_WITHOUT_VERTICES );

    ///////////////////////////////////////////////////////////////
    //
    //  2. make three partial edges.
    VDPartialEdge* ptrPartialEdges[3] = {rg_NULL, rg_NULL, rg_NULL};

    //     2.1 insert 3 partial edges of the edge into GPL(Global Partial edge List).
    rg_INT i=0;
	for ( i=0; i<3; i++ )  {
        ptrPartialEdges[i] = m_GlobalPartialedgeList.addTail( VDPartialEdge( m_GlobalPartialedgeList.getSize(), ptrEdge ) );

        ptrPartialEdges[i]->setNextPartialEdgeInLoop( ptrPartialEdges[i] );
        ptrPartialEdges[i]->setPreviousPartialEdgeInLoop( ptrPartialEdges[i] );
    }

    //     2.2 connect one partial edge to the edge.
    ptrEdge->setPartialEdge( ptrPartialEdges[0] );

    //     2.3 connect partial edges in Radial Cycle.
    ptrPartialEdges[0]->setNextPartialEdgeInRadialCycle( ptrPartialEdges[1] ); // big[0] -> big[1] ; false
    ptrPartialEdges[1]->setNextPartialEdgeInRadialCycle( ptrPartialEdges[2] ); // big[0] -> small  ; true
    ptrPartialEdges[2]->setNextPartialEdgeInRadialCycle( ptrPartialEdges[0] ); // small  -> big[1] ; true

    ptrPartialEdges[0]->isRightOrientationInLoop(rg_FALSE);
    ptrPartialEdges[1]->isRightOrientationInLoop(rg_TRUE);
    ptrPartialEdges[2]->isRightOrientationInLoop(rg_TRUE);


    ///////////////////////////////////////////////////////////////
    //
    //  3. make inner-loop of faceByBigBrothers.
    VDLoop* ptrInnerLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), 
                                                             faceByBigBrothers, rg_FALSE ) );
    faceByBigBrothers->addLoop( ptrInnerLoop );
    ptrInnerLoop->addPartialEdgeToLoop( ptrPartialEdges[0] );


    ///////////////////////////////////////////////////////////////
    //
    //  4. make faces and its outer-loop between small ball and each big ball
    VDFace* ptrFaceFromBig0ToSmall = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), 
                                                                       smallBrother->getCell(), 
                                                                       bigBrother[0]->getCell() ) );
    VDLoop* ptrLoopOfBig0ToSmall = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), 
                                                                     ptrFaceFromBig0ToSmall, rg_TRUE ) );
    ptrFaceFromBig0ToSmall->addLoop( ptrLoopOfBig0ToSmall );
    ptrLoopOfBig0ToSmall->addPartialEdgeToLoop( ptrPartialEdges[1] );

    //  define connection between VDCell and VDFace.
    smallBrother->getCell()->addBoundingFace( ptrFaceFromBig0ToSmall );
    bigBrother[0]->getCell()->addBoundingFace( ptrFaceFromBig0ToSmall );

    
    VDFace* ptrFaceFromSmallToBig1 = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), 
                                                                       bigBrother[1]->getCell(),
                                                                       smallBrother->getCell() ) );
    VDLoop* ptrLoopOfSmallToBig1 = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), 
                                                                     ptrFaceFromSmallToBig1, rg_TRUE ) );
    ptrFaceFromSmallToBig1->addLoop( ptrLoopOfSmallToBig1 );
    ptrLoopOfSmallToBig1->addPartialEdgeToLoop( ptrPartialEdges[2] );  
    //  define connection between VDCell and VDFace.
    bigBrother[1]->getCell()->addBoundingFace( ptrFaceFromSmallToBig1 );
    smallBrother->getCell()->addBoundingFace( ptrFaceFromSmallToBig1 );
    
}



void rg_SphereSetVoronoiDiagram::traceEdgesInSmallWorld(VDFace*         faceByBigBrothers,
                                                        VIDIC&          theVIDIC,     
                                                        rg_dList<Gate>& theGateStack )
{
    rg_dList<VDLoop*> loopsInSmallWorld;

    while ( theGateStack.getSize() > 0 )  {

        Gate currGate = theGateStack.popFront();

        findEndVertexBasedOnAngularDistanceInAcceleration( currGate ); 

        VDVertex* ptrEndVertex = exploreEdgeToBeMadeByCurrentGate( currGate, theVIDIC, theGateStack );

        if ( m_status != rg_TRUE )
            break;
    
        setLocalTopologyByEdgeInSmallWorld( faceByBigBrothers,
                                            currGate, 
                                            ptrEndVertex,
                                            loopsInSmallWorld );
    }

    if ( m_status == rg_TRUE )  {   
        arrangePartialEdgesInLoops( loopsInSmallWorld );
    }
}



void rg_SphereSetVoronoiDiagram::setLocalTopologyByEdgeInSmallWorld( 
                                                                VDFace*         faceByBigBrothers,
                                                                const Gate&     gate, 
                                                                VDVertex*       endVertex,
                                                                rg_dList<VDLoop*>& loopsInSmallWorld)
{
    ///////////////////////////////////////////////////////////////
    //
    //  make Voronoi edge.
    //
    VDEdge* ptrEdge = m_GlobalEdgeList.addTail( VDEdge( m_GlobalEdgeList.getSize(), 
                                                        gate.getStartVertex(), 
                                                        endVertex ) );
    ptrEdge->setEdgeType( gate.getEdgeType() );

    //  connect edge to vertex as incident edge of vertex.
    gate.getStartVertex()->setIncidentEdge( ptrEdge );
    endVertex->setIncidentEdge( ptrEdge );


    ///////////////////////////////////////////////////////////////
    //
    //  make three partial edges
    //
    VDPartialEdge* partialEdges[3];

    rg_INT i=0;
	for ( i=0; i<3; i++ ) {
        partialEdges[i] = m_GlobalPartialedgeList.addTail( VDPartialEdge( m_GlobalPartialedgeList.getSize(), ptrEdge ) );

        partialEdges[i]->setNextPartialEdgeInLoop( partialEdges[i] );
        partialEdges[i]->setPreviousPartialEdgeInLoop( partialEdges[i] );
    }

    //  connect one partial edge to the edge.
    ptrEdge->setPartialEdge( partialEdges[0] );

    //  connect partial edges in Radial Cycle.
    partialEdges[0]->setNextPartialEdgeInRadialCycle( partialEdges[1] );
    partialEdges[1]->setNextPartialEdgeInRadialCycle( partialEdges[2] );
    partialEdges[2]->setNextPartialEdgeInRadialCycle( partialEdges[0] );



    ///////////////////////////////////////////////////////////////
    //
    //  make faces and its loops.
    //
    VDCell* gateCell[4] = { gate.getGateBallCCW1()->getCell(), gate.getGateBallCCW2()->getCell(), 
                            gate.getGateBallCCW3()->getCell(), gate.getGateBallCCW1()->getCell() };

    //  partial edges[0] <- gate ball 0 & gate ball 1
    //  partial edges[1] <- gate ball 1 & gate ball 2
    //  partial edges[2] <- gate ball 2 & gate ball 0

    for ( i=0; i<3; i++ ) {
    
        VDFace* ptrFace = gateCell[i]->findFaceToShareWith( gateCell[i+1] );
        VDLoop* ptrLoop = rg_NULL;

        if ( ptrFace == faceByBigBrothers )  {

            if ( loopsInSmallWorld.getFirstEntity()->isOuterLoop() == rg_TRUE )  {
                ptrLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrFace, rg_FALSE ) );
                loopsInSmallWorld.addHead( ptrLoop );

                ptrFace->addLoop( ptrLoop );

                if ( ptrFace->getLeftCell() == gateCell[i] && ptrFace->getRightCell() == gateCell[i+1] )
                    partialEdges[i]->isRightOrientationInLoop(rg_TRUE);
                else
                    partialEdges[i]->isRightOrientationInLoop(rg_FALSE);
            }
            else  {
                ptrLoop = loopsInSmallWorld.getFirstEntity();

                if ( ptrFace->getLeftCell() == gateCell[i] && ptrFace->getRightCell() == gateCell[i+1] )
                    partialEdges[i]->isRightOrientationInLoop(rg_TRUE);
                else
                    partialEdges[i]->isRightOrientationInLoop(rg_FALSE);
            }
        }
        else  {        
            if ( ptrFace == rg_NULL )
            {
                ptrFace = m_GlobalFaceList.addTail( VDFace( m_GlobalFaceList.getSize(), gateCell[i+1], gateCell[i] ) );
                ptrLoop = m_GlobalLoopList.addTail( VDLoop( m_GlobalLoopList.getSize(), ptrFace, rg_TRUE ) );
                loopsInSmallWorld.addTail( ptrLoop );

                //  define connection between VDCell and VDFace.
                gateCell[i]->addBoundingFace( ptrFace );
                gateCell[i+1]->addBoundingFace( ptrFace );

                //  define connection between VDFace and VDLoop.
                ptrFace->addLoop( ptrLoop );

                //  define connection between VDLoop and VDPartialEdge.
                //partialEdges[i]->setLoop(ptrLoop);
                partialEdges[i]->isRightOrientationInLoop(rg_TRUE);
            }
            else
            {
                ptrLoop = ptrFace->getLoop( 0 );

                //  define connection between VDLoop and VDPartialEdge.
                //partialEdges[i]->setLoop(ptrLoop);

                if (    ptrFace->getLeftCell() == gateCell[i] 
                     && ptrFace->getRightCell() == gateCell[i+1] )
                    partialEdges[i]->isRightOrientationInLoop(rg_TRUE);
                else
                    partialEdges[i]->isRightOrientationInLoop(rg_FALSE);
            }
        }
        ptrLoop->addPartialEdgeToLoop( partialEdges[i] );
    }
}



void  rg_SphereSetVoronoiDiagram::arrangePartialEdgesInLoops( rg_dList<VDLoop*>& loopsInSmallWorld )
{
    VDLoop* currLoop = rg_NULL;

    loopsInSmallWorld.reset4Loop();
    while( loopsInSmallWorld.setNext4Loop() ) {
        currLoop = loopsInSmallWorld.getEntity();

        rg_dList<VDPartialEdge*> multipleOuterLoop;

        if ( currLoop->constructLoopCycle( multipleOuterLoop ) == rg_FALSE )  {
            if ( multipleOuterLoop.getSize() > 0 )
                makeAdditiveOuterLoopOfFace( currLoop->getFace(), multipleOuterLoop );
        }
    }
}



////////////////////////////////////////////////////////////////////////////////////////////////
//
//  VORONOI VERTEX COMPUTATION

rg_INT rg_SphereSetVoronoiDiagram::evaluate_Voronoi_Vertex(Sphere*  s1, Sphere* s2, Sphere* s3, Sphere* s4,
                                                           Sphere*& tangentSpheres)
{	
    //donguk's version (2004. 7. 5)
	Sphere* arrSphere[4] = {s1, s2, s3, s4};
	
	if((s1->getRadius() == s2->getRadius()) && (s2->getRadius() == s3->getRadius()) && (s3->getRadius()==s4->getRadius()))
	{
		return evaluate_Voronoi_Vertex_SameRadius( s1, s2, s3, s4, tangentSpheres);
	}
	else
	{
        Sphere balls[4] = { *s1, *s2, *s3, *s4 };

        sort4SpheresByRadiusInDescendingPowers( balls );
		
        // initilizing variables for calculating Vertex
		// 3 spheres are shrinked and translated by minimum sphere
        rg_Point3D centerOfSphere[4];
        int i=0;
		for (i=0; i<4; i++)  {
            arrSphere[i]      = &balls[i];
            centerOfSphere[i] = arrSphere[i]->getCenter();
        }
        
        double radius[4] = { arrSphere[0]->getRadius(), 
                             arrSphere[1]->getRadius(), 
                             arrSphere[2]->getRadius(), 
                             arrSphere[3]->getRadius() };

			
		double x1 = centerOfSphere[0].getX() - centerOfSphere[3].getX();
		double y1 = centerOfSphere[0].getY() - centerOfSphere[3].getY();
		double z1 = centerOfSphere[0].getZ() - centerOfSphere[3].getZ();
		double r1 = radius[0] - radius[3];
		//double r1 = arrSphere[0]->getRadius() - arrSphere[3]->getRadius();

		double x2 = centerOfSphere[1].getX() - centerOfSphere[3].getX();
		double y2 = centerOfSphere[1].getY() - centerOfSphere[3].getY();
		double z2 = centerOfSphere[1].getZ() - centerOfSphere[3].getZ();
		double r2 = radius[1] - radius[3];
		//double r2 = arrSphere[1]->getRadius() - arrSphere[3]->getRadius();

		double x3 = centerOfSphere[2].getX() - centerOfSphere[3].getX();
		double y3 = centerOfSphere[2].getY() - centerOfSphere[3].getY();
		double z3 = centerOfSphere[2].getZ() - centerOfSphere[3].getZ();
		double r3 = radius[2] - radius[3];

//		double determinant1 = y3*x2*z1+x3*y1*z2+x1*y2*z3-y2*x3*z1-x2*y1*z3-x1*y3*z2;
		double determinant1 = getDeterminant(x1,y1,z1,x2,y2,z2,x3,y3,z3);
		double determinant2 = getDeterminant(x1,y1,r1,x2,y2,r2,x3,y3,r3);
		double determinant3 = getDeterminant(x1,r1,z1,x2,r2,z2,x3,r3,z3);
		double determinant4 = getDeterminant(r1,y1,z1,r2,y2,z2,r3,y3,z3);

		//system of equations
		//eq1: x1*x + y1*y + z1*z = rhs1
		//eq2: x2*x + y2*y + z2*z = rhs2
		//eq3: x3*x + y3*y + z3*z = rhs3
		//eq4: x*x + y*y + z*z = r*r

		double rhs1 = (x1*x1+y1*y1+z1*z1-r1*r1)/2.0;
		double rhs2 = (x2*x2+y2*y2+z2*z2-r2*r2)/2.0;
		double rhs3 = (x3*x3+y3*y3+z3*z3-r3*r3)/2.0;

		double a_x=0., b_x=0., c_x=0.;
		double a_y=0., b_y=0., c_y=0.;
		double a_z=0., b_z=0., c_z=0.;
		double a_r=0., b_r=0., c_r=0.;

		double x = 0., y = 0., z = 0.;

		//finding non-null minor
		if( !rg_ZERO(determinant1) ) 
		{
			//a_r * r^2 + b_r * r + c_r = 0

			//x = a_x * r + b_x;
			//y = a_y * r + b_y;
			//z = a_z * r + b_z;

			a_x = (-z1*r2*y3+z2*r1*y3+z1*y2*r3+r2*z3*y1-z2*y1*r3-y2*z3*r1)/determinant1;
			b_x = (-z1*y2*rhs3+z1*rhs2*y3+z2*y1*rhs3+y2*z3*rhs1-z2*rhs1*y3-rhs2*z3*y1)/determinant1;
			a_y = (x3*z1*r2-z2*x3*r1-z1*x2*r3-x1*z3*r2+z2*x1*r3+z3*x2*r1)/determinant1;
			b_y = (z2*x3*rhs1-z3*x2*rhs1-z2*x1*rhs3-x3*z1*rhs2+x1*z3*rhs2+z1*x2*rhs3)/determinant1;
			a_z = (x2*y1*r3+r2*x1*y3-r2*x3*y1+y2*x3*r1-x2*r1*y3-y2*x1*r3)/determinant1;
			b_z = (-x2*y1*rhs3+x2*rhs1*y3-y2*x3*rhs1+y2*x1*rhs3-rhs2*x1*y3+rhs2*x3*y1)/determinant1;

			// r = (-b_r + sqrt(b_r^2 - 4. * a_r * c_r))/(2.* a_r);
			// r = (-b_r - sqrt(b_r^2 - 4. * a_r * c_r))/(2.* a_r);
			a_r = a_y*a_y+a_x*a_x+a_z*a_z-1.0;
			b_r = 2.0*b_x*a_x+2.0*b_y*a_y+2.0*b_z*a_z;
			c_r = b_x*b_x+b_z*b_z+b_y*b_y;

			double det = b_r*b_r - 4.*a_r*c_r; //b^2 - 4ac
			if(det < 0)
			{
				tangentSpheres = rg_NULL;
				return 0;
			}

			double radius1 = (-b_r + sqrt(det))/(2.* a_r);
			double radius2 = (-b_r - sqrt(det))/(2.* a_r);

			if( radius1 > 0 && radius2 > 0 )
			{
				tangentSpheres = new Sphere[2];

				x = a_x * radius1 + b_x;
				y = a_y * radius1 + b_y;
				z = a_z * radius1 + b_z;
				tangentSpheres[0].setSphere(rg_Point3D(x,y,z)+centerOfSphere[3], radius1-radius[3]);

				x = a_x * radius2 + b_x;
				y = a_y * radius2 + b_y;
				z = a_z * radius2 + b_z;
				tangentSpheres[1].setSphere(rg_Point3D(x,y,z)+centerOfSphere[3], radius2-radius[3]);

				return 2;
			}
			else if( radius1 > 0 )
			{
				tangentSpheres = new Sphere[1];

				x = a_x * radius1 + b_x;
				y = a_y * radius1 + b_y;
				z = a_z * radius1 + b_z;
				tangentSpheres[0].setSphere(rg_Point3D(x,y,z)+centerOfSphere[3], radius1-radius[3]);

				return 1;
			}
			else if( radius2 > 0 )
			{
				tangentSpheres = new Sphere[1];

				x = a_x * radius2 + b_x;
				y = a_y * radius2 + b_y;
				z = a_z * radius2 + b_z;
				tangentSpheres[0].setSphere(rg_Point3D(x,y,z)+centerOfSphere[3], radius2-radius[3]);

				return 1;
			}
			else
			{
				tangentSpheres = rg_NULL;
				return 0;
			}

		}
		else if( !rg_ZERO(determinant2) ) 
		{
			//a_z * z^2 + b_z * z + c_z = 0

			//x = a_x * z + b_x;
			//y = a_y * z + b_y;
			//r = a_r * z + b_r;

			a_x = (z1*r2*y3-z2*r1*y3-z1*y2*r3-r2*z3*y1+z2*y1*r3+y2*z3*r1)/determinant2;
			b_x = (r2*y1*rhs3-y2*r1*rhs3+rhs2*r1*y3-r2*rhs1*y3+y2*r3*rhs1-rhs2*y1*r3)/determinant2;
			a_y = (-x3*z1*r2+z2*x3*r1+z1*x2*r3+x1*z3*r2-z2*x1*r3-z3*x2*r1)/determinant2;
			b_y = (-x1*r2*rhs3+x3*r2*rhs1+x2*r1*rhs3+r3*x1*rhs2-x2*r3*rhs1-r1*x3*rhs2)/determinant2;
			a_r = (z3*x2*y1-z1*x2*y3+z2*x1*y3-z2*x3*y1+x3*z1*y2-x1*z3*y2)/determinant2;
			b_r = (-x2*y1*rhs3+x2*rhs1*y3-y2*x3*rhs1+y2*x1*rhs3-rhs2*x1*y3+rhs2*x3*y1)/determinant2;

			a_z = a_x*a_x+a_y*a_y+1.0-a_r*a_r;
			b_z = 2.0*b_x*a_x+2.0*b_y*a_y-2.0*b_r*a_r;
			c_z = b_x*b_x+b_y*b_y-b_r*b_r;

			double det = b_z*b_z - 4.0*a_z*c_z;
			if(det < 0)
			{
				tangentSpheres = rg_NULL;
				return 0;
			}

			double zValue1 = (-b_z + sqrt(det))/(2.0*a_z);
			double zValue2 = (-b_z - sqrt(det))/(2.0*a_z);

			double rad1 = a_r*zValue1 + b_r;
			double rad2 = a_r*zValue2 + b_r;
			if( rad1 > 0 && rad2 > 0) 
			{
				tangentSpheres = new Sphere [2];
				x = a_x*zValue1 + b_x;
				y = a_y*zValue1 + b_y;
				tangentSpheres[0].setSphere(rg_Point3D(x,y,zValue1)+centerOfSphere[3], rad1-radius[3]);

				x = a_x*zValue2 + b_x;
				y = a_y*zValue2 + b_y;
				tangentSpheres[1].setSphere(rg_Point3D(x,y,zValue2)+centerOfSphere[3], rad2-radius[3]);

				return 2;
			}
			else if( rad1 > 0 )
			{
				tangentSpheres = new Sphere [1];
				x = a_x*zValue1 + b_x;
				y = a_y*zValue1 + b_y;
				tangentSpheres[0].setSphere(rg_Point3D(x,y,zValue1)+centerOfSphere[3], rad1-radius[3]);

				return 1;
			}
			else if( rad2 > 0 )
			{
				tangentSpheres = new Sphere [1];
				x = a_x*zValue2 + b_x;
				y = a_y*zValue2 + b_y;
				tangentSpheres[0].setSphere(rg_Point3D(x,y,zValue2)+centerOfSphere[3], rad2-radius[3]);
				return 1;
			}
			else
			{
				tangentSpheres = rg_NULL;
				return 0;
			}
		
		
		}
		else if( !rg_ZERO(determinant3) ) 
		{
			//a_y * y^2 + b_y * y + c_y = 0

			//x = a_x * y + b_x;
			//r = a_r * y + b_r;
			//z = a_z * y + b_z;
			a_x = (z1*r2*y3-z2*r1*y3-z1*y2*r3-r2*z3*y1+z2*y1*r3+y2*z3*r1)/determinant3;
			b_x = (-z1*r2*rhs3+z1*r3*rhs2+r2*z3*rhs1+z2*r1*rhs3-rhs2*z3*r1-z2*rhs1*r3)/determinant3;
			a_r = (z3*x2*y1-z1*x2*y3+z2*x1*y3-z2*x3*y1+x3*z1*y2-x1*z3*y2)/determinant3;
			b_r = (z2*x3*rhs1-z3*x2*rhs1-z2*x1*rhs3-x3*z1*rhs2+x1*z3*rhs2+z1*x2*rhs3)/determinant3;
			a_z = (-x2*y1*r3+x2*r1*y3-y2*x3*r1+y2*x1*r3-r2*x1*y3+r2*x3*y1)/determinant3;
			b_z = (x1*r2*rhs3-r3*x1*rhs2-x2*r1*rhs3+r1*x3*rhs2-x3*r2*rhs1+x2*r3*rhs1)/determinant3;

			a_y = a_x*a_x+a_z*a_z+1.0-a_r*a_r;
			b_y = 2.0*b_x*a_x+2.0*b_z*a_z-2.0*b_r*a_r;
			c_y = b_x*b_x-b_r*b_r+b_z*b_z;

			double det = b_y*b_y - 4.0*a_y*c_y;
			if(det < 0)
			{
				tangentSpheres = rg_NULL;
				return 0;
			}

			double yValue1 = (-b_y + sqrt(det))/(2.0*a_y);
			double yValue2 = (-b_y - sqrt(det))/(2.0*a_y);

			double rad1 = a_r*yValue1 + b_r;
			double rad2 = a_r*yValue2 + b_r;
			if( rad1 > 0 && rad2 > 0) 
			{
				tangentSpheres = new Sphere [2];
				x = a_x*yValue1 + b_x;
				z = a_z*yValue1 + b_z;
				tangentSpheres[0].setSphere(rg_Point3D(x,yValue1,z)+centerOfSphere[3], rad1-radius[3]);

				x = a_x*yValue2 + b_x;
				z = a_z*yValue2 + b_z;
				tangentSpheres[1].setSphere(rg_Point3D(x,yValue2,z)+centerOfSphere[3], rad2-radius[3]);

				return 2;
			}
			else if( rad1 > 0 )
			{
				tangentSpheres = new Sphere [1];
				x = a_x*yValue1 + b_x;
				z = a_z*yValue1 + b_z;
				tangentSpheres[0].setSphere(rg_Point3D(x,yValue1,z)+centerOfSphere[3], rad1-radius[3]);

				return 1;
			}
			else if( rad2 > 0 )
			{
				tangentSpheres = new Sphere [1];
				x = a_x*yValue2 + b_x;
				z = a_z*yValue2 + b_z;
				tangentSpheres[0].setSphere(rg_Point3D(x,yValue2,z)+centerOfSphere[3], rad2-radius[3]);

				return 1;
			}
			else
			{
				tangentSpheres = rg_NULL;
				return 0;
			}
		
		}
		else if( !rg_ZERO(determinant4) ) 
		{
			//a_x * x^2 + b_x * x + c_x = 0

			//r = a_r * x + b_r;
			//y = a_y * x + b_y;
			//z = a_z * x + b_z;
			a_y = (-x3*z1*r2+z2*x3*r1+z1*x2*r3+x1*z3*r2-z2*x1*r3-z3*x2*r1)/determinant4;
			b_y = (z1*r2*rhs3-z1*r3*rhs2-z2*r1*rhs3+z2*rhs1*r3-r2*z3*rhs1+rhs2*z3*r1)/determinant4;
			a_r = (z3*x2*y1-z1*x2*y3+z2*x1*y3-z2*x3*y1+x3*z1*y2-x1*z3*y2)/determinant4;
			b_r = (-z1*y2*rhs3+z1*rhs2*y3+z2*y1*rhs3+y2*z3*rhs1-z2*rhs1*y3-rhs2*z3*y1)/determinant4;
			a_z = (-x2*y1*r3+x2*r1*y3-y2*x3*r1+y2*x1*r3-r2*x1*y3+r2*x3*y1)/determinant4;
			b_z = (y2*r1*rhs3-rhs2*r1*y3-r2*y1*rhs3+rhs2*y1*r3+r2*rhs1*y3-y2*r3*rhs1)/determinant4;

			a_x = a_z*a_z+a_y*a_y+1.0-a_r*a_r;
			b_x = 2.0*b_y*a_y-2.0*b_r*a_r+2.0*b_z*a_z;
			c_x = b_z*b_z+b_y*b_y-b_r*b_r;

			double det = b_x*b_x - 4.0*a_x*c_x;
			if(det < 0)
			{
				tangentSpheres = rg_NULL;
				return 0;
			}

			double xValue1 = (-b_x + sqrt(det))/(2.0*a_x);
			double xValue2 = (-b_x - sqrt(det))/(2.0*a_x);

			double rad1 = a_r*xValue1 + b_r;
			double rad2 = a_r*xValue2 + b_r;
			if( rad1 > 0 && rad2 > 0) 
			{
				tangentSpheres = new Sphere [2];
				z = a_z*xValue1 + b_z;
				y = a_y*xValue1 + b_y;
				tangentSpheres[0].setSphere(rg_Point3D(xValue1,y,z)+centerOfSphere[3], rad1-radius[3]);

				z = a_z*xValue2 + b_z;
				y = a_y*xValue2 + b_y;
				tangentSpheres[1].setSphere(rg_Point3D(xValue2,y,z)+centerOfSphere[3], rad2-radius[3]);

				return 2;
			}
			else if( rad1 > 0 )
			{
				tangentSpheres = new Sphere [1];
				z = a_z*xValue1 + b_z;
				y = a_y*xValue1 + b_y;
				tangentSpheres[0].setSphere(rg_Point3D(xValue1,y,z)+centerOfSphere[3], rad1-radius[3]);

				return 1;
			}
			else if( rad2 > 0 )
			{
				tangentSpheres = new Sphere [1];
				z = a_z*xValue2 + b_z;
				y = a_y*xValue2 + b_y;
				tangentSpheres[0].setSphere(rg_Point3D(xValue2,y,z)+centerOfSphere[3], rad2-radius[3]);
				return 1;
			}
			else
			{
				tangentSpheres = rg_NULL;
				return 0;
			}
		
		}
		else
		{
			tangentSpheres = rg_NULL;
			return 0;
		}
	}
	
	tangentSpheres = rg_NULL;
	return 0;
}




//  Modified by Youngsong Cho 2004-04-27
rg_INT rg_SphereSetVoronoiDiagram::evaluate_Voronoi_Vertex_SameRadius(Sphere* s1, Sphere* s2, Sphere* s3, Sphere* s4, 
                                                                      Sphere*& tangentSpheres)
{
	Sphere* arrSphere[4] = {s1, s2, s3, s4};
	
	// initilizing variables for calculating Vertex
	// 3 spheres are shrinked and translated by minimum sphere
    rg_Point3D centerOfSphere[4];
    int i=0;
	for (i=0; i<4; i++)  {
        centerOfSphere[i] = arrSphere[i]->getCenter();
    }

	double x1 = centerOfSphere[0].getX() - centerOfSphere[3].getX();
	double y1 = centerOfSphere[0].getY() - centerOfSphere[3].getY();
	double z1 = centerOfSphere[0].getZ() - centerOfSphere[3].getZ();
	double r1 = arrSphere[0]->getRadius() - arrSphere[3]->getRadius();

	double x2 = centerOfSphere[1].getX() - centerOfSphere[3].getX();
	double y2 = centerOfSphere[1].getY() - centerOfSphere[3].getY();
	double z2 = centerOfSphere[1].getZ() - centerOfSphere[3].getZ();
	double r2 = arrSphere[1]->getRadius() - arrSphere[3]->getRadius();

	double x3 = centerOfSphere[2].getX() - centerOfSphere[3].getX();
	double y3 = centerOfSphere[2].getY() - centerOfSphere[3].getY();
	double z3 = centerOfSphere[2].getZ() - centerOfSphere[3].getZ();
	double r3 = arrSphere[2]->getRadius() - arrSphere[3]->getRadius();

	double det_A =   ( x1 * ( ( y2 * z3 ) - ( z2 * y3 ) ) )
	 			   - ( y1 * ( ( x2 * z3 ) - ( z2 * x3 ) ) )
				   + ( z1 * ( ( x2 * y3 ) - ( y2 * x3 ) ) );
    /*
	double det_A = getDeterminent(x1, y1, z1, 
                                  x2, y2, z2, 
                                  x3, y3, z3);
    */
		
	if ( det_A == 0 )
	{
		tangentSpheres = rg_NULL;
		return 0;
	}

	// evalueate Final D-Sphere's X, Y, Z coordinates and radius
	double l1 = (x1*x1) + (y1*y1) + (z1*z1);
	double l2 = (x2*x2) + (y2*y2) + (z2*z2);
	double l3 = (x3*x3) + (y3*y3) + (z3*z3);
    /*
	double l1 = pow(x1,2) + pow(y1,2) + pow(z1,2);
	double l2 = pow(x2,2) + pow(y2,2) + pow(z2,2);
	double l3 = pow(x3,2) + pow(y3,2) + pow(z3,2);
    */
	double finalX = getDeterminant(l1, y1, z1, l2, y2, z2, l3, y3, z3)/(2*det_A) + centerOfSphere[3].getX();
	double finalY = getDeterminant(x1, l1, z1, x2, l2, z2, x3, l3, z3)/(2*det_A) + centerOfSphere[3].getY();
	double finalZ = getDeterminant(x1, y1, l1, x2, y2, l2, x3, y3, l3)/(2*det_A) + centerOfSphere[3].getZ();

	double finalR = sqrt(   ( (finalX - centerOfSphere[0].getX()) * (finalX - centerOfSphere[0].getX()) ) 
                          + ( (finalY - centerOfSphere[0].getY()) * (finalY - centerOfSphere[0].getY()) ) 
                          + ( (finalZ - centerOfSphere[0].getZ()) * (finalZ - centerOfSphere[0].getZ()) ) ) 
                    - arrSphere[3]->getRadius();

    tangentSpheres = new Sphere[1];

    tangentSpheres[0].setSphere( rg_Point3D(finalX, finalY, finalZ), finalR);

	return 1;
}





rg_REAL rg_SphereSetVoronoiDiagram::getDeterminant( const rg_REAL& mat11, const rg_REAL& mat12, const rg_REAL& mat13,
        					                        const rg_REAL& mat21, const rg_REAL& mat22, const rg_REAL& mat23,
		        			                        const rg_REAL& mat31, const rg_REAL& mat32, const rg_REAL& mat33 )
{
	rg_REAL det_Value =   mat11 * ( ( mat22 * mat33 ) - ( mat23 * mat32 ) )
	 				    - mat12 * ( ( mat21 * mat33 ) - ( mat23 * mat31 ) ) 
				        + mat13 * ( ( mat21 * mat32 ) - ( mat22 * mat31 ) );
	return det_Value;
}

    //  VORONOI EDGE COMPUTATION


rg_Point3D rg_SphereSetVoronoiDiagram::getTangentVectorOfEdgeAtVertex(const rg_Point3D& vertexPt,
																const Sphere& s1,
																const Sphere& s2,
																const Sphere& s3 )
{
	rg_Point3D v1 = (s1.getCenter() - vertexPt).getUnitVector();
	rg_Point3D v2 = (s2.getCenter() - vertexPt).getUnitVector();
	rg_Point3D v3 = (s3.getCenter() - vertexPt).getUnitVector();

	rg_REAL x1 = v1.getX();
	rg_REAL x2 = v2.getX();
	rg_REAL x3 = v3.getX();

	rg_REAL y1 = v1.getY();
	rg_REAL y2 = v2.getY();
	rg_REAL y3 = v3.getY();

	rg_REAL z1 = v1.getZ();
	rg_REAL z2 = v2.getZ();
	rg_REAL z3 = v3.getZ();

	rg_REAL x = 
		(y3*(z1 - z2) + y1*(z2 - z3) + 
        y2*(-z1 + z3))/
      (x1*y2 - x1*y3 - y2*z1 + 
        y3*z1 - x1*z2 + y1*z2 - 
        y3*z2 + 
        x3*(y1 - y2 - z1 + z2) + 
        x2*(-y1 + y3 + z1 - z3) + 
        x1*z3 - y1*z3 + y2*z3);
	rg_REAL y = 
		(x3*(-z1 + z2) + 
        x2*(z1 - z3) + x1*(-z2 + z3))
       /(x1*y2 - x1*y3 - y2*z1 + 
        y3*z1 - x1*z2 + y1*z2 - 
        y3*z2 + 
        x3*(y1 - y2 - z1 + z2) + 
        x2*(-y1 + y3 + z1 - z3) + 
        x1*z3 - y1*z3 + y2*z3);
	rg_REAL z = 
		-((x3*(-y1 + y2) + 
          x2*(y1 - y3) + 
          x1*(-y2 + y3))/
        ((x2 - x3 - y2 + y3)*
           (-x1 + x3 + z1 - z3) + 
          (-x1 + x3 + y1 - y3)*
           (-x2 + x3 + z2 - z3)));

	rg_Point3D tangentVector(x, y, z);
	tangentVector.normalize();
	return tangentVector;
}




rg_FLAG rg_SphereSetVoronoiDiagram::evaluateVoronoiEdge(       rg_RBzCurve3D*& curve,
                                                         const rg_Point3D&     startPt,
                                                         const rg_Point3D&     endPt,
                                                         const Sphere&         s1,
                                                         const Sphere&         s2,
                                                         const Sphere&         s3 )
{
	if (    rg_EQ( s1.getRadius(), s2.getRadius() )
		 && rg_EQ( s2.getRadius(), s3.getRadius() )
		 && rg_EQ( s3.getRadius(), s1.getRadius() ) )
	{
		curve = new rg_RBzCurve3D(rg_QUADRATIC); 

		curve->setCtrlPt(0, startPt);
		curve->setCtrlPt(1, (startPt+endPt)/2);
		curve->setCtrlPt(2, endPt);

		curve->setWeight(0, 1);
		curve->setWeight(1, 1);
		curve->setWeight(2, 1);

        return LINEAR_TYPE_EDGE;
    }

  	if ( startPt == endPt )
	{
		curve = new rg_RBzCurve3D(rg_QUADRATIC); 

		curve->setCtrlPt(0, startPt);
		curve->setCtrlPt(1, (startPt+endPt)/2);
		curve->setCtrlPt(2, endPt);

		curve->setWeight(0, 1);
		curve->setWeight(1, 1);
		curve->setWeight(2, 1);

        return LINEAR_TYPE_EDGE;
    }


	//two end points
//	rg_Point3D startPt = voronoiEdge->getStartVertex()->getPoint();
//	rg_Point3D endPt   = voronoiEdge->getEndVertex()->getPoint();

	//two tangent vectors
	rg_Point3D tangentVectorStart = getTangentVectorOfEdgeAtVertex(startPt, s1, s2, s3);
	rg_Point3D tangentVectorEnd = getTangentVectorOfEdgeAtVertex(endPt, s1, s2, s3);

	//to deternmine valid direction of tangent vector at start vertex
	//by checking the orientations of generators
	rg_TMatrix3D rotateStart2ZAxis;
	rotateStart2ZAxis.rotate(tangentVectorStart, rg_Point3D(0,0,1));

	//projected vectors of 3 vectors starting from start vertex end to generators
	rg_Point3D vectorStart1 = (s1.getCenter() - startPt).getUnitVector();
	rg_Point3D vectorStart2 = (s2.getCenter() - startPt).getUnitVector();
	rg_Point3D vectorStart3 = (s3.getCenter() - startPt).getUnitVector();

	rg_Point3D projVec1 = rotateStart2ZAxis*vectorStart1;
	rg_Point3D projVec2 = rotateStart2ZAxis*vectorStart2;
	rg_Point3D projVec3 = rotateStart2ZAxis*vectorStart3;

	rg_Point2D vec1, vec2;
	vec1.setPoint( (projVec2 - projVec1).getX(), (projVec2 - projVec1).getY() );
	vec2.setPoint( (projVec3 - projVec2).getX(), (projVec3 - projVec2).getY() );

	if( vec1*vec2 < 0 ) //opposite direction
	{
		tangentVectorStart = -tangentVectorStart;
	}

	//a passing point
	rg_Point3D passingPt;

	rg_Point3D vec21 = s2.getCenter()-s1.getCenter();
	rg_Point3D vec31 = s3.getCenter()-s1.getCenter();
	rg_Point3D planeNormal = vec21.crossProduct(vec31);

	rg_TMatrix3D transformMat;   //translate and then rotate
	transformMat.translate( -(s1.getCenter()) ); //translate 
	transformMat.rotate(planeNormal, rg_Point3D(0,0,1)); //rotate from planeNormal to z-axis
	
	rg_TMatrix3D invTransformMat;
	invTransformMat.rotate(rg_Point3D(0,0,1), planeNormal);
	invTransformMat.translate( s1.getCenter() );

	rg_Circle2D transformedCircle1((transformMat*s1.getCenter()).getX(),
								   (transformMat*s1.getCenter()).getY(),
								   s1.getRadius() );
	rg_Circle2D transformedCircle2((transformMat*s2.getCenter()).getX(),
								   (transformMat*s2.getCenter()).getY(),
								   s2.getRadius() );
	rg_Circle2D transformedCircle3((transformMat*s3.getCenter()).getX(),
								   (transformMat*s3.getCenter()).getY(),
								   s3.getRadius() );
	rg_INT numOfCC = 0;
	rg_Circle2D tangentCircle1, tangentCircle2;
	numOfCC = makeCircumcircle(transformedCircle1,transformedCircle2,transformedCircle3,
							   tangentCircle1, tangentCircle2);

	if( numOfCC == 1 )
	{
		passingPt.setPoint(tangentCircle1.getCenterPt().getX(),
						   tangentCircle1.getCenterPt().getY(),
						   0);
		passingPt = invTransformMat*passingPt;

	}
	else if( numOfCC == 2)	//\C7ذ\E1\C7ؾ\DF\C7ϴ\C2 \B9\AE\C1\A6! (\BE\C6\C1\F7\C0\BA \C7ذ\E1\C7\CF\C1\F6 \BE\CA\C0\BD-\B5\BF\BF\ED)
	{
		passingPt.setPoint(tangentCircle1.getCenterPt().getX(),
						   tangentCircle1.getCenterPt().getY(),
						   0);
		passingPt = invTransformMat*passingPt;
	}
	else // numOfCC == 0
	{
	}
	//we have found a passing point


 //another method for computing a passing point
 // 	passingPt = evaluate_Voronoi_Vertex(&s1,&s2,&s3,&Sphere(startPt,0)).getCenterPoint();

	curve = new rg_RBzCurve3D(rg_QUADRATIC); 

	rg_REAL area = (endPt-startPt).getUnitVector().crossProduct(tangentVectorStart).magnitude();

	if( rg_ZERO(area, resNeg4) )	//direction of start tangent == vector (start to end)
	{
		curve->setCtrlPt(0, startPt);
		curve->setCtrlPt(1, (startPt+endPt)/2);
		curve->setCtrlPt(2, endPt);

		curve->setWeight(0, 1);
		curve->setWeight(1, 1);
		curve->setWeight(2, 1);
	}
	else
	{
        //x-y plane\BF\A1 projection
		rg_TMatrix3D rotation;
        rg_Point3D baseVec = (endPt-startPt).crossProduct(tangentVectorStart);
        baseVec.normalize();
        const rg_FLAG bReverse=rg_EQ(baseVec.innerProduct(rg_Point3D(0,0,1)), -1);
        if( !bReverse )
	        rotation.rotate(baseVec,rg_Point3D(0,0,1));
        else
            rotation.rotateY(rg_PI);

		rg_TMatrix3D translation;
		translation.translate( -startPt );
		rg_TMatrix3D projection = rotation*translation;		//translate and rotate

		//rg_TMatrix3D rotation;
		//rotation.rotate((endPt-startPt).crossProduct(tangentVectorStart), rg_Point3D(0,0,1));
		//rg_TMatrix3D translation;
		//translation.translate( -startPt );
		//rg_TMatrix3D projection = rotation*translation;		//translate and rotate
		
		rg_Point3D newStartPt = projection*startPt;
		rg_Point3D newStartTangent = rotation*tangentVectorStart; //vector\B4\C2 translation\C7ϸ\E9 error!!!
		rg_Point3D newEndPt = projection*endPt;
		rg_Point3D newEndTangent = rotation*tangentVectorEnd; //vector\B4\C2 translation\C7ϸ\E9 error!!! (only rotation)
		rg_Point3D newPassingPt = projection*passingPt;

		//because of the numerical errors of rotation, we make z-value zero by force. 	
		newStartTangent.setZ(0);
		newEndPt.setZ(0);
		newEndTangent.setZ(0);
		newPassingPt.setZ(0);
		
		curve->makeConic(newStartPt, newStartTangent, newEndPt, newEndTangent, newPassingPt);

		//back transform of control points
		rg_TMatrix3D backTransform;
        if( !bReverse )
	        backTransform.rotate(rg_Point3D(0,0,1),baseVec);
        else
            backTransform.rotateY(-rg_PI);
	    backTransform.translate( startPt );

		//backTransform.rotate(rg_Point3D(0,0,1), (endPt-startPt).crossProduct(tangentVectorStart));
		//backTransform.translate( startPt );

		curve->setCtrlPt(0, backTransform*curve->getCtrlPt(0));
		curve->setCtrlPt(1, backTransform*curve->getCtrlPt(1));
		curve->setCtrlPt(2, backTransform*curve->getCtrlPt(2));

		rg_REAL crossProductTangentVector = newEndPt.crossProduct(newStartTangent).getZ();
		rg_REAL crossProductPassingPoint = newEndPt.crossProduct(newPassingPt).getZ();

		//if passing point is not on the desired segment, we shift the second weight.
		if( crossProductTangentVector*crossProductPassingPoint < 0 )
		{
			curve->setWeight( 1, -curve->getWeight(1) );
		}

		//FOR TEST
		//\B8\B8\BE\E0 orientation\C0\CC \C1\A6\B4\EB\B7\CE \C1־\EE\C1\F8\B4ٸ\E9, \C0\CC code\B4\C2 \BAҷ\C1\C1\F6\C1\F6 \BE\CA\C0\BB \B0\CD\C0̴\D9.
		//\C7\F6\C0\E7\B4\C2 orientation\C0\CC \C1־\EE\C1\F8\B4ٰ\ED \B0\A1\C1\A4\C7\D1 \BB\F3\C5\C2\C0̹Ƿ\CE \BEƷ\A1\C0\C7 code\B0\A1 \C7ʿ\E4\C7ϴ\D9.
		//hyperbola\C0\CF \B0\E6\BF\EC weight\B0\A1 \C0\BD\BC\F6\C0\CE \B0\E6\BF\EC\B4\C2 \B9߻\FD\C7\CF\C1\F6 \BEʴ´\D9.
		if( curve->getWeight(1) < -1 ) // hyperbola
		{
			curve->setWeight( 1, -curve->getWeight(1) );
		}
		//end of FOR TEST
		
	}

    return NONLINEAR_TYPE_EDGE;
}



//this function returns the number of circumcircles (0, 1, or 2)
rg_INT rg_SphereSetVoronoiDiagram::makeCircumcircle(const rg_Circle2D& circle1,
									 const rg_Circle2D& circle2,
									 const rg_Circle2D& circle3,
								     rg_Circle2D&	result1, 
								     rg_Circle2D&	result2)
{
	rg_INT numOfCC = 0;

	rg_ImplicitEquation tangentLine1(1);
	rg_ImplicitEquation tangentLine2(1);

	rg_Point2D z3;
	z3 = makeTangentLinesInWPlane(circle1, circle2, circle3, tangentLine1, tangentLine2);

	rg_REAL smallest = 0;
	if( z3 == circle1.getCenterPt() )
	{
		smallest = circle1.getRadius();
	}
	else if ( z3 == circle2.getCenterPt() )
	{
		smallest = circle2.getRadius();
	}
	else
	{
		smallest = circle3.getRadius();
	}

	//\C7ϳ\AA\C0\CE\C1\F6 \B5ΰ\B3\C0\CE\C1\F6 \C6Ǵ\DC (point location problem)
	rg_REAL sign1 = tangentLine1.evaluateImpEquation(0, 0);
	rg_REAL sign2 = tangentLine2.evaluateImpEquation(0, 0);

	if( sign1 * sign2 < 0 ) //alpha region
	{
		numOfCC = 1;
		if( sign1 > 0 )
		{
			result1 = transformW2Z(tangentLine1, z3);		//+ -
			result1.setRadius( result1.getRadius() - smallest );
		}
		else
		{
			result1 = transformW2Z(tangentLine2, z3);		//- +
			result1.setRadius( result1.getRadius() - smallest );
		}
	}
	else if( sign1 >0 && sign2 >0 ) //gamma region
	{
		numOfCC = 2;
		result1 = transformW2Z(tangentLine1, z3);			//+ +
		result1.setRadius( result1.getRadius() - smallest );

		result2 = transformW2Z(tangentLine2, z3);
		result2.setRadius( result2.getRadius() - smallest );
	}
	else if( sign1 < 0 && sign2 < 0 )						//- -
	{
		numOfCC = 0;
	}
	else if( rg_EQ(sign1, 0., resNeg15) && sign2 > 0 )						//0 +
	{
		numOfCC = 1;
		result1 = transformW2Z(tangentLine2, z3);
		result1.setRadius( result1.getRadius() - smallest );
	}
	else if( sign1 > 0 && rg_EQ(sign2, 0., resNeg15) )						//+ 0
	{
		numOfCC = 1;
		result1 = transformW2Z(tangentLine1, z3);
		result1.setRadius( result1.getRadius() - smallest );
	}
	else												//0 0
	{													//0 -
		numOfCC = 0;									//- 0
	}

	return numOfCC;
}

//this member function shrink radii of circles
//and c3Tilde becomes a point
void rg_SphereSetVoronoiDiagram::shrinkCircle(const		rg_Circle2D& c1,
							   const		rg_Circle2D& c2,
							   const		rg_Circle2D& c3,
							   rg_Circle2D&	c1Tilde,
							   rg_Circle2D&	c2Tilde,
							   rg_Point2D&		c3Tilde)
{
	if( c1.getRadius() < c2.getRadius() )
	{
		if( c1.getRadius() < c3.getRadius() )
		{
			c1Tilde.setCircle(c2.getCenterPt(), c2.getRadius() - c1.getRadius());
			c2Tilde.setCircle(c3.getCenterPt(), c3.getRadius() - c1.getRadius());
			c3Tilde = c1.getCenterPt();
		}
		else
		{
			c1Tilde.setCircle(c1.getCenterPt(), c1.getRadius() - c3.getRadius());
			c2Tilde.setCircle(c2.getCenterPt(), c2.getRadius() - c3.getRadius());
			c3Tilde = c3.getCenterPt();
		}
	}
	else
	{
		if( c2.getRadius() < c3.getRadius() )
		{
			c1Tilde.setCircle(c1.getCenterPt(), c1.getRadius() - c2.getRadius());
			c2Tilde.setCircle(c3.getCenterPt(), c3.getRadius() - c2.getRadius());
			c3Tilde = c2.getCenterPt();
		}
		else
		{
			c1Tilde.setCircle(c1.getCenterPt(), c1.getRadius() - c3.getRadius());
			c2Tilde.setCircle(c2.getCenterPt(), c2.getRadius() - c3.getRadius());
			c3Tilde = c3.getCenterPt();
		}
	}
}

rg_Point2D rg_SphereSetVoronoiDiagram::makeTangentLinesInWPlane(const rg_Circle2D&		c1, 
										      const rg_Circle2D&		c2,
											  const rg_Circle2D&		c3,
											  rg_ImplicitEquation&		result1,
											  rg_ImplicitEquation&		result2)
{
	rg_Circle2D c1Tilde, c2Tilde;
	rg_Point2D z3;

	//c3 has not the smallest radius between three circles
	//so we must make c3 to have smallest radius by swapping for convienience
	shrinkCircle(c1, c2, c3, c1Tilde, c2Tilde, z3);

	rg_Circle2D w1, w2;
	w1 = transformZ2W( c1Tilde, z3 );
	w2 = transformZ2W( c2Tilde, z3 );

	//need to be added generating two tangent lines
	makeTangentLinesOf2Circles(w1, w2, result1, result2);

	return z3;
}

//this function returns tangent lines 
void rg_SphereSetVoronoiDiagram::makeTangentLinesOf2Circles(const rg_Circle2D&	circle1,
											 const rg_Circle2D&	circle2,
											 rg_ImplicitEquation&	result1,
											 rg_ImplicitEquation&	result2)
{
	rg_Circle2D w1, w2;

	//the radius of w2 is less than that of w1
	if( circle1.getRadius() < circle2.getRadius() )
	{
		w1 = circle2;
		w2 = circle1;
	}
	else
	{
		w1 = circle1;
		w2 = circle2;
	}

	rg_Point2D c2cVector = w1.getCenterPt() - w2.getCenterPt();
	rg_REAL r = w1.getRadius() - w2.getRadius();
	rg_REAL length = c2cVector.magnitude();
	rg_REAL sine = r / length;
	rg_REAL cosine = sqrt( length*length - r*r ) / length;

	//rotate theta  /  -theta
	rg_Point2D normal1( c2cVector.getX() * cosine - c2cVector.getY() * sine ,
					 c2cVector.getX() * sine + c2cVector.getY() * cosine );
	rg_Point2D normal2( c2cVector.getX() * cosine + c2cVector.getY() * sine ,
					 c2cVector.getY() * cosine - c2cVector.getX() * sine );
	normal1 = normal1.getUnitVector();
	normal2 = normal2.getUnitVector();
	
	//rotate -PI/2  /  PI/2
	normal1.setPoint( normal1.getY() , -1 * normal1.getX() );
	normal2.setPoint( -1 * normal2.getY() , normal2.getX() );

	result1.setCoeff(1, 0, normal1.getX());
	result1.setCoeff(0, 1, normal1.getY());
	result1.setCoeff(0, 0, -1 * normal1.getX() * w2.getCenterPt().getX() 
						   - normal1.getY() * w2.getCenterPt().getY() 
						   + w2.getRadius()   );

	result2.setCoeff(1, 0, normal2.getX());
	result2.setCoeff(0, 1, normal2.getY());
	result2.setCoeff(0, 0, -1 * normal2.getX() * w2.getCenterPt().getX() 
						   - normal2.getY() * w2.getCenterPt().getY() 
						   + w2.getRadius()   );

}

//
// Z = 1 / w + z3
//
rg_Circle2D rg_SphereSetVoronoiDiagram::transformW2Z(const rg_ImplicitEquation& line,
								   const rg_Point2D& z3)
{
	//degenerate case must be considered
	//case c == 0;

	rg_REAL a = line.getCoeff(1, 0);
	rg_REAL b = line.getCoeff(0, 1);
	rg_REAL c = line.getCoeff(0, 0);

	rg_REAL tempX = -a/2./c + z3.getX();
	rg_REAL tempY =  b/2./c + z3.getY();
	rg_REAL rad = 1./2./c;

	return rg_Circle2D(tempX, tempY, rad);
}

//
// W = 1 / (z - z3)
//
rg_Circle2D rg_SphereSetVoronoiDiagram::transformZ2W(const rg_Circle2D&	cTilde,
								   const rg_Point2D&	c3Tilde)
{
	rg_REAL x1 = cTilde.getCenterPt().getX();
	rg_REAL y1 = cTilde.getCenterPt().getY();
	rg_REAL x3 = c3Tilde.getX();
	rg_REAL y3 = c3Tilde.getY();

	rg_REAL p1 = pow(x1 - x3, 2.) + pow(y1 - y3, 2.) - pow(cTilde.getRadius(), 2.);

	return rg_Circle2D( (x1 - x3)/p1 , -(y1 - y3)/p1, cTilde.getRadius() / p1 );
}



rg_INT rg_SphereSetVoronoiDiagram::makeCircumcircleOnCenterPlane(const Sphere& s1,
																 const Sphere& s2,
																 const Sphere& s3,
																 Sphere& result1,
																 Sphere& result2)
{
	rg_Point3D vec21 = s2.getCenter()-s1.getCenter();
	rg_Point3D vec31 = s3.getCenter()-s1.getCenter();
	rg_Point3D planeNormal = vec21.crossProduct(vec31);


    rg_Vector3D zeroVec(0.0, 0.0, 0.0);
    if ( planeNormal == zeroVec )  {
        planeNormal = vec21;
        planeNormal.normalize();

	    rg_TMatrix3D transMat;   //translate and then rotate
        if( !rg_EQ( planeNormal.innerProduct( rg_Point3D(0.0, 0.0, 1.0) ), -1) )
	        transMat.rotate(planeNormal, rg_Point3D(0,0,1)); //rotate from planeNormal to z-axis
        else
            transMat.rotateY(rg_PI);

        transMat.rotateY(rg_PI*0.5);
        planeNormal = transMat*planeNormal;
    }


	rg_TMatrix3D transformMat;   //translate and then rotate
	transformMat.translate( -(s1.getCenter()) ); //translate 
	transformMat.rotate(planeNormal, rg_Point3D(0,0,1)); //rotate from planeNormal to z-axis
	
	rg_TMatrix3D invTransformMat;
	invTransformMat.rotate(rg_Point3D(0,0,1), planeNormal);
	invTransformMat.translate( s1.getCenter() );

    
	rg_Circle2D transformedCircle1((transformMat*s1.getCenter()).getX(),
								   (transformMat*s1.getCenter()).getY(),
								   s1.getRadius() );
	rg_Circle2D transformedCircle2((transformMat*s2.getCenter()).getX(),
								   (transformMat*s2.getCenter()).getY(),
								   s2.getRadius() );
	rg_Circle2D transformedCircle3((transformMat*s3.getCenter()).getX(),
								   (transformMat*s3.getCenter()).getY(),
								   s3.getRadius() );
	rg_INT numOfCC = 0;
	rg_Circle2D tangentCircle1, tangentCircle2;
	numOfCC = makeCircumcircle(transformedCircle1,transformedCircle2,transformedCircle3,
							   tangentCircle1, tangentCircle2);

	result1.setCenter(invTransformMat*rg_Point3D(tangentCircle1.getCenterPt().getX(),tangentCircle1.getCenterPt().getY(),0));
	result2.setCenter(invTransformMat*rg_Point3D(tangentCircle2.getCenterPt().getX(),tangentCircle2.getCenterPt().getY(),0));
	result1.setRadius(tangentCircle1.getRadius());
	result2.setRadius(tangentCircle2.getRadius());

	return numOfCC;
}

//end of donguk
///////////////
///////////////////////////////////////////////////////////////////////////////
//
//  topological operators..


///////////////////////////////////////////////////////////////////////////////
//
//  geometric operators..
Sphere rg_SphereSetVoronoiDiagram::evaluateSmallestTangentSphereFor3Spheres(const rg_Point3D& normalOfPlaneBy3SphereCenters,
                                                                            const Sphere&     sphereInCCW1,
                                                                            const Sphere&     sphereInCCW2,
                                                                            const Sphere&     sphereInCCW3)
{
    rg_TMatrix3D transformMat;   //translate and then rotate
	transformMat.translate( -(sphereInCCW1.getCenter()) ); //translate 

    const rg_FLAG bReverse=rg_EQ(normalOfPlaneBy3SphereCenters.innerProduct(rg_Point3D(0,0,1)), -1);
    if( !bReverse )
	    transformMat.rotate(normalOfPlaneBy3SphereCenters,rg_Point3D(0,0,1));
    else
        transformMat.rotateY(rg_PI);

	rg_TMatrix3D invRotMat;
    if( !bReverse )
	    invRotMat.rotate(rg_Point3D(0,0,1),normalOfPlaneBy3SphereCenters);
    else
        invRotMat.rotateY(-rg_PI);
	invRotMat.translate( sphereInCCW1.getCenter() ); //translate 
	
	rg_Circle2D transformedCircle1((transformMat*sphereInCCW1.getCenter()).getX(),
								   (transformMat*sphereInCCW1.getCenter()).getY(),
								   sphereInCCW1.getRadius() );
	rg_Circle2D transformedCircle2((transformMat*sphereInCCW2.getCenter()).getX(),
								   (transformMat*sphereInCCW2.getCenter()).getY(),
								   sphereInCCW2.getRadius() );
	rg_Circle2D transformedCircle3((transformMat*sphereInCCW3.getCenter()).getX(),
								   (transformMat*sphereInCCW3.getCenter()).getY(),
								   sphereInCCW3.getRadius() );
	rg_INT numOfCC = 0;
	rg_Circle2D tangentCircle1, tangentCircle2;
	numOfCC = makeCircumcircle(transformedCircle1,transformedCircle2,transformedCircle3,
							   tangentCircle1, tangentCircle2);

	if(numOfCC == 1)
	{
        rg_Point3D centerOfSmallestTS( tangentCircle1.getCenterPt() );
        centerOfSmallestTS = invRotMat*centerOfSmallestTS;

        return Sphere( centerOfSmallestTS, tangentCircle1.getRadius());	
	}
	else if(numOfCC == 2)
	{
        rg_Point3D centerOfSmallestTS( tangentCircle1.getCenterPt() );
        centerOfSmallestTS = invRotMat*centerOfSmallestTS;

        return Sphere( centerOfSmallestTS, tangentCircle1.getRadius());	
	}
    else 
    {
        return Sphere();	
    }
}


void rg_SphereSetVoronoiDiagram::fineInternalVoids(const rg_Point3D& minPtOfBoundingBox,
                                                   const rg_Point3D& maxPtOfBoundingBox,
                                                         rg_INT&     numOfInternalVoids,
                                                         VDVertex**&   internalVoids)
{
    rg_dList<VDVertex*> internalVoidsInBoundingBox;

    rg_REAL minX = minPtOfBoundingBox.getX();
    rg_REAL minY = minPtOfBoundingBox.getY();
    rg_REAL minZ = minPtOfBoundingBox.getZ();

    rg_REAL maxX = maxPtOfBoundingBox.getX();
    rg_REAL maxY = maxPtOfBoundingBox.getY();
    rg_REAL maxZ = maxPtOfBoundingBox.getZ();

    VDVertex* currVertex = rg_NULL;
	m_GlobalVertexList.reset4Loop();
	while( m_GlobalVertexList.setNext4Loop() )  
    {
		currVertex = m_GlobalVertexList.getpEntity();

        if ( !currVertex->isConnectedWithUnboundedEdge() )
        {
            rg_REAL centerX = currVertex->getPoint().getX();
            rg_REAL centerY = currVertex->getPoint().getY();
            rg_REAL centerZ = currVertex->getPoint().getZ();
            rg_REAL radius  = currVertex->getRadiusOfTangentSphere();
        
            if (    rg_GT( centerX-radius, minX ) && rg_GT( centerY-radius, minY ) && rg_GT( centerZ-radius, minZ )
                 && rg_LT( centerX+radius, maxX ) && rg_LT( centerY+radius, maxY ) && rg_LT( centerZ+radius, maxZ ) )
            {
                internalVoidsInBoundingBox.addTail( currVertex );
            }
        }
    }

    numOfInternalVoids = internalVoidsInBoundingBox.getSize();
    internalVoids      = internalVoidsInBoundingBox.getArray();

    qsort( (void *) internalVoids, numOfInternalVoids, sizeof(VDVertex*), compareInternalVoids);
}



/*
rg_dList<VDFace*>* rg_SphereSetVoronoiDiagram::findInteractionInterface()
{
    rg_FLAG existPreviousTangibilityCheck = m_hasProbeTangibility;
    rg_REAL prevProbeRadius               = m_probeRadius;


    //rg_REAL dist = 13.;
    rg_REAL dist = 25.;
    if ( existPreviousTangibilityCheck == rg_FALSE || rg_NE(prevProbeRadius, dist) )
        defineProbeTangibility(dist);
    
    rg_dList<VDFace*>* facesWithTwoGeneratorInDifferentGroup = new rg_dList<VDFace*>;

    VDFace* currFace = rg_NULL;

    BallGenerator* rightGenerator = rg_NULL;
    BallGenerator* leftGenerator  = rg_NULL;

    ///////////////////////////////
    //  modified by Y.Cho 2004-10-20
//    PDBAtomAttribute* attributeOfRightGenerator = rg_NULL;
//    PDBAtomAttribute* attributeOfLeftGenerator  = rg_NULL;
    char chainIDOfRightGenerator;
    char chainIDOfLeftGenerator;

	m_GlobalFaceList.reset4Loop();
	while( m_GlobalFaceList.setNext4Loop() )  
    {
		currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isBounded() == rg_FALSE )
            continue;

        rightGenerator = (BallGenerator*) currFace->getRightCell()->getGenerator();
        leftGenerator  = (BallGenerator*) currFace->getLeftCell()->getGenerator();
        
//        attributeOfRightGenerator = (PDBAtomAttribute*) rightGenerator->getAttribute();
//        attributeOfLeftGenerator  = (PDBAtomAttribute*) leftGenerator->getAttribute();

//        if ( attributeOfRightGenerator->getGroupID() != attributeOfLeftGenerator->getGroupID() )
//        {
//            if ( currFace->isTangible() == FULLY_TANGIBLE )
//                continue;

//            if ( currFace->isTangible() == PARTIALLY_TANGIBLE )  {
//                continue;
//            }

//            facesWithTwoGeneratorInDifferentGroup->addTail( currFace );       
//        }
        chainIDOfRightGenerator = ((PDBAtom*)rightGenerator->getProperty())->getChainID();
        chainIDOfLeftGenerator  = ((PDBAtom*)leftGenerator->getProperty())->getChainID();
        if ( chainIDOfRightGenerator != chainIDOfLeftGenerator )
        {
            if ( currFace->isTangible() == FULLY_TANGIBLE )
                continue;

            if ( currFace->isTangible() == PARTIALLY_TANGIBLE )  {
                continue;
            }

            facesWithTwoGeneratorInDifferentGroup->addTail( currFace );       
        }
    }

    if ( existPreviousTangibilityCheck == rg_TRUE )  {
        if ( rg_NE(prevProbeRadius, dist) )
            defineProbeTangibility(prevProbeRadius);
    }
    else  {
        resetTangibility();
    }

    return facesWithTwoGeneratorInDifferentGroup;
}
*/


//  rg_REAL distancesBetweenCenters[4] 
//  rg_REAL distancesBetweenBalls[4]  
//      distances..[0] : sum of distance
//      distances..[1] : average of distance
//      distances..[2] : maxinum distance
//      distances..[3] : mininum distance
void rg_SphereSetVoronoiDiagram::analyzeDistanceBetweenGroups( rg_dList<VDFace*>* facesWithTwoGeneratorInDifferentGroup,
                                                               rg_REAL*           distancesBetweenCenters,
                                                               rg_REAL*           distancesBetweenBalls )
{
    rg_INT numOfFacesWithTwoGeneratorInDifferentGroup = facesWithTwoGeneratorInDifferentGroup->getSize();

    distancesBetweenCenters[0] = 0.0;
    distancesBetweenBalls[0]   = 0.0;
    distancesBetweenCenters[2] = -DBL_MAX;
    distancesBetweenBalls[2]   = -DBL_MAX;
    distancesBetweenCenters[3] = DBL_MAX;
    distancesBetweenBalls[3]   = DBL_MAX;


    Sphere rightBall;
    Sphere leftBall;

    rg_REAL distanceTwoCenter = 0.0;
    rg_REAL distanceTwoBall   = 0.0;

    VDFace* currFace = rg_NULL;
	facesWithTwoGeneratorInDifferentGroup->reset4Loop();
	while( facesWithTwoGeneratorInDifferentGroup->setNext4Loop() )  
    {
		currFace = facesWithTwoGeneratorInDifferentGroup->getEntity();

        rightBall = ((BallGenerator*) currFace->getRightCell()->getGenerator())->getBall();
        leftBall  = ((BallGenerator*) currFace->getLeftCell()->getGenerator())->getBall();

        distanceTwoCenter = rightBall.getCenter().distance( leftBall.getCenter() );
        distanceTwoBall   = distanceTwoCenter - rightBall.getRadius() - leftBall.getRadius();
        if ( rg_NEG( distanceTwoBall ) )
            distanceTwoBall = 0.0;

        distancesBetweenCenters[0] += distanceTwoCenter;

        if ( rg_GT( distanceTwoCenter, distancesBetweenCenters[2] ) )
            distancesBetweenCenters[2] = distanceTwoCenter;

        if ( rg_LT( distanceTwoCenter, distancesBetweenCenters[3] ) )
            distancesBetweenCenters[3] = distanceTwoCenter;


        distancesBetweenBalls[0]   += distanceTwoBall;

        if ( rg_GT( distanceTwoBall, distancesBetweenBalls[2] ) )
            distancesBetweenBalls[2] = distanceTwoBall;

        if ( rg_LT( distanceTwoBall, distancesBetweenBalls[3] ) )
            distancesBetweenBalls[3] = distanceTwoBall;

    }

    distancesBetweenCenters[1] = distancesBetweenCenters[0]/numOfFacesWithTwoGeneratorInDifferentGroup;
    distancesBetweenBalls[1]   = distancesBetweenBalls[0]/numOfFacesWithTwoGeneratorInDifferentGroup;

}

rg_dList<rg_NURBSplineSurface3D>* rg_SphereSetVoronoiDiagram::computeLinkBlends(const rg_REAL& radius)
{
	rg_RevolvedSurface blendSurf;

	rg_dList<rg_NURBSplineSurface3D>* linkBlends = new rg_dList<rg_NURBSplineSurface3D>;
	VDEdge* ptrEdge = NULL;

	m_GlobalEdgeList.reset4Loop();
	while(m_GlobalEdgeList.setNext4Loop())
	{
		ptrEdge = m_GlobalEdgeList.getpEntity();

		if(ptrEdge->getMinRadius() < radius)
		{
			rg_REAL maxRad = 0.;
			switch(ptrEdge->getMinPosition())
			{
			case START:
				//max side check
				if(ptrEdge->getEndVertex()->isOnInfinity() || 
				   ptrEdge->getEndVertex()->getRadiusOfTangentSphere() > radius)
				{
					//construct blend (end vertex)
					computeLinkBlendToEdge(ptrEdge, radius, END, blendSurf);
					linkBlends->add(blendSurf);
				}
				break;
				
			case END:
				if(ptrEdge->getStartVertex()->getRadiusOfTangentSphere() > radius)
				{
					//construct blend (start vertex)
					computeLinkBlendToEdge(ptrEdge, radius, START, blendSurf);
					linkBlends->add(blendSurf);
				}
				break;
			
			case MID:
				//minimum\C0\CC edge\C0\C7 \B0\A1\BF \C1\B8\C0\E7\C7\D2 \B0\E6\BF\EC, \BE\E7 \B9\E6\C7\E2(start/end)\B8\F0\B5ο\A1\BC\AD 
				//blending\C0\CC \BB\FD\B1\E6 \BC\F6 \C0ִ\D9.
				if(ptrEdge->getEndVertex()->isOnInfinity() || 
				   ptrEdge->getEndVertex()->getRadiusOfTangentSphere() > radius)
				{
					//construct blend (end vertex)
					computeLinkBlendToEdge(ptrEdge, radius, END, blendSurf);
					linkBlends->add(blendSurf);
				}
				if(ptrEdge->getStartVertex()->getRadiusOfTangentSphere() > radius)
				{
					//construct blend (start vertex)
					computeLinkBlendToEdge(ptrEdge, radius, START, blendSurf);
					linkBlends->add(blendSurf);
				}
				break;
			}
		}
	}
	return linkBlends;

}

rg_Point3D rg_SphereSetVoronoiDiagram::computeTangentProbe(VDEdge* theEdge,
														   const rg_REAL& radius,
														   const PositionOnEdge& position)
{
	VDCell* cells[3];
	theEdge->inquireIntoOnlyEdgeSharingCellsInCCW(cells);
	
	Sphere s1 = ((BallGenerator*)cells[0]->getGenerator())->getBall();
	Sphere s2 = ((BallGenerator*)cells[1]->getGenerator())->getBall();
	Sphere s3 = ((BallGenerator*)cells[2]->getGenerator())->getBall();

	rg_Point3D vec21 = s2.getCenter()-s1.getCenter();
	rg_Point3D vec31 = s3.getCenter()-s1.getCenter();
	rg_Point3D planeNormal = vec21.crossProduct(vec31);

	rg_TMatrix3D transformMat;   //translate and then rotate
	transformMat.translate( -(s1.getCenter()) ); //translate 
	transformMat.rotate(planeNormal, rg_Point3D(0,0,1)); //rotate from planeNormal to z-axis
	
	rg_TMatrix3D invTransformMat;
	invTransformMat.rotate(rg_Point3D(0,0,1), planeNormal);
	invTransformMat.translate( s1.getCenter() );

	//s1 is mapped to the origin
	//s2 is mapped to (x2,y2,0,r2)
	//s2 is mapped to (x3,y3,0,r3)

	//x1 and y1 is zero
	rg_REAL r1 = s1.getRadius();

	rg_REAL x2 = (transformMat*s2.getCenter()).getX();
	rg_REAL y2 = (transformMat*s2.getCenter()).getY();
	rg_REAL r2 = s2.getRadius();
	
	rg_REAL x3 = (transformMat*s3.getCenter()).getX();
	rg_REAL y3 = (transformMat*s3.getCenter()).getY();
	rg_REAL r3 = s3.getRadius();

	//mapped\B0\F8\B0\A3\BF\A1\BC\AD\C0\C7 probe\C0\C7 \C1߽\C9 (z=0)
	rg_REAL x = -(2.0*radius*r2*y3-x2*x2*y3+y2*x3*x3+2.0*y2*radius*r1+y2*y3*y3-2.0*y2*radius*r3-y2*
				r3*r3+y2*r1*r1-y2*y2*y3+r2*r2*y3-2.0*radius*r1*y3-r1*r1*y3)/(x2*y3-x3*y2)/2.0;

	rg_REAL y = (x2*x3*x3+2.0*x2*radius*r1+x2*y3*y3-2.0*x2*radius*r3-x2*r3*r3-x3*r1*r1+x2*r1*
				r1+2.0*x3*radius*r2-x3*x2*x2-x3*y2*y2+x3*r2*r2-2.0*x3*radius*r1)/(x2*y3-x3*y2)/2.0;

	if( (radius+r1)*(radius+r1)-x*x-y*y < 0 )
		return rg_Point3D(0,0,0);

	rg_REAL z = sqrt((radius+r1)*(radius+r1)-x*x-y*y);

	rg_Point3D center(x,y,z);
	rg_Point3D contactPt1 = center - center.getUnitVector()*radius;
	rg_Point3D contactPt2 = center + (rg_Point3D(x2,y2,0) - center).getUnitVector()*radius;
	rg_Point3D contactPt3 = center + (rg_Point3D(x3,y3,0) - center).getUnitVector()*radius;
	rg_Point3D vector1 = contactPt1 - contactPt2;
	rg_Point3D vector2 = contactPt2 - contactPt3;
	rg_REAL zValue = (vector1*vector2).getZ();

	switch(position)
	{
	case START:
		if( zValue > 0 )
		{
			z = -z;
		}
		break;
	case END:
		if( zValue < 0 )
		{
			z = -z;
		}
		break;
	}
	
	return invTransformMat*rg_Point3D(x,y,z);
}
						 

void rg_SphereSetVoronoiDiagram::computeLinkBlendToEdge(VDEdge* theEdge, 
														const rg_REAL& radius, 
														const PositionOnEdge& position, //START or END
														rg_RevolvedSurface& computedSurf)
{
	VDCell* cells[3];
	theEdge->inquireIntoOnlyEdgeSharingCellsInCCW(cells);
	
	Sphere s1 = ((BallGenerator*)cells[0]->getGenerator())->getBall();
	Sphere s2 = ((BallGenerator*)cells[1]->getGenerator())->getBall();
	Sphere s3 = ((BallGenerator*)cells[2]->getGenerator())->getBall();


	//center of probe at the original space
	rg_Point3D center = computeTangentProbe(theEdge, radius, position);

	rg_Point3D contactPt1;
	rg_Point3D contactPt2;
	rg_Point3D contactPt3;
	if(position == START) //\B9\E6\C7\E2\C0\CC \B9ٲ\F0 \B0\E6\BF\EC orientation\C0\BB \B9ٲ\E3\C1־\EE\BE\DF \B0\A1\BD\C3ȭ\BD\C3 normal\C0\CC \C0\CF\C1\A4\C7\D1 \B9\E6\C7\E2\C0\B8\B7\CE \B3\AA\BF´\D9.
	{
		contactPt1 = center + (s1.getCenter() - center).getUnitVector()*radius;
		contactPt3 = center + (s2.getCenter() - center).getUnitVector()*radius;
		contactPt2 = center + (s3.getCenter() - center).getUnitVector()*radius;
	}
	else
	{
		contactPt1 = center + (s1.getCenter() - center).getUnitVector()*radius;
		contactPt2 = center + (s2.getCenter() - center).getUnitVector()*radius;
		contactPt3 = center + (s3.getCenter() - center).getUnitVector()*radius;
	}

	rg_Point3D vec1 = contactPt1 - center;
	rg_Point3D vec2 = contactPt2 - center;
	rg_Point3D vec3 = contactPt3 - center;

	rg_TMatrix3D rotMat1;
	rotMat1.rotateArbitraryAxis( vec1.crossProduct(vec2), rg_PI/2 );
	rg_RBzCurve3D rbzCurve3(2);
	makeConicByProjection(contactPt1,rotMat1*vec1,contactPt2, rotMat1*vec2, center + (vec1+vec2).getUnitVector()*radius,rbzCurve3);

	rg_TMatrix3D rotMat2;
	rotMat2.rotateArbitraryAxis( vec2.crossProduct(vec3), rg_PI/2 );
	rg_RBzCurve3D rbzCurve1(2);
	makeConicByProjection(contactPt2, rotMat2*vec2, contactPt3, rotMat2*vec3, center + (vec2+vec3).getUnitVector()*radius,rbzCurve1);

	rg_TMatrix3D rotMat3;
	rotMat3.rotateArbitraryAxis( vec3.crossProduct(vec1), rg_PI/2 );
	rg_RBzCurve3D rbzCurve2(2);
	makeConicByProjection(contactPt3, rotMat3*vec3, contactPt1, rotMat3*vec1, center + (vec1+vec3).getUnitVector()*radius,rbzCurve2);

	computedSurf.setOrderOfSurface(3,3);
	computedSurf.setControlNetAndWeight(3,3);

	computedSurf.setPointOnControlNet(0,0,contactPt1);
	computedSurf.setPointOnControlNet(0,1,contactPt1);
	computedSurf.setPointOnControlNet(0,2,contactPt1);
	computedSurf.setWeight(0,0,1);
	computedSurf.setWeight(0,1,1);
	computedSurf.setWeight(0,2,1);

	computedSurf.setPointOnControlNet(1,0,rbzCurve3.getCtrlPt(1));
	computedSurf.setPointOnControlNet(1,1,((contactPt1+contactPt2+contactPt3)/3-center).getUnitVector()*radius*1.1+center);
	computedSurf.setPointOnControlNet(1,2,rbzCurve2.getCtrlPt(1));
	computedSurf.setWeight(1,0,rbzCurve3.getWeight(1));
	computedSurf.setWeight(1,1,0.2);
	computedSurf.setWeight(1,2,rbzCurve2.getWeight(1));

	computedSurf.setPointOnControlNet(2,0,contactPt2);
	computedSurf.setPointOnControlNet(2,1,rbzCurve1.getCtrlPt(1));
	computedSurf.setPointOnControlNet(2,2,contactPt3);
	computedSurf.setWeight(2,0,1);
	computedSurf.setWeight(2,1,rbzCurve1.getWeight(1));
	computedSurf.setWeight(2,2,1);

	rg_REAL knots[6] = {0,0,0,1,1,1};
	computedSurf.setKnotVectorOfU(6,knots);
	computedSurf.setKnotVectorOfV(6,knots);

}

rg_dList<rg_NURBSplineSurface3D>* rg_SphereSetVoronoiDiagram::computeRollingBlends(const rg_REAL& radius)
{
	int count = 0;

	rg_dList<rg_NURBSplineSurface3D>* rollingBlends = new rg_dList<rg_NURBSplineSurface3D>;

	Sphere lSphere, rSphere;
	VDFace* ptrFace = NULL;
	m_GlobalFaceList.reset4Loop();
	while(m_GlobalFaceList.setNext4Loop())
	{
		ptrFace = m_GlobalFaceList.getpEntity();
		lSphere = ((BallGenerator*)(ptrFace->getLeftCell()->getGenerator()))->getBall();
		rSphere = ((BallGenerator*)(ptrFace->getRightCell()->getGenerator()))->getBall();
		rg_REAL dist = lSphere.getCenter().distance(rSphere.getCenter())-lSphere.getRadius()-rSphere.getRadius();
		if(dist < radius*2)
		{
			computeRollingBlendToFace(ptrFace, radius, rollingBlends);
		}

//		if(count++ > 200)
//			break;
	}
	return rollingBlends;
}

void rg_SphereSetVoronoiDiagram::computeRollingBlendToFace(VDFace* ptrFace, 
														   const rg_REAL& radius, 
														   rg_dList<rg_NURBSplineSurface3D>* rollingBlends)
{
	//BASIC IDEA:
	//probe\C0\C7 center point\B8\A6 ccw\B9\E6\C7\E2\C0\B8\B7\CE ã\B4´\D9. 
	//ó\C0\BD \BD\C3\C0۵Ǵ\C2 center point\B4\C2 blending\C0\CC \BD\C3\C0۵Ǵ\C2 \C1\A1\C0\CC \B5ȴ\D9.
	//\B4\D9\C0\BD\C0\B8\B7\CE \B8\B8\B3\AA\B4\C2 center point\B4\C2 blending\C0\CC \B3\A1\B3\AA\B4\C2 point\B0\A1 \B5ȴ\D9.

	Sphere lSphere = ((BallGenerator*)(ptrFace->getLeftCell()->getGenerator()))->getBall();
	Sphere rSphere = ((BallGenerator*)(ptrFace->getRightCell()->getGenerator()))->getBall();

	rg_dList<rg_Point3D> centerList;
	VDLoop* loop = ptrFace->getLoop(0);
	VDPartialEdge* firstEdge = loop->getPartialEdge();

	VDPartialEdge* currEdge = firstEdge;
	VDEdge* original = NULL;

	rg_FLAG isStart = rg_FALSE;

	do
	{
		original = currEdge->getOriginalEdge();
		
		if( original->getMaxRadius() < radius )
		{//\BB\FD\B1\E2\C1\F6 \BEʴ´\D9.
			currEdge = currEdge->getNextPartialEdgeInLoop();
			continue;
		}
		if( original->getMinRadius() > radius )
		{//\C0\FCü\BF\A1\BC\AD \BB\FD\B1\E4\B4\D9.
			currEdge = currEdge->getNextPartialEdgeInLoop();
			continue;
		}

		switch(original->getMinPosition())
		{
		case START:
			if(currEdge->isRightOrientationInLoop())
			{
				//add to list (start of blend) --> end \C2ʿ\A1 \BB\FD\B1\E4\B4\D9.
				if( centerList.isEmpty() )
					isStart = rg_TRUE;
				centerList.add( computeTangentProbe(original,radius,END) );
			}
			else
			{
				//add to list (end of blend) --> end\B9\E6\C7⿡\BC\AD \B3\A1\B3\AD\B4\D9.
				centerList.add( computeTangentProbe(original,radius,END) );
			}
			break;
		case END:
			if(currEdge->isRightOrientationInLoop())
			{
				//end of blend --> start\B9\E6\C7⿡\BC\AD \B3\A1\B3\B2
				centerList.add( computeTangentProbe(original,radius,START) );
			}
			else
			{
				//start of blend --> start \B9\E6\C7⿡\BC\AD \BD\C3\C0۵\CA.
				if( centerList.isEmpty() )
					isStart = rg_TRUE;
				centerList.add( computeTangentProbe(original,radius,START) );
			}
			break;
		case MID:
			if(currEdge->isRightOrientationInLoop())
			{
				if(original->getStartVertex()->getRadiusOfTangentSphere() > radius)
				{
					//end of blend --> start \B9\E6\C7\E2
					centerList.add( computeTangentProbe(original,radius,START) );
				}
				if(original->getEndVertex()->isOnInfinity() ||
				   original->getEndVertex()->getRadiusOfTangentSphere() > radius)
				{
					//start of blend --> end \B9\E6\C7\E2
					if( centerList.isEmpty() )
						isStart = rg_TRUE;
					centerList.add( computeTangentProbe(original,radius,END) );
				}
			}
			else
			{
				//\BEƷ\A1\C0\C7 \BC\F8\BC\AD\B5\B5 \BF\B5\C7\E2\C0\BB \B9\CCģ\B4\D9.
				if(original->getEndVertex()->isOnInfinity() ||
				   original->getEndVertex()->getRadiusOfTangentSphere() > radius)
				{
					//end of blend --> end \B9\E6\C7\E2
					centerList.add( computeTangentProbe(original,radius,END) );
				}
				if(original->getStartVertex()->getRadiusOfTangentSphere() > radius)
				{
					//start of blend --> start\B9\E6\C7\E2
					if( centerList.isEmpty() )
						isStart = rg_TRUE;
					centerList.add( computeTangentProbe(original,radius,START) );
				}
			}
			break;
		}

		currEdge = currEdge->getNextPartialEdgeInLoop();
	}
	while( currEdge != firstEdge );

	rg_FLAG fullExist = rg_FALSE;
	if( centerList.isEmpty() )
	{
		if(original->getMinRadius() > radius)
			fullExist = rg_TRUE;
	}

	if( !centerList.isEmpty() && !isStart )
	{
		centerList.setHead( centerList.getTail() );
	}


	if( fullExist )
	{
		rg_TMatrix3D mat;
		mat.rotate(rg_Point3D(1,0,0), rSphere.getCenter() - lSphere.getCenter());
		mat.translate(lSphere.getCenter());

		rg_REAL r1 = lSphere.getRadius();
		rg_REAL r2 = rSphere.getRadius();
		rg_REAL x2 = (rSphere.getCenter()-lSphere.getCenter()).magnitude();

		rg_REAL x = (2.0*radius*r1+x2*x2-2.0*radius*r2-r2*r2+r1*r1)/x2/2.0;
		rg_REAL y = sqrt(-x*x+radius*radius+2.0*radius*r1+r1*r1);

		rg_Point3D center = mat*rg_Point3D(x,y,0);
		rg_Point3D contactPt1 = center + (lSphere.getCenter() - center).getUnitVector()*radius;
		rg_Point3D contactPt2 = center + (rSphere.getCenter() - center).getUnitVector()*radius;

		rg_Point3D vec1(contactPt1 - center);
		rg_Point3D vec2(contactPt2 - center);
		rg_Line3D axis(lSphere.getCenter(),rSphere.getCenter());

		if(y < radius)
		{
			rg_REAL dist = sqrt( radius*radius - y*y );
			rg_Point3D ptOnAxis1 = mat*rg_Point3D(x-dist,0,0);
			rg_Point3D ptOnAxis2 = mat*rg_Point3D(x+dist,0,0);

			rg_Outline3D circularCrv1;
			rg_Outline3D circularCrv2;

			circularCrv1.formArc(center, contactPt2, ptOnAxis2, vec2*vec1);
			circularCrv2.formArc(center, ptOnAxis1, contactPt1, vec2*vec1);

			rg_RevolvedSurface rollSurf1;
			rollSurf1.setGeneratrix(circularCrv1);
			rollSurf1.setAxis(axis);
			rollSurf1.setAngle(2*rg_PI);
			rollSurf1.makeSurface();
			rollingBlends->add(rollSurf1);

			rg_RevolvedSurface rollSurf2;
			rollSurf2.setGeneratrix(circularCrv2);
			rollSurf2.setAxis(axis);
			rollSurf2.setAngle(2*rg_PI);
			rollSurf2.makeSurface();
			rollingBlends->add(rollSurf2);
		}
		else
		{
			rg_Outline3D circularCrv;
			circularCrv.formArc(center,contactPt2,contactPt1,vec2*vec1);

			rg_RevolvedSurface rollSurf;
			rollSurf.setGeneratrix(circularCrv);
			rollSurf.setAxis(axis);
			rollSurf.setAngle(2*rg_PI);
			rollSurf.makeSurface();
			rollingBlends->add(rollSurf);
		}
	}
	else //no exist & partially exist
	{
		rg_Point3D startCenter, endCenter;
		rg_Point3D contactPt1, contactPt2;
		centerList.reset4Loop();
		while(centerList.setNext4Loop())
		{
			startCenter = centerList.getEntity();
			centerList.setNext4Loop();
			endCenter = centerList.getEntity();

			contactPt1 = startCenter + (lSphere.getCenter() - startCenter).getUnitVector()*radius;
			contactPt2 = startCenter + (rSphere.getCenter() - startCenter).getUnitVector()*radius;

			rg_Point3D vec1(contactPt1 - startCenter);
			rg_Point3D vec2(contactPt2 - startCenter);

			//axis of revolution
			rg_Line3D axis(lSphere.getCenter(),rSphere.getCenter());
			rg_REAL dist2axis = axis.getDistance(startCenter);

			if( dist2axis < radius ) //then, the surface has self intersection
			{
				rg_Point3D unitDir = axis.getUnitDirection();
				rg_Point3D target = startCenter - axis.getSP();
				rg_Point3D footOnAxis = target.innerProduct(unitDir)*unitDir + axis.getSP();

				rg_Point3D transVec = sqrt( radius*radius - dist2axis*dist2axis )*unitDir;

				rg_Point3D ptOnAxis1 = footOnAxis -transVec;
				rg_Point3D ptOnAxis2 = footOnAxis +transVec;

				//\B5\CE \B0\B3\C0\C7 \BF\F8ȣ\B7\CE \B3\AA\B4\A9\BE\EE\C1\F8\B4\D9.
				rg_Outline3D circularCrv1;
				rg_Outline3D circularCrv2;

				circularCrv1.formArc(startCenter,contactPt2,ptOnAxis2,vec2*vec1);
				circularCrv2.formArc(startCenter,ptOnAxis1,contactPt1,vec2*vec1);

				rg_REAL dist = startCenter.distance(endCenter);
				rg_REAL angle = 0.;

				rg_Point3D vectorS = startCenter-footOnAxis;
				rg_Point3D vectorE = endCenter-footOnAxis;

				if((vectorS*vectorE).innerProduct(unitDir) > 0)
				{
					angle = acos( 1. - (dist*dist)/(2.*dist2axis*dist2axis) );
				}
				else
				{
					angle = 2*rg_PI - acos( 1. - (dist*dist)/(2.*dist2axis*dist2axis) );
				}

				rg_RevolvedSurface rollSurf1;
				rollSurf1.setGeneratrix(circularCrv1);
				rollSurf1.setAxis(axis);
				rollSurf1.setAngle(angle);
				rollSurf1.makeSurface();
				rollingBlends->add(rollSurf1);

				rg_RevolvedSurface rollSurf2;
				rollSurf2.setGeneratrix(circularCrv2);
				rollSurf2.setAxis(axis);
				rollSurf2.setAngle(angle);
				rollSurf2.makeSurface();
				rollingBlends->add(rollSurf2);

			}
			else
			{
				rg_Outline3D circularCrv;
				circularCrv.formArc(startCenter,contactPt2,contactPt1,vec2*vec1);

				rg_REAL dist = startCenter.distance(endCenter);
				rg_REAL angle = 0.;

				rg_Point3D unitDir = axis.getUnitDirection();
				rg_Point3D target = startCenter - axis.getSP();
				rg_Point3D footOnAxis = target.innerProduct(unitDir)*unitDir + axis.getSP();
				rg_Point3D vectorS = startCenter-footOnAxis;
				rg_Point3D vectorE = endCenter-footOnAxis;

				if((vectorS*vectorE).innerProduct(unitDir) > 0)
				{
					angle = acos( 1. - (dist*dist)/(2.*dist2axis*dist2axis) );
				}
				else
				{
					angle = 2*rg_PI - acos( 1. - (dist*dist)/(2.*dist2axis*dist2axis) );
				}

				rg_RevolvedSurface rollSurf;

				rollSurf.setGeneratrix(circularCrv);
				rollSurf.setAxis(axis);
				rollSurf.setAngle(angle);
				rollSurf.makeSurface();

				rollingBlends->add(rollSurf);
			}
		}
	}
}

//passPt must lie on the desired curve segment.
void rg_SphereSetVoronoiDiagram::makeConicByProjection(const rg_Point3D& sPt, const rg_Point3D& tangentStart,
													   const rg_Point3D& ePt, const rg_Point3D& tangentEnd,
													   const rg_Point3D& passPt,
													   rg_RBzCurve3D& curve)
{
		//x-y plane\BF\A1 projection
		rg_TMatrix3D rotation;
		rotation.rotate((ePt-sPt).crossProduct(tangentStart), rg_Point3D(0,0,1));
		rg_TMatrix3D translation;
		translation.translate( -sPt );
		rg_TMatrix3D projection = rotation*translation;		//translate and rotate
		
		rg_Point3D newStartPt = projection*sPt;
		rg_Point3D newStartTangent = rotation*tangentStart; //vector\B4\C2 translation\C7ϸ\E9 error!!!
		rg_Point3D newEndPt = projection*ePt;
		rg_Point3D newEndTangent = rotation*tangentEnd; //vector\B4\C2 translation\C7ϸ\E9 error!!! (only rotation)
		rg_Point3D newPassingPt = projection*passPt;

		//because of the numerical errors of rotation, we make z-value zero by force. 	
		newStartTangent.setZ(0);
		newEndPt.setZ(0);
		newEndTangent.setZ(0);
		newPassingPt.setZ(0);
		
		curve.makeConic(newStartPt, newStartTangent, newEndPt, newEndTangent, newPassingPt);

		//back transform of control points
		rg_TMatrix3D backTransform;
		backTransform.rotate(rg_Point3D(0,0,1), (ePt-sPt).crossProduct(tangentStart));
		backTransform.translate( sPt );

		curve.setCtrlPt(0, backTransform*curve.getCtrlPt(0));
		curve.setCtrlPt(1, backTransform*curve.getCtrlPt(1));
		curve.setCtrlPt(2, backTransform*curve.getCtrlPt(2));
}

void rg_SphereSetVoronoiDiagram::defineProbeTangibility(const rg_REAL& radiusOfProbe)
{
    if ( m_hasProbeTangibility )
        resetTangibility();


    //  1. LINK BLEND 
    //      : link blend\B4\C2 voronoi edge\B8\A6 \C1\A4\C0\C7\C7ϴ\C2 3 generator\B5鿡 \C0\C7\C7\D8 \C1\A4\C0ǵ\CA. 
    //        voronoi edge\C0\C7 tangibility\B8\A6 \B0˻\F6\C7ϸ\E9 link blend\C0\C7 \C1\B8\C0\E7\C0\AF\B9\AB\B8\A6 \BE˼\F6 \C0\D6\C0\BD.
    //      1) fully-tangible edge    : probe\B0\A1 edge\C0\A7\C0\C7 \B8\F0\B5\E7 \C1\A1\B5\E9\C0\BB \C1\F6\B3\AF \BC\F6 \C0ִ\C2 edge.
    //      2) partially-tangible edge: probe\B0\A1 edge\C0\A7\C0\C7 \C0Ϻ\CE \C1\A1\B5鸸\C0\BB \C1\F6\B3\AF \BC\F6 \C0ִ\C2 edge.
    //      3) non-tangible edge      : probe\B0\A1 edge\C0\A7\C0\C7 \BE \C1\A1\B5\B5 \C1\F6\B3\AF \BC\F6 \BE\F8\B4\C2 edge.
    //
    defineProbeTangibleEdges( radiusOfProbe );

    //  2. ROLLING BLEND 
    //      : rolling blend\B4\C2 voronoi face\B8\A6 \C1\A4\C0\C7\C7ϴ\C2 2 generator\B5鿡 \C0\C7\C7\D8 \C1\A4\C0ǵ\CA.
    //        voronoi edge\C0\C7 tangibility\B8\A6 \B0˻\F6\C7ϸ\E9 rolling blend\C0\C7 \C1\B8\C0\E7\C0\AF\B9\AB\B8\A6 \BE˼\F6 \C0\D6\C0\BD.
    //      1) fully-tangible edge\B8\A6 \C1\A4\C0\C7\C7ϴ\C2 3 generator\B5\E9\C0\BA \BC\AD\B7\CE rolling blend\B8\A6 \B0\A1\C1\FA\BC\F6 \C0\D6\C0\BD.
    //      2) partially-tangible edge\B8\A6 \C1\A4\C0\C7\C7ϴ\C2 3 generator\B5\E9\C0\BA \C7׻\F3 \BC\AD\B7\CE rolling blend\B8\A6 \B0\A1\C1\FC.
    //      3) non-tangible edge\B8\A6 \C1\A4\C0\C7\C7ϴ\C2 3 generator\B5\E9\C0\BA \BC\AD\B7\CE rolling blend\B8\A6 \B0\A1\C1\FA\BC\F6 \BE\F8\C0\BD.
    defineProbeTangibleFaces( radiusOfProbe );

    m_hasProbeTangibility = rg_TRUE; 
    m_probeRadius         = radiusOfProbe;
}

void rg_SphereSetVoronoiDiagram::defineProbeTangibleEdges(const rg_REAL& radiusOfProbe)
{
	VDEdge* currEdge = rg_NULL;

    m_GlobalEdgeList.reset4Loop();
	while(m_GlobalEdgeList.setNext4Loop())
	{
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isOnInfinity() )
            currEdge->isTangible( FULLY_TANGIBLE );
    }

	m_GlobalEdgeList.reset4Loop();
	while(m_GlobalEdgeList.setNext4Loop())
	{
		currEdge = m_GlobalEdgeList.getpEntity();

        if (    currEdge->isBoundedEdge() 
             || currEdge->isTangible() != rg_UNKNOWN )
            continue;

        
        rg_dList< pair<VDEdge*, VDVertex*> > edgesToCheckTangibility;

        ///////////////////////////////////////////////////////////////////
        //
        //  initialization of edge-list to check tangibility.

        //  FULLY-TANGIBLE EDGE : min radius >= probe radius 
        if ( rg_LE(radiusOfProbe, currEdge->getMinRadius()) )
        {
            currEdge->isTangible( FULLY_TANGIBLE );

            VDVertex* startVertex = currEdge->getStartVertex();
            startVertex->isTangible(FULLY_TANGIBLE);

            //  define edges to be propagated
            for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
            {
                if ( currEdge != startVertex->getIncidentEdge( i ) )
                {
                    edgesToCheckTangibility.addTail( pair<VDEdge*, VDVertex*>(startVertex->getIncidentEdge( i ), startVertex) );
                }
            }
        }
        //  PARTIALLY-TANGIBLE EDGE : min radius < probe radius <= max radius.
        else if ( rg_LE(radiusOfProbe, currEdge->getMaxRadius()) )
        {
            currEdge->isTangible( PARTIALLY_TANGIBLE );
            currEdge->getStartVertex()->isTangible(NON_TANGIBLE);
        }
        //  NON-TANGIBLE EDGE : probe radius > max radius.
        else 
        {
            currEdge->isTangible( NON_TANGIBLE );
            currEdge->getStartVertex()->isTangible(NON_TANGIBLE);
        }



        ///////////////////////////////////////////////////////////////////
        //
        //  propagation of edge-list to check tangibility.
        pair<VDEdge*, VDVertex*> currEdgeToCheck;
        while ( edgesToCheckTangibility.getSize() > 0 )
        {
            currEdgeToCheck = edgesToCheckTangibility.popFront();

            if ( currEdgeToCheck.first->isTangible() != rg_UNKNOWN )
                continue;

            //  FULLY-TANGIBLE EDGE : min radius >= probe radius 
            if ( rg_LE(radiusOfProbe, currEdgeToCheck.first->getMinRadius()) )
            {
                currEdgeToCheck.first->isTangible( FULLY_TANGIBLE );

                //  define edges to be propagated
                VDVertex* vertexToBePropagated = rg_NULL;
                if ( currEdgeToCheck.first->getStartVertex() == currEdgeToCheck.second )
                    vertexToBePropagated = currEdgeToCheck.first->getEndVertex();
                else
                    vertexToBePropagated = currEdgeToCheck.first->getStartVertex();

                vertexToBePropagated->isTangible(FULLY_TANGIBLE);


                if ( !vertexToBePropagated->isOnInfinity() )
                {
                    for ( rg_INT i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
                    {
                        if ( currEdgeToCheck.first != vertexToBePropagated->getIncidentEdge( i ) )
                        {
                            edgesToCheckTangibility.addTail( pair<VDEdge*, VDVertex*>(vertexToBePropagated->getIncidentEdge( i ), vertexToBePropagated) );
                        }
                    }
                }

            }
            //  PARTIALLY-TANGIBLE EDGE : min radius < probe radius <= max radius.
            else if ( rg_LE(radiusOfProbe, currEdgeToCheck.first->getMaxRadius()) )
            {
                currEdgeToCheck.first->isTangible( PARTIALLY_TANGIBLE );
            }
            //  NON-TANGIBLE EDGE : probe radius > max radius.
            else 
            {
                currEdgeToCheck.first->isTangible( NON_TANGIBLE );

                if ( currEdgeToCheck.first->getStartVertex() == currEdgeToCheck.second )
                    currEdgeToCheck.first->getEndVertex()->isTangible(NON_TANGIBLE);
                else
                    currEdgeToCheck.first->getStartVertex()->isTangible(NON_TANGIBLE);
            }
        }        
    }

    m_GlobalEdgeList.reset4Loop();
	while(m_GlobalEdgeList.setNext4Loop())
	{
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isTangible() == rg_UNKNOWN )
            currEdge->isTangible( NON_TANGIBLE );
    }

    VDVertex* currVertex = rg_NULL;
    m_GlobalVertexList.reset4Loop();
	while(m_GlobalVertexList.setNext4Loop())
	{
		currVertex = m_GlobalVertexList.getpEntity();

        if ( currVertex->isTangible() == rg_UNKNOWN )
            currVertex->isTangible( NON_TANGIBLE );
    }

}

void rg_SphereSetVoronoiDiagram::defineProbeTangibleFaces(const rg_REAL& radiusOfProbe)
{
	VDEdge* currEdge = rg_NULL;

    VDPartialEdge* startPartEdge = rg_NULL;
    VDPartialEdge* currPartEdge  = rg_NULL;

    VDFace*        currFace  = rg_NULL;
    BallGenerator* rightBall = rg_NULL;
    BallGenerator* leftBall  = rg_NULL;

	m_GlobalEdgeList.reset4Loop();
	while(m_GlobalEdgeList.setNext4Loop())
	{
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isTangible() == FULLY_TANGIBLE )
        {
            startPartEdge = currEdge->getPartialEdge();
            currPartEdge  = startPartEdge;
            do
            {
                currFace = currPartEdge->getLoop()->getFace();
                if ( currFace->isOnInfinity() )
                {
                    currFace->isTangible( FULLY_TANGIBLE );

                    if ( currFace->getRightCell()->getID() != -1 )  {
                        rightBall = (BallGenerator*) (currFace->getRightCell()->getGenerator());
                        if ( rightBall->isTangible() != PARTIALLY_TANGIBLE )
                            rightBall->isTangible( FULLY_TANGIBLE );
                    }

                    if ( currFace->getLeftCell()->getID() != -1 )  {
                        leftBall = (BallGenerator*) (currFace->getLeftCell()->getGenerator());
                        if ( leftBall->isTangible() != PARTIALLY_TANGIBLE )
                            leftBall->isTangible( FULLY_TANGIBLE );                        
                    }
                }
                else
                {
                    rightBall = (BallGenerator*) (currFace->getRightCell()->getGenerator());
                    leftBall  = (BallGenerator*) (currFace->getLeftCell()->getGenerator());

                    Sphere rSphere = rightBall->getBall();
                    Sphere lSphere = leftBall->getBall();

		            rg_REAL dist = lSphere.getCenter().distance(rSphere.getCenter()) 
                                   - lSphere.getRadius() - rSphere.getRadius();


		            if( rg_LT( dist, radiusOfProbe*2.0) )
                    {
                        currFace->isTangible( PARTIALLY_TANGIBLE );

                        rightBall->isTangible( PARTIALLY_TANGIBLE );
                        leftBall->isTangible( PARTIALLY_TANGIBLE );
                    }
                    else
                    {
                        currFace->isTangible( FULLY_TANGIBLE );

                        if ( rightBall->isTangible() != PARTIALLY_TANGIBLE )
                            rightBall->isTangible( FULLY_TANGIBLE );

                        if ( leftBall->isTangible() != PARTIALLY_TANGIBLE )
                            leftBall->isTangible( FULLY_TANGIBLE );
                    }
                }

                currPartEdge = currPartEdge->getNextPartialEdgeInRadialCycle();

            } while ( startPartEdge != currPartEdge );
        }
        else if ( currEdge->isTangible() == PARTIALLY_TANGIBLE )
        {
            startPartEdge = currEdge->getPartialEdge();
            currPartEdge  = startPartEdge;
            do
            {
                currFace = currPartEdge->getLoop()->getFace();

                rightBall = (BallGenerator*) (currFace->getRightCell()->getGenerator());
                leftBall  = (BallGenerator*) (currFace->getLeftCell()->getGenerator());

                currFace->isTangible( PARTIALLY_TANGIBLE );

                rightBall->isTangible( PARTIALLY_TANGIBLE );
                leftBall->isTangible( PARTIALLY_TANGIBLE );

                currPartEdge = currPartEdge->getNextPartialEdgeInRadialCycle();

            } while ( startPartEdge != currPartEdge );
        }
        else
        {
        }
    }

    m_GlobalFaceList.reset4Loop();
	while(m_GlobalFaceList.setNext4Loop())
	{
		currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isTangible() == rg_UNKNOWN )
            currFace->isTangible( NON_TANGIBLE );
    }

    BallGenerator* currBall = rg_NULL;
    m_GlobalGeneratorList.reset4Loop();
	while(m_GlobalGeneratorList.setNext4Loop())
	{
		currBall = m_GlobalGeneratorList.getpEntity();

        if ( currBall->isTangible() == rg_UNKNOWN )
            currBall->isTangible( NON_TANGIBLE );
    }

}

void rg_SphereSetVoronoiDiagram::resetTangibility()
{
    m_hasProbeTangibility = rg_FALSE; 
    m_probeRadius         = 0.0;

    VDFace* currFace = rg_NULL;
    m_GlobalFaceList.reset4Loop();
	while(m_GlobalFaceList.setNext4Loop())
	{
		currFace = m_GlobalFaceList.getpEntity();
        
        currFace->isTangible( rg_UNKNOWN );
    }

    BallGenerator* currBall = rg_NULL;
    m_GlobalGeneratorList.reset4Loop();
	while(m_GlobalGeneratorList.setNext4Loop())
	{
		currBall = m_GlobalGeneratorList.getpEntity();
            
        currBall->isTangible( rg_UNKNOWN );
    }

    VDEdge* currEdge = rg_NULL;
    m_GlobalEdgeList.reset4Loop();
	while(m_GlobalEdgeList.setNext4Loop())
	{
		currEdge = m_GlobalEdgeList.getpEntity();
            
        currEdge->isTangible( rg_UNKNOWN );
    }

    VDVertex* currVertex = rg_NULL;
    m_GlobalVertexList.reset4Loop();
	while(m_GlobalVertexList.setNext4Loop())
	{
		currVertex = m_GlobalVertexList.getpEntity();
            
        currVertex->isTangible( rg_UNKNOWN );
    }

}


rg_dList< pair<VDEdge*, rg_FLAG> >* rg_SphereSetVoronoiDiagram::findSolventTrajectory()
{
    rg_REAL    x_dist = rg_ABS(m_maxPointOfBoundingBoxForBalls.getX() - m_minPointOfBoundingBoxForBalls.getX());

    rg_Point3D midPoint = (m_minPointOfBoundingBoxForBalls + m_maxPointOfBoundingBoxForBalls)*0.5;

    rg_Point3D startPoint = midPoint;
    startPoint.setX( startPoint.getX()- (1*x_dist) );

    rg_Point3D endPoint   = midPoint;
    endPoint.setX( endPoint.getX() + (1*x_dist) );

    startVertexOnSolventTrajectory = findClosestVertexToGivenPoint(startPoint);
    endVertexOnSolventTrajectory   = findClosestVertexToGivenPoint(endPoint);

    VDVertex* startVertex = startVertexOnSolventTrajectory;
    VDVertex* endVertex   = endVertexOnSolventTrajectory;

    rg_Point3D targetPoint = endVertex->getPoint();

    rg_Point3D idealDirection = targetPoint - startVertex->getPoint();

    rg_dList< PathInvestigator > solventTrajectory;
    rg_dList< VDVertex* >        verticesOnTrajectory;

    VDVertex* propagatingVertex = startVertex;
    verticesOnTrajectory.addTail( propagatingVertex );

    ///////////////////////////////////////////////////////////////////////////
    //
    // initialization for solvent trajectory.
    PathInvestigator onePath;
    onePath.setPropagatingVertex( propagatingVertex );


    rg_REAL angleOfEdges[4];
    rg_FLAG orientationOfEdges[4];
    rg_FLAG orderOfEdge[4];
    int i=0;
	for ( i=0; i<4; i++)
    {
        VDEdge* currEdge = propagatingVertex->getIncidentEdge(i);

        if ( currEdge->isBoundedEdge() && currEdge->isTangible() == FULLY_TANGIBLE )
        {
            if ( currEdge->getStartVertex() == propagatingVertex )
            {
                rg_Point3D edgeDirection = currEdge->getEndVertex()->getPoint() - currEdge->getStartVertex()->getPoint();
                angleOfEdges[i]       = idealDirection.angle(edgeDirection);
                orientationOfEdges[i] = rg_TRUE;
            }
            else
            {
                rg_Point3D edgeDirection = currEdge->getStartVertex()->getPoint() - currEdge->getEndVertex()->getPoint();
                angleOfEdges[i] = idealDirection.angle(edgeDirection);
                orientationOfEdges[i] = rg_FALSE;
            }
        }
        else
        {
            angleOfEdges[i] = DBL_MAX;
            orientationOfEdges[i] = -1;
        }
    }

    for ( i=0; i<4; i++ )
    {
        int order = 0;
        int j=0;
		for ( j=0; j<4; j++ )
        {
            if ( i==j )
                continue;

            if ( rg_LT( angleOfEdges[i], angleOfEdges[j] ) )
                continue;
            else
                order++;
        }
        orderOfEdge[i] = order;
    }

    for ( i=0; i<4; i++ )
    {
        if ( orderOfEdge[i] == 3 )
            continue;

        if ( orientationOfEdges[i] == -1 )
            continue;

        onePath.setPath(orderOfEdge[i], propagatingVertex->getIncidentEdge(i), orientationOfEdges[i]);
    }    
    onePath.setCurrPath(0);

    solventTrajectory.addTail( onePath );

    PathInvestigator* currInvestigator = solventTrajectory.getLastpEntity();

    ///////////////////////////////////////////////////////////////////////////
    //
    // main algorithm for solvent trajectory.    
    while ( currInvestigator->getPropagatingVertex() != endVertex )//&& solventTrajectory.getSize() <40 )
    {
        if ( currInvestigator->isValidInvestigator() == rg_FALSE )
        {
            solventTrajectory.killTail();
            verticesOnTrajectory.killTail();

            currInvestigator = solventTrajectory.getLastpEntity();
            currInvestigator->killCurrPath();

            continue;
        }

        PathInvestigator newInvestigator;

        if ( currInvestigator->makeNextInvestigator(targetPoint, verticesOnTrajectory, newInvestigator) == rg_TRUE )
        {
            solventTrajectory.addTail( newInvestigator );
            currInvestigator = solventTrajectory.getLastpEntity();
        }
        else
        {
            currInvestigator->killCurrPath();
        }
    }

    rg_dList< pair<VDEdge*, rg_FLAG> >* edgesOnSolventTrajectory = new rg_dList< pair<VDEdge*, rg_FLAG> >;


	solventTrajectory.reset4Loop();
	while(solventTrajectory.setNext4Loop())
	{
		currInvestigator = solventTrajectory.getpEntity();

        pair<VDEdge*, rg_FLAG> edgeToAdd;

        edgeToAdd.first  = currInvestigator->getCurrPath();
        edgeToAdd.second = currInvestigator->getOrientationOfCurrPath();

        edgesOnSolventTrajectory->addTail(edgeToAdd);
    }
    edgesOnSolventTrajectory->killTail();

    return edgesOnSolventTrajectory;
    
/*
    rg_dList< pair<VDEdge*, rg_FLAG> >* edgesOnSolventTrajectory = new rg_dList< pair<VDEdge*, rg_FLAG> >;

    VDEdge* currEdge = rg_NULL;

    VDEdge* prevPropagatingEdge = rg_NULL;

    const rg_INT maxNumOfEdges = 5;

    rg_INT countOfEdges  = 0;

    rg_FLAG bDirectionOfPrevEdge = rg_FALSE;
    rg_FLAG bSuccess             = rg_FALSE;

    VDVertex* propagatingVertex = rg_NULL;

	m_GlobalEdgeList.reset4Loop();
	while(m_GlobalEdgeList.setNext4Loop())
	{
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( !currEdge->isBoundedEdge() && currEdge->isTangible() == FULLY_TANGIBLE )
        {
            propagatingVertex   = currEdge->getStartVertex();
            prevPropagatingEdge = currEdge;

            bDirectionOfPrevEdge = rg_FALSE;
            bSuccess = rg_FALSE;
            do
            {
                int numFailure = 0;
                for ( int j=0; j<4; j++)
                {
                    VDEdge* tempEdge = propagatingVertex->getIncidentEdge(j);
                    if ( tempEdge != prevPropagatingEdge && tempEdge->isTangible() == FULLY_TANGIBLE )
                    {
                        //propagatingEdges[ countOfEdges++ ] = tempEdge;
                        pair<VDEdge*, rg_FLAG> edgeToAdd;

                        edgeToAdd.first = tempEdge;

                        prevPropagatingEdge = tempEdge;
                        if ( tempEdge->getStartVertex() == propagatingVertex )
                        {
                            propagatingVertex = tempEdge->getEndVertex();
                            edgeToAdd.second = rg_TRUE;
                        }
                        else
                        {
                            propagatingVertex = tempEdge->getStartVertex();
                            edgeToAdd.second = rg_FALSE;
                        }

                        edgesOnSolventTrajectory->addTail( edgeToAdd );

                        break;
                    }
                    else
                    {
                        numFailure++;
                    }
                }
                
                if ( numFailure < 4 )
                {
                    bSuccess = rg_TRUE;
                }
                else
                {
                    bSuccess = rg_FALSE;
                    break;
                }
            } while ( edgesOnSolventTrajectory->getSize() < maxNumOfEdges );   
            
            if ( bSuccess == rg_TRUE )
                break;
            else
            {
                edgesOnSolventTrajectory->removeAll();
                
                continue;
            }
        }
    }

//    for (i=0; i<maxNumOfEdges; i++)
//        edgesOnSolventTrajectory->addTail( propagatingEdges[i]);


    return edgesOnSolventTrajectory;

*/
}

///////////////////////////////////////////////////////////////////////////////
//
//  check Voronoi diagram
rg_FLAG rg_SphereSetVoronoiDiagram::isEmptySphere(const Sphere& aSphere, const rg_REAL& res)
{   
    VDCell* currCell = rg_NULL;

    rg_Point3D center = aSphere.getCenter();
    rg_REAL    radius = aSphere.getRadius();

	m_GlobalCellList.reset4Loop();
	while( m_GlobalCellList.setNext4Loop() )  
    {
        currCell = m_GlobalCellList.getpEntity();

        if ( currCell->getID() == -1 )
            continue;

        Sphere ball = ((BallGenerator*) currCell->getGenerator())->getBall();

        rg_REAL distanceBetweenCenters = center.distance( ball.getCenter() );
        rg_REAL sumOfRadii             = radius + ball.getRadius();
        rg_REAL intersection           = distanceBetweenCenters - sumOfRadii;

        if( !rg_ZERO( intersection, res ) )
            if ( !rg_POS( intersection ) )
                return rg_FALSE;
    }

    return rg_TRUE;
}

rg_FLAG rg_SphereSetVoronoiDiagram::isValid()
{
//	m_bucket.constructBucket( 0.5,
//                              m_minPointOfBoundingBoxForBalls,
//                              m_maxPointOfBoundingBoxForBalls,
//                              m_GlobalGeneratorList );

    rg_dList<BallGenerator*>* neighborBalls = rg_NULL;
    BallGenerator*            currBall      = rg_NULL;

    VDVertex* currVertex = rg_NULL;
	m_GlobalVertexList.reset4Loop();
	while( m_GlobalVertexList.setNext4Loop() )  
    {
		currVertex = m_GlobalVertexList.getpEntity();

        Sphere tangentSphere(currVertex->getPoint(), currVertex->getRadiusOfTangentSphere());
        neighborBalls = m_bucket.getCellsByMask( tangentSphere );
        
        neighborBalls->reset4Loop();
	    while( neighborBalls->setNext4Loop() )  
        {
            currBall = neighborBalls->getEntity();

            rg_REAL distanceBetweenCenters = tangentSphere.getCenter().distance( currBall->getCenter() );
            rg_REAL sumOfRadii             = tangentSphere.getRadius() + currBall->getRadius();
            rg_REAL intersection           = distanceBetweenCenters - sumOfRadii;

            if( !rg_ZERO( intersection, resNeg5 ) )  {
                if ( !rg_POS( intersection ) )  {
//                    m_bucket.removeBucketTable();
                    delete neighborBalls;
                    return rg_FALSE;
                }
            }
        }
        delete neighborBalls;
    }

    /*
    VDEdge* currEdge = rg_NULL;
	m_GlobalEdgeList.reset4Loop();
	while( m_GlobalEdgeList.setNext4Loop() )  
    {
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isOnInfinity() )
            continue;

    	rg_RBzCurve3D* curve       = (rg_RBzCurve3D*) currEdge->getEdgeEquation(); 
        rg_Point3D     pointOnEdge = curve->evaluatePt( 0.5 );

        VDCell* oneCell = currEdge->getPartialEdge()->getLoop()->getFace()->getRightCell();
        Sphere  oneBall = ((BallGenerator*) oneCell->getGenerator())->getBall();

        rg_REAL diatanceBetweenPointAndBall = pointOnEdge.distance( oneBall.getCenter() ) - oneBall.getRadius();

        if ( !isEmptySphere( Sphere(pointOnEdge, diatanceBetweenPointAndBall), resNeg3 ) )
            return rg_FALSE;
    }
    */
    //m_bucket.removeBucketTable();

    return rg_TRUE;
}



void rg_SphereSetVoronoiDiagram::perturbBallGenerators()
{
    BallGenerator* currBall = rg_NULL;

    rg_INT  numBalls = m_GlobalGeneratorList.getSize();
    rg_REAL  lowerBoundOfPerturbation = log10((double)numBalls)+3.0;

    rg_REAL unitOfPerturbation = pow(10.0, -lowerBoundOfPerturbation);

    m_GlobalGeneratorList.reset4Loop();
    while ( m_GlobalGeneratorList.setNext4Loop() )  {
        currBall = m_GlobalGeneratorList.getpEntity();

        rg_Point3D center = currBall->getCenter();
        rg_REAL x = center.getX() + 1.0*currBall->getCell()->getID()*unitOfPerturbation;
        rg_REAL y = center.getY() + 2.0*currBall->getCell()->getID()*unitOfPerturbation;
        rg_REAL z = center.getZ() + 3.0*currBall->getCell()->getID()*unitOfPerturbation;

        currBall->setCenter( rg_Point3D(x, y, z) );
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//  file I/O
/*
void rg_SphereSetVoronoiDiagram::readGenerators(const char* filename)
{
    rg_REAL minX =  DBL_MAX;
    rg_REAL minY =  DBL_MAX;
    rg_REAL minZ =  DBL_MAX;

    rg_REAL maxX = -DBL_MAX;
    rg_REAL maxY = -DBL_MAX;
    rg_REAL maxZ = -DBL_MAX;


    rg_REAL sumOfRadius = 0.0;

    
    std::ifstream fin(filename);
    
    rg_INT numGenerator;
    fin >> numGenerator;

    int i=0;
	for (i=0; i<numGenerator; i++)  {
        rg_REAL x;
        rg_REAL y;
        rg_REAL z;
        rg_REAL radius;
        char    atomType;
        
        fin >> x;
        fin >> y;
        fin >> z;
        fin >> radius;
        fin >> atomType;
        
        Attribute attribute;

		if(atomType == 'C')  
			attribute.setColor3f(C_COLOR);
		else if(atomType == 'H') 
			attribute.setColor3f(H_COLOR);
		else if(atomType == 'N')  
			attribute.setColor3f(N_COLOR);
		else if(atomType == 'O')  
			attribute.setColor3f(O_COLOR);
		else if(atomType == 'S') 
			attribute.setColor3f(S_COLOR);
		else if(atomType == 'P') 
			attribute.setColor3f(P_COLOR);
		else{}


        BallGenerator aBall(i, x, y, z, radius, &attribute);

        VDCell        aCell(i);

        BallGenerator* ptrABall = m_GlobalGeneratorList.addTail( aBall );
        VDCell*        ptrACell = m_GlobalCellList.addTail( aCell );

        connectCellAndGenerator( ptrACell, ptrABall );


		if( rg_GT( x + radius, maxX ) )
		    maxX = x + radius;
		if( rg_GT( y + radius, maxY ) )
		    maxY = y + radius;
		if( rg_GT( z + radius, maxZ ) )
		    maxZ = z + radius;

		if( rg_LT( x - radius, minX ) )
		    minX = x - radius;
		if( rg_LT( y - radius, minY ) )
		    minY = y - radius;
		if( rg_LT( z - radius, minZ ) )
		    minZ = z - radius;


        sumOfRadius += radius;
    }

    m_AvgRadiusOfBalls = sumOfRadius/numGenerator;

    m_minPointOfBoundingBoxForBalls.setPoint( minX, minY, minZ );
    m_maxPointOfBoundingBoxForBalls.setPoint( maxX, maxY, maxZ );
}



rg_INT rg_SphereSetVoronoiDiagram::readGeneratorsInPDB(const char* filename)
{
    rg_REAL minX =  DBL_MAX;
    rg_REAL minY =  DBL_MAX;
    rg_REAL minZ =  DBL_MAX;

    rg_REAL maxX = -DBL_MAX;
    rg_REAL maxY = -DBL_MAX;
    rg_REAL maxZ = -DBL_MAX;

    rg_REAL sumOfRadius = 0.0;


    std::ifstream fin;
    fin.open(filename);

	char	recordline[85];
	char*	seps = " \b,\t\n\v\r["; //\BF\A9\B9\E9 üũ
	char*	token = NULL;
	char	atomType;
	char	chainID;		//[\C1\A4\B9\CE\C3߰\A1] : to identify the group of chain

    rg_REAL x;
	rg_REAL y;
    rg_REAL z;
    rg_REAL radius;

    rg_INT  idOfAtom;

	int i=0;
	int posOfCoord;
    char atomSerial[6];
	char X_coordValueInCharType[8];
	char Y_coordValueInCharType[8];
	char Z_coordValueInCharType[8];
	int count = 0;

    char   charResSeq[5];
    rg_INT resSeq = -1;
    char   iCode  = ' ';



	while(fin.getline(recordline, 85))
	{
		token = strtok(recordline, seps); // space\B0\A1 \B3\AA\BFö\A7\B1\EE\C1\F6 token\BF\A1 \C0\FA\C0\E5
		if( !strcmp(token, "ATOM") )
		{
            resSeq = -1;
            iCode  = -1;
			i=0;
			for(posOfCoord=6; posOfCoord<=10; posOfCoord++)  {
				atomSerial[i++] = recordline[posOfCoord];
			}
			idOfAtom = atoi(atomSerial);

            
            atomType = recordline[13];
			chainID  = recordline[21];	//[\C1\A4\B9\CE\C3߰\A1] : to identify the group of chain
            
			i=0;
			for(posOfCoord=22; posOfCoord<=25; posOfCoord++)  {
				charResSeq[i++] = recordline[posOfCoord];
			}
			resSeq = atoi(charResSeq);
            
            if ( recordline[26] != ' ' )
            {
                iCode = recordline[26];
            }


            PDBAtomAttribute anPDAAtom;

            anPDAAtom.setAtomSerialNumber(idOfAtom);
            anPDAAtom.setResidueSequenceNumber(resSeq);
            anPDAAtom.setCodeForInsertionOfResidues(iCode);

			if(atomType == 'C')  {
				anPDAAtom.setPDBData(C_COLOR, "C", chainID);
                radius = C_RADIUS;
			}
			else if(atomType == 'H')  {
				anPDAAtom.setPDBData(H_COLOR, "H", chainID);
				radius = H_RADIUS;
			}
			else if(atomType == 'N')  {
				anPDAAtom.setPDBData(N_COLOR, "N", chainID);
				radius = N_RADIUS;
			}
			else if(atomType == 'O')  {
				anPDAAtom.setPDBData(O_COLOR, "O", chainID);
				radius = O_RADIUS;
			}
			else if(atomType == 'S')  {
				anPDAAtom.setPDBData(S_COLOR, "S", chainID);
				radius = S_RADIUS;
			}
			else if(atomType == 'P')  {
				anPDAAtom.setPDBData(P_COLOR, "P", chainID);
				radius = P_RADIUS;
			}
			else{}


			i=0;
			for(posOfCoord=30; posOfCoord<=37; posOfCoord++)  {
				X_coordValueInCharType[i++] = recordline[posOfCoord];
			}
			x = atof(X_coordValueInCharType);


			i=0;
			for(posOfCoord=38; posOfCoord<=45; posOfCoord++)  {
				Y_coordValueInCharType[i++] = recordline[posOfCoord];
			}
			y = atof(Y_coordValueInCharType);
			
			i=0;
			for(posOfCoord=46; posOfCoord<=53; posOfCoord++)  {
				Z_coordValueInCharType[i++] = recordline[posOfCoord];
			}
			z = atof(Z_coordValueInCharType);


            BallGenerator aBall( m_GlobalGeneratorList.getSize(), x, y, z, radius, (Attribute*) &anPDAAtom );
            VDCell        aCell( m_GlobalCellList.getSize() );

            BallGenerator* ptrABall = m_GlobalGeneratorList.addTail( aBall );
            VDCell*        ptrACell = m_GlobalCellList.addTail( aCell );

            connectCellAndGenerator( ptrACell, ptrABall );

			count++;



		    if( rg_GT( x + radius, maxX ) )
		        maxX = x + radius;
		    if( rg_GT( y + radius, maxY ) )
		        maxY = y + radius;
		    if( rg_GT( z + radius, maxZ ) )
		        maxZ = z + radius;

		    if( rg_LT( x - radius, minX ) )
		        minX = x - radius;
		    if( rg_LT( y - radius, minY ) )
		        minY = y - radius;
		    if( rg_LT( z - radius, minZ ) )
		        minZ = z - radius;

            sumOfRadius += radius;

		}
    }

    m_AvgRadiusOfBalls = sumOfRadius/count;

    m_minPointOfBoundingBoxForBalls.setPoint( minX, minY, minZ );
    m_maxPointOfBoundingBoxForBalls.setPoint( maxX, maxY, maxZ );

    fin.close();

    return rg_TRUE;
}



void rg_SphereSetVoronoiDiagram::fileOutDelaunayEdgeLength( ofstream& fout )
{
    VDFace* currFace = rg_NULL;


    fout << "//  0(Serial) 1(Serial) : 0(resSeq) 1(resSeq) : 0((ChainID) 1(ChainID) : 0(iCode) 1(iCode) : dist(0,1)" << endl << endl;

	m_GlobalFaceList.reset4Loop();
	while( m_GlobalFaceList.setNext4Loop() )  
    {
		currFace = m_GlobalFaceList.getpEntity();

        BallGenerator* rightGenerator = (BallGenerator*) currFace->getRightCell()->getGenerator();
        BallGenerator* leftGenerator  = (BallGenerator*) currFace->getLeftCell()->getGenerator();

        rg_Point3D rightBall = rightGenerator->getCenter();
        rg_Point3D leftBall  = leftGenerator->getCenter();

        ///////////////////////////////////////////////////////////
        //  modified by Y.Cho  2004-10-20
//        PDBAtomAttribute* rightPDBAttri = (PDBAtomAttribute*) ((BallGenerator*) (currFace->getRightCell()->getGenerator()))->getAttribute();
//        PDBAtomAttribute* leftPDBAttri  = (PDBAtomAttribute*) ((BallGenerator*) (currFace->getLeftCell()->getGenerator()))->getAttribute();
//
//        fout << rightPDBAttri->getAtomSerialNumber() << "\t" << leftPDBAttri->getAtomSerialNumber() << "\t:\t";
//        fout << rightPDBAttri->getResidueSequenceNumber() << "\t" << leftPDBAttri->getResidueSequenceNumber() << "\t:\t";
//
//        if ( rightPDBAttri->getGroupID() != ' ' )
//            fout << rightPDBAttri->getGroupID();
//        fout << "\t";
//        if ( leftPDBAttri->getGroupID() != ' ' )
//            fout << leftPDBAttri->getGroupID();
//        fout << "\t:\t";
//
//        if ( rightPDBAttri->getCodeForInsertionOfResidues() != -1 )
//            fout << rightPDBAttri->getCodeForInsertionOfResidues();
//        fout << "\t";
//        if ( leftPDBAttri->getCodeForInsertionOfResidues() != -1 )
//            fout << leftPDBAttri->getCodeForInsertionOfResidues();
//        fout << "\t:\t";
//
//        fout << rightBall.distance( leftBall ) << endl;
                
        PDBAtom* rightPDBAtom = (PDBAtom*)(rightGenerator->getProperty());
        PDBAtom* leftPDBAtom  = (PDBAtom*)(leftGenerator->getProperty());

        fout << rightPDBAtom->getAtomSerial() << "\t" << leftPDBAtom->getAtomSerial() << "\t:\t";
        fout << rightPDBAtom->getResidue()->getResidueSequqnceNumber() << "\t" << leftPDBAtom->getResidue()->getResidueSequqnceNumber() << "\t:\t";

        if ( rightPDBAtom->getChainID() != ' ' )
            fout << rightPDBAtom->getChainID();
        fout << "\t";
        if ( leftPDBAtom->getChainID() != ' ' )
            fout << leftPDBAtom->getChainID();
        fout << "\t:\t";

        if ( leftPDBAtom->getCodeForInsertionOfResidues() != -1 )
            fout << rightPDBAtom->getCodeForInsertionOfResidues();
        fout << "\t";
        if ( leftPDBAtom->getCodeForInsertionOfResidues() != -1 )
            fout << rightPDBAtom->getCodeForInsertionOfResidues();
        fout << "\t:\t";

        fout << rightBall.distance( leftBall ) << endl;


    }
}
*/



rg_FLAG rg_SphereSetVoronoiDiagram::reportConstructionEVDS(ofstream& fout)
{
    rg_FLAG bValid = isValid();

    if ( bValid )
        fout << "It is a valid EVD(S)." << endl << endl;
    else
        fout << "NOTE! It is NOT a valid EVD(S)." << endl << endl;

    fout << "# of Ball:   \t" << m_GlobalGeneratorList.getSize() << endl;
    fout << "# of Cell:   \t" << m_GlobalCellList.getSize() << endl;
    fout << "# of Face:   \t" << m_GlobalFaceList.getSize() << endl;
    fout << "# of Loop:   \t" << m_GlobalLoopList.getSize() << endl;
    fout << "# of Edge:   \t" << m_GlobalEdgeList.getSize() << endl;
    fout << "# of Vertex: \t" << m_GlobalVertexList.getSize() << endl;
    fout << endl;  

    /*
    fout << endl << "===  Time Analysis  ===" << endl;
    rg_REAL times[NUM_TIMES_FOR_ANALYSIS];
    for (rg_INT i=0; i<NUM_TIMES_FOR_ANALYSIS; i++)
        times[i] = (double)(computingTime[i])/CLOCKS_PER_SEC;
 

    fout << "   * Total time  : \t" << times[16] << endl;
    */


    fout << endl << "===  Time Analysis  ===" << endl;

    computingTime[1]  -= (timeToBeExtracted[0]+timeToBeExtracted[1]);
    computingTime[3]  -= timeToBeExtracted[0];
    computingTime[7]  -= timeToBeExtracted[1];
    computingTime[16] -= (timeToBeExtracted[0]+timeToBeExtracted[1]);
    rg_REAL times[NUM_TIMES_FOR_ANALYSIS];
    rg_INT i=0;
	for (i=0; i<NUM_TIMES_FOR_ANALYSIS; i++)
        times[i] = (double)(computingTime[i])/CLOCKS_PER_SEC;
 

    fout << "   * Total time  : \t" << times[16] << endl;
    fout << "   * BIG-WORLD -------------------- " << endl;
    fout << "      1. Initializing EVD(S)   : "     << times[0]  << "\t\t" << 100.*times[0]/times[16]  << " %" << endl;
    fout << "      2. finding end vertex    : "     << times[1]  << "\t\t" << 100.*times[1]/times[16]  << " %" << endl;
    fout << "         1) preparing Ve       :   "   << times[2]  << "  \t" << 100.*times[2]/times[16]  << " %" << endl;    
    fout << "         2) parabolic edge     :   "   << times[3]  << "  \t" << 100.*times[3]/times[16]  << " %" <<  endl;
    fout << "             a) in front FR    :     " << times[4]  << "\t"   << 100.*times[4]/times[16]  << " %" <<  endl;
    fout << "             b) from front FR  :     " << times[5]  << "\t"   << 100.*times[5]/times[16]  << " %" <<  endl;
    fout << "             c) in rear FR     :     " << times[6]  << "\t"   << 100.*times[6]/times[16]  << " %" <<  endl;
    fout << "         3) linear edge        :   "   << times[7]  << "  \t" << 100.*times[7]/times[16]  << " %" <<  endl;
    fout << "             a) in front FR    :     " << times[8]  << "\t"   << 100.*times[8]/times[16]  << " %" <<  endl;
    fout << "             b) from front FR  :     " << times[9]  << "\t"   << 100.*times[9]/times[16]  << " %" <<  endl;
    fout << "             c) in rear FR     :     " << times[10] << "\t"   << 100.*times[10]/times[16] << " %" <<  endl;
    fout << "         4) elliptic edge      :   "   << times[11] << "  \t" << 100.*times[11]/times[16] << " %" <<  endl;
    fout << "      3. exploring edge        : "     << times[12] << "\t\t" << 100.*times[12]/times[16] << " %" <<  endl;
    fout << "      4. local topology        : "     << times[13] << "\t\t" << 100.*times[13]/times[16] << " %" <<  endl;
    fout << "      5. global topology       : "     << times[14] << "\t\t" << 100.*times[14]/times[16] << " %" <<  endl;
    fout << endl;
    fout << "   * SMALL-WORLD -------------------- " << endl;
    fout << "     1. constructing small worlds : " << times[15] << "\t\t" << 100.*times[15]/times[16] << " %" << endl;
    fout << endl << endl;

    if ( m_bucket.isConstructed() )  {    
        fout << endl << "===  Bucket Analysis  ===" << endl;    
        fout << "   * element size                    : \t" << m_sizeOfBucketElement << endl;
    }

    /*
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //
    fout << endl;
    fout << "===  Candidate Balls Analysis  ===" << endl;   

    rg_INT  numEdges[3][2]   = { {0, 0}, {0, 0}, {0, 0} };
    rg_INT  numInfiniteEdges = 0;
            //  numEdges[0]     : parabolic edge
            //  numEdges[1]     : linear edge
            //  numEdges[2]     : elliptic edge
            //  numEdges[i][0]     : bounded edge
            //  numEdges[i][1]     : unbounded edge
    
    VDEdge* currEdge = rg_NULL;
    m_GlobalEdgeList.reset4Loop();
    while ( m_GlobalEdgeList.setNext4Loop() )  {
        currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isOnInfinity() )
            numInfiniteEdges++;
        else  {        
            if ( currEdge->getEdgeType() == PARABOLIC_OR_HYPERBOLIC_EDGE )  {
                if ( currEdge->isBoundedEdge() )
                    numEdges[0][0]++;
                else 
                    numEdges[0][1]++;
            }
            else if ( currEdge->getEdgeType() == LINEAR_EDGE )  {
                if ( currEdge->isBoundedEdge() )
                    numEdges[1][0]++;
                else 
                    numEdges[1][1]++;
            }
            else if ( currEdge->getEdgeType() == CIRCULAR_OR_ELLIPTIC_EDGE )  {
                if ( currEdge->isBoundedEdge() )
                    numEdges[2][0]++;
                else 
                    numEdges[2][1]++;
            }
            else  {
            }
        }
    }


    rg_INT sumNumEdge[3];
    sumNumEdge[0] = numEdges[0][0]+numEdges[1][0]+numEdges[2][0];
    sumNumEdge[1] = numEdges[0][0]+numEdges[1][0]+numEdges[2][0]+numEdges[0][1]+numEdges[1][1]+numEdges[2][1];
    sumNumEdge[2] = numEdges[0][0]+numEdges[1][0]+numEdges[2][0]+numEdges[0][1]+numEdges[1][1]+numEdges[2][1]+numInfiniteEdges;

    fout << "   * num. of edges: " << endl;
    fout << "                               \t" << "bounded" << "\t\t\t\t\t" << "unbounded" << "\t\t" << "sum" << endl;
    fout << "       1. non-elliptic edges : \t" << numEdges[0][0]+numEdges[1][0] << "\t" 
                                                << (100.0*(numEdges[0][0]+numEdges[1][0]))/sumNumEdge[0] << "\t"
                                                << (100.0*(numEdges[0][0]+numEdges[1][0]))/sumNumEdge[1] << "\t"
                                                << (100.0*(numEdges[0][0]+numEdges[1][0]))/sumNumEdge[2] << "\t"
                                                << numEdges[0][1]+numEdges[1][1] << "\t" 
                                                << numEdges[0][0]+numEdges[1][0]+numEdges[0][1]+numEdges[1][1] << endl;
    fout << "          1) parabolic edges : \t" << numEdges[0][0] << "\t" 
                                                << (100.0*(numEdges[0][0]))/sumNumEdge[0] << "\t" 
                                                << (100.0*(numEdges[0][0]))/sumNumEdge[1] << "\t" 
                                                << (100.0*(numEdges[0][0]))/sumNumEdge[2] << "\t" 
                                                << numEdges[0][1] << "\t" << numEdges[0][0]+numEdges[0][1] << endl;
    fout << "          2) linear edges    : \t" << numEdges[1][0] << "\t" 
                                                << (100.0*(numEdges[1][0]))/sumNumEdge[0] << "\t" 
                                                << (100.0*(numEdges[1][0]))/sumNumEdge[1] << "\t" 
                                                << (100.0*(numEdges[1][0]))/sumNumEdge[2] << "\t" 
                                                << numEdges[1][1] << "\t" << numEdges[1][0]+numEdges[1][1] << endl;
    fout << "       2. elliptic edges     : \t" << numEdges[2][0] << "\t" 
                                                << (100.0*(numEdges[2][0]))/sumNumEdge[0] << "\t" 
                                                << (100.0*(numEdges[2][0]))/sumNumEdge[1] << "\t" 
                                                << (100.0*(numEdges[2][0]))/sumNumEdge[2] << "\t" 
                                                << numEdges[2][1] << "\t" << numEdges[2][0]+numEdges[2][1] << endl;
    fout << "       3. sum                  \t" << sumNumEdge[0]  << "\t\t\t\t\t" 
                                                << sumNumEdge[1]-sumNumEdge[0] << "\t" 
                                                << sumNumEdge[1] << "\t" 
                                                << "(" << numInfiniteEdges << ")" << "\t" << m_GlobalEdgeList.getSize() << endl;
    fout << endl;


    rg_INT numTotalEdges = m_GlobalEdgeList.getSize();
    rg_INT numTotalBalls = m_GlobalGeneratorList.getSize();
    fout << "                   \t" << "total # balls" << "\t" << "avg. # balls" << "\t" << "# edges" << endl;
    fout << "   * full search : \t" << (numTotalBalls-3)*numTotalEdges << "\t" << (numTotalBalls-3) << endl;
    fout << endl;
    fout << "   * feasible search : \t" << endl;
    rg_INT numBallsInFR[3];
    numBallsInFR[0] =   m_numCandidateBallsOfNonElliptic[0][0][0]+m_numCandidateBallsOfNonElliptic[0][0][1]
                      + m_numCandidateBallsOfNonElliptic[0][1][0]+m_numCandidateBallsOfNonElliptic[0][1][1]
                      + m_numCandidateBallsOfNonElliptic[0][0][2]+m_numCandidateBallsOfNonElliptic[0][1][2];
    numBallsInFR[1] =   m_numCandidateBallsOfNonElliptic[0][0][0]+m_numCandidateBallsOfNonElliptic[0][0][1]
                      + m_numCandidateBallsOfNonElliptic[0][0][2];
    numBallsInFR[2] =   m_numCandidateBallsOfNonElliptic[0][1][0]+m_numCandidateBallsOfNonElliptic[0][1][1]
                      + m_numCandidateBallsOfNonElliptic[0][1][2];
    rg_INT numEndBallsInFR[3][2];
    numEndBallsInFR[0][0] =   m_numEndBallsInPosition[0][0]+m_numEndBallsInPosition[1][0]
                            + m_numEndBallsInPosition[0][1]+m_numEndBallsInPosition[1][1];
    
    numEndBallsInFR[0][1] =   m_numEndBallsInPosition[0][2]+m_numEndBallsInPosition[1][2];
    numEndBallsInFR[1][0] =   m_numEndBallsInPosition[0][0]+m_numEndBallsInPosition[0][1];
    numEndBallsInFR[1][1] =   m_numEndBallsInPosition[0][2];
    numEndBallsInFR[2][0] =   m_numEndBallsInPosition[1][0]+m_numEndBallsInPosition[1][1];
    numEndBallsInFR[2][1] =   m_numEndBallsInPosition[1][2];

    rg_INT numBallsInFFR[3];
    numBallsInFFR[0] = m_numIntersectingBallsWithFrontFR[0]+m_numIntersectingBallsWithFrontFR[1];
    numBallsInFFR[1] = m_numIntersectingBallsWithFrontFR[0];
    numBallsInFFR[2] = m_numIntersectingBallsWithFrontFR[1];
    fout << "                            : \t" << "balls"         << "\t" << "balls in FR  "                                            << "\t" << "A) end balls"        << "\t" << "B) end balls"         << "\t" << "C) balls "      << "\t" << "D) C/A"      << endl;
    fout << "                            : \t" << "in FR"         << "\t" << "/bounded edge"                                            << "\t" << " in front FR"        << "\t" << "  in rear FR"         << "\t" << "   for A "      << "\t" << " balls/edge" << endl;
    fout << "       1. non-elliptic edge : \t" << numBallsInFR[0] << "\t" << ((rg_REAL)numBallsInFR[0])/(numEdges[0][0]+numEdges[1][0]) << "\t" << numEndBallsInFR[0][0] << "\t" << numEndBallsInFR[0][1]  << "\t" << numBallsInFFR[0] << "\t" << ((rg_REAL)numBallsInFFR[0])/numEndBallsInFR[0][0] << endl;
    fout << "          1) parabolic edge : \t" << numBallsInFR[1] << "\t" << ((rg_REAL)numBallsInFR[0])/numEdges[0][0]                  << "\t" << numEndBallsInFR[1][0] << "\t" << numEndBallsInFR[1][1]  << "\t" << numBallsInFFR[1] << "\t" << ((rg_REAL)numBallsInFFR[1])/numEndBallsInFR[1][0] << endl;
    fout << "          2) linear    edge : \t" << numBallsInFR[2] << "\t" << ((rg_REAL)numBallsInFR[0])/numEdges[1][0]                  << "\t" << numEndBallsInFR[2][0] << "\t" << numEndBallsInFR[2][1]  << "\t" << numBallsInFFR[2] << "\t" << ((rg_REAL)numBallsInFFR[2])/numEndBallsInFR[2][0] << endl;
    fout << "       2. elliptic edge     : \t" << m_numCandidateBallsOfElliptic[0] << "\t" << ((rg_REAL)m_numCandidateBallsOfElliptic[0])/numEdges[2][0] << endl;
    fout << endl << endl;


    fout << "   * actual search   : \t" << endl;
    fout << "                              " << "\t" << "position of end balls" << endl;
    fout << "          num of end balls  : " << "\t" << "in front FR" << "\t" << "from front FR" << "\t" << "in rear FR" << "\t" << "not found" << endl;
    fout << "       0. non-elliptic edge : " << "\t" << m_numEdges[0][0]+m_numEdges[1][0] << "\t" << m_numEdges[0][1]+m_numEdges[1][1] << "\t" << m_numEdges[0][2]+m_numEdges[1][2] <<  "\t" << numEdges[0][1]+numEdges[1][1] << endl;
    fout << "          1) parabolic edge : " << "\t" << m_numEdges[0][0] << "\t" << m_numEdges[0][1] << "\t" << m_numEdges[0][2] << "\t" << numEdges[0][1] << endl;
    fout << "          2) linear    edge : " << "\t" << m_numEdges[1][0] << "\t" << m_numEdges[1][1] << "\t" << m_numEdges[1][2] << "\t" << numEdges[1][1] << endl;
    fout << "       1. non-elliptic edge : " << "\t" << m_numCandidateBallsOfNonElliptic[1][0][0]+m_numCandidateBallsOfNonElliptic[1][1][0]
                                             << "\t" << m_numCandidateBallsOfNonElliptic[1][0][1]+m_numCandidateBallsOfNonElliptic[1][1][1]
                                             << "\t" << m_numCandidateBallsOfNonElliptic[1][0][2]+m_numCandidateBallsOfNonElliptic[1][1][2]
                                             << "\t" << m_numBallsForUnboundedEdge[0]+m_numBallsForUnboundedEdge[1];
    fout                                     << "\t" << ((rg_REAL) m_numCandidateBallsOfNonElliptic[1][0][0]+m_numCandidateBallsOfNonElliptic[1][1][0])/(m_numEdges[0][0]+m_numEdges[1][0])
                                             << "\t" << ((rg_REAL) m_numCandidateBallsOfNonElliptic[1][0][1]+m_numCandidateBallsOfNonElliptic[1][1][1])/(m_numEdges[0][1]+m_numEdges[1][1])
                                             << "\t" << ((rg_REAL) m_numCandidateBallsOfNonElliptic[1][0][2]+m_numCandidateBallsOfNonElliptic[1][1][2])/(m_numEdges[0][2]+m_numEdges[1][2])
                                             << "\t" << ((rg_REAL) m_numBallsForUnboundedEdge[0]+m_numBallsForUnboundedEdge[1])/(numEdges[0][1]+numEdges[1][1]) << endl;
    fout << "          1) parabolic edge : " << "\t" << m_numCandidateBallsOfNonElliptic[1][0][0]
                                             << "\t" << m_numCandidateBallsOfNonElliptic[1][0][1]
                                             << "\t" << m_numCandidateBallsOfNonElliptic[1][0][2]
                                             << "\t" << m_numBallsForUnboundedEdge[0];
    fout                                     << "\t" << ((rg_REAL) m_numCandidateBallsOfNonElliptic[1][0][0])/m_numEdges[0][0]
                                             << "\t" << ((rg_REAL) m_numCandidateBallsOfNonElliptic[1][0][1])/m_numEdges[0][1]
                                             << "\t" << ((rg_REAL) m_numCandidateBallsOfNonElliptic[1][0][2])/m_numEdges[0][2]
                                             << "\t" << ((rg_REAL) m_numBallsForUnboundedEdge[0])/(numEdges[0][1]) << endl;
    fout << "          2) linear    edge : " << "\t" << m_numCandidateBallsOfNonElliptic[1][1][0]
                                             << "\t" << m_numCandidateBallsOfNonElliptic[1][1][1]
                                             << "\t" << m_numCandidateBallsOfNonElliptic[1][1][2]
                                             << "\t" << m_numBallsForUnboundedEdge[1];
    fout                                     << "\t" << ((rg_REAL) m_numCandidateBallsOfNonElliptic[1][1][0])/m_numEdges[1][0]
                                             << "\t" << ((rg_REAL) m_numCandidateBallsOfNonElliptic[1][1][1])/m_numEdges[1][1]
                                             << "\t" << ((rg_REAL) m_numCandidateBallsOfNonElliptic[1][1][2])/m_numEdges[1][2]
                                             << "\t" << ((rg_REAL) m_numBallsForUnboundedEdge[1])/(numEdges[1][1]) << endl;
    fout << "       2. elliptic edge     : " << "\t" << m_numCandidateBallsOfElliptic[1] << "\t" << ((rg_REAL)m_numCandidateBallsOfElliptic[1])/numEdges[2][0] << endl;



    
    if ( m_bucket.isConstructed() )  {    
        fout << endl << "===  Bucket Analysis  ===" << endl;    
        fout << "   * element size                    : \t" << m_sizeOfBucketElement << endl;
        

        rg_INT numXs = m_bucket.getNumElementInX();
        rg_INT numYs = m_bucket.getNumElementInY();
        rg_INT numZs = m_bucket.getNumElementInZ();
        rg_INT numElements = numXs*numYs*numZs;
        rg_INT numEmptyElements = m_bucket.getNumEmptyElements();

        fout << "   * num elements                    : \t" << numElements << "(  " << numXs << "\t" <<  numYs << "\t" <<  numZs << "  )" << endl;
        fout << "   * num empty elements              : \t" << numEmptyElements << "\t" << ((rg_REAL)numEmptyElements)/numElements << endl;
        fout << "   * avg. num of balls in an element : \t" << ((rg_REAL)m_GlobalGeneratorList.getSize())/numElements << endl;
        fout << "   * avg. num of balls in non-empty element : \t" << ((rg_REAL)m_GlobalGeneratorList.getSize())/(numElements-numEmptyElements) << endl;
        fout << "   * avg. diameter of sphere for FRONT FEASIBLE REGION" << endl;
        fout << "                               " << "\t" << "position of end balls" << endl;
        fout << "                               " << "\t" << "in front FR" << "\t" << "from front FR" << "\t" << "in rear FR" << endl;
        fout << "       1. non-elliptic edge  : " << "\t" << 2.0*(m_avgRadiusOfSphereForRf[0][0]+m_avgRadiusOfSphereForRf[1][0])/(m_numEdges[0][0]+m_numEdges[1][0])
                                                  << "\t" << 2.0*(m_avgRadiusOfSphereForRf[0][1]+m_avgRadiusOfSphereForRf[1][1])/(m_numEdges[0][1]+m_numEdges[1][1])
                                                  << "\t" << 2.0*(m_avgRadiusOfSphereForRf[0][2]+m_avgRadiusOfSphereForRf[1][2])/(m_numEdges[0][2]+m_numEdges[1][2]) << endl;
        fout << "          1) parabolic edge  : " << "\t" << 2.0*m_avgRadiusOfSphereForRf[0][0]/m_numEdges[0][0]
                                                  << "\t" << 2.0*m_avgRadiusOfSphereForRf[0][1]/m_numEdges[0][1]
                                                  << "\t" << 2.0*m_avgRadiusOfSphereForRf[0][2]/m_numEdges[0][2] << endl;
        fout << "          2) linear edge     : " << "\t" << 2.0*m_avgRadiusOfSphereForRf[1][0]/m_numEdges[1][0]
                                                  << "\t" << 2.0*m_avgRadiusOfSphereForRf[1][1]/m_numEdges[1][1]
                                                  << "\t" << 2.0*m_avgRadiusOfSphereForRf[1][2]/m_numEdges[1][2] << endl;
        fout << endl << endl;
        fout << "   * avg. num of searched bucket elements" << endl;
        fout << "                               " << "\t" << "position of end balls" << endl;
        fout << "                               " << "\t" << "in front FR" << "\t" << "from front FR" << "\t" << "in rear FR" << "\t" << "not found" << endl;
        fout << "       1. non-elliptic edge  : " << "\t" << (m_numSearchedBucketElements[0][0]+m_numSearchedBucketElements[1][0])
                                                  << "\t" << ((rg_REAL)(m_numSearchedBucketElements[0][0]+m_numSearchedBucketElements[1][0]))/(m_numEdges[0][0]+m_numEdges[1][0])
                                                  << "\t" << (m_numSearchedBucketElements[0][1]+m_numSearchedBucketElements[1][1])
                                                  << "\t" << ((rg_REAL)(m_numSearchedBucketElements[0][1]+m_numSearchedBucketElements[1][1]))/(m_numEdges[0][1]+m_numEdges[1][1])
                                                  << "\t" << (m_numSearchedBucketElements[0][2]+m_numSearchedBucketElements[1][2])
                                                  << "\t" << ((rg_REAL)(m_numSearchedBucketElements[0][2]+m_numSearchedBucketElements[1][2]))/(m_numEdges[0][2]+m_numEdges[1][2])
                                                  << "\t" << (m_numSearchedBucketElements[0][3]+m_numSearchedBucketElements[1][3])
                                                  << "\t" << ((rg_REAL)(m_numSearchedBucketElements[0][3]+m_numSearchedBucketElements[1][3]))/(numEdges[0][1]+numEdges[1][1]) << endl;
        fout << "          1) parabolic edge  : " << "\t" << m_numSearchedBucketElements[0][0]
                                                  << "\t" << ((rg_REAL)m_numSearchedBucketElements[0][0])/m_numEdges[0][0]
                                                  << "\t" << m_numSearchedBucketElements[0][1]
                                                  << "\t" << ((rg_REAL)m_numSearchedBucketElements[0][1])/m_numEdges[0][1]
                                                  << "\t" << m_numSearchedBucketElements[0][2]
                                                  << "\t" << ((rg_REAL)m_numSearchedBucketElements[0][2])/m_numEdges[0][2]
                                                  << "\t" << (m_numSearchedBucketElements[0][3])
                                                  << "\t" << ((rg_REAL)(m_numSearchedBucketElements[0][3]))/(numEdges[0][1]) << endl;
        fout << "          2) linear edge     : " << "\t" << m_numSearchedBucketElements[1][0]
                                                  << "\t" << ((rg_REAL)m_numSearchedBucketElements[1][0])/m_numEdges[1][0]
                                                  << "\t" << m_numSearchedBucketElements[1][1]
                                                  << "\t" << ((rg_REAL)m_numSearchedBucketElements[1][1])/m_numEdges[1][1]
                                                  << "\t" << m_numSearchedBucketElements[1][2]
                                                  << "\t" << ((rg_REAL)m_numSearchedBucketElements[1][2])/m_numEdges[1][2]
                                                  << "\t" << (m_numSearchedBucketElements[1][3])
                                                  << "\t" << ((rg_REAL)(m_numSearchedBucketElements[1][3]))/(numEdges[1][1]) << endl;
    }
    //
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //*/
    
    return bValid;
}



void rg_SphereSetVoronoiDiagram::reportVDSforQT(ofstream& fout)
{
    fout << "# of Ball:   \t" << m_GlobalGeneratorList.getSize() << endl;
    fout << "# of Cell:   \t" << m_GlobalCellList.getSize() << endl;
    fout << "# of Face:   \t" << m_GlobalFaceList.getSize() << endl;
    fout << "# of Loop:   \t" << m_GlobalLoopList.getSize() << endl;
    fout << "# of Edge:   \t" << m_GlobalEdgeList.getSize() << endl;
    fout << "# of Vertex: \t" << m_GlobalVertexList.getSize() << endl;

    fout << "VD-BALL:\t" << m_GlobalGeneratorList.getSize() << endl;
    BallGenerator* currBall = rg_NULL;
    m_GlobalGeneratorList.reset4Loop();
    while ( m_GlobalGeneratorList.setNext4Loop() )  {
        currBall = m_GlobalGeneratorList.getpEntity();

        rg_Point3D pt = currBall->getBall().getCenter();
        fout << currBall->getCell()->getID()+1 << "\t"; 
        fout << pt.getX() << "\t" << pt.getY() << "\t" << pt.getZ() << "\t"; 
        fout << currBall->getBall().getRadius() << endl;
    }
    fout << endl;

    fout << "VD-VERTEX:\t" << m_GlobalVertexList.getSize() << endl;
    VDVertex* currVertex = rg_NULL;
    m_GlobalVertexList.reset4Loop();
    while ( m_GlobalVertexList.setNext4Loop() )  {
        currVertex = m_GlobalVertexList.getpEntity();

        fout << currVertex->getID() << " :\t";

        if ( currVertex->isOnInfinity() )
            fout << "X" << "\t";
        else
            fout << "O" << "\t";

        VDCell*  cellsDefiningVertex[4];
        currVertex->inquireAllCellsToDefineThisVertex(cellsDefiningVertex);
        rg_INT i=0;
		for ( i=0; i<4; i++ )
            fout << cellsDefiningVertex[i]->getID()+1 << "\t";


        VDEdge** incidentEdge = currVertex->getAllIncidentEdges();
        for ( i=0; i<4; i++ )  {
            VDCell* cellsDefiningEdge[3];
            incidentEdge[i]->inquireIntoOnlyEdgeSharingCellsInCCW( cellsDefiningEdge );

            fout << "( " << incidentEdge[i]->getStartVertex()->getID() << " -> "
                 << incidentEdge[i]->getEndVertex()->getID() << " : ";

            fout << cellsDefiningEdge[0]->getID()+1 << "\t";
            fout << cellsDefiningEdge[1]->getID()+1 << "\t";
            fout << cellsDefiningEdge[2]->getID()+1 << " )\t";
        }

        fout << endl;
    }
    fout << endl;
}




////////////////////////////////////////////////////////////////////////////////
//
//  Convert Voronoi diagram of balls into Quasi-triangulation
//
    
void rg_SphereSetVoronoiDiagram::checkLinkageVD2QT()
{
    ofstream fout("check_linkVD2QT.txt");

    m_GlobalCellList.reset4Loop();
    while ( m_GlobalCellList.setNext4Loop() ) {
        VDCell*     currVCell = m_GlobalCellList.getpEntity();

        BetaVertex* betaVtx   = currVCell->getBetaVertex();

        if ( betaVtx == rg_NULL ) {
            fout << "VCELL\t" << currVCell->getID() << "\t" << "no b-vtx" << endl;            
        }
        else {
            Ball*          currBallOnBVtx  = betaVtx->getBallProperty();
            BallGenerator* currBallOnVCell = currVCell->getGenerator();

            if ( currBallOnBVtx!=rg_NULL && currBallOnVCell!=rg_NULL ) {
                if ( !(currBallOnBVtx->getGeometry() == currBallOnVCell->getBall()) ) {
                    fout << "VCELL\t" << currVCell->getID() << "\t" << "diff. balls" << endl;            
                }
            }
            else if ( (currBallOnBVtx==rg_NULL && currBallOnVCell!=rg_NULL) || (currBallOnBVtx!=rg_NULL && currBallOnVCell==rg_NULL)) {
                fout << "VCELL\t" << currVCell->getID() << "\t" << "diff. balls" << endl;            
            }
            else {
            }
        }
    }


    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() ) {
        VDFace*   currVFace = m_GlobalFaceList.getpEntity();
        BetaEdge* betaEdge  = currVFace->getBetaEdge();

        if ( betaEdge == rg_NULL ) {
            fout << "VFACE\t" << currVFace->getID() << "\t" << "no b-edge" << endl;            
        }
    }


    m_GlobalEdgeList.reset4Loop();
    while ( m_GlobalEdgeList.setNext4Loop() ) {
        VDEdge*   currVEdge = m_GlobalEdgeList.getpEntity();
        BetaFace* betaFace  = currVEdge->getBetaFace();

        if ( betaFace == rg_NULL ) {
            fout << "VEDGE\t" << currVEdge->getID() << "\t" << "no b-face" << endl;            
        }
    }


    m_GlobalVertexList.reset4Loop();
    while ( m_GlobalVertexList.setNext4Loop() ) {
        VDVertex* currVVtx = m_GlobalVertexList.getpEntity();
        BetaCell* betaCell = currVVtx->getBetaCell();

        if ( currVVtx == rg_NULL ) {
            fout << "VVTX \t" << currVVtx->getID() << "\t" << "no b-cell" << endl;            
        }
        else {
            Sphere minTangentSphereOnBCell = betaCell->getMinTangentSphere();
            Sphere minTangentSphereOnVVtx( currVVtx->getPoint(), currVVtx->getRadiusOfTangentSphere() );
            if ( !(minTangentSphereOnBCell == minTangentSphereOnVVtx) )  {
                fout << "VVTX \t" << currVVtx->getID() << "\t" << "diff. min. tangent sphere." << endl;            
            }
        }
    }
}


void rg_SphereSetVoronoiDiagram::convertToQuasiTriangulation(QuasiTriangulation& quasiTriangulation)
{
    //  1. Initialization
    //    1) make objects for vertices and tetrahedra in QT.
    //    2) initialize them 
    QTVertex**      verticesInQT   = initializeVerticesInQT(quasiTriangulation);
    QTTetrahedron** tetrahedraInQT = initializeTetrahedraInQT(quasiTriangulation);


    //  2. Construction of Topology between Vertices and Tetrahedra
    //    1) vertex     : set first tetrahedron
    //    2) tetrahedron: set 4 vertices and neighbors
    //
    //    A Voronoi vertex is mapped into a tetrahedron in QT which is in a big or small world.
    //    Therefore, multiplicity anomaly is handled by this function 
    //    but singularity and degeneracy anomalies are not.
    setTopologyBetweenVerticesAndTetrahedraInQT( verticesInQT, tetrahedraInQT );
    


    setGatesInQT( verticesInQT, tetrahedraInQT, quasiTriangulation);

    linkVDandQT(  verticesInQT, tetrahedraInQT);

    quasiTriangulation.isLinkedVD( rg_TRUE );

    delete [] tetrahedraInQT;
    delete [] verticesInQT;
}




QTVertex** rg_SphereSetVoronoiDiagram::initializeVerticesInQT(QuasiTriangulation& quasiTriangulation)
{
    QTVertex** verticesInQT = new QTVertex*[m_GlobalCellList.getSize()];

    verticesInQT[0] = quasiTriangulation.addVertex( QTVertex(0) );

    BallGenerator* currBall = rg_NULL;
    m_GlobalGeneratorList.reset4Loop();
    while ( m_GlobalGeneratorList.setNext4Loop() )  {
        currBall = m_GlobalGeneratorList.getpEntity();

        rg_INT pos        = currBall->getCell()->getID() + 1;

        Ball*  ball       = quasiTriangulation.addBall( Ball(currBall->getBall(), currBall->getProperty(), currBall->getIDFromInput() ) );
        verticesInQT[pos] = quasiTriangulation.addVertex( QTVertex(pos, ball) );
    }

    return verticesInQT;
}



QTTetrahedron** rg_SphereSetVoronoiDiagram::initializeTetrahedraInQT(QuasiTriangulation& quasiTriangulation)
{
    QTTetrahedron** tetrahedraInQT = new QTTetrahedron*[m_GlobalVertexList.getSize()];

    VDVertex* currVertex = rg_NULL;
    m_GlobalVertexList.reset4Loop();
    while ( m_GlobalVertexList.setNext4Loop() )  {
        currVertex = m_GlobalVertexList.getpEntity();

        rg_INT pos = currVertex->getID();
        Sphere emptyTangentSphere( currVertex->getPoint(), currVertex->getRadiusOfTangentSphere() );

        tetrahedraInQT[pos] = quasiTriangulation.addTetrahedron( QTTetrahedron( pos, emptyTangentSphere ) );
    }

    return tetrahedraInQT;
}



void rg_SphereSetVoronoiDiagram::setTopologyBetweenVerticesAndTetrahedraInQT(QTVertex**      verticesInQT, 
                                                                             QTTetrahedron** tetrahedraInQT)
{
    VDVertex* currVertex = rg_NULL;
    m_GlobalVertexList.reset4Loop();
    while ( m_GlobalVertexList.setNext4Loop() )  {
        currVertex = m_GlobalVertexList.getpEntity();

        VDEdge** incidentEdge = currVertex->getAllIncidentEdges();

        VDCell*  cellsDefiningVertex[4];
        currVertex->inquireAllCellsToDefineThisVertex(cellsDefiningVertex);


        VDCell* cellsDefiningEdge[4][3];
        VDCell* cellForEndBallOfEdge[4];

        rg_INT i=0;
		for ( i=0; i<4; i++ )  {
            incidentEdge[i]->inquireIntoOnlyEdgeSharingCellsInCCW( cellsDefiningEdge[i] );
			
			rg_INT j=0;
            for ( j=0; j<4; j++)  {
                if (    (cellsDefiningVertex[j] != cellsDefiningEdge[i][0])
                     && (cellsDefiningVertex[j] != cellsDefiningEdge[i][1])
                     && (cellsDefiningVertex[j] != cellsDefiningEdge[i][2]) ) {
                    cellForEndBallOfEdge[i] = cellsDefiningVertex[j];
                    break;
                }
            }
        }

        //  Ordering of Vertices in Tetrahedron of Quasi-Triangulation
        //    For four vertices, v_0, v_1, v_2, and v_3,
        //    v_1, v_2, and v_3 are arranged in the CCW order when they are viewed from v_0.
        
        rg_INT pos = currVertex->getID();
        if ( currVertex == incidentEdge[0]->getEndVertex() )  {

            rg_INT idOfVertex[4] = { cellForEndBallOfEdge[0]->getID() + 1, 
                                     cellsDefiningEdge[0][0]->getID() + 1,
                                     cellsDefiningEdge[0][1]->getID() + 1, 
                                     cellsDefiningEdge[0][2]->getID() + 1 };

            tetrahedraInQT[pos]->setVertices( verticesInQT[idOfVertex[0]], verticesInQT[idOfVertex[1]],
                                              verticesInQT[idOfVertex[2]], verticesInQT[idOfVertex[3]] );
            for ( i=0; i<4; i++ )  {
                if ( verticesInQT[ idOfVertex[i] ]->getFirstTetrahedron() == rg_NULL )
                    verticesInQT[ idOfVertex[i] ]->setFirstTetrahedron( tetrahedraInQT[pos] );
            }

            rg_INT edgeToNeighbor[3];
            for ( i=0; i<3; i++ )  {
                if ( cellsDefiningEdge[0][i] == cellForEndBallOfEdge[1] )
                    edgeToNeighbor[i] = 1;
                else if ( cellsDefiningEdge[0][i] == cellForEndBallOfEdge[2] )
                    edgeToNeighbor[i] = 2;
                else
                    edgeToNeighbor[i] = 3;
            }

            rg_INT posOfNeighbor = incidentEdge[0]->getStartVertex()->getID();
            tetrahedraInQT[pos]->setNeighbor( 0, tetrahedraInQT[ posOfNeighbor ] );

            for ( i=0; i<3; i++ )  {
                if ( currVertex == incidentEdge[edgeToNeighbor[i]]->getEndVertex() )
                    posOfNeighbor = incidentEdge[edgeToNeighbor[i]]->getStartVertex()->getID();
                else
                    posOfNeighbor = incidentEdge[edgeToNeighbor[i]]->getEndVertex()->getID();

                tetrahedraInQT[pos]->setNeighbor( i+1, tetrahedraInQT[ posOfNeighbor ] );
            }
        }
        else  {

            rg_INT idOfVertex[4] = { cellForEndBallOfEdge[0]->getID() + 1, 
                                     cellsDefiningEdge[0][2]->getID() + 1,
                                     cellsDefiningEdge[0][1]->getID() + 1, 
                                     cellsDefiningEdge[0][0]->getID() + 1 };

            tetrahedraInQT[pos]->setVertices( verticesInQT[idOfVertex[0]], verticesInQT[idOfVertex[1]],
                                              verticesInQT[idOfVertex[2]], verticesInQT[idOfVertex[3]] );
            for ( i=0; i<4; i++ )  {
                if ( verticesInQT[ idOfVertex[i] ]->getFirstTetrahedron() == rg_NULL )
                    verticesInQT[ idOfVertex[i] ]->setFirstTetrahedron( tetrahedraInQT[pos] );
            }

            
            rg_INT edgeToNeighbor[3];
            for ( i=0; i<3; i++ )  {
                if ( cellsDefiningEdge[0][2-i] == cellForEndBallOfEdge[1] )
                    edgeToNeighbor[i] = 1;
                else if ( cellsDefiningEdge[0][2-i] == cellForEndBallOfEdge[2] )
                    edgeToNeighbor[i] = 2;
                else
                    edgeToNeighbor[i] = 3;
            }

            rg_INT posOfNeighbor = incidentEdge[0]->getEndVertex()->getID();
            tetrahedraInQT[pos]->setNeighbor( 0, tetrahedraInQT[ posOfNeighbor ] );

            for ( i=0; i<3; i++ )  {
                if ( currVertex == incidentEdge[edgeToNeighbor[i]]->getEndVertex() )
                    posOfNeighbor = incidentEdge[edgeToNeighbor[i]]->getStartVertex()->getID();
                else
                    posOfNeighbor = incidentEdge[edgeToNeighbor[i]]->getEndVertex()->getID();

                tetrahedraInQT[pos]->setNeighbor( i+1, tetrahedraInQT[ posOfNeighbor ] );
            }
        }
    }
}



void rg_SphereSetVoronoiDiagram::setGatesInQT(QTVertex**          verticesInQT, 
                                              QTTetrahedron**     tetrahedraInQT,
                                              QuasiTriangulation& quasiTriangulation)
{
    rg_INT numVerticesInQT = m_GlobalCellList.getSize();

    rg_dList< pair<VDEdge*, QTTetrahedron*> > refListForDegeneracyAnomaly;


    rg_dList< pair< pair<QTVertex*, QTVertex*>, rg_INT > > bigBrothersList;

    VDFace* currFace = rg_NULL;
    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isThereInnerLoop() )  {
            //  find bigBrothers.
            QTVertex* bigBrother[2] = { verticesInQT[ currFace->getRightCell()->getID()+1 ],
                                        verticesInQT[ currFace->getLeftCell()->getID()+1 ]  };
            QTGate* currGate = quasiTriangulation.addGate( QTGate( bigBrother[0], bigBrother[1] ) );


            pair< pair<QTVertex*, QTVertex*>, rg_INT > oneBigBrother;
            oneBigBrother.first.first  = bigBrother[0];
            oneBigBrother.first.second = bigBrother[1];
            oneBigBrother.second         = 1;
            bigBrothersList.add( oneBigBrother );

            
            //  find bigTetrahedron.
            VDLoop* outerLoop = currFace->getOuterLoop();

            QTTetrahedron* bigTetrahedron  = rg_NULL;
            if ( outerLoop->isSingleEdge() ) {
                VDEdge* vEdgeForDegeneracy = outerLoop->getPartialEdge()->getOriginalEdge();

                refListForDegeneracyAnomaly.reset4Loop();
                while ( refListForDegeneracyAnomaly.setNext4Loop() ) {
                    pair<VDEdge*, QTTetrahedron*> currRef = refListForDegeneracyAnomaly.getEntity();

                    if ( currRef.first == vEdgeForDegeneracy ) {
                        bigTetrahedron = currRef.second;
                        break;
                    }
                }
            }
            else {
                rg_dList<QTTetrahedron*> incidentTetrahedraOfBigBrother1;
                bigBrother[0]->inquireIncidentTetrahedra(incidentTetrahedraOfBigBrother1);

                QTTetrahedron* currTetrahedron = rg_NULL;
                incidentTetrahedraOfBigBrother1.reset4Loop();
                while ( incidentTetrahedraOfBigBrother1.setNext4Loop() ) {
                    currTetrahedron = incidentTetrahedraOfBigBrother1.getEntity();

                    if ( currTetrahedron->isThere(bigBrother[1]) )  {
                        bigTetrahedron = currTetrahedron;
                        break;
                    }
                }
            }
            currGate->setBigTetrahedron( bigTetrahedron );
            
            //  find small tetrahedra.
            rg_dList<VDLoop*>* loops = currFace->getLoops();

            VDLoop* currLoop = rg_NULL;
            loops->reset4Loop();
            while ( loops->setNext4Loop() )  {
                currLoop = loops->getEntity();
            
                if ( currLoop->isOuterLoop() )
                    continue;

                VDEdge* boundingEdge = currLoop->getPartialEdge()->getOriginalEdge();

                //if ( boundingEdge->getStartVertex() == rg_NULL ){
                if ( currLoop->isSingleEdge() ) {
                    //  Degeneracy anomaly
                    VDCell* cellsDefiningEdge[3];
                    boundingEdge->inquireIntoOnlyEdgeSharingCellsInCCW( cellsDefiningEdge );

                    QTTetrahedron* degeneracyAnomaly = quasiTriangulation.addTetrahedron( QTTetrahedron() );
                    degeneracyAnomaly->setID( quasiTriangulation.getNumQTTetrahedron()-1 );
                    for ( rg_INT i=0; i<3; i++ ) {
                        QTVertex* vertexQT = verticesInQT[ cellsDefiningEdge[i]->getID()+1 ];
                        degeneracyAnomaly->setVertex( i, vertexQT );
                        
                        if ( vertexQT->getFirstTetrahedron() == rg_NULL )
                            vertexQT->setFirstTetrahedron(degeneracyAnomaly);
                    }
                    
                    currGate->addSmallTetrahedron( degeneracyAnomaly );

                    pair<VDEdge*, QTTetrahedron*> referForDegeneracy;
                    referForDegeneracy.first  = boundingEdge;
                    referForDegeneracy.second = degeneracyAnomaly;
                    refListForDegeneracyAnomaly.add( referForDegeneracy );
                }
                else  {
                    QTTetrahedron* smallTetrahedron = tetrahedraInQT[ boundingEdge->getStartVertex()->getID() ];

                    currGate->addSmallTetrahedron( smallTetrahedron );
                }
            }
        }
    }


    //  set depth of big brothers.
    //  Note that bigBrothersList must be sorted by radii of big-brothers in the descending order.
    if ( bigBrothersList.getSize() > 0 )  {
        
        if ( bigBrothersList.getSize() > 1 )  {
        
            rg_dNode< pair< pair<QTVertex*, QTVertex*>, rg_INT > >* firstNode = bigBrothersList.getFirstpNode();
            rg_dNode< pair< pair<QTVertex*, QTVertex*>, rg_INT > >* lastNode  = bigBrothersList.getLastpNode();
            rg_dNode< pair< pair<QTVertex*, QTVertex*>, rg_INT > >* currNode  = firstNode->getNext();
            rg_dNode< pair< pair<QTVertex*, QTVertex*>, rg_INT > >* prevNode  = rg_NULL;

            do {
                QTVertex* smallBrother[2] = {currNode->getEntity().first.first, currNode->getEntity().first.second};
                Sphere    oneSmallBall = smallBrother[0]->getBall();

                for ( prevNode = currNode->getPrev(); prevNode!=lastNode; prevNode=prevNode->getPrev())  {
                    QTVertex* bigBrother[2] = { prevNode->getEntity().first.first, prevNode->getEntity().first.second };
                    Sphere    bigBall[2]    = { bigBrother[0]->getBall(), bigBrother[1]->getBall() };

                    Sphere    tangentSphere[2];
                    rg_INT numTS = makeCircumcircleOnCenterPlane( bigBall[0], bigBall[1], oneSmallBall, 
                                                                  tangentSphere[0], tangentSphere[1] );
                    
                    if ( numTS == 2 )  {
                        currNode->getpEntity()->second = (prevNode->getpEntity()->second + 1);
                        break;
                    }
                }

                currNode = currNode->getNext();
            } while ( currNode != firstNode );
        }    
        

        pair< pair<QTVertex*, QTVertex*>, rg_INT >* currBrother = rg_NULL;
        bigBrothersList.reset4Loop();
        while ( bigBrothersList.setNext4Loop() )  {
            currBrother = bigBrothersList.getpEntity();

            QTVertex* bigBrother[2] = { currBrother->first.first, currBrother->first.second };
            bigBrother[0]->setID( bigBrother[0]->getID() - (currBrother->second*numVerticesInQT) );
            bigBrother[1]->setID( bigBrother[1]->getID() - (currBrother->second*numVerticesInQT) );
        }
    }    
}



void rg_SphereSetVoronoiDiagram::linkVDandQT( QTVertex** verticesInQT, QTTetrahedron** tetrahedraInQT)
{
    linkVCellAndQTVertex( verticesInQT   );
    linkVVertexAndQTCell( tetrahedraInQT );
}


void rg_SphereSetVoronoiDiagram::linkVCellAndQTVertex( QTVertex**      verticesInQT)
{
    m_GlobalCellList.reset4Loop();
    while ( m_GlobalCellList.setNext4Loop() ) {
        VDCell* currVCell = m_GlobalCellList.getpEntity();

        QTVertex* currQVtx = verticesInQT[ currVCell->getID()+1 ];

        currVCell->connectQVertex( currQVtx );
        currQVtx->connectVCell(    currVCell );
    }
}



void rg_SphereSetVoronoiDiagram::linkVVertexAndQTCell( QTTetrahedron** tetrahedraInQT)
{
    m_GlobalVertexList.reset4Loop();
    while ( m_GlobalVertexList.setNext4Loop() ) {
        VDVertex* currVVtx = m_GlobalVertexList.getpEntity();

        QTTetrahedron* currQCell = tetrahedraInQT[ currVVtx->getID() ];

        currVVtx->connectQCell( currQCell );
        currQCell->connectVVertex( currVVtx );
    }
}



void rg_SphereSetVoronoiDiagram::reportNumericalErrorInTangentSphereComputationOf4Balls(ofstream& fout)
{
    VDVertex* currVertex = rg_NULL;
    m_GlobalVertexList.reset4Loop();
    while ( m_GlobalVertexList.setNext4Loop() )  {
        currVertex = m_GlobalVertexList.getpEntity();

        if ( currVertex->isOnInfinity() )
            continue;

        VDCell* cellToDefineVertex[4];
        currVertex->inquireAllCellsToDefineThisVertex( cellToDefineVertex );

        rg_INT  ballIDWithMinRadius = -1;
        rg_REAL minRadius           = DBL_MAX;

        BallGenerator* ballGenerator[4];
        Sphere         ball[4];
        rg_INT i=0;
		for ( i=0; i<4; i++ )  {
            ballGenerator[i] = (BallGenerator*) cellToDefineVertex[i]->getGenerator();
            ball[i]          = ballGenerator[i]->getBall();

            if ( ball[i].getRadius() < minRadius )  {
                minRadius = ball[i].getRadius();
                ballIDWithMinRadius = i;
            }
        }

        fout << "v_" << currVertex->getID() << "\t:\t" 
             << "b_" << ballGenerator[0]->getID() << "\t"
             << "b_" << ballGenerator[1]->getID() << "\t"
             << "b_" << ballGenerator[2]->getID() << "\t"
             << "b_" << ballGenerator[3]->getID() << endl;
        
	    if( (ball[0].getRadius() == ball[1].getRadius()) && (ball[1].getRadius() == ball[2].getRadius()) && (ball[2].getRadius() == ball[3].getRadius()) )
	    {
            for ( i=0; i<4; i++ )  {
                rg_INT j=0;
				for ( j=0; j<4; j++ )  {
                    if ( j==i ) continue;

					rg_INT k=0;
                    for ( k=0; k<4; k++ )  {
                        if ( k==i || k==j ) continue;

                        for ( rg_INT l=0; l<4; l++ )  {
                            if ( l==i || l==j || l==k ) continue;

                            fout << "\t" 
                                 << "b_" << ballGenerator[i]->getID() << "\t"
                                 << "b_" << ballGenerator[j]->getID() << "\t"
                                 << "b_" << ballGenerator[k]->getID() << "\t"
                                 << "b_" << ballGenerator[l]->getID() << "\t";
                            
                            rg_INT numTS;
                            Sphere tangentSpheres[2];
                            numTS = computeSphereTangentTo4SpheresWithSameRadiusOutside( ball[i], ball[j], ball[k], ball[l], tangentSpheres);
                            fout << "( " << numTS << " )\t";

                            if ( numTS == 2 )  {
                                if ( tangentSpheres[1] < tangentSpheres[0] )  {
                                    Sphere temp = tangentSpheres[0];
                                    tangentSpheres[0] = tangentSpheres[1];
                                    tangentSpheres[1] = temp;
                                }
                            }

                            for ( rg_INT iTS = 0; iTS<numTS; iTS++ ) {
                                rg_Point3D center = tangentSpheres[iTS].getCenter();
                                fout << "|\t";
                                fout.width(15);
                                fout << center.getX() << "\t";
                                fout.width(15);
                                fout << center.getY() << "\t";
                                fout.width(15);
                                fout << center.getZ() << "\t";
                                fout.width(15);
                                fout << tangentSpheres[iTS].getRadius() << "\t";
                            }
                            fout << endl;
                        }
                    }
                }
            }
	    }
	    else
	    {
            for ( i=0; i<4; i++ )  {
                rg_INT j=0;
				for ( j=0; j<4; j++ )  {
                    if ( j==i ) continue;

                    rg_INT k=0;
					for ( k=0; k<4; k++ )  {
                        if ( k==i || k==j ) continue;

                        for ( rg_INT l=0; l<4; l++ )  {
                            if ( l==i || l==j || l==k ) continue;

                            if ( ball[l].getRadius() != minRadius )
                                continue;

                            fout << "\t" 
                                 << "b_" << ballGenerator[i]->getID() << "\t"
                                 << "b_" << ballGenerator[j]->getID() << "\t"
                                 << "b_" << ballGenerator[k]->getID() << "\t"
                                 << "b_" << ballGenerator[l]->getID() << "\t";
                            
                            rg_INT numTS;
                            Sphere tangentSpheres[2];
                            numTS = computeSphereTangentTo4SpheresWithDiffRadiusOutside( ball[i], ball[j], ball[k], ball[l], tangentSpheres);
                            fout << "( " << numTS << " )\t";

                            if ( numTS == 2 )  {
                                if ( tangentSpheres[1] < tangentSpheres[0] )  {
                                    Sphere temp = tangentSpheres[0];
                                    tangentSpheres[0] = tangentSpheres[1];
                                    tangentSpheres[1] = temp;
                                }
                            }

                            for ( rg_INT iTS = 0; iTS<numTS; iTS++ ) {
                                rg_Point3D center = tangentSpheres[iTS].getCenter();
                                fout << "|\t";
                                fout.width(15);
                                fout << center.getX() << "\t";
                                fout.width(15);
                                fout << center.getY() << "\t";
                                fout.width(15);
                                fout << center.getZ() << "\t";
                                fout.width(15);
                                fout << tangentSpheres[iTS].getRadius() << "\t";
                                //rg_Point3D center = tangentSpheres[iTS].getCenter();
                                //fout << "|\t"
                                //     << center.getX() << "\t"
                                //     << center.getY() << "\t"
                                //     << center.getZ() << "\t"
                                //     << tangentSpheres[iTS].getRadius() << "\t";
                            }
                            fout << endl;
                        }
                    }
                }
            }
        }       
    }
}




void rg_SphereSetVoronoiDiagram::analyzeSphereVoronoiDiagram(ofstream& fout)
{
    fout << "###  Numbers in Sphere Voronoi Diagram  ###" << endl;
    fout << endl;
    fout << "Num. balls"      << "\t" << getNumOfBalls()                 << endl;
    fout << "Num. cells"      << "\t" << getNumOfVoronoiCells()          << endl;
    fout << "Num. UnB. cells" << "\t" << getNumOfUnboundedVoronoiCells() << endl;
    fout << "Num. faces"      << "\t" << getNumOfVoronoiFaces()          << endl;
    fout << "Num. edges"      << "\t" << getNumOfVoronoiEdges()          << endl;
    fout << "Num. vertices"   << "\t" << getNumOfVoronoiVertices()       << endl;
    fout << endl;
    fout << endl;


}


void rg_SphereSetVoronoiDiagram::analyzeNumNeighborAtoms(ofstream& fout)
{
    fout << endl;
    fout << "###  Frequency of Num. of Neighbor Atoms  (num. faces/cell) ###" << endl;

    rg_INT numBall = m_GlobalGeneratorList.getSize();
    rg_INT maxNumNeighborAtoms = -1;
    rg_INT minNumNeighborAtoms = numBall;
    rg_INT sumNumNeighborAtoms = 0;

    rg_INT i=0;
    rg_INT* distriNumNeighborAtoms = new rg_INT[numBall];
    for ( i=0; i<numBall; i++)
        distriNumNeighborAtoms[i] = 0;
    
    rg_INT count               = 0;
    BallGenerator* currBall = rg_NULL;
    m_GlobalGeneratorList.reset4Loop();
    while ( m_GlobalGeneratorList.setNext4Loop() )  {
        currBall = m_GlobalGeneratorList.getpEntity();

        VDCell* currCell = currBall->getCell();
        if ( currCell->isBounded() == rg_FALSE )
            continue;

        rg_INT  numNeighborAtoms = currCell->getNumOfBoundingFaces();

        if ( numNeighborAtoms > maxNumNeighborAtoms )
            maxNumNeighborAtoms = numNeighborAtoms;
        if ( numNeighborAtoms < minNumNeighborAtoms )
            minNumNeighborAtoms = numNeighborAtoms;

        sumNumNeighborAtoms += numNeighborAtoms;

        distriNumNeighborAtoms[count] = numNeighborAtoms;

        count++;
    }

    rg_REAL mean = (rg_REAL)sumNumNeighborAtoms/count;
    rg_REAL variance = 0.0;
    for (i=0; i<numBall; i++ )  {
        variance += (distriNumNeighborAtoms[i]-mean)*(distriNumNeighborAtoms[i]-mean);
    }
    variance = variance/count;
    rg_REAL standard_deviation = sqrt(variance);


    fout << endl;
    fout << "---  Summary  ---" << endl;
    fout << "*** Without the consideration of atoms with unbounded Voronoi regions" << endl;
    fout << "num. of tested balls"           << "\t" << count << endl;
    fout << "Average Num. of Neighbor Atoms" << "\t" << mean << endl;
    fout << "    Variance                  " << "\t" << variance << endl;
    fout << "    Standard deviation        " << "\t" << standard_deviation << endl;
    fout << "Maxinum Num. of Neighbor Atoms" << "\t" << maxNumNeighborAtoms << endl;     
    fout << "Mininum Num. of Neighbor Atoms" << "\t" << minNumNeighborAtoms << endl;     
    fout << endl;
        

    maxNumNeighborAtoms = -1;
    minNumNeighborAtoms = numBall;
    sumNumNeighborAtoms = 0;

    for ( i=0; i<numBall; i++)
        distriNumNeighborAtoms[i] = 0;
    
    count               = 0;
    m_GlobalGeneratorList.reset4Loop();
    while ( m_GlobalGeneratorList.setNext4Loop() )  {
        currBall = m_GlobalGeneratorList.getpEntity();

        VDCell* currCell = currBall->getCell();

        rg_INT  numNeighborAtoms = currCell->getNumOfBoundingFaces();

        if ( numNeighborAtoms > maxNumNeighborAtoms )
            maxNumNeighborAtoms = numNeighborAtoms;
        if ( numNeighborAtoms < minNumNeighborAtoms )
            minNumNeighborAtoms = numNeighborAtoms;

        sumNumNeighborAtoms += numNeighborAtoms;

        distriNumNeighborAtoms[count] = numNeighborAtoms;

        count++;
    }

    mean = (rg_REAL)sumNumNeighborAtoms/count;
    variance = 0.0;
    for (i=0; i<numBall; i++ )  {
        variance += (distriNumNeighborAtoms[i]-mean)*(distriNumNeighborAtoms[i]-mean);
    }
    variance = variance/count;
    standard_deviation = sqrt(variance);


    fout << endl;
    fout << endl;
    fout << "*** With the consideration of atoms with unbounded Voronoi regions" << endl;
    fout << "num. of tested balls"           << "\t" << count << endl;
    fout << "Average Num. of Neighbor Atoms" << "\t" << mean << endl;
    fout << "    Variance                  " << "\t" << variance << endl;
    fout << "    Standard deviation        " << "\t" << standard_deviation << endl;
    fout << "Maxinum Num. of Neighbor Atoms" << "\t" << maxNumNeighborAtoms << endl;     
    fout << "Mininum Num. of Neighbor Atoms" << "\t" << minNumNeighborAtoms << endl;     
    fout << endl << endl;


    delete [] distriNumNeighborAtoms;
}



void rg_SphereSetVoronoiDiagram::analyzeFrequencyOfNumEdgesPerFace(ofstream& fout)
{
    fout << endl;
    fout << "###  Frequency of Num. of Edges per Face  (num. edges/face) ###" << endl;
    rg_INT count       = 0;
    rg_INT maxNumEdges = -1;
    rg_INT minNumEdges = m_GlobalEdgeList.getSize();
    rg_INT sumNumEdges = 0;

    rg_INT numFace = m_GlobalFaceList.getSize();
    rg_INT* distriNumEdgesPerFace = new rg_INT[numFace];
    rg_INT i=0;
    for (i=0; i<numFace; i++)
        distriNumEdgesPerFace[i] = 0;

    VDFace* currFace = rg_NULL;
    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;
        if ( currFace->isBounded() == rg_FALSE )
            continue;

        rg_INT numEdges = 0;
        rg_dList<VDLoop*>* loops = currFace->getLoops();
        loops->reset4Loop();
        while ( loops->setNext4Loop() )
            numEdges += loops->getEntity()->getNumOfBoundingEdges();

        if ( numEdges > maxNumEdges )
            maxNumEdges = numEdges;
        if ( numEdges < minNumEdges )
            minNumEdges = numEdges;

        sumNumEdges += numEdges;
        
        distriNumEdgesPerFace[currFace->getID()] = numEdges;
        //fout << count << "\t" << numEdges << endl;
        count++;
    }

    rg_REAL mean = (rg_REAL)sumNumEdges/count;
    rg_REAL variance = 0.0;
    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;

        variance += (distriNumEdgesPerFace[currFace->getID()]-mean)*(distriNumEdgesPerFace[currFace->getID()]-mean);
    }
    variance = variance/count;
    rg_REAL standard_deviation = sqrt(variance);

    fout << endl;
    fout << "---  Summary  ---" << endl;
    fout << "*** Without the consideration of unbounded Voronoi faces" << endl;
    fout << "num. of tested faces"           << "\t" << count << endl;
    fout << "Average Num. of Edges " << "\t" << mean << endl;
    fout << "    Variance          " << "\t" << variance << endl;
    fout << "    Standard deviation" << "\t" << standard_deviation << endl;
    fout << "Maxinum Num. of Edges " << "\t" << maxNumEdges << endl;     
    fout << "Mininum Num. of Edges " << "\t" << minNumEdges << endl;     
    fout << endl;


    count       = 0;
    maxNumEdges = -1;
    minNumEdges = m_GlobalEdgeList.getSize();
    sumNumEdges = 0;

    for (i=0; i<numFace; i++)
        distriNumEdgesPerFace[i] = 0;

    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;

        rg_INT numEdges = 0;
        rg_dList<VDLoop*>* loops = currFace->getLoops();
        loops->reset4Loop();
        while ( loops->setNext4Loop() )
            numEdges += loops->getEntity()->getNumOfBoundingEdges();

        if ( numEdges > maxNumEdges )
            maxNumEdges = numEdges;
        if ( numEdges < minNumEdges )
            minNumEdges = numEdges;

        sumNumEdges += numEdges;
        
        distriNumEdgesPerFace[currFace->getID()] = numEdges;
        //fout << count << "\t" << numEdges << endl;
        count++;
    }

    mean = (rg_REAL)sumNumEdges/count;
    variance = 0.0;
    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;

        variance += (distriNumEdgesPerFace[currFace->getID()]-mean)*(distriNumEdgesPerFace[currFace->getID()]-mean);
    }
    variance = variance/count;
    standard_deviation = sqrt(variance);

    fout << endl;
    fout << "*** With the consideration of unbounded Voronoi faces" << endl;
    fout << "num. of tested faces"           << "\t" << count << endl;
    fout << "Average Num. of Edges " << "\t" << mean << endl;
    fout << "    Variance          " << "\t" << variance << endl;
    fout << "    Standard deviation" << "\t" << standard_deviation << endl;
    fout << "Maxinum Num. of Edges " << "\t" << maxNumEdges << endl;     
    fout << "Mininum Num. of Edges " << "\t" << minNumEdges << endl;     
    fout << endl << endl;


    delete [] distriNumEdgesPerFace;
    /*
    fout << endl;
    fout << "###  Frequency of Num. of Edges per Face  (num. edges/face) ###" << endl;
    rg_INT count       = 1;
    rg_INT maxNumEdges = -1;
    rg_INT minNumEdges = m_GlobalEdgeList.getSize();
    rg_INT sumNumEdges = 0;

    VDFace* currFace = rg_NULL;
    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;

        rg_INT numEdges = 0;
        rg_dList<VDLoop*>* loops = currFace->getLoops();
        loops->reset4Loop();
        while ( loops->setNext4Loop() )
            numEdges += loops->getEntity()->getNumOfBoundingEdges();

        if ( numEdges > maxNumEdges )
            maxNumEdges = numEdges;
        if ( numEdges < minNumEdges )
            minNumEdges = numEdges;

        sumNumEdges += numEdges;
        
        fout << count << "\t" << numEdges << endl;
        count++;
    }

    fout << endl;
    fout << "---  Summary  ---" << endl;
    fout << "Average Num. of Edges" << "\t" << sumNumEdges/count << endl;
    fout << "Maxinum Num. of Edges" << "\t" << maxNumEdges << endl;     
    fout << "Mininum Num. of Edges" << "\t" << minNumEdges << endl;     
    fout << endl << endl;
    */
}



void rg_SphereSetVoronoiDiagram::analyzeDistanceNeighborBalls(ofstream& fout)
{
    fout << endl;
    fout << "###  Frequency of Distances between Two Neighbor Balls  ###" << endl;
    rg_INT count       = 1;
    rg_REAL maxDistanceBTWCenter = -DBL_MAX;
    rg_REAL minDistanceBTWCenter = DBL_MAX;

    rg_REAL maxDistanceBTWBall   = -DBL_MAX;
    rg_REAL minDistanceBTWBall   = DBL_MAX;
    
    rg_REAL sumDistanceBTWCenter = 0.0;
    rg_REAL sumDistanceBTWBall   = 0.0;

    rg_INT numFace = m_GlobalFaceList.getSize();
    rg_REAL* distriDistanceBTWCenter = new rg_REAL[numFace];
    rg_REAL* distriDistanceBTWBall   = new rg_REAL[numFace];
    rg_INT i=0;
    for (i=0; i<numFace; i++)  {
        distriDistanceBTWCenter[i] = 0;
        distriDistanceBTWBall[i]   = 0;
    }

    VDFace* currFace = rg_NULL;
    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;

        Sphere ball[2] = { currFace->getRightCell()->getGenerator()->getBall(),
                           currFace->getLeftCell()->getGenerator()->getBall()  };
        
        rg_REAL distBtwCenter = ball[0].getCenter().distance( ball[1].getCenter() );
        rg_REAL distBtwBall   = ball[0].distanceToSphere( ball[1] );

        if ( distBtwCenter > maxDistanceBTWCenter )
            maxDistanceBTWCenter = distBtwCenter;
        if ( distBtwCenter < minDistanceBTWCenter )
            minDistanceBTWCenter = distBtwCenter;

        if ( distBtwBall > maxDistanceBTWBall )
            maxDistanceBTWBall = distBtwBall;
        if ( distBtwBall < minDistanceBTWBall )
            minDistanceBTWBall = distBtwBall;

        sumDistanceBTWCenter += distBtwCenter;
        sumDistanceBTWBall   += distBtwBall;

        distriDistanceBTWCenter[currFace->getID()] = distBtwCenter;
        distriDistanceBTWBall[currFace->getID()]   = distBtwBall;

        //fout << count << "\t" << distBtwCenter << "\t" << distBtwBall << endl;
        count++;
    }

    rg_REAL meanBtwCenter = (rg_REAL)sumDistanceBTWCenter/count;
    rg_REAL meanBtwBall   = (rg_REAL)sumDistanceBTWBall/count;
    rg_REAL varianceBtwCenter = 0.0;
    rg_REAL varianceBtwBall = 0.0;
    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;

        varianceBtwCenter +=  (distriDistanceBTWCenter[currFace->getID()]-meanBtwCenter)
                             *(distriDistanceBTWCenter[currFace->getID()]-meanBtwCenter);
        varianceBtwBall   +=  (distriDistanceBTWBall[currFace->getID()]-meanBtwBall)
                             *(distriDistanceBTWBall[currFace->getID()]-meanBtwBall);
    }
    varianceBtwCenter = varianceBtwCenter/count;
    varianceBtwBall   = varianceBtwBall/count;
    rg_REAL standard_deviationBtwCenter = sqrt(varianceBtwCenter);
    rg_REAL standard_deviationBtwBall = sqrt(varianceBtwBall);

    fout << endl;
    fout << "---  Summary  ---" << endl;
    fout << "Average Distance Btw Ball's Center" << "\t" << meanBtwCenter << endl;
    fout << "    Variance                      " << "\t" << varianceBtwCenter << endl;
    fout << "    Standard deviation            " << "\t" << standard_deviationBtwCenter << endl;
    fout << "Maxinum Distance Btw Ball's Center" << "\t" << maxDistanceBTWCenter << endl;     
    fout << "Mininum Distance Btw Ball's Center" << "\t" << minDistanceBTWCenter << endl;     
    fout << endl;
    fout << "Average Distance Btw Balls" << "\t" << meanBtwBall << endl;
    fout << "    Variance              " << "\t" << varianceBtwBall << endl;
    fout << "    Standard deviation    " << "\t" << standard_deviationBtwBall << endl;
    fout << "Maxinum Distance Btw Balls" << "\t" << maxDistanceBTWBall << endl;     
    fout << "Mininum Distance Btw Balls" << "\t" << minDistanceBTWBall << endl;     
    fout << endl << endl;

    delete [] distriDistanceBTWCenter;
    delete [] distriDistanceBTWBall;
}



void rg_SphereSetVoronoiDiagram::analyzeDistanceBtwBoundedNeighborBalls(ofstream& fout)
{
    fout << endl;
    fout << "###  Frequency of Distances between Two Bounded Neighbor Balls  ###" << endl;
    rg_INT count       = 1;
    rg_REAL maxDistanceBTWCenter = -DBL_MAX;
    rg_REAL minDistanceBTWCenter = DBL_MAX;

    rg_REAL maxDistanceBTWBall   = -DBL_MAX;
    rg_REAL minDistanceBTWBall   = DBL_MAX;
    
    rg_REAL sumDistanceBTWCenter = 0.0;
    rg_REAL sumDistanceBTWBall   = 0.0;

    rg_INT numFace = m_GlobalFaceList.getSize();
    rg_REAL* distriDistanceBTWCenter = new rg_REAL[numFace];
    rg_REAL* distriDistanceBTWBall   = new rg_REAL[numFace];
    rg_INT i=0;
    for (i=0; i<numFace; i++)  {
        distriDistanceBTWCenter[i] = 0.0;
        distriDistanceBTWBall[i]   = 0.0;
    }

    VDFace* currFace = rg_NULL;
    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;
        if ( currFace->isBounded() == rg_FALSE )
            continue;

        Sphere ball[2] = { currFace->getRightCell()->getGenerator()->getBall(),
                           currFace->getLeftCell()->getGenerator()->getBall()  };
        
        rg_REAL distBtwCenter = ball[0].getCenter().distance( ball[1].getCenter() );
        rg_REAL distBtwBall   = ball[0].distanceToSphere( ball[1] );

        if ( distBtwCenter > maxDistanceBTWCenter )
            maxDistanceBTWCenter = distBtwCenter;
        if ( distBtwCenter < minDistanceBTWCenter )
            minDistanceBTWCenter = distBtwCenter;

        if ( distBtwBall > maxDistanceBTWBall )
            maxDistanceBTWBall = distBtwBall;
        if ( distBtwBall < minDistanceBTWBall )
            minDistanceBTWBall = distBtwBall;

        sumDistanceBTWCenter += distBtwCenter;
        sumDistanceBTWBall   += distBtwBall;

        distriDistanceBTWCenter[currFace->getID()] = distBtwCenter;
        distriDistanceBTWBall[currFace->getID()]   = distBtwBall;

        //fout << count << "\t" << distBtwCenter << "\t" << distBtwBall << endl;
        count++;
    }

    rg_REAL meanBtwCenter = (rg_REAL)sumDistanceBTWCenter/count;
    rg_REAL meanBtwBall   = (rg_REAL)sumDistanceBTWBall/count;
    rg_REAL varianceBtwCenter = 0.0;
    rg_REAL varianceBtwBall = 0.0;
    m_GlobalFaceList.reset4Loop();
    while ( m_GlobalFaceList.setNext4Loop() )  {
        currFace = m_GlobalFaceList.getpEntity();

        if ( currFace->isOnInfinity() == rg_TRUE )
            continue;
        if ( currFace->isBounded() == rg_FALSE )
            continue;

        varianceBtwCenter +=  (distriDistanceBTWCenter[currFace->getID()]-meanBtwCenter)
                             *(distriDistanceBTWCenter[currFace->getID()]-meanBtwCenter);
        varianceBtwBall   +=  (distriDistanceBTWBall[currFace->getID()]-meanBtwBall)
                             *(distriDistanceBTWBall[currFace->getID()]-meanBtwBall);
    }
    varianceBtwCenter = varianceBtwCenter/count;
    varianceBtwBall   = varianceBtwBall/count;
    rg_REAL standard_deviationBtwCenter = sqrt(varianceBtwCenter);
    rg_REAL standard_deviationBtwBall = sqrt(varianceBtwBall);

    fout << endl;
    fout << "---  Summary  ---" << endl;
    fout << "Average Distance Btw Ball's Center" << "\t" << meanBtwCenter << endl;
    fout << "    Variance                      " << "\t" << varianceBtwCenter << endl;
    fout << "    Standard deviation            " << "\t" << standard_deviationBtwCenter << endl;
    fout << "Maxinum Distance Btw Ball's Center" << "\t" << maxDistanceBTWCenter << endl;     
    fout << "Mininum Distance Btw Ball's Center" << "\t" << minDistanceBTWCenter << endl;     
    fout << endl;
    fout << "Average Distance Btw Balls" << "\t" << meanBtwBall << endl;
    fout << "    Variance              " << "\t" << varianceBtwBall << endl;
    fout << "    Standard deviation    " << "\t" << standard_deviationBtwBall << endl;
    fout << "Maxinum Distance Btw Balls" << "\t" << maxDistanceBTWBall << endl;     
    fout << "Mininum Distance Btw Balls" << "\t" << minDistanceBTWBall << endl;     
    fout << endl << endl;

    delete [] distriDistanceBTWCenter;
    delete [] distriDistanceBTWBall;
}


