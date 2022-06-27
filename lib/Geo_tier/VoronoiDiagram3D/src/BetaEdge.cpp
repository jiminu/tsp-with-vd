#include "BetaEdge.h"

#include "BetaVertex.h"
#include "BetaFace.h"
#include "BetaCell.h"

#include "FunctionsForVoronoiDiagram3D.h"
#include "VDFace.h"
using namespace V::GeometryTier;


BetaEdge::BetaEdge()
: m_startVertex(rg_NULL), m_endVertex(rg_NULL), m_firstFace(rg_NULL),
  m_visited(rg_FALSE)
{
    m_vFace = rg_NULL;
}



BetaEdge::BetaEdge(BetaVertex* startVertex, BetaVertex* endVertex)
: m_startVertex(startVertex), m_endVertex(endVertex), m_firstFace(rg_NULL),
  m_visited(rg_FALSE)
{
    m_vFace = rg_NULL;
}



BetaEdge::BetaEdge(const rg_INT& ID, BetaVertex* startVertex, BetaVertex* endVertex)
: TopologicalEntity(ID),
  m_startVertex(startVertex), m_endVertex(endVertex), m_firstFace(rg_NULL),
  m_visited(rg_FALSE)
{
    m_vFace = rg_NULL;
}



BetaEdge::BetaEdge(const rg_INT& ID, BetaVertex* startVertex, BetaVertex* endVertex, BetaFace* firstFace)
: TopologicalEntity(ID),
  m_startVertex(startVertex), m_endVertex(endVertex), m_firstFace(firstFace),
  m_visited(rg_FALSE)
{
    m_vFace = rg_NULL;
}



BetaEdge::BetaEdge(const BetaEdge& betaEdge)
: TopologicalEntity(betaEdge.m_ID),
  m_startVertex(betaEdge.m_startVertex), m_endVertex(betaEdge.m_endVertex), m_firstFace(betaEdge.m_firstFace),
  m_betaSpan(betaEdge.m_betaSpan),
  m_visited(betaEdge.m_visited)
{
    m_vFace = betaEdge.m_vFace;
}



BetaEdge::~BetaEdge()
{
    m_startVertex = rg_NULL;
    m_endVertex   = rg_NULL;
    m_firstFace   = rg_NULL;

    m_smallWorlds.removeAll();


    //if ( m_vFace != rg_NULL ) {
    //    m_vFace->disconnectBetaEdge(this);
    //}
}



rg_BOOL BetaEdge::isVirtual() const
{
    if ( m_startVertex->isVirtual() || m_endVertex->isVirtual() ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



BetaVertex* BetaEdge::getStartVertex() const
{
    return m_startVertex;
}



BetaVertex* BetaEdge::getEndVertex() const
{
    return m_endVertex;
}



BetaFace*   BetaEdge::getFirstFace() const
{
    return m_firstFace;
}




rg_BOOL BetaEdge::isGateToSmallWorlds() const
{
    if ( m_smallWorlds.getSize() > 0 ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_INT  BetaEdge::getNumOfSmallWorlds() const
{
    return m_smallWorlds.getSize();
}



rg_dList<BetaFace*>* BetaEdge::getSmallWorlds() 
{
    return &m_smallWorlds;
}



rg_INT      BetaEdge::getDepth() const
{
    rg_INT depth = m_startVertex->getDepth();
    if ( depth < m_endVertex->getDepth() ) {
        depth = m_endVertex->getDepth();
    }

    return depth;
}


rg_BetaSpan  BetaEdge::getBetaSpan() const
{
    return m_betaSpan;
}



rg_REAL BetaEdge::getMinValueForValidSimplex() const
{
    return m_betaSpan.getLowerValueOfInterval(1);
}



void BetaEdge::setStartVertex(BetaVertex* qtVertex)
{
    m_startVertex = qtVertex;
}



void BetaEdge::setEndVertex(BetaVertex* qtVertex)
{
    m_endVertex = qtVertex;
}



void BetaEdge::setFirstFace(BetaFace* firstFace)
{
    m_firstFace = firstFace;
}



void BetaEdge::addSmallWorld(BetaFace* smallWorld)
{
    m_smallWorlds.add( smallWorld );
}



void BetaEdge::setBetaSpan(const rg_BetaSpan& betaSpan)
{
    m_betaSpan = betaSpan;
}

void BetaEdge::shiftBetaSpan(const rg_REAL& delta)
{
	m_betaSpan.shiftBetaSpan(delta);
}

BetaEdge& BetaEdge::operator=(const BetaEdge& betaEdge)
{
    if ( this == &betaEdge ) {
        return *this;
    }

    m_ID          = betaEdge.m_ID;
    m_startVertex = betaEdge.m_startVertex;
    m_endVertex   = betaEdge.m_endVertex;
    m_firstFace   = betaEdge.m_firstFace;
    m_betaSpan    = betaEdge.m_betaSpan;

    m_visited     = betaEdge.m_visited;

    m_vFace = betaEdge.m_vFace;

    return *this;
}



rg_BOOL BetaEdge::isOnConvexHull() const
{
    rg_dList<BetaCell*> cellList;
    searchCellsInIntraWorld( cellList );

    BetaCell* currCell = rg_NULL;
    cellList.reset4Loop();
    while ( cellList.setNext4Loop() ) {
        currCell = cellList.getEntity();

        if ( currCell->isVirtual() ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_BOOL BetaEdge::isThere(const BetaVertex* const vertex) const
{
    if ( vertex == m_startVertex || vertex == m_endVertex ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



BetaCell* BetaEdge::getCCWNextCell(const BetaFace* const currFace) const
{
    if ( currFace->isThere( this ) ) {
        BetaCell* ccwNextCell = rg_NULL;
        if ( currFace->getEdgeOrientation( this ) ) {
            ccwNextCell = currFace->getRightCell();
        }
        else {
            ccwNextCell = currFace->getLeftCell();
        }

        return ccwNextCell;
    }
    else {
        return rg_NULL;
    }
}



BetaFace* BetaEdge::getCCWNextFace(const BetaFace* const currFace) const
{
    BetaCell* ccwNextCell = getCCWNextCell( currFace );

    if ( ccwNextCell != rg_NULL ) {
        BetaFace* incidentFace[2];
        ccwNextCell->findFacesIncidentToEdge( this, incidentFace);
     
        BetaFace* ccwNextFace = rg_NULL;
        if ( currFace == incidentFace[0] ) {
            ccwNextFace = incidentFace[1];
        }
        else {
            ccwNextFace = incidentFace[0];
        }

        return ccwNextFace;
    }
    else {
        return rg_NULL;
    }
}



BetaCell* BetaEdge::getCWNextCell(const BetaFace* const currFace) const
{
    if ( currFace->isThere( this ) ) {
        BetaCell* cwNextCell = rg_NULL;
        if ( currFace->getEdgeOrientation( this ) ) {
            cwNextCell = currFace->getLeftCell();
        }
        else {
            cwNextCell = currFace->getRightCell();
        }

        return cwNextCell;
    }
    else {
        return rg_NULL;
    }
}



BetaFace* BetaEdge::getCWNextFace(const BetaFace* const currFace) const
{
    BetaCell* cwNextCell = getCWNextCell( currFace );

    if ( cwNextCell != rg_NULL ) {
        BetaFace* incidentFace[2];
        cwNextCell->findFacesIncidentToEdge( this, incidentFace);
     
        BetaFace* cwNextFace = rg_NULL;
        if ( currFace == incidentFace[0] ) {
            cwNextFace = incidentFace[1];
        }
        else {
            cwNextFace = incidentFace[0];
        }

        return cwNextFace;
    }
    else {
        return rg_NULL;
    }
}



void BetaEdge::computeBetaSpan()
{
    Sphere minTangentSphere = computeMinSphereTangentToTwoBalls(m_startVertex->getBall(), m_endVertex->getBall());

    // construct beta-span in the big world.
    m_betaSpan = computeBetaSpan(minTangentSphere);


    BetaFace* currSmallWorld = rg_NULL;
    m_smallWorlds.reset4Loop();
    while ( m_smallWorlds.setNext4Loop() ) {
        currSmallWorld = m_smallWorlds.getEntity();

        // construct beta-span in the small world.
        rg_BetaSpan betaSpanInSmallWorld = computeBetaSpan(minTangentSphere, currSmallWorld);

        m_betaSpan.unify( betaSpanInSmallWorld );
    }
}



rg_BetaSpan BetaEdge::computeBetaSpan(const Sphere& minTangentSphere, BetaFace* world) const
{
    rg_REAL minOfCurrSimplex = minTangentSphere.getRadius();

    rg_REAL minOfHigherSimplex = REAL_PLUS_INFINITE;
    rg_REAL maxOfHigherSimplex = REAL_MINUS_INFINITE;

    
    rg_dList<BetaFace*> incidentFaces;
    searchFiniteFacesInIntraWorld( incidentFaces, world );
    
    BetaFace* currFace = rg_NULL;
    incidentFaces.reset4Loop();
    while ( incidentFaces.setNext4Loop() )  {
        currFace = incidentFaces.getEntity();

        Interval_REAL interval = currFace->getBetaSpan().getBetaIntervalOfSingularRegularState();

        rg_REAL minValue = interval.getLowerValue();
        rg_REAL maxValue = interval.getUpperValue();
        if ( minValue < minOfHigherSimplex ) {
            minOfHigherSimplex = minValue;
        }

        if ( maxValue > maxOfHigherSimplex ) {
            maxOfHigherSimplex = maxValue;
        }                
    }



    rg_BetaSpan betaSpan;
    if ( world == rg_NULL ) {
        // In the big world
        if ( isOnConvexHull() ) {
            if ( isAttached(minTangentSphere) ) {
                betaSpan.setNumBetaInterVal(2);
                betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
                betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  REAL_PLUS_INFINITE);
            }
            else {
                betaSpan.setNumBetaInterVal(3);
                betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
                betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
                betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  REAL_PLUS_INFINITE);
            }
        }
        else {
            if ( isAttached(minTangentSphere) ) {
                betaSpan.setNumBetaInterVal(3);
                betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
                betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
                betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);
            }
            else {
                betaSpan.setNumBetaInterVal(4);
                betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
                betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
                betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
                betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);

            }
        }
    }
    else {
        // In the small world
        if ( isAttached(minTangentSphere, world) ) {
            betaSpan.setNumBetaInterVal(3);
            betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
            betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
            betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);
        }
        else {
            betaSpan.setNumBetaInterVal(4);
            betaSpan.setBetaInterval(0, EXTRANEOUS_SIMPLEX, REAL_MINUS_INFINITE, minOfCurrSimplex);
            betaSpan.setBetaInterval(1, SINGULAR_SIMPLEX,   minOfCurrSimplex,    minOfHigherSimplex);
            betaSpan.setBetaInterval(2, REGULAR_SIMPLEX,    minOfHigherSimplex,  maxOfHigherSimplex);
            betaSpan.setBetaInterval(3, INTERIOR_SIMPLEX,   maxOfHigherSimplex,  REAL_PLUS_INFINITE);
        }
    }

    return betaSpan;
}



rg_BOOL BetaEdge::isAttached(const Sphere& minTangentSphere, BetaFace* world) const
{
    // check for whether two q-vertices of this q-edge 
    // define only a single q-edge in the whole quasi-triangulation or not.
    rg_dList<BetaVertex*> neighborVertexOfStartVertex;
    m_startVertex->searchVerticesInIntraWorld( neighborVertexOfStartVertex, world );

    rg_INT numEdgesWithSameVertices = 0;
    BetaVertex* neighborVertex = rg_NULL;
    neighborVertexOfStartVertex.reset4Loop();
    while ( neighborVertexOfStartVertex.setNext4Loop() ) {
        neighborVertex = neighborVertexOfStartVertex.getEntity();

        if ( m_endVertex == neighborVertex ) {
            numEdgesWithSameVertices++;
        }
    }

    //  When numEdgesWithSameVertices == 1, 
    //  two q-vertices of this q-edge 
    //  define only a single q-edge, this edge, in the whole quasi-triangulation.
    //  Therefore, the emptiness test of min. tangent sphere($\xi$) is needed only.
    if ( numEdgesWithSameVertices >= 1 ) {
        //  check the emptiness of minimum tangent sphere(minTangentSphere).
        //rg_dList<BetaFace*> incidentFaces;
        //searchFiniteFacesInIntraWorld( incidentFaces, world );
    
        //BetaFace* currFace = rg_NULL;
        //incidentFaces.reset4Loop();
        //while ( incidentFaces.setNext4Loop() )  {
        //    currFace = incidentFaces.getEntity();
    
        //    BetaVertex* mateVertex = currFace->findMateVertexOfEdge( this );
        //
        //    if ( minTangentSphere.isThereIntersectionWith( mateVertex->getBall() ) )  {
        //        return rg_TRUE;
        //    }                         
        //}

        neighborVertexOfStartVertex.reset4Loop();
        while ( neighborVertexOfStartVertex.setNext4Loop() ) {
            neighborVertex = neighborVertexOfStartVertex.getEntity();

            if ( m_endVertex == neighborVertex ) {
                continue;
            }

            if ( minTangentSphere.isThereIntersectionWith( neighborVertex->getBall() ) )  {
                return rg_TRUE;
            }                         
        }

        rg_dList<BetaVertex*> neighborVertexOfEndVertex;
        m_endVertex->searchVerticesInIntraWorld( neighborVertexOfEndVertex, world );

        neighborVertexOfEndVertex.reset4Loop();
        while ( neighborVertexOfEndVertex.setNext4Loop() ) {
            neighborVertex = neighborVertexOfEndVertex.getEntity();

            if ( m_startVertex == neighborVertex ) {
                continue;
            }

            if ( minTangentSphere.isThereIntersectionWith( neighborVertex->getBall() ) )  {
                return rg_TRUE;
            }                         
        }

    }

    //  When numEdgesWithSameVertices > 1, 
    //  two q-vertices of this q-edge define distinct multiple q-edges in the whole quasi-triangulation.
    //  In the primary structure (VD), this means two atoms define multiple v-faces.
    if ( numEdgesWithSameVertices > 1 ) {
        //  check for whether the center of minTangentSphere is on the Voronoi face or not.
        rg_dList<BetaCell*> incidentCells;
        searchFiniteCellsInIntraWorld( incidentCells, world );

        BetaCell* currCell = rg_NULL;
        incidentCells.reset4Loop();
        while ( incidentCells.setNext4Loop() ) {
            currCell = incidentCells.getEntity();

            rg_REAL signedVolume = currCell->computeSignedVolume();
            if ( signedVolume < 0.0 ) {
                return rg_TRUE;
            }
        }
    }

    return rg_FALSE;
}


///////////////////////////////////////////////////////////////////////////
//  Quasi-operator
//  
//  1. Q(e, Y | W): Query primitives for intra-world.
void BetaEdge::searchVerticesInIntraWorld( rg_dList<BetaVertex*>& vertexList, BetaFace* world ) const
{
    vertexList.add( m_startVertex );
    vertexList.add( m_endVertex );
}



void BetaEdge::searchEdgesInIntraWorld(    rg_dList<BetaEdge*>&   edgeList,   BetaFace* world ) const
{
    BetaFace* faceInCurrWorld = rg_NULL;
    // if world is null, find the biggest world of the vertex.
    if ( world == rg_NULL ) {
        faceInCurrWorld = m_firstFace;
    }
    else {
        if ( world->isThere( this ) ) {
            faceInCurrWorld = world;
        }
        else {
            return;
        }
    }


    // The current world of the vertex is an isolated face.
    if ( faceInCurrWorld->isIsolatedFace() ) {
        BetaEdge** edgeOnFace = faceInCurrWorld->getEdges();
        for ( rg_INT i_edge=0; i_edge<EIWDS_NUM_EDGE_ON_FACE; i_edge++ ) {
            if ( this != edgeOnFace[i_edge] ) {
                edgeList.add( edgeOnFace[i_edge] );
            }
        }
    }
    else {
        BetaFace* startFace = faceInCurrWorld;
        BetaFace* currFace  = startFace;

        do {
            BetaEdge** edgeOnFace = currFace->getEdges();
            for ( rg_INT i_edge=0; i_edge<EIWDS_NUM_EDGE_ON_FACE; i_edge++ ) {
                if ( edgeOnFace[i_edge] != this ) {
                    edgeList.add( edgeOnFace[i_edge] );
                }
            }

            currFace  = getCCWNextFace( currFace );

        } while ( currFace != startFace );
    }
}




void BetaEdge::searchFacesInIntraWorld(    rg_dList<BetaFace*>&   faceList,   BetaFace* world ) const
{
    BetaFace* faceInCurrWorld = rg_NULL;
    // if world is null, find the biggest world of the vertex.
    if ( world == rg_NULL ) {
        faceInCurrWorld = m_firstFace;
    }
    else {
        if ( world->isThere( this ) ) {
            faceInCurrWorld = world;
        }
        else {
            return;
        }
    }


    // The current world of the vertex is an isolated face.
    if ( faceInCurrWorld->isIsolatedFace() ) {
        faceList.add( faceInCurrWorld );
    }
    else {
        BetaFace* startFace = faceInCurrWorld;
        BetaFace* currFace  = startFace;

        do {
            faceList.add( currFace );

            currFace  = getCCWNextFace( currFace );

        } while ( currFace != startFace );
    }
}



void BetaEdge::searchCellsInIntraWorld(    rg_dList<BetaCell*>&   cellList,   BetaFace* world ) const
{
    BetaFace* faceInCurrWorld = rg_NULL;
    // if world is null, find the biggest world of the vertex.
    if ( world == rg_NULL ) {
        faceInCurrWorld = m_firstFace;
    }
    else {
        if ( world->isThere( this ) ) {
            faceInCurrWorld = world;
        }
        else {
            return;
        }
    }


    // The current world of the vertex is an isolated face.
    if ( faceInCurrWorld->isIsolatedFace() ) {
        return;
    }
    else {
        BetaFace* startFace = faceInCurrWorld;
        BetaFace* currFace  = startFace;

        do {
            BetaCell* currCell = getCCWNextCell(currFace);
            cellList.add( currCell );

            currFace  = getCCWNextFace( currFace );

        } while ( currFace != startFace );
    }
}



void BetaEdge::searchFiniteEdgesInIntraWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world ) const
{
    rg_dList<BetaEdge*> edgeList;
    searchEdgesInIntraWorld( edgeList, world );

    BetaEdge* currEdge = rg_NULL;
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        currEdge = edgeList.getEntity();

        if ( currEdge->isVirtual() == rg_FALSE ) {
            finiteEdgeList.add( currEdge );
        }
    }
}



void BetaEdge::searchFiniteFacesInIntraWorld(    rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world ) const
{
    rg_dList<BetaFace*> faceList;
    searchFacesInIntraWorld( faceList, world );

    BetaFace* currFace = rg_NULL;
    faceList.reset4Loop();
    while ( faceList.setNext4Loop() ) {
        currFace = faceList.getEntity();

        if ( currFace->isVirtual() == rg_FALSE ) {
            finiteFaceList.add( currFace );
        }
    }
}



void BetaEdge::searchFiniteCellsInIntraWorld(    rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world ) const
{
    rg_dList<BetaCell*> cellList;
    searchCellsInIntraWorld( cellList, world );

    BetaCell* currCell = rg_NULL;
    cellList.reset4Loop();
    while ( cellList.setNext4Loop() ) {
        currCell = cellList.getEntity();

        if ( currCell->isVirtual() == rg_FALSE ) {
            finiteCellList.add( currCell );
        }
    }
}



//  2. Query primitives for inter-world.
void BetaEdge::searchSmallWorlds( rg_dList<BetaFace*>& smallWorldList, BetaFace* world ) const
{
    rg_BOOL isOperationInBigWorld = rg_FALSE;
    if ( world == rg_NULL ) {
        isOperationInBigWorld = rg_TRUE;
    }
    else {
        if ( world->isThere(this) ) {
            // check for whether "world" is in the biggest world of this edge.
            rg_dList<BetaFace*> facesOnBiggestWorld;
            searchFacesInIntraWorld( facesOnBiggestWorld );

            facesOnBiggestWorld.reset4Loop();
            while ( facesOnBiggestWorld.setNext4Loop() ) {
                facesOnBiggestWorld.getEntity()->isVisited( rg_TRUE );
            }

            if ( world->isVisited() ) {
                isOperationInBigWorld = rg_TRUE;
            }

            facesOnBiggestWorld.reset4Loop();
            while ( facesOnBiggestWorld.setNext4Loop() ) {
                facesOnBiggestWorld.getEntity()->isVisited( rg_FALSE );
            }
        }
    }


    if ( isOperationInBigWorld ) {
        smallWorldList.append( m_smallWorlds );
    }
}



BetaFace* BetaEdge::searchBigWorld( BetaFace* world ) const
{
    rg_BOOL isOperationInSmallWorld = rg_FALSE;
    if ( world == rg_NULL ) {
        isOperationInSmallWorld = rg_FALSE;
    }
    else {
        if ( world->isThere(this) ) {
            // check for whether "world" is in the small world of this edge.
            rg_dList<BetaFace*> facesOnBiggestWorld;
            searchFacesInIntraWorld( facesOnBiggestWorld );

            facesOnBiggestWorld.reset4Loop();
            while ( facesOnBiggestWorld.setNext4Loop() ) {
                facesOnBiggestWorld.getEntity()->isVisited( rg_TRUE );
            }

            if ( !world->isVisited() ) {
                isOperationInSmallWorld = rg_TRUE;
            }

            facesOnBiggestWorld.reset4Loop();
            while ( facesOnBiggestWorld.setNext4Loop() ) {
                facesOnBiggestWorld.getEntity()->isVisited( rg_FALSE );
            }
        }
    }


    if ( isOperationInSmallWorld ) {
        return m_firstFace;
    }
    else {
        return rg_NULL;
    }

}



//  3. Q(e, Y): Query primitives for whole-world.
void BetaEdge::searchVerticesInWholeWorld( rg_dList<BetaVertex*>& vertexList ) const
{
    searchVerticesInIntraWorld( vertexList );
}



void BetaEdge::searchEdgesInWholeWorld(    rg_dList<BetaEdge*>&   edgeList   ) const
{
    rg_dList<BetaFace*> contributingWorlds;
    contributingWorlds.add( m_firstFace );

    searchSmallWorlds( contributingWorlds );

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        searchEdgesInIntraWorld( edgeList, currWorld );
    }
}



void BetaEdge::searchFacesInWholeWorld(    rg_dList<BetaFace*>&   faceList   ) const
{
    rg_dList<BetaFace*> contributingWorlds;
    contributingWorlds.add( m_firstFace );

    searchSmallWorlds( contributingWorlds );

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        searchFacesInIntraWorld( faceList, currWorld );
    }
}



void BetaEdge::searchCellsInWholeWorld(    rg_dList<BetaCell*>&   cellList   ) const
{
    rg_dList<BetaFace*> contributingWorlds;
    contributingWorlds.add( m_firstFace );

    searchSmallWorlds( contributingWorlds );

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        searchCellsInIntraWorld( cellList, currWorld );
    }
}



void BetaEdge::searchFiniteEdgesInWholeWorld(    rg_dList<BetaEdge*>&   finiteEdgeList   ) const
{
    rg_dList<BetaEdge*> edgeList;
    searchEdgesInWholeWorld( edgeList );

    BetaEdge* currEdge = rg_NULL;
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        currEdge = edgeList.getEntity();

        if ( currEdge->isVirtual() == rg_FALSE ) {
            finiteEdgeList.add( currEdge );
        }
    }
}



void BetaEdge::searchFiniteFacesInWholeWorld(    rg_dList<BetaFace*>&   finiteFaceList   ) const
{
    rg_dList<BetaFace*> faceList;
    searchFacesInWholeWorld( faceList );

    BetaFace* currFace = rg_NULL;
    faceList.reset4Loop();
    while ( faceList.setNext4Loop() ) {
        currFace = faceList.getEntity();

        if ( currFace->isVirtual() == rg_FALSE ) {
            finiteFaceList.add( currFace );
        }
    }
}



void BetaEdge::searchFiniteCellsInWholeWorld(    rg_dList<BetaCell*>&   finiteCellList   ) const
{
    rg_dList<BetaCell*> cellList;
    searchCellsInWholeWorld( cellList );

    BetaCell* currCell = rg_NULL;
    cellList.reset4Loop();
    while ( cellList.setNext4Loop() ) {
        currCell = cellList.getEntity();

        if ( currCell->isVirtual() == rg_FALSE ) {
            finiteCellList.add( currCell );
        }
    }
}




//
//  END of Quasi-operator
///////////////////////////////////////////////////////////////////////////};



///////////////////////////////////////////////////////////////////////////
//  Beta-operator
//  
void BetaEdge::searchEdgesInBetaComplex( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList ) const
{
    rg_dList<BetaEdge*> finiteEdgeList;
    searchFiniteEdgesInWholeWorld( finiteEdgeList );

    BetaEdge* currEdge = rg_NULL;
    finiteEdgeList.reset4Loop();
    while ( finiteEdgeList.setNext4Loop() ) {
        currEdge = finiteEdgeList.getEntity();

        if ( currEdge->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            continue;
        }
        else {
            betaEdgeList.add( currEdge );
        }
    }
}



void BetaEdge::searchFacesInBetaComplex( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const
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



void BetaEdge::searchCellsInBetaComplex( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList ) const
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



void BetaEdge::searchExtraneousEdgesInBetaComplex( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList ) const
{
    rg_dList<BetaEdge*> finiteEdgeList;
    searchFiniteEdgesInWholeWorld( finiteEdgeList );

    BetaEdge* currEdge = rg_NULL;
    finiteEdgeList.reset4Loop();
    while ( finiteEdgeList.setNext4Loop() ) {
        currEdge = finiteEdgeList.getEntity();

        if ( currEdge->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            betaEdgeList.add( currEdge );
        }
    }
}



void BetaEdge::searchExtraneousFacesInBetaComplex( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const
{
    rg_dList<BetaFace*> finiteFaceList;
    searchFiniteFacesInWholeWorld( finiteFaceList );

    BetaFace* currFace = rg_NULL;
    finiteFaceList.reset4Loop();
    while ( finiteFaceList.setNext4Loop() ) {
        currFace = finiteFaceList.getEntity();

        if ( currFace->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            betaFaceList.add( currFace );
        }
    }
}



void BetaEdge::searchExtraneousCellsInBetaComplex( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList ) const
{
    rg_dList<BetaCell*> finiteCellList;
    searchFiniteCellsInWholeWorld( finiteCellList );

    BetaCell* currCell = rg_NULL;
    finiteCellList.reset4Loop();
    while ( finiteCellList.setNext4Loop() ) {
        currCell = finiteCellList.getEntity();

        if ( currCell->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
            betaCellList.add( currCell );
        }
    }
}





void BetaEdge::searchEdgesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaEdge*>& betaEdgeList ) const
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



void BetaEdge::searchFacesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaFace*>& betaFaceList ) const
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


    if ( !isGateToSmallWorlds() || betaFaceList.getSize() <= 2 ) {
        return;
    }


    ///////////////////////////////////////////////////////////////////////////
    //  Some applications need the incident faces to be sorted.
    //
    //  1. generate plane for computing angle of beta-face.
    rg_Point3D startPt = m_startVertex->getBall().getCenter();
    rg_Point3D endPt   = m_endVertex->getBall().getCenter();

    rg_Point3D normal = endPt - startPt;
    normal.normalize();

    Plane basePlane;
    basePlane.definePlaneByNormalAndPassingPoint(normal, startPt);



    //  2. construct angle of beta-face.
    rg_INT numFaces  = betaFaceList.getSize();

    BetaFace**  faces        = betaFaceList.getArray();
    rg_Point3D* projectedVtx = new rg_Point3D[numFaces];
    rg_REAL*    angle        = new rg_REAL[numFaces];
    
    rg_INT i=0;
    for ( i=0; i<numFaces; i++ ) {
        BetaVertex* mateVtx = faces[i]->findMateVertexOfEdge( this );
        rg_Point3D  matePt  = mateVtx->getBall().getCenter();
        projectedVtx[i]     = basePlane.projectPointOnPlane( matePt );

        if ( i == 0) {
            angle[i] = 0.0;
        }
        else {
            angle[i] = basePlane.computeAngleInCCW( projectedVtx[0], startPt, projectedVtx[i] );
        }
    }


    //  3. sort beta-faces by angle.
    rg_BOOL swapped        = rg_TRUE;
    rg_INT  timesToCompare = numFaces;
    do {
		swapped = rg_FALSE;
		timesToCompare--;

		for ( i=0; i < timesToCompare; i++) {
			if ( angle[i] > angle[i+1] ) {
                rg_REAL tempAngle = angle[i];
                angle[i]   = angle[i+1];
                angle[i+1] = tempAngle;

                BetaFace* tempFace = faces[i];
                faces[i]   = faces[i+1];
                faces[i+1] = tempFace;

                swapped    = rg_TRUE;
			}
		}            
    } while( swapped );


    //  4. make betaFaceList for sorted face.
    betaFaceList.removeAll();
    for ( i=0; i<numFaces; i++ ) {
        betaFaceList.add( faces[i] );
    }


    delete [] faces;
    delete [] projectedVtx;
    delete [] angle;
}



void BetaEdge::searchCellsInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaCell*>& betaCellList ) const
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





void BetaEdge::searchEdgesInBetaShape( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList ) const
{
    rg_dList<BetaEdge*> finiteEdgeList;
    searchFiniteEdgesInWholeWorld( finiteEdgeList );

    BetaEdge* currEdge = rg_NULL;
    finiteEdgeList.reset4Loop();
    while ( finiteEdgeList.setNext4Loop() ) {
        currEdge = finiteEdgeList.getEntity();

        rg_INT state = currEdge->getBoundingState(beta);
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaEdgeList.add( currEdge );
        }
    }
}



void BetaEdge::searchFacesInBetaShape( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const
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
       

    if ( !isGateToSmallWorlds() || betaFaceList.getSize() <= 2 ) {
        return;
    }


    ///////////////////////////////////////////////////////////////////////////
    //  Some applications need the incident faces to be sorted.
    //
    //  1. generate plane for computing angle of beta-face.
    rg_Point3D startPt = m_startVertex->getBall().getCenter();
    rg_Point3D endPt   = m_endVertex->getBall().getCenter();

    rg_Point3D normal = endPt - startPt;
    normal.normalize();

    Plane basePlane;
    basePlane.definePlaneByNormalAndPassingPoint(normal, startPt);



    //  2. construct angle of beta-face.
    rg_INT numFaces  = betaFaceList.getSize();

    BetaFace**  faces        = betaFaceList.getArray();
    rg_Point3D* projectedVtx = new rg_Point3D[numFaces];
    rg_REAL*    angle        = new rg_REAL[numFaces];
    
    rg_INT i=0;
    for ( i=0; i<numFaces; i++ ) {
        BetaVertex* mateVtx = faces[i]->findMateVertexOfEdge( this );
        rg_Point3D  matePt  = mateVtx->getBall().getCenter();
        projectedVtx[i]     = basePlane.projectPointOnPlane( matePt );

        if ( i == 0) {
            angle[i] = 0.0;
        }
        else {
            angle[i] = basePlane.computeAngleInCCW( projectedVtx[0], startPt, projectedVtx[i] );
        }
    }


    //  3. sort beta-faces by angle.
    rg_BOOL swapped        = rg_TRUE;
    rg_INT  timesToCompare = numFaces;
    do {
		swapped = rg_FALSE;
		timesToCompare--;

		for ( i=0; i < timesToCompare; i++) {
			if ( angle[i] > angle[i+1] ) {
                rg_REAL tempAngle = angle[i];
                angle[i]   = angle[i+1];
                angle[i+1] = tempAngle;

                BetaFace* tempFace = faces[i];
                faces[i]   = faces[i+1];
                faces[i+1] = tempFace;

                swapped    = rg_TRUE;
			}
		}            
    } while( swapped );


    //  4. make betaFaceList for sorted face.
    betaFaceList.removeAll();
    for ( i=0; i<numFaces; i++ ) {
        betaFaceList.add( faces[i] );
    }


    delete [] faces;
    delete [] projectedVtx;
    delete [] angle;
}




//
//  END of Quasi-operator
///////////////////////////////////////////////////////////////////////////



void BetaEdge::connectVFace(VDFace* v_face)
{
    m_vFace = v_face;
}



void BetaEdge::disconnectVFace(VDFace* v_face)
{
    if ( m_vFace == v_face ) {
        m_vFace = rg_NULL;
    }
}

