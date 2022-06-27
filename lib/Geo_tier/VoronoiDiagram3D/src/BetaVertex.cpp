#include "BetaVertex.h"

#include "BetaEdge.h"
#include "BetaFace.h"
#include "BetaCell.h"
#include "VDCell.h"
using namespace V::GeometryTier;

#include <set>
using namespace std;

BetaVertex::BetaVertex()
: m_firstEdge(rg_NULL), m_firstCell(rg_NULL),
  m_visited(rg_FALSE), 
  m_ball(rg_NULL)
{
    m_depth = -1;
    m_vCell = rg_NULL;
}



BetaVertex::BetaVertex(const rg_INT& ID)
: TopologicalEntity(ID),
  m_firstEdge(rg_NULL), m_firstCell(rg_NULL),
  m_visited(rg_FALSE), 
  m_ball(rg_NULL)
{
    m_depth = -1;
    m_vCell = rg_NULL;
}



BetaVertex::BetaVertex(const rg_INT& ID, Ball* ball)
: TopologicalEntity(ID),
  m_firstEdge(rg_NULL), m_firstCell(rg_NULL),
  m_visited(rg_FALSE), 
  m_ball(ball)
{
    m_depth = -1;
    m_vCell = rg_NULL;
}



BetaVertex::BetaVertex(const BetaVertex& betaVtx)
: TopologicalEntity(betaVtx.m_ID),
  m_firstEdge(betaVtx.m_firstEdge), m_firstCell(betaVtx.m_firstCell),
  m_betaSpan(betaVtx.m_betaSpan),
  m_visited(betaVtx.m_visited), m_ball(betaVtx.m_ball)
{
    m_depth = betaVtx.m_depth;
    m_vCell = betaVtx.m_vCell;
}



BetaVertex::~BetaVertex()
{
    m_firstEdge = rg_NULL;
    m_firstCell = rg_NULL;
    m_ball      = rg_NULL;

    //if ( m_vCell != rg_NULL ) {
    //    m_vCell->disconnectBetaVertex(this);
    //}
}



rg_BOOL BetaVertex::isVirtual() const
{
    if ( m_ball == rg_NULL ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}


BetaEdge* BetaVertex::getFirstEdge() const
{
    return m_firstEdge;
}



BetaCell* BetaVertex::getFirstCell() const
{
    return m_firstCell;
}



rg_BetaSpan  BetaVertex::getBetaSpan() const
{
    return m_betaSpan;
}



rg_REAL BetaVertex::getMinValueForValidSimplex() const
{
    return m_betaSpan.getLowerValueOfInterval(1);
}


Ball*     BetaVertex::getBallProperty() const
{
    return m_ball;
}


Sphere         BetaVertex::getBall() const
{
    if ( m_ball != rg_NULL ) {
        return m_ball->getGeometry();
    }
    else {
        Sphere sphere;
        sphere.setRadius(REAL_MINUS_INFINITE);

        return sphere;
    }
}



void*          BetaVertex::getProperty()
{
    if ( m_ball != rg_NULL ) {
        return m_ball->getProperty();
    }
    else {
        return rg_NULL;
    }
}




void BetaVertex::setFirstEdge(BetaEdge* firstEdge)
{
    m_firstEdge = firstEdge;
}



void BetaVertex::setFirstCell(BetaCell* firstCell)
{
    m_firstCell = firstCell;
}



void BetaVertex::setBetaSpan(const rg_BetaSpan& betaSpan)
{
    m_betaSpan = betaSpan;
}



void BetaVertex::setBall(const Sphere& ball)
{
    if ( m_ball != rg_NULL ) {
        m_ball->setGeometry(ball);
    }
}



void BetaVertex::setProperty(void* property)
{
    if ( m_ball != rg_NULL ) {
        m_ball->setProperty(property);
    }
}



void BetaVertex::setBallProperty(Ball* ball)
{
    m_ball = ball;
}

void BetaVertex::shiftBetaSpan(const rg_REAL& delta)
{
	m_betaSpan.shiftBetaSpan(delta);
}

BetaVertex& BetaVertex::operator =(const BetaVertex& betaVtx)
{
    if ( this == &betaVtx ) {
        return *this;
    }

    m_ID        = betaVtx.m_ID;

    m_firstEdge = betaVtx.m_firstEdge;
    m_firstCell = betaVtx.m_firstCell;
    m_betaSpan  = betaVtx.m_betaSpan;
    m_depth     = betaVtx.m_depth;

    m_visited   = betaVtx.m_visited;
    m_ball      = betaVtx.m_ball;
        
    m_vCell = betaVtx.m_vCell;

    return *this;
}



rg_BOOL BetaVertex::isOnConvexHull() const
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




void BetaVertex::computeBetaSpan()
{
    rg_dList<BetaEdge*> incidentEdgeInWholeWorld;
    searchFiniteEdgesInWholeWorld( incidentEdgeInWholeWorld );

    rg_REAL minOfHigherSimplex = REAL_PLUS_INFINITE;
    rg_REAL maxOfHigherSimplex = REAL_MINUS_INFINITE;

    BetaEdge* currEdge = rg_NULL;
    incidentEdgeInWholeWorld.reset4Loop();
    while ( incidentEdgeInWholeWorld.setNext4Loop() ) {
        currEdge = incidentEdgeInWholeWorld.getEntity();

        Interval_REAL interval = currEdge->getBetaSpan().getBetaIntervalOfSingularRegularState();

        rg_REAL minValue = interval.getLowerValue();
        rg_REAL maxValue = interval.getUpperValue();
        if ( minValue < minOfHigherSimplex ) {
            minOfHigherSimplex = minValue;
        }

        if ( maxValue > maxOfHigherSimplex ) {
            maxOfHigherSimplex = maxValue;
        }                
    }


    if ( isOnConvexHull() ) {
        m_betaSpan.setNumBetaInterVal(2);
        m_betaSpan.setBetaInterval(0, SINGULAR_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
        m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,  minOfHigherSimplex,  REAL_PLUS_INFINITE);
    }
    else {
        m_betaSpan.setNumBetaInterVal(3);
        m_betaSpan.setBetaInterval(0, SINGULAR_SIMPLEX, REAL_MINUS_INFINITE, minOfHigherSimplex);
        m_betaSpan.setBetaInterval(1, REGULAR_SIMPLEX,  minOfHigherSimplex,  maxOfHigherSimplex);
        m_betaSpan.setBetaInterval(2, INTERIOR_SIMPLEX, maxOfHigherSimplex,  REAL_PLUS_INFINITE);
    }
}



///////////////////////////////////////////////////////////////////////////
//  Quasi-operator
//  
//  1. Q(v, Y | W): Query primitives for intra-world.
void BetaVertex::searchVerticesInIntraWorld( rg_dList<BetaVertex*>& vertexList, BetaFace* world ) const
{
    rg_dList<BetaEdge*> edgeList;
    searchEdgesInIntraWorld( edgeList, world );

    BetaEdge* currEdge = rg_NULL;
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        currEdge = edgeList.getEntity();

        if ( this == currEdge->getStartVertex() ) {
            vertexList.add( currEdge->getEndVertex() );
        }
        else {
            vertexList.add( currEdge->getStartVertex() );
        }
    }
}



void BetaVertex::searchEdgesInIntraWorld(    rg_dList<BetaEdge*>&   edgeList,   BetaFace* world ) const
{
    BetaFace* currWorld   = rg_NULL;
    // if world is null, find the biggest world of the vertex.
    if ( world == rg_NULL ) {
        currWorld = m_firstEdge->getFirstFace();
    }
    else {
        if ( world->isThere( this ) ) {
            currWorld = world;
        }
        else {
            return;
        }
    }

    // The current world of the vertex is an isolated face.
    if ( currWorld->isIsolatedFace() ) {
        BetaEdge* incidentEdge[2];
        currWorld->findEdgeIncidentToVertex(this, incidentEdge); 
        for ( rg_INT i=0; i<2; i++ ) {
            if ( !incidentEdge[i]->isVisited() ) {
                edgeList.add( incidentEdge[i] );
            }
        }
    }
    else {
        BetaCell* cellInWorld = currWorld->getRightCell();;


        rg_dList<BetaCell*> stackOfCellsToBeVisited;
        stackOfCellsToBeVisited.add( cellInWorld );

        rg_dList<BetaCell*> visitedCellList;


        while ( stackOfCellsToBeVisited.getSize() > 0 ) {
            BetaCell* currCell = stackOfCellsToBeVisited.popFront();

            if ( currCell->isVisited() ) {
                continue;
            }

            BetaFace**   faceOnCell   = currCell->getFaces();
            BetaVertex** vertexOnCell = currCell->getVertices();
            for ( rg_INT i_vtx=0; i_vtx<EIWDS_NUM_VERTEX_ON_CELL; i_vtx++ ) {
                if ( this == vertexOnCell[i_vtx] ) {
                    continue;
                }
        
                BetaCell*   neighbor = currCell->getMateCell( i_vtx );
                if ( neighbor->isVisited() ) {
                    continue;
                }
                else {
                    stackOfCellsToBeVisited.add( neighbor );
                }

                BetaEdge** edgeOnFace = faceOnCell[i_vtx]->getEdges();
                for ( rg_INT i_edge=0; i_edge<EIWDS_NUM_EDGE_ON_FACE; i_edge++ ) {
                    if ( edgeOnFace[i_edge]->isVisited() ) {
                        continue;
                    }
                    else {
                        if ( edgeOnFace[i_edge]->isThere( this ) ) {
                            edgeOnFace[i_edge]->isVisited( rg_TRUE );
                            edgeList.add( edgeOnFace[i_edge] );
                        }
                    }
                }
            }

            currCell->isVisited( rg_TRUE );
            visitedCellList.add( currCell );
        }

        visitedCellList.reset4Loop();
        while ( visitedCellList.setNext4Loop() ) {
            visitedCellList.getEntity()->isVisited( rg_FALSE );
        }

        edgeList.reset4Loop();
        while ( edgeList.setNext4Loop() ) {
            edgeList.getEntity()->isVisited( rg_FALSE );
        }
    }
}



void BetaVertex::searchFacesInIntraWorld(    rg_dList<BetaFace*>&   faceList,   BetaFace* world ) const
{
    BetaFace* currWorld   = rg_NULL;
    // if world is null, find the biggest world of the vertex.
    if ( world == rg_NULL ) {
        currWorld = m_firstEdge->getFirstFace();
    }
    else {
        if ( world->isThere( this ) ) {
            currWorld = world;
        }
        else {
            return;
        }
    }

    
    // The current world of the vertex is an isolated face.
    if ( currWorld->isIsolatedFace() ) {
        faceList.add( currWorld );
    }
    else {
        BetaCell* cellInWorld = currWorld->getRightCell();

        rg_dList<BetaCell*> stackOfCellsToBeVisited;
        stackOfCellsToBeVisited.add( cellInWorld );

        rg_dList<BetaCell*> visitedCellList;


        while ( stackOfCellsToBeVisited.getSize() > 0 ) {
            BetaCell* currCell = stackOfCellsToBeVisited.popFront();

            if ( currCell->isVisited() ) {
                continue;
            }


            BetaFace**   faceOnCell   = currCell->getFaces();
            BetaVertex** vertexOnCell = currCell->getVertices();
            for ( rg_INT i_vtx=0; i_vtx<EIWDS_NUM_VERTEX_ON_CELL; i_vtx++ ) {
                if ( this == vertexOnCell[i_vtx] ) {
                    continue;
                }

                BetaCell*   neighbor = currCell->getMateCell( i_vtx );

                if ( !neighbor->isVisited() ) {
                    faceList.add( faceOnCell[i_vtx] );
                    stackOfCellsToBeVisited.add( neighbor );
                }
            }
        
            currCell->isVisited( rg_TRUE );
            visitedCellList.add( currCell );
        }

        visitedCellList.reset4Loop();
        while ( visitedCellList.setNext4Loop() ) {
            visitedCellList.getEntity()->isVisited( rg_FALSE );
        }
    }
}



void BetaVertex::searchCellsInIntraWorld(    rg_dList<BetaCell*>&   cellList,   BetaFace* world ) const
{
    BetaFace* currWorld   = rg_NULL;
    // if world is null, find the biggest world of the vertex.
    if ( world == rg_NULL ) {
        currWorld = m_firstEdge->getFirstFace();
    }
    else {
        if ( world->isThere( this ) ) {
            currWorld = world;
        }
        else {
            return;
        }
    }

    
    // The current world of the vertex is an isolated face.
    if ( currWorld->isIsolatedFace() ) {
        return;
    }
    else {
        BetaCell* cellInWorld = currWorld->getRightCell();


        rg_dList<BetaCell*> stackOfCellsToBeVisited;
        cellInWorld->isVisited(rg_TRUE);
        stackOfCellsToBeVisited.add( cellInWorld );

        while ( stackOfCellsToBeVisited.getSize() > 0 ) {
            BetaCell* currCell = stackOfCellsToBeVisited.popFront();

            cellList.add( currCell );

            //BetaFace**   faceOnCell   = currCell->getFaces();
            BetaVertex** vertexOnCell = currCell->getVertices();
            for ( rg_INT i_vtx=0; i_vtx<EIWDS_NUM_VERTEX_ON_CELL; i_vtx++ ) {
                if ( this == vertexOnCell[i_vtx] ) {
                    continue;
                }

                BetaCell*   neighbor = currCell->getMateCell( i_vtx );

                if ( !neighbor->isVisited() ) {
                    neighbor->isVisited(rg_TRUE);
                    stackOfCellsToBeVisited.add( neighbor );
                }
            }
        }


        cellList.reset4Loop();
        while ( cellList.setNext4Loop() ) {
            cellList.getEntity()->isVisited( rg_FALSE );
        }
    }
}



void BetaVertex::searchFiniteVerticesInIntraWorld(    rg_dList<BetaVertex*>&   finiteVertexList,   BetaFace* world ) const
{
    rg_dList<BetaVertex*> vertexList;
    searchVerticesInIntraWorld( vertexList, world );

    BetaVertex* currVertex = rg_NULL;
    vertexList.reset4Loop();
    while ( vertexList.setNext4Loop() ) {
        currVertex = vertexList.getEntity();

        if ( currVertex->isVirtual() == rg_FALSE ) {
            finiteVertexList.add( currVertex );
        }
    }
}



void BetaVertex::searchFiniteEdgesInIntraWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world ) const
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



void BetaVertex::searchFiniteFacesInIntraWorld(    rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world ) const
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



void BetaVertex::searchFiniteCellsInIntraWorld(    rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world ) const
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



///////////////////////////////////////////////////////////////////////////};
//
//  2. Query primitives for inter-world.
void BetaVertex::searchSmallWorlds(      rg_dList<BetaFace*>& smallWorldList, BetaFace* world ) const
{
    BetaFace* currWorld = rg_NULL;
    // if world is null, find the biggest world of the vertex.
    if ( world == rg_NULL ) {
        BetaFace* faceOnBiggestWorld = m_firstEdge->getFirstFace();
        currWorld = faceOnBiggestWorld;
    }
    else {
        currWorld = world;
    }
    
    
    rg_dList<BetaFace*> facesOnCurrWorld;
    searchFiniteFacesInIntraWorld( facesOnCurrWorld, currWorld );

    facesOnCurrWorld.reset4Loop();
    while ( facesOnCurrWorld.setNext4Loop() ) {
        facesOnCurrWorld.getEntity()->isVisited( rg_TRUE );
    }

    rg_dList<BetaEdge*> edgeList;
    searchFiniteEdgesInIntraWorld( edgeList, currWorld );

    BetaEdge* currEdge = rg_NULL;
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        currEdge = edgeList.getEntity();

        if ( currEdge->getFirstFace()->isVisited() ) {
            rg_dList<BetaFace*>* smallWorlds = currEdge->getSmallWorlds();
            smallWorldList.append( *smallWorlds );
        }
    }

    facesOnCurrWorld.reset4Loop();
    while ( facesOnCurrWorld.setNext4Loop() ) {
        facesOnCurrWorld.getEntity()->isVisited( rg_FALSE );
    }
}



void BetaVertex::searchWholeSmallWorlds( rg_dList<BetaFace*>& smallWorldList, BetaFace* world ) const
{
    BetaFace* currWorld = rg_NULL;
    // if world is null, find the biggest world of the vertex.
    if ( world == rg_NULL ) {
        BetaFace* faceOnBiggestWorld = m_firstEdge->getFirstFace();
        currWorld = faceOnBiggestWorld;
    }
    else {
        currWorld = world;
    }
    

    rg_dList<BetaFace*> stackForWorlds;
    stackForWorlds.add( currWorld );

    while ( stackForWorlds.getSize() > 0 ) {
        BetaFace* searchingWorld = stackForWorlds.popFront();

        rg_dList<BetaFace*> facesOnCurrWorld;
        searchFiniteFacesInIntraWorld( facesOnCurrWorld, searchingWorld );

        facesOnCurrWorld.reset4Loop();
        while ( facesOnCurrWorld.setNext4Loop() ) {
            facesOnCurrWorld.getEntity()->isVisited( rg_TRUE );
        }

        rg_dList<BetaEdge*> edgeList;
        searchFiniteEdgesInIntraWorld( edgeList, searchingWorld );

        BetaEdge* currEdge = rg_NULL;
        edgeList.reset4Loop();
        while ( edgeList.setNext4Loop() ) {
            currEdge = edgeList.getEntity();

            if ( currEdge->getFirstFace()->isVisited() ) {
                if ( currEdge->isGateToSmallWorlds() ) {
                    rg_dList<BetaFace*>* smallWorlds = currEdge->getSmallWorlds();
                    smallWorldList.append( *smallWorlds );
                    stackForWorlds.append( *smallWorlds );
                }
            }
        }

        facesOnCurrWorld.reset4Loop();
        while ( facesOnCurrWorld.setNext4Loop() ) {
            facesOnCurrWorld.getEntity()->isVisited( rg_FALSE );
        }    
    }
}



void BetaVertex::searchVerticesInSmallWorld( rg_dList<BetaVertex*>& vertexList, BetaFace* world ) const
{
    rg_dList<BetaFace*> smallWorldList;
    searchSmallWorlds( smallWorldList, world );

    BetaFace* currWorld = rg_NULL;
    smallWorldList.reset4Loop();
    while ( smallWorldList.setNext4Loop() ) {
        currWorld = smallWorldList.getEntity();

        searchVerticesInIntraWorld( vertexList, currWorld );
    }
}



void BetaVertex::searchEdgesInSmallWorld(    rg_dList<BetaEdge*>&   edgeList,   BetaFace* world ) const
{
    rg_dList<BetaFace*> smallWorldList;
    searchSmallWorlds( smallWorldList, world );

    BetaFace* currWorld = rg_NULL;
    smallWorldList.reset4Loop();
    while ( smallWorldList.setNext4Loop() ) {
        currWorld = smallWorldList.getEntity();

        searchEdgesInIntraWorld( edgeList, currWorld );

        edgeList.reset4Loop();
        while ( edgeList.setNext4Loop() ) {
            edgeList.getEntity()->isVisited( rg_TRUE );
        }
    }
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        edgeList.getEntity()->isVisited( rg_FALSE );
    }
}



void BetaVertex::searchFacesInSmallWorld(    rg_dList<BetaFace*>&   faceList,   BetaFace* world ) const
{
    rg_dList<BetaFace*> smallWorldList;
    searchSmallWorlds( smallWorldList, world );

    BetaFace* currWorld = rg_NULL;
    smallWorldList.reset4Loop();
    while ( smallWorldList.setNext4Loop() ) {
        currWorld = smallWorldList.getEntity();

        searchFacesInIntraWorld( faceList, currWorld );
    }
}



void BetaVertex::searchCellsInSmallWorld(    rg_dList<BetaCell*>&   cellList,   BetaFace* world ) const
{
    rg_dList<BetaFace*> smallWorldList;
    searchSmallWorlds( smallWorldList, world );

    BetaFace* currWorld = rg_NULL;
    smallWorldList.reset4Loop();
    while ( smallWorldList.setNext4Loop() ) {
        currWorld = smallWorldList.getEntity();

        searchCellsInIntraWorld( cellList, currWorld );
    }
}



void BetaVertex::searchFiniteVerticesInSmallWorld( rg_dList<BetaVertex*>& finiteVertexList, BetaFace* world ) const
{
    rg_dList<BetaFace*> smallWorldList;
    searchSmallWorlds( smallWorldList, world );

    BetaFace* currWorld = rg_NULL;
    smallWorldList.reset4Loop();
    while ( smallWorldList.setNext4Loop() ) {
        currWorld = smallWorldList.getEntity();

        searchFiniteVerticesInIntraWorld( finiteVertexList, currWorld );
    }
}



void BetaVertex::searchFiniteEdgesInSmallWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world ) const
{
    rg_dList<BetaFace*> smallWorldList;
    searchSmallWorlds( smallWorldList, world );

    BetaFace* currWorld = rg_NULL;
    smallWorldList.reset4Loop();
    while ( smallWorldList.setNext4Loop() ) {
        currWorld = smallWorldList.getEntity();

        searchFiniteEdgesInIntraWorld( finiteEdgeList, currWorld );
    }
}



void BetaVertex::searchFiniteFacesInSmallWorld(    rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world ) const
{
    rg_dList<BetaFace*> smallWorldList;
    searchSmallWorlds( smallWorldList, world );

    BetaFace* currWorld = rg_NULL;
    smallWorldList.reset4Loop();
    while ( smallWorldList.setNext4Loop() ) {
        currWorld = smallWorldList.getEntity();

        searchFiniteFacesInIntraWorld( finiteFaceList, currWorld );
    }
}



void BetaVertex::searchFiniteCellsInSmallWorld(    rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world ) const
{
    rg_dList<BetaFace*> smallWorldList;
    searchSmallWorlds( smallWorldList, world );

    BetaFace* currWorld = rg_NULL;
    smallWorldList.reset4Loop();
    while ( smallWorldList.setNext4Loop() ) {
        currWorld = smallWorldList.getEntity();

        searchFiniteCellsInIntraWorld( finiteCellList, currWorld );
    }
}





///////////////////////////////////////////////////////////////////////////};
//
//  3. Q(f, Y): Query primitives for whole-world.
void BetaVertex::searchVerticesInWholeWorld( rg_dList<BetaVertex*>& vertexList ) const
{
    rg_dList<BetaEdge*> edgeList;
    searchEdgesInWholeWorld( edgeList );

    set<BetaVertex*> adjacentVertices;
    BetaEdge* currEdge = rg_NULL;
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        currEdge = edgeList.getEntity();

        if ( this == currEdge->getStartVertex() ) {
            adjacentVertices.insert( currEdge->getEndVertex() );
        }
        else {
            adjacentVertices.insert( currEdge->getStartVertex() );
        }
    }

    set<BetaVertex*>::iterator i_vtx;
    for ( i_vtx=adjacentVertices.begin(); i_vtx!=adjacentVertices.end(); i_vtx++ ) {
        BetaVertex* currVtx = *i_vtx;
        vertexList.add( currVtx );
    }

    //BetaEdge* currEdge = rg_NULL;
    //edgeList.reset4Loop();
    //while ( edgeList.setNext4Loop() ) {
    //    currEdge = edgeList.getEntity();

    //    if ( this == currEdge->getStartVertex() ) {
    //        vertexList.add( currEdge->getEndVertex() );
    //    }
    //    else {
    //        vertexList.add( currEdge->getStartVertex() );
    //    }
    //}
}



void BetaVertex::searchEdgesInWholeWorld(    rg_dList<BetaEdge*>&   edgeList   ) const
{
    rg_dList<BetaFace*> contributingWorlds;
    contributingWorlds.add( rg_NULL );

    searchWholeSmallWorlds( contributingWorlds );

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        rg_dList<BetaEdge*> edgesInCurrWorld;
        searchEdgesInIntraWorld( edgesInCurrWorld, currWorld);

        edgesInCurrWorld.reset4Loop();
        while ( edgesInCurrWorld.setNext4Loop() ) {
            edgesInCurrWorld.getEntity()->isVisited( rg_TRUE );
        }

        edgeList.mergeTail( edgesInCurrWorld );
    }

    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        edgeList.getEntity()->isVisited( rg_FALSE );
    }
}



void BetaVertex::searchFacesInWholeWorld(    rg_dList<BetaFace*>&   faceList   ) const
{
    searchFacesInIntraWorld( faceList );

    rg_dList<BetaFace*> contributingWorlds;
    searchWholeSmallWorlds( contributingWorlds );

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        searchFacesInIntraWorld( faceList, currWorld);
    }
}



void BetaVertex::searchCellsInWholeWorld(    rg_dList<BetaCell*>&   cellList   ) const
{
    searchCellsInIntraWorld( cellList );

    rg_dList<BetaFace*> contributingWorlds;
    searchWholeSmallWorlds( contributingWorlds );

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        searchCellsInIntraWorld( cellList, currWorld);
    }
}


void BetaVertex::searchFiniteVerticesInWholeWorld( rg_dList<BetaVertex*>& finiteVertexList ) const
{
    rg_dList<BetaVertex*> vertexList;
    searchVerticesInWholeWorld( vertexList );

    BetaVertex* currVertex = rg_NULL;
    vertexList.reset4Loop();
    while ( vertexList.setNext4Loop() ) {
        currVertex = vertexList.getEntity();

        if ( currVertex->isVirtual() == rg_FALSE ) {
            finiteVertexList.add( currVertex );
        }
    }
}



void BetaVertex::searchFiniteEdgesInWholeWorld(    rg_dList<BetaEdge*>&   finiteEdgeList   ) const
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



void BetaVertex::searchFiniteFacesInWholeWorld(    rg_dList<BetaFace*>&   finiteFaceList   ) const
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



void BetaVertex::searchFiniteCellsInWholeWorld(    rg_dList<BetaCell*>&   finiteCellList   ) const
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



void BetaVertex::searchFiniteVerticesInSmallerWorldViaGateEdge(BetaEdge* gateEdge, rg_dList<BetaVertex*>& finiteVertexList   ) const
{
    if ( !gateEdge->isGateToSmallWorlds() ) {
        return;
    }


    rg_dList<BetaFace*>* smallWorlds = gateEdge->getSmallWorlds();

    rg_dList<BetaFace*> contributingWorlds;
    contributingWorlds.appendTail( *smallWorlds );

    smallWorlds->reset4Loop();
    while ( smallWorlds->setNext4Loop() ) {
        BetaFace* world = smallWorlds->getEntity();
        searchWholeSmallWorlds( contributingWorlds, world );
    }

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        rg_dList<BetaVertex*> verticesInCurrWorld;
        searchVerticesInIntraWorld( verticesInCurrWorld, currWorld);

        verticesInCurrWorld.reset4Loop();
        while ( verticesInCurrWorld.setNext4Loop() ) {
            verticesInCurrWorld.getEntity()->isVisited( rg_TRUE );
        }

        finiteVertexList.mergeTail( verticesInCurrWorld );
    }

    finiteVertexList.reset4Loop();
    while ( finiteVertexList.setNext4Loop() ) {
        finiteVertexList.getEntity()->isVisited( rg_FALSE );
    }
}



void BetaVertex::searchFiniteEdgesInSmallerWorldViaGateEdge(   BetaEdge* gateEdge, rg_dList<BetaEdge*>&   finiteEdgeList   ) const
{
    if ( !gateEdge->isGateToSmallWorlds() ) {
        return;
    }


    rg_dList<BetaFace*>* smallWorlds = gateEdge->getSmallWorlds();

    rg_dList<BetaFace*> contributingWorlds;
    contributingWorlds.appendTail( *smallWorlds );

    smallWorlds->reset4Loop();
    while ( smallWorlds->setNext4Loop() ) {
        BetaFace* world = smallWorlds->getEntity();
        searchWholeSmallWorlds( contributingWorlds, world );
    }

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        rg_dList<BetaEdge*> edgesInCurrWorld;
        searchEdgesInIntraWorld( edgesInCurrWorld, currWorld);

        edgesInCurrWorld.reset4Loop();
        while ( edgesInCurrWorld.setNext4Loop() ) {
            edgesInCurrWorld.getEntity()->isVisited( rg_TRUE );
        }

        finiteEdgeList.mergeTail( edgesInCurrWorld );
    }

    finiteEdgeList.reset4Loop();
    while ( finiteEdgeList.setNext4Loop() ) {
        finiteEdgeList.getEntity()->isVisited( rg_FALSE );
    }
}



void BetaVertex::searchFiniteFacesInSmallerWorldViaGateEdge(   BetaEdge* gateEdge, rg_dList<BetaFace*>&   finiteFaceList   ) const
{
    if ( !gateEdge->isGateToSmallWorlds() ) {
        return;
    }


    rg_dList<BetaFace*>* smallWorlds = gateEdge->getSmallWorlds();

    rg_dList<BetaFace*> contributingWorlds;
    contributingWorlds.appendTail( *smallWorlds );

    smallWorlds->reset4Loop();
    while ( smallWorlds->setNext4Loop() ) {
        BetaFace* world = smallWorlds->getEntity();
        searchWholeSmallWorlds( contributingWorlds, world );
    }

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        rg_dList<BetaFace*> facesInCurrWorld;
        searchFacesInIntraWorld( facesInCurrWorld, currWorld);

        facesInCurrWorld.reset4Loop();
        while ( facesInCurrWorld.setNext4Loop() ) {
            facesInCurrWorld.getEntity()->isVisited( rg_TRUE );
        }

        finiteFaceList.mergeTail( facesInCurrWorld );
    }

    finiteFaceList.reset4Loop();
    while ( finiteFaceList.setNext4Loop() ) {
        finiteFaceList.getEntity()->isVisited( rg_FALSE );
    }
}



void BetaVertex::searchFiniteCellsInSmallerWorldViaGateEdge(   BetaEdge* gateEdge, rg_dList<BetaCell*>&   finiteCellList   ) const
{
    if ( !gateEdge->isGateToSmallWorlds() ) {
        return;
    }


    rg_dList<BetaFace*>* smallWorlds = gateEdge->getSmallWorlds();

    rg_dList<BetaFace*> contributingWorlds;
    contributingWorlds.appendTail( *smallWorlds );

    smallWorlds->reset4Loop();
    while ( smallWorlds->setNext4Loop() ) {
        BetaFace* world = smallWorlds->getEntity();
        searchWholeSmallWorlds( contributingWorlds, world );
    }

    while ( contributingWorlds.getSize() > 0 ) {
        BetaFace* currWorld = contributingWorlds.popFront();

        rg_dList<BetaCell*> cellsInCurrWorld;
        searchCellsInIntraWorld( cellsInCurrWorld, currWorld);

        cellsInCurrWorld.reset4Loop();
        while ( cellsInCurrWorld.setNext4Loop() ) {
            cellsInCurrWorld.getEntity()->isVisited( rg_TRUE );
        }

        finiteCellList.mergeTail( cellsInCurrWorld );
    }

    finiteCellList.reset4Loop();
    while ( finiteCellList.setNext4Loop() ) {
        finiteCellList.getEntity()->isVisited( rg_FALSE );
    }
}




//
//  END of Quasi-operator
///////////////////////////////////////////////////////////////////////////};



///////////////////////////////////////////////////////////////////////////
//  Beta-operator
//  
void BetaVertex::searchVerticesInBetaComplex(const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const
{
    rg_dList<BetaEdge*> incidentEdgeInBetaComplex;
    searchEdgesInBetaComplex(beta, incidentEdgeInBetaComplex );

    BetaEdge* currEdge = rg_NULL;
    incidentEdgeInBetaComplex.reset4Loop();
    while ( incidentEdgeInBetaComplex.setNext4Loop() ) {
        currEdge = incidentEdgeInBetaComplex.getEntity();

        if ( currEdge->getStartVertex() == this ) {
            betaVertexList.add( currEdge->getEndVertex() );
        }
        else {
            betaVertexList.add( currEdge->getStartVertex() );
        }
    }


    //  The following code has bugs
    //  because it does not check that incident edge is in the beta-complex or not
    //
    //BetaVertex* currVertex = rg_NULL;
    //finiteVertexList.reset4Loop();
    //while ( finiteVertexList.setNext4Loop() ) {
    //    currVertex = finiteVertexList.getEntity();

    //    if ( currVertex->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
    //        continue;
    //    }
    //    else {
    //        betaVertexList.add( currVertex );
    //    }
    //}
}



void BetaVertex::searchEdgesInBetaComplex(   const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList   ) const
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



void BetaVertex::searchFacesInBetaComplex(   const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList   ) const
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



void BetaVertex::searchCellsInBetaComplex(   const rg_REAL& beta, rg_dList<BetaCell*>&   betaCellList   ) const
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



void BetaVertex::searchExtraneousVerticesInBetaComplex(const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const
{
    //  There are no exterior vertics in beta-complex.
    //
    //rg_dList<BetaVertex*> finiteVertexList;
    //searchFiniteVerticesInWholeWorld( finiteVertexList );

    //BetaVertex* currVertex = rg_NULL;
    //finiteVertexList.reset4Loop();
    //while ( finiteVertexList.setNext4Loop() ) {
    //    currVertex = finiteVertexList.getEntity();

    //    if ( currVertex->getBoundingState(beta) == EXTRANEOUS_SIMPLEX ) {
    //        betaVertexList.add( currVertex );
    //    }
    //}
}



void BetaVertex::searchExtraneousEdgesInBetaComplex(   const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList   ) const
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



void BetaVertex::searchExtraneousFacesInBetaComplex(   const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList   ) const
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



void BetaVertex::searchExtraneousCellsInBetaComplex(   const rg_REAL& beta, rg_dList<BetaCell*>&   betaCellList   ) const
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



void BetaVertex::searchVerticesInBetaComplex(const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaVertex*>& betaVertexList ) const
{
    rg_dList<BetaVertex*> neighborVertexInBetaComplex;
    searchVerticesInBetaComplex(beta, neighborVertexInBetaComplex );

    
    BetaVertex* currVertex = rg_NULL;
    neighborVertexInBetaComplex.reset4Loop();
    while ( neighborVertexInBetaComplex.setNext4Loop() ) {
        currVertex = neighborVertexInBetaComplex.getEntity();

        if ( currVertex->getBoundingState(beta) == boundingState ) {
            betaVertexList.add( currVertex );
        }
    }


    //rg_dList<BetaVertex*> finiteVertexList;
    //searchFiniteVerticesInWholeWorld( finiteVertexList );

    //BetaVertex* currVertex = rg_NULL;
    //finiteVertexList.reset4Loop();
    //while ( finiteVertexList.setNext4Loop() ) {
    //    currVertex = finiteVertexList.getEntity();

    //    if ( currVertex->getBoundingState(beta) == boundingState ) {
    //        betaVertexList.add( currVertex );
    //    }
    //}
}


void BetaVertex::searchEdgesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaEdge*>& betaEdgeList ) const
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



void BetaVertex::searchFacesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaFace*>& betaFaceList ) const
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



void BetaVertex::searchCellsInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaCell*>& betaCellList ) const
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




void BetaVertex::searchVerticesInBetaShape( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const
{
    rg_dList<BetaEdge*> incidentEdgeInBetaShape;
    searchEdgesInBetaShape(beta, incidentEdgeInBetaShape );

    BetaEdge* currEdge = rg_NULL;
    incidentEdgeInBetaShape.reset4Loop();
    while ( incidentEdgeInBetaShape.setNext4Loop() ) {
        currEdge = incidentEdgeInBetaShape.getEntity();

        if ( currEdge->getStartVertex() == this ) {
            betaVertexList.add( currEdge->getEndVertex() );
        }
        else {
            betaVertexList.add( currEdge->getStartVertex() );
        }
    }

    //  The following code has bugs
    //  because it does not check that incident edge is in the beta-complex or not
    //
    //rg_dList<BetaVertex*> finiteVertexList;
    //searchFiniteVerticesInWholeWorld( finiteVertexList );

    //BetaVertex* currVertex = rg_NULL;
    //finiteVertexList.reset4Loop();
    //while ( finiteVertexList.setNext4Loop() ) {
    //    currVertex = finiteVertexList.getEntity();

    //    rg_INT state = currVertex->getBoundingState(beta);
    //    if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
    //        betaVertexList.add( currVertex );
    //    }
    //}
}



void BetaVertex::searchEdgesInBetaShape(    const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList   ) const
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



void BetaVertex::searchFacesInBetaShape(    const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList   ) const
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



rg_INT BetaVertex::countNumGroupsOfFaceConnectedExtraneousCellsInBetaComplex(  const rg_REAL& beta ) const
{
    rg_INT numGroups = 0;

    rg_dList<BetaCell*> incidentQuasiCells;
    searchCellsInWholeWorld( incidentQuasiCells );

    BetaCell* currCell = rg_NULL;    
    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        currCell = incidentQuasiCells.getEntity();

        if ( currCell->isVisited() ) {
            continue;
        }

        if ( currCell->getBoundingState(beta) != EXTRANEOUS_SIMPLEX ) {
            continue;
        }

        numGroups++;

        rg_dList<BetaCell*> searchStack;
        searchStack.add( currCell );
        while ( searchStack.getSize() > 0 ) {
            BetaCell* searchingCell = searchStack.popFront();

            if ( searchingCell->isVisited() ) {
                continue;
            }

            searchingCell->isVisited(rg_TRUE);

            BetaFace** boundingFace = searchingCell->getFaces();
            for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
                if ( boundingFace[i]->getBoundingState(beta) != EXTRANEOUS_SIMPLEX ) {
                    continue;
                }

                if ( boundingFace[i]->isThere( this ) ) {
                    BetaCell* neighbor = searchingCell->getNeighborCell(i);
                    searchStack.add( neighbor );
                }
            }
        }
    }


    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        incidentQuasiCells.getEntity()->isVisited(rg_FALSE);
    }


    return numGroups;
}



rg_INT BetaVertex::countNumGroupsOfFaceConnectedInteriorCellsInBetaComplex(  const rg_REAL& beta ) const
{
    rg_INT numGroups = 0;

    
    rg_dList<BetaCell*> incidentQuasiCells;
    searchCellsInWholeWorld( incidentQuasiCells );

    BetaCell* currCell = rg_NULL;    
    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        currCell = incidentQuasiCells.getEntity();

        if ( currCell->isVisited() ) {
            continue;
        }

        if ( currCell->getBoundingState(beta) != INTERIOR_SIMPLEX ) {
            continue;
        }

        numGroups++;

        rg_dList<BetaCell*> searchStack;
        searchStack.add( currCell );
        while ( searchStack.getSize() > 0 ) {
            BetaCell* searchingCell = searchStack.popFront();

            if ( searchingCell->isVisited() ) {
                continue;
            }

            searchingCell->isVisited(rg_TRUE);

            BetaFace** boundingFace = searchingCell->getFaces();
            for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
                if ( boundingFace[i]->getBoundingState(beta) != INTERIOR_SIMPLEX ) {
                    continue;
                }

                if ( boundingFace[i]->isThere( this ) ) {
                    BetaCell* neighbor = searchingCell->getNeighborCell(i);
                    searchStack.add( neighbor );
                }
            }
        }
    }


    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        incidentQuasiCells.getEntity()->isVisited(rg_FALSE);
    }


    return numGroups;
}



void BetaVertex::searchGroupsOfFaceConnectedExtraneousCellsInBetaComplex( const rg_REAL& beta, rg_dList< rg_dList<BetaCell*> >& groupsOfExtraneousCell ) const
{
    rg_dList<BetaCell*> incidentQuasiCells;
    searchCellsInWholeWorld( incidentQuasiCells );

    BetaCell* currCell = rg_NULL;    
    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        currCell = incidentQuasiCells.getEntity();

        if ( currCell->isVisited() ) {
            continue;
        }

        if ( currCell->getBoundingState(beta) != EXTRANEOUS_SIMPLEX ) {
            continue;
        }


        rg_dList<BetaCell*>* currCellGroup = groupsOfExtraneousCell.add( rg_dList<BetaCell*>() );

        rg_dList<BetaCell*> searchStack;
        searchStack.add( currCell );
        while ( searchStack.getSize() > 0 ) {
            BetaCell* searchingCell = searchStack.popFront();

            if ( searchingCell->isVisited() ) {
                continue;
            }

            currCellGroup->add( searchingCell );
            searchingCell->isVisited(rg_TRUE);

            BetaFace** boundingFace = searchingCell->getFaces();
            for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
                if ( boundingFace[i]->getBoundingState(beta) != EXTRANEOUS_SIMPLEX ) {
                    continue;
                }

                if ( boundingFace[i]->isThere( this ) ) {
                    BetaCell* neighbor = searchingCell->getNeighborCell(i);
                    searchStack.add( neighbor );
                }
            }
        }
    }


    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        incidentQuasiCells.getEntity()->isVisited(rg_FALSE);
    }
}



void BetaVertex::searchGroupsOfFaceConnectedInteriorCellsInBetaComplex(   const rg_REAL& beta, rg_dList< rg_dList<BetaCell*> >& groupsOfInteriorCell ) const
{
    rg_dList<BetaCell*> incidentQuasiCells;
    searchCellsInWholeWorld( incidentQuasiCells );

    BetaCell* currCell = rg_NULL;    
    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        currCell = incidentQuasiCells.getEntity();

        if ( currCell->isVisited() ) {
            continue;
        }

        if ( currCell->getBoundingState(beta) != INTERIOR_SIMPLEX ) {
            continue;
        }


        rg_dList<BetaCell*>* currCellGroup = groupsOfInteriorCell.add( rg_dList<BetaCell*>() );

        rg_dList<BetaCell*> searchStack;
        searchStack.add( currCell );
        while ( searchStack.getSize() > 0 ) {
            BetaCell* searchingCell = searchStack.popFront();

            if ( searchingCell->isVisited() ) {
                continue;
            }

            currCellGroup->add( searchingCell );
            searchingCell->isVisited(rg_TRUE);

            BetaFace** boundingFace = searchingCell->getFaces();
            for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
                if ( boundingFace[i]->getBoundingState(beta) != INTERIOR_SIMPLEX ) {
                    continue;
                }

                if ( boundingFace[i]->isThere( this ) ) {
                    BetaCell* neighbor = searchingCell->getNeighborCell(i);
                    searchStack.add( neighbor );
                }
            }
        }
    }


    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        incidentQuasiCells.getEntity()->isVisited(rg_FALSE);
    }
}



void   BetaVertex::searchGroupsOfFaceConnectedExtraneousCellsInGivenWorld(   const rg_REAL& beta, rg_dList< rg_dList<BetaCell*> >& groupsOfExtraneousCell, BetaFace* world ) const
{
    rg_dList<BetaCell*> cellListInGivenWorld;
    searchCellsInIntraWorld( cellListInGivenWorld, world );

    cellListInGivenWorld.reset4Loop();
    while ( cellListInGivenWorld.setNext4Loop() ) {
        BetaCell* currCell = cellListInGivenWorld.getEntity();

        if ( currCell->isVisited() ) {
            continue;
        }

        if ( currCell->getBoundingState(beta) != EXTRANEOUS_SIMPLEX ) {
            continue;
        }


        rg_dList<BetaCell*>* currCellGroup = groupsOfExtraneousCell.add( rg_dList<BetaCell*>() );

        rg_dList<BetaCell*> searchStack;
        searchStack.add( currCell );
        while ( searchStack.getSize() > 0 ) {
            BetaCell* searchingCell = searchStack.popFront();

            if ( searchingCell->isVisited() ) {
                continue;
            }

            currCellGroup->add( searchingCell );
            searchingCell->isVisited(rg_TRUE);

            BetaFace** boundingFace = searchingCell->getFaces();
            for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
                if ( boundingFace[i]->getBoundingState(beta) != EXTRANEOUS_SIMPLEX ) {
                    continue;
                }

                if ( boundingFace[i]->isThere( this ) ) {
                    BetaCell* neighbor = searchingCell->getNeighborCell(i);
                    searchStack.add( neighbor );
                }
            }
        }
    }

    cellListInGivenWorld.reset4Loop();
    while ( cellListInGivenWorld.setNext4Loop() ) {
        cellListInGivenWorld.getEntity()->isVisited(rg_FALSE);
    }
}

//
//  END of Beta-operator
///////////////////////////////////////////////////////////////////////////};




void BetaVertex::connectVCell(VDCell* v_cell)
{
    m_vCell = v_cell;
}



void BetaVertex::disconnectVCell(VDCell* v_cell)
{
    if ( m_vCell == v_cell ) {
        m_vCell = rg_NULL;
    }
}



