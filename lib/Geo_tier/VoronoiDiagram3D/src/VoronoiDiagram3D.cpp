#include "VoronoiDiagram3D.h"

#include "rg_Generator.h"
#include "rg_BallGenerator.h"
#include "FunctionsForVoronoiDiagram3D.h"
#include "rg_RelativeOp.h"
#include <float.h>
using namespace V::GeometryTier;



///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
VoronoiDiagram3D::VoronoiDiagram3D()
{
    m_status          = NORMAL_STATE;
    m_makeVoronoiFace = rg_FALSE;
}



VoronoiDiagram3D::~VoronoiDiagram3D()
{
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
rg_dList< VDCell >* VoronoiDiagram3D::getGlobalCellList() 
{
    return &m_GlobalCellList;
}

rg_dList< VDFace >* VoronoiDiagram3D::getGlobalFaceList() 
{
    return &m_GlobalFaceList;
}

rg_dList< VDLoop >* VoronoiDiagram3D::getGlobalLoopList() 
{
    return &m_GlobalLoopList;
}

rg_dList< VDPartialEdge >* VoronoiDiagram3D::getGlobalPartialEdgeList() 
{
    return &m_GlobalPartialedgeList;
}

rg_dList< VDEdge >* VoronoiDiagram3D::getGlobalEdgeList() 
{
    return &m_GlobalEdgeList;
}

rg_dList< VDVertex >* VoronoiDiagram3D::getGlobalVerticeList() 
{
    return &m_GlobalVertexList;
}


rg_INT VoronoiDiagram3D::getNumOfCells() const
{
    return m_GlobalCellList.getSize();
}

rg_INT VoronoiDiagram3D::getNumOfFaces() const
{
    return m_GlobalFaceList.getSize();
}

rg_INT VoronoiDiagram3D::getNumOfLoops() const
{
    return m_GlobalLoopList.getSize();
}

rg_INT VoronoiDiagram3D::getNumOfPartialEdges() const
{
    return m_GlobalPartialedgeList.getSize();
}

rg_INT VoronoiDiagram3D::getNumOfEdges() const
{
    return m_GlobalEdgeList.getSize();
}

rg_INT VoronoiDiagram3D::getNumOfVertices() const
{
    return m_GlobalVertexList.getSize();
}

rg_INT VoronoiDiagram3D::getNumOfDisconnectedGenerator()
{
    rg_INT numBallInDisconnectedCase = 0;

    VDCell* currCell = rg_NULL;
	m_GlobalCellList.reset4Loop();
	while( m_GlobalCellList.setNext4Loop() )  
    {
		currCell = m_GlobalCellList.getpEntity();
        if ( currCell->getNumOfBoundingFaces() == 0 )
            numBallInDisconnectedCase++;
    }

    return numBallInDisconnectedCase;
}

///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void VoronoiDiagram3D::duplicateTopology(const VoronoiDiagram3D& origin)
{
    //  map< origin_entity*, this_entity*>
    map<VDCell*, VDCell*>               cellMap;
    map<VDFace*, VDFace*>               faceMap;
    map<VDLoop*, VDLoop*>               loopMap;
    map<VDPartialEdge*, VDPartialEdge*> prEdgeMap;
    map<VDEdge*, VDEdge*>               edgeMap;
    map<VDVertex*, VDVertex*>           vertexMap;

    // make binary search tree with origin.entity as key
    makeMapForDuplication( origin, cellMap, faceMap, loopMap, prEdgeMap, edgeMap, vertexMap);

    // duplicate topology of each entity.
    duplicateTopologyOfCells(        origin, cellMap, faceMap);
    duplicateTopologyOfFaces(        origin, cellMap, faceMap, loopMap);
    duplicateTopologyOfLoops(        origin, faceMap, loopMap, prEdgeMap);
    duplicateTopologyOfPartialEdges( origin, loopMap, prEdgeMap, edgeMap);
    duplicateTopologyOfEdges(        origin, prEdgeMap, edgeMap, vertexMap);
    duplicateTopologyOfVertices(     origin,  edgeMap, vertexMap);
}



    
void VoronoiDiagram3D::clean()
{
    m_GlobalCellList.removeAll();
    m_GlobalFaceList.removeAll();
    m_GlobalLoopList.removeAll();
    m_GlobalPartialedgeList.removeAll();
    m_GlobalEdgeList.removeAll();
    m_GlobalVertexList.removeAll();

    m_bucket.removeBucketTable();
}



void VoronoiDiagram3D::cleanForRestart()
{
    m_GlobalCellList.reset4Loop();
    while ( m_GlobalCellList.setNext4Loop() ) {
        VDCell* currCell = m_GlobalCellList.getpEntity();

        currCell->getBoungindFaces()->removeAll();
        currCell->connectQVertex(    rg_NULL );
        currCell->connectBetaVertex( rg_NULL );
    }

    m_GlobalFaceList.removeAll();
    m_GlobalLoopList.removeAll();
    m_GlobalPartialedgeList.removeAll();
    m_GlobalEdgeList.removeAll();
    m_GlobalVertexList.removeAll();

    m_bucket.removeBucketTable();
}



void VoronoiDiagram3D::makeMapForDuplication( const VoronoiDiagram3D& origin, 
                                              map<VDCell*, VDCell*>&               cellMap,
                                              map<VDFace*, VDFace*>&               faceMap,
                                              map<VDLoop*, VDLoop*>&               loopMap,
                                              map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                              map<VDEdge*, VDEdge*>&               edgeMap,
                                              map<VDVertex*, VDVertex*>&           vertexMap)
{
    //  make map for V-cells
    origin.m_GlobalCellList.reset4Loop();
    while ( origin.m_GlobalCellList.setNext4Loop() ) {
        VDCell* currOriCell  = origin.m_GlobalCellList.getpEntity();
        VDCell* currThisCell = m_GlobalCellList.add( VDCell( currOriCell->getID() ) );

        cellMap.insert( make_pair( currOriCell, currThisCell ) );
    }


    //  make map for V-faces
    origin.m_GlobalFaceList.reset4Loop();
    while ( origin.m_GlobalFaceList.setNext4Loop() ) {
        VDFace* currOriFace  = origin.m_GlobalFaceList.getpEntity();
        VDFace* currThisFace = m_GlobalFaceList.add( VDFace( currOriFace->getID() ) );

        faceMap.insert( make_pair( currOriFace, currThisFace ) );
    }


    //  make map for V-loops
    origin.m_GlobalLoopList.reset4Loop();
    while ( origin.m_GlobalLoopList.setNext4Loop() ) {
        VDLoop* currOriLoop  = origin.m_GlobalLoopList.getpEntity();
        VDLoop* currThisLoop = m_GlobalLoopList.add( VDLoop( currOriLoop->getID(), rg_NULL ) );

        loopMap.insert( make_pair( currOriLoop, currThisLoop ) );
    }


    //  make map for V-partial_edges
    origin.m_GlobalPartialedgeList.reset4Loop();
    while ( origin.m_GlobalPartialedgeList.setNext4Loop() ) {
        VDPartialEdge* currOriPrEdge  = origin.m_GlobalPartialedgeList.getpEntity();
        VDPartialEdge* currThisPrEdge = m_GlobalPartialedgeList.add( VDPartialEdge( currOriPrEdge->getID(), rg_NULL ) );

        prEdgeMap.insert( make_pair( currOriPrEdge, currThisPrEdge ) );
    }


    //  make map for V-edges
    origin.m_GlobalEdgeList.reset4Loop();
    while ( origin.m_GlobalEdgeList.setNext4Loop() ) {
        VDEdge* currOriEdge  = origin.m_GlobalEdgeList.getpEntity();
        VDEdge* currThisEdge = m_GlobalEdgeList.add( VDEdge( currOriEdge->getID(), rg_NULL, rg_NULL ) );

        edgeMap.insert( make_pair( currOriEdge, currThisEdge ) );
    }


    //  make map for V-vertices
    origin.m_GlobalVertexList.reset4Loop();
    while ( origin.m_GlobalVertexList.setNext4Loop() ) {
        VDVertex* currOriVtx  = origin.m_GlobalVertexList.getpEntity();
        VDVertex* currThisVtx = m_GlobalVertexList.add( VDVertex( currOriVtx->getID() ) );

        vertexMap.insert( make_pair( currOriVtx, currThisVtx ) );
    }
}



void VoronoiDiagram3D::duplicateTopologyOfCells( const VoronoiDiagram3D& origin, 
                                   map<VDCell*, VDCell*>&               cellMap,
                                   map<VDFace*, VDFace*>&               faceMap)
{
    map<VDFace*, VDFace*>::iterator it_face;

    map<VDFace*, VDFace*>::iterator face_end   = faceMap.end();

    origin.m_GlobalCellList.reset4Loop();
    while ( origin.m_GlobalCellList.setNext4Loop() ) {
        VDCell* currOriCell  = origin.m_GlobalCellList.getpEntity();
        VDCell* currThisCell = cellMap.find( currOriCell )->second;

        rg_dList<VDFace*>* originBoundingFaces = currOriCell->getBoungindFaces();
        originBoundingFaces->reset4Loop();
        while ( originBoundingFaces->setNext4Loop() ) {
            VDFace* currOriFace  = originBoundingFaces->getEntity();

            it_face = faceMap.find( currOriFace );
            if ( it_face != face_end ) {
                currThisCell->addBoundingFace( it_face->second );
            }
        }

        currThisCell->connectBetaVertex( currOriCell->getBetaVertex() );
        currThisCell->connectQVertex(    currOriCell->getQVertex() );
    }
}



void VoronoiDiagram3D::duplicateTopologyOfFaces( const VoronoiDiagram3D& origin, 
                                   map<VDCell*, VDCell*>&               cellMap,
                                   map<VDFace*, VDFace*>&               faceMap,
                                   map<VDLoop*, VDLoop*>&               loopMap)
{
    map<VDCell*, VDCell*>::iterator it_cell;
    map<VDFace*, VDFace*>::iterator it_face;
    map<VDLoop*, VDLoop*>::iterator it_loop;

    map<VDCell*, VDCell*>::iterator cell_end   = cellMap.end();
    map<VDFace*, VDFace*>::iterator face_end   = faceMap.end();
    map<VDLoop*, VDLoop*>::iterator loop_end   = loopMap.end();

    origin.m_GlobalFaceList.reset4Loop();
    while ( origin.m_GlobalFaceList.setNext4Loop() ) {
        VDFace* currOriFace  = origin.m_GlobalFaceList.getpEntity();
        VDFace* currThisFace = faceMap.find( currOriFace )->second;

        VDCell* leftCell  = rg_NULL;
        it_cell = cellMap.find( currOriFace->getLeftCell() );
        if ( it_cell != cell_end ) {
            leftCell = it_cell->second;
        }

        VDCell* rightCell = rg_NULL;
        it_cell = cellMap.find( currOriFace->getRightCell() );
        if ( it_cell != cell_end ) {
            rightCell = it_cell->second;
        }

        currThisFace->setLeftCell(  leftCell );
        currThisFace->setRightCell( rightCell );

        rg_dList<VDLoop*>* originLoops = currOriFace->getLoops();
        originLoops->reset4Loop();
        while ( originLoops->setNext4Loop() ) {
            VDLoop* currOriLoop = originLoops->getEntity();
            VDLoop* currThisLoop = rg_NULL;
            it_loop = loopMap.find( currOriLoop );
            if ( it_loop != loop_end ) {
                currThisLoop = it_loop->second;
                currThisFace->addLoop( currThisLoop );
            }
        }

        currThisFace->connectBetaEdge( currOriFace->getBetaEdge() );
    }
}



void VoronoiDiagram3D::duplicateTopologyOfLoops( const VoronoiDiagram3D& origin, 
                                   map<VDFace*, VDFace*>&               faceMap,
                                   map<VDLoop*, VDLoop*>&               loopMap,
                                   map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap)
{
    map<VDFace*, VDFace*>::iterator               it_face;
    map<VDPartialEdge*, VDPartialEdge*>::iterator it_predge;

    map<VDFace*, VDFace*>::iterator               face_end   = faceMap.end();
    map<VDPartialEdge*, VDPartialEdge*>::iterator predge_end = prEdgeMap.end();

    origin.m_GlobalLoopList.reset4Loop();
    while ( origin.m_GlobalLoopList.setNext4Loop() ) {
        VDLoop* currOriLoop  = origin.m_GlobalLoopList.getpEntity();
        VDLoop* currThisLoop = loopMap.find( currOriLoop )->second;

        VDFace*        currFace   = rg_NULL;
        it_face = faceMap.find( currOriLoop->getFace() );
        if ( it_face != face_end ) {
            currFace = it_face->second;
        }

        VDPartialEdge* currPrEdge = rg_NULL;
        it_predge = prEdgeMap.find( currOriLoop->getPartialEdge() );
        if ( it_predge != predge_end ) {
            currPrEdge = it_predge->second;
        }

        currThisLoop->setLoop( currFace, currPrEdge, currOriLoop->isOuterLoop() );
    }
}



void VoronoiDiagram3D::duplicateTopologyOfPartialEdges( const VoronoiDiagram3D& origin, 
                                   map<VDLoop*, VDLoop*>&               loopMap,
                                   map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                   map<VDEdge*, VDEdge*>&               edgeMap)
{
    map<VDLoop*, VDLoop*>::iterator               it_loop;
    map<VDPartialEdge*, VDPartialEdge*>::iterator it_predge;
    map<VDEdge*, VDEdge*>::iterator               it_edge;

    map<VDLoop*, VDLoop*>::iterator               loop_end   = loopMap.end();
    map<VDPartialEdge*, VDPartialEdge*>::iterator predge_end = prEdgeMap.end();
    map<VDEdge*, VDEdge*>::iterator               edge_end   = edgeMap.end();

    origin.m_GlobalPartialedgeList.reset4Loop();
    while ( origin.m_GlobalPartialedgeList.setNext4Loop() ) {
        VDPartialEdge* currOriPrEdge  = origin.m_GlobalPartialedgeList.getpEntity();

        VDPartialEdge* currThisPrEdge     = rg_NULL;
        VDLoop*        currLoop           = rg_NULL;
        VDEdge*        currEdge           = rg_NULL;
        VDPartialEdge* nextPrEdgeInRCycle = rg_NULL;
        VDPartialEdge* nextPrEdgeInLoop   = rg_NULL;
        VDPartialEdge* prevPrEdgeInLoop   = rg_NULL;

        it_predge = prEdgeMap.find( currOriPrEdge );
        if ( it_predge != predge_end ) {
            currThisPrEdge = it_predge->second;
        }

        it_loop = loopMap.find( currOriPrEdge->getLoop() );
        if ( it_loop != loop_end ) {
            currLoop = it_loop->second;
        }

        it_edge = edgeMap.find( currOriPrEdge->getOriginalEdge() );
        if ( it_edge != edge_end ) {
            currEdge = it_edge->second;
        }

        it_predge = prEdgeMap.find( currOriPrEdge->getNextPartialEdgeInRadialCycle() );
        if ( it_predge != predge_end ) {
            nextPrEdgeInRCycle = it_predge->second;
        }

        it_predge = prEdgeMap.find( currOriPrEdge->getNextPartialEdgeInLoop() );
        if ( it_predge != predge_end ) {
            nextPrEdgeInLoop = it_predge->second;
        }

        it_predge = prEdgeMap.find( currOriPrEdge->getPreviousPartialEdgeInLoop() );
        if ( it_predge != predge_end ) {
            prevPrEdgeInLoop = it_predge->second;
        }

        currThisPrEdge->setPartialEdge( currLoop, currEdge, 
                                        nextPrEdgeInRCycle, nextPrEdgeInLoop, prevPrEdgeInLoop, 
                                        currOriPrEdge->isRightOrientationInLoop() );
    }
}



void VoronoiDiagram3D::duplicateTopologyOfEdges( const VoronoiDiagram3D& origin, 
                                   map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                   map<VDEdge*, VDEdge*>&               edgeMap,
                                   map<VDVertex*, VDVertex*>&           vertexMap)
{
    map<VDPartialEdge*, VDPartialEdge*>::iterator it_predge;
    map<VDEdge*, VDEdge*>::iterator               it_edge;
    map<VDVertex*, VDVertex*>::iterator           it_vtx;

    map<VDPartialEdge*, VDPartialEdge*>::iterator predge_end = prEdgeMap.end();
    map<VDEdge*, VDEdge*>::iterator               edge_end   = edgeMap.end();
    map<VDVertex*, VDVertex*>::iterator           vtx_end    = vertexMap.end();

    origin.m_GlobalEdgeList.reset4Loop();
    while ( origin.m_GlobalEdgeList.setNext4Loop() ) {
        VDEdge* currOriEdge  = origin.m_GlobalEdgeList.getpEntity();
        VDEdge* currThisEdge = edgeMap.find( currOriEdge )->second;

        if ( currOriEdge->getID() == 4 ) {
            int stop = 1;
        }


        VDVertex* startVertex = rg_NULL;
        it_vtx = vertexMap.find( currOriEdge->getStartVertex() );
        if ( it_vtx != vtx_end ) {
            startVertex = it_vtx->second;
        }

        VDVertex* endVertex   = rg_NULL;
        it_vtx = vertexMap.find( currOriEdge->getEndVertex() );
        if ( it_vtx != vtx_end ) {
            endVertex = it_vtx->second;
        }

        VDPartialEdge* prEdge = rg_NULL;
        it_predge = prEdgeMap.find( currOriEdge->getPartialEdge() );
        if ( it_predge != predge_end ) {
            prEdge = it_predge->second;
        }

        currThisEdge->setEdgeType( currOriEdge->getEdgeType() );
        currThisEdge->setStartVertex( startVertex );
        currThisEdge->setEndVertex( endVertex );
        currThisEdge->setPartialEdge( prEdge );
        if ( currThisEdge->isOnInfinity() ) {
            currThisEdge->setInfinity();
        }

        currThisEdge->connectBetaFace( currOriEdge->getBetaFace() );
    }

}



void VoronoiDiagram3D::duplicateTopologyOfVertices( const VoronoiDiagram3D& origin, 
                                   map<VDEdge*, VDEdge*>&               edgeMap,
                                   map<VDVertex*, VDVertex*>&           vertexMap)
{
    map<VDEdge*, VDEdge*>::iterator     it_edge;
    map<VDVertex*, VDVertex*>::iterator it_vtx;
    map<VDEdge*, VDEdge*>::iterator     edge_end = edgeMap.end();
    map<VDVertex*, VDVertex*>::iterator vtx_end  = vertexMap.end();


    origin.m_GlobalVertexList.reset4Loop();
    while ( origin.m_GlobalVertexList.setNext4Loop() ) {
        VDVertex* currOriVtx  = origin.m_GlobalVertexList.getpEntity();
        VDVertex* currThisVtx = rg_NULL;

        it_vtx = vertexMap.find( currOriVtx );
        if ( it_vtx != vtx_end ) {
            currThisVtx = it_vtx->second;
        }

        for ( rg_INT i=0; i<4; i++ ) {
            VDEdge*   incidentEdge = rg_NULL;
            it_edge = edgeMap.find( currOriVtx->getIncidentEdge(i) );
            if ( it_edge != edge_end ) {
                incidentEdge = it_edge->second;
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


///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
void VoronoiDiagram3D::fileOutGlobalEdgeList( ofstream& fout )
{
    VDEdge* currEdge = rg_NULL;

	m_GlobalEdgeList.reset4Loop();
	while( m_GlobalEdgeList.setNext4Loop() )  
    {
		currEdge = m_GlobalEdgeList.getpEntity();

        currEdge->fileOutVDEdge( fout );
    }
}

void VoronoiDiagram3D::fileOutEdgeLengthsOfDelaunayTetrahedron( ofstream& fout )
{
    VDVertex* currVertex = rg_NULL;


    fout << "//  0 1 2 3 : dist(0,1) dist(0,2) dist(0,3) dist(1,2) dist(1,3) dist(2,3)" << endl << endl;

	m_GlobalVertexList.reset4Loop();
	while( m_GlobalVertexList.setNext4Loop() )  
    {
		currVertex = m_GlobalVertexList.getpEntity();

        if ( currVertex->isOnInfinity() )
            continue;
        
        VDCell* cellsToDefineVertex[NUM_DEFINING_CELLS_OF_ONE_VERTEX];

        currVertex->inquireAllCellsToDefineThisVertex( cellsToDefineVertex );

        qsort( (void *) cellsToDefineVertex, NUM_DEFINING_CELLS_OF_ONE_VERTEX, sizeof(VDCell*), compareVDCellByID);

        rg_Point3D centersOfBallsToDefineVertex[NUM_DEFINING_CELLS_OF_ONE_VERTEX];
        rg_INT i=0;
		for ( i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        {
            centersOfBallsToDefineVertex[i] = ((BallGenerator*) (cellsToDefineVertex[i]->getGenerator()))->getCenter();
        }

        rg_REAL edgeLengths[6];
        // cell[0] : cell[1]
        edgeLengths[0] = centersOfBallsToDefineVertex[0].distance( centersOfBallsToDefineVertex[1] );
        // cell[0] : cell[2]
        edgeLengths[1] = centersOfBallsToDefineVertex[0].distance( centersOfBallsToDefineVertex[2] );
        // cell[0] : cell[3]
        edgeLengths[2] = centersOfBallsToDefineVertex[0].distance( centersOfBallsToDefineVertex[3] );
        // cell[1] : cell[2]
        edgeLengths[3] = centersOfBallsToDefineVertex[1].distance( centersOfBallsToDefineVertex[2] );
        // cell[1] : cell[3]
        edgeLengths[4] = centersOfBallsToDefineVertex[1].distance( centersOfBallsToDefineVertex[3] );
        // cell[2] : cell[3]
        edgeLengths[5] = centersOfBallsToDefineVertex[2].distance( centersOfBallsToDefineVertex[3] );

        for ( i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++ )
            fout << cellsToDefineVertex[i]->getID() << "\t";

        fout << " : " << "\t";
        for ( i=0; i<6; i++)
            fout << edgeLengths[i] << "\t";

        fout << endl;
    }
}

void VoronoiDiagram3D::writeCellList( ofstream& fout )
{
    fout << endl << endl;
    fout << "===================================================================" << endl;
    fout << "=  Voronoi Cells" << endl;
    fout << "=" << endl << endl;

    
    VDCell* currCell = rg_NULL;

	m_GlobalCellList.reset4Loop();
	while( m_GlobalCellList.setNext4Loop() )  
    {
		currCell = m_GlobalCellList.getpEntity();
        
        fout << "C_" << currCell->getID() << "\t";
        if ( currCell->isBounded() )
            fout << "(B)" << "\t";
        else
            fout << "(UB)" << "\t";

        fout << "G_" << currCell->getGenerator()->getID() << "\t";

        VDFace* currFace = rg_NULL;
        rg_dList<VDFace*>* faceList = currCell->getBoungindFaces();

        fout << faceList->getSize() << "\t";
        faceList->reset4Loop();
        while ( faceList->setNext4Loop() )  
        {
		    currFace = faceList->getEntity();

            fout << "F_" << currFace->getID() << "\t";
        }
        fout << endl;
    }
}

void VoronoiDiagram3D::writeFaceList( ofstream& fout )
{
    fout << endl << endl;
    fout << "===================================================================" << endl;
    fout << "=  Voronoi Faces" << endl;
    fout << "=" << endl << endl;

    VDFace* currFace = rg_NULL;

	m_GlobalFaceList.reset4Loop();
	while( m_GlobalFaceList.setNext4Loop() )  
    {
		currFace = m_GlobalFaceList.getpEntity();

        fout << "F_" << currFace->getID() << "\t";
        if ( currFace->isBounded() )
            fout << "(B)" << "\t";
        else
            fout << "(UB)" << "\t";

        fout << "C_" << currFace->getRightCell()->getID() << "\t";
        fout << "C_" << currFace->getLeftCell()->getID() << "\t";


        VDLoop* currLoop = rg_NULL;
        rg_dList<VDLoop*>* loopList = currFace->getLoops();
        fout << loopList->getSize() << "\t";
        loopList->reset4Loop();
        while ( loopList->setNext4Loop() )  
        {
		    currLoop = loopList->getEntity();

            fout << "L_" << currLoop->getID() << "\t";
        }
        fout << endl;
    }
}

void VoronoiDiagram3D::writeLoopList( ofstream& fout )
{
    fout << endl << endl;
    fout << "===================================================================" << endl;
    fout << "=  Loops" << endl;
    fout << "=" << endl << endl;

    VDLoop* currLoop = rg_NULL;

    m_GlobalLoopList.reset4Loop();
    while ( m_GlobalLoopList.setNext4Loop() )  
    {
		currLoop = m_GlobalLoopList.getpEntity();

        fout << "L_" << currLoop->getID() << "\t";
        fout << "F_" << currLoop->getFace()->getID() << "\t";

        VDPartialEdge* startPrEdge = currLoop->getPartialEdge();
        VDPartialEdge* currPrEdge  = startPrEdge;

        do 
        {
            fout << "PE_" << currPrEdge->getID();
            VDVertex* sVet = currPrEdge->getOriginalEdge()->getStartVertex();
            VDVertex* eVet = currPrEdge->getOriginalEdge()->getEndVertex();
            if ( currPrEdge->isRightOrientationInLoop() )

            {
                fout << "(+)\t";
                //fout << "(+){V_" << sVet->getID() << " - V_" << eVet->getID() << "}\t";
            }
            else
            {
                fout << "(-)\t";
                //fout << "(-){V_" << eVet->getID() << " - V_" << sVet->getID() << "}\t";
            }

            currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
        } while ( startPrEdge != currPrEdge );

        fout << endl << "\t\t";
        if ( currLoop->getFace()->isBounded() )
        {
            if ( startPrEdge->isRightOrientationInLoop() )
                fout << "V_" << startPrEdge->getOriginalEdge()->getStartVertex()->getID();
            else
                fout << "V_" << startPrEdge->getOriginalEdge()->getEndVertex()->getID();

            currPrEdge  = startPrEdge;
            do 
            {
                fout << " -> ";
                if ( currPrEdge->isRightOrientationInLoop() )
                    fout << "V_" << currPrEdge->getOriginalEdge()->getEndVertex()->getID();
                else
                    fout << "V_" << currPrEdge->getOriginalEdge()->getStartVertex()->getID();

                currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
            } while ( startPrEdge != currPrEdge );
        }
        else
        {
            while ( startPrEdge->getOriginalEdge()->isBoundedEdge() )
            {
                startPrEdge = startPrEdge->getNextPartialEdgeInLoop();
            }

            if ( !startPrEdge->getNextPartialEdgeInLoop()->getOriginalEdge()->isBoundedEdge() )
                startPrEdge = startPrEdge->getNextPartialEdgeInLoop();

            if ( startPrEdge->isRightOrientationInLoop() )
                fout << "V_" << startPrEdge->getOriginalEdge()->getStartVertex()->getID();
            else
                fout << "V_" << startPrEdge->getOriginalEdge()->getEndVertex()->getID();

            currPrEdge  = startPrEdge;
            do 
            {
                fout << " -> ";
                if ( currPrEdge->isRightOrientationInLoop() )
                    fout << "V_" << currPrEdge->getOriginalEdge()->getEndVertex()->getID();
                else
                    fout << "V_" << currPrEdge->getOriginalEdge()->getStartVertex()->getID();

                currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
            } while ( startPrEdge != currPrEdge );
        }
        fout << endl;
    }
}

void VoronoiDiagram3D::writePrEdgeList( ofstream& fout )
{
    fout << endl << endl;
    fout << "===================================================================" << endl;
    fout << "=  Partial Edges" << endl;
    fout << "=" << endl << endl;

    VDPartialEdge* currPrEdge = rg_NULL;

    m_GlobalPartialedgeList.reset4Loop();
    while ( m_GlobalPartialedgeList.setNext4Loop() )  
    {
		currPrEdge = m_GlobalPartialedgeList.getpEntity();

        fout << "PE_" << currPrEdge->getID() << "\t";
        if ( currPrEdge->isRightOrientationInLoop() )
            fout << "(+)\t";
        else
            fout << "(-)\t";
        fout << "E_" << currPrEdge->getOriginalEdge()->getID() << "\t";
        fout << "L_" << currPrEdge->getLoop()->getID() << "\t";
        fout << endl;
    }
}

void VoronoiDiagram3D::writeEdgeList( ofstream& fout )
{
    fout << endl << endl;
    fout << "===================================================================" << endl;
    fout << "=  Voronoi Edges" << endl;
    fout << "=" << endl << endl;

    VDEdge* currEdge = rg_NULL;

    m_GlobalEdgeList.reset4Loop();
    while ( m_GlobalEdgeList.setNext4Loop() )  
    {
		currEdge = m_GlobalEdgeList.getpEntity();

        fout << "E_" << currEdge->getID() << "\t";
        if ( currEdge->isBoundedEdge() )
            fout << "(B)" << "\t";
        else
            fout << "(UB)" << "\t";

        fout << "V_" << currEdge->getStartVertex()->getID() << "\t"; 
        fout << "V_" << currEdge->getEndVertex()->getID() << "\t"; 

        VDPartialEdge* startPrEdge = currEdge->getPartialEdge();
        VDPartialEdge* currPrEdge  = startPrEdge;
        do
        {
            fout << "PE_" << currPrEdge->getID() << "\t";

            currPrEdge = currPrEdge->getNextPartialEdgeInRadialCycle();
        } while ( startPrEdge != currPrEdge );

        VDCell* cellsInCCWSharingEdge[3];
        currEdge->inquireIntoOnlyEdgeSharingCellsInCCW( cellsInCCWSharingEdge );
        for ( int i=0; i<3; i++ )
        {
            fout << "C_" << cellsInCCWSharingEdge[i]->getID() << "\t";
        }


        fout << endl;
    }
}

void VoronoiDiagram3D::writeVertexList( ofstream& fout )
{
    fout << endl << endl;
    fout << "===================================================================" << endl;
    fout << "=  Voronoi Vertices" << endl;
    fout << "=" << endl << endl;

    VDVertex* currVertex = rg_NULL;

    m_GlobalVertexList.reset4Loop();
    while ( m_GlobalVertexList.setNext4Loop() )  
    {
		currVertex = m_GlobalVertexList.getpEntity();

        fout << "V_" << currVertex->getID() << "\t";
        if ( currVertex->isOnInfinity() )
        {
            fout << "(UB)" << "\t";
        }
        else
        {
            fout << "(B)" << "\t";

            rg_Point3D point( currVertex->getPoint() );
            fout << point.getX() << "\t" << point.getY() << "\t" << point.getZ() << "\t";

            for (int i=0; i<NUM_INCIDENT_EDGES_OF_ONE_VERTEX; i++)
            {
                fout << "E_" << currVertex->getIncidentEdge(i)->getID() << "\t";
            }
        }
        fout << endl;
    }
}

void VoronoiDiagram3D::writeGeneratorList( ofstream& fout )
{
}



///////////////////////////////////////////////////////////////////////////////
//
//  topological operators..
VDVertex* VoronoiDiagram3D::createVertex(const rg_FLAG& onInfinity, const rg_Point3D& point, const rg_REAL& radiusOfTS)
{
    VDVertex* newVtx = m_GlobalVertexList.add( VDVertex(m_GlobalVertexList.getSize(), onInfinity, point, radiusOfTS) );

    return newVtx;
}


    
VDEdge*   VoronoiDiagram3D::createEdge(VDVertex* startVtx, VDVertex* endVtx)
{
    VDEdge* newEdge = m_GlobalEdgeList.add( VDEdge( m_GlobalEdgeList.getSize(), startVtx, endVtx ) );

    VDPartialEdge* newPrEdge1 = m_GlobalPartialedgeList.add( VDPartialEdge( m_GlobalPartialedgeList.getSize(), newEdge ) );
    VDPartialEdge* newPrEdge2 = m_GlobalPartialedgeList.add( VDPartialEdge( m_GlobalPartialedgeList.getSize(), newEdge ) );
    VDPartialEdge* newPrEdge3 = m_GlobalPartialedgeList.add( VDPartialEdge( m_GlobalPartialedgeList.getSize(), newEdge ) );

    newPrEdge1->setNextPartialEdgeInRadialCycle( newPrEdge2 );
    newPrEdge2->setNextPartialEdgeInRadialCycle( newPrEdge3 );
    newPrEdge3->setNextPartialEdgeInRadialCycle( newPrEdge1 );

    newEdge->setPartialEdge( newPrEdge1 );

    return newEdge;
}



VDPartialEdge*  VoronoiDiagram3D::createPrEdge(VDEdge* oriEdge)
{
    VDPartialEdge* newPrEdge1 = m_GlobalPartialedgeList.add( VDPartialEdge( m_GlobalPartialedgeList.getSize(), oriEdge ) );
    VDPartialEdge* newPrEdge2 = m_GlobalPartialedgeList.add( VDPartialEdge( m_GlobalPartialedgeList.getSize(), oriEdge ) );
    VDPartialEdge* newPrEdge3 = m_GlobalPartialedgeList.add( VDPartialEdge( m_GlobalPartialedgeList.getSize(), oriEdge ) );

    newPrEdge1->setNextPartialEdgeInRadialCycle( newPrEdge2 );
    newPrEdge2->setNextPartialEdgeInRadialCycle( newPrEdge3 );
    newPrEdge3->setNextPartialEdgeInRadialCycle( newPrEdge1 );

    oriEdge->setPartialEdge( newPrEdge1 );

    return newPrEdge1;
}




VDCell*        VoronoiDiagram3D::createCell(  const VDCell& cell)
{
    return m_GlobalCellList.add( cell );
}



VDFace*        VoronoiDiagram3D::createFace(  const VDFace& face)
{
    return m_GlobalFaceList.add( face );
}



VDLoop*        VoronoiDiagram3D::createLoop(  const VDLoop& loop)
{
    return m_GlobalLoopList.add( loop );
}



VDPartialEdge* VoronoiDiagram3D::createPrEdge(const VDPartialEdge& prEdge)
{
    return m_GlobalPartialedgeList.add( prEdge );
}



VDEdge*        VoronoiDiagram3D::createEdge(  const VDEdge& edge)
{
    return m_GlobalEdgeList.add( edge );
}



VDVertex*      VoronoiDiagram3D::createVertex(const VDVertex& vertex)
{
    return m_GlobalVertexList.add( vertex );
}





void VoronoiDiagram3D::splitEdge(VDEdge* edge, VDVertex* splitVtx, VDEdge*& newEdge1, VDEdge*& newEdge2)
{
    //  sVtx    edge    eVtx        sVtx  edge=newEdge1  splitVtx  newEdge2  eVtx
    //    *--------------->*    =>    *---------------------->*--------------->*

    VDVertex* startVtx = edge->getStartVertex();
    VDVertex* endVtx   = edge->getEndVertex();
    VDPartialEdge* prEdge[3] = { edge->getPartialEdge(), prEdge[0]->getNextPartialEdgeInRadialCycle(), prEdge[1]->getNextPartialEdgeInRadialCycle()};

    //  create newEdge2
    newEdge2 = createEdge( splitVtx, endVtx );
    VDPartialEdge* newPrEdge[3] = { newEdge2->getPartialEdge(), newPrEdge[0]->getNextPartialEdgeInRadialCycle(), newPrEdge[1]->getNextPartialEdgeInRadialCycle() };


    if ( edge->getID() == 140 || edge->getID() == 23123 || edge->getID() == 23124 || edge->getID() == 22954 || edge->getID() == 289 || edge->getID() == 23155 || edge->getID() == 271 ) {
        int stop = 1;
    }


    //  set topology of prEdges of newEdge2
    rg_INT i=0;
    for ( i=0; i<3; i++ ) {
        newPrEdge[i]->setLoop( prEdge[i]->getLoop() );
        newPrEdge[i]->isRightOrientationInLoop( prEdge[i]->isRightOrientationInLoop() );

        VDPartialEdge* nextPrEdgeInLoop = rg_NULL;
        VDPartialEdge* prevPrEdgeInLoop = rg_NULL;

        if ( prEdge[i]->isRightOrientationInLoop() ) {
            nextPrEdgeInLoop = prEdge[i]->getNextPartialEdgeInLoop();
            prevPrEdgeInLoop = prEdge[i];
        }
        else {
            nextPrEdgeInLoop = prEdge[i];
            prevPrEdgeInLoop = prEdge[i]->getPreviousPartialEdgeInLoop();
        }

        newPrEdge[i]->setPreviousPartialEdgeInLoop( prevPrEdgeInLoop );
        newPrEdge[i]->setNextPartialEdgeInLoop(     nextPrEdgeInLoop );

        prevPrEdgeInLoop->setNextPartialEdgeInLoop(     newPrEdge[i] );
        nextPrEdgeInLoop->setPreviousPartialEdgeInLoop( newPrEdge[i] );
    }

    //  set topology of newEdge1
    newEdge1 = edge;
    newEdge1->setEndVertex( splitVtx );
    splitVtx->setIncidentEdge( 0, newEdge1 );

    //  set topology of newEdge2
    splitVtx->setIncidentEdge( 1, newEdge2 );
    for ( i=0; i<4; i++ ) {
        if ( endVtx->getIncidentEdge(i) == edge ) {
            endVtx->setIncidentEdge( i, newEdge2 );
            break;
        }
    }

    //  set linkage bet. VD and QT
    newEdge2->connectBetaFace( newEdge1->getBetaFace() );

}


VDVertex* VoronoiDiagram3D::findClosestVertexToGivenPoint( const rg_Point3D& givenPoint )
{
    // 시간이 없어 O(n) 방법을 사용함. vertex들을 전부 검색.
    // Youngsong Cho 2004. 06. 04.

    rg_REAL minDistance     = DBL_MAX;
    VDVertex* closestVertex = rg_NULL;

    VDVertex* currVertex = rg_NULL;

	m_GlobalVertexList.reset4Loop();
	while( m_GlobalVertexList.setNext4Loop() )  
    {
		currVertex = m_GlobalVertexList.getpEntity();

        if ( currVertex->isOnInfinity() || currVertex->isTangible() == NON_TANGIBLE )
            continue;

        rg_REAL dist = givenPoint.distance( currVertex->getPoint() );
        if ( rg_LT( dist, minDistance) )
        {
            minDistance   = dist;
            closestVertex = currVertex;
        }
    }

    return closestVertex;
}


///////////////////////////////////////////////////////////////////////////////
//
//  geometric operators..
rg_dList<BallGenerator*>* VoronoiDiagram3D::computeConvexHullOfGeneter()
{
    rg_dList<BallGenerator*>* convexHull = new rg_dList<BallGenerator*>;

    VDEdge* currEdge = rg_NULL;

	m_GlobalEdgeList.reset4Loop();
	while( m_GlobalEdgeList.setNext4Loop() )  
    {
		currEdge = m_GlobalEdgeList.getpEntity();

        if ( currEdge->isOnInfinity() )
            continue;

        if ( !currEdge->isBoundedEdge() )
        {
            rg_dList<VDCell*> onlyEdgeSharingCellsInCCW;
            
            currEdge->inquireIntoOnlyEdgeSharingCellsInCCW( onlyEdgeSharingCellsInCCW );

            onlyEdgeSharingCellsInCCW.reset4Loop();
	        while( onlyEdgeSharingCellsInCCW.setNext4Loop() )  
            {
                convexHull->addTail( onlyEdgeSharingCellsInCCW.getEntity()->getGenerator() );
            }
        }
    }

    return convexHull;
}

rg_dList<VDFace*>* VoronoiDiagram3D::obtainFacesWithIntersectingGenerators()
{
    rg_dList<VDFace*>* facesWithIntersectingGenerators = new rg_dList<VDFace*>;

    VDFace* currFace = rg_NULL;

	m_GlobalFaceList.reset4Loop();
	while( m_GlobalFaceList.setNext4Loop() )  
    {
		currFace = m_GlobalFaceList.getpEntity();
        
        if ( currFace->getLeftCell()->getID() == -1 || currFace->getRightCell()->getID() == -1 )
            continue;

        BallGenerator* generatorOfRightCell = currFace->getRightCell()->getGenerator();
        BallGenerator* generatorOfLeftCell  = currFace->getLeftCell()->getGenerator();

        if ( generatorOfRightCell->isThisGeneratorIntersectedWith(generatorOfLeftCell) )
        {
            facesWithIntersectingGenerators->addTail( currFace );
        }

    }

    return facesWithIntersectingGenerators;
}

rg_dList<VDEdge*>* VoronoiDiagram3D::obtainEdgesWithoutIntersectingGenerators()
{
    rg_dList<VDEdge*>* edgesWithoutIntersectingGenerators = new rg_dList<VDEdge*>;

    VDEdge* currEdge = rg_NULL;
    VDCell*    cellsSharedEdge[3]      = {rg_NULL, rg_NULL, rg_NULL};
    BallGenerator* generatorsSharedEdge[3] = {rg_NULL, rg_NULL, rg_NULL};

	m_GlobalEdgeList.reset4Loop();
	while( m_GlobalEdgeList.setNext4Loop() )  
    {
		currEdge = m_GlobalEdgeList.getpEntity();
        
        currEdge->inquireIntoOnlyEdgeSharingCellsInCCW( cellsSharedEdge );

        for ( rg_INT i=0; i<3; i++)
            generatorsSharedEdge[i] = cellsSharedEdge[i]->getGenerator();

        if (    !generatorsSharedEdge[0]->isThisGeneratorIntersectedWith( generatorsSharedEdge[1] ) 
             && !generatorsSharedEdge[1]->isThisGeneratorIntersectedWith( generatorsSharedEdge[2] ) 
             && !generatorsSharedEdge[2]->isThisGeneratorIntersectedWith( generatorsSharedEdge[0] ) )
        {
            edgesWithoutIntersectingGenerators->addTail( currEdge );
        }
    }

    return edgesWithoutIntersectingGenerators;
}


/*
rg_dList<VDCell*>* VoronoiDiagram3D::obtainCellsToIntersectWithSectionalYZPlane(const rg_REAL& xOfPlane)
{
    rg_dList<VDCell*>* cellsToIntersectWithSectionalYZPlane = new rg_dList<VDCell*>;

    VDCell* currCell = rg_NULL;

	m_GlobalCellList.reset4Loop();
	while( m_GlobalCellList.setNext4Loop() )  
    {
		currCell = m_GlobalCellList.getpEntity();

        rg_dList<VDFace*>* boundingFaces = currCell->getBoungindFaces();

        VDFace* currFace = rg_NULL;

	    boundingFaces->reset4Loop();
	    while( boundingFaces->setNext4Loop() )  
        {
		    currFace = boundingFaces->getEntity();

            if ( currFace->isThereIntersectionWithSectionalYZPlane( xOfPlane ) )
            {
                cellsToIntersectWithSectionalYZPlane->addTail( currCell );
                break;
            }
        }
    }

    return cellsToIntersectWithSectionalYZPlane;
}
*/



//void VoronoiDiagram3D::computeVoronoiFaces()
//{
//    VDFace* currFace = rg_NULL;
//
//	m_GlobalFaceList.reset4Loop();
//	while( m_GlobalFaceList.setNext4Loop() )  
//    {
//		currFace = m_GlobalFaceList.getpEntity();
//        
//        if ( currFace->getLeftCell()->getID() < 0 || currFace->getRightCell()->getID() < 0 )
//            continue;
//        
//        if ( currFace->isOnInfinity() )
//            continue;
//
//        if ( currFace->getFaceMesh() == rg_NULL )
//            currFace->makeMeshForVoronoiFace();
//    }
//    m_makeVoronoiFace = rg_TRUE;
//}
//
//
//
//void VoronoiDiagram3D::computeVoronoiFaces(rg_INT& numTriangles)
//{
//    numTriangles = 0;
//
//    VDFace* currFace = rg_NULL;
//	m_GlobalFaceList.reset4Loop();
//	while( m_GlobalFaceList.setNext4Loop() )  
//    {
//		currFace = m_GlobalFaceList.getpEntity();
//        
//        if ( currFace->getLeftCell()->getID() < 0 || currFace->getRightCell()->getID() < 0 )
//            continue;
//        
//        if ( currFace->isOnInfinity() || !currFace->isBounded())
//            continue;
//
//        if ( currFace->getFaceMesh() == rg_NULL )
//        {
//            currFace->makeMeshForVoronoiFace();
//            
//            rg_INT numTrianglesOnOneMesh = currFace->getFaceMesh()->getSize();
//
//            numTriangles += (numTrianglesOnOneMesh/6);
//        }
//    }
//    m_makeVoronoiFace = rg_TRUE;
//}

rg_INT VoronoiDiagram3D::approximateStorageCost()
{
    rg_INT sumOfBoundingFacesInCell = 0;
    rg_INT sumOfLoopsInFace         = 0;

    VDCell* currCell = rg_NULL;
	m_GlobalCellList.reset4Loop();
	while( m_GlobalCellList.setNext4Loop() )  
    {
		currCell = m_GlobalCellList.getpEntity();       
        sumOfBoundingFacesInCell += currCell->getNumOfBoundingFaces();
    }

    VDFace* currFace = rg_NULL;
	m_GlobalFaceList.reset4Loop();
	while( m_GlobalFaceList.setNext4Loop() )  
    {
		currFace = m_GlobalFaceList.getpEntity();
        sumOfLoopsInFace += currFace->getNumOfLoops();
    }

    rg_INT storageForCell   = (4 + 12 + (8*sumOfBoundingFacesInCell) + 8)*m_GlobalCellList.getSize() + 12;
    rg_INT storageForFace   = (4 + 8 + 12 + (8*sumOfLoopsInFace) + 8)*m_GlobalFaceList.getSize() + 12;
    rg_INT storageForLoop   = (10 + 8)*m_GlobalLoopList.getSize() + 12;
    rg_INT storageForPrEdge = (22+8)*m_GlobalPartialedgeList.getSize() + 12;
    rg_INT storageForEdge   = (4 + 12 + 8)*m_GlobalEdgeList.getSize() + 12;
    rg_INT storageForVertex = (24 + 16 + 8)*m_GlobalVertexList.getSize() + 12;


    long storageCostForEVDS = storageForCell + storageForFace + storageForLoop
                                + storageForPrEdge + storageForEdge + storageForVertex;
    
    return storageCostForEVDS;
}

