#include "BetaUniverse.h"

#include "VDCell.h"
#include "VDFace.h"
#include "VDEdge.h"
#include "VDVertex.h"
#include "rg_BallGenerator.h"
using namespace V::GeometryTier;

#include <time.h>
#include <stdlib.h> 

#include <list>
#include <set>
using namespace std;



BetaUniverse::BetaUniverse()
{
    m_numFiniteCells    = 0;
    m_numFiniteFaces    = 0;
    m_numFiniteEdges    = 0;
    m_numFiniteVertices = 0;

    m_sortedFiniteCells     = rg_NULL;
    m_sortedFiniteFaces     = rg_NULL;
    m_sortedFiniteEdges     = rg_NULL;
    m_sortedFiniteVertices  = rg_NULL;
	//m_simplexQueryingMachine = rg_NULL;
}



BetaUniverse::~BetaUniverse()
{
    clean();
}



// 
// rg_dList<Ball>* BetaUniverse::getBallList() 
// {
//     return &m_balls;
// }
// 
// 
// 
// rg_dList<BetaCell>* BetaUniverse::getCellList()
// {
//     return &m_cells;
// }
// 
// 
// 
// rg_dList<BetaFace>* BetaUniverse::getFaceList()
// {
//     return &m_faces;
// }
// 
// 
// 
// rg_dList<BetaEdge>* BetaUniverse::getEdgeList()
// {
//     return &m_edges;
// }
// 
// 
// 
// rg_dList<BetaVertex>* BetaUniverse::getVertexList()
// {
//     return &m_vertices;
// }
// 
// 
// 
// rg_INT          BetaUniverse::getNumBalls() const
// {
//     return m_balls.getSize();
// }
// 
// 
// 
// rg_INT BetaUniverse::getNumCells() const
// {
//     return m_numFiniteCells;
// }
// 
// 
// 
// rg_INT BetaUniverse::getNumFaces() const
// {
//     return m_numFiniteFaces;
// }
// 
// 
// 
// rg_INT BetaUniverse::getNumEdges() const
// {
//     return m_numFiniteEdges;
// }
// 
// 
// 
// rg_INT BetaUniverse::getNumVertices() const
// {
//     return m_numFiniteVertices;
// }
// 
// 
// 
// BetaCell**   BetaUniverse::getSortedCells() 
// {
//     return m_sortedFiniteCells;
// }
// 
// 
// 
// BetaFace**   BetaUniverse::getSortedFaces()
// {
//     return m_sortedFiniteFaces;
// }
// 
// 
// 
// BetaEdge**   BetaUniverse::getSortedEdges()
// {
//     return m_sortedFiniteEdges;
// }
// 
// 
// 
// BetaVertex** BetaUniverse::getSortedVertices()
// {
//     return m_sortedFiniteVertices;
// }
// 


void BetaUniverse::enlargeRadiiOfBallsNUpdateBetaSpan(const rg_REAL& deltaOfBallRadii)
{
	// enlarge the radii of each ball
	Ball* currBall = rg_NULL;
	m_balls.reset4Loop();
	while (m_balls.setNext4Loop())
	{
		currBall = m_balls.getpEntity();
		rg_REAL originalRadius = currBall->getGeometry().getRadius();
		Sphere  sphere = currBall->getGeometry();
		sphere.setRadius(originalRadius + deltaOfBallRadii);
		currBall->setGeometry(sphere);
	}

	// recompute beta span
	//updateBetaSpan(-deltaOfBallRadii);
	updateBetaSpan2(-deltaOfBallRadii);
}

void BetaUniverse::shrinkRadiiOfBallsNUpdateBetaSpan(const rg_REAL& deltaOfBallRadii)
{
	// shrink the radii of each ball
	Ball* currBall = rg_NULL;
	m_balls.reset4Loop();
	while (m_balls.setNext4Loop())
	{
		currBall = m_balls.getpEntity();
		rg_REAL originalRadius = currBall->getGeometry().getRadius();
		Sphere  sphere = currBall->getGeometry();
		sphere.setRadius(originalRadius - deltaOfBallRadii);
		currBall->setGeometry(sphere);
	}

	// recompute beta span
	//updateBetaSpan(deltaOfBallRadii);
	updateBetaSpan2(deltaOfBallRadii);
}

void BetaUniverse::resizeRadiiOfBallsNUpdateBetaSpan(const rg_REAL& deltaOfBallRadii)
{
	// resize the radii of each ball
	Ball* currBall = rg_NULL;
	m_balls.reset4Loop();
	while (m_balls.setNext4Loop())
	{
		currBall = m_balls.getpEntity();
		rg_REAL originalRadius = currBall->getGeometry().getRadius();
		Sphere  sphere = currBall->getGeometry();
		sphere.setRadius(originalRadius + deltaOfBallRadii);
		currBall->setGeometry(sphere);
	}

	// recompute beta span
	//updateBetaSpan(-deltaOfBallRadii);
	updateBetaSpan2(-deltaOfBallRadii);
}

void BetaUniverse::updateBetaSpan(const rg_REAL& deltaOfBallRadii)
{
	updateRadiiOfMinimumTangentSpheresForBetaCells(deltaOfBallRadii);
	computeBetaSpan();
	generateSortedFiniteSimplexes();
	rearrangeSimplexIDs();	
}

void BetaUniverse::updateRadiiOfMinimumTangentSpheresForBetaCells(const rg_REAL& deltaOfBallRadii)
{
    BetaCell* currCell = rg_NULL;
	m_cells.reset4Loop();
	while( m_cells.setNext4Loop() )  
	{
		currCell = m_cells.getpEntity();

        if(currCell->isVirtual())
		{
            continue;
        }
        currCell->updateRadiusOfMinimumTangentSphere(deltaOfBallRadii);
	}
}

void BetaUniverse::updateBetaSpan2(const rg_REAL& deltaOfBallRadii)
{
	shiftBetaSpan(deltaOfBallRadii);
// 	generateSortedFiniteSimplexes();
// 	rearrangeSimplexIDs();
}

void BetaUniverse::shiftBetaSpan(const rg_REAL& deltaOfBallRadii)
{
	// For beta cells
    BetaCell* currCell = rg_NULL;
	m_cells.reset4Loop();
	while( m_cells.setNext4Loop() )  {
		currCell = m_cells.getpEntity();

        if ( currCell->isVirtual() ) {
            continue;
        }

        currCell->shiftBetaSpan(deltaOfBallRadii);
	}	

	// For beta faces
    BetaFace* currFace = rg_NULL;
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		currFace = m_faces.getpEntity();

        if ( currFace->isVirtual() ) {
            continue;
        }

        currFace->shiftBetaSpan(deltaOfBallRadii);
	}

	// For beta edges
    BetaEdge* currEdge = rg_NULL;
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		currEdge = m_edges.getpEntity();

        if ( currEdge->isVirtual() ) {
            continue;
        }

        currEdge->shiftBetaSpan(deltaOfBallRadii);
	}

	// For beta vertices
    BetaVertex* currVertex = rg_NULL;
    m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		currVertex = m_vertices.getpEntity();

        if ( currVertex->isVirtual() ) {
            continue;
        }

        currVertex->shiftBetaSpan(deltaOfBallRadii);
	}
}

void BetaUniverse::enlargeRadiiOfBalls(const rg_REAL& deltaOfBallRadii)
{
	// enlarge the radii of each ball
	Ball* currBall = rg_NULL;
	m_balls.reset4Loop();
	while (m_balls.setNext4Loop())
	{
		currBall = m_balls.getpEntity();
		rg_REAL originalRadius = currBall->getGeometry().getRadius();
		Sphere  sphere = currBall->getGeometry();
		sphere.setRadius(originalRadius + deltaOfBallRadii);
		currBall->setGeometry(sphere);
	}
}

void BetaUniverse::shrinkRadiiOfBalls(const rg_REAL& deltaOfBallRadii)
{
	// shrink the radii of each ball
	Ball* currBall = rg_NULL;
	m_balls.reset4Loop();
	while (m_balls.setNext4Loop())
	{
		currBall = m_balls.getpEntity();
		rg_REAL originalRadius = currBall->getGeometry().getRadius();
		Sphere  sphere = currBall->getGeometry();
		sphere.setRadius(originalRadius - deltaOfBallRadii);
		currBall->setGeometry(sphere);
	}	
}

void BetaUniverse::construct(QuasiTriangulation& qtInIWDS) 
{    
    convertFromIWDSToEIWDS(qtInIWDS);
    
    computeBetaSpan();    
    
    generateSortedFiniteSimplexes();


    //  The ID of each real simplex is arranged from 0 by the order in the sorted array.
    //  for example,
    //  m_sortedFiniteCells[0]->getID() = 0;
    //  ...
    //  m_sortedFiniteCells[m_numFiniteCells-1]->getID() = m_numFiniteCells-1;
    //  The ID of each virtual simplex is arranged from m_numFiniteCells by the order in the list.

    //rearrangeSimplexIDs();
    //setDepthOfWorlds();

	// JKKIM ///////////////////////////////////////////
	//m_simplexQueryingMachine = new SimplexQueryingMachineInQTByKDTree(this);
	////////////////////////////////////////////////////
    


    /*
    //  for time analysis
    m_time4BC[0] = 0.0;
    m_time4BC[1] = 0.0;
    m_time4BC[2] = 0.0;
    m_time4BC[3] = 0.0;

    clock_t start, finish;

    start = clock();
        convertFromIWDSToEIWDS(qtInIWDS);
    finish = clock();
    m_time4BC[0] = (double)(finish-start)/CLOCKS_PER_SEC;

    
    start = clock();
        computeBetaSpan();    
    finish = clock();
    m_time4BC[1] = (double)(finish-start)/CLOCKS_PER_SEC;

    
    start = clock();
        generateSortedFiniteSimplexes();
    finish = clock();
    m_time4BC[2] = (double)(finish-start)/CLOCKS_PER_SEC;


    start = clock();
        rearrangeSimplexIDs();
    finish = clock();
    m_time4BC[3] = (double)(finish-start)/CLOCKS_PER_SEC;
    
    */
    
    
}



void BetaUniverse::clean()
{

    if ( m_sortedFiniteCells != rg_NULL ) {
        delete [] m_sortedFiniteCells;
    }
    if ( m_sortedFiniteFaces != rg_NULL ) {
        delete [] m_sortedFiniteFaces;
    }
    if ( m_sortedFiniteEdges != rg_NULL ) {
        delete [] m_sortedFiniteEdges;
    }
    if ( m_sortedFiniteVertices != rg_NULL ) {
        delete [] m_sortedFiniteVertices;
    }
//	if(m_simplexQueryingMachine != rg_NULL) {
//		delete m_simplexQueryingMachine;
//	}

    m_numFiniteCells    = 0;
    m_numFiniteFaces    = 0;
    m_numFiniteEdges    = 0;
    m_numFiniteVertices = 0;

    m_cells.removeAll();
    m_faces.removeAll();
    m_edges.removeAll();
    m_vertices.removeAll();
}



void  BetaUniverse::huntCellsInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const
{   
    startIndex = 0;
    endIndex   = 0;

    if ( beta < m_sortedFiniteCells[0]->getMinValueForValidSimplex() )  {
        startIndex = 0;
        endIndex   = -1;
    }
    else {
        for ( rg_INT i=1; i<m_numFiniteCells; i++)  {
            if ( beta < m_sortedFiniteCells[i]->getMinValueForValidSimplex() ) {
                break;
            }
            else {
                endIndex = i;
            }
        }
    }
}



void  BetaUniverse::huntFacesInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const
{
    startIndex = 0;
    endIndex   = 0;

    if ( beta < m_sortedFiniteFaces[0]->getMinValueForValidSimplex() )  {
        startIndex = 0;
        endIndex   = -1;
    }
    else {
        for ( rg_INT i=1; i<m_numFiniteFaces; i++)  {
            if ( beta < m_sortedFiniteFaces[i]->getMinValueForValidSimplex() ) {
                break;
            }
            else {
                endIndex = i;
            }
        }
    }
}



void  BetaUniverse::huntEdgesInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const
{
    startIndex = 0;
    endIndex   = 0;

    if ( beta < m_sortedFiniteEdges[0]->getMinValueForValidSimplex() )  {
        startIndex = 0;
        endIndex   = -1;
    }
    else {
        for ( rg_INT i=1; i<m_numFiniteEdges; i++)  {
            if ( beta < m_sortedFiniteEdges[i]->getMinValueForValidSimplex() ) {
                break;
            }
            else {
                endIndex = i;
            }
        }
    }
}



void  BetaUniverse::huntVerticesInBetaComplex( const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const
{
    startIndex = 0;
    endIndex   = m_numFiniteVertices-1;
}



void BetaUniverse::huntSignificantCellsInBetaComplex(const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList, const rg_REAL& tolerance /* = rg_MATH_RES  */)  const
{
    for ( rg_INT i=0; i<m_numFiniteCells; i++)  {
        if ( beta >= m_sortedFiniteCells[i]->getMinValueForValidSimplex() ) 
		{
			rg_REAL signedVolume 
				= m_sortedFiniteCells[i]->computeSignedVolume();

			// Only if the singed volume of tetrahedron is positive
			// add up the list
			if(signedVolume > tolerance)
				betaCellList.add( m_sortedFiniteCells[i] );
        }
        else {
            break;
        }
    }
}



void  BetaUniverse::huntCellsInBetaComplex(    const rg_REAL& beta, rg_dList<BetaCell*>&   betaCellList ) const
{
    for ( rg_INT i=0; i<m_numFiniteCells; i++)  {
        if ( beta >= m_sortedFiniteCells[i]->getMinValueForValidSimplex() ) {
            betaCellList.add( m_sortedFiniteCells[i] );
        }
        else {
            break;
        }
    }
}



void  BetaUniverse::huntFacesInBetaComplex(    const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList ) const
{
    for ( rg_INT i=0; i<m_numFiniteFaces; i++)  {
        if ( beta >= m_sortedFiniteFaces[i]->getMinValueForValidSimplex() ) {
            betaFaceList.add( m_sortedFiniteFaces[i] );
        }
        else {
            break;
        }
    }
}



void  BetaUniverse::huntEdgesInBetaComplex(    const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList ) const
{
    for ( rg_INT i=0; i<m_numFiniteEdges; i++)  {
        if ( beta >= m_sortedFiniteEdges[i]->getMinValueForValidSimplex() ) {
            betaEdgeList.add( m_sortedFiniteEdges[i] );
        }
        else {
            break;
        }
    }
}



void  BetaUniverse::huntVerticesInBetaComplex( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const
{
    for ( rg_INT i=0; i<m_numFiniteVertices; i++)  {
        if ( beta >= m_sortedFiniteVertices[i]->getMinValueForValidSimplex() ) {
            betaVertexList.add( m_sortedFiniteVertices[i] );
        }
        else {
            break;
        }
    }
}




void  BetaUniverse::huntFacesInBetaShape(    const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList ) const
{
    for ( rg_INT i=0; i<m_numFiniteFaces; i++)  {
        rg_INT state = m_sortedFiniteFaces[i]->getBoundingState(beta);
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaFaceList.add( m_sortedFiniteFaces[i] );
        }
    }
}



void  BetaUniverse::huntEdgesInBetaShape(    const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList ) const
{
    for ( rg_INT i=0; i<m_numFiniteEdges; i++)  {
        rg_INT state = m_sortedFiniteEdges[i]->getBoundingState(beta);
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaEdgeList.add( m_sortedFiniteEdges[i] );
        }
    }
}



void  BetaUniverse::huntVerticesInBetaShape( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const
{
    for ( rg_INT i=0; i<m_numFiniteVertices; i++)  {
        rg_INT state = m_sortedFiniteVertices[i]->getBoundingState(beta);
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaVertexList.add( m_sortedFiniteVertices[i] );
        }
    }
}



rg_INT BetaUniverse::countNumGroupsOfFaceConnectedExtraneousCellsInBetaComplex(  const rg_REAL& beta ) 
{
    rg_INT numGroups = 0;

    BetaCell* currCell = rg_NULL;    
    m_cells.reset4Loop();
    while ( m_cells.setNext4Loop() ) {
        currCell = m_cells.getpEntity();

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

                BetaCell* neighbor = searchingCell->getNeighborCell(i);
                searchStack.add( neighbor );
            }
        }
    }


    m_cells.reset4Loop();
    while ( m_cells.setNext4Loop() ) {
        m_cells.getpEntity()->isVisited(rg_FALSE);
    }


    return numGroups;
}



rg_INT BetaUniverse::countNumGroupsOfFaceConnectedInteriorCellsInBetaComplex(  const rg_REAL& beta ) const
{
    rg_INT numGroups = 0;

    
    rg_INT i_cell = 0;
    for ( i_cell = 0; i_cell<m_numFiniteCells; i_cell++) {
        BetaCell* currCell = m_sortedFiniteCells[i_cell];

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

                BetaCell* neighbor = searchingCell->getNeighborCell(i);
                searchStack.add( neighbor );
            }
        }
    }


    for ( i_cell = 0; i_cell<m_numFiniteCells; i_cell++) {
        m_sortedFiniteCells[i_cell]->isVisited(rg_FALSE);
    }


    return numGroups;
}



void BetaUniverse::searchGroupOfOutmostExtraneousCellsInBetaComplex(   const rg_REAL& beta, rg_dList<BetaCell*>& outmostExtraneousCells )
{
    BetaVertex* virtualVertex = findVirtualVertex();

    rg_dList<BetaCell*> stackToSearchExtraneousCells;
    virtualVertex->searchCellsInWholeWorld( stackToSearchExtraneousCells );

    while ( stackToSearchExtraneousCells.getSize() > 0 ) {
        BetaCell* searchingCell = stackToSearchExtraneousCells.popFront();

        if ( searchingCell->isVisited() ) {
            continue;
        }

        outmostExtraneousCells.add( searchingCell );
        searchingCell->isVisited(rg_TRUE);

        BetaFace** boundingFace = searchingCell->getFaces();
        for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
            if ( boundingFace[i]->getBoundingState(beta) != EXTRANEOUS_SIMPLEX ) {
                continue;
            }

            BetaCell* neighbor = searchingCell->getNeighborCell(i);
            stackToSearchExtraneousCells.add( neighbor );
        }
    }

    outmostExtraneousCells.reset4Loop();
    while ( outmostExtraneousCells.setNext4Loop() ) {
        outmostExtraneousCells.getEntity()->isVisited(rg_FALSE);
    }
}



void BetaUniverse::searchGroupsOfFaceConnectedExtraneousCellsInBetaComplex( 
                                                           const rg_REAL&                   beta, 
                                                           rg_dList< rg_dList<BetaCell*> >& groupsOfExtraneousCell ) 
{
    rg_dList<BetaCell*> outmostExtraneousCells;
    searchGroupOfOutmostExtraneousCellsInBetaComplex( beta, outmostExtraneousCells );
    outmostExtraneousCells.reset4Loop();
    while ( outmostExtraneousCells.setNext4Loop() ) {
        outmostExtraneousCells.getEntity()->isVisited(rg_TRUE);
    }

    rg_dList<BetaCell*>* groupOfOutmostExtraneousCells = groupsOfExtraneousCell.add( rg_dList<BetaCell*>() );
    groupOfOutmostExtraneousCells->mergeTail( outmostExtraneousCells );

    BetaCell* currCell = rg_NULL;    
    m_cells.reset4Loop();
    while ( m_cells.setNext4Loop() ) {
        currCell = m_cells.getpEntity();

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

                BetaCell* neighbor = searchingCell->getNeighborCell(i);
                searchStack.add( neighbor );
            }
        }
    }


    m_cells.reset4Loop();
    while ( m_cells.setNext4Loop() ) {
        m_cells.getpEntity()->isVisited(rg_FALSE);
    }
}



void BetaUniverse::searchGroupsOfFaceConnectedInteriorCellsInBetaComplex(   
                                                           const rg_REAL&                   beta, 
                                                           rg_dList< rg_dList<BetaCell*> >& groupsOfInteriorCell ) const
{
    rg_INT i_cell = 0;
    for ( i_cell = 0; i_cell<m_numFiniteCells; i_cell++) {
        BetaCell* currCell = m_sortedFiniteCells[i_cell];

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

                BetaCell* neighbor = searchingCell->getNeighborCell(i);
                searchStack.add( neighbor );
            }
        }
    }


    for ( i_cell = 0; i_cell<m_numFiniteCells; i_cell++) {
        m_sortedFiniteCells[i_cell]->isVisited(rg_FALSE);
    }
}



//////////////////////////////////////////////////////////////////////////////
//  by Y. Cho 20130319
rg_INT  BetaUniverse::findAll2AdjacencyAnomaly( list< pair<BetaCell*, BetaCell*> >& anomaliesWith2Adjacency ) const	
{
    set<BetaCell*> cellsInAnomaly;

    for ( rg_INT i=0; i<m_numFiniteCells; i++ ) {
        BetaCell* currCell = m_sortedFiniteCells[i];

        if ( cellsInAnomaly.find( currCell ) != cellsInAnomaly.end() ) {
            continue;
        }
        if ( currCell->isMultiplicity() != 2 ) {
            continue;
        }

        BetaCell* neighbor = currCell->getNeighborCellWithMultiplicity();

        cellsInAnomaly.insert( currCell );
        cellsInAnomaly.insert( neighbor );

        anomaliesWith2Adjacency.push_back( make_pair( currCell, neighbor ) );
    }

    return anomaliesWith2Adjacency.size();
}



rg_INT  BetaUniverse::findAll3AdjacencyAnomaly( list< pair<BetaCell*, BetaCell*> >& anomaliesWith3Adjacency ) const
{
    set<BetaCell*> cellsInAnomaly;

    for ( rg_INT i=0; i<m_numFiniteCells; i++ ) {
        BetaCell* currCell = m_sortedFiniteCells[i];

        if ( cellsInAnomaly.find( currCell ) != cellsInAnomaly.end() ) {
            continue;
        }
        if ( currCell->isMultiplicity() != 3 ) {
            continue;
        }

        BetaCell* neighbor = currCell->getNeighborCellWithMultiplicity();

        cellsInAnomaly.insert( currCell );
        cellsInAnomaly.insert( neighbor );

        anomaliesWith3Adjacency.push_back( make_pair( currCell, neighbor ) );
    }

    return anomaliesWith3Adjacency.size();
}



rg_INT  BetaUniverse::findAll4AdjacencyAnomaly( list< pair<BetaCell*, BetaCell*> >& anomaliesWith4Adjacency ) const
{
    set<BetaCell*> cellsInAnomaly;

    for ( rg_INT i=0; i<m_numFiniteCells; i++ ) {
        BetaCell* currCell = m_sortedFiniteCells[i];

        if ( cellsInAnomaly.find( currCell ) != cellsInAnomaly.end() ) {
            continue;
        }
        if ( currCell->isMultiplicity() != 4 ) {
            continue;
        }

        BetaCell* neighbor = currCell->getNeighborCellWithMultiplicity();

        cellsInAnomaly.insert( currCell );
        cellsInAnomaly.insert( neighbor );

        anomaliesWith4Adjacency.push_back( make_pair( currCell, neighbor ) );
    }

    return anomaliesWith4Adjacency.size();
}



rg_INT  BetaUniverse::findAllSmallWorlds(       list< BetaFace* >& smallWorlds ) const
{
    for ( rg_INT i=0; i<m_numFiniteEdges; i++ ) {
        BetaEdge* currEdge = m_sortedFiniteEdges[i];

        if ( currEdge->isGateToSmallWorlds() ) {
            rg_dList<BetaFace*>* smallWorldList = currEdge->getSmallWorlds();

            BetaFace* currSmallWorld = rg_NULL;
            smallWorldList->reset4Loop();
            while ( smallWorldList->setNext4Loop() ) {
                currSmallWorld = smallWorldList->getEntity();

                smallWorlds.push_back( currSmallWorld );
            }
        }
    }

    return smallWorlds.size();
}



rg_INT  BetaUniverse::findAllIsolatedTriangles( list< BetaFace* >& isolatedTriangles ) const
{
    for ( rg_INT i=0; i<m_numFiniteEdges; i++ ) {
        BetaEdge* currEdge = m_sortedFiniteEdges[i];

        if ( currEdge->isGateToSmallWorlds() ) {
            rg_dList<BetaFace*>* smallWorldList = currEdge->getSmallWorlds();

            BetaFace* currSmallWorld = rg_NULL;
            smallWorldList->reset4Loop();
            while ( smallWorldList->setNext4Loop() ) {
                currSmallWorld = smallWorldList->getEntity();

                if ( currSmallWorld->isIsolatedFace() ) {
                    isolatedTriangles.push_back( currSmallWorld );
                }
            }
        }

    }
        
    return isolatedTriangles.size();
}

    


rg_INT  BetaUniverse::findFaceConnectedCellsInWorld( const BetaFace* const world, list< BetaCell* >& cellsInWorld ) const
{
    if ( !world->isIsolatedFace() ) {
        set<BetaCell*> cellsInThisWorld;

        rg_dList<BetaCell*> cellStack;
        cellStack.pushBack( world->getLeftCell() );

	    while ( cellStack.getSize() > 0 ) {
		    BetaCell* popedCell = cellStack.popFront();
            cellsInThisWorld.insert( popedCell );

            BetaCell* neighbor[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
	        popedCell->getNeighborCells( neighbor );  
		    for ( rg_INT i=0; i<4; i++) {

                if ( cellsInThisWorld.find( neighbor[i] ) == cellsInThisWorld.end() ) {
				    cellStack.pushBack( neighbor[i] );
                }
		    }
        }


        rg_INT numCells = 0;
        set<BetaCell*>::iterator i_cell;
        for ( i_cell=cellsInThisWorld.begin(); i_cell!=cellsInThisWorld.end(); i_cell++ ) {
            BetaCell* currCell = *i_cell;
            if ( !currCell->isVirtual() ) {
                cellsInWorld.push_back( currCell );
            }
        }
    }

    return cellsInWorld.size();
}



BetaVertex* BetaUniverse::findVirtualVertex() const
{
    BetaVertex* virtualVertex = rg_NULL;

    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() ) {
        BetaVertex* currVtx = m_vertices.getpEntity();

        if ( currVtx->isVirtual() ) {
            virtualVertex = currVtx;
            break;
        }
    }

    return virtualVertex;
}



void BetaUniverse::computeBetaSpan()
{
    //computerBetaSpanOfCells();
    //computerBetaSpanOfFaces();
    //computerBetaSpanOfEdges();
    //computerBetaSpanOfVertices();


    m_time4BetaSpan[0] = 0.0;
    m_time4BetaSpan[1] = 0.0;
    m_time4BetaSpan[2] = 0.0;
    m_time4BetaSpan[3] = 0.0;

    clock_t start, finish;

    start = clock();
        computerBetaSpanOfCells();
    finish = clock();
    m_time4BetaSpan[0] = (double)(finish-start)/CLOCKS_PER_SEC;


    start = clock();
        computerBetaSpanOfFaces();
    finish = clock();
    m_time4BetaSpan[1] = (double)(finish-start)/CLOCKS_PER_SEC;


    start = clock();
        computerBetaSpanOfEdges();
    finish = clock();
    m_time4BetaSpan[2] = (double)(finish-start)/CLOCKS_PER_SEC;


    start = clock();
        computerBetaSpanOfVertices();
    finish = clock();
    m_time4BetaSpan[3] = (double)(finish-start)/CLOCKS_PER_SEC;

}



void BetaUniverse::computerBetaSpanOfCells()
{
    BetaCell* currCell = rg_NULL;
	m_cells.reset4Loop();
	while( m_cells.setNext4Loop() )  {
		currCell = m_cells.getpEntity();

        if ( currCell->isVirtual() ) {
            continue;
        }

        currCell->computeBetaSpan();
	}
}



void BetaUniverse::computerBetaSpanOfFaces()
{
    BetaFace* currFace = rg_NULL;
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		currFace = m_faces.getpEntity();

        if ( currFace->isVirtual() ) {
            continue;
        }

        currFace->computeBetaSpan();
	}
}



void BetaUniverse::computerBetaSpanOfEdges()
{
    BetaEdge* currEdge = rg_NULL;
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		currEdge = m_edges.getpEntity();

        if ( currEdge->isVirtual() ) {
            continue;
        }

        currEdge->computeBetaSpan();
	}
}



void BetaUniverse::computerBetaSpanOfVertices()
{
    /*
    rg_REAL     rMax  = -DBL_MAX;

	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		currVertex = m_vertices.getpEntity();

        rg_REAL radius = currVertex->getBall().getRadius();

        if ( radius > rMax ) {
            rMax = radius;
        }
    }
    */

    
    BetaVertex* currVertex = rg_NULL;
    m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		currVertex = m_vertices.getpEntity();

        if ( currVertex->isVirtual() ) {
            continue;
        }

        currVertex->computeBetaSpan();
	}
}



void BetaUniverse::setDepthOfWorlds()
{
    setDepthOfVerticesInWorlds();
}



void BetaUniverse::setDepthOfVerticesInWorlds()
{
    setDepthOfVerticesInRootWorld();

    setDepthOfVerticesInSmallWorlds();
}


            
void BetaUniverse::setDepthOfVerticesInRootWorld()
{
    BetaVertex* virtualVertex            = findVirtualVertex();
    BetaCell*   cellIncidentToVirtualVtx = virtualVertex->getFirstCell();

    //  collect cells in root world.
    set<BetaCell*> cellsInRootWorld;

    rg_dList<BetaCell*> cellStack;
    cellStack.pushBack( cellIncidentToVirtualVtx );

	while ( cellStack.getSize() > 0 ) {
		BetaCell* popedCell = cellStack.popFront();
        cellsInRootWorld.insert( popedCell );

        BetaCell* neighbor[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
	    popedCell->getNeighborCells( neighbor );  
		for ( rg_INT i=0; i<4; i++) {

            if ( cellsInRootWorld.find( neighbor[i] ) == cellsInRootWorld.end() ) {
				cellStack.pushBack( neighbor[i] );
            }
		}
    }


    //  set depth of vertices in root world.
    const rg_INT DEPTH_OF_ROOT = 0;
    rg_INT i=0;
    set<BetaCell*>::iterator i_cell;
    for ( i_cell=cellsInRootWorld.begin(); i_cell!=cellsInRootWorld.end(); i_cell++ ) {
        for ( i=0; i<4; i++ ) {
            (*i_cell)->getVertex(i)->setDepth( DEPTH_OF_ROOT );
        }
    }

}



void BetaUniverse::setDepthOfVerticesInSmallWorlds()
{
    //  collect first face in small worlds.
    list< pair<BetaEdge*,BetaFace*> > firstQFacesInWorlds;

    rg_INT i=0;
    for ( i=0; i<m_numFiniteEdges; i++ ) {
        if ( !m_sortedFiniteEdges[i]->isGateToSmallWorlds() ) {
            continue;
        }

        rg_dList<BetaFace*>* facesInSmallWorld = m_sortedFiniteEdges[i]->getSmallWorlds();
        facesInSmallWorld->reset4Loop();
        while ( facesInSmallWorld->setNext4Loop() ) {
            firstQFacesInWorlds.push_back( make_pair(m_sortedFiniteEdges[i], facesInSmallWorld->getEntity() ) );
        }
    }


    //  set depth of vertices in root world.
    while ( !firstQFacesInWorlds.empty() ) {
        BetaEdge* gateEdgeIntoWorld = firstQFacesInWorlds.front().first;
        BetaFace* seedFaceOfWorld   = firstQFacesInWorlds.front().second;

        firstQFacesInWorlds.pop_front();

        BetaVertex* gateVertex[2] = { gateEdgeIntoWorld->getStartVertex(), gateEdgeIntoWorld->getEndVertex()};

        //  define depth of this small world
        rg_INT depthOfThisSmallWorld = rg_UNKNOWN;
        if ( ( gateVertex[0]->getDepth() == rg_UNKNOWN) || (gateVertex[1]->getDepth() == rg_UNKNOWN) ) {
            firstQFacesInWorlds.push_back( make_pair(gateEdgeIntoWorld, seedFaceOfWorld ) );   
            continue;
        }
        else {
            depthOfThisSmallWorld = gateVertex[0]->getDepth();
            if ( depthOfThisSmallWorld < gateVertex[1]->getDepth() ) {
                depthOfThisSmallWorld = gateVertex[1]->getDepth();
            }
            depthOfThisSmallWorld++;
        }

        //  find vertices in this small world
        set<BetaVertex*> verticesInThisSmallWorld;

        if ( seedFaceOfWorld->isIsolatedFace() ) {
            rg_dList<BetaVertex*> vertexInSeedFace;
            seedFaceOfWorld->searchVerticesInWholeWorld( vertexInSeedFace );

            vertexInSeedFace.reset4Loop();
            while ( vertexInSeedFace.setNext4Loop() ) {
                BetaVertex* currVtx = vertexInSeedFace.getEntity();
                if ( (currVtx == gateVertex[0]) || (currVtx == gateVertex[1]) ) {
                    continue;
                }

                verticesInThisSmallWorld.insert( currVtx );
            }
        }
        else {
            set<BetaCell*> cellsInThisSmallWorld;

            rg_dList<BetaCell*> cellStack;
            cellStack.pushBack( seedFaceOfWorld->getLeftCell() );

	        while ( cellStack.getSize() > 0 ) {
		        BetaCell* popedCell = cellStack.popFront();
                cellsInThisSmallWorld.insert( popedCell );

                BetaCell* neighbor[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
	            popedCell->getNeighborCells( neighbor );  
		        for ( rg_INT i=0; i<4; i++) {

                    if ( cellsInThisSmallWorld.find( neighbor[i] ) == cellsInThisSmallWorld.end() ) {
				        cellStack.pushBack( neighbor[i] );
                    }
		        }
            }

            set<BetaCell*>::iterator i_cell;
            for ( i_cell=cellsInThisSmallWorld.begin(); i_cell!=cellsInThisSmallWorld.end(); i_cell++ ) {
                for ( i=0; i<4; i++ ) {
                    BetaVertex* currVtx = (*i_cell)->getVertex(i);
                    if ( (currVtx == gateVertex[0]) || (currVtx == gateVertex[1]) ) {
                        continue;
                    }
                    verticesInThisSmallWorld.insert( currVtx );
                }
            }
        }

        // set depth of vertics in this small world
        set<BetaVertex*>::iterator i_vtx;
        for ( i_vtx=verticesInThisSmallWorld.begin(); i_vtx!=verticesInThisSmallWorld.end(); i_vtx++ ) {
            (*i_vtx)->setDepth( depthOfThisSmallWorld );
        }
    }
}



void BetaUniverse::generateSortedFiniteSimplexes()
{
    generateSortedFiniteCells();
    generateSortedFiniteFaces();
    generateSortedFiniteEdges();
    generateSortedFiniteVertices();
}



int V::GeometryTier::compareBetaSpanOfCells(const void* ts1, const void* ts2)
{
    BetaCell* cell1 = *(BetaCell**) ts1;
    BetaCell* cell2 = *(BetaCell**) ts2;


    if ( cell1->getMinValueForValidSimplex() > cell2->getMinValueForValidSimplex() ) {
        return 1;
    }
    else if ( cell1->getMinValueForValidSimplex() < cell2->getMinValueForValidSimplex() ) {
        return -1;
    }
    else {
        return 0;
    }
}



int V::GeometryTier::compareBetaSpanOfFaces(const void* ts1, const void* ts2)
{
    BetaFace* face1 = *(BetaFace**) ts1;
    BetaFace* face2 = *(BetaFace**) ts2;


    if ( face1->getMinValueForValidSimplex() > face2->getMinValueForValidSimplex() ) {
        return 1;
    }
    else if ( face1->getMinValueForValidSimplex() < face2->getMinValueForValidSimplex() ) {
        return -1;
    }
    else {
        return 0;
    }
}



int V::GeometryTier::compareBetaSpanOfEdges(const void* ts1, const void* ts2)
{
    BetaEdge* edge1 = *(BetaEdge**) ts1;
    BetaEdge* edge2 = *(BetaEdge**) ts2;


    if ( edge1->getMinValueForValidSimplex() > edge2->getMinValueForValidSimplex() ) {
        return 1;
    }
    else if ( edge1->getMinValueForValidSimplex() < edge2->getMinValueForValidSimplex() ) {
        return -1;
    }
    else {
        return 0;
    }
}



int V::GeometryTier::compareBetaSpanOfVertices(const void* ts1, const void* ts2)
{
    BetaVertex* vertex1 = *(BetaVertex**) ts1;
    BetaVertex* vertex2 = *(BetaVertex**) ts2;


    if ( vertex1->getMinValueForValidSimplex() > vertex2->getMinValueForValidSimplex() ) {
        return 1;
    }
    else if ( vertex1->getMinValueForValidSimplex() < vertex2->getMinValueForValidSimplex() ) {
        return -1;
    }
    else {
        return 0;
    }
}



void BetaUniverse::generateSortedFiniteCells()
{
    m_numFiniteCells = 0;

    BetaCell* currCell = rg_NULL;
	m_cells.reset4Loop();
	while( m_cells.setNext4Loop() )  {
		currCell = m_cells.getpEntity();

        if ( currCell->isVirtual() ) {
            continue;
        }
        else {
            m_numFiniteCells++;
        }
	}


    m_sortedFiniteCells = new BetaCell*[m_numFiniteCells];
    rg_INT i = 0;
	m_cells.reset4Loop();
	while( m_cells.setNext4Loop() )  {
		currCell = m_cells.getpEntity();

        if ( currCell->isVirtual() ) {
            continue;
        }
        else {
            m_sortedFiniteCells[i] = currCell;
            i++;
        }
	}

    
    qsort( (void *) m_sortedFiniteCells, m_numFiniteCells, sizeof(BetaCell*), compareBetaSpanOfCells);
}



void BetaUniverse::generateSortedFiniteFaces()
{
    m_numFiniteFaces = 0;

    BetaFace* currFace = rg_NULL;
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		currFace = m_faces.getpEntity();

        if ( currFace->isVirtual() ) {
            continue;
        }
        else {
            m_numFiniteFaces++;
        }
	}

    m_sortedFiniteFaces = new BetaFace*[m_numFiniteFaces];
    rg_INT i = 0;
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		currFace = m_faces.getpEntity();

        if ( currFace->isVirtual() ) {
            continue;
        }
        else {
            m_sortedFiniteFaces[i] = currFace;
            i++;
        }
	}

    
    qsort( (void *) m_sortedFiniteFaces, m_numFiniteFaces, sizeof(BetaFace*), compareBetaSpanOfFaces);
}



void BetaUniverse::generateSortedFiniteEdges()
{
    m_numFiniteEdges = 0;

    BetaEdge* currEdge = rg_NULL;
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		currEdge = m_edges.getpEntity();

        if ( currEdge->isVirtual() ) {
            continue;
        }
        else {
            m_numFiniteEdges++;
        }
	}

    m_sortedFiniteEdges = new BetaEdge*[m_numFiniteEdges];
    rg_INT i = 0;
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		currEdge = m_edges.getpEntity();

        if ( currEdge->isVirtual() ) {
            continue;
        }
        else {
            m_sortedFiniteEdges[i] = currEdge;
            i++;
        }
	}


    qsort( (void *) m_sortedFiniteEdges, m_numFiniteEdges, sizeof(BetaEdge*), compareBetaSpanOfEdges);
}



void BetaUniverse::generateSortedFiniteVertices()
{
    m_numFiniteVertices = 0;

    BetaVertex* currVertex = rg_NULL;
	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		currVertex = m_vertices.getpEntity();

        if ( currVertex->isVirtual() ) {
            continue;
        }
        else {
            m_numFiniteVertices++;
        }
	}

    m_sortedFiniteVertices = new BetaVertex*[m_numFiniteVertices];
    rg_INT i = 0;
	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		currVertex = m_vertices.getpEntity();

        if ( currVertex->isVirtual() ) {
            continue;
        }
        else {
            m_sortedFiniteVertices[i] = currVertex;
            i++;
        }
	}


    qsort( (void *) m_sortedFiniteVertices, m_numFiniteVertices, sizeof(BetaVertex*), compareBetaSpanOfVertices);
}









void BetaUniverse::rearrangeSimplexIDs()
{
    rearrangeCellIDs();
    rearrangeFaceIDs();
    rearrangeEdgeIDs();
    rearrangeVertexIDs();
}



void BetaUniverse::rearrangeCellIDs()
{
	m_cells.reset4Loop();
	while( m_cells.setNext4Loop() )  {
		m_cells.getpEntity()->setID( -1 );
	}

    rg_INT i=0; 
    for ( i=0; i<m_numFiniteCells; i++ ) {
        m_sortedFiniteCells[i]->setID( i );
    }

    rg_INT ID = m_numFiniteCells;
    BetaCell* currCell = rg_NULL;
	m_cells.reset4Loop();
	while( m_cells.setNext4Loop() )  {
		currCell = m_cells.getpEntity();

        if ( currCell->getID() == -1 ) {
            currCell->setID( ID++ );
        }
	}
}



void BetaUniverse::rearrangeFaceIDs()
{
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		m_faces.getpEntity()->setID( -1 );
	}

    rg_INT i=0; 
    for ( i=0; i<m_numFiniteFaces; i++ ) {
        m_sortedFiniteFaces[i]->setID( i );
    }

    rg_INT ID = m_numFiniteFaces;
    BetaFace* currFace = rg_NULL;
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		currFace = m_faces.getpEntity();

        if ( currFace->getID() == -1 ) {
            currFace->setID( ID++ );
        }
	}
}



void BetaUniverse::rearrangeEdgeIDs()
{
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		m_edges.getpEntity()->setID( -1 );
	}

    rg_INT i=0; 
    for ( i=0; i<m_numFiniteEdges; i++ ) {
        m_sortedFiniteEdges[i]->setID( i );
    }

    rg_INT ID = m_numFiniteEdges;
    BetaEdge* currEdge = rg_NULL;
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		currEdge = m_edges.getpEntity();

        if ( currEdge->getID() == -1 ) {
            currEdge->setID( ID++ );
        }
	}
}



void BetaUniverse::rearrangeVertexIDs()
{
	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		m_vertices.getpEntity()->setID( -1 );
	}

    rg_INT i=0; 
    for ( i=0; i<m_numFiniteVertices; i++ ) {
        m_sortedFiniteVertices[i]->setID( i );
    }

    rg_INT ID = m_numFiniteVertices;
    BetaVertex* currVertex = rg_NULL;
	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		currVertex = m_vertices.getpEntity();

        if ( currVertex->getID() == -1 ) {
            currVertex->setID( ID++ );
        }
	}
}




///////////////////////////////////////////////////////////////////////////
//
//  Converter
//
//  1. IWDS -> eIWDS
void BetaUniverse::convertFromIWDSToEIWDS(QuasiTriangulation& qtInIWDS)
{
    BetaCell**   refTableForCells    = createCellsAndMakeRefTableForCells( qtInIWDS );
    BetaVertex** refTableForVertices = createVerticesAndMakeRefTableForVertices( qtInIWDS );

    //  set topology in intra-world.
    setTopologyBetweenCellsAndVertices( qtInIWDS, refTableForCells, refTableForVertices );

    
    createFacesAndSetTopologyBetweenCellsAndFaces( qtInIWDS, refTableForCells );


    createEdgesAndSetTopologyAmongFacesEdgesAndVertices( qtInIWDS, refTableForCells, refTableForVertices );


    BetaFace** refTableForIsolatedFaces = createIsolatedFacesAndMakeRefTableForFaces( qtInIWDS );

    
    //  set topology in inter-world.
    setTopologyBetweenGateEdgesAndSmallWorlds(qtInIWDS, refTableForCells, refTableForIsolatedFaces, refTableForVertices );
    
    
    if ( qtInIWDS.isLinkedVD() ) {
        linkVD( qtInIWDS, refTableForCells, refTableForVertices);
        //checkLinkVD2BC();
    }


    delete [] refTableForIsolatedFaces;
    delete [] refTableForVertices;
    delete [] refTableForCells;
}



BetaCell** BetaUniverse::createCellsAndMakeRefTableForCells(QuasiTriangulation& qtInIWDS)
{
    //  make reference table for cells of quasi-triangulation for eIWDS
    rg_INT     numCells         = qtInIWDS.getNumQTTetrahedron();
    BetaCell** refTableForCells = new BetaCell*[numCells];    
    for ( rg_INT i_cell=0; i_cell<numCells; i_cell++)  {
        refTableForCells[i_cell] = rg_NULL;
    }


    rg_dList<QTTetrahedron>* cellsInIWDS = qtInIWDS.getTetrahedra();

    QTTetrahedron* currCellInIWDS = rg_NULL;
    cellsInIWDS->reset4Loop();
    while ( cellsInIWDS->setNext4Loop() ) {
        currCellInIWDS = cellsInIWDS->getpEntity();

        if ( currCellInIWDS->isAnomaly() == DEGENERACY_ANOMALY ) {
            continue;
        }

        //  create cells of quasi-triangulation for eIWDS
        rg_INT    currCellID = currCellInIWDS->getID();
        BetaCell* currCell   = m_cells.add( BetaCell(currCellID) );
        currCell->setMinTangentSphere( currCellInIWDS->getEmptyTangentSphere() );


        //  set reference table for cells of quasi-triangulation for eIWDS
        refTableForCells[currCellID] = currCell;
    }

    return refTableForCells;
}



BetaVertex** BetaUniverse::createVerticesAndMakeRefTableForVertices(QuasiTriangulation& qtInIWDS)
{
    //  make reference table for cells of quasi-triangulation for eIWDS
    rg_INT       numVertices         = qtInIWDS.getNumQTVertices();
    BetaVertex** refTableForVertices = new BetaVertex*[numVertices];    
    for ( rg_INT i_vtx=0; i_vtx<numVertices; i_vtx++ )  {
        refTableForVertices[i_vtx] = rg_NULL;
    }


    rg_dList<QTVertex>* verticesInIWDS = qtInIWDS.getVertices();

    QTVertex* currVertexInIWDS = rg_NULL;
    verticesInIWDS->reset4Loop();
    while ( verticesInIWDS->setNext4Loop() ) {
        currVertexInIWDS = verticesInIWDS->getpEntity();

        //  create vertices of quasi-triangulation for eIWDS
        rg_INT      vertexID   = qtInIWDS.getAbsoluteIDOfVertex( currVertexInIWDS );

        Ball*       currBall   = rg_NULL;
        if ( !currVertexInIWDS->isVirtual() ) {
            currBall = m_balls.add( Ball(vertexID, currVertexInIWDS->getBall(), 
                                         currVertexInIWDS->getProperty(), currVertexInIWDS->getIDFromInput() ) );
        }
        BetaVertex* currVertex = m_vertices.add( BetaVertex( vertexID, currBall ) );


        //  set reference table for cells of quasi-triangulation for eIWDS
        refTableForVertices[vertexID] = currVertex;
    }

    return refTableForVertices;
}



void BetaUniverse::setTopologyBetweenCellsAndVertices( QuasiTriangulation& qtInIWDS, 
                                                                    BetaCell**     refTableForCells, 
                                                                    BetaVertex**   refTableForVertices )
{
    //  set topology from cell to vertex
    //  set bounding vertices of cell.
    rg_dList<QTTetrahedron>* cellsInIWDS = qtInIWDS.getTetrahedra();

    QTTetrahedron* currCellInIWDS = rg_NULL;
    cellsInIWDS->reset4Loop();
    while ( cellsInIWDS->setNext4Loop() ) {
        currCellInIWDS = cellsInIWDS->getpEntity();

        if ( currCellInIWDS->isAnomaly() == DEGENERACY_ANOMALY ) {
            continue;
        }
    
        rg_INT    cellID   = currCellInIWDS->getID();
        BetaCell* currCell = refTableForCells[ cellID ];

        QTVertex** vertexInIWDS = currCellInIWDS->getVertices();
        for ( rg_INT i_vtx=0; i_vtx<EIWDS_NUM_VERTEX_ON_CELL; i_vtx++ )  {
            rg_INT vertexID = qtInIWDS.getAbsoluteIDOfVertex( vertexInIWDS[i_vtx] );

            //  set bounding vertices of cell
            currCell->setVertex( i_vtx, refTableForVertices[ vertexID ] );
        }
    }


    //  set topology from vertex to cell
    //  set first cell of vertex (first cell is one of cells in the biggest world of vertex).
    rg_dList<QTVertex>* verticesInIWDS = qtInIWDS.getVertices();

    QTVertex* currVertexInIWDS = rg_NULL;
    verticesInIWDS->reset4Loop();
    while ( verticesInIWDS->setNext4Loop() ) {
        currVertexInIWDS = verticesInIWDS->getpEntity();

        rg_INT vertexID = qtInIWDS.getAbsoluteIDOfVertex( currVertexInIWDS );
        rg_INT cellID   = currVertexInIWDS->getFirstTetrahedron()->getID();

        BetaVertex* currVertex = refTableForVertices[ vertexID ];
        BetaCell*   firstCell  = refTableForCells[ cellID ];

        //  set first cell of vertex
        currVertex->setFirstCell( firstCell );
    }
}



void BetaUniverse::createFacesAndSetTopologyBetweenCellsAndFaces( QuasiTriangulation& qtInIWDS, 
                                                                  BetaCell** refTableForCells )
{
    rg_dList<QTTetrahedron>* cellsInIWDS = qtInIWDS.getTetrahedra();
    QTTetrahedron* currCellInIWDS = rg_NULL;
    cellsInIWDS->reset4Loop();
    while ( cellsInIWDS->setNext4Loop() ) {
        currCellInIWDS = cellsInIWDS->getpEntity();

        if ( currCellInIWDS->isAnomaly() == DEGENERACY_ANOMALY ) {
            continue;
        }
    
        rg_INT    cellID   = currCellInIWDS->getID();
        BetaCell* currCell = refTableForCells[ cellID ];

        QTTetrahedron** neighborInIWDS = currCellInIWDS->getNeighbors();
        for ( rg_INT i_cell=0; i_cell<EIWDS_NUM_FACE_ON_CELL; i_cell++ )  {
            rg_INT    neighborID = neighborInIWDS[i_cell]->getID();
            BetaCell* neighbor   = refTableForCells[ neighborID ];

            if ( currCell->getFace(i_cell) != rg_NULL ) {
                continue;
            }

            //  1. create faces of quasi-triangulation for eIWDS
            //  2. set topology from face to cells.
            //     'currCell' and 'neighbor' are the left cell and the right cell of 'currFace', respectively.
            BetaFace* currFace = m_faces.add( BetaFace( m_faces.getSize(), currCell, neighbor ) );

            
            //  3. set topology from cell to face. 
            currCell->setFace( i_cell, currFace );

            rg_INT matePosInNeighbor = currCellInIWDS->getMatePosInNeighbor(neighborInIWDS[i_cell], i_cell);
            neighbor->setFace( matePosInNeighbor, currFace );
        }
    }
}



void BetaUniverse::createEdgesAndSetTopologyAmongFacesEdgesAndVertices( QuasiTriangulation& qtInIWDS, 
                                                                                     BetaCell**     refTableForCells,
                                                                                     BetaVertex**   refTableForVertices )
{
    rg_dList<QTTetrahedron>* cellsInIWDS = qtInIWDS.getTetrahedra();

    QTTetrahedron* currCellInIWDS = rg_NULL;
    cellsInIWDS->reset4Loop();
    while ( cellsInIWDS->setNext4Loop() ) {
        currCellInIWDS = cellsInIWDS->getpEntity();

        if ( currCellInIWDS->isAnomaly() == DEGENERACY_ANOMALY ) {
            continue;
        }
    
        rg_INT    cellID   = currCellInIWDS->getID();
        BetaCell* currCell = refTableForCells[ cellID ];

        QTTetrahedron** neighborInIWDS = currCellInIWDS->getNeighbors();
        for ( rg_INT i_cell=0; i_cell<EIWDS_NUM_FACE_ON_CELL; i_cell++ )  {
            BetaFace* currFace = currCell->getFace( i_cell );
            if ( currFace->isVisited() ) {
                continue;
            }


            rg_INT    neighborID = neighborInIWDS[i_cell]->getID();
            BetaCell* neighbor   = refTableForCells[ neighborID ];

            QTFace    currFaceInIWDS(currCellInIWDS, currCellInIWDS->getVertex(i_cell));

            if ( currCell == currFace->getLeftCell() ) {
                currFaceInIWDS = currFaceInIWDS.getTwinFace();
            }

            rg_dList<QTVertex*> vertexList;
            qtInIWDS.inquireVerticesOfFaceInIntraWorld( currFaceInIWDS, vertexList);
            

            rg_INT i_vtx = 0;
            QTVertex*   vertexInIWDS[EIWDS_NUM_EDGE_ON_FACE+1];
            BetaVertex* vertex[EIWDS_NUM_EDGE_ON_FACE+1];
            vertexList.reset4Loop();
            while ( vertexList.setNext4Loop() ) {
                vertexInIWDS[i_vtx] = vertexList.getEntity();
                rg_INT vertexID     = qtInIWDS.getAbsoluteIDOfVertex( vertexInIWDS[i_vtx] );
                vertex[i_vtx]       = refTableForVertices[ vertexID ];
                i_vtx++;
            }
            vertexInIWDS[EIWDS_NUM_EDGE_ON_FACE] = vertexInIWDS[0];
            vertex[EIWDS_NUM_EDGE_ON_FACE]       = vertex[0];


            for ( i_vtx=0; i_vtx<EIWDS_NUM_EDGE_ON_FACE; i_vtx++ ) {
                BetaEdge* currEdge = currFace->findEdge( vertex[i_vtx], vertex[i_vtx+1] ); 

                if ( currEdge != rg_NULL ) {
                    continue;
                }


                //  1. create edge.
                //  2. set topology from edge to two extreme vertices.
                //  3. set topology from edge to face.
                currEdge = m_edges.add( BetaEdge( m_edges.getSize(), vertex[i_vtx], vertex[i_vtx+1], currFace ) );


                //  4. set topology from face to edge.
                QTEdge currEdgeInIWDS( currCellInIWDS, vertexInIWDS[i_vtx], vertexInIWDS[i_vtx+1] );

                setTopologyFromFaceToEdge( qtInIWDS, refTableForCells, refTableForVertices, currEdgeInIWDS, currEdge );
            }

            currFace->isVisited(rg_TRUE);
        }

        //currCell->isVisited( rg_TRUE );
    }

    m_faces.reset4Loop();
    while ( m_faces.setNext4Loop() ) {
        m_faces.getpEntity()->isVisited( rg_FALSE );
    }


    //  5. set topology from vertex to edge.
    setTopologyFromVerticesToEdges();
}




void BetaUniverse::setTopologyFromFaceToEdge( QuasiTriangulation& qtInIWDS, 
                                                           BetaCell**          refTableForCells,
                                                           BetaVertex**        refTableForVertices,
                                                           const QTEdge&       edgeInIWDS,
                                                           BetaEdge*           edgeInEIWDS )
{
    //  4. set topology from face to edge.
    QTVertex* startVertexInIWDS = edgeInIWDS.getStartVertex();
    QTVertex* endVertexInIWDS   = edgeInIWDS.getEndVertex();


    rg_dList<QTFace> faceListIncidentToEdgeInIWDS;
    qtInIWDS.inquireFacesOfEdgeInIntraWorld( edgeInIWDS, faceListIncidentToEdgeInIWDS );

    faceListIncidentToEdgeInIWDS.reset4Loop();
    while ( faceListIncidentToEdgeInIWDS.setNext4Loop() ) {
        QTFace incidentFaceInIWDS = faceListIncidentToEdgeInIWDS.getEntity();

//         QTTetrahedron* incidentCellInIWDS[2] = { incidentFaceInIWDS.getTetrahedron(),
//                                                  incidentFaceInIWDS.getIncidentTetrahedron() };
//         BetaCell*      incidentCell[2]       = { refTableForCells[ incidentCellInIWDS[0]->getID() ], 
//                                                  refTableForCells[ incidentCellInIWDS[1]->getID() ] };
// 
//         BetaFace*      incidentFace = incidentCell[0]->findSharingFaceWith( incidentCell[1] );
        QTTetrahedron* incidentCellInIWDS = incidentFaceInIWDS.getTetrahedron();
        QTVertex*      mateVertexInIWDS   = incidentFaceInIWDS.getMateVertex();

        BetaCell*      incidentCell = refTableForCells[ incidentCellInIWDS->getID() ];
        BetaVertex*    mateVertex   = refTableForVertices[ qtInIWDS.getAbsoluteIDOfVertex(mateVertexInIWDS) ];
        BetaFace*      incidentFace = incidentCell->getMateFace(mateVertex);

        //  'vertexOnFaceInIWDS' is a list of vertices to bound 'incidentFaceInIWDS'.
        //  the vertices in 'vertexOnFaceInIWDS' is arranged by the CCW orientation 
        //  viewing from the right cell of 'incidentFace' which is a face in eIWDS corresponding to 'incidentFaceInIWDS'.
        rg_dList<QTVertex*> vertexOnFaceInIWDS;
        if ( incidentCell == incidentFace->getRightCell() ) {
            qtInIWDS.inquireVerticesOfFaceInIntraWorld( incidentFaceInIWDS, vertexOnFaceInIWDS);
        }
        else {
            qtInIWDS.inquireVerticesOfFaceInIntraWorld( incidentFaceInIWDS.getTwinFace(), vertexOnFaceInIWDS);
        }


        QTVertex** vertexArray = vertexOnFaceInIWDS.getArray();
        if ( vertexArray[0] == startVertexInIWDS ) {
            //  edge orientation is forward 
            if ( vertexArray[1] == endVertexInIWDS ) {
                incidentFace->setEdge( 0, rg_TRUE, edgeInEIWDS );
            }

            //  edge orientation is backward 
            if ( vertexArray[2] == endVertexInIWDS ) {
                incidentFace->setEdge( 2, rg_FALSE, edgeInEIWDS );
            }
        }
        else if ( vertexArray[1] == startVertexInIWDS ) {
            //  edge orientation is forward 
            if ( vertexArray[2] == endVertexInIWDS ) {
                incidentFace->setEdge( 1, rg_TRUE, edgeInEIWDS );
            }

            //  edge orientation is backward 
            if ( vertexArray[0] == endVertexInIWDS ) {
                incidentFace->setEdge( 0, rg_FALSE, edgeInEIWDS );
            }
        }
        else { // if ( vertexArray[2] == startVertexInIWDS ) {
            //  edge orientation is forward 
            if ( vertexArray[0] == endVertexInIWDS ) {
                incidentFace->setEdge( 2, rg_TRUE, edgeInEIWDS );
            }

            //  edge orientation is backward 
            if ( vertexArray[1] == endVertexInIWDS ) {
                incidentFace->setEdge( 1, rg_FALSE, edgeInEIWDS );
            }
        }

        delete [] vertexArray;

    }

}



void BetaUniverse::setTopologyFromVerticesToEdges()
{
    BetaVertex* currVertex = rg_NULL;

    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() ) {
        currVertex = m_vertices.getpEntity();

        BetaCell* firstCell = currVertex->getFirstCell();
        //  the vertex contributes to an isolated triangular face.
        if ( firstCell == rg_NULL ) {
            continue;
        }

        BetaEdge* incidentEdge[3];
        firstCell->findEdgesIncidentToVertex( currVertex, incidentEdge );

        currVertex->setFirstEdge( incidentEdge[0] );
    }
}



BetaFace** BetaUniverse::createIsolatedFacesAndMakeRefTableForFaces( QuasiTriangulation& qtInIWDS )
{
    //  make reference table for cells of quasi-triangulation for eIWDS
    rg_INT numCells = qtInIWDS.getNumQTTetrahedron();
    BetaFace** refTableForIsolatedFaces   = new BetaFace*[numCells];    
    for ( rg_INT i_cell=0; i_cell<numCells; i_cell++)  {
        refTableForIsolatedFaces[i_cell] = rg_NULL;
    }


    rg_dList<QTTetrahedron>* cellsInIWDS = qtInIWDS.getTetrahedra();
    QTTetrahedron* currCellInIWDS = rg_NULL;
    cellsInIWDS->reset4Loop();
    while ( cellsInIWDS->setNext4Loop() ) {
        currCellInIWDS = cellsInIWDS->getpEntity();

        if ( currCellInIWDS->isAnomaly() != DEGENERACY_ANOMALY ) {
            continue;
        }

        //  create isolated faces of quasi-triangulation for eIWDS
        BetaFace* currFace = m_faces.add( BetaFace( m_faces.getSize() ) );

        //  set reference table for isolated faces of quasi-triangulation for eIWDS
        rg_INT currCellID = currCellInIWDS->getID();

        refTableForIsolatedFaces[currCellID] = currFace;
    }

    return refTableForIsolatedFaces;
}



void BetaUniverse::setTopologyBetweenGateEdgesAndSmallWorlds( QuasiTriangulation& qtInIWDS, 
                                                                           BetaCell**     refTableForCells,
                                                                           BetaFace**     refTableForIsolatedFaces,
                                                                           BetaVertex**   refTableForVertices )
{
    rg_dList<QTGate*> ptrGatesInIWDS;
    rg_dList<QTGate>* gatesInIWDS = qtInIWDS.getGates()->getGates();

    QTGate* currGateInIWDS = rg_NULL;
    gatesInIWDS->reset4Loop();
    while ( gatesInIWDS->setNext4Loop() )  {
        ptrGatesInIWDS.add( gatesInIWDS->getpEntity() );
    }


    rg_dList<BetaEdge*> edgesToBeRemoved;


    ptrGatesInIWDS.reset4Loop();
    while ( ptrGatesInIWDS.getSize() > 0 ) {
        QTGate* currGateInIWDS = ptrGatesInIWDS.popFront();

        // 1. set topology between a gate edge and its biggest world.
        QTTetrahedron*   cellInBigWorldOfIWDS  = currGateInIWDS->getBigTetrahedron();
        QTVertex*        vertexOfGateAtIWDS[2] = { currGateInIWDS->getVertex(0), currGateInIWDS->getVertex(1) };
        BetaVertex*      vertexOfGate[2]       = { refTableForVertices[ qtInIWDS.getAbsoluteIDOfVertex(vertexOfGateAtIWDS[0]) ], 
                                                   refTableForVertices[ qtInIWDS.getAbsoluteIDOfVertex(vertexOfGateAtIWDS[1]) ] };


        BetaEdge* gateEdge = rg_NULL;
        if ( cellInBigWorldOfIWDS->isAnomaly() != DEGENERACY_ANOMALY ) {
            // 1.1    The biggest world consists of more than two cells.
            // 1.1.1  find edge in eIWDS
            BetaCell*   cellInBigWorld = refTableForCells[ cellInBigWorldOfIWDS->getID() ];
            gateEdge  = cellInBigWorld->findEdge(vertexOfGate[0], vertexOfGate[1]); 

            // 1.1.2  set topology from edge to face ( set the first face of edge )
            BetaFace* faceInBigWorld[2] = {rg_NULL, rg_NULL};
            cellInBigWorld->findFacesIncidentToEdge(gateEdge, faceInBigWorld);

            gateEdge->setFirstFace( faceInBigWorld[0] );
        }
        else {
            // 1.2    The biggest world consists of a single isolated triangular face.
            // 1.2.1  find edge in eIWDS
            BetaFace*   faceInBigWorld = refTableForIsolatedFaces[ cellInBigWorldOfIWDS->getID() ];
            gateEdge  = faceInBigWorld->findEdge(vertexOfGate[0], vertexOfGate[1]);
            
            // 1.2.2  set topology from edge to face ( set the first face of edge )
            gateEdge->setFirstFace( faceInBigWorld );
        }
        
        
        // 2. set topology between a gate edge and its small worlds.
        rg_dList<QTTetrahedron*>* smallWorldsInIWDS = currGateInIWDS->getSmallTetrahedra();

        QTTetrahedron* cellInSmallWorldOfIWDS = rg_NULL;
        smallWorldsInIWDS->reset4Loop();
        while ( smallWorldsInIWDS->setNext4Loop() ) {
            cellInSmallWorldOfIWDS = smallWorldsInIWDS->getEntity();

            if ( cellInSmallWorldOfIWDS->isAnomaly() != DEGENERACY_ANOMALY ) {
                // set topology between a gate edge and a small world which consists of more than two cells.
                BetaEdge* gateEdgeInSmallWorld 
                     = setTopologyBetweenGateEdgeAndSmallWorldAsSetOfCells( 
                                                                vertexOfGate, gateEdge, cellInSmallWorldOfIWDS,
                                                                qtInIWDS, refTableForCells );


                edgesToBeRemoved.add( gateEdgeInSmallWorld );
            }
            else {
                // set topology between a gate edge and a small world which consists of an isolated face.
                setTopologyBetweenGateEdgeAndSmallWorldAsIsolatedFace( 
                                                                gateEdge, cellInSmallWorldOfIWDS, qtInIWDS, 
                                                                refTableForIsolatedFaces, refTableForVertices );

            }
        }
    }

    
    //  remove edges in small world corresponding to gate edges.
    edgesToBeRemoved.reset4Loop();
    while ( edgesToBeRemoved.setNext4Loop() ) {
        edgesToBeRemoved.getEntity()->isVisited( rg_TRUE );
    }

    rg_INT count = 0;
    rg_INT numEdgeToBeRemoved = edgesToBeRemoved.getSize();
    rg_dNode<BetaEdge>* currNode = m_edges.getFirstpNode();
    while ( count < numEdgeToBeRemoved ) {
        BetaEdge* currEdge = currNode->getpEntity();

        if ( currEdge->isVisited() )  {
            rg_dNode<BetaEdge>* nodeToBeRemoved = currNode;

            currNode = currNode->getPrev();

            m_edges.kill( nodeToBeRemoved );
            count++;
        }

        currNode = currNode->getNext();
    }
}



BetaEdge* BetaUniverse::setTopologyBetweenGateEdgeAndSmallWorldAsSetOfCells(
                                                                      BetaVertex**        vertexOfGate,
                                                                      BetaEdge*           gateEdge,
                                                                      QTTetrahedron*      cellInSmallWorldOfIWDS,
                                                                      QuasiTriangulation& qtInIWDS, 
                                                                      BetaCell**          refTableForCells )
{
    //  This small world consists of more than two cells 
    //  and its topology is completely connected.
    //  However, the gate edge in the small world is distinguished from the gate edge in big world.
    //  Hence, we will remove the gate edge in the small world which is returned by this function.


    //  1. find an edge(small gate edge) corresponding to the gate edge.
    BetaCell* cellInSmallWorld     = refTableForCells[ cellInSmallWorldOfIWDS->getID() ];
    BetaEdge* gateEdgeInSmallWorld = cellInSmallWorld->findEdge(vertexOfGate[0], vertexOfGate[1]); 


    //  2. replace the small gate edge with the gate edge in big world.
    BetaFace* faceIncidentToEdgeInCell[2] = {rg_NULL, rg_NULL};
    cellInSmallWorld->findFacesIncidentToEdge(gateEdgeInSmallWorld, faceIncidentToEdgeInCell);

    rg_dList<BetaFace*> faceaIncidentGateEdgeInSmallWorld;
    gateEdgeInSmallWorld->searchFacesInIntraWorld( faceaIncidentGateEdgeInSmallWorld, 
                                                   faceIncidentToEdgeInCell[0]        );

    if ( gateEdge->getStartVertex() == gateEdgeInSmallWorld->getStartVertex() ) {
        BetaFace* currFace = rg_NULL;
        faceaIncidentGateEdgeInSmallWorld.reset4Loop();
        while ( faceaIncidentGateEdgeInSmallWorld.setNext4Loop() ) {
            currFace = faceaIncidentGateEdgeInSmallWorld.getEntity();
        
            rg_INT posEdge = currFace->getPosOfEdge(gateEdgeInSmallWorld);

            currFace->setEdge( posEdge, gateEdge);
        }
    }
    else {
        BetaFace* currFace = rg_NULL;
        faceaIncidentGateEdgeInSmallWorld.reset4Loop();
        while ( faceaIncidentGateEdgeInSmallWorld.setNext4Loop() ) {
            currFace = faceaIncidentGateEdgeInSmallWorld.getEntity();
        
            rg_INT  posEdge = currFace->getPosOfEdge(gateEdgeInSmallWorld);
            rg_BOOL orientation = currFace->getEdgeOrientation( posEdge );
            currFace->setEdge( posEdge, !orientation, gateEdge);
        }
    }

    gateEdge->addSmallWorld( faceaIncidentGateEdgeInSmallWorld.getFirstEntity() );

    return gateEdgeInSmallWorld;
}



void BetaUniverse::setTopologyBetweenGateEdgeAndSmallWorldAsIsolatedFace( 
                                                                      BetaEdge*      gateEdge,
                                                                      QTTetrahedron*      cellInSmallWorldOfIWDS,
                                                                      QuasiTriangulation& qtInIWDS, 
                                                                      BetaFace**     refTableForIsolatedFaces,
                                                                      BetaVertex**   refTableForVertices )
{
    BetaFace* isolatedfaceInSmallWorld = refTableForIsolatedFaces[ cellInSmallWorldOfIWDS->getID() ];


    // find three vertices in eIWDS to bound the isolated face.
    rg_dList<QTVertex*> vertextOnIsolatedFaceInIWDS;
    qtInIWDS.inquireVerticesOfCellInIntraWorld(cellInSmallWorldOfIWDS, vertextOnIsolatedFaceInIWDS);

    rg_INT i_vtx = 0;
    BetaVertex* vertexOnIsolatedFace[3];
    vertextOnIsolatedFaceInIWDS.reset4Loop();
    while ( vertextOnIsolatedFaceInIWDS.setNext4Loop() ) {
        QTVertex* currVtxInIWDS = vertextOnIsolatedFaceInIWDS.getEntity();
        rg_INT    vtxID         = qtInIWDS.getAbsoluteIDOfVertex( currVtxInIWDS );

        vertexOnIsolatedFace[i_vtx] = refTableForVertices[ vtxID ];
        i_vtx++;
    }
    

    //  So far as, we do not define bounding edges of an isolated face.
    //  1. We set topology between a gate edge and isolated face.
    //  2. we create two edges on the isolated face 
    //     and set topology between the edges and isolated face.

    BetaVertex* startVtxOfGateEdge = gateEdge->getStartVertex();
    BetaVertex* endVtxOfGateEdge   = gateEdge->getEndVertex();

    if ( vertexOnIsolatedFace[0] == startVtxOfGateEdge ) {
        //  the orientation of the gate edge is FORWARD against the orientation of the isolated face.
        if ( vertexOnIsolatedFace[1] == endVtxOfGateEdge ) {
            // create two edges on the isolated face 
            // and set topology from edge to isolated face.
            BetaEdge* edgeOfIsolatedFace[2] = { rg_NULL, rg_NULL };
            edgeOfIsolatedFace[0] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[1], vertexOnIsolatedFace[2],
                                                                  isolatedfaceInSmallWorld ) );
            edgeOfIsolatedFace[1] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[2], vertexOnIsolatedFace[0],
                                                                  isolatedfaceInSmallWorld ) );

            // set topology from an isolated face to edges.
            isolatedfaceInSmallWorld->setEdge( 0, rg_TRUE, gateEdge );
            isolatedfaceInSmallWorld->setEdge( 1, rg_TRUE, edgeOfIsolatedFace[0] );
            isolatedfaceInSmallWorld->setEdge( 2, rg_TRUE, edgeOfIsolatedFace[1] );

            // set topology from vertex to edge.
            vertexOnIsolatedFace[2]->setFirstEdge( edgeOfIsolatedFace[0] );
        }
        //  the orientation of the gate edge is BACKWARD against the orientation of the isolated face.
        else { // if ( vertexOnIsolatedFace[2] == endVtxOfGateEdge ) {
            // create two edges on the isolated face 
            // and set topology from edge to isolated face.
            BetaEdge* edgeOfIsolatedFace[2] = { rg_NULL, rg_NULL };
            edgeOfIsolatedFace[0] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[0], vertexOnIsolatedFace[1],
                                                                  isolatedfaceInSmallWorld ) );
            edgeOfIsolatedFace[1] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[1], vertexOnIsolatedFace[2],
                                                                  isolatedfaceInSmallWorld ) );

            // set topology from an isolated face to edges.
            isolatedfaceInSmallWorld->setEdge( 0, rg_TRUE, edgeOfIsolatedFace[0] );
            isolatedfaceInSmallWorld->setEdge( 1, rg_TRUE, edgeOfIsolatedFace[1] );
            isolatedfaceInSmallWorld->setEdge( 2, rg_FALSE, gateEdge );

            // set topology from vertex to edge.
            vertexOnIsolatedFace[1]->setFirstEdge( edgeOfIsolatedFace[0] );
        }
    }
    else if ( vertexOnIsolatedFace[1] == startVtxOfGateEdge ) {
        if ( vertexOnIsolatedFace[2] == endVtxOfGateEdge ) {

            BetaEdge* edgeOfIsolatedFace[2] = { rg_NULL, rg_NULL };
            edgeOfIsolatedFace[0] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[0], vertexOnIsolatedFace[1],
                                                                  isolatedfaceInSmallWorld ) );
            edgeOfIsolatedFace[1] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[2], vertexOnIsolatedFace[0],
                                                                  isolatedfaceInSmallWorld ) );

            isolatedfaceInSmallWorld->setEdge( 0, rg_TRUE, edgeOfIsolatedFace[0] );
            isolatedfaceInSmallWorld->setEdge( 1, rg_TRUE, gateEdge );
            isolatedfaceInSmallWorld->setEdge( 2, rg_TRUE, edgeOfIsolatedFace[1] );

            // set topology from vertex to edge.
            vertexOnIsolatedFace[0]->setFirstEdge( edgeOfIsolatedFace[0] );
        }
        else { // if ( vertexOnIsolatedFace[0] == endVtxOfGateEdge ) {

            BetaEdge* edgeOfIsolatedFace[2] = { rg_NULL, rg_NULL };
            edgeOfIsolatedFace[0] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[1], vertexOnIsolatedFace[2],
                                                                  isolatedfaceInSmallWorld ) );
            edgeOfIsolatedFace[1] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[2], vertexOnIsolatedFace[0],
                                                                  isolatedfaceInSmallWorld ) );

            isolatedfaceInSmallWorld->setEdge( 0, rg_FALSE, gateEdge );
            isolatedfaceInSmallWorld->setEdge( 1, rg_TRUE, edgeOfIsolatedFace[0] );
            isolatedfaceInSmallWorld->setEdge( 2, rg_TRUE, edgeOfIsolatedFace[1] );

            // set topology from vertex to edge.
            vertexOnIsolatedFace[2]->setFirstEdge( edgeOfIsolatedFace[0] );
        }
    }
    else { // if ( vertexOnIsolatedFace[2] == startVtxOfGateEdge ) {
        if ( vertexOnIsolatedFace[0] == endVtxOfGateEdge ) {

            BetaEdge* edgeOfIsolatedFace[2] = { rg_NULL, rg_NULL };
            edgeOfIsolatedFace[0] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[0], vertexOnIsolatedFace[1],
                                                                  isolatedfaceInSmallWorld ) );
            edgeOfIsolatedFace[1] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[1], vertexOnIsolatedFace[2],
                                                                  isolatedfaceInSmallWorld ) );

            isolatedfaceInSmallWorld->setEdge( 0, rg_TRUE, edgeOfIsolatedFace[0] );
            isolatedfaceInSmallWorld->setEdge( 1, rg_TRUE, edgeOfIsolatedFace[1] );
            isolatedfaceInSmallWorld->setEdge( 2, rg_TRUE, gateEdge );

            // set topology from vertex to edge.
            vertexOnIsolatedFace[1]->setFirstEdge( edgeOfIsolatedFace[0] );
        }
        else { // if ( vertexOnIsolatedFace[1] == endVtxOfGateEdge ) {

            BetaEdge* edgeOfIsolatedFace[2] = { rg_NULL, rg_NULL };
            edgeOfIsolatedFace[0] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[0], vertexOnIsolatedFace[1],
                                                                  isolatedfaceInSmallWorld ) );
            edgeOfIsolatedFace[1] = m_edges.add( BetaEdge( m_edges.getSize(), 
                                                                  vertexOnIsolatedFace[2], vertexOnIsolatedFace[0],
                                                                  isolatedfaceInSmallWorld ) );

            isolatedfaceInSmallWorld->setEdge( 0, rg_TRUE, edgeOfIsolatedFace[0] );
            isolatedfaceInSmallWorld->setEdge( 1, rg_FALSE, gateEdge );
            isolatedfaceInSmallWorld->setEdge( 2, rg_TRUE, edgeOfIsolatedFace[1] );

            // set topology from vertex to edge.
            vertexOnIsolatedFace[0]->setFirstEdge( edgeOfIsolatedFace[0] );
        }
    }
    
    gateEdge->addSmallWorld( isolatedfaceInSmallWorld );
}



void BetaUniverse::linkVD( QuasiTriangulation& qtInIWDS, BetaCell** refTableForCells, BetaVertex** refTableForVertices)
{
    linkVVertexAndBetaCell( qtInIWDS, refTableForCells );
    linkVCellAndBetaVertex( qtInIWDS, refTableForVertices );
    linkVEdgeAndBetaFace();
    linkVFaceAndBetaEdge();
}



void BetaUniverse::linkVVertexAndBetaCell( QuasiTriangulation& qtInIWDS, BetaCell** refTableForCells )
{
    rg_dList<QTTetrahedron>* q_cells = qtInIWDS.getTetrahedra();

    q_cells->reset4Loop();
    while ( q_cells->setNext4Loop() ) {
        QTTetrahedron* currQCell = q_cells->getpEntity();

        if ( currQCell->isAnomaly() == DEGENERACY_ANOMALY ) {
            continue;
        }

        VDVertex* currVVtx     = currQCell->getVVertex();
        BetaCell* currBetaCell = refTableForCells[currQCell->getID()];

        currVVtx->connectBetaCell( currBetaCell );
        currBetaCell->connectVVertex( currVVtx );
    }
}



void BetaUniverse::linkVCellAndBetaVertex( QuasiTriangulation& qtInIWDS, BetaVertex** refTableForVertices )
{
    rg_dList<QTVertex>* q_vertices = qtInIWDS.getVertices();

    q_vertices->reset4Loop();
    while ( q_vertices->setNext4Loop() ) {
        QTVertex* currQVtx = q_vertices->getpEntity();

        VDCell*     currVCell   = currQVtx->getVCell();

        rg_INT      vertexID    = qtInIWDS.getAbsoluteIDOfVertex( currQVtx );
        BetaVertex* currBetaVtx = refTableForVertices[vertexID];

        currVCell->connectBetaVertex( currBetaVtx );
        currBetaVtx->connectVCell( currVCell );
    }
}


void BetaUniverse::linkVEdgeAndBetaFace()
{
    //  we already know the corresponding v-cell of each b-vertex.
    //  we find the corresponding v-edge of each b-face.
    //  Note that in the following codes, only v-edges with v-vtx are found.
    m_cells.reset4Loop();
    while ( m_cells.setNext4Loop() ) {
        BetaCell* currBetaCell = m_cells.getpEntity();
        VDVertex* currVDVtx    = currBetaCell->getVVertex();

        VDEdge*   v_edge[4]    = { rg_NULL, rg_NULL, rg_NULL, rg_NULL };
        VDCell*   v_cell[4][2] = { {rg_NULL, rg_NULL}, {rg_NULL, rg_NULL}, {rg_NULL, rg_NULL}, {rg_NULL, rg_NULL} };
        for ( rg_INT i=0; i<4; i++ ) {
            v_edge[i]    = currVDVtx->getIncidentEdge( i );
            v_cell[i][0] = v_edge[i]->getStartVCell();
            v_cell[i][1] = v_edge[i]->getEndVCell();

        }

        for ( rg_INT i=0; i<4; i++ ) {
            BetaFace* currBetaFace = currBetaCell->getFace(i);

            if ( currBetaFace->getVEdge() != rg_NULL ) {
                continue;
            }

            BetaVertex* startVtx = currBetaCell->getVertex(i);
            BetaVertex* endVtx   = currBetaCell->getNeighborCell(i)->getMateVertex( currBetaFace );

            for ( rg_INT j=0; j<4; j++ ) {
                if ( v_edge[j]->getBetaFace() != rg_NULL ) {
                    continue;
                }

                if (   (startVtx->getVCell() == v_cell[j][0] && endVtx->getVCell() == v_cell[j][1] )
                    || (startVtx->getVCell() == v_cell[j][1] && endVtx->getVCell() == v_cell[j][0] ) ) {
                    currBetaFace->connectVEdge( v_edge[j] );
                    v_edge[j]->connectBetaFace( currBetaFace );
                    break;
                }
            }
        }

    }



    //  Find isolated b-face corresponds to v-edge without v-vtx.
    rg_INT i=0;
    m_faces.reset4Loop();
    while ( m_faces.setNext4Loop() ) {
        BetaFace* currBetaFace = m_faces.getpEntity();

        if ( !currBetaFace->isIsolatedFace() ) {
            continue;
        }


        // Find gate edge whose small world contains the current isolated face.
        BetaEdge* gateEdgeToCurrFace = rg_NULL;
        for ( i=0; i<3; i++ ) {
            BetaEdge* currBetaEdge = currBetaFace->getEdge(i);

            if ( currBetaEdge->isGateToSmallWorlds() ) {
                if ( currBetaEdge->getSmallWorlds()->isInList( currBetaFace ) == rg_TRUE ) {
                    gateEdgeToCurrFace = currBetaEdge;
                    break;
                }
            }
        }


        //  make set of v-cells corresponding to b-vtx of the current isolated face.
        set<VDCell*>          setOfvcellForCurrBetaFace;
        rg_dList<BetaVertex*> boundingVtx;
        currBetaFace->searchVerticesInWholeWorld( boundingVtx );

        boundingVtx.reset4Loop();
        while ( boundingVtx.setNext4Loop() ) {
            setOfvcellForCurrBetaFace.insert( boundingVtx.getEntity()->getVCell() );
        }


        //  find v-faces corresponding to gate-edge.
        //  Note that there may be two more than v-faces defined by two v-cells
        VDCell* vcellATStartGateEdge = gateEdgeToCurrFace->getStartVertex()->getVCell(); 
        VDCell* vcellATEndGateEdge   = gateEdgeToCurrFace->getEndVertex()->getVCell(); 

        rg_dList<VDFace*> sharingVFaces;
        vcellATStartGateEdge->findFacesToShareWith(vcellATEndGateEdge, sharingVFaces);


        //  find v-edge corresponding to current b-face.
        sharingVFaces.reset4Loop();
        while ( sharingVFaces.setNext4Loop() ) {
            VDFace* currVFace = sharingVFaces.getEntity();

            rg_dList<VDEdge*> edgeList;
            currVFace->inquireBoundingEdges( edgeList );

            edgeList.reset4Loop();
            while ( edgeList.setNext4Loop() ) {
                VDEdge* currVEdge = edgeList.getEntity();

                if ( currVEdge->getStartVertex() != rg_NULL ) {
                    continue;
                }


                //  check for whether three v-cells to define v-edge is identical v-cells to define current b-face.
                rg_dList<VDCell*> cellList;
                currVEdge->inquireIncidentCells( cellList );

                rg_BOOL isCorrespondingVEdge = rg_TRUE;
                cellList.reset4Loop();
                while ( cellList.setNext4Loop() ) {
                    VDCell* currVCell = cellList.getEntity();

                    if ( setOfvcellForCurrBetaFace.find( currVCell ) == setOfvcellForCurrBetaFace.end() ) {
                        isCorrespondingVEdge = rg_FALSE;
                        break;
                    }
                }


                //  if identical, connect v-edge with b-face.
                if ( isCorrespondingVEdge ) {
                    currBetaFace->connectVEdge( currVEdge );
                    currVEdge->connectBetaFace( currBetaFace );
                    break;
                }
            }
        }
    }
}



void BetaUniverse::linkVFaceAndBetaEdge()
{
    m_edges.reset4Loop();
    while ( m_edges.setNext4Loop() ) {
        BetaEdge* currBetaEdge = m_edges.getpEntity();

        BetaFace* incidentBetaFace     = currBetaEdge->getFirstFace();
        VDEdge*   VEdgeForInciBetaFace = incidentBetaFace->getVEdge();


        //  find v-faces corresponding to gate-edge.
        //  Note that there may be two more than v-faces defined by two v-cells
        VDCell* vcellAtStartVtx = currBetaEdge->getStartVertex()->getVCell(); 
        VDCell* vcellAtEndVtx   = currBetaEdge->getEndVertex()->getVCell(); 

        rg_dList<VDFace*> sharingVFaces;
        vcellAtStartVtx->findFacesToShareWith(vcellAtEndVtx, sharingVFaces);

        if ( sharingVFaces.getSize() == 1 ) {
            VDFace* currVFace = sharingVFaces.getFirstEntity();
            currBetaEdge->connectVFace( currVFace );
            currVFace->connectBetaEdge( currBetaEdge );
        }
        else {
            sharingVFaces.reset4Loop();
            while ( sharingVFaces.setNext4Loop() ) {
                VDFace* currVFace = sharingVFaces.getEntity();

                if ( currVFace->getBetaEdge() != rg_NULL ) {
                    continue;
                }

                if ( VEdgeForInciBetaFace->isIncidentTo( currVFace ) ) {
                    currBetaEdge->connectVFace( currVFace );
                    currVFace->connectBetaEdge( currBetaEdge );
                }
            }
        }
    }
}




void BetaUniverse::checkLinkVD2BC()
{
    ofstream fout("check_linkVD2BC.txt");

    m_cells.reset4Loop();
    while ( m_cells.setNext4Loop() ) {
        BetaCell* currCell = m_cells.getpEntity();

        VDVertex* currVVtx = currCell->getVVertex();

        if ( currVVtx == rg_NULL ) {
            fout << "CELL\t" << currCell->getID() << "\t" << "no V-vtx" << endl;            
        }
        else {
            Sphere minTangentSphereOnBCell = currCell->getMinTangentSphere();
            Sphere minTangentSphereOnVVtx( currVVtx->getPoint(), currVVtx->getRadiusOfTangentSphere() );
            if ( !(minTangentSphereOnBCell == minTangentSphereOnVVtx) )  {
                fout << "CELL\t" << currCell->getID() << "\t" << "diff. min. tangent sphere." << endl;            
            }
        }
    }


    m_faces.reset4Loop();
    while ( m_faces.setNext4Loop() ) {
        BetaFace* currFace = m_faces.getpEntity();

        VDEdge*   currVEdge = currFace->getVEdge();


        if ( currVEdge == rg_NULL ) {
            fout << "FACE\t" << currFace->getID() << "\t" << "no V-edge" << endl;   
            continue;
        }
    }


    m_edges.reset4Loop();
    while ( m_edges.setNext4Loop() ) {
        BetaEdge* currEdge  = m_edges.getpEntity();
        VDFace*   currVFace = currEdge->getVFace();

        if ( currVFace == rg_NULL ) {
            fout << "EDGE\t" << currEdge->getID() << "\t" << "no V-face" << endl;            
        }
    }


    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() ) {
        BetaVertex* currVtx = m_vertices.getpEntity();
        VDCell*     currVCell = currVtx->getVCell();
        if ( currVCell == rg_NULL ) {
            fout << "VTX \t" << currVtx->getID() << "\t" << "no V-cell" << endl;            
        }
        else {
            Ball*          currBallOnBVtx  = currVtx->getBallProperty();
            BallGenerator* currBallOnVCell = currVCell->getGenerator();

            if ( currBallOnBVtx!=rg_NULL && currBallOnVCell!=rg_NULL ) {
                if ( !(currBallOnBVtx->getGeometry() == currBallOnVCell->getBall()) ) {
                    fout << "VTX \t" << currVtx->getID() << "\t" << "diff. balls" << endl;            
                }
            }
            else if ( (currBallOnBVtx==rg_NULL && currBallOnVCell!=rg_NULL) || (currBallOnBVtx!=rg_NULL && currBallOnVCell==rg_NULL)) {
                fout << "VTX \t" << currVtx->getID() << "\t" << "diff. balls" << endl;            
            }
            else {

            }
        }
    }
}

/*
///////////////////////////////////////////////////////////////////////////
//
//  Quasi-operator
//
//  1. Query primitives for intra-world.
//
//  1.1. Q(v, Y | W): Query primitives for intra-world.
void BetaUniverse::inquireVerticesOfVertexInIntraWorld(BetaVertex* givenVertex, rg_dList<BetaVertex*>& vertexList, BetaFace* world )
{
    givenVertex->inquireVerticesOfVertexInIntraWorld( vertexList, world );
}



void BetaUniverse::inquireEdgesOfVertexInIntraWorld(   BetaVertex* givenVertex, rg_dList<BetaEdge*>&   edgeList,   BetaFace* world )
{
    givenVertex->inquireEdgesOfVertexInIntraWorld( edgeList, world );
}



void BetaUniverse::inquireFacesOfVertexInIntraWorld(   BetaVertex* givenVertex, rg_dList<BetaFace*>&   faceList,   BetaFace* world )
{
    givenVertex->inquireFacesOfVertexInIntraWorld( faceList, world );
}



void BetaUniverse::inquireCellsOfVertexInIntraWorld(   BetaVertex* givenVertex, rg_dList<BetaCell*>&   cellList,   BetaFace* world )
{
    givenVertex->inquireCellsOfVertexInIntraWorld( cellList, world );
}




//  1.2. Q(e, Y | W): Query primitives for intra-world.
void BetaUniverse::inquireVerticesOfEdgeInIntraWorld( BetaEdge* givenEdge, rg_dList<BetaVertex*>& vertexList, BetaFace* world )
{
    givenEdge->inquireVerticesOfEdgeInIntraWorld( vertexList, world );
}



void BetaUniverse::inquireEdgesOfEdgeInIntraWorld(    BetaEdge* givenEdge, rg_dList<BetaEdge*>&   edgeList,   BetaFace* world )
{
    givenEdge->inquireEdgesOfEdgeInIntraWorld( edgeList, world );
}



void BetaUniverse::inquireFacesOfEdgeInIntraWorld(    BetaEdge* givenEdge, rg_dList<BetaFace*>&   faceList,   BetaFace* world )
{
    givenEdge->inquireFacesOfEdgeInIntraWorld( faceList, world );
}



void BetaUniverse::inquireCellsOfEdgeInIntraWorld(    BetaEdge* givenEdge, rg_dList<BetaCell*>&   cellList,   BetaFace* world )
{
    givenEdge->inquireCellsOfEdgeInIntraWorld( cellList, world );
}




//  1.3. Q(f, Y | W): Query primitives for intra-world.
void BetaUniverse::inquireVerticesOfFaceInIntraWorld( BetaFace* givenFace, rg_dList<BetaVertex*>& vertexList, BetaFace* world )
{
    givenFace->inquireVerticesOfFaceInIntraWorld( vertexList, world );
}



void BetaUniverse::inquireEdgesOfFaceInIntraWorld(    BetaFace* givenFace, rg_dList<BetaEdge*>&   edgeList,   BetaFace* world )
{
    givenFace->inquireEdgesOfFaceInIntraWorld( edgeList, world );
}



void BetaUniverse::inquireFacesOfFaceInIntraWorld(    BetaFace* givenFace, rg_dList<BetaFace*>&   faceList,   BetaFace* world )
{
    givenFace->inquireFacesOfFaceInIntraWorld( faceList, world );
}



void BetaUniverse::inquireCellsOfFaceInIntraWorld(    BetaFace* givenFace, rg_dList<BetaCell*>&   cellList,   BetaFace* world )
{
    givenFace->inquireCellsOfFaceInIntraWorld( cellList, world );
}




//  1.4. Q(c, Y | W): Query primitives for intra-world.
void BetaUniverse::inquireVerticesOfCellInIntraWorld(BetaCell* givenCell, rg_dList<BetaVertex*>& vertexList, BetaFace* world )
{
    givenCell->inquireVerticesOfCellInIntraWorld( vertexList, world );
}



void BetaUniverse::inquireEdgesOfCellInIntraWorld(   BetaCell* givenCell, rg_dList<BetaEdge*>&   edgeList,   BetaFace* world )
{
    givenCell->inquireEdgesOfCellInIntraWorld( edgeList, world );
}



void BetaUniverse::inquireFacesOfCellInIntraWorld(   BetaCell* givenCell, rg_dList<BetaFace*>&   faceList,   BetaFace* world )
{
    givenCell->inquireFacesOfCellInIntraWorld( faceList, world );
}



void BetaUniverse::inquireCellsOfCellInIntraWorld(   BetaCell* givenCell, rg_dList<BetaCell*>&   cellList,   BetaFace* world )
{
    givenCell->inquireCellsOfCellInIntraWorld( cellList, world );
}




//  2. Query primitives for inter-world.
//
void BetaUniverse::inquireSmallWorldsAtVertex(BetaVertex* givenVertex, rg_dList<BetaFace*>& smallWorldList, BetaFace* world )
{
    givenVertex->inquireSmallWorldsAtVertex( smallWorldList, world );
}



void BetaUniverse::inquireSmallWorldsAtEdge(  BetaEdge*   givenEdge,   rg_dList<BetaFace*>& smallWorldList, BetaFace* world )
{
    givenEdge->inquireSmallWorldsAtEdge( smallWorldList, world );
}



void BetaUniverse::inquireSmallWorldsAtFace(  BetaFace*   givenFace,   rg_dList<BetaFace*>& smallWorldList, BetaFace* world )
{
    givenFace->inquireSmallWorldsAtFace( smallWorldList, world );
}



void BetaUniverse::inquireSmallWorldsAtCell(  BetaCell*   givenCell,   rg_dList<BetaFace*>& smallWorldList, BetaFace* world )
{
    givenCell->inquireSmallWorldsAtCell( smallWorldList, world );
}




void BetaUniverse::inquireWholeSmallWorldsAtVertex(BetaVertex* givenVertex, rg_dList<BetaFace*>& smallWorldList, BetaFace* world )
{
    givenVertex->inquireWholeSmallWorldsAtVertex( smallWorldList, world );
}




//  3. Query primitives for whole-world.
//
//  3.1. Q(v, Y): Query primitives for whole-world.
void BetaUniverse::inquireVerticesOfVertex(BetaVertex* givenVertex, rg_dList<BetaVertex*>& vertexList )
{
    givenVertex->inquireVerticesOfVertex( vertexList );
}



void BetaUniverse::inquireEdgesOfVertex(   BetaVertex* givenVertex, rg_dList<BetaEdge*>&   edgeList   )
{
    givenVertex->inquireEdgesOfVertex( edgeList );
}



void BetaUniverse::inquireFacesOfVertex(   BetaVertex* givenVertex, rg_dList<BetaFace*>&   faceList   )
{
    givenVertex->inquireFacesOfVertex( faceList );
}



void BetaUniverse::inquireCellsOfVertex(   BetaVertex* givenVertex, rg_dList<BetaCell*>&   cellList   )
{
    givenVertex->inquireCellsOfVertex( cellList );
}




//  3.2. Q(e, Y): Query primitives for whole-world.
void BetaUniverse::inquireVerticesOfEdge( BetaEdge* givenEdge, rg_dList<BetaVertex*>& vertexList )
{
    givenEdge->inquireVerticesOfEdge( vertexList );
}



void BetaUniverse::inquireEdgesOfEdge(    BetaEdge* givenEdge, rg_dList<BetaEdge*>&   edgeList   )
{
    givenEdge->inquireEdgesOfEdge( edgeList );
}



void BetaUniverse::inquireFacesOfEdge(    BetaEdge* givenEdge, rg_dList<BetaFace*>&   faceList   )
{
    givenEdge->inquireFacesOfEdge( faceList );
}



void BetaUniverse::inquireCellsOfEdge(    BetaEdge* givenEdge, rg_dList<BetaCell*>&   cellList   )
{
    givenEdge->inquireCellsOfEdge( cellList );
}




//  3.3. Q(f, Y): Query primitives for whole-world.
void BetaUniverse::inquireVerticesOfFace( BetaFace* givenFace, rg_dList<BetaVertex*>& vertexList )
{
    givenFace->inquireVerticesOfFace( vertexList );
}



void BetaUniverse::inquireEdgesOfFace(    BetaFace* givenFace, rg_dList<BetaEdge*>&   edgeList   )
{
    givenFace->inquireEdgesOfFace( edgeList );
}



void BetaUniverse::inquireFacesOfFace(    BetaFace* givenFace, rg_dList<BetaFace*>&   faceList   )
{
    givenFace->inquireFacesOfFace( faceList );
}



void BetaUniverse::inquireCellsOfFace(    BetaFace* givenFace, rg_dList<BetaCell*>&   cellList   )
{
    givenFace->inquireCellsOfFace( cellList );
}




//  3.4. Q(c, Y): Query primitives for whole-world.
void BetaUniverse::inquireVerticesOfCell(BetaCell* givenCell, rg_dList<BetaVertex*>& vertexList )
{
    givenCell->inquireVerticesOfCell( vertexList );
}



void BetaUniverse::inquireEdgesOfCell(   BetaCell* givenCell, rg_dList<BetaEdge*>&   edgeList   )
{
    givenCell->inquireEdgesOfCell( edgeList );
}



void BetaUniverse::inquireFacesOfCell(   BetaCell* givenCell, rg_dList<BetaFace*>&   faceList   )
{
    givenCell->inquireFacesOfCell( faceList );
}



void BetaUniverse::inquireCellsOfCell(   BetaCell* givenCell, rg_dList<BetaCell*>&   cellList   )
{
    givenCell->inquireCellsOfCell( cellList );
}




//
//  END of Quasi-operator
//  
///////////////////////////////////////////////////////////////////////////
*/



// void BetaUniverse::reportQT(ofstream& fout)
// {
//     BetaCell* currCell = rg_NULL;
//     m_cells.reset4Loop();
//     while ( m_cells.setNext4Loop() ) {
//         currCell = m_cells.getpEntity();
// 
//         fout << "Cell " << currCell->getID() << " : ";
// 
//         fout << "vertices ";
// 	    BetaVertex** vertex = currCell->getVertices();
//         for( int i=0; i<4; i++ ) {
//             if ( vertex[i] == rg_NULL )
//                 fout << "null" << "\t";
//             else
//                 fout << vertex[i]->getID() << "\t";
//         }
// 
//         fout << "faces ";	
//         BetaFace** face = currCell->getFaces();
//         int j = 0;
//         for( j=0; j<4; j++ ) {
//             if ( face[j] == rg_NULL )
//                 fout << "null" << "\t";
//             else
//                 fout << face[j]->getID() << "\t";
//         }
// 
//         fout << "neighbor ";	
//         for( j=0; j<4; j++ ) {
//             if ( face[j] == rg_NULL )
//                 fout << "null" << "\t";
//             else {
//                 if ( face[j]->getLeftCell() == currCell )
//                     fout << face[j]->getRightCell()->getID() << "\t";
//                 else
//                     fout << face[j]->getLeftCell()->getID() << "\t";
//             }
//         }
//         fout << endl;
//     }
//     fout << endl;
// 
// 
//     BetaFace* currFace = rg_NULL;
//     m_faces.reset4Loop();
//     while ( m_faces.setNext4Loop() ) {
//         currFace = m_faces.getpEntity();
// 
//         fout << "Face " << currFace->getID() << " : ";
// 
//         if ( currFace->getLeftCell() == rg_NULL ) {
//             fout << "lc null \t";
//             fout << "rc null \t";
//         }
//         else {
//             fout << "lc " << currFace->getLeftCell()->getID() << "\t";
//             fout << "rc " << currFace->getRightCell()->getID() << "\t";
//         }
// 
//         fout << "vtx ";
//         BetaCell* rightCell = currFace->getRightCell();
//         if ( rightCell == rg_NULL ) {
//             fout << "null" << "\t" << "null" << "\t" << "null" << "\t";
//         }
//         else {
// 
//             rg_INT posFace = -1;
//             for ( int i=0; i<4; i++ ) {
//                 if ( currFace == rightCell->getFace(i) ) {
//                     posFace = i;
//                     break;
//                 }
//             }
// 
//             if ( posFace == 0 ) {
//                 fout << rightCell->getVertex(1)->getID() << "\t" 
//                      << rightCell->getVertex(2)->getID() << "\t" 
//                      << rightCell->getVertex(3)->getID() << "\t";
//             }
//             else if ( posFace == 1 ) {
//                 fout << rightCell->getVertex(0)->getID() << "\t" 
//                      << rightCell->getVertex(3)->getID() << "\t" 
//                      << rightCell->getVertex(2)->getID() << "\t";
//             }
//             else if ( posFace == 2 ) {
//                 fout << rightCell->getVertex(0)->getID() << "\t" 
//                      << rightCell->getVertex(1)->getID() << "\t" 
//                      << rightCell->getVertex(3)->getID() << "\t";
//             }
//             else { //if ( posFace == 3 ) {
//                 fout << rightCell->getVertex(0)->getID() << "\t" 
//                      << rightCell->getVertex(2)->getID() << "\t" 
//                      << rightCell->getVertex(1)->getID() << "\t";
//             }
// 
//         }
// 
//         fout << "edges ";
// 	    BetaEdge** edge = currFace->getEdges();
//         rg_BOOL*   ori  = currFace->getEdgeOrientations();
//         for ( int i=0; i<3; i++ ) {
//             if ( edge[i] == rg_NULL )
//                 fout << "null" << "\t";
//             else {
//                 fout << edge[i]->getID();
//                 if ( ori[i] ) {
//                     fout << "(+)\t";
//                 }
//                 else {
//                     fout << "(-)\t";
//                 }
//             }
//         }
//         fout << endl;
//     }
//     fout << endl;
// 
// 
//     BetaEdge* currEdge = rg_NULL;
//     m_edges.reset4Loop();
//     while ( m_edges.setNext4Loop() ) {
//         currEdge = m_edges.getpEntity();
// 
//         fout << "Edge " << currEdge->getID() << " : ";
// 
//         if ( currEdge->getStartVertex() == rg_NULL )
//             fout << "sv " << "null" << "\t";
//         else
//             fout << "sv " << currEdge->getStartVertex()->getID() << "\t";
// 
//         if ( currEdge->getEndVertex() == rg_NULL )
//             fout << "ev " << "null" << "\t";
//         else
//             fout << "ev " << currEdge->getEndVertex()->getID() << "\t";
// 
//         if ( currEdge->getFirstFace() == rg_NULL )
//             fout << "face " << "null" << "\t";
//         else
//             fout << "face " << currEdge->getFirstFace()->getID() << "\t";
// 
// 
//         fout << "small world ";
//         rg_dList<BetaFace*>* smallWorld = currEdge->getSmallWorlds();
//         BetaFace* currWorld = rg_NULL;
//         smallWorld->reset4Loop();
//         while ( smallWorld->setNext4Loop() ) {
//             currWorld = smallWorld->getEntity();
//             fout << currWorld->getID() << "\t";
//         }
//         fout << endl;
//     }
//     fout << endl;
// 
// 
//     BetaVertex* currVtx = rg_NULL;
//     m_vertices.reset4Loop();
//     while ( m_vertices.setNext4Loop() ) {
//         currVtx = m_vertices.getpEntity();
// 
//         fout << "Vertex " << currVtx->getID() << " : ";
// 
//         if ( currVtx->getFirstEdge() == rg_NULL )
//             fout << "edge " << "null" << "\t";
//         else
//             fout << "edge " << currVtx->getFirstEdge()->getID() << "\t";
//         
//         if ( currVtx->getFirstCell() == rg_NULL ) {
//             fout << "cell null \t";
//         }
//         else {
//             fout << "cell " << currVtx->getFirstCell()->getID() << "\t";
//         }
// 
//         Sphere ball = currVtx->getBall();
//         fout << ball.getCenter().getX() << "\t" << ball.getCenter().getY() << "\t"<< ball.getCenter().getZ() << "\t"<< ball.getRadius() << "\t";
//         fout << endl;
//     }
// 
// }
// 
// 
// 
// 
// void BetaUniverse::reportBetaSpan() 
// {
//     ofstream fout;
//     fout.open("test_beta_span.txt");
// 
//     BetaCell* currCell = rg_NULL;
//     m_cells.reset4Loop();
//     while ( m_cells.setNext4Loop() ) {
//         currCell = m_cells.getpEntity();
// 
//         if ( currCell->isVirtual() ) {
//             continue;
//         }
// 
//         fout << "Cell " << currCell->getID() << " : ";
//         BetaSpanCore span = currCell->getBetaSpan();
//         span.report( fout );
//     }
//     fout << endl;
// 
// 
//     BetaFace* currFace = rg_NULL;
//     m_faces.reset4Loop();
//     while ( m_faces.setNext4Loop() ) {
//         currFace = m_faces.getpEntity();
// 
//         if ( currFace->isVirtual() ) {
//             continue;
//         }
// 
//         fout << "Face " << currFace->getID() << " : ";
//         BetaSpanCore span = currFace->getBetaSpan();
//         span.report( fout );
//     }
//     fout << endl;
// 
// 
//     BetaEdge* currEdge = rg_NULL;
//     m_edges.reset4Loop();
//     while ( m_edges.setNext4Loop() ) {
//         currEdge = m_edges.getpEntity();
// 
//         if ( currEdge->isVirtual() ) {
//             continue;
//         }
// 
//         fout << "Edge " << currEdge->getID() << " : ";
//         BetaSpanCore span = currEdge->getBetaSpan();
//         span.report( fout );
//     }
//     fout << endl;
// 
// 
//     BetaVertex* currVtx = rg_NULL;
//     m_vertices.reset4Loop();
//     while ( m_vertices.setNext4Loop() ) {
//         currVtx = m_vertices.getpEntity();
// 
//         if ( currVtx->isVirtual() ) {
//             continue;
//         }
// 
//         fout << "Vertex " << currVtx->getID() << " : ";
//         BetaSpanCore span = currVtx->getBetaSpan();
//         span.report( fout );
//     }
// 
//     fout.close();
// }
// 
// 
// 
// 
// void BetaUniverse::testQuery_v_C(ofstream& fout)
// {
// //     fout << "* Q(v,C) test: ";
// //     fout << "# of vertices:" << m_vertices.getSize() << endl << endl;
// 
//     list<int> simplexInFailure;
// 
//     BetaVertex* currVtx = rg_NULL;
//     m_vertices.reset4Loop();
//     while ( m_vertices.setNext4Loop() ) {
//         currVtx = m_vertices.getpEntity();
// 
//         rg_dList<BetaCell*> cellList;
//         currVtx->searchCellsInWholeWorld(cellList);
// 
//         list<int> cell_localSearch;
//         cellList.reset4Loop();
//         while ( cellList.setNext4Loop() ) {
//             cell_localSearch.push_back( cellList.getEntity()->getID() );
//         }
// 
//         list<int> cell_globalSearch;
//         BetaCell* currCell = rg_NULL;
//         m_cells.reset4Loop();
//         while ( m_cells.setNext4Loop() ) {
//             currCell = m_cells.getpEntity();
// 
//             if ( currCell->isThere( currVtx ) ) {
//                 cell_globalSearch.push_back( currCell->getID() );
//             }
//         }
// 
// 
//         cell_localSearch.sort();
//         cell_globalSearch.sort();
// 
//         if ( !(cell_localSearch == cell_globalSearch) ) {
//             simplexInFailure.push_back( currVtx->getID() );
//         }
// 
// //         fout << "Vertex " << currVtx->getID() << ": Q(v,C)";
// // 
// //         if ( cell_localSearch == cell_globalSearch ) {
// //             fout << " SUCCESS!" << endl;
// //         }
// //         else {
// //             fout << " FAIL!" << endl;
// //             fout << "\tlocal search:\t";
// //             list< int >::iterator	it;
// // 	        for( it = cell_localSearch.begin(); it != cell_localSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //             fout << "\tglobal search:\t";
// // 	        for( it = cell_globalSearch.begin(); it != cell_globalSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //         }
// //         fout << endl;
//     }
// 
//     if ( simplexInFailure.size() == 0 ) {
//         fout << "success\t";
//     }
//     else {
//         fout << "fail\t";
// //         list< int >::iterator	it;
// // 	    for( it = simplexInFailure.begin(); it != simplexInFailure.end(); it++)  {
// //             fout << (*it) << "\t";
// //         }
// //         fout << endl;
//     }
// }
// 
// 
// 
// void BetaUniverse::testQuery_v_F(ofstream& fout)
// {
// //     fout << "* Q(v,F) test: ";
// //     fout << "****** Q(v,F) test ******" << endl;
// //     fout << "# of vertices:" << m_vertices.getSize() << endl << endl;
// 
//     list<int> simplexInFailure;
// 
//     BetaVertex* currVtx = rg_NULL;
//     m_vertices.reset4Loop();
//     while ( m_vertices.setNext4Loop() ) {
//         currVtx = m_vertices.getpEntity();
// 
//         rg_dList<BetaFace*> faceList;
//         currVtx->searchFacesInWholeWorld(faceList);
// 
//         list<int> face_localSearch;
//         faceList.reset4Loop();
//         while ( faceList.setNext4Loop() ) {
//             face_localSearch.push_back( faceList.getEntity()->getID() );
//         }
// 
//         list<int> face_globalSearch;
//         BetaFace* currFace = rg_NULL;
//         m_faces.reset4Loop();
//         while ( m_faces.setNext4Loop() ) {
//             currFace = m_faces.getpEntity();
// 
//             if ( currFace->isThere( currVtx ) ) {
//                 face_globalSearch.push_back( currFace->getID() );
//             }
//         }
// 
// 
//         face_localSearch.sort();
//         face_globalSearch.sort();
// 
//         if ( !(face_localSearch == face_globalSearch) ) {
//             simplexInFailure.push_back( currVtx->getID() );
//         }
// 
// //         fout << "Vertex " << currVtx->getID() << ": Q(v,F)";
// // 
// //         if ( face_localSearch == face_globalSearch ) {
// //             fout << " SUCCESS!" << endl;
// //         }
// //         else {
// //             fout << " FAIL!" << endl;
// //             fout << "\tlocal search:\t";
// //             list< int >::iterator	it;
// // 	           for( it = face_localSearch.begin(); it != face_localSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //             fout << "\tglobal search:\t";
// // 	           for( it = face_globalSearch.begin(); it != face_globalSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //         }
// //         fout << endl;
//     }
// 
// 
//     if ( simplexInFailure.size() == 0 ) {
//         fout << "success\t";
//     }
//     else {
//         fout << "fail\t";
// //         list< int >::iterator	it;
// // 	    for( it = simplexInFailure.begin(); it != simplexInFailure.end(); it++)  {
// //             fout << (*it) << "\t";
// //         }
// //         fout << endl;
//     }
// }
// 
// 
// 
// void BetaUniverse::testQuery_v_E(ofstream& fout)
// {
// 
// //     fout << "****** Q(v,E) test ******" << endl << endl;
// 
// //     fout << "* Q(v,E) test: ";
//     list<int> simplexInFailure;
// 
// 
//     BetaVertex* currVtx = rg_NULL;
//     m_vertices.reset4Loop();
//     while ( m_vertices.setNext4Loop() ) {
//         currVtx = m_vertices.getpEntity();
// 
//         rg_dList<BetaEdge*> edgeList;
//         currVtx->searchEdgesInWholeWorld(edgeList);
// 
//         list<int> edge_localSearch;
//         edgeList.reset4Loop();
//         while ( edgeList.setNext4Loop() ) {
//             edge_localSearch.push_back( edgeList.getEntity()->getID() );
//         }
// 
//         list<int> edge_globalSearch;
//         BetaEdge* currEdge = rg_NULL;
//         m_edges.reset4Loop();
//         while ( m_edges.setNext4Loop() ) {
//             currEdge = m_edges.getpEntity();
// 
//             if ( currEdge->isThere( currVtx ) ) {
//                 edge_globalSearch.push_back( currEdge->getID() );
//             }
//         }
// 
// 
//         edge_localSearch.sort();
//         edge_globalSearch.sort();
// 
//         if ( !(edge_localSearch == edge_globalSearch) ) {
//             simplexInFailure.push_back( currVtx->getID() );
//         }
// 
// //         fout << "Vertex " << currVtx->getID() << ": Q(v,E)";
// // 
// //         if ( edge_localSearch == edge_globalSearch ) {
// //             fout << " SUCCESS!" << endl;
// //         }
// //         else {
// //             fout << " FAIL!" << endl;
// //             fout << "\tlocal search:\t";
// //             list< int >::iterator	it;
// // 	        for( it = edge_localSearch.begin(); it != edge_localSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //             fout << "\tglobal search:\t";
// // 	        for( it = edge_globalSearch.begin(); it != edge_globalSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //         }
// //         fout << endl;
//     }
// 
// 
//     if ( simplexInFailure.size() == 0 ) {
//         fout << "success\t";
//     }
//     else {
//         fout << "fail\t";
// //         list< int >::iterator	it;
// // 	    for( it = simplexInFailure.begin(); it != simplexInFailure.end(); it++)  {
// //             fout << (*it) << "\t";
// //         }
// //         fout << endl;
//     }
// }
// 
// 
// 
// void BetaUniverse::testQuery_v_V(ofstream& fout)
// {
// //     fout << "* Q(v,V) test: ";
//     list<int> simplexInFailure;
// 
// 
// //     fout << "****** Q(v,V) test ******" << endl << endl;
//     BetaVertex* currVtx = rg_NULL;
//     m_vertices.reset4Loop();
//     while ( m_vertices.setNext4Loop() ) {
//         currVtx = m_vertices.getpEntity();
// 
//         rg_dList<BetaVertex*> vertexList;
//         currVtx->searchVerticesInWholeWorld(vertexList);
// 
//         list<int> vertex_localSearch;
//         vertexList.reset4Loop();
//         while ( vertexList.setNext4Loop() ) {
//             vertex_localSearch.push_back( vertexList.getEntity()->getID() );
//         }
// 
//         list<int> vertex_globalSearch;
//         BetaEdge* currEdge = rg_NULL;
//         m_edges.reset4Loop();
//         while ( m_edges.setNext4Loop() ) {
//             currEdge = m_edges.getpEntity();
// 
//             if ( currEdge->isThere( currVtx ) ) {
//                 if ( currVtx == currEdge->getStartVertex() ) {
//                     vertex_globalSearch.push_back( currEdge->getEndVertex()->getID() );
//                 }
//                 else {
//                     vertex_globalSearch.push_back( currEdge->getStartVertex()->getID() );
//                 }
//             }
//         }
// 
// 
//         vertex_localSearch.sort();
//         vertex_globalSearch.sort();
// 
// 
//         if ( !(vertex_localSearch == vertex_globalSearch) ) {
//             simplexInFailure.push_back( currVtx->getID() );
//         }
//         
// //         fout << "Vertex " << currVtx->getID() << ": Q(v,V)";
// // 
// //         if ( vertex_localSearch == vertex_globalSearch ) {
// //             fout << " SUCCESS!" << endl;
// //         }
// //         else {
// //             fout << " FAIL!" << endl;
// // 
// //             fout << "\tlocal search:\t";
// //             list< int >::iterator	it;
// // 	        for( it = vertex_localSearch.begin(); it != vertex_localSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //             fout << "\tglobal search:\t";
// // 	        for( it = vertex_globalSearch.begin(); it != vertex_globalSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //         }
// //         fout << endl;
//     }
// 
// 
//     if ( simplexInFailure.size() == 0 ) {
//         fout << "success\t";
//     }
//     else {
//         fout << "fail\t";
// //         list< int >::iterator	it;
// // 	    for( it = simplexInFailure.begin(); it != simplexInFailure.end(); it++)  {
// //             fout << (*it) << "\t";
// //         }
// //         fout << endl;
//     }
// }
// 
// 
// 
// 
// void BetaUniverse::testQuery_e_C(ofstream& fout)
// {
// //     fout << "* Q(e,C) test: ";
//     list<int> simplexInFailure;
// 
// 
// //     fout << "****** Q(e,C) test ******" << endl;
// //     fout << "# of edges:" << m_edges.getSize() << endl << endl;
// 
//     BetaEdge* currEdge = rg_NULL;
//     m_edges.reset4Loop();
//     while ( m_edges.setNext4Loop() ) {
//         currEdge = m_edges.getpEntity();
// 
//         rg_dList<BetaCell*> cellList;
//         currEdge->searchCellsInWholeWorld(cellList);
//         list<int> cell_localSearch;
//         cellList.reset4Loop();
//         while ( cellList.setNext4Loop() ) {
//             cell_localSearch.push_back( cellList.getEntity()->getID() );
//         }
// 
//         list<int> cell_wholeSearch;
//         BetaCell* currCell = rg_NULL;
//         m_cells.reset4Loop();
//         while ( m_cells.setNext4Loop() ) {
//             currCell = m_cells.getpEntity();
// 
//             if ( currCell->isThere( currEdge ) ) {
//                 cell_wholeSearch.push_back( currCell->getID() );
//             }
//         }
// 
// 
//         cell_localSearch.sort();
//         cell_wholeSearch.sort();
// 
//         if ( !(cell_localSearch == cell_wholeSearch) ) {
//             simplexInFailure.push_back( currEdge->getID() );
//         }
// 
// //         fout << "Edge " << currEdge->getID() << ": Q(e,C)";
// // 
// //         if ( cell_localSearch == cell_wholeSearch ) {
// //             fout << " SUCCESS!" << endl;
// //         }
// //         else {
// //             fout << " FAIL!" << endl;
// //             fout << "\tlocal search:\t";
// //             list< int >::iterator	it;
// // 	        for( it = cell_localSearch.begin(); it != cell_localSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //             fout << "\tglobal search:\t";
// // 	        for( it = cell_wholeSearch.begin(); it != cell_wholeSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //         }
// //         fout << endl;
// 
//     }
// 
// 
//     if ( simplexInFailure.size() == 0 ) {
//         fout << "success\t";
//     }
//     else {
//         fout << "fail\t";
// //         list< int >::iterator	it;
// // 	    for( it = simplexInFailure.begin(); it != simplexInFailure.end(); it++)  {
// //             fout << (*it) << "\t";
// //         }
// //         fout << endl;
//     }
// }
// 
// 
// 
// void BetaUniverse::testQuery_e_F(ofstream& fout)
// {
// //     fout << "* Q(e,F) test: ";
//     list<int> simplexInFailure;
// 
// //     fout << "****** Q(e,F) test ******" << endl;
// //     fout << "# of edges:" << m_edges.getSize() << endl << endl;
// 
//     BetaEdge* currEdge = rg_NULL;
//     m_edges.reset4Loop();
//     while ( m_edges.setNext4Loop() ) {
//         currEdge = m_edges.getpEntity();
// 
//         rg_dList<BetaFace*> faceList;
//         currEdge->searchFacesInWholeWorld(faceList);
//         list<int> face_localSearch;
//         faceList.reset4Loop();
//         while ( faceList.setNext4Loop() ) {
//             face_localSearch.push_back( faceList.getEntity()->getID() );
//         }
// 
//         list<int> face_wholeSearch;
//         BetaFace* currFace = rg_NULL;
//         m_faces.reset4Loop();
//         while ( m_faces.setNext4Loop() ) {
//             currFace = m_faces.getpEntity();
// 
//             if ( currFace->isThere( currEdge ) ) {
//                 face_wholeSearch.push_back( currFace->getID() );
//             }
//         }
// 
// 
//         face_localSearch.sort();
//         face_wholeSearch.sort();
// 
//         if ( !(face_localSearch == face_wholeSearch) ) {
//             simplexInFailure.push_back( currEdge->getID() );
//         }
// 
// //         fout << "Edge " << currEdge->getID() << ": Q(e,F)";
// // 
// //         if ( face_localSearch == face_wholeSearch ) {
// //             fout << " SUCCESS!" << endl;
// //         }
// //         else {
// //             fout << " FAIL!" << endl;
// //             fout << "\tlocal search:\t";
// //             list< int >::iterator	it;
// // 	        for( it = face_localSearch.begin(); it != face_localSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //             fout << "\tglobal search:\t";
// // 	        for( it = face_wholeSearch.begin(); it != face_wholeSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //         }
// //         fout << endl;
// 
//     }
// 
// 
//     if ( simplexInFailure.size() == 0 ) {
//         fout << "success\t";
//     }
//     else {
//         fout << "fail\t";
// //         list< int >::iterator	it;
// // 	    for( it = simplexInFailure.begin(); it != simplexInFailure.end(); it++)  {
// //             fout << (*it) << "\t";
// //         }
// //         fout << endl;
//     }
// }
// 
// 
// 
// void BetaUniverse::testQuery_e_E(ofstream& fout)
// {
// //     fout << "* Q(e,E) test: ";
//     list<int> simplexInFailure;
// 
// //     fout << "****** Q(e,E) test ******" << endl;
// //     fout << "# of edges:" << m_edges.getSize() << endl << endl;
// 
//     BetaEdge* currEdge = rg_NULL;
//     m_edges.reset4Loop();
//     while ( m_edges.setNext4Loop() ) {
//         currEdge = m_edges.getpEntity();
// 
//         rg_dList<BetaEdge*> edgeList;
//         currEdge->searchEdgesInWholeWorld(edgeList);
//         list<int> edge_localSearch;
//         edgeList.reset4Loop();
//         while ( edgeList.setNext4Loop() ) {
//             edge_localSearch.push_back( edgeList.getEntity()->getID() );
//         }
// 
//         list<int> edge_wholeSearch;
//         BetaFace* currFace = rg_NULL;
//         m_faces.reset4Loop();
//         while ( m_faces.setNext4Loop() ) {
//             currFace = m_faces.getpEntity();
// 
//             if ( currFace->isThere( currEdge ) ) {
//                 BetaEdge** edge = currFace->getEdges();
//                 for ( rg_INT i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++ ) {
//                     if ( currEdge != edge[i] ) {
//                         edge_wholeSearch.push_back( edge[i]->getID() );
//                     }
//                 }
//             }
//         }
// 
// 
//         edge_localSearch.sort();
//         edge_wholeSearch.sort();
// 
// 
//         if ( !(edge_localSearch == edge_wholeSearch) ) {
//             simplexInFailure.push_back( currEdge->getID() );
//         }
// 
// 
// //         fout << "Edge " << currEdge->getID() << ": Q(e,E)";
// // 
// //         if ( edge_localSearch == edge_wholeSearch ) {
// //             fout << " SUCCESS!" << endl;
// //         }
// //         else {
// //             fout << " FAIL!" << endl;
// //             fout << "\tlocal search:\t";
// //             list< int >::iterator	it;
// // 	        for( it = edge_localSearch.begin(); it != edge_localSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //             fout << "\tglobal search:\t";
// // 	        for( it = edge_wholeSearch.begin(); it != edge_wholeSearch.end(); it++)  {
// //                 fout << (*it) << "\t";
// //             }
// //             fout << endl;
// //         }
// //         fout << endl;
//     }
// 
// 
//     if ( simplexInFailure.size() == 0 ) {
//         fout << "success\t";
//     }
//     else {
//         fout << "fail\t";
// //         list< int >::iterator	it;
// // 	    for( it = simplexInFailure.begin(); it != simplexInFailure.end(); it++)  {
// //             fout << (*it) << "\t";
// //         }
// //         fout << endl;
//     }
// }
