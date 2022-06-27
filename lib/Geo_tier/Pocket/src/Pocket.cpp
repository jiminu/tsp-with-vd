#include "Pocket.h"
//#include "BetaUniverse.h"
#include "ManifoldizedBetaShape.h"
//#include "ConstForBetaComplex.h"
#include "FunctionsForPocket.h"

#include "rg_AugmentedBetaComplex.h"

//debug start
//#include "Particle.h"
#include "rg_Atom.h"
#include "Residue.h"
#include "StringFunctions.h"
//debug end

#include <fstream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////
//
// Constructor
Pocket::Pocket()
{
    m_innerBetaVal             = 0.0;
    m_outerBetaVal             = 0.0;

    m_pocketPrimitives      = rg_NULL;    
    m_sortedPocketPrimitives= rg_NULL;

	m_numOfPocketPrimitives = 0;
}



Pocket::~Pocket()
{
	removeAll();    
}




void Pocket::clear()
{
    removeAll();
}



/////////////////////////////////////////////////////////////////////////////////
//
// Get Functions

BetaUniverse* Pocket::getQusiTriangulationForBeta()
{
	return m_qusiTriangulationForBeta;
}


rg_REAL Pocket::getInnerBetaVal() const
{
    return m_innerBetaVal;    
}



rg_REAL Pocket::getOuterBetaVal() const
{
	return m_outerBetaVal;
}


ManifoldizedBetaShape* Pocket::getpInnerBetaShape()
{
    return &m_innerBetaShape;
}




PocketPrimitive* Pocket::getPocketPrimitives()
{
	return m_pocketPrimitives;
}



PocketPrimitive** Pocket::getSortedPocketPrimitives()
{
	return m_sortedPocketPrimitives;
}



rg_INT Pocket::getNumOfPocketPrimitives() const
{
	return m_numOfPocketPrimitives;
}



/////////////////////////////////////////////////////////////////////////////////
//
// Set Functions

void Pocket::setQuasiTriangulationForBeta(BetaUniverse* tempBetaComplex)
{
	m_qusiTriangulationForBeta = tempBetaComplex;
}




void Pocket::setOuterBetaVal(const rg_REAL& tempOuterBetaVal)
{
	m_outerBetaVal = tempOuterBetaVal;
}

void Pocket::setInnerBetaVal( const rg_REAL& innerBetaVal )
{
    m_innerBetaVal = innerBetaVal;
}


/////////////////////////////////////////////////////////////////////////////////
//
// Pocket Extraction

void Pocket::constructPocket( BetaUniverse* betaUniverse, const rg_REAL& outerBetaVal, const rg_REAL& innerBetaVal)
{
	m_qusiTriangulationForBeta	= betaUniverse;

    m_outerBetaVal		        = outerBetaVal;
    m_innerBetaVal              = innerBetaVal;
    if ( m_outerBetaVal < innerBetaVal ) {
        m_outerBetaVal = innerBetaVal;
    }
    m_innerBetaShape.setBetaValue( innerBetaVal );
	
	constructPocket();
}



void Pocket::constructPocket()
{    
	removeAll();

    if (    m_innerBetaShape.getBodies()->getSize() == 0
         || m_innerBetaShape.getBetaValue() != m_innerBetaVal ) {
        AugmentedBetaComplex augBC;
	    augBC.construct( m_innerBetaVal, m_qusiTriangulationForBeta, CONSTRUCT_WITH_REGULAR_BETA_SHAPE );
        augBC.extractBoundaryOfAugmentedBetaComplex( m_innerBetaShape );
    }
	
    
	markVerticesOnPocket();

	markFacesOnPocket();

	generatePocketPrimitiveByClusteringFaces();

    markVerticesToUnvisited();
    
	
    
    //shave
    shavePocketPrimitives();

    //markFaceOnShavedPocket
    markFacesOnCurrentPockets();

    generatePocketPrimitiveByClusteringFaces();
		
	
	
    
    //sortPocketPrimitives( COMPARE_WRT_AREA_OF_FACES );
    sortPocketPrimitives( COMPARE_WRT_NUM_OF_FACES );
}


void Pocket::markFacesOnCurrentPockets()
{
    for ( int i=0; i < m_numOfPocketPrimitives; ++i)	{							//for comparison to previous "raw pocket"
		rg_dList< MBSFace* >* facesOnThisPocketPrimitive
		                	     = m_pocketPrimitives[i].getMBSFaces();

        facesOnThisPocketPrimitive->reset4Loop();
        while ( facesOnThisPocketPrimitive->setNext4Loop() ) {
            MBSFace* currFaceOnCurrPocket = facesOnThisPocketPrimitive->getEntity();

            currFaceOnCurrPocket->isVisited( rg_TRUE );
        }
    }
}


void Pocket::markVerticesOnPocket()
{
    rg_dList< MBSBody* >* bodies = m_innerBetaShape.getBodies();
    bodies->reset4Loop();
    while ( bodies->setNext4Loop() ) {
        MBSBody* currBody = bodies->getEntity();

        MBSShell* exteriorShell = currBody->getExteriorShell();
        

        markVerticesOfVertexList( exteriorShell->getVertices() );

        rg_dList< MBSShell* >* interiorShells = currBody->getInteriorShells();
        interiorShells->reset4Loop();
        while ( interiorShells->setNext4Loop() ) {
            MBSShell* currInteriorShell = interiorShells->getEntity();

            markVerticesOfVertexList( currInteriorShell->getVertices() );
        }
    }

}

void Pocket::markVerticesOfVertexList( rg_dList< MBSVertex* >* vertices )
{    
	vertices->reset4Loop();
	while ( vertices->setNext4Loop() )		{
		MBSVertex* currVertex = vertices->getEntity();

        BetaVertex* currBetaVertex = currVertex->getOriginalBetaVertex();

		rg_REAL upperBoundOfRegularInterval = findUpperBoundOfRegularInterval( currBetaVertex );

		if ( upperBoundOfRegularInterval < m_outerBetaVal)
		{
			currVertex->isVisited( rg_TRUE );
		}
	}

}



rg_REAL Pocket::findUpperBoundOfRegularInterval( BetaVertex* givenVertex )
{
    rg_BetaSpan betaSpan = givenVertex->getBetaSpan();

    rg_REAL upperValueOfBoundingState = REAL_MINUS_INFINITE;

	rg_INT  numOfBetaInterval = betaSpan.getNumBetaInterval();
	rg_INT* boundingStateArray = betaSpan.getBoundingStates();
	rg_REAL* betaIntervalArray = betaSpan.getBetaIntervals();

	for (rg_INT i=0; i < numOfBetaInterval; i++ ) {
		if ( boundingStateArray[i] == REGULAR_SIMPLEX )
		{
			//interval of m_boundingState[i] : [m_betaInterval[i], m_betaInterval[i+1])
			if ( upperValueOfBoundingState < betaIntervalArray[i+1] )
			{
				upperValueOfBoundingState = betaIntervalArray[i+1];
			}
		}
	}

	if ( upperValueOfBoundingState == REAL_MINUS_INFINITE )
		return rg_MAX_REAL;

	return upperValueOfBoundingState;
}



void Pocket::markFacesOnPocket()
{
	rg_dList< MBSBody* >* bodies = m_innerBetaShape.getBodies();
    bodies->reset4Loop();
    while ( bodies->setNext4Loop() ) {
        MBSBody* currBody = bodies->getEntity();

        MBSShell* exteriorShell = currBody->getExteriorShell();        

        markFacesOfFaceList( exteriorShell->getFaces() );

        rg_dList< MBSShell* >* interiorShells = currBody->getInteriorShells();
        interiorShells->reset4Loop();
        while ( interiorShells->setNext4Loop() ) {
            MBSShell* currInteriorShell = interiorShells->getEntity();

            markFacesOfFaceList( currInteriorShell->getFaces() );
        }
    }
}


void Pocket::markFacesOfFaceList( rg_dList< MBSFace* >* faces )
{
	faces->reset4Loop();
	while ( faces->setNext4Loop() )		{
		MBSFace* currFace = faces->getEntity();

        rg_dList< MBSVertex* > boundingVertices;
        currFace->searchBoundingVertices( boundingVertices );

        rg_INT numOfVisitedVertices = 0;
        boundingVertices.reset4Loop();
        while ( boundingVertices.setNext4Loop() ) {
            MBSVertex* currBoundingVertex = boundingVertices.getEntity();

            if ( currBoundingVertex->isVisited() ) {
                numOfVisitedVertices++;
            }
        }

        if ( numOfVisitedVertices >= 3 ) {
            currFace->isVisited( rg_TRUE );
        }
    }
}




void Pocket::generatePocketPrimitiveByClusteringFaces()
{
    rg_dList< PocketPrimitive > pocketPrimitiveList;

    rg_dList< MBSBody* >* bodies = m_innerBetaShape.getBodies();
    bodies->reset4Loop();
    while ( bodies->setNext4Loop() ) {
        MBSBody* currBody = bodies->getEntity();

        MBSShell* exteriorShell = currBody->getExteriorShell();        

        clusterFaces( exteriorShell->getFaces(), pocketPrimitiveList );

        rg_dList< MBSShell* >* interiorShells = currBody->getInteriorShells();
        interiorShells->reset4Loop();
        while ( interiorShells->setNext4Loop() ) {
            MBSShell* currInteriorShell = interiorShells->getEntity();

            clusterFaces( currInteriorShell->getFaces(), pocketPrimitiveList );
        }
    }

    
    m_numOfPocketPrimitives = pocketPrimitiveList.getSize();


    if ( m_pocketPrimitives != rg_NULL ) {
        delete [] m_pocketPrimitives;
    }
	m_pocketPrimitives		= pocketPrimitiveList.getArray();
}


void Pocket::clusterFaces( rg_dList< MBSFace* >* faces, rg_dList< PocketPrimitive >& pocketPrimitiveList )
{

	faces->reset4Loop();
	while ( faces->setNext4Loop() )		{
		MBSFace* currFace = faces->getEntity();

        if ( currFace->isVisited() == rg_FALSE ) {
			continue;
        }

		PocketPrimitive* currPP = pocketPrimitiveList.add( PocketPrimitive() );

		rg_dList< MBSFace* > searchStack;
		searchStack.pushBack( currFace );
		while ( searchStack.getSize() > 0 )	{
			MBSFace* searchingFace = searchStack.popBack();

			currPP->addMBSFaceToThisPrimitive( searchingFace );
			searchingFace->isVisited( rg_FALSE );

			rg_dList< MBSFace* > incidentFaces;
			searchingFace->searchIncidentFaces( incidentFaces );
			incidentFaces.reset4Loop();
			while ( incidentFaces.setNext4Loop() )	{
				MBSFace* currIncidentFace = incidentFaces.getEntity();

                if ( currIncidentFace->isVisited() == rg_FALSE ) {
					continue;
				}

				searchStack.pushBack( currIncidentFace );
			}
		}
	}	
}


void Pocket::markVerticesToUnvisited()
{
    rg_dList< MBSBody* >* bodies = m_innerBetaShape.getBodies();
    bodies->reset4Loop();
    while ( bodies->setNext4Loop() ) {
        MBSBody* currBody = bodies->getEntity();

        MBSShell* exteriorShell = currBody->getExteriorShell();
        

        rg_dList< MBSVertex*>* vertices = exteriorShell->getVertices();
        vertices->reset4Loop();
        while ( vertices->setNext4Loop() ) {
            vertices->getEntity()->isVisited( rg_FALSE );
        }

        rg_dList< MBSShell* >* interiorShells = currBody->getInteriorShells();
        interiorShells->reset4Loop();
        while ( interiorShells->setNext4Loop() ) {
            MBSShell* currInteriorShell = interiorShells->getEntity();

            rg_dList< MBSVertex* >* verticesOnInterior = currInteriorShell->getVertices();
            verticesOnInterior->reset4Loop();
            while ( verticesOnInterior->setNext4Loop() ) {
                verticesOnInterior->getEntity()->isVisited( rg_FALSE );
            }            
        }
    }
}



void Pocket::removeAll()
{
    //m_qusiTriangulationForBeta = rg_NULL;
	//m_innerBetaShape = rg_NULL;
	//m_innerBetaShape.removeAll();
    
    //m_outerBetaVal = DBL_MAX;


    if ( m_pocketPrimitives != rg_NULL )
	{
        delete [] m_pocketPrimitives;
		m_pocketPrimitives = rg_NULL;
	}
	
    if (m_sortedPocketPrimitives != rg_NULL)
    {
        delete [] m_sortedPocketPrimitives;
		m_sortedPocketPrimitives = rg_NULL;
    }

	m_numOfPocketPrimitives = 0;
}



void Pocket::sortPocketPrimitives( const rg_INT& compareWRT )
{
	if(m_sortedPocketPrimitives == rg_NULL)
	{
		m_sortedPocketPrimitives = new PocketPrimitive*[ m_numOfPocketPrimitives ];
		for (rg_INT i=0; i<m_numOfPocketPrimitives; i++)
		{
			m_sortedPocketPrimitives[i] = &( m_pocketPrimitives[i] );
		}
	}
	

    switch( compareWRT )
    {
    case COMPARE_WRT_NUM_OF_FACES :
        qsort( (void *) m_sortedPocketPrimitives, m_numOfPocketPrimitives, sizeof(PocketPrimitive*), comparePocketPrimitivesByNumFaces);
        break;

    case COMPARE_WRT_AREA_OF_FACES :
        qsort( (void *) m_sortedPocketPrimitives, m_numOfPocketPrimitives, sizeof(PocketPrimitive*), comparePocketPrimitivesByAreaOfFaces);
        break;

    default:
        qsort( (void *) m_sortedPocketPrimitives, m_numOfPocketPrimitives, sizeof(PocketPrimitive*), comparePocketPrimitivesByAreaOfFaces);
        break;
    }
}






void Pocket::shavePocketPrimitives()
{
    shrinkPocketPrimitives();
    
    //removeMinorEdgeConnectedFaceGroup();
    
    enlargePocketPrimitives();
}



void Pocket::shrinkPocketPrimitives()
{
    for ( int i=0; i < m_numOfPocketPrimitives; ++i)	{							//for comparison to previous "raw pocket"
	//for ( int i=0; i < 10 && i < m_numOfPocketPrimitives; ++i)	{							//for comparison to previous "raw pocket"
		//PocketPrimitive* currPocketPrimitive = m_sortedPocketPrimitives[i];
        PocketPrimitive* currPocketPrimitive = &m_pocketPrimitives[i];

		shrinkPocketPrimitive( currPocketPrimitive );
	}
}


void Pocket::shrinkPocketPrimitive( PocketPrimitive* currPocketPrimitive )
{
    rg_dList< MBSFace* >* facesOnThisPocketPrimitive
			 = currPocketPrimitive->getMBSFaces();

    if ( facesOnThisPocketPrimitive->getSize() == 0 ) {
        return;
    }

	facesOnThisPocketPrimitive->reset4Loop();
	while ( facesOnThisPocketPrimitive->setNext4Loop() )	{
		MBSFace* currFace = facesOnThisPocketPrimitive->getEntity();

		currFace->isVisited( rg_TRUE );
	}	
		

    ////
    // shirnk할 face들을 찾아서 nodesToBeDiscarded에 넣는다.
	rg_dList< rg_dNode< MBSFace* >* > nodesToBeDiscarded;

	rg_dNode< MBSFace* >* startNode = facesOnThisPocketPrimitive->getFirstpNode();
	rg_dNode< MBSFace* >* currNode  = startNode;

	do 
	{
		MBSFace* currFace = currNode->getEntity();
        
        /*
		if ( isThisFaceOnBoundaryOfPocket( currFace ) )
		{				
			nodesToBeDiscarded.add( currNode );									
		}
        */

        rg_FLAG isThisFaceToBeDiscarded = rg_FALSE;

        rg_dList< MBSVertex* > boundingVertices;
        currFace->searchBoundingVertices( boundingVertices );

        boundingVertices.reset4Loop();
        while ( boundingVertices.setNext4Loop() ) {
            MBSVertex* currBoundingVertex = boundingVertices.getEntity();

            rg_dList< MBSFace* > incidentFaces;
            currBoundingVertex->searchIncidentFaces( incidentFaces );

            incidentFaces.reset4Loop();
            while ( incidentFaces.setNext4Loop() ) {
                if ( incidentFaces.getEntity()->isVisited() == rg_FALSE ) {
                    isThisFaceToBeDiscarded = rg_TRUE;
                    break;
                }
            }

            if ( isThisFaceToBeDiscarded == rg_TRUE ) {
                break;
            }
        }

        if ( isThisFaceToBeDiscarded == rg_TRUE ) {
            nodesToBeDiscarded.add( currNode );
        }

		currNode = currNode->getNext();
	} while( currNode != startNode );
    //
    ////

    ////
    // nodesToBeDiscarded에 있는 face 들을 제거 한다.
	nodesToBeDiscarded.reset4Loop();
	while ( nodesToBeDiscarded.setNext4Loop() )	{
		rg_dNode< MBSFace* >* currNode = nodesToBeDiscarded.getEntity();
		MBSFace*			  currFace = currNode->getEntity();

        currFace->isVisited( rg_FALSE );
		facesOnThisPocketPrimitive->kill( currNode );
	}
    //
    ////
	

	facesOnThisPocketPrimitive->reset4Loop();
	while ( facesOnThisPocketPrimitive->setNext4Loop() )	{
		MBSFace* currFace = facesOnThisPocketPrimitive->getEntity();

		currFace->isVisited( rg_FALSE );
	}
}


rg_BOOL Pocket::isThisFaceOnBoundaryOfPocket( MBSFace* currFace )
{
	rg_BOOL resultIsThisFaceOnBoundaryOfPocket = rg_FALSE;

	rg_dList< MBSFace* > incidentFaces; 
	//findIncidentRegularFacesOfGivenFace( currFace, incidentFaces );
	
	//shirink vertex-incident face
	//currFace->getIncidentFaceList( incidentFaces );
	findRegularFacesBoundedByVertexOfGivenFace( currFace, incidentFaces );

	incidentFaces.reset4Loop();
	while ( incidentFaces.setNext4Loop() )	{
		MBSFace* currIncidentFace = incidentFaces.getEntity();

        if ( currIncidentFace->isVisited() == rg_FALSE ) {
			resultIsThisFaceOnBoundaryOfPocket = rg_TRUE;
			break;
		}
	}

	return resultIsThisFaceOnBoundaryOfPocket;	
}


// geometry comparison
// manifoldized beta-shape 에서 주어진 face와 vertex-incident관계에 있는 face를 뽑아 오려면,
// 먼저 vertex의 geometry가 같은(같은 Ball*을 pointing하고 있는) vertex들을
// 먼저 찾고(navigation을 통해) 그 다음 찾아진 vertex들로부터 incident한 faces를 찾아야 한다.
void Pocket::findRegularFacesBoundedByVertexOfGivenFace( MBSFace* givenFace, rg_dList< MBSFace* >& resultIncidentFaces )
{
	rg_dList< MBSVertex* > boundingVertices;
	givenFace->searchBoundingVertices( boundingVertices );

	boundingVertices.reset4Loop();
	while ( boundingVertices.setNext4Loop() )	{
		MBSVertex* currVertex = boundingVertices.getEntity();

        // artificial한 것까지 포함해서 1차 neighbor만 뽑는다.
        rg_dList< MBSFace* > faceListincidentToCurrVertex;
		currVertex->searchIncidentFaces( faceListincidentToCurrVertex );

        faceListincidentToCurrVertex.reset4Loop();
		while ( faceListincidentToCurrVertex.setNext4Loop() )	{
			MBSFace* currIncidentFace = faceListincidentToCurrVertex.getEntity();

            if ( currIncidentFace == givenFace )    {
				continue;
            }
			
            resultIncidentFaces.addWithoutSame( currIncidentFace );	
		}
        //


        /*
        rg_dList< MBSVertex* > verticesHavingIdenticalGeometry;
		findAllVerticesHavingIdenticalGeometry( currVertex, verticesHavingIdenticalGeometry);

		verticesHavingIdenticalGeometry.reset4Loop();
		while ( verticesHavingIdenticalGeometry.setNext4Loop() )	{
			MBSVertex* currIdenticalVertex = verticesHavingIdenticalGeometry.getEntity();

			rg_dList< MBSFace* > faceListincidentToCurrVertex;
			currIdenticalVertex->searchIncidentFaces( faceListincidentToCurrVertex );

            
            
			faceListincidentToCurrVertex.reset4Loop();
			while ( faceListincidentToCurrVertex.setNext4Loop() )	{
				MBSFace* currIncidentFace = faceListincidentToCurrVertex.getEntity();

                if ( currIncidentFace == givenFace )    {
					continue;
                }
				
                resultIncidentFaces.addWithoutSame( currIncidentFace );	
			}
		}
        */
	}
}




void Pocket::findAllVerticesHavingIdenticalGeometry( MBSVertex* givenVertex, rg_dList< MBSVertex* >& resultIdenticalGeometryVertices )
{  
    //  incident한 vertex중에 geometry가 같은 것을 navigation
	rg_dList< MBSVertex* > searchStack;
	searchStack.pushBack( givenVertex );
	while ( searchStack.getSize() > 0 )	{
		MBSVertex* searchingVertex = searchStack.popBack();

		resultIdenticalGeometryVertices.add( searchingVertex );
		searchingVertex->isVisited( rg_TRUE );

		rg_dList< MBSVertex* > incidentVertex;
		searchingVertex->searchAdjacentVertices( incidentVertex );

		incidentVertex.reset4Loop();
		while ( incidentVertex.setNext4Loop() )	{
			MBSVertex* currIncidentVertex = incidentVertex.getEntity();

			if ( currIncidentVertex->isVisited() ) {
				continue;
			}

			//condition
			Ball* geometryOfGivenVertex        = givenVertex->getOriginalBetaVertex()->getBallProperty();
			Ball* geometryOfCurrIncidentVertex = currIncidentVertex->getOriginalBetaVertex()->getBallProperty();
			if ( geometryOfGivenVertex == geometryOfCurrIncidentVertex )
			{
				searchStack.pushBack( currIncidentVertex );
			}			
		}
	}

	resultIdenticalGeometryVertices.reset4Loop();
	while ( resultIdenticalGeometryVertices.setNext4Loop() )	{
		MBSVertex* currResultIncidentVertex = resultIdenticalGeometryVertices.getEntity();

		currResultIncidentVertex->isVisited( rg_FALSE );
	}   
}



void Pocket::removeMinorEdgeConnectedFaceGroup()
{
    for ( int i=0; i < m_numOfPocketPrimitives; ++i)	{							//for comparison to previous "raw pocket"
	//for ( int i=0; i < 10 && i < m_numOfPocketPrimitives; ++i)	{							//for comparison to previous "raw pocket"
		//PocketPrimitive* currPocketPrimitive = m_sortedPocketPrimitives[i];
        PocketPrimitive* currPocketPrimitive = &m_pocketPrimitives[i];

		removeMinorEdgeConnectedFaceGroupForEachPocketPrimitive( currPocketPrimitive );
	}
}


void Pocket::removeMinorEdgeConnectedFaceGroupForEachPocketPrimitive( PocketPrimitive* currPocketPrimitive )
{
    rg_dList< MBSFace* >* facesOnThisPocketPrimitive
			 = currPocketPrimitive->getMBSFaces();

	if ( facesOnThisPocketPrimitive->getSize() == 0 )
		return;

	
    facesOnThisPocketPrimitive->reset4Loop();
	while ( facesOnThisPocketPrimitive->setNext4Loop() )	{
		MBSFace* currFace = facesOnThisPocketPrimitive->getEntity();
		currFace->isVisited( rg_TRUE );
	}

    // edge-connected face group으로 나눈다.
	rg_dList< rg_dList< MBSFace* > > edgeConnectedFaceGroup;

	facesOnThisPocketPrimitive->reset4Loop();
	while ( facesOnThisPocketPrimitive->setNext4Loop() )		{
		MBSFace* currFace = facesOnThisPocketPrimitive->getEntity();

		if ( currFace->isVisited() == rg_FALSE )
			continue;

		rg_dList< MBSFace* >* currFaceGroup = edgeConnectedFaceGroup.add( rg_dList< MBSFace* >() );

		rg_dList< MBSFace* > searchStack;
		searchStack.pushBack( currFace );
		while ( searchStack.getSize() > 0 )	{
			MBSFace* searchingFace = searchStack.popBack();

			currFaceGroup->add( searchingFace );
			searchingFace->isVisited( rg_FALSE );  //처음에 visit check했던 모든 face가 결국 false로 check 된다.
			//isThisFaceOnPocket[searchingFace->getID()] = rg_FALSE;

			
			rg_dList< MBSFace* > incidentFaces;
			//findIncidentRegularFacesOfGivenFace( searchingFace, incidentFaces );
			searchingFace->searchIncidentFaces( incidentFaces );

			incidentFaces.reset4Loop();
			while ( incidentFaces.setNext4Loop() )	{
				MBSFace* currIncidentFace = incidentFaces.getEntity();

				if ( currIncidentFace->isVisited() == rg_FALSE )					
				{
					continue;
				}
				
				searchStack.pushBack( currIncidentFace );
			}
		}		
	}
	

	// 가장 큰 edge-connected face group을 찾는다.
	rg_REAL				   largestAreaOfGroup = -1;
	rg_dList< MBSFace* >*  largestEdgeConnectedFaceGroup = rg_NULL;	

	edgeConnectedFaceGroup.reset4Loop();
	while (edgeConnectedFaceGroup.setNext4Loop() )	{
		rg_dList< MBSFace* >* currEdgeConnectedFaceGroup = edgeConnectedFaceGroup.getpEntity();

		rg_REAL areaOfCurrGroup = computeSumOfAreaOnGivenFaceList( currEdgeConnectedFaceGroup );

		if ( largestAreaOfGroup < areaOfCurrGroup )
		{
			largestAreaOfGroup = areaOfCurrGroup;
			largestEdgeConnectedFaceGroup = currEdgeConnectedFaceGroup;
		}
	}

    // 가장 큰 edge-connected face group만 남기고 나머지는 지운다.
	facesOnThisPocketPrimitive->removeAll();
	facesOnThisPocketPrimitive->append( *largestEdgeConnectedFaceGroup );
}


rg_REAL Pocket::computeSumOfAreaOnGivenFaceList( rg_dList< MBSFace* >* currFaceList )
{
	rg_REAL sumOfArea = 0;

	currFaceList->reset4Loop();
	while ( currFaceList->setNext4Loop() )	{
		MBSFace* currFace = currFaceList->getEntity();

		sumOfArea += computeAreaOfMBSFace( currFace );
	}

	return sumOfArea;
}



rg_REAL Pocket::computeAreaOfMBSFace( MBSFace* givenBetaFace )
{
    rg_dList< MBSVertex* > boundingVertices;
    givenBetaFace->searchBoundingVertices( boundingVertices );

    boundingVertices.reset4Loop();
    Ball* boundingSpheres[3];
    
	rg_INT i=0;
	boundingVertices.reset4Loop();
	while ( boundingVertices.setNext4Loop() )	{
        MBSVertex* currVertex = boundingVertices.getEntity();

        boundingSpheres[i++] = currVertex->getOriginalBetaVertex()->getBallProperty();
    }

    if ( boundingSpheres[0] == boundingSpheres[1] 
      || boundingSpheres[0] == boundingSpheres[2] 
      || boundingSpheres[1] == boundingSpheres[2] 
      )                                             //area = 0
      return 0;

    rg_Point3D vector1 = boundingSpheres[1]->getGeometry().getCenter()
                       - boundingSpheres[0]->getGeometry().getCenter();
    rg_Point3D vector2 = boundingSpheres[2]->getGeometry().getCenter()
                       - boundingSpheres[0]->getGeometry().getCenter();

    rg_Point3D crossProductOfVectors = vector1.crossProduct( vector2 );
    return crossProductOfVectors.magnitude() / 2.0;
}


void Pocket::enlargePocketPrimitives()
{
    for ( int i=0; i < m_numOfPocketPrimitives; ++i)	{							//for comparison to previous "raw pocket"
	//for ( int i=0; i < 10 && i < m_numOfPocketPrimitives; ++i)	{							//for comparison to previous "raw pocket"
		//PocketPrimitive* currPocketPrimitive = m_sortedPocketPrimitives[i];
        PocketPrimitive* currPocketPrimitive = &m_pocketPrimitives[i];

		enlargePocketPrimitive( currPocketPrimitive );
	}
}


void Pocket::enlargePocketPrimitive( PocketPrimitive* currPocketPrimitive )
{
    rg_dList< MBSFace* >* faceOnThisPocket = currPocketPrimitive->getMBSFaces();

    // 현재 pocket에 해당하는 face에 모두 visit check 해둔다.
    faceOnThisPocket->reset4Loop();
	while ( faceOnThisPocket->setNext4Loop() )	{
		MBSFace* currFace = faceOnThisPocket->getEntity();
		currFace->isVisited( rg_TRUE );
	}
	

	rg_dList< MBSFace* > facesToBeAdded;
	
	faceOnThisPocket->reset4Loop();
	while ( faceOnThisPocket->setNext4Loop() )	{
		MBSFace* currFace = faceOnThisPocket->getEntity();

		rg_dList< MBSFace* > incidentFaces;
		//findIncidentRegularFacesOfGivenFace( currFace, incidentFaces );
		findRegularFacesBoundedByVertexOfGivenFace( currFace, incidentFaces );

		incidentFaces.reset4Loop();
		while ( incidentFaces.setNext4Loop() )	{
			MBSFace* currIncidentFace = incidentFaces.getEntity();

        
			if ( currIncidentFace->isVisited() == rg_FALSE )
			{
				currIncidentFace->isVisited( rg_TRUE );
				facesToBeAdded.add( currIncidentFace );
			}
		}
	}

	facesToBeAdded.reset4Loop();
	while ( facesToBeAdded.setNext4Loop() )	{
		MBSFace* currFaceToBeAdded = facesToBeAdded.getEntity();

		faceOnThisPocket->add( currFaceToBeAdded );
	}
	


	faceOnThisPocket->reset4Loop();
	while ( faceOnThisPocket->setNext4Loop() )	{
		MBSFace* currFace = faceOnThisPocket->getEntity();
		currFace->isVisited( rg_FALSE );
	}
}




void Pocket::reportAtomsOfPockets( const char* filePath, const int& numOfPockets )
{
	ofstream fout(filePath);

	rg_INT validNumOfPockets;
	if ( m_numOfPocketPrimitives < numOfPockets )
		validNumOfPockets = m_numOfPocketPrimitives;
	else
		validNumOfPockets = numOfPockets;

	cout << "AT_ID" << "\t" << "AT_NAME" << "\t" << "RE_ID" << "\t" << "RE_NAME" << endl;
	fout << "AT_ID" << "\t" << "AT_NAME" << "\t" << "RE_ID" << "\t" << "RE_NAME" << endl;

	for ( rg_INT i=0; i<validNumOfPockets; i++)
	{
		cout << "Pocket " << i+1 << "\tPocketThreshold = " << m_outerBetaVal << endl;
		fout << "Pocket " << i+1 << "\tPocketThreshold = " << m_outerBetaVal << endl;

		rg_dList<void*>* atomList 
					= m_sortedPocketPrimitives[i]->getBallPropertyList();

		atomList->reset4Loop();
		while ( atomList->setNext4Loop() )
		{
			//should be replaced
			Atom* currAtom = (Atom*)atomList->getEntity();

            
            string resName = currAtom->getResidue()->getResidueName();
            string resSeq  = StringFunctions::convertIntegerToString( currAtom->getResidue()->getSequenceNumber() );

			cout << currAtom->getSerialFromInputFile() << "\t";
			cout << currAtom->getAtomNameFromInputFile() << "\t";
			cout << resSeq << "\t";
			cout << resName << endl;
			//should be replaced


			fout << currAtom->getSerialFromInputFile() << "\t";
			fout << currAtom->getAtomNameFromInputFile() << "\t";
			fout << resSeq << "\t";
			fout << resName << endl;
			//should be replaced

		}
        cout << endl;
		fout << endl;
	}
	fout.close();
}









rg_INT  Pocket::computeNumOfVoids()
{
    BetaUniverse* qtforBeta = m_qusiTriangulationForBeta;
    rg_REAL betaValue = m_innerBetaShape.getBetaValue();

	rg_dList< BetaCell >* cellsOfQTforBeta = qtforBeta->getCellList();

	cellsOfQTforBeta->reset4Loop();
	while ( cellsOfQTforBeta->setNext4Loop() )
		cellsOfQTforBeta->getpEntity()->isVisited( rg_FALSE );
	
	rg_INT numOfVoids = 0;	

	cellsOfQTforBeta->reset4Loop();
	while ( cellsOfQTforBeta->setNext4Loop() )
	{
		BetaCell* currCell = cellsOfQTforBeta->getpEntity();

		if ( currCell->getBoundingState( betaValue ) != EXTRANEOUS_SIMPLEX )
			continue;

		if ( currCell->isVisited() )
			continue;

		rg_dList<BetaCell*> cellStack;
		cellStack.pushBack( currCell );

		rg_BOOL isThisGroupTouchedConvexHull = rg_FALSE;

		while ( cellStack.getSize() > 0 )
		{
			BetaCell* popedCell = cellStack.popBack();
			popedCell->isVisited( rg_TRUE );

			if ( popedCell->isVirtual() )
			{
				isThisGroupTouchedConvexHull = rg_TRUE;
			}

			BetaFace** boundingFace = popedCell->getFaces();

			for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++)
			{
				if ( boundingFace[i]->getBoundingState( betaValue ) != EXTRANEOUS_SIMPLEX )
					continue;

				BetaCell* nextCell = rg_NULL;
				if ( boundingFace[i]->getLeftCell() == popedCell )
					nextCell = boundingFace[i]->getRightCell();
				else
					nextCell = boundingFace[i]->getLeftCell();

				if ( nextCell->getBoundingState( betaValue ) == EXTRANEOUS_SIMPLEX
				  && nextCell->isVisited() == rg_FALSE )				
				{
					cellStack.pushBack( nextCell );
				}
			}
		}

		if ( isThisGroupTouchedConvexHull == rg_FALSE )	
			numOfVoids++;
	}

	cellsOfQTforBeta->reset4Loop();
	while ( cellsOfQTforBeta->setNext4Loop() )
		cellsOfQTforBeta->getpEntity()->isVisited( rg_FALSE );

	return numOfVoids;	
}

