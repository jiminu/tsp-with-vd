#include "rg_AugmentedBetaComplex.h"

#include <time.h>

#include <set>
#include <utility>
using namespace std;


///////////////////////////////////////////////////////////////////////////
//   
//  Constructor ...
AugmentedBetaComplex::AugmentedBetaComplex()
: m_betaValue(-1.0)
{
    r_numExtraneousCellGroupAtVertex = rg_NULL;
    r_numVerticesDifferentNumConesInputAndTT_D_V = 0;
}



AugmentedBetaComplex::~AugmentedBetaComplex()
{
    if ( r_numExtraneousCellGroupAtVertex != rg_NULL )
		delete [] r_numExtraneousCellGroupAtVertex;

}




///////////////////////////////////////////////////////////////////////////
//   
//  GET functions
BetaUniverse* AugmentedBetaComplex::getBetaUniverse() const
{
    return m_betaUniverse;
}


rg_REAL AugmentedBetaComplex::getBetaValue() const
{
    return m_betaValue;
}




rg_INT AugmentedBetaComplex::getNumAugmentedBCCells() const
{
    return m_augmentedBCCells.getSize();
}



rg_INT AugmentedBetaComplex::getNumAugmentedBCFaces() const
{
    return m_augmentedBCFaces.getSize();
}



rg_INT AugmentedBetaComplex::getNumAugmentedBCEdges() const
{
    return m_augmentedBCEdges.getSize();
}



rg_INT AugmentedBetaComplex::getNumAugmentedBCVertices() const
{
    return m_augmentedBCVertices.getSize();
}





rg_INT AugmentedBetaComplex::getNumArtificialAugmentedBCCells() const
{
    return m_artificialAugmentedBCCells.getSize();
}



rg_INT AugmentedBetaComplex::getNumArtificialAugmentedBCFaces() const
{
    return m_artificialAugmentedBCFaces.getSize();
}



rg_INT AugmentedBetaComplex::getNumArtificialAugmentedBCEdges() const
{
    return m_artificialAugmentedBCEdges.getSize();
}



rg_INT AugmentedBetaComplex::getNumArtificialAugmentedBCVertices() const
{
    return m_artificialAugmentedBCVertices.getSize();
}




rg_INT AugmentedBetaComplex::getNumAllAugmentedBCCells() const
{
    return m_augmentedBCCells.getSize() + m_artificialAugmentedBCCells.getSize();
}



rg_INT AugmentedBetaComplex::getNumAllAugmentedBCFaces() const
{
    return m_augmentedBCFaces.getSize() + m_artificialAugmentedBCFaces.getSize();
}



rg_INT AugmentedBetaComplex::getNumAllAugmentedBCEdges() const
{
    return m_augmentedBCEdges.getSize() + m_artificialAugmentedBCEdges.getSize();
}



rg_INT AugmentedBetaComplex::getNumAllAugmentedBCVertices() const
{
    return m_augmentedBCVertices.getSize() + m_artificialAugmentedBCVertices.getSize();
}



rg_dList<AugmentedBCCell>* AugmentedBetaComplex::getAugmentedBCCells()
{
    return &m_augmentedBCCells;
}



rg_dList<AugmentedBCFace>* AugmentedBetaComplex::getAugmentedBCFaces()
{
    return &m_augmentedBCFaces;
}



rg_dList<AugmentedBCEdge>* AugmentedBetaComplex::getAugmentedBCEdges()
{
    return &m_augmentedBCEdges;
}



rg_dList<AugmentedBCVertex>* AugmentedBetaComplex::getAugmentedBCVertices()
{
    return &m_augmentedBCVertices;
}



rg_dList<AugmentedBCCell>* AugmentedBetaComplex::getArtificialAugmentedBCCells()
{
    return &m_artificialAugmentedBCCells;
}



rg_dList<AugmentedBCFace>* AugmentedBetaComplex::getArtificialAugmentedBCFaces()
{
    return &m_artificialAugmentedBCFaces;
}



rg_dList<AugmentedBCEdge>* AugmentedBetaComplex::getArtificialAugmentedBCEdges()
{
    return &m_artificialAugmentedBCEdges;
}



rg_dList<AugmentedBCVertex>* AugmentedBetaComplex::getArtificialAugmentedBCVertices()
{
    return &m_artificialAugmentedBCVertices;
}




rg_dList<Ball>* AugmentedBetaComplex::getBalls()
{
    return &m_balls;
}






///////////////////////////////////////////////////////////////////////////
//   
//  SET functions
void AugmentedBetaComplex::setBetaUniverse( BetaUniverse* betaUniverse )
{
    m_betaUniverse = betaUniverse;
}

void AugmentedBetaComplex::setBetaValue(const rg_REAL& value)
{
    m_betaValue = value;
}




// MAKE augmented beta-complex
// void AugmentedBetaComplex::construct( const rg_REAL& betaValue, BetaUniverse* betaUniverse, const rg_INT& constructOption )
// {
//     if ( betaUniverse->getNumVertices() == 0 )    //BetaComplex가 계산되어 있지 않으면
//         return;
//     
//     if ( m_augmentedBCVertices.getSize() != 0 )
//         clearMyself(); 
// 
// 
// 
//     setBetaUniverse( betaUniverse );
// 	setBetaValue( betaValue );
// 
//     duplicateOriginalQuasiTriangulationForBetaToAugmentedBetaComplex(betaUniverse);
// 
// 
//     if ( constructOption == CONSTRUCT_WITH_WHOLE_BETA_SHAPE ) {
//         //do nothing.    
//     }
//     else if ( constructOption == CONSTRUCT_WITH_REGULAR_BETA_SHAPE ) {
//         setSingularEdgesToFalseExterior();
// 
//         setSingularFacesToFalseExterior();    
//     }
//     
// 
//     indexExteriorRegionOfEdgesAndFaces();
// 
// 
//     //For time checking
//     clock_t start;    
//     start  = clock();
// 	r_numExtraneousCellGroupAtVertex = new rg_INT[betaUniverse->getVertexList()->getSize()+1];
// 	
// 
//     rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices = 
//                                         new rg_dList<NonManifoldCorner*>[ m_augmentedBCVertices.getSize() * 2 ];
//     rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges = 
//                                         new rg_dList<NonManifoldCorner*>[ m_augmentedBCEdges.getSize() * 2 ];
// 
//     rg_dList< NonManifoldCorner > nonManifoldCornerList;
//     searchNonManifoldConfigurationsForAllVertices(ptToNMCornerOfNMVertices, nonManifoldCornerList);
//     searchNonManifoldConfigurationsForAllEdges(   ptToNMCornerOfNMEdges,    nonManifoldCornerList);
// 
// 	reportNumNMCorners(nonManifoldCornerList);
// 
//     detachAndRemoveExtraneousSimplexes();
// 
//     augment(ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges, nonManifoldCornerList);
//         
//     
//     if ( ptToNMCornerOfNMVertices != rg_NULL )
//         delete [] ptToNMCornerOfNMVertices;
//     if ( ptToNMCornerOfNMEdges != rg_NULL )
//         delete [] ptToNMCornerOfNMEdges;
//     
//     //For time checking
//     clock_t finish = clock();
// 
//     //reportResultOfManifoldization( (double)(finish - start) / CLOCKS_PER_SEC );  
//     
// }









///////////////////////////////////////////////////////////////////////////
//   
//  MAKE augmented beta-complex
//
//  The following member function is fixed by Youngsong Cho 2011-04-23
//void AugmentedBetaComplex::construct( const rg_REAL& betaValue, BetaUniverse* betaUniverse, const rg_INT& constructOption )
//{
//    if (    betaUniverse->getNumVertices() == 0 
//         || betaUniverse->getNumEdges() == 0 ) {
//        return;
//    }
//    
//    if ( m_augmentedBCVertices.getSize() != 0 )
//        clearMyself(); 
//
//
//
//    setBetaUniverse( betaUniverse );
//	setBetaValue( betaValue );
//
//    
//    duplicateOriginalQuasiTriangulationForBetaToAugmentedBetaComplex(betaUniverse);
//
//    
//    setEdgeAsSingularIfEdgeIncidentToExtraneousFacesOnly();
//
//    
//    if (      constructOption == CONSTRUCT_WITH_WHOLE_BETA_SHAPE ) {
//        //do nothing.    
//    }
//    else if ( constructOption == CONSTRUCT_WITH_REGULAR_BETA_SHAPE ) {
//        setSingularEdgesToFalseExterior();
//
//        setSingularFacesToFalseExterior();
//        
//        setSingularVerticesToFalseExterior();
//    }
//    
//
//    
//    indexExteriorRegionOfEdgesAndFaces();
//
//
//
//
//
//    rg_INT maxIndexOfQTVertex = betaUniverse->getVertexList()->getSize() + 1;
//    r_numExtraneousCellGroupAtVertex = new rg_INT[maxIndexOfQTVertex];
//    for ( rg_INT i=0; i<maxIndexOfQTVertex; i++ ) {
//        r_numExtraneousCellGroupAtVertex[i] = 0;
//    }
//
//
//
//
//    rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices = 
//                                        new rg_dList<NonManifoldCorner*>[ m_augmentedBCVertices.getSize() ];
//    rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges = 
//                                        new rg_dList<NonManifoldCorner*>[ m_augmentedBCEdges.getSize() ];
//
//    rg_dList< NonManifoldCorner > nonManifoldCornerList;
//    searchNonManifoldConfigurationsForAllVertices(ptToNMCornerOfNMVertices, nonManifoldCornerList);
//    searchNonManifoldConfigurationsForAllEdges(   ptToNMCornerOfNMEdges,    nonManifoldCornerList);
//
//	
//    detachAndRemoveExtraneousSimplexes();
//
//
//    augment(ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges, nonManifoldCornerList);
//        
//    
//    if ( ptToNMCornerOfNMVertices != rg_NULL )
//        delete [] ptToNMCornerOfNMVertices;
//    if ( ptToNMCornerOfNMEdges != rg_NULL )
//        delete [] ptToNMCornerOfNMEdges;   
//    
//}




void AugmentedBetaComplex::construct( const rg_REAL& betaValue, BetaUniverse* betaUniverse, const rg_INT& constructOption )
{
    if ( betaUniverse->getNumVertices() == 0 || betaUniverse->getNumEdges() == 0 ) {
        return;
    }
    
    if ( m_augmentedBCVertices.getSize() != 0 ) {
        clearMyself(); 
    }


    m_betaUniverse = betaUniverse;
	m_betaValue    = betaValue;

    
    //  prepare to augment beta-complex.
    duplicateBetaComplexToAugmentedBetaComplex();
    setEdgeAsSingularIfEdgeIncidentToExtraneousFacesOnly();
    
    if (      constructOption == CONSTRUCT_WITH_WHOLE_BETA_SHAPE ) {
        //do nothing.    
    }
    else if ( constructOption == CONSTRUCT_WITH_REGULAR_BETA_SHAPE ) {
        setSingularEdgesToFalseExterior();
        setSingularFacesToFalseExterior();
        setSingularVerticesToFalseExterior();
    }
    
    indexExteriorRegionOfEdgesAndFaces();


    //  search non-manifold configurations
    rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices = new rg_dList<NonManifoldCorner*>[ m_augmentedBCVertices.getSize() ];
    rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges    = new rg_dList<NonManifoldCorner*>[ m_augmentedBCEdges.getSize() ];

    rg_dList< NonManifoldCorner > nonManifoldCornerList;
    searchNonManifoldConfigurationsForAllVertices(ptToNMCornerOfNMVertices, nonManifoldCornerList);
    searchNonManifoldConfigurationsForAllEdges(   ptToNMCornerOfNMEdges,    nonManifoldCornerList);


	//  detach exterior simplexes
    detachAndRemoveExtraneousSimplexes();


    //  augment
    augment(ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges, nonManifoldCornerList);
        
    
    if ( ptToNMCornerOfNMVertices != rg_NULL )
        delete [] ptToNMCornerOfNMVertices;
    if ( ptToNMCornerOfNMEdges != rg_NULL )
        delete [] ptToNMCornerOfNMEdges;   
    
}




//  Youngsong Cho 2011-04-22 
void AugmentedBetaComplex::duplicateBetaComplexToAugmentedBetaComplex()
{
    map<BetaVertex*, AugmentedBCVertex*> mapBetaVertexToAugmentedBetaVertex;
    map<BetaEdge*,   AugmentedBCEdge*>   mapBetaEdgeToAugmentedBetaEdge;
    map<BetaFace*,   AugmentedBCFace*>   mapBetaFaceToAugmentedBetaFace;
    map<BetaCell*,   AugmentedBCCell*>   mapBetaCellToAugmentedBetaCell;

    createAugmentedBetaSimplexesAndRegisterToMap(  mapBetaVertexToAugmentedBetaVertex,
                                                   mapBetaEdgeToAugmentedBetaEdge,
                                                   mapBetaFaceToAugmentedBetaFace,
                                                   mapBetaCellToAugmentedBetaCell );


    duplicateBetaVerticesToAugmentedBCVertices(    mapBetaVertexToAugmentedBetaVertex,
                                                   mapBetaEdgeToAugmentedBetaEdge,
                                                   mapBetaCellToAugmentedBetaCell );

    duplicateBetaEdgesToAugmentedBCEdges(          mapBetaVertexToAugmentedBetaVertex,
                                                   mapBetaEdgeToAugmentedBetaEdge,
                                                   mapBetaFaceToAugmentedBetaFace);

    duplicateBetaCellsToAugmentedBCCells(          mapBetaVertexToAugmentedBetaVertex,
                                                   mapBetaFaceToAugmentedBetaFace,
                                                   mapBetaCellToAugmentedBetaCell );
    
    duplicateBetaFacesToAugmentedBCFaces(          mapBetaEdgeToAugmentedBetaEdge,
                                                   mapBetaFaceToAugmentedBetaFace,
                                                   mapBetaCellToAugmentedBetaCell );

}


void AugmentedBetaComplex::createAugmentedBetaSimplexesAndRegisterToMap( map<BetaVertex*, AugmentedBCVertex*>& mapBetaVertexToAugmentedBetaVertex,
                                                                         map<BetaEdge*,   AugmentedBCEdge*>&   mapBetaEdgeToAugmentedBetaEdge,
                                                                         map<BetaFace*,   AugmentedBCFace*>&   mapBetaFaceToAugmentedBetaFace,
                                                                         map<BetaCell*,   AugmentedBCCell*>&   mapBetaCellToAugmentedBetaCell )
{
    rg_INT i=0;
    rg_dList<BetaVertex>* betaVertices = m_betaUniverse->getVertexList();
    betaVertices->reset4Loop();
    while ( betaVertices->setNext4Loop() ) {
        BetaVertex*        currBetaVtx          = betaVertices->getpEntity();
        AugmentedBCVertex* currAugmentedBetaVtx = m_augmentedBCVertices.add( AugmentedBCVertex( i++ ) );            
        //AugmentedBCVertex* currAugmentedBetaVtx = m_augmentedBCVertices.add( AugmentedBCVertex( currBetaVtx->getID() ) );            

        currAugmentedBetaVtx->setOriginalBetaVertex( currBetaVtx );

        mapBetaVertexToAugmentedBetaVertex.insert( make_pair(currBetaVtx, currAugmentedBetaVtx) ); 
    }


    i=0;
    rg_dList<BetaEdge>* betaEdges = m_betaUniverse->getEdgeList();
    betaEdges->reset4Loop();
    while ( betaEdges->setNext4Loop() ) {
        BetaEdge*        currBetaEdge          = betaEdges->getpEntity();
        AugmentedBCEdge* currAugmentedBetaEdge = m_augmentedBCEdges.add( AugmentedBCEdge( i++ ) );            
        //AugmentedBCEdge* currAugmentedBetaEdge = m_augmentedBCEdges.add( AugmentedBCEdge( currBetaEdge->getID() ) );            

        mapBetaEdgeToAugmentedBetaEdge.insert( make_pair(currBetaEdge, currAugmentedBetaEdge) ); 
    }


    i=0;
    rg_dList<BetaFace>* betaFaces = m_betaUniverse->getFaceList();
    betaFaces->reset4Loop();
    while ( betaFaces->setNext4Loop() ) {
        BetaFace*        currBetaFace          = betaFaces->getpEntity();
        AugmentedBCFace* currAugmentedBetaFace = m_augmentedBCFaces.add( AugmentedBCFace( i++ ) );            
        //AugmentedBCFace* currAugmentedBetaFace = m_augmentedBCFaces.add( AugmentedBCFace( currBetaFace->getID() ) );            

        mapBetaFaceToAugmentedBetaFace.insert( make_pair(currBetaFace, currAugmentedBetaFace) ); 
    }


    i=0;
    rg_dList<BetaCell>* betaCells = m_betaUniverse->getCellList();
    betaCells->reset4Loop();
    while ( betaCells->setNext4Loop() ) {
        BetaCell*        currBetaCell          = betaCells->getpEntity();
        AugmentedBCCell* currAugmentedBetaCell = m_augmentedBCCells.add( AugmentedBCCell( i++ ) );            
        //AugmentedBCCell* currAugmentedBetaCell = m_augmentedBCCells.add( AugmentedBCCell( currBetaCell->getID() ) );            

        mapBetaCellToAugmentedBetaCell.insert( make_pair(currBetaCell, currAugmentedBetaCell) ); 
    }
}



void AugmentedBetaComplex::duplicateBetaVerticesToAugmentedBCVertices( const map<BetaVertex*, AugmentedBCVertex*>& mapBetaVertexToAugmentedBetaVertex,
                                                                       const map<BetaEdge*,   AugmentedBCEdge*>&   mapBetaEdgeToAugmentedBetaEdge,
                                                                       const map<BetaCell*,   AugmentedBCCell*>&   mapBetaCellToAugmentedBetaCell )
{
    rg_dList<BetaVertex>* betaVertices = m_betaUniverse->getVertexList();
    betaVertices->reset4Loop();
    while ( betaVertices->setNext4Loop() ) {
        BetaVertex*        currBetaVtx          = betaVertices->getpEntity();
        AugmentedBCVertex* currAugmentedBetaVtx = mapBetaVertexToAugmentedBetaVertex.find( currBetaVtx )->second;            

        currAugmentedBetaVtx->setFirstEdge( mapBetaEdgeToAugmentedBetaEdge.find( currBetaVtx->getFirstEdge() )->second );
        if ( currBetaVtx->getFirstCell() == rg_NULL ) {
            currAugmentedBetaVtx->setFirstCell( rg_NULL );
        }
        else {
            currAugmentedBetaVtx->setFirstCell( mapBetaCellToAugmentedBetaCell.find( currBetaVtx->getFirstCell() )->second );
        }
        

		if ( !currBetaVtx->isVirtual() ) {
			currAugmentedBetaVtx->setBetaSpan( currBetaVtx->getBetaSpan() );			

			Ball* ptOfBall = m_balls.add( *(currBetaVtx->getBallProperty()) );
			currAugmentedBetaVtx->setBallProperty( ptOfBall ); 
		}        
		
		currAugmentedBetaVtx->setBoundingState( currBetaVtx->getBoundingState( m_betaValue ) );		
        currAugmentedBetaVtx->setOriginalBetaVertex( currBetaVtx );
        currAugmentedBetaVtx->setIsArtificial( rg_FALSE );		
		currAugmentedBetaVtx->isVisited( rg_FALSE );
    }    
}



void AugmentedBetaComplex::duplicateBetaEdgesToAugmentedBCEdges( const map<BetaVertex*, AugmentedBCVertex*>& mapBetaVertexToAugmentedBetaVertex,
                                                                 const map<BetaEdge*,   AugmentedBCEdge*>&   mapBetaEdgeToAugmentedBetaEdge,
                                                                 const map<BetaFace*,   AugmentedBCFace*>&   mapBetaFaceToAugmentedBetaFace)
{
    rg_dList<BetaEdge>* betaEdges = m_betaUniverse->getEdgeList();

    betaEdges->reset4Loop();
    while( betaEdges->setNext4Loop() )  {
        BetaEdge*        currBetaEdge        = betaEdges->getpEntity();;
        AugmentedBCEdge* currAugmentedBCEdge = mapBetaEdgeToAugmentedBetaEdge.find( currBetaEdge )->second;


        currAugmentedBCEdge->setStartVertex( mapBetaVertexToAugmentedBetaVertex.find( currBetaEdge->getStartVertex() )->second );
        currAugmentedBCEdge->setEndVertex(   mapBetaVertexToAugmentedBetaVertex.find( currBetaEdge->getEndVertex() )->second   );
        currAugmentedBCEdge->setFirstFace(   mapBetaFaceToAugmentedBetaFace.find(     currBetaEdge->getFirstFace() )->second   );
		
		rg_dList<BetaFace*>* smallWorlds = currBetaEdge->getSmallWorlds();
		BetaFace*			 currSmallWorldFace = rg_NULL;
		smallWorlds->reset4Loop();
		while ( smallWorlds->setNext4Loop() ) {
			currSmallWorldFace = smallWorlds->getEntity();
            currAugmentedBCEdge->getSmallWorlds()->add( mapBetaFaceToAugmentedBetaFace.find( currSmallWorldFace )->second );
		}

		if ( !currAugmentedBCEdge->isVirtual() ) {
			currAugmentedBCEdge->setBetaSpan( currBetaEdge->getBetaSpan() );
		}
		
		currAugmentedBCEdge->setBoundingState( currBetaEdge->getBoundingState( m_betaValue ) );
		currAugmentedBCEdge->setIsArtificial( rg_FALSE );
		currAugmentedBCEdge->isVisited( rg_FALSE );
    }
}



void AugmentedBetaComplex::duplicateBetaFacesToAugmentedBCFaces( const map<BetaEdge*,   AugmentedBCEdge*>&   mapBetaEdgeToAugmentedBetaEdge,
                                                                 const map<BetaFace*,   AugmentedBCFace*>&   mapBetaFaceToAugmentedBetaFace,
                                                                 const map<BetaCell*,   AugmentedBCCell*>&   mapBetaCellToAugmentedBetaCell )
{
    rg_dList<BetaFace>* betaFaces = m_betaUniverse->getFaceList();

    betaFaces->reset4Loop();
    while( betaFaces->setNext4Loop() )   {
        BetaFace*        currBetaFace        = betaFaces->getpEntity();        
        AugmentedBCFace* currAugmentedBCFace = mapBetaFaceToAugmentedBetaFace.find( currBetaFace )->second;

        if ( currBetaFace->getRightCell() == rg_NULL ) { //is currBetaFace an isolated face?
            currAugmentedBCFace->setRightCell( rg_NULL );
        }
        else {
            currAugmentedBCFace->setRightCell( mapBetaCellToAugmentedBetaCell.find( currBetaFace->getRightCell() )->second );
        }
        
        if ( currBetaFace->getLeftCell() == rg_NULL ) {
            currAugmentedBCFace->setLeftCell( rg_NULL );
        }
        else {
            currAugmentedBCFace->setLeftCell(  mapBetaCellToAugmentedBetaCell.find( currBetaFace->getLeftCell() )->second );
        }

        

		for (rg_INT i = 0; i < EIWDS_NUM_EDGE_ON_FACE; i++)        {
            currAugmentedBCFace->setEdge(i, currBetaFace->getEdgeOrientation(i), 
                                            mapBetaEdgeToAugmentedBetaEdge.find( currBetaFace->getEdge(i) )->second );
        }

		if ( !currAugmentedBCFace->isVirtual() ) {
			currAugmentedBCFace->setBetaSpan( currBetaFace->getBetaSpan() );
		}
        
        currAugmentedBCFace->setBoundingState( currBetaFace->getBoundingState( m_betaValue ) );
		currAugmentedBCFace->setIsArtificial( rg_FALSE );
		currAugmentedBCFace->isVisited( rg_FALSE );

        if ( currAugmentedBCFace->getID() == 26 )
            int stop = 0;

        //reverse orientations if leftCell is "Virtual" or "Extraneous" <-regular face의 경우 항상 right cell이 rg_NULL이도록 setting
        AugmentedBCCell* leftCell = (AugmentedBCCell*)currAugmentedBCFace->getLeftCell();
        
        if ( leftCell != rg_NULL ) {
            rg_INT leftCellBoundingState= leftCell->getBoundingState();
		    if (    leftCellBoundingState == EXTRANEOUS_SIMPLEX 
                 || leftCellBoundingState == FALSE_EXTERIOR_SIMPLEX ) {
                currAugmentedBCFace->reverse();            
            }     
        }


        if ( currAugmentedBCFace->getID() == 26 )
            int stop = 0;

	}
}



void AugmentedBetaComplex::duplicateBetaCellsToAugmentedBCCells( const map<BetaVertex*, AugmentedBCVertex*>& mapBetaVertexToAugmentedBetaVertex,
                                                                 const map<BetaFace*,   AugmentedBCFace*>&   mapBetaFaceToAugmentedBetaFace,
                                                                 const map<BetaCell*,   AugmentedBCCell*>&   mapBetaCellToAugmentedBetaCell )
{
    rg_dList<BetaCell>* betaCells = m_betaUniverse->getCellList();

    betaCells->reset4Loop();
    while ( betaCells->setNext4Loop() ) {
        BetaCell*        currBetaCell        = betaCells->getpEntity();
        AugmentedBCCell* currAugmentedBCCell = mapBetaCellToAugmentedBetaCell.find( currBetaCell )->second;

		rg_INT i = 0;
        for (i = 0; i < EIWDS_NUM_FACE_ON_CELL; i++)  {
            currAugmentedBCCell->setFace(i, mapBetaFaceToAugmentedBetaFace.find( currBetaCell->getFace(i) )->second );
        }

		for (i = 0; i < EIWDS_NUM_VERTEX_ON_CELL; i++)  {
            currAugmentedBCCell->setVertex(i, mapBetaVertexToAugmentedBetaVertex.find( currBetaCell->getVertex(i) )->second );
        }

		if ( !currAugmentedBCCell->isVirtual() ) {
			currAugmentedBCCell->setBetaSpan( currBetaCell->getBetaSpan() );
			currAugmentedBCCell->setMinTangentSphere( currBetaCell->getMinTangentSphere() );
		}

        currAugmentedBCCell->setBoundingState( currBetaCell->getBoundingState( m_betaValue ) );
		currAugmentedBCCell->setIsArtificial( rg_FALSE );
		currAugmentedBCCell->isVisited( rg_FALSE );
    }
}



void AugmentedBetaComplex::writeDuplicate()
{
    ofstream fout;
    fout.open("duplicate_AugmentedBC.txt");

    fout << "Augmented beta-vertices: " << endl;
    m_augmentedBCVertices.reset4Loop();
    while ( m_augmentedBCVertices.setNext4Loop() ) {
        AugmentedBCVertex* currVtx = m_augmentedBCVertices.getpEntity();

        fout << currVtx->getID() << "\t";
        if ( currVtx->getFirstEdge() != rg_NULL )
            fout << currVtx->getFirstEdge()->getID() << "\t";
        if ( currVtx->getFirstCell() != rg_NULL )
            fout << currVtx->getFirstCell()->getID() << "\t";

        fout << endl;
    }


    fout << "Augmented beta-edges: " << endl;
    m_augmentedBCEdges.reset4Loop();
    while ( m_augmentedBCEdges.setNext4Loop() ) {
        AugmentedBCEdge* currEdge = m_augmentedBCEdges.getpEntity();

        fout << currEdge->getID() << "\t";
        fout << currEdge->getIndexOfExteriorRegion() << "\t";
        fout << currEdge->getStartVertex()->getID() << "\t";
        fout << currEdge->getEndVertex()->getID() << "\t";
        fout << currEdge->getFirstFace()->getID() << "\t";

        fout << endl;
    }


    fout << "Augmented beta-faces: " << endl;
    m_augmentedBCFaces.reset4Loop();
    while ( m_augmentedBCFaces.setNext4Loop() ) {
        AugmentedBCFace* currFace = m_augmentedBCFaces.getpEntity();

        fout << currFace->getID() << "\t";
        if ( currFace->getLeftCell() != rg_NULL )
            fout << currFace->getLeftCell()->getID() << "(" << ((AugmentedBCCell*) currFace->getLeftCell())->getBoundingState() << ")\t";
        if ( currFace->getRightCell() != rg_NULL )
            fout << currFace->getRightCell()->getID() << "(" << ((AugmentedBCCell*) currFace->getRightCell())->getBoundingState() << ")\t";
        fout << currFace->getEdge(0)->getID() << "\t";
        fout << currFace->getEdge(1)->getID() << "\t";
        fout << currFace->getEdge(2)->getID() << "\t";
        fout << currFace->getIndexOfExteriorRegionToLeftCell() << "\t";
        fout << currFace->getIndexOfExteriorRegionToRightCell() << "\t";

        fout << endl;
    }


    fout << "Augmented beta-cells: " << endl;
    m_augmentedBCCells.reset4Loop();
    while ( m_augmentedBCCells.setNext4Loop() ) {
        AugmentedBCCell* currCell = m_augmentedBCCells.getpEntity();

        fout << currCell->getID() << "\t";
        fout << currCell->getVertex(0)->getID() << "\t";
        fout << currCell->getVertex(1)->getID() << "\t";
        fout << currCell->getVertex(2)->getID() << "\t";
        fout << currCell->getVertex(3)->getID() << "\t";

        fout << currCell->getFace(0)->getID() << "\t";
        fout << currCell->getFace(1)->getID() << "\t";
        fout << currCell->getFace(2)->getID() << "\t";
        fout << currCell->getFace(3)->getID() << "\t";

        fout << endl;
    }



    fout << endl << endl;
    fout << "------------ARTIFICIAL AUGMENTED BETA SIMPLEXES ------------------------" << endl;
    fout << endl << endl;


    fout << "Artificial augmented beta-vertices: " << endl;
    m_artificialAugmentedBCVertices.reset4Loop();
    while ( m_artificialAugmentedBCVertices.setNext4Loop() ) {
        AugmentedBCVertex* currVtx = m_artificialAugmentedBCVertices.getpEntity();

        fout << currVtx->getID() << "\t";
        if ( currVtx->getFirstEdge() != rg_NULL )
            fout << currVtx->getFirstEdge()->getID() << "\t";
        if ( currVtx->getFirstCell() != rg_NULL )
            fout << currVtx->getFirstCell()->getID() << "\t";

        fout << endl;
    }


    fout << "Artificial augmented beta-edges: " << endl;
    m_artificialAugmentedBCEdges.reset4Loop();
    while ( m_artificialAugmentedBCEdges.setNext4Loop() ) {
        AugmentedBCEdge* currEdge = m_artificialAugmentedBCEdges.getpEntity();

        fout << currEdge->getID() << "\t";
        fout << currEdge->getIndexOfExteriorRegion() << "\t";
        fout << currEdge->getStartVertex()->getID() << "\t";
        fout << currEdge->getEndVertex()->getID() << "\t";
        fout << currEdge->getFirstFace()->getID() << "\t";

        fout << endl;
    }


    fout << "Artificial augmented beta-faces: " << endl;
    m_artificialAugmentedBCFaces.reset4Loop();
    while ( m_artificialAugmentedBCFaces.setNext4Loop() ) {
        AugmentedBCFace* currFace = m_artificialAugmentedBCFaces.getpEntity();

        fout << currFace->getID() << "\t";
        if ( currFace->getLeftCell() != rg_NULL )
            fout << currFace->getLeftCell()->getID() << "(" << ((AugmentedBCCell*) currFace->getLeftCell())->getBoundingState() << ")\t";
        if ( currFace->getRightCell() != rg_NULL )
            fout << currFace->getRightCell()->getID() << "(" << ((AugmentedBCCell*) currFace->getRightCell())->getBoundingState() << ")\t";
        fout << currFace->getEdge(0)->getID() << "\t";
        fout << currFace->getEdge(1)->getID() << "\t";
        fout << currFace->getEdge(2)->getID() << "\t";
        fout << currFace->getIndexOfExteriorRegionToLeftCell() << "\t";
        fout << currFace->getIndexOfExteriorRegionToRightCell() << "\t";

        fout << endl;
    }


    fout << "Artificial augmented beta-cells: " << endl;
    m_artificialAugmentedBCCells.reset4Loop();
    while ( m_artificialAugmentedBCCells.setNext4Loop() ) {
        AugmentedBCCell* currCell = m_artificialAugmentedBCCells.getpEntity();

        fout << currCell->getID() << "\t";
        fout << currCell->getVertex(0)->getID() << "\t";
        fout << currCell->getVertex(1)->getID() << "\t";
        fout << currCell->getVertex(2)->getID() << "\t";
        fout << currCell->getVertex(3)->getID() << "\t";

        fout << currCell->getFace(0)->getID() << "\t";
        fout << currCell->getFace(1)->getID() << "\t";
        fout << currCell->getFace(2)->getID() << "\t";
        fout << currCell->getFace(3)->getID() << "\t";

        fout << endl;
    }
}



void AugmentedBetaComplex::writeNonManifoldCorner( rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices, 
                                                   rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges )
{
    ofstream fout;
    fout.open("non-manifold corners.txt");

    rg_INT numNMCornerOfNMVertices = m_augmentedBCVertices.getSize();
    rg_INT numNMCornerOfNMEdges    = m_augmentedBCEdges.getSize();

    fout << "Non-manifold corners at non-manifold vertices" << endl;

    rg_INT i=0; 
    for ( i=0; i<numNMCornerOfNMVertices; i++ ) {
        if ( ptToNMCornerOfNMVertices[i].getSize()==0 ) {
            continue;
        }

        fout << i << "_non-manifold vertex" << endl;
        rg_INT ID = 0;
        ptToNMCornerOfNMVertices[i].reset4Loop();
        while ( ptToNMCornerOfNMVertices[i].setNext4Loop() ) {
            NonManifoldCorner* nmConner = ptToNMCornerOfNMVertices[i].getEntity();


            fout << "\t" << ID++ << "\t" << nmConner->getTypeOfNonManifoldConfiguration() << "\t";

            PointerToSimplex* nmSimplex[3] = {nmConner->getpNMSimplex(), nmConner->getpFirstSimplex(), nmConner->getpSecondSimplex()};

            for ( rg_INT j=0; j<3; j++ ) {
                switch ( nmSimplex[j]->getTypeOfSimplex() ) {
                    case VERTEX_SIMPLEX:
                        fout << "v_" << ((AugmentedBCVertex*)nmSimplex[j]->getSimplex())->getID() << "\t";
                        break;
                    case EDGE_SIMPLEX:
                        fout << "e_" << ((AugmentedBCEdge*)nmSimplex[j]->getSimplex())->getID() << "\t";
                        break;
                    case FACE_SIMPLEX:
                        fout << "f_" << ((AugmentedBCFace*)nmSimplex[j]->getSimplex())->getID() << "\t";
                        break;
                    case CELL_SIMPLEX:
                        fout << "c_" << ((AugmentedBCCell*)nmSimplex[j]->getSimplex())->getID() << "\t";
                        break;
                }
            }
            fout << endl;
        }
    }


    fout << "Non-manifold corners at non-manifold edges" << endl;

    for ( i=0; i<numNMCornerOfNMEdges; i++ ) {
        if ( ptToNMCornerOfNMEdges[i].getSize()==0 ) {
            continue;
        }

        fout << i << "_non-manifold edge" << endl;
        rg_INT ID = 0;
        ptToNMCornerOfNMEdges[i].reset4Loop();
        while ( ptToNMCornerOfNMEdges[i].setNext4Loop() ) {
            NonManifoldCorner* nmConner = ptToNMCornerOfNMEdges[i].getEntity();


            fout << "\t" << ID++ << "\t" << nmConner->getTypeOfNonManifoldConfiguration() << "\t";

            PointerToSimplex* nmSimplex[3] = {nmConner->getpNMSimplex(), nmConner->getpFirstSimplex(), nmConner->getpSecondSimplex()};

            for ( rg_INT j=0; j<3; j++ ) {
                switch ( nmSimplex[j]->getTypeOfSimplex() ) {
                    case VERTEX_SIMPLEX:
                        fout << "v_" << ((AugmentedBCVertex*)nmSimplex[j]->getSimplex())->getID() << "\t";
                        break;
                    case EDGE_SIMPLEX:
                        fout << "e_" << ((AugmentedBCEdge*)nmSimplex[j]->getSimplex())->getID() << "\t";
                        break;
                    case FACE_SIMPLEX:
                        fout << "f_" << ((AugmentedBCFace*)nmSimplex[j]->getSimplex())->getID() << "\t";
                        break;
                    case CELL_SIMPLEX:
                        fout << "c_" << ((AugmentedBCCell*)nmSimplex[j]->getSimplex())->getID() << "\t";
                        break;
                }
            }
            fout << endl;
        }
    }
}



void AugmentedBetaComplex::adjustIDOfAugmentedBetaComplex()
{
	adjustIDOfAugmentedBCVertices();
	adjustIDOfAugmentedBCEdges();
	adjustIDOfAugmentedBCFaces();
	adjustIDOfAugmentedBCCells();

	adjustIDOfBalls();
}


void AugmentedBetaComplex::adjustIDOfAugmentedBCVertices()
{
    int i = 0;
    AugmentedBCVertex* currVertex = rg_NULL;
    m_augmentedBCVertices.reset4Loop();
    while (m_augmentedBCVertices.setNext4Loop() )
    {
        currVertex = m_augmentedBCVertices.getpEntity();
        currVertex->setID(i);
        i++;        
    }

    m_artificialAugmentedBCVertices.reset4Loop();
    while (m_artificialAugmentedBCVertices.setNext4Loop() )
    {
        currVertex = m_artificialAugmentedBCVertices.getpEntity();
        currVertex->setID(i);
        i++;
    }    
}



void AugmentedBetaComplex::adjustIDOfAugmentedBCEdges()
{
    int i = 0;
    AugmentedBCEdge* currEdge = rg_NULL;
    m_augmentedBCEdges.reset4Loop();
    while (m_augmentedBCEdges.setNext4Loop() )
    {
        currEdge = m_augmentedBCEdges.getpEntity();
        currEdge->setID(i);
        i++;
    }

    m_artificialAugmentedBCEdges.reset4Loop();
    while (m_artificialAugmentedBCEdges.setNext4Loop() )
    {
        currEdge = m_artificialAugmentedBCEdges.getpEntity();
        currEdge->setID(i);
        i++;
    }
}



void AugmentedBetaComplex::adjustIDOfAugmentedBCFaces()
{
    int i = 0;
    AugmentedBCFace* currFace = rg_NULL;
    m_augmentedBCFaces.reset4Loop();
    while (m_augmentedBCFaces.setNext4Loop() )
    {
        currFace = m_augmentedBCFaces.getpEntity();
        currFace->setID(i);
        i++;
    }

    m_artificialAugmentedBCFaces.reset4Loop();
    while (m_artificialAugmentedBCFaces.setNext4Loop() )
    {
        currFace = m_artificialAugmentedBCFaces.getpEntity();
        currFace->setID(i);
        i++;
    }
}



void AugmentedBetaComplex::adjustIDOfAugmentedBCCells()
{
    int i = 0;
    AugmentedBCCell* currCell = rg_NULL;
    m_augmentedBCCells.reset4Loop();
    while (m_augmentedBCCells.setNext4Loop() )
    {
        currCell = m_augmentedBCCells.getpEntity();  
        currCell->setID(i);
        i++;
    }

    m_artificialAugmentedBCCells.reset4Loop();
    while (m_artificialAugmentedBCCells.setNext4Loop() )
    {
        currCell = m_artificialAugmentedBCCells.getpEntity();
        currCell->setID(i);
        i++;
    }
}


void AugmentedBetaComplex::adjustIDOfBalls()
{
	int i = 0;
    Ball* currBall = rg_NULL;
    m_balls.reset4Loop();
    while (m_balls.setNext4Loop() )
    {
        currBall = m_balls.getpEntity();  
        currBall->setID(i);
        i++;
    }
}





void AugmentedBetaComplex::adjustIDOfQTForBeta( BetaUniverse* betaUniverse )
{
	adjustIDOfQTForBetaVertices(betaUniverse);
	adjustIDOfQTForBetaEdges(betaUniverse);
	adjustIDOfQTForBetaFaces(betaUniverse);
	adjustIDOfQTForBetaCells(betaUniverse);

	adjustIDOfBallsInQTForBeta(betaUniverse);
}


void AugmentedBetaComplex::adjustIDOfQTForBetaVertices( BetaUniverse* betaUniverse )
{
    int i = 0;
    
    rg_dList< BetaVertex >* vertexList = betaUniverse->getVertexList();
    vertexList->reset4Loop();
    while (vertexList->setNext4Loop() )
    {
        BetaVertex* currVertex = vertexList->getpEntity();
        currVertex->setID(i);
        i++;        
    }
}



void AugmentedBetaComplex::adjustIDOfQTForBetaEdges( BetaUniverse* betaUniverse )
{
    int i = 0;
    
    rg_dList< BetaEdge >* edgeList = betaUniverse->getEdgeList();
    edgeList->reset4Loop();
    while (edgeList->setNext4Loop() )
    {
        BetaEdge* currEdge = edgeList->getpEntity();
        currEdge->setID(i);
        i++;        
    }
}



void AugmentedBetaComplex::adjustIDOfQTForBetaFaces( BetaUniverse* betaUniverse )
{
    int i = 0;
    
    rg_dList< BetaFace >* faceList = betaUniverse->getFaceList();
    faceList->reset4Loop();
    while (faceList->setNext4Loop() )
    {
        BetaFace* currFace = faceList->getpEntity();
        currFace->setID(i);
        i++;        
    }
}



void AugmentedBetaComplex::adjustIDOfQTForBetaCells( BetaUniverse* betaUniverse )
{
    int i = 0;
    
    rg_dList< BetaCell >* cellList = betaUniverse->getCellList();
    cellList->reset4Loop();
    while (cellList->setNext4Loop() )
    {
        BetaCell* currCell = cellList->getpEntity();
        currCell->setID(i);
        i++;        
    }
}


void AugmentedBetaComplex::adjustIDOfBallsInQTForBeta( BetaUniverse* betaUniverse )
{
	int i = 0;
    
    rg_dList< Ball >* ballList = betaUniverse->getBallList();
    ballList->reset4Loop();
    while (ballList->setNext4Loop() )
    {
        Ball* currBall = ballList->getpEntity();
        currBall->setID(i);
        i++;        
    }
}






void AugmentedBetaComplex::duplicateOriginalQuasiTriangulationForBetaToAugmentedBetaComplex( BetaUniverse* betaUniverse )
{
    adjustIDOfQTForBeta( betaUniverse );

    rg_INT numOfBCVertices   = betaUniverse->getVertexList()->getSize();
    rg_INT numOfBCEdges      = betaUniverse->getEdgeList()->getSize();
    rg_INT numOfBCFaces      = betaUniverse->getFaceList()->getSize();
    rg_INT numOfBCCells		 = betaUniverse->getCellList()->getSize();

    AugmentedBCVertex**   ptOfAugVertices	= new AugmentedBCVertex*[numOfBCVertices];
    AugmentedBCEdge**     ptOfAugEdges		= new AugmentedBCEdge*[numOfBCEdges];
    AugmentedBCFace**     ptOfAugFaces		= new AugmentedBCFace*[numOfBCFaces];
    AugmentedBCCell**	  ptOfAugCells		= new AugmentedBCCell*[numOfBCCells];

    rg_INT i;
    for (i = 0; i < numOfBCVertices ; i++)    {
        AugmentedBCVertex currAugmentedBCVertex(i);
        ptOfAugVertices[i] = m_augmentedBCVertices.add( currAugmentedBCVertex );
    }

    for (i = 0; i < numOfBCEdges ; i++)    {
        AugmentedBCEdge currAugmentedBCEdge(i);
        ptOfAugEdges[i] = m_augmentedBCEdges.add( currAugmentedBCEdge );
    }

    for (i = 0; i < numOfBCFaces ; i++)    {
        AugmentedBCFace currAugmentedBCFace(i);
        ptOfAugFaces[i] = m_augmentedBCFaces.add( currAugmentedBCFace );
    }

    for (i = 0; i < numOfBCCells ; i++)    {
        AugmentedBCCell currAugmentedBCCell(i);
        ptOfAugCells[i] = m_augmentedBCCells.add( currAugmentedBCCell );
    }


    duplicateOriginalBetaVerticesToAugmentedBCVertices( betaUniverse->getVertexList(), ptOfAugVertices, ptOfAugEdges, ptOfAugCells );
    duplicateOrigianlBetaEdgesToAugmentedBCEdges(       betaUniverse->getEdgeList(), ptOfAugVertices, ptOfAugEdges, ptOfAugFaces);
    duplicateOriginalBetaCellsToAugmentedBCCells(		betaUniverse->getCellList(), ptOfAugVertices, ptOfAugFaces, ptOfAugCells);
    duplicateOriginalBetaFacesToAugmentedBCFaces(       betaUniverse->getFaceList(), ptOfAugEdges, ptOfAugFaces, ptOfAugCells);  //cell 정보가 있어야 left cell을 extraneous가 아니도록 할 수 있다.
	
    if ( ptOfAugVertices != rg_NULL )
        delete [] ptOfAugVertices;
    if ( ptOfAugEdges != rg_NULL )
        delete [] ptOfAugEdges;
    if ( ptOfAugFaces != rg_NULL )
        delete [] ptOfAugFaces;  
    if ( ptOfAugCells != rg_NULL )
        delete [] ptOfAugCells;


    //betaUniverse->rearrangeSimplexIDs();
}


void AugmentedBetaComplex::duplicateOriginalBetaVerticesToAugmentedBCVertices(rg_dList<BetaVertex>* betaVertices,
																AugmentedBCVertex** ptOfAugVertices,
																AugmentedBCEdge**   ptOfAugEdges,
																AugmentedBCCell**   ptOfAugCells)
{
    BetaVertex*  currBetaVertex   = rg_NULL;
    rg_INT       currBetaVertexID = 0;
    AugmentedBCVertex* currAugmentedBCVertex = rg_NULL;
    betaVertices->reset4Loop();
    while( betaVertices->setNext4Loop() )
    {
        currBetaVertex   = betaVertices->getpEntity();
        currBetaVertexID = currBetaVertex->getID();

        currAugmentedBCVertex = ptOfAugVertices[currBetaVertexID];

        currAugmentedBCVertex->setFirstEdge( ptOfAugEdges[ currBetaVertex->getFirstEdge()->getID() ]);
        if ( currBetaVertex->getFirstCell() == rg_NULL ) {
            currAugmentedBCVertex->setFirstCell( rg_NULL );
        }
        else {
            currAugmentedBCVertex->setFirstCell( ptOfAugCells[ currBetaVertex->getFirstCell()->getID() ]);
        }
        

		if ( !currBetaVertex->isVirtual() )
		{
			currAugmentedBCVertex->setBetaSpan( currBetaVertex->getBetaSpan() );			

			Ball* ptOfBall = m_balls.add( *(currBetaVertex->getBallProperty()) );
			currAugmentedBCVertex->setBallProperty( ptOfBall ); 
		}        
		
		currAugmentedBCVertex->setBoundingState( currBetaVertex->getBoundingState( m_betaValue ) );
		
        currAugmentedBCVertex->setOriginalBetaVertex( currBetaVertex );
        currAugmentedBCVertex->setIsArtificial( rg_FALSE );
		
		currAugmentedBCVertex->isVisited( rg_FALSE );
    }    
}



void AugmentedBetaComplex::duplicateOrigianlBetaEdgesToAugmentedBCEdges( rg_dList<BetaEdge>* betaEdges, 
															   AugmentedBCVertex** ptOfAugVertices, 
															   AugmentedBCEdge**   ptOfAugEdges, 
															   AugmentedBCFace**   ptOfAugFaces)
{
    BetaEdge*  currBetaEdge   = rg_NULL;
    rg_INT     currBetaEdgeID = 0;
    AugmentedBCEdge* currAugmentedBCEdge = rg_NULL;
    betaEdges->reset4Loop();
    while( betaEdges->setNext4Loop() )
    {
        currBetaEdge = betaEdges->getpEntity();
        currBetaEdgeID = currBetaEdge->getID();


        currAugmentedBCEdge = ptOfAugEdges[ currBetaEdgeID ];

        currAugmentedBCEdge->setStartVertex( ptOfAugVertices[ currBetaEdge->getStartVertex()->getID() ] );
        currAugmentedBCEdge->setEndVertex(   ptOfAugVertices[ currBetaEdge->getEndVertex()->getID() ] );
        currAugmentedBCEdge->setFirstFace(   ptOfAugFaces[ currBetaEdge->getFirstFace()->getID() ] );
		
		rg_dList<BetaFace*>* smallWorlds = currBetaEdge->getSmallWorlds();
		BetaFace*			 currSmallWorldFace = rg_NULL;
		smallWorlds->reset4Loop();
		while ( smallWorlds->setNext4Loop() )
		{
			currSmallWorldFace = smallWorlds->getEntity();
			currAugmentedBCEdge->getSmallWorlds()->add( ptOfAugFaces[ currSmallWorldFace->getID() ] );
		}

		if ( !currAugmentedBCEdge->isVirtual() )		
		{
			currAugmentedBCEdge->setBetaSpan( currBetaEdge->getBetaSpan() );
		}
		
		currAugmentedBCEdge->setBoundingState( currBetaEdge->getBoundingState( m_betaValue ) );
		currAugmentedBCEdge->setIsArtificial( rg_FALSE );

		currAugmentedBCEdge->isVisited( rg_FALSE );
    }
}




void AugmentedBetaComplex::duplicateOriginalBetaFacesToAugmentedBCFaces( rg_dList<BetaFace>* betaFaces, 
                                                     AugmentedBCEdge**   ptOfAugEdges, 
                                                     AugmentedBCFace**   ptOfAugFaces, 
                                                     AugmentedBCCell**   ptOfAugCells)
{    
    BetaFace* currBetaFace   = rg_NULL;
    rg_INT    currBetaFaceID = 0;

    AugmentedBCFace* currAugmentedBCFace = rg_NULL;
    betaFaces->reset4Loop();
    while( betaFaces->setNext4Loop() )
    {
        currBetaFace   = betaFaces->getpEntity();
        currBetaFaceID = currBetaFace->getID();    
        
        currAugmentedBCFace = ptOfAugFaces[ currBetaFaceID ];

        if ( currBetaFace->getRightCell() == rg_NULL ) { //is currBetaFace an isolated face?
            currAugmentedBCFace->setRightCell( rg_NULL );
        }
        else {
            currAugmentedBCFace->setRightCell( ptOfAugCells[ currBetaFace->getRightCell()->getID() ] );
        }
        
        if ( currBetaFace->getLeftCell() == rg_NULL ) {
            currAugmentedBCFace->setLeftCell( rg_NULL );
        }
        else {
            currAugmentedBCFace->setLeftCell(  ptOfAugCells[ currBetaFace->getLeftCell()->getID() ] );
        }

        

		for (rg_INT i = 0; i < EIWDS_NUM_EDGE_ON_FACE; i++)        {
            currAugmentedBCFace->setEdge(i, currBetaFace->getEdgeOrientation(i), 
                                  ptOfAugEdges[ currBetaFace->getEdge(i)->getID() ] );
        }

		if ( !currAugmentedBCFace->isVirtual() )		
		{
			currAugmentedBCFace->setBetaSpan( currBetaFace->getBetaSpan() );
		}
        
        currAugmentedBCFace->setBoundingState( currBetaFace->getBoundingState( m_betaValue ) );

		
		currAugmentedBCFace->setIsArtificial( rg_FALSE );

		currAugmentedBCFace->isVisited( rg_FALSE );

        //reverse orientations if leftCell is "Virtual" or "Extraneous" <-regular face의 경우 항상 right cell이 rg_NULL이도록 setting
        AugmentedBCCell* leftCell = (AugmentedBCCell*)currAugmentedBCFace->getLeftCell();
        
        if ( leftCell != rg_NULL ) {
            rg_INT leftCellBoundingState= leftCell->getBoundingState();
		    if (    leftCellBoundingState == EXTRANEOUS_SIMPLEX 
                 || leftCellBoundingState == FALSE_EXTERIOR_SIMPLEX ) {
                currAugmentedBCFace->reverse();            
            }     
        }
	}
}


void AugmentedBetaComplex::duplicateOriginalBetaCellsToAugmentedBCCells( rg_dList<BetaCell>* betaCells, 
													 AugmentedBCVertex** ptOfAugVertices, 
                                                     AugmentedBCFace**   ptOfAugFaces, 
                                                     AugmentedBCCell**   ptOfAugCells)
{
    BetaCell*  currBetaCell   = rg_NULL;
    rg_INT     currBetaCellID = 0;
    AugmentedBCCell* currAugmentedBCCell = rg_NULL;
    betaCells->reset4Loop();
    while ( betaCells->setNext4Loop() )
    {
        currBetaCell   = betaCells->getpEntity();
        currBetaCellID = currBetaCell->getID();

        currAugmentedBCCell = ptOfAugCells[ currBetaCellID ];

		rg_INT i = 0;
        for (i = 0; i < EIWDS_NUM_FACE_ON_CELL; i++)  {
            currAugmentedBCCell->setFace(i, ptOfAugFaces[ currBetaCell->getFace(i)->getID() ] );
        }

		for (i = 0; i < EIWDS_NUM_VERTEX_ON_CELL; i++)  {
            currAugmentedBCCell->setVertex(i, ptOfAugVertices[ currBetaCell->getVertex(i)->getID() ] );
        }

		if ( !currAugmentedBCCell->isVirtual() )		
		{
			currAugmentedBCCell->setBetaSpan( currBetaCell->getBetaSpan() );
			currAugmentedBCCell->setMinTangentSphere( currBetaCell->getMinTangentSphere() );
		}

        currAugmentedBCCell->setBoundingState( currBetaCell->getBoundingState( m_betaValue ) );
		currAugmentedBCCell->setIsArtificial( rg_FALSE );

		currAugmentedBCCell->isVisited( rg_FALSE );
    }
}


void AugmentedBetaComplex::indexExteriorRegionOfEdgesAndFaces()
{
    rg_INT indexOfExteriorRegion = 0;	

    
	m_augmentedBCCells.reset4Loop();
	while ( m_augmentedBCCells.setNext4Loop() )
	{
		AugmentedBCCell* currCell = m_augmentedBCCells.getpEntity();

        rg_INT currCellBoundingState = currCell->getBoundingState();
		if (    currCellBoundingState != EXTRANEOUS_SIMPLEX
             && currCellBoundingState != FALSE_EXTERIOR_SIMPLEX )
			continue;

		if ( currCell->isVisited() )
			continue;

		rg_dList<AugmentedBCCell*> cellStack;
		cellStack.pushBack( currCell );

		while ( cellStack.getSize() > 0 )
		{
			AugmentedBCCell* popedCell = cellStack.popBack();

            if ( popedCell->isVisited() ) {
                continue;
            }

			popedCell->isVisited( rg_TRUE );

            BetaFace** boundingFace = popedCell->getFaces();

                                                      
			rg_INT i=0;
            for ( i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
                // popedCell의 boundingEdge 중 singular Edge가 있다면
                // indexOfExterior Region을 입력.
                BetaEdge** boundingEdge = boundingFace[i]->getEdges();

                for ( rg_INT j=0; j<EIWDS_NUM_EDGE_ON_FACE; j++ ) {
                    if ( boundingEdge[j]->isVisited() ) {
                        continue;
                    }

                    boundingEdge[j]->isVisited( rg_TRUE );

                    if ( ((AugmentedBCEdge*) boundingEdge[j])->getBoundingState() == SINGULAR_SIMPLEX ) {
                        ((AugmentedBCEdge*) boundingEdge[j])->setIndexOfExteriorRegion( indexOfExteriorRegion );
                    }
                }


                // popedCell의 boundingFace 중
                // regular  face : indexOfExteriorToRight Region = indexOfExteriorRegion
                // singular face : popedCeel이 singularFace의 leftCell이면 indexOfExteriorToLeft에 입력                
                //                                           rightCell이면 indexExteriorToRight에 입력

                /* singular face는 2번 marking 되어야 하므로 아래 조건 없애야 한다.
                if ( boundingFace[i]->isVisited() ) {
                    continue;
                }

                boundingFace[i]->isVisited( rg_TRUE );
                */

                rg_INT boundingStateOfCurrBoundingFace = ((AugmentedBCFace*) boundingFace[i])->getBoundingState();
                if (    boundingStateOfCurrBoundingFace == SINGULAR_SIMPLEX
                     || boundingStateOfCurrBoundingFace == REGULAR_SIMPLEX ) {
                    if ( boundingFace[i]->getRightCell() == popedCell ) {
                        ((AugmentedBCFace*) boundingFace[i])->setIndexOfExteriorRegionToRightCell( indexOfExteriorRegion );
                    }

                    if ( boundingFace[i]->getLeftCell() == popedCell ) {
                        ((AugmentedBCFace*) boundingFace[i])->setIndexOfExteriorRegionToLeftCell( indexOfExteriorRegion );
                    }
                    
                }
            }
            
            // find nextCells and push-back it to stack
			for ( i=0; i<EIWDS_NUM_FACE_ON_CELL; i++) {
                rg_INT boundingStateOfCurrBoundingFace = ((AugmentedBCFace*) boundingFace[i])->getBoundingState();
				if (    boundingStateOfCurrBoundingFace != EXTRANEOUS_SIMPLEX
                     && boundingStateOfCurrBoundingFace != FALSE_EXTERIOR_SIMPLEX )
					continue;

				AugmentedBCCell* nextCell = rg_NULL;
                if ( boundingFace[i]->getLeftCell() == popedCell ) {
					nextCell = (AugmentedBCCell*) boundingFace[i]->getRightCell();
                }
                else {
					nextCell = (AugmentedBCCell*) boundingFace[i]->getLeftCell();
                }

                rg_INT nextCellBoundingState = nextCell->getBoundingState();
				if (    nextCellBoundingState == EXTRANEOUS_SIMPLEX
                     || nextCellBoundingState == FALSE_EXTERIOR_SIMPLEX ) {
                    if ( nextCell->isVisited() == rg_FALSE ) {
					    cellStack.pushBack( nextCell );
                    }
				}
			}
		}

		indexOfExteriorRegion++;
	}



	m_augmentedBCEdges.reset4Loop();
    while ( m_augmentedBCEdges.setNext4Loop() ) {
		m_augmentedBCEdges.getpEntity()->isVisited( rg_FALSE );
    }

    /*
    m_augmentedBCFaces.reset4Loop();
    while ( m_augmentedBCFaces.setNext4Loop() ) {
		m_augmentedBCFaces.getpEntity()->isVisited( rg_FALSE );
    }
    */

    m_augmentedBCCells.reset4Loop();
    while ( m_augmentedBCCells.setNext4Loop() ) {
		m_augmentedBCCells.getpEntity()->isVisited( rg_FALSE );
    }
	
}


void AugmentedBetaComplex::searchNonManifoldConfigurationsForAllVertices( rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMVertices, 
                                                                          rg_dList< NonManifoldCorner >& nonManifoldCornerList)  const
{
	AugmentedBCVertex* currVertex;
    m_augmentedBCVertices.reset4Loop();
    
    while ( m_augmentedBCVertices.setNext4Loop() )
    {
        currVertex = m_augmentedBCVertices.getpEntity();
        if ( currVertex->getBoundingState() == REGULAR_SIMPLEX )  //Extraneous, singular, or interior vertices are ignored.
            generateNonManifoldVertexConfiguration(currVertex, ptToNMCornerOfNMVertices, nonManifoldCornerList);
    }
}




void AugmentedBetaComplex::searchNonManifoldConfigurationsForAllEdges( rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMEdges, 
                                                                       rg_dList< NonManifoldCorner >& nonManifoldCornerList)  
{
    AugmentedBCEdge* currEdge;
    m_augmentedBCEdges.reset4Loop();
    while ( m_augmentedBCEdges.setNext4Loop() )
    {
        currEdge = m_augmentedBCEdges.getpEntity();

        if ( currEdge->getBoundingState() == REGULAR_SIMPLEX )  ////Extraneous, singular, or interior edges are ignored.
            generateNonManifoldEdgeConfiguration( currEdge, ptToNMCornerOfNMEdges, nonManifoldCornerList );
    }
}





void AugmentedBetaComplex::generateNonManifoldVertexConfiguration(AugmentedBCVertex*             givenVertex, 
                                                                  rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMVertices, 
                                                                  rg_dList< NonManifoldCorner >& nonManifoldCornerList) const
{
    rg_dList< BetaEdge* > singularEdgeList;
    findSingularEdgesIncidentToVertex( givenVertex, singularEdgeList );
    //givenVertex->searchEdgesInBetaComplex( m_betaValue, SINGULAR_SIMPLEX, singularEdgeList);

	rg_dList< BetaFace* > singularFaceList;
    findSingularFacesIncidentToVertex( givenVertex, singularFaceList );
	//givenVertex->searchFacesInBetaComplex( m_betaValue, SINGULAR_SIMPLEX, singularFaceList);

	rg_dList< rg_dList< BetaCell* > > interiorCellGroupList;
	//givenVertex->searchGroupsOfFaceConnectedInteriorCellsInBetaComplex( m_betaValue, interiorCellGroupList );
    searchGroupsOfFaceConnectedInteriorCellsInBetaComplex( givenVertex, interiorCellGroupList );
    // bounding state를 m_betaValue를 통해 확인하는 것이 아니라 해당 simplex의 getBoundingState()라는 함수로 확인하는 함수를
    // this에서 가지고 있도록 하였다.

	rg_dList< rg_dList< BetaCell* > > extraneousCellGroupList;
	//givenVertex->searchGroupsOfFaceConnectedExtraneousCellsInBetaComplex( m_betaValue, extraneousCellGroupList );
    searchGroupsOfFaceConnectedExtraneousCellsInBetaComplex( givenVertex, extraneousCellGroupList );


    rg_INT numOfSingularEdges        = singularEdgeList.getSize();
    rg_INT numOfSingularFaces        = singularFaceList.getSize();
	rg_INT numOfInteriorCellGroups   = interiorCellGroupList.getSize();
	rg_INT numOfVolumetricComponents = numOfSingularEdges + numOfSingularFaces + numOfInteriorCellGroups;

	rg_INT numOfExtraneousCellGroups = extraneousCellGroupList.getSize();


    //  ycho 2011-04-23
    //r_numExtraneousCellGroupAtVertex[givenVertex->getID()] = numOfExtraneousCellGroups;

    if ( numOfVolumetricComponents == 1 
        && numOfExtraneousCellGroups == 1)   
        return;                                 //currVertex is a manifold vertex

	
	
	
	generateNMVertexCornersInEachExtraneousCellGroup( givenVertex, extraneousCellGroupList, ptToNMCornerOfNMVertices, nonManifoldCornerList );

    //TT_D_V configuration
    if ( numOfExtraneousCellGroups > 1 )
    {
        NonManifoldCorner tempNMConf( givenVertex );
        NonManifoldCorner* tempPtNMConf = nonManifoldCornerList.addTail( tempNMConf );
		tempPtNMConf->getpFirstSimplex()->setTypeOfSimplex( numOfExtraneousCellGroups - 1 ); //for counting the number of TT^D_V configuration.
    }
}







void AugmentedBetaComplex::findSingularEdgesIncidentToVertex( AugmentedBCVertex*     givenVertex, 
                                                              rg_dList< BetaEdge* >& singularEdgeList ) const
{
    rg_dList< BetaEdge* > incidentEdges;
    givenVertex->searchFiniteEdgesInWholeWorld( incidentEdges );

    incidentEdges.reset4Loop();
    while ( incidentEdges.setNext4Loop() )   {
        AugmentedBCEdge* currEdge = (AugmentedBCEdge*) incidentEdges.getEntity();

        if ( currEdge->getBoundingState() == SINGULAR_SIMPLEX ) {
            singularEdgeList.add( (BetaEdge*) currEdge );
        }
    }
}




void AugmentedBetaComplex::findSingularFacesIncidentToVertex( AugmentedBCVertex*     givenVertex, 
                                                              rg_dList< BetaFace* >& singularFaceList  ) const
{
    rg_dList< BetaFace* > incidentFaces;
    givenVertex->searchFiniteFacesInWholeWorld( incidentFaces );

    incidentFaces.reset4Loop();
    while ( incidentFaces.setNext4Loop() )   {
        AugmentedBCFace* currFace = (AugmentedBCFace*) incidentFaces.getEntity();

        if ( currFace->getBoundingState() == SINGULAR_SIMPLEX ) {
            singularFaceList.add( (BetaFace*) currFace );
        }
    }
}






void AugmentedBetaComplex::generateNMVertexCornersInEachExtraneousCellGroup( AugmentedBCVertex*               givenVertex, 
                                                                             rg_dList< rg_dList<BetaCell*> >& extraneousCellGroupList, 
																			 rg_dList<NonManifoldCorner*>*    ptToNMCornerOfNMVertices,
																			 rg_dList< NonManifoldCorner >&   nonManifoldCornerList ) const
{
	extraneousCellGroupList.reset4Loop();
	while( extraneousCellGroupList.setNext4Loop() )
	{
		rg_dList< BetaCell* >* currExtraneousCellGroup = extraneousCellGroupList.getpEntity();

		rg_dList< BetaEdge* > singularEdgesInCurrExtraneousCellGroup;
		findSingularEdgesIncidentToVertexInCellGroup(givenVertex, currExtraneousCellGroup, singularEdgesInCurrExtraneousCellGroup);

		rg_dList< BetaFace* > betaShapeFacesInCurrExtraneousCellGroup;
		findBetaShapeFacesIncidentToVertexInCellGroup(givenVertex, currExtraneousCellGroup, betaShapeFacesInCurrExtraneousCellGroup);

		//find edge-connected face group		
		rg_dList< rg_dList < BetaFace* > > edgeConnectedFaceGroupList;
		separateFacesIntoEdgeConnectedComponent( givenVertex, betaShapeFacesInCurrExtraneousCellGroup, edgeConnectedFaceGroupList );

		
		///////////////////////////////////////////////
		//2. Generate nm-corners
		//2.1 EE_V configuration
		generateNMVertexCornersWithSingularEdges( givenVertex, singularEdgesInCurrExtraneousCellGroup, ptToNMCornerOfNMVertices, nonManifoldCornerList );
		
		//2.2 EF_V or ET_V configuration
		if ( singularEdgesInCurrExtraneousCellGroup.getSize() > 0 
		  && edgeConnectedFaceGroupList.getSize() > 0 )
		{
			PointerToSimplex tempFirstSimplex( (AugmentedBCEdge*)singularEdgesInCurrExtraneousCellGroup.getLastEntity() );

			PointerToSimplex tempSecondSimplex;			

			BetaFace* regularFace = findRegularFaceInGroupOfFaceList( edgeConnectedFaceGroupList );
			if ( regularFace == rg_NULL ) //EF_V
			{
				tempSecondSimplex.setPointerToSimplex( FACE_SIMPLEX, (AugmentedBCFace*)edgeConnectedFaceGroupList.getFirstpEntity()->getFirstEntity());		
			}
			else //ET_V
			{
				tempSecondSimplex.setPointerToSimplex( CELL_SIMPLEX, (AugmentedBCFace*)regularFace  );	//cell을 넣어주면, resolve할 때 다른 extraneous영역에 있는 regular face를 가져올 수 있으므로 face를 직접 넣어준다.				
			}
    
            NonManifoldCorner  tempNMConf( givenVertex, tempFirstSimplex, tempSecondSimplex );
            NonManifoldCorner* tempPtrNMConf = nonManifoldCornerList.pushFront( tempNMConf );
        
            ptToNMCornerOfNMVertices[givenVertex->getID()].add( tempPtrNMConf );
		}

		//2.3 FF_V, FT_V, or TT_V configuration
		generateNMVertexCornersWithEdgeConnectedFaceGroup(  givenVertex, edgeConnectedFaceGroupList, ptToNMCornerOfNMVertices, nonManifoldCornerList );
	}
}




BetaFace* AugmentedBetaComplex::findRegularFaceInGroupOfFaceList( rg_dList< rg_dList < BetaFace* > >& faceGroupList ) const
{
	faceGroupList.reset4Loop();
	while ( faceGroupList.setNext4Loop() )
	{
		rg_dList < BetaFace* >* currFaceGroup = faceGroupList.getpEntity();

		currFaceGroup->reset4Loop();
		while ( currFaceGroup->setNext4Loop() )
		{
			BetaFace* currFace = currFaceGroup->getEntity();

			if ( ((AugmentedBCFace*)currFace)->getBoundingState() == REGULAR_SIMPLEX )
				return currFace;
		}
	}
	return rg_NULL;
}



BetaFace* AugmentedBetaComplex::findRegularFaceInFaceList(rg_dList < BetaFace* >* faceList) const
{
	faceList->reset4Loop();
	while ( faceList->setNext4Loop() )
	{
		BetaFace* currFace = faceList->getEntity();

		if ( ((AugmentedBCFace*)currFace)->getBoundingState() == REGULAR_SIMPLEX )
			return currFace;	
	}
	return rg_NULL;
}



//2.3 FF_V, FT_V, or TT_V configuration
void AugmentedBetaComplex::generateNMVertexCornersWithEdgeConnectedFaceGroup(  AugmentedBCVertex*               givenVertex, 
                                                                               rg_dList< rg_dList<BetaFace*> >& edgeConnectedFaceList, 
                                                                              rg_dList<NonManifoldCorner*>*     ptToNMCornerOfNMVertices,
                                                                              rg_dList< NonManifoldCorner >&    nonManifoldCornerList ) const
{
	BetaFace* currGroupFace = rg_NULL;
    BetaFace* nextGroupFace = rg_NULL;
    
    rg_INT numOfGroups = edgeConnectedFaceList.getSize();
    rg_INT i=0;

    edgeConnectedFaceList.reset4Loop();
    while( edgeConnectedFaceList.setNext4Loop() && i < numOfGroups - 1)
    {
		PointerToSimplex tempNMSimplex(VERTEX_SIMPLEX, (void*)givenVertex );
        PointerToSimplex tempSimplexOfCurrFace;
		PointerToSimplex tempSimplexOfNextFace;
		
		//singular face는 방향이 없으므로 regular face를 뽑는다.(regular face가 없을 경우 singular face를 뽑는다.)
        currGroupFace = findRegularFaceInFaceList( edgeConnectedFaceList.getpEntity() );
		if ( currGroupFace == rg_NULL )	//currSimplex = Face
		{
			tempSimplexOfCurrFace.setPointerToSimplex( FACE_SIMPLEX, (AugmentedBCFace*)edgeConnectedFaceList.getpEntity()->getFirstEntity());
		}
		else	//currSimplex = Cell
		{
			tempSimplexOfCurrFace.setPointerToSimplex( CELL_SIMPLEX, (AugmentedBCFace*)currGroupFace );
		}
			
        nextGroupFace = findRegularFaceInFaceList( edgeConnectedFaceList.getpNextEntity() );
		if ( nextGroupFace == rg_NULL ) //nextSimplex = Face
		{
			tempSimplexOfNextFace.setPointerToSimplex( FACE_SIMPLEX, (AugmentedBCFace*)edgeConnectedFaceList.getpNextEntity()->getFirstEntity() );
		}
		else	//nextSimplex = Cell
		{
			tempSimplexOfNextFace.setPointerToSimplex( CELL_SIMPLEX, (AugmentedBCFace*)nextGroupFace );
		}	

		NonManifoldCorner  tempNMConf;
		tempNMConf.setNMSimplex( tempNMSimplex );
		
		tempNMConf.setFirstSimplex( tempSimplexOfCurrFace );
		tempNMConf.setSecondSimplex( tempSimplexOfNextFace) ;			
		

        NonManifoldCorner* tempPtNMConf = nonManifoldCornerList.pushFront( tempNMConf );

        ptToNMCornerOfNMVertices[givenVertex->getID()].add( tempPtNMConf );

        i++;
    }

}


//EE_V
void AugmentedBetaComplex::generateNMVertexCornersWithSingularEdges( AugmentedBCVertex*             givenVertex, 
                                                                     rg_dList<BetaEdge*>&           singularEdgeList, 
                                                                     rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMVertices,
                                                                     rg_dList< NonManifoldCorner >& nonManifoldCornerList ) const
{
    AugmentedBCEdge* currSingularEdge = rg_NULL;
    AugmentedBCEdge* nextSingularEdge = rg_NULL;
    
    rg_INT numOfSingularEdges = singularEdgeList.getSize();
    rg_INT i=0;

    singularEdgeList.reset4Loop();
    while( singularEdgeList.setNext4Loop() && i < numOfSingularEdges - 1)
    {
        currSingularEdge = (AugmentedBCEdge*)singularEdgeList.getEntity();
        nextSingularEdge = (AugmentedBCEdge*)singularEdgeList.getNextEntity();

        NonManifoldCorner  tempNMConf( givenVertex, currSingularEdge, nextSingularEdge );
        NonManifoldCorner* tempPtNMConf = nonManifoldCornerList.pushFront( tempNMConf );

        ptToNMCornerOfNMVertices[givenVertex->getID()].add( tempPtNMConf );

        i++;
    }
}



void AugmentedBetaComplex::separateFacesIntoEdgeConnectedComponent( AugmentedBCVertex*                  givenVertex, 
                                                                    rg_dList < BetaFace* >&             facesToBeSeparatedByEdgeConnection, 
                                                                    rg_dList< rg_dList < BetaFace* > >& edgeConnectedFaceGroupList ) const
{
	while ( facesToBeSeparatedByEdgeConnection.getSize() > 0 )
	{
		rg_dNode< BetaFace* >* firstNode = facesToBeSeparatedByEdgeConnection.getFirstpNode();
		BetaFace* firstFace				 = firstNode->getEntity();

		rg_dList< rg_dNode< BetaFace* >* > nodeToBeKilled;

		rg_dList<BetaFace*>* currEdgeConnectedFaceGroup = 
								edgeConnectedFaceGroupList.add( rg_dList<BetaFace*>() );

		currEdgeConnectedFaceGroup->add( firstFace );
		nodeToBeKilled.add( firstNode );

		if ( facesToBeSeparatedByEdgeConnection.getSize() == 1 )
		{
			facesToBeSeparatedByEdgeConnection.kill( firstNode );
			break;
		}

		rg_dNode< BetaFace* >* currNode = facesToBeSeparatedByEdgeConnection.getSecondpNode();

		do 
		{
			BetaFace* currFace = currNode->getEntity();

			if ( isTheseFacesConnectedWithEdgesWhichBoundedByGivenVertex( givenVertex, firstFace, currFace ) )
			{
				currEdgeConnectedFaceGroup->add( currFace );
				nodeToBeKilled.add( currNode );
			}

			currNode = currNode->getNext();
		} while(currNode != firstNode);

		nodeToBeKilled.reset4Loop();
		while ( nodeToBeKilled.setNext4Loop() )
		{
			rg_dNode< BetaFace* >* currNode = nodeToBeKilled.getEntity();
			facesToBeSeparatedByEdgeConnection.kill( currNode );
		}
	}
}





rg_BOOL AugmentedBetaComplex::isTheseFacesConnectedWithEdgesWhichBoundedByGivenVertex(AugmentedBCVertex* givenVertex, BetaFace* sourceFace, BetaFace* targetFace) const
{
    rg_dList< BetaFace* > VisitedFacesList;
    
    sourceFace->isVisited( rg_TRUE );
    VisitedFacesList.add( sourceFace );
    
    BetaEdge* edgesBoundedByVertex[2];
    sourceFace->findEdgeIncidentToVertex( givenVertex, edgesBoundedByVertex);

    rg_dList< BetaEdge* > traverseEdgesStack;
    traverseEdgesStack.add( edgesBoundedByVertex[0] );
    traverseEdgesStack.add( edgesBoundedByVertex[1] );

    BetaEdge* currTraverseEdge = rg_NULL;
    while( traverseEdgesStack.getSize() > 0 )
    {
        currTraverseEdge = traverseEdgesStack.popBack();

        rg_dList< BetaFace* > incidentFacesInBetaComplex;
        //currTraverseEdge->searchFacesInBetaComplex( m_betaValue, incidentFacesInBetaComplex ); //searchFacesInBetaShape을 이용할 경우, TT_D_V인 case를 TT_V로 잘 못 detect할 수 있다.
        searchFacesIncidentToEdgeInBetaComplex( currTraverseEdge, incidentFacesInBetaComplex );

        BetaFace* currFace = rg_NULL;
        incidentFacesInBetaComplex.reset4Loop();
        while( incidentFacesInBetaComplex.setNext4Loop() )
        {
            currFace = incidentFacesInBetaComplex.getEntity();

            if ( currFace->isVisited() )
                continue;

            currFace->isVisited( rg_TRUE );
            VisitedFacesList.add( currFace );
            
            if ( currFace == targetFace )
                break;

            BetaEdge* primitivesOfNextEdge[2];
            currFace->findEdgeIncidentToVertex( givenVertex, primitivesOfNextEdge);

            if ( primitivesOfNextEdge[0] != currTraverseEdge )
                traverseEdgesStack.add( primitivesOfNextEdge[0] );
            else if ( primitivesOfNextEdge[1] != currTraverseEdge )
                traverseEdgesStack.add( primitivesOfNextEdge[1] );
        }

        if ( targetFace->isVisited() == rg_TRUE )
            break;
    }

    rg_BOOL isTheseFacesConnectedWithEdges;
    if ( targetFace->isVisited() )
        isTheseFacesConnectedWithEdges = rg_TRUE;
    else //( targetFace->isVisited() == rg_FALSE )
        isTheseFacesConnectedWithEdges = rg_FALSE;
    
    VisitedFacesList.reset4Loop();
    while( VisitedFacesList.setNext4Loop() )    {
        VisitedFacesList.getEntity()->isVisited( rg_FALSE );
    }

    return isTheseFacesConnectedWithEdges;
}




void AugmentedBetaComplex::findSingularEdgesIncidentToVertexInCellGroup(AugmentedBCVertex*     givenVertex, 
                                                                        rg_dList< BetaCell* >* givenCellGroup, 
                                                                        rg_dList< BetaEdge* >& singularEdgesIncidentToGivenVertex) const
{
	givenCellGroup->reset4Loop();
	while ( givenCellGroup->setNext4Loop() )
	{
		BetaCell* currCell = givenCellGroup->getEntity();
		
		rg_dList< BetaEdge* > boundingEdgesOfCurrCell;
		currCell->searchEdgesInWholeWorld( boundingEdgesOfCurrCell );

		boundingEdgesOfCurrCell.reset4Loop();
		while ( boundingEdgesOfCurrCell.setNext4Loop() )
		{
			BetaEdge* currEdge = boundingEdgesOfCurrCell.getEntity();

			if ( ((AugmentedBCEdge*)currEdge)->getBoundingState() == SINGULAR_SIMPLEX
				&& currEdge->isThere( givenVertex ) )
				singularEdgesIncidentToGivenVertex.addWithoutSame( currEdge );
		}
	}
}


void AugmentedBetaComplex::findBetaShapeFacesIncidentToVertexInCellGroup(AugmentedBCVertex*     givenVertex, 
                                                                         rg_dList< BetaCell* >* givenCellGroup, 
																		 rg_dList< BetaFace* >& betaShapeFacesIncidentToGivenVertex) const
{
	givenCellGroup->reset4Loop();
	while ( givenCellGroup->setNext4Loop() )
	{
		BetaCell* currCell = givenCellGroup->getEntity();
		
		rg_dList< BetaFace* > boundingFacesOfCurrCell;
		currCell->searchFacesInWholeWorld( boundingFacesOfCurrCell );

		boundingFacesOfCurrCell.reset4Loop();
		while ( boundingFacesOfCurrCell.setNext4Loop() )
		{
			BetaFace* currFace = boundingFacesOfCurrCell.getEntity();

			if ( ((AugmentedBCFace*)currFace)->getBoundingState() == SINGULAR_SIMPLEX
				|| ((AugmentedBCFace*)currFace)->getBoundingState() == REGULAR_SIMPLEX )
			{
				if ( currFace->isThere( givenVertex ) )
				betaShapeFacesIncidentToGivenVertex.addWithoutSame( currFace );
			}
		}
	}
}



















void AugmentedBetaComplex::generateNonManifoldEdgeConfiguration(AugmentedBCEdge*               givenEdge, 
                                                                rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMEdges,
                                                                rg_dList< NonManifoldCorner >& nonManifoldCornerList)
{
	rg_dList< BetaFace* > facesIncidentToEdgeInBetaShape;
	//givenEdge->searchFacesInBetaShape( m_betaValue, facesIncidentToEdgeInBetaShape );
    searchFacesIncidentToEdgeInBetaShape( givenEdge, facesIncidentToEdgeInBetaShape );
    if ( givenEdge->isGateToSmallWorlds() && facesIncidentToEdgeInBetaShape.getSize() > 2 ) {
        sortFacesCCWOfEdgeByGeometricComparison(givenEdge, facesIncidentToEdgeInBetaShape );
    }
    
    
    AugmentedBCFace* currFace = rg_NULL;
    rg_INT numOfSingularFaces = 0;
    rg_INT numOfRegularFaces  = 0;

    rg_FLAG isThereEllipticQTFace = rg_FALSE;
	facesIncidentToEdgeInBetaShape.reset4Loop();
	while ( facesIncidentToEdgeInBetaShape.setNext4Loop() )
	{
		currFace = (AugmentedBCFace*) facesIncidentToEdgeInBetaShape.getEntity();

        if ( currFace->getBoundingState() == SINGULAR_SIMPLEX ) {
			numOfSingularFaces++;
            
            if ( currFace->isSingularEllipticFace(m_betaValue) ) {
                isThereEllipticQTFace = rg_TRUE;
            }
        }
        else if ( currFace->getBoundingState() == REGULAR_SIMPLEX ) {
			numOfRegularFaces++;
        }
	}

    if ( numOfSingularFaces == 0 && numOfRegularFaces == 2 )
        return;

    if ( numOfSingularFaces <= 1 && numOfRegularFaces == 0)
        return;    
    
    if ( isThereEllipticQTFace ) {
        sortFacesCCWOfEdgeByGeometricComparison(givenEdge, facesIncidentToEdgeInBetaShape );
    }
    
    //FF_E, FT_E, TF_E, TT_E
    generateNonManifoldEdgeConfigurationByOrderOfFacesWRTGivenEdge(givenEdge, facesIncidentToEdgeInBetaShape,
                                                    ptToNMCornerOfNMEdges,
                                                    nonManifoldCornerList);
    
}


//FF_E, FT_E, TF_E, TT_E
void AugmentedBetaComplex::generateNonManifoldEdgeConfigurationByOrderOfFacesWRTGivenEdge(
                                                    AugmentedBCEdge*                givenEdge, 
                                                    rg_dList< BetaFace* >&          facesIncidentToEdgeInBetaShape,
                                                    rg_dList<NonManifoldCorner*>*   ptToNMCornerOfNMEdges,
                                                    rg_dList< NonManifoldCorner >&  nonManifoldCornerList) 
{
    if ( givenEdge->getID() == 4443 ) {
        int stop = 1;
    }
    rg_dNode< BetaFace* >* startNode = facesIncidentToEdgeInBetaShape.getFirstpNode();
    rg_dNode< BetaFace* >* currNode  = startNode;
    rg_dNode< BetaFace* >* initialFaceNode = rg_NULL; //REGULAR_SIMPLEX이고 edge방향이 true인 face를 initial로 하면 안된다.
                                                      //=> 이렇게 NM-corner를 만들면 augmenting한 결과가 givenEdge주변을 완전히 둘러싸게 되어 interior edge가 된다.

    do {
        AugmentedBCFace* currFace = (AugmentedBCFace*)currNode->getEntity();

        if (    currFace->getBoundingState() == REGULAR_SIMPLEX 
             && currFace->getEdgeOrientation( givenEdge ) == rg_TRUE )  {        
            currNode = currNode->getNext();
            continue;
        }
        else {
            //  From initialFaceNode, we traverse the interior of augmented beta-complex.
            initialFaceNode = currNode;
            break;
        }
    } while( currNode != startNode );


    rg_dList< AugmentedBCFace* > facesToGenerateNMCorner;

    currNode  = initialFaceNode;
    do {
        AugmentedBCFace* currFace = (AugmentedBCFace*)currNode->getEntity();

            
        //  Some interior tetrahedra can be incident to this edge.
        //  They compose more than one groups which are only connected in this edge.
        //  The interior of one of groups is from initialFaceNode to currNode.
        
        /*
        //  For elliptic q-face
        if ( currFace->getBoundingState() == SINGULAR_SIMPLEX ) {
            if ( leftCell == INTERIOR && rightCell == INTERIOR ) {
                1. copy currFace, copy_currFace,
                2. currFace->setBoundingState( SINGULAR_SIMPLEX );
                3. copy_currFace->setBoundingState( INTERIOR_SIMPLEX );
            }
        }
        if ( nextFace->getBoundingState() == SINGULAR_SIMPLEX ) {
            if ( leftCell == INTERIOR && rightCell == INTERIOR ) {
                1. copy nextFace, copy_nextFace,
                2. nextFace->setBoundingState( SINGULAR_SIMPLEX );
                3. copy_nextFace->setBoundingState( INTERIOR_SIMPLEX );
            }
        }
        */
        
        if ( currFace->isSingularEllipticFace( m_betaValue ) ) {
            AugmentedBCFace* createdSingularEllipticFace = createSingularEllipticFace( givenEdge, currFace );
        
            currFace->setBoundingState( INTERIOR_SIMPLEX );
            fixBoundingStateOfBoundingEdgesOfFace( currFace );
            fixBoundingStateOfBoundingVerticesOfFace( currFace );

            facesToGenerateNMCorner.add( createdSingularEllipticFace );              
        }
        else
        {
            facesToGenerateNMCorner.add( currFace );
        } 
        
        
        currNode = currNode->getNext();
    } while( currNode != initialFaceNode );


    rg_INT i = 0;
    rg_INT numOfNMCorner = facesToGenerateNMCorner.getSize() - 1; // -1 하는 이유, face개수 만큼 NM-Coner를 만들면,
                                                                  // augmenting한 결과가 givenEdge주변을 완전히 둘러싸게 되어 interior edge가 된다.
    
    facesToGenerateNMCorner.reset4Loop();
    for ( i = 0; i< numOfNMCorner; ++i )    {
        facesToGenerateNMCorner.setNext4Loop();

        AugmentedBCFace* currFace = facesToGenerateNMCorner.getEntity();
        AugmentedBCFace* nextFace = facesToGenerateNMCorner.getNextEntity();

        if (    currFace->getBoundingState()              == REGULAR_SIMPLEX 
             && currFace->getEdgeOrientation( givenEdge ) == rg_FALSE ) {
            continue;
        }
        
        
        PointerToSimplex nmSimplex(EDGE_SIMPLEX, (void*)givenEdge );

        PointerToSimplex simplexOfCurrFace(FACE_SIMPLEX, currFace);
		if ( currFace->getBoundingState() == REGULAR_SIMPLEX )
			simplexOfCurrFace.setTypeOfSimplex( CELL_SIMPLEX );
        
        PointerToSimplex simplexOfNextFace(FACE_SIMPLEX, nextFace);
		if ( nextFace->getBoundingState() == REGULAR_SIMPLEX )
			simplexOfNextFace.setTypeOfSimplex( CELL_SIMPLEX );
        
        
        NonManifoldCorner* tempPtNMConf = nonManifoldCornerList.pushFront( 
                                            NonManifoldCorner( nmSimplex, simplexOfCurrFace, simplexOfNextFace ) );

        ptToNMCornerOfNMEdges[givenEdge->getID()].add( tempPtNMConf );
    }
}


void AugmentedBetaComplex::fixBoundingStateOfBoundingEdgesOfFace( AugmentedBCFace* givenFace )
{
    BetaEdge** boundingEdge = givenFace->getEdges();

    for ( int i=0; i< EIWDS_NUM_EDGE_ON_FACE; i++ ) {
        rg_dList< BetaFace* > incidentFaces;
        boundingEdge[i]->searchFacesInWholeWorld( incidentFaces );

        rg_BOOL isThisEdgeIncidentToInteriorFaceOnly = rg_TRUE;
        incidentFaces.reset4Loop();
        while ( incidentFaces.setNext4Loop() )  {
            AugmentedBCFace* currFace = (AugmentedBCFace*) incidentFaces.getEntity();

            if ( currFace->getBoundingState() != INTERIOR_SIMPLEX ) {
                isThisEdgeIncidentToInteriorFaceOnly = rg_FALSE;
                break;
            }
        }

        if ( isThisEdgeIncidentToInteriorFaceOnly == rg_TRUE )  {
            ((AugmentedBCEdge*)boundingEdge[i])->setBoundingState( INTERIOR_SIMPLEX );
        }
    }
}




void AugmentedBetaComplex::fixBoundingStateOfBoundingVerticesOfFace( AugmentedBCFace* givenFace )
{
    rg_dList< BetaVertex* > boundingVertices;
    givenFace->searchFiniteVerticesInWholeWorld( boundingVertices );

    boundingVertices.reset4Loop();
    while ( boundingVertices.setNext4Loop() ) {
        AugmentedBCVertex* currVertex = (AugmentedBCVertex*) boundingVertices.getEntity();

        rg_dList< BetaEdge* > incidentEdges;
        currVertex->searchFiniteEdgesInWholeWorld( incidentEdges );

        rg_BOOL isThisVertexIncidentToInteriorEdgeOnly = rg_TRUE;
        incidentEdges.reset4Loop();

        while ( incidentEdges.setNext4Loop() ) {
            AugmentedBCEdge* currEdge = (AugmentedBCEdge*) incidentEdges.getEntity();

            if ( currEdge->getBoundingState() != INTERIOR_SIMPLEX ) {
                isThisVertexIncidentToInteriorEdgeOnly = rg_FALSE;
                break;
            }
        }

        if ( isThisVertexIncidentToInteriorEdgeOnly == rg_TRUE ) {
            currVertex->setBoundingState( INTERIOR_SIMPLEX );
        }
    }
}






AugmentedBCFace* AugmentedBetaComplex::createSingularEllipticFace( AugmentedBCEdge* givenEdge, 
                                                                   AugmentedBCFace* singularEllipticFace )
{
    //  singularEllipticFace is a singular beta-face whose left and right cells are interior
    //  newFace : ( givenEdge, newEdge[0], newEdge[1] )
    //            ( givenEdge->m_startVertex, givenEdge->m_endVertex, newVertex )
        

    ////////////////////////////////////////////////////////////////////////////////
    //  1. create simplices( vertex, edge, face )
    AugmentedBCVertex* startVertex = (AugmentedBCVertex*) givenEdge->getStartVertex();
    AugmentedBCVertex* endVertex   = (AugmentedBCVertex*) givenEdge->getEndVertex();

    // (1.1) duplicate mate-vertex
    AugmentedBCVertex* mateVertex = (AugmentedBCVertex*) singularEllipticFace->findMateVertexOfEdge( givenEdge );
    //AugmentedBCVertex* newVertex  = m_augmentedBCVertices.add( *mateVertex );
    
    rg_INT maxIndexOfVertex = getMaxIndexOfAugmentedBCVertices();
    AugmentedBCVertex* newVertex  = m_artificialAugmentedBCVertices.add( *mateVertex );
    
    newVertex->setID( maxIndexOfVertex + 1 );    
    newVertex->setBoundingState( REGULAR_SIMPLEX );
    newVertex->setIsArtificial( rg_TRUE );


    // (1.2) create bounding edges except givenEdge
    AugmentedBCEdge* newEdge[2] = { rg_NULL, rg_NULL };

    rg_INT maxIndexOfEdge = getMaxIndexOfAugmentedBCEdges();
    newEdge[0] = m_artificialAugmentedBCEdges.add( AugmentedBCEdge( endVertex, newVertex ) );
    newEdge[0]->setID( maxIndexOfEdge + 1 );
    newEdge[0]->setBoundingState( REGULAR_SIMPLEX );
    newEdge[0]->setIsArtificial( rg_TRUE );
    
    newEdge[1] = m_artificialAugmentedBCEdges.add( AugmentedBCEdge( newVertex, startVertex ) );
    newEdge[1]->setID( maxIndexOfEdge + 2 );
    newEdge[1]->setBoundingState( REGULAR_SIMPLEX );
    newEdge[1]->setIsArtificial( rg_TRUE );


    // (1.3) create duplicated singularEllipticFace
    rg_INT maxIndexOfFace = getMaxIndexOfAugmentedBCFaces();
    AugmentedBCFace* newFace = m_artificialAugmentedBCFaces.add(  AugmentedBCFace( givenEdge, newEdge[0], newEdge[1] )  );

    newFace->setEdgeOrientation(0, rg_TRUE);
    newFace->setEdgeOrientation(1, rg_TRUE);
    newFace->setEdgeOrientation(2, rg_TRUE);
    
    newFace->setID( maxIndexOfFace + 1 );
    newFace->setBoundingState( SINGULAR_SIMPLEX );
    newFace->setIsArtificial( rg_TRUE );
    

    ////////////////////////////////////////////////////////////////////////////////
    //  2. build topology lower-dimensional simplex to higher-dimensional simplex
    newVertex->setFirstEdge( newEdge[0] );
    newVertex->setFirstCell( rg_NULL );

    newEdge[0]->setFirstFace( newFace );
    newEdge[1]->setFirstFace( newFace );    

    newFace->setLeftCell(  rg_NULL );
    newFace->setRightCell( rg_NULL );


    ////////////////////////////////////////////////////////////////////////////////
    //  3. set exterior index by regular faces of left cell and right cell
    BetaCell* leftCell = singularEllipticFace->getLeftCell();
    BetaFace* facesIncidentToGivenEdgeInLeftCell[2];    
    leftCell->findFacesIncidentToEdge( givenEdge, facesIncidentToGivenEdgeInLeftCell );
    
    BetaFace* regularFaceOfLeftCell = rg_NULL;
    if ( facesIncidentToGivenEdgeInLeftCell[0] == singularEllipticFace ) {
        regularFaceOfLeftCell = facesIncidentToGivenEdgeInLeftCell[1];
    }
    else {
        regularFaceOfLeftCell = facesIncidentToGivenEdgeInLeftCell[0];
    }

    newFace->setIndexOfExteriorRegionToLeftCell( ((AugmentedBCFace*) regularFaceOfLeftCell)->getIndexOfExteriorRegionToRightCell() );
                                                                                         //->getIndexOfExteriorRegionToRightCell()  regular face이므로 right방향으르 가져오는 것이 맞다.
    
    BetaCell* rightCell = singularEllipticFace->getRightCell();
    BetaFace* facesIncidentToGivenEdgeInRightCell[2];    
    rightCell->findFacesIncidentToEdge( givenEdge, facesIncidentToGivenEdgeInRightCell );
    
    BetaFace* regularFaceOfRightCell = rg_NULL;
    if ( facesIncidentToGivenEdgeInRightCell[0] == singularEllipticFace ) {
        regularFaceOfRightCell = facesIncidentToGivenEdgeInRightCell[1];
    }
    else {
        regularFaceOfRightCell = facesIncidentToGivenEdgeInRightCell[0];
    }

    newFace->setIndexOfExteriorRegionToRightCell( ((AugmentedBCFace*) regularFaceOfRightCell)->getIndexOfExteriorRegionToRightCell() );

        
    
    return newFace;
}


void AugmentedBetaComplex::createPtOfSingularFaceOrCellBoundedByRegularFace( AugmentedBCFace* givenFace, PointerToSimplex& ptToSingularFaceOrCell ) const
{
    if( givenFace->getBoundingState() == REGULAR_SIMPLEX )
    {
        ptToSingularFaceOrCell.setTypeOfSimplex( CELL_SIMPLEX );
        ptToSingularFaceOrCell.setSimplex( (void*)(givenFace->getLeftCell() ) );
    }
    else //if ( givenFace->getBoundingState() == SINGULAR_SIMPLEX )
    {
        ptToSingularFaceOrCell.setTypeOfSimplex( FACE_SIMPLEX );
        ptToSingularFaceOrCell.setSimplex( (void*)givenFace );
    }    
}











void AugmentedBetaComplex::detachAndRemoveExtraneousSimplexes()
{
    // regular edge의 incident face가 모두 extraneous이면 해당 edge를 extraneous로 바꾸는 작업을 한다.
    //setEdgeAsExtraneousIfEdgeIncidentToExtraneousFacesOnly();

    detachExtraneousEdgesFromAllVertices(); 

	detachExtraneousCellsFromAllVertices();

    detachExtraneousFacesFromAllEdges();

    detachExtraneousCellsFromAllFaces();
    

    removeExtraneousCells();

    removeExtraneousFaces();

    removeExtraneousEdges();
    
    removeExtraneousVertices();

}

void AugmentedBetaComplex::setEdgeAsSingularIfEdgeIncidentToExtraneousFacesOnly()
{
    AugmentedBCEdge* currEdge     = rg_NULL;
    rg_INT   currEdgeType = -1;
    m_augmentedBCEdges.reset4Loop();
    while ( m_augmentedBCEdges.setNext4Loop() )
    {
        currEdge     = m_augmentedBCEdges.getpEntity();
        currEdgeType = currEdge->getBoundingState();

        if ( currEdgeType == REGULAR_SIMPLEX ) {            
            rg_INT boundingStateOfFirstFace = ((AugmentedBCFace*)currEdge->getFirstFace())->getBoundingState();
            if (    boundingStateOfFirstFace == EXTRANEOUS_SIMPLEX 
                 || boundingStateOfFirstFace == FALSE_EXTERIOR_SIMPLEX )
            {
                rg_dList< BetaFace* > incidentFaceList;
                currEdge->searchFiniteFacesInWholeWorld( incidentFaceList );

                rg_BOOL isThisEdgeIncidentToExtraneousFacesOnly = rg_TRUE;

                incidentFaceList.reset4Loop();
                while (incidentFaceList.setNext4Loop() )
                {
                    AugmentedBCFace* currIncidentFace = (AugmentedBCFace*) incidentFaceList.getEntity();

                    rg_INT currIncidentFaceBoundingState = currIncidentFace->getBoundingState();
                    if (    currIncidentFaceBoundingState != EXTRANEOUS_SIMPLEX 
                         && currIncidentFaceBoundingState != FALSE_EXTERIOR_SIMPLEX ) {
                        isThisEdgeIncidentToExtraneousFacesOnly = rg_FALSE;
                    }
                }

                ///////////////////////////////////////////////////////////////////
                // regular edge인데 incident한 face가 모두 extraneous인 경우가 있다.
                // 이 경우 일단 이 regular edge를 extraneous로 바꾼다.                
                //
                if ( isThisEdgeIncidentToExtraneousFacesOnly ) {
                    currEdge->setBoundingState( SINGULAR_SIMPLEX );
                }


                //
                ////////////////////////////////////////////////////////////////////////
            }
        }
    }
}

void AugmentedBetaComplex::detachExtraneousEdgesFromAllVertices()
{
    AugmentedBCVertex* currVertex = rg_NULL;
    rg_INT     currVertexType = -1;

    m_augmentedBCVertices.reset4Loop();
    while ( m_augmentedBCVertices.setNext4Loop() )
    {
        currVertex     = m_augmentedBCVertices.getpEntity();
        currVertexType = currVertex->getBoundingState();        

        if ( currVertexType == SINGULAR_SIMPLEX ) {
            currVertex->setFirstEdge( rg_NULL );
        }
        else if (    currVertexType == EXTRANEOUS_SIMPLEX
                  || currVertexType == FALSE_EXTERIOR_SIMPLEX
                  || currVertexType == INTERIOR_SIMPLEX ) {
            continue;
        }
        else { // if currVertexType == REGULAR_SIMPLEX
            rg_INT boundingStateOfFirstEdge = ((AugmentedBCEdge*)currVertex->getFirstEdge())->getBoundingState();
            if (    boundingStateOfFirstEdge == EXTRANEOUS_SIMPLEX
                 || boundingStateOfFirstEdge == FALSE_EXTERIOR_SIMPLEX) {
                rg_dList< BetaEdge* > incidentEdgeList;
			    currVertex->searchFiniteEdgesInWholeWorld( incidentEdgeList );

                currVertex->setFirstEdge( rg_NULL );
                incidentEdgeList.reset4Loop();
                while (incidentEdgeList.setNext4Loop() )
                {
                    AugmentedBCEdge* currIncidentEdge = (AugmentedBCEdge*) incidentEdgeList.getEntity();

                    rg_INT currIncidentEdgeBoundingState = currIncidentEdge->getBoundingState();
                    if (    currIncidentEdgeBoundingState != EXTRANEOUS_SIMPLEX 
                         && currIncidentEdgeBoundingState != FALSE_EXTERIOR_SIMPLEX )
                    {
                        currVertex->setFirstEdge( currIncidentEdge );
                        break;
                    }
                }
           }
        }
    }
}



void AugmentedBetaComplex::detachExtraneousCellsFromAllVertices()
{
    AugmentedBCVertex* currVertex = rg_NULL;
    rg_INT     currVertexType = -1;

    m_augmentedBCVertices.reset4Loop();
    while ( m_augmentedBCVertices.setNext4Loop() )
    {
        currVertex     = m_augmentedBCVertices.getpEntity();
        currVertexType = currVertex->getBoundingState();        

        if ( currVertexType == SINGULAR_SIMPLEX )   {
            currVertex->setFirstCell( rg_NULL );
            continue;
        }
        else if (    currVertexType == EXTRANEOUS_SIMPLEX
                  || currVertexType == FALSE_EXTERIOR_SIMPLEX
                  || currVertexType == INTERIOR_SIMPLEX ) {
            continue;
        }
        else    { // currVertexType == REGULAR_SIMPLEX
            if ( currVertex->getFirstCell() == rg_NULL )    {
                continue;
            }
            
            rg_INT boundingStateOfFirstCell = ((AugmentedBCCell*)currVertex->getFirstCell())->getBoundingState();
            if (    boundingStateOfFirstCell == EXTRANEOUS_SIMPLEX 
                 || boundingStateOfFirstCell == FALSE_EXTERIOR_SIMPLEX ) {            
                rg_dList< BetaCell* > incidentCellList;
			    currVertex->searchFiniteCellsInWholeWorld( incidentCellList );	

                currVertex->setFirstCell( rg_NULL );
                incidentCellList.reset4Loop();
                while (incidentCellList.setNext4Loop() )
                {
                    AugmentedBCCell* currIncidentCell = (AugmentedBCCell*) incidentCellList.getEntity();

                    rg_INT currIncidentCellBoundingState = currIncidentCell->getBoundingState();
                    if (    currIncidentCellBoundingState != EXTRANEOUS_SIMPLEX 
                         && currIncidentCellBoundingState != FALSE_EXTERIOR_SIMPLEX ) {
                        currVertex->setFirstCell( currIncidentCell );
                        break;
                    }
                }
            }
        }    
    }
}


void AugmentedBetaComplex::detachExtraneousFacesFromAllEdges()
{
    AugmentedBCEdge* currEdge     = rg_NULL;
    rg_INT   currEdgeType = -1;
    m_augmentedBCEdges.reset4Loop();
    while ( m_augmentedBCEdges.setNext4Loop() )
    {
        currEdge     = m_augmentedBCEdges.getpEntity();
        currEdgeType = currEdge->getBoundingState();


        if ( currEdgeType == SINGULAR_SIMPLEX ) {
            currEdge->setFirstFace( rg_NULL );
        }
        else if (    currEdgeType == EXTRANEOUS_SIMPLEX
                  || currEdgeType == FALSE_EXTERIOR_SIMPLEX
                  || currEdgeType == INTERIOR_SIMPLEX ) {
               continue;
        }
        else { //if currEdgeType == REGULAR_SIMPLEX
            rg_INT boundingStateOfFirstFace = ((AugmentedBCFace*)currEdge->getFirstFace())->getBoundingState();
            if (    boundingStateOfFirstFace == EXTRANEOUS_SIMPLEX 
                 || boundingStateOfFirstFace == FALSE_EXTERIOR_SIMPLEX )
            {
                rg_dList< BetaFace* > incidentFaceList;
                currEdge->searchFiniteFacesInWholeWorld( incidentFaceList );

                currEdge->setFirstFace( rg_NULL );
                incidentFaceList.reset4Loop();
                while (incidentFaceList.setNext4Loop() )
                {
                    AugmentedBCFace* currIncidentFace = (AugmentedBCFace*) incidentFaceList.getEntity();

                    rg_INT currIncidentFaceBoundingState = currIncidentFace->getBoundingState();
                    if (    currIncidentFaceBoundingState != EXTRANEOUS_SIMPLEX 
                         && currIncidentFaceBoundingState != FALSE_EXTERIOR_SIMPLEX ) {
                        currEdge->setFirstFace( currIncidentFace );
                        break;
                    }
                }                
            }
        }
    }
}
            

void AugmentedBetaComplex::detachExtraneousCellsFromAllFaces()
{
    m_augmentedBCFaces.reset4Loop();
    while ( m_augmentedBCFaces.setNext4Loop() )
    {
		AugmentedBCFace* currFace	  = m_augmentedBCFaces.getpEntity();
        //rg_INT			 currFaceType = currFace->getBoundingState();

        AugmentedBCCell* leftCell  = (AugmentedBCCell*) currFace->getLeftCell();
        AugmentedBCCell* rightCell = (AugmentedBCCell*) currFace->getRightCell();

        if ( leftCell != rg_NULL )  {
            rg_INT leftCellBoundingState = leftCell->getBoundingState();
            if (    leftCellBoundingState == EXTRANEOUS_SIMPLEX 
                 || leftCellBoundingState == FALSE_EXTERIOR_SIMPLEX ) {
                currFace->setLeftCell(  rg_NULL );
            }
        }

        if ( rightCell != rg_NULL )  {
            rg_INT rightCellBoundingState = rightCell->getBoundingState();
            if (    rightCellBoundingState == EXTRANEOUS_SIMPLEX
                 || rightCellBoundingState == FALSE_EXTERIOR_SIMPLEX ) {
                currFace->setRightCell(  rg_NULL );
            }
        }        
    }
        
	/*
    //  2010-02-11 위의 내용으로 수정(이창희, 조영송)
    AugmentedBCFace* currFace     = rg_NULL;
    rg_INT   currFaceType = -1;
    m_augmentedBCFaces.reset4Loop();
    while ( m_augmentedBCFaces.setNext4Loop() )
    {
        currFace     = m_augmentedBCFaces.getpEntity();
        currFaceType = currFace->getBoundingState();

        if ( currFaceType == INTERIOR_SIMPLEX 
            || currFaceType == EXTRANEOUS_SIMPLEX )
            continue;
        else if ( currFace->getBoundingState() == SINGULAR_SIMPLEX )
        {
            currFace->setLeftCell(  rg_NULL );
            currFace->setRightCell( rg_NULL);
        }
        else if ( currFace->getLeftCell() == rg_NULL )
        {
            currFace->setLeftCell( currFace->getRightCell() );
            currFace->setRightCell( rg_NULL );

            currFace->reverseEdgeOrientations();
        }
        else if ( ((AugmentedBCCell*)currFace->getLeftCell())->getBoundingState() == EXTRANEOUS_SIMPLEX )
        {
            currFace->setLeftCell( currFace->getRightCell() );
            currFace->setRightCell( rg_NULL );

            currFace->reverseEdgeOrientations();
        }
        else if ( currFace->getRightCell() == rg_NULL )
            currFace->setRightCell( rg_NULL );
        else if ( ((AugmentedBCCell*)currFace->getRightCell())->getBoundingState() == EXTRANEOUS_SIMPLEX )
                currFace->setRightCell( rg_NULL );            
        
    }
	*/
}






void AugmentedBetaComplex::removeExtraneousCells()
{
    rg_dNode< AugmentedBCCell >* startNode = m_augmentedBCCells.getFirstpNode();
    rg_dNode< AugmentedBCCell >* currNode  = startNode;

    rg_dList< rg_dNode< AugmentedBCCell >* > toBeKilledNodeList;
    
    AugmentedBCCell* currCell =  rg_NULL;
    do {
        currCell = currNode->getpEntity();

        rg_INT currCellBoundingState = currCell->getBoundingState();
        if (    currCellBoundingState == EXTRANEOUS_SIMPLEX 
             || currCellBoundingState == FALSE_EXTERIOR_SIMPLEX )
            toBeKilledNodeList.add( currNode );        

        currNode = currNode->getNext();
    	
    } while( currNode != startNode );

    rg_dNode< AugmentedBCCell >* currToBeKilledNode = rg_NULL;
    toBeKilledNodeList.reset4Loop();
    while ( toBeKilledNodeList.setNext4Loop() )
    {
        currToBeKilledNode = toBeKilledNodeList.getEntity();
        m_augmentedBCCells.kill( currToBeKilledNode );
    }
}




void AugmentedBetaComplex::removeExtraneousFaces()
{
    rg_dNode< AugmentedBCFace >* startNode = m_augmentedBCFaces.getFirstpNode();
    rg_dNode< AugmentedBCFace >* currNode  = startNode;

    rg_dList< rg_dNode< AugmentedBCFace >* > toBeKilledNodeList;

    AugmentedBCFace* currFace =  rg_NULL;

    do {
        currFace = currNode->getpEntity();
       
        rg_INT currFaceBoundingState = currFace->getBoundingState();
        if (    currFaceBoundingState == EXTRANEOUS_SIMPLEX 
             || currFaceBoundingState == FALSE_EXTERIOR_SIMPLEX )
            toBeKilledNodeList.add( currNode );
        
        currNode = currNode->getNext();
    	
    } while( currNode != startNode );

    rg_dNode< AugmentedBCFace >* currToBeKilledNode = rg_NULL;
    toBeKilledNodeList.reset4Loop();
    while ( toBeKilledNodeList.setNext4Loop() )
    {
        currToBeKilledNode = toBeKilledNodeList.getEntity();
        m_augmentedBCFaces.kill( currToBeKilledNode );
    }
}




void AugmentedBetaComplex::removeExtraneousEdges()
{
    rg_dNode< AugmentedBCEdge >* startNode = m_augmentedBCEdges.getFirstpNode();
    rg_dNode< AugmentedBCEdge >* currNode  = startNode;

    rg_dList< rg_dNode< AugmentedBCEdge >* > toBeKilledNodeList;

    AugmentedBCEdge* currEdge =  rg_NULL;

    do {
        currEdge = currNode->getpEntity();

        rg_INT currEdgeBoundingState = currEdge->getBoundingState();
        if (    currEdgeBoundingState == EXTRANEOUS_SIMPLEX
             || currEdgeBoundingState == FALSE_EXTERIOR_SIMPLEX )        
            toBeKilledNodeList.add( currNode );
        
        currNode = currNode->getNext();
    	
    } while( currNode != startNode );

    rg_dNode< AugmentedBCEdge >* currToBeKilledNode = rg_NULL;
    toBeKilledNodeList.reset4Loop();
    while ( toBeKilledNodeList.setNext4Loop() )
    {
        currToBeKilledNode = toBeKilledNodeList.getEntity();
        m_augmentedBCEdges.kill( currToBeKilledNode );
    }
}




void AugmentedBetaComplex::removeExtraneousVertices()
{
    rg_dNode< AugmentedBCVertex >* startNode = m_augmentedBCVertices.getFirstpNode();
    rg_dNode< AugmentedBCVertex >* currNode  = startNode;

    rg_dList< rg_dNode< AugmentedBCVertex >* > toBeKilledNodeList;

    AugmentedBCVertex* currVertex =  rg_NULL;

    do {
        currVertex = currNode->getpEntity();

        rg_INT currVertexBoundingState = currVertex->getBoundingState();
        if (    currVertexBoundingState == EXTRANEOUS_SIMPLEX 
             || currVertexBoundingState == FALSE_EXTERIOR_SIMPLEX )        
            toBeKilledNodeList.add( currNode );
        
        currNode = currNode->getNext();
    	
    } while( currNode != startNode );

    rg_dNode< AugmentedBCVertex >* currToBeKilledNode = rg_NULL;
    toBeKilledNodeList.reset4Loop();
    while ( toBeKilledNodeList.setNext4Loop() )
    {
        currToBeKilledNode = toBeKilledNodeList.getEntity();
        m_augmentedBCVertices.kill( currToBeKilledNode );
    }
}




// void AugmentedBetaComplex::augment(rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges,
//                                                                rg_dList< NonManifoldCorner >& nonManifoldCornerList)
// {
// 	rg_BOOL isFirstTT_DV = rg_TRUE;
// 
//     while( nonManifoldCornerList.getSize() > 0 )
//     {
// 		NonManifoldCorner currNonManifoldConf = nonManifoldCornerList.popFront();
//         
// 
//         deletePTToNMCornerOfGivenNMCorner( &currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges,
//                                                                nonManifoldCornerList  );        
//         
//         rg_INT currNonManifoldConfType = currNonManifoldConf.getTypeOfNonManifoldConfiguration();
//         NonManifoldCorner raisedNonManifoldConf;
// 
// 
//         switch( currNonManifoldConfType )
//         {
//         case EE_V_CONFIGURATION:
//             raisedNonManifoldConf = resolve_EE_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices );
//             nonManifoldCornerList.pushFront( raisedNonManifoldConf );
//             break;
// 
//         case EF_V_CONFIGURATION:
//             raisedNonManifoldConf = resolve_EF_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
//             nonManifoldCornerList.pushFront( raisedNonManifoldConf );
//             break;
// 
//         case ET_V_CONFIGURATION:
//             raisedNonManifoldConf = resolve_ET_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices );
//             nonManifoldCornerList.pushFront( raisedNonManifoldConf );
//             break;
// 
//         case FF_V_CONFIGURATION:
//             raisedNonManifoldConf = resolve_FF_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
//             nonManifoldCornerList.pushFront( raisedNonManifoldConf );
//             break;
// 
//         case FT_V_CONFIGURATION:
//             raisedNonManifoldConf = resolve_FT_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
//             nonManifoldCornerList.pushFront( raisedNonManifoldConf );
//             break;
// 
//         case TT_V_CONFIGURATION:
//             raisedNonManifoldConf = resolve_TT_V_Configuration( currNonManifoldConf );
//             nonManifoldCornerList.pushFront( raisedNonManifoldConf );
//             break;
// 
//         case FF_E_CONFIGURATION:
//             resolve_FF_E_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
//             break;
// 
//         case FT_E_CONFIGURATION:
//             resolve_FT_E_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
//             break;
// 
//         case TF_E_CONFIGURATION:
//             resolve_TF_E_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
//             break;
// 
//         case TT_E_CONFIGURATION:
//             resolve_TT_E_Configuration( currNonManifoldConf );
//             break;
// 
//         case TT_D_V_CONFIGURATION:
// 			if ( isFirstTT_DV == rg_TRUE )
// 			{
// 				isFirstTT_DV = rg_FALSE;
// 				checkComparingNumExtraneousCellGroupAndNumRegularFaceGroup();
// 			}
// 
//             resolve_TT_D_V_Configuration( currNonManifoldConf );
//             break;        
//             
//         default:
//             break;
//         }
//     }
// }




void AugmentedBetaComplex::augment(rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges,
                                                               rg_dList< NonManifoldCorner >& nonManifoldCornerList)
{
    rg_BOOL isFirstTT_DV = rg_TRUE;
    
	while( nonManifoldCornerList.getSize() > 0 )
    {
		NonManifoldCorner currNonManifoldConf = nonManifoldCornerList.popFront();
        

        deletePTToNMCornerOfGivenNMCorner( &currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges,
                                                               nonManifoldCornerList  );        
        
        rg_INT currNonManifoldConfType = currNonManifoldConf.getTypeOfNonManifoldConfiguration();
        NonManifoldCorner raisedNonManifoldConf;


        switch( currNonManifoldConfType )
        {
            // nm-simplex: VERTEX
        case EE_V_CONFIGURATION:
            raisedNonManifoldConf = resolve_EE_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices );
            nonManifoldCornerList.pushFront( raisedNonManifoldConf );
            break;

        case EF_V_CONFIGURATION:
            raisedNonManifoldConf = resolve_EF_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
            nonManifoldCornerList.pushFront( raisedNonManifoldConf );
            break;

        case ET_V_CONFIGURATION:
            raisedNonManifoldConf = resolve_ET_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices );
            nonManifoldCornerList.pushFront( raisedNonManifoldConf );
            break;

        case FF_V_CONFIGURATION:
            raisedNonManifoldConf = resolve_FF_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
            nonManifoldCornerList.pushFront( raisedNonManifoldConf );
            break;

        case FT_V_CONFIGURATION:
            raisedNonManifoldConf = resolve_FT_V_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
            nonManifoldCornerList.pushFront( raisedNonManifoldConf );
            break;

        case TT_V_CONFIGURATION:
            raisedNonManifoldConf = resolve_TT_V_Configuration( currNonManifoldConf );
            nonManifoldCornerList.pushFront( raisedNonManifoldConf );
            break;




            // nm-simplex: EDGE
        case FF_E_CONFIGURATION:
            resolve_FF_E_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
            break;

        case FT_E_CONFIGURATION:
            resolve_FT_E_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
            break;

        case TF_E_CONFIGURATION:
            resolve_TF_E_Configuration( currNonManifoldConf, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
            break;

        case TT_E_CONFIGURATION:
            resolve_TT_E_Configuration( currNonManifoldConf );
            break;

        case TT_D_V_CONFIGURATION:	
            if ( isFirstTT_DV == rg_TRUE )
			{
				isFirstTT_DV = rg_FALSE;
				//checkComparingNumExtraneousCellGroupAndNumRegularFaceGroup();
			}
            resolve_TT_D_V_Configuration( currNonManifoldConf );
            break;        
            
        default:
            break;
        }
    }
}




void AugmentedBetaComplex::deletePTToNMCornerOfGivenNMCorner( NonManifoldCorner* givenNMCorner, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges,
                                                               rg_dList< NonManifoldCorner >& nonManifoldCornerList  )
{
    PointerToSimplex* nmSimplexOfFirstNMConf = givenNMCorner->getpNMSimplex();

    rg_dNode< NonManifoldCorner* >* toBeDeletedNode = rg_NULL;

    if ( nmSimplexOfFirstNMConf->getTypeOfSimplex() == VERTEX_SIMPLEX )
    {
        AugmentedBCVertex* nmVertex = (AugmentedBCVertex*)( nmSimplexOfFirstNMConf->getSimplex() );

        rg_INT idOfNMVertex = nmVertex->getID();
        ptToNMCornerOfNMVertices[ idOfNMVertex ].reset4Loop();
        NonManifoldCorner* currNMConf = rg_NULL;

        while( ptToNMCornerOfNMVertices[ idOfNMVertex ].setNext4Loop() )
        {
            currNMConf = ptToNMCornerOfNMVertices[ idOfNMVertex ].getEntity();

            if ( currNMConf == givenNMCorner )
            {
                toBeDeletedNode = ptToNMCornerOfNMVertices[ idOfNMVertex ].getCurrentpNode();
                break;
            }
        }

        if ( toBeDeletedNode != rg_NULL )
            ptToNMCornerOfNMVertices[ idOfNMVertex ].kill( toBeDeletedNode );
    }
    else //if ( nmSimplexOfFirstNMConf->getTypeOfSimplex() == EDGE_SIMPLEX )
    {
        AugmentedBCEdge* nmEdge = (AugmentedBCEdge*)( nmSimplexOfFirstNMConf->getSimplex() );

        if ( nmEdge->isArtificial() == rg_TRUE )
            return;

        rg_INT idOfNMEdge = nmEdge->getID();
        ptToNMCornerOfNMEdges[ idOfNMEdge ].reset4Loop();
        NonManifoldCorner* currNMConf = rg_NULL;

        while( ptToNMCornerOfNMEdges[ idOfNMEdge ].setNext4Loop() )
        {
            currNMConf = ptToNMCornerOfNMEdges[ idOfNMEdge ].getEntity();

            if ( currNMConf == givenNMCorner )
            {
                toBeDeletedNode = ptToNMCornerOfNMEdges[ idOfNMEdge ].getCurrentpNode();
                break;
            }
        }

        if ( toBeDeletedNode != rg_NULL )
            ptToNMCornerOfNMEdges[ idOfNMEdge ].kill( toBeDeletedNode );
    }
}


NonManifoldCorner AugmentedBetaComplex::resolve_EE_V_Configuration(const NonManifoldCorner& givenNMConf, 
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices )
{
    AugmentedBCVertex* nmVertex      = (AugmentedBCVertex*)( givenNMConf.getNMSimplex().getSimplex() );
    AugmentedBCEdge* givenFirstEdge  = (AugmentedBCEdge*)(   givenNMConf.getFirstSimplex().getSimplex() );
    AugmentedBCEdge* givenSecondEdge = (AugmentedBCEdge*)(   givenNMConf.getSecondSimplex().getSimplex() );

    AugmentedBCEdge* raisedNMEdge = rg_NULL;
    AugmentedBCFace* raisedFace1  = rg_NULL;
    AugmentedBCFace* raisedFace2  = rg_NULL;


    rg_INT indexOfCommonExteriorRegion = givenFirstEdge->getIndexOfExteriorRegion();


    O_EE_V_operator( nmVertex, givenFirstEdge, givenSecondEdge, raisedNMEdge, raisedFace1, raisedFace2 );


    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenFirstEdge,  raisedFace1, ptToNMCornerOfNMVertices );
    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenSecondEdge, raisedFace2, ptToNMCornerOfNMVertices );
    

    setIndexOfExteriorRegionOfSingularFace( indexOfCommonExteriorRegion, raisedFace1 );
    setIndexOfExteriorRegionOfSingularFace( indexOfCommonExteriorRegion, raisedFace2 );
    

    NonManifoldCorner raisedConf( raisedNMEdge, raisedFace1, raisedFace2 );
    return raisedConf;
}


NonManifoldCorner AugmentedBetaComplex::resolve_EF_V_Configuration(const NonManifoldCorner& givenNMConf, 
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges )
{
    AugmentedBCVertex* nmVertex = (AugmentedBCVertex*)( givenNMConf.getNMSimplex().getSimplex() );    
    
    AugmentedBCEdge* givenEdge = (AugmentedBCEdge*) givenNMConf.getEdgeSimplex();
    AugmentedBCFace* givenFace = (AugmentedBCFace*) givenNMConf.getFaceSimplex();
    
    AugmentedBCEdge* raisedNMEdge = rg_NULL;
    AugmentedBCFace* raisedFace   = rg_NULL;
    AugmentedBCCell* raisedTetra  = rg_NULL;


    rg_INT indexOfCommonExteriorRegion = givenEdge->getIndexOfExteriorRegion();


    if ( givenFace->getIndexOfExteriorRegionToRightCell() != indexOfCommonExteriorRegion ) {
        givenFace->reverse();
    }


    O_EF_V_operator( nmVertex, givenEdge, givenFace, raisedNMEdge, raisedFace, raisedTetra );


    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenEdge, raisedFace, ptToNMCornerOfNMVertices );
    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenFace, raisedTetra, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
    

    setIndexOfExteriorRegionOfSingularFace( indexOfCommonExteriorRegion, raisedFace  );
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, givenFace, raisedTetra );


    NonManifoldCorner raisedConf( raisedNMEdge, raisedFace, givenFace );
	raisedConf.getpSecondSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
    return raisedConf;
}


NonManifoldCorner AugmentedBetaComplex::resolve_ET_V_Configuration( const NonManifoldCorner& givenNMConf, 
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices )
{
    AugmentedBCVertex* nmVertex   = (AugmentedBCVertex*)( givenNMConf.getNMSimplex().getSimplex() );    
    
    AugmentedBCEdge* givenEdge = (AugmentedBCEdge*) givenNMConf.getEdgeSimplex();
    AugmentedBCFace* givenFace = (AugmentedBCFace*) givenNMConf.getCellSimplex();
    
    AugmentedBCEdge*   raisedNMEdge = rg_NULL;
    AugmentedBCFace*   raisedFace   = rg_NULL;
    AugmentedBCCell*   raisedTetra  = rg_NULL;


    rg_INT indexOfCommonExteriorRegion = givenEdge->getIndexOfExteriorRegion();


    AugmentedBCFace* regularFaceOnExactExteriorRegion = findRegularFaceTouchedExactExterior( nmVertex, givenFace, indexOfCommonExteriorRegion );


    O_ET_V_operator( nmVertex, givenEdge, regularFaceOnExactExteriorRegion, raisedNMEdge, raisedFace, raisedTetra );


    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenEdge,  raisedFace, ptToNMCornerOfNMVertices );
    

    setIndexOfExteriorRegionOfSingularFace( indexOfCommonExteriorRegion, raisedFace );
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, regularFaceOnExactExteriorRegion, raisedTetra );


    NonManifoldCorner raisedConf( raisedNMEdge, raisedFace, regularFaceOnExactExteriorRegion );
	raisedConf.getpSecondSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
    return raisedConf;
}



NonManifoldCorner AugmentedBetaComplex::resolve_FF_V_Configuration(const NonManifoldCorner& givenNMConf, 
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges )
{
    AugmentedBCVertex* nmVertex = (AugmentedBCVertex*)( givenNMConf.getNMSimplex().getSimplex() );    
    AugmentedBCFace* givenFace1 = (AugmentedBCFace*)( givenNMConf.getFirstSimplex().getSimplex() );
    AugmentedBCFace* givenFace2 = (AugmentedBCFace*)( givenNMConf.getSecondSimplex().getSimplex() );

    AugmentedBCEdge* raisedNMEdge = rg_NULL;
    AugmentedBCCell* raisedTetra1 = rg_NULL;
    AugmentedBCCell* raisedTetra2 = rg_NULL;

    rg_INT indexOfCommonExteriorRegion = findIndexOfCommonExteriorRegion( givenFace1, givenFace2 );

    if ( givenFace1->getIndexOfExteriorRegionToRightCell() != indexOfCommonExteriorRegion ) {
        givenFace1->reverse();
    }

    if ( givenFace2->getIndexOfExteriorRegionToRightCell() != indexOfCommonExteriorRegion ) {
        givenFace2->reverse();
    }


    O_FF_V_operator( nmVertex, givenFace1, givenFace2, raisedNMEdge, raisedTetra1, raisedTetra2 );


    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenFace1, raisedTetra1, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenFace2, raisedTetra2, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
    

    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, givenFace1, raisedTetra1 );
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, givenFace2, raisedTetra2 );


    NonManifoldCorner raisedConf( raisedNMEdge, givenFace1, givenFace2 );
	raisedConf.getpFirstSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
	raisedConf.getpSecondSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
    return raisedConf;
}


NonManifoldCorner AugmentedBetaComplex::resolve_FT_V_Configuration(const NonManifoldCorner& givenNMConf, 
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges )
{
    AugmentedBCVertex* nmVertex   = (AugmentedBCVertex*)( givenNMConf.getNMSimplex().getSimplex() );

    AugmentedBCFace* givenFace        = (AugmentedBCFace*) givenNMConf.getFaceSimplex();
    AugmentedBCFace* givenRegularFace = (AugmentedBCFace*) givenNMConf.getCellSimplex();
    
    AugmentedBCEdge* raisedNMEdge = rg_NULL;
    AugmentedBCCell* raisedTetra1 = rg_NULL;
    AugmentedBCCell* raisedTetra2 = rg_NULL;

    rg_INT indexOfCommonExteriorRegion = findIndexOfCommonExteriorRegion( givenFace, givenRegularFace );


    if ( givenFace->getIndexOfExteriorRegionToRightCell() != indexOfCommonExteriorRegion ) {
        givenFace->reverse();
    }
    AugmentedBCFace* regularFaceOnExactExteriorRegion = findRegularFaceTouchedExactExterior( nmVertex, givenRegularFace, indexOfCommonExteriorRegion );



    O_FT_V_operator( nmVertex, givenFace, regularFaceOnExactExteriorRegion, raisedNMEdge, raisedTetra1, raisedTetra2 );


    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenFace, raisedTetra1, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
    

    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, givenFace, raisedTetra1 );
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, regularFaceOnExactExteriorRegion, raisedTetra2 );



    NonManifoldCorner raisedConf( raisedNMEdge, givenFace, regularFaceOnExactExteriorRegion );
	raisedConf.getpFirstSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
	raisedConf.getpSecondSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
    return raisedConf;
}


NonManifoldCorner AugmentedBetaComplex::resolve_TT_V_Configuration(const NonManifoldCorner& givenNMConf )
{
    AugmentedBCVertex* nmVertex			 = (AugmentedBCVertex*)( givenNMConf.getNMSimplex().getSimplex() );
    AugmentedBCFace*   givenRegularFace1 = (AugmentedBCFace*)( givenNMConf.getFirstSimplex().getSimplex() );
    AugmentedBCFace*   givenRegularFace2 = (AugmentedBCFace*)( givenNMConf.getSecondSimplex().getSimplex() );

    AugmentedBCEdge* raisedNMEdge = rg_NULL;
    AugmentedBCCell* raisedTetra1 = rg_NULL;
    AugmentedBCCell* raisedTetra2 = rg_NULL;

    
    rg_INT indexOfCommonExteriorRegion = findIndexOfCommonExteriorRegion( givenRegularFace1, givenRegularFace2 );


    AugmentedBCFace* firstRegularFaceOnExactExteriorRegion = findRegularFaceTouchedExactExterior( nmVertex, givenRegularFace1, indexOfCommonExteriorRegion );
    AugmentedBCFace* secondRegularFaceOnExactExteriorRegion = findRegularFaceTouchedExactExterior( nmVertex, givenRegularFace2, indexOfCommonExteriorRegion );



    O_TT_V_operator( nmVertex, firstRegularFaceOnExactExteriorRegion, secondRegularFaceOnExactExteriorRegion, raisedNMEdge, raisedTetra1, raisedTetra2 );

    
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, firstRegularFaceOnExactExteriorRegion, raisedTetra1 );
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, secondRegularFaceOnExactExteriorRegion, raisedTetra2 );


    NonManifoldCorner raisedConf( raisedNMEdge, firstRegularFaceOnExactExteriorRegion, secondRegularFaceOnExactExteriorRegion );
	raisedConf.getpFirstSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
	raisedConf.getpSecondSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
    return raisedConf;
}






NonManifoldCorner AugmentedBetaComplex::resolve_FF_E_Configuration(const NonManifoldCorner& givenNMConf, 
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges )
{
    AugmentedBCEdge* nmEdge     = (AugmentedBCEdge*)( givenNMConf.getNMSimplex().getSimplex() );
    AugmentedBCFace* givenFace1 = (AugmentedBCFace*)( givenNMConf.getFirstSimplex().getSimplex() );
    AugmentedBCFace* givenFace2 = (AugmentedBCFace*)( givenNMConf.getSecondSimplex().getSimplex() );

    AugmentedBCFace* raisedNMFace = rg_NULL;
    AugmentedBCCell* raisedTetra1 = rg_NULL;
    AugmentedBCCell* raisedTetra2 = rg_NULL;

    


    if ( givenFace1->getEdgeOrientation( nmEdge ) != rg_TRUE ) {
        givenFace1->reverse();
    }

    if ( givenFace2->getEdgeOrientation( nmEdge ) != rg_FALSE ) {
        givenFace2->reverse();
    }
    
    rg_INT indexOfCommonExteriorRegion = givenFace1->getIndexOfExteriorRegionToRightCell(); 



    O_FF_E_operator( nmEdge, givenFace1, givenFace2, raisedNMFace, raisedTetra1, raisedTetra2 );

    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenFace1, raisedTetra1, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );
    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenFace2, raisedTetra2, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );    


    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, givenFace1, raisedTetra1 );
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, givenFace2, raisedTetra2 );


    NonManifoldCorner raisedConf( raisedNMFace, givenFace1, givenFace2 );
	raisedConf.getpFirstSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
	raisedConf.getpSecondSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
    return raisedConf;
}



NonManifoldCorner AugmentedBetaComplex::resolve_FT_E_Configuration(const NonManifoldCorner& givenNMConf, 
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges)
{
    AugmentedBCEdge* nmEdge			  = (AugmentedBCEdge*)( givenNMConf.getNMSimplex().getSimplex() );
    AugmentedBCFace* givenFace		  = (AugmentedBCFace*)( givenNMConf.getFirstSimplex().getSimplex() );
    AugmentedBCFace* givenRegularFace = (AugmentedBCFace*)( givenNMConf.getSecondSimplex().getSimplex() );

    
    AugmentedBCFace* raisedNMFace = rg_NULL;
    AugmentedBCCell* raisedTetra1 = rg_NULL;
    AugmentedBCCell* raisedTetra2 = rg_NULL;

    
    AugmentedBCFace* regularFaceBoudedByGivenEdge = findRegularFaceBoundedByGivenEdge( nmEdge, givenRegularFace );	
    AugmentedBCFace* validRegularFace = findRegularFaceOfGivenCellByBackwardOrientation( regularFaceBoudedByGivenEdge, nmEdge );

    if (    givenFace->getEdgeOrientation( nmEdge) 
         == validRegularFace->getEdgeOrientation( nmEdge) ) {
        givenFace->reverse();   //givenRegularFace에 방향을 맞춘다.
    }

    rg_INT indexOfCommonExteriorRegion = validRegularFace->getIndexOfExteriorRegionToRightCell();


    O_FT_E_operator( nmEdge, givenFace, validRegularFace, raisedNMFace, raisedTetra1, raisedTetra2 );


    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenFace, raisedTetra1, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );        


    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, givenFace, raisedTetra1 );
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, validRegularFace, raisedTetra2 );


    NonManifoldCorner raisedConf( raisedNMFace, givenFace, validRegularFace );
	raisedConf.getpFirstSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
	raisedConf.getpSecondSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
    return raisedConf;
}



NonManifoldCorner AugmentedBetaComplex::resolve_TF_E_Configuration(const NonManifoldCorner& givenNMConf, 
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                                                rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges)
{
    AugmentedBCEdge* nmEdge			  = (AugmentedBCEdge*)( givenNMConf.getNMSimplex().getSimplex() );
    AugmentedBCFace* givenRegularFace = (AugmentedBCFace*)( givenNMConf.getFirstSimplex().getSimplex() );
    AugmentedBCFace* givenFace		  = (AugmentedBCFace*)( givenNMConf.getSecondSimplex().getSimplex() );    

    AugmentedBCFace* raisedNMFace        = rg_NULL;
    AugmentedBCCell* raisedTetra1 = rg_NULL;
    AugmentedBCCell* raisedTetra2 = rg_NULL;

    
    
    AugmentedBCFace* regularFaceBoudedByGivenEdge = findRegularFaceBoundedByGivenEdge( nmEdge, givenRegularFace );	
    AugmentedBCFace* validRegularFace = findRegularFaceOfGivenCellByForwardOrientation( regularFaceBoudedByGivenEdge, nmEdge );

    if (    givenFace->getEdgeOrientation( nmEdge )
         == validRegularFace->getEdgeOrientation( nmEdge) ) {
        givenFace->reverse();
    }


    rg_INT indexOfCommonExteriorRegion = validRegularFace->getIndexOfExteriorRegionToRightCell();


    
    

    O_TF_E_operator( nmEdge, validRegularFace, givenFace, raisedNMFace, raisedTetra1, raisedTetra2 );

    updateNMCornersForLowerDimensionSimplexOfGivenSimplex( givenFace, raisedTetra2, ptToNMCornerOfNMVertices, ptToNMCornerOfNMEdges );


    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, validRegularFace, raisedTetra1 );
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, givenFace, raisedTetra2 );


    NonManifoldCorner raisedConf( raisedNMFace, validRegularFace, givenFace );
	raisedConf.getpFirstSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
	raisedConf.getpSecondSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
    return raisedConf;
}


NonManifoldCorner AugmentedBetaComplex::resolve_TT_E_Configuration( const NonManifoldCorner& givenNMConf )
{
    AugmentedBCEdge* nmEdge			   = (AugmentedBCEdge*)( givenNMConf.getNMSimplex().getSimplex() );
    AugmentedBCFace* givenRegularFace1 = (AugmentedBCFace*)( givenNMConf.getFirstSimplex().getSimplex() );
    AugmentedBCFace* givenRegularFace2 = (AugmentedBCFace*)( givenNMConf.getSecondSimplex().getSimplex() );

    AugmentedBCFace* raisedNMFace = rg_NULL;
    AugmentedBCCell* raisedTetra1 = rg_NULL;
    AugmentedBCCell* raisedTetra2 = rg_NULL;
    
    
    AugmentedBCFace* regularFaceBoudedByGivenEdge1 = findRegularFaceBoundedByGivenEdge( nmEdge, givenRegularFace1 );
	AugmentedBCFace* regularFaceOfGivenTet1 = findRegularFaceOfGivenCellByForwardOrientation( regularFaceBoudedByGivenEdge1, nmEdge );

	AugmentedBCFace* regularFaceBoudedByGivenEdge2 = findRegularFaceBoundedByGivenEdge( nmEdge, givenRegularFace2 );	
    AugmentedBCFace* regularFaceOfGivenTet2 = findRegularFaceOfGivenCellByBackwardOrientation( regularFaceBoudedByGivenEdge2, nmEdge );

    
    rg_INT indexOfCommonExteriorRegion = regularFaceOfGivenTet1->getIndexOfExteriorRegionToRightCell();

        
    O_TT_E_operator( nmEdge, regularFaceOfGivenTet1, regularFaceOfGivenTet2, raisedNMFace, raisedTetra1, raisedTetra2 );


    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, regularFaceOfGivenTet1, raisedTetra1 );
    setIndexOfExteriorRegionOfCellExceptGivenFace( indexOfCommonExteriorRegion, regularFaceOfGivenTet2, raisedTetra2 );


    NonManifoldCorner raisedConf( raisedNMFace, regularFaceOfGivenTet1, regularFaceOfGivenTet2 );
	raisedConf.getpFirstSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
	raisedConf.getpSecondSimplex()->setTypeOfSimplex( CELL_SIMPLEX );
    return raisedConf;
}


NonManifoldCorner AugmentedBetaComplex::resolve_TT_D_V_Configuration( const NonManifoldCorner& givenNMConf )
{
    AugmentedBCVertex* nmVertex = (AugmentedBCVertex*)( givenNMConf.getNMSimplex().getSimplex() );

	O_TT_D_V_operator( nmVertex );	
    
    NonManifoldCorner raisedConf( nmVertex );
    return raisedConf;
}


void AugmentedBetaComplex::O_TT_D_V_operator( AugmentedBCVertex* nmVertex )
{
	rg_dList< rg_dList<AugmentedBCFace*> > faceGroupList;
    groupFacesIncidentToVertexByEdge( nmVertex, faceGroupList);

    if ( faceGroupList.getSize() <=1 )
    {                                   //TT_T conf의 detect가 정확하다면 return되면 안된다.
        NonManifoldCorner raisedConf( nmVertex );
        return;
    }        

    rg_INT numOfFaceGroups = faceGroupList.getSize();
    rg_INT i=0;
    rg_dList<AugmentedBCFace*>* currFaceGroup = rg_NULL;
    faceGroupList.reset4Loop();
    while ( faceGroupList.setNext4Loop() && i < numOfFaceGroups - 1 )
    {
        currFaceGroup = faceGroupList.getpEntity();

        addArtificialCellsByGivenVertexAndFaceGroup( nmVertex, currFaceGroup );

        i++;
    }
}




rg_INT AugmentedBetaComplex::findIndexOfCommonExteriorRegion( AugmentedBCFace* givenFace1, AugmentedBCFace* givenFace2 )
{
    rg_INT leftExtOfFace1  = givenFace1->getIndexOfExteriorRegionToLeftCell();
    rg_INT rightExtOfFace1 = givenFace1->getIndexOfExteriorRegionToRightCell();

    rg_INT leftExtOfFace2  = givenFace2->getIndexOfExteriorRegionToLeftCell();
    rg_INT rightExtOfFace2 = givenFace2->getIndexOfExteriorRegionToRightCell();

    if ( leftExtOfFace1 != -1 ) {
        if ( leftExtOfFace1 == leftExtOfFace2 ) {
            return leftExtOfFace1;
        }
        else if ( leftExtOfFace1 == rightExtOfFace2 ) {
            return leftExtOfFace1;
        }
    }

    
    if ( rightExtOfFace1 != -1 ) {
        if ( rightExtOfFace1 == leftExtOfFace2 ) {
            return rightExtOfFace1;
        }
        else if ( rightExtOfFace1 == rightExtOfFace2 ) {
            return rightExtOfFace1;
        }
    }

    
    return -1;
}




void AugmentedBetaComplex::setIndexOfExteriorRegionOfSingularFace( const rg_INT& indexOfExteriorRegion, AugmentedBCFace* singularFace )
{
    singularFace->setIndexOfExteriorRegionToLeftCell( indexOfExteriorRegion );
    singularFace->setIndexOfExteriorRegionToRightCell( indexOfExteriorRegion );
}



void AugmentedBetaComplex::setIndexOfExteriorRegionOfCellExceptGivenFace( const rg_INT& indexOfExteriorRegion, AugmentedBCFace* exceptThisFace, AugmentedBCCell* givenCell )
{
    BetaFace** boundingFaceOfGivenCell = givenCell->getFaces();
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++) {
        if ( boundingFaceOfGivenCell[i] == exceptThisFace ) {
            continue;
        }
        
        //((AugmentedBCFace*) boundingFaceOfGivenCell[i])->setIndexOfExteriorRegionToLeftCell( indexOfExteriorRegion );
        ((AugmentedBCFace*) boundingFaceOfGivenCell[i])->setIndexOfExteriorRegionToRightCell( indexOfExteriorRegion );
    }    
}


void AugmentedBetaComplex::checkComparingNumExtraneousCellGroupAndNumRegularFaceGroup()
{
	r_numVerticesDifferentNumConesInputAndTT_D_V = 0;

	ofstream fout("numOfExtraneousAndRegular.txt");

	m_augmentedBCVertices.reset4Loop();
	while ( m_augmentedBCVertices.setNext4Loop() )
	{
		AugmentedBCVertex* currVertex = m_augmentedBCVertices.getpEntity();

		if ( currVertex->getBoundingState() != REGULAR_SIMPLEX )
			continue;

        rg_INT numOfExteriorSpace = 0;
        if ( currVertex->getFirstEdge() == rg_NULL ) {
            numOfExteriorSpace = 1;
        }
        else if ( ((AugmentedBCEdge* )currVertex->getFirstEdge())->getBoundingState() == SINGULAR_SIMPLEX ) {
            numOfExteriorSpace = 1;
        }
        else if ( currVertex->getFirstEdge()->getFirstFace() == rg_NULL ) {
            numOfExteriorSpace = 1;
        }
        else if ( ((AugmentedBCFace*) currVertex->getFirstEdge()->getFirstFace())->getBoundingState()
                       == SINGULAR_SIMPLEX ) {
            numOfExteriorSpace = 1;
        }
        else {
            rg_dList< rg_dList<AugmentedBCFace*> > faceGroupList;
		    groupFacesIncidentToVertexByEdge( currVertex, faceGroupList);

            numOfExteriorSpace = faceGroupList.getSize();
        }

		

		if ( r_numExtraneousCellGroupAtVertex[currVertex->getID()] > 1
			|| numOfExteriorSpace > 1 )
			fout << "*";

		if ( r_numExtraneousCellGroupAtVertex[currVertex->getID()] != numOfExteriorSpace )
		{
			r_numVerticesDifferentNumConesInputAndTT_D_V++;
			fout << "!!";	
		}
		
			
		fout << "Vertex ID: " << currVertex->getID() << "\t" << "numOfExtraneousCellGroup: " << r_numExtraneousCellGroupAtVertex[currVertex->getID()] << "\t"
			 << "numOfRegularFaceGroup: " << numOfExteriorSpace << endl;

	}
	fout.close();
}




void AugmentedBetaComplex::reportNumNMCorners(rg_dList< NonManifoldCorner >& nmCorners)
{
	ofstream fout("numOfNMCorners.txt");

	rg_INT numOfNMCorner[11];
	for (rg_INT i=0;i<11;i++)
		numOfNMCorner[i]=0;

	nmCorners.reset4Loop();
	while ( nmCorners.setNext4Loop() )
	{
		NonManifoldCorner* currCorner = nmCorners.getpEntity();
		rg_INT			   currConf   = currCorner->getTypeOfNonManifoldConfiguration();

		if ( currConf == TT_D_V_CONFIGURATION)
			numOfNMCorner[TT_D_V_CONFIGURATION] += currCorner->getpFirstSimplex()->getTypeOfSimplex();
		else
			numOfNMCorner[ currConf ]++;
	}

	fout << "EE_V" << '\t' << numOfNMCorner[EE_V_CONFIGURATION] << endl;
	fout << "EF_V" << '\t' << numOfNMCorner[EF_V_CONFIGURATION] << endl;
	fout << "ET_V" << '\t' << numOfNMCorner[ET_V_CONFIGURATION] << endl;
	fout << "FF_V" << '\t' << numOfNMCorner[FF_V_CONFIGURATION] << endl;
	fout << "FT_V" << '\t' << numOfNMCorner[FT_V_CONFIGURATION] << endl;
	fout << "TT_V" << '\t' << numOfNMCorner[TT_V_CONFIGURATION] << endl;
	fout << "FF_E" << '\t' << numOfNMCorner[FF_E_CONFIGURATION] << endl;
	fout << "FT_E" << '\t' << numOfNMCorner[FT_E_CONFIGURATION]+numOfNMCorner[TF_E_CONFIGURATION] << endl;
	fout << "TT_E" << '\t' << numOfNMCorner[TT_E_CONFIGURATION] << endl;
	fout << "TT_D_V" << '\t' << numOfNMCorner[TT_D_V_CONFIGURATION] << endl;

	fout.close();
}



void AugmentedBetaComplex::updateNMCornersForLowerDimensionSimplexOfGivenSimplex(AugmentedBCEdge* givenEdge, AugmentedBCFace* raisedFace, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices)
{
    AugmentedBCVertex* lowerDimensionSimplex[2] = {(AugmentedBCVertex*)givenEdge->getStartVertex(), (AugmentedBCVertex*)givenEdge->getEndVertex()};

    rg_dList<NonManifoldCorner*>* currPTToNMCornerOfNMVertices = rg_NULL;

    rg_INT i;
    for (i=0; i<2; i++)
    {
        if ( lowerDimensionSimplex[i]->isArtificial() == rg_TRUE )
            continue;

        currPTToNMCornerOfNMVertices = &( ptToNMCornerOfNMVertices[ lowerDimensionSimplex[i]->getID() ] );

        replaceGivenSimplexToRaisedSimplex( givenEdge, raisedFace, currPTToNMCornerOfNMVertices );
    }
}



void AugmentedBetaComplex::updateNMCornersForLowerDimensionSimplexOfGivenSimplex(AugmentedBCFace* givenFace, AugmentedBCCell* raisedTetra, 
                                                                                           rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges)
{
	rg_dList<BetaVertex*> boundingVertices;
    givenFace->searchVerticesInWholeWorld( boundingVertices );    

    rg_dList<NonManifoldCorner*>* currPTToNMCornerOfNMVertices = rg_NULL;

	boundingVertices.reset4Loop();
	while ( boundingVertices.setNext4Loop() )
	{
		AugmentedBCVertex* currAugVertex = (AugmentedBCVertex*)boundingVertices.getEntity();
        if ( currAugVertex->isArtificial() == rg_TRUE )
            continue;

        currPTToNMCornerOfNMVertices = &( ptToNMCornerOfNMVertices[ currAugVertex->getID() ] );

        replaceGivenSimplexToRaisedSimplex( givenFace, raisedTetra, currPTToNMCornerOfNMVertices );
    }

	BetaEdge** boundingEdge = givenFace->getEdges();

	rg_dList<NonManifoldCorner*>* currPTToNMCornerOfNMEdges = rg_NULL;
	rg_INT i;
    for (i=0; i<EIWDS_NUM_EDGE_ON_FACE; i++)
    {
		AugmentedBCEdge* currAugEdge = (AugmentedBCEdge*)boundingEdge[i];
        if ( currAugEdge->isArtificial() == rg_TRUE )
            continue;

        currPTToNMCornerOfNMEdges = &( ptToNMCornerOfNMEdges[ currAugEdge->getID() ] );

        replaceGivenSimplexToRaisedSimplex( givenFace, raisedTetra, currPTToNMCornerOfNMEdges);
    }
}


void AugmentedBetaComplex::replaceGivenSimplexToRaisedSimplex( void* original, void* substitute, rg_dList<NonManifoldCorner*>* currPTToNMCornerOfNMSimplex)
{
    currPTToNMCornerOfNMSimplex->reset4Loop();
    while ( currPTToNMCornerOfNMSimplex->setNext4Loop() )
    {
        NonManifoldCorner* currNMConfiguration = currPTToNMCornerOfNMSimplex->getEntity();

        PointerToSimplex* firstSimplexOfCurrNMConf  = currNMConfiguration->getpFirstSimplex();
        PointerToSimplex* secondSimplexOfCurrNMConf = currNMConfiguration->getpSecondSimplex();

        if ( firstSimplexOfCurrNMConf->getSimplex() == original )
        {
            rg_INT originalTypeOfSimplex = firstSimplexOfCurrNMConf->getTypeOfSimplex();
            if ( originalTypeOfSimplex == EDGE_SIMPLEX )
			{
                firstSimplexOfCurrNMConf->setTypeOfSimplex( FACE_SIMPLEX );
				firstSimplexOfCurrNMConf->setSimplex( substitute );
			}
            else if ( originalTypeOfSimplex == FACE_SIMPLEX )
			{
                firstSimplexOfCurrNMConf->setTypeOfSimplex( CELL_SIMPLEX );
			}
        }
        else if ( secondSimplexOfCurrNMConf->getSimplex() == original )
        {
            rg_INT originalTypeOfSimplex = secondSimplexOfCurrNMConf->getTypeOfSimplex();
            if ( originalTypeOfSimplex == EDGE_SIMPLEX )
			{
                secondSimplexOfCurrNMConf->setTypeOfSimplex( FACE_SIMPLEX );
				secondSimplexOfCurrNMConf->setSimplex( substitute );
			}
            else if ( originalTypeOfSimplex == FACE_SIMPLEX )
			{
                secondSimplexOfCurrNMConf->setTypeOfSimplex( CELL_SIMPLEX );
			}
        }
    }
}



void AugmentedBetaComplex::O_EE_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCEdge* edgeSimplex1, AugmentedBCEdge* edgeSimplex2, AugmentedBCEdge*& raisedNMSimplex, AugmentedBCFace*& raisedSimplex1, AugmentedBCFace*& raisedSimplex2)
{   
    AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex);
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCFace* artificialFace1 = addArtificialFaceBy2Edges( artificialEdge, edgeSimplex1, givenVertex );
    AugmentedBCFace* artificialFace2 = addArtificialFaceBy2Edges( edgeSimplex2, artificialEdge, givenVertex );    
    
    artificialFace1->setBoundingState( SINGULAR_SIMPLEX );
    artificialFace2->setBoundingState( SINGULAR_SIMPLEX );

    raisedNMSimplex = artificialEdge;
    raisedSimplex1 = artificialFace1;
    raisedSimplex2 = artificialFace2;
}



void AugmentedBetaComplex::O_EF_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCEdge* givenEdge, AugmentedBCFace* givenFace, AugmentedBCEdge*& raisedNMSimplex, AugmentedBCFace*& raisedSimplex1, AugmentedBCCell*& raisedSimplex2)
{
    AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCFace* artificialFace = addArtificialFaceBy2Edges( artificialEdge, givenEdge, givenVertex );
    artificialFace->setBoundingState( SINGULAR_SIMPLEX );

    AugmentedBCCell* artificialCell = addArtificialCellBy1EdgeAnd1Face( artificialEdge, givenFace, givenVertex );

    raisedNMSimplex = artificialEdge;
    raisedSimplex1  = artificialFace;
    raisedSimplex2  = artificialCell;
}


void AugmentedBetaComplex::O_ET_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCEdge* givenEdge, AugmentedBCFace* givenRegularFace, AugmentedBCEdge*& raisedNMEdge, AugmentedBCFace*& raisedFace, AugmentedBCCell*& raisedTetra)
{
    O_EF_V_operator( givenVertex, givenEdge, givenRegularFace, raisedNMEdge, raisedFace, raisedTetra );

    /*
    AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCFace* artificialFace = addArtificialFaceBy2Edges( artificialEdge, givenEdge, givenVertex );
    artificialFace->setBoundingState( SINGULAR_SIMPLEX );

    AugmentedBCCell* artificialCell = addArtificialCellBy1EdgeAnd1Face( artificialEdge, givenRegularFace, givenVertex );
    
    raisedNMEdge = artificialEdge;
    raisedFace   = artificialFace;
    raisedTetra  = artificialCell;    
    */
}


AugmentedBCFace* AugmentedBetaComplex::findRegularFaceTouchedExactExterior( AugmentedBCVertex* givenVertex, 
																		    AugmentedBCFace* givenFace,
                                                                            const rg_INT& indexOfExteriorRegion )
{
	if ( givenFace->getBoundingState() == REGULAR_SIMPLEX
		&& givenFace->isThere( givenVertex )
        && givenFace->getIndexOfExteriorRegionToRightCell() == indexOfExteriorRegion )
	{	
		return givenFace;
	}

	BetaCell* firstCell = givenFace->getLeftCell();

    rg_dList< BetaCell* > visitedCellList;

    rg_dList< BetaCell* > cellStack;
    cellStack.pushBack( firstCell );

    while ( cellStack.getSize() > 0 ) {
        BetaCell* popedCell = cellStack.popBack();

        if ( popedCell->isVisited() ) {
            continue;
        }

        popedCell->isVisited( rg_TRUE );
        visitedCellList.add( popedCell );

        BetaFace** boudingFace = popedCell->getFaces();

        for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ ) {
            if ( boudingFace[i]->isThere( givenVertex) == rg_FALSE ) {
                continue;
            }


            if (    ((AugmentedBCFace*) boudingFace[i])->getBoundingState() == REGULAR_SIMPLEX		        
                 && ((AugmentedBCFace*) boudingFace[i])->getIndexOfExteriorRegionToRightCell() == indexOfExteriorRegion )
	        {	
                visitedCellList.reset4Loop();
                while ( visitedCellList.setNext4Loop() ) {
		            visitedCellList.getEntity()->isVisited( rg_FALSE );
                }

		        return (AugmentedBCFace*)boudingFace[i];
	        }


            BetaCell* nextCell = rg_NULL;
            if ( boudingFace[i]->getLeftCell() == popedCell ) {
                nextCell = boudingFace[i]->getRightCell();
            }
            else {
                nextCell = boudingFace[i]->getLeftCell();
            }

            if ( nextCell == rg_NULL ) {
                continue;
            }

            if ( nextCell->isVisited() == rg_FALSE ) {
                cellStack.pushBack( nextCell );
            }
        }
    }


	

	visitedCellList.reset4Loop();
    while ( visitedCellList.setNext4Loop() ) {
		visitedCellList.getEntity()->isVisited( rg_FALSE );
    }

	return rg_NULL;
}




void AugmentedBetaComplex::O_FF_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCFace* givenFirstFace, AugmentedBCFace* givenSecondFace, AugmentedBCEdge*& raisedNMEdge, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2)
{
    AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCCell* artificialCell1 = addArtificialCellBy1EdgeAnd1Face( artificialEdge, givenFirstFace, givenVertex );
    AugmentedBCCell* artificialCell2 = addArtificialCellBy1EdgeAnd1Face( artificialEdge, givenSecondFace, givenVertex );

    raisedNMEdge = artificialEdge;
    raisedTetra1 = artificialCell1;
    raisedTetra2 = artificialCell2;
}


void AugmentedBetaComplex::O_FT_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCFace* givenFace, AugmentedBCFace* givenRegularFace, AugmentedBCEdge*& raisedNMEdge, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2)
{
    O_FF_V_operator( givenVertex, givenFace, givenRegularFace, raisedNMEdge, raisedTetra1, raisedTetra2 );


    /*
	AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCCell* artificialCell1 = addArtificialCellBy1EdgeAnd1Face( artificialEdge, givenFace, givenVertex );
    AugmentedBCCell* artificialCell2 = addArtificialCellBy1EdgeAnd1Face( artificialEdge, givenRegularFace, givenVertex );

    raisedNMEdge = artificialEdge;
    raisedTetra1 = artificialCell1;
    raisedTetra2 = artificialCell2;
    */
}


void AugmentedBetaComplex::O_TT_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCFace* givenFirstRegularFace, AugmentedBCFace* givenSecondRegularFace, AugmentedBCEdge*& raisedNMEdge, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2)
{
    O_FF_V_operator( givenVertex, givenFirstRegularFace, givenSecondRegularFace, raisedNMEdge, raisedTetra1, raisedTetra2 );


    /*
	AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCCell* artificialCell1 = addArtificialCellBy1EdgeAnd1Face( artificialEdge, givenFirstRegularFace, givenVertex );
    AugmentedBCCell* artificialCell2 = addArtificialCellBy1EdgeAnd1Face( artificialEdge, givenSecondRegularFace, givenVertex );    
    
    raisedNMEdge = artificialEdge;
    raisedTetra1 = artificialCell1;
    raisedTetra2 = artificialCell2;
    */
}



void AugmentedBetaComplex::O_FF_E_operator(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenFirstFace, AugmentedBCFace* givenSecondFace, AugmentedBCFace*& raisedNMFace, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2)
{
    AugmentedBCVertex* givenVertex = (AugmentedBCVertex*)givenEdge->getStartVertex();  //임의로 start vertex로 하였다.(end vertex로 하여도 상관 없음)
    AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCFace* artificialFace = addArtificialFaceBy2Edges( givenEdge, artificialEdge, givenVertex );
    AugmentedBCCell* artificialTet1 = addArtificialCellByRegularAndSingularFaces( givenFirstFace, artificialFace,  givenEdge);
    AugmentedBCCell* artificialTet2 = addArtificialCellByRegularAndSingularFaces( givenSecondFace, artificialFace, givenEdge);
    
    artificialFace->setRightCell( artificialTet1 );

    raisedNMFace = artificialFace;
    raisedTetra1 = artificialTet1;
    raisedTetra2 = artificialTet2;
}



void AugmentedBetaComplex::O_FT_E_operator(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenFace, AugmentedBCFace* givenRegularFace, AugmentedBCFace*& raisedNMFace, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2)
{
    O_FF_E_operator( givenEdge, givenFace, givenRegularFace, raisedNMFace, raisedTetra1, raisedTetra2 );


    /*
	AugmentedBCVertex* givenVertex = (AugmentedBCVertex*)givenEdge->getStartVertex();  //임의로 start vertex로 하였다.(end vertex로 하여도 상관 없음)
    AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCFace* artificialFace = addArtificialFaceBy2Edges( givenEdge, artificialEdge, givenVertex );
    AugmentedBCCell* artificialTet1 = addArtificialCellByRegularAndSingularFaces( givenFace, artificialFace,  givenEdge);
    AugmentedBCCell* artificialTet2 = addArtificialCellByRegularAndSingularFaces( givenRegularFace, artificialFace, givenEdge);
    
    artificialFace->setRightCell( artificialTet1 );
    
    raisedNMFace = artificialFace;
    raisedTetra1 = artificialTet1;
    raisedTetra2 = artificialTet2;
    */
}


void AugmentedBetaComplex::O_TF_E_operator(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenRegularFace, AugmentedBCFace* givenFace, AugmentedBCFace*& raisedNMFace, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2)
{
    O_FF_E_operator( givenEdge, givenRegularFace, givenFace, raisedNMFace, raisedTetra1, raisedTetra2 );



    /*
    AugmentedBCVertex* givenVertex = (AugmentedBCVertex*)givenEdge->getStartVertex();  //임의로 start vertex로 하였다.(end vertex로 하여도 상관 없음)
    AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCFace*        artificialFace = addArtificialFaceBy2Edges( givenEdge, artificialEdge, givenVertex );
    AugmentedBCCell* artificialTet1 = addArtificialCellByRegularAndSingularFaces( givenRegularFace, artificialFace, givenEdge);
    AugmentedBCCell* artificialTet2 = addArtificialCellByRegularAndSingularFaces( givenFace, artificialFace,  givenEdge);
    
    artificialFace->setRightCell( artificialTet1 );

    raisedNMFace = artificialFace;
    raisedTetra1 = artificialTet1;
    raisedTetra2 = artificialTet2;
    */
}




AugmentedBCFace* AugmentedBetaComplex::findRegularFaceBoundedByGivenEdge( AugmentedBCEdge* givenEdge, AugmentedBCFace* givenRegularFace )
{
	if ( givenRegularFace->isThere( givenEdge ) 
		&& givenRegularFace->getBoundingState() == REGULAR_SIMPLEX )
		return givenRegularFace;

	rg_dList< BetaCell*> visitedCellList;
	rg_dList< AugmentedBCFace* > facesToBeTraced;

	BetaCell* cellIncidentToGivenFace[2] = { givenRegularFace->getLeftCell(), givenRegularFace->getRightCell() };
	for ( rg_INT i=0; i<2; i++ )	{
		if ( cellIncidentToGivenFace[i] == rg_NULL )
			continue;

		cellIncidentToGivenFace[i]->isVisited(rg_TRUE);
		visitedCellList.add( cellIncidentToGivenFace[i] );

		BetaFace** boundingFace = cellIncidentToGivenFace[i]->getFaces();
		for ( rg_INT face_j=0; face_j<EIWDS_NUM_FACE_ON_CELL; face_j++)	{
			if ( boundingFace[face_j]->isThere( givenEdge )
				&& ((AugmentedBCFace*)boundingFace[face_j])->getBoundingState() == REGULAR_SIMPLEX )
			{
				visitedCellList.reset4Loop();
				while ( visitedCellList.setNext4Loop() )
					visitedCellList.getEntity()->isVisited( rg_FALSE );
				
				return (AugmentedBCFace*)boundingFace[face_j];
			}

			if ( boundingFace[face_j]->isThere( givenEdge->getStartVertex() )
			  || boundingFace[face_j]->isThere( givenEdge->getEndVertex() ) )
			{
				if ( ((AugmentedBCFace*)boundingFace[face_j])->getBoundingState() == INTERIOR_SIMPLEX )
				{
					facesToBeTraced.add( (AugmentedBCFace*)boundingFace[face_j] );
				}
			}
		}
	}

	while ( facesToBeTraced.getSize() > 0 )
	{
		AugmentedBCFace* currFace = facesToBeTraced.popFront();

		BetaCell* nextCell = rg_NULL;
		if ( currFace->getLeftCell()->isVisited() )
			nextCell = currFace->getRightCell();
		else
			nextCell = currFace->getLeftCell();

		if ( nextCell == rg_NULL || nextCell->isVisited() ) //nextCell == rg_NULL 조건은 필요없다.
			continue;

		nextCell->isVisited( rg_TRUE );
		visitedCellList.add( nextCell );

		BetaFace** boundingFaceOfNextCell = nextCell->getFaces();
		for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )
		{
			AugmentedBCFace* currFace = (AugmentedBCFace*)boundingFaceOfNextCell[i];
			if ( currFace->isThere( givenEdge )
				&& currFace->getBoundingState() == REGULAR_SIMPLEX )
			{
				visitedCellList.reset4Loop();
				while ( visitedCellList.setNext4Loop() )
					visitedCellList.getEntity()->isVisited( rg_FALSE );
				
				return currFace;
			}

			if ( currFace->isThere( givenEdge->getStartVertex() )
			  || currFace->isThere( givenEdge->getEndVertex() ) )
			{
				if ( currFace->getBoundingState() == INTERIOR_SIMPLEX )
				{
					facesToBeTraced.add( currFace );
				}
			}
		}

	}

	return rg_NULL;
}








void AugmentedBetaComplex::O_TT_E_operator(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenRegularFace1, AugmentedBCFace* givenRegularFace2, AugmentedBCFace*& raisedNMFace, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2)
{
    O_FF_E_operator( givenEdge, givenRegularFace1, givenRegularFace2, raisedNMFace, raisedTetra1, raisedTetra2 );
    

    /*
    AugmentedBCVertex* givenVertex = (AugmentedBCVertex*)givenEdge->getStartVertex();  //임의로 start vertex로 하였다.(end vertex로 하여도 상관 없음)
    AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCFace*        artificialFace = addArtificialFaceBy2Edges( givenEdge, artificialEdge, givenVertex );
    AugmentedBCCell* artificialTet1 = addArtificialCellByRegularAndSingularFaces( givenRegularFace1, artificialFace, givenEdge);
    AugmentedBCCell* artificialTet2 = addArtificialCellByRegularAndSingularFaces( givenRegularFace2, artificialFace, givenEdge);     

    artificialFace->setRightCell( artificialTet1 );

    raisedNMFace = artificialFace;
    raisedTetra1 = artificialTet1;
    raisedTetra2 = artificialTet2;
    */
}





void AugmentedBetaComplex::addArtificialCellsByGivenVertexAndFaceGroup(AugmentedBCVertex* givenVertex, rg_dList<AugmentedBCFace*>* faceGroup)
{
    if ( faceGroup->getSize() < 3) //왜 2가 아닌가? ->regular face가 3개는 있어야 loop을 이룬다.
        return;

    AugmentedBCVertex* artificialVertex = replicateGivenVertex( givenVertex );
    AugmentedBCEdge*   artificialEdge   = addArtificialEdge( givenVertex, artificialVertex );
    artificialVertex->setFirstEdge( artificialEdge );

    AugmentedBCFace* firstRegularFace = faceGroup->getFirstEntity();
    AugmentedBCCell* firstArtificialTetr = addArtificialCellBy1EdgeAnd1Face( artificialEdge, firstRegularFace, givenVertex );

    AugmentedBCFace* secondRegularFace = faceGroup->getSecondEntity();
    AugmentedBCFace* lastRegularFace   = faceGroup->getLastEntity();

    AugmentedBCEdge* currCommonEdge = rg_NULL;
    AugmentedBCFace* currBridgeFace = rg_NULL;
    findBridgeFaceByIncidentFaceAndCellNotThisBoundingFace(givenVertex, secondRegularFace, firstArtificialTetr, firstRegularFace,
                                                                    currCommonEdge, currBridgeFace);

    AugmentedBCEdge* lastCommonEdge = rg_NULL;
    AugmentedBCFace* lastBridgeFace = rg_NULL;
    findBridgeFaceByIncidentFaceAndCellNotThisBoundingFace(givenVertex, lastRegularFace, firstArtificialTetr, firstRegularFace,
                                                                    lastCommonEdge, lastBridgeFace);

    AugmentedBCCell* prevCell = firstArtificialTetr;
    AugmentedBCFace* currRegaulrFace = rg_NULL;
    AugmentedBCFace* nextRegaulrFace = rg_NULL;
    AugmentedBCCell* currArtificialTetr = rg_NULL;
    faceGroup->reset4Loop();
    faceGroup->setNext4Loop(); //두번째 regular face부터 시작하기 위해
    while (faceGroup->setNext4Loop() )
    {
        currRegaulrFace = faceGroup->getEntity();

        if ( currRegaulrFace == lastRegularFace )
            break;

        currArtificialTetr = addArtificialCellByRegularAndSingularFaces( currRegaulrFace, currBridgeFace, currCommonEdge );
        currBridgeFace->setRightCell( prevCell );

        prevCell = currArtificialTetr;

        nextRegaulrFace = faceGroup->getNextEntity();

        findBridgeFaceByIncidentFaceAndCellNotThisBoundingFace(givenVertex, nextRegaulrFace, currArtificialTetr, currRegaulrFace,
                                                                     currCommonEdge, currBridgeFace );
    }

    AugmentedBCEdge* e1 = lastRegularFace->findMateEdgeOfVertex( givenVertex );
    AugmentedBCEdge* e2 = lastBridgeFace->findMateEdgeOfVertex( givenVertex );
    AugmentedBCEdge* e3 = currBridgeFace->findMateEdgeOfVertex(givenVertex );

    AugmentedBCFace* lastArtificialFace = addArtificialFaceBy3Edges(e1, e2, e3);

    
    AugmentedBCCell tempCell(currRegaulrFace, currBridgeFace, lastBridgeFace, lastArtificialFace);
    tempCell.setIsArtificial( rg_TRUE );
    tempCell.setBoundingState( INTERIOR_SIMPLEX );
    
    tempCell.setID( getMaxIndexOfAugmentedBCCells() + 1 );
    

    AugmentedBCCell* newCell = m_artificialAugmentedBCCells.add( tempCell );

    currRegaulrFace->setRightCell( newCell);
    currBridgeFace->setRightCell( newCell );
    lastBridgeFace->setRightCell( newCell );
    lastArtificialFace->setLeftCell( newCell );

    currRegaulrFace->setBoundingState( INTERIOR_SIMPLEX );
    currBridgeFace->setBoundingState( INTERIOR_SIMPLEX );
    lastBridgeFace->setBoundingState( INTERIOR_SIMPLEX );
    lastArtificialFace->setBoundingState( REGULAR_SIMPLEX );

    artificialEdge->setBoundingState( INTERIOR_SIMPLEX );
    
    setClassOfAugmentedBCEdgesToInteriorIfEdgesBoundRegularFaces(givenVertex, faceGroup);
}



void AugmentedBetaComplex::findBridgeFaceByIncidentFaceAndCellNotThisBoundingFace(AugmentedBCVertex* givenVertex, AugmentedBCFace* incidentFace, AugmentedBCCell* givenCell, AugmentedBCFace* noThisBoundingFace,
                                                                    AugmentedBCEdge*& outputCommonEdge, AugmentedBCFace*& outputBoundingFace)
{
    for ( rg_INT i=0; i<3; i++)    
    {
        BetaEdge* currEdge = incidentFace->getEdge(i);
        /*
        if ( givenCell->isThere( currEdge ) == rg_TRUE ) {
            if ( currEdge->isThere( givenVertex ) == rg_TRUE ) { //multiplicity anomaly때문에 반드시 필요한 조건이다.
                outputCommonEdge = (AugmentedBCEdge*)currEdge;
                break;
            }
        }
        */
        if (    givenCell->isThere( currEdge ) 
             && currEdge->isThere( givenVertex )
             && noThisBoundingFace->isThere( currEdge ) ) { //multiplicity anomaly때문에 반드시 필요한 조건이다.
                outputCommonEdge = (AugmentedBCEdge*)currEdge;
            break;            
        }
    }

    AugmentedBCFace* boundingFacesincidentToGivenEdge[2];
    givenCell->findFacesIncidentToEdge( outputCommonEdge, boundingFacesincidentToGivenEdge );

    if ( boundingFacesincidentToGivenEdge[0] == noThisBoundingFace )
        outputBoundingFace = boundingFacesincidentToGivenEdge[1];
    else
        outputBoundingFace = boundingFacesincidentToGivenEdge[0];    
}


void AugmentedBetaComplex::setClassOfAugmentedBCEdgesToInteriorIfEdgesBoundRegularFaces(AugmentedBCVertex* givenVertex, rg_dList<AugmentedBCFace*>* faceGroup)
{
    AugmentedBCFace* currAugmentedBCFace = rg_NULL;
    faceGroup->reset4Loop();
    while ( faceGroup->setNext4Loop() )
    {
        currAugmentedBCFace = faceGroup->getEntity();

        BetaEdge* edgesBoundedByGivenVertex[2];
        currAugmentedBCFace->findEdgeIncidentToVertex( givenVertex, edgesBoundedByGivenVertex );

        rg_INT i=0;
        for (i=0; i<2; i++)
        {
            if ( ((AugmentedBCEdge*)edgesBoundedByGivenVertex[i])->getBoundingState() == REGULAR_SIMPLEX )
            {
                ((AugmentedBCEdge*)edgesBoundedByGivenVertex[i])->setBoundingState( INTERIOR_SIMPLEX );
            }
        }
    }
}



AugmentedBCFace* AugmentedBetaComplex::findRegularFaceOfGivenCellByBackwardOrientation(AugmentedBCFace* givenRegularFace, AugmentedBCEdge* givenEdge)
{
    AugmentedBCFace* startFace = givenRegularFace;
	 
    AugmentedBCFace* currFace  = startFace;
	AugmentedBCCell* currCell = rg_NULL;
    do  {
        if ( currFace->getEdgeOrientation(givenEdge) == rg_FALSE )  {
            currCell = (AugmentedBCCell*)currFace->getRightCell();

	
            if ( currCell == rg_NULL) 
            {
                return currFace;
            }
            
			currFace = (AugmentedBCFace*)currCell->findNextIncidentFaceOfEdge(givenEdge, currFace);
        }
        else if ( currFace->getEdgeOrientation(givenEdge) == rg_TRUE )  {
            currCell = (AugmentedBCCell*)currFace->getLeftCell();
            
            /*
            if ( currCell == rg_NULL )
            {
                return currFace;
            }
            */

            currFace = currCell->findNextIncidentFaceOfEdge(givenEdge, currFace);
        }
        else  {
            break;
        }
    } while ( currFace != startFace );

    return rg_NULL;
}




AugmentedBCFace* AugmentedBetaComplex::findRegularFaceOfGivenCellByForwardOrientation(AugmentedBCFace* givenRegularFace, AugmentedBCEdge* givenEdge)
{
    AugmentedBCFace* startFace = givenRegularFace;
	
    AugmentedBCFace* currFace  = startFace;
	AugmentedBCCell* currCell = rg_NULL;
    do  {
        if ( currFace->getEdgeOrientation(givenEdge) == rg_TRUE )  {
            currCell = (AugmentedBCCell*)currFace->getRightCell();

			if ( currCell == rg_NULL) 
            {
                return currFace;
            }
            
			currFace = currCell->findNextIncidentFaceOfEdge(givenEdge, currFace);
        }
        else if ( currFace->getEdgeOrientation(givenEdge) == rg_FALSE )  {
            currCell = (AugmentedBCCell*)currFace->getLeftCell();
            
            /*
            if ( currCell == rg_NULL )
            {
                return currFace;
            }
            */
            
			currFace = currCell->findNextIncidentFaceOfEdge(givenEdge, currFace);
        }
        else  {
            break;
        }
    } while ( currFace != startFace );

    return rg_NULL;
}





AugmentedBCCell* AugmentedBetaComplex::addArtificialCellByRegularAndSingularFaces(AugmentedBCFace* regularFace, AugmentedBCFace* singularFace, AugmentedBCEdge* commonEdge)
{
	AugmentedBCVertex* mateVertexOfRegularFace = (AugmentedBCVertex*)regularFace->findMateVertexOfEdge( commonEdge );
    AugmentedBCVertex* mateVertexOfSingularFace = (AugmentedBCVertex*)singularFace->findMateVertexOfEdge( commonEdge );

    AugmentedBCEdge* newEdge = addArtificialEdge( mateVertexOfSingularFace, mateVertexOfRegularFace );

    AugmentedBCVertex* startVertexOfCommonEdge = (AugmentedBCVertex*)commonEdge->getStartVertex();
    AugmentedBCEdge* edgeOfFace1IncidentToStartVertexOfCommonEdge = (AugmentedBCEdge*)regularFace->findIncidentEdgeOfVertexNotGivenEdge( startVertexOfCommonEdge, commonEdge );
    AugmentedBCEdge* edgeOfFace2IncidentToStartVertexOfCommonEdge = (AugmentedBCEdge*)singularFace->findIncidentEdgeOfVertexNotGivenEdge( startVertexOfCommonEdge, commonEdge );

    AugmentedBCFace* artificialFace1 = addArtificialFaceBy3Edges(newEdge, edgeOfFace1IncidentToStartVertexOfCommonEdge, edgeOfFace2IncidentToStartVertexOfCommonEdge);

    AugmentedBCVertex* endVertexOfCommonEdge = (AugmentedBCVertex*)commonEdge->getEndVertex();
    AugmentedBCEdge* edgeOfFace1IncidentToEndVertexOfCommonEdge = regularFace->findIncidentEdgeOfVertexNotGivenEdge( endVertexOfCommonEdge, commonEdge );
    AugmentedBCEdge* edgeOfFace2IncidentToEndVertexOfCommonEdge = singularFace->findIncidentEdgeOfVertexNotGivenEdge( endVertexOfCommonEdge, commonEdge );

    AugmentedBCFace* artificialFace2 = addArtificialFaceBy3Edges(newEdge, edgeOfFace2IncidentToEndVertexOfCommonEdge, edgeOfFace1IncidentToEndVertexOfCommonEdge);

    //regularFace 에 방향을 맞춘다.
    if (    regularFace->getEdgeOrientation( commonEdge )
         != singularFace->getEdgeOrientation( commonEdge) ) {
        singularFace->reverse();
    }
    
    if (    regularFace->getEdgeOrientation( edgeOfFace1IncidentToStartVertexOfCommonEdge )
         != artificialFace1->getEdgeOrientation( edgeOfFace1IncidentToStartVertexOfCommonEdge )) {
        artificialFace1->reverse();
    }

    if (    regularFace->getEdgeOrientation( edgeOfFace1IncidentToEndVertexOfCommonEdge )
         != artificialFace2->getEdgeOrientation( edgeOfFace1IncidentToEndVertexOfCommonEdge ) )  {
        artificialFace2->reverse();
    }

    
    AugmentedBCCell tempCell(regularFace, singularFace, artificialFace1, artificialFace2);
    tempCell.setIsArtificial( rg_TRUE );    
    tempCell.setBoundingState( INTERIOR_SIMPLEX );

    tempCell.setID( getMaxIndexOfAugmentedBCCells() + 1 );
    

    AugmentedBCCell* newCell = m_artificialAugmentedBCCells.add( tempCell );
    
    if ( regularFace->getBoundingState() == SINGULAR_SIMPLEX )
    {
        regularFace->reverse();
        regularFace->setLeftCell( newCell);
        //regularFace->reverseEdgeOrientations();
        regularFace->setBoundingState( REGULAR_SIMPLEX );
    }
    else if ( regularFace->getBoundingState() == REGULAR_SIMPLEX )
    {
        regularFace->setRightCell( newCell );
        regularFace->setBoundingState( INTERIOR_SIMPLEX );
    }
    
    singularFace->setLeftCell( newCell ); //regularFace에 방향을 맞췄기 때문.
    artificialFace1->setLeftCell( newCell );
    artificialFace2->setLeftCell( newCell );    
   
    
    if ( singularFace->getBoundingState() == SINGULAR_SIMPLEX )
        singularFace->setBoundingState( REGULAR_SIMPLEX );
    else if ( singularFace->getBoundingState() == REGULAR_SIMPLEX )
        singularFace->setBoundingState( INTERIOR_SIMPLEX );

    artificialFace1->setBoundingState( REGULAR_SIMPLEX );
    artificialFace2->setBoundingState( REGULAR_SIMPLEX );

    newEdge->setBoundingState( REGULAR_SIMPLEX );
    newEdge->setFirstFace( artificialFace1 );
    
    return newCell;
}


AugmentedBCCell* AugmentedBetaComplex::addArtificialCellBy1EdgeAnd1Face(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenFace, AugmentedBCVertex* commonVertex)
{
    BetaEdge* betaEdgesBoundedByCommonVertex[2];
    givenFace->findEdgeIncidentToVertex( (BetaVertex*)commonVertex, betaEdgesBoundedByCommonVertex );

	AugmentedBCEdge* edgeBoundedByCommonVertex[2]
		= {(AugmentedBCEdge*)betaEdgesBoundedByCommonVertex[0], (AugmentedBCEdge*)betaEdgesBoundedByCommonVertex[1]};

    
    AugmentedBCFace* artificialFace1 = addArtificialFaceBy2Edges( givenEdge, edgeBoundedByCommonVertex[0], commonVertex );
    AugmentedBCFace* artificialFace2 = addArtificialFaceBy2Edges( edgeBoundedByCommonVertex[1], givenEdge, commonVertex );

    AugmentedBCEdge* edgeNotBoundedByCommonVertex = givenFace->findMateEdgeOfVertex( commonVertex );
    AugmentedBCEdge* edgeOfFace1NotBoundedByCommonVertex = artificialFace1->findMateEdgeOfVertex( commonVertex );
    AugmentedBCEdge* edgeOfFace2NotBoundedByCommonVertex = artificialFace2->findMateEdgeOfVertex( commonVertex );
    
    AugmentedBCFace* artificialFace3 = addArtificialFaceBy3Edges(edgeNotBoundedByCommonVertex, edgeOfFace1NotBoundedByCommonVertex, edgeOfFace2NotBoundedByCommonVertex);

    if (    givenFace->getEdgeOrientation( edgeBoundedByCommonVertex[0] ) 
         != artificialFace1->getEdgeOrientation( edgeBoundedByCommonVertex[0] )) {
        artificialFace1->reverse();
    }

    if (    givenFace->getEdgeOrientation( edgeBoundedByCommonVertex[1] ) 
         != artificialFace2->getEdgeOrientation( edgeBoundedByCommonVertex[1] )) {
        artificialFace2->reverse();
    }

    if (    givenFace->getEdgeOrientation( edgeNotBoundedByCommonVertex )
         != artificialFace3->getEdgeOrientation( edgeNotBoundedByCommonVertex )) {
        artificialFace3->reverse();
    }

    AugmentedBCCell tempCell(givenFace, artificialFace1, artificialFace2, artificialFace3);
    tempCell.setIsArtificial( rg_TRUE );    
    tempCell.setBoundingState( INTERIOR_SIMPLEX );

    tempCell.setID( getMaxIndexOfAugmentedBCCells() + 1 );
    

    AugmentedBCCell* newCell = m_artificialAugmentedBCCells.add( tempCell );
    
    if ( givenFace->getBoundingState() == SINGULAR_SIMPLEX )
    {
        givenFace->reverse();
        givenFace->setLeftCell( newCell);
        //givenFace->reverseEdgeOrientations();
        givenFace->setBoundingState( REGULAR_SIMPLEX );
    }
    else if ( givenFace->getBoundingState() == REGULAR_SIMPLEX )
    {
        givenFace->setRightCell( newCell );    
        givenFace->setBoundingState( INTERIOR_SIMPLEX );
    }
    
    artificialFace1->setLeftCell( newCell );
    artificialFace2->setLeftCell( newCell );
    artificialFace3->setLeftCell( newCell );

    artificialFace1->setBoundingState( REGULAR_SIMPLEX );
    artificialFace2->setBoundingState( REGULAR_SIMPLEX );
    artificialFace3->setBoundingState( REGULAR_SIMPLEX );

    return newCell;
}

AugmentedBCVertex* AugmentedBetaComplex::replicateGivenVertex(AugmentedBCVertex* givenVertex)
{
    AugmentedBCVertex tempArtificialVertex(*givenVertex);
    tempArtificialVertex.setIsArtificial(rg_TRUE);    
	tempArtificialVertex.setOriginalBetaVertex(givenVertex->getOriginalBetaVertex());

    tempArtificialVertex.setID( getMaxIndexOfAugmentedBCVertices() + 1 );
        
    AugmentedBCVertex* artificialVertex = m_artificialAugmentedBCVertices.add( tempArtificialVertex );

    return artificialVertex;
}

AugmentedBCEdge* AugmentedBetaComplex::addArtificialEdge(AugmentedBCVertex* givenVertex, AugmentedBCVertex* artificialVertex)
{
    AugmentedBCEdge tempArtificialEdge(givenVertex, artificialVertex);
    tempArtificialEdge.setIsArtificial( rg_TRUE );    
    tempArtificialEdge.setBoundingState( SINGULAR_SIMPLEX );

    tempArtificialEdge.setID( getMaxIndexOfAugmentedBCEdges() + 1 );
    
	AugmentedBCEdge* artificialEdge = m_artificialAugmentedBCEdges.add( tempArtificialEdge );
    artificialVertex->setFirstEdge( artificialEdge );
    givenVertex->setFirstEdge( artificialEdge );

    return artificialEdge;
}


AugmentedBCFace* AugmentedBetaComplex::addArtificialFaceBy2Edges(AugmentedBCEdge* e1, AugmentedBCEdge* e2, AugmentedBCVertex* commonVertex)
{    
    AugmentedBCVertex* v1 = e1->getBoundingVertexNotGivenVertex(commonVertex);
    AugmentedBCVertex* v2 = e2->getBoundingVertexNotGivenVertex(commonVertex);

    AugmentedBCEdge* artificialEdge = addArtificialEdge( v2, v1 );

    AugmentedBCFace tempNewFace(e1, e2, artificialEdge);
    tempNewFace.setIsArtificial( rg_TRUE );

    tempNewFace.setID( getMaxIndexOfAugmentedBCFaces() + 1 );
    
    AugmentedBCFace* artificialFace = m_artificialAugmentedBCFaces.add( tempNewFace );

    artificialEdge->setFirstFace( artificialFace );
    e1->setFirstFace( artificialFace );
    e2->setFirstFace( artificialFace );

    artificialEdge->setBoundingState( REGULAR_SIMPLEX );
    e1->setBoundingState( REGULAR_SIMPLEX );
    e2->setBoundingState( REGULAR_SIMPLEX );

    artificialFace->setBoundingState( SINGULAR_SIMPLEX );

    return artificialFace;
}


AugmentedBCFace* AugmentedBetaComplex::addArtificialFaceBy3Edges(AugmentedBCEdge* e1, AugmentedBCEdge* e2, AugmentedBCEdge* e3)
{
    AugmentedBCFace tempFace(e1, e2, e3);
    tempFace.setIsArtificial( rg_TRUE );
    
    tempFace.setID( getMaxIndexOfAugmentedBCFaces() + 1 );
    
    AugmentedBCFace* newFace = m_artificialAugmentedBCFaces.add( tempFace );

    e1->setFirstFace( newFace );
    e2->setFirstFace( newFace );
    e3->setFirstFace( newFace );

    return newFace;
}



AugmentedBCFace* AugmentedBetaComplex::findIncidentRegularFaceOfGivenCellAndVertex(AugmentedBCCell* givenCell, AugmentedBCVertex* givenVertex)
{
    rg_dList<BetaCell*> cellList;    

    rg_dList< pair<BetaFace*, BetaCell*> > facesToBeTraced;

    pair< BetaFace*, BetaCell* > currFace;
    BetaFace** boundingFaceOfFirstCell = givenCell->getFaces();
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {

        if ( boundingFaceOfFirstCell[i]->isThere( (BetaVertex*)givenVertex ) )  {
            currFace.first  = boundingFaceOfFirstCell[i];

            if ( ((AugmentedBCFace*)currFace.first)->getBoundingState() == REGULAR_SIMPLEX )
                return (AugmentedBCFace*)currFace.first;

            if ( givenCell == (AugmentedBCCell*)boundingFaceOfFirstCell[i]->getLeftCell() )  
                currFace.second = boundingFaceOfFirstCell[i]->getRightCell();
            else  
                currFace.second = boundingFaceOfFirstCell[i]->getLeftCell();

            facesToBeTraced.add( currFace );
        }
    }

    cellList.add( givenCell );
    while ( facesToBeTraced.getSize() > 0 )  {
        currFace = facesToBeTraced.popFront();

        if ( currFace.second == rg_NULL )
            continue;

        if ( cellList.isInList( currFace.second ) == rg_FALSE )  {        
            cellList.add( currFace.second );

            BetaFace** boundingFace = currFace.second->getFaces();
            for ( rg_INT i=0; i<4; i++ )  {
                if ( boundingFace[i] == currFace.first )
                    continue;

                if ( boundingFace[i]->isThere( givenVertex ) )  {

                    pair<BetaFace*, BetaCell*> insertingFace;
                    insertingFace.first  = boundingFace[i];

                    if ( ((AugmentedBCFace*)insertingFace.first)->getBoundingState() == REGULAR_SIMPLEX )
                    {
                        cellList.removeAll();
                        return (AugmentedBCFace*)insertingFace.first;
                    }

                    if ( currFace.second == boundingFace[i]->getLeftCell() )  
                        insertingFace.second = boundingFace[i]->getRightCell();
                    else  
                        insertingFace.second = boundingFace[i]->getLeftCell();

                    facesToBeTraced.add( insertingFace );
                }
            }
        }
    }

    return rg_NULL;
}



void AugmentedBetaComplex::extractBoundaryOfAugmentedBetaComplex( ManifoldizedBetaShape& manifoldizedBS)
{
    if ( m_augmentedBCVertices.getSize() == 0 ) {
        return;
    }

    if ( manifoldizedBS.getBodies()->getSize() != 0 )
        manifoldizedBS.removeAll();
    
	manifoldizedBS.setBetaValue( m_betaValue );

    //adjustIDOfAugmentedBetaComplex();    
    
    
    
    MBSVertex** ptMBSVertices = rg_NULL;
    MBSEdge**   ptMBSEdges    = rg_NULL;
    MBSFace**   ptMBSFaces    = rg_NULL;    
    initializePointers( ptMBSVertices, ptMBSEdges, ptMBSFaces );

    
    rg_dList< MBSVertex* >  mbsVertices;
    rg_dList< MBSEdge* >    mbsEdges;
    rg_dList< MBSFace* >    mbsFaces;
    duplicateBetaShapeSimplexesOfAugBCtoMBS( ptMBSVertices, ptMBSEdges, ptMBSFaces,
                                             mbsVertices, mbsEdges, mbsFaces );


    buildTopologyOfMBSEdges( ptMBSVertices, ptMBSEdges, ptMBSFaces );


    rg_dList< MBSShell* > mbsShells;
    divideMBSSimplexesIntoShells( mbsVertices, mbsEdges, mbsFaces, mbsShells );

    
    divideMBSShellsIntoBodies( mbsShells, manifoldizedBS );


    divideExteriorShellFromInteriorShells( manifoldizedBS );


    if ( ptMBSVertices != rg_NULL )
        delete [] ptMBSVertices;
    
    if ( ptMBSEdges != rg_NULL )
        delete [] ptMBSEdges;
    
    if ( ptMBSFaces != rg_NULL )
        delete [] ptMBSFaces;
}



rg_INT AugmentedBetaComplex::getMaxIndexOfAugmentedBCVertices()
{
    rg_INT maxIndexOfAugmentedBCVertex = 0;
    if ( m_artificialAugmentedBCVertices.getSize() == 0) {
        if ( m_augmentedBCVertices.getSize() == 0 ) {
            maxIndexOfAugmentedBCVertex = -1;
        }
        else {
            maxIndexOfAugmentedBCVertex = m_augmentedBCVertices.getTail()->getpEntity()->getID();
        }
    }
    else {
        maxIndexOfAugmentedBCVertex = m_artificialAugmentedBCVertices.getTail()->getpEntity()->getID();
    }

    return maxIndexOfAugmentedBCVertex;
}


rg_INT AugmentedBetaComplex::getMaxIndexOfAugmentedBCEdges()
{
    rg_INT maxIndexOfAugmentedBCEdge = 0;
    if ( m_artificialAugmentedBCEdges.getSize() == 0) {
        if ( m_augmentedBCEdges.getSize() == 0 ) {
            maxIndexOfAugmentedBCEdge = -1;
        }
        else {
            maxIndexOfAugmentedBCEdge = m_augmentedBCEdges.getTail()->getpEntity()->getID();
        }
    }
    else {
        maxIndexOfAugmentedBCEdge = m_artificialAugmentedBCEdges.getTail()->getpEntity()->getID();
    }

    return maxIndexOfAugmentedBCEdge;
}

rg_INT AugmentedBetaComplex::getMaxIndexOfAugmentedBCFaces()
{
    rg_INT maxIndexOfAugmentedBCFace = 0;
    if ( m_artificialAugmentedBCFaces.getSize() == 0) {
        if ( m_augmentedBCFaces.getSize() == 0 ) {
            maxIndexOfAugmentedBCFace = -1;
        }
        else {
            maxIndexOfAugmentedBCFace = m_augmentedBCFaces.getTail()->getpEntity()->getID();
        }
    }
    else {
        maxIndexOfAugmentedBCFace = m_artificialAugmentedBCFaces.getTail()->getpEntity()->getID();
    }

    return maxIndexOfAugmentedBCFace;
}


rg_INT AugmentedBetaComplex::getMaxIndexOfAugmentedBCCells()
{
    rg_INT maxIndexOfAugmentedBCCell = 0;
    if ( m_artificialAugmentedBCCells.getSize() == 0) {
        if ( m_augmentedBCCells.getSize() == 0 ) {
            maxIndexOfAugmentedBCCell = -1;
        }
        else {
            maxIndexOfAugmentedBCCell = m_augmentedBCCells.getTail()->getpEntity()->getID();
        }
    }
    else {
        maxIndexOfAugmentedBCCell = m_artificialAugmentedBCCells.getTail()->getpEntity()->getID();
    }

    return maxIndexOfAugmentedBCCell;
}



void AugmentedBetaComplex::initializePointers( MBSVertex**& ptMBSVertices, 
                                    MBSEdge**& ptMBSEdges, MBSFace**& ptMBSFaces )
{
    // 1. create
    rg_INT maxIndexOfAugmentedBCVertex = getMaxIndexOfAugmentedBCVertices();
    rg_INT maxIndexOfAugmentedBCEdge   = getMaxIndexOfAugmentedBCEdges();
    rg_INT maxIndexOfAugmentedBCFace   = getMaxIndexOfAugmentedBCFaces();    

    if ( maxIndexOfAugmentedBCVertex + 1 != 0 ) {
        ptMBSVertices = new MBSVertex*[ maxIndexOfAugmentedBCVertex + 1 ];
    }
        
    if ( maxIndexOfAugmentedBCEdge + 1 != 0 ) {
        ptMBSEdges = new MBSEdge*[   maxIndexOfAugmentedBCEdge + 1 ];
    }

    if ( maxIndexOfAugmentedBCFace + 1 != 0 ) {
        ptMBSFaces = new MBSFace*[   maxIndexOfAugmentedBCFace + 1 ];
    }


    // 2. set pointers to MBS-simplex rg_NULL
    rg_INT i = 0;

    rg_INT numOfBCVertices = getMaxIndexOfAugmentedBCVertices() + 1;
    for ( i=0; i < numOfBCVertices; i++)
        ptMBSVertices[i] = rg_NULL;

    rg_INT numOfBCEdges = getMaxIndexOfAugmentedBCEdges() + 1;
    for ( i=0; i < numOfBCEdges; i++)
        ptMBSEdges[i] = rg_NULL;

    rg_INT numOfBCFaces = getMaxIndexOfAugmentedBCFaces() + 1;
    for ( i=0; i < numOfBCFaces; i++)
        ptMBSFaces[i] = rg_NULL;
}



void AugmentedBetaComplex::duplicateBetaShapeSimplexesOfAugBCtoMBS( MBSVertex** ptMBSVertices, MBSEdge** ptMBSEdges, MBSFace** ptMBSFaces, 
                                                                    rg_dList< MBSVertex* >& mbsVertices, rg_dList< MBSEdge* >& mbsEdges, rg_dList< MBSFace* >& mbsFaces )
{
    duplicateVerticesOfBetaShapeInAugBCtoMBS( ptMBSVertices, mbsVertices );
    duplicateEdgesOfBetaShapeInAugBCtoMBS( ptMBSVertices, ptMBSEdges, mbsEdges );
    duplicateFacesOfBetaShapeInAugBCtoMBS( ptMBSEdges, ptMBSFaces, mbsFaces );
}


void AugmentedBetaComplex::duplicateVerticesOfBetaShapeInAugBCtoMBS( MBSVertex** ptMBSVertices, rg_dList< MBSVertex* >& mbsVertices )
{
    m_augmentedBCVertices.reset4Loop();
    while( m_augmentedBCVertices.setNext4Loop() )
    {
        AugmentedBCVertex* currAugmentedBCVertex = m_augmentedBCVertices.getpEntity();

        rg_INT boundingStateOfCurrVertex = currAugmentedBCVertex->getBoundingState();

        if (   boundingStateOfCurrVertex == REGULAR_SIMPLEX 
            || boundingStateOfCurrVertex == SINGULAR_SIMPLEX ) {
            rg_INT currAugmentedBCVertexID = currAugmentedBCVertex->getID();

            MBSVertex* tempMBSVertex= new MBSVertex( currAugmentedBCVertexID );
            
            mbsVertices.add( tempMBSVertex );
            ptMBSVertices[ currAugmentedBCVertexID ] = tempMBSVertex;

            tempMBSVertex->setOriginalBetaVertex( currAugmentedBCVertex->getOriginalBetaVertex() );
            tempMBSVertex->isArtificial( currAugmentedBCVertex->isArtificial() );
            tempMBSVertex->isVisited( rg_FALSE );
        }
    }


    m_artificialAugmentedBCVertices.reset4Loop();
    while( m_artificialAugmentedBCVertices.setNext4Loop() )
    {
        AugmentedBCVertex* currAugmentedBCVertex = m_artificialAugmentedBCVertices.getpEntity();

        rg_INT boundingStateOfCurrVertex = currAugmentedBCVertex->getBoundingState();

        if (   boundingStateOfCurrVertex == REGULAR_SIMPLEX 
            || boundingStateOfCurrVertex == SINGULAR_SIMPLEX ) {
            rg_INT currAugmentedBCVertexID = currAugmentedBCVertex->getID();

            MBSVertex* tempMBSVertex = new MBSVertex( currAugmentedBCVertexID);            
            
            mbsVertices.add( tempMBSVertex );
            ptMBSVertices[ currAugmentedBCVertexID ] = tempMBSVertex;

            tempMBSVertex->setOriginalBetaVertex( currAugmentedBCVertex->getOriginalBetaVertex() );
            tempMBSVertex->isArtificial( currAugmentedBCVertex->isArtificial() );
            tempMBSVertex->isVisited( rg_FALSE );
        }
    }


}





void AugmentedBetaComplex::duplicateEdgesOfBetaShapeInAugBCtoMBS( MBSVertex** ptMBSVertices, MBSEdge** ptMBSEdges,
                                                                  rg_dList< MBSEdge* >& mbsEdges )
{
    m_augmentedBCEdges.reset4Loop();
    while( m_augmentedBCEdges.setNext4Loop() )
    {
        AugmentedBCEdge* currAugmentedBCEdge = m_augmentedBCEdges.getpEntity();

        rg_INT boundingStateOfCurrEdge = currAugmentedBCEdge->getBoundingState();

        if (    boundingStateOfCurrEdge == REGULAR_SIMPLEX
             || boundingStateOfCurrEdge == SINGULAR_SIMPLEX) {
            rg_INT currAugmentedBCEdgeID = currAugmentedBCEdge->getID();
            MBSVertex* startVertex = ptMBSVertices[ currAugmentedBCEdge->getStartVertex()->getID() ];
            MBSVertex* endVertex   = ptMBSVertices[ currAugmentedBCEdge->getEndVertex()->getID() ];
            
            MBSEdge* tempMBSEdge = new MBSEdge( currAugmentedBCEdgeID, startVertex, endVertex);
            
            mbsEdges.add( tempMBSEdge );
            ptMBSEdges[ currAugmentedBCEdgeID ] = tempMBSEdge;
            
            //tempMBSEdge->setOriginalBetaEdge( currAugmentedBCEdge->getOriginalBetaEdge() );
            tempMBSEdge->isArtificial( currAugmentedBCEdge->isArtificial() );
            tempMBSEdge->isVisited( rg_FALSE );

            startVertex->setFirstEdge( tempMBSEdge );
            endVertex->setFirstEdge( tempMBSEdge );
        }
    }


    m_artificialAugmentedBCEdges.reset4Loop();
    while( m_artificialAugmentedBCEdges.setNext4Loop() )
    {
        AugmentedBCEdge*currAugmentedBCEdge = m_artificialAugmentedBCEdges.getpEntity();

        rg_INT boundingStateOfCurrEdge = currAugmentedBCEdge->getBoundingState();

        if (    boundingStateOfCurrEdge == REGULAR_SIMPLEX
             || boundingStateOfCurrEdge == SINGULAR_SIMPLEX) {
            rg_INT currAugmentedBCEdgeID = currAugmentedBCEdge->getID();
            MBSVertex* startVertex = ptMBSVertices[ currAugmentedBCEdge->getStartVertex()->getID() ];
            MBSVertex* endVertex   = ptMBSVertices[ currAugmentedBCEdge->getEndVertex()->getID() ];
            
            MBSEdge* tempMBSEdge = new MBSEdge( currAugmentedBCEdgeID, startVertex, endVertex);
            
            mbsEdges.add( tempMBSEdge );
            ptMBSEdges[ currAugmentedBCEdgeID ] = tempMBSEdge;

            //tempMBSEdge->setOriginalBetaEdge( currAugmentedBCEdge->getOriginalBetaEdge() );
            tempMBSEdge->isArtificial( currAugmentedBCEdge->isArtificial() );
            tempMBSEdge->isVisited( rg_FALSE );
            
            startVertex->setFirstEdge( ptMBSEdges[ currAugmentedBCEdgeID ] );
            endVertex->setFirstEdge( ptMBSEdges[ currAugmentedBCEdgeID ] );
        }
    }
}



void AugmentedBetaComplex::duplicateFacesOfBetaShapeInAugBCtoMBS( MBSEdge** ptMBSEdges, MBSFace** ptMBSFaces, 
                                                                  rg_dList< MBSFace* >& mbsFaces )
{
    AugmentedBCFace*   currAugmentedBCFace = rg_NULL;
    rg_INT     currAugmentedBCFaceID;
    
    m_augmentedBCFaces.reset4Loop();
    while( m_augmentedBCFaces.setNext4Loop() )
    {
        currAugmentedBCFace = m_augmentedBCFaces.getpEntity();

        rg_INT boundingStateOfCurrFace = currAugmentedBCFace->getBoundingState();

        if (    boundingStateOfCurrFace == REGULAR_SIMPLEX
             || boundingStateOfCurrFace == SINGULAR_SIMPLEX) {
            currAugmentedBCFaceID = currAugmentedBCFace->getID();
            MBSEdge* firstEdge = ptMBSEdges[ currAugmentedBCFace->getEdge(0)->getID() ];
                        
            MBSFace* tempMBSFace = new MBSFace( currAugmentedBCFaceID );
            tempMBSFace->setFirstEdge( firstEdge );
            //tempMBSFace->setOriginalBetaFace( currAugmentedBCFace->getOriginalBetaFace() );
            tempMBSFace->isArtificial( currAugmentedBCFace->isArtificial() );
            tempMBSFace->isVisited( rg_FALSE );

            mbsFaces.add( tempMBSFace );
            ptMBSFaces[ currAugmentedBCFaceID ] = tempMBSFace;            
        }
    }


    m_artificialAugmentedBCFaces.reset4Loop();
    while( m_artificialAugmentedBCFaces.setNext4Loop() )
    {
        currAugmentedBCFace = m_artificialAugmentedBCFaces.getpEntity();

        rg_INT boundingStateOfCurrFace = currAugmentedBCFace->getBoundingState();

        if (    boundingStateOfCurrFace == REGULAR_SIMPLEX
             || boundingStateOfCurrFace == SINGULAR_SIMPLEX) {
            currAugmentedBCFaceID = currAugmentedBCFace->getID();
            MBSEdge* firstEdge = ptMBSEdges[ currAugmentedBCFace->getEdge(0)->getID() ];
                        
            MBSFace* tempMBSFace = new MBSFace( currAugmentedBCFaceID );
            tempMBSFace->setFirstEdge( firstEdge );
            //tempMBSFace->setOriginalBetaFace( currAugmentedBCFace->getOriginalBetaFace() );
            tempMBSFace->isArtificial( currAugmentedBCFace->isArtificial() );
            tempMBSFace->isVisited( rg_FALSE );

            mbsFaces.add( tempMBSFace );
            ptMBSFaces[ currAugmentedBCFaceID ] = tempMBSFace;  
        }
    }

}


void AugmentedBetaComplex::buildTopologyOfMBSEdges( MBSVertex** ptMBSVertices, MBSEdge** ptMBSEdges, MBSFace** ptMBSFaces )
{
    m_augmentedBCEdges.reset4Loop();
    while( m_augmentedBCEdges.setNext4Loop() )
    {
        AugmentedBCEdge* currAugmentedBCEdge = m_augmentedBCEdges.getpEntity();


        if ( currAugmentedBCEdge->getBoundingState() == REGULAR_SIMPLEX ) {
            buildTopologyOfMBSEdge( currAugmentedBCEdge, ptMBSVertices, ptMBSEdges, ptMBSFaces );
        }
    }


    m_artificialAugmentedBCEdges.reset4Loop();
    while( m_artificialAugmentedBCEdges.setNext4Loop() )
    {
        AugmentedBCEdge*currAugmentedBCEdge = m_artificialAugmentedBCEdges.getpEntity();

        if ( currAugmentedBCEdge->getBoundingState() == REGULAR_SIMPLEX ) {
            buildTopologyOfMBSEdge( currAugmentedBCEdge, ptMBSVertices, ptMBSEdges, ptMBSFaces );
        }
    }
}



void AugmentedBetaComplex::buildTopologyOfMBSEdge( AugmentedBCEdge* currAugmentedBCEdge, MBSVertex** ptMBSVertices, MBSEdge** ptMBSEdges, MBSFace** ptMBSFaces )
{
    if ( currAugmentedBCEdge->getFirstFace() == rg_NULL ) {
        return;
    }

    rg_INT currAugmentedBCEdgeID = currAugmentedBCEdge->getID();

    rg_INT boundingStateOfFirstFaceOfCurrRegularEdge 
            = ((AugmentedBCFace*) currAugmentedBCEdge->getFirstFace())->getBoundingState();

    if ( boundingStateOfFirstFaceOfCurrRegularEdge == SINGULAR_SIMPLEX ) { //currEdge bounds SINGULAR face only
        AugmentedBCFace* singularFace = (AugmentedBCFace*) currAugmentedBCEdge->getFirstFace();

        if ( singularFace->getEdgeOrientation( currAugmentedBCEdge ) == rg_TRUE ) { //==> singularFace is a right-face of currEdge
            ptMBSEdges[ currAugmentedBCEdgeID ]->setLeftFace( ptMBSFaces[ singularFace->getID() ] );

            BetaEdge** boundingEdges = singularFace->getEdges();        
            for ( rg_INT i=0; i<3; i++ )
            {
			    AugmentedBCEdge* currEdge = (AugmentedBCEdge*)boundingEdges[i];
                if ( currEdge == currAugmentedBCEdge ) {
                    continue;
                }
                else if ( currEdge->isThere( currAugmentedBCEdge->getStartVertex() ) ) {
                    ptMBSEdges[ currAugmentedBCEdgeID ]->setLeftLeg( ptMBSEdges[ currEdge->getID() ] );
                }
                else if ( currEdge->isThere( currAugmentedBCEdge->getEndVertex() ) ) {
                    ptMBSEdges[ currAugmentedBCEdgeID ]->setLeftHand( ptMBSEdges[ currEdge->getID() ] );
                }
            }
        }
        else { // if( singularFace->getEdgeOrientation( currAugmentedBCEdge ) == rg_FALSE ) ==> singularFace is a right-face of currEdge
            ptMBSEdges[ currAugmentedBCEdgeID ]->setRightFace( ptMBSFaces[ singularFace->getID() ] );    

            BetaEdge** boundingEdges = singularFace->getEdges();        
            for ( rg_INT i=0; i<3; i++ )
            {
                if ( boundingEdges[i] == currAugmentedBCEdge )
                    continue;
                else if ( boundingEdges[i]->isThere( currAugmentedBCEdge->getStartVertex() ) )
                    ptMBSEdges[ currAugmentedBCEdgeID ]->setRightLeg( ptMBSEdges[ boundingEdges[i]->getID() ] );
                else if ( boundingEdges[i]->isThere( currAugmentedBCEdge->getEndVertex() ) )
                    ptMBSEdges[ currAugmentedBCEdgeID ]->setRightHand( ptMBSEdges[ boundingEdges[i]->getID() ] );
            }
        }
    }
    else { //currEdge bounds REGULAR or INTERIOR face --> currEdge bounds ordinary regular faces
        rg_dList< AugmentedBCFace* > incidentFaces;
        findAugRegularFacesIncidentToAugEdge( currAugmentedBCEdge, incidentFaces);
    
        AugmentedBCFace* leftAugmentedBCFace = rg_NULL;
        AugmentedBCFace* rightAugmentedBCFace = rg_NULL;

        incidentFaces.reset4Loop();
        AugmentedBCFace* currIncidentFace = rg_NULL;
        while( incidentFaces.setNext4Loop() )
        {
            currIncidentFace = incidentFaces.getEntity();

            if ( currIncidentFace->getEdgeOrientation( currAugmentedBCEdge ) == rg_TRUE )
                leftAugmentedBCFace = currIncidentFace;
            else if ( currIncidentFace->getEdgeOrientation( currAugmentedBCEdge ) == rg_FALSE )
                rightAugmentedBCFace = currIncidentFace;                
        }

        ptMBSEdges[ currAugmentedBCEdgeID ]->setLeftFace( ptMBSFaces[ leftAugmentedBCFace->getID() ] );
        ptMBSEdges[ currAugmentedBCEdgeID ]->setRightFace( ptMBSFaces[ rightAugmentedBCFace->getID() ] );

        rg_INT i=0;
        BetaEdge** boundingEdges = leftAugmentedBCFace->getEdges();        
        for ( i=0; i<3; i++ )
        {
			AugmentedBCEdge* currEdge = (AugmentedBCEdge*)boundingEdges[i];
            if ( currEdge == currAugmentedBCEdge )
                continue;
            else if ( currEdge->isThere( currAugmentedBCEdge->getStartVertex() ) )
                ptMBSEdges[ currAugmentedBCEdgeID ]->setLeftLeg( ptMBSEdges[ currEdge->getID() ] );
            else if ( currEdge->isThere( currAugmentedBCEdge->getEndVertex() ) )
                ptMBSEdges[ currAugmentedBCEdgeID ]->setLeftHand( ptMBSEdges[ currEdge->getID() ] );
        }

        boundingEdges = rightAugmentedBCFace->getEdges();        
        for ( i=0; i<3; i++ )
        {
            if ( boundingEdges[i] == currAugmentedBCEdge )
                continue;
            else if ( boundingEdges[i]->isThere( currAugmentedBCEdge->getStartVertex() ) )
                ptMBSEdges[ currAugmentedBCEdgeID ]->setRightLeg( ptMBSEdges[ boundingEdges[i]->getID() ] );
            else if ( boundingEdges[i]->isThere( currAugmentedBCEdge->getEndVertex() ) )
                ptMBSEdges[ currAugmentedBCEdgeID ]->setRightHand( ptMBSEdges[ boundingEdges[i]->getID() ] );
        }
    }
}





void AugmentedBetaComplex::divideMBSSimplexesIntoShells( rg_dList< MBSVertex* >& mbsVertices, rg_dList< MBSEdge* >& mbsEdges, rg_dList< MBSFace* >& mbsFaces, 
                                                         rg_dList< MBSShell* >& mbsShells )
{
    rg_INT indexOfShell = 0;

    mbsFaces.reset4Loop();
    while ( mbsFaces.setNext4Loop() ) {
        MBSFace* currFace = mbsFaces.getEntity();

        if ( currFace->isVisited() ) {
            continue;
        }

        MBSShell* currShell = new MBSShell(indexOfShell);
        indexOfShell++;
        mbsShells.add( currShell );

        rg_dList< MBSFace* > faceStack;
        faceStack.pushBack( currFace );

        while ( faceStack.getSize() > 0 ) {
            MBSFace* popedFace = faceStack.popBack();

            if ( popedFace->isVisited() ) {
                continue;
            }

            popedFace->isVisited( rg_TRUE );
            currShell->addFace( popedFace );



            rg_dList< MBSEdge* > boundingEdges;
            popedFace->searchBoundingEdges( boundingEdges );
            boundingEdges.reset4Loop();
            while ( boundingEdges.setNext4Loop() ) {
                MBSEdge* currBoundingEdge = boundingEdges.getEntity();

                if ( currBoundingEdge->isVisited() ) {
                    continue;
                }

                currBoundingEdge->isVisited( rg_TRUE );
                currShell->addEdge( currBoundingEdge );
            }



            rg_dList< MBSVertex* > boundingVertices;
            popedFace->searchBoundingVertices( boundingVertices );
            boundingVertices.reset4Loop();
            while ( boundingVertices.setNext4Loop() ) {
                MBSVertex* currBoundingVertex = boundingVertices.getEntity();

                if ( currBoundingVertex->isVisited() ) {
                    continue;
                }

                currBoundingVertex->isVisited( rg_TRUE );
                currShell->addVertex( currBoundingVertex );

                currBoundingVertex->setShell( currShell );
            }



            rg_dList< MBSFace* > incidentFaces;
            popedFace->searchIncidentFaces( incidentFaces );

            incidentFaces.reset4Loop();
            while ( incidentFaces.setNext4Loop() ) {
                MBSFace* currIncidentFace = incidentFaces.getEntity();

                if ( currIncidentFace->isVisited() ) {
                    continue;
                }

                faceStack.pushBack( currIncidentFace );
            }
        }
    }


    mbsEdges.reset4Loop();
    while ( mbsEdges.setNext4Loop() ) {
        MBSEdge* currEdge = mbsEdges.getEntity();

        if ( currEdge->isVisited() ) {
            continue;
        }

        //currEdge is a stand-alone edge
        MBSShell* currShell = new MBSShell(indexOfShell);
        indexOfShell++;
        mbsShells.add( currShell );

        currEdge->isVisited( rg_TRUE );
        currShell->addEdge( currEdge );

        MBSVertex* startVertex = currEdge->getStartVertex();
        currShell->addVertex( startVertex );
        startVertex->isVisited( rg_TRUE );
        startVertex->setShell( currShell );


        MBSVertex* endVertex = currEdge->getEndVertex();
        currShell->addVertex( endVertex );
        endVertex->isVisited( rg_TRUE );
        endVertex->setShell( currShell );
    }



    mbsVertices.reset4Loop();
    while ( mbsVertices.setNext4Loop() ) {
        MBSVertex* currVertex = mbsVertices.getEntity();

        if ( currVertex->isVisited() ) {
            continue;
        }

        //currVertex is a stand-alone vertex
        MBSShell* currShell = new MBSShell(indexOfShell);
        indexOfShell++;
        mbsShells.add( currShell );

        currVertex->isVisited( rg_TRUE );
        currShell->addVertex( currVertex );
        currVertex->setShell( currShell );
    }



    
    mbsVertices.reset4Loop();
    while ( mbsVertices.setNext4Loop() ) {
        mbsVertices.getEntity()->isVisited( rg_FALSE );
    }

    mbsEdges.reset4Loop();
    while ( mbsEdges.setNext4Loop() ) {
        mbsEdges.getEntity()->isVisited( rg_FALSE );
    }

    mbsFaces.reset4Loop();
    while( mbsFaces.setNext4Loop() ) {
        mbsFaces.getEntity()->isVisited( rg_FALSE );
    }
}



void AugmentedBetaComplex::divideMBSShellsIntoBodies( rg_dList< MBSShell* > mbsShells, ManifoldizedBetaShape& manifoldizedBS )
{
    rg_INT maxIndexOfAugmentedBCVertex = getMaxIndexOfAugmentedBCVertices();
    
    rg_INT i=0;
    rg_INT* mbsVerticesBelongToThisBody = new rg_INT[ maxIndexOfAugmentedBCVertex + 1 ]; //MBSVertex or AugmentedBCVertex의 ID를 넣으면, body의 ID가 나온다.
    for ( i=0; i < maxIndexOfAugmentedBCVertex + 1; i++ ) {
        mbsVerticesBelongToThisBody[i] = -1;
    }

    rg_INT numOfBodies = markToVerticesWhichBodyTheyBelongAndReturnNumberOfBodies( mbsVerticesBelongToThisBody );

    MBSBody** ptToMBSBody = rg_NULL;
    if ( numOfBodies != 0 ) {
        ptToMBSBody = new MBSBody*[ numOfBodies ];
    }

    for ( i=0; i < numOfBodies; i++ ) {
        MBSBody* currBody = new MBSBody(i);

        manifoldizedBS.addBody( currBody );
        ptToMBSBody[i] = currBody;
    }

    rg_INT IndexOfSingularBody = numOfBodies;
    mbsShells.reset4Loop();
    while ( mbsShells.setNext4Loop() ) {
        MBSShell* currShell = mbsShells.getEntity();

        rg_INT indexOfBody = mbsVerticesBelongToThisBody[ currShell->getVertices()->getFirstEntity()->getID() ];
        if ( indexOfBody == -1 ) { //currShell bounds singular-body
            //add body            
            MBSBody* tempBody = new MBSBody(IndexOfSingularBody);
            IndexOfSingularBody++;
            manifoldizedBS.addBody( tempBody );
            tempBody->setExteriorShell( currShell );
            currShell->setBody( tempBody );
        }
        else {
            ptToMBSBody[ indexOfBody ]->addInteriorShell( currShell );
            currShell->setBody( ptToMBSBody[ indexOfBody ] );
        }
    }

    
    
    if ( ptToMBSBody != rg_NULL )
        delete [] ptToMBSBody;

    delete [] mbsVerticesBelongToThisBody;
}



rg_INT AugmentedBetaComplex::markToVerticesWhichBodyTheyBelongAndReturnNumberOfBodies( rg_INT* mbsVerticesBelongToThisBody )
{
    rg_INT indexOfBody = 0;

    m_augmentedBCCells.reset4Loop();
    while ( m_augmentedBCCells.setNext4Loop() ) {
        AugmentedBCCell* currCell = m_augmentedBCCells.getpEntity();

        if ( currCell->isVisited() ) {
            continue;
        }

        rg_dList< AugmentedBCCell* > cellStack;
        cellStack.pushBack( currCell );

        while ( cellStack.getSize() > 0 ) {
            AugmentedBCCell* popedCell = cellStack.popBack();

            if ( popedCell->isVisited() ) {
                continue;
            }

            popedCell->isVisited( rg_TRUE );

            rg_INT i=0;
            for ( i=0; i < EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
                mbsVerticesBelongToThisBody[ popedCell->getVertex(i)->getID() ]
                                                                        = indexOfBody;
            }

            
            for ( i=0; i < EIWDS_NUM_FACE_ON_CELL; i++ ) {
                AugmentedBCFace* currFace = (AugmentedBCFace*)( popedCell->getFace(i) );

                if ( currFace->getBoundingState() != INTERIOR_SIMPLEX ) {
                    continue;
                }

                AugmentedBCCell* incidentCell = rg_NULL;
                if ( currFace->getLeftCell() == popedCell ) {
                    incidentCell = (AugmentedBCCell*) currFace->getRightCell();
                }
                else { //if ( currFace->getRightCell() == popedCell )
                    incidentCell = (AugmentedBCCell*) currFace->getLeftCell();
                }

                if ( incidentCell->isVisited() == rg_FALSE ) {
                    cellStack.pushBack( incidentCell );
                }
            }
        }

        indexOfBody++;
    }




    m_artificialAugmentedBCCells.reset4Loop();
    while ( m_artificialAugmentedBCCells.setNext4Loop() ) {
        AugmentedBCCell* currCell = m_artificialAugmentedBCCells.getpEntity();

        if ( currCell->isVisited() ) {
            continue;
        }

        rg_dList< AugmentedBCCell* > cellStack;
        cellStack.pushBack( currCell );

        while ( cellStack.getSize() > 0 ) {
            AugmentedBCCell* popedCell = cellStack.popBack();

            if ( popedCell->isVisited() ) {
                continue;
            }

            popedCell->isVisited( rg_TRUE );

            rg_INT i=0;
            for ( i=0; i < EIWDS_NUM_VERTEX_ON_CELL; i++ ) {
                mbsVerticesBelongToThisBody[ popedCell->getVertex(i)->getID() ]
                                                                        = indexOfBody;
            }

            
            for ( i=0; i < EIWDS_NUM_FACE_ON_CELL; i++ ) {
                AugmentedBCFace* currFace = (AugmentedBCFace*)( popedCell->getFace(i) );

                if ( currFace->getBoundingState() != INTERIOR_SIMPLEX ) {
                    continue;
                }

                AugmentedBCCell* incidentCell = rg_NULL;
                if ( currFace->getLeftCell() == popedCell ) {
                    incidentCell = (AugmentedBCCell*) currFace->getRightCell();
                }
                else { //if ( currFace->getRightCell() == popedCell )
                    incidentCell = (AugmentedBCCell*) currFace->getLeftCell();
                }

                if ( incidentCell->isVisited() == rg_FALSE ) {
                    cellStack.pushBack( incidentCell );
                }
            }
        }

        indexOfBody++;
    }



    m_augmentedBCCells.reset4Loop();
    while ( m_augmentedBCCells.setNext4Loop() ) {
        m_augmentedBCCells.getpEntity()->isVisited( rg_FALSE );
    }
    m_artificialAugmentedBCCells.reset4Loop();
    while ( m_artificialAugmentedBCCells.setNext4Loop() ) {
        m_artificialAugmentedBCCells.getpEntity()->isVisited( rg_FALSE );
    }

    return indexOfBody;
}


void AugmentedBetaComplex::divideExteriorShellFromInteriorShells( ManifoldizedBetaShape& manifoldizedBS )
{
    if ( m_betaUniverse->getVertexList()->getSize() == 0 ) {
        return;
    }


    rg_dList< BetaVertex* > verticesBoundExteriorShell;
    searchVerticesWhichBoundExteriorShell( verticesBoundExteriorShell );
    
    set< BetaVertex* > betaVerticesToTouchExteriorShell;

    verticesBoundExteriorShell.reset4Loop();
    while ( verticesBoundExteriorShell.setNext4Loop() ) {
        BetaVertex* currVertex = verticesBoundExteriorShell.getEntity();

        rg_INT numOfExteriorCones = currVertex->countNumGroupsOfFaceConnectedExtraneousCellsInBetaComplex( m_betaValue );

        if ( numOfExteriorCones == 1) {
            betaVerticesToTouchExteriorShell.insert(currVertex);
        }
    }


    rg_dList< MBSBody* >* bodies = manifoldizedBS.getBodies();
    bodies->reset4Loop();
    while ( bodies->setNext4Loop() ) {
        MBSBody* currBody = bodies->getEntity();

        rg_dList< MBSShell* >* shellsInCurrBody = currBody->getInteriorShells();
        
        if ( shellsInCurrBody->getSize() == 0 ) { //currBody is a singular body
            //  do nothing.
        }
        else if ( shellsInCurrBody->getSize() == 1 ) { //interiorShell이 1개 이상은 있다. exteriorShell은 singular body인 경우를 제외하고는 항상 rg_NULL이다.
            currBody->setExteriorShell( shellsInCurrBody->getFirstEntity() );
            currBody->getInteriorShells()->removeAll();
        }
        else {
            rg_dNode< MBSShell* >* exteriorShellNode          = rg_NULL;
            rg_dNode< MBSShell* >* shellNodeWithMaxNumOfFaces = rg_NULL;
            rg_INT maxNumOfFaces = 0;

            rg_dNode< MBSShell* >* currNode = shellsInCurrBody->getFirstpNode();
            for ( rg_INT i=0; i<shellsInCurrBody->getSize(); i++, currNode=currNode->getNext() ) {
                MBSShell* currShell = currNode->getEntity();

                rg_dList< MBSVertex* >* verticesOnCurrShell = currShell->getVertices();
                verticesOnCurrShell->reset4Loop();
                while ( verticesOnCurrShell->setNext4Loop() ) {
                    MBSVertex* currVtx = verticesOnCurrShell->getEntity();

                    if ( betaVerticesToTouchExteriorShell.find( currVtx->getOriginalBetaVertex() ) != betaVerticesToTouchExteriorShell.end() ) {
                        exteriorShellNode = currNode;
                        break;
                    }
                }

                rg_INT numOfFacesOnCurrShell = currShell->getFaces()->getSize();
                if ( numOfFacesOnCurrShell > maxNumOfFaces ) {
                    shellNodeWithMaxNumOfFaces = currNode;
                    maxNumOfFaces              = numOfFacesOnCurrShell;
                }
            }


            if ( exteriorShellNode != rg_NULL ) {
                currBody->setExteriorShell( exteriorShellNode->getEntity() );
                shellsInCurrBody->kill( exteriorShellNode );        
            }
            else {
                currBody->setExteriorShell( shellNodeWithMaxNumOfFaces->getEntity() );
                shellsInCurrBody->kill( shellNodeWithMaxNumOfFaces );
            }
        }
    }



    //  The following cods are fixed by Youngsong Cho 2011-04-23.
    //
    //if ( m_betaUniverse->getVertexList()->getSize() == 0 ) {
    //    return;
    //}

    //rg_INT  maxIndexOfMBSVertex = getMaxIndexOfAugmentedBCVertices();
    //rg_INT* isThisVertexTouchedOnExteriorShell = new rg_INT[ maxIndexOfMBSVertex + 1 ];
   
    ////status of "isThisVertexTouchedOnExteriorShell"
    //const rg_INT TOUCHES_EXTERIOR_SHELL         = 1;
    //const rg_INT DOES_NOT_TOUCH_EXTERIOR_SHELL  = 0;
    //const rg_INT UNKNWON                        = -1;   //Because, This vertex is the nm-vertex of TT^D_V
    //

    //for ( rg_INT i=0; i < maxIndexOfMBSVertex + 1; i++) {
    //    isThisVertexTouchedOnExteriorShell[i] = DOES_NOT_TOUCH_EXTERIOR_SHELL;
    //}

    //rg_dList< BetaVertex* > verticesBoundExteriorShell;
    //searchVerticesWhichBoundExteriorShell( verticesBoundExteriorShell );

    //
    //verticesBoundExteriorShell.reset4Loop();
    //while ( verticesBoundExteriorShell.setNext4Loop() ) {
    //    BetaVertex* currVertex = verticesBoundExteriorShell.getEntity();

    //    rg_INT numOfExteriorCones
    //                = currVertex->countNumGroupsOfFaceConnectedExtraneousCellsInBetaComplex( m_betaValue );

    //    if ( numOfExteriorCones == 1) {
    //        isThisVertexTouchedOnExteriorShell[ currVertex->getID() ] = TOUCHES_EXTERIOR_SHELL;
    //    }
    //    else if ( numOfExteriorCones > 1 ) {
    //        isThisVertexTouchedOnExteriorShell[ currVertex->getID() ] = UNKNWON;
    //    }
    //}


    //rg_dList< MBSBody* >* bodies = manifoldizedBS.getBodies();
    //bodies->reset4Loop();
    //while ( bodies->setNext4Loop() ) {
    //    MBSBody* currBody = bodies->getEntity();

    //    rg_dList< MBSShell* >* currInteriorShells = currBody->getInteriorShells();
    //    
    //    if ( currInteriorShells->getSize() == 0 ) { //currBody is a singular body
    //        continue;
    //    }

    //    if ( currInteriorShells->getSize() == 1 ) { //interiorShell이 1개 이상은 있다. exteriorShell은 singular body인 경우를 제외하고는 항상 rg_NULL이다.
    //        currBody->setExteriorShell( currInteriorShells->getFirstEntity() );
    //        currBody->getInteriorShells()->removeAll();

    //        continue;
    //    }

    //    
    //    
    //    rg_dNode< MBSShell* >* startNode = currInteriorShells->getFirstpNode();
    //    rg_dNode< MBSShell* >* currNode  = startNode;

    //    rg_dNode< MBSShell* >* exteriorShellNode             = rg_NULL;
    //    rg_dNode< MBSShell* >* exteriorShellNodeByNumOfFaces = rg_NULL;
    //    rg_INT maxNumOfFaces = 0;
    //
    //    do {
    //        MBSShell* currShell = currNode->getEntity();
    //        
    //        rg_dList< MBSVertex* >* verticesOnCurrShell = currShell->getVertices();
    //        verticesOnCurrShell->reset4Loop();
    //        while ( verticesOnCurrShell->setNext4Loop() ) {
    //            MBSVertex* currVertexonCurrShell = verticesOnCurrShell->getEntity();

    //            if ( isThisVertexTouchedOnExteriorShell[ currVertexonCurrShell->getID() ]
    //                         == TOUCHES_EXTERIOR_SHELL ) {
    //                exteriorShellNode = currNode;
    //                break;
    //            }
    //        }
    //        
    //        rg_INT numOfFacesOnCurrShell = currShell->getFaces()->getSize();
    //        if ( numOfFacesOnCurrShell > maxNumOfFaces ) {
    //            exteriorShellNodeByNumOfFaces = currNode;
    //            maxNumOfFaces                 = numOfFacesOnCurrShell;
    //        }
    //        

    //        currNode = currNode->getNext();
    //	    
    //    } while( currNode != startNode );

    //    if ( exteriorShellNode != rg_NULL ) {
    //        currBody->setExteriorShell( exteriorShellNode->getEntity() );
    //        currInteriorShells->kill( exteriorShellNode );        
    //    }
    //    else {
    //        currBody->setExteriorShell( exteriorShellNodeByNumOfFaces->getEntity() );
    //        currInteriorShells->kill( exteriorShellNodeByNumOfFaces );
    //    }
    //}

    //delete [] isThisVertexTouchedOnExteriorShell;
}



void AugmentedBetaComplex::searchVerticesWhichBoundExteriorShell( rg_dList< BetaVertex* >& verticesBoundExteriorShell )
{
    rg_dList< BetaCell >* betaCells = m_betaUniverse->getCellList();
    BetaCell* firstVirtualCell = rg_NULL;
    betaCells->reset4Loop();
    while ( betaCells->setNext4Loop() ) {
        BetaCell* currCell = betaCells->getpEntity();

        if ( currCell->isVirtual() ) {
            firstVirtualCell = currCell;
            break;
        }
    }

    rg_dList< BetaCell* > cellStack;
    cellStack.pushBack( firstVirtualCell );

    while ( cellStack.getSize() > 0 ) {
        BetaCell* popedCell = cellStack.popBack();

        if ( popedCell->isVisited() ) {
            continue;
        }

        popedCell->isVisited( rg_TRUE );
        for ( rg_INT i_vtx=0; i_vtx<EIWDS_NUM_VERTEX_ON_CELL; i_vtx++ ) {
            BetaVertex* currVertex = popedCell->getVertex( i_vtx );

            if ( currVertex->isVisited() ) {
                continue;
            }

            rg_INT currVertexBoundingState = currVertex->getBoundingState( m_betaValue );
            if (    currVertexBoundingState == EXTRANEOUS_SIMPLEX ) {
                continue;
            }

            currVertex->isVisited( rg_TRUE );
            verticesBoundExteriorShell.add( currVertex );

        }

        for ( rg_INT i_face=0; i_face<EIWDS_NUM_FACE_ON_CELL; i_face++ ) {
            BetaFace* currFace = popedCell->getFace( i_face );

            rg_INT currFaceBoundingState = currFace->getBoundingState( m_betaValue );
            if (    currFaceBoundingState != EXTRANEOUS_SIMPLEX ) {
                continue;
            }

            BetaCell* nextCell = rg_NULL;
            if ( currFace->getLeftCell() == popedCell ) {
                nextCell = currFace->getRightCell();
            }
            else {
                nextCell = currFace->getLeftCell();
            }

            if (    nextCell->getBoundingState( m_betaValue ) ==  EXTRANEOUS_SIMPLEX
                 && nextCell->isVisited() == rg_FALSE ) {
                cellStack.pushBack( nextCell );
            }
        }
    }

    verticesBoundExteriorShell.reset4Loop();
    while ( verticesBoundExteriorShell.setNext4Loop() ) {
        verticesBoundExteriorShell.getEntity()->isVisited( rg_FALSE );
    }

    betaCells->reset4Loop();
    while ( betaCells->setNext4Loop() ) {
        betaCells->getpEntity()->isVisited( rg_FALSE );
    }
}




/*
rg_REAL AugmentedBetaComplex::findUpperValueOfBoundingState( const rg_BetaSpan& betaSpan, rg_INT givenBoundingState)
{
	rg_REAL upperValueOfBoundingState = REAL_MINUS_INFINITE;

	rg_INT  numOfBetaInterval = betaSpan.getNumBetaInterval();
	rg_INT* boundingStateArray = betaSpan.getBoundingStates();
	rg_REAL* betaIntervalArray = betaSpan.getBetaIntervals();

	for (rg_INT i=0; i < numOfBetaInterval; i++ ) {
		if ( boundingStateArray[i] == givenBoundingState )
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
*/



void AugmentedBetaComplex::findAugRegularFacesIncidentToAugEdge( AugmentedBCEdge* givenAugEdge, rg_dList< AugmentedBCFace* >& incidentFaces )
{
    AugmentedBCFace* startFace = (AugmentedBCFace*)givenAugEdge->getFirstFace();
	AugmentedBCFace* currFace  = startFace;

	while ( givenAugEdge->getCCWNextFace(currFace) != rg_NULL )	{
        currFace = (AugmentedBCFace*)givenAugEdge->getCCWNextFace(currFace);
	}

	if ( currFace->getBoundingState() == REGULAR_SIMPLEX )
		incidentFaces.add( currFace );

	currFace  = startFace;
	while ( givenAugEdge->getCWNextFace(currFace) != rg_NULL ) {
        currFace = (AugmentedBCFace*)givenAugEdge->getCWNextFace(currFace);
	}

	if ( currFace->getBoundingState() == REGULAR_SIMPLEX )
		incidentFaces.add( currFace );
}


void AugmentedBetaComplex::groupFacesIncidentToVertexByEdge(AugmentedBCVertex* givenVertex, 
                                                                rg_dList< rg_dList<AugmentedBCFace*> >& faceGroupList)
{
    rg_dList<AugmentedBCFace*> faceListIncidentToVertex;
    findAugRegularFacesIncidentToAugVertex(givenVertex, faceListIncidentToVertex);

	AugmentedBCFace* currFace = rg_NULL;
	faceListIncidentToVertex.reset4Loop();
	while ( faceListIncidentToVertex.setNext4Loop() )
	{
		currFace = faceListIncidentToVertex.getEntity();

		if ( currFace->isVisited() )
			continue;

        rg_dList<AugmentedBCFace*>* ptrRegularFaceGroup = faceGroupList.add( rg_dList<AugmentedBCFace*>() );

        currFace->isVisited( rg_TRUE );
		BetaFace* startFace = currFace;

        do {
            ptrRegularFaceGroup->add( currFace );


            BetaEdge* ccwEdgeAtGivenVertex = rg_NULL;

            BetaEdge** boundingEdge = currFace->getEdges();
            rg_INT j;
            for( j=0; j<EIWDS_NUM_EDGE_ON_FACE; j++ )  {
                if ( currFace->getEdgeOrientation(j) == rg_TRUE )  {
                    if ( boundingEdge[j]->getEndVertex() == givenVertex ) {
                        ccwEdgeAtGivenVertex = boundingEdge[j];
						break;
                    }
                }
                else {
                    if ( boundingEdge[j]->getStartVertex() == givenVertex ) {
                        ccwEdgeAtGivenVertex = boundingEdge[j];
						break;
                    }
                }
            }

            
            rg_dList<AugmentedBCFace*> faceListIncidentToEdge;
            findAugRegularFacesIncidentToAugEdge((AugmentedBCEdge*)ccwEdgeAtGivenVertex, faceListIncidentToEdge);
        

            if ( currFace == faceListIncidentToEdge.getFirstEntity() )
                currFace = faceListIncidentToEdge.getLastEntity();
            else
                currFace = faceListIncidentToEdge.getFirstEntity();

            currFace->isVisited( rg_TRUE );
            
        } while ( currFace != startFace );
    }

	faceListIncidentToVertex.reset4Loop();
	while ( faceListIncidentToVertex.setNext4Loop() )
		faceListIncidentToVertex.getEntity()->isVisited( rg_FALSE );
}


void AugmentedBetaComplex::findAugRegularFacesIncidentToAugVertex( AugmentedBCVertex* givenAugVertex, rg_dList< AugmentedBCFace* >& incidentFaces )
{
	rg_dList<AugmentedBCCell*> cellsIncidentToVertex;
	findAugInteriorCellsIncidentToAugVertex( givenAugVertex,  cellsIncidentToVertex );

	cellsIncidentToVertex.reset4Loop();
	while ( cellsIncidentToVertex.setNext4Loop() )
	{
		AugmentedBCCell* currCell = cellsIncidentToVertex.getEntity();

		BetaFace** faceOnCell = currCell->getFaces();
		for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++)
		{
			AugmentedBCFace* currAugFace = (AugmentedBCFace*)faceOnCell[i];

			if ( currAugFace->isThere( (BetaVertex*)givenAugVertex ) )  
			{
				if ( currAugFace->getBoundingState() == REGULAR_SIMPLEX )
					incidentFaces.addWithoutSame( currAugFace );
			}
		}
	}
}


void AugmentedBetaComplex::findAugInteriorCellsIncidentToAugVertex( AugmentedBCVertex* givenAugVertex, rg_dList< AugmentedBCCell* >& cellsIncidentToVertex )
{
	//AugmentedBCFace* firstFace = (AugmentedBCFace*)givenAugVertex->getFirstEdge()->getFirstFace();

    AugmentedBCEdge* firstEdge = (AugmentedBCEdge*) givenAugVertex->getFirstEdge();
    if ( firstEdge == rg_NULL ) { //stand-alone vertex
        return; 
    }

    AugmentedBCFace* firstFace = (AugmentedBCFace* ) firstEdge->getFirstFace();
    if ( firstFace == rg_NULL ) { // givenAugVertex bounds an edge which is a stand-alone edge
        return;
    }


	AugmentedBCCell* firstCell = rg_NULL;
	if ( firstFace->getRightCell() != rg_NULL )
		firstCell = (AugmentedBCCell*)firstFace->getRightCell();
	else
		firstCell = (AugmentedBCCell*)firstFace->getLeftCell();

	
	cellsIncidentToVertex.add( (AugmentedBCCell*)firstCell );
    firstCell->isVisited(rg_TRUE);

    rg_dList< pair<BetaFace*, BetaCell*> > facesToBeTraced;

    pair< BetaFace*, BetaCell* > currFace;
    BetaFace** boundingFace = firstCell->getFaces();
    for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {

        if ( boundingFace[i]->isThere( (BetaVertex*)givenAugVertex ) )  {
            currFace.first  = boundingFace[i];

            if ( firstCell == boundingFace[i]->getRightCell() )  
                currFace.second = boundingFace[i]->getLeftCell();
            else  
                currFace.second = boundingFace[i]->getRightCell();

            facesToBeTraced.add( currFace );
        }
    }


    while ( facesToBeTraced.getSize() > 0 )  {
        currFace = facesToBeTraced.popFront();

        if ( currFace.second == rg_NULL )
            continue;

        if ( currFace.second->isVisited() == rg_FALSE )  {        
            cellsIncidentToVertex.add( (AugmentedBCCell*)currFace.second );
            currFace.second->isVisited(rg_TRUE);

            boundingFace = currFace.second->getFaces();
            for ( rg_INT i=0; i<EIWDS_NUM_FACE_ON_CELL; i++ )  {
                if ( boundingFace[i] == currFace.first )
                    continue;

                if ( boundingFace[i]->isThere( (BetaVertex*)givenAugVertex ) )  {

                    pair<BetaFace*, BetaCell*> insertingFace;
                    insertingFace.first  = boundingFace[i];
                    if ( currFace.second == boundingFace[i]->getRightCell() )  
                        insertingFace.second = boundingFace[i]->getLeftCell();
                    else  
                        insertingFace.second = boundingFace[i]->getRightCell();

                    facesToBeTraced.add( insertingFace );
                }
            }
        }
    }

	cellsIncidentToVertex.reset4Loop();
	while ( cellsIncidentToVertex.setNext4Loop() )
		cellsIncidentToVertex.getEntity()->isVisited( rg_FALSE );
	
}




// 
// void AugmentedBetaComplex::constructMBSWithNoGenus( BetaComplex*  betaComplex, 
//                                                        ManifoldizedBetaShape& manifoldizedBetaShape )
// {
//     double defaultUpperBoundOfBetaValue   = 6.0;
//     double defaultMinimumRangeOfBetaValue = 0.2;
// 
//     constructMBSWithNoGenus( betaComplex, defaultUpperBoundOfBetaValue, defaultMinimumRangeOfBetaValue, manifoldizedBetaShape);
// }
// 
// 
// void AugmentedBetaComplex::constructMBSWithNoGenus( BetaComplex*  betaComplex, 
//                                                        const rg_REAL& givenUpperBoundOfBetaValue,
//                                                        const rg_REAL& minimumRangeOfBetaValue,
//                                                        ManifoldizedBetaShape& manifoldizedBetaShape )
// {
//     constructAugmentedBetaComplex( *betaComplex );
//     extractBoundaryOfAugmentedBetaComplex( manifoldizedBetaShape );
// 
//     rg_INT initialNumOfShell    = manifoldizedBetaShape.computeNumOfShells();
//     rg_INT initialEulerEquation = manifoldizedBetaShape.getVertices()->getSize()
//                                 - manifoldizedBetaShape.getEdges()->getSize()
//                                 + manifoldizedBetaShape.getFaces()->getSize();
// 
//     if ( initialEulerEquation == 2 * initialNumOfShell ) //there is no genus
//     {
//         return;        
//     }
//     
//     rg_REAL lowerBoundOfBetaValue = betaComplex->getBetaValue();
//     rg_REAL upperBoundOfBetaValue = givenUpperBoundOfBetaValue;    
// 
//     while ( upperBoundOfBetaValue - lowerBoundOfBetaValue > minimumRangeOfBetaValue )
//     {
//         this->removeAll();
//         manifoldizedBetaShape.removeAll();
// 
//         rg_REAL middleBetaValue = ( lowerBoundOfBetaValue + upperBoundOfBetaValue ) / 2;
//         
//         betaComplex->setBetaValue( middleBetaValue );
//         constructAugmentedBetaComplex( *betaComplex );
//         extractBoundaryOfAugmentedBetaComplex( manifoldizedBetaShape );
// 
//         rg_INT numOfShell    = manifoldizedBetaShape.computeNumOfShells();
//         rg_INT eulerEquation = manifoldizedBetaShape.getVertices()->getSize()
//                              - manifoldizedBetaShape.getEdges()->getSize()
//                              + manifoldizedBetaShape.getFaces()->getSize();
// 
//         if ( eulerEquation == 2 * numOfShell ) //there is no genus
//         {
//             upperBoundOfBetaValue = middleBetaValue;        
//         }
//         else
//         {
//             lowerBoundOfBetaValue = middleBetaValue;
//         }
//     }
// 
//     rg_INT resultNumOfShell    = manifoldizedBetaShape.computeNumOfShells();
//     rg_INT resultEulerEquation = manifoldizedBetaShape.getVertices()->getSize()
//                                - manifoldizedBetaShape.getEdges()->getSize()
//                                + manifoldizedBetaShape.getFaces()->getSize();
// 
//     if ( resultEulerEquation == 2 * resultNumOfShell ) //there is no genus
//     {
//         return;        
//     }
//     else
//     {
//         betaComplex->setBetaValue( upperBoundOfBetaValue );
//         constructAugmentedBetaComplex( *betaComplex );
//         extractBoundaryOfAugmentedBetaComplex( manifoldizedBetaShape );
// 
//         return;
//     }
// }


void AugmentedBetaComplex::clearMyself()
{
    m_betaValue = -1;

	m_augmentedBCCells.removeAll();
	m_augmentedBCFaces.removeAll();
	m_augmentedBCEdges.removeAll();
	m_augmentedBCVertices.removeAll();

   	m_artificialAugmentedBCCells.removeAll();
	m_artificialAugmentedBCFaces.removeAll();
	m_artificialAugmentedBCEdges.removeAll();
	m_artificialAugmentedBCVertices.removeAll();	

    m_balls.removeAll();
}



rg_INT AugmentedBetaComplex::getNumVerticesOutputConeIsNotOne()
{
	rg_INT numVerticesOutputConeIsNotOne = 0;

	m_augmentedBCVertices.reset4Loop();
	while ( m_augmentedBCVertices.setNext4Loop() )
	{
		AugmentedBCVertex* currVertex = m_augmentedBCVertices.getpEntity();

		if ( currVertex->getBoundingState() != REGULAR_SIMPLEX )
			continue;

        if ( currVertex->getFirstEdge() == rg_NULL ) {
            continue;
        }

        if ( ((AugmentedBCEdge* )currVertex->getFirstEdge())->getBoundingState() == SINGULAR_SIMPLEX ) {
            continue;   //currVertex bounds SINGULAR edge only.
        }
        
        if ( currVertex->getFirstEdge()->getFirstFace() == rg_NULL ) {
            continue;
        }

        if ( ((AugmentedBCFace*) currVertex->getFirstEdge()->getFirstFace())->getBoundingState()
                       == SINGULAR_SIMPLEX ) {
            continue;   //currVertex bounds SINGULAR face only
        }
        

		rg_dList< rg_dList<AugmentedBCFace*> > faceGroupList;
		groupFacesIncidentToVertexByEdge( currVertex, faceGroupList);

		if ( faceGroupList.getSize() > 1 )
			numVerticesOutputConeIsNotOne++;		
	}
	

    m_artificialAugmentedBCVertices.reset4Loop();
	while ( m_artificialAugmentedBCVertices.setNext4Loop() )
	{
		AugmentedBCVertex* currVertex = m_artificialAugmentedBCVertices.getpEntity();

		if ( currVertex->getBoundingState() != REGULAR_SIMPLEX )
			continue;


        if ( ((AugmentedBCEdge* )currVertex->getFirstEdge())->getBoundingState() == SINGULAR_SIMPLEX ) {
            continue;   //currVertex bounds SINGULAR edge only.
        }
        
        if ( ((AugmentedBCFace*) currVertex->getFirstEdge()->getFirstFace())->getBoundingState()
                       == SINGULAR_SIMPLEX ) {
            continue;   //currVertex bounds SINGULAR face only
        }
        

		rg_dList< rg_dList<AugmentedBCFace*> > faceGroupList;
		groupFacesIncidentToVertexByEdge( currVertex, faceGroupList);

		if ( faceGroupList.getSize() > 1 )
			numVerticesOutputConeIsNotOne++;		
	}


	return numVerticesOutputConeIsNotOne;
}


rg_INT AugmentedBetaComplex::getNumVerticesDifferentNumConesInputAndTT_D_V()
{
	return r_numVerticesDifferentNumConesInputAndTT_D_V;
}




//betaUniverse에 있는 simplex들 중에 지우고 싶은 것은 extraneous로 설정했을 때,
//generate NMConf at vertex 를 하면 아래 query를 이용했을 경우 문제가 된다.
//그래서 QuasitriangulationForBeta에 있는 query를 가져와서 
//boundingState를 get하는 부분을 수정하였다.
void AugmentedBetaComplex::searchGroupsOfFaceConnectedExtraneousCellsInBetaComplex( AugmentedBCVertex* givenVertex, 
                                                                                   rg_dList< rg_dList<BetaCell*> >& groupsOfExtraneousCell ) const
{
    rg_dList<BetaCell*> incidentQuasiCells;
    givenVertex->searchCellsInWholeWorld( incidentQuasiCells );

    BetaCell* currCell = rg_NULL;    
    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        currCell = incidentQuasiCells.getEntity();

        if ( currCell->isVisited() ) {
            continue;
        }
    
        rg_INT currCellBoundingState = ((AugmentedBCCell*) currCell)->getBoundingState();
        if (    currCellBoundingState != EXTRANEOUS_SIMPLEX 
             && currCellBoundingState != FALSE_EXTERIOR_SIMPLEX ) {
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
                rg_INT boundingStateOfBoundingFace = ((AugmentedBCFace*)boundingFace[i])->getBoundingState();
                if (    boundingStateOfBoundingFace != EXTRANEOUS_SIMPLEX
                     && boundingStateOfBoundingFace != FALSE_EXTERIOR_SIMPLEX ) {
                    continue;
                }

                if ( boundingFace[i]->isThere( givenVertex ) ) {
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



void AugmentedBetaComplex::searchGroupsOfFaceConnectedInteriorCellsInBetaComplex( AugmentedBCVertex* givenVertex,  rg_dList< rg_dList<BetaCell*> >& groupsOfInteriorCell ) const
{
    rg_dList<BetaCell*> incidentQuasiCells;
    givenVertex->searchCellsInWholeWorld( incidentQuasiCells );

    BetaCell* currCell = rg_NULL;    
    incidentQuasiCells.reset4Loop();
    while ( incidentQuasiCells.setNext4Loop() ) {
        currCell = incidentQuasiCells.getEntity();

        if ( currCell->isVisited() ) {
            continue;
        }

        if ( ((AugmentedBCCell*)currCell)->getBoundingState() != INTERIOR_SIMPLEX ) {
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
                if ( ((AugmentedBCFace*) boundingFace[i])->getBoundingState() != INTERIOR_SIMPLEX ) {
                    continue;
                }

                if ( boundingFace[i]->isThere( givenVertex ) ) {
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





void AugmentedBetaComplex::searchFacesIncidentToEdgeInBetaComplex( BetaEdge* givenEdge, 
                                                                   rg_dList< BetaFace* >& incidentFacesInBetaComplex ) const
{
    rg_dList< BetaFace* > incidentFaces;
    givenEdge->searchFiniteFacesInWholeWorld( incidentFaces );

    incidentFaces.reset4Loop();
    while ( incidentFaces.setNext4Loop() )   {
        AugmentedBCFace* currFace = (AugmentedBCFace*) incidentFaces.getEntity();

        rg_INT currFaceBoundingState = currFace->getBoundingState();
        if (    currFaceBoundingState != EXTRANEOUS_SIMPLEX
             && currFaceBoundingState != FALSE_EXTERIOR_SIMPLEX ) {
            incidentFacesInBetaComplex.add( (BetaFace*) currFace );
        }
    }
}




void AugmentedBetaComplex::searchFacesIncidentToEdgeInBetaShape( AugmentedBCEdge*       givenEdge, 
                                                                 rg_dList< BetaFace* >& facesIncidentToEdgeInBetaShape ) const
{
    rg_dList<BetaFace*> finiteFaceList;
    givenEdge->searchFiniteFacesInWholeWorld( finiteFaceList );

    finiteFaceList.reset4Loop();
    while ( finiteFaceList.setNext4Loop() ) {
        AugmentedBCFace* currFace = (AugmentedBCFace*) finiteFaceList.getEntity();

        rg_INT state = currFace->getBoundingState();
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            facesIncidentToEdgeInBetaShape.add( currFace );
        }
    }
}

void AugmentedBetaComplex::sortFacesCCWOfEdgeByGeometricComparison( AugmentedBCEdge*       givenEdge, 
                                                                    rg_dList< BetaFace* >& facesIncidentToEdge ) const
{
    ///////////////////////////////////////////////////////////////////////////
    //  Some applications need the incident faces to be sorted.
    //
    //  1. generate plane for computing angle of beta-face.
    rg_Point3D startPt = givenEdge->getStartVertex()->getBall().getCenter();
    rg_Point3D endPt   = givenEdge->getEndVertex()->getBall().getCenter();

    rg_Point3D normal = endPt - startPt;
    normal.normalize();

    Plane basePlane;
    basePlane.definePlaneByNormalAndPassingPoint(normal, startPt);



    //  2. compute angle of beta-face.
    rg_INT numFaces  = facesIncidentToEdge.getSize();

    BetaFace**  faces        = facesIncidentToEdge.getArray();
    rg_Point3D* projectedVtx = new rg_Point3D[numFaces];
    rg_REAL*    angle        = new rg_REAL[numFaces];
    
    rg_INT i=0;
    for ( i=0; i<numFaces; i++ ) {
        BetaVertex* mateVtx = faces[i]->findMateVertexOfEdge( givenEdge );
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


    //  4. make facesIncidentToEdgeInBetaShape for sorted face.
    facesIncidentToEdge.removeAll();
    for ( i=0; i<numFaces; i++ ) {
        facesIncidentToEdge.add( faces[i] );
    }


    delete [] faces;
    delete [] projectedVtx;
    delete [] angle;
}








void AugmentedBetaComplex::printBetaCell(ofstream& fout, BetaCell* betaCell, const rg_REAL& beta) 
{
    fout << "c" << betaCell->getID();
    fout  << "(" << getBoundingState( betaCell->getBoundingState(beta) ) << ")\t";
    int i=0;
    for ( i=0; i<4; i++ ) {
	    fout << "v" << betaCell->getVertex(i)->getID() << "\t";
    }
    fout << endl;
    for ( i=0; i<4; i++ ) {
        fout << "\t";
        printBetaFace(fout, betaCell->getFace(i), beta);
    }

}

void AugmentedBetaComplex::printBetaEdge(ofstream& fout, BetaEdge* betaEdge, const rg_REAL& beta) 
{
    fout << "e" << betaEdge->getID();
    fout  << "(" << getBoundingState( betaEdge->getBoundingState(beta) ) << ")\t";
    
    fout << "v" << betaEdge->getStartVertex()->getID();
    fout  << "(" << getBoundingState( betaEdge->getStartVertex()->getBoundingState(beta) ) << ")\t";
    
    fout << "v" << betaEdge->getEndVertex()->getID();
    fout  << "(" << getBoundingState( betaEdge->getEndVertex()->getBoundingState(beta) ) << ")\t";
    
    fout << endl;
    
    fout << "firstFace: ";
    printBetaFace( fout, betaEdge->getFirstFace(), beta);
    
}

void AugmentedBetaComplex::printBetaFace(ofstream& fout, BetaFace* betaFace, const rg_REAL& beta) 
{
    fout << "f" << betaFace->getID();
    fout  << "(" << getBoundingState( betaFace->getBoundingState(beta) ) << ")\t";

    fout << "lc" << betaFace->getLeftCell()->getID();
    fout  << "(" << getBoundingState( betaFace->getLeftCell()->getBoundingState(beta) ) << ")\t";
    
    fout << "rc" <<  betaFace->getRightCell()->getID() << "(";
    fout << getBoundingState( betaFace->getRightCell()->getBoundingState(beta) ) << ")\t";

    for ( int i=0; i<3; i++ ) {
        BetaEdge* currEdge = betaFace->getEdge(i);
        fout << "e" << currEdge->getID();
        fout << "(" << getBoundingState( currEdge->getBoundingState(beta) ) << ")";
        if ( betaFace->getEdgeOrientation(i) ) {
            fout << "+ (v" << currEdge->getStartVertex()->getID() << " -> v";
            fout << currEdge->getEndVertex()->getID() << ")";
        }
        else {
            fout << "- (v" << currEdge->getEndVertex()->getID() << " -> v";
            fout << currEdge->getStartVertex()->getID() << ")";
        }
        fout << "\t";
    }
    fout << endl;
}

char AugmentedBetaComplex::getBoundingState(const rg_INT& boundingState) 
{
    char state = 'I';
    switch ( boundingState )  {
    case FALSE_EXTERIOR_SIMPLEX:
        state = 'F';
        break;
    case EXTRANEOUS_SIMPLEX:
        state = 'E';
        break;
    case SINGULAR_SIMPLEX:
        state = 'S';
        break;
    case REGULAR_SIMPLEX:
        state = 'R';
        break;
    case INTERIOR_SIMPLEX:
        state = 'I';
        break;
    }

    return state;
}



void AugmentedBetaComplex::setSingularEdgesToFalseExterior()
{
    m_augmentedBCEdges.reset4Loop();
    while ( m_augmentedBCEdges.setNext4Loop() ) {
        AugmentedBCEdge* currEdge = m_augmentedBCEdges.getpEntity();

        setSingularEdgeToFalseExterior( currEdge );
    }
}

void AugmentedBetaComplex::setSingularEdgeToFalseExterior( AugmentedBCEdge* currEdge )
{
    if ( currEdge->getBoundingState() == SINGULAR_SIMPLEX ) {
        currEdge->setBoundingState( FALSE_EXTERIOR_SIMPLEX );
        currEdge->setBetaSpan( rg_BetaSpan() );            

        AugmentedBCVertex* vertex[2] = { (AugmentedBCVertex*) currEdge->getStartVertex(), 
                                         (AugmentedBCVertex*) currEdge->getEndVertex() };
        for (int i=0; i<2; i++) {
            rg_dList< BetaCell* > incidentCells;
            vertex[i]->searchCellsInWholeWorld( incidentCells );

            rg_BOOL isThereOnlyExtraneous = rg_TRUE;
            incidentCells.reset4Loop();
            while ( incidentCells.setNext4Loop() )  {
                AugmentedBCCell* currIncidentCell = (AugmentedBCCell*) incidentCells.getEntity();

                rg_INT currIncidentCellBoundingState = currIncidentCell->getBoundingState();
                if (    currIncidentCellBoundingState != EXTRANEOUS_SIMPLEX
                     && currIncidentCellBoundingState != FALSE_EXTERIOR_SIMPLEX )   {
                    isThereOnlyExtraneous = rg_FALSE;
                    break;
                }
            }

            if ( isThereOnlyExtraneous )    {
                vertex[i]->setBoundingState( FALSE_EXTERIOR_SIMPLEX );
                vertex[i]->setBetaSpan( rg_BetaSpan() );
            }
        }
    }
}



void AugmentedBetaComplex::setSingularFacesToFalseExterior()
{
    m_augmentedBCFaces.reset4Loop();
    while ( m_augmentedBCFaces.setNext4Loop() ) {
        AugmentedBCFace* currFace = m_augmentedBCFaces.getpEntity();

        if ( currFace->isSingularEllipticFace( m_betaValue ) ) {
            currFace->setBoundingState( INTERIOR_SIMPLEX );
            fixBoundingStateOfBoundingEdgesOfFace( currFace );
            fixBoundingStateOfBoundingVerticesOfFace( currFace );
        }
        else if ( currFace->getBoundingState() == SINGULAR_SIMPLEX ) {
            currFace->setBoundingState( FALSE_EXTERIOR_SIMPLEX );

            BetaEdge** boundingEdges = currFace->getEdges();
            for ( int i=0; i<3; i++)    {
                rg_dList< BetaCell* > incidentCells;
                boundingEdges[i]->searchCellsInWholeWorld( incidentCells );

                rg_BOOL isThereOnlyExtraneousCell = rg_TRUE;
                incidentCells.reset4Loop();
                while ( incidentCells.setNext4Loop() )  {
                    AugmentedBCCell* currIncidentCell = (AugmentedBCCell*) incidentCells.getEntity();

                    rg_INT currIncidentCellBoundingSate = currIncidentCell->getBoundingState();
                    if (   currIncidentCellBoundingSate != EXTRANEOUS_SIMPLEX 
                        && currIncidentCellBoundingSate != FALSE_EXTERIOR_SIMPLEX ) {
                        isThereOnlyExtraneousCell = rg_FALSE;
                        break;
                    }
                }

                if ( isThereOnlyExtraneousCell )    {
                    ((AugmentedBCEdge*)boundingEdges[i])->setBoundingState( FALSE_EXTERIOR_SIMPLEX );
                }
            }

            rg_dList< BetaVertex* > boundingVertices;
            currFace->searchVerticesInWholeWorld( boundingVertices );

            boundingVertices.reset4Loop();
            while ( boundingVertices.setNext4Loop() )   {
                //BetaVertex* currBoundingVertex = boundingVertices.getEntity();
                AugmentedBCVertex* currBoundingVertex = (AugmentedBCVertex*)boundingVertices.getEntity();

                rg_dList< BetaCell* > incidentCells;
                currBoundingVertex->searchCellsInWholeWorld( incidentCells );

                rg_BOOL isThereOnlyExtraneous = rg_TRUE;
                incidentCells.reset4Loop();
                while ( incidentCells.setNext4Loop() )  {
                    AugmentedBCCell* currIncidentCell = (AugmentedBCCell*) incidentCells.getEntity();

                    rg_INT currIncidentCellBoundingState = currIncidentCell->getBoundingState();
                    if (    currIncidentCellBoundingState != EXTRANEOUS_SIMPLEX
                         && currIncidentCellBoundingState != FALSE_EXTERIOR_SIMPLEX )   {
                        isThereOnlyExtraneous = rg_FALSE;
                        break;
                    }
                }

                if ( isThereOnlyExtraneous )    {
                    //((AugmentedBCVertex*)currBoundingVertex)->setBoundingState(FALSE_EXTERIOR_SIMPLEX);
                    currBoundingVertex->setBoundingState(FALSE_EXTERIOR_SIMPLEX);
                }
            }
        }
    }
}


void AugmentedBetaComplex::setSingularVerticesToFalseExterior()
{
    m_augmentedBCVertices.reset4Loop();
    while ( m_augmentedBCVertices.setNext4Loop() ) {
        AugmentedBCVertex* currVertex = m_augmentedBCVertices.getpEntity();

        if ( currVertex->getBoundingState() == SINGULAR_SIMPLEX ) {
            currVertex->setBoundingState( FALSE_EXTERIOR_SIMPLEX );            
        }
    }
}
