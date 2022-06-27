// SimplexQueryingMachineInQTByKDTree.cpp: implementation of the SimplexQueryingMachineInQTByKDTree class.
//
//////////////////////////////////////////////////////////////////////

#include "SimplexQueryingMachineInQTByKDTree.h"
#include "BetaUniverse.h"
using namespace V::GeometryTier;


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

template<> int node_less<BetaCell*>::axis = 0;
template<> int node_less<BetaFace*>::axis = 0;
template<> int node_less<BetaEdge*>::axis = 0;
template<> int node_less<BetaVertex*>::axis = 0;


SimplexQueryingMachineInQTByKDTree::SimplexQueryingMachineInQTByKDTree()
{

}

// SimplexQueryingMachineInQTByKDTree::SimplexQueryingMachineInQTByKDTree( rg_dList<BetaCell>& cellsInBetaUniverse, rg_dList<BetaFace>& facesInBetaUniverse,   rg_dList<BetaEdge>& edgesInBetaUniverse, rg_dList<BetaVertex>& verticesInBetaUniverse )
// {
// 	makeKDTreeForSearchingSimplex(cellsInBetaUniverse, facesInBetaUniverse,
// 								  edgesInBetaUniverse, verticesInBetaUniverse);
// }

SimplexQueryingMachineInQTByKDTree::SimplexQueryingMachineInQTByKDTree( BetaUniverse* betaUniverse )
{
	if(betaUniverse != rg_NULL) {
		makeKDTreeForSearchingSimplex(betaUniverse);
	}
}


SimplexQueryingMachineInQTByKDTree::~SimplexQueryingMachineInQTByKDTree()
{

}




// void SimplexQueryingMachineInQTByKDTree::makeKDTreeForSearchingSimplex(rg_dList<BetaCell>&  cellsInBetaUniverse, rg_dList<BetaFace>& facesInBetaUniverse,
// 									   rg_dList<BetaEdge>& edgesInBetaUniverse, rg_dList<BetaVertex>& verticesInBetaUniverse)
// {
// 	makeKDTreeOfCells(cellsInBetaUniverse);
// 	makeKDTreeOfFaces(facesInBetaUniverse);
// 	makeKDTreeOfEdges(edgesInBetaUniverse);
// 	makeKDTreeOfVertices(verticesInBetaUniverse);
// }

void SimplexQueryingMachineInQTByKDTree::makeKDTreeForSearchingSimplex( BetaUniverse* betaUniverse )
{
	if(betaUniverse != rg_NULL) {
		makeKDTreeOfCells(betaUniverse->getCellList());
 		makeKDTreeOfFaces(betaUniverse->getFaceList());
 		makeKDTreeOfEdges(betaUniverse->getEdgeList());
 		makeKDTreeOfVertices(betaUniverse->getVertexList());
	}	
}


void SimplexQueryingMachineInQTByKDTree::makeKDTreeOfCells(rg_dList<BetaCell>*  cellsInBetaUniverse)
{
	rg_dList<BetaCell*>	cellList;
	rg_dList<double*>	valueList;

	getBetaCellsForKDTree(cellsInBetaUniverse, cellList, valueList);

	m_cellsInKDTree.setAllEntitiesFromList(cellList, valueList, 3);
}


void SimplexQueryingMachineInQTByKDTree::getBetaCellsForKDTree( rg_dList<BetaCell>*  cellsInBetaUniverse, rg_dList<BetaCell*>& cellList, rg_dList<double*>& valueList )
{
	cellList.removeAll();
	valueList.removeAll();

	cellsInBetaUniverse->reset4Loop();
	while(cellsInBetaUniverse->setNext4Loop()) {
		BetaCell* currCell = cellsInBetaUniverse->getpEntity();

		if ( currCell->isVirtual() ) {
            continue;
        }

		rg_BetaSpan currBetaSpan = currCell->getBetaSpan();
        double* values;
        values = new double[3];

		rg_REAL* betaIntervals = currBetaSpan.getBetaIntervals();
       
		values[0] = betaIntervals[1];
		values[1] = betaIntervals[1];
		values[2] = betaIntervals[1];
        
		cellList.add(currCell);
        valueList.add(values);
	}
}


void SimplexQueryingMachineInQTByKDTree::makeKDTreeOfFaces(rg_dList<BetaFace>* facesInBetaUniverse)
{
	rg_dList<BetaFace*>	faceList;
	rg_dList<double*>	valueList;

	getBetaFacesForKDTreeAndSetAnomalyCase(facesInBetaUniverse, faceList, valueList);

	m_facesInKDTree.setAllEntitiesFromList(faceList, valueList, 3);
}


void SimplexQueryingMachineInQTByKDTree::getBetaFacesForKDTreeAndSetAnomalyCase( rg_dList<BetaFace>* facesInBetaUniverse, rg_dList<BetaFace*>& faceList, rg_dList<double*>& valueList )
{
	faceList.removeAll();
	valueList.removeAll();

	facesInBetaUniverse->reset4Loop();
	while(facesInBetaUniverse->setNext4Loop()) {
		BetaFace* currFace = facesInBetaUniverse->getpEntity();

		if ( currFace->isVirtual() ) {
            continue;
        }

// after decide name of function, modify.
// 		if(currFace->isXXXXXX()) {
// 			m_ellipticFacesOfAnomalyCase.add(currFace);
// 
// 			continue;
// 		}


		if(currFace->getBetaSpan().getNumBetaInterval() == 5) {
			m_ellipticFacesOfAnomalyCase.add(currFace);
			continue;
		}
		else if(currFace->getBetaSpan().getNumBetaInterval() == 4){
			if(currFace->getBetaSpan().getBoundingState(1) == REGULAR_SIMPLEX) {
				m_ellipticFacesOfAnomalyCase.add(currFace);
				continue;
			}
		}


		rg_BetaSpan currBetaSpan = currFace->getBetaSpan();
        double* values;
        values = new double[3];

		rg_REAL* betaIntervals = currBetaSpan.getBetaIntervals();

		values[0] = betaIntervals[1];

		if(currBetaSpan.getBoundingState(1) == SINGULAR_SIMPLEX) {
			values[1] = betaIntervals[2];
			values[2] = betaIntervals[3];
		}
		else {
			values[1] = betaIntervals[1];
			values[2] = betaIntervals[2];
		}


		faceList.add(currFace);
        valueList.add(values);
	}
}


void SimplexQueryingMachineInQTByKDTree::makeKDTreeOfEdges(rg_dList<BetaEdge>* edgesInBetaUniverse)
{
	rg_dList<BetaEdge*>	edgeList;
	rg_dList<double*>	valueList;

	getBetaEdgesForKDTreeAndAddGateEdges(edgesInBetaUniverse, edgeList, valueList);

	m_edgesInKDTree.setAllEntitiesFromList(edgeList, valueList, 3);
}


void SimplexQueryingMachineInQTByKDTree::getBetaEdgesForKDTreeAndAddGateEdges( rg_dList<BetaEdge>* edgesInBetaUniverse, rg_dList<BetaEdge*>& edgeList, rg_dList<double*>& valueList )
{
	edgeList.removeAll();
	valueList.removeAll();

	edgesInBetaUniverse->reset4Loop();
	while(edgesInBetaUniverse->setNext4Loop()) {
		BetaEdge* currEdge = edgesInBetaUniverse->getpEntity();

		if ( currEdge->isVirtual() ) {
            continue;
        }

		if(currEdge->isGateToSmallWorlds()) {
			m_gateEedges.add(currEdge);

			continue;
		}

		rg_BetaSpan currBetaSpan = currEdge->getBetaSpan();
        double* values;
        values = new double[3];

		rg_REAL* betaIntervals = currBetaSpan.getBetaIntervals();

        values[0] = betaIntervals[1];

		if(currBetaSpan.getBoundingState(1) == SINGULAR_SIMPLEX) {
			values[1] = betaIntervals[2];
			values[2] = betaIntervals[3];
		}
		else {
			values[1] = betaIntervals[1];
			values[2] = betaIntervals[2];
		}


		edgeList.add(currEdge);
        valueList.add(values);
	}
}


void SimplexQueryingMachineInQTByKDTree::makeKDTreeOfVertices(rg_dList<BetaVertex>* verticesInBetaUniverse)
{
	rg_dList<BetaVertex*>	vertexList;
	rg_dList<double*>		valueList;

	getBetaVerticesForKDTree(verticesInBetaUniverse, vertexList, valueList);

	m_verticesInKDTree.setAllEntitiesFromList(vertexList, valueList, 3);
}


void SimplexQueryingMachineInQTByKDTree::getBetaVerticesForKDTree( rg_dList<BetaVertex>* verticesInBetaUniverse, rg_dList<BetaVertex*>& vertexList, rg_dList<double*>& valueList )
{
	vertexList.removeAll();
	valueList.removeAll();

	verticesInBetaUniverse->reset4Loop();
	while(verticesInBetaUniverse->setNext4Loop()) {
		BetaVertex* currVertex = verticesInBetaUniverse->getpEntity();

		if ( currVertex->isVirtual() ) {
            continue;
        }

		rg_BetaSpan currBetaSpan = currVertex->getBetaSpan();
        double* values;
        values = new double[3];

		rg_REAL* betaIntervals = currBetaSpan.getBetaIntervals();
        
		values[0] = betaIntervals[0];
		values[1] = betaIntervals[1];
		values[2] = betaIntervals[2];


		vertexList.add(currVertex);
        valueList.add(values);
	}
}



void SimplexQueryingMachineInQTByKDTree::huntCellsForGivenBoundingState( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaCell*>& betaCellList )
{
	if(boundingState == SINGULAR_SIMPLEX || boundingState == REGULAR_SIMPLEX) {
		return;
	}

	KDTreeInterval* range;
	// KDTreeInterval range[3];
	
	makeRangeForGivenBoundingState(beta, boundingState, range);

	m_cellsInKDTree.getEntitiesInGivenRange(range, m_cellsInKDTree.getRootNode(), betaCellList);

	delete[] range;
	
}


void SimplexQueryingMachineInQTByKDTree::huntFacesForGivenBoundingState( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaFace*>& betaFaceList )
{
	KDTreeInterval* range;
	
	makeRangeForGivenBoundingState(beta, boundingState, range);

	m_facesInKDTree.getEntitiesInGivenRange(range, m_facesInKDTree.getRootNode(), betaFaceList);

	m_ellipticFacesOfAnomalyCase.reset4Loop();
	while(m_ellipticFacesOfAnomalyCase.setNext4Loop()) {
		BetaFace* currFace = m_ellipticFacesOfAnomalyCase.getEntity();
		
		rg_INT state = currFace->getBoundingState(beta);        
        if ( state == boundingState) {
			betaFaceList.add( currFace );
        }
    }

	delete[] range;
}


void SimplexQueryingMachineInQTByKDTree::huntEdgesForGivenBoundingState( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaEdge*>& betaEdgeList )
{
	KDTreeInterval* range;
	
	makeRangeForGivenBoundingState(beta, boundingState, range);

	m_edgesInKDTree.getEntitiesInGivenRange(range, m_edgesInKDTree.getRootNode(), betaEdgeList);

	m_gateEedges.reset4Loop();
	while(m_gateEedges.setNext4Loop()) {
		BetaEdge* currEdge = m_gateEedges.getEntity();
		
		rg_INT state = currEdge->getBoundingState(beta);
		if ( state == boundingState) {
			betaEdgeList.add( currEdge );
        }
    }

	delete[] range;
}


void SimplexQueryingMachineInQTByKDTree::huntVerticesForGivenBoundingState( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaVertex*>& betaVertexList )
{
	if(boundingState == EXTERIOR_SIMPLEX) {
		return;
	}
	
	KDTreeInterval* range;
	
	makeRangeForGivenBoundingState(beta, boundingState, range);

	m_verticesInKDTree.getEntitiesInGivenRange(range, m_verticesInKDTree.getRootNode(), betaVertexList);

	delete[] range;
}


void SimplexQueryingMachineInQTByKDTree::huntCellsForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaCell*>& betaCellList )
{
	if(fromBoundingState == toBoundingState) {
		return;
	}

	KDTreeInterval* range;
	
	makeRangeForGivenChangedState(fromBeta, toBeta, fromBoundingState, toBoundingState, range);

	m_cellsInKDTree.getEntitiesInGivenRange(range, m_cellsInKDTree.getRootNode(), betaCellList);

	delete[] range;
}


void SimplexQueryingMachineInQTByKDTree::huntFacesForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaFace*>& betaFaceList )
{
	if(fromBoundingState == toBoundingState) {
		return;
	}

	KDTreeInterval* range;
	
	makeRangeForGivenChangedState(fromBeta, toBeta, fromBoundingState, toBoundingState, range);

	m_facesInKDTree.getEntitiesInGivenRange(range, m_facesInKDTree.getRootNode(), betaFaceList);

	m_ellipticFacesOfAnomalyCase.reset4Loop();
	while(m_ellipticFacesOfAnomalyCase.setNext4Loop()) {
		BetaFace* currFace = m_ellipticFacesOfAnomalyCase.getEntity();
		
		rg_INT stateForFromBeta = currFace->getBoundingState(fromBeta);        
		rg_INT stateForToBeta	= currFace->getBoundingState(toBeta);        
        if ( stateForFromBeta == fromBoundingState && stateForToBeta == toBoundingState) {
			betaFaceList.add( currFace );
        }
    }

	delete[] range;
}


void SimplexQueryingMachineInQTByKDTree::huntEdgesForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaEdge*>& betaEdgeList )
{
	if(fromBoundingState == toBoundingState) {
		return;
	}

	KDTreeInterval* range;

	makeRangeForGivenChangedState(fromBeta, toBeta, fromBoundingState, toBoundingState, range);

	m_edgesInKDTree.getEntitiesInGivenRange(range, m_edgesInKDTree.getRootNode(), betaEdgeList);

	m_gateEedges.reset4Loop();
	while(m_gateEedges.setNext4Loop()) {
		BetaEdge* currEdge = m_gateEedges.getEntity();
		
		rg_INT stateForFromBeta = currEdge->getBoundingState(fromBeta);        
		rg_INT stateForToBeta	= currEdge->getBoundingState(toBeta);        
        if ( stateForFromBeta == fromBoundingState && stateForToBeta == toBoundingState) {
			betaEdgeList.add( currEdge );
        }
    }

	delete[] range;
}


void SimplexQueryingMachineInQTByKDTree::huntVerticesForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaVertex*>& betaVertexList )
{
	if(fromBoundingState == toBoundingState) {
		return;
	}

	KDTreeInterval* range;

	makeRangeForGivenChangedState(fromBeta, toBeta, fromBoundingState, toBoundingState, range);

	m_verticesInKDTree.getEntitiesInGivenRange(range, m_verticesInKDTree.getRootNode(), betaVertexList);

	delete[] range;
}


void SimplexQueryingMachineInQTByKDTree::makeRangeForGivenBoundingState( const rg_REAL& beta, const rg_INT& boundingState, KDTreeInterval*& range)
{
	range = new KDTreeInterval[3];

	switch(boundingState)
	{
		case EXTERIOR_SIMPLEX: // (b, inf]x(b, inf]x(b, inf]
			range[0].setLowerValue(beta);
			range[0].setStatusOfOpenClose(OPEN_CLOSE);    

			range[1].setLowerValue(beta);
			range[1].setStatusOfOpenClose(OPEN_CLOSE);    

			range[2].setLowerValue(beta);
			range[2].setStatusOfOpenClose(OPEN_CLOSE);    

			break;
		case SINGULAR_SIMPLEX: // [-inf, b]x(b, inf]x(b, inf]
			range[0].setUpperValue(beta);
			range[0].setStatusOfOpenClose(CLOSE_CLOSE);    

			range[1].setLowerValue(beta);
			range[1].setStatusOfOpenClose(OPEN_CLOSE);    

			range[2].setLowerValue(beta);
			range[2].setStatusOfOpenClose(OPEN_CLOSE);  

			break;
		case REGULAR_SIMPLEX: // [-inf, b]x[-inf, b]x(b, inf]
			range[0].setUpperValue(beta);
			range[0].setStatusOfOpenClose(CLOSE_CLOSE);    

			range[1].setUpperValue(beta);
			range[1].setStatusOfOpenClose(CLOSE_CLOSE);    

			range[2].setLowerValue(beta);
			range[2].setStatusOfOpenClose(OPEN_CLOSE);  

			break;
		case INTERIOR_SIMPLEX: // [-inf, b]x[-inf, b]x[-inf, b]
			range[0].setUpperValue(beta);
			range[0].setStatusOfOpenClose(CLOSE_CLOSE);    

			range[1].setUpperValue(beta);
			range[1].setStatusOfOpenClose(CLOSE_CLOSE);    

			range[2].setUpperValue(beta);
			range[2].setStatusOfOpenClose(CLOSE_CLOSE);  

			break;
		default:
			break;
	}
}


void SimplexQueryingMachineInQTByKDTree::makeRangeForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, KDTreeInterval*& range )
{
	range = new KDTreeInterval[3];

	switch(fromBoundingState)
	{
		case EXTERIOR_SIMPLEX:
			switch(toBoundingState)
			{
				//increasing beta
				case INTERIOR_SIMPLEX: // (b1, b2]x(b1, b2]x(b1, b2]
					range[0].setLowerValue(fromBeta);
					range[0].setUpperValue(toBeta);
					range[0].setStatusOfOpenClose(OPEN_CLOSE);    

					range[1].setLowerValue(fromBeta);
					range[1].setUpperValue(toBeta);
					range[1].setStatusOfOpenClose(OPEN_CLOSE);  
					
					range[2].setLowerValue(fromBeta);
					range[2].setUpperValue(toBeta);
					range[2].setStatusOfOpenClose(OPEN_CLOSE);  

					break;
				case REGULAR_SIMPLEX: // (b1, b2]x(b1, b2]x(b2, inf]
					range[0].setLowerValue(fromBeta);
					range[0].setUpperValue(toBeta);
					range[0].setStatusOfOpenClose(OPEN_CLOSE);    

					range[1].setLowerValue(fromBeta);
					range[1].setUpperValue(toBeta);
					range[1].setStatusOfOpenClose(OPEN_CLOSE);  
					
					range[2].setLowerValue(toBeta);					
					range[2].setStatusOfOpenClose(OPEN_CLOSE);    

					break;
				case SINGULAR_SIMPLEX: // (b1, b2]x(b2, inf]x(b2, inf]
					range[0].setLowerValue(fromBeta);
					range[0].setUpperValue(toBeta);
					range[0].setStatusOfOpenClose(OPEN_CLOSE);    

					range[1].setLowerValue(toBeta);					
					range[1].setStatusOfOpenClose(OPEN_CLOSE);  
					
					range[2].setLowerValue(toBeta);					
					range[2].setStatusOfOpenClose(OPEN_CLOSE);    

					break;
			}
			break;
		case SINGULAR_SIMPLEX:
			switch(toBoundingState)
			{
				//increasing beta
				case INTERIOR_SIMPLEX: // [-inf, b1]x(b1, b2]x(b1, b2]
					range[0].setUpperValue(fromBeta);
					range[0].setStatusOfOpenClose(CLOSE_CLOSE);    

					range[1].setLowerValue(fromBeta);
					range[1].setUpperValue(toBeta);
					range[1].setStatusOfOpenClose(OPEN_CLOSE);  
					
					range[2].setLowerValue(fromBeta);
					range[2].setUpperValue(toBeta);
					range[2].setStatusOfOpenClose(OPEN_CLOSE);     

					break;
				case REGULAR_SIMPLEX: // [-inf, b1]x(b1, b2]x(b2, inf]
					range[0].setUpperValue(fromBeta);
					range[0].setStatusOfOpenClose(CLOSE_CLOSE);    

					range[1].setLowerValue(fromBeta);					
					range[1].setUpperValue(toBeta);		
					range[1].setStatusOfOpenClose(OPEN_CLOSE);  
					
					range[2].setLowerValue(toBeta);					
					range[2].setStatusOfOpenClose(OPEN_CLOSE);      

					break;		
					
				//decreasing beta	
				case EXTERIOR_SIMPLEX: // (b3, b1]x(b1, inf]x(b1, inf]
					range[0].setLowerValue(toBeta);
					range[0].setUpperValue(fromBeta);
					range[0].setStatusOfOpenClose(OPEN_CLOSE);    

					range[1].setLowerValue(fromBeta);					
					range[1].setStatusOfOpenClose(OPEN_CLOSE);  
					
					range[2].setLowerValue(fromBeta);					
					range[2].setStatusOfOpenClose(OPEN_CLOSE);      

					break;			
			}

			break;
		case REGULAR_SIMPLEX:
			switch(toBoundingState)
			{
				//increasing beta
				case INTERIOR_SIMPLEX: // [-inf, b1]x[-inf, b1]x(b1, b2]
					range[0].setUpperValue(fromBeta);
					range[0].setStatusOfOpenClose(CLOSE_CLOSE);    

					range[1].setUpperValue(fromBeta);				
					range[1].setStatusOfOpenClose(OPEN_CLOSE);  
					
					range[2].setLowerValue(fromBeta);	
					range[2].setUpperValue(toBeta);
					range[2].setStatusOfOpenClose(OPEN_CLOSE);      

					break;			
					
				//decreasing beta	
				case EXTERIOR_SIMPLEX: // (b3, b1]x(b3, b1]x(b1, inf]
					range[0].setLowerValue(toBeta);
					range[0].setUpperValue(fromBeta);
					range[0].setStatusOfOpenClose(OPEN_CLOSE);    

					range[1].setLowerValue(toBeta);
					range[1].setUpperValue(fromBeta);
					range[1].setStatusOfOpenClose(OPEN_CLOSE);    
					
					range[2].setLowerValue(fromBeta);					
					range[2].setStatusOfOpenClose(OPEN_CLOSE);      

					break;	
				case SINGULAR_SIMPLEX: // [-inf, b3]x(b3, b1]x(b1, inf]					
					range[0].setUpperValue(toBeta);
					range[0].setStatusOfOpenClose(CLOSE_CLOSE);    

					range[1].setLowerValue(toBeta);
					range[1].setUpperValue(fromBeta);
					range[1].setStatusOfOpenClose(OPEN_CLOSE);    
					
					range[2].setLowerValue(fromBeta);					
					range[2].setStatusOfOpenClose(OPEN_CLOSE);      

					break;	
			}


			break;
		case INTERIOR_SIMPLEX:
			switch(toBoundingState)
			{
				//decreasing beta
				case EXTERIOR_SIMPLEX: // (b3, b1]x(b3, b1]x(b3, b1]
					range[0].setLowerValue(toBeta);
					range[0].setUpperValue(fromBeta);
					range[0].setStatusOfOpenClose(OPEN_CLOSE);    

					range[1].setLowerValue(toBeta);
					range[1].setUpperValue(fromBeta);
					range[1].setStatusOfOpenClose(OPEN_CLOSE);  
					
					range[2].setLowerValue(toBeta);
					range[2].setUpperValue(fromBeta);
					range[2].setStatusOfOpenClose(OPEN_CLOSE);  

					break;
				case REGULAR_SIMPLEX: // [-inf, b3]x[-inf, b3]x(b3, b1]					
					range[0].setUpperValue(toBeta);
					range[0].setStatusOfOpenClose(CLOSE_CLOSE);    
					
					range[1].setUpperValue(toBeta);
					range[1].setStatusOfOpenClose(CLOSE_CLOSE);  
					
					range[2].setLowerValue(toBeta);	
					range[2].setUpperValue(fromBeta);
					range[2].setStatusOfOpenClose(OPEN_CLOSE);    

					break;
				case SINGULAR_SIMPLEX: // [-inf, b3]x(b3, b1]x(b3, b1]					
					range[0].setUpperValue(toBeta);
					range[0].setStatusOfOpenClose(CLOSE_CLOSE);    

					range[1].setLowerValue(toBeta);	
					range[1].setUpperValue(fromBeta);
					range[1].setStatusOfOpenClose(OPEN_CLOSE);  
					
					range[2].setLowerValue(toBeta);	
					range[2].setUpperValue(fromBeta);
					range[2].setStatusOfOpenClose(OPEN_CLOSE);   

					break;
			}

			break;
		default:
			break;
	}
}


// HQ1
void SimplexQueryingMachineInQTByKDTree::huntFacesInSqueezedHull( rg_dList<BetaFace*>& betaFaceList )
{
	huntFacesForGivenBoundingState(DBL_MAX, REGULAR_SIMPLEX, betaFaceList);
}


void SimplexQueryingMachineInQTByKDTree::huntEdgesInSqueezedHull( rg_dList<BetaEdge*>& betaEdgeList )
{
	huntEdgesForGivenBoundingState(DBL_MAX, REGULAR_SIMPLEX, betaEdgeList);
}


void SimplexQueryingMachineInQTByKDTree::huntVerticesInSqueezedHull( rg_dList<BetaVertex*>& betaVertexList )
{
	huntVerticesForGivenBoundingState(DBL_MAX, REGULAR_SIMPLEX, betaVertexList);
}


// HQ2
void SimplexQueryingMachineInQTByKDTree::huntCellsInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList )
{
	KDTreeInterval* range = new KDTreeInterval[3];
	
	// [-inf, b]x[-inf, inf]x[-inf,inf]

 	range[0].setUpperValue(beta);    
    
	m_cellsInKDTree.getEntitiesInGivenRange(range, m_cellsInKDTree.getRootNode(), betaCellList);

	delete[] range;
}



void SimplexQueryingMachineInQTByKDTree::huntFacesInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList )
{
	KDTreeInterval* range = new KDTreeInterval[3];
	// betacomplex range
	// [-inf, b]x[-inf, inf]x[-inf,inf]

 	range[0].setUpperValue(beta);
        

	m_facesInKDTree.getEntitiesInGivenRange(range, m_facesInKDTree.getRootNode(), betaFaceList);

	m_ellipticFacesOfAnomalyCase.reset4Loop();
	while(m_ellipticFacesOfAnomalyCase.setNext4Loop()) {
		BetaFace* currFace = m_ellipticFacesOfAnomalyCase.getEntity();
		
		rg_INT state = currFace->getBoundingState(beta);        
        if ( state == EXTERIOR_SIMPLEX) {
			continue;
        }
		
		betaFaceList.add( currFace );
    }

	delete[] range;
}

void SimplexQueryingMachineInQTByKDTree::huntEdgesInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList )
{
	KDTreeInterval* range = new KDTreeInterval[3];

	// [-inf, b]x[-inf, inf]x[-inf,inf]

 	range[0].setUpperValue(beta);
        

	m_edgesInKDTree.getEntitiesInGivenRange(range, m_edgesInKDTree.getRootNode(), betaEdgeList);
	
	m_gateEedges.reset4Loop();
	while(m_gateEedges.setNext4Loop()) {
		BetaEdge* currEdge = m_gateEedges.getEntity();
		
		rg_INT state = currEdge->getBoundingState(beta);
		if ( state == EXTERIOR_SIMPLEX) {
			continue;
        }
        
        betaEdgeList.add( currEdge );
        
    }

	delete[] range;
}

void SimplexQueryingMachineInQTByKDTree::huntVerticesInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList )
{
	KDTreeInterval* range = new KDTreeInterval[3];

	// [-inf, b]x[-inf, inf]x[-inf,inf]
	
 	range[0].setUpperValue(beta);
    

	m_verticesInKDTree.getEntitiesInGivenRange(range, m_verticesInKDTree.getRootNode(), betaVertexList);

	delete[] range;
}



// HQ3
void SimplexQueryingMachineInQTByKDTree::huntFacesInBetaShapeInKDTree( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList )
{
	KDTreeInterval* range = new KDTreeInterval[3];

	// [-inf, b]x[-inf, inf]x(inf, b]

 	range[0].setUpperValue(beta);    
    
    range[2].setLowerValue(beta);
    range[2].setStatusOfOpenClose(OPEN_CLOSE);

	m_facesInKDTree.getEntitiesInGivenRange(range, m_facesInKDTree.getRootNode(), betaFaceList);

	m_ellipticFacesOfAnomalyCase.reset4Loop();
	while(m_ellipticFacesOfAnomalyCase.setNext4Loop()) {
		BetaFace* currFace = m_ellipticFacesOfAnomalyCase.getEntity();
		
		rg_INT state = currFace->getBoundingState(beta);        
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaFaceList.add( currFace );
        }
    }

	delete[] range;
}

void SimplexQueryingMachineInQTByKDTree::huntEdgesInBetaShapeInKDTree( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList )
{
	KDTreeInterval* range = new KDTreeInterval[3];

	// [-inf, b]x[-inf, inf]x(inf, b]

 	range[0].setUpperValue(beta);    
    
    range[2].setLowerValue(beta);
    range[2].setStatusOfOpenClose(OPEN_CLOSE);

	m_edgesInKDTree.getEntitiesInGivenRange(range, m_edgesInKDTree.getRootNode(), betaEdgeList);
	
	m_gateEedges.reset4Loop();
	while(m_gateEedges.setNext4Loop()) {
		BetaEdge* currEdge = m_gateEedges.getEntity();
		
		rg_INT state = currEdge->getBoundingState(beta);        
        if ( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
            betaEdgeList.add( currEdge );
        }
    }

	delete[] range;
}

void SimplexQueryingMachineInQTByKDTree::huntVerticesInBetaShapeInKDTree( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList )
{
	KDTreeInterval* range = new KDTreeInterval[3];

	// [-inf, b]x[-inf, inf]x(inf, b]

 	range[0].setUpperValue(beta);    
    
    range[2].setLowerValue(beta);
    range[2].setStatusOfOpenClose(OPEN_CLOSE);

	m_verticesInKDTree.getEntitiesInGivenRange(range, m_verticesInKDTree.getRootNode(), betaVertexList);

	delete[] range;
}


// HQ4
void SimplexQueryingMachineInQTByKDTree::huntFacesInRegularShapeInKDTree( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList )
{
	huntFacesForGivenBoundingState(beta, REGULAR_SIMPLEX, betaFaceList);
}


void SimplexQueryingMachineInQTByKDTree::huntEdgesInRegularShapeInKDTree( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList )
{
	huntEdgesForGivenBoundingState(beta, REGULAR_SIMPLEX, betaEdgeList);
}


void SimplexQueryingMachineInQTByKDTree::huntVerticesInRegularShapeInKDTree( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList )
{
	huntVerticesForGivenBoundingState(beta, REGULAR_SIMPLEX, betaVertexList);
}


// HQ5
void SimplexQueryingMachineInQTByKDTree::huntCellsInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaCell*>& betaCellListWithCellsForFromBeta )
{
	rg_dList<BetaCell*> marginalaCellsForE2I;
	rg_dList<BetaCell*> marginalaCellsForE2R;
	rg_dList<BetaCell*> marginalaCellsForE2S;

	huntCellsForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, INTERIOR_SIMPLEX, marginalaCellsForE2I);
	huntCellsForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, REGULAR_SIMPLEX, marginalaCellsForE2R);
	huntCellsForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, SINGULAR_SIMPLEX, marginalaCellsForE2S);

	betaCellListWithCellsForFromBeta.append(marginalaCellsForE2I);
	betaCellListWithCellsForFromBeta.append(marginalaCellsForE2R);
	betaCellListWithCellsForFromBeta.append(marginalaCellsForE2S);
}


void SimplexQueryingMachineInQTByKDTree::huntFacesInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaFace*>& betaFaceListWithFacesForFromBeta )
{
	rg_dList<BetaFace*> marginalaFacesForE2I;
	rg_dList<BetaFace*> marginalaFacesForE2R;
	rg_dList<BetaFace*> marginalaFacesForE2S;

	huntFacesForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, INTERIOR_SIMPLEX, marginalaFacesForE2I);
	huntFacesForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, REGULAR_SIMPLEX, marginalaFacesForE2R);
	huntFacesForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, SINGULAR_SIMPLEX, marginalaFacesForE2S);

	betaFaceListWithFacesForFromBeta.append(marginalaFacesForE2I);
	betaFaceListWithFacesForFromBeta.append(marginalaFacesForE2R);
	betaFaceListWithFacesForFromBeta.append(marginalaFacesForE2S);
}


void SimplexQueryingMachineInQTByKDTree::huntEdgesInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaEdge*>& betaEdgeListWithEdgesForFromBeta )
{
	rg_dList<BetaEdge*> marginalaEdgesForE2I;
	rg_dList<BetaEdge*> marginalaEdgesForE2R;
	rg_dList<BetaEdge*> marginalaEdgesForE2S;

	huntEdgesForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, INTERIOR_SIMPLEX, marginalaEdgesForE2I);
	huntEdgesForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, REGULAR_SIMPLEX, marginalaEdgesForE2R);
	huntEdgesForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, SINGULAR_SIMPLEX, marginalaEdgesForE2S);

	betaEdgeListWithEdgesForFromBeta.append(marginalaEdgesForE2I);
	betaEdgeListWithEdgesForFromBeta.append(marginalaEdgesForE2R);
	betaEdgeListWithEdgesForFromBeta.append(marginalaEdgesForE2S);
}


void SimplexQueryingMachineInQTByKDTree::huntVerticesInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaVertex*>& betaVertexListWithVerticesForFromBeta )
{
	rg_dList<BetaVertex*> marginalaVerticesForE2I;
	rg_dList<BetaVertex*> marginalaVerticesForE2R;
	rg_dList<BetaVertex*> marginalaVerticesForE2S;

	huntVerticesForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, INTERIOR_SIMPLEX, marginalaVerticesForE2I);
	huntVerticesForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, REGULAR_SIMPLEX, marginalaVerticesForE2R);
	huntVerticesForGivenChangedState(fromBeta, toBeta, EXTERIOR_SIMPLEX, SINGULAR_SIMPLEX, marginalaVerticesForE2S);

	betaVertexListWithVerticesForFromBeta.append(marginalaVerticesForE2I);
	betaVertexListWithVerticesForFromBeta.append(marginalaVerticesForE2R);
	betaVertexListWithVerticesForFromBeta.append(marginalaVerticesForE2S);
}