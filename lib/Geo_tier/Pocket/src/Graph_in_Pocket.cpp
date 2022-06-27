#include "Graph_in_Pocket.h"

#include <float.h>

#include "rg_RelativeOp.h"

/////////////////////////////////////////////////////////////////////////////////
//
// Constructor
Graph_in_Pocket::Graph_in_Pocket()
{
}



Graph_in_Pocket::Graph_in_Pocket(rg_REAL** adjacencyMatrix, const rg_INT& n)  //n by n matrix
{
    setGraph(adjacencyMatrix, n);
}


Graph_in_Pocket::Graph_in_Pocket(rg_dList< GEdge >* edgeList, 
                rg_dList< GVertex >* vertexList)
{
    //m_EdgeList  = *EdgeList;
    //m_vertexList = *vertexList;

    m_edgeList.duplicateList(*edgeList);
    m_vertexList.duplicateList(*vertexList);
}

                
Graph_in_Pocket::Graph_in_Pocket(const Graph_in_Pocket& graph)
{    
    /*
    rg_dList< GVertex >* vertexList;
    vertexList = graph.getVertexList();
    rg_INT numberOfVertices = vertexList.getSize();

    rg_REAL** adjacencyMatrix = new rg_REAL*[numberOfVertices];
    for ( rg_INT i=0; i < numberOfVertices; i++)  {
        adjacencyMatrix[i] = new rg_REAL[numberOfVertices];
    }

    i = 0;
    vertexList.reset4Loop();
    while (vertexList.setNext4Loop())
    {
        GVertex* currVertex = vertexList.getFirstpEntity();
        rg_dList< GEdge* >* adjacentEdges = currVertex->getAdjacentEdges();
        rg_INT numberOfAdjacentEdges = adjacentEdges->size();

        for ( int j=0; j < numberOfAdjacentEdges; j++)
        {
            adjacencyMatrix[i][j] = (adjacentEdges+j).getCost();
        }

        i++;
    }

    setGraph(adjacencyMatrix, numberOfVertices); 
    */
}



Graph_in_Pocket::~Graph_in_Pocket()
{
}

/////////////////////////////////////////////////////////////////////////////////
//
// Get Functions
rg_dList< GEdge >*    Graph_in_Pocket::getEdgeList()
{
    return &m_edgeList;
}



rg_dList< GVertex >*   Graph_in_Pocket::getVertexList()
{
    return &m_vertexList;
}

/////////////////////////////////////////////////////////////////////////////////
//
// Set Functions
void Graph_in_Pocket::setEdgeList(rg_dList< GEdge >* edgeList)
{
    m_edgeList.duplicateList(*edgeList);
    
    //m_EdgeList = *EdgeList;
}



void Graph_in_Pocket::setVertexList(rg_dList< GVertex >* vertexList)
{
    m_vertexList.duplicateList(*vertexList);
    
    //m_vertexList = *vertexList;
}





GEdge* Graph_in_Pocket::addEdge( const GEdge& edge )
{
    return m_edgeList.add( edge );
}



GVertex* Graph_in_Pocket::addVertex( const GVertex& vertex)
{
    return m_vertexList.add( vertex );
}


/////////////////////////////////////////////////////////////////////////////////
//
// Operator Overloadings
//Graph& operator =(const Graph& graph)







void Graph_in_Pocket::setGraph(rg_REAL** adjacencyMatrix, const rg_INT& n) //n by n matrix
{
    rg_INT i = 0;
    for (i = 0; i < n; i++)    {
        GVertex tempVertex(i);
        m_vertexList.add(tempVertex);
    }

    rg_dNode< GVertex >* ithVertex = m_vertexList.getFirstpNode();
    rg_dNode< GVertex >* jthVertex = m_vertexList.getFirstpNode();

    for ( i = 0; i < n; i++)    {
        for (rg_INT j = (i+1) ; j < n; j++)   {
            
            if ( rg_GT(adjacencyMatrix[i][j], 0.0) )   {
                GEdge tempEdge( (i*n + j) );
                tempEdge.setStartVertex(ithVertex->getpEntity());
                tempEdge.setEndVertex(jthVertex->getpEntity());
                tempEdge.setCost( adjacencyMatrix[i][j] );

                GEdge* newEdge = m_edgeList.add(tempEdge);
                ithVertex->getEntity().addAdjacentEdge( newEdge );
                jthVertex->getEntity().addAdjacentEdge( newEdge );

                jthVertex = jthVertex->getNext();
            }

            
        }
        ithVertex = ithVertex->getNext();
    }
}







void Graph_in_Pocket::setUserDataForVerticesInGraph(void** vertexData, const rg_INT& n)    //n vertex data
{
    m_vertexList.reset4Loop();

    rg_INT i = 0;
    while( m_vertexList.setNext4Loop() )
    {
        GVertex* currVertex = m_vertexList.getpEntity();
        currVertex->setVertexData(vertexData[i]);
        i++;
    }
}






void Graph_in_Pocket::setUserDataForEdgesInGraph(void*** EdgeData, const rg_INT& n)     //n by n matrix
{
    m_edgeList.reset4Loop();

    while ( m_edgeList.setNext4Loop() )
    {
        GEdge* currEdge    = m_edgeList.getpEntity();
        rg_INT rowIndex    = currEdge->getStartVertex()->getID();
        rg_INT columnIndex = currEdge->getEndVertex()->getID();

        currEdge->setEdgeData(EdgeData[rowIndex][columnIndex]);
    }
}






void Graph_in_Pocket::convertEdgePathToVertexPath(     rg_dList< GEdge* >* EdgeList, 
                                        rg_dList< GVertex* >* vertexList)
{    
    GEdge* firstEdge = EdgeList->getFirstEntity();
    GEdge* secondEdge = EdgeList->getSecondEntity();
    GVertex* prevVertex = rg_NULL;

    if ( firstEdge->getStartVertex() == secondEdge->getStartVertex()
        || firstEdge->getStartVertex() == secondEdge->getEndVertex())
    {
        prevVertex = firstEdge->getEndVertex();
    }
    else if ( firstEdge->getEndVertex() == secondEdge->getStartVertex()
             || firstEdge->getEndVertex() == secondEdge->getEndVertex())
    {
        prevVertex = firstEdge->getStartVertex();
    }

    vertexList->add(prevVertex);

    EdgeList->reset4Loop();
    while (EdgeList->setNext4Loop())
    {
        GEdge* currEdge = EdgeList->getEntity();
        if ( prevVertex == currEdge->getStartVertex() )
        {
            prevVertex = currEdge->getEndVertex();
            vertexList->add( prevVertex );
        }
        else if ( prevVertex == currEdge->getEndVertex() )
        {
            prevVertex = currEdge->getStartVertex();
            vertexList->add( prevVertex );
        }
    }
}









void Graph_in_Pocket::convertEdgePathToVertexDataPath( rg_dList< GEdge* >* EdgeList, 
                                        rg_dList<void*>* vertexDataList)
{
    rg_dList< GVertex* >* vertexList = rg_NULL;
    convertEdgePathToVertexPath( EdgeList, vertexList);

    vertexList->reset4Loop();
    while (vertexList->setNext4Loop() )
    {
        GVertex* currVertex = vertexList->getEntity();
        vertexDataList->add( currVertex->getVertexData() );
    }
}




void Graph_in_Pocket::convertEdgePathToEdgeDataPath(  rg_dList< GEdge* >* EdgeList, 
                                        rg_dList<void*>* EdgeDataList)
{
    EdgeList->reset4Loop();
    while (EdgeList->setNext4Loop() )
    {
        GEdge* currEdge = EdgeList->getEntity();
        EdgeDataList->add( currEdge->getEdgeData() );
    }
}





void Graph_in_Pocket::convertEdgePathToIndexPath(    rg_dList< GEdge* >* EdgeList, 
                                        rg_dList<rg_INT>* indexList)
{
    rg_dList< GVertex* >* vertexList = rg_NULL;
    convertEdgePathToVertexPath( EdgeList, vertexList);

    vertexList->reset4Loop();
    while (vertexList->setNext4Loop() )
    {
        GVertex* currVertex = vertexList->getEntity();
        indexList->add( currVertex->getID() );
    }
}






/////////////////////////////////////////////////////////////////////////////////
//
// Dijkstra

void Graph_in_Pocket::runDijkstra(GVertex* source, 
                    GVertex* destination, 
                    rg_dList< GEdge* >* shortestPath)
{
    executeDijkstraAlgorithm(source, destination, shortestPath);
}





void Graph_in_Pocket::runDijkstra(void* sourceVOID, 
                    void* destinationVOID, 
                    rg_dList< GEdge* >* shortestPath)
{
    GVertex* source      = rg_NULL;
    GVertex* destination = rg_NULL;

    settingSourceAndDestination(sourceVOID, destinationVOID, 
                                    source, destination);
    executeDijkstraAlgorithm(source, destination, shortestPath);
}






void Graph_in_Pocket::runDijkstra(const rg_INT& sourceIndex, 
                    const rg_INT& destinationIndex, 
                    rg_dList< GEdge* >* shortestPath)
{
    GVertex* source      = rg_NULL;
    GVertex* destination = rg_NULL;

    settingSourceAndDestination(sourceIndex, destinationIndex, 
                                    source, destination);
    executeDijkstraAlgorithm(source, destination, shortestPath);
}




/*
void Graph::runDijkstra(GVertex* source)
{
    rg_dList< GVertex* > unvisitedVertexList;  //unvisited vertex list
    rg_INT maxIndexOfVertex = findMaxIndexOfVertex();
    
    rg_dNode< GVertex* >** vertexIndexArray = new rg_dNode< GVertex* >*[ maxIndexOfVertex + 1 ];
                                //unvisited vertex list의 Container를 pointing하는 array ( Vertex ID 기반)
                                //예를 들어 vertexIndexArray[17]에는 Vertex ID가 17인 GVertex* 를 containing 하고 있는
                                // unvisitedVertexList의 vertex 주소가 들어있다.

    initialize(source, &unvisitedVertexList, vertexIndexArray);

    calculateCostFromSourceToAllOtherVertices(source, 
                                            &unvisitedVertexList, vertexIndexArray);
    
    if (vertexIndexArray != rg_NULL)
    {
        delete [] vertexIndexArray;
    }
}
*/


void Graph_in_Pocket::settingSourceAndDestination(void * sourceVOID, void* destinationVOID, 
                                        GVertex*& source, GVertex*& destination)
{
    m_vertexList.reset4Loop();
    while( m_vertexList.setNext4Loop() )
    {
        GVertex* currVertex = m_vertexList.getpEntity();
        if ( currVertex->getVertexData() == sourceVOID)
        {
            source = currVertex;
        }

        if ( currVertex->getVertexData() == destinationVOID)
        {
            destination = currVertex;
        }
    }
}




void Graph_in_Pocket::settingSourceAndDestination(const rg_INT& sourceIndex, const rg_INT& destinationIndex, 
                                        GVertex*& source, GVertex*& destination)
{
    m_vertexList.reset4Loop();
    while( m_vertexList.setNext4Loop() )
    {
        GVertex* currVertex = m_vertexList.getpEntity();
        if ( currVertex->getID() == sourceIndex)
        {
            source = currVertex;
        }

        if ( currVertex->getID() == destinationIndex)
        {
            destination = currVertex;
        }
    }
}




void Graph_in_Pocket::executeDijkstraAlgorithm(GVertex* source, GVertex* destination, 
                                    rg_dList< GEdge* >* shortestPath)
{
	PriorityQueueForGVertex unvisitedVertices( m_vertexList.getSize() );	

    rg_INT maxIndexOfVertex = findMaxIndexOfVertex();
    PQNode** vertexIndexArray = new PQNode*[ maxIndexOfVertex + 1 ];
                                //unvisited vertex list의 Container를 pointing하는 array ( Vertex ID 기반)
                                //예를 들어 vertexIndexArray[17]에는 Vertex ID가 17인 GVertex* 를 containing 하고 있는
                                // unvisitedVertexList의 vertex 주소가 들어있다.

    initialize(source, &unvisitedVertices, vertexIndexArray);

    calculateCostFromSourceToDestination(source, destination, 
                                            &unvisitedVertices, vertexIndexArray);
    defineShortestPath(source, destination, shortestPath);

    if (vertexIndexArray != rg_NULL)
    {
        delete [] vertexIndexArray;
    }
}



rg_INT Graph_in_Pocket::findMaxIndexOfVertex()
{
    rg_INT maxIndex = -1;

    GVertex*    currGVertex = rg_NULL;
    rg_INT      currGVertexID = -1;
    m_vertexList.reset4Loop();
    while ( m_vertexList.setNext4Loop() )
    {
        currGVertex = m_vertexList.getpEntity();
        currGVertexID = currGVertex->getID();

        if ( currGVertexID > maxIndex )
            maxIndex = currGVertexID;
    }

    return maxIndex;
}


void Graph_in_Pocket::initialize(GVertex* source,
					   PriorityQueueForGVertex* unvisitedVertices, 
					   PQNode** vertexIndexArray)
{
    m_vertexList.reset4Loop();
    while( m_vertexList.setNext4Loop() )
    {
        GVertex* currVertex = m_vertexList.getpEntity();
        currVertex->setPreEdge(rg_NULL);
        
        PQNode* addressOfAddedVertex;
        if ( currVertex == source )
        {
			addressOfAddedVertex = unvisitedVertices->insertSourcePQNode( currVertex );
        }
        else
        {
            addressOfAddedVertex = unvisitedVertices->insertPQNodeWithoutBubbling( currVertex );
        }
        
        vertexIndexArray[currVertex->getID()] = addressOfAddedVertex;     
    }

	unvisitedVertices->setNumOfNodes( unvisitedVertices->getNumOfNodes() + 1);
}





void Graph_in_Pocket::calculateCostFromSourceToDestination(GVertex* source, GVertex* destination, 
                                                    PriorityQueueForGVertex* unvisitedVertices, 
                                                    PQNode** vertexIndexArray)
{
    while( unvisitedVertices->getNumOfNodes() != 0 )
    {
        GVertex* nearestVertex = unvisitedVertices->getMinPQNode()->getGVertex();
        if ( nearestVertex == rg_NULL)
            break;        
        
        rg_dList< GEdge* >* adjacentEdges = nearestVertex->getAdjacentEdges();
        
		adjacentEdges->reset4Loop();
		while ( adjacentEdges->setNext4Loop() )
		{
			GEdge* currAdjacentEdgeEntity = adjacentEdges->getEntity();

            PQNode* incidentVertexIterator = vertexIndexArray[currAdjacentEdgeEntity->getEndVertex()->getID()];
            if ( incidentVertexIterator == rg_NULL )	//visited vertex는 rg_NULL로 setting.
				continue;   //poped vertex는 다시 계산할 필요가 없다.			
            
			
			GVertex* incidentVertex = incidentVertexIterator->getGVertex();
			if ( incidentVertex == rg_NULL ) //poped vertex
				continue;
            
			rg_REAL accumCost = unvisitedVertices->getMinPQNode()->getKeyValue() + currAdjacentEdgeEntity->getCost();
			
			if ( incidentVertexIterator->getKeyValue() > accumCost )
			{
				unvisitedVertices->modifyKeyOfPQNode( accumCost, vertexIndexArray[ incidentVertex->getID() ], vertexIndexArray );
				incidentVertex->setPreEdge( currAdjacentEdgeEntity );
			}
        }
        
		vertexIndexArray[ nearestVertex->getID() ] =  rg_NULL;
		unvisitedVertices->killMinPQNode( vertexIndexArray );

        //We only interested in a shortest path between source to destination.
        if (nearestVertex == destination)
            break;
    }    
}




void Graph_in_Pocket::calculateCostFromSourceToAllOtherVertices(GVertex* source, 
														PriorityQueueForGVertex* unvisitedVertices, 
														PQNode** vertexIndexArray)
{
    while( unvisitedVertices->getNumOfNodes() != 0 )
    {
        GVertex* nearestVertex = unvisitedVertices->getMinPQNode()->getGVertex();
        if ( nearestVertex == rg_NULL)
            break;        
        
        rg_dList< GEdge* >* adjacentEdges = nearestVertex->getAdjacentEdges();
        
		adjacentEdges->reset4Loop();
		while ( adjacentEdges->setNext4Loop() )
		{
			GEdge* currAdjacentEdgeEntity = adjacentEdges->getEntity();

            PQNode* incidentVertexIterator = vertexIndexArray[currAdjacentEdgeEntity->getEndVertex()->getID()];
            if ( incidentVertexIterator == rg_NULL )	//visited vertex는 rg_NULL로 setting.
            {
                continue;   //poped vertex는 다시 계산할 필요가 없다.
			}
            
			
			GVertex* incidentVertex = incidentVertexIterator->getGVertex();
            
			rg_REAL accumCost = unvisitedVertices->getMinPQNode()->getKeyValue() + currAdjacentEdgeEntity->getCost();

			if ( incidentVertexIterator->getKeyValue() > accumCost )
			{
				unvisitedVertices->modifyKeyOfPQNode( accumCost, vertexIndexArray[ incidentVertex->getID() ], vertexIndexArray );
				incidentVertex->setPreEdge( currAdjacentEdgeEntity );
			}
			
        }
        
		vertexIndexArray[ nearestVertex->getID() ] =  rg_NULL;
		unvisitedVertices->killMinPQNode( vertexIndexArray );
    }    
}








void Graph_in_Pocket::defineShortestPath(GVertex* source, GVertex* destination,
                               rg_dList< GEdge* >* shortestPath)
{
    GVertex* currentVertex = destination;
	while( !(currentVertex == source) )
	{
        GEdge* preEdge = currentVertex->getPreEdge();

        if ( preEdge == rg_NULL )          //shortest path가 없다는 뜻.
        {
            return;
        }

        currentVertex = preEdge->getStartVertex();
                
	    shortestPath->addHead(preEdge);       
	}
}