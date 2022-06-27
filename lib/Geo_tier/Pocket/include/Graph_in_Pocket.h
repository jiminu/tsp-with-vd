/////////////////////////////////////////////////////////////////////
//                           HISTORY
//
//  [20.03.31, CYSONG] 
//     This class is temporally renamed as Graph_in_Pocket from Graph.
//     Because the old name is collied when Graph in Mesh is moved to VUFS. 
//
//
/////////////////////////////////////////////////////////////////////

#ifndef _GRAPH_IN_POCKET_H_
#define _GRAPH_IN_POCKET_H_

#include "GEdge.h"
#include "GVertex.h"

#include "rg_dList.h"

#include "PriorityQueueForGVertex.h"

//a graph represented by adjacency list
class Graph_in_Pocket  
{
private:
	rg_dList< GVertex >   m_vertexList;
    rg_dList< GEdge >     m_edgeList;	


public:
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Constructor
	Graph_in_Pocket();
    Graph_in_Pocket(rg_REAL** adjacencyMatrix, const rg_INT& n);  //n by n matrix
    

    Graph_in_Pocket(rg_dList< GEdge >* edgeList, rg_dList< GVertex >* vertexList);

                    
	Graph_in_Pocket(const Graph_in_Pocket& graph);
    ~Graph_in_Pocket();

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Get Functions
	rg_dList< GEdge >*     getEdgeList();
    rg_dList< GVertex >*   getVertexList();
	
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Set Functions
	void                   setEdgeList(rg_dList< GEdge >* edgeList);
    void                   setVertexList(rg_dList< GVertex >* vertexList);

    GEdge*                 addEdge( const GEdge& edge );
    GVertex*               addVertex( const GVertex& vertex);
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Operator Overloadings
    //Graph& operator =(const Graph& graph);



    void setGraph(rg_REAL** adjacencyMatrix, const rg_INT& n);
    
    void setUserDataForVerticesInGraph(void** vertexData, const rg_INT& n);    //n vertex data
    void setUserDataForEdgesInGraph(void*** edgeData, const rg_INT& n);   //n by n matrix



    /////////////////////////////////////////////////////////////////////////////////
	//
	// Dijkstra    
    void runDijkstra(GVertex* source, 
                        GVertex* destination, 
                        rg_dList< GEdge* >* shortestPath);
    void runDijkstra(void* sourceVOID, 
                        void* destinationVOID, 
                        rg_dList< GEdge* >* shortestPath);
    void runDijkstra(const rg_INT& sourceIndex, 
                        const rg_INT& destinationIndex, 
                        rg_dList< GEdge* >* shortestPath);

    //void runDijkstra(GVertex* source);
    


	void convertEdgePathToVertexPath(     rg_dList< GEdge* >* edgeList, 
                                        rg_dList< GVertex* >* vertexList);
    void convertEdgePathToVertexDataPath( rg_dList< GEdge* >* edgeList, 
                                        rg_dList<void*>* vertexDataList);
	void convertEdgePathToEdgeDataPath(  rg_dList< GEdge* >* edgeList, 
                                        rg_dList<void*>* edgeDataList);
	void convertEdgePathToIndexPath(    rg_dList< GEdge* >* edgeList, 
                                        rg_dList<rg_INT>* indexList);

private:

    void settingSourceAndDestination(void * sourceVOID, void* destinationVOID, 
                                        GVertex*& source, GVertex*& destination);    
    void settingSourceAndDestination(const rg_INT& sourceIndex, const rg_INT& destinationIndex, 
                                        GVertex*& source, GVertex*& destination);
    

    void executeDijkstraAlgorithm(GVertex* source, GVertex* destination, 
                                    rg_dList< GEdge* >* shortestPath);
        rg_INT findMaxIndexOfVertex();
    
        void initialize(GVertex* source, 
                        PriorityQueueForGVertex* unvisitedVertices, 
                        PQNode** vertexIndexArray);

        void calculateCostFromSourceToDestination(GVertex* source, GVertex* destination, 
                                                    PriorityQueueForGVertex* unvisitedVertices, 
                                                    PQNode** vertexIndexArray);
        void calculateCostFromSourceToAllOtherVertices(GVertex* source, 
                                                    PriorityQueueForGVertex* unvisitedVertices, 
                                                    PQNode** vertexIndexArray);
            void updateCostOfAdjacentVertex(GVertex* nearestVertex,
                                             GEdge* adjacentEdge,
                                             GVertex* incidentVertex,
                                             rg_dNode< GVertex* >* incidentVertexIterator,
                                             rg_dList< GVertex* >* unvisitedVertexList, 
                                             PQNode** vertexIndexArray);
            
        void defineShortestPath(GVertex* source, GVertex* destination,
                               rg_dList< GEdge* >* shortestPath);
};

#endif
