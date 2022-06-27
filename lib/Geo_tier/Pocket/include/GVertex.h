#ifndef _GVertex_H_
#define _GVertex_H_

#include "rg_Const.h"
#include "rg_dList.h"


//#include <rg_dList>
//using namespace std;


class GEdge;

class GVertex
{
private:
	rg_INT		        	    m_ID;
	rg_dList< GEdge* >	        m_adjacentEdges;
	
	void*		        	    m_vertexData; 

    GEdge*     		            m_prevEdge;


public:
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Constructor
	GVertex();
    GVertex(const rg_INT& id);
    GVertex(const rg_INT& id, void* vertexData);
	//GVertex(const rg_INT id, 
    //                rg_dList< GEdge* >* adjacentEdges, 
    //                const rg_REAL& accumulatedCostFromSource,
    //                GEdge* prevEdge,
    //                void* vertexData);
	GVertex(const GVertex& GVertex);
	~GVertex();

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Get Functions
	rg_INT				        getID() const;
	rg_dList< GEdge* >*	        getAdjacentEdges();
    GEdge*                      getPreEdge();
	void*					    getVertexData();

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Set Functions
	void 				        setID(const rg_INT& id);
	void                	    setAdjacentEdges(rg_dList< GEdge *> adjacentEdges);
    void                        setPreEdge(GEdge* prevEdge);
	void					    setVertexData(void* vertexData);

    GEdge**                     addAdjacentEdge(GEdge* arc);
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Operator Overloadings
	GVertex& operator =(const GVertex& GVertex);

};

#endif
