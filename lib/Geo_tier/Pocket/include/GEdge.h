#ifndef _GEdge_H_
#define _GEdge_H_

#include "rg_Const.h"

class GVertex;

class GEdge
{
private:
    rg_INT          m_ID;

	GVertex*	    m_startVertex;
	GVertex*	    m_endVertex;

	rg_REAL     	m_cost;

	void*   		m_edgeData;


public:
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Constructor
	GEdge();
    GEdge(const rg_INT& id);
    GEdge(const rg_INT& id, void* edgeData);
	GEdge(const rg_INT& id,
                    GVertex* startVertex,
                    GVertex* endVertex,
                    const rg_REAL& cost,
                    void* edgeData);
    GEdge(const GEdge& GEdge);
	~GEdge();

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Get Functions
	rg_INT				getID() const;
	GVertex*	        getStartVertex();
    GVertex*	        getEndVertex();
    rg_REAL		        getCost() const;
	void*				getEdgeData();

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Set Functions
	void				    setID(const rg_INT& id);
	void	                setStartVertex(GVertex* startVertex);
    void    	            setEndVertex(GVertex* endVertex);
    void		            setCost(const rg_REAL& cost);
	void					setEdgeData(void* edgeData);

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Operator Overloadings
	GEdge& operator =(const GEdge& GEdge);

};

#endif
