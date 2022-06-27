#ifndef _PQNODE_H
#define _PQNODE_H

#include "rg_Const.h"
#include "GVertex.h"

class PQNode
{
private:
    rg_REAL		m_keyValue; //accum distance from source
	GVertex*	m_GVertex;
	

public:
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Constructor
	PQNode();
    PQNode(const PQNode& pqNode);
	~PQNode();

	void clearMyself();

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Get Functions
	rg_REAL		        getKeyValue() const;	
    GVertex*	        getGVertex();
    
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Set Functions
	void			    setKeyValue(const rg_REAL& keyValue);
	void                setGVertex(GVertex* gVertex);

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Operator Overloadings
	PQNode& operator =(const PQNode& pqNode);

};

#endif
