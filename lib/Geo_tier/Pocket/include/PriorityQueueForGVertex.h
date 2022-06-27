#ifndef PRIORITY_QUEUE_FOR_GVERTEX_H_
#define PRIORITY_QUEUE_FOR_GVERTEX_H_

#include "PQNode.h"

const rg_INT DEFAULT_MAX_NUM_OF_NODES = 100;

class PriorityQueueForGVertex
{
private:
    PQNode*		m_PQNodes;
	rg_INT		m_numOfNodes;

	rg_INT		m_maxNumOfNodes;	

public:
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Constructor
	PriorityQueueForGVertex();
	PriorityQueueForGVertex(const rg_INT& maxNumOfNodes);
    PriorityQueueForGVertex(const PriorityQueueForGVertex& pq);
	~PriorityQueueForGVertex();

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Get Functions
	PQNode*		        getPQNodes() const;	
    rg_INT				getNumOfNodes() const;
	rg_INT				getMaxNumOfNodes() const;
    
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Set Functions
	void			    setPQNodes(PQNode* pqNodes);
	void                setNumOfNodes(const rg_INT& numOfNodes);
	void				setMaxNumOfNodes(const rg_INT& maxNumOfNodes);

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Operator Overloadings
	PriorityQueueForGVertex& operator =(const PriorityQueueForGVertex& pq);



	PQNode* insertPQNodeWithoutBubbling(GVertex* gVertex);

	PQNode* insertSourcePQNode(GVertex* gVertex);

	PQNode* getMinPQNode();
	void    killMinPQNode();
	void	killMinPQNode(PQNode** pPQNodeOfGVertex);

	void modifyKeyOfPQNode( const rg_REAL& key, PQNode* pqNode );
	void modifyKeyOfPQNode( const rg_REAL& key, PQNode* pqNode, PQNode** pPQNodeOfGVertex );
		void upHeapBubbling( PQNode* pqNode );
		void upHeapBubbling( PQNode* pqNode, PQNode** pPQNodeOfGVertex );
		void downHeapBubbling( PQNode* pqNode );
		void downHeapBubbling( PQNode* pqNode, PQNode** pPQNodeOfGVertex );


};

#endif