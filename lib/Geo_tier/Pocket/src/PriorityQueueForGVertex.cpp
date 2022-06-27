#include "PriorityQueueForGVertex.h"


PriorityQueueForGVertex::PriorityQueueForGVertex()
{
	m_PQNodes		= new PQNode[DEFAULT_MAX_NUM_OF_NODES+1];	//0번째는 비워 둠.
	m_numOfNodes	= 0;
	m_maxNumOfNodes = DEFAULT_MAX_NUM_OF_NODES;
}



PriorityQueueForGVertex::PriorityQueueForGVertex(const PriorityQueueForGVertex& pq)
{
	m_PQNodes		= pq.m_PQNodes;
	m_numOfNodes	= pq.m_numOfNodes;
	m_maxNumOfNodes = pq.m_maxNumOfNodes;
}

PriorityQueueForGVertex::PriorityQueueForGVertex( const rg_INT& maxNumOfNodes )
{
	m_PQNodes		= new PQNode[maxNumOfNodes+1];			//0번째는 비워 둠.
	m_numOfNodes	= 0;
	m_maxNumOfNodes = maxNumOfNodes;
}

PriorityQueueForGVertex::~PriorityQueueForGVertex()
{
	if ( m_PQNodes != rg_NULL )
		delete [] m_PQNodes;
}


/////////////////////////////////////////////////////////////////////////////////
//
// Get Functions
PQNode* PriorityQueueForGVertex::getPQNodes() const
{
	return m_PQNodes;
}

	
rg_INT  PriorityQueueForGVertex::getNumOfNodes() const
{
	return m_numOfNodes;
}


rg_INT  PriorityQueueForGVertex::getMaxNumOfNodes() const
{
	return m_maxNumOfNodes;
}


/////////////////////////////////////////////////////////////////////////////////
//
// Set Functions
void PriorityQueueForGVertex::setPQNodes(PQNode* pqNodes)
{
	m_PQNodes = pqNodes;
}


void PriorityQueueForGVertex::setNumOfNodes(const rg_INT& numOfNodes)
{
	m_numOfNodes = numOfNodes;
}


void PriorityQueueForGVertex::setMaxNumOfNodes(const rg_INT& maxNumOfNodes)
{
	m_maxNumOfNodes = maxNumOfNodes;
}


/////////////////////////////////////////////////////////////////////////////////
//
// Operator Overloadings
PriorityQueueForGVertex& PriorityQueueForGVertex::operator =(const PriorityQueueForGVertex& pq)
{
	if ( this == &pq)
		return *this;

	m_PQNodes		= pq.m_PQNodes;
	m_numOfNodes	= pq.m_numOfNodes;
	m_maxNumOfNodes = pq.m_maxNumOfNodes;

	return *this;
}




//For Dijkstra's Algorithm
PQNode* PriorityQueueForGVertex::insertPQNodeWithoutBubbling(GVertex* gVertex)
{
	m_PQNodes[m_numOfNodes + 2].setGVertex( gVertex );
	m_PQNodes[m_numOfNodes + 2].setKeyValue( rg_MAX_REAL );

	m_numOfNodes++;

	return &(m_PQNodes[m_numOfNodes + 1]);
}


//For Dijkstra's Algorithm
PQNode* PriorityQueueForGVertex::insertSourcePQNode(GVertex* gVertex)
{
	m_PQNodes[1].setGVertex( gVertex );
	m_PQNodes[1].setKeyValue( 0 );
	
	return &(m_PQNodes[1]);	
}


PQNode* PriorityQueueForGVertex::getMinPQNode()
{
	return &(m_PQNodes[1]);
}

void PriorityQueueForGVertex::killMinPQNode()
{
	m_PQNodes[1] = m_PQNodes[m_numOfNodes];
	m_PQNodes[m_numOfNodes].clearMyself();
	m_numOfNodes--;

	downHeapBubbling( &(m_PQNodes[1]) );
}

void PriorityQueueForGVertex::killMinPQNode( PQNode** pPQNodeOfGVertex )
{
	pPQNodeOfGVertex[ m_PQNodes[1].getGVertex()->getID() ] = pPQNodeOfGVertex[ m_PQNodes[m_numOfNodes].getGVertex()->getID() ];
	pPQNodeOfGVertex[ m_PQNodes[m_numOfNodes].getGVertex()->getID() ] = rg_NULL;

	m_PQNodes[1] = m_PQNodes[m_numOfNodes];
	m_PQNodes[m_numOfNodes].clearMyself();
	m_numOfNodes--;

	downHeapBubbling( &(m_PQNodes[1]), pPQNodeOfGVertex );	
}

void PriorityQueueForGVertex::modifyKeyOfPQNode( const rg_REAL& key, PQNode* pqNode )
{
	pqNode->setKeyValue( key );

	upHeapBubbling(   pqNode );
	downHeapBubbling( pqNode );
}

void PriorityQueueForGVertex::modifyKeyOfPQNode( const rg_REAL& key, PQNode* pqNode, PQNode** pPQNodeOfGVertex )
{
	pqNode->setKeyValue( key );

	upHeapBubbling(   pqNode, pPQNodeOfGVertex );
	downHeapBubbling( pqNode, pPQNodeOfGVertex );	
}


void PriorityQueueForGVertex::upHeapBubbling( PQNode* pqNode )
{
	rg_INT indexOfNode = (int)(pqNode-m_PQNodes);

	while ( indexOfNode > 1)
	{
		if ( m_PQNodes[indexOfNode].getKeyValue() < m_PQNodes[indexOfNode/2].getKeyValue()  )
		{
			PQNode tempNode          = m_PQNodes[indexOfNode];
			m_PQNodes[indexOfNode]   = m_PQNodes[indexOfNode/2];
			m_PQNodes[indexOfNode/2] = tempNode;
		}
		else
		{
			break;
		}

		indexOfNode = indexOfNode / 2;
	}	
}

void PriorityQueueForGVertex::upHeapBubbling( PQNode* pqNode, PQNode** pPQNodeOfGVertex )
{
	rg_INT indexOfNode = (int)(pqNode-m_PQNodes);

	while ( indexOfNode > 1)
	{
		if ( m_PQNodes[indexOfNode].getKeyValue() < m_PQNodes[indexOfNode/2].getKeyValue()  )
		{
			PQNode* pTempNode = pPQNodeOfGVertex[ m_PQNodes[indexOfNode].getGVertex()->getID() ];
			pPQNodeOfGVertex[ m_PQNodes[indexOfNode].getGVertex()->getID() ]
				= pPQNodeOfGVertex[ m_PQNodes[indexOfNode/2].getGVertex()->getID() ];
			pPQNodeOfGVertex[ m_PQNodes[indexOfNode/2].getGVertex()->getID() ] = pTempNode;

			PQNode tempNode          = m_PQNodes[indexOfNode];
			m_PQNodes[indexOfNode]   = m_PQNodes[indexOfNode/2];
			m_PQNodes[indexOfNode/2] = tempNode;
		}
		else
		{
			break;
		}

		indexOfNode = indexOfNode / 2;
	}		
}
void PriorityQueueForGVertex::downHeapBubbling( PQNode* pqNode )
{
	rg_INT indexOfNode = (int)(pqNode-m_PQNodes);

	while ( indexOfNode*2 <= m_numOfNodes )
	{
		if ( m_PQNodes[indexOfNode*2 + 1].getGVertex() == rg_NULL ) //마지막 node의 sibling이 없을 경
		{
			indexOfNode = indexOfNode*2;
		}
		else 
		{
			if ( m_PQNodes[indexOfNode*2].getKeyValue() < m_PQNodes[indexOfNode*2 + 1].getKeyValue() )
			{
				indexOfNode = indexOfNode*2;
			}
			else
			{
				indexOfNode = indexOfNode*2 + 1;				
			}
		}

		PQNode* childNode = &m_PQNodes[indexOfNode];
		if ( m_PQNodes[indexOfNode/2].getKeyValue() > childNode->getKeyValue()  )
		{
			PQNode tempNode          = m_PQNodes[indexOfNode/2];
			m_PQNodes[indexOfNode/2]   = *childNode;
			*childNode = tempNode;
		}
		else
		{
			break;
		}		
	}	
}

void PriorityQueueForGVertex::downHeapBubbling( PQNode* pqNode, PQNode** pPQNodeOfGVertex )
{
	rg_INT indexOfNode = (int)(pqNode-m_PQNodes);

	while ( indexOfNode*2 <= m_numOfNodes )
	{
		if ( indexOfNode*2 == m_numOfNodes ) //마지막 node의 sibling이 없을 경
		{
			indexOfNode = indexOfNode*2;
		}
		else 
		{
			if ( m_PQNodes[indexOfNode*2].getKeyValue() < m_PQNodes[indexOfNode*2 + 1].getKeyValue() )
			{
				indexOfNode = indexOfNode*2;
			}
			else
			{
				indexOfNode = indexOfNode*2 + 1;				
			}
		}

		PQNode* childNode = &m_PQNodes[indexOfNode];
		if ( m_PQNodes[indexOfNode/2].getKeyValue() > childNode->getKeyValue()  )
		{
			PQNode* pTempNode = pPQNodeOfGVertex[ m_PQNodes[indexOfNode/2].getGVertex()->getID() ];
			pPQNodeOfGVertex[ m_PQNodes[indexOfNode/2].getGVertex()->getID() ]
				= pPQNodeOfGVertex[ childNode->getGVertex()->getID() ];
			pPQNodeOfGVertex[ childNode->getGVertex()->getID() ] = pTempNode;

			PQNode tempNode          = m_PQNodes[indexOfNode/2];
			m_PQNodes[indexOfNode/2]   = *childNode;
			*childNode = tempNode;
		}
		else
		{
			break;
		}		
	}		
}