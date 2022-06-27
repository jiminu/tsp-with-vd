#include "PQNode.h"


PQNode::PQNode()
{
	m_keyValue = rg_MAX_REAL;
	m_GVertex  = rg_NULL;
}

PQNode::PQNode( const PQNode& pqNode )
{
	m_keyValue = pqNode.m_keyValue;
	m_GVertex  = pqNode.m_GVertex;
}

rg_REAL PQNode::getKeyValue() const
{
	return m_keyValue;	
}

GVertex* PQNode::getGVertex()
{
	return m_GVertex;	
}

void PQNode::setKeyValue( const rg_REAL& keyValue )
{
	m_keyValue = keyValue;	
}

void PQNode::setGVertex( GVertex* gVertex )
{
	m_GVertex = gVertex;	
}

PQNode& PQNode::operator=( const PQNode& pqNode )
{
	if ( this == &pqNode )
        return *this;

    m_keyValue = pqNode.m_keyValue;
	m_GVertex  = pqNode.m_GVertex;

    return *this;	
}

PQNode::~PQNode()
{
}



void PQNode::clearMyself()
{
	m_keyValue = rg_MAX_REAL;
	m_GVertex  = rg_NULL;
}
