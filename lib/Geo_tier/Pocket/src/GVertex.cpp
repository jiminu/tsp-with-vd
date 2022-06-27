#include "GVertex.h"
#include <float.h>

GVertex::GVertex()
: m_ID(-1)
{
    m_vertexData  = rg_NULL;
    m_prevEdge    = rg_NULL;
}




GVertex::GVertex(const rg_INT& id)
: m_ID(id)
{
    m_vertexData  = rg_NULL;
    m_prevEdge    = rg_NULL;
}





GVertex::GVertex(const rg_INT& id, void* vertexData)
: m_ID(id)
{
    m_vertexData = vertexData;
    m_prevEdge   = rg_NULL;    
}





//GVertex::GVertex(const rg_INT id, 
//                    rg_dList< GEdge* >* adjacentEdges, 
//                    const rg_REAL& accumulatedCostFromSource,
//                    GEdge* prevEdge,
//                    void* vertexData)
//{
//    m_ID                        = id;
//    m_adjacentEdges              = *adjacentEdges;
//    m_accumulatedCostFromSource = accumulatedCostFromSource;
//    m_prevEdge                   = prevEdge;
//    m_vertexData                  = vertexData;
//}

	


GVertex::GVertex(const GVertex& GVertex)
{
    m_ID                        = GVertex.m_ID;
    m_adjacentEdges             = GVertex.m_adjacentEdges;
    m_prevEdge                  = GVertex.m_prevEdge;
    m_vertexData                = GVertex.m_vertexData;
}



GVertex::~GVertex()
{
}



/////////////////////////////////////////////////////////////////////////////////
//
// Get Functions	
rg_INT GVertex::getID() const
{
    return m_ID;
}


rg_dList< GEdge* >*	GVertex::getAdjacentEdges()
{
    return &m_adjacentEdges;
}



GEdge* GVertex::getPreEdge()
{
    return m_prevEdge;
}



void* GVertex::getVertexData()
{
    return m_vertexData;
}






/////////////////////////////////////////////////////////////////////////////////
//
// Set Functions
void GVertex::setID(const rg_INT& id)
{
    m_ID = id;
}

void GVertex::setAdjacentEdges(rg_dList< GEdge *> adjacentEdges)
{
    m_adjacentEdges = adjacentEdges;
}




void GVertex::setPreEdge(GEdge* prevEdge)
{
    m_prevEdge = prevEdge;
}

void GVertex::setVertexData(void* vertexData)
{
    m_vertexData = vertexData;
}



GEdge** GVertex::addAdjacentEdge(GEdge* arc)
{
    return m_adjacentEdges.add( arc );
}


/////////////////////////////////////////////////////////////////////////////////
//
// Operator Overloadings
GVertex& GVertex::operator =(const GVertex& GVertex)
{
    if ( this == &GVertex )
        return *this;

    m_ID                         = GVertex.m_ID;
    m_adjacentEdges              = GVertex.m_adjacentEdges;
    m_prevEdge                   = GVertex.m_prevEdge;
    m_vertexData                 = GVertex.m_vertexData;

    return *this;
}
