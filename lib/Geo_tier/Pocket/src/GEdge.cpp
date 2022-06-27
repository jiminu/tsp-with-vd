#include "GEdge.h"

GEdge::GEdge()
: m_ID(-1), m_cost(0.0)
{
    m_startVertex = rg_NULL;
    m_endVertex   = rg_NULL;
    m_edgeData   = rg_NULL;
}



GEdge::GEdge(const rg_INT& id)
: m_ID(id), m_cost(0.0)
{
    m_startVertex = rg_NULL;
    m_endVertex   = rg_NULL;
    m_edgeData   = rg_NULL;
}




GEdge::GEdge(const rg_INT& id, void* edgeData)
: m_ID(id), m_cost(0.0)
{
    m_edgeData   = edgeData;

    m_startVertex = rg_NULL;
    m_endVertex   = rg_NULL;    
}



GEdge::GEdge(const rg_INT& id,
                GVertex* startVertex,
                GVertex* endVertex,
                const rg_REAL& cost,
                void* edgeData)
{
    m_ID        = id;
    m_startVertex = startVertex;
    m_endVertex   = endVertex;
    m_cost      = cost;
    m_edgeData   = edgeData;
}


GEdge::GEdge(const GEdge& GEdge)
{
    m_ID        = GEdge.m_ID;
    m_startVertex = GEdge.m_startVertex;
    m_endVertex   = GEdge.m_endVertex;
    m_cost      = GEdge.m_cost;
    m_edgeData   = GEdge.m_edgeData;
}


GEdge::~GEdge()
{
}


/////////////////////////////////////////////////////////////////////////////////
//
// Get Functions
rg_INT GEdge::getID() const
{
    return m_ID;
}



GVertex* GEdge::getStartVertex()
{
    return m_startVertex;
}



GVertex* GEdge::getEndVertex()
{
    return m_endVertex;
}



rg_REAL GEdge::getCost() const
{
    return m_cost;
}



void* GEdge::getEdgeData()
{
    return m_edgeData;
}

/////////////////////////////////////////////////////////////////////////////////
//
// Set Functions
void GEdge:: setID(const rg_INT& id)
{
    m_ID = id;
}



void GEdge::setStartVertex(GVertex* startVertex)
{
    m_startVertex = startVertex;
}




void GEdge::setEndVertex(GVertex* endVertex)
{
    m_endVertex   = endVertex;
}




void GEdge::setCost(const rg_REAL& cost)
{
    m_cost      = cost;
}





void GEdge::setEdgeData(void* edgeData)
{
    m_edgeData   = edgeData;
}


/////////////////////////////////////////////////////////////////////////////////
//
// Operator Overloadings
GEdge& GEdge::operator =(const GEdge& GEdge)
{
    if ( this == &GEdge)
        return *this;

    m_ID        = GEdge.m_ID;
    m_startVertex = GEdge.m_startVertex;
    m_endVertex   = GEdge.m_endVertex;
    m_cost      = GEdge.m_cost;
    m_edgeData   = GEdge.m_edgeData;

    return *this;
}
