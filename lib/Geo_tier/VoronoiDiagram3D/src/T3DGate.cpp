#include "T3DGate.h"
using namespace V::GeometryTier;


T3DGate::T3DGate()
: m_tetrahedron(rg_NULL), m_mateVertexPos(-1), m_indexOfPriorEdge(rg_UNKNOWN)
{
    for (rg_INT i=0; i<3; i++)  {
        m_vertex[i]                   = rg_NULL;
        m_vertexForNextTetrahedron[i] = rg_NULL;
        m_isValidVertex[i]            = rg_FALSE;
    }
}



T3DGate::T3DGate(T3DTetrahedron* tetrahedron, const rg_INT& mateVertexPos)
: m_tetrahedron(tetrahedron), m_mateVertexPos(mateVertexPos), m_indexOfPriorEdge(rg_UNKNOWN)
{
    m_tetrahedron->getVerticesOnFace(m_mateVertexPos, m_vertex);
    for (rg_INT i=0; i<3; i++)  {
        m_vertexForNextTetrahedron[i] = rg_NULL;
        m_isValidVertex[i]            = rg_FALSE;
    }
}



T3DGate::T3DGate(const T3DGate& gate)
: m_tetrahedron(gate.m_tetrahedron), m_mateVertexPos(gate.m_mateVertexPos),
  m_indexOfPriorEdge(gate.m_indexOfPriorEdge)
{
    for (rg_INT i=0; i<3; i++)  {
        m_vertex[i]                   = gate.m_vertex[i];
        m_vertexForNextTetrahedron[i] = gate.m_vertexForNextTetrahedron[i];
        m_isValidVertex[i]            = gate.m_isValidVertex[i];
    }
}



T3DGate::~T3DGate()
{
}




T3DTetrahedron* T3DGate::getTetrahedron() const
{
    return m_tetrahedron;
}



rg_INT          T3DGate::getMateVertexPos() const
{
    return m_mateVertexPos;
}



T3DVertex*      T3DGate::getVertex(const rg_INT& i) const
{
    return m_vertex[i];
}



T3DVertex**     T3DGate::getVertices() 
{
    return m_vertex;
}



T3DVertex*      T3DGate::getVertexForNextTetrahedron(const rg_INT& i) const
{
    return m_vertexForNextTetrahedron[i];
}



T3DVertex**     T3DGate::getVerticesForNextTetrahedron() 
{
    return m_vertexForNextTetrahedron;
}



T3DTetrahedron* T3DGate::getIncidentNeighbor() const
{
    return m_tetrahedron->getNeighbor(m_mateVertexPos);
}



rg_INT T3DGate::getIndexOfPriorEdge() const
{
    return m_indexOfPriorEdge;
}



void T3DGate::getVerticesOfPriorEdge(T3DVertex*& startVertex, T3DVertex*& endVertex)
{
    switch ( m_indexOfPriorEdge )  {
        case 0:
            startVertex = m_vertex[0];
            endVertex   = m_vertex[1];
            break;
        case 1:
            startVertex = m_vertex[1];
            endVertex   = m_vertex[2];
            break;
        case 2:
            startVertex = m_vertex[2];
            endVertex   = m_vertex[0];
            break;
        default:
            startVertex = rg_NULL;
            endVertex   = rg_NULL;
            break;
    }
}



rg_FLAG         T3DGate::isValidVertexForNextTetrahedron(const rg_INT& i) const
{
    return m_isValidVertex[i];
}



rg_FLAG T3DGate::compareGateByVertices( const T3DGate& gate )
{
    T3DVertex* comparedVertex[5];
    comparedVertex[0] = gate.m_vertex[2];
    comparedVertex[1] = gate.m_vertex[1];
    comparedVertex[2] = gate.m_vertex[0];
    comparedVertex[3] = gate.m_vertex[2];
    comparedVertex[4] = gate.m_vertex[1];

    for (rg_INT i=0; i<3; i++)  {
        if (    m_vertex[0] == comparedVertex[i] 
             && m_vertex[1] == comparedVertex[i+1] 
             && m_vertex[2] == comparedVertex[i+2])  {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_FLAG T3DGate::compareGateByVertices( T3DVertex** vertexOnGate )
{
    T3DVertex* comparedVertex[5];
    comparedVertex[0] = vertexOnGate[2];
    comparedVertex[1] = vertexOnGate[1];
    comparedVertex[2] = vertexOnGate[0];
    comparedVertex[3] = vertexOnGate[2];
    comparedVertex[4] = vertexOnGate[1];

    for (rg_INT i=0; i<3; i++)  {
        if (    m_vertex[0] == comparedVertex[i] 
             && m_vertex[1] == comparedVertex[i+1] 
             && m_vertex[2] == comparedVertex[i+2])  {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



void T3DGate::setTetrahedron(T3DTetrahedron* tetrahedron)
{
    m_tetrahedron = tetrahedron;
    m_tetrahedron->getVerticesOnFace(m_mateVertexPos, m_vertex);
}



void T3DGate::setMateVertexPos(const rg_INT& pos)
{
    m_mateVertexPos = pos;
    m_tetrahedron->getVerticesOnFace(m_mateVertexPos, m_vertex);
}



void T3DGate::setVertexForNextTetrahedron(const rg_INT& i, T3DVertex* vertex)
{
    m_isValidVertex[i]            = rg_TRUE;
    m_vertexForNextTetrahedron[i] = vertex;
}



void T3DGate::setGate(T3DTetrahedron* tetrahedron, const rg_INT& mateVertexPos)
{
    m_tetrahedron   = tetrahedron;
    m_mateVertexPos = mateVertexPos;

    m_tetrahedron->getVerticesOnFace(m_mateVertexPos, m_vertex);
    for (rg_INT i=0; i<3; i++)  {
        m_vertexForNextTetrahedron[i] = rg_NULL;
        m_isValidVertex[i]            = rg_FALSE;
    }
}



void T3DGate::setIndexOfPriorEdge(const rg_INT& i)
{
    m_indexOfPriorEdge = i;
}



T3DGate& T3DGate::operator =(const T3DGate& gate)
{
    if ( this == &gate )
        return *this;

    m_tetrahedron   = gate.m_tetrahedron;
    m_mateVertexPos = gate.m_mateVertexPos;
    m_indexOfPriorEdge = gate.m_indexOfPriorEdge;

    for (rg_INT i=0; i<3; i++)  {
        m_vertex[i]                   = gate.m_vertex[i];
        m_vertexForNextTetrahedron[i] = gate.m_vertexForNextTetrahedron[i];
        m_isValidVertex[i]            = gate.m_isValidVertex[i];
    }
    
    return *this;
}



