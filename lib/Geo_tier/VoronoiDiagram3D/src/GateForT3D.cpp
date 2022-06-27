#include "GateForT3D.h"
using namespace V::GeometryTier;

GateForT3D::GateForT3D()
: m_vertex(rg_NULL), m_mateVertex(rg_NULL), m_tetrahedron(rg_NULL)
{
}



GateForT3D::GateForT3D(T3DVertex* vertex, T3DVertex* mateVertex, T3DTetrahedron* tetrahedron)
: m_vertex(vertex), m_mateVertex(mateVertex), m_tetrahedron(tetrahedron)
{
}



GateForT3D::GateForT3D(const GateForT3D& gate)
: m_vertex(gate.m_vertex), m_mateVertex(gate.m_mateVertex), m_tetrahedron(gate.m_tetrahedron)
{
}



GateForT3D::~GateForT3D()
{
}




T3DVertex* GateForT3D::getVertex() const
{
    return m_vertex;
}



T3DVertex* GateForT3D::getMateVertex() const
{
    return m_mateVertex;
}



T3DTetrahedron* GateForT3D::getTetrahedron() const
{
    return m_tetrahedron;
}




void GateForT3D::setVertex(T3DVertex* vertex)
{
    m_vertex = vertex;
}



void GateForT3D::setMateVertex(T3DVertex* mateVertex)
{
    m_mateVertex = mateVertex;
}



void GateForT3D::setTetrahedron(T3DTetrahedron* tetrahedron)
{
    m_tetrahedron = tetrahedron;
}




GateForT3D& GateForT3D::operator =(const GateForT3D& gate)
{
    if ( this == &gate )
        return *this;

    m_vertex      = gate.m_vertex;
    m_mateVertex  = gate.m_mateVertex;
    m_tetrahedron = gate.m_tetrahedron;

    return *this;
}



