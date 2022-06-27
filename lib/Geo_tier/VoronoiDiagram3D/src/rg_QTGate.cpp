#include "rg_QTGate.h"

#include "rg_QTVertex.h"
#include "rg_QTTetrahedron.h"
using namespace V::GeometryTier;

QTGate::QTGate()
{
    m_vertex[0] = rg_NULL;
    m_vertex[1] = rg_NULL;

    m_bigTetrahedron = rg_NULL;
}



QTGate::QTGate(QTVertex* vertex1, QTVertex* vertex2)
{
    m_vertex[0] = vertex1;
    m_vertex[1] = vertex2;

    m_bigTetrahedron = rg_NULL;
}



QTGate::QTGate(const QTGate& gate)
{
    m_vertex[0] = gate.m_vertex[0];
    m_vertex[1] = gate.m_vertex[1];

    m_bigTetrahedron = gate.m_bigTetrahedron;
    m_smallTetrahedra.duplicateList(gate.m_smallTetrahedra);
}



QTGate::~QTGate()
{
}




QTVertex* QTGate::getVertex(const rg_INDEX& i)
{
    if ( i<0 || i>=2 )
        return rg_NULL;
    else
        return m_vertex[i];
}



QTVertex** QTGate::getVertices()
{
    return m_vertex;
}



QTTetrahedron* QTGate::getBigTetrahedron()
{
    return m_bigTetrahedron;
}



rg_dList<QTTetrahedron*>* QTGate::getSmallTetrahedra()
{
    return &m_smallTetrahedra;
}




void QTGate::setVertex(const rg_INDEX& i, QTVertex* vertex)
{
    if ( i<0 || i>=2 )
        return;
    else
        m_vertex[i] = vertex;
}



void QTGate::setVertices(QTVertex* vertex1, QTVertex* vertex2)
{
    m_vertex[0] = vertex1;
    m_vertex[1] = vertex2;
}



void QTGate::setVertices(QTVertex** vertices)
{
    m_vertex[0] = vertices[0];
    m_vertex[1] = vertices[1];
}



void QTGate::setBigTetrahedron(QTTetrahedron* bigTetrahedron)
{
    m_bigTetrahedron = bigTetrahedron;
}



void QTGate::addSmallTetrahedron(QTTetrahedron* smallTetrahedron)
{
    m_smallTetrahedra.add( smallTetrahedron );
}




QTGate& QTGate::operator =(const QTGate& gate)
{
    if ( this == &gate )
        return *this;

    m_vertex[0] = gate.m_vertex[0];
    m_vertex[1] = gate.m_vertex[1];

    m_bigTetrahedron = gate.m_bigTetrahedron;
    m_smallTetrahedra.duplicateList(gate.m_smallTetrahedra);

    return *this;
}



rg_FLAG QTGate::operator==(const QTGate& gate)
{
    if (    (m_vertex[0]==gate.m_vertex[0] && m_vertex[1]==gate.m_vertex[1])
         || (m_vertex[0]==gate.m_vertex[1] && m_vertex[1]==gate.m_vertex[0])  )
        return rg_TRUE;
    else
        return rg_FALSE;
}

