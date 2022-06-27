#include "T3DFace.h"
#include "T3DTetrahedron.h"
using namespace V::GeometryTier;



T3DFace::T3DFace()
: m_tetrahedron(rg_NULL), m_mateVertexPos(-1)
{
}



T3DFace::T3DFace(T3DTetrahedron* tetrahedron, const rg_INT& mateVertexPos)
: m_tetrahedron(tetrahedron), m_mateVertexPos(mateVertexPos)
{
}



T3DFace::T3DFace(const T3DFace& face)
{
    m_tetrahedron   = face.m_tetrahedron;
    m_mateVertexPos = face.m_mateVertexPos;
}



T3DFace::~T3DFace()
{
}




T3DTetrahedron* T3DFace::getTetrahedron() const
{
    return m_tetrahedron;
}



rg_INT T3DFace::getMateVertexPos() const
{
    return m_mateVertexPos;
}



T3DTetrahedron* T3DFace::getIncidentNeighbor() const
{
    return m_tetrahedron->getNeighbor(m_mateVertexPos);
}



void T3DFace::getBoundingVertices(T3DVertex** vertices) const
{
    switch( m_mateVertexPos )  {
        case 0:
            vertices[0] = m_tetrahedron->getVertex(1);
            vertices[1] = m_tetrahedron->getVertex(3);
            vertices[2] = m_tetrahedron->getVertex(2);
            break;
        case 1:
            vertices[0] = m_tetrahedron->getVertex(0);
            vertices[1] = m_tetrahedron->getVertex(2);
            vertices[2] = m_tetrahedron->getVertex(3);
            break;
        case 2:
            vertices[0] = m_tetrahedron->getVertex(0);
            vertices[1] = m_tetrahedron->getVertex(3);
            vertices[2] = m_tetrahedron->getVertex(1);
            break;
        case 3:
            vertices[0] = m_tetrahedron->getVertex(0);
            vertices[1] = m_tetrahedron->getVertex(1);
            vertices[2] = m_tetrahedron->getVertex(2);
            break;
        default:
            break;
    }
}




void T3DFace::setTetrahedron(T3DTetrahedron* tetrahedron)
{
    m_tetrahedron = tetrahedron;
}



void T3DFace::setMateVertexPos(const rg_INT& mateVertexPos)
{
    m_mateVertexPos = mateVertexPos;
}



void T3DFace::setFace(T3DTetrahedron* tetrahedron, const rg_INT& mateVertexPos)
{
    m_tetrahedron   = tetrahedron;
    m_mateVertexPos = mateVertexPos;
}



rg_FLAG T3DFace::compareFaceByVertices( const T3DFace& face )
{
    T3DVertex* thisVertices[3];
    this->getBoundingVertices(thisVertices);

    T3DVertex* tempVertices[3];
    face.getBoundingVertices(tempVertices);
    T3DVertex* targetVertices[3];
    targetVertices[0] = tempVertices[2];
    targetVertices[1] = tempVertices[1];
    targetVertices[2] = tempVertices[0];
    targetVertices[3] = tempVertices[2];
    targetVertices[4] = tempVertices[1];

    for (rg_INT i=0; i<3; i++)  {
        if (    thisVertices[0] == targetVertices[i] 
             && thisVertices[1] == targetVertices[i+1] 
             && thisVertices[2] == targetVertices[i+2])  {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}




T3DFace& T3DFace::operator =(const T3DFace& face)
{
    if ( this == &face )
        return *this;

    m_tetrahedron   = face.m_tetrahedron;
    m_mateVertexPos = face.m_mateVertexPos;

    return *this;
}



