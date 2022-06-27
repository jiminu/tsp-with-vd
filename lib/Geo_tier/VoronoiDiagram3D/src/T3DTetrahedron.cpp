#include "T3DTetrahedron.h"
#include "T3DVertex.h"
using namespace V::GeometryTier;



T3DTetrahedron::T3DTetrahedron()
: m_check( rg_FALSE )
{
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i]   = rg_NULL;
        m_neighbor[i] = rg_NULL;
    }
}



T3DTetrahedron::T3DTetrahedron(const rg_INT& ID)
: TopologicalEntity(ID), m_check( rg_FALSE )
{
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i]   = rg_NULL;
        m_neighbor[i] = rg_NULL;
    }
}



T3DTetrahedron::T3DTetrahedron(const rg_INT& ID,
                               T3DVertex* vertex1, T3DVertex* vertex2, 
                               T3DVertex* vertex3, T3DVertex* vertex4)
: TopologicalEntity(ID), m_check( rg_FALSE )
{
    m_vertex[0] = vertex1;
    m_vertex[1] = vertex2;
    m_vertex[2] = vertex3;
    m_vertex[3] = vertex4;

    for( rg_INT i=0; i<4; i++ )  {
        m_neighbor[i] = rg_NULL;
    }
}




T3DTetrahedron::T3DTetrahedron(const T3DTetrahedron& tetrahedron)
{
	m_ID = tetrahedron.m_ID;
    m_check = tetrahedron.m_check;
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i]   = tetrahedron.m_vertex[i];
        m_neighbor[i] = tetrahedron.m_neighbor[i];
    }
}



T3DTetrahedron::~T3DTetrahedron()
{
}




T3DVertex*       T3DTetrahedron::getVertex(const rg_INDEX& i) const
{
    if ( i<0 || i>=4)
        return rg_NULL;
    else
        return m_vertex[i];
}



T3DVertex**      T3DTetrahedron::getVertices()
{
    return m_vertex;
}



T3DTetrahedron*  T3DTetrahedron::getNeighbor(const rg_INDEX& i) const
{
    if ( i<0 || i>=4)
        return rg_NULL;
    else
        return m_neighbor[i];
}



T3DTetrahedron** T3DTetrahedron::getNeighbors()
{
    return m_neighbor;
}



T3DTetrahedron*  T3DTetrahedron::getMateTetrahedron(const rg_INT& pos) const
{
    if ( pos<0 || pos>=4 )
        return rg_NULL;
    else
        return m_neighbor[pos];
}



T3DTetrahedron*  T3DTetrahedron::getMateTetrahedron(T3DVertex* vertex) const
{
    for (rg_INT i=0; i<4; i++ )  {
        if ( vertex == m_vertex[i] )
            return m_neighbor[i];
    }

    return rg_NULL;
}



T3DVertex* T3DTetrahedron::getMateVertex(T3DTetrahedron* neighbor) const
{
    for ( rg_INT i=0; i<4; i++ )  {
        if ( neighbor == m_neighbor[i] )
            return m_vertex[i];
    }

    return rg_NULL;
}



rg_INT T3DTetrahedron::getMateVertexPos(T3DTetrahedron* neighbor) const
{
    for ( rg_INT i=0; i<4; i++ )  {
        if ( neighbor == m_neighbor[i] )
            return i;
    }

    return -1;
}



rg_INT           T3DTetrahedron::getPosOfVertex(T3DVertex* vertex) const
{
    for (rg_INT i=0; i<4; i++ )  {
        if ( vertex == m_vertex[i] )
            return i;
    }

    return -1;
}



rg_INT T3DTetrahedron::getMateVertexPosOfIncidentFace(T3DTetrahedron* neighbor, const rg_INT& mateVertexPos)
{
    T3DVertex* incidentVertex[3];
    rg_INT     vIndex = 0;
    rg_INT i=0;
	for ( i=0; i<4; i++ )  {
        if ( i == mateVertexPos )
            continue;

        incidentVertex[vIndex++] = neighbor->m_vertex[i];
    }

    vIndex = 0;
    rg_INT  thisMateVertexPos = -1;

    for ( i=0; i<4; i++ )  {
        rg_FLAG isIncident = rg_FALSE;
        for ( rg_INT j=0; j<3; j++)  {        
            if ( m_vertex[i] == incidentVertex[j] )  {
                isIncident = rg_TRUE;
                break;
            }
        }

        if ( isIncident == rg_FALSE )  {
            return i;
        }
    }

    return thisMateVertexPos;
}



void T3DTetrahedron::getVerticesOnFace(const rg_INT& mateVertexPos, T3DVertex** vertices)
{
    switch( mateVertexPos )  {
        case 0:
            vertices[0] = m_vertex[1];
            vertices[1] = m_vertex[3];
            vertices[2] = m_vertex[2];
            break;
        case 1:
            vertices[0] = m_vertex[0];
            vertices[1] = m_vertex[2];
            vertices[2] = m_vertex[3];
            break;
        case 2:
            vertices[0] = m_vertex[0];
            vertices[1] = m_vertex[3];
            vertices[2] = m_vertex[1];
            break;
        case 3:
            vertices[0] = m_vertex[0];
            vertices[1] = m_vertex[1];
            vertices[2] = m_vertex[2];
            break;
        default:
            break;
    }
}



T3DFace T3DTetrahedron::getFace(const rg_INT& mateVertexPos)
{
    return T3DFace(this, mateVertexPos);
}



rg_INT T3DTetrahedron::locateNeighborPos(T3DTetrahedron* neighbor)
{
    for ( rg_INT i=0; i<4; i++)  {
        if ( neighbor == m_neighbor[i] )
            return i;
    }

    return -1;
}



rg_FLAG T3DTetrahedron::locateFace(T3DFace& faceToLocate, const rg_INT& vertexID1, const rg_INT& vertexID2, const rg_INT& vertexID3)
{
    rg_FLAG isLocated[4] = {rg_FALSE, rg_FALSE, rg_FALSE, rg_FALSE};
    rg_INT i=0;
	for ( i=0; i<4; i++ )  {
        rg_INT ID = m_vertex[i]->getID();
        if ( ID == vertexID1 || ID == vertexID2 || ID == vertexID3 )  {
            isLocated[i] = rg_TRUE;
        }
    }

    rg_INT mateVertexID = -1;
    rg_INT determinantForLocation = 0;
    for ( i=0; i<4; i++ )  {
        if ( isLocated[i] == rg_TRUE )
            determinantForLocation++;
        else
            mateVertexID = i;
    }

    if ( determinantForLocation == 3 )  {
        faceToLocate.setFace(this, mateVertexID);
        return rg_TRUE;
    }
    else
        return rg_FALSE;
}



rg_FLAG T3DTetrahedron::isThereThisVertex(T3DVertex* vertex) const
{
    for (rg_INT i=0; i<4; i++ )  {
        if ( vertex == m_vertex[i] )
            return rg_TRUE;
    }

    return rg_FALSE;
}



rg_FLAG T3DTetrahedron::isThereThisEdge(T3DVertex* vertex1, T3DVertex* vertex2) const
{
    if ( isThereThisVertex(vertex1) && isThereThisVertex(vertex2) )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_FLAG T3DTetrahedron::isThereThisNeighbor(T3DTetrahedron* neighbor) const
{
    for (rg_INT i=0; i<4; i++ )  {
        if ( neighbor == m_neighbor[i] )
            return rg_TRUE;
    }

    return rg_FALSE;
}



rg_FLAG T3DTetrahedron::isOnBoundaryOfConvexHull() const
{
    if (    ( m_neighbor[0] == rg_NULL )
         || ( m_neighbor[1] == rg_NULL )
         || ( m_neighbor[2] == rg_NULL )
         || ( m_neighbor[3] == rg_NULL ) )
         return rg_TRUE;
    else
        return rg_FALSE;
}



rg_FLAG T3DTetrahedron::isChecked() const
{
    return m_check;
}



void T3DTetrahedron::setCheck(const rg_FLAG& check)
{
    m_check = check;
}



void T3DTetrahedron::setVertex(const rg_INDEX& i, T3DVertex* vertex)
{
    if ( i<0 || i>=4)
        return;
    else
        m_vertex[i] = vertex;
}



void T3DTetrahedron::setVertices(T3DVertex* vertex1, T3DVertex* vertex2, 
                                 T3DVertex* vertex3, T3DVertex* vertex4)
{
    m_vertex[0] = vertex1;
    m_vertex[1] = vertex2;
    m_vertex[2] = vertex3;
    m_vertex[3] = vertex4;    
}



void T3DTetrahedron::setVertices(T3DVertex** vertices)
{
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i] = vertices[i];
    }
}



void T3DTetrahedron::setNeighbor(const rg_INDEX& i, T3DTetrahedron* neighbor)
{
    if ( i<0 || i>=4)
        return;
    else
        m_neighbor[i] = neighbor;
}



void T3DTetrahedron::setNeighbor(T3DVertex* vertex, T3DTetrahedron* neighbor)
{
    for ( rg_INT i=0; i<4; i++ )  {
        if ( m_vertex[i] == vertex )  {
            m_neighbor[i] = neighbor;
            break;
        }
    }
}



void T3DTetrahedron::setNeighbors(T3DTetrahedron* neighbor1, T3DTetrahedron* neighbor2, 
                                  T3DTetrahedron* neighbor3, T3DTetrahedron* neighbor4)
{
    m_neighbor[0] = neighbor1;
    m_neighbor[1] = neighbor2;
    m_neighbor[2] = neighbor3;
    m_neighbor[3] = neighbor4;
}



void T3DTetrahedron::setNeighbors(T3DTetrahedron** neighbors)
{
    for( rg_INT i=0; i<4; i++ )  {
        m_neighbor[i] = neighbors[i];
    }
}




T3DTetrahedron& T3DTetrahedron::operator =(const T3DTetrahedron& tetrahedron)
{
    if ( this == &tetrahedron )
        return *this;

    m_ID = tetrahedron.m_ID;
    m_check = tetrahedron.m_check;
    for( rg_INT i=0; i<4; i++ )  {
        m_vertex[i]   = tetrahedron.m_vertex[i];
        m_neighbor[i] = tetrahedron.m_neighbor[i];
    }

    return *this;
}



