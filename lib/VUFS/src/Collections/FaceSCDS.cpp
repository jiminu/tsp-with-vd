#include "FaceSCDS.h"

FaceSCDS::FaceSCDS()
: m_visited(rg_FALSE)
{
    for ( rg_INT i=0; i<3; i++ ) {
        m_vertex[i]   = rg_NULL;
        m_neighbor[i] = rg_NULL;
    }
}



FaceSCDS::FaceSCDS(const rg_INT& ID)
: TopologicalEntity(ID), m_visited(rg_FALSE)
{
    for ( rg_INT i=0; i<3; i++ ) {
        m_vertex[i]   = rg_NULL;
        m_neighbor[i] = rg_NULL;
    }
}



FaceSCDS::FaceSCDS(const FaceSCDS& face)
: TopologicalEntity(face), m_visited(face.m_visited)
{
    for ( rg_INT i=0; i<3; i++ ) {
        m_vertex[i]   = face.m_vertex[i];
        m_neighbor[i] = face.m_neighbor[i];
    }
}



FaceSCDS::~FaceSCDS()
{
}



VertexSCDS* FaceSCDS::getVertex(const rg_INDEX& i) const
{
    if ( i >= 0 && i < 3 ) {
        return m_vertex[i];
    }
    else {
        return rg_NULL;
    }
}



FaceSCDS* FaceSCDS::getNeighbor(const rg_INDEX& i) const
{
    if ( i >= 0 && i < 3 ) {
        return m_neighbor[i];
    }
    else {
        return rg_NULL;
    }
}



//VertexSCDS** FaceSCDS::getVertices()
//{
//    return m_vertex;
//}
//
//
//
//FaceSCDS** FaceSCDS::getNeighbors()
//{
//    return m_neighbor;
//}



void FaceSCDS::setVertex(  const rg_INDEX& i, VertexSCDS* vertex)
{
    if ( i >= 0 && i < 3 ) {
        m_vertex[i] = vertex;
    }
}



void FaceSCDS::setNeighbor(const rg_INDEX& i, FaceSCDS*   neighbor)
{
    if ( i >= 0 && i < 3 ) {
        m_neighbor[i] = neighbor;
    }
}




void FaceSCDS::setVertices(VertexSCDS* vertex1,   VertexSCDS* vertex2,   VertexSCDS* vertex3)
{
    m_vertex[0] = vertex1;
    m_vertex[1] = vertex2;
    m_vertex[2] = vertex3;
}



void FaceSCDS::setNeighbors(FaceSCDS*  neighbor1, FaceSCDS*   neighbor2, FaceSCDS*   neighbor3)
{
    m_neighbor[0] = neighbor1;
    m_neighbor[1] = neighbor2;
    m_neighbor[2] = neighbor3;
}




FaceSCDS& FaceSCDS::operator =(const FaceSCDS& face)
{
    if ( this == &face ) {
        return *this;
    }

    TopologicalEntity::operator =(face);
    
    for ( rg_INT i=0; i<3; i++ ) {
        m_vertex[i]   = face.m_vertex[i];
        m_neighbor[i] = face.m_neighbor[i];
    }

    m_visited = face.m_visited;

    return *this;
}



FaceSCDS* FaceSCDS::getMateFace(const rg_INT& pos) const
{
    if ( pos >= 0 && pos < 3 ) {
        return m_neighbor[pos];
    }
    else {
        return rg_NULL;
    }
}



FaceSCDS* FaceSCDS::getMateFace(VertexSCDS* vertex) const
{
    for (rg_INT i=0; i<3; i++ )  {
        if ( vertex == m_vertex[i] ) {
            return m_neighbor[i];
        }
    }

    return rg_NULL;
}



rg_INT FaceSCDS::getPosOfVertex(VertexSCDS* vertex) const
{
    for (rg_INT i=0; i<3; i++ )  {
        if ( vertex == m_vertex[i] ) {
            return i;
        }
    }

    return -1;
}



rg_INT FaceSCDS::getPosOfNeighbor(FaceSCDS* neighbor) const
{
    for (rg_INT i=0; i<3; i++ )  {
        if ( neighbor == m_neighbor[i] ) {
            return i;
        }
    }

    return -1;
}


rg_INT FaceSCDS::getThisPosInNeighbor(const rg_INT& neighborPos, FaceSCDS* neighbor) const
{
    if ( m_neighbor[neighborPos] != neighbor ) {
        return -1;
    }

    VertexSCDS* sharingVtx[2] = {rg_NULL, rg_NULL};
    switch ( neighborPos ) {
        case 0:
            sharingVtx[0] = m_vertex[1];
            sharingVtx[1] = m_vertex[2];
            break;
        case 1:
            sharingVtx[0] = m_vertex[0];
            sharingVtx[1] = m_vertex[2];
            break;
        case 2:
            sharingVtx[0] = m_vertex[0];
            sharingVtx[1] = m_vertex[1];
            break;
        default:
            return -1;
    }

    rg_INT posInNeighbor = -1;
    for ( rg_INT i=0; i<3; i++ ) {
        if (    neighbor->m_vertex[i] != sharingVtx[0]
             && neighbor->m_vertex[i] != sharingVtx[1] ) {
             posInNeighbor = i;
             break;
        }
    }

    return posInNeighbor;
}


rg_BOOL FaceSCDS::isThere(VertexSCDS* vertex) const
{
    for (rg_INT i=0; i<3; i++ )  {
        if ( vertex == m_vertex[i] ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_BOOL FaceSCDS::isThere(FaceSCDS*   neighbor) const
{
    for (rg_INT i=0; i<3; i++ )  {
        if ( neighbor == m_neighbor[i] ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



void FaceSCDS::getBoundingVertices(rg_dList<VertexSCDS*> vertexList) const
{
    for (rg_INT i=0; i<3; i++ )  {
        vertexList.add( m_vertex[i] );
    }
}



void FaceSCDS::getIncidentNeighbors(rg_dList<FaceSCDS*>  neighborList) const
{
    for (rg_INT i=0; i<3; i++ )  {
        neighborList.add( m_neighbor[i] );
    }
}





