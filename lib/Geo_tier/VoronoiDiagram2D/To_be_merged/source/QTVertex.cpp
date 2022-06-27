#include "QTVertex.h"
#include "QTEdge.h"
#include "QTFace.h"
using namespace BULL2D::GeometryTier;



QTVertex::QTVertex(void)
{
    m_ID        = -1;
    m_firstEdge = NULL;
    m_circle    = NULL;
}



QTVertex::QTVertex( const int& ID )
{
    m_ID        = ID;
    m_firstEdge = NULL;
    m_circle    = NULL;
}



QTVertex::QTVertex( const int& ID, QTEdge* const firstEdge )
{
    m_ID        = ID;
    m_firstEdge = firstEdge;
    m_circle    = NULL;
}



QTVertex::QTVertex( const int& ID, rg_Circle2D* const circle )
{
    m_ID        = ID;
    m_firstEdge = NULL;
    m_circle    = circle;
}



QTVertex::QTVertex( const int& ID, rg_Circle2D* const circle, QTEdge* const firstEdge )
{
    m_ID        = ID;
    m_firstEdge = firstEdge;
    m_circle    = circle;
}



QTVertex::QTVertex( const QTVertex& vertex )
{
    m_ID        = vertex.m_ID;
    m_firstEdge = vertex.m_firstEdge;
    m_circle    = vertex.m_circle;
}



QTVertex::~QTVertex(void)
{
}



QTVertex QTVertex::operator=( const QTVertex& vertex )
{
    if( this == &vertex ) {
        return *this;
    }

    m_ID        = vertex.m_ID;
    m_firstEdge = vertex.m_firstEdge;
    m_circle    = vertex.m_circle;

    return *this;
}



bool QTVertex::operator==( const QTVertex& vertex ) const
{
    if( this == &vertex ) {
        return true;
    }

    if( m_ID        == vertex.m_ID        &&
        m_firstEdge == vertex.m_firstEdge &&
        m_circle    == vertex.m_circle )
    {
        return true;
    }
    else {
        return false;
    }
}



void QTVertex::getIncidentEdges( list<QTEdge*>& incidentEdges ) const
{
    QTEdge* currEdge = m_firstEdge;

    do {
        incidentEdges.push_back( currEdge );

        if( this == currEdge->getStartVertex() )
        {
            currEdge = currEdge->getLeftLeg();
        }
        else
        {
            currEdge = currEdge->getRightHand();
        }
    } while( currEdge != m_firstEdge );
}



void QTVertex::getIncidentFaces( list<QTFace*>& incidentFaces ) const
{
    QTEdge* currEdge = m_firstEdge;

    do {
        if( this == currEdge->getStartVertex() )
        {
            incidentFaces.push_back( currEdge->getLeftFace() );
            currEdge = currEdge->getLeftLeg();
        }
        else
        {
            incidentFaces.push_back( currEdge->getRightFace() );
            currEdge = currEdge->getRightHand();
        }
    } while( currEdge != m_firstEdge );
}



void QTVertex::getAdjacentVertices( list<QTVertex*>& adjacentVertices ) const
{
    QTEdge* currEdge = m_firstEdge;

    do {
        if( this == currEdge->getStartVertex() )
        {
            adjacentVertices.push_back( currEdge->getEndVertex() );
            currEdge = currEdge->getLeftLeg();
        }
        else
        {
            adjacentVertices.push_back( currEdge->getStartVertex() );
            currEdge = currEdge->getRightHand();
        }
    } while( currEdge != m_firstEdge );
}



bool QTVertex::isInfinite() const
{
    if( m_circle == NULL ) {
        return true;
    }
    else {
        return false;
    }
}
