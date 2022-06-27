#include "QTFace.h"
#include "QTVertex.h"
#include "QTEdge.h"
using namespace BULL2D::GeometryTier;



QTFace::QTFace(void)
{
    m_ID        = -1;
    m_firstEdge = NULL;
}



QTFace::QTFace( const int& ID )
{
    m_ID        = ID;
    m_firstEdge = NULL;
}



QTFace::QTFace( const int& ID, QTEdge* const firstEdge )
{
    m_ID        = ID;
    m_firstEdge = NULL;
}



QTFace::QTFace( const QTFace& face )
{
    m_ID        = face.m_ID;
    m_firstEdge = face.m_firstEdge;
}



QTFace::~QTFace(void)
{
}

QTFace& QTFace::operator=( const QTFace& face )
{
    if( this == &face ) {
        return *this;
    }

    m_ID        = face.m_ID;
    m_firstEdge = face.m_firstEdge;

    return *this;
}

bool QTFace::operator==( const QTFace& face ) const
{
    if( this == &face ) {
        return true;
    }

    if( m_ID        == face.m_ID &&
        m_firstEdge == face.m_firstEdge )
    {
        return true;
    }
    else {
        return false;
    }
}



void QTFace::getBoundaryEdges( list<QTEdge*>& boundaryEdges ) const
{
    QTEdge* currEdge = m_firstEdge;

    do {
        boundaryEdges.push_back( currEdge );

        if( this == currEdge->getLeftFace() )
        {
            currEdge = currEdge->getLeftHand();
        }
        else
        {
            currEdge = currEdge->getRightLeg();
        }
    } while( currEdge != m_firstEdge );
}



void QTFace::getBoundaryVertices( list<QTVertex*>& boundaryVertices ) const
{
    QTEdge* currEdge = m_firstEdge;

    do {
        if( this == currEdge->getLeftFace() )
        {
            boundaryVertices.push_back( currEdge->getStartVertex() );
            currEdge = currEdge->getLeftHand();
        }
        else
        {
            boundaryVertices.push_back( currEdge->getEndVertex() );
            currEdge = currEdge->getRightLeg();
        }
    } while( currEdge != m_firstEdge );
}



void QTFace::getAdjacentFaces( list<QTFace*>& adjacentFaces ) const
{
    QTEdge* currEdge = m_firstEdge;

    do {
        if( this == currEdge->getLeftFace() )
        {
            adjacentFaces.push_back( currEdge->getRightFace() );
            currEdge = currEdge->getLeftHand();
        }
        else
        {
            adjacentFaces.push_back( currEdge->getLeftFace() );
            currEdge = currEdge->getRightLeg();
        }
    } while( currEdge != m_firstEdge );
}



bool QTFace::isInfinite() const
{
    list<QTVertex*> boundaryVertices;
    getBoundaryVertices( boundaryVertices );

    bool ret = false;

    list<QTVertex*>::iterator i_vertex = boundaryVertices.begin();
    for( ; i_vertex != boundaryVertices.end() ; i_vertex++ ) {
        QTVertex* currVertex = *i_vertex;
        if( currVertex->isInfinite() ) {
            ret = true;
            break;
        }
    }

    return ret;
}
