#include "QTEdge.h"
#include "QTFace.h"
#include "QTVertex.h"
using namespace BULL2D::GeometryTier;



QTEdge::QTEdge(void)
{
    m_ID            = -1;
    m_startVertex   = NULL;
    m_endVertex     = NULL;
    m_leftFace      = NULL;
    m_rightFace     = NULL;
    m_leftHand      = NULL;
    m_rightHand     = NULL;
    m_leftLeg       = NULL;
    m_rightLeg      = NULL;
}



QTEdge::QTEdge( const int& ID )
{
    m_ID            = ID;
    m_startVertex   = NULL;
    m_endVertex     = NULL;
    m_leftFace      = NULL;
    m_rightFace     = NULL;
    m_leftHand      = NULL;
    m_rightHand     = NULL;
    m_leftLeg       = NULL;
    m_rightLeg      = NULL;
}



QTEdge::QTEdge( const int& ID, QTVertex* const startVertex, QTVertex* const endVertex )
{
    m_ID            = ID;
    m_startVertex   = startVertex;
    m_endVertex     = endVertex;
    m_leftFace      = NULL;
    m_rightFace     = NULL;
    m_leftHand      = NULL;
    m_rightHand     = NULL;
    m_leftLeg       = NULL;
    m_rightLeg      = NULL;
}



QTEdge::QTEdge( const QTEdge& edge )
{
    m_ID            = edge.m_ID;
    m_startVertex   = edge.m_startVertex;
    m_endVertex     = edge.m_endVertex;
    m_leftFace      = edge.m_leftFace;
    m_rightFace     = edge.m_rightFace;
    m_leftHand      = edge.m_leftHand;
    m_rightHand     = edge.m_rightHand;
    m_leftLeg       = edge.m_leftLeg;
    m_rightLeg      = edge.m_rightLeg;
}



QTEdge::~QTEdge(void)
{
}



QTEdge& QTEdge::operator=( const QTEdge& edge )
{
    if( this == &edge ) {
        return *this;
    }

    m_ID            = edge.m_ID;
    m_startVertex   = edge.m_startVertex;
    m_endVertex     = edge.m_endVertex;
    m_leftFace      = edge.m_leftFace;
    m_rightFace     = edge.m_rightFace;
    m_leftHand      = edge.m_leftHand;
    m_rightHand     = edge.m_rightHand;
    m_leftLeg       = edge.m_leftLeg;
    m_rightLeg      = edge.m_rightLeg;

    return *this;
}



bool QTEdge::operator==( const QTEdge& edge ) const
{
    if( this == &edge ) {
        return true;
    }

    if( m_ID            == edge.m_ID            &&
        m_startVertex   == edge.m_startVertex   &&
        m_endVertex     == edge.m_endVertex     &&
        m_leftFace      == edge.m_leftFace      &&
        m_rightFace     == edge.m_rightFace     &&
        m_leftHand      == edge.m_leftHand      &&
        m_rightHand     == edge.m_rightHand     &&
        m_leftLeg       == edge.m_leftLeg       &&
        m_rightLeg      == edge.m_rightLeg )
    {
        return true;
    }
    else
    {
        return false;
    }
}



void QTEdge::getIncidentVertices( list<QTVertex*> incidentVertices )
{
    incidentVertices.push_back( m_startVertex );
    incidentVertices.push_back( m_endVertex );
}



void QTEdge::getIncidentFaces( list<QTFace*> incidentFaces )
{
    incidentFaces.push_back( m_leftFace );
    incidentFaces.push_back( m_rightFace );
}



void QTEdge::getAdjacentEdges( list<QTEdge*> adjacentEdges )
{
    adjacentEdges.push_back( m_rightHand );
    adjacentEdges.push_back( m_leftHand );
    adjacentEdges.push_back( m_leftLeg );
    adjacentEdges.push_back( m_rightLeg );
}



bool QTEdge::isGoingToInfinite() const
{
    if( m_startVertex->isInfinite() || m_endVertex->isInfinite() ) {
        return true;
    }
    else {
        return false;
    }
}
