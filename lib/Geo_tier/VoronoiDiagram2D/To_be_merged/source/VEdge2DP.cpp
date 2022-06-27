#include "VEdge2DP.h"
#include "VVertex2DP.h"
#include "VFace2DP.h"

using namespace BULL2D::GeometryTier;

VEdge2DP::VEdge2DP()
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
    m_status        = WHITE_E_P;
}



VEdge2DP::VEdge2DP( const int& ID )
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
    m_status        = WHITE_E_P;
}



VEdge2DP::VEdge2DP( const int& ID, VVertex2DP* const startVertex, VVertex2DP* const endVertex )
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
    m_status        = WHITE_E_P;
}



VEdge2DP::VEdge2DP( const VEdge2DP& edge )
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
    m_status        = edge.m_status;
}



VEdge2DP::~VEdge2DP()
{
}



VEdge2DP& VEdge2DP::operator=( const VEdge2DP& edge )
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
    m_status        = edge.m_status;

    return *this;
}



VFace2DP* VEdge2DP::getMateFace( VVertex2DP* vertex ) const
{
    VFace2DP* mateFace = NULL;

    if( vertex == m_startVertex )
    {
        if( m_endVertex == m_leftHand->m_startVertex )
        {
            mateFace = m_leftHand->m_rightFace;
        }
        else
        {
            mateFace = m_leftHand->m_leftFace;
        }
    }
    else if( vertex == m_endVertex )
    {
        if( m_startVertex == m_leftLeg->m_startVertex )
        {
            mateFace = m_leftLeg->m_leftFace;
        }
        else
        {
            mateFace = m_leftLeg->m_rightFace;
        }
    }
    else
    {
        mateFace = NULL;
    }

    return mateFace;
}



void VEdge2DP::getIncidentVVertices( list<VVertex2DP*>& incidentVVerticesList ) const
{
    incidentVVerticesList.push_back( m_startVertex );
    incidentVVerticesList.push_back( m_endVertex );
}



void VEdge2DP::getIncidentVFaces( list<VFace2DP*>& incidentVFacesList ) const
{
    incidentVFacesList.push_back( m_leftFace );
    incidentVFacesList.push_back( m_rightFace );
}



void VEdge2DP::getAdjacentVEdges( list<VEdge2DP*>& adjacentVEdgesList ) const
{
    adjacentVEdgesList.push_back( m_rightHand );
    adjacentVEdgesList.push_back( m_leftHand );
    adjacentVEdgesList.push_back( m_leftLeg );
    adjacentVEdgesList.push_back( m_rightLeg );
}



VVertex2DP* VEdge2DP::getOppositeVVertex( VVertex2DP* vertex ) const
{
    if( vertex == m_startVertex )
    {
        return m_endVertex;
    }
    else if( vertex == m_endVertex )
    {
        return m_startVertex;
    }
    else
    {
        return NULL;
    }
}



VFace2DP* VEdge2DP::getOppositeVFace( VFace2DP* face ) const
{
    if( face == m_leftFace )
    {
        return m_rightFace;
    }
    else if( face == m_rightFace )
    {
        return m_leftFace;
    }
    else
    {
        return NULL;
    }
}



void VEdge2DP::setTopology( VVertex2DP* const startVertex, VVertex2DP* const endVertex, 
                           VFace2DP* const leftFace, VFace2DP* const rightFace, 
                           VEdge2DP* const leftHand, VEdge2DP* const rightHand, 
                           VEdge2DP* const leftLeg, VEdge2DP* const rightLeg )
{
    m_startVertex   = startVertex;
    m_endVertex     = endVertex;
    m_leftFace      = leftFace;
    m_rightFace     = rightFace;
    m_leftHand      = leftHand;
    m_rightHand     = rightHand;
    m_leftLeg       = leftLeg;
    m_rightLeg      = rightLeg;
}



void VEdge2DP::setTopology( VFace2DP* const leftFace, VFace2DP* const rightFace, 
                           VEdge2DP* const leftHand, VEdge2DP* const rightHand, 
                           VEdge2DP* const leftLeg, VEdge2DP* const rightLeg )
{
    m_leftFace      = leftFace;
    m_rightFace     = rightFace;
    m_leftHand      = leftHand;
    m_rightHand     = rightHand;
    m_leftLeg       = leftLeg;
    m_rightLeg      = rightLeg;
}



void VEdge2DP::setTopology( VEdge2DP* const leftHand, VEdge2DP* const rightHand, 
                           VEdge2DP* const leftLeg, VEdge2DP* const rightLeg )
{
    m_leftHand      = leftHand;
    m_rightHand     = rightHand;
    m_leftLeg       = leftLeg;
    m_rightLeg      = rightLeg;
}




bool VEdge2DP::isInfinite() const
{
    if( m_leftFace->isInfinite() || m_rightFace->isInfinite() ) {
        return true;
    }
    else {
        return false;
    }
}



bool VEdge2DP::isBounded() const
{
    if( isInfinite() ) {
        return false;
    }


    if( m_startVertex->isInfinite() || m_endVertex->isInfinite() ) {
        return false;
    }
    else {
        return true;
    }
}



bool VEdge2DP::isUnbounded() const
{
    if( isInfinite() ) {
        return false;
    }


    if( m_startVertex->isInfinite() || m_endVertex->isInfinite() ) {
        return true;
    }
    else {
        return false;
    }
}



bool VEdge2DP::isMemberOf( const VFace2DP* const face ) const
{
    if( m_leftFace == face || m_rightFace == face )
    {
        return true;
    }
    else
    {
        return false;
    }
}



bool VEdge2DP::isAdjacentTo( const VEdge2DP* const edge ) const
{
    if( m_leftHand == edge || m_rightHand == edge || m_leftLeg == edge || m_rightLeg == edge )
    {
        return true;
    }
    else
    {
        return false;
    }
}



bool VEdge2DP::isIncidentTo( const VVertex2DP* const vertex ) const
{
    if( m_startVertex == vertex || m_endVertex == vertex )
    {
        return true;
    }
    else
    {
        return false;
    }
}


