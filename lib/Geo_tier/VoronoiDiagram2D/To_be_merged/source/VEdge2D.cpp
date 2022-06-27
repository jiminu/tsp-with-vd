#include "VEdge2D.h"
#include "VVertex2D.h"
#include "VFace2D.h"
using namespace BULL2D::GeometryTier;



VEdge2D::VEdge2D()
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
    m_status        = WHITE_E;
    m_isAnomalyTestDone = false;
    m_isAlreadyCandidateForFlippingInPhantomRemoval = false;
}



VEdge2D::VEdge2D( const int& ID )
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
    m_status        = WHITE_E;
    m_isAnomalyTestDone = false;
    m_isAlreadyCandidateForFlippingInPhantomRemoval = false;
}



VEdge2D::VEdge2D( const int& ID, VVertex2D* const startVertex, VVertex2D* const endVertex )
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
    m_status        = WHITE_E;
    m_isAnomalyTestDone = false;
    m_isAlreadyCandidateForFlippingInPhantomRemoval = false;
}



VEdge2D::VEdge2D( const VEdge2D& edge )
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
    m_isAnomalyTestDone = edge.m_isAnomalyTestDone;
    m_isAlreadyCandidateForFlippingInPhantomRemoval = edge.m_isAlreadyCandidateForFlippingInPhantomRemoval;
}



VEdge2D::~VEdge2D()
{
}



void VEdge2D::setGeometry( const rg_Point2D& sp, const rg_Point2D& tvs, const rg_Point2D& ep, const rg_Point2D& tve, const rg_Point2D& passPt )
{
    m_geometry.makeRQBezier(sp, tvs, ep, tve, passPt);

    if( m_geometry.getWeight(1) < 0 )
    {
        m_geometry.setWeight( 1, m_geometry.getWeight(1) * -1 );
    }
}



void VEdge2D::setGeometryLine()
{
    rg_Point2D ctrlPt0 = m_startVertex->getLocation();
    rg_Point2D ctrlPt2 = m_endVertex->getLocation();
    rg_Point2D ctrlPt1 = (ctrlPt0 + ctrlPt2)/2.;

    m_geometry.setCtrlPt( 0, ctrlPt0 );
    m_geometry.setCtrlPt( 1, ctrlPt1 );
    m_geometry.setCtrlPt( 2, ctrlPt2 );
    m_geometry.setWeight(0, 1.);
    m_geometry.setWeight(1, 1.);
    m_geometry.setWeight(2, 1.);
}



VEdge2D& VEdge2D::operator=( const VEdge2D& edge )
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
    m_isAnomalyTestDone = edge.m_isAnomalyTestDone;
    m_isAlreadyCandidateForFlippingInPhantomRemoval = edge.m_isAlreadyCandidateForFlippingInPhantomRemoval;

    return *this;
}




bool VEdge2D::operator==( const VEdge2D& edge ) const
{
    if( this == &edge ) {
        return true;
    }

    if( m_rightHand == edge.m_rightHand &&
        m_leftHand == edge.m_leftHand &&
        m_rightLeg == edge.m_rightLeg &&
        m_leftLeg == edge.m_leftLeg &&
        m_leftFace == edge.m_leftFace &&
        m_rightFace == edge.m_rightFace &&
        m_startVertex == edge.m_startVertex &&
        m_endVertex == edge.m_endVertex &&
        m_ID == edge.m_ID &&
        m_status == edge.m_status &&
        m_isAnomalyTestDone == edge.m_isAnomalyTestDone &&
        m_isAlreadyCandidateForFlippingInPhantomRemoval == edge.m_isAlreadyCandidateForFlippingInPhantomRemoval )
    {
        return true;
    }
    else {
        return false;
    }
}



void VEdge2D::getIncidentVVertices( list<VVertex2D*>& incidentVVerticesList ) const
{
    incidentVVerticesList.push_back( m_startVertex );
    incidentVVerticesList.push_back( m_endVertex );
}



void VEdge2D::getIncidentVFaces( list<VFace2D*>& incidentVFacesList ) const
{
    incidentVFacesList.push_back( m_leftFace );
    incidentVFacesList.push_back( m_rightFace );
}



void VEdge2D::getAdjacentVEdges( list<VEdge2D*>& adjacentVEdgesList ) const
{
    adjacentVEdgesList.push_back( m_rightHand );
    adjacentVEdgesList.push_back( m_leftHand );
    adjacentVEdgesList.push_back( m_leftLeg );
    adjacentVEdgesList.push_back( m_rightLeg );
}




VVertex2D* VEdge2D::getOppositVVertex( VVertex2D* vertex ) const
{
    if( vertex == m_startVertex ) {
        return m_endVertex;
    }
    else if( vertex == m_endVertex ) {
        return m_startVertex;
    }
    else {
        return NULL;
    }
}



VFace2D* VEdge2D::getOppositeVFace( VFace2D* face ) const
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



void VEdge2D::setTopology( VVertex2D* const startVertex, VVertex2D* const endVertex, 
                           VFace2D* const leftFace, VFace2D* const rightFace, 
                           VEdge2D* const leftHand, VEdge2D* const rightHand, 
                           VEdge2D* const leftLeg, VEdge2D* const rightLeg )
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



void VEdge2D::setTopology( VFace2D* const leftFace, VFace2D* const rightFace, 
                           VEdge2D* const leftHand, VEdge2D* const rightHand, 
                           VEdge2D* const leftLeg, VEdge2D* const rightLeg )
{
    m_leftFace      = leftFace;
    m_rightFace     = rightFace;
    m_leftHand      = leftHand;
    m_rightHand     = rightHand;
    m_leftLeg       = leftLeg;
    m_rightLeg      = rightLeg;
}



void VEdge2D::setTopology( VEdge2D* const leftHand, VEdge2D* const rightHand, 
                           VEdge2D* const leftLeg, VEdge2D* const rightLeg )
{
    m_leftHand      = leftHand;
    m_rightHand     = rightHand;
    m_leftLeg       = leftLeg;
    m_rightLeg      = rightLeg;
}



bool VEdge2D::isInfinite() const
{
    if( m_leftFace->isInfinite() || m_rightFace->isInfinite() ) {
        return true;
    }
    else {
        return false;
    }
}



bool VEdge2D::isBounded() const
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



bool VEdge2D::isUnBounded() const
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



void VEdge2D::flip()
{
    VEdge2D*    RH = m_rightHand;
    VEdge2D*    LH = m_leftHand;
    VEdge2D*    RL = m_rightLeg;
    VEdge2D*    LL = m_leftLeg;

    VVertex2D*  SV = m_startVertex;
    VVertex2D*  EV = m_endVertex;

    VFace2D*    RF = m_rightFace;
    VFace2D*    LF = m_leftFace;

    //move first edge on two faces
    if( this == m_leftFace->getFirstVEdge() )	{
        if( LH != NULL) {
            m_leftFace->setFirstVEdge( LH );
        }
        else if(LL != NULL) {
            m_leftFace->setFirstVEdge( LL );
        }
        else {
            return;
        }
    }

    if( this == m_rightFace->getFirstVEdge() ) {
        if(RH != NULL) {
            m_rightFace->setFirstVEdge( RH );
        }
        else if(RL != NULL) {
            m_rightFace->setFirstVEdge( RL );
        }
        else {
            return; //error
        }
    }

    //move first edge on two vertices (start/end)
    SV->setFirstVEdge( this );
    EV->setFirstVEdge( this );

    //change vertex(start or end) of incident edges
    //left hand and right leg does not change start or end vertex
    if(RH->getStartVertex() == EV)
    {
        RH->setStartVertex( SV );

        RH->setRightLeg( RL );
        RH->setLeftLeg( this );
    }
    else //end
    {
        RH->setEndVertex( SV );

        RH->setRightHand( this );
        RH->setLeftHand( RL );
    }

    if(LL->getStartVertex() == SV)
    {
        LL->setStartVertex( EV );

        LL->setRightLeg( LH );
        LL->setLeftLeg( this );
    }
    else
    {
        LL->setEndVertex( EV );

        LL->setRightHand( this );
        LL->setLeftHand( LH );
    }

    if( RL->getStartVertex() == SV )
    {
        RL->setRightLeg( this );
        RL->setLeftLeg( RH );
    }
    else
    {
        RL->setRightHand( RH );
        RL->setLeftHand( this );
    }

    if( LH->getStartVertex() == EV )
    {
        LH->setRightLeg( this );
        LH->setLeftLeg( LL );
    }
    else
    {
        LH->setRightHand( LL );
        LH->setLeftHand( this );
    }

    //change right and left face
    if( LF == LH->getLeftFace() ) {
        m_rightFace = LH->getRightFace();
    }
    else {
        m_rightFace = LH->getLeftFace();
    }

    if( RF == RL->getRightFace() ) {
        m_leftFace = RL->getLeftFace();
    }
    else {
        m_leftFace = RL->getRightFace();
    }

    //update topology of two legs and two hands
    m_rightHand = LH;
    m_leftHand  = LL;
    m_rightLeg  = RH;
    m_leftLeg   = RL;
}



VFace2D* VEdge2D::getMateFace( VVertex2D* vertex ) const
{
    VFace2D* mateFace = rg_NULL;

    if ( vertex == getStartVertex() ) {
        if ( getEndVertex() == getLeftHand()->getStartVertex() ) {
            mateFace = getLeftHand()->getRightFace();
        }
        else {
            mateFace = getLeftHand()->getLeftFace();
        }
    }
    else if ( vertex == getEndVertex() ) {
        if ( getStartVertex() == getLeftLeg()->getStartVertex() ) {
            mateFace = getLeftLeg()->getLeftFace();
        }
        else {
            mateFace = getLeftLeg()->getRightFace();
        }
    }
    else {
        mateFace = rg_NULL;
    }

    return mateFace;
}



void VEdge2D::get4GeneratorsDefineEdgeEquationAndStartAndEndVertices( list<Generator2D*>& generatorList )
{
    generatorList.push_back( getRightFace()->getGenerator() );
    generatorList.push_back( getMateFace( m_startVertex )->getGenerator() );
    generatorList.push_back( getLeftFace()->getGenerator() );
    generatorList.push_back( getMateFace( m_endVertex )->getGenerator() );    
}



bool VEdge2D::isMemberOf( const VFace2D* const face ) const
{
//     bool isMemberOfFace = false;
// 
//     if( m_leftFace == face || m_rightFace == face )
//     {
//         isMemberOfFace = true;
//     }
// 
//     return isMemberOfFace;

    if( m_leftFace == face || m_rightFace == face )
    {
        return true;
    }
    else
    {
        return false;
    }
}



bool VEdge2D::isAdjacentTo( const VEdge2D* const edge ) const
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



bool VEdge2D::isIncidentTo( const VVertex2D* const vertex ) const
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


