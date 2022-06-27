#include "MBSVertex.h"
#include "MBSEdge.h"
#include "MBSFace.h"

MBSEdge::MBSEdge()
: TopologicalEntity(),
  m_leftFace(rg_NULL), m_rightFace(rg_NULL),
  m_leftHand(rg_NULL), m_rightHand(rg_NULL),
  m_leftLeg(rg_NULL), m_rightLeg(rg_NULL),
  m_startVertex(rg_NULL), m_endVertex(rg_NULL),
  m_originalBetaEdge(rg_NULL), m_isArtificial(rg_FALSE),
  m_visited(rg_FALSE)
{
}



MBSEdge::MBSEdge(const rg_INT& ID)
: TopologicalEntity(ID),
  m_leftFace(rg_NULL), m_rightFace(rg_NULL),
  m_leftHand(rg_NULL), m_rightHand(rg_NULL),
  m_leftLeg(rg_NULL), m_rightLeg(rg_NULL),
  m_startVertex(rg_NULL), m_endVertex(rg_NULL),
  m_originalBetaEdge(rg_NULL), m_isArtificial(rg_FALSE),
  m_visited(rg_FALSE)
{
}



MBSEdge::MBSEdge(const rg_INT& ID, MBSVertex* startVertex, MBSVertex* endVertex)
: TopologicalEntity(ID),
  m_leftFace(rg_NULL), m_rightFace(rg_NULL),
  m_leftHand(rg_NULL), m_rightHand(rg_NULL),
  m_leftLeg(rg_NULL), m_rightLeg(rg_NULL),
  m_startVertex(startVertex), m_endVertex(endVertex),
  m_originalBetaEdge(rg_NULL), m_isArtificial(rg_FALSE),
  m_visited(rg_FALSE)
{
}



MBSEdge::MBSEdge(const MBSEdge& edge)
: TopologicalEntity(edge),
  m_leftFace(edge.m_leftFace), m_rightFace(edge.m_rightFace),
  m_leftHand(edge.m_leftHand), m_rightHand(edge.m_rightHand),
  m_leftLeg(edge.m_leftLeg), m_rightLeg(edge.m_rightLeg),
  m_startVertex(edge.m_startVertex), m_endVertex(edge.m_endVertex),
  m_originalBetaEdge(edge.m_originalBetaEdge), m_isArtificial(edge.m_isArtificial),
  m_visited(edge.m_visited)
{
}



MBSEdge::~MBSEdge()
{
}




MBSFace*     MBSEdge::getLeftFace() const
{
    return m_leftFace;
}



MBSFace*     MBSEdge::getRightFace() const
{
    return m_rightFace;
}



MBSEdge*     MBSEdge::getLeftHand() const
{
    return m_leftHand;
}



MBSEdge*     MBSEdge::getRightHand() const
{
    return m_rightHand;
}



MBSEdge*     MBSEdge::getLeftLeg() const
{
    return m_leftLeg;
}



MBSEdge*     MBSEdge::getRightLeg() const
{
    return m_rightLeg;
}



MBSVertex*   MBSEdge::getStartVertex() const
{
    return m_startVertex;
}



MBSVertex*   MBSEdge::getEndVertex() const
{
    return m_endVertex;
}



BetaEdge*   MBSEdge::getOriginalBetaEdge() const
{
    return m_originalBetaEdge;
}

rg_BOOL     MBSEdge::isArtificial() const
{
    return m_isArtificial;
}



MBSShell*    MBSEdge::getShell() const
{
    return m_startVertex->getShell();
}



void        MBSEdge::setLeftFace(MBSFace* leftFace)
{
    m_leftFace = leftFace;
}



void        MBSEdge::setRightFace(MBSFace* rightFace)
{
    m_rightFace = rightFace;
}



void        MBSEdge::setLeftHand(MBSEdge* leftHand)
{
    m_leftHand = leftHand;
}



void        MBSEdge::setRightHand(MBSEdge* rightHand)
{
    m_rightHand = rightHand;
}



void        MBSEdge::setLeftLeg(MBSEdge* leftLeg)
{
    m_leftLeg = leftLeg;
}



void        MBSEdge::setRightLeg(MBSEdge* rightLeg)
{
    m_rightLeg = rightLeg;
}



void        MBSEdge::setStartVertex(MBSVertex* startVertex)
{
    m_startVertex = startVertex;
}



void        MBSEdge::setEndVertex(MBSVertex* endVertex)
{
    m_endVertex = endVertex;
}



void        MBSEdge::setOriginalBetaEdge( BetaEdge* originalBetaEdge )
{
    m_originalBetaEdge = originalBetaEdge;
}


void        MBSEdge::isArtificial( const rg_BOOL& isArtificial )
{
    m_isArtificial = isArtificial;
}



MBSEdge&     MBSEdge::operator =(const MBSEdge& edge)
{
    if ( this == & edge ) {
        return *this;
    }

    TopologicalEntity::operator =(edge);

    m_leftFace    = edge.m_leftFace;
    m_rightFace   = edge.m_rightFace;

    m_leftHand    = edge.m_leftHand;
    m_rightHand   = edge.m_rightHand;
    m_leftLeg     = edge.m_leftLeg;
    m_rightLeg    = edge.m_rightLeg;

    m_startVertex = edge.m_startVertex;
    m_endVertex   = edge.m_endVertex;

    m_originalBetaEdge  = edge.m_originalBetaEdge;
    m_isArtificial      = edge.m_isArtificial;
    
    m_visited     = edge.m_visited;

    return *this;
}





void MBSEdge::searchBoundingVertices( rg_dList<MBSVertex*>& vertexList )
{
    if ( m_startVertex != rg_NULL ) {
        vertexList.add( m_startVertex );
    }

    if ( m_endVertex != rg_NULL ) {
        vertexList.add( m_endVertex );
    }
}



rg_FLAG MBSEdge::searchAdjacentEdges( rg_dList<MBSEdge*>& edgeList )
{
    MBSEdge* currEdge  = m_leftLeg;

    do {
        edgeList.add( currEdge );

        if ( m_startVertex == currEdge->getStartVertex() ) {
            currEdge = currEdge->getLeftLeg();
        }
        else if ( m_startVertex == currEdge->getEndVertex() ) {
            currEdge = currEdge->getRightHand();
        }
        else {
            return rg_FALSE;
        }
    } while (currEdge != this);


    currEdge  = m_rightHand;

    do {
        edgeList.add( currEdge );

        if ( m_endVertex == currEdge->getStartVertex() ) {
            currEdge = currEdge->getLeftLeg();
        }
        else if ( m_endVertex == currEdge->getEndVertex() ) {
            currEdge = currEdge->getRightHand();
        }
        else {
            return rg_FALSE;
        }
    } while (currEdge != this);


    return rg_TRUE;
}



void MBSEdge::searchAdjacentFaces( rg_dList<MBSFace*>& faceList )
{
    if ( m_leftFace != rg_NULL ) {
        faceList.add( m_leftFace );
    }

    if ( m_rightFace != rg_NULL ) {
        faceList.add( m_rightFace );
    }    
}
