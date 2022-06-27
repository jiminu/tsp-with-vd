#include "MBSVertex.h"
#include "MBSEdge.h"
#include "MBSFace.h"



MBSVertex::MBSVertex()
: TopologicalEntity(), 
  m_shell(rg_NULL), m_firstEdge(rg_NULL), m_visited(rg_FALSE),
  m_originalBetaVertex(rg_NULL), m_isArtificial(rg_FALSE)
{
}



MBSVertex::MBSVertex(const rg_INT& ID)
: TopologicalEntity(ID), 
  m_shell(rg_NULL), m_firstEdge(rg_NULL), m_visited(rg_FALSE),
  m_originalBetaVertex(rg_NULL), m_isArtificial(rg_FALSE)
{
}




MBSVertex::MBSVertex(const MBSVertex& vertex)
: TopologicalEntity(vertex), 
  m_shell(vertex.m_shell), m_firstEdge(vertex.m_firstEdge),m_visited(vertex.m_visited),
  m_originalBetaVertex(vertex.m_originalBetaVertex), m_isArtificial(vertex.m_isArtificial)
{
}



MBSVertex::~MBSVertex()
{
}






MBSEdge* MBSVertex::getFirstEdge() const
{
    return m_firstEdge;
}



BetaVertex* MBSVertex::getOriginalBetaVertex() const
{
    return m_originalBetaVertex;
}

rg_BOOL MBSVertex::isArtificial() const
{
    return m_isArtificial;
}



MBSShell* MBSVertex::getShell() const
{
    return m_shell;
}






void MBSVertex::setShell( MBSShell* shell )
{
    m_shell = shell;
}


void MBSVertex::setFirstEdge(MBSEdge* firstEdge)
{   
    m_firstEdge = firstEdge;
}



void MBSVertex::setOriginalBetaVertex(BetaVertex* originalBetaVertex)
{
    m_originalBetaVertex = originalBetaVertex;
}

void MBSVertex::isArtificial( const rg_BOOL& isArtificial )
{
    m_isArtificial = isArtificial;
}


MBSVertex& MBSVertex::operator =(const MBSVertex& vertex)
{
    if ( this == &vertex ) {
        return *this;
    }

    TopologicalEntity::operator =(vertex);

    m_shell     = vertex.m_shell;
    m_firstEdge = vertex.m_firstEdge;
    
    m_originalBetaVertex = vertex.m_originalBetaVertex;
    m_isArtificial       = vertex.m_isArtificial;

    m_visited   = vertex.m_visited;

    return *this;
}



rg_FLAG MBSVertex::searchAdjacentVertices(rg_dList<MBSVertex*>& vertexList) const
{
    MBSEdge* startEdge = m_firstEdge;
	MBSEdge* currEdge  = startEdge;

    do {
        if ( this == currEdge->getStartVertex() ) {
            vertexList.add( currEdge->getEndVertex() );
            currEdge = currEdge->getLeftLeg();
        }
        else if ( this == currEdge->getEndVertex() ) {
            vertexList.add( currEdge->getStartVertex() );
            currEdge = currEdge->getRightHand();
        }
        else {
            return rg_FALSE;
        }
    } while (startEdge != currEdge);

    return rg_TRUE;
}



rg_FLAG MBSVertex::searchIncidentEdges(rg_dList<MBSEdge*>& edgeList) const
{
    MBSEdge* startEdge = m_firstEdge;
	MBSEdge* currEdge  = startEdge;

    do {
        edgeList.add( currEdge );
        if ( this == currEdge->getStartVertex() )
            currEdge = currEdge->getLeftLeg();
        else if ( this == currEdge->getEndVertex() )
            currEdge = currEdge->getRightHand();
        else
            return rg_FALSE;
    } while (startEdge != currEdge);

    return rg_TRUE;
}



rg_FLAG MBSVertex::searchIncidentFaces(rg_dList<MBSFace*>& faceList) const
{
    MBSEdge* startEdge = m_firstEdge;
	MBSEdge* currEdge  = startEdge;

    do {
        if ( this == currEdge->getStartVertex() ) {
            if ( currEdge->getLeftFace() != rg_NULL ) {
                faceList.add( currEdge->getLeftFace() );
            }
            currEdge = currEdge->getLeftLeg();
        }
        else if ( this == currEdge->getEndVertex() ) {
            if ( currEdge->getRightFace() != rg_NULL ) {
                faceList.add( currEdge->getRightFace() );
            }
            currEdge = currEdge->getRightHand();
        }
        else {
            return rg_FALSE;
        }
    } while (startEdge != currEdge);

    return rg_TRUE;
}


