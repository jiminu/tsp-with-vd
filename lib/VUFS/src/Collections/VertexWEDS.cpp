#include "VertexWEDS.h"
#include "EdgeWEDS.h"
#include "FaceWEDS.h"



VertexWEDS::VertexWEDS()
: m_firstEdge(rg_NULL), m_visited(rg_FALSE)
{
}



VertexWEDS::VertexWEDS(const rg_INT& ID)
: TopologicalEntity(ID), m_firstEdge(rg_NULL), m_visited(rg_FALSE)
{
}



VertexWEDS::VertexWEDS(const rg_INT& ID, EdgeWEDS* firstEdge)
: TopologicalEntity(ID), m_firstEdge(firstEdge), m_visited(rg_FALSE)
{
}



VertexWEDS::VertexWEDS(const VertexWEDS& vertex)
: TopologicalEntity(vertex), 
  m_firstEdge(vertex.m_firstEdge), m_visited(vertex.m_visited)
{
}



VertexWEDS::~VertexWEDS()
{
}





VertexWEDS& VertexWEDS::operator=(const VertexWEDS& vertex)
{
    if ( this == &vertex ) {
        return *this;
    }

    TopologicalEntity::operator =(vertex);
    m_firstEdge = vertex.m_firstEdge;
    m_visited   = vertex.m_visited;

    return *this;
}




//  Topological operators
rg_INT VertexWEDS::getNeighborVertices(  rg_dList<VertexWEDS*>& vertexList ) const
{
	EdgeWEDS* startEdge = m_firstEdge;
	EdgeWEDS* currEdge  = startEdge;

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
            break;
        }
    } while (startEdge != currEdge);

    return vertexList.getSize();
}



rg_INT VertexWEDS::getIncidentEdges( rg_dList<EdgeWEDS*>& edgeList ) const
{
	EdgeWEDS* startEdge = m_firstEdge;
	EdgeWEDS* currEdge  = startEdge;

    do {
        edgeList.add( currEdge );
        if ( this == currEdge->getStartVertex() ) {
            currEdge = currEdge->getLeftLeg();
        }
        else if ( this == currEdge->getEndVertex() ) {
            currEdge = currEdge->getRightHand();
        }
        else {
            break;
        }
    } while (startEdge != currEdge);

    return edgeList.getSize();
}



rg_INT VertexWEDS::getIncidentFaces( rg_dList<FaceWEDS*>& faceList ) const
{
	EdgeWEDS* startEdge = m_firstEdge;
	EdgeWEDS* currEdge  = startEdge;

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
            break;
        }
    } while (startEdge != currEdge);

    return faceList.getSize();
}




rg_INT VertexWEDS::getEdgesInStar( rg_dList<EdgeWEDS*>& edgeList ) const
{
    rg_dList<EdgeWEDS*> edgesIncidentToVertex;
    getIncidentEdges(edgesIncidentToVertex);

    rg_dList<EdgeWEDS*> edgesInShellOfVertex;
    getEdgesInShell(edgesInShellOfVertex);


    edgeList.append( edgesIncidentToVertex );
    edgeList.append( edgesInShellOfVertex );

    return edgeList.getSize();
}



rg_INT VertexWEDS::getEdgesInShell( rg_dList<EdgeWEDS*>& edgeList ) const
{
    rg_dList<EdgeWEDS*> incidentEdgeList;
    getIncidentEdges(incidentEdgeList);

    rg_INT numIncidentEdges = incidentEdgeList.getSize();
    rg_dNode< EdgeWEDS* >* currNode = incidentEdgeList.getFirstpNode();

    for ( rg_INT i=0; i<numIncidentEdges; i++, currNode=currNode->getNext() ) {
        EdgeWEDS* currEdge = currNode->getEntity();
        EdgeWEDS* nextEdge = currNode->getNext()->getEntity();

        EdgeWEDS* coverEdge = rg_NULL;
        if ( this == currEdge->getStartVertex() ) {
            coverEdge = currEdge->getLeftHand();
        }
        else {
            coverEdge = currEdge->getRightLeg();
        }

        if ( this == nextEdge->getStartVertex() ) {
            if ( coverEdge != nextEdge->getRightHand() ) {
                coverEdge = rg_NULL;
            }
        }
        else {
            if ( coverEdge != nextEdge->getLeftLeg() ) {
                coverEdge = rg_NULL;
            }
        }

        if ( coverEdge != rg_NULL ) {
            edgeList.add( coverEdge );
        }
    }

    return edgeList.getSize();
}




EdgeWEDS* VertexWEDS::findConnectingEdge( VertexWEDS* vertex ) const
{
	EdgeWEDS* connectingEdge = rg_NULL;
	EdgeWEDS* startEdge      = m_firstEdge;
	EdgeWEDS* currEdge       = startEdge;

    do {
        if ( this == currEdge->getStartVertex() ) {
            if ( vertex == currEdge->getEndVertex() ) {
                connectingEdge = currEdge;
                break;
            }
            currEdge = currEdge->getLeftLeg();
        }
        else { //if ( this == currEdge->getEndVertex() ) {
            if ( vertex == currEdge->getStartVertex() ) {
                connectingEdge = currEdge;
                break;
            }

            currEdge = currEdge->getRightHand();
        }
    } while (startEdge != currEdge);

    return connectingEdge;
}



rg_INT VertexWEDS::findSharingFace( VertexWEDS* vertex, rg_dList<FaceWEDS*>& faceList ) const
{
    rg_dList<FaceWEDS*> facesIncidentToThisVtx;
    getIncidentFaces(facesIncidentToThisVtx);

    facesIncidentToThisVtx.reset4Loop();
    while ( facesIncidentToThisVtx.setNext4Loop() ) {
        FaceWEDS* currFace = facesIncidentToThisVtx.getEntity();

        if ( currFace->isThere( vertex ) ) {
            faceList.add( currFace );
        }
    }

    return faceList.getSize();
}



