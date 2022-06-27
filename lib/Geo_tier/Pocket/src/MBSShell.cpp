#include "MBSShell.h"
#include "MBSFace.h"
#include "MBSEdge.h"
#include "MBSVertex.h"



MBSShell::MBSShell()
: TopologicalEntity(),
  m_body(rg_NULL),
  m_visited(rg_FALSE)
{
}



MBSShell::MBSShell(const rg_INT& ID)
: TopologicalEntity(ID),
  m_body(rg_NULL),
  m_visited(rg_FALSE)
{
}



MBSShell::MBSShell(const rg_INT& ID, MBSBody* body)
: TopologicalEntity(ID),
  m_body(body),
  m_visited(rg_FALSE)
{
}



MBSShell::MBSShell(const MBSShell& shell)
: TopologicalEntity(shell.m_ID),
  m_body(shell.m_body),
  m_visited(shell.m_visited)
{
    m_faces.duplicateList(shell.m_faces);
    m_edges.duplicateList(shell.m_edges);
    m_vertices.duplicateList(shell.m_vertices);
}



MBSShell::~MBSShell()
{    
    rg_INT i = 0;
    rg_INT numFaces = m_faces.getSize();
    rg_dNode<MBSFace*>* faceNode = m_faces.getFirstpNode();
    for ( i=0; i<numFaces; i++, faceNode=faceNode->getNext() ) {
        MBSFace* currFace = faceNode->getEntity();
        if ( currFace != rg_NULL ) {
            delete currFace;
            
            faceNode->setEntity( rg_NULL );
        }
    }
    m_faces.removeAll();



    rg_INT numEdges = m_edges.getSize();
    rg_dNode<MBSEdge*>* edgeNode = m_edges.getFirstpNode();
    for ( i=0; i<numEdges; i++, edgeNode=edgeNode->getNext() ) {
        MBSEdge* currEdge = edgeNode->getEntity();
        if ( currEdge != rg_NULL ) {
            delete currEdge;
            
            edgeNode->setEntity( rg_NULL );
        }
    }
    m_edges.removeAll();


    rg_INT numVertices = m_vertices.getSize();
    rg_dNode<MBSVertex*>* vertexNode = m_vertices.getFirstpNode();
    for ( i=0; i<numVertices; i++, vertexNode=vertexNode->getNext() ) {
        MBSVertex* currVtx = vertexNode->getEntity();
        if ( currVtx != rg_NULL ) {
            delete currVtx;
            
            vertexNode->setEntity( rg_NULL );
        }
    }
    m_vertices.removeAll();


//     MBSFace* currFace = rg_NULL;
//     m_faces.reset4Loop();
//     while ( m_faces.setNext4Loop() ) {
//         currFace = m_faces.getEntity();
// 
//         if ( currFace != rg_NULL ) {
//             delete currFace;
//             currFace = rg_NULL;
//         }
//     }
// 
//     MBSLoop* currLoop = rg_NULL;
//     m_loops.reset4Loop();
//     while ( m_loops.setNext4Loop() ) {
//         currLoop = m_loops.getEntity();
// 
//         if ( currLoop != rg_NULL ) {
//             delete currLoop;
//             currLoop = rg_NULL;
//         }
//     }
// 
//     MBSEdge* currEdge = rg_NULL;
//     m_edges.reset4Loop();
//     while ( m_edges.setNext4Loop() ) {
//         currEdge = m_edges.getEntity();
// 
//         if ( currEdge != rg_NULL ) {
//             delete currEdge;
//             currEdge = rg_NULL;
//         }
//     }
// 
//     MBSVertex* currVtx = rg_NULL;
//     m_vertices.reset4Loop();
//     while ( m_vertices.setNext4Loop() ) {
//         currVtx = m_vertices.getEntity();
// 
//         if ( currVtx != rg_NULL ) {
//             delete currVtx;
//             currVtx = rg_NULL;
//         }
//     }
}




MBSBody*  MBSShell::getBody() const
{
    return m_body;
}




rg_INT MBSShell::getNumFaces() const
{
    return m_faces.getSize();
}





rg_INT MBSShell::getNumEdges() const
{
    return m_edges.getSize();
}



rg_INT MBSShell::getNumVertices() const
{
    return m_vertices.getSize();
}




rg_dList<MBSFace*>*  MBSShell::getFaces() 
{
    return &m_faces;
}





rg_dList<MBSEdge*>*  MBSShell::getEdges() 
{
    return &m_edges;
}



rg_dList<MBSVertex*>* MBSShell::getVertices() 
{
    return &m_vertices;
}




void MBSShell::setBody(MBSBody* o_body)
{
    m_body = o_body;
}



void MBSShell::addFace(MBSFace* o_face)
{
    m_faces.add( o_face );
}





void MBSShell::addEdge(MBSEdge* o_edge)
{
    m_edges.add( o_edge );
}



void MBSShell::addVertex(MBSVertex* o_vertex)
{
    m_vertices.add( o_vertex );
}




MBSShell& MBSShell::operator =(const MBSShell& shell)
{
    if ( this == &shell ) {
        return *this;
    }

    TopologicalEntity::operator =(shell);

    m_body = shell.m_body;

    m_faces.duplicateList(shell.m_faces);
    m_edges.duplicateList(shell.m_edges);
    m_vertices.duplicateList(shell.m_vertices);

    m_visited = shell.m_visited;

    return *this;
}

