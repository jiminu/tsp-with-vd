#include "SurfaceModel.h"


#include <set>
#include <vector>
using namespace std;




SurfaceModel::SurfaceModel()
: m_number_of_bodies(0),
  m_number_of_shells(0),
  m_number_of_faces(0),
  m_number_of_loops(0),
  m_number_of_edges(0),
  m_number_of_vertices(0)
{
}



SurfaceModel::SurfaceModel(const SurfaceModel& surfaceModel)
{
}



SurfaceModel::~SurfaceModel()
{
    clear();
}



void SurfaceModel::clear()
{
    for (list<Body*>::iterator i_body=m_bodies.begin(); i_body!=m_bodies.end(); ++i_body ) {
        delete (*i_body);
    }

    for (list<Shell*>::iterator i_shell=m_shells.begin(); i_shell!=m_shells.end(); ++i_shell ) {
        delete (*i_shell);
    }

    for (list<Face*>::iterator i_face=m_faces.begin(); i_face!=m_faces.end(); ++i_face ) {
        delete (*i_face);
    }

    for (list<Loop*>::iterator i_loop=m_loops.begin(); i_loop!=m_loops.end(); ++i_loop ) {
        delete (*i_loop);
    }

    for (list<Edge*>::iterator i_edge=m_edges.begin(); i_edge!=m_edges.end(); ++i_edge ) {
        delete (*i_edge);
    }

    for (list<Vertex*>::iterator i_vx=m_vertices.begin(); i_vx!=m_vertices.end(); ++i_vx ) {
        delete (*i_vx);
    }


    m_bodies.clear();
    m_shells.clear();
    m_faces.clear();
    m_loops.clear();
    m_edges.clear();
    m_vertices.clear();

    m_number_of_bodies   = 0;
    m_number_of_shells   = 0;
    m_number_of_faces    = 0;
    m_number_of_loops    = 0;
    m_number_of_edges    = 0;
    m_number_of_vertices = 0;
}

    


void SurfaceModel::release_deleted_entities()
{
    list<Vertex*>::iterator i_vtx=m_vertices.begin();
    while ( i_vtx!=m_vertices.end() ) {
        if ( (*i_vtx)->will_be_deleted() ) {
            delete (*i_vtx);
            --m_number_of_vertices;
            i_vtx = m_vertices.erase( i_vtx );
        }
        else {
            ++i_vtx;
        }
    }

    list<Edge*>::iterator i_edge=m_edges.begin();
    while ( i_edge!=m_edges.end() ) {
        if ( (*i_edge)->will_be_deleted() ) {
            delete (*i_edge);
            --m_number_of_edges;
            i_edge = m_edges.erase( i_edge );
        }
        else {
            ++i_edge;
        }
    }

    list<Loop*>::iterator i_loop=m_loops.begin();
    while ( i_loop!=m_loops.end() ) {
        if ( (*i_loop)->will_be_deleted() ) {
            delete (*i_loop);
            --m_number_of_loops;
            i_loop = m_loops.erase( i_loop );
        }
        else {
            ++i_loop;
        }
    }


    list<Face*>::iterator i_face=m_faces.begin();
    while ( i_face!=m_faces.end() ) {
        if ( (*i_face)->will_be_deleted() ) {
            delete (*i_face);
            --m_number_of_faces;
            i_face = m_faces.erase( i_face );
        }
        else {
            ++i_face;
        }
    }

    list<Shell*>::iterator i_shell=m_shells.begin(); 
    while ( i_shell!=m_shells.end() ) {
        if ( (*i_shell)->will_be_deleted() ) {
            delete (*i_shell);
            --m_number_of_shells;
            i_shell = m_shells.erase( i_shell );
        }
        else {
            ++i_shell;
        }
    }

    list<Body*>::iterator i_body=m_bodies.begin();
    while ( i_body!=m_bodies.end() ) {
        if ( (*i_body)->will_be_deleted() ) {
            delete (*i_body);
            --m_number_of_bodies;
            i_body = m_bodies.erase( i_body );
        }
        else {
            ++i_body;
        }
    }
}



SurfaceModel& SurfaceModel::operator =(const SurfaceModel& surfaceModel)
{
    if ( this != &surfaceModel ) {
    }

    return *this;
}



const list< SurfaceModel::Body* >&    SurfaceModel::get_all_bodies() const
{
    return m_bodies;
}



const list< SurfaceModel::Shell* >&   SurfaceModel::get_all_shells() const
{
    return m_shells;
}



const list< SurfaceModel::Face* >&    SurfaceModel::get_all_faces() const
{
    return m_faces;
}



const list< SurfaceModel::Loop* >&    SurfaceModel::get_all_loops() const
{
    return m_loops;
}



const list< SurfaceModel::Edge* >&    SurfaceModel::get_all_edges() const
{
    return m_edges;
}



const list< SurfaceModel::Vertex* >&  SurfaceModel::get_all_vertices() const
{
    return m_vertices;
}



unsigned int            SurfaceModel::number_of_bodies() const
{
    return m_number_of_bodies;
}



unsigned int            SurfaceModel::number_of_shells() const
{
    return m_number_of_shells;
}



unsigned int            SurfaceModel::number_of_faces() const
{
    return m_number_of_faces;
}



unsigned int            SurfaceModel::number_of_loops() const
{
    return m_number_of_loops;
}



unsigned int            SurfaceModel::number_of_edges() const
{
    return m_number_of_edges;
}



unsigned int            SurfaceModel::number_of_vertices() const
{
    return m_number_of_vertices;
}



SurfaceModel::Body*     SurfaceModel::create_body()
{
    int ID = 0;
    if (!m_bodies.empty()) {
        ID = m_bodies.back()->getID() + 1;
    }
    m_bodies.push_back( new Body(ID) );
    ++m_number_of_bodies;

    return m_bodies.back();
}



SurfaceModel::Shell*    SurfaceModel::create_shell()
{
    int ID = 0;
    if (!m_shells.empty()) {
        ID = m_shells.back()->getID() + 1;
    }
    m_shells.push_back( new Shell(ID) );
    ++m_number_of_shells;

    return m_shells.back();
}



SurfaceModel::Face*     SurfaceModel::create_face()
{
    int ID = 0;
    if (!m_faces.empty()) {
        ID = m_faces.back()->getID() + 1;
    }
    m_faces.push_back( new Face(ID) );
    ++m_number_of_faces;

    return m_faces.back();
}



SurfaceModel::Loop*     SurfaceModel::create_loop()
{
    int ID = 0;
    if (!m_loops.empty()) {
        ID = m_loops.back()->getID() + 1;
    }
    m_loops.push_back( new Loop(ID) );
    ++m_number_of_loops;

    return m_loops.back();
}



SurfaceModel::Edge*     SurfaceModel::create_edge()
{
    int ID = 0;
    if (!m_edges.empty()) {
        ID = m_edges.back()->getID() + 1;
    }
    m_edges.push_back( new Edge(ID) );
    ++m_number_of_edges;

    return m_edges.back();
}



SurfaceModel::Vertex*   SurfaceModel::create_vertex()
{
    int ID = 0;
    if (!m_vertices.empty()) {
        ID = m_vertices.back()->getID() + 1;
    }
    m_vertices.push_back( new Vertex(ID) );
    ++m_number_of_vertices;

    return m_vertices.back();
}



void                    SurfaceModel::concatenate_body(   const SurfaceModel::Body*       const body )
{
    m_bodies.push_back( const_cast<Body*>(body) ); 
    ++m_number_of_bodies;
}



void                    SurfaceModel::concatenate_shell(  const SurfaceModel::Shell*      const shell )
{
    m_shells.push_back( const_cast<Shell*>(shell) ); 
    ++m_number_of_shells;
}



void                    SurfaceModel::concatenate_face(   const SurfaceModel::Face*       const face )
{
    m_faces.push_back( const_cast<Face*>(face) ); 
    ++m_number_of_faces;
}



void                    SurfaceModel::concatenate_loop(   const SurfaceModel::Loop*       const loop )
{
    m_loops.push_back( const_cast<Loop*>(loop) ); 
    ++m_number_of_loops;
}



void                    SurfaceModel::concatenate_edge(   const SurfaceModel::Edge*       const edge )
{
    m_edges.push_back( const_cast<Edge*>(edge) ); 
    ++m_number_of_edges;
}



void                    SurfaceModel::concatenate_vertex( const SurfaceModel::Vertex*     const vertex )
{
    m_vertices.push_back( const_cast<Vertex*>(vertex) ); 
    ++m_number_of_vertices;
}



void                    SurfaceModel::remove_body(        const SurfaceModel::Body*       const body )
{
    if ( body != rg_NULL ) {
        m_bodies.remove( const_cast<Body*>(body) );
        --m_number_of_bodies;
        delete body;
    }
}



void                    SurfaceModel::remove_shell(       const SurfaceModel::Shell*      const shell )
{
    if ( shell != rg_NULL ) {
        m_shells.remove( const_cast<Shell*>( shell ) );
        --m_number_of_shells;
        delete shell;
    }
}



void                    SurfaceModel::remove_face(        const SurfaceModel::Face*       const face )
{
    if ( face != rg_NULL ) {
        m_faces.remove( const_cast<Face*>( face ) );
        --m_number_of_faces;
        delete face;
    }
}



void                    SurfaceModel::remove_loop(        const SurfaceModel::Loop*       const loop )
{
    if ( loop != rg_NULL ) {
        m_loops.remove( const_cast<Loop*>( loop ) );
        --m_number_of_loops;
        delete loop;
    }
}



void                    SurfaceModel::remove_edge(        const SurfaceModel::Edge*       const edge )
{
    if ( edge != rg_NULL ) {
        m_edges.remove( const_cast<Edge*>( edge ) );
        --m_number_of_edges;
        delete edge;
    }
}



void                    SurfaceModel::remove_vertex(      const SurfaceModel::Vertex*     const vertex )
{
    if ( vertex != rg_NULL ) {
        m_vertices.remove( const_cast<Vertex*>( vertex ) );
        --m_number_of_vertices;
        delete vertex;
    }
}

    

void  SurfaceModel::splice(SurfaceModel::Shell* shellToExtend, SurfaceModel::Shell* shellToRemove)
{
    if ( shellToExtend == shellToRemove ) {
        return;
    }


    list< SurfaceModel::Face* > facesToBeMoved;
    shellToRemove->find_bounding_faces( facesToBeMoved );

    list< SurfaceModel::Face* >::iterator i_face;
    for ( i_face=facesToBeMoved.begin(); i_face!=facesToBeMoved.end(); ++i_face ) {
        SurfaceModel::Face* currFace = *i_face;

        currFace->set_shell( shellToExtend );
        shellToExtend->add_face( currFace );
    }


    SurfaceModel::Body*  bodyToRemove = shellToRemove->body();
    remove_shell( shellToRemove );
    remove_body(  bodyToRemove );
}



bool SurfaceModel::divide_edge_by_inserting_vertex(SurfaceModel::Edge* edge, SurfaceModel::Vertex* splitVertex, SurfaceModel::Edge*& outputEdge1, SurfaceModel::Edge*& outputEdge2)
{
    int edgeID = m_edges.back()->getID() + 1;

    Vertex* start = edge->start_vertex();
    Vertex* end   = edge->end_vertex();

    if ( start != rg_NULL && end != rg_NULL ) {
        //  newEdges[0] = this : vertex[0] -> vertex[1]
        //  newEdges[1]        : vertex[1] -> vertex[2]
        Vertex* vertex[3]  = {start, splitVertex, end};
        Edge*   newEdge[2] = {rg_NULL, rg_NULL};

        newEdge[0] = edge;
        newEdge[0]->set_end_vertex( vertex[1] );

        newEdge[1] = create_edge();
        newEdge[1]->setID( edgeID );
        newEdge[1]->set_start_vertex( vertex[1] );
        newEdge[1]->set_end_vertex( vertex[2] );

        newEdge[1]->set_right_face( edge->right_face() );
        newEdge[1]->set_left_face(  edge->left_face() );
        
        newEdge[1]->connect_left_hand_edge(  edge->left_hand_edge() );
        newEdge[1]->connect_right_hand_edge( edge->right_hand_edge() );
        newEdge[1]->connect_left_leg_edge(   edge );
        newEdge[1]->connect_right_leg_edge(  edge );


        vertex[1]->set_first_edge( newEdge[1] );
        vertex[2]->set_first_edge( newEdge[1] );


        outputEdge1 = newEdge[0];
        outputEdge2 = newEdge[1];


        return true;
    }
    else {
        return false;
    }
}

    

bool                    SurfaceModel::divide_face_into_triangular_faces_by_inserting_vertex( 
                                        SurfaceModel::Face* face, SurfaceModel::Vertex* vertex, 
                                        list<SurfaceModel::Edge*>& newEdges, list<SurfaceModel::Face*>& newFaces)
{
    if ( !face->no_interior_loop() ) {
        return false;
    }


    int i=0;
    int edgeID = m_edges.back()->getID();
    int faceID = m_faces.back()->getID();
    int numBoundingVertices = face->number_of_bounding_vertices();


    //  collect bounding edges and vertices of face
    vector<SurfaceModel::Vertex*> boundingVertex(numBoundingVertices+2, rg_NULL);
    vector<SurfaceModel::Edge*>   boundingEdge(  numBoundingVertices+2, rg_NULL);
    {
        int pos = 1;
        SurfaceModel::Edge* firstEdge = face->first_edge();
        SurfaceModel::Edge* currEdge  = firstEdge;

        do {
            boundingEdge[pos]   = currEdge;
            if ( face == currEdge->left_face() ) {
                boundingVertex[pos] = currEdge->start_vertex();
                currEdge = currEdge->left_hand_edge();
            }
            else {
                boundingVertex[pos] = currEdge->end_vertex();
                currEdge = currEdge->right_leg_edge();
            }
            ++pos;
        } while ( currEdge != firstEdge );
        boundingVertex.front() = boundingVertex[numBoundingVertices];
        boundingVertex.back()  = boundingVertex[1];
        boundingEdge.front()   = boundingEdge[numBoundingVertices];
        boundingEdge.back()    = boundingEdge[1];
    }


    //  make new edges and faces
    vector<SurfaceModel::Edge*> newEdge(numBoundingVertices+2, rg_NULL);
    vector<SurfaceModel::Face*> newFace(numBoundingVertices+2, rg_NULL);
    {
        for ( i=1; i<=numBoundingVertices; ++i ) {
            newEdge[i] = create_edge();
            newEdge[i]->setID( ++edgeID );
            newEdge[i]->set_start_vertex( vertex );
            newEdge[i]->set_end_vertex(   boundingVertex[i] );
        }
        newEdge.front() = newEdge[numBoundingVertices];
        newEdge.back()  = newEdge[1];


        newFace[1] = face;
        for ( i=2; i<=numBoundingVertices; ++i ) {
            newFace[i] = create_face();
            newFace[i]->setID( ++faceID );

            SurfaceModel::Loop* newLoop = create_loop();
            newLoop->setID( faceID );

            newFace[i]->set_exterior_loop( newLoop );
            newLoop->set_face( newFace[i] );

            newFace[i]->set_shell( face->shell() );
            face->shell()->add_face( newFace[i] );
        }
        newFace.front() = newFace[numBoundingVertices];
        newFace.back()  = newFace[1];
    }


    //  set topology among new entities and from new ones to bounding ones.
    for ( i=1; i<=numBoundingVertices; ++i ) {
        newEdge[i]->set_left_face(       newFace[i]        );
        newEdge[i]->set_right_face(      newFace[i-1]      );

        newEdge[i]->set_left_hand_edge(  boundingEdge[i]   );
        newEdge[i]->set_right_hand_edge( boundingEdge[i-1] );
        newEdge[i]->set_left_leg_edge(   newEdge[i+1]      );
        newEdge[i]->set_right_leg_edge(  newEdge[i-1]      );

        newFace[i]->set_first_edge(                  boundingEdge[i] );
        newFace[i]->exterior_loop()->set_first_edge( boundingEdge[i] );
    }


    //  set topology from bounding edge to new entities
    for ( i=1; i<=numBoundingVertices; ++i ) {
        if ( face == boundingEdge[i]->left_face() ) {
            boundingEdge[i]->set_left_face(      newFace[i]   );
            boundingEdge[i]->set_left_hand_edge( newEdge[i+1] );
            boundingEdge[i]->set_left_leg_edge(  newEdge[i]   );
        }
        else {
            boundingEdge[i]->set_right_face(      newFace[i]   );
            boundingEdge[i]->set_right_hand_edge( newEdge[i] );
            boundingEdge[i]->set_right_leg_edge(  newEdge[i+1]   );
        }
    }  

    //  set first edge of inserted vertex
    vertex->set_first_edge( newEdge[1] );



    //  make output arguments.
    for ( i=1; i<=numBoundingVertices; ++i ) {
        newEdges.push_back( newEdge[i] );
        newFaces.push_back( newFace[i] );
    }

    return true;
}



bool                    SurfaceModel::divide_face_by_inserting_edge( 
                                        SurfaceModel::Face* face, SurfaceModel::Vertex* vtx1OfEdge, SurfaceModel::Vertex* vtx2OfEdge, 
                                        SurfaceModel::Edge* &newEdge, SurfaceModel::Face* &newFace1, SurfaceModel::Face* &newFace2)
{
    //      splitFace
    //              vtx2OnSplitFace                 
    //      *---------------*---------------*       
    //      |               |               |       
    //      |    newFace1   |   newFace2    |        
    //      |               |               |        
    //      |               |               |       
    //      *---------------*---------------*       
    //              vtx1OnSplitFace                         
    //
    int faceID = m_faces.back()->getID()+1;
    int loopID = faceID;
    int edgeID = m_edges.back()->getID()+1;
    
    if ( face->no_interior_loop() ) {
        //  create new entities and rearrange them 
        Face*   newFace[2] = { face, create_face() };
        newFace[1]->setID( faceID );

        Loop*   newLoop[2] = { face->exterior_loop(), create_loop() };
        newLoop[1]->setID( loopID );

        newEdge = create_edge();
        newEdge->setID( edgeID );

        Vertex* vtxOfNewEdge[2] = { vtx1OfEdge, vtx2OfEdge };

        Edge* edgeIncidentToVtx[2][2] = { {rg_NULL, rg_NULL}, {rg_NULL, rg_NULL} };
        face->find_bounding_edges_incident_to_vertex( vtxOfNewEdge[0], edgeIncidentToVtx[0][0], edgeIncidentToVtx[0][1]);
        face->find_bounding_edges_incident_to_vertex( vtxOfNewEdge[1], edgeIncidentToVtx[1][0], edgeIncidentToVtx[1][1]);


        //  set topology between vertices and edges;
        newEdge->set_start_vertex( vtxOfNewEdge[0] );
        newEdge->set_end_vertex(   vtxOfNewEdge[1] );


        //  set topology among edges;
        newEdge->connect_left_hand_edge(  edgeIncidentToVtx[1][1] );
        newEdge->connect_right_hand_edge( edgeIncidentToVtx[1][0] );
        newEdge->connect_left_leg_edge(   edgeIncidentToVtx[0][0] );
        newEdge->connect_right_leg_edge(  edgeIncidentToVtx[0][1] );


        //  set topology between edges and faces;
        newEdge->set_left_face(  newFace[0] );
        newEdge->set_right_face( newFace[1] );

        Edge* currEdge = newEdge;
        do {
            Edge* nextEdge = rg_NULL;
            if ( newFace[1] == currEdge->right_face() ) {
                nextEdge = currEdge->right_hand_edge();
                if ( currEdge->end_vertex() == nextEdge->start_vertex() ) {
                    nextEdge->set_right_face( newFace[1] );
                }
                else {
                    nextEdge->set_left_face(  newFace[1] );
                }
            }
            else {
                nextEdge = currEdge->left_leg_edge();
                if ( currEdge->start_vertex() == nextEdge->start_vertex() ) {
                    nextEdge->set_right_face( newFace[1] );
                }
                else {
                    nextEdge->set_left_face(  newFace[1] );
                }
            }

            currEdge = nextEdge;
        } while ( currEdge != newEdge );


        //  set topology between edges and loops;
        newLoop[0]->set_first_edge( newEdge );
        newLoop[1]->set_first_edge( newEdge );
        newLoop[0]->set_face( newFace[0] );
        newLoop[1]->set_face( newFace[1] );



        //  set topology between loops and faces;
        newFace[0]->set_first_edge( newEdge );
        newFace[1]->set_first_edge( newEdge );
        newFace[0]->set_exterior_loop( newLoop[0] );
        newFace[1]->set_exterior_loop( newLoop[1] );

        //  set topology between faces and shell;
        newFace[1]->set_shell( face->shell() );
        face->shell()->add_face( newFace[1] );

        newFace1 = newFace[0];
        newFace2 = newFace[1];

        return true;
    }
    else {
        return false;
    }
}


    

bool SurfaceModel::merge_two_faces_by_removing_edge( SurfaceModel::Face* face1, SurfaceModel::Face* face2 )
{
    list<SurfaceModel::Edge*> sharingEdges;
    int numSharingEdges = face1->find_edges_sharing_with_adjacent_face( face2, sharingEdges );
    if ( numSharingEdges == 0 ) {
        return false;
    }


    //  check that two faces intersect in a single edge-chain.
    map<SurfaceModel::Vertex*, vector<SurfaceModel::Edge*> > vertex2edgeMap;
    map<SurfaceModel::Vertex*, vector<SurfaceModel::Edge*> >::iterator it_vtx2edge;
    list<SurfaceModel::Edge*>::iterator i_edge;
    for ( i_edge=sharingEdges.begin(); i_edge!=sharingEdges.end(); ++i_edge ) {
        Edge*   currEdge = *i_edge;
        Vertex* vertex[2] = { currEdge->start_vertex(), currEdge->end_vertex() };

        for ( int i=0; i<2; ++i ) {
            it_vtx2edge = vertex2edgeMap.find( vertex[i] );
            if ( it_vtx2edge == vertex2edgeMap.end() ) {
                vertex2edgeMap.insert( make_pair( vertex[i], vector<SurfaceModel::Edge*>() ) );
                it_vtx2edge = vertex2edgeMap.find( vertex[i] );
            }
            it_vtx2edge->second.push_back( currEdge );
        }
    }

    int numVerticesWithDegree1 = 0;
    for ( it_vtx2edge=vertex2edgeMap.begin(); it_vtx2edge!=vertex2edgeMap.end(); ++it_vtx2edge) {
        if ( it_vtx2edge->second.size() == 1 ) {
            ++numVerticesWithDegree1;
        }
    }

    if ( numVerticesWithDegree1 != 2 ) {
        return false;
    }

    bool    doTwoFacesMeetInAnotherPositionExceptSharingEdges = false;
    list<SurfaceModel::Vertex*> boundingVerticesOfFace2;
    face2->find_bounding_vertices( boundingVerticesOfFace2 );
    for ( list<SurfaceModel::Vertex*>::iterator it_v=boundingVerticesOfFace2.begin(); it_v!=boundingVerticesOfFace2.end(); ++it_v ) {
        SurfaceModel::Vertex* currVtxOnFace2 = *it_v;
        if ( vertex2edgeMap.find( currVtxOnFace2 ) != vertex2edgeMap.end() ) {
            continue;
        }

        if ( face1->is_incident_to( currVtxOnFace2 ) ) {
            doTwoFacesMeetInAnotherPositionExceptSharingEdges = true;
            break;
        }
    }
    if ( doTwoFacesMeetInAnotherPositionExceptSharingEdges ) {
        return false;
    }



    /////////////////////////////////////////////////////////////
    //
    //  merge two faces which share just a single edge-chain.
    //
    list<SurfaceModel::Edge*> boundingEdgesOfFace2;
    face2->find_bounding_edges( boundingEdgesOfFace2 );
    for ( i_edge=boundingEdgesOfFace2.begin(); i_edge!=boundingEdgesOfFace2.end(); ++i_edge ) {
        Edge*   currEdge = *i_edge;
        if ( currEdge->right_face() == face2 ) {
            currEdge->set_right_face( face1 );
        }
        else {
            currEdge->set_left_face( face1 );
        }
    }
    face2->shell()->remove_face( face2 );



    for ( it_vtx2edge=vertex2edgeMap.begin(); it_vtx2edge!=vertex2edgeMap.end(); ++it_vtx2edge) {
        if ( it_vtx2edge->second.size() != 1 ) {
            continue;
        }

        Vertex* currVtx     = it_vtx2edge->first;
        Edge*   sharingEdge = it_vtx2edge->second[0];

        Edge*   incidentEdge[2] = {rg_NULL, rg_NULL};
        if ( currVtx->is_end_vertex_of( sharingEdge ) ) {
            incidentEdge[0] = sharingEdge->right_hand_edge();
            incidentEdge[1] = sharingEdge->left_hand_edge();
        }
        else {
            incidentEdge[0] = sharingEdge->left_leg_edge();
            incidentEdge[1] = sharingEdge->right_leg_edge();
        }

        if ( currVtx->is_end_vertex_of( incidentEdge[0] ) ) {
            if ( face1 == incidentEdge[0]->left_face() ) {
                incidentEdge[0]->connect_left_hand_edge( incidentEdge[1] );
            }
            else {
                incidentEdge[0]->connect_right_hand_edge( incidentEdge[1] );
            }
        }
        else {
            if ( face1 == incidentEdge[0]->left_face() ) {
                incidentEdge[0]->connect_left_leg_edge( incidentEdge[1] );
            }
            else {
                incidentEdge[0]->connect_right_leg_edge( incidentEdge[1] );
            }
        }

        currVtx->set_first_edge( incidentEdge[0] );
        face1->set_first_edge( incidentEdge[0] );
        face1->exterior_loop()->set_first_edge( incidentEdge[0] );
    }



    for ( i_edge=sharingEdges.begin(); i_edge!=sharingEdges.end(); ++i_edge ) {
        (*i_edge)->will_be_deleted( true );
    }

    for ( it_vtx2edge=vertex2edgeMap.begin(); it_vtx2edge!=vertex2edgeMap.end(); ++it_vtx2edge) {
        if ( it_vtx2edge->second.size() != 1 ) {
            it_vtx2edge->first->will_be_deleted( true );
        }
    }

    face2->will_be_deleted( true );
    face2->exterior_loop()->will_be_deleted( true );

    return true;
}



void SurfaceModel::duplicate(const SurfaceModel& surfaceModel)
{
    clear();

    map<Body*,   Body*>   bodyMap;
    map<Shell*,  Shell*>  shellMap;
    map<Face*,   Face*>   faceMap;
    map<Loop*,   Loop*>   loopMap;
    map<Edge*,   Edge*>   edgeMap;
    map<Vertex*, Vertex*> vertexMap;
    make_entity_map_for_duplication(    surfaceModel, bodyMap, shellMap, faceMap, loopMap, edgeMap, vertexMap );

    duplicate_topology(                 surfaceModel, bodyMap, shellMap, faceMap, loopMap, edgeMap, vertexMap );
}


void            SurfaceModel::make_entity_map_for_duplication( const SurfaceModel&          surfaceModel, 
                                                               map<SurfaceModel::Body*,     SurfaceModel::Body*>&    bodyMap,
                                                               map<SurfaceModel::Shell*,    SurfaceModel::Shell*>&   shellMap,
                                                               map<SurfaceModel::Face*,     SurfaceModel::Face*>&    faceMap,
                                                               map<SurfaceModel::Loop*,     SurfaceModel::Loop*>&    loopMap,
                                                               map<SurfaceModel::Edge*,     SurfaceModel::Edge*>&    edgeMap,
                                                               map<SurfaceModel::Vertex*,   SurfaceModel::Vertex*>&  vertexMap)
{
    //  make map for bodies
    list<Body*>::const_iterator i_body;
    for ( i_body=surfaceModel.m_bodies.begin(); i_body!=surfaceModel.m_bodies.end(); ++i_body ) {
        Body* originalBody   = const_cast<Body*>( *i_body );
        Body* duplicatedBody = this->create_body();
        duplicatedBody->setID( originalBody->getID() );

        bodyMap.insert( make_pair( originalBody, duplicatedBody ) );
    }
    bodyMap.insert( make_pair( (Body*)rg_NULL, (Body*)rg_NULL ) );


    // make map for shells
    list<Shell*>::const_iterator i_shell;
    for ( i_shell=surfaceModel.m_shells.begin(); i_shell!=surfaceModel.m_shells.end(); ++i_shell ) {
        Shell* originalShell   = const_cast<Shell*>( *i_shell );
        Shell* duplicatedShell = this->create_shell();
        duplicatedShell->setID( originalShell->getID() );

        shellMap.insert( make_pair( originalShell, duplicatedShell ) );
    }
    shellMap.insert( make_pair( (Shell*)rg_NULL, (Shell*)rg_NULL ) );


    // make map for faces
    list<Face*>::const_iterator i_face;
    for ( i_face=surfaceModel.m_faces.begin(); i_face!=surfaceModel.m_faces.end(); ++i_face ) {
        Face* originalFace   = const_cast<Face*>( *i_face );
        Face* duplicatedFace = this->create_face();
        duplicatedFace->setID( originalFace->getID() );

        faceMap.insert( make_pair( originalFace, duplicatedFace ) );
    }
    faceMap.insert( make_pair( (Face*)rg_NULL, (Face*)rg_NULL ) );


    //  make map for loops
    list<Loop*>::const_iterator i_loop;
    for ( i_loop=surfaceModel.m_loops.begin(); i_loop!=surfaceModel.m_loops.end(); ++i_loop ) {
        Loop* originalLoop = const_cast<Loop*>( *i_loop );
        Loop* duplicatedLoop = this->create_loop();
        duplicatedLoop->setID( originalLoop->getID() );

        loopMap.insert( make_pair( originalLoop, duplicatedLoop ) );
    }
    loopMap.insert( make_pair( (Loop*)rg_NULL, (Loop*)rg_NULL ) );


    // make map for edges
    list<Edge*>::const_iterator i_edge;
    for ( i_edge=surfaceModel.m_edges.begin(); i_edge!=surfaceModel.m_edges.end(); ++i_edge ) {
        Edge* originalEdge = const_cast<Edge*>( *i_edge );
        Edge* duplicatedEdge = this->create_edge();
        duplicatedEdge->setID( originalEdge->getID() );

        edgeMap.insert( make_pair( originalEdge, duplicatedEdge ) );
    }
    edgeMap.insert( make_pair( (Edge*)rg_NULL, (Edge*)rg_NULL ) );


    // make map for vertices
    list<Vertex*>::const_iterator i_vtx;
    for ( i_vtx=surfaceModel.m_vertices.begin(); i_vtx!=surfaceModel.m_vertices.end(); ++i_vtx ) {
        Vertex* originalVtx = const_cast<Vertex*>( *i_vtx );
        Vertex* duplicatedVtx = this->create_vertex();
        duplicatedVtx->setID( originalVtx->getID() );

        vertexMap.insert( make_pair( originalVtx, duplicatedVtx ) );
    }
    vertexMap.insert( make_pair( (Vertex*)rg_NULL, (Vertex*)rg_NULL ) );

}



void            SurfaceModel::duplicate_topology(        const SurfaceModel&                surfaceModel, 
                                                         const map<SurfaceModel::Body*,     SurfaceModel::Body*>&    bodyMap,
                                                         const map<SurfaceModel::Shell*,    SurfaceModel::Shell*>&   shellMap,
                                                         const map<SurfaceModel::Face*,     SurfaceModel::Face*>&    faceMap,
                                                         const map<SurfaceModel::Loop*,     SurfaceModel::Loop*>&    loopMap,
                                                         const map<SurfaceModel::Edge*,     SurfaceModel::Edge*>&    edgeMap,
                                                         const map<SurfaceModel::Vertex*,   SurfaceModel::Vertex*>&  vertexMap)
{
    list<Body*>::const_iterator i_body;
    for ( i_body=surfaceModel.m_bodies.begin(); i_body!=surfaceModel.m_bodies.end(); ++i_body ) {
        Body* original   = const_cast<Body*>( *i_body );
        Body* duplicated = bodyMap.find( original )->second;

        duplicated->set_exterior_shell(  shellMap.find( original->exterior_shell() )->second  );

        const list<Shell*>& interiorShellOfOriginal = original->interior_shell();
        for (list<Shell*>::const_iterator i_shell=interiorShellOfOriginal.begin(); i_shell!=interiorShellOfOriginal.end(); ++i_shell ) {
            duplicated->add_interior_shell(  shellMap.find( *i_shell )->second  );
        }
    }


    list<Shell*>::const_iterator i_shell;
    for ( i_shell=surfaceModel.m_shells.begin(); i_shell!=surfaceModel.m_shells.end(); ++i_shell ) {
        Shell* original   = const_cast<Shell*>( *i_shell );
        Shell* duplicated = shellMap.find( original )->second;

        duplicated->set_body( bodyMap.find( original->body() )->second );

        const list<Face* >& boundingFacesOfOriginal = original->faces();
        for ( list<Face*>::const_iterator i_face=boundingFacesOfOriginal.begin(); i_face!=boundingFacesOfOriginal.end(); ++i_face ) {
            duplicated->add_face( faceMap.find( *i_face )->second );
        }
    }


    list<Face*>::const_iterator i_face;
    for ( i_face=surfaceModel.m_faces.begin(); i_face!=surfaceModel.m_faces.end(); ++i_face ) {
        Face* original   = const_cast<Face*>( *i_face );
        Face* duplicated = faceMap.find( original )->second;

        duplicated->set_first_edge(    edgeMap.find(  original->first_edge()    )->second );
        duplicated->set_shell(         shellMap.find( original->shell()         )->second );
        duplicated->set_exterior_loop( loopMap.find(  original->exterior_loop() )->second );

        const list<Loop*>& interiorLoopsOfOriginal = original->interior_loops();
        for ( list<Loop*>::const_iterator i_loop=interiorLoopsOfOriginal.begin(); i_loop!=interiorLoopsOfOriginal.end(); ++i_loop ) {
            duplicated->add_interior_loop( loopMap.find( *i_loop )->second );
        }
    }


    list<Loop*>::const_iterator i_loop;
    for ( i_loop=surfaceModel.m_loops.begin(); i_loop!=surfaceModel.m_loops.end(); ++i_loop ) {
        Loop* original   = const_cast<Loop*>( *i_loop );
        Loop* duplicated = loopMap.find( original )->second;

        duplicated->set_face(       faceMap.find( original->face() )->second );
        duplicated->set_first_edge( edgeMap.find( original->first_edge() )->second );
    }


    list<Edge*>::const_iterator i_edge;
    for ( i_edge=surfaceModel.m_edges.begin(); i_edge!=surfaceModel.m_edges.end(); ++i_edge ) {
        Edge* original   = const_cast<Edge*>( *i_edge );
        Edge* duplicated = edgeMap.find( original )->second;

        duplicated->set_edge( vertexMap.find( original->start_vertex()    )->second,
                              vertexMap.find( original->end_vertex()      )->second,
                              faceMap.find(   original->right_face()      )->second,
                              faceMap.find(   original->left_face()       )->second,
                              edgeMap.find(   original->right_hand_edge() )->second,
                              edgeMap.find(   original->left_hand_edge()  )->second,
                              edgeMap.find(   original->right_leg_edge()  )->second,
                              edgeMap.find(   original->left_leg_edge()   )->second );
    }


    list<Vertex*>::const_iterator i_vtx;
    for ( i_vtx=surfaceModel.m_vertices.begin(); i_vtx!=surfaceModel.m_vertices.end(); ++i_vtx ) {
        Vertex* original   = const_cast<Vertex*>( *i_vtx );
        Vertex* duplicated = vertexMap.find( original )->second;

        duplicated->set_first_edge( edgeMap.find( original->first_edge() )->second );
    }
}





///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Body 
//
SurfaceModel::Body::Body()
: TopologicalEntity(),
  m_exteriorShell(rg_NULL)
{
}



SurfaceModel::Body::Body(const rg_INT& ID)
: TopologicalEntity(ID),
  m_exteriorShell(rg_NULL)
{
}



SurfaceModel::Body::Body(const Body& body)
: TopologicalEntity( body ),
  m_exteriorShell(   body.m_exteriorShell ),
  m_interiorShell(   body.m_interiorShell )
{
}



SurfaceModel::Body::~Body()
{
}



SurfaceModel::Body& SurfaceModel::Body::operator =(const Body& body)
{
    if ( this != &body ) {
        TopologicalEntity::operator=(body);
        m_exteriorShell = body.m_exteriorShell;
        m_interiorShell = body.m_interiorShell;
    }

    return *this;
}



SurfaceModel::Shell*                 SurfaceModel::Body::exterior_shell() const
{
    return m_exteriorShell;
}



const list< SurfaceModel::Shell* >&  SurfaceModel::Body::interior_shell() const
{
    return m_interiorShell;
}



void                    SurfaceModel::Body::set_exterior_shell(        const SurfaceModel::Shell* const shell)
{
    m_exteriorShell = const_cast<SurfaceModel::Shell*>( shell );
}



void                    SurfaceModel::Body::add_interior_shell(    const SurfaceModel::Shell* const shell)
{
    m_interiorShell.push_back( const_cast<SurfaceModel::Shell*>( shell ) );
}



void                    SurfaceModel::Body::remove_interior_shell( const SurfaceModel::Shell* const shell)
{
    m_interiorShell.remove( const_cast<SurfaceModel::Shell*>( shell ) );
}



rg_INT                  SurfaceModel::Body::find_bounding_shells(         list<SurfaceModel::Shell*>&  boundingShells) const
{
    boundingShells.push_back( m_exteriorShell );
    boundingShells.insert( boundingShells.end(), m_interiorShell.begin(), m_interiorShell.end() );

    return boundingShells.size();
}



rg_INT                  SurfaceModel::Body::find_bounding_faces(          list<SurfaceModel::Face*>&   boundingFaces) const
{
    list<SurfaceModel::Shell*> boundingShells;
    find_bounding_shells( boundingShells );

    rg_INT numFaces = 0;
    list<SurfaceModel::Shell*>::iterator i_shell;
    for ( i_shell=boundingShells.begin(); i_shell!=boundingShells.end(); ++i_shell ) {
        numFaces += (*i_shell)->find_bounding_faces( boundingFaces );
    }

    return numFaces;
}



rg_INT                  SurfaceModel::Body::find_loops_of_bounding_faces( list<SurfaceModel::Loop*>&   loopsOfBoundingFaces) const
{
    list<SurfaceModel::Shell*> boundingShells;
    find_bounding_shells( boundingShells );

    rg_INT numLoops = 0;
    list<SurfaceModel::Shell*>::iterator i_shell;
    for ( i_shell=boundingShells.begin(); i_shell!=boundingShells.end(); ++i_shell ) {
        numLoops += (*i_shell)->find_loops_of_bounding_faces( loopsOfBoundingFaces );
    }

    return numLoops;
}



rg_INT                  SurfaceModel::Body::find_bounding_edges(          list<SurfaceModel::Edge*>&   boundingEdges) const
{
    list<SurfaceModel::Shell*> boundingShells;
    find_bounding_shells( boundingShells );

    rg_INT numEdges = 0;
    list<SurfaceModel::Shell*>::iterator i_shell;
    for ( i_shell=boundingShells.begin(); i_shell!=boundingShells.end(); ++i_shell ) {
        numEdges += (*i_shell)->find_bounding_edges( boundingEdges );
    }

    return numEdges;
}



rg_INT                  SurfaceModel::Body::find_bounding_vertices(       list<SurfaceModel::Vertex*>& boundingVertices) const
{
    list<SurfaceModel::Shell*> boundingShells;
    find_bounding_shells( boundingShells );

    rg_INT numVertices = 0;
    list<SurfaceModel::Shell*>::iterator i_shell;
    for ( i_shell=boundingShells.begin(); i_shell!=boundingShells.end(); ++i_shell ) {
        numVertices += (*i_shell)->find_bounding_vertices( boundingVertices );
    }

    return numVertices;
}




///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Shell 
//
SurfaceModel::Shell::Shell()
: TopologicalEntity(),
  m_body( rg_NULL ),
  m_number_of_faces(0)
{
}



SurfaceModel::Shell::Shell(const rg_INT& ID)
: TopologicalEntity(ID),
  m_body( rg_NULL ),
  m_number_of_faces(0) 
{
}



SurfaceModel::Shell::Shell(const SurfaceModel::Shell& shell)
: TopologicalEntity( shell ),
  m_body(            shell.m_body ),
  m_faces(           shell.m_faces ),
  m_number_of_faces( shell.m_number_of_faces ) 
{
}



SurfaceModel::Shell::~Shell()
{
}



SurfaceModel::Shell& SurfaceModel::Shell::operator =(const SurfaceModel::Shell& shell)
{
    if ( this != &shell ) {
        TopologicalEntity::operator=(shell);
        m_body              = shell.m_body;
        m_faces             = shell.m_faces;
        m_number_of_faces   = shell.m_number_of_faces;
    }

    return *this;
}



SurfaceModel::Body*     SurfaceModel::Shell::body() const
{
    return m_body;
}



void                    SurfaceModel::Shell::set_body(const SurfaceModel::Body* const attachedBody)
{
    m_body = const_cast<SurfaceModel::Body*>( attachedBody );
}



const list< SurfaceModel::Face* >& SurfaceModel::Shell::faces() const
{
    return m_faces;
}



rg_INT                  SurfaceModel::Shell::number_of_faces() const
{
    return m_number_of_faces;
}



rg_INT                  SurfaceModel::Shell::number_of_edges() const
{
    list<SurfaceModel::Edge*> boundingEdges;
    return find_bounding_edges( boundingEdges );
}



rg_INT                  SurfaceModel::Shell::number_of_vertices() const
{
    list<SurfaceModel::Vertex*> boundingVertices;
    return find_bounding_vertices( boundingVertices );
}



void                    SurfaceModel::Shell::add_face(    const SurfaceModel::Face*  const face)
{
    m_faces.push_back( const_cast<SurfaceModel::Face*>( face ) );
    ++m_number_of_faces;
}



void                    SurfaceModel::Shell::remove_face( const SurfaceModel::Face*  const face)
{
    m_faces.remove( const_cast<SurfaceModel::Face*>( face ) );
    --m_number_of_faces;
}



bool                    SurfaceModel::Shell::is_exterior_shell() const
{
    return ( this == m_body->exterior_shell() ) ? true : false;
}


    
rg_INT      SurfaceModel::Shell::find_bounding_faces(          list<SurfaceModel::Face*>&   boundingFaces) const
{
    boundingFaces.insert( boundingFaces.end(), m_faces.begin(), m_faces.end() );
    return m_number_of_faces;
}



rg_INT      SurfaceModel::Shell::find_loops_of_bounding_faces( list<SurfaceModel::Loop*>&   loopsOfBoundingFaces) const
{
    int numLoops = 0;
    list<SurfaceModel::Face*>::const_iterator i_face;
    for ( i_face=m_faces.begin(); i_face!=m_faces.end(); ++i_face ) {
        numLoops += (*i_face)->find_bounding_loops( loopsOfBoundingFaces );
    }

    return numLoops;
}



rg_INT      SurfaceModel::Shell::find_bounding_edges(          list<SurfaceModel::Edge*>&   boundingEdges) const
{
    set<SurfaceModel::Edge*> setOfBoundingEdges;

    list<SurfaceModel::Face*>::const_iterator i_face;
    for ( i_face=m_faces.begin(); i_face!=m_faces.end(); ++i_face ) {
        SurfaceModel::Face* currFace = (*i_face);

        list<SurfaceModel::Edge*> boundingEdgesOfCurrFace;
        currFace->find_bounding_edges( boundingEdgesOfCurrFace );

        setOfBoundingEdges.insert( boundingEdgesOfCurrFace.begin(), boundingEdgesOfCurrFace.end() );
    }

    boundingEdges.insert( boundingEdges.end(), setOfBoundingEdges.begin(), setOfBoundingEdges.end() );

    return boundingEdges.size();
}



rg_INT      SurfaceModel::Shell::find_bounding_vertices(       list<SurfaceModel::Vertex*>& boundingVertices) const
{
    set<SurfaceModel::Vertex*> setOfBoundingVertices;

    list<SurfaceModel::Face*>::const_iterator i_face;
    for ( i_face=m_faces.begin(); i_face!=m_faces.end(); ++i_face ) {
        SurfaceModel::Face* currFace = (*i_face);

        list<SurfaceModel::Vertex*> boundingVtxsOfCurrFace;
        currFace->find_bounding_vertices( boundingVtxsOfCurrFace );

        setOfBoundingVertices.insert( boundingVtxsOfCurrFace.begin(), boundingVtxsOfCurrFace.end() );
    }

    boundingVertices.insert( boundingVertices.end(), setOfBoundingVertices.begin(), setOfBoundingVertices.end() );

    return boundingVertices.size();
}




///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Face 
//
SurfaceModel::Face::Face()
: WingedEdgeDataStructure::Face(),
  m_shell(rg_NULL),
  m_exteriorLoop(rg_NULL)
{
}



SurfaceModel::Face::Face(const rg_INT& ID)
: WingedEdgeDataStructure::Face(ID),
  m_shell(rg_NULL),
  m_exteriorLoop(rg_NULL)
{
}



SurfaceModel::Face::Face(const SurfaceModel::Face& face)
: WingedEdgeDataStructure::Face( face),
  m_shell(          face.m_shell ),
  m_exteriorLoop(   face.m_exteriorLoop ),
  m_interiorLoops(  face.m_interiorLoops )
{
}



SurfaceModel::Face::~Face()
{
}



SurfaceModel::Face& SurfaceModel::Face::operator =(const SurfaceModel::Face& face)
{
    if ( this != &face ) {
        WingedEdgeDataStructure::Face::operator=(face);

        m_shell         = face.m_shell;
        m_exteriorLoop  = face.m_exteriorLoop;
        m_interiorLoops = face.m_interiorLoops;
    }

    return *this;
}



SurfaceModel::Shell*                SurfaceModel::Face::shell() const
{
    return m_shell;
}



SurfaceModel::Loop*                 SurfaceModel::Face::exterior_loop() const
{
    return m_exteriorLoop;
}



const list< SurfaceModel::Loop* >&  SurfaceModel::Face::interior_loops() const
{
    return m_interiorLoops;
}



void    SurfaceModel::Face::set_shell(                const SurfaceModel::Shell* const attachedShell) 
{
    m_shell = const_cast<SurfaceModel::Shell*>( attachedShell );
}



void    SurfaceModel::Face::set_exterior_loop(        const SurfaceModel::Loop*  const exteriorLoop) 
{
    m_exteriorLoop = const_cast<SurfaceModel::Loop*>( exteriorLoop );
}



void    SurfaceModel::Face::add_interior_loop(    const SurfaceModel::Loop*  const interiorLoop) 
{
    m_interiorLoops.push_back( const_cast<SurfaceModel::Loop*>(interiorLoop) );
}



void    SurfaceModel::Face::remove_interior_loop( const SurfaceModel::Loop*  const interiorLoop) 
{
    m_interiorLoops.remove( const_cast<SurfaceModel::Loop*>(interiorLoop) );
}



bool    SurfaceModel::Face::no_interior_loop() const
{
    if ( m_interiorLoops.empty() ) {
        return true;
    }
    else {
        return false;
    }
}



rg_INT  SurfaceModel::Face::number_of_loops() const
{
    return 1+m_interiorLoops.size();
}



SurfaceModel::Edge*   SurfaceModel::Face::first_edge() const
{
    return m_exteriorLoop->first_edge();
}

    

SurfaceModel::Edge*   SurfaceModel::Face::find_edge(const SurfaceModel::Vertex* const vertex1, const SurfaceModel::Vertex* const vertex2 ) const
{
    return (SurfaceModel::Edge*) WingedEdgeDataStructure::Face::find_edge(vertex1, vertex2);
}



rg_INT  SurfaceModel::Face::number_of_bounding_vertices() const
{
    list<SurfaceModel::Vertex*> boundingVertices;
    return find_bounding_vertices( boundingVertices );
}



rg_INT  SurfaceModel::Face::number_of_bounding_edges() const
{
    list<SurfaceModel::Edge*> boundingEdges;
    return find_bounding_edges( boundingEdges );
}



rg_INT  SurfaceModel::Face::number_of_adjacent_faces() const
{
    list<SurfaceModel::Face*> adjacentFaces;
    return find_adjacent_faces( adjacentFaces );
}



rg_INT  SurfaceModel::Face::find_bounding_vertices( list<SurfaceModel::Vertex*>& boundingVertices) const
{
    list<SurfaceModel::Loop*> all_loops;
    find_bounding_loops( all_loops );

    list<SurfaceModel::Loop*>::iterator i_loop;
    for ( i_loop=all_loops.begin(); i_loop!=all_loops.end(); ++i_loop ) {
        SurfaceModel::Loop* currLoop = *i_loop;
        currLoop->find_bounding_vertices( boundingVertices );
    }

    return boundingVertices.size();
}



rg_INT  SurfaceModel::Face::find_bounding_edges(    list<SurfaceModel::Edge*>&   boundingEdges) const
{
    list<SurfaceModel::Loop*> all_loops;
    find_bounding_loops( all_loops );

    list<SurfaceModel::Loop*>::iterator i_loop;
    for ( i_loop=all_loops.begin(); i_loop!=all_loops.end(); ++i_loop ) {
        SurfaceModel::Loop* currLoop = *i_loop;
        currLoop->find_bounding_edges( boundingEdges );
    }

    return boundingEdges.size();
}



rg_INT  SurfaceModel::Face::find_bounding_loops(    list<SurfaceModel::Loop*>&   boundingLoops) const
{
    boundingLoops.push_back( m_exteriorLoop );
    boundingLoops.insert( boundingLoops.end(), m_interiorLoops.begin(), m_interiorLoops.end() );

    return boundingLoops.size();
}



rg_INT  SurfaceModel::Face::find_adjacent_faces(    list<SurfaceModel::Face*>&   adjacentFaces) const
{
    list<SurfaceModel::Loop*> all_loops;
    find_bounding_loops( all_loops );

    list<SurfaceModel::Loop*>::iterator i_loop;
    for ( i_loop=all_loops.begin(); i_loop!=all_loops.end(); ++i_loop ) {
        SurfaceModel::Loop* currLoop = *i_loop;
        currLoop->find_incident_faces( adjacentFaces );
    }

    set<SurfaceModel::Face*> setOfAdjacentFaces;
    setOfAdjacentFaces.insert( adjacentFaces.begin(), adjacentFaces.end() );

    adjacentFaces.clear();
    adjacentFaces.insert( adjacentFaces.begin(), setOfAdjacentFaces.begin(), setOfAdjacentFaces.end() );


    return adjacentFaces.size();
}



rg_INT  SurfaceModel::Face::find_bounding_edges_incident_to_vertex(const SurfaceModel::Vertex* const vertex, SurfaceModel::Edge*& prevEdge, SurfaceModel::Edge*& nextEdge) const
{
    vector<SurfaceModel::Edge*> boundingEdgesIncidentToVertex;

    list<SurfaceModel::Edge*> incidentEdges;
    vertex->find_incident_edges( incidentEdges );

    for (list<SurfaceModel::Edge*>::iterator i_edge=incidentEdges.begin(); i_edge!=incidentEdges.end(); ++i_edge ) {
        SurfaceModel::Edge* currEdge = *i_edge;

        if ( currEdge->is_member_of( this ) ) {
            boundingEdgesIncidentToVertex.push_back( currEdge );
        }
    }

    if ( vertex == boundingEdgesIncidentToVertex[0]->start_vertex() ) {
        if ( this == boundingEdgesIncidentToVertex[0]->left_face() ) {
            prevEdge = boundingEdgesIncidentToVertex[1];
            nextEdge = boundingEdgesIncidentToVertex[0];
        }
        else {
            prevEdge = boundingEdgesIncidentToVertex[0];
            nextEdge = boundingEdgesIncidentToVertex[1];
        }
    }
    else if ( vertex == boundingEdgesIncidentToVertex[0]->end_vertex() ) {
        if ( this == boundingEdgesIncidentToVertex[0]->left_face() ) {
            prevEdge = boundingEdgesIncidentToVertex[0];
            nextEdge = boundingEdgesIncidentToVertex[1];
        }
        else {
            prevEdge = boundingEdgesIncidentToVertex[1];
            nextEdge = boundingEdgesIncidentToVertex[0];
        }
    }
    else {
        prevEdge = rg_NULL;
        nextEdge = rg_NULL;
    }

    return boundingEdgesIncidentToVertex.size();
}


    
rg_INT  SurfaceModel::Face::find_edges_sharing_with_adjacent_face( const SurfaceModel::Face* const face, list<SurfaceModel::Edge*>& sharingEdges) const
{
    list<SurfaceModel::Edge*> boundingEdges;
    find_bounding_edges( boundingEdges );

    int numSharingEdges = 0;
    list<SurfaceModel::Edge*>::iterator i_edge;
    for ( i_edge=boundingEdges.begin(); i_edge!=boundingEdges.end(); ++i_edge ) {
        if ( face->is_incident_to( *i_edge ) ) {
            sharingEdges.push_back( *i_edge );
            ++numSharingEdges;
        }
    }

    return numSharingEdges;
}



///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Loop 
//
SurfaceModel::Loop::Loop()
: TopologicalEntity()
{
}



SurfaceModel::Loop::Loop(const rg_INT& ID)
: TopologicalEntity(ID)
{
}



SurfaceModel::Loop::Loop(const SurfaceModel::Loop& loop)
: TopologicalEntity( loop),
  m_face(            loop.m_face ),
  m_firstEdge(       loop.m_firstEdge )
{
}



SurfaceModel::Loop::~Loop()
{
}



SurfaceModel::Loop& SurfaceModel::Loop::operator =(const SurfaceModel::Loop& loop)
{
    if ( this != &loop ) {
        TopologicalEntity::operator=(loop);
        m_face      = loop.m_face;
        m_firstEdge = loop.m_firstEdge;
    }

    return *this;
}



SurfaceModel::Face*     SurfaceModel::Loop::face() const
{
    return m_face;
}



SurfaceModel::Edge*     SurfaceModel::Loop::first_edge() const
{
    return m_firstEdge;
}



void                    SurfaceModel::Loop::set_face(const SurfaceModel::Face* const attachedFace)
{
    m_face = const_cast<SurfaceModel::Face*>( attachedFace );
}



void                    SurfaceModel::Loop::set_first_edge(const SurfaceModel::Edge* const firstEdge)
{
    m_firstEdge = const_cast<SurfaceModel::Edge*>( firstEdge );
}


    
bool    SurfaceModel::Loop::is_exterior_loop() const
{
    if ( this == m_face->exterior_loop() ) {
        return true;
    }
    else {
        return false;
    }
}



rg_INT  SurfaceModel::Loop::find_bounding_vertices( list<SurfaceModel::Vertex*>& boundingVertices) const
{
    SurfaceModel::Edge* startEdge = m_firstEdge;
    SurfaceModel::Edge* currEdge  = startEdge;

    unsigned int numBoundingVertices = 0;

    do  {
        if ( m_face == currEdge->left_face() )  {
            boundingVertices.push_back( currEdge->end_vertex() );
            ++numBoundingVertices;

            currEdge = currEdge->left_hand_edge();
        }
        else if ( m_face == currEdge->right_face() )  {
            boundingVertices.push_back( currEdge->start_vertex() );
            ++numBoundingVertices;

            currEdge = currEdge->right_leg_edge();
        }
        else  {
            break;
        }
    } while ( currEdge != rg_NULL && currEdge != startEdge );   

    return numBoundingVertices;
}



rg_INT  SurfaceModel::Loop::find_bounding_edges(    list<SurfaceModel::Edge*>&   boundingEdges) const
{
    SurfaceModel::Edge* startEdge = m_firstEdge;
    SurfaceModel::Edge* currEdge  = startEdge;

    unsigned int numBoundingEdges = 0;
    do  {
        boundingEdges.push_back( currEdge );
        ++numBoundingEdges;

        if ( m_face == currEdge->left_face() ) {
            currEdge = currEdge->left_hand_edge();
        }
        else if ( m_face == currEdge->right_face() )  {
            currEdge = currEdge->right_leg_edge();
        }
        else {
            break;
        }
    } while ( currEdge != rg_NULL && currEdge != startEdge );   

    return numBoundingEdges;
}


    
rg_INT  SurfaceModel::Loop::find_incident_faces(    list<SurfaceModel::Face*>&   incidentFaces) const
{
    SurfaceModel::Edge* startEdge = m_firstEdge;
    SurfaceModel::Edge* currEdge  = startEdge;

    unsigned int numAdjacentFaces = 0;
    do  {
        SurfaceModel::Face* rightFace = currEdge->left_face();
        SurfaceModel::Face* leftFace  = currEdge->right_face();

        if ( m_face == leftFace )  {
            if ( rightFace != rg_NULL ) {
                incidentFaces.push_back( rightFace );
                ++numAdjacentFaces;
            }

            currEdge = currEdge->left_hand_edge();
        }
        else if ( m_face == rightFace )  {
            if ( leftFace != rg_NULL ) {
                incidentFaces.push_back( leftFace );
                ++numAdjacentFaces;
            }

            currEdge = currEdge->right_leg_edge();
        }
        else  {
            break;
        }
    } while ( currEdge != rg_NULL && currEdge != startEdge );   

    return numAdjacentFaces;

}


///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Edge 
//

SurfaceModel::Edge::Edge()
: WingedEdgeDataStructure::Edge()
{
}



SurfaceModel::Edge::Edge(const rg_INT& ID)
: WingedEdgeDataStructure::Edge(ID)
{
}



SurfaceModel::Edge::Edge(const SurfaceModel::Edge& edge)
: WingedEdgeDataStructure::Edge(edge)
{
}



SurfaceModel::Edge::~Edge()
{
}



SurfaceModel::Edge& SurfaceModel::Edge::operator =(const SurfaceModel::Edge& edge)
{
    if ( this != &edge ) {
        WingedEdgeDataStructure::Edge::operator=(edge);
    }

    return *this;
}


    
SurfaceModel::Shell*  SurfaceModel::Edge::shell() const
{
    if ( left_face() != rg_NULL ) {
        return left_face()->shell();
    }
    else if ( right_face() != rg_NULL ) {
        return right_face()->shell();
    }
    else {
        Shell* vtxShell = start_vertex()->shell();
        if ( vtxShell == rg_NULL ) {
            return end_vertex()->shell();
        }
        else {
            return vtxShell;
        }
    }
}



SurfaceModel::Vertex* SurfaceModel::Edge::start_vertex() const
{
    return (SurfaceModel::Vertex*)( WingedEdgeDataStructure::Edge::start_vertex() );
}



SurfaceModel::Vertex* SurfaceModel::Edge::end_vertex() const
{
    return (SurfaceModel::Vertex*)( WingedEdgeDataStructure::Edge::end_vertex() );
}



SurfaceModel::Face*   SurfaceModel::Edge::right_face() const
{
    return (SurfaceModel::Face*)( WingedEdgeDataStructure::Edge::right_face() );
}



SurfaceModel::Face*   SurfaceModel::Edge::left_face() const
{
    return (SurfaceModel::Face*)( WingedEdgeDataStructure::Edge::left_face() );
}



SurfaceModel::Edge*   SurfaceModel::Edge::right_hand_edge() const
{
    return (SurfaceModel::Edge*)( WingedEdgeDataStructure::Edge::right_hand_edge() );
}



SurfaceModel::Edge*   SurfaceModel::Edge::left_hand_edge() const
{
    return (SurfaceModel::Edge*)( WingedEdgeDataStructure::Edge::left_hand_edge() );
}



SurfaceModel::Edge*   SurfaceModel::Edge::right_leg_edge() const
{
    return (SurfaceModel::Edge*)( WingedEdgeDataStructure::Edge::right_leg_edge() );
}



SurfaceModel::Edge*   SurfaceModel::Edge::left_leg_edge() const
{
    return (SurfaceModel::Edge*)( WingedEdgeDataStructure::Edge::left_leg_edge() );
}



SurfaceModel::Vertex* SurfaceModel::Edge::opposite_vertex(const SurfaceModel::Vertex* const vertex) const
{
    return (SurfaceModel::Vertex*)( WingedEdgeDataStructure::Edge::opposite_vertex( vertex ) );
}

    

SurfaceModel::Face*   SurfaceModel::Edge::opposite_face(const SurfaceModel::Face* const face) const
{
    return (SurfaceModel::Face*)( WingedEdgeDataStructure::Edge::opposite_face( face ) );
}



rg_INT  SurfaceModel::Edge::find_bounding_vertices( list<SurfaceModel::Vertex*>& boundingVertices) const
{
    list< WingedEdgeDataStructure::Vertex* > verticesInWEDS;
    rg_INT numBoundingVertices = WingedEdgeDataStructure::Edge::find_bounding_vertices( verticesInWEDS );

    list< WingedEdgeDataStructure::Vertex* >::iterator i_vtx;
    for ( i_vtx=verticesInWEDS.begin(); i_vtx!=verticesInWEDS.end(); ++i_vtx ) {
        boundingVertices.push_back( (SurfaceModel::Vertex*)(*i_vtx) );
    }

    return numBoundingVertices;
}



rg_INT  SurfaceModel::Edge::find_adjacent_edges(    list<SurfaceModel::Edge*>&   adjacentEdges) const
{
    list< WingedEdgeDataStructure::Edge* > edgesInWEDS;
    rg_INT numAdjacentEdges = WingedEdgeDataStructure::Edge::find_adjacent_edges( edgesInWEDS );

    list< WingedEdgeDataStructure::Edge* >::iterator i_edge;
    for ( i_edge=edgesInWEDS.begin(); i_edge!=edgesInWEDS.end(); ++i_edge ) {
        adjacentEdges.push_back( (SurfaceModel::Edge*)(*i_edge) );
    }

    return numAdjacentEdges;
}



rg_INT  SurfaceModel::Edge::find_incident_faces(    list<SurfaceModel::Face*>&   incidentFaces) const
{
    list< WingedEdgeDataStructure::Face* > facesInWEDS;
    rg_INT numIncidentFaces = WingedEdgeDataStructure::Edge::find_incident_faces( facesInWEDS );

    list< WingedEdgeDataStructure::Face* >::iterator i_face;
    for ( i_face=facesInWEDS.begin(); i_face!=facesInWEDS.end(); ++i_face ) {
        incidentFaces.push_back( (SurfaceModel::Face*)(*i_face) );
    }

    return numIncidentFaces;
}


    
rg_INT  SurfaceModel::Edge::find_edges_in_star(     list<SurfaceModel::Edge*>&   edgesInStar) const
{
    list< WingedEdgeDataStructure::Edge* > edgesInStarInWEDS;
    rg_INT numEdgesInStar = WingedEdgeDataStructure::Edge::find_edges_in_star( edgesInStarInWEDS );

    list< WingedEdgeDataStructure::Edge* >::iterator i_edge;
    for ( i_edge=edgesInStarInWEDS.begin(); i_edge!=edgesInStarInWEDS.end(); ++i_edge ) {
        edgesInStar.push_back( (SurfaceModel::Edge*)(*i_edge) );
    }

    return numEdgesInStar;

}


///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Vertex 
//

SurfaceModel::Vertex::Vertex()
: WingedEdgeDataStructure::Vertex()
{
}



SurfaceModel::Vertex::Vertex(const rg_INT& ID)
: WingedEdgeDataStructure::Vertex(ID)
{
}



SurfaceModel::Vertex::Vertex(const SurfaceModel::Vertex& vertex)
: WingedEdgeDataStructure::Vertex(vertex)
{
}



SurfaceModel::Vertex::~Vertex()
{
}



SurfaceModel::Vertex& SurfaceModel::Vertex::operator =(const Vertex& vertex)
{
    if ( this != &vertex ) {
        WingedEdgeDataStructure::Vertex::operator=(vertex);
    }

    return *this;
}



SurfaceModel::Shell*  SurfaceModel::Vertex::shell() const
{
    SurfaceModel::Edge* firstEdge = first_edge();
    if ( firstEdge == rg_NULL ) {
        return rg_NULL;
    }

    SurfaceModel::Shell* containingShell = rg_NULL;
    if ( firstEdge->left_face() != rg_NULL ) {
        containingShell = firstEdge->left_face()->shell();
    }
    else if ( firstEdge->right_face() != rg_NULL ) {
        containingShell = firstEdge->right_face()->shell();
    }
    else {
        SurfaceModel::Vertex* oppositeVtx = firstEdge->opposite_vertex( this );

        list<SurfaceModel::Face*> incidentFaces;
        int numFaces = oppositeVtx->find_incident_faces( incidentFaces );

        if ( numFaces != 0 ) {
            containingShell = incidentFaces.front()->shell();
        }
    }

    return containingShell;
}



SurfaceModel::Edge*   SurfaceModel::Vertex::first_edge() const
{
    return (SurfaceModel::Edge*)( WingedEdgeDataStructure::Vertex::first_edge() );
}



rg_INT  SurfaceModel::Vertex::find_adjacent_vertices( list<SurfaceModel::Vertex*>& adjacentVertices) const
{
    list< WingedEdgeDataStructure::Vertex* > verticesInWEDS;
    rg_INT numAdjacentVertices = WingedEdgeDataStructure::Vertex::find_adjacent_vertices( verticesInWEDS );

    list< WingedEdgeDataStructure::Vertex* >::iterator i_vtx;
    for ( i_vtx=verticesInWEDS.begin(); i_vtx!=verticesInWEDS.end(); ++i_vtx ) {
        adjacentVertices.push_back( (SurfaceModel::Vertex*)(*i_vtx) );
    }

    return numAdjacentVertices;
}



rg_INT  SurfaceModel::Vertex::find_incident_edges(    list<SurfaceModel::Edge*>&   incidentEdges) const
{
    list< WingedEdgeDataStructure::Edge* > edgesInWEDS;
    rg_INT numIncidentEdges = WingedEdgeDataStructure::Vertex::find_incident_edges( edgesInWEDS );

    list< WingedEdgeDataStructure::Edge* >::iterator i_edge;
    for ( i_edge=edgesInWEDS.begin(); i_edge!=edgesInWEDS.end(); ++i_edge ) {
        incidentEdges.push_back( (SurfaceModel::Edge*)(*i_edge) );
    }

    return numIncidentEdges;
}



rg_INT  SurfaceModel::Vertex::find_incident_faces(    list<SurfaceModel::Face*>&   incidentFaces) const
{
    list< WingedEdgeDataStructure::Face* > facesInWEDS;
    rg_INT numIncidentFaces = WingedEdgeDataStructure::Vertex::find_incident_faces( facesInWEDS );

    list< WingedEdgeDataStructure::Face* >::iterator i_face;
    for ( i_face=facesInWEDS.begin(); i_face!=facesInWEDS.end(); ++i_face ) {
        incidentFaces.push_back( (SurfaceModel::Face*)(*i_face) );
    }

    return numIncidentFaces;
}


    
rg_INT  SurfaceModel::Vertex::find_edges_in_star(     list<SurfaceModel::Edge*>&   edgesInStar) const
{
    list<WingedEdgeDataStructure::Edge*> edgesInStarWEDS;
    rg_INT  numEdgesInStar = WingedEdgeDataStructure::Vertex::find_edges_in_star( edgesInStarWEDS );

    list<WingedEdgeDataStructure::Edge*>::iterator i_edge;
    for ( i_edge=edgesInStarWEDS.begin(); i_edge!=edgesInStarWEDS.end(); ++i_edge ) {
        edgesInStar.push_back( (SurfaceModel::Edge*)(*i_edge) );
    }

    return numEdgesInStar;
}



rg_INT  SurfaceModel::Vertex::find_edges_in_shell(    list<SurfaceModel::Edge*>&   edgesInShell) const
{
    list<WingedEdgeDataStructure::Edge*> edgesInShellWEDS;
    rg_INT  numEdgesInShell = WingedEdgeDataStructure::Vertex::find_edges_in_shell( edgesInShellWEDS );

    list<WingedEdgeDataStructure::Edge*>::iterator i_edge;
    for ( i_edge=edgesInShellWEDS.begin(); i_edge!=edgesInShellWEDS.end(); ++i_edge ) {
        edgesInShell.push_back( (SurfaceModel::Edge*)(*i_edge) );
    }

    return numEdgesInShell;
}

