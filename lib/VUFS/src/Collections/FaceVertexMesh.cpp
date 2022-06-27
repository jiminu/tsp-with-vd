#include "FaceVertexMesh.h"
#include "OperationsForPolygonMeshModel.h"

#include <algorithm>
#include <map>
using namespace std;



FaceVertexMesh::FaceVertexMesh()
: m_number_of_faces(0), 
  m_number_of_vertices(0)
{
}



FaceVertexMesh::FaceVertexMesh(const FaceVertexMesh& FVmesh)
{
}



FaceVertexMesh::~FaceVertexMesh()
{
    clear();
}



FaceVertexMesh& FaceVertexMesh::operator =(const FaceVertexMesh& FVmesh)
{
    return *this;
}



void            FaceVertexMesh::clear()
{
    for (list<Face*>::iterator i_face=m_faces.begin(); i_face!=m_faces.end(); ++i_face ) {
        delete (*i_face);
    }

    for (list<Vertex*>::iterator i_vx=m_vertices.begin(); i_vx!=m_vertices.end(); ++i_vx ) {
        delete (*i_vx);
    }

    m_faces.clear();
    m_vertices.clear();

    m_number_of_faces    = 0;
    m_number_of_vertices = 0;
}



unsigned int    FaceVertexMesh::number_of_faces() const
{
    return m_number_of_faces;
}



unsigned int    FaceVertexMesh::number_of_vertices() const
{
    return m_number_of_vertices;
}



const list<FaceVertexMesh::Face*>&     FaceVertexMesh::get_faces() const
{
    return m_faces;
}



const list<FaceVertexMesh::Vertex*>&   FaceVertexMesh::get_vertices() const
{
    return m_vertices;
}



FaceVertexMesh::Face*           FaceVertexMesh::create_face()
{
    m_faces.push_back( new Face() );
    ++m_number_of_faces;

    return m_faces.back();
}



FaceVertexMesh::Vertex*         FaceVertexMesh::create_vertex()
{
    m_vertices.push_back( new Vertex() );
    ++m_number_of_vertices;

    return m_vertices.back();
}



void            FaceVertexMesh::remove_face(const FaceVertexMesh::Face* const face)
{
    if ( face != rg_NULL ) {
        m_faces.remove( const_cast<Face*>( face ) );
        --m_number_of_faces;
        delete face;
    }
}



void            FaceVertexMesh::remove_vertex(const FaceVertexMesh::Vertex* const vertex)
{
    if ( vertex != rg_NULL ) {
        m_vertices.remove( const_cast<Vertex*>( vertex ) );
        --m_number_of_vertices;
        delete vertex;
    }
}



FaceVertexMesh::Vertex*  FaceVertexMesh::find(const rg_Point3D& point) const
{
    Vertex* correspondingVertex = rg_NULL;
    list<Vertex*>::const_iterator i_vtx;
    for ( i_vtx=m_vertices.begin(); i_vtx!=m_vertices.end(); ++i_vtx ) {
        Vertex* currVtx = *i_vtx;

        if ( currVtx->has_same_coordinate( point ) ) {
            correspondingVertex = currVtx;
            break;
        }
    }

    return correspondingVertex;
}



void            FaceVertexMesh::check_face_with_same_vertex( list<FaceVertexMesh::Face*>& facesWithSameVertex ) const
{
    for (list<Face*>::const_iterator i_face=m_faces.begin(); i_face!=m_faces.end(); ++i_face ) {
        if ( (*i_face)->has_same_vertex() ) {
            facesWithSameVertex.push_back( (*i_face) );
        }
    }
}

    

void            FaceVertexMesh::check_non_manifold_vertices( list<FaceVertexMesh::Vertex*>& nonManifoldVertices ) const
{
    for (list<Vertex*>::const_iterator i_vtx=m_vertices.begin(); i_vtx!=m_vertices.end(); ++i_vtx ) {
        Vertex* currVtx = (*i_vtx);

        list<Face*> incidentFaces;
        int numIncidentFaces = currVtx->find_incident_faces( incidentFaces );


        int     countTracingFace = 0;
        Face*   startFace = incidentFaces.front();
        Face*   currFace  = startFace;

        do {
            ++countTracingFace;

            Edge edge[2];
            bool    bFound = currFace->find_edges_incident_to_vertex( currVtx, edge[0], edge[1] );

            currFace = edge[0].find_twin_face();

        } while ( currFace!=rg_NULL && currFace!=startFace );

        if ( countTracingFace != numIncidentFaces || currFace==rg_NULL ) {
            nonManifoldVertices.push_back( currVtx );
        }
    }
}



bool            FaceVertexMesh::convert_to(PolygonMeshModel& polygonMeshModel) const
{
    polygonMeshModel.clear();


    map<FaceVertexMesh::Face*, PolygonMeshModel::Face*>     faceMap;    
    map<FaceVertexMesh::Vertex*, PolygonMeshModel::Vertex*> vertexMap;

    make_map_from_FaceVertexMesh_to_PolygonMesh( polygonMeshModel, faceMap, vertexMap );

    construct_shells_of_PolygonMesh( polygonMeshModel, faceMap, vertexMap );

    construct_bodies_of_PolygonMesh( polygonMeshModel );

    return true;
}



void            FaceVertexMesh::make_map_from_FaceVertexMesh_to_PolygonMesh(
                                    PolygonMeshModel&                                        polygonMeshModel,
                                    map<FaceVertexMesh::Face*, PolygonMeshModel::Face*>&     faceMap,
                                    map<FaceVertexMesh::Vertex*, PolygonMeshModel::Vertex*>& vertexMap  ) const
{
    list<FaceVertexMesh::Face*>::const_iterator i_face;
    for ( i_face=m_faces.begin(); i_face!=m_faces.end(); ++i_face ) {
        FaceVertexMesh::Face*   faceFVMesh = *i_face;
        PolygonMeshModel::Face* facePMesh  = polygonMeshModel.create( PolygonMeshModel::Face() );

        facePMesh->setID(      faceFVMesh->getID() );
        facePMesh->set_normal( faceFVMesh->normal() );

        faceMap.insert( make_pair( faceFVMesh, facePMesh ) );
    }
    faceMap.insert( make_pair( (FaceVertexMesh::Face*)rg_NULL, (PolygonMeshModel::Face*)rg_NULL ) );



    list<FaceVertexMesh::Vertex*>::const_iterator i_vtx;
    for ( i_vtx=m_vertices.begin(); i_vtx!=m_vertices.end(); ++i_vtx ) {
        FaceVertexMesh::Vertex*   vtxFVMesh = *i_vtx;
        PolygonMeshModel::Vertex* vtxPMesh  = polygonMeshModel.create( PolygonMeshModel::Vertex() );

        vtxPMesh->setID(          vtxFVMesh->getID() );
        vtxPMesh->set_coordinate( vtxFVMesh->coordinate() );

        vertexMap.insert( make_pair( vtxFVMesh, vtxPMesh ) );
    }
    vertexMap.insert( make_pair( (FaceVertexMesh::Vertex*)rg_NULL, (PolygonMeshModel::Vertex*)rg_NULL ) );

}




void            FaceVertexMesh::construct_shells_of_PolygonMesh(
                                PolygonMeshModel&                                        polygonMeshModel,
                                map<FaceVertexMesh::Face*, PolygonMeshModel::Face*>&     faceMap,
                                map<FaceVertexMesh::Vertex*, PolygonMeshModel::Vertex*>& vertexMap  ) const
{
    int shellID = 0;
    int edgeID  = 0;

    set<FaceVertexMesh::Face*> visitedFacesInFVMesh;

    list<FaceVertexMesh::Face*>::const_iterator i_face;
    for ( i_face=m_faces.begin(); i_face!=m_faces.end(); ++i_face ) {
        FaceVertexMesh::Face*   faceFVMesh = *i_face;

        if ( visitedFacesInFVMesh.find( faceFVMesh ) != visitedFacesInFVMesh.end() ) {
            continue;
        }


        PolygonMeshModel::Shell* newShell = polygonMeshModel.create( PolygonMeshModel::Shell() );
        newShell->setID( ++shellID );


        list<FaceVertexMesh::Face*> faceStack;
        faceStack.push_back( faceFVMesh );

        while ( !faceStack.empty() ) {
            FaceVertexMesh::Face* currFaceFVM = faceStack.front();
            faceStack.pop_front();

            if ( currFaceFVM == rg_NULL ) {
                continue;
            }

            if ( visitedFacesInFVMesh.find( currFaceFVM ) != visitedFacesInFVMesh.end() ) {
                continue;
            }
            else {
                visitedFacesInFVMesh.insert( currFaceFVM );
            }



            PolygonMeshModel::Face* currFacePM  = faceMap.find( currFaceFVM )->second;


            //  set topology of shell <-> face
            currFacePM->set_shell( newShell );
            newShell->add_face( currFacePM );


            //  create and find bounding edge of currFacePM
            //  and set topology of edge <-> vertex and edge -> face
            PolygonMeshModel::Edge* boundingEdgePM[3] = {rg_NULL, rg_NULL, rg_NULL};
            for ( int pos=0; pos<3; ++pos ) {
                FaceVertexMesh::Edge  boundingEdgeFVM = currFaceFVM->bounding_edge( pos );
                FaceVertexMesh::Face* adjFaceFVM      = currFaceFVM->adjacent_face( pos );
                if ( adjFaceFVM == rg_NULL ) {
                    int stop = 1;
                }


                PolygonMeshModel::Vertex* startVtxPM  = vertexMap.find( boundingEdgeFVM.start_vertex() )->second;
                PolygonMeshModel::Vertex* endVtxPM    = vertexMap.find( boundingEdgeFVM.end_vertex() )->second;

                if ( visitedFacesInFVMesh.find( adjFaceFVM ) == visitedFacesInFVMesh.end() ) {
                    //  create edge of currFacePM in PolygonMeshModel
                    boundingEdgePM[pos] = polygonMeshModel.create( PolygonMeshModel::Edge() );
                    boundingEdgePM[pos]->setID( ++edgeID );
                    boundingEdgePM[pos]->set_start_vertex( startVtxPM );
                    boundingEdgePM[pos]->set_end_vertex(   endVtxPM   );
                    startVtxPM->set_first_edge( boundingEdgePM[pos] );
                    endVtxPM->set_first_edge( boundingEdgePM[pos] );

                    boundingEdgePM[pos]->set_left_face(    currFacePM );

                    //  push back adjFaceFVM for next search
                    faceStack.push_back( adjFaceFVM );                  
                }
                else {
                    //  find edge of currFacePM from adjacent face.
                    PolygonMeshModel::Face* adjFacePM  = faceMap.find( adjFaceFVM )->second;

                    boundingEdgePM[pos] = adjFacePM->find_edge( startVtxPM, endVtxPM ); 
                    boundingEdgePM[pos]->set_right_face(    currFacePM );
                }
            }


            //   set topology of edge <-> edge
            if ( currFacePM == boundingEdgePM[0]->left_face() ) {
                boundingEdgePM[0]->set_left_hand_edge(  boundingEdgePM[1] );
                boundingEdgePM[0]->set_left_leg_edge(   boundingEdgePM[2] );
            }
            else if ( currFacePM == boundingEdgePM[0]->right_face() ) {
                boundingEdgePM[0]->set_right_leg_edge(  boundingEdgePM[1] );
                boundingEdgePM[0]->set_right_hand_edge( boundingEdgePM[2] );
            }

            if ( currFacePM == boundingEdgePM[1]->left_face() ) {
                boundingEdgePM[1]->set_left_hand_edge(  boundingEdgePM[2] );
                boundingEdgePM[1]->set_left_leg_edge(   boundingEdgePM[0] );
            }
            else if ( currFacePM == boundingEdgePM[1]->right_face() ) {
                boundingEdgePM[1]->set_right_leg_edge(  boundingEdgePM[2] );
                boundingEdgePM[1]->set_right_hand_edge( boundingEdgePM[0] );
            }

            if ( currFacePM == boundingEdgePM[2]->left_face() ) {
                boundingEdgePM[2]->set_left_hand_edge(  boundingEdgePM[0] );
                boundingEdgePM[2]->set_left_leg_edge(   boundingEdgePM[1] );
            }
            else if ( currFacePM == boundingEdgePM[2]->right_face() ) {
                boundingEdgePM[2]->set_right_leg_edge(  boundingEdgePM[0] );
                boundingEdgePM[2]->set_right_hand_edge( boundingEdgePM[1] );
            }

            //   create loop and set topology of loop -> edge and loop <-> face
            PolygonMeshModel::Loop* loopOfCurrFacePM = polygonMeshModel.create( PolygonMeshModel::Loop() );
            loopOfCurrFacePM->setID( currFacePM->getID() );            
            loopOfCurrFacePM->set_first_edge( boundingEdgePM[0] );

            loopOfCurrFacePM->set_face( currFacePM );
            currFacePM->set_exterior_loop( loopOfCurrFacePM );

            currFacePM->set_first_edge( boundingEdgePM[0] );



        }
    }
}



void            FaceVertexMesh::construct_bodies_of_PolygonMesh(
                                PolygonMeshModel&                                        polygonMeshModel ) const
{
    int bodyID = 0;
    if ( polygonMeshModel.number_of_shells() == 1 ) {
        PolygonMeshModel::Shell* oneShell = (PolygonMeshModel::Shell*) polygonMeshModel.get_all_shells().front();
        PolygonMeshModel::Body*  newBody  = polygonMeshModel.create( PolygonMeshModel::Body() );

        newBody->setID( ++bodyID );
        newBody->set_exterior_shell( oneShell );

        oneShell->set_body( newBody );
    }
    else {
        //  When there are many shells, 
        //  we assume that there are no intersection among shells.
        //  The following codes are incorrect and must be fixed.

        list<PolygonMeshModel::Shell*> all_shells;
        polygonMeshModel.all_shells( all_shells );

        list<PolygonMeshModel::Shell*>::iterator i_shell;
        for ( i_shell=all_shells.begin(); i_shell!=all_shells.end(); ++i_shell ) {
            PolygonMeshModel::Shell* currShell = *i_shell;

            PolygonMeshModel::Body*  newBody  = polygonMeshModel.create( PolygonMeshModel::Body() );

            newBody->setID( ++bodyID );
            newBody->set_exterior_shell( currShell );

            currShell->set_body( newBody );
        }


        //vector< OperationsForPolygonMeshModel::ShellWithVolume > sortedShells( all_shells.size() );
        //int pos = 0;

        //list<PolygonMeshModel::Shell*>::iterator i_shell;
        //for ( i_shell=all_shells.begin(); i_shell!=all_shells.end(); ++i_shell, ++pos ) {
        //    PolygonMeshModel::Shell* currShell = *i_shell;

        //    sortedShells[ pos ].m_shell = currShell;
        //    sortedShells[ pos ].m_volume = currShell->bounding_box().computeVolume();
        //}
        //sort( sortedShells.begin(), sortedShells.end(), OperationsForPolygonMeshModel::GreaterThanShellWithVolume );



        //map< PolygonMeshModel::Shell*, list<PolygonMeshModel::Shell*> > shellHierarchy;


    }
}




///////////////////////////////////////////////////////////////////////////////
//
// class FaceVertexMesh::Face 
//
FaceVertexMesh::Face::Face()
: TopologicalEntity()
{
    m_bounding_vertex[0] = rg_NULL;
    m_bounding_vertex[1] = rg_NULL;
    m_bounding_vertex[2] = rg_NULL;
}



FaceVertexMesh::Face::Face(const rg_INT& ID)
: TopologicalEntity( ID )
{
    m_bounding_vertex[0] = rg_NULL;
    m_bounding_vertex[1] = rg_NULL;
    m_bounding_vertex[2] = rg_NULL;
}



FaceVertexMesh::Face::Face(const FaceVertexMesh::Face& face)
: TopologicalEntity( face )
{
    m_bounding_vertex[0] = face.m_bounding_vertex[0];
    m_bounding_vertex[1] = face.m_bounding_vertex[1];
    m_bounding_vertex[2] = face.m_bounding_vertex[2];
}



FaceVertexMesh::Face::~Face()
{
}



FaceVertexMesh::Face& FaceVertexMesh::Face::operator =(const FaceVertexMesh::Face& face)
{
    if ( this != &face ) {
        TopologicalEntity::operator=( face );
        m_bounding_vertex[0] = face.m_bounding_vertex[0];
        m_bounding_vertex[1] = face.m_bounding_vertex[1];
        m_bounding_vertex[2] = face.m_bounding_vertex[2];
    }

    return *this;
}



FaceVertexMesh::Vertex*  FaceVertexMesh::Face::bounding_vertex(const int& pos) const
{
    return ( pos>=0 && pos<3 ) ? m_bounding_vertex[pos]: rg_NULL;
}

    

rg_Point3D      FaceVertexMesh::Face::normal() const
{
    return m_normal;
}



void  FaceVertexMesh::Face::set_bounding_vertex(const int& pos, const FaceVertexMesh::Vertex* const vertex)
{
    if ( pos>=0 && pos<3 ) {
        m_bounding_vertex[pos] = const_cast<FaceVertexMesh::Vertex*>( vertex );
    }
}


    
void    FaceVertexMesh::Face::set_normal(const rg_Point3D& normal)
{
    m_normal = normal;
}



bool    FaceVertexMesh::Face::is_incident_to(const FaceVertexMesh::Vertex* const vertex) const
{
    bool is_incident = false;
    for ( int pos=0; pos<3; ++pos ) {
        if ( m_bounding_vertex[pos] == vertex ) {
            is_incident = true;
            break;
        }
    }

    return is_incident;
}

    

bool    FaceVertexMesh::Face::is_incident_to(const Edge& edge) const
{
    if ( this->is_incident_to( edge.start_vertex() ) && this->is_incident_to( edge.end_vertex() ) ) {
        return true;
    }
    else {
        return false;
    }
}

    
bool            FaceVertexMesh::Face::has_same_vertex() const
{
    if (    (m_bounding_vertex[0] == m_bounding_vertex[1])
         || (m_bounding_vertex[1] == m_bounding_vertex[2])
         || (m_bounding_vertex[2] == m_bounding_vertex[0]) ) 
    {
        return true;
    }
    else {
        return false;
    }
}


FaceVertexMesh::Edge            FaceVertexMesh::Face::bounding_edge(const int& pos) const
{
    if ( pos == 0 ) {
        return Edge( this, m_bounding_vertex[0], m_bounding_vertex[1] );
    }
    else if ( pos == 1 ) {
        return Edge( this, m_bounding_vertex[1], m_bounding_vertex[2] );
    }
    else if ( pos == 2 ) {
        return Edge( this, m_bounding_vertex[2], m_bounding_vertex[0] );
    }
    else {
        return Edge();
    }
}

    
    
FaceVertexMesh::Face*           FaceVertexMesh::Face::adjacent_face(const int& pos) const
{
    Edge edge = bounding_edge( pos );
    return edge.find_twin_face();
}



int             FaceVertexMesh::Face::find_edge_position(const FaceVertexMesh::Edge& edge) const
{
    if ( this != edge.face() ) {
        return -1;
    }

    int pos = -1;
    if ( m_bounding_vertex[0] == edge.start_vertex() && m_bounding_vertex[1] == edge.end_vertex()) {
        pos = 0;
    }
    else if ( m_bounding_vertex[1] == edge.start_vertex() && m_bounding_vertex[2] == edge.end_vertex()) {
        pos = 1;
    }
    else if ( m_bounding_vertex[2] == edge.start_vertex() && m_bounding_vertex[0] == edge.end_vertex()) {
        pos = 2;
    }

    return pos;
}

    

int             FaceVertexMesh::Face::find_vertex_position(const FaceVertexMesh::Vertex* const vertex) const
{
    int pos = -1;
    if ( m_bounding_vertex[0] == vertex ) {
        pos = 0;
    }
    else if ( m_bounding_vertex[1] == vertex ) {
        pos = 1;
    }
    else if ( m_bounding_vertex[2] == vertex ) {
        pos = 2;
    }

    return pos;
}



unsigned int    FaceVertexMesh::Face::number_of_bounding_vertices() const
{
    return 3;
}



unsigned int    FaceVertexMesh::Face::find_bounding_vertices(list<FaceVertexMesh::Vertex*>& boundingVertices) const
{
    boundingVertices.push_back( m_bounding_vertex[0] );
    boundingVertices.push_back( m_bounding_vertex[1] );
    boundingVertices.push_back( m_bounding_vertex[2] );

    return 3;
}



unsigned int    FaceVertexMesh::Face::find_bounding_edges(    list<FaceVertexMesh::Edge>&    boundingEdges    ) const
{
    for ( int pos=0; pos<3; ++pos ) {
        boundingEdges.push_back( bounding_edge( pos ) );
    }

    return 3;
}


    
unsigned int    FaceVertexMesh::Face::find_adjacent_faces(    list<FaceVertexMesh::Face*>& adjacentFaces      ) const
{
    for ( int pos=0; pos<3; ++pos ) {
        Edge edge = bounding_edge( pos );

        Face* adjFace = edge.find_twin_face();
        if ( adjFace != rg_NULL ) {
            adjacentFaces.push_back( adjFace );
        }
    }

    return adjacentFaces.size();
}

    

bool            FaceVertexMesh::Face::find_edges_incident_to_vertex(const FaceVertexMesh::Vertex* const vertex, FaceVertexMesh::Edge& inEdge, FaceVertexMesh::Edge& outEdge) const
{
    bool    isFound = true;
    int vtxPos = find_vertex_position( vertex );
    
    if ( vtxPos == 0 ) {
        inEdge.set( this, m_bounding_vertex[2], m_bounding_vertex[0] );
        outEdge.set( this, m_bounding_vertex[0], m_bounding_vertex[1] );
    }
    else if ( vtxPos == 1 ) {
        inEdge.set( this, m_bounding_vertex[0], m_bounding_vertex[1] );
        outEdge.set( this, m_bounding_vertex[1], m_bounding_vertex[2] );
    }
    else if ( vtxPos == 2 ) {
        inEdge.set( this, m_bounding_vertex[1], m_bounding_vertex[2] );
        outEdge.set( this, m_bounding_vertex[2], m_bounding_vertex[0] );
    }
    else {
        isFound = false;
    }

    return isFound;
}



///////////////////////////////////////////////////////////////////////////////
//
// class FaceVertexMesh::Edge 
//
FaceVertexMesh::Edge::Edge()
: m_face(         rg_NULL ), 
  m_start_vertex( rg_NULL ), 
  m_end_vertex(   rg_NULL )
{
}


    
FaceVertexMesh::Edge::Edge(const FaceVertexMesh::Face* const face, const FaceVertexMesh::Vertex* const startVertex, const FaceVertexMesh::Vertex* const endVertex )
: m_face(         const_cast<Face*>(face) ), 
  m_start_vertex( const_cast<Vertex*>(startVertex) ), 
  m_end_vertex(   const_cast<Vertex*>(endVertex)   )
{
}



FaceVertexMesh::Edge::Edge(const FaceVertexMesh::Edge& edge)
: m_face(         edge.m_face ), 
  m_start_vertex( edge.m_start_vertex ), 
  m_end_vertex(   edge.m_end_vertex )
{
}



FaceVertexMesh::Edge::~Edge()
{
}



FaceVertexMesh::Edge& FaceVertexMesh::Edge::operator =(const FaceVertexMesh::Edge& edge)
{
    if ( this != &edge ) {
        m_face         = edge.m_face;
        m_start_vertex = edge.m_start_vertex;
        m_end_vertex   = edge.m_end_vertex;
    }

    return *this;
}



FaceVertexMesh::Face*    FaceVertexMesh::Edge::face() const
{
    return m_face;
}



FaceVertexMesh::Vertex*  FaceVertexMesh::Edge::start_vertex() const
{
    return m_start_vertex;
}



FaceVertexMesh::Vertex*  FaceVertexMesh::Edge::end_vertex() const
{
    return m_end_vertex;
}



void                     FaceVertexMesh::Edge::set_face( const FaceVertexMesh::Face* const face )
{
    m_face = const_cast<Face*>( face );
}



void                     FaceVertexMesh::Edge::set_start_vertex( const FaceVertexMesh::Vertex* const vertex )
{
    m_start_vertex = const_cast<Vertex*>( vertex );
}



void                     FaceVertexMesh::Edge::set_end_vertex( const FaceVertexMesh::Vertex* const vertex )
{
    m_end_vertex = const_cast<Vertex*>( vertex );
}



void                     FaceVertexMesh::Edge::set( const FaceVertexMesh::Face* const   face,
                                                    const FaceVertexMesh::Vertex* const startVertex, 
                                                    const FaceVertexMesh::Vertex* const endVertex )
{
    m_face = const_cast<Face*>(           face );
    m_start_vertex = const_cast<Vertex*>( startVertex );
    m_end_vertex = const_cast<Vertex*>(   endVertex );
}


    
FaceVertexMesh::Face*       FaceVertexMesh::Edge::find_twin_face() const
{
    if ( m_face == rg_NULL ) {
        return rg_NULL;
    }

    list<Face*> incidentFaces;
    unsigned int numIncidentFaces = this->find_incident_faces( incidentFaces );

    bool doesExistThisFace = false;
    for( list<Face*>::iterator i_face=incidentFaces.begin(); i_face!=incidentFaces.end(); ++i_face ) {
        if ( m_face == (*i_face) ) {
            doesExistThisFace = true;
            incidentFaces.erase( i_face );
            break;
        }
    }

    if ( doesExistThisFace && incidentFaces.size()==1 ) {
        return incidentFaces.front();
    }
    else {
        return rg_NULL;
    }
}

    

unsigned int FaceVertexMesh::Edge::find_incident_faces( list<FaceVertexMesh::Face*>& incidentFaces ) const
{
    list<Face*> facesIncidentToStartVtx;
    m_start_vertex->find_incident_faces( facesIncidentToStartVtx );

    for ( list<Face*>::iterator i_face=facesIncidentToStartVtx.begin(); i_face!=facesIncidentToStartVtx.end(); ++i_face ) {
        Face* currFace = *i_face;

        if ( currFace->is_incident_to( m_end_vertex ) ) {
            incidentFaces.push_back( currFace );
        }
    }

    return incidentFaces.size();
}

///////////////////////////////////////////////////////////////////////////////
//
// class FaceVertexMesh::Vertex 
//
FaceVertexMesh::Vertex::Vertex()
: TopologicalEntity()
{
}



FaceVertexMesh::Vertex::Vertex(const rg_INT& ID)
: TopologicalEntity( ID )
{
}



FaceVertexMesh::Vertex::Vertex(const FaceVertexMesh::Vertex& vertex)
: TopologicalEntity( vertex ),
  m_coordinate( vertex.m_coordinate ),
  m_incident_faces( vertex.m_incident_faces )
{
}



FaceVertexMesh::Vertex::~Vertex()
{
}



FaceVertexMesh::Vertex& FaceVertexMesh::Vertex::operator =(const FaceVertexMesh::Vertex& vertex)
{
    if ( this != &vertex ) {
        TopologicalEntity::operator=( vertex );
        m_coordinate     = vertex.m_coordinate;
        m_incident_faces = vertex.m_incident_faces;
    }

    return *this;
}



rg_Point3D      FaceVertexMesh::Vertex::coordinate() const
{
    return m_coordinate;
}



void            FaceVertexMesh::Vertex::set_coordinate( const rg_Point3D& coordinate )
{
    m_coordinate = coordinate;
}



void            FaceVertexMesh::Vertex::add_incident_face( const FaceVertexMesh::Face* const face )
{
    m_incident_faces.insert( const_cast<FaceVertexMesh::Face*>( face ) );
}



void            FaceVertexMesh::Vertex::remove_incident_face( const FaceVertexMesh::Face* const face )
{
    m_incident_faces.erase( const_cast<FaceVertexMesh::Face*>( face ) );
}



bool            FaceVertexMesh::Vertex::has_same_coordinate(const rg_Point3D& point, const rg_REAL& res)
{
    return ( m_coordinate.isEqual( point, res ) ) ? true : false;
}



bool            FaceVertexMesh::Vertex::is_equal(const Vertex* const vertex, const rg_REAL& res)
{
    return ( m_coordinate.isEqual( vertex->m_coordinate, res ) ) ? true : false;
}



bool            FaceVertexMesh::Vertex::is_member_of( const FaceVertexMesh::Face* const face ) const
{
    return face->is_incident_to( this );
}


    
unsigned int    FaceVertexMesh::Vertex::number_of_incident_faces() const
{
    return m_incident_faces.size();
}



unsigned int    FaceVertexMesh::Vertex::find_incident_faces(list<FaceVertexMesh::Face*>& incidentFaces) const
{
    incidentFaces.insert( incidentFaces.end(), m_incident_faces.begin(), m_incident_faces.end() );

    return m_incident_faces.size();
}




