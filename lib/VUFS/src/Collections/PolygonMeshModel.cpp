#include "PolygonMeshModel.h"

#include "Plane.h"
#include "rg_TMatrix3D.h"
#include "Ratcliff_Triangulate2D.h"
#include "SimplificationByGarlan.h"
#include "rg_GeoFunc.h"


#include <set>
#include <algorithm>
using namespace std;


PolygonMeshModel::PolygonMeshModel()
{
}


    
PolygonMeshModel::PolygonMeshModel(const PolygonMeshModel& polygonMeshModel)
{
    duplicate( polygonMeshModel );
}



PolygonMeshModel::~PolygonMeshModel()
{
    clear();
}


    
PolygonMeshModel& PolygonMeshModel::operator =(const PolygonMeshModel& polygonMeshModel)
{
    if ( this != &polygonMeshModel ) {
        duplicate( polygonMeshModel );
    }

    return *this;
}



SurfaceModel::Body*     PolygonMeshModel::create_body()
{
    int ID = 0;
    if (!m_bodies.empty()) {
        ID = m_bodies.back()->getID() + 1;
    }
    m_bodies.push_back(new PolygonMeshModel::Body(ID));
    ++m_number_of_bodies;

    return m_bodies.back();
}



SurfaceModel::Shell*    PolygonMeshModel::create_shell()
{
    int ID = 0;
    if (!m_shells.empty()) {
        ID = m_shells.back()->getID() + 1;
    }
    m_shells.push_back(new PolygonMeshModel::Shell(ID));
    ++m_number_of_shells;

    return m_shells.back();
}



SurfaceModel::Face*     PolygonMeshModel::create_face()
{
    int ID = 0;
    if (!m_faces.empty()) {
        ID = m_faces.back()->getID() + 1;
    }
    m_faces.push_back(new PolygonMeshModel::Face(ID));
    ++m_number_of_faces;

    return m_faces.back();
}



SurfaceModel::Loop*     PolygonMeshModel::create_loop()
{
    int ID = 0;
    if (!m_loops.empty()) {
        ID = m_loops.back()->getID() + 1;
    }
    m_loops.push_back(new PolygonMeshModel::Loop(ID));
    ++m_number_of_loops;

    return m_loops.back();
}



SurfaceModel::Edge*     PolygonMeshModel::create_edge()
{
    int ID = 0;
    if (!m_edges.empty()) {
        ID = m_edges.back()->getID() + 1;
    }
    m_edges.push_back(new PolygonMeshModel::Edge(ID));
    ++m_number_of_edges;

    return m_edges.back();
}



SurfaceModel::Vertex*   PolygonMeshModel::create_vertex()
{
    int ID = 0;
    if (!m_vertices.empty()) {
        ID = m_vertices.back()->getID() + 1;
    }
    m_vertices.push_back(new PolygonMeshModel::Vertex(ID));
    ++m_number_of_vertices;

    return m_vertices.back();
}

    

PolygonMeshModel::Body*     PolygonMeshModel::create(const PolygonMeshModel::Body&   body)
{
    PolygonMeshModel::Body* newBody = new PolygonMeshModel::Body( body );
    SurfaceModel::concatenate_body( newBody );

    return newBody;
}



PolygonMeshModel::Shell*    PolygonMeshModel::create(const PolygonMeshModel::Shell&  shell)
{
    PolygonMeshModel::Shell* newShell = new PolygonMeshModel::Shell( shell );
    SurfaceModel::concatenate_shell( newShell );

    return newShell;
}



PolygonMeshModel::Face*     PolygonMeshModel::create(const PolygonMeshModel::Face&   face)
{
    PolygonMeshModel::Face* newFace = new PolygonMeshModel::Face( face );
    SurfaceModel::concatenate_face( newFace );

    return newFace;
}



PolygonMeshModel::Loop*     PolygonMeshModel::create(const PolygonMeshModel::Loop&   loop)
{
    PolygonMeshModel::Loop* newLoop = new PolygonMeshModel::Loop( loop );
    SurfaceModel::concatenate_loop( newLoop );

    return newLoop;
}



PolygonMeshModel::Edge*     PolygonMeshModel::create(const PolygonMeshModel::Edge&   edge)
{
    PolygonMeshModel::Edge* newEdge = new PolygonMeshModel::Edge( edge );
    SurfaceModel::concatenate_edge( newEdge );

    return newEdge;
}



PolygonMeshModel::Vertex*   PolygonMeshModel::create(const PolygonMeshModel::Vertex& vertex)
{
    PolygonMeshModel::Vertex* newVertex = new PolygonMeshModel::Vertex( vertex );
    SurfaceModel::concatenate_vertex( newVertex );

    return newVertex;
}



PolygonMeshModel::Vertex*  PolygonMeshModel::find(const rg_Point3D& point) const
{
    Vertex* correspondingVertex = rg_NULL;

    const list< SurfaceModel::Vertex* >& all_vertices = SurfaceModel::get_all_vertices();

    list<SurfaceModel::Vertex*>::const_iterator i_vtx;
    for ( i_vtx=all_vertices.begin(); i_vtx!=all_vertices.end(); ++i_vtx ) {
        PolygonMeshModel::Vertex* currVtx = (PolygonMeshModel::Vertex*)(*i_vtx);

        if ( currVtx->coordinate().isEqual( point, rg_SYSTEM_RES ) ) {
            correspondingVertex = currVtx;
            break;
        }
    }

    return correspondingVertex;
}



AxisAlignedBox          PolygonMeshModel::computeBoundingBox() const
{
    AxisAlignedBox boundingBox;

    for ( list<SurfaceModel::Vertex*>::const_iterator i_vtx=m_vertices.begin(); i_vtx!=m_vertices.end(); ++i_vtx ) {
        PolygonMeshModel::Vertex* currVtx = (PolygonMeshModel::Vertex*)(*i_vtx);

        rg_Point3D coordinate = currVtx->coordinate();

        boundingBox.update( coordinate );
    }

    return boundingBox;
}


unsigned int    PolygonMeshModel::all_bodies(   list<PolygonMeshModel::Body*>&    all_bodies) const
{
    all_bodies.clear();

    unsigned int numBodies = 0;
    const list< SurfaceModel::Body* >& bodies = SurfaceModel::get_all_bodies();

    list< SurfaceModel::Body* >::const_iterator i_body;
    for ( i_body=bodies.begin(); i_body!=bodies.end(); ++i_body ) {
        all_bodies.push_back( (PolygonMeshModel::Body*)(*i_body) );
        ++numBodies;
    }

    return numBodies;
}



unsigned int    PolygonMeshModel::all_shells(   list<PolygonMeshModel::Shell*>&   all_shells) const
{
    all_shells.clear();

    unsigned int numShells = 0;
    const list< SurfaceModel::Shell* >& shells = SurfaceModel::get_all_shells();

    list< SurfaceModel::Shell* >::const_iterator i_shell;
    for ( i_shell=shells.begin(); i_shell!=shells.end(); ++i_shell ) {
        all_shells.push_back( (PolygonMeshModel::Shell*)(*i_shell) );
        ++numShells;
    }

    return numShells;
}



unsigned int    PolygonMeshModel::all_faces(    list<PolygonMeshModel::Face*>&    all_faces) const
{
    all_faces.clear();

    unsigned int numFaces = 0;
    const list< SurfaceModel::Face* >& faces = SurfaceModel::get_all_faces();

    list< SurfaceModel::Face* >::const_iterator i_face;
    for ( i_face=faces.begin(); i_face!=faces.end(); ++i_face ) {
        all_faces.push_back( (PolygonMeshModel::Face*)(*i_face) );
        ++numFaces;
    }

    return numFaces;
}



unsigned int    PolygonMeshModel::all_loops(    list<PolygonMeshModel::Loop*>&    all_loops) const
{
    all_loops.clear();

    unsigned int numLoops = 0;
    const list< SurfaceModel::Loop* >& loops = SurfaceModel::get_all_loops();

    list< SurfaceModel::Loop* >::const_iterator i_loop;
    for ( i_loop=loops.begin(); i_loop!=loops.end(); ++i_loop ) {
        all_loops.push_back( (PolygonMeshModel::Loop*)(*i_loop) );
        ++numLoops;
    }

    return numLoops;
}



unsigned int    PolygonMeshModel::all_edges(    list<PolygonMeshModel::Edge*>&    all_edges) const
{
    all_edges.clear();

    unsigned int numEdges = 0;
    const list< SurfaceModel::Edge* >& edges = SurfaceModel::get_all_edges();

    list< SurfaceModel::Edge* >::const_iterator i_edge;
    for ( i_edge=edges.begin(); i_edge!=edges.end(); ++i_edge ) {
        all_edges.push_back( (PolygonMeshModel::Edge*)(*i_edge) );
        ++numEdges;
    }

    return numEdges;
}



unsigned int    PolygonMeshModel::all_vertices( list<PolygonMeshModel::Vertex*>&  all_vertices) const
{
    all_vertices.clear();

    unsigned int numVertices = 0;
    const list< SurfaceModel::Vertex* >& vertices = SurfaceModel::get_all_vertices();

    list< SurfaceModel::Vertex* >::const_iterator i_vtx;
    for ( i_vtx=vertices.begin(); i_vtx!=vertices.end(); ++i_vtx ) {
        all_vertices.push_back( (PolygonMeshModel::Vertex*)(*i_vtx) );
        ++numVertices;
    }

    return numVertices;
}

    

bool             PolygonMeshModel::divide_edge_by_inserting_vertex( PolygonMeshModel::Edge* edge, PolygonMeshModel::Vertex* splitVertex, PolygonMeshModel::Edge*& outputEdge1, PolygonMeshModel::Edge*& outputEdge2)
{
    SurfaceModel::Edge* outputEdge1SM  = (SurfaceModel::Edge*)outputEdge1; 
    SurfaceModel::Edge* outputEdge2SM  = (SurfaceModel::Edge*)outputEdge2; 

    bool result = SurfaceModel::divide_edge_by_inserting_vertex( edge, splitVertex, outputEdge1SM, outputEdge2SM );

    if ( result ) {
        outputEdge1  = (PolygonMeshModel::Edge*)outputEdge1SM;
        outputEdge2  = (PolygonMeshModel::Edge*)outputEdge2SM;
    }

    return result;
}


    
bool             PolygonMeshModel::divide_face_into_triangular_faces_by_inserting_vertex( 
                                        PolygonMeshModel::Face* face, PolygonMeshModel::Vertex* vertex, 
                                        list<PolygonMeshModel::Edge*>& newEdges, list<PolygonMeshModel::Face*>& newFaces)
{
    list<SurfaceModel::Edge*> newEdgesSM;
    list<SurfaceModel::Face*> newFacesSM;

    bool result = SurfaceModel::divide_face_into_triangular_faces_by_inserting_vertex( face, vertex, newEdgesSM, newFacesSM);

    for ( list<SurfaceModel::Edge*>::iterator i_edge=newEdgesSM.begin(); i_edge!=newEdgesSM.end(); ++i_edge ) {
        newEdges.push_back( (PolygonMeshModel::Edge*) *i_edge );
    }
    for ( list<SurfaceModel::Face*>::iterator i_face=newFacesSM.begin(); i_face!=newFacesSM.end(); ++i_face ) {
        PolygonMeshModel::Face* currFace = (PolygonMeshModel::Face*)(*i_face);

        vector<PolygonMeshModel::Vertex*> vertex;
        currFace->find_bounding_vertices( vertex );

        Plane facePlane( vertex[0]->coordinate(), vertex[1]->coordinate(), vertex[2]->coordinate() );

        currFace->set_normal( facePlane.getNormal() );

        newFaces.push_back( (PolygonMeshModel::Face*) *i_face );
    }

    return result;
}



bool             PolygonMeshModel::divide_face_by_inserting_edge( 
                                        PolygonMeshModel::Face* face, PolygonMeshModel::Vertex* vtx1OfEdge, PolygonMeshModel::Vertex* vtx2OfEdge, 
                                        PolygonMeshModel::Edge*& newEdge, PolygonMeshModel::Face*& newFace1, PolygonMeshModel::Face*& newFace2)
{
    SurfaceModel::Edge* newEdgeSM  = (SurfaceModel::Edge*)newEdge; 
    SurfaceModel::Face* newFace1SM = (SurfaceModel::Face*)newFace1; 
    SurfaceModel::Face* newFace2SM = (SurfaceModel::Face*)newFace2; 

    bool result = SurfaceModel::divide_face_by_inserting_edge(face, vtx1OfEdge, vtx2OfEdge, newEdgeSM, newFace1SM, newFace2SM);

    if ( result ) {
        newEdge  = (PolygonMeshModel::Edge*)newEdgeSM;
        newFace1 = (PolygonMeshModel::Face*)newFace1SM;
        newFace2 = (PolygonMeshModel::Face*)newFace2SM;

        newFace2->set_normal( newFace1->normal() );
        newFace2->set_face_group_ID( newFace1->face_group_ID() );
    }

    return result;
}


    
bool PolygonMeshModel::merge_two_faces_by_removing_edge( PolygonMeshModel::Face* face1, PolygonMeshModel::Face* face2 )
{
    rg_Point3D normal[2] = { face1->normal(), face2->normal() };
    bool areTwoFacesMerged = SurfaceModel::merge_two_faces_by_removing_edge( face1, face2 );

    if ( areTwoFacesMerged ) {
        face1->set_normal( (normal[0] + normal[1])/2.0 );
    }

    return areTwoFacesMerged;
}



//PolygonMeshModel::ResultOfEdgeCollapse PolygonMeshModel::collapse_edge( PolygonMeshModel::Edge* edgeToCollapse, const rg_Point3D& vertexPosition )
void PolygonMeshModel::collapse_edge( PolygonMeshModel::Edge* edgeToCollapse, const rg_Point3D& vertexPosition )
{
    //if ( edgeToCollapse->is_on_boundary() ) {
    //    return EDGE_IS_ON_BOUNDARY;
    //}

    //if ( does_nonmanifold_happen_by_edge_collapse( edgeToCollapse ) ) {
    //    return NONMANIFOLD_HAPPEN;
    //}

    //if ( does_folding_face_happen_by_edge_collapse( edgeToCollapse, vertexPosition ) ) {
    //    return FOLDING_FACE_HAPPEN;
    //}



    ////////////////////////////////////////////////////////////////////////
    //
    //  Collapse edge
    //

    Vertex* start_vertex = edgeToCollapse->start_vertex();
    Vertex* end_vertex   = edgeToCollapse->end_vertex();

    Face*   right_face   = edgeToCollapse->right_face();
    Face*   left_face    = edgeToCollapse->left_face();
    Edge*   right_hand   = edgeToCollapse->right_hand_edge();
    Edge*   left_hand    = edgeToCollapse->left_hand_edge();
    Edge*   right_leg    = edgeToCollapse->right_leg_edge();
    Edge*   left_leg     = edgeToCollapse->left_leg_edge();


    ///////////////////////////////////////////////////////////
    //  fix topology caused by edge-collapse

    //  set topology of vertices influenced by edge-collapse
    //      end vertex, right vertex, left vertex
    Vertex* right_vertex = right_hand->opposite_vertex( end_vertex );
    Vertex* left_vertex  = left_hand->opposite_vertex(  end_vertex );
    {
        end_vertex->set_first_edge(   right_hand );
        right_vertex->set_first_edge( right_hand );
        left_vertex->set_first_edge(  left_hand  );
    }


    //  set topology of faces influenced by edge-collapse
    //      right leg face, left leg face
    Face* right_leg_face = right_leg->opposite_face( right_face );
    Face* left_leg_face  = left_leg->opposite_face(  left_face  );
    {
        right_leg_face->set_first_edge( right_hand );
        left_leg_face->set_first_edge(  left_hand  );
        right_leg_face->exterior_loop()->set_first_edge( right_hand );
        left_leg_face->exterior_loop()->set_first_edge(  left_hand  );
    }


    //  set topology of edges influenced by edge-collapse
    //      right-hand, left-hand, edges of right-leg and left-leg
    Edge*   edge_for_right_leg[2] = {rg_NULL, rg_NULL};
    Edge*   edge_for_left_leg[2]  = {rg_NULL, rg_NULL};
    if ( start_vertex->is_start_vertex_of( right_leg ) ) {
        edge_for_right_leg[0] = right_leg->right_hand_edge();
        edge_for_right_leg[1] = right_leg->right_leg_edge();
    }
    else {
        edge_for_right_leg[0] = right_leg->left_leg_edge();
        edge_for_right_leg[1] = right_leg->left_hand_edge();
    }

    if ( start_vertex->is_start_vertex_of( left_leg ) ) {
        edge_for_left_leg[0] = left_leg->left_hand_edge();
        edge_for_left_leg[1] = left_leg->left_leg_edge();
    }
    else {
        edge_for_left_leg[0] = left_leg->right_leg_edge();
        edge_for_left_leg[1] = left_leg->right_hand_edge();
    }


    if ( edge_for_right_leg[0] == edge_for_left_leg[0] ) {
        //  connect righ-hand with edge_for_right_leg[0]
        if ( right_vertex->is_end_vertex_of( right_hand ) ) {
            right_hand->connect_right_hand_edge( edge_for_right_leg[0] );
            right_hand->set_right_face( right_leg_face );
        }
        else {
            right_hand->connect_left_leg_edge(   edge_for_right_leg[0] );
            right_hand->set_left_face( right_leg_face );
        }

        //  connect left-hand with edge_for_left_leg[0]
        if ( left_vertex->is_end_vertex_of( left_hand ) ) {
            left_hand->connect_left_hand_edge( edge_for_left_leg[0] );
            left_hand->set_left_face( left_leg_face );
        }
        else {
            left_hand->connect_right_leg_edge( edge_for_left_leg[0] );
            left_hand->set_right_face( left_leg_face );
        }

        //  connect righ-hand with left-hand
        if ( end_vertex->is_end_vertex_of( right_hand ) ) {
            right_hand->connect_left_hand_edge(  left_hand );
        }
        else {
            right_hand->connect_right_leg_edge(  left_hand );
        }
    }
    else {
        list<Edge*> edgesIncidentToStartVtx;
        start_vertex->find_incident_edges( edgesIncidentToStartVtx );
        for ( list<Edge*>::iterator i_edge=edgesIncidentToStartVtx.begin(); i_edge!=edgesIncidentToStartVtx.end(); ++i_edge ) {
            Edge* edge = (*i_edge);
            if ( start_vertex->is_end_vertex_of( edge ) ) {
                edge->set_end_vertex( end_vertex );
            }
            else {
                edge->set_start_vertex( end_vertex );
            }
        }

        //  connect righ-hand with edges for right leg
        if ( right_vertex->is_end_vertex_of( right_hand ) ) {
            right_hand->connect_right_hand_edge( edge_for_right_leg[0] );
            right_hand->connect_right_leg_edge(  edge_for_right_leg[1] );
            right_hand->set_right_face( right_leg_face );
        }
        else {
            right_hand->connect_left_leg_edge(   edge_for_right_leg[0] );
            right_hand->connect_left_hand_edge(  edge_for_right_leg[1] );
            right_hand->set_left_face( right_leg_face );
        }


        //  connect left-hand with edge_for_left_leg[0]
        if ( left_vertex->is_end_vertex_of( left_hand ) ) {
            left_hand->connect_left_hand_edge(  edge_for_left_leg[0] );
            left_hand->connect_left_leg_edge(   edge_for_left_leg[1] );
            left_hand->set_left_face(           left_leg_face );
        }
        else {
            left_hand->connect_right_leg_edge(  edge_for_left_leg[0] );
            left_hand->connect_right_hand_edge( edge_for_left_leg[1] );
            left_hand->set_right_face(          left_leg_face );
        }
    }


    Shell* shell = right_face->shell();
    shell->remove_face( right_face );
    shell->remove_face( left_face  );


    ///////////////////////////////////////////////////////////
    //  set coordinate of end vertex
    end_vertex->set_coordinate( vertexPosition );

    list<PolygonMeshModel::Face*> facesIncidentToEnd;
    end_vertex->find_incident_faces( facesIncidentToEnd );
    for ( list<PolygonMeshModel::Face*>::iterator i_face=facesIncidentToEnd.begin(); i_face!=facesIncidentToEnd.end(); ++i_face ) {
        (*i_face)->compute_normal();
    }

    ///////////////////////////////////////////////////////////
    //  remove  start vertex, right-leg edge, left-leg edge, 
    //          right face, left face.
    start_vertex->will_be_deleted(true);
    right_leg->will_be_deleted(true);
    left_leg->will_be_deleted(true);
    edgeToCollapse->will_be_deleted(true);

    right_face->exterior_loop()->will_be_deleted(true);
    left_face->exterior_loop()->will_be_deleted(true);
    right_face->will_be_deleted(true);
    left_face->will_be_deleted(true);

    //remove_vertex( start_vertex     );
    //remove_edge(   right_leg        );
    //remove_edge(   left_leg         );
    //remove_edge(   edgeToCollapse   );

    //remove_loop(   right_face->exterior_loop() );
    //remove_loop(   left_face->exterior_loop() );
    //remove_face(   right_face );
    //remove_face(   left_face  );

        
    
    //return EDGE_IS_COLLAPSED;
}



void PolygonMeshModel::run_CatmallClark_subdivision_once()
{
}

    
void PolygonMeshModel::run_Loop_subdivision(const int& iteration)
{
    for ( int i=0; i<iteration; ++i ) {
        run_Loop_subdivision_once();
    }
}

void PolygonMeshModel::run_Loop_subdivision_once()
{
    list<Shell*> shellsInModel;
    all_shells( shellsInModel );

    for ( list<Shell*>::iterator i_shell=shellsInModel.begin(); i_shell!=shellsInModel.end(); ++i_shell ) {
        Shell* currShell = *i_shell;

        run_Loop_subdivision_once_for_single_shell( currShell );
    }
}



void PolygonMeshModel::triangulate()
{
    list<Face*> facesInModel;
    all_faces( facesInModel );


    for ( list<Face*>::iterator i_face=facesInModel.begin(); i_face!=facesInModel.end(); ++i_face ) {
        triangulate_single_face( *i_face );
    }
}



void  PolygonMeshModel::simplify(const double& limitOfQEM)
{
    SimplificationByGarlan simplifier;
    simplifier.run_simplification( this, limitOfQEM );
}

    

bool  PolygonMeshModel::triangulate_single_face( PolygonMeshModel::Face* face )
{

    list<Vertex*> boundingVertices;
    int numBoundingVertices = face->find_bounding_vertices( boundingVertices );

    if ( numBoundingVertices <= 3 ) {
        return false;
    }



    rg_Point3D originToTranslate = boundingVertices.front()->coordinate();
    rg_Point3D faceNormal        = face->normal();

	rg_TMatrix3D trMatrix;
    {
	    trMatrix.translate(-originToTranslate);
        rg_FLAG bReverse = rg_EQ( faceNormal.innerProduct( rg_Point3D(0.0, 0.0, 1.0) ), -1);
        if( !bReverse )
	        trMatrix.rotate(faceNormal, rg_Point3D(0.0, 0.0, 1.0));
        else
            trMatrix.rotateY(rg_PI);
    }


    map<int, Vertex*> vertexIDMap;
    vector<Ratcliff::Vector2d> points2D( boundingVertices.size() );
    int i=0;
    list<Vertex*>::iterator i_vtx;
    for ( i_vtx=boundingVertices.begin(); i_vtx!=boundingVertices.end(); ++i_vtx, ++i ) {
        Vertex* currVtx = *i_vtx;

        vertexIDMap.insert( make_pair( currVtx->getID(), currVtx ) );

        rg_Point3D  point = currVtx->coordinate();
        rg_Point3D  transformedPoint = trMatrix*point;
        points2D[i].set( currVtx->getID(), transformedPoint.getX(), transformedPoint.getY() );
    }


    list<Ratcliff::Triangle2d> triangles;
    bool  isTriangulationCompleted = Ratcliff::Triangulate2d::triangulate( points2D, triangles );
    if ( !isTriangulationCompleted ) {
        return false;
    }


    while ( !triangles.empty() ) {
        Ratcliff::Triangle2d currTriangle = triangles.front();
        triangles.pop_front();

        int     vertexID[3] = { currTriangle.point(0).ID(), 
                                currTriangle.point(1).ID(), 
                                currTriangle.point(2).ID() };
        Vertex* vertex[3]   = { vertexIDMap.find( vertexID[0] )->second, 
                                vertexIDMap.find( vertexID[1] )->second, 
                                vertexIDMap.find( vertexID[2] )->second };

        Edge*   edge[3]     = { face->find_edge( vertex[0], vertex[1] ), 
                                face->find_edge( vertex[1], vertex[2] ), 
                                face->find_edge( vertex[2], vertex[0] ) };

        Edge* newEdge = rg_NULL;
        Face* newFace[2] = {rg_NULL, rg_NULL};

        if ( edge[0]!=rg_NULL && edge[1]!=rg_NULL && edge[2]==rg_NULL ) {
            divide_face_by_inserting_edge( face, vertex[0], vertex[2], newEdge, newFace[0], newFace[1] );
        }
        else if ( edge[0]!=rg_NULL && edge[1]==rg_NULL && edge[2]!=rg_NULL ) {
            divide_face_by_inserting_edge( face, vertex[2], vertex[1], newEdge, newFace[0], newFace[1] );
        }
        else if ( edge[0]==rg_NULL && edge[1]!=rg_NULL && edge[2]!=rg_NULL ) {
            divide_face_by_inserting_edge( face, vertex[1], vertex[0], newEdge, newFace[0], newFace[1] );
        }
        else if ( edge[0]!=rg_NULL && edge[1]!=rg_NULL && edge[2]!=rg_NULL ) {
            break;
        }
        else {
            triangles.push_back( currTriangle );
        }
    } 

        
    return true;
}



void PolygonMeshModel::duplicate(const PolygonMeshModel& polygonMeshModel)
{
    clear();

    map<SurfaceModel::Body*,   SurfaceModel::Body*>   bodyMap;
    map<SurfaceModel::Shell*,  SurfaceModel::Shell*>  shellMap;
    map<SurfaceModel::Face*,   SurfaceModel::Face*>   faceMap;
    map<SurfaceModel::Loop*,   SurfaceModel::Loop*>   loopMap;
    map<SurfaceModel::Edge*,   SurfaceModel::Edge*>   edgeMap;
    map<SurfaceModel::Vertex*, SurfaceModel::Vertex*> vertexMap;

    SurfaceModel::make_entity_map_for_duplication(    polygonMeshModel, bodyMap, shellMap, faceMap, loopMap, edgeMap, vertexMap );

    SurfaceModel::duplicate_topology(                 polygonMeshModel, bodyMap, shellMap, faceMap, loopMap, edgeMap, vertexMap );

    PolygonMeshModel::duplicate_property(             polygonMeshModel, bodyMap, shellMap, faceMap, loopMap, edgeMap, vertexMap );
}
        


void PolygonMeshModel::duplicate_property( 
                                 const PolygonMeshModel&                                   polygonMeshModel, 
                                 const map<SurfaceModel::Body*,   SurfaceModel::Body*>&    bodyMap,
                                 const map<SurfaceModel::Shell*,  SurfaceModel::Shell*>&   shellMap,
                                 const map<SurfaceModel::Face*,   SurfaceModel::Face*>&    faceMap,
                                 const map<SurfaceModel::Loop*,   SurfaceModel::Loop*>&    loopMap,
                                 const map<SurfaceModel::Edge*,   SurfaceModel::Edge*>&    edgeMap,
                                 const map<SurfaceModel::Vertex*, SurfaceModel::Vertex*>&  vertexMap)
{
    map<SurfaceModel::Face*,   SurfaceModel::Face*>::const_iterator i_face;
    for ( i_face=faceMap.begin(); i_face!=faceMap.end(); ++i_face ) {
        PolygonMeshModel::Face* original   = (PolygonMeshModel::Face*) (i_face->first);
        PolygonMeshModel::Face* duplicated = (PolygonMeshModel::Face*) (i_face->second);

        if ( original == rg_NULL ) {
            continue;
        }

        duplicated->set_normal( original->normal() );
    }


    map<SurfaceModel::Vertex*,   SurfaceModel::Vertex*>::const_iterator i_vtx;
    for ( i_vtx=vertexMap.begin(); i_vtx!=vertexMap.end(); ++i_vtx ) {
        PolygonMeshModel::Vertex* original   = (PolygonMeshModel::Vertex*) (i_vtx->first);
        PolygonMeshModel::Vertex* duplicated = (PolygonMeshModel::Vertex*) (i_vtx->second);

        if ( original == rg_NULL ) {
            continue;
        }

        duplicated->set_coordinate( original->coordinate() );
    }
}




void  PolygonMeshModel::run_Loop_subdivision_once_for_single_shell(PolygonMeshModel::Shell* shell)
{
    map<PolygonMeshModel::Vertex*, PolygonMeshModel::Edge*> oddVertexMap;
    map<PolygonMeshModel::Vertex*, rg_Point3D>              evenVertexMap;
    compute_odd_and_even_vertices( shell, oddVertexMap, evenVertexMap );


    //  divide edges in shell by inserting odd vertices
    map<PolygonMeshModel::Vertex*, PolygonMeshModel::Edge*>::iterator i_odd;
    for ( i_odd=oddVertexMap.begin(); i_odd!=oddVertexMap.end(); ++i_odd ) {
        PolygonMeshModel::Vertex* oddVertex = i_odd->first;
        PolygonMeshModel::Edge*   edgeToInsertOdd = i_odd->second;

        PolygonMeshModel::Edge*   newEdge[2] = {rg_NULL, rg_NULL};
        divide_edge_by_inserting_vertex( edgeToInsertOdd, oddVertex, newEdge[0], newEdge[1]);
    }


    //  divide faces in shell by inserting edges to connect two odd vertices
    list<PolygonMeshModel::Face*> facesInShell;
    shell->find_bounding_faces( facesInShell );

    list<PolygonMeshModel::Face*>::iterator i_face;
    for ( i_face=facesInShell.begin(); i_face!=facesInShell.end(); ++i_face ) {
        PolygonMeshModel::Face* currFace = *i_face;

        list<PolygonMeshModel::Vertex*> boundingVertices;
        currFace->find_bounding_vertices( boundingVertices );

        int pos = 0;
        vector<PolygonMeshModel::Vertex*> oddVertexInFace(4, rg_NULL);
        for ( list<PolygonMeshModel::Vertex*>::iterator i_vtx=boundingVertices.begin(); i_vtx!=boundingVertices.end(); ++i_vtx ) {
            PolygonMeshModel::Vertex* currVtx = *i_vtx;
            if ( oddVertexMap.find( currVtx ) != oddVertexMap.end() ) {
                oddVertexInFace[pos] = currVtx;
                ++pos;
            }
        }
        oddVertexInFace[3] = oddVertexInFace[0];

        PolygonMeshModel::Face* faceToDivide = currFace;
        for ( pos=0; pos<3; ++pos ) {
            PolygonMeshModel::Vertex* start = oddVertexInFace[pos];
            PolygonMeshModel::Vertex* end   = oddVertexInFace[pos+1];

            PolygonMeshModel::Edge*   newEdge    = rg_NULL;
            PolygonMeshModel::Face*   newFace[2] = {rg_NULL, rg_NULL};

            divide_face_by_inserting_edge( faceToDivide, start, end, newEdge, newFace[0], newFace[1] );

            faceToDivide = newFace[0];
        }
    }


    //  modify the position of even vertices
    map<PolygonMeshModel::Vertex*, rg_Point3D>::iterator i_even;
    for ( i_even=evenVertexMap.begin(); i_even!=evenVertexMap.end(); ++i_even ) {
        PolygonMeshModel::Vertex* evenVertex = i_even->first;
        evenVertex->set_coordinate( i_even->second );
    }
}



void  PolygonMeshModel::compute_odd_and_even_vertices(
                            PolygonMeshModel::Shell* shell,
                            map<PolygonMeshModel::Vertex*, PolygonMeshModel::Edge*>& oddVertexMap,
                            map<PolygonMeshModel::Vertex*, rg_Point3D>&              evenVertexMap )
{
    //  compute position of odd vertices
    int vertexID = number_of_vertices();

    list<PolygonMeshModel::Edge*> edgesInShell;
    shell->find_bounding_edges( edgesInShell );

    list<PolygonMeshModel::Edge*>::iterator i_edge;
    for ( i_edge=edgesInShell.begin(); i_edge!=edgesInShell.end(); ++i_edge ) {
        PolygonMeshModel::Edge* currEdge = *i_edge;
        PolygonMeshModel::Vertex* startVtx = currEdge->start_vertex();
        PolygonMeshModel::Vertex* endVtx   = currEdge->end_vertex();

        rg_Point3D startPoint = startVtx->coordinate();
        rg_Point3D endPoint   = endVtx->coordinate();

        rg_Point3D positionOfOddVertex;
        if ( currEdge->is_on_boundary() ) {
            positionOfOddVertex = 0.5*(startPoint+endPoint);
        }
        else {
            rg_Point3D leftPoint  = currEdge->left_hand_edge()->opposite_vertex( endVtx )->coordinate();
            rg_Point3D rightPoint = currEdge->right_hand_edge()->opposite_vertex( endVtx )->coordinate();

            positionOfOddVertex = 0.375*(startPoint+endPoint) + 0.125*(leftPoint+rightPoint);
            //positionOfOddVertex = (3.0/8.0)*(startPoint+endPoint) + (1.0/8.0)*(leftPoint+rightPoint);
        }

        PolygonMeshModel::Vertex* oddVertex = create( PolygonMeshModel::Vertex( ++vertexID ) );
        oddVertex->set_coordinate( positionOfOddVertex );

        oddVertexMap.insert( make_pair(oddVertex, currEdge) );
    }


    //  compute position of even vertices
    list<PolygonMeshModel::Vertex*> verticesInShell;
    shell->find_bounding_vertices( verticesInShell );
    list<PolygonMeshModel::Vertex*>::iterator i_vtx;
    for ( i_vtx=verticesInShell.begin(); i_vtx!=verticesInShell.end(); ++i_vtx ) {
        PolygonMeshModel::Vertex* currVtx = *i_vtx;

        rg_Point3D positionOfEvenVertex;
        if ( currVtx->is_on_boundary() ) {
            list<PolygonMeshModel::Edge*> incidentEdges;
            currVtx->find_incident_edges( incidentEdges );

            rg_Point3D sumPointOfAdjVerticesOnBoundary(0.0, 0.0, 0.0);
            list<PolygonMeshModel::Edge*>::iterator i_incie;
            for ( i_incie=incidentEdges.begin(); i_incie!=incidentEdges.end(); ++i_incie ) {
                if ( (*i_incie)->is_on_boundary() ) {
                    sumPointOfAdjVerticesOnBoundary += (*i_incie)->opposite_vertex( currVtx )->coordinate();
                }
            }

            positionOfEvenVertex = 0.125*sumPointOfAdjVerticesOnBoundary + 0.75*currVtx->coordinate();
        }
        else {
            list<PolygonMeshModel::Vertex*> adjacentVertices;
            int numAdjVertices = currVtx->find_adjacent_vertices( adjacentVertices );

            rg_Point3D sumPointOfAdjVertices(0.0, 0.0, 0.0);
            list<PolygonMeshModel::Vertex*>::iterator i_adjvtx;
            for ( i_adjvtx=adjacentVertices.begin(); i_adjvtx!=adjacentVertices.end(); ++i_adjvtx ) {
                sumPointOfAdjVertices += (*i_adjvtx)->coordinate();
            }
        
            double beta = 0.0;
            if ( numAdjVertices == 3 ) {
                beta = 0.1875; 
                // beta = 3.0/16.0;
            }
            else {
                beta = 0.375/numAdjVertices; 
                // beta = 3.0/(8.0*numAdjVertices);
            }

            positionOfEvenVertex = beta*sumPointOfAdjVertices + (1.0-beta*numAdjVertices)*currVtx->coordinate();
        }

        evenVertexMap.insert( make_pair( currVtx, positionOfEvenVertex ) );
    }
}


    
bool    PolygonMeshModel::does_nonmanifold_happen_by_edge_collapse( Edge* edgeToCollapse ) const
{
    PolygonMeshModel::Vertex* startVtx = edgeToCollapse->start_vertex();
    PolygonMeshModel::Vertex* endVtx   = edgeToCollapse->end_vertex();

    list<PolygonMeshModel::Vertex*> adjacentVertices[2];
    startVtx->find_adjacent_vertices( adjacentVertices[0] );
    endVtx->find_adjacent_vertices(   adjacentVertices[1] );

    vector<PolygonMeshModel::Vertex*> vertices[2];
    vertices[0].insert( vertices[0].begin(), adjacentVertices[0].begin(), adjacentVertices[0].end() );
    vertices[1].insert( vertices[1].begin(), adjacentVertices[1].begin(), adjacentVertices[1].end() );

    sort( vertices[0].begin(), vertices[0].end() );
    sort( vertices[1].begin(), vertices[1].end() );

    int size = max( vertices[0].size(), vertices[1].size() );
    vector<PolygonMeshModel::Vertex*> common_adjacent_vertices(size);
    vector<PolygonMeshModel::Vertex*>::iterator result;
    result = set_intersection( vertices[0].begin(), vertices[0].end(), 
                               vertices[1].begin(), vertices[1].end(), 
                               common_adjacent_vertices.begin() );


    bool doesNonmanifoldHappenByEdgeCollapse = false;
    PolygonMeshModel::Vertex* rightVtx = edgeToCollapse->right_hand_edge()->opposite_vertex( endVtx );
    PolygonMeshModel::Vertex* leftVtx  = edgeToCollapse->left_hand_edge()->opposite_vertex( endVtx );

    vector<PolygonMeshModel::Vertex*>::iterator i_vtx;
    for ( i_vtx=common_adjacent_vertices.begin(); i_vtx!=result; ++i_vtx ) {
        PolygonMeshModel::Vertex* commonVertex = *i_vtx;

        if ( commonVertex!=rightVtx && commonVertex!=leftVtx ) {
            doesNonmanifoldHappenByEdgeCollapse = true;
            break;
        }
    }

    return doesNonmanifoldHappenByEdgeCollapse;
}


    
bool    PolygonMeshModel::does_folding_face_happen_by_edge_collapse( Edge* edgeToCollapse, const rg_Point3D& vertexPosition ) const
{
    double limitOfPlaneChange = rg_PI/4.0;

    PolygonMeshModel::Vertex* startVtx  = edgeToCollapse->start_vertex();
    PolygonMeshModel::Edge*   rightHand = edgeToCollapse->right_hand_edge();
    PolygonMeshModel::Edge*   leftHand  = edgeToCollapse->left_hand_edge();

    rg_Point3D startPoint = startVtx->coordinate();

    list<PolygonMeshModel::Edge*> edgesInShellOfVertex;
    startVtx->find_edges_in_shell( edgesInShellOfVertex );

    list<PolygonMeshModel::Edge*>::iterator i_edge;
    for ( i_edge=edgesInShellOfVertex.begin(); i_edge!=edgesInShellOfVertex.end(); ++i_edge ) {
        PolygonMeshModel::Edge* currEdge = *i_edge;

        if ( currEdge == rightHand || currEdge == leftHand ) {
            continue;
        }

        rg_Point3D point[2];
        if ( currEdge->left_hand_edge()->is_incident_to( startVtx ) ) {
            point[0] = currEdge->start_vertex()->coordinate();
            point[1] = currEdge->end_vertex()->coordinate();
        }
        else {
            point[1] = currEdge->start_vertex()->coordinate();
            point[0] = currEdge->end_vertex()->coordinate();
        }
        Plane planeBefore( point[0], point[1], startPoint );
        Plane planeAfter(  point[0], point[1], vertexPosition );
        double angleVariance = planeBefore.getNormal().angle( planeAfter.getNormal() );
        if ( angleVariance > limitOfPlaneChange ) {
            return true;
        }
    }



    PolygonMeshModel::Vertex* endVtx   = edgeToCollapse->end_vertex();
    PolygonMeshModel::Edge*   rightLeg = edgeToCollapse->right_leg_edge();
    PolygonMeshModel::Edge*   leftLeg  = edgeToCollapse->left_leg_edge();

    rg_Point3D endPoint = endVtx->coordinate();
    edgesInShellOfVertex.clear();
    edgeToCollapse->end_vertex()->find_edges_in_shell(   edgesInShellOfVertex );

    for ( i_edge=edgesInShellOfVertex.begin(); i_edge!=edgesInShellOfVertex.end(); ++i_edge ) {
        PolygonMeshModel::Edge* currEdge = *i_edge;

        if ( currEdge == rightLeg || currEdge == leftLeg ) {
            continue;
        }

        rg_Point3D point[2];
        if ( currEdge->left_hand_edge()->is_incident_to( endVtx ) ) {
            point[0] = currEdge->start_vertex()->coordinate();
            point[1] = currEdge->end_vertex()->coordinate();
        }
        else {
            point[1] = currEdge->start_vertex()->coordinate();
            point[0] = currEdge->end_vertex()->coordinate();
        }
        Plane planeBefore( point[0], point[1], endPoint );
        Plane planeAfter(  point[0], point[1], vertexPosition );
        double angleVariance = planeBefore.getNormal().angle( planeAfter.getNormal() );

        if ( angleVariance > limitOfPlaneChange ) {
            return true;
        }
    }


    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Body 
//
PolygonMeshModel::Body::Body()
: SurfaceModel::Body()
{
}



PolygonMeshModel::Body::Body(const rg_INT& ID)
: SurfaceModel::Body(ID)
{
}



PolygonMeshModel::Body::Body(const Body& body)
: SurfaceModel::Body( body )
{
}



PolygonMeshModel::Body::~Body()
{
}



PolygonMeshModel::Body& PolygonMeshModel::Body::operator =(const Body& body)
{
    if ( this != &body ) {
        SurfaceModel::Body::operator=(body);
    }

    return *this;
}



PolygonMeshModel::Shell*    PolygonMeshModel::Body::exterior_shell() const
{
    return (PolygonMeshModel::Shell*)SurfaceModel::Body::exterior_shell();
}



unsigned int                PolygonMeshModel::Body::interior_shell(list< PolygonMeshModel::Shell* >& interiorShells) const
{
    interiorShells.clear();
    unsigned int numInteriorShells = 0;
    const list< SurfaceModel::Shell* >& interior_shells = SurfaceModel::Body::interior_shell();

    list< SurfaceModel::Shell* >::const_iterator i_shell;
    for ( i_shell=interior_shells.begin(); i_shell!=interior_shells.end(); ++i_shell ) {
        interiorShells.push_back( (PolygonMeshModel::Shell*)(*i_shell) );
        ++numInteriorShells;
    }

    return numInteriorShells;
}


///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Shell 
//
PolygonMeshModel::Shell::Shell()
: SurfaceModel::Shell()
, m_visible(true)
{
}



PolygonMeshModel::Shell::Shell(const rg_INT& ID)
: SurfaceModel::Shell(ID)
, m_visible(true)
{
}



PolygonMeshModel::Shell::Shell(const Shell& shell)
: SurfaceModel::Shell(shell)
, m_visible(shell.m_visible)
{
}



PolygonMeshModel::Shell::~Shell()
{
}



PolygonMeshModel::Shell& PolygonMeshModel::Shell::operator =(const Shell& shell)
{
    if ( this != &shell ) {
        SurfaceModel::Shell::operator=( shell );
        m_visible = shell.m_visible;
    }

    return *this;
}



unsigned int                 PolygonMeshModel::Shell::find_bounding_faces( list<PolygonMeshModel::Face*>&   boundingFaces) const
{
    boundingFaces.clear();
    unsigned int numBoundingFaces = 0;

    const list< SurfaceModel::Face* >& boundingFacesSM = SurfaceModel::Shell::faces();
    
    list< SurfaceModel::Face* >::const_iterator i_face;
    for ( i_face=boundingFacesSM.begin(); i_face!=boundingFacesSM.end(); ++i_face ) {
        boundingFaces.push_back( (PolygonMeshModel::Face*)(*i_face) );
        ++numBoundingFaces;
    }

    return numBoundingFaces;
}



unsigned int                 PolygonMeshModel::Shell::find_bounding_edges( list<PolygonMeshModel::Edge*>&   boundingEdges) const
{
    boundingEdges.clear();

    list<SurfaceModel::Edge*>   boundingEdgesSM;
    SurfaceModel::Shell::find_bounding_edges( boundingEdgesSM );

    list< SurfaceModel::Edge* >::const_iterator i_edge;
    for ( i_edge=boundingEdgesSM.begin(); i_edge!=boundingEdgesSM.end(); ++i_edge ) {
        boundingEdges.push_back( (PolygonMeshModel::Edge*) (*i_edge) );
    }

    return boundingEdges.size();
}


    

unsigned int                 PolygonMeshModel::Shell::find_bounding_vertices( list<PolygonMeshModel::Vertex*>& boundingVertices) const
{
    boundingVertices.clear();

    list<SurfaceModel::Vertex*>   boundingVerticesSM;
    SurfaceModel::Shell::find_bounding_vertices( boundingVerticesSM );

    list< SurfaceModel::Vertex* >::const_iterator i_vtx;
    for ( i_vtx=boundingVerticesSM.begin(); i_vtx!=boundingVerticesSM.end(); ++i_vtx ) {
        boundingVertices.push_back( (PolygonMeshModel::Vertex*) (*i_vtx) );
    }

    return boundingVertices.size();
}


    
AxisAlignedBox              PolygonMeshModel::Shell::bounding_box() const
{
    AxisAlignedBox boundingBox;

    list<PolygonMeshModel::Vertex*> boundingVertices;
    find_bounding_vertices( boundingVertices );

    
    list<PolygonMeshModel::Vertex*>::iterator i_vtx;
    for ( i_vtx=boundingVertices.begin(); i_vtx!=boundingVertices.end(); ++i_vtx ) {
        boundingBox.update( (*i_vtx)->coordinate() );
    }

    return boundingBox;
}

    

double                      PolygonMeshModel::Shell::volume() const
{
    double shellVolume = 0.0;

    AxisAlignedBox boundingBox = bounding_box();
    boundingBox.enlarge(2.0);

    Plane  bottomPlane( rg_Point3D(0.0, 0.0, 1.0), boundingBox.getMinPoint() );

    list<PolygonMeshModel::Face*> facesInShell;
    find_bounding_faces( facesInShell );

    list<PolygonMeshModel::Face*>::iterator i_face;
    for ( i_face=facesInShell.begin(); i_face!=facesInShell.end(); ++i_face ) {
        PolygonMeshModel::Face* currFace = *i_face;

        vector<PolygonMeshModel::Vertex*> vertex;
        currFace->find_bounding_vertices( vertex );


        shellVolume += rg_GeoFunc::computeSignedVolumeOfCutAwayTriangularPrism(
                                        bottomPlane, 
                                        vertex[0]->coordinate(), vertex[1]->coordinate(), vertex[2]->coordinate());

    }

    return shellVolume;
}



double                      PolygonMeshModel::Shell::area() const
{
    double shellArea = 0.0;

    list<PolygonMeshModel::Face*> facesInShell;
    find_bounding_faces(facesInShell);

    list<PolygonMeshModel::Face*>::iterator i_face;
    for (i_face = facesInShell.begin(); i_face != facesInShell.end(); ++i_face) {
        PolygonMeshModel::Face* currFace = *i_face;

        vector<PolygonMeshModel::Vertex*> vertex;
        currFace->find_bounding_vertices(vertex);

        shellArea += rg_GeoFunc::computeTriangleArea(vertex[0]->coordinate(), vertex[1]->coordinate(), vertex[2]->coordinate());

    }

    return shellArea;
}


///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Face 
//
PolygonMeshModel::Face::Face()
: SurfaceModel::Face(),
  m_faceGroupID(-1)
{
}



PolygonMeshModel::Face::Face(const rg_INT& ID)
: SurfaceModel::Face( ID ),
  m_faceGroupID(-1)
{
}



PolygonMeshModel::Face::Face(const Face& face)
: SurfaceModel::Face( face ),
  m_normal( face.m_normal ),
  m_faceGroupID( face.m_faceGroupID )
{
}



PolygonMeshModel::Face::~Face()
{
}



PolygonMeshModel::Face& PolygonMeshModel::Face::operator =(const Face& face)
{
    if ( this != &face ) {
        SurfaceModel::Face::operator=( face );

        m_normal        = face.m_normal;
        m_faceGroupID   = face.m_faceGroupID;
    }

    return *this;
}



rg_Point3D  PolygonMeshModel::Face::normal() const
{
    return m_normal;
}



void        PolygonMeshModel::Face::set_normal(const rg_Point3D& normal)
{
    m_normal = normal;
}

    

PolygonMeshModel::Edge*       PolygonMeshModel::Face::find_edge(const PolygonMeshModel::Vertex* const vertex1, const PolygonMeshModel::Vertex* const vertex2 ) const
{
    return (PolygonMeshModel::Edge*) SurfaceModel::Face::find_edge( vertex1, vertex2 );
}


    
unsigned int  PolygonMeshModel::Face::find_bounding_vertices( list<PolygonMeshModel::Vertex*>& boundingVertices) const
{
    boundingVertices.clear();

    list<SurfaceModel::Vertex*> boundingVerticesSM;
    unsigned int numBoundingVertices = SurfaceModel::Face::find_bounding_vertices( boundingVerticesSM );

    list<SurfaceModel::Vertex*>::iterator i_vtx;
    for ( i_vtx=boundingVerticesSM.begin(); i_vtx!=boundingVerticesSM.end(); ++i_vtx ) {
        boundingVertices.push_back( (PolygonMeshModel::Vertex*)(*i_vtx) );
    }

    return numBoundingVertices;
}

    

unsigned int  PolygonMeshModel::Face::find_bounding_vertices( vector<PolygonMeshModel::Vertex*>& boundingVertices) const
{
    boundingVertices.clear();

    list<SurfaceModel::Vertex*> boundingVerticesSM;
    unsigned int numBoundingVertices = SurfaceModel::Face::find_bounding_vertices( boundingVerticesSM );

    list<SurfaceModel::Vertex*>::iterator i_vtx;
    for ( i_vtx=boundingVerticesSM.begin(); i_vtx!=boundingVerticesSM.end(); ++i_vtx ) {
        boundingVertices.push_back( (PolygonMeshModel::Vertex*)(*i_vtx) );
    }

    return numBoundingVertices;
}



unsigned int  PolygonMeshModel::Face::find_bounding_edges(    list<PolygonMeshModel::Edge*>&   boundingEdges) const
{
    boundingEdges.clear();

    list<SurfaceModel::Edge*> boundingEdgesSM;
    rg_INT numBoundingEdges = SurfaceModel::Face::find_bounding_edges( boundingEdgesSM );

    list<SurfaceModel::Edge*>::iterator i_edge;
    for ( i_edge=boundingEdgesSM.begin(); i_edge!=boundingEdgesSM.end(); ++i_edge ) {
        boundingEdges.push_back( (PolygonMeshModel::Edge*)(*i_edge) );
    }

    return numBoundingEdges;
}
    

    
unsigned int    PolygonMeshModel::Face::find_edges_sharing_with_adjacent_face( const PolygonMeshModel::Face* const face, list<PolygonMeshModel::Edge*>& sharingEdges) const
{
    list<SurfaceModel::Edge*> sharingEdgesSM;
    rg_INT numSharingEdges = SurfaceModel::Face::find_edges_sharing_with_adjacent_face( face, sharingEdgesSM );

    list<SurfaceModel::Edge*>::iterator i_edge;
    for ( i_edge=sharingEdgesSM.begin(); i_edge!=sharingEdgesSM.end(); ++i_edge ) {
        sharingEdges.push_back( (PolygonMeshModel::Edge*)(*i_edge) );
    }

    return numSharingEdges;
}



void            PolygonMeshModel::Face::compute_normal()
{
    vector<PolygonMeshModel::Vertex*> boundingVertex;
    find_bounding_vertices( boundingVertex );

    Plane facePlane( boundingVertex[0]->coordinate(), boundingVertex[1]->coordinate(), boundingVertex[2]->coordinate() );

    m_normal = facePlane.getNormal();
}



///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Loop 
//
PolygonMeshModel::Loop::Loop()
: SurfaceModel::Loop()
{
}



PolygonMeshModel::Loop::Loop(const rg_INT& ID)
: SurfaceModel::Loop( ID )
{
}



PolygonMeshModel::Loop::Loop(const Loop& loop)
: SurfaceModel::Loop( loop )
{
}



PolygonMeshModel::Loop::~Loop()
{
}



PolygonMeshModel::Loop& PolygonMeshModel::Loop::operator =(const Loop& loop)
{
    if ( this != &loop ) {
        SurfaceModel::Loop::operator=( loop );
    }

    return *this;
}





///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Edge 
//
PolygonMeshModel::Edge::Edge()
: SurfaceModel::Edge(),
  m_boundary_of_faceGroup(false)
{
}



PolygonMeshModel::Edge::Edge(const rg_INT& ID)
: SurfaceModel::Edge( ID ),
  m_boundary_of_faceGroup(false)
{
}



PolygonMeshModel::Edge::Edge(const Edge& edge)
: SurfaceModel::Edge( edge ),
  m_boundary_of_faceGroup( edge.m_boundary_of_faceGroup )
{
}



PolygonMeshModel::Edge::~Edge()
{
}



PolygonMeshModel::Edge& PolygonMeshModel::Edge::operator =(const Edge& edge)
{
    if ( this != &edge ) {
        SurfaceModel::Edge::operator=( edge );
        m_boundary_of_faceGroup = edge.m_boundary_of_faceGroup;
    }

    return *this;
}

    

PolygonMeshModel::Shell*  PolygonMeshModel::Edge::shell() const
{
    return (PolygonMeshModel::Shell*)(SurfaceModel::Edge::shell());
}


    
PolygonMeshModel::Vertex* PolygonMeshModel::Edge::opposite_vertex(const PolygonMeshModel::Vertex* const vertex) const
{
    return (PolygonMeshModel::Vertex*)( SurfaceModel::Edge::opposite_vertex( vertex ) );
}

PolygonMeshModel::Face*   PolygonMeshModel::Edge::opposite_face(const PolygonMeshModel::Face* const face) const
{
    return (PolygonMeshModel::Face*)( SurfaceModel::Edge::opposite_face( face ) );
}



rg_INT  PolygonMeshModel::Edge::find_edges_in_star(     list<PolygonMeshModel::Edge*>&   edgesInStar) const
{
    list<SurfaceModel::Edge*> edgesInStarSM;
    int numEdgesInStar = SurfaceModel::Edge::find_edges_in_star( edgesInStarSM );

    list<SurfaceModel::Edge*>::iterator i_edge;
    for ( i_edge=edgesInStarSM.begin(); i_edge!=edgesInStarSM.end(); ++i_edge ) {
        edgesInStar.push_back( (PolygonMeshModel::Edge*)(*i_edge) );
    }

    return numEdgesInStar;
}



rg_INT  PolygonMeshModel::Edge::find_incident_faces( list<PolygonMeshModel::Face*>& incidentFaces ) const
{
    list<SurfaceModel::Face*> incidentFacesSM;
    int numIncidentFaces = SurfaceModel::Edge::find_incident_faces( incidentFacesSM );

    list<SurfaceModel::Face*>::iterator i_face;
    for ( i_face=incidentFacesSM.begin(); i_face!=incidentFacesSM.end(); ++i_face ) {
        incidentFaces.push_back( (PolygonMeshModel::Face*)(*i_face) );
    }

    return numIncidentFaces;
}



///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Vertex
//
PolygonMeshModel::Vertex::Vertex()
: SurfaceModel::Vertex()
{
}



PolygonMeshModel::Vertex::Vertex(const rg_INT& ID)
: SurfaceModel::Vertex( ID )
{
}



PolygonMeshModel::Vertex::Vertex(const Vertex& vertex)
: SurfaceModel::Vertex( vertex ),
  m_coordinate( vertex.m_coordinate )
{
}



PolygonMeshModel::Vertex::~Vertex()
{
}



PolygonMeshModel::Vertex& PolygonMeshModel::Vertex::operator =(const Vertex& vertex)
{
    if ( this != &vertex ) {
        SurfaceModel::Vertex::operator=( vertex );

        m_coordinate = vertex.m_coordinate;
    }

    return *this;
}



rg_Point3D  PolygonMeshModel::Vertex::coordinate() const
{
    return m_coordinate;
}



void        PolygonMeshModel::Vertex::set_coordinate(const rg_Point3D& coordinate)
{
    m_coordinate = coordinate;
}


    
PolygonMeshModel::Shell*  PolygonMeshModel::Vertex::shell() const
{
    return (PolygonMeshModel::Shell*)(SurfaceModel::Vertex::shell());
}


        
    
rg_Point3D  PolygonMeshModel::Vertex::compute_normal() const
{
    list<PolygonMeshModel::Face*> incidentFaces;
    find_incident_faces( incidentFaces );

    rg_Point3D normal(0.0, 0.0, 0.0);
    for ( list<PolygonMeshModel::Face*>::iterator i_face=incidentFaces.begin(); i_face!=incidentFaces.end(); ++i_face ) {
        normal += (*i_face)->normal();
    }
    normal.normalize();

    return normal;
}



rg_INT      PolygonMeshModel::Vertex::find_adjacent_vertices( list<PolygonMeshModel::Vertex*>& adjacentVertices) const
{
    list<SurfaceModel::Vertex*> adjacentVerticesSM;
    int numAdjacentVertices = SurfaceModel::Vertex::find_adjacent_vertices( adjacentVerticesSM );

    list<SurfaceModel::Vertex*>::iterator i_vtx;
    for ( i_vtx=adjacentVerticesSM.begin(); i_vtx!=adjacentVerticesSM.end(); ++i_vtx ) {
        adjacentVertices.push_back( (PolygonMeshModel::Vertex*)(*i_vtx) );
    }

    return numAdjacentVertices;
}



rg_INT      PolygonMeshModel::Vertex::find_incident_edges(    list<PolygonMeshModel::Edge*>&   incidentEdges) const
{
    list<SurfaceModel::Edge*> incidentEdgesSM;
    int numIncidentEdges = SurfaceModel::Vertex::find_incident_edges( incidentEdgesSM );

    list<SurfaceModel::Edge*>::iterator i_edge;
    for ( i_edge=incidentEdgesSM.begin(); i_edge!=incidentEdgesSM.end(); ++i_edge ) {
        incidentEdges.push_back( (PolygonMeshModel::Edge*)(*i_edge) );
    }

    return numIncidentEdges;

}



rg_INT      PolygonMeshModel::Vertex::find_incident_faces(    list<PolygonMeshModel::Face*>&   incidentFaces) const
{
    list<SurfaceModel::Face*> incidentFacesSM;
    int numIncidentFaces = SurfaceModel::Vertex::find_incident_faces( incidentFacesSM );

    list<SurfaceModel::Face*>::iterator i_face;
    for ( i_face=incidentFacesSM.begin(); i_face!=incidentFacesSM.end(); ++i_face ) {
        incidentFaces.push_back( (PolygonMeshModel::Face*)(*i_face) );
    }

    return numIncidentFaces;
}


    
rg_INT      PolygonMeshModel::Vertex::find_edges_in_star(     list<PolygonMeshModel::Edge*>&   edgesInStar) const
{
    list<SurfaceModel::Edge*> edgesInStarSM;
    int numEdgesInStar = SurfaceModel::Vertex::find_edges_in_star( edgesInStarSM );

    list<SurfaceModel::Edge*>::iterator i_edge;
    for ( i_edge=edgesInStarSM.begin(); i_edge!=edgesInStarSM.end(); ++i_edge ) {
        edgesInStar.push_back( (PolygonMeshModel::Edge*)(*i_edge) );
    }

    return numEdgesInStar;
}



rg_INT      PolygonMeshModel::Vertex::find_edges_in_shell(    list<PolygonMeshModel::Edge*>&   edgesInShell) const
{
    list<SurfaceModel::Edge*> edgesInShellSM;
    int numEdgesInShell = SurfaceModel::Vertex::find_edges_in_shell( edgesInShellSM );

    list<SurfaceModel::Edge*>::iterator i_edge;
    for ( i_edge=edgesInShellSM.begin(); i_edge!=edgesInShellSM.end(); ++i_edge ) {
        edgesInShell.push_back( (PolygonMeshModel::Edge*)(*i_edge) );
    }

    return numEdgesInShell;
}



