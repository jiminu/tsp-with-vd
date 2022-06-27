#ifndef POLYGONMESHMODEL_H
#define POLYGONMESHMODEL_H

#include "SurfaceModel.h"
#include "rg_Point3D.h"

#include "AxisAlignedBox.h"


class PolygonMeshModel : public SurfaceModel
{
//  Face, Edge, Vertex는 winged-edge DS에 묶여있음.
public:
    class Body;
    class Shell;
    class Face;
    class Loop;
    class Edge;
    class Vertex;

    enum  ResultOfEdgeCollapse { EDGE_IS_COLLAPSED, NONMANIFOLD_HAPPEN, FOLDING_FACE_HAPPEN, EDGE_IS_ON_BOUNDARY };

public:
    PolygonMeshModel();
    PolygonMeshModel(const PolygonMeshModel& polygonMeshModel);
    ~PolygonMeshModel();

    PolygonMeshModel& operator =(const PolygonMeshModel& polygonMeshModel);

    virtual SurfaceModel::Body*         create_body();
    virtual SurfaceModel::Shell*        create_shell();
    virtual SurfaceModel::Face*         create_face();
    virtual SurfaceModel::Loop*         create_loop();
    virtual SurfaceModel::Edge*         create_edge();
    virtual SurfaceModel::Vertex*       create_vertex();

    PolygonMeshModel::Body*     create(const PolygonMeshModel::Body&   body);
    PolygonMeshModel::Shell*    create(const PolygonMeshModel::Shell&  shell);
    PolygonMeshModel::Face*     create(const PolygonMeshModel::Face&   face);
    PolygonMeshModel::Loop*     create(const PolygonMeshModel::Loop&   loop);
    PolygonMeshModel::Edge*     create(const PolygonMeshModel::Edge&   edge);
    PolygonMeshModel::Vertex*   create(const PolygonMeshModel::Vertex& vertex);


    Vertex*                     find(const rg_Point3D& point) const;

    AxisAlignedBox              computeBoundingBox() const;


    unsigned int                all_bodies(   list<Body*>&    all_bodies) const;
    unsigned int                all_shells(   list<Shell*>&   all_shells) const;
    unsigned int                all_faces(    list<Face*>&    all_faces) const;
    unsigned int                all_loops(    list<Loop*>&    all_loops) const;
    unsigned int                all_edges(    list<Edge*>&    all_edges) const;
    unsigned int                all_vertices( list<Vertex*>&  all_vertices) const;


    bool                        divide_edge_by_inserting_vertex( Edge* edge, Vertex* splitVertex, Edge*& outputEdge1, Edge*& outputEdge2);
    bool                        divide_face_into_triangular_faces_by_inserting_vertex( 
                                           Face* face, Vertex* vertex, 
                                           list<Edge*>& newEdges, list<Face*>& newFaces);
    bool                        divide_face_by_inserting_edge( 
                                            Face* face, Vertex* vtx1OfEdge, Vertex* vtx2OfEdge, 
                                            Edge*& newEdge, Face*& newFace1, Face*& newFace2);


    bool                        merge_two_faces_by_removing_edge( Face* face1, Face* face2 );



    /////////////////////////////////////////////////////////////////
    //
    //  use for only triangular mesh model.
    void                        collapse_edge( Edge* edgeToCollapse, const rg_Point3D& vertexPosition );
    bool                        does_nonmanifold_happen_by_edge_collapse( Edge* edgeToCollapse ) const;
    bool                        does_folding_face_happen_by_edge_collapse( Edge* edgeToCollapse, const rg_Point3D& vertexPosition ) const;
    //
    //////////////////////////////////////////////////////////////////

    void                        run_Loop_subdivision(const int& iteration);
    void                        run_Loop_subdivision_once_for_single_shell(PolygonMeshModel::Shell* shell);

    void                        triangulate();
    void                        simplify(const double& limitOfQEM); //  Garland's mesh decimation algorithm


protected:
    void duplicate(const PolygonMeshModel& polygonMeshModel);
        void duplicate_property( const PolygonMeshModel&                                   polygonMeshModel, 
                                 const map<SurfaceModel::Body*,   SurfaceModel::Body*>&    bodyMap,
                                 const map<SurfaceModel::Shell*,  SurfaceModel::Shell*>&   shellMap,
                                 const map<SurfaceModel::Face*,   SurfaceModel::Face*>&    faceMap,
                                 const map<SurfaceModel::Loop*,   SurfaceModel::Loop*>&    loopMap,
                                 const map<SurfaceModel::Edge*,   SurfaceModel::Edge*>&    edgeMap,
                                 const map<SurfaceModel::Vertex*, SurfaceModel::Vertex*>&  vertexMap);

    void  run_CatmallClark_subdivision_once();
    void  run_Loop_subdivision_once();
            void    compute_odd_and_even_vertices(
                            PolygonMeshModel::Shell* shell,
                            map<PolygonMeshModel::Vertex*, PolygonMeshModel::Edge*>& oddVertexMap,
                            map<PolygonMeshModel::Vertex*, rg_Point3D>&              evenVertexMap );


    bool  triangulate_single_face( PolygonMeshModel::Face* face );




};




///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Body 
//
class PolygonMeshModel::Body : public SurfaceModel::Body
{
public:
    Body();
    Body(const rg_INT& ID);
    Body(const Body& body);
    virtual ~Body();

    Body& operator =(const Body& body);

    PolygonMeshModel::Shell*    exterior_shell() const;
    unsigned int                interior_shell(list< PolygonMeshModel::Shell* >& interiorShells) const;

};



///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Shell 
//
class PolygonMeshModel::Shell : public SurfaceModel::Shell
{
private:
    bool    m_visible;

public:
    Shell();
    Shell(const rg_INT& ID);
    Shell(const Shell& shell);
    virtual ~Shell();

    Shell& operator =(const Shell& shell);

    PolygonMeshModel::Body*     body() const;

    bool                        is_visible() const;
    void                        is_visible(const bool& visibleOrNot);

    unsigned int                find_bounding_faces(    list<PolygonMeshModel::Face*>&   boundingFaces) const;
    unsigned int                find_bounding_edges(    list<PolygonMeshModel::Edge*>&   boundingEdges) const;
    unsigned int                find_bounding_vertices( list<PolygonMeshModel::Vertex*>& boundingVertices) const;


    AxisAlignedBox              bounding_box() const;

    double                      volume() const;
    double                      area() const;
};



///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Face 
//
class PolygonMeshModel::Face : public SurfaceModel::Face
{
private:
    typedef PolygonMeshModel::Edge   Edge;
    typedef PolygonMeshModel::Vertex Vertex;

    rg_Point3D  m_normal;
    int         m_faceGroupID;

public:
    Face();
    Face(const rg_INT& ID);
    Face(const Face& face);
    virtual ~Face();

    Face& operator =(const Face& face);

    PolygonMeshModel::Edge*   first_edge() const;

    rg_Point3D  normal() const;
    void        set_normal(const rg_Point3D& normal);

    int         face_group_ID() const;
    void        set_face_group_ID(const int& groupID);
    Edge*       find_edge(const Vertex* const vertex1, const Vertex* const vertex2 ) const;

    PolygonMeshModel::Shell*  shell() const;

    unsigned int  find_bounding_vertices( list<PolygonMeshModel::Vertex*>& boundingVertices) const;
    unsigned int  find_bounding_vertices( vector<PolygonMeshModel::Vertex*>& boundingVertices) const;
    unsigned int  find_bounding_edges(    list<PolygonMeshModel::Edge*>&   boundingEdges) const;

    unsigned int  find_edges_sharing_with_adjacent_face( const PolygonMeshModel::Face* const face, list<PolygonMeshModel::Edge*>& sharingEdges) const;


    Plane           plane() const;

    void            compute_normal();
};



///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Loop 
//
class PolygonMeshModel::Loop : public SurfaceModel::Loop
{
public:
    Loop();
    Loop(const rg_INT& ID);
    Loop(const Loop& loop);
    virtual ~Loop();

    Loop& operator =(const Loop& loop);
};



///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Edge 
//
class PolygonMeshModel::Edge : public SurfaceModel::Edge
{
private:
    bool    m_boundary_of_faceGroup;

public:
    Edge();
    Edge(const rg_INT& ID);
    Edge(const Edge& edge);
    virtual ~Edge();

    Edge& operator =(const Edge& edge);

    PolygonMeshModel::Shell*  shell() const;

    bool    is_boundary_of_faceGroup() const;
    void    is_boundary_of_faceGroup(const bool& boundaryOrNot);

    PolygonMeshModel::Vertex* start_vertex() const;
    PolygonMeshModel::Vertex* end_vertex() const;
    PolygonMeshModel::Face*   right_face() const;
    PolygonMeshModel::Face*   left_face() const;
    PolygonMeshModel::Edge*   right_hand_edge() const;
    PolygonMeshModel::Edge*   left_hand_edge() const;
    PolygonMeshModel::Edge*   right_leg_edge() const;
    PolygonMeshModel::Edge*   left_leg_edge() const;

    PolygonMeshModel::Vertex* opposite_vertex(const PolygonMeshModel::Vertex* const vertex) const;
    PolygonMeshModel::Face*   opposite_face(const PolygonMeshModel::Face* const face) const;

    rg_INT  find_edges_in_star(     list<PolygonMeshModel::Edge*>&   edgesInStar) const;
    rg_INT  find_incident_faces(    list<PolygonMeshModel::Face*>&   incidentFaces) const;

};



///////////////////////////////////////////////////////////////////////////////
//
// class PolygonMeshModel::Vertex 
//
class PolygonMeshModel::Vertex : public SurfaceModel::Vertex
{
private:
    rg_Point3D  m_coordinate;

public:
    Vertex();
    Vertex(const rg_INT& ID);
    Vertex(const Vertex& vertex);
    virtual ~Vertex();

    Vertex& operator =(const Vertex& vertex);

    rg_Point3D  coordinate() const;
    void        set_coordinate(const rg_Point3D& coordinate);

    PolygonMeshModel::Shell*  shell() const;

    rg_Point3D  compute_normal() const;

    rg_INT      find_adjacent_vertices( list<PolygonMeshModel::Vertex*>& adjacentVertices) const;
    rg_INT      find_incident_edges(    list<PolygonMeshModel::Edge*>&   incidentEdges) const;
    rg_INT      find_incident_faces(    list<PolygonMeshModel::Face*>&   incidentFaces) const;

    rg_INT      find_edges_in_star(     list<PolygonMeshModel::Edge*>&   edgesInStar) const;
    rg_INT      find_edges_in_shell(    list<PolygonMeshModel::Edge*>&   edgesInShell) const;
};



inline  PolygonMeshModel::Body*     PolygonMeshModel::Shell::body() const           { return (PolygonMeshModel::Body*)SurfaceModel::Shell::body(); }
inline  bool                        PolygonMeshModel::Shell::is_visible() const     { return m_visible; };
inline  void                        PolygonMeshModel::Shell::is_visible(const bool& visibleOrNot) { m_visible = visibleOrNot; };

inline  PolygonMeshModel::Edge*     PolygonMeshModel::Face::first_edge() const      { return (PolygonMeshModel::Edge*)SurfaceModel::Face::first_edge();     }
inline  int                         PolygonMeshModel::Face::face_group_ID() const   { return m_faceGroupID; }
inline  void                        PolygonMeshModel::Face::set_face_group_ID(const int& groupID) { m_faceGroupID = groupID; }
inline  PolygonMeshModel::Shell*    PolygonMeshModel::Face::shell() const           { return (PolygonMeshModel::Shell*)SurfaceModel::Face::shell(); }
inline  Plane                       PolygonMeshModel::Face::plane() const           { return Plane( m_normal, first_edge()->start_vertex()->coordinate() );     }

inline  PolygonMeshModel::Vertex*   PolygonMeshModel::Edge::start_vertex() const    { return (PolygonMeshModel::Vertex*)SurfaceModel::Edge::start_vertex();     }
inline  PolygonMeshModel::Vertex*   PolygonMeshModel::Edge::end_vertex() const      { return (PolygonMeshModel::Vertex*)SurfaceModel::Edge::end_vertex();       }
inline  PolygonMeshModel::Face*     PolygonMeshModel::Edge::right_face() const      { return (PolygonMeshModel::Face*)  SurfaceModel::Edge::right_face();       }
inline  PolygonMeshModel::Face*     PolygonMeshModel::Edge::left_face() const       { return (PolygonMeshModel::Face*)  SurfaceModel::Edge::left_face();        }
inline  PolygonMeshModel::Edge*     PolygonMeshModel::Edge::right_hand_edge() const { return (PolygonMeshModel::Edge*)  SurfaceModel::Edge::right_hand_edge();  }
inline  PolygonMeshModel::Edge*     PolygonMeshModel::Edge::left_hand_edge() const  { return (PolygonMeshModel::Edge*)  SurfaceModel::Edge::left_hand_edge();   }
inline  PolygonMeshModel::Edge*     PolygonMeshModel::Edge::right_leg_edge() const  { return (PolygonMeshModel::Edge*)  SurfaceModel::Edge::right_leg_edge();   }
inline  PolygonMeshModel::Edge*     PolygonMeshModel::Edge::left_leg_edge() const   { return (PolygonMeshModel::Edge*)  SurfaceModel::Edge::left_leg_edge();    }
inline  bool                        PolygonMeshModel::Edge::is_boundary_of_faceGroup() const { return m_boundary_of_faceGroup; }
inline  void                        PolygonMeshModel::Edge::is_boundary_of_faceGroup(const bool& boundaryOrNot) { m_boundary_of_faceGroup = boundaryOrNot; }

#endif

