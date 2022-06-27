#ifndef SURFACEMODEL_H
#define SURFACEMODEL_H

#include "rg_Const.h"
#include "TopologicalEntity.h"

#include "WingedEdgeDataStructure.h"

#include <list>
#include <vector>
#include <map>
using namespace std;



class SurfaceModel
{
public:
    class Body;
    class Shell;
    class Face;
    class Loop;
    class Edge;
    class Vertex;

protected:
    list< Body* >   m_bodies;
    list< Shell* >  m_shells;
    list< Face* >   m_faces;
    list< Loop* >   m_loops;
    list< Edge* >   m_edges;
    list< Vertex* > m_vertices;

    unsigned int    m_number_of_bodies;
    unsigned int    m_number_of_shells;
    unsigned int    m_number_of_faces;
    unsigned int    m_number_of_loops;
    unsigned int    m_number_of_edges;
    unsigned int    m_number_of_vertices;

public:
    SurfaceModel();
    SurfaceModel(const SurfaceModel& surfaceModel);
    ~SurfaceModel();

    void clear();
    void release_deleted_entities();

    SurfaceModel& operator =(const SurfaceModel& surfaceModel);

    const list< Body* >&    get_all_bodies() const;
    const list< Shell* >&   get_all_shells() const;
    const list< Face* >&    get_all_faces() const;
    const list< Loop* >&    get_all_loops() const;
    const list< Edge* >&    get_all_edges() const;
    const list< Vertex* >&  get_all_vertices() const;

    unsigned int            number_of_bodies() const;
    unsigned int            number_of_shells() const;
    unsigned int            number_of_faces() const;
    unsigned int            number_of_loops() const;
    unsigned int            number_of_edges() const;
    unsigned int            number_of_vertices() const;

    virtual SurfaceModel::Body*     create_body();
    virtual SurfaceModel::Shell*    create_shell();
    virtual SurfaceModel::Face*     create_face();
    virtual SurfaceModel::Loop*     create_loop();
    virtual SurfaceModel::Edge*     create_edge();
    virtual SurfaceModel::Vertex*   create_vertex();

    void                    concatenate_body(   const SurfaceModel::Body*       const body );
    void                    concatenate_shell(  const SurfaceModel::Shell*      const shell );
    void                    concatenate_face(   const SurfaceModel::Face*       const face );
    void                    concatenate_loop(   const SurfaceModel::Loop*       const loop );
    void                    concatenate_edge(   const SurfaceModel::Edge*       const edge );
    void                    concatenate_vertex( const SurfaceModel::Vertex*     const vertex );

    void                    remove_body(        const SurfaceModel::Body*       const body );
    void                    remove_shell(       const SurfaceModel::Shell*      const shell );
    void                    remove_face(        const SurfaceModel::Face*       const face );
    void                    remove_loop(        const SurfaceModel::Loop*       const loop );
    void                    remove_edge(        const SurfaceModel::Edge*       const edge );
    void                    remove_vertex(      const SurfaceModel::Vertex*     const vertex );

    void                    splice(SurfaceModel::Shell* shellToExtend, SurfaceModel::Shell* shellToRemove);

    bool                    divide_edge_by_inserting_vertex( Edge* edge, Vertex* splitVertex, Edge*& outputEdge1, Edge*& outputEdge2);
    bool                    divide_face_into_triangular_faces_by_inserting_vertex( 
                                        Face* face, Vertex* vertex, 
                                        list<Edge*>& newEdges, list<Face*>& newFaces);
    bool                    divide_face_by_inserting_edge( 
                                        Face* face, Vertex* vtx1OfEdge, Vertex* vtx2OfEdge, 
                                        Edge* &newEdge, Face* &newFace1, Face* &newFace2);


    bool                    merge_two_faces_by_removing_edge( Face* face1, Face* face2 );


protected:
    void duplicate(const SurfaceModel& surfaceModel);

        void            make_entity_map_for_duplication( const SurfaceModel&                surfaceModel, 
                                                               map<Body*,        Body*>&    bodyMap,
                                                               map<Shell*,       Shell*>&   shellMap,
                                                               map<Face*,        Face*>&    faceMap,
                                                               map<Loop*,        Loop*>&    loopMap,
                                                               map<Edge*,        Edge*>&    edgeMap,
                                                               map<Vertex*,      Vertex*>&  vertexMap);
        void            duplicate_topology(              const SurfaceModel&                surfaceModel, 
                                                         const map<Body*,        Body*>&    bodyMap,
                                                         const map<Shell*,       Shell*>&   shellMap,
                                                         const map<Face*,        Face*>&    faceMap,
                                                         const map<Loop*,        Loop*>&    loopMap,
                                                         const map<Edge*,        Edge*>&    edgeMap,
                                                         const map<Vertex*,      Vertex*>&  vertexMap);

};


///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Body 
//
class SurfaceModel::Body : public TopologicalEntity
{
private:
    SurfaceModel::Shell*            m_exteriorShell;
    list< SurfaceModel::Shell* >    m_interiorShell;

public:
    Body();
    Body(const rg_INT& ID);
    Body(const Body& body);
    virtual ~Body();

    Body& operator =(const Body& body);

    SurfaceModel::Shell*                 exterior_shell() const;
    const list< SurfaceModel::Shell* >&  interior_shell() const;

    void                    set_exterior_shell(        const SurfaceModel::Shell* const shell);
    void                    add_interior_shell(    const SurfaceModel::Shell* const shell);
    void                    remove_interior_shell( const SurfaceModel::Shell* const shell);

    rg_INT                  find_bounding_shells(         list<SurfaceModel::Shell*>&  boundingShells) const;
    rg_INT                  find_bounding_faces(          list<SurfaceModel::Face*>&   boundingFaces) const;
    rg_INT                  find_loops_of_bounding_faces( list<SurfaceModel::Loop*>&   loopsOfBoundingFaces) const;
    rg_INT                  find_bounding_edges(          list<SurfaceModel::Edge*>&   boundingEdges) const;
    rg_INT                  find_bounding_vertices(       list<SurfaceModel::Vertex*>& boundingVertices) const;

};

///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Shell 
//
class SurfaceModel::Shell : public TopologicalEntity
{
private:
    SurfaceModel::Body*             m_body;
    list< SurfaceModel::Face* >     m_faces;

    unsigned int                    m_number_of_faces;

public:
    Shell();
    Shell(const rg_INT& ID);
    Shell(const Shell& shell);
    virtual ~Shell();

    Shell& operator =(const Shell& shell);

    SurfaceModel::Body*     body() const;
    void                    set_body(const SurfaceModel::Body* const attachedBody);

    const list< SurfaceModel::Face* >& faces() const;

    rg_INT                  number_of_faces() const;
    rg_INT                  number_of_edges() const;
    rg_INT                  number_of_vertices() const;

    void                    add_face(    const SurfaceModel::Face*  const face);
    void                    remove_face( const SurfaceModel::Face*  const face);

    bool                    is_exterior_shell() const;

    rg_INT                  find_bounding_faces(          list<SurfaceModel::Face*>&   boundingFaces) const;
    rg_INT                  find_loops_of_bounding_faces( list<SurfaceModel::Loop*>&   loopsOfBoundingFaces) const;
    rg_INT                  find_bounding_edges(          list<SurfaceModel::Edge*>&   boundingEdges) const;
    rg_INT                  find_bounding_vertices(       list<SurfaceModel::Vertex*>& boundingVertices) const;

};


///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Face 
//
class SurfaceModel::Face : public WingedEdgeDataStructure::Face
{
private:
    SurfaceModel::Shell*            m_shell;
    SurfaceModel::Loop*             m_exteriorLoop;
    list< SurfaceModel::Loop* >     m_interiorLoops;

public:
    Face();
    Face(const rg_INT& ID);
    Face(const Face& face);
    virtual ~Face();

    Face& operator =(const Face& face);

    SurfaceModel::Shell*                shell() const;
    SurfaceModel::Loop*                 exterior_loop() const;
    const list< SurfaceModel::Loop* >&  interior_loops() const;

    void    set_shell(                const SurfaceModel::Shell* const attachedShell);
    void    set_exterior_loop(        const SurfaceModel::Loop*  const exteriorLoop);
    void    add_interior_loop(    const SurfaceModel::Loop*  const interiorLoop);
    void    remove_interior_loop( const SurfaceModel::Loop*  const interiorLoop);

    bool    no_interior_loop() const;

    rg_INT  number_of_loops() const;
    rg_INT  number_of_bounding_vertices() const;
    rg_INT  number_of_bounding_edges() const;
    rg_INT  number_of_adjacent_faces() const;


    SurfaceModel::Edge*   first_edge() const;
    SurfaceModel::Edge*   find_edge(const SurfaceModel::Vertex* const vertex1, const SurfaceModel::Vertex* const vertex2 ) const;

    rg_INT  find_bounding_vertices( list<SurfaceModel::Vertex*>& boundingVertices) const;
    rg_INT  find_bounding_edges(    list<SurfaceModel::Edge*>&   boundingEdges) const;
    rg_INT  find_bounding_loops(    list<SurfaceModel::Loop*>&   boundingLoops) const;
    rg_INT  find_adjacent_faces(    list<SurfaceModel::Face*>&   adjacentFaces) const;

    rg_INT  find_bounding_edges_incident_to_vertex(const SurfaceModel::Vertex* const vertex, SurfaceModel::Edge*& prevEdge, SurfaceModel::Edge*& nextEdge) const;

    rg_INT  find_edges_sharing_with_adjacent_face( const SurfaceModel::Face* const face, list<SurfaceModel::Edge*>& sharingEdges) const;

};

///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Loop 
//
class SurfaceModel::Loop : public TopologicalEntity
{
private:
    SurfaceModel::Face*     m_face;
    SurfaceModel::Edge*     m_firstEdge;

public:
    Loop();
    Loop(const rg_INT& ID);
    Loop(const Loop& loop);
    virtual ~Loop();

    Loop& operator =(const Loop& loop);

    SurfaceModel::Face*     face() const;
    SurfaceModel::Edge*     first_edge() const;

    void                    set_face(const SurfaceModel::Face* const attachedFace);
    void                    set_first_edge(const SurfaceModel::Edge* const firstEdge);

    bool    is_exterior_loop() const;

    rg_INT  find_bounding_vertices( list<SurfaceModel::Vertex*>& boundingVertices) const;
    rg_INT  find_bounding_edges(    list<SurfaceModel::Edge*>&   boundingEdges) const;
    rg_INT  find_incident_faces(    list<SurfaceModel::Face*>&   incidentFaces) const;
};

///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Edge 
//
class SurfaceModel::Edge : public WingedEdgeDataStructure::Edge
{
public:
    Edge();
    Edge(const rg_INT& ID);
    Edge(const Edge& edge);
    virtual ~Edge();

    Edge& operator =(const Edge& edge);

    SurfaceModel::Shell*  shell() const;

    SurfaceModel::Vertex* start_vertex() const;
    SurfaceModel::Vertex* end_vertex() const;   
    SurfaceModel::Face*   right_face() const;
    SurfaceModel::Face*   left_face() const;
    SurfaceModel::Edge*   right_hand_edge() const;
    SurfaceModel::Edge*   left_hand_edge() const;
    SurfaceModel::Edge*   right_leg_edge() const;
    SurfaceModel::Edge*   left_leg_edge() const;

    SurfaceModel::Vertex* opposite_vertex(const SurfaceModel::Vertex* const vertex) const;
    SurfaceModel::Face*   opposite_face(const SurfaceModel::Face* const face) const;

    rg_INT  find_bounding_vertices( list<SurfaceModel::Vertex*>& boundingVertices) const;
    rg_INT  find_adjacent_edges(    list<SurfaceModel::Edge*>&   adjacentEdges) const;
    rg_INT  find_incident_faces(    list<SurfaceModel::Face*>&   incidentFaces) const;

    rg_INT  find_edges_in_star(     list<SurfaceModel::Edge*>&   edgesInStar) const;

};

///////////////////////////////////////////////////////////////////////////////
//
// class SurfaceModel::Vertex 
//
class SurfaceModel::Vertex : public WingedEdgeDataStructure::Vertex
{
public:
    Vertex();
    Vertex(const rg_INT& ID);
    Vertex(const Vertex& vertex);
    virtual ~Vertex();

    Vertex& operator =(const Vertex& vertex);

    SurfaceModel::Shell*  shell() const;

    SurfaceModel::Edge*   first_edge() const;

    rg_INT  find_adjacent_vertices( list<SurfaceModel::Vertex*>& adjacentVertices) const;
    rg_INT  find_incident_edges(    list<SurfaceModel::Edge*>&   incidentEdges) const;
    rg_INT  find_incident_faces(    list<SurfaceModel::Face*>&   incidentFaces) const;

    rg_INT  find_edges_in_star(     list<SurfaceModel::Edge*>&   edgesInStar) const;
    rg_INT  find_edges_in_shell(    list<SurfaceModel::Edge*>&   edgesInShell) const;

};

#endif


