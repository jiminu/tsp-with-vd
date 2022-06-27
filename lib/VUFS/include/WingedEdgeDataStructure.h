#ifndef _WINGEDEDGEDATASTRUCTURE_H
#define _WINGEDEDGEDATASTRUCTURE_H

#include "rg_Const.h"
#include "TopologicalEntity.h"

#include <list>
#include <map>
using namespace std;

class WingedEdgeDataStructure
{
public:
    class Face;
    class Edge;
    class Vertex;

private:
    list<Face*>     m_faces;
    list<Edge*>     m_edges;
    list<Vertex*>   m_vertices;

    unsigned int    m_number_of_faces;
    unsigned int    m_number_of_edges;
    unsigned int    m_number_of_vertices;

public:
    WingedEdgeDataStructure();
    WingedEdgeDataStructure(const WingedEdgeDataStructure& weds);
    ~WingedEdgeDataStructure();

    void    clear();

    WingedEdgeDataStructure& operator =(const WingedEdgeDataStructure& weds);

    const list<Face*>&          get_all_faces() const;
    const list<Edge*>&          get_all_edges() const;
    const list<Vertex*>&        get_all_vertices() const;

    rg_INT          number_of_faces() const;
    rg_INT          number_of_edges() const;
    rg_INT          number_of_vertices() const;

    Face*           create_face(        const Face&        face );
    Edge*           create_edge(        const Edge&        edge );
    Vertex*         create_vertex(      const Vertex&      vertex );
    virtual Face*   create_face();
    virtual Edge*   create_edge();
    virtual Vertex* create_vertex();

    void            concatenate_face(   const Face*        const face );
    void            concatenate_edge(   const Edge*        const edge );
    void            concatenate_vertex( const Vertex*      const vertex );

    void            remove_face(        const Face*        const face );
    void            remove_edge(        const Edge*        const edge );
    void            remove_vertex(      const Vertex*      const vertex );


protected:
    void            duplicate(const WingedEdgeDataStructure& reds);
    void            make_entity_map_for_duplication( const WingedEdgeDataStructure& weds, 
                                                     map<Face*,        Face*>&          faceMap,
                                                     map<Edge*,        Edge*>&          edgeMap,
                                                     map<Vertex*,      Vertex*>&        vertexMap);
    void            duplicate_topology(              const WingedEdgeDataStructure& weds, 
                                                     const map<Face*,        Face*>&    faceMap,
                                                     const map<Edge*,        Edge*>&    edgeMap,
                                                     const map<Vertex*,      Vertex*>&  vertexMap);

};



///////////////////////////////////////////////////////////////////////////////
//
// class WingedEdgeDataStructureCore::Face : public TopologicalEntity
//

class WingedEdgeDataStructure::Face : public TopologicalEntity
{
protected:
    WingedEdgeDataStructure::Edge* m_firstEdge;

public:
    Face();
    Face(const rg_INT& ID);
    Face(const Face& face);
    virtual ~Face();

    Face&   operator =(const WingedEdgeDataStructure::Face& face);

    WingedEdgeDataStructure::Edge*   first_edge() const;
    void    set_first_edge(const WingedEdgeDataStructure::Edge* const edge);

    bool    is_incident_to(const WingedEdgeDataStructure::Vertex* const vertex) const;
    bool    is_incident_to(const WingedEdgeDataStructure::Edge* const edge) const;
    bool    is_adjacent_to(const WingedEdgeDataStructure::Face* const face) const;

    WingedEdgeDataStructure::Edge*   find_edge(const WingedEdgeDataStructure::Vertex* const vertex1, const WingedEdgeDataStructure::Vertex* const vertex2 ) const;

    rg_INT  number_of_bounding_vertices() const;
    rg_INT  number_of_bounding_edges() const;
    rg_INT  number_of_adjacent_faces() const;

    rg_INT  find_bounding_vertices( list<WingedEdgeDataStructure::Vertex*>& boundingVertices) const;
    rg_INT  find_bounding_edges(    list<WingedEdgeDataStructure::Edge*>&   boundingEdges) const;
    rg_INT  find_adjacent_faces(    list<WingedEdgeDataStructure::Face*>&   adjacentFaces) const;
};



///////////////////////////////////////////////////////////////////////////////
//
// class WingedEdgeDataStructureCore::Edge : public TopologicalEntity
//

class WingedEdgeDataStructure::Edge : public TopologicalEntity
{
protected:
    WingedEdgeDataStructure::Vertex* m_startVertex;
    WingedEdgeDataStructure::Vertex* m_endVertex;
    
    WingedEdgeDataStructure::Face*   m_rightFace;
    WingedEdgeDataStructure::Face*   m_leftFace;

    WingedEdgeDataStructure::Edge*   m_rightHand;
    WingedEdgeDataStructure::Edge*   m_leftHand;
    WingedEdgeDataStructure::Edge*   m_rightLeg;
    WingedEdgeDataStructure::Edge*   m_leftLeg;

public:
    Edge();
    Edge(const rg_INT& ID);
    Edge(const WingedEdgeDataStructure::Edge& edge);
    virtual ~Edge();

    Edge& operator =(const WingedEdgeDataStructure::Edge& edge);

    WingedEdgeDataStructure::Vertex* start_vertex() const;
    WingedEdgeDataStructure::Vertex* end_vertex() const;   
    WingedEdgeDataStructure::Face*   right_face() const;
    WingedEdgeDataStructure::Face*   left_face() const;
    WingedEdgeDataStructure::Edge*   right_hand_edge() const;
    WingedEdgeDataStructure::Edge*   left_hand_edge() const;
    WingedEdgeDataStructure::Edge*   right_leg_edge() const;
    WingedEdgeDataStructure::Edge*   left_leg_edge() const;

    WingedEdgeDataStructure::Vertex* opposite_vertex(const WingedEdgeDataStructure::Vertex* const vertex) const;
    WingedEdgeDataStructure::Face*   opposite_face(const WingedEdgeDataStructure::Face* const face) const;

    void    set_start_vertex(    const WingedEdgeDataStructure::Vertex* const vertex);
    void    set_end_vertex(      const WingedEdgeDataStructure::Vertex* const vertex);   
    void    set_right_face(      const WingedEdgeDataStructure::Face* const face);
    void    set_left_face(       const WingedEdgeDataStructure::Face* const face);
    void    set_right_hand_edge( const WingedEdgeDataStructure::Edge* const edge);
    void    set_left_hand_edge(  const WingedEdgeDataStructure::Edge* const edge);
    void    set_right_leg_edge(  const WingedEdgeDataStructure::Edge* const edge);
    void    set_left_leg_edge(   const WingedEdgeDataStructure::Edge* const edge);

    void    set_edge(  const WingedEdgeDataStructure::Vertex* const startVertex ,          
                       const WingedEdgeDataStructure::Vertex* const endVertex,
                       const WingedEdgeDataStructure::Face*   const rightFace       = rg_NULL,     
                       const WingedEdgeDataStructure::Face*   const leftFace        = rg_NULL,
                       const WingedEdgeDataStructure::Edge*   const rightHandEdge   = rg_NULL, 
                       const WingedEdgeDataStructure::Edge*   const leftHandEdge    = rg_NULL,
                       const WingedEdgeDataStructure::Edge*   const rightLegEdge    = rg_NULL,  
                       const WingedEdgeDataStructure::Edge*   const leftLegEdge     = rg_NULL );

    bool    is_start_vertex(const WingedEdgeDataStructure::Vertex* const vertex) const;
    bool    is_end_vertex(  const WingedEdgeDataStructure::Vertex* const vertex) const;

    bool    is_incident_to(const WingedEdgeDataStructure::Vertex* const vertex) const;
    bool    is_adjacent_to(const WingedEdgeDataStructure::Edge* const edge) const;
    bool    is_member_of(  const WingedEdgeDataStructure::Face* const face) const;

    rg_INT  find_bounding_vertices( list<WingedEdgeDataStructure::Vertex*>& boundingVertices) const;
    rg_INT  find_adjacent_edges(    list<WingedEdgeDataStructure::Edge*>&   adjacentEdges) const;
    rg_INT  find_incident_faces(    list<WingedEdgeDataStructure::Face*>&   incidentFaces) const;

    bool    is_on_boundary() const;

    void    reverse();
    void    flip();

    void    connect_right_hand_edge( WingedEdgeDataStructure::Edge* const rightHand);
    void    connect_left_hand_edge(  WingedEdgeDataStructure::Edge* const leftHand);
    void    connect_right_leg_edge(  WingedEdgeDataStructure::Edge* const rightLeg);
    void    connect_left_leg_edge(   WingedEdgeDataStructure::Edge* const leftLeg);


    rg_INT  find_edges_in_star(     list<WingedEdgeDataStructure::Edge*>&   edgesInStar) const;

};



///////////////////////////////////////////////////////////////////////////////
//
// class WingedEdgeDataStructureCore::Vertex : public TopologicalEntity
//

class WingedEdgeDataStructure::Vertex : public TopologicalEntity
{
protected:
    WingedEdgeDataStructure::Edge* m_firstEdge;

public:
    Vertex();
    Vertex(const rg_INT& ID);
    Vertex(const WingedEdgeDataStructure::Vertex& vertex);
    virtual ~Vertex();

    Vertex& operator =(const WingedEdgeDataStructure::Vertex& vertex);

    WingedEdgeDataStructure::Edge*   first_edge() const;
    void    set_first_edge(const WingedEdgeDataStructure::Edge* const edge);

    bool    is_adjacent_to(const WingedEdgeDataStructure::Vertex* const vertex) const;
    bool    is_member_of(  const WingedEdgeDataStructure::Edge* const edge) const;
    bool    is_member_of(  const WingedEdgeDataStructure::Face* const face) const;

    bool    is_start_vertex_of(const WingedEdgeDataStructure::Edge* const edge) const;
    bool    is_end_vertex_of(  const WingedEdgeDataStructure::Edge* const edge) const;

    bool    is_on_boundary() const;

    rg_INT  number_of_adjacent_vertices() const;
    rg_INT  number_of_incident_edges() const;
    rg_INT  number_of_incident_faces() const;

    rg_INT  find_adjacent_vertices( list<WingedEdgeDataStructure::Vertex*>& adjacentVertices) const;
    rg_INT  find_incident_edges(    list<WingedEdgeDataStructure::Edge*>&   incidentEdges) const;
    rg_INT  find_incident_faces(    list<WingedEdgeDataStructure::Face*>&   incidentFaces) const;

    rg_INT  find_edges_in_star(     list<WingedEdgeDataStructure::Edge*>&   edgesInStar) const;
    rg_INT  find_edges_in_shell(    list<WingedEdgeDataStructure::Edge*>&   edgesInShell) const;
};





///////////////////////////////////////////////////////////////////////////////
//
// inline functions for class WingedEdgeDataStructureCore 
//
inline  const list<WingedEdgeDataStructure::Face*>&   WingedEdgeDataStructure::get_all_faces() const    { return m_faces; }
inline  const list<WingedEdgeDataStructure::Edge*>&   WingedEdgeDataStructure::get_all_edges() const    { return m_edges; }
inline  const list<WingedEdgeDataStructure::Vertex*>& WingedEdgeDataStructure::get_all_vertices() const { return m_vertices; }

inline  rg_INT WingedEdgeDataStructure::number_of_faces() const     { return m_number_of_faces; }
inline  rg_INT WingedEdgeDataStructure::number_of_edges() const     { return m_number_of_edges; }
inline  rg_INT WingedEdgeDataStructure::number_of_vertices() const  { return m_number_of_vertices; }




///////////////////////////////////////////////////////////////////////////////
//
// inline functions for class WingedEdgeDataStructureCore::Face 
//
inline  WingedEdgeDataStructure::Edge*   WingedEdgeDataStructure::Face::first_edge() const { return m_firstEdge; }

inline  void    WingedEdgeDataStructure::Face::set_first_edge(const WingedEdgeDataStructure::Edge* const edge) {
    m_firstEdge = const_cast<WingedEdgeDataStructure::Edge*>(edge);
}



///////////////////////////////////////////////////////////////////////////////
//
// inline functions for class WingedEdgeDataStructureCore::Edge 
//
inline  WingedEdgeDataStructure::Vertex* WingedEdgeDataStructure::Edge::start_vertex() const    { return m_startVertex; }
inline  WingedEdgeDataStructure::Vertex* WingedEdgeDataStructure::Edge::end_vertex() const      { return m_endVertex;   }
inline  WingedEdgeDataStructure::Face*   WingedEdgeDataStructure::Edge::right_face() const      { return m_rightFace;   }
inline  WingedEdgeDataStructure::Face*   WingedEdgeDataStructure::Edge::left_face() const       { return m_leftFace;    }
inline  WingedEdgeDataStructure::Edge*   WingedEdgeDataStructure::Edge::right_hand_edge() const { return m_rightHand;   }
inline  WingedEdgeDataStructure::Edge*   WingedEdgeDataStructure::Edge::left_hand_edge() const  { return m_leftHand;    }
inline  WingedEdgeDataStructure::Edge*   WingedEdgeDataStructure::Edge::right_leg_edge() const  { return m_rightLeg;    }
inline  WingedEdgeDataStructure::Edge*   WingedEdgeDataStructure::Edge::left_leg_edge() const   { return m_leftLeg;     }


inline  void    WingedEdgeDataStructure::Edge::set_start_vertex(    const WingedEdgeDataStructure::Vertex* const vertex) {
    m_startVertex = const_cast<WingedEdgeDataStructure::Vertex*>(vertex);
}

inline  void    WingedEdgeDataStructure::Edge::set_end_vertex(      const WingedEdgeDataStructure::Vertex* const vertex) {
    m_endVertex = const_cast<WingedEdgeDataStructure::Vertex*>(vertex);
}

inline  void    WingedEdgeDataStructure::Edge::set_right_face(      const WingedEdgeDataStructure::Face* const face) {
    m_rightFace = const_cast<WingedEdgeDataStructure::Face*>(face);
}

inline  void    WingedEdgeDataStructure::Edge::set_left_face(       const WingedEdgeDataStructure::Face* const face) {
    m_leftFace = const_cast<WingedEdgeDataStructure::Face*>(face);
}

inline  void    WingedEdgeDataStructure::Edge::set_right_hand_edge( const WingedEdgeDataStructure::Edge* const edge) {
    m_rightHand = const_cast<WingedEdgeDataStructure::Edge*>(edge);
}

inline  void    WingedEdgeDataStructure::Edge::set_left_hand_edge(  const WingedEdgeDataStructure::Edge* const edge) {
    m_leftHand = const_cast<WingedEdgeDataStructure::Edge*>(edge);
}

inline  void    WingedEdgeDataStructure::Edge::set_right_leg_edge(  const WingedEdgeDataStructure::Edge* const edge) {
    m_rightLeg = const_cast<WingedEdgeDataStructure::Edge*>(edge);
}

inline  void    WingedEdgeDataStructure::Edge::set_left_leg_edge(   const WingedEdgeDataStructure::Edge* const edge) {
    m_leftLeg = const_cast<WingedEdgeDataStructure::Edge*>(edge);
}


inline  bool    WingedEdgeDataStructure::Edge::is_start_vertex(const WingedEdgeDataStructure::Vertex* const vertex) const {
    return (vertex==m_startVertex) ? true : false;
}

inline  bool    WingedEdgeDataStructure::Edge::is_end_vertex(  const WingedEdgeDataStructure::Vertex* const vertex) const {
    return (vertex==m_endVertex) ? true : false;
}

inline  bool    WingedEdgeDataStructure::Edge::is_incident_to(const WingedEdgeDataStructure::Vertex* const vertex) const  {
    return (vertex==m_startVertex || vertex==m_endVertex ) ? true : false;
}

inline  bool    WingedEdgeDataStructure::Edge::is_adjacent_to(const WingedEdgeDataStructure::Edge* const edge) const {
    return (edge->is_incident_to(m_startVertex) || edge->is_incident_to(m_endVertex) ) ? true : false;
}

inline  bool    WingedEdgeDataStructure::Edge::is_member_of(  const WingedEdgeDataStructure::Face* const face) const {
    return (face==m_rightFace || face==m_leftFace) ? true : false;
}

inline  bool    WingedEdgeDataStructure::Edge::is_on_boundary() const {
    return ( m_rightFace==rg_NULL || m_leftFace==rg_NULL ) ? true : false;
}


///////////////////////////////////////////////////////////////////////////////
//
// inline functions for class WingedEdgeDataStructureCore::Vertex 
//
inline  WingedEdgeDataStructure::Edge*   WingedEdgeDataStructure::Vertex::first_edge() const { return m_firstEdge; }

inline  void    WingedEdgeDataStructure::Vertex::set_first_edge(const WingedEdgeDataStructure::Edge* const edge) {
    m_firstEdge = const_cast<WingedEdgeDataStructure::Edge*>(edge);
}


inline  bool    WingedEdgeDataStructure::Vertex::is_member_of(  const WingedEdgeDataStructure::Edge* const edge) const {
    return edge->is_incident_to( this );
}

inline  bool    WingedEdgeDataStructure::Vertex::is_member_of(  const WingedEdgeDataStructure::Face* const face) const {
    return face->is_incident_to( this );
}

inline  bool    WingedEdgeDataStructure::Vertex::is_start_vertex_of(const WingedEdgeDataStructure::Edge* const edge) const {
    return ( this == edge->start_vertex() ) ? true : false;
}

inline  bool    WingedEdgeDataStructure::Vertex::is_end_vertex_of(  const WingedEdgeDataStructure::Edge* const edge) const {
    return ( this == edge->end_vertex() ) ? true : false;
}


#endif

