#ifndef _SIMPLIFICATIONBYGARLAN_H
#define _SIMPLIFICATIONBYGARLAN_H


#include "PolygonMeshModel.h"


#include <map>
#include <list>
using namespace std;


class SimplificationByGarlan
{
public:
    class   Edge;
    class   Vertex;

private:
    PolygonMeshModel*   m_model;
    double              m_limitOfQEM;


public:
    SimplificationByGarlan();
    ~SimplificationByGarlan();

    void    run_simplification( PolygonMeshModel*   model, const double& limitOfQEM );
    void    simplify_shell( PolygonMeshModel::Shell* shell);

    void    set_simplification(PolygonMeshModel* model, const double& limitOfQEM);

private:
    void    initialize( PolygonMeshModel::Shell* shell,
                            map<PolygonMeshModel::Edge*, SimplificationByGarlan::Edge>& edgeMap,
                            map<PolygonMeshModel::Vertex*, SimplificationByGarlan::Vertex>& vertexMap,
                            list<SimplificationByGarlan::Edge*>& edgesSortedByCost);

        void    collect_edges_influenced_by_edge_collapse(
                            SimplificationByGarlan::Edge* edgeToCollapse,
                            const map<PolygonMeshModel::Edge*, SimplificationByGarlan::Edge>& edgeMap,
                            list<SimplificationByGarlan::Edge*>& edgesInfluencedByEdgeCollapse);

        void    remove_influenced_edges_from_sorted_list(
                            SimplificationByGarlan::Edge* edgeToCollapse,
                            const map<PolygonMeshModel::Edge*, SimplificationByGarlan::Edge>& edgeMap,
                            list<SimplificationByGarlan::Edge*>& edgesSortedByCost, 
                            list<SimplificationByGarlan::Edge*>& edgesInfluencedByEdgeCollapse);

        void    update_influenced_edges(
                            const map<PolygonMeshModel::Vertex*, SimplificationByGarlan::Vertex>& vertexMap,
                            list<SimplificationByGarlan::Edge*>& edgesInfluencedByEdgeCollapse);

        void    insert_influenced_edges_into_sorted_list(
                            list<SimplificationByGarlan::Edge*>& edgesSortedByCost, 
                            list<SimplificationByGarlan::Edge*>& edgesInfluencedByEdgeCollapse);
};



///////////////////////////////////////////////////////////////////////////////
//
// class SimplificationByGarlan::Edge 
//
class SimplificationByGarlan::Edge
{
private:
    SimplificationByGarlan::Vertex* m_vertex[2];
    PolygonMeshModel::Edge*         m_pedge;
    double                          m_cost;
	rg_Point3D                      m_position4newVertex;

public:
    Edge();
    Edge(SimplificationByGarlan::Vertex* start, SimplificationByGarlan::Vertex* end, PolygonMeshModel::Edge* pedge);
    Edge(const Edge& edge);
    ~Edge();

    Edge& operator =(const Edge& edge);

    SimplificationByGarlan::Vertex*     start_vertex() const;
    SimplificationByGarlan::Vertex*     end_vertex() const;
    PolygonMeshModel::Edge*             pedge() const;
    double                              cost() const;
	rg_Point3D                          position_for_new_vertex() const;

    void    set_start_vertex(   SimplificationByGarlan::Vertex* vertex );
    void    set_end_vertex(     SimplificationByGarlan::Vertex* vertex );
    void    set_pedge(      PolygonMeshModel::Edge*         pedge);

    void    initialize();

    static  bool LessThan( Edge* e1, Edge* e2 );
};



///////////////////////////////////////////////////////////////////////////////
//
// class SimplificationByGarlan::Vertex 
//
class SimplificationByGarlan::Vertex
{
private:
    PolygonMeshModel::Vertex*   m_pvertex;
    double                      m_Q[10];    //  Quadric error metric

public:
    Vertex();
    Vertex(PolygonMeshModel::Vertex* pvertex);
    Vertex(const Vertex& vertex);
    ~Vertex();

    Vertex& operator =(const Vertex& vertex);

    PolygonMeshModel::Vertex*   pvertex() const;
    double                      Q(const int& i);

    void    set_pvertex(PolygonMeshModel::Vertex* pvertex);
    void    set_Q(const int& i, const double& value);
    void    set_Q(SimplificationByGarlan::Edge* edge);

    void    initialize();

};

inline  void                            SimplificationByGarlan::set_simplification(PolygonMeshModel* model, const double& limitOfQEM) {
    m_model = model;
    m_limitOfQEM = limitOfQEM; 
}

inline  SimplificationByGarlan::Vertex* SimplificationByGarlan::Edge::start_vertex() const              { return m_vertex[0]; }
inline  SimplificationByGarlan::Vertex* SimplificationByGarlan::Edge::end_vertex() const                { return m_vertex[1]; }
inline  PolygonMeshModel::Edge*         SimplificationByGarlan::Edge::pedge() const                     { return m_pedge; }
inline  double                          SimplificationByGarlan::Edge::cost() const                      { return m_cost; }
inline  rg_Point3D                      SimplificationByGarlan::Edge::position_for_new_vertex() const   { return m_position4newVertex; }
inline  void                            SimplificationByGarlan::Edge::set_start_vertex(SimplificationByGarlan::Vertex* vertex)  { m_vertex[0] = vertex; }
inline  void                            SimplificationByGarlan::Edge::set_end_vertex(SimplificationByGarlan::Vertex* vertex)    { m_vertex[1] = vertex; }
inline  void                            SimplificationByGarlan::Edge::set_pedge(PolygonMeshModel::Edge* pedge)              { m_pedge = pedge; }

inline  PolygonMeshModel::Vertex*       SimplificationByGarlan::Vertex::pvertex() const     { return m_pvertex; }
inline  double                          SimplificationByGarlan::Vertex::Q(const int& i)     { return m_Q[i]; }
inline  void                            SimplificationByGarlan::Vertex::set_pvertex(PolygonMeshModel::Vertex* pvertex)  { m_pvertex = pvertex; }
inline  void                            SimplificationByGarlan::Vertex::set_Q(const int& i, const double& value)        { m_Q[i] = value; }


#endif


