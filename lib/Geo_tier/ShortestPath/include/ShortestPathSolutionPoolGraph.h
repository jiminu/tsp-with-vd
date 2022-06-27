#ifndef _SHORTEST_PATH_SOLUTION_POOL_GRAPH_
#define _SHORTEST_PATH_SOLUTION_POOL_GRAPH_


#include "VertexForSolutionPoolGraph.h"
#include "EdgeForSolutionPoolGraph.h"
#include "rg_Circle2D.H"
#include <unordered_map>
using namespace std;

class ShortestPathSolutionPoolGraph
{
private:
    list<VertexForSolutionPoolGraph> m_Vertices;
    list<EdgeForSolutionPoolGraph>   m_Edges;

    VertexForSolutionPoolGraph* m_StarVtx;
    VertexForSolutionPoolGraph* m_EndVtx;

public:
    //constructor
    ShortestPathSolutionPoolGraph();
    ShortestPathSolutionPoolGraph(const ShortestPathSolutionPoolGraph& SPSPG);

    //destructor;
    ~ShortestPathSolutionPoolGraph();

    //equal operation
    ShortestPathSolutionPoolGraph& operator=(const ShortestPathSolutionPoolGraph& SPSPG);

    //setter
    void set_vertices(const list<VertexForSolutionPoolGraph>& vertices);
    void set_edges(const list<EdgeForSolutionPoolGraph>& edges);
    inline void set_start_vtx(VertexForSolutionPoolGraph* startVtx) { m_StarVtx = startVtx; };
    inline void set_end_vtx(VertexForSolutionPoolGraph* endVtx)     { m_EndVtx = endVtx; };
    
    //getter
    void get_vertices(list<VertexForSolutionPoolGraph*>& vertices) const;
    void get_edges(list<EdgeForSolutionPoolGraph*>& edges) const;

    inline VertexForSolutionPoolGraph* get_start_vtx() const { return m_StarVtx; };
    inline VertexForSolutionPoolGraph* get_end_vtx() const   { return m_EndVtx; };

    //function
    VertexForSolutionPoolGraph* create_vertex(const VertexForSolutionPoolGraph& vtx);
    EdgeForSolutionPoolGraph*   create_edge(const EdgeForSolutionPoolGraph& edge);

private:
    void copy_from(const ShortestPathSolutionPoolGraph& SPSPG);
    void make_link_btw_Edges(const list<EdgeForSolutionPoolGraph*>& prevEdges, unordered_map<EdgeForSolutionPoolGraph*, EdgeForSolutionPoolGraph*>& EdgeLinkerFromPrevToCurr);
    void make_link_btw_Verices(const list<VertexForSolutionPoolGraph*>& prevVertices, unordered_map<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>& VertexLinkerFromPrevToCurr);
    void modify_Edge_topology(const unordered_map<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>& VertexLinkerFromPrevToCurr);
    void modify_Vertex_topology(const unordered_map<EdgeForSolutionPoolGraph*, EdgeForSolutionPoolGraph*>& EdgeLinkerFromPrevToCurr) ;
    void modify_graph_start_N_end_vertex(VertexForSolutionPoolGraph* oldStartVtx, VertexForSolutionPoolGraph* oldEndVtx, const unordered_map<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>& VertexLinkerFromPrevToCurr);

    void clear();
};


#endif