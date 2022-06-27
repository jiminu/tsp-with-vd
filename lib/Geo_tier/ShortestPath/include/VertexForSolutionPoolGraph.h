#ifndef _VERTEX_FOR_SOLUTION_POOL_GRAPH_
#define _VERTEX_FOR_SOLUTION_POOL_GRAPH_

#include "rg_Point2D.h"
#include "rg_Circle2D.h"
#include <list>
using namespace std;

class EdgeForSolutionPoolGraph;

class VertexForSolutionPoolGraph
{
private:
    list<EdgeForSolutionPoolGraph*> m_TangentLineSegmentEdges;
    list<EdgeForSolutionPoolGraph*> m_ArcEdges;

    rg_Point2D                      m_TangentPointToDisk;
    double                          m_AccumulateLengthFromSourceVtx;

    double                          m_Angle;
    EdgeForSolutionPoolGraph*       m_PrevEdge;

    void*                           m_userData;

public:
    class AngleLess {
    public:
        bool    operator()(VertexForSolutionPoolGraph* vertex1, VertexForSolutionPoolGraph* vertex2) {
            return (vertex1->get_angle() < vertex2->get_angle()) ? true : false;
        }
    };

    //constructor
    VertexForSolutionPoolGraph();
    VertexForSolutionPoolGraph(const VertexForSolutionPoolGraph& vertexForSPG);

    //destructor;
    ~VertexForSolutionPoolGraph();

    //equal operation
    VertexForSolutionPoolGraph& operator=(const VertexForSolutionPoolGraph& vertexForSPG);


    //setter
    inline void set_tangent_point_to_disk(const rg_Point2D& tangentPoint) { m_TangentPointToDisk = tangentPoint; };
    inline void set_accumulate_length_from_source_vtx(const double& accumulateLength) { m_AccumulateLengthFromSourceVtx = accumulateLength; };
    inline void set_angle(const double& angle) { m_Angle = angle; };
    inline void set_prev_edge(EdgeForSolutionPoolGraph* edge)  { m_PrevEdge = edge; };

    void set_tangent_line_segment_edges(const list<EdgeForSolutionPoolGraph*>& lineEdges);
    void set_arc_edges(const list<EdgeForSolutionPoolGraph*>& arcEdges);

    void    set_user_data(void* userData) { m_userData = userData; }
    void    accumulate_length_from_source_vtx(const double& accumulateLength) { m_AccumulateLengthFromSourceVtx += accumulateLength; };

    //getter
    inline const rg_Point2D&  get_tangent_point_to_disk()   const { return m_TangentPointToDisk; };
    inline double get_accumulate_length_from_source_vtx() const { return m_AccumulateLengthFromSourceVtx; };
    inline double get_angle() const { return m_Angle; };

    inline void   get_line_edges(list<EdgeForSolutionPoolGraph*>& lineEdges) const { lineEdges = m_TangentLineSegmentEdges; };
    inline void   get_arc_edges(list<EdgeForSolutionPoolGraph*>& arcEdges) const { arcEdges = m_ArcEdges; };
    void          get_all_edges(list<EdgeForSolutionPoolGraph*>& edges) const;
    inline EdgeForSolutionPoolGraph* get_prev_edge()    const { return m_PrevEdge; };

    void*         user_data() const { return m_userData; }

    //function
    inline void add_line_edge(EdgeForSolutionPoolGraph* lineEdge) { m_TangentLineSegmentEdges.push_back(lineEdge); };
    inline void add_arc_edge(EdgeForSolutionPoolGraph* arcEdge) { m_ArcEdges.push_back(arcEdge); };
    void erase_line_edge(EdgeForSolutionPoolGraph* lineEdge) ;
    void erase_arc_edge(EdgeForSolutionPoolGraph* arcEdge) ;

    void clear();

    //operator
    static bool are_two_same_point(VertexForSolutionPoolGraph* first, VertexForSolutionPoolGraph* second);

private:
    void copy_from(const VertexForSolutionPoolGraph& vertexForSPG);

};




#endif