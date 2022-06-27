#ifndef _EDGE_FOR_SOLUTION_POOL_GRAPH_
#define _EDGE_FOR_SOLUTION_POOL_GRAPH_

#include "rg_Point2D.h"
#include "rg_Circle2D.h"
#include "rg_ImplicitEquation.h"
#include "Arc2D.h"
#include <list>
using namespace std;

class VertexForSolutionPoolGraph;

enum EdgeTypeOfSPG {LINE_SEGMENT_SPG, ARC_SPG};

class EdgeForSolutionPoolGraph
{
private:
    VertexForSolutionPoolGraph*     m_StartVtx;
    VertexForSolutionPoolGraph*     m_EndVtx;

    EdgeTypeOfSPG       m_EdgeType;
    rg_ImplicitEquation m_ImplicitEquation;
    double              m_EdgeLength;

    void*                           m_userData;

public:
    //constructor
    EdgeForSolutionPoolGraph();
    EdgeForSolutionPoolGraph(const EdgeForSolutionPoolGraph& edgeForSPG);

    //destructor;
    ~EdgeForSolutionPoolGraph();

    //equal operation
    EdgeForSolutionPoolGraph& operator=(const EdgeForSolutionPoolGraph& edgeForSPG);
    
    //setter
    inline void set_start_vtx(VertexForSolutionPoolGraph* startVtx) { m_StartVtx = startVtx; };
    inline void set_end_vtx(VertexForSolutionPoolGraph* endVtx)     { m_EndVtx   = endVtx; };

    inline void set_edge_type(const EdgeTypeOfSPG& edgeType)                      { m_EdgeType = edgeType; };
    inline void set_implicit_equation(const rg_ImplicitEquation& implictEquation) { m_ImplicitEquation = implictEquation; };
    inline void set_edge_length(const double& edgeLength)                         { m_EdgeLength = edgeLength; };

    void    set_user_data(void* userData) { m_userData = userData; }

    //getter
    inline VertexForSolutionPoolGraph*  get_start_vtx() const { return m_StartVtx; };
    inline VertexForSolutionPoolGraph*  get_end_vtx()   const { return m_EndVtx; };
    inline EdgeTypeOfSPG                get_edge_type() const { return m_EdgeType;};
    inline const rg_ImplicitEquation&   get_implicit_equation() const { return m_ImplicitEquation; };
    inline double get_edge_length() const { return m_EdgeLength; };

    void*         user_data() const { return m_userData; }

    VertexForSolutionPoolGraph* get_opposite_vtx(VertexForSolutionPoolGraph* currVtx) const;
    Arc2D get_arc() const;

    //function
    void clear();

private:
    void copy_from(const EdgeForSolutionPoolGraph& edgeForSPG);

};


#endif