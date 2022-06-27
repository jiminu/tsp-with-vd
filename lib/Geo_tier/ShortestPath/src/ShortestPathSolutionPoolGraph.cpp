#include "ShortestPathSolutionPoolGraph.h"

ShortestPathSolutionPoolGraph::ShortestPathSolutionPoolGraph()
{
    m_StarVtx = NULL;
    m_EndVtx  = NULL;
}


ShortestPathSolutionPoolGraph::ShortestPathSolutionPoolGraph(const ShortestPathSolutionPoolGraph& SPSPG)
{
    copy_from(SPSPG);
}


ShortestPathSolutionPoolGraph::~ShortestPathSolutionPoolGraph()
{

}


ShortestPathSolutionPoolGraph& ShortestPathSolutionPoolGraph::operator=(const ShortestPathSolutionPoolGraph& SPSPG)
{
    if (this == &SPSPG)
    {
        return *this;
    }

    clear();
    copy_from(SPSPG);

    return *this;
}



void ShortestPathSolutionPoolGraph::copy_from(const ShortestPathSolutionPoolGraph& SPSPG)
{
    clear();

    unordered_map<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*> VertexLinkerFromPrevToCurr;
    unordered_map<EdgeForSolutionPoolGraph*, EdgeForSolutionPoolGraph*> EdgeLinkerFromPrevToCurr;

    list<VertexForSolutionPoolGraph*> prevVertices;
    list<EdgeForSolutionPoolGraph*> prevEdges;

    SPSPG.get_vertices(prevVertices);
    SPSPG.get_edges(prevEdges);

    //1. make link for vertex    
    make_link_btw_Verices(prevVertices, VertexLinkerFromPrevToCurr);

    //2. make link for edge
    make_link_btw_Edges(prevEdges, EdgeLinkerFromPrevToCurr);

    //3. change vertex topology from prev to curr
    modify_Vertex_topology(EdgeLinkerFromPrevToCurr);

    //4. change edge topology from prev to curr
    modify_Edge_topology(VertexLinkerFromPrevToCurr);

    //5. change start and end vtx
    modify_graph_start_N_end_vertex(SPSPG.get_start_vtx(), SPSPG.get_end_vtx(), VertexLinkerFromPrevToCurr);
}

void ShortestPathSolutionPoolGraph::make_link_btw_Edges(const list<EdgeForSolutionPoolGraph*>& prevEdges, unordered_map<EdgeForSolutionPoolGraph*, EdgeForSolutionPoolGraph*>& EdgeLinkerFromPrevToCurr)
{
    for (list<EdgeForSolutionPoolGraph*>::const_iterator it_Edge = prevEdges.begin(); it_Edge != prevEdges.end(); it_Edge++)
    {
        EdgeForSolutionPoolGraph* prevEdge = *it_Edge;
        EdgeForSolutionPoolGraph* newEdge = create_edge(*prevEdge);
        
        EdgeLinkerFromPrevToCurr.insert(make_pair(prevEdge, newEdge));
    }
}

void ShortestPathSolutionPoolGraph::make_link_btw_Verices(const list<VertexForSolutionPoolGraph*>& prevVertices, unordered_map<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>& VertexLinkerFromPrevToCurr)
{
    for (list<VertexForSolutionPoolGraph*>::const_iterator it_Vertex = prevVertices.begin(); it_Vertex != prevVertices.end(); it_Vertex++)
    {
        VertexForSolutionPoolGraph* prevVertex = *it_Vertex;
        VertexForSolutionPoolGraph* newVertex   = create_vertex(*prevVertex);

        VertexLinkerFromPrevToCurr.insert(make_pair(prevVertex, newVertex));
    }
}

void ShortestPathSolutionPoolGraph::modify_Edge_topology(const unordered_map<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>& VertexLinkerFromPrevToCurr)
{
    for (list<EdgeForSolutionPoolGraph>::iterator it_Edge = m_Edges.begin(); it_Edge != m_Edges.end(); it_Edge++)
    {
        EdgeForSolutionPoolGraph* newEdge = &(*it_Edge);

        VertexForSolutionPoolGraph* newStartVtx = VertexLinkerFromPrevToCurr.at(newEdge->get_start_vtx());
        VertexForSolutionPoolGraph* newEndVtx   = VertexLinkerFromPrevToCurr.at(newEdge->get_end_vtx());

        newEdge->set_start_vtx(newStartVtx);
        newEdge->set_end_vtx(newEndVtx);
    }
}

void ShortestPathSolutionPoolGraph::modify_Vertex_topology(const unordered_map<EdgeForSolutionPoolGraph*, EdgeForSolutionPoolGraph*>& EdgeLinkerFromPrevToCurr)
{
    for (list<VertexForSolutionPoolGraph>::iterator it_Vertex = m_Vertices.begin(); it_Vertex != m_Vertices.end(); it_Vertex++)
    {
        VertexForSolutionPoolGraph* newVertex = &(*it_Vertex);

        list<EdgeForSolutionPoolGraph*> oldLineEdges;
        list<EdgeForSolutionPoolGraph*> oldArcEdges;

        newVertex->get_line_edges(oldLineEdges);
        newVertex->get_arc_edges(oldArcEdges);
        
        // 1. line edge
        list<EdgeForSolutionPoolGraph*> newLineEdges;
        
        for (list<EdgeForSolutionPoolGraph*>::const_iterator it_OldEdge = oldLineEdges.begin(); it_OldEdge != oldLineEdges.end(); it_OldEdge++)
        {
            EdgeForSolutionPoolGraph* oldEdge = *it_OldEdge;

            newLineEdges.push_back(EdgeLinkerFromPrevToCurr.at(oldEdge));
        }

        // 2. arc edge
        list<EdgeForSolutionPoolGraph*> newArcEdges;

        for (list<EdgeForSolutionPoolGraph*>::const_iterator it_OldEdge = oldArcEdges.begin(); it_OldEdge != oldArcEdges.end(); it_OldEdge++)
        {
            EdgeForSolutionPoolGraph* oldEdge = *it_OldEdge;

            newArcEdges.push_back(EdgeLinkerFromPrevToCurr.at(oldEdge));
        }

        // 3. prev edge
        EdgeForSolutionPoolGraph* newPrevEdge = NULL;

        if (newVertex->get_prev_edge() != NULL)
        {
            newPrevEdge = EdgeLinkerFromPrevToCurr.at(newVertex->get_prev_edge());
        }
        
        // 4. set new edges
        newVertex->set_tangent_line_segment_edges(newLineEdges);
        newVertex->set_arc_edges(newArcEdges);
        newVertex->set_prev_edge(newPrevEdge);
    }
}


void ShortestPathSolutionPoolGraph::modify_graph_start_N_end_vertex(VertexForSolutionPoolGraph* oldStartVtx, VertexForSolutionPoolGraph* oldEndVtx, const unordered_map<VertexForSolutionPoolGraph*, VertexForSolutionPoolGraph*>& VertexLinkerFromPrevToCurr)
{
    if (oldStartVtx != NULL)
    {
        m_StarVtx = VertexLinkerFromPrevToCurr.at(oldStartVtx);
    }

    if (oldEndVtx != NULL)
    {
        m_EndVtx = VertexLinkerFromPrevToCurr.at(oldEndVtx);
    }
}



void ShortestPathSolutionPoolGraph::clear()
{
    m_Vertices.clear();
    m_Edges.clear();

    m_StarVtx = NULL;
    m_EndVtx  = NULL;
}


void ShortestPathSolutionPoolGraph::set_vertices(const list<VertexForSolutionPoolGraph>& vertices)
{
    m_Vertices.clear();
    m_Vertices = vertices;
}


void ShortestPathSolutionPoolGraph::set_edges(const list<EdgeForSolutionPoolGraph>& edges)
{
    m_Edges.clear();
    m_Edges = edges;
}


void ShortestPathSolutionPoolGraph::get_vertices(list<VertexForSolutionPoolGraph*>& vertices) const
{
    for (list<VertexForSolutionPoolGraph>::const_iterator it_vtx = m_Vertices.begin(); it_vtx != m_Vertices.end(); it_vtx++)
    {
        VertexForSolutionPoolGraph* currVtx = const_cast<VertexForSolutionPoolGraph*>(&(*it_vtx));
        vertices.push_back(currVtx);
    }
}


void ShortestPathSolutionPoolGraph::get_edges(list<EdgeForSolutionPoolGraph*>& edges) const
{
    for (list<EdgeForSolutionPoolGraph>::const_iterator it_edge = m_Edges.begin(); it_edge != m_Edges.end(); it_edge++)
    {
        EdgeForSolutionPoolGraph* currEdge = const_cast<EdgeForSolutionPoolGraph*>(&(*it_edge));
        edges.push_back(currEdge);
    }
}


VertexForSolutionPoolGraph* ShortestPathSolutionPoolGraph::create_vertex(const VertexForSolutionPoolGraph& vtx)
{
    m_Vertices.push_back(vtx);
    return (&(*m_Vertices.rbegin()));
}


EdgeForSolutionPoolGraph* ShortestPathSolutionPoolGraph::create_edge(const EdgeForSolutionPoolGraph& edge)
{
    m_Edges.push_back(edge);
    return (&(*m_Edges.rbegin()));
}



