#include "VertexForSolutionPoolGraph.h"

VertexForSolutionPoolGraph::VertexForSolutionPoolGraph()
{
    m_AccumulateLengthFromSourceVtx = DBL_MAX;
    m_Angle = -1.0;
    m_PrevEdge = NULL;
}


VertexForSolutionPoolGraph::VertexForSolutionPoolGraph(const VertexForSolutionPoolGraph& vertexForSPG)
{
    copy_from(vertexForSPG);
}


VertexForSolutionPoolGraph::~VertexForSolutionPoolGraph()
{
}


VertexForSolutionPoolGraph& VertexForSolutionPoolGraph::operator=(const VertexForSolutionPoolGraph& vertexForSPG)
{
    if (this == &vertexForSPG)
    {
        return *this;
    }

    clear();
    copy_from(vertexForSPG);

    return *this;
}


void VertexForSolutionPoolGraph::set_arc_edges(const list<EdgeForSolutionPoolGraph*>& arcEdges)
{
    m_ArcEdges.clear();
    m_ArcEdges = arcEdges;
}



void VertexForSolutionPoolGraph::set_tangent_line_segment_edges(const list<EdgeForSolutionPoolGraph*>& lineEdges)
{
    m_TangentLineSegmentEdges.clear();
    m_TangentLineSegmentEdges = lineEdges;
}



void VertexForSolutionPoolGraph::copy_from(const VertexForSolutionPoolGraph& vertexForSPG)
{
    m_ArcEdges.clear();
    m_TangentLineSegmentEdges.clear();

    m_TangentLineSegmentEdges       = vertexForSPG.m_TangentLineSegmentEdges;
    m_ArcEdges                      = vertexForSPG.m_ArcEdges;
    m_TangentPointToDisk            = vertexForSPG.m_TangentPointToDisk;
    m_AccumulateLengthFromSourceVtx = vertexForSPG.m_AccumulateLengthFromSourceVtx;
    m_Angle                         = vertexForSPG.m_Angle;
    m_PrevEdge                      = vertexForSPG.m_PrevEdge;
}


void VertexForSolutionPoolGraph::erase_line_edge(EdgeForSolutionPoolGraph * lineEdge)
{
    list<EdgeForSolutionPoolGraph*>::iterator it_TargetEdge = find(m_TangentLineSegmentEdges.begin(), m_TangentLineSegmentEdges.end(), lineEdge);

    if (it_TargetEdge != m_TangentLineSegmentEdges.end())
    {
        m_TangentLineSegmentEdges.erase(it_TargetEdge);
    }
}


void VertexForSolutionPoolGraph::erase_arc_edge(EdgeForSolutionPoolGraph * arcEdge)
{
    list<EdgeForSolutionPoolGraph*>::iterator it_TargetEdge = find(m_ArcEdges.begin(), m_ArcEdges.end(), arcEdge);

    if (it_TargetEdge != m_ArcEdges.end())
    {
        m_ArcEdges.erase(it_TargetEdge);
    }
}


void VertexForSolutionPoolGraph::clear()
{
    m_ArcEdges.clear();
    m_TangentLineSegmentEdges.clear();

    m_TangentPointToDisk            = rg_Point2D(DBL_MAX,DBL_MAX);
    m_AccumulateLengthFromSourceVtx = DBL_MAX;
    m_Angle = -1.0;
    m_PrevEdge = NULL;
}



bool VertexForSolutionPoolGraph::are_two_same_point(VertexForSolutionPoolGraph* first, VertexForSolutionPoolGraph* second)
{
    if (first->get_tangent_point_to_disk() == second->get_tangent_point_to_disk())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void VertexForSolutionPoolGraph::get_all_edges(list<EdgeForSolutionPoolGraph*>& edges) const
{
    for (list<EdgeForSolutionPoolGraph*>::const_iterator it_TangentLine = m_TangentLineSegmentEdges.begin();
        it_TangentLine != m_TangentLineSegmentEdges.end();
        it_TangentLine++)
    {

        edges.push_back(*it_TangentLine);
    }

    for (list<EdgeForSolutionPoolGraph*>::const_iterator it_Arc = m_ArcEdges.begin();
        it_Arc != m_ArcEdges.end();
        it_Arc++)
    {

        edges.push_back(*it_Arc);
    }
}
