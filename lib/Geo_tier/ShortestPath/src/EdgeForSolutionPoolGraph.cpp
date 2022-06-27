#include "EdgeForSolutionPoolGraph.h"
#include "VertexForSolutionPoolGraph.h"

EdgeForSolutionPoolGraph::EdgeForSolutionPoolGraph()
{
    m_StartVtx = NULL;
    m_EndVtx   = NULL;
    m_EdgeType = LINE_SEGMENT_SPG;
    m_ImplicitEquation = rg_ImplicitEquation();
    m_EdgeLength = 0.0;
}

EdgeForSolutionPoolGraph::EdgeForSolutionPoolGraph(const EdgeForSolutionPoolGraph& edgeForSPG)
{
    copy_from(edgeForSPG); 
}

EdgeForSolutionPoolGraph::~EdgeForSolutionPoolGraph()
{

}

EdgeForSolutionPoolGraph& EdgeForSolutionPoolGraph::operator=(const EdgeForSolutionPoolGraph& edgeForSPG)
{
    if (this == &edgeForSPG)
    {
        return *this;
    }

    clear();
    copy_from(edgeForSPG);

    return *this;
}

void EdgeForSolutionPoolGraph::clear()
{
    m_StartVtx = NULL;
    m_EndVtx   = NULL;
    m_EdgeType = LINE_SEGMENT_SPG;
    m_ImplicitEquation = rg_ImplicitEquation();
    m_EdgeLength = 0.0;
}

void EdgeForSolutionPoolGraph::copy_from(const EdgeForSolutionPoolGraph& edgeForSPG)
{
    m_StartVtx   = edgeForSPG.m_StartVtx;
    m_EndVtx     = edgeForSPG.m_EndVtx;
    m_EdgeType   = edgeForSPG.m_EdgeType;
    m_ImplicitEquation = edgeForSPG.m_ImplicitEquation;
    m_EdgeLength = edgeForSPG.m_EdgeLength;
}



Arc2D EdgeForSolutionPoolGraph::get_arc() const
{
    Arc2D arcOfEdge;

    rg_Point2D startPt = get_start_vtx()->get_tangent_point_to_disk();
    rg_Point2D endPt   = get_end_vtx()->get_tangent_point_to_disk();

    switch (m_EdgeType)
    {
    case ARC_SPG:
    {
        // x2 + y2 + Ax + By + C = 0
        // (x - a)2 + (y - b)2 = r2
        // a = -A/2, b = -B/2, r = sqrt( (A2+B2)/4 - C )
        double A = m_ImplicitEquation.getCoeff(1, 0);
        double B = m_ImplicitEquation.getCoeff(0, 1);
        double C = m_ImplicitEquation.getCoeff(0, 0);

        double a = -A / 2.0;
        double b = -B / 2.0;
        double r = sqrt((pow(A, 2) + pow(B, 2)) / 4.0 - C);

        rg_Point2D center(a, b);
        rg_Point2D vecCenterToStart = startPt - center;
        rg_Point2D vecCenterToEnd   = endPt   - center;

        double theta = angleFromVec1toVec2(vecCenterToStart, vecCenterToEnd);

        arcOfEdge = Arc2D(rg_Circle2D(a, b, r), startPt, endPt);

        /*
        if (theta > rg_PI)
        {
            arcOfEdge = Arc2D(rg_Circle2D(a, b, r), endPt, startPt);
        }
        else
        {
            arcOfEdge = Arc2D(rg_Circle2D(a, b, r), startPt, endPt);
        }
        */
        break;
    }

    case LINE_SEGMENT_SPG:
    {
        arcOfEdge = Arc2D(rg_Circle2D(DBL_MAX, DBL_MAX, DBL_MAX), startPt, endPt);
        break;
    }

    default:
        break;
    }

    return arcOfEdge;
}



VertexForSolutionPoolGraph* EdgeForSolutionPoolGraph::get_opposite_vtx(VertexForSolutionPoolGraph* currVtx) const
{
    if (currVtx == m_StartVtx)
    {
        return m_EndVtx;
    }
    else if (currVtx == m_EndVtx)
    {
        return m_StartVtx;
    }
    else
    {
        return NULL;
    }
}
