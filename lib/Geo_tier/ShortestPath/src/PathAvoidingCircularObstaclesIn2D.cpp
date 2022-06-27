#include "PathAvoidingCircularObstaclesIn2D.h"
#include "EdgeForSolutionPoolGraph.h"
#include "VertexForSolutionPoolGraph.h"



PathAvoidingCircularObstaclesIn2D::PathAvoidingCircularObstaclesIn2D()
{
}



PathAvoidingCircularObstaclesIn2D::PathAvoidingCircularObstaclesIn2D(const PathAvoidingCircularObstaclesIn2D& path)
: m_start_point(path.m_start_point)
, m_end_point  (path.m_end_point)
, m_path       (path.m_path)
, m_orientation(path.m_orientation)
, m_distance_knot(path.m_distance_knot)
{
}



PathAvoidingCircularObstaclesIn2D::~PathAvoidingCircularObstaclesIn2D()
{
}



void    PathAvoidingCircularObstaclesIn2D::clear()
{
    m_start_point = rg_Point2D();
    m_end_point   = rg_Point2D();
    m_path.clear();
    m_orientation.clear();
    m_distance_knot.clear();

}



PathAvoidingCircularObstaclesIn2D& PathAvoidingCircularObstaclesIn2D::operator =(const PathAvoidingCircularObstaclesIn2D& path)
{
    if (this != &path) {
        m_start_point   = path.m_start_point;
        m_end_point     = path.m_end_point;
        m_path          = path.m_path;
        m_orientation   = path.m_orientation;
        m_distance_knot = path.m_distance_knot;
    }

    return *this;
}



void    PathAvoidingCircularObstaclesIn2D::construct(const rg_Point2D& startPoint, const rg_Point2D& endPoint, const list<EdgeForSolutionPoolGraph*>& path)
{
    clear();
    m_start_point = startPoint;
    m_end_point   = endPoint;
    m_path.insert(m_path.begin(), path.begin(), path.end());

    int numEdgeInPath = m_path.size();
    m_orientation.resize(numEdgeInPath);
    m_distance_knot.resize(numEdgeInPath+1);


    int i = 0;
    //  determine the orientation of edges in solution pool
    VertexForSolutionPoolGraph* prevVtx = m_path[0]->get_start_vtx();
    rg_Point2D startingPt = prevVtx->get_tangent_point_to_disk();
    if (!(startingPt == m_start_point)) {
        prevVtx = m_path[0]->get_end_vtx();
        startingPt = prevVtx->get_tangent_point_to_disk();
    }

    for (i = 0; i < numEdgeInPath; ++i) {
        m_orientation[i] = (prevVtx == m_path[i]->get_start_vtx()) ? true : false;

        prevVtx = m_path[i]->get_opposite_vtx(prevVtx);
    }

    //  determine the m_distance_knot values 
    m_distance_knot[0] = 0.0;
    for (i = 0; i < numEdgeInPath; ++i) {
        double edgeLength = m_path[i]->get_edge_length();

        m_distance_knot[i+1] = m_distance_knot[i] + edgeLength;
    }

}



rg_Point2D  PathAvoidingCircularObstaclesIn2D::evaluate_point_on_path(const double& distanceFromStartPt) const
{
    rg_Point2D pointOnPath;
    if (rg_LE(distanceFromStartPt, 0.0)) {
        pointOnPath = m_start_point;
    }
    else if (rg_GE(distanceFromStartPt, m_distance_knot.back()) ) {
        pointOnPath = m_end_point;
    }
    else {
        int edgeIndex = find_index_of_edge_on_path_with_distance(distanceFromStartPt);

        EdgeForSolutionPoolGraph* edgeOnPath = m_path[edgeIndex];
        if (edgeOnPath->get_edge_type() == LINE_SEGMENT_SPG) {
            rg_Point2D  start;
            rg_Point2D  end;
            if (m_orientation[edgeIndex] == true) {
                start = edgeOnPath->get_start_vtx()->get_tangent_point_to_disk();
                end   = edgeOnPath->get_end_vtx()->get_tangent_point_to_disk();
            }
            else {
                start = edgeOnPath->get_end_vtx()->get_tangent_point_to_disk();
                end   = edgeOnPath->get_start_vtx()->get_tangent_point_to_disk();
            }
            rg_Point2D  direction = (end - start).getUnitVector();

            double      difference = distanceFromStartPt - m_distance_knot[edgeIndex];

            pointOnPath = start + (difference*direction);
        }
        else if (edgeOnPath->get_edge_type() == ARC_SPG) {
            Arc2D arc = edgeOnPath->get_arc();

            double  lengthDifference = 0.0;
            if (m_orientation[edgeIndex] == true) {
                lengthDifference = distanceFromStartPt - m_distance_knot[edgeIndex];
            }
            else {
                lengthDifference = m_distance_knot[edgeIndex+1] - distanceFromStartPt;
            }
            double angle = lengthDifference / arc.getRadius();

            rg_Point2D uVector  = (arc.getStartPoint() - arc.getCenterPt()).getUnitVector();
            rg_Point2D vVector(-uVector.getY(), uVector.getX());
            pointOnPath = arc.getCenterPt() + arc.getRadius()*(cos(angle)*uVector + sin(angle)*vVector);
        }
    }

    return pointOnPath;
}



void        PathAvoidingCircularObstaclesIn2D::evaluate_points_on_path(vector<rg_Point2D>& pointsOnPath, const int& numPoints) const
{
    double pathLength = length_of_path();
    double distInterval = pathLength / (numPoints - 1);
    
    evaluate_points_on_path(pointsOnPath, distInterval);
}



void        PathAvoidingCircularObstaclesIn2D::evaluate_points_on_path(vector<rg_Point2D>& pointsOnPath, const double& distInterval) const
{
    list<rg_Point2D> ptsOnPath;

    double pathLength = length_of_path();
    for (double dist = 0.0; dist < pathLength; dist += distInterval) {
        rg_Point2D point = evaluate_point_on_path(dist);
        ptsOnPath.push_back(point);
    }
    ptsOnPath.push_back(m_end_point);

    pointsOnPath.insert(pointsOnPath.begin(), ptsOnPath.begin(), ptsOnPath.end());
}



unsigned int    PathAvoidingCircularObstaclesIn2D::find_index_of_edge_on_path_with_distance(const double& distanceFromStartPt) const
{
    int numEdgeInPath = m_path.size();
    unsigned int edgeIndex = 0;
    if (distanceFromStartPt < m_distance_knot[0]) {
        edgeIndex = 0;
    }
    else if (distanceFromStartPt >= m_distance_knot.back()) {
        edgeIndex = numEdgeInPath-1;
    }
    else {
        for (int i = 0; i < numEdgeInPath; ++i) {
            if (distanceFromStartPt >= m_distance_knot[i] && distanceFromStartPt < m_distance_knot[i + 1]) {
                edgeIndex = i;
                break;
            }
        }
    }

    return edgeIndex;
}
