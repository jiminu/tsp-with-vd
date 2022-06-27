#ifndef _PATHAVOIDINGCIRCULAROBSTACLESIN2D_H
#define _PATHAVOIDINGCIRCULAROBSTACLESIN2D_H

#include "rg_Point2D.h"

#include <vector>
#include <list>
using namespace std;

class EdgeForSolutionPoolGraph;


class PathAvoidingCircularObstaclesIn2D
{
private:
    rg_Point2D                        m_start_point;
    rg_Point2D                        m_end_point;

    vector<EdgeForSolutionPoolGraph*> m_path;
    vector<bool>                      m_orientation;
    vector<double>                    m_distance_knot;


public:
    PathAvoidingCircularObstaclesIn2D();
    PathAvoidingCircularObstaclesIn2D(const PathAvoidingCircularObstaclesIn2D& path);
    ~PathAvoidingCircularObstaclesIn2D();

    void    clear();
    PathAvoidingCircularObstaclesIn2D& operator =(const PathAvoidingCircularObstaclesIn2D& path);

    void    construct(const rg_Point2D& startPoint, const rg_Point2D& endPoint, const list<EdgeForSolutionPoolGraph*>& path);

    double      length_of_path() const { return m_distance_knot.back(); }

    rg_Point2D  evaluate_point_on_path( const double& distanceFromStartPt) const;
    void        evaluate_points_on_path(vector<rg_Point2D>& pointsOnPath, const int& numPoints = 101) const;
    void        evaluate_points_on_path(vector<rg_Point2D>& pointsOnPath, const double& distInterval = 1.0) const;

private:
    unsigned int    find_index_of_edge_on_path_with_distance(const double& distanceFromStartPt) const;
};

#endif


