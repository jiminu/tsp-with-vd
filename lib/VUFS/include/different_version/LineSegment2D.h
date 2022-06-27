#ifndef LINESEGMENT2D_H
#define LINESEGMENT2D_H

#include "rg_Point2D.h"

class LineSegment2D
{
private:
    rg_Point2D  m_start_point;
    rg_Point2D  m_end_point;

public:
    enum PointPos { NOT_ON_LINESEGMENT, AT_START_POINT, AT_MID_POINT, AT_END_POINT, AT_LINESEGMENT };
    enum IntersectionPos { NON_INTERSECTION, INTERSECTION_AT_START_POINT, INTERSECTION_AT_MID_POINT, INTERSECTION_AT_END_POINT };
public:
    LineSegment2D();
    LineSegment2D(const rg_Point2D& startPoint, const rg_Point2D endPoint);
    LineSegment2D(const LineSegment2D& lineSegment);
    ~LineSegment2D();

    LineSegment2D& operator =(const LineSegment2D& lineSegment);

    rg_Point2D  start_point() const;
    rg_Point2D  end_point() const;

    void        set_start_point(const rg_Point2D& startPoint);
    void        set_end_point(const rg_Point2D& endPoint);

    double      squared_length() const;
    double      length() const;

    double      signed_distance(const rg_Point2D& point) const;
    double      evaluate_parameter(const rg_Point2D& point) const;

    rg_Point2D  direction() const;
    rg_Point2D  normal() const;

    rg_Point2D  evaluate_point(const double& param) const;
    rg_Point2D  project_point(const rg_Point2D& point) const;

    bool        is_parallel(const LineSegment2D& lineSegment) const;
    bool        is_overlapped(const LineSegment2D& lineSegment) const;


    bool        does_contain_in_corresponding_line(const rg_Point2D& point) const;
    bool        does_contain(const rg_Point2D& point) const;
    bool        does_contain(const rg_Point2D& point, PointPos& position) const;
    bool        does_intersect(const LineSegment2D& lineSegment) const;
    bool        intersect(const LineSegment2D& lineSegment, rg_Point2D& intersection, PointPos& position) const;
};



inline  rg_Point2D  LineSegment2D::start_point() const { return m_start_point; }
inline  rg_Point2D  LineSegment2D::end_point() const { return m_end_point; }
inline  void        LineSegment2D::set_start_point(const rg_Point2D& startPoint) { m_start_point = startPoint; }
inline  void        LineSegment2D::set_end_point(const rg_Point2D& endPoint) { m_end_point = endPoint; }

#endif // !
