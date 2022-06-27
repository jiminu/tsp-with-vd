#include "LineSegment2D.h"
#include "rg_RelativeOp.h"
#include "rg_GeoFunc.h"

#include <cmath>
using namespace std;


LineSegment2D::LineSegment2D()
{
}



LineSegment2D::LineSegment2D(const rg_Point2D& startPoint, const rg_Point2D endPoint)
:  m_start_point(startPoint)
,  m_end_point(endPoint)
{
}



LineSegment2D::LineSegment2D(const LineSegment2D& lineSegment)
: m_start_point(lineSegment.m_start_point)
, m_end_point(lineSegment.m_end_point)
{
}



LineSegment2D::~LineSegment2D()
{
}



LineSegment2D& LineSegment2D::operator =(const LineSegment2D& lineSegment)
{
    if (this != &lineSegment) {
        m_start_point = lineSegment.m_start_point;
        m_end_point = lineSegment.m_end_point;
    }

    return *this;
}


double      LineSegment2D::squared_length() const
{
    return (   (m_end_point.getX() - m_start_point.getX())*(m_end_point.getX() - m_start_point.getX()) 
             + (m_end_point.getY() - m_start_point.getY())*(m_end_point.getY() - m_start_point.getY()));
}



double      LineSegment2D::length() const
{   
    return sqrt(squared_length());
}


double      LineSegment2D::signed_distance(const rg_Point2D& point) const
{
    rg_Point2D   spPtVector(point - m_start_point);
    rg_Point2D   epPtVector(point - m_end_point);
    rg_Point2D   lineVector(m_end_point - m_start_point);

    rg_REAL spPtLen = spPtVector.magnitude();
    rg_REAL epPtLen = epPtVector.magnitude();
    rg_REAL lineLen = lineVector.magnitude();

    if (rg_GE(epPtLen*epPtLen, (spPtLen*spPtLen) + (lineLen*lineLen))) {
        return spPtLen;
    }
    else if (rg_GE((spPtLen*spPtLen), (epPtLen*epPtLen) + (lineLen*lineLen))) {
        return epPtLen;
    }
    else{
        return (lineVector*spPtVector/lineLen);
    }
}


double      LineSegment2D::evaluate_parameter(const rg_Point2D& point) const
{
    double    parameter = DBL_MAX;

    if (point == m_start_point) {
        parameter = 0.0;
    }
    else if (point == m_end_point) {
        parameter = 1.0;
    }
    else {
        if (rg_EQ(m_start_point.getX(), m_end_point.getX())) {
            parameter = (point.getY() - m_start_point.getY()) / (m_end_point.getY() - m_start_point.getY());
        }
        else if (rg_EQ(m_start_point.getY(), m_end_point.getY())) {
            parameter = (point.getX() - m_start_point.getX()) / (m_end_point.getX() - m_start_point.getX());

        }
        else {
            double t[2];
            t[0] = (point.getX() - m_start_point.getX()) / (m_end_point.getX() - m_start_point.getX());
            t[1] = (point.getY() - m_start_point.getY()) / (m_end_point.getY() - m_start_point.getY());

            if (rg_EQ(t[0], t[1]) ) {
                parameter = t[0];
            }
        }
    }

    return parameter;
}


rg_Point2D  LineSegment2D::direction() const
{
    rg_Point2D direct = (m_end_point - m_start_point).getUnitVector();
    return direct;
}



rg_Point2D  LineSegment2D::normal() const
{
    rg_Point2D direct = direction();
    return rg_Point2D(-direct.getY(), direct.getX());
}



rg_Point2D  LineSegment2D::evaluate_point(const double& param) const
{
    return m_start_point + param * (m_end_point - m_start_point);
}



rg_Point2D  LineSegment2D::project_point(const rg_Point2D& point) const
{
    double dist = signed_distance(point);
    rg_Point2D normal = this->normal();

    return point - (dist*normal);
}



bool        LineSegment2D::is_parallel(const LineSegment2D& lineSegment) const
{
    rg_Point2D vector[2] = { (m_end_point - m_start_point).getUnitVector(), 
                             (lineSegment.m_end_point - lineSegment.m_start_point).getUnitVector() };

    return (rg_ZERO(vector[0] * vector[1])) ? true : false;
}


bool        LineSegment2D::is_overlapped(const LineSegment2D& lineSegment) const
{
    bool overlapped = false;
    if (is_parallel(lineSegment) && does_contain_in_corresponding_line(lineSegment.m_start_point)) {
        double param[2][2] = { { evaluate_parameter(lineSegment.m_start_point), evaluate_parameter(lineSegment.m_end_point) },
                               { lineSegment.evaluate_parameter(m_start_point), lineSegment.evaluate_parameter(m_end_point) } };

        if (    (param[0][0] >= 0.0 && param[0][0] <= 1.0) 
             || (param[0][1] >= 0.0 && param[0][1] <= 1.0)
             || (param[1][0] >= 0.0 && param[1][0] <= 1.0)
             || (param[1][1] >= 0.0 && param[1][1] <= 1.0) ) {
            overlapped = true;
        }
    }

    return overlapped ;
}



bool        LineSegment2D::does_contain_in_corresponding_line(const rg_Point2D& point) const
{
    rg_Line2D line(m_start_point, m_end_point);
    double distance = line.signed_distance(point);

    return (rg_ZERO(distance)) ? true : false;
}



bool        LineSegment2D::does_contain(const rg_Point2D& point) const
{
    PointPos position;
    return does_contain( point, position);
}



bool        LineSegment2D::does_contain(const rg_Point2D& point, LineSegment2D::PointPos& position) const
{   
    position = NOT_ON_LINESEGMENT;
    bool    doeslineSegmentContainPoint = false;

    if (point == m_start_point) {
        position = AT_START_POINT;
        doeslineSegmentContainPoint = true;
    }
    else if (point == m_end_point) {
        position = AT_END_POINT;
        doeslineSegmentContainPoint = true;
    }
    else {
        if (rg_EQ(m_start_point.getX(), m_end_point.getX())) {
            double t = (point.getY() - m_start_point.getY()) / (m_end_point.getY() - m_start_point.getY());

            if ((t > 0.0) && (t < 1.0)) {
                position = AT_MID_POINT;
                doeslineSegmentContainPoint = true;
            }
        }
        else if (rg_EQ(m_start_point.getY(), m_end_point.getY())) {
            double t = (point.getX() - m_start_point.getX()) / (m_end_point.getX() - m_start_point.getX());

            if ((t > 0.0) && (t < 1.0)) {
                position = AT_MID_POINT;
                doeslineSegmentContainPoint = true;
            }
        }
        else {
            double t[2];
            t[0] = (point.getX() - m_start_point.getX()) / (m_end_point.getX() - m_start_point.getX());
            t[1] = (point.getY() - m_start_point.getY()) / (m_end_point.getY() - m_start_point.getY());

            if (rg_EQ(t[0], t[1]) && (t[0] > 0.0) && (t[0] < 1.0)) {
                position = AT_MID_POINT;
                doeslineSegmentContainPoint = true;
            }
        }
    }

    return doeslineSegmentContainPoint;
}



bool        LineSegment2D::does_intersect(const LineSegment2D& lineSegment) const
{
    rg_Point2D intersection;
    LineSegment2D::PointPos position;
    return intersect(lineSegment, intersection, position);
}



bool        LineSegment2D::intersect(const LineSegment2D& lineSegment, rg_Point2D& intersection, LineSegment2D::PointPos& position) const
{
    bool    doesIntersect = false;
    bool bTwoLinesAreParallel;
    double parameter[2];
    intersection = rg_GeoFunc::compute_intersection_between_two_lines(
                                                    m_start_point, m_end_point, 
                                                    lineSegment.m_start_point, lineSegment.m_end_point, 
                                                    parameter[0], parameter[1], bTwoLinesAreParallel);

    if (bTwoLinesAreParallel) {
        if (is_overlapped(lineSegment)) {
            position = AT_LINESEGMENT;
            doesIntersect = true;
        }
    }
    else {
        if (rg_GE(parameter[1], 0.0) && rg_LE(parameter[1], 1.0)) {
            if (rg_EQ(parameter[0], 0.0)) {
                position = AT_START_POINT;
            }
            else if (rg_EQ(parameter[0], 1.0)) {
                position = AT_END_POINT;
            }
            else if ((parameter[0] > 0.0) && (parameter[0] < 1.0)) {
                position = AT_MID_POINT;
            }
            else {
                position = NOT_ON_LINESEGMENT;
            }
        }
        else {
            position = NOT_ON_LINESEGMENT;
        }

        doesIntersect = (position == NOT_ON_LINESEGMENT) ? false : true;
    }

    return doesIntersect;
}



