//********************************************************************
//
//	  FILENAME    : rg_Line2D.cpp
//	  
//    DESCRIPoint2DION : 
//           This consists of the implementation of class rg_Line2D.      
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 7 Apr 1998    
//
//    HISTORY:
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************
#include "rg_Line2D.h"
#include "rg_RelativeOp.h"
#include "rg_GeoFunc.h"
#include "rg_TMatrix2D.h"

rg_Line2D::rg_Line2D()
{
}

rg_Line2D::rg_Line2D(const rg_Point2D &sPt, const rg_Point2D &ePt)
{
	m_StartPoint = sPt;
	m_EndPoint   = ePt;
}

rg_Line2D::rg_Line2D(const rg_Line2D &line)
{
	m_StartPoint = line.m_StartPoint;
	m_EndPoint   = line.m_EndPoint;
}

rg_Line2D::~rg_Line2D()
{
}

////////////////////////////////////////////////////////////////////////
// Access elements



/*
rg_REAL rg_Line2D::getTangent() const
{
	return (endP.getY()-startP.getY())/(endP.getX()-startP.getX());
}
*/
          

bool rg_Line2D::does_contain_in_line_segment(const rg_Point2D& point) const
{
    bool    doeslineSegmentContainPoint = false;

    if (point == m_StartPoint) {
        doeslineSegmentContainPoint = true;
    }
    else if (point == m_EndPoint) {
        doeslineSegmentContainPoint = true;
    }
    else {
        if (rg_EQ(m_StartPoint.getX(), m_EndPoint.getX())) {
            double t = (point.getY() - m_StartPoint.getY()) / (m_EndPoint.getY() - m_StartPoint.getY());

            if ((t > 0.0) && (t < 1.0)) {
                doeslineSegmentContainPoint = true;
            }
        }
        else if (rg_EQ(m_StartPoint.getY(), m_EndPoint.getY())) {
            double t = (point.getX() - m_StartPoint.getX()) / (m_EndPoint.getX() - m_StartPoint.getX());

            if ((t > 0.0) && (t < 1.0)) {
                doeslineSegmentContainPoint = true;
            }
        }
        else {
            double t[2];
            t[0] = (point.getX() - m_StartPoint.getX()) / (m_EndPoint.getX() - m_StartPoint.getX());
            t[1] = (point.getY() - m_StartPoint.getY()) / (m_EndPoint.getY() - m_StartPoint.getY());

            if (rg_EQ(t[0], t[1]) && (t[0] > 0.0) && (t[0] < 1.0)) {
                doeslineSegmentContainPoint = true;
            }
        }
    }

    return doeslineSegmentContainPoint;
}


////////////////////////////////////////////////////////////////////////
// Operations

void rg_Line2D::rotate_this_line_about_point(const rg_REAL& angle, const rg_Point2D & point)
{
	rg_TMatrix2D rotationMat;
	rotationMat.rotate(angle, point);

	m_StartPoint = rotationMat * m_StartPoint;
	m_EndPoint   = rotationMat * m_EndPoint;
}

//void rg_Line2D::get_coefficients_of_implicit_form_of_line_equation(rg_REAL& coefficientOfX, rg_REAL& coefficientOfY, rg_REAL& coefficientOfConstant) const
//{
//    rg_REAL x1 = sp.getX();
//    rg_REAL y1 = sp.getY();
//    rg_REAL x2 = ep.getX();
//    rg_REAL y2 = ep.getY();
//
//    coefficientOfX        = y2-y1;
//    coefficientOfY        = x1-x2;
//    coefficientOfConstant = x2*y1 - x1*y2;
//}

rg_Point2D rg_Line2D::project(const rg_Point2D& givenPoint, rg_REAL& parameterOfFoot) const
{
    rg_Point2D vecStart2End   = m_EndPoint - m_StartPoint;
    rg_Point2D vecStart2Given = givenPoint - m_StartPoint;
    parameterOfFoot           = vecStart2Given.operator%(vecStart2End);
   
    double denominator = vecStart2End.magnitudeSquare();
    parameterOfFoot    = parameterOfFoot / denominator;

    rg_Point2D footOnEntireLine = m_StartPoint + parameterOfFoot * vecStart2End;

    return footOnEntireLine;
}


rg_Point2D rg_Line2D::compute_footprint_of_point_onto_line_segment(const rg_Point2D& givenPoint, rg_REAL& parameterOfFoot) const
{
    rg_Point2D footOnLineSegment = compute_perpendicular_footprint_of_point_onto_entire_line(givenPoint, parameterOfFoot);

    if (parameterOfFoot < 0.0)
    {
        parameterOfFoot = 0.0;
        footOnLineSegment = m_StartPoint;
    }
    else
    {
        if (parameterOfFoot >= 1.0)
        {
            parameterOfFoot = 1.0;
            footOnLineSegment = m_EndPoint;
        }
    }

    return footOnLineSegment;
}

bool rg_Line2D::compute_footprint_of_point_onto_line_segment(const rg_Point2D& givenPoint, rg_Point2D& footOnLineSegment) const
{
    rg_REAL parameterOfFoot;
    footOnLineSegment = compute_footprint_of_point_onto_line_segment(givenPoint, parameterOfFoot);
    return true;
}

rg_Point2D rg_Line2D::compute_perpendicular_footprint_of_point_onto_entire_line(const rg_Point2D& givenPoint, rg_REAL& parameterOfFoot) const
{
    //rg_Point2D vecStart2End = ep - sp;
    //rg_Point2D vecStart2Given = givenPoint - sp;
    //parameterOfFoot = vecStart2Given.operator%(vecStart2End);
    //double denominator = vecStart2End.magnitudeSquare();
    //parameterOfFoot = parameterOfFoot / denominator;

    //rg_Point2D footOnEntireLine = sp + parameterOfFoot * vecStart2End;

    //return footOnEntireLine;

    return project(givenPoint, parameterOfFoot);
}


bool rg_Line2D::compute_perpendicular_footprint_of_point_onto_entire_line(const rg_Point2D& givenPoint, rg_Point2D& footOnEntireLine) const
{
    rg_REAL parameterOfFoot = DBL_MAX;
    footOnEntireLine = compute_perpendicular_footprint_of_point_onto_entire_line(givenPoint, parameterOfFoot);
    return true;
}


rg_Point2D rg_Line2D::compute_intersection_with_line(const rg_Line2D& line, bool& bTwoLinesAreParallel) const
{
    return rg_GeoFunc::compute_intersection_between_two_lines(m_StartPoint, m_EndPoint, line.m_StartPoint, line.m_EndPoint, bTwoLinesAreParallel);
}


rg_Point2D rg_Line2D::compute_intersection_btw_two_line_segments(const rg_Line2D& lineSegment, bool& doesIntersect, IntersectionPointPosition& position) const
{
    doesIntersect = false;

    bool bTwoLinesAreParallel;
    double parameter[2];
    rg_Point2D intersectionPt = rg_GeoFunc::compute_intersection_between_two_lines(
                                                    m_StartPoint, m_EndPoint, 
                                                    lineSegment.m_StartPoint, lineSegment.m_EndPoint, 
                                                    parameter[0], parameter[1], bTwoLinesAreParallel);

    if (bTwoLinesAreParallel) 
    {
        if (is_overlapped(lineSegment)) 
        {
            position      = ON_A_WHOLE_LINESEGMENT;
            doesIntersect = true;
        }
        else
        {
            position      = NONE_INTERSECTION;
            doesIntersect = false;
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
                position = NOT_ON_LINESEGMENT_BUT_ON_EXTENDED_LINE;
            }
        }
        else {
            position = NOT_ON_LINESEGMENT_BUT_ON_EXTENDED_LINE;
        }

        doesIntersect = (position == NOT_ON_LINESEGMENT_BUT_ON_EXTENDED_LINE) ? false : true;
    }

    return intersectionPt;
}


rg_Point2D rg_Line2D::compute_tangent_vector_at_this_point_for_marching_direction(const rg_Point2D& givenPoint, const bool& bLeftFace, const rg_Point2D& marchingPoint)
{
    rg_Point2D tangentVectorOnLine = evaluateVector();
    rg_Point2D guideVec = marchingPoint - givenPoint;
    //rg_Point2D tangentVectorOnLine;

    switch (bLeftFace)
    {
    case true:
    {
        // cross product cannot be zero.
        if (rg_NEG(guideVec.operator*(tangentVectorOnLine)))
            tangentVectorOnLine = -1.0 * tangentVectorOnLine;
    }
    break;
    case false:
    {
        if (rg_POS(guideVec.operator*(tangentVectorOnLine)))
            tangentVectorOnLine = -1.0 * tangentVectorOnLine;
    }
    break;
    }

    return tangentVectorOnLine;
}



bool rg_Line2D::is_overlapped(const rg_Line2D& lineSegment) const /** * \brief convert a line passing through two points into the implicit form ( \f$Ax + By + C = 0 (A^2 + B^2 = 1)\f$ ). */ 
{
    bool overlapped = false;

    if (is_parallel_to(lineSegment) && does_contain(lineSegment.m_StartPoint)) {
        double param[2][2] = { { get_parameter_of_point_on_line_segment(lineSegment.m_StartPoint), get_parameter_of_point_on_line_segment(lineSegment.m_EndPoint) },
                               { lineSegment.get_parameter_of_point_on_line_segment(m_StartPoint), lineSegment.get_parameter_of_point_on_line_segment(m_EndPoint) } };

        if (    (param[0][0] >= 0.0 && param[0][0] <= 1.0) 
             || (param[0][1] >= 0.0 && param[0][1] <= 1.0)
             || (param[1][0] >= 0.0 && param[1][0] <= 1.0)
             || (param[1][1] >= 0.0 && param[1][1] <= 1.0) ) 
        {
            overlapped = true;
        }
    }

    return overlapped ;
}


bool rg_Line2D::does_intersect_as_line_segment(const rg_Line2D& lineSegment) const
{
    bool doesIntersect                 = false;
    IntersectionPointPosition position = NONE_INTERSECTION;

    compute_intersection_btw_two_line_segments(lineSegment, doesIntersect, position);
    
    return doesIntersect;
}



bool rg_Line2D::Is_perpendicular_footprint_of_point_on_line_segment( const rg_Point2D& givenPoint )
{
    rg_Point2D footprint;
    compute_perpendicular_footprint_of_point_onto_entire_line( givenPoint, footprint );
    return does_contain_in_line_segment( footprint );
}


////////////////////////////////////////////////////////////////////////
// Operator Overloading
rg_FLAG rg_Line2D::operator==(const rg_Line2D &line) const
{
    if ( (m_StartPoint == line.m_StartPoint && m_EndPoint == line.m_EndPoint)
         || (m_StartPoint == line.m_EndPoint && m_EndPoint == line.m_StartPoint) )
    {
        return rg_TRUE;
    }
    else 
        return rg_FALSE;
}

rg_Line2D& rg_Line2D::operator =(const rg_Line2D& line)
{
    if ( this == &line )
        return *this;

    m_StartPoint = line.m_StartPoint;
    m_EndPoint = line.m_EndPoint;

    return *this;
}


bool  rg_Line2D::convertIntoImplicitForm(double& A, double& B, double& C)
{
    double x1 = m_StartPoint.getX();
    double y1 = m_StartPoint.getY();
    double x2 = m_EndPoint.getX();
    double y2 = m_EndPoint.getY();

    //  make the line equation passing through two points.
    //
    //  Two point equation of a line
    //   y - y1     y2 - y1
    //  -------- = ---------
    //   x - x1     x2 - x1
    //
    //  General form of the line equation
    //  Ax + By + C = 0 (A^2 + B^2 = 1)

    bool isConversionComplete = true;

    if (m_StartPoint.isEqual(m_EndPoint)) {
        isConversionComplete = false;
    }
    else {
        if (x1 == x2) {
            A = 1.0;
            B = 0.0;
            C = -x1;
        }
        else if (y1 == y2) {
            A = 0.0;
            B = 1.0;
            C = -y1;
        }
        else {
            double Aprime = (y2 - y1) / (x2 - x1);
            double Bprime = -1.0;
            double Cprime = y1 - (Aprime*x1);
            double normalization = sqrt(Aprime*Aprime + Bprime*Bprime);
            A = Aprime / normalization;
            B = Bprime / normalization;
            C = Cprime / normalization;
        }

        isConversionComplete = true;
    }

    return isConversionComplete;
}
