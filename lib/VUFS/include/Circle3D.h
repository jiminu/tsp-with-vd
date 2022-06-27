////////////////////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : Circle3D.h
//	  
//    DESCRIPTION : 
//           This is the implementation of the class Circle3D
//           This class defines a circle on a plane of 3D space.
//           The plane is defined by the center point of a circle and normal vector.
//
//    AUTHOR      : Ryu, Jooghyun
//    START DATE  : Nov. 10, 2008    
//
//           Copyright ¨Ï 2008 by Voronoi Diagram Research Center, Hanyang University
//
/////////////////////////////////////////////////////////////////////////////////////

#ifndef _CIRCLE3D_H
#define _CIRCLE3D_H

#include "rg_Point3D.h"
#include "Plane.h"
#include "rg_dList.h"

#include "LineSegment3D.h"


#include <list>
using namespace std;


class SquaredDistFuncLineSegCircle;


class Circle3D
{
protected:
    rg_Point3D m_center;
    rg_REAL    m_radius;

    rg_Point3D m_normal;

public:
    Circle3D();
    Circle3D(const rg_Point3D& center, const rg_REAL& radius, const rg_Point3D& normal);
    Circle3D(const Circle3D& circle);
    virtual ~Circle3D();

    inline rg_Point3D getCenter() const { return m_center; }
    inline rg_REAL    getRadius() const { return m_radius; }
    inline rg_Point3D getNormal() const { return m_normal; }

    void       setCenter(const rg_Point3D& center);
    void       setRadius(const rg_REAL& radius);
    void       setNormal(const rg_Point3D& normal);
    void       setCircle(const rg_Point3D& center, const rg_REAL& radius, const rg_Point3D& normal);

    Circle3D&  operator =(const Circle3D& circle);
    rg_BOOL    operator==(const Circle3D& circle) const;
    rg_BOOL    operator!=(const Circle3D& circle) const;

    rg_BOOL    isOnCircle(const rg_Point3D& point, const rg_REAL tolerance = resNeg6) const;
	rg_FLAG    isInsideCircle(const rg_Point3D& targetPt, const rg_REAL tolerance = resNeg6) const;

    void       reverse();

    inline Plane getPlaneContainingThisCircle() const { return Plane(m_normal, m_center); }

    rg_INT     intersect(const LineSegment3D& lineSegment, rg_Point3D* intersectPts) const;
    //  return: the number of intersections between a circle in 3D and plane 
    //          NOTE: If the intersection between a circle in 3D and plane is circle itself,
    //                the function return INT_MAX.
    rg_INT     intersect(const Plane& plane,     rg_Point3D* intersectPts, const rg_REAL tolerance = resNeg6) const;
    rg_INT     intersect(const Circle3D& circle, rg_Point3D* intersectPts, const rg_REAL tolerance = resNeg6) const;

    rg_REAL    evaluateCircumference() const;
    rg_REAL    evaluateAngle(const rg_Point3D& pt1OnCircle, const rg_Point3D& pt2OnCircle) const;
    rg_Point3D evaluateMiddlePoint(const rg_Point3D& point1, const rg_Point3D& point2) const;
    void       evaluatePointsOnCircle(const rg_REAL& unitLength, const rg_Point3D& startPoint,
                                        rg_dList<rg_Point3D>& pointsOnCircle) const;
    void       evaluatePointsOnCircle(list<rg_Point3D>& pointsOnCircle, const int& numPoints = 12) const;


	virtual rg_REAL computeMinDistanceFromLineSegment(const LineSegment3D& lineSegment, rg_Point3D& minPtOnCircle, rg_Point3D& minPtOnLineSeg, PosOnLineSegOrArc& pos) const;
	rg_REAL computeMinDistFromPoint(const rg_Point3D& targetPt, rg_Point3D& minPtOnCircle) const;
	rg_INT  computeIntersectionWith(const LineSegment3D& linesegment, rg_Point3D intersectPts[], rg_REAL param[]) const;

	rg_INT locatePointsWithLocalMinimaFromLineSegment(const LineSegment3D&  lineSegment,
		                                              rg_dList<rg_Point3D>& pointsOnCircle,
													  rg_dList<rg_Point3D>& pointsOnLineSegment,
													  rg_dList<rg_REAL>&    paramsOnLineSegment) const;

	rg_INT locatePointsWithLocalMinimaFromLineSegment(const LineSegment3D&  lineSegment,
		                                              rg_dList<rg_Point3D>& pointsOnCircle,
		                                              rg_dList<rg_Point3D>& pointsOnLineSegment,
													  rg_dList<rg_REAL>&    paramsOnLineSegment,
													  SquaredDistFuncLineSegCircle* sqdDistFunc) const;

	rg_INT locatePointsWithLocalMinimaFromLineSegment(const LineSegment3D&    lineSegment,
		                                              rg_dList<rg_Point3D>&   pointsOnCircle,
		                                              rg_dList<rg_Point3D>&   pointsOnLineSegment,
													  rg_dList<PosOnLineSegOrArc>& positionsOnCircle,
													  rg_dList<PosOnLineSegOrArc>& positionsOnLineSegment,
		                                              rg_dList<rg_REAL>&      paramsOnLineSegment,
													  SquaredDistFuncLineSegCircle* sqdDistFunc) const;

	rg_INT locatePointsWithLocalMinimaFromLineSegment(const LineSegment3D& lineSegment, 
		                                              rg_Point3D*&         pointsOnCircle,
											          rg_Point3D*&         pointsOnLineSegment,
											          rg_REAL*&            paramsOnLineSegment,
											          PosOnLineSegOrArc*&       positionsOnLineSegment,
											          rg_REAL*&            minDistances) const;


	rg_INT locatePointsWithLocalMinimaMaximaFromLineSegment(const LineSegment3D& lineSegment, 
												            rg_Point3D*& pointsOnCircle, 
												            rg_Point3D*& pointsOnLineSegment, 
												            rg_REAL*&    paramsOnLineSegment,
												            rg_REAL*&    minDistances) const;

	rg_INT locatePointsWithLocalMinimaMaximaFromLineSegment(const LineSegment3D&  lineSegment, 
															rg_dList<rg_Point3D>& pointsOnCircle, 
															rg_dList<rg_Point3D>& pointsOnLineSegment, 
															rg_dList<rg_REAL>&    paramsOnLineSegment,
												            rg_dList<rg_REAL>&    minDistances) const;


};

#endif


