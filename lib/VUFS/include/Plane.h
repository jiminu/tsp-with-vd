#ifndef _PLANE_H
#define _PLANE_H

#include "rg_Point3D.h"
#include "Line3D.h"
#include "LineSegment3D.h"



//  Coordinate representation of the plane in 3D space
//  ax + by + cz = d 
//    m_normal (a, b, c) is a unit vector.
//    d is a distance between the plane and origin.

// History ---------------------------------------------
// Inserted by Joonghyun
// void        projectLineSegmentOnPlane()
// rg_FLAG computeIntersectionWithLineSegment()
// rg_FLAG isParallelTo()
// rg_FLAG isCoincidentWith()
// Plane getReversedPlane()

class Plane
{
private:
	rg_Point3D m_normal;
    rg_REAL    m_distanceFromOrigin;

public:
    Plane();
    Plane(const rg_Point3D& normal, const rg_REAL& distanceFromOrigin);
    Plane(const rg_Point3D& normal, const rg_Point3D& passingPoint);
    Plane(const rg_Point3D& passingPoint1, const rg_Point3D& passingPoint2, const rg_Point3D& passingPoint3);
    Plane(const Plane& plane);
    ~Plane();

    rg_Point3D  getNormal() const;
    rg_REAL     getDistanceFromOrigin() const;
	Plane       getReversedPlane() const;

    void        setNormal(const rg_Point3D& normal);
    void        setDistanceFromOrigin(const rg_REAL& distanceFromOrigin);
    void        setPlane(const rg_Point3D& normal, const rg_REAL& distanceFromOrigin);
    void        reverseNormal();

    void        definePlaneByNormalAndPassingPoint(const rg_Point3D& normal, const rg_Point3D& passingPoint);
    void        definePlaneByThreePoints(const rg_Point3D& passingPoint1, const rg_Point3D& passingPoint2, const rg_Point3D& passingPoint3);

    rg_REAL     distanceFromPoint(const rg_Point3D& point) const;

	rg_REAL     computeAngleInCCW(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3) const;
    rg_REAL     computeProjectedAngleInCCW(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3) const;
    rg_REAL     computeRelativeAngleInCCW(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3);

    rg_Point3D  reflectPoint(const rg_Point3D& point) const;
    rg_Point3D  projectPointOnPlane(const rg_Point3D& point) const;
	rg_FLAG     computeIntersectionWithPlane(const Plane& plane, Line3D& line) const;

    Plane&  operator =(const Plane& plane);
    rg_BOOL operator==(const Plane& plane) const;
    rg_BOOL operator!=(const Plane& plane) const;

	rg_FLAG isOnNormalSide(const rg_Point3D& point) const;	
	rg_FLAG isOnOppositeNormalSide(const rg_Point3D& point) const;
	rg_FLAG isOnThisPlane(const rg_Point3D& point) const;
	rg_FLAG isOnNormalSideOrOnThisPlane(const rg_Point3D& point) const;
	rg_FLAG isOnOppositeNormalSideOrOnThisPlane(const rg_Point3D& point) const;

    rg_BOOL     isThereIntersectionWith(const Plane& plane) const;
    rg_BOOL     intersect(const Plane& plane, Line3D& line) const;
    rg_INT      intersect(const LineSegment3D& lineSegment, rg_Point3D& intersectionPoint) const;
    rg_INT      intersectLineSegment(const rg_Point3D& pt1, const rg_Point3D& pt2, rg_Point3D& intersectionPoint) const;
	rg_Point3D  intersectLineSegment(const rg_Point3D& pt1, const rg_Point3D& pt2) const;

	void        projectLineSegmentOnPlane(const LineSegment3D& targetLineSegment, LineSegment3D& projection) const;
	rg_FLAG computeIntersectionWithLineSegment(const rg_Point3D& pt1, const rg_Point3D& pt2, rg_Point3D& intersecPt) const;	
	rg_FLAG computeIntersectionWithLineSegment(const LineSegment3D& linesegment, rg_Point3D& intersecPt, rg_REAL& intersecParam, PosOnLineSegOrArc& pos) const;	
	rg_FLAG computeIntersectionWithLineSegment(const LineSegment3D& linesegment, rg_Point3D& intersecPt, rg_REAL& intersecParam) const;	
	rg_FLAG computeIntersectionWithLineSegment(const LineSegment3D& linesegment, rg_REAL& intersecParam) const;	

	rg_FLAG isParallelTo(const Plane& plane) const;
	rg_FLAG isCoincidentWith(const Plane& plane) const;	

};

#endif
