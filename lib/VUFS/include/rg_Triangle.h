#ifndef RG_TRIANGLE_H
#define RG_TRIANGLE_H

#include "rg_Point3D.h"
#include "LineSegment3D.h"
#include "Circle3D.h"
#include "Sphere.h"

#include <list>
using namespace std;



class rg_Triangle
{
protected:
    rg_Point3D m_point[3];

public:
    enum    Position { IN_LOWER_HALF_SPACE, ON_HYPERPLANE, IN_UPPER_HALF_SPACE };

public:
    rg_Triangle();
    rg_Triangle(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3);
    rg_Triangle(const rg_Triangle& triangle);
    ~rg_Triangle();

    inline rg_Point3D point(const rg_INT& i) const { return m_point[i]; }
    inline rg_Point3D getPoint(const rg_INT& i) const { return m_point[i]; }
    inline void       setPoint(const rg_INT& i, const rg_Point3D& point) { m_point[i] = point; }
    inline void       setTriangle(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3) { m_point[0] = pt1; m_point[1] = pt2; m_point[2] = pt3; }

    rg_Triangle& operator =(const rg_Triangle& triangle);
    inline  rg_Point3D   operator[](const int& i) const { return m_point[i]; }

    rg_BOOL     isPointOnTriangle(const rg_Point3D& point) const;
    rg_BOOL     isCircleInTriangle(const Circle3D& circle) const;

    bool        doesIntersectWith(const Sphere& sphere) const;

    double      area() const;
    rg_Point3D  computeNormalVector() const;

    rg_INT      intersect(const LineSegment3D& lineSegment, rg_Point3D& intersectionPoint) const;

    static double  area(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3);


    Position    positionOnHyperplane( const Plane& hyperplane, const rg_REAL& res=rg_MATH_RES );

    bool        intersectWithLowerHalfSpace( const Plane& hyperplane, list<rg_Point3D>& polygon, const rg_REAL& res=rg_MATH_RES );
    bool        intersectWithLowerHalfSpace( const Plane& hyperplane, list<rg_Point3D>& polygon, LineSegment3D& lineSegToIntersect, const rg_REAL& res=rg_MATH_RES );

    bool        intersectWithPlane( const Plane& plane, LineSegment3D& lineSegToIntersect, const rg_REAL& res=rg_MATH_RES );
    
    
    rg_REAL     computeMinDistFromPoint(const rg_Point3D& targetPt) const;
};

#endif 

