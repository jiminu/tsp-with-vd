#ifndef _RG_TETRAHEDRON_H
#define _RG_TETRAHEDRON_H

#include "rg_Point3D.h"
#include "LineSegment3D.h"
#include "rg_Triangle.h"
#include "Sphere.h"

#include <list>
using namespace std;



class rg_Tetrahedron
{
private:
    rg_Point3D m_point[4];
    //  m_point[1], m_point[2], and m_point[3] has CCW orientation viewing from m_point[0].

public:
    enum    Position { IN_LOWER_HALF_SPACE, ON_HYPERPLANE, IN_UPPER_HALF_SPACE };

public:
    rg_Tetrahedron();
    rg_Tetrahedron(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4);
    rg_Tetrahedron(const rg_Tetrahedron& tetrahedron);
    ~rg_Tetrahedron();

    rg_Point3D getPoint(const rg_INT& i) const;
    void       setPoint(const rg_INT& i, const rg_Point3D& pt);

    rg_Tetrahedron& operator =(const rg_Tetrahedron& tetrahedron);

    rg_BOOL     isPointInTetrahedron(const rg_Point3D& point) const;

    bool        isLineSegmentInTetrahedron(const LineSegment3D& lineSegment) const;
    bool        isTriangleInTetrahedron(const rg_Triangle& triangle) const;

    bool        doesIntersectWith(const LineSegment3D& lineSegment) const;
    bool        doesIntersectWith(const rg_Triangle& triangle) const;
    bool        doesIntersectWith(const Sphere& sphere) const;

    rg_Triangle triangle(const rg_INT& i) const;

    rg_REAL     computeSignedVolume() const;


    Position    positionOnHyperplane( const Plane& hyperplane, const rg_REAL& res=rg_MATH_RES );

    bool        hasIntersectionWithHyperplane( const Plane& hyperplane, const rg_REAL& res=rg_MATH_RES );
    bool        isInOpenLowerHalfSpace(        const Plane& hyperplane, const rg_REAL& res=rg_MATH_RES );
    bool        isInOpenUpperHalfSpace(        const Plane& hyperplane, const rg_REAL& res=rg_MATH_RES );
    bool        isInClosedLowerHalfSpace(      const Plane& hyperplane, const rg_REAL& res=rg_MATH_RES );
    bool        isInClosedUpperHalfSpace(      const Plane& hyperplane, const rg_REAL& res=rg_MATH_RES );

    bool        intersectWithPlane( const Plane& hyperplane, list<rg_Point3D>& polygon, const rg_REAL& res=rg_MATH_RES );


    static rg_REAL computeSignedVolume(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4);
    

private:
    rg_BOOL checkOrientation() const;
};

#endif


