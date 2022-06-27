#include "TriangleWithVertexNormal.h"


TriangleWithVertexNormal::TriangleWithVertexNormal()
{
}



TriangleWithVertexNormal::TriangleWithVertexNormal(const rg_Point3D& point1, const rg_Point3D& point2, const rg_Point3D& point3,
                                const rg_Point3D& normal1, const rg_Point3D& normal2, const rg_Point3D& normal3)
: rg_Triangle(point1, point2, point3)
{
    m_normal[0] = normal1;
    m_normal[1] = normal2;
    m_normal[2] = normal3;
}



TriangleWithVertexNormal::TriangleWithVertexNormal(const TriangleWithVertexNormal& triangle)
: rg_Triangle(triangle)
{
    m_normal[0] = triangle.m_normal[0];
    m_normal[1] = triangle.m_normal[1];
    m_normal[2] = triangle.m_normal[2];
}



TriangleWithVertexNormal::~TriangleWithVertexNormal()
{
}



TriangleWithVertexNormal& TriangleWithVertexNormal::operator =(const TriangleWithVertexNormal& triangle)
{
    if (this != &triangle) {
        rg_Triangle::operator=(triangle);
        m_normal[0] = triangle.m_normal[0];
        m_normal[1] = triangle.m_normal[1];
        m_normal[2] = triangle.m_normal[2];
    }

    return *this;
}



rg_Point3D  TriangleWithVertexNormal::normal(const int& i) const
{
    return m_normal[i];
}



void        TriangleWithVertexNormal::set_normal(const int& i, const rg_Point3D& normal)
{
    m_normal[i] = normal;
}



void        TriangleWithVertexNormal::set_triangle(const rg_Point3D& point1, const rg_Point3D& point2, const rg_Point3D& point3,
    const rg_Point3D& normal1, const rg_Point3D& normal2, const rg_Point3D& normal3)
{
    setTriangle(point1, point2, point3);
    m_normal[0] = normal1;
    m_normal[1] = normal2;
    m_normal[2] = normal3;
}



bool        TriangleWithVertexNormal::intersectWithLowerHalfSpace(const Plane& hyperplane, vector<TriangleWithVertexNormal>& region, const rg_REAL& res )
{
    int     numPtsInOpenLowerHalfSpace = 0;
    int     numPtsInOpenUpperHalfSpace = 0;
    int     numPtsOnHyperPlane = 0;

    bool    onHyperplane[3] = { false, false, false };
    bool    onUpperHalfSpace[3] = { false, false, false };
    bool    onLowerHalfSpace[3] = { false, false, false };

    int i = 0;
    for (i = 0; i<3; ++i) {
        double distance = hyperplane.distanceFromPoint(m_point[i]);
        if (rg_ZERO(distance, res)) {
            ++numPtsOnHyperPlane;
            onHyperplane[i] = true;
        }
        else if (rg_POS(distance, res)) {
            ++numPtsInOpenUpperHalfSpace;
            onUpperHalfSpace[i] = true;
        }
        else {
            ++numPtsInOpenLowerHalfSpace;
            onLowerHalfSpace[i] = true;
        }
    }


    bool    hasIntersectionWithLowerHalfSpace = false;
    if (numPtsInOpenUpperHalfSpace + numPtsOnHyperPlane == 3) {
        hasIntersectionWithLowerHalfSpace = false;
    }
    else if (numPtsInOpenLowerHalfSpace + numPtsOnHyperPlane == 3) {
        hasIntersectionWithLowerHalfSpace = true;

        region.push_back(*this);
    }
    else if (numPtsInOpenUpperHalfSpace == 1 && numPtsInOpenLowerHalfSpace == 2) {
        hasIntersectionWithLowerHalfSpace = true;

        if (onUpperHalfSpace[0]) {
            rg_Point3D intersection[2];
            hyperplane.intersect(LineSegment3D(m_point[0], m_point[1]), intersection[0]);
            hyperplane.intersect(LineSegment3D(m_point[2], m_point[0]), intersection[1]);
            double param[2] = { m_point[0].distance(intersection[0]) / m_point[0].distance(m_point[1]),
                                m_point[0].distance(intersection[1]) / m_point[0].distance(m_point[2]) };
            rg_Point3D normalAtIntersection[2] = { (1.0 - param[0])*m_normal[0] + param[0] * m_normal[1], 
                                                   (1.0 - param[1])*m_normal[0] + param[1] * m_normal[2] };

            // intersection[0] -> m_point[1] -> m_point[2] -> intersection[1]
            region.resize(2);
            region[0].set_triangle(intersection[0], m_point[1], m_point[2], normalAtIntersection[0], m_normal[1], m_normal[2]);
            region[1].set_triangle(intersection[0], m_point[2], intersection[1], normalAtIntersection[0], m_normal[2], normalAtIntersection[1]);
        }
        else if (onUpperHalfSpace[1]) {
            rg_Point3D intersection[2];
            hyperplane.intersect(LineSegment3D(m_point[1], m_point[2]), intersection[0]);
            hyperplane.intersect(LineSegment3D(m_point[0], m_point[1]), intersection[1]);
            double param[2] = { m_point[1].distance(intersection[0]) / m_point[1].distance(m_point[2]),
                                m_point[1].distance(intersection[1]) / m_point[1].distance(m_point[0]) };
            rg_Point3D normalAtIntersection[2] = { (1.0 - param[0])*m_normal[1] + param[0] * m_normal[2],
                                                   (1.0 - param[1])*m_normal[1] + param[1] * m_normal[0] };

            // intersection[0] -> m_point[2] -> m_point[0] -> intersection[1]
            region.resize(2);
            region[0].set_triangle(intersection[0], m_point[2], m_point[0], normalAtIntersection[0], m_normal[2], m_normal[0]);
            region[1].set_triangle(intersection[0], m_point[0], intersection[1], normalAtIntersection[0], m_normal[0], normalAtIntersection[1]);
        }
        else { //if ( onUpperHalfSpace[2] ) {
            rg_Point3D intersection[2];
            hyperplane.intersect(LineSegment3D(m_point[2], m_point[0]), intersection[0]);
            hyperplane.intersect(LineSegment3D(m_point[1], m_point[2]), intersection[1]);
            double param[2] = { m_point[2].distance(intersection[0]) / m_point[2].distance(m_point[0]),
                                m_point[2].distance(intersection[1]) / m_point[2].distance(m_point[1]) };
            rg_Point3D normalAtIntersection[2] = { (1.0 - param[0])*m_normal[2] + param[0] * m_normal[0],
                                                   (1.0 - param[1])*m_normal[2] + param[1] * m_normal[1] };

            // intersection[0] -> m_point[0] -> m_point[1] -> intersection[1]
            region.resize(2);
            region[0].set_triangle(intersection[0], m_point[0], m_point[1], normalAtIntersection[0], m_normal[0], m_normal[1]);
            region[1].set_triangle(intersection[0], m_point[1], intersection[1], normalAtIntersection[0], m_normal[1], normalAtIntersection[1]);
        }
    }
    else if (numPtsInOpenUpperHalfSpace == 2 && numPtsInOpenLowerHalfSpace == 1) {
        hasIntersectionWithLowerHalfSpace = true;

        if (onLowerHalfSpace[0]) {
            rg_Point3D intersection[2];
            hyperplane.intersect(LineSegment3D(m_point[0], m_point[1]), intersection[0]);
            hyperplane.intersect(LineSegment3D(m_point[2], m_point[0]), intersection[1]);
            double param[2] = { m_point[0].distance(intersection[0]) / m_point[0].distance(m_point[1]),
                                m_point[0].distance(intersection[1]) / m_point[0].distance(m_point[2]) };
            rg_Point3D normalAtIntersection[2] = { (1.0 - param[0])*m_normal[0] + param[0] * m_normal[1],
                                                   (1.0 - param[1])*m_normal[0] + param[1] * m_normal[2] };

            // m_point[0] -> intersection[0] -> intersection[1]
            region.resize(1);
            region[0].set_triangle(m_point[0], intersection[0], intersection[1], m_normal[0], normalAtIntersection[0], normalAtIntersection[1]);
        }
        else if (onLowerHalfSpace[1]) {
            rg_Point3D intersection[2];
            hyperplane.intersect(LineSegment3D(m_point[1], m_point[2]), intersection[0]);
            hyperplane.intersect(LineSegment3D(m_point[0], m_point[1]), intersection[1]);
            double param[2] = { m_point[1].distance(intersection[0]) / m_point[1].distance(m_point[2]),
                                m_point[1].distance(intersection[1]) / m_point[1].distance(m_point[0]) };
            rg_Point3D normalAtIntersection[2] = { (1.0 - param[0])*m_normal[1] + param[0] * m_normal[2],
                                                   (1.0 - param[1])*m_normal[1] + param[1] * m_normal[0] };

            // m_point[1] -> intersection[0] -> intersection[1]
            region.resize(1);
            region[0].set_triangle(m_point[1], intersection[0], intersection[1], m_normal[1], normalAtIntersection[0], normalAtIntersection[1]);
        }
        else { //if ( onLowerHalfSpace[2] ) {
            rg_Point3D intersection[2];
            hyperplane.intersect(LineSegment3D(m_point[2], m_point[0]), intersection[0]);
            hyperplane.intersect(LineSegment3D(m_point[1], m_point[2]), intersection[1]);
            double param[2] = { m_point[2].distance(intersection[0]) / m_point[2].distance(m_point[0]),
                                m_point[2].distance(intersection[1]) / m_point[2].distance(m_point[1]) };
            rg_Point3D normalAtIntersection[2] = { (1.0 - param[0])*m_normal[2] + param[0] * m_normal[0],
                                                   (1.0 - param[1])*m_normal[2] + param[1] * m_normal[1] };

            // m_point[2] -> intersection[0] -> intersection[1]
            region.resize(1);
            region[0].set_triangle(m_point[2], intersection[0], intersection[1], m_normal[2], normalAtIntersection[0], normalAtIntersection[1]);
        }
    }
    else {
        hasIntersectionWithLowerHalfSpace = true;
    }

    return hasIntersectionWithLowerHalfSpace;
}




