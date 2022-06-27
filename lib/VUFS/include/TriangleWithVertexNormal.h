#ifndef TRIANGLE_WITH_VERTEX_NORMAL_H
#define TRIANGLE_WITH_VERTEX_NORMAL_H

#include "rg_Triangle.h"

#include <vector>
using namespace std;


class TriangleWithVertexNormal : public rg_Triangle
{
private:
    rg_Point3D m_normal[3];

public:
    TriangleWithVertexNormal();
    TriangleWithVertexNormal(const rg_Point3D& point1, const rg_Point3D& point2, const rg_Point3D& point3,
                             const rg_Point3D& normal1, const rg_Point3D& normal2, const rg_Point3D& normal3);
    TriangleWithVertexNormal(const TriangleWithVertexNormal& triangle);
    ~TriangleWithVertexNormal();

    TriangleWithVertexNormal& operator =(const TriangleWithVertexNormal& triangle);

    rg_Point3D  normal(const int& i) const;
    void        set_normal(const int& i, const rg_Point3D& normal);
    void        set_triangle(const rg_Point3D& point1, const rg_Point3D& point2, const rg_Point3D& point3,
                             const rg_Point3D& normal1, const rg_Point3D& normal2, const rg_Point3D& normal3);

    bool        intersectWithLowerHalfSpace(const Plane& hyperplane, vector<TriangleWithVertexNormal>& region, const rg_REAL& res = rg_MATH_RES);


};

#endif


