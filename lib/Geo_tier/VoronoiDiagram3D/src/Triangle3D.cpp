#include "Triangle3D.h"
using namespace V::GeometryTier;

#include <float.h>


Triangle3D::Triangle3D()
{
}



Triangle3D::Triangle3D(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3)
{
    m_point[0] = pt1;
    m_point[1] = pt2;
    m_point[2] = pt3;
}



Triangle3D::Triangle3D(const Triangle3D& triangle)
{
    m_point[0] = triangle.m_point[0];
    m_point[1] = triangle.m_point[1];
    m_point[2] = triangle.m_point[2];
}



Triangle3D::~Triangle3D()
{

}




rg_Point3D  Triangle3D::getPoint(const rg_INT& i) const
{
    if ( i<0 || i>3 ) {
        return rg_Point3D(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    }
    else {
        return m_point[i];
    }
}



rg_Point3D* Triangle3D::getPoints()
{
    return m_point;
}




void        Triangle3D::setPoint(const rg_INT& i, const rg_Point3D& pt)
{
    if ( i<0 || i>3 ) {
        return ;
    }
    else {
        m_point[i] = pt;
    }
}



void        Triangle3D::setPoints(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3)
{
    m_point[0] = pt1;
    m_point[1] = pt2;
    m_point[2] = pt3;
}




Triangle3D& Triangle3D::operator =(const Triangle3D& triangle)
{
    if ( this == &triangle ) {
        return *this;
    }

    m_point[0] = triangle.m_point[0];
    m_point[1] = triangle.m_point[1];
    m_point[2] = triangle.m_point[2];

    return *this;
}



Plane Triangle3D::getPlane() const
{
    return Plane(m_point[0], m_point[1], m_point[2]);
}

