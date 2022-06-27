#include <float.h>
#include "Triangle.h"
using namespace V::GeometryTier;


Triangle::Triangle()
{
    for( rg_INT i=0; i<3; i++)  {
        m_vertex[i].setPoint(-DBL_MAX, -DBL_MAX, -DBL_MAX);
        m_normalOnVertex[i].setPoint(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    }
}

Triangle::Triangle(const Triangle& aTriangle)
{
    for( rg_INT i=0; i<3; i++)  {
        m_vertex[i]         = aTriangle.m_vertex[i];
        m_normalOnVertex[i] = aTriangle.m_normalOnVertex[i];
    }
}

Triangle::~Triangle()
{}


rg_Point3D Triangle::getVertex(const rg_INT& i) const
{
    if ( i<0 || i>2 )
        return rg_Point3D(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    else
        return m_vertex[i];
}

rg_Point3D Triangle::getNormalOnVertex(const rg_INT& i) const
{
    if ( i<0 || i>2 )
        return rg_Point3D(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    else
        return m_normalOnVertex[i];
}


void Triangle::setVertex(const rg_INT& i, const rg_Point3D& aVertex)
{
    if ( i<0 || i>2 )
        return;
    else
        m_vertex[i] = aVertex;

}

void Triangle::setNormalOnVertex(const rg_INT& i, const rg_Point3D& aNormal)
{
    if ( i<0 || i>2 )
        return;
    else
        m_normalOnVertex[i] = aNormal;
}


Triangle& Triangle::operator =(const Triangle& aTriangle)
{
    if ( this == &aTriangle )
        return *this;

    for( rg_INT i=0; i<3; i++)  {
        m_vertex[i]         = aTriangle.m_vertex[i];
        m_normalOnVertex[i] = aTriangle.m_normalOnVertex[i];
    }

    return *this;
}

