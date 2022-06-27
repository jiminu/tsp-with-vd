#include <float.h>
#include "TriangleOnVDFace.h"
using namespace V::GeometryTier;


TriangleOnVDFace::TriangleOnVDFace()
: m_distanceToGenerator(0.0)
{
    for( rg_INT i=0; i<3; i++)  {
        m_distance[i] = 0.0;
    }
}

TriangleOnVDFace::TriangleOnVDFace(const TriangleOnVDFace& aTriangle)
: Triangle(aTriangle)
{
    m_distanceToGenerator = aTriangle.m_distanceToGenerator;
    for( rg_INT i=0; i<3; i++)  {
        m_distance[i] = aTriangle.m_distance[i];
    }
}

TriangleOnVDFace::~TriangleOnVDFace()
{}


rg_REAL TriangleOnVDFace::getDistanceToGenerator() const
{
    return m_distanceToGenerator;
}

rg_REAL TriangleOnVDFace::getDistance(const rg_INT& i) const
{
    if ( i<0 || i>2 )
        return -DBL_MAX;
    else
        return m_distance[i];
}

void    TriangleOnVDFace::setDistanceToGenerator(const rg_REAL& distance)
{
    m_distanceToGenerator = distance;
}

void    TriangleOnVDFace::setDistance(const rg_INT& i, const rg_REAL& distance)
{
    if ( i<0 || i>2 )
        return;
    else
        m_distance[i] = distance;
}

TriangleOnVDFace& TriangleOnVDFace::operator =(const TriangleOnVDFace& aTriangle)
{
    if ( this == &aTriangle )
        return *this;

    for( rg_INT i=0; i<3; i++)  {
        m_vertex[i]         = aTriangle.m_vertex[i];
        m_normalOnVertex[i] = aTriangle.m_normalOnVertex[i];
        m_distance[i]       = aTriangle.m_distance[i];
    }

    m_distanceToGenerator = aTriangle.m_distanceToGenerator;

    return *this;
}

