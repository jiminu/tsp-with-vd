#include "RightCircularCone.h"
using namespace V::GeometryTier;


RightCircularCone::RightCircularCone()
{
}



RightCircularCone::RightCircularCone(const rg_Point3D& apex, const Circle3D& baseCircle)
: m_apex(apex), m_baseCircle(baseCircle)
{
    Plane   basePlane = m_baseCircle.getPlaneContainingThisCircle();
    rg_REAL dist      = basePlane.distanceFromPoint(m_apex);

    if ( rg_NEG(dist) ) {
        m_baseCircle.reverse();
    }
}



RightCircularCone::RightCircularCone(const RightCircularCone& cone)
: m_apex(cone.m_apex), m_baseCircle(cone.m_baseCircle)
{
}



RightCircularCone::~RightCircularCone()
{
}




rg_Point3D RightCircularCone::getApex() const
{
    return m_apex;
}



Circle3D   RightCircularCone::getBaseCircle() const
{
    return m_baseCircle;
}



rg_REAL    RightCircularCone::getHeight() const
{
    Plane   basePlane = m_baseCircle.getPlaneContainingThisCircle();
    rg_REAL height    = basePlane.distanceFromPoint(m_apex);

    return height;
}



rg_REAL    RightCircularCone::getSlantHeight() const
{
    Plane   basePlane = m_baseCircle.getPlaneContainingThisCircle();
    rg_REAL height    = basePlane.distanceFromPoint(m_apex);
    rg_REAL radius    = m_baseCircle.getRadius();

    rg_REAL slantHeight = sqrt( radius*radius + height*height );

    return slantHeight;
}



rg_REAL    RightCircularCone::getHalfVertexAngle() const
{
    Plane   basePlane   = m_baseCircle.getPlaneContainingThisCircle();
    rg_REAL radius      = m_baseCircle.getRadius();

    rg_REAL height      = basePlane.distanceFromPoint(m_apex);
    rg_REAL slantHeight = sqrt( radius*radius + height*height );

    rg_REAL halfVertexAngle = acos( height/slantHeight );

    return halfVertexAngle;
}



void       RightCircularCone::setApex(const rg_Point3D& apex)
{
    m_apex       = apex;
}



void       RightCircularCone::setBaseCircle(const Circle3D& baseCircle)
{
    m_baseCircle = baseCircle;

    Plane   basePlane = m_baseCircle.getPlaneContainingThisCircle();
    rg_REAL dist      = basePlane.distanceFromPoint(m_apex);

    if ( rg_NEG(dist) ) {
        m_baseCircle.reverse();
    }
}



void       RightCircularCone::setCone(const rg_Point3D& apex, const Circle3D& baseCircle)
{
    m_apex       = apex;
    m_baseCircle = baseCircle;

    Plane   basePlane = m_baseCircle.getPlaneContainingThisCircle();
    rg_REAL dist      = basePlane.distanceFromPoint(m_apex);

    if ( rg_NEG(dist) ) {
        m_baseCircle.reverse();
    }
}




RightCircularCone& RightCircularCone::operator =(const RightCircularCone& cone)
{
    if ( this == &cone ) {
        return *this;
    }

    m_apex       = cone.m_apex;
    m_baseCircle = cone.m_baseCircle;

    return *this;
}



void       RightCircularCone::makeFromSphericalTriangle( const Sphere&     sphere, 
                                                         const rg_Point3D& pt1OnSphere, 
                                                         const rg_Point3D& pt2OnSphere, 
                                                         const rg_Point3D& pt3OnSphere )
{
    m_apex       = sphere.getCenter();

    Plane    basePlane( pt1OnSphere, pt2OnSphere, pt3OnSphere );
    
    Circle3D intersectionCircle;
    sphere.intersect( basePlane, intersectionCircle );


    m_baseCircle = intersectionCircle;

    rg_REAL dist = basePlane.distanceFromPoint(m_apex);

    if ( rg_NEG(dist) ) {
        m_baseCircle.reverse();
    }
}



rg_BOOL    RightCircularCone::isOnBaseCircle(const rg_Point3D& point) const
{
    return m_baseCircle.isOnCircle(point);
}



rg_BOOL    RightCircularCone::isContainedInCone(const rg_Point3D& point) const
{
    if ( isContainedInInfiniteCone(point) ) {
        Plane   basePlane = m_baseCircle.getPlaneContainingThisCircle();
        rg_REAL dist      = basePlane.distanceFromPoint(point);

        if ( !rg_NEG( dist ) ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_BOOL    RightCircularCone::isContainedInInfiniteCone(const rg_Point3D& point) const
{
    rg_REAL halfVertexAngle = getHalfVertexAngle();

    rg_Point3D vecToCenter = -m_baseCircle.getNormal();
    rg_Point3D vecToPoint  = point - m_apex;

    rg_REAL angle = vecToCenter.angle(vecToPoint);

    if ( rg_GE(halfVertexAngle, angle) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}







rg_BOOL    RightCircularCone::isInsideOfInfiniteCone(const rg_Point3D& point) const
{
    rg_REAL halfVertexAngle = getHalfVertexAngle();

    rg_Point3D vecToCenter = -m_baseCircle.getNormal();
    rg_Point3D vecToPoint  = point - m_apex;

    rg_REAL angle = vecToCenter.angle(vecToPoint);

    if ( rg_GT(halfVertexAngle, angle) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}


