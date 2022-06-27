#ifndef _RIGHTCIRCULARCONE_H
#define _RIGHTCIRCULARCONE_H

#include "rg_Const.h"
#include "rg_Point3D.h"
#include "Circle3D.h"
#include "Sphere.h"

namespace V {

namespace GeometryTier {


class RightCircularCone
{
private:
	rg_Point3D m_apex;
	Circle3D   m_baseCircle;

public:
    RightCircularCone();
    RightCircularCone(const rg_Point3D& apex, const Circle3D& baseCircle);
    RightCircularCone(const RightCircularCone& cone);
    ~RightCircularCone();

    rg_Point3D getApex() const;
    Circle3D   getBaseCircle() const;
    rg_REAL    getHeight() const;
    rg_REAL    getSlantHeight() const;
    rg_REAL    getHalfVertexAngle() const;

    void       setApex(const rg_Point3D& apex);
    void       setBaseCircle(const Circle3D& baseCircle);
    void       setCone(const rg_Point3D& apex, const Circle3D& baseCircle);

    RightCircularCone& operator =(const RightCircularCone& cone);

    void       makeFromSphericalTriangle( const Sphere&     sphere, 
                                          const rg_Point3D& pt1OnSphere, 
                                          const rg_Point3D& pt2OnSphere, 
                                          const rg_Point3D& pt3OnSphere);

    rg_BOOL    isOnBaseCircle(const rg_Point3D& point) const;
    rg_BOOL    isContainedInCone(const rg_Point3D& point) const;
    rg_BOOL    isContainedInInfiniteCone(const rg_Point3D& point) const;
    rg_BOOL    isInsideOfInfiniteCone(const rg_Point3D& point) const;

};


} // namespace GeometryTier

} // namespace V


#endif

