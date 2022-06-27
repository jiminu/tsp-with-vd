#ifndef _TRIANGLE3D_H
#define _TRIANGLE3D_H

#include "rg_Point3D.h"
#include "Plane.h"


namespace V {

namespace GeometryTier {


class Triangle3D
{
private:
	rg_Point3D m_point[3];

public:
    Triangle3D();
    Triangle3D(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3);
    Triangle3D(const Triangle3D& triangle);
    ~Triangle3D();

    rg_Point3D  getPoint(const rg_INT& i) const;
    rg_Point3D* getPoints();

    void        setPoint(const rg_INT& i, const rg_Point3D& pt);
    void        setPoints(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3);

    Triangle3D& operator =(const Triangle3D& triangle);

    Plane       getPlane() const;
};

} // namespace GeometryTier

} // namespace V


#endif
