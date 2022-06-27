#ifndef _AXISALIGNEDBOX_H
#define _AXISALIGNEDBOX_H

#include "rg_Point3D.h"
#include "Sphere.h"

class AxisAlignedBox
{
private:
	rg_Point3D m_minPoint;
	rg_Point3D m_maxPoint;

public:
	AxisAlignedBox();
	AxisAlignedBox(const rg_Point3D& minPt, const rg_Point3D& maxPt);
	AxisAlignedBox(const AxisAlignedBox& box);
	~AxisAlignedBox();

    void        clear();

    rg_Point3D  getMinPoint() const;
    rg_Point3D  getMaxPoint() const;

    rg_Point3D  getSize() const;
    rg_Point3D  getCenter() const;

    void        setMinPoint(const rg_Point3D& minPt);
    void        setMaxPoint(const rg_Point3D& maxPt);
    void        setAxisAlignedBox(const rg_Point3D& minPt, const rg_Point3D& maxPt);

    rg_REAL     computeLengthOfDiagonal() const;
    rg_REAL     computeAreaOfDiagonalFace() const;
    rg_REAL     computeVolume() const;

    AxisAlignedBox& operator =(const AxisAlignedBox& box);

    void       enlarge(const double& ratio);

    void       update(const rg_Point3D& point);
    void       update(const Sphere&     sphere);

    bool       doesInclude(const rg_Point3D& pt) const;
    bool       doesInclude(const Sphere& sphere) const;
};

#endif

