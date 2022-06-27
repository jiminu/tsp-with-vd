#ifndef _MESHFORSPHERE_H
#define _MESHFORSPHERE_H

#include "rg_Const.h"
#include "SolidObject.h"
#include "Sphere.h"


namespace V {

namespace GeometryTier {


class MeshForSphere : public SolidObject
{
private:
    rg_INT m_resolution;
    Sphere m_sphere;

public:
    MeshForSphere();
    MeshForSphere(const rg_INT& resolution, const Sphere& sphere);
    //MeshForSphere(const MeshForSphere& sphereMesh);
    ~MeshForSphere();

    Sphere getSphere() const;

    void   setSphere(const Sphere& sphere);
    void   setResolution(const rg_INT& resolution);
    //MeshForSphere& operator =(const MeshForSphere& sphereMesh);
};

} // namespace GeometryTier

} // namespace V


#endif

