#include "rg_Frustum3D.h"

rg_Frustum3D::rg_Frustum3D(const unsigned rg_INT& tID)
: rg_Geometry(FRUSTUM,tID), rg_PrimitiveCylinder3D(FRUSTUM,tID), rg_PrimitiveFrustum3D(FRUSTUM,tID), rg_Cylinder3D(FRUSTUM,tID)
{
}

rg_Frustum3D::rg_Frustum3D(const GeometryType& tType, const unsigned rg_INT& tID)
: rg_Geometry(tType,tID), rg_PrimitiveCylinder3D(tType,tID), rg_PrimitiveFrustum3D(tType,tID), rg_Cylinder3D(tType,tID)
{
}

rg_Frustum3D::rg_Frustum3D(const rg_Frustum3D& temp)
: rg_Geometry(temp), rg_PrimitiveCylinder3D(temp), rg_PrimitiveFrustum3D(temp), rg_Cylinder3D(temp)
{
}

rg_Frustum3D::~rg_Frustum3D()
{
}

void   rg_Frustum3D::setFrustum(const rg_Frustum3D& temp)
{
    rg_Cylinder3D::setCylinder(temp);
    rg_PrimitiveFrustum3D::setTopRadius(temp.getTopRadius());
}

rg_Frustum3D& rg_Frustum3D::operator=(const rg_Frustum3D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }

    rg_Frustum3D::setFrustum(temp);
    return *this;
}


