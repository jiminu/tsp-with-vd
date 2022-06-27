#include "rg_PrimitiveFrustum3D.h"



rg_PrimitiveFrustum3D::rg_PrimitiveFrustum3D(const unsigned rg_INT& tID)
: rg_Geometry( PRIMITIVE_FRUSTUM ), rg_PrimitiveCylinder3D(PRIMITIVE_FRUSTUM,tID), topRadius(0.0)
{
}

rg_PrimitiveFrustum3D::rg_PrimitiveFrustum3D(const GeometryType& tType,const unsigned rg_INT& tID)
: rg_Geometry( tType ), rg_PrimitiveCylinder3D(tType,tID), topRadius(0.0)
{
}

rg_PrimitiveFrustum3D::rg_PrimitiveFrustum3D(const rg_PrimitiveFrustum3D& temp)
: rg_Geometry(temp), rg_PrimitiveCylinder3D( temp ), topRadius(temp.topRadius)
{
}

rg_PrimitiveFrustum3D::~rg_PrimitiveFrustum3D()
{
}

rg_REAL   rg_PrimitiveFrustum3D::getBottomRadius() const
{
    return rg_PrimitiveCylinder3D::getRadius();
}

rg_REAL   rg_PrimitiveFrustum3D::getTopRadius() const
{
    return topRadius;
}

void   rg_PrimitiveFrustum3D::setBottomRadius(const rg_REAL& tBottomRadius)
{
    rg_PrimitiveCylinder3D::setRadius( tBottomRadius);
}
    
void   rg_PrimitiveFrustum3D::setTopRadius(const rg_REAL& tTopRadius)
{
    topRadius=tTopRadius;
}

void   rg_PrimitiveFrustum3D::setPrimitiveFrustum(const rg_PrimitiveFrustum3D& temp)
{
    rg_PrimitiveCylinder3D::setPrimitiveCylinder(temp);
    topRadius=temp.topRadius;
}

rg_PrimitiveFrustum3D& rg_PrimitiveFrustum3D::operator=(const rg_PrimitiveFrustum3D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }

    rg_PrimitiveFrustum3D::setPrimitiveFrustum(temp);

    return *this;
}


