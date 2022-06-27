#include "rg_Sphere3D.h"

rg_Sphere3D::rg_Sphere3D(const rg_Point3D&  tCenter,
                         const rg_REAL&     tRadius)
:center(tCenter), radius(tRadius)
{
}

rg_Sphere3D::rg_Sphere3D(const rg_Sphere3D&   temp)
:center(temp.center), radius(temp.radius)
{
}

//  get functions
rg_Point3D rg_Sphere3D::getCenter() const
{
    return center;
}

rg_REAL    rg_Sphere3D::getRadius() const
{
    return radius;
}

//  set functions
void rg_Sphere3D::setCenter(const rg_Point3D& tCenter)
{
    center=tCenter;
}

void rg_Sphere3D::setRadius(const rg_REAL& tRadius)
{
    radius=tRadius;
}

void rg_Sphere3D::setAll( const rg_Sphere3D& temp)
{
    center=temp.center;
    radius=temp.radius;
}


//  operator
rg_Sphere3D& rg_Sphere3D::operator=(const rg_Sphere3D& temp)
{
    if ( this == &temp)
    {
        return *this;
    }
    center=temp.center;
    radius=temp.radius;
    return *this;
}


