#ifndef _RG_SPHERE3D_H
#define _RG_SPHERE3D_H

#include "rg_Const.h"
#include "rg_Point3D.h"

class rg_Sphere3D
{
protected:
    rg_Point3D center;
    rg_REAL    radius;
public:
//    rg_Sphere3D();
    rg_Sphere3D(const rg_Point3D&  tCenter=rg_Point3D(0.0,0.0,0.0),
                const rg_REAL&     tRadius=0.0);
    rg_Sphere3D(const rg_Sphere3D&   temp);
//  get functions
    rg_Point3D getCenter() const;
    rg_REAL    getRadius() const;

//  set functions
    void setCenter(const rg_Point3D& tCenter);
    void setRadius(const rg_REAL& tRadius);

    void setAll( const rg_Sphere3D& temp);

//  operator
    rg_Sphere3D& operator=(const rg_Sphere3D& temp);
};

#endif



