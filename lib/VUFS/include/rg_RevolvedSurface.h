#ifndef _RG_REVOLVEDSURFACE3D_H
#define _RG_REVOLVEDSURFACE3D_H

#include "rg_Line3D.h"
#include "rg_Outline3D.h"
#include "rg_NURBSplineSurface3D.h"

class rg_RevolvedSurface: public rg_NURBSplineSurface3D 
{
private:
    rg_Outline3D generatrix;

    rg_Line3D   axis;
    rg_REAL     angle;

public:
    rg_RevolvedSurface();
    rg_RevolvedSurface(const rg_Outline3D& tGeneratix);
    rg_RevolvedSurface(const rg_Outline3D& tGeneratix,
                    const rg_Line3D&   tAxis,
                    const rg_REAL&     tAngle);
    rg_RevolvedSurface(const rg_RevolvedSurface& temp);
    ~rg_RevolvedSurface();
    // get function
    rg_Outline3D getGeneratrix() const;
    rg_Line3D   getAxis() const;
    rg_REAL     getAngle() const;

    // set function
    void   setGeneratrix(const rg_Outline3D& tGeneratix );
    void   setAxis(const rg_Line3D& tAxis);
    void   setAngle(const rg_REAL& tAngle);

    void   makeSurface();
    
    // operator overloading
    rg_RevolvedSurface& operator=(const rg_RevolvedSurface& temp);
};

#endif


