#ifndef _RG_OUTLINE3D_H
#define _RG_OUTLINE3D_H

#include "rg_Point3D.h"
#include "rg_NURBSplineCurve3D.h"

class rg_Outline3D: public rg_NURBSplineCurve3D
{
protected:
    rg_Point3D* passingPts;
    rg_INT    numOfPassingPts;

public:
    rg_Outline3D();
    rg_Outline3D(const rg_Outline3D& temp);
    ~rg_Outline3D();

    // get function
    rg_Point3D* getPassingPts() const;
    rg_Point3D  getPassingPt(const rg_INT& index) const;
    rg_INT    getNumOfPassingPts() const;

    // set function
    void   setPassingPts( rg_Point3D* tPassingPts,
                          const rg_INT&  tNumOfPassingPts);
    void   setPassingPt( const rg_INT& index,
                         const rg_Point3D& pt );

    // make curve
    void fitCurve( const rg_INT& tOrder,
                   const rg_Point3D& startTangent=rg_Point3D(0.0,0.0,0.0),
                   const rg_Point3D& endTangent=rg_Point3D(0.0,0.0,0.0)   );

    // operator overloading
    rg_Outline3D& operator=(const rg_Outline3D& temp);

};

#endif 


