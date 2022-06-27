#ifndef _RG_POLYLINE3D_H
#define _RG_POLYLINE3D_H

#include "rg_Const.h"
#include "rg_Curve.h"
#include "rg_Point3D.h"

class rg_Polyline3D: public rg_Curve
{
protected:
    rg_INT      numOfPts;
    rg_Point3D*  pts;

public:
    rg_Polyline3D();
    rg_Polyline3D(const rg_INT            &tNumOfPts,
                  const rg_Point3D* const  tPts);
    rg_Polyline3D(const rg_Polyline3D&     temp);

    ~rg_Polyline3D();
    void removeAll();

    rg_INT        getNumOfPts() const;
    rg_Point3D*   getPts()  const;
    rg_Point3D*   accessPts() const;

    void          setNumOfPts( const rg_INT& tNumOfPts);
    void          setPts(const rg_Point3D* const tPts);
    void          setPt(const rg_INT     &index,
                        const rg_Point3D &tPt);

    void          setAll(const rg_Polyline3D& temp);
    
    rg_Polyline3D& operator=(const rg_Polyline3D& temp);
    


};
#endif


