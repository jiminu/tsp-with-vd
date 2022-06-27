#ifndef _RG_POLYLINE2D_H
#define _RG_POLYLINE2D_H

#include "rg_Const.h"
#include "rg_Curve.h"
#include "rg_Point2D.h"

class rg_Polyline2D: public rg_Curve
{
protected:
    rg_INT      numOfPts;
    rg_Point2D*  pts;

public:
    rg_Polyline2D();
    rg_Polyline2D(const rg_INT            &tNumOfPts,
                  const rg_Point2D* const  tPts);
    rg_Polyline2D(const rg_Polyline2D&     temp);

    ~rg_Polyline2D();
    void removeAll();

    rg_INT        getNumOfPts() const;
    rg_Point2D    getPt(const rg_INT &index)  const;
    rg_Point2D*   getPts() const;
    rg_Point2D*   accessPts() const;


    void          setNumOfPts( const rg_INT& tNumOfPts);
    void          setPts(const rg_Point2D* const tPts);
    void          setPt(const rg_INT     &index,
                        const rg_Point2D &tPt);

    void          setAll(const rg_Polyline2D& temp);
    
    rg_Polyline2D& operator=(const rg_Polyline2D& temp);

    rg_REAL        evaluateDistance(const rg_Point2D& pt) const;    
    // FAQ
    rg_FLAG        isNull() const;
    rg_FLAG        isOverlapped(const rg_Point2D& pt, const rg_REAL& tolerance ) const;

};
#endif


