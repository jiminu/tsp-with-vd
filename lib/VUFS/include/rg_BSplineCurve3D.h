//********************************************************************
//
//	  FILENAME    : rg_BSplineCurve3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_BSplineCurve3D 
//           which define B-Spline rg_Curve and its property. 
//                          
//	  CLASS NAME  : rg_BSplineCurve3D
//
//    BASE CLASS  : rg_Curve      
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun 1996    
//
//    HISTORY     :
//          BY Young-Song Cho.  13 Jul. 1997
//              rg_INT rg_BSplineCurve3D::getNumOfKnotSpan() const
//          BY Young-Song Cho.  22 Apr. 1999
//	            rg_FLAG isPlanarCurve(rg_Point3D& normal) const;
//          By Taeboom Jang    1999. 10.14 
//              rg_Point3D  getStartPoint() const;
//              rg_Point3D  getEndPoint() const;
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_BSPLINECURVE3D_H
#define _RG_BSPLINECURVE3D_H

#include "rg_Curve.h"

#include "rg_Const.h"
//#include "DefConst.h"
#include "rg_Point3D.h"
#include "rg_dList.h"
#include "rg_BoundingBox3D.h"
#include "rg_BoundingBox2D.h"
#include "rg_Polyline3D.h"
#include "rg_Polyline2D.h"
#include "rg_TMatrix3D.h"

class rg_BSplineCurve3D : public rg_Curve
{
protected:
    rg_INT    numOfCtrlPts;
    rg_Point3D* ctrlPts;
    rg_ORDER  order;

public:
////    Constructor & Destructor
    rg_BSplineCurve3D();
    rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                    const rg_Planarity    &newPlanarity );
    rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                    const rg_Planarity    &newPlanarity,
                    const rg_ORDER        &newOrder );
    rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                    const rg_Planarity    &newPlanarity,
                    const rg_INT          &num );
    rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                    const rg_Planarity    &newPlanarity,
                    const rg_INT          &num, 
                    rg_Point3D*             newControlP );
    rg_BSplineCurve3D( const rg_ORDER &newOrder );
    rg_BSplineCurve3D( const rg_INT &num );
    rg_BSplineCurve3D( const rg_INT &num, 
                    rg_Point3D*    newControlP );
    rg_BSplineCurve3D( const rg_ORDER &newOrder,
                    const rg_INT   &num, 
                    rg_Point3D*       newControlP );
    ////  Constructor      : March 13 1997
    rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                    const rg_Planarity    &newPlanarity,
                    const rg_ORDER        &newOrder,
                    const rg_INT          &num, 
                    rg_Point3D*             newControlP );

    ////  Copy Constructor : March 13 1997
    rg_BSplineCurve3D( const rg_BSplineCurve3D &curve);

    virtual ~rg_BSplineCurve3D();

////    Get Functions.
    rg_ORDER        getOrder() const;
    rg_INT          getNumOfCtrlPts() const;
    rg_Point3D      getCtrlPt( const rg_INT &i ) const;
    rg_Point3D*     getCtrlPts() const;
	rg_REAL         getKnotValue( const rg_INT &i ) const;
    rg_INT          getNumOfKnotSpan() const;
    rg_Point3D*     accessCtrlPts() const;

    rg_Point3D  getStartPoint() const;
    rg_Point3D  getEndPoint() const;

////    Set Functions.
    void         setOrder( const rg_ORDER &newOrder );
    void         setNumOfCtrlPts( const rg_INT &numOfCtrlPt );
    void         setCtrlPt( const rg_INT   &i, 
                                     const rg_Point3D &newPt );
    void         setCtrlPts( const rg_INT         &numOfCtrlPt, 
                                   rg_Point3D*     newCtrlPlygn );

////    Calculations & Operations
    virtual rg_REAL   evaluateBasisFunc( const rg_INT  &index, 
                                      const rg_REAL &param,
                                      const rg_INT  &Order ) const;
    virtual rg_Point3D  evaluatePt( const rg_REAL &u ) const;
    virtual rg_Point3D* evaluatePtsInEvenParameter( const rg_INT &numfPtOnCurve ) const;
    rg_Polyline2D makePolyline2DInEvenParameter(const rg_INT &numOfPts) const;

    virtual rg_Polyline3D makePolyline3DInEvenParameter(const rg_INT &numOfPts) const;
	bool    isPlanarCurve(rg_Point3D& normal) const;
    void    transform(const rg_TMatrix3D& transform);
    rg_Point3D  evaluateCenterOfCtrlPts() const;
    void        reverseTrace();

    rg_REAL     evaluateCurveLengthByCtrlPts() const;

////    Derivative
    rg_BSplineCurve3D makeDerivative();

    rg_REAL getCurvatureInXY( const rg_REAL &u );
    rg_REAL getCurvatureInXZ( const rg_REAL &u );
    rg_REAL getCurvatureInYZ( const rg_REAL &u );

	rg_INT   findIndexOfNearestCtrlPt(const rg_Point3D& pt) const;
	rg_BoundingBox3D makeBoundingBox3D() const;
    rg_BoundingBox2D makeBoundingBox2D() const;

////    Operator Overloading
    rg_BSplineCurve3D& operator =(const rg_BSplineCurve3D &curve);

////    FAQ functions
    rg_FLAG isNull() const;

    void   removeAll();
};

#endif


