//********************************************************************
//
//	  FILENAME    : BSplineSurface.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_BSplineSurface3D 
//           which define B-Spline rg_Surface and its property. 
//                          
//	  CLASS NAME  : rg_BSplineSurface3D
//
//    BASE CLASS  : rg_Surface
//      
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun 1996    
//
//    History:
//           1. insert new function related wit knots
//                 rg_REAL*  getKnotVectorOfU() const;
//                 rg_REAL*  getKnotVecotrOfV() const;
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_BSPLINESURFACE3D_H
#define _RG_BSPLINESURFACE3D_H

#include "rg_Const.h"
//#include "DefConst.h"
#include "rg_Point3D.h"

#include "rg_Surface.h"
#include "rg_BSplineCurve3D.h"

class rg_BSplineSurface3D : public rg_Surface
{
protected:
    rg_INT	    rowOfControlNet;
    rg_INT	    columnOfControlNet;
    rg_Point3D** control_net;
    rg_ORDER   u_order;
    rg_ORDER   v_order;

public:
////	Constructor & Destructor.
	rg_BSplineSurface3D();
	rg_BSplineSurface3D( const unsigned rg_INT &newID, 
		              const rg_Planarity    &newPlanarity );
	rg_BSplineSurface3D( const unsigned rg_INT &newID, 
		              const rg_Planarity    &newPlanarity, 
					  const rg_INT          &row, 
					  const rg_INT          &col ); 
	rg_BSplineSurface3D( const unsigned rg_INT &newID, 
		              const rg_Planarity    &newPlanarity, 
					  const rg_ORDER        &uOrder, 
					  const rg_ORDER        &vOrder );
	rg_BSplineSurface3D( const rg_INT &row, 
		              const rg_INT &col );
	rg_BSplineSurface3D( const rg_ORDER &uOrder, 
		              const rg_ORDER &vOrder );
	rg_BSplineSurface3D( const rg_INT   &row, 
		              const rg_INT   &col,
					  const rg_ORDER &uOrder, 
					  const rg_ORDER &vOrder );
	rg_BSplineSurface3D( const unsigned rg_INT &newID, 
		              const rg_Planarity    &newPlanarity, 
					  const rg_INT          &row, 
					  const rg_INT          &col, 
					  const rg_ORDER        &uOrder, 
					  const rg_ORDER        &vOrder );
    ////  Constructor      : March 13 1997
  	rg_BSplineSurface3D( const unsigned rg_INT &newID, 
		              const rg_Planarity    &newPlanarity, 
					  const rg_INT          &row, 
					  const rg_INT          &col, 
					  const rg_ORDER        &uOrder, 
					  const rg_ORDER        &vOrder, 
                      rg_Point3D**             ctrlNet);
    ////  Copy Constructor : March 13 1997
    rg_BSplineSurface3D( const rg_BSplineSurface3D &surface );
	
    virtual ~rg_BSplineSurface3D();

////	Get Functions.
	rg_INT   isValidSurface() const;
	rg_INT	  getRowOfControlNet() const;
	rg_INT	  getColumnOfControlNet() const;	
	rg_ORDER getOrderOfU() const;
	rg_ORDER getOrderOfV() const;
	rg_Point3D getPointOnControlNet( const rg_INT &row, 
		                        const rg_INT &col ) const;
    rg_Point3D** getControlNet() const;
	rg_REAL    getKnotValueOfU( const rg_INT &kIndex ) const;
	rg_REAL    getKnotValueOfV( const rg_INT &kIndex ) const;
    rg_REAL*   getKnotVectorOfU() const;
    rg_REAL*   getKnotVecotrOfV() const;

//    rg_BSplineCurve3D getUIsoparametricCurve() const;

////	Set Functions.
	void    setControlNet( const rg_INT &row, 
                           const rg_INT &col );
	void    setControlNet( const rg_INT &row, 
                           const rg_INT &col, 
                           rg_Point3D**   controlNet );
	void	setOrderOfU( const rg_ORDER &uOrder );
	void	setOrderOfV( const rg_ORDER &vOrder );
	void	setOrderOfSurface( const rg_ORDER &uOrder, 
                               const rg_ORDER &vOrder );
	void	setPointOnControlNet( const rg_INT   &row, 
                                  const rg_INT   &col, 
                                  const rg_Point3D &point);

////	Operating & Calculating.
	virtual rg_REAL evaluateBasisFuncU( const rg_INDEX     &index, 
		                             const rg_PARAMETER &u,
                                     const rg_ORDER     &uOrder );
	virtual rg_REAL evaluateBasisFuncV( const rg_INDEX     &index, 
		                             const rg_PARAMETER &v,
                                     const rg_ORDER     &vOrder );
	virtual rg_Point3D evaluatePt( const rg_PARAMETER &u, 
                              const rg_PARAMETER &v );


////    Derivative.
    rg_BSplineSurface3D derivativeSurfaceOfU();
    rg_BSplineSurface3D derivativeSurfaceOfV();
    rg_BSplineSurface3D derivativeSurfaceOfUV();

////    Operator Overloading.
    rg_BSplineSurface3D& operator =(const rg_BSplineSurface3D &surface);

////    File-Out Function.
	//void	     fileOut( const char* fileName );   

};

#endif

