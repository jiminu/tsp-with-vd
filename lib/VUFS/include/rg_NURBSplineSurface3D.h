//********************************************************************
//
//	  FILENAME    : rg_NURBSplineSurface3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_NURBSplineSurface3D 
//           which define Non-uniform Rational B-Spline rg_Surface and its property. 
//                          
//	  CLASS NAME  : rg_NURBSplineSurface3D
//
//    BASE CLASS  : rg_NUBSplineSurface3D
//
//    AUTHOR      : Young-Song Cho
//
//    HISTORY     :
//
//          BY Dong-Gyou Lee 8 Apr. 1998
//                  bool makeTabulatedCylinder( const rg_NURBSplineCurve3D& directionalCurve,
//                                              const rg_Point3D& endPoint )        
//
//    START DATE  : 21 Jun 1996    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_NURBSPLINESURFACE3D_H
#define _RG_NURBSPLINESURFACE3D_H

#include "rg_NUBSplineSurface3D.h"
#include "rg_RBzSurface3D.h"

#include "rg_Const.h"
#include "rg_Point3D.h"

#include "rg_NURBSplineCurve3D.h"

class rg_NURBSplineSurface3D : public rg_NUBSplineSurface3D
{
protected:
      rg_REAL** weight_vectors;  

public:
////    Constructor & Destructor.----------------------------------------------
    rg_NURBSplineSurface3D();
    rg_NURBSplineSurface3D( const rg_INT   &row, 
                         const rg_INT   &col, 
                         const rg_ORDER &uOrder, 
                         const rg_ORDER &vOrder );
    rg_NURBSplineSurface3D( const unsigned rg_INT &newID, 
                         const rg_Planarity    &newPlanarity, 
                         const rg_INT          &row, 
                         const rg_INT          &col, 
                         const rg_ORDER        &uOrder, 
                         const rg_ORDER        &vOrder,
                         rg_Point3D**             ctrlNet,
                         rg_REAL*               uKnotVector,
                         rg_REAL*               vKnotVector,
                         rg_REAL**              weightVector);
    rg_NURBSplineSurface3D( const rg_NURBSplineSurface3D &surface );

    virtual ~rg_NURBSplineSurface3D();

////    Get Functions.---------------------------------------------------------
    rg_REAL   getWeight(const rg_INDEX &i, const rg_INDEX &j) const;
    rg_REAL** getWeightVector() const;

    rg_NURBSplineCurve3D getUIsoparamCurve(const rg_PARAMETER &u);
    rg_NURBSplineCurve3D getVIsoparamCurve(const rg_PARAMETER &v);

	void getPiecewiseSurfaceInPowerForm(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_REAL**** & polyCoeffOfWeight) const;

////    Set Functions.---------------------------------------------------------
    void setControlNetAndWeight(const rg_INDEX& tRow,
                                const rg_INDEX& tCol );
                                        
    void setWeight(const rg_INDEX &i, const rg_INDEX &j, const rg_REAL &weight);
    void setWeightVector(rg_REAL** weightVector);
    void setSurface(const rg_RBzSurface3D& surface);

////    Operating & Calculating.-----------------------------------------------
    virtual rg_Point3D evaluatePt(const rg_PARAMETER &u, const rg_PARAMETER &v);

////    Derivative.------------------------------------------------------------
    rg_Point3D derivativeSurfaceOfU(const rg_PARAMETER &u, const rg_PARAMETER &v);
    rg_Point3D derivativeSurfaceOfV(const rg_PARAMETER &u, const rg_PARAMETER &v);
//    rg_Point3D derivativeSurfaceOfUV(const rg_REAL &u, const rg_REAL &v);

////    rg_Surface Construction Technique.----------------------------------------
    void skinnedSurface(const rg_INT                     &n,
                        const rg_NURBSplineCurve3D* const sectionCurves,
                        rg_REAL*                          parameterization = rg_NULL);

    rg_NURBSplineCurve3D* vDirectionalCurveInterpolation4Skinning(
                        const rg_INT                     &n,
                        const rg_NURBSplineCurve3D* const sectionCurves,
						rg_REAL*                          parameterization = rg_NULL);

////    Operator Overloading.
    rg_NURBSplineSurface3D& operator =(const rg_NURBSplineSurface3D &surface);

//	representation of special surface using NURBS surface form.	             
//	BY Dong-Gyou Lee 8. Apr. 1998
	rg_FLAG makeTabulatedCylinder( const rg_NURBSplineCurve3D& directionalCurve,
		                        const rg_Point3D& endPoint );        

};

#endif


