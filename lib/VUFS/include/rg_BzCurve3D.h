//********************************************************************
//
//	  FILENAME    : rg_BzCurve3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class Bezier rg_Curve in 3-D
//           which define a Bezier rg_Curve and its property. 
//                          
//	  CLASS NAME  : rg_BzCurve3D
//
//    BASE CLASS  : rg_Curve
//      
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//
//    HISTORY     : 	
//	      1.  By Dong-Gyou Lee 18 Mar. 1998
//                  void powerToBezierCurve( const rg_DEGREE& dgr,
//                                           const rg_REAL[] paramValues,
//                                           const rg_Matrix& powerCoeff )
//
//	      2.  By Taeboom Jang 1999. 8.27
//                  rg_dList<rg_Point3D> intersectWithPlaneForCubic(const rg_Point3D& plane) const;
//    START DATE  : 9 Jul. 1997    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_BZCURVE3D_H
#define _RG_BZCURVE3D_H

#include "rg_Const.h"
#include "rg_Point3D.h"
#include "rg_Plane3D.h"
#include "rg_ComplexNumber.h"
#include "rg_Matrix.h"
#include "rg_Curve.h"
#include "rg_ListByPtr.h"
#include "rg_Polynomial.h"
#include "rg_dList.h"
#include "rg_BzCurve2D.h"

class rg_BzCurve3D : public rg_Curve
{
protected:
    rg_DEGREE degree;   
    rg_Point3D* ctrlPts;

public:
    ////  Constructor & Destructor
    rg_BzCurve3D();
    rg_BzCurve3D(const rg_DEGREE &dgr);
    rg_BzCurve3D(const rg_DEGREE &dgr, const rg_Point3D* ctrlPts);
    rg_BzCurve3D(const rg_BzCurve3D &curve);
    virtual ~rg_BzCurve3D();

    ////  Access elements
    rg_DEGREE  getDegree() const;
	rg_DEGREE  getOrder() const;
    rg_Point3D   getCtrlPt(const rg_INDEX &i) const;
	rg_BzCurve2D evaluateBzCurve2D() const;

    //return pointer to the copy of control points.
    rg_Point3D*  getCtrlPts() const;

    void    setDegree(const rg_DEGREE &dgr);
	void    setOrder(const rg_DEGREE &dgr);
    void    setCtrlPt(const rg_INDEX &i, const rg_Point3D &pt);
    void    setCtrlPts(const rg_INT &numOfCtrlPt, const rg_Point3D* ctrlPts);
    void    setCurve(const rg_BzCurve3D &curve);

    //  Operations
    virtual rg_Point3D  evaluatePt(const rg_PARAMETER &u) const;
    rg_BzCurve3D      makeDerivative();
    void            raiseDegree(const rg_DEGREE& raisingTimes);

    rg_Point3D**     deCasteljau(const rg_PARAMETER &u);
    rg_REAL        bernstein(const rg_INDEX &i, const rg_DEGREE &dgr, const rg_PARAMETER &t) const;
	rg_Polynomial bernstein(const rg_DEGREE& n, const rg_INDEX& i) const;
    rg_Polynomial* convertBzCurve2Polynomial() const;
    rg_REAL        factorial(const rg_REAL &n) const;

    //  Intersection 
    rg_sListByPtr* intersectOfCubicBezierAndPlane(const rg_Plane3D &plane);
    rg_dList<rg_Point3D> intersectWithPlaneForCubic(const rg_Plane3D& plane) const;
	//	Conversion between power basis & bezier form.
	//	added by Dong-Gyou Lee 18 Mar. 1998
   	void powerToBezierCurve( const rg_DEGREE& dgr,
			                 const rg_REAL paramValues[],
							 const rg_Matrix& powerCoeff );


};

#endif


