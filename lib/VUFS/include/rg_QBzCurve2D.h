/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_QBzCurve2D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_QBzCurve2D 
//           which define quadratic Bezier rg_Curve and its property. 
//                          
//	  CLASS NAME  : rg_QBzCurve2D
//
//    BASE CLASS  : rg_BzCurve2D
//      
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jun 1997    
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _RG_QBZCURVE2D_H
#define _RG_QBZCURVE2D_H

#include "rg_Point2D.h"
#include "rg_Point3D.h"
#include "rg_BzCurve2D.h"
#include "rg_dList.h"

class rg_QBzCurve2D : public rg_BzCurve2D
{
public:
	rg_QBzCurve2D();
	rg_QBzCurve2D(const rg_Point2D *ctrlpt);
	rg_QBzCurve2D(const rg_Point2D &p1, const rg_Point2D &p2, const rg_Point2D &p3);
	rg_QBzCurve2D(const rg_Point3D &p1, const rg_Point3D &p2, const rg_Point3D &p3);
	rg_QBzCurve2D(const rg_QBzCurve2D &curve);
	~rg_QBzCurve2D();

	//Operations
	rg_INT     isOriginExteriorPoint();
	rg_INT     ctrlPtsMonotoneInPolarAngle();
	rg_INT     ctrlPtsIncreasingInPolarAngle();
	rg_INT     ctrlPtsDecreasingInPolarAngle();
	rg_REAL  isOriginOnCurve();
	rg_INT     isNotConjugateTangentVector();
//	rg_REAL *getInflectionParameter();
	rg_REAL  getOriginParameter();
	rg_dList<rg_REAL> getParameterOntheAxis();
	rg_REAL *getConjugateTangentParameters();
	rg_REAL  parameterValueOfPt(const rg_Point2D &point);
	rg_REAL  parameterValueOfPt(const rg_REAL &x, const rg_REAL &y);

	//Access elements
	void     setCtrlPts(const rg_Point2D &p1, const rg_Point2D &p2, const rg_Point2D &p3);
	void     setCtrlPts(const rg_Point2D *point);
};

#endif  // rg_QBzCurve2D class


