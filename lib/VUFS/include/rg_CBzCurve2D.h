/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_CBzCurve2D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_CBzCurve2D 
//           which define cubic Bezier rg_Curve and its property. 
//                          
//	  CLASS NAME  : rg_CBzCurve2D
//
//    BASE CLASS  : rg_BzCurve2D
//      
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jun 1997    
//
//           Copyright ⓒ 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _RG_CBZCURVE2D_H
#define _RG_CBZCURVE2D_H

#include "rg_Const.h"
#include "rg_Point2D.h"
#include "rg_BzCurve2D.h"
#include "rg_QBzCurve2D.h"
#include "rg_dList.h"
//#include "rg_RQBzCurve2D.h"

class rg_RQBzCurve2D; // 왜 이것을 반드시 선언해 주어야 할까? 위의 include만 시켜 주면 안돼나?

class rg_CBzCurve2D : public rg_BzCurve2D 
{
public:
	rg_CBzCurve2D();
	rg_CBzCurve2D(const rg_Point2D *ctrlpt);
	rg_CBzCurve2D(const rg_Point2D &p0, const rg_Point2D &p1, 
					   const rg_Point2D &p2, const rg_Point2D &p3);
	rg_CBzCurve2D(const rg_CBzCurve2D &curve);
	~rg_CBzCurve2D();

	//Operations
	rg_INT        typeOfCurve() const;
	rg_dList<rg_REAL> approximateRQBzCurves(const rg_REAL &t0, 
					 				         const rg_REAL &t1,
											 rg_dList<rg_RQBzCurve2D> &rqBzCurveList) const;
	rg_RQBzCurve2D makeOneRQBzCurve(const rg_REAL &t0,
							   const rg_REAL &t1) const;
	rg_dList<rg_CBzCurve2D> subdivideTwoCurve(const rg_REAL &t);
	void subdivideTwoCurve(const rg_REAL &t, 
						   rg_CBzCurve2D &curve1,
						   rg_CBzCurve2D &curve2);
	rg_INT     isNotSelfIntersection() const;
    rg_REAL* getInflectionParameter() const;

	//Access elements
	void   setCtrlPts(const rg_Point2D &b0, const rg_Point2D &b1, 
						const rg_Point2D &b2, const rg_Point2D &b3);
	void   setCtrlPts(const rg_Point2D *point);
};

#endif  // rg_CBzCurve2D class


