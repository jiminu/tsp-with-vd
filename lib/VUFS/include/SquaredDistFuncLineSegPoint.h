////////////////////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : SquaredDistFuncLineSegPoint.h
//	  
//    DESCRIPTION : 
//           This is the implementation of the class SquaredDistFuncLineSegPoint
//           and is used for representing the squared distance function defined between a line segment and a point.
//
//    AUTHOR      : Ryu, Jooghyun
//    START DATE  : Jan. 14, 2010
//
//           Copyright ¨Ï 2010 by Voronoi Diagram Research Center, Hanyang University
//
/////////////////////////////////////////////////////////////////////////////////////


#ifndef _SQUARED_DISTFUNC_LINESEG_POINT_H_
#define _SQUARED_DISTFUNC_LINESEG_POINT_H_

#include "DistFunc.h"
#include "rg_QuadraticPolynomial.h"
#include "LineSegment3D.h"
#include "rg_dList.h"

const rg_INT NUM_COEFF_SDIST_FUNC_BTW_LINESEG_POINT = 3;

class SquaredDistFuncLineSegPoint : public DistFunc
{
private:
	rg_QuadraticPolynomial m_distFunc;

public:
	SquaredDistFuncLineSegPoint();
	SquaredDistFuncLineSegPoint(const SquaredDistFuncLineSegPoint& sDistFunc);
	SquaredDistFuncLineSegPoint(const LineSegment3D& linesegment, 
		                        const rg_Point3D& point, 
								rg_REAL lower = 0.0,
								rg_REAL upper = 1.0);
	~SquaredDistFuncLineSegPoint();

	rg_INT   computeLocalMinimaInsideInterval(rg_dList<rg_REAL>& localMinima);
	rg_INT   computeLocalExremaWithoutConsideringBoundary();
	rg_INT   computeLocalMinima(rg_dList<rg_REAL>& localMinima);
	rg_INT   computeLocalExrema();

	LocalPtStatus getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point) const;
	//LocalPtStatus getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point, rg_REAL& fVal) const;
	LocalPtStatus getStatusOfNearestExtremum(const rg_REAL& point) const;

	rg_REAL  evaluate(const rg_REAL& point) const;

	//rg_INT computeInflectionPoints(rg_dList<rg_REAL>& inflectionPts) const;

	SquaredDistFuncLineSegPoint& operator =(const SquaredDistFuncLineSegPoint& distFunc);
};

#endif

