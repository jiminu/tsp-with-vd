////////////////////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : SquaredDistFuncLineSegCircle.h
//	  
//    DESCRIPTION : 
//           This is the implementation of the class DistFuncLineSegCircle
//           and is used for representing the squared distance function defined between a line segment and a circle.
//
//    AUTHOR      : Ryu, Jooghyun
//    START DATE  : Dec. 28, 2009
//
//           Copyright ¨Ï 2009 by Voronoi Diagram Research Center, Hanyang University
//
/////////////////////////////////////////////////////////////////////////////////////


#ifndef _SQUARED_DISTFUNC_LINESEG_CIRCLE_H_
#define _SQUARED_DISTFUNC_LINESEG_CIRCLE_H_

#include "DistFunc.h"
#include "Circle3D.h"
#include "rg_QuarticPolynomial.h"
#include "rg_CubicPolynomial.h"
#include "rg_QuadraticPolynomial.h"

const rg_INT  NUM_COEFF_SDIST_FUNC_BTW_LINESEG_CIRCLE =  7;

class SquaredDistFuncLineSegCircle : public DistFunc
{
protected:
	rg_REAL        m_coeff[NUM_COEFF_SDIST_FUNC_BTW_LINESEG_CIRCLE];
	rg_REAL        m_radius;

public:
	SquaredDistFuncLineSegCircle();
	SquaredDistFuncLineSegCircle(const SquaredDistFuncLineSegCircle& sDistFunc);
	SquaredDistFuncLineSegCircle(const LineSegment3D& lineSegment, 
		                         const Circle3D& circle,
								 rg_REAL lower = 0.0,
								 rg_REAL upper = 1.0);
	~SquaredDistFuncLineSegCircle();

	rg_INT computeLocalMinimaInsideInterval(rg_dList<rg_REAL>& localMinima);
	rg_INT computeLocalExremaInsideInterval(rg_dList<rg_REAL>& localExtrema, rg_dList<LocalPtStatus>& localStatuses) const;

	LocalPtStatus getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point) const;
	//LocalPtStatus getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point, rg_REAL& fVal) const;
	LocalPtStatus getStatusOfNearestExtremum(const rg_REAL& point) const;

	rg_INT computeLocalMinima(rg_dList<rg_REAL>& localMinima);
	rg_INT computeLocalExrema(rg_dList<rg_REAL>& localExtrema, rg_dList<LocalPtStatus>& localStatuses) const;
	void   computeCriticalPoints(rg_dList<rg_REAL>& roots) const;
	void   sortLocalExtremaAccordingToParam(rg_dList<rg_REAL>& unsortedLocalExtrema, rg_dList<LocalPtStatus>& unsortedLocalPtStatuses);
	void computeSqdDerivativeForSqdDistFunc(rg_QuarticPolynomial& derivative) const;
	void computeDerivativeForSqdDistFunc(rg_Polynomial& derivative, rg_DEGREE& degree) const;
	LocalPtStatus isThisPtLocalMinOrMax(const rg_REAL& point, const LocalPtStatus& status) const;
	LocalPtStatus isThisPtLocalMinOrMax(const rg_REAL& point) const;	
	rg_REAL evaluate(const rg_REAL& point) const;
	//rg_REAL evaluateDerivative(const rg_REAL& point) const;
	//rg_REAL evaluateDerivativeOfDerivative(const rg_REAL& point) const;

	rg_INT computeInflectionPoints(rg_dList<rg_REAL>& inflectionPts) const;
	rg_INT computeRootOfThirdDerivative(rg_dList<rg_REAL>& root) const;

	SquaredDistFuncLineSegCircle& operator=(const SquaredDistFuncLineSegCircle& sDistFunc);
};

#endif

