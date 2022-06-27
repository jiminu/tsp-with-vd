////////////////////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : DistFunc.h
//	  
//    DESCRIPTION : 
//           This is the implementation of the class DistFunc
//           and is used for representing the distance function defined between geometric entity.
//
//    AUTHOR      : Ryu, Jooghyun
//    START DATE  : Jan. 18, 2010
//
//           Copyright ¨Ï 2009 by Voronoi Diagram Research Center, Hanyang University
//
/////////////////////////////////////////////////////////////////////////////////////

#ifndef _DISTFUNC_H_
#define _DISTFUNC_H_

#include "rg_RelativeOp.h"
#include "rg_dList.h"
#include "float.h"

enum InvolvedEntityType {UNDEFINED, LINESEGMENT_POINT, LINESEGMENT_CIRCLE, LINESEGMENT_ARC};

enum LocalPtStatus {LOCAL_UNKNOWN, LOCAL_MIN, LOCAL_MAX, LOCAL_EXTREMUM, LOCAL_NONE};

class DistFunc
{
protected:
	InvolvedEntityType m_type;
	rg_REAL        m_interval[ 2 ];
	rg_REAL*       m_sortedLocalExtrema;
	LocalPtStatus* m_statusOfSortedLocalExtrema;
	rg_INT         m_numOfLocalExtrema;
	
	rg_FLAG        isThisPtIn(const rg_REAL& point, rg_dList<rg_REAL>& ptList) const;
	
public:
	DistFunc();
	DistFunc(const DistFunc& distFunc);
	DistFunc(const InvolvedEntityType& type, 
		     const rg_REAL& lower, 
			 const rg_REAL& upper);
	virtual ~DistFunc();

	InvolvedEntityType getType() const;
	void getInterval(rg_REAL interval[]);
// 	virtual rg_REAL evaluate(const rg_REAL& point) const;
// 	virtual rg_REAL evaluateDerivative(const rg_REAL& point) const;
// 	virtual rg_REAL evaluateDerivativeOfDerivative(const rg_REAL& point) const;
// 	virtual rg_INT  computeLocalMinimaWithoutConsideringBoundary(rg_dList<rg_REAL>& localMinima);
// 	virtual rg_INT  computeLocalMinima(rg_dList<rg_REAL>& localMinima);
// 
// 	virtual LocalPtStatus  getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point) const;
// 	virtual LocalPtStatus  getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point, rg_REAL& fVal) const;
// 	virtual LocalPtStatus  getStatusOfNearestExtremum(const rg_REAL& point) const;
	
	void setType(const InvolvedEntityType& type);
	void setInterval(const rg_REAL& lower, const rg_REAL& upper);

	DistFunc& operator =(const DistFunc& distFunc);
};

#endif

