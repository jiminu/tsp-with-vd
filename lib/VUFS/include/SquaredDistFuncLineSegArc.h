////////////////////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : SquaredDistFuncLineSegCircle.h
//	  
//    DESCRIPTION : 
//           This is the implementation of the class DistFuncLineSegArc
//           and is used for representing the squared distance function defined between a line segment and an arc.
//
//    AUTHOR      : Ryu, Jooghyun
//    START DATE  : Jan. 17, 2010
//
//           Copyright ¨Ï 2009 by Voronoi Diagram Research Center, Hanyang University
//
/////////////////////////////////////////////////////////////////////////////////////

#ifndef _SQUARED_DISTFUNC_LINESEG_ARC_H_
#define _SQUARED_DISTFUNC_LINESEG_ARC_H_

#include "SquaredDistFuncLineSegCircle.h"
#include "SquaredDistFuncLineSegPoint.h"
#include "Arc3D.h"
#include "TrisectorOfSpaceByArc.h"

class SquaredDistFuncLineSegArc : public DistFunc
{
private:
	rg_INT      m_numOfIntervals;
	DistFunc**  m_sDistFunc;	
	rg_REAL*    m_intervalPoint;
	
public:
	SquaredDistFuncLineSegArc();
	SquaredDistFuncLineSegArc(const SquaredDistFuncLineSegArc& sDistFunc);
	SquaredDistFuncLineSegArc(const LineSegment3D& lineSegment, 
		                      const Arc3D& arc, 
							  const rg_REAL& lower = 0.0,
							  const rg_REAL& upper = 1.0);
	~SquaredDistFuncLineSegArc();

	void setIntervalPts_sDistFuncs(const rg_INT& numOfIntervals, 
		                           const rg_REAL* intervalPoint,
								   DistFunc** sDistFunc);

	rg_INT computeLocalMinima(rg_dList<rg_REAL>& localMinima);		
		//LocalPtStatus isThisPtLocalMinOrMax(const rg_REAL& point, const rg_REAL& delta_neighborhood) const;		
		LocalPtStatus isThisIntervalPtLocalMinOrMax(const rg_REAL& point) const;
		rg_REAL computeDeltaNeighborhood(rg_dList<rg_REAL>& localMinima) const;
		LocalPtStatus isThisPtLocalMinOrMax(const rg_REAL& point) const;
		LocalPtStatus isThisIntervalPtLocalMinOrMax(const rg_REAL& point, const rg_FLAG& isStartOrEndIntervalPt) const;
	rg_REAL evaluate(const rg_REAL& point) const;
	//rg_REAL evaluateDerivative(const rg_REAL& point) const;
	//rg_REAL evaluateDerivativeOfDerivative(const rg_REAL& point) const;

	rg_INT computeInflectionPoints(rg_dList<rg_REAL>& inflectionPts) const;
	rg_INT computeRootOfThirdDerivative(rg_dList<rg_REAL>& root) const;

	SquaredDistFuncLineSegArc& operator=(const SquaredDistFuncLineSegArc& sDistFunc);
};

#endif
