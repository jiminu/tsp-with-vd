////////////////////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : LineSegment3D.h
//	  
//    DESCRIPTION : 
//           This is the implementation of the class LineSegment3D
//           and is used for representing a line.
//
//    AUTHOR      : Ryu, Jooghyun
//    START DATE  : Nov. 10, 2008   
//
//           Copyright ¨Ï 2009 by Voronoi Diagram Research Center, Hanyang University
//
/////////////////////////////////////////////////////////////////////////////////////

#ifndef _LINESEGMENT3D_H_
#define _LINESEGMENT3D_H_

#include "Line3D.h"
#include "rg_RelativeOp.h"

enum PosOnLineSegOrArc {START_PT, END_PT, NONEXTREME_PT, UNKNOWN_POS};

class LineSegment3D : public Line3D
{
private:
	rg_Point3D m_sPoint;
	rg_Point3D m_ePoint;

public:
	LineSegment3D();
	LineSegment3D(const rg_Point3D& sPt, const rg_Point3D& ePt);
	LineSegment3D(const LineSegment3D& lineSegment);
	~LineSegment3D();

	rg_Point3D getStartPt() const;
	rg_Point3D getEndPt()   const;
	rg_REAL    getLength()  const;

    void       setPoints(const rg_Point3D& sPt, const rg_Point3D& ePt);

    LineSegment3D& operator =(const LineSegment3D& lineSegment);

    rg_BOOL isOnLineSegment(const rg_Point3D& point, const rg_REAL& res = rg_MATH_RES) const;


	rg_Point3D evaluatePt(const rg_REAL& param) const;
	rg_Point3D computeProjectedPointForPoint3D(const rg_Point3D& target, rg_REAL& param) const;

    rg_REAL computeMinDistFromPoint(const rg_Point3D& targetPt) const;
	rg_REAL computeMinDistFromPoint(const rg_Point3D& targetPt, 
									rg_Point3D& minPt, 
									rg_REAL& minPtParam,
									PosOnLineSegOrArc& pos) const;
	rg_REAL computeMinDistFromPoint(const rg_Point3D& targetPt, 
		                            rg_Point3D& minPt, 
		                            rg_REAL& minPtParam) const;
	rg_REAL computeMaxDistFromPoint(const rg_Point3D& targetPt, 
		                            rg_Point3D& maxPt,
									PosOnLineSegOrArc& pos) const;
};

#endif

