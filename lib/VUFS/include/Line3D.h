#ifndef _LINE3D_H_
#define _LINE3D_H_

#include "rg_Point3D.h"

// parametric representation L(t) of an infinite line
// with a direction vector V(v_x, v_y, v_z)
//      a passing point    P0(x0, y0, z0)
//
//  L(t) = P0 + tV

class Line3D
{
private:
	rg_Point3D m_directionVec;
	rg_Point3D m_passingPt;

public:
	Line3D();
	Line3D(const rg_Point3D& directionVec, const rg_Point3D& passingPt);
	Line3D(const Line3D& line);
	~Line3D();

	rg_Point3D getDirVec() const;
	rg_Point3D getPassPt() const;

	void       set(const rg_Point3D& directionVec, const rg_Point3D& passingPt);
	void       setDirVec(const rg_Point3D& directionVec);
	void       setPassPt(const rg_Point3D& passingPt);

	Line3D& operator=(const Line3D& line);

    rg_BOOL    isOnLine(const rg_Point3D& point, const rg_REAL& res = rg_MATH_RES) const;

	rg_Point3D computeProjectedPointForPoint3D(const rg_Point3D& point) const;
	rg_Point3D computeProjectedPointForPoint3D(const rg_Point3D& point, rg_REAL& dist) const;

	rg_REAL    computeDistanceFromPoint3D(  const rg_Point3D& point) const;
	rg_REAL    computeMinDistFromPoint(     const rg_Point3D& targetPt, 
		                                    rg_Point3D&       minPt, 
									        rg_REAL&          minPtParam) const;

};

#endif

