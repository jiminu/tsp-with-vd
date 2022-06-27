#ifndef _RG_BZINTERSECTOR_H
#define _RG_BZINTERSECTOR_H

#include "rg_Const.h"
#include "rg_BzCurve2D.h"
#include "rg_RQBzCurve2D.h"
#include "rg_Polynomial.h"
#include "rg_ComplexNumber.h"
#include "rg_Point2D.h"
#include "rg_Line.h"
#include "rg_dList.h"

class rg_BzIntersector
{
public:
	rg_BzIntersector();
	~rg_BzIntersector();

	rg_dList<rg_Point2D>  intersectBzCurveVsBzCurve(const rg_BzCurve2D &curve_s, 
                                              const rg_BzCurve2D &curve_t, 
                                              rg_dList<rg_REAL*> &seedParam4TwoCurve);

	rg_dList<rg_Point2D>  intersectBzCurveVsBzCurve(const rg_BzCurve2D &curve_s, 
                                              const rg_BzCurve2D &curve_t, 
                                              rg_dList<rg_REAL*> &seedParam4TwoCurve,
											  rg_REAL & time);

	rg_dList<rg_REAL*>  makeSeed(rg_dList<rg_RQBzCurve2D> &rqcurve_s, 
                             rg_dList<rg_REAL> subParam_s,
		                     rg_dList<rg_RQBzCurve2D> &rqcurve_t, 
                             rg_dList<rg_REAL> subParam_t);
	rg_Point2D          iterationWithSeed(const rg_BzCurve2D &curve_s, const rg_REAL &param_s,
		                               const rg_BzCurve2D &curve_t, const rg_REAL &param_t);
	rg_dList<rg_RQBzCurve2D> approximateRQBzCurves(const rg_BzCurve2D &curve, 
                                           rg_dList<rg_REAL> &param);
	rg_RQBzCurve2D        makeOneRQBzCurve(const rg_BzCurve2D &curve, 
                                      const rg_REAL &t0, 
                                      const rg_REAL &t1);
	rg_dList<rg_REAL>    findCharParam(const rg_BzCurve2D &curve);
    void             makeSimpleParam(const rg_BzCurve2D& curve, 
                                     const rg_REAL& t0, 
                                     const rg_REAL& t1,
                                     rg_dList<rg_REAL>& charParam);
    rg_REAL*          getConjugateTangentParam(const rg_BzCurve2D& curve, 
                                              const rg_REAL& t0,
                                              const rg_REAL& t1);

};

#endif // rg_BzIntersector class


