#ifndef _RG_NURBSCURVEINTERSECTOR_H
#define _RG_NURBSCURVEINTERSECTOR_H

#include "rg_Const.h"
#include "rg_BzCurve2D.h"
#include "rg_RQBzCurve2D.h"
#include "rg_Polynomial.h"
#include "rg_ComplexNumber.h"
#include "rg_Point2D.h"
#include "rg_Line.h"
#include "rg_dList.h"

#include "rg_NUBSplineCurve3D.h"
#include "rg_NURBSplineCurve3D.h"

#include "rg_BzIntersector.h"

class rg_NURBSCurveIntersector
{
public:
	rg_NURBSCurveIntersector();
	~rg_NURBSCurveIntersector();
/*

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Decompose the given NURBS curve into Rational Bezier curves and 
// then apply Cocktail Algorithm to each pair of Bezier rg_Curve !!(filtering process not yet implemented)


	rg_dList<rg_Point2D> intersectBSplineCurveVsBSplineCurveUsingCurveDecomposition(rg_NURBSplineCurve3D& curve_s,
																			  rg_NURBSplineCurve3D& curve_t,
																			  rg_dList<rg_REAL*> & seedParamList4TwoCurve,
																			  rg_REAL& time);
	rg_dList<rg_Point3D>    findCharPointUsingCurveDecomposition(rg_NURBSplineCurve3D &curve);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// make the each curve segment of given NURBS curve polynomial form 
// by making the basis functions polynomial and
// then apply Cocktail Algorithm(find inflection point from polynomial form curve)

	rg_dList<rg_Point2D> intersectNURBSplineCurves(rg_NURBSplineCurve3D& curve_s, 
								        	 rg_NURBSplineCurve3D& curve_t,
										     rg_dList<rg_REAL*> & parameters, rg_REAL & time);

	rg_dList<rg_REAL>    findCharParam(rg_NURBSplineCurve3D &curve);

	rg_dList<rg_RQBzCurve2D> approximateRQBzCurves(rg_NURBSplineCurve3D &curve, rg_dList<rg_REAL> &param);

	rg_RQBzCurve2D        makeOneRQBzCurve(rg_NURBSplineCurve3D &curve, const rg_REAL &t0, const rg_REAL &t1);

	rg_dList<rg_REAL*>   makeSeed(rg_dList<rg_RQBzCurve2D> &rqcurve_s, rg_dList<rg_REAL> subParam_s, rg_dList<rg_RQBzCurve2D> &rqcurve_t, rg_dList<rg_REAL> subParam_t);

	rg_Point2D          iterationWithSeed( rg_NURBSplineCurve3D &curve_s, const rg_REAL &param_s,
		                                rg_NURBSplineCurve3D &curve_t, const rg_REAL &param_t);
    rg_Point2D          iterationWithParameter( rg_NURBSplineCurve3D &curve_s, rg_REAL &param_s,
					    		             rg_NURBSplineCurve3D &curve_t, rg_REAL &param_t);


    void             makeSimpleParam(rg_NURBSplineCurve3D& curve, const rg_REAL& t0, const rg_REAL& t1, rg_dList<rg_REAL>& charParam);

    rg_REAL*          getConjugateTangentParam(rg_NURBSplineCurve3D& curve, const rg_REAL& t0, const rg_REAL& t1);

	/*
	rg_dList<rg_Point2D>  intersectBSplineCurveVsBSplineCurve(const  &curve_s,
														const  &curve_t, 
														rg_dList<rg_REAL*> &seedParam4TwoCurve);
											  */

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Decompose the given NURBS curve into Rational Bezier curves and 
// then apply Cocktail Algorithm to entire NURBS curve only one time !!
// (obtain the info. about inflection and derivative from each Bezier curves)


	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


};

#endif


