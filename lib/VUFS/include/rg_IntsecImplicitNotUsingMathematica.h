#ifndef _INTSECIMPLICITNOTUSINGMATHEMATICA_H
#define _INTSECIMPLICITNOTUSINGMATHEMATICA_H
		  
#include "rg_Const.h"
#include "rg_CBzCurve2D.h"
#include "rg_Polynomial.h"
#include "rg_ImplicitEquation.h"
#include "rg_dList.h"
#include "rg_Point2D.h"
#include "rg_ComplexNumber.h"
#include "rg_Polynomial.h"

class rg_IntsecImplicitNotUsingMathematica
{
public:
	rg_IntsecImplicitNotUsingMathematica();
	~rg_IntsecImplicitNotUsingMathematica();
	rg_dList<rg_Point2D>   intersectBzCurveVsBzCurve1(const rg_BzCurve2D &curve_s, 
											   const rg_BzCurve2D &curve_t, 
											   rg_REAL &time);

	rg_dList<rg_Point2D>   intersectBzCurveVsBzCurve2(const rg_BzCurve2D &curve_s, 
											    const rg_BzCurve2D &curve_t, 
											    rg_REAL &time);

	rg_ImplicitEquation implicitize1(const rg_BzCurve2D &curve);
	rg_ImplicitEquation implicitize2(const rg_BzCurve2D &curve);

	rg_REAL           inversion(const rg_BzCurve2D &curve, const rg_Point2D &point);
	rg_REAL           inversionUsingLinearEquation(const rg_BzCurve2D &curve, const rg_Point2D &point);
};

#endif // rg_IntsecImplicitNotUsingMathematica class


