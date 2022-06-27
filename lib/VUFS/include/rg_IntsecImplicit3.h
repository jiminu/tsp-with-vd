#ifndef _RG_INTSECIMPLICT3_H
#define _RG_INTSECIMPLICT3_H
		  
#include "rg_Const.h"
#include "rg_CBzCurve2D.h"
#include "rg_Polynomial.h"
#include "rg_ImplicitEquation.h"
#include "rg_dList.h"
#include "rg_Point2D.h"
#include "rg_ComplexNumber.h"
#include "rg_Polynomial.h"

class rg_IntsecImplicit3
{
public:
	rg_IntsecImplicit3();
	~rg_IntsecImplicit3();
	rg_dList<rg_Point2D>   intersectBzCurveVsBzCurve(const rg_BzCurve2D &curve_s, 
											   const rg_BzCurve2D &curve_t, 
											   rg_REAL &time);
	rg_ImplicitEquation implicitize(const rg_BzCurve2D &curve);
	rg_REAL           inversion(const rg_BzCurve2D &curve, const rg_Point2D &point);
};

#endif // rg_IntsecImplicit3 class

