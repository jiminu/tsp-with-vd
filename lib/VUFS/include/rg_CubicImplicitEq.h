#ifndef _RG_CUBICIMPLICITEQ_H
#define _RG_CUBICIMPLICITEQ_H
#include "rg_Const.h"
class rg_CubicImplicitEq
{
public:
	rg_REAL a;
	rg_REAL b;
	rg_REAL c;
	rg_REAL d;
	rg_REAL e;
	rg_REAL f;
	rg_REAL g;
	rg_REAL h;
	rg_REAL i;
	rg_REAL j;

public:
	rg_CubicImplicitEq();
	rg_CubicImplicitEq(const rg_REAL &A, const rg_REAL &B, const rg_REAL &C, 
				 const rg_REAL &D, const rg_REAL &E, const rg_REAL &F,
				 const rg_REAL &G, const rg_REAL &H, const rg_REAL &I, const rg_REAL &J);
	rg_CubicImplicitEq(const rg_CubicImplicitEq &eqation);
	~rg_CubicImplicitEq();
};

#endif


