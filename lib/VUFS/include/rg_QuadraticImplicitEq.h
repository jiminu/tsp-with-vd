#ifndef _RG_QUADRATICIMPLICITEQ_H
#define _RG_QUADRATICIMPLICITEQ_H
#include "rg_Const.h"
class rg_QuadraticImplicitEq
{
public:
	rg_REAL a;
	rg_REAL b;
	rg_REAL c;
	rg_REAL d;
	rg_REAL e;
	rg_REAL f;

public:
	rg_QuadraticImplicitEq();
	rg_QuadraticImplicitEq(const rg_REAL &A, const rg_REAL &B, const rg_REAL &C, const rg_REAL &D, const rg_REAL &E, const rg_REAL &F);
	rg_QuadraticImplicitEq(const rg_QuadraticImplicitEq &eqation);
	~rg_QuadraticImplicitEq();
};

#endif


