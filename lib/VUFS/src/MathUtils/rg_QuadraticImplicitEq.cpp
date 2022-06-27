#include "rg_QuadraticImplicitEq.h"

rg_QuadraticImplicitEq::rg_QuadraticImplicitEq()
{
}

rg_QuadraticImplicitEq::rg_QuadraticImplicitEq(const rg_REAL &A, const rg_REAL &B, const rg_REAL &C, const rg_REAL &D, const rg_REAL &E, const rg_REAL &F)
{
	a = A;
	b = B;
	c = C;
	d = D;
	e = E;
	f = F;
}

rg_QuadraticImplicitEq::rg_QuadraticImplicitEq(const rg_QuadraticImplicitEq &equation)
{
	a = equation.a;
	b = equation.b;
	c = equation.c;
	d = equation.d;
	e = equation.e;
	f = equation.f;
}

rg_QuadraticImplicitEq::~rg_QuadraticImplicitEq()
{
}


