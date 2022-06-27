#include "rg_CubicImplicitEq.h"

rg_CubicImplicitEq::rg_CubicImplicitEq()
{
}

rg_CubicImplicitEq::rg_CubicImplicitEq(const rg_REAL &A, const rg_REAL &B, const rg_REAL &C, 
								 const rg_REAL &D, const rg_REAL &E, const rg_REAL &F,
								 const rg_REAL &G, const rg_REAL &H, const rg_REAL &I, const rg_REAL &J)
{
	a = A;
	b = B;
	c = C;
	d = D;
	e = E;
	f = F;
	g = G;
	h = H;
	i = I;
	j = J;
}

rg_CubicImplicitEq::rg_CubicImplicitEq(const rg_CubicImplicitEq &equation)
{
	a = equation.a;
	b = equation.b;
	c = equation.c;
	d = equation.d;
	e = equation.e;
	f = equation.f;
	g = equation.g;
	h = equation.h;
	i = equation.i;
	j = equation.j;
}

rg_CubicImplicitEq::~rg_CubicImplicitEq()
{
}


