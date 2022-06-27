#include <math.h>
#include "rg_LinearPolynomial.h"
#include "rg_MathFunc.h"
#include "rg_RelativeOp.h"

rg_LinearPolynomial::rg_LinearPolynomial()
{
	degree = rg_LINEAR;
}

rg_LinearPolynomial::rg_LinearPolynomial(rg_REAL* coeff)
: rg_Polynomial(rg_LINEAR, coeff)
{
}

rg_LinearPolynomial::rg_LinearPolynomial(const rg_LinearPolynomial& lPolynomial)
: rg_Polynomial(lPolynomial)
{
}

rg_LinearPolynomial::~rg_LinearPolynomial()
{
	
}

void rg_LinearPolynomial::setPolynomial(rg_REAL *coeff)
{
	degree = rg_LINEAR;
	rg_DEGREE order = degree + 1;
	
	if(coefficient == rg_NULL)
		coefficient = new rg_REAL[order];
	
	for(rg_INT i = 0;i < order;i++)
		coefficient[ i ] = coeff[ i ];
}

rg_ComplexNumber* rg_LinearPolynomial::solve()
{
	rg_ComplexNumber* r = new rg_ComplexNumber[rg_LINEAR];
	r[ 0 ].setRealNumber(-coefficient[ 0 ] / coefficient[ 1 ]);
	return r;
}


