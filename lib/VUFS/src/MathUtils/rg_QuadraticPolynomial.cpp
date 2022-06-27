#include <math.h>
#include "rg_QuadraticPolynomial.h"
#include "rg_MathFunc.h"
#include "rg_RelativeOp.h"

rg_QuadraticPolynomial::rg_QuadraticPolynomial()
{
	degree = rg_QUADRATIC;
}

rg_QuadraticPolynomial::rg_QuadraticPolynomial(rg_REAL* coeff)
: rg_Polynomial(rg_QUADRATIC, coeff)
{
}

rg_QuadraticPolynomial::rg_QuadraticPolynomial(const rg_QuadraticPolynomial& qPolynomial)
: rg_Polynomial(qPolynomial)
{
}

rg_QuadraticPolynomial::~rg_QuadraticPolynomial()
{

}

void rg_QuadraticPolynomial::setPolynomial(rg_REAL *coeff)
{
	degree = rg_QUADRATIC;
	rg_DEGREE order = degree + 1;

	if(coefficient == rg_NULL)
		coefficient = new rg_REAL[order];

	for(rg_INT i = 0;i < order;i++)
		coefficient[ i ] = coeff[ i ];
}

rg_ComplexNumber* rg_QuadraticPolynomial::solve()
{	
	if(!rg_ZERO(coefficient[ 2 ]))
	{
		rg_ComplexNumber* r = new rg_ComplexNumber[rg_QUADRATIC];
		rg_REAL p = coefficient[ 1 ] / (2.0 * coefficient[ 2 ]);
		rg_REAL q = coefficient[ 0 ] / coefficient[ 2 ];
		
		rg_REAL determinant =  p * p - q;
		// determinant >= 0
		if(rg_NNEG(determinant))
		{
			r[ 0 ].setRealNumber(-p + sqrt(determinant));
			r[ 1 ].setRealNumber(-p - sqrt(determinant));
		}
		// determinant < 0
		else
		{
			r[ 0 ].setRealNumber(-p);
			r[ 0 ].setImaginaryNumber(sqrt(-determinant));
			r[ 1 ].setRealNumber(-p);
			r[ 1 ].setImaginaryNumber(-sqrt(-determinant));
		}
		return r;
	}
	// coefficient[ 2 ] == 0
	else
	{
		rg_ComplexNumber* r = new rg_ComplexNumber[rg_LINEAR];
		r[ 0 ].setRealNumber(- coefficient[ 0 ] / coefficient[ 1 ]);
		return r;
	}
}


