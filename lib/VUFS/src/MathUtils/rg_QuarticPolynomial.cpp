/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : QuatricPolynomial.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class QuatricPolynomial 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jan 1997    
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include <math.h>
#include "rg_CubicPolynomial.h"
#include "rg_QuarticPolynomial.h"
#include "rg_MathFunc.h"
#include "rg_RelativeOp.h"

rg_QuarticPolynomial::rg_QuarticPolynomial()
{
	degree = 4;
}

// If a cubic polynomial, for example, is 2x^3 - 4x^2 + x^1 + 3, then 
//
//            coefficient[0] =  3, 
//            coefficient[1] =  1,
//            coefficient[2] = -4,
//        and coefficient[3] =  2.

rg_QuarticPolynomial::rg_QuarticPolynomial(rg_REAL *coeff)
: rg_Polynomial(4, coeff)
{
}

rg_QuarticPolynomial::rg_QuarticPolynomial(const rg_QuarticPolynomial &eq)
: rg_Polynomial(eq.degree, eq.coefficient)
{
}

rg_QuarticPolynomial::~rg_QuarticPolynomial()
{
}

void rg_QuarticPolynomial::setQuarticPolynomial(rg_REAL *coeff)
{
	for(rg_INT i = 0; i < 4+1; i++)
		coefficient[i] = coeff[i];

}

rg_ComplexNumber *rg_QuarticPolynomial::solve()
{
	rg_REAL *a = new rg_REAL[4];
	if( rg_NZERO(coefficient[4]) )
	{
		a[3] = coefficient[3]/coefficient[4];
		a[2] = coefficient[2]/coefficient[4];
		a[1] = coefficient[1]/coefficient[4];
		a[0] = coefficient[0]/coefficient[4];
	}

	rg_REAL p = -3.0/8.0*pow(a[3], 2) + a[2];
	rg_REAL q =  1.0/8.0*pow(a[3], 3) - 1.0/2.0*a[3]*a[2] + a[1];
	rg_REAL r = -3.0/256.0*pow(a[3], 4) + 1.0/16.0*pow(a[3], 2)*a[2] - 1.0/4.0*a[3]*a[1] + a[0];

	// Handling special cases
	// q == 0 (Biquadratic eq)
	if( rg_ZERO( q ) )
	{
		rg_ComplexNumber *x = new rg_ComplexNumber[4];
		rg_ComplexNumber determin_over_2(p * p / 4.0 - r, 0.0);
		rg_ComplexNumber cCoeff_over_4_times_qCoeff(coefficient[3] / (4.0 * coefficient[4]), 0.0);
		rg_ComplexNumber minus_p_over_2(-p/2.0, 0.0);

		rg_ComplexNumber *r = rg_MathFunc::root(determin_over_2, 2);
		rg_ComplexNumber y[ 2 ] = {minus_p_over_2 + r[ 0 ], minus_p_over_2 + r[ 1 ]};
		delete [] r;

		r = rg_MathFunc::root(y[ 0 ], 2);
		x[ 0 ] =  r[ 0 ] - cCoeff_over_4_times_qCoeff;
		x[ 1 ] =  r[ 1 ] - cCoeff_over_4_times_qCoeff;
		delete [] r;

		r = rg_MathFunc::root(y[ 1 ], 2);
		x[ 2 ] =  r[ 0 ] - cCoeff_over_4_times_qCoeff;
		x[ 3 ] =  r[ 1 ] - cCoeff_over_4_times_qCoeff;
		delete [] r;

		return x;
	}
	// r == 0
	if( rg_ZERO( r ) )
	{
		rg_ComplexNumber *x = new rg_ComplexNumber[4];
		rg_REAL *coeff3 = new rg_REAL[4];
		coeff3[0] =  q;
		coeff3[1] =  p;
		coeff3[2] =  0.0;
		coeff3[3] =  1.0;
		rg_CubicPolynomial cubicEq(coeff3);
		rg_ComplexNumber *y = cubicEq.solve();

		x[ 0 ] = rg_ComplexNumber(0.0, 0.0) - rg_ComplexNumber(0.25 * a[3], 0.0);
		for(rg_INT i = 0;i < 3;i++)
			x[i + 1] = y[ i ] - rg_ComplexNumber(0.25 * a[3], 0.0);
		
		delete [] y;
		return x;
	}
	/////////////////////////

	rg_REAL *coeff3 = new rg_REAL[4];
	coeff3[0] =  4*p*r - pow(q, 2);
	coeff3[1] = -4*r;
	coeff3[2] = -p;
	coeff3[3] =  1.0;

	rg_CubicPolynomial cubicEq(coeff3);
	rg_ComplexNumber *u = cubicEq.solve();
	delete []coeff3;

	rg_ComplexNumber *z      = rg_MathFunc::root(u[0]-p, 2);
	rg_ComplexNumber *coeff2 = new rg_ComplexNumber[2];

	rg_INT i;
	rg_ComplexNumber *x = new rg_ComplexNumber[4];
	for(i = 0; i < 2; i++)
	{
		coeff2[0] = -pow(-1., i) * z[0];
		coeff2[1] = u[0]/2.0 + pow(-1., i) * q / 2.0 / z[0];

		rg_ComplexNumber *tempx = solveQuadraticEq(coeff2);
		x[i*2]   = tempx[0] - a[3]/4.0;
		x[i*2+1] = tempx[1] - a[3]/4.0;
		delete []tempx;
	}

    delete[] coeff2;
	delete[] u;
	delete[] z;
//  delete[] a;
/*
	for(i = 0; i < 4; i++)
		TRACE(" x[%d] = %5.3f + %5.3fi\n", i, x[i].getRealNumber(), x[i].getImaginaryNumber());

	// This is a part of testing result values.
	rg_INT okCount = 0;
	for(i = 0; i < 4; i++)
	{
		rg_ComplexNumber vTest(x[i]*x[i]*x[i]*x[i] + a[3]*x[i]*x[i]*x[i] + a[2]*x[i]*x[i] + a[1]*x[i] + a[0]);
		if(vTest.isZero()) okCount++;
	}

	if(okCount == 4)
		TRACE("value test is good :quartic case\n");
	else 
		return rg_NULL;
*/
    rg_INT okCount = 0;
	for(i = 0; i < 4; i++)
	{
		rg_ComplexNumber vTest(x[i]*x[i]*x[i]*x[i] + a[3]*x[i]*x[i]*x[i] + a[2]*x[i]*x[i] + a[1]*x[i] + a[0]);
		if(vTest.isZero()) okCount++;
	}

	//if(okCount != 4) return rg_NULL;

	return x;
}


rg_ComplexNumber *rg_QuarticPolynomial::solve_OLD()
{
	rg_REAL *a = new rg_REAL[4];
	if (rg_NZERO(coefficient[4]))
	{
		a[3] = coefficient[3] / coefficient[4];
		a[2] = coefficient[2] / coefficient[4];
		a[1] = coefficient[1] / coefficient[4];
		a[0] = coefficient[0] / coefficient[4];
	}

	rg_REAL p = -3.0 / 8.0*pow(a[3], 2) + a[2];
	rg_REAL q = 1.0 / 8.0*pow(a[3], 3) - 1.0 / 2.0*a[3] * a[2] + a[1];
	rg_REAL r = -3.0 / 256.0*pow(a[3], 4) + 1.0 / 16.0*pow(a[3], 2)*a[2] - 1.0 / 4.0*a[3] * a[1] + a[0];

	rg_REAL *coeff3 = new rg_REAL[4];
	coeff3[0] = 4 * p*r - pow(q, 2);
	coeff3[1] = -4 * r;
	coeff3[2] = -p;
	coeff3[3] = 1.0;

	rg_CubicPolynomial cubicEq(coeff3);
	rg_ComplexNumber *u = cubicEq.solve();
	delete[]coeff3;

	rg_ComplexNumber *z = rg_MathFunc::root(u[0] - p, 2);
	rg_ComplexNumber *coeff2 = new rg_ComplexNumber[2];

	rg_ComplexNumber *x = new rg_ComplexNumber[4];
	rg_INT i = 0;
	for (i = 0; i < 2; i++)
	{
		coeff2[0] = -pow(-1.0, i) * z[0];
		coeff2[1] = u[0] / 2.0 + pow(-1.0, i) * q / 2.0 / z[0];

		rg_ComplexNumber *tempx = solveQuadraticEq(coeff2);
		x[i * 2] = tempx[0] - a[3] / 4.0;
		x[i * 2 + 1] = tempx[1] - a[3] / 4.0;
		delete[]tempx;
	}

	delete[] coeff2;
	delete[] u;
	delete[] z;
	//  delete[] a;
	/*
	for(i = 0; i < 4; i++)
	TRACE(" x[%d] = %5.3f + %5.3fi\n", i, x[i].getRealNumber(), x[i].getImaginaryNumber());

	// This is a part of testing result values.
	rg_INT okCount = 0;
	for(i = 0; i < 4; i++)
	{
	rg_ComplexNumber vTest(x[i]*x[i]*x[i]*x[i] + a[3]*x[i]*x[i]*x[i] + a[2]*x[i]*x[i] + a[1]*x[i] + a[0]);
	if(vTest.isZero()) okCount++;
	}

	if(okCount == 4)
	TRACE("value test is good :quartic case\n");
	else
	return rg_NULL;
	*/
	rg_INT okCount = 0;
	for (i = 0; i < 4; i++)
	{
		rg_ComplexNumber vTest(x[i] * x[i] * x[i] * x[i] + a[3] * x[i] * x[i] * x[i] + a[2] * x[i] * x[i] + a[1] * x[i] + a[0]);
		if (vTest.isZero()) okCount++;
	}

	if(okCount != 4) return rg_NULL;

	return x;
}


rg_ComplexNumber *rg_QuarticPolynomial::solveQuadraticEq(rg_ComplexNumber *coeff2)
{
	rg_ComplexNumber *x = rg_MathFunc::root(coeff2[0]*coeff2[0] - 4*coeff2[1], 2);
	
	for(rg_INT i=0; i<2; i++)
		x[i] = (-coeff2[0] + x[i])/2.0;

/*
	// This is a part of testing result values.
	rg_INT okCount = 0;
	for(i = 0; i < 2; i++)
	{
		rg_ComplexNumber vTest(x[i]*x[i] + coeff2[0]*x[i] + coeff2[1]);
		if(vTest.isZero()) okCount++;
	}

	if(okCount == 2)
		TRACE("value test is good\n");
*/
	return x;
}


