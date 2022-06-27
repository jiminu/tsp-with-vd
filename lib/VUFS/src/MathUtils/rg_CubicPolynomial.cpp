/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_CubicPolynomial.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_CubicPolynomial 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jan 1997    
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include <math.h>
#include "rg_CubicPolynomial.h"
#include "rg_MathFunc.h"
#include "rg_RelativeOp.h"

rg_CubicPolynomial::rg_CubicPolynomial()
{
    degree = rg_CUBIC;
}

// If a cubic polynomial, for example, is 2x^3 - 4x^2 + x^1 + 3, then 
//
//            coefficient[0] =  3, 
//            coefficient[1] =  1,
//            coefficient[2] = -4,
//        and coefficient[3] =  2.

rg_CubicPolynomial::rg_CubicPolynomial(rg_REAL *coeff)
: rg_Polynomial(rg_CUBIC, coeff)
{
}

rg_CubicPolynomial::rg_CubicPolynomial(const rg_CubicPolynomial &eq)
: rg_Polynomial(eq.degree, eq.coefficient)
{
}

rg_CubicPolynomial::~rg_CubicPolynomial()
{
}

void rg_CubicPolynomial::setPolynomial(rg_REAL *coeff)
{
	for(rg_INDEX i = 0; i < rg_CUBIC+1; i++)
		coefficient[i] = coeff[i];
}

rg_ComplexNumber* rg_CubicPolynomial::solve()
{
	rg_REAL *a = new rg_REAL[rg_CUBIC];
	if( rg_NZERO(coefficient[3]) ) 
	{
		a[2] = coefficient[2]/coefficient[3];
		a[1] = coefficient[1]/coefficient[3];
		a[0] = coefficient[0]/coefficient[3];
	}

	rg_REAL p = a[1] - pow(a[2], 2) / 3.0;
	rg_REAL q = 2 * pow(a[2], 3) / 27.0 - a[2] * a[1] / 3.0 + a[0];

	rg_REAL A = -q / 2.0;
	rg_REAL B = pow(p, 3) / 27.0 + pow(q, 2) / 4.0;

	// w = -1/2 + rg_MathFunc::root(3)/2 i
	rg_ComplexNumber omega(-0.5, sqrt(3.)/2.0);
	rg_ComplexNumber z1;
	rg_ComplexNumber z2;

	if( rg_POS(B, resNeg19) ) 
	{
		z1.setRealNumber( rg_MathFunc::cubicRoot(A+sqrt(B)) );
		z2.setRealNumber( rg_MathFunc::cubicRoot(A-sqrt(B)) );
	}
	
	else 
	{
		if( rg_ZERO(A, resNeg19) )
		{
            rg_REAL cubicRootOfB=-rg_MathFunc::cubicRoot(B);
            rg_REAL temp=sqrt(cubicRootOfB);
			z1.setImaginaryNumber( temp );
			z2.setImaginaryNumber(-temp );
		}

		else 
		{
            rg_REAL rootOfB_=sqrt(-B);
			rg_ComplexNumber other1(A, rootOfB_);
			rg_ComplexNumber *root1 = rg_MathFunc::cubicRoot(other1);
	
			rg_ComplexNumber other2(A, -rootOfB_);
			rg_ComplexNumber *root2 = rg_MathFunc::cubicRoot(other2);

			rg_ComplexNumber copyP(-p/3, 0);
			for(rg_INDEX i=0; i<3; i++)
			{
				for(rg_INDEX j=0; j<3; j++)
				{
					if(root1[i]*root2[j] == copyP) 
					{
						z1 = root1[i];
						z2 = root2[j];
					}
				}
			}
			delete []root1;
			delete []root2;
		}
	}

	rg_ComplexNumber *x = new rg_ComplexNumber[rg_CUBIC];

	x[0] = z1 + z2 - a[2] / 3.0;
	x[1] = omega * z1 + omega * omega * z2 - a[2] / 3.0;
	x[2] = omega * omega * z1 + omega * z2 - a[2] / 3.0;

    delete[] a;
    // We display solutions of a cubic equation.
/*
	for(rg_INDEX i = 0; i < 3; i++)
		TRACE(" x[%d] = %5.3f + %5.3fi\n", i, x[i].getRealNumber(), x[i].getImaginaryNumber());

	// This is a part of testing result values.
	rg_INT okCount = 0;
	for(i = 0; i < 3; i++)
	{
		rg_ComplexNumber vTest(x[i]*x[i]*x[i] + a[2]*x[i]*x[i] + a[1]*x[i] + a[0]);
		if(vTest.isZero()) okCount++;
	}

	if(okCount == 3)
		TRACE("value test is good :cubic case\n");
*/
    return x;
}


