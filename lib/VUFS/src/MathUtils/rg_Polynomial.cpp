/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_Polynomial.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Polynomial 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jan 1997    
//
//    HISTORY     :
//          1. By Young-Song Cho,  13 Jul. 1997
//              make : rg_REAL  rg_Polynomial::getCoefficient(const rg_INDEX &i) const
//              make : rg_REAL* rg_Polynomial::getCoefficient() const
//              make : void  rg_Polynomial::setCoefficient(const rg_INDEX &i, 
//                                                      const rg_REAL &coeff)
//          2. By Hyung-Joo Lee,  15 Jan. 1998
//              make : void rg_Polynomial::setCoefficient( const rg_REAL* coeff )
//              make : rg_Polynomial rg_Polynomial::makeDerivative() const
//
//          3. By Hyung-Joo Lee,  16 Jan. 1998
//              make : void rg_Polynomial::updateDegree()
//                     - Check the coefficient of the highest degree term and
//                      , if that is zero, decrease the degree
//
//          4. By Jung-Hyun Ryu,  Feb 26. 1998
//              update : rg_REAL rg_Polynomial::getCoefficient(const rg_INDEX &i) const
//              
//          5. By Jung-Hyun Ryu,  Mar 6. 1998
//              make   : rg_Polynomial  operator ^(const rg_DEGREE index) const
//              make   : rg_REAL& rg_Polynomial::operator [](const rg_DEGREE index)
//
//          6. By Jung-Hyun Ryu,  July 14. 1998
//				make   : void		reconstruct();
//				make   : rg_Polynomial  operator /(const rg_Polynomial &operand) const;
//				make   : rg_Polynomial  operator /(const rg_REAL &n) const;
//
//          6. By Jung-Hyun Ryu,  July 15. 1998
//				make   : rg_REAL evaluatePolynomial(const rg_REAL& valueForVariable) const;
//
//          7. By Jung-Hyun Ryu,  July 28. 1998
//			    Update : void setDegree(const rg_DEGREE &dg);
//				Update : rg_Polynomial& operator =(const rg_Polynomial &polynomial);
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include <iostream>

#include <math.h>
#include "rg_Polynomial.h"
#include "rg_MathFunc.h"
#include "rg_RelativeOp.h"
#include "PolynomialSolver.h"

rg_Polynomial::rg_Polynomial()
{
	degree      = -1;
	coefficient = rg_NULL;
}

// If a cubic polynomial, for example, is 2x^3 - 4x^2 + x^1 + 3, then
//            coefficient[0] =  3, 
//            coefficient[1] =  1,
//            coefficient[2] = -4,
//        and coefficient[3] =  2.

rg_Polynomial::rg_Polynomial(const rg_DEGREE &dg)
{
	degree      = dg;
	coefficient = new rg_REAL[degree+1];
	for(rg_INT i=0; i<degree+1; i++)
		coefficient[i] = 0.0;
}

rg_Polynomial::rg_Polynomial(const rg_DEGREE &dg, rg_REAL *coeff)
{
	degree      = dg;
	coefficient = new rg_REAL[degree+1];
	for(rg_INT i=0; i<degree+1; i++)
		coefficient[i] = coeff[i];
}

rg_Polynomial::rg_Polynomial(const rg_Polynomial &eq)
{
	degree      = eq.degree;
	coefficient = new rg_REAL[degree+1];
	for(rg_INT i=0; i<degree+1; i++)
		coefficient[i] = eq.coefficient[i];
}

rg_Polynomial::~rg_Polynomial()
{
	if(coefficient != rg_NULL) delete []coefficient;
}

void rg_Polynomial::setPolynomial(const rg_DEGREE &dg, const rg_REAL *coeff)
{
	if(coefficient != rg_NULL) delete[] coefficient;
	
	degree      = dg;
	coefficient = new rg_REAL[degree+1];	
	for(rg_INT i=0; i<degree+1; i++)
		coefficient[i] = coeff[i];

}

void rg_Polynomial::setDegree(const rg_DEGREE &dg)
{
	if(coefficient != rg_NULL)
	{
		delete[] coefficient;
	}

	degree = dg;
	coefficient = new rg_REAL[degree+1];
	for(rg_INT i=0; i<degree+1; i++)
		coefficient[i] = 0.0;
/*
	degree      = dg;
	coefficient = new rg_REAL[degree+1];	
	for(rg_INT i=0; i<degree+1; i++)
		coefficient[i] = 0.0;
		*/
}

void rg_Polynomial::setCoefficient(const rg_INDEX &i, const rg_REAL &coeff)
{
    coefficient[i] = coeff;
}

// inserted in histrory by Hyung-Joo Lee
// function to find the first derivative of the polynomial
// The code for setCoefficient( const rg_REAL* coeff ) start here
void rg_Polynomial::setCoefficient( const rg_REAL* coeff )
{
    if ( coefficient )
        delete [] coefficient;

	coefficient = new rg_REAL[degree + 1];

	for( rg_INT i = 0; i < degree + 1; i++ )
	{
		coefficient[i] = coeff[i];
	}
}

// reduce the storage by deleting the all-zero-coefficient terms
// (note that this case may happen during the operation among the Polynomials!!)
void rg_Polynomial::reconstruct()
{
	rg_DEGREE newDegree = -1;
	for(rg_INT i = 0;i <= degree;i++)
	{
		if(!rg_ZERO(coefficient[ i ], rg_SYSTEM_RES))
			newDegree = i;
	}

	if(degree != newDegree)
	{
		degree = newDegree;
		rg_REAL* newCoeff = new rg_REAL[degree + 1];
		for(rg_INT i = 0;i <= degree;i++)
		{
			newCoeff[ i ] = coefficient[ i ];
		}
		setCoefficient(newCoeff);
		delete [] newCoeff;
	}
	/*
	rg_DEGREE newDegree = degree;
	rg_FLAG isConsecutive = rg_TRUE;

	rg_INT i = 0;
	for(i = degree;i >= 1;i--)
	{
		if(rg_ZERO(coefficient[ i ]) && isConsecutive)
			newDegree--;
		isConsecutive = rg_FALSE;
	}

	rg_Polynomial temp(newDegree);

	for(i = 0;i <= newDegree;i++)
		temp[ i ] = coefficient[ i ];

	delete[] coefficient;
	coefficient = rg_NULL;
	degree = 0;

	(*this) = temp;
	*/
}

// The code for setCoefficient( const rg_REAL* coeff ) start here

rg_DEGREE rg_Polynomial::getDegree() const
{
	return degree;
}

rg_REAL rg_Polynomial::getCoefficient(const rg_INDEX &i) const
{
	if((0 <= i) && (i <= degree))
	{
		return coefficient[i];
	}
	else
	{
		return 0.0;
	}
}

rg_REAL* rg_Polynomial::getCoefficient() const
{
    return coefficient;
}

rg_REAL rg_Polynomial::evaluatePolynomial(const rg_REAL& valueForVariable) const
{
	rg_REAL result = 0.0;
    /*
	for(rg_INDEX i = 0;i < degree + 1;i++)
		result += coefficient[ i ] * pow(valueForVariable, i);
    */
    // update by Taeboom, Jang
    for( rg_INDEX i=0; i < degree+1; i++ )
    {
        result = result* valueForVariable + coefficient[degree-i];
    }

	return result;
}

rg_Polynomial rg_Polynomial::operator *( const rg_REAL &n ) const
{
	rg_Polynomial result(degree);

	for(rg_INT i = 0; i < degree + 1; i++)
		result.coefficient[i] = n*coefficient[i];

	//result.reconstruct();

	return result;
}

rg_Polynomial rg_Polynomial::operator /(const rg_REAL &n) const
{
	rg_Polynomial result(degree);

	for(rg_INT i = 0; i < degree + 1; i++)
		result.coefficient[i] = coefficient[i]/n;

	//result.reconstruct();

	return result;	
}


rg_Polynomial operator *( const rg_REAL &n, const rg_Polynomial &polynomial )
{
	return polynomial * n;
}

rg_Polynomial rg_Polynomial::operator *(const rg_Polynomial &polynomial) const
{
	rg_Polynomial result(degree + polynomial.degree);

	for(rg_INT i = 0; i < degree + 1; i++)
	{
		for(rg_INT j = 0; j < polynomial.degree + 1; j++)
		{
			result.coefficient[i + j] += coefficient[i]*polynomial.coefficient[j];
		}
	}

	//result.reconstruct();

	return result;
}

rg_Polynomial rg_Polynomial::operator /(const rg_Polynomial &operand) const
{
	if(operand.degree > degree)
		return rg_Polynomial(0);
	
	else
	{
		if(operand.degree == 0)
		{
			return (*this) / operand.coefficient[0];
		}
		else
		{
			rg_Polynomial result, accumulated = (*this);

			rg_Polynomial temp(accumulated.degree - operand.degree);
			temp[accumulated.degree - operand.degree] = 1.0;
			result = temp = (accumulated.coefficient[accumulated.degree] / operand.coefficient[operand.degree])*temp;

			rg_Polynomial remainder =(accumulated - temp * operand);
			remainder.reconstruct();

			while(remainder.degree >= operand.degree )
			{
				accumulated = accumulated - temp * operand;
				accumulated.reconstruct();
				temp.setDegree(accumulated.degree - operand.degree);
				temp[accumulated.degree - operand.degree] = 1.0;
				temp = (accumulated.coefficient[accumulated.degree] / operand.coefficient[operand.degree])*temp;

				result = result + temp;
				remainder=(accumulated - temp * operand);
				remainder.reconstruct();
			}
				
			result.reconstruct();
			return result;
		}
	}
}

rg_Polynomial rg_Polynomial::operator ^(const rg_DEGREE index) const
{
	rg_Polynomial temp(0);
	temp[ 0 ] = 1.0;
	for(rg_INT i = 0;i < index;i++)
		temp = temp * (*this);
	
	return temp;
}

rg_Polynomial rg_Polynomial::operator +(const rg_Polynomial &polynomial) const
{
	rg_DEGREE maxDegree;
	if(degree > polynomial.degree)
	{
		maxDegree = degree;
	}
	else 
	{
		maxDegree = polynomial.degree;
	}

	rg_Polynomial result(maxDegree);

	for(rg_INT i = 0; i < maxDegree + 1; i++)
	{
		if( i > degree)
		{
			result.coefficient[i] = polynomial.coefficient[i];
		}
		else if( i > polynomial.degree)
		{
			result.coefficient[i] = coefficient[i];
		}
		else
		{
			result.coefficient[i] = coefficient[i] + polynomial.coefficient[i];
		}
	}

	//result.reconstruct();

	return result;
}

rg_Polynomial rg_Polynomial::operator -(const rg_Polynomial &polynomial) const
{
	return (*this) + (-1)*polynomial;
}

rg_Polynomial& rg_Polynomial::operator =(const rg_Polynomial &polynomial)
{
    if (this != &polynomial) {
        setDegree(polynomial.degree);

        for (rg_INT i = 0; i < degree + 1; i++) {
            coefficient[i] = polynomial.coefficient[i];
        }
    }

	return *(this);
}

rg_FLAG rg_Polynomial::operator ==(const rg_Polynomial &operand)
{
	rg_FLAG isEqual = rg_TRUE;

	for(rg_INDEX i = 0;i <= degree;i++)
	{
		if(coefficient[ i ] != operand.coefficient[ i ])
			isEqual = rg_FALSE;
	}

	if((degree == operand.degree) && isEqual)
		return rg_TRUE;
	else
		return rg_FALSE;
}

rg_REAL& rg_Polynomial::operator [](const rg_DEGREE index)
{
	return coefficient[ index ];
}

rg_Polynomial rg_Polynomial::power(const rg_INT &k) const
{
	if( k == 0 )
	{
		rg_REAL *var = new rg_REAL[1];
		var[0]      = 1.0;
		return rg_Polynomial(0, var);
	}

	else if( k == 1)
	{
		return *this;
	}

	else 
	{
		rg_Polynomial result((*this)*(*this));
		for(rg_INT i = 1; i <= k - 2; i++)
		{
			result = result * (*this);
		}
		return result;
	}
}

rg_ComplexNumber* rg_Polynomial::solve()
{
	static dcomplex  p[MAXCOEFF],// coefficient vector of polynomial   
				     pred[MAXCOEFF],// coeff. vector of deflated polynom. 
				     root[MAXCOEFF];// vector of determined roots         
	static unsigned rg_INT flag = rg_FALSE;
                     // flag = rg_TRUE  => complex coefficients      
                     // flag = rg_FALSE => real coefficients         
	unsigned char null(dcomplex *p, dcomplex *pred, rg_INT *n, dcomplex *root,
                   rg_REAL *maxerr, unsigned char flag);

	rg_ComplexNumber *rootValue = new rg_ComplexNumber[degree];
	
	rg_REAL maxerr;
	rg_INT i = 0;
	for(i = 0; i <= degree; i++)
		p[i].r = coefficient[i];

	rg_INT dg = degree;
	null(p, pred, &dg, root, &maxerr, (unsigned char)flag);
	
	for(i = 0; i < degree; i++)
	{
		rootValue[i].setRealNumber(root[i].r);
		rootValue[i].setImaginaryNumber(root[i].i);
	}

	return rootValue;
}

rg_ComplexNumber* rg_Polynomial::solve(const rg_REAL& RHS)
{
	static dcomplex  p[MAXCOEFF],// coefficient vector of polynomial   
				     pred[MAXCOEFF],// coeff. vector of deflated polynom. 
				     root[MAXCOEFF];// vector of determined roots         
	static unsigned rg_INT flag = rg_FALSE;
                     // flag = rg_TRUE  => complex coefficients      
                     // flag = rg_FALSE => real coefficients         
	unsigned char null(dcomplex *p, dcomplex *pred, rg_INT *n, dcomplex *root,
                   rg_REAL *maxerr, unsigned char flag);

	rg_ComplexNumber *rootValue = new rg_ComplexNumber[degree];
	
	rg_REAL maxerr;
	rg_INT i = 0;
	for(i = 0; i <= degree; i++)
		p[i].r = coefficient[i];

    p[0].r-=RHS;
	rg_INT dg = degree;
	null(p, pred, &dg, root, &maxerr, (unsigned char)flag);
	
	for(i = 0; i < degree; i++)
	{
		rootValue[i].setRealNumber(root[i].r);
		rootValue[i].setImaginaryNumber(root[i].i);
	}

	return rootValue;
}
// inserted in histrory by Hyung-Joo Lee
// function to find the first derivative of the polynomial
// The code for makeDerivative() start here
rg_Polynomial rg_Polynomial::makeDerivative() const
{
	rg_Polynomial derivative;
	
	if( degree == 0 )
	{
		derivative.setDegree( 0 );
		derivative.coefficient = new rg_REAL[1];
		derivative.setCoefficient( 0, 0.0 );

		return derivative;
	}
	
	derivative.setDegree( degree - 1 );

	rg_REAL* temp = new rg_REAL[degree];
	for( rg_INT i = 0; i < degree; i++ )
	{
		temp[i] = (i + 1) * getCoefficient(i + 1);
	}

	derivative.setCoefficient( temp );
	delete[] temp;

	return derivative;
}
// The code for makeDerivative() end here

// inserted in histrory by Hyung-Joo Lee
// function to find the first derivative of the polynomial
// The code for setUpdateDegree() start here
void rg_Polynomial::updateDegree()
{
	if( degree > 0 )
	{
		while( rg_ZERO( getCoefficient(degree), 0.000001 ) )
		{
			degree = degree - 1;
			if( degree == 0 )
			{
				break;
			}
		}
		rg_REAL* tempCoeff = new rg_REAL[degree+1];
		for( rg_INT i = 0; i <= degree; i++ )
		{
			tempCoeff[i] = coefficient[i];
		}
		setPolynomial( degree, tempCoeff );
		delete[] tempCoeff;
	}
}
// The code for setUpdateDegree() end here

void rg_Polynomial::displayCoefficient()
{
	for(rg_INT i = 0;i <= degree;i++);
		//cout << 'a' << i << '=' << coefficient[ i ] << " ";
}


