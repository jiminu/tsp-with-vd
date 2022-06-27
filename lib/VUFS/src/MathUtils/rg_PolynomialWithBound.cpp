#include "rg_PolynomialWithBound.h"

#include <iostream>


// Constructors & Destructor

rg_PolynomialWithBound::rg_PolynomialWithBound() : rg_Polynomial()
{
	lowerBound = lowerInf;
	upperBound = upperInf;
}

rg_PolynomialWithBound::rg_PolynomialWithBound(const rg_DEGREE& dg) : rg_Polynomial(dg)
{
	lowerBound = lowerInf;
	upperBound = upperInf;
}

rg_PolynomialWithBound::rg_PolynomialWithBound(const rg_DEGREE &dg, rg_REAL *coeff) 
 : rg_Polynomial(dg, coeff)
{
	lowerBound = lowerInf;
	upperBound = upperInf;
}

rg_PolynomialWithBound::rg_PolynomialWithBound(const rg_REAL lower, const rg_REAL upper)
 : rg_Polynomial()
{
	lowerBound = lower;
	upperBound = upper;
}

rg_PolynomialWithBound::rg_PolynomialWithBound(const rg_DEGREE &dg, rg_REAL *coeff, const rg_REAL lower, const rg_REAL upper)
 : rg_Polynomial(dg, coeff)
{
	lowerBound = lower;
	upperBound = upper;
}

rg_PolynomialWithBound::rg_PolynomialWithBound(const rg_PolynomialWithBound& SourceObj)
 : rg_Polynomial(SourceObj.getDegree(), SourceObj.getCoefficient())
{
	lowerBound = SourceObj.lowerBound;
	upperBound = SourceObj.upperBound;
}

rg_PolynomialWithBound::rg_PolynomialWithBound(const rg_Polynomial& SourceObj) : rg_Polynomial( SourceObj )
{
}

rg_PolynomialWithBound::~rg_PolynomialWithBound()
{
}


// Member functions

void rg_PolynomialWithBound::setBound(const rg_REAL lower, const rg_REAL upper)
{
	lowerBound = lower;
	upperBound = upper;
}

rg_REAL* rg_PolynomialWithBound::getBound() const
{
	rg_REAL* bounds = new rg_REAL[ 2 ];

	bounds[ 0 ] = lowerBound;
	bounds[ 1 ] = upperBound;

	return bounds;
}

rg_REAL rg_PolynomialWithBound::getLowerBound() const
{
	return lowerBound;
}

rg_REAL rg_PolynomialWithBound::getUpperBound() const
{
	return upperBound;
}

rg_PolynomialWithBound& rg_PolynomialWithBound::operator =(const rg_PolynomialWithBound &SourceObj)
{
    if (this != &SourceObj) {
        rg_Polynomial::operator =(SourceObj);

        lowerBound = SourceObj.lowerBound;
        upperBound = SourceObj.upperBound;
    }

	return *(this);
}

/*
rg_PolynomialWithBound rg_PolynomialWithBound::operator  *(const rg_PolynomialWithBound &operand  ) const
{
		
	(*this) * operand;	
}


rg_PolynomialWithBound rg_PolynomialWithBound::operator  /(const rg_PolynomialWithBound &operand  ) const
{
}

rg_PolynomialWithBound rg_PolynomialWithBound::operator  +(const rg_PolynomialWithBound &operand  ) const
{
}

rg_PolynomialWithBound rg_PolynomialWithBound::operator  -(const rg_PolynomialWithBound &operand  ) const
{
}

rg_FLAG                rg_PolynomialWithBound::operator ==(const rg_PolynomialWithBound &operand  )
{
}
*/

/*
bool rg_PolynomialWithBound::solve(rg_ComplexNumber* rg_MathFunc::root)
{

	if()
	{
		return true;
	}
	else
	{
		return false;
	}
}*/

rg_PolynomialWithBound rg_PolynomialWithBound::makeDerivative() const
{
	rg_PolynomialWithBound derivative;
	
	if( degree == 0 )
	{
		derivative.setDegree( 0 );
		derivative.coefficient = new rg_REAL[1];
		derivative.setCoefficient( 0, 0.0 );
		derivative.setBound(lowerBound, upperBound);

		return derivative;
	}
	
	derivative.setDegree( degree - 1 );
	derivative.setBound(lowerBound, upperBound);

	rg_REAL* temp = new rg_REAL[degree];
	for( rg_INT i = 0; i < degree; i++ )
	{
		temp[i] = (i + 1) * getCoefficient(i + 1);
	}

	derivative.setCoefficient( temp );
	delete[] temp;

	return derivative;
}

void rg_PolynomialWithBound::display()
{
	for(rg_INT i = 0;i <= degree;i++);
		//cout << 'a' << i << '=' << coefficient[ i ] << " ";
	//cout << '[' << lowerBound << ", " << upperBound << ']';
}


