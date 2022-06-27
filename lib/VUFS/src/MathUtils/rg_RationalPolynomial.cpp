/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_RationalPolynomial.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_RationalPolynomial 
//           which define the polynomial with bound and its property. 
//                          
//	  CLASS NAME  : rg_RationalPolynomial
//
//    BASE CLASS  : None
//      
//    AUTHOR      : Ryu, JungHyun
//    START DATE  : 7 Sep. 1998    
//
//    Copyright ¨Ï 1998 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////


#include <iostream>

#include "rg_RationalPolynomial.h"
#include <math.h>
#include "rg_MathFunc.h"
#include "rg_RelativeOp.h"
#include "PolynomialSolver.h"


// Constructors & Destructor

rg_RationalPolynomial::rg_RationalPolynomial()
{
}

rg_RationalPolynomial::rg_RationalPolynomial(const rg_DEGREE& numerDegree,  const rg_DEGREE& denomiDegree)
{
	numerator.setDegree(numerDegree);
	denominator.setDegree(denomiDegree);
}

rg_RationalPolynomial::rg_RationalPolynomial(const rg_DEGREE& numerDegree,  const rg_REAL* numerCoeff,
	                                   const rg_DEGREE& denomiDegree, const rg_REAL* denomiCoeff)
{
	numerator.setPolynomial(numerDegree, numerCoeff);
	denominator.setPolynomial(denomiDegree, denomiCoeff);
}

rg_RationalPolynomial::rg_RationalPolynomial(const rg_Polynomial numeratorP, const rg_Polynomial denominatorP)
{
	numerator = numeratorP;
	denominator = denominatorP;
}

rg_RationalPolynomial::rg_RationalPolynomial(const rg_RationalPolynomial& rp)
{
	numerator = rp.numerator;
	denominator = rp.denominator;
}

rg_RationalPolynomial::~rg_RationalPolynomial()
{
}


// member functions

void rg_RationalPolynomial::setRationalPolynomial(const rg_DEGREE& numerDegree,  const rg_REAL* numerCoeff,
	                                           const rg_DEGREE& denomiDegree, const rg_REAL* denomiCoeff)
{
	numerator.setPolynomial(numerDegree, numerCoeff);
	denominator.setPolynomial(denomiDegree, denomiCoeff);
}

void rg_RationalPolynomial::setRationalPolynomial(const rg_Polynomial& numer, const rg_Polynomial& denomi)
{
	setNumerator(numer);
	setDenominator(denomi);
}

void rg_RationalPolynomial::setNumerator(const rg_Polynomial& numer)
{
	numerator = numer;
}

void rg_RationalPolynomial::setDenominator(const rg_Polynomial& denomi)
{
	denominator = denomi;
}

void rg_RationalPolynomial::setNumeraor(const rg_DEGREE& numerDegree,  const rg_REAL* numerCoeff)
{
	numerator.setPolynomial(numerDegree, numerCoeff);
}

void rg_RationalPolynomial::setDenominator(const rg_DEGREE& denomiDegree, const rg_REAL* denomiCoeff)
{
	denominator.setPolynomial(denomiDegree, denomiCoeff);
}

void rg_RationalPolynomial::setDegree(const rg_DEGREE& numerDegree,  const rg_DEGREE& denomiDegree)
{
	numerator.setDegree(numerDegree);
	denominator.setDegree(denomiDegree);
}

void rg_RationalPolynomial::setNumerDegree(const rg_DEGREE& numerDegree)
{
	numerator.setDegree(numerDegree);
}

void rg_RationalPolynomial::setDenomiDegree(const rg_DEGREE& denomiDegree)
{
	denominator.setDegree(denomiDegree);
}

void rg_RationalPolynomial::setCoefficient(const rg_REAL* numerCoeff, const rg_REAL* denomiCoeff)
{
	numerator.setCoefficient(numerCoeff);
	denominator.setCoefficient(denomiCoeff);
}

void rg_RationalPolynomial::setNumerCoefficient(const rg_INDEX& i, const rg_REAL& numerCoeff,
	                                         const rg_INDEX& j, const rg_REAL& denomiCoeff)
{
	numerator.setCoefficient(i, numerCoeff);
	denominator.setCoefficient(j, denomiCoeff);
}

void rg_RationalPolynomial::setNumerCoefficient(const rg_REAL* numerCoeff)
{
	numerator.setCoefficient(numerCoeff);
}

void rg_RationalPolynomial::setNumerCoefficient(const rg_INDEX& i, const rg_REAL& numerCoeff)
{
	numerator.setCoefficient(i, numerCoeff);
}

void rg_RationalPolynomial::setDenomiCoefficient(const rg_REAL* denomiCoeff)
{
	denominator.setCoefficient(denomiCoeff);
}

void rg_RationalPolynomial::setDenomiCoefficient(const rg_INDEX& j, const rg_REAL& denomiCoeff)
{
	denominator.setCoefficient(j, denomiCoeff);
}

rg_Polynomial rg_RationalPolynomial::getNumerator() const
{
	return numerator;
}

rg_Polynomial rg_RationalPolynomial::getDenominator() const
{
	return denominator;
}

rg_DEGREE* rg_RationalPolynomial::getDegree() const
{
	rg_DEGREE* degree = new rg_DEGREE [ 2 ];
	degree[ 0 ] = numerator.getDegree();
	degree[ 1 ] = denominator.getDegree();

	return degree;
}

rg_DEGREE rg_RationalPolynomial::getNumerDegree() const
{
	return numerator.getDegree();
}

rg_DEGREE rg_RationalPolynomial::getDenomiDegree() const
{
	return denominator.getDegree();
}

rg_REAL**  rg_RationalPolynomial::getCoefficient() const
{
	rg_REAL** coefficient = new rg_REAL* [ 2 ];

	coefficient[ 0 ] = numerator.getCoefficient();
	coefficient[ 1 ] = denominator.getCoefficient();

	return coefficient;
}

rg_REAL* rg_RationalPolynomial::getCoefficient(const rg_INDEX& i, const rg_INDEX& j) const
{
	rg_REAL* coefficient = new rg_REAL[ 2 ];

	coefficient[ 0 ] = numerator.getCoefficient( i );
	coefficient[ 1 ] = denominator.getCoefficient( j );

	return coefficient;
}

rg_REAL* rg_RationalPolynomial::getNumerCoefficient() const
{
	return numerator.getCoefficient();
}

rg_REAL rg_RationalPolynomial::getNumerCoefficient(const rg_INDEX& i) const
{
	return numerator.getCoefficient( i );
}

rg_REAL* rg_RationalPolynomial::getDenomiCoefficient() const
{
	return denominator.getCoefficient();
}

rg_REAL rg_RationalPolynomial::getDenomiCoefficient(const rg_INDEX& j) const
{
	return denominator.getCoefficient( j );
}


rg_REAL rg_RationalPolynomial::evaluatePolynomial(const rg_REAL& valueForVariable) const
{
	return numerator.evaluatePolynomial(valueForVariable) / denominator.evaluatePolynomial(valueForVariable);
}

rg_RationalPolynomial rg_RationalPolynomial::operator *( const rg_REAL &n ) const
{
	rg_Polynomial numeratorP = numerator * n;

	rg_RationalPolynomial result(numeratorP, denominator);

	return result;
}

rg_RationalPolynomial rg_RationalPolynomial::operator /(const rg_REAL &n) const
{
	rg_Polynomial numeratorP = numerator / n;

	rg_RationalPolynomial result(numeratorP, denominator);

	return result;
}


rg_RationalPolynomial operator *( const rg_REAL &n, const rg_RationalPolynomial &operand )
{
	return operand * n;
}

rg_RationalPolynomial rg_RationalPolynomial::operator *(const rg_RationalPolynomial &operand) const
{
	rg_Polynomial numeratorP = numerator * operand.numerator;
	rg_Polynomial denominatorP = denominator * operand.denominator;
	numeratorP.reconstruct();
	denominatorP.reconstruct();

	rg_RationalPolynomial result(numeratorP, denominatorP);

	return result;
}


rg_RationalPolynomial rg_RationalPolynomial::operator +(const rg_RationalPolynomial &operand) const
{
	rg_Polynomial numeratorP = (numerator * operand.denominator) + (denominator * operand.numerator);
	rg_Polynomial denominatorP = denominator * operand.denominator;
	numeratorP.reconstruct();
	denominatorP.reconstruct();

	rg_RationalPolynomial result(numeratorP, denominatorP);

	return result;
}

rg_RationalPolynomial rg_RationalPolynomial::operator -(const rg_RationalPolynomial &operand) const
{
	return (*this) + (-1)*operand;
}

rg_RationalPolynomial& rg_RationalPolynomial::operator =(const rg_RationalPolynomial &operand)
{
	(*this).numerator = operand.numerator;
	(*this).denominator = operand.denominator;

	return *(this);
}

rg_FLAG rg_RationalPolynomial::operator ==(const rg_RationalPolynomial &operand)
{
	if((numerator == operand.numerator) && (denominator == operand.denominator))
		return rg_TRUE;
	else
		return rg_FALSE;
}

rg_ComplexNumber* rg_RationalPolynomial::solve()
{
    //////////////////////////////////////////////////////////////////////////
	// Solvign the numerator of rational polynomial
	//////////////////////////////////////////////////////////////////////////

	static dcomplex  p[MAXCOEFF],// coefficient vector of polynomial   
				     pred[MAXCOEFF],// coeff. vector of deflated polynom. 
				     root[MAXCOEFF];// vector of determined roots         
	static unsigned rg_INT flag = rg_FALSE;
                     // flag = rg_TRUE  => complex coefficients      
                     // flag = rg_FALSE => real coefficients         
	unsigned char null(dcomplex *p, dcomplex *pred, rg_INT *n, dcomplex *root,
                   rg_REAL *maxerr, unsigned char flag);

	rg_DEGREE degree = numerator.getDegree();

	rg_ComplexNumber *rootValue = new rg_ComplexNumber[degree];
	
	rg_REAL maxerr;
	rg_INT i = 0;
	for(i = 0; i <= degree; i++)
		p[i].r = numerator.getCoefficient( i );

	rg_INT dg = degree;
	null(p, pred, &dg, root, &maxerr, (unsigned char)flag);
	
	for(i = 0; i < degree; i++)
	{
		rootValue[i].setRealNumber(root[i].r);
		rootValue[i].setImaginaryNumber(root[i].i);
	}
    //////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////
	// Solvign the denominator of rational polynomial
	//////////////////////////////////////////////////////////////////////////

	static dcomplex  p2[MAXCOEFF],// coefficient vector of polynomial   
				     pred2[MAXCOEFF],// coeff. vector of deflated polynom. 
				     root2[MAXCOEFF];// vector of determined roots         
	static unsigned rg_INT flag2 = rg_FALSE;
                     // flag = rg_TRUE  => complex coefficients      
                     // flag = rg_FALSE => real coefficients         
	unsigned char null(dcomplex *p2, dcomplex *pred2, rg_INT *n2, dcomplex *root2,
                   rg_REAL *maxerr2, unsigned char flag2);

	rg_DEGREE degree2 = denominator.getDegree();

	rg_ComplexNumber *rootValue2 = new rg_ComplexNumber[degree2];
	
	rg_REAL maxerr2;
	for(i = 0; i <= degree2; i++)
		p2[i].r = denominator.getCoefficient( i );

	rg_INT dg2 = degree2;
	null(p2, pred2, &dg2, root2, &maxerr2, (unsigned char)flag2);
	
	for(i = 0; i < degree2; i++)
	{
		rootValue2[i].setRealNumber(root2[i].r);
		rootValue2[i].setImaginaryNumber(root2[i].i);
	}

	//////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////

	rg_ComplexNumber *solution = new rg_ComplexNumber[degree];
	rg_INT count = 0;
	bool IsZeroOfDenomi = false;

	for(i = 0;i < degree;i++)
	{
		for(rg_INT j = 0;j < degree2;j++)
		{
			if(rootValue[ i ] == rootValue2[ j ])
				IsZeroOfDenomi = true;
		}

		if(!IsZeroOfDenomi)
			solution[count++] = rootValue[ i ];
	}

	return solution;
}


rg_RationalPolynomial rg_RationalPolynomial::makeDerivative() const
{
	rg_RationalPolynomial result;
	rg_Polynomial newNumerator, newDenominator;

	newDenominator = denominator.power( 2 );
	newNumerator = (numerator.makeDerivative() * denominator) - (numerator * denominator.makeDerivative());
	newNumerator.reconstruct();
	newDenominator.reconstruct();

	result.setRationalPolynomial(newNumerator, newDenominator);

	return result;
}

void rg_RationalPolynomial::display()
{
	rg_DEGREE denoDg = denominator.getDegree();
	rg_DEGREE nuDg = numerator.getDegree();
	rg_DEGREE mxDrg = (nuDg < denoDg) ? denoDg : nuDg;

	numerator.displayCoefficient();
	//cout << '\n';
	for(rg_INT i = 0;i < mxDrg + 1;i++)
		//cout << "-----";
	//cout << '\n';
	denominator.displayCoefficient();
	//cout << '\n';
}


