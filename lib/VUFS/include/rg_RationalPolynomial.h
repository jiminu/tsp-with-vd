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

#ifndef _RG_RATIONALPOLYNOMIAL_H
#define _RG_RATIONALPOLYNOMIAL_H

#include "rg_Polynomial.h"
#include "rg_Const.h"

class rg_RationalPolynomial
{
private:
	rg_Polynomial numerator;
	rg_Polynomial denominator;

public:
	rg_RationalPolynomial();
	rg_RationalPolynomial(const rg_DEGREE& numerDegree,  const rg_DEGREE& denomiDegree);
	rg_RationalPolynomial(const rg_DEGREE& numerDegree,  const rg_REAL* numerCoeff,
		               const rg_DEGREE& denomiDegree, const rg_REAL* denomiCoeff);
	rg_RationalPolynomial(const rg_Polynomial numeratorP, const rg_Polynomial denominatorP);
	rg_RationalPolynomial(const rg_RationalPolynomial& rp);
	~rg_RationalPolynomial();

public:
	void setRationalPolynomial(const rg_DEGREE& numerDegree,  const rg_REAL* numerCoeff,
		                       const rg_DEGREE& denomiDegree, const rg_REAL* denomiCoeff);
    void setRationalPolynomial(const rg_Polynomial& numer, const rg_Polynomial& denomi);
	void setNumerator(const rg_Polynomial& numer);
	void setDenominator(const rg_Polynomial& denomi);
	void setNumeraor(const rg_DEGREE& numerDegree,  const rg_REAL* numerCoeff);
	void setDenominator(const rg_DEGREE& denomiDegree, const rg_REAL* denomiCoeff);

	void setDegree(const rg_DEGREE& numerDegree,  const rg_DEGREE& denomiDegree);
	void setNumerDegree(const rg_DEGREE& numerDegree);
	void setDenomiDegree(const rg_DEGREE& denomiDegree);

	void setCoefficient(const rg_REAL* numerCoeff, const rg_REAL* denomiCoeff);
	void setNumerCoefficient(const rg_INDEX& i, const rg_REAL& numerCoeff,
		                     const rg_INDEX& j, const rg_REAL& denomiCoeff);
	void setNumerCoefficient(const rg_REAL* numerCoeff);
	void setNumerCoefficient(const rg_INDEX& i, const rg_REAL& numerCoeff);
    void setDenomiCoefficient(const rg_REAL* denomiCoeff);
	void setDenomiCoefficient(const rg_INDEX& j, const rg_REAL& denomiCoeff);

	rg_Polynomial getNumerator() const;
	rg_Polynomial getDenominator() const;
	rg_DEGREE* getDegree() const;
	rg_DEGREE  getNumerDegree() const;
	rg_DEGREE  getDenomiDegree() const;

	rg_REAL**  getCoefficient() const;
    rg_REAL*   getCoefficient(const rg_INDEX& i, const rg_INDEX& j) const;
	rg_REAL*   getNumerCoefficient() const;
	rg_REAL    getNumerCoefficient(const rg_INDEX& i) const;
	rg_REAL*   getDenomiCoefficient() const;
	rg_REAL    getDenomiCoefficient(const rg_INDEX& j) const;

	rg_REAL    evaluatePolynomial(const rg_REAL& valueForVariable) const;

    rg_RationalPolynomial  operator *(const rg_RationalPolynomial &operand) const;
	rg_RationalPolynomial  operator +(const rg_RationalPolynomial &operand) const;
	rg_RationalPolynomial  operator -(const rg_RationalPolynomial &operand) const;
	rg_RationalPolynomial  operator *(const rg_REAL &n) const;
	rg_RationalPolynomial  operator /(const rg_REAL &n) const;
	rg_RationalPolynomial& operator =(const rg_RationalPolynomial &operand);
	rg_FLAG        operator ==(const rg_RationalPolynomial &operand);

	friend rg_RationalPolynomial operator *(const rg_REAL &n, const rg_RationalPolynomial &operand);

	rg_ComplexNumber* solve();


	rg_RationalPolynomial makeDerivative() const;

	void display();
};

#endif


