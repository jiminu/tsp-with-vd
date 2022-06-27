/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_PolynomialWithBound.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_PolynomialWithBound 
//           which define the polynomial with bound and its property. 
//                          
//	  CLASS NAME  : rg_PolynomialWithBound
//
//    BASE CLASS  : rg_Polynomial
//      
//    AUTHOR      : Ryu, JungHyun
//    START DATE  : 30 Aug. 1998    
//
//    Copyright ¨Ï 1998 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _RG_POLYNOMIALWITHBOUND_H
#define _RG_POLYNOMIALWITHBOUND_H

#include "rg_Polynomial.h"
#include "rg_Const.h"

const rg_REAL lowerInf = -10000;
const rg_REAL upperInf =  10000;

class rg_PolynomialWithBound : public rg_Polynomial
{
private:
	rg_REAL upperBound;
	rg_REAL lowerBound;

public:
	rg_PolynomialWithBound();
	rg_PolynomialWithBound(const rg_DEGREE& dg);
	rg_PolynomialWithBound(const rg_DEGREE &dg, rg_REAL *coeff);
	rg_PolynomialWithBound(const rg_REAL lower, const rg_REAL upper);
	rg_PolynomialWithBound(const rg_DEGREE &dg, rg_REAL *coeff, 
		                const rg_REAL lower, const rg_REAL upper);
	rg_PolynomialWithBound(const rg_PolynomialWithBound& SourceObj);
	rg_PolynomialWithBound(const rg_Polynomial& SourceObj);
	~rg_PolynomialWithBound();

	rg_PolynomialWithBound& operator  =(const rg_PolynomialWithBound &SourceObj);
    /*rg_PolynomialWithBound  operator  *(const rg_PolynomialWithBound &operand  ) const;
    rg_PolynomialWithBound  operator  /(const rg_PolynomialWithBound &operand  ) const;
	rg_PolynomialWithBound  operator  +(const rg_PolynomialWithBound &operand  ) const;
	rg_PolynomialWithBound  operator  -(const rg_PolynomialWithBound &operand  ) const;
	rg_FLAG                 operator ==(const rg_PolynomialWithBound &operand  );*/

public:
	void setBound(const rg_REAL lower, const rg_REAL upper);

	rg_REAL* getBound() const;
	rg_REAL  getLowerBound() const;
	rg_REAL  getUpperBound() const;

	//bool solve(rg_ComplexNumber* rg_MathFunc::root);

	rg_PolynomialWithBound makeDerivative() const;

	void display(); //test
};

#endif


