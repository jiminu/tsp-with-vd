/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_QuadraticPolynomial.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_QuadraticPolynomial 
//           which represents a quadratic polynomial and its property. 
//                          
//	  CLASS NAME  : rg_QuadraticPolynomial
//
//    BASE CLASS  : rg_Polynomial
//      
//    AUTHOR      : Ryu, Joonghyun
//    START DATE  : Jan. 6. 2010   
//
//           Copyright ¨Ï 2010 by VDRC at Hanyang University
//
/////////////////////////////////////////////////////////////////////


#ifndef _QUADRATIC_POLYNOMIAL_H_
#define _QUADRATIC_POLYNOMIAL_H_

#include "rg_Const.h"
#include "rg_ComplexNumber.h"
#include "rg_Polynomial.h"

class rg_QuadraticPolynomial : public rg_Polynomial
{
public:
	rg_QuadraticPolynomial();
	rg_QuadraticPolynomial(rg_REAL* coeff);
	rg_QuadraticPolynomial(const rg_QuadraticPolynomial& qPolynomial);
	~rg_QuadraticPolynomial();

	void              setPolynomial(rg_REAL *coeff);
	rg_ComplexNumber *solve();
};

#endif


