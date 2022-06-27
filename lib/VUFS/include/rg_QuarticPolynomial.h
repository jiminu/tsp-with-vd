/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_QuarticPolynomial.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_QuarticPolynomial 
//           which define quartic polynomial and its property. 
//                          
//	  CLASS NAME  : rg_QuarticPolynomial
//
//    BASE CLASS  : rg_Polynomial
//      
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 18 Jul. 1996    
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _RG_QUARTICPOLYNOMIAL_H
#define _RG_QUARTICPOLYNOMIAL_H

#include "rg_Const.h"
#include "rg_ComplexNumber.h"
#include "rg_Polynomial.h"

// This is a class about a quartic polynomial.
// a0 * x^0 + a1 * x^1 + a2 * x^2 + a3 * x^3 + a4 * x^4, 
// member data : degree = 4, 
//               coefficient[0] = a0,
//               coefficient[1] = a1,
//               coefficient[2] = a2,
//               coefficient[3] = a3, 
//           and coefficient[4] = a4.

class rg_QuarticPolynomial : public rg_Polynomial
{
public:
	rg_QuarticPolynomial();
	rg_QuarticPolynomial(rg_REAL *coeff);
	rg_QuarticPolynomial(const rg_QuarticPolynomial &eq);
	~rg_QuarticPolynomial();

	void           setQuarticPolynomial(rg_REAL *coeff);
	rg_ComplexNumber *solve();

	rg_ComplexNumber *solve_OLD();

private:
	rg_ComplexNumber *solveQuadraticEq(rg_ComplexNumber *coeff2);
};

#endif  // rg_QuarticPolynomial class


