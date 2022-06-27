/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_LinearPolynomial.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_LinearPolynomial 
//           which represents a linear polynomial and its property. 
//                          
//	  CLASS NAME  : rg_LinearPolynomial
//
//    BASE CLASS  : rg_Polynomial
//      
//    AUTHOR      : Ryu, Joonghyun
//    START DATE  : Jan. 7. 2010   
//
//           Copyright ¨Ï 2010 by VDRC at Hanyang University
//
/////////////////////////////////////////////////////////////////////


#ifndef _LINEAR_POLYNOMIAL_H_
#define _LINEAR_POLYNOMIAL_H_

#include "rg_Const.h"
#include "rg_ComplexNumber.h"
#include "rg_Polynomial.h"

class rg_LinearPolynomial : public rg_Polynomial
{
public:
	rg_LinearPolynomial();
	rg_LinearPolynomial(rg_REAL* coeff);
	rg_LinearPolynomial(const rg_LinearPolynomial& qPolynomial);
	~rg_LinearPolynomial();
	
	void              setPolynomial(rg_REAL *coeff);
	rg_ComplexNumber *solve();
};

#endif


