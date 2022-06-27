/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_CubicPolynomial.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_CubicPolynomial 
//           which define cubic polynomial and its property. 
//                          
//	  CLASS NAME  : rg_CubicPolynomial
//
//    BASE CLASS  : rg_Polynomial
//      
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 18 Jul. 1996    
//
//           Copyright ⓒ 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _RG_CUBICPOLYNOMIAL_H
#define _RG_CUBICPOLYNOMIAL_H

#include "rg_Const.h"
#include "rg_ComplexNumber.h"
#include "rg_Polynomial.h"

// This is a class about a cubic polynomial.
// a0 * x^0 + a1 * x^1 + a2 * x^2 + a3 * x^3, 
// member data : degree = 3, 
//               coefficient[0] = a0,
//               coefficient[1] = a1,
//               coefficient[2] = a2,
//               coefficient[3] = a3.

// rg_Polynomial class로부터 상속을 받는데, 
// 그 class의 member data가 protected로 되어 있다고 가정한다. 
// 추후에 생각해 볼 여지가 있는 부분이다. 

class rg_CubicPolynomial : public rg_Polynomial
{
public:
	rg_CubicPolynomial();
	rg_CubicPolynomial(rg_REAL *coeff);
	rg_CubicPolynomial(const rg_CubicPolynomial &eq);
	~rg_CubicPolynomial();

	void           setPolynomial(rg_REAL *coeff);
	rg_ComplexNumber *solve();

};

#endif  // rg_CubicPolynomial class


