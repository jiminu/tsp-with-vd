/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_RPolynomialCurve2D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_RPolynomialCurve2D 
//           which define rational polynomialCurve rg_Curve whose dimesion is 2 and its property. 
//                          
//	  CLASS NAME  : rg_RPolynomialCurve2D
//
//    BASE CLASS  : rg_PolynomialCurve2D
//      
//    AUTHOR      : Taeboom Jang
//    START DATE  : 2000.3.18    
//
//    HISTORY     :
//
//           Copyright 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _RG_RPOLYNOMIALCURVE3D_H
#define _RG_RPOLYNOMIALCURVE3D_H

#include "rg_Const.h"
#include "rg_PolynomialCurve2D.h"
#include "rg_Point2D.h"

class rg_RPolynomialCurve2D : public rg_PolynomialCurve2D
{
private:
    rg_REAL* coefficientsOfDenominator;
public:
	rg_RPolynomialCurve2D();
	rg_RPolynomialCurve2D(const rg_DEGREE &tDegree);
	rg_RPolynomialCurve2D(const rg_DEGREE &tDegree, const rg_Point2D* tCoefficients);
	rg_RPolynomialCurve2D(const rg_PolynomialCurve2D &tCurve);
  	rg_RPolynomialCurve2D(const rg_RPolynomialCurve2D &temp);

	~rg_RPolynomialCurve2D();

	//Operations
    void removeAll();
	rg_Point2D    evaluatePt(const rg_PARAMETER& t) const;
    rg_RPolynomialCurve2D getDerivative() const;


	//Access elements
    rg_REAL     getCoefficientOfDenominator(const rg_INT &index) const ;
    rg_REAL*    getCoefficientsOfDenominator();
	
	void     setCoefficientOfDenominator(const rg_INDEX& i, const rg_REAL& tCoefficient);
	void     setCoefficientsOfDenominator(const rg_REAL* tCoefficients);
	void     setCurve(const rg_RPolynomialCurve2D& temp);
	rg_RPolynomialCurve2D& operator=(const rg_RPolynomialCurve2D& temp);
};

#endif


