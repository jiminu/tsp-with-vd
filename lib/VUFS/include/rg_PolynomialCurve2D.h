/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_PolynomialCurve2D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_PolynomialCurve2D 
//           which define PolynomialCurve rg_Curve whose dimesion is 2 and its property. 
//                          
//	  CLASS NAME  : rg_PolynomialCurve2D
//
//    BASE CLASS  : rg_Curve
//      
//    AUTHOR      : Taeboom Jang
//    START DATE  : 2000.3.18    
//
//    HISTORY     :
//
//           Copyright 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _RG_POLYNOMIALCURVE2D_H
#define _RG_POLYNOMIALCURVE2D_H

#include "rg_Const.h"
#include "rg_Curve.h"
#include "rg_Point2D.h"

class rg_PolynomialCurve2D : public rg_Curve
{
private:
	rg_DEGREE   degree;
	rg_Point2D* coefficients;
public:


	rg_PolynomialCurve2D();
	rg_PolynomialCurve2D(const rg_DEGREE &tDegree);
	rg_PolynomialCurve2D(const rg_DEGREE &tDegree, const rg_Point2D* tCoefficients);
	rg_PolynomialCurve2D(const rg_PolynomialCurve2D &temp);
	~rg_PolynomialCurve2D();

    virtual void removeAll();
	//Operations
	rg_Point2D    evaluatePt(const rg_PARAMETER& t) const;


	//Access elements
	rg_DEGREE   getDegree() const;
	rg_Point2D  getCoefficient(const rg_DEGREE& i) const;
	rg_Point2D* getCoefficients() const;
	
	void     setDegree(const rg_DEGREE& i);
	void     setCoefficient(const rg_INDEX& i, const rg_Point2D& point);
	void     setCoefficients(const rg_Point2D* point);
	void     setCurve(const rg_PolynomialCurve2D& curve);
	rg_PolynomialCurve2D& operator=(const rg_PolynomialCurve2D& curve);
};

#endif


