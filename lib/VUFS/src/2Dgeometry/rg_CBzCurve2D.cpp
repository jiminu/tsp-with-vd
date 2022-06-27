/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_CBzCurve2D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_CBzCurve2D 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jun 1997    
//
//           Copyright ⓒ 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include <math.h>
#include "rg_CBzCurve2D.h"
#include "rg_RelativeOp.h"
#include "rg_Point2D.h"
#include "rg_Line.h"
#include "rg_IntersectFunc.h"
#include "rg_MathFunc.h"
#include "rg_RQBzCurve2D.h"

rg_CBzCurve2D::rg_CBzCurve2D() :rg_BzCurve2D(rg_CUBIC)
{
}

rg_CBzCurve2D::rg_CBzCurve2D(const rg_Point2D *ctrlpt) :rg_BzCurve2D(rg_CUBIC, ctrlpt)
{
}

rg_CBzCurve2D::rg_CBzCurve2D(const rg_CBzCurve2D &curve) :rg_BzCurve2D(curve)
{
	degree = rg_CUBIC;
}

rg_CBzCurve2D::rg_CBzCurve2D(const rg_Point2D &p0, const rg_Point2D &p1, 
				   const rg_Point2D &p2, const rg_Point2D &p3)
{
	rg_Point2D *pt = new rg_Point2D[degree+1];
	pt[0] = p0;
	pt[1] = p1;
	pt[2] = p2;
	pt[3] = p3;

	degree = rg_CUBIC;
	setCtrlPts(pt);
}

rg_CBzCurve2D::~rg_CBzCurve2D()
{
}

rg_INT rg_CBzCurve2D::typeOfCurve() const
{
	rg_QBzCurve2D QCurve;
	QCurve.setCurve(makeDerivative());
	if(QCurve.isOriginExteriorPoint())
	{
		if(QCurve.ctrlPtsMonotoneInPolarAngle()) return 1;
		else return 2;
	}
	else 
	{
		if(QCurve.isOriginOnCurve()) return 5;
		else if(QCurve.isNotConjugateTangentVector()) return 3; 
		else if(isNotSelfIntersection()) return 41;
		else return 42;
	}
}

rg_INT rg_CBzCurve2D::isNotSelfIntersection() const
{
	rg_Point2D X(rg_IntersectFunc::intersectLineVsLine(rg_Line<rg_Point2D>(ctrlPts[0], ctrlPts[1]), 
		                          rg_Line<rg_Point2D>(ctrlPts[2], ctrlPts[3]) ));

	rg_Point2D a1(ctrlPts[1]-ctrlPts[0]);
	rg_Point2D a2(X           -ctrlPts[0]);
	rg_REAL alpha = a1.getX()/a2.getX();

	rg_Point2D b1(ctrlPts[3]-ctrlPts[2]);
	rg_Point2D b2(ctrlPts[3]-X);
	rg_REAL beta = b1.getX()/b2.getX();

	if( rg_LE(alpha, 1.0) || rg_LE(beta, 1.0) ) return 1;
	else 
	{
		if(rg_GT( (alpha-4.0/3.0)*(beta-4.0/3.0), 4.0/9.0) ) return 0;
		else return 1;
	}
}

rg_dList<rg_REAL> rg_CBzCurve2D::approximateRQBzCurves(const rg_REAL &t0, 
											  const rg_REAL &t1,
											  rg_dList<rg_RQBzCurve2D> &rqBzCurveList) const
{
	rg_QBzCurve2D Hodograph;
	Hodograph.setCurve(makeDerivative());
	rg_dList<rg_REAL> param;
	param.add(t0);
	rg_REAL *t;

	switch(typeOfCurve())
	{
		case 1:	// no inflection case
			t = rg_NULL;
			rqBzCurveList.add(makeOneRQBzCurve(t0, t1));
			break;

		case 2: // one inflection case
			t        = getInflectionParameter();
			param.add(t[0]);
			
			rqBzCurveList.add(makeOneRQBzCurve(t0, t0+(t1-t0)*t[0]));
			rqBzCurveList.add(makeOneRQBzCurve(t0+(t1-t0)*t[0], t1));
			break;

		case 3:  // two inflection case
			t           = getInflectionParameter();
			param.add(t[0]);
			param.add(t[1]);

			rqBzCurveList.add(makeOneRQBzCurve(t0, t0+(t1-t0)*t[0]));
			rqBzCurveList.add(makeOneRQBzCurve(t0+(t1-t0)*t[0], t0+(t1-t0)*t[1]));
			rqBzCurveList.add(makeOneRQBzCurve(t0+(t1-t0)*t[1], t1));
			break;		

		case 41: // conjugate tangent case 
			t           = Hodograph.getConjugateTangentParameters();
			param.add(t[0]);
			param.add(t[1]);

			rqBzCurveList.add(makeOneRQBzCurve(t0, t0+(t1-t0)*t[0]));
			rqBzCurveList.add(makeOneRQBzCurve(t0+(t1-t0)*t[0], t0+(t1-t0)*t[1]));
			rqBzCurveList.add(makeOneRQBzCurve(t0+(t1-t0)*t[1], t1));
			break;		

		case 42: // conjugate tangent case with selfintersection
			t           = Hodograph.getConjugateTangentParameters();
			param.add(t[0]);
			param.add(t[1]);

			rqBzCurveList.add(makeOneRQBzCurve(t0, t0+(t1-t0)*t[0]));
			rqBzCurveList.add(makeOneRQBzCurve(t0+(t1-t0)*t[0], t0+(t1-t0)*t[1]));
			rqBzCurveList.add(makeOneRQBzCurve(t0+(t1-t0)*t[1], t1));
			break;		

		case 5: // cusp case
			t    = new rg_REAL[1];
			t[0] = Hodograph.getOriginParameter();
			param.add(t[0]);
			
			rqBzCurveList.add(makeOneRQBzCurve(t0, t0+(t1-t0)*t[0]));
			rqBzCurveList.add(makeOneRQBzCurve(t0+(t1-t0)*t[0], t1));
            delete[] t;

			break;

		default:
			break;
	};

	param.add(t1);
	return param;
}

rg_RQBzCurve2D rg_CBzCurve2D::makeOneRQBzCurve(const rg_REAL &t0,
								     const rg_REAL &t1) const 
{
	rg_QBzCurve2D devCurve;
	devCurve.setCurve(makeDerivative());
	rg_Point2D          *pt          = new rg_Point2D[2];
	rg_Point2D          *pt_prime    = new rg_Point2D[2];

	pt[0]       = evaluatePt(t0);
	pt_prime[0] = pt[0] + devCurve.evaluatePt(t0);

	pt[1]       = evaluatePt(t1);
	pt_prime[1] = pt[1] + devCurve.evaluatePt(t1);

	rg_RQBzCurve2D rqcurve;
	rqcurve.makeRQBezier(pt[0], pt_prime[0], 
		                 pt[1], pt_prime[1],
						 rg_Point2D(evaluatePt((t0+t1)/2.0)));

	delete []pt;
	delete []pt_prime;
	
	return rqcurve;
}

rg_dList<rg_CBzCurve2D> rg_CBzCurve2D::subdivideTwoCurve(const rg_REAL &t)
{
	rg_dList<rg_CBzCurve2D>  curve;
	rg_Point2D  **b = deCasteljau(t);

	curve.add(rg_CBzCurve2D(b[0][0], b[1][0], b[2][0], b[3][0]));
	curve.add(rg_CBzCurve2D(b[3][0], b[2][1], b[1][2], b[0][3]));
	
	return curve;
}

void rg_CBzCurve2D::subdivideTwoCurve(const rg_REAL &t, 
								 rg_CBzCurve2D &curve1,
								 rg_CBzCurve2D &curve2)
{
	rg_Point2D  **b = deCasteljau(t);

	curve1.setCurve(rg_CBzCurve2D(b[0][0], b[1][0], b[2][0], b[3][0]));
	curve2.setCurve(rg_CBzCurve2D(b[3][0], b[2][1], b[1][2], b[0][3]));
	
}

rg_REAL *rg_CBzCurve2D::getInflectionParameter() const
{
   	rg_QBzCurve2D QCurve;
	QCurve.setCurve(makeDerivative());

    rg_Point2D *ctrlPt = QCurve.getCtrlPts();

	rg_Point2D a(   ctrlPt[0] - 2*ctrlPt[1] + ctrlPt[2] );
	rg_Point2D b(2*(ctrlPt[1] -   ctrlPt[0])            );
	rg_Point2D c(   ctrlPt[0]                           );

	// y(t) = k * x(t)를 만족하는 k 값을 찾는다. 
	rg_REAL  p =    b.getX()*b.getX() - 4*a.getX()*c.getX();
	rg_REAL  q = -2*b.getY()*b.getX() + 4*a.getY()*c.getX() + 4*a.getX()*c.getY();
	rg_REAL  r =    b.getY()*b.getY() - 4*a.getY()*c.getY();
	rg_REAL *k = rg_MathFunc::solveQuadraticEq(p, q, r);

	rg_REAL *t = new rg_REAL[2];
	t[0] = -(b.getY() - k[0]*b.getX())/(2*(a.getY()-k[0]*a.getX()));
	t[1] = -(b.getY() - k[1]*b.getX())/(2*(a.getY()-k[1]*a.getX()));

	// k값을 대입한 parameter의 값이 0에서 1사이의 값인지를 
	// 판단해서 return 한다. 
//	if(rg_BTORexclusive(0.0, -t[0], 1.0) && rg_BTORexclusive(0.0, -t[1], 1.0)) 
//	{
//		t[0] = -t[0];
//		t[1] = -t[1];
//	}
	if( rg_GT(t[0], t[1]) )
	{
		rg_REAL temp = t[0];
		t[0] = t[1];
		t[1] = temp;
	}

	if(rg_BTORexclusive(0.0, t[0], 1.0)) 
	{
		if(rg_BTORexclusive(0.0, t[1], 1.0)) return t;
		else return &t[0];
	}
	else 
	{
		if(rg_BTORexclusive(0.0, t[1], 1.0)) return &t[1];
		else return rg_NULL;
	}
}

void rg_CBzCurve2D::setCtrlPts(const rg_Point2D &b0, const rg_Point2D &b1, 
							const rg_Point2D &b2, const rg_Point2D &b3)
{
	ctrlPts[0] = b0;
	ctrlPts[1] = b1;
	ctrlPts[2] = b2;
	ctrlPts[3] = b3;
}

void rg_CBzCurve2D::setCtrlPts(const rg_Point2D *point)
{
	setDegree(rg_CUBIC);
	rg_BzCurve2D::setCtrlPts(point);
}


