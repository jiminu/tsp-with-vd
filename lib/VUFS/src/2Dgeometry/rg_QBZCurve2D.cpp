/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_QBzCurve2D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_QBzCurve2D 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jun 1997    
//
//           Copyright ⓒ 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include <math.h>
#include "rg_QBzCurve2D.h"
#include "rg_IntersectFunc.h"
#include "rg_MathFunc.h"
#include "rg_Line.h"
#include "sortFunc.h"

rg_QBzCurve2D::rg_QBzCurve2D() :rg_BzCurve2D(rg_QUADRATIC)
{
}

rg_QBzCurve2D::rg_QBzCurve2D(const rg_Point2D &p1, const rg_Point2D &p2, const rg_Point2D &p3)
{
	rg_Point2D *pt  = new rg_Point2D[rg_QUADRATIC+1];
	pt[0] = p1;
	pt[1] = p2;
	pt[2] = p3;

	//rg_BzCurve2D::rg_BzCurve2D(rg_QUADRATIC, pt);
	degree = rg_QUADRATIC;
	setCtrlPts(pt);
}

rg_QBzCurve2D::rg_QBzCurve2D(const rg_Point3D &p1, const rg_Point3D &p2, const rg_Point3D &p3)
{
	rg_Point2D *pt  = new rg_Point2D[rg_QUADRATIC];
	pt[0].setPoint(p1.getX(), p1.getY());
	pt[1].setPoint(p2.getX(), p2.getY());
	pt[2].setPoint(p3.getX(), p3.getY());

	//rg_BzCurve2D::rg_BzCurve2D(rg_QUADRATIC, pt);
	degree = rg_QUADRATIC;
	setCtrlPts(pt);

}

rg_QBzCurve2D::rg_QBzCurve2D(const rg_Point2D *ctrlpt) :rg_BzCurve2D(rg_QUADRATIC, ctrlpt)
{
}

rg_QBzCurve2D::rg_QBzCurve2D(const rg_QBzCurve2D &curve) :rg_BzCurve2D(curve)
{
	degree = rg_QUADRATIC;
}

rg_QBzCurve2D::~rg_QBzCurve2D()
{
}

rg_INT rg_QBzCurve2D::isOriginExteriorPoint()
{
	rg_Point2D ZR(0.0, 0.0);
	if( (ctrlPts[1]-ctrlPts[0])*(ctrlPts[2]-ctrlPts[1]) < 0)
	{
		if(   (ctrlPts[1]-ctrlPts[0])*(ZR-ctrlPts[0]) < 0 
		   && (ctrlPts[2]-ctrlPts[1])*(ZR-ctrlPts[1]) < 0
		   && (ctrlPts[0]-ctrlPts[2])*(ZR-ctrlPts[2]) < 0 ) 
				return 0; //interior point
		else return 1;
	}
	else 
	{
		if(   (ctrlPts[1]-ctrlPts[0])*(ZR-ctrlPts[0]) > 0 
		   && (ctrlPts[2]-ctrlPts[1])*(ZR-ctrlPts[1]) > 0
		   && (ctrlPts[0]-ctrlPts[2])*(ZR-ctrlPts[2]) > 0 ) 
				return 0; //interior point
		else return 1;
	}
}

rg_INT rg_QBzCurve2D::ctrlPtsMonotoneInPolarAngle()
{
	if(ctrlPtsIncreasingInPolarAngle() || ctrlPtsDecreasingInPolarAngle())
		return 1;

	else 
		return 0;
}	

rg_INT rg_QBzCurve2D::ctrlPtsIncreasingInPolarAngle()
{
	if( ctrlPts[0]*ctrlPts[1] > 0)
	{
		if( ctrlPts[1]*ctrlPts[2] > 0) return 1;
		else return 0;
	}
	else return 0;
}

rg_INT rg_QBzCurve2D::ctrlPtsDecreasingInPolarAngle()
{
//	if(crossProduct(ctrlPts[0], ctrlPts[1]) < 0)
//	{
//		if(crossProduct(ctrlPts[1], ctrlPts[2]) < 0) return 1;
	if( ctrlPts[0]*ctrlPts[1] < 0)
	{
		if( ctrlPts[1]*ctrlPts[2] < 0) return 1;
		else return 0;
	}
	else return 0;
}

rg_REAL rg_QBzCurve2D::isOriginOnCurve()
{
	return parameterValueOfPt(rg_Point2D(0.0, 0.0));
}

rg_REAL rg_QBzCurve2D::getOriginParameter()
{
	return isOriginOnCurve();
}

rg_REAL rg_QBzCurve2D::parameterValueOfPt(const rg_Point2D &point)
{
	rg_Point2D a(   ctrlPts[0] - 2*ctrlPts[1] + ctrlPts[2] );
	rg_Point2D b(2*(ctrlPts[1] -   ctrlPts[0])               );
	rg_Point2D c(   ctrlPts[0] -   point                       );

	rg_Point2D discriminant(  pow(b.getX(),2)- 4*a.getX()*c.getX(),
		                   pow(b.getY(),2)- 4*a.getY()*c.getY()  );
	if( rg_POS(discriminant.getX()) && rg_POS(discriminant.getY()) )
	{
		rg_REAL *tx = rg_MathFunc::solveQuadraticEq(a.getX(), b.getX(), c.getX());
        rg_REAL *ty = rg_MathFunc::solveQuadraticEq(a.getY(), b.getY(), c.getY());

		if( rg_BTORexclusive(0.0, tx[0], 1.0) )
		{
			if( rg_EQ(tx[0], ty[0]) || rg_EQ(tx[0], ty[1]) )
			{
				return tx[0];
			}
			
			else  return 0.0;
		}
		
		else if( rg_BTORexclusive(0.0, tx[1], 1.0) )
		{
			if( rg_EQ(tx[1], ty[0]) || rg_EQ(tx[1], ty[1]) )
			{
				return tx[1];
			}

			else  return 0.0;
		}

		else  return 0.0;
	}
	else return 0;
}

/*
rg_REAL rg_QBzCurve2D::parameterValueOfPt(const rg_Point2D &point)
{
	rg_Point2D a(   ctrlPts[0] - 2*ctrlPts[1] + ctrlPts[2] );
	rg_Point2D b(2*(ctrlPts[1] -   ctrlPts[0])               );
	rg_Point2D c(   ctrlPts[0] -   point                       );

	rg_Point2D discriminant(b*b - 4*a*c);
	if( rg_POS(discriminant.getX()) && rg_POS(discriminant.getY()) )
	{
		rg_Point2D *t = rg_MathFunc::solveQuadraticEq(a, b, c);
		if(    rg_BTORexclusive(0.0, t[0].getX(), 1.0) &&  rg_BTORexclusive(0.0, t[0].getY(), 1.0) 
			&& rg_EQ(t[0].getX(), t[0].getY()) ) 
			return t[0].getX();
		
		else if( rg_BTORexclusive(0.0, t[1].getX(), 1.0) &&  rg_BTORexclusive(0.0, t[1].getY(), 1.0) 
			  && rg_EQ(t[1].getX(), t[1].getY()) ) 
			return t[1].getX();

		else if( rg_BTORexclusive(0.0, t[0].getX(), 1.0) &&  rg_BTORexclusive(0.0, t[1].getY(), 1.0) 
			  && rg_EQ(t[0].getX(), t[1].getY()) ) 
			return t[0].getX();

		else if( rg_BTORexclusive(0.0, t[1].getX(), 1.0) &&  rg_BTORexclusive(0.0, t[0].getY(), 1.0) 
			  && rg_EQ(t[1].getX(), t[0].getY()) ) 
			return t[1].getX();

		else return 0;
	}
	else return 0;
}
*/

rg_REAL rg_QBzCurve2D::parameterValueOfPt(const rg_REAL &x, const rg_REAL &y)
{
	return parameterValueOfPt(rg_Point2D(x, y));
}

rg_INT rg_QBzCurve2D::isNotConjugateTangentVector()
{
	rg_Point2D ZR(0.0, 0.0);

	rg_Point2D *H_t = rg_IntersectFunc::intersectQBezierVsLine(*this, rg_Line<rg_Point2D>(ZR, ctrlPts[0]));
	rg_Point2D  H_otherIntersect;

	if( H_t[0] == ctrlPts[0] ) H_otherIntersect = H_t[1];
	else H_otherIntersect = H_t[0];

	rg_Point2D comparison(ctrlPts[0] - H_otherIntersect);
	if( rg_LT(comparison.magnitude(), ctrlPts[0].magnitude()) ) return 1;
	else return 0;
}

rg_dList<rg_REAL> rg_QBzCurve2D::getParameterOntheAxis()
{
	rg_Point2D a(   ctrlPts[0] - 2*ctrlPts[1] + ctrlPts[2] );
	rg_Point2D b(2*(ctrlPts[1] -   ctrlPts[0])               );
	rg_Point2D c(   ctrlPts[0]                                 );
	
	rg_REAL *tx = rg_MathFunc::solveQuadraticEq(a.getX(), b.getX(), c.getX());
	rg_REAL *ty = rg_MathFunc::solveQuadraticEq(a.getY(), b.getY(), c.getY());

	rg_REAL *tm = new rg_REAL[4];
	tm[0] = tx[0];
	tm[1] = tx[1];
	tm[2] = ty[0];
	tm[3] = ty[1];
	
	QuickSort(tm, 0, 3);	
	rg_dList <rg_REAL> t;
	for(rg_INT i = 0; i < 4; i++)
		if(rg_BTORexclusive(0.0, tm[i], 1.0)) t.add(tm[i]);

//	delete []tm;
	return t;
}

rg_REAL *rg_QBzCurve2D::getConjugateTangentParameters()
{
	rg_Point2D  ZR(0.0, 0.0);
	rg_Point2D  H_otherIntersect;
	rg_REAL  *param = new rg_REAL[2];

	rg_Point2D *t1 = rg_IntersectFunc::intersectQBezierVsLine(*this, rg_Line<rg_Point2D>(ZR, ctrlPts[2]));
	if( t1[0] == ctrlPts[2] ) param[0] = parameterValueOfPt(t1[1]);
	else param[0] = parameterValueOfPt(t1[0]);

	rg_Point2D *t0 = rg_IntersectFunc::intersectQBezierVsLine(*this, rg_Line<rg_Point2D>(ZR, ctrlPts[0]));
	if( t0[0] == ctrlPts[0] ) param[1] = parameterValueOfPt(t0[1]);
	else param[1] = parameterValueOfPt(t0[0]);	
	
	return param;
}

/*
// This function must be changed the name
rg_REAL *rg_QBzCurve2D::getInflectionParameter()
{
	rg_Point2D a(   ctrlPts[0] - 2*ctrlPts[1] + ctrlPts[2] );
	rg_Point2D b(2*(ctrlPts[1] -   ctrlPts[0])               );
	rg_Point2D c(   ctrlPts[0]                                 );

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
*/

void rg_QBzCurve2D::setCtrlPts(const rg_Point2D &p1, const rg_Point2D &p2, const rg_Point2D &p3)
{
	ctrlPts[0] = p1;
	ctrlPts[1] = p2;
	ctrlPts[2] = p3;
}

void rg_QBzCurve2D::setCtrlPts(const rg_Point2D *point)
{
	setDegree(rg_QUADRATIC);
	rg_BzCurve2D::setCtrlPts(point);
}


