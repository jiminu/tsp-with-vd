#include "rg_IntsecImplicit4.h"
#include "rg_Polynomial.h"
#include "rg_RelativeOp.h"

// inserted by Jung-Hyun Ryu
#include "rg_IntersectFunc.h"

//#include <fstream>
#include <time.h>

rg_IntsecImplicit4::rg_IntsecImplicit4()
{

}

rg_IntsecImplicit4::~rg_IntsecImplicit4()
{

}

rg_dList<rg_Point2D> rg_IntsecImplicit4::intersectBzCurveVsBzCurve(const rg_BzCurve2D &curve_s, const rg_BzCurve2D &curve_t, rg_REAL &time)
{
	rg_dList<rg_Point2D> intersectPointList;
	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime = clock(); 
	//for(rg_INT i=0; i<20; i++)
	//{

		intersectPointList.removeAll();
		rg_ImplicitEquation  impEq   = implicitize(curve_t);

		//rg_ImplicitEquation  impEq    = rg_IntersectFunc::implicitizeNotUsingMathematica1(curve_t);//test
		//rg_ImplicitEquation  impEq    = rg_IntersectFunc::implicitizeNotUsingMathematica2(curve_t);//test

		rg_Polynomial*       cs      = curve_s.convertBzCurve2Polynomial();
		rg_Polynomial        poly    = impEq.substitute(cs[0], cs[1]);
		delete[] cs;
		rg_ComplexNumber*    param_s = poly.solve();

	// 	rg_dList<rg_Point2D> intersectPointList;
		for(rg_INT i = 0; i < poly.getDegree(); i++)
		{
			if(param_s[i].isPureRealNumber() && rg_BTOR(0.0, param_s[i].getRealNumber(), 1.0) )
			{
				rg_Point2D point_s = curve_s.evaluatePt(param_s[i].getRealNumber());
				rg_REAL  param_t = inversion(curve_t, point_s); 
				rg_Point2D point_t = curve_t.evaluatePt(param_t);
				if(rg_BTOR(0.0, param_t, 1.0))
				{
					intersectPointList.add(point_s);
				}
			}
		}

		delete[] param_s;
	//}

	/*
	/////////////////////////////////////////// test

	rg_ImplicitEquation  impEq1   = rg_IntersectFunc::implicitizeNotUsingMathematica1(curve_t);//test
	rg_ImplicitEquation  impEq2   = rg_IntersectFunc::implicitizeNotUsingMathematica2(curve_t);//test
		

	/////////////////////////////////////////// test
	ofstream fout("sylvester.dat", ios::app);


	for(rg_INT j=0; j < impEq1.getDegree()+1; j++)
	{
		for( rg_INDEX k=0; k < impEq1.getDegree()+1-j; k++)
		{
			fout.width(5);
			fout<<impEq1.getCoeff(j, k)<<"x^"<<j<<"y^"<<k<< "  ";
		}
		cout<<endl;
	}

	ofstream fout2("bezout.dat", ios::app);
	for(rg_INT l=0; l < impEq2.getDegree()+1; l++)
	{
		for(rg_INT m=0; m < impEq2.getDegree()+1-l; m++)
		{
			fout2.width(5);
			fout2<<impEq2.getCoeff(l, m)<<"x^"<<l<<"y^"<<m<< "  ";
		}
		cout<<endl;
	}
	/////////////////////////////////////////// test
*/
	clock_t EndTime = clock();
	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	
 	return intersectPointList;
}

rg_dList<rg_Point2D> rg_IntsecImplicit4::intersectBzCurveVsBzCurve( const rg_BzCurve2D &curve_s, 
														   const rg_BzCurve2D &curve_t, 
														   rg_Polynomial  &poly )
{
	rg_dList<rg_Point2D> intersectPointList;
	// This part is to find a seed
	// record the time to find the seed using implicitization approach

	for(rg_INT i=0; i<20; i++)
	{

		intersectPointList.removeAll();
		rg_ImplicitEquation  impEq   = implicitize(curve_t);
		rg_Polynomial*       cs      = curve_s.convertBzCurve2Polynomial();
		poly    = impEq.substitute(cs[0], cs[1]);
		delete[] cs;
		rg_ComplexNumber*    param_s = poly.solve();

	// 	rg_dList<rg_Point2D> intersectPointList;
		for(rg_INT i = 0; i < poly.getDegree(); i++)
		{
			if(param_s[i].isPureRealNumber() && rg_BTOR(0.0, param_s[i].getRealNumber(), 1.0) )
			{
				rg_Point2D point_s = curve_s.evaluatePt(param_s[i].getRealNumber());
				rg_REAL  param_t = inversion(curve_t, point_s); 
				rg_Point2D point_t = curve_t.evaluatePt(param_t);
				if(rg_BTOR(0.0, param_t, 1.0))
				{
					intersectPointList.add(point_s);
				}
			}
		}

		delete[] param_s;
	}

 	return intersectPointList;
}

  

rg_ImplicitEquation rg_IntsecImplicit4::implicitize(const rg_BzCurve2D &curve)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL x3 = curve.getCtrlPt(3).getX();
	rg_REAL x4 = curve.getCtrlPt(4).getX();
	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL y3 = curve.getCtrlPt(3).getY();
	rg_REAL y4 = curve.getCtrlPt(4).getY();

//  If
//		    
//	x(t) = a4 t^4 + a3 t^3 + a2 t^2 + a1 t + a0
//
//	y(t) = b4 t^4 + b3 t^3 + b2 t^2 + b1 t + b0
//  
//  then the resultant of p(x, t) and q(y, t) is like the following.

	rg_REAL a4 = x4 - 4*x3 +  6*x2 -  4*x1 +   x0;
	rg_REAL a3 =      4*x3 - 12*x2 + 12*x1 - 4*x0;
	rg_REAL a2 =              6*x2 - 12*x1 + 6*x0;
	rg_REAL a1 =                      4*x1 - 4*x0;
	rg_REAL a0 =                               x0;

	rg_REAL b4 = y4 - 4*y3 +  6*y2 -  4*y1 +   y0;
	rg_REAL b3 =      4*y3 - 12*y2 + 12*y1 - 4*y0;
	rg_REAL b2 =              6*y2 - 12*y1 + 6*y0;
	rg_REAL b1 =                      4*y1 - 4*y0;
	rg_REAL b0 =                               y0;

  //			|															                   |
  //			|	m0 x + m1 y + m2	n0 x + n1 y + n2	o0 x + o1 y + o2  p0 x + p1 y + p2 |
  //			|															                   |
  // R(P, Q) =	|	n0 x + n1 y + n2	q0 x + q1 y + q2	r0 x + r1 y + r2  s0 x + s1 y + s2 |
  //			|															                   |
  //			|	o0 x + o1 y + o2	r0 x + r1 y + r2	t0 x + t1 y + t2  u0 x + u1 y + u2 |
  //			|															                   |
  //            |   p0 x + p1 y + p2    s0 x + s1 y + s2    u0 x + u1 y + u2  v0 x + v1 y + v2 |
  //            |                                                                              |
  
	rg_REAL m0 =   0.0;
	rg_REAL m1 =   0.0;
	rg_REAL m2 =   a4*b3 - a3*b4;

  	rg_REAL n0 =   0.0;
	rg_REAL n1 =   0.0;
	rg_REAL n2 =   a4*b2 - a2*b4;

  	rg_REAL o0 =   0.0;
	rg_REAL o1 =   0.0;
	rg_REAL o2 =   a4*b1 - a1*b4;

  	rg_REAL p0 =   b4;
	rg_REAL p1 = - a4;
	rg_REAL p2 =   a4*b0 - a0*b4;

   	rg_REAL q0 =   0.0;
	rg_REAL q1 =   0.0;
	rg_REAL q2 =   a4*b1 + a3*b2 - a2*b3 - a1*b4;

   	rg_REAL r0 =   b4;
	rg_REAL r1 = - a4;
	rg_REAL r2 =   a4*b0 + a3*b1 - a1*b3 - a0*b4;

	rg_REAL s0 =   b3;
	rg_REAL s1 = - a3;
	rg_REAL s2 =   a3*b0 - a0*b3;

   	rg_REAL t0 =   b3;
	rg_REAL t1 = - a3;
	rg_REAL t2 =   a3*b0 + a2*b1 - a1*b2 - a0*b3;

	rg_REAL u0 =   b2;
	rg_REAL u1 = - a2;
	rg_REAL u2 =   a2*b0 - a0*b2;

	rg_REAL v0 =   b1;
	rg_REAL v1 = - a1;
	rg_REAL v2 =   a1*b0 - a0*b1;

// If the form of the implicitization is 
//    k0 
//  + k1  x       + k2      y
//  + k3  x^2     + k4  x   y   + k5      y^2
//  + k6  x^3     + k7  x^2 y   + k8  x   y^2  + k9      y^3
//  + k10 x^4     + k11 x^3 y   + k12 x^2 y^2  + k13 x   y^3 + k14 y^4

    rg_REAL *k = new rg_REAL[15];
    
    k[0]  = p2*p2*r2*r2 - 2*o2*p2*r2*s2 + o2*o2*s2*s2 - p2*p2*q2*t2 + 2*n2*p2*s2*t2 - 
  m2*s2*s2*t2 + 2*o2*p2*q2*u2 - 2*n2*p2*r2*u2 - 2*n2*o2*s2*u2 + 
  2*m2*r2*s2*u2 + n2*n2*u2*u2 - m2*q2*u2*u2 - o2*o2*q2*v2 + 2*n2*o2*r2*v2 - 
  m2*r2*r2*v2 - n2*n2*t2*v2 + m2*q2*t2*v2;
    
    k[1]  = 2*p2*p2*r0*r2 + 2*p0*p2*r2*r2 - 
  2*o2*p2*r2*s0 - 2*o2*p2*r0*s2 - 2*o2*p0*r2*s2 - 2*o0*p2*r2*s2 + 
  2*o2*o2*s0*s2 + 2*o0*o2*s2*s2 - p2*p2*q2*t0 + 2*n2*p2*s2*t0 - 
  m2*s2*s2*t0 - p2*p2*q0*t2 - 2*p0*p2*q2*t2 + 2*n2*p2*s0*t2 + 
  2*n2*p0*s2*t2 + 2*n0*p2*s2*t2 - 2*m2*s0*s2*t2 - m0*s2*s2*t2 + 
  2*o2*p2*q2*u0 - 2*n2*p2*r2*u0 - 2*n2*o2*s2*u0 + 2*m2*r2*s2*u0 + 
  2*o2*p2*q0*u2 + 2*o2*p0*q2*u2 + 2*o0*p2*q2*u2 - 2*n2*p2*r0*u2 - 
  2*n2*p0*r2*u2 - 2*n0*p2*r2*u2 - 2*n2*o2*s0*u2 + 2*m2*r2*s0*u2 - 
  2*n2*o0*s2*u2 - 2*n0*o2*s2*u2 + 2*m2*r0*s2*u2 + 2*m0*r2*s2*u2 + 
  2*n2*n2*u0*u2 - 2*m2*q2*u0*u2 + 2*n0*n2*u2*u2 - m2*q0*u2*u2 - 
  m0*q2*u2*u2 - o2*o2*q2*v0 + 2*n2*o2*r2*v0 - m2*r2*r2*v0 - 
  n2*n2*t2*v0 + m2*q2*t2*v0 - o2*o2*q0*v2 - 2*o0*o2*q2*v2 + 
  2*n2*o2*r0*v2 + 2*n2*o0*r2*v2 + 2*n0*o2*r2*v2 - 2*m2*r0*r2*v2 - 
  m0*r2*r2*v2 - n2*n2*t0*v2 + m2*q2*t0*v2 - 2*n0*n2*t2*v2 + 
  m2*q0*t2*v2 + m0*q2*t2*v2;
    
    k[2]  = 2*p2*p2*r1*r2 + 2*p1*p2*r2*r2 - 
  2*o2*p2*r2*s1 - 2*o2*p2*r1*s2 - 2*o2*p1*r2*s2 - 2*o1*p2*r2*s2 + 
  2*o2*o2*s1*s2 + 2*o1*o2*s2*s2 - p2*p2*q2*t1 + 2*n2*p2*s2*t1 - 
  m2*s2*s2*t1 - p2*p2*q1*t2 - 2*p1*p2*q2*t2 + 2*n2*p2*s1*t2 + 
  2*n2*p1*s2*t2 + 2*n1*p2*s2*t2 - 2*m2*s1*s2*t2 - m1*s2*s2*t2 + 
  2*o2*p2*q2*u1 - 2*n2*p2*r2*u1 - 2*n2*o2*s2*u1 + 2*m2*r2*s2*u1 + 
  2*o2*p2*q1*u2 + 2*o2*p1*q2*u2 + 2*o1*p2*q2*u2 - 2*n2*p2*r1*u2 - 
  2*n2*p1*r2*u2 - 2*n1*p2*r2*u2 - 2*n2*o2*s1*u2 + 2*m2*r2*s1*u2 - 
  2*n2*o1*s2*u2 - 2*n1*o2*s2*u2 + 2*m2*r1*s2*u2 + 2*m1*r2*s2*u2 + 
  2*n2*n2*u1*u2 - 2*m2*q2*u1*u2 + 2*n1*n2*u2*u2 - m2*q1*u2*u2 - 
  m1*q2*u2*u2 - o2*o2*q2*v1 + 2*n2*o2*r2*v1 - m2*r2*r2*v1 - 
  n2*n2*t2*v1 + m2*q2*t2*v1 - o2*o2*q1*v2 - 2*o1*o2*q2*v2 + 
  2*n2*o2*r1*v2 + 2*n2*o1*r2*v2 + 2*n1*o2*r2*v2 - 2*m2*r1*r2*v2 - 
  m1*r2*r2*v2 - n2*n2*t1*v2 + m2*q2*t1*v2 - 2*n1*n2*t2*v2 + 
  m2*q1*t2*v2 + m1*q2*t2*v2;
    
    k[3]  = p2*p2*r0*r0 + 4*p0*p2*r0*r2 + 
  p0*p0*r2*r2 - 2*o2*p2*r0*s0 - 2*o2*p0*r2*s0 - 2*o0*p2*r2*s0 + 
  o2*o2*s0*s0 - 2*o2*p0*r0*s2 - 2*o0*p2*r0*s2 - 2*o0*p0*r2*s2 + 
  4*o0*o2*s0*s2 + o0*o0*s2*s2 - p2*p2*q0*t0 - 2*p0*p2*q2*t0 + 
  2*n2*p2*s0*t0 + 2*n2*p0*s2*t0 + 2*n0*p2*s2*t0 - 
  2*m2*s0*s2*t0 - m0*s2*s2*t0 - 2*p0*p2*q0*t2 - p0*p0*q2*t2 + 
  2*n2*p0*s0*t2 + 2*n0*p2*s0*t2 - m2*s0*s0*t2 + 
  2*n0*p0*s2*t2 - 2*m0*s0*s2*t2 + 2*o2*p2*q0*u0 + 
  2*o2*p0*q2*u0 + 2*o0*p2*q2*u0 - 2*n2*p2*r0*u0 - 
  2*n2*p0*r2*u0 - 2*n0*p2*r2*u0 - 2*n2*o2*s0*u0 + 
  2*m2*r2*s0*u0 - 2*n2*o0*s2*u0 - 2*n0*o2*s2*u0 + 
  2*m2*r0*s2*u0 + 2*m0*r2*s2*u0 + n2*n2*u0*u0 - m2*q2*u0*u0 + 
  2*o2*p0*q0*u2 + 2*o0*p2*q0*u2 + 2*o0*p0*q2*u2 - 
  2*n2*p0*r0*u2 - 2*n0*p2*r0*u2 - 2*n0*p0*r2*u2 - 
  2*n2*o0*s0*u2 - 2*n0*o2*s0*u2 + 2*m2*r0*s0*u2 + 
  2*m0*r2*s0*u2 - 2*n0*o0*s2*u2 + 2*m0*r0*s2*u2 + 
  4*n0*n2*u0*u2 - 2*m2*q0*u0*u2 - 2*m0*q2*u0*u2 + n0*n0*u2*u2 - 
  m0*q0*u2*u2 - o2*o2*q0*v0 - 2*o0*o2*q2*v0 + 2*n2*o2*r0*v0 + 
  2*n2*o0*r2*v0 + 2*n0*o2*r2*v0 - 2*m2*r0*r2*v0 - 
  m0*r2*r2*v0 - n2*n2*t0*v0 + m2*q2*t0*v0 - 2*n0*n2*t2*v0 + 
  m2*q0*t2*v0 + m0*q2*t2*v0 - 2*o0*o2*q0*v2 - o0*o0*q2*v2 + 
  2*n2*o0*r0*v2 + 2*n0*o2*r0*v2 - m2*r0*r0*v2 + 
  2*n0*o0*r2*v2 - 2*m0*r0*r2*v2 - 2*n0*n2*t0*v2 + 
  m2*q0*t0*v2 + m0*q2*t0*v2 - n0*n0*t2*v2 + m0*q0*t2*v2;
    
    k[4]  = 2*p2*p2*r0*r1 + 4*p1*p2*r0*r2 + 
  4*p0*p2*r1*r2 + 2*p0*p1*r2*r2 - 2*o2*p2*r1*s0 - 
  2*o2*p1*r2*s0 - 2*o1*p2*r2*s0 - 2*o2*p2*r0*s1 - 
  2*o2*p0*r2*s1 - 2*o0*p2*r2*s1 + 2*o2*o2*s0*s1 - 
  2*o2*p1*r0*s2 - 2*o1*p2*r0*s2 - 2*o2*p0*r1*s2 - 
  2*o0*p2*r1*s2 - 2*o1*p0*r2*s2 - 2*o0*p1*r2*s2 + 
  4*o1*o2*s0*s2 + 4*o0*o2*s1*s2 + 2*o0*o1*s2*s2 - p2*p2*q1*t0 - 
  2*p1*p2*q2*t0 + 2*n2*p2*s1*t0 + 2*n2*p1*s2*t0 + 
  2*n1*p2*s2*t0 - 2*m2*s1*s2*t0 - m1*s2*s2*t0 - p2*p2*q0*t1 - 
  2*p0*p2*q2*t1 + 2*n2*p2*s0*t1 + 2*n2*p0*s2*t1 + 
  2*n0*p2*s2*t1 - 2*m2*s0*s2*t1 - m0*s2*s2*t1 - 
  2*p1*p2*q0*t2 - 2*p0*p2*q1*t2 - 2*p0*p1*q2*t2 + 
  2*n2*p1*s0*t2 + 2*n1*p2*s0*t2 + 2*n2*p0*s1*t2 + 
  2*n0*p2*s1*t2 - 2*m2*s0*s1*t2 + 2*n1*p0*s2*t2 + 
  2*n0*p1*s2*t2 - 2*m1*s0*s2*t2 - 2*m0*s1*s2*t2 + 
  2*o2*p2*q1*u0 + 2*o2*p1*q2*u0 + 2*o1*p2*q2*u0 - 
  2*n2*p2*r1*u0 - 2*n2*p1*r2*u0 - 2*n1*p2*r2*u0 - 
  2*n2*o2*s1*u0 + 2*m2*r2*s1*u0 - 2*n2*o1*s2*u0 - 
  2*n1*o2*s2*u0 + 2*m2*r1*s2*u0 + 2*m1*r2*s2*u0 + 
  2*o2*p2*q0*u1 + 2*o2*p0*q2*u1 + 2*o0*p2*q2*u1 - 
  2*n2*p2*r0*u1 - 2*n2*p0*r2*u1 - 2*n0*p2*r2*u1 - 
  2*n2*o2*s0*u1 + 2*m2*r2*s0*u1 - 2*n2*o0*s2*u1 - 
  2*n0*o2*s2*u1 + 2*m2*r0*s2*u1 + 2*m0*r2*s2*u1 + 
  2*n2*n2*u0*u1 - 2*m2*q2*u0*u1 + 2*o2*p1*q0*u2 + 
  2*o1*p2*q0*u2 + 2*o2*p0*q1*u2 + 2*o0*p2*q1*u2 + 
  2*o1*p0*q2*u2 + 2*o0*p1*q2*u2 - 2*n2*p1*r0*u2 - 
  2*n1*p2*r0*u2 - 2*n2*p0*r1*u2 - 2*n0*p2*r1*u2 - 
  2*n1*p0*r2*u2 - 2*n0*p1*r2*u2 - 2*n2*o1*s0*u2 - 
  2*n1*o2*s0*u2 + 2*m2*r1*s0*u2 + 2*m1*r2*s0*u2 - 
  2*n2*o0*s1*u2 - 2*n0*o2*s1*u2 + 2*m2*r0*s1*u2 + 
  2*m0*r2*s1*u2 - 2*n1*o0*s2*u2 - 2*n0*o1*s2*u2 + 
  2*m1*r0*s2*u2 + 2*m0*r1*s2*u2 + 4*n1*n2*u0*u2 - 
  2*m2*q1*u0*u2 - 2*m1*q2*u0*u2 + 4*n0*n2*u1*u2 - 
  2*m2*q0*u1*u2 - 2*m0*q2*u1*u2 + 2*n0*n1*u2*u2 - m1*q0*u2*u2 - 
  m0*q1*u2*u2 - o2*o2*q1*v0 - 2*o1*o2*q2*v0 + 2*n2*o2*r1*v0 + 
  2*n2*o1*r2*v0 + 2*n1*o2*r2*v0 - 2*m2*r1*r2*v0 - 
  m1*r2*r2*v0 - n2*n2*t1*v0 + m2*q2*t1*v0 - 2*n1*n2*t2*v0 + 
  m2*q1*t2*v0 + m1*q2*t2*v0 - o2*o2*q0*v1 - 2*o0*o2*q2*v1 + 
  2*n2*o2*r0*v1 + 2*n2*o0*r2*v1 + 2*n0*o2*r2*v1 - 
  2*m2*r0*r2*v1 - m0*r2*r2*v1 - n2*n2*t0*v1 + m2*q2*t0*v1 - 
  2*n0*n2*t2*v1 + m2*q0*t2*v1 + m0*q2*t2*v1 - 2*o1*o2*q0*v2 - 
  2*o0*o2*q1*v2 - 2*o0*o1*q2*v2 + 2*n2*o1*r0*v2 + 
  2*n1*o2*r0*v2 + 2*n2*o0*r1*v2 + 2*n0*o2*r1*v2 - 
  2*m2*r0*r1*v2 + 2*n1*o0*r2*v2 + 2*n0*o1*r2*v2 - 
  2*m1*r0*r2*v2 - 2*m0*r1*r2*v2 - 2*n1*n2*t0*v2 + 
  m2*q1*t0*v2 + m1*q2*t0*v2 - 2*n0*n2*t1*v2 + m2*q0*t1*v2 + 
  m0*q2*t1*v2 - 2*n0*n1*t2*v2 + m1*q0*t2*v2 + m0*q1*t2*v2;
    
    k[5]  = p2*p2*r1*r1 + 4*p1*p2*r1*r2 + 
  p1*p1*r2*r2 - 2*o2*p2*r1*s1 - 2*o2*p1*r2*s1 - 2*o1*p2*r2*s1 + 
  o2*o2*s1*s1 - 2*o2*p1*r1*s2 - 2*o1*p2*r1*s2 - 2*o1*p1*r2*s2 + 
  4*o1*o2*s1*s2 + o1*o1*s2*s2 - p2*p2*q1*t1 - 2*p1*p2*q2*t1 + 
  2*n2*p2*s1*t1 + 2*n2*p1*s2*t1 + 2*n1*p2*s2*t1 - 
  2*m2*s1*s2*t1 - m1*s2*s2*t1 - 2*p1*p2*q1*t2 - p1*p1*q2*t2 + 
  2*n2*p1*s1*t2 + 2*n1*p2*s1*t2 - m2*s1*s1*t2 + 
  2*n1*p1*s2*t2 - 2*m1*s1*s2*t2 + 2*o2*p2*q1*u1 + 
  2*o2*p1*q2*u1 + 2*o1*p2*q2*u1 - 2*n2*p2*r1*u1 - 
  2*n2*p1*r2*u1 - 2*n1*p2*r2*u1 - 2*n2*o2*s1*u1 + 
  2*m2*r2*s1*u1 - 2*n2*o1*s2*u1 - 2*n1*o2*s2*u1 + 
  2*m2*r1*s2*u1 + 2*m1*r2*s2*u1 + n2*n2*u1*u1 - m2*q2*u1*u1 + 
  2*o2*p1*q1*u2 + 2*o1*p2*q1*u2 + 2*o1*p1*q2*u2 - 
  2*n2*p1*r1*u2 - 2*n1*p2*r1*u2 - 2*n1*p1*r2*u2 - 
  2*n2*o1*s1*u2 - 2*n1*o2*s1*u2 + 2*m2*r1*s1*u2 + 
  2*m1*r2*s1*u2 - 2*n1*o1*s2*u2 + 2*m1*r1*s2*u2 + 
  4*n1*n2*u1*u2 - 2*m2*q1*u1*u2 - 2*m1*q2*u1*u2 + n1*n1*u2*u2 - 
  m1*q1*u2*u2 - o2*o2*q1*v1 - 2*o1*o2*q2*v1 + 2*n2*o2*r1*v1 + 
  2*n2*o1*r2*v1 + 2*n1*o2*r2*v1 - 2*m2*r1*r2*v1 - 
  m1*r2*r2*v1 - n2*n2*t1*v1 + m2*q2*t1*v1 - 2*n1*n2*t2*v1 + 
  m2*q1*t2*v1 + m1*q2*t2*v1 - 2*o1*o2*q1*v2 - o1*o1*q2*v2 + 
  2*n2*o1*r1*v2 + 2*n1*o2*r1*v2 - m2*r1*r1*v2 + 
  2*n1*o1*r2*v2 - 2*m1*r1*r2*v2 - 2*n1*n2*t1*v2 + 
  m2*q1*t1*v2 + m1*q2*t1*v2 - n1*n1*t2*v2 + m1*q1*t2*v2;
    
    k[6]  = 2*p0*p2*r0*r0 + 2*p0*p0*r0*r2 - 2*o2*p0*r0*s0 - 
  2*o0*p2*r0*s0 - 2*o0*p0*r2*s0 + 2*o0*o2*s0*s0 - 
  2*o0*p0*r0*s2 + 2*o0*o0*s0*s2 - 2*p0*p2*q0*t0 - p0*p0*q2*t0 + 
  2*n2*p0*s0*t0 + 2*n0*p2*s0*t0 - m2*s0*s0*t0 + 
  2*n0*p0*s2*t0 - 2*m0*s0*s2*t0 - p0*p0*q0*t2 + 
  2*n0*p0*s0*t2 - m0*s0*s0*t2 + 2*o2*p0*q0*u0 + 
  2*o0*p2*q0*u0 + 2*o0*p0*q2*u0 - 2*n2*p0*r0*u0 - 
  2*n0*p2*r0*u0 - 2*n0*p0*r2*u0 - 2*n2*o0*s0*u0 - 
  2*n0*o2*s0*u0 + 2*m2*r0*s0*u0 + 2*m0*r2*s0*u0 - 
  2*n0*o0*s2*u0 + 2*m0*r0*s2*u0 + 2*n0*n2*u0*u0 - m2*q0*u0*u0 - 
  m0*q2*u0*u0 + 2*o0*p0*q0*u2 - 2*n0*p0*r0*u2 - 
  2*n0*o0*s0*u2 + 2*m0*r0*s0*u2 + 2*n0*n0*u0*u2 - 
  2*m0*q0*u0*u2 - 2*o0*o2*q0*v0 - o0*o0*q2*v0 + 
  2*n2*o0*r0*v0 + 2*n0*o2*r0*v0 - m2*r0*r0*v0 + 
  2*n0*o0*r2*v0 - 2*m0*r0*r2*v0 - 2*n0*n2*t0*v0 + 
  m2*q0*t0*v0 + m0*q2*t0*v0 - n0*n0*t2*v0 + m0*q0*t2*v0 - 
  o0*o0*q0*v2 + 2*n0*o0*r0*v2 - m0*r0*r0*v2 - n0*n0*t0*v2 + 
  m0*q0*t0*v2;
    
    k[7]  = 2*p1*p2*r0*r0 + 4*p0*p2*r0*r1 + 4*p0*p1*r0*r2 + 
  2*p0*p0*r1*r2 - 2*o2*p1*r0*s0 - 2*o1*p2*r0*s0 - 
  2*o2*p0*r1*s0 - 2*o0*p2*r1*s0 - 2*o1*p0*r2*s0 - 
  2*o0*p1*r2*s0 + 2*o1*o2*s0*s0 - 2*o2*p0*r0*s1 - 
  2*o0*p2*r0*s1 - 2*o0*p0*r2*s1 + 4*o0*o2*s0*s1 - 
  2*o1*p0*r0*s2 - 2*o0*p1*r0*s2 - 2*o0*p0*r1*s2 + 
  4*o0*o1*s0*s2 + 2*o0*o0*s1*s2 - 2*p1*p2*q0*t0 - 
  2*p0*p2*q1*t0 - 2*p0*p1*q2*t0 + 2*n2*p1*s0*t0 + 
  2*n1*p2*s0*t0 + 2*n2*p0*s1*t0 + 2*n0*p2*s1*t0 - 
  2*m2*s0*s1*t0 + 2*n1*p0*s2*t0 + 2*n0*p1*s2*t0 - 
  2*m1*s0*s2*t0 - 2*m0*s1*s2*t0 - 2*p0*p2*q0*t1 - 
  p0*p0*q2*t1 + 2*n2*p0*s0*t1 + 2*n0*p2*s0*t1 - 
  m2*s0*s0*t1 + 2*n0*p0*s2*t1 - 2*m0*s0*s2*t1 - 
  2*p0*p1*q0*t2 - p0*p0*q1*t2 + 2*n1*p0*s0*t2 + 
  2*n0*p1*s0*t2 - m1*s0*s0*t2 + 2*n0*p0*s1*t2 - 
  2*m0*s0*s1*t2 + 2*o2*p1*q0*u0 + 2*o1*p2*q0*u0 + 
  2*o2*p0*q1*u0 + 2*o0*p2*q1*u0 + 2*o1*p0*q2*u0 + 
  2*o0*p1*q2*u0 - 2*n2*p1*r0*u0 - 2*n1*p2*r0*u0 - 
  2*n2*p0*r1*u0 - 2*n0*p2*r1*u0 - 2*n1*p0*r2*u0 - 
  2*n0*p1*r2*u0 - 2*n2*o1*s0*u0 - 2*n1*o2*s0*u0 + 
  2*m2*r1*s0*u0 + 2*m1*r2*s0*u0 - 2*n2*o0*s1*u0 - 
  2*n0*o2*s1*u0 + 2*m2*r0*s1*u0 + 2*m0*r2*s1*u0 - 
  2*n1*o0*s2*u0 - 2*n0*o1*s2*u0 + 2*m1*r0*s2*u0 + 
  2*m0*r1*s2*u0 + 2*n1*n2*u0*u0 - m2*q1*u0*u0 - 
  m1*q2*u0*u0 + 2*o2*p0*q0*u1 + 2*o0*p2*q0*u1 + 
  2*o0*p0*q2*u1 - 2*n2*p0*r0*u1 - 2*n0*p2*r0*u1 - 
  2*n0*p0*r2*u1 - 2*n2*o0*s0*u1 - 2*n0*o2*s0*u1 + 
  2*m2*r0*s0*u1 + 2*m0*r2*s0*u1 - 2*n0*o0*s2*u1 + 
  2*m0*r0*s2*u1 + 4*n0*n2*u0*u1 - 2*m2*q0*u0*u1 - 
  2*m0*q2*u0*u1 + 2*o1*p0*q0*u2 + 2*o0*p1*q0*u2 + 
  2*o0*p0*q1*u2 - 2*n1*p0*r0*u2 - 2*n0*p1*r0*u2 - 
  2*n0*p0*r1*u2 - 2*n1*o0*s0*u2 - 2*n0*o1*s0*u2 + 
  2*m1*r0*s0*u2 + 2*m0*r1*s0*u2 - 2*n0*o0*s1*u2 + 
  2*m0*r0*s1*u2 + 4*n0*n1*u0*u2 - 2*m1*q0*u0*u2 - 
  2*m0*q1*u0*u2 + 2*n0*n0*u1*u2 - 2*m0*q0*u1*u2 - 
  2*o1*o2*q0*v0 - 2*o0*o2*q1*v0 - 2*o0*o1*q2*v0 + 
  2*n2*o1*r0*v0 + 2*n1*o2*r0*v0 + 2*n2*o0*r1*v0 + 
  2*n0*o2*r1*v0 - 2*m2*r0*r1*v0 + 2*n1*o0*r2*v0 + 
  2*n0*o1*r2*v0 - 2*m1*r0*r2*v0 - 2*m0*r1*r2*v0 - 
  2*n1*n2*t0*v0 + m2*q1*t0*v0 + m1*q2*t0*v0 - 
  2*n0*n2*t1*v0 + m2*q0*t1*v0 + m0*q2*t1*v0 - 
  2*n0*n1*t2*v0 + m1*q0*t2*v0 + m0*q1*t2*v0 - 
  2*o0*o2*q0*v1 - o0*o0*q2*v1 + 2*n2*o0*r0*v1 + 
  2*n0*o2*r0*v1 - m2*r0*r0*v1 + 2*n0*o0*r2*v1 - 
  2*m0*r0*r2*v1 - 2*n0*n2*t0*v1 + m2*q0*t0*v1 + 
  m0*q2*t0*v1 - n0*n0*t2*v1 + m0*q0*t2*v1 - 
  2*o0*o1*q0*v2 - o0*o0*q1*v2 + 2*n1*o0*r0*v2 + 
  2*n0*o1*r0*v2 - m1*r0*r0*v2 + 2*n0*o0*r1*v2 - 
  2*m0*r0*r1*v2 - 2*n0*n1*t0*v2 + m1*q0*t0*v2 + 
  m0*q1*t0*v2 - n0*n0*t1*v2 + m0*q0*t1*v2;
    
    k[8]  = 4*p1*p2*r0*r1 + 2*p0*p2*r1*r1 + 2*p1*p1*r0*r2 + 
  4*p0*p1*r1*r2 - 2*o2*p1*r1*s0 - 2*o1*p2*r1*s0 - 
  2*o1*p1*r2*s0 - 2*o2*p1*r0*s1 - 2*o1*p2*r0*s1 - 
  2*o2*p0*r1*s1 - 2*o0*p2*r1*s1 - 2*o1*p0*r2*s1 - 
  2*o0*p1*r2*s1 + 4*o1*o2*s0*s1 + 2*o0*o2*s1*s1 - 
  2*o1*p1*r0*s2 - 2*o1*p0*r1*s2 - 2*o0*p1*r1*s2 + 
  2*o1*o1*s0*s2 + 4*o0*o1*s1*s2 - 2*p1*p2*q1*t0 - 
  p1*p1*q2*t0 + 2*n2*p1*s1*t0 + 2*n1*p2*s1*t0 - 
  m2*s1*s1*t0 + 2*n1*p1*s2*t0 - 2*m1*s1*s2*t0 - 
  2*p1*p2*q0*t1 - 2*p0*p2*q1*t1 - 2*p0*p1*q2*t1 + 
  2*n2*p1*s0*t1 + 2*n1*p2*s0*t1 + 2*n2*p0*s1*t1 + 
  2*n0*p2*s1*t1 - 2*m2*s0*s1*t1 + 2*n1*p0*s2*t1 + 
  2*n0*p1*s2*t1 - 2*m1*s0*s2*t1 - 2*m0*s1*s2*t1 - 
  p1*p1*q0*t2 - 2*p0*p1*q1*t2 + 2*n1*p1*s0*t2 + 
  2*n1*p0*s1*t2 + 2*n0*p1*s1*t2 - 2*m1*s0*s1*t2 - 
  m0*s1*s1*t2 + 2*o2*p1*q1*u0 + 2*o1*p2*q1*u0 + 
  2*o1*p1*q2*u0 - 2*n2*p1*r1*u0 - 2*n1*p2*r1*u0 - 
  2*n1*p1*r2*u0 - 2*n2*o1*s1*u0 - 2*n1*o2*s1*u0 + 
  2*m2*r1*s1*u0 + 2*m1*r2*s1*u0 - 2*n1*o1*s2*u0 + 
  2*m1*r1*s2*u0 + 2*o2*p1*q0*u1 + 2*o1*p2*q0*u1 + 
  2*o2*p0*q1*u1 + 2*o0*p2*q1*u1 + 2*o1*p0*q2*u1 + 
  2*o0*p1*q2*u1 - 2*n2*p1*r0*u1 - 2*n1*p2*r0*u1 - 
  2*n2*p0*r1*u1 - 2*n0*p2*r1*u1 - 2*n1*p0*r2*u1 - 
  2*n0*p1*r2*u1 - 2*n2*o1*s0*u1 - 2*n1*o2*s0*u1 + 
  2*m2*r1*s0*u1 + 2*m1*r2*s0*u1 - 2*n2*o0*s1*u1 - 
  2*n0*o2*s1*u1 + 2*m2*r0*s1*u1 + 2*m0*r2*s1*u1 - 
  2*n1*o0*s2*u1 - 2*n0*o1*s2*u1 + 2*m1*r0*s2*u1 + 
  2*m0*r1*s2*u1 + 4*n1*n2*u0*u1 - 2*m2*q1*u0*u1 - 
  2*m1*q2*u0*u1 + 2*n0*n2*u1*u1 - m2*q0*u1*u1 - 
  m0*q2*u1*u1 + 2*o1*p1*q0*u2 + 2*o1*p0*q1*u2 + 
  2*o0*p1*q1*u2 - 2*n1*p1*r0*u2 - 2*n1*p0*r1*u2 - 
  2*n0*p1*r1*u2 - 2*n1*o1*s0*u2 + 2*m1*r1*s0*u2 - 
  2*n1*o0*s1*u2 - 2*n0*o1*s1*u2 + 2*m1*r0*s1*u2 + 
  2*m0*r1*s1*u2 + 2*n1*n1*u0*u2 - 2*m1*q1*u0*u2 + 
  4*n0*n1*u1*u2 - 2*m1*q0*u1*u2 - 2*m0*q1*u1*u2 - 
  2*o1*o2*q1*v0 - o1*o1*q2*v0 + 2*n2*o1*r1*v0 + 
  2*n1*o2*r1*v0 - m2*r1*r1*v0 + 2*n1*o1*r2*v0 - 
  2*m1*r1*r2*v0 - 2*n1*n2*t1*v0 + m2*q1*t1*v0 + 
  m1*q2*t1*v0 - n1*n1*t2*v0 + m1*q1*t2*v0 - 
  2*o1*o2*q0*v1 - 2*o0*o2*q1*v1 - 2*o0*o1*q2*v1 + 
  2*n2*o1*r0*v1 + 2*n1*o2*r0*v1 + 2*n2*o0*r1*v1 + 
  2*n0*o2*r1*v1 - 2*m2*r0*r1*v1 + 2*n1*o0*r2*v1 + 
  2*n0*o1*r2*v1 - 2*m1*r0*r2*v1 - 2*m0*r1*r2*v1 - 
  2*n1*n2*t0*v1 + m2*q1*t0*v1 + m1*q2*t0*v1 - 
  2*n0*n2*t1*v1 + m2*q0*t1*v1 + m0*q2*t1*v1 - 
  2*n0*n1*t2*v1 + m1*q0*t2*v1 + m0*q1*t2*v1 - 
  o1*o1*q0*v2 - 2*o0*o1*q1*v2 + 2*n1*o1*r0*v2 + 
  2*n1*o0*r1*v2 + 2*n0*o1*r1*v2 - 2*m1*r0*r1*v2 - 
  m0*r1*r1*v2 - n1*n1*t0*v2 + m1*q1*t0*v2 - 
  2*n0*n1*t1*v2 + m1*q0*t1*v2 + m0*q1*t1*v2;
    
    k[9]  = 2*p1*p2*r1*r1 + 2*p1*p1*r1*r2 - 2*o2*p1*r1*s1 - 
  2*o1*p2*r1*s1 - 2*o1*p1*r2*s1 + 2*o1*o2*s1*s1 - 
  2*o1*p1*r1*s2 + 2*o1*o1*s1*s2 - 2*p1*p2*q1*t1 - p1*p1*q2*t1 + 
  2*n2*p1*s1*t1 + 2*n1*p2*s1*t1 - m2*s1*s1*t1 + 
  2*n1*p1*s2*t1 - 2*m1*s1*s2*t1 - p1*p1*q1*t2 + 
  2*n1*p1*s1*t2 - m1*s1*s1*t2 + 2*o2*p1*q1*u1 + 
  2*o1*p2*q1*u1 + 2*o1*p1*q2*u1 - 2*n2*p1*r1*u1 - 
  2*n1*p2*r1*u1 - 2*n1*p1*r2*u1 - 2*n2*o1*s1*u1 - 
  2*n1*o2*s1*u1 + 2*m2*r1*s1*u1 + 2*m1*r2*s1*u1 - 
  2*n1*o1*s2*u1 + 2*m1*r1*s2*u1 + 2*n1*n2*u1*u1 - m2*q1*u1*u1 - 
  m1*q2*u1*u1 + 2*o1*p1*q1*u2 - 2*n1*p1*r1*u2 - 
  2*n1*o1*s1*u2 + 2*m1*r1*s1*u2 + 2*n1*n1*u1*u2 - 
  2*m1*q1*u1*u2 - 2*o1*o2*q1*v1 - o1*o1*q2*v1 + 
  2*n2*o1*r1*v1 + 2*n1*o2*r1*v1 - m2*r1*r1*v1 + 
  2*n1*o1*r2*v1 - 2*m1*r1*r2*v1 - 2*n1*n2*t1*v1 + 
  m2*q1*t1*v1 + m1*q2*t1*v1 - n1*n1*t2*v1 + m1*q1*t2*v1 - 
  o1*o1*q1*v2 + 2*n1*o1*r1*v2 - m1*r1*r1*v2 - n1*n1*t1*v2 + 
  m1*q1*t1*v2;
    
    k[10] = p0*p0*r0*r0 - 2*o0*p0*r0*s0 + o0*o0*s0*s0 - 
  p0*p0*q0*t0 + 2*n0*p0*s0*t0 - m0*s0*s0*t0 + 2*o0*p0*q0*u0 - 
  2*n0*p0*r0*u0 - 2*n0*o0*s0*u0 + 2*m0*r0*s0*u0 + n0*n0*u0*u0 - 
  m0*q0*u0*u0 - o0*o0*q0*v0 + 2*n0*o0*r0*v0 - m0*r0*r0*v0 - 
  n0*n0*t0*v0 + m0*q0*t0*v0;
    
    k[11] = 2*p0*p1*r0*r0 + 2*p0*p0*r0*r1 - 2*o1*p0*r0*s0 - 
  2*o0*p1*r0*s0 - 2*o0*p0*r1*s0 + 2*o0*o1*s0*s0 - 
  2*o0*p0*r0*s1 + 2*o0*o0*s0*s1 - 2*p0*p1*q0*t0 - 
  p0*p0*q1*t0 + 2*n1*p0*s0*t0 + 2*n0*p1*s0*t0 - 
  m1*s0*s0*t0 + 2*n0*p0*s1*t0 - 2*m0*s0*s1*t0 - 
  p0*p0*q0*t1 + 2*n0*p0*s0*t1 - m0*s0*s0*t1 + 
  2*o1*p0*q0*u0 + 2*o0*p1*q0*u0 + 2*o0*p0*q1*u0 - 
  2*n1*p0*r0*u0 - 2*n0*p1*r0*u0 - 2*n0*p0*r1*u0 - 
  2*n1*o0*s0*u0 - 2*n0*o1*s0*u0 + 2*m1*r0*s0*u0 + 
  2*m0*r1*s0*u0 - 2*n0*o0*s1*u0 + 2*m0*r0*s1*u0 + 
  2*n0*n1*u0*u0 - m1*q0*u0*u0 - m0*q1*u0*u0 + 
  2*o0*p0*q0*u1 - 2*n0*p0*r0*u1 - 2*n0*o0*s0*u1 + 
  2*m0*r0*s0*u1 + 2*n0*n0*u0*u1 - 2*m0*q0*u0*u1 - 
  2*o0*o1*q0*v0 - o0*o0*q1*v0 + 2*n1*o0*r0*v0 + 
  2*n0*o1*r0*v0 - m1*r0*r0*v0 + 2*n0*o0*r1*v0 - 
  2*m0*r0*r1*v0 - 2*n0*n1*t0*v0 + m1*q0*t0*v0 + 
  m0*q1*t0*v0 - n0*n0*t1*v0 + m0*q0*t1*v0 - 
  o0*o0*q0*v1 + 2*n0*o0*r0*v1 - m0*r0*r0*v1 - 
  n0*n0*t0*v1 + m0*q0*t0*v1;
    
    k[12] = p1*p1*r0*r0 + 4*p0*p1*r0*r1 + p0*p0*r1*r1 - 
  2*o1*p1*r0*s0 - 2*o1*p0*r1*s0 - 2*o0*p1*r1*s0 + 
  o1*o1*s0*s0 - 2*o1*p0*r0*s1 - 2*o0*p1*r0*s1 - 
  2*o0*p0*r1*s1 + 4*o0*o1*s0*s1 + o0*o0*s1*s1 - 
  p1*p1*q0*t0 - 2*p0*p1*q1*t0 + 2*n1*p1*s0*t0 + 
  2*n1*p0*s1*t0 + 2*n0*p1*s1*t0 - 2*m1*s0*s1*t0 - 
  m0*s1*s1*t0 - 2*p0*p1*q0*t1 - p0*p0*q1*t1 + 
  2*n1*p0*s0*t1 + 2*n0*p1*s0*t1 - m1*s0*s0*t1 + 
  2*n0*p0*s1*t1 - 2*m0*s0*s1*t1 + 2*o1*p1*q0*u0 + 
  2*o1*p0*q1*u0 + 2*o0*p1*q1*u0 - 2*n1*p1*r0*u0 - 
  2*n1*p0*r1*u0 - 2*n0*p1*r1*u0 - 2*n1*o1*s0*u0 + 
  2*m1*r1*s0*u0 - 2*n1*o0*s1*u0 - 2*n0*o1*s1*u0 + 
  2*m1*r0*s1*u0 + 2*m0*r1*s1*u0 + n1*n1*u0*u0 - 
  m1*q1*u0*u0 + 2*o1*p0*q0*u1 + 2*o0*p1*q0*u1 + 
  2*o0*p0*q1*u1 - 2*n1*p0*r0*u1 - 2*n0*p1*r0*u1 - 
  2*n0*p0*r1*u1 - 2*n1*o0*s0*u1 - 2*n0*o1*s0*u1 + 
  2*m1*r0*s0*u1 + 2*m0*r1*s0*u1 - 2*n0*o0*s1*u1 + 
  2*m0*r0*s1*u1 + 4*n0*n1*u0*u1 - 2*m1*q0*u0*u1 - 
  2*m0*q1*u0*u1 + n0*n0*u1*u1 - m0*q0*u1*u1 - 
  o1*o1*q0*v0 - 2*o0*o1*q1*v0 + 2*n1*o1*r0*v0 + 
  2*n1*o0*r1*v0 + 2*n0*o1*r1*v0 - 2*m1*r0*r1*v0 - 
  m0*r1*r1*v0 - n1*n1*t0*v0 + m1*q1*t0*v0 - 
  2*n0*n1*t1*v0 + m1*q0*t1*v0 + m0*q1*t1*v0 - 
  2*o0*o1*q0*v1 - o0*o0*q1*v1 + 2*n1*o0*r0*v1 + 
  2*n0*o1*r0*v1 - m1*r0*r0*v1 + 2*n0*o0*r1*v1 - 
  2*m0*r0*r1*v1 - 2*n0*n1*t0*v1 + m1*q0*t0*v1 + 
  m0*q1*t0*v1 - n0*n0*t1*v1 + m0*q0*t1*v1;
    
    k[13] = 2*p1*p1*r0*r1 + 2*p0*p1*r1*r1 - 
  2*o1*p1*r1*s0 - 2*o1*p1*r0*s1 - 2*o1*p0*r1*s1 - 
  2*o0*p1*r1*s1 + 2*o1*o1*s0*s1 + 2*o0*o1*s1*s1 - 
  p1*p1*q1*t0 + 2*n1*p1*s1*t0 - m1*s1*s1*t0 - 
  p1*p1*q0*t1 - 2*p0*p1*q1*t1 + 2*n1*p1*s0*t1 + 
  2*n1*p0*s1*t1 + 2*n0*p1*s1*t1 - 2*m1*s0*s1*t1 - 
  m0*s1*s1*t1 + 2*o1*p1*q1*u0 - 2*n1*p1*r1*u0 - 
  2*n1*o1*s1*u0 + 2*m1*r1*s1*u0 + 2*o1*p1*q0*u1 + 
  2*o1*p0*q1*u1 + 2*o0*p1*q1*u1 - 2*n1*p1*r0*u1 - 
  2*n1*p0*r1*u1 - 2*n0*p1*r1*u1 - 2*n1*o1*s0*u1 + 
  2*m1*r1*s0*u1 - 2*n1*o0*s1*u1 - 2*n0*o1*s1*u1 + 
  2*m1*r0*s1*u1 + 2*m0*r1*s1*u1 + 2*n1*n1*u0*u1 - 
  2*m1*q1*u0*u1 + 2*n0*n1*u1*u1 - m1*q0*u1*u1 - 
  m0*q1*u1*u1 - o1*o1*q1*v0 + 2*n1*o1*r1*v0 - 
  m1*r1*r1*v0 - n1*n1*t1*v0 + m1*q1*t1*v0 - 
  o1*o1*q0*v1 - 2*o0*o1*q1*v1 + 2*n1*o1*r0*v1 + 
  2*n1*o0*r1*v1 + 2*n0*o1*r1*v1 - 2*m1*r0*r1*v1 - 
  m0*r1*r1*v1 - n1*n1*t0*v1 + m1*q1*t0*v1 - 
  2*n0*n1*t1*v1 + m1*q0*t1*v1 + m0*q1*t1*v1;
    
    k[14] = p1*p1*r1*r1 - 2*o1*p1*r1*s1 + o1*o1*s1*s1 - p1*p1*q1*t1 + 
  2*n1*p1*s1*t1 - m1*s1*s1*t1 + 2*o1*p1*q1*u1 - 
  2*n1*p1*r1*u1 - 2*n1*o1*s1*u1 + 2*m1*r1*s1*u1 + n1*n1*u1*u1 - 
  m1*q1*u1*u1 - o1*o1*q1*v1 + 2*n1*o1*r1*v1 - m1*r1*r1*v1 - 
  n1*n1*t1*v1 + m1*q1*t1*v1;
   
	return 	rg_ImplicitEquation(4, k);
}

rg_REAL rg_IntsecImplicit4::inversion(const rg_BzCurve2D &curve, const rg_Point2D &point)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL x3 = curve.getCtrlPt(3).getX();
	rg_REAL x4 = curve.getCtrlPt(4).getX();
	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL y3 = curve.getCtrlPt(3).getY();
	rg_REAL y4 = curve.getCtrlPt(4).getY();

//  If
//		    
//	x(t) = a4 t^4 + a3 t^3 + a2 t^2 + a1 t + a0
//
//	y(t) = b4 t^4 + b3 t^3 + b2 t^2 + b1 t + b0
//  
//  then the resultant of p(x, t) and q(y, t) is like the following.

	rg_REAL a4 = x4 - 4*x3 +  6*x2 -  4*x1 +   x0;
	rg_REAL a3 =      4*x3 - 12*x2 + 12*x1 - 4*x0;
	rg_REAL a2 =              6*x2 - 12*x1 + 6*x0;
	rg_REAL a1 =                      4*x1 - 4*x0;
	rg_REAL a0 =                               x0;

	rg_REAL b4 = y4 - 4*y3 +  6*y2 -  4*y1 +   y0;
	rg_REAL b3 =      4*y3 - 12*y2 + 12*y1 - 4*y0;
	rg_REAL b2 =              6*y2 - 12*y1 + 6*y0;
	rg_REAL b1 =                      4*y1 - 4*y0;
	rg_REAL b0 =                               y0;

  //			|															                   |
  //			|	m0 x + m1 y + m2	n0 x + n1 y + n2	o0 x + o1 y + o2  p0 x + p1 y + p2 |
  //			|															                   |
  // R(P, Q) =	|	n0 x + n1 y + n2	q0 x + q1 y + q2	r0 x + r1 y + r2  s0 x + s1 y + s2 |
  //			|															                   |
  //			|	o0 x + o1 y + o2	r0 x + r1 y + r2	t0 x + t1 y + t2  u0 x + u1 y + u2 |
  //			|															                   |
  //            |   p0 x + p1 y + p2    s0 x + s1 y + s2    u0 x + u1 y + u2  v0 x + v1 y + v2 |
  //            |                                                                              |
  
	rg_REAL m0 =   0.0;
	rg_REAL m1 =   0.0;
	rg_REAL m2 =   a4*b3 - a3*b4;

  	rg_REAL n0 =   0.0;
	rg_REAL n1 =   0.0;
	rg_REAL n2 =   a4*b2 - a2*b4;

  	rg_REAL o0 =   0.0;
	rg_REAL o1 =   0.0;
	rg_REAL o2 =   a4*b1 - a1*b4;

  	rg_REAL p0 =   b4;
	rg_REAL p1 = - a4;
	rg_REAL p2 =   a4*b0 - a0*b4;

   	rg_REAL q0 =   0.0;
	rg_REAL q1 =   0.0;
	rg_REAL q2 =   a4*b1 + a3*b2 - a2*b3 - a1*b4;

   	rg_REAL r0 =   b4;
	rg_REAL r1 = - a4;
	rg_REAL r2 =   a4*b0 + a3*b1 - a1*b3 - a0*b4;

	rg_REAL s0 =   b3;
	rg_REAL s1 = - a3;
	rg_REAL s2 =   a3*b0 - a0*b3;

   	rg_REAL t0 =   b3;
	rg_REAL t1 = - a3;
	rg_REAL t2 =   a3*b0 + a2*b1 - a1*b2 - a0*b3;

	rg_REAL u0 =   b2;
	rg_REAL u1 = - a2;
	rg_REAL u2 =   a2*b0 - a0*b2;

	rg_REAL v0 =   b1;
	rg_REAL v1 = - a1;
	rg_REAL v2 =   a1*b0 - a0*b1;

// 고쳐야할 부분
	rg_REAL x  =  point.getX();
	rg_REAL y  =  point.getY();

	rg_REAL t_3 = -(p2*r2*r2) + o2*r2*s2 + p2*q2*t2 - n2*s2*t2 - o2*q2*u2 + n2*r2*u2 - 
  2*p2*r0*r2*x - p0*r2*r2*x + o2*r2*s0*x + o2*r0*s2*x + o0*r2*s2*x + 
  p2*q2*t0*x - n2*s2*t0*x + p2*q0*t2*x + p0*q2*t2*x - n2*s0*t2*x - 
  n0*s2*t2*x - o2*q2*u0*x + n2*r2*u0*x - o2*q0*u2*x - o0*q2*u2*x + 
  n2*r0*u2*x + n0*r2*u2*x - p2*r0*r0*x*x - 2*p0*r0*r2*x*x + o2*r0*s0*x*x + 
  o0*r2*s0*x*x + o0*r0*s2*x*x + p2*q0*t0*x*x + p0*q2*t0*x*x - n2*s0*t0*x*x - 
  n0*s2*t0*x*x + p0*q0*t2*x*x - n0*s0*t2*x*x - o2*q0*u0*x*x - o0*q2*u0*x*x + 
  n2*r0*u0*x*x + n0*r2*u0*x*x - o0*q0*u2*x*x + n0*r0*u2*x*x - p0*r0*r0*x*x*x + 
  o0*r0*s0*x*x*x + p0*q0*t0*x*x*x - n0*s0*t0*x*x*x - o0*q0*u0*x*x*x + n0*r0*u0*x*x*x - 
  2*p2*r1*r2*y - p1*r2*r2*y + o2*r2*s1*y + o2*r1*s2*y + o1*r2*s2*y + 
  p2*q2*t1*y - n2*s2*t1*y + p2*q1*t2*y + p1*q2*t2*y - n2*s1*t2*y - 
  n1*s2*t2*y - o2*q2*u1*y + n2*r2*u1*y - o2*q1*u2*y - o1*q2*u2*y + 
  n2*r1*u2*y + n1*r2*u2*y - 2*p2*r0*r1*x*y - 2*p1*r0*r2*x*y - 
  2*p0*r1*r2*x*y + o2*r1*s0*x*y + o1*r2*s0*x*y + o2*r0*s1*x*y + 
  o0*r2*s1*x*y + o1*r0*s2*x*y + o0*r1*s2*x*y + p2*q1*t0*x*y + p1*q2*t0*x*y - 
  n2*s1*t0*x*y - n1*s2*t0*x*y + p2*q0*t1*x*y + p0*q2*t1*x*y - n2*s0*t1*x*y - 
  n0*s2*t1*x*y + p1*q0*t2*x*y + p0*q1*t2*x*y - n1*s0*t2*x*y - n0*s1*t2*x*y - 
  o2*q1*u0*x*y - o1*q2*u0*x*y + n2*r1*u0*x*y + n1*r2*u0*x*y - o2*q0*u1*x*y - 
  o0*q2*u1*x*y + n2*r0*u1*x*y + n0*r2*u1*x*y - o1*q0*u2*x*y - o0*q1*u2*x*y + 
  n1*r0*u2*x*y + n0*r1*u2*x*y - p1*r0*r0*x*x*y - 2*p0*r0*r1*x*x*y + 
  o1*r0*s0*x*x*y + o0*r1*s0*x*x*y + o0*r0*s1*x*x*y + p1*q0*t0*x*x*y + 
  p0*q1*t0*x*x*y - n1*s0*t0*x*x*y - n0*s1*t0*x*x*y + p0*q0*t1*x*x*y - 
  n0*s0*t1*x*x*y - o1*q0*u0*x*x*y - o0*q1*u0*x*x*y + n1*r0*u0*x*x*y + 
  n0*r1*u0*x*x*y - o0*q0*u1*x*x*y + n0*r0*u1*x*x*y - p2*r1*r1*y*y - 
  2*p1*r1*r2*y*y + o2*r1*s1*y*y + o1*r2*s1*y*y + o1*r1*s2*y*y + 
  p2*q1*t1*y*y + p1*q2*t1*y*y - n2*s1*t1*y*y - n1*s2*t1*y*y + p1*q1*t2*y*y - 
  n1*s1*t2*y*y - o2*q1*u1*y*y - o1*q2*u1*y*y + n2*r1*u1*y*y + n1*r2*u1*y*y - 
  o1*q1*u2*y*y + n1*r1*u2*y*y - 2*p1*r0*r1*x*y*y - p0*r1*r1*x*y*y + 
  o1*r1*s0*x*y*y + o1*r0*s1*x*y*y + o0*r1*s1*x*y*y + p1*q1*t0*x*y*y - 
  n1*s1*t0*x*y*y + p1*q0*t1*x*y*y + p0*q1*t1*x*y*y - n1*s0*t1*x*y*y - 
  n0*s1*t1*x*y*y - o1*q1*u0*x*y*y + n1*r1*u0*x*y*y - o1*q0*u1*x*y*y - 
  o0*q1*u1*x*y*y + n1*r0*u1*x*y*y + n0*r1*u1*x*y*y - p1*r1*r1*y*y*y + 
  o1*r1*s1*y*y*y + p1*q1*t1*y*y*y - n1*s1*t1*y*y*y - o1*q1*u1*y*y*y + n1*r1*u1*y*y*y;

	rg_REAL t_2  = -(o2*p2*r2) + o2*o2*s2 + n2*p2*t2 - m2*s2*t2 - n2*o2*u2 + m2*r2*u2 - 
  o2*p2*r0*x - o2*p0*r2*x - o0*p2*r2*x + o2*o2*s0*x + 2*o0*o2*s2*x + 
  n2*p2*t0*x - m2*s2*t0*x + n2*p0*t2*x + n0*p2*t2*x - m2*s0*t2*x - 
  m0*s2*t2*x - n2*o2*u0*x + m2*r2*u0*x - n2*o0*u2*x - n0*o2*u2*x + 
  m2*r0*u2*x + m0*r2*u2*x - o2*p0*r0*x*x - o0*p2*r0*x*x - o0*p0*r2*x*x + 
  2*o0*o2*s0*x*x + o0*o0*s2*x*x + n2*p0*t0*x*x + n0*p2*t0*x*x - m2*s0*t0*x*x - 
  m0*s2*t0*x*x + n0*p0*t2*x*x - m0*s0*t2*x*x - n2*o0*u0*x*x - n0*o2*u0*x*x + 
  m2*r0*u0*x*x + m0*r2*u0*x*x - n0*o0*u2*x*x + m0*r0*u2*x*x - o0*p0*r0*x*x*x + 
  o0*o0*s0*x*x*x + n0*p0*t0*x*x*x - m0*s0*t0*x*x*x - n0*o0*u0*x*x*x + m0*r0*u0*x*x*x - 
  o2*p2*r1*y - o2*p1*r2*y - o1*p2*r2*y + o2*o2*s1*y + 2*o1*o2*s2*y + 
  n2*p2*t1*y - m2*s2*t1*y + n2*p1*t2*y + n1*p2*t2*y - m2*s1*t2*y - 
  m1*s2*t2*y - n2*o2*u1*y + m2*r2*u1*y - n2*o1*u2*y - n1*o2*u2*y + 
  m2*r1*u2*y + m1*r2*u2*y - o2*p1*r0*x*y - o1*p2*r0*x*y - o2*p0*r1*x*y - 
  o0*p2*r1*x*y - o1*p0*r2*x*y - o0*p1*r2*x*y + 2*o1*o2*s0*x*y + 
  2*o0*o2*s1*x*y + 2*o0*o1*s2*x*y + n2*p1*t0*x*y + n1*p2*t0*x*y - 
  m2*s1*t0*x*y - m1*s2*t0*x*y + n2*p0*t1*x*y + n0*p2*t1*x*y - m2*s0*t1*x*y - 
  m0*s2*t1*x*y + n1*p0*t2*x*y + n0*p1*t2*x*y - m1*s0*t2*x*y - m0*s1*t2*x*y - 
  n2*o1*u0*x*y - n1*o2*u0*x*y + m2*r1*u0*x*y + m1*r2*u0*x*y - n2*o0*u1*x*y - 
  n0*o2*u1*x*y + m2*r0*u1*x*y + m0*r2*u1*x*y - n1*o0*u2*x*y - n0*o1*u2*x*y + 
  m1*r0*u2*x*y + m0*r1*u2*x*y - o1*p0*r0*x*x*y - o0*p1*r0*x*x*y - 
  o0*p0*r1*x*x*y + 2*o0*o1*s0*x*x*y + o0*o0*s1*x*x*y + n1*p0*t0*x*x*y + 
  n0*p1*t0*x*x*y - m1*s0*t0*x*x*y - m0*s1*t0*x*x*y + n0*p0*t1*x*x*y - 
  m0*s0*t1*x*x*y - n1*o0*u0*x*x*y - n0*o1*u0*x*x*y + m1*r0*u0*x*x*y + 
  m0*r1*u0*x*x*y - n0*o0*u1*x*x*y + m0*r0*u1*x*x*y - o2*p1*r1*y*y - 
  o1*p2*r1*y*y - o1*p1*r2*y*y + 2*o1*o2*s1*y*y + o1*o1*s2*y*y + n2*p1*t1*y*y + 
  n1*p2*t1*y*y - m2*s1*t1*y*y - m1*s2*t1*y*y + n1*p1*t2*y*y - m1*s1*t2*y*y - 
  n2*o1*u1*y*y - n1*o2*u1*y*y + m2*r1*u1*y*y + m1*r2*u1*y*y - n1*o1*u2*y*y + 
  m1*r1*u2*y*y - o1*p1*r0*x*y*y - o1*p0*r1*x*y*y - o0*p1*r1*x*y*y + 
  o1*o1*s0*x*y*y + 2*o0*o1*s1*x*y*y + n1*p1*t0*x*y*y - m1*s1*t0*x*y*y + 
  n1*p0*t1*x*y*y + n0*p1*t1*x*y*y - m1*s0*t1*x*y*y - m0*s1*t1*x*y*y - 
  n1*o1*u0*x*y*y + m1*r1*u0*x*y*y - n1*o0*u1*x*y*y - n0*o1*u1*x*y*y + 
  m1*r0*u1*x*y*y + m0*r1*u1*x*y*y - o1*p1*r1*y*y*y + o1*o1*s1*y*y*y + 
  n1*p1*t1*y*y*y - m1*s1*t1*y*y*y - n1*o1*u1*y*y*y + m1*r1*u1*y*y*y;

	return fabs(t_3/t_2);
//  return t_3/t_2;
}


