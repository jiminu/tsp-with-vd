#include "rg_IntsecImplicit5.h"
#include "rg_Polynomial.h"
#include "rg_RelativeOp.h"

// inserted by Jung-Hyun Ryu
#include "rg_IntersectFunc.h"

//#include <fstream>
#include <time.h>

rg_IntsecImplicit5::rg_IntsecImplicit5()
{

}

rg_IntsecImplicit5::~rg_IntsecImplicit5()
{

}

rg_dList<rg_Point2D> rg_IntsecImplicit5::intersectBzCurveVsBzCurve(const rg_BzCurve2D &curve_s, const rg_BzCurve2D &curve_t, rg_REAL &time)
{
	rg_dList<rg_Point2D> intersectPointList;
	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime = clock(); 
	//for(rg_INT i=0; i<20; i++)
	//{

	intersectPointList.removeAll();
    rg_ImplicitEquation  impEq   = implicitize(curve_t);
	//rg_ImplicitEquation  impEq   = rg_IntersectFunc::implicitizeNotUsingMathematica1(curve_t);//test

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
 
/*
rg_dList<rg_Point2D> rg_IntsecImplicit5::intersectBzCurveVsBzCurve(rg_BzCurve2D &curve_s, rg_BzCurve2D &curve_t)
{
	rg_dList<rg_Point2D> intersectPointList;
	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime = clock(); 
	for(rg_INT i=0; i<50; i++)
	{

//	intersectPointList.ClearAll();
    rg_ImplicitEquation  impEq   = rg_IntersectFunc::implicitize(curve_t);
    rg_Polynomial*       cs      = curve_s.convertBzCurve2Polynomial();
    rg_Polynomial        poly    = impEq.substitute(cs[0], cs[1]);
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

	}

	clock_t EndTime = clock();
	rg_REAL ComputingTime = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	ofstream fout("out.dat", ios::app);
	fout << " ---- Imp5 Alg. ----" << endl;
    fout << " degree of C1 : " << curve_s.getDegree();
    fout << ",  degree of C2 : " << curve_t.getDegree() << endl;
	fout << " # of intersection points  : " << intersectPointList.getSize() << endl;
	fout << " computation time  : " << ComputingTime << endl ;


 	return intersectPointList;
}  
*/

rg_ImplicitEquation rg_IntsecImplicit5::implicitize(const rg_BzCurve2D &curve)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL x3 = curve.getCtrlPt(3).getX();
	rg_REAL x4 = curve.getCtrlPt(4).getX();
    rg_REAL x5 = curve.getCtrlPt(5).getX();

	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL y3 = curve.getCtrlPt(3).getY();
	rg_REAL y4 = curve.getCtrlPt(4).getY();
    rg_REAL y5 = curve.getCtrlPt(5).getY();

//  If
//		    
//	x(t) = a5 t^5 + a4 t^4 + a3 t^3 + a2 t^2 + a1 t + a0
//
//	y(t) = b5 t^5 + b4 t^4 + b3 t^3 + b2 t^2 + b1 t + b0
//  
//  then the resultant of p(x, t) and q(y, t) is like the following.

    rg_REAL a5 = x5- 5*x4 + 10*x3 - 10*x2 +  5*x1 -    x0;
	rg_REAL a4 =     5*x4 - 20*x3 + 30*x2 - 20*x1 +  5*x0;
	rg_REAL a3 =            10*x3 - 30*x2 + 30*x1 - 10*x0;
	rg_REAL a2 =                    10*x2 - 20*x1 + 10*x0;
	rg_REAL a1 =                             5*x1 -  5*x0;
	rg_REAL a0 =                                       x0;

    rg_REAL b5 = y5- 5*y4 + 10*y3 - 10*y2 +  5*y1 -    y0;
	rg_REAL b4 =     5*y4 - 20*y3 + 30*y2 - 20*y1 +  5*y0;
	rg_REAL b3 =            10*y3 - 30*y2 + 30*y1 - 10*y0;
	rg_REAL b2 =                    10*y2 - 20*y1 + 10*y0;
	rg_REAL b1 =                             5*y1 -  5*y0;
	rg_REAL b0 =                                       y0;

  //			|														                                      |
  //			|	i0 x + i1 y + i2  j0 x + j1 y + j2	k0 x + k1 y + k2  l0 x + l1 y + l2  m0 x + m1 y + m2  |
  //			|														                                      |
  // R(P, Q) =	|	j0 x + j1 y + j2  n0 x + n1 y + n2	o0 x + o1 y + o2  p0 x + p1 y + p2  q0 x + q1 y + q2  |
  //			|														                                      |
  //			|	k0 x + k1 y + k2  o0 x + o1 y + o2	r0 x + r1 y + r2  s0 x + s1 y + s2  t0 x + t1 y + t2  |
  //			|														                                      |
  //            |   l0 x + l1 y + l2  p0 x + p1 y + p2  s0 x + s1 y + s2  u0 x + u1 y + u2  v0 x + v1 y + v2  |
  //            |                                                                                             |
  //            |   m0 x + m1 y + m2  q0 x + q1 y + q2  t0 x + t1 y + t2  v0 x + v1 y + v2  w0 x + w1 y + w2  |
   
	rg_REAL i0 =   0.0;
	rg_REAL i1 =   0.0;
	rg_REAL i2 =   a5*b4 - a4*b5;

	rg_REAL j0 =   0.0;
	rg_REAL j1 =   0.0;
	rg_REAL j2 =   a5*b3 - a3*b5;

	rg_REAL k0 =   0.0;
	rg_REAL k1 =   0.0;
	rg_REAL k2 =   a5*b2 - a2*b5;

	rg_REAL l0 =   0.0;
	rg_REAL l1 =   0.0;
	rg_REAL l2 =   a5*b1 - a1*b5;

	rg_REAL m0 =   b5;
	rg_REAL m1 = - a5;
	rg_REAL m2 =   a5*b0 - a0*b5;

	rg_REAL n0 =   0.0;
	rg_REAL n1 =   0.0;
	rg_REAL n2 =   a5*b2 - a2*b5 + a4*b3 - a3*b4;

	rg_REAL o0 =   0.0;
	rg_REAL o1 =   0.0;
	rg_REAL o2 =   a5*b1 - a1*b5 + a4*b2 - a2*b4;

    rg_REAL p0 =   b5;
	rg_REAL p1 = - a5;
	rg_REAL p2 =   a5*b0 - a0*b5 + a4*b1 - a1*b4;

  	rg_REAL q0 =   b4;
	rg_REAL q1 = - a4;
	rg_REAL q2 =   a4*b0 - a0*b4;

  	rg_REAL r0 =   b5;
	rg_REAL r1 = - a5;
	rg_REAL r2 =   a5*b0 - a0*b5 + a4*b1 - a1*b4 + a3*b2 - a2*b3;

  	rg_REAL s0 =   b4;
	rg_REAL s1 = - a4;
	rg_REAL s2 =   a4*b0 - a0*b4 + a3*b1 - a1*b3;

   	rg_REAL t0 =   b3;
	rg_REAL t1 = - a3;
	rg_REAL t2 =   a3*b0 - a0*b3;

   	rg_REAL u0 =   b3;
	rg_REAL u1 = - a3;
	rg_REAL u2 =   a3*b0 - a0*b3 + a2*b1 - a1*b2;

	rg_REAL v0 =   b2;
	rg_REAL v1 = - a2;
	rg_REAL v2 =   a2*b0 - a0*b2;

   	rg_REAL w0 =   b1;
	rg_REAL w1 = - a1;
	rg_REAL w2 =   a1*b0 - a0*b1;

// If the form of the implicitization is 
//    k0 
//  + k1  x       + k2      y
//  + k3  x^2     + k4  x   y   + k5      y^2
//  + k6  x^3     + k7  x^2 y   + k8  x   y^2  + k9      y^3
//  + k10 x^4     + k11 x^3 y   + k12 x^2 y^2  + k13 x   y^3 + k14     y^4
//  + k15 x^5     + k16 x^4 y   + k17 x^3 y^2  + k18 x^2 y^3 + k19 x   y^4 + k20 y^5

    rg_REAL *k = new rg_REAL[21];
    
    k[0]  = m2*m2*p2*p2*r2 - 2*l2*m2*p2*q2*r2 + l2*l2*q2*q2*r2 - 2*m2*m2*o2*p2*s2 + 
  2*l2*m2*o2*q2*s2 + 2*k2*m2*p2*q2*s2 - 2*k2*l2*q2*q2*s2 + m2*m2*n2*s2*s2 - 
  2*j2*m2*q2*s2*s2 + i2*q2*q2*s2*s2 + 2*l2*m2*o2*p2*t2 - 2*k2*m2*p2*p2*t2 - 
  2*l2*l2*o2*q2*t2 + 2*k2*l2*p2*q2*t2 - 2*l2*m2*n2*s2*t2 + 2*j2*m2*p2*s2*t2 + 
  2*j2*l2*q2*s2*t2 - 2*i2*p2*q2*s2*t2 + l2*l2*n2*t2*t2 - 2*j2*l2*p2*t2*t2 + 
  i2*p2*p2*t2*t2 + m2*m2*o2*o2*u2 - 2*k2*m2*o2*q2*u2 + k2*k2*q2*q2*u2 - 
  m2*m2*n2*r2*u2 + 2*j2*m2*q2*r2*u2 - i2*q2*q2*r2*u2 + 2*k2*m2*n2*t2*u2 - 
  2*j2*m2*o2*t2*u2 - 2*j2*k2*q2*t2*u2 + 2*i2*o2*q2*t2*u2 + j2*j2*t2*t2*u2 - 
  i2*n2*t2*t2*u2 - 2*l2*m2*o2*o2*v2 + 2*k2*m2*o2*p2*v2 + 2*k2*l2*o2*q2*v2 - 
  2*k2*k2*p2*q2*v2 + 2*l2*m2*n2*r2*v2 - 2*j2*m2*p2*r2*v2 - 2*j2*l2*q2*r2*v2 + 
  2*i2*p2*q2*r2*v2 - 2*k2*m2*n2*s2*v2 + 2*j2*m2*o2*s2*v2 + 2*j2*k2*q2*s2*v2 - 
  2*i2*o2*q2*s2*v2 - 2*k2*l2*n2*t2*v2 + 2*j2*l2*o2*t2*v2 + 2*j2*k2*p2*t2*v2 - 
  2*i2*o2*p2*t2*v2 - 2*j2*j2*s2*t2*v2 + 2*i2*n2*s2*t2*v2 + k2*k2*n2*v2*v2 - 
  2*j2*k2*o2*v2*v2 + i2*o2*o2*v2*v2 + j2*j2*r2*v2*v2 - i2*n2*r2*v2*v2 + 
  l2*l2*o2*o2*w2 - 2*k2*l2*o2*p2*w2 + k2*k2*p2*p2*w2 - l2*l2*n2*r2*w2 + 
  2*j2*l2*p2*r2*w2 - i2*p2*p2*r2*w2 + 2*k2*l2*n2*s2*w2 - 2*j2*l2*o2*s2*w2 - 
  2*j2*k2*p2*s2*w2 + 2*i2*o2*p2*s2*w2 + j2*j2*s2*s2*w2 - i2*n2*s2*s2*w2 - 
  k2*k2*n2*u2*w2 + 2*j2*k2*o2*u2*w2 - i2*o2*o2*u2*w2 - j2*j2*r2*u2*w2 + 
  i2*n2*r2*u2*w2;
    
    k[1]  = m2*m2*p2*p2*r0 - 2*l2*m2*p2*q2*r0 + l2*l2*q2*q2*r0 + 
  2*m2*m2*p0*p2*r2 + 2*m0*m2*p2*p2*r2 - 2*l2*m2*p2*q0*r2 - 
  2*l2*m2*p0*q2*r2 - 2*l2*m0*p2*q2*r2 - 2*l0*m2*p2*q2*r2 + 
  2*l2*l2*q0*q2*r2 + 2*l0*l2*q2*q2*r2 - 2*m2*m2*o2*p2*s0 + 
  2*l2*m2*o2*q2*s0 + 2*k2*m2*p2*q2*s0 - 2*k2*l2*q2*q2*s0 - 
  2*m2*m2*o2*p0*s2 - 2*m2*m2*o0*p2*s2 - 4*m0*m2*o2*p2*s2 + 
  2*l2*m2*o2*q0*s2 + 2*k2*m2*p2*q0*s2 + 2*l2*m2*o0*q2*s2 + 
  2*l2*m0*o2*q2*s2 + 2*l0*m2*o2*q2*s2 + 2*k2*m2*p0*q2*s2 + 
  2*k2*m0*p2*q2*s2 + 2*k0*m2*p2*q2*s2 - 4*k2*l2*q0*q2*s2 - 
  2*k2*l0*q2*q2*s2 - 2*k0*l2*q2*q2*s2 + 2*m2*m2*n2*s0*s2 - 
  4*j2*m2*q2*s0*s2 + 2*i2*q2*q2*s0*s2 + m2*m2*n0*s2*s2 + 
  2*m0*m2*n2*s2*s2 - 2*j2*m2*q0*s2*s2 - 2*j2*m0*q2*s2*s2 - 
  2*j0*m2*q2*s2*s2 + 2*i2*q0*q2*s2*s2 + i0*q2*q2*s2*s2 + 
  2*l2*m2*o2*p2*t0 - 2*k2*m2*p2*p2*t0 - 2*l2*l2*o2*q2*t0 + 
  2*k2*l2*p2*q2*t0 - 2*l2*m2*n2*s2*t0 + 2*j2*m2*p2*s2*t0 + 
  2*j2*l2*q2*s2*t0 - 2*i2*p2*q2*s2*t0 + 2*l2*m2*o2*p0*t2 + 
  2*l2*m2*o0*p2*t2 + 2*l2*m0*o2*p2*t2 + 2*l0*m2*o2*p2*t2 - 
  4*k2*m2*p0*p2*t2 - 2*k2*m0*p2*p2*t2 - 2*k0*m2*p2*p2*t2 - 
  2*l2*l2*o2*q0*t2 + 2*k2*l2*p2*q0*t2 - 2*l2*l2*o0*q2*t2 - 
  4*l0*l2*o2*q2*t2 + 2*k2*l2*p0*q2*t2 + 2*k2*l0*p2*q2*t2 + 
  2*k0*l2*p2*q2*t2 - 2*l2*m2*n2*s0*t2 + 2*j2*m2*p2*s0*t2 + 
  2*j2*l2*q2*s0*t2 - 2*i2*p2*q2*s0*t2 - 2*l2*m2*n0*s2*t2 - 
  2*l2*m0*n2*s2*t2 - 2*l0*m2*n2*s2*t2 + 2*j2*m2*p0*s2*t2 + 
  2*j2*m0*p2*s2*t2 + 2*j0*m2*p2*s2*t2 + 2*j2*l2*q0*s2*t2 - 
  2*i2*p2*q0*s2*t2 + 2*j2*l0*q2*s2*t2 + 2*j0*l2*q2*s2*t2 - 
  2*i2*p0*q2*s2*t2 - 2*i0*p2*q2*s2*t2 + 2*l2*l2*n2*t0*t2 - 
  4*j2*l2*p2*t0*t2 + 2*i2*p2*p2*t0*t2 + l2*l2*n0*t2*t2 + 
  2*l0*l2*n2*t2*t2 - 2*j2*l2*p0*t2*t2 - 2*j2*l0*p2*t2*t2 - 
  2*j0*l2*p2*t2*t2 + 2*i2*p0*p2*t2*t2 + i0*p2*p2*t2*t2 + m2*m2*o2*o2*u0 - 
  2*k2*m2*o2*q2*u0 + k2*k2*q2*q2*u0 - m2*m2*n2*r2*u0 + 
  2*j2*m2*q2*r2*u0 - i2*q2*q2*r2*u0 + 2*k2*m2*n2*t2*u0 - 
  2*j2*m2*o2*t2*u0 - 2*j2*k2*q2*t2*u0 + 2*i2*o2*q2*t2*u0 + 
  j2*j2*t2*t2*u0 - i2*n2*t2*t2*u0 + 2*m2*m2*o0*o2*u2 + 2*m0*m2*o2*o2*u2 - 
  2*k2*m2*o2*q0*u2 - 2*k2*m2*o0*q2*u2 - 2*k2*m0*o2*q2*u2 - 
  2*k0*m2*o2*q2*u2 + 2*k2*k2*q0*q2*u2 + 2*k0*k2*q2*q2*u2 - 
  m2*m2*n2*r0*u2 + 2*j2*m2*q2*r0*u2 - i2*q2*q2*r0*u2 - m2*m2*n0*r2*u2 - 
  2*m0*m2*n2*r2*u2 + 2*j2*m2*q0*r2*u2 + 2*j2*m0*q2*r2*u2 + 
  2*j0*m2*q2*r2*u2 - 2*i2*q0*q2*r2*u2 - i0*q2*q2*r2*u2 + 
  2*k2*m2*n2*t0*u2 - 2*j2*m2*o2*t0*u2 - 2*j2*k2*q2*t0*u2 + 
  2*i2*o2*q2*t0*u2 + 2*k2*m2*n0*t2*u2 + 2*k2*m0*n2*t2*u2 + 
  2*k0*m2*n2*t2*u2 - 2*j2*m2*o0*t2*u2 - 2*j2*m0*o2*t2*u2 - 
  2*j0*m2*o2*t2*u2 - 2*j2*k2*q0*t2*u2 + 2*i2*o2*q0*t2*u2 - 
  2*j2*k0*q2*t2*u2 - 2*j0*k2*q2*t2*u2 + 2*i2*o0*q2*t2*u2 + 
  2*i0*o2*q2*t2*u2 + 2*j2*j2*t0*t2*u2 - 2*i2*n2*t0*t2*u2 + 
  2*j0*j2*t2*t2*u2 - i2*n0*t2*t2*u2 - i0*n2*t2*t2*u2 - 2*l2*m2*o2*o2*v0 + 
  2*k2*m2*o2*p2*v0 + 2*k2*l2*o2*q2*v0 - 2*k2*k2*p2*q2*v0 + 
  2*l2*m2*n2*r2*v0 - 2*j2*m2*p2*r2*v0 - 2*j2*l2*q2*r2*v0 + 
  2*i2*p2*q2*r2*v0 - 2*k2*m2*n2*s2*v0 + 2*j2*m2*o2*s2*v0 + 
  2*j2*k2*q2*s2*v0 - 2*i2*o2*q2*s2*v0 - 2*k2*l2*n2*t2*v0 + 
  2*j2*l2*o2*t2*v0 + 2*j2*k2*p2*t2*v0 - 2*i2*o2*p2*t2*v0 - 
  2*j2*j2*s2*t2*v0 + 2*i2*n2*s2*t2*v0 - 4*l2*m2*o0*o2*v2 - 
  2*l2*m0*o2*o2*v2 - 2*l0*m2*o2*o2*v2 + 2*k2*m2*o2*p0*v2 + 
  2*k2*m2*o0*p2*v2 + 2*k2*m0*o2*p2*v2 + 2*k0*m2*o2*p2*v2 + 
  2*k2*l2*o2*q0*v2 - 2*k2*k2*p2*q0*v2 + 2*k2*l2*o0*q2*v2 + 
  2*k2*l0*o2*q2*v2 + 2*k0*l2*o2*q2*v2 - 2*k2*k2*p0*q2*v2 - 
  4*k0*k2*p2*q2*v2 + 2*l2*m2*n2*r0*v2 - 2*j2*m2*p2*r0*v2 - 
  2*j2*l2*q2*r0*v2 + 2*i2*p2*q2*r0*v2 + 2*l2*m2*n0*r2*v2 + 
  2*l2*m0*n2*r2*v2 + 2*l0*m2*n2*r2*v2 - 2*j2*m2*p0*r2*v2 - 
  2*j2*m0*p2*r2*v2 - 2*j0*m2*p2*r2*v2 - 2*j2*l2*q0*r2*v2 + 
  2*i2*p2*q0*r2*v2 - 2*j2*l0*q2*r2*v2 - 2*j0*l2*q2*r2*v2 + 
  2*i2*p0*q2*r2*v2 + 2*i0*p2*q2*r2*v2 - 2*k2*m2*n2*s0*v2 + 
  2*j2*m2*o2*s0*v2 + 2*j2*k2*q2*s0*v2 - 2*i2*o2*q2*s0*v2 - 
  2*k2*m2*n0*s2*v2 - 2*k2*m0*n2*s2*v2 - 2*k0*m2*n2*s2*v2 + 
  2*j2*m2*o0*s2*v2 + 2*j2*m0*o2*s2*v2 + 2*j0*m2*o2*s2*v2 + 
  2*j2*k2*q0*s2*v2 - 2*i2*o2*q0*s2*v2 + 2*j2*k0*q2*s2*v2 + 
  2*j0*k2*q2*s2*v2 - 2*i2*o0*q2*s2*v2 - 2*i0*o2*q2*s2*v2 - 
  2*k2*l2*n2*t0*v2 + 2*j2*l2*o2*t0*v2 + 2*j2*k2*p2*t0*v2 - 
  2*i2*o2*p2*t0*v2 - 2*j2*j2*s2*t0*v2 + 2*i2*n2*s2*t0*v2 - 
  2*k2*l2*n0*t2*v2 - 2*k2*l0*n2*t2*v2 - 2*k0*l2*n2*t2*v2 + 
  2*j2*l2*o0*t2*v2 + 2*j2*l0*o2*t2*v2 + 2*j0*l2*o2*t2*v2 + 
  2*j2*k2*p0*t2*v2 - 2*i2*o2*p0*t2*v2 + 2*j2*k0*p2*t2*v2 + 
  2*j0*k2*p2*t2*v2 - 2*i2*o0*p2*t2*v2 - 2*i0*o2*p2*t2*v2 - 
  2*j2*j2*s0*t2*v2 + 2*i2*n2*s0*t2*v2 - 4*j0*j2*s2*t2*v2 + 
  2*i2*n0*s2*t2*v2 + 2*i0*n2*s2*t2*v2 + 2*k2*k2*n2*v0*v2 - 
  4*j2*k2*o2*v0*v2 + 2*i2*o2*o2*v0*v2 + 2*j2*j2*r2*v0*v2 - 
  2*i2*n2*r2*v0*v2 + k2*k2*n0*v2*v2 + 2*k0*k2*n2*v2*v2 - 
  2*j2*k2*o0*v2*v2 - 2*j2*k0*o2*v2*v2 - 2*j0*k2*o2*v2*v2 + 
  2*i2*o0*o2*v2*v2 + i0*o2*o2*v2*v2 + j2*j2*r0*v2*v2 - i2*n2*r0*v2*v2 + 
  2*j0*j2*r2*v2*v2 - i2*n0*r2*v2*v2 - i0*n2*r2*v2*v2 + l2*l2*o2*o2*w0 - 
  2*k2*l2*o2*p2*w0 + k2*k2*p2*p2*w0 - l2*l2*n2*r2*w0 + 
  2*j2*l2*p2*r2*w0 - i2*p2*p2*r2*w0 + 2*k2*l2*n2*s2*w0 - 
  2*j2*l2*o2*s2*w0 - 2*j2*k2*p2*s2*w0 + 2*i2*o2*p2*s2*w0 + 
  j2*j2*s2*s2*w0 - i2*n2*s2*s2*w0 - k2*k2*n2*u2*w0 + 2*j2*k2*o2*u2*w0 - 
  i2*o2*o2*u2*w0 - j2*j2*r2*u2*w0 + i2*n2*r2*u2*w0 + 2*l2*l2*o0*o2*w2 + 
  2*l0*l2*o2*o2*w2 - 2*k2*l2*o2*p0*w2 - 2*k2*l2*o0*p2*w2 - 
  2*k2*l0*o2*p2*w2 - 2*k0*l2*o2*p2*w2 + 2*k2*k2*p0*p2*w2 + 
  2*k0*k2*p2*p2*w2 - l2*l2*n2*r0*w2 + 2*j2*l2*p2*r0*w2 - 
  i2*p2*p2*r0*w2 - l2*l2*n0*r2*w2 - 2*l0*l2*n2*r2*w2 + 
  2*j2*l2*p0*r2*w2 + 2*j2*l0*p2*r2*w2 + 2*j0*l2*p2*r2*w2 - 
  2*i2*p0*p2*r2*w2 - i0*p2*p2*r2*w2 + 2*k2*l2*n2*s0*w2 - 
  2*j2*l2*o2*s0*w2 - 2*j2*k2*p2*s0*w2 + 2*i2*o2*p2*s0*w2 + 
  2*k2*l2*n0*s2*w2 + 2*k2*l0*n2*s2*w2 + 2*k0*l2*n2*s2*w2 - 
  2*j2*l2*o0*s2*w2 - 2*j2*l0*o2*s2*w2 - 2*j0*l2*o2*s2*w2 - 
  2*j2*k2*p0*s2*w2 + 2*i2*o2*p0*s2*w2 - 2*j2*k0*p2*s2*w2 - 
  2*j0*k2*p2*s2*w2 + 2*i2*o0*p2*s2*w2 + 2*i0*o2*p2*s2*w2 + 
  2*j2*j2*s0*s2*w2 - 2*i2*n2*s0*s2*w2 + 2*j0*j2*s2*s2*w2 - 
  i2*n0*s2*s2*w2 - i0*n2*s2*s2*w2 - k2*k2*n2*u0*w2 + 2*j2*k2*o2*u0*w2 - 
  i2*o2*o2*u0*w2 - j2*j2*r2*u0*w2 + i2*n2*r2*u0*w2 - k2*k2*n0*u2*w2 - 
  2*k0*k2*n2*u2*w2 + 2*j2*k2*o0*u2*w2 + 2*j2*k0*o2*u2*w2 + 
  2*j0*k2*o2*u2*w2 - 2*i2*o0*o2*u2*w2 - i0*o2*o2*u2*w2 - 
  j2*j2*r0*u2*w2 + i2*n2*r0*u2*w2 - 2*j0*j2*r2*u2*w2 + 
  i2*n0*r2*u2*w2 + i0*n2*r2*u2*w2;
    
    k[2]  = m2*m2*p2*p2*r1 - 
  2*l2*m2*p2*q2*r1 + l2*l2*q2*q2*r1 + 2*m2*m2*p1*p2*r2 + 
  2*m1*m2*p2*p2*r2 - 2*l2*m2*p2*q1*r2 - 2*l2*m2*p1*q2*r2 - 
  2*l2*m1*p2*q2*r2 - 2*l1*m2*p2*q2*r2 + 2*l2*l2*q1*q2*r2 + 
  2*l1*l2*q2*q2*r2 - 2*m2*m2*o2*p2*s1 + 2*l2*m2*o2*q2*s1 + 
  2*k2*m2*p2*q2*s1 - 2*k2*l2*q2*q2*s1 - 2*m2*m2*o2*p1*s2 - 
  2*m2*m2*o1*p2*s2 - 4*m1*m2*o2*p2*s2 + 2*l2*m2*o2*q1*s2 + 
  2*k2*m2*p2*q1*s2 + 2*l2*m2*o1*q2*s2 + 2*l2*m1*o2*q2*s2 + 
  2*l1*m2*o2*q2*s2 + 2*k2*m2*p1*q2*s2 + 2*k2*m1*p2*q2*s2 + 
  2*k1*m2*p2*q2*s2 - 4*k2*l2*q1*q2*s2 - 2*k2*l1*q2*q2*s2 - 
  2*k1*l2*q2*q2*s2 + 2*m2*m2*n2*s1*s2 - 4*j2*m2*q2*s1*s2 + 
  2*i2*q2*q2*s1*s2 + m2*m2*n1*s2*s2 + 2*m1*m2*n2*s2*s2 - 
  2*j2*m2*q1*s2*s2 - 2*j2*m1*q2*s2*s2 - 2*j1*m2*q2*s2*s2 + 
  2*i2*q1*q2*s2*s2 + i1*q2*q2*s2*s2 + 2*l2*m2*o2*p2*t1 - 
  2*k2*m2*p2*p2*t1 - 2*l2*l2*o2*q2*t1 + 2*k2*l2*p2*q2*t1 - 
  2*l2*m2*n2*s2*t1 + 2*j2*m2*p2*s2*t1 + 2*j2*l2*q2*s2*t1 - 
  2*i2*p2*q2*s2*t1 + 2*l2*m2*o2*p1*t2 + 2*l2*m2*o1*p2*t2 + 
  2*l2*m1*o2*p2*t2 + 2*l1*m2*o2*p2*t2 - 4*k2*m2*p1*p2*t2 - 
  2*k2*m1*p2*p2*t2 - 2*k1*m2*p2*p2*t2 - 2*l2*l2*o2*q1*t2 + 
  2*k2*l2*p2*q1*t2 - 2*l2*l2*o1*q2*t2 - 4*l1*l2*o2*q2*t2 + 
  2*k2*l2*p1*q2*t2 + 2*k2*l1*p2*q2*t2 + 2*k1*l2*p2*q2*t2 - 
  2*l2*m2*n2*s1*t2 + 2*j2*m2*p2*s1*t2 + 2*j2*l2*q2*s1*t2 - 
  2*i2*p2*q2*s1*t2 - 2*l2*m2*n1*s2*t2 - 2*l2*m1*n2*s2*t2 - 
  2*l1*m2*n2*s2*t2 + 2*j2*m2*p1*s2*t2 + 2*j2*m1*p2*s2*t2 + 
  2*j1*m2*p2*s2*t2 + 2*j2*l2*q1*s2*t2 - 2*i2*p2*q1*s2*t2 + 
  2*j2*l1*q2*s2*t2 + 2*j1*l2*q2*s2*t2 - 2*i2*p1*q2*s2*t2 - 
  2*i1*p2*q2*s2*t2 + 2*l2*l2*n2*t1*t2 - 4*j2*l2*p2*t1*t2 + 
  2*i2*p2*p2*t1*t2 + l2*l2*n1*t2*t2 + 2*l1*l2*n2*t2*t2 - 
  2*j2*l2*p1*t2*t2 - 2*j2*l1*p2*t2*t2 - 2*j1*l2*p2*t2*t2 + 
  2*i2*p1*p2*t2*t2 + i1*p2*p2*t2*t2 + m2*m2*o2*o2*u1 - 2*k2*m2*o2*q2*u1 + 
  k2*k2*q2*q2*u1 - m2*m2*n2*r2*u1 + 2*j2*m2*q2*r2*u1 - i2*q2*q2*r2*u1 + 
  2*k2*m2*n2*t2*u1 - 2*j2*m2*o2*t2*u1 - 2*j2*k2*q2*t2*u1 + 
  2*i2*o2*q2*t2*u1 + j2*j2*t2*t2*u1 - i2*n2*t2*t2*u1 + 2*m2*m2*o1*o2*u2 + 
  2*m1*m2*o2*o2*u2 - 2*k2*m2*o2*q1*u2 - 2*k2*m2*o1*q2*u2 - 
  2*k2*m1*o2*q2*u2 - 2*k1*m2*o2*q2*u2 + 2*k2*k2*q1*q2*u2 + 
  2*k1*k2*q2*q2*u2 - m2*m2*n2*r1*u2 + 2*j2*m2*q2*r1*u2 - 
  i2*q2*q2*r1*u2 - m2*m2*n1*r2*u2 - 2*m1*m2*n2*r2*u2 + 
  2*j2*m2*q1*r2*u2 + 2*j2*m1*q2*r2*u2 + 2*j1*m2*q2*r2*u2 - 
  2*i2*q1*q2*r2*u2 - i1*q2*q2*r2*u2 + 2*k2*m2*n2*t1*u2 - 
  2*j2*m2*o2*t1*u2 - 2*j2*k2*q2*t1*u2 + 2*i2*o2*q2*t1*u2 + 
  2*k2*m2*n1*t2*u2 + 2*k2*m1*n2*t2*u2 + 2*k1*m2*n2*t2*u2 - 
  2*j2*m2*o1*t2*u2 - 2*j2*m1*o2*t2*u2 - 2*j1*m2*o2*t2*u2 - 
  2*j2*k2*q1*t2*u2 + 2*i2*o2*q1*t2*u2 - 2*j2*k1*q2*t2*u2 - 
  2*j1*k2*q2*t2*u2 + 2*i2*o1*q2*t2*u2 + 2*i1*o2*q2*t2*u2 + 
  2*j2*j2*t1*t2*u2 - 2*i2*n2*t1*t2*u2 + 2*j1*j2*t2*t2*u2 - 
  i2*n1*t2*t2*u2 - i1*n2*t2*t2*u2 - 2*l2*m2*o2*o2*v1 + 
  2*k2*m2*o2*p2*v1 + 2*k2*l2*o2*q2*v1 - 2*k2*k2*p2*q2*v1 + 
  2*l2*m2*n2*r2*v1 - 2*j2*m2*p2*r2*v1 - 2*j2*l2*q2*r2*v1 + 
  2*i2*p2*q2*r2*v1 - 2*k2*m2*n2*s2*v1 + 2*j2*m2*o2*s2*v1 + 
  2*j2*k2*q2*s2*v1 - 2*i2*o2*q2*s2*v1 - 2*k2*l2*n2*t2*v1 + 
  2*j2*l2*o2*t2*v1 + 2*j2*k2*p2*t2*v1 - 2*i2*o2*p2*t2*v1 - 
  2*j2*j2*s2*t2*v1 + 2*i2*n2*s2*t2*v1 - 4*l2*m2*o1*o2*v2 - 
  2*l2*m1*o2*o2*v2 - 2*l1*m2*o2*o2*v2 + 2*k2*m2*o2*p1*v2 + 
  2*k2*m2*o1*p2*v2 + 2*k2*m1*o2*p2*v2 + 2*k1*m2*o2*p2*v2 + 
  2*k2*l2*o2*q1*v2 - 2*k2*k2*p2*q1*v2 + 2*k2*l2*o1*q2*v2 + 
  2*k2*l1*o2*q2*v2 + 2*k1*l2*o2*q2*v2 - 2*k2*k2*p1*q2*v2 - 
  4*k1*k2*p2*q2*v2 + 2*l2*m2*n2*r1*v2 - 2*j2*m2*p2*r1*v2 - 
  2*j2*l2*q2*r1*v2 + 2*i2*p2*q2*r1*v2 + 2*l2*m2*n1*r2*v2 + 
  2*l2*m1*n2*r2*v2 + 2*l1*m2*n2*r2*v2 - 2*j2*m2*p1*r2*v2 - 
  2*j2*m1*p2*r2*v2 - 2*j1*m2*p2*r2*v2 - 2*j2*l2*q1*r2*v2 + 
  2*i2*p2*q1*r2*v2 - 2*j2*l1*q2*r2*v2 - 2*j1*l2*q2*r2*v2 + 
  2*i2*p1*q2*r2*v2 + 2*i1*p2*q2*r2*v2 - 2*k2*m2*n2*s1*v2 + 
  2*j2*m2*o2*s1*v2 + 2*j2*k2*q2*s1*v2 - 2*i2*o2*q2*s1*v2 - 
  2*k2*m2*n1*s2*v2 - 2*k2*m1*n2*s2*v2 - 2*k1*m2*n2*s2*v2 + 
  2*j2*m2*o1*s2*v2 + 2*j2*m1*o2*s2*v2 + 2*j1*m2*o2*s2*v2 + 
  2*j2*k2*q1*s2*v2 - 2*i2*o2*q1*s2*v2 + 2*j2*k1*q2*s2*v2 + 
  2*j1*k2*q2*s2*v2 - 2*i2*o1*q2*s2*v2 - 2*i1*o2*q2*s2*v2 - 
  2*k2*l2*n2*t1*v2 + 2*j2*l2*o2*t1*v2 + 2*j2*k2*p2*t1*v2 - 
  2*i2*o2*p2*t1*v2 - 2*j2*j2*s2*t1*v2 + 2*i2*n2*s2*t1*v2 - 
  2*k2*l2*n1*t2*v2 - 2*k2*l1*n2*t2*v2 - 2*k1*l2*n2*t2*v2 + 
  2*j2*l2*o1*t2*v2 + 2*j2*l1*o2*t2*v2 + 2*j1*l2*o2*t2*v2 + 
  2*j2*k2*p1*t2*v2 - 2*i2*o2*p1*t2*v2 + 2*j2*k1*p2*t2*v2 + 
  2*j1*k2*p2*t2*v2 - 2*i2*o1*p2*t2*v2 - 2*i1*o2*p2*t2*v2 - 
  2*j2*j2*s1*t2*v2 + 2*i2*n2*s1*t2*v2 - 4*j1*j2*s2*t2*v2 + 
  2*i2*n1*s2*t2*v2 + 2*i1*n2*s2*t2*v2 + 2*k2*k2*n2*v1*v2 - 
  4*j2*k2*o2*v1*v2 + 2*i2*o2*o2*v1*v2 + 2*j2*j2*r2*v1*v2 - 
  2*i2*n2*r2*v1*v2 + k2*k2*n1*v2*v2 + 2*k1*k2*n2*v2*v2 - 
  2*j2*k2*o1*v2*v2 - 2*j2*k1*o2*v2*v2 - 2*j1*k2*o2*v2*v2 + 
  2*i2*o1*o2*v2*v2 + i1*o2*o2*v2*v2 + j2*j2*r1*v2*v2 - i2*n2*r1*v2*v2 + 
  2*j1*j2*r2*v2*v2 - i2*n1*r2*v2*v2 - i1*n2*r2*v2*v2 + l2*l2*o2*o2*w1 - 
  2*k2*l2*o2*p2*w1 + k2*k2*p2*p2*w1 - l2*l2*n2*r2*w1 + 
  2*j2*l2*p2*r2*w1 - i2*p2*p2*r2*w1 + 2*k2*l2*n2*s2*w1 - 
  2*j2*l2*o2*s2*w1 - 2*j2*k2*p2*s2*w1 + 2*i2*o2*p2*s2*w1 + 
  j2*j2*s2*s2*w1 - i2*n2*s2*s2*w1 - k2*k2*n2*u2*w1 + 2*j2*k2*o2*u2*w1 - 
  i2*o2*o2*u2*w1 - j2*j2*r2*u2*w1 + i2*n2*r2*u2*w1 + 2*l2*l2*o1*o2*w2 + 
  2*l1*l2*o2*o2*w2 - 2*k2*l2*o2*p1*w2 - 2*k2*l2*o1*p2*w2 - 
  2*k2*l1*o2*p2*w2 - 2*k1*l2*o2*p2*w2 + 2*k2*k2*p1*p2*w2 + 
  2*k1*k2*p2*p2*w2 - l2*l2*n2*r1*w2 + 2*j2*l2*p2*r1*w2 - 
  i2*p2*p2*r1*w2 - l2*l2*n1*r2*w2 - 2*l1*l2*n2*r2*w2 + 
  2*j2*l2*p1*r2*w2 + 2*j2*l1*p2*r2*w2 + 2*j1*l2*p2*r2*w2 - 
  2*i2*p1*p2*r2*w2 - i1*p2*p2*r2*w2 + 2*k2*l2*n2*s1*w2 - 
  2*j2*l2*o2*s1*w2 - 2*j2*k2*p2*s1*w2 + 2*i2*o2*p2*s1*w2 + 
  2*k2*l2*n1*s2*w2 + 2*k2*l1*n2*s2*w2 + 2*k1*l2*n2*s2*w2 - 
  2*j2*l2*o1*s2*w2 - 2*j2*l1*o2*s2*w2 - 2*j1*l2*o2*s2*w2 - 
  2*j2*k2*p1*s2*w2 + 2*i2*o2*p1*s2*w2 - 2*j2*k1*p2*s2*w2 - 
  2*j1*k2*p2*s2*w2 + 2*i2*o1*p2*s2*w2 + 2*i1*o2*p2*s2*w2 + 
  2*j2*j2*s1*s2*w2 - 2*i2*n2*s1*s2*w2 + 2*j1*j2*s2*s2*w2 - 
  i2*n1*s2*s2*w2 - i1*n2*s2*s2*w2 - k2*k2*n2*u1*w2 + 2*j2*k2*o2*u1*w2 - 
  i2*o2*o2*u1*w2 - j2*j2*r2*u1*w2 + i2*n2*r2*u1*w2 - k2*k2*n1*u2*w2 - 
  2*k1*k2*n2*u2*w2 + 2*j2*k2*o1*u2*w2 + 2*j2*k1*o2*u2*w2 + 
  2*j1*k2*o2*u2*w2 - 2*i2*o1*o2*u2*w2 - i1*o2*o2*u2*w2 - 
  j2*j2*r1*u2*w2 + i2*n2*r1*u2*w2 - 2*j1*j2*r2*u2*w2 + 
  i2*n1*r2*u2*w2 + i1*n2*r2*u2*w2;
    
    k[3]  = 2*m2*m2*p0*p2*r0 + 
  2*m0*m2*p2*p2*r0 - 2*l2*m2*p2*q0*r0 - 2*l2*m2*p0*q2*r0 - 
  2*l2*m0*p2*q2*r0 - 2*l0*m2*p2*q2*r0 + 2*l2*l2*q0*q2*r0 + 
  2*l0*l2*q2*q2*r0 + m2*m2*p0*p0*r2 + 4*m0*m2*p0*p2*r2 + 
  m0*m0*p2*p2*r2 - 2*l2*m2*p0*q0*r2 - 2*l2*m0*p2*q0*r2 - 
  2*l0*m2*p2*q0*r2 + l2*l2*q0*q0*r2 - 2*l2*m0*p0*q2*r2 - 
  2*l0*m2*p0*q2*r2 - 2*l0*m0*p2*q2*r2 + 4*l0*l2*q0*q2*r2 + 
  l0*l0*q2*q2*r2 - 2*m2*m2*o2*p0*s0 - 2*m2*m2*o0*p2*s0 - 
  4*m0*m2*o2*p2*s0 + 2*l2*m2*o2*q0*s0 + 2*k2*m2*p2*q0*s0 + 
  2*l2*m2*o0*q2*s0 + 2*l2*m0*o2*q2*s0 + 2*l0*m2*o2*q2*s0 + 
  2*k2*m2*p0*q2*s0 + 2*k2*m0*p2*q2*s0 + 2*k0*m2*p2*q2*s0 - 
  4*k2*l2*q0*q2*s0 - 2*k2*l0*q2*q2*s0 - 2*k0*l2*q2*q2*s0 + 
  m2*m2*n2*s0*s0 - 2*j2*m2*q2*s0*s0 + i2*q2*q2*s0*s0 - 
  2*m2*m2*o0*p0*s2 - 4*m0*m2*o2*p0*s2 - 4*m0*m2*o0*p2*s2 - 
  2*m0*m0*o2*p2*s2 + 2*l2*m2*o0*q0*s2 + 2*l2*m0*o2*q0*s2 + 
  2*l0*m2*o2*q0*s2 + 2*k2*m2*p0*q0*s2 + 2*k2*m0*p2*q0*s2 + 
  2*k0*m2*p2*q0*s2 - 2*k2*l2*q0*q0*s2 + 2*l2*m0*o0*q2*s2 + 
  2*l0*m2*o0*q2*s2 + 2*l0*m0*o2*q2*s2 + 2*k2*m0*p0*q2*s2 + 
  2*k0*m2*p0*q2*s2 + 2*k0*m0*p2*q2*s2 - 4*k2*l0*q0*q2*s2 - 
  4*k0*l2*q0*q2*s2 - 2*k0*l0*q2*q2*s2 + 2*m2*m2*n0*s0*s2 + 
  4*m0*m2*n2*s0*s2 - 4*j2*m2*q0*s0*s2 - 4*j2*m0*q2*s0*s2 - 
  4*j0*m2*q2*s0*s2 + 4*i2*q0*q2*s0*s2 + 2*i0*q2*q2*s0*s2 + 
  2*m0*m2*n0*s2*s2 + m0*m0*n2*s2*s2 - 2*j2*m0*q0*s2*s2 - 
  2*j0*m2*q0*s2*s2 + i2*q0*q0*s2*s2 - 2*j0*m0*q2*s2*s2 + 
  2*i0*q0*q2*s2*s2 + 2*l2*m2*o2*p0*t0 + 2*l2*m2*o0*p2*t0 + 
  2*l2*m0*o2*p2*t0 + 2*l0*m2*o2*p2*t0 - 4*k2*m2*p0*p2*t0 - 
  2*k2*m0*p2*p2*t0 - 2*k0*m2*p2*p2*t0 - 2*l2*l2*o2*q0*t0 + 
  2*k2*l2*p2*q0*t0 - 2*l2*l2*o0*q2*t0 - 4*l0*l2*o2*q2*t0 + 
  2*k2*l2*p0*q2*t0 + 2*k2*l0*p2*q2*t0 + 2*k0*l2*p2*q2*t0 - 
  2*l2*m2*n2*s0*t0 + 2*j2*m2*p2*s0*t0 + 2*j2*l2*q2*s0*t0 - 
  2*i2*p2*q2*s0*t0 - 2*l2*m2*n0*s2*t0 - 2*l2*m0*n2*s2*t0 - 
  2*l0*m2*n2*s2*t0 + 2*j2*m2*p0*s2*t0 + 2*j2*m0*p2*s2*t0 + 
  2*j0*m2*p2*s2*t0 + 2*j2*l2*q0*s2*t0 - 2*i2*p2*q0*s2*t0 + 
  2*j2*l0*q2*s2*t0 + 2*j0*l2*q2*s2*t0 - 2*i2*p0*q2*s2*t0 - 
  2*i0*p2*q2*s2*t0 + l2*l2*n2*t0*t0 - 2*j2*l2*p2*t0*t0 + 
  i2*p2*p2*t0*t0 + 2*l2*m2*o0*p0*t2 + 2*l2*m0*o2*p0*t2 + 
  2*l0*m2*o2*p0*t2 - 2*k2*m2*p0*p0*t2 + 2*l2*m0*o0*p2*t2 + 
  2*l0*m2*o0*p2*t2 + 2*l0*m0*o2*p2*t2 - 4*k2*m0*p0*p2*t2 - 
  4*k0*m2*p0*p2*t2 - 2*k0*m0*p2*p2*t2 - 2*l2*l2*o0*q0*t2 - 
  4*l0*l2*o2*q0*t2 + 2*k2*l2*p0*q0*t2 + 2*k2*l0*p2*q0*t2 + 
  2*k0*l2*p2*q0*t2 - 4*l0*l2*o0*q2*t2 - 2*l0*l0*o2*q2*t2 + 
  2*k2*l0*p0*q2*t2 + 2*k0*l2*p0*q2*t2 + 2*k0*l0*p2*q2*t2 - 
  2*l2*m2*n0*s0*t2 - 2*l2*m0*n2*s0*t2 - 2*l0*m2*n2*s0*t2 + 
  2*j2*m2*p0*s0*t2 + 2*j2*m0*p2*s0*t2 + 2*j0*m2*p2*s0*t2 + 
  2*j2*l2*q0*s0*t2 - 2*i2*p2*q0*s0*t2 + 2*j2*l0*q2*s0*t2 + 
  2*j0*l2*q2*s0*t2 - 2*i2*p0*q2*s0*t2 - 2*i0*p2*q2*s0*t2 - 
  2*l2*m0*n0*s2*t2 - 2*l0*m2*n0*s2*t2 - 2*l0*m0*n2*s2*t2 + 
  2*j2*m0*p0*s2*t2 + 2*j0*m2*p0*s2*t2 + 2*j0*m0*p2*s2*t2 + 
  2*j2*l0*q0*s2*t2 + 2*j0*l2*q0*s2*t2 - 2*i2*p0*q0*s2*t2 - 
  2*i0*p2*q0*s2*t2 + 2*j0*l0*q2*s2*t2 - 2*i0*p0*q2*s2*t2 + 
  2*l2*l2*n0*t0*t2 + 4*l0*l2*n2*t0*t2 - 4*j2*l2*p0*t0*t2 - 
  4*j2*l0*p2*t0*t2 - 4*j0*l2*p2*t0*t2 + 4*i2*p0*p2*t0*t2 + 
  2*i0*p2*p2*t0*t2 + 2*l0*l2*n0*t2*t2 + l0*l0*n2*t2*t2 - 
  2*j2*l0*p0*t2*t2 - 2*j0*l2*p0*t2*t2 + i2*p0*p0*t2*t2 - 
  2*j0*l0*p2*t2*t2 + 2*i0*p0*p2*t2*t2 + 2*m2*m2*o0*o2*u0 + 
  2*m0*m2*o2*o2*u0 - 2*k2*m2*o2*q0*u0 - 2*k2*m2*o0*q2*u0 - 
  2*k2*m0*o2*q2*u0 - 2*k0*m2*o2*q2*u0 + 2*k2*k2*q0*q2*u0 + 
  2*k0*k2*q2*q2*u0 - m2*m2*n2*r0*u0 + 2*j2*m2*q2*r0*u0 - 
  i2*q2*q2*r0*u0 - m2*m2*n0*r2*u0 - 2*m0*m2*n2*r2*u0 + 
  2*j2*m2*q0*r2*u0 + 2*j2*m0*q2*r2*u0 + 2*j0*m2*q2*r2*u0 - 
  2*i2*q0*q2*r2*u0 - i0*q2*q2*r2*u0 + 2*k2*m2*n2*t0*u0 - 
  2*j2*m2*o2*t0*u0 - 2*j2*k2*q2*t0*u0 + 2*i2*o2*q2*t0*u0 + 
  2*k2*m2*n0*t2*u0 + 2*k2*m0*n2*t2*u0 + 2*k0*m2*n2*t2*u0 - 
  2*j2*m2*o0*t2*u0 - 2*j2*m0*o2*t2*u0 - 2*j0*m2*o2*t2*u0 - 
  2*j2*k2*q0*t2*u0 + 2*i2*o2*q0*t2*u0 - 2*j2*k0*q2*t2*u0 - 
  2*j0*k2*q2*t2*u0 + 2*i2*o0*q2*t2*u0 + 2*i0*o2*q2*t2*u0 + 
  2*j2*j2*t0*t2*u0 - 2*i2*n2*t0*t2*u0 + 2*j0*j2*t2*t2*u0 - 
  i2*n0*t2*t2*u0 - i0*n2*t2*t2*u0 + m2*m2*o0*o0*u2 + 
  4*m0*m2*o0*o2*u2 + m0*m0*o2*o2*u2 - 2*k2*m2*o0*q0*u2 - 
  2*k2*m0*o2*q0*u2 - 2*k0*m2*o2*q0*u2 + k2*k2*q0*q0*u2 - 
  2*k2*m0*o0*q2*u2 - 2*k0*m2*o0*q2*u2 - 2*k0*m0*o2*q2*u2 + 
  4*k0*k2*q0*q2*u2 + k0*k0*q2*q2*u2 - m2*m2*n0*r0*u2 - 
  2*m0*m2*n2*r0*u2 + 2*j2*m2*q0*r0*u2 + 2*j2*m0*q2*r0*u2 + 
  2*j0*m2*q2*r0*u2 - 2*i2*q0*q2*r0*u2 - i0*q2*q2*r0*u2 - 
  2*m0*m2*n0*r2*u2 - m0*m0*n2*r2*u2 + 2*j2*m0*q0*r2*u2 + 
  2*j0*m2*q0*r2*u2 - i2*q0*q0*r2*u2 + 2*j0*m0*q2*r2*u2 - 
  2*i0*q0*q2*r2*u2 + 2*k2*m2*n0*t0*u2 + 2*k2*m0*n2*t0*u2 + 
  2*k0*m2*n2*t0*u2 - 2*j2*m2*o0*t0*u2 - 2*j2*m0*o2*t0*u2 - 
  2*j0*m2*o2*t0*u2 - 2*j2*k2*q0*t0*u2 + 2*i2*o2*q0*t0*u2 - 
  2*j2*k0*q2*t0*u2 - 2*j0*k2*q2*t0*u2 + 2*i2*o0*q2*t0*u2 + 
  2*i0*o2*q2*t0*u2 + j2*j2*t0*t0*u2 - i2*n2*t0*t0*u2 + 
  2*k2*m0*n0*t2*u2 + 2*k0*m2*n0*t2*u2 + 2*k0*m0*n2*t2*u2 - 
  2*j2*m0*o0*t2*u2 - 2*j0*m2*o0*t2*u2 - 2*j0*m0*o2*t2*u2 - 
  2*j2*k0*q0*t2*u2 - 2*j0*k2*q0*t2*u2 + 2*i2*o0*q0*t2*u2 + 
  2*i0*o2*q0*t2*u2 - 2*j0*k0*q2*t2*u2 + 2*i0*o0*q2*t2*u2 + 
  4*j0*j2*t0*t2*u2 - 2*i2*n0*t0*t2*u2 - 2*i0*n2*t0*t2*u2 + 
  j0*j0*t2*t2*u2 - i0*n0*t2*t2*u2 - 4*l2*m2*o0*o2*v0 - 
  2*l2*m0*o2*o2*v0 - 2*l0*m2*o2*o2*v0 + 2*k2*m2*o2*p0*v0 + 
  2*k2*m2*o0*p2*v0 + 2*k2*m0*o2*p2*v0 + 2*k0*m2*o2*p2*v0 + 
  2*k2*l2*o2*q0*v0 - 2*k2*k2*p2*q0*v0 + 2*k2*l2*o0*q2*v0 + 
  2*k2*l0*o2*q2*v0 + 2*k0*l2*o2*q2*v0 - 2*k2*k2*p0*q2*v0 - 
  4*k0*k2*p2*q2*v0 + 2*l2*m2*n2*r0*v0 - 2*j2*m2*p2*r0*v0 - 
  2*j2*l2*q2*r0*v0 + 2*i2*p2*q2*r0*v0 + 2*l2*m2*n0*r2*v0 + 
  2*l2*m0*n2*r2*v0 + 2*l0*m2*n2*r2*v0 - 2*j2*m2*p0*r2*v0 - 
  2*j2*m0*p2*r2*v0 - 2*j0*m2*p2*r2*v0 - 2*j2*l2*q0*r2*v0 + 
  2*i2*p2*q0*r2*v0 - 2*j2*l0*q2*r2*v0 - 2*j0*l2*q2*r2*v0 + 
  2*i2*p0*q2*r2*v0 + 2*i0*p2*q2*r2*v0 - 2*k2*m2*n2*s0*v0 + 
  2*j2*m2*o2*s0*v0 + 2*j2*k2*q2*s0*v0 - 2*i2*o2*q2*s0*v0 - 
  2*k2*m2*n0*s2*v0 - 2*k2*m0*n2*s2*v0 - 2*k0*m2*n2*s2*v0 + 
  2*j2*m2*o0*s2*v0 + 2*j2*m0*o2*s2*v0 + 2*j0*m2*o2*s2*v0 + 
  2*j2*k2*q0*s2*v0 - 2*i2*o2*q0*s2*v0 + 2*j2*k0*q2*s2*v0 + 
  2*j0*k2*q2*s2*v0 - 2*i2*o0*q2*s2*v0 - 2*i0*o2*q2*s2*v0 - 
  2*k2*l2*n2*t0*v0 + 2*j2*l2*o2*t0*v0 + 2*j2*k2*p2*t0*v0 - 
  2*i2*o2*p2*t0*v0 - 2*j2*j2*s2*t0*v0 + 2*i2*n2*s2*t0*v0 - 
  2*k2*l2*n0*t2*v0 - 2*k2*l0*n2*t2*v0 - 2*k0*l2*n2*t2*v0 + 
  2*j2*l2*o0*t2*v0 + 2*j2*l0*o2*t2*v0 + 2*j0*l2*o2*t2*v0 + 
  2*j2*k2*p0*t2*v0 - 2*i2*o2*p0*t2*v0 + 2*j2*k0*p2*t2*v0 + 
  2*j0*k2*p2*t2*v0 - 2*i2*o0*p2*t2*v0 - 2*i0*o2*p2*t2*v0 - 
  2*j2*j2*s0*t2*v0 + 2*i2*n2*s0*t2*v0 - 4*j0*j2*s2*t2*v0 + 
  2*i2*n0*s2*t2*v0 + 2*i0*n2*s2*t2*v0 + k2*k2*n2*v0*v0 - 
  2*j2*k2*o2*v0*v0 + i2*o2*o2*v0*v0 + j2*j2*r2*v0*v0 - 
  i2*n2*r2*v0*v0 - 2*l2*m2*o0*o0*v2 - 4*l2*m0*o0*o2*v2 - 
  4*l0*m2*o0*o2*v2 - 2*l0*m0*o2*o2*v2 + 2*k2*m2*o0*p0*v2 + 
  2*k2*m0*o2*p0*v2 + 2*k0*m2*o2*p0*v2 + 2*k2*m0*o0*p2*v2 + 
  2*k0*m2*o0*p2*v2 + 2*k0*m0*o2*p2*v2 + 2*k2*l2*o0*q0*v2 + 
  2*k2*l0*o2*q0*v2 + 2*k0*l2*o2*q0*v2 - 2*k2*k2*p0*q0*v2 - 
  4*k0*k2*p2*q0*v2 + 2*k2*l0*o0*q2*v2 + 2*k0*l2*o0*q2*v2 + 
  2*k0*l0*o2*q2*v2 - 4*k0*k2*p0*q2*v2 - 2*k0*k0*p2*q2*v2 + 
  2*l2*m2*n0*r0*v2 + 2*l2*m0*n2*r0*v2 + 2*l0*m2*n2*r0*v2 - 
  2*j2*m2*p0*r0*v2 - 2*j2*m0*p2*r0*v2 - 2*j0*m2*p2*r0*v2 - 
  2*j2*l2*q0*r0*v2 + 2*i2*p2*q0*r0*v2 - 2*j2*l0*q2*r0*v2 - 
  2*j0*l2*q2*r0*v2 + 2*i2*p0*q2*r0*v2 + 2*i0*p2*q2*r0*v2 + 
  2*l2*m0*n0*r2*v2 + 2*l0*m2*n0*r2*v2 + 2*l0*m0*n2*r2*v2 - 
  2*j2*m0*p0*r2*v2 - 2*j0*m2*p0*r2*v2 - 2*j0*m0*p2*r2*v2 - 
  2*j2*l0*q0*r2*v2 - 2*j0*l2*q0*r2*v2 + 2*i2*p0*q0*r2*v2 + 
  2*i0*p2*q0*r2*v2 - 2*j0*l0*q2*r2*v2 + 2*i0*p0*q2*r2*v2 - 
  2*k2*m2*n0*s0*v2 - 2*k2*m0*n2*s0*v2 - 2*k0*m2*n2*s0*v2 + 
  2*j2*m2*o0*s0*v2 + 2*j2*m0*o2*s0*v2 + 2*j0*m2*o2*s0*v2 + 
  2*j2*k2*q0*s0*v2 - 2*i2*o2*q0*s0*v2 + 2*j2*k0*q2*s0*v2 + 
  2*j0*k2*q2*s0*v2 - 2*i2*o0*q2*s0*v2 - 2*i0*o2*q2*s0*v2 - 
  2*k2*m0*n0*s2*v2 - 2*k0*m2*n0*s2*v2 - 2*k0*m0*n2*s2*v2 + 
  2*j2*m0*o0*s2*v2 + 2*j0*m2*o0*s2*v2 + 2*j0*m0*o2*s2*v2 + 
  2*j2*k0*q0*s2*v2 + 2*j0*k2*q0*s2*v2 - 2*i2*o0*q0*s2*v2 - 
  2*i0*o2*q0*s2*v2 + 2*j0*k0*q2*s2*v2 - 2*i0*o0*q2*s2*v2 - 
  2*k2*l2*n0*t0*v2 - 2*k2*l0*n2*t0*v2 - 2*k0*l2*n2*t0*v2 + 
  2*j2*l2*o0*t0*v2 + 2*j2*l0*o2*t0*v2 + 2*j0*l2*o2*t0*v2 + 
  2*j2*k2*p0*t0*v2 - 2*i2*o2*p0*t0*v2 + 2*j2*k0*p2*t0*v2 + 
  2*j0*k2*p2*t0*v2 - 2*i2*o0*p2*t0*v2 - 2*i0*o2*p2*t0*v2 - 
  2*j2*j2*s0*t0*v2 + 2*i2*n2*s0*t0*v2 - 4*j0*j2*s2*t0*v2 + 
  2*i2*n0*s2*t0*v2 + 2*i0*n2*s2*t0*v2 - 2*k2*l0*n0*t2*v2 - 
  2*k0*l2*n0*t2*v2 - 2*k0*l0*n2*t2*v2 + 2*j2*l0*o0*t2*v2 + 
  2*j0*l2*o0*t2*v2 + 2*j0*l0*o2*t2*v2 + 2*j2*k0*p0*t2*v2 + 
  2*j0*k2*p0*t2*v2 - 2*i2*o0*p0*t2*v2 - 2*i0*o2*p0*t2*v2 + 
  2*j0*k0*p2*t2*v2 - 2*i0*o0*p2*t2*v2 - 4*j0*j2*s0*t2*v2 + 
  2*i2*n0*s0*t2*v2 + 2*i0*n2*s0*t2*v2 - 2*j0*j0*s2*t2*v2 + 
  2*i0*n0*s2*t2*v2 + 2*k2*k2*n0*v0*v2 + 4*k0*k2*n2*v0*v2 - 
  4*j2*k2*o0*v0*v2 - 4*j2*k0*o2*v0*v2 - 4*j0*k2*o2*v0*v2 + 
  4*i2*o0*o2*v0*v2 + 2*i0*o2*o2*v0*v2 + 2*j2*j2*r0*v0*v2 - 
  2*i2*n2*r0*v0*v2 + 4*j0*j2*r2*v0*v2 - 2*i2*n0*r2*v0*v2 - 
  2*i0*n2*r2*v0*v2 + 2*k0*k2*n0*v2*v2 + k0*k0*n2*v2*v2 - 
  2*j2*k0*o0*v2*v2 - 2*j0*k2*o0*v2*v2 + i2*o0*o0*v2*v2 - 
  2*j0*k0*o2*v2*v2 + 2*i0*o0*o2*v2*v2 + 2*j0*j2*r0*v2*v2 - 
  i2*n0*r0*v2*v2 - i0*n2*r0*v2*v2 + j0*j0*r2*v2*v2 - 
  i0*n0*r2*v2*v2 + 2*l2*l2*o0*o2*w0 + 2*l0*l2*o2*o2*w0 - 
  2*k2*l2*o2*p0*w0 - 2*k2*l2*o0*p2*w0 - 2*k2*l0*o2*p2*w0 - 
  2*k0*l2*o2*p2*w0 + 2*k2*k2*p0*p2*w0 + 2*k0*k2*p2*p2*w0 - 
  l2*l2*n2*r0*w0 + 2*j2*l2*p2*r0*w0 - i2*p2*p2*r0*w0 - 
  l2*l2*n0*r2*w0 - 2*l0*l2*n2*r2*w0 + 2*j2*l2*p0*r2*w0 + 
  2*j2*l0*p2*r2*w0 + 2*j0*l2*p2*r2*w0 - 2*i2*p0*p2*r2*w0 - 
  i0*p2*p2*r2*w0 + 2*k2*l2*n2*s0*w0 - 2*j2*l2*o2*s0*w0 - 
  2*j2*k2*p2*s0*w0 + 2*i2*o2*p2*s0*w0 + 2*k2*l2*n0*s2*w0 + 
  2*k2*l0*n2*s2*w0 + 2*k0*l2*n2*s2*w0 - 2*j2*l2*o0*s2*w0 - 
  2*j2*l0*o2*s2*w0 - 2*j0*l2*o2*s2*w0 - 2*j2*k2*p0*s2*w0 + 
  2*i2*o2*p0*s2*w0 - 2*j2*k0*p2*s2*w0 - 2*j0*k2*p2*s2*w0 + 
  2*i2*o0*p2*s2*w0 + 2*i0*o2*p2*s2*w0 + 2*j2*j2*s0*s2*w0 - 
  2*i2*n2*s0*s2*w0 + 2*j0*j2*s2*s2*w0 - i2*n0*s2*s2*w0 - 
  i0*n2*s2*s2*w0 - k2*k2*n2*u0*w0 + 2*j2*k2*o2*u0*w0 - 
  i2*o2*o2*u0*w0 - j2*j2*r2*u0*w0 + i2*n2*r2*u0*w0 - 
  k2*k2*n0*u2*w0 - 2*k0*k2*n2*u2*w0 + 2*j2*k2*o0*u2*w0 + 
  2*j2*k0*o2*u2*w0 + 2*j0*k2*o2*u2*w0 - 2*i2*o0*o2*u2*w0 - 
  i0*o2*o2*u2*w0 - j2*j2*r0*u2*w0 + i2*n2*r0*u2*w0 - 
  2*j0*j2*r2*u2*w0 + i2*n0*r2*u2*w0 + i0*n2*r2*u2*w0 + 
  l2*l2*o0*o0*w2 + 4*l0*l2*o0*o2*w2 + l0*l0*o2*o2*w2 - 
  2*k2*l2*o0*p0*w2 - 2*k2*l0*o2*p0*w2 - 2*k0*l2*o2*p0*w2 + 
  k2*k2*p0*p0*w2 - 2*k2*l0*o0*p2*w2 - 2*k0*l2*o0*p2*w2 - 
  2*k0*l0*o2*p2*w2 + 4*k0*k2*p0*p2*w2 + k0*k0*p2*p2*w2 - 
  l2*l2*n0*r0*w2 - 2*l0*l2*n2*r0*w2 + 2*j2*l2*p0*r0*w2 + 
  2*j2*l0*p2*r0*w2 + 2*j0*l2*p2*r0*w2 - 2*i2*p0*p2*r0*w2 - 
  i0*p2*p2*r0*w2 - 2*l0*l2*n0*r2*w2 - l0*l0*n2*r2*w2 + 
  2*j2*l0*p0*r2*w2 + 2*j0*l2*p0*r2*w2 - i2*p0*p0*r2*w2 + 
  2*j0*l0*p2*r2*w2 - 2*i0*p0*p2*r2*w2 + 2*k2*l2*n0*s0*w2 + 
  2*k2*l0*n2*s0*w2 + 2*k0*l2*n2*s0*w2 - 2*j2*l2*o0*s0*w2 - 
  2*j2*l0*o2*s0*w2 - 2*j0*l2*o2*s0*w2 - 2*j2*k2*p0*s0*w2 + 
  2*i2*o2*p0*s0*w2 - 2*j2*k0*p2*s0*w2 - 2*j0*k2*p2*s0*w2 + 
  2*i2*o0*p2*s0*w2 + 2*i0*o2*p2*s0*w2 + j2*j2*s0*s0*w2 - 
  i2*n2*s0*s0*w2 + 2*k2*l0*n0*s2*w2 + 2*k0*l2*n0*s2*w2 + 
  2*k0*l0*n2*s2*w2 - 2*j2*l0*o0*s2*w2 - 2*j0*l2*o0*s2*w2 - 
  2*j0*l0*o2*s2*w2 - 2*j2*k0*p0*s2*w2 - 2*j0*k2*p0*s2*w2 + 
  2*i2*o0*p0*s2*w2 + 2*i0*o2*p0*s2*w2 - 2*j0*k0*p2*s2*w2 + 
  2*i0*o0*p2*s2*w2 + 4*j0*j2*s0*s2*w2 - 2*i2*n0*s0*s2*w2 - 
  2*i0*n2*s0*s2*w2 + j0*j0*s2*s2*w2 - i0*n0*s2*s2*w2 - 
  k2*k2*n0*u0*w2 - 2*k0*k2*n2*u0*w2 + 2*j2*k2*o0*u0*w2 + 
  2*j2*k0*o2*u0*w2 + 2*j0*k2*o2*u0*w2 - 2*i2*o0*o2*u0*w2 - 
  i0*o2*o2*u0*w2 - j2*j2*r0*u0*w2 + i2*n2*r0*u0*w2 - 
  2*j0*j2*r2*u0*w2 + i2*n0*r2*u0*w2 + i0*n2*r2*u0*w2 - 
  2*k0*k2*n0*u2*w2 - k0*k0*n2*u2*w2 + 2*j2*k0*o0*u2*w2 + 
  2*j0*k2*o0*u2*w2 - i2*o0*o0*u2*w2 + 2*j0*k0*o2*u2*w2 - 
  2*i0*o0*o2*u2*w2 - 2*j0*j2*r0*u2*w2 + i2*n0*r0*u2*w2 + 
  i0*n2*r0*u2*w2 - j0*j0*r2*u2*w2 + i0*n0*r2*u2*w2;
    
    k[4]  = 2*m2*m2*p1*p2*r0 + 
  2*m1*m2*p2*p2*r0 - 2*l2*m2*p2*q1*r0 - 2*l2*m2*p1*q2*r0 - 
  2*l2*m1*p2*q2*r0 - 2*l1*m2*p2*q2*r0 + 2*l2*l2*q1*q2*r0 + 
  2*l1*l2*q2*q2*r0 + 2*m2*m2*p0*p2*r1 + 2*m0*m2*p2*p2*r1 - 
  2*l2*m2*p2*q0*r1 - 2*l2*m2*p0*q2*r1 - 2*l2*m0*p2*q2*r1 - 
  2*l0*m2*p2*q2*r1 + 2*l2*l2*q0*q2*r1 + 2*l0*l2*q2*q2*r1 + 
  2*m2*m2*p0*p1*r2 + 4*m1*m2*p0*p2*r2 + 4*m0*m2*p1*p2*r2 + 
  2*m0*m1*p2*p2*r2 - 2*l2*m2*p1*q0*r2 - 2*l2*m1*p2*q0*r2 - 
  2*l1*m2*p2*q0*r2 - 2*l2*m2*p0*q1*r2 - 2*l2*m0*p2*q1*r2 - 
  2*l0*m2*p2*q1*r2 + 2*l2*l2*q0*q1*r2 - 2*l2*m1*p0*q2*r2 - 
  2*l1*m2*p0*q2*r2 - 2*l2*m0*p1*q2*r2 - 2*l0*m2*p1*q2*r2 - 
  2*l1*m0*p2*q2*r2 - 2*l0*m1*p2*q2*r2 + 4*l1*l2*q0*q2*r2 + 
  4*l0*l2*q1*q2*r2 + 2*l0*l1*q2*q2*r2 - 2*m2*m2*o2*p1*s0 - 
  2*m2*m2*o1*p2*s0 - 4*m1*m2*o2*p2*s0 + 2*l2*m2*o2*q1*s0 + 
  2*k2*m2*p2*q1*s0 + 2*l2*m2*o1*q2*s0 + 2*l2*m1*o2*q2*s0 + 
  2*l1*m2*o2*q2*s0 + 2*k2*m2*p1*q2*s0 + 2*k2*m1*p2*q2*s0 + 
  2*k1*m2*p2*q2*s0 - 4*k2*l2*q1*q2*s0 - 2*k2*l1*q2*q2*s0 - 
  2*k1*l2*q2*q2*s0 - 2*m2*m2*o2*p0*s1 - 2*m2*m2*o0*p2*s1 - 
  4*m0*m2*o2*p2*s1 + 2*l2*m2*o2*q0*s1 + 2*k2*m2*p2*q0*s1 + 
  2*l2*m2*o0*q2*s1 + 2*l2*m0*o2*q2*s1 + 2*l0*m2*o2*q2*s1 + 
  2*k2*m2*p0*q2*s1 + 2*k2*m0*p2*q2*s1 + 2*k0*m2*p2*q2*s1 - 
  4*k2*l2*q0*q2*s1 - 2*k2*l0*q2*q2*s1 - 2*k0*l2*q2*q2*s1 + 
  2*m2*m2*n2*s0*s1 - 4*j2*m2*q2*s0*s1 + 2*i2*q2*q2*s0*s1 - 
  2*m2*m2*o1*p0*s2 - 4*m1*m2*o2*p0*s2 - 2*m2*m2*o0*p1*s2 - 
  4*m0*m2*o2*p1*s2 - 4*m1*m2*o0*p2*s2 - 4*m0*m2*o1*p2*s2 - 
  4*m0*m1*o2*p2*s2 + 2*l2*m2*o1*q0*s2 + 2*l2*m1*o2*q0*s2 + 
  2*l1*m2*o2*q0*s2 + 2*k2*m2*p1*q0*s2 + 2*k2*m1*p2*q0*s2 + 
  2*k1*m2*p2*q0*s2 + 2*l2*m2*o0*q1*s2 + 2*l2*m0*o2*q1*s2 + 
  2*l0*m2*o2*q1*s2 + 2*k2*m2*p0*q1*s2 + 2*k2*m0*p2*q1*s2 + 
  2*k0*m2*p2*q1*s2 - 4*k2*l2*q0*q1*s2 + 2*l2*m1*o0*q2*s2 + 
  2*l1*m2*o0*q2*s2 + 2*l2*m0*o1*q2*s2 + 2*l0*m2*o1*q2*s2 + 
  2*l1*m0*o2*q2*s2 + 2*l0*m1*o2*q2*s2 + 2*k2*m1*p0*q2*s2 + 
  2*k1*m2*p0*q2*s2 + 2*k2*m0*p1*q2*s2 + 2*k0*m2*p1*q2*s2 + 
  2*k1*m0*p2*q2*s2 + 2*k0*m1*p2*q2*s2 - 4*k2*l1*q0*q2*s2 - 
  4*k1*l2*q0*q2*s2 - 4*k2*l0*q1*q2*s2 - 4*k0*l2*q1*q2*s2 - 
  2*k1*l0*q2*q2*s2 - 2*k0*l1*q2*q2*s2 + 2*m2*m2*n1*s0*s2 + 
  4*m1*m2*n2*s0*s2 - 4*j2*m2*q1*s0*s2 - 4*j2*m1*q2*s0*s2 - 
  4*j1*m2*q2*s0*s2 + 4*i2*q1*q2*s0*s2 + 2*i1*q2*q2*s0*s2 + 
  2*m2*m2*n0*s1*s2 + 4*m0*m2*n2*s1*s2 - 4*j2*m2*q0*s1*s2 - 
  4*j2*m0*q2*s1*s2 - 4*j0*m2*q2*s1*s2 + 4*i2*q0*q2*s1*s2 + 
  2*i0*q2*q2*s1*s2 + 2*m1*m2*n0*s2*s2 + 2*m0*m2*n1*s2*s2 + 
  2*m0*m1*n2*s2*s2 - 2*j2*m1*q0*s2*s2 - 2*j1*m2*q0*s2*s2 - 
  2*j2*m0*q1*s2*s2 - 2*j0*m2*q1*s2*s2 + 2*i2*q0*q1*s2*s2 - 
  2*j1*m0*q2*s2*s2 - 2*j0*m1*q2*s2*s2 + 2*i1*q0*q2*s2*s2 + 
  2*i0*q1*q2*s2*s2 + 2*l2*m2*o2*p1*t0 + 2*l2*m2*o1*p2*t0 + 
  2*l2*m1*o2*p2*t0 + 2*l1*m2*o2*p2*t0 - 4*k2*m2*p1*p2*t0 - 
  2*k2*m1*p2*p2*t0 - 2*k1*m2*p2*p2*t0 - 2*l2*l2*o2*q1*t0 + 
  2*k2*l2*p2*q1*t0 - 2*l2*l2*o1*q2*t0 - 4*l1*l2*o2*q2*t0 + 
  2*k2*l2*p1*q2*t0 + 2*k2*l1*p2*q2*t0 + 2*k1*l2*p2*q2*t0 - 
  2*l2*m2*n2*s1*t0 + 2*j2*m2*p2*s1*t0 + 2*j2*l2*q2*s1*t0 - 
  2*i2*p2*q2*s1*t0 - 2*l2*m2*n1*s2*t0 - 2*l2*m1*n2*s2*t0 - 
  2*l1*m2*n2*s2*t0 + 2*j2*m2*p1*s2*t0 + 2*j2*m1*p2*s2*t0 + 
  2*j1*m2*p2*s2*t0 + 2*j2*l2*q1*s2*t0 - 2*i2*p2*q1*s2*t0 + 
  2*j2*l1*q2*s2*t0 + 2*j1*l2*q2*s2*t0 - 2*i2*p1*q2*s2*t0 - 
  2*i1*p2*q2*s2*t0 + 2*l2*m2*o2*p0*t1 + 2*l2*m2*o0*p2*t1 + 
  2*l2*m0*o2*p2*t1 + 2*l0*m2*o2*p2*t1 - 4*k2*m2*p0*p2*t1 - 
  2*k2*m0*p2*p2*t1 - 2*k0*m2*p2*p2*t1 - 2*l2*l2*o2*q0*t1 + 
  2*k2*l2*p2*q0*t1 - 2*l2*l2*o0*q2*t1 - 4*l0*l2*o2*q2*t1 + 
  2*k2*l2*p0*q2*t1 + 2*k2*l0*p2*q2*t1 + 2*k0*l2*p2*q2*t1 - 
  2*l2*m2*n2*s0*t1 + 2*j2*m2*p2*s0*t1 + 2*j2*l2*q2*s0*t1 - 
  2*i2*p2*q2*s0*t1 - 2*l2*m2*n0*s2*t1 - 2*l2*m0*n2*s2*t1 - 
  2*l0*m2*n2*s2*t1 + 2*j2*m2*p0*s2*t1 + 2*j2*m0*p2*s2*t1 + 
  2*j0*m2*p2*s2*t1 + 2*j2*l2*q0*s2*t1 - 2*i2*p2*q0*s2*t1 + 
  2*j2*l0*q2*s2*t1 + 2*j0*l2*q2*s2*t1 - 2*i2*p0*q2*s2*t1 - 
  2*i0*p2*q2*s2*t1 + 2*l2*l2*n2*t0*t1 - 4*j2*l2*p2*t0*t1 + 
  2*i2*p2*p2*t0*t1 + 2*l2*m2*o1*p0*t2 + 2*l2*m1*o2*p0*t2 + 
  2*l1*m2*o2*p0*t2 + 2*l2*m2*o0*p1*t2 + 2*l2*m0*o2*p1*t2 + 
  2*l0*m2*o2*p1*t2 - 4*k2*m2*p0*p1*t2 + 2*l2*m1*o0*p2*t2 + 
  2*l1*m2*o0*p2*t2 + 2*l2*m0*o1*p2*t2 + 2*l0*m2*o1*p2*t2 + 
  2*l1*m0*o2*p2*t2 + 2*l0*m1*o2*p2*t2 - 4*k2*m1*p0*p2*t2 - 
  4*k1*m2*p0*p2*t2 - 4*k2*m0*p1*p2*t2 - 4*k0*m2*p1*p2*t2 - 
  2*k1*m0*p2*p2*t2 - 2*k0*m1*p2*p2*t2 - 2*l2*l2*o1*q0*t2 - 
  4*l1*l2*o2*q0*t2 + 2*k2*l2*p1*q0*t2 + 2*k2*l1*p2*q0*t2 + 
  2*k1*l2*p2*q0*t2 - 2*l2*l2*o0*q1*t2 - 4*l0*l2*o2*q1*t2 + 
  2*k2*l2*p0*q1*t2 + 2*k2*l0*p2*q1*t2 + 2*k0*l2*p2*q1*t2 - 
  4*l1*l2*o0*q2*t2 - 4*l0*l2*o1*q2*t2 - 4*l0*l1*o2*q2*t2 + 
  2*k2*l1*p0*q2*t2 + 2*k1*l2*p0*q2*t2 + 2*k2*l0*p1*q2*t2 + 
  2*k0*l2*p1*q2*t2 + 2*k1*l0*p2*q2*t2 + 2*k0*l1*p2*q2*t2 - 
  2*l2*m2*n1*s0*t2 - 2*l2*m1*n2*s0*t2 - 2*l1*m2*n2*s0*t2 + 
  2*j2*m2*p1*s0*t2 + 2*j2*m1*p2*s0*t2 + 2*j1*m2*p2*s0*t2 + 
  2*j2*l2*q1*s0*t2 - 2*i2*p2*q1*s0*t2 + 2*j2*l1*q2*s0*t2 + 
  2*j1*l2*q2*s0*t2 - 2*i2*p1*q2*s0*t2 - 2*i1*p2*q2*s0*t2 - 
  2*l2*m2*n0*s1*t2 - 2*l2*m0*n2*s1*t2 - 2*l0*m2*n2*s1*t2 + 
  2*j2*m2*p0*s1*t2 + 2*j2*m0*p2*s1*t2 + 2*j0*m2*p2*s1*t2 + 
  2*j2*l2*q0*s1*t2 - 2*i2*p2*q0*s1*t2 + 2*j2*l0*q2*s1*t2 + 
  2*j0*l2*q2*s1*t2 - 2*i2*p0*q2*s1*t2 - 2*i0*p2*q2*s1*t2 - 
  2*l2*m1*n0*s2*t2 - 2*l1*m2*n0*s2*t2 - 2*l2*m0*n1*s2*t2 - 
  2*l0*m2*n1*s2*t2 - 2*l1*m0*n2*s2*t2 - 2*l0*m1*n2*s2*t2 + 
  2*j2*m1*p0*s2*t2 + 2*j1*m2*p0*s2*t2 + 2*j2*m0*p1*s2*t2 + 
  2*j0*m2*p1*s2*t2 + 2*j1*m0*p2*s2*t2 + 2*j0*m1*p2*s2*t2 + 
  2*j2*l1*q0*s2*t2 + 2*j1*l2*q0*s2*t2 - 2*i2*p1*q0*s2*t2 - 
  2*i1*p2*q0*s2*t2 + 2*j2*l0*q1*s2*t2 + 2*j0*l2*q1*s2*t2 - 
  2*i2*p0*q1*s2*t2 - 2*i0*p2*q1*s2*t2 + 2*j1*l0*q2*s2*t2 + 
  2*j0*l1*q2*s2*t2 - 2*i1*p0*q2*s2*t2 - 2*i0*p1*q2*s2*t2 + 
  2*l2*l2*n1*t0*t2 + 4*l1*l2*n2*t0*t2 - 4*j2*l2*p1*t0*t2 - 
  4*j2*l1*p2*t0*t2 - 4*j1*l2*p2*t0*t2 + 4*i2*p1*p2*t0*t2 + 
  2*i1*p2*p2*t0*t2 + 2*l2*l2*n0*t1*t2 + 4*l0*l2*n2*t1*t2 - 
  4*j2*l2*p0*t1*t2 - 4*j2*l0*p2*t1*t2 - 4*j0*l2*p2*t1*t2 + 
  4*i2*p0*p2*t1*t2 + 2*i0*p2*p2*t1*t2 + 2*l1*l2*n0*t2*t2 + 
  2*l0*l2*n1*t2*t2 + 2*l0*l1*n2*t2*t2 - 2*j2*l1*p0*t2*t2 - 
  2*j1*l2*p0*t2*t2 - 2*j2*l0*p1*t2*t2 - 2*j0*l2*p1*t2*t2 + 
  2*i2*p0*p1*t2*t2 - 2*j1*l0*p2*t2*t2 - 2*j0*l1*p2*t2*t2 + 
  2*i1*p0*p2*t2*t2 + 2*i0*p1*p2*t2*t2 + 2*m2*m2*o1*o2*u0 + 
  2*m1*m2*o2*o2*u0 - 2*k2*m2*o2*q1*u0 - 2*k2*m2*o1*q2*u0 - 
  2*k2*m1*o2*q2*u0 - 2*k1*m2*o2*q2*u0 + 2*k2*k2*q1*q2*u0 + 
  2*k1*k2*q2*q2*u0 - m2*m2*n2*r1*u0 + 2*j2*m2*q2*r1*u0 - 
  i2*q2*q2*r1*u0 - m2*m2*n1*r2*u0 - 2*m1*m2*n2*r2*u0 + 
  2*j2*m2*q1*r2*u0 + 2*j2*m1*q2*r2*u0 + 2*j1*m2*q2*r2*u0 - 
  2*i2*q1*q2*r2*u0 - i1*q2*q2*r2*u0 + 2*k2*m2*n2*t1*u0 - 
  2*j2*m2*o2*t1*u0 - 2*j2*k2*q2*t1*u0 + 2*i2*o2*q2*t1*u0 + 
  2*k2*m2*n1*t2*u0 + 2*k2*m1*n2*t2*u0 + 2*k1*m2*n2*t2*u0 - 
  2*j2*m2*o1*t2*u0 - 2*j2*m1*o2*t2*u0 - 2*j1*m2*o2*t2*u0 - 
  2*j2*k2*q1*t2*u0 + 2*i2*o2*q1*t2*u0 - 2*j2*k1*q2*t2*u0 - 
  2*j1*k2*q2*t2*u0 + 2*i2*o1*q2*t2*u0 + 2*i1*o2*q2*t2*u0 + 
  2*j2*j2*t1*t2*u0 - 2*i2*n2*t1*t2*u0 + 2*j1*j2*t2*t2*u0 - 
  i2*n1*t2*t2*u0 - i1*n2*t2*t2*u0 + 2*m2*m2*o0*o2*u1 + 
  2*m0*m2*o2*o2*u1 - 2*k2*m2*o2*q0*u1 - 2*k2*m2*o0*q2*u1 - 
  2*k2*m0*o2*q2*u1 - 2*k0*m2*o2*q2*u1 + 2*k2*k2*q0*q2*u1 + 
  2*k0*k2*q2*q2*u1 - m2*m2*n2*r0*u1 + 2*j2*m2*q2*r0*u1 - 
  i2*q2*q2*r0*u1 - m2*m2*n0*r2*u1 - 2*m0*m2*n2*r2*u1 + 
  2*j2*m2*q0*r2*u1 + 2*j2*m0*q2*r2*u1 + 2*j0*m2*q2*r2*u1 - 
  2*i2*q0*q2*r2*u1 - i0*q2*q2*r2*u1 + 2*k2*m2*n2*t0*u1 - 
  2*j2*m2*o2*t0*u1 - 2*j2*k2*q2*t0*u1 + 2*i2*o2*q2*t0*u1 + 
  2*k2*m2*n0*t2*u1 + 2*k2*m0*n2*t2*u1 + 2*k0*m2*n2*t2*u1 - 
  2*j2*m2*o0*t2*u1 - 2*j2*m0*o2*t2*u1 - 2*j0*m2*o2*t2*u1 - 
  2*j2*k2*q0*t2*u1 + 2*i2*o2*q0*t2*u1 - 2*j2*k0*q2*t2*u1 - 
  2*j0*k2*q2*t2*u1 + 2*i2*o0*q2*t2*u1 + 2*i0*o2*q2*t2*u1 + 
  2*j2*j2*t0*t2*u1 - 2*i2*n2*t0*t2*u1 + 2*j0*j2*t2*t2*u1 - 
  i2*n0*t2*t2*u1 - i0*n2*t2*t2*u1 + 2*m2*m2*o0*o1*u2 + 
  4*m1*m2*o0*o2*u2 + 4*m0*m2*o1*o2*u2 + 2*m0*m1*o2*o2*u2 - 
  2*k2*m2*o1*q0*u2 - 2*k2*m1*o2*q0*u2 - 2*k1*m2*o2*q0*u2 - 
  2*k2*m2*o0*q1*u2 - 2*k2*m0*o2*q1*u2 - 2*k0*m2*o2*q1*u2 + 
  2*k2*k2*q0*q1*u2 - 2*k2*m1*o0*q2*u2 - 2*k1*m2*o0*q2*u2 - 
  2*k2*m0*o1*q2*u2 - 2*k0*m2*o1*q2*u2 - 2*k1*m0*o2*q2*u2 - 
  2*k0*m1*o2*q2*u2 + 4*k1*k2*q0*q2*u2 + 4*k0*k2*q1*q2*u2 + 
  2*k0*k1*q2*q2*u2 - m2*m2*n1*r0*u2 - 2*m1*m2*n2*r0*u2 + 
  2*j2*m2*q1*r0*u2 + 2*j2*m1*q2*r0*u2 + 2*j1*m2*q2*r0*u2 - 
  2*i2*q1*q2*r0*u2 - i1*q2*q2*r0*u2 - m2*m2*n0*r1*u2 - 
  2*m0*m2*n2*r1*u2 + 2*j2*m2*q0*r1*u2 + 2*j2*m0*q2*r1*u2 + 
  2*j0*m2*q2*r1*u2 - 2*i2*q0*q2*r1*u2 - i0*q2*q2*r1*u2 - 
  2*m1*m2*n0*r2*u2 - 2*m0*m2*n1*r2*u2 - 2*m0*m1*n2*r2*u2 + 
  2*j2*m1*q0*r2*u2 + 2*j1*m2*q0*r2*u2 + 2*j2*m0*q1*r2*u2 + 
  2*j0*m2*q1*r2*u2 - 2*i2*q0*q1*r2*u2 + 2*j1*m0*q2*r2*u2 + 
  2*j0*m1*q2*r2*u2 - 2*i1*q0*q2*r2*u2 - 2*i0*q1*q2*r2*u2 + 
  2*k2*m2*n1*t0*u2 + 2*k2*m1*n2*t0*u2 + 2*k1*m2*n2*t0*u2 - 
  2*j2*m2*o1*t0*u2 - 2*j2*m1*o2*t0*u2 - 2*j1*m2*o2*t0*u2 - 
  2*j2*k2*q1*t0*u2 + 2*i2*o2*q1*t0*u2 - 2*j2*k1*q2*t0*u2 - 
  2*j1*k2*q2*t0*u2 + 2*i2*o1*q2*t0*u2 + 2*i1*o2*q2*t0*u2 + 
  2*k2*m2*n0*t1*u2 + 2*k2*m0*n2*t1*u2 + 2*k0*m2*n2*t1*u2 - 
  2*j2*m2*o0*t1*u2 - 2*j2*m0*o2*t1*u2 - 2*j0*m2*o2*t1*u2 - 
  2*j2*k2*q0*t1*u2 + 2*i2*o2*q0*t1*u2 - 2*j2*k0*q2*t1*u2 - 
  2*j0*k2*q2*t1*u2 + 2*i2*o0*q2*t1*u2 + 2*i0*o2*q2*t1*u2 + 
  2*j2*j2*t0*t1*u2 - 2*i2*n2*t0*t1*u2 + 2*k2*m1*n0*t2*u2 + 
  2*k1*m2*n0*t2*u2 + 2*k2*m0*n1*t2*u2 + 2*k0*m2*n1*t2*u2 + 
  2*k1*m0*n2*t2*u2 + 2*k0*m1*n2*t2*u2 - 2*j2*m1*o0*t2*u2 - 
  2*j1*m2*o0*t2*u2 - 2*j2*m0*o1*t2*u2 - 2*j0*m2*o1*t2*u2 - 
  2*j1*m0*o2*t2*u2 - 2*j0*m1*o2*t2*u2 - 2*j2*k1*q0*t2*u2 - 
  2*j1*k2*q0*t2*u2 + 2*i2*o1*q0*t2*u2 + 2*i1*o2*q0*t2*u2 - 
  2*j2*k0*q1*t2*u2 - 2*j0*k2*q1*t2*u2 + 2*i2*o0*q1*t2*u2 + 
  2*i0*o2*q1*t2*u2 - 2*j1*k0*q2*t2*u2 - 2*j0*k1*q2*t2*u2 + 
  2*i1*o0*q2*t2*u2 + 2*i0*o1*q2*t2*u2 + 4*j1*j2*t0*t2*u2 - 
  2*i2*n1*t0*t2*u2 - 2*i1*n2*t0*t2*u2 + 4*j0*j2*t1*t2*u2 - 
  2*i2*n0*t1*t2*u2 - 2*i0*n2*t1*t2*u2 + 2*j0*j1*t2*t2*u2 - 
  i1*n0*t2*t2*u2 - i0*n1*t2*t2*u2 - 4*l2*m2*o1*o2*v0 - 
  2*l2*m1*o2*o2*v0 - 2*l1*m2*o2*o2*v0 + 2*k2*m2*o2*p1*v0 + 
  2*k2*m2*o1*p2*v0 + 2*k2*m1*o2*p2*v0 + 2*k1*m2*o2*p2*v0 + 
  2*k2*l2*o2*q1*v0 - 2*k2*k2*p2*q1*v0 + 2*k2*l2*o1*q2*v0 + 
  2*k2*l1*o2*q2*v0 + 2*k1*l2*o2*q2*v0 - 2*k2*k2*p1*q2*v0 - 
  4*k1*k2*p2*q2*v0 + 2*l2*m2*n2*r1*v0 - 2*j2*m2*p2*r1*v0 - 
  2*j2*l2*q2*r1*v0 + 2*i2*p2*q2*r1*v0 + 2*l2*m2*n1*r2*v0 + 
  2*l2*m1*n2*r2*v0 + 2*l1*m2*n2*r2*v0 - 2*j2*m2*p1*r2*v0 - 
  2*j2*m1*p2*r2*v0 - 2*j1*m2*p2*r2*v0 - 2*j2*l2*q1*r2*v0 + 
  2*i2*p2*q1*r2*v0 - 2*j2*l1*q2*r2*v0 - 2*j1*l2*q2*r2*v0 + 
  2*i2*p1*q2*r2*v0 + 2*i1*p2*q2*r2*v0 - 2*k2*m2*n2*s1*v0 + 
  2*j2*m2*o2*s1*v0 + 2*j2*k2*q2*s1*v0 - 2*i2*o2*q2*s1*v0 - 
  2*k2*m2*n1*s2*v0 - 2*k2*m1*n2*s2*v0 - 2*k1*m2*n2*s2*v0 + 
  2*j2*m2*o1*s2*v0 + 2*j2*m1*o2*s2*v0 + 2*j1*m2*o2*s2*v0 + 
  2*j2*k2*q1*s2*v0 - 2*i2*o2*q1*s2*v0 + 2*j2*k1*q2*s2*v0 + 
  2*j1*k2*q2*s2*v0 - 2*i2*o1*q2*s2*v0 - 2*i1*o2*q2*s2*v0 - 
  2*k2*l2*n2*t1*v0 + 2*j2*l2*o2*t1*v0 + 2*j2*k2*p2*t1*v0 - 
  2*i2*o2*p2*t1*v0 - 2*j2*j2*s2*t1*v0 + 2*i2*n2*s2*t1*v0 - 
  2*k2*l2*n1*t2*v0 - 2*k2*l1*n2*t2*v0 - 2*k1*l2*n2*t2*v0 + 
  2*j2*l2*o1*t2*v0 + 2*j2*l1*o2*t2*v0 + 2*j1*l2*o2*t2*v0 + 
  2*j2*k2*p1*t2*v0 - 2*i2*o2*p1*t2*v0 + 2*j2*k1*p2*t2*v0 + 
  2*j1*k2*p2*t2*v0 - 2*i2*o1*p2*t2*v0 - 2*i1*o2*p2*t2*v0 - 
  2*j2*j2*s1*t2*v0 + 2*i2*n2*s1*t2*v0 - 4*j1*j2*s2*t2*v0 + 
  2*i2*n1*s2*t2*v0 + 2*i1*n2*s2*t2*v0 - 4*l2*m2*o0*o2*v1 - 
  2*l2*m0*o2*o2*v1 - 2*l0*m2*o2*o2*v1 + 2*k2*m2*o2*p0*v1 + 
  2*k2*m2*o0*p2*v1 + 2*k2*m0*o2*p2*v1 + 2*k0*m2*o2*p2*v1 + 
  2*k2*l2*o2*q0*v1 - 2*k2*k2*p2*q0*v1 + 2*k2*l2*o0*q2*v1 + 
  2*k2*l0*o2*q2*v1 + 2*k0*l2*o2*q2*v1 - 2*k2*k2*p0*q2*v1 - 
  4*k0*k2*p2*q2*v1 + 2*l2*m2*n2*r0*v1 - 2*j2*m2*p2*r0*v1 - 
  2*j2*l2*q2*r0*v1 + 2*i2*p2*q2*r0*v1 + 2*l2*m2*n0*r2*v1 + 
  2*l2*m0*n2*r2*v1 + 2*l0*m2*n2*r2*v1 - 2*j2*m2*p0*r2*v1 - 
  2*j2*m0*p2*r2*v1 - 2*j0*m2*p2*r2*v1 - 2*j2*l2*q0*r2*v1 + 
  2*i2*p2*q0*r2*v1 - 2*j2*l0*q2*r2*v1 - 2*j0*l2*q2*r2*v1 + 
  2*i2*p0*q2*r2*v1 + 2*i0*p2*q2*r2*v1 - 2*k2*m2*n2*s0*v1 + 
  2*j2*m2*o2*s0*v1 + 2*j2*k2*q2*s0*v1 - 2*i2*o2*q2*s0*v1 - 
  2*k2*m2*n0*s2*v1 - 2*k2*m0*n2*s2*v1 - 2*k0*m2*n2*s2*v1 + 
  2*j2*m2*o0*s2*v1 + 2*j2*m0*o2*s2*v1 + 2*j0*m2*o2*s2*v1 + 
  2*j2*k2*q0*s2*v1 - 2*i2*o2*q0*s2*v1 + 2*j2*k0*q2*s2*v1 + 
  2*j0*k2*q2*s2*v1 - 2*i2*o0*q2*s2*v1 - 2*i0*o2*q2*s2*v1 - 
  2*k2*l2*n2*t0*v1 + 2*j2*l2*o2*t0*v1 + 2*j2*k2*p2*t0*v1 - 
  2*i2*o2*p2*t0*v1 - 2*j2*j2*s2*t0*v1 + 2*i2*n2*s2*t0*v1 - 
  2*k2*l2*n0*t2*v1 - 2*k2*l0*n2*t2*v1 - 2*k0*l2*n2*t2*v1 + 
  2*j2*l2*o0*t2*v1 + 2*j2*l0*o2*t2*v1 + 2*j0*l2*o2*t2*v1 + 
  2*j2*k2*p0*t2*v1 - 2*i2*o2*p0*t2*v1 + 2*j2*k0*p2*t2*v1 + 
  2*j0*k2*p2*t2*v1 - 2*i2*o0*p2*t2*v1 - 2*i0*o2*p2*t2*v1 - 
  2*j2*j2*s0*t2*v1 + 2*i2*n2*s0*t2*v1 - 4*j0*j2*s2*t2*v1 + 
  2*i2*n0*s2*t2*v1 + 2*i0*n2*s2*t2*v1 + 2*k2*k2*n2*v0*v1 - 
  4*j2*k2*o2*v0*v1 + 2*i2*o2*o2*v0*v1 + 2*j2*j2*r2*v0*v1 - 
  2*i2*n2*r2*v0*v1 - 4*l2*m2*o0*o1*v2 - 4*l2*m1*o0*o2*v2 - 
  4*l1*m2*o0*o2*v2 - 4*l2*m0*o1*o2*v2 - 4*l0*m2*o1*o2*v2 - 
  2*l1*m0*o2*o2*v2 - 2*l0*m1*o2*o2*v2 + 2*k2*m2*o1*p0*v2 + 
  2*k2*m1*o2*p0*v2 + 2*k1*m2*o2*p0*v2 + 2*k2*m2*o0*p1*v2 + 
  2*k2*m0*o2*p1*v2 + 2*k0*m2*o2*p1*v2 + 2*k2*m1*o0*p2*v2 + 
  2*k1*m2*o0*p2*v2 + 2*k2*m0*o1*p2*v2 + 2*k0*m2*o1*p2*v2 + 
  2*k1*m0*o2*p2*v2 + 2*k0*m1*o2*p2*v2 + 2*k2*l2*o1*q0*v2 + 
  2*k2*l1*o2*q0*v2 + 2*k1*l2*o2*q0*v2 - 2*k2*k2*p1*q0*v2 - 
  4*k1*k2*p2*q0*v2 + 2*k2*l2*o0*q1*v2 + 2*k2*l0*o2*q1*v2 + 
  2*k0*l2*o2*q1*v2 - 2*k2*k2*p0*q1*v2 - 4*k0*k2*p2*q1*v2 + 
  2*k2*l1*o0*q2*v2 + 2*k1*l2*o0*q2*v2 + 2*k2*l0*o1*q2*v2 + 
  2*k0*l2*o1*q2*v2 + 2*k1*l0*o2*q2*v2 + 2*k0*l1*o2*q2*v2 - 
  4*k1*k2*p0*q2*v2 - 4*k0*k2*p1*q2*v2 - 4*k0*k1*p2*q2*v2 + 
  2*l2*m2*n1*r0*v2 + 2*l2*m1*n2*r0*v2 + 2*l1*m2*n2*r0*v2 - 
  2*j2*m2*p1*r0*v2 - 2*j2*m1*p2*r0*v2 - 2*j1*m2*p2*r0*v2 - 
  2*j2*l2*q1*r0*v2 + 2*i2*p2*q1*r0*v2 - 2*j2*l1*q2*r0*v2 - 
  2*j1*l2*q2*r0*v2 + 2*i2*p1*q2*r0*v2 + 2*i1*p2*q2*r0*v2 + 
  2*l2*m2*n0*r1*v2 + 2*l2*m0*n2*r1*v2 + 2*l0*m2*n2*r1*v2 - 
  2*j2*m2*p0*r1*v2 - 2*j2*m0*p2*r1*v2 - 2*j0*m2*p2*r1*v2 - 
  2*j2*l2*q0*r1*v2 + 2*i2*p2*q0*r1*v2 - 2*j2*l0*q2*r1*v2 - 
  2*j0*l2*q2*r1*v2 + 2*i2*p0*q2*r1*v2 + 2*i0*p2*q2*r1*v2 + 
  2*l2*m1*n0*r2*v2 + 2*l1*m2*n0*r2*v2 + 2*l2*m0*n1*r2*v2 + 
  2*l0*m2*n1*r2*v2 + 2*l1*m0*n2*r2*v2 + 2*l0*m1*n2*r2*v2 - 
  2*j2*m1*p0*r2*v2 - 2*j1*m2*p0*r2*v2 - 2*j2*m0*p1*r2*v2 - 
  2*j0*m2*p1*r2*v2 - 2*j1*m0*p2*r2*v2 - 2*j0*m1*p2*r2*v2 - 
  2*j2*l1*q0*r2*v2 - 2*j1*l2*q0*r2*v2 + 2*i2*p1*q0*r2*v2 + 
  2*i1*p2*q0*r2*v2 - 2*j2*l0*q1*r2*v2 - 2*j0*l2*q1*r2*v2 + 
  2*i2*p0*q1*r2*v2 + 2*i0*p2*q1*r2*v2 - 2*j1*l0*q2*r2*v2 - 
  2*j0*l1*q2*r2*v2 + 2*i1*p0*q2*r2*v2 + 2*i0*p1*q2*r2*v2 - 
  2*k2*m2*n1*s0*v2 - 2*k2*m1*n2*s0*v2 - 2*k1*m2*n2*s0*v2 + 
  2*j2*m2*o1*s0*v2 + 2*j2*m1*o2*s0*v2 + 2*j1*m2*o2*s0*v2 + 
  2*j2*k2*q1*s0*v2 - 2*i2*o2*q1*s0*v2 + 2*j2*k1*q2*s0*v2 + 
  2*j1*k2*q2*s0*v2 - 2*i2*o1*q2*s0*v2 - 2*i1*o2*q2*s0*v2 - 
  2*k2*m2*n0*s1*v2 - 2*k2*m0*n2*s1*v2 - 2*k0*m2*n2*s1*v2 + 
  2*j2*m2*o0*s1*v2 + 2*j2*m0*o2*s1*v2 + 2*j0*m2*o2*s1*v2 + 
  2*j2*k2*q0*s1*v2 - 2*i2*o2*q0*s1*v2 + 2*j2*k0*q2*s1*v2 + 
  2*j0*k2*q2*s1*v2 - 2*i2*o0*q2*s1*v2 - 2*i0*o2*q2*s1*v2 - 
  2*k2*m1*n0*s2*v2 - 2*k1*m2*n0*s2*v2 - 2*k2*m0*n1*s2*v2 - 
  2*k0*m2*n1*s2*v2 - 2*k1*m0*n2*s2*v2 - 2*k0*m1*n2*s2*v2 + 
  2*j2*m1*o0*s2*v2 + 2*j1*m2*o0*s2*v2 + 2*j2*m0*o1*s2*v2 + 
  2*j0*m2*o1*s2*v2 + 2*j1*m0*o2*s2*v2 + 2*j0*m1*o2*s2*v2 + 
  2*j2*k1*q0*s2*v2 + 2*j1*k2*q0*s2*v2 - 2*i2*o1*q0*s2*v2 - 
  2*i1*o2*q0*s2*v2 + 2*j2*k0*q1*s2*v2 + 2*j0*k2*q1*s2*v2 - 
  2*i2*o0*q1*s2*v2 - 2*i0*o2*q1*s2*v2 + 2*j1*k0*q2*s2*v2 + 
  2*j0*k1*q2*s2*v2 - 2*i1*o0*q2*s2*v2 - 2*i0*o1*q2*s2*v2 - 
  2*k2*l2*n1*t0*v2 - 2*k2*l1*n2*t0*v2 - 2*k1*l2*n2*t0*v2 + 
  2*j2*l2*o1*t0*v2 + 2*j2*l1*o2*t0*v2 + 2*j1*l2*o2*t0*v2 + 
  2*j2*k2*p1*t0*v2 - 2*i2*o2*p1*t0*v2 + 2*j2*k1*p2*t0*v2 + 
  2*j1*k2*p2*t0*v2 - 2*i2*o1*p2*t0*v2 - 2*i1*o2*p2*t0*v2 - 
  2*j2*j2*s1*t0*v2 + 2*i2*n2*s1*t0*v2 - 4*j1*j2*s2*t0*v2 + 
  2*i2*n1*s2*t0*v2 + 2*i1*n2*s2*t0*v2 - 2*k2*l2*n0*t1*v2 - 
  2*k2*l0*n2*t1*v2 - 2*k0*l2*n2*t1*v2 + 2*j2*l2*o0*t1*v2 + 
  2*j2*l0*o2*t1*v2 + 2*j0*l2*o2*t1*v2 + 2*j2*k2*p0*t1*v2 - 
  2*i2*o2*p0*t1*v2 + 2*j2*k0*p2*t1*v2 + 2*j0*k2*p2*t1*v2 - 
  2*i2*o0*p2*t1*v2 - 2*i0*o2*p2*t1*v2 - 2*j2*j2*s0*t1*v2 + 
  2*i2*n2*s0*t1*v2 - 4*j0*j2*s2*t1*v2 + 2*i2*n0*s2*t1*v2 + 
  2*i0*n2*s2*t1*v2 - 2*k2*l1*n0*t2*v2 - 2*k1*l2*n0*t2*v2 - 
  2*k2*l0*n1*t2*v2 - 2*k0*l2*n1*t2*v2 - 2*k1*l0*n2*t2*v2 - 
  2*k0*l1*n2*t2*v2 + 2*j2*l1*o0*t2*v2 + 2*j1*l2*o0*t2*v2 + 
  2*j2*l0*o1*t2*v2 + 2*j0*l2*o1*t2*v2 + 2*j1*l0*o2*t2*v2 + 
  2*j0*l1*o2*t2*v2 + 2*j2*k1*p0*t2*v2 + 2*j1*k2*p0*t2*v2 - 
  2*i2*o1*p0*t2*v2 - 2*i1*o2*p0*t2*v2 + 2*j2*k0*p1*t2*v2 + 
  2*j0*k2*p1*t2*v2 - 2*i2*o0*p1*t2*v2 - 2*i0*o2*p1*t2*v2 + 
  2*j1*k0*p2*t2*v2 + 2*j0*k1*p2*t2*v2 - 2*i1*o0*p2*t2*v2 - 
  2*i0*o1*p2*t2*v2 - 4*j1*j2*s0*t2*v2 + 2*i2*n1*s0*t2*v2 + 
  2*i1*n2*s0*t2*v2 - 4*j0*j2*s1*t2*v2 + 2*i2*n0*s1*t2*v2 + 
  2*i0*n2*s1*t2*v2 - 4*j0*j1*s2*t2*v2 + 2*i1*n0*s2*t2*v2 + 
  2*i0*n1*s2*t2*v2 + 2*k2*k2*n1*v0*v2 + 4*k1*k2*n2*v0*v2 - 
  4*j2*k2*o1*v0*v2 - 4*j2*k1*o2*v0*v2 - 4*j1*k2*o2*v0*v2 + 
  4*i2*o1*o2*v0*v2 + 2*i1*o2*o2*v0*v2 + 2*j2*j2*r1*v0*v2 - 
  2*i2*n2*r1*v0*v2 + 4*j1*j2*r2*v0*v2 - 2*i2*n1*r2*v0*v2 - 
  2*i1*n2*r2*v0*v2 + 2*k2*k2*n0*v1*v2 + 4*k0*k2*n2*v1*v2 - 
  4*j2*k2*o0*v1*v2 - 4*j2*k0*o2*v1*v2 - 4*j0*k2*o2*v1*v2 + 
  4*i2*o0*o2*v1*v2 + 2*i0*o2*o2*v1*v2 + 2*j2*j2*r0*v1*v2 - 
  2*i2*n2*r0*v1*v2 + 4*j0*j2*r2*v1*v2 - 2*i2*n0*r2*v1*v2 - 
  2*i0*n2*r2*v1*v2 + 2*k1*k2*n0*v2*v2 + 2*k0*k2*n1*v2*v2 + 
  2*k0*k1*n2*v2*v2 - 2*j2*k1*o0*v2*v2 - 2*j1*k2*o0*v2*v2 - 
  2*j2*k0*o1*v2*v2 - 2*j0*k2*o1*v2*v2 + 2*i2*o0*o1*v2*v2 - 
  2*j1*k0*o2*v2*v2 - 2*j0*k1*o2*v2*v2 + 2*i1*o0*o2*v2*v2 + 
  2*i0*o1*o2*v2*v2 + 2*j1*j2*r0*v2*v2 - i2*n1*r0*v2*v2 - 
  i1*n2*r0*v2*v2 + 2*j0*j2*r1*v2*v2 - i2*n0*r1*v2*v2 - 
  i0*n2*r1*v2*v2 + 2*j0*j1*r2*v2*v2 - i1*n0*r2*v2*v2 - 
  i0*n1*r2*v2*v2 + 2*l2*l2*o1*o2*w0 + 2*l1*l2*o2*o2*w0 - 
  2*k2*l2*o2*p1*w0 - 2*k2*l2*o1*p2*w0 - 2*k2*l1*o2*p2*w0 - 
  2*k1*l2*o2*p2*w0 + 2*k2*k2*p1*p2*w0 + 2*k1*k2*p2*p2*w0 - 
  l2*l2*n2*r1*w0 + 2*j2*l2*p2*r1*w0 - i2*p2*p2*r1*w0 - 
  l2*l2*n1*r2*w0 - 2*l1*l2*n2*r2*w0 + 2*j2*l2*p1*r2*w0 + 
  2*j2*l1*p2*r2*w0 + 2*j1*l2*p2*r2*w0 - 2*i2*p1*p2*r2*w0 - 
  i1*p2*p2*r2*w0 + 2*k2*l2*n2*s1*w0 - 2*j2*l2*o2*s1*w0 - 
  2*j2*k2*p2*s1*w0 + 2*i2*o2*p2*s1*w0 + 2*k2*l2*n1*s2*w0 + 
  2*k2*l1*n2*s2*w0 + 2*k1*l2*n2*s2*w0 - 2*j2*l2*o1*s2*w0 - 
  2*j2*l1*o2*s2*w0 - 2*j1*l2*o2*s2*w0 - 2*j2*k2*p1*s2*w0 + 
  2*i2*o2*p1*s2*w0 - 2*j2*k1*p2*s2*w0 - 2*j1*k2*p2*s2*w0 + 
  2*i2*o1*p2*s2*w0 + 2*i1*o2*p2*s2*w0 + 2*j2*j2*s1*s2*w0 - 
  2*i2*n2*s1*s2*w0 + 2*j1*j2*s2*s2*w0 - i2*n1*s2*s2*w0 - 
  i1*n2*s2*s2*w0 - k2*k2*n2*u1*w0 + 2*j2*k2*o2*u1*w0 - 
  i2*o2*o2*u1*w0 - j2*j2*r2*u1*w0 + i2*n2*r2*u1*w0 - 
  k2*k2*n1*u2*w0 - 2*k1*k2*n2*u2*w0 + 2*j2*k2*o1*u2*w0 + 
  2*j2*k1*o2*u2*w0 + 2*j1*k2*o2*u2*w0 - 2*i2*o1*o2*u2*w0 - 
  i1*o2*o2*u2*w0 - j2*j2*r1*u2*w0 + i2*n2*r1*u2*w0 - 
  2*j1*j2*r2*u2*w0 + i2*n1*r2*u2*w0 + i1*n2*r2*u2*w0 + 
  2*l2*l2*o0*o2*w1 + 2*l0*l2*o2*o2*w1 - 2*k2*l2*o2*p0*w1 - 
  2*k2*l2*o0*p2*w1 - 2*k2*l0*o2*p2*w1 - 2*k0*l2*o2*p2*w1 + 
  2*k2*k2*p0*p2*w1 + 2*k0*k2*p2*p2*w1 - l2*l2*n2*r0*w1 + 
  2*j2*l2*p2*r0*w1 - i2*p2*p2*r0*w1 - l2*l2*n0*r2*w1 - 
  2*l0*l2*n2*r2*w1 + 2*j2*l2*p0*r2*w1 + 2*j2*l0*p2*r2*w1 + 
  2*j0*l2*p2*r2*w1 - 2*i2*p0*p2*r2*w1 - i0*p2*p2*r2*w1 + 
  2*k2*l2*n2*s0*w1 - 2*j2*l2*o2*s0*w1 - 2*j2*k2*p2*s0*w1 + 
  2*i2*o2*p2*s0*w1 + 2*k2*l2*n0*s2*w1 + 2*k2*l0*n2*s2*w1 + 
  2*k0*l2*n2*s2*w1 - 2*j2*l2*o0*s2*w1 - 2*j2*l0*o2*s2*w1 - 
  2*j0*l2*o2*s2*w1 - 2*j2*k2*p0*s2*w1 + 2*i2*o2*p0*s2*w1 - 
  2*j2*k0*p2*s2*w1 - 2*j0*k2*p2*s2*w1 + 2*i2*o0*p2*s2*w1 + 
  2*i0*o2*p2*s2*w1 + 2*j2*j2*s0*s2*w1 - 2*i2*n2*s0*s2*w1 + 
  2*j0*j2*s2*s2*w1 - i2*n0*s2*s2*w1 - i0*n2*s2*s2*w1 - 
  k2*k2*n2*u0*w1 + 2*j2*k2*o2*u0*w1 - i2*o2*o2*u0*w1 - 
  j2*j2*r2*u0*w1 + i2*n2*r2*u0*w1 - k2*k2*n0*u2*w1 - 
  2*k0*k2*n2*u2*w1 + 2*j2*k2*o0*u2*w1 + 2*j2*k0*o2*u2*w1 + 
  2*j0*k2*o2*u2*w1 - 2*i2*o0*o2*u2*w1 - i0*o2*o2*u2*w1 - 
  j2*j2*r0*u2*w1 + i2*n2*r0*u2*w1 - 2*j0*j2*r2*u2*w1 + 
  i2*n0*r2*u2*w1 + i0*n2*r2*u2*w1 + 2*l2*l2*o0*o1*w2 + 
  4*l1*l2*o0*o2*w2 + 4*l0*l2*o1*o2*w2 + 2*l0*l1*o2*o2*w2 - 
  2*k2*l2*o1*p0*w2 - 2*k2*l1*o2*p0*w2 - 2*k1*l2*o2*p0*w2 - 
  2*k2*l2*o0*p1*w2 - 2*k2*l0*o2*p1*w2 - 2*k0*l2*o2*p1*w2 + 
  2*k2*k2*p0*p1*w2 - 2*k2*l1*o0*p2*w2 - 2*k1*l2*o0*p2*w2 - 
  2*k2*l0*o1*p2*w2 - 2*k0*l2*o1*p2*w2 - 2*k1*l0*o2*p2*w2 - 
  2*k0*l1*o2*p2*w2 + 4*k1*k2*p0*p2*w2 + 4*k0*k2*p1*p2*w2 + 
  2*k0*k1*p2*p2*w2 - l2*l2*n1*r0*w2 - 2*l1*l2*n2*r0*w2 + 
  2*j2*l2*p1*r0*w2 + 2*j2*l1*p2*r0*w2 + 2*j1*l2*p2*r0*w2 - 
  2*i2*p1*p2*r0*w2 - i1*p2*p2*r0*w2 - l2*l2*n0*r1*w2 - 
  2*l0*l2*n2*r1*w2 + 2*j2*l2*p0*r1*w2 + 2*j2*l0*p2*r1*w2 + 
  2*j0*l2*p2*r1*w2 - 2*i2*p0*p2*r1*w2 - i0*p2*p2*r1*w2 - 
  2*l1*l2*n0*r2*w2 - 2*l0*l2*n1*r2*w2 - 2*l0*l1*n2*r2*w2 + 
  2*j2*l1*p0*r2*w2 + 2*j1*l2*p0*r2*w2 + 2*j2*l0*p1*r2*w2 + 
  2*j0*l2*p1*r2*w2 - 2*i2*p0*p1*r2*w2 + 2*j1*l0*p2*r2*w2 + 
  2*j0*l1*p2*r2*w2 - 2*i1*p0*p2*r2*w2 - 2*i0*p1*p2*r2*w2 + 
  2*k2*l2*n1*s0*w2 + 2*k2*l1*n2*s0*w2 + 2*k1*l2*n2*s0*w2 - 
  2*j2*l2*o1*s0*w2 - 2*j2*l1*o2*s0*w2 - 2*j1*l2*o2*s0*w2 - 
  2*j2*k2*p1*s0*w2 + 2*i2*o2*p1*s0*w2 - 2*j2*k1*p2*s0*w2 - 
  2*j1*k2*p2*s0*w2 + 2*i2*o1*p2*s0*w2 + 2*i1*o2*p2*s0*w2 + 
  2*k2*l2*n0*s1*w2 + 2*k2*l0*n2*s1*w2 + 2*k0*l2*n2*s1*w2 - 
  2*j2*l2*o0*s1*w2 - 2*j2*l0*o2*s1*w2 - 2*j0*l2*o2*s1*w2 - 
  2*j2*k2*p0*s1*w2 + 2*i2*o2*p0*s1*w2 - 2*j2*k0*p2*s1*w2 - 
  2*j0*k2*p2*s1*w2 + 2*i2*o0*p2*s1*w2 + 2*i0*o2*p2*s1*w2 + 
  2*j2*j2*s0*s1*w2 - 2*i2*n2*s0*s1*w2 + 2*k2*l1*n0*s2*w2 + 
  2*k1*l2*n0*s2*w2 + 2*k2*l0*n1*s2*w2 + 2*k0*l2*n1*s2*w2 + 
  2*k1*l0*n2*s2*w2 + 2*k0*l1*n2*s2*w2 - 2*j2*l1*o0*s2*w2 - 
  2*j1*l2*o0*s2*w2 - 2*j2*l0*o1*s2*w2 - 2*j0*l2*o1*s2*w2 - 
  2*j1*l0*o2*s2*w2 - 2*j0*l1*o2*s2*w2 - 2*j2*k1*p0*s2*w2 - 
  2*j1*k2*p0*s2*w2 + 2*i2*o1*p0*s2*w2 + 2*i1*o2*p0*s2*w2 - 
  2*j2*k0*p1*s2*w2 - 2*j0*k2*p1*s2*w2 + 2*i2*o0*p1*s2*w2 + 
  2*i0*o2*p1*s2*w2 - 2*j1*k0*p2*s2*w2 - 2*j0*k1*p2*s2*w2 + 
  2*i1*o0*p2*s2*w2 + 2*i0*o1*p2*s2*w2 + 4*j1*j2*s0*s2*w2 - 
  2*i2*n1*s0*s2*w2 - 2*i1*n2*s0*s2*w2 + 4*j0*j2*s1*s2*w2 - 
  2*i2*n0*s1*s2*w2 - 2*i0*n2*s1*s2*w2 + 2*j0*j1*s2*s2*w2 - 
  i1*n0*s2*s2*w2 - i0*n1*s2*s2*w2 - k2*k2*n1*u0*w2 - 
  2*k1*k2*n2*u0*w2 + 2*j2*k2*o1*u0*w2 + 2*j2*k1*o2*u0*w2 + 
  2*j1*k2*o2*u0*w2 - 2*i2*o1*o2*u0*w2 - i1*o2*o2*u0*w2 - 
  j2*j2*r1*u0*w2 + i2*n2*r1*u0*w2 - 2*j1*j2*r2*u0*w2 + 
  i2*n1*r2*u0*w2 + i1*n2*r2*u0*w2 - k2*k2*n0*u1*w2 - 
  2*k0*k2*n2*u1*w2 + 2*j2*k2*o0*u1*w2 + 2*j2*k0*o2*u1*w2 + 
  2*j0*k2*o2*u1*w2 - 2*i2*o0*o2*u1*w2 - i0*o2*o2*u1*w2 - 
  j2*j2*r0*u1*w2 + i2*n2*r0*u1*w2 - 2*j0*j2*r2*u1*w2 + 
  i2*n0*r2*u1*w2 + i0*n2*r2*u1*w2 - 2*k1*k2*n0*u2*w2 - 
  2*k0*k2*n1*u2*w2 - 2*k0*k1*n2*u2*w2 + 2*j2*k1*o0*u2*w2 + 
  2*j1*k2*o0*u2*w2 + 2*j2*k0*o1*u2*w2 + 2*j0*k2*o1*u2*w2 - 
  2*i2*o0*o1*u2*w2 + 2*j1*k0*o2*u2*w2 + 2*j0*k1*o2*u2*w2 - 
  2*i1*o0*o2*u2*w2 - 2*i0*o1*o2*u2*w2 - 2*j1*j2*r0*u2*w2 + 
  i2*n1*r0*u2*w2 + i1*n2*r0*u2*w2 - 2*j0*j2*r1*u2*w2 + 
  i2*n0*r1*u2*w2 + i0*n2*r1*u2*w2 - 2*j0*j1*r2*u2*w2 + 
  i1*n0*r2*u2*w2 + i0*n1*r2*u2*w2;
    
    k[5]  = 2*m2*m2*p1*p2*r1 + 
  2*m1*m2*p2*p2*r1 - 2*l2*m2*p2*q1*r1 - 2*l2*m2*p1*q2*r1 - 
  2*l2*m1*p2*q2*r1 - 2*l1*m2*p2*q2*r1 + 2*l2*l2*q1*q2*r1 + 
  2*l1*l2*q2*q2*r1 + m2*m2*p1*p1*r2 + 4*m1*m2*p1*p2*r2 + 
  m1*m1*p2*p2*r2 - 2*l2*m2*p1*q1*r2 - 2*l2*m1*p2*q1*r2 - 
  2*l1*m2*p2*q1*r2 + l2*l2*q1*q1*r2 - 2*l2*m1*p1*q2*r2 - 
  2*l1*m2*p1*q2*r2 - 2*l1*m1*p2*q2*r2 + 4*l1*l2*q1*q2*r2 + 
  l1*l1*q2*q2*r2 - 2*m2*m2*o2*p1*s1 - 2*m2*m2*o1*p2*s1 - 
  4*m1*m2*o2*p2*s1 + 2*l2*m2*o2*q1*s1 + 2*k2*m2*p2*q1*s1 + 
  2*l2*m2*o1*q2*s1 + 2*l2*m1*o2*q2*s1 + 2*l1*m2*o2*q2*s1 + 
  2*k2*m2*p1*q2*s1 + 2*k2*m1*p2*q2*s1 + 2*k1*m2*p2*q2*s1 - 
  4*k2*l2*q1*q2*s1 - 2*k2*l1*q2*q2*s1 - 2*k1*l2*q2*q2*s1 + 
  m2*m2*n2*s1*s1 - 2*j2*m2*q2*s1*s1 + i2*q2*q2*s1*s1 - 
  2*m2*m2*o1*p1*s2 - 4*m1*m2*o2*p1*s2 - 4*m1*m2*o1*p2*s2 - 
  2*m1*m1*o2*p2*s2 + 2*l2*m2*o1*q1*s2 + 2*l2*m1*o2*q1*s2 + 
  2*l1*m2*o2*q1*s2 + 2*k2*m2*p1*q1*s2 + 2*k2*m1*p2*q1*s2 + 
  2*k1*m2*p2*q1*s2 - 2*k2*l2*q1*q1*s2 + 2*l2*m1*o1*q2*s2 + 
  2*l1*m2*o1*q2*s2 + 2*l1*m1*o2*q2*s2 + 2*k2*m1*p1*q2*s2 + 
  2*k1*m2*p1*q2*s2 + 2*k1*m1*p2*q2*s2 - 4*k2*l1*q1*q2*s2 - 
  4*k1*l2*q1*q2*s2 - 2*k1*l1*q2*q2*s2 + 2*m2*m2*n1*s1*s2 + 
  4*m1*m2*n2*s1*s2 - 4*j2*m2*q1*s1*s2 - 4*j2*m1*q2*s1*s2 - 
  4*j1*m2*q2*s1*s2 + 4*i2*q1*q2*s1*s2 + 2*i1*q2*q2*s1*s2 + 
  2*m1*m2*n1*s2*s2 + m1*m1*n2*s2*s2 - 2*j2*m1*q1*s2*s2 - 
  2*j1*m2*q1*s2*s2 + i2*q1*q1*s2*s2 - 2*j1*m1*q2*s2*s2 + 
  2*i1*q1*q2*s2*s2 + 2*l2*m2*o2*p1*t1 + 2*l2*m2*o1*p2*t1 + 
  2*l2*m1*o2*p2*t1 + 2*l1*m2*o2*p2*t1 - 4*k2*m2*p1*p2*t1 - 
  2*k2*m1*p2*p2*t1 - 2*k1*m2*p2*p2*t1 - 2*l2*l2*o2*q1*t1 + 
  2*k2*l2*p2*q1*t1 - 2*l2*l2*o1*q2*t1 - 4*l1*l2*o2*q2*t1 + 
  2*k2*l2*p1*q2*t1 + 2*k2*l1*p2*q2*t1 + 2*k1*l2*p2*q2*t1 - 
  2*l2*m2*n2*s1*t1 + 2*j2*m2*p2*s1*t1 + 2*j2*l2*q2*s1*t1 - 
  2*i2*p2*q2*s1*t1 - 2*l2*m2*n1*s2*t1 - 2*l2*m1*n2*s2*t1 - 
  2*l1*m2*n2*s2*t1 + 2*j2*m2*p1*s2*t1 + 2*j2*m1*p2*s2*t1 + 
  2*j1*m2*p2*s2*t1 + 2*j2*l2*q1*s2*t1 - 2*i2*p2*q1*s2*t1 + 
  2*j2*l1*q2*s2*t1 + 2*j1*l2*q2*s2*t1 - 2*i2*p1*q2*s2*t1 - 
  2*i1*p2*q2*s2*t1 + l2*l2*n2*t1*t1 - 2*j2*l2*p2*t1*t1 + 
  i2*p2*p2*t1*t1 + 2*l2*m2*o1*p1*t2 + 2*l2*m1*o2*p1*t2 + 
  2*l1*m2*o2*p1*t2 - 2*k2*m2*p1*p1*t2 + 2*l2*m1*o1*p2*t2 + 
  2*l1*m2*o1*p2*t2 + 2*l1*m1*o2*p2*t2 - 4*k2*m1*p1*p2*t2 - 
  4*k1*m2*p1*p2*t2 - 2*k1*m1*p2*p2*t2 - 2*l2*l2*o1*q1*t2 - 
  4*l1*l2*o2*q1*t2 + 2*k2*l2*p1*q1*t2 + 2*k2*l1*p2*q1*t2 + 
  2*k1*l2*p2*q1*t2 - 4*l1*l2*o1*q2*t2 - 2*l1*l1*o2*q2*t2 + 
  2*k2*l1*p1*q2*t2 + 2*k1*l2*p1*q2*t2 + 2*k1*l1*p2*q2*t2 - 
  2*l2*m2*n1*s1*t2 - 2*l2*m1*n2*s1*t2 - 2*l1*m2*n2*s1*t2 + 
  2*j2*m2*p1*s1*t2 + 2*j2*m1*p2*s1*t2 + 2*j1*m2*p2*s1*t2 + 
  2*j2*l2*q1*s1*t2 - 2*i2*p2*q1*s1*t2 + 2*j2*l1*q2*s1*t2 + 
  2*j1*l2*q2*s1*t2 - 2*i2*p1*q2*s1*t2 - 2*i1*p2*q2*s1*t2 - 
  2*l2*m1*n1*s2*t2 - 2*l1*m2*n1*s2*t2 - 2*l1*m1*n2*s2*t2 + 
  2*j2*m1*p1*s2*t2 + 2*j1*m2*p1*s2*t2 + 2*j1*m1*p2*s2*t2 + 
  2*j2*l1*q1*s2*t2 + 2*j1*l2*q1*s2*t2 - 2*i2*p1*q1*s2*t2 - 
  2*i1*p2*q1*s2*t2 + 2*j1*l1*q2*s2*t2 - 2*i1*p1*q2*s2*t2 + 
  2*l2*l2*n1*t1*t2 + 4*l1*l2*n2*t1*t2 - 4*j2*l2*p1*t1*t2 - 
  4*j2*l1*p2*t1*t2 - 4*j1*l2*p2*t1*t2 + 4*i2*p1*p2*t1*t2 + 
  2*i1*p2*p2*t1*t2 + 2*l1*l2*n1*t2*t2 + l1*l1*n2*t2*t2 - 
  2*j2*l1*p1*t2*t2 - 2*j1*l2*p1*t2*t2 + i2*p1*p1*t2*t2 - 
  2*j1*l1*p2*t2*t2 + 2*i1*p1*p2*t2*t2 + 2*m2*m2*o1*o2*u1 + 
  2*m1*m2*o2*o2*u1 - 2*k2*m2*o2*q1*u1 - 2*k2*m2*o1*q2*u1 - 
  2*k2*m1*o2*q2*u1 - 2*k1*m2*o2*q2*u1 + 2*k2*k2*q1*q2*u1 + 
  2*k1*k2*q2*q2*u1 - m2*m2*n2*r1*u1 + 2*j2*m2*q2*r1*u1 - 
  i2*q2*q2*r1*u1 - m2*m2*n1*r2*u1 - 2*m1*m2*n2*r2*u1 + 
  2*j2*m2*q1*r2*u1 + 2*j2*m1*q2*r2*u1 + 2*j1*m2*q2*r2*u1 - 
  2*i2*q1*q2*r2*u1 - i1*q2*q2*r2*u1 + 2*k2*m2*n2*t1*u1 - 
  2*j2*m2*o2*t1*u1 - 2*j2*k2*q2*t1*u1 + 2*i2*o2*q2*t1*u1 + 
  2*k2*m2*n1*t2*u1 + 2*k2*m1*n2*t2*u1 + 2*k1*m2*n2*t2*u1 - 
  2*j2*m2*o1*t2*u1 - 2*j2*m1*o2*t2*u1 - 2*j1*m2*o2*t2*u1 - 
  2*j2*k2*q1*t2*u1 + 2*i2*o2*q1*t2*u1 - 2*j2*k1*q2*t2*u1 - 
  2*j1*k2*q2*t2*u1 + 2*i2*o1*q2*t2*u1 + 2*i1*o2*q2*t2*u1 + 
  2*j2*j2*t1*t2*u1 - 2*i2*n2*t1*t2*u1 + 2*j1*j2*t2*t2*u1 - 
  i2*n1*t2*t2*u1 - i1*n2*t2*t2*u1 + m2*m2*o1*o1*u2 + 
  4*m1*m2*o1*o2*u2 + m1*m1*o2*o2*u2 - 2*k2*m2*o1*q1*u2 - 
  2*k2*m1*o2*q1*u2 - 2*k1*m2*o2*q1*u2 + k2*k2*q1*q1*u2 - 
  2*k2*m1*o1*q2*u2 - 2*k1*m2*o1*q2*u2 - 2*k1*m1*o2*q2*u2 + 
  4*k1*k2*q1*q2*u2 + k1*k1*q2*q2*u2 - m2*m2*n1*r1*u2 - 
  2*m1*m2*n2*r1*u2 + 2*j2*m2*q1*r1*u2 + 2*j2*m1*q2*r1*u2 + 
  2*j1*m2*q2*r1*u2 - 2*i2*q1*q2*r1*u2 - i1*q2*q2*r1*u2 - 
  2*m1*m2*n1*r2*u2 - m1*m1*n2*r2*u2 + 2*j2*m1*q1*r2*u2 + 
  2*j1*m2*q1*r2*u2 - i2*q1*q1*r2*u2 + 2*j1*m1*q2*r2*u2 - 
  2*i1*q1*q2*r2*u2 + 2*k2*m2*n1*t1*u2 + 2*k2*m1*n2*t1*u2 + 
  2*k1*m2*n2*t1*u2 - 2*j2*m2*o1*t1*u2 - 2*j2*m1*o2*t1*u2 - 
  2*j1*m2*o2*t1*u2 - 2*j2*k2*q1*t1*u2 + 2*i2*o2*q1*t1*u2 - 
  2*j2*k1*q2*t1*u2 - 2*j1*k2*q2*t1*u2 + 2*i2*o1*q2*t1*u2 + 
  2*i1*o2*q2*t1*u2 + j2*j2*t1*t1*u2 - i2*n2*t1*t1*u2 + 
  2*k2*m1*n1*t2*u2 + 2*k1*m2*n1*t2*u2 + 2*k1*m1*n2*t2*u2 - 
  2*j2*m1*o1*t2*u2 - 2*j1*m2*o1*t2*u2 - 2*j1*m1*o2*t2*u2 - 
  2*j2*k1*q1*t2*u2 - 2*j1*k2*q1*t2*u2 + 2*i2*o1*q1*t2*u2 + 
  2*i1*o2*q1*t2*u2 - 2*j1*k1*q2*t2*u2 + 2*i1*o1*q2*t2*u2 + 
  4*j1*j2*t1*t2*u2 - 2*i2*n1*t1*t2*u2 - 2*i1*n2*t1*t2*u2 + 
  j1*j1*t2*t2*u2 - i1*n1*t2*t2*u2 - 4*l2*m2*o1*o2*v1 - 
  2*l2*m1*o2*o2*v1 - 2*l1*m2*o2*o2*v1 + 2*k2*m2*o2*p1*v1 + 
  2*k2*m2*o1*p2*v1 + 2*k2*m1*o2*p2*v1 + 2*k1*m2*o2*p2*v1 + 
  2*k2*l2*o2*q1*v1 - 2*k2*k2*p2*q1*v1 + 2*k2*l2*o1*q2*v1 + 
  2*k2*l1*o2*q2*v1 + 2*k1*l2*o2*q2*v1 - 2*k2*k2*p1*q2*v1 - 
  4*k1*k2*p2*q2*v1 + 2*l2*m2*n2*r1*v1 - 2*j2*m2*p2*r1*v1 - 
  2*j2*l2*q2*r1*v1 + 2*i2*p2*q2*r1*v1 + 2*l2*m2*n1*r2*v1 + 
  2*l2*m1*n2*r2*v1 + 2*l1*m2*n2*r2*v1 - 2*j2*m2*p1*r2*v1 - 
  2*j2*m1*p2*r2*v1 - 2*j1*m2*p2*r2*v1 - 2*j2*l2*q1*r2*v1 + 
  2*i2*p2*q1*r2*v1 - 2*j2*l1*q2*r2*v1 - 2*j1*l2*q2*r2*v1 + 
  2*i2*p1*q2*r2*v1 + 2*i1*p2*q2*r2*v1 - 2*k2*m2*n2*s1*v1 + 
  2*j2*m2*o2*s1*v1 + 2*j2*k2*q2*s1*v1 - 2*i2*o2*q2*s1*v1 - 
  2*k2*m2*n1*s2*v1 - 2*k2*m1*n2*s2*v1 - 2*k1*m2*n2*s2*v1 + 
  2*j2*m2*o1*s2*v1 + 2*j2*m1*o2*s2*v1 + 2*j1*m2*o2*s2*v1 + 
  2*j2*k2*q1*s2*v1 - 2*i2*o2*q1*s2*v1 + 2*j2*k1*q2*s2*v1 + 
  2*j1*k2*q2*s2*v1 - 2*i2*o1*q2*s2*v1 - 2*i1*o2*q2*s2*v1 - 
  2*k2*l2*n2*t1*v1 + 2*j2*l2*o2*t1*v1 + 2*j2*k2*p2*t1*v1 - 
  2*i2*o2*p2*t1*v1 - 2*j2*j2*s2*t1*v1 + 2*i2*n2*s2*t1*v1 - 
  2*k2*l2*n1*t2*v1 - 2*k2*l1*n2*t2*v1 - 2*k1*l2*n2*t2*v1 + 
  2*j2*l2*o1*t2*v1 + 2*j2*l1*o2*t2*v1 + 2*j1*l2*o2*t2*v1 + 
  2*j2*k2*p1*t2*v1 - 2*i2*o2*p1*t2*v1 + 2*j2*k1*p2*t2*v1 + 
  2*j1*k2*p2*t2*v1 - 2*i2*o1*p2*t2*v1 - 2*i1*o2*p2*t2*v1 - 
  2*j2*j2*s1*t2*v1 + 2*i2*n2*s1*t2*v1 - 4*j1*j2*s2*t2*v1 + 
  2*i2*n1*s2*t2*v1 + 2*i1*n2*s2*t2*v1 + k2*k2*n2*v1*v1 - 
  2*j2*k2*o2*v1*v1 + i2*o2*o2*v1*v1 + j2*j2*r2*v1*v1 - 
  i2*n2*r2*v1*v1 - 2*l2*m2*o1*o1*v2 - 4*l2*m1*o1*o2*v2 - 
  4*l1*m2*o1*o2*v2 - 2*l1*m1*o2*o2*v2 + 2*k2*m2*o1*p1*v2 + 
  2*k2*m1*o2*p1*v2 + 2*k1*m2*o2*p1*v2 + 2*k2*m1*o1*p2*v2 + 
  2*k1*m2*o1*p2*v2 + 2*k1*m1*o2*p2*v2 + 2*k2*l2*o1*q1*v2 + 
  2*k2*l1*o2*q1*v2 + 2*k1*l2*o2*q1*v2 - 2*k2*k2*p1*q1*v2 - 
  4*k1*k2*p2*q1*v2 + 2*k2*l1*o1*q2*v2 + 2*k1*l2*o1*q2*v2 + 
  2*k1*l1*o2*q2*v2 - 4*k1*k2*p1*q2*v2 - 2*k1*k1*p2*q2*v2 + 
  2*l2*m2*n1*r1*v2 + 2*l2*m1*n2*r1*v2 + 2*l1*m2*n2*r1*v2 - 
  2*j2*m2*p1*r1*v2 - 2*j2*m1*p2*r1*v2 - 2*j1*m2*p2*r1*v2 - 
  2*j2*l2*q1*r1*v2 + 2*i2*p2*q1*r1*v2 - 2*j2*l1*q2*r1*v2 - 
  2*j1*l2*q2*r1*v2 + 2*i2*p1*q2*r1*v2 + 2*i1*p2*q2*r1*v2 + 
  2*l2*m1*n1*r2*v2 + 2*l1*m2*n1*r2*v2 + 2*l1*m1*n2*r2*v2 - 
  2*j2*m1*p1*r2*v2 - 2*j1*m2*p1*r2*v2 - 2*j1*m1*p2*r2*v2 - 
  2*j2*l1*q1*r2*v2 - 2*j1*l2*q1*r2*v2 + 2*i2*p1*q1*r2*v2 + 
  2*i1*p2*q1*r2*v2 - 2*j1*l1*q2*r2*v2 + 2*i1*p1*q2*r2*v2 - 
  2*k2*m2*n1*s1*v2 - 2*k2*m1*n2*s1*v2 - 2*k1*m2*n2*s1*v2 + 
  2*j2*m2*o1*s1*v2 + 2*j2*m1*o2*s1*v2 + 2*j1*m2*o2*s1*v2 + 
  2*j2*k2*q1*s1*v2 - 2*i2*o2*q1*s1*v2 + 2*j2*k1*q2*s1*v2 + 
  2*j1*k2*q2*s1*v2 - 2*i2*o1*q2*s1*v2 - 2*i1*o2*q2*s1*v2 - 
  2*k2*m1*n1*s2*v2 - 2*k1*m2*n1*s2*v2 - 2*k1*m1*n2*s2*v2 + 
  2*j2*m1*o1*s2*v2 + 2*j1*m2*o1*s2*v2 + 2*j1*m1*o2*s2*v2 + 
  2*j2*k1*q1*s2*v2 + 2*j1*k2*q1*s2*v2 - 2*i2*o1*q1*s2*v2 - 
  2*i1*o2*q1*s2*v2 + 2*j1*k1*q2*s2*v2 - 2*i1*o1*q2*s2*v2 - 
  2*k2*l2*n1*t1*v2 - 2*k2*l1*n2*t1*v2 - 2*k1*l2*n2*t1*v2 + 
  2*j2*l2*o1*t1*v2 + 2*j2*l1*o2*t1*v2 + 2*j1*l2*o2*t1*v2 + 
  2*j2*k2*p1*t1*v2 - 2*i2*o2*p1*t1*v2 + 2*j2*k1*p2*t1*v2 + 
  2*j1*k2*p2*t1*v2 - 2*i2*o1*p2*t1*v2 - 2*i1*o2*p2*t1*v2 - 
  2*j2*j2*s1*t1*v2 + 2*i2*n2*s1*t1*v2 - 4*j1*j2*s2*t1*v2 + 
  2*i2*n1*s2*t1*v2 + 2*i1*n2*s2*t1*v2 - 2*k2*l1*n1*t2*v2 - 
  2*k1*l2*n1*t2*v2 - 2*k1*l1*n2*t2*v2 + 2*j2*l1*o1*t2*v2 + 
  2*j1*l2*o1*t2*v2 + 2*j1*l1*o2*t2*v2 + 2*j2*k1*p1*t2*v2 + 
  2*j1*k2*p1*t2*v2 - 2*i2*o1*p1*t2*v2 - 2*i1*o2*p1*t2*v2 + 
  2*j1*k1*p2*t2*v2 - 2*i1*o1*p2*t2*v2 - 4*j1*j2*s1*t2*v2 + 
  2*i2*n1*s1*t2*v2 + 2*i1*n2*s1*t2*v2 - 2*j1*j1*s2*t2*v2 + 
  2*i1*n1*s2*t2*v2 + 2*k2*k2*n1*v1*v2 + 4*k1*k2*n2*v1*v2 - 
  4*j2*k2*o1*v1*v2 - 4*j2*k1*o2*v1*v2 - 4*j1*k2*o2*v1*v2 + 
  4*i2*o1*o2*v1*v2 + 2*i1*o2*o2*v1*v2 + 2*j2*j2*r1*v1*v2 - 
  2*i2*n2*r1*v1*v2 + 4*j1*j2*r2*v1*v2 - 2*i2*n1*r2*v1*v2 - 
  2*i1*n2*r2*v1*v2 + 2*k1*k2*n1*v2*v2 + k1*k1*n2*v2*v2 - 
  2*j2*k1*o1*v2*v2 - 2*j1*k2*o1*v2*v2 + i2*o1*o1*v2*v2 - 
  2*j1*k1*o2*v2*v2 + 2*i1*o1*o2*v2*v2 + 2*j1*j2*r1*v2*v2 - 
  i2*n1*r1*v2*v2 - i1*n2*r1*v2*v2 + j1*j1*r2*v2*v2 - 
  i1*n1*r2*v2*v2 + 2*l2*l2*o1*o2*w1 + 2*l1*l2*o2*o2*w1 - 
  2*k2*l2*o2*p1*w1 - 2*k2*l2*o1*p2*w1 - 2*k2*l1*o2*p2*w1 - 
  2*k1*l2*o2*p2*w1 + 2*k2*k2*p1*p2*w1 + 2*k1*k2*p2*p2*w1 - 
  l2*l2*n2*r1*w1 + 2*j2*l2*p2*r1*w1 - i2*p2*p2*r1*w1 - 
  l2*l2*n1*r2*w1 - 2*l1*l2*n2*r2*w1 + 2*j2*l2*p1*r2*w1 + 
  2*j2*l1*p2*r2*w1 + 2*j1*l2*p2*r2*w1 - 2*i2*p1*p2*r2*w1 - 
  i1*p2*p2*r2*w1 + 2*k2*l2*n2*s1*w1 - 2*j2*l2*o2*s1*w1 - 
  2*j2*k2*p2*s1*w1 + 2*i2*o2*p2*s1*w1 + 2*k2*l2*n1*s2*w1 + 
  2*k2*l1*n2*s2*w1 + 2*k1*l2*n2*s2*w1 - 2*j2*l2*o1*s2*w1 - 
  2*j2*l1*o2*s2*w1 - 2*j1*l2*o2*s2*w1 - 2*j2*k2*p1*s2*w1 + 
  2*i2*o2*p1*s2*w1 - 2*j2*k1*p2*s2*w1 - 2*j1*k2*p2*s2*w1 + 
  2*i2*o1*p2*s2*w1 + 2*i1*o2*p2*s2*w1 + 2*j2*j2*s1*s2*w1 - 
  2*i2*n2*s1*s2*w1 + 2*j1*j2*s2*s2*w1 - i2*n1*s2*s2*w1 - 
  i1*n2*s2*s2*w1 - k2*k2*n2*u1*w1 + 2*j2*k2*o2*u1*w1 - 
  i2*o2*o2*u1*w1 - j2*j2*r2*u1*w1 + i2*n2*r2*u1*w1 - 
  k2*k2*n1*u2*w1 - 2*k1*k2*n2*u2*w1 + 2*j2*k2*o1*u2*w1 + 
  2*j2*k1*o2*u2*w1 + 2*j1*k2*o2*u2*w1 - 2*i2*o1*o2*u2*w1 - 
  i1*o2*o2*u2*w1 - j2*j2*r1*u2*w1 + i2*n2*r1*u2*w1 - 
  2*j1*j2*r2*u2*w1 + i2*n1*r2*u2*w1 + i1*n2*r2*u2*w1 + 
  l2*l2*o1*o1*w2 + 4*l1*l2*o1*o2*w2 + l1*l1*o2*o2*w2 - 
  2*k2*l2*o1*p1*w2 - 2*k2*l1*o2*p1*w2 - 2*k1*l2*o2*p1*w2 + 
  k2*k2*p1*p1*w2 - 2*k2*l1*o1*p2*w2 - 2*k1*l2*o1*p2*w2 - 
  2*k1*l1*o2*p2*w2 + 4*k1*k2*p1*p2*w2 + k1*k1*p2*p2*w2 - 
  l2*l2*n1*r1*w2 - 2*l1*l2*n2*r1*w2 + 2*j2*l2*p1*r1*w2 + 
  2*j2*l1*p2*r1*w2 + 2*j1*l2*p2*r1*w2 - 2*i2*p1*p2*r1*w2 - 
  i1*p2*p2*r1*w2 - 2*l1*l2*n1*r2*w2 - l1*l1*n2*r2*w2 + 
  2*j2*l1*p1*r2*w2 + 2*j1*l2*p1*r2*w2 - i2*p1*p1*r2*w2 + 
  2*j1*l1*p2*r2*w2 - 2*i1*p1*p2*r2*w2 + 2*k2*l2*n1*s1*w2 + 
  2*k2*l1*n2*s1*w2 + 2*k1*l2*n2*s1*w2 - 2*j2*l2*o1*s1*w2 - 
  2*j2*l1*o2*s1*w2 - 2*j1*l2*o2*s1*w2 - 2*j2*k2*p1*s1*w2 + 
  2*i2*o2*p1*s1*w2 - 2*j2*k1*p2*s1*w2 - 2*j1*k2*p2*s1*w2 + 
  2*i2*o1*p2*s1*w2 + 2*i1*o2*p2*s1*w2 + j2*j2*s1*s1*w2 - 
  i2*n2*s1*s1*w2 + 2*k2*l1*n1*s2*w2 + 2*k1*l2*n1*s2*w2 + 
  2*k1*l1*n2*s2*w2 - 2*j2*l1*o1*s2*w2 - 2*j1*l2*o1*s2*w2 - 
  2*j1*l1*o2*s2*w2 - 2*j2*k1*p1*s2*w2 - 2*j1*k2*p1*s2*w2 + 
  2*i2*o1*p1*s2*w2 + 2*i1*o2*p1*s2*w2 - 2*j1*k1*p2*s2*w2 + 
  2*i1*o1*p2*s2*w2 + 4*j1*j2*s1*s2*w2 - 2*i2*n1*s1*s2*w2 - 
  2*i1*n2*s1*s2*w2 + j1*j1*s2*s2*w2 - i1*n1*s2*s2*w2 - 
  k2*k2*n1*u1*w2 - 2*k1*k2*n2*u1*w2 + 2*j2*k2*o1*u1*w2 + 
  2*j2*k1*o2*u1*w2 + 2*j1*k2*o2*u1*w2 - 2*i2*o1*o2*u1*w2 - 
  i1*o2*o2*u1*w2 - j2*j2*r1*u1*w2 + i2*n2*r1*u1*w2 - 
  2*j1*j2*r2*u1*w2 + i2*n1*r2*u1*w2 + i1*n2*r2*u1*w2 - 
  2*k1*k2*n1*u2*w2 - k1*k1*n2*u2*w2 + 2*j2*k1*o1*u2*w2 + 
  2*j1*k2*o1*u2*w2 - i2*o1*o1*u2*w2 + 2*j1*k1*o2*u2*w2 - 
  2*i1*o1*o2*u2*w2 - 2*j1*j2*r1*u2*w2 + i2*n1*r1*u2*w2 + 
  i1*n2*r1*u2*w2 - j1*j1*r2*u2*w2 + i1*n1*r2*u2*w2;
    
    k[6]  =   m2*m2*p0*p0*r0 + 4*m0*m2*p0*p2*r0 + m0*m0*p2*p2*r0 - 
  2*l2*m2*p0*q0*r0 - 2*l2*m0*p2*q0*r0 - 2*l0*m2*p2*q0*r0 + 
  l2*l2*q0*q0*r0 - 2*l2*m0*p0*q2*r0 - 2*l0*m2*p0*q2*r0 - 
  2*l0*m0*p2*q2*r0 + 4*l0*l2*q0*q2*r0 + l0*l0*q2*q2*r0 + 
  2*m0*m2*p0*p0*r2 + 2*m0*m0*p0*p2*r2 - 2*l2*m0*p0*q0*r2 - 
  2*l0*m2*p0*q0*r2 - 2*l0*m0*p2*q0*r2 + 2*l0*l2*q0*q0*r2 - 
  2*l0*m0*p0*q2*r2 + 2*l0*l0*q0*q2*r2 - 2*m2*m2*o0*p0*s0 - 
  4*m0*m2*o2*p0*s0 - 4*m0*m2*o0*p2*s0 - 2*m0*m0*o2*p2*s0 + 
  2*l2*m2*o0*q0*s0 + 2*l2*m0*o2*q0*s0 + 2*l0*m2*o2*q0*s0 + 
  2*k2*m2*p0*q0*s0 + 2*k2*m0*p2*q0*s0 + 2*k0*m2*p2*q0*s0 - 
  2*k2*l2*q0*q0*s0 + 2*l2*m0*o0*q2*s0 + 2*l0*m2*o0*q2*s0 + 
  2*l0*m0*o2*q2*s0 + 2*k2*m0*p0*q2*s0 + 2*k0*m2*p0*q2*s0 + 
  2*k0*m0*p2*q2*s0 - 4*k2*l0*q0*q2*s0 - 4*k0*l2*q0*q2*s0 - 
  2*k0*l0*q2*q2*s0 + m2*m2*n0*s0*s0 + 2*m0*m2*n2*s0*s0 - 
  2*j2*m2*q0*s0*s0 - 2*j2*m0*q2*s0*s0 - 2*j0*m2*q2*s0*s0 + 
  2*i2*q0*q2*s0*s0 + i0*q2*q2*s0*s0 - 4*m0*m2*o0*p0*s2 - 
  2*m0*m0*o2*p0*s2 - 2*m0*m0*o0*p2*s2 + 2*l2*m0*o0*q0*s2 + 
  2*l0*m2*o0*q0*s2 + 2*l0*m0*o2*q0*s2 + 2*k2*m0*p0*q0*s2 + 
  2*k0*m2*p0*q0*s2 + 2*k0*m0*p2*q0*s2 - 2*k2*l0*q0*q0*s2 - 
  2*k0*l2*q0*q0*s2 + 2*l0*m0*o0*q2*s2 + 2*k0*m0*p0*q2*s2 - 
  4*k0*l0*q0*q2*s2 + 4*m0*m2*n0*s0*s2 + 2*m0*m0*n2*s0*s2 - 
  4*j2*m0*q0*s0*s2 - 4*j0*m2*q0*s0*s2 + 2*i2*q0*q0*s0*s2 - 
  4*j0*m0*q2*s0*s2 + 4*i0*q0*q2*s0*s2 + m0*m0*n0*s2*s2 - 
  2*j0*m0*q0*s2*s2 + i0*q0*q0*s2*s2 + 2*l2*m2*o0*p0*t0 + 
  2*l2*m0*o2*p0*t0 + 2*l0*m2*o2*p0*t0 - 2*k2*m2*p0*p0*t0 + 
  2*l2*m0*o0*p2*t0 + 2*l0*m2*o0*p2*t0 + 2*l0*m0*o2*p2*t0 - 
  4*k2*m0*p0*p2*t0 - 4*k0*m2*p0*p2*t0 - 2*k0*m0*p2*p2*t0 - 
  2*l2*l2*o0*q0*t0 - 4*l0*l2*o2*q0*t0 + 2*k2*l2*p0*q0*t0 + 
  2*k2*l0*p2*q0*t0 + 2*k0*l2*p2*q0*t0 - 4*l0*l2*o0*q2*t0 - 
  2*l0*l0*o2*q2*t0 + 2*k2*l0*p0*q2*t0 + 2*k0*l2*p0*q2*t0 + 
  2*k0*l0*p2*q2*t0 - 2*l2*m2*n0*s0*t0 - 2*l2*m0*n2*s0*t0 - 
  2*l0*m2*n2*s0*t0 + 2*j2*m2*p0*s0*t0 + 2*j2*m0*p2*s0*t0 + 
  2*j0*m2*p2*s0*t0 + 2*j2*l2*q0*s0*t0 - 2*i2*p2*q0*s0*t0 + 
  2*j2*l0*q2*s0*t0 + 2*j0*l2*q2*s0*t0 - 2*i2*p0*q2*s0*t0 - 
  2*i0*p2*q2*s0*t0 - 2*l2*m0*n0*s2*t0 - 2*l0*m2*n0*s2*t0 - 
  2*l0*m0*n2*s2*t0 + 2*j2*m0*p0*s2*t0 + 2*j0*m2*p0*s2*t0 + 
  2*j0*m0*p2*s2*t0 + 2*j2*l0*q0*s2*t0 + 2*j0*l2*q0*s2*t0 - 
  2*i2*p0*q0*s2*t0 - 2*i0*p2*q0*s2*t0 + 2*j0*l0*q2*s2*t0 - 
  2*i0*p0*q2*s2*t0 + l2*l2*n0*t0*t0 + 2*l0*l2*n2*t0*t0 - 
  2*j2*l2*p0*t0*t0 - 2*j2*l0*p2*t0*t0 - 2*j0*l2*p2*t0*t0 + 
  2*i2*p0*p2*t0*t0 + i0*p2*p2*t0*t0 + 2*l2*m0*o0*p0*t2 + 
  2*l0*m2*o0*p0*t2 + 2*l0*m0*o2*p0*t2 - 2*k2*m0*p0*p0*t2 - 
  2*k0*m2*p0*p0*t2 + 2*l0*m0*o0*p2*t2 - 4*k0*m0*p0*p2*t2 - 
  4*l0*l2*o0*q0*t2 - 2*l0*l0*o2*q0*t2 + 2*k2*l0*p0*q0*t2 + 
  2*k0*l2*p0*q0*t2 + 2*k0*l0*p2*q0*t2 - 2*l0*l0*o0*q2*t2 + 
  2*k0*l0*p0*q2*t2 - 2*l2*m0*n0*s0*t2 - 2*l0*m2*n0*s0*t2 - 
  2*l0*m0*n2*s0*t2 + 2*j2*m0*p0*s0*t2 + 2*j0*m2*p0*s0*t2 + 
  2*j0*m0*p2*s0*t2 + 2*j2*l0*q0*s0*t2 + 2*j0*l2*q0*s0*t2 - 
  2*i2*p0*q0*s0*t2 - 2*i0*p2*q0*s0*t2 + 2*j0*l0*q2*s0*t2 - 
  2*i0*p0*q2*s0*t2 - 2*l0*m0*n0*s2*t2 + 2*j0*m0*p0*s2*t2 + 
  2*j0*l0*q0*s2*t2 - 2*i0*p0*q0*s2*t2 + 4*l0*l2*n0*t0*t2 + 
  2*l0*l0*n2*t0*t2 - 4*j2*l0*p0*t0*t2 - 4*j0*l2*p0*t0*t2 + 
  2*i2*p0*p0*t0*t2 - 4*j0*l0*p2*t0*t2 + 4*i0*p0*p2*t0*t2 + 
  l0*l0*n0*t2*t2 - 2*j0*l0*p0*t2*t2 + i0*p0*p0*t2*t2 + 
  m2*m2*o0*o0*u0 + 4*m0*m2*o0*o2*u0 + m0*m0*o2*o2*u0 - 
  2*k2*m2*o0*q0*u0 - 2*k2*m0*o2*q0*u0 - 2*k0*m2*o2*q0*u0 + 
  k2*k2*q0*q0*u0 - 2*k2*m0*o0*q2*u0 - 2*k0*m2*o0*q2*u0 - 
  2*k0*m0*o2*q2*u0 + 4*k0*k2*q0*q2*u0 + k0*k0*q2*q2*u0 - 
  m2*m2*n0*r0*u0 - 2*m0*m2*n2*r0*u0 + 2*j2*m2*q0*r0*u0 + 
  2*j2*m0*q2*r0*u0 + 2*j0*m2*q2*r0*u0 - 2*i2*q0*q2*r0*u0 - 
  i0*q2*q2*r0*u0 - 2*m0*m2*n0*r2*u0 - m0*m0*n2*r2*u0 + 
  2*j2*m0*q0*r2*u0 + 2*j0*m2*q0*r2*u0 - i2*q0*q0*r2*u0 + 
  2*j0*m0*q2*r2*u0 - 2*i0*q0*q2*r2*u0 + 2*k2*m2*n0*t0*u0 + 
  2*k2*m0*n2*t0*u0 + 2*k0*m2*n2*t0*u0 - 2*j2*m2*o0*t0*u0 - 
  2*j2*m0*o2*t0*u0 - 2*j0*m2*o2*t0*u0 - 2*j2*k2*q0*t0*u0 + 
  2*i2*o2*q0*t0*u0 - 2*j2*k0*q2*t0*u0 - 2*j0*k2*q2*t0*u0 + 
  2*i2*o0*q2*t0*u0 + 2*i0*o2*q2*t0*u0 + j2*j2*t0*t0*u0 - 
  i2*n2*t0*t0*u0 + 2*k2*m0*n0*t2*u0 + 2*k0*m2*n0*t2*u0 + 
  2*k0*m0*n2*t2*u0 - 2*j2*m0*o0*t2*u0 - 2*j0*m2*o0*t2*u0 - 
  2*j0*m0*o2*t2*u0 - 2*j2*k0*q0*t2*u0 - 2*j0*k2*q0*t2*u0 + 
  2*i2*o0*q0*t2*u0 + 2*i0*o2*q0*t2*u0 - 2*j0*k0*q2*t2*u0 + 
  2*i0*o0*q2*t2*u0 + 4*j0*j2*t0*t2*u0 - 2*i2*n0*t0*t2*u0 - 
  2*i0*n2*t0*t2*u0 + j0*j0*t2*t2*u0 - i0*n0*t2*t2*u0 + 
  2*m0*m2*o0*o0*u2 + 2*m0*m0*o0*o2*u2 - 2*k2*m0*o0*q0*u2 - 
  2*k0*m2*o0*q0*u2 - 2*k0*m0*o2*q0*u2 + 2*k0*k2*q0*q0*u2 - 
  2*k0*m0*o0*q2*u2 + 2*k0*k0*q0*q2*u2 - 2*m0*m2*n0*r0*u2 - 
  m0*m0*n2*r0*u2 + 2*j2*m0*q0*r0*u2 + 2*j0*m2*q0*r0*u2 - 
  i2*q0*q0*r0*u2 + 2*j0*m0*q2*r0*u2 - 2*i0*q0*q2*r0*u2 - 
  m0*m0*n0*r2*u2 + 2*j0*m0*q0*r2*u2 - i0*q0*q0*r2*u2 + 
  2*k2*m0*n0*t0*u2 + 2*k0*m2*n0*t0*u2 + 2*k0*m0*n2*t0*u2 - 
  2*j2*m0*o0*t0*u2 - 2*j0*m2*o0*t0*u2 - 2*j0*m0*o2*t0*u2 - 
  2*j2*k0*q0*t0*u2 - 2*j0*k2*q0*t0*u2 + 2*i2*o0*q0*t0*u2 + 
  2*i0*o2*q0*t0*u2 - 2*j0*k0*q2*t0*u2 + 2*i0*o0*q2*t0*u2 + 
  2*j0*j2*t0*t0*u2 - i2*n0*t0*t0*u2 - i0*n2*t0*t0*u2 + 
  2*k0*m0*n0*t2*u2 - 2*j0*m0*o0*t2*u2 - 2*j0*k0*q0*t2*u2 + 
  2*i0*o0*q0*t2*u2 + 2*j0*j0*t0*t2*u2 - 2*i0*n0*t0*t2*u2 - 
  2*l2*m2*o0*o0*v0 - 4*l2*m0*o0*o2*v0 - 4*l0*m2*o0*o2*v0 - 
  2*l0*m0*o2*o2*v0 + 2*k2*m2*o0*p0*v0 + 2*k2*m0*o2*p0*v0 + 
  2*k0*m2*o2*p0*v0 + 2*k2*m0*o0*p2*v0 + 2*k0*m2*o0*p2*v0 + 
  2*k0*m0*o2*p2*v0 + 2*k2*l2*o0*q0*v0 + 2*k2*l0*o2*q0*v0 + 
  2*k0*l2*o2*q0*v0 - 2*k2*k2*p0*q0*v0 - 4*k0*k2*p2*q0*v0 + 
  2*k2*l0*o0*q2*v0 + 2*k0*l2*o0*q2*v0 + 2*k0*l0*o2*q2*v0 - 
  4*k0*k2*p0*q2*v0 - 2*k0*k0*p2*q2*v0 + 2*l2*m2*n0*r0*v0 + 
  2*l2*m0*n2*r0*v0 + 2*l0*m2*n2*r0*v0 - 2*j2*m2*p0*r0*v0 - 
  2*j2*m0*p2*r0*v0 - 2*j0*m2*p2*r0*v0 - 2*j2*l2*q0*r0*v0 + 
  2*i2*p2*q0*r0*v0 - 2*j2*l0*q2*r0*v0 - 2*j0*l2*q2*r0*v0 + 
  2*i2*p0*q2*r0*v0 + 2*i0*p2*q2*r0*v0 + 2*l2*m0*n0*r2*v0 + 
  2*l0*m2*n0*r2*v0 + 2*l0*m0*n2*r2*v0 - 2*j2*m0*p0*r2*v0 - 
  2*j0*m2*p0*r2*v0 - 2*j0*m0*p2*r2*v0 - 2*j2*l0*q0*r2*v0 - 
  2*j0*l2*q0*r2*v0 + 2*i2*p0*q0*r2*v0 + 2*i0*p2*q0*r2*v0 - 
  2*j0*l0*q2*r2*v0 + 2*i0*p0*q2*r2*v0 - 2*k2*m2*n0*s0*v0 - 
  2*k2*m0*n2*s0*v0 - 2*k0*m2*n2*s0*v0 + 2*j2*m2*o0*s0*v0 + 
  2*j2*m0*o2*s0*v0 + 2*j0*m2*o2*s0*v0 + 2*j2*k2*q0*s0*v0 - 
  2*i2*o2*q0*s0*v0 + 2*j2*k0*q2*s0*v0 + 2*j0*k2*q2*s0*v0 - 
  2*i2*o0*q2*s0*v0 - 2*i0*o2*q2*s0*v0 - 2*k2*m0*n0*s2*v0 - 
  2*k0*m2*n0*s2*v0 - 2*k0*m0*n2*s2*v0 + 2*j2*m0*o0*s2*v0 + 
  2*j0*m2*o0*s2*v0 + 2*j0*m0*o2*s2*v0 + 2*j2*k0*q0*s2*v0 + 
  2*j0*k2*q0*s2*v0 - 2*i2*o0*q0*s2*v0 - 2*i0*o2*q0*s2*v0 + 
  2*j0*k0*q2*s2*v0 - 2*i0*o0*q2*s2*v0 - 2*k2*l2*n0*t0*v0 - 
  2*k2*l0*n2*t0*v0 - 2*k0*l2*n2*t0*v0 + 2*j2*l2*o0*t0*v0 + 
  2*j2*l0*o2*t0*v0 + 2*j0*l2*o2*t0*v0 + 2*j2*k2*p0*t0*v0 - 
  2*i2*o2*p0*t0*v0 + 2*j2*k0*p2*t0*v0 + 2*j0*k2*p2*t0*v0 - 
  2*i2*o0*p2*t0*v0 - 2*i0*o2*p2*t0*v0 - 2*j2*j2*s0*t0*v0 + 
  2*i2*n2*s0*t0*v0 - 4*j0*j2*s2*t0*v0 + 2*i2*n0*s2*t0*v0 + 
  2*i0*n2*s2*t0*v0 - 2*k2*l0*n0*t2*v0 - 2*k0*l2*n0*t2*v0 - 
  2*k0*l0*n2*t2*v0 + 2*j2*l0*o0*t2*v0 + 2*j0*l2*o0*t2*v0 + 
  2*j0*l0*o2*t2*v0 + 2*j2*k0*p0*t2*v0 + 2*j0*k2*p0*t2*v0 - 
  2*i2*o0*p0*t2*v0 - 2*i0*o2*p0*t2*v0 + 2*j0*k0*p2*t2*v0 - 
  2*i0*o0*p2*t2*v0 - 4*j0*j2*s0*t2*v0 + 2*i2*n0*s0*t2*v0 + 
  2*i0*n2*s0*t2*v0 - 2*j0*j0*s2*t2*v0 + 2*i0*n0*s2*t2*v0 + 
  k2*k2*n0*v0*v0 + 2*k0*k2*n2*v0*v0 - 2*j2*k2*o0*v0*v0 - 
  2*j2*k0*o2*v0*v0 - 2*j0*k2*o2*v0*v0 + 2*i2*o0*o2*v0*v0 + 
  i0*o2*o2*v0*v0 + j2*j2*r0*v0*v0 - i2*n2*r0*v0*v0 + 
  2*j0*j2*r2*v0*v0 - i2*n0*r2*v0*v0 - i0*n2*r2*v0*v0 - 
  2*l2*m0*o0*o0*v2 - 2*l0*m2*o0*o0*v2 - 4*l0*m0*o0*o2*v2 + 
  2*k2*m0*o0*p0*v2 + 2*k0*m2*o0*p0*v2 + 2*k0*m0*o2*p0*v2 + 
  2*k0*m0*o0*p2*v2 + 2*k2*l0*o0*q0*v2 + 2*k0*l2*o0*q0*v2 + 
  2*k0*l0*o2*q0*v2 - 4*k0*k2*p0*q0*v2 - 2*k0*k0*p2*q0*v2 + 
  2*k0*l0*o0*q2*v2 - 2*k0*k0*p0*q2*v2 + 2*l2*m0*n0*r0*v2 + 
  2*l0*m2*n0*r0*v2 + 2*l0*m0*n2*r0*v2 - 2*j2*m0*p0*r0*v2 - 
  2*j0*m2*p0*r0*v2 - 2*j0*m0*p2*r0*v2 - 2*j2*l0*q0*r0*v2 - 
  2*j0*l2*q0*r0*v2 + 2*i2*p0*q0*r0*v2 + 2*i0*p2*q0*r0*v2 - 
  2*j0*l0*q2*r0*v2 + 2*i0*p0*q2*r0*v2 + 2*l0*m0*n0*r2*v2 - 
  2*j0*m0*p0*r2*v2 - 2*j0*l0*q0*r2*v2 + 2*i0*p0*q0*r2*v2 - 
  2*k2*m0*n0*s0*v2 - 2*k0*m2*n0*s0*v2 - 2*k0*m0*n2*s0*v2 + 
  2*j2*m0*o0*s0*v2 + 2*j0*m2*o0*s0*v2 + 2*j0*m0*o2*s0*v2 + 
  2*j2*k0*q0*s0*v2 + 2*j0*k2*q0*s0*v2 - 2*i2*o0*q0*s0*v2 - 
  2*i0*o2*q0*s0*v2 + 2*j0*k0*q2*s0*v2 - 2*i0*o0*q2*s0*v2 - 
  2*k0*m0*n0*s2*v2 + 2*j0*m0*o0*s2*v2 + 2*j0*k0*q0*s2*v2 - 
  2*i0*o0*q0*s2*v2 - 2*k2*l0*n0*t0*v2 - 2*k0*l2*n0*t0*v2 - 
  2*k0*l0*n2*t0*v2 + 2*j2*l0*o0*t0*v2 + 2*j0*l2*o0*t0*v2 + 
  2*j0*l0*o2*t0*v2 + 2*j2*k0*p0*t0*v2 + 2*j0*k2*p0*t0*v2 - 
  2*i2*o0*p0*t0*v2 - 2*i0*o2*p0*t0*v2 + 2*j0*k0*p2*t0*v2 - 
  2*i0*o0*p2*t0*v2 - 4*j0*j2*s0*t0*v2 + 2*i2*n0*s0*t0*v2 + 
  2*i0*n2*s0*t0*v2 - 2*j0*j0*s2*t0*v2 + 2*i0*n0*s2*t0*v2 - 
  2*k0*l0*n0*t2*v2 + 2*j0*l0*o0*t2*v2 + 2*j0*k0*p0*t2*v2 - 
  2*i0*o0*p0*t2*v2 - 2*j0*j0*s0*t2*v2 + 2*i0*n0*s0*t2*v2 + 
  4*k0*k2*n0*v0*v2 + 2*k0*k0*n2*v0*v2 - 4*j2*k0*o0*v0*v2 - 
  4*j0*k2*o0*v0*v2 + 2*i2*o0*o0*v0*v2 - 4*j0*k0*o2*v0*v2 + 
  4*i0*o0*o2*v0*v2 + 4*j0*j2*r0*v0*v2 - 2*i2*n0*r0*v0*v2 - 
  2*i0*n2*r0*v0*v2 + 2*j0*j0*r2*v0*v2 - 2*i0*n0*r2*v0*v2 + 
  k0*k0*n0*v2*v2 - 2*j0*k0*o0*v2*v2 + i0*o0*o0*v2*v2 + 
  j0*j0*r0*v2*v2 - i0*n0*r0*v2*v2 + l2*l2*o0*o0*w0 + 
  4*l0*l2*o0*o2*w0 + l0*l0*o2*o2*w0 - 2*k2*l2*o0*p0*w0 - 
  2*k2*l0*o2*p0*w0 - 2*k0*l2*o2*p0*w0 + k2*k2*p0*p0*w0 - 
  2*k2*l0*o0*p2*w0 - 2*k0*l2*o0*p2*w0 - 2*k0*l0*o2*p2*w0 + 
  4*k0*k2*p0*p2*w0 + k0*k0*p2*p2*w0 - l2*l2*n0*r0*w0 - 
  2*l0*l2*n2*r0*w0 + 2*j2*l2*p0*r0*w0 + 2*j2*l0*p2*r0*w0 + 
  2*j0*l2*p2*r0*w0 - 2*i2*p0*p2*r0*w0 - i0*p2*p2*r0*w0 - 
  2*l0*l2*n0*r2*w0 - l0*l0*n2*r2*w0 + 2*j2*l0*p0*r2*w0 + 
  2*j0*l2*p0*r2*w0 - i2*p0*p0*r2*w0 + 2*j0*l0*p2*r2*w0 - 
  2*i0*p0*p2*r2*w0 + 2*k2*l2*n0*s0*w0 + 2*k2*l0*n2*s0*w0 + 
  2*k0*l2*n2*s0*w0 - 2*j2*l2*o0*s0*w0 - 2*j2*l0*o2*s0*w0 - 
  2*j0*l2*o2*s0*w0 - 2*j2*k2*p0*s0*w0 + 2*i2*o2*p0*s0*w0 - 
  2*j2*k0*p2*s0*w0 - 2*j0*k2*p2*s0*w0 + 2*i2*o0*p2*s0*w0 + 
  2*i0*o2*p2*s0*w0 + j2*j2*s0*s0*w0 - i2*n2*s0*s0*w0 + 
  2*k2*l0*n0*s2*w0 + 2*k0*l2*n0*s2*w0 + 2*k0*l0*n2*s2*w0 - 
  2*j2*l0*o0*s2*w0 - 2*j0*l2*o0*s2*w0 - 2*j0*l0*o2*s2*w0 - 
  2*j2*k0*p0*s2*w0 - 2*j0*k2*p0*s2*w0 + 2*i2*o0*p0*s2*w0 + 
  2*i0*o2*p0*s2*w0 - 2*j0*k0*p2*s2*w0 + 2*i0*o0*p2*s2*w0 + 
  4*j0*j2*s0*s2*w0 - 2*i2*n0*s0*s2*w0 - 2*i0*n2*s0*s2*w0 + 
  j0*j0*s2*s2*w0 - i0*n0*s2*s2*w0 - k2*k2*n0*u0*w0 - 
  2*k0*k2*n2*u0*w0 + 2*j2*k2*o0*u0*w0 + 2*j2*k0*o2*u0*w0 + 
  2*j0*k2*o2*u0*w0 - 2*i2*o0*o2*u0*w0 - i0*o2*o2*u0*w0 - 
  j2*j2*r0*u0*w0 + i2*n2*r0*u0*w0 - 2*j0*j2*r2*u0*w0 + 
  i2*n0*r2*u0*w0 + i0*n2*r2*u0*w0 - 2*k0*k2*n0*u2*w0 - 
  k0*k0*n2*u2*w0 + 2*j2*k0*o0*u2*w0 + 2*j0*k2*o0*u2*w0 - 
  i2*o0*o0*u2*w0 + 2*j0*k0*o2*u2*w0 - 2*i0*o0*o2*u2*w0 - 
  2*j0*j2*r0*u2*w0 + i2*n0*r0*u2*w0 + i0*n2*r0*u2*w0 - 
  j0*j0*r2*u2*w0 + i0*n0*r2*u2*w0 + 2*l0*l2*o0*o0*w2 + 
  2*l0*l0*o0*o2*w2 - 2*k2*l0*o0*p0*w2 - 2*k0*l2*o0*p0*w2 - 
  2*k0*l0*o2*p0*w2 + 2*k0*k2*p0*p0*w2 - 2*k0*l0*o0*p2*w2 + 
  2*k0*k0*p0*p2*w2 - 2*l0*l2*n0*r0*w2 - l0*l0*n2*r0*w2 + 
  2*j2*l0*p0*r0*w2 + 2*j0*l2*p0*r0*w2 - i2*p0*p0*r0*w2 + 
  2*j0*l0*p2*r0*w2 - 2*i0*p0*p2*r0*w2 - l0*l0*n0*r2*w2 + 
  2*j0*l0*p0*r2*w2 - i0*p0*p0*r2*w2 + 2*k2*l0*n0*s0*w2 + 
  2*k0*l2*n0*s0*w2 + 2*k0*l0*n2*s0*w2 - 2*j2*l0*o0*s0*w2 - 
  2*j0*l2*o0*s0*w2 - 2*j0*l0*o2*s0*w2 - 2*j2*k0*p0*s0*w2 - 
  2*j0*k2*p0*s0*w2 + 2*i2*o0*p0*s0*w2 + 2*i0*o2*p0*s0*w2 - 
  2*j0*k0*p2*s0*w2 + 2*i0*o0*p2*s0*w2 + 2*j0*j2*s0*s0*w2 - 
  i2*n0*s0*s0*w2 - i0*n2*s0*s0*w2 + 2*k0*l0*n0*s2*w2 - 
  2*j0*l0*o0*s2*w2 - 2*j0*k0*p0*s2*w2 + 2*i0*o0*p0*s2*w2 + 
  2*j0*j0*s0*s2*w2 - 2*i0*n0*s0*s2*w2 - 2*k0*k2*n0*u0*w2 - 
  k0*k0*n2*u0*w2 + 2*j2*k0*o0*u0*w2 + 2*j0*k2*o0*u0*w2 - 
  i2*o0*o0*u0*w2 + 2*j0*k0*o2*u0*w2 - 2*i0*o0*o2*u0*w2 - 
  2*j0*j2*r0*u0*w2 + i2*n0*r0*u0*w2 + i0*n2*r0*u0*w2 - 
  j0*j0*r2*u0*w2 + i0*n0*r2*u0*w2 - k0*k0*n0*u2*w2 + 
  2*j0*k0*o0*u2*w2 - i0*o0*o0*u2*w2 - j0*j0*r0*u2*w2 + 
  i0*n0*r0*u2*w2;
    
    k[7]  = 2*m2*m2*p0*p1*r0 + 
  4*m1*m2*p0*p2*r0 + 4*m0*m2*p1*p2*r0 + 2*m0*m1*p2*p2*r0 - 
  2*l2*m2*p1*q0*r0 - 2*l2*m1*p2*q0*r0 - 2*l1*m2*p2*q0*r0 - 
  2*l2*m2*p0*q1*r0 - 2*l2*m0*p2*q1*r0 - 2*l0*m2*p2*q1*r0 + 
  2*l2*l2*q0*q1*r0 - 2*l2*m1*p0*q2*r0 - 2*l1*m2*p0*q2*r0 - 
  2*l2*m0*p1*q2*r0 - 2*l0*m2*p1*q2*r0 - 2*l1*m0*p2*q2*r0 - 
  2*l0*m1*p2*q2*r0 + 4*l1*l2*q0*q2*r0 + 4*l0*l2*q1*q2*r0 + 
  2*l0*l1*q2*q2*r0 + m2*m2*p0*p0*r1 + 4*m0*m2*p0*p2*r1 + 
  m0*m0*p2*p2*r1 - 2*l2*m2*p0*q0*r1 - 2*l2*m0*p2*q0*r1 - 
  2*l0*m2*p2*q0*r1 + l2*l2*q0*q0*r1 - 2*l2*m0*p0*q2*r1 - 
  2*l0*m2*p0*q2*r1 - 2*l0*m0*p2*q2*r1 + 4*l0*l2*q0*q2*r1 + 
  l0*l0*q2*q2*r1 + 2*m1*m2*p0*p0*r2 + 4*m0*m2*p0*p1*r2 + 
  4*m0*m1*p0*p2*r2 + 2*m0*m0*p1*p2*r2 - 2*l2*m1*p0*q0*r2 - 
  2*l1*m2*p0*q0*r2 - 2*l2*m0*p1*q0*r2 - 2*l0*m2*p1*q0*r2 - 
  2*l1*m0*p2*q0*r2 - 2*l0*m1*p2*q0*r2 + 2*l1*l2*q0*q0*r2 - 
  2*l2*m0*p0*q1*r2 - 2*l0*m2*p0*q1*r2 - 2*l0*m0*p2*q1*r2 + 
  4*l0*l2*q0*q1*r2 - 2*l1*m0*p0*q2*r2 - 2*l0*m1*p0*q2*r2 - 
  2*l0*m0*p1*q2*r2 + 4*l0*l1*q0*q2*r2 + 2*l0*l0*q1*q2*r2 - 
  2*m2*m2*o1*p0*s0 - 4*m1*m2*o2*p0*s0 - 2*m2*m2*o0*p1*s0 - 
  4*m0*m2*o2*p1*s0 - 4*m1*m2*o0*p2*s0 - 4*m0*m2*o1*p2*s0 - 
  4*m0*m1*o2*p2*s0 + 2*l2*m2*o1*q0*s0 + 2*l2*m1*o2*q0*s0 + 
  2*l1*m2*o2*q0*s0 + 2*k2*m2*p1*q0*s0 + 2*k2*m1*p2*q0*s0 + 
  2*k1*m2*p2*q0*s0 + 2*l2*m2*o0*q1*s0 + 2*l2*m0*o2*q1*s0 + 
  2*l0*m2*o2*q1*s0 + 2*k2*m2*p0*q1*s0 + 2*k2*m0*p2*q1*s0 + 
  2*k0*m2*p2*q1*s0 - 4*k2*l2*q0*q1*s0 + 2*l2*m1*o0*q2*s0 + 
  2*l1*m2*o0*q2*s0 + 2*l2*m0*o1*q2*s0 + 2*l0*m2*o1*q2*s0 + 
  2*l1*m0*o2*q2*s0 + 2*l0*m1*o2*q2*s0 + 2*k2*m1*p0*q2*s0 + 
  2*k1*m2*p0*q2*s0 + 2*k2*m0*p1*q2*s0 + 2*k0*m2*p1*q2*s0 + 
  2*k1*m0*p2*q2*s0 + 2*k0*m1*p2*q2*s0 - 4*k2*l1*q0*q2*s0 - 
  4*k1*l2*q0*q2*s0 - 4*k2*l0*q1*q2*s0 - 4*k0*l2*q1*q2*s0 - 
  2*k1*l0*q2*q2*s0 - 2*k0*l1*q2*q2*s0 + m2*m2*n1*s0*s0 + 
  2*m1*m2*n2*s0*s0 - 2*j2*m2*q1*s0*s0 - 2*j2*m1*q2*s0*s0 - 
  2*j1*m2*q2*s0*s0 + 2*i2*q1*q2*s0*s0 + i1*q2*q2*s0*s0 - 
  2*m2*m2*o0*p0*s1 - 4*m0*m2*o2*p0*s1 - 4*m0*m2*o0*p2*s1 - 
  2*m0*m0*o2*p2*s1 + 2*l2*m2*o0*q0*s1 + 2*l2*m0*o2*q0*s1 + 
  2*l0*m2*o2*q0*s1 + 2*k2*m2*p0*q0*s1 + 2*k2*m0*p2*q0*s1 + 
  2*k0*m2*p2*q0*s1 - 2*k2*l2*q0*q0*s1 + 2*l2*m0*o0*q2*s1 + 
  2*l0*m2*o0*q2*s1 + 2*l0*m0*o2*q2*s1 + 2*k2*m0*p0*q2*s1 + 
  2*k0*m2*p0*q2*s1 + 2*k0*m0*p2*q2*s1 - 4*k2*l0*q0*q2*s1 - 
  4*k0*l2*q0*q2*s1 - 2*k0*l0*q2*q2*s1 + 2*m2*m2*n0*s0*s1 + 
  4*m0*m2*n2*s0*s1 - 4*j2*m2*q0*s0*s1 - 4*j2*m0*q2*s0*s1 - 
  4*j0*m2*q2*s0*s1 + 4*i2*q0*q2*s0*s1 + 2*i0*q2*q2*s0*s1 - 
  4*m1*m2*o0*p0*s2 - 4*m0*m2*o1*p0*s2 - 4*m0*m1*o2*p0*s2 - 
  4*m0*m2*o0*p1*s2 - 2*m0*m0*o2*p1*s2 - 4*m0*m1*o0*p2*s2 - 
  2*m0*m0*o1*p2*s2 + 2*l2*m1*o0*q0*s2 + 2*l1*m2*o0*q0*s2 + 
  2*l2*m0*o1*q0*s2 + 2*l0*m2*o1*q0*s2 + 2*l1*m0*o2*q0*s2 + 
  2*l0*m1*o2*q0*s2 + 2*k2*m1*p0*q0*s2 + 2*k1*m2*p0*q0*s2 + 
  2*k2*m0*p1*q0*s2 + 2*k0*m2*p1*q0*s2 + 2*k1*m0*p2*q0*s2 + 
  2*k0*m1*p2*q0*s2 - 2*k2*l1*q0*q0*s2 - 2*k1*l2*q0*q0*s2 + 
  2*l2*m0*o0*q1*s2 + 2*l0*m2*o0*q1*s2 + 2*l0*m0*o2*q1*s2 + 
  2*k2*m0*p0*q1*s2 + 2*k0*m2*p0*q1*s2 + 2*k0*m0*p2*q1*s2 - 
  4*k2*l0*q0*q1*s2 - 4*k0*l2*q0*q1*s2 + 2*l1*m0*o0*q2*s2 + 
  2*l0*m1*o0*q2*s2 + 2*l0*m0*o1*q2*s2 + 2*k1*m0*p0*q2*s2 + 
  2*k0*m1*p0*q2*s2 + 2*k0*m0*p1*q2*s2 - 4*k1*l0*q0*q2*s2 - 
  4*k0*l1*q0*q2*s2 - 4*k0*l0*q1*q2*s2 + 4*m1*m2*n0*s0*s2 + 
  4*m0*m2*n1*s0*s2 + 4*m0*m1*n2*s0*s2 - 4*j2*m1*q0*s0*s2 - 
  4*j1*m2*q0*s0*s2 - 4*j2*m0*q1*s0*s2 - 4*j0*m2*q1*s0*s2 + 
  4*i2*q0*q1*s0*s2 - 4*j1*m0*q2*s0*s2 - 4*j0*m1*q2*s0*s2 + 
  4*i1*q0*q2*s0*s2 + 4*i0*q1*q2*s0*s2 + 4*m0*m2*n0*s1*s2 + 
  2*m0*m0*n2*s1*s2 - 4*j2*m0*q0*s1*s2 - 4*j0*m2*q0*s1*s2 + 
  2*i2*q0*q0*s1*s2 - 4*j0*m0*q2*s1*s2 + 4*i0*q0*q2*s1*s2 + 
  2*m0*m1*n0*s2*s2 + m0*m0*n1*s2*s2 - 2*j1*m0*q0*s2*s2 - 
  2*j0*m1*q0*s2*s2 + i1*q0*q0*s2*s2 - 2*j0*m0*q1*s2*s2 + 
  2*i0*q0*q1*s2*s2 + 2*l2*m2*o1*p0*t0 + 2*l2*m1*o2*p0*t0 + 
  2*l1*m2*o2*p0*t0 + 2*l2*m2*o0*p1*t0 + 2*l2*m0*o2*p1*t0 + 
  2*l0*m2*o2*p1*t0 - 4*k2*m2*p0*p1*t0 + 2*l2*m1*o0*p2*t0 + 
  2*l1*m2*o0*p2*t0 + 2*l2*m0*o1*p2*t0 + 2*l0*m2*o1*p2*t0 + 
  2*l1*m0*o2*p2*t0 + 2*l0*m1*o2*p2*t0 - 4*k2*m1*p0*p2*t0 - 
  4*k1*m2*p0*p2*t0 - 4*k2*m0*p1*p2*t0 - 4*k0*m2*p1*p2*t0 - 
  2*k1*m0*p2*p2*t0 - 2*k0*m1*p2*p2*t0 - 2*l2*l2*o1*q0*t0 - 
  4*l1*l2*o2*q0*t0 + 2*k2*l2*p1*q0*t0 + 2*k2*l1*p2*q0*t0 + 
  2*k1*l2*p2*q0*t0 - 2*l2*l2*o0*q1*t0 - 4*l0*l2*o2*q1*t0 + 
  2*k2*l2*p0*q1*t0 + 2*k2*l0*p2*q1*t0 + 2*k0*l2*p2*q1*t0 - 
  4*l1*l2*o0*q2*t0 - 4*l0*l2*o1*q2*t0 - 4*l0*l1*o2*q2*t0 + 
  2*k2*l1*p0*q2*t0 + 2*k1*l2*p0*q2*t0 + 2*k2*l0*p1*q2*t0 + 
  2*k0*l2*p1*q2*t0 + 2*k1*l0*p2*q2*t0 + 2*k0*l1*p2*q2*t0 - 
  2*l2*m2*n1*s0*t0 - 2*l2*m1*n2*s0*t0 - 2*l1*m2*n2*s0*t0 + 
  2*j2*m2*p1*s0*t0 + 2*j2*m1*p2*s0*t0 + 2*j1*m2*p2*s0*t0 + 
  2*j2*l2*q1*s0*t0 - 2*i2*p2*q1*s0*t0 + 2*j2*l1*q2*s0*t0 + 
  2*j1*l2*q2*s0*t0 - 2*i2*p1*q2*s0*t0 - 2*i1*p2*q2*s0*t0 - 
  2*l2*m2*n0*s1*t0 - 2*l2*m0*n2*s1*t0 - 2*l0*m2*n2*s1*t0 + 
  2*j2*m2*p0*s1*t0 + 2*j2*m0*p2*s1*t0 + 2*j0*m2*p2*s1*t0 + 
  2*j2*l2*q0*s1*t0 - 2*i2*p2*q0*s1*t0 + 2*j2*l0*q2*s1*t0 + 
  2*j0*l2*q2*s1*t0 - 2*i2*p0*q2*s1*t0 - 2*i0*p2*q2*s1*t0 - 
  2*l2*m1*n0*s2*t0 - 2*l1*m2*n0*s2*t0 - 2*l2*m0*n1*s2*t0 - 
  2*l0*m2*n1*s2*t0 - 2*l1*m0*n2*s2*t0 - 2*l0*m1*n2*s2*t0 + 
  2*j2*m1*p0*s2*t0 + 2*j1*m2*p0*s2*t0 + 2*j2*m0*p1*s2*t0 + 
  2*j0*m2*p1*s2*t0 + 2*j1*m0*p2*s2*t0 + 2*j0*m1*p2*s2*t0 + 
  2*j2*l1*q0*s2*t0 + 2*j1*l2*q0*s2*t0 - 2*i2*p1*q0*s2*t0 - 
  2*i1*p2*q0*s2*t0 + 2*j2*l0*q1*s2*t0 + 2*j0*l2*q1*s2*t0 - 
  2*i2*p0*q1*s2*t0 - 2*i0*p2*q1*s2*t0 + 2*j1*l0*q2*s2*t0 + 
  2*j0*l1*q2*s2*t0 - 2*i1*p0*q2*s2*t0 - 2*i0*p1*q2*s2*t0 + 
  l2*l2*n1*t0*t0 + 2*l1*l2*n2*t0*t0 - 2*j2*l2*p1*t0*t0 - 
  2*j2*l1*p2*t0*t0 - 2*j1*l2*p2*t0*t0 + 2*i2*p1*p2*t0*t0 + 
  i1*p2*p2*t0*t0 + 2*l2*m2*o0*p0*t1 + 2*l2*m0*o2*p0*t1 + 
  2*l0*m2*o2*p0*t1 - 2*k2*m2*p0*p0*t1 + 2*l2*m0*o0*p2*t1 + 
  2*l0*m2*o0*p2*t1 + 2*l0*m0*o2*p2*t1 - 4*k2*m0*p0*p2*t1 - 
  4*k0*m2*p0*p2*t1 - 2*k0*m0*p2*p2*t1 - 2*l2*l2*o0*q0*t1 - 
  4*l0*l2*o2*q0*t1 + 2*k2*l2*p0*q0*t1 + 2*k2*l0*p2*q0*t1 + 
  2*k0*l2*p2*q0*t1 - 4*l0*l2*o0*q2*t1 - 2*l0*l0*o2*q2*t1 + 
  2*k2*l0*p0*q2*t1 + 2*k0*l2*p0*q2*t1 + 2*k0*l0*p2*q2*t1 - 
  2*l2*m2*n0*s0*t1 - 2*l2*m0*n2*s0*t1 - 2*l0*m2*n2*s0*t1 + 
  2*j2*m2*p0*s0*t1 + 2*j2*m0*p2*s0*t1 + 2*j0*m2*p2*s0*t1 + 
  2*j2*l2*q0*s0*t1 - 2*i2*p2*q0*s0*t1 + 2*j2*l0*q2*s0*t1 + 
  2*j0*l2*q2*s0*t1 - 2*i2*p0*q2*s0*t1 - 2*i0*p2*q2*s0*t1 - 
  2*l2*m0*n0*s2*t1 - 2*l0*m2*n0*s2*t1 - 2*l0*m0*n2*s2*t1 + 
  2*j2*m0*p0*s2*t1 + 2*j0*m2*p0*s2*t1 + 2*j0*m0*p2*s2*t1 + 
  2*j2*l0*q0*s2*t1 + 2*j0*l2*q0*s2*t1 - 2*i2*p0*q0*s2*t1 - 
  2*i0*p2*q0*s2*t1 + 2*j0*l0*q2*s2*t1 - 2*i0*p0*q2*s2*t1 + 
  2*l2*l2*n0*t0*t1 + 4*l0*l2*n2*t0*t1 - 4*j2*l2*p0*t0*t1 - 
  4*j2*l0*p2*t0*t1 - 4*j0*l2*p2*t0*t1 + 4*i2*p0*p2*t0*t1 + 
  2*i0*p2*p2*t0*t1 + 2*l2*m1*o0*p0*t2 + 2*l1*m2*o0*p0*t2 + 
  2*l2*m0*o1*p0*t2 + 2*l0*m2*o1*p0*t2 + 2*l1*m0*o2*p0*t2 + 
  2*l0*m1*o2*p0*t2 - 2*k2*m1*p0*p0*t2 - 2*k1*m2*p0*p0*t2 + 
  2*l2*m0*o0*p1*t2 + 2*l0*m2*o0*p1*t2 + 2*l0*m0*o2*p1*t2 - 
  4*k2*m0*p0*p1*t2 - 4*k0*m2*p0*p1*t2 + 2*l1*m0*o0*p2*t2 + 
  2*l0*m1*o0*p2*t2 + 2*l0*m0*o1*p2*t2 - 4*k1*m0*p0*p2*t2 - 
  4*k0*m1*p0*p2*t2 - 4*k0*m0*p1*p2*t2 - 4*l1*l2*o0*q0*t2 - 
  4*l0*l2*o1*q0*t2 - 4*l0*l1*o2*q0*t2 + 2*k2*l1*p0*q0*t2 + 
  2*k1*l2*p0*q0*t2 + 2*k2*l0*p1*q0*t2 + 2*k0*l2*p1*q0*t2 + 
  2*k1*l0*p2*q0*t2 + 2*k0*l1*p2*q0*t2 - 4*l0*l2*o0*q1*t2 - 
  2*l0*l0*o2*q1*t2 + 2*k2*l0*p0*q1*t2 + 2*k0*l2*p0*q1*t2 + 
  2*k0*l0*p2*q1*t2 - 4*l0*l1*o0*q2*t2 - 2*l0*l0*o1*q2*t2 + 
  2*k1*l0*p0*q2*t2 + 2*k0*l1*p0*q2*t2 + 2*k0*l0*p1*q2*t2 - 
  2*l2*m1*n0*s0*t2 - 2*l1*m2*n0*s0*t2 - 2*l2*m0*n1*s0*t2 - 
  2*l0*m2*n1*s0*t2 - 2*l1*m0*n2*s0*t2 - 2*l0*m1*n2*s0*t2 + 
  2*j2*m1*p0*s0*t2 + 2*j1*m2*p0*s0*t2 + 2*j2*m0*p1*s0*t2 + 
  2*j0*m2*p1*s0*t2 + 2*j1*m0*p2*s0*t2 + 2*j0*m1*p2*s0*t2 + 
  2*j2*l1*q0*s0*t2 + 2*j1*l2*q0*s0*t2 - 2*i2*p1*q0*s0*t2 - 
  2*i1*p2*q0*s0*t2 + 2*j2*l0*q1*s0*t2 + 2*j0*l2*q1*s0*t2 - 
  2*i2*p0*q1*s0*t2 - 2*i0*p2*q1*s0*t2 + 2*j1*l0*q2*s0*t2 + 
  2*j0*l1*q2*s0*t2 - 2*i1*p0*q2*s0*t2 - 2*i0*p1*q2*s0*t2 - 
  2*l2*m0*n0*s1*t2 - 2*l0*m2*n0*s1*t2 - 2*l0*m0*n2*s1*t2 + 
  2*j2*m0*p0*s1*t2 + 2*j0*m2*p0*s1*t2 + 2*j0*m0*p2*s1*t2 + 
  2*j2*l0*q0*s1*t2 + 2*j0*l2*q0*s1*t2 - 2*i2*p0*q0*s1*t2 - 
  2*i0*p2*q0*s1*t2 + 2*j0*l0*q2*s1*t2 - 2*i0*p0*q2*s1*t2 - 
  2*l1*m0*n0*s2*t2 - 2*l0*m1*n0*s2*t2 - 2*l0*m0*n1*s2*t2 + 
  2*j1*m0*p0*s2*t2 + 2*j0*m1*p0*s2*t2 + 2*j0*m0*p1*s2*t2 + 
  2*j1*l0*q0*s2*t2 + 2*j0*l1*q0*s2*t2 - 2*i1*p0*q0*s2*t2 - 
  2*i0*p1*q0*s2*t2 + 2*j0*l0*q1*s2*t2 - 2*i0*p0*q1*s2*t2 + 
  4*l1*l2*n0*t0*t2 + 4*l0*l2*n1*t0*t2 + 4*l0*l1*n2*t0*t2 - 
  4*j2*l1*p0*t0*t2 - 4*j1*l2*p0*t0*t2 - 4*j2*l0*p1*t0*t2 - 
  4*j0*l2*p1*t0*t2 + 4*i2*p0*p1*t0*t2 - 4*j1*l0*p2*t0*t2 - 
  4*j0*l1*p2*t0*t2 + 4*i1*p0*p2*t0*t2 + 4*i0*p1*p2*t0*t2 + 
  4*l0*l2*n0*t1*t2 + 2*l0*l0*n2*t1*t2 - 4*j2*l0*p0*t1*t2 - 
  4*j0*l2*p0*t1*t2 + 2*i2*p0*p0*t1*t2 - 4*j0*l0*p2*t1*t2 + 
  4*i0*p0*p2*t1*t2 + 2*l0*l1*n0*t2*t2 + l0*l0*n1*t2*t2 - 
  2*j1*l0*p0*t2*t2 - 2*j0*l1*p0*t2*t2 + i1*p0*p0*t2*t2 - 
  2*j0*l0*p1*t2*t2 + 2*i0*p0*p1*t2*t2 + 2*m2*m2*o0*o1*u0 + 
  4*m1*m2*o0*o2*u0 + 4*m0*m2*o1*o2*u0 + 2*m0*m1*o2*o2*u0 - 
  2*k2*m2*o1*q0*u0 - 2*k2*m1*o2*q0*u0 - 2*k1*m2*o2*q0*u0 - 
  2*k2*m2*o0*q1*u0 - 2*k2*m0*o2*q1*u0 - 2*k0*m2*o2*q1*u0 + 
  2*k2*k2*q0*q1*u0 - 2*k2*m1*o0*q2*u0 - 2*k1*m2*o0*q2*u0 - 
  2*k2*m0*o1*q2*u0 - 2*k0*m2*o1*q2*u0 - 2*k1*m0*o2*q2*u0 - 
  2*k0*m1*o2*q2*u0 + 4*k1*k2*q0*q2*u0 + 4*k0*k2*q1*q2*u0 + 
  2*k0*k1*q2*q2*u0 - m2*m2*n1*r0*u0 - 2*m1*m2*n2*r0*u0 + 
  2*j2*m2*q1*r0*u0 + 2*j2*m1*q2*r0*u0 + 2*j1*m2*q2*r0*u0 - 
  2*i2*q1*q2*r0*u0 - i1*q2*q2*r0*u0 - m2*m2*n0*r1*u0 - 
  2*m0*m2*n2*r1*u0 + 2*j2*m2*q0*r1*u0 + 2*j2*m0*q2*r1*u0 + 
  2*j0*m2*q2*r1*u0 - 2*i2*q0*q2*r1*u0 - i0*q2*q2*r1*u0 - 
  2*m1*m2*n0*r2*u0 - 2*m0*m2*n1*r2*u0 - 2*m0*m1*n2*r2*u0 + 
  2*j2*m1*q0*r2*u0 + 2*j1*m2*q0*r2*u0 + 2*j2*m0*q1*r2*u0 + 
  2*j0*m2*q1*r2*u0 - 2*i2*q0*q1*r2*u0 + 2*j1*m0*q2*r2*u0 + 
  2*j0*m1*q2*r2*u0 - 2*i1*q0*q2*r2*u0 - 2*i0*q1*q2*r2*u0 + 
  2*k2*m2*n1*t0*u0 + 2*k2*m1*n2*t0*u0 + 2*k1*m2*n2*t0*u0 - 
  2*j2*m2*o1*t0*u0 - 2*j2*m1*o2*t0*u0 - 2*j1*m2*o2*t0*u0 - 
  2*j2*k2*q1*t0*u0 + 2*i2*o2*q1*t0*u0 - 2*j2*k1*q2*t0*u0 - 
  2*j1*k2*q2*t0*u0 + 2*i2*o1*q2*t0*u0 + 2*i1*o2*q2*t0*u0 + 
  2*k2*m2*n0*t1*u0 + 2*k2*m0*n2*t1*u0 + 2*k0*m2*n2*t1*u0 - 
  2*j2*m2*o0*t1*u0 - 2*j2*m0*o2*t1*u0 - 2*j0*m2*o2*t1*u0 - 
  2*j2*k2*q0*t1*u0 + 2*i2*o2*q0*t1*u0 - 2*j2*k0*q2*t1*u0 - 
  2*j0*k2*q2*t1*u0 + 2*i2*o0*q2*t1*u0 + 2*i0*o2*q2*t1*u0 + 
  2*j2*j2*t0*t1*u0 - 2*i2*n2*t0*t1*u0 + 2*k2*m1*n0*t2*u0 + 
  2*k1*m2*n0*t2*u0 + 2*k2*m0*n1*t2*u0 + 2*k0*m2*n1*t2*u0 + 
  2*k1*m0*n2*t2*u0 + 2*k0*m1*n2*t2*u0 - 2*j2*m1*o0*t2*u0 - 
  2*j1*m2*o0*t2*u0 - 2*j2*m0*o1*t2*u0 - 2*j0*m2*o1*t2*u0 - 
  2*j1*m0*o2*t2*u0 - 2*j0*m1*o2*t2*u0 - 2*j2*k1*q0*t2*u0 - 
  2*j1*k2*q0*t2*u0 + 2*i2*o1*q0*t2*u0 + 2*i1*o2*q0*t2*u0 - 
  2*j2*k0*q1*t2*u0 - 2*j0*k2*q1*t2*u0 + 2*i2*o0*q1*t2*u0 + 
  2*i0*o2*q1*t2*u0 - 2*j1*k0*q2*t2*u0 - 2*j0*k1*q2*t2*u0 + 
  2*i1*o0*q2*t2*u0 + 2*i0*o1*q2*t2*u0 + 4*j1*j2*t0*t2*u0 - 
  2*i2*n1*t0*t2*u0 - 2*i1*n2*t0*t2*u0 + 4*j0*j2*t1*t2*u0 - 
  2*i2*n0*t1*t2*u0 - 2*i0*n2*t1*t2*u0 + 2*j0*j1*t2*t2*u0 - 
  i1*n0*t2*t2*u0 - i0*n1*t2*t2*u0 + m2*m2*o0*o0*u1 + 
  4*m0*m2*o0*o2*u1 + m0*m0*o2*o2*u1 - 2*k2*m2*o0*q0*u1 - 
  2*k2*m0*o2*q0*u1 - 2*k0*m2*o2*q0*u1 + k2*k2*q0*q0*u1 - 
  2*k2*m0*o0*q2*u1 - 2*k0*m2*o0*q2*u1 - 2*k0*m0*o2*q2*u1 + 
  4*k0*k2*q0*q2*u1 + k0*k0*q2*q2*u1 - m2*m2*n0*r0*u1 - 
  2*m0*m2*n2*r0*u1 + 2*j2*m2*q0*r0*u1 + 2*j2*m0*q2*r0*u1 + 
  2*j0*m2*q2*r0*u1 - 2*i2*q0*q2*r0*u1 - i0*q2*q2*r0*u1 - 
  2*m0*m2*n0*r2*u1 - m0*m0*n2*r2*u1 + 2*j2*m0*q0*r2*u1 + 
  2*j0*m2*q0*r2*u1 - i2*q0*q0*r2*u1 + 2*j0*m0*q2*r2*u1 - 
  2*i0*q0*q2*r2*u1 + 2*k2*m2*n0*t0*u1 + 2*k2*m0*n2*t0*u1 + 
  2*k0*m2*n2*t0*u1 - 2*j2*m2*o0*t0*u1 - 2*j2*m0*o2*t0*u1 - 
  2*j0*m2*o2*t0*u1 - 2*j2*k2*q0*t0*u1 + 2*i2*o2*q0*t0*u1 - 
  2*j2*k0*q2*t0*u1 - 2*j0*k2*q2*t0*u1 + 2*i2*o0*q2*t0*u1 + 
  2*i0*o2*q2*t0*u1 + j2*j2*t0*t0*u1 - i2*n2*t0*t0*u1 + 
  2*k2*m0*n0*t2*u1 + 2*k0*m2*n0*t2*u1 + 2*k0*m0*n2*t2*u1 - 
  2*j2*m0*o0*t2*u1 - 2*j0*m2*o0*t2*u1 - 2*j0*m0*o2*t2*u1 - 
  2*j2*k0*q0*t2*u1 - 2*j0*k2*q0*t2*u1 + 2*i2*o0*q0*t2*u1 + 
  2*i0*o2*q0*t2*u1 - 2*j0*k0*q2*t2*u1 + 2*i0*o0*q2*t2*u1 + 
  4*j0*j2*t0*t2*u1 - 2*i2*n0*t0*t2*u1 - 2*i0*n2*t0*t2*u1 + 
  j0*j0*t2*t2*u1 - i0*n0*t2*t2*u1 + 2*m1*m2*o0*o0*u2 + 
  4*m0*m2*o0*o1*u2 + 4*m0*m1*o0*o2*u2 + 2*m0*m0*o1*o2*u2 - 
  2*k2*m1*o0*q0*u2 - 2*k1*m2*o0*q0*u2 - 2*k2*m0*o1*q0*u2 - 
  2*k0*m2*o1*q0*u2 - 2*k1*m0*o2*q0*u2 - 2*k0*m1*o2*q0*u2 + 
  2*k1*k2*q0*q0*u2 - 2*k2*m0*o0*q1*u2 - 2*k0*m2*o0*q1*u2 - 
  2*k0*m0*o2*q1*u2 + 4*k0*k2*q0*q1*u2 - 2*k1*m0*o0*q2*u2 - 
  2*k0*m1*o0*q2*u2 - 2*k0*m0*o1*q2*u2 + 4*k0*k1*q0*q2*u2 + 
  2*k0*k0*q1*q2*u2 - 2*m1*m2*n0*r0*u2 - 2*m0*m2*n1*r0*u2 - 
  2*m0*m1*n2*r0*u2 + 2*j2*m1*q0*r0*u2 + 2*j1*m2*q0*r0*u2 + 
  2*j2*m0*q1*r0*u2 + 2*j0*m2*q1*r0*u2 - 2*i2*q0*q1*r0*u2 + 
  2*j1*m0*q2*r0*u2 + 2*j0*m1*q2*r0*u2 - 2*i1*q0*q2*r0*u2 - 
  2*i0*q1*q2*r0*u2 - 2*m0*m2*n0*r1*u2 - m0*m0*n2*r1*u2 + 
  2*j2*m0*q0*r1*u2 + 2*j0*m2*q0*r1*u2 - i2*q0*q0*r1*u2 + 
  2*j0*m0*q2*r1*u2 - 2*i0*q0*q2*r1*u2 - 2*m0*m1*n0*r2*u2 - 
  m0*m0*n1*r2*u2 + 2*j1*m0*q0*r2*u2 + 2*j0*m1*q0*r2*u2 - 
  i1*q0*q0*r2*u2 + 2*j0*m0*q1*r2*u2 - 2*i0*q0*q1*r2*u2 + 
  2*k2*m1*n0*t0*u2 + 2*k1*m2*n0*t0*u2 + 2*k2*m0*n1*t0*u2 + 
  2*k0*m2*n1*t0*u2 + 2*k1*m0*n2*t0*u2 + 2*k0*m1*n2*t0*u2 - 
  2*j2*m1*o0*t0*u2 - 2*j1*m2*o0*t0*u2 - 2*j2*m0*o1*t0*u2 - 
  2*j0*m2*o1*t0*u2 - 2*j1*m0*o2*t0*u2 - 2*j0*m1*o2*t0*u2 - 
  2*j2*k1*q0*t0*u2 - 2*j1*k2*q0*t0*u2 + 2*i2*o1*q0*t0*u2 + 
  2*i1*o2*q0*t0*u2 - 2*j2*k0*q1*t0*u2 - 2*j0*k2*q1*t0*u2 + 
  2*i2*o0*q1*t0*u2 + 2*i0*o2*q1*t0*u2 - 2*j1*k0*q2*t0*u2 - 
  2*j0*k1*q2*t0*u2 + 2*i1*o0*q2*t0*u2 + 2*i0*o1*q2*t0*u2 + 
  2*j1*j2*t0*t0*u2 - i2*n1*t0*t0*u2 - i1*n2*t0*t0*u2 + 
  2*k2*m0*n0*t1*u2 + 2*k0*m2*n0*t1*u2 + 2*k0*m0*n2*t1*u2 - 
  2*j2*m0*o0*t1*u2 - 2*j0*m2*o0*t1*u2 - 2*j0*m0*o2*t1*u2 - 
  2*j2*k0*q0*t1*u2 - 2*j0*k2*q0*t1*u2 + 2*i2*o0*q0*t1*u2 + 
  2*i0*o2*q0*t1*u2 - 2*j0*k0*q2*t1*u2 + 2*i0*o0*q2*t1*u2 + 
  4*j0*j2*t0*t1*u2 - 2*i2*n0*t0*t1*u2 - 2*i0*n2*t0*t1*u2 + 
  2*k1*m0*n0*t2*u2 + 2*k0*m1*n0*t2*u2 + 2*k0*m0*n1*t2*u2 - 
  2*j1*m0*o0*t2*u2 - 2*j0*m1*o0*t2*u2 - 2*j0*m0*o1*t2*u2 - 
  2*j1*k0*q0*t2*u2 - 2*j0*k1*q0*t2*u2 + 2*i1*o0*q0*t2*u2 + 
  2*i0*o1*q0*t2*u2 - 2*j0*k0*q1*t2*u2 + 2*i0*o0*q1*t2*u2 + 
  4*j0*j1*t0*t2*u2 - 2*i1*n0*t0*t2*u2 - 2*i0*n1*t0*t2*u2 + 
  2*j0*j0*t1*t2*u2 - 2*i0*n0*t1*t2*u2 - 4*l2*m2*o0*o1*v0 - 
  4*l2*m1*o0*o2*v0 - 4*l1*m2*o0*o2*v0 - 4*l2*m0*o1*o2*v0 - 
  4*l0*m2*o1*o2*v0 - 2*l1*m0*o2*o2*v0 - 2*l0*m1*o2*o2*v0 + 
  2*k2*m2*o1*p0*v0 + 2*k2*m1*o2*p0*v0 + 2*k1*m2*o2*p0*v0 + 
  2*k2*m2*o0*p1*v0 + 2*k2*m0*o2*p1*v0 + 2*k0*m2*o2*p1*v0 + 
  2*k2*m1*o0*p2*v0 + 2*k1*m2*o0*p2*v0 + 2*k2*m0*o1*p2*v0 + 
  2*k0*m2*o1*p2*v0 + 2*k1*m0*o2*p2*v0 + 2*k0*m1*o2*p2*v0 + 
  2*k2*l2*o1*q0*v0 + 2*k2*l1*o2*q0*v0 + 2*k1*l2*o2*q0*v0 - 
  2*k2*k2*p1*q0*v0 - 4*k1*k2*p2*q0*v0 + 2*k2*l2*o0*q1*v0 + 
  2*k2*l0*o2*q1*v0 + 2*k0*l2*o2*q1*v0 - 2*k2*k2*p0*q1*v0 - 
  4*k0*k2*p2*q1*v0 + 2*k2*l1*o0*q2*v0 + 2*k1*l2*o0*q2*v0 + 
  2*k2*l0*o1*q2*v0 + 2*k0*l2*o1*q2*v0 + 2*k1*l0*o2*q2*v0 + 
  2*k0*l1*o2*q2*v0 - 4*k1*k2*p0*q2*v0 - 4*k0*k2*p1*q2*v0 - 
  4*k0*k1*p2*q2*v0 + 2*l2*m2*n1*r0*v0 + 2*l2*m1*n2*r0*v0 + 
  2*l1*m2*n2*r0*v0 - 2*j2*m2*p1*r0*v0 - 2*j2*m1*p2*r0*v0 - 
  2*j1*m2*p2*r0*v0 - 2*j2*l2*q1*r0*v0 + 2*i2*p2*q1*r0*v0 - 
  2*j2*l1*q2*r0*v0 - 2*j1*l2*q2*r0*v0 + 2*i2*p1*q2*r0*v0 + 
  2*i1*p2*q2*r0*v0 + 2*l2*m2*n0*r1*v0 + 2*l2*m0*n2*r1*v0 + 
  2*l0*m2*n2*r1*v0 - 2*j2*m2*p0*r1*v0 - 2*j2*m0*p2*r1*v0 - 
  2*j0*m2*p2*r1*v0 - 2*j2*l2*q0*r1*v0 + 2*i2*p2*q0*r1*v0 - 
  2*j2*l0*q2*r1*v0 - 2*j0*l2*q2*r1*v0 + 2*i2*p0*q2*r1*v0 + 
  2*i0*p2*q2*r1*v0 + 2*l2*m1*n0*r2*v0 + 2*l1*m2*n0*r2*v0 + 
  2*l2*m0*n1*r2*v0 + 2*l0*m2*n1*r2*v0 + 2*l1*m0*n2*r2*v0 + 
  2*l0*m1*n2*r2*v0 - 2*j2*m1*p0*r2*v0 - 2*j1*m2*p0*r2*v0 - 
  2*j2*m0*p1*r2*v0 - 2*j0*m2*p1*r2*v0 - 2*j1*m0*p2*r2*v0 - 
  2*j0*m1*p2*r2*v0 - 2*j2*l1*q0*r2*v0 - 2*j1*l2*q0*r2*v0 + 
  2*i2*p1*q0*r2*v0 + 2*i1*p2*q0*r2*v0 - 2*j2*l0*q1*r2*v0 - 
  2*j0*l2*q1*r2*v0 + 2*i2*p0*q1*r2*v0 + 2*i0*p2*q1*r2*v0 - 
  2*j1*l0*q2*r2*v0 - 2*j0*l1*q2*r2*v0 + 2*i1*p0*q2*r2*v0 + 
  2*i0*p1*q2*r2*v0 - 2*k2*m2*n1*s0*v0 - 2*k2*m1*n2*s0*v0 - 
  2*k1*m2*n2*s0*v0 + 2*j2*m2*o1*s0*v0 + 2*j2*m1*o2*s0*v0 + 
  2*j1*m2*o2*s0*v0 + 2*j2*k2*q1*s0*v0 - 2*i2*o2*q1*s0*v0 + 
  2*j2*k1*q2*s0*v0 + 2*j1*k2*q2*s0*v0 - 2*i2*o1*q2*s0*v0 - 
  2*i1*o2*q2*s0*v0 - 2*k2*m2*n0*s1*v0 - 2*k2*m0*n2*s1*v0 - 
  2*k0*m2*n2*s1*v0 + 2*j2*m2*o0*s1*v0 + 2*j2*m0*o2*s1*v0 + 
  2*j0*m2*o2*s1*v0 + 2*j2*k2*q0*s1*v0 - 2*i2*o2*q0*s1*v0 + 
  2*j2*k0*q2*s1*v0 + 2*j0*k2*q2*s1*v0 - 2*i2*o0*q2*s1*v0 - 
  2*i0*o2*q2*s1*v0 - 2*k2*m1*n0*s2*v0 - 2*k1*m2*n0*s2*v0 - 
  2*k2*m0*n1*s2*v0 - 2*k0*m2*n1*s2*v0 - 2*k1*m0*n2*s2*v0 - 
  2*k0*m1*n2*s2*v0 + 2*j2*m1*o0*s2*v0 + 2*j1*m2*o0*s2*v0 + 
  2*j2*m0*o1*s2*v0 + 2*j0*m2*o1*s2*v0 + 2*j1*m0*o2*s2*v0 + 
  2*j0*m1*o2*s2*v0 + 2*j2*k1*q0*s2*v0 + 2*j1*k2*q0*s2*v0 - 
  2*i2*o1*q0*s2*v0 - 2*i1*o2*q0*s2*v0 + 2*j2*k0*q1*s2*v0 + 
  2*j0*k2*q1*s2*v0 - 2*i2*o0*q1*s2*v0 - 2*i0*o2*q1*s2*v0 + 
  2*j1*k0*q2*s2*v0 + 2*j0*k1*q2*s2*v0 - 2*i1*o0*q2*s2*v0 - 
  2*i0*o1*q2*s2*v0 - 2*k2*l2*n1*t0*v0 - 2*k2*l1*n2*t0*v0 - 
  2*k1*l2*n2*t0*v0 + 2*j2*l2*o1*t0*v0 + 2*j2*l1*o2*t0*v0 + 
  2*j1*l2*o2*t0*v0 + 2*j2*k2*p1*t0*v0 - 2*i2*o2*p1*t0*v0 + 
  2*j2*k1*p2*t0*v0 + 2*j1*k2*p2*t0*v0 - 2*i2*o1*p2*t0*v0 - 
  2*i1*o2*p2*t0*v0 - 2*j2*j2*s1*t0*v0 + 2*i2*n2*s1*t0*v0 - 
  4*j1*j2*s2*t0*v0 + 2*i2*n1*s2*t0*v0 + 2*i1*n2*s2*t0*v0 - 
  2*k2*l2*n0*t1*v0 - 2*k2*l0*n2*t1*v0 - 2*k0*l2*n2*t1*v0 + 
  2*j2*l2*o0*t1*v0 + 2*j2*l0*o2*t1*v0 + 2*j0*l2*o2*t1*v0 + 
  2*j2*k2*p0*t1*v0 - 2*i2*o2*p0*t1*v0 + 2*j2*k0*p2*t1*v0 + 
  2*j0*k2*p2*t1*v0 - 2*i2*o0*p2*t1*v0 - 2*i0*o2*p2*t1*v0 - 
  2*j2*j2*s0*t1*v0 + 2*i2*n2*s0*t1*v0 - 4*j0*j2*s2*t1*v0 + 
  2*i2*n0*s2*t1*v0 + 2*i0*n2*s2*t1*v0 - 2*k2*l1*n0*t2*v0 - 
  2*k1*l2*n0*t2*v0 - 2*k2*l0*n1*t2*v0 - 2*k0*l2*n1*t2*v0 - 
  2*k1*l0*n2*t2*v0 - 2*k0*l1*n2*t2*v0 + 2*j2*l1*o0*t2*v0 + 
  2*j1*l2*o0*t2*v0 + 2*j2*l0*o1*t2*v0 + 2*j0*l2*o1*t2*v0 + 
  2*j1*l0*o2*t2*v0 + 2*j0*l1*o2*t2*v0 + 2*j2*k1*p0*t2*v0 + 
  2*j1*k2*p0*t2*v0 - 2*i2*o1*p0*t2*v0 - 2*i1*o2*p0*t2*v0 + 
  2*j2*k0*p1*t2*v0 + 2*j0*k2*p1*t2*v0 - 2*i2*o0*p1*t2*v0 - 
  2*i0*o2*p1*t2*v0 + 2*j1*k0*p2*t2*v0 + 2*j0*k1*p2*t2*v0 - 
  2*i1*o0*p2*t2*v0 - 2*i0*o1*p2*t2*v0 - 4*j1*j2*s0*t2*v0 + 
  2*i2*n1*s0*t2*v0 + 2*i1*n2*s0*t2*v0 - 4*j0*j2*s1*t2*v0 + 
  2*i2*n0*s1*t2*v0 + 2*i0*n2*s1*t2*v0 - 4*j0*j1*s2*t2*v0 + 
  2*i1*n0*s2*t2*v0 + 2*i0*n1*s2*t2*v0 + k2*k2*n1*v0*v0 + 
  2*k1*k2*n2*v0*v0 - 2*j2*k2*o1*v0*v0 - 2*j2*k1*o2*v0*v0 - 
  2*j1*k2*o2*v0*v0 + 2*i2*o1*o2*v0*v0 + i1*o2*o2*v0*v0 + 
  j2*j2*r1*v0*v0 - i2*n2*r1*v0*v0 + 2*j1*j2*r2*v0*v0 - 
  i2*n1*r2*v0*v0 - i1*n2*r2*v0*v0 - 2*l2*m2*o0*o0*v1 - 
  4*l2*m0*o0*o2*v1 - 4*l0*m2*o0*o2*v1 - 2*l0*m0*o2*o2*v1 + 
  2*k2*m2*o0*p0*v1 + 2*k2*m0*o2*p0*v1 + 2*k0*m2*o2*p0*v1 + 
  2*k2*m0*o0*p2*v1 + 2*k0*m2*o0*p2*v1 + 2*k0*m0*o2*p2*v1 + 
  2*k2*l2*o0*q0*v1 + 2*k2*l0*o2*q0*v1 + 2*k0*l2*o2*q0*v1 - 
  2*k2*k2*p0*q0*v1 - 4*k0*k2*p2*q0*v1 + 2*k2*l0*o0*q2*v1 + 
  2*k0*l2*o0*q2*v1 + 2*k0*l0*o2*q2*v1 - 4*k0*k2*p0*q2*v1 - 
  2*k0*k0*p2*q2*v1 + 2*l2*m2*n0*r0*v1 + 2*l2*m0*n2*r0*v1 + 
  2*l0*m2*n2*r0*v1 - 2*j2*m2*p0*r0*v1 - 2*j2*m0*p2*r0*v1 - 
  2*j0*m2*p2*r0*v1 - 2*j2*l2*q0*r0*v1 + 2*i2*p2*q0*r0*v1 - 
  2*j2*l0*q2*r0*v1 - 2*j0*l2*q2*r0*v1 + 2*i2*p0*q2*r0*v1 + 
  2*i0*p2*q2*r0*v1 + 2*l2*m0*n0*r2*v1 + 2*l0*m2*n0*r2*v1 + 
  2*l0*m0*n2*r2*v1 - 2*j2*m0*p0*r2*v1 - 2*j0*m2*p0*r2*v1 - 
  2*j0*m0*p2*r2*v1 - 2*j2*l0*q0*r2*v1 - 2*j0*l2*q0*r2*v1 + 
  2*i2*p0*q0*r2*v1 + 2*i0*p2*q0*r2*v1 - 2*j0*l0*q2*r2*v1 + 
  2*i0*p0*q2*r2*v1 - 2*k2*m2*n0*s0*v1 - 2*k2*m0*n2*s0*v1 - 
  2*k0*m2*n2*s0*v1 + 2*j2*m2*o0*s0*v1 + 2*j2*m0*o2*s0*v1 + 
  2*j0*m2*o2*s0*v1 + 2*j2*k2*q0*s0*v1 - 2*i2*o2*q0*s0*v1 + 
  2*j2*k0*q2*s0*v1 + 2*j0*k2*q2*s0*v1 - 2*i2*o0*q2*s0*v1 - 
  2*i0*o2*q2*s0*v1 - 2*k2*m0*n0*s2*v1 - 2*k0*m2*n0*s2*v1 - 
  2*k0*m0*n2*s2*v1 + 2*j2*m0*o0*s2*v1 + 2*j0*m2*o0*s2*v1 + 
  2*j0*m0*o2*s2*v1 + 2*j2*k0*q0*s2*v1 + 2*j0*k2*q0*s2*v1 - 
  2*i2*o0*q0*s2*v1 - 2*i0*o2*q0*s2*v1 + 2*j0*k0*q2*s2*v1 - 
  2*i0*o0*q2*s2*v1 - 2*k2*l2*n0*t0*v1 - 2*k2*l0*n2*t0*v1 - 
  2*k0*l2*n2*t0*v1 + 2*j2*l2*o0*t0*v1 + 2*j2*l0*o2*t0*v1 + 
  2*j0*l2*o2*t0*v1 + 2*j2*k2*p0*t0*v1 - 2*i2*o2*p0*t0*v1 + 
  2*j2*k0*p2*t0*v1 + 2*j0*k2*p2*t0*v1 - 2*i2*o0*p2*t0*v1 - 
  2*i0*o2*p2*t0*v1 - 2*j2*j2*s0*t0*v1 + 2*i2*n2*s0*t0*v1 - 
  4*j0*j2*s2*t0*v1 + 2*i2*n0*s2*t0*v1 + 2*i0*n2*s2*t0*v1 - 
  2*k2*l0*n0*t2*v1 - 2*k0*l2*n0*t2*v1 - 2*k0*l0*n2*t2*v1 + 
  2*j2*l0*o0*t2*v1 + 2*j0*l2*o0*t2*v1 + 2*j0*l0*o2*t2*v1 + 
  2*j2*k0*p0*t2*v1 + 2*j0*k2*p0*t2*v1 - 2*i2*o0*p0*t2*v1 - 
  2*i0*o2*p0*t2*v1 + 2*j0*k0*p2*t2*v1 - 2*i0*o0*p2*t2*v1 - 
  4*j0*j2*s0*t2*v1 + 2*i2*n0*s0*t2*v1 + 2*i0*n2*s0*t2*v1 - 
  2*j0*j0*s2*t2*v1 + 2*i0*n0*s2*t2*v1 + 2*k2*k2*n0*v0*v1 + 
  4*k0*k2*n2*v0*v1 - 4*j2*k2*o0*v0*v1 - 4*j2*k0*o2*v0*v1 - 
  4*j0*k2*o2*v0*v1 + 4*i2*o0*o2*v0*v1 + 2*i0*o2*o2*v0*v1 + 
  2*j2*j2*r0*v0*v1 - 2*i2*n2*r0*v0*v1 + 4*j0*j2*r2*v0*v1 - 
  2*i2*n0*r2*v0*v1 - 2*i0*n2*r2*v0*v1 - 2*l2*m1*o0*o0*v2 - 
  2*l1*m2*o0*o0*v2 - 4*l2*m0*o0*o1*v2 - 4*l0*m2*o0*o1*v2 - 
  4*l1*m0*o0*o2*v2 - 4*l0*m1*o0*o2*v2 - 4*l0*m0*o1*o2*v2 + 
  2*k2*m1*o0*p0*v2 + 2*k1*m2*o0*p0*v2 + 2*k2*m0*o1*p0*v2 + 
  2*k0*m2*o1*p0*v2 + 2*k1*m0*o2*p0*v2 + 2*k0*m1*o2*p0*v2 + 
  2*k2*m0*o0*p1*v2 + 2*k0*m2*o0*p1*v2 + 2*k0*m0*o2*p1*v2 + 
  2*k1*m0*o0*p2*v2 + 2*k0*m1*o0*p2*v2 + 2*k0*m0*o1*p2*v2 + 
  2*k2*l1*o0*q0*v2 + 2*k1*l2*o0*q0*v2 + 2*k2*l0*o1*q0*v2 + 
  2*k0*l2*o1*q0*v2 + 2*k1*l0*o2*q0*v2 + 2*k0*l1*o2*q0*v2 - 
  4*k1*k2*p0*q0*v2 - 4*k0*k2*p1*q0*v2 - 4*k0*k1*p2*q0*v2 + 
  2*k2*l0*o0*q1*v2 + 2*k0*l2*o0*q1*v2 + 2*k0*l0*o2*q1*v2 - 
  4*k0*k2*p0*q1*v2 - 2*k0*k0*p2*q1*v2 + 2*k1*l0*o0*q2*v2 + 
  2*k0*l1*o0*q2*v2 + 2*k0*l0*o1*q2*v2 - 4*k0*k1*p0*q2*v2 - 
  2*k0*k0*p1*q2*v2 + 2*l2*m1*n0*r0*v2 + 2*l1*m2*n0*r0*v2 + 
  2*l2*m0*n1*r0*v2 + 2*l0*m2*n1*r0*v2 + 2*l1*m0*n2*r0*v2 + 
  2*l0*m1*n2*r0*v2 - 2*j2*m1*p0*r0*v2 - 2*j1*m2*p0*r0*v2 - 
  2*j2*m0*p1*r0*v2 - 2*j0*m2*p1*r0*v2 - 2*j1*m0*p2*r0*v2 - 
  2*j0*m1*p2*r0*v2 - 2*j2*l1*q0*r0*v2 - 2*j1*l2*q0*r0*v2 + 
  2*i2*p1*q0*r0*v2 + 2*i1*p2*q0*r0*v2 - 2*j2*l0*q1*r0*v2 - 
  2*j0*l2*q1*r0*v2 + 2*i2*p0*q1*r0*v2 + 2*i0*p2*q1*r0*v2 - 
  2*j1*l0*q2*r0*v2 - 2*j0*l1*q2*r0*v2 + 2*i1*p0*q2*r0*v2 + 
  2*i0*p1*q2*r0*v2 + 2*l2*m0*n0*r1*v2 + 2*l0*m2*n0*r1*v2 + 
  2*l0*m0*n2*r1*v2 - 2*j2*m0*p0*r1*v2 - 2*j0*m2*p0*r1*v2 - 
  2*j0*m0*p2*r1*v2 - 2*j2*l0*q0*r1*v2 - 2*j0*l2*q0*r1*v2 + 
  2*i2*p0*q0*r1*v2 + 2*i0*p2*q0*r1*v2 - 2*j0*l0*q2*r1*v2 + 
  2*i0*p0*q2*r1*v2 + 2*l1*m0*n0*r2*v2 + 2*l0*m1*n0*r2*v2 + 
  2*l0*m0*n1*r2*v2 - 2*j1*m0*p0*r2*v2 - 2*j0*m1*p0*r2*v2 - 
  2*j0*m0*p1*r2*v2 - 2*j1*l0*q0*r2*v2 - 2*j0*l1*q0*r2*v2 + 
  2*i1*p0*q0*r2*v2 + 2*i0*p1*q0*r2*v2 - 2*j0*l0*q1*r2*v2 + 
  2*i0*p0*q1*r2*v2 - 2*k2*m1*n0*s0*v2 - 2*k1*m2*n0*s0*v2 - 
  2*k2*m0*n1*s0*v2 - 2*k0*m2*n1*s0*v2 - 2*k1*m0*n2*s0*v2 - 
  2*k0*m1*n2*s0*v2 + 2*j2*m1*o0*s0*v2 + 2*j1*m2*o0*s0*v2 + 
  2*j2*m0*o1*s0*v2 + 2*j0*m2*o1*s0*v2 + 2*j1*m0*o2*s0*v2 + 
  2*j0*m1*o2*s0*v2 + 2*j2*k1*q0*s0*v2 + 2*j1*k2*q0*s0*v2 - 
  2*i2*o1*q0*s0*v2 - 2*i1*o2*q0*s0*v2 + 2*j2*k0*q1*s0*v2 + 
  2*j0*k2*q1*s0*v2 - 2*i2*o0*q1*s0*v2 - 2*i0*o2*q1*s0*v2 + 
  2*j1*k0*q2*s0*v2 + 2*j0*k1*q2*s0*v2 - 2*i1*o0*q2*s0*v2 - 
  2*i0*o1*q2*s0*v2 - 2*k2*m0*n0*s1*v2 - 2*k0*m2*n0*s1*v2 - 
  2*k0*m0*n2*s1*v2 + 2*j2*m0*o0*s1*v2 + 2*j0*m2*o0*s1*v2 + 
  2*j0*m0*o2*s1*v2 + 2*j2*k0*q0*s1*v2 + 2*j0*k2*q0*s1*v2 - 
  2*i2*o0*q0*s1*v2 - 2*i0*o2*q0*s1*v2 + 2*j0*k0*q2*s1*v2 - 
  2*i0*o0*q2*s1*v2 - 2*k1*m0*n0*s2*v2 - 2*k0*m1*n0*s2*v2 - 
  2*k0*m0*n1*s2*v2 + 2*j1*m0*o0*s2*v2 + 2*j0*m1*o0*s2*v2 + 
  2*j0*m0*o1*s2*v2 + 2*j1*k0*q0*s2*v2 + 2*j0*k1*q0*s2*v2 - 
  2*i1*o0*q0*s2*v2 - 2*i0*o1*q0*s2*v2 + 2*j0*k0*q1*s2*v2 - 
  2*i0*o0*q1*s2*v2 - 2*k2*l1*n0*t0*v2 - 2*k1*l2*n0*t0*v2 - 
  2*k2*l0*n1*t0*v2 - 2*k0*l2*n1*t0*v2 - 2*k1*l0*n2*t0*v2 - 
  2*k0*l1*n2*t0*v2 + 2*j2*l1*o0*t0*v2 + 2*j1*l2*o0*t0*v2 + 
  2*j2*l0*o1*t0*v2 + 2*j0*l2*o1*t0*v2 + 2*j1*l0*o2*t0*v2 + 
  2*j0*l1*o2*t0*v2 + 2*j2*k1*p0*t0*v2 + 2*j1*k2*p0*t0*v2 - 
  2*i2*o1*p0*t0*v2 - 2*i1*o2*p0*t0*v2 + 2*j2*k0*p1*t0*v2 + 
  2*j0*k2*p1*t0*v2 - 2*i2*o0*p1*t0*v2 - 2*i0*o2*p1*t0*v2 + 
  2*j1*k0*p2*t0*v2 + 2*j0*k1*p2*t0*v2 - 2*i1*o0*p2*t0*v2 - 
  2*i0*o1*p2*t0*v2 - 4*j1*j2*s0*t0*v2 + 2*i2*n1*s0*t0*v2 + 
  2*i1*n2*s0*t0*v2 - 4*j0*j2*s1*t0*v2 + 2*i2*n0*s1*t0*v2 + 
  2*i0*n2*s1*t0*v2 - 4*j0*j1*s2*t0*v2 + 2*i1*n0*s2*t0*v2 + 
  2*i0*n1*s2*t0*v2 - 2*k2*l0*n0*t1*v2 - 2*k0*l2*n0*t1*v2 - 
  2*k0*l0*n2*t1*v2 + 2*j2*l0*o0*t1*v2 + 2*j0*l2*o0*t1*v2 + 
  2*j0*l0*o2*t1*v2 + 2*j2*k0*p0*t1*v2 + 2*j0*k2*p0*t1*v2 - 
  2*i2*o0*p0*t1*v2 - 2*i0*o2*p0*t1*v2 + 2*j0*k0*p2*t1*v2 - 
  2*i0*o0*p2*t1*v2 - 4*j0*j2*s0*t1*v2 + 2*i2*n0*s0*t1*v2 + 
  2*i0*n2*s0*t1*v2 - 2*j0*j0*s2*t1*v2 + 2*i0*n0*s2*t1*v2 - 
  2*k1*l0*n0*t2*v2 - 2*k0*l1*n0*t2*v2 - 2*k0*l0*n1*t2*v2 + 
  2*j1*l0*o0*t2*v2 + 2*j0*l1*o0*t2*v2 + 2*j0*l0*o1*t2*v2 + 
  2*j1*k0*p0*t2*v2 + 2*j0*k1*p0*t2*v2 - 2*i1*o0*p0*t2*v2 - 
  2*i0*o1*p0*t2*v2 + 2*j0*k0*p1*t2*v2 - 2*i0*o0*p1*t2*v2 - 
  4*j0*j1*s0*t2*v2 + 2*i1*n0*s0*t2*v2 + 2*i0*n1*s0*t2*v2 - 
  2*j0*j0*s1*t2*v2 + 2*i0*n0*s1*t2*v2 + 4*k1*k2*n0*v0*v2 + 
  4*k0*k2*n1*v0*v2 + 4*k0*k1*n2*v0*v2 - 4*j2*k1*o0*v0*v2 - 
  4*j1*k2*o0*v0*v2 - 4*j2*k0*o1*v0*v2 - 4*j0*k2*o1*v0*v2 + 
  4*i2*o0*o1*v0*v2 - 4*j1*k0*o2*v0*v2 - 4*j0*k1*o2*v0*v2 + 
  4*i1*o0*o2*v0*v2 + 4*i0*o1*o2*v0*v2 + 4*j1*j2*r0*v0*v2 - 
  2*i2*n1*r0*v0*v2 - 2*i1*n2*r0*v0*v2 + 4*j0*j2*r1*v0*v2 - 
  2*i2*n0*r1*v0*v2 - 2*i0*n2*r1*v0*v2 + 4*j0*j1*r2*v0*v2 - 
  2*i1*n0*r2*v0*v2 - 2*i0*n1*r2*v0*v2 + 4*k0*k2*n0*v1*v2 + 
  2*k0*k0*n2*v1*v2 - 4*j2*k0*o0*v1*v2 - 4*j0*k2*o0*v1*v2 + 
  2*i2*o0*o0*v1*v2 - 4*j0*k0*o2*v1*v2 + 4*i0*o0*o2*v1*v2 + 
  4*j0*j2*r0*v1*v2 - 2*i2*n0*r0*v1*v2 - 2*i0*n2*r0*v1*v2 + 
  2*j0*j0*r2*v1*v2 - 2*i0*n0*r2*v1*v2 + 2*k0*k1*n0*v2*v2 + 
  k0*k0*n1*v2*v2 - 2*j1*k0*o0*v2*v2 - 2*j0*k1*o0*v2*v2 + 
  i1*o0*o0*v2*v2 - 2*j0*k0*o1*v2*v2 + 2*i0*o0*o1*v2*v2 + 
  2*j0*j1*r0*v2*v2 - i1*n0*r0*v2*v2 - i0*n1*r0*v2*v2 + 
  j0*j0*r1*v2*v2 - i0*n0*r1*v2*v2 + 2*l2*l2*o0*o1*w0 + 
  4*l1*l2*o0*o2*w0 + 4*l0*l2*o1*o2*w0 + 2*l0*l1*o2*o2*w0 - 
  2*k2*l2*o1*p0*w0 - 2*k2*l1*o2*p0*w0 - 2*k1*l2*o2*p0*w0 - 
  2*k2*l2*o0*p1*w0 - 2*k2*l0*o2*p1*w0 - 2*k0*l2*o2*p1*w0 + 
  2*k2*k2*p0*p1*w0 - 2*k2*l1*o0*p2*w0 - 2*k1*l2*o0*p2*w0 - 
  2*k2*l0*o1*p2*w0 - 2*k0*l2*o1*p2*w0 - 2*k1*l0*o2*p2*w0 - 
  2*k0*l1*o2*p2*w0 + 4*k1*k2*p0*p2*w0 + 4*k0*k2*p1*p2*w0 + 
  2*k0*k1*p2*p2*w0 - l2*l2*n1*r0*w0 - 2*l1*l2*n2*r0*w0 + 
  2*j2*l2*p1*r0*w0 + 2*j2*l1*p2*r0*w0 + 2*j1*l2*p2*r0*w0 - 
  2*i2*p1*p2*r0*w0 - i1*p2*p2*r0*w0 - l2*l2*n0*r1*w0 - 
  2*l0*l2*n2*r1*w0 + 2*j2*l2*p0*r1*w0 + 2*j2*l0*p2*r1*w0 + 
  2*j0*l2*p2*r1*w0 - 2*i2*p0*p2*r1*w0 - i0*p2*p2*r1*w0 - 
  2*l1*l2*n0*r2*w0 - 2*l0*l2*n1*r2*w0 - 2*l0*l1*n2*r2*w0 + 
  2*j2*l1*p0*r2*w0 + 2*j1*l2*p0*r2*w0 + 2*j2*l0*p1*r2*w0 + 
  2*j0*l2*p1*r2*w0 - 2*i2*p0*p1*r2*w0 + 2*j1*l0*p2*r2*w0 + 
  2*j0*l1*p2*r2*w0 - 2*i1*p0*p2*r2*w0 - 2*i0*p1*p2*r2*w0 + 
  2*k2*l2*n1*s0*w0 + 2*k2*l1*n2*s0*w0 + 2*k1*l2*n2*s0*w0 - 
  2*j2*l2*o1*s0*w0 - 2*j2*l1*o2*s0*w0 - 2*j1*l2*o2*s0*w0 - 
  2*j2*k2*p1*s0*w0 + 2*i2*o2*p1*s0*w0 - 2*j2*k1*p2*s0*w0 - 
  2*j1*k2*p2*s0*w0 + 2*i2*o1*p2*s0*w0 + 2*i1*o2*p2*s0*w0 + 
  2*k2*l2*n0*s1*w0 + 2*k2*l0*n2*s1*w0 + 2*k0*l2*n2*s1*w0 - 
  2*j2*l2*o0*s1*w0 - 2*j2*l0*o2*s1*w0 - 2*j0*l2*o2*s1*w0 - 
  2*j2*k2*p0*s1*w0 + 2*i2*o2*p0*s1*w0 - 2*j2*k0*p2*s1*w0 - 
  2*j0*k2*p2*s1*w0 + 2*i2*o0*p2*s1*w0 + 2*i0*o2*p2*s1*w0 + 
  2*j2*j2*s0*s1*w0 - 2*i2*n2*s0*s1*w0 + 2*k2*l1*n0*s2*w0 + 
  2*k1*l2*n0*s2*w0 + 2*k2*l0*n1*s2*w0 + 2*k0*l2*n1*s2*w0 + 
  2*k1*l0*n2*s2*w0 + 2*k0*l1*n2*s2*w0 - 2*j2*l1*o0*s2*w0 - 
  2*j1*l2*o0*s2*w0 - 2*j2*l0*o1*s2*w0 - 2*j0*l2*o1*s2*w0 - 
  2*j1*l0*o2*s2*w0 - 2*j0*l1*o2*s2*w0 - 2*j2*k1*p0*s2*w0 - 
  2*j1*k2*p0*s2*w0 + 2*i2*o1*p0*s2*w0 + 2*i1*o2*p0*s2*w0 - 
  2*j2*k0*p1*s2*w0 - 2*j0*k2*p1*s2*w0 + 2*i2*o0*p1*s2*w0 + 
  2*i0*o2*p1*s2*w0 - 2*j1*k0*p2*s2*w0 - 2*j0*k1*p2*s2*w0 + 
  2*i1*o0*p2*s2*w0 + 2*i0*o1*p2*s2*w0 + 4*j1*j2*s0*s2*w0 - 
  2*i2*n1*s0*s2*w0 - 2*i1*n2*s0*s2*w0 + 4*j0*j2*s1*s2*w0 - 
  2*i2*n0*s1*s2*w0 - 2*i0*n2*s1*s2*w0 + 2*j0*j1*s2*s2*w0 - 
  i1*n0*s2*s2*w0 - i0*n1*s2*s2*w0 - k2*k2*n1*u0*w0 - 
  2*k1*k2*n2*u0*w0 + 2*j2*k2*o1*u0*w0 + 2*j2*k1*o2*u0*w0 + 
  2*j1*k2*o2*u0*w0 - 2*i2*o1*o2*u0*w0 - i1*o2*o2*u0*w0 - 
  j2*j2*r1*u0*w0 + i2*n2*r1*u0*w0 - 2*j1*j2*r2*u0*w0 + 
  i2*n1*r2*u0*w0 + i1*n2*r2*u0*w0 - k2*k2*n0*u1*w0 - 
  2*k0*k2*n2*u1*w0 + 2*j2*k2*o0*u1*w0 + 2*j2*k0*o2*u1*w0 + 
  2*j0*k2*o2*u1*w0 - 2*i2*o0*o2*u1*w0 - i0*o2*o2*u1*w0 - 
  j2*j2*r0*u1*w0 + i2*n2*r0*u1*w0 - 2*j0*j2*r2*u1*w0 + 
  i2*n0*r2*u1*w0 + i0*n2*r2*u1*w0 - 2*k1*k2*n0*u2*w0 - 
  2*k0*k2*n1*u2*w0 - 2*k0*k1*n2*u2*w0 + 2*j2*k1*o0*u2*w0 + 
  2*j1*k2*o0*u2*w0 + 2*j2*k0*o1*u2*w0 + 2*j0*k2*o1*u2*w0 - 
  2*i2*o0*o1*u2*w0 + 2*j1*k0*o2*u2*w0 + 2*j0*k1*o2*u2*w0 - 
  2*i1*o0*o2*u2*w0 - 2*i0*o1*o2*u2*w0 - 2*j1*j2*r0*u2*w0 + 
  i2*n1*r0*u2*w0 + i1*n2*r0*u2*w0 - 2*j0*j2*r1*u2*w0 + 
  i2*n0*r1*u2*w0 + i0*n2*r1*u2*w0 - 2*j0*j1*r2*u2*w0 + 
  i1*n0*r2*u2*w0 + i0*n1*r2*u2*w0 + l2*l2*o0*o0*w1 + 
  4*l0*l2*o0*o2*w1 + l0*l0*o2*o2*w1 - 2*k2*l2*o0*p0*w1 - 
  2*k2*l0*o2*p0*w1 - 2*k0*l2*o2*p0*w1 + k2*k2*p0*p0*w1 - 
  2*k2*l0*o0*p2*w1 - 2*k0*l2*o0*p2*w1 - 2*k0*l0*o2*p2*w1 + 
  4*k0*k2*p0*p2*w1 + k0*k0*p2*p2*w1 - l2*l2*n0*r0*w1 - 
  2*l0*l2*n2*r0*w1 + 2*j2*l2*p0*r0*w1 + 2*j2*l0*p2*r0*w1 + 
  2*j0*l2*p2*r0*w1 - 2*i2*p0*p2*r0*w1 - i0*p2*p2*r0*w1 - 
  2*l0*l2*n0*r2*w1 - l0*l0*n2*r2*w1 + 2*j2*l0*p0*r2*w1 + 
  2*j0*l2*p0*r2*w1 - i2*p0*p0*r2*w1 + 2*j0*l0*p2*r2*w1 - 
  2*i0*p0*p2*r2*w1 + 2*k2*l2*n0*s0*w1 + 2*k2*l0*n2*s0*w1 + 
  2*k0*l2*n2*s0*w1 - 2*j2*l2*o0*s0*w1 - 2*j2*l0*o2*s0*w1 - 
  2*j0*l2*o2*s0*w1 - 2*j2*k2*p0*s0*w1 + 2*i2*o2*p0*s0*w1 - 
  2*j2*k0*p2*s0*w1 - 2*j0*k2*p2*s0*w1 + 2*i2*o0*p2*s0*w1 + 
  2*i0*o2*p2*s0*w1 + j2*j2*s0*s0*w1 - i2*n2*s0*s0*w1 + 
  2*k2*l0*n0*s2*w1 + 2*k0*l2*n0*s2*w1 + 2*k0*l0*n2*s2*w1 - 
  2*j2*l0*o0*s2*w1 - 2*j0*l2*o0*s2*w1 - 2*j0*l0*o2*s2*w1 - 
  2*j2*k0*p0*s2*w1 - 2*j0*k2*p0*s2*w1 + 2*i2*o0*p0*s2*w1 + 
  2*i0*o2*p0*s2*w1 - 2*j0*k0*p2*s2*w1 + 2*i0*o0*p2*s2*w1 + 
  4*j0*j2*s0*s2*w1 - 2*i2*n0*s0*s2*w1 - 2*i0*n2*s0*s2*w1 + 
  j0*j0*s2*s2*w1 - i0*n0*s2*s2*w1 - k2*k2*n0*u0*w1 - 
  2*k0*k2*n2*u0*w1 + 2*j2*k2*o0*u0*w1 + 2*j2*k0*o2*u0*w1 + 
  2*j0*k2*o2*u0*w1 - 2*i2*o0*o2*u0*w1 - i0*o2*o2*u0*w1 - 
  j2*j2*r0*u0*w1 + i2*n2*r0*u0*w1 - 2*j0*j2*r2*u0*w1 + 
  i2*n0*r2*u0*w1 + i0*n2*r2*u0*w1 - 2*k0*k2*n0*u2*w1 - 
  k0*k0*n2*u2*w1 + 2*j2*k0*o0*u2*w1 + 2*j0*k2*o0*u2*w1 - 
  i2*o0*o0*u2*w1 + 2*j0*k0*o2*u2*w1 - 2*i0*o0*o2*u2*w1 - 
  2*j0*j2*r0*u2*w1 + i2*n0*r0*u2*w1 + i0*n2*r0*u2*w1 - 
  j0*j0*r2*u2*w1 + i0*n0*r2*u2*w1 + 2*l1*l2*o0*o0*w2 + 
  4*l0*l2*o0*o1*w2 + 4*l0*l1*o0*o2*w2 + 2*l0*l0*o1*o2*w2 - 
  2*k2*l1*o0*p0*w2 - 2*k1*l2*o0*p0*w2 - 2*k2*l0*o1*p0*w2 - 
  2*k0*l2*o1*p0*w2 - 2*k1*l0*o2*p0*w2 - 2*k0*l1*o2*p0*w2 + 
  2*k1*k2*p0*p0*w2 - 2*k2*l0*o0*p1*w2 - 2*k0*l2*o0*p1*w2 - 
  2*k0*l0*o2*p1*w2 + 4*k0*k2*p0*p1*w2 - 2*k1*l0*o0*p2*w2 - 
  2*k0*l1*o0*p2*w2 - 2*k0*l0*o1*p2*w2 + 4*k0*k1*p0*p2*w2 + 
  2*k0*k0*p1*p2*w2 - 2*l1*l2*n0*r0*w2 - 2*l0*l2*n1*r0*w2 - 
  2*l0*l1*n2*r0*w2 + 2*j2*l1*p0*r0*w2 + 2*j1*l2*p0*r0*w2 + 
  2*j2*l0*p1*r0*w2 + 2*j0*l2*p1*r0*w2 - 2*i2*p0*p1*r0*w2 + 
  2*j1*l0*p2*r0*w2 + 2*j0*l1*p2*r0*w2 - 2*i1*p0*p2*r0*w2 - 
  2*i0*p1*p2*r0*w2 - 2*l0*l2*n0*r1*w2 - l0*l0*n2*r1*w2 + 
  2*j2*l0*p0*r1*w2 + 2*j0*l2*p0*r1*w2 - i2*p0*p0*r1*w2 + 
  2*j0*l0*p2*r1*w2 - 2*i0*p0*p2*r1*w2 - 2*l0*l1*n0*r2*w2 - 
  l0*l0*n1*r2*w2 + 2*j1*l0*p0*r2*w2 + 2*j0*l1*p0*r2*w2 - 
  i1*p0*p0*r2*w2 + 2*j0*l0*p1*r2*w2 - 2*i0*p0*p1*r2*w2 + 
  2*k2*l1*n0*s0*w2 + 2*k1*l2*n0*s0*w2 + 2*k2*l0*n1*s0*w2 + 
  2*k0*l2*n1*s0*w2 + 2*k1*l0*n2*s0*w2 + 2*k0*l1*n2*s0*w2 - 
  2*j2*l1*o0*s0*w2 - 2*j1*l2*o0*s0*w2 - 2*j2*l0*o1*s0*w2 - 
  2*j0*l2*o1*s0*w2 - 2*j1*l0*o2*s0*w2 - 2*j0*l1*o2*s0*w2 - 
  2*j2*k1*p0*s0*w2 - 2*j1*k2*p0*s0*w2 + 2*i2*o1*p0*s0*w2 + 
  2*i1*o2*p0*s0*w2 - 2*j2*k0*p1*s0*w2 - 2*j0*k2*p1*s0*w2 + 
  2*i2*o0*p1*s0*w2 + 2*i0*o2*p1*s0*w2 - 2*j1*k0*p2*s0*w2 - 
  2*j0*k1*p2*s0*w2 + 2*i1*o0*p2*s0*w2 + 2*i0*o1*p2*s0*w2 + 
  2*j1*j2*s0*s0*w2 - i2*n1*s0*s0*w2 - i1*n2*s0*s0*w2 + 
  2*k2*l0*n0*s1*w2 + 2*k0*l2*n0*s1*w2 + 2*k0*l0*n2*s1*w2 - 
  2*j2*l0*o0*s1*w2 - 2*j0*l2*o0*s1*w2 - 2*j0*l0*o2*s1*w2 - 
  2*j2*k0*p0*s1*w2 - 2*j0*k2*p0*s1*w2 + 2*i2*o0*p0*s1*w2 + 
  2*i0*o2*p0*s1*w2 - 2*j0*k0*p2*s1*w2 + 2*i0*o0*p2*s1*w2 + 
  4*j0*j2*s0*s1*w2 - 2*i2*n0*s0*s1*w2 - 2*i0*n2*s0*s1*w2 + 
  2*k1*l0*n0*s2*w2 + 2*k0*l1*n0*s2*w2 + 2*k0*l0*n1*s2*w2 - 
  2*j1*l0*o0*s2*w2 - 2*j0*l1*o0*s2*w2 - 2*j0*l0*o1*s2*w2 - 
  2*j1*k0*p0*s2*w2 - 2*j0*k1*p0*s2*w2 + 2*i1*o0*p0*s2*w2 + 
  2*i0*o1*p0*s2*w2 - 2*j0*k0*p1*s2*w2 + 2*i0*o0*p1*s2*w2 + 
  4*j0*j1*s0*s2*w2 - 2*i1*n0*s0*s2*w2 - 2*i0*n1*s0*s2*w2 + 
  2*j0*j0*s1*s2*w2 - 2*i0*n0*s1*s2*w2 - 2*k1*k2*n0*u0*w2 - 
  2*k0*k2*n1*u0*w2 - 2*k0*k1*n2*u0*w2 + 2*j2*k1*o0*u0*w2 + 
  2*j1*k2*o0*u0*w2 + 2*j2*k0*o1*u0*w2 + 2*j0*k2*o1*u0*w2 - 
  2*i2*o0*o1*u0*w2 + 2*j1*k0*o2*u0*w2 + 2*j0*k1*o2*u0*w2 - 
  2*i1*o0*o2*u0*w2 - 2*i0*o1*o2*u0*w2 - 2*j1*j2*r0*u0*w2 + 
  i2*n1*r0*u0*w2 + i1*n2*r0*u0*w2 - 2*j0*j2*r1*u0*w2 + 
  i2*n0*r1*u0*w2 + i0*n2*r1*u0*w2 - 2*j0*j1*r2*u0*w2 + 
  i1*n0*r2*u0*w2 + i0*n1*r2*u0*w2 - 2*k0*k2*n0*u1*w2 - 
  k0*k0*n2*u1*w2 + 2*j2*k0*o0*u1*w2 + 2*j0*k2*o0*u1*w2 - 
  i2*o0*o0*u1*w2 + 2*j0*k0*o2*u1*w2 - 2*i0*o0*o2*u1*w2 - 
  2*j0*j2*r0*u1*w2 + i2*n0*r0*u1*w2 + i0*n2*r0*u1*w2 - 
  j0*j0*r2*u1*w2 + i0*n0*r2*u1*w2 - 2*k0*k1*n0*u2*w2 - 
  k0*k0*n1*u2*w2 + 2*j1*k0*o0*u2*w2 + 2*j0*k1*o0*u2*w2 - 
  i1*o0*o0*u2*w2 + 2*j0*k0*o1*u2*w2 - 2*i0*o0*o1*u2*w2 - 
  2*j0*j1*r0*u2*w2 + i1*n0*r0*u2*w2 + i0*n1*r0*u2*w2 - 
  j0*j0*r1*u2*w2 + i0*n0*r1*u2*w2;
    
    k[8]  =   m2*m2*p1*p1*r0 + 4*m1*m2*p1*p2*r0 + m1*m1*p2*p2*r0 - 
  2*l2*m2*p1*q1*r0 - 2*l2*m1*p2*q1*r0 - 2*l1*m2*p2*q1*r0 + 
  l2*l2*q1*q1*r0 - 2*l2*m1*p1*q2*r0 - 2*l1*m2*p1*q2*r0 - 
  2*l1*m1*p2*q2*r0 + 4*l1*l2*q1*q2*r0 + l1*l1*q2*q2*r0 + 
  2*m2*m2*p0*p1*r1 + 4*m1*m2*p0*p2*r1 + 4*m0*m2*p1*p2*r1 + 
  2*m0*m1*p2*p2*r1 - 2*l2*m2*p1*q0*r1 - 2*l2*m1*p2*q0*r1 - 
  2*l1*m2*p2*q0*r1 - 2*l2*m2*p0*q1*r1 - 2*l2*m0*p2*q1*r1 - 
  2*l0*m2*p2*q1*r1 + 2*l2*l2*q0*q1*r1 - 2*l2*m1*p0*q2*r1 - 
  2*l1*m2*p0*q2*r1 - 2*l2*m0*p1*q2*r1 - 2*l0*m2*p1*q2*r1 - 
  2*l1*m0*p2*q2*r1 - 2*l0*m1*p2*q2*r1 + 4*l1*l2*q0*q2*r1 + 
  4*l0*l2*q1*q2*r1 + 2*l0*l1*q2*q2*r1 + 4*m1*m2*p0*p1*r2 + 
  2*m0*m2*p1*p1*r2 + 2*m1*m1*p0*p2*r2 + 4*m0*m1*p1*p2*r2 - 
  2*l2*m1*p1*q0*r2 - 2*l1*m2*p1*q0*r2 - 2*l1*m1*p2*q0*r2 - 
  2*l2*m1*p0*q1*r2 - 2*l1*m2*p0*q1*r2 - 2*l2*m0*p1*q1*r2 - 
  2*l0*m2*p1*q1*r2 - 2*l1*m0*p2*q1*r2 - 2*l0*m1*p2*q1*r2 + 
  4*l1*l2*q0*q1*r2 + 2*l0*l2*q1*q1*r2 - 2*l1*m1*p0*q2*r2 - 
  2*l1*m0*p1*q2*r2 - 2*l0*m1*p1*q2*r2 + 2*l1*l1*q0*q2*r2 + 
  4*l0*l1*q1*q2*r2 - 2*m2*m2*o1*p1*s0 - 4*m1*m2*o2*p1*s0 - 
  4*m1*m2*o1*p2*s0 - 2*m1*m1*o2*p2*s0 + 2*l2*m2*o1*q1*s0 + 
  2*l2*m1*o2*q1*s0 + 2*l1*m2*o2*q1*s0 + 2*k2*m2*p1*q1*s0 + 
  2*k2*m1*p2*q1*s0 + 2*k1*m2*p2*q1*s0 - 2*k2*l2*q1*q1*s0 + 
  2*l2*m1*o1*q2*s0 + 2*l1*m2*o1*q2*s0 + 2*l1*m1*o2*q2*s0 + 
  2*k2*m1*p1*q2*s0 + 2*k1*m2*p1*q2*s0 + 2*k1*m1*p2*q2*s0 - 
  4*k2*l1*q1*q2*s0 - 4*k1*l2*q1*q2*s0 - 2*k1*l1*q2*q2*s0 - 
  2*m2*m2*o1*p0*s1 - 4*m1*m2*o2*p0*s1 - 2*m2*m2*o0*p1*s1 - 
  4*m0*m2*o2*p1*s1 - 4*m1*m2*o0*p2*s1 - 4*m0*m2*o1*p2*s1 - 
  4*m0*m1*o2*p2*s1 + 2*l2*m2*o1*q0*s1 + 2*l2*m1*o2*q0*s1 + 
  2*l1*m2*o2*q0*s1 + 2*k2*m2*p1*q0*s1 + 2*k2*m1*p2*q0*s1 + 
  2*k1*m2*p2*q0*s1 + 2*l2*m2*o0*q1*s1 + 2*l2*m0*o2*q1*s1 + 
  2*l0*m2*o2*q1*s1 + 2*k2*m2*p0*q1*s1 + 2*k2*m0*p2*q1*s1 + 
  2*k0*m2*p2*q1*s1 - 4*k2*l2*q0*q1*s1 + 2*l2*m1*o0*q2*s1 + 
  2*l1*m2*o0*q2*s1 + 2*l2*m0*o1*q2*s1 + 2*l0*m2*o1*q2*s1 + 
  2*l1*m0*o2*q2*s1 + 2*l0*m1*o2*q2*s1 + 2*k2*m1*p0*q2*s1 + 
  2*k1*m2*p0*q2*s1 + 2*k2*m0*p1*q2*s1 + 2*k0*m2*p1*q2*s1 + 
  2*k1*m0*p2*q2*s1 + 2*k0*m1*p2*q2*s1 - 4*k2*l1*q0*q2*s1 - 
  4*k1*l2*q0*q2*s1 - 4*k2*l0*q1*q2*s1 - 4*k0*l2*q1*q2*s1 - 
  2*k1*l0*q2*q2*s1 - 2*k0*l1*q2*q2*s1 + 2*m2*m2*n1*s0*s1 + 
  4*m1*m2*n2*s0*s1 - 4*j2*m2*q1*s0*s1 - 4*j2*m1*q2*s0*s1 - 
  4*j1*m2*q2*s0*s1 + 4*i2*q1*q2*s0*s1 + 2*i1*q2*q2*s0*s1 + 
  m2*m2*n0*s1*s1 + 2*m0*m2*n2*s1*s1 - 2*j2*m2*q0*s1*s1 - 
  2*j2*m0*q2*s1*s1 - 2*j0*m2*q2*s1*s1 + 2*i2*q0*q2*s1*s1 + 
  i0*q2*q2*s1*s1 - 4*m1*m2*o1*p0*s2 - 2*m1*m1*o2*p0*s2 - 
  4*m1*m2*o0*p1*s2 - 4*m0*m2*o1*p1*s2 - 4*m0*m1*o2*p1*s2 - 
  2*m1*m1*o0*p2*s2 - 4*m0*m1*o1*p2*s2 + 2*l2*m1*o1*q0*s2 + 
  2*l1*m2*o1*q0*s2 + 2*l1*m1*o2*q0*s2 + 2*k2*m1*p1*q0*s2 + 
  2*k1*m2*p1*q0*s2 + 2*k1*m1*p2*q0*s2 + 2*l2*m1*o0*q1*s2 + 
  2*l1*m2*o0*q1*s2 + 2*l2*m0*o1*q1*s2 + 2*l0*m2*o1*q1*s2 + 
  2*l1*m0*o2*q1*s2 + 2*l0*m1*o2*q1*s2 + 2*k2*m1*p0*q1*s2 + 
  2*k1*m2*p0*q1*s2 + 2*k2*m0*p1*q1*s2 + 2*k0*m2*p1*q1*s2 + 
  2*k1*m0*p2*q1*s2 + 2*k0*m1*p2*q1*s2 - 4*k2*l1*q0*q1*s2 - 
  4*k1*l2*q0*q1*s2 - 2*k2*l0*q1*q1*s2 - 2*k0*l2*q1*q1*s2 + 
  2*l1*m1*o0*q2*s2 + 2*l1*m0*o1*q2*s2 + 2*l0*m1*o1*q2*s2 + 
  2*k1*m1*p0*q2*s2 + 2*k1*m0*p1*q2*s2 + 2*k0*m1*p1*q2*s2 - 
  4*k1*l1*q0*q2*s2 - 4*k1*l0*q1*q2*s2 - 4*k0*l1*q1*q2*s2 + 
  4*m1*m2*n1*s0*s2 + 2*m1*m1*n2*s0*s2 - 4*j2*m1*q1*s0*s2 - 
  4*j1*m2*q1*s0*s2 + 2*i2*q1*q1*s0*s2 - 4*j1*m1*q2*s0*s2 + 
  4*i1*q1*q2*s0*s2 + 4*m1*m2*n0*s1*s2 + 4*m0*m2*n1*s1*s2 + 
  4*m0*m1*n2*s1*s2 - 4*j2*m1*q0*s1*s2 - 4*j1*m2*q0*s1*s2 - 
  4*j2*m0*q1*s1*s2 - 4*j0*m2*q1*s1*s2 + 4*i2*q0*q1*s1*s2 - 
  4*j1*m0*q2*s1*s2 - 4*j0*m1*q2*s1*s2 + 4*i1*q0*q2*s1*s2 + 
  4*i0*q1*q2*s1*s2 + m1*m1*n0*s2*s2 + 2*m0*m1*n1*s2*s2 - 
  2*j1*m1*q0*s2*s2 - 2*j1*m0*q1*s2*s2 - 2*j0*m1*q1*s2*s2 + 
  2*i1*q0*q1*s2*s2 + i0*q1*q1*s2*s2 + 2*l2*m2*o1*p1*t0 + 
  2*l2*m1*o2*p1*t0 + 2*l1*m2*o2*p1*t0 - 2*k2*m2*p1*p1*t0 + 
  2*l2*m1*o1*p2*t0 + 2*l1*m2*o1*p2*t0 + 2*l1*m1*o2*p2*t0 - 
  4*k2*m1*p1*p2*t0 - 4*k1*m2*p1*p2*t0 - 2*k1*m1*p2*p2*t0 - 
  2*l2*l2*o1*q1*t0 - 4*l1*l2*o2*q1*t0 + 2*k2*l2*p1*q1*t0 + 
  2*k2*l1*p2*q1*t0 + 2*k1*l2*p2*q1*t0 - 4*l1*l2*o1*q2*t0 - 
  2*l1*l1*o2*q2*t0 + 2*k2*l1*p1*q2*t0 + 2*k1*l2*p1*q2*t0 + 
  2*k1*l1*p2*q2*t0 - 2*l2*m2*n1*s1*t0 - 2*l2*m1*n2*s1*t0 - 
  2*l1*m2*n2*s1*t0 + 2*j2*m2*p1*s1*t0 + 2*j2*m1*p2*s1*t0 + 
  2*j1*m2*p2*s1*t0 + 2*j2*l2*q1*s1*t0 - 2*i2*p2*q1*s1*t0 + 
  2*j2*l1*q2*s1*t0 + 2*j1*l2*q2*s1*t0 - 2*i2*p1*q2*s1*t0 - 
  2*i1*p2*q2*s1*t0 - 2*l2*m1*n1*s2*t0 - 2*l1*m2*n1*s2*t0 - 
  2*l1*m1*n2*s2*t0 + 2*j2*m1*p1*s2*t0 + 2*j1*m2*p1*s2*t0 + 
  2*j1*m1*p2*s2*t0 + 2*j2*l1*q1*s2*t0 + 2*j1*l2*q1*s2*t0 - 
  2*i2*p1*q1*s2*t0 - 2*i1*p2*q1*s2*t0 + 2*j1*l1*q2*s2*t0 - 
  2*i1*p1*q2*s2*t0 + 2*l2*m2*o1*p0*t1 + 2*l2*m1*o2*p0*t1 + 
  2*l1*m2*o2*p0*t1 + 2*l2*m2*o0*p1*t1 + 2*l2*m0*o2*p1*t1 + 
  2*l0*m2*o2*p1*t1 - 4*k2*m2*p0*p1*t1 + 2*l2*m1*o0*p2*t1 + 
  2*l1*m2*o0*p2*t1 + 2*l2*m0*o1*p2*t1 + 2*l0*m2*o1*p2*t1 + 
  2*l1*m0*o2*p2*t1 + 2*l0*m1*o2*p2*t1 - 4*k2*m1*p0*p2*t1 - 
  4*k1*m2*p0*p2*t1 - 4*k2*m0*p1*p2*t1 - 4*k0*m2*p1*p2*t1 - 
  2*k1*m0*p2*p2*t1 - 2*k0*m1*p2*p2*t1 - 2*l2*l2*o1*q0*t1 - 
  4*l1*l2*o2*q0*t1 + 2*k2*l2*p1*q0*t1 + 2*k2*l1*p2*q0*t1 + 
  2*k1*l2*p2*q0*t1 - 2*l2*l2*o0*q1*t1 - 4*l0*l2*o2*q1*t1 + 
  2*k2*l2*p0*q1*t1 + 2*k2*l0*p2*q1*t1 + 2*k0*l2*p2*q1*t1 - 
  4*l1*l2*o0*q2*t1 - 4*l0*l2*o1*q2*t1 - 4*l0*l1*o2*q2*t1 + 
  2*k2*l1*p0*q2*t1 + 2*k1*l2*p0*q2*t1 + 2*k2*l0*p1*q2*t1 + 
  2*k0*l2*p1*q2*t1 + 2*k1*l0*p2*q2*t1 + 2*k0*l1*p2*q2*t1 - 
  2*l2*m2*n1*s0*t1 - 2*l2*m1*n2*s0*t1 - 2*l1*m2*n2*s0*t1 + 
  2*j2*m2*p1*s0*t1 + 2*j2*m1*p2*s0*t1 + 2*j1*m2*p2*s0*t1 + 
  2*j2*l2*q1*s0*t1 - 2*i2*p2*q1*s0*t1 + 2*j2*l1*q2*s0*t1 + 
  2*j1*l2*q2*s0*t1 - 2*i2*p1*q2*s0*t1 - 2*i1*p2*q2*s0*t1 - 
  2*l2*m2*n0*s1*t1 - 2*l2*m0*n2*s1*t1 - 2*l0*m2*n2*s1*t1 + 
  2*j2*m2*p0*s1*t1 + 2*j2*m0*p2*s1*t1 + 2*j0*m2*p2*s1*t1 + 
  2*j2*l2*q0*s1*t1 - 2*i2*p2*q0*s1*t1 + 2*j2*l0*q2*s1*t1 + 
  2*j0*l2*q2*s1*t1 - 2*i2*p0*q2*s1*t1 - 2*i0*p2*q2*s1*t1 - 
  2*l2*m1*n0*s2*t1 - 2*l1*m2*n0*s2*t1 - 2*l2*m0*n1*s2*t1 - 
  2*l0*m2*n1*s2*t1 - 2*l1*m0*n2*s2*t1 - 2*l0*m1*n2*s2*t1 + 
  2*j2*m1*p0*s2*t1 + 2*j1*m2*p0*s2*t1 + 2*j2*m0*p1*s2*t1 + 
  2*j0*m2*p1*s2*t1 + 2*j1*m0*p2*s2*t1 + 2*j0*m1*p2*s2*t1 + 
  2*j2*l1*q0*s2*t1 + 2*j1*l2*q0*s2*t1 - 2*i2*p1*q0*s2*t1 - 
  2*i1*p2*q0*s2*t1 + 2*j2*l0*q1*s2*t1 + 2*j0*l2*q1*s2*t1 - 
  2*i2*p0*q1*s2*t1 - 2*i0*p2*q1*s2*t1 + 2*j1*l0*q2*s2*t1 + 
  2*j0*l1*q2*s2*t1 - 2*i1*p0*q2*s2*t1 - 2*i0*p1*q2*s2*t1 + 
  2*l2*l2*n1*t0*t1 + 4*l1*l2*n2*t0*t1 - 4*j2*l2*p1*t0*t1 - 
  4*j2*l1*p2*t0*t1 - 4*j1*l2*p2*t0*t1 + 4*i2*p1*p2*t0*t1 + 
  2*i1*p2*p2*t0*t1 + l2*l2*n0*t1*t1 + 2*l0*l2*n2*t1*t1 - 
  2*j2*l2*p0*t1*t1 - 2*j2*l0*p2*t1*t1 - 2*j0*l2*p2*t1*t1 + 
  2*i2*p0*p2*t1*t1 + i0*p2*p2*t1*t1 + 2*l2*m1*o1*p0*t2 + 
  2*l1*m2*o1*p0*t2 + 2*l1*m1*o2*p0*t2 + 2*l2*m1*o0*p1*t2 + 
  2*l1*m2*o0*p1*t2 + 2*l2*m0*o1*p1*t2 + 2*l0*m2*o1*p1*t2 + 
  2*l1*m0*o2*p1*t2 + 2*l0*m1*o2*p1*t2 - 4*k2*m1*p0*p1*t2 - 
  4*k1*m2*p0*p1*t2 - 2*k2*m0*p1*p1*t2 - 2*k0*m2*p1*p1*t2 + 
  2*l1*m1*o0*p2*t2 + 2*l1*m0*o1*p2*t2 + 2*l0*m1*o1*p2*t2 - 
  4*k1*m1*p0*p2*t2 - 4*k1*m0*p1*p2*t2 - 4*k0*m1*p1*p2*t2 - 
  4*l1*l2*o1*q0*t2 - 2*l1*l1*o2*q0*t2 + 2*k2*l1*p1*q0*t2 + 
  2*k1*l2*p1*q0*t2 + 2*k1*l1*p2*q0*t2 - 4*l1*l2*o0*q1*t2 - 
  4*l0*l2*o1*q1*t2 - 4*l0*l1*o2*q1*t2 + 2*k2*l1*p0*q1*t2 + 
  2*k1*l2*p0*q1*t2 + 2*k2*l0*p1*q1*t2 + 2*k0*l2*p1*q1*t2 + 
  2*k1*l0*p2*q1*t2 + 2*k0*l1*p2*q1*t2 - 2*l1*l1*o0*q2*t2 - 
  4*l0*l1*o1*q2*t2 + 2*k1*l1*p0*q2*t2 + 2*k1*l0*p1*q2*t2 + 
  2*k0*l1*p1*q2*t2 - 2*l2*m1*n1*s0*t2 - 2*l1*m2*n1*s0*t2 - 
  2*l1*m1*n2*s0*t2 + 2*j2*m1*p1*s0*t2 + 2*j1*m2*p1*s0*t2 + 
  2*j1*m1*p2*s0*t2 + 2*j2*l1*q1*s0*t2 + 2*j1*l2*q1*s0*t2 - 
  2*i2*p1*q1*s0*t2 - 2*i1*p2*q1*s0*t2 + 2*j1*l1*q2*s0*t2 - 
  2*i1*p1*q2*s0*t2 - 2*l2*m1*n0*s1*t2 - 2*l1*m2*n0*s1*t2 - 
  2*l2*m0*n1*s1*t2 - 2*l0*m2*n1*s1*t2 - 2*l1*m0*n2*s1*t2 - 
  2*l0*m1*n2*s1*t2 + 2*j2*m1*p0*s1*t2 + 2*j1*m2*p0*s1*t2 + 
  2*j2*m0*p1*s1*t2 + 2*j0*m2*p1*s1*t2 + 2*j1*m0*p2*s1*t2 + 
  2*j0*m1*p2*s1*t2 + 2*j2*l1*q0*s1*t2 + 2*j1*l2*q0*s1*t2 - 
  2*i2*p1*q0*s1*t2 - 2*i1*p2*q0*s1*t2 + 2*j2*l0*q1*s1*t2 + 
  2*j0*l2*q1*s1*t2 - 2*i2*p0*q1*s1*t2 - 2*i0*p2*q1*s1*t2 + 
  2*j1*l0*q2*s1*t2 + 2*j0*l1*q2*s1*t2 - 2*i1*p0*q2*s1*t2 - 
  2*i0*p1*q2*s1*t2 - 2*l1*m1*n0*s2*t2 - 2*l1*m0*n1*s2*t2 - 
  2*l0*m1*n1*s2*t2 + 2*j1*m1*p0*s2*t2 + 2*j1*m0*p1*s2*t2 + 
  2*j0*m1*p1*s2*t2 + 2*j1*l1*q0*s2*t2 - 2*i1*p1*q0*s2*t2 + 
  2*j1*l0*q1*s2*t2 + 2*j0*l1*q1*s2*t2 - 2*i1*p0*q1*s2*t2 - 
  2*i0*p1*q1*s2*t2 + 4*l1*l2*n1*t0*t2 + 2*l1*l1*n2*t0*t2 - 
  4*j2*l1*p1*t0*t2 - 4*j1*l2*p1*t0*t2 + 2*i2*p1*p1*t0*t2 - 
  4*j1*l1*p2*t0*t2 + 4*i1*p1*p2*t0*t2 + 4*l1*l2*n0*t1*t2 + 
  4*l0*l2*n1*t1*t2 + 4*l0*l1*n2*t1*t2 - 4*j2*l1*p0*t1*t2 - 
  4*j1*l2*p0*t1*t2 - 4*j2*l0*p1*t1*t2 - 4*j0*l2*p1*t1*t2 + 
  4*i2*p0*p1*t1*t2 - 4*j1*l0*p2*t1*t2 - 4*j0*l1*p2*t1*t2 + 
  4*i1*p0*p2*t1*t2 + 4*i0*p1*p2*t1*t2 + l1*l1*n0*t2*t2 + 
  2*l0*l1*n1*t2*t2 - 2*j1*l1*p0*t2*t2 - 2*j1*l0*p1*t2*t2 - 
  2*j0*l1*p1*t2*t2 + 2*i1*p0*p1*t2*t2 + i0*p1*p1*t2*t2 + 
  m2*m2*o1*o1*u0 + 4*m1*m2*o1*o2*u0 + m1*m1*o2*o2*u0 - 
  2*k2*m2*o1*q1*u0 - 2*k2*m1*o2*q1*u0 - 2*k1*m2*o2*q1*u0 + 
  k2*k2*q1*q1*u0 - 2*k2*m1*o1*q2*u0 - 2*k1*m2*o1*q2*u0 - 
  2*k1*m1*o2*q2*u0 + 4*k1*k2*q1*q2*u0 + k1*k1*q2*q2*u0 - 
  m2*m2*n1*r1*u0 - 2*m1*m2*n2*r1*u0 + 2*j2*m2*q1*r1*u0 + 
  2*j2*m1*q2*r1*u0 + 2*j1*m2*q2*r1*u0 - 2*i2*q1*q2*r1*u0 - 
  i1*q2*q2*r1*u0 - 2*m1*m2*n1*r2*u0 - m1*m1*n2*r2*u0 + 
  2*j2*m1*q1*r2*u0 + 2*j1*m2*q1*r2*u0 - i2*q1*q1*r2*u0 + 
  2*j1*m1*q2*r2*u0 - 2*i1*q1*q2*r2*u0 + 2*k2*m2*n1*t1*u0 + 
  2*k2*m1*n2*t1*u0 + 2*k1*m2*n2*t1*u0 - 2*j2*m2*o1*t1*u0 - 
  2*j2*m1*o2*t1*u0 - 2*j1*m2*o2*t1*u0 - 2*j2*k2*q1*t1*u0 + 
  2*i2*o2*q1*t1*u0 - 2*j2*k1*q2*t1*u0 - 2*j1*k2*q2*t1*u0 + 
  2*i2*o1*q2*t1*u0 + 2*i1*o2*q2*t1*u0 + j2*j2*t1*t1*u0 - 
  i2*n2*t1*t1*u0 + 2*k2*m1*n1*t2*u0 + 2*k1*m2*n1*t2*u0 + 
  2*k1*m1*n2*t2*u0 - 2*j2*m1*o1*t2*u0 - 2*j1*m2*o1*t2*u0 - 
  2*j1*m1*o2*t2*u0 - 2*j2*k1*q1*t2*u0 - 2*j1*k2*q1*t2*u0 + 
  2*i2*o1*q1*t2*u0 + 2*i1*o2*q1*t2*u0 - 2*j1*k1*q2*t2*u0 + 
  2*i1*o1*q2*t2*u0 + 4*j1*j2*t1*t2*u0 - 2*i2*n1*t1*t2*u0 - 
  2*i1*n2*t1*t2*u0 + j1*j1*t2*t2*u0 - i1*n1*t2*t2*u0 + 
  2*m2*m2*o0*o1*u1 + 4*m1*m2*o0*o2*u1 + 4*m0*m2*o1*o2*u1 + 
  2*m0*m1*o2*o2*u1 - 2*k2*m2*o1*q0*u1 - 2*k2*m1*o2*q0*u1 - 
  2*k1*m2*o2*q0*u1 - 2*k2*m2*o0*q1*u1 - 2*k2*m0*o2*q1*u1 - 
  2*k0*m2*o2*q1*u1 + 2*k2*k2*q0*q1*u1 - 2*k2*m1*o0*q2*u1 - 
  2*k1*m2*o0*q2*u1 - 2*k2*m0*o1*q2*u1 - 2*k0*m2*o1*q2*u1 - 
  2*k1*m0*o2*q2*u1 - 2*k0*m1*o2*q2*u1 + 4*k1*k2*q0*q2*u1 + 
  4*k0*k2*q1*q2*u1 + 2*k0*k1*q2*q2*u1 - m2*m2*n1*r0*u1 - 
  2*m1*m2*n2*r0*u1 + 2*j2*m2*q1*r0*u1 + 2*j2*m1*q2*r0*u1 + 
  2*j1*m2*q2*r0*u1 - 2*i2*q1*q2*r0*u1 - i1*q2*q2*r0*u1 - 
  m2*m2*n0*r1*u1 - 2*m0*m2*n2*r1*u1 + 2*j2*m2*q0*r1*u1 + 
  2*j2*m0*q2*r1*u1 + 2*j0*m2*q2*r1*u1 - 2*i2*q0*q2*r1*u1 - 
  i0*q2*q2*r1*u1 - 2*m1*m2*n0*r2*u1 - 2*m0*m2*n1*r2*u1 - 
  2*m0*m1*n2*r2*u1 + 2*j2*m1*q0*r2*u1 + 2*j1*m2*q0*r2*u1 + 
  2*j2*m0*q1*r2*u1 + 2*j0*m2*q1*r2*u1 - 2*i2*q0*q1*r2*u1 + 
  2*j1*m0*q2*r2*u1 + 2*j0*m1*q2*r2*u1 - 2*i1*q0*q2*r2*u1 - 
  2*i0*q1*q2*r2*u1 + 2*k2*m2*n1*t0*u1 + 2*k2*m1*n2*t0*u1 + 
  2*k1*m2*n2*t0*u1 - 2*j2*m2*o1*t0*u1 - 2*j2*m1*o2*t0*u1 - 
  2*j1*m2*o2*t0*u1 - 2*j2*k2*q1*t0*u1 + 2*i2*o2*q1*t0*u1 - 
  2*j2*k1*q2*t0*u1 - 2*j1*k2*q2*t0*u1 + 2*i2*o1*q2*t0*u1 + 
  2*i1*o2*q2*t0*u1 + 2*k2*m2*n0*t1*u1 + 2*k2*m0*n2*t1*u1 + 
  2*k0*m2*n2*t1*u1 - 2*j2*m2*o0*t1*u1 - 2*j2*m0*o2*t1*u1 - 
  2*j0*m2*o2*t1*u1 - 2*j2*k2*q0*t1*u1 + 2*i2*o2*q0*t1*u1 - 
  2*j2*k0*q2*t1*u1 - 2*j0*k2*q2*t1*u1 + 2*i2*o0*q2*t1*u1 + 
  2*i0*o2*q2*t1*u1 + 2*j2*j2*t0*t1*u1 - 2*i2*n2*t0*t1*u1 + 
  2*k2*m1*n0*t2*u1 + 2*k1*m2*n0*t2*u1 + 2*k2*m0*n1*t2*u1 + 
  2*k0*m2*n1*t2*u1 + 2*k1*m0*n2*t2*u1 + 2*k0*m1*n2*t2*u1 - 
  2*j2*m1*o0*t2*u1 - 2*j1*m2*o0*t2*u1 - 2*j2*m0*o1*t2*u1 - 
  2*j0*m2*o1*t2*u1 - 2*j1*m0*o2*t2*u1 - 2*j0*m1*o2*t2*u1 - 
  2*j2*k1*q0*t2*u1 - 2*j1*k2*q0*t2*u1 + 2*i2*o1*q0*t2*u1 + 
  2*i1*o2*q0*t2*u1 - 2*j2*k0*q1*t2*u1 - 2*j0*k2*q1*t2*u1 + 
  2*i2*o0*q1*t2*u1 + 2*i0*o2*q1*t2*u1 - 2*j1*k0*q2*t2*u1 - 
  2*j0*k1*q2*t2*u1 + 2*i1*o0*q2*t2*u1 + 2*i0*o1*q2*t2*u1 + 
  4*j1*j2*t0*t2*u1 - 2*i2*n1*t0*t2*u1 - 2*i1*n2*t0*t2*u1 + 
  4*j0*j2*t1*t2*u1 - 2*i2*n0*t1*t2*u1 - 2*i0*n2*t1*t2*u1 + 
  2*j0*j1*t2*t2*u1 - i1*n0*t2*t2*u1 - i0*n1*t2*t2*u1 + 
  4*m1*m2*o0*o1*u2 + 2*m0*m2*o1*o1*u2 + 2*m1*m1*o0*o2*u2 + 
  4*m0*m1*o1*o2*u2 - 2*k2*m1*o1*q0*u2 - 2*k1*m2*o1*q0*u2 - 
  2*k1*m1*o2*q0*u2 - 2*k2*m1*o0*q1*u2 - 2*k1*m2*o0*q1*u2 - 
  2*k2*m0*o1*q1*u2 - 2*k0*m2*o1*q1*u2 - 2*k1*m0*o2*q1*u2 - 
  2*k0*m1*o2*q1*u2 + 4*k1*k2*q0*q1*u2 + 2*k0*k2*q1*q1*u2 - 
  2*k1*m1*o0*q2*u2 - 2*k1*m0*o1*q2*u2 - 2*k0*m1*o1*q2*u2 + 
  2*k1*k1*q0*q2*u2 + 4*k0*k1*q1*q2*u2 - 2*m1*m2*n1*r0*u2 - 
  m1*m1*n2*r0*u2 + 2*j2*m1*q1*r0*u2 + 2*j1*m2*q1*r0*u2 - 
  i2*q1*q1*r0*u2 + 2*j1*m1*q2*r0*u2 - 2*i1*q1*q2*r0*u2 - 
  2*m1*m2*n0*r1*u2 - 2*m0*m2*n1*r1*u2 - 2*m0*m1*n2*r1*u2 + 
  2*j2*m1*q0*r1*u2 + 2*j1*m2*q0*r1*u2 + 2*j2*m0*q1*r1*u2 + 
  2*j0*m2*q1*r1*u2 - 2*i2*q0*q1*r1*u2 + 2*j1*m0*q2*r1*u2 + 
  2*j0*m1*q2*r1*u2 - 2*i1*q0*q2*r1*u2 - 2*i0*q1*q2*r1*u2 - 
  m1*m1*n0*r2*u2 - 2*m0*m1*n1*r2*u2 + 2*j1*m1*q0*r2*u2 + 
  2*j1*m0*q1*r2*u2 + 2*j0*m1*q1*r2*u2 - 2*i1*q0*q1*r2*u2 - 
  i0*q1*q1*r2*u2 + 2*k2*m1*n1*t0*u2 + 2*k1*m2*n1*t0*u2 + 
  2*k1*m1*n2*t0*u2 - 2*j2*m1*o1*t0*u2 - 2*j1*m2*o1*t0*u2 - 
  2*j1*m1*o2*t0*u2 - 2*j2*k1*q1*t0*u2 - 2*j1*k2*q1*t0*u2 + 
  2*i2*o1*q1*t0*u2 + 2*i1*o2*q1*t0*u2 - 2*j1*k1*q2*t0*u2 + 
  2*i1*o1*q2*t0*u2 + 2*k2*m1*n0*t1*u2 + 2*k1*m2*n0*t1*u2 + 
  2*k2*m0*n1*t1*u2 + 2*k0*m2*n1*t1*u2 + 2*k1*m0*n2*t1*u2 + 
  2*k0*m1*n2*t1*u2 - 2*j2*m1*o0*t1*u2 - 2*j1*m2*o0*t1*u2 - 
  2*j2*m0*o1*t1*u2 - 2*j0*m2*o1*t1*u2 - 2*j1*m0*o2*t1*u2 - 
  2*j0*m1*o2*t1*u2 - 2*j2*k1*q0*t1*u2 - 2*j1*k2*q0*t1*u2 + 
  2*i2*o1*q0*t1*u2 + 2*i1*o2*q0*t1*u2 - 2*j2*k0*q1*t1*u2 - 
  2*j0*k2*q1*t1*u2 + 2*i2*o0*q1*t1*u2 + 2*i0*o2*q1*t1*u2 - 
  2*j1*k0*q2*t1*u2 - 2*j0*k1*q2*t1*u2 + 2*i1*o0*q2*t1*u2 + 
  2*i0*o1*q2*t1*u2 + 4*j1*j2*t0*t1*u2 - 2*i2*n1*t0*t1*u2 - 
  2*i1*n2*t0*t1*u2 + 2*j0*j2*t1*t1*u2 - i2*n0*t1*t1*u2 - 
  i0*n2*t1*t1*u2 + 2*k1*m1*n0*t2*u2 + 2*k1*m0*n1*t2*u2 + 
  2*k0*m1*n1*t2*u2 - 2*j1*m1*o0*t2*u2 - 2*j1*m0*o1*t2*u2 - 
  2*j0*m1*o1*t2*u2 - 2*j1*k1*q0*t2*u2 + 2*i1*o1*q0*t2*u2 - 
  2*j1*k0*q1*t2*u2 - 2*j0*k1*q1*t2*u2 + 2*i1*o0*q1*t2*u2 + 
  2*i0*o1*q1*t2*u2 + 2*j1*j1*t0*t2*u2 - 2*i1*n1*t0*t2*u2 + 
  4*j0*j1*t1*t2*u2 - 2*i1*n0*t1*t2*u2 - 2*i0*n1*t1*t2*u2 - 
  2*l2*m2*o1*o1*v0 - 4*l2*m1*o1*o2*v0 - 4*l1*m2*o1*o2*v0 - 
  2*l1*m1*o2*o2*v0 + 2*k2*m2*o1*p1*v0 + 2*k2*m1*o2*p1*v0 + 
  2*k1*m2*o2*p1*v0 + 2*k2*m1*o1*p2*v0 + 2*k1*m2*o1*p2*v0 + 
  2*k1*m1*o2*p2*v0 + 2*k2*l2*o1*q1*v0 + 2*k2*l1*o2*q1*v0 + 
  2*k1*l2*o2*q1*v0 - 2*k2*k2*p1*q1*v0 - 4*k1*k2*p2*q1*v0 + 
  2*k2*l1*o1*q2*v0 + 2*k1*l2*o1*q2*v0 + 2*k1*l1*o2*q2*v0 - 
  4*k1*k2*p1*q2*v0 - 2*k1*k1*p2*q2*v0 + 2*l2*m2*n1*r1*v0 + 
  2*l2*m1*n2*r1*v0 + 2*l1*m2*n2*r1*v0 - 2*j2*m2*p1*r1*v0 - 
  2*j2*m1*p2*r1*v0 - 2*j1*m2*p2*r1*v0 - 2*j2*l2*q1*r1*v0 + 
  2*i2*p2*q1*r1*v0 - 2*j2*l1*q2*r1*v0 - 2*j1*l2*q2*r1*v0 + 
  2*i2*p1*q2*r1*v0 + 2*i1*p2*q2*r1*v0 + 2*l2*m1*n1*r2*v0 + 
  2*l1*m2*n1*r2*v0 + 2*l1*m1*n2*r2*v0 - 2*j2*m1*p1*r2*v0 - 
  2*j1*m2*p1*r2*v0 - 2*j1*m1*p2*r2*v0 - 2*j2*l1*q1*r2*v0 - 
  2*j1*l2*q1*r2*v0 + 2*i2*p1*q1*r2*v0 + 2*i1*p2*q1*r2*v0 - 
  2*j1*l1*q2*r2*v0 + 2*i1*p1*q2*r2*v0 - 2*k2*m2*n1*s1*v0 - 
  2*k2*m1*n2*s1*v0 - 2*k1*m2*n2*s1*v0 + 2*j2*m2*o1*s1*v0 + 
  2*j2*m1*o2*s1*v0 + 2*j1*m2*o2*s1*v0 + 2*j2*k2*q1*s1*v0 - 
  2*i2*o2*q1*s1*v0 + 2*j2*k1*q2*s1*v0 + 2*j1*k2*q2*s1*v0 - 
  2*i2*o1*q2*s1*v0 - 2*i1*o2*q2*s1*v0 - 2*k2*m1*n1*s2*v0 - 
  2*k1*m2*n1*s2*v0 - 2*k1*m1*n2*s2*v0 + 2*j2*m1*o1*s2*v0 + 
  2*j1*m2*o1*s2*v0 + 2*j1*m1*o2*s2*v0 + 2*j2*k1*q1*s2*v0 + 
  2*j1*k2*q1*s2*v0 - 2*i2*o1*q1*s2*v0 - 2*i1*o2*q1*s2*v0 + 
  2*j1*k1*q2*s2*v0 - 2*i1*o1*q2*s2*v0 - 2*k2*l2*n1*t1*v0 - 
  2*k2*l1*n2*t1*v0 - 2*k1*l2*n2*t1*v0 + 2*j2*l2*o1*t1*v0 + 
  2*j2*l1*o2*t1*v0 + 2*j1*l2*o2*t1*v0 + 2*j2*k2*p1*t1*v0 - 
  2*i2*o2*p1*t1*v0 + 2*j2*k1*p2*t1*v0 + 2*j1*k2*p2*t1*v0 - 
  2*i2*o1*p2*t1*v0 - 2*i1*o2*p2*t1*v0 - 2*j2*j2*s1*t1*v0 + 
  2*i2*n2*s1*t1*v0 - 4*j1*j2*s2*t1*v0 + 2*i2*n1*s2*t1*v0 + 
  2*i1*n2*s2*t1*v0 - 2*k2*l1*n1*t2*v0 - 2*k1*l2*n1*t2*v0 - 
  2*k1*l1*n2*t2*v0 + 2*j2*l1*o1*t2*v0 + 2*j1*l2*o1*t2*v0 + 
  2*j1*l1*o2*t2*v0 + 2*j2*k1*p1*t2*v0 + 2*j1*k2*p1*t2*v0 - 
  2*i2*o1*p1*t2*v0 - 2*i1*o2*p1*t2*v0 + 2*j1*k1*p2*t2*v0 - 
  2*i1*o1*p2*t2*v0 - 4*j1*j2*s1*t2*v0 + 2*i2*n1*s1*t2*v0 + 
  2*i1*n2*s1*t2*v0 - 2*j1*j1*s2*t2*v0 + 2*i1*n1*s2*t2*v0 - 
  4*l2*m2*o0*o1*v1 - 4*l2*m1*o0*o2*v1 - 4*l1*m2*o0*o2*v1 - 
  4*l2*m0*o1*o2*v1 - 4*l0*m2*o1*o2*v1 - 2*l1*m0*o2*o2*v1 - 
  2*l0*m1*o2*o2*v1 + 2*k2*m2*o1*p0*v1 + 2*k2*m1*o2*p0*v1 + 
  2*k1*m2*o2*p0*v1 + 2*k2*m2*o0*p1*v1 + 2*k2*m0*o2*p1*v1 + 
  2*k0*m2*o2*p1*v1 + 2*k2*m1*o0*p2*v1 + 2*k1*m2*o0*p2*v1 + 
  2*k2*m0*o1*p2*v1 + 2*k0*m2*o1*p2*v1 + 2*k1*m0*o2*p2*v1 + 
  2*k0*m1*o2*p2*v1 + 2*k2*l2*o1*q0*v1 + 2*k2*l1*o2*q0*v1 + 
  2*k1*l2*o2*q0*v1 - 2*k2*k2*p1*q0*v1 - 4*k1*k2*p2*q0*v1 + 
  2*k2*l2*o0*q1*v1 + 2*k2*l0*o2*q1*v1 + 2*k0*l2*o2*q1*v1 - 
  2*k2*k2*p0*q1*v1 - 4*k0*k2*p2*q1*v1 + 2*k2*l1*o0*q2*v1 + 
  2*k1*l2*o0*q2*v1 + 2*k2*l0*o1*q2*v1 + 2*k0*l2*o1*q2*v1 + 
  2*k1*l0*o2*q2*v1 + 2*k0*l1*o2*q2*v1 - 4*k1*k2*p0*q2*v1 - 
  4*k0*k2*p1*q2*v1 - 4*k0*k1*p2*q2*v1 + 2*l2*m2*n1*r0*v1 + 
  2*l2*m1*n2*r0*v1 + 2*l1*m2*n2*r0*v1 - 2*j2*m2*p1*r0*v1 - 
  2*j2*m1*p2*r0*v1 - 2*j1*m2*p2*r0*v1 - 2*j2*l2*q1*r0*v1 + 
  2*i2*p2*q1*r0*v1 - 2*j2*l1*q2*r0*v1 - 2*j1*l2*q2*r0*v1 + 
  2*i2*p1*q2*r0*v1 + 2*i1*p2*q2*r0*v1 + 2*l2*m2*n0*r1*v1 + 
  2*l2*m0*n2*r1*v1 + 2*l0*m2*n2*r1*v1 - 2*j2*m2*p0*r1*v1 - 
  2*j2*m0*p2*r1*v1 - 2*j0*m2*p2*r1*v1 - 2*j2*l2*q0*r1*v1 + 
  2*i2*p2*q0*r1*v1 - 2*j2*l0*q2*r1*v1 - 2*j0*l2*q2*r1*v1 + 
  2*i2*p0*q2*r1*v1 + 2*i0*p2*q2*r1*v1 + 2*l2*m1*n0*r2*v1 + 
  2*l1*m2*n0*r2*v1 + 2*l2*m0*n1*r2*v1 + 2*l0*m2*n1*r2*v1 + 
  2*l1*m0*n2*r2*v1 + 2*l0*m1*n2*r2*v1 - 2*j2*m1*p0*r2*v1 - 
  2*j1*m2*p0*r2*v1 - 2*j2*m0*p1*r2*v1 - 2*j0*m2*p1*r2*v1 - 
  2*j1*m0*p2*r2*v1 - 2*j0*m1*p2*r2*v1 - 2*j2*l1*q0*r2*v1 - 
  2*j1*l2*q0*r2*v1 + 2*i2*p1*q0*r2*v1 + 2*i1*p2*q0*r2*v1 - 
  2*j2*l0*q1*r2*v1 - 2*j0*l2*q1*r2*v1 + 2*i2*p0*q1*r2*v1 + 
  2*i0*p2*q1*r2*v1 - 2*j1*l0*q2*r2*v1 - 2*j0*l1*q2*r2*v1 + 
  2*i1*p0*q2*r2*v1 + 2*i0*p1*q2*r2*v1 - 2*k2*m2*n1*s0*v1 - 
  2*k2*m1*n2*s0*v1 - 2*k1*m2*n2*s0*v1 + 2*j2*m2*o1*s0*v1 + 
  2*j2*m1*o2*s0*v1 + 2*j1*m2*o2*s0*v1 + 2*j2*k2*q1*s0*v1 - 
  2*i2*o2*q1*s0*v1 + 2*j2*k1*q2*s0*v1 + 2*j1*k2*q2*s0*v1 - 
  2*i2*o1*q2*s0*v1 - 2*i1*o2*q2*s0*v1 - 2*k2*m2*n0*s1*v1 - 
  2*k2*m0*n2*s1*v1 - 2*k0*m2*n2*s1*v1 + 2*j2*m2*o0*s1*v1 + 
  2*j2*m0*o2*s1*v1 + 2*j0*m2*o2*s1*v1 + 2*j2*k2*q0*s1*v1 - 
  2*i2*o2*q0*s1*v1 + 2*j2*k0*q2*s1*v1 + 2*j0*k2*q2*s1*v1 - 
  2*i2*o0*q2*s1*v1 - 2*i0*o2*q2*s1*v1 - 2*k2*m1*n0*s2*v1 - 
  2*k1*m2*n0*s2*v1 - 2*k2*m0*n1*s2*v1 - 2*k0*m2*n1*s2*v1 - 
  2*k1*m0*n2*s2*v1 - 2*k0*m1*n2*s2*v1 + 2*j2*m1*o0*s2*v1 + 
  2*j1*m2*o0*s2*v1 + 2*j2*m0*o1*s2*v1 + 2*j0*m2*o1*s2*v1 + 
  2*j1*m0*o2*s2*v1 + 2*j0*m1*o2*s2*v1 + 2*j2*k1*q0*s2*v1 + 
  2*j1*k2*q0*s2*v1 - 2*i2*o1*q0*s2*v1 - 2*i1*o2*q0*s2*v1 + 
  2*j2*k0*q1*s2*v1 + 2*j0*k2*q1*s2*v1 - 2*i2*o0*q1*s2*v1 - 
  2*i0*o2*q1*s2*v1 + 2*j1*k0*q2*s2*v1 + 2*j0*k1*q2*s2*v1 - 
  2*i1*o0*q2*s2*v1 - 2*i0*o1*q2*s2*v1 - 2*k2*l2*n1*t0*v1 - 
  2*k2*l1*n2*t0*v1 - 2*k1*l2*n2*t0*v1 + 2*j2*l2*o1*t0*v1 + 
  2*j2*l1*o2*t0*v1 + 2*j1*l2*o2*t0*v1 + 2*j2*k2*p1*t0*v1 - 
  2*i2*o2*p1*t0*v1 + 2*j2*k1*p2*t0*v1 + 2*j1*k2*p2*t0*v1 - 
  2*i2*o1*p2*t0*v1 - 2*i1*o2*p2*t0*v1 - 2*j2*j2*s1*t0*v1 + 
  2*i2*n2*s1*t0*v1 - 4*j1*j2*s2*t0*v1 + 2*i2*n1*s2*t0*v1 + 
  2*i1*n2*s2*t0*v1 - 2*k2*l2*n0*t1*v1 - 2*k2*l0*n2*t1*v1 - 
  2*k0*l2*n2*t1*v1 + 2*j2*l2*o0*t1*v1 + 2*j2*l0*o2*t1*v1 + 
  2*j0*l2*o2*t1*v1 + 2*j2*k2*p0*t1*v1 - 2*i2*o2*p0*t1*v1 + 
  2*j2*k0*p2*t1*v1 + 2*j0*k2*p2*t1*v1 - 2*i2*o0*p2*t1*v1 - 
  2*i0*o2*p2*t1*v1 - 2*j2*j2*s0*t1*v1 + 2*i2*n2*s0*t1*v1 - 
  4*j0*j2*s2*t1*v1 + 2*i2*n0*s2*t1*v1 + 2*i0*n2*s2*t1*v1 - 
  2*k2*l1*n0*t2*v1 - 2*k1*l2*n0*t2*v1 - 2*k2*l0*n1*t2*v1 - 
  2*k0*l2*n1*t2*v1 - 2*k1*l0*n2*t2*v1 - 2*k0*l1*n2*t2*v1 + 
  2*j2*l1*o0*t2*v1 + 2*j1*l2*o0*t2*v1 + 2*j2*l0*o1*t2*v1 + 
  2*j0*l2*o1*t2*v1 + 2*j1*l0*o2*t2*v1 + 2*j0*l1*o2*t2*v1 + 
  2*j2*k1*p0*t2*v1 + 2*j1*k2*p0*t2*v1 - 2*i2*o1*p0*t2*v1 - 
  2*i1*o2*p0*t2*v1 + 2*j2*k0*p1*t2*v1 + 2*j0*k2*p1*t2*v1 - 
  2*i2*o0*p1*t2*v1 - 2*i0*o2*p1*t2*v1 + 2*j1*k0*p2*t2*v1 + 
  2*j0*k1*p2*t2*v1 - 2*i1*o0*p2*t2*v1 - 2*i0*o1*p2*t2*v1 - 
  4*j1*j2*s0*t2*v1 + 2*i2*n1*s0*t2*v1 + 2*i1*n2*s0*t2*v1 - 
  4*j0*j2*s1*t2*v1 + 2*i2*n0*s1*t2*v1 + 2*i0*n2*s1*t2*v1 - 
  4*j0*j1*s2*t2*v1 + 2*i1*n0*s2*t2*v1 + 2*i0*n1*s2*t2*v1 + 
  2*k2*k2*n1*v0*v1 + 4*k1*k2*n2*v0*v1 - 4*j2*k2*o1*v0*v1 - 
  4*j2*k1*o2*v0*v1 - 4*j1*k2*o2*v0*v1 + 4*i2*o1*o2*v0*v1 + 
  2*i1*o2*o2*v0*v1 + 2*j2*j2*r1*v0*v1 - 2*i2*n2*r1*v0*v1 + 
  4*j1*j2*r2*v0*v1 - 2*i2*n1*r2*v0*v1 - 2*i1*n2*r2*v0*v1 + 
  k2*k2*n0*v1*v1 + 2*k0*k2*n2*v1*v1 - 2*j2*k2*o0*v1*v1 - 
  2*j2*k0*o2*v1*v1 - 2*j0*k2*o2*v1*v1 + 2*i2*o0*o2*v1*v1 + 
  i0*o2*o2*v1*v1 + j2*j2*r0*v1*v1 - i2*n2*r0*v1*v1 + 
  2*j0*j2*r2*v1*v1 - i2*n0*r2*v1*v1 - i0*n2*r2*v1*v1 - 
  4*l2*m1*o0*o1*v2 - 4*l1*m2*o0*o1*v2 - 2*l2*m0*o1*o1*v2 - 
  2*l0*m2*o1*o1*v2 - 4*l1*m1*o0*o2*v2 - 4*l1*m0*o1*o2*v2 - 
  4*l0*m1*o1*o2*v2 + 2*k2*m1*o1*p0*v2 + 2*k1*m2*o1*p0*v2 + 
  2*k1*m1*o2*p0*v2 + 2*k2*m1*o0*p1*v2 + 2*k1*m2*o0*p1*v2 + 
  2*k2*m0*o1*p1*v2 + 2*k0*m2*o1*p1*v2 + 2*k1*m0*o2*p1*v2 + 
  2*k0*m1*o2*p1*v2 + 2*k1*m1*o0*p2*v2 + 2*k1*m0*o1*p2*v2 + 
  2*k0*m1*o1*p2*v2 + 2*k2*l1*o1*q0*v2 + 2*k1*l2*o1*q0*v2 + 
  2*k1*l1*o2*q0*v2 - 4*k1*k2*p1*q0*v2 - 2*k1*k1*p2*q0*v2 + 
  2*k2*l1*o0*q1*v2 + 2*k1*l2*o0*q1*v2 + 2*k2*l0*o1*q1*v2 + 
  2*k0*l2*o1*q1*v2 + 2*k1*l0*o2*q1*v2 + 2*k0*l1*o2*q1*v2 - 
  4*k1*k2*p0*q1*v2 - 4*k0*k2*p1*q1*v2 - 4*k0*k1*p2*q1*v2 + 
  2*k1*l1*o0*q2*v2 + 2*k1*l0*o1*q2*v2 + 2*k0*l1*o1*q2*v2 - 
  2*k1*k1*p0*q2*v2 - 4*k0*k1*p1*q2*v2 + 2*l2*m1*n1*r0*v2 + 
  2*l1*m2*n1*r0*v2 + 2*l1*m1*n2*r0*v2 - 2*j2*m1*p1*r0*v2 - 
  2*j1*m2*p1*r0*v2 - 2*j1*m1*p2*r0*v2 - 2*j2*l1*q1*r0*v2 - 
  2*j1*l2*q1*r0*v2 + 2*i2*p1*q1*r0*v2 + 2*i1*p2*q1*r0*v2 - 
  2*j1*l1*q2*r0*v2 + 2*i1*p1*q2*r0*v2 + 2*l2*m1*n0*r1*v2 + 
  2*l1*m2*n0*r1*v2 + 2*l2*m0*n1*r1*v2 + 2*l0*m2*n1*r1*v2 + 
  2*l1*m0*n2*r1*v2 + 2*l0*m1*n2*r1*v2 - 2*j2*m1*p0*r1*v2 - 
  2*j1*m2*p0*r1*v2 - 2*j2*m0*p1*r1*v2 - 2*j0*m2*p1*r1*v2 - 
  2*j1*m0*p2*r1*v2 - 2*j0*m1*p2*r1*v2 - 2*j2*l1*q0*r1*v2 - 
  2*j1*l2*q0*r1*v2 + 2*i2*p1*q0*r1*v2 + 2*i1*p2*q0*r1*v2 - 
  2*j2*l0*q1*r1*v2 - 2*j0*l2*q1*r1*v2 + 2*i2*p0*q1*r1*v2 + 
  2*i0*p2*q1*r1*v2 - 2*j1*l0*q2*r1*v2 - 2*j0*l1*q2*r1*v2 + 
  2*i1*p0*q2*r1*v2 + 2*i0*p1*q2*r1*v2 + 2*l1*m1*n0*r2*v2 + 
  2*l1*m0*n1*r2*v2 + 2*l0*m1*n1*r2*v2 - 2*j1*m1*p0*r2*v2 - 
  2*j1*m0*p1*r2*v2 - 2*j0*m1*p1*r2*v2 - 2*j1*l1*q0*r2*v2 + 
  2*i1*p1*q0*r2*v2 - 2*j1*l0*q1*r2*v2 - 2*j0*l1*q1*r2*v2 + 
  2*i1*p0*q1*r2*v2 + 2*i0*p1*q1*r2*v2 - 2*k2*m1*n1*s0*v2 - 
  2*k1*m2*n1*s0*v2 - 2*k1*m1*n2*s0*v2 + 2*j2*m1*o1*s0*v2 + 
  2*j1*m2*o1*s0*v2 + 2*j1*m1*o2*s0*v2 + 2*j2*k1*q1*s0*v2 + 
  2*j1*k2*q1*s0*v2 - 2*i2*o1*q1*s0*v2 - 2*i1*o2*q1*s0*v2 + 
  2*j1*k1*q2*s0*v2 - 2*i1*o1*q2*s0*v2 - 2*k2*m1*n0*s1*v2 - 
  2*k1*m2*n0*s1*v2 - 2*k2*m0*n1*s1*v2 - 2*k0*m2*n1*s1*v2 - 
  2*k1*m0*n2*s1*v2 - 2*k0*m1*n2*s1*v2 + 2*j2*m1*o0*s1*v2 + 
  2*j1*m2*o0*s1*v2 + 2*j2*m0*o1*s1*v2 + 2*j0*m2*o1*s1*v2 + 
  2*j1*m0*o2*s1*v2 + 2*j0*m1*o2*s1*v2 + 2*j2*k1*q0*s1*v2 + 
  2*j1*k2*q0*s1*v2 - 2*i2*o1*q0*s1*v2 - 2*i1*o2*q0*s1*v2 + 
  2*j2*k0*q1*s1*v2 + 2*j0*k2*q1*s1*v2 - 2*i2*o0*q1*s1*v2 - 
  2*i0*o2*q1*s1*v2 + 2*j1*k0*q2*s1*v2 + 2*j0*k1*q2*s1*v2 - 
  2*i1*o0*q2*s1*v2 - 2*i0*o1*q2*s1*v2 - 2*k1*m1*n0*s2*v2 - 
  2*k1*m0*n1*s2*v2 - 2*k0*m1*n1*s2*v2 + 2*j1*m1*o0*s2*v2 + 
  2*j1*m0*o1*s2*v2 + 2*j0*m1*o1*s2*v2 + 2*j1*k1*q0*s2*v2 - 
  2*i1*o1*q0*s2*v2 + 2*j1*k0*q1*s2*v2 + 2*j0*k1*q1*s2*v2 - 
  2*i1*o0*q1*s2*v2 - 2*i0*o1*q1*s2*v2 - 2*k2*l1*n1*t0*v2 - 
  2*k1*l2*n1*t0*v2 - 2*k1*l1*n2*t0*v2 + 2*j2*l1*o1*t0*v2 + 
  2*j1*l2*o1*t0*v2 + 2*j1*l1*o2*t0*v2 + 2*j2*k1*p1*t0*v2 + 
  2*j1*k2*p1*t0*v2 - 2*i2*o1*p1*t0*v2 - 2*i1*o2*p1*t0*v2 + 
  2*j1*k1*p2*t0*v2 - 2*i1*o1*p2*t0*v2 - 4*j1*j2*s1*t0*v2 + 
  2*i2*n1*s1*t0*v2 + 2*i1*n2*s1*t0*v2 - 2*j1*j1*s2*t0*v2 + 
  2*i1*n1*s2*t0*v2 - 2*k2*l1*n0*t1*v2 - 2*k1*l2*n0*t1*v2 - 
  2*k2*l0*n1*t1*v2 - 2*k0*l2*n1*t1*v2 - 2*k1*l0*n2*t1*v2 - 
  2*k0*l1*n2*t1*v2 + 2*j2*l1*o0*t1*v2 + 2*j1*l2*o0*t1*v2 + 
  2*j2*l0*o1*t1*v2 + 2*j0*l2*o1*t1*v2 + 2*j1*l0*o2*t1*v2 + 
  2*j0*l1*o2*t1*v2 + 2*j2*k1*p0*t1*v2 + 2*j1*k2*p0*t1*v2 - 
  2*i2*o1*p0*t1*v2 - 2*i1*o2*p0*t1*v2 + 2*j2*k0*p1*t1*v2 + 
  2*j0*k2*p1*t1*v2 - 2*i2*o0*p1*t1*v2 - 2*i0*o2*p1*t1*v2 + 
  2*j1*k0*p2*t1*v2 + 2*j0*k1*p2*t1*v2 - 2*i1*o0*p2*t1*v2 - 
  2*i0*o1*p2*t1*v2 - 4*j1*j2*s0*t1*v2 + 2*i2*n1*s0*t1*v2 + 
  2*i1*n2*s0*t1*v2 - 4*j0*j2*s1*t1*v2 + 2*i2*n0*s1*t1*v2 + 
  2*i0*n2*s1*t1*v2 - 4*j0*j1*s2*t1*v2 + 2*i1*n0*s2*t1*v2 + 
  2*i0*n1*s2*t1*v2 - 2*k1*l1*n0*t2*v2 - 2*k1*l0*n1*t2*v2 - 
  2*k0*l1*n1*t2*v2 + 2*j1*l1*o0*t2*v2 + 2*j1*l0*o1*t2*v2 + 
  2*j0*l1*o1*t2*v2 + 2*j1*k1*p0*t2*v2 - 2*i1*o1*p0*t2*v2 + 
  2*j1*k0*p1*t2*v2 + 2*j0*k1*p1*t2*v2 - 2*i1*o0*p1*t2*v2 - 
  2*i0*o1*p1*t2*v2 - 2*j1*j1*s0*t2*v2 + 2*i1*n1*s0*t2*v2 - 
  4*j0*j1*s1*t2*v2 + 2*i1*n0*s1*t2*v2 + 2*i0*n1*s1*t2*v2 + 
  4*k1*k2*n1*v0*v2 + 2*k1*k1*n2*v0*v2 - 4*j2*k1*o1*v0*v2 - 
  4*j1*k2*o1*v0*v2 + 2*i2*o1*o1*v0*v2 - 4*j1*k1*o2*v0*v2 + 
  4*i1*o1*o2*v0*v2 + 4*j1*j2*r1*v0*v2 - 2*i2*n1*r1*v0*v2 - 
  2*i1*n2*r1*v0*v2 + 2*j1*j1*r2*v0*v2 - 2*i1*n1*r2*v0*v2 + 
  4*k1*k2*n0*v1*v2 + 4*k0*k2*n1*v1*v2 + 4*k0*k1*n2*v1*v2 - 
  4*j2*k1*o0*v1*v2 - 4*j1*k2*o0*v1*v2 - 4*j2*k0*o1*v1*v2 - 
  4*j0*k2*o1*v1*v2 + 4*i2*o0*o1*v1*v2 - 4*j1*k0*o2*v1*v2 - 
  4*j0*k1*o2*v1*v2 + 4*i1*o0*o2*v1*v2 + 4*i0*o1*o2*v1*v2 + 
  4*j1*j2*r0*v1*v2 - 2*i2*n1*r0*v1*v2 - 2*i1*n2*r0*v1*v2 + 
  4*j0*j2*r1*v1*v2 - 2*i2*n0*r1*v1*v2 - 2*i0*n2*r1*v1*v2 + 
  4*j0*j1*r2*v1*v2 - 2*i1*n0*r2*v1*v2 - 2*i0*n1*r2*v1*v2 + 
  k1*k1*n0*v2*v2 + 2*k0*k1*n1*v2*v2 - 2*j1*k1*o0*v2*v2 - 
  2*j1*k0*o1*v2*v2 - 2*j0*k1*o1*v2*v2 + 2*i1*o0*o1*v2*v2 + 
  i0*o1*o1*v2*v2 + j1*j1*r0*v2*v2 - i1*n1*r0*v2*v2 + 
  2*j0*j1*r1*v2*v2 - i1*n0*r1*v2*v2 - i0*n1*r1*v2*v2 + 
  l2*l2*o1*o1*w0 + 4*l1*l2*o1*o2*w0 + l1*l1*o2*o2*w0 - 
  2*k2*l2*o1*p1*w0 - 2*k2*l1*o2*p1*w0 - 2*k1*l2*o2*p1*w0 + 
  k2*k2*p1*p1*w0 - 2*k2*l1*o1*p2*w0 - 2*k1*l2*o1*p2*w0 - 
  2*k1*l1*o2*p2*w0 + 4*k1*k2*p1*p2*w0 + k1*k1*p2*p2*w0 - 
  l2*l2*n1*r1*w0 - 2*l1*l2*n2*r1*w0 + 2*j2*l2*p1*r1*w0 + 
  2*j2*l1*p2*r1*w0 + 2*j1*l2*p2*r1*w0 - 2*i2*p1*p2*r1*w0 - 
  i1*p2*p2*r1*w0 - 2*l1*l2*n1*r2*w0 - l1*l1*n2*r2*w0 + 
  2*j2*l1*p1*r2*w0 + 2*j1*l2*p1*r2*w0 - i2*p1*p1*r2*w0 + 
  2*j1*l1*p2*r2*w0 - 2*i1*p1*p2*r2*w0 + 2*k2*l2*n1*s1*w0 + 
  2*k2*l1*n2*s1*w0 + 2*k1*l2*n2*s1*w0 - 2*j2*l2*o1*s1*w0 - 
  2*j2*l1*o2*s1*w0 - 2*j1*l2*o2*s1*w0 - 2*j2*k2*p1*s1*w0 + 
  2*i2*o2*p1*s1*w0 - 2*j2*k1*p2*s1*w0 - 2*j1*k2*p2*s1*w0 + 
  2*i2*o1*p2*s1*w0 + 2*i1*o2*p2*s1*w0 + j2*j2*s1*s1*w0 - 
  i2*n2*s1*s1*w0 + 2*k2*l1*n1*s2*w0 + 2*k1*l2*n1*s2*w0 + 
  2*k1*l1*n2*s2*w0 - 2*j2*l1*o1*s2*w0 - 2*j1*l2*o1*s2*w0 - 
  2*j1*l1*o2*s2*w0 - 2*j2*k1*p1*s2*w0 - 2*j1*k2*p1*s2*w0 + 
  2*i2*o1*p1*s2*w0 + 2*i1*o2*p1*s2*w0 - 2*j1*k1*p2*s2*w0 + 
  2*i1*o1*p2*s2*w0 + 4*j1*j2*s1*s2*w0 - 2*i2*n1*s1*s2*w0 - 
  2*i1*n2*s1*s2*w0 + j1*j1*s2*s2*w0 - i1*n1*s2*s2*w0 - 
  k2*k2*n1*u1*w0 - 2*k1*k2*n2*u1*w0 + 2*j2*k2*o1*u1*w0 + 
  2*j2*k1*o2*u1*w0 + 2*j1*k2*o2*u1*w0 - 2*i2*o1*o2*u1*w0 - 
  i1*o2*o2*u1*w0 - j2*j2*r1*u1*w0 + i2*n2*r1*u1*w0 - 
  2*j1*j2*r2*u1*w0 + i2*n1*r2*u1*w0 + i1*n2*r2*u1*w0 - 
  2*k1*k2*n1*u2*w0 - k1*k1*n2*u2*w0 + 2*j2*k1*o1*u2*w0 + 
  2*j1*k2*o1*u2*w0 - i2*o1*o1*u2*w0 + 2*j1*k1*o2*u2*w0 - 
  2*i1*o1*o2*u2*w0 - 2*j1*j2*r1*u2*w0 + i2*n1*r1*u2*w0 + 
  i1*n2*r1*u2*w0 - j1*j1*r2*u2*w0 + i1*n1*r2*u2*w0 + 
  2*l2*l2*o0*o1*w1 + 4*l1*l2*o0*o2*w1 + 4*l0*l2*o1*o2*w1 + 
  2*l0*l1*o2*o2*w1 - 2*k2*l2*o1*p0*w1 - 2*k2*l1*o2*p0*w1 - 
  2*k1*l2*o2*p0*w1 - 2*k2*l2*o0*p1*w1 - 2*k2*l0*o2*p1*w1 - 
  2*k0*l2*o2*p1*w1 + 2*k2*k2*p0*p1*w1 - 2*k2*l1*o0*p2*w1 - 
  2*k1*l2*o0*p2*w1 - 2*k2*l0*o1*p2*w1 - 2*k0*l2*o1*p2*w1 - 
  2*k1*l0*o2*p2*w1 - 2*k0*l1*o2*p2*w1 + 4*k1*k2*p0*p2*w1 + 
  4*k0*k2*p1*p2*w1 + 2*k0*k1*p2*p2*w1 - l2*l2*n1*r0*w1 - 
  2*l1*l2*n2*r0*w1 + 2*j2*l2*p1*r0*w1 + 2*j2*l1*p2*r0*w1 + 
  2*j1*l2*p2*r0*w1 - 2*i2*p1*p2*r0*w1 - i1*p2*p2*r0*w1 - 
  l2*l2*n0*r1*w1 - 2*l0*l2*n2*r1*w1 + 2*j2*l2*p0*r1*w1 + 
  2*j2*l0*p2*r1*w1 + 2*j0*l2*p2*r1*w1 - 2*i2*p0*p2*r1*w1 - 
  i0*p2*p2*r1*w1 - 2*l1*l2*n0*r2*w1 - 2*l0*l2*n1*r2*w1 - 
  2*l0*l1*n2*r2*w1 + 2*j2*l1*p0*r2*w1 + 2*j1*l2*p0*r2*w1 + 
  2*j2*l0*p1*r2*w1 + 2*j0*l2*p1*r2*w1 - 2*i2*p0*p1*r2*w1 + 
  2*j1*l0*p2*r2*w1 + 2*j0*l1*p2*r2*w1 - 2*i1*p0*p2*r2*w1 - 
  2*i0*p1*p2*r2*w1 + 2*k2*l2*n1*s0*w1 + 2*k2*l1*n2*s0*w1 + 
  2*k1*l2*n2*s0*w1 - 2*j2*l2*o1*s0*w1 - 2*j2*l1*o2*s0*w1 - 
  2*j1*l2*o2*s0*w1 - 2*j2*k2*p1*s0*w1 + 2*i2*o2*p1*s0*w1 - 
  2*j2*k1*p2*s0*w1 - 2*j1*k2*p2*s0*w1 + 2*i2*o1*p2*s0*w1 + 
  2*i1*o2*p2*s0*w1 + 2*k2*l2*n0*s1*w1 + 2*k2*l0*n2*s1*w1 + 
  2*k0*l2*n2*s1*w1 - 2*j2*l2*o0*s1*w1 - 2*j2*l0*o2*s1*w1 - 
  2*j0*l2*o2*s1*w1 - 2*j2*k2*p0*s1*w1 + 2*i2*o2*p0*s1*w1 - 
  2*j2*k0*p2*s1*w1 - 2*j0*k2*p2*s1*w1 + 2*i2*o0*p2*s1*w1 + 
  2*i0*o2*p2*s1*w1 + 2*j2*j2*s0*s1*w1 - 2*i2*n2*s0*s1*w1 + 
  2*k2*l1*n0*s2*w1 + 2*k1*l2*n0*s2*w1 + 2*k2*l0*n1*s2*w1 + 
  2*k0*l2*n1*s2*w1 + 2*k1*l0*n2*s2*w1 + 2*k0*l1*n2*s2*w1 - 
  2*j2*l1*o0*s2*w1 - 2*j1*l2*o0*s2*w1 - 2*j2*l0*o1*s2*w1 - 
  2*j0*l2*o1*s2*w1 - 2*j1*l0*o2*s2*w1 - 2*j0*l1*o2*s2*w1 - 
  2*j2*k1*p0*s2*w1 - 2*j1*k2*p0*s2*w1 + 2*i2*o1*p0*s2*w1 + 
  2*i1*o2*p0*s2*w1 - 2*j2*k0*p1*s2*w1 - 2*j0*k2*p1*s2*w1 + 
  2*i2*o0*p1*s2*w1 + 2*i0*o2*p1*s2*w1 - 2*j1*k0*p2*s2*w1 - 
  2*j0*k1*p2*s2*w1 + 2*i1*o0*p2*s2*w1 + 2*i0*o1*p2*s2*w1 + 
  4*j1*j2*s0*s2*w1 - 2*i2*n1*s0*s2*w1 - 2*i1*n2*s0*s2*w1 + 
  4*j0*j2*s1*s2*w1 - 2*i2*n0*s1*s2*w1 - 2*i0*n2*s1*s2*w1 + 
  2*j0*j1*s2*s2*w1 - i1*n0*s2*s2*w1 - i0*n1*s2*s2*w1 - 
  k2*k2*n1*u0*w1 - 2*k1*k2*n2*u0*w1 + 2*j2*k2*o1*u0*w1 + 
  2*j2*k1*o2*u0*w1 + 2*j1*k2*o2*u0*w1 - 2*i2*o1*o2*u0*w1 - 
  i1*o2*o2*u0*w1 - j2*j2*r1*u0*w1 + i2*n2*r1*u0*w1 - 
  2*j1*j2*r2*u0*w1 + i2*n1*r2*u0*w1 + i1*n2*r2*u0*w1 - 
  k2*k2*n0*u1*w1 - 2*k0*k2*n2*u1*w1 + 2*j2*k2*o0*u1*w1 + 
  2*j2*k0*o2*u1*w1 + 2*j0*k2*o2*u1*w1 - 2*i2*o0*o2*u1*w1 - 
  i0*o2*o2*u1*w1 - j2*j2*r0*u1*w1 + i2*n2*r0*u1*w1 - 
  2*j0*j2*r2*u1*w1 + i2*n0*r2*u1*w1 + i0*n2*r2*u1*w1 - 
  2*k1*k2*n0*u2*w1 - 2*k0*k2*n1*u2*w1 - 2*k0*k1*n2*u2*w1 + 
  2*j2*k1*o0*u2*w1 + 2*j1*k2*o0*u2*w1 + 2*j2*k0*o1*u2*w1 + 
  2*j0*k2*o1*u2*w1 - 2*i2*o0*o1*u2*w1 + 2*j1*k0*o2*u2*w1 + 
  2*j0*k1*o2*u2*w1 - 2*i1*o0*o2*u2*w1 - 2*i0*o1*o2*u2*w1 - 
  2*j1*j2*r0*u2*w1 + i2*n1*r0*u2*w1 + i1*n2*r0*u2*w1 - 
  2*j0*j2*r1*u2*w1 + i2*n0*r1*u2*w1 + i0*n2*r1*u2*w1 - 
  2*j0*j1*r2*u2*w1 + i1*n0*r2*u2*w1 + i0*n1*r2*u2*w1 + 
  4*l1*l2*o0*o1*w2 + 2*l0*l2*o1*o1*w2 + 2*l1*l1*o0*o2*w2 + 
  4*l0*l1*o1*o2*w2 - 2*k2*l1*o1*p0*w2 - 2*k1*l2*o1*p0*w2 - 
  2*k1*l1*o2*p0*w2 - 2*k2*l1*o0*p1*w2 - 2*k1*l2*o0*p1*w2 - 
  2*k2*l0*o1*p1*w2 - 2*k0*l2*o1*p1*w2 - 2*k1*l0*o2*p1*w2 - 
  2*k0*l1*o2*p1*w2 + 4*k1*k2*p0*p1*w2 + 2*k0*k2*p1*p1*w2 - 
  2*k1*l1*o0*p2*w2 - 2*k1*l0*o1*p2*w2 - 2*k0*l1*o1*p2*w2 + 
  2*k1*k1*p0*p2*w2 + 4*k0*k1*p1*p2*w2 - 2*l1*l2*n1*r0*w2 - 
  l1*l1*n2*r0*w2 + 2*j2*l1*p1*r0*w2 + 2*j1*l2*p1*r0*w2 - 
  i2*p1*p1*r0*w2 + 2*j1*l1*p2*r0*w2 - 2*i1*p1*p2*r0*w2 - 
  2*l1*l2*n0*r1*w2 - 2*l0*l2*n1*r1*w2 - 2*l0*l1*n2*r1*w2 + 
  2*j2*l1*p0*r1*w2 + 2*j1*l2*p0*r1*w2 + 2*j2*l0*p1*r1*w2 + 
  2*j0*l2*p1*r1*w2 - 2*i2*p0*p1*r1*w2 + 2*j1*l0*p2*r1*w2 + 
  2*j0*l1*p2*r1*w2 - 2*i1*p0*p2*r1*w2 - 2*i0*p1*p2*r1*w2 - 
  l1*l1*n0*r2*w2 - 2*l0*l1*n1*r2*w2 + 2*j1*l1*p0*r2*w2 + 
  2*j1*l0*p1*r2*w2 + 2*j0*l1*p1*r2*w2 - 2*i1*p0*p1*r2*w2 - 
  i0*p1*p1*r2*w2 + 2*k2*l1*n1*s0*w2 + 2*k1*l2*n1*s0*w2 + 
  2*k1*l1*n2*s0*w2 - 2*j2*l1*o1*s0*w2 - 2*j1*l2*o1*s0*w2 - 
  2*j1*l1*o2*s0*w2 - 2*j2*k1*p1*s0*w2 - 2*j1*k2*p1*s0*w2 + 
  2*i2*o1*p1*s0*w2 + 2*i1*o2*p1*s0*w2 - 2*j1*k1*p2*s0*w2 + 
  2*i1*o1*p2*s0*w2 + 2*k2*l1*n0*s1*w2 + 2*k1*l2*n0*s1*w2 + 
  2*k2*l0*n1*s1*w2 + 2*k0*l2*n1*s1*w2 + 2*k1*l0*n2*s1*w2 + 
  2*k0*l1*n2*s1*w2 - 2*j2*l1*o0*s1*w2 - 2*j1*l2*o0*s1*w2 - 
  2*j2*l0*o1*s1*w2 - 2*j0*l2*o1*s1*w2 - 2*j1*l0*o2*s1*w2 - 
  2*j0*l1*o2*s1*w2 - 2*j2*k1*p0*s1*w2 - 2*j1*k2*p0*s1*w2 + 
  2*i2*o1*p0*s1*w2 + 2*i1*o2*p0*s1*w2 - 2*j2*k0*p1*s1*w2 - 
  2*j0*k2*p1*s1*w2 + 2*i2*o0*p1*s1*w2 + 2*i0*o2*p1*s1*w2 - 
  2*j1*k0*p2*s1*w2 - 2*j0*k1*p2*s1*w2 + 2*i1*o0*p2*s1*w2 + 
  2*i0*o1*p2*s1*w2 + 4*j1*j2*s0*s1*w2 - 2*i2*n1*s0*s1*w2 - 
  2*i1*n2*s0*s1*w2 + 2*j0*j2*s1*s1*w2 - i2*n0*s1*s1*w2 - 
  i0*n2*s1*s1*w2 + 2*k1*l1*n0*s2*w2 + 2*k1*l0*n1*s2*w2 + 
  2*k0*l1*n1*s2*w2 - 2*j1*l1*o0*s2*w2 - 2*j1*l0*o1*s2*w2 - 
  2*j0*l1*o1*s2*w2 - 2*j1*k1*p0*s2*w2 + 2*i1*o1*p0*s2*w2 - 
  2*j1*k0*p1*s2*w2 - 2*j0*k1*p1*s2*w2 + 2*i1*o0*p1*s2*w2 + 
  2*i0*o1*p1*s2*w2 + 2*j1*j1*s0*s2*w2 - 2*i1*n1*s0*s2*w2 + 
  4*j0*j1*s1*s2*w2 - 2*i1*n0*s1*s2*w2 - 2*i0*n1*s1*s2*w2 - 
  2*k1*k2*n1*u0*w2 - k1*k1*n2*u0*w2 + 2*j2*k1*o1*u0*w2 + 
  2*j1*k2*o1*u0*w2 - i2*o1*o1*u0*w2 + 2*j1*k1*o2*u0*w2 - 
  2*i1*o1*o2*u0*w2 - 2*j1*j2*r1*u0*w2 + i2*n1*r1*u0*w2 + 
  i1*n2*r1*u0*w2 - j1*j1*r2*u0*w2 + i1*n1*r2*u0*w2 - 
  2*k1*k2*n0*u1*w2 - 2*k0*k2*n1*u1*w2 - 2*k0*k1*n2*u1*w2 + 
  2*j2*k1*o0*u1*w2 + 2*j1*k2*o0*u1*w2 + 2*j2*k0*o1*u1*w2 + 
  2*j0*k2*o1*u1*w2 - 2*i2*o0*o1*u1*w2 + 2*j1*k0*o2*u1*w2 + 
  2*j0*k1*o2*u1*w2 - 2*i1*o0*o2*u1*w2 - 2*i0*o1*o2*u1*w2 - 
  2*j1*j2*r0*u1*w2 + i2*n1*r0*u1*w2 + i1*n2*r0*u1*w2 - 
  2*j0*j2*r1*u1*w2 + i2*n0*r1*u1*w2 + i0*n2*r1*u1*w2 - 
  2*j0*j1*r2*u1*w2 + i1*n0*r2*u1*w2 + i0*n1*r2*u1*w2 - 
  k1*k1*n0*u2*w2 - 2*k0*k1*n1*u2*w2 + 2*j1*k1*o0*u2*w2 + 
  2*j1*k0*o1*u2*w2 + 2*j0*k1*o1*u2*w2 - 2*i1*o0*o1*u2*w2 - 
  i0*o1*o1*u2*w2 - j1*j1*r0*u2*w2 + i1*n1*r0*u2*w2 - 
  2*j0*j1*r1*u2*w2 + i1*n0*r1*u2*w2 + i0*n1*r1*u2*w2;
    
    k[9]  = m2*m2*p1*p1*r1 + 4*m1*m2*p1*p2*r1 + m1*m1*p2*p2*r1 - 
  2*l2*m2*p1*q1*r1 - 2*l2*m1*p2*q1*r1 - 2*l1*m2*p2*q1*r1 + 
  l2*l2*q1*q1*r1 - 2*l2*m1*p1*q2*r1 - 2*l1*m2*p1*q2*r1 - 
  2*l1*m1*p2*q2*r1 + 4*l1*l2*q1*q2*r1 + l1*l1*q2*q2*r1 + 
  2*m1*m2*p1*p1*r2 + 2*m1*m1*p1*p2*r2 - 2*l2*m1*p1*q1*r2 - 
  2*l1*m2*p1*q1*r2 - 2*l1*m1*p2*q1*r2 + 2*l1*l2*q1*q1*r2 - 
  2*l1*m1*p1*q2*r2 + 2*l1*l1*q1*q2*r2 - 2*m2*m2*o1*p1*s1 - 
  4*m1*m2*o2*p1*s1 - 4*m1*m2*o1*p2*s1 - 2*m1*m1*o2*p2*s1 + 
  2*l2*m2*o1*q1*s1 + 2*l2*m1*o2*q1*s1 + 2*l1*m2*o2*q1*s1 + 
  2*k2*m2*p1*q1*s1 + 2*k2*m1*p2*q1*s1 + 2*k1*m2*p2*q1*s1 - 
  2*k2*l2*q1*q1*s1 + 2*l2*m1*o1*q2*s1 + 2*l1*m2*o1*q2*s1 + 
  2*l1*m1*o2*q2*s1 + 2*k2*m1*p1*q2*s1 + 2*k1*m2*p1*q2*s1 + 
  2*k1*m1*p2*q2*s1 - 4*k2*l1*q1*q2*s1 - 4*k1*l2*q1*q2*s1 - 
  2*k1*l1*q2*q2*s1 + m2*m2*n1*s1*s1 + 2*m1*m2*n2*s1*s1 - 
  2*j2*m2*q1*s1*s1 - 2*j2*m1*q2*s1*s1 - 2*j1*m2*q2*s1*s1 + 
  2*i2*q1*q2*s1*s1 + i1*q2*q2*s1*s1 - 4*m1*m2*o1*p1*s2 - 
  2*m1*m1*o2*p1*s2 - 2*m1*m1*o1*p2*s2 + 2*l2*m1*o1*q1*s2 + 
  2*l1*m2*o1*q1*s2 + 2*l1*m1*o2*q1*s2 + 2*k2*m1*p1*q1*s2 + 
  2*k1*m2*p1*q1*s2 + 2*k1*m1*p2*q1*s2 - 2*k2*l1*q1*q1*s2 - 
  2*k1*l2*q1*q1*s2 + 2*l1*m1*o1*q2*s2 + 2*k1*m1*p1*q2*s2 - 
  4*k1*l1*q1*q2*s2 + 4*m1*m2*n1*s1*s2 + 2*m1*m1*n2*s1*s2 - 
  4*j2*m1*q1*s1*s2 - 4*j1*m2*q1*s1*s2 + 2*i2*q1*q1*s1*s2 - 
  4*j1*m1*q2*s1*s2 + 4*i1*q1*q2*s1*s2 + m1*m1*n1*s2*s2 - 
  2*j1*m1*q1*s2*s2 + i1*q1*q1*s2*s2 + 2*l2*m2*o1*p1*t1 + 
  2*l2*m1*o2*p1*t1 + 2*l1*m2*o2*p1*t1 - 2*k2*m2*p1*p1*t1 + 
  2*l2*m1*o1*p2*t1 + 2*l1*m2*o1*p2*t1 + 2*l1*m1*o2*p2*t1 - 
  4*k2*m1*p1*p2*t1 - 4*k1*m2*p1*p2*t1 - 2*k1*m1*p2*p2*t1 - 
  2*l2*l2*o1*q1*t1 - 4*l1*l2*o2*q1*t1 + 2*k2*l2*p1*q1*t1 + 
  2*k2*l1*p2*q1*t1 + 2*k1*l2*p2*q1*t1 - 4*l1*l2*o1*q2*t1 - 
  2*l1*l1*o2*q2*t1 + 2*k2*l1*p1*q2*t1 + 2*k1*l2*p1*q2*t1 + 
  2*k1*l1*p2*q2*t1 - 2*l2*m2*n1*s1*t1 - 2*l2*m1*n2*s1*t1 - 
  2*l1*m2*n2*s1*t1 + 2*j2*m2*p1*s1*t1 + 2*j2*m1*p2*s1*t1 + 
  2*j1*m2*p2*s1*t1 + 2*j2*l2*q1*s1*t1 - 2*i2*p2*q1*s1*t1 + 
  2*j2*l1*q2*s1*t1 + 2*j1*l2*q2*s1*t1 - 2*i2*p1*q2*s1*t1 - 
  2*i1*p2*q2*s1*t1 - 2*l2*m1*n1*s2*t1 - 2*l1*m2*n1*s2*t1 - 
  2*l1*m1*n2*s2*t1 + 2*j2*m1*p1*s2*t1 + 2*j1*m2*p1*s2*t1 + 
  2*j1*m1*p2*s2*t1 + 2*j2*l1*q1*s2*t1 + 2*j1*l2*q1*s2*t1 - 
  2*i2*p1*q1*s2*t1 - 2*i1*p2*q1*s2*t1 + 2*j1*l1*q2*s2*t1 - 
  2*i1*p1*q2*s2*t1 + l2*l2*n1*t1*t1 + 2*l1*l2*n2*t1*t1 - 
  2*j2*l2*p1*t1*t1 - 2*j2*l1*p2*t1*t1 - 2*j1*l2*p2*t1*t1 + 
  2*i2*p1*p2*t1*t1 + i1*p2*p2*t1*t1 + 2*l2*m1*o1*p1*t2 + 
  2*l1*m2*o1*p1*t2 + 2*l1*m1*o2*p1*t2 - 2*k2*m1*p1*p1*t2 - 
  2*k1*m2*p1*p1*t2 + 2*l1*m1*o1*p2*t2 - 4*k1*m1*p1*p2*t2 - 
  4*l1*l2*o1*q1*t2 - 2*l1*l1*o2*q1*t2 + 2*k2*l1*p1*q1*t2 + 
  2*k1*l2*p1*q1*t2 + 2*k1*l1*p2*q1*t2 - 2*l1*l1*o1*q2*t2 + 
  2*k1*l1*p1*q2*t2 - 2*l2*m1*n1*s1*t2 - 2*l1*m2*n1*s1*t2 - 
  2*l1*m1*n2*s1*t2 + 2*j2*m1*p1*s1*t2 + 2*j1*m2*p1*s1*t2 + 
  2*j1*m1*p2*s1*t2 + 2*j2*l1*q1*s1*t2 + 2*j1*l2*q1*s1*t2 - 
  2*i2*p1*q1*s1*t2 - 2*i1*p2*q1*s1*t2 + 2*j1*l1*q2*s1*t2 - 
  2*i1*p1*q2*s1*t2 - 2*l1*m1*n1*s2*t2 + 2*j1*m1*p1*s2*t2 + 
  2*j1*l1*q1*s2*t2 - 2*i1*p1*q1*s2*t2 + 4*l1*l2*n1*t1*t2 + 
  2*l1*l1*n2*t1*t2 - 4*j2*l1*p1*t1*t2 - 4*j1*l2*p1*t1*t2 + 
  2*i2*p1*p1*t1*t2 - 4*j1*l1*p2*t1*t2 + 4*i1*p1*p2*t1*t2 + 
  l1*l1*n1*t2*t2 - 2*j1*l1*p1*t2*t2 + i1*p1*p1*t2*t2 + 
  m2*m2*o1*o1*u1 + 4*m1*m2*o1*o2*u1 + m1*m1*o2*o2*u1 - 
  2*k2*m2*o1*q1*u1 - 2*k2*m1*o2*q1*u1 - 2*k1*m2*o2*q1*u1 + 
  k2*k2*q1*q1*u1 - 2*k2*m1*o1*q2*u1 - 2*k1*m2*o1*q2*u1 - 
  2*k1*m1*o2*q2*u1 + 4*k1*k2*q1*q2*u1 + k1*k1*q2*q2*u1 - 
  m2*m2*n1*r1*u1 - 2*m1*m2*n2*r1*u1 + 2*j2*m2*q1*r1*u1 + 
  2*j2*m1*q2*r1*u1 + 2*j1*m2*q2*r1*u1 - 2*i2*q1*q2*r1*u1 - 
  i1*q2*q2*r1*u1 - 2*m1*m2*n1*r2*u1 - m1*m1*n2*r2*u1 + 
  2*j2*m1*q1*r2*u1 + 2*j1*m2*q1*r2*u1 - i2*q1*q1*r2*u1 + 
  2*j1*m1*q2*r2*u1 - 2*i1*q1*q2*r2*u1 + 2*k2*m2*n1*t1*u1 + 
  2*k2*m1*n2*t1*u1 + 2*k1*m2*n2*t1*u1 - 2*j2*m2*o1*t1*u1 - 
  2*j2*m1*o2*t1*u1 - 2*j1*m2*o2*t1*u1 - 2*j2*k2*q1*t1*u1 + 
  2*i2*o2*q1*t1*u1 - 2*j2*k1*q2*t1*u1 - 2*j1*k2*q2*t1*u1 + 
  2*i2*o1*q2*t1*u1 + 2*i1*o2*q2*t1*u1 + j2*j2*t1*t1*u1 - 
  i2*n2*t1*t1*u1 + 2*k2*m1*n1*t2*u1 + 2*k1*m2*n1*t2*u1 + 
  2*k1*m1*n2*t2*u1 - 2*j2*m1*o1*t2*u1 - 2*j1*m2*o1*t2*u1 - 
  2*j1*m1*o2*t2*u1 - 2*j2*k1*q1*t2*u1 - 2*j1*k2*q1*t2*u1 + 
  2*i2*o1*q1*t2*u1 + 2*i1*o2*q1*t2*u1 - 2*j1*k1*q2*t2*u1 + 
  2*i1*o1*q2*t2*u1 + 4*j1*j2*t1*t2*u1 - 2*i2*n1*t1*t2*u1 - 
  2*i1*n2*t1*t2*u1 + j1*j1*t2*t2*u1 - i1*n1*t2*t2*u1 + 
  2*m1*m2*o1*o1*u2 + 2*m1*m1*o1*o2*u2 - 2*k2*m1*o1*q1*u2 - 
  2*k1*m2*o1*q1*u2 - 2*k1*m1*o2*q1*u2 + 2*k1*k2*q1*q1*u2 - 
  2*k1*m1*o1*q2*u2 + 2*k1*k1*q1*q2*u2 - 2*m1*m2*n1*r1*u2 - 
  m1*m1*n2*r1*u2 + 2*j2*m1*q1*r1*u2 + 2*j1*m2*q1*r1*u2 - 
  i2*q1*q1*r1*u2 + 2*j1*m1*q2*r1*u2 - 2*i1*q1*q2*r1*u2 - 
  m1*m1*n1*r2*u2 + 2*j1*m1*q1*r2*u2 - i1*q1*q1*r2*u2 + 
  2*k2*m1*n1*t1*u2 + 2*k1*m2*n1*t1*u2 + 2*k1*m1*n2*t1*u2 - 
  2*j2*m1*o1*t1*u2 - 2*j1*m2*o1*t1*u2 - 2*j1*m1*o2*t1*u2 - 
  2*j2*k1*q1*t1*u2 - 2*j1*k2*q1*t1*u2 + 2*i2*o1*q1*t1*u2 + 
  2*i1*o2*q1*t1*u2 - 2*j1*k1*q2*t1*u2 + 2*i1*o1*q2*t1*u2 + 
  2*j1*j2*t1*t1*u2 - i2*n1*t1*t1*u2 - i1*n2*t1*t1*u2 + 
  2*k1*m1*n1*t2*u2 - 2*j1*m1*o1*t2*u2 - 2*j1*k1*q1*t2*u2 + 
  2*i1*o1*q1*t2*u2 + 2*j1*j1*t1*t2*u2 - 2*i1*n1*t1*t2*u2 - 
  2*l2*m2*o1*o1*v1 - 4*l2*m1*o1*o2*v1 - 4*l1*m2*o1*o2*v1 - 
  2*l1*m1*o2*o2*v1 + 2*k2*m2*o1*p1*v1 + 2*k2*m1*o2*p1*v1 + 
  2*k1*m2*o2*p1*v1 + 2*k2*m1*o1*p2*v1 + 2*k1*m2*o1*p2*v1 + 
  2*k1*m1*o2*p2*v1 + 2*k2*l2*o1*q1*v1 + 2*k2*l1*o2*q1*v1 + 
  2*k1*l2*o2*q1*v1 - 2*k2*k2*p1*q1*v1 - 4*k1*k2*p2*q1*v1 + 
  2*k2*l1*o1*q2*v1 + 2*k1*l2*o1*q2*v1 + 2*k1*l1*o2*q2*v1 - 
  4*k1*k2*p1*q2*v1 - 2*k1*k1*p2*q2*v1 + 2*l2*m2*n1*r1*v1 + 
  2*l2*m1*n2*r1*v1 + 2*l1*m2*n2*r1*v1 - 2*j2*m2*p1*r1*v1 - 
  2*j2*m1*p2*r1*v1 - 2*j1*m2*p2*r1*v1 - 2*j2*l2*q1*r1*v1 + 
  2*i2*p2*q1*r1*v1 - 2*j2*l1*q2*r1*v1 - 2*j1*l2*q2*r1*v1 + 
  2*i2*p1*q2*r1*v1 + 2*i1*p2*q2*r1*v1 + 2*l2*m1*n1*r2*v1 + 
  2*l1*m2*n1*r2*v1 + 2*l1*m1*n2*r2*v1 - 2*j2*m1*p1*r2*v1 - 
  2*j1*m2*p1*r2*v1 - 2*j1*m1*p2*r2*v1 - 2*j2*l1*q1*r2*v1 - 
  2*j1*l2*q1*r2*v1 + 2*i2*p1*q1*r2*v1 + 2*i1*p2*q1*r2*v1 - 
  2*j1*l1*q2*r2*v1 + 2*i1*p1*q2*r2*v1 - 2*k2*m2*n1*s1*v1 - 
  2*k2*m1*n2*s1*v1 - 2*k1*m2*n2*s1*v1 + 2*j2*m2*o1*s1*v1 + 
  2*j2*m1*o2*s1*v1 + 2*j1*m2*o2*s1*v1 + 2*j2*k2*q1*s1*v1 - 
  2*i2*o2*q1*s1*v1 + 2*j2*k1*q2*s1*v1 + 2*j1*k2*q2*s1*v1 - 
  2*i2*o1*q2*s1*v1 - 2*i1*o2*q2*s1*v1 - 2*k2*m1*n1*s2*v1 - 
  2*k1*m2*n1*s2*v1 - 2*k1*m1*n2*s2*v1 + 2*j2*m1*o1*s2*v1 + 
  2*j1*m2*o1*s2*v1 + 2*j1*m1*o2*s2*v1 + 2*j2*k1*q1*s2*v1 + 
  2*j1*k2*q1*s2*v1 - 2*i2*o1*q1*s2*v1 - 2*i1*o2*q1*s2*v1 + 
  2*j1*k1*q2*s2*v1 - 2*i1*o1*q2*s2*v1 - 2*k2*l2*n1*t1*v1 - 
  2*k2*l1*n2*t1*v1 - 2*k1*l2*n2*t1*v1 + 2*j2*l2*o1*t1*v1 + 
  2*j2*l1*o2*t1*v1 + 2*j1*l2*o2*t1*v1 + 2*j2*k2*p1*t1*v1 - 
  2*i2*o2*p1*t1*v1 + 2*j2*k1*p2*t1*v1 + 2*j1*k2*p2*t1*v1 - 
  2*i2*o1*p2*t1*v1 - 2*i1*o2*p2*t1*v1 - 2*j2*j2*s1*t1*v1 + 
  2*i2*n2*s1*t1*v1 - 4*j1*j2*s2*t1*v1 + 2*i2*n1*s2*t1*v1 + 
  2*i1*n2*s2*t1*v1 - 2*k2*l1*n1*t2*v1 - 2*k1*l2*n1*t2*v1 - 
  2*k1*l1*n2*t2*v1 + 2*j2*l1*o1*t2*v1 + 2*j1*l2*o1*t2*v1 + 
  2*j1*l1*o2*t2*v1 + 2*j2*k1*p1*t2*v1 + 2*j1*k2*p1*t2*v1 - 
  2*i2*o1*p1*t2*v1 - 2*i1*o2*p1*t2*v1 + 2*j1*k1*p2*t2*v1 - 
  2*i1*o1*p2*t2*v1 - 4*j1*j2*s1*t2*v1 + 2*i2*n1*s1*t2*v1 + 
  2*i1*n2*s1*t2*v1 - 2*j1*j1*s2*t2*v1 + 2*i1*n1*s2*t2*v1 + 
  k2*k2*n1*v1*v1 + 2*k1*k2*n2*v1*v1 - 2*j2*k2*o1*v1*v1 - 
  2*j2*k1*o2*v1*v1 - 2*j1*k2*o2*v1*v1 + 2*i2*o1*o2*v1*v1 + 
  i1*o2*o2*v1*v1 + j2*j2*r1*v1*v1 - i2*n2*r1*v1*v1 + 
  2*j1*j2*r2*v1*v1 - i2*n1*r2*v1*v1 - i1*n2*r2*v1*v1 - 
  2*l2*m1*o1*o1*v2 - 2*l1*m2*o1*o1*v2 - 4*l1*m1*o1*o2*v2 + 
  2*k2*m1*o1*p1*v2 + 2*k1*m2*o1*p1*v2 + 2*k1*m1*o2*p1*v2 + 
  2*k1*m1*o1*p2*v2 + 2*k2*l1*o1*q1*v2 + 2*k1*l2*o1*q1*v2 + 
  2*k1*l1*o2*q1*v2 - 4*k1*k2*p1*q1*v2 - 2*k1*k1*p2*q1*v2 + 
  2*k1*l1*o1*q2*v2 - 2*k1*k1*p1*q2*v2 + 2*l2*m1*n1*r1*v2 + 
  2*l1*m2*n1*r1*v2 + 2*l1*m1*n2*r1*v2 - 2*j2*m1*p1*r1*v2 - 
  2*j1*m2*p1*r1*v2 - 2*j1*m1*p2*r1*v2 - 2*j2*l1*q1*r1*v2 - 
  2*j1*l2*q1*r1*v2 + 2*i2*p1*q1*r1*v2 + 2*i1*p2*q1*r1*v2 - 
  2*j1*l1*q2*r1*v2 + 2*i1*p1*q2*r1*v2 + 2*l1*m1*n1*r2*v2 - 
  2*j1*m1*p1*r2*v2 - 2*j1*l1*q1*r2*v2 + 2*i1*p1*q1*r2*v2 - 
  2*k2*m1*n1*s1*v2 - 2*k1*m2*n1*s1*v2 - 2*k1*m1*n2*s1*v2 + 
  2*j2*m1*o1*s1*v2 + 2*j1*m2*o1*s1*v2 + 2*j1*m1*o2*s1*v2 + 
  2*j2*k1*q1*s1*v2 + 2*j1*k2*q1*s1*v2 - 2*i2*o1*q1*s1*v2 - 
  2*i1*o2*q1*s1*v2 + 2*j1*k1*q2*s1*v2 - 2*i1*o1*q2*s1*v2 - 
  2*k1*m1*n1*s2*v2 + 2*j1*m1*o1*s2*v2 + 2*j1*k1*q1*s2*v2 - 
  2*i1*o1*q1*s2*v2 - 2*k2*l1*n1*t1*v2 - 2*k1*l2*n1*t1*v2 - 
  2*k1*l1*n2*t1*v2 + 2*j2*l1*o1*t1*v2 + 2*j1*l2*o1*t1*v2 + 
  2*j1*l1*o2*t1*v2 + 2*j2*k1*p1*t1*v2 + 2*j1*k2*p1*t1*v2 - 
  2*i2*o1*p1*t1*v2 - 2*i1*o2*p1*t1*v2 + 2*j1*k1*p2*t1*v2 - 
  2*i1*o1*p2*t1*v2 - 4*j1*j2*s1*t1*v2 + 2*i2*n1*s1*t1*v2 + 
  2*i1*n2*s1*t1*v2 - 2*j1*j1*s2*t1*v2 + 2*i1*n1*s2*t1*v2 - 
  2*k1*l1*n1*t2*v2 + 2*j1*l1*o1*t2*v2 + 2*j1*k1*p1*t2*v2 - 
  2*i1*o1*p1*t2*v2 - 2*j1*j1*s1*t2*v2 + 2*i1*n1*s1*t2*v2 + 
  4*k1*k2*n1*v1*v2 + 2*k1*k1*n2*v1*v2 - 4*j2*k1*o1*v1*v2 - 
  4*j1*k2*o1*v1*v2 + 2*i2*o1*o1*v1*v2 - 4*j1*k1*o2*v1*v2 + 
  4*i1*o1*o2*v1*v2 + 4*j1*j2*r1*v1*v2 - 2*i2*n1*r1*v1*v2 - 
  2*i1*n2*r1*v1*v2 + 2*j1*j1*r2*v1*v2 - 2*i1*n1*r2*v1*v2 + 
  k1*k1*n1*v2*v2 - 2*j1*k1*o1*v2*v2 + i1*o1*o1*v2*v2 + 
  j1*j1*r1*v2*v2 - i1*n1*r1*v2*v2 + l2*l2*o1*o1*w1 + 
  4*l1*l2*o1*o2*w1 + l1*l1*o2*o2*w1 - 2*k2*l2*o1*p1*w1 - 
  2*k2*l1*o2*p1*w1 - 2*k1*l2*o2*p1*w1 + k2*k2*p1*p1*w1 - 
  2*k2*l1*o1*p2*w1 - 2*k1*l2*o1*p2*w1 - 2*k1*l1*o2*p2*w1 + 
  4*k1*k2*p1*p2*w1 + k1*k1*p2*p2*w1 - l2*l2*n1*r1*w1 - 
  2*l1*l2*n2*r1*w1 + 2*j2*l2*p1*r1*w1 + 2*j2*l1*p2*r1*w1 + 
  2*j1*l2*p2*r1*w1 - 2*i2*p1*p2*r1*w1 - i1*p2*p2*r1*w1 - 
  2*l1*l2*n1*r2*w1 - l1*l1*n2*r2*w1 + 2*j2*l1*p1*r2*w1 + 
  2*j1*l2*p1*r2*w1 - i2*p1*p1*r2*w1 + 2*j1*l1*p2*r2*w1 - 
  2*i1*p1*p2*r2*w1 + 2*k2*l2*n1*s1*w1 + 2*k2*l1*n2*s1*w1 + 
  2*k1*l2*n2*s1*w1 - 2*j2*l2*o1*s1*w1 - 2*j2*l1*o2*s1*w1 - 
  2*j1*l2*o2*s1*w1 - 2*j2*k2*p1*s1*w1 + 2*i2*o2*p1*s1*w1 - 
  2*j2*k1*p2*s1*w1 - 2*j1*k2*p2*s1*w1 + 2*i2*o1*p2*s1*w1 + 
  2*i1*o2*p2*s1*w1 + j2*j2*s1*s1*w1 - i2*n2*s1*s1*w1 + 
  2*k2*l1*n1*s2*w1 + 2*k1*l2*n1*s2*w1 + 2*k1*l1*n2*s2*w1 - 
  2*j2*l1*o1*s2*w1 - 2*j1*l2*o1*s2*w1 - 2*j1*l1*o2*s2*w1 - 
  2*j2*k1*p1*s2*w1 - 2*j1*k2*p1*s2*w1 + 2*i2*o1*p1*s2*w1 + 
  2*i1*o2*p1*s2*w1 - 2*j1*k1*p2*s2*w1 + 2*i1*o1*p2*s2*w1 + 
  4*j1*j2*s1*s2*w1 - 2*i2*n1*s1*s2*w1 - 2*i1*n2*s1*s2*w1 + 
  j1*j1*s2*s2*w1 - i1*n1*s2*s2*w1 - k2*k2*n1*u1*w1 - 
  2*k1*k2*n2*u1*w1 + 2*j2*k2*o1*u1*w1 + 2*j2*k1*o2*u1*w1 + 
  2*j1*k2*o2*u1*w1 - 2*i2*o1*o2*u1*w1 - i1*o2*o2*u1*w1 - 
  j2*j2*r1*u1*w1 + i2*n2*r1*u1*w1 - 2*j1*j2*r2*u1*w1 + 
  i2*n1*r2*u1*w1 + i1*n2*r2*u1*w1 - 2*k1*k2*n1*u2*w1 - 
  k1*k1*n2*u2*w1 + 2*j2*k1*o1*u2*w1 + 2*j1*k2*o1*u2*w1 - 
  i2*o1*o1*u2*w1 + 2*j1*k1*o2*u2*w1 - 2*i1*o1*o2*u2*w1 - 
  2*j1*j2*r1*u2*w1 + i2*n1*r1*u2*w1 + i1*n2*r1*u2*w1 - 
  j1*j1*r2*u2*w1 + i1*n1*r2*u2*w1 + 2*l1*l2*o1*o1*w2 + 
  2*l1*l1*o1*o2*w2 - 2*k2*l1*o1*p1*w2 - 2*k1*l2*o1*p1*w2 - 
  2*k1*l1*o2*p1*w2 + 2*k1*k2*p1*p1*w2 - 2*k1*l1*o1*p2*w2 + 
  2*k1*k1*p1*p2*w2 - 2*l1*l2*n1*r1*w2 - l1*l1*n2*r1*w2 + 
  2*j2*l1*p1*r1*w2 + 2*j1*l2*p1*r1*w2 - i2*p1*p1*r1*w2 + 
  2*j1*l1*p2*r1*w2 - 2*i1*p1*p2*r1*w2 - l1*l1*n1*r2*w2 + 
  2*j1*l1*p1*r2*w2 - i1*p1*p1*r2*w2 + 2*k2*l1*n1*s1*w2 + 
  2*k1*l2*n1*s1*w2 + 2*k1*l1*n2*s1*w2 - 2*j2*l1*o1*s1*w2 - 
  2*j1*l2*o1*s1*w2 - 2*j1*l1*o2*s1*w2 - 2*j2*k1*p1*s1*w2 - 
  2*j1*k2*p1*s1*w2 + 2*i2*o1*p1*s1*w2 + 2*i1*o2*p1*s1*w2 - 
  2*j1*k1*p2*s1*w2 + 2*i1*o1*p2*s1*w2 + 2*j1*j2*s1*s1*w2 - 
  i2*n1*s1*s1*w2 - i1*n2*s1*s1*w2 + 2*k1*l1*n1*s2*w2 - 
  2*j1*l1*o1*s2*w2 - 2*j1*k1*p1*s2*w2 + 2*i1*o1*p1*s2*w2 + 
  2*j1*j1*s1*s2*w2 - 2*i1*n1*s1*s2*w2 - 2*k1*k2*n1*u1*w2 - 
  k1*k1*n2*u1*w2 + 2*j2*k1*o1*u1*w2 + 2*j1*k2*o1*u1*w2 - 
  i2*o1*o1*u1*w2 + 2*j1*k1*o2*u1*w2 - 2*i1*o1*o2*u1*w2 - 
  2*j1*j2*r1*u1*w2 + i2*n1*r1*u1*w2 + i1*n2*r1*u1*w2 - 
  j1*j1*r2*u1*w2 + i1*n1*r2*u1*w2 - k1*k1*n1*u2*w2 + 
  2*j1*k1*o1*u2*w2 - i1*o1*o1*u2*w2 - j1*j1*r1*u2*w2 + 
  i1*n1*r1*u2*w2;
    
    k[10] = 2*m0*m2*p0*p0*r0 + 2*m0*m0*p0*p2*r0 - 
  2*l2*m0*p0*q0*r0 - 2*l0*m2*p0*q0*r0 - 2*l0*m0*p2*q0*r0 + 
  2*l0*l2*q0*q0*r0 - 2*l0*m0*p0*q2*r0 + 2*l0*l0*q0*q2*r0 + 
  m0*m0*p0*p0*r2 - 2*l0*m0*p0*q0*r2 + l0*l0*q0*q0*r2 - 
  4*m0*m2*o0*p0*s0 - 2*m0*m0*o2*p0*s0 - 2*m0*m0*o0*p2*s0 + 
  2*l2*m0*o0*q0*s0 + 2*l0*m2*o0*q0*s0 + 2*l0*m0*o2*q0*s0 + 
  2*k2*m0*p0*q0*s0 + 2*k0*m2*p0*q0*s0 + 2*k0*m0*p2*q0*s0 - 
  2*k2*l0*q0*q0*s0 - 2*k0*l2*q0*q0*s0 + 2*l0*m0*o0*q2*s0 + 
  2*k0*m0*p0*q2*s0 - 4*k0*l0*q0*q2*s0 + 2*m0*m2*n0*s0*s0 + 
  m0*m0*n2*s0*s0 - 2*j2*m0*q0*s0*s0 - 2*j0*m2*q0*s0*s0 + 
  i2*q0*q0*s0*s0 - 2*j0*m0*q2*s0*s0 + 2*i0*q0*q2*s0*s0 - 
  2*m0*m0*o0*p0*s2 + 2*l0*m0*o0*q0*s2 + 2*k0*m0*p0*q0*s2 - 
  2*k0*l0*q0*q0*s2 + 2*m0*m0*n0*s0*s2 - 4*j0*m0*q0*s0*s2 + 
  2*i0*q0*q0*s0*s2 + 2*l2*m0*o0*p0*t0 + 2*l0*m2*o0*p0*t0 + 
  2*l0*m0*o2*p0*t0 - 2*k2*m0*p0*p0*t0 - 2*k0*m2*p0*p0*t0 + 
  2*l0*m0*o0*p2*t0 - 4*k0*m0*p0*p2*t0 - 4*l0*l2*o0*q0*t0 - 
  2*l0*l0*o2*q0*t0 + 2*k2*l0*p0*q0*t0 + 2*k0*l2*p0*q0*t0 + 
  2*k0*l0*p2*q0*t0 - 2*l0*l0*o0*q2*t0 + 2*k0*l0*p0*q2*t0 - 
  2*l2*m0*n0*s0*t0 - 2*l0*m2*n0*s0*t0 - 2*l0*m0*n2*s0*t0 + 
  2*j2*m0*p0*s0*t0 + 2*j0*m2*p0*s0*t0 + 2*j0*m0*p2*s0*t0 + 
  2*j2*l0*q0*s0*t0 + 2*j0*l2*q0*s0*t0 - 2*i2*p0*q0*s0*t0 - 
  2*i0*p2*q0*s0*t0 + 2*j0*l0*q2*s0*t0 - 2*i0*p0*q2*s0*t0 - 
  2*l0*m0*n0*s2*t0 + 2*j0*m0*p0*s2*t0 + 2*j0*l0*q0*s2*t0 - 
  2*i0*p0*q0*s2*t0 + 2*l0*l2*n0*t0*t0 + l0*l0*n2*t0*t0 - 
  2*j2*l0*p0*t0*t0 - 2*j0*l2*p0*t0*t0 + i2*p0*p0*t0*t0 - 
  2*j0*l0*p2*t0*t0 + 2*i0*p0*p2*t0*t0 + 2*l0*m0*o0*p0*t2 - 
  2*k0*m0*p0*p0*t2 - 2*l0*l0*o0*q0*t2 + 2*k0*l0*p0*q0*t2 - 
  2*l0*m0*n0*s0*t2 + 2*j0*m0*p0*s0*t2 + 2*j0*l0*q0*s0*t2 - 
  2*i0*p0*q0*s0*t2 + 2*l0*l0*n0*t0*t2 - 4*j0*l0*p0*t0*t2 + 
  2*i0*p0*p0*t0*t2 + 2*m0*m2*o0*o0*u0 + 2*m0*m0*o0*o2*u0 - 
  2*k2*m0*o0*q0*u0 - 2*k0*m2*o0*q0*u0 - 2*k0*m0*o2*q0*u0 + 
  2*k0*k2*q0*q0*u0 - 2*k0*m0*o0*q2*u0 + 2*k0*k0*q0*q2*u0 - 
  2*m0*m2*n0*r0*u0 - m0*m0*n2*r0*u0 + 2*j2*m0*q0*r0*u0 + 
  2*j0*m2*q0*r0*u0 - i2*q0*q0*r0*u0 + 2*j0*m0*q2*r0*u0 - 
  2*i0*q0*q2*r0*u0 - m0*m0*n0*r2*u0 + 2*j0*m0*q0*r2*u0 - 
  i0*q0*q0*r2*u0 + 2*k2*m0*n0*t0*u0 + 2*k0*m2*n0*t0*u0 + 
  2*k0*m0*n2*t0*u0 - 2*j2*m0*o0*t0*u0 - 2*j0*m2*o0*t0*u0 - 
  2*j0*m0*o2*t0*u0 - 2*j2*k0*q0*t0*u0 - 2*j0*k2*q0*t0*u0 + 
  2*i2*o0*q0*t0*u0 + 2*i0*o2*q0*t0*u0 - 2*j0*k0*q2*t0*u0 + 
  2*i0*o0*q2*t0*u0 + 2*j0*j2*t0*t0*u0 - i2*n0*t0*t0*u0 - 
  i0*n2*t0*t0*u0 + 2*k0*m0*n0*t2*u0 - 2*j0*m0*o0*t2*u0 - 
  2*j0*k0*q0*t2*u0 + 2*i0*o0*q0*t2*u0 + 2*j0*j0*t0*t2*u0 - 
  2*i0*n0*t0*t2*u0 + m0*m0*o0*o0*u2 - 2*k0*m0*o0*q0*u2 + 
  k0*k0*q0*q0*u2 - m0*m0*n0*r0*u2 + 2*j0*m0*q0*r0*u2 - 
  i0*q0*q0*r0*u2 + 2*k0*m0*n0*t0*u2 - 2*j0*m0*o0*t0*u2 - 
  2*j0*k0*q0*t0*u2 + 2*i0*o0*q0*t0*u2 + j0*j0*t0*t0*u2 - 
  i0*n0*t0*t0*u2 - 2*l2*m0*o0*o0*v0 - 2*l0*m2*o0*o0*v0 - 
  4*l0*m0*o0*o2*v0 + 2*k2*m0*o0*p0*v0 + 2*k0*m2*o0*p0*v0 + 
  2*k0*m0*o2*p0*v0 + 2*k0*m0*o0*p2*v0 + 2*k2*l0*o0*q0*v0 + 
  2*k0*l2*o0*q0*v0 + 2*k0*l0*o2*q0*v0 - 4*k0*k2*p0*q0*v0 - 
  2*k0*k0*p2*q0*v0 + 2*k0*l0*o0*q2*v0 - 2*k0*k0*p0*q2*v0 + 
  2*l2*m0*n0*r0*v0 + 2*l0*m2*n0*r0*v0 + 2*l0*m0*n2*r0*v0 - 
  2*j2*m0*p0*r0*v0 - 2*j0*m2*p0*r0*v0 - 2*j0*m0*p2*r0*v0 - 
  2*j2*l0*q0*r0*v0 - 2*j0*l2*q0*r0*v0 + 2*i2*p0*q0*r0*v0 + 
  2*i0*p2*q0*r0*v0 - 2*j0*l0*q2*r0*v0 + 2*i0*p0*q2*r0*v0 + 
  2*l0*m0*n0*r2*v0 - 2*j0*m0*p0*r2*v0 - 2*j0*l0*q0*r2*v0 + 
  2*i0*p0*q0*r2*v0 - 2*k2*m0*n0*s0*v0 - 2*k0*m2*n0*s0*v0 - 
  2*k0*m0*n2*s0*v0 + 2*j2*m0*o0*s0*v0 + 2*j0*m2*o0*s0*v0 + 
  2*j0*m0*o2*s0*v0 + 2*j2*k0*q0*s0*v0 + 2*j0*k2*q0*s0*v0 - 
  2*i2*o0*q0*s0*v0 - 2*i0*o2*q0*s0*v0 + 2*j0*k0*q2*s0*v0 - 
  2*i0*o0*q2*s0*v0 - 2*k0*m0*n0*s2*v0 + 2*j0*m0*o0*s2*v0 + 
  2*j0*k0*q0*s2*v0 - 2*i0*o0*q0*s2*v0 - 2*k2*l0*n0*t0*v0 - 
  2*k0*l2*n0*t0*v0 - 2*k0*l0*n2*t0*v0 + 2*j2*l0*o0*t0*v0 + 
  2*j0*l2*o0*t0*v0 + 2*j0*l0*o2*t0*v0 + 2*j2*k0*p0*t0*v0 + 
  2*j0*k2*p0*t0*v0 - 2*i2*o0*p0*t0*v0 - 2*i0*o2*p0*t0*v0 + 
  2*j0*k0*p2*t0*v0 - 2*i0*o0*p2*t0*v0 - 4*j0*j2*s0*t0*v0 + 
  2*i2*n0*s0*t0*v0 + 2*i0*n2*s0*t0*v0 - 2*j0*j0*s2*t0*v0 + 
  2*i0*n0*s2*t0*v0 - 2*k0*l0*n0*t2*v0 + 2*j0*l0*o0*t2*v0 + 
  2*j0*k0*p0*t2*v0 - 2*i0*o0*p0*t2*v0 - 2*j0*j0*s0*t2*v0 + 
  2*i0*n0*s0*t2*v0 + 2*k0*k2*n0*v0*v0 + k0*k0*n2*v0*v0 - 
  2*j2*k0*o0*v0*v0 - 2*j0*k2*o0*v0*v0 + i2*o0*o0*v0*v0 - 
  2*j0*k0*o2*v0*v0 + 2*i0*o0*o2*v0*v0 + 2*j0*j2*r0*v0*v0 - 
  i2*n0*r0*v0*v0 - i0*n2*r0*v0*v0 + j0*j0*r2*v0*v0 - 
  i0*n0*r2*v0*v0 - 2*l0*m0*o0*o0*v2 + 2*k0*m0*o0*p0*v2 + 
  2*k0*l0*o0*q0*v2 - 2*k0*k0*p0*q0*v2 + 2*l0*m0*n0*r0*v2 - 
  2*j0*m0*p0*r0*v2 - 2*j0*l0*q0*r0*v2 + 2*i0*p0*q0*r0*v2 - 
  2*k0*m0*n0*s0*v2 + 2*j0*m0*o0*s0*v2 + 2*j0*k0*q0*s0*v2 - 
  2*i0*o0*q0*s0*v2 - 2*k0*l0*n0*t0*v2 + 2*j0*l0*o0*t0*v2 + 
  2*j0*k0*p0*t0*v2 - 2*i0*o0*p0*t0*v2 - 2*j0*j0*s0*t0*v2 + 
  2*i0*n0*s0*t0*v2 + 2*k0*k0*n0*v0*v2 - 4*j0*k0*o0*v0*v2 + 
  2*i0*o0*o0*v0*v2 + 2*j0*j0*r0*v0*v2 - 2*i0*n0*r0*v0*v2 + 
  2*l0*l2*o0*o0*w0 + 2*l0*l0*o0*o2*w0 - 2*k2*l0*o0*p0*w0 - 
  2*k0*l2*o0*p0*w0 - 2*k0*l0*o2*p0*w0 + 2*k0*k2*p0*p0*w0 - 
  2*k0*l0*o0*p2*w0 + 2*k0*k0*p0*p2*w0 - 2*l0*l2*n0*r0*w0 - 
  l0*l0*n2*r0*w0 + 2*j2*l0*p0*r0*w0 + 2*j0*l2*p0*r0*w0 - 
  i2*p0*p0*r0*w0 + 2*j0*l0*p2*r0*w0 - 2*i0*p0*p2*r0*w0 - 
  l0*l0*n0*r2*w0 + 2*j0*l0*p0*r2*w0 - i0*p0*p0*r2*w0 + 
  2*k2*l0*n0*s0*w0 + 2*k0*l2*n0*s0*w0 + 2*k0*l0*n2*s0*w0 - 
  2*j2*l0*o0*s0*w0 - 2*j0*l2*o0*s0*w0 - 2*j0*l0*o2*s0*w0 - 
  2*j2*k0*p0*s0*w0 - 2*j0*k2*p0*s0*w0 + 2*i2*o0*p0*s0*w0 + 
  2*i0*o2*p0*s0*w0 - 2*j0*k0*p2*s0*w0 + 2*i0*o0*p2*s0*w0 + 
  2*j0*j2*s0*s0*w0 - i2*n0*s0*s0*w0 - i0*n2*s0*s0*w0 + 
  2*k0*l0*n0*s2*w0 - 2*j0*l0*o0*s2*w0 - 2*j0*k0*p0*s2*w0 + 
  2*i0*o0*p0*s2*w0 + 2*j0*j0*s0*s2*w0 - 2*i0*n0*s0*s2*w0 - 
  2*k0*k2*n0*u0*w0 - k0*k0*n2*u0*w0 + 2*j2*k0*o0*u0*w0 + 
  2*j0*k2*o0*u0*w0 - i2*o0*o0*u0*w0 + 2*j0*k0*o2*u0*w0 - 
  2*i0*o0*o2*u0*w0 - 2*j0*j2*r0*u0*w0 + i2*n0*r0*u0*w0 + 
  i0*n2*r0*u0*w0 - j0*j0*r2*u0*w0 + i0*n0*r2*u0*w0 - 
  k0*k0*n0*u2*w0 + 2*j0*k0*o0*u2*w0 - i0*o0*o0*u2*w0 - 
  j0*j0*r0*u2*w0 + i0*n0*r0*u2*w0 + l0*l0*o0*o0*w2 - 
  2*k0*l0*o0*p0*w2 + k0*k0*p0*p0*w2 - l0*l0*n0*r0*w2 + 
  2*j0*l0*p0*r0*w2 - i0*p0*p0*r0*w2 + 2*k0*l0*n0*s0*w2 - 
  2*j0*l0*o0*s0*w2 - 2*j0*k0*p0*s0*w2 + 2*i0*o0*p0*s0*w2 + 
  j0*j0*s0*s0*w2 - i0*n0*s0*s0*w2 - k0*k0*n0*u0*w2 + 
  2*j0*k0*o0*u0*w2 - i0*o0*o0*u0*w2 - j0*j0*r0*u0*w2 + 
  i0*n0*r0*u0*w2;
    
    k[11] = 2*m1*m2*p0*p0*r0 + 
  4*m0*m2*p0*p1*r0 + 4*m0*m1*p0*p2*r0 + 2*m0*m0*p1*p2*r0 - 
  2*l2*m1*p0*q0*r0 - 2*l1*m2*p0*q0*r0 - 2*l2*m0*p1*q0*r0 - 
  2*l0*m2*p1*q0*r0 - 2*l1*m0*p2*q0*r0 - 2*l0*m1*p2*q0*r0 + 
  2*l1*l2*q0*q0*r0 - 2*l2*m0*p0*q1*r0 - 2*l0*m2*p0*q1*r0 - 
  2*l0*m0*p2*q1*r0 + 4*l0*l2*q0*q1*r0 - 2*l1*m0*p0*q2*r0 - 
  2*l0*m1*p0*q2*r0 - 2*l0*m0*p1*q2*r0 + 4*l0*l1*q0*q2*r0 + 
  2*l0*l0*q1*q2*r0 + 2*m0*m2*p0*p0*r1 + 2*m0*m0*p0*p2*r1 - 
  2*l2*m0*p0*q0*r1 - 2*l0*m2*p0*q0*r1 - 2*l0*m0*p2*q0*r1 + 
  2*l0*l2*q0*q0*r1 - 2*l0*m0*p0*q2*r1 + 2*l0*l0*q0*q2*r1 + 
  2*m0*m1*p0*p0*r2 + 2*m0*m0*p0*p1*r2 - 2*l1*m0*p0*q0*r2 - 
  2*l0*m1*p0*q0*r2 - 2*l0*m0*p1*q0*r2 + 2*l0*l1*q0*q0*r2 - 
  2*l0*m0*p0*q1*r2 + 2*l0*l0*q0*q1*r2 - 4*m1*m2*o0*p0*s0 - 
  4*m0*m2*o1*p0*s0 - 4*m0*m1*o2*p0*s0 - 4*m0*m2*o0*p1*s0 - 
  2*m0*m0*o2*p1*s0 - 4*m0*m1*o0*p2*s0 - 2*m0*m0*o1*p2*s0 + 
  2*l2*m1*o0*q0*s0 + 2*l1*m2*o0*q0*s0 + 2*l2*m0*o1*q0*s0 + 
  2*l0*m2*o1*q0*s0 + 2*l1*m0*o2*q0*s0 + 2*l0*m1*o2*q0*s0 + 
  2*k2*m1*p0*q0*s0 + 2*k1*m2*p0*q0*s0 + 2*k2*m0*p1*q0*s0 + 
  2*k0*m2*p1*q0*s0 + 2*k1*m0*p2*q0*s0 + 2*k0*m1*p2*q0*s0 - 
  2*k2*l1*q0*q0*s0 - 2*k1*l2*q0*q0*s0 + 2*l2*m0*o0*q1*s0 + 
  2*l0*m2*o0*q1*s0 + 2*l0*m0*o2*q1*s0 + 2*k2*m0*p0*q1*s0 + 
  2*k0*m2*p0*q1*s0 + 2*k0*m0*p2*q1*s0 - 4*k2*l0*q0*q1*s0 - 
  4*k0*l2*q0*q1*s0 + 2*l1*m0*o0*q2*s0 + 2*l0*m1*o0*q2*s0 + 
  2*l0*m0*o1*q2*s0 + 2*k1*m0*p0*q2*s0 + 2*k0*m1*p0*q2*s0 + 
  2*k0*m0*p1*q2*s0 - 4*k1*l0*q0*q2*s0 - 4*k0*l1*q0*q2*s0 - 
  4*k0*l0*q1*q2*s0 + 2*m1*m2*n0*s0*s0 + 2*m0*m2*n1*s0*s0 + 
  2*m0*m1*n2*s0*s0 - 2*j2*m1*q0*s0*s0 - 2*j1*m2*q0*s0*s0 - 
  2*j2*m0*q1*s0*s0 - 2*j0*m2*q1*s0*s0 + 2*i2*q0*q1*s0*s0 - 
  2*j1*m0*q2*s0*s0 - 2*j0*m1*q2*s0*s0 + 2*i1*q0*q2*s0*s0 + 
  2*i0*q1*q2*s0*s0 - 4*m0*m2*o0*p0*s1 - 2*m0*m0*o2*p0*s1 - 
  2*m0*m0*o0*p2*s1 + 2*l2*m0*o0*q0*s1 + 2*l0*m2*o0*q0*s1 + 
  2*l0*m0*o2*q0*s1 + 2*k2*m0*p0*q0*s1 + 2*k0*m2*p0*q0*s1 + 
  2*k0*m0*p2*q0*s1 - 2*k2*l0*q0*q0*s1 - 2*k0*l2*q0*q0*s1 + 
  2*l0*m0*o0*q2*s1 + 2*k0*m0*p0*q2*s1 - 4*k0*l0*q0*q2*s1 + 
  4*m0*m2*n0*s0*s1 + 2*m0*m0*n2*s0*s1 - 4*j2*m0*q0*s0*s1 - 
  4*j0*m2*q0*s0*s1 + 2*i2*q0*q0*s0*s1 - 4*j0*m0*q2*s0*s1 + 
  4*i0*q0*q2*s0*s1 - 4*m0*m1*o0*p0*s2 - 2*m0*m0*o1*p0*s2 - 
  2*m0*m0*o0*p1*s2 + 2*l1*m0*o0*q0*s2 + 2*l0*m1*o0*q0*s2 + 
  2*l0*m0*o1*q0*s2 + 2*k1*m0*p0*q0*s2 + 2*k0*m1*p0*q0*s2 + 
  2*k0*m0*p1*q0*s2 - 2*k1*l0*q0*q0*s2 - 2*k0*l1*q0*q0*s2 + 
  2*l0*m0*o0*q1*s2 + 2*k0*m0*p0*q1*s2 - 4*k0*l0*q0*q1*s2 + 
  4*m0*m1*n0*s0*s2 + 2*m0*m0*n1*s0*s2 - 4*j1*m0*q0*s0*s2 - 
  4*j0*m1*q0*s0*s2 + 2*i1*q0*q0*s0*s2 - 4*j0*m0*q1*s0*s2 + 
  4*i0*q0*q1*s0*s2 + 2*m0*m0*n0*s1*s2 - 4*j0*m0*q0*s1*s2 + 
  2*i0*q0*q0*s1*s2 + 2*l2*m1*o0*p0*t0 + 2*l1*m2*o0*p0*t0 + 
  2*l2*m0*o1*p0*t0 + 2*l0*m2*o1*p0*t0 + 2*l1*m0*o2*p0*t0 + 
  2*l0*m1*o2*p0*t0 - 2*k2*m1*p0*p0*t0 - 2*k1*m2*p0*p0*t0 + 
  2*l2*m0*o0*p1*t0 + 2*l0*m2*o0*p1*t0 + 2*l0*m0*o2*p1*t0 - 
  4*k2*m0*p0*p1*t0 - 4*k0*m2*p0*p1*t0 + 2*l1*m0*o0*p2*t0 + 
  2*l0*m1*o0*p2*t0 + 2*l0*m0*o1*p2*t0 - 4*k1*m0*p0*p2*t0 - 
  4*k0*m1*p0*p2*t0 - 4*k0*m0*p1*p2*t0 - 4*l1*l2*o0*q0*t0 - 
  4*l0*l2*o1*q0*t0 - 4*l0*l1*o2*q0*t0 + 2*k2*l1*p0*q0*t0 + 
  2*k1*l2*p0*q0*t0 + 2*k2*l0*p1*q0*t0 + 2*k0*l2*p1*q0*t0 + 
  2*k1*l0*p2*q0*t0 + 2*k0*l1*p2*q0*t0 - 4*l0*l2*o0*q1*t0 - 
  2*l0*l0*o2*q1*t0 + 2*k2*l0*p0*q1*t0 + 2*k0*l2*p0*q1*t0 + 
  2*k0*l0*p2*q1*t0 - 4*l0*l1*o0*q2*t0 - 2*l0*l0*o1*q2*t0 + 
  2*k1*l0*p0*q2*t0 + 2*k0*l1*p0*q2*t0 + 2*k0*l0*p1*q2*t0 - 
  2*l2*m1*n0*s0*t0 - 2*l1*m2*n0*s0*t0 - 2*l2*m0*n1*s0*t0 - 
  2*l0*m2*n1*s0*t0 - 2*l1*m0*n2*s0*t0 - 2*l0*m1*n2*s0*t0 + 
  2*j2*m1*p0*s0*t0 + 2*j1*m2*p0*s0*t0 + 2*j2*m0*p1*s0*t0 + 
  2*j0*m2*p1*s0*t0 + 2*j1*m0*p2*s0*t0 + 2*j0*m1*p2*s0*t0 + 
  2*j2*l1*q0*s0*t0 + 2*j1*l2*q0*s0*t0 - 2*i2*p1*q0*s0*t0 - 
  2*i1*p2*q0*s0*t0 + 2*j2*l0*q1*s0*t0 + 2*j0*l2*q1*s0*t0 - 
  2*i2*p0*q1*s0*t0 - 2*i0*p2*q1*s0*t0 + 2*j1*l0*q2*s0*t0 + 
  2*j0*l1*q2*s0*t0 - 2*i1*p0*q2*s0*t0 - 2*i0*p1*q2*s0*t0 - 
  2*l2*m0*n0*s1*t0 - 2*l0*m2*n0*s1*t0 - 2*l0*m0*n2*s1*t0 + 
  2*j2*m0*p0*s1*t0 + 2*j0*m2*p0*s1*t0 + 2*j0*m0*p2*s1*t0 + 
  2*j2*l0*q0*s1*t0 + 2*j0*l2*q0*s1*t0 - 2*i2*p0*q0*s1*t0 - 
  2*i0*p2*q0*s1*t0 + 2*j0*l0*q2*s1*t0 - 2*i0*p0*q2*s1*t0 - 
  2*l1*m0*n0*s2*t0 - 2*l0*m1*n0*s2*t0 - 2*l0*m0*n1*s2*t0 + 
  2*j1*m0*p0*s2*t0 + 2*j0*m1*p0*s2*t0 + 2*j0*m0*p1*s2*t0 + 
  2*j1*l0*q0*s2*t0 + 2*j0*l1*q0*s2*t0 - 2*i1*p0*q0*s2*t0 - 
  2*i0*p1*q0*s2*t0 + 2*j0*l0*q1*s2*t0 - 2*i0*p0*q1*s2*t0 + 
  2*l1*l2*n0*t0*t0 + 2*l0*l2*n1*t0*t0 + 2*l0*l1*n2*t0*t0 - 
  2*j2*l1*p0*t0*t0 - 2*j1*l2*p0*t0*t0 - 2*j2*l0*p1*t0*t0 - 
  2*j0*l2*p1*t0*t0 + 2*i2*p0*p1*t0*t0 - 2*j1*l0*p2*t0*t0 - 
  2*j0*l1*p2*t0*t0 + 2*i1*p0*p2*t0*t0 + 2*i0*p1*p2*t0*t0 + 
  2*l2*m0*o0*p0*t1 + 2*l0*m2*o0*p0*t1 + 2*l0*m0*o2*p0*t1 - 
  2*k2*m0*p0*p0*t1 - 2*k0*m2*p0*p0*t1 + 2*l0*m0*o0*p2*t1 - 
  4*k0*m0*p0*p2*t1 - 4*l0*l2*o0*q0*t1 - 2*l0*l0*o2*q0*t1 + 
  2*k2*l0*p0*q0*t1 + 2*k0*l2*p0*q0*t1 + 2*k0*l0*p2*q0*t1 - 
  2*l0*l0*o0*q2*t1 + 2*k0*l0*p0*q2*t1 - 2*l2*m0*n0*s0*t1 - 
  2*l0*m2*n0*s0*t1 - 2*l0*m0*n2*s0*t1 + 2*j2*m0*p0*s0*t1 + 
  2*j0*m2*p0*s0*t1 + 2*j0*m0*p2*s0*t1 + 2*j2*l0*q0*s0*t1 + 
  2*j0*l2*q0*s0*t1 - 2*i2*p0*q0*s0*t1 - 2*i0*p2*q0*s0*t1 + 
  2*j0*l0*q2*s0*t1 - 2*i0*p0*q2*s0*t1 - 2*l0*m0*n0*s2*t1 + 
  2*j0*m0*p0*s2*t1 + 2*j0*l0*q0*s2*t1 - 2*i0*p0*q0*s2*t1 + 
  4*l0*l2*n0*t0*t1 + 2*l0*l0*n2*t0*t1 - 4*j2*l0*p0*t0*t1 - 
  4*j0*l2*p0*t0*t1 + 2*i2*p0*p0*t0*t1 - 4*j0*l0*p2*t0*t1 + 
  4*i0*p0*p2*t0*t1 + 2*l1*m0*o0*p0*t2 + 2*l0*m1*o0*p0*t2 + 
  2*l0*m0*o1*p0*t2 - 2*k1*m0*p0*p0*t2 - 2*k0*m1*p0*p0*t2 + 
  2*l0*m0*o0*p1*t2 - 4*k0*m0*p0*p1*t2 - 4*l0*l1*o0*q0*t2 - 
  2*l0*l0*o1*q0*t2 + 2*k1*l0*p0*q0*t2 + 2*k0*l1*p0*q0*t2 + 
  2*k0*l0*p1*q0*t2 - 2*l0*l0*o0*q1*t2 + 2*k0*l0*p0*q1*t2 - 
  2*l1*m0*n0*s0*t2 - 2*l0*m1*n0*s0*t2 - 2*l0*m0*n1*s0*t2 + 
  2*j1*m0*p0*s0*t2 + 2*j0*m1*p0*s0*t2 + 2*j0*m0*p1*s0*t2 + 
  2*j1*l0*q0*s0*t2 + 2*j0*l1*q0*s0*t2 - 2*i1*p0*q0*s0*t2 - 
  2*i0*p1*q0*s0*t2 + 2*j0*l0*q1*s0*t2 - 2*i0*p0*q1*s0*t2 - 
  2*l0*m0*n0*s1*t2 + 2*j0*m0*p0*s1*t2 + 2*j0*l0*q0*s1*t2 - 
  2*i0*p0*q0*s1*t2 + 4*l0*l1*n0*t0*t2 + 2*l0*l0*n1*t0*t2 - 
  4*j1*l0*p0*t0*t2 - 4*j0*l1*p0*t0*t2 + 2*i1*p0*p0*t0*t2 - 
  4*j0*l0*p1*t0*t2 + 4*i0*p0*p1*t0*t2 + 2*l0*l0*n0*t1*t2 - 
  4*j0*l0*p0*t1*t2 + 2*i0*p0*p0*t1*t2 + 2*m1*m2*o0*o0*u0 + 
  4*m0*m2*o0*o1*u0 + 4*m0*m1*o0*o2*u0 + 2*m0*m0*o1*o2*u0 - 
  2*k2*m1*o0*q0*u0 - 2*k1*m2*o0*q0*u0 - 2*k2*m0*o1*q0*u0 - 
  2*k0*m2*o1*q0*u0 - 2*k1*m0*o2*q0*u0 - 2*k0*m1*o2*q0*u0 + 
  2*k1*k2*q0*q0*u0 - 2*k2*m0*o0*q1*u0 - 2*k0*m2*o0*q1*u0 - 
  2*k0*m0*o2*q1*u0 + 4*k0*k2*q0*q1*u0 - 2*k1*m0*o0*q2*u0 - 
  2*k0*m1*o0*q2*u0 - 2*k0*m0*o1*q2*u0 + 4*k0*k1*q0*q2*u0 + 
  2*k0*k0*q1*q2*u0 - 2*m1*m2*n0*r0*u0 - 2*m0*m2*n1*r0*u0 - 
  2*m0*m1*n2*r0*u0 + 2*j2*m1*q0*r0*u0 + 2*j1*m2*q0*r0*u0 + 
  2*j2*m0*q1*r0*u0 + 2*j0*m2*q1*r0*u0 - 2*i2*q0*q1*r0*u0 + 
  2*j1*m0*q2*r0*u0 + 2*j0*m1*q2*r0*u0 - 2*i1*q0*q2*r0*u0 - 
  2*i0*q1*q2*r0*u0 - 2*m0*m2*n0*r1*u0 - m0*m0*n2*r1*u0 + 
  2*j2*m0*q0*r1*u0 + 2*j0*m2*q0*r1*u0 - i2*q0*q0*r1*u0 + 
  2*j0*m0*q2*r1*u0 - 2*i0*q0*q2*r1*u0 - 2*m0*m1*n0*r2*u0 - 
  m0*m0*n1*r2*u0 + 2*j1*m0*q0*r2*u0 + 2*j0*m1*q0*r2*u0 - 
  i1*q0*q0*r2*u0 + 2*j0*m0*q1*r2*u0 - 2*i0*q0*q1*r2*u0 + 
  2*k2*m1*n0*t0*u0 + 2*k1*m2*n0*t0*u0 + 2*k2*m0*n1*t0*u0 + 
  2*k0*m2*n1*t0*u0 + 2*k1*m0*n2*t0*u0 + 2*k0*m1*n2*t0*u0 - 
  2*j2*m1*o0*t0*u0 - 2*j1*m2*o0*t0*u0 - 2*j2*m0*o1*t0*u0 - 
  2*j0*m2*o1*t0*u0 - 2*j1*m0*o2*t0*u0 - 2*j0*m1*o2*t0*u0 - 
  2*j2*k1*q0*t0*u0 - 2*j1*k2*q0*t0*u0 + 2*i2*o1*q0*t0*u0 + 
  2*i1*o2*q0*t0*u0 - 2*j2*k0*q1*t0*u0 - 2*j0*k2*q1*t0*u0 + 
  2*i2*o0*q1*t0*u0 + 2*i0*o2*q1*t0*u0 - 2*j1*k0*q2*t0*u0 - 
  2*j0*k1*q2*t0*u0 + 2*i1*o0*q2*t0*u0 + 2*i0*o1*q2*t0*u0 + 
  2*j1*j2*t0*t0*u0 - i2*n1*t0*t0*u0 - i1*n2*t0*t0*u0 + 
  2*k2*m0*n0*t1*u0 + 2*k0*m2*n0*t1*u0 + 2*k0*m0*n2*t1*u0 - 
  2*j2*m0*o0*t1*u0 - 2*j0*m2*o0*t1*u0 - 2*j0*m0*o2*t1*u0 - 
  2*j2*k0*q0*t1*u0 - 2*j0*k2*q0*t1*u0 + 2*i2*o0*q0*t1*u0 + 
  2*i0*o2*q0*t1*u0 - 2*j0*k0*q2*t1*u0 + 2*i0*o0*q2*t1*u0 + 
  4*j0*j2*t0*t1*u0 - 2*i2*n0*t0*t1*u0 - 2*i0*n2*t0*t1*u0 + 
  2*k1*m0*n0*t2*u0 + 2*k0*m1*n0*t2*u0 + 2*k0*m0*n1*t2*u0 - 
  2*j1*m0*o0*t2*u0 - 2*j0*m1*o0*t2*u0 - 2*j0*m0*o1*t2*u0 - 
  2*j1*k0*q0*t2*u0 - 2*j0*k1*q0*t2*u0 + 2*i1*o0*q0*t2*u0 + 
  2*i0*o1*q0*t2*u0 - 2*j0*k0*q1*t2*u0 + 2*i0*o0*q1*t2*u0 + 
  4*j0*j1*t0*t2*u0 - 2*i1*n0*t0*t2*u0 - 2*i0*n1*t0*t2*u0 + 
  2*j0*j0*t1*t2*u0 - 2*i0*n0*t1*t2*u0 + 2*m0*m2*o0*o0*u1 + 
  2*m0*m0*o0*o2*u1 - 2*k2*m0*o0*q0*u1 - 2*k0*m2*o0*q0*u1 - 
  2*k0*m0*o2*q0*u1 + 2*k0*k2*q0*q0*u1 - 2*k0*m0*o0*q2*u1 + 
  2*k0*k0*q0*q2*u1 - 2*m0*m2*n0*r0*u1 - m0*m0*n2*r0*u1 + 
  2*j2*m0*q0*r0*u1 + 2*j0*m2*q0*r0*u1 - i2*q0*q0*r0*u1 + 
  2*j0*m0*q2*r0*u1 - 2*i0*q0*q2*r0*u1 - m0*m0*n0*r2*u1 + 
  2*j0*m0*q0*r2*u1 - i0*q0*q0*r2*u1 + 2*k2*m0*n0*t0*u1 + 
  2*k0*m2*n0*t0*u1 + 2*k0*m0*n2*t0*u1 - 2*j2*m0*o0*t0*u1 - 
  2*j0*m2*o0*t0*u1 - 2*j0*m0*o2*t0*u1 - 2*j2*k0*q0*t0*u1 - 
  2*j0*k2*q0*t0*u1 + 2*i2*o0*q0*t0*u1 + 2*i0*o2*q0*t0*u1 - 
  2*j0*k0*q2*t0*u1 + 2*i0*o0*q2*t0*u1 + 2*j0*j2*t0*t0*u1 - 
  i2*n0*t0*t0*u1 - i0*n2*t0*t0*u1 + 2*k0*m0*n0*t2*u1 - 
  2*j0*m0*o0*t2*u1 - 2*j0*k0*q0*t2*u1 + 2*i0*o0*q0*t2*u1 + 
  2*j0*j0*t0*t2*u1 - 2*i0*n0*t0*t2*u1 + 2*m0*m1*o0*o0*u2 + 
  2*m0*m0*o0*o1*u2 - 2*k1*m0*o0*q0*u2 - 2*k0*m1*o0*q0*u2 - 
  2*k0*m0*o1*q0*u2 + 2*k0*k1*q0*q0*u2 - 2*k0*m0*o0*q1*u2 + 
  2*k0*k0*q0*q1*u2 - 2*m0*m1*n0*r0*u2 - m0*m0*n1*r0*u2 + 
  2*j1*m0*q0*r0*u2 + 2*j0*m1*q0*r0*u2 - i1*q0*q0*r0*u2 + 
  2*j0*m0*q1*r0*u2 - 2*i0*q0*q1*r0*u2 - m0*m0*n0*r1*u2 + 
  2*j0*m0*q0*r1*u2 - i0*q0*q0*r1*u2 + 2*k1*m0*n0*t0*u2 + 
  2*k0*m1*n0*t0*u2 + 2*k0*m0*n1*t0*u2 - 2*j1*m0*o0*t0*u2 - 
  2*j0*m1*o0*t0*u2 - 2*j0*m0*o1*t0*u2 - 2*j1*k0*q0*t0*u2 - 
  2*j0*k1*q0*t0*u2 + 2*i1*o0*q0*t0*u2 + 2*i0*o1*q0*t0*u2 - 
  2*j0*k0*q1*t0*u2 + 2*i0*o0*q1*t0*u2 + 2*j0*j1*t0*t0*u2 - 
  i1*n0*t0*t0*u2 - i0*n1*t0*t0*u2 + 2*k0*m0*n0*t1*u2 - 
  2*j0*m0*o0*t1*u2 - 2*j0*k0*q0*t1*u2 + 2*i0*o0*q0*t1*u2 + 
  2*j0*j0*t0*t1*u2 - 2*i0*n0*t0*t1*u2 - 2*l2*m1*o0*o0*v0 - 
  2*l1*m2*o0*o0*v0 - 4*l2*m0*o0*o1*v0 - 4*l0*m2*o0*o1*v0 - 
  4*l1*m0*o0*o2*v0 - 4*l0*m1*o0*o2*v0 - 4*l0*m0*o1*o2*v0 + 
  2*k2*m1*o0*p0*v0 + 2*k1*m2*o0*p0*v0 + 2*k2*m0*o1*p0*v0 + 
  2*k0*m2*o1*p0*v0 + 2*k1*m0*o2*p0*v0 + 2*k0*m1*o2*p0*v0 + 
  2*k2*m0*o0*p1*v0 + 2*k0*m2*o0*p1*v0 + 2*k0*m0*o2*p1*v0 + 
  2*k1*m0*o0*p2*v0 + 2*k0*m1*o0*p2*v0 + 2*k0*m0*o1*p2*v0 + 
  2*k2*l1*o0*q0*v0 + 2*k1*l2*o0*q0*v0 + 2*k2*l0*o1*q0*v0 + 
  2*k0*l2*o1*q0*v0 + 2*k1*l0*o2*q0*v0 + 2*k0*l1*o2*q0*v0 - 
  4*k1*k2*p0*q0*v0 - 4*k0*k2*p1*q0*v0 - 4*k0*k1*p2*q0*v0 + 
  2*k2*l0*o0*q1*v0 + 2*k0*l2*o0*q1*v0 + 2*k0*l0*o2*q1*v0 - 
  4*k0*k2*p0*q1*v0 - 2*k0*k0*p2*q1*v0 + 2*k1*l0*o0*q2*v0 + 
  2*k0*l1*o0*q2*v0 + 2*k0*l0*o1*q2*v0 - 4*k0*k1*p0*q2*v0 - 
  2*k0*k0*p1*q2*v0 + 2*l2*m1*n0*r0*v0 + 2*l1*m2*n0*r0*v0 + 
  2*l2*m0*n1*r0*v0 + 2*l0*m2*n1*r0*v0 + 2*l1*m0*n2*r0*v0 + 
  2*l0*m1*n2*r0*v0 - 2*j2*m1*p0*r0*v0 - 2*j1*m2*p0*r0*v0 - 
  2*j2*m0*p1*r0*v0 - 2*j0*m2*p1*r0*v0 - 2*j1*m0*p2*r0*v0 - 
  2*j0*m1*p2*r0*v0 - 2*j2*l1*q0*r0*v0 - 2*j1*l2*q0*r0*v0 + 
  2*i2*p1*q0*r0*v0 + 2*i1*p2*q0*r0*v0 - 2*j2*l0*q1*r0*v0 - 
  2*j0*l2*q1*r0*v0 + 2*i2*p0*q1*r0*v0 + 2*i0*p2*q1*r0*v0 - 
  2*j1*l0*q2*r0*v0 - 2*j0*l1*q2*r0*v0 + 2*i1*p0*q2*r0*v0 + 
  2*i0*p1*q2*r0*v0 + 2*l2*m0*n0*r1*v0 + 2*l0*m2*n0*r1*v0 + 
  2*l0*m0*n2*r1*v0 - 2*j2*m0*p0*r1*v0 - 2*j0*m2*p0*r1*v0 - 
  2*j0*m0*p2*r1*v0 - 2*j2*l0*q0*r1*v0 - 2*j0*l2*q0*r1*v0 + 
  2*i2*p0*q0*r1*v0 + 2*i0*p2*q0*r1*v0 - 2*j0*l0*q2*r1*v0 + 
  2*i0*p0*q2*r1*v0 + 2*l1*m0*n0*r2*v0 + 2*l0*m1*n0*r2*v0 + 
  2*l0*m0*n1*r2*v0 - 2*j1*m0*p0*r2*v0 - 2*j0*m1*p0*r2*v0 - 
  2*j0*m0*p1*r2*v0 - 2*j1*l0*q0*r2*v0 - 2*j0*l1*q0*r2*v0 + 
  2*i1*p0*q0*r2*v0 + 2*i0*p1*q0*r2*v0 - 2*j0*l0*q1*r2*v0 + 
  2*i0*p0*q1*r2*v0 - 2*k2*m1*n0*s0*v0 - 2*k1*m2*n0*s0*v0 - 
  2*k2*m0*n1*s0*v0 - 2*k0*m2*n1*s0*v0 - 2*k1*m0*n2*s0*v0 - 
  2*k0*m1*n2*s0*v0 + 2*j2*m1*o0*s0*v0 + 2*j1*m2*o0*s0*v0 + 
  2*j2*m0*o1*s0*v0 + 2*j0*m2*o1*s0*v0 + 2*j1*m0*o2*s0*v0 + 
  2*j0*m1*o2*s0*v0 + 2*j2*k1*q0*s0*v0 + 2*j1*k2*q0*s0*v0 - 
  2*i2*o1*q0*s0*v0 - 2*i1*o2*q0*s0*v0 + 2*j2*k0*q1*s0*v0 + 
  2*j0*k2*q1*s0*v0 - 2*i2*o0*q1*s0*v0 - 2*i0*o2*q1*s0*v0 + 
  2*j1*k0*q2*s0*v0 + 2*j0*k1*q2*s0*v0 - 2*i1*o0*q2*s0*v0 - 
  2*i0*o1*q2*s0*v0 - 2*k2*m0*n0*s1*v0 - 2*k0*m2*n0*s1*v0 - 
  2*k0*m0*n2*s1*v0 + 2*j2*m0*o0*s1*v0 + 2*j0*m2*o0*s1*v0 + 
  2*j0*m0*o2*s1*v0 + 2*j2*k0*q0*s1*v0 + 2*j0*k2*q0*s1*v0 - 
  2*i2*o0*q0*s1*v0 - 2*i0*o2*q0*s1*v0 + 2*j0*k0*q2*s1*v0 - 
  2*i0*o0*q2*s1*v0 - 2*k1*m0*n0*s2*v0 - 2*k0*m1*n0*s2*v0 - 
  2*k0*m0*n1*s2*v0 + 2*j1*m0*o0*s2*v0 + 2*j0*m1*o0*s2*v0 + 
  2*j0*m0*o1*s2*v0 + 2*j1*k0*q0*s2*v0 + 2*j0*k1*q0*s2*v0 - 
  2*i1*o0*q0*s2*v0 - 2*i0*o1*q0*s2*v0 + 2*j0*k0*q1*s2*v0 - 
  2*i0*o0*q1*s2*v0 - 2*k2*l1*n0*t0*v0 - 2*k1*l2*n0*t0*v0 - 
  2*k2*l0*n1*t0*v0 - 2*k0*l2*n1*t0*v0 - 2*k1*l0*n2*t0*v0 - 
  2*k0*l1*n2*t0*v0 + 2*j2*l1*o0*t0*v0 + 2*j1*l2*o0*t0*v0 + 
  2*j2*l0*o1*t0*v0 + 2*j0*l2*o1*t0*v0 + 2*j1*l0*o2*t0*v0 + 
  2*j0*l1*o2*t0*v0 + 2*j2*k1*p0*t0*v0 + 2*j1*k2*p0*t0*v0 - 
  2*i2*o1*p0*t0*v0 - 2*i1*o2*p0*t0*v0 + 2*j2*k0*p1*t0*v0 + 
  2*j0*k2*p1*t0*v0 - 2*i2*o0*p1*t0*v0 - 2*i0*o2*p1*t0*v0 + 
  2*j1*k0*p2*t0*v0 + 2*j0*k1*p2*t0*v0 - 2*i1*o0*p2*t0*v0 - 
  2*i0*o1*p2*t0*v0 - 4*j1*j2*s0*t0*v0 + 2*i2*n1*s0*t0*v0 + 
  2*i1*n2*s0*t0*v0 - 4*j0*j2*s1*t0*v0 + 2*i2*n0*s1*t0*v0 + 
  2*i0*n2*s1*t0*v0 - 4*j0*j1*s2*t0*v0 + 2*i1*n0*s2*t0*v0 + 
  2*i0*n1*s2*t0*v0 - 2*k2*l0*n0*t1*v0 - 2*k0*l2*n0*t1*v0 - 
  2*k0*l0*n2*t1*v0 + 2*j2*l0*o0*t1*v0 + 2*j0*l2*o0*t1*v0 + 
  2*j0*l0*o2*t1*v0 + 2*j2*k0*p0*t1*v0 + 2*j0*k2*p0*t1*v0 - 
  2*i2*o0*p0*t1*v0 - 2*i0*o2*p0*t1*v0 + 2*j0*k0*p2*t1*v0 - 
  2*i0*o0*p2*t1*v0 - 4*j0*j2*s0*t1*v0 + 2*i2*n0*s0*t1*v0 + 
  2*i0*n2*s0*t1*v0 - 2*j0*j0*s2*t1*v0 + 2*i0*n0*s2*t1*v0 - 
  2*k1*l0*n0*t2*v0 - 2*k0*l1*n0*t2*v0 - 2*k0*l0*n1*t2*v0 + 
  2*j1*l0*o0*t2*v0 + 2*j0*l1*o0*t2*v0 + 2*j0*l0*o1*t2*v0 + 
  2*j1*k0*p0*t2*v0 + 2*j0*k1*p0*t2*v0 - 2*i1*o0*p0*t2*v0 - 
  2*i0*o1*p0*t2*v0 + 2*j0*k0*p1*t2*v0 - 2*i0*o0*p1*t2*v0 - 
  4*j0*j1*s0*t2*v0 + 2*i1*n0*s0*t2*v0 + 2*i0*n1*s0*t2*v0 - 
  2*j0*j0*s1*t2*v0 + 2*i0*n0*s1*t2*v0 + 2*k1*k2*n0*v0*v0 + 
  2*k0*k2*n1*v0*v0 + 2*k0*k1*n2*v0*v0 - 2*j2*k1*o0*v0*v0 - 
  2*j1*k2*o0*v0*v0 - 2*j2*k0*o1*v0*v0 - 2*j0*k2*o1*v0*v0 + 
  2*i2*o0*o1*v0*v0 - 2*j1*k0*o2*v0*v0 - 2*j0*k1*o2*v0*v0 + 
  2*i1*o0*o2*v0*v0 + 2*i0*o1*o2*v0*v0 + 2*j1*j2*r0*v0*v0 - 
  i2*n1*r0*v0*v0 - i1*n2*r0*v0*v0 + 2*j0*j2*r1*v0*v0 - 
  i2*n0*r1*v0*v0 - i0*n2*r1*v0*v0 + 2*j0*j1*r2*v0*v0 - 
  i1*n0*r2*v0*v0 - i0*n1*r2*v0*v0 - 2*l2*m0*o0*o0*v1 - 
  2*l0*m2*o0*o0*v1 - 4*l0*m0*o0*o2*v1 + 2*k2*m0*o0*p0*v1 + 
  2*k0*m2*o0*p0*v1 + 2*k0*m0*o2*p0*v1 + 2*k0*m0*o0*p2*v1 + 
  2*k2*l0*o0*q0*v1 + 2*k0*l2*o0*q0*v1 + 2*k0*l0*o2*q0*v1 - 
  4*k0*k2*p0*q0*v1 - 2*k0*k0*p2*q0*v1 + 2*k0*l0*o0*q2*v1 - 
  2*k0*k0*p0*q2*v1 + 2*l2*m0*n0*r0*v1 + 2*l0*m2*n0*r0*v1 + 
  2*l0*m0*n2*r0*v1 - 2*j2*m0*p0*r0*v1 - 2*j0*m2*p0*r0*v1 - 
  2*j0*m0*p2*r0*v1 - 2*j2*l0*q0*r0*v1 - 2*j0*l2*q0*r0*v1 + 
  2*i2*p0*q0*r0*v1 + 2*i0*p2*q0*r0*v1 - 2*j0*l0*q2*r0*v1 + 
  2*i0*p0*q2*r0*v1 + 2*l0*m0*n0*r2*v1 - 2*j0*m0*p0*r2*v1 - 
  2*j0*l0*q0*r2*v1 + 2*i0*p0*q0*r2*v1 - 2*k2*m0*n0*s0*v1 - 
  2*k0*m2*n0*s0*v1 - 2*k0*m0*n2*s0*v1 + 2*j2*m0*o0*s0*v1 + 
  2*j0*m2*o0*s0*v1 + 2*j0*m0*o2*s0*v1 + 2*j2*k0*q0*s0*v1 + 
  2*j0*k2*q0*s0*v1 - 2*i2*o0*q0*s0*v1 - 2*i0*o2*q0*s0*v1 + 
  2*j0*k0*q2*s0*v1 - 2*i0*o0*q2*s0*v1 - 2*k0*m0*n0*s2*v1 + 
  2*j0*m0*o0*s2*v1 + 2*j0*k0*q0*s2*v1 - 2*i0*o0*q0*s2*v1 - 
  2*k2*l0*n0*t0*v1 - 2*k0*l2*n0*t0*v1 - 2*k0*l0*n2*t0*v1 + 
  2*j2*l0*o0*t0*v1 + 2*j0*l2*o0*t0*v1 + 2*j0*l0*o2*t0*v1 + 
  2*j2*k0*p0*t0*v1 + 2*j0*k2*p0*t0*v1 - 2*i2*o0*p0*t0*v1 - 
  2*i0*o2*p0*t0*v1 + 2*j0*k0*p2*t0*v1 - 2*i0*o0*p2*t0*v1 - 
  4*j0*j2*s0*t0*v1 + 2*i2*n0*s0*t0*v1 + 2*i0*n2*s0*t0*v1 - 
  2*j0*j0*s2*t0*v1 + 2*i0*n0*s2*t0*v1 - 2*k0*l0*n0*t2*v1 + 
  2*j0*l0*o0*t2*v1 + 2*j0*k0*p0*t2*v1 - 2*i0*o0*p0*t2*v1 - 
  2*j0*j0*s0*t2*v1 + 2*i0*n0*s0*t2*v1 + 4*k0*k2*n0*v0*v1 + 
  2*k0*k0*n2*v0*v1 - 4*j2*k0*o0*v0*v1 - 4*j0*k2*o0*v0*v1 + 
  2*i2*o0*o0*v0*v1 - 4*j0*k0*o2*v0*v1 + 4*i0*o0*o2*v0*v1 + 
  4*j0*j2*r0*v0*v1 - 2*i2*n0*r0*v0*v1 - 2*i0*n2*r0*v0*v1 + 
  2*j0*j0*r2*v0*v1 - 2*i0*n0*r2*v0*v1 - 2*l1*m0*o0*o0*v2 - 
  2*l0*m1*o0*o0*v2 - 4*l0*m0*o0*o1*v2 + 2*k1*m0*o0*p0*v2 + 
  2*k0*m1*o0*p0*v2 + 2*k0*m0*o1*p0*v2 + 2*k0*m0*o0*p1*v2 + 
  2*k1*l0*o0*q0*v2 + 2*k0*l1*o0*q0*v2 + 2*k0*l0*o1*q0*v2 - 
  4*k0*k1*p0*q0*v2 - 2*k0*k0*p1*q0*v2 + 2*k0*l0*o0*q1*v2 - 
  2*k0*k0*p0*q1*v2 + 2*l1*m0*n0*r0*v2 + 2*l0*m1*n0*r0*v2 + 
  2*l0*m0*n1*r0*v2 - 2*j1*m0*p0*r0*v2 - 2*j0*m1*p0*r0*v2 - 
  2*j0*m0*p1*r0*v2 - 2*j1*l0*q0*r0*v2 - 2*j0*l1*q0*r0*v2 + 
  2*i1*p0*q0*r0*v2 + 2*i0*p1*q0*r0*v2 - 2*j0*l0*q1*r0*v2 + 
  2*i0*p0*q1*r0*v2 + 2*l0*m0*n0*r1*v2 - 2*j0*m0*p0*r1*v2 - 
  2*j0*l0*q0*r1*v2 + 2*i0*p0*q0*r1*v2 - 2*k1*m0*n0*s0*v2 - 
  2*k0*m1*n0*s0*v2 - 2*k0*m0*n1*s0*v2 + 2*j1*m0*o0*s0*v2 + 
  2*j0*m1*o0*s0*v2 + 2*j0*m0*o1*s0*v2 + 2*j1*k0*q0*s0*v2 + 
  2*j0*k1*q0*s0*v2 - 2*i1*o0*q0*s0*v2 - 2*i0*o1*q0*s0*v2 + 
  2*j0*k0*q1*s0*v2 - 2*i0*o0*q1*s0*v2 - 2*k0*m0*n0*s1*v2 + 
  2*j0*m0*o0*s1*v2 + 2*j0*k0*q0*s1*v2 - 2*i0*o0*q0*s1*v2 - 
  2*k1*l0*n0*t0*v2 - 2*k0*l1*n0*t0*v2 - 2*k0*l0*n1*t0*v2 + 
  2*j1*l0*o0*t0*v2 + 2*j0*l1*o0*t0*v2 + 2*j0*l0*o1*t0*v2 + 
  2*j1*k0*p0*t0*v2 + 2*j0*k1*p0*t0*v2 - 2*i1*o0*p0*t0*v2 - 
  2*i0*o1*p0*t0*v2 + 2*j0*k0*p1*t0*v2 - 2*i0*o0*p1*t0*v2 - 
  4*j0*j1*s0*t0*v2 + 2*i1*n0*s0*t0*v2 + 2*i0*n1*s0*t0*v2 - 
  2*j0*j0*s1*t0*v2 + 2*i0*n0*s1*t0*v2 - 2*k0*l0*n0*t1*v2 + 
  2*j0*l0*o0*t1*v2 + 2*j0*k0*p0*t1*v2 - 2*i0*o0*p0*t1*v2 - 
  2*j0*j0*s0*t1*v2 + 2*i0*n0*s0*t1*v2 + 4*k0*k1*n0*v0*v2 + 
  2*k0*k0*n1*v0*v2 - 4*j1*k0*o0*v0*v2 - 4*j0*k1*o0*v0*v2 + 
  2*i1*o0*o0*v0*v2 - 4*j0*k0*o1*v0*v2 + 4*i0*o0*o1*v0*v2 + 
  4*j0*j1*r0*v0*v2 - 2*i1*n0*r0*v0*v2 - 2*i0*n1*r0*v0*v2 + 
  2*j0*j0*r1*v0*v2 - 2*i0*n0*r1*v0*v2 + 2*k0*k0*n0*v1*v2 - 
  4*j0*k0*o0*v1*v2 + 2*i0*o0*o0*v1*v2 + 2*j0*j0*r0*v1*v2 - 
  2*i0*n0*r0*v1*v2 + 2*l1*l2*o0*o0*w0 + 4*l0*l2*o0*o1*w0 + 
  4*l0*l1*o0*o2*w0 + 2*l0*l0*o1*o2*w0 - 2*k2*l1*o0*p0*w0 - 
  2*k1*l2*o0*p0*w0 - 2*k2*l0*o1*p0*w0 - 2*k0*l2*o1*p0*w0 - 
  2*k1*l0*o2*p0*w0 - 2*k0*l1*o2*p0*w0 + 2*k1*k2*p0*p0*w0 - 
  2*k2*l0*o0*p1*w0 - 2*k0*l2*o0*p1*w0 - 2*k0*l0*o2*p1*w0 + 
  4*k0*k2*p0*p1*w0 - 2*k1*l0*o0*p2*w0 - 2*k0*l1*o0*p2*w0 - 
  2*k0*l0*o1*p2*w0 + 4*k0*k1*p0*p2*w0 + 2*k0*k0*p1*p2*w0 - 
  2*l1*l2*n0*r0*w0 - 2*l0*l2*n1*r0*w0 - 2*l0*l1*n2*r0*w0 + 
  2*j2*l1*p0*r0*w0 + 2*j1*l2*p0*r0*w0 + 2*j2*l0*p1*r0*w0 + 
  2*j0*l2*p1*r0*w0 - 2*i2*p0*p1*r0*w0 + 2*j1*l0*p2*r0*w0 + 
  2*j0*l1*p2*r0*w0 - 2*i1*p0*p2*r0*w0 - 2*i0*p1*p2*r0*w0 - 
  2*l0*l2*n0*r1*w0 - l0*l0*n2*r1*w0 + 2*j2*l0*p0*r1*w0 + 
  2*j0*l2*p0*r1*w0 - i2*p0*p0*r1*w0 + 2*j0*l0*p2*r1*w0 - 
  2*i0*p0*p2*r1*w0 - 2*l0*l1*n0*r2*w0 - l0*l0*n1*r2*w0 + 
  2*j1*l0*p0*r2*w0 + 2*j0*l1*p0*r2*w0 - i1*p0*p0*r2*w0 + 
  2*j0*l0*p1*r2*w0 - 2*i0*p0*p1*r2*w0 + 2*k2*l1*n0*s0*w0 + 
  2*k1*l2*n0*s0*w0 + 2*k2*l0*n1*s0*w0 + 2*k0*l2*n1*s0*w0 + 
  2*k1*l0*n2*s0*w0 + 2*k0*l1*n2*s0*w0 - 2*j2*l1*o0*s0*w0 - 
  2*j1*l2*o0*s0*w0 - 2*j2*l0*o1*s0*w0 - 2*j0*l2*o1*s0*w0 - 
  2*j1*l0*o2*s0*w0 - 2*j0*l1*o2*s0*w0 - 2*j2*k1*p0*s0*w0 - 
  2*j1*k2*p0*s0*w0 + 2*i2*o1*p0*s0*w0 + 2*i1*o2*p0*s0*w0 - 
  2*j2*k0*p1*s0*w0 - 2*j0*k2*p1*s0*w0 + 2*i2*o0*p1*s0*w0 + 
  2*i0*o2*p1*s0*w0 - 2*j1*k0*p2*s0*w0 - 2*j0*k1*p2*s0*w0 + 
  2*i1*o0*p2*s0*w0 + 2*i0*o1*p2*s0*w0 + 2*j1*j2*s0*s0*w0 - 
  i2*n1*s0*s0*w0 - i1*n2*s0*s0*w0 + 2*k2*l0*n0*s1*w0 + 
  2*k0*l2*n0*s1*w0 + 2*k0*l0*n2*s1*w0 - 2*j2*l0*o0*s1*w0 - 
  2*j0*l2*o0*s1*w0 - 2*j0*l0*o2*s1*w0 - 2*j2*k0*p0*s1*w0 - 
  2*j0*k2*p0*s1*w0 + 2*i2*o0*p0*s1*w0 + 2*i0*o2*p0*s1*w0 - 
  2*j0*k0*p2*s1*w0 + 2*i0*o0*p2*s1*w0 + 4*j0*j2*s0*s1*w0 - 
  2*i2*n0*s0*s1*w0 - 2*i0*n2*s0*s1*w0 + 2*k1*l0*n0*s2*w0 + 
  2*k0*l1*n0*s2*w0 + 2*k0*l0*n1*s2*w0 - 2*j1*l0*o0*s2*w0 - 
  2*j0*l1*o0*s2*w0 - 2*j0*l0*o1*s2*w0 - 2*j1*k0*p0*s2*w0 - 
  2*j0*k1*p0*s2*w0 + 2*i1*o0*p0*s2*w0 + 2*i0*o1*p0*s2*w0 - 
  2*j0*k0*p1*s2*w0 + 2*i0*o0*p1*s2*w0 + 4*j0*j1*s0*s2*w0 - 
  2*i1*n0*s0*s2*w0 - 2*i0*n1*s0*s2*w0 + 2*j0*j0*s1*s2*w0 - 
  2*i0*n0*s1*s2*w0 - 2*k1*k2*n0*u0*w0 - 2*k0*k2*n1*u0*w0 - 
  2*k0*k1*n2*u0*w0 + 2*j2*k1*o0*u0*w0 + 2*j1*k2*o0*u0*w0 + 
  2*j2*k0*o1*u0*w0 + 2*j0*k2*o1*u0*w0 - 2*i2*o0*o1*u0*w0 + 
  2*j1*k0*o2*u0*w0 + 2*j0*k1*o2*u0*w0 - 2*i1*o0*o2*u0*w0 - 
  2*i0*o1*o2*u0*w0 - 2*j1*j2*r0*u0*w0 + i2*n1*r0*u0*w0 + 
  i1*n2*r0*u0*w0 - 2*j0*j2*r1*u0*w0 + i2*n0*r1*u0*w0 + 
  i0*n2*r1*u0*w0 - 2*j0*j1*r2*u0*w0 + i1*n0*r2*u0*w0 + 
  i0*n1*r2*u0*w0 - 2*k0*k2*n0*u1*w0 - k0*k0*n2*u1*w0 + 
  2*j2*k0*o0*u1*w0 + 2*j0*k2*o0*u1*w0 - i2*o0*o0*u1*w0 + 
  2*j0*k0*o2*u1*w0 - 2*i0*o0*o2*u1*w0 - 2*j0*j2*r0*u1*w0 + 
  i2*n0*r0*u1*w0 + i0*n2*r0*u1*w0 - j0*j0*r2*u1*w0 + 
  i0*n0*r2*u1*w0 - 2*k0*k1*n0*u2*w0 - k0*k0*n1*u2*w0 + 
  2*j1*k0*o0*u2*w0 + 2*j0*k1*o0*u2*w0 - i1*o0*o0*u2*w0 + 
  2*j0*k0*o1*u2*w0 - 2*i0*o0*o1*u2*w0 - 2*j0*j1*r0*u2*w0 + 
  i1*n0*r0*u2*w0 + i0*n1*r0*u2*w0 - j0*j0*r1*u2*w0 + 
  i0*n0*r1*u2*w0 + 2*l0*l2*o0*o0*w1 + 2*l0*l0*o0*o2*w1 - 
  2*k2*l0*o0*p0*w1 - 2*k0*l2*o0*p0*w1 - 2*k0*l0*o2*p0*w1 + 
  2*k0*k2*p0*p0*w1 - 2*k0*l0*o0*p2*w1 + 2*k0*k0*p0*p2*w1 - 
  2*l0*l2*n0*r0*w1 - l0*l0*n2*r0*w1 + 2*j2*l0*p0*r0*w1 + 
  2*j0*l2*p0*r0*w1 - i2*p0*p0*r0*w1 + 2*j0*l0*p2*r0*w1 - 
  2*i0*p0*p2*r0*w1 - l0*l0*n0*r2*w1 + 2*j0*l0*p0*r2*w1 - 
  i0*p0*p0*r2*w1 + 2*k2*l0*n0*s0*w1 + 2*k0*l2*n0*s0*w1 + 
  2*k0*l0*n2*s0*w1 - 2*j2*l0*o0*s0*w1 - 2*j0*l2*o0*s0*w1 - 
  2*j0*l0*o2*s0*w1 - 2*j2*k0*p0*s0*w1 - 2*j0*k2*p0*s0*w1 + 
  2*i2*o0*p0*s0*w1 + 2*i0*o2*p0*s0*w1 - 2*j0*k0*p2*s0*w1 + 
  2*i0*o0*p2*s0*w1 + 2*j0*j2*s0*s0*w1 - i2*n0*s0*s0*w1 - 
  i0*n2*s0*s0*w1 + 2*k0*l0*n0*s2*w1 - 2*j0*l0*o0*s2*w1 - 
  2*j0*k0*p0*s2*w1 + 2*i0*o0*p0*s2*w1 + 2*j0*j0*s0*s2*w1 - 
  2*i0*n0*s0*s2*w1 - 2*k0*k2*n0*u0*w1 - k0*k0*n2*u0*w1 + 
  2*j2*k0*o0*u0*w1 + 2*j0*k2*o0*u0*w1 - i2*o0*o0*u0*w1 + 
  2*j0*k0*o2*u0*w1 - 2*i0*o0*o2*u0*w1 - 2*j0*j2*r0*u0*w1 + 
  i2*n0*r0*u0*w1 + i0*n2*r0*u0*w1 - j0*j0*r2*u0*w1 + 
  i0*n0*r2*u0*w1 - k0*k0*n0*u2*w1 + 2*j0*k0*o0*u2*w1 - 
  i0*o0*o0*u2*w1 - j0*j0*r0*u2*w1 + i0*n0*r0*u2*w1 + 
  2*l0*l1*o0*o0*w2 + 2*l0*l0*o0*o1*w2 - 2*k1*l0*o0*p0*w2 - 
  2*k0*l1*o0*p0*w2 - 2*k0*l0*o1*p0*w2 + 2*k0*k1*p0*p0*w2 - 
  2*k0*l0*o0*p1*w2 + 2*k0*k0*p0*p1*w2 - 2*l0*l1*n0*r0*w2 - 
  l0*l0*n1*r0*w2 + 2*j1*l0*p0*r0*w2 + 2*j0*l1*p0*r0*w2 - 
  i1*p0*p0*r0*w2 + 2*j0*l0*p1*r0*w2 - 2*i0*p0*p1*r0*w2 - 
  l0*l0*n0*r1*w2 + 2*j0*l0*p0*r1*w2 - i0*p0*p0*r1*w2 + 
  2*k1*l0*n0*s0*w2 + 2*k0*l1*n0*s0*w2 + 2*k0*l0*n1*s0*w2 - 
  2*j1*l0*o0*s0*w2 - 2*j0*l1*o0*s0*w2 - 2*j0*l0*o1*s0*w2 - 
  2*j1*k0*p0*s0*w2 - 2*j0*k1*p0*s0*w2 + 2*i1*o0*p0*s0*w2 + 
  2*i0*o1*p0*s0*w2 - 2*j0*k0*p1*s0*w2 + 2*i0*o0*p1*s0*w2 + 
  2*j0*j1*s0*s0*w2 - i1*n0*s0*s0*w2 - i0*n1*s0*s0*w2 + 
  2*k0*l0*n0*s1*w2 - 2*j0*l0*o0*s1*w2 - 2*j0*k0*p0*s1*w2 + 
  2*i0*o0*p0*s1*w2 + 2*j0*j0*s0*s1*w2 - 2*i0*n0*s0*s1*w2 - 
  2*k0*k1*n0*u0*w2 - k0*k0*n1*u0*w2 + 2*j1*k0*o0*u0*w2 + 
  2*j0*k1*o0*u0*w2 - i1*o0*o0*u0*w2 + 2*j0*k0*o1*u0*w2 - 
  2*i0*o0*o1*u0*w2 - 2*j0*j1*r0*u0*w2 + i1*n0*r0*u0*w2 + 
  i0*n1*r0*u0*w2 - j0*j0*r1*u0*w2 + i0*n0*r1*u0*w2 - 
  k0*k0*n0*u1*w2 + 2*j0*k0*o0*u1*w2 - i0*o0*o0*u1*w2 - 
  j0*j0*r0*u1*w2 + i0*n0*r0*u1*w2;
    
    k[12] = 4*m1*m2*p0*p1*r0 + 2*m0*m2*p1*p1*r0 + 
  2*m1*m1*p0*p2*r0 + 4*m0*m1*p1*p2*r0 - 
  2*l2*m1*p1*q0*r0 - 2*l1*m2*p1*q0*r0 - 
  2*l1*m1*p2*q0*r0 - 2*l2*m1*p0*q1*r0 - 
  2*l1*m2*p0*q1*r0 - 2*l2*m0*p1*q1*r0 - 
  2*l0*m2*p1*q1*r0 - 2*l1*m0*p2*q1*r0 - 
  2*l0*m1*p2*q1*r0 + 4*l1*l2*q0*q1*r0 + 
  2*l0*l2*q1*q1*r0 - 2*l1*m1*p0*q2*r0 - 
  2*l1*m0*p1*q2*r0 - 2*l0*m1*p1*q2*r0 + 
  2*l1*l1*q0*q2*r0 + 4*l0*l1*q1*q2*r0 + 
  2*m1*m2*p0*p0*r1 + 4*m0*m2*p0*p1*r1 + 
  4*m0*m1*p0*p2*r1 + 2*m0*m0*p1*p2*r1 - 
  2*l2*m1*p0*q0*r1 - 2*l1*m2*p0*q0*r1 - 
  2*l2*m0*p1*q0*r1 - 2*l0*m2*p1*q0*r1 - 
  2*l1*m0*p2*q0*r1 - 2*l0*m1*p2*q0*r1 + 
  2*l1*l2*q0*q0*r1 - 2*l2*m0*p0*q1*r1 - 
  2*l0*m2*p0*q1*r1 - 2*l0*m0*p2*q1*r1 + 
  4*l0*l2*q0*q1*r1 - 2*l1*m0*p0*q2*r1 - 
  2*l0*m1*p0*q2*r1 - 2*l0*m0*p1*q2*r1 + 
  4*l0*l1*q0*q2*r1 + 2*l0*l0*q1*q2*r1 + m1*m1*p0*p0*r2 + 
  4*m0*m1*p0*p1*r2 + m0*m0*p1*p1*r2 - 
  2*l1*m1*p0*q0*r2 - 2*l1*m0*p1*q0*r2 - 
  2*l0*m1*p1*q0*r2 + l1*l1*q0*q0*r2 - 
  2*l1*m0*p0*q1*r2 - 2*l0*m1*p0*q1*r2 - 
  2*l0*m0*p1*q1*r2 + 4*l0*l1*q0*q1*r2 + 
  l0*l0*q1*q1*r2 - 4*m1*m2*o1*p0*s0 - 2*m1*m1*o2*p0*s0 - 
  4*m1*m2*o0*p1*s0 - 4*m0*m2*o1*p1*s0 - 
  4*m0*m1*o2*p1*s0 - 2*m1*m1*o0*p2*s0 - 
  4*m0*m1*o1*p2*s0 + 2*l2*m1*o1*q0*s0 + 
  2*l1*m2*o1*q0*s0 + 2*l1*m1*o2*q0*s0 + 
  2*k2*m1*p1*q0*s0 + 2*k1*m2*p1*q0*s0 + 
  2*k1*m1*p2*q0*s0 + 2*l2*m1*o0*q1*s0 + 
  2*l1*m2*o0*q1*s0 + 2*l2*m0*o1*q1*s0 + 
  2*l0*m2*o1*q1*s0 + 2*l1*m0*o2*q1*s0 + 
  2*l0*m1*o2*q1*s0 + 2*k2*m1*p0*q1*s0 + 
  2*k1*m2*p0*q1*s0 + 2*k2*m0*p1*q1*s0 + 
  2*k0*m2*p1*q1*s0 + 2*k1*m0*p2*q1*s0 + 
  2*k0*m1*p2*q1*s0 - 4*k2*l1*q0*q1*s0 - 
  4*k1*l2*q0*q1*s0 - 2*k2*l0*q1*q1*s0 - 
  2*k0*l2*q1*q1*s0 + 2*l1*m1*o0*q2*s0 + 
  2*l1*m0*o1*q2*s0 + 2*l0*m1*o1*q2*s0 + 
  2*k1*m1*p0*q2*s0 + 2*k1*m0*p1*q2*s0 + 
  2*k0*m1*p1*q2*s0 - 4*k1*l1*q0*q2*s0 - 
  4*k1*l0*q1*q2*s0 - 4*k0*l1*q1*q2*s0 + 
  2*m1*m2*n1*s0*s0 + m1*m1*n2*s0*s0 - 2*j2*m1*q1*s0*s0 - 
  2*j1*m2*q1*s0*s0 + i2*q1*q1*s0*s0 - 2*j1*m1*q2*s0*s0 + 
  2*i1*q1*q2*s0*s0 - 4*m1*m2*o0*p0*s1 - 
  4*m0*m2*o1*p0*s1 - 4*m0*m1*o2*p0*s1 - 
  4*m0*m2*o0*p1*s1 - 2*m0*m0*o2*p1*s1 - 
  4*m0*m1*o0*p2*s1 - 2*m0*m0*o1*p2*s1 + 
  2*l2*m1*o0*q0*s1 + 2*l1*m2*o0*q0*s1 + 
  2*l2*m0*o1*q0*s1 + 2*l0*m2*o1*q0*s1 + 
  2*l1*m0*o2*q0*s1 + 2*l0*m1*o2*q0*s1 + 
  2*k2*m1*p0*q0*s1 + 2*k1*m2*p0*q0*s1 + 
  2*k2*m0*p1*q0*s1 + 2*k0*m2*p1*q0*s1 + 
  2*k1*m0*p2*q0*s1 + 2*k0*m1*p2*q0*s1 - 
  2*k2*l1*q0*q0*s1 - 2*k1*l2*q0*q0*s1 + 
  2*l2*m0*o0*q1*s1 + 2*l0*m2*o0*q1*s1 + 
  2*l0*m0*o2*q1*s1 + 2*k2*m0*p0*q1*s1 + 
  2*k0*m2*p0*q1*s1 + 2*k0*m0*p2*q1*s1 - 
  4*k2*l0*q0*q1*s1 - 4*k0*l2*q0*q1*s1 + 
  2*l1*m0*o0*q2*s1 + 2*l0*m1*o0*q2*s1 + 
  2*l0*m0*o1*q2*s1 + 2*k1*m0*p0*q2*s1 + 
  2*k0*m1*p0*q2*s1 + 2*k0*m0*p1*q2*s1 - 
  4*k1*l0*q0*q2*s1 - 4*k0*l1*q0*q2*s1 - 
  4*k0*l0*q1*q2*s1 + 4*m1*m2*n0*s0*s1 + 
  4*m0*m2*n1*s0*s1 + 4*m0*m1*n2*s0*s1 - 
  4*j2*m1*q0*s0*s1 - 4*j1*m2*q0*s0*s1 - 
  4*j2*m0*q1*s0*s1 - 4*j0*m2*q1*s0*s1 + 
  4*i2*q0*q1*s0*s1 - 4*j1*m0*q2*s0*s1 - 
  4*j0*m1*q2*s0*s1 + 4*i1*q0*q2*s0*s1 + 
  4*i0*q1*q2*s0*s1 + 2*m0*m2*n0*s1*s1 + m0*m0*n2*s1*s1 - 
  2*j2*m0*q0*s1*s1 - 2*j0*m2*q0*s1*s1 + i2*q0*q0*s1*s1 - 
  2*j0*m0*q2*s1*s1 + 2*i0*q0*q2*s1*s1 - 
  2*m1*m1*o0*p0*s2 - 4*m0*m1*o1*p0*s2 - 
  4*m0*m1*o0*p1*s2 - 2*m0*m0*o1*p1*s2 + 
  2*l1*m1*o0*q0*s2 + 2*l1*m0*o1*q0*s2 + 
  2*l0*m1*o1*q0*s2 + 2*k1*m1*p0*q0*s2 + 
  2*k1*m0*p1*q0*s2 + 2*k0*m1*p1*q0*s2 - 
  2*k1*l1*q0*q0*s2 + 2*l1*m0*o0*q1*s2 + 
  2*l0*m1*o0*q1*s2 + 2*l0*m0*o1*q1*s2 + 
  2*k1*m0*p0*q1*s2 + 2*k0*m1*p0*q1*s2 + 
  2*k0*m0*p1*q1*s2 - 4*k1*l0*q0*q1*s2 - 
  4*k0*l1*q0*q1*s2 - 2*k0*l0*q1*q1*s2 + 
  2*m1*m1*n0*s0*s2 + 4*m0*m1*n1*s0*s2 - 
  4*j1*m1*q0*s0*s2 - 4*j1*m0*q1*s0*s2 - 
  4*j0*m1*q1*s0*s2 + 4*i1*q0*q1*s0*s2 + 
  2*i0*q1*q1*s0*s2 + 4*m0*m1*n0*s1*s2 + 
  2*m0*m0*n1*s1*s2 - 4*j1*m0*q0*s1*s2 - 
  4*j0*m1*q0*s1*s2 + 2*i1*q0*q0*s1*s2 - 
  4*j0*m0*q1*s1*s2 + 4*i0*q0*q1*s1*s2 + 
  2*l2*m1*o1*p0*t0 + 2*l1*m2*o1*p0*t0 + 
  2*l1*m1*o2*p0*t0 + 2*l2*m1*o0*p1*t0 + 
  2*l1*m2*o0*p1*t0 + 2*l2*m0*o1*p1*t0 + 
  2*l0*m2*o1*p1*t0 + 2*l1*m0*o2*p1*t0 + 
  2*l0*m1*o2*p1*t0 - 4*k2*m1*p0*p1*t0 - 
  4*k1*m2*p0*p1*t0 - 2*k2*m0*p1*p1*t0 - 
  2*k0*m2*p1*p1*t0 + 2*l1*m1*o0*p2*t0 + 
  2*l1*m0*o1*p2*t0 + 2*l0*m1*o1*p2*t0 - 
  4*k1*m1*p0*p2*t0 - 4*k1*m0*p1*p2*t0 - 
  4*k0*m1*p1*p2*t0 - 4*l1*l2*o1*q0*t0 - 
  2*l1*l1*o2*q0*t0 + 2*k2*l1*p1*q0*t0 + 
  2*k1*l2*p1*q0*t0 + 2*k1*l1*p2*q0*t0 - 
  4*l1*l2*o0*q1*t0 - 4*l0*l2*o1*q1*t0 - 
  4*l0*l1*o2*q1*t0 + 2*k2*l1*p0*q1*t0 + 
  2*k1*l2*p0*q1*t0 + 2*k2*l0*p1*q1*t0 + 
  2*k0*l2*p1*q1*t0 + 2*k1*l0*p2*q1*t0 + 
  2*k0*l1*p2*q1*t0 - 2*l1*l1*o0*q2*t0 - 
  4*l0*l1*o1*q2*t0 + 2*k1*l1*p0*q2*t0 + 
  2*k1*l0*p1*q2*t0 + 2*k0*l1*p1*q2*t0 - 
  2*l2*m1*n1*s0*t0 - 2*l1*m2*n1*s0*t0 - 
  2*l1*m1*n2*s0*t0 + 2*j2*m1*p1*s0*t0 + 
  2*j1*m2*p1*s0*t0 + 2*j1*m1*p2*s0*t0 + 
  2*j2*l1*q1*s0*t0 + 2*j1*l2*q1*s0*t0 - 
  2*i2*p1*q1*s0*t0 - 2*i1*p2*q1*s0*t0 + 
  2*j1*l1*q2*s0*t0 - 2*i1*p1*q2*s0*t0 - 
  2*l2*m1*n0*s1*t0 - 2*l1*m2*n0*s1*t0 - 
  2*l2*m0*n1*s1*t0 - 2*l0*m2*n1*s1*t0 - 
  2*l1*m0*n2*s1*t0 - 2*l0*m1*n2*s1*t0 + 
  2*j2*m1*p0*s1*t0 + 2*j1*m2*p0*s1*t0 + 
  2*j2*m0*p1*s1*t0 + 2*j0*m2*p1*s1*t0 + 
  2*j1*m0*p2*s1*t0 + 2*j0*m1*p2*s1*t0 + 
  2*j2*l1*q0*s1*t0 + 2*j1*l2*q0*s1*t0 - 
  2*i2*p1*q0*s1*t0 - 2*i1*p2*q0*s1*t0 + 
  2*j2*l0*q1*s1*t0 + 2*j0*l2*q1*s1*t0 - 
  2*i2*p0*q1*s1*t0 - 2*i0*p2*q1*s1*t0 + 
  2*j1*l0*q2*s1*t0 + 2*j0*l1*q2*s1*t0 - 
  2*i1*p0*q2*s1*t0 - 2*i0*p1*q2*s1*t0 - 
  2*l1*m1*n0*s2*t0 - 2*l1*m0*n1*s2*t0 - 
  2*l0*m1*n1*s2*t0 + 2*j1*m1*p0*s2*t0 + 
  2*j1*m0*p1*s2*t0 + 2*j0*m1*p1*s2*t0 + 
  2*j1*l1*q0*s2*t0 - 2*i1*p1*q0*s2*t0 + 
  2*j1*l0*q1*s2*t0 + 2*j0*l1*q1*s2*t0 - 
  2*i1*p0*q1*s2*t0 - 2*i0*p1*q1*s2*t0 + 
  2*l1*l2*n1*t0*t0 + l1*l1*n2*t0*t0 - 2*j2*l1*p1*t0*t0 - 
  2*j1*l2*p1*t0*t0 + i2*p1*p1*t0*t0 - 2*j1*l1*p2*t0*t0 + 
  2*i1*p1*p2*t0*t0 + 2*l2*m1*o0*p0*t1 + 
  2*l1*m2*o0*p0*t1 + 2*l2*m0*o1*p0*t1 + 
  2*l0*m2*o1*p0*t1 + 2*l1*m0*o2*p0*t1 + 
  2*l0*m1*o2*p0*t1 - 2*k2*m1*p0*p0*t1 - 
  2*k1*m2*p0*p0*t1 + 2*l2*m0*o0*p1*t1 + 
  2*l0*m2*o0*p1*t1 + 2*l0*m0*o2*p1*t1 - 
  4*k2*m0*p0*p1*t1 - 4*k0*m2*p0*p1*t1 + 
  2*l1*m0*o0*p2*t1 + 2*l0*m1*o0*p2*t1 + 
  2*l0*m0*o1*p2*t1 - 4*k1*m0*p0*p2*t1 - 
  4*k0*m1*p0*p2*t1 - 4*k0*m0*p1*p2*t1 - 
  4*l1*l2*o0*q0*t1 - 4*l0*l2*o1*q0*t1 - 
  4*l0*l1*o2*q0*t1 + 2*k2*l1*p0*q0*t1 + 
  2*k1*l2*p0*q0*t1 + 2*k2*l0*p1*q0*t1 + 
  2*k0*l2*p1*q0*t1 + 2*k1*l0*p2*q0*t1 + 
  2*k0*l1*p2*q0*t1 - 4*l0*l2*o0*q1*t1 - 
  2*l0*l0*o2*q1*t1 + 2*k2*l0*p0*q1*t1 + 
  2*k0*l2*p0*q1*t1 + 2*k0*l0*p2*q1*t1 - 
  4*l0*l1*o0*q2*t1 - 2*l0*l0*o1*q2*t1 + 
  2*k1*l0*p0*q2*t1 + 2*k0*l1*p0*q2*t1 + 
  2*k0*l0*p1*q2*t1 - 2*l2*m1*n0*s0*t1 - 
  2*l1*m2*n0*s0*t1 - 2*l2*m0*n1*s0*t1 - 
  2*l0*m2*n1*s0*t1 - 2*l1*m0*n2*s0*t1 - 
  2*l0*m1*n2*s0*t1 + 2*j2*m1*p0*s0*t1 + 
  2*j1*m2*p0*s0*t1 + 2*j2*m0*p1*s0*t1 + 
  2*j0*m2*p1*s0*t1 + 2*j1*m0*p2*s0*t1 + 
  2*j0*m1*p2*s0*t1 + 2*j2*l1*q0*s0*t1 + 
  2*j1*l2*q0*s0*t1 - 2*i2*p1*q0*s0*t1 - 
  2*i1*p2*q0*s0*t1 + 2*j2*l0*q1*s0*t1 + 
  2*j0*l2*q1*s0*t1 - 2*i2*p0*q1*s0*t1 - 
  2*i0*p2*q1*s0*t1 + 2*j1*l0*q2*s0*t1 + 
  2*j0*l1*q2*s0*t1 - 2*i1*p0*q2*s0*t1 - 
  2*i0*p1*q2*s0*t1 - 2*l2*m0*n0*s1*t1 - 
  2*l0*m2*n0*s1*t1 - 2*l0*m0*n2*s1*t1 + 
  2*j2*m0*p0*s1*t1 + 2*j0*m2*p0*s1*t1 + 
  2*j0*m0*p2*s1*t1 + 2*j2*l0*q0*s1*t1 + 
  2*j0*l2*q0*s1*t1 - 2*i2*p0*q0*s1*t1 - 
  2*i0*p2*q0*s1*t1 + 2*j0*l0*q2*s1*t1 - 
  2*i0*p0*q2*s1*t1 - 2*l1*m0*n0*s2*t1 - 
  2*l0*m1*n0*s2*t1 - 2*l0*m0*n1*s2*t1 + 
  2*j1*m0*p0*s2*t1 + 2*j0*m1*p0*s2*t1 + 
  2*j0*m0*p1*s2*t1 + 2*j1*l0*q0*s2*t1 + 
  2*j0*l1*q0*s2*t1 - 2*i1*p0*q0*s2*t1 - 
  2*i0*p1*q0*s2*t1 + 2*j0*l0*q1*s2*t1 - 
  2*i0*p0*q1*s2*t1 + 4*l1*l2*n0*t0*t1 + 
  4*l0*l2*n1*t0*t1 + 4*l0*l1*n2*t0*t1 - 
  4*j2*l1*p0*t0*t1 - 4*j1*l2*p0*t0*t1 - 
  4*j2*l0*p1*t0*t1 - 4*j0*l2*p1*t0*t1 + 
  4*i2*p0*p1*t0*t1 - 4*j1*l0*p2*t0*t1 - 
  4*j0*l1*p2*t0*t1 + 4*i1*p0*p2*t0*t1 + 
  4*i0*p1*p2*t0*t1 + 2*l0*l2*n0*t1*t1 + l0*l0*n2*t1*t1 - 
  2*j2*l0*p0*t1*t1 - 2*j0*l2*p0*t1*t1 + i2*p0*p0*t1*t1 - 
  2*j0*l0*p2*t1*t1 + 2*i0*p0*p2*t1*t1 + 
  2*l1*m1*o0*p0*t2 + 2*l1*m0*o1*p0*t2 + 
  2*l0*m1*o1*p0*t2 - 2*k1*m1*p0*p0*t2 + 
  2*l1*m0*o0*p1*t2 + 2*l0*m1*o0*p1*t2 + 
  2*l0*m0*o1*p1*t2 - 4*k1*m0*p0*p1*t2 - 
  4*k0*m1*p0*p1*t2 - 2*k0*m0*p1*p1*t2 - 
  2*l1*l1*o0*q0*t2 - 4*l0*l1*o1*q0*t2 + 
  2*k1*l1*p0*q0*t2 + 2*k1*l0*p1*q0*t2 + 
  2*k0*l1*p1*q0*t2 - 4*l0*l1*o0*q1*t2 - 
  2*l0*l0*o1*q1*t2 + 2*k1*l0*p0*q1*t2 + 
  2*k0*l1*p0*q1*t2 + 2*k0*l0*p1*q1*t2 - 
  2*l1*m1*n0*s0*t2 - 2*l1*m0*n1*s0*t2 - 
  2*l0*m1*n1*s0*t2 + 2*j1*m1*p0*s0*t2 + 
  2*j1*m0*p1*s0*t2 + 2*j0*m1*p1*s0*t2 + 
  2*j1*l1*q0*s0*t2 - 2*i1*p1*q0*s0*t2 + 
  2*j1*l0*q1*s0*t2 + 2*j0*l1*q1*s0*t2 - 
  2*i1*p0*q1*s0*t2 - 2*i0*p1*q1*s0*t2 - 
  2*l1*m0*n0*s1*t2 - 2*l0*m1*n0*s1*t2 - 
  2*l0*m0*n1*s1*t2 + 2*j1*m0*p0*s1*t2 + 
  2*j0*m1*p0*s1*t2 + 2*j0*m0*p1*s1*t2 + 
  2*j1*l0*q0*s1*t2 + 2*j0*l1*q0*s1*t2 - 
  2*i1*p0*q0*s1*t2 - 2*i0*p1*q0*s1*t2 + 
  2*j0*l0*q1*s1*t2 - 2*i0*p0*q1*s1*t2 + 
  2*l1*l1*n0*t0*t2 + 4*l0*l1*n1*t0*t2 - 
  4*j1*l1*p0*t0*t2 - 4*j1*l0*p1*t0*t2 - 
  4*j0*l1*p1*t0*t2 + 4*i1*p0*p1*t0*t2 + 
  2*i0*p1*p1*t0*t2 + 4*l0*l1*n0*t1*t2 + 
  2*l0*l0*n1*t1*t2 - 4*j1*l0*p0*t1*t2 - 
  4*j0*l1*p0*t1*t2 + 2*i1*p0*p0*t1*t2 - 
  4*j0*l0*p1*t1*t2 + 4*i0*p0*p1*t1*t2 + 
  4*m1*m2*o0*o1*u0 + 2*m0*m2*o1*o1*u0 + 
  2*m1*m1*o0*o2*u0 + 4*m0*m1*o1*o2*u0 - 
  2*k2*m1*o1*q0*u0 - 2*k1*m2*o1*q0*u0 - 
  2*k1*m1*o2*q0*u0 - 2*k2*m1*o0*q1*u0 - 
  2*k1*m2*o0*q1*u0 - 2*k2*m0*o1*q1*u0 - 
  2*k0*m2*o1*q1*u0 - 2*k1*m0*o2*q1*u0 - 
  2*k0*m1*o2*q1*u0 + 4*k1*k2*q0*q1*u0 + 
  2*k0*k2*q1*q1*u0 - 2*k1*m1*o0*q2*u0 - 
  2*k1*m0*o1*q2*u0 - 2*k0*m1*o1*q2*u0 + 
  2*k1*k1*q0*q2*u0 + 4*k0*k1*q1*q2*u0 - 
  2*m1*m2*n1*r0*u0 - m1*m1*n2*r0*u0 + 
  2*j2*m1*q1*r0*u0 + 2*j1*m2*q1*r0*u0 - 
  i2*q1*q1*r0*u0 + 2*j1*m1*q2*r0*u0 - 
  2*i1*q1*q2*r0*u0 - 2*m1*m2*n0*r1*u0 - 
  2*m0*m2*n1*r1*u0 - 2*m0*m1*n2*r1*u0 + 
  2*j2*m1*q0*r1*u0 + 2*j1*m2*q0*r1*u0 + 
  2*j2*m0*q1*r1*u0 + 2*j0*m2*q1*r1*u0 - 
  2*i2*q0*q1*r1*u0 + 2*j1*m0*q2*r1*u0 + 
  2*j0*m1*q2*r1*u0 - 2*i1*q0*q2*r1*u0 - 
  2*i0*q1*q2*r1*u0 - m1*m1*n0*r2*u0 - 
  2*m0*m1*n1*r2*u0 + 2*j1*m1*q0*r2*u0 + 
  2*j1*m0*q1*r2*u0 + 2*j0*m1*q1*r2*u0 - 
  2*i1*q0*q1*r2*u0 - i0*q1*q1*r2*u0 + 
  2*k2*m1*n1*t0*u0 + 2*k1*m2*n1*t0*u0 + 
  2*k1*m1*n2*t0*u0 - 2*j2*m1*o1*t0*u0 - 
  2*j1*m2*o1*t0*u0 - 2*j1*m1*o2*t0*u0 - 
  2*j2*k1*q1*t0*u0 - 2*j1*k2*q1*t0*u0 + 
  2*i2*o1*q1*t0*u0 + 2*i1*o2*q1*t0*u0 - 
  2*j1*k1*q2*t0*u0 + 2*i1*o1*q2*t0*u0 + 
  2*k2*m1*n0*t1*u0 + 2*k1*m2*n0*t1*u0 + 
  2*k2*m0*n1*t1*u0 + 2*k0*m2*n1*t1*u0 + 
  2*k1*m0*n2*t1*u0 + 2*k0*m1*n2*t1*u0 - 
  2*j2*m1*o0*t1*u0 - 2*j1*m2*o0*t1*u0 - 
  2*j2*m0*o1*t1*u0 - 2*j0*m2*o1*t1*u0 - 
  2*j1*m0*o2*t1*u0 - 2*j0*m1*o2*t1*u0 - 
  2*j2*k1*q0*t1*u0 - 2*j1*k2*q0*t1*u0 + 
  2*i2*o1*q0*t1*u0 + 2*i1*o2*q0*t1*u0 - 
  2*j2*k0*q1*t1*u0 - 2*j0*k2*q1*t1*u0 + 
  2*i2*o0*q1*t1*u0 + 2*i0*o2*q1*t1*u0 - 
  2*j1*k0*q2*t1*u0 - 2*j0*k1*q2*t1*u0 + 
  2*i1*o0*q2*t1*u0 + 2*i0*o1*q2*t1*u0 + 
  4*j1*j2*t0*t1*u0 - 2*i2*n1*t0*t1*u0 - 
  2*i1*n2*t0*t1*u0 + 2*j0*j2*t1*t1*u0 - 
  i2*n0*t1*t1*u0 - i0*n2*t1*t1*u0 + 2*k1*m1*n0*t2*u0 + 
  2*k1*m0*n1*t2*u0 + 2*k0*m1*n1*t2*u0 - 
  2*j1*m1*o0*t2*u0 - 2*j1*m0*o1*t2*u0 - 
  2*j0*m1*o1*t2*u0 - 2*j1*k1*q0*t2*u0 + 
  2*i1*o1*q0*t2*u0 - 2*j1*k0*q1*t2*u0 - 
  2*j0*k1*q1*t2*u0 + 2*i1*o0*q1*t2*u0 + 
  2*i0*o1*q1*t2*u0 + 2*j1*j1*t0*t2*u0 - 
  2*i1*n1*t0*t2*u0 + 4*j0*j1*t1*t2*u0 - 
  2*i1*n0*t1*t2*u0 - 2*i0*n1*t1*t2*u0 + 
  2*m1*m2*o0*o0*u1 + 4*m0*m2*o0*o1*u1 + 
  4*m0*m1*o0*o2*u1 + 2*m0*m0*o1*o2*u1 - 
  2*k2*m1*o0*q0*u1 - 2*k1*m2*o0*q0*u1 - 
  2*k2*m0*o1*q0*u1 - 2*k0*m2*o1*q0*u1 - 
  2*k1*m0*o2*q0*u1 - 2*k0*m1*o2*q0*u1 + 
  2*k1*k2*q0*q0*u1 - 2*k2*m0*o0*q1*u1 - 
  2*k0*m2*o0*q1*u1 - 2*k0*m0*o2*q1*u1 + 
  4*k0*k2*q0*q1*u1 - 2*k1*m0*o0*q2*u1 - 
  2*k0*m1*o0*q2*u1 - 2*k0*m0*o1*q2*u1 + 
  4*k0*k1*q0*q2*u1 + 2*k0*k0*q1*q2*u1 - 
  2*m1*m2*n0*r0*u1 - 2*m0*m2*n1*r0*u1 - 
  2*m0*m1*n2*r0*u1 + 2*j2*m1*q0*r0*u1 + 
  2*j1*m2*q0*r0*u1 + 2*j2*m0*q1*r0*u1 + 
  2*j0*m2*q1*r0*u1 - 2*i2*q0*q1*r0*u1 + 
  2*j1*m0*q2*r0*u1 + 2*j0*m1*q2*r0*u1 - 
  2*i1*q0*q2*r0*u1 - 2*i0*q1*q2*r0*u1 - 
  2*m0*m2*n0*r1*u1 - m0*m0*n2*r1*u1 + 
  2*j2*m0*q0*r1*u1 + 2*j0*m2*q0*r1*u1 - 
  i2*q0*q0*r1*u1 + 2*j0*m0*q2*r1*u1 - 
  2*i0*q0*q2*r1*u1 - 2*m0*m1*n0*r2*u1 - 
  m0*m0*n1*r2*u1 + 2*j1*m0*q0*r2*u1 + 
  2*j0*m1*q0*r2*u1 - i1*q0*q0*r2*u1 + 
  2*j0*m0*q1*r2*u1 - 2*i0*q0*q1*r2*u1 + 
  2*k2*m1*n0*t0*u1 + 2*k1*m2*n0*t0*u1 + 
  2*k2*m0*n1*t0*u1 + 2*k0*m2*n1*t0*u1 + 
  2*k1*m0*n2*t0*u1 + 2*k0*m1*n2*t0*u1 - 
  2*j2*m1*o0*t0*u1 - 2*j1*m2*o0*t0*u1 - 
  2*j2*m0*o1*t0*u1 - 2*j0*m2*o1*t0*u1 - 
  2*j1*m0*o2*t0*u1 - 2*j0*m1*o2*t0*u1 - 
  2*j2*k1*q0*t0*u1 - 2*j1*k2*q0*t0*u1 + 
  2*i2*o1*q0*t0*u1 + 2*i1*o2*q0*t0*u1 - 
  2*j2*k0*q1*t0*u1 - 2*j0*k2*q1*t0*u1 + 
  2*i2*o0*q1*t0*u1 + 2*i0*o2*q1*t0*u1 - 
  2*j1*k0*q2*t0*u1 - 2*j0*k1*q2*t0*u1 + 
  2*i1*o0*q2*t0*u1 + 2*i0*o1*q2*t0*u1 + 
  2*j1*j2*t0*t0*u1 - i2*n1*t0*t0*u1 - i1*n2*t0*t0*u1 + 
  2*k2*m0*n0*t1*u1 + 2*k0*m2*n0*t1*u1 + 
  2*k0*m0*n2*t1*u1 - 2*j2*m0*o0*t1*u1 - 
  2*j0*m2*o0*t1*u1 - 2*j0*m0*o2*t1*u1 - 
  2*j2*k0*q0*t1*u1 - 2*j0*k2*q0*t1*u1 + 
  2*i2*o0*q0*t1*u1 + 2*i0*o2*q0*t1*u1 - 
  2*j0*k0*q2*t1*u1 + 2*i0*o0*q2*t1*u1 + 
  4*j0*j2*t0*t1*u1 - 2*i2*n0*t0*t1*u1 - 
  2*i0*n2*t0*t1*u1 + 2*k1*m0*n0*t2*u1 + 
  2*k0*m1*n0*t2*u1 + 2*k0*m0*n1*t2*u1 - 
  2*j1*m0*o0*t2*u1 - 2*j0*m1*o0*t2*u1 - 
  2*j0*m0*o1*t2*u1 - 2*j1*k0*q0*t2*u1 - 
  2*j0*k1*q0*t2*u1 + 2*i1*o0*q0*t2*u1 + 
  2*i0*o1*q0*t2*u1 - 2*j0*k0*q1*t2*u1 + 
  2*i0*o0*q1*t2*u1 + 4*j0*j1*t0*t2*u1 - 
  2*i1*n0*t0*t2*u1 - 2*i0*n1*t0*t2*u1 + 
  2*j0*j0*t1*t2*u1 - 2*i0*n0*t1*t2*u1 + m1*m1*o0*o0*u2 + 
  4*m0*m1*o0*o1*u2 + m0*m0*o1*o1*u2 - 
  2*k1*m1*o0*q0*u2 - 2*k1*m0*o1*q0*u2 - 
  2*k0*m1*o1*q0*u2 + k1*k1*q0*q0*u2 - 
  2*k1*m0*o0*q1*u2 - 2*k0*m1*o0*q1*u2 - 
  2*k0*m0*o1*q1*u2 + 4*k0*k1*q0*q1*u2 + 
  k0*k0*q1*q1*u2 - m1*m1*n0*r0*u2 - 2*m0*m1*n1*r0*u2 + 
  2*j1*m1*q0*r0*u2 + 2*j1*m0*q1*r0*u2 + 
  2*j0*m1*q1*r0*u2 - 2*i1*q0*q1*r0*u2 - 
  i0*q1*q1*r0*u2 - 2*m0*m1*n0*r1*u2 - m0*m0*n1*r1*u2 + 
  2*j1*m0*q0*r1*u2 + 2*j0*m1*q0*r1*u2 - 
  i1*q0*q0*r1*u2 + 2*j0*m0*q1*r1*u2 - 
  2*i0*q0*q1*r1*u2 + 2*k1*m1*n0*t0*u2 + 
  2*k1*m0*n1*t0*u2 + 2*k0*m1*n1*t0*u2 - 
  2*j1*m1*o0*t0*u2 - 2*j1*m0*o1*t0*u2 - 
  2*j0*m1*o1*t0*u2 - 2*j1*k1*q0*t0*u2 + 
  2*i1*o1*q0*t0*u2 - 2*j1*k0*q1*t0*u2 - 
  2*j0*k1*q1*t0*u2 + 2*i1*o0*q1*t0*u2 + 
  2*i0*o1*q1*t0*u2 + j1*j1*t0*t0*u2 - i1*n1*t0*t0*u2 + 
  2*k1*m0*n0*t1*u2 + 2*k0*m1*n0*t1*u2 + 
  2*k0*m0*n1*t1*u2 - 2*j1*m0*o0*t1*u2 - 
  2*j0*m1*o0*t1*u2 - 2*j0*m0*o1*t1*u2 - 
  2*j1*k0*q0*t1*u2 - 2*j0*k1*q0*t1*u2 + 
  2*i1*o0*q0*t1*u2 + 2*i0*o1*q0*t1*u2 - 
  2*j0*k0*q1*t1*u2 + 2*i0*o0*q1*t1*u2 + 
  4*j0*j1*t0*t1*u2 - 2*i1*n0*t0*t1*u2 - 
  2*i0*n1*t0*t1*u2 + j0*j0*t1*t1*u2 - i0*n0*t1*t1*u2 - 
  4*l2*m1*o0*o1*v0 - 4*l1*m2*o0*o1*v0 - 
  2*l2*m0*o1*o1*v0 - 2*l0*m2*o1*o1*v0 - 
  4*l1*m1*o0*o2*v0 - 4*l1*m0*o1*o2*v0 - 
  4*l0*m1*o1*o2*v0 + 2*k2*m1*o1*p0*v0 + 
  2*k1*m2*o1*p0*v0 + 2*k1*m1*o2*p0*v0 + 
  2*k2*m1*o0*p1*v0 + 2*k1*m2*o0*p1*v0 + 
  2*k2*m0*o1*p1*v0 + 2*k0*m2*o1*p1*v0 + 
  2*k1*m0*o2*p1*v0 + 2*k0*m1*o2*p1*v0 + 
  2*k1*m1*o0*p2*v0 + 2*k1*m0*o1*p2*v0 + 
  2*k0*m1*o1*p2*v0 + 2*k2*l1*o1*q0*v0 + 
  2*k1*l2*o1*q0*v0 + 2*k1*l1*o2*q0*v0 - 
  4*k1*k2*p1*q0*v0 - 2*k1*k1*p2*q0*v0 + 
  2*k2*l1*o0*q1*v0 + 2*k1*l2*o0*q1*v0 + 
  2*k2*l0*o1*q1*v0 + 2*k0*l2*o1*q1*v0 + 
  2*k1*l0*o2*q1*v0 + 2*k0*l1*o2*q1*v0 - 
  4*k1*k2*p0*q1*v0 - 4*k0*k2*p1*q1*v0 - 
  4*k0*k1*p2*q1*v0 + 2*k1*l1*o0*q2*v0 + 
  2*k1*l0*o1*q2*v0 + 2*k0*l1*o1*q2*v0 - 
  2*k1*k1*p0*q2*v0 - 4*k0*k1*p1*q2*v0 + 
  2*l2*m1*n1*r0*v0 + 2*l1*m2*n1*r0*v0 + 
  2*l1*m1*n2*r0*v0 - 2*j2*m1*p1*r0*v0 - 
  2*j1*m2*p1*r0*v0 - 2*j1*m1*p2*r0*v0 - 
  2*j2*l1*q1*r0*v0 - 2*j1*l2*q1*r0*v0 + 
  2*i2*p1*q1*r0*v0 + 2*i1*p2*q1*r0*v0 - 
  2*j1*l1*q2*r0*v0 + 2*i1*p1*q2*r0*v0 + 
  2*l2*m1*n0*r1*v0 + 2*l1*m2*n0*r1*v0 + 
  2*l2*m0*n1*r1*v0 + 2*l0*m2*n1*r1*v0 + 
  2*l1*m0*n2*r1*v0 + 2*l0*m1*n2*r1*v0 - 
  2*j2*m1*p0*r1*v0 - 2*j1*m2*p0*r1*v0 - 
  2*j2*m0*p1*r1*v0 - 2*j0*m2*p1*r1*v0 - 
  2*j1*m0*p2*r1*v0 - 2*j0*m1*p2*r1*v0 - 
  2*j2*l1*q0*r1*v0 - 2*j1*l2*q0*r1*v0 + 
  2*i2*p1*q0*r1*v0 + 2*i1*p2*q0*r1*v0 - 
  2*j2*l0*q1*r1*v0 - 2*j0*l2*q1*r1*v0 + 
  2*i2*p0*q1*r1*v0 + 2*i0*p2*q1*r1*v0 - 
  2*j1*l0*q2*r1*v0 - 2*j0*l1*q2*r1*v0 + 
  2*i1*p0*q2*r1*v0 + 2*i0*p1*q2*r1*v0 + 
  2*l1*m1*n0*r2*v0 + 2*l1*m0*n1*r2*v0 + 
  2*l0*m1*n1*r2*v0 - 2*j1*m1*p0*r2*v0 - 
  2*j1*m0*p1*r2*v0 - 2*j0*m1*p1*r2*v0 - 
  2*j1*l1*q0*r2*v0 + 2*i1*p1*q0*r2*v0 - 
  2*j1*l0*q1*r2*v0 - 2*j0*l1*q1*r2*v0 + 
  2*i1*p0*q1*r2*v0 + 2*i0*p1*q1*r2*v0 - 
  2*k2*m1*n1*s0*v0 - 2*k1*m2*n1*s0*v0 - 
  2*k1*m1*n2*s0*v0 + 2*j2*m1*o1*s0*v0 + 
  2*j1*m2*o1*s0*v0 + 2*j1*m1*o2*s0*v0 + 
  2*j2*k1*q1*s0*v0 + 2*j1*k2*q1*s0*v0 - 
  2*i2*o1*q1*s0*v0 - 2*i1*o2*q1*s0*v0 + 
  2*j1*k1*q2*s0*v0 - 2*i1*o1*q2*s0*v0 - 
  2*k2*m1*n0*s1*v0 - 2*k1*m2*n0*s1*v0 - 
  2*k2*m0*n1*s1*v0 - 2*k0*m2*n1*s1*v0 - 
  2*k1*m0*n2*s1*v0 - 2*k0*m1*n2*s1*v0 + 
  2*j2*m1*o0*s1*v0 + 2*j1*m2*o0*s1*v0 + 
  2*j2*m0*o1*s1*v0 + 2*j0*m2*o1*s1*v0 + 
  2*j1*m0*o2*s1*v0 + 2*j0*m1*o2*s1*v0 + 
  2*j2*k1*q0*s1*v0 + 2*j1*k2*q0*s1*v0 - 
  2*i2*o1*q0*s1*v0 - 2*i1*o2*q0*s1*v0 + 
  2*j2*k0*q1*s1*v0 + 2*j0*k2*q1*s1*v0 - 
  2*i2*o0*q1*s1*v0 - 2*i0*o2*q1*s1*v0 + 
  2*j1*k0*q2*s1*v0 + 2*j0*k1*q2*s1*v0 - 
  2*i1*o0*q2*s1*v0 - 2*i0*o1*q2*s1*v0 - 
  2*k1*m1*n0*s2*v0 - 2*k1*m0*n1*s2*v0 - 
  2*k0*m1*n1*s2*v0 + 2*j1*m1*o0*s2*v0 + 
  2*j1*m0*o1*s2*v0 + 2*j0*m1*o1*s2*v0 + 
  2*j1*k1*q0*s2*v0 - 2*i1*o1*q0*s2*v0 + 
  2*j1*k0*q1*s2*v0 + 2*j0*k1*q1*s2*v0 - 
  2*i1*o0*q1*s2*v0 - 2*i0*o1*q1*s2*v0 - 
  2*k2*l1*n1*t0*v0 - 2*k1*l2*n1*t0*v0 - 
  2*k1*l1*n2*t0*v0 + 2*j2*l1*o1*t0*v0 + 
  2*j1*l2*o1*t0*v0 + 2*j1*l1*o2*t0*v0 + 
  2*j2*k1*p1*t0*v0 + 2*j1*k2*p1*t0*v0 - 
  2*i2*o1*p1*t0*v0 - 2*i1*o2*p1*t0*v0 + 
  2*j1*k1*p2*t0*v0 - 2*i1*o1*p2*t0*v0 - 
  4*j1*j2*s1*t0*v0 + 2*i2*n1*s1*t0*v0 + 
  2*i1*n2*s1*t0*v0 - 2*j1*j1*s2*t0*v0 + 
  2*i1*n1*s2*t0*v0 - 2*k2*l1*n0*t1*v0 - 
  2*k1*l2*n0*t1*v0 - 2*k2*l0*n1*t1*v0 - 
  2*k0*l2*n1*t1*v0 - 2*k1*l0*n2*t1*v0 - 
  2*k0*l1*n2*t1*v0 + 2*j2*l1*o0*t1*v0 + 
  2*j1*l2*o0*t1*v0 + 2*j2*l0*o1*t1*v0 + 
  2*j0*l2*o1*t1*v0 + 2*j1*l0*o2*t1*v0 + 
  2*j0*l1*o2*t1*v0 + 2*j2*k1*p0*t1*v0 + 
  2*j1*k2*p0*t1*v0 - 2*i2*o1*p0*t1*v0 - 
  2*i1*o2*p0*t1*v0 + 2*j2*k0*p1*t1*v0 + 
  2*j0*k2*p1*t1*v0 - 2*i2*o0*p1*t1*v0 - 
  2*i0*o2*p1*t1*v0 + 2*j1*k0*p2*t1*v0 + 
  2*j0*k1*p2*t1*v0 - 2*i1*o0*p2*t1*v0 - 
  2*i0*o1*p2*t1*v0 - 4*j1*j2*s0*t1*v0 + 
  2*i2*n1*s0*t1*v0 + 2*i1*n2*s0*t1*v0 - 
  4*j0*j2*s1*t1*v0 + 2*i2*n0*s1*t1*v0 + 
  2*i0*n2*s1*t1*v0 - 4*j0*j1*s2*t1*v0 + 
  2*i1*n0*s2*t1*v0 + 2*i0*n1*s2*t1*v0 - 
  2*k1*l1*n0*t2*v0 - 2*k1*l0*n1*t2*v0 - 
  2*k0*l1*n1*t2*v0 + 2*j1*l1*o0*t2*v0 + 
  2*j1*l0*o1*t2*v0 + 2*j0*l1*o1*t2*v0 + 
  2*j1*k1*p0*t2*v0 - 2*i1*o1*p0*t2*v0 + 
  2*j1*k0*p1*t2*v0 + 2*j0*k1*p1*t2*v0 - 
  2*i1*o0*p1*t2*v0 - 2*i0*o1*p1*t2*v0 - 
  2*j1*j1*s0*t2*v0 + 2*i1*n1*s0*t2*v0 - 
  4*j0*j1*s1*t2*v0 + 2*i1*n0*s1*t2*v0 + 
  2*i0*n1*s1*t2*v0 + 2*k1*k2*n1*v0*v0 + k1*k1*n2*v0*v0 - 
  2*j2*k1*o1*v0*v0 - 2*j1*k2*o1*v0*v0 + i2*o1*o1*v0*v0 - 
  2*j1*k1*o2*v0*v0 + 2*i1*o1*o2*v0*v0 + 
  2*j1*j2*r1*v0*v0 - i2*n1*r1*v0*v0 - i1*n2*r1*v0*v0 + 
  j1*j1*r2*v0*v0 - i1*n1*r2*v0*v0 - 2*l2*m1*o0*o0*v1 - 
  2*l1*m2*o0*o0*v1 - 4*l2*m0*o0*o1*v1 - 
  4*l0*m2*o0*o1*v1 - 4*l1*m0*o0*o2*v1 - 
  4*l0*m1*o0*o2*v1 - 4*l0*m0*o1*o2*v1 + 
  2*k2*m1*o0*p0*v1 + 2*k1*m2*o0*p0*v1 + 
  2*k2*m0*o1*p0*v1 + 2*k0*m2*o1*p0*v1 + 
  2*k1*m0*o2*p0*v1 + 2*k0*m1*o2*p0*v1 + 
  2*k2*m0*o0*p1*v1 + 2*k0*m2*o0*p1*v1 + 
  2*k0*m0*o2*p1*v1 + 2*k1*m0*o0*p2*v1 + 
  2*k0*m1*o0*p2*v1 + 2*k0*m0*o1*p2*v1 + 
  2*k2*l1*o0*q0*v1 + 2*k1*l2*o0*q0*v1 + 
  2*k2*l0*o1*q0*v1 + 2*k0*l2*o1*q0*v1 + 
  2*k1*l0*o2*q0*v1 + 2*k0*l1*o2*q0*v1 - 
  4*k1*k2*p0*q0*v1 - 4*k0*k2*p1*q0*v1 - 
  4*k0*k1*p2*q0*v1 + 2*k2*l0*o0*q1*v1 + 
  2*k0*l2*o0*q1*v1 + 2*k0*l0*o2*q1*v1 - 
  4*k0*k2*p0*q1*v1 - 2*k0*k0*p2*q1*v1 + 
  2*k1*l0*o0*q2*v1 + 2*k0*l1*o0*q2*v1 + 
  2*k0*l0*o1*q2*v1 - 4*k0*k1*p0*q2*v1 - 
  2*k0*k0*p1*q2*v1 + 2*l2*m1*n0*r0*v1 + 
  2*l1*m2*n0*r0*v1 + 2*l2*m0*n1*r0*v1 + 
  2*l0*m2*n1*r0*v1 + 2*l1*m0*n2*r0*v1 + 
  2*l0*m1*n2*r0*v1 - 2*j2*m1*p0*r0*v1 - 
  2*j1*m2*p0*r0*v1 - 2*j2*m0*p1*r0*v1 - 
  2*j0*m2*p1*r0*v1 - 2*j1*m0*p2*r0*v1 - 
  2*j0*m1*p2*r0*v1 - 2*j2*l1*q0*r0*v1 - 
  2*j1*l2*q0*r0*v1 + 2*i2*p1*q0*r0*v1 + 
  2*i1*p2*q0*r0*v1 - 2*j2*l0*q1*r0*v1 - 
  2*j0*l2*q1*r0*v1 + 2*i2*p0*q1*r0*v1 + 
  2*i0*p2*q1*r0*v1 - 2*j1*l0*q2*r0*v1 - 
  2*j0*l1*q2*r0*v1 + 2*i1*p0*q2*r0*v1 + 
  2*i0*p1*q2*r0*v1 + 2*l2*m0*n0*r1*v1 + 
  2*l0*m2*n0*r1*v1 + 2*l0*m0*n2*r1*v1 - 
  2*j2*m0*p0*r1*v1 - 2*j0*m2*p0*r1*v1 - 
  2*j0*m0*p2*r1*v1 - 2*j2*l0*q0*r1*v1 - 
  2*j0*l2*q0*r1*v1 + 2*i2*p0*q0*r1*v1 + 
  2*i0*p2*q0*r1*v1 - 2*j0*l0*q2*r1*v1 + 
  2*i0*p0*q2*r1*v1 + 2*l1*m0*n0*r2*v1 + 
  2*l0*m1*n0*r2*v1 + 2*l0*m0*n1*r2*v1 - 
  2*j1*m0*p0*r2*v1 - 2*j0*m1*p0*r2*v1 - 
  2*j0*m0*p1*r2*v1 - 2*j1*l0*q0*r2*v1 - 
  2*j0*l1*q0*r2*v1 + 2*i1*p0*q0*r2*v1 + 
  2*i0*p1*q0*r2*v1 - 2*j0*l0*q1*r2*v1 + 
  2*i0*p0*q1*r2*v1 - 2*k2*m1*n0*s0*v1 - 
  2*k1*m2*n0*s0*v1 - 2*k2*m0*n1*s0*v1 - 
  2*k0*m2*n1*s0*v1 - 2*k1*m0*n2*s0*v1 - 
  2*k0*m1*n2*s0*v1 + 2*j2*m1*o0*s0*v1 + 
  2*j1*m2*o0*s0*v1 + 2*j2*m0*o1*s0*v1 + 
  2*j0*m2*o1*s0*v1 + 2*j1*m0*o2*s0*v1 + 
  2*j0*m1*o2*s0*v1 + 2*j2*k1*q0*s0*v1 + 
  2*j1*k2*q0*s0*v1 - 2*i2*o1*q0*s0*v1 - 
  2*i1*o2*q0*s0*v1 + 2*j2*k0*q1*s0*v1 + 
  2*j0*k2*q1*s0*v1 - 2*i2*o0*q1*s0*v1 - 
  2*i0*o2*q1*s0*v1 + 2*j1*k0*q2*s0*v1 + 
  2*j0*k1*q2*s0*v1 - 2*i1*o0*q2*s0*v1 - 
  2*i0*o1*q2*s0*v1 - 2*k2*m0*n0*s1*v1 - 
  2*k0*m2*n0*s1*v1 - 2*k0*m0*n2*s1*v1 + 
  2*j2*m0*o0*s1*v1 + 2*j0*m2*o0*s1*v1 + 
  2*j0*m0*o2*s1*v1 + 2*j2*k0*q0*s1*v1 + 
  2*j0*k2*q0*s1*v1 - 2*i2*o0*q0*s1*v1 - 
  2*i0*o2*q0*s1*v1 + 2*j0*k0*q2*s1*v1 - 
  2*i0*o0*q2*s1*v1 - 2*k1*m0*n0*s2*v1 - 
  2*k0*m1*n0*s2*v1 - 2*k0*m0*n1*s2*v1 + 
  2*j1*m0*o0*s2*v1 + 2*j0*m1*o0*s2*v1 + 
  2*j0*m0*o1*s2*v1 + 2*j1*k0*q0*s2*v1 + 
  2*j0*k1*q0*s2*v1 - 2*i1*o0*q0*s2*v1 - 
  2*i0*o1*q0*s2*v1 + 2*j0*k0*q1*s2*v1 - 
  2*i0*o0*q1*s2*v1 - 2*k2*l1*n0*t0*v1 - 
  2*k1*l2*n0*t0*v1 - 2*k2*l0*n1*t0*v1 - 
  2*k0*l2*n1*t0*v1 - 2*k1*l0*n2*t0*v1 - 
  2*k0*l1*n2*t0*v1 + 2*j2*l1*o0*t0*v1 + 
  2*j1*l2*o0*t0*v1 + 2*j2*l0*o1*t0*v1 + 
  2*j0*l2*o1*t0*v1 + 2*j1*l0*o2*t0*v1 + 
  2*j0*l1*o2*t0*v1 + 2*j2*k1*p0*t0*v1 + 
  2*j1*k2*p0*t0*v1 - 2*i2*o1*p0*t0*v1 - 
  2*i1*o2*p0*t0*v1 + 2*j2*k0*p1*t0*v1 + 
  2*j0*k2*p1*t0*v1 - 2*i2*o0*p1*t0*v1 - 
  2*i0*o2*p1*t0*v1 + 2*j1*k0*p2*t0*v1 + 
  2*j0*k1*p2*t0*v1 - 2*i1*o0*p2*t0*v1 - 
  2*i0*o1*p2*t0*v1 - 4*j1*j2*s0*t0*v1 + 
  2*i2*n1*s0*t0*v1 + 2*i1*n2*s0*t0*v1 - 
  4*j0*j2*s1*t0*v1 + 2*i2*n0*s1*t0*v1 + 
  2*i0*n2*s1*t0*v1 - 4*j0*j1*s2*t0*v1 + 
  2*i1*n0*s2*t0*v1 + 2*i0*n1*s2*t0*v1 - 
  2*k2*l0*n0*t1*v1 - 2*k0*l2*n0*t1*v1 - 
  2*k0*l0*n2*t1*v1 + 2*j2*l0*o0*t1*v1 + 
  2*j0*l2*o0*t1*v1 + 2*j0*l0*o2*t1*v1 + 
  2*j2*k0*p0*t1*v1 + 2*j0*k2*p0*t1*v1 - 
  2*i2*o0*p0*t1*v1 - 2*i0*o2*p0*t1*v1 + 
  2*j0*k0*p2*t1*v1 - 2*i0*o0*p2*t1*v1 - 
  4*j0*j2*s0*t1*v1 + 2*i2*n0*s0*t1*v1 + 
  2*i0*n2*s0*t1*v1 - 2*j0*j0*s2*t1*v1 + 
  2*i0*n0*s2*t1*v1 - 2*k1*l0*n0*t2*v1 - 
  2*k0*l1*n0*t2*v1 - 2*k0*l0*n1*t2*v1 + 
  2*j1*l0*o0*t2*v1 + 2*j0*l1*o0*t2*v1 + 
  2*j0*l0*o1*t2*v1 + 2*j1*k0*p0*t2*v1 + 
  2*j0*k1*p0*t2*v1 - 2*i1*o0*p0*t2*v1 - 
  2*i0*o1*p0*t2*v1 + 2*j0*k0*p1*t2*v1 - 
  2*i0*o0*p1*t2*v1 - 4*j0*j1*s0*t2*v1 + 
  2*i1*n0*s0*t2*v1 + 2*i0*n1*s0*t2*v1 - 
  2*j0*j0*s1*t2*v1 + 2*i0*n0*s1*t2*v1 + 
  4*k1*k2*n0*v0*v1 + 4*k0*k2*n1*v0*v1 + 
  4*k0*k1*n2*v0*v1 - 4*j2*k1*o0*v0*v1 - 
  4*j1*k2*o0*v0*v1 - 4*j2*k0*o1*v0*v1 - 
  4*j0*k2*o1*v0*v1 + 4*i2*o0*o1*v0*v1 - 
  4*j1*k0*o2*v0*v1 - 4*j0*k1*o2*v0*v1 + 
  4*i1*o0*o2*v0*v1 + 4*i0*o1*o2*v0*v1 + 
  4*j1*j2*r0*v0*v1 - 2*i2*n1*r0*v0*v1 - 
  2*i1*n2*r0*v0*v1 + 4*j0*j2*r1*v0*v1 - 
  2*i2*n0*r1*v0*v1 - 2*i0*n2*r1*v0*v1 + 
  4*j0*j1*r2*v0*v1 - 2*i1*n0*r2*v0*v1 - 
  2*i0*n1*r2*v0*v1 + 2*k0*k2*n0*v1*v1 + k0*k0*n2*v1*v1 - 
  2*j2*k0*o0*v1*v1 - 2*j0*k2*o0*v1*v1 + i2*o0*o0*v1*v1 - 
  2*j0*k0*o2*v1*v1 + 2*i0*o0*o2*v1*v1 + 
  2*j0*j2*r0*v1*v1 - i2*n0*r0*v1*v1 - i0*n2*r0*v1*v1 + 
  j0*j0*r2*v1*v1 - i0*n0*r2*v1*v1 - 2*l1*m1*o0*o0*v2 - 
  4*l1*m0*o0*o1*v2 - 4*l0*m1*o0*o1*v2 - 
  2*l0*m0*o1*o1*v2 + 2*k1*m1*o0*p0*v2 + 
  2*k1*m0*o1*p0*v2 + 2*k0*m1*o1*p0*v2 + 
  2*k1*m0*o0*p1*v2 + 2*k0*m1*o0*p1*v2 + 
  2*k0*m0*o1*p1*v2 + 2*k1*l1*o0*q0*v2 + 
  2*k1*l0*o1*q0*v2 + 2*k0*l1*o1*q0*v2 - 
  2*k1*k1*p0*q0*v2 - 4*k0*k1*p1*q0*v2 + 
  2*k1*l0*o0*q1*v2 + 2*k0*l1*o0*q1*v2 + 
  2*k0*l0*o1*q1*v2 - 4*k0*k1*p0*q1*v2 - 
  2*k0*k0*p1*q1*v2 + 2*l1*m1*n0*r0*v2 + 
  2*l1*m0*n1*r0*v2 + 2*l0*m1*n1*r0*v2 - 
  2*j1*m1*p0*r0*v2 - 2*j1*m0*p1*r0*v2 - 
  2*j0*m1*p1*r0*v2 - 2*j1*l1*q0*r0*v2 + 
  2*i1*p1*q0*r0*v2 - 2*j1*l0*q1*r0*v2 - 
  2*j0*l1*q1*r0*v2 + 2*i1*p0*q1*r0*v2 + 
  2*i0*p1*q1*r0*v2 + 2*l1*m0*n0*r1*v2 + 
  2*l0*m1*n0*r1*v2 + 2*l0*m0*n1*r1*v2 - 
  2*j1*m0*p0*r1*v2 - 2*j0*m1*p0*r1*v2 - 
  2*j0*m0*p1*r1*v2 - 2*j1*l0*q0*r1*v2 - 
  2*j0*l1*q0*r1*v2 + 2*i1*p0*q0*r1*v2 + 
  2*i0*p1*q0*r1*v2 - 2*j0*l0*q1*r1*v2 + 
  2*i0*p0*q1*r1*v2 - 2*k1*m1*n0*s0*v2 - 
  2*k1*m0*n1*s0*v2 - 2*k0*m1*n1*s0*v2 + 
  2*j1*m1*o0*s0*v2 + 2*j1*m0*o1*s0*v2 + 
  2*j0*m1*o1*s0*v2 + 2*j1*k1*q0*s0*v2 - 
  2*i1*o1*q0*s0*v2 + 2*j1*k0*q1*s0*v2 + 
  2*j0*k1*q1*s0*v2 - 2*i1*o0*q1*s0*v2 - 
  2*i0*o1*q1*s0*v2 - 2*k1*m0*n0*s1*v2 - 
  2*k0*m1*n0*s1*v2 - 2*k0*m0*n1*s1*v2 + 
  2*j1*m0*o0*s1*v2 + 2*j0*m1*o0*s1*v2 + 
  2*j0*m0*o1*s1*v2 + 2*j1*k0*q0*s1*v2 + 
  2*j0*k1*q0*s1*v2 - 2*i1*o0*q0*s1*v2 - 
  2*i0*o1*q0*s1*v2 + 2*j0*k0*q1*s1*v2 - 
  2*i0*o0*q1*s1*v2 - 2*k1*l1*n0*t0*v2 - 
  2*k1*l0*n1*t0*v2 - 2*k0*l1*n1*t0*v2 + 
  2*j1*l1*o0*t0*v2 + 2*j1*l0*o1*t0*v2 + 
  2*j0*l1*o1*t0*v2 + 2*j1*k1*p0*t0*v2 - 
  2*i1*o1*p0*t0*v2 + 2*j1*k0*p1*t0*v2 + 
  2*j0*k1*p1*t0*v2 - 2*i1*o0*p1*t0*v2 - 
  2*i0*o1*p1*t0*v2 - 2*j1*j1*s0*t0*v2 + 
  2*i1*n1*s0*t0*v2 - 4*j0*j1*s1*t0*v2 + 
  2*i1*n0*s1*t0*v2 + 2*i0*n1*s1*t0*v2 - 
  2*k1*l0*n0*t1*v2 - 2*k0*l1*n0*t1*v2 - 
  2*k0*l0*n1*t1*v2 + 2*j1*l0*o0*t1*v2 + 
  2*j0*l1*o0*t1*v2 + 2*j0*l0*o1*t1*v2 + 
  2*j1*k0*p0*t1*v2 + 2*j0*k1*p0*t1*v2 - 
  2*i1*o0*p0*t1*v2 - 2*i0*o1*p0*t1*v2 + 
  2*j0*k0*p1*t1*v2 - 2*i0*o0*p1*t1*v2 - 
  4*j0*j1*s0*t1*v2 + 2*i1*n0*s0*t1*v2 + 
  2*i0*n1*s0*t1*v2 - 2*j0*j0*s1*t1*v2 + 
  2*i0*n0*s1*t1*v2 + 2*k1*k1*n0*v0*v2 + 
  4*k0*k1*n1*v0*v2 - 4*j1*k1*o0*v0*v2 - 
  4*j1*k0*o1*v0*v2 - 4*j0*k1*o1*v0*v2 + 
  4*i1*o0*o1*v0*v2 + 2*i0*o1*o1*v0*v2 + 
  2*j1*j1*r0*v0*v2 - 2*i1*n1*r0*v0*v2 + 
  4*j0*j1*r1*v0*v2 - 2*i1*n0*r1*v0*v2 - 
  2*i0*n1*r1*v0*v2 + 4*k0*k1*n0*v1*v2 + 
  2*k0*k0*n1*v1*v2 - 4*j1*k0*o0*v1*v2 - 
  4*j0*k1*o0*v1*v2 + 2*i1*o0*o0*v1*v2 - 
  4*j0*k0*o1*v1*v2 + 4*i0*o0*o1*v1*v2 + 
  4*j0*j1*r0*v1*v2 - 2*i1*n0*r0*v1*v2 - 
  2*i0*n1*r0*v1*v2 + 2*j0*j0*r1*v1*v2 - 
  2*i0*n0*r1*v1*v2 + 4*l1*l2*o0*o1*w0 + 
  2*l0*l2*o1*o1*w0 + 2*l1*l1*o0*o2*w0 + 
  4*l0*l1*o1*o2*w0 - 2*k2*l1*o1*p0*w0 - 
  2*k1*l2*o1*p0*w0 - 2*k1*l1*o2*p0*w0 - 
  2*k2*l1*o0*p1*w0 - 2*k1*l2*o0*p1*w0 - 
  2*k2*l0*o1*p1*w0 - 2*k0*l2*o1*p1*w0 - 
  2*k1*l0*o2*p1*w0 - 2*k0*l1*o2*p1*w0 + 
  4*k1*k2*p0*p1*w0 + 2*k0*k2*p1*p1*w0 - 
  2*k1*l1*o0*p2*w0 - 2*k1*l0*o1*p2*w0 - 
  2*k0*l1*o1*p2*w0 + 2*k1*k1*p0*p2*w0 + 
  4*k0*k1*p1*p2*w0 - 2*l1*l2*n1*r0*w0 - 
  l1*l1*n2*r0*w0 + 2*j2*l1*p1*r0*w0 + 
  2*j1*l2*p1*r0*w0 - i2*p1*p1*r0*w0 + 
  2*j1*l1*p2*r0*w0 - 2*i1*p1*p2*r0*w0 - 
  2*l1*l2*n0*r1*w0 - 2*l0*l2*n1*r1*w0 - 
  2*l0*l1*n2*r1*w0 + 2*j2*l1*p0*r1*w0 + 
  2*j1*l2*p0*r1*w0 + 2*j2*l0*p1*r1*w0 + 
  2*j0*l2*p1*r1*w0 - 2*i2*p0*p1*r1*w0 + 
  2*j1*l0*p2*r1*w0 + 2*j0*l1*p2*r1*w0 - 
  2*i1*p0*p2*r1*w0 - 2*i0*p1*p2*r1*w0 - 
  l1*l1*n0*r2*w0 - 2*l0*l1*n1*r2*w0 + 
  2*j1*l1*p0*r2*w0 + 2*j1*l0*p1*r2*w0 + 
  2*j0*l1*p1*r2*w0 - 2*i1*p0*p1*r2*w0 - 
  i0*p1*p1*r2*w0 + 2*k2*l1*n1*s0*w0 + 
  2*k1*l2*n1*s0*w0 + 2*k1*l1*n2*s0*w0 - 
  2*j2*l1*o1*s0*w0 - 2*j1*l2*o1*s0*w0 - 
  2*j1*l1*o2*s0*w0 - 2*j2*k1*p1*s0*w0 - 
  2*j1*k2*p1*s0*w0 + 2*i2*o1*p1*s0*w0 + 
  2*i1*o2*p1*s0*w0 - 2*j1*k1*p2*s0*w0 + 
  2*i1*o1*p2*s0*w0 + 2*k2*l1*n0*s1*w0 + 
  2*k1*l2*n0*s1*w0 + 2*k2*l0*n1*s1*w0 + 
  2*k0*l2*n1*s1*w0 + 2*k1*l0*n2*s1*w0 + 
  2*k0*l1*n2*s1*w0 - 2*j2*l1*o0*s1*w0 - 
  2*j1*l2*o0*s1*w0 - 2*j2*l0*o1*s1*w0 - 
  2*j0*l2*o1*s1*w0 - 2*j1*l0*o2*s1*w0 - 
  2*j0*l1*o2*s1*w0 - 2*j2*k1*p0*s1*w0 - 
  2*j1*k2*p0*s1*w0 + 2*i2*o1*p0*s1*w0 + 
  2*i1*o2*p0*s1*w0 - 2*j2*k0*p1*s1*w0 - 
  2*j0*k2*p1*s1*w0 + 2*i2*o0*p1*s1*w0 + 
  2*i0*o2*p1*s1*w0 - 2*j1*k0*p2*s1*w0 - 
  2*j0*k1*p2*s1*w0 + 2*i1*o0*p2*s1*w0 + 
  2*i0*o1*p2*s1*w0 + 4*j1*j2*s0*s1*w0 - 
  2*i2*n1*s0*s1*w0 - 2*i1*n2*s0*s1*w0 + 
  2*j0*j2*s1*s1*w0 - i2*n0*s1*s1*w0 - i0*n2*s1*s1*w0 + 
  2*k1*l1*n0*s2*w0 + 2*k1*l0*n1*s2*w0 + 
  2*k0*l1*n1*s2*w0 - 2*j1*l1*o0*s2*w0 - 
  2*j1*l0*o1*s2*w0 - 2*j0*l1*o1*s2*w0 - 
  2*j1*k1*p0*s2*w0 + 2*i1*o1*p0*s2*w0 - 
  2*j1*k0*p1*s2*w0 - 2*j0*k1*p1*s2*w0 + 
  2*i1*o0*p1*s2*w0 + 2*i0*o1*p1*s2*w0 + 
  2*j1*j1*s0*s2*w0 - 2*i1*n1*s0*s2*w0 + 
  4*j0*j1*s1*s2*w0 - 2*i1*n0*s1*s2*w0 - 
  2*i0*n1*s1*s2*w0 - 2*k1*k2*n1*u0*w0 - 
  k1*k1*n2*u0*w0 + 2*j2*k1*o1*u0*w0 + 
  2*j1*k2*o1*u0*w0 - i2*o1*o1*u0*w0 + 
  2*j1*k1*o2*u0*w0 - 2*i1*o1*o2*u0*w0 - 
  2*j1*j2*r1*u0*w0 + i2*n1*r1*u0*w0 + 
  i1*n2*r1*u0*w0 - j1*j1*r2*u0*w0 + i1*n1*r2*u0*w0 - 
  2*k1*k2*n0*u1*w0 - 2*k0*k2*n1*u1*w0 - 
  2*k0*k1*n2*u1*w0 + 2*j2*k1*o0*u1*w0 + 
  2*j1*k2*o0*u1*w0 + 2*j2*k0*o1*u1*w0 + 
  2*j0*k2*o1*u1*w0 - 2*i2*o0*o1*u1*w0 + 
  2*j1*k0*o2*u1*w0 + 2*j0*k1*o2*u1*w0 - 
  2*i1*o0*o2*u1*w0 - 2*i0*o1*o2*u1*w0 - 
  2*j1*j2*r0*u1*w0 + i2*n1*r0*u1*w0 + 
  i1*n2*r0*u1*w0 - 2*j0*j2*r1*u1*w0 + 
  i2*n0*r1*u1*w0 + i0*n2*r1*u1*w0 - 
  2*j0*j1*r2*u1*w0 + i1*n0*r2*u1*w0 + 
  i0*n1*r2*u1*w0 - k1*k1*n0*u2*w0 - 2*k0*k1*n1*u2*w0 + 
  2*j1*k1*o0*u2*w0 + 2*j1*k0*o1*u2*w0 + 
  2*j0*k1*o1*u2*w0 - 2*i1*o0*o1*u2*w0 - 
  i0*o1*o1*u2*w0 - j1*j1*r0*u2*w0 + i1*n1*r0*u2*w0 - 
  2*j0*j1*r1*u2*w0 + i1*n0*r1*u2*w0 + 
  i0*n1*r1*u2*w0 + 2*l1*l2*o0*o0*w1 + 
  4*l0*l2*o0*o1*w1 + 4*l0*l1*o0*o2*w1 + 
  2*l0*l0*o1*o2*w1 - 2*k2*l1*o0*p0*w1 - 
  2*k1*l2*o0*p0*w1 - 2*k2*l0*o1*p0*w1 - 
  2*k0*l2*o1*p0*w1 - 2*k1*l0*o2*p0*w1 - 
  2*k0*l1*o2*p0*w1 + 2*k1*k2*p0*p0*w1 - 
  2*k2*l0*o0*p1*w1 - 2*k0*l2*o0*p1*w1 - 
  2*k0*l0*o2*p1*w1 + 4*k0*k2*p0*p1*w1 - 
  2*k1*l0*o0*p2*w1 - 2*k0*l1*o0*p2*w1 - 
  2*k0*l0*o1*p2*w1 + 4*k0*k1*p0*p2*w1 + 
  2*k0*k0*p1*p2*w1 - 2*l1*l2*n0*r0*w1 - 
  2*l0*l2*n1*r0*w1 - 2*l0*l1*n2*r0*w1 + 
  2*j2*l1*p0*r0*w1 + 2*j1*l2*p0*r0*w1 + 
  2*j2*l0*p1*r0*w1 + 2*j0*l2*p1*r0*w1 - 
  2*i2*p0*p1*r0*w1 + 2*j1*l0*p2*r0*w1 + 
  2*j0*l1*p2*r0*w1 - 2*i1*p0*p2*r0*w1 - 
  2*i0*p1*p2*r0*w1 - 2*l0*l2*n0*r1*w1 - 
  l0*l0*n2*r1*w1 + 2*j2*l0*p0*r1*w1 + 
  2*j0*l2*p0*r1*w1 - i2*p0*p0*r1*w1 + 
  2*j0*l0*p2*r1*w1 - 2*i0*p0*p2*r1*w1 - 
  2*l0*l1*n0*r2*w1 - l0*l0*n1*r2*w1 + 
  2*j1*l0*p0*r2*w1 + 2*j0*l1*p0*r2*w1 - 
  i1*p0*p0*r2*w1 + 2*j0*l0*p1*r2*w1 - 
  2*i0*p0*p1*r2*w1 + 2*k2*l1*n0*s0*w1 + 
  2*k1*l2*n0*s0*w1 + 2*k2*l0*n1*s0*w1 + 
  2*k0*l2*n1*s0*w1 + 2*k1*l0*n2*s0*w1 + 
  2*k0*l1*n2*s0*w1 - 2*j2*l1*o0*s0*w1 - 
  2*j1*l2*o0*s0*w1 - 2*j2*l0*o1*s0*w1 - 
  2*j0*l2*o1*s0*w1 - 2*j1*l0*o2*s0*w1 - 
  2*j0*l1*o2*s0*w1 - 2*j2*k1*p0*s0*w1 - 
  2*j1*k2*p0*s0*w1 + 2*i2*o1*p0*s0*w1 + 
  2*i1*o2*p0*s0*w1 - 2*j2*k0*p1*s0*w1 - 
  2*j0*k2*p1*s0*w1 + 2*i2*o0*p1*s0*w1 + 
  2*i0*o2*p1*s0*w1 - 2*j1*k0*p2*s0*w1 - 
  2*j0*k1*p2*s0*w1 + 2*i1*o0*p2*s0*w1 + 
  2*i0*o1*p2*s0*w1 + 2*j1*j2*s0*s0*w1 - 
  i2*n1*s0*s0*w1 - i1*n2*s0*s0*w1 + 2*k2*l0*n0*s1*w1 + 
  2*k0*l2*n0*s1*w1 + 2*k0*l0*n2*s1*w1 - 
  2*j2*l0*o0*s1*w1 - 2*j0*l2*o0*s1*w1 - 
  2*j0*l0*o2*s1*w1 - 2*j2*k0*p0*s1*w1 - 
  2*j0*k2*p0*s1*w1 + 2*i2*o0*p0*s1*w1 + 
  2*i0*o2*p0*s1*w1 - 2*j0*k0*p2*s1*w1 + 
  2*i0*o0*p2*s1*w1 + 4*j0*j2*s0*s1*w1 - 
  2*i2*n0*s0*s1*w1 - 2*i0*n2*s0*s1*w1 + 
  2*k1*l0*n0*s2*w1 + 2*k0*l1*n0*s2*w1 + 
  2*k0*l0*n1*s2*w1 - 2*j1*l0*o0*s2*w1 - 
  2*j0*l1*o0*s2*w1 - 2*j0*l0*o1*s2*w1 - 
  2*j1*k0*p0*s2*w1 - 2*j0*k1*p0*s2*w1 + 
  2*i1*o0*p0*s2*w1 + 2*i0*o1*p0*s2*w1 - 
  2*j0*k0*p1*s2*w1 + 2*i0*o0*p1*s2*w1 + 
  4*j0*j1*s0*s2*w1 - 2*i1*n0*s0*s2*w1 - 
  2*i0*n1*s0*s2*w1 + 2*j0*j0*s1*s2*w1 - 
  2*i0*n0*s1*s2*w1 - 2*k1*k2*n0*u0*w1 - 
  2*k0*k2*n1*u0*w1 - 2*k0*k1*n2*u0*w1 + 
  2*j2*k1*o0*u0*w1 + 2*j1*k2*o0*u0*w1 + 
  2*j2*k0*o1*u0*w1 + 2*j0*k2*o1*u0*w1 - 
  2*i2*o0*o1*u0*w1 + 2*j1*k0*o2*u0*w1 + 
  2*j0*k1*o2*u0*w1 - 2*i1*o0*o2*u0*w1 - 
  2*i0*o1*o2*u0*w1 - 2*j1*j2*r0*u0*w1 + 
  i2*n1*r0*u0*w1 + i1*n2*r0*u0*w1 - 
  2*j0*j2*r1*u0*w1 + i2*n0*r1*u0*w1 + 
  i0*n2*r1*u0*w1 - 2*j0*j1*r2*u0*w1 + 
  i1*n0*r2*u0*w1 + i0*n1*r2*u0*w1 - 
  2*k0*k2*n0*u1*w1 - k0*k0*n2*u1*w1 + 
  2*j2*k0*o0*u1*w1 + 2*j0*k2*o0*u1*w1 - 
  i2*o0*o0*u1*w1 + 2*j0*k0*o2*u1*w1 - 
  2*i0*o0*o2*u1*w1 - 2*j0*j2*r0*u1*w1 + 
  i2*n0*r0*u1*w1 + i0*n2*r0*u1*w1 - j0*j0*r2*u1*w1 + 
  i0*n0*r2*u1*w1 - 2*k0*k1*n0*u2*w1 - k0*k0*n1*u2*w1 + 
  2*j1*k0*o0*u2*w1 + 2*j0*k1*o0*u2*w1 - 
  i1*o0*o0*u2*w1 + 2*j0*k0*o1*u2*w1 - 
  2*i0*o0*o1*u2*w1 - 2*j0*j1*r0*u2*w1 + 
  i1*n0*r0*u2*w1 + i0*n1*r0*u2*w1 - j0*j0*r1*u2*w1 + 
  i0*n0*r1*u2*w1 + l1*l1*o0*o0*w2 + 4*l0*l1*o0*o1*w2 + 
  l0*l0*o1*o1*w2 - 2*k1*l1*o0*p0*w2 - 
  2*k1*l0*o1*p0*w2 - 2*k0*l1*o1*p0*w2 + 
  k1*k1*p0*p0*w2 - 2*k1*l0*o0*p1*w2 - 
  2*k0*l1*o0*p1*w2 - 2*k0*l0*o1*p1*w2 + 
  4*k0*k1*p0*p1*w2 + k0*k0*p1*p1*w2 - l1*l1*n0*r0*w2 - 
  2*l0*l1*n1*r0*w2 + 2*j1*l1*p0*r0*w2 + 
  2*j1*l0*p1*r0*w2 + 2*j0*l1*p1*r0*w2 - 
  2*i1*p0*p1*r0*w2 - i0*p1*p1*r0*w2 - 
  2*l0*l1*n0*r1*w2 - l0*l0*n1*r1*w2 + 
  2*j1*l0*p0*r1*w2 + 2*j0*l1*p0*r1*w2 - 
  i1*p0*p0*r1*w2 + 2*j0*l0*p1*r1*w2 - 
  2*i0*p0*p1*r1*w2 + 2*k1*l1*n0*s0*w2 + 
  2*k1*l0*n1*s0*w2 + 2*k0*l1*n1*s0*w2 - 
  2*j1*l1*o0*s0*w2 - 2*j1*l0*o1*s0*w2 - 
  2*j0*l1*o1*s0*w2 - 2*j1*k1*p0*s0*w2 + 
  2*i1*o1*p0*s0*w2 - 2*j1*k0*p1*s0*w2 - 
  2*j0*k1*p1*s0*w2 + 2*i1*o0*p1*s0*w2 + 
  2*i0*o1*p1*s0*w2 + j1*j1*s0*s0*w2 - i1*n1*s0*s0*w2 + 
  2*k1*l0*n0*s1*w2 + 2*k0*l1*n0*s1*w2 + 
  2*k0*l0*n1*s1*w2 - 2*j1*l0*o0*s1*w2 - 
  2*j0*l1*o0*s1*w2 - 2*j0*l0*o1*s1*w2 - 
  2*j1*k0*p0*s1*w2 - 2*j0*k1*p0*s1*w2 + 
  2*i1*o0*p0*s1*w2 + 2*i0*o1*p0*s1*w2 - 
  2*j0*k0*p1*s1*w2 + 2*i0*o0*p1*s1*w2 + 
  4*j0*j1*s0*s1*w2 - 2*i1*n0*s0*s1*w2 - 
  2*i0*n1*s0*s1*w2 + j0*j0*s1*s1*w2 - i0*n0*s1*s1*w2 - 
  k1*k1*n0*u0*w2 - 2*k0*k1*n1*u0*w2 + 
  2*j1*k1*o0*u0*w2 + 2*j1*k0*o1*u0*w2 + 
  2*j0*k1*o1*u0*w2 - 2*i1*o0*o1*u0*w2 - 
  i0*o1*o1*u0*w2 - j1*j1*r0*u0*w2 + i1*n1*r0*u0*w2 - 
  2*j0*j1*r1*u0*w2 + i1*n0*r1*u0*w2 + 
  i0*n1*r1*u0*w2 - 2*k0*k1*n0*u1*w2 - k0*k0*n1*u1*w2 + 
  2*j1*k0*o0*u1*w2 + 2*j0*k1*o0*u1*w2 - 
  i1*o0*o0*u1*w2 + 2*j0*k0*o1*u1*w2 - 
  2*i0*o0*o1*u1*w2 - 2*j0*j1*r0*u1*w2 + 
  i1*n0*r0*u1*w2 + i0*n1*r0*u1*w2 - j0*j0*r1*u1*w2 + 
  i0*n0*r1*u1*w2;
    
    k[13] = 2*m1*m2*p1*p1*r0 + 2*m1*m1*p1*p2*r0 - 
  2*l2*m1*p1*q1*r0 - 2*l1*m2*p1*q1*r0 - 2*l1*m1*p2*q1*r0 + 
  2*l1*l2*q1*q1*r0 - 2*l1*m1*p1*q2*r0 + 2*l1*l1*q1*q2*r0 + 
  4*m1*m2*p0*p1*r1 + 2*m0*m2*p1*p1*r1 + 2*m1*m1*p0*p2*r1 + 
  4*m0*m1*p1*p2*r1 - 2*l2*m1*p1*q0*r1 - 2*l1*m2*p1*q0*r1 - 
  2*l1*m1*p2*q0*r1 - 2*l2*m1*p0*q1*r1 - 2*l1*m2*p0*q1*r1 - 
  2*l2*m0*p1*q1*r1 - 2*l0*m2*p1*q1*r1 - 2*l1*m0*p2*q1*r1 - 
  2*l0*m1*p2*q1*r1 + 4*l1*l2*q0*q1*r1 + 2*l0*l2*q1*q1*r1 - 
  2*l1*m1*p0*q2*r1 - 2*l1*m0*p1*q2*r1 - 2*l0*m1*p1*q2*r1 + 
  2*l1*l1*q0*q2*r1 + 4*l0*l1*q1*q2*r1 + 2*m1*m1*p0*p1*r2 + 
  2*m0*m1*p1*p1*r2 - 2*l1*m1*p1*q0*r2 - 2*l1*m1*p0*q1*r2 - 
  2*l1*m0*p1*q1*r2 - 2*l0*m1*p1*q1*r2 + 2*l1*l1*q0*q1*r2 + 
  2*l0*l1*q1*q1*r2 - 4*m1*m2*o1*p1*s0 - 2*m1*m1*o2*p1*s0 - 
  2*m1*m1*o1*p2*s0 + 2*l2*m1*o1*q1*s0 + 2*l1*m2*o1*q1*s0 + 
  2*l1*m1*o2*q1*s0 + 2*k2*m1*p1*q1*s0 + 2*k1*m2*p1*q1*s0 + 
  2*k1*m1*p2*q1*s0 - 2*k2*l1*q1*q1*s0 - 2*k1*l2*q1*q1*s0 + 
  2*l1*m1*o1*q2*s0 + 2*k1*m1*p1*q2*s0 - 4*k1*l1*q1*q2*s0 - 
  4*m1*m2*o1*p0*s1 - 2*m1*m1*o2*p0*s1 - 4*m1*m2*o0*p1*s1 - 
  4*m0*m2*o1*p1*s1 - 4*m0*m1*o2*p1*s1 - 2*m1*m1*o0*p2*s1 - 
  4*m0*m1*o1*p2*s1 + 2*l2*m1*o1*q0*s1 + 2*l1*m2*o1*q0*s1 + 
  2*l1*m1*o2*q0*s1 + 2*k2*m1*p1*q0*s1 + 2*k1*m2*p1*q0*s1 + 
  2*k1*m1*p2*q0*s1 + 2*l2*m1*o0*q1*s1 + 2*l1*m2*o0*q1*s1 + 
  2*l2*m0*o1*q1*s1 + 2*l0*m2*o1*q1*s1 + 2*l1*m0*o2*q1*s1 + 
  2*l0*m1*o2*q1*s1 + 2*k2*m1*p0*q1*s1 + 2*k1*m2*p0*q1*s1 + 
  2*k2*m0*p1*q1*s1 + 2*k0*m2*p1*q1*s1 + 2*k1*m0*p2*q1*s1 + 
  2*k0*m1*p2*q1*s1 - 4*k2*l1*q0*q1*s1 - 4*k1*l2*q0*q1*s1 - 
  2*k2*l0*q1*q1*s1 - 2*k0*l2*q1*q1*s1 + 2*l1*m1*o0*q2*s1 + 
  2*l1*m0*o1*q2*s1 + 2*l0*m1*o1*q2*s1 + 2*k1*m1*p0*q2*s1 + 
  2*k1*m0*p1*q2*s1 + 2*k0*m1*p1*q2*s1 - 4*k1*l1*q0*q2*s1 - 
  4*k1*l0*q1*q2*s1 - 4*k0*l1*q1*q2*s1 + 4*m1*m2*n1*s0*s1 + 
  2*m1*m1*n2*s0*s1 - 4*j2*m1*q1*s0*s1 - 4*j1*m2*q1*s0*s1 + 
  2*i2*q1*q1*s0*s1 - 4*j1*m1*q2*s0*s1 + 4*i1*q1*q2*s0*s1 + 
  2*m1*m2*n0*s1*s1 + 2*m0*m2*n1*s1*s1 + 2*m0*m1*n2*s1*s1 - 
  2*j2*m1*q0*s1*s1 - 2*j1*m2*q0*s1*s1 - 2*j2*m0*q1*s1*s1 - 
  2*j0*m2*q1*s1*s1 + 2*i2*q0*q1*s1*s1 - 2*j1*m0*q2*s1*s1 - 
  2*j0*m1*q2*s1*s1 + 2*i1*q0*q2*s1*s1 + 2*i0*q1*q2*s1*s1 - 
  2*m1*m1*o1*p0*s2 - 2*m1*m1*o0*p1*s2 - 4*m0*m1*o1*p1*s2 + 
  2*l1*m1*o1*q0*s2 + 2*k1*m1*p1*q0*s2 + 2*l1*m1*o0*q1*s2 + 
  2*l1*m0*o1*q1*s2 + 2*l0*m1*o1*q1*s2 + 2*k1*m1*p0*q1*s2 + 
  2*k1*m0*p1*q1*s2 + 2*k0*m1*p1*q1*s2 - 4*k1*l1*q0*q1*s2 - 
  2*k1*l0*q1*q1*s2 - 2*k0*l1*q1*q1*s2 + 2*m1*m1*n1*s0*s2 - 
  4*j1*m1*q1*s0*s2 + 2*i1*q1*q1*s0*s2 + 2*m1*m1*n0*s1*s2 + 
  4*m0*m1*n1*s1*s2 - 4*j1*m1*q0*s1*s2 - 4*j1*m0*q1*s1*s2 - 
  4*j0*m1*q1*s1*s2 + 4*i1*q0*q1*s1*s2 + 2*i0*q1*q1*s1*s2 + 
  2*l2*m1*o1*p1*t0 + 2*l1*m2*o1*p1*t0 + 2*l1*m1*o2*p1*t0 - 
  2*k2*m1*p1*p1*t0 - 2*k1*m2*p1*p1*t0 + 2*l1*m1*o1*p2*t0 - 
  4*k1*m1*p1*p2*t0 - 4*l1*l2*o1*q1*t0 - 2*l1*l1*o2*q1*t0 + 
  2*k2*l1*p1*q1*t0 + 2*k1*l2*p1*q1*t0 + 2*k1*l1*p2*q1*t0 - 
  2*l1*l1*o1*q2*t0 + 2*k1*l1*p1*q2*t0 - 2*l2*m1*n1*s1*t0 - 
  2*l1*m2*n1*s1*t0 - 2*l1*m1*n2*s1*t0 + 2*j2*m1*p1*s1*t0 + 
  2*j1*m2*p1*s1*t0 + 2*j1*m1*p2*s1*t0 + 2*j2*l1*q1*s1*t0 + 
  2*j1*l2*q1*s1*t0 - 2*i2*p1*q1*s1*t0 - 2*i1*p2*q1*s1*t0 + 
  2*j1*l1*q2*s1*t0 - 2*i1*p1*q2*s1*t0 - 2*l1*m1*n1*s2*t0 + 
  2*j1*m1*p1*s2*t0 + 2*j1*l1*q1*s2*t0 - 2*i1*p1*q1*s2*t0 + 
  2*l2*m1*o1*p0*t1 + 2*l1*m2*o1*p0*t1 + 2*l1*m1*o2*p0*t1 + 
  2*l2*m1*o0*p1*t1 + 2*l1*m2*o0*p1*t1 + 2*l2*m0*o1*p1*t1 + 
  2*l0*m2*o1*p1*t1 + 2*l1*m0*o2*p1*t1 + 2*l0*m1*o2*p1*t1 - 
  4*k2*m1*p0*p1*t1 - 4*k1*m2*p0*p1*t1 - 2*k2*m0*p1*p1*t1 - 
  2*k0*m2*p1*p1*t1 + 2*l1*m1*o0*p2*t1 + 2*l1*m0*o1*p2*t1 + 
  2*l0*m1*o1*p2*t1 - 4*k1*m1*p0*p2*t1 - 4*k1*m0*p1*p2*t1 - 
  4*k0*m1*p1*p2*t1 - 4*l1*l2*o1*q0*t1 - 2*l1*l1*o2*q0*t1 + 
  2*k2*l1*p1*q0*t1 + 2*k1*l2*p1*q0*t1 + 2*k1*l1*p2*q0*t1 - 
  4*l1*l2*o0*q1*t1 - 4*l0*l2*o1*q1*t1 - 4*l0*l1*o2*q1*t1 + 
  2*k2*l1*p0*q1*t1 + 2*k1*l2*p0*q1*t1 + 2*k2*l0*p1*q1*t1 + 
  2*k0*l2*p1*q1*t1 + 2*k1*l0*p2*q1*t1 + 2*k0*l1*p2*q1*t1 - 
  2*l1*l1*o0*q2*t1 - 4*l0*l1*o1*q2*t1 + 2*k1*l1*p0*q2*t1 + 
  2*k1*l0*p1*q2*t1 + 2*k0*l1*p1*q2*t1 - 2*l2*m1*n1*s0*t1 - 
  2*l1*m2*n1*s0*t1 - 2*l1*m1*n2*s0*t1 + 2*j2*m1*p1*s0*t1 + 
  2*j1*m2*p1*s0*t1 + 2*j1*m1*p2*s0*t1 + 2*j2*l1*q1*s0*t1 + 
  2*j1*l2*q1*s0*t1 - 2*i2*p1*q1*s0*t1 - 2*i1*p2*q1*s0*t1 + 
  2*j1*l1*q2*s0*t1 - 2*i1*p1*q2*s0*t1 - 2*l2*m1*n0*s1*t1 - 
  2*l1*m2*n0*s1*t1 - 2*l2*m0*n1*s1*t1 - 2*l0*m2*n1*s1*t1 - 
  2*l1*m0*n2*s1*t1 - 2*l0*m1*n2*s1*t1 + 2*j2*m1*p0*s1*t1 + 
  2*j1*m2*p0*s1*t1 + 2*j2*m0*p1*s1*t1 + 2*j0*m2*p1*s1*t1 + 
  2*j1*m0*p2*s1*t1 + 2*j0*m1*p2*s1*t1 + 2*j2*l1*q0*s1*t1 + 
  2*j1*l2*q0*s1*t1 - 2*i2*p1*q0*s1*t1 - 2*i1*p2*q0*s1*t1 + 
  2*j2*l0*q1*s1*t1 + 2*j0*l2*q1*s1*t1 - 2*i2*p0*q1*s1*t1 - 
  2*i0*p2*q1*s1*t1 + 2*j1*l0*q2*s1*t1 + 2*j0*l1*q2*s1*t1 - 
  2*i1*p0*q2*s1*t1 - 2*i0*p1*q2*s1*t1 - 2*l1*m1*n0*s2*t1 - 
  2*l1*m0*n1*s2*t1 - 2*l0*m1*n1*s2*t1 + 2*j1*m1*p0*s2*t1 + 
  2*j1*m0*p1*s2*t1 + 2*j0*m1*p1*s2*t1 + 2*j1*l1*q0*s2*t1 - 
  2*i1*p1*q0*s2*t1 + 2*j1*l0*q1*s2*t1 + 2*j0*l1*q1*s2*t1 - 
  2*i1*p0*q1*s2*t1 - 2*i0*p1*q1*s2*t1 + 4*l1*l2*n1*t0*t1 + 
  2*l1*l1*n2*t0*t1 - 4*j2*l1*p1*t0*t1 - 4*j1*l2*p1*t0*t1 + 
  2*i2*p1*p1*t0*t1 - 4*j1*l1*p2*t0*t1 + 4*i1*p1*p2*t0*t1 + 
  2*l1*l2*n0*t1*t1 + 2*l0*l2*n1*t1*t1 + 2*l0*l1*n2*t1*t1 - 
  2*j2*l1*p0*t1*t1 - 2*j1*l2*p0*t1*t1 - 2*j2*l0*p1*t1*t1 - 
  2*j0*l2*p1*t1*t1 + 2*i2*p0*p1*t1*t1 - 2*j1*l0*p2*t1*t1 - 
  2*j0*l1*p2*t1*t1 + 2*i1*p0*p2*t1*t1 + 2*i0*p1*p2*t1*t1 + 
  2*l1*m1*o1*p0*t2 + 2*l1*m1*o0*p1*t2 + 2*l1*m0*o1*p1*t2 + 
  2*l0*m1*o1*p1*t2 - 4*k1*m1*p0*p1*t2 - 2*k1*m0*p1*p1*t2 - 
  2*k0*m1*p1*p1*t2 - 2*l1*l1*o1*q0*t2 + 2*k1*l1*p1*q0*t2 - 
  2*l1*l1*o0*q1*t2 - 4*l0*l1*o1*q1*t2 + 2*k1*l1*p0*q1*t2 + 
  2*k1*l0*p1*q1*t2 + 2*k0*l1*p1*q1*t2 - 2*l1*m1*n1*s0*t2 + 
  2*j1*m1*p1*s0*t2 + 2*j1*l1*q1*s0*t2 - 2*i1*p1*q1*s0*t2 - 
  2*l1*m1*n0*s1*t2 - 2*l1*m0*n1*s1*t2 - 2*l0*m1*n1*s1*t2 + 
  2*j1*m1*p0*s1*t2 + 2*j1*m0*p1*s1*t2 + 2*j0*m1*p1*s1*t2 + 
  2*j1*l1*q0*s1*t2 - 2*i1*p1*q0*s1*t2 + 2*j1*l0*q1*s1*t2 + 
  2*j0*l1*q1*s1*t2 - 2*i1*p0*q1*s1*t2 - 2*i0*p1*q1*s1*t2 + 
  2*l1*l1*n1*t0*t2 - 4*j1*l1*p1*t0*t2 + 2*i1*p1*p1*t0*t2 + 
  2*l1*l1*n0*t1*t2 + 4*l0*l1*n1*t1*t2 - 4*j1*l1*p0*t1*t2 - 
  4*j1*l0*p1*t1*t2 - 4*j0*l1*p1*t1*t2 + 4*i1*p0*p1*t1*t2 + 
  2*i0*p1*p1*t1*t2 + 2*m1*m2*o1*o1*u0 + 2*m1*m1*o1*o2*u0 - 
  2*k2*m1*o1*q1*u0 - 2*k1*m2*o1*q1*u0 - 2*k1*m1*o2*q1*u0 + 
  2*k1*k2*q1*q1*u0 - 2*k1*m1*o1*q2*u0 + 2*k1*k1*q1*q2*u0 - 
  2*m1*m2*n1*r1*u0 - m1*m1*n2*r1*u0 + 2*j2*m1*q1*r1*u0 + 
  2*j1*m2*q1*r1*u0 - i2*q1*q1*r1*u0 + 2*j1*m1*q2*r1*u0 - 
  2*i1*q1*q2*r1*u0 - m1*m1*n1*r2*u0 + 2*j1*m1*q1*r2*u0 - 
  i1*q1*q1*r2*u0 + 2*k2*m1*n1*t1*u0 + 2*k1*m2*n1*t1*u0 + 
  2*k1*m1*n2*t1*u0 - 2*j2*m1*o1*t1*u0 - 2*j1*m2*o1*t1*u0 - 
  2*j1*m1*o2*t1*u0 - 2*j2*k1*q1*t1*u0 - 2*j1*k2*q1*t1*u0 + 
  2*i2*o1*q1*t1*u0 + 2*i1*o2*q1*t1*u0 - 2*j1*k1*q2*t1*u0 + 
  2*i1*o1*q2*t1*u0 + 2*j1*j2*t1*t1*u0 - i2*n1*t1*t1*u0 - 
  i1*n2*t1*t1*u0 + 2*k1*m1*n1*t2*u0 - 2*j1*m1*o1*t2*u0 - 
  2*j1*k1*q1*t2*u0 + 2*i1*o1*q1*t2*u0 + 2*j1*j1*t1*t2*u0 - 
  2*i1*n1*t1*t2*u0 + 4*m1*m2*o0*o1*u1 + 2*m0*m2*o1*o1*u1 + 
  2*m1*m1*o0*o2*u1 + 4*m0*m1*o1*o2*u1 - 2*k2*m1*o1*q0*u1 - 
  2*k1*m2*o1*q0*u1 - 2*k1*m1*o2*q0*u1 - 2*k2*m1*o0*q1*u1 - 
  2*k1*m2*o0*q1*u1 - 2*k2*m0*o1*q1*u1 - 2*k0*m2*o1*q1*u1 - 
  2*k1*m0*o2*q1*u1 - 2*k0*m1*o2*q1*u1 + 4*k1*k2*q0*q1*u1 + 
  2*k0*k2*q1*q1*u1 - 2*k1*m1*o0*q2*u1 - 2*k1*m0*o1*q2*u1 - 
  2*k0*m1*o1*q2*u1 + 2*k1*k1*q0*q2*u1 + 4*k0*k1*q1*q2*u1 - 
  2*m1*m2*n1*r0*u1 - m1*m1*n2*r0*u1 + 2*j2*m1*q1*r0*u1 + 
  2*j1*m2*q1*r0*u1 - i2*q1*q1*r0*u1 + 2*j1*m1*q2*r0*u1 - 
  2*i1*q1*q2*r0*u1 - 2*m1*m2*n0*r1*u1 - 2*m0*m2*n1*r1*u1 - 
  2*m0*m1*n2*r1*u1 + 2*j2*m1*q0*r1*u1 + 2*j1*m2*q0*r1*u1 + 
  2*j2*m0*q1*r1*u1 + 2*j0*m2*q1*r1*u1 - 2*i2*q0*q1*r1*u1 + 
  2*j1*m0*q2*r1*u1 + 2*j0*m1*q2*r1*u1 - 2*i1*q0*q2*r1*u1 - 
  2*i0*q1*q2*r1*u1 - m1*m1*n0*r2*u1 - 2*m0*m1*n1*r2*u1 + 
  2*j1*m1*q0*r2*u1 + 2*j1*m0*q1*r2*u1 + 2*j0*m1*q1*r2*u1 - 
  2*i1*q0*q1*r2*u1 - i0*q1*q1*r2*u1 + 2*k2*m1*n1*t0*u1 + 
  2*k1*m2*n1*t0*u1 + 2*k1*m1*n2*t0*u1 - 2*j2*m1*o1*t0*u1 - 
  2*j1*m2*o1*t0*u1 - 2*j1*m1*o2*t0*u1 - 2*j2*k1*q1*t0*u1 - 
  2*j1*k2*q1*t0*u1 + 2*i2*o1*q1*t0*u1 + 2*i1*o2*q1*t0*u1 - 
  2*j1*k1*q2*t0*u1 + 2*i1*o1*q2*t0*u1 + 2*k2*m1*n0*t1*u1 + 
  2*k1*m2*n0*t1*u1 + 2*k2*m0*n1*t1*u1 + 2*k0*m2*n1*t1*u1 + 
  2*k1*m0*n2*t1*u1 + 2*k0*m1*n2*t1*u1 - 2*j2*m1*o0*t1*u1 - 
  2*j1*m2*o0*t1*u1 - 2*j2*m0*o1*t1*u1 - 2*j0*m2*o1*t1*u1 - 
  2*j1*m0*o2*t1*u1 - 2*j0*m1*o2*t1*u1 - 2*j2*k1*q0*t1*u1 - 
  2*j1*k2*q0*t1*u1 + 2*i2*o1*q0*t1*u1 + 2*i1*o2*q0*t1*u1 - 
  2*j2*k0*q1*t1*u1 - 2*j0*k2*q1*t1*u1 + 2*i2*o0*q1*t1*u1 + 
  2*i0*o2*q1*t1*u1 - 2*j1*k0*q2*t1*u1 - 2*j0*k1*q2*t1*u1 + 
  2*i1*o0*q2*t1*u1 + 2*i0*o1*q2*t1*u1 + 4*j1*j2*t0*t1*u1 - 
  2*i2*n1*t0*t1*u1 - 2*i1*n2*t0*t1*u1 + 2*j0*j2*t1*t1*u1 - 
  i2*n0*t1*t1*u1 - i0*n2*t1*t1*u1 + 2*k1*m1*n0*t2*u1 + 
  2*k1*m0*n1*t2*u1 + 2*k0*m1*n1*t2*u1 - 2*j1*m1*o0*t2*u1 - 
  2*j1*m0*o1*t2*u1 - 2*j0*m1*o1*t2*u1 - 2*j1*k1*q0*t2*u1 + 
  2*i1*o1*q0*t2*u1 - 2*j1*k0*q1*t2*u1 - 2*j0*k1*q1*t2*u1 + 
  2*i1*o0*q1*t2*u1 + 2*i0*o1*q1*t2*u1 + 2*j1*j1*t0*t2*u1 - 
  2*i1*n1*t0*t2*u1 + 4*j0*j1*t1*t2*u1 - 2*i1*n0*t1*t2*u1 - 
  2*i0*n1*t1*t2*u1 + 2*m1*m1*o0*o1*u2 + 2*m0*m1*o1*o1*u2 - 
  2*k1*m1*o1*q0*u2 - 2*k1*m1*o0*q1*u2 - 2*k1*m0*o1*q1*u2 - 
  2*k0*m1*o1*q1*u2 + 2*k1*k1*q0*q1*u2 + 2*k0*k1*q1*q1*u2 - 
  m1*m1*n1*r0*u2 + 2*j1*m1*q1*r0*u2 - i1*q1*q1*r0*u2 - 
  m1*m1*n0*r1*u2 - 2*m0*m1*n1*r1*u2 + 2*j1*m1*q0*r1*u2 + 
  2*j1*m0*q1*r1*u2 + 2*j0*m1*q1*r1*u2 - 2*i1*q0*q1*r1*u2 - 
  i0*q1*q1*r1*u2 + 2*k1*m1*n1*t0*u2 - 2*j1*m1*o1*t0*u2 - 
  2*j1*k1*q1*t0*u2 + 2*i1*o1*q1*t0*u2 + 2*k1*m1*n0*t1*u2 + 
  2*k1*m0*n1*t1*u2 + 2*k0*m1*n1*t1*u2 - 2*j1*m1*o0*t1*u2 - 
  2*j1*m0*o1*t1*u2 - 2*j0*m1*o1*t1*u2 - 2*j1*k1*q0*t1*u2 + 
  2*i1*o1*q0*t1*u2 - 2*j1*k0*q1*t1*u2 - 2*j0*k1*q1*t1*u2 + 
  2*i1*o0*q1*t1*u2 + 2*i0*o1*q1*t1*u2 + 2*j1*j1*t0*t1*u2 - 
  2*i1*n1*t0*t1*u2 + 2*j0*j1*t1*t1*u2 - i1*n0*t1*t1*u2 - 
  i0*n1*t1*t1*u2 - 2*l2*m1*o1*o1*v0 - 2*l1*m2*o1*o1*v0 - 
  4*l1*m1*o1*o2*v0 + 2*k2*m1*o1*p1*v0 + 2*k1*m2*o1*p1*v0 + 
  2*k1*m1*o2*p1*v0 + 2*k1*m1*o1*p2*v0 + 2*k2*l1*o1*q1*v0 + 
  2*k1*l2*o1*q1*v0 + 2*k1*l1*o2*q1*v0 - 4*k1*k2*p1*q1*v0 - 
  2*k1*k1*p2*q1*v0 + 2*k1*l1*o1*q2*v0 - 2*k1*k1*p1*q2*v0 + 
  2*l2*m1*n1*r1*v0 + 2*l1*m2*n1*r1*v0 + 2*l1*m1*n2*r1*v0 - 
  2*j2*m1*p1*r1*v0 - 2*j1*m2*p1*r1*v0 - 2*j1*m1*p2*r1*v0 - 
  2*j2*l1*q1*r1*v0 - 2*j1*l2*q1*r1*v0 + 2*i2*p1*q1*r1*v0 + 
  2*i1*p2*q1*r1*v0 - 2*j1*l1*q2*r1*v0 + 2*i1*p1*q2*r1*v0 + 
  2*l1*m1*n1*r2*v0 - 2*j1*m1*p1*r2*v0 - 2*j1*l1*q1*r2*v0 + 
  2*i1*p1*q1*r2*v0 - 2*k2*m1*n1*s1*v0 - 2*k1*m2*n1*s1*v0 - 
  2*k1*m1*n2*s1*v0 + 2*j2*m1*o1*s1*v0 + 2*j1*m2*o1*s1*v0 + 
  2*j1*m1*o2*s1*v0 + 2*j2*k1*q1*s1*v0 + 2*j1*k2*q1*s1*v0 - 
  2*i2*o1*q1*s1*v0 - 2*i1*o2*q1*s1*v0 + 2*j1*k1*q2*s1*v0 - 
  2*i1*o1*q2*s1*v0 - 2*k1*m1*n1*s2*v0 + 2*j1*m1*o1*s2*v0 + 
  2*j1*k1*q1*s2*v0 - 2*i1*o1*q1*s2*v0 - 2*k2*l1*n1*t1*v0 - 
  2*k1*l2*n1*t1*v0 - 2*k1*l1*n2*t1*v0 + 2*j2*l1*o1*t1*v0 + 
  2*j1*l2*o1*t1*v0 + 2*j1*l1*o2*t1*v0 + 2*j2*k1*p1*t1*v0 + 
  2*j1*k2*p1*t1*v0 - 2*i2*o1*p1*t1*v0 - 2*i1*o2*p1*t1*v0 + 
  2*j1*k1*p2*t1*v0 - 2*i1*o1*p2*t1*v0 - 4*j1*j2*s1*t1*v0 + 
  2*i2*n1*s1*t1*v0 + 2*i1*n2*s1*t1*v0 - 2*j1*j1*s2*t1*v0 + 
  2*i1*n1*s2*t1*v0 - 2*k1*l1*n1*t2*v0 + 2*j1*l1*o1*t2*v0 + 
  2*j1*k1*p1*t2*v0 - 2*i1*o1*p1*t2*v0 - 2*j1*j1*s1*t2*v0 + 
  2*i1*n1*s1*t2*v0 - 4*l2*m1*o0*o1*v1 - 4*l1*m2*o0*o1*v1 - 
  2*l2*m0*o1*o1*v1 - 2*l0*m2*o1*o1*v1 - 4*l1*m1*o0*o2*v1 - 
  4*l1*m0*o1*o2*v1 - 4*l0*m1*o1*o2*v1 + 2*k2*m1*o1*p0*v1 + 
  2*k1*m2*o1*p0*v1 + 2*k1*m1*o2*p0*v1 + 2*k2*m1*o0*p1*v1 + 
  2*k1*m2*o0*p1*v1 + 2*k2*m0*o1*p1*v1 + 2*k0*m2*o1*p1*v1 + 
  2*k1*m0*o2*p1*v1 + 2*k0*m1*o2*p1*v1 + 2*k1*m1*o0*p2*v1 + 
  2*k1*m0*o1*p2*v1 + 2*k0*m1*o1*p2*v1 + 2*k2*l1*o1*q0*v1 + 
  2*k1*l2*o1*q0*v1 + 2*k1*l1*o2*q0*v1 - 4*k1*k2*p1*q0*v1 - 
  2*k1*k1*p2*q0*v1 + 2*k2*l1*o0*q1*v1 + 2*k1*l2*o0*q1*v1 + 
  2*k2*l0*o1*q1*v1 + 2*k0*l2*o1*q1*v1 + 2*k1*l0*o2*q1*v1 + 
  2*k0*l1*o2*q1*v1 - 4*k1*k2*p0*q1*v1 - 4*k0*k2*p1*q1*v1 - 
  4*k0*k1*p2*q1*v1 + 2*k1*l1*o0*q2*v1 + 2*k1*l0*o1*q2*v1 + 
  2*k0*l1*o1*q2*v1 - 2*k1*k1*p0*q2*v1 - 4*k0*k1*p1*q2*v1 + 
  2*l2*m1*n1*r0*v1 + 2*l1*m2*n1*r0*v1 + 2*l1*m1*n2*r0*v1 - 
  2*j2*m1*p1*r0*v1 - 2*j1*m2*p1*r0*v1 - 2*j1*m1*p2*r0*v1 - 
  2*j2*l1*q1*r0*v1 - 2*j1*l2*q1*r0*v1 + 2*i2*p1*q1*r0*v1 + 
  2*i1*p2*q1*r0*v1 - 2*j1*l1*q2*r0*v1 + 2*i1*p1*q2*r0*v1 + 
  2*l2*m1*n0*r1*v1 + 2*l1*m2*n0*r1*v1 + 2*l2*m0*n1*r1*v1 + 
  2*l0*m2*n1*r1*v1 + 2*l1*m0*n2*r1*v1 + 2*l0*m1*n2*r1*v1 - 
  2*j2*m1*p0*r1*v1 - 2*j1*m2*p0*r1*v1 - 2*j2*m0*p1*r1*v1 - 
  2*j0*m2*p1*r1*v1 - 2*j1*m0*p2*r1*v1 - 2*j0*m1*p2*r1*v1 - 
  2*j2*l1*q0*r1*v1 - 2*j1*l2*q0*r1*v1 + 2*i2*p1*q0*r1*v1 + 
  2*i1*p2*q0*r1*v1 - 2*j2*l0*q1*r1*v1 - 2*j0*l2*q1*r1*v1 + 
  2*i2*p0*q1*r1*v1 + 2*i0*p2*q1*r1*v1 - 2*j1*l0*q2*r1*v1 - 
  2*j0*l1*q2*r1*v1 + 2*i1*p0*q2*r1*v1 + 2*i0*p1*q2*r1*v1 + 
  2*l1*m1*n0*r2*v1 + 2*l1*m0*n1*r2*v1 + 2*l0*m1*n1*r2*v1 - 
  2*j1*m1*p0*r2*v1 - 2*j1*m0*p1*r2*v1 - 2*j0*m1*p1*r2*v1 - 
  2*j1*l1*q0*r2*v1 + 2*i1*p1*q0*r2*v1 - 2*j1*l0*q1*r2*v1 - 
  2*j0*l1*q1*r2*v1 + 2*i1*p0*q1*r2*v1 + 2*i0*p1*q1*r2*v1 - 
  2*k2*m1*n1*s0*v1 - 2*k1*m2*n1*s0*v1 - 2*k1*m1*n2*s0*v1 + 
  2*j2*m1*o1*s0*v1 + 2*j1*m2*o1*s0*v1 + 2*j1*m1*o2*s0*v1 + 
  2*j2*k1*q1*s0*v1 + 2*j1*k2*q1*s0*v1 - 2*i2*o1*q1*s0*v1 - 
  2*i1*o2*q1*s0*v1 + 2*j1*k1*q2*s0*v1 - 2*i1*o1*q2*s0*v1 - 
  2*k2*m1*n0*s1*v1 - 2*k1*m2*n0*s1*v1 - 2*k2*m0*n1*s1*v1 - 
  2*k0*m2*n1*s1*v1 - 2*k1*m0*n2*s1*v1 - 2*k0*m1*n2*s1*v1 + 
  2*j2*m1*o0*s1*v1 + 2*j1*m2*o0*s1*v1 + 2*j2*m0*o1*s1*v1 + 
  2*j0*m2*o1*s1*v1 + 2*j1*m0*o2*s1*v1 + 2*j0*m1*o2*s1*v1 + 
  2*j2*k1*q0*s1*v1 + 2*j1*k2*q0*s1*v1 - 2*i2*o1*q0*s1*v1 - 
  2*i1*o2*q0*s1*v1 + 2*j2*k0*q1*s1*v1 + 2*j0*k2*q1*s1*v1 - 
  2*i2*o0*q1*s1*v1 - 2*i0*o2*q1*s1*v1 + 2*j1*k0*q2*s1*v1 + 
  2*j0*k1*q2*s1*v1 - 2*i1*o0*q2*s1*v1 - 2*i0*o1*q2*s1*v1 - 
  2*k1*m1*n0*s2*v1 - 2*k1*m0*n1*s2*v1 - 2*k0*m1*n1*s2*v1 + 
  2*j1*m1*o0*s2*v1 + 2*j1*m0*o1*s2*v1 + 2*j0*m1*o1*s2*v1 + 
  2*j1*k1*q0*s2*v1 - 2*i1*o1*q0*s2*v1 + 2*j1*k0*q1*s2*v1 + 
  2*j0*k1*q1*s2*v1 - 2*i1*o0*q1*s2*v1 - 2*i0*o1*q1*s2*v1 - 
  2*k2*l1*n1*t0*v1 - 2*k1*l2*n1*t0*v1 - 2*k1*l1*n2*t0*v1 + 
  2*j2*l1*o1*t0*v1 + 2*j1*l2*o1*t0*v1 + 2*j1*l1*o2*t0*v1 + 
  2*j2*k1*p1*t0*v1 + 2*j1*k2*p1*t0*v1 - 2*i2*o1*p1*t0*v1 - 
  2*i1*o2*p1*t0*v1 + 2*j1*k1*p2*t0*v1 - 2*i1*o1*p2*t0*v1 - 
  4*j1*j2*s1*t0*v1 + 2*i2*n1*s1*t0*v1 + 2*i1*n2*s1*t0*v1 - 
  2*j1*j1*s2*t0*v1 + 2*i1*n1*s2*t0*v1 - 2*k2*l1*n0*t1*v1 - 
  2*k1*l2*n0*t1*v1 - 2*k2*l0*n1*t1*v1 - 2*k0*l2*n1*t1*v1 - 
  2*k1*l0*n2*t1*v1 - 2*k0*l1*n2*t1*v1 + 2*j2*l1*o0*t1*v1 + 
  2*j1*l2*o0*t1*v1 + 2*j2*l0*o1*t1*v1 + 2*j0*l2*o1*t1*v1 + 
  2*j1*l0*o2*t1*v1 + 2*j0*l1*o2*t1*v1 + 2*j2*k1*p0*t1*v1 + 
  2*j1*k2*p0*t1*v1 - 2*i2*o1*p0*t1*v1 - 2*i1*o2*p0*t1*v1 + 
  2*j2*k0*p1*t1*v1 + 2*j0*k2*p1*t1*v1 - 2*i2*o0*p1*t1*v1 - 
  2*i0*o2*p1*t1*v1 + 2*j1*k0*p2*t1*v1 + 2*j0*k1*p2*t1*v1 - 
  2*i1*o0*p2*t1*v1 - 2*i0*o1*p2*t1*v1 - 4*j1*j2*s0*t1*v1 + 
  2*i2*n1*s0*t1*v1 + 2*i1*n2*s0*t1*v1 - 4*j0*j2*s1*t1*v1 + 
  2*i2*n0*s1*t1*v1 + 2*i0*n2*s1*t1*v1 - 4*j0*j1*s2*t1*v1 + 
  2*i1*n0*s2*t1*v1 + 2*i0*n1*s2*t1*v1 - 2*k1*l1*n0*t2*v1 - 
  2*k1*l0*n1*t2*v1 - 2*k0*l1*n1*t2*v1 + 2*j1*l1*o0*t2*v1 + 
  2*j1*l0*o1*t2*v1 + 2*j0*l1*o1*t2*v1 + 2*j1*k1*p0*t2*v1 - 
  2*i1*o1*p0*t2*v1 + 2*j1*k0*p1*t2*v1 + 2*j0*k1*p1*t2*v1 - 
  2*i1*o0*p1*t2*v1 - 2*i0*o1*p1*t2*v1 - 2*j1*j1*s0*t2*v1 + 
  2*i1*n1*s0*t2*v1 - 4*j0*j1*s1*t2*v1 + 2*i1*n0*s1*t2*v1 + 
  2*i0*n1*s1*t2*v1 + 4*k1*k2*n1*v0*v1 + 2*k1*k1*n2*v0*v1 - 
  4*j2*k1*o1*v0*v1 - 4*j1*k2*o1*v0*v1 + 2*i2*o1*o1*v0*v1 - 
  4*j1*k1*o2*v0*v1 + 4*i1*o1*o2*v0*v1 + 4*j1*j2*r1*v0*v1 - 
  2*i2*n1*r1*v0*v1 - 2*i1*n2*r1*v0*v1 + 2*j1*j1*r2*v0*v1 - 
  2*i1*n1*r2*v0*v1 + 2*k1*k2*n0*v1*v1 + 2*k0*k2*n1*v1*v1 + 
  2*k0*k1*n2*v1*v1 - 2*j2*k1*o0*v1*v1 - 2*j1*k2*o0*v1*v1 - 
  2*j2*k0*o1*v1*v1 - 2*j0*k2*o1*v1*v1 + 2*i2*o0*o1*v1*v1 - 
  2*j1*k0*o2*v1*v1 - 2*j0*k1*o2*v1*v1 + 2*i1*o0*o2*v1*v1 + 
  2*i0*o1*o2*v1*v1 + 2*j1*j2*r0*v1*v1 - i2*n1*r0*v1*v1 - 
  i1*n2*r0*v1*v1 + 2*j0*j2*r1*v1*v1 - i2*n0*r1*v1*v1 - 
  i0*n2*r1*v1*v1 + 2*j0*j1*r2*v1*v1 - i1*n0*r2*v1*v1 - 
  i0*n1*r2*v1*v1 - 4*l1*m1*o0*o1*v2 - 2*l1*m0*o1*o1*v2 - 
  2*l0*m1*o1*o1*v2 + 2*k1*m1*o1*p0*v2 + 2*k1*m1*o0*p1*v2 + 
  2*k1*m0*o1*p1*v2 + 2*k0*m1*o1*p1*v2 + 2*k1*l1*o1*q0*v2 - 
  2*k1*k1*p1*q0*v2 + 2*k1*l1*o0*q1*v2 + 2*k1*l0*o1*q1*v2 + 
  2*k0*l1*o1*q1*v2 - 2*k1*k1*p0*q1*v2 - 4*k0*k1*p1*q1*v2 + 
  2*l1*m1*n1*r0*v2 - 2*j1*m1*p1*r0*v2 - 2*j1*l1*q1*r0*v2 + 
  2*i1*p1*q1*r0*v2 + 2*l1*m1*n0*r1*v2 + 2*l1*m0*n1*r1*v2 + 
  2*l0*m1*n1*r1*v2 - 2*j1*m1*p0*r1*v2 - 2*j1*m0*p1*r1*v2 - 
  2*j0*m1*p1*r1*v2 - 2*j1*l1*q0*r1*v2 + 2*i1*p1*q0*r1*v2 - 
  2*j1*l0*q1*r1*v2 - 2*j0*l1*q1*r1*v2 + 2*i1*p0*q1*r1*v2 + 
  2*i0*p1*q1*r1*v2 - 2*k1*m1*n1*s0*v2 + 2*j1*m1*o1*s0*v2 + 
  2*j1*k1*q1*s0*v2 - 2*i1*o1*q1*s0*v2 - 2*k1*m1*n0*s1*v2 - 
  2*k1*m0*n1*s1*v2 - 2*k0*m1*n1*s1*v2 + 2*j1*m1*o0*s1*v2 + 
  2*j1*m0*o1*s1*v2 + 2*j0*m1*o1*s1*v2 + 2*j1*k1*q0*s1*v2 - 
  2*i1*o1*q0*s1*v2 + 2*j1*k0*q1*s1*v2 + 2*j0*k1*q1*s1*v2 - 
  2*i1*o0*q1*s1*v2 - 2*i0*o1*q1*s1*v2 - 2*k1*l1*n1*t0*v2 + 
  2*j1*l1*o1*t0*v2 + 2*j1*k1*p1*t0*v2 - 2*i1*o1*p1*t0*v2 - 
  2*j1*j1*s1*t0*v2 + 2*i1*n1*s1*t0*v2 - 2*k1*l1*n0*t1*v2 - 
  2*k1*l0*n1*t1*v2 - 2*k0*l1*n1*t1*v2 + 2*j1*l1*o0*t1*v2 + 
  2*j1*l0*o1*t1*v2 + 2*j0*l1*o1*t1*v2 + 2*j1*k1*p0*t1*v2 - 
  2*i1*o1*p0*t1*v2 + 2*j1*k0*p1*t1*v2 + 2*j0*k1*p1*t1*v2 - 
  2*i1*o0*p1*t1*v2 - 2*i0*o1*p1*t1*v2 - 2*j1*j1*s0*t1*v2 + 
  2*i1*n1*s0*t1*v2 - 4*j0*j1*s1*t1*v2 + 2*i1*n0*s1*t1*v2 + 
  2*i0*n1*s1*t1*v2 + 2*k1*k1*n1*v0*v2 - 4*j1*k1*o1*v0*v2 + 
  2*i1*o1*o1*v0*v2 + 2*j1*j1*r1*v0*v2 - 2*i1*n1*r1*v0*v2 + 
  2*k1*k1*n0*v1*v2 + 4*k0*k1*n1*v1*v2 - 4*j1*k1*o0*v1*v2 - 
  4*j1*k0*o1*v1*v2 - 4*j0*k1*o1*v1*v2 + 4*i1*o0*o1*v1*v2 + 
  2*i0*o1*o1*v1*v2 + 2*j1*j1*r0*v1*v2 - 2*i1*n1*r0*v1*v2 + 
  4*j0*j1*r1*v1*v2 - 2*i1*n0*r1*v1*v2 - 2*i0*n1*r1*v1*v2 + 
  2*l1*l2*o1*o1*w0 + 2*l1*l1*o1*o2*w0 - 2*k2*l1*o1*p1*w0 - 
  2*k1*l2*o1*p1*w0 - 2*k1*l1*o2*p1*w0 + 2*k1*k2*p1*p1*w0 - 
  2*k1*l1*o1*p2*w0 + 2*k1*k1*p1*p2*w0 - 2*l1*l2*n1*r1*w0 - 
  l1*l1*n2*r1*w0 + 2*j2*l1*p1*r1*w0 + 2*j1*l2*p1*r1*w0 - 
  i2*p1*p1*r1*w0 + 2*j1*l1*p2*r1*w0 - 2*i1*p1*p2*r1*w0 - 
  l1*l1*n1*r2*w0 + 2*j1*l1*p1*r2*w0 - i1*p1*p1*r2*w0 + 
  2*k2*l1*n1*s1*w0 + 2*k1*l2*n1*s1*w0 + 2*k1*l1*n2*s1*w0 - 
  2*j2*l1*o1*s1*w0 - 2*j1*l2*o1*s1*w0 - 2*j1*l1*o2*s1*w0 - 
  2*j2*k1*p1*s1*w0 - 2*j1*k2*p1*s1*w0 + 2*i2*o1*p1*s1*w0 + 
  2*i1*o2*p1*s1*w0 - 2*j1*k1*p2*s1*w0 + 2*i1*o1*p2*s1*w0 + 
  2*j1*j2*s1*s1*w0 - i2*n1*s1*s1*w0 - i1*n2*s1*s1*w0 + 
  2*k1*l1*n1*s2*w0 - 2*j1*l1*o1*s2*w0 - 2*j1*k1*p1*s2*w0 + 
  2*i1*o1*p1*s2*w0 + 2*j1*j1*s1*s2*w0 - 2*i1*n1*s1*s2*w0 - 
  2*k1*k2*n1*u1*w0 - k1*k1*n2*u1*w0 + 2*j2*k1*o1*u1*w0 + 
  2*j1*k2*o1*u1*w0 - i2*o1*o1*u1*w0 + 2*j1*k1*o2*u1*w0 - 
  2*i1*o1*o2*u1*w0 - 2*j1*j2*r1*u1*w0 + i2*n1*r1*u1*w0 + 
  i1*n2*r1*u1*w0 - j1*j1*r2*u1*w0 + i1*n1*r2*u1*w0 - 
  k1*k1*n1*u2*w0 + 2*j1*k1*o1*u2*w0 - i1*o1*o1*u2*w0 - 
  j1*j1*r1*u2*w0 + i1*n1*r1*u2*w0 + 4*l1*l2*o0*o1*w1 + 
  2*l0*l2*o1*o1*w1 + 2*l1*l1*o0*o2*w1 + 4*l0*l1*o1*o2*w1 - 
  2*k2*l1*o1*p0*w1 - 2*k1*l2*o1*p0*w1 - 2*k1*l1*o2*p0*w1 - 
  2*k2*l1*o0*p1*w1 - 2*k1*l2*o0*p1*w1 - 2*k2*l0*o1*p1*w1 - 
  2*k0*l2*o1*p1*w1 - 2*k1*l0*o2*p1*w1 - 2*k0*l1*o2*p1*w1 + 
  4*k1*k2*p0*p1*w1 + 2*k0*k2*p1*p1*w1 - 2*k1*l1*o0*p2*w1 - 
  2*k1*l0*o1*p2*w1 - 2*k0*l1*o1*p2*w1 + 2*k1*k1*p0*p2*w1 + 
  4*k0*k1*p1*p2*w1 - 2*l1*l2*n1*r0*w1 - l1*l1*n2*r0*w1 + 
  2*j2*l1*p1*r0*w1 + 2*j1*l2*p1*r0*w1 - i2*p1*p1*r0*w1 + 
  2*j1*l1*p2*r0*w1 - 2*i1*p1*p2*r0*w1 - 2*l1*l2*n0*r1*w1 - 
  2*l0*l2*n1*r1*w1 - 2*l0*l1*n2*r1*w1 + 2*j2*l1*p0*r1*w1 + 
  2*j1*l2*p0*r1*w1 + 2*j2*l0*p1*r1*w1 + 2*j0*l2*p1*r1*w1 - 
  2*i2*p0*p1*r1*w1 + 2*j1*l0*p2*r1*w1 + 2*j0*l1*p2*r1*w1 - 
  2*i1*p0*p2*r1*w1 - 2*i0*p1*p2*r1*w1 - l1*l1*n0*r2*w1 - 
  2*l0*l1*n1*r2*w1 + 2*j1*l1*p0*r2*w1 + 2*j1*l0*p1*r2*w1 + 
  2*j0*l1*p1*r2*w1 - 2*i1*p0*p1*r2*w1 - i0*p1*p1*r2*w1 + 
  2*k2*l1*n1*s0*w1 + 2*k1*l2*n1*s0*w1 + 2*k1*l1*n2*s0*w1 - 
  2*j2*l1*o1*s0*w1 - 2*j1*l2*o1*s0*w1 - 2*j1*l1*o2*s0*w1 - 
  2*j2*k1*p1*s0*w1 - 2*j1*k2*p1*s0*w1 + 2*i2*o1*p1*s0*w1 + 
  2*i1*o2*p1*s0*w1 - 2*j1*k1*p2*s0*w1 + 2*i1*o1*p2*s0*w1 + 
  2*k2*l1*n0*s1*w1 + 2*k1*l2*n0*s1*w1 + 2*k2*l0*n1*s1*w1 + 
  2*k0*l2*n1*s1*w1 + 2*k1*l0*n2*s1*w1 + 2*k0*l1*n2*s1*w1 - 
  2*j2*l1*o0*s1*w1 - 2*j1*l2*o0*s1*w1 - 2*j2*l0*o1*s1*w1 - 
  2*j0*l2*o1*s1*w1 - 2*j1*l0*o2*s1*w1 - 2*j0*l1*o2*s1*w1 - 
  2*j2*k1*p0*s1*w1 - 2*j1*k2*p0*s1*w1 + 2*i2*o1*p0*s1*w1 + 
  2*i1*o2*p0*s1*w1 - 2*j2*k0*p1*s1*w1 - 2*j0*k2*p1*s1*w1 + 
  2*i2*o0*p1*s1*w1 + 2*i0*o2*p1*s1*w1 - 2*j1*k0*p2*s1*w1 - 
  2*j0*k1*p2*s1*w1 + 2*i1*o0*p2*s1*w1 + 2*i0*o1*p2*s1*w1 + 
  4*j1*j2*s0*s1*w1 - 2*i2*n1*s0*s1*w1 - 2*i1*n2*s0*s1*w1 + 
  2*j0*j2*s1*s1*w1 - i2*n0*s1*s1*w1 - i0*n2*s1*s1*w1 + 
  2*k1*l1*n0*s2*w1 + 2*k1*l0*n1*s2*w1 + 2*k0*l1*n1*s2*w1 - 
  2*j1*l1*o0*s2*w1 - 2*j1*l0*o1*s2*w1 - 2*j0*l1*o1*s2*w1 - 
  2*j1*k1*p0*s2*w1 + 2*i1*o1*p0*s2*w1 - 2*j1*k0*p1*s2*w1 - 
  2*j0*k1*p1*s2*w1 + 2*i1*o0*p1*s2*w1 + 2*i0*o1*p1*s2*w1 + 
  2*j1*j1*s0*s2*w1 - 2*i1*n1*s0*s2*w1 + 4*j0*j1*s1*s2*w1 - 
  2*i1*n0*s1*s2*w1 - 2*i0*n1*s1*s2*w1 - 2*k1*k2*n1*u0*w1 - 
  k1*k1*n2*u0*w1 + 2*j2*k1*o1*u0*w1 + 2*j1*k2*o1*u0*w1 - 
  i2*o1*o1*u0*w1 + 2*j1*k1*o2*u0*w1 - 2*i1*o1*o2*u0*w1 - 
  2*j1*j2*r1*u0*w1 + i2*n1*r1*u0*w1 + i1*n2*r1*u0*w1 - 
  j1*j1*r2*u0*w1 + i1*n1*r2*u0*w1 - 2*k1*k2*n0*u1*w1 - 
  2*k0*k2*n1*u1*w1 - 2*k0*k1*n2*u1*w1 + 2*j2*k1*o0*u1*w1 + 
  2*j1*k2*o0*u1*w1 + 2*j2*k0*o1*u1*w1 + 2*j0*k2*o1*u1*w1 - 
  2*i2*o0*o1*u1*w1 + 2*j1*k0*o2*u1*w1 + 2*j0*k1*o2*u1*w1 - 
  2*i1*o0*o2*u1*w1 - 2*i0*o1*o2*u1*w1 - 2*j1*j2*r0*u1*w1 + 
  i2*n1*r0*u1*w1 + i1*n2*r0*u1*w1 - 2*j0*j2*r1*u1*w1 + 
  i2*n0*r1*u1*w1 + i0*n2*r1*u1*w1 - 2*j0*j1*r2*u1*w1 + 
  i1*n0*r2*u1*w1 + i0*n1*r2*u1*w1 - k1*k1*n0*u2*w1 - 
  2*k0*k1*n1*u2*w1 + 2*j1*k1*o0*u2*w1 + 2*j1*k0*o1*u2*w1 + 
  2*j0*k1*o1*u2*w1 - 2*i1*o0*o1*u2*w1 - i0*o1*o1*u2*w1 - 
  j1*j1*r0*u2*w1 + i1*n1*r0*u2*w1 - 2*j0*j1*r1*u2*w1 + 
  i1*n0*r1*u2*w1 + i0*n1*r1*u2*w1 + 2*l1*l1*o0*o1*w2 + 
  2*l0*l1*o1*o1*w2 - 2*k1*l1*o1*p0*w2 - 2*k1*l1*o0*p1*w2 - 
  2*k1*l0*o1*p1*w2 - 2*k0*l1*o1*p1*w2 + 2*k1*k1*p0*p1*w2 + 
  2*k0*k1*p1*p1*w2 - l1*l1*n1*r0*w2 + 2*j1*l1*p1*r0*w2 - 
  i1*p1*p1*r0*w2 - l1*l1*n0*r1*w2 - 2*l0*l1*n1*r1*w2 + 
  2*j1*l1*p0*r1*w2 + 2*j1*l0*p1*r1*w2 + 2*j0*l1*p1*r1*w2 - 
  2*i1*p0*p1*r1*w2 - i0*p1*p1*r1*w2 + 2*k1*l1*n1*s0*w2 - 
  2*j1*l1*o1*s0*w2 - 2*j1*k1*p1*s0*w2 + 2*i1*o1*p1*s0*w2 + 
  2*k1*l1*n0*s1*w2 + 2*k1*l0*n1*s1*w2 + 2*k0*l1*n1*s1*w2 - 
  2*j1*l1*o0*s1*w2 - 2*j1*l0*o1*s1*w2 - 2*j0*l1*o1*s1*w2 - 
  2*j1*k1*p0*s1*w2 + 2*i1*o1*p0*s1*w2 - 2*j1*k0*p1*s1*w2 - 
  2*j0*k1*p1*s1*w2 + 2*i1*o0*p1*s1*w2 + 2*i0*o1*p1*s1*w2 + 
  2*j1*j1*s0*s1*w2 - 2*i1*n1*s0*s1*w2 + 2*j0*j1*s1*s1*w2 - 
  i1*n0*s1*s1*w2 - i0*n1*s1*s1*w2 - k1*k1*n1*u0*w2 + 
  2*j1*k1*o1*u0*w2 - i1*o1*o1*u0*w2 - j1*j1*r1*u0*w2 + 
  i1*n1*r1*u0*w2 - k1*k1*n0*u1*w2 - 2*k0*k1*n1*u1*w2 + 
  2*j1*k1*o0*u1*w2 + 2*j1*k0*o1*u1*w2 + 2*j0*k1*o1*u1*w2 - 
  2*i1*o0*o1*u1*w2 - i0*o1*o1*u1*w2 - j1*j1*r0*u1*w2 + 
  i1*n1*r0*u1*w2 - 2*j0*j1*r1*u1*w2 + i1*n0*r1*u1*w2 + 
  i0*n1*r1*u1*w2;
    
    k[14] = 2*m1*m2*p1*p1*r1 + 2*m1*m1*p1*p2*r1 - 
  2*l2*m1*p1*q1*r1 - 2*l1*m2*p1*q1*r1 - 2*l1*m1*p2*q1*r1 + 
  2*l1*l2*q1*q1*r1 - 2*l1*m1*p1*q2*r1 + 2*l1*l1*q1*q2*r1 + 
  m1*m1*p1*p1*r2 - 2*l1*m1*p1*q1*r2 + l1*l1*q1*q1*r2 - 
  4*m1*m2*o1*p1*s1 - 2*m1*m1*o2*p1*s1 - 2*m1*m1*o1*p2*s1 + 
  2*l2*m1*o1*q1*s1 + 2*l1*m2*o1*q1*s1 + 2*l1*m1*o2*q1*s1 + 
  2*k2*m1*p1*q1*s1 + 2*k1*m2*p1*q1*s1 + 2*k1*m1*p2*q1*s1 - 
  2*k2*l1*q1*q1*s1 - 2*k1*l2*q1*q1*s1 + 2*l1*m1*o1*q2*s1 + 
  2*k1*m1*p1*q2*s1 - 4*k1*l1*q1*q2*s1 + 2*m1*m2*n1*s1*s1 + 
  m1*m1*n2*s1*s1 - 2*j2*m1*q1*s1*s1 - 2*j1*m2*q1*s1*s1 + 
  i2*q1*q1*s1*s1 - 2*j1*m1*q2*s1*s1 + 2*i1*q1*q2*s1*s1 - 
  2*m1*m1*o1*p1*s2 + 2*l1*m1*o1*q1*s2 + 2*k1*m1*p1*q1*s2 - 
  2*k1*l1*q1*q1*s2 + 2*m1*m1*n1*s1*s2 - 4*j1*m1*q1*s1*s2 + 
  2*i1*q1*q1*s1*s2 + 2*l2*m1*o1*p1*t1 + 2*l1*m2*o1*p1*t1 + 
  2*l1*m1*o2*p1*t1 - 2*k2*m1*p1*p1*t1 - 2*k1*m2*p1*p1*t1 + 
  2*l1*m1*o1*p2*t1 - 4*k1*m1*p1*p2*t1 - 4*l1*l2*o1*q1*t1 - 
  2*l1*l1*o2*q1*t1 + 2*k2*l1*p1*q1*t1 + 2*k1*l2*p1*q1*t1 + 
  2*k1*l1*p2*q1*t1 - 2*l1*l1*o1*q2*t1 + 2*k1*l1*p1*q2*t1 - 
  2*l2*m1*n1*s1*t1 - 2*l1*m2*n1*s1*t1 - 2*l1*m1*n2*s1*t1 + 
  2*j2*m1*p1*s1*t1 + 2*j1*m2*p1*s1*t1 + 2*j1*m1*p2*s1*t1 + 
  2*j2*l1*q1*s1*t1 + 2*j1*l2*q1*s1*t1 - 2*i2*p1*q1*s1*t1 - 
  2*i1*p2*q1*s1*t1 + 2*j1*l1*q2*s1*t1 - 2*i1*p1*q2*s1*t1 - 
  2*l1*m1*n1*s2*t1 + 2*j1*m1*p1*s2*t1 + 2*j1*l1*q1*s2*t1 - 
  2*i1*p1*q1*s2*t1 + 2*l1*l2*n1*t1*t1 + l1*l1*n2*t1*t1 - 
  2*j2*l1*p1*t1*t1 - 2*j1*l2*p1*t1*t1 + i2*p1*p1*t1*t1 - 
  2*j1*l1*p2*t1*t1 + 2*i1*p1*p2*t1*t1 + 2*l1*m1*o1*p1*t2 - 
  2*k1*m1*p1*p1*t2 - 2*l1*l1*o1*q1*t2 + 2*k1*l1*p1*q1*t2 - 
  2*l1*m1*n1*s1*t2 + 2*j1*m1*p1*s1*t2 + 2*j1*l1*q1*s1*t2 - 
  2*i1*p1*q1*s1*t2 + 2*l1*l1*n1*t1*t2 - 4*j1*l1*p1*t1*t2 + 
  2*i1*p1*p1*t1*t2 + 2*m1*m2*o1*o1*u1 + 2*m1*m1*o1*o2*u1 - 
  2*k2*m1*o1*q1*u1 - 2*k1*m2*o1*q1*u1 - 2*k1*m1*o2*q1*u1 + 
  2*k1*k2*q1*q1*u1 - 2*k1*m1*o1*q2*u1 + 2*k1*k1*q1*q2*u1 - 
  2*m1*m2*n1*r1*u1 - m1*m1*n2*r1*u1 + 2*j2*m1*q1*r1*u1 + 
  2*j1*m2*q1*r1*u1 - i2*q1*q1*r1*u1 + 2*j1*m1*q2*r1*u1 - 
  2*i1*q1*q2*r1*u1 - m1*m1*n1*r2*u1 + 2*j1*m1*q1*r2*u1 - 
  i1*q1*q1*r2*u1 + 2*k2*m1*n1*t1*u1 + 2*k1*m2*n1*t1*u1 + 
  2*k1*m1*n2*t1*u1 - 2*j2*m1*o1*t1*u1 - 2*j1*m2*o1*t1*u1 - 
  2*j1*m1*o2*t1*u1 - 2*j2*k1*q1*t1*u1 - 2*j1*k2*q1*t1*u1 + 
  2*i2*o1*q1*t1*u1 + 2*i1*o2*q1*t1*u1 - 2*j1*k1*q2*t1*u1 + 
  2*i1*o1*q2*t1*u1 + 2*j1*j2*t1*t1*u1 - i2*n1*t1*t1*u1 - 
  i1*n2*t1*t1*u1 + 2*k1*m1*n1*t2*u1 - 2*j1*m1*o1*t2*u1 - 
  2*j1*k1*q1*t2*u1 + 2*i1*o1*q1*t2*u1 + 2*j1*j1*t1*t2*u1 - 
  2*i1*n1*t1*t2*u1 + m1*m1*o1*o1*u2 - 2*k1*m1*o1*q1*u2 + 
  k1*k1*q1*q1*u2 - m1*m1*n1*r1*u2 + 2*j1*m1*q1*r1*u2 - 
  i1*q1*q1*r1*u2 + 2*k1*m1*n1*t1*u2 - 2*j1*m1*o1*t1*u2 - 
  2*j1*k1*q1*t1*u2 + 2*i1*o1*q1*t1*u2 + j1*j1*t1*t1*u2 - 
  i1*n1*t1*t1*u2 - 2*l2*m1*o1*o1*v1 - 2*l1*m2*o1*o1*v1 - 
  4*l1*m1*o1*o2*v1 + 2*k2*m1*o1*p1*v1 + 2*k1*m2*o1*p1*v1 + 
  2*k1*m1*o2*p1*v1 + 2*k1*m1*o1*p2*v1 + 2*k2*l1*o1*q1*v1 + 
  2*k1*l2*o1*q1*v1 + 2*k1*l1*o2*q1*v1 - 4*k1*k2*p1*q1*v1 - 
  2*k1*k1*p2*q1*v1 + 2*k1*l1*o1*q2*v1 - 2*k1*k1*p1*q2*v1 + 
  2*l2*m1*n1*r1*v1 + 2*l1*m2*n1*r1*v1 + 2*l1*m1*n2*r1*v1 - 
  2*j2*m1*p1*r1*v1 - 2*j1*m2*p1*r1*v1 - 2*j1*m1*p2*r1*v1 - 
  2*j2*l1*q1*r1*v1 - 2*j1*l2*q1*r1*v1 + 2*i2*p1*q1*r1*v1 + 
  2*i1*p2*q1*r1*v1 - 2*j1*l1*q2*r1*v1 + 2*i1*p1*q2*r1*v1 + 
  2*l1*m1*n1*r2*v1 - 2*j1*m1*p1*r2*v1 - 2*j1*l1*q1*r2*v1 + 
  2*i1*p1*q1*r2*v1 - 2*k2*m1*n1*s1*v1 - 2*k1*m2*n1*s1*v1 - 
  2*k1*m1*n2*s1*v1 + 2*j2*m1*o1*s1*v1 + 2*j1*m2*o1*s1*v1 + 
  2*j1*m1*o2*s1*v1 + 2*j2*k1*q1*s1*v1 + 2*j1*k2*q1*s1*v1 - 
  2*i2*o1*q1*s1*v1 - 2*i1*o2*q1*s1*v1 + 2*j1*k1*q2*s1*v1 - 
  2*i1*o1*q2*s1*v1 - 2*k1*m1*n1*s2*v1 + 2*j1*m1*o1*s2*v1 + 
  2*j1*k1*q1*s2*v1 - 2*i1*o1*q1*s2*v1 - 2*k2*l1*n1*t1*v1 - 
  2*k1*l2*n1*t1*v1 - 2*k1*l1*n2*t1*v1 + 2*j2*l1*o1*t1*v1 + 
  2*j1*l2*o1*t1*v1 + 2*j1*l1*o2*t1*v1 + 2*j2*k1*p1*t1*v1 + 
  2*j1*k2*p1*t1*v1 - 2*i2*o1*p1*t1*v1 - 2*i1*o2*p1*t1*v1 + 
  2*j1*k1*p2*t1*v1 - 2*i1*o1*p2*t1*v1 - 4*j1*j2*s1*t1*v1 + 
  2*i2*n1*s1*t1*v1 + 2*i1*n2*s1*t1*v1 - 2*j1*j1*s2*t1*v1 + 
  2*i1*n1*s2*t1*v1 - 2*k1*l1*n1*t2*v1 + 2*j1*l1*o1*t2*v1 + 
  2*j1*k1*p1*t2*v1 - 2*i1*o1*p1*t2*v1 - 2*j1*j1*s1*t2*v1 + 
  2*i1*n1*s1*t2*v1 + 2*k1*k2*n1*v1*v1 + k1*k1*n2*v1*v1 - 
  2*j2*k1*o1*v1*v1 - 2*j1*k2*o1*v1*v1 + i2*o1*o1*v1*v1 - 
  2*j1*k1*o2*v1*v1 + 2*i1*o1*o2*v1*v1 + 2*j1*j2*r1*v1*v1 - 
  i2*n1*r1*v1*v1 - i1*n2*r1*v1*v1 + j1*j1*r2*v1*v1 - 
  i1*n1*r2*v1*v1 - 2*l1*m1*o1*o1*v2 + 2*k1*m1*o1*p1*v2 + 
  2*k1*l1*o1*q1*v2 - 2*k1*k1*p1*q1*v2 + 2*l1*m1*n1*r1*v2 - 
  2*j1*m1*p1*r1*v2 - 2*j1*l1*q1*r1*v2 + 2*i1*p1*q1*r1*v2 - 
  2*k1*m1*n1*s1*v2 + 2*j1*m1*o1*s1*v2 + 2*j1*k1*q1*s1*v2 - 
  2*i1*o1*q1*s1*v2 - 2*k1*l1*n1*t1*v2 + 2*j1*l1*o1*t1*v2 + 
  2*j1*k1*p1*t1*v2 - 2*i1*o1*p1*t1*v2 - 2*j1*j1*s1*t1*v2 + 
  2*i1*n1*s1*t1*v2 + 2*k1*k1*n1*v1*v2 - 4*j1*k1*o1*v1*v2 + 
  2*i1*o1*o1*v1*v2 + 2*j1*j1*r1*v1*v2 - 2*i1*n1*r1*v1*v2 + 
  2*l1*l2*o1*o1*w1 + 2*l1*l1*o1*o2*w1 - 2*k2*l1*o1*p1*w1 - 
  2*k1*l2*o1*p1*w1 - 2*k1*l1*o2*p1*w1 + 2*k1*k2*p1*p1*w1 - 
  2*k1*l1*o1*p2*w1 + 2*k1*k1*p1*p2*w1 - 2*l1*l2*n1*r1*w1 - 
  l1*l1*n2*r1*w1 + 2*j2*l1*p1*r1*w1 + 2*j1*l2*p1*r1*w1 - 
  i2*p1*p1*r1*w1 + 2*j1*l1*p2*r1*w1 - 2*i1*p1*p2*r1*w1 - 
  l1*l1*n1*r2*w1 + 2*j1*l1*p1*r2*w1 - i1*p1*p1*r2*w1 + 
  2*k2*l1*n1*s1*w1 + 2*k1*l2*n1*s1*w1 + 2*k1*l1*n2*s1*w1 - 
  2*j2*l1*o1*s1*w1 - 2*j1*l2*o1*s1*w1 - 2*j1*l1*o2*s1*w1 - 
  2*j2*k1*p1*s1*w1 - 2*j1*k2*p1*s1*w1 + 2*i2*o1*p1*s1*w1 + 
  2*i1*o2*p1*s1*w1 - 2*j1*k1*p2*s1*w1 + 2*i1*o1*p2*s1*w1 + 
  2*j1*j2*s1*s1*w1 - i2*n1*s1*s1*w1 - i1*n2*s1*s1*w1 + 
  2*k1*l1*n1*s2*w1 - 2*j1*l1*o1*s2*w1 - 2*j1*k1*p1*s2*w1 + 
  2*i1*o1*p1*s2*w1 + 2*j1*j1*s1*s2*w1 - 2*i1*n1*s1*s2*w1 - 
  2*k1*k2*n1*u1*w1 - k1*k1*n2*u1*w1 + 2*j2*k1*o1*u1*w1 + 
  2*j1*k2*o1*u1*w1 - i2*o1*o1*u1*w1 + 2*j1*k1*o2*u1*w1 - 
  2*i1*o1*o2*u1*w1 - 2*j1*j2*r1*u1*w1 + i2*n1*r1*u1*w1 + 
  i1*n2*r1*u1*w1 - j1*j1*r2*u1*w1 + i1*n1*r2*u1*w1 - 
  k1*k1*n1*u2*w1 + 2*j1*k1*o1*u2*w1 - i1*o1*o1*u2*w1 - 
  j1*j1*r1*u2*w1 + i1*n1*r1*u2*w1 + l1*l1*o1*o1*w2 - 
  2*k1*l1*o1*p1*w2 + k1*k1*p1*p1*w2 - l1*l1*n1*r1*w2 + 
  2*j1*l1*p1*r1*w2 - i1*p1*p1*r1*w2 + 2*k1*l1*n1*s1*w2 - 
  2*j1*l1*o1*s1*w2 - 2*j1*k1*p1*s1*w2 + 2*i1*o1*p1*s1*w2 + 
  j1*j1*s1*s1*w2 - i1*n1*s1*s1*w2 - k1*k1*n1*u1*w2 + 
  2*j1*k1*o1*u1*w2 - i1*o1*o1*u1*w2 - j1*j1*r1*u1*w2 + 
  i1*n1*r1*u1*w2;
   
    k[15] = m0*m0*p0*p0*r0 - 2*l0*m0*p0*q0*r0 + 
  l0*l0*q0*q0*r0 - 2*m0*m0*o0*p0*s0 + 2*l0*m0*o0*q0*s0 + 
  2*k0*m0*p0*q0*s0 - 2*k0*l0*q0*q0*s0 + m0*m0*n0*s0*s0 - 
  2*j0*m0*q0*s0*s0 + i0*q0*q0*s0*s0 + 2*l0*m0*o0*p0*t0 - 
  2*k0*m0*p0*p0*t0 - 2*l0*l0*o0*q0*t0 + 2*k0*l0*p0*q0*t0 - 
  2*l0*m0*n0*s0*t0 + 2*j0*m0*p0*s0*t0 + 2*j0*l0*q0*s0*t0 - 
  2*i0*p0*q0*s0*t0 + l0*l0*n0*t0*t0 - 2*j0*l0*p0*t0*t0 + 
  i0*p0*p0*t0*t0 + m0*m0*o0*o0*u0 - 2*k0*m0*o0*q0*u0 + 
  k0*k0*q0*q0*u0 - m0*m0*n0*r0*u0 + 2*j0*m0*q0*r0*u0 - 
  i0*q0*q0*r0*u0 + 2*k0*m0*n0*t0*u0 - 2*j0*m0*o0*t0*u0 - 
  2*j0*k0*q0*t0*u0 + 2*i0*o0*q0*t0*u0 + j0*j0*t0*t0*u0 - 
  i0*n0*t0*t0*u0 - 2*l0*m0*o0*o0*v0 + 2*k0*m0*o0*p0*v0 + 
  2*k0*l0*o0*q0*v0 - 2*k0*k0*p0*q0*v0 + 2*l0*m0*n0*r0*v0 - 
  2*j0*m0*p0*r0*v0 - 2*j0*l0*q0*r0*v0 + 2*i0*p0*q0*r0*v0 - 
  2*k0*m0*n0*s0*v0 + 2*j0*m0*o0*s0*v0 + 2*j0*k0*q0*s0*v0 - 
  2*i0*o0*q0*s0*v0 - 2*k0*l0*n0*t0*v0 + 2*j0*l0*o0*t0*v0 + 
  2*j0*k0*p0*t0*v0 - 2*i0*o0*p0*t0*v0 - 2*j0*j0*s0*t0*v0 + 
  2*i0*n0*s0*t0*v0 + k0*k0*n0*v0*v0 - 2*j0*k0*o0*v0*v0 + 
  i0*o0*o0*v0*v0 + j0*j0*r0*v0*v0 - i0*n0*r0*v0*v0 + 
  l0*l0*o0*o0*w0 - 2*k0*l0*o0*p0*w0 + k0*k0*p0*p0*w0 - 
  l0*l0*n0*r0*w0 + 2*j0*l0*p0*r0*w0 - i0*p0*p0*r0*w0 + 
  2*k0*l0*n0*s0*w0 - 2*j0*l0*o0*s0*w0 - 2*j0*k0*p0*s0*w0 + 
  2*i0*o0*p0*s0*w0 + j0*j0*s0*s0*w0 - i0*n0*s0*s0*w0 - 
  k0*k0*n0*u0*w0 + 2*j0*k0*o0*u0*w0 - i0*o0*o0*u0*w0 - 
  j0*j0*r0*u0*w0 + i0*n0*r0*u0*w0;
   
    k[16] = 2*m0*m1*p0*p0*r0 + 
  2*m0*m0*p0*p1*r0 - 2*l1*m0*p0*q0*r0 - 2*l0*m1*p0*q0*r0 - 
  2*l0*m0*p1*q0*r0 + 2*l0*l1*q0*q0*r0 - 2*l0*m0*p0*q1*r0 + 
  2*l0*l0*q0*q1*r0 + m0*m0*p0*p0*r1 - 2*l0*m0*p0*q0*r1 + 
  l0*l0*q0*q0*r1 - 4*m0*m1*o0*p0*s0 - 2*m0*m0*o1*p0*s0 - 
  2*m0*m0*o0*p1*s0 + 2*l1*m0*o0*q0*s0 + 2*l0*m1*o0*q0*s0 + 
  2*l0*m0*o1*q0*s0 + 2*k1*m0*p0*q0*s0 + 2*k0*m1*p0*q0*s0 + 
  2*k0*m0*p1*q0*s0 - 2*k1*l0*q0*q0*s0 - 2*k0*l1*q0*q0*s0 + 
  2*l0*m0*o0*q1*s0 + 2*k0*m0*p0*q1*s0 - 4*k0*l0*q0*q1*s0 + 
  2*m0*m1*n0*s0*s0 + m0*m0*n1*s0*s0 - 2*j1*m0*q0*s0*s0 - 
  2*j0*m1*q0*s0*s0 + i1*q0*q0*s0*s0 - 2*j0*m0*q1*s0*s0 + 
  2*i0*q0*q1*s0*s0 - 2*m0*m0*o0*p0*s1 + 2*l0*m0*o0*q0*s1 + 
  2*k0*m0*p0*q0*s1 - 2*k0*l0*q0*q0*s1 + 2*m0*m0*n0*s0*s1 - 
  4*j0*m0*q0*s0*s1 + 2*i0*q0*q0*s0*s1 + 2*l1*m0*o0*p0*t0 + 
  2*l0*m1*o0*p0*t0 + 2*l0*m0*o1*p0*t0 - 2*k1*m0*p0*p0*t0 - 
  2*k0*m1*p0*p0*t0 + 2*l0*m0*o0*p1*t0 - 4*k0*m0*p0*p1*t0 - 
  4*l0*l1*o0*q0*t0 - 2*l0*l0*o1*q0*t0 + 2*k1*l0*p0*q0*t0 + 
  2*k0*l1*p0*q0*t0 + 2*k0*l0*p1*q0*t0 - 2*l0*l0*o0*q1*t0 + 
  2*k0*l0*p0*q1*t0 - 2*l1*m0*n0*s0*t0 - 2*l0*m1*n0*s0*t0 - 
  2*l0*m0*n1*s0*t0 + 2*j1*m0*p0*s0*t0 + 2*j0*m1*p0*s0*t0 + 
  2*j0*m0*p1*s0*t0 + 2*j1*l0*q0*s0*t0 + 2*j0*l1*q0*s0*t0 - 
  2*i1*p0*q0*s0*t0 - 2*i0*p1*q0*s0*t0 + 2*j0*l0*q1*s0*t0 - 
  2*i0*p0*q1*s0*t0 - 2*l0*m0*n0*s1*t0 + 2*j0*m0*p0*s1*t0 + 
  2*j0*l0*q0*s1*t0 - 2*i0*p0*q0*s1*t0 + 2*l0*l1*n0*t0*t0 + 
  l0*l0*n1*t0*t0 - 2*j1*l0*p0*t0*t0 - 2*j0*l1*p0*t0*t0 + 
  i1*p0*p0*t0*t0 - 2*j0*l0*p1*t0*t0 + 2*i0*p0*p1*t0*t0 + 
  2*l0*m0*o0*p0*t1 - 2*k0*m0*p0*p0*t1 - 2*l0*l0*o0*q0*t1 + 
  2*k0*l0*p0*q0*t1 - 2*l0*m0*n0*s0*t1 + 2*j0*m0*p0*s0*t1 + 
  2*j0*l0*q0*s0*t1 - 2*i0*p0*q0*s0*t1 + 2*l0*l0*n0*t0*t1 - 
  4*j0*l0*p0*t0*t1 + 2*i0*p0*p0*t0*t1 + 2*m0*m1*o0*o0*u0 + 
  2*m0*m0*o0*o1*u0 - 2*k1*m0*o0*q0*u0 - 2*k0*m1*o0*q0*u0 - 
  2*k0*m0*o1*q0*u0 + 2*k0*k1*q0*q0*u0 - 2*k0*m0*o0*q1*u0 + 
  2*k0*k0*q0*q1*u0 - 2*m0*m1*n0*r0*u0 - m0*m0*n1*r0*u0 + 
  2*j1*m0*q0*r0*u0 + 2*j0*m1*q0*r0*u0 - i1*q0*q0*r0*u0 + 
  2*j0*m0*q1*r0*u0 - 2*i0*q0*q1*r0*u0 - m0*m0*n0*r1*u0 + 
  2*j0*m0*q0*r1*u0 - i0*q0*q0*r1*u0 + 2*k1*m0*n0*t0*u0 + 
  2*k0*m1*n0*t0*u0 + 2*k0*m0*n1*t0*u0 - 2*j1*m0*o0*t0*u0 - 
  2*j0*m1*o0*t0*u0 - 2*j0*m0*o1*t0*u0 - 2*j1*k0*q0*t0*u0 - 
  2*j0*k1*q0*t0*u0 + 2*i1*o0*q0*t0*u0 + 2*i0*o1*q0*t0*u0 - 
  2*j0*k0*q1*t0*u0 + 2*i0*o0*q1*t0*u0 + 2*j0*j1*t0*t0*u0 - 
  i1*n0*t0*t0*u0 - i0*n1*t0*t0*u0 + 2*k0*m0*n0*t1*u0 - 
  2*j0*m0*o0*t1*u0 - 2*j0*k0*q0*t1*u0 + 2*i0*o0*q0*t1*u0 + 
  2*j0*j0*t0*t1*u0 - 2*i0*n0*t0*t1*u0 + m0*m0*o0*o0*u1 - 
  2*k0*m0*o0*q0*u1 + k0*k0*q0*q0*u1 - m0*m0*n0*r0*u1 + 
  2*j0*m0*q0*r0*u1 - i0*q0*q0*r0*u1 + 2*k0*m0*n0*t0*u1 - 
  2*j0*m0*o0*t0*u1 - 2*j0*k0*q0*t0*u1 + 2*i0*o0*q0*t0*u1 + 
  j0*j0*t0*t0*u1 - i0*n0*t0*t0*u1 - 2*l1*m0*o0*o0*v0 - 
  2*l0*m1*o0*o0*v0 - 4*l0*m0*o0*o1*v0 + 2*k1*m0*o0*p0*v0 + 
  2*k0*m1*o0*p0*v0 + 2*k0*m0*o1*p0*v0 + 2*k0*m0*o0*p1*v0 + 
  2*k1*l0*o0*q0*v0 + 2*k0*l1*o0*q0*v0 + 2*k0*l0*o1*q0*v0 - 
  4*k0*k1*p0*q0*v0 - 2*k0*k0*p1*q0*v0 + 2*k0*l0*o0*q1*v0 - 
  2*k0*k0*p0*q1*v0 + 2*l1*m0*n0*r0*v0 + 2*l0*m1*n0*r0*v0 + 
  2*l0*m0*n1*r0*v0 - 2*j1*m0*p0*r0*v0 - 2*j0*m1*p0*r0*v0 - 
  2*j0*m0*p1*r0*v0 - 2*j1*l0*q0*r0*v0 - 2*j0*l1*q0*r0*v0 + 
  2*i1*p0*q0*r0*v0 + 2*i0*p1*q0*r0*v0 - 2*j0*l0*q1*r0*v0 + 
  2*i0*p0*q1*r0*v0 + 2*l0*m0*n0*r1*v0 - 2*j0*m0*p0*r1*v0 - 
  2*j0*l0*q0*r1*v0 + 2*i0*p0*q0*r1*v0 - 2*k1*m0*n0*s0*v0 - 
  2*k0*m1*n0*s0*v0 - 2*k0*m0*n1*s0*v0 + 2*j1*m0*o0*s0*v0 + 
  2*j0*m1*o0*s0*v0 + 2*j0*m0*o1*s0*v0 + 2*j1*k0*q0*s0*v0 + 
  2*j0*k1*q0*s0*v0 - 2*i1*o0*q0*s0*v0 - 2*i0*o1*q0*s0*v0 + 
  2*j0*k0*q1*s0*v0 - 2*i0*o0*q1*s0*v0 - 2*k0*m0*n0*s1*v0 + 
  2*j0*m0*o0*s1*v0 + 2*j0*k0*q0*s1*v0 - 2*i0*o0*q0*s1*v0 - 
  2*k1*l0*n0*t0*v0 - 2*k0*l1*n0*t0*v0 - 2*k0*l0*n1*t0*v0 + 
  2*j1*l0*o0*t0*v0 + 2*j0*l1*o0*t0*v0 + 2*j0*l0*o1*t0*v0 + 
  2*j1*k0*p0*t0*v0 + 2*j0*k1*p0*t0*v0 - 2*i1*o0*p0*t0*v0 - 
  2*i0*o1*p0*t0*v0 + 2*j0*k0*p1*t0*v0 - 2*i0*o0*p1*t0*v0 - 
  4*j0*j1*s0*t0*v0 + 2*i1*n0*s0*t0*v0 + 2*i0*n1*s0*t0*v0 - 
  2*j0*j0*s1*t0*v0 + 2*i0*n0*s1*t0*v0 - 2*k0*l0*n0*t1*v0 + 
  2*j0*l0*o0*t1*v0 + 2*j0*k0*p0*t1*v0 - 2*i0*o0*p0*t1*v0 - 
  2*j0*j0*s0*t1*v0 + 2*i0*n0*s0*t1*v0 + 2*k0*k1*n0*v0*v0 + 
  k0*k0*n1*v0*v0 - 2*j1*k0*o0*v0*v0 - 2*j0*k1*o0*v0*v0 + 
  i1*o0*o0*v0*v0 - 2*j0*k0*o1*v0*v0 + 2*i0*o0*o1*v0*v0 + 
  2*j0*j1*r0*v0*v0 - i1*n0*r0*v0*v0 - i0*n1*r0*v0*v0 + 
  j0*j0*r1*v0*v0 - i0*n0*r1*v0*v0 - 2*l0*m0*o0*o0*v1 + 
  2*k0*m0*o0*p0*v1 + 2*k0*l0*o0*q0*v1 - 2*k0*k0*p0*q0*v1 + 
  2*l0*m0*n0*r0*v1 - 2*j0*m0*p0*r0*v1 - 2*j0*l0*q0*r0*v1 + 
  2*i0*p0*q0*r0*v1 - 2*k0*m0*n0*s0*v1 + 2*j0*m0*o0*s0*v1 + 
  2*j0*k0*q0*s0*v1 - 2*i0*o0*q0*s0*v1 - 2*k0*l0*n0*t0*v1 + 
  2*j0*l0*o0*t0*v1 + 2*j0*k0*p0*t0*v1 - 2*i0*o0*p0*t0*v1 - 
  2*j0*j0*s0*t0*v1 + 2*i0*n0*s0*t0*v1 + 2*k0*k0*n0*v0*v1 - 
  4*j0*k0*o0*v0*v1 + 2*i0*o0*o0*v0*v1 + 2*j0*j0*r0*v0*v1 - 
  2*i0*n0*r0*v0*v1 + 2*l0*l1*o0*o0*w0 + 2*l0*l0*o0*o1*w0 - 
  2*k1*l0*o0*p0*w0 - 2*k0*l1*o0*p0*w0 - 2*k0*l0*o1*p0*w0 + 
  2*k0*k1*p0*p0*w0 - 2*k0*l0*o0*p1*w0 + 2*k0*k0*p0*p1*w0 - 
  2*l0*l1*n0*r0*w0 - l0*l0*n1*r0*w0 + 2*j1*l0*p0*r0*w0 + 
  2*j0*l1*p0*r0*w0 - i1*p0*p0*r0*w0 + 2*j0*l0*p1*r0*w0 - 
  2*i0*p0*p1*r0*w0 - l0*l0*n0*r1*w0 + 2*j0*l0*p0*r1*w0 - 
  i0*p0*p0*r1*w0 + 2*k1*l0*n0*s0*w0 + 2*k0*l1*n0*s0*w0 + 
  2*k0*l0*n1*s0*w0 - 2*j1*l0*o0*s0*w0 - 2*j0*l1*o0*s0*w0 - 
  2*j0*l0*o1*s0*w0 - 2*j1*k0*p0*s0*w0 - 2*j0*k1*p0*s0*w0 + 
  2*i1*o0*p0*s0*w0 + 2*i0*o1*p0*s0*w0 - 2*j0*k0*p1*s0*w0 + 
  2*i0*o0*p1*s0*w0 + 2*j0*j1*s0*s0*w0 - i1*n0*s0*s0*w0 - 
  i0*n1*s0*s0*w0 + 2*k0*l0*n0*s1*w0 - 2*j0*l0*o0*s1*w0 - 
  2*j0*k0*p0*s1*w0 + 2*i0*o0*p0*s1*w0 + 2*j0*j0*s0*s1*w0 - 
  2*i0*n0*s0*s1*w0 - 2*k0*k1*n0*u0*w0 - k0*k0*n1*u0*w0 + 
  2*j1*k0*o0*u0*w0 + 2*j0*k1*o0*u0*w0 - i1*o0*o0*u0*w0 + 
  2*j0*k0*o1*u0*w0 - 2*i0*o0*o1*u0*w0 - 2*j0*j1*r0*u0*w0 + 
  i1*n0*r0*u0*w0 + i0*n1*r0*u0*w0 - j0*j0*r1*u0*w0 + 
  i0*n0*r1*u0*w0 - k0*k0*n0*u1*w0 + 2*j0*k0*o0*u1*w0 - 
  i0*o0*o0*u1*w0 - j0*j0*r0*u1*w0 + i0*n0*r0*u1*w0 + 
  l0*l0*o0*o0*w1 - 2*k0*l0*o0*p0*w1 + k0*k0*p0*p0*w1 - 
  l0*l0*n0*r0*w1 + 2*j0*l0*p0*r0*w1 - i0*p0*p0*r0*w1 + 
  2*k0*l0*n0*s0*w1 - 2*j0*l0*o0*s0*w1 - 2*j0*k0*p0*s0*w1 + 
  2*i0*o0*p0*s0*w1 + j0*j0*s0*s0*w1 - i0*n0*s0*s0*w1 - 
  k0*k0*n0*u0*w1 + 2*j0*k0*o0*u0*w1 - i0*o0*o0*u0*w1 - 
  j0*j0*r0*u0*w1 + i0*n0*r0*u0*w1;
   
    k[17] = m1*m1*p0*p0*r0 + 4*m0*m1*p0*p1*r0 + 
  m0*m0*p1*p1*r0 - 2*l1*m1*p0*q0*r0 - 
  2*l1*m0*p1*q0*r0 - 2*l0*m1*p1*q0*r0 + 
  l1*l1*q0*q0*r0 - 2*l1*m0*p0*q1*r0 - 
  2*l0*m1*p0*q1*r0 - 2*l0*m0*p1*q1*r0 + 
  4*l0*l1*q0*q1*r0 + l0*l0*q1*q1*r0 + 2*m0*m1*p0*p0*r1 + 
  2*m0*m0*p0*p1*r1 - 2*l1*m0*p0*q0*r1 - 
  2*l0*m1*p0*q0*r1 - 2*l0*m0*p1*q0*r1 + 
  2*l0*l1*q0*q0*r1 - 2*l0*m0*p0*q1*r1 + 
  2*l0*l0*q0*q1*r1 - 2*m1*m1*o0*p0*s0 - 
  4*m0*m1*o1*p0*s0 - 4*m0*m1*o0*p1*s0 - 
  2*m0*m0*o1*p1*s0 + 2*l1*m1*o0*q0*s0 + 
  2*l1*m0*o1*q0*s0 + 2*l0*m1*o1*q0*s0 + 
  2*k1*m1*p0*q0*s0 + 2*k1*m0*p1*q0*s0 + 
  2*k0*m1*p1*q0*s0 - 2*k1*l1*q0*q0*s0 + 
  2*l1*m0*o0*q1*s0 + 2*l0*m1*o0*q1*s0 + 
  2*l0*m0*o1*q1*s0 + 2*k1*m0*p0*q1*s0 + 
  2*k0*m1*p0*q1*s0 + 2*k0*m0*p1*q1*s0 - 
  4*k1*l0*q0*q1*s0 - 4*k0*l1*q0*q1*s0 - 
  2*k0*l0*q1*q1*s0 + m1*m1*n0*s0*s0 + 2*m0*m1*n1*s0*s0 - 
  2*j1*m1*q0*s0*s0 - 2*j1*m0*q1*s0*s0 - 
  2*j0*m1*q1*s0*s0 + 2*i1*q0*q1*s0*s0 + i0*q1*q1*s0*s0 - 
  4*m0*m1*o0*p0*s1 - 2*m0*m0*o1*p0*s1 - 
  2*m0*m0*o0*p1*s1 + 2*l1*m0*o0*q0*s1 + 
  2*l0*m1*o0*q0*s1 + 2*l0*m0*o1*q0*s1 + 
  2*k1*m0*p0*q0*s1 + 2*k0*m1*p0*q0*s1 + 
  2*k0*m0*p1*q0*s1 - 2*k1*l0*q0*q0*s1 - 
  2*k0*l1*q0*q0*s1 + 2*l0*m0*o0*q1*s1 + 
  2*k0*m0*p0*q1*s1 - 4*k0*l0*q0*q1*s1 + 
  4*m0*m1*n0*s0*s1 + 2*m0*m0*n1*s0*s1 - 
  4*j1*m0*q0*s0*s1 - 4*j0*m1*q0*s0*s1 + 
  2*i1*q0*q0*s0*s1 - 4*j0*m0*q1*s0*s1 + 
  4*i0*q0*q1*s0*s1 + m0*m0*n0*s1*s1 - 2*j0*m0*q0*s1*s1 + 
  i0*q0*q0*s1*s1 + 2*l1*m1*o0*p0*t0 + 
  2*l1*m0*o1*p0*t0 + 2*l0*m1*o1*p0*t0 - 
  2*k1*m1*p0*p0*t0 + 2*l1*m0*o0*p1*t0 + 
  2*l0*m1*o0*p1*t0 + 2*l0*m0*o1*p1*t0 - 
  4*k1*m0*p0*p1*t0 - 4*k0*m1*p0*p1*t0 - 
  2*k0*m0*p1*p1*t0 - 2*l1*l1*o0*q0*t0 - 
  4*l0*l1*o1*q0*t0 + 2*k1*l1*p0*q0*t0 + 
  2*k1*l0*p1*q0*t0 + 2*k0*l1*p1*q0*t0 - 
  4*l0*l1*o0*q1*t0 - 2*l0*l0*o1*q1*t0 + 
  2*k1*l0*p0*q1*t0 + 2*k0*l1*p0*q1*t0 + 
  2*k0*l0*p1*q1*t0 - 2*l1*m1*n0*s0*t0 - 
  2*l1*m0*n1*s0*t0 - 2*l0*m1*n1*s0*t0 + 
  2*j1*m1*p0*s0*t0 + 2*j1*m0*p1*s0*t0 + 
  2*j0*m1*p1*s0*t0 + 2*j1*l1*q0*s0*t0 - 
  2*i1*p1*q0*s0*t0 + 2*j1*l0*q1*s0*t0 + 
  2*j0*l1*q1*s0*t0 - 2*i1*p0*q1*s0*t0 - 
  2*i0*p1*q1*s0*t0 - 2*l1*m0*n0*s1*t0 - 
  2*l0*m1*n0*s1*t0 - 2*l0*m0*n1*s1*t0 + 
  2*j1*m0*p0*s1*t0 + 2*j0*m1*p0*s1*t0 + 
  2*j0*m0*p1*s1*t0 + 2*j1*l0*q0*s1*t0 + 
  2*j0*l1*q0*s1*t0 - 2*i1*p0*q0*s1*t0 - 
  2*i0*p1*q0*s1*t0 + 2*j0*l0*q1*s1*t0 - 
  2*i0*p0*q1*s1*t0 + l1*l1*n0*t0*t0 + 2*l0*l1*n1*t0*t0 - 
  2*j1*l1*p0*t0*t0 - 2*j1*l0*p1*t0*t0 - 
  2*j0*l1*p1*t0*t0 + 2*i1*p0*p1*t0*t0 + i0*p1*p1*t0*t0 + 
  2*l1*m0*o0*p0*t1 + 2*l0*m1*o0*p0*t1 + 
  2*l0*m0*o1*p0*t1 - 2*k1*m0*p0*p0*t1 - 
  2*k0*m1*p0*p0*t1 + 2*l0*m0*o0*p1*t1 - 
  4*k0*m0*p0*p1*t1 - 4*l0*l1*o0*q0*t1 - 
  2*l0*l0*o1*q0*t1 + 2*k1*l0*p0*q0*t1 + 
  2*k0*l1*p0*q0*t1 + 2*k0*l0*p1*q0*t1 - 
  2*l0*l0*o0*q1*t1 + 2*k0*l0*p0*q1*t1 - 
  2*l1*m0*n0*s0*t1 - 2*l0*m1*n0*s0*t1 - 
  2*l0*m0*n1*s0*t1 + 2*j1*m0*p0*s0*t1 + 
  2*j0*m1*p0*s0*t1 + 2*j0*m0*p1*s0*t1 + 
  2*j1*l0*q0*s0*t1 + 2*j0*l1*q0*s0*t1 - 
  2*i1*p0*q0*s0*t1 - 2*i0*p1*q0*s0*t1 + 
  2*j0*l0*q1*s0*t1 - 2*i0*p0*q1*s0*t1 - 
  2*l0*m0*n0*s1*t1 + 2*j0*m0*p0*s1*t1 + 
  2*j0*l0*q0*s1*t1 - 2*i0*p0*q0*s1*t1 + 
  4*l0*l1*n0*t0*t1 + 2*l0*l0*n1*t0*t1 - 
  4*j1*l0*p0*t0*t1 - 4*j0*l1*p0*t0*t1 + 
  2*i1*p0*p0*t0*t1 - 4*j0*l0*p1*t0*t1 + 
  4*i0*p0*p1*t0*t1 + l0*l0*n0*t1*t1 - 2*j0*l0*p0*t1*t1 + 
  i0*p0*p0*t1*t1 + m1*m1*o0*o0*u0 + 4*m0*m1*o0*o1*u0 + 
  m0*m0*o1*o1*u0 - 2*k1*m1*o0*q0*u0 - 
  2*k1*m0*o1*q0*u0 - 2*k0*m1*o1*q0*u0 + 
  k1*k1*q0*q0*u0 - 2*k1*m0*o0*q1*u0 - 
  2*k0*m1*o0*q1*u0 - 2*k0*m0*o1*q1*u0 + 
  4*k0*k1*q0*q1*u0 + k0*k0*q1*q1*u0 - m1*m1*n0*r0*u0 - 
  2*m0*m1*n1*r0*u0 + 2*j1*m1*q0*r0*u0 + 
  2*j1*m0*q1*r0*u0 + 2*j0*m1*q1*r0*u0 - 
  2*i1*q0*q1*r0*u0 - i0*q1*q1*r0*u0 - 
  2*m0*m1*n0*r1*u0 - m0*m0*n1*r1*u0 + 
  2*j1*m0*q0*r1*u0 + 2*j0*m1*q0*r1*u0 - 
  i1*q0*q0*r1*u0 + 2*j0*m0*q1*r1*u0 - 
  2*i0*q0*q1*r1*u0 + 2*k1*m1*n0*t0*u0 + 
  2*k1*m0*n1*t0*u0 + 2*k0*m1*n1*t0*u0 - 
  2*j1*m1*o0*t0*u0 - 2*j1*m0*o1*t0*u0 - 
  2*j0*m1*o1*t0*u0 - 2*j1*k1*q0*t0*u0 + 
  2*i1*o1*q0*t0*u0 - 2*j1*k0*q1*t0*u0 - 
  2*j0*k1*q1*t0*u0 + 2*i1*o0*q1*t0*u0 + 
  2*i0*o1*q1*t0*u0 + j1*j1*t0*t0*u0 - i1*n1*t0*t0*u0 + 
  2*k1*m0*n0*t1*u0 + 2*k0*m1*n0*t1*u0 + 
  2*k0*m0*n1*t1*u0 - 2*j1*m0*o0*t1*u0 - 
  2*j0*m1*o0*t1*u0 - 2*j0*m0*o1*t1*u0 - 
  2*j1*k0*q0*t1*u0 - 2*j0*k1*q0*t1*u0 + 
  2*i1*o0*q0*t1*u0 + 2*i0*o1*q0*t1*u0 - 
  2*j0*k0*q1*t1*u0 + 2*i0*o0*q1*t1*u0 + 
  4*j0*j1*t0*t1*u0 - 2*i1*n0*t0*t1*u0 - 
  2*i0*n1*t0*t1*u0 + j0*j0*t1*t1*u0 - i0*n0*t1*t1*u0 + 
  2*m0*m1*o0*o0*u1 + 2*m0*m0*o0*o1*u1 - 
  2*k1*m0*o0*q0*u1 - 2*k0*m1*o0*q0*u1 - 
  2*k0*m0*o1*q0*u1 + 2*k0*k1*q0*q0*u1 - 
  2*k0*m0*o0*q1*u1 + 2*k0*k0*q0*q1*u1 - 
  2*m0*m1*n0*r0*u1 - m0*m0*n1*r0*u1 + 
  2*j1*m0*q0*r0*u1 + 2*j0*m1*q0*r0*u1 - 
  i1*q0*q0*r0*u1 + 2*j0*m0*q1*r0*u1 - 
  2*i0*q0*q1*r0*u1 - m0*m0*n0*r1*u1 + 
  2*j0*m0*q0*r1*u1 - i0*q0*q0*r1*u1 + 
  2*k1*m0*n0*t0*u1 + 2*k0*m1*n0*t0*u1 + 
  2*k0*m0*n1*t0*u1 - 2*j1*m0*o0*t0*u1 - 
  2*j0*m1*o0*t0*u1 - 2*j0*m0*o1*t0*u1 - 
  2*j1*k0*q0*t0*u1 - 2*j0*k1*q0*t0*u1 + 
  2*i1*o0*q0*t0*u1 + 2*i0*o1*q0*t0*u1 - 
  2*j0*k0*q1*t0*u1 + 2*i0*o0*q1*t0*u1 + 
  2*j0*j1*t0*t0*u1 - i1*n0*t0*t0*u1 - i0*n1*t0*t0*u1 + 
  2*k0*m0*n0*t1*u1 - 2*j0*m0*o0*t1*u1 - 
  2*j0*k0*q0*t1*u1 + 2*i0*o0*q0*t1*u1 + 
  2*j0*j0*t0*t1*u1 - 2*i0*n0*t0*t1*u1 - 
  2*l1*m1*o0*o0*v0 - 4*l1*m0*o0*o1*v0 - 
  4*l0*m1*o0*o1*v0 - 2*l0*m0*o1*o1*v0 + 
  2*k1*m1*o0*p0*v0 + 2*k1*m0*o1*p0*v0 + 
  2*k0*m1*o1*p0*v0 + 2*k1*m0*o0*p1*v0 + 
  2*k0*m1*o0*p1*v0 + 2*k0*m0*o1*p1*v0 + 
  2*k1*l1*o0*q0*v0 + 2*k1*l0*o1*q0*v0 + 
  2*k0*l1*o1*q0*v0 - 2*k1*k1*p0*q0*v0 - 
  4*k0*k1*p1*q0*v0 + 2*k1*l0*o0*q1*v0 + 
  2*k0*l1*o0*q1*v0 + 2*k0*l0*o1*q1*v0 - 
  4*k0*k1*p0*q1*v0 - 2*k0*k0*p1*q1*v0 + 
  2*l1*m1*n0*r0*v0 + 2*l1*m0*n1*r0*v0 + 
  2*l0*m1*n1*r0*v0 - 2*j1*m1*p0*r0*v0 - 
  2*j1*m0*p1*r0*v0 - 2*j0*m1*p1*r0*v0 - 
  2*j1*l1*q0*r0*v0 + 2*i1*p1*q0*r0*v0 - 
  2*j1*l0*q1*r0*v0 - 2*j0*l1*q1*r0*v0 + 
  2*i1*p0*q1*r0*v0 + 2*i0*p1*q1*r0*v0 + 
  2*l1*m0*n0*r1*v0 + 2*l0*m1*n0*r1*v0 + 
  2*l0*m0*n1*r1*v0 - 2*j1*m0*p0*r1*v0 - 
  2*j0*m1*p0*r1*v0 - 2*j0*m0*p1*r1*v0 - 
  2*j1*l0*q0*r1*v0 - 2*j0*l1*q0*r1*v0 + 
  2*i1*p0*q0*r1*v0 + 2*i0*p1*q0*r1*v0 - 
  2*j0*l0*q1*r1*v0 + 2*i0*p0*q1*r1*v0 - 
  2*k1*m1*n0*s0*v0 - 2*k1*m0*n1*s0*v0 - 
  2*k0*m1*n1*s0*v0 + 2*j1*m1*o0*s0*v0 + 
  2*j1*m0*o1*s0*v0 + 2*j0*m1*o1*s0*v0 + 
  2*j1*k1*q0*s0*v0 - 2*i1*o1*q0*s0*v0 + 
  2*j1*k0*q1*s0*v0 + 2*j0*k1*q1*s0*v0 - 
  2*i1*o0*q1*s0*v0 - 2*i0*o1*q1*s0*v0 - 
  2*k1*m0*n0*s1*v0 - 2*k0*m1*n0*s1*v0 - 
  2*k0*m0*n1*s1*v0 + 2*j1*m0*o0*s1*v0 + 
  2*j0*m1*o0*s1*v0 + 2*j0*m0*o1*s1*v0 + 
  2*j1*k0*q0*s1*v0 + 2*j0*k1*q0*s1*v0 - 
  2*i1*o0*q0*s1*v0 - 2*i0*o1*q0*s1*v0 + 
  2*j0*k0*q1*s1*v0 - 2*i0*o0*q1*s1*v0 - 
  2*k1*l1*n0*t0*v0 - 2*k1*l0*n1*t0*v0 - 
  2*k0*l1*n1*t0*v0 + 2*j1*l1*o0*t0*v0 + 
  2*j1*l0*o1*t0*v0 + 2*j0*l1*o1*t0*v0 + 
  2*j1*k1*p0*t0*v0 - 2*i1*o1*p0*t0*v0 + 
  2*j1*k0*p1*t0*v0 + 2*j0*k1*p1*t0*v0 - 
  2*i1*o0*p1*t0*v0 - 2*i0*o1*p1*t0*v0 - 
  2*j1*j1*s0*t0*v0 + 2*i1*n1*s0*t0*v0 - 
  4*j0*j1*s1*t0*v0 + 2*i1*n0*s1*t0*v0 + 
  2*i0*n1*s1*t0*v0 - 2*k1*l0*n0*t1*v0 - 
  2*k0*l1*n0*t1*v0 - 2*k0*l0*n1*t1*v0 + 
  2*j1*l0*o0*t1*v0 + 2*j0*l1*o0*t1*v0 + 
  2*j0*l0*o1*t1*v0 + 2*j1*k0*p0*t1*v0 + 
  2*j0*k1*p0*t1*v0 - 2*i1*o0*p0*t1*v0 - 
  2*i0*o1*p0*t1*v0 + 2*j0*k0*p1*t1*v0 - 
  2*i0*o0*p1*t1*v0 - 4*j0*j1*s0*t1*v0 + 
  2*i1*n0*s0*t1*v0 + 2*i0*n1*s0*t1*v0 - 
  2*j0*j0*s1*t1*v0 + 2*i0*n0*s1*t1*v0 + k1*k1*n0*v0*v0 + 
  2*k0*k1*n1*v0*v0 - 2*j1*k1*o0*v0*v0 - 
  2*j1*k0*o1*v0*v0 - 2*j0*k1*o1*v0*v0 + 
  2*i1*o0*o1*v0*v0 + i0*o1*o1*v0*v0 + j1*j1*r0*v0*v0 - 
  i1*n1*r0*v0*v0 + 2*j0*j1*r1*v0*v0 - i1*n0*r1*v0*v0 - 
  i0*n1*r1*v0*v0 - 2*l1*m0*o0*o0*v1 - 2*l0*m1*o0*o0*v1 - 
  4*l0*m0*o0*o1*v1 + 2*k1*m0*o0*p0*v1 + 
  2*k0*m1*o0*p0*v1 + 2*k0*m0*o1*p0*v1 + 
  2*k0*m0*o0*p1*v1 + 2*k1*l0*o0*q0*v1 + 
  2*k0*l1*o0*q0*v1 + 2*k0*l0*o1*q0*v1 - 
  4*k0*k1*p0*q0*v1 - 2*k0*k0*p1*q0*v1 + 
  2*k0*l0*o0*q1*v1 - 2*k0*k0*p0*q1*v1 + 
  2*l1*m0*n0*r0*v1 + 2*l0*m1*n0*r0*v1 + 
  2*l0*m0*n1*r0*v1 - 2*j1*m0*p0*r0*v1 - 
  2*j0*m1*p0*r0*v1 - 2*j0*m0*p1*r0*v1 - 
  2*j1*l0*q0*r0*v1 - 2*j0*l1*q0*r0*v1 + 
  2*i1*p0*q0*r0*v1 + 2*i0*p1*q0*r0*v1 - 
  2*j0*l0*q1*r0*v1 + 2*i0*p0*q1*r0*v1 + 
  2*l0*m0*n0*r1*v1 - 2*j0*m0*p0*r1*v1 - 
  2*j0*l0*q0*r1*v1 + 2*i0*p0*q0*r1*v1 - 
  2*k1*m0*n0*s0*v1 - 2*k0*m1*n0*s0*v1 - 
  2*k0*m0*n1*s0*v1 + 2*j1*m0*o0*s0*v1 + 
  2*j0*m1*o0*s0*v1 + 2*j0*m0*o1*s0*v1 + 
  2*j1*k0*q0*s0*v1 + 2*j0*k1*q0*s0*v1 - 
  2*i1*o0*q0*s0*v1 - 2*i0*o1*q0*s0*v1 + 
  2*j0*k0*q1*s0*v1 - 2*i0*o0*q1*s0*v1 - 
  2*k0*m0*n0*s1*v1 + 2*j0*m0*o0*s1*v1 + 
  2*j0*k0*q0*s1*v1 - 2*i0*o0*q0*s1*v1 - 
  2*k1*l0*n0*t0*v1 - 2*k0*l1*n0*t0*v1 - 
  2*k0*l0*n1*t0*v1 + 2*j1*l0*o0*t0*v1 + 
  2*j0*l1*o0*t0*v1 + 2*j0*l0*o1*t0*v1 + 
  2*j1*k0*p0*t0*v1 + 2*j0*k1*p0*t0*v1 - 
  2*i1*o0*p0*t0*v1 - 2*i0*o1*p0*t0*v1 + 
  2*j0*k0*p1*t0*v1 - 2*i0*o0*p1*t0*v1 - 
  4*j0*j1*s0*t0*v1 + 2*i1*n0*s0*t0*v1 + 
  2*i0*n1*s0*t0*v1 - 2*j0*j0*s1*t0*v1 + 
  2*i0*n0*s1*t0*v1 - 2*k0*l0*n0*t1*v1 + 
  2*j0*l0*o0*t1*v1 + 2*j0*k0*p0*t1*v1 - 
  2*i0*o0*p0*t1*v1 - 2*j0*j0*s0*t1*v1 + 
  2*i0*n0*s0*t1*v1 + 4*k0*k1*n0*v0*v1 + 
  2*k0*k0*n1*v0*v1 - 4*j1*k0*o0*v0*v1 - 
  4*j0*k1*o0*v0*v1 + 2*i1*o0*o0*v0*v1 - 
  4*j0*k0*o1*v0*v1 + 4*i0*o0*o1*v0*v1 + 
  4*j0*j1*r0*v0*v1 - 2*i1*n0*r0*v0*v1 - 
  2*i0*n1*r0*v0*v1 + 2*j0*j0*r1*v0*v1 - 
  2*i0*n0*r1*v0*v1 + k0*k0*n0*v1*v1 - 2*j0*k0*o0*v1*v1 + 
  i0*o0*o0*v1*v1 + j0*j0*r0*v1*v1 - i0*n0*r0*v1*v1 + 
  l1*l1*o0*o0*w0 + 4*l0*l1*o0*o1*w0 + l0*l0*o1*o1*w0 - 
  2*k1*l1*o0*p0*w0 - 2*k1*l0*o1*p0*w0 - 
  2*k0*l1*o1*p0*w0 + k1*k1*p0*p0*w0 - 
  2*k1*l0*o0*p1*w0 - 2*k0*l1*o0*p1*w0 - 
  2*k0*l0*o1*p1*w0 + 4*k0*k1*p0*p1*w0 + 
  k0*k0*p1*p1*w0 - l1*l1*n0*r0*w0 - 2*l0*l1*n1*r0*w0 + 
  2*j1*l1*p0*r0*w0 + 2*j1*l0*p1*r0*w0 + 
  2*j0*l1*p1*r0*w0 - 2*i1*p0*p1*r0*w0 - 
  i0*p1*p1*r0*w0 - 2*l0*l1*n0*r1*w0 - l0*l0*n1*r1*w0 + 
  2*j1*l0*p0*r1*w0 + 2*j0*l1*p0*r1*w0 - 
  i1*p0*p0*r1*w0 + 2*j0*l0*p1*r1*w0 - 
  2*i0*p0*p1*r1*w0 + 2*k1*l1*n0*s0*w0 + 
  2*k1*l0*n1*s0*w0 + 2*k0*l1*n1*s0*w0 - 
  2*j1*l1*o0*s0*w0 - 2*j1*l0*o1*s0*w0 - 
  2*j0*l1*o1*s0*w0 - 2*j1*k1*p0*s0*w0 + 
  2*i1*o1*p0*s0*w0 - 2*j1*k0*p1*s0*w0 - 
  2*j0*k1*p1*s0*w0 + 2*i1*o0*p1*s0*w0 + 
  2*i0*o1*p1*s0*w0 + j1*j1*s0*s0*w0 - i1*n1*s0*s0*w0 + 
  2*k1*l0*n0*s1*w0 + 2*k0*l1*n0*s1*w0 + 
  2*k0*l0*n1*s1*w0 - 2*j1*l0*o0*s1*w0 - 
  2*j0*l1*o0*s1*w0 - 2*j0*l0*o1*s1*w0 - 
  2*j1*k0*p0*s1*w0 - 2*j0*k1*p0*s1*w0 + 
  2*i1*o0*p0*s1*w0 + 2*i0*o1*p0*s1*w0 - 
  2*j0*k0*p1*s1*w0 + 2*i0*o0*p1*s1*w0 + 
  4*j0*j1*s0*s1*w0 - 2*i1*n0*s0*s1*w0 - 
  2*i0*n1*s0*s1*w0 + j0*j0*s1*s1*w0 - i0*n0*s1*s1*w0 - 
  k1*k1*n0*u0*w0 - 2*k0*k1*n1*u0*w0 + 
  2*j1*k1*o0*u0*w0 + 2*j1*k0*o1*u0*w0 + 
  2*j0*k1*o1*u0*w0 - 2*i1*o0*o1*u0*w0 - 
  i0*o1*o1*u0*w0 - j1*j1*r0*u0*w0 + i1*n1*r0*u0*w0 - 
  2*j0*j1*r1*u0*w0 + i1*n0*r1*u0*w0 + 
  i0*n1*r1*u0*w0 - 2*k0*k1*n0*u1*w0 - k0*k0*n1*u1*w0 + 
  2*j1*k0*o0*u1*w0 + 2*j0*k1*o0*u1*w0 - 
  i1*o0*o0*u1*w0 + 2*j0*k0*o1*u1*w0 - 
  2*i0*o0*o1*u1*w0 - 2*j0*j1*r0*u1*w0 + 
  i1*n0*r0*u1*w0 + i0*n1*r0*u1*w0 - j0*j0*r1*u1*w0 + 
  i0*n0*r1*u1*w0 + 2*l0*l1*o0*o0*w1 + 
  2*l0*l0*o0*o1*w1 - 2*k1*l0*o0*p0*w1 - 
  2*k0*l1*o0*p0*w1 - 2*k0*l0*o1*p0*w1 + 
  2*k0*k1*p0*p0*w1 - 2*k0*l0*o0*p1*w1 + 
  2*k0*k0*p0*p1*w1 - 2*l0*l1*n0*r0*w1 - 
  l0*l0*n1*r0*w1 + 2*j1*l0*p0*r0*w1 + 
  2*j0*l1*p0*r0*w1 - i1*p0*p0*r0*w1 + 
  2*j0*l0*p1*r0*w1 - 2*i0*p0*p1*r0*w1 - 
  l0*l0*n0*r1*w1 + 2*j0*l0*p0*r1*w1 - i0*p0*p0*r1*w1 + 
  2*k1*l0*n0*s0*w1 + 2*k0*l1*n0*s0*w1 + 
  2*k0*l0*n1*s0*w1 - 2*j1*l0*o0*s0*w1 - 
  2*j0*l1*o0*s0*w1 - 2*j0*l0*o1*s0*w1 - 
  2*j1*k0*p0*s0*w1 - 2*j0*k1*p0*s0*w1 + 
  2*i1*o0*p0*s0*w1 + 2*i0*o1*p0*s0*w1 - 
  2*j0*k0*p1*s0*w1 + 2*i0*o0*p1*s0*w1 + 
  2*j0*j1*s0*s0*w1 - i1*n0*s0*s0*w1 - i0*n1*s0*s0*w1 + 
  2*k0*l0*n0*s1*w1 - 2*j0*l0*o0*s1*w1 - 
  2*j0*k0*p0*s1*w1 + 2*i0*o0*p0*s1*w1 + 
  2*j0*j0*s0*s1*w1 - 2*i0*n0*s0*s1*w1 - 
  2*k0*k1*n0*u0*w1 - k0*k0*n1*u0*w1 + 
  2*j1*k0*o0*u0*w1 + 2*j0*k1*o0*u0*w1 - 
  i1*o0*o0*u0*w1 + 2*j0*k0*o1*u0*w1 - 
  2*i0*o0*o1*u0*w1 - 2*j0*j1*r0*u0*w1 + 
  i1*n0*r0*u0*w1 + i0*n1*r0*u0*w1 - j0*j0*r1*u0*w1 + 
  i0*n0*r1*u0*w1 - k0*k0*n0*u1*w1 + 2*j0*k0*o0*u1*w1 - 
  i0*o0*o0*u1*w1 - j0*j0*r0*u1*w1 + i0*n0*r0*u1*w1;
   
    k[18] = 2*m1*m1*p0*p1*r0 + 2*m0*m1*p1*p1*r0 - 
  2*l1*m1*p1*q0*r0 - 2*l1*m1*p0*q1*r0 - 
  2*l1*m0*p1*q1*r0 - 2*l0*m1*p1*q1*r0 + 
  2*l1*l1*q0*q1*r0 + 2*l0*l1*q1*q1*r0 + m1*m1*p0*p0*r1 + 
  4*m0*m1*p0*p1*r1 + m0*m0*p1*p1*r1 - 
  2*l1*m1*p0*q0*r1 - 2*l1*m0*p1*q0*r1 - 
  2*l0*m1*p1*q0*r1 + l1*l1*q0*q0*r1 - 
  2*l1*m0*p0*q1*r1 - 2*l0*m1*p0*q1*r1 - 
  2*l0*m0*p1*q1*r1 + 4*l0*l1*q0*q1*r1 + 
  l0*l0*q1*q1*r1 - 2*m1*m1*o1*p0*s0 - 2*m1*m1*o0*p1*s0 - 
  4*m0*m1*o1*p1*s0 + 2*l1*m1*o1*q0*s0 + 
  2*k1*m1*p1*q0*s0 + 2*l1*m1*o0*q1*s0 + 
  2*l1*m0*o1*q1*s0 + 2*l0*m1*o1*q1*s0 + 
  2*k1*m1*p0*q1*s0 + 2*k1*m0*p1*q1*s0 + 
  2*k0*m1*p1*q1*s0 - 4*k1*l1*q0*q1*s0 - 
  2*k1*l0*q1*q1*s0 - 2*k0*l1*q1*q1*s0 + m1*m1*n1*s0*s0 - 
  2*j1*m1*q1*s0*s0 + i1*q1*q1*s0*s0 - 2*m1*m1*o0*p0*s1 - 
  4*m0*m1*o1*p0*s1 - 4*m0*m1*o0*p1*s1 - 
  2*m0*m0*o1*p1*s1 + 2*l1*m1*o0*q0*s1 + 
  2*l1*m0*o1*q0*s1 + 2*l0*m1*o1*q0*s1 + 
  2*k1*m1*p0*q0*s1 + 2*k1*m0*p1*q0*s1 + 
  2*k0*m1*p1*q0*s1 - 2*k1*l1*q0*q0*s1 + 
  2*l1*m0*o0*q1*s1 + 2*l0*m1*o0*q1*s1 + 
  2*l0*m0*o1*q1*s1 + 2*k1*m0*p0*q1*s1 + 
  2*k0*m1*p0*q1*s1 + 2*k0*m0*p1*q1*s1 - 
  4*k1*l0*q0*q1*s1 - 4*k0*l1*q0*q1*s1 - 
  2*k0*l0*q1*q1*s1 + 2*m1*m1*n0*s0*s1 + 
  4*m0*m1*n1*s0*s1 - 4*j1*m1*q0*s0*s1 - 
  4*j1*m0*q1*s0*s1 - 4*j0*m1*q1*s0*s1 + 
  4*i1*q0*q1*s0*s1 + 2*i0*q1*q1*s0*s1 + 
  2*m0*m1*n0*s1*s1 + m0*m0*n1*s1*s1 - 2*j1*m0*q0*s1*s1 - 
  2*j0*m1*q0*s1*s1 + i1*q0*q0*s1*s1 - 2*j0*m0*q1*s1*s1 + 
  2*i0*q0*q1*s1*s1 + 2*l1*m1*o1*p0*t0 + 
  2*l1*m1*o0*p1*t0 + 2*l1*m0*o1*p1*t0 + 
  2*l0*m1*o1*p1*t0 - 4*k1*m1*p0*p1*t0 - 
  2*k1*m0*p1*p1*t0 - 2*k0*m1*p1*p1*t0 - 
  2*l1*l1*o1*q0*t0 + 2*k1*l1*p1*q0*t0 - 
  2*l1*l1*o0*q1*t0 - 4*l0*l1*o1*q1*t0 + 
  2*k1*l1*p0*q1*t0 + 2*k1*l0*p1*q1*t0 + 
  2*k0*l1*p1*q1*t0 - 2*l1*m1*n1*s0*t0 + 
  2*j1*m1*p1*s0*t0 + 2*j1*l1*q1*s0*t0 - 
  2*i1*p1*q1*s0*t0 - 2*l1*m1*n0*s1*t0 - 
  2*l1*m0*n1*s1*t0 - 2*l0*m1*n1*s1*t0 + 
  2*j1*m1*p0*s1*t0 + 2*j1*m0*p1*s1*t0 + 
  2*j0*m1*p1*s1*t0 + 2*j1*l1*q0*s1*t0 - 
  2*i1*p1*q0*s1*t0 + 2*j1*l0*q1*s1*t0 + 
  2*j0*l1*q1*s1*t0 - 2*i1*p0*q1*s1*t0 - 
  2*i0*p1*q1*s1*t0 + l1*l1*n1*t0*t0 - 2*j1*l1*p1*t0*t0 + 
  i1*p1*p1*t0*t0 + 2*l1*m1*o0*p0*t1 + 
  2*l1*m0*o1*p0*t1 + 2*l0*m1*o1*p0*t1 - 
  2*k1*m1*p0*p0*t1 + 2*l1*m0*o0*p1*t1 + 
  2*l0*m1*o0*p1*t1 + 2*l0*m0*o1*p1*t1 - 
  4*k1*m0*p0*p1*t1 - 4*k0*m1*p0*p1*t1 - 
  2*k0*m0*p1*p1*t1 - 2*l1*l1*o0*q0*t1 - 
  4*l0*l1*o1*q0*t1 + 2*k1*l1*p0*q0*t1 + 
  2*k1*l0*p1*q0*t1 + 2*k0*l1*p1*q0*t1 - 
  4*l0*l1*o0*q1*t1 - 2*l0*l0*o1*q1*t1 + 
  2*k1*l0*p0*q1*t1 + 2*k0*l1*p0*q1*t1 + 
  2*k0*l0*p1*q1*t1 - 2*l1*m1*n0*s0*t1 - 
  2*l1*m0*n1*s0*t1 - 2*l0*m1*n1*s0*t1 + 
  2*j1*m1*p0*s0*t1 + 2*j1*m0*p1*s0*t1 + 
  2*j0*m1*p1*s0*t1 + 2*j1*l1*q0*s0*t1 - 
  2*i1*p1*q0*s0*t1 + 2*j1*l0*q1*s0*t1 + 
  2*j0*l1*q1*s0*t1 - 2*i1*p0*q1*s0*t1 - 
  2*i0*p1*q1*s0*t1 - 2*l1*m0*n0*s1*t1 - 
  2*l0*m1*n0*s1*t1 - 2*l0*m0*n1*s1*t1 + 
  2*j1*m0*p0*s1*t1 + 2*j0*m1*p0*s1*t1 + 
  2*j0*m0*p1*s1*t1 + 2*j1*l0*q0*s1*t1 + 
  2*j0*l1*q0*s1*t1 - 2*i1*p0*q0*s1*t1 - 
  2*i0*p1*q0*s1*t1 + 2*j0*l0*q1*s1*t1 - 
  2*i0*p0*q1*s1*t1 + 2*l1*l1*n0*t0*t1 + 
  4*l0*l1*n1*t0*t1 - 4*j1*l1*p0*t0*t1 - 
  4*j1*l0*p1*t0*t1 - 4*j0*l1*p1*t0*t1 + 
  4*i1*p0*p1*t0*t1 + 2*i0*p1*p1*t0*t1 + 
  2*l0*l1*n0*t1*t1 + l0*l0*n1*t1*t1 - 2*j1*l0*p0*t1*t1 - 
  2*j0*l1*p0*t1*t1 + i1*p0*p0*t1*t1 - 2*j0*l0*p1*t1*t1 + 
  2*i0*p0*p1*t1*t1 + 2*m1*m1*o0*o1*u0 + 
  2*m0*m1*o1*o1*u0 - 2*k1*m1*o1*q0*u0 - 
  2*k1*m1*o0*q1*u0 - 2*k1*m0*o1*q1*u0 - 
  2*k0*m1*o1*q1*u0 + 2*k1*k1*q0*q1*u0 + 
  2*k0*k1*q1*q1*u0 - m1*m1*n1*r0*u0 + 
  2*j1*m1*q1*r0*u0 - i1*q1*q1*r0*u0 - m1*m1*n0*r1*u0 - 
  2*m0*m1*n1*r1*u0 + 2*j1*m1*q0*r1*u0 + 
  2*j1*m0*q1*r1*u0 + 2*j0*m1*q1*r1*u0 - 
  2*i1*q0*q1*r1*u0 - i0*q1*q1*r1*u0 + 
  2*k1*m1*n1*t0*u0 - 2*j1*m1*o1*t0*u0 - 
  2*j1*k1*q1*t0*u0 + 2*i1*o1*q1*t0*u0 + 
  2*k1*m1*n0*t1*u0 + 2*k1*m0*n1*t1*u0 + 
  2*k0*m1*n1*t1*u0 - 2*j1*m1*o0*t1*u0 - 
  2*j1*m0*o1*t1*u0 - 2*j0*m1*o1*t1*u0 - 
  2*j1*k1*q0*t1*u0 + 2*i1*o1*q0*t1*u0 - 
  2*j1*k0*q1*t1*u0 - 2*j0*k1*q1*t1*u0 + 
  2*i1*o0*q1*t1*u0 + 2*i0*o1*q1*t1*u0 + 
  2*j1*j1*t0*t1*u0 - 2*i1*n1*t0*t1*u0 + 
  2*j0*j1*t1*t1*u0 - i1*n0*t1*t1*u0 - i0*n1*t1*t1*u0 + 
  m1*m1*o0*o0*u1 + 4*m0*m1*o0*o1*u1 + m0*m0*o1*o1*u1 - 
  2*k1*m1*o0*q0*u1 - 2*k1*m0*o1*q0*u1 - 
  2*k0*m1*o1*q0*u1 + k1*k1*q0*q0*u1 - 
  2*k1*m0*o0*q1*u1 - 2*k0*m1*o0*q1*u1 - 
  2*k0*m0*o1*q1*u1 + 4*k0*k1*q0*q1*u1 + 
  k0*k0*q1*q1*u1 - m1*m1*n0*r0*u1 - 2*m0*m1*n1*r0*u1 + 
  2*j1*m1*q0*r0*u1 + 2*j1*m0*q1*r0*u1 + 
  2*j0*m1*q1*r0*u1 - 2*i1*q0*q1*r0*u1 - 
  i0*q1*q1*r0*u1 - 2*m0*m1*n0*r1*u1 - m0*m0*n1*r1*u1 + 
  2*j1*m0*q0*r1*u1 + 2*j0*m1*q0*r1*u1 - 
  i1*q0*q0*r1*u1 + 2*j0*m0*q1*r1*u1 - 
  2*i0*q0*q1*r1*u1 + 2*k1*m1*n0*t0*u1 + 
  2*k1*m0*n1*t0*u1 + 2*k0*m1*n1*t0*u1 - 
  2*j1*m1*o0*t0*u1 - 2*j1*m0*o1*t0*u1 - 
  2*j0*m1*o1*t0*u1 - 2*j1*k1*q0*t0*u1 + 
  2*i1*o1*q0*t0*u1 - 2*j1*k0*q1*t0*u1 - 
  2*j0*k1*q1*t0*u1 + 2*i1*o0*q1*t0*u1 + 
  2*i0*o1*q1*t0*u1 + j1*j1*t0*t0*u1 - i1*n1*t0*t0*u1 + 
  2*k1*m0*n0*t1*u1 + 2*k0*m1*n0*t1*u1 + 
  2*k0*m0*n1*t1*u1 - 2*j1*m0*o0*t1*u1 - 
  2*j0*m1*o0*t1*u1 - 2*j0*m0*o1*t1*u1 - 
  2*j1*k0*q0*t1*u1 - 2*j0*k1*q0*t1*u1 + 
  2*i1*o0*q0*t1*u1 + 2*i0*o1*q0*t1*u1 - 
  2*j0*k0*q1*t1*u1 + 2*i0*o0*q1*t1*u1 + 
  4*j0*j1*t0*t1*u1 - 2*i1*n0*t0*t1*u1 - 
  2*i0*n1*t0*t1*u1 + j0*j0*t1*t1*u1 - i0*n0*t1*t1*u1 - 
  4*l1*m1*o0*o1*v0 - 2*l1*m0*o1*o1*v0 - 
  2*l0*m1*o1*o1*v0 + 2*k1*m1*o1*p0*v0 + 
  2*k1*m1*o0*p1*v0 + 2*k1*m0*o1*p1*v0 + 
  2*k0*m1*o1*p1*v0 + 2*k1*l1*o1*q0*v0 - 
  2*k1*k1*p1*q0*v0 + 2*k1*l1*o0*q1*v0 + 
  2*k1*l0*o1*q1*v0 + 2*k0*l1*o1*q1*v0 - 
  2*k1*k1*p0*q1*v0 - 4*k0*k1*p1*q1*v0 + 
  2*l1*m1*n1*r0*v0 - 2*j1*m1*p1*r0*v0 - 
  2*j1*l1*q1*r0*v0 + 2*i1*p1*q1*r0*v0 + 
  2*l1*m1*n0*r1*v0 + 2*l1*m0*n1*r1*v0 + 
  2*l0*m1*n1*r1*v0 - 2*j1*m1*p0*r1*v0 - 
  2*j1*m0*p1*r1*v0 - 2*j0*m1*p1*r1*v0 - 
  2*j1*l1*q0*r1*v0 + 2*i1*p1*q0*r1*v0 - 
  2*j1*l0*q1*r1*v0 - 2*j0*l1*q1*r1*v0 + 
  2*i1*p0*q1*r1*v0 + 2*i0*p1*q1*r1*v0 - 
  2*k1*m1*n1*s0*v0 + 2*j1*m1*o1*s0*v0 + 
  2*j1*k1*q1*s0*v0 - 2*i1*o1*q1*s0*v0 - 
  2*k1*m1*n0*s1*v0 - 2*k1*m0*n1*s1*v0 - 
  2*k0*m1*n1*s1*v0 + 2*j1*m1*o0*s1*v0 + 
  2*j1*m0*o1*s1*v0 + 2*j0*m1*o1*s1*v0 + 
  2*j1*k1*q0*s1*v0 - 2*i1*o1*q0*s1*v0 + 
  2*j1*k0*q1*s1*v0 + 2*j0*k1*q1*s1*v0 - 
  2*i1*o0*q1*s1*v0 - 2*i0*o1*q1*s1*v0 - 
  2*k1*l1*n1*t0*v0 + 2*j1*l1*o1*t0*v0 + 
  2*j1*k1*p1*t0*v0 - 2*i1*o1*p1*t0*v0 - 
  2*j1*j1*s1*t0*v0 + 2*i1*n1*s1*t0*v0 - 
  2*k1*l1*n0*t1*v0 - 2*k1*l0*n1*t1*v0 - 
  2*k0*l1*n1*t1*v0 + 2*j1*l1*o0*t1*v0 + 
  2*j1*l0*o1*t1*v0 + 2*j0*l1*o1*t1*v0 + 
  2*j1*k1*p0*t1*v0 - 2*i1*o1*p0*t1*v0 + 
  2*j1*k0*p1*t1*v0 + 2*j0*k1*p1*t1*v0 - 
  2*i1*o0*p1*t1*v0 - 2*i0*o1*p1*t1*v0 - 
  2*j1*j1*s0*t1*v0 + 2*i1*n1*s0*t1*v0 - 
  4*j0*j1*s1*t1*v0 + 2*i1*n0*s1*t1*v0 + 
  2*i0*n1*s1*t1*v0 + k1*k1*n1*v0*v0 - 2*j1*k1*o1*v0*v0 + 
  i1*o1*o1*v0*v0 + j1*j1*r1*v0*v0 - i1*n1*r1*v0*v0 - 
  2*l1*m1*o0*o0*v1 - 4*l1*m0*o0*o1*v1 - 
  4*l0*m1*o0*o1*v1 - 2*l0*m0*o1*o1*v1 + 
  2*k1*m1*o0*p0*v1 + 2*k1*m0*o1*p0*v1 + 
  2*k0*m1*o1*p0*v1 + 2*k1*m0*o0*p1*v1 + 
  2*k0*m1*o0*p1*v1 + 2*k0*m0*o1*p1*v1 + 
  2*k1*l1*o0*q0*v1 + 2*k1*l0*o1*q0*v1 + 
  2*k0*l1*o1*q0*v1 - 2*k1*k1*p0*q0*v1 - 
  4*k0*k1*p1*q0*v1 + 2*k1*l0*o0*q1*v1 + 
  2*k0*l1*o0*q1*v1 + 2*k0*l0*o1*q1*v1 - 
  4*k0*k1*p0*q1*v1 - 2*k0*k0*p1*q1*v1 + 
  2*l1*m1*n0*r0*v1 + 2*l1*m0*n1*r0*v1 + 
  2*l0*m1*n1*r0*v1 - 2*j1*m1*p0*r0*v1 - 
  2*j1*m0*p1*r0*v1 - 2*j0*m1*p1*r0*v1 - 
  2*j1*l1*q0*r0*v1 + 2*i1*p1*q0*r0*v1 - 
  2*j1*l0*q1*r0*v1 - 2*j0*l1*q1*r0*v1 + 
  2*i1*p0*q1*r0*v1 + 2*i0*p1*q1*r0*v1 + 
  2*l1*m0*n0*r1*v1 + 2*l0*m1*n0*r1*v1 + 
  2*l0*m0*n1*r1*v1 - 2*j1*m0*p0*r1*v1 - 
  2*j0*m1*p0*r1*v1 - 2*j0*m0*p1*r1*v1 - 
  2*j1*l0*q0*r1*v1 - 2*j0*l1*q0*r1*v1 + 
  2*i1*p0*q0*r1*v1 + 2*i0*p1*q0*r1*v1 - 
  2*j0*l0*q1*r1*v1 + 2*i0*p0*q1*r1*v1 - 
  2*k1*m1*n0*s0*v1 - 2*k1*m0*n1*s0*v1 - 
  2*k0*m1*n1*s0*v1 + 2*j1*m1*o0*s0*v1 + 
  2*j1*m0*o1*s0*v1 + 2*j0*m1*o1*s0*v1 + 
  2*j1*k1*q0*s0*v1 - 2*i1*o1*q0*s0*v1 + 
  2*j1*k0*q1*s0*v1 + 2*j0*k1*q1*s0*v1 - 
  2*i1*o0*q1*s0*v1 - 2*i0*o1*q1*s0*v1 - 
  2*k1*m0*n0*s1*v1 - 2*k0*m1*n0*s1*v1 - 
  2*k0*m0*n1*s1*v1 + 2*j1*m0*o0*s1*v1 + 
  2*j0*m1*o0*s1*v1 + 2*j0*m0*o1*s1*v1 + 
  2*j1*k0*q0*s1*v1 + 2*j0*k1*q0*s1*v1 - 
  2*i1*o0*q0*s1*v1 - 2*i0*o1*q0*s1*v1 + 
  2*j0*k0*q1*s1*v1 - 2*i0*o0*q1*s1*v1 - 
  2*k1*l1*n0*t0*v1 - 2*k1*l0*n1*t0*v1 - 
  2*k0*l1*n1*t0*v1 + 2*j1*l1*o0*t0*v1 + 
  2*j1*l0*o1*t0*v1 + 2*j0*l1*o1*t0*v1 + 
  2*j1*k1*p0*t0*v1 - 2*i1*o1*p0*t0*v1 + 
  2*j1*k0*p1*t0*v1 + 2*j0*k1*p1*t0*v1 - 
  2*i1*o0*p1*t0*v1 - 2*i0*o1*p1*t0*v1 - 
  2*j1*j1*s0*t0*v1 + 2*i1*n1*s0*t0*v1 - 
  4*j0*j1*s1*t0*v1 + 2*i1*n0*s1*t0*v1 + 
  2*i0*n1*s1*t0*v1 - 2*k1*l0*n0*t1*v1 - 
  2*k0*l1*n0*t1*v1 - 2*k0*l0*n1*t1*v1 + 
  2*j1*l0*o0*t1*v1 + 2*j0*l1*o0*t1*v1 + 
  2*j0*l0*o1*t1*v1 + 2*j1*k0*p0*t1*v1 + 
  2*j0*k1*p0*t1*v1 - 2*i1*o0*p0*t1*v1 - 
  2*i0*o1*p0*t1*v1 + 2*j0*k0*p1*t1*v1 - 
  2*i0*o0*p1*t1*v1 - 4*j0*j1*s0*t1*v1 + 
  2*i1*n0*s0*t1*v1 + 2*i0*n1*s0*t1*v1 - 
  2*j0*j0*s1*t1*v1 + 2*i0*n0*s1*t1*v1 + 
  2*k1*k1*n0*v0*v1 + 4*k0*k1*n1*v0*v1 - 
  4*j1*k1*o0*v0*v1 - 4*j1*k0*o1*v0*v1 - 
  4*j0*k1*o1*v0*v1 + 4*i1*o0*o1*v0*v1 + 
  2*i0*o1*o1*v0*v1 + 2*j1*j1*r0*v0*v1 - 
  2*i1*n1*r0*v0*v1 + 4*j0*j1*r1*v0*v1 - 
  2*i1*n0*r1*v0*v1 - 2*i0*n1*r1*v0*v1 + 
  2*k0*k1*n0*v1*v1 + k0*k0*n1*v1*v1 - 2*j1*k0*o0*v1*v1 - 
  2*j0*k1*o0*v1*v1 + i1*o0*o0*v1*v1 - 2*j0*k0*o1*v1*v1 + 
  2*i0*o0*o1*v1*v1 + 2*j0*j1*r0*v1*v1 - i1*n0*r0*v1*v1 - 
  i0*n1*r0*v1*v1 + j0*j0*r1*v1*v1 - i0*n0*r1*v1*v1 + 
  2*l1*l1*o0*o1*w0 + 2*l0*l1*o1*o1*w0 - 
  2*k1*l1*o1*p0*w0 - 2*k1*l1*o0*p1*w0 - 
  2*k1*l0*o1*p1*w0 - 2*k0*l1*o1*p1*w0 + 
  2*k1*k1*p0*p1*w0 + 2*k0*k1*p1*p1*w0 - l1*l1*n1*r0*w0 + 
  2*j1*l1*p1*r0*w0 - i1*p1*p1*r0*w0 - l1*l1*n0*r1*w0 - 
  2*l0*l1*n1*r1*w0 + 2*j1*l1*p0*r1*w0 + 
  2*j1*l0*p1*r1*w0 + 2*j0*l1*p1*r1*w0 - 
  2*i1*p0*p1*r1*w0 - i0*p1*p1*r1*w0 + 
  2*k1*l1*n1*s0*w0 - 2*j1*l1*o1*s0*w0 - 
  2*j1*k1*p1*s0*w0 + 2*i1*o1*p1*s0*w0 + 
  2*k1*l1*n0*s1*w0 + 2*k1*l0*n1*s1*w0 + 
  2*k0*l1*n1*s1*w0 - 2*j1*l1*o0*s1*w0 - 
  2*j1*l0*o1*s1*w0 - 2*j0*l1*o1*s1*w0 - 
  2*j1*k1*p0*s1*w0 + 2*i1*o1*p0*s1*w0 - 
  2*j1*k0*p1*s1*w0 - 2*j0*k1*p1*s1*w0 + 
  2*i1*o0*p1*s1*w0 + 2*i0*o1*p1*s1*w0 + 
  2*j1*j1*s0*s1*w0 - 2*i1*n1*s0*s1*w0 + 
  2*j0*j1*s1*s1*w0 - i1*n0*s1*s1*w0 - i0*n1*s1*s1*w0 - 
  k1*k1*n1*u0*w0 + 2*j1*k1*o1*u0*w0 - i1*o1*o1*u0*w0 - 
  j1*j1*r1*u0*w0 + i1*n1*r1*u0*w0 - k1*k1*n0*u1*w0 - 
  2*k0*k1*n1*u1*w0 + 2*j1*k1*o0*u1*w0 + 
  2*j1*k0*o1*u1*w0 + 2*j0*k1*o1*u1*w0 - 
  2*i1*o0*o1*u1*w0 - i0*o1*o1*u1*w0 - j1*j1*r0*u1*w0 + 
  i1*n1*r0*u1*w0 - 2*j0*j1*r1*u1*w0 + 
  i1*n0*r1*u1*w0 + i0*n1*r1*u1*w0 + l1*l1*o0*o0*w1 + 
  4*l0*l1*o0*o1*w1 + l0*l0*o1*o1*w1 - 
  2*k1*l1*o0*p0*w1 - 2*k1*l0*o1*p0*w1 - 
  2*k0*l1*o1*p0*w1 + k1*k1*p0*p0*w1 - 
  2*k1*l0*o0*p1*w1 - 2*k0*l1*o0*p1*w1 - 
  2*k0*l0*o1*p1*w1 + 4*k0*k1*p0*p1*w1 + 
  k0*k0*p1*p1*w1 - l1*l1*n0*r0*w1 - 2*l0*l1*n1*r0*w1 + 
  2*j1*l1*p0*r0*w1 + 2*j1*l0*p1*r0*w1 + 
  2*j0*l1*p1*r0*w1 - 2*i1*p0*p1*r0*w1 - 
  i0*p1*p1*r0*w1 - 2*l0*l1*n0*r1*w1 - l0*l0*n1*r1*w1 + 
  2*j1*l0*p0*r1*w1 + 2*j0*l1*p0*r1*w1 - 
  i1*p0*p0*r1*w1 + 2*j0*l0*p1*r1*w1 - 
  2*i0*p0*p1*r1*w1 + 2*k1*l1*n0*s0*w1 + 
  2*k1*l0*n1*s0*w1 + 2*k0*l1*n1*s0*w1 - 
  2*j1*l1*o0*s0*w1 - 2*j1*l0*o1*s0*w1 - 
  2*j0*l1*o1*s0*w1 - 2*j1*k1*p0*s0*w1 + 
  2*i1*o1*p0*s0*w1 - 2*j1*k0*p1*s0*w1 - 
  2*j0*k1*p1*s0*w1 + 2*i1*o0*p1*s0*w1 + 
  2*i0*o1*p1*s0*w1 + j1*j1*s0*s0*w1 - i1*n1*s0*s0*w1 + 
  2*k1*l0*n0*s1*w1 + 2*k0*l1*n0*s1*w1 + 
  2*k0*l0*n1*s1*w1 - 2*j1*l0*o0*s1*w1 - 
  2*j0*l1*o0*s1*w1 - 2*j0*l0*o1*s1*w1 - 
  2*j1*k0*p0*s1*w1 - 2*j0*k1*p0*s1*w1 + 
  2*i1*o0*p0*s1*w1 + 2*i0*o1*p0*s1*w1 - 
  2*j0*k0*p1*s1*w1 + 2*i0*o0*p1*s1*w1 + 
  4*j0*j1*s0*s1*w1 - 2*i1*n0*s0*s1*w1 - 
  2*i0*n1*s0*s1*w1 + j0*j0*s1*s1*w1 - i0*n0*s1*s1*w1 - 
  k1*k1*n0*u0*w1 - 2*k0*k1*n1*u0*w1 + 
  2*j1*k1*o0*u0*w1 + 2*j1*k0*o1*u0*w1 + 
  2*j0*k1*o1*u0*w1 - 2*i1*o0*o1*u0*w1 - 
  i0*o1*o1*u0*w1 - j1*j1*r0*u0*w1 + i1*n1*r0*u0*w1 - 
  2*j0*j1*r1*u0*w1 + i1*n0*r1*u0*w1 + 
  i0*n1*r1*u0*w1 - 2*k0*k1*n0*u1*w1 - k0*k0*n1*u1*w1 + 
  2*j1*k0*o0*u1*w1 + 2*j0*k1*o0*u1*w1 - 
  i1*o0*o0*u1*w1 + 2*j0*k0*o1*u1*w1 - 
  2*i0*o0*o1*u1*w1 - 2*j0*j1*r0*u1*w1 + 
  i1*n0*r0*u1*w1 + i0*n1*r0*u1*w1 - j0*j0*r1*u1*w1 + 
  i0*n0*r1*u1*w1;
   
    k[19] = m1*m1*p1*p1*r0 - 2*l1*m1*p1*q1*r0 + 
  l1*l1*q1*q1*r0 + 2*m1*m1*p0*p1*r1 + 2*m0*m1*p1*p1*r1 - 
  2*l1*m1*p1*q0*r1 - 2*l1*m1*p0*q1*r1 - 2*l1*m0*p1*q1*r1 - 
  2*l0*m1*p1*q1*r1 + 2*l1*l1*q0*q1*r1 + 2*l0*l1*q1*q1*r1 - 
  2*m1*m1*o1*p1*s0 + 2*l1*m1*o1*q1*s0 + 2*k1*m1*p1*q1*s0 - 
  2*k1*l1*q1*q1*s0 - 2*m1*m1*o1*p0*s1 - 2*m1*m1*o0*p1*s1 - 
  4*m0*m1*o1*p1*s1 + 2*l1*m1*o1*q0*s1 + 2*k1*m1*p1*q0*s1 + 
  2*l1*m1*o0*q1*s1 + 2*l1*m0*o1*q1*s1 + 2*l0*m1*o1*q1*s1 + 
  2*k1*m1*p0*q1*s1 + 2*k1*m0*p1*q1*s1 + 2*k0*m1*p1*q1*s1 - 
  4*k1*l1*q0*q1*s1 - 2*k1*l0*q1*q1*s1 - 2*k0*l1*q1*q1*s1 + 
  2*m1*m1*n1*s0*s1 - 4*j1*m1*q1*s0*s1 + 2*i1*q1*q1*s0*s1 + 
  m1*m1*n0*s1*s1 + 2*m0*m1*n1*s1*s1 - 2*j1*m1*q0*s1*s1 - 
  2*j1*m0*q1*s1*s1 - 2*j0*m1*q1*s1*s1 + 2*i1*q0*q1*s1*s1 + 
  i0*q1*q1*s1*s1 + 2*l1*m1*o1*p1*t0 - 2*k1*m1*p1*p1*t0 - 
  2*l1*l1*o1*q1*t0 + 2*k1*l1*p1*q1*t0 - 2*l1*m1*n1*s1*t0 + 
  2*j1*m1*p1*s1*t0 + 2*j1*l1*q1*s1*t0 - 2*i1*p1*q1*s1*t0 + 
  2*l1*m1*o1*p0*t1 + 2*l1*m1*o0*p1*t1 + 2*l1*m0*o1*p1*t1 + 
  2*l0*m1*o1*p1*t1 - 4*k1*m1*p0*p1*t1 - 2*k1*m0*p1*p1*t1 - 
  2*k0*m1*p1*p1*t1 - 2*l1*l1*o1*q0*t1 + 2*k1*l1*p1*q0*t1 - 
  2*l1*l1*o0*q1*t1 - 4*l0*l1*o1*q1*t1 + 2*k1*l1*p0*q1*t1 + 
  2*k1*l0*p1*q1*t1 + 2*k0*l1*p1*q1*t1 - 2*l1*m1*n1*s0*t1 + 
  2*j1*m1*p1*s0*t1 + 2*j1*l1*q1*s0*t1 - 2*i1*p1*q1*s0*t1 - 
  2*l1*m1*n0*s1*t1 - 2*l1*m0*n1*s1*t1 - 2*l0*m1*n1*s1*t1 + 
  2*j1*m1*p0*s1*t1 + 2*j1*m0*p1*s1*t1 + 2*j0*m1*p1*s1*t1 + 
  2*j1*l1*q0*s1*t1 - 2*i1*p1*q0*s1*t1 + 2*j1*l0*q1*s1*t1 + 
  2*j0*l1*q1*s1*t1 - 2*i1*p0*q1*s1*t1 - 2*i0*p1*q1*s1*t1 + 
  2*l1*l1*n1*t0*t1 - 4*j1*l1*p1*t0*t1 + 2*i1*p1*p1*t0*t1 + 
  l1*l1*n0*t1*t1 + 2*l0*l1*n1*t1*t1 - 2*j1*l1*p0*t1*t1 - 
  2*j1*l0*p1*t1*t1 - 2*j0*l1*p1*t1*t1 + 2*i1*p0*p1*t1*t1 + 
  i0*p1*p1*t1*t1 + m1*m1*o1*o1*u0 - 2*k1*m1*o1*q1*u0 + 
  k1*k1*q1*q1*u0 - m1*m1*n1*r1*u0 + 2*j1*m1*q1*r1*u0 - 
  i1*q1*q1*r1*u0 + 2*k1*m1*n1*t1*u0 - 2*j1*m1*o1*t1*u0 - 
  2*j1*k1*q1*t1*u0 + 2*i1*o1*q1*t1*u0 + j1*j1*t1*t1*u0 - 
  i1*n1*t1*t1*u0 + 2*m1*m1*o0*o1*u1 + 2*m0*m1*o1*o1*u1 - 
  2*k1*m1*o1*q0*u1 - 2*k1*m1*o0*q1*u1 - 2*k1*m0*o1*q1*u1 - 
  2*k0*m1*o1*q1*u1 + 2*k1*k1*q0*q1*u1 + 2*k0*k1*q1*q1*u1 - 
  m1*m1*n1*r0*u1 + 2*j1*m1*q1*r0*u1 - i1*q1*q1*r0*u1 - 
  m1*m1*n0*r1*u1 - 2*m0*m1*n1*r1*u1 + 2*j1*m1*q0*r1*u1 + 
  2*j1*m0*q1*r1*u1 + 2*j0*m1*q1*r1*u1 - 2*i1*q0*q1*r1*u1 - 
  i0*q1*q1*r1*u1 + 2*k1*m1*n1*t0*u1 - 2*j1*m1*o1*t0*u1 - 
  2*j1*k1*q1*t0*u1 + 2*i1*o1*q1*t0*u1 + 2*k1*m1*n0*t1*u1 + 
  2*k1*m0*n1*t1*u1 + 2*k0*m1*n1*t1*u1 - 2*j1*m1*o0*t1*u1 - 
  2*j1*m0*o1*t1*u1 - 2*j0*m1*o1*t1*u1 - 2*j1*k1*q0*t1*u1 + 
  2*i1*o1*q0*t1*u1 - 2*j1*k0*q1*t1*u1 - 2*j0*k1*q1*t1*u1 + 
  2*i1*o0*q1*t1*u1 + 2*i0*o1*q1*t1*u1 + 2*j1*j1*t0*t1*u1 - 
  2*i1*n1*t0*t1*u1 + 2*j0*j1*t1*t1*u1 - i1*n0*t1*t1*u1 - 
  i0*n1*t1*t1*u1 - 2*l1*m1*o1*o1*v0 + 2*k1*m1*o1*p1*v0 + 
  2*k1*l1*o1*q1*v0 - 2*k1*k1*p1*q1*v0 + 2*l1*m1*n1*r1*v0 - 
  2*j1*m1*p1*r1*v0 - 2*j1*l1*q1*r1*v0 + 2*i1*p1*q1*r1*v0 - 
  2*k1*m1*n1*s1*v0 + 2*j1*m1*o1*s1*v0 + 2*j1*k1*q1*s1*v0 - 
  2*i1*o1*q1*s1*v0 - 2*k1*l1*n1*t1*v0 + 2*j1*l1*o1*t1*v0 + 
  2*j1*k1*p1*t1*v0 - 2*i1*o1*p1*t1*v0 - 2*j1*j1*s1*t1*v0 + 
  2*i1*n1*s1*t1*v0 - 4*l1*m1*o0*o1*v1 - 2*l1*m0*o1*o1*v1 - 
  2*l0*m1*o1*o1*v1 + 2*k1*m1*o1*p0*v1 + 2*k1*m1*o0*p1*v1 + 
  2*k1*m0*o1*p1*v1 + 2*k0*m1*o1*p1*v1 + 2*k1*l1*o1*q0*v1 - 
  2*k1*k1*p1*q0*v1 + 2*k1*l1*o0*q1*v1 + 2*k1*l0*o1*q1*v1 + 
  2*k0*l1*o1*q1*v1 - 2*k1*k1*p0*q1*v1 - 4*k0*k1*p1*q1*v1 + 
  2*l1*m1*n1*r0*v1 - 2*j1*m1*p1*r0*v1 - 2*j1*l1*q1*r0*v1 + 
  2*i1*p1*q1*r0*v1 + 2*l1*m1*n0*r1*v1 + 2*l1*m0*n1*r1*v1 + 
  2*l0*m1*n1*r1*v1 - 2*j1*m1*p0*r1*v1 - 2*j1*m0*p1*r1*v1 - 
  2*j0*m1*p1*r1*v1 - 2*j1*l1*q0*r1*v1 + 2*i1*p1*q0*r1*v1 - 
  2*j1*l0*q1*r1*v1 - 2*j0*l1*q1*r1*v1 + 2*i1*p0*q1*r1*v1 + 
  2*i0*p1*q1*r1*v1 - 2*k1*m1*n1*s0*v1 + 2*j1*m1*o1*s0*v1 + 
  2*j1*k1*q1*s0*v1 - 2*i1*o1*q1*s0*v1 - 2*k1*m1*n0*s1*v1 - 
  2*k1*m0*n1*s1*v1 - 2*k0*m1*n1*s1*v1 + 2*j1*m1*o0*s1*v1 + 
  2*j1*m0*o1*s1*v1 + 2*j0*m1*o1*s1*v1 + 2*j1*k1*q0*s1*v1 - 
  2*i1*o1*q0*s1*v1 + 2*j1*k0*q1*s1*v1 + 2*j0*k1*q1*s1*v1 - 
  2*i1*o0*q1*s1*v1 - 2*i0*o1*q1*s1*v1 - 2*k1*l1*n1*t0*v1 + 
  2*j1*l1*o1*t0*v1 + 2*j1*k1*p1*t0*v1 - 2*i1*o1*p1*t0*v1 - 
  2*j1*j1*s1*t0*v1 + 2*i1*n1*s1*t0*v1 - 2*k1*l1*n0*t1*v1 - 
  2*k1*l0*n1*t1*v1 - 2*k0*l1*n1*t1*v1 + 2*j1*l1*o0*t1*v1 + 
  2*j1*l0*o1*t1*v1 + 2*j0*l1*o1*t1*v1 + 2*j1*k1*p0*t1*v1 - 
  2*i1*o1*p0*t1*v1 + 2*j1*k0*p1*t1*v1 + 2*j0*k1*p1*t1*v1 - 
  2*i1*o0*p1*t1*v1 - 2*i0*o1*p1*t1*v1 - 2*j1*j1*s0*t1*v1 + 
  2*i1*n1*s0*t1*v1 - 4*j0*j1*s1*t1*v1 + 2*i1*n0*s1*t1*v1 + 
  2*i0*n1*s1*t1*v1 + 2*k1*k1*n1*v0*v1 - 4*j1*k1*o1*v0*v1 + 
  2*i1*o1*o1*v0*v1 + 2*j1*j1*r1*v0*v1 - 2*i1*n1*r1*v0*v1 + 
  k1*k1*n0*v1*v1 + 2*k0*k1*n1*v1*v1 - 2*j1*k1*o0*v1*v1 - 
  2*j1*k0*o1*v1*v1 - 2*j0*k1*o1*v1*v1 + 2*i1*o0*o1*v1*v1 + 
  i0*o1*o1*v1*v1 + j1*j1*r0*v1*v1 - i1*n1*r0*v1*v1 + 
  2*j0*j1*r1*v1*v1 - i1*n0*r1*v1*v1 - i0*n1*r1*v1*v1 + 
  l1*l1*o1*o1*w0 - 2*k1*l1*o1*p1*w0 + k1*k1*p1*p1*w0 - 
  l1*l1*n1*r1*w0 + 2*j1*l1*p1*r1*w0 - i1*p1*p1*r1*w0 + 
  2*k1*l1*n1*s1*w0 - 2*j1*l1*o1*s1*w0 - 2*j1*k1*p1*s1*w0 + 
  2*i1*o1*p1*s1*w0 + j1*j1*s1*s1*w0 - i1*n1*s1*s1*w0 - 
  k1*k1*n1*u1*w0 + 2*j1*k1*o1*u1*w0 - i1*o1*o1*u1*w0 - 
  j1*j1*r1*u1*w0 + i1*n1*r1*u1*w0 + 2*l1*l1*o0*o1*w1 + 
  2*l0*l1*o1*o1*w1 - 2*k1*l1*o1*p0*w1 - 2*k1*l1*o0*p1*w1 - 
  2*k1*l0*o1*p1*w1 - 2*k0*l1*o1*p1*w1 + 2*k1*k1*p0*p1*w1 + 
  2*k0*k1*p1*p1*w1 - l1*l1*n1*r0*w1 + 2*j1*l1*p1*r0*w1 - 
  i1*p1*p1*r0*w1 - l1*l1*n0*r1*w1 - 2*l0*l1*n1*r1*w1 + 
  2*j1*l1*p0*r1*w1 + 2*j1*l0*p1*r1*w1 + 2*j0*l1*p1*r1*w1 - 
  2*i1*p0*p1*r1*w1 - i0*p1*p1*r1*w1 + 2*k1*l1*n1*s0*w1 - 
  2*j1*l1*o1*s0*w1 - 2*j1*k1*p1*s0*w1 + 2*i1*o1*p1*s0*w1 + 
  2*k1*l1*n0*s1*w1 + 2*k1*l0*n1*s1*w1 + 2*k0*l1*n1*s1*w1 - 
  2*j1*l1*o0*s1*w1 - 2*j1*l0*o1*s1*w1 - 2*j0*l1*o1*s1*w1 - 
  2*j1*k1*p0*s1*w1 + 2*i1*o1*p0*s1*w1 - 2*j1*k0*p1*s1*w1 - 
  2*j0*k1*p1*s1*w1 + 2*i1*o0*p1*s1*w1 + 2*i0*o1*p1*s1*w1 + 
  2*j1*j1*s0*s1*w1 - 2*i1*n1*s0*s1*w1 + 2*j0*j1*s1*s1*w1 - 
  i1*n0*s1*s1*w1 - i0*n1*s1*s1*w1 - k1*k1*n1*u0*w1 + 
  2*j1*k1*o1*u0*w1 - i1*o1*o1*u0*w1 - j1*j1*r1*u0*w1 + 
  i1*n1*r1*u0*w1 - k1*k1*n0*u1*w1 - 2*k0*k1*n1*u1*w1 + 
  2*j1*k1*o0*u1*w1 + 2*j1*k0*o1*u1*w1 + 2*j0*k1*o1*u1*w1 - 
  2*i1*o0*o1*u1*w1 - i0*o1*o1*u1*w1 - j1*j1*r0*u1*w1 + 
  i1*n1*r0*u1*w1 - 2*j0*j1*r1*u1*w1 + i1*n0*r1*u1*w1 + 
  i0*n1*r1*u1*w1;
   
    k[20] = m1*m1*p1*p1*r1 - 2*l1*m1*p1*q1*r1 + 
  l1*l1*q1*q1*r1 - 2*m1*m1*o1*p1*s1 + 2*l1*m1*o1*q1*s1 + 
  2*k1*m1*p1*q1*s1 - 2*k1*l1*q1*q1*s1 + m1*m1*n1*s1*s1 - 
  2*j1*m1*q1*s1*s1 + i1*q1*q1*s1*s1 + 2*l1*m1*o1*p1*t1 - 
  2*k1*m1*p1*p1*t1 - 2*l1*l1*o1*q1*t1 + 2*k1*l1*p1*q1*t1 - 
  2*l1*m1*n1*s1*t1 + 2*j1*m1*p1*s1*t1 + 2*j1*l1*q1*s1*t1 - 
  2*i1*p1*q1*s1*t1 + l1*l1*n1*t1*t1 - 2*j1*l1*p1*t1*t1 + 
  i1*p1*p1*t1*t1 + m1*m1*o1*o1*u1 - 2*k1*m1*o1*q1*u1 + 
  k1*k1*q1*q1*u1 - m1*m1*n1*r1*u1 + 2*j1*m1*q1*r1*u1 - 
  i1*q1*q1*r1*u1 + 2*k1*m1*n1*t1*u1 - 2*j1*m1*o1*t1*u1 - 
  2*j1*k1*q1*t1*u1 + 2*i1*o1*q1*t1*u1 + j1*j1*t1*t1*u1 - 
  i1*n1*t1*t1*u1 - 2*l1*m1*o1*o1*v1 + 2*k1*m1*o1*p1*v1 + 
  2*k1*l1*o1*q1*v1 - 2*k1*k1*p1*q1*v1 + 2*l1*m1*n1*r1*v1 - 
  2*j1*m1*p1*r1*v1 - 2*j1*l1*q1*r1*v1 + 2*i1*p1*q1*r1*v1 - 
  2*k1*m1*n1*s1*v1 + 2*j1*m1*o1*s1*v1 + 2*j1*k1*q1*s1*v1 - 
  2*i1*o1*q1*s1*v1 - 2*k1*l1*n1*t1*v1 + 2*j1*l1*o1*t1*v1 + 
  2*j1*k1*p1*t1*v1 - 2*i1*o1*p1*t1*v1 - 2*j1*j1*s1*t1*v1 + 
  2*i1*n1*s1*t1*v1 + k1*k1*n1*v1*v1 - 2*j1*k1*o1*v1*v1 + 
  i1*o1*o1*v1*v1 + j1*j1*r1*v1*v1 - i1*n1*r1*v1*v1 + 
  l1*l1*o1*o1*w1 - 2*k1*l1*o1*p1*w1 + k1*k1*p1*p1*w1 - 
  l1*l1*n1*r1*w1 + 2*j1*l1*p1*r1*w1 - i1*p1*p1*r1*w1 + 
  2*k1*l1*n1*s1*w1 - 2*j1*l1*o1*s1*w1 - 2*j1*k1*p1*s1*w1 + 
  2*i1*o1*p1*s1*w1 + j1*j1*s1*s1*w1 - i1*n1*s1*s1*w1 - 
  k1*k1*n1*u1*w1 + 2*j1*k1*o1*u1*w1 - i1*o1*o1*u1*w1 - 
  j1*j1*r1*u1*w1 + i1*n1*r1*u1*w1;
   
	return 	rg_ImplicitEquation(5, k);
}

rg_REAL rg_IntsecImplicit5::inversion(const rg_BzCurve2D &curve, const rg_Point2D &point)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL x3 = curve.getCtrlPt(3).getX();
	rg_REAL x4 = curve.getCtrlPt(4).getX();
    rg_REAL x5 = curve.getCtrlPt(5).getX();

	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL y3 = curve.getCtrlPt(3).getY();
	rg_REAL y4 = curve.getCtrlPt(4).getY();
    rg_REAL y5 = curve.getCtrlPt(5).getY();

//  If
//		    
//	x(t) = a5 t^5 + a4 t^4 + a3 t^3 + a2 t^2 + a1 t + a0
//
//	y(t) = b5 t^5 + b4 t^4 + b3 t^3 + b2 t^2 + b1 t + b0
//  
//  then the resultant of p(x, t) and q(y, t) is like the following.

    rg_REAL a5 = x5- 5*x4 + 10*x3 - 10*x2 +  5*x1 -    x0;
	rg_REAL a4 =     5*x4 - 20*x3 + 30*x2 - 20*x1 +  5*x0;
	rg_REAL a3 =            10*x3 - 30*x2 + 30*x1 - 10*x0;
	rg_REAL a2 =                    10*x2 - 20*x1 + 10*x0;
	rg_REAL a1 =                             5*x1 -  5*x0;
	rg_REAL a0 =                                       x0;

    rg_REAL b5 = y5- 5*y4 + 10*y3 - 10*y2 +  5*y1 -    y0;
	rg_REAL b4 =     5*y4 - 20*y3 + 30*y2 - 20*y1 +  5*y0;
	rg_REAL b3 =            10*y3 - 30*y2 + 30*y1 - 10*y0;
	rg_REAL b2 =                    10*y2 - 20*y1 + 10*y0;
	rg_REAL b1 =                             5*y1 -  5*y0;
	rg_REAL b0 =                                       y0;

  //			|														                                      |
  //			|	i0 x + i1 y + i2  j0 x + j1 y + j2	k0 x + k1 y + k2  l0 x + l1 y + l2  m0 x + m1 y + m2  |
  //			|														                                      |
  // R(P, Q) =	|	j0 x + j1 y + j2  n0 x + n1 y + n2	o0 x + o1 y + o2  p0 x + p1 y + p2  q0 x + q1 y + q2  |
  //			|														                                      |
  //			|	k0 x + k1 y + k2  o0 x + o1 y + o2	r0 x + r1 y + r2  s0 x + s1 y + s2  t0 x + t1 y + t2  |
  //			|														                                      |
  //            |   l0 x + l1 y + l2  p0 x + p1 y + p2  s0 x + s1 y + s2  u0 x + u1 y + u2  v0 x + v1 y + v2  |
  //            |                                                                                             |
  //            |   m0 x + m1 y + m2  q0 x + q1 y + q2  t0 x + t1 y + t2  v0 x + v1 y + v2  w0 x + w1 y + w2  |
   
	rg_REAL i0 =   0.0;
	rg_REAL i1 =   0.0;
	rg_REAL i2 =   a5*b4 - a4*b5;

	rg_REAL j0 =   0.0;
	rg_REAL j1 =   0.0;
	rg_REAL j2 =   a5*b3 - a3*b5;

	rg_REAL k0 =   0.0;
	rg_REAL k1 =   0.0;
	rg_REAL k2 =   a5*b2 - a2*b5;

	rg_REAL l0 =   0.0;
	rg_REAL l1 =   0.0;
	rg_REAL l2 =   a5*b1 - a1*b5;

	rg_REAL m0 =   b5;
	rg_REAL m1 = - a5;
	rg_REAL m2 =   a5*b0 - a0*b5;

	rg_REAL n0 =   0.0;
	rg_REAL n1 =   0.0;
	rg_REAL n2 =   a5*b2 - a2*b5 + a4*b3 - a3*b4;

	rg_REAL o0 =   0.0;
	rg_REAL o1 =   0.0;
	rg_REAL o2 =   a5*b1 - a1*b5 + a4*b2 - a2*b4;

    rg_REAL p0 =   b5;
	rg_REAL p1 = - a5;
	rg_REAL p2 =   a5*b0 - a0*b5 + a4*b1 - a1*b4;

  	rg_REAL q0 =   b4;
	rg_REAL q1 = - a4;
	rg_REAL q2 =   a4*b0 - a0*b4;

  	rg_REAL r0 =   b5;
	rg_REAL r1 = - a5;
	rg_REAL r2 =   a5*b0 - a0*b5 + a4*b1 - a1*b4 + a3*b2 - a2*b3;

  	rg_REAL s0 =   b4;
	rg_REAL s1 = - a4;
	rg_REAL s2 =   a4*b0 - a0*b4 + a3*b1 - a1*b3;

   	rg_REAL t0 =   b3;
	rg_REAL t1 = - a3;
	rg_REAL t2 =   a3*b0 - a0*b3;

   	rg_REAL u0 =   b3;
	rg_REAL u1 = - a3;
	rg_REAL u2 =   a3*b0 - a0*b3 + a2*b1 - a1*b2;

	rg_REAL v0 =   b2;
	rg_REAL v1 = - a2;
	rg_REAL v2 =   a2*b0 - a0*b2;

   	rg_REAL w0 =   b1;
	rg_REAL w1 = - a1;
	rg_REAL w2 =   a1*b0 - a0*b1;

// 고쳐야할 부분
	rg_REAL x  =  point.getX();
	rg_REAL y  =  point.getY();

	rg_REAL t_4 =m2*p2*p2*r2 - l2*p2*q2*r2 - 2*m2*o2*p2*s2 + l2*o2*q2*s2 + k2*p2*q2*s2 + 
  m2*n2*s2*s2 - j2*q2*s2*s2 + l2*o2*p2*t2 - k2*p2*p2*t2 - l2*n2*s2*t2 + 
  j2*p2*s2*t2 + m2*o2*o2*u2 - k2*o2*q2*u2 - m2*n2*r2*u2 + j2*q2*r2*u2 + 
  k2*n2*t2*u2 - j2*o2*t2*u2 - l2*o2*o2*v2 + k2*o2*p2*v2 + l2*n2*r2*v2 - 
  j2*p2*r2*v2 - k2*n2*s2*v2 + j2*o2*s2*v2 + m2*p2*p2*r0*x - l2*p2*q2*r0*x + 
  2*m2*p0*p2*r2*x + m0*p2*p2*r2*x - l2*p2*q0*r2*x - l2*p0*q2*r2*x - 
  l0*p2*q2*r2*x - 2*m2*o2*p2*s0*x + l2*o2*q2*s0*x + k2*p2*q2*s0*x - 
  2*m2*o2*p0*s2*x - 2*m2*o0*p2*s2*x - 2*m0*o2*p2*s2*x + l2*o2*q0*s2*x + 
  k2*p2*q0*s2*x + l2*o0*q2*s2*x + l0*o2*q2*s2*x + k2*p0*q2*s2*x + 
  k0*p2*q2*s2*x + 2*m2*n2*s0*s2*x - 2*j2*q2*s0*s2*x + m2*n0*s2*s2*x + 
  m0*n2*s2*s2*x - j2*q0*s2*s2*x - j0*q2*s2*s2*x + l2*o2*p2*t0*x - k2*p2*p2*t0*x - 
  l2*n2*s2*t0*x + j2*p2*s2*t0*x + l2*o2*p0*t2*x + l2*o0*p2*t2*x + 
  l0*o2*p2*t2*x - 2*k2*p0*p2*t2*x - k0*p2*p2*t2*x - l2*n2*s0*t2*x + 
  j2*p2*s0*t2*x - l2*n0*s2*t2*x - l0*n2*s2*t2*x + j2*p0*s2*t2*x + 
  j0*p2*s2*t2*x + m2*o2*o2*u0*x - k2*o2*q2*u0*x - m2*n2*r2*u0*x + 
  j2*q2*r2*u0*x + k2*n2*t2*u0*x - j2*o2*t2*u0*x + 2*m2*o0*o2*u2*x + 
  m0*o2*o2*u2*x - k2*o2*q0*u2*x - k2*o0*q2*u2*x - k0*o2*q2*u2*x - 
  m2*n2*r0*u2*x + j2*q2*r0*u2*x - m2*n0*r2*u2*x - m0*n2*r2*u2*x + 
  j2*q0*r2*u2*x + j0*q2*r2*u2*x + k2*n2*t0*u2*x - j2*o2*t0*u2*x + 
  k2*n0*t2*u2*x + k0*n2*t2*u2*x - j2*o0*t2*u2*x - j0*o2*t2*u2*x - 
  l2*o2*o2*v0*x + k2*o2*p2*v0*x + l2*n2*r2*v0*x - j2*p2*r2*v0*x - 
  k2*n2*s2*v0*x + j2*o2*s2*v0*x - 2*l2*o0*o2*v2*x - l0*o2*o2*v2*x + 
  k2*o2*p0*v2*x + k2*o0*p2*v2*x + k0*o2*p2*v2*x + l2*n2*r0*v2*x - 
  j2*p2*r0*v2*x + l2*n0*r2*v2*x + l0*n2*r2*v2*x - j2*p0*r2*v2*x - 
  j0*p2*r2*v2*x - k2*n2*s0*v2*x + j2*o2*s0*v2*x - k2*n0*s2*v2*x - 
  k0*n2*s2*v2*x + j2*o0*s2*v2*x + j0*o2*s2*v2*x + 2*m2*p0*p2*r0*x*x + 
  m0*p2*p2*r0*x*x - l2*p2*q0*r0*x*x - l2*p0*q2*r0*x*x - l0*p2*q2*r0*x*x + 
  m2*p0*p0*r2*x*x + 2*m0*p0*p2*r2*x*x - l2*p0*q0*r2*x*x - l0*p2*q0*r2*x*x - 
  l0*p0*q2*r2*x*x - 2*m2*o2*p0*s0*x*x - 2*m2*o0*p2*s0*x*x - 
  2*m0*o2*p2*s0*x*x + l2*o2*q0*s0*x*x + k2*p2*q0*s0*x*x + l2*o0*q2*s0*x*x + 
  l0*o2*q2*s0*x*x + k2*p0*q2*s0*x*x + k0*p2*q2*s0*x*x + m2*n2*s0*s0*x*x - 
  j2*q2*s0*s0*x*x - 2*m2*o0*p0*s2*x*x - 2*m0*o2*p0*s2*x*x - 
  2*m0*o0*p2*s2*x*x + l2*o0*q0*s2*x*x + l0*o2*q0*s2*x*x + k2*p0*q0*s2*x*x + 
  k0*p2*q0*s2*x*x + l0*o0*q2*s2*x*x + k0*p0*q2*s2*x*x + 2*m2*n0*s0*s2*x*x + 
  2*m0*n2*s0*s2*x*x - 2*j2*q0*s0*s2*x*x - 2*j0*q2*s0*s2*x*x + 
  m0*n0*s2*s2*x*x - j0*q0*s2*s2*x*x + l2*o2*p0*t0*x*x + l2*o0*p2*t0*x*x + 
  l0*o2*p2*t0*x*x - 2*k2*p0*p2*t0*x*x - k0*p2*p2*t0*x*x - l2*n2*s0*t0*x*x + 
  j2*p2*s0*t0*x*x - l2*n0*s2*t0*x*x - l0*n2*s2*t0*x*x + j2*p0*s2*t0*x*x + 
  j0*p2*s2*t0*x*x + l2*o0*p0*t2*x*x + l0*o2*p0*t2*x*x - k2*p0*p0*t2*x*x + 
  l0*o0*p2*t2*x*x - 2*k0*p0*p2*t2*x*x - l2*n0*s0*t2*x*x - l0*n2*s0*t2*x*x + 
  j2*p0*s0*t2*x*x + j0*p2*s0*t2*x*x - l0*n0*s2*t2*x*x + j0*p0*s2*t2*x*x + 
  2*m2*o0*o2*u0*x*x + m0*o2*o2*u0*x*x - k2*o2*q0*u0*x*x - k2*o0*q2*u0*x*x - 
  k0*o2*q2*u0*x*x - m2*n2*r0*u0*x*x + j2*q2*r0*u0*x*x - m2*n0*r2*u0*x*x - 
  m0*n2*r2*u0*x*x + j2*q0*r2*u0*x*x + j0*q2*r2*u0*x*x + k2*n2*t0*u0*x*x - 
  j2*o2*t0*u0*x*x + k2*n0*t2*u0*x*x + k0*n2*t2*u0*x*x - j2*o0*t2*u0*x*x - 
  j0*o2*t2*u0*x*x + m2*o0*o0*u2*x*x + 2*m0*o0*o2*u2*x*x - k2*o0*q0*u2*x*x - 
  k0*o2*q0*u2*x*x - k0*o0*q2*u2*x*x - m2*n0*r0*u2*x*x - m0*n2*r0*u2*x*x + 
  j2*q0*r0*u2*x*x + j0*q2*r0*u2*x*x - m0*n0*r2*u2*x*x + j0*q0*r2*u2*x*x + 
  k2*n0*t0*u2*x*x + k0*n2*t0*u2*x*x - j2*o0*t0*u2*x*x - j0*o2*t0*u2*x*x + 
  k0*n0*t2*u2*x*x - j0*o0*t2*u2*x*x - 2*l2*o0*o2*v0*x*x - l0*o2*o2*v0*x*x + 
  k2*o2*p0*v0*x*x + k2*o0*p2*v0*x*x + k0*o2*p2*v0*x*x + l2*n2*r0*v0*x*x - 
  j2*p2*r0*v0*x*x + l2*n0*r2*v0*x*x + l0*n2*r2*v0*x*x - j2*p0*r2*v0*x*x - 
  j0*p2*r2*v0*x*x - k2*n2*s0*v0*x*x + j2*o2*s0*v0*x*x - k2*n0*s2*v0*x*x - 
  k0*n2*s2*v0*x*x + j2*o0*s2*v0*x*x + j0*o2*s2*v0*x*x - l2*o0*o0*v2*x*x - 
  2*l0*o0*o2*v2*x*x + k2*o0*p0*v2*x*x + k0*o2*p0*v2*x*x + k0*o0*p2*v2*x*x + 
  l2*n0*r0*v2*x*x + l0*n2*r0*v2*x*x - j2*p0*r0*v2*x*x - j0*p2*r0*v2*x*x + 
  l0*n0*r2*v2*x*x - j0*p0*r2*v2*x*x - k2*n0*s0*v2*x*x - k0*n2*s0*v2*x*x + 
  j2*o0*s0*v2*x*x + j0*o2*s0*v2*x*x - k0*n0*s2*v2*x*x + j0*o0*s2*v2*x*x + 
  m2*p0*p0*r0*x*x*x + 2*m0*p0*p2*r0*x*x*x - l2*p0*q0*r0*x*x*x - l0*p2*q0*r0*x*x*x - 
  l0*p0*q2*r0*x*x*x + m0*p0*p0*r2*x*x*x - l0*p0*q0*r2*x*x*x - 2*m2*o0*p0*s0*x*x*x - 
  2*m0*o2*p0*s0*x*x*x - 2*m0*o0*p2*s0*x*x*x + l2*o0*q0*s0*x*x*x + l0*o2*q0*s0*x*x*x + 
  k2*p0*q0*s0*x*x*x + k0*p2*q0*s0*x*x*x + l0*o0*q2*s0*x*x*x + k0*p0*q2*s0*x*x*x + 
  m2*n0*s0*s0*x*x*x + m0*n2*s0*s0*x*x*x - j2*q0*s0*s0*x*x*x - j0*q2*s0*s0*x*x*x - 
  2*m0*o0*p0*s2*x*x*x + l0*o0*q0*s2*x*x*x + k0*p0*q0*s2*x*x*x + 2*m0*n0*s0*s2*x*x*x - 
  2*j0*q0*s0*s2*x*x*x + l2*o0*p0*t0*x*x*x + l0*o2*p0*t0*x*x*x - k2*p0*p0*t0*x*x*x + 
  l0*o0*p2*t0*x*x*x - 2*k0*p0*p2*t0*x*x*x - l2*n0*s0*t0*x*x*x - l0*n2*s0*t0*x*x*x + 
  j2*p0*s0*t0*x*x*x + j0*p2*s0*t0*x*x*x - l0*n0*s2*t0*x*x*x + j0*p0*s2*t0*x*x*x + 
  l0*o0*p0*t2*x*x*x - k0*p0*p0*t2*x*x*x - l0*n0*s0*t2*x*x*x + j0*p0*s0*t2*x*x*x + 
  m2*o0*o0*u0*x*x*x + 2*m0*o0*o2*u0*x*x*x - k2*o0*q0*u0*x*x*x - k0*o2*q0*u0*x*x*x - 
  k0*o0*q2*u0*x*x*x - m2*n0*r0*u0*x*x*x - m0*n2*r0*u0*x*x*x + j2*q0*r0*u0*x*x*x + 
  j0*q2*r0*u0*x*x*x - m0*n0*r2*u0*x*x*x + j0*q0*r2*u0*x*x*x + k2*n0*t0*u0*x*x*x + 
  k0*n2*t0*u0*x*x*x - j2*o0*t0*u0*x*x*x - j0*o2*t0*u0*x*x*x + k0*n0*t2*u0*x*x*x - 
  j0*o0*t2*u0*x*x*x + m0*o0*o0*u2*x*x*x - k0*o0*q0*u2*x*x*x - m0*n0*r0*u2*x*x*x + 
  j0*q0*r0*u2*x*x*x + k0*n0*t0*u2*x*x*x - j0*o0*t0*u2*x*x*x - l2*o0*o0*v0*x*x*x - 
  2*l0*o0*o2*v0*x*x*x + k2*o0*p0*v0*x*x*x + k0*o2*p0*v0*x*x*x + k0*o0*p2*v0*x*x*x + 
  l2*n0*r0*v0*x*x*x + l0*n2*r0*v0*x*x*x - j2*p0*r0*v0*x*x*x - j0*p2*r0*v0*x*x*x + 
  l0*n0*r2*v0*x*x*x - j0*p0*r2*v0*x*x*x - k2*n0*s0*v0*x*x*x - k0*n2*s0*v0*x*x*x + 
  j2*o0*s0*v0*x*x*x + j0*o2*s0*v0*x*x*x - k0*n0*s2*v0*x*x*x + j0*o0*s2*v0*x*x*x - 
  l0*o0*o0*v2*x*x*x + k0*o0*p0*v2*x*x*x + l0*n0*r0*v2*x*x*x - j0*p0*r0*v2*x*x*x - 
  k0*n0*s0*v2*x*x*x + j0*o0*s0*v2*x*x*x + m0*p0*p0*r0*x*x*x*x - l0*p0*q0*r0*x*x*x*x - 
  2*m0*o0*p0*s0*x*x*x*x + l0*o0*q0*s0*x*x*x*x + k0*p0*q0*s0*x*x*x*x + m0*n0*s0*s0*x*x*x*x - 
  j0*q0*s0*s0*x*x*x*x + l0*o0*p0*t0*x*x*x*x - k0*p0*p0*t0*x*x*x*x - l0*n0*s0*t0*x*x*x*x + 
  j0*p0*s0*t0*x*x*x*x + m0*o0*o0*u0*x*x*x*x - k0*o0*q0*u0*x*x*x*x - m0*n0*r0*u0*x*x*x*x + 
  j0*q0*r0*u0*x*x*x*x + k0*n0*t0*u0*x*x*x*x - j0*o0*t0*u0*x*x*x*x - l0*o0*o0*v0*x*x*x*x + 
  k0*o0*p0*v0*x*x*x*x + l0*n0*r0*v0*x*x*x*x - j0*p0*r0*v0*x*x*x*x - k0*n0*s0*v0*x*x*x*x + 
  j0*o0*s0*v0*x*x*x*x + m2*p2*p2*r1*y - l2*p2*q2*r1*y + 2*m2*p1*p2*r2*y + 
  m1*p2*p2*r2*y - l2*p2*q1*r2*y - l2*p1*q2*r2*y - l1*p2*q2*r2*y - 
  2*m2*o2*p2*s1*y + l2*o2*q2*s1*y + k2*p2*q2*s1*y - 2*m2*o2*p1*s2*y - 
  2*m2*o1*p2*s2*y - 2*m1*o2*p2*s2*y + l2*o2*q1*s2*y + k2*p2*q1*s2*y + 
  l2*o1*q2*s2*y + l1*o2*q2*s2*y + k2*p1*q2*s2*y + k1*p2*q2*s2*y + 
  2*m2*n2*s1*s2*y - 2*j2*q2*s1*s2*y + m2*n1*s2*s2*y + m1*n2*s2*s2*y - 
  j2*q1*s2*s2*y - j1*q2*s2*s2*y + l2*o2*p2*t1*y - k2*p2*p2*t1*y - 
  l2*n2*s2*t1*y + j2*p2*s2*t1*y + l2*o2*p1*t2*y + l2*o1*p2*t2*y + 
  l1*o2*p2*t2*y - 2*k2*p1*p2*t2*y - k1*p2*p2*t2*y - l2*n2*s1*t2*y + 
  j2*p2*s1*t2*y - l2*n1*s2*t2*y - l1*n2*s2*t2*y + j2*p1*s2*t2*y + 
  j1*p2*s2*t2*y + m2*o2*o2*u1*y - k2*o2*q2*u1*y - m2*n2*r2*u1*y + 
  j2*q2*r2*u1*y + k2*n2*t2*u1*y - j2*o2*t2*u1*y + 2*m2*o1*o2*u2*y + 
  m1*o2*o2*u2*y - k2*o2*q1*u2*y - k2*o1*q2*u2*y - k1*o2*q2*u2*y - 
  m2*n2*r1*u2*y + j2*q2*r1*u2*y - m2*n1*r2*u2*y - m1*n2*r2*u2*y + 
  j2*q1*r2*u2*y + j1*q2*r2*u2*y + k2*n2*t1*u2*y - j2*o2*t1*u2*y + 
  k2*n1*t2*u2*y + k1*n2*t2*u2*y - j2*o1*t2*u2*y - j1*o2*t2*u2*y - 
  l2*o2*o2*v1*y + k2*o2*p2*v1*y + l2*n2*r2*v1*y - j2*p2*r2*v1*y - 
  k2*n2*s2*v1*y + j2*o2*s2*v1*y - 2*l2*o1*o2*v2*y - l1*o2*o2*v2*y + 
  k2*o2*p1*v2*y + k2*o1*p2*v2*y + k1*o2*p2*v2*y + l2*n2*r1*v2*y - 
  j2*p2*r1*v2*y + l2*n1*r2*v2*y + l1*n2*r2*v2*y - j2*p1*r2*v2*y - 
  j1*p2*r2*v2*y - k2*n2*s1*v2*y + j2*o2*s1*v2*y - k2*n1*s2*v2*y - 
  k1*n2*s2*v2*y + j2*o1*s2*v2*y + j1*o2*s2*v2*y + 2*m2*p1*p2*r0*x*y + 
  m1*p2*p2*r0*x*y - l2*p2*q1*r0*x*y - l2*p1*q2*r0*x*y - l1*p2*q2*r0*x*y + 
  2*m2*p0*p2*r1*x*y + m0*p2*p2*r1*x*y - l2*p2*q0*r1*x*y - l2*p0*q2*r1*x*y - 
  l0*p2*q2*r1*x*y + 2*m2*p0*p1*r2*x*y + 2*m1*p0*p2*r2*x*y + 
  2*m0*p1*p2*r2*x*y - l2*p1*q0*r2*x*y - l1*p2*q0*r2*x*y - l2*p0*q1*r2*x*y - 
  l0*p2*q1*r2*x*y - l1*p0*q2*r2*x*y - l0*p1*q2*r2*x*y - 2*m2*o2*p1*s0*x*y - 
  2*m2*o1*p2*s0*x*y - 2*m1*o2*p2*s0*x*y + l2*o2*q1*s0*x*y + k2*p2*q1*s0*x*y + 
  l2*o1*q2*s0*x*y + l1*o2*q2*s0*x*y + k2*p1*q2*s0*x*y + k1*p2*q2*s0*x*y - 
  2*m2*o2*p0*s1*x*y - 2*m2*o0*p2*s1*x*y - 2*m0*o2*p2*s1*x*y + 
  l2*o2*q0*s1*x*y + k2*p2*q0*s1*x*y + l2*o0*q2*s1*x*y + l0*o2*q2*s1*x*y + 
  k2*p0*q2*s1*x*y + k0*p2*q2*s1*x*y + 2*m2*n2*s0*s1*x*y - 2*j2*q2*s0*s1*x*y - 
  2*m2*o1*p0*s2*x*y - 2*m1*o2*p0*s2*x*y - 2*m2*o0*p1*s2*x*y - 
  2*m0*o2*p1*s2*x*y - 2*m1*o0*p2*s2*x*y - 2*m0*o1*p2*s2*x*y + 
  l2*o1*q0*s2*x*y + l1*o2*q0*s2*x*y + k2*p1*q0*s2*x*y + k1*p2*q0*s2*x*y + 
  l2*o0*q1*s2*x*y + l0*o2*q1*s2*x*y + k2*p0*q1*s2*x*y + k0*p2*q1*s2*x*y + 
  l1*o0*q2*s2*x*y + l0*o1*q2*s2*x*y + k1*p0*q2*s2*x*y + k0*p1*q2*s2*x*y + 
  2*m2*n1*s0*s2*x*y + 2*m1*n2*s0*s2*x*y - 2*j2*q1*s0*s2*x*y - 
  2*j1*q2*s0*s2*x*y + 2*m2*n0*s1*s2*x*y + 2*m0*n2*s1*s2*x*y - 
  2*j2*q0*s1*s2*x*y - 2*j0*q2*s1*s2*x*y + m1*n0*s2*s2*x*y + m0*n1*s2*s2*x*y - 
  j1*q0*s2*s2*x*y - j0*q1*s2*s2*x*y + l2*o2*p1*t0*x*y + l2*o1*p2*t0*x*y + 
  l1*o2*p2*t0*x*y - 2*k2*p1*p2*t0*x*y - k1*p2*p2*t0*x*y - l2*n2*s1*t0*x*y + 
  j2*p2*s1*t0*x*y - l2*n1*s2*t0*x*y - l1*n2*s2*t0*x*y + j2*p1*s2*t0*x*y + 
  j1*p2*s2*t0*x*y + l2*o2*p0*t1*x*y + l2*o0*p2*t1*x*y + l0*o2*p2*t1*x*y - 
  2*k2*p0*p2*t1*x*y - k0*p2*p2*t1*x*y - l2*n2*s0*t1*x*y + j2*p2*s0*t1*x*y - 
  l2*n0*s2*t1*x*y - l0*n2*s2*t1*x*y + j2*p0*s2*t1*x*y + j0*p2*s2*t1*x*y + 
  l2*o1*p0*t2*x*y + l1*o2*p0*t2*x*y + l2*o0*p1*t2*x*y + l0*o2*p1*t2*x*y - 
  2*k2*p0*p1*t2*x*y + l1*o0*p2*t2*x*y + l0*o1*p2*t2*x*y - 2*k1*p0*p2*t2*x*y - 
  2*k0*p1*p2*t2*x*y - l2*n1*s0*t2*x*y - l1*n2*s0*t2*x*y + j2*p1*s0*t2*x*y + 
  j1*p2*s0*t2*x*y - l2*n0*s1*t2*x*y - l0*n2*s1*t2*x*y + j2*p0*s1*t2*x*y + 
  j0*p2*s1*t2*x*y - l1*n0*s2*t2*x*y - l0*n1*s2*t2*x*y + j1*p0*s2*t2*x*y + 
  j0*p1*s2*t2*x*y + 2*m2*o1*o2*u0*x*y + m1*o2*o2*u0*x*y - k2*o2*q1*u0*x*y - 
  k2*o1*q2*u0*x*y - k1*o2*q2*u0*x*y - m2*n2*r1*u0*x*y + j2*q2*r1*u0*x*y - 
  m2*n1*r2*u0*x*y - m1*n2*r2*u0*x*y + j2*q1*r2*u0*x*y + j1*q2*r2*u0*x*y + 
  k2*n2*t1*u0*x*y - j2*o2*t1*u0*x*y + k2*n1*t2*u0*x*y + k1*n2*t2*u0*x*y - 
  j2*o1*t2*u0*x*y - j1*o2*t2*u0*x*y + 2*m2*o0*o2*u1*x*y + m0*o2*o2*u1*x*y - 
  k2*o2*q0*u1*x*y - k2*o0*q2*u1*x*y - k0*o2*q2*u1*x*y - m2*n2*r0*u1*x*y + 
  j2*q2*r0*u1*x*y - m2*n0*r2*u1*x*y - m0*n2*r2*u1*x*y + j2*q0*r2*u1*x*y + 
  j0*q2*r2*u1*x*y + k2*n2*t0*u1*x*y - j2*o2*t0*u1*x*y + k2*n0*t2*u1*x*y + 
  k0*n2*t2*u1*x*y - j2*o0*t2*u1*x*y - j0*o2*t2*u1*x*y + 2*m2*o0*o1*u2*x*y + 
  2*m1*o0*o2*u2*x*y + 2*m0*o1*o2*u2*x*y - k2*o1*q0*u2*x*y - k1*o2*q0*u2*x*y - 
  k2*o0*q1*u2*x*y - k0*o2*q1*u2*x*y - k1*o0*q2*u2*x*y - k0*o1*q2*u2*x*y - 
  m2*n1*r0*u2*x*y - m1*n2*r0*u2*x*y + j2*q1*r0*u2*x*y + j1*q2*r0*u2*x*y - 
  m2*n0*r1*u2*x*y - m0*n2*r1*u2*x*y + j2*q0*r1*u2*x*y + j0*q2*r1*u2*x*y - 
  m1*n0*r2*u2*x*y - m0*n1*r2*u2*x*y + j1*q0*r2*u2*x*y + j0*q1*r2*u2*x*y + 
  k2*n1*t0*u2*x*y + k1*n2*t0*u2*x*y - j2*o1*t0*u2*x*y - j1*o2*t0*u2*x*y + 
  k2*n0*t1*u2*x*y + k0*n2*t1*u2*x*y - j2*o0*t1*u2*x*y - j0*o2*t1*u2*x*y + 
  k1*n0*t2*u2*x*y + k0*n1*t2*u2*x*y - j1*o0*t2*u2*x*y - j0*o1*t2*u2*x*y - 
  2*l2*o1*o2*v0*x*y - l1*o2*o2*v0*x*y + k2*o2*p1*v0*x*y + k2*o1*p2*v0*x*y + 
  k1*o2*p2*v0*x*y + l2*n2*r1*v0*x*y - j2*p2*r1*v0*x*y + l2*n1*r2*v0*x*y + 
  l1*n2*r2*v0*x*y - j2*p1*r2*v0*x*y - j1*p2*r2*v0*x*y - k2*n2*s1*v0*x*y + 
  j2*o2*s1*v0*x*y - k2*n1*s2*v0*x*y - k1*n2*s2*v0*x*y + j2*o1*s2*v0*x*y + 
  j1*o2*s2*v0*x*y - 2*l2*o0*o2*v1*x*y - l0*o2*o2*v1*x*y + k2*o2*p0*v1*x*y + 
  k2*o0*p2*v1*x*y + k0*o2*p2*v1*x*y + l2*n2*r0*v1*x*y - j2*p2*r0*v1*x*y + 
  l2*n0*r2*v1*x*y + l0*n2*r2*v1*x*y - j2*p0*r2*v1*x*y - j0*p2*r2*v1*x*y - 
  k2*n2*s0*v1*x*y + j2*o2*s0*v1*x*y - k2*n0*s2*v1*x*y - k0*n2*s2*v1*x*y + 
  j2*o0*s2*v1*x*y + j0*o2*s2*v1*x*y - 2*l2*o0*o1*v2*x*y - 2*l1*o0*o2*v2*x*y - 
  2*l0*o1*o2*v2*x*y + k2*o1*p0*v2*x*y + k1*o2*p0*v2*x*y + k2*o0*p1*v2*x*y + 
  k0*o2*p1*v2*x*y + k1*o0*p2*v2*x*y + k0*o1*p2*v2*x*y + l2*n1*r0*v2*x*y + 
  l1*n2*r0*v2*x*y - j2*p1*r0*v2*x*y - j1*p2*r0*v2*x*y + l2*n0*r1*v2*x*y + 
  l0*n2*r1*v2*x*y - j2*p0*r1*v2*x*y - j0*p2*r1*v2*x*y + l1*n0*r2*v2*x*y + 
  l0*n1*r2*v2*x*y - j1*p0*r2*v2*x*y - j0*p1*r2*v2*x*y - k2*n1*s0*v2*x*y - 
  k1*n2*s0*v2*x*y + j2*o1*s0*v2*x*y + j1*o2*s0*v2*x*y - k2*n0*s1*v2*x*y - 
  k0*n2*s1*v2*x*y + j2*o0*s1*v2*x*y + j0*o2*s1*v2*x*y - k1*n0*s2*v2*x*y - 
  k0*n1*s2*v2*x*y + j1*o0*s2*v2*x*y + j0*o1*s2*v2*x*y + 2*m2*p0*p1*r0*x*x*y + 
  2*m1*p0*p2*r0*x*x*y + 2*m0*p1*p2*r0*x*x*y - l2*p1*q0*r0*x*x*y - 
  l1*p2*q0*r0*x*x*y - l2*p0*q1*r0*x*x*y - l0*p2*q1*r0*x*x*y - 
  l1*p0*q2*r0*x*x*y - l0*p1*q2*r0*x*x*y + m2*p0*p0*r1*x*x*y + 
  2*m0*p0*p2*r1*x*x*y - l2*p0*q0*r1*x*x*y - l0*p2*q0*r1*x*x*y - 
  l0*p0*q2*r1*x*x*y + m1*p0*p0*r2*x*x*y + 2*m0*p0*p1*r2*x*x*y - 
  l1*p0*q0*r2*x*x*y - l0*p1*q0*r2*x*x*y - l0*p0*q1*r2*x*x*y - 
  2*m2*o1*p0*s0*x*x*y - 2*m1*o2*p0*s0*x*x*y - 2*m2*o0*p1*s0*x*x*y - 
  2*m0*o2*p1*s0*x*x*y - 2*m1*o0*p2*s0*x*x*y - 2*m0*o1*p2*s0*x*x*y + 
  l2*o1*q0*s0*x*x*y + l1*o2*q0*s0*x*x*y + k2*p1*q0*s0*x*x*y + 
  k1*p2*q0*s0*x*x*y + l2*o0*q1*s0*x*x*y + l0*o2*q1*s0*x*x*y + 
  k2*p0*q1*s0*x*x*y + k0*p2*q1*s0*x*x*y + l1*o0*q2*s0*x*x*y + 
  l0*o1*q2*s0*x*x*y + k1*p0*q2*s0*x*x*y + k0*p1*q2*s0*x*x*y + 
  m2*n1*s0*s0*x*x*y + m1*n2*s0*s0*x*x*y - j2*q1*s0*s0*x*x*y - j1*q2*s0*s0*x*x*y - 
  2*m2*o0*p0*s1*x*x*y - 2*m0*o2*p0*s1*x*x*y - 2*m0*o0*p2*s1*x*x*y + 
  l2*o0*q0*s1*x*x*y + l0*o2*q0*s1*x*x*y + k2*p0*q0*s1*x*x*y + 
  k0*p2*q0*s1*x*x*y + l0*o0*q2*s1*x*x*y + k0*p0*q2*s1*x*x*y + 
  2*m2*n0*s0*s1*x*x*y + 2*m0*n2*s0*s1*x*x*y - 2*j2*q0*s0*s1*x*x*y - 
  2*j0*q2*s0*s1*x*x*y - 2*m1*o0*p0*s2*x*x*y - 2*m0*o1*p0*s2*x*x*y - 
  2*m0*o0*p1*s2*x*x*y + l1*o0*q0*s2*x*x*y + l0*o1*q0*s2*x*x*y + 
  k1*p0*q0*s2*x*x*y + k0*p1*q0*s2*x*x*y + l0*o0*q1*s2*x*x*y + 
  k0*p0*q1*s2*x*x*y + 2*m1*n0*s0*s2*x*x*y + 2*m0*n1*s0*s2*x*x*y - 
  2*j1*q0*s0*s2*x*x*y - 2*j0*q1*s0*s2*x*x*y + 2*m0*n0*s1*s2*x*x*y - 
  2*j0*q0*s1*s2*x*x*y + l2*o1*p0*t0*x*x*y + l1*o2*p0*t0*x*x*y + 
  l2*o0*p1*t0*x*x*y + l0*o2*p1*t0*x*x*y - 2*k2*p0*p1*t0*x*x*y + 
  l1*o0*p2*t0*x*x*y + l0*o1*p2*t0*x*x*y - 2*k1*p0*p2*t0*x*x*y - 
  2*k0*p1*p2*t0*x*x*y - l2*n1*s0*t0*x*x*y - l1*n2*s0*t0*x*x*y + 
  j2*p1*s0*t0*x*x*y + j1*p2*s0*t0*x*x*y - l2*n0*s1*t0*x*x*y - 
  l0*n2*s1*t0*x*x*y + j2*p0*s1*t0*x*x*y + j0*p2*s1*t0*x*x*y - 
  l1*n0*s2*t0*x*x*y - l0*n1*s2*t0*x*x*y + j1*p0*s2*t0*x*x*y + 
  j0*p1*s2*t0*x*x*y + l2*o0*p0*t1*x*x*y + l0*o2*p0*t1*x*x*y - 
  k2*p0*p0*t1*x*x*y + l0*o0*p2*t1*x*x*y - 2*k0*p0*p2*t1*x*x*y - 
  l2*n0*s0*t1*x*x*y - l0*n2*s0*t1*x*x*y + j2*p0*s0*t1*x*x*y + 
  j0*p2*s0*t1*x*x*y - l0*n0*s2*t1*x*x*y + j0*p0*s2*t1*x*x*y + 
  l1*o0*p0*t2*x*x*y + l0*o1*p0*t2*x*x*y - k1*p0*p0*t2*x*x*y + 
  l0*o0*p1*t2*x*x*y - 2*k0*p0*p1*t2*x*x*y - l1*n0*s0*t2*x*x*y - 
  l0*n1*s0*t2*x*x*y + j1*p0*s0*t2*x*x*y + j0*p1*s0*t2*x*x*y - 
  l0*n0*s1*t2*x*x*y + j0*p0*s1*t2*x*x*y + 2*m2*o0*o1*u0*x*x*y + 
  2*m1*o0*o2*u0*x*x*y + 2*m0*o1*o2*u0*x*x*y - k2*o1*q0*u0*x*x*y - 
  k1*o2*q0*u0*x*x*y - k2*o0*q1*u0*x*x*y - k0*o2*q1*u0*x*x*y - 
  k1*o0*q2*u0*x*x*y - k0*o1*q2*u0*x*x*y - m2*n1*r0*u0*x*x*y - 
  m1*n2*r0*u0*x*x*y + j2*q1*r0*u0*x*x*y + j1*q2*r0*u0*x*x*y - 
  m2*n0*r1*u0*x*x*y - m0*n2*r1*u0*x*x*y + j2*q0*r1*u0*x*x*y + 
  j0*q2*r1*u0*x*x*y - m1*n0*r2*u0*x*x*y - m0*n1*r2*u0*x*x*y + 
  j1*q0*r2*u0*x*x*y + j0*q1*r2*u0*x*x*y + k2*n1*t0*u0*x*x*y + 
  k1*n2*t0*u0*x*x*y - j2*o1*t0*u0*x*x*y - j1*o2*t0*u0*x*x*y + 
  k2*n0*t1*u0*x*x*y + k0*n2*t1*u0*x*x*y - j2*o0*t1*u0*x*x*y - 
  j0*o2*t1*u0*x*x*y + k1*n0*t2*u0*x*x*y + k0*n1*t2*u0*x*x*y - 
  j1*o0*t2*u0*x*x*y - j0*o1*t2*u0*x*x*y + m2*o0*o0*u1*x*x*y + 
  2*m0*o0*o2*u1*x*x*y - k2*o0*q0*u1*x*x*y - k0*o2*q0*u1*x*x*y - 
  k0*o0*q2*u1*x*x*y - m2*n0*r0*u1*x*x*y - m0*n2*r0*u1*x*x*y + 
  j2*q0*r0*u1*x*x*y + j0*q2*r0*u1*x*x*y - m0*n0*r2*u1*x*x*y + 
  j0*q0*r2*u1*x*x*y + k2*n0*t0*u1*x*x*y + k0*n2*t0*u1*x*x*y - 
  j2*o0*t0*u1*x*x*y - j0*o2*t0*u1*x*x*y + k0*n0*t2*u1*x*x*y - 
  j0*o0*t2*u1*x*x*y + m1*o0*o0*u2*x*x*y + 2*m0*o0*o1*u2*x*x*y - 
  k1*o0*q0*u2*x*x*y - k0*o1*q0*u2*x*x*y - k0*o0*q1*u2*x*x*y - 
  m1*n0*r0*u2*x*x*y - m0*n1*r0*u2*x*x*y + j1*q0*r0*u2*x*x*y + 
  j0*q1*r0*u2*x*x*y - m0*n0*r1*u2*x*x*y + j0*q0*r1*u2*x*x*y + 
  k1*n0*t0*u2*x*x*y + k0*n1*t0*u2*x*x*y - j1*o0*t0*u2*x*x*y - 
  j0*o1*t0*u2*x*x*y + k0*n0*t1*u2*x*x*y - j0*o0*t1*u2*x*x*y - 
  2*l2*o0*o1*v0*x*x*y - 2*l1*o0*o2*v0*x*x*y - 2*l0*o1*o2*v0*x*x*y + 
  k2*o1*p0*v0*x*x*y + k1*o2*p0*v0*x*x*y + k2*o0*p1*v0*x*x*y + 
  k0*o2*p1*v0*x*x*y + k1*o0*p2*v0*x*x*y + k0*o1*p2*v0*x*x*y + 
  l2*n1*r0*v0*x*x*y + l1*n2*r0*v0*x*x*y - j2*p1*r0*v0*x*x*y - 
  j1*p2*r0*v0*x*x*y + l2*n0*r1*v0*x*x*y + l0*n2*r1*v0*x*x*y - 
  j2*p0*r1*v0*x*x*y - j0*p2*r1*v0*x*x*y + l1*n0*r2*v0*x*x*y + 
  l0*n1*r2*v0*x*x*y - j1*p0*r2*v0*x*x*y - j0*p1*r2*v0*x*x*y - 
  k2*n1*s0*v0*x*x*y - k1*n2*s0*v0*x*x*y + j2*o1*s0*v0*x*x*y + 
  j1*o2*s0*v0*x*x*y - k2*n0*s1*v0*x*x*y - k0*n2*s1*v0*x*x*y + 
  j2*o0*s1*v0*x*x*y + j0*o2*s1*v0*x*x*y - k1*n0*s2*v0*x*x*y - 
  k0*n1*s2*v0*x*x*y + j1*o0*s2*v0*x*x*y + j0*o1*s2*v0*x*x*y - 
  l2*o0*o0*v1*x*x*y - 2*l0*o0*o2*v1*x*x*y + k2*o0*p0*v1*x*x*y + 
  k0*o2*p0*v1*x*x*y + k0*o0*p2*v1*x*x*y + l2*n0*r0*v1*x*x*y + 
  l0*n2*r0*v1*x*x*y - j2*p0*r0*v1*x*x*y - j0*p2*r0*v1*x*x*y + 
  l0*n0*r2*v1*x*x*y - j0*p0*r2*v1*x*x*y - k2*n0*s0*v1*x*x*y - 
  k0*n2*s0*v1*x*x*y + j2*o0*s0*v1*x*x*y + j0*o2*s0*v1*x*x*y - 
  k0*n0*s2*v1*x*x*y + j0*o0*s2*v1*x*x*y - l1*o0*o0*v2*x*x*y - 
  2*l0*o0*o1*v2*x*x*y + k1*o0*p0*v2*x*x*y + k0*o1*p0*v2*x*x*y + 
  k0*o0*p1*v2*x*x*y + l1*n0*r0*v2*x*x*y + l0*n1*r0*v2*x*x*y - 
  j1*p0*r0*v2*x*x*y - j0*p1*r0*v2*x*x*y + l0*n0*r1*v2*x*x*y - 
  j0*p0*r1*v2*x*x*y - k1*n0*s0*v2*x*x*y - k0*n1*s0*v2*x*x*y + 
  j1*o0*s0*v2*x*x*y + j0*o1*s0*v2*x*x*y - k0*n0*s1*v2*x*x*y + 
  j0*o0*s1*v2*x*x*y + m1*p0*p0*r0*x*x*x*y + 2*m0*p0*p1*r0*x*x*x*y - 
  l1*p0*q0*r0*x*x*x*y - l0*p1*q0*r0*x*x*x*y - l0*p0*q1*r0*x*x*x*y + 
  m0*p0*p0*r1*x*x*x*y - l0*p0*q0*r1*x*x*x*y - 2*m1*o0*p0*s0*x*x*x*y - 
  2*m0*o1*p0*s0*x*x*x*y - 2*m0*o0*p1*s0*x*x*x*y + l1*o0*q0*s0*x*x*x*y + 
  l0*o1*q0*s0*x*x*x*y + k1*p0*q0*s0*x*x*x*y + k0*p1*q0*s0*x*x*x*y + 
  l0*o0*q1*s0*x*x*x*y + k0*p0*q1*s0*x*x*x*y + m1*n0*s0*s0*x*x*x*y + 
  m0*n1*s0*s0*x*x*x*y - j1*q0*s0*s0*x*x*x*y - j0*q1*s0*s0*x*x*x*y - 
  2*m0*o0*p0*s1*x*x*x*y + l0*o0*q0*s1*x*x*x*y + k0*p0*q0*s1*x*x*x*y + 
  2*m0*n0*s0*s1*x*x*x*y - 2*j0*q0*s0*s1*x*x*x*y + l1*o0*p0*t0*x*x*x*y + 
  l0*o1*p0*t0*x*x*x*y - k1*p0*p0*t0*x*x*x*y + l0*o0*p1*t0*x*x*x*y - 
  2*k0*p0*p1*t0*x*x*x*y - l1*n0*s0*t0*x*x*x*y - l0*n1*s0*t0*x*x*x*y + 
  j1*p0*s0*t0*x*x*x*y + j0*p1*s0*t0*x*x*x*y - l0*n0*s1*t0*x*x*x*y + 
  j0*p0*s1*t0*x*x*x*y + l0*o0*p0*t1*x*x*x*y - k0*p0*p0*t1*x*x*x*y - 
  l0*n0*s0*t1*x*x*x*y + j0*p0*s0*t1*x*x*x*y + m1*o0*o0*u0*x*x*x*y + 
  2*m0*o0*o1*u0*x*x*x*y - k1*o0*q0*u0*x*x*x*y - k0*o1*q0*u0*x*x*x*y - 
  k0*o0*q1*u0*x*x*x*y - m1*n0*r0*u0*x*x*x*y - m0*n1*r0*u0*x*x*x*y + 
  j1*q0*r0*u0*x*x*x*y + j0*q1*r0*u0*x*x*x*y - m0*n0*r1*u0*x*x*x*y + 
  j0*q0*r1*u0*x*x*x*y + k1*n0*t0*u0*x*x*x*y + k0*n1*t0*u0*x*x*x*y - 
  j1*o0*t0*u0*x*x*x*y - j0*o1*t0*u0*x*x*x*y + k0*n0*t1*u0*x*x*x*y - 
  j0*o0*t1*u0*x*x*x*y + m0*o0*o0*u1*x*x*x*y - k0*o0*q0*u1*x*x*x*y - 
  m0*n0*r0*u1*x*x*x*y + j0*q0*r0*u1*x*x*x*y + k0*n0*t0*u1*x*x*x*y - 
  j0*o0*t0*u1*x*x*x*y - l1*o0*o0*v0*x*x*x*y - 2*l0*o0*o1*v0*x*x*x*y + 
  k1*o0*p0*v0*x*x*x*y + k0*o1*p0*v0*x*x*x*y + k0*o0*p1*v0*x*x*x*y + 
  l1*n0*r0*v0*x*x*x*y + l0*n1*r0*v0*x*x*x*y - j1*p0*r0*v0*x*x*x*y - 
  j0*p1*r0*v0*x*x*x*y + l0*n0*r1*v0*x*x*x*y - j0*p0*r1*v0*x*x*x*y - 
  k1*n0*s0*v0*x*x*x*y - k0*n1*s0*v0*x*x*x*y + j1*o0*s0*v0*x*x*x*y + 
  j0*o1*s0*v0*x*x*x*y - k0*n0*s1*v0*x*x*x*y + j0*o0*s1*v0*x*x*x*y - 
  l0*o0*o0*v1*x*x*x*y + k0*o0*p0*v1*x*x*x*y + l0*n0*r0*v1*x*x*x*y - 
  j0*p0*r0*v1*x*x*x*y - k0*n0*s0*v1*x*x*x*y + j0*o0*s0*v1*x*x*x*y + 
  2*m2*p1*p2*r1*y*y + m1*p2*p2*r1*y*y - l2*p2*q1*r1*y*y - l2*p1*q2*r1*y*y - 
  l1*p2*q2*r1*y*y + m2*p1*p1*r2*y*y + 2*m1*p1*p2*r2*y*y - l2*p1*q1*r2*y*y - 
  l1*p2*q1*r2*y*y - l1*p1*q2*r2*y*y - 2*m2*o2*p1*s1*y*y - 2*m2*o1*p2*s1*y*y - 
  2*m1*o2*p2*s1*y*y + l2*o2*q1*s1*y*y + k2*p2*q1*s1*y*y + l2*o1*q2*s1*y*y + 
  l1*o2*q2*s1*y*y + k2*p1*q2*s1*y*y + k1*p2*q2*s1*y*y + m2*n2*s1*s1*y*y - 
  j2*q2*s1*s1*y*y - 2*m2*o1*p1*s2*y*y - 2*m1*o2*p1*s2*y*y - 
  2*m1*o1*p2*s2*y*y + l2*o1*q1*s2*y*y + l1*o2*q1*s2*y*y + k2*p1*q1*s2*y*y + 
  k1*p2*q1*s2*y*y + l1*o1*q2*s2*y*y + k1*p1*q2*s2*y*y + 2*m2*n1*s1*s2*y*y + 
  2*m1*n2*s1*s2*y*y - 2*j2*q1*s1*s2*y*y - 2*j1*q2*s1*s2*y*y + 
  m1*n1*s2*s2*y*y - j1*q1*s2*s2*y*y + l2*o2*p1*t1*y*y + l2*o1*p2*t1*y*y + 
  l1*o2*p2*t1*y*y - 2*k2*p1*p2*t1*y*y - k1*p2*p2*t1*y*y - l2*n2*s1*t1*y*y + 
  j2*p2*s1*t1*y*y - l2*n1*s2*t1*y*y - l1*n2*s2*t1*y*y + j2*p1*s2*t1*y*y + 
  j1*p2*s2*t1*y*y + l2*o1*p1*t2*y*y + l1*o2*p1*t2*y*y - k2*p1*p1*t2*y*y + 
  l1*o1*p2*t2*y*y - 2*k1*p1*p2*t2*y*y - l2*n1*s1*t2*y*y - l1*n2*s1*t2*y*y + 
  j2*p1*s1*t2*y*y + j1*p2*s1*t2*y*y - l1*n1*s2*t2*y*y + j1*p1*s2*t2*y*y + 
  2*m2*o1*o2*u1*y*y + m1*o2*o2*u1*y*y - k2*o2*q1*u1*y*y - k2*o1*q2*u1*y*y - 
  k1*o2*q2*u1*y*y - m2*n2*r1*u1*y*y + j2*q2*r1*u1*y*y - m2*n1*r2*u1*y*y - 
  m1*n2*r2*u1*y*y + j2*q1*r2*u1*y*y + j1*q2*r2*u1*y*y + k2*n2*t1*u1*y*y - 
  j2*o2*t1*u1*y*y + k2*n1*t2*u1*y*y + k1*n2*t2*u1*y*y - j2*o1*t2*u1*y*y - 
  j1*o2*t2*u1*y*y + m2*o1*o1*u2*y*y + 2*m1*o1*o2*u2*y*y - k2*o1*q1*u2*y*y - 
  k1*o2*q1*u2*y*y - k1*o1*q2*u2*y*y - m2*n1*r1*u2*y*y - m1*n2*r1*u2*y*y + 
  j2*q1*r1*u2*y*y + j1*q2*r1*u2*y*y - m1*n1*r2*u2*y*y + j1*q1*r2*u2*y*y + 
  k2*n1*t1*u2*y*y + k1*n2*t1*u2*y*y - j2*o1*t1*u2*y*y - j1*o2*t1*u2*y*y + 
  k1*n1*t2*u2*y*y - j1*o1*t2*u2*y*y - 2*l2*o1*o2*v1*y*y - l1*o2*o2*v1*y*y + 
  k2*o2*p1*v1*y*y + k2*o1*p2*v1*y*y + k1*o2*p2*v1*y*y + l2*n2*r1*v1*y*y - 
  j2*p2*r1*v1*y*y + l2*n1*r2*v1*y*y + l1*n2*r2*v1*y*y - j2*p1*r2*v1*y*y - 
  j1*p2*r2*v1*y*y - k2*n2*s1*v1*y*y + j2*o2*s1*v1*y*y - k2*n1*s2*v1*y*y - 
  k1*n2*s2*v1*y*y + j2*o1*s2*v1*y*y + j1*o2*s2*v1*y*y - l2*o1*o1*v2*y*y - 
  2*l1*o1*o2*v2*y*y + k2*o1*p1*v2*y*y + k1*o2*p1*v2*y*y + k1*o1*p2*v2*y*y + 
  l2*n1*r1*v2*y*y + l1*n2*r1*v2*y*y - j2*p1*r1*v2*y*y - j1*p2*r1*v2*y*y + 
  l1*n1*r2*v2*y*y - j1*p1*r2*v2*y*y - k2*n1*s1*v2*y*y - k1*n2*s1*v2*y*y + 
  j2*o1*s1*v2*y*y + j1*o2*s1*v2*y*y - k1*n1*s2*v2*y*y + j1*o1*s2*v2*y*y + 
  m2*p1*p1*r0*x*y*y + 2*m1*p1*p2*r0*x*y*y - l2*p1*q1*r0*x*y*y - 
  l1*p2*q1*r0*x*y*y - l1*p1*q2*r0*x*y*y + 2*m2*p0*p1*r1*x*y*y + 
  2*m1*p0*p2*r1*x*y*y + 2*m0*p1*p2*r1*x*y*y - l2*p1*q0*r1*x*y*y - 
  l1*p2*q0*r1*x*y*y - l2*p0*q1*r1*x*y*y - l0*p2*q1*r1*x*y*y - 
  l1*p0*q2*r1*x*y*y - l0*p1*q2*r1*x*y*y + 2*m1*p0*p1*r2*x*y*y + 
  m0*p1*p1*r2*x*y*y - l1*p1*q0*r2*x*y*y - l1*p0*q1*r2*x*y*y - 
  l0*p1*q1*r2*x*y*y - 2*m2*o1*p1*s0*x*y*y - 2*m1*o2*p1*s0*x*y*y - 
  2*m1*o1*p2*s0*x*y*y + l2*o1*q1*s0*x*y*y + l1*o2*q1*s0*x*y*y + 
  k2*p1*q1*s0*x*y*y + k1*p2*q1*s0*x*y*y + l1*o1*q2*s0*x*y*y + 
  k1*p1*q2*s0*x*y*y - 2*m2*o1*p0*s1*x*y*y - 2*m1*o2*p0*s1*x*y*y - 
  2*m2*o0*p1*s1*x*y*y - 2*m0*o2*p1*s1*x*y*y - 2*m1*o0*p2*s1*x*y*y - 
  2*m0*o1*p2*s1*x*y*y + l2*o1*q0*s1*x*y*y + l1*o2*q0*s1*x*y*y + 
  k2*p1*q0*s1*x*y*y + k1*p2*q0*s1*x*y*y + l2*o0*q1*s1*x*y*y + 
  l0*o2*q1*s1*x*y*y + k2*p0*q1*s1*x*y*y + k0*p2*q1*s1*x*y*y + 
  l1*o0*q2*s1*x*y*y + l0*o1*q2*s1*x*y*y + k1*p0*q2*s1*x*y*y + 
  k0*p1*q2*s1*x*y*y + 2*m2*n1*s0*s1*x*y*y + 2*m1*n2*s0*s1*x*y*y - 
  2*j2*q1*s0*s1*x*y*y - 2*j1*q2*s0*s1*x*y*y + m2*n0*s1*s1*x*y*y + 
  m0*n2*s1*s1*x*y*y - j2*q0*s1*s1*x*y*y - j0*q2*s1*s1*x*y*y - 
  2*m1*o1*p0*s2*x*y*y - 2*m1*o0*p1*s2*x*y*y - 2*m0*o1*p1*s2*x*y*y + 
  l1*o1*q0*s2*x*y*y + k1*p1*q0*s2*x*y*y + l1*o0*q1*s2*x*y*y + 
  l0*o1*q1*s2*x*y*y + k1*p0*q1*s2*x*y*y + k0*p1*q1*s2*x*y*y + 
  2*m1*n1*s0*s2*x*y*y - 2*j1*q1*s0*s2*x*y*y + 2*m1*n0*s1*s2*x*y*y + 
  2*m0*n1*s1*s2*x*y*y - 2*j1*q0*s1*s2*x*y*y - 2*j0*q1*s1*s2*x*y*y + 
  l2*o1*p1*t0*x*y*y + l1*o2*p1*t0*x*y*y - k2*p1*p1*t0*x*y*y + 
  l1*o1*p2*t0*x*y*y - 2*k1*p1*p2*t0*x*y*y - l2*n1*s1*t0*x*y*y - 
  l1*n2*s1*t0*x*y*y + j2*p1*s1*t0*x*y*y + j1*p2*s1*t0*x*y*y - 
  l1*n1*s2*t0*x*y*y + j1*p1*s2*t0*x*y*y + l2*o1*p0*t1*x*y*y + 
  l1*o2*p0*t1*x*y*y + l2*o0*p1*t1*x*y*y + l0*o2*p1*t1*x*y*y - 
  2*k2*p0*p1*t1*x*y*y + l1*o0*p2*t1*x*y*y + l0*o1*p2*t1*x*y*y - 
  2*k1*p0*p2*t1*x*y*y - 2*k0*p1*p2*t1*x*y*y - l2*n1*s0*t1*x*y*y - 
  l1*n2*s0*t1*x*y*y + j2*p1*s0*t1*x*y*y + j1*p2*s0*t1*x*y*y - 
  l2*n0*s1*t1*x*y*y - l0*n2*s1*t1*x*y*y + j2*p0*s1*t1*x*y*y + 
  j0*p2*s1*t1*x*y*y - l1*n0*s2*t1*x*y*y - l0*n1*s2*t1*x*y*y + 
  j1*p0*s2*t1*x*y*y + j0*p1*s2*t1*x*y*y + l1*o1*p0*t2*x*y*y + 
  l1*o0*p1*t2*x*y*y + l0*o1*p1*t2*x*y*y - 2*k1*p0*p1*t2*x*y*y - 
  k0*p1*p1*t2*x*y*y - l1*n1*s0*t2*x*y*y + j1*p1*s0*t2*x*y*y - 
  l1*n0*s1*t2*x*y*y - l0*n1*s1*t2*x*y*y + j1*p0*s1*t2*x*y*y + 
  j0*p1*s1*t2*x*y*y + m2*o1*o1*u0*x*y*y + 2*m1*o1*o2*u0*x*y*y - 
  k2*o1*q1*u0*x*y*y - k1*o2*q1*u0*x*y*y - k1*o1*q2*u0*x*y*y - 
  m2*n1*r1*u0*x*y*y - m1*n2*r1*u0*x*y*y + j2*q1*r1*u0*x*y*y + 
  j1*q2*r1*u0*x*y*y - m1*n1*r2*u0*x*y*y + j1*q1*r2*u0*x*y*y + 
  k2*n1*t1*u0*x*y*y + k1*n2*t1*u0*x*y*y - j2*o1*t1*u0*x*y*y - 
  j1*o2*t1*u0*x*y*y + k1*n1*t2*u0*x*y*y - j1*o1*t2*u0*x*y*y + 
  2*m2*o0*o1*u1*x*y*y + 2*m1*o0*o2*u1*x*y*y + 2*m0*o1*o2*u1*x*y*y - 
  k2*o1*q0*u1*x*y*y - k1*o2*q0*u1*x*y*y - k2*o0*q1*u1*x*y*y - 
  k0*o2*q1*u1*x*y*y - k1*o0*q2*u1*x*y*y - k0*o1*q2*u1*x*y*y - 
  m2*n1*r0*u1*x*y*y - m1*n2*r0*u1*x*y*y + j2*q1*r0*u1*x*y*y + 
  j1*q2*r0*u1*x*y*y - m2*n0*r1*u1*x*y*y - m0*n2*r1*u1*x*y*y + 
  j2*q0*r1*u1*x*y*y + j0*q2*r1*u1*x*y*y - m1*n0*r2*u1*x*y*y - 
  m0*n1*r2*u1*x*y*y + j1*q0*r2*u1*x*y*y + j0*q1*r2*u1*x*y*y + 
  k2*n1*t0*u1*x*y*y + k1*n2*t0*u1*x*y*y - j2*o1*t0*u1*x*y*y - 
  j1*o2*t0*u1*x*y*y + k2*n0*t1*u1*x*y*y + k0*n2*t1*u1*x*y*y - 
  j2*o0*t1*u1*x*y*y - j0*o2*t1*u1*x*y*y + k1*n0*t2*u1*x*y*y + 
  k0*n1*t2*u1*x*y*y - j1*o0*t2*u1*x*y*y - j0*o1*t2*u1*x*y*y + 
  2*m1*o0*o1*u2*x*y*y + m0*o1*o1*u2*x*y*y - k1*o1*q0*u2*x*y*y - 
  k1*o0*q1*u2*x*y*y - k0*o1*q1*u2*x*y*y - m1*n1*r0*u2*x*y*y + 
  j1*q1*r0*u2*x*y*y - m1*n0*r1*u2*x*y*y - m0*n1*r1*u2*x*y*y + 
  j1*q0*r1*u2*x*y*y + j0*q1*r1*u2*x*y*y + k1*n1*t0*u2*x*y*y - 
  j1*o1*t0*u2*x*y*y + k1*n0*t1*u2*x*y*y + k0*n1*t1*u2*x*y*y - 
  j1*o0*t1*u2*x*y*y - j0*o1*t1*u2*x*y*y - l2*o1*o1*v0*x*y*y - 
  2*l1*o1*o2*v0*x*y*y + k2*o1*p1*v0*x*y*y + k1*o2*p1*v0*x*y*y + 
  k1*o1*p2*v0*x*y*y + l2*n1*r1*v0*x*y*y + l1*n2*r1*v0*x*y*y - 
  j2*p1*r1*v0*x*y*y - j1*p2*r1*v0*x*y*y + l1*n1*r2*v0*x*y*y - 
  j1*p1*r2*v0*x*y*y - k2*n1*s1*v0*x*y*y - k1*n2*s1*v0*x*y*y + 
  j2*o1*s1*v0*x*y*y + j1*o2*s1*v0*x*y*y - k1*n1*s2*v0*x*y*y + 
  j1*o1*s2*v0*x*y*y - 2*l2*o0*o1*v1*x*y*y - 2*l1*o0*o2*v1*x*y*y - 
  2*l0*o1*o2*v1*x*y*y + k2*o1*p0*v1*x*y*y + k1*o2*p0*v1*x*y*y + 
  k2*o0*p1*v1*x*y*y + k0*o2*p1*v1*x*y*y + k1*o0*p2*v1*x*y*y + 
  k0*o1*p2*v1*x*y*y + l2*n1*r0*v1*x*y*y + l1*n2*r0*v1*x*y*y - 
  j2*p1*r0*v1*x*y*y - j1*p2*r0*v1*x*y*y + l2*n0*r1*v1*x*y*y + 
  l0*n2*r1*v1*x*y*y - j2*p0*r1*v1*x*y*y - j0*p2*r1*v1*x*y*y + 
  l1*n0*r2*v1*x*y*y + l0*n1*r2*v1*x*y*y - j1*p0*r2*v1*x*y*y - 
  j0*p1*r2*v1*x*y*y - k2*n1*s0*v1*x*y*y - k1*n2*s0*v1*x*y*y + 
  j2*o1*s0*v1*x*y*y + j1*o2*s0*v1*x*y*y - k2*n0*s1*v1*x*y*y - 
  k0*n2*s1*v1*x*y*y + j2*o0*s1*v1*x*y*y + j0*o2*s1*v1*x*y*y - 
  k1*n0*s2*v1*x*y*y - k0*n1*s2*v1*x*y*y + j1*o0*s2*v1*x*y*y + 
  j0*o1*s2*v1*x*y*y - 2*l1*o0*o1*v2*x*y*y - l0*o1*o1*v2*x*y*y + 
  k1*o1*p0*v2*x*y*y + k1*o0*p1*v2*x*y*y + k0*o1*p1*v2*x*y*y + 
  l1*n1*r0*v2*x*y*y - j1*p1*r0*v2*x*y*y + l1*n0*r1*v2*x*y*y + 
  l0*n1*r1*v2*x*y*y - j1*p0*r1*v2*x*y*y - j0*p1*r1*v2*x*y*y - 
  k1*n1*s0*v2*x*y*y + j1*o1*s0*v2*x*y*y - k1*n0*s1*v2*x*y*y - 
  k0*n1*s1*v2*x*y*y + j1*o0*s1*v2*x*y*y + j0*o1*s1*v2*x*y*y + 
  2*m1*p0*p1*r0*x*x*y*y + m0*p1*p1*r0*x*x*y*y - l1*p1*q0*r0*x*x*y*y - 
  l1*p0*q1*r0*x*x*y*y - l0*p1*q1*r0*x*x*y*y + m1*p0*p0*r1*x*x*y*y + 
  2*m0*p0*p1*r1*x*x*y*y - l1*p0*q0*r1*x*x*y*y - l0*p1*q0*r1*x*x*y*y - 
  l0*p0*q1*r1*x*x*y*y - 2*m1*o1*p0*s0*x*x*y*y - 2*m1*o0*p1*s0*x*x*y*y - 
  2*m0*o1*p1*s0*x*x*y*y + l1*o1*q0*s0*x*x*y*y + k1*p1*q0*s0*x*x*y*y + 
  l1*o0*q1*s0*x*x*y*y + l0*o1*q1*s0*x*x*y*y + k1*p0*q1*s0*x*x*y*y + 
  k0*p1*q1*s0*x*x*y*y + m1*n1*s0*s0*x*x*y*y - j1*q1*s0*s0*x*x*y*y - 
  2*m1*o0*p0*s1*x*x*y*y - 2*m0*o1*p0*s1*x*x*y*y - 2*m0*o0*p1*s1*x*x*y*y + 
  l1*o0*q0*s1*x*x*y*y + l0*o1*q0*s1*x*x*y*y + k1*p0*q0*s1*x*x*y*y + 
  k0*p1*q0*s1*x*x*y*y + l0*o0*q1*s1*x*x*y*y + k0*p0*q1*s1*x*x*y*y + 
  2*m1*n0*s0*s1*x*x*y*y + 2*m0*n1*s0*s1*x*x*y*y - 2*j1*q0*s0*s1*x*x*y*y - 
  2*j0*q1*s0*s1*x*x*y*y + m0*n0*s1*s1*x*x*y*y - j0*q0*s1*s1*x*x*y*y + 
  l1*o1*p0*t0*x*x*y*y + l1*o0*p1*t0*x*x*y*y + l0*o1*p1*t0*x*x*y*y - 
  2*k1*p0*p1*t0*x*x*y*y - k0*p1*p1*t0*x*x*y*y - l1*n1*s0*t0*x*x*y*y + 
  j1*p1*s0*t0*x*x*y*y - l1*n0*s1*t0*x*x*y*y - l0*n1*s1*t0*x*x*y*y + 
  j1*p0*s1*t0*x*x*y*y + j0*p1*s1*t0*x*x*y*y + l1*o0*p0*t1*x*x*y*y + 
  l0*o1*p0*t1*x*x*y*y - k1*p0*p0*t1*x*x*y*y + l0*o0*p1*t1*x*x*y*y - 
  2*k0*p0*p1*t1*x*x*y*y - l1*n0*s0*t1*x*x*y*y - l0*n1*s0*t1*x*x*y*y + 
  j1*p0*s0*t1*x*x*y*y + j0*p1*s0*t1*x*x*y*y - l0*n0*s1*t1*x*x*y*y + 
  j0*p0*s1*t1*x*x*y*y + 2*m1*o0*o1*u0*x*x*y*y + m0*o1*o1*u0*x*x*y*y - 
  k1*o1*q0*u0*x*x*y*y - k1*o0*q1*u0*x*x*y*y - k0*o1*q1*u0*x*x*y*y - 
  m1*n1*r0*u0*x*x*y*y + j1*q1*r0*u0*x*x*y*y - m1*n0*r1*u0*x*x*y*y - 
  m0*n1*r1*u0*x*x*y*y + j1*q0*r1*u0*x*x*y*y + j0*q1*r1*u0*x*x*y*y + 
  k1*n1*t0*u0*x*x*y*y - j1*o1*t0*u0*x*x*y*y + k1*n0*t1*u0*x*x*y*y + 
  k0*n1*t1*u0*x*x*y*y - j1*o0*t1*u0*x*x*y*y - j0*o1*t1*u0*x*x*y*y + 
  m1*o0*o0*u1*x*x*y*y + 2*m0*o0*o1*u1*x*x*y*y - k1*o0*q0*u1*x*x*y*y - 
  k0*o1*q0*u1*x*x*y*y - k0*o0*q1*u1*x*x*y*y - m1*n0*r0*u1*x*x*y*y - 
  m0*n1*r0*u1*x*x*y*y + j1*q0*r0*u1*x*x*y*y + j0*q1*r0*u1*x*x*y*y - 
  m0*n0*r1*u1*x*x*y*y + j0*q0*r1*u1*x*x*y*y + k1*n0*t0*u1*x*x*y*y + 
  k0*n1*t0*u1*x*x*y*y - j1*o0*t0*u1*x*x*y*y - j0*o1*t0*u1*x*x*y*y + 
  k0*n0*t1*u1*x*x*y*y - j0*o0*t1*u1*x*x*y*y - 2*l1*o0*o1*v0*x*x*y*y - 
  l0*o1*o1*v0*x*x*y*y + k1*o1*p0*v0*x*x*y*y + k1*o0*p1*v0*x*x*y*y + 
  k0*o1*p1*v0*x*x*y*y + l1*n1*r0*v0*x*x*y*y - j1*p1*r0*v0*x*x*y*y + 
  l1*n0*r1*v0*x*x*y*y + l0*n1*r1*v0*x*x*y*y - j1*p0*r1*v0*x*x*y*y - 
  j0*p1*r1*v0*x*x*y*y - k1*n1*s0*v0*x*x*y*y + j1*o1*s0*v0*x*x*y*y - 
  k1*n0*s1*v0*x*x*y*y - k0*n1*s1*v0*x*x*y*y + j1*o0*s1*v0*x*x*y*y + 
  j0*o1*s1*v0*x*x*y*y - l1*o0*o0*v1*x*x*y*y - 2*l0*o0*o1*v1*x*x*y*y + 
  k1*o0*p0*v1*x*x*y*y + k0*o1*p0*v1*x*x*y*y + k0*o0*p1*v1*x*x*y*y + 
  l1*n0*r0*v1*x*x*y*y + l0*n1*r0*v1*x*x*y*y - j1*p0*r0*v1*x*x*y*y - 
  j0*p1*r0*v1*x*x*y*y + l0*n0*r1*v1*x*x*y*y - j0*p0*r1*v1*x*x*y*y - 
  k1*n0*s0*v1*x*x*y*y - k0*n1*s0*v1*x*x*y*y + j1*o0*s0*v1*x*x*y*y + 
  j0*o1*s0*v1*x*x*y*y - k0*n0*s1*v1*x*x*y*y + j0*o0*s1*v1*x*x*y*y + 
  m2*p1*p1*r1*y*y*y + 2*m1*p1*p2*r1*y*y*y - l2*p1*q1*r1*y*y*y - l1*p2*q1*r1*y*y*y - 
  l1*p1*q2*r1*y*y*y + m1*p1*p1*r2*y*y*y - l1*p1*q1*r2*y*y*y - 2*m2*o1*p1*s1*y*y*y - 
  2*m1*o2*p1*s1*y*y*y - 2*m1*o1*p2*s1*y*y*y + l2*o1*q1*s1*y*y*y + l1*o2*q1*s1*y*y*y + 
  k2*p1*q1*s1*y*y*y + k1*p2*q1*s1*y*y*y + l1*o1*q2*s1*y*y*y + k1*p1*q2*s1*y*y*y + 
  m2*n1*s1*s1*y*y*y + m1*n2*s1*s1*y*y*y - j2*q1*s1*s1*y*y*y - j1*q2*s1*s1*y*y*y - 
  2*m1*o1*p1*s2*y*y*y + l1*o1*q1*s2*y*y*y + k1*p1*q1*s2*y*y*y + 2*m1*n1*s1*s2*y*y*y - 
  2*j1*q1*s1*s2*y*y*y + l2*o1*p1*t1*y*y*y + l1*o2*p1*t1*y*y*y - k2*p1*p1*t1*y*y*y + 
  l1*o1*p2*t1*y*y*y - 2*k1*p1*p2*t1*y*y*y - l2*n1*s1*t1*y*y*y - l1*n2*s1*t1*y*y*y + 
  j2*p1*s1*t1*y*y*y + j1*p2*s1*t1*y*y*y - l1*n1*s2*t1*y*y*y + j1*p1*s2*t1*y*y*y + 
  l1*o1*p1*t2*y*y*y - k1*p1*p1*t2*y*y*y - l1*n1*s1*t2*y*y*y + j1*p1*s1*t2*y*y*y + 
  m2*o1*o1*u1*y*y*y + 2*m1*o1*o2*u1*y*y*y - k2*o1*q1*u1*y*y*y - k1*o2*q1*u1*y*y*y - 
  k1*o1*q2*u1*y*y*y - m2*n1*r1*u1*y*y*y - m1*n2*r1*u1*y*y*y + j2*q1*r1*u1*y*y*y + 
  j1*q2*r1*u1*y*y*y - m1*n1*r2*u1*y*y*y + j1*q1*r2*u1*y*y*y + k2*n1*t1*u1*y*y*y + 
  k1*n2*t1*u1*y*y*y - j2*o1*t1*u1*y*y*y - j1*o2*t1*u1*y*y*y + k1*n1*t2*u1*y*y*y - 
  j1*o1*t2*u1*y*y*y + m1*o1*o1*u2*y*y*y - k1*o1*q1*u2*y*y*y - m1*n1*r1*u2*y*y*y + 
  j1*q1*r1*u2*y*y*y + k1*n1*t1*u2*y*y*y - j1*o1*t1*u2*y*y*y - l2*o1*o1*v1*y*y*y - 
  2*l1*o1*o2*v1*y*y*y + k2*o1*p1*v1*y*y*y + k1*o2*p1*v1*y*y*y + k1*o1*p2*v1*y*y*y + 
  l2*n1*r1*v1*y*y*y + l1*n2*r1*v1*y*y*y - j2*p1*r1*v1*y*y*y - j1*p2*r1*v1*y*y*y + 
  l1*n1*r2*v1*y*y*y - j1*p1*r2*v1*y*y*y - k2*n1*s1*v1*y*y*y - k1*n2*s1*v1*y*y*y + 
  j2*o1*s1*v1*y*y*y + j1*o2*s1*v1*y*y*y - k1*n1*s2*v1*y*y*y + j1*o1*s2*v1*y*y*y - 
  l1*o1*o1*v2*y*y*y + k1*o1*p1*v2*y*y*y + l1*n1*r1*v2*y*y*y - j1*p1*r1*v2*y*y*y - 
  k1*n1*s1*v2*y*y*y + j1*o1*s1*v2*y*y*y + m1*p1*p1*r0*x*y*y*y - l1*p1*q1*r0*x*y*y*y + 
  2*m1*p0*p1*r1*x*y*y*y + m0*p1*p1*r1*x*y*y*y - l1*p1*q0*r1*x*y*y*y - 
  l1*p0*q1*r1*x*y*y*y - l0*p1*q1*r1*x*y*y*y - 2*m1*o1*p1*s0*x*y*y*y + 
  l1*o1*q1*s0*x*y*y*y + k1*p1*q1*s0*x*y*y*y - 2*m1*o1*p0*s1*x*y*y*y - 
  2*m1*o0*p1*s1*x*y*y*y - 2*m0*o1*p1*s1*x*y*y*y + l1*o1*q0*s1*x*y*y*y + 
  k1*p1*q0*s1*x*y*y*y + l1*o0*q1*s1*x*y*y*y + l0*o1*q1*s1*x*y*y*y + 
  k1*p0*q1*s1*x*y*y*y + k0*p1*q1*s1*x*y*y*y + 2*m1*n1*s0*s1*x*y*y*y - 
  2*j1*q1*s0*s1*x*y*y*y + m1*n0*s1*s1*x*y*y*y + m0*n1*s1*s1*x*y*y*y - 
  j1*q0*s1*s1*x*y*y*y - j0*q1*s1*s1*x*y*y*y + l1*o1*p1*t0*x*y*y*y - 
  k1*p1*p1*t0*x*y*y*y - l1*n1*s1*t0*x*y*y*y + j1*p1*s1*t0*x*y*y*y + 
  l1*o1*p0*t1*x*y*y*y + l1*o0*p1*t1*x*y*y*y + l0*o1*p1*t1*x*y*y*y - 
  2*k1*p0*p1*t1*x*y*y*y - k0*p1*p1*t1*x*y*y*y - l1*n1*s0*t1*x*y*y*y + 
  j1*p1*s0*t1*x*y*y*y - l1*n0*s1*t1*x*y*y*y - l0*n1*s1*t1*x*y*y*y + 
  j1*p0*s1*t1*x*y*y*y + j0*p1*s1*t1*x*y*y*y + m1*o1*o1*u0*x*y*y*y - 
  k1*o1*q1*u0*x*y*y*y - m1*n1*r1*u0*x*y*y*y + j1*q1*r1*u0*x*y*y*y + 
  k1*n1*t1*u0*x*y*y*y - j1*o1*t1*u0*x*y*y*y + 2*m1*o0*o1*u1*x*y*y*y + 
  m0*o1*o1*u1*x*y*y*y - k1*o1*q0*u1*x*y*y*y - k1*o0*q1*u1*x*y*y*y - 
  k0*o1*q1*u1*x*y*y*y - m1*n1*r0*u1*x*y*y*y + j1*q1*r0*u1*x*y*y*y - 
  m1*n0*r1*u1*x*y*y*y - m0*n1*r1*u1*x*y*y*y + j1*q0*r1*u1*x*y*y*y + 
  j0*q1*r1*u1*x*y*y*y + k1*n1*t0*u1*x*y*y*y - j1*o1*t0*u1*x*y*y*y + 
  k1*n0*t1*u1*x*y*y*y + k0*n1*t1*u1*x*y*y*y - j1*o0*t1*u1*x*y*y*y - 
  j0*o1*t1*u1*x*y*y*y - l1*o1*o1*v0*x*y*y*y + k1*o1*p1*v0*x*y*y*y + 
  l1*n1*r1*v0*x*y*y*y - j1*p1*r1*v0*x*y*y*y - k1*n1*s1*v0*x*y*y*y + 
  j1*o1*s1*v0*x*y*y*y - 2*l1*o0*o1*v1*x*y*y*y - l0*o1*o1*v1*x*y*y*y + 
  k1*o1*p0*v1*x*y*y*y + k1*o0*p1*v1*x*y*y*y + k0*o1*p1*v1*x*y*y*y + 
  l1*n1*r0*v1*x*y*y*y - j1*p1*r0*v1*x*y*y*y + l1*n0*r1*v1*x*y*y*y + 
  l0*n1*r1*v1*x*y*y*y - j1*p0*r1*v1*x*y*y*y - j0*p1*r1*v1*x*y*y*y - 
  k1*n1*s0*v1*x*y*y*y + j1*o1*s0*v1*x*y*y*y - k1*n0*s1*v1*x*y*y*y - 
  k0*n1*s1*v1*x*y*y*y + j1*o0*s1*v1*x*y*y*y + j0*o1*s1*v1*x*y*y*y + 
  m1*p1*p1*r1*y*y*y*y - l1*p1*q1*r1*y*y*y*y - 2*m1*o1*p1*s1*y*y*y*y + l1*o1*q1*s1*y*y*y*y + 
  k1*p1*q1*s1*y*y*y*y + m1*n1*s1*s1*y*y*y*y - j1*q1*s1*s1*y*y*y*y + l1*o1*p1*t1*y*y*y*y - 
  k1*p1*p1*t1*y*y*y*y - l1*n1*s1*t1*y*y*y*y + j1*p1*s1*t1*y*y*y*y + m1*o1*o1*u1*y*y*y*y - 
  k1*o1*q1*u1*y*y*y*y - m1*n1*r1*u1*y*y*y*y + j1*q1*r1*u1*y*y*y*y + k1*n1*t1*u1*y*y*y*y - 
  j1*o1*t1*u1*y*y*y*y - l1*o1*o1*v1*y*y*y*y + k1*o1*p1*v1*y*y*y*y + l1*n1*r1*v1*y*y*y*y - 
  j1*p1*r1*v1*y*y*y*y - k1*n1*s1*v1*y*y*y*y + j1*o1*s1*v1*y*y*y*y;

	rg_REAL t_3  = l2*m2*p2*r2 - l2*l2*q2*r2 - l2*m2*o2*s2 - k2*m2*p2*s2 + 2*k2*l2*q2*s2 + 
  j2*m2*s2*s2 - i2*q2*s2*s2 + l2*l2*o2*t2 - k2*l2*p2*t2 - j2*l2*s2*t2 + 
  i2*p2*s2*t2 + k2*m2*o2*u2 - k2*k2*q2*u2 - j2*m2*r2*u2 + i2*q2*r2*u2 + 
  j2*k2*t2*u2 - i2*o2*t2*u2 - k2*l2*o2*v2 + k2*k2*p2*v2 + j2*l2*r2*v2 - 
  i2*p2*r2*v2 - j2*k2*s2*v2 + i2*o2*s2*v2 + l2*m2*p2*r0*x - l2*l2*q2*r0*x + 
  l2*m2*p0*r2*x + l2*m0*p2*r2*x + l0*m2*p2*r2*x - l2*l2*q0*r2*x - 
  2*l0*l2*q2*r2*x - l2*m2*o2*s0*x - k2*m2*p2*s0*x + 2*k2*l2*q2*s0*x - 
  l2*m2*o0*s2*x - l2*m0*o2*s2*x - l0*m2*o2*s2*x - k2*m2*p0*s2*x - 
  k2*m0*p2*s2*x - k0*m2*p2*s2*x + 2*k2*l2*q0*s2*x + 2*k2*l0*q2*s2*x + 
  2*k0*l2*q2*s2*x + 2*j2*m2*s0*s2*x - 2*i2*q2*s0*s2*x + j2*m0*s2*s2*x + 
  j0*m2*s2*s2*x - i2*q0*s2*s2*x - i0*q2*s2*s2*x + l2*l2*o2*t0*x - k2*l2*p2*t0*x - 
  j2*l2*s2*t0*x + i2*p2*s2*t0*x + l2*l2*o0*t2*x + 2*l0*l2*o2*t2*x - 
  k2*l2*p0*t2*x - k2*l0*p2*t2*x - k0*l2*p2*t2*x - j2*l2*s0*t2*x + 
  i2*p2*s0*t2*x - j2*l0*s2*t2*x - j0*l2*s2*t2*x + i2*p0*s2*t2*x + 
  i0*p2*s2*t2*x + k2*m2*o2*u0*x - k2*k2*q2*u0*x - j2*m2*r2*u0*x + 
  i2*q2*r2*u0*x + j2*k2*t2*u0*x - i2*o2*t2*u0*x + k2*m2*o0*u2*x + 
  k2*m0*o2*u2*x + k0*m2*o2*u2*x - k2*k2*q0*u2*x - 2*k0*k2*q2*u2*x - 
  j2*m2*r0*u2*x + i2*q2*r0*u2*x - j2*m0*r2*u2*x - j0*m2*r2*u2*x + 
  i2*q0*r2*u2*x + i0*q2*r2*u2*x + j2*k2*t0*u2*x - i2*o2*t0*u2*x + 
  j2*k0*t2*u2*x + j0*k2*t2*u2*x - i2*o0*t2*u2*x - i0*o2*t2*u2*x - 
  k2*l2*o2*v0*x + k2*k2*p2*v0*x + j2*l2*r2*v0*x - i2*p2*r2*v0*x - 
  j2*k2*s2*v0*x + i2*o2*s2*v0*x - k2*l2*o0*v2*x - k2*l0*o2*v2*x - 
  k0*l2*o2*v2*x + k2*k2*p0*v2*x + 2*k0*k2*p2*v2*x + j2*l2*r0*v2*x - 
  i2*p2*r0*v2*x + j2*l0*r2*v2*x + j0*l2*r2*v2*x - i2*p0*r2*v2*x - 
  i0*p2*r2*v2*x - j2*k2*s0*v2*x + i2*o2*s0*v2*x - j2*k0*s2*v2*x - 
  j0*k2*s2*v2*x + i2*o0*s2*v2*x + i0*o2*s2*v2*x + l2*m2*p0*r0*x*x + 
  l2*m0*p2*r0*x*x + l0*m2*p2*r0*x*x - l2*l2*q0*r0*x*x - 2*l0*l2*q2*r0*x*x + 
  l2*m0*p0*r2*x*x + l0*m2*p0*r2*x*x + l0*m0*p2*r2*x*x - 2*l0*l2*q0*r2*x*x - 
  l0*l0*q2*r2*x*x - l2*m2*o0*s0*x*x - l2*m0*o2*s0*x*x - l0*m2*o2*s0*x*x - 
  k2*m2*p0*s0*x*x - k2*m0*p2*s0*x*x - k0*m2*p2*s0*x*x + 2*k2*l2*q0*s0*x*x + 
  2*k2*l0*q2*s0*x*x + 2*k0*l2*q2*s0*x*x + j2*m2*s0*s0*x*x - i2*q2*s0*s0*x*x - 
  l2*m0*o0*s2*x*x - l0*m2*o0*s2*x*x - l0*m0*o2*s2*x*x - k2*m0*p0*s2*x*x - 
  k0*m2*p0*s2*x*x - k0*m0*p2*s2*x*x + 2*k2*l0*q0*s2*x*x + 2*k0*l2*q0*s2*x*x + 
  2*k0*l0*q2*s2*x*x + 2*j2*m0*s0*s2*x*x + 2*j0*m2*s0*s2*x*x - 
  2*i2*q0*s0*s2*x*x - 2*i0*q2*s0*s2*x*x + j0*m0*s2*s2*x*x - i0*q0*s2*s2*x*x + 
  l2*l2*o0*t0*x*x + 2*l0*l2*o2*t0*x*x - k2*l2*p0*t0*x*x - k2*l0*p2*t0*x*x - 
  k0*l2*p2*t0*x*x - j2*l2*s0*t0*x*x + i2*p2*s0*t0*x*x - j2*l0*s2*t0*x*x - 
  j0*l2*s2*t0*x*x + i2*p0*s2*t0*x*x + i0*p2*s2*t0*x*x + 2*l0*l2*o0*t2*x*x + 
  l0*l0*o2*t2*x*x - k2*l0*p0*t2*x*x - k0*l2*p0*t2*x*x - k0*l0*p2*t2*x*x - 
  j2*l0*s0*t2*x*x - j0*l2*s0*t2*x*x + i2*p0*s0*t2*x*x + i0*p2*s0*t2*x*x - 
  j0*l0*s2*t2*x*x + i0*p0*s2*t2*x*x + k2*m2*o0*u0*x*x + k2*m0*o2*u0*x*x + 
  k0*m2*o2*u0*x*x - k2*k2*q0*u0*x*x - 2*k0*k2*q2*u0*x*x - j2*m2*r0*u0*x*x + 
  i2*q2*r0*u0*x*x - j2*m0*r2*u0*x*x - j0*m2*r2*u0*x*x + i2*q0*r2*u0*x*x + 
  i0*q2*r2*u0*x*x + j2*k2*t0*u0*x*x - i2*o2*t0*u0*x*x + j2*k0*t2*u0*x*x + 
  j0*k2*t2*u0*x*x - i2*o0*t2*u0*x*x - i0*o2*t2*u0*x*x + k2*m0*o0*u2*x*x + 
  k0*m2*o0*u2*x*x + k0*m0*o2*u2*x*x - 2*k0*k2*q0*u2*x*x - k0*k0*q2*u2*x*x - 
  j2*m0*r0*u2*x*x - j0*m2*r0*u2*x*x + i2*q0*r0*u2*x*x + i0*q2*r0*u2*x*x - 
  j0*m0*r2*u2*x*x + i0*q0*r2*u2*x*x + j2*k0*t0*u2*x*x + j0*k2*t0*u2*x*x - 
  i2*o0*t0*u2*x*x - i0*o2*t0*u2*x*x + j0*k0*t2*u2*x*x - i0*o0*t2*u2*x*x - 
  k2*l2*o0*v0*x*x - k2*l0*o2*v0*x*x - k0*l2*o2*v0*x*x + k2*k2*p0*v0*x*x + 
  2*k0*k2*p2*v0*x*x + j2*l2*r0*v0*x*x - i2*p2*r0*v0*x*x + j2*l0*r2*v0*x*x + 
  j0*l2*r2*v0*x*x - i2*p0*r2*v0*x*x - i0*p2*r2*v0*x*x - j2*k2*s0*v0*x*x + 
  i2*o2*s0*v0*x*x - j2*k0*s2*v0*x*x - j0*k2*s2*v0*x*x + i2*o0*s2*v0*x*x + 
  i0*o2*s2*v0*x*x - k2*l0*o0*v2*x*x - k0*l2*o0*v2*x*x - k0*l0*o2*v2*x*x + 
  2*k0*k2*p0*v2*x*x + k0*k0*p2*v2*x*x + j2*l0*r0*v2*x*x + j0*l2*r0*v2*x*x - 
  i2*p0*r0*v2*x*x - i0*p2*r0*v2*x*x + j0*l0*r2*v2*x*x - i0*p0*r2*v2*x*x - 
  j2*k0*s0*v2*x*x - j0*k2*s0*v2*x*x + i2*o0*s0*v2*x*x + i0*o2*s0*v2*x*x - 
  j0*k0*s2*v2*x*x + i0*o0*s2*v2*x*x + l2*m0*p0*r0*x*x*x + l0*m2*p0*r0*x*x*x + 
  l0*m0*p2*r0*x*x*x - 2*l0*l2*q0*r0*x*x*x - l0*l0*q2*r0*x*x*x + l0*m0*p0*r2*x*x*x - 
  l0*l0*q0*r2*x*x*x - l2*m0*o0*s0*x*x*x - l0*m2*o0*s0*x*x*x - l0*m0*o2*s0*x*x*x - 
  k2*m0*p0*s0*x*x*x - k0*m2*p0*s0*x*x*x - k0*m0*p2*s0*x*x*x + 2*k2*l0*q0*s0*x*x*x + 
  2*k0*l2*q0*s0*x*x*x + 2*k0*l0*q2*s0*x*x*x + j2*m0*s0*s0*x*x*x + j0*m2*s0*s0*x*x*x - 
  i2*q0*s0*s0*x*x*x - i0*q2*s0*s0*x*x*x - l0*m0*o0*s2*x*x*x - k0*m0*p0*s2*x*x*x + 
  2*k0*l0*q0*s2*x*x*x + 2*j0*m0*s0*s2*x*x*x - 2*i0*q0*s0*s2*x*x*x + 
  2*l0*l2*o0*t0*x*x*x + l0*l0*o2*t0*x*x*x - k2*l0*p0*t0*x*x*x - k0*l2*p0*t0*x*x*x - 
  k0*l0*p2*t0*x*x*x - j2*l0*s0*t0*x*x*x - j0*l2*s0*t0*x*x*x + i2*p0*s0*t0*x*x*x + 
  i0*p2*s0*t0*x*x*x - j0*l0*s2*t0*x*x*x + i0*p0*s2*t0*x*x*x + l0*l0*o0*t2*x*x*x - 
  k0*l0*p0*t2*x*x*x - j0*l0*s0*t2*x*x*x + i0*p0*s0*t2*x*x*x + k2*m0*o0*u0*x*x*x + 
  k0*m2*o0*u0*x*x*x + k0*m0*o2*u0*x*x*x - 2*k0*k2*q0*u0*x*x*x - k0*k0*q2*u0*x*x*x - 
  j2*m0*r0*u0*x*x*x - j0*m2*r0*u0*x*x*x + i2*q0*r0*u0*x*x*x + i0*q2*r0*u0*x*x*x - 
  j0*m0*r2*u0*x*x*x + i0*q0*r2*u0*x*x*x + j2*k0*t0*u0*x*x*x + j0*k2*t0*u0*x*x*x - 
  i2*o0*t0*u0*x*x*x - i0*o2*t0*u0*x*x*x + j0*k0*t2*u0*x*x*x - i0*o0*t2*u0*x*x*x + 
  k0*m0*o0*u2*x*x*x - k0*k0*q0*u2*x*x*x - j0*m0*r0*u2*x*x*x + i0*q0*r0*u2*x*x*x + 
  j0*k0*t0*u2*x*x*x - i0*o0*t0*u2*x*x*x - k2*l0*o0*v0*x*x*x - k0*l2*o0*v0*x*x*x - 
  k0*l0*o2*v0*x*x*x + 2*k0*k2*p0*v0*x*x*x + k0*k0*p2*v0*x*x*x + j2*l0*r0*v0*x*x*x + 
  j0*l2*r0*v0*x*x*x - i2*p0*r0*v0*x*x*x - i0*p2*r0*v0*x*x*x + j0*l0*r2*v0*x*x*x - 
  i0*p0*r2*v0*x*x*x - j2*k0*s0*v0*x*x*x - j0*k2*s0*v0*x*x*x + i2*o0*s0*v0*x*x*x + 
  i0*o2*s0*v0*x*x*x - j0*k0*s2*v0*x*x*x + i0*o0*s2*v0*x*x*x - k0*l0*o0*v2*x*x*x + 
  k0*k0*p0*v2*x*x*x + j0*l0*r0*v2*x*x*x - i0*p0*r0*v2*x*x*x - j0*k0*s0*v2*x*x*x + 
  i0*o0*s0*v2*x*x*x + l0*m0*p0*r0*x*x*x*x - l0*l0*q0*r0*x*x*x*x - l0*m0*o0*s0*x*x*x*x - 
  k0*m0*p0*s0*x*x*x*x + 2*k0*l0*q0*s0*x*x*x*x + j0*m0*s0*s0*x*x*x*x - i0*q0*s0*s0*x*x*x*x + 
  l0*l0*o0*t0*x*x*x*x - k0*l0*p0*t0*x*x*x*x - j0*l0*s0*t0*x*x*x*x + i0*p0*s0*t0*x*x*x*x + 
  k0*m0*o0*u0*x*x*x*x - k0*k0*q0*u0*x*x*x*x - j0*m0*r0*u0*x*x*x*x + i0*q0*r0*u0*x*x*x*x + 
  j0*k0*t0*u0*x*x*x*x - i0*o0*t0*u0*x*x*x*x - k0*l0*o0*v0*x*x*x*x + k0*k0*p0*v0*x*x*x*x + 
  j0*l0*r0*v0*x*x*x*x - i0*p0*r0*v0*x*x*x*x - j0*k0*s0*v0*x*x*x*x + i0*o0*s0*v0*x*x*x*x + 
  l2*m2*p2*r1*y - l2*l2*q2*r1*y + l2*m2*p1*r2*y + l2*m1*p2*r2*y + 
  l1*m2*p2*r2*y - l2*l2*q1*r2*y - 2*l1*l2*q2*r2*y - l2*m2*o2*s1*y - 
  k2*m2*p2*s1*y + 2*k2*l2*q2*s1*y - l2*m2*o1*s2*y - l2*m1*o2*s2*y - 
  l1*m2*o2*s2*y - k2*m2*p1*s2*y - k2*m1*p2*s2*y - k1*m2*p2*s2*y + 
  2*k2*l2*q1*s2*y + 2*k2*l1*q2*s2*y + 2*k1*l2*q2*s2*y + 2*j2*m2*s1*s2*y - 
  2*i2*q2*s1*s2*y + j2*m1*s2*s2*y + j1*m2*s2*s2*y - i2*q1*s2*s2*y - 
  i1*q2*s2*s2*y + l2*l2*o2*t1*y - k2*l2*p2*t1*y - j2*l2*s2*t1*y + 
  i2*p2*s2*t1*y + l2*l2*o1*t2*y + 2*l1*l2*o2*t2*y - k2*l2*p1*t2*y - 
  k2*l1*p2*t2*y - k1*l2*p2*t2*y - j2*l2*s1*t2*y + i2*p2*s1*t2*y - 
  j2*l1*s2*t2*y - j1*l2*s2*t2*y + i2*p1*s2*t2*y + i1*p2*s2*t2*y + 
  k2*m2*o2*u1*y - k2*k2*q2*u1*y - j2*m2*r2*u1*y + i2*q2*r2*u1*y + 
  j2*k2*t2*u1*y - i2*o2*t2*u1*y + k2*m2*o1*u2*y + k2*m1*o2*u2*y + 
  k1*m2*o2*u2*y - k2*k2*q1*u2*y - 2*k1*k2*q2*u2*y - j2*m2*r1*u2*y + 
  i2*q2*r1*u2*y - j2*m1*r2*u2*y - j1*m2*r2*u2*y + i2*q1*r2*u2*y + 
  i1*q2*r2*u2*y + j2*k2*t1*u2*y - i2*o2*t1*u2*y + j2*k1*t2*u2*y + 
  j1*k2*t2*u2*y - i2*o1*t2*u2*y - i1*o2*t2*u2*y - k2*l2*o2*v1*y + 
  k2*k2*p2*v1*y + j2*l2*r2*v1*y - i2*p2*r2*v1*y - j2*k2*s2*v1*y + 
  i2*o2*s2*v1*y - k2*l2*o1*v2*y - k2*l1*o2*v2*y - k1*l2*o2*v2*y + 
  k2*k2*p1*v2*y + 2*k1*k2*p2*v2*y + j2*l2*r1*v2*y - i2*p2*r1*v2*y + 
  j2*l1*r2*v2*y + j1*l2*r2*v2*y - i2*p1*r2*v2*y - i1*p2*r2*v2*y - 
  j2*k2*s1*v2*y + i2*o2*s1*v2*y - j2*k1*s2*v2*y - j1*k2*s2*v2*y + 
  i2*o1*s2*v2*y + i1*o2*s2*v2*y + l2*m2*p1*r0*x*y + l2*m1*p2*r0*x*y + 
  l1*m2*p2*r0*x*y - l2*l2*q1*r0*x*y - 2*l1*l2*q2*r0*x*y + l2*m2*p0*r1*x*y + 
  l2*m0*p2*r1*x*y + l0*m2*p2*r1*x*y - l2*l2*q0*r1*x*y - 2*l0*l2*q2*r1*x*y + 
  l2*m1*p0*r2*x*y + l1*m2*p0*r2*x*y + l2*m0*p1*r2*x*y + l0*m2*p1*r2*x*y + 
  l1*m0*p2*r2*x*y + l0*m1*p2*r2*x*y - 2*l1*l2*q0*r2*x*y - 2*l0*l2*q1*r2*x*y - 
  2*l0*l1*q2*r2*x*y - l2*m2*o1*s0*x*y - l2*m1*o2*s0*x*y - l1*m2*o2*s0*x*y - 
  k2*m2*p1*s0*x*y - k2*m1*p2*s0*x*y - k1*m2*p2*s0*x*y + 2*k2*l2*q1*s0*x*y + 
  2*k2*l1*q2*s0*x*y + 2*k1*l2*q2*s0*x*y - l2*m2*o0*s1*x*y - l2*m0*o2*s1*x*y - 
  l0*m2*o2*s1*x*y - k2*m2*p0*s1*x*y - k2*m0*p2*s1*x*y - k0*m2*p2*s1*x*y + 
  2*k2*l2*q0*s1*x*y + 2*k2*l0*q2*s1*x*y + 2*k0*l2*q2*s1*x*y + 
  2*j2*m2*s0*s1*x*y - 2*i2*q2*s0*s1*x*y - l2*m1*o0*s2*x*y - l1*m2*o0*s2*x*y - 
  l2*m0*o1*s2*x*y - l0*m2*o1*s2*x*y - l1*m0*o2*s2*x*y - l0*m1*o2*s2*x*y - 
  k2*m1*p0*s2*x*y - k1*m2*p0*s2*x*y - k2*m0*p1*s2*x*y - k0*m2*p1*s2*x*y - 
  k1*m0*p2*s2*x*y - k0*m1*p2*s2*x*y + 2*k2*l1*q0*s2*x*y + 2*k1*l2*q0*s2*x*y + 
  2*k2*l0*q1*s2*x*y + 2*k0*l2*q1*s2*x*y + 2*k1*l0*q2*s2*x*y + 
  2*k0*l1*q2*s2*x*y + 2*j2*m1*s0*s2*x*y + 2*j1*m2*s0*s2*x*y - 
  2*i2*q1*s0*s2*x*y - 2*i1*q2*s0*s2*x*y + 2*j2*m0*s1*s2*x*y + 
  2*j0*m2*s1*s2*x*y - 2*i2*q0*s1*s2*x*y - 2*i0*q2*s1*s2*x*y + 
  j1*m0*s2*s2*x*y + j0*m1*s2*s2*x*y - i1*q0*s2*s2*x*y - i0*q1*s2*s2*x*y + 
  l2*l2*o1*t0*x*y + 2*l1*l2*o2*t0*x*y - k2*l2*p1*t0*x*y - k2*l1*p2*t0*x*y - 
  k1*l2*p2*t0*x*y - j2*l2*s1*t0*x*y + i2*p2*s1*t0*x*y - j2*l1*s2*t0*x*y - 
  j1*l2*s2*t0*x*y + i2*p1*s2*t0*x*y + i1*p2*s2*t0*x*y + l2*l2*o0*t1*x*y + 
  2*l0*l2*o2*t1*x*y - k2*l2*p0*t1*x*y - k2*l0*p2*t1*x*y - k0*l2*p2*t1*x*y - 
  j2*l2*s0*t1*x*y + i2*p2*s0*t1*x*y - j2*l0*s2*t1*x*y - j0*l2*s2*t1*x*y + 
  i2*p0*s2*t1*x*y + i0*p2*s2*t1*x*y + 2*l1*l2*o0*t2*x*y + 2*l0*l2*o1*t2*x*y + 
  2*l0*l1*o2*t2*x*y - k2*l1*p0*t2*x*y - k1*l2*p0*t2*x*y - k2*l0*p1*t2*x*y - 
  k0*l2*p1*t2*x*y - k1*l0*p2*t2*x*y - k0*l1*p2*t2*x*y - j2*l1*s0*t2*x*y - 
  j1*l2*s0*t2*x*y + i2*p1*s0*t2*x*y + i1*p2*s0*t2*x*y - j2*l0*s1*t2*x*y - 
  j0*l2*s1*t2*x*y + i2*p0*s1*t2*x*y + i0*p2*s1*t2*x*y - j1*l0*s2*t2*x*y - 
  j0*l1*s2*t2*x*y + i1*p0*s2*t2*x*y + i0*p1*s2*t2*x*y + k2*m2*o1*u0*x*y + 
  k2*m1*o2*u0*x*y + k1*m2*o2*u0*x*y - k2*k2*q1*u0*x*y - 2*k1*k2*q2*u0*x*y - 
  j2*m2*r1*u0*x*y + i2*q2*r1*u0*x*y - j2*m1*r2*u0*x*y - j1*m2*r2*u0*x*y + 
  i2*q1*r2*u0*x*y + i1*q2*r2*u0*x*y + j2*k2*t1*u0*x*y - i2*o2*t1*u0*x*y + 
  j2*k1*t2*u0*x*y + j1*k2*t2*u0*x*y - i2*o1*t2*u0*x*y - i1*o2*t2*u0*x*y + 
  k2*m2*o0*u1*x*y + k2*m0*o2*u1*x*y + k0*m2*o2*u1*x*y - k2*k2*q0*u1*x*y - 
  2*k0*k2*q2*u1*x*y - j2*m2*r0*u1*x*y + i2*q2*r0*u1*x*y - j2*m0*r2*u1*x*y - 
  j0*m2*r2*u1*x*y + i2*q0*r2*u1*x*y + i0*q2*r2*u1*x*y + j2*k2*t0*u1*x*y - 
  i2*o2*t0*u1*x*y + j2*k0*t2*u1*x*y + j0*k2*t2*u1*x*y - i2*o0*t2*u1*x*y - 
  i0*o2*t2*u1*x*y + k2*m1*o0*u2*x*y + k1*m2*o0*u2*x*y + k2*m0*o1*u2*x*y + 
  k0*m2*o1*u2*x*y + k1*m0*o2*u2*x*y + k0*m1*o2*u2*x*y - 2*k1*k2*q0*u2*x*y - 
  2*k0*k2*q1*u2*x*y - 2*k0*k1*q2*u2*x*y - j2*m1*r0*u2*x*y - j1*m2*r0*u2*x*y + 
  i2*q1*r0*u2*x*y + i1*q2*r0*u2*x*y - j2*m0*r1*u2*x*y - j0*m2*r1*u2*x*y + 
  i2*q0*r1*u2*x*y + i0*q2*r1*u2*x*y - j1*m0*r2*u2*x*y - j0*m1*r2*u2*x*y + 
  i1*q0*r2*u2*x*y + i0*q1*r2*u2*x*y + j2*k1*t0*u2*x*y + j1*k2*t0*u2*x*y - 
  i2*o1*t0*u2*x*y - i1*o2*t0*u2*x*y + j2*k0*t1*u2*x*y + j0*k2*t1*u2*x*y - 
  i2*o0*t1*u2*x*y - i0*o2*t1*u2*x*y + j1*k0*t2*u2*x*y + j0*k1*t2*u2*x*y - 
  i1*o0*t2*u2*x*y - i0*o1*t2*u2*x*y - k2*l2*o1*v0*x*y - k2*l1*o2*v0*x*y - 
  k1*l2*o2*v0*x*y + k2*k2*p1*v0*x*y + 2*k1*k2*p2*v0*x*y + j2*l2*r1*v0*x*y - 
  i2*p2*r1*v0*x*y + j2*l1*r2*v0*x*y + j1*l2*r2*v0*x*y - i2*p1*r2*v0*x*y - 
  i1*p2*r2*v0*x*y - j2*k2*s1*v0*x*y + i2*o2*s1*v0*x*y - j2*k1*s2*v0*x*y - 
  j1*k2*s2*v0*x*y + i2*o1*s2*v0*x*y + i1*o2*s2*v0*x*y - k2*l2*o0*v1*x*y - 
  k2*l0*o2*v1*x*y - k0*l2*o2*v1*x*y + k2*k2*p0*v1*x*y + 2*k0*k2*p2*v1*x*y + 
  j2*l2*r0*v1*x*y - i2*p2*r0*v1*x*y + j2*l0*r2*v1*x*y + j0*l2*r2*v1*x*y - 
  i2*p0*r2*v1*x*y - i0*p2*r2*v1*x*y - j2*k2*s0*v1*x*y + i2*o2*s0*v1*x*y - 
  j2*k0*s2*v1*x*y - j0*k2*s2*v1*x*y + i2*o0*s2*v1*x*y + i0*o2*s2*v1*x*y - 
  k2*l1*o0*v2*x*y - k1*l2*o0*v2*x*y - k2*l0*o1*v2*x*y - k0*l2*o1*v2*x*y - 
  k1*l0*o2*v2*x*y - k0*l1*o2*v2*x*y + 2*k1*k2*p0*v2*x*y + 2*k0*k2*p1*v2*x*y + 
  2*k0*k1*p2*v2*x*y + j2*l1*r0*v2*x*y + j1*l2*r0*v2*x*y - i2*p1*r0*v2*x*y - 
  i1*p2*r0*v2*x*y + j2*l0*r1*v2*x*y + j0*l2*r1*v2*x*y - i2*p0*r1*v2*x*y - 
  i0*p2*r1*v2*x*y + j1*l0*r2*v2*x*y + j0*l1*r2*v2*x*y - i1*p0*r2*v2*x*y - 
  i0*p1*r2*v2*x*y - j2*k1*s0*v2*x*y - j1*k2*s0*v2*x*y + i2*o1*s0*v2*x*y + 
  i1*o2*s0*v2*x*y - j2*k0*s1*v2*x*y - j0*k2*s1*v2*x*y + i2*o0*s1*v2*x*y + 
  i0*o2*s1*v2*x*y - j1*k0*s2*v2*x*y - j0*k1*s2*v2*x*y + i1*o0*s2*v2*x*y + 
  i0*o1*s2*v2*x*y + l2*m1*p0*r0*x*x*y + l1*m2*p0*r0*x*x*y + 
  l2*m0*p1*r0*x*x*y + l0*m2*p1*r0*x*x*y + l1*m0*p2*r0*x*x*y + 
  l0*m1*p2*r0*x*x*y - 2*l1*l2*q0*r0*x*x*y - 2*l0*l2*q1*r0*x*x*y - 
  2*l0*l1*q2*r0*x*x*y + l2*m0*p0*r1*x*x*y + l0*m2*p0*r1*x*x*y + 
  l0*m0*p2*r1*x*x*y - 2*l0*l2*q0*r1*x*x*y - l0*l0*q2*r1*x*x*y + 
  l1*m0*p0*r2*x*x*y + l0*m1*p0*r2*x*x*y + l0*m0*p1*r2*x*x*y - 
  2*l0*l1*q0*r2*x*x*y - l0*l0*q1*r2*x*x*y - l2*m1*o0*s0*x*x*y - 
  l1*m2*o0*s0*x*x*y - l2*m0*o1*s0*x*x*y - l0*m2*o1*s0*x*x*y - 
  l1*m0*o2*s0*x*x*y - l0*m1*o2*s0*x*x*y - k2*m1*p0*s0*x*x*y - 
  k1*m2*p0*s0*x*x*y - k2*m0*p1*s0*x*x*y - k0*m2*p1*s0*x*x*y - 
  k1*m0*p2*s0*x*x*y - k0*m1*p2*s0*x*x*y + 2*k2*l1*q0*s0*x*x*y + 
  2*k1*l2*q0*s0*x*x*y + 2*k2*l0*q1*s0*x*x*y + 2*k0*l2*q1*s0*x*x*y + 
  2*k1*l0*q2*s0*x*x*y + 2*k0*l1*q2*s0*x*x*y + j2*m1*s0*s0*x*x*y + 
  j1*m2*s0*s0*x*x*y - i2*q1*s0*s0*x*x*y - i1*q2*s0*s0*x*x*y - 
  l2*m0*o0*s1*x*x*y - l0*m2*o0*s1*x*x*y - l0*m0*o2*s1*x*x*y - 
  k2*m0*p0*s1*x*x*y - k0*m2*p0*s1*x*x*y - k0*m0*p2*s1*x*x*y + 
  2*k2*l0*q0*s1*x*x*y + 2*k0*l2*q0*s1*x*x*y + 2*k0*l0*q2*s1*x*x*y + 
  2*j2*m0*s0*s1*x*x*y + 2*j0*m2*s0*s1*x*x*y - 2*i2*q0*s0*s1*x*x*y - 
  2*i0*q2*s0*s1*x*x*y - l1*m0*o0*s2*x*x*y - l0*m1*o0*s2*x*x*y - 
  l0*m0*o1*s2*x*x*y - k1*m0*p0*s2*x*x*y - k0*m1*p0*s2*x*x*y - 
  k0*m0*p1*s2*x*x*y + 2*k1*l0*q0*s2*x*x*y + 2*k0*l1*q0*s2*x*x*y + 
  2*k0*l0*q1*s2*x*x*y + 2*j1*m0*s0*s2*x*x*y + 2*j0*m1*s0*s2*x*x*y - 
  2*i1*q0*s0*s2*x*x*y - 2*i0*q1*s0*s2*x*x*y + 2*j0*m0*s1*s2*x*x*y - 
  2*i0*q0*s1*s2*x*x*y + 2*l1*l2*o0*t0*x*x*y + 2*l0*l2*o1*t0*x*x*y + 
  2*l0*l1*o2*t0*x*x*y - k2*l1*p0*t0*x*x*y - k1*l2*p0*t0*x*x*y - 
  k2*l0*p1*t0*x*x*y - k0*l2*p1*t0*x*x*y - k1*l0*p2*t0*x*x*y - 
  k0*l1*p2*t0*x*x*y - j2*l1*s0*t0*x*x*y - j1*l2*s0*t0*x*x*y + 
  i2*p1*s0*t0*x*x*y + i1*p2*s0*t0*x*x*y - j2*l0*s1*t0*x*x*y - 
  j0*l2*s1*t0*x*x*y + i2*p0*s1*t0*x*x*y + i0*p2*s1*t0*x*x*y - 
  j1*l0*s2*t0*x*x*y - j0*l1*s2*t0*x*x*y + i1*p0*s2*t0*x*x*y + 
  i0*p1*s2*t0*x*x*y + 2*l0*l2*o0*t1*x*x*y + l0*l0*o2*t1*x*x*y - 
  k2*l0*p0*t1*x*x*y - k0*l2*p0*t1*x*x*y - k0*l0*p2*t1*x*x*y - 
  j2*l0*s0*t1*x*x*y - j0*l2*s0*t1*x*x*y + i2*p0*s0*t1*x*x*y + 
  i0*p2*s0*t1*x*x*y - j0*l0*s2*t1*x*x*y + i0*p0*s2*t1*x*x*y + 
  2*l0*l1*o0*t2*x*x*y + l0*l0*o1*t2*x*x*y - k1*l0*p0*t2*x*x*y - 
  k0*l1*p0*t2*x*x*y - k0*l0*p1*t2*x*x*y - j1*l0*s0*t2*x*x*y - 
  j0*l1*s0*t2*x*x*y + i1*p0*s0*t2*x*x*y + i0*p1*s0*t2*x*x*y - 
  j0*l0*s1*t2*x*x*y + i0*p0*s1*t2*x*x*y + k2*m1*o0*u0*x*x*y + 
  k1*m2*o0*u0*x*x*y + k2*m0*o1*u0*x*x*y + k0*m2*o1*u0*x*x*y + 
  k1*m0*o2*u0*x*x*y + k0*m1*o2*u0*x*x*y - 2*k1*k2*q0*u0*x*x*y - 
  2*k0*k2*q1*u0*x*x*y - 2*k0*k1*q2*u0*x*x*y - j2*m1*r0*u0*x*x*y - 
  j1*m2*r0*u0*x*x*y + i2*q1*r0*u0*x*x*y + i1*q2*r0*u0*x*x*y - 
  j2*m0*r1*u0*x*x*y - j0*m2*r1*u0*x*x*y + i2*q0*r1*u0*x*x*y + 
  i0*q2*r1*u0*x*x*y - j1*m0*r2*u0*x*x*y - j0*m1*r2*u0*x*x*y + 
  i1*q0*r2*u0*x*x*y + i0*q1*r2*u0*x*x*y + j2*k1*t0*u0*x*x*y + 
  j1*k2*t0*u0*x*x*y - i2*o1*t0*u0*x*x*y - i1*o2*t0*u0*x*x*y + 
  j2*k0*t1*u0*x*x*y + j0*k2*t1*u0*x*x*y - i2*o0*t1*u0*x*x*y - 
  i0*o2*t1*u0*x*x*y + j1*k0*t2*u0*x*x*y + j0*k1*t2*u0*x*x*y - 
  i1*o0*t2*u0*x*x*y - i0*o1*t2*u0*x*x*y + k2*m0*o0*u1*x*x*y + 
  k0*m2*o0*u1*x*x*y + k0*m0*o2*u1*x*x*y - 2*k0*k2*q0*u1*x*x*y - 
  k0*k0*q2*u1*x*x*y - j2*m0*r0*u1*x*x*y - j0*m2*r0*u1*x*x*y + 
  i2*q0*r0*u1*x*x*y + i0*q2*r0*u1*x*x*y - j0*m0*r2*u1*x*x*y + 
  i0*q0*r2*u1*x*x*y + j2*k0*t0*u1*x*x*y + j0*k2*t0*u1*x*x*y - 
  i2*o0*t0*u1*x*x*y - i0*o2*t0*u1*x*x*y + j0*k0*t2*u1*x*x*y - 
  i0*o0*t2*u1*x*x*y + k1*m0*o0*u2*x*x*y + k0*m1*o0*u2*x*x*y + 
  k0*m0*o1*u2*x*x*y - 2*k0*k1*q0*u2*x*x*y - k0*k0*q1*u2*x*x*y - 
  j1*m0*r0*u2*x*x*y - j0*m1*r0*u2*x*x*y + i1*q0*r0*u2*x*x*y + 
  i0*q1*r0*u2*x*x*y - j0*m0*r1*u2*x*x*y + i0*q0*r1*u2*x*x*y + 
  j1*k0*t0*u2*x*x*y + j0*k1*t0*u2*x*x*y - i1*o0*t0*u2*x*x*y - 
  i0*o1*t0*u2*x*x*y + j0*k0*t1*u2*x*x*y - i0*o0*t1*u2*x*x*y - 
  k2*l1*o0*v0*x*x*y - k1*l2*o0*v0*x*x*y - k2*l0*o1*v0*x*x*y - 
  k0*l2*o1*v0*x*x*y - k1*l0*o2*v0*x*x*y - k0*l1*o2*v0*x*x*y + 
  2*k1*k2*p0*v0*x*x*y + 2*k0*k2*p1*v0*x*x*y + 2*k0*k1*p2*v0*x*x*y + 
  j2*l1*r0*v0*x*x*y + j1*l2*r0*v0*x*x*y - i2*p1*r0*v0*x*x*y - 
  i1*p2*r0*v0*x*x*y + j2*l0*r1*v0*x*x*y + j0*l2*r1*v0*x*x*y - 
  i2*p0*r1*v0*x*x*y - i0*p2*r1*v0*x*x*y + j1*l0*r2*v0*x*x*y + 
  j0*l1*r2*v0*x*x*y - i1*p0*r2*v0*x*x*y - i0*p1*r2*v0*x*x*y - 
  j2*k1*s0*v0*x*x*y - j1*k2*s0*v0*x*x*y + i2*o1*s0*v0*x*x*y + 
  i1*o2*s0*v0*x*x*y - j2*k0*s1*v0*x*x*y - j0*k2*s1*v0*x*x*y + 
  i2*o0*s1*v0*x*x*y + i0*o2*s1*v0*x*x*y - j1*k0*s2*v0*x*x*y - 
  j0*k1*s2*v0*x*x*y + i1*o0*s2*v0*x*x*y + i0*o1*s2*v0*x*x*y - 
  k2*l0*o0*v1*x*x*y - k0*l2*o0*v1*x*x*y - k0*l0*o2*v1*x*x*y + 
  2*k0*k2*p0*v1*x*x*y + k0*k0*p2*v1*x*x*y + j2*l0*r0*v1*x*x*y + 
  j0*l2*r0*v1*x*x*y - i2*p0*r0*v1*x*x*y - i0*p2*r0*v1*x*x*y + 
  j0*l0*r2*v1*x*x*y - i0*p0*r2*v1*x*x*y - j2*k0*s0*v1*x*x*y - 
  j0*k2*s0*v1*x*x*y + i2*o0*s0*v1*x*x*y + i0*o2*s0*v1*x*x*y - 
  j0*k0*s2*v1*x*x*y + i0*o0*s2*v1*x*x*y - k1*l0*o0*v2*x*x*y - 
  k0*l1*o0*v2*x*x*y - k0*l0*o1*v2*x*x*y + 2*k0*k1*p0*v2*x*x*y + 
  k0*k0*p1*v2*x*x*y + j1*l0*r0*v2*x*x*y + j0*l1*r0*v2*x*x*y - 
  i1*p0*r0*v2*x*x*y - i0*p1*r0*v2*x*x*y + j0*l0*r1*v2*x*x*y - 
  i0*p0*r1*v2*x*x*y - j1*k0*s0*v2*x*x*y - j0*k1*s0*v2*x*x*y + 
  i1*o0*s0*v2*x*x*y + i0*o1*s0*v2*x*x*y - j0*k0*s1*v2*x*x*y + 
  i0*o0*s1*v2*x*x*y + l1*m0*p0*r0*x*x*x*y + l0*m1*p0*r0*x*x*x*y + 
  l0*m0*p1*r0*x*x*x*y - 2*l0*l1*q0*r0*x*x*x*y - l0*l0*q1*r0*x*x*x*y + 
  l0*m0*p0*r1*x*x*x*y - l0*l0*q0*r1*x*x*x*y - l1*m0*o0*s0*x*x*x*y - 
  l0*m1*o0*s0*x*x*x*y - l0*m0*o1*s0*x*x*x*y - k1*m0*p0*s0*x*x*x*y - 
  k0*m1*p0*s0*x*x*x*y - k0*m0*p1*s0*x*x*x*y + 2*k1*l0*q0*s0*x*x*x*y + 
  2*k0*l1*q0*s0*x*x*x*y + 2*k0*l0*q1*s0*x*x*x*y + j1*m0*s0*s0*x*x*x*y + 
  j0*m1*s0*s0*x*x*x*y - i1*q0*s0*s0*x*x*x*y - i0*q1*s0*s0*x*x*x*y - 
  l0*m0*o0*s1*x*x*x*y - k0*m0*p0*s1*x*x*x*y + 2*k0*l0*q0*s1*x*x*x*y + 
  2*j0*m0*s0*s1*x*x*x*y - 2*i0*q0*s0*s1*x*x*x*y + 2*l0*l1*o0*t0*x*x*x*y + 
  l0*l0*o1*t0*x*x*x*y - k1*l0*p0*t0*x*x*x*y - k0*l1*p0*t0*x*x*x*y - 
  k0*l0*p1*t0*x*x*x*y - j1*l0*s0*t0*x*x*x*y - j0*l1*s0*t0*x*x*x*y + 
  i1*p0*s0*t0*x*x*x*y + i0*p1*s0*t0*x*x*x*y - j0*l0*s1*t0*x*x*x*y + 
  i0*p0*s1*t0*x*x*x*y + l0*l0*o0*t1*x*x*x*y - k0*l0*p0*t1*x*x*x*y - 
  j0*l0*s0*t1*x*x*x*y + i0*p0*s0*t1*x*x*x*y + k1*m0*o0*u0*x*x*x*y + 
  k0*m1*o0*u0*x*x*x*y + k0*m0*o1*u0*x*x*x*y - 2*k0*k1*q0*u0*x*x*x*y - 
  k0*k0*q1*u0*x*x*x*y - j1*m0*r0*u0*x*x*x*y - j0*m1*r0*u0*x*x*x*y + 
  i1*q0*r0*u0*x*x*x*y + i0*q1*r0*u0*x*x*x*y - j0*m0*r1*u0*x*x*x*y + 
  i0*q0*r1*u0*x*x*x*y + j1*k0*t0*u0*x*x*x*y + j0*k1*t0*u0*x*x*x*y - 
  i1*o0*t0*u0*x*x*x*y - i0*o1*t0*u0*x*x*x*y + j0*k0*t1*u0*x*x*x*y - 
  i0*o0*t1*u0*x*x*x*y + k0*m0*o0*u1*x*x*x*y - k0*k0*q0*u1*x*x*x*y - 
  j0*m0*r0*u1*x*x*x*y + i0*q0*r0*u1*x*x*x*y + j0*k0*t0*u1*x*x*x*y - 
  i0*o0*t0*u1*x*x*x*y - k1*l0*o0*v0*x*x*x*y - k0*l1*o0*v0*x*x*x*y - 
  k0*l0*o1*v0*x*x*x*y + 2*k0*k1*p0*v0*x*x*x*y + k0*k0*p1*v0*x*x*x*y + 
  j1*l0*r0*v0*x*x*x*y + j0*l1*r0*v0*x*x*x*y - i1*p0*r0*v0*x*x*x*y - 
  i0*p1*r0*v0*x*x*x*y + j0*l0*r1*v0*x*x*x*y - i0*p0*r1*v0*x*x*x*y - 
  j1*k0*s0*v0*x*x*x*y - j0*k1*s0*v0*x*x*x*y + i1*o0*s0*v0*x*x*x*y + 
  i0*o1*s0*v0*x*x*x*y - j0*k0*s1*v0*x*x*x*y + i0*o0*s1*v0*x*x*x*y - 
  k0*l0*o0*v1*x*x*x*y + k0*k0*p0*v1*x*x*x*y + j0*l0*r0*v1*x*x*x*y - 
  i0*p0*r0*v1*x*x*x*y - j0*k0*s0*v1*x*x*x*y + i0*o0*s0*v1*x*x*x*y + 
  l2*m2*p1*r1*y*y + l2*m1*p2*r1*y*y + l1*m2*p2*r1*y*y - l2*l2*q1*r1*y*y - 
  2*l1*l2*q2*r1*y*y + l2*m1*p1*r2*y*y + l1*m2*p1*r2*y*y + l1*m1*p2*r2*y*y - 
  2*l1*l2*q1*r2*y*y - l1*l1*q2*r2*y*y - l2*m2*o1*s1*y*y - l2*m1*o2*s1*y*y - 
  l1*m2*o2*s1*y*y - k2*m2*p1*s1*y*y - k2*m1*p2*s1*y*y - k1*m2*p2*s1*y*y + 
  2*k2*l2*q1*s1*y*y + 2*k2*l1*q2*s1*y*y + 2*k1*l2*q2*s1*y*y + 
  j2*m2*s1*s1*y*y - i2*q2*s1*s1*y*y - l2*m1*o1*s2*y*y - l1*m2*o1*s2*y*y - 
  l1*m1*o2*s2*y*y - k2*m1*p1*s2*y*y - k1*m2*p1*s2*y*y - k1*m1*p2*s2*y*y + 
  2*k2*l1*q1*s2*y*y + 2*k1*l2*q1*s2*y*y + 2*k1*l1*q2*s2*y*y + 
  2*j2*m1*s1*s2*y*y + 2*j1*m2*s1*s2*y*y - 2*i2*q1*s1*s2*y*y - 
  2*i1*q2*s1*s2*y*y + j1*m1*s2*s2*y*y - i1*q1*s2*s2*y*y + l2*l2*o1*t1*y*y + 
  2*l1*l2*o2*t1*y*y - k2*l2*p1*t1*y*y - k2*l1*p2*t1*y*y - k1*l2*p2*t1*y*y - 
  j2*l2*s1*t1*y*y + i2*p2*s1*t1*y*y - j2*l1*s2*t1*y*y - j1*l2*s2*t1*y*y + 
  i2*p1*s2*t1*y*y + i1*p2*s2*t1*y*y + 2*l1*l2*o1*t2*y*y + l1*l1*o2*t2*y*y - 
  k2*l1*p1*t2*y*y - k1*l2*p1*t2*y*y - k1*l1*p2*t2*y*y - j2*l1*s1*t2*y*y - 
  j1*l2*s1*t2*y*y + i2*p1*s1*t2*y*y + i1*p2*s1*t2*y*y - j1*l1*s2*t2*y*y + 
  i1*p1*s2*t2*y*y + k2*m2*o1*u1*y*y + k2*m1*o2*u1*y*y + k1*m2*o2*u1*y*y - 
  k2*k2*q1*u1*y*y - 2*k1*k2*q2*u1*y*y - j2*m2*r1*u1*y*y + i2*q2*r1*u1*y*y - 
  j2*m1*r2*u1*y*y - j1*m2*r2*u1*y*y + i2*q1*r2*u1*y*y + i1*q2*r2*u1*y*y + 
  j2*k2*t1*u1*y*y - i2*o2*t1*u1*y*y + j2*k1*t2*u1*y*y + j1*k2*t2*u1*y*y - 
  i2*o1*t2*u1*y*y - i1*o2*t2*u1*y*y + k2*m1*o1*u2*y*y + k1*m2*o1*u2*y*y + 
  k1*m1*o2*u2*y*y - 2*k1*k2*q1*u2*y*y - k1*k1*q2*u2*y*y - j2*m1*r1*u2*y*y - 
  j1*m2*r1*u2*y*y + i2*q1*r1*u2*y*y + i1*q2*r1*u2*y*y - j1*m1*r2*u2*y*y + 
  i1*q1*r2*u2*y*y + j2*k1*t1*u2*y*y + j1*k2*t1*u2*y*y - i2*o1*t1*u2*y*y - 
  i1*o2*t1*u2*y*y + j1*k1*t2*u2*y*y - i1*o1*t2*u2*y*y - k2*l2*o1*v1*y*y - 
  k2*l1*o2*v1*y*y - k1*l2*o2*v1*y*y + k2*k2*p1*v1*y*y + 2*k1*k2*p2*v1*y*y + 
  j2*l2*r1*v1*y*y - i2*p2*r1*v1*y*y + j2*l1*r2*v1*y*y + j1*l2*r2*v1*y*y - 
  i2*p1*r2*v1*y*y - i1*p2*r2*v1*y*y - j2*k2*s1*v1*y*y + i2*o2*s1*v1*y*y - 
  j2*k1*s2*v1*y*y - j1*k2*s2*v1*y*y + i2*o1*s2*v1*y*y + i1*o2*s2*v1*y*y - 
  k2*l1*o1*v2*y*y - k1*l2*o1*v2*y*y - k1*l1*o2*v2*y*y + 2*k1*k2*p1*v2*y*y + 
  k1*k1*p2*v2*y*y + j2*l1*r1*v2*y*y + j1*l2*r1*v2*y*y - i2*p1*r1*v2*y*y - 
  i1*p2*r1*v2*y*y + j1*l1*r2*v2*y*y - i1*p1*r2*v2*y*y - j2*k1*s1*v2*y*y - 
  j1*k2*s1*v2*y*y + i2*o1*s1*v2*y*y + i1*o2*s1*v2*y*y - j1*k1*s2*v2*y*y + 
  i1*o1*s2*v2*y*y + l2*m1*p1*r0*x*y*y + l1*m2*p1*r0*x*y*y + 
  l1*m1*p2*r0*x*y*y - 2*l1*l2*q1*r0*x*y*y - l1*l1*q2*r0*x*y*y + 
  l2*m1*p0*r1*x*y*y + l1*m2*p0*r1*x*y*y + l2*m0*p1*r1*x*y*y + 
  l0*m2*p1*r1*x*y*y + l1*m0*p2*r1*x*y*y + l0*m1*p2*r1*x*y*y - 
  2*l1*l2*q0*r1*x*y*y - 2*l0*l2*q1*r1*x*y*y - 2*l0*l1*q2*r1*x*y*y + 
  l1*m1*p0*r2*x*y*y + l1*m0*p1*r2*x*y*y + l0*m1*p1*r2*x*y*y - 
  l1*l1*q0*r2*x*y*y - 2*l0*l1*q1*r2*x*y*y - l2*m1*o1*s0*x*y*y - 
  l1*m2*o1*s0*x*y*y - l1*m1*o2*s0*x*y*y - k2*m1*p1*s0*x*y*y - 
  k1*m2*p1*s0*x*y*y - k1*m1*p2*s0*x*y*y + 2*k2*l1*q1*s0*x*y*y + 
  2*k1*l2*q1*s0*x*y*y + 2*k1*l1*q2*s0*x*y*y - l2*m1*o0*s1*x*y*y - 
  l1*m2*o0*s1*x*y*y - l2*m0*o1*s1*x*y*y - l0*m2*o1*s1*x*y*y - 
  l1*m0*o2*s1*x*y*y - l0*m1*o2*s1*x*y*y - k2*m1*p0*s1*x*y*y - 
  k1*m2*p0*s1*x*y*y - k2*m0*p1*s1*x*y*y - k0*m2*p1*s1*x*y*y - 
  k1*m0*p2*s1*x*y*y - k0*m1*p2*s1*x*y*y + 2*k2*l1*q0*s1*x*y*y + 
  2*k1*l2*q0*s1*x*y*y + 2*k2*l0*q1*s1*x*y*y + 2*k0*l2*q1*s1*x*y*y + 
  2*k1*l0*q2*s1*x*y*y + 2*k0*l1*q2*s1*x*y*y + 2*j2*m1*s0*s1*x*y*y + 
  2*j1*m2*s0*s1*x*y*y - 2*i2*q1*s0*s1*x*y*y - 2*i1*q2*s0*s1*x*y*y + 
  j2*m0*s1*s1*x*y*y + j0*m2*s1*s1*x*y*y - i2*q0*s1*s1*x*y*y - i0*q2*s1*s1*x*y*y - 
  l1*m1*o0*s2*x*y*y - l1*m0*o1*s2*x*y*y - l0*m1*o1*s2*x*y*y - 
  k1*m1*p0*s2*x*y*y - k1*m0*p1*s2*x*y*y - k0*m1*p1*s2*x*y*y + 
  2*k1*l1*q0*s2*x*y*y + 2*k1*l0*q1*s2*x*y*y + 2*k0*l1*q1*s2*x*y*y + 
  2*j1*m1*s0*s2*x*y*y - 2*i1*q1*s0*s2*x*y*y + 2*j1*m0*s1*s2*x*y*y + 
  2*j0*m1*s1*s2*x*y*y - 2*i1*q0*s1*s2*x*y*y - 2*i0*q1*s1*s2*x*y*y + 
  2*l1*l2*o1*t0*x*y*y + l1*l1*o2*t0*x*y*y - k2*l1*p1*t0*x*y*y - 
  k1*l2*p1*t0*x*y*y - k1*l1*p2*t0*x*y*y - j2*l1*s1*t0*x*y*y - 
  j1*l2*s1*t0*x*y*y + i2*p1*s1*t0*x*y*y + i1*p2*s1*t0*x*y*y - 
  j1*l1*s2*t0*x*y*y + i1*p1*s2*t0*x*y*y + 2*l1*l2*o0*t1*x*y*y + 
  2*l0*l2*o1*t1*x*y*y + 2*l0*l1*o2*t1*x*y*y - k2*l1*p0*t1*x*y*y - 
  k1*l2*p0*t1*x*y*y - k2*l0*p1*t1*x*y*y - k0*l2*p1*t1*x*y*y - 
  k1*l0*p2*t1*x*y*y - k0*l1*p2*t1*x*y*y - j2*l1*s0*t1*x*y*y - 
  j1*l2*s0*t1*x*y*y + i2*p1*s0*t1*x*y*y + i1*p2*s0*t1*x*y*y - 
  j2*l0*s1*t1*x*y*y - j0*l2*s1*t1*x*y*y + i2*p0*s1*t1*x*y*y + 
  i0*p2*s1*t1*x*y*y - j1*l0*s2*t1*x*y*y - j0*l1*s2*t1*x*y*y + 
  i1*p0*s2*t1*x*y*y + i0*p1*s2*t1*x*y*y + l1*l1*o0*t2*x*y*y + 
  2*l0*l1*o1*t2*x*y*y - k1*l1*p0*t2*x*y*y - k1*l0*p1*t2*x*y*y - 
  k0*l1*p1*t2*x*y*y - j1*l1*s0*t2*x*y*y + i1*p1*s0*t2*x*y*y - 
  j1*l0*s1*t2*x*y*y - j0*l1*s1*t2*x*y*y + i1*p0*s1*t2*x*y*y + 
  i0*p1*s1*t2*x*y*y + k2*m1*o1*u0*x*y*y + k1*m2*o1*u0*x*y*y + 
  k1*m1*o2*u0*x*y*y - 2*k1*k2*q1*u0*x*y*y - k1*k1*q2*u0*x*y*y - 
  j2*m1*r1*u0*x*y*y - j1*m2*r1*u0*x*y*y + i2*q1*r1*u0*x*y*y + 
  i1*q2*r1*u0*x*y*y - j1*m1*r2*u0*x*y*y + i1*q1*r2*u0*x*y*y + 
  j2*k1*t1*u0*x*y*y + j1*k2*t1*u0*x*y*y - i2*o1*t1*u0*x*y*y - 
  i1*o2*t1*u0*x*y*y + j1*k1*t2*u0*x*y*y - i1*o1*t2*u0*x*y*y + 
  k2*m1*o0*u1*x*y*y + k1*m2*o0*u1*x*y*y + k2*m0*o1*u1*x*y*y + 
  k0*m2*o1*u1*x*y*y + k1*m0*o2*u1*x*y*y + k0*m1*o2*u1*x*y*y - 
  2*k1*k2*q0*u1*x*y*y - 2*k0*k2*q1*u1*x*y*y - 2*k0*k1*q2*u1*x*y*y - 
  j2*m1*r0*u1*x*y*y - j1*m2*r0*u1*x*y*y + i2*q1*r0*u1*x*y*y + 
  i1*q2*r0*u1*x*y*y - j2*m0*r1*u1*x*y*y - j0*m2*r1*u1*x*y*y + 
  i2*q0*r1*u1*x*y*y + i0*q2*r1*u1*x*y*y - j1*m0*r2*u1*x*y*y - 
  j0*m1*r2*u1*x*y*y + i1*q0*r2*u1*x*y*y + i0*q1*r2*u1*x*y*y + 
  j2*k1*t0*u1*x*y*y + j1*k2*t0*u1*x*y*y - i2*o1*t0*u1*x*y*y - 
  i1*o2*t0*u1*x*y*y + j2*k0*t1*u1*x*y*y + j0*k2*t1*u1*x*y*y - 
  i2*o0*t1*u1*x*y*y - i0*o2*t1*u1*x*y*y + j1*k0*t2*u1*x*y*y + 
  j0*k1*t2*u1*x*y*y - i1*o0*t2*u1*x*y*y - i0*o1*t2*u1*x*y*y + 
  k1*m1*o0*u2*x*y*y + k1*m0*o1*u2*x*y*y + k0*m1*o1*u2*x*y*y - 
  k1*k1*q0*u2*x*y*y - 2*k0*k1*q1*u2*x*y*y - j1*m1*r0*u2*x*y*y + 
  i1*q1*r0*u2*x*y*y - j1*m0*r1*u2*x*y*y - j0*m1*r1*u2*x*y*y + 
  i1*q0*r1*u2*x*y*y + i0*q1*r1*u2*x*y*y + j1*k1*t0*u2*x*y*y - 
  i1*o1*t0*u2*x*y*y + j1*k0*t1*u2*x*y*y + j0*k1*t1*u2*x*y*y - 
  i1*o0*t1*u2*x*y*y - i0*o1*t1*u2*x*y*y - k2*l1*o1*v0*x*y*y - 
  k1*l2*o1*v0*x*y*y - k1*l1*o2*v0*x*y*y + 2*k1*k2*p1*v0*x*y*y + 
  k1*k1*p2*v0*x*y*y + j2*l1*r1*v0*x*y*y + j1*l2*r1*v0*x*y*y - 
  i2*p1*r1*v0*x*y*y - i1*p2*r1*v0*x*y*y + j1*l1*r2*v0*x*y*y - 
  i1*p1*r2*v0*x*y*y - j2*k1*s1*v0*x*y*y - j1*k2*s1*v0*x*y*y + 
  i2*o1*s1*v0*x*y*y + i1*o2*s1*v0*x*y*y - j1*k1*s2*v0*x*y*y + 
  i1*o1*s2*v0*x*y*y - k2*l1*o0*v1*x*y*y - k1*l2*o0*v1*x*y*y - 
  k2*l0*o1*v1*x*y*y - k0*l2*o1*v1*x*y*y - k1*l0*o2*v1*x*y*y - 
  k0*l1*o2*v1*x*y*y + 2*k1*k2*p0*v1*x*y*y + 2*k0*k2*p1*v1*x*y*y + 
  2*k0*k1*p2*v1*x*y*y + j2*l1*r0*v1*x*y*y + j1*l2*r0*v1*x*y*y - 
  i2*p1*r0*v1*x*y*y - i1*p2*r0*v1*x*y*y + j2*l0*r1*v1*x*y*y + 
  j0*l2*r1*v1*x*y*y - i2*p0*r1*v1*x*y*y - i0*p2*r1*v1*x*y*y + 
  j1*l0*r2*v1*x*y*y + j0*l1*r2*v1*x*y*y - i1*p0*r2*v1*x*y*y - 
  i0*p1*r2*v1*x*y*y - j2*k1*s0*v1*x*y*y - j1*k2*s0*v1*x*y*y + 
  i2*o1*s0*v1*x*y*y + i1*o2*s0*v1*x*y*y - j2*k0*s1*v1*x*y*y - 
  j0*k2*s1*v1*x*y*y + i2*o0*s1*v1*x*y*y + i0*o2*s1*v1*x*y*y - 
  j1*k0*s2*v1*x*y*y - j0*k1*s2*v1*x*y*y + i1*o0*s2*v1*x*y*y + 
  i0*o1*s2*v1*x*y*y - k1*l1*o0*v2*x*y*y - k1*l0*o1*v2*x*y*y - 
  k0*l1*o1*v2*x*y*y + k1*k1*p0*v2*x*y*y + 2*k0*k1*p1*v2*x*y*y + 
  j1*l1*r0*v2*x*y*y - i1*p1*r0*v2*x*y*y + j1*l0*r1*v2*x*y*y + 
  j0*l1*r1*v2*x*y*y - i1*p0*r1*v2*x*y*y - i0*p1*r1*v2*x*y*y - 
  j1*k1*s0*v2*x*y*y + i1*o1*s0*v2*x*y*y - j1*k0*s1*v2*x*y*y - 
  j0*k1*s1*v2*x*y*y + i1*o0*s1*v2*x*y*y + i0*o1*s1*v2*x*y*y + 
  l1*m1*p0*r0*x*x*y*y + l1*m0*p1*r0*x*x*y*y + l0*m1*p1*r0*x*x*y*y - 
  l1*l1*q0*r0*x*x*y*y - 2*l0*l1*q1*r0*x*x*y*y + l1*m0*p0*r1*x*x*y*y + 
  l0*m1*p0*r1*x*x*y*y + l0*m0*p1*r1*x*x*y*y - 2*l0*l1*q0*r1*x*x*y*y - 
  l0*l0*q1*r1*x*x*y*y - l1*m1*o0*s0*x*x*y*y - l1*m0*o1*s0*x*x*y*y - 
  l0*m1*o1*s0*x*x*y*y - k1*m1*p0*s0*x*x*y*y - k1*m0*p1*s0*x*x*y*y - 
  k0*m1*p1*s0*x*x*y*y + 2*k1*l1*q0*s0*x*x*y*y + 2*k1*l0*q1*s0*x*x*y*y + 
  2*k0*l1*q1*s0*x*x*y*y + j1*m1*s0*s0*x*x*y*y - i1*q1*s0*s0*x*x*y*y - 
  l1*m0*o0*s1*x*x*y*y - l0*m1*o0*s1*x*x*y*y - l0*m0*o1*s1*x*x*y*y - 
  k1*m0*p0*s1*x*x*y*y - k0*m1*p0*s1*x*x*y*y - k0*m0*p1*s1*x*x*y*y + 
  2*k1*l0*q0*s1*x*x*y*y + 2*k0*l1*q0*s1*x*x*y*y + 2*k0*l0*q1*s1*x*x*y*y + 
  2*j1*m0*s0*s1*x*x*y*y + 2*j0*m1*s0*s1*x*x*y*y - 2*i1*q0*s0*s1*x*x*y*y - 
  2*i0*q1*s0*s1*x*x*y*y + j0*m0*s1*s1*x*x*y*y - i0*q0*s1*s1*x*x*y*y + 
  l1*l1*o0*t0*x*x*y*y + 2*l0*l1*o1*t0*x*x*y*y - k1*l1*p0*t0*x*x*y*y - 
  k1*l0*p1*t0*x*x*y*y - k0*l1*p1*t0*x*x*y*y - j1*l1*s0*t0*x*x*y*y + 
  i1*p1*s0*t0*x*x*y*y - j1*l0*s1*t0*x*x*y*y - j0*l1*s1*t0*x*x*y*y + 
  i1*p0*s1*t0*x*x*y*y + i0*p1*s1*t0*x*x*y*y + 2*l0*l1*o0*t1*x*x*y*y + 
  l0*l0*o1*t1*x*x*y*y - k1*l0*p0*t1*x*x*y*y - k0*l1*p0*t1*x*x*y*y - 
  k0*l0*p1*t1*x*x*y*y - j1*l0*s0*t1*x*x*y*y - j0*l1*s0*t1*x*x*y*y + 
  i1*p0*s0*t1*x*x*y*y + i0*p1*s0*t1*x*x*y*y - j0*l0*s1*t1*x*x*y*y + 
  i0*p0*s1*t1*x*x*y*y + k1*m1*o0*u0*x*x*y*y + k1*m0*o1*u0*x*x*y*y + 
  k0*m1*o1*u0*x*x*y*y - k1*k1*q0*u0*x*x*y*y - 2*k0*k1*q1*u0*x*x*y*y - 
  j1*m1*r0*u0*x*x*y*y + i1*q1*r0*u0*x*x*y*y - j1*m0*r1*u0*x*x*y*y - 
  j0*m1*r1*u0*x*x*y*y + i1*q0*r1*u0*x*x*y*y + i0*q1*r1*u0*x*x*y*y + 
  j1*k1*t0*u0*x*x*y*y - i1*o1*t0*u0*x*x*y*y + j1*k0*t1*u0*x*x*y*y + 
  j0*k1*t1*u0*x*x*y*y - i1*o0*t1*u0*x*x*y*y - i0*o1*t1*u0*x*x*y*y + 
  k1*m0*o0*u1*x*x*y*y + k0*m1*o0*u1*x*x*y*y + k0*m0*o1*u1*x*x*y*y - 
  2*k0*k1*q0*u1*x*x*y*y - k0*k0*q1*u1*x*x*y*y - j1*m0*r0*u1*x*x*y*y - 
  j0*m1*r0*u1*x*x*y*y + i1*q0*r0*u1*x*x*y*y + i0*q1*r0*u1*x*x*y*y - 
  j0*m0*r1*u1*x*x*y*y + i0*q0*r1*u1*x*x*y*y + j1*k0*t0*u1*x*x*y*y + 
  j0*k1*t0*u1*x*x*y*y - i1*o0*t0*u1*x*x*y*y - i0*o1*t0*u1*x*x*y*y + 
  j0*k0*t1*u1*x*x*y*y - i0*o0*t1*u1*x*x*y*y - k1*l1*o0*v0*x*x*y*y - 
  k1*l0*o1*v0*x*x*y*y - k0*l1*o1*v0*x*x*y*y + k1*k1*p0*v0*x*x*y*y + 
  2*k0*k1*p1*v0*x*x*y*y + j1*l1*r0*v0*x*x*y*y - i1*p1*r0*v0*x*x*y*y + 
  j1*l0*r1*v0*x*x*y*y + j0*l1*r1*v0*x*x*y*y - i1*p0*r1*v0*x*x*y*y - 
  i0*p1*r1*v0*x*x*y*y - j1*k1*s0*v0*x*x*y*y + i1*o1*s0*v0*x*x*y*y - 
  j1*k0*s1*v0*x*x*y*y - j0*k1*s1*v0*x*x*y*y + i1*o0*s1*v0*x*x*y*y + 
  i0*o1*s1*v0*x*x*y*y - k1*l0*o0*v1*x*x*y*y - k0*l1*o0*v1*x*x*y*y - 
  k0*l0*o1*v1*x*x*y*y + 2*k0*k1*p0*v1*x*x*y*y + k0*k0*p1*v1*x*x*y*y + 
  j1*l0*r0*v1*x*x*y*y + j0*l1*r0*v1*x*x*y*y - i1*p0*r0*v1*x*x*y*y - 
  i0*p1*r0*v1*x*x*y*y + j0*l0*r1*v1*x*x*y*y - i0*p0*r1*v1*x*x*y*y - 
  j1*k0*s0*v1*x*x*y*y - j0*k1*s0*v1*x*x*y*y + i1*o0*s0*v1*x*x*y*y + 
  i0*o1*s0*v1*x*x*y*y - j0*k0*s1*v1*x*x*y*y + i0*o0*s1*v1*x*x*y*y + 
  l2*m1*p1*r1*y*y*y + l1*m2*p1*r1*y*y*y + l1*m1*p2*r1*y*y*y - 2*l1*l2*q1*r1*y*y*y - 
  l1*l1*q2*r1*y*y*y + l1*m1*p1*r2*y*y*y - l1*l1*q1*r2*y*y*y - l2*m1*o1*s1*y*y*y - 
  l1*m2*o1*s1*y*y*y - l1*m1*o2*s1*y*y*y - k2*m1*p1*s1*y*y*y - k1*m2*p1*s1*y*y*y - 
  k1*m1*p2*s1*y*y*y + 2*k2*l1*q1*s1*y*y*y + 2*k1*l2*q1*s1*y*y*y + 
  2*k1*l1*q2*s1*y*y*y + j2*m1*s1*s1*y*y*y + j1*m2*s1*s1*y*y*y - i2*q1*s1*s1*y*y*y - 
  i1*q2*s1*s1*y*y*y - l1*m1*o1*s2*y*y*y - k1*m1*p1*s2*y*y*y + 2*k1*l1*q1*s2*y*y*y + 
  2*j1*m1*s1*s2*y*y*y - 2*i1*q1*s1*s2*y*y*y + 2*l1*l2*o1*t1*y*y*y + 
  l1*l1*o2*t1*y*y*y - k2*l1*p1*t1*y*y*y - k1*l2*p1*t1*y*y*y - k1*l1*p2*t1*y*y*y - 
  j2*l1*s1*t1*y*y*y - j1*l2*s1*t1*y*y*y + i2*p1*s1*t1*y*y*y + i1*p2*s1*t1*y*y*y - 
  j1*l1*s2*t1*y*y*y + i1*p1*s2*t1*y*y*y + l1*l1*o1*t2*y*y*y - k1*l1*p1*t2*y*y*y - 
  j1*l1*s1*t2*y*y*y + i1*p1*s1*t2*y*y*y + k2*m1*o1*u1*y*y*y + k1*m2*o1*u1*y*y*y + 
  k1*m1*o2*u1*y*y*y - 2*k1*k2*q1*u1*y*y*y - k1*k1*q2*u1*y*y*y - j2*m1*r1*u1*y*y*y - 
  j1*m2*r1*u1*y*y*y + i2*q1*r1*u1*y*y*y + i1*q2*r1*u1*y*y*y - j1*m1*r2*u1*y*y*y + 
  i1*q1*r2*u1*y*y*y + j2*k1*t1*u1*y*y*y + j1*k2*t1*u1*y*y*y - i2*o1*t1*u1*y*y*y - 
  i1*o2*t1*u1*y*y*y + j1*k1*t2*u1*y*y*y - i1*o1*t2*u1*y*y*y + k1*m1*o1*u2*y*y*y - 
  k1*k1*q1*u2*y*y*y - j1*m1*r1*u2*y*y*y + i1*q1*r1*u2*y*y*y + j1*k1*t1*u2*y*y*y - 
  i1*o1*t1*u2*y*y*y - k2*l1*o1*v1*y*y*y - k1*l2*o1*v1*y*y*y - k1*l1*o2*v1*y*y*y + 
  2*k1*k2*p1*v1*y*y*y + k1*k1*p2*v1*y*y*y + j2*l1*r1*v1*y*y*y + j1*l2*r1*v1*y*y*y - 
  i2*p1*r1*v1*y*y*y - i1*p2*r1*v1*y*y*y + j1*l1*r2*v1*y*y*y - i1*p1*r2*v1*y*y*y - 
  j2*k1*s1*v1*y*y*y - j1*k2*s1*v1*y*y*y + i2*o1*s1*v1*y*y*y + i1*o2*s1*v1*y*y*y - 
  j1*k1*s2*v1*y*y*y + i1*o1*s2*v1*y*y*y - k1*l1*o1*v2*y*y*y + k1*k1*p1*v2*y*y*y + 
  j1*l1*r1*v2*y*y*y - i1*p1*r1*v2*y*y*y - j1*k1*s1*v2*y*y*y + i1*o1*s1*v2*y*y*y + 
  l1*m1*p1*r0*x*y*y*y - l1*l1*q1*r0*x*y*y*y + l1*m1*p0*r1*x*y*y*y + 
  l1*m0*p1*r1*x*y*y*y + l0*m1*p1*r1*x*y*y*y - l1*l1*q0*r1*x*y*y*y - 
  2*l0*l1*q1*r1*x*y*y*y - l1*m1*o1*s0*x*y*y*y - k1*m1*p1*s0*x*y*y*y + 
  2*k1*l1*q1*s0*x*y*y*y - l1*m1*o0*s1*x*y*y*y - l1*m0*o1*s1*x*y*y*y - 
  l0*m1*o1*s1*x*y*y*y - k1*m1*p0*s1*x*y*y*y - k1*m0*p1*s1*x*y*y*y - 
  k0*m1*p1*s1*x*y*y*y + 2*k1*l1*q0*s1*x*y*y*y + 2*k1*l0*q1*s1*x*y*y*y + 
  2*k0*l1*q1*s1*x*y*y*y + 2*j1*m1*s0*s1*x*y*y*y - 2*i1*q1*s0*s1*x*y*y*y + 
  j1*m0*s1*s1*x*y*y*y + j0*m1*s1*s1*x*y*y*y - i1*q0*s1*s1*x*y*y*y - i0*q1*s1*s1*x*y*y*y + 
  l1*l1*o1*t0*x*y*y*y - k1*l1*p1*t0*x*y*y*y - j1*l1*s1*t0*x*y*y*y + 
  i1*p1*s1*t0*x*y*y*y + l1*l1*o0*t1*x*y*y*y + 2*l0*l1*o1*t1*x*y*y*y - 
  k1*l1*p0*t1*x*y*y*y - k1*l0*p1*t1*x*y*y*y - k0*l1*p1*t1*x*y*y*y - 
  j1*l1*s0*t1*x*y*y*y + i1*p1*s0*t1*x*y*y*y - j1*l0*s1*t1*x*y*y*y - 
  j0*l1*s1*t1*x*y*y*y + i1*p0*s1*t1*x*y*y*y + i0*p1*s1*t1*x*y*y*y + 
  k1*m1*o1*u0*x*y*y*y - k1*k1*q1*u0*x*y*y*y - j1*m1*r1*u0*x*y*y*y + 
  i1*q1*r1*u0*x*y*y*y + j1*k1*t1*u0*x*y*y*y - i1*o1*t1*u0*x*y*y*y + 
  k1*m1*o0*u1*x*y*y*y + k1*m0*o1*u1*x*y*y*y + k0*m1*o1*u1*x*y*y*y - 
  k1*k1*q0*u1*x*y*y*y - 2*k0*k1*q1*u1*x*y*y*y - j1*m1*r0*u1*x*y*y*y + 
  i1*q1*r0*u1*x*y*y*y - j1*m0*r1*u1*x*y*y*y - j0*m1*r1*u1*x*y*y*y + 
  i1*q0*r1*u1*x*y*y*y + i0*q1*r1*u1*x*y*y*y + j1*k1*t0*u1*x*y*y*y - 
  i1*o1*t0*u1*x*y*y*y + j1*k0*t1*u1*x*y*y*y + j0*k1*t1*u1*x*y*y*y - 
  i1*o0*t1*u1*x*y*y*y - i0*o1*t1*u1*x*y*y*y - k1*l1*o1*v0*x*y*y*y + 
  k1*k1*p1*v0*x*y*y*y + j1*l1*r1*v0*x*y*y*y - i1*p1*r1*v0*x*y*y*y - 
  j1*k1*s1*v0*x*y*y*y + i1*o1*s1*v0*x*y*y*y - k1*l1*o0*v1*x*y*y*y - 
  k1*l0*o1*v1*x*y*y*y - k0*l1*o1*v1*x*y*y*y + k1*k1*p0*v1*x*y*y*y + 
  2*k0*k1*p1*v1*x*y*y*y + j1*l1*r0*v1*x*y*y*y - i1*p1*r0*v1*x*y*y*y + 
  j1*l0*r1*v1*x*y*y*y + j0*l1*r1*v1*x*y*y*y - i1*p0*r1*v1*x*y*y*y - 
  i0*p1*r1*v1*x*y*y*y - j1*k1*s0*v1*x*y*y*y + i1*o1*s0*v1*x*y*y*y - 
  j1*k0*s1*v1*x*y*y*y - j0*k1*s1*v1*x*y*y*y + i1*o0*s1*v1*x*y*y*y + 
  i0*o1*s1*v1*x*y*y*y + l1*m1*p1*r1*y*y*y*y - l1*l1*q1*r1*y*y*y*y - l1*m1*o1*s1*y*y*y*y - 
  k1*m1*p1*s1*y*y*y*y + 2*k1*l1*q1*s1*y*y*y*y + j1*m1*s1*s1*y*y*y*y - i1*q1*s1*s1*y*y*y*y + 
  l1*l1*o1*t1*y*y*y*y - k1*l1*p1*t1*y*y*y*y - j1*l1*s1*t1*y*y*y*y + i1*p1*s1*t1*y*y*y*y + 
  k1*m1*o1*u1*y*y*y*y - k1*k1*q1*u1*y*y*y*y - j1*m1*r1*u1*y*y*y*y + i1*q1*r1*u1*y*y*y*y + 
  j1*k1*t1*u1*y*y*y*y - i1*o1*t1*u1*y*y*y*y - k1*l1*o1*v1*y*y*y*y + k1*k1*p1*v1*y*y*y*y + 
  j1*l1*r1*v1*y*y*y*y - i1*p1*r1*v1*y*y*y*y - j1*k1*s1*v1*y*y*y*y + i1*o1*s1*v1*y*y*y*y;

	return fabs(t_4/t_3);
}


