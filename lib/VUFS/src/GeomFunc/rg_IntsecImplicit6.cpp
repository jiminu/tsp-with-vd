#include "rg_IntsecImplicit6.h"
#include "rg_Polynomial.h"
#include "rg_RelativeOp.h"

// inserted by Jung-Hyun Ryu
#include "rg_IntersectFunc.h"

//#include <fstream>
#include <time.h>

rg_IntsecImplicit6::rg_IntsecImplicit6()
{

}

rg_IntsecImplicit6::~rg_IntsecImplicit6()
{

}

rg_dList<rg_Point2D> rg_IntsecImplicit6::intersectBzCurveVsBzCurve(const rg_BzCurve2D &curve_s, const rg_BzCurve2D &curve_t, rg_REAL &time)
{
	rg_dList<rg_Point2D> intersectPointList;
	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	//clock_t StartTime = clock(); 
	//for(rg_INT i=0; i<20; i++)
	//{

	intersectPointList.removeAll();
	clock_t StartTime = clock(); 
	//for(rg_INT o = 0;o < 19;o++) rg_IntersectFunc::implicitize(curve_t);
	rg_ImplicitEquation  impEq   = implicitize(curve_t);
	clock_t EndTime = clock();
	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	//ofstream fout("timeOut6.dat", ios::app);
	//fout << time << '\n';

	//rg_ImplicitEquation  impEq   = rg_IntersectFunc::implicitizeNotUsingMathematica2(curve_t);//test

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
//	}

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


//	clock_t EndTime = clock();
//  time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	
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

rg_ImplicitEquation rg_IntsecImplicit6::implicitize(const rg_BzCurve2D &curve)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL x3 = curve.getCtrlPt(3).getX();
	rg_REAL x4 = curve.getCtrlPt(4).getX();
    rg_REAL x5 = curve.getCtrlPt(5).getX();
	rg_REAL x6 = curve.getCtrlPt(6).getX();

	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL y3 = curve.getCtrlPt(3).getY();
	rg_REAL y4 = curve.getCtrlPt(4).getY();
    rg_REAL y5 = curve.getCtrlPt(5).getY();
	rg_REAL y6 = curve.getCtrlPt(6).getY();

//  If
//		    
//	x(t) = a6 t^6 + a5 t^5 + a4 t^4 + a3 t^3 + a2 t^2 + a1 t + a0
//
//	y(t) = b6 t^6 + b5 t^5 + b4 t^4 + b3 t^3 + b2 t^2 + b1 t + b0
//  
//  then the resultant of p(x, t) and q(y, t) is like the following.

	rg_REAL a6 =    x0 - 6*x1 + 15*x2 - 20*x3 + 15*x4 - 6*x5 + x6;
    rg_REAL a5 = -6*x0 + 30*x1 - 60*x2 + 60*x3 - 30*x4 + 6*x5;
	rg_REAL a4 = 15*x0 - 60*x1 + 90*x2 - 60*x3 + 15*x4;
	rg_REAL a3 =-20*x0 + 60*x1 - 60*x2 + 20*x3;
	rg_REAL a2 = 15*x0 - 30*x1 + 15*x2;
	rg_REAL a1 = -6*x0 + 6*x1;
	rg_REAL a0 =    x0;

	rg_REAL b6 = y0 - 6*y1 + 15*y2 - 20*y3 + 15*y4 - 6*y5 + y6;
    rg_REAL b5 = -6*y0 + 30*y1 - 60*y2 + 60*y3 - 30*y4 + 6*y5;
	rg_REAL b4 = 15*y0 - 60*y1 + 90*y2 - 60*y3 + 15*y4;
	rg_REAL b3 = -20*y0 + 60*y1 - 60*y2 + 20*y3;
	rg_REAL b2 = 15*y0 - 30*y1 + 15*y2;
	rg_REAL b1 = -6*y0 + 6*y1;
	rg_REAL b0 =    y0;

  //			|														                                                       |
  //			|	i0 x + i1 y + i2  j0 x + j1 y + j2	k0 x + k1 y + k2  l0 x + l1 y + l2  m0 x + m1 y + m2  n0 x + n1 y + n2 |
  //			|														                                                       |
  // R(P, Q) =	|	j0 x + j1 y + j2  o0 x + o1 y + o2	p0 x + p1 y + p2  q0 x + q1 y + q2  r0 x + r1 y + r2  s0 x + s1 y + s2 |
  //			|														                                                       |
  //			|	k0 x + k1 y + k2  p0 x + p1 y + p2	t0 x + t1 y + t2  u0 x + u1 y + u2  v0 x + v1 y + v2  w0 x + w1 y + w2 |
  //			|														                                                       |
  //            |   l0 x + l1 y + l2  q0 x + q1 y + q2  u0 x + u1 y + u2  z0 x + z1 y + z2  c0 x + c1 y + c2  d0 x + d1 y + d2 |
  //            |                                                                                                              |
  //            |   m0 x + m1 y + m2  r0 x + r1 y + r2  v0 x + v1 y + v2  c0 x + c1 y + c2  e0 x + e1 y + e2  f0 x + f1 y + f2 |
  //            |                                                                                                              |
  //            |   n0 x + n1 y + n2  s0 x + s1 y + s2  w0 x + w1 y + w2  d0 x + d1 y + d2  f0 x + f1 y + f2  g0 x + g1 y + g2 |
  //            |                                                                                                              |

	rg_REAL i0 =   0.0;
	rg_REAL i1 =   0.0;
	rg_REAL i2 =   a6*b5 - a5*b6;

	rg_REAL j0 =   0.0;
	rg_REAL j1 =   0.0;
	rg_REAL j2 =   a6*b4 - a4*b6;

	rg_REAL k0 =   0.0;
	rg_REAL k1 =   0.0;
	rg_REAL k2 =   a6*b3 - a3*b6;

	rg_REAL l0 =   0.0;
	rg_REAL l1 =   0.0;
	rg_REAL l2 =   a6*b2 - a2*b6;

	rg_REAL m0 =   0.0;
	rg_REAL m1 =   0.0;
	rg_REAL m2 =   a6*b1 - a1*b6;

	rg_REAL n0 =    b6;
	rg_REAL n1 =   -a6;
	rg_REAL n2 =   a6*b0 - a0*b6;

	rg_REAL o0 =   0.0;
	rg_REAL o1 =   0.0;
	rg_REAL o2 =   a5*b4 - a4*b5 + a6*b3 - a3*b6;

    rg_REAL p0 =   0.0;
	rg_REAL p1 =   0.0;
	rg_REAL p2 =   a5*b3 - a3*b5 + a6*b2 - a2*b6;

  	rg_REAL q0 =   0.0;
	rg_REAL q1 =   0.0;
	rg_REAL q2 =   a5*b2 - a2*b5 + a6*b1 - a1*b6;

  	rg_REAL r0 =   b6;
	rg_REAL r1 = - a6;
	rg_REAL r2 =   a5*b1 - a1*b5 + a6*b0 - a0*b6;

  	rg_REAL s0 =   b5;
	rg_REAL s1 = - a5;
	rg_REAL s2 =   a5*b0 - a0*b5;

   	rg_REAL t0 =   0.0;
	rg_REAL t1 =   0.0;
	rg_REAL t2 =   a4*b3 - a3*b4 + a5*b2 - a2*b5 + a6*b1 - a1*b6;

   	rg_REAL u0 =   b6;
	rg_REAL u1 = - a6;
	rg_REAL u2 =   a4*b2 - a2*b4 + a5*b1 - a1*b5 + a6*b0 - a0*b6;

	rg_REAL v0 =   b5;
	rg_REAL v1 = - a5;
	rg_REAL v2 =   a4*b1 - a1*b4 + a5*b0 - a0*b5;

   	rg_REAL w0 =   b4;
	rg_REAL w1 = - a4;
	rg_REAL w2 =   a4*b0 - a0*b4;

   	rg_REAL z0 =   b5;
	rg_REAL z1 = - a5;
	rg_REAL z2 =   a3*b2 - a2*b3 + a4*b1 - a1*b4 + a5*b0 - a0*b5;

   	rg_REAL c0 =   b4;
	rg_REAL c1 = - a4;
	rg_REAL c2 =   a3*b1 - a1*b3 + a4*b0 - a0*b4;

   	rg_REAL d0 =   b3;
	rg_REAL d1 = - a3;
	rg_REAL d2 =   a3*b0 - a0*b3;

   	rg_REAL e0 =   b3;
	rg_REAL e1 = - a3;
	rg_REAL e2 =   a2*b1 - a1*b2 + a3*b0 - a0*b3;

   	rg_REAL f0 =   b2;
	rg_REAL f1 = - a2;
	rg_REAL f2 =   a2*b0 - a0*b2;

   	rg_REAL g0 =   b1;
	rg_REAL g1 = - a1;
	rg_REAL g2 =   a1*b0 - a0*b1;

// If the form of the implicitization is 
//    k0 
//  + k1  x       + k2      y
//  + k3  x^2     + k4  x   y   + k5      y^2
//  + k6  x^3     + k7  x^2 y   + k8  x   y^2  + k9      y^3
//  + k10 x^4     + k11 x^3 y   + k12 x^2 y^2  + k13 x   y^3 + k14     y^4
//  + k15 x^5     + k16 x^4 y   + k17 x^3 y^2  + k18 x^2 y^3 + k19 x   y^4 + k20 y^5
//  + k21 x^6     + k22 x^5 y   + k23 x^4 y^2  + k24 x^3 y^3 + k25 x^2 y^4 + k26 x   y^5 + k27 y^6

    rg_REAL *k = new rg_REAL[28];
    
    k[0]  = pow(d2,2)*e2*pow(k2,2)*o2 - 2*c2*d2*f2*pow(k2,2)*o2 + 
   pow(c2,2)*g2*pow(k2,2)*o2 - 2*pow(d2,2)*e2*j2*k2*p2 + 
   4*c2*d2*f2*j2*k2*p2 - 2*pow(c2,2)*g2*j2*k2*p2 + 
   pow(d2,2)*e2*i2*pow(p2,2) - 2*c2*d2*f2*i2*pow(p2,2) + 
   pow(c2,2)*g2*i2*pow(p2,2) - pow(f2,2)*pow(l2,2)*pow(p2,2) + 
   e2*g2*pow(l2,2)*pow(p2,2) + 2*d2*f2*l2*m2*pow(p2,2) - 
   2*c2*g2*l2*m2*pow(p2,2) - pow(d2,2)*pow(m2,2)*pow(p2,2) - 
   2*d2*e2*l2*n2*pow(p2,2) + 2*c2*f2*l2*n2*pow(p2,2) + 
   2*c2*d2*m2*n2*pow(p2,2) - pow(c2,2)*pow(n2,2)*pow(p2,2) + 
   2*pow(f2,2)*k2*l2*p2*q2 - 2*e2*g2*k2*l2*p2*q2 - 2*d2*f2*k2*m2*p2*q2 + 
   2*c2*g2*k2*m2*p2*q2 + 2*d2*e2*k2*n2*p2*q2 - 2*c2*f2*k2*n2*p2*q2 - 
   pow(f2,2)*pow(k2,2)*pow(q2,2) + e2*g2*pow(k2,2)*pow(q2,2) - 
   2*d2*f2*k2*l2*p2*r2 + 2*c2*g2*k2*l2*p2*r2 + 2*pow(d2,2)*k2*m2*p2*r2 - 
   2*c2*d2*k2*n2*p2*r2 + 2*d2*f2*pow(k2,2)*q2*r2 - 
   2*c2*g2*pow(k2,2)*q2*r2 - pow(d2,2)*pow(k2,2)*pow(r2,2) + 
   2*d2*e2*k2*l2*p2*s2 - 2*c2*f2*k2*l2*p2*s2 - 2*c2*d2*k2*m2*p2*s2 + 
   2*pow(c2,2)*k2*n2*p2*s2 - 2*d2*e2*pow(k2,2)*q2*s2 + 
   2*c2*f2*pow(k2,2)*q2*s2 + 2*c2*d2*pow(k2,2)*r2*s2 - 
   pow(c2,2)*pow(k2,2)*pow(s2,2) + pow(d2,2)*e2*pow(j2,2)*t2 - 
   2*c2*d2*f2*pow(j2,2)*t2 + pow(c2,2)*g2*pow(j2,2)*t2 - 
   pow(d2,2)*e2*i2*o2*t2 + 2*c2*d2*f2*i2*o2*t2 - pow(c2,2)*g2*i2*o2*t2 + 
   pow(f2,2)*pow(l2,2)*o2*t2 - e2*g2*pow(l2,2)*o2*t2 - 
   2*d2*f2*l2*m2*o2*t2 + 2*c2*g2*l2*m2*o2*t2 + 
   pow(d2,2)*pow(m2,2)*o2*t2 + 2*d2*e2*l2*n2*o2*t2 - 
   2*c2*f2*l2*n2*o2*t2 - 2*c2*d2*m2*n2*o2*t2 + 
   pow(c2,2)*pow(n2,2)*o2*t2 - 2*pow(f2,2)*j2*l2*q2*t2 + 
   2*e2*g2*j2*l2*q2*t2 + 2*d2*f2*j2*m2*q2*t2 - 2*c2*g2*j2*m2*q2*t2 - 
   2*d2*e2*j2*n2*q2*t2 + 2*c2*f2*j2*n2*q2*t2 + 
   pow(f2,2)*i2*pow(q2,2)*t2 - e2*g2*i2*pow(q2,2)*t2 + 
   g2*pow(m2,2)*pow(q2,2)*t2 - 2*f2*m2*n2*pow(q2,2)*t2 + 
   e2*pow(n2,2)*pow(q2,2)*t2 + 2*d2*f2*j2*l2*r2*t2 - 
   2*c2*g2*j2*l2*r2*t2 - 2*pow(d2,2)*j2*m2*r2*t2 + 2*c2*d2*j2*n2*r2*t2 - 
   2*d2*f2*i2*q2*r2*t2 + 2*c2*g2*i2*q2*r2*t2 - 2*g2*l2*m2*q2*r2*t2 + 
   2*f2*l2*n2*q2*r2*t2 + 2*d2*m2*n2*q2*r2*t2 - 2*c2*pow(n2,2)*q2*r2*t2 + 
   pow(d2,2)*i2*pow(r2,2)*t2 + g2*pow(l2,2)*pow(r2,2)*t2 - 
   2*d2*l2*n2*pow(r2,2)*t2 - 2*d2*e2*j2*l2*s2*t2 + 2*c2*f2*j2*l2*s2*t2 + 
   2*c2*d2*j2*m2*s2*t2 - 2*pow(c2,2)*j2*n2*s2*t2 + 2*d2*e2*i2*q2*s2*t2 - 
   2*c2*f2*i2*q2*s2*t2 + 2*f2*l2*m2*q2*s2*t2 - 2*d2*pow(m2,2)*q2*s2*t2 - 
   2*e2*l2*n2*q2*s2*t2 + 2*c2*m2*n2*q2*s2*t2 - 2*c2*d2*i2*r2*s2*t2 - 
   2*f2*pow(l2,2)*r2*s2*t2 + 2*d2*l2*m2*r2*s2*t2 + 2*c2*l2*n2*r2*s2*t2 + 
   pow(c2,2)*i2*pow(s2,2)*t2 + e2*pow(l2,2)*pow(s2,2)*t2 - 
   2*c2*l2*m2*pow(s2,2)*t2 - 2*pow(f2,2)*k2*l2*o2*u2 + 
   2*e2*g2*k2*l2*o2*u2 + 2*d2*f2*k2*m2*o2*u2 - 2*c2*g2*k2*m2*o2*u2 - 
   2*d2*e2*k2*n2*o2*u2 + 2*c2*f2*k2*n2*o2*u2 + 2*pow(f2,2)*j2*l2*p2*u2 - 
   2*e2*g2*j2*l2*p2*u2 - 2*d2*f2*j2*m2*p2*u2 + 2*c2*g2*j2*m2*p2*u2 + 
   2*d2*e2*j2*n2*p2*u2 - 2*c2*f2*j2*n2*p2*u2 + 2*pow(f2,2)*j2*k2*q2*u2 - 
   2*e2*g2*j2*k2*q2*u2 - 2*pow(f2,2)*i2*p2*q2*u2 + 2*e2*g2*i2*p2*q2*u2 - 
   2*g2*pow(m2,2)*p2*q2*u2 + 4*f2*m2*n2*p2*q2*u2 - 
   2*e2*pow(n2,2)*p2*q2*u2 - 2*d2*f2*j2*k2*r2*u2 + 2*c2*g2*j2*k2*r2*u2 + 
   2*d2*f2*i2*p2*r2*u2 - 2*c2*g2*i2*p2*r2*u2 + 2*g2*l2*m2*p2*r2*u2 - 
   2*f2*l2*n2*p2*r2*u2 - 2*d2*m2*n2*p2*r2*u2 + 2*c2*pow(n2,2)*p2*r2*u2 + 
   2*g2*k2*m2*q2*r2*u2 - 2*f2*k2*n2*q2*r2*u2 - 2*g2*k2*l2*pow(r2,2)*u2 + 
   2*d2*k2*n2*pow(r2,2)*u2 + 2*d2*e2*j2*k2*s2*u2 - 2*c2*f2*j2*k2*s2*u2 - 
   2*d2*e2*i2*p2*s2*u2 + 2*c2*f2*i2*p2*s2*u2 - 2*f2*l2*m2*p2*s2*u2 + 
   2*d2*pow(m2,2)*p2*s2*u2 + 2*e2*l2*n2*p2*s2*u2 - 2*c2*m2*n2*p2*s2*u2 - 
   2*f2*k2*m2*q2*s2*u2 + 2*e2*k2*n2*q2*s2*u2 + 4*f2*k2*l2*r2*s2*u2 - 
   2*d2*k2*m2*r2*s2*u2 - 2*c2*k2*n2*r2*s2*u2 - 2*e2*k2*l2*pow(s2,2)*u2 + 
   2*c2*k2*m2*pow(s2,2)*u2 - pow(f2,2)*pow(j2,2)*pow(u2,2) + 
   e2*g2*pow(j2,2)*pow(u2,2) + pow(f2,2)*i2*o2*pow(u2,2) - 
   e2*g2*i2*o2*pow(u2,2) + g2*pow(m2,2)*o2*pow(u2,2) - 
   2*f2*m2*n2*o2*pow(u2,2) + e2*pow(n2,2)*o2*pow(u2,2) - 
   2*g2*j2*m2*r2*pow(u2,2) + 2*f2*j2*n2*r2*pow(u2,2) + 
   g2*i2*pow(r2,2)*pow(u2,2) - pow(n2,2)*pow(r2,2)*pow(u2,2) + 
   2*f2*j2*m2*s2*pow(u2,2) - 2*e2*j2*n2*s2*pow(u2,2) - 
   2*f2*i2*r2*s2*pow(u2,2) + 2*m2*n2*r2*s2*pow(u2,2) + 
   e2*i2*pow(s2,2)*pow(u2,2) - pow(m2,2)*pow(s2,2)*pow(u2,2) + 
   2*d2*f2*k2*l2*o2*v2 - 2*c2*g2*k2*l2*o2*v2 - 2*pow(d2,2)*k2*m2*o2*v2 + 
   2*c2*d2*k2*n2*o2*v2 - 2*d2*f2*j2*l2*p2*v2 + 2*c2*g2*j2*l2*p2*v2 + 
   2*pow(d2,2)*j2*m2*p2*v2 - 2*c2*d2*j2*n2*p2*v2 - 2*d2*f2*j2*k2*q2*v2 + 
   2*c2*g2*j2*k2*q2*v2 + 2*d2*f2*i2*p2*q2*v2 - 2*c2*g2*i2*p2*q2*v2 + 
   2*g2*l2*m2*p2*q2*v2 - 2*f2*l2*n2*p2*q2*v2 - 2*d2*m2*n2*p2*q2*v2 + 
   2*c2*pow(n2,2)*p2*q2*v2 - 2*g2*k2*m2*pow(q2,2)*v2 + 
   2*f2*k2*n2*pow(q2,2)*v2 + 2*pow(d2,2)*j2*k2*r2*v2 - 
   2*pow(d2,2)*i2*p2*r2*v2 - 2*g2*pow(l2,2)*p2*r2*v2 + 
   4*d2*l2*n2*p2*r2*v2 + 2*g2*k2*l2*q2*r2*v2 - 2*d2*k2*n2*q2*r2*v2 - 
   2*c2*d2*j2*k2*s2*v2 + 2*c2*d2*i2*p2*s2*v2 + 2*f2*pow(l2,2)*p2*s2*v2 - 
   2*d2*l2*m2*p2*s2*v2 - 2*c2*l2*n2*p2*s2*v2 - 2*f2*k2*l2*q2*s2*v2 + 
   4*d2*k2*m2*q2*s2*v2 - 2*c2*k2*n2*q2*s2*v2 - 2*d2*k2*l2*r2*s2*v2 + 
   2*c2*k2*l2*pow(s2,2)*v2 + 2*d2*f2*pow(j2,2)*u2*v2 - 
   2*c2*g2*pow(j2,2)*u2*v2 - 2*d2*f2*i2*o2*u2*v2 + 2*c2*g2*i2*o2*u2*v2 - 
   2*g2*l2*m2*o2*u2*v2 + 2*f2*l2*n2*o2*u2*v2 + 2*d2*m2*n2*o2*u2*v2 - 
   2*c2*pow(n2,2)*o2*u2*v2 + 2*g2*j2*m2*q2*u2*v2 - 2*f2*j2*n2*q2*u2*v2 + 
   2*g2*j2*l2*r2*u2*v2 - 2*d2*j2*n2*r2*u2*v2 - 2*g2*i2*q2*r2*u2*v2 + 
   2*pow(n2,2)*q2*r2*u2*v2 - 2*f2*j2*l2*s2*u2*v2 - 2*d2*j2*m2*s2*u2*v2 + 
   4*c2*j2*n2*s2*u2*v2 + 2*f2*i2*q2*s2*u2*v2 - 2*m2*n2*q2*s2*u2*v2 + 
   2*d2*i2*r2*s2*u2*v2 - 2*l2*n2*r2*s2*u2*v2 - 2*c2*i2*pow(s2,2)*u2*v2 + 
   2*l2*m2*pow(s2,2)*u2*v2 - pow(d2,2)*pow(j2,2)*pow(v2,2) + 
   pow(d2,2)*i2*o2*pow(v2,2) + g2*pow(l2,2)*o2*pow(v2,2) - 
   2*d2*l2*n2*o2*pow(v2,2) - 2*g2*j2*l2*q2*pow(v2,2) + 
   2*d2*j2*n2*q2*pow(v2,2) + g2*i2*pow(q2,2)*pow(v2,2) - 
   pow(n2,2)*pow(q2,2)*pow(v2,2) + 2*d2*j2*l2*s2*pow(v2,2) - 
   2*d2*i2*q2*s2*pow(v2,2) + 2*l2*n2*q2*s2*pow(v2,2) - 
   pow(l2,2)*pow(s2,2)*pow(v2,2) - 2*d2*e2*k2*l2*o2*w2 + 
   2*c2*f2*k2*l2*o2*w2 + 2*c2*d2*k2*m2*o2*w2 - 2*pow(c2,2)*k2*n2*o2*w2 + 
   2*d2*e2*j2*l2*p2*w2 - 2*c2*f2*j2*l2*p2*w2 - 2*c2*d2*j2*m2*p2*w2 + 
   2*pow(c2,2)*j2*n2*p2*w2 + 2*d2*e2*j2*k2*q2*w2 - 2*c2*f2*j2*k2*q2*w2 - 
   2*d2*e2*i2*p2*q2*w2 + 2*c2*f2*i2*p2*q2*w2 - 2*f2*l2*m2*p2*q2*w2 + 
   2*d2*pow(m2,2)*p2*q2*w2 + 2*e2*l2*n2*p2*q2*w2 - 2*c2*m2*n2*p2*q2*w2 + 
   2*f2*k2*m2*pow(q2,2)*w2 - 2*e2*k2*n2*pow(q2,2)*w2 - 
   2*c2*d2*j2*k2*r2*w2 + 2*c2*d2*i2*p2*r2*w2 + 2*f2*pow(l2,2)*p2*r2*w2 - 
   2*d2*l2*m2*p2*r2*w2 - 2*c2*l2*n2*p2*r2*w2 - 2*f2*k2*l2*q2*r2*w2 - 
   2*d2*k2*m2*q2*r2*w2 + 4*c2*k2*n2*q2*r2*w2 + 2*d2*k2*l2*pow(r2,2)*w2 + 
   2*pow(c2,2)*j2*k2*s2*w2 - 2*pow(c2,2)*i2*p2*s2*w2 - 
   2*e2*pow(l2,2)*p2*s2*w2 + 4*c2*l2*m2*p2*s2*w2 + 2*e2*k2*l2*q2*s2*w2 - 
   2*c2*k2*m2*q2*s2*w2 - 2*c2*k2*l2*r2*s2*w2 - 2*d2*e2*pow(j2,2)*u2*w2 + 
   2*c2*f2*pow(j2,2)*u2*w2 + 2*d2*e2*i2*o2*u2*w2 - 2*c2*f2*i2*o2*u2*w2 + 
   2*f2*l2*m2*o2*u2*w2 - 2*d2*pow(m2,2)*o2*u2*w2 - 2*e2*l2*n2*o2*u2*w2 + 
   2*c2*m2*n2*o2*u2*w2 - 2*f2*j2*m2*q2*u2*w2 + 2*e2*j2*n2*q2*u2*w2 - 
   2*f2*j2*l2*r2*u2*w2 + 4*d2*j2*m2*r2*u2*w2 - 2*c2*j2*n2*r2*u2*w2 + 
   2*f2*i2*q2*r2*u2*w2 - 2*m2*n2*q2*r2*u2*w2 - 2*d2*i2*pow(r2,2)*u2*w2 + 
   2*l2*n2*pow(r2,2)*u2*w2 + 2*e2*j2*l2*s2*u2*w2 - 2*c2*j2*m2*s2*u2*w2 - 
   2*e2*i2*q2*s2*u2*w2 + 2*pow(m2,2)*q2*s2*u2*w2 + 2*c2*i2*r2*s2*u2*w2 - 
   2*l2*m2*r2*s2*u2*w2 + 2*c2*d2*pow(j2,2)*v2*w2 - 2*c2*d2*i2*o2*v2*w2 - 
   2*f2*pow(l2,2)*o2*v2*w2 + 2*d2*l2*m2*o2*v2*w2 + 2*c2*l2*n2*o2*v2*w2 + 
   4*f2*j2*l2*q2*v2*w2 - 2*d2*j2*m2*q2*v2*w2 - 2*c2*j2*n2*q2*v2*w2 - 
   2*f2*i2*pow(q2,2)*v2*w2 + 2*m2*n2*pow(q2,2)*v2*w2 - 
   2*d2*j2*l2*r2*v2*w2 + 2*d2*i2*q2*r2*v2*w2 - 2*l2*n2*q2*r2*v2*w2 - 
   2*c2*j2*l2*s2*v2*w2 + 2*c2*i2*q2*s2*v2*w2 - 2*l2*m2*q2*s2*v2*w2 + 
   2*pow(l2,2)*r2*s2*v2*w2 - pow(c2,2)*pow(j2,2)*pow(w2,2) + 
   pow(c2,2)*i2*o2*pow(w2,2) + e2*pow(l2,2)*o2*pow(w2,2) - 
   2*c2*l2*m2*o2*pow(w2,2) - 2*e2*j2*l2*q2*pow(w2,2) + 
   2*c2*j2*m2*q2*pow(w2,2) + e2*i2*pow(q2,2)*pow(w2,2) - 
   pow(m2,2)*pow(q2,2)*pow(w2,2) + 2*c2*j2*l2*r2*pow(w2,2) - 
   2*c2*i2*q2*r2*pow(w2,2) + 2*l2*m2*q2*r2*pow(w2,2) - 
   pow(l2,2)*pow(r2,2)*pow(w2,2) + pow(f2,2)*pow(k2,2)*o2*z2 - 
   e2*g2*pow(k2,2)*o2*z2 - 2*pow(f2,2)*j2*k2*p2*z2 + 
   2*e2*g2*j2*k2*p2*z2 + pow(f2,2)*i2*pow(p2,2)*z2 - 
   e2*g2*i2*pow(p2,2)*z2 + g2*pow(m2,2)*pow(p2,2)*z2 - 
   2*f2*m2*n2*pow(p2,2)*z2 + e2*pow(n2,2)*pow(p2,2)*z2 - 
   2*g2*k2*m2*p2*r2*z2 + 2*f2*k2*n2*p2*r2*z2 + 
   g2*pow(k2,2)*pow(r2,2)*z2 + 2*f2*k2*m2*p2*s2*z2 - 
   2*e2*k2*n2*p2*s2*z2 - 2*f2*pow(k2,2)*r2*s2*z2 + 
   e2*pow(k2,2)*pow(s2,2)*z2 + pow(f2,2)*pow(j2,2)*t2*z2 - 
   e2*g2*pow(j2,2)*t2*z2 - pow(f2,2)*i2*o2*t2*z2 + e2*g2*i2*o2*t2*z2 - 
   g2*pow(m2,2)*o2*t2*z2 + 2*f2*m2*n2*o2*t2*z2 - e2*pow(n2,2)*o2*t2*z2 + 
   2*g2*j2*m2*r2*t2*z2 - 2*f2*j2*n2*r2*t2*z2 - g2*i2*pow(r2,2)*t2*z2 + 
   pow(n2,2)*pow(r2,2)*t2*z2 - 2*f2*j2*m2*s2*t2*z2 + 
   2*e2*j2*n2*s2*t2*z2 + 2*f2*i2*r2*s2*t2*z2 - 2*m2*n2*r2*s2*t2*z2 - 
   e2*i2*pow(s2,2)*t2*z2 + pow(m2,2)*pow(s2,2)*t2*z2 + 
   2*g2*k2*m2*o2*v2*z2 - 2*f2*k2*n2*o2*v2*z2 - 2*g2*j2*m2*p2*v2*z2 + 
   2*f2*j2*n2*p2*v2*z2 - 2*g2*j2*k2*r2*v2*z2 + 2*g2*i2*p2*r2*v2*z2 - 
   2*pow(n2,2)*p2*r2*v2*z2 + 2*f2*j2*k2*s2*v2*z2 - 2*f2*i2*p2*s2*v2*z2 + 
   2*m2*n2*p2*s2*v2*z2 + 2*k2*n2*r2*s2*v2*z2 - 2*k2*m2*pow(s2,2)*v2*z2 + 
   g2*pow(j2,2)*pow(v2,2)*z2 - g2*i2*o2*pow(v2,2)*z2 + 
   pow(n2,2)*o2*pow(v2,2)*z2 - 2*j2*n2*s2*pow(v2,2)*z2 + 
   i2*pow(s2,2)*pow(v2,2)*z2 - 2*f2*k2*m2*o2*w2*z2 + 
   2*e2*k2*n2*o2*w2*z2 + 2*f2*j2*m2*p2*w2*z2 - 2*e2*j2*n2*p2*w2*z2 + 
   2*f2*j2*k2*r2*w2*z2 - 2*f2*i2*p2*r2*w2*z2 + 2*m2*n2*p2*r2*w2*z2 - 
   2*k2*n2*pow(r2,2)*w2*z2 - 2*e2*j2*k2*s2*w2*z2 + 2*e2*i2*p2*s2*w2*z2 - 
   2*pow(m2,2)*p2*s2*w2*z2 + 2*k2*m2*r2*s2*w2*z2 - 
   2*f2*pow(j2,2)*v2*w2*z2 + 2*f2*i2*o2*v2*w2*z2 - 2*m2*n2*o2*v2*w2*z2 + 
   2*j2*n2*r2*v2*w2*z2 + 2*j2*m2*s2*v2*w2*z2 - 2*i2*r2*s2*v2*w2*z2 + 
   e2*pow(j2,2)*pow(w2,2)*z2 - e2*i2*o2*pow(w2,2)*z2 + 
   pow(m2,2)*o2*pow(w2,2)*z2 - 2*j2*m2*r2*pow(w2,2)*z2 + 
   i2*pow(r2,2)*pow(w2,2)*z2
;
    
    k[1]  = pow(d2,2)*e0*pow(k2,2)*o2 + 2*d0*d2*e2*pow(k2,2)*o2 - 
   2*c2*d2*f0*pow(k2,2)*o2 - 2*c2*d0*f2*pow(k2,2)*o2 - 
   2*c0*d2*f2*pow(k2,2)*o2 + pow(c2,2)*g0*pow(k2,2)*o2 + 
   2*c0*c2*g2*pow(k2,2)*o2 - 2*pow(d2,2)*e0*j2*k2*p2 - 
   4*d0*d2*e2*j2*k2*p2 + 4*c2*d2*f0*j2*k2*p2 + 4*c2*d0*f2*j2*k2*p2 + 
   4*c0*d2*f2*j2*k2*p2 - 2*pow(c2,2)*g0*j2*k2*p2 - 4*c0*c2*g2*j2*k2*p2 + 
   pow(d2,2)*e0*i2*pow(p2,2) + 2*d0*d2*e2*i2*pow(p2,2) - 
   2*c2*d2*f0*i2*pow(p2,2) - 2*c2*d0*f2*i2*pow(p2,2) - 
   2*c0*d2*f2*i2*pow(p2,2) + pow(c2,2)*g0*i2*pow(p2,2) + 
   2*c0*c2*g2*i2*pow(p2,2) - 2*f0*f2*pow(l2,2)*pow(p2,2) + 
   e2*g0*pow(l2,2)*pow(p2,2) + e0*g2*pow(l2,2)*pow(p2,2) + 
   2*d2*f0*l2*m2*pow(p2,2) + 2*d0*f2*l2*m2*pow(p2,2) - 
   2*c2*g0*l2*m2*pow(p2,2) - 2*c0*g2*l2*m2*pow(p2,2) - 
   2*d0*d2*pow(m2,2)*pow(p2,2) - 2*d2*e2*l2*n0*pow(p2,2) + 
   2*c2*f2*l2*n0*pow(p2,2) + 2*c2*d2*m2*n0*pow(p2,2) - 
   2*d2*e0*l2*n2*pow(p2,2) - 2*d0*e2*l2*n2*pow(p2,2) + 
   2*c2*f0*l2*n2*pow(p2,2) + 2*c0*f2*l2*n2*pow(p2,2) + 
   2*c2*d0*m2*n2*pow(p2,2) + 2*c0*d2*m2*n2*pow(p2,2) - 
   2*pow(c2,2)*n0*n2*pow(p2,2) - 2*c0*c2*pow(n2,2)*pow(p2,2) + 
   4*f0*f2*k2*l2*p2*q2 - 2*e2*g0*k2*l2*p2*q2 - 2*e0*g2*k2*l2*p2*q2 - 
   2*d2*f0*k2*m2*p2*q2 - 2*d0*f2*k2*m2*p2*q2 + 2*c2*g0*k2*m2*p2*q2 + 
   2*c0*g2*k2*m2*p2*q2 + 2*d2*e2*k2*n0*p2*q2 - 2*c2*f2*k2*n0*p2*q2 + 
   2*d2*e0*k2*n2*p2*q2 + 2*d0*e2*k2*n2*p2*q2 - 2*c2*f0*k2*n2*p2*q2 - 
   2*c0*f2*k2*n2*p2*q2 - 2*f0*f2*pow(k2,2)*pow(q2,2) + 
   e2*g0*pow(k2,2)*pow(q2,2) + e0*g2*pow(k2,2)*pow(q2,2) - 
   2*d2*f2*k2*l2*p2*r0 + 2*c2*g2*k2*l2*p2*r0 + 2*pow(d2,2)*k2*m2*p2*r0 - 
   2*c2*d2*k2*n2*p2*r0 + 2*d2*f2*pow(k2,2)*q2*r0 - 
   2*c2*g2*pow(k2,2)*q2*r0 - 2*d2*f0*k2*l2*p2*r2 - 2*d0*f2*k2*l2*p2*r2 + 
   2*c2*g0*k2*l2*p2*r2 + 2*c0*g2*k2*l2*p2*r2 + 4*d0*d2*k2*m2*p2*r2 - 
   2*c2*d2*k2*n0*p2*r2 - 2*c2*d0*k2*n2*p2*r2 - 2*c0*d2*k2*n2*p2*r2 + 
   2*d2*f0*pow(k2,2)*q2*r2 + 2*d0*f2*pow(k2,2)*q2*r2 - 
   2*c2*g0*pow(k2,2)*q2*r2 - 2*c0*g2*pow(k2,2)*q2*r2 - 
   2*pow(d2,2)*pow(k2,2)*r0*r2 - 2*d0*d2*pow(k2,2)*pow(r2,2) + 
   2*d2*e2*k2*l2*p2*s0 - 2*c2*f2*k2*l2*p2*s0 - 2*c2*d2*k2*m2*p2*s0 + 
   2*pow(c2,2)*k2*n2*p2*s0 - 2*d2*e2*pow(k2,2)*q2*s0 + 
   2*c2*f2*pow(k2,2)*q2*s0 + 2*c2*d2*pow(k2,2)*r2*s0 + 
   2*d2*e0*k2*l2*p2*s2 + 2*d0*e2*k2*l2*p2*s2 - 2*c2*f0*k2*l2*p2*s2 - 
   2*c0*f2*k2*l2*p2*s2 - 2*c2*d0*k2*m2*p2*s2 - 2*c0*d2*k2*m2*p2*s2 + 
   2*pow(c2,2)*k2*n0*p2*s2 + 4*c0*c2*k2*n2*p2*s2 - 
   2*d2*e0*pow(k2,2)*q2*s2 - 2*d0*e2*pow(k2,2)*q2*s2 + 
   2*c2*f0*pow(k2,2)*q2*s2 + 2*c0*f2*pow(k2,2)*q2*s2 + 
   2*c2*d2*pow(k2,2)*r0*s2 + 2*c2*d0*pow(k2,2)*r2*s2 + 
   2*c0*d2*pow(k2,2)*r2*s2 - 2*pow(c2,2)*pow(k2,2)*s0*s2 - 
   2*c0*c2*pow(k2,2)*pow(s2,2) + pow(d2,2)*e0*pow(j2,2)*t2 + 
   2*d0*d2*e2*pow(j2,2)*t2 - 2*c2*d2*f0*pow(j2,2)*t2 - 
   2*c2*d0*f2*pow(j2,2)*t2 - 2*c0*d2*f2*pow(j2,2)*t2 + 
   pow(c2,2)*g0*pow(j2,2)*t2 + 2*c0*c2*g2*pow(j2,2)*t2 - 
   pow(d2,2)*e0*i2*o2*t2 - 2*d0*d2*e2*i2*o2*t2 + 2*c2*d2*f0*i2*o2*t2 + 
   2*c2*d0*f2*i2*o2*t2 + 2*c0*d2*f2*i2*o2*t2 - pow(c2,2)*g0*i2*o2*t2 - 
   2*c0*c2*g2*i2*o2*t2 + 2*f0*f2*pow(l2,2)*o2*t2 - 
   e2*g0*pow(l2,2)*o2*t2 - e0*g2*pow(l2,2)*o2*t2 - 2*d2*f0*l2*m2*o2*t2 - 
   2*d0*f2*l2*m2*o2*t2 + 2*c2*g0*l2*m2*o2*t2 + 2*c0*g2*l2*m2*o2*t2 + 
   2*d0*d2*pow(m2,2)*o2*t2 + 2*d2*e2*l2*n0*o2*t2 - 2*c2*f2*l2*n0*o2*t2 - 
   2*c2*d2*m2*n0*o2*t2 + 2*d2*e0*l2*n2*o2*t2 + 2*d0*e2*l2*n2*o2*t2 - 
   2*c2*f0*l2*n2*o2*t2 - 2*c0*f2*l2*n2*o2*t2 - 2*c2*d0*m2*n2*o2*t2 - 
   2*c0*d2*m2*n2*o2*t2 + 2*pow(c2,2)*n0*n2*o2*t2 + 
   2*c0*c2*pow(n2,2)*o2*t2 - 4*f0*f2*j2*l2*q2*t2 + 2*e2*g0*j2*l2*q2*t2 + 
   2*e0*g2*j2*l2*q2*t2 + 2*d2*f0*j2*m2*q2*t2 + 2*d0*f2*j2*m2*q2*t2 - 
   2*c2*g0*j2*m2*q2*t2 - 2*c0*g2*j2*m2*q2*t2 - 2*d2*e2*j2*n0*q2*t2 + 
   2*c2*f2*j2*n0*q2*t2 - 2*d2*e0*j2*n2*q2*t2 - 2*d0*e2*j2*n2*q2*t2 + 
   2*c2*f0*j2*n2*q2*t2 + 2*c0*f2*j2*n2*q2*t2 + 2*f0*f2*i2*pow(q2,2)*t2 - 
   e2*g0*i2*pow(q2,2)*t2 - e0*g2*i2*pow(q2,2)*t2 + 
   g0*pow(m2,2)*pow(q2,2)*t2 - 2*f2*m2*n0*pow(q2,2)*t2 - 
   2*f0*m2*n2*pow(q2,2)*t2 + 2*e2*n0*n2*pow(q2,2)*t2 + 
   e0*pow(n2,2)*pow(q2,2)*t2 + 2*d2*f2*j2*l2*r0*t2 - 
   2*c2*g2*j2*l2*r0*t2 - 2*pow(d2,2)*j2*m2*r0*t2 + 2*c2*d2*j2*n2*r0*t2 - 
   2*d2*f2*i2*q2*r0*t2 + 2*c2*g2*i2*q2*r0*t2 - 2*g2*l2*m2*q2*r0*t2 + 
   2*f2*l2*n2*q2*r0*t2 + 2*d2*m2*n2*q2*r0*t2 - 2*c2*pow(n2,2)*q2*r0*t2 + 
   2*d2*f0*j2*l2*r2*t2 + 2*d0*f2*j2*l2*r2*t2 - 2*c2*g0*j2*l2*r2*t2 - 
   2*c0*g2*j2*l2*r2*t2 - 4*d0*d2*j2*m2*r2*t2 + 2*c2*d2*j2*n0*r2*t2 + 
   2*c2*d0*j2*n2*r2*t2 + 2*c0*d2*j2*n2*r2*t2 - 2*d2*f0*i2*q2*r2*t2 - 
   2*d0*f2*i2*q2*r2*t2 + 2*c2*g0*i2*q2*r2*t2 + 2*c0*g2*i2*q2*r2*t2 - 
   2*g0*l2*m2*q2*r2*t2 + 2*f2*l2*n0*q2*r2*t2 + 2*d2*m2*n0*q2*r2*t2 + 
   2*f0*l2*n2*q2*r2*t2 + 2*d0*m2*n2*q2*r2*t2 - 4*c2*n0*n2*q2*r2*t2 - 
   2*c0*pow(n2,2)*q2*r2*t2 + 2*pow(d2,2)*i2*r0*r2*t2 + 
   2*g2*pow(l2,2)*r0*r2*t2 - 4*d2*l2*n2*r0*r2*t2 + 
   2*d0*d2*i2*pow(r2,2)*t2 + g0*pow(l2,2)*pow(r2,2)*t2 - 
   2*d2*l2*n0*pow(r2,2)*t2 - 2*d0*l2*n2*pow(r2,2)*t2 - 
   2*d2*e2*j2*l2*s0*t2 + 2*c2*f2*j2*l2*s0*t2 + 2*c2*d2*j2*m2*s0*t2 - 
   2*pow(c2,2)*j2*n2*s0*t2 + 2*d2*e2*i2*q2*s0*t2 - 2*c2*f2*i2*q2*s0*t2 + 
   2*f2*l2*m2*q2*s0*t2 - 2*d2*pow(m2,2)*q2*s0*t2 - 2*e2*l2*n2*q2*s0*t2 + 
   2*c2*m2*n2*q2*s0*t2 - 2*c2*d2*i2*r2*s0*t2 - 2*f2*pow(l2,2)*r2*s0*t2 + 
   2*d2*l2*m2*r2*s0*t2 + 2*c2*l2*n2*r2*s0*t2 - 2*d2*e0*j2*l2*s2*t2 - 
   2*d0*e2*j2*l2*s2*t2 + 2*c2*f0*j2*l2*s2*t2 + 2*c0*f2*j2*l2*s2*t2 + 
   2*c2*d0*j2*m2*s2*t2 + 2*c0*d2*j2*m2*s2*t2 - 2*pow(c2,2)*j2*n0*s2*t2 - 
   4*c0*c2*j2*n2*s2*t2 + 2*d2*e0*i2*q2*s2*t2 + 2*d0*e2*i2*q2*s2*t2 - 
   2*c2*f0*i2*q2*s2*t2 - 2*c0*f2*i2*q2*s2*t2 + 2*f0*l2*m2*q2*s2*t2 - 
   2*d0*pow(m2,2)*q2*s2*t2 - 2*e2*l2*n0*q2*s2*t2 + 2*c2*m2*n0*q2*s2*t2 - 
   2*e0*l2*n2*q2*s2*t2 + 2*c0*m2*n2*q2*s2*t2 - 2*c2*d2*i2*r0*s2*t2 - 
   2*f2*pow(l2,2)*r0*s2*t2 + 2*d2*l2*m2*r0*s2*t2 + 2*c2*l2*n2*r0*s2*t2 - 
   2*c2*d0*i2*r2*s2*t2 - 2*c0*d2*i2*r2*s2*t2 - 2*f0*pow(l2,2)*r2*s2*t2 + 
   2*d0*l2*m2*r2*s2*t2 + 2*c2*l2*n0*r2*s2*t2 + 2*c0*l2*n2*r2*s2*t2 + 
   2*pow(c2,2)*i2*s0*s2*t2 + 2*e2*pow(l2,2)*s0*s2*t2 - 
   4*c2*l2*m2*s0*s2*t2 + 2*c0*c2*i2*pow(s2,2)*t2 + 
   e0*pow(l2,2)*pow(s2,2)*t2 - 2*c0*l2*m2*pow(s2,2)*t2 - 
   2*pow(f2,2)*k2*l2*o2*u0 + 2*e2*g2*k2*l2*o2*u0 + 2*d2*f2*k2*m2*o2*u0 - 
   2*c2*g2*k2*m2*o2*u0 - 2*d2*e2*k2*n2*o2*u0 + 2*c2*f2*k2*n2*o2*u0 + 
   2*pow(f2,2)*j2*l2*p2*u0 - 2*e2*g2*j2*l2*p2*u0 - 2*d2*f2*j2*m2*p2*u0 + 
   2*c2*g2*j2*m2*p2*u0 + 2*d2*e2*j2*n2*p2*u0 - 2*c2*f2*j2*n2*p2*u0 + 
   2*pow(f2,2)*j2*k2*q2*u0 - 2*e2*g2*j2*k2*q2*u0 - 
   2*pow(f2,2)*i2*p2*q2*u0 + 2*e2*g2*i2*p2*q2*u0 - 
   2*g2*pow(m2,2)*p2*q2*u0 + 4*f2*m2*n2*p2*q2*u0 - 
   2*e2*pow(n2,2)*p2*q2*u0 - 2*d2*f2*j2*k2*r2*u0 + 2*c2*g2*j2*k2*r2*u0 + 
   2*d2*f2*i2*p2*r2*u0 - 2*c2*g2*i2*p2*r2*u0 + 2*g2*l2*m2*p2*r2*u0 - 
   2*f2*l2*n2*p2*r2*u0 - 2*d2*m2*n2*p2*r2*u0 + 2*c2*pow(n2,2)*p2*r2*u0 + 
   2*g2*k2*m2*q2*r2*u0 - 2*f2*k2*n2*q2*r2*u0 - 2*g2*k2*l2*pow(r2,2)*u0 + 
   2*d2*k2*n2*pow(r2,2)*u0 + 2*d2*e2*j2*k2*s2*u0 - 2*c2*f2*j2*k2*s2*u0 - 
   2*d2*e2*i2*p2*s2*u0 + 2*c2*f2*i2*p2*s2*u0 - 2*f2*l2*m2*p2*s2*u0 + 
   2*d2*pow(m2,2)*p2*s2*u0 + 2*e2*l2*n2*p2*s2*u0 - 2*c2*m2*n2*p2*s2*u0 - 
   2*f2*k2*m2*q2*s2*u0 + 2*e2*k2*n2*q2*s2*u0 + 4*f2*k2*l2*r2*s2*u0 - 
   2*d2*k2*m2*r2*s2*u0 - 2*c2*k2*n2*r2*s2*u0 - 2*e2*k2*l2*pow(s2,2)*u0 + 
   2*c2*k2*m2*pow(s2,2)*u0 - 4*f0*f2*k2*l2*o2*u2 + 2*e2*g0*k2*l2*o2*u2 + 
   2*e0*g2*k2*l2*o2*u2 + 2*d2*f0*k2*m2*o2*u2 + 2*d0*f2*k2*m2*o2*u2 - 
   2*c2*g0*k2*m2*o2*u2 - 2*c0*g2*k2*m2*o2*u2 - 2*d2*e2*k2*n0*o2*u2 + 
   2*c2*f2*k2*n0*o2*u2 - 2*d2*e0*k2*n2*o2*u2 - 2*d0*e2*k2*n2*o2*u2 + 
   2*c2*f0*k2*n2*o2*u2 + 2*c0*f2*k2*n2*o2*u2 + 4*f0*f2*j2*l2*p2*u2 - 
   2*e2*g0*j2*l2*p2*u2 - 2*e0*g2*j2*l2*p2*u2 - 2*d2*f0*j2*m2*p2*u2 - 
   2*d0*f2*j2*m2*p2*u2 + 2*c2*g0*j2*m2*p2*u2 + 2*c0*g2*j2*m2*p2*u2 + 
   2*d2*e2*j2*n0*p2*u2 - 2*c2*f2*j2*n0*p2*u2 + 2*d2*e0*j2*n2*p2*u2 + 
   2*d0*e2*j2*n2*p2*u2 - 2*c2*f0*j2*n2*p2*u2 - 2*c0*f2*j2*n2*p2*u2 + 
   4*f0*f2*j2*k2*q2*u2 - 2*e2*g0*j2*k2*q2*u2 - 2*e0*g2*j2*k2*q2*u2 - 
   4*f0*f2*i2*p2*q2*u2 + 2*e2*g0*i2*p2*q2*u2 + 2*e0*g2*i2*p2*q2*u2 - 
   2*g0*pow(m2,2)*p2*q2*u2 + 4*f2*m2*n0*p2*q2*u2 + 4*f0*m2*n2*p2*q2*u2 - 
   4*e2*n0*n2*p2*q2*u2 - 2*e0*pow(n2,2)*p2*q2*u2 - 2*d2*f2*j2*k2*r0*u2 + 
   2*c2*g2*j2*k2*r0*u2 + 2*d2*f2*i2*p2*r0*u2 - 2*c2*g2*i2*p2*r0*u2 + 
   2*g2*l2*m2*p2*r0*u2 - 2*f2*l2*n2*p2*r0*u2 - 2*d2*m2*n2*p2*r0*u2 + 
   2*c2*pow(n2,2)*p2*r0*u2 + 2*g2*k2*m2*q2*r0*u2 - 2*f2*k2*n2*q2*r0*u2 - 
   2*d2*f0*j2*k2*r2*u2 - 2*d0*f2*j2*k2*r2*u2 + 2*c2*g0*j2*k2*r2*u2 + 
   2*c0*g2*j2*k2*r2*u2 + 2*d2*f0*i2*p2*r2*u2 + 2*d0*f2*i2*p2*r2*u2 - 
   2*c2*g0*i2*p2*r2*u2 - 2*c0*g2*i2*p2*r2*u2 + 2*g0*l2*m2*p2*r2*u2 - 
   2*f2*l2*n0*p2*r2*u2 - 2*d2*m2*n0*p2*r2*u2 - 2*f0*l2*n2*p2*r2*u2 - 
   2*d0*m2*n2*p2*r2*u2 + 4*c2*n0*n2*p2*r2*u2 + 2*c0*pow(n2,2)*p2*r2*u2 + 
   2*g0*k2*m2*q2*r2*u2 - 2*f2*k2*n0*q2*r2*u2 - 2*f0*k2*n2*q2*r2*u2 - 
   4*g2*k2*l2*r0*r2*u2 + 4*d2*k2*n2*r0*r2*u2 - 2*g0*k2*l2*pow(r2,2)*u2 + 
   2*d2*k2*n0*pow(r2,2)*u2 + 2*d0*k2*n2*pow(r2,2)*u2 + 
   2*d2*e2*j2*k2*s0*u2 - 2*c2*f2*j2*k2*s0*u2 - 2*d2*e2*i2*p2*s0*u2 + 
   2*c2*f2*i2*p2*s0*u2 - 2*f2*l2*m2*p2*s0*u2 + 2*d2*pow(m2,2)*p2*s0*u2 + 
   2*e2*l2*n2*p2*s0*u2 - 2*c2*m2*n2*p2*s0*u2 - 2*f2*k2*m2*q2*s0*u2 + 
   2*e2*k2*n2*q2*s0*u2 + 4*f2*k2*l2*r2*s0*u2 - 2*d2*k2*m2*r2*s0*u2 - 
   2*c2*k2*n2*r2*s0*u2 + 2*d2*e0*j2*k2*s2*u2 + 2*d0*e2*j2*k2*s2*u2 - 
   2*c2*f0*j2*k2*s2*u2 - 2*c0*f2*j2*k2*s2*u2 - 2*d2*e0*i2*p2*s2*u2 - 
   2*d0*e2*i2*p2*s2*u2 + 2*c2*f0*i2*p2*s2*u2 + 2*c0*f2*i2*p2*s2*u2 - 
   2*f0*l2*m2*p2*s2*u2 + 2*d0*pow(m2,2)*p2*s2*u2 + 2*e2*l2*n0*p2*s2*u2 - 
   2*c2*m2*n0*p2*s2*u2 + 2*e0*l2*n2*p2*s2*u2 - 2*c0*m2*n2*p2*s2*u2 - 
   2*f0*k2*m2*q2*s2*u2 + 2*e2*k2*n0*q2*s2*u2 + 2*e0*k2*n2*q2*s2*u2 + 
   4*f2*k2*l2*r0*s2*u2 - 2*d2*k2*m2*r0*s2*u2 - 2*c2*k2*n2*r0*s2*u2 + 
   4*f0*k2*l2*r2*s2*u2 - 2*d0*k2*m2*r2*s2*u2 - 2*c2*k2*n0*r2*s2*u2 - 
   2*c0*k2*n2*r2*s2*u2 - 4*e2*k2*l2*s0*s2*u2 + 4*c2*k2*m2*s0*s2*u2 - 
   2*e0*k2*l2*pow(s2,2)*u2 + 2*c0*k2*m2*pow(s2,2)*u2 - 
   2*pow(f2,2)*pow(j2,2)*u0*u2 + 2*e2*g2*pow(j2,2)*u0*u2 + 
   2*pow(f2,2)*i2*o2*u0*u2 - 2*e2*g2*i2*o2*u0*u2 + 
   2*g2*pow(m2,2)*o2*u0*u2 - 4*f2*m2*n2*o2*u0*u2 + 
   2*e2*pow(n2,2)*o2*u0*u2 - 4*g2*j2*m2*r2*u0*u2 + 4*f2*j2*n2*r2*u0*u2 + 
   2*g2*i2*pow(r2,2)*u0*u2 - 2*pow(n2,2)*pow(r2,2)*u0*u2 + 
   4*f2*j2*m2*s2*u0*u2 - 4*e2*j2*n2*s2*u0*u2 - 4*f2*i2*r2*s2*u0*u2 + 
   4*m2*n2*r2*s2*u0*u2 + 2*e2*i2*pow(s2,2)*u0*u2 - 
   2*pow(m2,2)*pow(s2,2)*u0*u2 - 2*f0*f2*pow(j2,2)*pow(u2,2) + 
   e2*g0*pow(j2,2)*pow(u2,2) + e0*g2*pow(j2,2)*pow(u2,2) + 
   2*f0*f2*i2*o2*pow(u2,2) - e2*g0*i2*o2*pow(u2,2) - 
   e0*g2*i2*o2*pow(u2,2) + g0*pow(m2,2)*o2*pow(u2,2) - 
   2*f2*m2*n0*o2*pow(u2,2) - 2*f0*m2*n2*o2*pow(u2,2) + 
   2*e2*n0*n2*o2*pow(u2,2) + e0*pow(n2,2)*o2*pow(u2,2) - 
   2*g2*j2*m2*r0*pow(u2,2) + 2*f2*j2*n2*r0*pow(u2,2) - 
   2*g0*j2*m2*r2*pow(u2,2) + 2*f2*j2*n0*r2*pow(u2,2) + 
   2*f0*j2*n2*r2*pow(u2,2) + 2*g2*i2*r0*r2*pow(u2,2) - 
   2*pow(n2,2)*r0*r2*pow(u2,2) + g0*i2*pow(r2,2)*pow(u2,2) - 
   2*n0*n2*pow(r2,2)*pow(u2,2) + 2*f2*j2*m2*s0*pow(u2,2) - 
   2*e2*j2*n2*s0*pow(u2,2) - 2*f2*i2*r2*s0*pow(u2,2) + 
   2*m2*n2*r2*s0*pow(u2,2) + 2*f0*j2*m2*s2*pow(u2,2) - 
   2*e2*j2*n0*s2*pow(u2,2) - 2*e0*j2*n2*s2*pow(u2,2) - 
   2*f2*i2*r0*s2*pow(u2,2) + 2*m2*n2*r0*s2*pow(u2,2) - 
   2*f0*i2*r2*s2*pow(u2,2) + 2*m2*n0*r2*s2*pow(u2,2) + 
   2*e2*i2*s0*s2*pow(u2,2) - 2*pow(m2,2)*s0*s2*pow(u2,2) + 
   e0*i2*pow(s2,2)*pow(u2,2) + 2*d2*f2*k2*l2*o2*v0 - 
   2*c2*g2*k2*l2*o2*v0 - 2*pow(d2,2)*k2*m2*o2*v0 + 2*c2*d2*k2*n2*o2*v0 - 
   2*d2*f2*j2*l2*p2*v0 + 2*c2*g2*j2*l2*p2*v0 + 2*pow(d2,2)*j2*m2*p2*v0 - 
   2*c2*d2*j2*n2*p2*v0 - 2*d2*f2*j2*k2*q2*v0 + 2*c2*g2*j2*k2*q2*v0 + 
   2*d2*f2*i2*p2*q2*v0 - 2*c2*g2*i2*p2*q2*v0 + 2*g2*l2*m2*p2*q2*v0 - 
   2*f2*l2*n2*p2*q2*v0 - 2*d2*m2*n2*p2*q2*v0 + 2*c2*pow(n2,2)*p2*q2*v0 - 
   2*g2*k2*m2*pow(q2,2)*v0 + 2*f2*k2*n2*pow(q2,2)*v0 + 
   2*pow(d2,2)*j2*k2*r2*v0 - 2*pow(d2,2)*i2*p2*r2*v0 - 
   2*g2*pow(l2,2)*p2*r2*v0 + 4*d2*l2*n2*p2*r2*v0 + 2*g2*k2*l2*q2*r2*v0 - 
   2*d2*k2*n2*q2*r2*v0 - 2*c2*d2*j2*k2*s2*v0 + 2*c2*d2*i2*p2*s2*v0 + 
   2*f2*pow(l2,2)*p2*s2*v0 - 2*d2*l2*m2*p2*s2*v0 - 2*c2*l2*n2*p2*s2*v0 - 
   2*f2*k2*l2*q2*s2*v0 + 4*d2*k2*m2*q2*s2*v0 - 2*c2*k2*n2*q2*s2*v0 - 
   2*d2*k2*l2*r2*s2*v0 + 2*c2*k2*l2*pow(s2,2)*v0 + 
   2*d2*f2*pow(j2,2)*u2*v0 - 2*c2*g2*pow(j2,2)*u2*v0 - 
   2*d2*f2*i2*o2*u2*v0 + 2*c2*g2*i2*o2*u2*v0 - 2*g2*l2*m2*o2*u2*v0 + 
   2*f2*l2*n2*o2*u2*v0 + 2*d2*m2*n2*o2*u2*v0 - 2*c2*pow(n2,2)*o2*u2*v0 + 
   2*g2*j2*m2*q2*u2*v0 - 2*f2*j2*n2*q2*u2*v0 + 2*g2*j2*l2*r2*u2*v0 - 
   2*d2*j2*n2*r2*u2*v0 - 2*g2*i2*q2*r2*u2*v0 + 2*pow(n2,2)*q2*r2*u2*v0 - 
   2*f2*j2*l2*s2*u2*v0 - 2*d2*j2*m2*s2*u2*v0 + 4*c2*j2*n2*s2*u2*v0 + 
   2*f2*i2*q2*s2*u2*v0 - 2*m2*n2*q2*s2*u2*v0 + 2*d2*i2*r2*s2*u2*v0 - 
   2*l2*n2*r2*s2*u2*v0 - 2*c2*i2*pow(s2,2)*u2*v0 + 
   2*l2*m2*pow(s2,2)*u2*v0 + 2*d2*f0*k2*l2*o2*v2 + 2*d0*f2*k2*l2*o2*v2 - 
   2*c2*g0*k2*l2*o2*v2 - 2*c0*g2*k2*l2*o2*v2 - 4*d0*d2*k2*m2*o2*v2 + 
   2*c2*d2*k2*n0*o2*v2 + 2*c2*d0*k2*n2*o2*v2 + 2*c0*d2*k2*n2*o2*v2 - 
   2*d2*f0*j2*l2*p2*v2 - 2*d0*f2*j2*l2*p2*v2 + 2*c2*g0*j2*l2*p2*v2 + 
   2*c0*g2*j2*l2*p2*v2 + 4*d0*d2*j2*m2*p2*v2 - 2*c2*d2*j2*n0*p2*v2 - 
   2*c2*d0*j2*n2*p2*v2 - 2*c0*d2*j2*n2*p2*v2 - 2*d2*f0*j2*k2*q2*v2 - 
   2*d0*f2*j2*k2*q2*v2 + 2*c2*g0*j2*k2*q2*v2 + 2*c0*g2*j2*k2*q2*v2 + 
   2*d2*f0*i2*p2*q2*v2 + 2*d0*f2*i2*p2*q2*v2 - 2*c2*g0*i2*p2*q2*v2 - 
   2*c0*g2*i2*p2*q2*v2 + 2*g0*l2*m2*p2*q2*v2 - 2*f2*l2*n0*p2*q2*v2 - 
   2*d2*m2*n0*p2*q2*v2 - 2*f0*l2*n2*p2*q2*v2 - 2*d0*m2*n2*p2*q2*v2 + 
   4*c2*n0*n2*p2*q2*v2 + 2*c0*pow(n2,2)*p2*q2*v2 - 
   2*g0*k2*m2*pow(q2,2)*v2 + 2*f2*k2*n0*pow(q2,2)*v2 + 
   2*f0*k2*n2*pow(q2,2)*v2 + 2*pow(d2,2)*j2*k2*r0*v2 - 
   2*pow(d2,2)*i2*p2*r0*v2 - 2*g2*pow(l2,2)*p2*r0*v2 + 
   4*d2*l2*n2*p2*r0*v2 + 2*g2*k2*l2*q2*r0*v2 - 2*d2*k2*n2*q2*r0*v2 + 
   4*d0*d2*j2*k2*r2*v2 - 4*d0*d2*i2*p2*r2*v2 - 2*g0*pow(l2,2)*p2*r2*v2 + 
   4*d2*l2*n0*p2*r2*v2 + 4*d0*l2*n2*p2*r2*v2 + 2*g0*k2*l2*q2*r2*v2 - 
   2*d2*k2*n0*q2*r2*v2 - 2*d0*k2*n2*q2*r2*v2 - 2*c2*d2*j2*k2*s0*v2 + 
   2*c2*d2*i2*p2*s0*v2 + 2*f2*pow(l2,2)*p2*s0*v2 - 2*d2*l2*m2*p2*s0*v2 - 
   2*c2*l2*n2*p2*s0*v2 - 2*f2*k2*l2*q2*s0*v2 + 4*d2*k2*m2*q2*s0*v2 - 
   2*c2*k2*n2*q2*s0*v2 - 2*d2*k2*l2*r2*s0*v2 - 2*c2*d0*j2*k2*s2*v2 - 
   2*c0*d2*j2*k2*s2*v2 + 2*c2*d0*i2*p2*s2*v2 + 2*c0*d2*i2*p2*s2*v2 + 
   2*f0*pow(l2,2)*p2*s2*v2 - 2*d0*l2*m2*p2*s2*v2 - 2*c2*l2*n0*p2*s2*v2 - 
   2*c0*l2*n2*p2*s2*v2 - 2*f0*k2*l2*q2*s2*v2 + 4*d0*k2*m2*q2*s2*v2 - 
   2*c2*k2*n0*q2*s2*v2 - 2*c0*k2*n2*q2*s2*v2 - 2*d2*k2*l2*r0*s2*v2 - 
   2*d0*k2*l2*r2*s2*v2 + 4*c2*k2*l2*s0*s2*v2 + 2*c0*k2*l2*pow(s2,2)*v2 + 
   2*d2*f2*pow(j2,2)*u0*v2 - 2*c2*g2*pow(j2,2)*u0*v2 - 
   2*d2*f2*i2*o2*u0*v2 + 2*c2*g2*i2*o2*u0*v2 - 2*g2*l2*m2*o2*u0*v2 + 
   2*f2*l2*n2*o2*u0*v2 + 2*d2*m2*n2*o2*u0*v2 - 2*c2*pow(n2,2)*o2*u0*v2 + 
   2*g2*j2*m2*q2*u0*v2 - 2*f2*j2*n2*q2*u0*v2 + 2*g2*j2*l2*r2*u0*v2 - 
   2*d2*j2*n2*r2*u0*v2 - 2*g2*i2*q2*r2*u0*v2 + 2*pow(n2,2)*q2*r2*u0*v2 - 
   2*f2*j2*l2*s2*u0*v2 - 2*d2*j2*m2*s2*u0*v2 + 4*c2*j2*n2*s2*u0*v2 + 
   2*f2*i2*q2*s2*u0*v2 - 2*m2*n2*q2*s2*u0*v2 + 2*d2*i2*r2*s2*u0*v2 - 
   2*l2*n2*r2*s2*u0*v2 - 2*c2*i2*pow(s2,2)*u0*v2 + 
   2*l2*m2*pow(s2,2)*u0*v2 + 2*d2*f0*pow(j2,2)*u2*v2 + 
   2*d0*f2*pow(j2,2)*u2*v2 - 2*c2*g0*pow(j2,2)*u2*v2 - 
   2*c0*g2*pow(j2,2)*u2*v2 - 2*d2*f0*i2*o2*u2*v2 - 2*d0*f2*i2*o2*u2*v2 + 
   2*c2*g0*i2*o2*u2*v2 + 2*c0*g2*i2*o2*u2*v2 - 2*g0*l2*m2*o2*u2*v2 + 
   2*f2*l2*n0*o2*u2*v2 + 2*d2*m2*n0*o2*u2*v2 + 2*f0*l2*n2*o2*u2*v2 + 
   2*d0*m2*n2*o2*u2*v2 - 4*c2*n0*n2*o2*u2*v2 - 2*c0*pow(n2,2)*o2*u2*v2 + 
   2*g0*j2*m2*q2*u2*v2 - 2*f2*j2*n0*q2*u2*v2 - 2*f0*j2*n2*q2*u2*v2 + 
   2*g2*j2*l2*r0*u2*v2 - 2*d2*j2*n2*r0*u2*v2 - 2*g2*i2*q2*r0*u2*v2 + 
   2*pow(n2,2)*q2*r0*u2*v2 + 2*g0*j2*l2*r2*u2*v2 - 2*d2*j2*n0*r2*u2*v2 - 
   2*d0*j2*n2*r2*u2*v2 - 2*g0*i2*q2*r2*u2*v2 + 4*n0*n2*q2*r2*u2*v2 - 
   2*f2*j2*l2*s0*u2*v2 - 2*d2*j2*m2*s0*u2*v2 + 4*c2*j2*n2*s0*u2*v2 + 
   2*f2*i2*q2*s0*u2*v2 - 2*m2*n2*q2*s0*u2*v2 + 2*d2*i2*r2*s0*u2*v2 - 
   2*l2*n2*r2*s0*u2*v2 - 2*f0*j2*l2*s2*u2*v2 - 2*d0*j2*m2*s2*u2*v2 + 
   4*c2*j2*n0*s2*u2*v2 + 4*c0*j2*n2*s2*u2*v2 + 2*f0*i2*q2*s2*u2*v2 - 
   2*m2*n0*q2*s2*u2*v2 + 2*d2*i2*r0*s2*u2*v2 - 2*l2*n2*r0*s2*u2*v2 + 
   2*d0*i2*r2*s2*u2*v2 - 2*l2*n0*r2*s2*u2*v2 - 4*c2*i2*s0*s2*u2*v2 + 
   4*l2*m2*s0*s2*u2*v2 - 2*c0*i2*pow(s2,2)*u2*v2 - 
   2*pow(d2,2)*pow(j2,2)*v0*v2 + 2*pow(d2,2)*i2*o2*v0*v2 + 
   2*g2*pow(l2,2)*o2*v0*v2 - 4*d2*l2*n2*o2*v0*v2 - 4*g2*j2*l2*q2*v0*v2 + 
   4*d2*j2*n2*q2*v0*v2 + 2*g2*i2*pow(q2,2)*v0*v2 - 
   2*pow(n2,2)*pow(q2,2)*v0*v2 + 4*d2*j2*l2*s2*v0*v2 - 
   4*d2*i2*q2*s2*v0*v2 + 4*l2*n2*q2*s2*v0*v2 - 
   2*pow(l2,2)*pow(s2,2)*v0*v2 - 2*d0*d2*pow(j2,2)*pow(v2,2) + 
   2*d0*d2*i2*o2*pow(v2,2) + g0*pow(l2,2)*o2*pow(v2,2) - 
   2*d2*l2*n0*o2*pow(v2,2) - 2*d0*l2*n2*o2*pow(v2,2) - 
   2*g0*j2*l2*q2*pow(v2,2) + 2*d2*j2*n0*q2*pow(v2,2) + 
   2*d0*j2*n2*q2*pow(v2,2) + g0*i2*pow(q2,2)*pow(v2,2) - 
   2*n0*n2*pow(q2,2)*pow(v2,2) + 2*d2*j2*l2*s0*pow(v2,2) - 
   2*d2*i2*q2*s0*pow(v2,2) + 2*l2*n2*q2*s0*pow(v2,2) + 
   2*d0*j2*l2*s2*pow(v2,2) - 2*d0*i2*q2*s2*pow(v2,2) + 
   2*l2*n0*q2*s2*pow(v2,2) - 2*pow(l2,2)*s0*s2*pow(v2,2) - 
   2*d2*e2*k2*l2*o2*w0 + 2*c2*f2*k2*l2*o2*w0 + 2*c2*d2*k2*m2*o2*w0 - 
   2*pow(c2,2)*k2*n2*o2*w0 + 2*d2*e2*j2*l2*p2*w0 - 2*c2*f2*j2*l2*p2*w0 - 
   2*c2*d2*j2*m2*p2*w0 + 2*pow(c2,2)*j2*n2*p2*w0 + 2*d2*e2*j2*k2*q2*w0 - 
   2*c2*f2*j2*k2*q2*w0 - 2*d2*e2*i2*p2*q2*w0 + 2*c2*f2*i2*p2*q2*w0 - 
   2*f2*l2*m2*p2*q2*w0 + 2*d2*pow(m2,2)*p2*q2*w0 + 2*e2*l2*n2*p2*q2*w0 - 
   2*c2*m2*n2*p2*q2*w0 + 2*f2*k2*m2*pow(q2,2)*w0 - 
   2*e2*k2*n2*pow(q2,2)*w0 - 2*c2*d2*j2*k2*r2*w0 + 2*c2*d2*i2*p2*r2*w0 + 
   2*f2*pow(l2,2)*p2*r2*w0 - 2*d2*l2*m2*p2*r2*w0 - 2*c2*l2*n2*p2*r2*w0 - 
   2*f2*k2*l2*q2*r2*w0 - 2*d2*k2*m2*q2*r2*w0 + 4*c2*k2*n2*q2*r2*w0 + 
   2*d2*k2*l2*pow(r2,2)*w0 + 2*pow(c2,2)*j2*k2*s2*w0 - 
   2*pow(c2,2)*i2*p2*s2*w0 - 2*e2*pow(l2,2)*p2*s2*w0 + 
   4*c2*l2*m2*p2*s2*w0 + 2*e2*k2*l2*q2*s2*w0 - 2*c2*k2*m2*q2*s2*w0 - 
   2*c2*k2*l2*r2*s2*w0 - 2*d2*e2*pow(j2,2)*u2*w0 + 
   2*c2*f2*pow(j2,2)*u2*w0 + 2*d2*e2*i2*o2*u2*w0 - 2*c2*f2*i2*o2*u2*w0 + 
   2*f2*l2*m2*o2*u2*w0 - 2*d2*pow(m2,2)*o2*u2*w0 - 2*e2*l2*n2*o2*u2*w0 + 
   2*c2*m2*n2*o2*u2*w0 - 2*f2*j2*m2*q2*u2*w0 + 2*e2*j2*n2*q2*u2*w0 - 
   2*f2*j2*l2*r2*u2*w0 + 4*d2*j2*m2*r2*u2*w0 - 2*c2*j2*n2*r2*u2*w0 + 
   2*f2*i2*q2*r2*u2*w0 - 2*m2*n2*q2*r2*u2*w0 - 2*d2*i2*pow(r2,2)*u2*w0 + 
   2*l2*n2*pow(r2,2)*u2*w0 + 2*e2*j2*l2*s2*u2*w0 - 2*c2*j2*m2*s2*u2*w0 - 
   2*e2*i2*q2*s2*u2*w0 + 2*pow(m2,2)*q2*s2*u2*w0 + 2*c2*i2*r2*s2*u2*w0 - 
   2*l2*m2*r2*s2*u2*w0 + 2*c2*d2*pow(j2,2)*v2*w0 - 2*c2*d2*i2*o2*v2*w0 - 
   2*f2*pow(l2,2)*o2*v2*w0 + 2*d2*l2*m2*o2*v2*w0 + 2*c2*l2*n2*o2*v2*w0 + 
   4*f2*j2*l2*q2*v2*w0 - 2*d2*j2*m2*q2*v2*w0 - 2*c2*j2*n2*q2*v2*w0 - 
   2*f2*i2*pow(q2,2)*v2*w0 + 2*m2*n2*pow(q2,2)*v2*w0 - 
   2*d2*j2*l2*r2*v2*w0 + 2*d2*i2*q2*r2*v2*w0 - 2*l2*n2*q2*r2*v2*w0 - 
   2*c2*j2*l2*s2*v2*w0 + 2*c2*i2*q2*s2*v2*w0 - 2*l2*m2*q2*s2*v2*w0 + 
   2*pow(l2,2)*r2*s2*v2*w0 - 2*d2*e0*k2*l2*o2*w2 - 2*d0*e2*k2*l2*o2*w2 + 
   2*c2*f0*k2*l2*o2*w2 + 2*c0*f2*k2*l2*o2*w2 + 2*c2*d0*k2*m2*o2*w2 + 
   2*c0*d2*k2*m2*o2*w2 - 2*pow(c2,2)*k2*n0*o2*w2 - 4*c0*c2*k2*n2*o2*w2 + 
   2*d2*e0*j2*l2*p2*w2 + 2*d0*e2*j2*l2*p2*w2 - 2*c2*f0*j2*l2*p2*w2 - 
   2*c0*f2*j2*l2*p2*w2 - 2*c2*d0*j2*m2*p2*w2 - 2*c0*d2*j2*m2*p2*w2 + 
   2*pow(c2,2)*j2*n0*p2*w2 + 4*c0*c2*j2*n2*p2*w2 + 2*d2*e0*j2*k2*q2*w2 + 
   2*d0*e2*j2*k2*q2*w2 - 2*c2*f0*j2*k2*q2*w2 - 2*c0*f2*j2*k2*q2*w2 - 
   2*d2*e0*i2*p2*q2*w2 - 2*d0*e2*i2*p2*q2*w2 + 2*c2*f0*i2*p2*q2*w2 + 
   2*c0*f2*i2*p2*q2*w2 - 2*f0*l2*m2*p2*q2*w2 + 2*d0*pow(m2,2)*p2*q2*w2 + 
   2*e2*l2*n0*p2*q2*w2 - 2*c2*m2*n0*p2*q2*w2 + 2*e0*l2*n2*p2*q2*w2 - 
   2*c0*m2*n2*p2*q2*w2 + 2*f0*k2*m2*pow(q2,2)*w2 - 
   2*e2*k2*n0*pow(q2,2)*w2 - 2*e0*k2*n2*pow(q2,2)*w2 - 
   2*c2*d2*j2*k2*r0*w2 + 2*c2*d2*i2*p2*r0*w2 + 2*f2*pow(l2,2)*p2*r0*w2 - 
   2*d2*l2*m2*p2*r0*w2 - 2*c2*l2*n2*p2*r0*w2 - 2*f2*k2*l2*q2*r0*w2 - 
   2*d2*k2*m2*q2*r0*w2 + 4*c2*k2*n2*q2*r0*w2 - 2*c2*d0*j2*k2*r2*w2 - 
   2*c0*d2*j2*k2*r2*w2 + 2*c2*d0*i2*p2*r2*w2 + 2*c0*d2*i2*p2*r2*w2 + 
   2*f0*pow(l2,2)*p2*r2*w2 - 2*d0*l2*m2*p2*r2*w2 - 2*c2*l2*n0*p2*r2*w2 - 
   2*c0*l2*n2*p2*r2*w2 - 2*f0*k2*l2*q2*r2*w2 - 2*d0*k2*m2*q2*r2*w2 + 
   4*c2*k2*n0*q2*r2*w2 + 4*c0*k2*n2*q2*r2*w2 + 4*d2*k2*l2*r0*r2*w2 + 
   2*d0*k2*l2*pow(r2,2)*w2 + 2*pow(c2,2)*j2*k2*s0*w2 - 
   2*pow(c2,2)*i2*p2*s0*w2 - 2*e2*pow(l2,2)*p2*s0*w2 + 
   4*c2*l2*m2*p2*s0*w2 + 2*e2*k2*l2*q2*s0*w2 - 2*c2*k2*m2*q2*s0*w2 - 
   2*c2*k2*l2*r2*s0*w2 + 4*c0*c2*j2*k2*s2*w2 - 4*c0*c2*i2*p2*s2*w2 - 
   2*e0*pow(l2,2)*p2*s2*w2 + 4*c0*l2*m2*p2*s2*w2 + 2*e0*k2*l2*q2*s2*w2 - 
   2*c0*k2*m2*q2*s2*w2 - 2*c2*k2*l2*r0*s2*w2 - 2*c0*k2*l2*r2*s2*w2 - 
   2*d2*e2*pow(j2,2)*u0*w2 + 2*c2*f2*pow(j2,2)*u0*w2 + 
   2*d2*e2*i2*o2*u0*w2 - 2*c2*f2*i2*o2*u0*w2 + 2*f2*l2*m2*o2*u0*w2 - 
   2*d2*pow(m2,2)*o2*u0*w2 - 2*e2*l2*n2*o2*u0*w2 + 2*c2*m2*n2*o2*u0*w2 - 
   2*f2*j2*m2*q2*u0*w2 + 2*e2*j2*n2*q2*u0*w2 - 2*f2*j2*l2*r2*u0*w2 + 
   4*d2*j2*m2*r2*u0*w2 - 2*c2*j2*n2*r2*u0*w2 + 2*f2*i2*q2*r2*u0*w2 - 
   2*m2*n2*q2*r2*u0*w2 - 2*d2*i2*pow(r2,2)*u0*w2 + 
   2*l2*n2*pow(r2,2)*u0*w2 + 2*e2*j2*l2*s2*u0*w2 - 2*c2*j2*m2*s2*u0*w2 - 
   2*e2*i2*q2*s2*u0*w2 + 2*pow(m2,2)*q2*s2*u0*w2 + 2*c2*i2*r2*s2*u0*w2 - 
   2*l2*m2*r2*s2*u0*w2 - 2*d2*e0*pow(j2,2)*u2*w2 - 
   2*d0*e2*pow(j2,2)*u2*w2 + 2*c2*f0*pow(j2,2)*u2*w2 + 
   2*c0*f2*pow(j2,2)*u2*w2 + 2*d2*e0*i2*o2*u2*w2 + 2*d0*e2*i2*o2*u2*w2 - 
   2*c2*f0*i2*o2*u2*w2 - 2*c0*f2*i2*o2*u2*w2 + 2*f0*l2*m2*o2*u2*w2 - 
   2*d0*pow(m2,2)*o2*u2*w2 - 2*e2*l2*n0*o2*u2*w2 + 2*c2*m2*n0*o2*u2*w2 - 
   2*e0*l2*n2*o2*u2*w2 + 2*c0*m2*n2*o2*u2*w2 - 2*f0*j2*m2*q2*u2*w2 + 
   2*e2*j2*n0*q2*u2*w2 + 2*e0*j2*n2*q2*u2*w2 - 2*f2*j2*l2*r0*u2*w2 + 
   4*d2*j2*m2*r0*u2*w2 - 2*c2*j2*n2*r0*u2*w2 + 2*f2*i2*q2*r0*u2*w2 - 
   2*m2*n2*q2*r0*u2*w2 - 2*f0*j2*l2*r2*u2*w2 + 4*d0*j2*m2*r2*u2*w2 - 
   2*c2*j2*n0*r2*u2*w2 - 2*c0*j2*n2*r2*u2*w2 + 2*f0*i2*q2*r2*u2*w2 - 
   2*m2*n0*q2*r2*u2*w2 - 4*d2*i2*r0*r2*u2*w2 + 4*l2*n2*r0*r2*u2*w2 - 
   2*d0*i2*pow(r2,2)*u2*w2 + 2*l2*n0*pow(r2,2)*u2*w2 + 
   2*e2*j2*l2*s0*u2*w2 - 2*c2*j2*m2*s0*u2*w2 - 2*e2*i2*q2*s0*u2*w2 + 
   2*pow(m2,2)*q2*s0*u2*w2 + 2*c2*i2*r2*s0*u2*w2 - 2*l2*m2*r2*s0*u2*w2 + 
   2*e0*j2*l2*s2*u2*w2 - 2*c0*j2*m2*s2*u2*w2 - 2*e0*i2*q2*s2*u2*w2 + 
   2*c2*i2*r0*s2*u2*w2 - 2*l2*m2*r0*s2*u2*w2 + 2*c0*i2*r2*s2*u2*w2 + 
   2*c2*d2*pow(j2,2)*v0*w2 - 2*c2*d2*i2*o2*v0*w2 - 
   2*f2*pow(l2,2)*o2*v0*w2 + 2*d2*l2*m2*o2*v0*w2 + 2*c2*l2*n2*o2*v0*w2 + 
   4*f2*j2*l2*q2*v0*w2 - 2*d2*j2*m2*q2*v0*w2 - 2*c2*j2*n2*q2*v0*w2 - 
   2*f2*i2*pow(q2,2)*v0*w2 + 2*m2*n2*pow(q2,2)*v0*w2 - 
   2*d2*j2*l2*r2*v0*w2 + 2*d2*i2*q2*r2*v0*w2 - 2*l2*n2*q2*r2*v0*w2 - 
   2*c2*j2*l2*s2*v0*w2 + 2*c2*i2*q2*s2*v0*w2 - 2*l2*m2*q2*s2*v0*w2 + 
   2*pow(l2,2)*r2*s2*v0*w2 + 2*c2*d0*pow(j2,2)*v2*w2 + 
   2*c0*d2*pow(j2,2)*v2*w2 - 2*c2*d0*i2*o2*v2*w2 - 2*c0*d2*i2*o2*v2*w2 - 
   2*f0*pow(l2,2)*o2*v2*w2 + 2*d0*l2*m2*o2*v2*w2 + 2*c2*l2*n0*o2*v2*w2 + 
   2*c0*l2*n2*o2*v2*w2 + 4*f0*j2*l2*q2*v2*w2 - 2*d0*j2*m2*q2*v2*w2 - 
   2*c2*j2*n0*q2*v2*w2 - 2*c0*j2*n2*q2*v2*w2 - 2*f0*i2*pow(q2,2)*v2*w2 + 
   2*m2*n0*pow(q2,2)*v2*w2 - 2*d2*j2*l2*r0*v2*w2 + 2*d2*i2*q2*r0*v2*w2 - 
   2*l2*n2*q2*r0*v2*w2 - 2*d0*j2*l2*r2*v2*w2 + 2*d0*i2*q2*r2*v2*w2 - 
   2*l2*n0*q2*r2*v2*w2 - 2*c2*j2*l2*s0*v2*w2 + 2*c2*i2*q2*s0*v2*w2 - 
   2*l2*m2*q2*s0*v2*w2 + 2*pow(l2,2)*r2*s0*v2*w2 - 2*c0*j2*l2*s2*v2*w2 + 
   2*c0*i2*q2*s2*v2*w2 + 2*pow(l2,2)*r0*s2*v2*w2 - 
   2*pow(c2,2)*pow(j2,2)*w0*w2 + 2*pow(c2,2)*i2*o2*w0*w2 + 
   2*e2*pow(l2,2)*o2*w0*w2 - 4*c2*l2*m2*o2*w0*w2 - 4*e2*j2*l2*q2*w0*w2 + 
   4*c2*j2*m2*q2*w0*w2 + 2*e2*i2*pow(q2,2)*w0*w2 - 
   2*pow(m2,2)*pow(q2,2)*w0*w2 + 4*c2*j2*l2*r2*w0*w2 - 
   4*c2*i2*q2*r2*w0*w2 + 4*l2*m2*q2*r2*w0*w2 - 
   2*pow(l2,2)*pow(r2,2)*w0*w2 - 2*c0*c2*pow(j2,2)*pow(w2,2) + 
   2*c0*c2*i2*o2*pow(w2,2) + e0*pow(l2,2)*o2*pow(w2,2) - 
   2*c0*l2*m2*o2*pow(w2,2) - 2*e0*j2*l2*q2*pow(w2,2) + 
   2*c0*j2*m2*q2*pow(w2,2) + e0*i2*pow(q2,2)*pow(w2,2) + 
   2*c2*j2*l2*r0*pow(w2,2) - 2*c2*i2*q2*r0*pow(w2,2) + 
   2*l2*m2*q2*r0*pow(w2,2) + 2*c0*j2*l2*r2*pow(w2,2) - 
   2*c0*i2*q2*r2*pow(w2,2) - 2*pow(l2,2)*r0*r2*pow(w2,2) + 
   pow(f2,2)*pow(k2,2)*o2*z0 - e2*g2*pow(k2,2)*o2*z0 - 
   2*pow(f2,2)*j2*k2*p2*z0 + 2*e2*g2*j2*k2*p2*z0 + 
   pow(f2,2)*i2*pow(p2,2)*z0 - e2*g2*i2*pow(p2,2)*z0 + 
   g2*pow(m2,2)*pow(p2,2)*z0 - 2*f2*m2*n2*pow(p2,2)*z0 + 
   e2*pow(n2,2)*pow(p2,2)*z0 - 2*g2*k2*m2*p2*r2*z0 + 
   2*f2*k2*n2*p2*r2*z0 + g2*pow(k2,2)*pow(r2,2)*z0 + 
   2*f2*k2*m2*p2*s2*z0 - 2*e2*k2*n2*p2*s2*z0 - 2*f2*pow(k2,2)*r2*s2*z0 + 
   e2*pow(k2,2)*pow(s2,2)*z0 + pow(f2,2)*pow(j2,2)*t2*z0 - 
   e2*g2*pow(j2,2)*t2*z0 - pow(f2,2)*i2*o2*t2*z0 + e2*g2*i2*o2*t2*z0 - 
   g2*pow(m2,2)*o2*t2*z0 + 2*f2*m2*n2*o2*t2*z0 - e2*pow(n2,2)*o2*t2*z0 + 
   2*g2*j2*m2*r2*t2*z0 - 2*f2*j2*n2*r2*t2*z0 - g2*i2*pow(r2,2)*t2*z0 + 
   pow(n2,2)*pow(r2,2)*t2*z0 - 2*f2*j2*m2*s2*t2*z0 + 
   2*e2*j2*n2*s2*t2*z0 + 2*f2*i2*r2*s2*t2*z0 - 2*m2*n2*r2*s2*t2*z0 - 
   e2*i2*pow(s2,2)*t2*z0 + pow(m2,2)*pow(s2,2)*t2*z0 + 
   2*g2*k2*m2*o2*v2*z0 - 2*f2*k2*n2*o2*v2*z0 - 2*g2*j2*m2*p2*v2*z0 + 
   2*f2*j2*n2*p2*v2*z0 - 2*g2*j2*k2*r2*v2*z0 + 2*g2*i2*p2*r2*v2*z0 - 
   2*pow(n2,2)*p2*r2*v2*z0 + 2*f2*j2*k2*s2*v2*z0 - 2*f2*i2*p2*s2*v2*z0 + 
   2*m2*n2*p2*s2*v2*z0 + 2*k2*n2*r2*s2*v2*z0 - 2*k2*m2*pow(s2,2)*v2*z0 + 
   g2*pow(j2,2)*pow(v2,2)*z0 - g2*i2*o2*pow(v2,2)*z0 + 
   pow(n2,2)*o2*pow(v2,2)*z0 - 2*j2*n2*s2*pow(v2,2)*z0 + 
   i2*pow(s2,2)*pow(v2,2)*z0 - 2*f2*k2*m2*o2*w2*z0 + 
   2*e2*k2*n2*o2*w2*z0 + 2*f2*j2*m2*p2*w2*z0 - 2*e2*j2*n2*p2*w2*z0 + 
   2*f2*j2*k2*r2*w2*z0 - 2*f2*i2*p2*r2*w2*z0 + 2*m2*n2*p2*r2*w2*z0 - 
   2*k2*n2*pow(r2,2)*w2*z0 - 2*e2*j2*k2*s2*w2*z0 + 2*e2*i2*p2*s2*w2*z0 - 
   2*pow(m2,2)*p2*s2*w2*z0 + 2*k2*m2*r2*s2*w2*z0 - 
   2*f2*pow(j2,2)*v2*w2*z0 + 2*f2*i2*o2*v2*w2*z0 - 2*m2*n2*o2*v2*w2*z0 + 
   2*j2*n2*r2*v2*w2*z0 + 2*j2*m2*s2*v2*w2*z0 - 2*i2*r2*s2*v2*w2*z0 + 
   e2*pow(j2,2)*pow(w2,2)*z0 - e2*i2*o2*pow(w2,2)*z0 + 
   pow(m2,2)*o2*pow(w2,2)*z0 - 2*j2*m2*r2*pow(w2,2)*z0 + 
   i2*pow(r2,2)*pow(w2,2)*z0 + 2*f0*f2*pow(k2,2)*o2*z2 - 
   e2*g0*pow(k2,2)*o2*z2 - e0*g2*pow(k2,2)*o2*z2 - 4*f0*f2*j2*k2*p2*z2 + 
   2*e2*g0*j2*k2*p2*z2 + 2*e0*g2*j2*k2*p2*z2 + 2*f0*f2*i2*pow(p2,2)*z2 - 
   e2*g0*i2*pow(p2,2)*z2 - e0*g2*i2*pow(p2,2)*z2 + 
   g0*pow(m2,2)*pow(p2,2)*z2 - 2*f2*m2*n0*pow(p2,2)*z2 - 
   2*f0*m2*n2*pow(p2,2)*z2 + 2*e2*n0*n2*pow(p2,2)*z2 + 
   e0*pow(n2,2)*pow(p2,2)*z2 - 2*g2*k2*m2*p2*r0*z2 + 
   2*f2*k2*n2*p2*r0*z2 - 2*g0*k2*m2*p2*r2*z2 + 2*f2*k2*n0*p2*r2*z2 + 
   2*f0*k2*n2*p2*r2*z2 + 2*g2*pow(k2,2)*r0*r2*z2 + 
   g0*pow(k2,2)*pow(r2,2)*z2 + 2*f2*k2*m2*p2*s0*z2 - 
   2*e2*k2*n2*p2*s0*z2 - 2*f2*pow(k2,2)*r2*s0*z2 + 2*f0*k2*m2*p2*s2*z2 - 
   2*e2*k2*n0*p2*s2*z2 - 2*e0*k2*n2*p2*s2*z2 - 2*f2*pow(k2,2)*r0*s2*z2 - 
   2*f0*pow(k2,2)*r2*s2*z2 + 2*e2*pow(k2,2)*s0*s2*z2 + 
   e0*pow(k2,2)*pow(s2,2)*z2 + 2*f0*f2*pow(j2,2)*t2*z2 - 
   e2*g0*pow(j2,2)*t2*z2 - e0*g2*pow(j2,2)*t2*z2 - 2*f0*f2*i2*o2*t2*z2 + 
   e2*g0*i2*o2*t2*z2 + e0*g2*i2*o2*t2*z2 - g0*pow(m2,2)*o2*t2*z2 + 
   2*f2*m2*n0*o2*t2*z2 + 2*f0*m2*n2*o2*t2*z2 - 2*e2*n0*n2*o2*t2*z2 - 
   e0*pow(n2,2)*o2*t2*z2 + 2*g2*j2*m2*r0*t2*z2 - 2*f2*j2*n2*r0*t2*z2 + 
   2*g0*j2*m2*r2*t2*z2 - 2*f2*j2*n0*r2*t2*z2 - 2*f0*j2*n2*r2*t2*z2 - 
   2*g2*i2*r0*r2*t2*z2 + 2*pow(n2,2)*r0*r2*t2*z2 - 
   g0*i2*pow(r2,2)*t2*z2 + 2*n0*n2*pow(r2,2)*t2*z2 - 
   2*f2*j2*m2*s0*t2*z2 + 2*e2*j2*n2*s0*t2*z2 + 2*f2*i2*r2*s0*t2*z2 - 
   2*m2*n2*r2*s0*t2*z2 - 2*f0*j2*m2*s2*t2*z2 + 2*e2*j2*n0*s2*t2*z2 + 
   2*e0*j2*n2*s2*t2*z2 + 2*f2*i2*r0*s2*t2*z2 - 2*m2*n2*r0*s2*t2*z2 + 
   2*f0*i2*r2*s2*t2*z2 - 2*m2*n0*r2*s2*t2*z2 - 2*e2*i2*s0*s2*t2*z2 + 
   2*pow(m2,2)*s0*s2*t2*z2 - e0*i2*pow(s2,2)*t2*z2 + 
   2*g2*k2*m2*o2*v0*z2 - 2*f2*k2*n2*o2*v0*z2 - 2*g2*j2*m2*p2*v0*z2 + 
   2*f2*j2*n2*p2*v0*z2 - 2*g2*j2*k2*r2*v0*z2 + 2*g2*i2*p2*r2*v0*z2 - 
   2*pow(n2,2)*p2*r2*v0*z2 + 2*f2*j2*k2*s2*v0*z2 - 2*f2*i2*p2*s2*v0*z2 + 
   2*m2*n2*p2*s2*v0*z2 + 2*k2*n2*r2*s2*v0*z2 - 2*k2*m2*pow(s2,2)*v0*z2 + 
   2*g0*k2*m2*o2*v2*z2 - 2*f2*k2*n0*o2*v2*z2 - 2*f0*k2*n2*o2*v2*z2 - 
   2*g0*j2*m2*p2*v2*z2 + 2*f2*j2*n0*p2*v2*z2 + 2*f0*j2*n2*p2*v2*z2 - 
   2*g2*j2*k2*r0*v2*z2 + 2*g2*i2*p2*r0*v2*z2 - 2*pow(n2,2)*p2*r0*v2*z2 - 
   2*g0*j2*k2*r2*v2*z2 + 2*g0*i2*p2*r2*v2*z2 - 4*n0*n2*p2*r2*v2*z2 + 
   2*f2*j2*k2*s0*v2*z2 - 2*f2*i2*p2*s0*v2*z2 + 2*m2*n2*p2*s0*v2*z2 + 
   2*k2*n2*r2*s0*v2*z2 + 2*f0*j2*k2*s2*v2*z2 - 2*f0*i2*p2*s2*v2*z2 + 
   2*m2*n0*p2*s2*v2*z2 + 2*k2*n2*r0*s2*v2*z2 + 2*k2*n0*r2*s2*v2*z2 - 
   4*k2*m2*s0*s2*v2*z2 + 2*g2*pow(j2,2)*v0*v2*z2 - 2*g2*i2*o2*v0*v2*z2 + 
   2*pow(n2,2)*o2*v0*v2*z2 - 4*j2*n2*s2*v0*v2*z2 + 
   2*i2*pow(s2,2)*v0*v2*z2 + g0*pow(j2,2)*pow(v2,2)*z2 - 
   g0*i2*o2*pow(v2,2)*z2 + 2*n0*n2*o2*pow(v2,2)*z2 - 
   2*j2*n2*s0*pow(v2,2)*z2 - 2*j2*n0*s2*pow(v2,2)*z2 + 
   2*i2*s0*s2*pow(v2,2)*z2 - 2*f2*k2*m2*o2*w0*z2 + 2*e2*k2*n2*o2*w0*z2 + 
   2*f2*j2*m2*p2*w0*z2 - 2*e2*j2*n2*p2*w0*z2 + 2*f2*j2*k2*r2*w0*z2 - 
   2*f2*i2*p2*r2*w0*z2 + 2*m2*n2*p2*r2*w0*z2 - 2*k2*n2*pow(r2,2)*w0*z2 - 
   2*e2*j2*k2*s2*w0*z2 + 2*e2*i2*p2*s2*w0*z2 - 2*pow(m2,2)*p2*s2*w0*z2 + 
   2*k2*m2*r2*s2*w0*z2 - 2*f2*pow(j2,2)*v2*w0*z2 + 2*f2*i2*o2*v2*w0*z2 - 
   2*m2*n2*o2*v2*w0*z2 + 2*j2*n2*r2*v2*w0*z2 + 2*j2*m2*s2*v2*w0*z2 - 
   2*i2*r2*s2*v2*w0*z2 - 2*f0*k2*m2*o2*w2*z2 + 2*e2*k2*n0*o2*w2*z2 + 
   2*e0*k2*n2*o2*w2*z2 + 2*f0*j2*m2*p2*w2*z2 - 2*e2*j2*n0*p2*w2*z2 - 
   2*e0*j2*n2*p2*w2*z2 + 2*f2*j2*k2*r0*w2*z2 - 2*f2*i2*p2*r0*w2*z2 + 
   2*m2*n2*p2*r0*w2*z2 + 2*f0*j2*k2*r2*w2*z2 - 2*f0*i2*p2*r2*w2*z2 + 
   2*m2*n0*p2*r2*w2*z2 - 4*k2*n2*r0*r2*w2*z2 - 2*k2*n0*pow(r2,2)*w2*z2 - 
   2*e2*j2*k2*s0*w2*z2 + 2*e2*i2*p2*s0*w2*z2 - 2*pow(m2,2)*p2*s0*w2*z2 + 
   2*k2*m2*r2*s0*w2*z2 - 2*e0*j2*k2*s2*w2*z2 + 2*e0*i2*p2*s2*w2*z2 + 
   2*k2*m2*r0*s2*w2*z2 - 2*f2*pow(j2,2)*v0*w2*z2 + 2*f2*i2*o2*v0*w2*z2 - 
   2*m2*n2*o2*v0*w2*z2 + 2*j2*n2*r2*v0*w2*z2 + 2*j2*m2*s2*v0*w2*z2 - 
   2*i2*r2*s2*v0*w2*z2 - 2*f0*pow(j2,2)*v2*w2*z2 + 2*f0*i2*o2*v2*w2*z2 - 
   2*m2*n0*o2*v2*w2*z2 + 2*j2*n2*r0*v2*w2*z2 + 2*j2*n0*r2*v2*w2*z2 + 
   2*j2*m2*s0*v2*w2*z2 - 2*i2*r2*s0*v2*w2*z2 - 2*i2*r0*s2*v2*w2*z2 + 
   2*e2*pow(j2,2)*w0*w2*z2 - 2*e2*i2*o2*w0*w2*z2 + 
   2*pow(m2,2)*o2*w0*w2*z2 - 4*j2*m2*r2*w0*w2*z2 + 
   2*i2*pow(r2,2)*w0*w2*z2 + e0*pow(j2,2)*pow(w2,2)*z2 - 
   e0*i2*o2*pow(w2,2)*z2 - 2*j2*m2*r0*pow(w2,2)*z2 + 
   2*i2*r0*r2*pow(w2,2)*z2
;
    
    k[2]  = pow(d2,2)*e1*pow(k2,2)*o2 + 2*d1*d2*e2*pow(k2,2)*o2 - 
   2*c2*d2*f1*pow(k2,2)*o2 - 2*c2*d1*f2*pow(k2,2)*o2 - 
   2*c1*d2*f2*pow(k2,2)*o2 + pow(c2,2)*g1*pow(k2,2)*o2 + 
   2*c1*c2*g2*pow(k2,2)*o2 - 2*pow(d2,2)*e1*j2*k2*p2 - 
   4*d1*d2*e2*j2*k2*p2 + 4*c2*d2*f1*j2*k2*p2 + 4*c2*d1*f2*j2*k2*p2 + 
   4*c1*d2*f2*j2*k2*p2 - 2*pow(c2,2)*g1*j2*k2*p2 - 4*c1*c2*g2*j2*k2*p2 + 
   pow(d2,2)*e1*i2*pow(p2,2) + 2*d1*d2*e2*i2*pow(p2,2) - 
   2*c2*d2*f1*i2*pow(p2,2) - 2*c2*d1*f2*i2*pow(p2,2) - 
   2*c1*d2*f2*i2*pow(p2,2) + pow(c2,2)*g1*i2*pow(p2,2) + 
   2*c1*c2*g2*i2*pow(p2,2) - 2*f1*f2*pow(l2,2)*pow(p2,2) + 
   e2*g1*pow(l2,2)*pow(p2,2) + e1*g2*pow(l2,2)*pow(p2,2) + 
   2*d2*f1*l2*m2*pow(p2,2) + 2*d1*f2*l2*m2*pow(p2,2) - 
   2*c2*g1*l2*m2*pow(p2,2) - 2*c1*g2*l2*m2*pow(p2,2) - 
   2*d1*d2*pow(m2,2)*pow(p2,2) - 2*d2*e2*l2*n1*pow(p2,2) + 
   2*c2*f2*l2*n1*pow(p2,2) + 2*c2*d2*m2*n1*pow(p2,2) - 
   2*d2*e1*l2*n2*pow(p2,2) - 2*d1*e2*l2*n2*pow(p2,2) + 
   2*c2*f1*l2*n2*pow(p2,2) + 2*c1*f2*l2*n2*pow(p2,2) + 
   2*c2*d1*m2*n2*pow(p2,2) + 2*c1*d2*m2*n2*pow(p2,2) - 
   2*pow(c2,2)*n1*n2*pow(p2,2) - 2*c1*c2*pow(n2,2)*pow(p2,2) + 
   4*f1*f2*k2*l2*p2*q2 - 2*e2*g1*k2*l2*p2*q2 - 2*e1*g2*k2*l2*p2*q2 - 
   2*d2*f1*k2*m2*p2*q2 - 2*d1*f2*k2*m2*p2*q2 + 2*c2*g1*k2*m2*p2*q2 + 
   2*c1*g2*k2*m2*p2*q2 + 2*d2*e2*k2*n1*p2*q2 - 2*c2*f2*k2*n1*p2*q2 + 
   2*d2*e1*k2*n2*p2*q2 + 2*d1*e2*k2*n2*p2*q2 - 2*c2*f1*k2*n2*p2*q2 - 
   2*c1*f2*k2*n2*p2*q2 - 2*f1*f2*pow(k2,2)*pow(q2,2) + 
   e2*g1*pow(k2,2)*pow(q2,2) + e1*g2*pow(k2,2)*pow(q2,2) - 
   2*d2*f2*k2*l2*p2*r1 + 2*c2*g2*k2*l2*p2*r1 + 2*pow(d2,2)*k2*m2*p2*r1 - 
   2*c2*d2*k2*n2*p2*r1 + 2*d2*f2*pow(k2,2)*q2*r1 - 
   2*c2*g2*pow(k2,2)*q2*r1 - 2*d2*f1*k2*l2*p2*r2 - 2*d1*f2*k2*l2*p2*r2 + 
   2*c2*g1*k2*l2*p2*r2 + 2*c1*g2*k2*l2*p2*r2 + 4*d1*d2*k2*m2*p2*r2 - 
   2*c2*d2*k2*n1*p2*r2 - 2*c2*d1*k2*n2*p2*r2 - 2*c1*d2*k2*n2*p2*r2 + 
   2*d2*f1*pow(k2,2)*q2*r2 + 2*d1*f2*pow(k2,2)*q2*r2 - 
   2*c2*g1*pow(k2,2)*q2*r2 - 2*c1*g2*pow(k2,2)*q2*r2 - 
   2*pow(d2,2)*pow(k2,2)*r1*r2 - 2*d1*d2*pow(k2,2)*pow(r2,2) + 
   2*d2*e2*k2*l2*p2*s1 - 2*c2*f2*k2*l2*p2*s1 - 2*c2*d2*k2*m2*p2*s1 + 
   2*pow(c2,2)*k2*n2*p2*s1 - 2*d2*e2*pow(k2,2)*q2*s1 + 
   2*c2*f2*pow(k2,2)*q2*s1 + 2*c2*d2*pow(k2,2)*r2*s1 + 
   2*d2*e1*k2*l2*p2*s2 + 2*d1*e2*k2*l2*p2*s2 - 2*c2*f1*k2*l2*p2*s2 - 
   2*c1*f2*k2*l2*p2*s2 - 2*c2*d1*k2*m2*p2*s2 - 2*c1*d2*k2*m2*p2*s2 + 
   2*pow(c2,2)*k2*n1*p2*s2 + 4*c1*c2*k2*n2*p2*s2 - 
   2*d2*e1*pow(k2,2)*q2*s2 - 2*d1*e2*pow(k2,2)*q2*s2 + 
   2*c2*f1*pow(k2,2)*q2*s2 + 2*c1*f2*pow(k2,2)*q2*s2 + 
   2*c2*d2*pow(k2,2)*r1*s2 + 2*c2*d1*pow(k2,2)*r2*s2 + 
   2*c1*d2*pow(k2,2)*r2*s2 - 2*pow(c2,2)*pow(k2,2)*s1*s2 - 
   2*c1*c2*pow(k2,2)*pow(s2,2) + pow(d2,2)*e1*pow(j2,2)*t2 + 
   2*d1*d2*e2*pow(j2,2)*t2 - 2*c2*d2*f1*pow(j2,2)*t2 - 
   2*c2*d1*f2*pow(j2,2)*t2 - 2*c1*d2*f2*pow(j2,2)*t2 + 
   pow(c2,2)*g1*pow(j2,2)*t2 + 2*c1*c2*g2*pow(j2,2)*t2 - 
   pow(d2,2)*e1*i2*o2*t2 - 2*d1*d2*e2*i2*o2*t2 + 2*c2*d2*f1*i2*o2*t2 + 
   2*c2*d1*f2*i2*o2*t2 + 2*c1*d2*f2*i2*o2*t2 - pow(c2,2)*g1*i2*o2*t2 - 
   2*c1*c2*g2*i2*o2*t2 + 2*f1*f2*pow(l2,2)*o2*t2 - 
   e2*g1*pow(l2,2)*o2*t2 - e1*g2*pow(l2,2)*o2*t2 - 2*d2*f1*l2*m2*o2*t2 - 
   2*d1*f2*l2*m2*o2*t2 + 2*c2*g1*l2*m2*o2*t2 + 2*c1*g2*l2*m2*o2*t2 + 
   2*d1*d2*pow(m2,2)*o2*t2 + 2*d2*e2*l2*n1*o2*t2 - 2*c2*f2*l2*n1*o2*t2 - 
   2*c2*d2*m2*n1*o2*t2 + 2*d2*e1*l2*n2*o2*t2 + 2*d1*e2*l2*n2*o2*t2 - 
   2*c2*f1*l2*n2*o2*t2 - 2*c1*f2*l2*n2*o2*t2 - 2*c2*d1*m2*n2*o2*t2 - 
   2*c1*d2*m2*n2*o2*t2 + 2*pow(c2,2)*n1*n2*o2*t2 + 
   2*c1*c2*pow(n2,2)*o2*t2 - 4*f1*f2*j2*l2*q2*t2 + 2*e2*g1*j2*l2*q2*t2 + 
   2*e1*g2*j2*l2*q2*t2 + 2*d2*f1*j2*m2*q2*t2 + 2*d1*f2*j2*m2*q2*t2 - 
   2*c2*g1*j2*m2*q2*t2 - 2*c1*g2*j2*m2*q2*t2 - 2*d2*e2*j2*n1*q2*t2 + 
   2*c2*f2*j2*n1*q2*t2 - 2*d2*e1*j2*n2*q2*t2 - 2*d1*e2*j2*n2*q2*t2 + 
   2*c2*f1*j2*n2*q2*t2 + 2*c1*f2*j2*n2*q2*t2 + 2*f1*f2*i2*pow(q2,2)*t2 - 
   e2*g1*i2*pow(q2,2)*t2 - e1*g2*i2*pow(q2,2)*t2 + 
   g1*pow(m2,2)*pow(q2,2)*t2 - 2*f2*m2*n1*pow(q2,2)*t2 - 
   2*f1*m2*n2*pow(q2,2)*t2 + 2*e2*n1*n2*pow(q2,2)*t2 + 
   e1*pow(n2,2)*pow(q2,2)*t2 + 2*d2*f2*j2*l2*r1*t2 - 
   2*c2*g2*j2*l2*r1*t2 - 2*pow(d2,2)*j2*m2*r1*t2 + 2*c2*d2*j2*n2*r1*t2 - 
   2*d2*f2*i2*q2*r1*t2 + 2*c2*g2*i2*q2*r1*t2 - 2*g2*l2*m2*q2*r1*t2 + 
   2*f2*l2*n2*q2*r1*t2 + 2*d2*m2*n2*q2*r1*t2 - 2*c2*pow(n2,2)*q2*r1*t2 + 
   2*d2*f1*j2*l2*r2*t2 + 2*d1*f2*j2*l2*r2*t2 - 2*c2*g1*j2*l2*r2*t2 - 
   2*c1*g2*j2*l2*r2*t2 - 4*d1*d2*j2*m2*r2*t2 + 2*c2*d2*j2*n1*r2*t2 + 
   2*c2*d1*j2*n2*r2*t2 + 2*c1*d2*j2*n2*r2*t2 - 2*d2*f1*i2*q2*r2*t2 - 
   2*d1*f2*i2*q2*r2*t2 + 2*c2*g1*i2*q2*r2*t2 + 2*c1*g2*i2*q2*r2*t2 - 
   2*g1*l2*m2*q2*r2*t2 + 2*f2*l2*n1*q2*r2*t2 + 2*d2*m2*n1*q2*r2*t2 + 
   2*f1*l2*n2*q2*r2*t2 + 2*d1*m2*n2*q2*r2*t2 - 4*c2*n1*n2*q2*r2*t2 - 
   2*c1*pow(n2,2)*q2*r2*t2 + 2*pow(d2,2)*i2*r1*r2*t2 + 
   2*g2*pow(l2,2)*r1*r2*t2 - 4*d2*l2*n2*r1*r2*t2 + 
   2*d1*d2*i2*pow(r2,2)*t2 + g1*pow(l2,2)*pow(r2,2)*t2 - 
   2*d2*l2*n1*pow(r2,2)*t2 - 2*d1*l2*n2*pow(r2,2)*t2 - 
   2*d2*e2*j2*l2*s1*t2 + 2*c2*f2*j2*l2*s1*t2 + 2*c2*d2*j2*m2*s1*t2 - 
   2*pow(c2,2)*j2*n2*s1*t2 + 2*d2*e2*i2*q2*s1*t2 - 2*c2*f2*i2*q2*s1*t2 + 
   2*f2*l2*m2*q2*s1*t2 - 2*d2*pow(m2,2)*q2*s1*t2 - 2*e2*l2*n2*q2*s1*t2 + 
   2*c2*m2*n2*q2*s1*t2 - 2*c2*d2*i2*r2*s1*t2 - 2*f2*pow(l2,2)*r2*s1*t2 + 
   2*d2*l2*m2*r2*s1*t2 + 2*c2*l2*n2*r2*s1*t2 - 2*d2*e1*j2*l2*s2*t2 - 
   2*d1*e2*j2*l2*s2*t2 + 2*c2*f1*j2*l2*s2*t2 + 2*c1*f2*j2*l2*s2*t2 + 
   2*c2*d1*j2*m2*s2*t2 + 2*c1*d2*j2*m2*s2*t2 - 2*pow(c2,2)*j2*n1*s2*t2 - 
   4*c1*c2*j2*n2*s2*t2 + 2*d2*e1*i2*q2*s2*t2 + 2*d1*e2*i2*q2*s2*t2 - 
   2*c2*f1*i2*q2*s2*t2 - 2*c1*f2*i2*q2*s2*t2 + 2*f1*l2*m2*q2*s2*t2 - 
   2*d1*pow(m2,2)*q2*s2*t2 - 2*e2*l2*n1*q2*s2*t2 + 2*c2*m2*n1*q2*s2*t2 - 
   2*e1*l2*n2*q2*s2*t2 + 2*c1*m2*n2*q2*s2*t2 - 2*c2*d2*i2*r1*s2*t2 - 
   2*f2*pow(l2,2)*r1*s2*t2 + 2*d2*l2*m2*r1*s2*t2 + 2*c2*l2*n2*r1*s2*t2 - 
   2*c2*d1*i2*r2*s2*t2 - 2*c1*d2*i2*r2*s2*t2 - 2*f1*pow(l2,2)*r2*s2*t2 + 
   2*d1*l2*m2*r2*s2*t2 + 2*c2*l2*n1*r2*s2*t2 + 2*c1*l2*n2*r2*s2*t2 + 
   2*pow(c2,2)*i2*s1*s2*t2 + 2*e2*pow(l2,2)*s1*s2*t2 - 
   4*c2*l2*m2*s1*s2*t2 + 2*c1*c2*i2*pow(s2,2)*t2 + 
   e1*pow(l2,2)*pow(s2,2)*t2 - 2*c1*l2*m2*pow(s2,2)*t2 - 
   2*pow(f2,2)*k2*l2*o2*u1 + 2*e2*g2*k2*l2*o2*u1 + 2*d2*f2*k2*m2*o2*u1 - 
   2*c2*g2*k2*m2*o2*u1 - 2*d2*e2*k2*n2*o2*u1 + 2*c2*f2*k2*n2*o2*u1 + 
   2*pow(f2,2)*j2*l2*p2*u1 - 2*e2*g2*j2*l2*p2*u1 - 2*d2*f2*j2*m2*p2*u1 + 
   2*c2*g2*j2*m2*p2*u1 + 2*d2*e2*j2*n2*p2*u1 - 2*c2*f2*j2*n2*p2*u1 + 
   2*pow(f2,2)*j2*k2*q2*u1 - 2*e2*g2*j2*k2*q2*u1 - 
   2*pow(f2,2)*i2*p2*q2*u1 + 2*e2*g2*i2*p2*q2*u1 - 
   2*g2*pow(m2,2)*p2*q2*u1 + 4*f2*m2*n2*p2*q2*u1 - 
   2*e2*pow(n2,2)*p2*q2*u1 - 2*d2*f2*j2*k2*r2*u1 + 2*c2*g2*j2*k2*r2*u1 + 
   2*d2*f2*i2*p2*r2*u1 - 2*c2*g2*i2*p2*r2*u1 + 2*g2*l2*m2*p2*r2*u1 - 
   2*f2*l2*n2*p2*r2*u1 - 2*d2*m2*n2*p2*r2*u1 + 2*c2*pow(n2,2)*p2*r2*u1 + 
   2*g2*k2*m2*q2*r2*u1 - 2*f2*k2*n2*q2*r2*u1 - 2*g2*k2*l2*pow(r2,2)*u1 + 
   2*d2*k2*n2*pow(r2,2)*u1 + 2*d2*e2*j2*k2*s2*u1 - 2*c2*f2*j2*k2*s2*u1 - 
   2*d2*e2*i2*p2*s2*u1 + 2*c2*f2*i2*p2*s2*u1 - 2*f2*l2*m2*p2*s2*u1 + 
   2*d2*pow(m2,2)*p2*s2*u1 + 2*e2*l2*n2*p2*s2*u1 - 2*c2*m2*n2*p2*s2*u1 - 
   2*f2*k2*m2*q2*s2*u1 + 2*e2*k2*n2*q2*s2*u1 + 4*f2*k2*l2*r2*s2*u1 - 
   2*d2*k2*m2*r2*s2*u1 - 2*c2*k2*n2*r2*s2*u1 - 2*e2*k2*l2*pow(s2,2)*u1 + 
   2*c2*k2*m2*pow(s2,2)*u1 - 4*f1*f2*k2*l2*o2*u2 + 2*e2*g1*k2*l2*o2*u2 + 
   2*e1*g2*k2*l2*o2*u2 + 2*d2*f1*k2*m2*o2*u2 + 2*d1*f2*k2*m2*o2*u2 - 
   2*c2*g1*k2*m2*o2*u2 - 2*c1*g2*k2*m2*o2*u2 - 2*d2*e2*k2*n1*o2*u2 + 
   2*c2*f2*k2*n1*o2*u2 - 2*d2*e1*k2*n2*o2*u2 - 2*d1*e2*k2*n2*o2*u2 + 
   2*c2*f1*k2*n2*o2*u2 + 2*c1*f2*k2*n2*o2*u2 + 4*f1*f2*j2*l2*p2*u2 - 
   2*e2*g1*j2*l2*p2*u2 - 2*e1*g2*j2*l2*p2*u2 - 2*d2*f1*j2*m2*p2*u2 - 
   2*d1*f2*j2*m2*p2*u2 + 2*c2*g1*j2*m2*p2*u2 + 2*c1*g2*j2*m2*p2*u2 + 
   2*d2*e2*j2*n1*p2*u2 - 2*c2*f2*j2*n1*p2*u2 + 2*d2*e1*j2*n2*p2*u2 + 
   2*d1*e2*j2*n2*p2*u2 - 2*c2*f1*j2*n2*p2*u2 - 2*c1*f2*j2*n2*p2*u2 + 
   4*f1*f2*j2*k2*q2*u2 - 2*e2*g1*j2*k2*q2*u2 - 2*e1*g2*j2*k2*q2*u2 - 
   4*f1*f2*i2*p2*q2*u2 + 2*e2*g1*i2*p2*q2*u2 + 2*e1*g2*i2*p2*q2*u2 - 
   2*g1*pow(m2,2)*p2*q2*u2 + 4*f2*m2*n1*p2*q2*u2 + 4*f1*m2*n2*p2*q2*u2 - 
   4*e2*n1*n2*p2*q2*u2 - 2*e1*pow(n2,2)*p2*q2*u2 - 2*d2*f2*j2*k2*r1*u2 + 
   2*c2*g2*j2*k2*r1*u2 + 2*d2*f2*i2*p2*r1*u2 - 2*c2*g2*i2*p2*r1*u2 + 
   2*g2*l2*m2*p2*r1*u2 - 2*f2*l2*n2*p2*r1*u2 - 2*d2*m2*n2*p2*r1*u2 + 
   2*c2*pow(n2,2)*p2*r1*u2 + 2*g2*k2*m2*q2*r1*u2 - 2*f2*k2*n2*q2*r1*u2 - 
   2*d2*f1*j2*k2*r2*u2 - 2*d1*f2*j2*k2*r2*u2 + 2*c2*g1*j2*k2*r2*u2 + 
   2*c1*g2*j2*k2*r2*u2 + 2*d2*f1*i2*p2*r2*u2 + 2*d1*f2*i2*p2*r2*u2 - 
   2*c2*g1*i2*p2*r2*u2 - 2*c1*g2*i2*p2*r2*u2 + 2*g1*l2*m2*p2*r2*u2 - 
   2*f2*l2*n1*p2*r2*u2 - 2*d2*m2*n1*p2*r2*u2 - 2*f1*l2*n2*p2*r2*u2 - 
   2*d1*m2*n2*p2*r2*u2 + 4*c2*n1*n2*p2*r2*u2 + 2*c1*pow(n2,2)*p2*r2*u2 + 
   2*g1*k2*m2*q2*r2*u2 - 2*f2*k2*n1*q2*r2*u2 - 2*f1*k2*n2*q2*r2*u2 - 
   4*g2*k2*l2*r1*r2*u2 + 4*d2*k2*n2*r1*r2*u2 - 2*g1*k2*l2*pow(r2,2)*u2 + 
   2*d2*k2*n1*pow(r2,2)*u2 + 2*d1*k2*n2*pow(r2,2)*u2 + 
   2*d2*e2*j2*k2*s1*u2 - 2*c2*f2*j2*k2*s1*u2 - 2*d2*e2*i2*p2*s1*u2 + 
   2*c2*f2*i2*p2*s1*u2 - 2*f2*l2*m2*p2*s1*u2 + 2*d2*pow(m2,2)*p2*s1*u2 + 
   2*e2*l2*n2*p2*s1*u2 - 2*c2*m2*n2*p2*s1*u2 - 2*f2*k2*m2*q2*s1*u2 + 
   2*e2*k2*n2*q2*s1*u2 + 4*f2*k2*l2*r2*s1*u2 - 2*d2*k2*m2*r2*s1*u2 - 
   2*c2*k2*n2*r2*s1*u2 + 2*d2*e1*j2*k2*s2*u2 + 2*d1*e2*j2*k2*s2*u2 - 
   2*c2*f1*j2*k2*s2*u2 - 2*c1*f2*j2*k2*s2*u2 - 2*d2*e1*i2*p2*s2*u2 - 
   2*d1*e2*i2*p2*s2*u2 + 2*c2*f1*i2*p2*s2*u2 + 2*c1*f2*i2*p2*s2*u2 - 
   2*f1*l2*m2*p2*s2*u2 + 2*d1*pow(m2,2)*p2*s2*u2 + 2*e2*l2*n1*p2*s2*u2 - 
   2*c2*m2*n1*p2*s2*u2 + 2*e1*l2*n2*p2*s2*u2 - 2*c1*m2*n2*p2*s2*u2 - 
   2*f1*k2*m2*q2*s2*u2 + 2*e2*k2*n1*q2*s2*u2 + 2*e1*k2*n2*q2*s2*u2 + 
   4*f2*k2*l2*r1*s2*u2 - 2*d2*k2*m2*r1*s2*u2 - 2*c2*k2*n2*r1*s2*u2 + 
   4*f1*k2*l2*r2*s2*u2 - 2*d1*k2*m2*r2*s2*u2 - 2*c2*k2*n1*r2*s2*u2 - 
   2*c1*k2*n2*r2*s2*u2 - 4*e2*k2*l2*s1*s2*u2 + 4*c2*k2*m2*s1*s2*u2 - 
   2*e1*k2*l2*pow(s2,2)*u2 + 2*c1*k2*m2*pow(s2,2)*u2 - 
   2*pow(f2,2)*pow(j2,2)*u1*u2 + 2*e2*g2*pow(j2,2)*u1*u2 + 
   2*pow(f2,2)*i2*o2*u1*u2 - 2*e2*g2*i2*o2*u1*u2 + 
   2*g2*pow(m2,2)*o2*u1*u2 - 4*f2*m2*n2*o2*u1*u2 + 
   2*e2*pow(n2,2)*o2*u1*u2 - 4*g2*j2*m2*r2*u1*u2 + 4*f2*j2*n2*r2*u1*u2 + 
   2*g2*i2*pow(r2,2)*u1*u2 - 2*pow(n2,2)*pow(r2,2)*u1*u2 + 
   4*f2*j2*m2*s2*u1*u2 - 4*e2*j2*n2*s2*u1*u2 - 4*f2*i2*r2*s2*u1*u2 + 
   4*m2*n2*r2*s2*u1*u2 + 2*e2*i2*pow(s2,2)*u1*u2 - 
   2*pow(m2,2)*pow(s2,2)*u1*u2 - 2*f1*f2*pow(j2,2)*pow(u2,2) + 
   e2*g1*pow(j2,2)*pow(u2,2) + e1*g2*pow(j2,2)*pow(u2,2) + 
   2*f1*f2*i2*o2*pow(u2,2) - e2*g1*i2*o2*pow(u2,2) - 
   e1*g2*i2*o2*pow(u2,2) + g1*pow(m2,2)*o2*pow(u2,2) - 
   2*f2*m2*n1*o2*pow(u2,2) - 2*f1*m2*n2*o2*pow(u2,2) + 
   2*e2*n1*n2*o2*pow(u2,2) + e1*pow(n2,2)*o2*pow(u2,2) - 
   2*g2*j2*m2*r1*pow(u2,2) + 2*f2*j2*n2*r1*pow(u2,2) - 
   2*g1*j2*m2*r2*pow(u2,2) + 2*f2*j2*n1*r2*pow(u2,2) + 
   2*f1*j2*n2*r2*pow(u2,2) + 2*g2*i2*r1*r2*pow(u2,2) - 
   2*pow(n2,2)*r1*r2*pow(u2,2) + g1*i2*pow(r2,2)*pow(u2,2) - 
   2*n1*n2*pow(r2,2)*pow(u2,2) + 2*f2*j2*m2*s1*pow(u2,2) - 
   2*e2*j2*n2*s1*pow(u2,2) - 2*f2*i2*r2*s1*pow(u2,2) + 
   2*m2*n2*r2*s1*pow(u2,2) + 2*f1*j2*m2*s2*pow(u2,2) - 
   2*e2*j2*n1*s2*pow(u2,2) - 2*e1*j2*n2*s2*pow(u2,2) - 
   2*f2*i2*r1*s2*pow(u2,2) + 2*m2*n2*r1*s2*pow(u2,2) - 
   2*f1*i2*r2*s2*pow(u2,2) + 2*m2*n1*r2*s2*pow(u2,2) + 
   2*e2*i2*s1*s2*pow(u2,2) - 2*pow(m2,2)*s1*s2*pow(u2,2) + 
   e1*i2*pow(s2,2)*pow(u2,2) + 2*d2*f2*k2*l2*o2*v1 - 
   2*c2*g2*k2*l2*o2*v1 - 2*pow(d2,2)*k2*m2*o2*v1 + 2*c2*d2*k2*n2*o2*v1 - 
   2*d2*f2*j2*l2*p2*v1 + 2*c2*g2*j2*l2*p2*v1 + 2*pow(d2,2)*j2*m2*p2*v1 - 
   2*c2*d2*j2*n2*p2*v1 - 2*d2*f2*j2*k2*q2*v1 + 2*c2*g2*j2*k2*q2*v1 + 
   2*d2*f2*i2*p2*q2*v1 - 2*c2*g2*i2*p2*q2*v1 + 2*g2*l2*m2*p2*q2*v1 - 
   2*f2*l2*n2*p2*q2*v1 - 2*d2*m2*n2*p2*q2*v1 + 2*c2*pow(n2,2)*p2*q2*v1 - 
   2*g2*k2*m2*pow(q2,2)*v1 + 2*f2*k2*n2*pow(q2,2)*v1 + 
   2*pow(d2,2)*j2*k2*r2*v1 - 2*pow(d2,2)*i2*p2*r2*v1 - 
   2*g2*pow(l2,2)*p2*r2*v1 + 4*d2*l2*n2*p2*r2*v1 + 2*g2*k2*l2*q2*r2*v1 - 
   2*d2*k2*n2*q2*r2*v1 - 2*c2*d2*j2*k2*s2*v1 + 2*c2*d2*i2*p2*s2*v1 + 
   2*f2*pow(l2,2)*p2*s2*v1 - 2*d2*l2*m2*p2*s2*v1 - 2*c2*l2*n2*p2*s2*v1 - 
   2*f2*k2*l2*q2*s2*v1 + 4*d2*k2*m2*q2*s2*v1 - 2*c2*k2*n2*q2*s2*v1 - 
   2*d2*k2*l2*r2*s2*v1 + 2*c2*k2*l2*pow(s2,2)*v1 + 
   2*d2*f2*pow(j2,2)*u2*v1 - 2*c2*g2*pow(j2,2)*u2*v1 - 
   2*d2*f2*i2*o2*u2*v1 + 2*c2*g2*i2*o2*u2*v1 - 2*g2*l2*m2*o2*u2*v1 + 
   2*f2*l2*n2*o2*u2*v1 + 2*d2*m2*n2*o2*u2*v1 - 2*c2*pow(n2,2)*o2*u2*v1 + 
   2*g2*j2*m2*q2*u2*v1 - 2*f2*j2*n2*q2*u2*v1 + 2*g2*j2*l2*r2*u2*v1 - 
   2*d2*j2*n2*r2*u2*v1 - 2*g2*i2*q2*r2*u2*v1 + 2*pow(n2,2)*q2*r2*u2*v1 - 
   2*f2*j2*l2*s2*u2*v1 - 2*d2*j2*m2*s2*u2*v1 + 4*c2*j2*n2*s2*u2*v1 + 
   2*f2*i2*q2*s2*u2*v1 - 2*m2*n2*q2*s2*u2*v1 + 2*d2*i2*r2*s2*u2*v1 - 
   2*l2*n2*r2*s2*u2*v1 - 2*c2*i2*pow(s2,2)*u2*v1 + 
   2*l2*m2*pow(s2,2)*u2*v1 + 2*d2*f1*k2*l2*o2*v2 + 2*d1*f2*k2*l2*o2*v2 - 
   2*c2*g1*k2*l2*o2*v2 - 2*c1*g2*k2*l2*o2*v2 - 4*d1*d2*k2*m2*o2*v2 + 
   2*c2*d2*k2*n1*o2*v2 + 2*c2*d1*k2*n2*o2*v2 + 2*c1*d2*k2*n2*o2*v2 - 
   2*d2*f1*j2*l2*p2*v2 - 2*d1*f2*j2*l2*p2*v2 + 2*c2*g1*j2*l2*p2*v2 + 
   2*c1*g2*j2*l2*p2*v2 + 4*d1*d2*j2*m2*p2*v2 - 2*c2*d2*j2*n1*p2*v2 - 
   2*c2*d1*j2*n2*p2*v2 - 2*c1*d2*j2*n2*p2*v2 - 2*d2*f1*j2*k2*q2*v2 - 
   2*d1*f2*j2*k2*q2*v2 + 2*c2*g1*j2*k2*q2*v2 + 2*c1*g2*j2*k2*q2*v2 + 
   2*d2*f1*i2*p2*q2*v2 + 2*d1*f2*i2*p2*q2*v2 - 2*c2*g1*i2*p2*q2*v2 - 
   2*c1*g2*i2*p2*q2*v2 + 2*g1*l2*m2*p2*q2*v2 - 2*f2*l2*n1*p2*q2*v2 - 
   2*d2*m2*n1*p2*q2*v2 - 2*f1*l2*n2*p2*q2*v2 - 2*d1*m2*n2*p2*q2*v2 + 
   4*c2*n1*n2*p2*q2*v2 + 2*c1*pow(n2,2)*p2*q2*v2 - 
   2*g1*k2*m2*pow(q2,2)*v2 + 2*f2*k2*n1*pow(q2,2)*v2 + 
   2*f1*k2*n2*pow(q2,2)*v2 + 2*pow(d2,2)*j2*k2*r1*v2 - 
   2*pow(d2,2)*i2*p2*r1*v2 - 2*g2*pow(l2,2)*p2*r1*v2 + 
   4*d2*l2*n2*p2*r1*v2 + 2*g2*k2*l2*q2*r1*v2 - 2*d2*k2*n2*q2*r1*v2 + 
   4*d1*d2*j2*k2*r2*v2 - 4*d1*d2*i2*p2*r2*v2 - 2*g1*pow(l2,2)*p2*r2*v2 + 
   4*d2*l2*n1*p2*r2*v2 + 4*d1*l2*n2*p2*r2*v2 + 2*g1*k2*l2*q2*r2*v2 - 
   2*d2*k2*n1*q2*r2*v2 - 2*d1*k2*n2*q2*r2*v2 - 2*c2*d2*j2*k2*s1*v2 + 
   2*c2*d2*i2*p2*s1*v2 + 2*f2*pow(l2,2)*p2*s1*v2 - 2*d2*l2*m2*p2*s1*v2 - 
   2*c2*l2*n2*p2*s1*v2 - 2*f2*k2*l2*q2*s1*v2 + 4*d2*k2*m2*q2*s1*v2 - 
   2*c2*k2*n2*q2*s1*v2 - 2*d2*k2*l2*r2*s1*v2 - 2*c2*d1*j2*k2*s2*v2 - 
   2*c1*d2*j2*k2*s2*v2 + 2*c2*d1*i2*p2*s2*v2 + 2*c1*d2*i2*p2*s2*v2 + 
   2*f1*pow(l2,2)*p2*s2*v2 - 2*d1*l2*m2*p2*s2*v2 - 2*c2*l2*n1*p2*s2*v2 - 
   2*c1*l2*n2*p2*s2*v2 - 2*f1*k2*l2*q2*s2*v2 + 4*d1*k2*m2*q2*s2*v2 - 
   2*c2*k2*n1*q2*s2*v2 - 2*c1*k2*n2*q2*s2*v2 - 2*d2*k2*l2*r1*s2*v2 - 
   2*d1*k2*l2*r2*s2*v2 + 4*c2*k2*l2*s1*s2*v2 + 2*c1*k2*l2*pow(s2,2)*v2 + 
   2*d2*f2*pow(j2,2)*u1*v2 - 2*c2*g2*pow(j2,2)*u1*v2 - 
   2*d2*f2*i2*o2*u1*v2 + 2*c2*g2*i2*o2*u1*v2 - 2*g2*l2*m2*o2*u1*v2 + 
   2*f2*l2*n2*o2*u1*v2 + 2*d2*m2*n2*o2*u1*v2 - 2*c2*pow(n2,2)*o2*u1*v2 + 
   2*g2*j2*m2*q2*u1*v2 - 2*f2*j2*n2*q2*u1*v2 + 2*g2*j2*l2*r2*u1*v2 - 
   2*d2*j2*n2*r2*u1*v2 - 2*g2*i2*q2*r2*u1*v2 + 2*pow(n2,2)*q2*r2*u1*v2 - 
   2*f2*j2*l2*s2*u1*v2 - 2*d2*j2*m2*s2*u1*v2 + 4*c2*j2*n2*s2*u1*v2 + 
   2*f2*i2*q2*s2*u1*v2 - 2*m2*n2*q2*s2*u1*v2 + 2*d2*i2*r2*s2*u1*v2 - 
   2*l2*n2*r2*s2*u1*v2 - 2*c2*i2*pow(s2,2)*u1*v2 + 
   2*l2*m2*pow(s2,2)*u1*v2 + 2*d2*f1*pow(j2,2)*u2*v2 + 
   2*d1*f2*pow(j2,2)*u2*v2 - 2*c2*g1*pow(j2,2)*u2*v2 - 
   2*c1*g2*pow(j2,2)*u2*v2 - 2*d2*f1*i2*o2*u2*v2 - 2*d1*f2*i2*o2*u2*v2 + 
   2*c2*g1*i2*o2*u2*v2 + 2*c1*g2*i2*o2*u2*v2 - 2*g1*l2*m2*o2*u2*v2 + 
   2*f2*l2*n1*o2*u2*v2 + 2*d2*m2*n1*o2*u2*v2 + 2*f1*l2*n2*o2*u2*v2 + 
   2*d1*m2*n2*o2*u2*v2 - 4*c2*n1*n2*o2*u2*v2 - 2*c1*pow(n2,2)*o2*u2*v2 + 
   2*g1*j2*m2*q2*u2*v2 - 2*f2*j2*n1*q2*u2*v2 - 2*f1*j2*n2*q2*u2*v2 + 
   2*g2*j2*l2*r1*u2*v2 - 2*d2*j2*n2*r1*u2*v2 - 2*g2*i2*q2*r1*u2*v2 + 
   2*pow(n2,2)*q2*r1*u2*v2 + 2*g1*j2*l2*r2*u2*v2 - 2*d2*j2*n1*r2*u2*v2 - 
   2*d1*j2*n2*r2*u2*v2 - 2*g1*i2*q2*r2*u2*v2 + 4*n1*n2*q2*r2*u2*v2 - 
   2*f2*j2*l2*s1*u2*v2 - 2*d2*j2*m2*s1*u2*v2 + 4*c2*j2*n2*s1*u2*v2 + 
   2*f2*i2*q2*s1*u2*v2 - 2*m2*n2*q2*s1*u2*v2 + 2*d2*i2*r2*s1*u2*v2 - 
   2*l2*n2*r2*s1*u2*v2 - 2*f1*j2*l2*s2*u2*v2 - 2*d1*j2*m2*s2*u2*v2 + 
   4*c2*j2*n1*s2*u2*v2 + 4*c1*j2*n2*s2*u2*v2 + 2*f1*i2*q2*s2*u2*v2 - 
   2*m2*n1*q2*s2*u2*v2 + 2*d2*i2*r1*s2*u2*v2 - 2*l2*n2*r1*s2*u2*v2 + 
   2*d1*i2*r2*s2*u2*v2 - 2*l2*n1*r2*s2*u2*v2 - 4*c2*i2*s1*s2*u2*v2 + 
   4*l2*m2*s1*s2*u2*v2 - 2*c1*i2*pow(s2,2)*u2*v2 - 
   2*pow(d2,2)*pow(j2,2)*v1*v2 + 2*pow(d2,2)*i2*o2*v1*v2 + 
   2*g2*pow(l2,2)*o2*v1*v2 - 4*d2*l2*n2*o2*v1*v2 - 4*g2*j2*l2*q2*v1*v2 + 
   4*d2*j2*n2*q2*v1*v2 + 2*g2*i2*pow(q2,2)*v1*v2 - 
   2*pow(n2,2)*pow(q2,2)*v1*v2 + 4*d2*j2*l2*s2*v1*v2 - 
   4*d2*i2*q2*s2*v1*v2 + 4*l2*n2*q2*s2*v1*v2 - 
   2*pow(l2,2)*pow(s2,2)*v1*v2 - 2*d1*d2*pow(j2,2)*pow(v2,2) + 
   2*d1*d2*i2*o2*pow(v2,2) + g1*pow(l2,2)*o2*pow(v2,2) - 
   2*d2*l2*n1*o2*pow(v2,2) - 2*d1*l2*n2*o2*pow(v2,2) - 
   2*g1*j2*l2*q2*pow(v2,2) + 2*d2*j2*n1*q2*pow(v2,2) + 
   2*d1*j2*n2*q2*pow(v2,2) + g1*i2*pow(q2,2)*pow(v2,2) - 
   2*n1*n2*pow(q2,2)*pow(v2,2) + 2*d2*j2*l2*s1*pow(v2,2) - 
   2*d2*i2*q2*s1*pow(v2,2) + 2*l2*n2*q2*s1*pow(v2,2) + 
   2*d1*j2*l2*s2*pow(v2,2) - 2*d1*i2*q2*s2*pow(v2,2) + 
   2*l2*n1*q2*s2*pow(v2,2) - 2*pow(l2,2)*s1*s2*pow(v2,2) - 
   2*d2*e2*k2*l2*o2*w1 + 2*c2*f2*k2*l2*o2*w1 + 2*c2*d2*k2*m2*o2*w1 - 
   2*pow(c2,2)*k2*n2*o2*w1 + 2*d2*e2*j2*l2*p2*w1 - 2*c2*f2*j2*l2*p2*w1 - 
   2*c2*d2*j2*m2*p2*w1 + 2*pow(c2,2)*j2*n2*p2*w1 + 2*d2*e2*j2*k2*q2*w1 - 
   2*c2*f2*j2*k2*q2*w1 - 2*d2*e2*i2*p2*q2*w1 + 2*c2*f2*i2*p2*q2*w1 - 
   2*f2*l2*m2*p2*q2*w1 + 2*d2*pow(m2,2)*p2*q2*w1 + 2*e2*l2*n2*p2*q2*w1 - 
   2*c2*m2*n2*p2*q2*w1 + 2*f2*k2*m2*pow(q2,2)*w1 - 
   2*e2*k2*n2*pow(q2,2)*w1 - 2*c2*d2*j2*k2*r2*w1 + 2*c2*d2*i2*p2*r2*w1 + 
   2*f2*pow(l2,2)*p2*r2*w1 - 2*d2*l2*m2*p2*r2*w1 - 2*c2*l2*n2*p2*r2*w1 - 
   2*f2*k2*l2*q2*r2*w1 - 2*d2*k2*m2*q2*r2*w1 + 4*c2*k2*n2*q2*r2*w1 + 
   2*d2*k2*l2*pow(r2,2)*w1 + 2*pow(c2,2)*j2*k2*s2*w1 - 
   2*pow(c2,2)*i2*p2*s2*w1 - 2*e2*pow(l2,2)*p2*s2*w1 + 
   4*c2*l2*m2*p2*s2*w1 + 2*e2*k2*l2*q2*s2*w1 - 2*c2*k2*m2*q2*s2*w1 - 
   2*c2*k2*l2*r2*s2*w1 - 2*d2*e2*pow(j2,2)*u2*w1 + 
   2*c2*f2*pow(j2,2)*u2*w1 + 2*d2*e2*i2*o2*u2*w1 - 2*c2*f2*i2*o2*u2*w1 + 
   2*f2*l2*m2*o2*u2*w1 - 2*d2*pow(m2,2)*o2*u2*w1 - 2*e2*l2*n2*o2*u2*w1 + 
   2*c2*m2*n2*o2*u2*w1 - 2*f2*j2*m2*q2*u2*w1 + 2*e2*j2*n2*q2*u2*w1 - 
   2*f2*j2*l2*r2*u2*w1 + 4*d2*j2*m2*r2*u2*w1 - 2*c2*j2*n2*r2*u2*w1 + 
   2*f2*i2*q2*r2*u2*w1 - 2*m2*n2*q2*r2*u2*w1 - 2*d2*i2*pow(r2,2)*u2*w1 + 
   2*l2*n2*pow(r2,2)*u2*w1 + 2*e2*j2*l2*s2*u2*w1 - 2*c2*j2*m2*s2*u2*w1 - 
   2*e2*i2*q2*s2*u2*w1 + 2*pow(m2,2)*q2*s2*u2*w1 + 2*c2*i2*r2*s2*u2*w1 - 
   2*l2*m2*r2*s2*u2*w1 + 2*c2*d2*pow(j2,2)*v2*w1 - 2*c2*d2*i2*o2*v2*w1 - 
   2*f2*pow(l2,2)*o2*v2*w1 + 2*d2*l2*m2*o2*v2*w1 + 2*c2*l2*n2*o2*v2*w1 + 
   4*f2*j2*l2*q2*v2*w1 - 2*d2*j2*m2*q2*v2*w1 - 2*c2*j2*n2*q2*v2*w1 - 
   2*f2*i2*pow(q2,2)*v2*w1 + 2*m2*n2*pow(q2,2)*v2*w1 - 
   2*d2*j2*l2*r2*v2*w1 + 2*d2*i2*q2*r2*v2*w1 - 2*l2*n2*q2*r2*v2*w1 - 
   2*c2*j2*l2*s2*v2*w1 + 2*c2*i2*q2*s2*v2*w1 - 2*l2*m2*q2*s2*v2*w1 + 
   2*pow(l2,2)*r2*s2*v2*w1 - 2*d2*e1*k2*l2*o2*w2 - 2*d1*e2*k2*l2*o2*w2 + 
   2*c2*f1*k2*l2*o2*w2 + 2*c1*f2*k2*l2*o2*w2 + 2*c2*d1*k2*m2*o2*w2 + 
   2*c1*d2*k2*m2*o2*w2 - 2*pow(c2,2)*k2*n1*o2*w2 - 4*c1*c2*k2*n2*o2*w2 + 
   2*d2*e1*j2*l2*p2*w2 + 2*d1*e2*j2*l2*p2*w2 - 2*c2*f1*j2*l2*p2*w2 - 
   2*c1*f2*j2*l2*p2*w2 - 2*c2*d1*j2*m2*p2*w2 - 2*c1*d2*j2*m2*p2*w2 + 
   2*pow(c2,2)*j2*n1*p2*w2 + 4*c1*c2*j2*n2*p2*w2 + 2*d2*e1*j2*k2*q2*w2 + 
   2*d1*e2*j2*k2*q2*w2 - 2*c2*f1*j2*k2*q2*w2 - 2*c1*f2*j2*k2*q2*w2 - 
   2*d2*e1*i2*p2*q2*w2 - 2*d1*e2*i2*p2*q2*w2 + 2*c2*f1*i2*p2*q2*w2 + 
   2*c1*f2*i2*p2*q2*w2 - 2*f1*l2*m2*p2*q2*w2 + 2*d1*pow(m2,2)*p2*q2*w2 + 
   2*e2*l2*n1*p2*q2*w2 - 2*c2*m2*n1*p2*q2*w2 + 2*e1*l2*n2*p2*q2*w2 - 
   2*c1*m2*n2*p2*q2*w2 + 2*f1*k2*m2*pow(q2,2)*w2 - 
   2*e2*k2*n1*pow(q2,2)*w2 - 2*e1*k2*n2*pow(q2,2)*w2 - 
   2*c2*d2*j2*k2*r1*w2 + 2*c2*d2*i2*p2*r1*w2 + 2*f2*pow(l2,2)*p2*r1*w2 - 
   2*d2*l2*m2*p2*r1*w2 - 2*c2*l2*n2*p2*r1*w2 - 2*f2*k2*l2*q2*r1*w2 - 
   2*d2*k2*m2*q2*r1*w2 + 4*c2*k2*n2*q2*r1*w2 - 2*c2*d1*j2*k2*r2*w2 - 
   2*c1*d2*j2*k2*r2*w2 + 2*c2*d1*i2*p2*r2*w2 + 2*c1*d2*i2*p2*r2*w2 + 
   2*f1*pow(l2,2)*p2*r2*w2 - 2*d1*l2*m2*p2*r2*w2 - 2*c2*l2*n1*p2*r2*w2 - 
   2*c1*l2*n2*p2*r2*w2 - 2*f1*k2*l2*q2*r2*w2 - 2*d1*k2*m2*q2*r2*w2 + 
   4*c2*k2*n1*q2*r2*w2 + 4*c1*k2*n2*q2*r2*w2 + 4*d2*k2*l2*r1*r2*w2 + 
   2*d1*k2*l2*pow(r2,2)*w2 + 2*pow(c2,2)*j2*k2*s1*w2 - 
   2*pow(c2,2)*i2*p2*s1*w2 - 2*e2*pow(l2,2)*p2*s1*w2 + 
   4*c2*l2*m2*p2*s1*w2 + 2*e2*k2*l2*q2*s1*w2 - 2*c2*k2*m2*q2*s1*w2 - 
   2*c2*k2*l2*r2*s1*w2 + 4*c1*c2*j2*k2*s2*w2 - 4*c1*c2*i2*p2*s2*w2 - 
   2*e1*pow(l2,2)*p2*s2*w2 + 4*c1*l2*m2*p2*s2*w2 + 2*e1*k2*l2*q2*s2*w2 - 
   2*c1*k2*m2*q2*s2*w2 - 2*c2*k2*l2*r1*s2*w2 - 2*c1*k2*l2*r2*s2*w2 - 
   2*d2*e2*pow(j2,2)*u1*w2 + 2*c2*f2*pow(j2,2)*u1*w2 + 
   2*d2*e2*i2*o2*u1*w2 - 2*c2*f2*i2*o2*u1*w2 + 2*f2*l2*m2*o2*u1*w2 - 
   2*d2*pow(m2,2)*o2*u1*w2 - 2*e2*l2*n2*o2*u1*w2 + 2*c2*m2*n2*o2*u1*w2 - 
   2*f2*j2*m2*q2*u1*w2 + 2*e2*j2*n2*q2*u1*w2 - 2*f2*j2*l2*r2*u1*w2 + 
   4*d2*j2*m2*r2*u1*w2 - 2*c2*j2*n2*r2*u1*w2 + 2*f2*i2*q2*r2*u1*w2 - 
   2*m2*n2*q2*r2*u1*w2 - 2*d2*i2*pow(r2,2)*u1*w2 + 
   2*l2*n2*pow(r2,2)*u1*w2 + 2*e2*j2*l2*s2*u1*w2 - 2*c2*j2*m2*s2*u1*w2 - 
   2*e2*i2*q2*s2*u1*w2 + 2*pow(m2,2)*q2*s2*u1*w2 + 2*c2*i2*r2*s2*u1*w2 - 
   2*l2*m2*r2*s2*u1*w2 - 2*d2*e1*pow(j2,2)*u2*w2 - 
   2*d1*e2*pow(j2,2)*u2*w2 + 2*c2*f1*pow(j2,2)*u2*w2 + 
   2*c1*f2*pow(j2,2)*u2*w2 + 2*d2*e1*i2*o2*u2*w2 + 2*d1*e2*i2*o2*u2*w2 - 
   2*c2*f1*i2*o2*u2*w2 - 2*c1*f2*i2*o2*u2*w2 + 2*f1*l2*m2*o2*u2*w2 - 
   2*d1*pow(m2,2)*o2*u2*w2 - 2*e2*l2*n1*o2*u2*w2 + 2*c2*m2*n1*o2*u2*w2 - 
   2*e1*l2*n2*o2*u2*w2 + 2*c1*m2*n2*o2*u2*w2 - 2*f1*j2*m2*q2*u2*w2 + 
   2*e2*j2*n1*q2*u2*w2 + 2*e1*j2*n2*q2*u2*w2 - 2*f2*j2*l2*r1*u2*w2 + 
   4*d2*j2*m2*r1*u2*w2 - 2*c2*j2*n2*r1*u2*w2 + 2*f2*i2*q2*r1*u2*w2 - 
   2*m2*n2*q2*r1*u2*w2 - 2*f1*j2*l2*r2*u2*w2 + 4*d1*j2*m2*r2*u2*w2 - 
   2*c2*j2*n1*r2*u2*w2 - 2*c1*j2*n2*r2*u2*w2 + 2*f1*i2*q2*r2*u2*w2 - 
   2*m2*n1*q2*r2*u2*w2 - 4*d2*i2*r1*r2*u2*w2 + 4*l2*n2*r1*r2*u2*w2 - 
   2*d1*i2*pow(r2,2)*u2*w2 + 2*l2*n1*pow(r2,2)*u2*w2 + 
   2*e2*j2*l2*s1*u2*w2 - 2*c2*j2*m2*s1*u2*w2 - 2*e2*i2*q2*s1*u2*w2 + 
   2*pow(m2,2)*q2*s1*u2*w2 + 2*c2*i2*r2*s1*u2*w2 - 2*l2*m2*r2*s1*u2*w2 + 
   2*e1*j2*l2*s2*u2*w2 - 2*c1*j2*m2*s2*u2*w2 - 2*e1*i2*q2*s2*u2*w2 + 
   2*c2*i2*r1*s2*u2*w2 - 2*l2*m2*r1*s2*u2*w2 + 2*c1*i2*r2*s2*u2*w2 + 
   2*c2*d2*pow(j2,2)*v1*w2 - 2*c2*d2*i2*o2*v1*w2 - 
   2*f2*pow(l2,2)*o2*v1*w2 + 2*d2*l2*m2*o2*v1*w2 + 2*c2*l2*n2*o2*v1*w2 + 
   4*f2*j2*l2*q2*v1*w2 - 2*d2*j2*m2*q2*v1*w2 - 2*c2*j2*n2*q2*v1*w2 - 
   2*f2*i2*pow(q2,2)*v1*w2 + 2*m2*n2*pow(q2,2)*v1*w2 - 
   2*d2*j2*l2*r2*v1*w2 + 2*d2*i2*q2*r2*v1*w2 - 2*l2*n2*q2*r2*v1*w2 - 
   2*c2*j2*l2*s2*v1*w2 + 2*c2*i2*q2*s2*v1*w2 - 2*l2*m2*q2*s2*v1*w2 + 
   2*pow(l2,2)*r2*s2*v1*w2 + 2*c2*d1*pow(j2,2)*v2*w2 + 
   2*c1*d2*pow(j2,2)*v2*w2 - 2*c2*d1*i2*o2*v2*w2 - 2*c1*d2*i2*o2*v2*w2 - 
   2*f1*pow(l2,2)*o2*v2*w2 + 2*d1*l2*m2*o2*v2*w2 + 2*c2*l2*n1*o2*v2*w2 + 
   2*c1*l2*n2*o2*v2*w2 + 4*f1*j2*l2*q2*v2*w2 - 2*d1*j2*m2*q2*v2*w2 - 
   2*c2*j2*n1*q2*v2*w2 - 2*c1*j2*n2*q2*v2*w2 - 2*f1*i2*pow(q2,2)*v2*w2 + 
   2*m2*n1*pow(q2,2)*v2*w2 - 2*d2*j2*l2*r1*v2*w2 + 2*d2*i2*q2*r1*v2*w2 - 
   2*l2*n2*q2*r1*v2*w2 - 2*d1*j2*l2*r2*v2*w2 + 2*d1*i2*q2*r2*v2*w2 - 
   2*l2*n1*q2*r2*v2*w2 - 2*c2*j2*l2*s1*v2*w2 + 2*c2*i2*q2*s1*v2*w2 - 
   2*l2*m2*q2*s1*v2*w2 + 2*pow(l2,2)*r2*s1*v2*w2 - 2*c1*j2*l2*s2*v2*w2 + 
   2*c1*i2*q2*s2*v2*w2 + 2*pow(l2,2)*r1*s2*v2*w2 - 
   2*pow(c2,2)*pow(j2,2)*w1*w2 + 2*pow(c2,2)*i2*o2*w1*w2 + 
   2*e2*pow(l2,2)*o2*w1*w2 - 4*c2*l2*m2*o2*w1*w2 - 4*e2*j2*l2*q2*w1*w2 + 
   4*c2*j2*m2*q2*w1*w2 + 2*e2*i2*pow(q2,2)*w1*w2 - 
   2*pow(m2,2)*pow(q2,2)*w1*w2 + 4*c2*j2*l2*r2*w1*w2 - 
   4*c2*i2*q2*r2*w1*w2 + 4*l2*m2*q2*r2*w1*w2 - 
   2*pow(l2,2)*pow(r2,2)*w1*w2 - 2*c1*c2*pow(j2,2)*pow(w2,2) + 
   2*c1*c2*i2*o2*pow(w2,2) + e1*pow(l2,2)*o2*pow(w2,2) - 
   2*c1*l2*m2*o2*pow(w2,2) - 2*e1*j2*l2*q2*pow(w2,2) + 
   2*c1*j2*m2*q2*pow(w2,2) + e1*i2*pow(q2,2)*pow(w2,2) + 
   2*c2*j2*l2*r1*pow(w2,2) - 2*c2*i2*q2*r1*pow(w2,2) + 
   2*l2*m2*q2*r1*pow(w2,2) + 2*c1*j2*l2*r2*pow(w2,2) - 
   2*c1*i2*q2*r2*pow(w2,2) - 2*pow(l2,2)*r1*r2*pow(w2,2) + 
   pow(f2,2)*pow(k2,2)*o2*z1 - e2*g2*pow(k2,2)*o2*z1 - 
   2*pow(f2,2)*j2*k2*p2*z1 + 2*e2*g2*j2*k2*p2*z1 + 
   pow(f2,2)*i2*pow(p2,2)*z1 - e2*g2*i2*pow(p2,2)*z1 + 
   g2*pow(m2,2)*pow(p2,2)*z1 - 2*f2*m2*n2*pow(p2,2)*z1 + 
   e2*pow(n2,2)*pow(p2,2)*z1 - 2*g2*k2*m2*p2*r2*z1 + 
   2*f2*k2*n2*p2*r2*z1 + g2*pow(k2,2)*pow(r2,2)*z1 + 
   2*f2*k2*m2*p2*s2*z1 - 2*e2*k2*n2*p2*s2*z1 - 2*f2*pow(k2,2)*r2*s2*z1 + 
   e2*pow(k2,2)*pow(s2,2)*z1 + pow(f2,2)*pow(j2,2)*t2*z1 - 
   e2*g2*pow(j2,2)*t2*z1 - pow(f2,2)*i2*o2*t2*z1 + e2*g2*i2*o2*t2*z1 - 
   g2*pow(m2,2)*o2*t2*z1 + 2*f2*m2*n2*o2*t2*z1 - e2*pow(n2,2)*o2*t2*z1 + 
   2*g2*j2*m2*r2*t2*z1 - 2*f2*j2*n2*r2*t2*z1 - g2*i2*pow(r2,2)*t2*z1 + 
   pow(n2,2)*pow(r2,2)*t2*z1 - 2*f2*j2*m2*s2*t2*z1 + 
   2*e2*j2*n2*s2*t2*z1 + 2*f2*i2*r2*s2*t2*z1 - 2*m2*n2*r2*s2*t2*z1 - 
   e2*i2*pow(s2,2)*t2*z1 + pow(m2,2)*pow(s2,2)*t2*z1 + 
   2*g2*k2*m2*o2*v2*z1 - 2*f2*k2*n2*o2*v2*z1 - 2*g2*j2*m2*p2*v2*z1 + 
   2*f2*j2*n2*p2*v2*z1 - 2*g2*j2*k2*r2*v2*z1 + 2*g2*i2*p2*r2*v2*z1 - 
   2*pow(n2,2)*p2*r2*v2*z1 + 2*f2*j2*k2*s2*v2*z1 - 2*f2*i2*p2*s2*v2*z1 + 
   2*m2*n2*p2*s2*v2*z1 + 2*k2*n2*r2*s2*v2*z1 - 2*k2*m2*pow(s2,2)*v2*z1 + 
   g2*pow(j2,2)*pow(v2,2)*z1 - g2*i2*o2*pow(v2,2)*z1 + 
   pow(n2,2)*o2*pow(v2,2)*z1 - 2*j2*n2*s2*pow(v2,2)*z1 + 
   i2*pow(s2,2)*pow(v2,2)*z1 - 2*f2*k2*m2*o2*w2*z1 + 
   2*e2*k2*n2*o2*w2*z1 + 2*f2*j2*m2*p2*w2*z1 - 2*e2*j2*n2*p2*w2*z1 + 
   2*f2*j2*k2*r2*w2*z1 - 2*f2*i2*p2*r2*w2*z1 + 2*m2*n2*p2*r2*w2*z1 - 
   2*k2*n2*pow(r2,2)*w2*z1 - 2*e2*j2*k2*s2*w2*z1 + 2*e2*i2*p2*s2*w2*z1 - 
   2*pow(m2,2)*p2*s2*w2*z1 + 2*k2*m2*r2*s2*w2*z1 - 
   2*f2*pow(j2,2)*v2*w2*z1 + 2*f2*i2*o2*v2*w2*z1 - 2*m2*n2*o2*v2*w2*z1 + 
   2*j2*n2*r2*v2*w2*z1 + 2*j2*m2*s2*v2*w2*z1 - 2*i2*r2*s2*v2*w2*z1 + 
   e2*pow(j2,2)*pow(w2,2)*z1 - e2*i2*o2*pow(w2,2)*z1 + 
   pow(m2,2)*o2*pow(w2,2)*z1 - 2*j2*m2*r2*pow(w2,2)*z1 + 
   i2*pow(r2,2)*pow(w2,2)*z1 + 2*f1*f2*pow(k2,2)*o2*z2 - 
   e2*g1*pow(k2,2)*o2*z2 - e1*g2*pow(k2,2)*o2*z2 - 4*f1*f2*j2*k2*p2*z2 + 
   2*e2*g1*j2*k2*p2*z2 + 2*e1*g2*j2*k2*p2*z2 + 2*f1*f2*i2*pow(p2,2)*z2 - 
   e2*g1*i2*pow(p2,2)*z2 - e1*g2*i2*pow(p2,2)*z2 + 
   g1*pow(m2,2)*pow(p2,2)*z2 - 2*f2*m2*n1*pow(p2,2)*z2 - 
   2*f1*m2*n2*pow(p2,2)*z2 + 2*e2*n1*n2*pow(p2,2)*z2 + 
   e1*pow(n2,2)*pow(p2,2)*z2 - 2*g2*k2*m2*p2*r1*z2 + 
   2*f2*k2*n2*p2*r1*z2 - 2*g1*k2*m2*p2*r2*z2 + 2*f2*k2*n1*p2*r2*z2 + 
   2*f1*k2*n2*p2*r2*z2 + 2*g2*pow(k2,2)*r1*r2*z2 + 
   g1*pow(k2,2)*pow(r2,2)*z2 + 2*f2*k2*m2*p2*s1*z2 - 
   2*e2*k2*n2*p2*s1*z2 - 2*f2*pow(k2,2)*r2*s1*z2 + 2*f1*k2*m2*p2*s2*z2 - 
   2*e2*k2*n1*p2*s2*z2 - 2*e1*k2*n2*p2*s2*z2 - 2*f2*pow(k2,2)*r1*s2*z2 - 
   2*f1*pow(k2,2)*r2*s2*z2 + 2*e2*pow(k2,2)*s1*s2*z2 + 
   e1*pow(k2,2)*pow(s2,2)*z2 + 2*f1*f2*pow(j2,2)*t2*z2 - 
   e2*g1*pow(j2,2)*t2*z2 - e1*g2*pow(j2,2)*t2*z2 - 2*f1*f2*i2*o2*t2*z2 + 
   e2*g1*i2*o2*t2*z2 + e1*g2*i2*o2*t2*z2 - g1*pow(m2,2)*o2*t2*z2 + 
   2*f2*m2*n1*o2*t2*z2 + 2*f1*m2*n2*o2*t2*z2 - 2*e2*n1*n2*o2*t2*z2 - 
   e1*pow(n2,2)*o2*t2*z2 + 2*g2*j2*m2*r1*t2*z2 - 2*f2*j2*n2*r1*t2*z2 + 
   2*g1*j2*m2*r2*t2*z2 - 2*f2*j2*n1*r2*t2*z2 - 2*f1*j2*n2*r2*t2*z2 - 
   2*g2*i2*r1*r2*t2*z2 + 2*pow(n2,2)*r1*r2*t2*z2 - 
   g1*i2*pow(r2,2)*t2*z2 + 2*n1*n2*pow(r2,2)*t2*z2 - 
   2*f2*j2*m2*s1*t2*z2 + 2*e2*j2*n2*s1*t2*z2 + 2*f2*i2*r2*s1*t2*z2 - 
   2*m2*n2*r2*s1*t2*z2 - 2*f1*j2*m2*s2*t2*z2 + 2*e2*j2*n1*s2*t2*z2 + 
   2*e1*j2*n2*s2*t2*z2 + 2*f2*i2*r1*s2*t2*z2 - 2*m2*n2*r1*s2*t2*z2 + 
   2*f1*i2*r2*s2*t2*z2 - 2*m2*n1*r2*s2*t2*z2 - 2*e2*i2*s1*s2*t2*z2 + 
   2*pow(m2,2)*s1*s2*t2*z2 - e1*i2*pow(s2,2)*t2*z2 + 
   2*g2*k2*m2*o2*v1*z2 - 2*f2*k2*n2*o2*v1*z2 - 2*g2*j2*m2*p2*v1*z2 + 
   2*f2*j2*n2*p2*v1*z2 - 2*g2*j2*k2*r2*v1*z2 + 2*g2*i2*p2*r2*v1*z2 - 
   2*pow(n2,2)*p2*r2*v1*z2 + 2*f2*j2*k2*s2*v1*z2 - 2*f2*i2*p2*s2*v1*z2 + 
   2*m2*n2*p2*s2*v1*z2 + 2*k2*n2*r2*s2*v1*z2 - 2*k2*m2*pow(s2,2)*v1*z2 + 
   2*g1*k2*m2*o2*v2*z2 - 2*f2*k2*n1*o2*v2*z2 - 2*f1*k2*n2*o2*v2*z2 - 
   2*g1*j2*m2*p2*v2*z2 + 2*f2*j2*n1*p2*v2*z2 + 2*f1*j2*n2*p2*v2*z2 - 
   2*g2*j2*k2*r1*v2*z2 + 2*g2*i2*p2*r1*v2*z2 - 2*pow(n2,2)*p2*r1*v2*z2 - 
   2*g1*j2*k2*r2*v2*z2 + 2*g1*i2*p2*r2*v2*z2 - 4*n1*n2*p2*r2*v2*z2 + 
   2*f2*j2*k2*s1*v2*z2 - 2*f2*i2*p2*s1*v2*z2 + 2*m2*n2*p2*s1*v2*z2 + 
   2*k2*n2*r2*s1*v2*z2 + 2*f1*j2*k2*s2*v2*z2 - 2*f1*i2*p2*s2*v2*z2 + 
   2*m2*n1*p2*s2*v2*z2 + 2*k2*n2*r1*s2*v2*z2 + 2*k2*n1*r2*s2*v2*z2 - 
   4*k2*m2*s1*s2*v2*z2 + 2*g2*pow(j2,2)*v1*v2*z2 - 2*g2*i2*o2*v1*v2*z2 + 
   2*pow(n2,2)*o2*v1*v2*z2 - 4*j2*n2*s2*v1*v2*z2 + 
   2*i2*pow(s2,2)*v1*v2*z2 + g1*pow(j2,2)*pow(v2,2)*z2 - 
   g1*i2*o2*pow(v2,2)*z2 + 2*n1*n2*o2*pow(v2,2)*z2 - 
   2*j2*n2*s1*pow(v2,2)*z2 - 2*j2*n1*s2*pow(v2,2)*z2 + 
   2*i2*s1*s2*pow(v2,2)*z2 - 2*f2*k2*m2*o2*w1*z2 + 2*e2*k2*n2*o2*w1*z2 + 
   2*f2*j2*m2*p2*w1*z2 - 2*e2*j2*n2*p2*w1*z2 + 2*f2*j2*k2*r2*w1*z2 - 
   2*f2*i2*p2*r2*w1*z2 + 2*m2*n2*p2*r2*w1*z2 - 2*k2*n2*pow(r2,2)*w1*z2 - 
   2*e2*j2*k2*s2*w1*z2 + 2*e2*i2*p2*s2*w1*z2 - 2*pow(m2,2)*p2*s2*w1*z2 + 
   2*k2*m2*r2*s2*w1*z2 - 2*f2*pow(j2,2)*v2*w1*z2 + 2*f2*i2*o2*v2*w1*z2 - 
   2*m2*n2*o2*v2*w1*z2 + 2*j2*n2*r2*v2*w1*z2 + 2*j2*m2*s2*v2*w1*z2 - 
   2*i2*r2*s2*v2*w1*z2 - 2*f1*k2*m2*o2*w2*z2 + 2*e2*k2*n1*o2*w2*z2 + 
   2*e1*k2*n2*o2*w2*z2 + 2*f1*j2*m2*p2*w2*z2 - 2*e2*j2*n1*p2*w2*z2 - 
   2*e1*j2*n2*p2*w2*z2 + 2*f2*j2*k2*r1*w2*z2 - 2*f2*i2*p2*r1*w2*z2 + 
   2*m2*n2*p2*r1*w2*z2 + 2*f1*j2*k2*r2*w2*z2 - 2*f1*i2*p2*r2*w2*z2 + 
   2*m2*n1*p2*r2*w2*z2 - 4*k2*n2*r1*r2*w2*z2 - 2*k2*n1*pow(r2,2)*w2*z2 - 
   2*e2*j2*k2*s1*w2*z2 + 2*e2*i2*p2*s1*w2*z2 - 2*pow(m2,2)*p2*s1*w2*z2 + 
   2*k2*m2*r2*s1*w2*z2 - 2*e1*j2*k2*s2*w2*z2 + 2*e1*i2*p2*s2*w2*z2 + 
   2*k2*m2*r1*s2*w2*z2 - 2*f2*pow(j2,2)*v1*w2*z2 + 2*f2*i2*o2*v1*w2*z2 - 
   2*m2*n2*o2*v1*w2*z2 + 2*j2*n2*r2*v1*w2*z2 + 2*j2*m2*s2*v1*w2*z2 - 
   2*i2*r2*s2*v1*w2*z2 - 2*f1*pow(j2,2)*v2*w2*z2 + 2*f1*i2*o2*v2*w2*z2 - 
   2*m2*n1*o2*v2*w2*z2 + 2*j2*n2*r1*v2*w2*z2 + 2*j2*n1*r2*v2*w2*z2 + 
   2*j2*m2*s1*v2*w2*z2 - 2*i2*r2*s1*v2*w2*z2 - 2*i2*r1*s2*v2*w2*z2 + 
   2*e2*pow(j2,2)*w1*w2*z2 - 2*e2*i2*o2*w1*w2*z2 + 
   2*pow(m2,2)*o2*w1*w2*z2 - 4*j2*m2*r2*w1*w2*z2 + 
   2*i2*pow(r2,2)*w1*w2*z2 + e1*pow(j2,2)*pow(w2,2)*z2 - 
   e1*i2*o2*pow(w2,2)*z2 - 2*j2*m2*r1*pow(w2,2)*z2 + 
   2*i2*r1*r2*pow(w2,2)*z2
;
    
	k[3]  = 2*d0*d2*e0*pow(k2,2)*o2 + pow(d0,2)*e2*pow(k2,2)*o2 - 
   2*c2*d0*f0*pow(k2,2)*o2 - 2*c0*d2*f0*pow(k2,2)*o2 - 
   2*c0*d0*f2*pow(k2,2)*o2 + 2*c0*c2*g0*pow(k2,2)*o2 + 
   pow(c0,2)*g2*pow(k2,2)*o2 - 4*d0*d2*e0*j2*k2*p2 - 
   2*pow(d0,2)*e2*j2*k2*p2 + 4*c2*d0*f0*j2*k2*p2 + 4*c0*d2*f0*j2*k2*p2 + 
   4*c0*d0*f2*j2*k2*p2 - 4*c0*c2*g0*j2*k2*p2 - 2*pow(c0,2)*g2*j2*k2*p2 + 
   2*d0*d2*e0*i2*pow(p2,2) + pow(d0,2)*e2*i2*pow(p2,2) - 
   2*c2*d0*f0*i2*pow(p2,2) - 2*c0*d2*f0*i2*pow(p2,2) - 
   2*c0*d0*f2*i2*pow(p2,2) + 2*c0*c2*g0*i2*pow(p2,2) + 
   pow(c0,2)*g2*i2*pow(p2,2) - pow(f0,2)*pow(l2,2)*pow(p2,2) + 
   e0*g0*pow(l2,2)*pow(p2,2) + 2*d0*f0*l2*m2*pow(p2,2) - 
   2*c0*g0*l2*m2*pow(p2,2) - pow(d0,2)*pow(m2,2)*pow(p2,2) - 
   2*d2*e0*l2*n0*pow(p2,2) - 2*d0*e2*l2*n0*pow(p2,2) + 
   2*c2*f0*l2*n0*pow(p2,2) + 2*c0*f2*l2*n0*pow(p2,2) + 
   2*c2*d0*m2*n0*pow(p2,2) + 2*c0*d2*m2*n0*pow(p2,2) - 
   pow(c2,2)*pow(n0,2)*pow(p2,2) - 2*d0*e0*l2*n2*pow(p2,2) + 
   2*c0*f0*l2*n2*pow(p2,2) + 2*c0*d0*m2*n2*pow(p2,2) - 
   4*c0*c2*n0*n2*pow(p2,2) - pow(c0,2)*pow(n2,2)*pow(p2,2) + 
   2*pow(f0,2)*k2*l2*p2*q2 - 2*e0*g0*k2*l2*p2*q2 - 2*d0*f0*k2*m2*p2*q2 + 
   2*c0*g0*k2*m2*p2*q2 + 2*d2*e0*k2*n0*p2*q2 + 2*d0*e2*k2*n0*p2*q2 - 
   2*c2*f0*k2*n0*p2*q2 - 2*c0*f2*k2*n0*p2*q2 + 2*d0*e0*k2*n2*p2*q2 - 
   2*c0*f0*k2*n2*p2*q2 - pow(f0,2)*pow(k2,2)*pow(q2,2) + 
   e0*g0*pow(k2,2)*pow(q2,2) - 2*d2*f0*k2*l2*p2*r0 - 
   2*d0*f2*k2*l2*p2*r0 + 2*c2*g0*k2*l2*p2*r0 + 2*c0*g2*k2*l2*p2*r0 + 
   4*d0*d2*k2*m2*p2*r0 - 2*c2*d2*k2*n0*p2*r0 - 2*c2*d0*k2*n2*p2*r0 - 
   2*c0*d2*k2*n2*p2*r0 + 2*d2*f0*pow(k2,2)*q2*r0 + 
   2*d0*f2*pow(k2,2)*q2*r0 - 2*c2*g0*pow(k2,2)*q2*r0 - 
   2*c0*g2*pow(k2,2)*q2*r0 - pow(d2,2)*pow(k2,2)*pow(r0,2) - 
   2*d0*f0*k2*l2*p2*r2 + 2*c0*g0*k2*l2*p2*r2 + 2*pow(d0,2)*k2*m2*p2*r2 - 
   2*c2*d0*k2*n0*p2*r2 - 2*c0*d2*k2*n0*p2*r2 - 2*c0*d0*k2*n2*p2*r2 + 
   2*d0*f0*pow(k2,2)*q2*r2 - 2*c0*g0*pow(k2,2)*q2*r2 - 
   4*d0*d2*pow(k2,2)*r0*r2 - pow(d0,2)*pow(k2,2)*pow(r2,2) + 
   2*d2*e0*k2*l2*p2*s0 + 2*d0*e2*k2*l2*p2*s0 - 2*c2*f0*k2*l2*p2*s0 - 
   2*c0*f2*k2*l2*p2*s0 - 2*c2*d0*k2*m2*p2*s0 - 2*c0*d2*k2*m2*p2*s0 + 
   2*pow(c2,2)*k2*n0*p2*s0 + 4*c0*c2*k2*n2*p2*s0 - 
   2*d2*e0*pow(k2,2)*q2*s0 - 2*d0*e2*pow(k2,2)*q2*s0 + 
   2*c2*f0*pow(k2,2)*q2*s0 + 2*c0*f2*pow(k2,2)*q2*s0 + 
   2*c2*d2*pow(k2,2)*r0*s0 + 2*c2*d0*pow(k2,2)*r2*s0 + 
   2*c0*d2*pow(k2,2)*r2*s0 - pow(c2,2)*pow(k2,2)*pow(s0,2) + 
   2*d0*e0*k2*l2*p2*s2 - 2*c0*f0*k2*l2*p2*s2 - 2*c0*d0*k2*m2*p2*s2 + 
   4*c0*c2*k2*n0*p2*s2 + 2*pow(c0,2)*k2*n2*p2*s2 - 
   2*d0*e0*pow(k2,2)*q2*s2 + 2*c0*f0*pow(k2,2)*q2*s2 + 
   2*c2*d0*pow(k2,2)*r0*s2 + 2*c0*d2*pow(k2,2)*r0*s2 + 
   2*c0*d0*pow(k2,2)*r2*s2 - 4*c0*c2*pow(k2,2)*s0*s2 - 
   pow(c0,2)*pow(k2,2)*pow(s2,2) + 2*d0*d2*e0*pow(j2,2)*t2 + 
   pow(d0,2)*e2*pow(j2,2)*t2 - 2*c2*d0*f0*pow(j2,2)*t2 - 
   2*c0*d2*f0*pow(j2,2)*t2 - 2*c0*d0*f2*pow(j2,2)*t2 + 
   2*c0*c2*g0*pow(j2,2)*t2 + pow(c0,2)*g2*pow(j2,2)*t2 - 
   2*d0*d2*e0*i2*o2*t2 - pow(d0,2)*e2*i2*o2*t2 + 2*c2*d0*f0*i2*o2*t2 + 
   2*c0*d2*f0*i2*o2*t2 + 2*c0*d0*f2*i2*o2*t2 - 2*c0*c2*g0*i2*o2*t2 - 
   pow(c0,2)*g2*i2*o2*t2 + pow(f0,2)*pow(l2,2)*o2*t2 - 
   e0*g0*pow(l2,2)*o2*t2 - 2*d0*f0*l2*m2*o2*t2 + 2*c0*g0*l2*m2*o2*t2 + 
   pow(d0,2)*pow(m2,2)*o2*t2 + 2*d2*e0*l2*n0*o2*t2 + 
   2*d0*e2*l2*n0*o2*t2 - 2*c2*f0*l2*n0*o2*t2 - 2*c0*f2*l2*n0*o2*t2 - 
   2*c2*d0*m2*n0*o2*t2 - 2*c0*d2*m2*n0*o2*t2 + 
   pow(c2,2)*pow(n0,2)*o2*t2 + 2*d0*e0*l2*n2*o2*t2 - 
   2*c0*f0*l2*n2*o2*t2 - 2*c0*d0*m2*n2*o2*t2 + 4*c0*c2*n0*n2*o2*t2 + 
   pow(c0,2)*pow(n2,2)*o2*t2 - 2*pow(f0,2)*j2*l2*q2*t2 + 
   2*e0*g0*j2*l2*q2*t2 + 2*d0*f0*j2*m2*q2*t2 - 2*c0*g0*j2*m2*q2*t2 - 
   2*d2*e0*j2*n0*q2*t2 - 2*d0*e2*j2*n0*q2*t2 + 2*c2*f0*j2*n0*q2*t2 + 
   2*c0*f2*j2*n0*q2*t2 - 2*d0*e0*j2*n2*q2*t2 + 2*c0*f0*j2*n2*q2*t2 + 
   pow(f0,2)*i2*pow(q2,2)*t2 - e0*g0*i2*pow(q2,2)*t2 - 
   2*f0*m2*n0*pow(q2,2)*t2 + e2*pow(n0,2)*pow(q2,2)*t2 + 
   2*e0*n0*n2*pow(q2,2)*t2 + 2*d2*f0*j2*l2*r0*t2 + 2*d0*f2*j2*l2*r0*t2 - 
   2*c2*g0*j2*l2*r0*t2 - 2*c0*g2*j2*l2*r0*t2 - 4*d0*d2*j2*m2*r0*t2 + 
   2*c2*d2*j2*n0*r0*t2 + 2*c2*d0*j2*n2*r0*t2 + 2*c0*d2*j2*n2*r0*t2 - 
   2*d2*f0*i2*q2*r0*t2 - 2*d0*f2*i2*q2*r0*t2 + 2*c2*g0*i2*q2*r0*t2 + 
   2*c0*g2*i2*q2*r0*t2 - 2*g0*l2*m2*q2*r0*t2 + 2*f2*l2*n0*q2*r0*t2 + 
   2*d2*m2*n0*q2*r0*t2 + 2*f0*l2*n2*q2*r0*t2 + 2*d0*m2*n2*q2*r0*t2 - 
   4*c2*n0*n2*q2*r0*t2 - 2*c0*pow(n2,2)*q2*r0*t2 + 
   pow(d2,2)*i2*pow(r0,2)*t2 + g2*pow(l2,2)*pow(r0,2)*t2 - 
   2*d2*l2*n2*pow(r0,2)*t2 + 2*d0*f0*j2*l2*r2*t2 - 2*c0*g0*j2*l2*r2*t2 - 
   2*pow(d0,2)*j2*m2*r2*t2 + 2*c2*d0*j2*n0*r2*t2 + 2*c0*d2*j2*n0*r2*t2 + 
   2*c0*d0*j2*n2*r2*t2 - 2*d0*f0*i2*q2*r2*t2 + 2*c0*g0*i2*q2*r2*t2 + 
   2*f0*l2*n0*q2*r2*t2 + 2*d0*m2*n0*q2*r2*t2 - 2*c2*pow(n0,2)*q2*r2*t2 - 
   4*c0*n0*n2*q2*r2*t2 + 4*d0*d2*i2*r0*r2*t2 + 2*g0*pow(l2,2)*r0*r2*t2 - 
   4*d2*l2*n0*r0*r2*t2 - 4*d0*l2*n2*r0*r2*t2 + 
   pow(d0,2)*i2*pow(r2,2)*t2 - 2*d0*l2*n0*pow(r2,2)*t2 - 
   2*d2*e0*j2*l2*s0*t2 - 2*d0*e2*j2*l2*s0*t2 + 2*c2*f0*j2*l2*s0*t2 + 
   2*c0*f2*j2*l2*s0*t2 + 2*c2*d0*j2*m2*s0*t2 + 2*c0*d2*j2*m2*s0*t2 - 
   2*pow(c2,2)*j2*n0*s0*t2 - 4*c0*c2*j2*n2*s0*t2 + 2*d2*e0*i2*q2*s0*t2 + 
   2*d0*e2*i2*q2*s0*t2 - 2*c2*f0*i2*q2*s0*t2 - 2*c0*f2*i2*q2*s0*t2 + 
   2*f0*l2*m2*q2*s0*t2 - 2*d0*pow(m2,2)*q2*s0*t2 - 2*e2*l2*n0*q2*s0*t2 + 
   2*c2*m2*n0*q2*s0*t2 - 2*e0*l2*n2*q2*s0*t2 + 2*c0*m2*n2*q2*s0*t2 - 
   2*c2*d2*i2*r0*s0*t2 - 2*f2*pow(l2,2)*r0*s0*t2 + 2*d2*l2*m2*r0*s0*t2 + 
   2*c2*l2*n2*r0*s0*t2 - 2*c2*d0*i2*r2*s0*t2 - 2*c0*d2*i2*r2*s0*t2 - 
   2*f0*pow(l2,2)*r2*s0*t2 + 2*d0*l2*m2*r2*s0*t2 + 2*c2*l2*n0*r2*s0*t2 + 
   2*c0*l2*n2*r2*s0*t2 + pow(c2,2)*i2*pow(s0,2)*t2 + 
   e2*pow(l2,2)*pow(s0,2)*t2 - 2*c2*l2*m2*pow(s0,2)*t2 - 
   2*d0*e0*j2*l2*s2*t2 + 2*c0*f0*j2*l2*s2*t2 + 2*c0*d0*j2*m2*s2*t2 - 
   4*c0*c2*j2*n0*s2*t2 - 2*pow(c0,2)*j2*n2*s2*t2 + 2*d0*e0*i2*q2*s2*t2 - 
   2*c0*f0*i2*q2*s2*t2 - 2*e0*l2*n0*q2*s2*t2 + 2*c0*m2*n0*q2*s2*t2 - 
   2*c2*d0*i2*r0*s2*t2 - 2*c0*d2*i2*r0*s2*t2 - 2*f0*pow(l2,2)*r0*s2*t2 + 
   2*d0*l2*m2*r0*s2*t2 + 2*c2*l2*n0*r0*s2*t2 + 2*c0*l2*n2*r0*s2*t2 - 
   2*c0*d0*i2*r2*s2*t2 + 2*c0*l2*n0*r2*s2*t2 + 4*c0*c2*i2*s0*s2*t2 + 
   2*e0*pow(l2,2)*s0*s2*t2 - 4*c0*l2*m2*s0*s2*t2 + 
   pow(c0,2)*i2*pow(s2,2)*t2 - 4*f0*f2*k2*l2*o2*u0 + 
   2*e2*g0*k2*l2*o2*u0 + 2*e0*g2*k2*l2*o2*u0 + 2*d2*f0*k2*m2*o2*u0 + 
   2*d0*f2*k2*m2*o2*u0 - 2*c2*g0*k2*m2*o2*u0 - 2*c0*g2*k2*m2*o2*u0 - 
   2*d2*e2*k2*n0*o2*u0 + 2*c2*f2*k2*n0*o2*u0 - 2*d2*e0*k2*n2*o2*u0 - 
   2*d0*e2*k2*n2*o2*u0 + 2*c2*f0*k2*n2*o2*u0 + 2*c0*f2*k2*n2*o2*u0 + 
   4*f0*f2*j2*l2*p2*u0 - 2*e2*g0*j2*l2*p2*u0 - 2*e0*g2*j2*l2*p2*u0 - 
   2*d2*f0*j2*m2*p2*u0 - 2*d0*f2*j2*m2*p2*u0 + 2*c2*g0*j2*m2*p2*u0 + 
   2*c0*g2*j2*m2*p2*u0 + 2*d2*e2*j2*n0*p2*u0 - 2*c2*f2*j2*n0*p2*u0 + 
   2*d2*e0*j2*n2*p2*u0 + 2*d0*e2*j2*n2*p2*u0 - 2*c2*f0*j2*n2*p2*u0 - 
   2*c0*f2*j2*n2*p2*u0 + 4*f0*f2*j2*k2*q2*u0 - 2*e2*g0*j2*k2*q2*u0 - 
   2*e0*g2*j2*k2*q2*u0 - 4*f0*f2*i2*p2*q2*u0 + 2*e2*g0*i2*p2*q2*u0 + 
   2*e0*g2*i2*p2*q2*u0 - 2*g0*pow(m2,2)*p2*q2*u0 + 4*f2*m2*n0*p2*q2*u0 + 
   4*f0*m2*n2*p2*q2*u0 - 4*e2*n0*n2*p2*q2*u0 - 2*e0*pow(n2,2)*p2*q2*u0 - 
   2*d2*f2*j2*k2*r0*u0 + 2*c2*g2*j2*k2*r0*u0 + 2*d2*f2*i2*p2*r0*u0 - 
   2*c2*g2*i2*p2*r0*u0 + 2*g2*l2*m2*p2*r0*u0 - 2*f2*l2*n2*p2*r0*u0 - 
   2*d2*m2*n2*p2*r0*u0 + 2*c2*pow(n2,2)*p2*r0*u0 + 2*g2*k2*m2*q2*r0*u0 - 
   2*f2*k2*n2*q2*r0*u0 - 2*d2*f0*j2*k2*r2*u0 - 2*d0*f2*j2*k2*r2*u0 + 
   2*c2*g0*j2*k2*r2*u0 + 2*c0*g2*j2*k2*r2*u0 + 2*d2*f0*i2*p2*r2*u0 + 
   2*d0*f2*i2*p2*r2*u0 - 2*c2*g0*i2*p2*r2*u0 - 2*c0*g2*i2*p2*r2*u0 + 
   2*g0*l2*m2*p2*r2*u0 - 2*f2*l2*n0*p2*r2*u0 - 2*d2*m2*n0*p2*r2*u0 - 
   2*f0*l2*n2*p2*r2*u0 - 2*d0*m2*n2*p2*r2*u0 + 4*c2*n0*n2*p2*r2*u0 + 
   2*c0*pow(n2,2)*p2*r2*u0 + 2*g0*k2*m2*q2*r2*u0 - 2*f2*k2*n0*q2*r2*u0 - 
   2*f0*k2*n2*q2*r2*u0 - 4*g2*k2*l2*r0*r2*u0 + 4*d2*k2*n2*r0*r2*u0 - 
   2*g0*k2*l2*pow(r2,2)*u0 + 2*d2*k2*n0*pow(r2,2)*u0 + 
   2*d0*k2*n2*pow(r2,2)*u0 + 2*d2*e2*j2*k2*s0*u0 - 2*c2*f2*j2*k2*s0*u0 - 
   2*d2*e2*i2*p2*s0*u0 + 2*c2*f2*i2*p2*s0*u0 - 2*f2*l2*m2*p2*s0*u0 + 
   2*d2*pow(m2,2)*p2*s0*u0 + 2*e2*l2*n2*p2*s0*u0 - 2*c2*m2*n2*p2*s0*u0 - 
   2*f2*k2*m2*q2*s0*u0 + 2*e2*k2*n2*q2*s0*u0 + 4*f2*k2*l2*r2*s0*u0 - 
   2*d2*k2*m2*r2*s0*u0 - 2*c2*k2*n2*r2*s0*u0 + 2*d2*e0*j2*k2*s2*u0 + 
   2*d0*e2*j2*k2*s2*u0 - 2*c2*f0*j2*k2*s2*u0 - 2*c0*f2*j2*k2*s2*u0 - 
   2*d2*e0*i2*p2*s2*u0 - 2*d0*e2*i2*p2*s2*u0 + 2*c2*f0*i2*p2*s2*u0 + 
   2*c0*f2*i2*p2*s2*u0 - 2*f0*l2*m2*p2*s2*u0 + 2*d0*pow(m2,2)*p2*s2*u0 + 
   2*e2*l2*n0*p2*s2*u0 - 2*c2*m2*n0*p2*s2*u0 + 2*e0*l2*n2*p2*s2*u0 - 
   2*c0*m2*n2*p2*s2*u0 - 2*f0*k2*m2*q2*s2*u0 + 2*e2*k2*n0*q2*s2*u0 + 
   2*e0*k2*n2*q2*s2*u0 + 4*f2*k2*l2*r0*s2*u0 - 2*d2*k2*m2*r0*s2*u0 - 
   2*c2*k2*n2*r0*s2*u0 + 4*f0*k2*l2*r2*s2*u0 - 2*d0*k2*m2*r2*s2*u0 - 
   2*c2*k2*n0*r2*s2*u0 - 2*c0*k2*n2*r2*s2*u0 - 4*e2*k2*l2*s0*s2*u0 + 
   4*c2*k2*m2*s0*s2*u0 - 2*e0*k2*l2*pow(s2,2)*u0 + 
   2*c0*k2*m2*pow(s2,2)*u0 - pow(f2,2)*pow(j2,2)*pow(u0,2) + 
   e2*g2*pow(j2,2)*pow(u0,2) + pow(f2,2)*i2*o2*pow(u0,2) - 
   e2*g2*i2*o2*pow(u0,2) + g2*pow(m2,2)*o2*pow(u0,2) - 
   2*f2*m2*n2*o2*pow(u0,2) + e2*pow(n2,2)*o2*pow(u0,2) - 
   2*g2*j2*m2*r2*pow(u0,2) + 2*f2*j2*n2*r2*pow(u0,2) + 
   g2*i2*pow(r2,2)*pow(u0,2) - pow(n2,2)*pow(r2,2)*pow(u0,2) + 
   2*f2*j2*m2*s2*pow(u0,2) - 2*e2*j2*n2*s2*pow(u0,2) - 
   2*f2*i2*r2*s2*pow(u0,2) + 2*m2*n2*r2*s2*pow(u0,2) + 
   e2*i2*pow(s2,2)*pow(u0,2) - pow(m2,2)*pow(s2,2)*pow(u0,2) - 
   2*pow(f0,2)*k2*l2*o2*u2 + 2*e0*g0*k2*l2*o2*u2 + 2*d0*f0*k2*m2*o2*u2 - 
   2*c0*g0*k2*m2*o2*u2 - 2*d2*e0*k2*n0*o2*u2 - 2*d0*e2*k2*n0*o2*u2 + 
   2*c2*f0*k2*n0*o2*u2 + 2*c0*f2*k2*n0*o2*u2 - 2*d0*e0*k2*n2*o2*u2 + 
   2*c0*f0*k2*n2*o2*u2 + 2*pow(f0,2)*j2*l2*p2*u2 - 2*e0*g0*j2*l2*p2*u2 - 
   2*d0*f0*j2*m2*p2*u2 + 2*c0*g0*j2*m2*p2*u2 + 2*d2*e0*j2*n0*p2*u2 + 
   2*d0*e2*j2*n0*p2*u2 - 2*c2*f0*j2*n0*p2*u2 - 2*c0*f2*j2*n0*p2*u2 + 
   2*d0*e0*j2*n2*p2*u2 - 2*c0*f0*j2*n2*p2*u2 + 2*pow(f0,2)*j2*k2*q2*u2 - 
   2*e0*g0*j2*k2*q2*u2 - 2*pow(f0,2)*i2*p2*q2*u2 + 2*e0*g0*i2*p2*q2*u2 + 
   4*f0*m2*n0*p2*q2*u2 - 2*e2*pow(n0,2)*p2*q2*u2 - 4*e0*n0*n2*p2*q2*u2 - 
   2*d2*f0*j2*k2*r0*u2 - 2*d0*f2*j2*k2*r0*u2 + 2*c2*g0*j2*k2*r0*u2 + 
   2*c0*g2*j2*k2*r0*u2 + 2*d2*f0*i2*p2*r0*u2 + 2*d0*f2*i2*p2*r0*u2 - 
   2*c2*g0*i2*p2*r0*u2 - 2*c0*g2*i2*p2*r0*u2 + 2*g0*l2*m2*p2*r0*u2 - 
   2*f2*l2*n0*p2*r0*u2 - 2*d2*m2*n0*p2*r0*u2 - 2*f0*l2*n2*p2*r0*u2 - 
   2*d0*m2*n2*p2*r0*u2 + 4*c2*n0*n2*p2*r0*u2 + 2*c0*pow(n2,2)*p2*r0*u2 + 
   2*g0*k2*m2*q2*r0*u2 - 2*f2*k2*n0*q2*r0*u2 - 2*f0*k2*n2*q2*r0*u2 - 
   2*g2*k2*l2*pow(r0,2)*u2 + 2*d2*k2*n2*pow(r0,2)*u2 - 
   2*d0*f0*j2*k2*r2*u2 + 2*c0*g0*j2*k2*r2*u2 + 2*d0*f0*i2*p2*r2*u2 - 
   2*c0*g0*i2*p2*r2*u2 - 2*f0*l2*n0*p2*r2*u2 - 2*d0*m2*n0*p2*r2*u2 + 
   2*c2*pow(n0,2)*p2*r2*u2 + 4*c0*n0*n2*p2*r2*u2 - 2*f0*k2*n0*q2*r2*u2 - 
   4*g0*k2*l2*r0*r2*u2 + 4*d2*k2*n0*r0*r2*u2 + 4*d0*k2*n2*r0*r2*u2 + 
   2*d0*k2*n0*pow(r2,2)*u2 + 2*d2*e0*j2*k2*s0*u2 + 2*d0*e2*j2*k2*s0*u2 - 
   2*c2*f0*j2*k2*s0*u2 - 2*c0*f2*j2*k2*s0*u2 - 2*d2*e0*i2*p2*s0*u2 - 
   2*d0*e2*i2*p2*s0*u2 + 2*c2*f0*i2*p2*s0*u2 + 2*c0*f2*i2*p2*s0*u2 - 
   2*f0*l2*m2*p2*s0*u2 + 2*d0*pow(m2,2)*p2*s0*u2 + 2*e2*l2*n0*p2*s0*u2 - 
   2*c2*m2*n0*p2*s0*u2 + 2*e0*l2*n2*p2*s0*u2 - 2*c0*m2*n2*p2*s0*u2 - 
   2*f0*k2*m2*q2*s0*u2 + 2*e2*k2*n0*q2*s0*u2 + 2*e0*k2*n2*q2*s0*u2 + 
   4*f2*k2*l2*r0*s0*u2 - 2*d2*k2*m2*r0*s0*u2 - 2*c2*k2*n2*r0*s0*u2 + 
   4*f0*k2*l2*r2*s0*u2 - 2*d0*k2*m2*r2*s0*u2 - 2*c2*k2*n0*r2*s0*u2 - 
   2*c0*k2*n2*r2*s0*u2 - 2*e2*k2*l2*pow(s0,2)*u2 + 
   2*c2*k2*m2*pow(s0,2)*u2 + 2*d0*e0*j2*k2*s2*u2 - 2*c0*f0*j2*k2*s2*u2 - 
   2*d0*e0*i2*p2*s2*u2 + 2*c0*f0*i2*p2*s2*u2 + 2*e0*l2*n0*p2*s2*u2 - 
   2*c0*m2*n0*p2*s2*u2 + 2*e0*k2*n0*q2*s2*u2 + 4*f0*k2*l2*r0*s2*u2 - 
   2*d0*k2*m2*r0*s2*u2 - 2*c2*k2*n0*r0*s2*u2 - 2*c0*k2*n2*r0*s2*u2 - 
   2*c0*k2*n0*r2*s2*u2 - 4*e0*k2*l2*s0*s2*u2 + 4*c0*k2*m2*s0*s2*u2 - 
   4*f0*f2*pow(j2,2)*u0*u2 + 2*e2*g0*pow(j2,2)*u0*u2 + 
   2*e0*g2*pow(j2,2)*u0*u2 + 4*f0*f2*i2*o2*u0*u2 - 2*e2*g0*i2*o2*u0*u2 - 
   2*e0*g2*i2*o2*u0*u2 + 2*g0*pow(m2,2)*o2*u0*u2 - 4*f2*m2*n0*o2*u0*u2 - 
   4*f0*m2*n2*o2*u0*u2 + 4*e2*n0*n2*o2*u0*u2 + 2*e0*pow(n2,2)*o2*u0*u2 - 
   4*g2*j2*m2*r0*u0*u2 + 4*f2*j2*n2*r0*u0*u2 - 4*g0*j2*m2*r2*u0*u2 + 
   4*f2*j2*n0*r2*u0*u2 + 4*f0*j2*n2*r2*u0*u2 + 4*g2*i2*r0*r2*u0*u2 - 
   4*pow(n2,2)*r0*r2*u0*u2 + 2*g0*i2*pow(r2,2)*u0*u2 - 
   4*n0*n2*pow(r2,2)*u0*u2 + 4*f2*j2*m2*s0*u0*u2 - 4*e2*j2*n2*s0*u0*u2 - 
   4*f2*i2*r2*s0*u0*u2 + 4*m2*n2*r2*s0*u0*u2 + 4*f0*j2*m2*s2*u0*u2 - 
   4*e2*j2*n0*s2*u0*u2 - 4*e0*j2*n2*s2*u0*u2 - 4*f2*i2*r0*s2*u0*u2 + 
   4*m2*n2*r0*s2*u0*u2 - 4*f0*i2*r2*s2*u0*u2 + 4*m2*n0*r2*s2*u0*u2 + 
   4*e2*i2*s0*s2*u0*u2 - 4*pow(m2,2)*s0*s2*u0*u2 + 
   2*e0*i2*pow(s2,2)*u0*u2 - pow(f0,2)*pow(j2,2)*pow(u2,2) + 
   e0*g0*pow(j2,2)*pow(u2,2) + pow(f0,2)*i2*o2*pow(u2,2) - 
   e0*g0*i2*o2*pow(u2,2) - 2*f0*m2*n0*o2*pow(u2,2) + 
   e2*pow(n0,2)*o2*pow(u2,2) + 2*e0*n0*n2*o2*pow(u2,2) - 
   2*g0*j2*m2*r0*pow(u2,2) + 2*f2*j2*n0*r0*pow(u2,2) + 
   2*f0*j2*n2*r0*pow(u2,2) + g2*i2*pow(r0,2)*pow(u2,2) - 
   pow(n2,2)*pow(r0,2)*pow(u2,2) + 2*f0*j2*n0*r2*pow(u2,2) + 
   2*g0*i2*r0*r2*pow(u2,2) - 4*n0*n2*r0*r2*pow(u2,2) - 
   pow(n0,2)*pow(r2,2)*pow(u2,2) + 2*f0*j2*m2*s0*pow(u2,2) - 
   2*e2*j2*n0*s0*pow(u2,2) - 2*e0*j2*n2*s0*pow(u2,2) - 
   2*f2*i2*r0*s0*pow(u2,2) + 2*m2*n2*r0*s0*pow(u2,2) - 
   2*f0*i2*r2*s0*pow(u2,2) + 2*m2*n0*r2*s0*pow(u2,2) + 
   e2*i2*pow(s0,2)*pow(u2,2) - pow(m2,2)*pow(s0,2)*pow(u2,2) - 
   2*e0*j2*n0*s2*pow(u2,2) - 2*f0*i2*r0*s2*pow(u2,2) + 
   2*m2*n0*r0*s2*pow(u2,2) + 2*e0*i2*s0*s2*pow(u2,2) + 
   2*d2*f0*k2*l2*o2*v0 + 2*d0*f2*k2*l2*o2*v0 - 2*c2*g0*k2*l2*o2*v0 - 
   2*c0*g2*k2*l2*o2*v0 - 4*d0*d2*k2*m2*o2*v0 + 2*c2*d2*k2*n0*o2*v0 + 
   2*c2*d0*k2*n2*o2*v0 + 2*c0*d2*k2*n2*o2*v0 - 2*d2*f0*j2*l2*p2*v0 - 
   2*d0*f2*j2*l2*p2*v0 + 2*c2*g0*j2*l2*p2*v0 + 2*c0*g2*j2*l2*p2*v0 + 
   4*d0*d2*j2*m2*p2*v0 - 2*c2*d2*j2*n0*p2*v0 - 2*c2*d0*j2*n2*p2*v0 - 
   2*c0*d2*j2*n2*p2*v0 - 2*d2*f0*j2*k2*q2*v0 - 2*d0*f2*j2*k2*q2*v0 + 
   2*c2*g0*j2*k2*q2*v0 + 2*c0*g2*j2*k2*q2*v0 + 2*d2*f0*i2*p2*q2*v0 + 
   2*d0*f2*i2*p2*q2*v0 - 2*c2*g0*i2*p2*q2*v0 - 2*c0*g2*i2*p2*q2*v0 + 
   2*g0*l2*m2*p2*q2*v0 - 2*f2*l2*n0*p2*q2*v0 - 2*d2*m2*n0*p2*q2*v0 - 
   2*f0*l2*n2*p2*q2*v0 - 2*d0*m2*n2*p2*q2*v0 + 4*c2*n0*n2*p2*q2*v0 + 
   2*c0*pow(n2,2)*p2*q2*v0 - 2*g0*k2*m2*pow(q2,2)*v0 + 
   2*f2*k2*n0*pow(q2,2)*v0 + 2*f0*k2*n2*pow(q2,2)*v0 + 
   2*pow(d2,2)*j2*k2*r0*v0 - 2*pow(d2,2)*i2*p2*r0*v0 - 
   2*g2*pow(l2,2)*p2*r0*v0 + 4*d2*l2*n2*p2*r0*v0 + 2*g2*k2*l2*q2*r0*v0 - 
   2*d2*k2*n2*q2*r0*v0 + 4*d0*d2*j2*k2*r2*v0 - 4*d0*d2*i2*p2*r2*v0 - 
   2*g0*pow(l2,2)*p2*r2*v0 + 4*d2*l2*n0*p2*r2*v0 + 4*d0*l2*n2*p2*r2*v0 + 
   2*g0*k2*l2*q2*r2*v0 - 2*d2*k2*n0*q2*r2*v0 - 2*d0*k2*n2*q2*r2*v0 - 
   2*c2*d2*j2*k2*s0*v0 + 2*c2*d2*i2*p2*s0*v0 + 2*f2*pow(l2,2)*p2*s0*v0 - 
   2*d2*l2*m2*p2*s0*v0 - 2*c2*l2*n2*p2*s0*v0 - 2*f2*k2*l2*q2*s0*v0 + 
   4*d2*k2*m2*q2*s0*v0 - 2*c2*k2*n2*q2*s0*v0 - 2*d2*k2*l2*r2*s0*v0 - 
   2*c2*d0*j2*k2*s2*v0 - 2*c0*d2*j2*k2*s2*v0 + 2*c2*d0*i2*p2*s2*v0 + 
   2*c0*d2*i2*p2*s2*v0 + 2*f0*pow(l2,2)*p2*s2*v0 - 2*d0*l2*m2*p2*s2*v0 - 
   2*c2*l2*n0*p2*s2*v0 - 2*c0*l2*n2*p2*s2*v0 - 2*f0*k2*l2*q2*s2*v0 + 
   4*d0*k2*m2*q2*s2*v0 - 2*c2*k2*n0*q2*s2*v0 - 2*c0*k2*n2*q2*s2*v0 - 
   2*d2*k2*l2*r0*s2*v0 - 2*d0*k2*l2*r2*s2*v0 + 4*c2*k2*l2*s0*s2*v0 + 
   2*c0*k2*l2*pow(s2,2)*v0 + 2*d2*f2*pow(j2,2)*u0*v0 - 
   2*c2*g2*pow(j2,2)*u0*v0 - 2*d2*f2*i2*o2*u0*v0 + 2*c2*g2*i2*o2*u0*v0 - 
   2*g2*l2*m2*o2*u0*v0 + 2*f2*l2*n2*o2*u0*v0 + 2*d2*m2*n2*o2*u0*v0 - 
   2*c2*pow(n2,2)*o2*u0*v0 + 2*g2*j2*m2*q2*u0*v0 - 2*f2*j2*n2*q2*u0*v0 + 
   2*g2*j2*l2*r2*u0*v0 - 2*d2*j2*n2*r2*u0*v0 - 2*g2*i2*q2*r2*u0*v0 + 
   2*pow(n2,2)*q2*r2*u0*v0 - 2*f2*j2*l2*s2*u0*v0 - 2*d2*j2*m2*s2*u0*v0 + 
   4*c2*j2*n2*s2*u0*v0 + 2*f2*i2*q2*s2*u0*v0 - 2*m2*n2*q2*s2*u0*v0 + 
   2*d2*i2*r2*s2*u0*v0 - 2*l2*n2*r2*s2*u0*v0 - 2*c2*i2*pow(s2,2)*u0*v0 + 
   2*l2*m2*pow(s2,2)*u0*v0 + 2*d2*f0*pow(j2,2)*u2*v0 + 
   2*d0*f2*pow(j2,2)*u2*v0 - 2*c2*g0*pow(j2,2)*u2*v0 - 
   2*c0*g2*pow(j2,2)*u2*v0 - 2*d2*f0*i2*o2*u2*v0 - 2*d0*f2*i2*o2*u2*v0 + 
   2*c2*g0*i2*o2*u2*v0 + 2*c0*g2*i2*o2*u2*v0 - 2*g0*l2*m2*o2*u2*v0 + 
   2*f2*l2*n0*o2*u2*v0 + 2*d2*m2*n0*o2*u2*v0 + 2*f0*l2*n2*o2*u2*v0 + 
   2*d0*m2*n2*o2*u2*v0 - 4*c2*n0*n2*o2*u2*v0 - 2*c0*pow(n2,2)*o2*u2*v0 + 
   2*g0*j2*m2*q2*u2*v0 - 2*f2*j2*n0*q2*u2*v0 - 2*f0*j2*n2*q2*u2*v0 + 
   2*g2*j2*l2*r0*u2*v0 - 2*d2*j2*n2*r0*u2*v0 - 2*g2*i2*q2*r0*u2*v0 + 
   2*pow(n2,2)*q2*r0*u2*v0 + 2*g0*j2*l2*r2*u2*v0 - 2*d2*j2*n0*r2*u2*v0 - 
   2*d0*j2*n2*r2*u2*v0 - 2*g0*i2*q2*r2*u2*v0 + 4*n0*n2*q2*r2*u2*v0 - 
   2*f2*j2*l2*s0*u2*v0 - 2*d2*j2*m2*s0*u2*v0 + 4*c2*j2*n2*s0*u2*v0 + 
   2*f2*i2*q2*s0*u2*v0 - 2*m2*n2*q2*s0*u2*v0 + 2*d2*i2*r2*s0*u2*v0 - 
   2*l2*n2*r2*s0*u2*v0 - 2*f0*j2*l2*s2*u2*v0 - 2*d0*j2*m2*s2*u2*v0 + 
   4*c2*j2*n0*s2*u2*v0 + 4*c0*j2*n2*s2*u2*v0 + 2*f0*i2*q2*s2*u2*v0 - 
   2*m2*n0*q2*s2*u2*v0 + 2*d2*i2*r0*s2*u2*v0 - 2*l2*n2*r0*s2*u2*v0 + 
   2*d0*i2*r2*s2*u2*v0 - 2*l2*n0*r2*s2*u2*v0 - 4*c2*i2*s0*s2*u2*v0 + 
   4*l2*m2*s0*s2*u2*v0 - 2*c0*i2*pow(s2,2)*u2*v0 - 
   pow(d2,2)*pow(j2,2)*pow(v0,2) + pow(d2,2)*i2*o2*pow(v0,2) + 
   g2*pow(l2,2)*o2*pow(v0,2) - 2*d2*l2*n2*o2*pow(v0,2) - 
   2*g2*j2*l2*q2*pow(v0,2) + 2*d2*j2*n2*q2*pow(v0,2) + 
   g2*i2*pow(q2,2)*pow(v0,2) - pow(n2,2)*pow(q2,2)*pow(v0,2) + 
   2*d2*j2*l2*s2*pow(v0,2) - 2*d2*i2*q2*s2*pow(v0,2) + 
   2*l2*n2*q2*s2*pow(v0,2) - pow(l2,2)*pow(s2,2)*pow(v0,2) + 
   2*d0*f0*k2*l2*o2*v2 - 2*c0*g0*k2*l2*o2*v2 - 2*pow(d0,2)*k2*m2*o2*v2 + 
   2*c2*d0*k2*n0*o2*v2 + 2*c0*d2*k2*n0*o2*v2 + 2*c0*d0*k2*n2*o2*v2 - 
   2*d0*f0*j2*l2*p2*v2 + 2*c0*g0*j2*l2*p2*v2 + 2*pow(d0,2)*j2*m2*p2*v2 - 
   2*c2*d0*j2*n0*p2*v2 - 2*c0*d2*j2*n0*p2*v2 - 2*c0*d0*j2*n2*p2*v2 - 
   2*d0*f0*j2*k2*q2*v2 + 2*c0*g0*j2*k2*q2*v2 + 2*d0*f0*i2*p2*q2*v2 - 
   2*c0*g0*i2*p2*q2*v2 - 2*f0*l2*n0*p2*q2*v2 - 2*d0*m2*n0*p2*q2*v2 + 
   2*c2*pow(n0,2)*p2*q2*v2 + 4*c0*n0*n2*p2*q2*v2 + 
   2*f0*k2*n0*pow(q2,2)*v2 + 4*d0*d2*j2*k2*r0*v2 - 4*d0*d2*i2*p2*r0*v2 - 
   2*g0*pow(l2,2)*p2*r0*v2 + 4*d2*l2*n0*p2*r0*v2 + 4*d0*l2*n2*p2*r0*v2 + 
   2*g0*k2*l2*q2*r0*v2 - 2*d2*k2*n0*q2*r0*v2 - 2*d0*k2*n2*q2*r0*v2 + 
   2*pow(d0,2)*j2*k2*r2*v2 - 2*pow(d0,2)*i2*p2*r2*v2 + 
   4*d0*l2*n0*p2*r2*v2 - 2*d0*k2*n0*q2*r2*v2 - 2*c2*d0*j2*k2*s0*v2 - 
   2*c0*d2*j2*k2*s0*v2 + 2*c2*d0*i2*p2*s0*v2 + 2*c0*d2*i2*p2*s0*v2 + 
   2*f0*pow(l2,2)*p2*s0*v2 - 2*d0*l2*m2*p2*s0*v2 - 2*c2*l2*n0*p2*s0*v2 - 
   2*c0*l2*n2*p2*s0*v2 - 2*f0*k2*l2*q2*s0*v2 + 4*d0*k2*m2*q2*s0*v2 - 
   2*c2*k2*n0*q2*s0*v2 - 2*c0*k2*n2*q2*s0*v2 - 2*d2*k2*l2*r0*s0*v2 - 
   2*d0*k2*l2*r2*s0*v2 + 2*c2*k2*l2*pow(s0,2)*v2 - 2*c0*d0*j2*k2*s2*v2 + 
   2*c0*d0*i2*p2*s2*v2 - 2*c0*l2*n0*p2*s2*v2 - 2*c0*k2*n0*q2*s2*v2 - 
   2*d0*k2*l2*r0*s2*v2 + 4*c0*k2*l2*s0*s2*v2 + 2*d2*f0*pow(j2,2)*u0*v2 + 
   2*d0*f2*pow(j2,2)*u0*v2 - 2*c2*g0*pow(j2,2)*u0*v2 - 
   2*c0*g2*pow(j2,2)*u0*v2 - 2*d2*f0*i2*o2*u0*v2 - 2*d0*f2*i2*o2*u0*v2 + 
   2*c2*g0*i2*o2*u0*v2 + 2*c0*g2*i2*o2*u0*v2 - 2*g0*l2*m2*o2*u0*v2 + 
   2*f2*l2*n0*o2*u0*v2 + 2*d2*m2*n0*o2*u0*v2 + 2*f0*l2*n2*o2*u0*v2 + 
   2*d0*m2*n2*o2*u0*v2 - 4*c2*n0*n2*o2*u0*v2 - 2*c0*pow(n2,2)*o2*u0*v2 + 
   2*g0*j2*m2*q2*u0*v2 - 2*f2*j2*n0*q2*u0*v2 - 2*f0*j2*n2*q2*u0*v2 + 
   2*g2*j2*l2*r0*u0*v2 - 2*d2*j2*n2*r0*u0*v2 - 2*g2*i2*q2*r0*u0*v2 + 
   2*pow(n2,2)*q2*r0*u0*v2 + 2*g0*j2*l2*r2*u0*v2 - 2*d2*j2*n0*r2*u0*v2 - 
   2*d0*j2*n2*r2*u0*v2 - 2*g0*i2*q2*r2*u0*v2 + 4*n0*n2*q2*r2*u0*v2 - 
   2*f2*j2*l2*s0*u0*v2 - 2*d2*j2*m2*s0*u0*v2 + 4*c2*j2*n2*s0*u0*v2 + 
   2*f2*i2*q2*s0*u0*v2 - 2*m2*n2*q2*s0*u0*v2 + 2*d2*i2*r2*s0*u0*v2 - 
   2*l2*n2*r2*s0*u0*v2 - 2*f0*j2*l2*s2*u0*v2 - 2*d0*j2*m2*s2*u0*v2 + 
   4*c2*j2*n0*s2*u0*v2 + 4*c0*j2*n2*s2*u0*v2 + 2*f0*i2*q2*s2*u0*v2 - 
   2*m2*n0*q2*s2*u0*v2 + 2*d2*i2*r0*s2*u0*v2 - 2*l2*n2*r0*s2*u0*v2 + 
   2*d0*i2*r2*s2*u0*v2 - 2*l2*n0*r2*s2*u0*v2 - 4*c2*i2*s0*s2*u0*v2 + 
   4*l2*m2*s0*s2*u0*v2 - 2*c0*i2*pow(s2,2)*u0*v2 + 
   2*d0*f0*pow(j2,2)*u2*v2 - 2*c0*g0*pow(j2,2)*u2*v2 - 
   2*d0*f0*i2*o2*u2*v2 + 2*c0*g0*i2*o2*u2*v2 + 2*f0*l2*n0*o2*u2*v2 + 
   2*d0*m2*n0*o2*u2*v2 - 2*c2*pow(n0,2)*o2*u2*v2 - 4*c0*n0*n2*o2*u2*v2 - 
   2*f0*j2*n0*q2*u2*v2 + 2*g0*j2*l2*r0*u2*v2 - 2*d2*j2*n0*r0*u2*v2 - 
   2*d0*j2*n2*r0*u2*v2 - 2*g0*i2*q2*r0*u2*v2 + 4*n0*n2*q2*r0*u2*v2 - 
   2*d0*j2*n0*r2*u2*v2 + 2*pow(n0,2)*q2*r2*u2*v2 - 2*f0*j2*l2*s0*u2*v2 - 
   2*d0*j2*m2*s0*u2*v2 + 4*c2*j2*n0*s0*u2*v2 + 4*c0*j2*n2*s0*u2*v2 + 
   2*f0*i2*q2*s0*u2*v2 - 2*m2*n0*q2*s0*u2*v2 + 2*d2*i2*r0*s0*u2*v2 - 
   2*l2*n2*r0*s0*u2*v2 + 2*d0*i2*r2*s0*u2*v2 - 2*l2*n0*r2*s0*u2*v2 - 
   2*c2*i2*pow(s0,2)*u2*v2 + 2*l2*m2*pow(s0,2)*u2*v2 + 
   4*c0*j2*n0*s2*u2*v2 + 2*d0*i2*r0*s2*u2*v2 - 2*l2*n0*r0*s2*u2*v2 - 
   4*c0*i2*s0*s2*u2*v2 - 4*d0*d2*pow(j2,2)*v0*v2 + 4*d0*d2*i2*o2*v0*v2 + 
   2*g0*pow(l2,2)*o2*v0*v2 - 4*d2*l2*n0*o2*v0*v2 - 4*d0*l2*n2*o2*v0*v2 - 
   4*g0*j2*l2*q2*v0*v2 + 4*d2*j2*n0*q2*v0*v2 + 4*d0*j2*n2*q2*v0*v2 + 
   2*g0*i2*pow(q2,2)*v0*v2 - 4*n0*n2*pow(q2,2)*v0*v2 + 
   4*d2*j2*l2*s0*v0*v2 - 4*d2*i2*q2*s0*v0*v2 + 4*l2*n2*q2*s0*v0*v2 + 
   4*d0*j2*l2*s2*v0*v2 - 4*d0*i2*q2*s2*v0*v2 + 4*l2*n0*q2*s2*v0*v2 - 
   4*pow(l2,2)*s0*s2*v0*v2 - pow(d0,2)*pow(j2,2)*pow(v2,2) + 
   pow(d0,2)*i2*o2*pow(v2,2) - 2*d0*l2*n0*o2*pow(v2,2) + 
   2*d0*j2*n0*q2*pow(v2,2) - pow(n0,2)*pow(q2,2)*pow(v2,2) + 
   2*d0*j2*l2*s0*pow(v2,2) - 2*d0*i2*q2*s0*pow(v2,2) + 
   2*l2*n0*q2*s0*pow(v2,2) - pow(l2,2)*pow(s0,2)*pow(v2,2) - 
   2*d2*e0*k2*l2*o2*w0 - 2*d0*e2*k2*l2*o2*w0 + 2*c2*f0*k2*l2*o2*w0 + 
   2*c0*f2*k2*l2*o2*w0 + 2*c2*d0*k2*m2*o2*w0 + 2*c0*d2*k2*m2*o2*w0 - 
   2*pow(c2,2)*k2*n0*o2*w0 - 4*c0*c2*k2*n2*o2*w0 + 2*d2*e0*j2*l2*p2*w0 + 
   2*d0*e2*j2*l2*p2*w0 - 2*c2*f0*j2*l2*p2*w0 - 2*c0*f2*j2*l2*p2*w0 - 
   2*c2*d0*j2*m2*p2*w0 - 2*c0*d2*j2*m2*p2*w0 + 2*pow(c2,2)*j2*n0*p2*w0 + 
   4*c0*c2*j2*n2*p2*w0 + 2*d2*e0*j2*k2*q2*w0 + 2*d0*e2*j2*k2*q2*w0 - 
   2*c2*f0*j2*k2*q2*w0 - 2*c0*f2*j2*k2*q2*w0 - 2*d2*e0*i2*p2*q2*w0 - 
   2*d0*e2*i2*p2*q2*w0 + 2*c2*f0*i2*p2*q2*w0 + 2*c0*f2*i2*p2*q2*w0 - 
   2*f0*l2*m2*p2*q2*w0 + 2*d0*pow(m2,2)*p2*q2*w0 + 2*e2*l2*n0*p2*q2*w0 - 
   2*c2*m2*n0*p2*q2*w0 + 2*e0*l2*n2*p2*q2*w0 - 2*c0*m2*n2*p2*q2*w0 + 
   2*f0*k2*m2*pow(q2,2)*w0 - 2*e2*k2*n0*pow(q2,2)*w0 - 
   2*e0*k2*n2*pow(q2,2)*w0 - 2*c2*d2*j2*k2*r0*w0 + 2*c2*d2*i2*p2*r0*w0 + 
   2*f2*pow(l2,2)*p2*r0*w0 - 2*d2*l2*m2*p2*r0*w0 - 2*c2*l2*n2*p2*r0*w0 - 
   2*f2*k2*l2*q2*r0*w0 - 2*d2*k2*m2*q2*r0*w0 + 4*c2*k2*n2*q2*r0*w0 - 
   2*c2*d0*j2*k2*r2*w0 - 2*c0*d2*j2*k2*r2*w0 + 2*c2*d0*i2*p2*r2*w0 + 
   2*c0*d2*i2*p2*r2*w0 + 2*f0*pow(l2,2)*p2*r2*w0 - 2*d0*l2*m2*p2*r2*w0 - 
   2*c2*l2*n0*p2*r2*w0 - 2*c0*l2*n2*p2*r2*w0 - 2*f0*k2*l2*q2*r2*w0 - 
   2*d0*k2*m2*q2*r2*w0 + 4*c2*k2*n0*q2*r2*w0 + 4*c0*k2*n2*q2*r2*w0 + 
   4*d2*k2*l2*r0*r2*w0 + 2*d0*k2*l2*pow(r2,2)*w0 + 
   2*pow(c2,2)*j2*k2*s0*w0 - 2*pow(c2,2)*i2*p2*s0*w0 - 
   2*e2*pow(l2,2)*p2*s0*w0 + 4*c2*l2*m2*p2*s0*w0 + 2*e2*k2*l2*q2*s0*w0 - 
   2*c2*k2*m2*q2*s0*w0 - 2*c2*k2*l2*r2*s0*w0 + 4*c0*c2*j2*k2*s2*w0 - 
   4*c0*c2*i2*p2*s2*w0 - 2*e0*pow(l2,2)*p2*s2*w0 + 4*c0*l2*m2*p2*s2*w0 + 
   2*e0*k2*l2*q2*s2*w0 - 2*c0*k2*m2*q2*s2*w0 - 2*c2*k2*l2*r0*s2*w0 - 
   2*c0*k2*l2*r2*s2*w0 - 2*d2*e2*pow(j2,2)*u0*w0 + 
   2*c2*f2*pow(j2,2)*u0*w0 + 2*d2*e2*i2*o2*u0*w0 - 2*c2*f2*i2*o2*u0*w0 + 
   2*f2*l2*m2*o2*u0*w0 - 2*d2*pow(m2,2)*o2*u0*w0 - 2*e2*l2*n2*o2*u0*w0 + 
   2*c2*m2*n2*o2*u0*w0 - 2*f2*j2*m2*q2*u0*w0 + 2*e2*j2*n2*q2*u0*w0 - 
   2*f2*j2*l2*r2*u0*w0 + 4*d2*j2*m2*r2*u0*w0 - 2*c2*j2*n2*r2*u0*w0 + 
   2*f2*i2*q2*r2*u0*w0 - 2*m2*n2*q2*r2*u0*w0 - 2*d2*i2*pow(r2,2)*u0*w0 + 
   2*l2*n2*pow(r2,2)*u0*w0 + 2*e2*j2*l2*s2*u0*w0 - 2*c2*j2*m2*s2*u0*w0 - 
   2*e2*i2*q2*s2*u0*w0 + 2*pow(m2,2)*q2*s2*u0*w0 + 2*c2*i2*r2*s2*u0*w0 - 
   2*l2*m2*r2*s2*u0*w0 - 2*d2*e0*pow(j2,2)*u2*w0 - 
   2*d0*e2*pow(j2,2)*u2*w0 + 2*c2*f0*pow(j2,2)*u2*w0 + 
   2*c0*f2*pow(j2,2)*u2*w0 + 2*d2*e0*i2*o2*u2*w0 + 2*d0*e2*i2*o2*u2*w0 - 
   2*c2*f0*i2*o2*u2*w0 - 2*c0*f2*i2*o2*u2*w0 + 2*f0*l2*m2*o2*u2*w0 - 
   2*d0*pow(m2,2)*o2*u2*w0 - 2*e2*l2*n0*o2*u2*w0 + 2*c2*m2*n0*o2*u2*w0 - 
   2*e0*l2*n2*o2*u2*w0 + 2*c0*m2*n2*o2*u2*w0 - 2*f0*j2*m2*q2*u2*w0 + 
   2*e2*j2*n0*q2*u2*w0 + 2*e0*j2*n2*q2*u2*w0 - 2*f2*j2*l2*r0*u2*w0 + 
   4*d2*j2*m2*r0*u2*w0 - 2*c2*j2*n2*r0*u2*w0 + 2*f2*i2*q2*r0*u2*w0 - 
   2*m2*n2*q2*r0*u2*w0 - 2*f0*j2*l2*r2*u2*w0 + 4*d0*j2*m2*r2*u2*w0 - 
   2*c2*j2*n0*r2*u2*w0 - 2*c0*j2*n2*r2*u2*w0 + 2*f0*i2*q2*r2*u2*w0 - 
   2*m2*n0*q2*r2*u2*w0 - 4*d2*i2*r0*r2*u2*w0 + 4*l2*n2*r0*r2*u2*w0 - 
   2*d0*i2*pow(r2,2)*u2*w0 + 2*l2*n0*pow(r2,2)*u2*w0 + 
   2*e2*j2*l2*s0*u2*w0 - 2*c2*j2*m2*s0*u2*w0 - 2*e2*i2*q2*s0*u2*w0 + 
   2*pow(m2,2)*q2*s0*u2*w0 + 2*c2*i2*r2*s0*u2*w0 - 2*l2*m2*r2*s0*u2*w0 + 
   2*e0*j2*l2*s2*u2*w0 - 2*c0*j2*m2*s2*u2*w0 - 2*e0*i2*q2*s2*u2*w0 + 
   2*c2*i2*r0*s2*u2*w0 - 2*l2*m2*r0*s2*u2*w0 + 2*c0*i2*r2*s2*u2*w0 + 
   2*c2*d2*pow(j2,2)*v0*w0 - 2*c2*d2*i2*o2*v0*w0 - 
   2*f2*pow(l2,2)*o2*v0*w0 + 2*d2*l2*m2*o2*v0*w0 + 2*c2*l2*n2*o2*v0*w0 + 
   4*f2*j2*l2*q2*v0*w0 - 2*d2*j2*m2*q2*v0*w0 - 2*c2*j2*n2*q2*v0*w0 - 
   2*f2*i2*pow(q2,2)*v0*w0 + 2*m2*n2*pow(q2,2)*v0*w0 - 
   2*d2*j2*l2*r2*v0*w0 + 2*d2*i2*q2*r2*v0*w0 - 2*l2*n2*q2*r2*v0*w0 - 
   2*c2*j2*l2*s2*v0*w0 + 2*c2*i2*q2*s2*v0*w0 - 2*l2*m2*q2*s2*v0*w0 + 
   2*pow(l2,2)*r2*s2*v0*w0 + 2*c2*d0*pow(j2,2)*v2*w0 + 
   2*c0*d2*pow(j2,2)*v2*w0 - 2*c2*d0*i2*o2*v2*w0 - 2*c0*d2*i2*o2*v2*w0 - 
   2*f0*pow(l2,2)*o2*v2*w0 + 2*d0*l2*m2*o2*v2*w0 + 2*c2*l2*n0*o2*v2*w0 + 
   2*c0*l2*n2*o2*v2*w0 + 4*f0*j2*l2*q2*v2*w0 - 2*d0*j2*m2*q2*v2*w0 - 
   2*c2*j2*n0*q2*v2*w0 - 2*c0*j2*n2*q2*v2*w0 - 2*f0*i2*pow(q2,2)*v2*w0 + 
   2*m2*n0*pow(q2,2)*v2*w0 - 2*d2*j2*l2*r0*v2*w0 + 2*d2*i2*q2*r0*v2*w0 - 
   2*l2*n2*q2*r0*v2*w0 - 2*d0*j2*l2*r2*v2*w0 + 2*d0*i2*q2*r2*v2*w0 - 
   2*l2*n0*q2*r2*v2*w0 - 2*c2*j2*l2*s0*v2*w0 + 2*c2*i2*q2*s0*v2*w0 - 
   2*l2*m2*q2*s0*v2*w0 + 2*pow(l2,2)*r2*s0*v2*w0 - 2*c0*j2*l2*s2*v2*w0 + 
   2*c0*i2*q2*s2*v2*w0 + 2*pow(l2,2)*r0*s2*v2*w0 - 
   pow(c2,2)*pow(j2,2)*pow(w0,2) + pow(c2,2)*i2*o2*pow(w0,2) + 
   e2*pow(l2,2)*o2*pow(w0,2) - 2*c2*l2*m2*o2*pow(w0,2) - 
   2*e2*j2*l2*q2*pow(w0,2) + 2*c2*j2*m2*q2*pow(w0,2) + 
   e2*i2*pow(q2,2)*pow(w0,2) - pow(m2,2)*pow(q2,2)*pow(w0,2) + 
   2*c2*j2*l2*r2*pow(w0,2) - 2*c2*i2*q2*r2*pow(w0,2) + 
   2*l2*m2*q2*r2*pow(w0,2) - pow(l2,2)*pow(r2,2)*pow(w0,2) - 
   2*d0*e0*k2*l2*o2*w2 + 2*c0*f0*k2*l2*o2*w2 + 2*c0*d0*k2*m2*o2*w2 - 
   4*c0*c2*k2*n0*o2*w2 - 2*pow(c0,2)*k2*n2*o2*w2 + 2*d0*e0*j2*l2*p2*w2 - 
   2*c0*f0*j2*l2*p2*w2 - 2*c0*d0*j2*m2*p2*w2 + 4*c0*c2*j2*n0*p2*w2 + 
   2*pow(c0,2)*j2*n2*p2*w2 + 2*d0*e0*j2*k2*q2*w2 - 2*c0*f0*j2*k2*q2*w2 - 
   2*d0*e0*i2*p2*q2*w2 + 2*c0*f0*i2*p2*q2*w2 + 2*e0*l2*n0*p2*q2*w2 - 
   2*c0*m2*n0*p2*q2*w2 - 2*e0*k2*n0*pow(q2,2)*w2 - 2*c2*d0*j2*k2*r0*w2 - 
   2*c0*d2*j2*k2*r0*w2 + 2*c2*d0*i2*p2*r0*w2 + 2*c0*d2*i2*p2*r0*w2 + 
   2*f0*pow(l2,2)*p2*r0*w2 - 2*d0*l2*m2*p2*r0*w2 - 2*c2*l2*n0*p2*r0*w2 - 
   2*c0*l2*n2*p2*r0*w2 - 2*f0*k2*l2*q2*r0*w2 - 2*d0*k2*m2*q2*r0*w2 + 
   4*c2*k2*n0*q2*r0*w2 + 4*c0*k2*n2*q2*r0*w2 + 2*d2*k2*l2*pow(r0,2)*w2 - 
   2*c0*d0*j2*k2*r2*w2 + 2*c0*d0*i2*p2*r2*w2 - 2*c0*l2*n0*p2*r2*w2 + 
   4*c0*k2*n0*q2*r2*w2 + 4*d0*k2*l2*r0*r2*w2 + 4*c0*c2*j2*k2*s0*w2 - 
   4*c0*c2*i2*p2*s0*w2 - 2*e0*pow(l2,2)*p2*s0*w2 + 4*c0*l2*m2*p2*s0*w2 + 
   2*e0*k2*l2*q2*s0*w2 - 2*c0*k2*m2*q2*s0*w2 - 2*c2*k2*l2*r0*s0*w2 - 
   2*c0*k2*l2*r2*s0*w2 + 2*pow(c0,2)*j2*k2*s2*w2 - 
   2*pow(c0,2)*i2*p2*s2*w2 - 2*c0*k2*l2*r0*s2*w2 - 
   2*d2*e0*pow(j2,2)*u0*w2 - 2*d0*e2*pow(j2,2)*u0*w2 + 
   2*c2*f0*pow(j2,2)*u0*w2 + 2*c0*f2*pow(j2,2)*u0*w2 + 
   2*d2*e0*i2*o2*u0*w2 + 2*d0*e2*i2*o2*u0*w2 - 2*c2*f0*i2*o2*u0*w2 - 
   2*c0*f2*i2*o2*u0*w2 + 2*f0*l2*m2*o2*u0*w2 - 2*d0*pow(m2,2)*o2*u0*w2 - 
   2*e2*l2*n0*o2*u0*w2 + 2*c2*m2*n0*o2*u0*w2 - 2*e0*l2*n2*o2*u0*w2 + 
   2*c0*m2*n2*o2*u0*w2 - 2*f0*j2*m2*q2*u0*w2 + 2*e2*j2*n0*q2*u0*w2 + 
   2*e0*j2*n2*q2*u0*w2 - 2*f2*j2*l2*r0*u0*w2 + 4*d2*j2*m2*r0*u0*w2 - 
   2*c2*j2*n2*r0*u0*w2 + 2*f2*i2*q2*r0*u0*w2 - 2*m2*n2*q2*r0*u0*w2 - 
   2*f0*j2*l2*r2*u0*w2 + 4*d0*j2*m2*r2*u0*w2 - 2*c2*j2*n0*r2*u0*w2 - 
   2*c0*j2*n2*r2*u0*w2 + 2*f0*i2*q2*r2*u0*w2 - 2*m2*n0*q2*r2*u0*w2 - 
   4*d2*i2*r0*r2*u0*w2 + 4*l2*n2*r0*r2*u0*w2 - 2*d0*i2*pow(r2,2)*u0*w2 + 
   2*l2*n0*pow(r2,2)*u0*w2 + 2*e2*j2*l2*s0*u0*w2 - 2*c2*j2*m2*s0*u0*w2 - 
   2*e2*i2*q2*s0*u0*w2 + 2*pow(m2,2)*q2*s0*u0*w2 + 2*c2*i2*r2*s0*u0*w2 - 
   2*l2*m2*r2*s0*u0*w2 + 2*e0*j2*l2*s2*u0*w2 - 2*c0*j2*m2*s2*u0*w2 - 
   2*e0*i2*q2*s2*u0*w2 + 2*c2*i2*r0*s2*u0*w2 - 2*l2*m2*r0*s2*u0*w2 + 
   2*c0*i2*r2*s2*u0*w2 - 2*d0*e0*pow(j2,2)*u2*w2 + 
   2*c0*f0*pow(j2,2)*u2*w2 + 2*d0*e0*i2*o2*u2*w2 - 2*c0*f0*i2*o2*u2*w2 - 
   2*e0*l2*n0*o2*u2*w2 + 2*c0*m2*n0*o2*u2*w2 + 2*e0*j2*n0*q2*u2*w2 - 
   2*f0*j2*l2*r0*u2*w2 + 4*d0*j2*m2*r0*u2*w2 - 2*c2*j2*n0*r0*u2*w2 - 
   2*c0*j2*n2*r0*u2*w2 + 2*f0*i2*q2*r0*u2*w2 - 2*m2*n0*q2*r0*u2*w2 - 
   2*d2*i2*pow(r0,2)*u2*w2 + 2*l2*n2*pow(r0,2)*u2*w2 - 
   2*c0*j2*n0*r2*u2*w2 - 4*d0*i2*r0*r2*u2*w2 + 4*l2*n0*r0*r2*u2*w2 + 
   2*e0*j2*l2*s0*u2*w2 - 2*c0*j2*m2*s0*u2*w2 - 2*e0*i2*q2*s0*u2*w2 + 
   2*c2*i2*r0*s0*u2*w2 - 2*l2*m2*r0*s0*u2*w2 + 2*c0*i2*r2*s0*u2*w2 + 
   2*c0*i2*r0*s2*u2*w2 + 2*c2*d0*pow(j2,2)*v0*w2 + 
   2*c0*d2*pow(j2,2)*v0*w2 - 2*c2*d0*i2*o2*v0*w2 - 2*c0*d2*i2*o2*v0*w2 - 
   2*f0*pow(l2,2)*o2*v0*w2 + 2*d0*l2*m2*o2*v0*w2 + 2*c2*l2*n0*o2*v0*w2 + 
   2*c0*l2*n2*o2*v0*w2 + 4*f0*j2*l2*q2*v0*w2 - 2*d0*j2*m2*q2*v0*w2 - 
   2*c2*j2*n0*q2*v0*w2 - 2*c0*j2*n2*q2*v0*w2 - 2*f0*i2*pow(q2,2)*v0*w2 + 
   2*m2*n0*pow(q2,2)*v0*w2 - 2*d2*j2*l2*r0*v0*w2 + 2*d2*i2*q2*r0*v0*w2 - 
   2*l2*n2*q2*r0*v0*w2 - 2*d0*j2*l2*r2*v0*w2 + 2*d0*i2*q2*r2*v0*w2 - 
   2*l2*n0*q2*r2*v0*w2 - 2*c2*j2*l2*s0*v0*w2 + 2*c2*i2*q2*s0*v0*w2 - 
   2*l2*m2*q2*s0*v0*w2 + 2*pow(l2,2)*r2*s0*v0*w2 - 2*c0*j2*l2*s2*v0*w2 + 
   2*c0*i2*q2*s2*v0*w2 + 2*pow(l2,2)*r0*s2*v0*w2 + 
   2*c0*d0*pow(j2,2)*v2*w2 - 2*c0*d0*i2*o2*v2*w2 + 2*c0*l2*n0*o2*v2*w2 - 
   2*c0*j2*n0*q2*v2*w2 - 2*d0*j2*l2*r0*v2*w2 + 2*d0*i2*q2*r0*v2*w2 - 
   2*l2*n0*q2*r0*v2*w2 - 2*c0*j2*l2*s0*v2*w2 + 2*c0*i2*q2*s0*v2*w2 + 
   2*pow(l2,2)*r0*s0*v2*w2 - 4*c0*c2*pow(j2,2)*w0*w2 + 
   4*c0*c2*i2*o2*w0*w2 + 2*e0*pow(l2,2)*o2*w0*w2 - 4*c0*l2*m2*o2*w0*w2 - 
   4*e0*j2*l2*q2*w0*w2 + 4*c0*j2*m2*q2*w0*w2 + 2*e0*i2*pow(q2,2)*w0*w2 + 
   4*c2*j2*l2*r0*w0*w2 - 4*c2*i2*q2*r0*w0*w2 + 4*l2*m2*q2*r0*w0*w2 + 
   4*c0*j2*l2*r2*w0*w2 - 4*c0*i2*q2*r2*w0*w2 - 4*pow(l2,2)*r0*r2*w0*w2 - 
   pow(c0,2)*pow(j2,2)*pow(w2,2) + pow(c0,2)*i2*o2*pow(w2,2) + 
   2*c0*j2*l2*r0*pow(w2,2) - 2*c0*i2*q2*r0*pow(w2,2) - 
   pow(l2,2)*pow(r0,2)*pow(w2,2) + 2*f0*f2*pow(k2,2)*o2*z0 - 
   e2*g0*pow(k2,2)*o2*z0 - e0*g2*pow(k2,2)*o2*z0 - 4*f0*f2*j2*k2*p2*z0 + 
   2*e2*g0*j2*k2*p2*z0 + 2*e0*g2*j2*k2*p2*z0 + 2*f0*f2*i2*pow(p2,2)*z0 - 
   e2*g0*i2*pow(p2,2)*z0 - e0*g2*i2*pow(p2,2)*z0 + 
   g0*pow(m2,2)*pow(p2,2)*z0 - 2*f2*m2*n0*pow(p2,2)*z0 - 
   2*f0*m2*n2*pow(p2,2)*z0 + 2*e2*n0*n2*pow(p2,2)*z0 + 
   e0*pow(n2,2)*pow(p2,2)*z0 - 2*g2*k2*m2*p2*r0*z0 + 
   2*f2*k2*n2*p2*r0*z0 - 2*g0*k2*m2*p2*r2*z0 + 2*f2*k2*n0*p2*r2*z0 + 
   2*f0*k2*n2*p2*r2*z0 + 2*g2*pow(k2,2)*r0*r2*z0 + 
   g0*pow(k2,2)*pow(r2,2)*z0 + 2*f2*k2*m2*p2*s0*z0 - 
   2*e2*k2*n2*p2*s0*z0 - 2*f2*pow(k2,2)*r2*s0*z0 + 2*f0*k2*m2*p2*s2*z0 - 
   2*e2*k2*n0*p2*s2*z0 - 2*e0*k2*n2*p2*s2*z0 - 2*f2*pow(k2,2)*r0*s2*z0 - 
   2*f0*pow(k2,2)*r2*s2*z0 + 2*e2*pow(k2,2)*s0*s2*z0 + 
   e0*pow(k2,2)*pow(s2,2)*z0 + 2*f0*f2*pow(j2,2)*t2*z0 - 
   e2*g0*pow(j2,2)*t2*z0 - e0*g2*pow(j2,2)*t2*z0 - 2*f0*f2*i2*o2*t2*z0 + 
   e2*g0*i2*o2*t2*z0 + e0*g2*i2*o2*t2*z0 - g0*pow(m2,2)*o2*t2*z0 + 
   2*f2*m2*n0*o2*t2*z0 + 2*f0*m2*n2*o2*t2*z0 - 2*e2*n0*n2*o2*t2*z0 - 
   e0*pow(n2,2)*o2*t2*z0 + 2*g2*j2*m2*r0*t2*z0 - 2*f2*j2*n2*r0*t2*z0 + 
   2*g0*j2*m2*r2*t2*z0 - 2*f2*j2*n0*r2*t2*z0 - 2*f0*j2*n2*r2*t2*z0 - 
   2*g2*i2*r0*r2*t2*z0 + 2*pow(n2,2)*r0*r2*t2*z0 - 
   g0*i2*pow(r2,2)*t2*z0 + 2*n0*n2*pow(r2,2)*t2*z0 - 
   2*f2*j2*m2*s0*t2*z0 + 2*e2*j2*n2*s0*t2*z0 + 2*f2*i2*r2*s0*t2*z0 - 
   2*m2*n2*r2*s0*t2*z0 - 2*f0*j2*m2*s2*t2*z0 + 2*e2*j2*n0*s2*t2*z0 + 
   2*e0*j2*n2*s2*t2*z0 + 2*f2*i2*r0*s2*t2*z0 - 2*m2*n2*r0*s2*t2*z0 + 
   2*f0*i2*r2*s2*t2*z0 - 2*m2*n0*r2*s2*t2*z0 - 2*e2*i2*s0*s2*t2*z0 + 
   2*pow(m2,2)*s0*s2*t2*z0 - e0*i2*pow(s2,2)*t2*z0 + 
   2*g2*k2*m2*o2*v0*z0 - 2*f2*k2*n2*o2*v0*z0 - 2*g2*j2*m2*p2*v0*z0 + 
   2*f2*j2*n2*p2*v0*z0 - 2*g2*j2*k2*r2*v0*z0 + 2*g2*i2*p2*r2*v0*z0 - 
   2*pow(n2,2)*p2*r2*v0*z0 + 2*f2*j2*k2*s2*v0*z0 - 2*f2*i2*p2*s2*v0*z0 + 
   2*m2*n2*p2*s2*v0*z0 + 2*k2*n2*r2*s2*v0*z0 - 2*k2*m2*pow(s2,2)*v0*z0 + 
   2*g0*k2*m2*o2*v2*z0 - 2*f2*k2*n0*o2*v2*z0 - 2*f0*k2*n2*o2*v2*z0 - 
   2*g0*j2*m2*p2*v2*z0 + 2*f2*j2*n0*p2*v2*z0 + 2*f0*j2*n2*p2*v2*z0 - 
   2*g2*j2*k2*r0*v2*z0 + 2*g2*i2*p2*r0*v2*z0 - 2*pow(n2,2)*p2*r0*v2*z0 - 
   2*g0*j2*k2*r2*v2*z0 + 2*g0*i2*p2*r2*v2*z0 - 4*n0*n2*p2*r2*v2*z0 + 
   2*f2*j2*k2*s0*v2*z0 - 2*f2*i2*p2*s0*v2*z0 + 2*m2*n2*p2*s0*v2*z0 + 
   2*k2*n2*r2*s0*v2*z0 + 2*f0*j2*k2*s2*v2*z0 - 2*f0*i2*p2*s2*v2*z0 + 
   2*m2*n0*p2*s2*v2*z0 + 2*k2*n2*r0*s2*v2*z0 + 2*k2*n0*r2*s2*v2*z0 - 
   4*k2*m2*s0*s2*v2*z0 + 2*g2*pow(j2,2)*v0*v2*z0 - 2*g2*i2*o2*v0*v2*z0 + 
   2*pow(n2,2)*o2*v0*v2*z0 - 4*j2*n2*s2*v0*v2*z0 + 
   2*i2*pow(s2,2)*v0*v2*z0 + g0*pow(j2,2)*pow(v2,2)*z0 - 
   g0*i2*o2*pow(v2,2)*z0 + 2*n0*n2*o2*pow(v2,2)*z0 - 
   2*j2*n2*s0*pow(v2,2)*z0 - 2*j2*n0*s2*pow(v2,2)*z0 + 
   2*i2*s0*s2*pow(v2,2)*z0 - 2*f2*k2*m2*o2*w0*z0 + 2*e2*k2*n2*o2*w0*z0 + 
   2*f2*j2*m2*p2*w0*z0 - 2*e2*j2*n2*p2*w0*z0 + 2*f2*j2*k2*r2*w0*z0 - 
   2*f2*i2*p2*r2*w0*z0 + 2*m2*n2*p2*r2*w0*z0 - 2*k2*n2*pow(r2,2)*w0*z0 - 
   2*e2*j2*k2*s2*w0*z0 + 2*e2*i2*p2*s2*w0*z0 - 2*pow(m2,2)*p2*s2*w0*z0 + 
   2*k2*m2*r2*s2*w0*z0 - 2*f2*pow(j2,2)*v2*w0*z0 + 2*f2*i2*o2*v2*w0*z0 - 
   2*m2*n2*o2*v2*w0*z0 + 2*j2*n2*r2*v2*w0*z0 + 2*j2*m2*s2*v2*w0*z0 - 
   2*i2*r2*s2*v2*w0*z0 - 2*f0*k2*m2*o2*w2*z0 + 2*e2*k2*n0*o2*w2*z0 + 
   2*e0*k2*n2*o2*w2*z0 + 2*f0*j2*m2*p2*w2*z0 - 2*e2*j2*n0*p2*w2*z0 - 
   2*e0*j2*n2*p2*w2*z0 + 2*f2*j2*k2*r0*w2*z0 - 2*f2*i2*p2*r0*w2*z0 + 
   2*m2*n2*p2*r0*w2*z0 + 2*f0*j2*k2*r2*w2*z0 - 2*f0*i2*p2*r2*w2*z0 + 
   2*m2*n0*p2*r2*w2*z0 - 4*k2*n2*r0*r2*w2*z0 - 2*k2*n0*pow(r2,2)*w2*z0 - 
   2*e2*j2*k2*s0*w2*z0 + 2*e2*i2*p2*s0*w2*z0 - 2*pow(m2,2)*p2*s0*w2*z0 + 
   2*k2*m2*r2*s0*w2*z0 - 2*e0*j2*k2*s2*w2*z0 + 2*e0*i2*p2*s2*w2*z0 + 
   2*k2*m2*r0*s2*w2*z0 - 2*f2*pow(j2,2)*v0*w2*z0 + 2*f2*i2*o2*v0*w2*z0 - 
   2*m2*n2*o2*v0*w2*z0 + 2*j2*n2*r2*v0*w2*z0 + 2*j2*m2*s2*v0*w2*z0 - 
   2*i2*r2*s2*v0*w2*z0 - 2*f0*pow(j2,2)*v2*w2*z0 + 2*f0*i2*o2*v2*w2*z0 - 
   2*m2*n0*o2*v2*w2*z0 + 2*j2*n2*r0*v2*w2*z0 + 2*j2*n0*r2*v2*w2*z0 + 
   2*j2*m2*s0*v2*w2*z0 - 2*i2*r2*s0*v2*w2*z0 - 2*i2*r0*s2*v2*w2*z0 + 
   2*e2*pow(j2,2)*w0*w2*z0 - 2*e2*i2*o2*w0*w2*z0 + 
   2*pow(m2,2)*o2*w0*w2*z0 - 4*j2*m2*r2*w0*w2*z0 + 
   2*i2*pow(r2,2)*w0*w2*z0 + e0*pow(j2,2)*pow(w2,2)*z0 - 
   e0*i2*o2*pow(w2,2)*z0 - 2*j2*m2*r0*pow(w2,2)*z0 + 
   2*i2*r0*r2*pow(w2,2)*z0 + pow(f0,2)*pow(k2,2)*o2*z2 - 
   e0*g0*pow(k2,2)*o2*z2 - 2*pow(f0,2)*j2*k2*p2*z2 + 
   2*e0*g0*j2*k2*p2*z2 + pow(f0,2)*i2*pow(p2,2)*z2 - 
   e0*g0*i2*pow(p2,2)*z2 - 2*f0*m2*n0*pow(p2,2)*z2 + 
   e2*pow(n0,2)*pow(p2,2)*z2 + 2*e0*n0*n2*pow(p2,2)*z2 - 
   2*g0*k2*m2*p2*r0*z2 + 2*f2*k2*n0*p2*r0*z2 + 2*f0*k2*n2*p2*r0*z2 + 
   g2*pow(k2,2)*pow(r0,2)*z2 + 2*f0*k2*n0*p2*r2*z2 + 
   2*g0*pow(k2,2)*r0*r2*z2 + 2*f0*k2*m2*p2*s0*z2 - 2*e2*k2*n0*p2*s0*z2 - 
   2*e0*k2*n2*p2*s0*z2 - 2*f2*pow(k2,2)*r0*s0*z2 - 
   2*f0*pow(k2,2)*r2*s0*z2 + e2*pow(k2,2)*pow(s0,2)*z2 - 
   2*e0*k2*n0*p2*s2*z2 - 2*f0*pow(k2,2)*r0*s2*z2 + 
   2*e0*pow(k2,2)*s0*s2*z2 + pow(f0,2)*pow(j2,2)*t2*z2 - 
   e0*g0*pow(j2,2)*t2*z2 - pow(f0,2)*i2*o2*t2*z2 + e0*g0*i2*o2*t2*z2 + 
   2*f0*m2*n0*o2*t2*z2 - e2*pow(n0,2)*o2*t2*z2 - 2*e0*n0*n2*o2*t2*z2 + 
   2*g0*j2*m2*r0*t2*z2 - 2*f2*j2*n0*r0*t2*z2 - 2*f0*j2*n2*r0*t2*z2 - 
   g2*i2*pow(r0,2)*t2*z2 + pow(n2,2)*pow(r0,2)*t2*z2 - 
   2*f0*j2*n0*r2*t2*z2 - 2*g0*i2*r0*r2*t2*z2 + 4*n0*n2*r0*r2*t2*z2 + 
   pow(n0,2)*pow(r2,2)*t2*z2 - 2*f0*j2*m2*s0*t2*z2 + 
   2*e2*j2*n0*s0*t2*z2 + 2*e0*j2*n2*s0*t2*z2 + 2*f2*i2*r0*s0*t2*z2 - 
   2*m2*n2*r0*s0*t2*z2 + 2*f0*i2*r2*s0*t2*z2 - 2*m2*n0*r2*s0*t2*z2 - 
   e2*i2*pow(s0,2)*t2*z2 + pow(m2,2)*pow(s0,2)*t2*z2 + 
   2*e0*j2*n0*s2*t2*z2 + 2*f0*i2*r0*s2*t2*z2 - 2*m2*n0*r0*s2*t2*z2 - 
   2*e0*i2*s0*s2*t2*z2 + 2*g0*k2*m2*o2*v0*z2 - 2*f2*k2*n0*o2*v0*z2 - 
   2*f0*k2*n2*o2*v0*z2 - 2*g0*j2*m2*p2*v0*z2 + 2*f2*j2*n0*p2*v0*z2 + 
   2*f0*j2*n2*p2*v0*z2 - 2*g2*j2*k2*r0*v0*z2 + 2*g2*i2*p2*r0*v0*z2 - 
   2*pow(n2,2)*p2*r0*v0*z2 - 2*g0*j2*k2*r2*v0*z2 + 2*g0*i2*p2*r2*v0*z2 - 
   4*n0*n2*p2*r2*v0*z2 + 2*f2*j2*k2*s0*v0*z2 - 2*f2*i2*p2*s0*v0*z2 + 
   2*m2*n2*p2*s0*v0*z2 + 2*k2*n2*r2*s0*v0*z2 + 2*f0*j2*k2*s2*v0*z2 - 
   2*f0*i2*p2*s2*v0*z2 + 2*m2*n0*p2*s2*v0*z2 + 2*k2*n2*r0*s2*v0*z2 + 
   2*k2*n0*r2*s2*v0*z2 - 4*k2*m2*s0*s2*v0*z2 + 
   g2*pow(j2,2)*pow(v0,2)*z2 - g2*i2*o2*pow(v0,2)*z2 + 
   pow(n2,2)*o2*pow(v0,2)*z2 - 2*j2*n2*s2*pow(v0,2)*z2 + 
   i2*pow(s2,2)*pow(v0,2)*z2 - 2*f0*k2*n0*o2*v2*z2 + 
   2*f0*j2*n0*p2*v2*z2 - 2*g0*j2*k2*r0*v2*z2 + 2*g0*i2*p2*r0*v2*z2 - 
   4*n0*n2*p2*r0*v2*z2 - 2*pow(n0,2)*p2*r2*v2*z2 + 2*f0*j2*k2*s0*v2*z2 - 
   2*f0*i2*p2*s0*v2*z2 + 2*m2*n0*p2*s0*v2*z2 + 2*k2*n2*r0*s0*v2*z2 + 
   2*k2*n0*r2*s0*v2*z2 - 2*k2*m2*pow(s0,2)*v2*z2 + 2*k2*n0*r0*s2*v2*z2 + 
   2*g0*pow(j2,2)*v0*v2*z2 - 2*g0*i2*o2*v0*v2*z2 + 4*n0*n2*o2*v0*v2*z2 - 
   4*j2*n2*s0*v0*v2*z2 - 4*j2*n0*s2*v0*v2*z2 + 4*i2*s0*s2*v0*v2*z2 + 
   pow(n0,2)*o2*pow(v2,2)*z2 - 2*j2*n0*s0*pow(v2,2)*z2 + 
   i2*pow(s0,2)*pow(v2,2)*z2 - 2*f0*k2*m2*o2*w0*z2 + 
   2*e2*k2*n0*o2*w0*z2 + 2*e0*k2*n2*o2*w0*z2 + 2*f0*j2*m2*p2*w0*z2 - 
   2*e2*j2*n0*p2*w0*z2 - 2*e0*j2*n2*p2*w0*z2 + 2*f2*j2*k2*r0*w0*z2 - 
   2*f2*i2*p2*r0*w0*z2 + 2*m2*n2*p2*r0*w0*z2 + 2*f0*j2*k2*r2*w0*z2 - 
   2*f0*i2*p2*r2*w0*z2 + 2*m2*n0*p2*r2*w0*z2 - 4*k2*n2*r0*r2*w0*z2 - 
   2*k2*n0*pow(r2,2)*w0*z2 - 2*e2*j2*k2*s0*w0*z2 + 2*e2*i2*p2*s0*w0*z2 - 
   2*pow(m2,2)*p2*s0*w0*z2 + 2*k2*m2*r2*s0*w0*z2 - 2*e0*j2*k2*s2*w0*z2 + 
   2*e0*i2*p2*s2*w0*z2 + 2*k2*m2*r0*s2*w0*z2 - 2*f2*pow(j2,2)*v0*w0*z2 + 
   2*f2*i2*o2*v0*w0*z2 - 2*m2*n2*o2*v0*w0*z2 + 2*j2*n2*r2*v0*w0*z2 + 
   2*j2*m2*s2*v0*w0*z2 - 2*i2*r2*s2*v0*w0*z2 - 2*f0*pow(j2,2)*v2*w0*z2 + 
   2*f0*i2*o2*v2*w0*z2 - 2*m2*n0*o2*v2*w0*z2 + 2*j2*n2*r0*v2*w0*z2 + 
   2*j2*n0*r2*v2*w0*z2 + 2*j2*m2*s0*v2*w0*z2 - 2*i2*r2*s0*v2*w0*z2 - 
   2*i2*r0*s2*v2*w0*z2 + e2*pow(j2,2)*pow(w0,2)*z2 - 
   e2*i2*o2*pow(w0,2)*z2 + pow(m2,2)*o2*pow(w0,2)*z2 - 
   2*j2*m2*r2*pow(w0,2)*z2 + i2*pow(r2,2)*pow(w0,2)*z2 + 
   2*e0*k2*n0*o2*w2*z2 - 2*e0*j2*n0*p2*w2*z2 + 2*f0*j2*k2*r0*w2*z2 - 
   2*f0*i2*p2*r0*w2*z2 + 2*m2*n0*p2*r0*w2*z2 - 2*k2*n2*pow(r0,2)*w2*z2 - 
   4*k2*n0*r0*r2*w2*z2 - 2*e0*j2*k2*s0*w2*z2 + 2*e0*i2*p2*s0*w2*z2 + 
   2*k2*m2*r0*s0*w2*z2 - 2*f0*pow(j2,2)*v0*w2*z2 + 2*f0*i2*o2*v0*w2*z2 - 
   2*m2*n0*o2*v0*w2*z2 + 2*j2*n2*r0*v0*w2*z2 + 2*j2*n0*r2*v0*w2*z2 + 
   2*j2*m2*s0*v0*w2*z2 - 2*i2*r2*s0*v0*w2*z2 - 2*i2*r0*s2*v0*w2*z2 + 
   2*j2*n0*r0*v2*w2*z2 - 2*i2*r0*s0*v2*w2*z2 + 2*e0*pow(j2,2)*w0*w2*z2 - 
   2*e0*i2*o2*w0*w2*z2 - 4*j2*m2*r0*w0*w2*z2 + 4*i2*r0*r2*w0*w2*z2 + 
   i2*pow(r0,2)*pow(w2,2)*z2
;

	k[4]  = 2*d1*d2*e0*pow(k2,2)*o2 + 2*d0*d2*e1*pow(k2,2)*o2 + 
   2*d0*d1*e2*pow(k2,2)*o2 - 2*c2*d1*f0*pow(k2,2)*o2 - 
   2*c1*d2*f0*pow(k2,2)*o2 - 2*c2*d0*f1*pow(k2,2)*o2 - 
   2*c0*d2*f1*pow(k2,2)*o2 - 2*c1*d0*f2*pow(k2,2)*o2 - 
   2*c0*d1*f2*pow(k2,2)*o2 + 2*c1*c2*g0*pow(k2,2)*o2 + 
   2*c0*c2*g1*pow(k2,2)*o2 + 2*c0*c1*g2*pow(k2,2)*o2 - 
   4*d1*d2*e0*j2*k2*p2 - 4*d0*d2*e1*j2*k2*p2 - 4*d0*d1*e2*j2*k2*p2 + 
   4*c2*d1*f0*j2*k2*p2 + 4*c1*d2*f0*j2*k2*p2 + 4*c2*d0*f1*j2*k2*p2 + 
   4*c0*d2*f1*j2*k2*p2 + 4*c1*d0*f2*j2*k2*p2 + 4*c0*d1*f2*j2*k2*p2 - 
   4*c1*c2*g0*j2*k2*p2 - 4*c0*c2*g1*j2*k2*p2 - 4*c0*c1*g2*j2*k2*p2 + 
   2*d1*d2*e0*i2*pow(p2,2) + 2*d0*d2*e1*i2*pow(p2,2) + 
   2*d0*d1*e2*i2*pow(p2,2) - 2*c2*d1*f0*i2*pow(p2,2) - 
   2*c1*d2*f0*i2*pow(p2,2) - 2*c2*d0*f1*i2*pow(p2,2) - 
   2*c0*d2*f1*i2*pow(p2,2) - 2*c1*d0*f2*i2*pow(p2,2) - 
   2*c0*d1*f2*i2*pow(p2,2) + 2*c1*c2*g0*i2*pow(p2,2) + 
   2*c0*c2*g1*i2*pow(p2,2) + 2*c0*c1*g2*i2*pow(p2,2) - 
   2*f0*f1*pow(l2,2)*pow(p2,2) + e1*g0*pow(l2,2)*pow(p2,2) + 
   e0*g1*pow(l2,2)*pow(p2,2) + 2*d1*f0*l2*m2*pow(p2,2) + 
   2*d0*f1*l2*m2*pow(p2,2) - 2*c1*g0*l2*m2*pow(p2,2) - 
   2*c0*g1*l2*m2*pow(p2,2) - 2*d0*d1*pow(m2,2)*pow(p2,2) - 
   2*d2*e1*l2*n0*pow(p2,2) - 2*d1*e2*l2*n0*pow(p2,2) + 
   2*c2*f1*l2*n0*pow(p2,2) + 2*c1*f2*l2*n0*pow(p2,2) + 
   2*c2*d1*m2*n0*pow(p2,2) + 2*c1*d2*m2*n0*pow(p2,2) - 
   2*d2*e0*l2*n1*pow(p2,2) - 2*d0*e2*l2*n1*pow(p2,2) + 
   2*c2*f0*l2*n1*pow(p2,2) + 2*c0*f2*l2*n1*pow(p2,2) + 
   2*c2*d0*m2*n1*pow(p2,2) + 2*c0*d2*m2*n1*pow(p2,2) - 
   2*pow(c2,2)*n0*n1*pow(p2,2) - 2*d1*e0*l2*n2*pow(p2,2) - 
   2*d0*e1*l2*n2*pow(p2,2) + 2*c1*f0*l2*n2*pow(p2,2) + 
   2*c0*f1*l2*n2*pow(p2,2) + 2*c1*d0*m2*n2*pow(p2,2) + 
   2*c0*d1*m2*n2*pow(p2,2) - 4*c1*c2*n0*n2*pow(p2,2) - 
   4*c0*c2*n1*n2*pow(p2,2) - 2*c0*c1*pow(n2,2)*pow(p2,2) + 
   4*f0*f1*k2*l2*p2*q2 - 2*e1*g0*k2*l2*p2*q2 - 2*e0*g1*k2*l2*p2*q2 - 
   2*d1*f0*k2*m2*p2*q2 - 2*d0*f1*k2*m2*p2*q2 + 2*c1*g0*k2*m2*p2*q2 + 
   2*c0*g1*k2*m2*p2*q2 + 2*d2*e1*k2*n0*p2*q2 + 2*d1*e2*k2*n0*p2*q2 - 
   2*c2*f1*k2*n0*p2*q2 - 2*c1*f2*k2*n0*p2*q2 + 2*d2*e0*k2*n1*p2*q2 + 
   2*d0*e2*k2*n1*p2*q2 - 2*c2*f0*k2*n1*p2*q2 - 2*c0*f2*k2*n1*p2*q2 + 
   2*d1*e0*k2*n2*p2*q2 + 2*d0*e1*k2*n2*p2*q2 - 2*c1*f0*k2*n2*p2*q2 - 
   2*c0*f1*k2*n2*p2*q2 - 2*f0*f1*pow(k2,2)*pow(q2,2) + 
   e1*g0*pow(k2,2)*pow(q2,2) + e0*g1*pow(k2,2)*pow(q2,2) - 
   2*d2*f1*k2*l2*p2*r0 - 2*d1*f2*k2*l2*p2*r0 + 2*c2*g1*k2*l2*p2*r0 + 
   2*c1*g2*k2*l2*p2*r0 + 4*d1*d2*k2*m2*p2*r0 - 2*c2*d2*k2*n1*p2*r0 - 
   2*c2*d1*k2*n2*p2*r0 - 2*c1*d2*k2*n2*p2*r0 + 2*d2*f1*pow(k2,2)*q2*r0 + 
   2*d1*f2*pow(k2,2)*q2*r0 - 2*c2*g1*pow(k2,2)*q2*r0 - 
   2*c1*g2*pow(k2,2)*q2*r0 - 2*d2*f0*k2*l2*p2*r1 - 2*d0*f2*k2*l2*p2*r1 + 
   2*c2*g0*k2*l2*p2*r1 + 2*c0*g2*k2*l2*p2*r1 + 4*d0*d2*k2*m2*p2*r1 - 
   2*c2*d2*k2*n0*p2*r1 - 2*c2*d0*k2*n2*p2*r1 - 2*c0*d2*k2*n2*p2*r1 + 
   2*d2*f0*pow(k2,2)*q2*r1 + 2*d0*f2*pow(k2,2)*q2*r1 - 
   2*c2*g0*pow(k2,2)*q2*r1 - 2*c0*g2*pow(k2,2)*q2*r1 - 
   2*pow(d2,2)*pow(k2,2)*r0*r1 - 2*d1*f0*k2*l2*p2*r2 - 
   2*d0*f1*k2*l2*p2*r2 + 2*c1*g0*k2*l2*p2*r2 + 2*c0*g1*k2*l2*p2*r2 + 
   4*d0*d1*k2*m2*p2*r2 - 2*c2*d1*k2*n0*p2*r2 - 2*c1*d2*k2*n0*p2*r2 - 
   2*c2*d0*k2*n1*p2*r2 - 2*c0*d2*k2*n1*p2*r2 - 2*c1*d0*k2*n2*p2*r2 - 
   2*c0*d1*k2*n2*p2*r2 + 2*d1*f0*pow(k2,2)*q2*r2 + 
   2*d0*f1*pow(k2,2)*q2*r2 - 2*c1*g0*pow(k2,2)*q2*r2 - 
   2*c0*g1*pow(k2,2)*q2*r2 - 4*d1*d2*pow(k2,2)*r0*r2 - 
   4*d0*d2*pow(k2,2)*r1*r2 - 2*d0*d1*pow(k2,2)*pow(r2,2) + 
   2*d2*e1*k2*l2*p2*s0 + 2*d1*e2*k2*l2*p2*s0 - 2*c2*f1*k2*l2*p2*s0 - 
   2*c1*f2*k2*l2*p2*s0 - 2*c2*d1*k2*m2*p2*s0 - 2*c1*d2*k2*m2*p2*s0 + 
   2*pow(c2,2)*k2*n1*p2*s0 + 4*c1*c2*k2*n2*p2*s0 - 
   2*d2*e1*pow(k2,2)*q2*s0 - 2*d1*e2*pow(k2,2)*q2*s0 + 
   2*c2*f1*pow(k2,2)*q2*s0 + 2*c1*f2*pow(k2,2)*q2*s0 + 
   2*c2*d2*pow(k2,2)*r1*s0 + 2*c2*d1*pow(k2,2)*r2*s0 + 
   2*c1*d2*pow(k2,2)*r2*s0 + 2*d2*e0*k2*l2*p2*s1 + 2*d0*e2*k2*l2*p2*s1 - 
   2*c2*f0*k2*l2*p2*s1 - 2*c0*f2*k2*l2*p2*s1 - 2*c2*d0*k2*m2*p2*s1 - 
   2*c0*d2*k2*m2*p2*s1 + 2*pow(c2,2)*k2*n0*p2*s1 + 4*c0*c2*k2*n2*p2*s1 - 
   2*d2*e0*pow(k2,2)*q2*s1 - 2*d0*e2*pow(k2,2)*q2*s1 + 
   2*c2*f0*pow(k2,2)*q2*s1 + 2*c0*f2*pow(k2,2)*q2*s1 + 
   2*c2*d2*pow(k2,2)*r0*s1 + 2*c2*d0*pow(k2,2)*r2*s1 + 
   2*c0*d2*pow(k2,2)*r2*s1 - 2*pow(c2,2)*pow(k2,2)*s0*s1 + 
   2*d1*e0*k2*l2*p2*s2 + 2*d0*e1*k2*l2*p2*s2 - 2*c1*f0*k2*l2*p2*s2 - 
   2*c0*f1*k2*l2*p2*s2 - 2*c1*d0*k2*m2*p2*s2 - 2*c0*d1*k2*m2*p2*s2 + 
   4*c1*c2*k2*n0*p2*s2 + 4*c0*c2*k2*n1*p2*s2 + 4*c0*c1*k2*n2*p2*s2 - 
   2*d1*e0*pow(k2,2)*q2*s2 - 2*d0*e1*pow(k2,2)*q2*s2 + 
   2*c1*f0*pow(k2,2)*q2*s2 + 2*c0*f1*pow(k2,2)*q2*s2 + 
   2*c2*d1*pow(k2,2)*r0*s2 + 2*c1*d2*pow(k2,2)*r0*s2 + 
   2*c2*d0*pow(k2,2)*r1*s2 + 2*c0*d2*pow(k2,2)*r1*s2 + 
   2*c1*d0*pow(k2,2)*r2*s2 + 2*c0*d1*pow(k2,2)*r2*s2 - 
   4*c1*c2*pow(k2,2)*s0*s2 - 4*c0*c2*pow(k2,2)*s1*s2 - 
   2*c0*c1*pow(k2,2)*pow(s2,2) + 2*d1*d2*e0*pow(j2,2)*t2 + 
   2*d0*d2*e1*pow(j2,2)*t2 + 2*d0*d1*e2*pow(j2,2)*t2 - 
   2*c2*d1*f0*pow(j2,2)*t2 - 2*c1*d2*f0*pow(j2,2)*t2 - 
   2*c2*d0*f1*pow(j2,2)*t2 - 2*c0*d2*f1*pow(j2,2)*t2 - 
   2*c1*d0*f2*pow(j2,2)*t2 - 2*c0*d1*f2*pow(j2,2)*t2 + 
   2*c1*c2*g0*pow(j2,2)*t2 + 2*c0*c2*g1*pow(j2,2)*t2 + 
   2*c0*c1*g2*pow(j2,2)*t2 - 2*d1*d2*e0*i2*o2*t2 - 2*d0*d2*e1*i2*o2*t2 - 
   2*d0*d1*e2*i2*o2*t2 + 2*c2*d1*f0*i2*o2*t2 + 2*c1*d2*f0*i2*o2*t2 + 
   2*c2*d0*f1*i2*o2*t2 + 2*c0*d2*f1*i2*o2*t2 + 2*c1*d0*f2*i2*o2*t2 + 
   2*c0*d1*f2*i2*o2*t2 - 2*c1*c2*g0*i2*o2*t2 - 2*c0*c2*g1*i2*o2*t2 - 
   2*c0*c1*g2*i2*o2*t2 + 2*f0*f1*pow(l2,2)*o2*t2 - 
   e1*g0*pow(l2,2)*o2*t2 - e0*g1*pow(l2,2)*o2*t2 - 2*d1*f0*l2*m2*o2*t2 - 
   2*d0*f1*l2*m2*o2*t2 + 2*c1*g0*l2*m2*o2*t2 + 2*c0*g1*l2*m2*o2*t2 + 
   2*d0*d1*pow(m2,2)*o2*t2 + 2*d2*e1*l2*n0*o2*t2 + 2*d1*e2*l2*n0*o2*t2 - 
   2*c2*f1*l2*n0*o2*t2 - 2*c1*f2*l2*n0*o2*t2 - 2*c2*d1*m2*n0*o2*t2 - 
   2*c1*d2*m2*n0*o2*t2 + 2*d2*e0*l2*n1*o2*t2 + 2*d0*e2*l2*n1*o2*t2 - 
   2*c2*f0*l2*n1*o2*t2 - 2*c0*f2*l2*n1*o2*t2 - 2*c2*d0*m2*n1*o2*t2 - 
   2*c0*d2*m2*n1*o2*t2 + 2*pow(c2,2)*n0*n1*o2*t2 + 2*d1*e0*l2*n2*o2*t2 + 
   2*d0*e1*l2*n2*o2*t2 - 2*c1*f0*l2*n2*o2*t2 - 2*c0*f1*l2*n2*o2*t2 - 
   2*c1*d0*m2*n2*o2*t2 - 2*c0*d1*m2*n2*o2*t2 + 4*c1*c2*n0*n2*o2*t2 + 
   4*c0*c2*n1*n2*o2*t2 + 2*c0*c1*pow(n2,2)*o2*t2 - 4*f0*f1*j2*l2*q2*t2 + 
   2*e1*g0*j2*l2*q2*t2 + 2*e0*g1*j2*l2*q2*t2 + 2*d1*f0*j2*m2*q2*t2 + 
   2*d0*f1*j2*m2*q2*t2 - 2*c1*g0*j2*m2*q2*t2 - 2*c0*g1*j2*m2*q2*t2 - 
   2*d2*e1*j2*n0*q2*t2 - 2*d1*e2*j2*n0*q2*t2 + 2*c2*f1*j2*n0*q2*t2 + 
   2*c1*f2*j2*n0*q2*t2 - 2*d2*e0*j2*n1*q2*t2 - 2*d0*e2*j2*n1*q2*t2 + 
   2*c2*f0*j2*n1*q2*t2 + 2*c0*f2*j2*n1*q2*t2 - 2*d1*e0*j2*n2*q2*t2 - 
   2*d0*e1*j2*n2*q2*t2 + 2*c1*f0*j2*n2*q2*t2 + 2*c0*f1*j2*n2*q2*t2 + 
   2*f0*f1*i2*pow(q2,2)*t2 - e1*g0*i2*pow(q2,2)*t2 - 
   e0*g1*i2*pow(q2,2)*t2 - 2*f1*m2*n0*pow(q2,2)*t2 - 
   2*f0*m2*n1*pow(q2,2)*t2 + 2*e2*n0*n1*pow(q2,2)*t2 + 
   2*e1*n0*n2*pow(q2,2)*t2 + 2*e0*n1*n2*pow(q2,2)*t2 + 
   2*d2*f1*j2*l2*r0*t2 + 2*d1*f2*j2*l2*r0*t2 - 2*c2*g1*j2*l2*r0*t2 - 
   2*c1*g2*j2*l2*r0*t2 - 4*d1*d2*j2*m2*r0*t2 + 2*c2*d2*j2*n1*r0*t2 + 
   2*c2*d1*j2*n2*r0*t2 + 2*c1*d2*j2*n2*r0*t2 - 2*d2*f1*i2*q2*r0*t2 - 
   2*d1*f2*i2*q2*r0*t2 + 2*c2*g1*i2*q2*r0*t2 + 2*c1*g2*i2*q2*r0*t2 - 
   2*g1*l2*m2*q2*r0*t2 + 2*f2*l2*n1*q2*r0*t2 + 2*d2*m2*n1*q2*r0*t2 + 
   2*f1*l2*n2*q2*r0*t2 + 2*d1*m2*n2*q2*r0*t2 - 4*c2*n1*n2*q2*r0*t2 - 
   2*c1*pow(n2,2)*q2*r0*t2 + 2*d2*f0*j2*l2*r1*t2 + 2*d0*f2*j2*l2*r1*t2 - 
   2*c2*g0*j2*l2*r1*t2 - 2*c0*g2*j2*l2*r1*t2 - 4*d0*d2*j2*m2*r1*t2 + 
   2*c2*d2*j2*n0*r1*t2 + 2*c2*d0*j2*n2*r1*t2 + 2*c0*d2*j2*n2*r1*t2 - 
   2*d2*f0*i2*q2*r1*t2 - 2*d0*f2*i2*q2*r1*t2 + 2*c2*g0*i2*q2*r1*t2 + 
   2*c0*g2*i2*q2*r1*t2 - 2*g0*l2*m2*q2*r1*t2 + 2*f2*l2*n0*q2*r1*t2 + 
   2*d2*m2*n0*q2*r1*t2 + 2*f0*l2*n2*q2*r1*t2 + 2*d0*m2*n2*q2*r1*t2 - 
   4*c2*n0*n2*q2*r1*t2 - 2*c0*pow(n2,2)*q2*r1*t2 + 
   2*pow(d2,2)*i2*r0*r1*t2 + 2*g2*pow(l2,2)*r0*r1*t2 - 
   4*d2*l2*n2*r0*r1*t2 + 2*d1*f0*j2*l2*r2*t2 + 2*d0*f1*j2*l2*r2*t2 - 
   2*c1*g0*j2*l2*r2*t2 - 2*c0*g1*j2*l2*r2*t2 - 4*d0*d1*j2*m2*r2*t2 + 
   2*c2*d1*j2*n0*r2*t2 + 2*c1*d2*j2*n0*r2*t2 + 2*c2*d0*j2*n1*r2*t2 + 
   2*c0*d2*j2*n1*r2*t2 + 2*c1*d0*j2*n2*r2*t2 + 2*c0*d1*j2*n2*r2*t2 - 
   2*d1*f0*i2*q2*r2*t2 - 2*d0*f1*i2*q2*r2*t2 + 2*c1*g0*i2*q2*r2*t2 + 
   2*c0*g1*i2*q2*r2*t2 + 2*f1*l2*n0*q2*r2*t2 + 2*d1*m2*n0*q2*r2*t2 + 
   2*f0*l2*n1*q2*r2*t2 + 2*d0*m2*n1*q2*r2*t2 - 4*c2*n0*n1*q2*r2*t2 - 
   4*c1*n0*n2*q2*r2*t2 - 4*c0*n1*n2*q2*r2*t2 + 4*d1*d2*i2*r0*r2*t2 + 
   2*g1*pow(l2,2)*r0*r2*t2 - 4*d2*l2*n1*r0*r2*t2 - 4*d1*l2*n2*r0*r2*t2 + 
   4*d0*d2*i2*r1*r2*t2 + 2*g0*pow(l2,2)*r1*r2*t2 - 4*d2*l2*n0*r1*r2*t2 - 
   4*d0*l2*n2*r1*r2*t2 + 2*d0*d1*i2*pow(r2,2)*t2 - 
   2*d1*l2*n0*pow(r2,2)*t2 - 2*d0*l2*n1*pow(r2,2)*t2 - 
   2*d2*e1*j2*l2*s0*t2 - 2*d1*e2*j2*l2*s0*t2 + 2*c2*f1*j2*l2*s0*t2 + 
   2*c1*f2*j2*l2*s0*t2 + 2*c2*d1*j2*m2*s0*t2 + 2*c1*d2*j2*m2*s0*t2 - 
   2*pow(c2,2)*j2*n1*s0*t2 - 4*c1*c2*j2*n2*s0*t2 + 2*d2*e1*i2*q2*s0*t2 + 
   2*d1*e2*i2*q2*s0*t2 - 2*c2*f1*i2*q2*s0*t2 - 2*c1*f2*i2*q2*s0*t2 + 
   2*f1*l2*m2*q2*s0*t2 - 2*d1*pow(m2,2)*q2*s0*t2 - 2*e2*l2*n1*q2*s0*t2 + 
   2*c2*m2*n1*q2*s0*t2 - 2*e1*l2*n2*q2*s0*t2 + 2*c1*m2*n2*q2*s0*t2 - 
   2*c2*d2*i2*r1*s0*t2 - 2*f2*pow(l2,2)*r1*s0*t2 + 2*d2*l2*m2*r1*s0*t2 + 
   2*c2*l2*n2*r1*s0*t2 - 2*c2*d1*i2*r2*s0*t2 - 2*c1*d2*i2*r2*s0*t2 - 
   2*f1*pow(l2,2)*r2*s0*t2 + 2*d1*l2*m2*r2*s0*t2 + 2*c2*l2*n1*r2*s0*t2 + 
   2*c1*l2*n2*r2*s0*t2 - 2*d2*e0*j2*l2*s1*t2 - 2*d0*e2*j2*l2*s1*t2 + 
   2*c2*f0*j2*l2*s1*t2 + 2*c0*f2*j2*l2*s1*t2 + 2*c2*d0*j2*m2*s1*t2 + 
   2*c0*d2*j2*m2*s1*t2 - 2*pow(c2,2)*j2*n0*s1*t2 - 4*c0*c2*j2*n2*s1*t2 + 
   2*d2*e0*i2*q2*s1*t2 + 2*d0*e2*i2*q2*s1*t2 - 2*c2*f0*i2*q2*s1*t2 - 
   2*c0*f2*i2*q2*s1*t2 + 2*f0*l2*m2*q2*s1*t2 - 2*d0*pow(m2,2)*q2*s1*t2 - 
   2*e2*l2*n0*q2*s1*t2 + 2*c2*m2*n0*q2*s1*t2 - 2*e0*l2*n2*q2*s1*t2 + 
   2*c0*m2*n2*q2*s1*t2 - 2*c2*d2*i2*r0*s1*t2 - 2*f2*pow(l2,2)*r0*s1*t2 + 
   2*d2*l2*m2*r0*s1*t2 + 2*c2*l2*n2*r0*s1*t2 - 2*c2*d0*i2*r2*s1*t2 - 
   2*c0*d2*i2*r2*s1*t2 - 2*f0*pow(l2,2)*r2*s1*t2 + 2*d0*l2*m2*r2*s1*t2 + 
   2*c2*l2*n0*r2*s1*t2 + 2*c0*l2*n2*r2*s1*t2 + 2*pow(c2,2)*i2*s0*s1*t2 + 
   2*e2*pow(l2,2)*s0*s1*t2 - 4*c2*l2*m2*s0*s1*t2 - 2*d1*e0*j2*l2*s2*t2 - 
   2*d0*e1*j2*l2*s2*t2 + 2*c1*f0*j2*l2*s2*t2 + 2*c0*f1*j2*l2*s2*t2 + 
   2*c1*d0*j2*m2*s2*t2 + 2*c0*d1*j2*m2*s2*t2 - 4*c1*c2*j2*n0*s2*t2 - 
   4*c0*c2*j2*n1*s2*t2 - 4*c0*c1*j2*n2*s2*t2 + 2*d1*e0*i2*q2*s2*t2 + 
   2*d0*e1*i2*q2*s2*t2 - 2*c1*f0*i2*q2*s2*t2 - 2*c0*f1*i2*q2*s2*t2 - 
   2*e1*l2*n0*q2*s2*t2 + 2*c1*m2*n0*q2*s2*t2 - 2*e0*l2*n1*q2*s2*t2 + 
   2*c0*m2*n1*q2*s2*t2 - 2*c2*d1*i2*r0*s2*t2 - 2*c1*d2*i2*r0*s2*t2 - 
   2*f1*pow(l2,2)*r0*s2*t2 + 2*d1*l2*m2*r0*s2*t2 + 2*c2*l2*n1*r0*s2*t2 + 
   2*c1*l2*n2*r0*s2*t2 - 2*c2*d0*i2*r1*s2*t2 - 2*c0*d2*i2*r1*s2*t2 - 
   2*f0*pow(l2,2)*r1*s2*t2 + 2*d0*l2*m2*r1*s2*t2 + 2*c2*l2*n0*r1*s2*t2 + 
   2*c0*l2*n2*r1*s2*t2 - 2*c1*d0*i2*r2*s2*t2 - 2*c0*d1*i2*r2*s2*t2 + 
   2*c1*l2*n0*r2*s2*t2 + 2*c0*l2*n1*r2*s2*t2 + 4*c1*c2*i2*s0*s2*t2 + 
   2*e1*pow(l2,2)*s0*s2*t2 - 4*c1*l2*m2*s0*s2*t2 + 4*c0*c2*i2*s1*s2*t2 + 
   2*e0*pow(l2,2)*s1*s2*t2 - 4*c0*l2*m2*s1*s2*t2 + 
   2*c0*c1*i2*pow(s2,2)*t2 - 4*f1*f2*k2*l2*o2*u0 + 2*e2*g1*k2*l2*o2*u0 + 
   2*e1*g2*k2*l2*o2*u0 + 2*d2*f1*k2*m2*o2*u0 + 2*d1*f2*k2*m2*o2*u0 - 
   2*c2*g1*k2*m2*o2*u0 - 2*c1*g2*k2*m2*o2*u0 - 2*d2*e2*k2*n1*o2*u0 + 
   2*c2*f2*k2*n1*o2*u0 - 2*d2*e1*k2*n2*o2*u0 - 2*d1*e2*k2*n2*o2*u0 + 
   2*c2*f1*k2*n2*o2*u0 + 2*c1*f2*k2*n2*o2*u0 + 4*f1*f2*j2*l2*p2*u0 - 
   2*e2*g1*j2*l2*p2*u0 - 2*e1*g2*j2*l2*p2*u0 - 2*d2*f1*j2*m2*p2*u0 - 
   2*d1*f2*j2*m2*p2*u0 + 2*c2*g1*j2*m2*p2*u0 + 2*c1*g2*j2*m2*p2*u0 + 
   2*d2*e2*j2*n1*p2*u0 - 2*c2*f2*j2*n1*p2*u0 + 2*d2*e1*j2*n2*p2*u0 + 
   2*d1*e2*j2*n2*p2*u0 - 2*c2*f1*j2*n2*p2*u0 - 2*c1*f2*j2*n2*p2*u0 + 
   4*f1*f2*j2*k2*q2*u0 - 2*e2*g1*j2*k2*q2*u0 - 2*e1*g2*j2*k2*q2*u0 - 
   4*f1*f2*i2*p2*q2*u0 + 2*e2*g1*i2*p2*q2*u0 + 2*e1*g2*i2*p2*q2*u0 - 
   2*g1*pow(m2,2)*p2*q2*u0 + 4*f2*m2*n1*p2*q2*u0 + 4*f1*m2*n2*p2*q2*u0 - 
   4*e2*n1*n2*p2*q2*u0 - 2*e1*pow(n2,2)*p2*q2*u0 - 2*d2*f2*j2*k2*r1*u0 + 
   2*c2*g2*j2*k2*r1*u0 + 2*d2*f2*i2*p2*r1*u0 - 2*c2*g2*i2*p2*r1*u0 + 
   2*g2*l2*m2*p2*r1*u0 - 2*f2*l2*n2*p2*r1*u0 - 2*d2*m2*n2*p2*r1*u0 + 
   2*c2*pow(n2,2)*p2*r1*u0 + 2*g2*k2*m2*q2*r1*u0 - 2*f2*k2*n2*q2*r1*u0 - 
   2*d2*f1*j2*k2*r2*u0 - 2*d1*f2*j2*k2*r2*u0 + 2*c2*g1*j2*k2*r2*u0 + 
   2*c1*g2*j2*k2*r2*u0 + 2*d2*f1*i2*p2*r2*u0 + 2*d1*f2*i2*p2*r2*u0 - 
   2*c2*g1*i2*p2*r2*u0 - 2*c1*g2*i2*p2*r2*u0 + 2*g1*l2*m2*p2*r2*u0 - 
   2*f2*l2*n1*p2*r2*u0 - 2*d2*m2*n1*p2*r2*u0 - 2*f1*l2*n2*p2*r2*u0 - 
   2*d1*m2*n2*p2*r2*u0 + 4*c2*n1*n2*p2*r2*u0 + 2*c1*pow(n2,2)*p2*r2*u0 + 
   2*g1*k2*m2*q2*r2*u0 - 2*f2*k2*n1*q2*r2*u0 - 2*f1*k2*n2*q2*r2*u0 - 
   4*g2*k2*l2*r1*r2*u0 + 4*d2*k2*n2*r1*r2*u0 - 2*g1*k2*l2*pow(r2,2)*u0 + 
   2*d2*k2*n1*pow(r2,2)*u0 + 2*d1*k2*n2*pow(r2,2)*u0 + 
   2*d2*e2*j2*k2*s1*u0 - 2*c2*f2*j2*k2*s1*u0 - 2*d2*e2*i2*p2*s1*u0 + 
   2*c2*f2*i2*p2*s1*u0 - 2*f2*l2*m2*p2*s1*u0 + 2*d2*pow(m2,2)*p2*s1*u0 + 
   2*e2*l2*n2*p2*s1*u0 - 2*c2*m2*n2*p2*s1*u0 - 2*f2*k2*m2*q2*s1*u0 + 
   2*e2*k2*n2*q2*s1*u0 + 4*f2*k2*l2*r2*s1*u0 - 2*d2*k2*m2*r2*s1*u0 - 
   2*c2*k2*n2*r2*s1*u0 + 2*d2*e1*j2*k2*s2*u0 + 2*d1*e2*j2*k2*s2*u0 - 
   2*c2*f1*j2*k2*s2*u0 - 2*c1*f2*j2*k2*s2*u0 - 2*d2*e1*i2*p2*s2*u0 - 
   2*d1*e2*i2*p2*s2*u0 + 2*c2*f1*i2*p2*s2*u0 + 2*c1*f2*i2*p2*s2*u0 - 
   2*f1*l2*m2*p2*s2*u0 + 2*d1*pow(m2,2)*p2*s2*u0 + 2*e2*l2*n1*p2*s2*u0 - 
   2*c2*m2*n1*p2*s2*u0 + 2*e1*l2*n2*p2*s2*u0 - 2*c1*m2*n2*p2*s2*u0 - 
   2*f1*k2*m2*q2*s2*u0 + 2*e2*k2*n1*q2*s2*u0 + 2*e1*k2*n2*q2*s2*u0 + 
   4*f2*k2*l2*r1*s2*u0 - 2*d2*k2*m2*r1*s2*u0 - 2*c2*k2*n2*r1*s2*u0 + 
   4*f1*k2*l2*r2*s2*u0 - 2*d1*k2*m2*r2*s2*u0 - 2*c2*k2*n1*r2*s2*u0 - 
   2*c1*k2*n2*r2*s2*u0 - 4*e2*k2*l2*s1*s2*u0 + 4*c2*k2*m2*s1*s2*u0 - 
   2*e1*k2*l2*pow(s2,2)*u0 + 2*c1*k2*m2*pow(s2,2)*u0 - 
   4*f0*f2*k2*l2*o2*u1 + 2*e2*g0*k2*l2*o2*u1 + 2*e0*g2*k2*l2*o2*u1 + 
   2*d2*f0*k2*m2*o2*u1 + 2*d0*f2*k2*m2*o2*u1 - 2*c2*g0*k2*m2*o2*u1 - 
   2*c0*g2*k2*m2*o2*u1 - 2*d2*e2*k2*n0*o2*u1 + 2*c2*f2*k2*n0*o2*u1 - 
   2*d2*e0*k2*n2*o2*u1 - 2*d0*e2*k2*n2*o2*u1 + 2*c2*f0*k2*n2*o2*u1 + 
   2*c0*f2*k2*n2*o2*u1 + 4*f0*f2*j2*l2*p2*u1 - 2*e2*g0*j2*l2*p2*u1 - 
   2*e0*g2*j2*l2*p2*u1 - 2*d2*f0*j2*m2*p2*u1 - 2*d0*f2*j2*m2*p2*u1 + 
   2*c2*g0*j2*m2*p2*u1 + 2*c0*g2*j2*m2*p2*u1 + 2*d2*e2*j2*n0*p2*u1 - 
   2*c2*f2*j2*n0*p2*u1 + 2*d2*e0*j2*n2*p2*u1 + 2*d0*e2*j2*n2*p2*u1 - 
   2*c2*f0*j2*n2*p2*u1 - 2*c0*f2*j2*n2*p2*u1 + 4*f0*f2*j2*k2*q2*u1 - 
   2*e2*g0*j2*k2*q2*u1 - 2*e0*g2*j2*k2*q2*u1 - 4*f0*f2*i2*p2*q2*u1 + 
   2*e2*g0*i2*p2*q2*u1 + 2*e0*g2*i2*p2*q2*u1 - 2*g0*pow(m2,2)*p2*q2*u1 + 
   4*f2*m2*n0*p2*q2*u1 + 4*f0*m2*n2*p2*q2*u1 - 4*e2*n0*n2*p2*q2*u1 - 
   2*e0*pow(n2,2)*p2*q2*u1 - 2*d2*f2*j2*k2*r0*u1 + 2*c2*g2*j2*k2*r0*u1 + 
   2*d2*f2*i2*p2*r0*u1 - 2*c2*g2*i2*p2*r0*u1 + 2*g2*l2*m2*p2*r0*u1 - 
   2*f2*l2*n2*p2*r0*u1 - 2*d2*m2*n2*p2*r0*u1 + 2*c2*pow(n2,2)*p2*r0*u1 + 
   2*g2*k2*m2*q2*r0*u1 - 2*f2*k2*n2*q2*r0*u1 - 2*d2*f0*j2*k2*r2*u1 - 
   2*d0*f2*j2*k2*r2*u1 + 2*c2*g0*j2*k2*r2*u1 + 2*c0*g2*j2*k2*r2*u1 + 
   2*d2*f0*i2*p2*r2*u1 + 2*d0*f2*i2*p2*r2*u1 - 2*c2*g0*i2*p2*r2*u1 - 
   2*c0*g2*i2*p2*r2*u1 + 2*g0*l2*m2*p2*r2*u1 - 2*f2*l2*n0*p2*r2*u1 - 
   2*d2*m2*n0*p2*r2*u1 - 2*f0*l2*n2*p2*r2*u1 - 2*d0*m2*n2*p2*r2*u1 + 
   4*c2*n0*n2*p2*r2*u1 + 2*c0*pow(n2,2)*p2*r2*u1 + 2*g0*k2*m2*q2*r2*u1 - 
   2*f2*k2*n0*q2*r2*u1 - 2*f0*k2*n2*q2*r2*u1 - 4*g2*k2*l2*r0*r2*u1 + 
   4*d2*k2*n2*r0*r2*u1 - 2*g0*k2*l2*pow(r2,2)*u1 + 
   2*d2*k2*n0*pow(r2,2)*u1 + 2*d0*k2*n2*pow(r2,2)*u1 + 
   2*d2*e2*j2*k2*s0*u1 - 2*c2*f2*j2*k2*s0*u1 - 2*d2*e2*i2*p2*s0*u1 + 
   2*c2*f2*i2*p2*s0*u1 - 2*f2*l2*m2*p2*s0*u1 + 2*d2*pow(m2,2)*p2*s0*u1 + 
   2*e2*l2*n2*p2*s0*u1 - 2*c2*m2*n2*p2*s0*u1 - 2*f2*k2*m2*q2*s0*u1 + 
   2*e2*k2*n2*q2*s0*u1 + 4*f2*k2*l2*r2*s0*u1 - 2*d2*k2*m2*r2*s0*u1 - 
   2*c2*k2*n2*r2*s0*u1 + 2*d2*e0*j2*k2*s2*u1 + 2*d0*e2*j2*k2*s2*u1 - 
   2*c2*f0*j2*k2*s2*u1 - 2*c0*f2*j2*k2*s2*u1 - 2*d2*e0*i2*p2*s2*u1 - 
   2*d0*e2*i2*p2*s2*u1 + 2*c2*f0*i2*p2*s2*u1 + 2*c0*f2*i2*p2*s2*u1 - 
   2*f0*l2*m2*p2*s2*u1 + 2*d0*pow(m2,2)*p2*s2*u1 + 2*e2*l2*n0*p2*s2*u1 - 
   2*c2*m2*n0*p2*s2*u1 + 2*e0*l2*n2*p2*s2*u1 - 2*c0*m2*n2*p2*s2*u1 - 
   2*f0*k2*m2*q2*s2*u1 + 2*e2*k2*n0*q2*s2*u1 + 2*e0*k2*n2*q2*s2*u1 + 
   4*f2*k2*l2*r0*s2*u1 - 2*d2*k2*m2*r0*s2*u1 - 2*c2*k2*n2*r0*s2*u1 + 
   4*f0*k2*l2*r2*s2*u1 - 2*d0*k2*m2*r2*s2*u1 - 2*c2*k2*n0*r2*s2*u1 - 
   2*c0*k2*n2*r2*s2*u1 - 4*e2*k2*l2*s0*s2*u1 + 4*c2*k2*m2*s0*s2*u1 - 
   2*e0*k2*l2*pow(s2,2)*u1 + 2*c0*k2*m2*pow(s2,2)*u1 - 
   2*pow(f2,2)*pow(j2,2)*u0*u1 + 2*e2*g2*pow(j2,2)*u0*u1 + 
   2*pow(f2,2)*i2*o2*u0*u1 - 2*e2*g2*i2*o2*u0*u1 + 
   2*g2*pow(m2,2)*o2*u0*u1 - 4*f2*m2*n2*o2*u0*u1 + 
   2*e2*pow(n2,2)*o2*u0*u1 - 4*g2*j2*m2*r2*u0*u1 + 4*f2*j2*n2*r2*u0*u1 + 
   2*g2*i2*pow(r2,2)*u0*u1 - 2*pow(n2,2)*pow(r2,2)*u0*u1 + 
   4*f2*j2*m2*s2*u0*u1 - 4*e2*j2*n2*s2*u0*u1 - 4*f2*i2*r2*s2*u0*u1 + 
   4*m2*n2*r2*s2*u0*u1 + 2*e2*i2*pow(s2,2)*u0*u1 - 
   2*pow(m2,2)*pow(s2,2)*u0*u1 - 4*f0*f1*k2*l2*o2*u2 + 
   2*e1*g0*k2*l2*o2*u2 + 2*e0*g1*k2*l2*o2*u2 + 2*d1*f0*k2*m2*o2*u2 + 
   2*d0*f1*k2*m2*o2*u2 - 2*c1*g0*k2*m2*o2*u2 - 2*c0*g1*k2*m2*o2*u2 - 
   2*d2*e1*k2*n0*o2*u2 - 2*d1*e2*k2*n0*o2*u2 + 2*c2*f1*k2*n0*o2*u2 + 
   2*c1*f2*k2*n0*o2*u2 - 2*d2*e0*k2*n1*o2*u2 - 2*d0*e2*k2*n1*o2*u2 + 
   2*c2*f0*k2*n1*o2*u2 + 2*c0*f2*k2*n1*o2*u2 - 2*d1*e0*k2*n2*o2*u2 - 
   2*d0*e1*k2*n2*o2*u2 + 2*c1*f0*k2*n2*o2*u2 + 2*c0*f1*k2*n2*o2*u2 + 
   4*f0*f1*j2*l2*p2*u2 - 2*e1*g0*j2*l2*p2*u2 - 2*e0*g1*j2*l2*p2*u2 - 
   2*d1*f0*j2*m2*p2*u2 - 2*d0*f1*j2*m2*p2*u2 + 2*c1*g0*j2*m2*p2*u2 + 
   2*c0*g1*j2*m2*p2*u2 + 2*d2*e1*j2*n0*p2*u2 + 2*d1*e2*j2*n0*p2*u2 - 
   2*c2*f1*j2*n0*p2*u2 - 2*c1*f2*j2*n0*p2*u2 + 2*d2*e0*j2*n1*p2*u2 + 
   2*d0*e2*j2*n1*p2*u2 - 2*c2*f0*j2*n1*p2*u2 - 2*c0*f2*j2*n1*p2*u2 + 
   2*d1*e0*j2*n2*p2*u2 + 2*d0*e1*j2*n2*p2*u2 - 2*c1*f0*j2*n2*p2*u2 - 
   2*c0*f1*j2*n2*p2*u2 + 4*f0*f1*j2*k2*q2*u2 - 2*e1*g0*j2*k2*q2*u2 - 
   2*e0*g1*j2*k2*q2*u2 - 4*f0*f1*i2*p2*q2*u2 + 2*e1*g0*i2*p2*q2*u2 + 
   2*e0*g1*i2*p2*q2*u2 + 4*f1*m2*n0*p2*q2*u2 + 4*f0*m2*n1*p2*q2*u2 - 
   4*e2*n0*n1*p2*q2*u2 - 4*e1*n0*n2*p2*q2*u2 - 4*e0*n1*n2*p2*q2*u2 - 
   2*d2*f1*j2*k2*r0*u2 - 2*d1*f2*j2*k2*r0*u2 + 2*c2*g1*j2*k2*r0*u2 + 
   2*c1*g2*j2*k2*r0*u2 + 2*d2*f1*i2*p2*r0*u2 + 2*d1*f2*i2*p2*r0*u2 - 
   2*c2*g1*i2*p2*r0*u2 - 2*c1*g2*i2*p2*r0*u2 + 2*g1*l2*m2*p2*r0*u2 - 
   2*f2*l2*n1*p2*r0*u2 - 2*d2*m2*n1*p2*r0*u2 - 2*f1*l2*n2*p2*r0*u2 - 
   2*d1*m2*n2*p2*r0*u2 + 4*c2*n1*n2*p2*r0*u2 + 2*c1*pow(n2,2)*p2*r0*u2 + 
   2*g1*k2*m2*q2*r0*u2 - 2*f2*k2*n1*q2*r0*u2 - 2*f1*k2*n2*q2*r0*u2 - 
   2*d2*f0*j2*k2*r1*u2 - 2*d0*f2*j2*k2*r1*u2 + 2*c2*g0*j2*k2*r1*u2 + 
   2*c0*g2*j2*k2*r1*u2 + 2*d2*f0*i2*p2*r1*u2 + 2*d0*f2*i2*p2*r1*u2 - 
   2*c2*g0*i2*p2*r1*u2 - 2*c0*g2*i2*p2*r1*u2 + 2*g0*l2*m2*p2*r1*u2 - 
   2*f2*l2*n0*p2*r1*u2 - 2*d2*m2*n0*p2*r1*u2 - 2*f0*l2*n2*p2*r1*u2 - 
   2*d0*m2*n2*p2*r1*u2 + 4*c2*n0*n2*p2*r1*u2 + 2*c0*pow(n2,2)*p2*r1*u2 + 
   2*g0*k2*m2*q2*r1*u2 - 2*f2*k2*n0*q2*r1*u2 - 2*f0*k2*n2*q2*r1*u2 - 
   4*g2*k2*l2*r0*r1*u2 + 4*d2*k2*n2*r0*r1*u2 - 2*d1*f0*j2*k2*r2*u2 - 
   2*d0*f1*j2*k2*r2*u2 + 2*c1*g0*j2*k2*r2*u2 + 2*c0*g1*j2*k2*r2*u2 + 
   2*d1*f0*i2*p2*r2*u2 + 2*d0*f1*i2*p2*r2*u2 - 2*c1*g0*i2*p2*r2*u2 - 
   2*c0*g1*i2*p2*r2*u2 - 2*f1*l2*n0*p2*r2*u2 - 2*d1*m2*n0*p2*r2*u2 - 
   2*f0*l2*n1*p2*r2*u2 - 2*d0*m2*n1*p2*r2*u2 + 4*c2*n0*n1*p2*r2*u2 + 
   4*c1*n0*n2*p2*r2*u2 + 4*c0*n1*n2*p2*r2*u2 - 2*f1*k2*n0*q2*r2*u2 - 
   2*f0*k2*n1*q2*r2*u2 - 4*g1*k2*l2*r0*r2*u2 + 4*d2*k2*n1*r0*r2*u2 + 
   4*d1*k2*n2*r0*r2*u2 - 4*g0*k2*l2*r1*r2*u2 + 4*d2*k2*n0*r1*r2*u2 + 
   4*d0*k2*n2*r1*r2*u2 + 2*d1*k2*n0*pow(r2,2)*u2 + 
   2*d0*k2*n1*pow(r2,2)*u2 + 2*d2*e1*j2*k2*s0*u2 + 2*d1*e2*j2*k2*s0*u2 - 
   2*c2*f1*j2*k2*s0*u2 - 2*c1*f2*j2*k2*s0*u2 - 2*d2*e1*i2*p2*s0*u2 - 
   2*d1*e2*i2*p2*s0*u2 + 2*c2*f1*i2*p2*s0*u2 + 2*c1*f2*i2*p2*s0*u2 - 
   2*f1*l2*m2*p2*s0*u2 + 2*d1*pow(m2,2)*p2*s0*u2 + 2*e2*l2*n1*p2*s0*u2 - 
   2*c2*m2*n1*p2*s0*u2 + 2*e1*l2*n2*p2*s0*u2 - 2*c1*m2*n2*p2*s0*u2 - 
   2*f1*k2*m2*q2*s0*u2 + 2*e2*k2*n1*q2*s0*u2 + 2*e1*k2*n2*q2*s0*u2 + 
   4*f2*k2*l2*r1*s0*u2 - 2*d2*k2*m2*r1*s0*u2 - 2*c2*k2*n2*r1*s0*u2 + 
   4*f1*k2*l2*r2*s0*u2 - 2*d1*k2*m2*r2*s0*u2 - 2*c2*k2*n1*r2*s0*u2 - 
   2*c1*k2*n2*r2*s0*u2 + 2*d2*e0*j2*k2*s1*u2 + 2*d0*e2*j2*k2*s1*u2 - 
   2*c2*f0*j2*k2*s1*u2 - 2*c0*f2*j2*k2*s1*u2 - 2*d2*e0*i2*p2*s1*u2 - 
   2*d0*e2*i2*p2*s1*u2 + 2*c2*f0*i2*p2*s1*u2 + 2*c0*f2*i2*p2*s1*u2 - 
   2*f0*l2*m2*p2*s1*u2 + 2*d0*pow(m2,2)*p2*s1*u2 + 2*e2*l2*n0*p2*s1*u2 - 
   2*c2*m2*n0*p2*s1*u2 + 2*e0*l2*n2*p2*s1*u2 - 2*c0*m2*n2*p2*s1*u2 - 
   2*f0*k2*m2*q2*s1*u2 + 2*e2*k2*n0*q2*s1*u2 + 2*e0*k2*n2*q2*s1*u2 + 
   4*f2*k2*l2*r0*s1*u2 - 2*d2*k2*m2*r0*s1*u2 - 2*c2*k2*n2*r0*s1*u2 + 
   4*f0*k2*l2*r2*s1*u2 - 2*d0*k2*m2*r2*s1*u2 - 2*c2*k2*n0*r2*s1*u2 - 
   2*c0*k2*n2*r2*s1*u2 - 4*e2*k2*l2*s0*s1*u2 + 4*c2*k2*m2*s0*s1*u2 + 
   2*d1*e0*j2*k2*s2*u2 + 2*d0*e1*j2*k2*s2*u2 - 2*c1*f0*j2*k2*s2*u2 - 
   2*c0*f1*j2*k2*s2*u2 - 2*d1*e0*i2*p2*s2*u2 - 2*d0*e1*i2*p2*s2*u2 + 
   2*c1*f0*i2*p2*s2*u2 + 2*c0*f1*i2*p2*s2*u2 + 2*e1*l2*n0*p2*s2*u2 - 
   2*c1*m2*n0*p2*s2*u2 + 2*e0*l2*n1*p2*s2*u2 - 2*c0*m2*n1*p2*s2*u2 + 
   2*e1*k2*n0*q2*s2*u2 + 2*e0*k2*n1*q2*s2*u2 + 4*f1*k2*l2*r0*s2*u2 - 
   2*d1*k2*m2*r0*s2*u2 - 2*c2*k2*n1*r0*s2*u2 - 2*c1*k2*n2*r0*s2*u2 + 
   4*f0*k2*l2*r1*s2*u2 - 2*d0*k2*m2*r1*s2*u2 - 2*c2*k2*n0*r1*s2*u2 - 
   2*c0*k2*n2*r1*s2*u2 - 2*c1*k2*n0*r2*s2*u2 - 2*c0*k2*n1*r2*s2*u2 - 
   4*e1*k2*l2*s0*s2*u2 + 4*c1*k2*m2*s0*s2*u2 - 4*e0*k2*l2*s1*s2*u2 + 
   4*c0*k2*m2*s1*s2*u2 - 4*f1*f2*pow(j2,2)*u0*u2 + 
   2*e2*g1*pow(j2,2)*u0*u2 + 2*e1*g2*pow(j2,2)*u0*u2 + 
   4*f1*f2*i2*o2*u0*u2 - 2*e2*g1*i2*o2*u0*u2 - 2*e1*g2*i2*o2*u0*u2 + 
   2*g1*pow(m2,2)*o2*u0*u2 - 4*f2*m2*n1*o2*u0*u2 - 4*f1*m2*n2*o2*u0*u2 + 
   4*e2*n1*n2*o2*u0*u2 + 2*e1*pow(n2,2)*o2*u0*u2 - 4*g2*j2*m2*r1*u0*u2 + 
   4*f2*j2*n2*r1*u0*u2 - 4*g1*j2*m2*r2*u0*u2 + 4*f2*j2*n1*r2*u0*u2 + 
   4*f1*j2*n2*r2*u0*u2 + 4*g2*i2*r1*r2*u0*u2 - 4*pow(n2,2)*r1*r2*u0*u2 + 
   2*g1*i2*pow(r2,2)*u0*u2 - 4*n1*n2*pow(r2,2)*u0*u2 + 
   4*f2*j2*m2*s1*u0*u2 - 4*e2*j2*n2*s1*u0*u2 - 4*f2*i2*r2*s1*u0*u2 + 
   4*m2*n2*r2*s1*u0*u2 + 4*f1*j2*m2*s2*u0*u2 - 4*e2*j2*n1*s2*u0*u2 - 
   4*e1*j2*n2*s2*u0*u2 - 4*f2*i2*r1*s2*u0*u2 + 4*m2*n2*r1*s2*u0*u2 - 
   4*f1*i2*r2*s2*u0*u2 + 4*m2*n1*r2*s2*u0*u2 + 4*e2*i2*s1*s2*u0*u2 - 
   4*pow(m2,2)*s1*s2*u0*u2 + 2*e1*i2*pow(s2,2)*u0*u2 - 
   4*f0*f2*pow(j2,2)*u1*u2 + 2*e2*g0*pow(j2,2)*u1*u2 + 
   2*e0*g2*pow(j2,2)*u1*u2 + 4*f0*f2*i2*o2*u1*u2 - 2*e2*g0*i2*o2*u1*u2 - 
   2*e0*g2*i2*o2*u1*u2 + 2*g0*pow(m2,2)*o2*u1*u2 - 4*f2*m2*n0*o2*u1*u2 - 
   4*f0*m2*n2*o2*u1*u2 + 4*e2*n0*n2*o2*u1*u2 + 2*e0*pow(n2,2)*o2*u1*u2 - 
   4*g2*j2*m2*r0*u1*u2 + 4*f2*j2*n2*r0*u1*u2 - 4*g0*j2*m2*r2*u1*u2 + 
   4*f2*j2*n0*r2*u1*u2 + 4*f0*j2*n2*r2*u1*u2 + 4*g2*i2*r0*r2*u1*u2 - 
   4*pow(n2,2)*r0*r2*u1*u2 + 2*g0*i2*pow(r2,2)*u1*u2 - 
   4*n0*n2*pow(r2,2)*u1*u2 + 4*f2*j2*m2*s0*u1*u2 - 4*e2*j2*n2*s0*u1*u2 - 
   4*f2*i2*r2*s0*u1*u2 + 4*m2*n2*r2*s0*u1*u2 + 4*f0*j2*m2*s2*u1*u2 - 
   4*e2*j2*n0*s2*u1*u2 - 4*e0*j2*n2*s2*u1*u2 - 4*f2*i2*r0*s2*u1*u2 + 
   4*m2*n2*r0*s2*u1*u2 - 4*f0*i2*r2*s2*u1*u2 + 4*m2*n0*r2*s2*u1*u2 + 
   4*e2*i2*s0*s2*u1*u2 - 4*pow(m2,2)*s0*s2*u1*u2 + 
   2*e0*i2*pow(s2,2)*u1*u2 - 2*f0*f1*pow(j2,2)*pow(u2,2) + 
   e1*g0*pow(j2,2)*pow(u2,2) + e0*g1*pow(j2,2)*pow(u2,2) + 
   2*f0*f1*i2*o2*pow(u2,2) - e1*g0*i2*o2*pow(u2,2) - 
   e0*g1*i2*o2*pow(u2,2) - 2*f1*m2*n0*o2*pow(u2,2) - 
   2*f0*m2*n1*o2*pow(u2,2) + 2*e2*n0*n1*o2*pow(u2,2) + 
   2*e1*n0*n2*o2*pow(u2,2) + 2*e0*n1*n2*o2*pow(u2,2) - 
   2*g1*j2*m2*r0*pow(u2,2) + 2*f2*j2*n1*r0*pow(u2,2) + 
   2*f1*j2*n2*r0*pow(u2,2) - 2*g0*j2*m2*r1*pow(u2,2) + 
   2*f2*j2*n0*r1*pow(u2,2) + 2*f0*j2*n2*r1*pow(u2,2) + 
   2*g2*i2*r0*r1*pow(u2,2) - 2*pow(n2,2)*r0*r1*pow(u2,2) + 
   2*f1*j2*n0*r2*pow(u2,2) + 2*f0*j2*n1*r2*pow(u2,2) + 
   2*g1*i2*r0*r2*pow(u2,2) - 4*n1*n2*r0*r2*pow(u2,2) + 
   2*g0*i2*r1*r2*pow(u2,2) - 4*n0*n2*r1*r2*pow(u2,2) - 
   2*n0*n1*pow(r2,2)*pow(u2,2) + 2*f1*j2*m2*s0*pow(u2,2) - 
   2*e2*j2*n1*s0*pow(u2,2) - 2*e1*j2*n2*s0*pow(u2,2) - 
   2*f2*i2*r1*s0*pow(u2,2) + 2*m2*n2*r1*s0*pow(u2,2) - 
   2*f1*i2*r2*s0*pow(u2,2) + 2*m2*n1*r2*s0*pow(u2,2) + 
   2*f0*j2*m2*s1*pow(u2,2) - 2*e2*j2*n0*s1*pow(u2,2) - 
   2*e0*j2*n2*s1*pow(u2,2) - 2*f2*i2*r0*s1*pow(u2,2) + 
   2*m2*n2*r0*s1*pow(u2,2) - 2*f0*i2*r2*s1*pow(u2,2) + 
   2*m2*n0*r2*s1*pow(u2,2) + 2*e2*i2*s0*s1*pow(u2,2) - 
   2*pow(m2,2)*s0*s1*pow(u2,2) - 2*e1*j2*n0*s2*pow(u2,2) - 
   2*e0*j2*n1*s2*pow(u2,2) - 2*f1*i2*r0*s2*pow(u2,2) + 
   2*m2*n1*r0*s2*pow(u2,2) - 2*f0*i2*r1*s2*pow(u2,2) + 
   2*m2*n0*r1*s2*pow(u2,2) + 2*e1*i2*s0*s2*pow(u2,2) + 
   2*e0*i2*s1*s2*pow(u2,2) + 2*d2*f1*k2*l2*o2*v0 + 2*d1*f2*k2*l2*o2*v0 - 
   2*c2*g1*k2*l2*o2*v0 - 2*c1*g2*k2*l2*o2*v0 - 4*d1*d2*k2*m2*o2*v0 + 
   2*c2*d2*k2*n1*o2*v0 + 2*c2*d1*k2*n2*o2*v0 + 2*c1*d2*k2*n2*o2*v0 - 
   2*d2*f1*j2*l2*p2*v0 - 2*d1*f2*j2*l2*p2*v0 + 2*c2*g1*j2*l2*p2*v0 + 
   2*c1*g2*j2*l2*p2*v0 + 4*d1*d2*j2*m2*p2*v0 - 2*c2*d2*j2*n1*p2*v0 - 
   2*c2*d1*j2*n2*p2*v0 - 2*c1*d2*j2*n2*p2*v0 - 2*d2*f1*j2*k2*q2*v0 - 
   2*d1*f2*j2*k2*q2*v0 + 2*c2*g1*j2*k2*q2*v0 + 2*c1*g2*j2*k2*q2*v0 + 
   2*d2*f1*i2*p2*q2*v0 + 2*d1*f2*i2*p2*q2*v0 - 2*c2*g1*i2*p2*q2*v0 - 
   2*c1*g2*i2*p2*q2*v0 + 2*g1*l2*m2*p2*q2*v0 - 2*f2*l2*n1*p2*q2*v0 - 
   2*d2*m2*n1*p2*q2*v0 - 2*f1*l2*n2*p2*q2*v0 - 2*d1*m2*n2*p2*q2*v0 + 
   4*c2*n1*n2*p2*q2*v0 + 2*c1*pow(n2,2)*p2*q2*v0 - 
   2*g1*k2*m2*pow(q2,2)*v0 + 2*f2*k2*n1*pow(q2,2)*v0 + 
   2*f1*k2*n2*pow(q2,2)*v0 + 2*pow(d2,2)*j2*k2*r1*v0 - 
   2*pow(d2,2)*i2*p2*r1*v0 - 2*g2*pow(l2,2)*p2*r1*v0 + 
   4*d2*l2*n2*p2*r1*v0 + 2*g2*k2*l2*q2*r1*v0 - 2*d2*k2*n2*q2*r1*v0 + 
   4*d1*d2*j2*k2*r2*v0 - 4*d1*d2*i2*p2*r2*v0 - 2*g1*pow(l2,2)*p2*r2*v0 + 
   4*d2*l2*n1*p2*r2*v0 + 4*d1*l2*n2*p2*r2*v0 + 2*g1*k2*l2*q2*r2*v0 - 
   2*d2*k2*n1*q2*r2*v0 - 2*d1*k2*n2*q2*r2*v0 - 2*c2*d2*j2*k2*s1*v0 + 
   2*c2*d2*i2*p2*s1*v0 + 2*f2*pow(l2,2)*p2*s1*v0 - 2*d2*l2*m2*p2*s1*v0 - 
   2*c2*l2*n2*p2*s1*v0 - 2*f2*k2*l2*q2*s1*v0 + 4*d2*k2*m2*q2*s1*v0 - 
   2*c2*k2*n2*q2*s1*v0 - 2*d2*k2*l2*r2*s1*v0 - 2*c2*d1*j2*k2*s2*v0 - 
   2*c1*d2*j2*k2*s2*v0 + 2*c2*d1*i2*p2*s2*v0 + 2*c1*d2*i2*p2*s2*v0 + 
   2*f1*pow(l2,2)*p2*s2*v0 - 2*d1*l2*m2*p2*s2*v0 - 2*c2*l2*n1*p2*s2*v0 - 
   2*c1*l2*n2*p2*s2*v0 - 2*f1*k2*l2*q2*s2*v0 + 4*d1*k2*m2*q2*s2*v0 - 
   2*c2*k2*n1*q2*s2*v0 - 2*c1*k2*n2*q2*s2*v0 - 2*d2*k2*l2*r1*s2*v0 - 
   2*d1*k2*l2*r2*s2*v0 + 4*c2*k2*l2*s1*s2*v0 + 2*c1*k2*l2*pow(s2,2)*v0 + 
   2*d2*f2*pow(j2,2)*u1*v0 - 2*c2*g2*pow(j2,2)*u1*v0 - 
   2*d2*f2*i2*o2*u1*v0 + 2*c2*g2*i2*o2*u1*v0 - 2*g2*l2*m2*o2*u1*v0 + 
   2*f2*l2*n2*o2*u1*v0 + 2*d2*m2*n2*o2*u1*v0 - 2*c2*pow(n2,2)*o2*u1*v0 + 
   2*g2*j2*m2*q2*u1*v0 - 2*f2*j2*n2*q2*u1*v0 + 2*g2*j2*l2*r2*u1*v0 - 
   2*d2*j2*n2*r2*u1*v0 - 2*g2*i2*q2*r2*u1*v0 + 2*pow(n2,2)*q2*r2*u1*v0 - 
   2*f2*j2*l2*s2*u1*v0 - 2*d2*j2*m2*s2*u1*v0 + 4*c2*j2*n2*s2*u1*v0 + 
   2*f2*i2*q2*s2*u1*v0 - 2*m2*n2*q2*s2*u1*v0 + 2*d2*i2*r2*s2*u1*v0 - 
   2*l2*n2*r2*s2*u1*v0 - 2*c2*i2*pow(s2,2)*u1*v0 + 
   2*l2*m2*pow(s2,2)*u1*v0 + 2*d2*f1*pow(j2,2)*u2*v0 + 
   2*d1*f2*pow(j2,2)*u2*v0 - 2*c2*g1*pow(j2,2)*u2*v0 - 
   2*c1*g2*pow(j2,2)*u2*v0 - 2*d2*f1*i2*o2*u2*v0 - 2*d1*f2*i2*o2*u2*v0 + 
   2*c2*g1*i2*o2*u2*v0 + 2*c1*g2*i2*o2*u2*v0 - 2*g1*l2*m2*o2*u2*v0 + 
   2*f2*l2*n1*o2*u2*v0 + 2*d2*m2*n1*o2*u2*v0 + 2*f1*l2*n2*o2*u2*v0 + 
   2*d1*m2*n2*o2*u2*v0 - 4*c2*n1*n2*o2*u2*v0 - 2*c1*pow(n2,2)*o2*u2*v0 + 
   2*g1*j2*m2*q2*u2*v0 - 2*f2*j2*n1*q2*u2*v0 - 2*f1*j2*n2*q2*u2*v0 + 
   2*g2*j2*l2*r1*u2*v0 - 2*d2*j2*n2*r1*u2*v0 - 2*g2*i2*q2*r1*u2*v0 + 
   2*pow(n2,2)*q2*r1*u2*v0 + 2*g1*j2*l2*r2*u2*v0 - 2*d2*j2*n1*r2*u2*v0 - 
   2*d1*j2*n2*r2*u2*v0 - 2*g1*i2*q2*r2*u2*v0 + 4*n1*n2*q2*r2*u2*v0 - 
   2*f2*j2*l2*s1*u2*v0 - 2*d2*j2*m2*s1*u2*v0 + 4*c2*j2*n2*s1*u2*v0 + 
   2*f2*i2*q2*s1*u2*v0 - 2*m2*n2*q2*s1*u2*v0 + 2*d2*i2*r2*s1*u2*v0 - 
   2*l2*n2*r2*s1*u2*v0 - 2*f1*j2*l2*s2*u2*v0 - 2*d1*j2*m2*s2*u2*v0 + 
   4*c2*j2*n1*s2*u2*v0 + 4*c1*j2*n2*s2*u2*v0 + 2*f1*i2*q2*s2*u2*v0 - 
   2*m2*n1*q2*s2*u2*v0 + 2*d2*i2*r1*s2*u2*v0 - 2*l2*n2*r1*s2*u2*v0 + 
   2*d1*i2*r2*s2*u2*v0 - 2*l2*n1*r2*s2*u2*v0 - 4*c2*i2*s1*s2*u2*v0 + 
   4*l2*m2*s1*s2*u2*v0 - 2*c1*i2*pow(s2,2)*u2*v0 + 2*d2*f0*k2*l2*o2*v1 + 
   2*d0*f2*k2*l2*o2*v1 - 2*c2*g0*k2*l2*o2*v1 - 2*c0*g2*k2*l2*o2*v1 - 
   4*d0*d2*k2*m2*o2*v1 + 2*c2*d2*k2*n0*o2*v1 + 2*c2*d0*k2*n2*o2*v1 + 
   2*c0*d2*k2*n2*o2*v1 - 2*d2*f0*j2*l2*p2*v1 - 2*d0*f2*j2*l2*p2*v1 + 
   2*c2*g0*j2*l2*p2*v1 + 2*c0*g2*j2*l2*p2*v1 + 4*d0*d2*j2*m2*p2*v1 - 
   2*c2*d2*j2*n0*p2*v1 - 2*c2*d0*j2*n2*p2*v1 - 2*c0*d2*j2*n2*p2*v1 - 
   2*d2*f0*j2*k2*q2*v1 - 2*d0*f2*j2*k2*q2*v1 + 2*c2*g0*j2*k2*q2*v1 + 
   2*c0*g2*j2*k2*q2*v1 + 2*d2*f0*i2*p2*q2*v1 + 2*d0*f2*i2*p2*q2*v1 - 
   2*c2*g0*i2*p2*q2*v1 - 2*c0*g2*i2*p2*q2*v1 + 2*g0*l2*m2*p2*q2*v1 - 
   2*f2*l2*n0*p2*q2*v1 - 2*d2*m2*n0*p2*q2*v1 - 2*f0*l2*n2*p2*q2*v1 - 
   2*d0*m2*n2*p2*q2*v1 + 4*c2*n0*n2*p2*q2*v1 + 2*c0*pow(n2,2)*p2*q2*v1 - 
   2*g0*k2*m2*pow(q2,2)*v1 + 2*f2*k2*n0*pow(q2,2)*v1 + 
   2*f0*k2*n2*pow(q2,2)*v1 + 2*pow(d2,2)*j2*k2*r0*v1 - 
   2*pow(d2,2)*i2*p2*r0*v1 - 2*g2*pow(l2,2)*p2*r0*v1 + 
   4*d2*l2*n2*p2*r0*v1 + 2*g2*k2*l2*q2*r0*v1 - 2*d2*k2*n2*q2*r0*v1 + 
   4*d0*d2*j2*k2*r2*v1 - 4*d0*d2*i2*p2*r2*v1 - 2*g0*pow(l2,2)*p2*r2*v1 + 
   4*d2*l2*n0*p2*r2*v1 + 4*d0*l2*n2*p2*r2*v1 + 2*g0*k2*l2*q2*r2*v1 - 
   2*d2*k2*n0*q2*r2*v1 - 2*d0*k2*n2*q2*r2*v1 - 2*c2*d2*j2*k2*s0*v1 + 
   2*c2*d2*i2*p2*s0*v1 + 2*f2*pow(l2,2)*p2*s0*v1 - 2*d2*l2*m2*p2*s0*v1 - 
   2*c2*l2*n2*p2*s0*v1 - 2*f2*k2*l2*q2*s0*v1 + 4*d2*k2*m2*q2*s0*v1 - 
   2*c2*k2*n2*q2*s0*v1 - 2*d2*k2*l2*r2*s0*v1 - 2*c2*d0*j2*k2*s2*v1 - 
   2*c0*d2*j2*k2*s2*v1 + 2*c2*d0*i2*p2*s2*v1 + 2*c0*d2*i2*p2*s2*v1 + 
   2*f0*pow(l2,2)*p2*s2*v1 - 2*d0*l2*m2*p2*s2*v1 - 2*c2*l2*n0*p2*s2*v1 - 
   2*c0*l2*n2*p2*s2*v1 - 2*f0*k2*l2*q2*s2*v1 + 4*d0*k2*m2*q2*s2*v1 - 
   2*c2*k2*n0*q2*s2*v1 - 2*c0*k2*n2*q2*s2*v1 - 2*d2*k2*l2*r0*s2*v1 - 
   2*d0*k2*l2*r2*s2*v1 + 4*c2*k2*l2*s0*s2*v1 + 2*c0*k2*l2*pow(s2,2)*v1 + 
   2*d2*f2*pow(j2,2)*u0*v1 - 2*c2*g2*pow(j2,2)*u0*v1 - 
   2*d2*f2*i2*o2*u0*v1 + 2*c2*g2*i2*o2*u0*v1 - 2*g2*l2*m2*o2*u0*v1 + 
   2*f2*l2*n2*o2*u0*v1 + 2*d2*m2*n2*o2*u0*v1 - 2*c2*pow(n2,2)*o2*u0*v1 + 
   2*g2*j2*m2*q2*u0*v1 - 2*f2*j2*n2*q2*u0*v1 + 2*g2*j2*l2*r2*u0*v1 - 
   2*d2*j2*n2*r2*u0*v1 - 2*g2*i2*q2*r2*u0*v1 + 2*pow(n2,2)*q2*r2*u0*v1 - 
   2*f2*j2*l2*s2*u0*v1 - 2*d2*j2*m2*s2*u0*v1 + 4*c2*j2*n2*s2*u0*v1 + 
   2*f2*i2*q2*s2*u0*v1 - 2*m2*n2*q2*s2*u0*v1 + 2*d2*i2*r2*s2*u0*v1 - 
   2*l2*n2*r2*s2*u0*v1 - 2*c2*i2*pow(s2,2)*u0*v1 + 
   2*l2*m2*pow(s2,2)*u0*v1 + 2*d2*f0*pow(j2,2)*u2*v1 + 
   2*d0*f2*pow(j2,2)*u2*v1 - 2*c2*g0*pow(j2,2)*u2*v1 - 
   2*c0*g2*pow(j2,2)*u2*v1 - 2*d2*f0*i2*o2*u2*v1 - 2*d0*f2*i2*o2*u2*v1 + 
   2*c2*g0*i2*o2*u2*v1 + 2*c0*g2*i2*o2*u2*v1 - 2*g0*l2*m2*o2*u2*v1 + 
   2*f2*l2*n0*o2*u2*v1 + 2*d2*m2*n0*o2*u2*v1 + 2*f0*l2*n2*o2*u2*v1 + 
   2*d0*m2*n2*o2*u2*v1 - 4*c2*n0*n2*o2*u2*v1 - 2*c0*pow(n2,2)*o2*u2*v1 + 
   2*g0*j2*m2*q2*u2*v1 - 2*f2*j2*n0*q2*u2*v1 - 2*f0*j2*n2*q2*u2*v1 + 
   2*g2*j2*l2*r0*u2*v1 - 2*d2*j2*n2*r0*u2*v1 - 2*g2*i2*q2*r0*u2*v1 + 
   2*pow(n2,2)*q2*r0*u2*v1 + 2*g0*j2*l2*r2*u2*v1 - 2*d2*j2*n0*r2*u2*v1 - 
   2*d0*j2*n2*r2*u2*v1 - 2*g0*i2*q2*r2*u2*v1 + 4*n0*n2*q2*r2*u2*v1 - 
   2*f2*j2*l2*s0*u2*v1 - 2*d2*j2*m2*s0*u2*v1 + 4*c2*j2*n2*s0*u2*v1 + 
   2*f2*i2*q2*s0*u2*v1 - 2*m2*n2*q2*s0*u2*v1 + 2*d2*i2*r2*s0*u2*v1 - 
   2*l2*n2*r2*s0*u2*v1 - 2*f0*j2*l2*s2*u2*v1 - 2*d0*j2*m2*s2*u2*v1 + 
   4*c2*j2*n0*s2*u2*v1 + 4*c0*j2*n2*s2*u2*v1 + 2*f0*i2*q2*s2*u2*v1 - 
   2*m2*n0*q2*s2*u2*v1 + 2*d2*i2*r0*s2*u2*v1 - 2*l2*n2*r0*s2*u2*v1 + 
   2*d0*i2*r2*s2*u2*v1 - 2*l2*n0*r2*s2*u2*v1 - 4*c2*i2*s0*s2*u2*v1 + 
   4*l2*m2*s0*s2*u2*v1 - 2*c0*i2*pow(s2,2)*u2*v1 - 
   2*pow(d2,2)*pow(j2,2)*v0*v1 + 2*pow(d2,2)*i2*o2*v0*v1 + 
   2*g2*pow(l2,2)*o2*v0*v1 - 4*d2*l2*n2*o2*v0*v1 - 4*g2*j2*l2*q2*v0*v1 + 
   4*d2*j2*n2*q2*v0*v1 + 2*g2*i2*pow(q2,2)*v0*v1 - 
   2*pow(n2,2)*pow(q2,2)*v0*v1 + 4*d2*j2*l2*s2*v0*v1 - 
   4*d2*i2*q2*s2*v0*v1 + 4*l2*n2*q2*s2*v0*v1 - 
   2*pow(l2,2)*pow(s2,2)*v0*v1 + 2*d1*f0*k2*l2*o2*v2 + 
   2*d0*f1*k2*l2*o2*v2 - 2*c1*g0*k2*l2*o2*v2 - 2*c0*g1*k2*l2*o2*v2 - 
   4*d0*d1*k2*m2*o2*v2 + 2*c2*d1*k2*n0*o2*v2 + 2*c1*d2*k2*n0*o2*v2 + 
   2*c2*d0*k2*n1*o2*v2 + 2*c0*d2*k2*n1*o2*v2 + 2*c1*d0*k2*n2*o2*v2 + 
   2*c0*d1*k2*n2*o2*v2 - 2*d1*f0*j2*l2*p2*v2 - 2*d0*f1*j2*l2*p2*v2 + 
   2*c1*g0*j2*l2*p2*v2 + 2*c0*g1*j2*l2*p2*v2 + 4*d0*d1*j2*m2*p2*v2 - 
   2*c2*d1*j2*n0*p2*v2 - 2*c1*d2*j2*n0*p2*v2 - 2*c2*d0*j2*n1*p2*v2 - 
   2*c0*d2*j2*n1*p2*v2 - 2*c1*d0*j2*n2*p2*v2 - 2*c0*d1*j2*n2*p2*v2 - 
   2*d1*f0*j2*k2*q2*v2 - 2*d0*f1*j2*k2*q2*v2 + 2*c1*g0*j2*k2*q2*v2 + 
   2*c0*g1*j2*k2*q2*v2 + 2*d1*f0*i2*p2*q2*v2 + 2*d0*f1*i2*p2*q2*v2 - 
   2*c1*g0*i2*p2*q2*v2 - 2*c0*g1*i2*p2*q2*v2 - 2*f1*l2*n0*p2*q2*v2 - 
   2*d1*m2*n0*p2*q2*v2 - 2*f0*l2*n1*p2*q2*v2 - 2*d0*m2*n1*p2*q2*v2 + 
   4*c2*n0*n1*p2*q2*v2 + 4*c1*n0*n2*p2*q2*v2 + 4*c0*n1*n2*p2*q2*v2 + 
   2*f1*k2*n0*pow(q2,2)*v2 + 2*f0*k2*n1*pow(q2,2)*v2 + 
   4*d1*d2*j2*k2*r0*v2 - 4*d1*d2*i2*p2*r0*v2 - 2*g1*pow(l2,2)*p2*r0*v2 + 
   4*d2*l2*n1*p2*r0*v2 + 4*d1*l2*n2*p2*r0*v2 + 2*g1*k2*l2*q2*r0*v2 - 
   2*d2*k2*n1*q2*r0*v2 - 2*d1*k2*n2*q2*r0*v2 + 4*d0*d2*j2*k2*r1*v2 - 
   4*d0*d2*i2*p2*r1*v2 - 2*g0*pow(l2,2)*p2*r1*v2 + 4*d2*l2*n0*p2*r1*v2 + 
   4*d0*l2*n2*p2*r1*v2 + 2*g0*k2*l2*q2*r1*v2 - 2*d2*k2*n0*q2*r1*v2 - 
   2*d0*k2*n2*q2*r1*v2 + 4*d0*d1*j2*k2*r2*v2 - 4*d0*d1*i2*p2*r2*v2 + 
   4*d1*l2*n0*p2*r2*v2 + 4*d0*l2*n1*p2*r2*v2 - 2*d1*k2*n0*q2*r2*v2 - 
   2*d0*k2*n1*q2*r2*v2 - 2*c2*d1*j2*k2*s0*v2 - 2*c1*d2*j2*k2*s0*v2 + 
   2*c2*d1*i2*p2*s0*v2 + 2*c1*d2*i2*p2*s0*v2 + 2*f1*pow(l2,2)*p2*s0*v2 - 
   2*d1*l2*m2*p2*s0*v2 - 2*c2*l2*n1*p2*s0*v2 - 2*c1*l2*n2*p2*s0*v2 - 
   2*f1*k2*l2*q2*s0*v2 + 4*d1*k2*m2*q2*s0*v2 - 2*c2*k2*n1*q2*s0*v2 - 
   2*c1*k2*n2*q2*s0*v2 - 2*d2*k2*l2*r1*s0*v2 - 2*d1*k2*l2*r2*s0*v2 - 
   2*c2*d0*j2*k2*s1*v2 - 2*c0*d2*j2*k2*s1*v2 + 2*c2*d0*i2*p2*s1*v2 + 
   2*c0*d2*i2*p2*s1*v2 + 2*f0*pow(l2,2)*p2*s1*v2 - 2*d0*l2*m2*p2*s1*v2 - 
   2*c2*l2*n0*p2*s1*v2 - 2*c0*l2*n2*p2*s1*v2 - 2*f0*k2*l2*q2*s1*v2 + 
   4*d0*k2*m2*q2*s1*v2 - 2*c2*k2*n0*q2*s1*v2 - 2*c0*k2*n2*q2*s1*v2 - 
   2*d2*k2*l2*r0*s1*v2 - 2*d0*k2*l2*r2*s1*v2 + 4*c2*k2*l2*s0*s1*v2 - 
   2*c1*d0*j2*k2*s2*v2 - 2*c0*d1*j2*k2*s2*v2 + 2*c1*d0*i2*p2*s2*v2 + 
   2*c0*d1*i2*p2*s2*v2 - 2*c1*l2*n0*p2*s2*v2 - 2*c0*l2*n1*p2*s2*v2 - 
   2*c1*k2*n0*q2*s2*v2 - 2*c0*k2*n1*q2*s2*v2 - 2*d1*k2*l2*r0*s2*v2 - 
   2*d0*k2*l2*r1*s2*v2 + 4*c1*k2*l2*s0*s2*v2 + 4*c0*k2*l2*s1*s2*v2 + 
   2*d2*f1*pow(j2,2)*u0*v2 + 2*d1*f2*pow(j2,2)*u0*v2 - 
   2*c2*g1*pow(j2,2)*u0*v2 - 2*c1*g2*pow(j2,2)*u0*v2 - 
   2*d2*f1*i2*o2*u0*v2 - 2*d1*f2*i2*o2*u0*v2 + 2*c2*g1*i2*o2*u0*v2 + 
   2*c1*g2*i2*o2*u0*v2 - 2*g1*l2*m2*o2*u0*v2 + 2*f2*l2*n1*o2*u0*v2 + 
   2*d2*m2*n1*o2*u0*v2 + 2*f1*l2*n2*o2*u0*v2 + 2*d1*m2*n2*o2*u0*v2 - 
   4*c2*n1*n2*o2*u0*v2 - 2*c1*pow(n2,2)*o2*u0*v2 + 2*g1*j2*m2*q2*u0*v2 - 
   2*f2*j2*n1*q2*u0*v2 - 2*f1*j2*n2*q2*u0*v2 + 2*g2*j2*l2*r1*u0*v2 - 
   2*d2*j2*n2*r1*u0*v2 - 2*g2*i2*q2*r1*u0*v2 + 2*pow(n2,2)*q2*r1*u0*v2 + 
   2*g1*j2*l2*r2*u0*v2 - 2*d2*j2*n1*r2*u0*v2 - 2*d1*j2*n2*r2*u0*v2 - 
   2*g1*i2*q2*r2*u0*v2 + 4*n1*n2*q2*r2*u0*v2 - 2*f2*j2*l2*s1*u0*v2 - 
   2*d2*j2*m2*s1*u0*v2 + 4*c2*j2*n2*s1*u0*v2 + 2*f2*i2*q2*s1*u0*v2 - 
   2*m2*n2*q2*s1*u0*v2 + 2*d2*i2*r2*s1*u0*v2 - 2*l2*n2*r2*s1*u0*v2 - 
   2*f1*j2*l2*s2*u0*v2 - 2*d1*j2*m2*s2*u0*v2 + 4*c2*j2*n1*s2*u0*v2 + 
   4*c1*j2*n2*s2*u0*v2 + 2*f1*i2*q2*s2*u0*v2 - 2*m2*n1*q2*s2*u0*v2 + 
   2*d2*i2*r1*s2*u0*v2 - 2*l2*n2*r1*s2*u0*v2 + 2*d1*i2*r2*s2*u0*v2 - 
   2*l2*n1*r2*s2*u0*v2 - 4*c2*i2*s1*s2*u0*v2 + 4*l2*m2*s1*s2*u0*v2 - 
   2*c1*i2*pow(s2,2)*u0*v2 + 2*d2*f0*pow(j2,2)*u1*v2 + 
   2*d0*f2*pow(j2,2)*u1*v2 - 2*c2*g0*pow(j2,2)*u1*v2 - 
   2*c0*g2*pow(j2,2)*u1*v2 - 2*d2*f0*i2*o2*u1*v2 - 2*d0*f2*i2*o2*u1*v2 + 
   2*c2*g0*i2*o2*u1*v2 + 2*c0*g2*i2*o2*u1*v2 - 2*g0*l2*m2*o2*u1*v2 + 
   2*f2*l2*n0*o2*u1*v2 + 2*d2*m2*n0*o2*u1*v2 + 2*f0*l2*n2*o2*u1*v2 + 
   2*d0*m2*n2*o2*u1*v2 - 4*c2*n0*n2*o2*u1*v2 - 2*c0*pow(n2,2)*o2*u1*v2 + 
   2*g0*j2*m2*q2*u1*v2 - 2*f2*j2*n0*q2*u1*v2 - 2*f0*j2*n2*q2*u1*v2 + 
   2*g2*j2*l2*r0*u1*v2 - 2*d2*j2*n2*r0*u1*v2 - 2*g2*i2*q2*r0*u1*v2 + 
   2*pow(n2,2)*q2*r0*u1*v2 + 2*g0*j2*l2*r2*u1*v2 - 2*d2*j2*n0*r2*u1*v2 - 
   2*d0*j2*n2*r2*u1*v2 - 2*g0*i2*q2*r2*u1*v2 + 4*n0*n2*q2*r2*u1*v2 - 
   2*f2*j2*l2*s0*u1*v2 - 2*d2*j2*m2*s0*u1*v2 + 4*c2*j2*n2*s0*u1*v2 + 
   2*f2*i2*q2*s0*u1*v2 - 2*m2*n2*q2*s0*u1*v2 + 2*d2*i2*r2*s0*u1*v2 - 
   2*l2*n2*r2*s0*u1*v2 - 2*f0*j2*l2*s2*u1*v2 - 2*d0*j2*m2*s2*u1*v2 + 
   4*c2*j2*n0*s2*u1*v2 + 4*c0*j2*n2*s2*u1*v2 + 2*f0*i2*q2*s2*u1*v2 - 
   2*m2*n0*q2*s2*u1*v2 + 2*d2*i2*r0*s2*u1*v2 - 2*l2*n2*r0*s2*u1*v2 + 
   2*d0*i2*r2*s2*u1*v2 - 2*l2*n0*r2*s2*u1*v2 - 4*c2*i2*s0*s2*u1*v2 + 
   4*l2*m2*s0*s2*u1*v2 - 2*c0*i2*pow(s2,2)*u1*v2 + 
   2*d1*f0*pow(j2,2)*u2*v2 + 2*d0*f1*pow(j2,2)*u2*v2 - 
   2*c1*g0*pow(j2,2)*u2*v2 - 2*c0*g1*pow(j2,2)*u2*v2 - 
   2*d1*f0*i2*o2*u2*v2 - 2*d0*f1*i2*o2*u2*v2 + 2*c1*g0*i2*o2*u2*v2 + 
   2*c0*g1*i2*o2*u2*v2 + 2*f1*l2*n0*o2*u2*v2 + 2*d1*m2*n0*o2*u2*v2 + 
   2*f0*l2*n1*o2*u2*v2 + 2*d0*m2*n1*o2*u2*v2 - 4*c2*n0*n1*o2*u2*v2 - 
   4*c1*n0*n2*o2*u2*v2 - 4*c0*n1*n2*o2*u2*v2 - 2*f1*j2*n0*q2*u2*v2 - 
   2*f0*j2*n1*q2*u2*v2 + 2*g1*j2*l2*r0*u2*v2 - 2*d2*j2*n1*r0*u2*v2 - 
   2*d1*j2*n2*r0*u2*v2 - 2*g1*i2*q2*r0*u2*v2 + 4*n1*n2*q2*r0*u2*v2 + 
   2*g0*j2*l2*r1*u2*v2 - 2*d2*j2*n0*r1*u2*v2 - 2*d0*j2*n2*r1*u2*v2 - 
   2*g0*i2*q2*r1*u2*v2 + 4*n0*n2*q2*r1*u2*v2 - 2*d1*j2*n0*r2*u2*v2 - 
   2*d0*j2*n1*r2*u2*v2 + 4*n0*n1*q2*r2*u2*v2 - 2*f1*j2*l2*s0*u2*v2 - 
   2*d1*j2*m2*s0*u2*v2 + 4*c2*j2*n1*s0*u2*v2 + 4*c1*j2*n2*s0*u2*v2 + 
   2*f1*i2*q2*s0*u2*v2 - 2*m2*n1*q2*s0*u2*v2 + 2*d2*i2*r1*s0*u2*v2 - 
   2*l2*n2*r1*s0*u2*v2 + 2*d1*i2*r2*s0*u2*v2 - 2*l2*n1*r2*s0*u2*v2 - 
   2*f0*j2*l2*s1*u2*v2 - 2*d0*j2*m2*s1*u2*v2 + 4*c2*j2*n0*s1*u2*v2 + 
   4*c0*j2*n2*s1*u2*v2 + 2*f0*i2*q2*s1*u2*v2 - 2*m2*n0*q2*s1*u2*v2 + 
   2*d2*i2*r0*s1*u2*v2 - 2*l2*n2*r0*s1*u2*v2 + 2*d0*i2*r2*s1*u2*v2 - 
   2*l2*n0*r2*s1*u2*v2 - 4*c2*i2*s0*s1*u2*v2 + 4*l2*m2*s0*s1*u2*v2 + 
   4*c1*j2*n0*s2*u2*v2 + 4*c0*j2*n1*s2*u2*v2 + 2*d1*i2*r0*s2*u2*v2 - 
   2*l2*n1*r0*s2*u2*v2 + 2*d0*i2*r1*s2*u2*v2 - 2*l2*n0*r1*s2*u2*v2 - 
   4*c1*i2*s0*s2*u2*v2 - 4*c0*i2*s1*s2*u2*v2 - 4*d1*d2*pow(j2,2)*v0*v2 + 
   4*d1*d2*i2*o2*v0*v2 + 2*g1*pow(l2,2)*o2*v0*v2 - 4*d2*l2*n1*o2*v0*v2 - 
   4*d1*l2*n2*o2*v0*v2 - 4*g1*j2*l2*q2*v0*v2 + 4*d2*j2*n1*q2*v0*v2 + 
   4*d1*j2*n2*q2*v0*v2 + 2*g1*i2*pow(q2,2)*v0*v2 - 
   4*n1*n2*pow(q2,2)*v0*v2 + 4*d2*j2*l2*s1*v0*v2 - 4*d2*i2*q2*s1*v0*v2 + 
   4*l2*n2*q2*s1*v0*v2 + 4*d1*j2*l2*s2*v0*v2 - 4*d1*i2*q2*s2*v0*v2 + 
   4*l2*n1*q2*s2*v0*v2 - 4*pow(l2,2)*s1*s2*v0*v2 - 
   4*d0*d2*pow(j2,2)*v1*v2 + 4*d0*d2*i2*o2*v1*v2 + 
   2*g0*pow(l2,2)*o2*v1*v2 - 4*d2*l2*n0*o2*v1*v2 - 4*d0*l2*n2*o2*v1*v2 - 
   4*g0*j2*l2*q2*v1*v2 + 4*d2*j2*n0*q2*v1*v2 + 4*d0*j2*n2*q2*v1*v2 + 
   2*g0*i2*pow(q2,2)*v1*v2 - 4*n0*n2*pow(q2,2)*v1*v2 + 
   4*d2*j2*l2*s0*v1*v2 - 4*d2*i2*q2*s0*v1*v2 + 4*l2*n2*q2*s0*v1*v2 + 
   4*d0*j2*l2*s2*v1*v2 - 4*d0*i2*q2*s2*v1*v2 + 4*l2*n0*q2*s2*v1*v2 - 
   4*pow(l2,2)*s0*s2*v1*v2 - 2*d0*d1*pow(j2,2)*pow(v2,2) + 
   2*d0*d1*i2*o2*pow(v2,2) - 2*d1*l2*n0*o2*pow(v2,2) - 
   2*d0*l2*n1*o2*pow(v2,2) + 2*d1*j2*n0*q2*pow(v2,2) + 
   2*d0*j2*n1*q2*pow(v2,2) - 2*n0*n1*pow(q2,2)*pow(v2,2) + 
   2*d1*j2*l2*s0*pow(v2,2) - 2*d1*i2*q2*s0*pow(v2,2) + 
   2*l2*n1*q2*s0*pow(v2,2) + 2*d0*j2*l2*s1*pow(v2,2) - 
   2*d0*i2*q2*s1*pow(v2,2) + 2*l2*n0*q2*s1*pow(v2,2) - 
   2*pow(l2,2)*s0*s1*pow(v2,2) - 2*d2*e1*k2*l2*o2*w0 - 
   2*d1*e2*k2*l2*o2*w0 + 2*c2*f1*k2*l2*o2*w0 + 2*c1*f2*k2*l2*o2*w0 + 
   2*c2*d1*k2*m2*o2*w0 + 2*c1*d2*k2*m2*o2*w0 - 2*pow(c2,2)*k2*n1*o2*w0 - 
   4*c1*c2*k2*n2*o2*w0 + 2*d2*e1*j2*l2*p2*w0 + 2*d1*e2*j2*l2*p2*w0 - 
   2*c2*f1*j2*l2*p2*w0 - 2*c1*f2*j2*l2*p2*w0 - 2*c2*d1*j2*m2*p2*w0 - 
   2*c1*d2*j2*m2*p2*w0 + 2*pow(c2,2)*j2*n1*p2*w0 + 4*c1*c2*j2*n2*p2*w0 + 
   2*d2*e1*j2*k2*q2*w0 + 2*d1*e2*j2*k2*q2*w0 - 2*c2*f1*j2*k2*q2*w0 - 
   2*c1*f2*j2*k2*q2*w0 - 2*d2*e1*i2*p2*q2*w0 - 2*d1*e2*i2*p2*q2*w0 + 
   2*c2*f1*i2*p2*q2*w0 + 2*c1*f2*i2*p2*q2*w0 - 2*f1*l2*m2*p2*q2*w0 + 
   2*d1*pow(m2,2)*p2*q2*w0 + 2*e2*l2*n1*p2*q2*w0 - 2*c2*m2*n1*p2*q2*w0 + 
   2*e1*l2*n2*p2*q2*w0 - 2*c1*m2*n2*p2*q2*w0 + 2*f1*k2*m2*pow(q2,2)*w0 - 
   2*e2*k2*n1*pow(q2,2)*w0 - 2*e1*k2*n2*pow(q2,2)*w0 - 
   2*c2*d2*j2*k2*r1*w0 + 2*c2*d2*i2*p2*r1*w0 + 2*f2*pow(l2,2)*p2*r1*w0 - 
   2*d2*l2*m2*p2*r1*w0 - 2*c2*l2*n2*p2*r1*w0 - 2*f2*k2*l2*q2*r1*w0 - 
   2*d2*k2*m2*q2*r1*w0 + 4*c2*k2*n2*q2*r1*w0 - 2*c2*d1*j2*k2*r2*w0 - 
   2*c1*d2*j2*k2*r2*w0 + 2*c2*d1*i2*p2*r2*w0 + 2*c1*d2*i2*p2*r2*w0 + 
   2*f1*pow(l2,2)*p2*r2*w0 - 2*d1*l2*m2*p2*r2*w0 - 2*c2*l2*n1*p2*r2*w0 - 
   2*c1*l2*n2*p2*r2*w0 - 2*f1*k2*l2*q2*r2*w0 - 2*d1*k2*m2*q2*r2*w0 + 
   4*c2*k2*n1*q2*r2*w0 + 4*c1*k2*n2*q2*r2*w0 + 4*d2*k2*l2*r1*r2*w0 + 
   2*d1*k2*l2*pow(r2,2)*w0 + 2*pow(c2,2)*j2*k2*s1*w0 - 
   2*pow(c2,2)*i2*p2*s1*w0 - 2*e2*pow(l2,2)*p2*s1*w0 + 
   4*c2*l2*m2*p2*s1*w0 + 2*e2*k2*l2*q2*s1*w0 - 2*c2*k2*m2*q2*s1*w0 - 
   2*c2*k2*l2*r2*s1*w0 + 4*c1*c2*j2*k2*s2*w0 - 4*c1*c2*i2*p2*s2*w0 - 
   2*e1*pow(l2,2)*p2*s2*w0 + 4*c1*l2*m2*p2*s2*w0 + 2*e1*k2*l2*q2*s2*w0 - 
   2*c1*k2*m2*q2*s2*w0 - 2*c2*k2*l2*r1*s2*w0 - 2*c1*k2*l2*r2*s2*w0 - 
   2*d2*e2*pow(j2,2)*u1*w0 + 2*c2*f2*pow(j2,2)*u1*w0 + 
   2*d2*e2*i2*o2*u1*w0 - 2*c2*f2*i2*o2*u1*w0 + 2*f2*l2*m2*o2*u1*w0 - 
   2*d2*pow(m2,2)*o2*u1*w0 - 2*e2*l2*n2*o2*u1*w0 + 2*c2*m2*n2*o2*u1*w0 - 
   2*f2*j2*m2*q2*u1*w0 + 2*e2*j2*n2*q2*u1*w0 - 2*f2*j2*l2*r2*u1*w0 + 
   4*d2*j2*m2*r2*u1*w0 - 2*c2*j2*n2*r2*u1*w0 + 2*f2*i2*q2*r2*u1*w0 - 
   2*m2*n2*q2*r2*u1*w0 - 2*d2*i2*pow(r2,2)*u1*w0 + 
   2*l2*n2*pow(r2,2)*u1*w0 + 2*e2*j2*l2*s2*u1*w0 - 2*c2*j2*m2*s2*u1*w0 - 
   2*e2*i2*q2*s2*u1*w0 + 2*pow(m2,2)*q2*s2*u1*w0 + 2*c2*i2*r2*s2*u1*w0 - 
   2*l2*m2*r2*s2*u1*w0 - 2*d2*e1*pow(j2,2)*u2*w0 - 
   2*d1*e2*pow(j2,2)*u2*w0 + 2*c2*f1*pow(j2,2)*u2*w0 + 
   2*c1*f2*pow(j2,2)*u2*w0 + 2*d2*e1*i2*o2*u2*w0 + 2*d1*e2*i2*o2*u2*w0 - 
   2*c2*f1*i2*o2*u2*w0 - 2*c1*f2*i2*o2*u2*w0 + 2*f1*l2*m2*o2*u2*w0 - 
   2*d1*pow(m2,2)*o2*u2*w0 - 2*e2*l2*n1*o2*u2*w0 + 2*c2*m2*n1*o2*u2*w0 - 
   2*e1*l2*n2*o2*u2*w0 + 2*c1*m2*n2*o2*u2*w0 - 2*f1*j2*m2*q2*u2*w0 + 
   2*e2*j2*n1*q2*u2*w0 + 2*e1*j2*n2*q2*u2*w0 - 2*f2*j2*l2*r1*u2*w0 + 
   4*d2*j2*m2*r1*u2*w0 - 2*c2*j2*n2*r1*u2*w0 + 2*f2*i2*q2*r1*u2*w0 - 
   2*m2*n2*q2*r1*u2*w0 - 2*f1*j2*l2*r2*u2*w0 + 4*d1*j2*m2*r2*u2*w0 - 
   2*c2*j2*n1*r2*u2*w0 - 2*c1*j2*n2*r2*u2*w0 + 2*f1*i2*q2*r2*u2*w0 - 
   2*m2*n1*q2*r2*u2*w0 - 4*d2*i2*r1*r2*u2*w0 + 4*l2*n2*r1*r2*u2*w0 - 
   2*d1*i2*pow(r2,2)*u2*w0 + 2*l2*n1*pow(r2,2)*u2*w0 + 
   2*e2*j2*l2*s1*u2*w0 - 2*c2*j2*m2*s1*u2*w0 - 2*e2*i2*q2*s1*u2*w0 + 
   2*pow(m2,2)*q2*s1*u2*w0 + 2*c2*i2*r2*s1*u2*w0 - 2*l2*m2*r2*s1*u2*w0 + 
   2*e1*j2*l2*s2*u2*w0 - 2*c1*j2*m2*s2*u2*w0 - 2*e1*i2*q2*s2*u2*w0 + 
   2*c2*i2*r1*s2*u2*w0 - 2*l2*m2*r1*s2*u2*w0 + 2*c1*i2*r2*s2*u2*w0 + 
   2*c2*d2*pow(j2,2)*v1*w0 - 2*c2*d2*i2*o2*v1*w0 - 
   2*f2*pow(l2,2)*o2*v1*w0 + 2*d2*l2*m2*o2*v1*w0 + 2*c2*l2*n2*o2*v1*w0 + 
   4*f2*j2*l2*q2*v1*w0 - 2*d2*j2*m2*q2*v1*w0 - 2*c2*j2*n2*q2*v1*w0 - 
   2*f2*i2*pow(q2,2)*v1*w0 + 2*m2*n2*pow(q2,2)*v1*w0 - 
   2*d2*j2*l2*r2*v1*w0 + 2*d2*i2*q2*r2*v1*w0 - 2*l2*n2*q2*r2*v1*w0 - 
   2*c2*j2*l2*s2*v1*w0 + 2*c2*i2*q2*s2*v1*w0 - 2*l2*m2*q2*s2*v1*w0 + 
   2*pow(l2,2)*r2*s2*v1*w0 + 2*c2*d1*pow(j2,2)*v2*w0 + 
   2*c1*d2*pow(j2,2)*v2*w0 - 2*c2*d1*i2*o2*v2*w0 - 2*c1*d2*i2*o2*v2*w0 - 
   2*f1*pow(l2,2)*o2*v2*w0 + 2*d1*l2*m2*o2*v2*w0 + 2*c2*l2*n1*o2*v2*w0 + 
   2*c1*l2*n2*o2*v2*w0 + 4*f1*j2*l2*q2*v2*w0 - 2*d1*j2*m2*q2*v2*w0 - 
   2*c2*j2*n1*q2*v2*w0 - 2*c1*j2*n2*q2*v2*w0 - 2*f1*i2*pow(q2,2)*v2*w0 + 
   2*m2*n1*pow(q2,2)*v2*w0 - 2*d2*j2*l2*r1*v2*w0 + 2*d2*i2*q2*r1*v2*w0 - 
   2*l2*n2*q2*r1*v2*w0 - 2*d1*j2*l2*r2*v2*w0 + 2*d1*i2*q2*r2*v2*w0 - 
   2*l2*n1*q2*r2*v2*w0 - 2*c2*j2*l2*s1*v2*w0 + 2*c2*i2*q2*s1*v2*w0 - 
   2*l2*m2*q2*s1*v2*w0 + 2*pow(l2,2)*r2*s1*v2*w0 - 2*c1*j2*l2*s2*v2*w0 + 
   2*c1*i2*q2*s2*v2*w0 + 2*pow(l2,2)*r1*s2*v2*w0 - 2*d2*e0*k2*l2*o2*w1 - 
   2*d0*e2*k2*l2*o2*w1 + 2*c2*f0*k2*l2*o2*w1 + 2*c0*f2*k2*l2*o2*w1 + 
   2*c2*d0*k2*m2*o2*w1 + 2*c0*d2*k2*m2*o2*w1 - 2*pow(c2,2)*k2*n0*o2*w1 - 
   4*c0*c2*k2*n2*o2*w1 + 2*d2*e0*j2*l2*p2*w1 + 2*d0*e2*j2*l2*p2*w1 - 
   2*c2*f0*j2*l2*p2*w1 - 2*c0*f2*j2*l2*p2*w1 - 2*c2*d0*j2*m2*p2*w1 - 
   2*c0*d2*j2*m2*p2*w1 + 2*pow(c2,2)*j2*n0*p2*w1 + 4*c0*c2*j2*n2*p2*w1 + 
   2*d2*e0*j2*k2*q2*w1 + 2*d0*e2*j2*k2*q2*w1 - 2*c2*f0*j2*k2*q2*w1 - 
   2*c0*f2*j2*k2*q2*w1 - 2*d2*e0*i2*p2*q2*w1 - 2*d0*e2*i2*p2*q2*w1 + 
   2*c2*f0*i2*p2*q2*w1 + 2*c0*f2*i2*p2*q2*w1 - 2*f0*l2*m2*p2*q2*w1 + 
   2*d0*pow(m2,2)*p2*q2*w1 + 2*e2*l2*n0*p2*q2*w1 - 2*c2*m2*n0*p2*q2*w1 + 
   2*e0*l2*n2*p2*q2*w1 - 2*c0*m2*n2*p2*q2*w1 + 2*f0*k2*m2*pow(q2,2)*w1 - 
   2*e2*k2*n0*pow(q2,2)*w1 - 2*e0*k2*n2*pow(q2,2)*w1 - 
   2*c2*d2*j2*k2*r0*w1 + 2*c2*d2*i2*p2*r0*w1 + 2*f2*pow(l2,2)*p2*r0*w1 - 
   2*d2*l2*m2*p2*r0*w1 - 2*c2*l2*n2*p2*r0*w1 - 2*f2*k2*l2*q2*r0*w1 - 
   2*d2*k2*m2*q2*r0*w1 + 4*c2*k2*n2*q2*r0*w1 - 2*c2*d0*j2*k2*r2*w1 - 
   2*c0*d2*j2*k2*r2*w1 + 2*c2*d0*i2*p2*r2*w1 + 2*c0*d2*i2*p2*r2*w1 + 
   2*f0*pow(l2,2)*p2*r2*w1 - 2*d0*l2*m2*p2*r2*w1 - 2*c2*l2*n0*p2*r2*w1 - 
   2*c0*l2*n2*p2*r2*w1 - 2*f0*k2*l2*q2*r2*w1 - 2*d0*k2*m2*q2*r2*w1 + 
   4*c2*k2*n0*q2*r2*w1 + 4*c0*k2*n2*q2*r2*w1 + 4*d2*k2*l2*r0*r2*w1 + 
   2*d0*k2*l2*pow(r2,2)*w1 + 2*pow(c2,2)*j2*k2*s0*w1 - 
   2*pow(c2,2)*i2*p2*s0*w1 - 2*e2*pow(l2,2)*p2*s0*w1 + 
   4*c2*l2*m2*p2*s0*w1 + 2*e2*k2*l2*q2*s0*w1 - 2*c2*k2*m2*q2*s0*w1 - 
   2*c2*k2*l2*r2*s0*w1 + 4*c0*c2*j2*k2*s2*w1 - 4*c0*c2*i2*p2*s2*w1 - 
   2*e0*pow(l2,2)*p2*s2*w1 + 4*c0*l2*m2*p2*s2*w1 + 2*e0*k2*l2*q2*s2*w1 - 
   2*c0*k2*m2*q2*s2*w1 - 2*c2*k2*l2*r0*s2*w1 - 2*c0*k2*l2*r2*s2*w1 - 
   2*d2*e2*pow(j2,2)*u0*w1 + 2*c2*f2*pow(j2,2)*u0*w1 + 
   2*d2*e2*i2*o2*u0*w1 - 2*c2*f2*i2*o2*u0*w1 + 2*f2*l2*m2*o2*u0*w1 - 
   2*d2*pow(m2,2)*o2*u0*w1 - 2*e2*l2*n2*o2*u0*w1 + 2*c2*m2*n2*o2*u0*w1 - 
   2*f2*j2*m2*q2*u0*w1 + 2*e2*j2*n2*q2*u0*w1 - 2*f2*j2*l2*r2*u0*w1 + 
   4*d2*j2*m2*r2*u0*w1 - 2*c2*j2*n2*r2*u0*w1 + 2*f2*i2*q2*r2*u0*w1 - 
   2*m2*n2*q2*r2*u0*w1 - 2*d2*i2*pow(r2,2)*u0*w1 + 
   2*l2*n2*pow(r2,2)*u0*w1 + 2*e2*j2*l2*s2*u0*w1 - 2*c2*j2*m2*s2*u0*w1 - 
   2*e2*i2*q2*s2*u0*w1 + 2*pow(m2,2)*q2*s2*u0*w1 + 2*c2*i2*r2*s2*u0*w1 - 
   2*l2*m2*r2*s2*u0*w1 - 2*d2*e0*pow(j2,2)*u2*w1 - 
   2*d0*e2*pow(j2,2)*u2*w1 + 2*c2*f0*pow(j2,2)*u2*w1 + 
   2*c0*f2*pow(j2,2)*u2*w1 + 2*d2*e0*i2*o2*u2*w1 + 2*d0*e2*i2*o2*u2*w1 - 
   2*c2*f0*i2*o2*u2*w1 - 2*c0*f2*i2*o2*u2*w1 + 2*f0*l2*m2*o2*u2*w1 - 
   2*d0*pow(m2,2)*o2*u2*w1 - 2*e2*l2*n0*o2*u2*w1 + 2*c2*m2*n0*o2*u2*w1 - 
   2*e0*l2*n2*o2*u2*w1 + 2*c0*m2*n2*o2*u2*w1 - 2*f0*j2*m2*q2*u2*w1 + 
   2*e2*j2*n0*q2*u2*w1 + 2*e0*j2*n2*q2*u2*w1 - 2*f2*j2*l2*r0*u2*w1 + 
   4*d2*j2*m2*r0*u2*w1 - 2*c2*j2*n2*r0*u2*w1 + 2*f2*i2*q2*r0*u2*w1 - 
   2*m2*n2*q2*r0*u2*w1 - 2*f0*j2*l2*r2*u2*w1 + 4*d0*j2*m2*r2*u2*w1 - 
   2*c2*j2*n0*r2*u2*w1 - 2*c0*j2*n2*r2*u2*w1 + 2*f0*i2*q2*r2*u2*w1 - 
   2*m2*n0*q2*r2*u2*w1 - 4*d2*i2*r0*r2*u2*w1 + 4*l2*n2*r0*r2*u2*w1 - 
   2*d0*i2*pow(r2,2)*u2*w1 + 2*l2*n0*pow(r2,2)*u2*w1 + 
   2*e2*j2*l2*s0*u2*w1 - 2*c2*j2*m2*s0*u2*w1 - 2*e2*i2*q2*s0*u2*w1 + 
   2*pow(m2,2)*q2*s0*u2*w1 + 2*c2*i2*r2*s0*u2*w1 - 2*l2*m2*r2*s0*u2*w1 + 
   2*e0*j2*l2*s2*u2*w1 - 2*c0*j2*m2*s2*u2*w1 - 2*e0*i2*q2*s2*u2*w1 + 
   2*c2*i2*r0*s2*u2*w1 - 2*l2*m2*r0*s2*u2*w1 + 2*c0*i2*r2*s2*u2*w1 + 
   2*c2*d2*pow(j2,2)*v0*w1 - 2*c2*d2*i2*o2*v0*w1 - 
   2*f2*pow(l2,2)*o2*v0*w1 + 2*d2*l2*m2*o2*v0*w1 + 2*c2*l2*n2*o2*v0*w1 + 
   4*f2*j2*l2*q2*v0*w1 - 2*d2*j2*m2*q2*v0*w1 - 2*c2*j2*n2*q2*v0*w1 - 
   2*f2*i2*pow(q2,2)*v0*w1 + 2*m2*n2*pow(q2,2)*v0*w1 - 
   2*d2*j2*l2*r2*v0*w1 + 2*d2*i2*q2*r2*v0*w1 - 2*l2*n2*q2*r2*v0*w1 - 
   2*c2*j2*l2*s2*v0*w1 + 2*c2*i2*q2*s2*v0*w1 - 2*l2*m2*q2*s2*v0*w1 + 
   2*pow(l2,2)*r2*s2*v0*w1 + 2*c2*d0*pow(j2,2)*v2*w1 + 
   2*c0*d2*pow(j2,2)*v2*w1 - 2*c2*d0*i2*o2*v2*w1 - 2*c0*d2*i2*o2*v2*w1 - 
   2*f0*pow(l2,2)*o2*v2*w1 + 2*d0*l2*m2*o2*v2*w1 + 2*c2*l2*n0*o2*v2*w1 + 
   2*c0*l2*n2*o2*v2*w1 + 4*f0*j2*l2*q2*v2*w1 - 2*d0*j2*m2*q2*v2*w1 - 
   2*c2*j2*n0*q2*v2*w1 - 2*c0*j2*n2*q2*v2*w1 - 2*f0*i2*pow(q2,2)*v2*w1 + 
   2*m2*n0*pow(q2,2)*v2*w1 - 2*d2*j2*l2*r0*v2*w1 + 2*d2*i2*q2*r0*v2*w1 - 
   2*l2*n2*q2*r0*v2*w1 - 2*d0*j2*l2*r2*v2*w1 + 2*d0*i2*q2*r2*v2*w1 - 
   2*l2*n0*q2*r2*v2*w1 - 2*c2*j2*l2*s0*v2*w1 + 2*c2*i2*q2*s0*v2*w1 - 
   2*l2*m2*q2*s0*v2*w1 + 2*pow(l2,2)*r2*s0*v2*w1 - 2*c0*j2*l2*s2*v2*w1 + 
   2*c0*i2*q2*s2*v2*w1 + 2*pow(l2,2)*r0*s2*v2*w1 - 
   2*pow(c2,2)*pow(j2,2)*w0*w1 + 2*pow(c2,2)*i2*o2*w0*w1 + 
   2*e2*pow(l2,2)*o2*w0*w1 - 4*c2*l2*m2*o2*w0*w1 - 4*e2*j2*l2*q2*w0*w1 + 
   4*c2*j2*m2*q2*w0*w1 + 2*e2*i2*pow(q2,2)*w0*w1 - 
   2*pow(m2,2)*pow(q2,2)*w0*w1 + 4*c2*j2*l2*r2*w0*w1 - 
   4*c2*i2*q2*r2*w0*w1 + 4*l2*m2*q2*r2*w0*w1 - 
   2*pow(l2,2)*pow(r2,2)*w0*w1 - 2*d1*e0*k2*l2*o2*w2 - 
   2*d0*e1*k2*l2*o2*w2 + 2*c1*f0*k2*l2*o2*w2 + 2*c0*f1*k2*l2*o2*w2 + 
   2*c1*d0*k2*m2*o2*w2 + 2*c0*d1*k2*m2*o2*w2 - 4*c1*c2*k2*n0*o2*w2 - 
   4*c0*c2*k2*n1*o2*w2 - 4*c0*c1*k2*n2*o2*w2 + 2*d1*e0*j2*l2*p2*w2 + 
   2*d0*e1*j2*l2*p2*w2 - 2*c1*f0*j2*l2*p2*w2 - 2*c0*f1*j2*l2*p2*w2 - 
   2*c1*d0*j2*m2*p2*w2 - 2*c0*d1*j2*m2*p2*w2 + 4*c1*c2*j2*n0*p2*w2 + 
   4*c0*c2*j2*n1*p2*w2 + 4*c0*c1*j2*n2*p2*w2 + 2*d1*e0*j2*k2*q2*w2 + 
   2*d0*e1*j2*k2*q2*w2 - 2*c1*f0*j2*k2*q2*w2 - 2*c0*f1*j2*k2*q2*w2 - 
   2*d1*e0*i2*p2*q2*w2 - 2*d0*e1*i2*p2*q2*w2 + 2*c1*f0*i2*p2*q2*w2 + 
   2*c0*f1*i2*p2*q2*w2 + 2*e1*l2*n0*p2*q2*w2 - 2*c1*m2*n0*p2*q2*w2 + 
   2*e0*l2*n1*p2*q2*w2 - 2*c0*m2*n1*p2*q2*w2 - 2*e1*k2*n0*pow(q2,2)*w2 - 
   2*e0*k2*n1*pow(q2,2)*w2 - 2*c2*d1*j2*k2*r0*w2 - 2*c1*d2*j2*k2*r0*w2 + 
   2*c2*d1*i2*p2*r0*w2 + 2*c1*d2*i2*p2*r0*w2 + 2*f1*pow(l2,2)*p2*r0*w2 - 
   2*d1*l2*m2*p2*r0*w2 - 2*c2*l2*n1*p2*r0*w2 - 2*c1*l2*n2*p2*r0*w2 - 
   2*f1*k2*l2*q2*r0*w2 - 2*d1*k2*m2*q2*r0*w2 + 4*c2*k2*n1*q2*r0*w2 + 
   4*c1*k2*n2*q2*r0*w2 - 2*c2*d0*j2*k2*r1*w2 - 2*c0*d2*j2*k2*r1*w2 + 
   2*c2*d0*i2*p2*r1*w2 + 2*c0*d2*i2*p2*r1*w2 + 2*f0*pow(l2,2)*p2*r1*w2 - 
   2*d0*l2*m2*p2*r1*w2 - 2*c2*l2*n0*p2*r1*w2 - 2*c0*l2*n2*p2*r1*w2 - 
   2*f0*k2*l2*q2*r1*w2 - 2*d0*k2*m2*q2*r1*w2 + 4*c2*k2*n0*q2*r1*w2 + 
   4*c0*k2*n2*q2*r1*w2 + 4*d2*k2*l2*r0*r1*w2 - 2*c1*d0*j2*k2*r2*w2 - 
   2*c0*d1*j2*k2*r2*w2 + 2*c1*d0*i2*p2*r2*w2 + 2*c0*d1*i2*p2*r2*w2 - 
   2*c1*l2*n0*p2*r2*w2 - 2*c0*l2*n1*p2*r2*w2 + 4*c1*k2*n0*q2*r2*w2 + 
   4*c0*k2*n1*q2*r2*w2 + 4*d1*k2*l2*r0*r2*w2 + 4*d0*k2*l2*r1*r2*w2 + 
   4*c1*c2*j2*k2*s0*w2 - 4*c1*c2*i2*p2*s0*w2 - 2*e1*pow(l2,2)*p2*s0*w2 + 
   4*c1*l2*m2*p2*s0*w2 + 2*e1*k2*l2*q2*s0*w2 - 2*c1*k2*m2*q2*s0*w2 - 
   2*c2*k2*l2*r1*s0*w2 - 2*c1*k2*l2*r2*s0*w2 + 4*c0*c2*j2*k2*s1*w2 - 
   4*c0*c2*i2*p2*s1*w2 - 2*e0*pow(l2,2)*p2*s1*w2 + 4*c0*l2*m2*p2*s1*w2 + 
   2*e0*k2*l2*q2*s1*w2 - 2*c0*k2*m2*q2*s1*w2 - 2*c2*k2*l2*r0*s1*w2 - 
   2*c0*k2*l2*r2*s1*w2 + 4*c0*c1*j2*k2*s2*w2 - 4*c0*c1*i2*p2*s2*w2 - 
   2*c1*k2*l2*r0*s2*w2 - 2*c0*k2*l2*r1*s2*w2 - 2*d2*e1*pow(j2,2)*u0*w2 - 
   2*d1*e2*pow(j2,2)*u0*w2 + 2*c2*f1*pow(j2,2)*u0*w2 + 
   2*c1*f2*pow(j2,2)*u0*w2 + 2*d2*e1*i2*o2*u0*w2 + 2*d1*e2*i2*o2*u0*w2 - 
   2*c2*f1*i2*o2*u0*w2 - 2*c1*f2*i2*o2*u0*w2 + 2*f1*l2*m2*o2*u0*w2 - 
   2*d1*pow(m2,2)*o2*u0*w2 - 2*e2*l2*n1*o2*u0*w2 + 2*c2*m2*n1*o2*u0*w2 - 
   2*e1*l2*n2*o2*u0*w2 + 2*c1*m2*n2*o2*u0*w2 - 2*f1*j2*m2*q2*u0*w2 + 
   2*e2*j2*n1*q2*u0*w2 + 2*e1*j2*n2*q2*u0*w2 - 2*f2*j2*l2*r1*u0*w2 + 
   4*d2*j2*m2*r1*u0*w2 - 2*c2*j2*n2*r1*u0*w2 + 2*f2*i2*q2*r1*u0*w2 - 
   2*m2*n2*q2*r1*u0*w2 - 2*f1*j2*l2*r2*u0*w2 + 4*d1*j2*m2*r2*u0*w2 - 
   2*c2*j2*n1*r2*u0*w2 - 2*c1*j2*n2*r2*u0*w2 + 2*f1*i2*q2*r2*u0*w2 - 
   2*m2*n1*q2*r2*u0*w2 - 4*d2*i2*r1*r2*u0*w2 + 4*l2*n2*r1*r2*u0*w2 - 
   2*d1*i2*pow(r2,2)*u0*w2 + 2*l2*n1*pow(r2,2)*u0*w2 + 
   2*e2*j2*l2*s1*u0*w2 - 2*c2*j2*m2*s1*u0*w2 - 2*e2*i2*q2*s1*u0*w2 + 
   2*pow(m2,2)*q2*s1*u0*w2 + 2*c2*i2*r2*s1*u0*w2 - 2*l2*m2*r2*s1*u0*w2 + 
   2*e1*j2*l2*s2*u0*w2 - 2*c1*j2*m2*s2*u0*w2 - 2*e1*i2*q2*s2*u0*w2 + 
   2*c2*i2*r1*s2*u0*w2 - 2*l2*m2*r1*s2*u0*w2 + 2*c1*i2*r2*s2*u0*w2 - 
   2*d2*e0*pow(j2,2)*u1*w2 - 2*d0*e2*pow(j2,2)*u1*w2 + 
   2*c2*f0*pow(j2,2)*u1*w2 + 2*c0*f2*pow(j2,2)*u1*w2 + 
   2*d2*e0*i2*o2*u1*w2 + 2*d0*e2*i2*o2*u1*w2 - 2*c2*f0*i2*o2*u1*w2 - 
   2*c0*f2*i2*o2*u1*w2 + 2*f0*l2*m2*o2*u1*w2 - 2*d0*pow(m2,2)*o2*u1*w2 - 
   2*e2*l2*n0*o2*u1*w2 + 2*c2*m2*n0*o2*u1*w2 - 2*e0*l2*n2*o2*u1*w2 + 
   2*c0*m2*n2*o2*u1*w2 - 2*f0*j2*m2*q2*u1*w2 + 2*e2*j2*n0*q2*u1*w2 + 
   2*e0*j2*n2*q2*u1*w2 - 2*f2*j2*l2*r0*u1*w2 + 4*d2*j2*m2*r0*u1*w2 - 
   2*c2*j2*n2*r0*u1*w2 + 2*f2*i2*q2*r0*u1*w2 - 2*m2*n2*q2*r0*u1*w2 - 
   2*f0*j2*l2*r2*u1*w2 + 4*d0*j2*m2*r2*u1*w2 - 2*c2*j2*n0*r2*u1*w2 - 
   2*c0*j2*n2*r2*u1*w2 + 2*f0*i2*q2*r2*u1*w2 - 2*m2*n0*q2*r2*u1*w2 - 
   4*d2*i2*r0*r2*u1*w2 + 4*l2*n2*r0*r2*u1*w2 - 2*d0*i2*pow(r2,2)*u1*w2 + 
   2*l2*n0*pow(r2,2)*u1*w2 + 2*e2*j2*l2*s0*u1*w2 - 2*c2*j2*m2*s0*u1*w2 - 
   2*e2*i2*q2*s0*u1*w2 + 2*pow(m2,2)*q2*s0*u1*w2 + 2*c2*i2*r2*s0*u1*w2 - 
   2*l2*m2*r2*s0*u1*w2 + 2*e0*j2*l2*s2*u1*w2 - 2*c0*j2*m2*s2*u1*w2 - 
   2*e0*i2*q2*s2*u1*w2 + 2*c2*i2*r0*s2*u1*w2 - 2*l2*m2*r0*s2*u1*w2 + 
   2*c0*i2*r2*s2*u1*w2 - 2*d1*e0*pow(j2,2)*u2*w2 - 
   2*d0*e1*pow(j2,2)*u2*w2 + 2*c1*f0*pow(j2,2)*u2*w2 + 
   2*c0*f1*pow(j2,2)*u2*w2 + 2*d1*e0*i2*o2*u2*w2 + 2*d0*e1*i2*o2*u2*w2 - 
   2*c1*f0*i2*o2*u2*w2 - 2*c0*f1*i2*o2*u2*w2 - 2*e1*l2*n0*o2*u2*w2 + 
   2*c1*m2*n0*o2*u2*w2 - 2*e0*l2*n1*o2*u2*w2 + 2*c0*m2*n1*o2*u2*w2 + 
   2*e1*j2*n0*q2*u2*w2 + 2*e0*j2*n1*q2*u2*w2 - 2*f1*j2*l2*r0*u2*w2 + 
   4*d1*j2*m2*r0*u2*w2 - 2*c2*j2*n1*r0*u2*w2 - 2*c1*j2*n2*r0*u2*w2 + 
   2*f1*i2*q2*r0*u2*w2 - 2*m2*n1*q2*r0*u2*w2 - 2*f0*j2*l2*r1*u2*w2 + 
   4*d0*j2*m2*r1*u2*w2 - 2*c2*j2*n0*r1*u2*w2 - 2*c0*j2*n2*r1*u2*w2 + 
   2*f0*i2*q2*r1*u2*w2 - 2*m2*n0*q2*r1*u2*w2 - 4*d2*i2*r0*r1*u2*w2 + 
   4*l2*n2*r0*r1*u2*w2 - 2*c1*j2*n0*r2*u2*w2 - 2*c0*j2*n1*r2*u2*w2 - 
   4*d1*i2*r0*r2*u2*w2 + 4*l2*n1*r0*r2*u2*w2 - 4*d0*i2*r1*r2*u2*w2 + 
   4*l2*n0*r1*r2*u2*w2 + 2*e1*j2*l2*s0*u2*w2 - 2*c1*j2*m2*s0*u2*w2 - 
   2*e1*i2*q2*s0*u2*w2 + 2*c2*i2*r1*s0*u2*w2 - 2*l2*m2*r1*s0*u2*w2 + 
   2*c1*i2*r2*s0*u2*w2 + 2*e0*j2*l2*s1*u2*w2 - 2*c0*j2*m2*s1*u2*w2 - 
   2*e0*i2*q2*s1*u2*w2 + 2*c2*i2*r0*s1*u2*w2 - 2*l2*m2*r0*s1*u2*w2 + 
   2*c0*i2*r2*s1*u2*w2 + 2*c1*i2*r0*s2*u2*w2 + 2*c0*i2*r1*s2*u2*w2 + 
   2*c2*d1*pow(j2,2)*v0*w2 + 2*c1*d2*pow(j2,2)*v0*w2 - 
   2*c2*d1*i2*o2*v0*w2 - 2*c1*d2*i2*o2*v0*w2 - 2*f1*pow(l2,2)*o2*v0*w2 + 
   2*d1*l2*m2*o2*v0*w2 + 2*c2*l2*n1*o2*v0*w2 + 2*c1*l2*n2*o2*v0*w2 + 
   4*f1*j2*l2*q2*v0*w2 - 2*d1*j2*m2*q2*v0*w2 - 2*c2*j2*n1*q2*v0*w2 - 
   2*c1*j2*n2*q2*v0*w2 - 2*f1*i2*pow(q2,2)*v0*w2 + 
   2*m2*n1*pow(q2,2)*v0*w2 - 2*d2*j2*l2*r1*v0*w2 + 2*d2*i2*q2*r1*v0*w2 - 
   2*l2*n2*q2*r1*v0*w2 - 2*d1*j2*l2*r2*v0*w2 + 2*d1*i2*q2*r2*v0*w2 - 
   2*l2*n1*q2*r2*v0*w2 - 2*c2*j2*l2*s1*v0*w2 + 2*c2*i2*q2*s1*v0*w2 - 
   2*l2*m2*q2*s1*v0*w2 + 2*pow(l2,2)*r2*s1*v0*w2 - 2*c1*j2*l2*s2*v0*w2 + 
   2*c1*i2*q2*s2*v0*w2 + 2*pow(l2,2)*r1*s2*v0*w2 + 
   2*c2*d0*pow(j2,2)*v1*w2 + 2*c0*d2*pow(j2,2)*v1*w2 - 
   2*c2*d0*i2*o2*v1*w2 - 2*c0*d2*i2*o2*v1*w2 - 2*f0*pow(l2,2)*o2*v1*w2 + 
   2*d0*l2*m2*o2*v1*w2 + 2*c2*l2*n0*o2*v1*w2 + 2*c0*l2*n2*o2*v1*w2 + 
   4*f0*j2*l2*q2*v1*w2 - 2*d0*j2*m2*q2*v1*w2 - 2*c2*j2*n0*q2*v1*w2 - 
   2*c0*j2*n2*q2*v1*w2 - 2*f0*i2*pow(q2,2)*v1*w2 + 
   2*m2*n0*pow(q2,2)*v1*w2 - 2*d2*j2*l2*r0*v1*w2 + 2*d2*i2*q2*r0*v1*w2 - 
   2*l2*n2*q2*r0*v1*w2 - 2*d0*j2*l2*r2*v1*w2 + 2*d0*i2*q2*r2*v1*w2 - 
   2*l2*n0*q2*r2*v1*w2 - 2*c2*j2*l2*s0*v1*w2 + 2*c2*i2*q2*s0*v1*w2 - 
   2*l2*m2*q2*s0*v1*w2 + 2*pow(l2,2)*r2*s0*v1*w2 - 2*c0*j2*l2*s2*v1*w2 + 
   2*c0*i2*q2*s2*v1*w2 + 2*pow(l2,2)*r0*s2*v1*w2 + 
   2*c1*d0*pow(j2,2)*v2*w2 + 2*c0*d1*pow(j2,2)*v2*w2 - 
   2*c1*d0*i2*o2*v2*w2 - 2*c0*d1*i2*o2*v2*w2 + 2*c1*l2*n0*o2*v2*w2 + 
   2*c0*l2*n1*o2*v2*w2 - 2*c1*j2*n0*q2*v2*w2 - 2*c0*j2*n1*q2*v2*w2 - 
   2*d1*j2*l2*r0*v2*w2 + 2*d1*i2*q2*r0*v2*w2 - 2*l2*n1*q2*r0*v2*w2 - 
   2*d0*j2*l2*r1*v2*w2 + 2*d0*i2*q2*r1*v2*w2 - 2*l2*n0*q2*r1*v2*w2 - 
   2*c1*j2*l2*s0*v2*w2 + 2*c1*i2*q2*s0*v2*w2 + 2*pow(l2,2)*r1*s0*v2*w2 - 
   2*c0*j2*l2*s1*v2*w2 + 2*c0*i2*q2*s1*v2*w2 + 2*pow(l2,2)*r0*s1*v2*w2 - 
   4*c1*c2*pow(j2,2)*w0*w2 + 4*c1*c2*i2*o2*w0*w2 + 
   2*e1*pow(l2,2)*o2*w0*w2 - 4*c1*l2*m2*o2*w0*w2 - 4*e1*j2*l2*q2*w0*w2 + 
   4*c1*j2*m2*q2*w0*w2 + 2*e1*i2*pow(q2,2)*w0*w2 + 4*c2*j2*l2*r1*w0*w2 - 
   4*c2*i2*q2*r1*w0*w2 + 4*l2*m2*q2*r1*w0*w2 + 4*c1*j2*l2*r2*w0*w2 - 
   4*c1*i2*q2*r2*w0*w2 - 4*pow(l2,2)*r1*r2*w0*w2 - 
   4*c0*c2*pow(j2,2)*w1*w2 + 4*c0*c2*i2*o2*w1*w2 + 
   2*e0*pow(l2,2)*o2*w1*w2 - 4*c0*l2*m2*o2*w1*w2 - 4*e0*j2*l2*q2*w1*w2 + 
   4*c0*j2*m2*q2*w1*w2 + 2*e0*i2*pow(q2,2)*w1*w2 + 4*c2*j2*l2*r0*w1*w2 - 
   4*c2*i2*q2*r0*w1*w2 + 4*l2*m2*q2*r0*w1*w2 + 4*c0*j2*l2*r2*w1*w2 - 
   4*c0*i2*q2*r2*w1*w2 - 4*pow(l2,2)*r0*r2*w1*w2 - 
   2*c0*c1*pow(j2,2)*pow(w2,2) + 2*c0*c1*i2*o2*pow(w2,2) + 
   2*c1*j2*l2*r0*pow(w2,2) - 2*c1*i2*q2*r0*pow(w2,2) + 
   2*c0*j2*l2*r1*pow(w2,2) - 2*c0*i2*q2*r1*pow(w2,2) - 
   2*pow(l2,2)*r0*r1*pow(w2,2) + 2*f1*f2*pow(k2,2)*o2*z0 - 
   e2*g1*pow(k2,2)*o2*z0 - e1*g2*pow(k2,2)*o2*z0 - 4*f1*f2*j2*k2*p2*z0 + 
   2*e2*g1*j2*k2*p2*z0 + 2*e1*g2*j2*k2*p2*z0 + 2*f1*f2*i2*pow(p2,2)*z0 - 
   e2*g1*i2*pow(p2,2)*z0 - e1*g2*i2*pow(p2,2)*z0 + 
   g1*pow(m2,2)*pow(p2,2)*z0 - 2*f2*m2*n1*pow(p2,2)*z0 - 
   2*f1*m2*n2*pow(p2,2)*z0 + 2*e2*n1*n2*pow(p2,2)*z0 + 
   e1*pow(n2,2)*pow(p2,2)*z0 - 2*g2*k2*m2*p2*r1*z0 + 
   2*f2*k2*n2*p2*r1*z0 - 2*g1*k2*m2*p2*r2*z0 + 2*f2*k2*n1*p2*r2*z0 + 
   2*f1*k2*n2*p2*r2*z0 + 2*g2*pow(k2,2)*r1*r2*z0 + 
   g1*pow(k2,2)*pow(r2,2)*z0 + 2*f2*k2*m2*p2*s1*z0 - 
   2*e2*k2*n2*p2*s1*z0 - 2*f2*pow(k2,2)*r2*s1*z0 + 2*f1*k2*m2*p2*s2*z0 - 
   2*e2*k2*n1*p2*s2*z0 - 2*e1*k2*n2*p2*s2*z0 - 2*f2*pow(k2,2)*r1*s2*z0 - 
   2*f1*pow(k2,2)*r2*s2*z0 + 2*e2*pow(k2,2)*s1*s2*z0 + 
   e1*pow(k2,2)*pow(s2,2)*z0 + 2*f1*f2*pow(j2,2)*t2*z0 - 
   e2*g1*pow(j2,2)*t2*z0 - e1*g2*pow(j2,2)*t2*z0 - 2*f1*f2*i2*o2*t2*z0 + 
   e2*g1*i2*o2*t2*z0 + e1*g2*i2*o2*t2*z0 - g1*pow(m2,2)*o2*t2*z0 + 
   2*f2*m2*n1*o2*t2*z0 + 2*f1*m2*n2*o2*t2*z0 - 2*e2*n1*n2*o2*t2*z0 - 
   e1*pow(n2,2)*o2*t2*z0 + 2*g2*j2*m2*r1*t2*z0 - 2*f2*j2*n2*r1*t2*z0 + 
   2*g1*j2*m2*r2*t2*z0 - 2*f2*j2*n1*r2*t2*z0 - 2*f1*j2*n2*r2*t2*z0 - 
   2*g2*i2*r1*r2*t2*z0 + 2*pow(n2,2)*r1*r2*t2*z0 - 
   g1*i2*pow(r2,2)*t2*z0 + 2*n1*n2*pow(r2,2)*t2*z0 - 
   2*f2*j2*m2*s1*t2*z0 + 2*e2*j2*n2*s1*t2*z0 + 2*f2*i2*r2*s1*t2*z0 - 
   2*m2*n2*r2*s1*t2*z0 - 2*f1*j2*m2*s2*t2*z0 + 2*e2*j2*n1*s2*t2*z0 + 
   2*e1*j2*n2*s2*t2*z0 + 2*f2*i2*r1*s2*t2*z0 - 2*m2*n2*r1*s2*t2*z0 + 
   2*f1*i2*r2*s2*t2*z0 - 2*m2*n1*r2*s2*t2*z0 - 2*e2*i2*s1*s2*t2*z0 + 
   2*pow(m2,2)*s1*s2*t2*z0 - e1*i2*pow(s2,2)*t2*z0 + 
   2*g2*k2*m2*o2*v1*z0 - 2*f2*k2*n2*o2*v1*z0 - 2*g2*j2*m2*p2*v1*z0 + 
   2*f2*j2*n2*p2*v1*z0 - 2*g2*j2*k2*r2*v1*z0 + 2*g2*i2*p2*r2*v1*z0 - 
   2*pow(n2,2)*p2*r2*v1*z0 + 2*f2*j2*k2*s2*v1*z0 - 2*f2*i2*p2*s2*v1*z0 + 
   2*m2*n2*p2*s2*v1*z0 + 2*k2*n2*r2*s2*v1*z0 - 2*k2*m2*pow(s2,2)*v1*z0 + 
   2*g1*k2*m2*o2*v2*z0 - 2*f2*k2*n1*o2*v2*z0 - 2*f1*k2*n2*o2*v2*z0 - 
   2*g1*j2*m2*p2*v2*z0 + 2*f2*j2*n1*p2*v2*z0 + 2*f1*j2*n2*p2*v2*z0 - 
   2*g2*j2*k2*r1*v2*z0 + 2*g2*i2*p2*r1*v2*z0 - 2*pow(n2,2)*p2*r1*v2*z0 - 
   2*g1*j2*k2*r2*v2*z0 + 2*g1*i2*p2*r2*v2*z0 - 4*n1*n2*p2*r2*v2*z0 + 
   2*f2*j2*k2*s1*v2*z0 - 2*f2*i2*p2*s1*v2*z0 + 2*m2*n2*p2*s1*v2*z0 + 
   2*k2*n2*r2*s1*v2*z0 + 2*f1*j2*k2*s2*v2*z0 - 2*f1*i2*p2*s2*v2*z0 + 
   2*m2*n1*p2*s2*v2*z0 + 2*k2*n2*r1*s2*v2*z0 + 2*k2*n1*r2*s2*v2*z0 - 
   4*k2*m2*s1*s2*v2*z0 + 2*g2*pow(j2,2)*v1*v2*z0 - 2*g2*i2*o2*v1*v2*z0 + 
   2*pow(n2,2)*o2*v1*v2*z0 - 4*j2*n2*s2*v1*v2*z0 + 
   2*i2*pow(s2,2)*v1*v2*z0 + g1*pow(j2,2)*pow(v2,2)*z0 - 
   g1*i2*o2*pow(v2,2)*z0 + 2*n1*n2*o2*pow(v2,2)*z0 - 
   2*j2*n2*s1*pow(v2,2)*z0 - 2*j2*n1*s2*pow(v2,2)*z0 + 
   2*i2*s1*s2*pow(v2,2)*z0 - 2*f2*k2*m2*o2*w1*z0 + 2*e2*k2*n2*o2*w1*z0 + 
   2*f2*j2*m2*p2*w1*z0 - 2*e2*j2*n2*p2*w1*z0 + 2*f2*j2*k2*r2*w1*z0 - 
   2*f2*i2*p2*r2*w1*z0 + 2*m2*n2*p2*r2*w1*z0 - 2*k2*n2*pow(r2,2)*w1*z0 - 
   2*e2*j2*k2*s2*w1*z0 + 2*e2*i2*p2*s2*w1*z0 - 2*pow(m2,2)*p2*s2*w1*z0 + 
   2*k2*m2*r2*s2*w1*z0 - 2*f2*pow(j2,2)*v2*w1*z0 + 2*f2*i2*o2*v2*w1*z0 - 
   2*m2*n2*o2*v2*w1*z0 + 2*j2*n2*r2*v2*w1*z0 + 2*j2*m2*s2*v2*w1*z0 - 
   2*i2*r2*s2*v2*w1*z0 - 2*f1*k2*m2*o2*w2*z0 + 2*e2*k2*n1*o2*w2*z0 + 
   2*e1*k2*n2*o2*w2*z0 + 2*f1*j2*m2*p2*w2*z0 - 2*e2*j2*n1*p2*w2*z0 - 
   2*e1*j2*n2*p2*w2*z0 + 2*f2*j2*k2*r1*w2*z0 - 2*f2*i2*p2*r1*w2*z0 + 
   2*m2*n2*p2*r1*w2*z0 + 2*f1*j2*k2*r2*w2*z0 - 2*f1*i2*p2*r2*w2*z0 + 
   2*m2*n1*p2*r2*w2*z0 - 4*k2*n2*r1*r2*w2*z0 - 2*k2*n1*pow(r2,2)*w2*z0 - 
   2*e2*j2*k2*s1*w2*z0 + 2*e2*i2*p2*s1*w2*z0 - 2*pow(m2,2)*p2*s1*w2*z0 + 
   2*k2*m2*r2*s1*w2*z0 - 2*e1*j2*k2*s2*w2*z0 + 2*e1*i2*p2*s2*w2*z0 + 
   2*k2*m2*r1*s2*w2*z0 - 2*f2*pow(j2,2)*v1*w2*z0 + 2*f2*i2*o2*v1*w2*z0 - 
   2*m2*n2*o2*v1*w2*z0 + 2*j2*n2*r2*v1*w2*z0 + 2*j2*m2*s2*v1*w2*z0 - 
   2*i2*r2*s2*v1*w2*z0 - 2*f1*pow(j2,2)*v2*w2*z0 + 2*f1*i2*o2*v2*w2*z0 - 
   2*m2*n1*o2*v2*w2*z0 + 2*j2*n2*r1*v2*w2*z0 + 2*j2*n1*r2*v2*w2*z0 + 
   2*j2*m2*s1*v2*w2*z0 - 2*i2*r2*s1*v2*w2*z0 - 2*i2*r1*s2*v2*w2*z0 + 
   2*e2*pow(j2,2)*w1*w2*z0 - 2*e2*i2*o2*w1*w2*z0 + 
   2*pow(m2,2)*o2*w1*w2*z0 - 4*j2*m2*r2*w1*w2*z0 + 
   2*i2*pow(r2,2)*w1*w2*z0 + e1*pow(j2,2)*pow(w2,2)*z0 - 
   e1*i2*o2*pow(w2,2)*z0 - 2*j2*m2*r1*pow(w2,2)*z0 + 
   2*i2*r1*r2*pow(w2,2)*z0 + 2*f0*f2*pow(k2,2)*o2*z1 - 
   e2*g0*pow(k2,2)*o2*z1 - e0*g2*pow(k2,2)*o2*z1 - 4*f0*f2*j2*k2*p2*z1 + 
   2*e2*g0*j2*k2*p2*z1 + 2*e0*g2*j2*k2*p2*z1 + 2*f0*f2*i2*pow(p2,2)*z1 - 
   e2*g0*i2*pow(p2,2)*z1 - e0*g2*i2*pow(p2,2)*z1 + 
   g0*pow(m2,2)*pow(p2,2)*z1 - 2*f2*m2*n0*pow(p2,2)*z1 - 
   2*f0*m2*n2*pow(p2,2)*z1 + 2*e2*n0*n2*pow(p2,2)*z1 + 
   e0*pow(n2,2)*pow(p2,2)*z1 - 2*g2*k2*m2*p2*r0*z1 + 
   2*f2*k2*n2*p2*r0*z1 - 2*g0*k2*m2*p2*r2*z1 + 2*f2*k2*n0*p2*r2*z1 + 
   2*f0*k2*n2*p2*r2*z1 + 2*g2*pow(k2,2)*r0*r2*z1 + 
   g0*pow(k2,2)*pow(r2,2)*z1 + 2*f2*k2*m2*p2*s0*z1 - 
   2*e2*k2*n2*p2*s0*z1 - 2*f2*pow(k2,2)*r2*s0*z1 + 2*f0*k2*m2*p2*s2*z1 - 
   2*e2*k2*n0*p2*s2*z1 - 2*e0*k2*n2*p2*s2*z1 - 2*f2*pow(k2,2)*r0*s2*z1 - 
   2*f0*pow(k2,2)*r2*s2*z1 + 2*e2*pow(k2,2)*s0*s2*z1 + 
   e0*pow(k2,2)*pow(s2,2)*z1 + 2*f0*f2*pow(j2,2)*t2*z1 - 
   e2*g0*pow(j2,2)*t2*z1 - e0*g2*pow(j2,2)*t2*z1 - 2*f0*f2*i2*o2*t2*z1 + 
   e2*g0*i2*o2*t2*z1 + e0*g2*i2*o2*t2*z1 - g0*pow(m2,2)*o2*t2*z1 + 
   2*f2*m2*n0*o2*t2*z1 + 2*f0*m2*n2*o2*t2*z1 - 2*e2*n0*n2*o2*t2*z1 - 
   e0*pow(n2,2)*o2*t2*z1 + 2*g2*j2*m2*r0*t2*z1 - 2*f2*j2*n2*r0*t2*z1 + 
   2*g0*j2*m2*r2*t2*z1 - 2*f2*j2*n0*r2*t2*z1 - 2*f0*j2*n2*r2*t2*z1 - 
   2*g2*i2*r0*r2*t2*z1 + 2*pow(n2,2)*r0*r2*t2*z1 - 
   g0*i2*pow(r2,2)*t2*z1 + 2*n0*n2*pow(r2,2)*t2*z1 - 
   2*f2*j2*m2*s0*t2*z1 + 2*e2*j2*n2*s0*t2*z1 + 2*f2*i2*r2*s0*t2*z1 - 
   2*m2*n2*r2*s0*t2*z1 - 2*f0*j2*m2*s2*t2*z1 + 2*e2*j2*n0*s2*t2*z1 + 
   2*e0*j2*n2*s2*t2*z1 + 2*f2*i2*r0*s2*t2*z1 - 2*m2*n2*r0*s2*t2*z1 + 
   2*f0*i2*r2*s2*t2*z1 - 2*m2*n0*r2*s2*t2*z1 - 2*e2*i2*s0*s2*t2*z1 + 
   2*pow(m2,2)*s0*s2*t2*z1 - e0*i2*pow(s2,2)*t2*z1 + 
   2*g2*k2*m2*o2*v0*z1 - 2*f2*k2*n2*o2*v0*z1 - 2*g2*j2*m2*p2*v0*z1 + 
   2*f2*j2*n2*p2*v0*z1 - 2*g2*j2*k2*r2*v0*z1 + 2*g2*i2*p2*r2*v0*z1 - 
   2*pow(n2,2)*p2*r2*v0*z1 + 2*f2*j2*k2*s2*v0*z1 - 2*f2*i2*p2*s2*v0*z1 + 
   2*m2*n2*p2*s2*v0*z1 + 2*k2*n2*r2*s2*v0*z1 - 2*k2*m2*pow(s2,2)*v0*z1 + 
   2*g0*k2*m2*o2*v2*z1 - 2*f2*k2*n0*o2*v2*z1 - 2*f0*k2*n2*o2*v2*z1 - 
   2*g0*j2*m2*p2*v2*z1 + 2*f2*j2*n0*p2*v2*z1 + 2*f0*j2*n2*p2*v2*z1 - 
   2*g2*j2*k2*r0*v2*z1 + 2*g2*i2*p2*r0*v2*z1 - 2*pow(n2,2)*p2*r0*v2*z1 - 
   2*g0*j2*k2*r2*v2*z1 + 2*g0*i2*p2*r2*v2*z1 - 4*n0*n2*p2*r2*v2*z1 + 
   2*f2*j2*k2*s0*v2*z1 - 2*f2*i2*p2*s0*v2*z1 + 2*m2*n2*p2*s0*v2*z1 + 
   2*k2*n2*r2*s0*v2*z1 + 2*f0*j2*k2*s2*v2*z1 - 2*f0*i2*p2*s2*v2*z1 + 
   2*m2*n0*p2*s2*v2*z1 + 2*k2*n2*r0*s2*v2*z1 + 2*k2*n0*r2*s2*v2*z1 - 
   4*k2*m2*s0*s2*v2*z1 + 2*g2*pow(j2,2)*v0*v2*z1 - 2*g2*i2*o2*v0*v2*z1 + 
   2*pow(n2,2)*o2*v0*v2*z1 - 4*j2*n2*s2*v0*v2*z1 + 
   2*i2*pow(s2,2)*v0*v2*z1 + g0*pow(j2,2)*pow(v2,2)*z1 - 
   g0*i2*o2*pow(v2,2)*z1 + 2*n0*n2*o2*pow(v2,2)*z1 - 
   2*j2*n2*s0*pow(v2,2)*z1 - 2*j2*n0*s2*pow(v2,2)*z1 + 
   2*i2*s0*s2*pow(v2,2)*z1 - 2*f2*k2*m2*o2*w0*z1 + 2*e2*k2*n2*o2*w0*z1 + 
   2*f2*j2*m2*p2*w0*z1 - 2*e2*j2*n2*p2*w0*z1 + 2*f2*j2*k2*r2*w0*z1 - 
   2*f2*i2*p2*r2*w0*z1 + 2*m2*n2*p2*r2*w0*z1 - 2*k2*n2*pow(r2,2)*w0*z1 - 
   2*e2*j2*k2*s2*w0*z1 + 2*e2*i2*p2*s2*w0*z1 - 2*pow(m2,2)*p2*s2*w0*z1 + 
   2*k2*m2*r2*s2*w0*z1 - 2*f2*pow(j2,2)*v2*w0*z1 + 2*f2*i2*o2*v2*w0*z1 - 
   2*m2*n2*o2*v2*w0*z1 + 2*j2*n2*r2*v2*w0*z1 + 2*j2*m2*s2*v2*w0*z1 - 
   2*i2*r2*s2*v2*w0*z1 - 2*f0*k2*m2*o2*w2*z1 + 2*e2*k2*n0*o2*w2*z1 + 
   2*e0*k2*n2*o2*w2*z1 + 2*f0*j2*m2*p2*w2*z1 - 2*e2*j2*n0*p2*w2*z1 - 
   2*e0*j2*n2*p2*w2*z1 + 2*f2*j2*k2*r0*w2*z1 - 2*f2*i2*p2*r0*w2*z1 + 
   2*m2*n2*p2*r0*w2*z1 + 2*f0*j2*k2*r2*w2*z1 - 2*f0*i2*p2*r2*w2*z1 + 
   2*m2*n0*p2*r2*w2*z1 - 4*k2*n2*r0*r2*w2*z1 - 2*k2*n0*pow(r2,2)*w2*z1 - 
   2*e2*j2*k2*s0*w2*z1 + 2*e2*i2*p2*s0*w2*z1 - 2*pow(m2,2)*p2*s0*w2*z1 + 
   2*k2*m2*r2*s0*w2*z1 - 2*e0*j2*k2*s2*w2*z1 + 2*e0*i2*p2*s2*w2*z1 + 
   2*k2*m2*r0*s2*w2*z1 - 2*f2*pow(j2,2)*v0*w2*z1 + 2*f2*i2*o2*v0*w2*z1 - 
   2*m2*n2*o2*v0*w2*z1 + 2*j2*n2*r2*v0*w2*z1 + 2*j2*m2*s2*v0*w2*z1 - 
   2*i2*r2*s2*v0*w2*z1 - 2*f0*pow(j2,2)*v2*w2*z1 + 2*f0*i2*o2*v2*w2*z1 - 
   2*m2*n0*o2*v2*w2*z1 + 2*j2*n2*r0*v2*w2*z1 + 2*j2*n0*r2*v2*w2*z1 + 
   2*j2*m2*s0*v2*w2*z1 - 2*i2*r2*s0*v2*w2*z1 - 2*i2*r0*s2*v2*w2*z1 + 
   2*e2*pow(j2,2)*w0*w2*z1 - 2*e2*i2*o2*w0*w2*z1 + 
   2*pow(m2,2)*o2*w0*w2*z1 - 4*j2*m2*r2*w0*w2*z1 + 
   2*i2*pow(r2,2)*w0*w2*z1 + e0*pow(j2,2)*pow(w2,2)*z1 - 
   e0*i2*o2*pow(w2,2)*z1 - 2*j2*m2*r0*pow(w2,2)*z1 + 
   2*i2*r0*r2*pow(w2,2)*z1 + 2*f0*f1*pow(k2,2)*o2*z2 - 
   e1*g0*pow(k2,2)*o2*z2 - e0*g1*pow(k2,2)*o2*z2 - 4*f0*f1*j2*k2*p2*z2 + 
   2*e1*g0*j2*k2*p2*z2 + 2*e0*g1*j2*k2*p2*z2 + 2*f0*f1*i2*pow(p2,2)*z2 - 
   e1*g0*i2*pow(p2,2)*z2 - e0*g1*i2*pow(p2,2)*z2 - 
   2*f1*m2*n0*pow(p2,2)*z2 - 2*f0*m2*n1*pow(p2,2)*z2 + 
   2*e2*n0*n1*pow(p2,2)*z2 + 2*e1*n0*n2*pow(p2,2)*z2 + 
   2*e0*n1*n2*pow(p2,2)*z2 - 2*g1*k2*m2*p2*r0*z2 + 2*f2*k2*n1*p2*r0*z2 + 
   2*f1*k2*n2*p2*r0*z2 - 2*g0*k2*m2*p2*r1*z2 + 2*f2*k2*n0*p2*r1*z2 + 
   2*f0*k2*n2*p2*r1*z2 + 2*g2*pow(k2,2)*r0*r1*z2 + 2*f1*k2*n0*p2*r2*z2 + 
   2*f0*k2*n1*p2*r2*z2 + 2*g1*pow(k2,2)*r0*r2*z2 + 
   2*g0*pow(k2,2)*r1*r2*z2 + 2*f1*k2*m2*p2*s0*z2 - 2*e2*k2*n1*p2*s0*z2 - 
   2*e1*k2*n2*p2*s0*z2 - 2*f2*pow(k2,2)*r1*s0*z2 - 
   2*f1*pow(k2,2)*r2*s0*z2 + 2*f0*k2*m2*p2*s1*z2 - 2*e2*k2*n0*p2*s1*z2 - 
   2*e0*k2*n2*p2*s1*z2 - 2*f2*pow(k2,2)*r0*s1*z2 - 
   2*f0*pow(k2,2)*r2*s1*z2 + 2*e2*pow(k2,2)*s0*s1*z2 - 
   2*e1*k2*n0*p2*s2*z2 - 2*e0*k2*n1*p2*s2*z2 - 2*f1*pow(k2,2)*r0*s2*z2 - 
   2*f0*pow(k2,2)*r1*s2*z2 + 2*e1*pow(k2,2)*s0*s2*z2 + 
   2*e0*pow(k2,2)*s1*s2*z2 + 2*f0*f1*pow(j2,2)*t2*z2 - 
   e1*g0*pow(j2,2)*t2*z2 - e0*g1*pow(j2,2)*t2*z2 - 2*f0*f1*i2*o2*t2*z2 + 
   e1*g0*i2*o2*t2*z2 + e0*g1*i2*o2*t2*z2 + 2*f1*m2*n0*o2*t2*z2 + 
   2*f0*m2*n1*o2*t2*z2 - 2*e2*n0*n1*o2*t2*z2 - 2*e1*n0*n2*o2*t2*z2 - 
   2*e0*n1*n2*o2*t2*z2 + 2*g1*j2*m2*r0*t2*z2 - 2*f2*j2*n1*r0*t2*z2 - 
   2*f1*j2*n2*r0*t2*z2 + 2*g0*j2*m2*r1*t2*z2 - 2*f2*j2*n0*r1*t2*z2 - 
   2*f0*j2*n2*r1*t2*z2 - 2*g2*i2*r0*r1*t2*z2 + 2*pow(n2,2)*r0*r1*t2*z2 - 
   2*f1*j2*n0*r2*t2*z2 - 2*f0*j2*n1*r2*t2*z2 - 2*g1*i2*r0*r2*t2*z2 + 
   4*n1*n2*r0*r2*t2*z2 - 2*g0*i2*r1*r2*t2*z2 + 4*n0*n2*r1*r2*t2*z2 + 
   2*n0*n1*pow(r2,2)*t2*z2 - 2*f1*j2*m2*s0*t2*z2 + 2*e2*j2*n1*s0*t2*z2 + 
   2*e1*j2*n2*s0*t2*z2 + 2*f2*i2*r1*s0*t2*z2 - 2*m2*n2*r1*s0*t2*z2 + 
   2*f1*i2*r2*s0*t2*z2 - 2*m2*n1*r2*s0*t2*z2 - 2*f0*j2*m2*s1*t2*z2 + 
   2*e2*j2*n0*s1*t2*z2 + 2*e0*j2*n2*s1*t2*z2 + 2*f2*i2*r0*s1*t2*z2 - 
   2*m2*n2*r0*s1*t2*z2 + 2*f0*i2*r2*s1*t2*z2 - 2*m2*n0*r2*s1*t2*z2 - 
   2*e2*i2*s0*s1*t2*z2 + 2*pow(m2,2)*s0*s1*t2*z2 + 2*e1*j2*n0*s2*t2*z2 + 
   2*e0*j2*n1*s2*t2*z2 + 2*f1*i2*r0*s2*t2*z2 - 2*m2*n1*r0*s2*t2*z2 + 
   2*f0*i2*r1*s2*t2*z2 - 2*m2*n0*r1*s2*t2*z2 - 2*e1*i2*s0*s2*t2*z2 - 
   2*e0*i2*s1*s2*t2*z2 + 2*g1*k2*m2*o2*v0*z2 - 2*f2*k2*n1*o2*v0*z2 - 
   2*f1*k2*n2*o2*v0*z2 - 2*g1*j2*m2*p2*v0*z2 + 2*f2*j2*n1*p2*v0*z2 + 
   2*f1*j2*n2*p2*v0*z2 - 2*g2*j2*k2*r1*v0*z2 + 2*g2*i2*p2*r1*v0*z2 - 
   2*pow(n2,2)*p2*r1*v0*z2 - 2*g1*j2*k2*r2*v0*z2 + 2*g1*i2*p2*r2*v0*z2 - 
   4*n1*n2*p2*r2*v0*z2 + 2*f2*j2*k2*s1*v0*z2 - 2*f2*i2*p2*s1*v0*z2 + 
   2*m2*n2*p2*s1*v0*z2 + 2*k2*n2*r2*s1*v0*z2 + 2*f1*j2*k2*s2*v0*z2 - 
   2*f1*i2*p2*s2*v0*z2 + 2*m2*n1*p2*s2*v0*z2 + 2*k2*n2*r1*s2*v0*z2 + 
   2*k2*n1*r2*s2*v0*z2 - 4*k2*m2*s1*s2*v0*z2 + 2*g0*k2*m2*o2*v1*z2 - 
   2*f2*k2*n0*o2*v1*z2 - 2*f0*k2*n2*o2*v1*z2 - 2*g0*j2*m2*p2*v1*z2 + 
   2*f2*j2*n0*p2*v1*z2 + 2*f0*j2*n2*p2*v1*z2 - 2*g2*j2*k2*r0*v1*z2 + 
   2*g2*i2*p2*r0*v1*z2 - 2*pow(n2,2)*p2*r0*v1*z2 - 2*g0*j2*k2*r2*v1*z2 + 
   2*g0*i2*p2*r2*v1*z2 - 4*n0*n2*p2*r2*v1*z2 + 2*f2*j2*k2*s0*v1*z2 - 
   2*f2*i2*p2*s0*v1*z2 + 2*m2*n2*p2*s0*v1*z2 + 2*k2*n2*r2*s0*v1*z2 + 
   2*f0*j2*k2*s2*v1*z2 - 2*f0*i2*p2*s2*v1*z2 + 2*m2*n0*p2*s2*v1*z2 + 
   2*k2*n2*r0*s2*v1*z2 + 2*k2*n0*r2*s2*v1*z2 - 4*k2*m2*s0*s2*v1*z2 + 
   2*g2*pow(j2,2)*v0*v1*z2 - 2*g2*i2*o2*v0*v1*z2 + 
   2*pow(n2,2)*o2*v0*v1*z2 - 4*j2*n2*s2*v0*v1*z2 + 
   2*i2*pow(s2,2)*v0*v1*z2 - 2*f1*k2*n0*o2*v2*z2 - 2*f0*k2*n1*o2*v2*z2 + 
   2*f1*j2*n0*p2*v2*z2 + 2*f0*j2*n1*p2*v2*z2 - 2*g1*j2*k2*r0*v2*z2 + 
   2*g1*i2*p2*r0*v2*z2 - 4*n1*n2*p2*r0*v2*z2 - 2*g0*j2*k2*r1*v2*z2 + 
   2*g0*i2*p2*r1*v2*z2 - 4*n0*n2*p2*r1*v2*z2 - 4*n0*n1*p2*r2*v2*z2 + 
   2*f1*j2*k2*s0*v2*z2 - 2*f1*i2*p2*s0*v2*z2 + 2*m2*n1*p2*s0*v2*z2 + 
   2*k2*n2*r1*s0*v2*z2 + 2*k2*n1*r2*s0*v2*z2 + 2*f0*j2*k2*s1*v2*z2 - 
   2*f0*i2*p2*s1*v2*z2 + 2*m2*n0*p2*s1*v2*z2 + 2*k2*n2*r0*s1*v2*z2 + 
   2*k2*n0*r2*s1*v2*z2 - 4*k2*m2*s0*s1*v2*z2 + 2*k2*n1*r0*s2*v2*z2 + 
   2*k2*n0*r1*s2*v2*z2 + 2*g1*pow(j2,2)*v0*v2*z2 - 2*g1*i2*o2*v0*v2*z2 + 
   4*n1*n2*o2*v0*v2*z2 - 4*j2*n2*s1*v0*v2*z2 - 4*j2*n1*s2*v0*v2*z2 + 
   4*i2*s1*s2*v0*v2*z2 + 2*g0*pow(j2,2)*v1*v2*z2 - 2*g0*i2*o2*v1*v2*z2 + 
   4*n0*n2*o2*v1*v2*z2 - 4*j2*n2*s0*v1*v2*z2 - 4*j2*n0*s2*v1*v2*z2 + 
   4*i2*s0*s2*v1*v2*z2 + 2*n0*n1*o2*pow(v2,2)*z2 - 
   2*j2*n1*s0*pow(v2,2)*z2 - 2*j2*n0*s1*pow(v2,2)*z2 + 
   2*i2*s0*s1*pow(v2,2)*z2 - 2*f1*k2*m2*o2*w0*z2 + 2*e2*k2*n1*o2*w0*z2 + 
   2*e1*k2*n2*o2*w0*z2 + 2*f1*j2*m2*p2*w0*z2 - 2*e2*j2*n1*p2*w0*z2 - 
   2*e1*j2*n2*p2*w0*z2 + 2*f2*j2*k2*r1*w0*z2 - 2*f2*i2*p2*r1*w0*z2 + 
   2*m2*n2*p2*r1*w0*z2 + 2*f1*j2*k2*r2*w0*z2 - 2*f1*i2*p2*r2*w0*z2 + 
   2*m2*n1*p2*r2*w0*z2 - 4*k2*n2*r1*r2*w0*z2 - 2*k2*n1*pow(r2,2)*w0*z2 - 
   2*e2*j2*k2*s1*w0*z2 + 2*e2*i2*p2*s1*w0*z2 - 2*pow(m2,2)*p2*s1*w0*z2 + 
   2*k2*m2*r2*s1*w0*z2 - 2*e1*j2*k2*s2*w0*z2 + 2*e1*i2*p2*s2*w0*z2 + 
   2*k2*m2*r1*s2*w0*z2 - 2*f2*pow(j2,2)*v1*w0*z2 + 2*f2*i2*o2*v1*w0*z2 - 
   2*m2*n2*o2*v1*w0*z2 + 2*j2*n2*r2*v1*w0*z2 + 2*j2*m2*s2*v1*w0*z2 - 
   2*i2*r2*s2*v1*w0*z2 - 2*f1*pow(j2,2)*v2*w0*z2 + 2*f1*i2*o2*v2*w0*z2 - 
   2*m2*n1*o2*v2*w0*z2 + 2*j2*n2*r1*v2*w0*z2 + 2*j2*n1*r2*v2*w0*z2 + 
   2*j2*m2*s1*v2*w0*z2 - 2*i2*r2*s1*v2*w0*z2 - 2*i2*r1*s2*v2*w0*z2 - 
   2*f0*k2*m2*o2*w1*z2 + 2*e2*k2*n0*o2*w1*z2 + 2*e0*k2*n2*o2*w1*z2 + 
   2*f0*j2*m2*p2*w1*z2 - 2*e2*j2*n0*p2*w1*z2 - 2*e0*j2*n2*p2*w1*z2 + 
   2*f2*j2*k2*r0*w1*z2 - 2*f2*i2*p2*r0*w1*z2 + 2*m2*n2*p2*r0*w1*z2 + 
   2*f0*j2*k2*r2*w1*z2 - 2*f0*i2*p2*r2*w1*z2 + 2*m2*n0*p2*r2*w1*z2 - 
   4*k2*n2*r0*r2*w1*z2 - 2*k2*n0*pow(r2,2)*w1*z2 - 2*e2*j2*k2*s0*w1*z2 + 
   2*e2*i2*p2*s0*w1*z2 - 2*pow(m2,2)*p2*s0*w1*z2 + 2*k2*m2*r2*s0*w1*z2 - 
   2*e0*j2*k2*s2*w1*z2 + 2*e0*i2*p2*s2*w1*z2 + 2*k2*m2*r0*s2*w1*z2 - 
   2*f2*pow(j2,2)*v0*w1*z2 + 2*f2*i2*o2*v0*w1*z2 - 2*m2*n2*o2*v0*w1*z2 + 
   2*j2*n2*r2*v0*w1*z2 + 2*j2*m2*s2*v0*w1*z2 - 2*i2*r2*s2*v0*w1*z2 - 
   2*f0*pow(j2,2)*v2*w1*z2 + 2*f0*i2*o2*v2*w1*z2 - 2*m2*n0*o2*v2*w1*z2 + 
   2*j2*n2*r0*v2*w1*z2 + 2*j2*n0*r2*v2*w1*z2 + 2*j2*m2*s0*v2*w1*z2 - 
   2*i2*r2*s0*v2*w1*z2 - 2*i2*r0*s2*v2*w1*z2 + 2*e2*pow(j2,2)*w0*w1*z2 - 
   2*e2*i2*o2*w0*w1*z2 + 2*pow(m2,2)*o2*w0*w1*z2 - 4*j2*m2*r2*w0*w1*z2 + 
   2*i2*pow(r2,2)*w0*w1*z2 + 2*e1*k2*n0*o2*w2*z2 + 2*e0*k2*n1*o2*w2*z2 - 
   2*e1*j2*n0*p2*w2*z2 - 2*e0*j2*n1*p2*w2*z2 + 2*f1*j2*k2*r0*w2*z2 - 
   2*f1*i2*p2*r0*w2*z2 + 2*m2*n1*p2*r0*w2*z2 + 2*f0*j2*k2*r1*w2*z2 - 
   2*f0*i2*p2*r1*w2*z2 + 2*m2*n0*p2*r1*w2*z2 - 4*k2*n2*r0*r1*w2*z2 - 
   4*k2*n1*r0*r2*w2*z2 - 4*k2*n0*r1*r2*w2*z2 - 2*e1*j2*k2*s0*w2*z2 + 
   2*e1*i2*p2*s0*w2*z2 + 2*k2*m2*r1*s0*w2*z2 - 2*e0*j2*k2*s1*w2*z2 + 
   2*e0*i2*p2*s1*w2*z2 + 2*k2*m2*r0*s1*w2*z2 - 2*f1*pow(j2,2)*v0*w2*z2 + 
   2*f1*i2*o2*v0*w2*z2 - 2*m2*n1*o2*v0*w2*z2 + 2*j2*n2*r1*v0*w2*z2 + 
   2*j2*n1*r2*v0*w2*z2 + 2*j2*m2*s1*v0*w2*z2 - 2*i2*r2*s1*v0*w2*z2 - 
   2*i2*r1*s2*v0*w2*z2 - 2*f0*pow(j2,2)*v1*w2*z2 + 2*f0*i2*o2*v1*w2*z2 - 
   2*m2*n0*o2*v1*w2*z2 + 2*j2*n2*r0*v1*w2*z2 + 2*j2*n0*r2*v1*w2*z2 + 
   2*j2*m2*s0*v1*w2*z2 - 2*i2*r2*s0*v1*w2*z2 - 2*i2*r0*s2*v1*w2*z2 + 
   2*j2*n1*r0*v2*w2*z2 + 2*j2*n0*r1*v2*w2*z2 - 2*i2*r1*s0*v2*w2*z2 - 
   2*i2*r0*s1*v2*w2*z2 + 2*e1*pow(j2,2)*w0*w2*z2 - 2*e1*i2*o2*w0*w2*z2 - 
   4*j2*m2*r1*w0*w2*z2 + 4*i2*r1*r2*w0*w2*z2 + 2*e0*pow(j2,2)*w1*w2*z2 - 
   2*e0*i2*o2*w1*w2*z2 - 4*j2*m2*r0*w1*w2*z2 + 4*i2*r0*r2*w1*w2*z2 + 
   2*i2*r0*r1*pow(w2,2)*z2
;

	k[5]  = 2*d1*d2*e1*pow(k2,2)*o2 + pow(d1,2)*e2*pow(k2,2)*o2 - 
   2*c2*d1*f1*pow(k2,2)*o2 - 2*c1*d2*f1*pow(k2,2)*o2 - 
   2*c1*d1*f2*pow(k2,2)*o2 + 2*c1*c2*g1*pow(k2,2)*o2 + 
   pow(c1,2)*g2*pow(k2,2)*o2 - 4*d1*d2*e1*j2*k2*p2 - 
   2*pow(d1,2)*e2*j2*k2*p2 + 4*c2*d1*f1*j2*k2*p2 + 4*c1*d2*f1*j2*k2*p2 + 
   4*c1*d1*f2*j2*k2*p2 - 4*c1*c2*g1*j2*k2*p2 - 2*pow(c1,2)*g2*j2*k2*p2 + 
   2*d1*d2*e1*i2*pow(p2,2) + pow(d1,2)*e2*i2*pow(p2,2) - 
   2*c2*d1*f1*i2*pow(p2,2) - 2*c1*d2*f1*i2*pow(p2,2) - 
   2*c1*d1*f2*i2*pow(p2,2) + 2*c1*c2*g1*i2*pow(p2,2) + 
   pow(c1,2)*g2*i2*pow(p2,2) - pow(f1,2)*pow(l2,2)*pow(p2,2) + 
   e1*g1*pow(l2,2)*pow(p2,2) + 2*d1*f1*l2*m2*pow(p2,2) - 
   2*c1*g1*l2*m2*pow(p2,2) - pow(d1,2)*pow(m2,2)*pow(p2,2) - 
   2*d2*e1*l2*n1*pow(p2,2) - 2*d1*e2*l2*n1*pow(p2,2) + 
   2*c2*f1*l2*n1*pow(p2,2) + 2*c1*f2*l2*n1*pow(p2,2) + 
   2*c2*d1*m2*n1*pow(p2,2) + 2*c1*d2*m2*n1*pow(p2,2) - 
   pow(c2,2)*pow(n1,2)*pow(p2,2) - 2*d1*e1*l2*n2*pow(p2,2) + 
   2*c1*f1*l2*n2*pow(p2,2) + 2*c1*d1*m2*n2*pow(p2,2) - 
   4*c1*c2*n1*n2*pow(p2,2) - pow(c1,2)*pow(n2,2)*pow(p2,2) + 
   2*pow(f1,2)*k2*l2*p2*q2 - 2*e1*g1*k2*l2*p2*q2 - 2*d1*f1*k2*m2*p2*q2 + 
   2*c1*g1*k2*m2*p2*q2 + 2*d2*e1*k2*n1*p2*q2 + 2*d1*e2*k2*n1*p2*q2 - 
   2*c2*f1*k2*n1*p2*q2 - 2*c1*f2*k2*n1*p2*q2 + 2*d1*e1*k2*n2*p2*q2 - 
   2*c1*f1*k2*n2*p2*q2 - pow(f1,2)*pow(k2,2)*pow(q2,2) + 
   e1*g1*pow(k2,2)*pow(q2,2) - 2*d2*f1*k2*l2*p2*r1 - 
   2*d1*f2*k2*l2*p2*r1 + 2*c2*g1*k2*l2*p2*r1 + 2*c1*g2*k2*l2*p2*r1 + 
   4*d1*d2*k2*m2*p2*r1 - 2*c2*d2*k2*n1*p2*r1 - 2*c2*d1*k2*n2*p2*r1 - 
   2*c1*d2*k2*n2*p2*r1 + 2*d2*f1*pow(k2,2)*q2*r1 + 
   2*d1*f2*pow(k2,2)*q2*r1 - 2*c2*g1*pow(k2,2)*q2*r1 - 
   2*c1*g2*pow(k2,2)*q2*r1 - pow(d2,2)*pow(k2,2)*pow(r1,2) - 
   2*d1*f1*k2*l2*p2*r2 + 2*c1*g1*k2*l2*p2*r2 + 2*pow(d1,2)*k2*m2*p2*r2 - 
   2*c2*d1*k2*n1*p2*r2 - 2*c1*d2*k2*n1*p2*r2 - 2*c1*d1*k2*n2*p2*r2 + 
   2*d1*f1*pow(k2,2)*q2*r2 - 2*c1*g1*pow(k2,2)*q2*r2 - 
   4*d1*d2*pow(k2,2)*r1*r2 - pow(d1,2)*pow(k2,2)*pow(r2,2) + 
   2*d2*e1*k2*l2*p2*s1 + 2*d1*e2*k2*l2*p2*s1 - 2*c2*f1*k2*l2*p2*s1 - 
   2*c1*f2*k2*l2*p2*s1 - 2*c2*d1*k2*m2*p2*s1 - 2*c1*d2*k2*m2*p2*s1 + 
   2*pow(c2,2)*k2*n1*p2*s1 + 4*c1*c2*k2*n2*p2*s1 - 
   2*d2*e1*pow(k2,2)*q2*s1 - 2*d1*e2*pow(k2,2)*q2*s1 + 
   2*c2*f1*pow(k2,2)*q2*s1 + 2*c1*f2*pow(k2,2)*q2*s1 + 
   2*c2*d2*pow(k2,2)*r1*s1 + 2*c2*d1*pow(k2,2)*r2*s1 + 
   2*c1*d2*pow(k2,2)*r2*s1 - pow(c2,2)*pow(k2,2)*pow(s1,2) + 
   2*d1*e1*k2*l2*p2*s2 - 2*c1*f1*k2*l2*p2*s2 - 2*c1*d1*k2*m2*p2*s2 + 
   4*c1*c2*k2*n1*p2*s2 + 2*pow(c1,2)*k2*n2*p2*s2 - 
   2*d1*e1*pow(k2,2)*q2*s2 + 2*c1*f1*pow(k2,2)*q2*s2 + 
   2*c2*d1*pow(k2,2)*r1*s2 + 2*c1*d2*pow(k2,2)*r1*s2 + 
   2*c1*d1*pow(k2,2)*r2*s2 - 4*c1*c2*pow(k2,2)*s1*s2 - 
   pow(c1,2)*pow(k2,2)*pow(s2,2) + 2*d1*d2*e1*pow(j2,2)*t2 + 
   pow(d1,2)*e2*pow(j2,2)*t2 - 2*c2*d1*f1*pow(j2,2)*t2 - 
   2*c1*d2*f1*pow(j2,2)*t2 - 2*c1*d1*f2*pow(j2,2)*t2 + 
   2*c1*c2*g1*pow(j2,2)*t2 + pow(c1,2)*g2*pow(j2,2)*t2 - 
   2*d1*d2*e1*i2*o2*t2 - pow(d1,2)*e2*i2*o2*t2 + 2*c2*d1*f1*i2*o2*t2 + 
   2*c1*d2*f1*i2*o2*t2 + 2*c1*d1*f2*i2*o2*t2 - 2*c1*c2*g1*i2*o2*t2 - 
   pow(c1,2)*g2*i2*o2*t2 + pow(f1,2)*pow(l2,2)*o2*t2 - 
   e1*g1*pow(l2,2)*o2*t2 - 2*d1*f1*l2*m2*o2*t2 + 2*c1*g1*l2*m2*o2*t2 + 
   pow(d1,2)*pow(m2,2)*o2*t2 + 2*d2*e1*l2*n1*o2*t2 + 
   2*d1*e2*l2*n1*o2*t2 - 2*c2*f1*l2*n1*o2*t2 - 2*c1*f2*l2*n1*o2*t2 - 
   2*c2*d1*m2*n1*o2*t2 - 2*c1*d2*m2*n1*o2*t2 + 
   pow(c2,2)*pow(n1,2)*o2*t2 + 2*d1*e1*l2*n2*o2*t2 - 
   2*c1*f1*l2*n2*o2*t2 - 2*c1*d1*m2*n2*o2*t2 + 4*c1*c2*n1*n2*o2*t2 + 
   pow(c1,2)*pow(n2,2)*o2*t2 - 2*pow(f1,2)*j2*l2*q2*t2 + 
   2*e1*g1*j2*l2*q2*t2 + 2*d1*f1*j2*m2*q2*t2 - 2*c1*g1*j2*m2*q2*t2 - 
   2*d2*e1*j2*n1*q2*t2 - 2*d1*e2*j2*n1*q2*t2 + 2*c2*f1*j2*n1*q2*t2 + 
   2*c1*f2*j2*n1*q2*t2 - 2*d1*e1*j2*n2*q2*t2 + 2*c1*f1*j2*n2*q2*t2 + 
   pow(f1,2)*i2*pow(q2,2)*t2 - e1*g1*i2*pow(q2,2)*t2 - 
   2*f1*m2*n1*pow(q2,2)*t2 + e2*pow(n1,2)*pow(q2,2)*t2 + 
   2*e1*n1*n2*pow(q2,2)*t2 + 2*d2*f1*j2*l2*r1*t2 + 2*d1*f2*j2*l2*r1*t2 - 
   2*c2*g1*j2*l2*r1*t2 - 2*c1*g2*j2*l2*r1*t2 - 4*d1*d2*j2*m2*r1*t2 + 
   2*c2*d2*j2*n1*r1*t2 + 2*c2*d1*j2*n2*r1*t2 + 2*c1*d2*j2*n2*r1*t2 - 
   2*d2*f1*i2*q2*r1*t2 - 2*d1*f2*i2*q2*r1*t2 + 2*c2*g1*i2*q2*r1*t2 + 
   2*c1*g2*i2*q2*r1*t2 - 2*g1*l2*m2*q2*r1*t2 + 2*f2*l2*n1*q2*r1*t2 + 
   2*d2*m2*n1*q2*r1*t2 + 2*f1*l2*n2*q2*r1*t2 + 2*d1*m2*n2*q2*r1*t2 - 
   4*c2*n1*n2*q2*r1*t2 - 2*c1*pow(n2,2)*q2*r1*t2 + 
   pow(d2,2)*i2*pow(r1,2)*t2 + g2*pow(l2,2)*pow(r1,2)*t2 - 
   2*d2*l2*n2*pow(r1,2)*t2 + 2*d1*f1*j2*l2*r2*t2 - 2*c1*g1*j2*l2*r2*t2 - 
   2*pow(d1,2)*j2*m2*r2*t2 + 2*c2*d1*j2*n1*r2*t2 + 2*c1*d2*j2*n1*r2*t2 + 
   2*c1*d1*j2*n2*r2*t2 - 2*d1*f1*i2*q2*r2*t2 + 2*c1*g1*i2*q2*r2*t2 + 
   2*f1*l2*n1*q2*r2*t2 + 2*d1*m2*n1*q2*r2*t2 - 2*c2*pow(n1,2)*q2*r2*t2 - 
   4*c1*n1*n2*q2*r2*t2 + 4*d1*d2*i2*r1*r2*t2 + 2*g1*pow(l2,2)*r1*r2*t2 - 
   4*d2*l2*n1*r1*r2*t2 - 4*d1*l2*n2*r1*r2*t2 + 
   pow(d1,2)*i2*pow(r2,2)*t2 - 2*d1*l2*n1*pow(r2,2)*t2 - 
   2*d2*e1*j2*l2*s1*t2 - 2*d1*e2*j2*l2*s1*t2 + 2*c2*f1*j2*l2*s1*t2 + 
   2*c1*f2*j2*l2*s1*t2 + 2*c2*d1*j2*m2*s1*t2 + 2*c1*d2*j2*m2*s1*t2 - 
   2*pow(c2,2)*j2*n1*s1*t2 - 4*c1*c2*j2*n2*s1*t2 + 2*d2*e1*i2*q2*s1*t2 + 
   2*d1*e2*i2*q2*s1*t2 - 2*c2*f1*i2*q2*s1*t2 - 2*c1*f2*i2*q2*s1*t2 + 
   2*f1*l2*m2*q2*s1*t2 - 2*d1*pow(m2,2)*q2*s1*t2 - 2*e2*l2*n1*q2*s1*t2 + 
   2*c2*m2*n1*q2*s1*t2 - 2*e1*l2*n2*q2*s1*t2 + 2*c1*m2*n2*q2*s1*t2 - 
   2*c2*d2*i2*r1*s1*t2 - 2*f2*pow(l2,2)*r1*s1*t2 + 2*d2*l2*m2*r1*s1*t2 + 
   2*c2*l2*n2*r1*s1*t2 - 2*c2*d1*i2*r2*s1*t2 - 2*c1*d2*i2*r2*s1*t2 - 
   2*f1*pow(l2,2)*r2*s1*t2 + 2*d1*l2*m2*r2*s1*t2 + 2*c2*l2*n1*r2*s1*t2 + 
   2*c1*l2*n2*r2*s1*t2 + pow(c2,2)*i2*pow(s1,2)*t2 + 
   e2*pow(l2,2)*pow(s1,2)*t2 - 2*c2*l2*m2*pow(s1,2)*t2 - 
   2*d1*e1*j2*l2*s2*t2 + 2*c1*f1*j2*l2*s2*t2 + 2*c1*d1*j2*m2*s2*t2 - 
   4*c1*c2*j2*n1*s2*t2 - 2*pow(c1,2)*j2*n2*s2*t2 + 2*d1*e1*i2*q2*s2*t2 - 
   2*c1*f1*i2*q2*s2*t2 - 2*e1*l2*n1*q2*s2*t2 + 2*c1*m2*n1*q2*s2*t2 - 
   2*c2*d1*i2*r1*s2*t2 - 2*c1*d2*i2*r1*s2*t2 - 2*f1*pow(l2,2)*r1*s2*t2 + 
   2*d1*l2*m2*r1*s2*t2 + 2*c2*l2*n1*r1*s2*t2 + 2*c1*l2*n2*r1*s2*t2 - 
   2*c1*d1*i2*r2*s2*t2 + 2*c1*l2*n1*r2*s2*t2 + 4*c1*c2*i2*s1*s2*t2 + 
   2*e1*pow(l2,2)*s1*s2*t2 - 4*c1*l2*m2*s1*s2*t2 + 
   pow(c1,2)*i2*pow(s2,2)*t2 - 4*f1*f2*k2*l2*o2*u1 + 
   2*e2*g1*k2*l2*o2*u1 + 2*e1*g2*k2*l2*o2*u1 + 2*d2*f1*k2*m2*o2*u1 + 
   2*d1*f2*k2*m2*o2*u1 - 2*c2*g1*k2*m2*o2*u1 - 2*c1*g2*k2*m2*o2*u1 - 
   2*d2*e2*k2*n1*o2*u1 + 2*c2*f2*k2*n1*o2*u1 - 2*d2*e1*k2*n2*o2*u1 - 
   2*d1*e2*k2*n2*o2*u1 + 2*c2*f1*k2*n2*o2*u1 + 2*c1*f2*k2*n2*o2*u1 + 
   4*f1*f2*j2*l2*p2*u1 - 2*e2*g1*j2*l2*p2*u1 - 2*e1*g2*j2*l2*p2*u1 - 
   2*d2*f1*j2*m2*p2*u1 - 2*d1*f2*j2*m2*p2*u1 + 2*c2*g1*j2*m2*p2*u1 + 
   2*c1*g2*j2*m2*p2*u1 + 2*d2*e2*j2*n1*p2*u1 - 2*c2*f2*j2*n1*p2*u1 + 
   2*d2*e1*j2*n2*p2*u1 + 2*d1*e2*j2*n2*p2*u1 - 2*c2*f1*j2*n2*p2*u1 - 
   2*c1*f2*j2*n2*p2*u1 + 4*f1*f2*j2*k2*q2*u1 - 2*e2*g1*j2*k2*q2*u1 - 
   2*e1*g2*j2*k2*q2*u1 - 4*f1*f2*i2*p2*q2*u1 + 2*e2*g1*i2*p2*q2*u1 + 
   2*e1*g2*i2*p2*q2*u1 - 2*g1*pow(m2,2)*p2*q2*u1 + 4*f2*m2*n1*p2*q2*u1 + 
   4*f1*m2*n2*p2*q2*u1 - 4*e2*n1*n2*p2*q2*u1 - 2*e1*pow(n2,2)*p2*q2*u1 - 
   2*d2*f2*j2*k2*r1*u1 + 2*c2*g2*j2*k2*r1*u1 + 2*d2*f2*i2*p2*r1*u1 - 
   2*c2*g2*i2*p2*r1*u1 + 2*g2*l2*m2*p2*r1*u1 - 2*f2*l2*n2*p2*r1*u1 - 
   2*d2*m2*n2*p2*r1*u1 + 2*c2*pow(n2,2)*p2*r1*u1 + 2*g2*k2*m2*q2*r1*u1 - 
   2*f2*k2*n2*q2*r1*u1 - 2*d2*f1*j2*k2*r2*u1 - 2*d1*f2*j2*k2*r2*u1 + 
   2*c2*g1*j2*k2*r2*u1 + 2*c1*g2*j2*k2*r2*u1 + 2*d2*f1*i2*p2*r2*u1 + 
   2*d1*f2*i2*p2*r2*u1 - 2*c2*g1*i2*p2*r2*u1 - 2*c1*g2*i2*p2*r2*u1 + 
   2*g1*l2*m2*p2*r2*u1 - 2*f2*l2*n1*p2*r2*u1 - 2*d2*m2*n1*p2*r2*u1 - 
   2*f1*l2*n2*p2*r2*u1 - 2*d1*m2*n2*p2*r2*u1 + 4*c2*n1*n2*p2*r2*u1 + 
   2*c1*pow(n2,2)*p2*r2*u1 + 2*g1*k2*m2*q2*r2*u1 - 2*f2*k2*n1*q2*r2*u1 - 
   2*f1*k2*n2*q2*r2*u1 - 4*g2*k2*l2*r1*r2*u1 + 4*d2*k2*n2*r1*r2*u1 - 
   2*g1*k2*l2*pow(r2,2)*u1 + 2*d2*k2*n1*pow(r2,2)*u1 + 
   2*d1*k2*n2*pow(r2,2)*u1 + 2*d2*e2*j2*k2*s1*u1 - 2*c2*f2*j2*k2*s1*u1 - 
   2*d2*e2*i2*p2*s1*u1 + 2*c2*f2*i2*p2*s1*u1 - 2*f2*l2*m2*p2*s1*u1 + 
   2*d2*pow(m2,2)*p2*s1*u1 + 2*e2*l2*n2*p2*s1*u1 - 2*c2*m2*n2*p2*s1*u1 - 
   2*f2*k2*m2*q2*s1*u1 + 2*e2*k2*n2*q2*s1*u1 + 4*f2*k2*l2*r2*s1*u1 - 
   2*d2*k2*m2*r2*s1*u1 - 2*c2*k2*n2*r2*s1*u1 + 2*d2*e1*j2*k2*s2*u1 + 
   2*d1*e2*j2*k2*s2*u1 - 2*c2*f1*j2*k2*s2*u1 - 2*c1*f2*j2*k2*s2*u1 - 
   2*d2*e1*i2*p2*s2*u1 - 2*d1*e2*i2*p2*s2*u1 + 2*c2*f1*i2*p2*s2*u1 + 
   2*c1*f2*i2*p2*s2*u1 - 2*f1*l2*m2*p2*s2*u1 + 2*d1*pow(m2,2)*p2*s2*u1 + 
   2*e2*l2*n1*p2*s2*u1 - 2*c2*m2*n1*p2*s2*u1 + 2*e1*l2*n2*p2*s2*u1 - 
   2*c1*m2*n2*p2*s2*u1 - 2*f1*k2*m2*q2*s2*u1 + 2*e2*k2*n1*q2*s2*u1 + 
   2*e1*k2*n2*q2*s2*u1 + 4*f2*k2*l2*r1*s2*u1 - 2*d2*k2*m2*r1*s2*u1 - 
   2*c2*k2*n2*r1*s2*u1 + 4*f1*k2*l2*r2*s2*u1 - 2*d1*k2*m2*r2*s2*u1 - 
   2*c2*k2*n1*r2*s2*u1 - 2*c1*k2*n2*r2*s2*u1 - 4*e2*k2*l2*s1*s2*u1 + 
   4*c2*k2*m2*s1*s2*u1 - 2*e1*k2*l2*pow(s2,2)*u1 + 
   2*c1*k2*m2*pow(s2,2)*u1 - pow(f2,2)*pow(j2,2)*pow(u1,2) + 
   e2*g2*pow(j2,2)*pow(u1,2) + pow(f2,2)*i2*o2*pow(u1,2) - 
   e2*g2*i2*o2*pow(u1,2) + g2*pow(m2,2)*o2*pow(u1,2) - 
   2*f2*m2*n2*o2*pow(u1,2) + e2*pow(n2,2)*o2*pow(u1,2) - 
   2*g2*j2*m2*r2*pow(u1,2) + 2*f2*j2*n2*r2*pow(u1,2) + 
   g2*i2*pow(r2,2)*pow(u1,2) - pow(n2,2)*pow(r2,2)*pow(u1,2) + 
   2*f2*j2*m2*s2*pow(u1,2) - 2*e2*j2*n2*s2*pow(u1,2) - 
   2*f2*i2*r2*s2*pow(u1,2) + 2*m2*n2*r2*s2*pow(u1,2) + 
   e2*i2*pow(s2,2)*pow(u1,2) - pow(m2,2)*pow(s2,2)*pow(u1,2) - 
   2*pow(f1,2)*k2*l2*o2*u2 + 2*e1*g1*k2*l2*o2*u2 + 2*d1*f1*k2*m2*o2*u2 - 
   2*c1*g1*k2*m2*o2*u2 - 2*d2*e1*k2*n1*o2*u2 - 2*d1*e2*k2*n1*o2*u2 + 
   2*c2*f1*k2*n1*o2*u2 + 2*c1*f2*k2*n1*o2*u2 - 2*d1*e1*k2*n2*o2*u2 + 
   2*c1*f1*k2*n2*o2*u2 + 2*pow(f1,2)*j2*l2*p2*u2 - 2*e1*g1*j2*l2*p2*u2 - 
   2*d1*f1*j2*m2*p2*u2 + 2*c1*g1*j2*m2*p2*u2 + 2*d2*e1*j2*n1*p2*u2 + 
   2*d1*e2*j2*n1*p2*u2 - 2*c2*f1*j2*n1*p2*u2 - 2*c1*f2*j2*n1*p2*u2 + 
   2*d1*e1*j2*n2*p2*u2 - 2*c1*f1*j2*n2*p2*u2 + 2*pow(f1,2)*j2*k2*q2*u2 - 
   2*e1*g1*j2*k2*q2*u2 - 2*pow(f1,2)*i2*p2*q2*u2 + 2*e1*g1*i2*p2*q2*u2 + 
   4*f1*m2*n1*p2*q2*u2 - 2*e2*pow(n1,2)*p2*q2*u2 - 4*e1*n1*n2*p2*q2*u2 - 
   2*d2*f1*j2*k2*r1*u2 - 2*d1*f2*j2*k2*r1*u2 + 2*c2*g1*j2*k2*r1*u2 + 
   2*c1*g2*j2*k2*r1*u2 + 2*d2*f1*i2*p2*r1*u2 + 2*d1*f2*i2*p2*r1*u2 - 
   2*c2*g1*i2*p2*r1*u2 - 2*c1*g2*i2*p2*r1*u2 + 2*g1*l2*m2*p2*r1*u2 - 
   2*f2*l2*n1*p2*r1*u2 - 2*d2*m2*n1*p2*r1*u2 - 2*f1*l2*n2*p2*r1*u2 - 
   2*d1*m2*n2*p2*r1*u2 + 4*c2*n1*n2*p2*r1*u2 + 2*c1*pow(n2,2)*p2*r1*u2 + 
   2*g1*k2*m2*q2*r1*u2 - 2*f2*k2*n1*q2*r1*u2 - 2*f1*k2*n2*q2*r1*u2 - 
   2*g2*k2*l2*pow(r1,2)*u2 + 2*d2*k2*n2*pow(r1,2)*u2 - 
   2*d1*f1*j2*k2*r2*u2 + 2*c1*g1*j2*k2*r2*u2 + 2*d1*f1*i2*p2*r2*u2 - 
   2*c1*g1*i2*p2*r2*u2 - 2*f1*l2*n1*p2*r2*u2 - 2*d1*m2*n1*p2*r2*u2 + 
   2*c2*pow(n1,2)*p2*r2*u2 + 4*c1*n1*n2*p2*r2*u2 - 2*f1*k2*n1*q2*r2*u2 - 
   4*g1*k2*l2*r1*r2*u2 + 4*d2*k2*n1*r1*r2*u2 + 4*d1*k2*n2*r1*r2*u2 + 
   2*d1*k2*n1*pow(r2,2)*u2 + 2*d2*e1*j2*k2*s1*u2 + 2*d1*e2*j2*k2*s1*u2 - 
   2*c2*f1*j2*k2*s1*u2 - 2*c1*f2*j2*k2*s1*u2 - 2*d2*e1*i2*p2*s1*u2 - 
   2*d1*e2*i2*p2*s1*u2 + 2*c2*f1*i2*p2*s1*u2 + 2*c1*f2*i2*p2*s1*u2 - 
   2*f1*l2*m2*p2*s1*u2 + 2*d1*pow(m2,2)*p2*s1*u2 + 2*e2*l2*n1*p2*s1*u2 - 
   2*c2*m2*n1*p2*s1*u2 + 2*e1*l2*n2*p2*s1*u2 - 2*c1*m2*n2*p2*s1*u2 - 
   2*f1*k2*m2*q2*s1*u2 + 2*e2*k2*n1*q2*s1*u2 + 2*e1*k2*n2*q2*s1*u2 + 
   4*f2*k2*l2*r1*s1*u2 - 2*d2*k2*m2*r1*s1*u2 - 2*c2*k2*n2*r1*s1*u2 + 
   4*f1*k2*l2*r2*s1*u2 - 2*d1*k2*m2*r2*s1*u2 - 2*c2*k2*n1*r2*s1*u2 - 
   2*c1*k2*n2*r2*s1*u2 - 2*e2*k2*l2*pow(s1,2)*u2 + 
   2*c2*k2*m2*pow(s1,2)*u2 + 2*d1*e1*j2*k2*s2*u2 - 2*c1*f1*j2*k2*s2*u2 - 
   2*d1*e1*i2*p2*s2*u2 + 2*c1*f1*i2*p2*s2*u2 + 2*e1*l2*n1*p2*s2*u2 - 
   2*c1*m2*n1*p2*s2*u2 + 2*e1*k2*n1*q2*s2*u2 + 4*f1*k2*l2*r1*s2*u2 - 
   2*d1*k2*m2*r1*s2*u2 - 2*c2*k2*n1*r1*s2*u2 - 2*c1*k2*n2*r1*s2*u2 - 
   2*c1*k2*n1*r2*s2*u2 - 4*e1*k2*l2*s1*s2*u2 + 4*c1*k2*m2*s1*s2*u2 - 
   4*f1*f2*pow(j2,2)*u1*u2 + 2*e2*g1*pow(j2,2)*u1*u2 + 
   2*e1*g2*pow(j2,2)*u1*u2 + 4*f1*f2*i2*o2*u1*u2 - 2*e2*g1*i2*o2*u1*u2 - 
   2*e1*g2*i2*o2*u1*u2 + 2*g1*pow(m2,2)*o2*u1*u2 - 4*f2*m2*n1*o2*u1*u2 - 
   4*f1*m2*n2*o2*u1*u2 + 4*e2*n1*n2*o2*u1*u2 + 2*e1*pow(n2,2)*o2*u1*u2 - 
   4*g2*j2*m2*r1*u1*u2 + 4*f2*j2*n2*r1*u1*u2 - 4*g1*j2*m2*r2*u1*u2 + 
   4*f2*j2*n1*r2*u1*u2 + 4*f1*j2*n2*r2*u1*u2 + 4*g2*i2*r1*r2*u1*u2 - 
   4*pow(n2,2)*r1*r2*u1*u2 + 2*g1*i2*pow(r2,2)*u1*u2 - 
   4*n1*n2*pow(r2,2)*u1*u2 + 4*f2*j2*m2*s1*u1*u2 - 4*e2*j2*n2*s1*u1*u2 - 
   4*f2*i2*r2*s1*u1*u2 + 4*m2*n2*r2*s1*u1*u2 + 4*f1*j2*m2*s2*u1*u2 - 
   4*e2*j2*n1*s2*u1*u2 - 4*e1*j2*n2*s2*u1*u2 - 4*f2*i2*r1*s2*u1*u2 + 
   4*m2*n2*r1*s2*u1*u2 - 4*f1*i2*r2*s2*u1*u2 + 4*m2*n1*r2*s2*u1*u2 + 
   4*e2*i2*s1*s2*u1*u2 - 4*pow(m2,2)*s1*s2*u1*u2 + 
   2*e1*i2*pow(s2,2)*u1*u2 - pow(f1,2)*pow(j2,2)*pow(u2,2) + 
   e1*g1*pow(j2,2)*pow(u2,2) + pow(f1,2)*i2*o2*pow(u2,2) - 
   e1*g1*i2*o2*pow(u2,2) - 2*f1*m2*n1*o2*pow(u2,2) + 
   e2*pow(n1,2)*o2*pow(u2,2) + 2*e1*n1*n2*o2*pow(u2,2) - 
   2*g1*j2*m2*r1*pow(u2,2) + 2*f2*j2*n1*r1*pow(u2,2) + 
   2*f1*j2*n2*r1*pow(u2,2) + g2*i2*pow(r1,2)*pow(u2,2) - 
   pow(n2,2)*pow(r1,2)*pow(u2,2) + 2*f1*j2*n1*r2*pow(u2,2) + 
   2*g1*i2*r1*r2*pow(u2,2) - 4*n1*n2*r1*r2*pow(u2,2) - 
   pow(n1,2)*pow(r2,2)*pow(u2,2) + 2*f1*j2*m2*s1*pow(u2,2) - 
   2*e2*j2*n1*s1*pow(u2,2) - 2*e1*j2*n2*s1*pow(u2,2) - 
   2*f2*i2*r1*s1*pow(u2,2) + 2*m2*n2*r1*s1*pow(u2,2) - 
   2*f1*i2*r2*s1*pow(u2,2) + 2*m2*n1*r2*s1*pow(u2,2) + 
   e2*i2*pow(s1,2)*pow(u2,2) - pow(m2,2)*pow(s1,2)*pow(u2,2) - 
   2*e1*j2*n1*s2*pow(u2,2) - 2*f1*i2*r1*s2*pow(u2,2) + 
   2*m2*n1*r1*s2*pow(u2,2) + 2*e1*i2*s1*s2*pow(u2,2) + 
   2*d2*f1*k2*l2*o2*v1 + 2*d1*f2*k2*l2*o2*v1 - 2*c2*g1*k2*l2*o2*v1 - 
   2*c1*g2*k2*l2*o2*v1 - 4*d1*d2*k2*m2*o2*v1 + 2*c2*d2*k2*n1*o2*v1 + 
   2*c2*d1*k2*n2*o2*v1 + 2*c1*d2*k2*n2*o2*v1 - 2*d2*f1*j2*l2*p2*v1 - 
   2*d1*f2*j2*l2*p2*v1 + 2*c2*g1*j2*l2*p2*v1 + 2*c1*g2*j2*l2*p2*v1 + 
   4*d1*d2*j2*m2*p2*v1 - 2*c2*d2*j2*n1*p2*v1 - 2*c2*d1*j2*n2*p2*v1 - 
   2*c1*d2*j2*n2*p2*v1 - 2*d2*f1*j2*k2*q2*v1 - 2*d1*f2*j2*k2*q2*v1 + 
   2*c2*g1*j2*k2*q2*v1 + 2*c1*g2*j2*k2*q2*v1 + 2*d2*f1*i2*p2*q2*v1 + 
   2*d1*f2*i2*p2*q2*v1 - 2*c2*g1*i2*p2*q2*v1 - 2*c1*g2*i2*p2*q2*v1 + 
   2*g1*l2*m2*p2*q2*v1 - 2*f2*l2*n1*p2*q2*v1 - 2*d2*m2*n1*p2*q2*v1 - 
   2*f1*l2*n2*p2*q2*v1 - 2*d1*m2*n2*p2*q2*v1 + 4*c2*n1*n2*p2*q2*v1 + 
   2*c1*pow(n2,2)*p2*q2*v1 - 2*g1*k2*m2*pow(q2,2)*v1 + 
   2*f2*k2*n1*pow(q2,2)*v1 + 2*f1*k2*n2*pow(q2,2)*v1 + 
   2*pow(d2,2)*j2*k2*r1*v1 - 2*pow(d2,2)*i2*p2*r1*v1 - 
   2*g2*pow(l2,2)*p2*r1*v1 + 4*d2*l2*n2*p2*r1*v1 + 2*g2*k2*l2*q2*r1*v1 - 
   2*d2*k2*n2*q2*r1*v1 + 4*d1*d2*j2*k2*r2*v1 - 4*d1*d2*i2*p2*r2*v1 - 
   2*g1*pow(l2,2)*p2*r2*v1 + 4*d2*l2*n1*p2*r2*v1 + 4*d1*l2*n2*p2*r2*v1 + 
   2*g1*k2*l2*q2*r2*v1 - 2*d2*k2*n1*q2*r2*v1 - 2*d1*k2*n2*q2*r2*v1 - 
   2*c2*d2*j2*k2*s1*v1 + 2*c2*d2*i2*p2*s1*v1 + 2*f2*pow(l2,2)*p2*s1*v1 - 
   2*d2*l2*m2*p2*s1*v1 - 2*c2*l2*n2*p2*s1*v1 - 2*f2*k2*l2*q2*s1*v1 + 
   4*d2*k2*m2*q2*s1*v1 - 2*c2*k2*n2*q2*s1*v1 - 2*d2*k2*l2*r2*s1*v1 - 
   2*c2*d1*j2*k2*s2*v1 - 2*c1*d2*j2*k2*s2*v1 + 2*c2*d1*i2*p2*s2*v1 + 
   2*c1*d2*i2*p2*s2*v1 + 2*f1*pow(l2,2)*p2*s2*v1 - 2*d1*l2*m2*p2*s2*v1 - 
   2*c2*l2*n1*p2*s2*v1 - 2*c1*l2*n2*p2*s2*v1 - 2*f1*k2*l2*q2*s2*v1 + 
   4*d1*k2*m2*q2*s2*v1 - 2*c2*k2*n1*q2*s2*v1 - 2*c1*k2*n2*q2*s2*v1 - 
   2*d2*k2*l2*r1*s2*v1 - 2*d1*k2*l2*r2*s2*v1 + 4*c2*k2*l2*s1*s2*v1 + 
   2*c1*k2*l2*pow(s2,2)*v1 + 2*d2*f2*pow(j2,2)*u1*v1 - 
   2*c2*g2*pow(j2,2)*u1*v1 - 2*d2*f2*i2*o2*u1*v1 + 2*c2*g2*i2*o2*u1*v1 - 
   2*g2*l2*m2*o2*u1*v1 + 2*f2*l2*n2*o2*u1*v1 + 2*d2*m2*n2*o2*u1*v1 - 
   2*c2*pow(n2,2)*o2*u1*v1 + 2*g2*j2*m2*q2*u1*v1 - 2*f2*j2*n2*q2*u1*v1 + 
   2*g2*j2*l2*r2*u1*v1 - 2*d2*j2*n2*r2*u1*v1 - 2*g2*i2*q2*r2*u1*v1 + 
   2*pow(n2,2)*q2*r2*u1*v1 - 2*f2*j2*l2*s2*u1*v1 - 2*d2*j2*m2*s2*u1*v1 + 
   4*c2*j2*n2*s2*u1*v1 + 2*f2*i2*q2*s2*u1*v1 - 2*m2*n2*q2*s2*u1*v1 + 
   2*d2*i2*r2*s2*u1*v1 - 2*l2*n2*r2*s2*u1*v1 - 2*c2*i2*pow(s2,2)*u1*v1 + 
   2*l2*m2*pow(s2,2)*u1*v1 + 2*d2*f1*pow(j2,2)*u2*v1 + 
   2*d1*f2*pow(j2,2)*u2*v1 - 2*c2*g1*pow(j2,2)*u2*v1 - 
   2*c1*g2*pow(j2,2)*u2*v1 - 2*d2*f1*i2*o2*u2*v1 - 2*d1*f2*i2*o2*u2*v1 + 
   2*c2*g1*i2*o2*u2*v1 + 2*c1*g2*i2*o2*u2*v1 - 2*g1*l2*m2*o2*u2*v1 + 
   2*f2*l2*n1*o2*u2*v1 + 2*d2*m2*n1*o2*u2*v1 + 2*f1*l2*n2*o2*u2*v1 + 
   2*d1*m2*n2*o2*u2*v1 - 4*c2*n1*n2*o2*u2*v1 - 2*c1*pow(n2,2)*o2*u2*v1 + 
   2*g1*j2*m2*q2*u2*v1 - 2*f2*j2*n1*q2*u2*v1 - 2*f1*j2*n2*q2*u2*v1 + 
   2*g2*j2*l2*r1*u2*v1 - 2*d2*j2*n2*r1*u2*v1 - 2*g2*i2*q2*r1*u2*v1 + 
   2*pow(n2,2)*q2*r1*u2*v1 + 2*g1*j2*l2*r2*u2*v1 - 2*d2*j2*n1*r2*u2*v1 - 
   2*d1*j2*n2*r2*u2*v1 - 2*g1*i2*q2*r2*u2*v1 + 4*n1*n2*q2*r2*u2*v1 - 
   2*f2*j2*l2*s1*u2*v1 - 2*d2*j2*m2*s1*u2*v1 + 4*c2*j2*n2*s1*u2*v1 + 
   2*f2*i2*q2*s1*u2*v1 - 2*m2*n2*q2*s1*u2*v1 + 2*d2*i2*r2*s1*u2*v1 - 
   2*l2*n2*r2*s1*u2*v1 - 2*f1*j2*l2*s2*u2*v1 - 2*d1*j2*m2*s2*u2*v1 + 
   4*c2*j2*n1*s2*u2*v1 + 4*c1*j2*n2*s2*u2*v1 + 2*f1*i2*q2*s2*u2*v1 - 
   2*m2*n1*q2*s2*u2*v1 + 2*d2*i2*r1*s2*u2*v1 - 2*l2*n2*r1*s2*u2*v1 + 
   2*d1*i2*r2*s2*u2*v1 - 2*l2*n1*r2*s2*u2*v1 - 4*c2*i2*s1*s2*u2*v1 + 
   4*l2*m2*s1*s2*u2*v1 - 2*c1*i2*pow(s2,2)*u2*v1 - 
   pow(d2,2)*pow(j2,2)*pow(v1,2) + pow(d2,2)*i2*o2*pow(v1,2) + 
   g2*pow(l2,2)*o2*pow(v1,2) - 2*d2*l2*n2*o2*pow(v1,2) - 
   2*g2*j2*l2*q2*pow(v1,2) + 2*d2*j2*n2*q2*pow(v1,2) + 
   g2*i2*pow(q2,2)*pow(v1,2) - pow(n2,2)*pow(q2,2)*pow(v1,2) + 
   2*d2*j2*l2*s2*pow(v1,2) - 2*d2*i2*q2*s2*pow(v1,2) + 
   2*l2*n2*q2*s2*pow(v1,2) - pow(l2,2)*pow(s2,2)*pow(v1,2) + 
   2*d1*f1*k2*l2*o2*v2 - 2*c1*g1*k2*l2*o2*v2 - 2*pow(d1,2)*k2*m2*o2*v2 + 
   2*c2*d1*k2*n1*o2*v2 + 2*c1*d2*k2*n1*o2*v2 + 2*c1*d1*k2*n2*o2*v2 - 
   2*d1*f1*j2*l2*p2*v2 + 2*c1*g1*j2*l2*p2*v2 + 2*pow(d1,2)*j2*m2*p2*v2 - 
   2*c2*d1*j2*n1*p2*v2 - 2*c1*d2*j2*n1*p2*v2 - 2*c1*d1*j2*n2*p2*v2 - 
   2*d1*f1*j2*k2*q2*v2 + 2*c1*g1*j2*k2*q2*v2 + 2*d1*f1*i2*p2*q2*v2 - 
   2*c1*g1*i2*p2*q2*v2 - 2*f1*l2*n1*p2*q2*v2 - 2*d1*m2*n1*p2*q2*v2 + 
   2*c2*pow(n1,2)*p2*q2*v2 + 4*c1*n1*n2*p2*q2*v2 + 
   2*f1*k2*n1*pow(q2,2)*v2 + 4*d1*d2*j2*k2*r1*v2 - 4*d1*d2*i2*p2*r1*v2 - 
   2*g1*pow(l2,2)*p2*r1*v2 + 4*d2*l2*n1*p2*r1*v2 + 4*d1*l2*n2*p2*r1*v2 + 
   2*g1*k2*l2*q2*r1*v2 - 2*d2*k2*n1*q2*r1*v2 - 2*d1*k2*n2*q2*r1*v2 + 
   2*pow(d1,2)*j2*k2*r2*v2 - 2*pow(d1,2)*i2*p2*r2*v2 + 
   4*d1*l2*n1*p2*r2*v2 - 2*d1*k2*n1*q2*r2*v2 - 2*c2*d1*j2*k2*s1*v2 - 
   2*c1*d2*j2*k2*s1*v2 + 2*c2*d1*i2*p2*s1*v2 + 2*c1*d2*i2*p2*s1*v2 + 
   2*f1*pow(l2,2)*p2*s1*v2 - 2*d1*l2*m2*p2*s1*v2 - 2*c2*l2*n1*p2*s1*v2 - 
   2*c1*l2*n2*p2*s1*v2 - 2*f1*k2*l2*q2*s1*v2 + 4*d1*k2*m2*q2*s1*v2 - 
   2*c2*k2*n1*q2*s1*v2 - 2*c1*k2*n2*q2*s1*v2 - 2*d2*k2*l2*r1*s1*v2 - 
   2*d1*k2*l2*r2*s1*v2 + 2*c2*k2*l2*pow(s1,2)*v2 - 2*c1*d1*j2*k2*s2*v2 + 
   2*c1*d1*i2*p2*s2*v2 - 2*c1*l2*n1*p2*s2*v2 - 2*c1*k2*n1*q2*s2*v2 - 
   2*d1*k2*l2*r1*s2*v2 + 4*c1*k2*l2*s1*s2*v2 + 2*d2*f1*pow(j2,2)*u1*v2 + 
   2*d1*f2*pow(j2,2)*u1*v2 - 2*c2*g1*pow(j2,2)*u1*v2 - 
   2*c1*g2*pow(j2,2)*u1*v2 - 2*d2*f1*i2*o2*u1*v2 - 2*d1*f2*i2*o2*u1*v2 + 
   2*c2*g1*i2*o2*u1*v2 + 2*c1*g2*i2*o2*u1*v2 - 2*g1*l2*m2*o2*u1*v2 + 
   2*f2*l2*n1*o2*u1*v2 + 2*d2*m2*n1*o2*u1*v2 + 2*f1*l2*n2*o2*u1*v2 + 
   2*d1*m2*n2*o2*u1*v2 - 4*c2*n1*n2*o2*u1*v2 - 2*c1*pow(n2,2)*o2*u1*v2 + 
   2*g1*j2*m2*q2*u1*v2 - 2*f2*j2*n1*q2*u1*v2 - 2*f1*j2*n2*q2*u1*v2 + 
   2*g2*j2*l2*r1*u1*v2 - 2*d2*j2*n2*r1*u1*v2 - 2*g2*i2*q2*r1*u1*v2 + 
   2*pow(n2,2)*q2*r1*u1*v2 + 2*g1*j2*l2*r2*u1*v2 - 2*d2*j2*n1*r2*u1*v2 - 
   2*d1*j2*n2*r2*u1*v2 - 2*g1*i2*q2*r2*u1*v2 + 4*n1*n2*q2*r2*u1*v2 - 
   2*f2*j2*l2*s1*u1*v2 - 2*d2*j2*m2*s1*u1*v2 + 4*c2*j2*n2*s1*u1*v2 + 
   2*f2*i2*q2*s1*u1*v2 - 2*m2*n2*q2*s1*u1*v2 + 2*d2*i2*r2*s1*u1*v2 - 
   2*l2*n2*r2*s1*u1*v2 - 2*f1*j2*l2*s2*u1*v2 - 2*d1*j2*m2*s2*u1*v2 + 
   4*c2*j2*n1*s2*u1*v2 + 4*c1*j2*n2*s2*u1*v2 + 2*f1*i2*q2*s2*u1*v2 - 
   2*m2*n1*q2*s2*u1*v2 + 2*d2*i2*r1*s2*u1*v2 - 2*l2*n2*r1*s2*u1*v2 + 
   2*d1*i2*r2*s2*u1*v2 - 2*l2*n1*r2*s2*u1*v2 - 4*c2*i2*s1*s2*u1*v2 + 
   4*l2*m2*s1*s2*u1*v2 - 2*c1*i2*pow(s2,2)*u1*v2 + 
   2*d1*f1*pow(j2,2)*u2*v2 - 2*c1*g1*pow(j2,2)*u2*v2 - 
   2*d1*f1*i2*o2*u2*v2 + 2*c1*g1*i2*o2*u2*v2 + 2*f1*l2*n1*o2*u2*v2 + 
   2*d1*m2*n1*o2*u2*v2 - 2*c2*pow(n1,2)*o2*u2*v2 - 4*c1*n1*n2*o2*u2*v2 - 
   2*f1*j2*n1*q2*u2*v2 + 2*g1*j2*l2*r1*u2*v2 - 2*d2*j2*n1*r1*u2*v2 - 
   2*d1*j2*n2*r1*u2*v2 - 2*g1*i2*q2*r1*u2*v2 + 4*n1*n2*q2*r1*u2*v2 - 
   2*d1*j2*n1*r2*u2*v2 + 2*pow(n1,2)*q2*r2*u2*v2 - 2*f1*j2*l2*s1*u2*v2 - 
   2*d1*j2*m2*s1*u2*v2 + 4*c2*j2*n1*s1*u2*v2 + 4*c1*j2*n2*s1*u2*v2 + 
   2*f1*i2*q2*s1*u2*v2 - 2*m2*n1*q2*s1*u2*v2 + 2*d2*i2*r1*s1*u2*v2 - 
   2*l2*n2*r1*s1*u2*v2 + 2*d1*i2*r2*s1*u2*v2 - 2*l2*n1*r2*s1*u2*v2 - 
   2*c2*i2*pow(s1,2)*u2*v2 + 2*l2*m2*pow(s1,2)*u2*v2 + 
   4*c1*j2*n1*s2*u2*v2 + 2*d1*i2*r1*s2*u2*v2 - 2*l2*n1*r1*s2*u2*v2 - 
   4*c1*i2*s1*s2*u2*v2 - 4*d1*d2*pow(j2,2)*v1*v2 + 4*d1*d2*i2*o2*v1*v2 + 
   2*g1*pow(l2,2)*o2*v1*v2 - 4*d2*l2*n1*o2*v1*v2 - 4*d1*l2*n2*o2*v1*v2 - 
   4*g1*j2*l2*q2*v1*v2 + 4*d2*j2*n1*q2*v1*v2 + 4*d1*j2*n2*q2*v1*v2 + 
   2*g1*i2*pow(q2,2)*v1*v2 - 4*n1*n2*pow(q2,2)*v1*v2 + 
   4*d2*j2*l2*s1*v1*v2 - 4*d2*i2*q2*s1*v1*v2 + 4*l2*n2*q2*s1*v1*v2 + 
   4*d1*j2*l2*s2*v1*v2 - 4*d1*i2*q2*s2*v1*v2 + 4*l2*n1*q2*s2*v1*v2 - 
   4*pow(l2,2)*s1*s2*v1*v2 - pow(d1,2)*pow(j2,2)*pow(v2,2) + 
   pow(d1,2)*i2*o2*pow(v2,2) - 2*d1*l2*n1*o2*pow(v2,2) + 
   2*d1*j2*n1*q2*pow(v2,2) - pow(n1,2)*pow(q2,2)*pow(v2,2) + 
   2*d1*j2*l2*s1*pow(v2,2) - 2*d1*i2*q2*s1*pow(v2,2) + 
   2*l2*n1*q2*s1*pow(v2,2) - pow(l2,2)*pow(s1,2)*pow(v2,2) - 
   2*d2*e1*k2*l2*o2*w1 - 2*d1*e2*k2*l2*o2*w1 + 2*c2*f1*k2*l2*o2*w1 + 
   2*c1*f2*k2*l2*o2*w1 + 2*c2*d1*k2*m2*o2*w1 + 2*c1*d2*k2*m2*o2*w1 - 
   2*pow(c2,2)*k2*n1*o2*w1 - 4*c1*c2*k2*n2*o2*w1 + 2*d2*e1*j2*l2*p2*w1 + 
   2*d1*e2*j2*l2*p2*w1 - 2*c2*f1*j2*l2*p2*w1 - 2*c1*f2*j2*l2*p2*w1 - 
   2*c2*d1*j2*m2*p2*w1 - 2*c1*d2*j2*m2*p2*w1 + 2*pow(c2,2)*j2*n1*p2*w1 + 
   4*c1*c2*j2*n2*p2*w1 + 2*d2*e1*j2*k2*q2*w1 + 2*d1*e2*j2*k2*q2*w1 - 
   2*c2*f1*j2*k2*q2*w1 - 2*c1*f2*j2*k2*q2*w1 - 2*d2*e1*i2*p2*q2*w1 - 
   2*d1*e2*i2*p2*q2*w1 + 2*c2*f1*i2*p2*q2*w1 + 2*c1*f2*i2*p2*q2*w1 - 
   2*f1*l2*m2*p2*q2*w1 + 2*d1*pow(m2,2)*p2*q2*w1 + 2*e2*l2*n1*p2*q2*w1 - 
   2*c2*m2*n1*p2*q2*w1 + 2*e1*l2*n2*p2*q2*w1 - 2*c1*m2*n2*p2*q2*w1 + 
   2*f1*k2*m2*pow(q2,2)*w1 - 2*e2*k2*n1*pow(q2,2)*w1 - 
   2*e1*k2*n2*pow(q2,2)*w1 - 2*c2*d2*j2*k2*r1*w1 + 2*c2*d2*i2*p2*r1*w1 + 
   2*f2*pow(l2,2)*p2*r1*w1 - 2*d2*l2*m2*p2*r1*w1 - 2*c2*l2*n2*p2*r1*w1 - 
   2*f2*k2*l2*q2*r1*w1 - 2*d2*k2*m2*q2*r1*w1 + 4*c2*k2*n2*q2*r1*w1 - 
   2*c2*d1*j2*k2*r2*w1 - 2*c1*d2*j2*k2*r2*w1 + 2*c2*d1*i2*p2*r2*w1 + 
   2*c1*d2*i2*p2*r2*w1 + 2*f1*pow(l2,2)*p2*r2*w1 - 2*d1*l2*m2*p2*r2*w1 - 
   2*c2*l2*n1*p2*r2*w1 - 2*c1*l2*n2*p2*r2*w1 - 2*f1*k2*l2*q2*r2*w1 - 
   2*d1*k2*m2*q2*r2*w1 + 4*c2*k2*n1*q2*r2*w1 + 4*c1*k2*n2*q2*r2*w1 + 
   4*d2*k2*l2*r1*r2*w1 + 2*d1*k2*l2*pow(r2,2)*w1 + 
   2*pow(c2,2)*j2*k2*s1*w1 - 2*pow(c2,2)*i2*p2*s1*w1 - 
   2*e2*pow(l2,2)*p2*s1*w1 + 4*c2*l2*m2*p2*s1*w1 + 2*e2*k2*l2*q2*s1*w1 - 
   2*c2*k2*m2*q2*s1*w1 - 2*c2*k2*l2*r2*s1*w1 + 4*c1*c2*j2*k2*s2*w1 - 
   4*c1*c2*i2*p2*s2*w1 - 2*e1*pow(l2,2)*p2*s2*w1 + 4*c1*l2*m2*p2*s2*w1 + 
   2*e1*k2*l2*q2*s2*w1 - 2*c1*k2*m2*q2*s2*w1 - 2*c2*k2*l2*r1*s2*w1 - 
   2*c1*k2*l2*r2*s2*w1 - 2*d2*e2*pow(j2,2)*u1*w1 + 
   2*c2*f2*pow(j2,2)*u1*w1 + 2*d2*e2*i2*o2*u1*w1 - 2*c2*f2*i2*o2*u1*w1 + 
   2*f2*l2*m2*o2*u1*w1 - 2*d2*pow(m2,2)*o2*u1*w1 - 2*e2*l2*n2*o2*u1*w1 + 
   2*c2*m2*n2*o2*u1*w1 - 2*f2*j2*m2*q2*u1*w1 + 2*e2*j2*n2*q2*u1*w1 - 
   2*f2*j2*l2*r2*u1*w1 + 4*d2*j2*m2*r2*u1*w1 - 2*c2*j2*n2*r2*u1*w1 + 
   2*f2*i2*q2*r2*u1*w1 - 2*m2*n2*q2*r2*u1*w1 - 2*d2*i2*pow(r2,2)*u1*w1 + 
   2*l2*n2*pow(r2,2)*u1*w1 + 2*e2*j2*l2*s2*u1*w1 - 2*c2*j2*m2*s2*u1*w1 - 
   2*e2*i2*q2*s2*u1*w1 + 2*pow(m2,2)*q2*s2*u1*w1 + 2*c2*i2*r2*s2*u1*w1 - 
   2*l2*m2*r2*s2*u1*w1 - 2*d2*e1*pow(j2,2)*u2*w1 - 
   2*d1*e2*pow(j2,2)*u2*w1 + 2*c2*f1*pow(j2,2)*u2*w1 + 
   2*c1*f2*pow(j2,2)*u2*w1 + 2*d2*e1*i2*o2*u2*w1 + 2*d1*e2*i2*o2*u2*w1 - 
   2*c2*f1*i2*o2*u2*w1 - 2*c1*f2*i2*o2*u2*w1 + 2*f1*l2*m2*o2*u2*w1 - 
   2*d1*pow(m2,2)*o2*u2*w1 - 2*e2*l2*n1*o2*u2*w1 + 2*c2*m2*n1*o2*u2*w1 - 
   2*e1*l2*n2*o2*u2*w1 + 2*c1*m2*n2*o2*u2*w1 - 2*f1*j2*m2*q2*u2*w1 + 
   2*e2*j2*n1*q2*u2*w1 + 2*e1*j2*n2*q2*u2*w1 - 2*f2*j2*l2*r1*u2*w1 + 
   4*d2*j2*m2*r1*u2*w1 - 2*c2*j2*n2*r1*u2*w1 + 2*f2*i2*q2*r1*u2*w1 - 
   2*m2*n2*q2*r1*u2*w1 - 2*f1*j2*l2*r2*u2*w1 + 4*d1*j2*m2*r2*u2*w1 - 
   2*c2*j2*n1*r2*u2*w1 - 2*c1*j2*n2*r2*u2*w1 + 2*f1*i2*q2*r2*u2*w1 - 
   2*m2*n1*q2*r2*u2*w1 - 4*d2*i2*r1*r2*u2*w1 + 4*l2*n2*r1*r2*u2*w1 - 
   2*d1*i2*pow(r2,2)*u2*w1 + 2*l2*n1*pow(r2,2)*u2*w1 + 
   2*e2*j2*l2*s1*u2*w1 - 2*c2*j2*m2*s1*u2*w1 - 2*e2*i2*q2*s1*u2*w1 + 
   2*pow(m2,2)*q2*s1*u2*w1 + 2*c2*i2*r2*s1*u2*w1 - 2*l2*m2*r2*s1*u2*w1 + 
   2*e1*j2*l2*s2*u2*w1 - 2*c1*j2*m2*s2*u2*w1 - 2*e1*i2*q2*s2*u2*w1 + 
   2*c2*i2*r1*s2*u2*w1 - 2*l2*m2*r1*s2*u2*w1 + 2*c1*i2*r2*s2*u2*w1 + 
   2*c2*d2*pow(j2,2)*v1*w1 - 2*c2*d2*i2*o2*v1*w1 - 
   2*f2*pow(l2,2)*o2*v1*w1 + 2*d2*l2*m2*o2*v1*w1 + 2*c2*l2*n2*o2*v1*w1 + 
   4*f2*j2*l2*q2*v1*w1 - 2*d2*j2*m2*q2*v1*w1 - 2*c2*j2*n2*q2*v1*w1 - 
   2*f2*i2*pow(q2,2)*v1*w1 + 2*m2*n2*pow(q2,2)*v1*w1 - 
   2*d2*j2*l2*r2*v1*w1 + 2*d2*i2*q2*r2*v1*w1 - 2*l2*n2*q2*r2*v1*w1 - 
   2*c2*j2*l2*s2*v1*w1 + 2*c2*i2*q2*s2*v1*w1 - 2*l2*m2*q2*s2*v1*w1 + 
   2*pow(l2,2)*r2*s2*v1*w1 + 2*c2*d1*pow(j2,2)*v2*w1 + 
   2*c1*d2*pow(j2,2)*v2*w1 - 2*c2*d1*i2*o2*v2*w1 - 2*c1*d2*i2*o2*v2*w1 - 
   2*f1*pow(l2,2)*o2*v2*w1 + 2*d1*l2*m2*o2*v2*w1 + 2*c2*l2*n1*o2*v2*w1 + 
   2*c1*l2*n2*o2*v2*w1 + 4*f1*j2*l2*q2*v2*w1 - 2*d1*j2*m2*q2*v2*w1 - 
   2*c2*j2*n1*q2*v2*w1 - 2*c1*j2*n2*q2*v2*w1 - 2*f1*i2*pow(q2,2)*v2*w1 + 
   2*m2*n1*pow(q2,2)*v2*w1 - 2*d2*j2*l2*r1*v2*w1 + 2*d2*i2*q2*r1*v2*w1 - 
   2*l2*n2*q2*r1*v2*w1 - 2*d1*j2*l2*r2*v2*w1 + 2*d1*i2*q2*r2*v2*w1 - 
   2*l2*n1*q2*r2*v2*w1 - 2*c2*j2*l2*s1*v2*w1 + 2*c2*i2*q2*s1*v2*w1 - 
   2*l2*m2*q2*s1*v2*w1 + 2*pow(l2,2)*r2*s1*v2*w1 - 2*c1*j2*l2*s2*v2*w1 + 
   2*c1*i2*q2*s2*v2*w1 + 2*pow(l2,2)*r1*s2*v2*w1 - 
   pow(c2,2)*pow(j2,2)*pow(w1,2) + pow(c2,2)*i2*o2*pow(w1,2) + 
   e2*pow(l2,2)*o2*pow(w1,2) - 2*c2*l2*m2*o2*pow(w1,2) - 
   2*e2*j2*l2*q2*pow(w1,2) + 2*c2*j2*m2*q2*pow(w1,2) + 
   e2*i2*pow(q2,2)*pow(w1,2) - pow(m2,2)*pow(q2,2)*pow(w1,2) + 
   2*c2*j2*l2*r2*pow(w1,2) - 2*c2*i2*q2*r2*pow(w1,2) + 
   2*l2*m2*q2*r2*pow(w1,2) - pow(l2,2)*pow(r2,2)*pow(w1,2) - 
   2*d1*e1*k2*l2*o2*w2 + 2*c1*f1*k2*l2*o2*w2 + 2*c1*d1*k2*m2*o2*w2 - 
   4*c1*c2*k2*n1*o2*w2 - 2*pow(c1,2)*k2*n2*o2*w2 + 2*d1*e1*j2*l2*p2*w2 - 
   2*c1*f1*j2*l2*p2*w2 - 2*c1*d1*j2*m2*p2*w2 + 4*c1*c2*j2*n1*p2*w2 + 
   2*pow(c1,2)*j2*n2*p2*w2 + 2*d1*e1*j2*k2*q2*w2 - 2*c1*f1*j2*k2*q2*w2 - 
   2*d1*e1*i2*p2*q2*w2 + 2*c1*f1*i2*p2*q2*w2 + 2*e1*l2*n1*p2*q2*w2 - 
   2*c1*m2*n1*p2*q2*w2 - 2*e1*k2*n1*pow(q2,2)*w2 - 2*c2*d1*j2*k2*r1*w2 - 
   2*c1*d2*j2*k2*r1*w2 + 2*c2*d1*i2*p2*r1*w2 + 2*c1*d2*i2*p2*r1*w2 + 
   2*f1*pow(l2,2)*p2*r1*w2 - 2*d1*l2*m2*p2*r1*w2 - 2*c2*l2*n1*p2*r1*w2 - 
   2*c1*l2*n2*p2*r1*w2 - 2*f1*k2*l2*q2*r1*w2 - 2*d1*k2*m2*q2*r1*w2 + 
   4*c2*k2*n1*q2*r1*w2 + 4*c1*k2*n2*q2*r1*w2 + 2*d2*k2*l2*pow(r1,2)*w2 - 
   2*c1*d1*j2*k2*r2*w2 + 2*c1*d1*i2*p2*r2*w2 - 2*c1*l2*n1*p2*r2*w2 + 
   4*c1*k2*n1*q2*r2*w2 + 4*d1*k2*l2*r1*r2*w2 + 4*c1*c2*j2*k2*s1*w2 - 
   4*c1*c2*i2*p2*s1*w2 - 2*e1*pow(l2,2)*p2*s1*w2 + 4*c1*l2*m2*p2*s1*w2 + 
   2*e1*k2*l2*q2*s1*w2 - 2*c1*k2*m2*q2*s1*w2 - 2*c2*k2*l2*r1*s1*w2 - 
   2*c1*k2*l2*r2*s1*w2 + 2*pow(c1,2)*j2*k2*s2*w2 - 
   2*pow(c1,2)*i2*p2*s2*w2 - 2*c1*k2*l2*r1*s2*w2 - 
   2*d2*e1*pow(j2,2)*u1*w2 - 2*d1*e2*pow(j2,2)*u1*w2 + 
   2*c2*f1*pow(j2,2)*u1*w2 + 2*c1*f2*pow(j2,2)*u1*w2 + 
   2*d2*e1*i2*o2*u1*w2 + 2*d1*e2*i2*o2*u1*w2 - 2*c2*f1*i2*o2*u1*w2 - 
   2*c1*f2*i2*o2*u1*w2 + 2*f1*l2*m2*o2*u1*w2 - 2*d1*pow(m2,2)*o2*u1*w2 - 
   2*e2*l2*n1*o2*u1*w2 + 2*c2*m2*n1*o2*u1*w2 - 2*e1*l2*n2*o2*u1*w2 + 
   2*c1*m2*n2*o2*u1*w2 - 2*f1*j2*m2*q2*u1*w2 + 2*e2*j2*n1*q2*u1*w2 + 
   2*e1*j2*n2*q2*u1*w2 - 2*f2*j2*l2*r1*u1*w2 + 4*d2*j2*m2*r1*u1*w2 - 
   2*c2*j2*n2*r1*u1*w2 + 2*f2*i2*q2*r1*u1*w2 - 2*m2*n2*q2*r1*u1*w2 - 
   2*f1*j2*l2*r2*u1*w2 + 4*d1*j2*m2*r2*u1*w2 - 2*c2*j2*n1*r2*u1*w2 - 
   2*c1*j2*n2*r2*u1*w2 + 2*f1*i2*q2*r2*u1*w2 - 2*m2*n1*q2*r2*u1*w2 - 
   4*d2*i2*r1*r2*u1*w2 + 4*l2*n2*r1*r2*u1*w2 - 2*d1*i2*pow(r2,2)*u1*w2 + 
   2*l2*n1*pow(r2,2)*u1*w2 + 2*e2*j2*l2*s1*u1*w2 - 2*c2*j2*m2*s1*u1*w2 - 
   2*e2*i2*q2*s1*u1*w2 + 2*pow(m2,2)*q2*s1*u1*w2 + 2*c2*i2*r2*s1*u1*w2 - 
   2*l2*m2*r2*s1*u1*w2 + 2*e1*j2*l2*s2*u1*w2 - 2*c1*j2*m2*s2*u1*w2 - 
   2*e1*i2*q2*s2*u1*w2 + 2*c2*i2*r1*s2*u1*w2 - 2*l2*m2*r1*s2*u1*w2 + 
   2*c1*i2*r2*s2*u1*w2 - 2*d1*e1*pow(j2,2)*u2*w2 + 
   2*c1*f1*pow(j2,2)*u2*w2 + 2*d1*e1*i2*o2*u2*w2 - 2*c1*f1*i2*o2*u2*w2 - 
   2*e1*l2*n1*o2*u2*w2 + 2*c1*m2*n1*o2*u2*w2 + 2*e1*j2*n1*q2*u2*w2 - 
   2*f1*j2*l2*r1*u2*w2 + 4*d1*j2*m2*r1*u2*w2 - 2*c2*j2*n1*r1*u2*w2 - 
   2*c1*j2*n2*r1*u2*w2 + 2*f1*i2*q2*r1*u2*w2 - 2*m2*n1*q2*r1*u2*w2 - 
   2*d2*i2*pow(r1,2)*u2*w2 + 2*l2*n2*pow(r1,2)*u2*w2 - 
   2*c1*j2*n1*r2*u2*w2 - 4*d1*i2*r1*r2*u2*w2 + 4*l2*n1*r1*r2*u2*w2 + 
   2*e1*j2*l2*s1*u2*w2 - 2*c1*j2*m2*s1*u2*w2 - 2*e1*i2*q2*s1*u2*w2 + 
   2*c2*i2*r1*s1*u2*w2 - 2*l2*m2*r1*s1*u2*w2 + 2*c1*i2*r2*s1*u2*w2 + 
   2*c1*i2*r1*s2*u2*w2 + 2*c2*d1*pow(j2,2)*v1*w2 + 
   2*c1*d2*pow(j2,2)*v1*w2 - 2*c2*d1*i2*o2*v1*w2 - 2*c1*d2*i2*o2*v1*w2 - 
   2*f1*pow(l2,2)*o2*v1*w2 + 2*d1*l2*m2*o2*v1*w2 + 2*c2*l2*n1*o2*v1*w2 + 
   2*c1*l2*n2*o2*v1*w2 + 4*f1*j2*l2*q2*v1*w2 - 2*d1*j2*m2*q2*v1*w2 - 
   2*c2*j2*n1*q2*v1*w2 - 2*c1*j2*n2*q2*v1*w2 - 2*f1*i2*pow(q2,2)*v1*w2 + 
   2*m2*n1*pow(q2,2)*v1*w2 - 2*d2*j2*l2*r1*v1*w2 + 2*d2*i2*q2*r1*v1*w2 - 
   2*l2*n2*q2*r1*v1*w2 - 2*d1*j2*l2*r2*v1*w2 + 2*d1*i2*q2*r2*v1*w2 - 
   2*l2*n1*q2*r2*v1*w2 - 2*c2*j2*l2*s1*v1*w2 + 2*c2*i2*q2*s1*v1*w2 - 
   2*l2*m2*q2*s1*v1*w2 + 2*pow(l2,2)*r2*s1*v1*w2 - 2*c1*j2*l2*s2*v1*w2 + 
   2*c1*i2*q2*s2*v1*w2 + 2*pow(l2,2)*r1*s2*v1*w2 + 
   2*c1*d1*pow(j2,2)*v2*w2 - 2*c1*d1*i2*o2*v2*w2 + 2*c1*l2*n1*o2*v2*w2 - 
   2*c1*j2*n1*q2*v2*w2 - 2*d1*j2*l2*r1*v2*w2 + 2*d1*i2*q2*r1*v2*w2 - 
   2*l2*n1*q2*r1*v2*w2 - 2*c1*j2*l2*s1*v2*w2 + 2*c1*i2*q2*s1*v2*w2 + 
   2*pow(l2,2)*r1*s1*v2*w2 - 4*c1*c2*pow(j2,2)*w1*w2 + 
   4*c1*c2*i2*o2*w1*w2 + 2*e1*pow(l2,2)*o2*w1*w2 - 4*c1*l2*m2*o2*w1*w2 - 
   4*e1*j2*l2*q2*w1*w2 + 4*c1*j2*m2*q2*w1*w2 + 2*e1*i2*pow(q2,2)*w1*w2 + 
   4*c2*j2*l2*r1*w1*w2 - 4*c2*i2*q2*r1*w1*w2 + 4*l2*m2*q2*r1*w1*w2 + 
   4*c1*j2*l2*r2*w1*w2 - 4*c1*i2*q2*r2*w1*w2 - 4*pow(l2,2)*r1*r2*w1*w2 - 
   pow(c1,2)*pow(j2,2)*pow(w2,2) + pow(c1,2)*i2*o2*pow(w2,2) + 
   2*c1*j2*l2*r1*pow(w2,2) - 2*c1*i2*q2*r1*pow(w2,2) - 
   pow(l2,2)*pow(r1,2)*pow(w2,2) + 2*f1*f2*pow(k2,2)*o2*z1 - 
   e2*g1*pow(k2,2)*o2*z1 - e1*g2*pow(k2,2)*o2*z1 - 4*f1*f2*j2*k2*p2*z1 + 
   2*e2*g1*j2*k2*p2*z1 + 2*e1*g2*j2*k2*p2*z1 + 2*f1*f2*i2*pow(p2,2)*z1 - 
   e2*g1*i2*pow(p2,2)*z1 - e1*g2*i2*pow(p2,2)*z1 + 
   g1*pow(m2,2)*pow(p2,2)*z1 - 2*f2*m2*n1*pow(p2,2)*z1 - 
   2*f1*m2*n2*pow(p2,2)*z1 + 2*e2*n1*n2*pow(p2,2)*z1 + 
   e1*pow(n2,2)*pow(p2,2)*z1 - 2*g2*k2*m2*p2*r1*z1 + 
   2*f2*k2*n2*p2*r1*z1 - 2*g1*k2*m2*p2*r2*z1 + 2*f2*k2*n1*p2*r2*z1 + 
   2*f1*k2*n2*p2*r2*z1 + 2*g2*pow(k2,2)*r1*r2*z1 + 
   g1*pow(k2,2)*pow(r2,2)*z1 + 2*f2*k2*m2*p2*s1*z1 - 
   2*e2*k2*n2*p2*s1*z1 - 2*f2*pow(k2,2)*r2*s1*z1 + 2*f1*k2*m2*p2*s2*z1 - 
   2*e2*k2*n1*p2*s2*z1 - 2*e1*k2*n2*p2*s2*z1 - 2*f2*pow(k2,2)*r1*s2*z1 - 
   2*f1*pow(k2,2)*r2*s2*z1 + 2*e2*pow(k2,2)*s1*s2*z1 + 
   e1*pow(k2,2)*pow(s2,2)*z1 + 2*f1*f2*pow(j2,2)*t2*z1 - 
   e2*g1*pow(j2,2)*t2*z1 - e1*g2*pow(j2,2)*t2*z1 - 2*f1*f2*i2*o2*t2*z1 + 
   e2*g1*i2*o2*t2*z1 + e1*g2*i2*o2*t2*z1 - g1*pow(m2,2)*o2*t2*z1 + 
   2*f2*m2*n1*o2*t2*z1 + 2*f1*m2*n2*o2*t2*z1 - 2*e2*n1*n2*o2*t2*z1 - 
   e1*pow(n2,2)*o2*t2*z1 + 2*g2*j2*m2*r1*t2*z1 - 2*f2*j2*n2*r1*t2*z1 + 
   2*g1*j2*m2*r2*t2*z1 - 2*f2*j2*n1*r2*t2*z1 - 2*f1*j2*n2*r2*t2*z1 - 
   2*g2*i2*r1*r2*t2*z1 + 2*pow(n2,2)*r1*r2*t2*z1 - 
   g1*i2*pow(r2,2)*t2*z1 + 2*n1*n2*pow(r2,2)*t2*z1 - 
   2*f2*j2*m2*s1*t2*z1 + 2*e2*j2*n2*s1*t2*z1 + 2*f2*i2*r2*s1*t2*z1 - 
   2*m2*n2*r2*s1*t2*z1 - 2*f1*j2*m2*s2*t2*z1 + 2*e2*j2*n1*s2*t2*z1 + 
   2*e1*j2*n2*s2*t2*z1 + 2*f2*i2*r1*s2*t2*z1 - 2*m2*n2*r1*s2*t2*z1 + 
   2*f1*i2*r2*s2*t2*z1 - 2*m2*n1*r2*s2*t2*z1 - 2*e2*i2*s1*s2*t2*z1 + 
   2*pow(m2,2)*s1*s2*t2*z1 - e1*i2*pow(s2,2)*t2*z1 + 
   2*g2*k2*m2*o2*v1*z1 - 2*f2*k2*n2*o2*v1*z1 - 2*g2*j2*m2*p2*v1*z1 + 
   2*f2*j2*n2*p2*v1*z1 - 2*g2*j2*k2*r2*v1*z1 + 2*g2*i2*p2*r2*v1*z1 - 
   2*pow(n2,2)*p2*r2*v1*z1 + 2*f2*j2*k2*s2*v1*z1 - 2*f2*i2*p2*s2*v1*z1 + 
   2*m2*n2*p2*s2*v1*z1 + 2*k2*n2*r2*s2*v1*z1 - 2*k2*m2*pow(s2,2)*v1*z1 + 
   2*g1*k2*m2*o2*v2*z1 - 2*f2*k2*n1*o2*v2*z1 - 2*f1*k2*n2*o2*v2*z1 - 
   2*g1*j2*m2*p2*v2*z1 + 2*f2*j2*n1*p2*v2*z1 + 2*f1*j2*n2*p2*v2*z1 - 
   2*g2*j2*k2*r1*v2*z1 + 2*g2*i2*p2*r1*v2*z1 - 2*pow(n2,2)*p2*r1*v2*z1 - 
   2*g1*j2*k2*r2*v2*z1 + 2*g1*i2*p2*r2*v2*z1 - 4*n1*n2*p2*r2*v2*z1 + 
   2*f2*j2*k2*s1*v2*z1 - 2*f2*i2*p2*s1*v2*z1 + 2*m2*n2*p2*s1*v2*z1 + 
   2*k2*n2*r2*s1*v2*z1 + 2*f1*j2*k2*s2*v2*z1 - 2*f1*i2*p2*s2*v2*z1 + 
   2*m2*n1*p2*s2*v2*z1 + 2*k2*n2*r1*s2*v2*z1 + 2*k2*n1*r2*s2*v2*z1 - 
   4*k2*m2*s1*s2*v2*z1 + 2*g2*pow(j2,2)*v1*v2*z1 - 2*g2*i2*o2*v1*v2*z1 + 
   2*pow(n2,2)*o2*v1*v2*z1 - 4*j2*n2*s2*v1*v2*z1 + 
   2*i2*pow(s2,2)*v1*v2*z1 + g1*pow(j2,2)*pow(v2,2)*z1 - 
   g1*i2*o2*pow(v2,2)*z1 + 2*n1*n2*o2*pow(v2,2)*z1 - 
   2*j2*n2*s1*pow(v2,2)*z1 - 2*j2*n1*s2*pow(v2,2)*z1 + 
   2*i2*s1*s2*pow(v2,2)*z1 - 2*f2*k2*m2*o2*w1*z1 + 2*e2*k2*n2*o2*w1*z1 + 
   2*f2*j2*m2*p2*w1*z1 - 2*e2*j2*n2*p2*w1*z1 + 2*f2*j2*k2*r2*w1*z1 - 
   2*f2*i2*p2*r2*w1*z1 + 2*m2*n2*p2*r2*w1*z1 - 2*k2*n2*pow(r2,2)*w1*z1 - 
   2*e2*j2*k2*s2*w1*z1 + 2*e2*i2*p2*s2*w1*z1 - 2*pow(m2,2)*p2*s2*w1*z1 + 
   2*k2*m2*r2*s2*w1*z1 - 2*f2*pow(j2,2)*v2*w1*z1 + 2*f2*i2*o2*v2*w1*z1 - 
   2*m2*n2*o2*v2*w1*z1 + 2*j2*n2*r2*v2*w1*z1 + 2*j2*m2*s2*v2*w1*z1 - 
   2*i2*r2*s2*v2*w1*z1 - 2*f1*k2*m2*o2*w2*z1 + 2*e2*k2*n1*o2*w2*z1 + 
   2*e1*k2*n2*o2*w2*z1 + 2*f1*j2*m2*p2*w2*z1 - 2*e2*j2*n1*p2*w2*z1 - 
   2*e1*j2*n2*p2*w2*z1 + 2*f2*j2*k2*r1*w2*z1 - 2*f2*i2*p2*r1*w2*z1 + 
   2*m2*n2*p2*r1*w2*z1 + 2*f1*j2*k2*r2*w2*z1 - 2*f1*i2*p2*r2*w2*z1 + 
   2*m2*n1*p2*r2*w2*z1 - 4*k2*n2*r1*r2*w2*z1 - 2*k2*n1*pow(r2,2)*w2*z1 - 
   2*e2*j2*k2*s1*w2*z1 + 2*e2*i2*p2*s1*w2*z1 - 2*pow(m2,2)*p2*s1*w2*z1 + 
   2*k2*m2*r2*s1*w2*z1 - 2*e1*j2*k2*s2*w2*z1 + 2*e1*i2*p2*s2*w2*z1 + 
   2*k2*m2*r1*s2*w2*z1 - 2*f2*pow(j2,2)*v1*w2*z1 + 2*f2*i2*o2*v1*w2*z1 - 
   2*m2*n2*o2*v1*w2*z1 + 2*j2*n2*r2*v1*w2*z1 + 2*j2*m2*s2*v1*w2*z1 - 
   2*i2*r2*s2*v1*w2*z1 - 2*f1*pow(j2,2)*v2*w2*z1 + 2*f1*i2*o2*v2*w2*z1 - 
   2*m2*n1*o2*v2*w2*z1 + 2*j2*n2*r1*v2*w2*z1 + 2*j2*n1*r2*v2*w2*z1 + 
   2*j2*m2*s1*v2*w2*z1 - 2*i2*r2*s1*v2*w2*z1 - 2*i2*r1*s2*v2*w2*z1 + 
   2*e2*pow(j2,2)*w1*w2*z1 - 2*e2*i2*o2*w1*w2*z1 + 
   2*pow(m2,2)*o2*w1*w2*z1 - 4*j2*m2*r2*w1*w2*z1 + 
   2*i2*pow(r2,2)*w1*w2*z1 + e1*pow(j2,2)*pow(w2,2)*z1 - 
   e1*i2*o2*pow(w2,2)*z1 - 2*j2*m2*r1*pow(w2,2)*z1 + 
   2*i2*r1*r2*pow(w2,2)*z1 + pow(f1,2)*pow(k2,2)*o2*z2 - 
   e1*g1*pow(k2,2)*o2*z2 - 2*pow(f1,2)*j2*k2*p2*z2 + 
   2*e1*g1*j2*k2*p2*z2 + pow(f1,2)*i2*pow(p2,2)*z2 - 
   e1*g1*i2*pow(p2,2)*z2 - 2*f1*m2*n1*pow(p2,2)*z2 + 
   e2*pow(n1,2)*pow(p2,2)*z2 + 2*e1*n1*n2*pow(p2,2)*z2 - 
   2*g1*k2*m2*p2*r1*z2 + 2*f2*k2*n1*p2*r1*z2 + 2*f1*k2*n2*p2*r1*z2 + 
   g2*pow(k2,2)*pow(r1,2)*z2 + 2*f1*k2*n1*p2*r2*z2 + 
   2*g1*pow(k2,2)*r1*r2*z2 + 2*f1*k2*m2*p2*s1*z2 - 2*e2*k2*n1*p2*s1*z2 - 
   2*e1*k2*n2*p2*s1*z2 - 2*f2*pow(k2,2)*r1*s1*z2 - 
   2*f1*pow(k2,2)*r2*s1*z2 + e2*pow(k2,2)*pow(s1,2)*z2 - 
   2*e1*k2*n1*p2*s2*z2 - 2*f1*pow(k2,2)*r1*s2*z2 + 
   2*e1*pow(k2,2)*s1*s2*z2 + pow(f1,2)*pow(j2,2)*t2*z2 - 
   e1*g1*pow(j2,2)*t2*z2 - pow(f1,2)*i2*o2*t2*z2 + e1*g1*i2*o2*t2*z2 + 
   2*f1*m2*n1*o2*t2*z2 - e2*pow(n1,2)*o2*t2*z2 - 2*e1*n1*n2*o2*t2*z2 + 
   2*g1*j2*m2*r1*t2*z2 - 2*f2*j2*n1*r1*t2*z2 - 2*f1*j2*n2*r1*t2*z2 - 
   g2*i2*pow(r1,2)*t2*z2 + pow(n2,2)*pow(r1,2)*t2*z2 - 
   2*f1*j2*n1*r2*t2*z2 - 2*g1*i2*r1*r2*t2*z2 + 4*n1*n2*r1*r2*t2*z2 + 
   pow(n1,2)*pow(r2,2)*t2*z2 - 2*f1*j2*m2*s1*t2*z2 + 
   2*e2*j2*n1*s1*t2*z2 + 2*e1*j2*n2*s1*t2*z2 + 2*f2*i2*r1*s1*t2*z2 - 
   2*m2*n2*r1*s1*t2*z2 + 2*f1*i2*r2*s1*t2*z2 - 2*m2*n1*r2*s1*t2*z2 - 
   e2*i2*pow(s1,2)*t2*z2 + pow(m2,2)*pow(s1,2)*t2*z2 + 
   2*e1*j2*n1*s2*t2*z2 + 2*f1*i2*r1*s2*t2*z2 - 2*m2*n1*r1*s2*t2*z2 - 
   2*e1*i2*s1*s2*t2*z2 + 2*g1*k2*m2*o2*v1*z2 - 2*f2*k2*n1*o2*v1*z2 - 
   2*f1*k2*n2*o2*v1*z2 - 2*g1*j2*m2*p2*v1*z2 + 2*f2*j2*n1*p2*v1*z2 + 
   2*f1*j2*n2*p2*v1*z2 - 2*g2*j2*k2*r1*v1*z2 + 2*g2*i2*p2*r1*v1*z2 - 
   2*pow(n2,2)*p2*r1*v1*z2 - 2*g1*j2*k2*r2*v1*z2 + 2*g1*i2*p2*r2*v1*z2 - 
   4*n1*n2*p2*r2*v1*z2 + 2*f2*j2*k2*s1*v1*z2 - 2*f2*i2*p2*s1*v1*z2 + 
   2*m2*n2*p2*s1*v1*z2 + 2*k2*n2*r2*s1*v1*z2 + 2*f1*j2*k2*s2*v1*z2 - 
   2*f1*i2*p2*s2*v1*z2 + 2*m2*n1*p2*s2*v1*z2 + 2*k2*n2*r1*s2*v1*z2 + 
   2*k2*n1*r2*s2*v1*z2 - 4*k2*m2*s1*s2*v1*z2 + 
   g2*pow(j2,2)*pow(v1,2)*z2 - g2*i2*o2*pow(v1,2)*z2 + 
   pow(n2,2)*o2*pow(v1,2)*z2 - 2*j2*n2*s2*pow(v1,2)*z2 + 
   i2*pow(s2,2)*pow(v1,2)*z2 - 2*f1*k2*n1*o2*v2*z2 + 
   2*f1*j2*n1*p2*v2*z2 - 2*g1*j2*k2*r1*v2*z2 + 2*g1*i2*p2*r1*v2*z2 - 
   4*n1*n2*p2*r1*v2*z2 - 2*pow(n1,2)*p2*r2*v2*z2 + 2*f1*j2*k2*s1*v2*z2 - 
   2*f1*i2*p2*s1*v2*z2 + 2*m2*n1*p2*s1*v2*z2 + 2*k2*n2*r1*s1*v2*z2 + 
   2*k2*n1*r2*s1*v2*z2 - 2*k2*m2*pow(s1,2)*v2*z2 + 2*k2*n1*r1*s2*v2*z2 + 
   2*g1*pow(j2,2)*v1*v2*z2 - 2*g1*i2*o2*v1*v2*z2 + 4*n1*n2*o2*v1*v2*z2 - 
   4*j2*n2*s1*v1*v2*z2 - 4*j2*n1*s2*v1*v2*z2 + 4*i2*s1*s2*v1*v2*z2 + 
   pow(n1,2)*o2*pow(v2,2)*z2 - 2*j2*n1*s1*pow(v2,2)*z2 + 
   i2*pow(s1,2)*pow(v2,2)*z2 - 2*f1*k2*m2*o2*w1*z2 + 
   2*e2*k2*n1*o2*w1*z2 + 2*e1*k2*n2*o2*w1*z2 + 2*f1*j2*m2*p2*w1*z2 - 
   2*e2*j2*n1*p2*w1*z2 - 2*e1*j2*n2*p2*w1*z2 + 2*f2*j2*k2*r1*w1*z2 - 
   2*f2*i2*p2*r1*w1*z2 + 2*m2*n2*p2*r1*w1*z2 + 2*f1*j2*k2*r2*w1*z2 - 
   2*f1*i2*p2*r2*w1*z2 + 2*m2*n1*p2*r2*w1*z2 - 4*k2*n2*r1*r2*w1*z2 - 
   2*k2*n1*pow(r2,2)*w1*z2 - 2*e2*j2*k2*s1*w1*z2 + 2*e2*i2*p2*s1*w1*z2 - 
   2*pow(m2,2)*p2*s1*w1*z2 + 2*k2*m2*r2*s1*w1*z2 - 2*e1*j2*k2*s2*w1*z2 + 
   2*e1*i2*p2*s2*w1*z2 + 2*k2*m2*r1*s2*w1*z2 - 2*f2*pow(j2,2)*v1*w1*z2 + 
   2*f2*i2*o2*v1*w1*z2 - 2*m2*n2*o2*v1*w1*z2 + 2*j2*n2*r2*v1*w1*z2 + 
   2*j2*m2*s2*v1*w1*z2 - 2*i2*r2*s2*v1*w1*z2 - 2*f1*pow(j2,2)*v2*w1*z2 + 
   2*f1*i2*o2*v2*w1*z2 - 2*m2*n1*o2*v2*w1*z2 + 2*j2*n2*r1*v2*w1*z2 + 
   2*j2*n1*r2*v2*w1*z2 + 2*j2*m2*s1*v2*w1*z2 - 2*i2*r2*s1*v2*w1*z2 - 
   2*i2*r1*s2*v2*w1*z2 + e2*pow(j2,2)*pow(w1,2)*z2 - 
   e2*i2*o2*pow(w1,2)*z2 + pow(m2,2)*o2*pow(w1,2)*z2 - 
   2*j2*m2*r2*pow(w1,2)*z2 + i2*pow(r2,2)*pow(w1,2)*z2 + 
   2*e1*k2*n1*o2*w2*z2 - 2*e1*j2*n1*p2*w2*z2 + 2*f1*j2*k2*r1*w2*z2 - 
   2*f1*i2*p2*r1*w2*z2 + 2*m2*n1*p2*r1*w2*z2 - 2*k2*n2*pow(r1,2)*w2*z2 - 
   4*k2*n1*r1*r2*w2*z2 - 2*e1*j2*k2*s1*w2*z2 + 2*e1*i2*p2*s1*w2*z2 + 
   2*k2*m2*r1*s1*w2*z2 - 2*f1*pow(j2,2)*v1*w2*z2 + 2*f1*i2*o2*v1*w2*z2 - 
   2*m2*n1*o2*v1*w2*z2 + 2*j2*n2*r1*v1*w2*z2 + 2*j2*n1*r2*v1*w2*z2 + 
   2*j2*m2*s1*v1*w2*z2 - 2*i2*r2*s1*v1*w2*z2 - 2*i2*r1*s2*v1*w2*z2 + 
   2*j2*n1*r1*v2*w2*z2 - 2*i2*r1*s1*v2*w2*z2 + 2*e1*pow(j2,2)*w1*w2*z2 - 
   2*e1*i2*o2*w1*w2*z2 - 4*j2*m2*r1*w1*w2*z2 + 4*i2*r1*r2*w1*w2*z2 + 
   i2*pow(r1,2)*pow(w2,2)*z2
;

   	k[6]  = pow(d0,2)*e0*pow(k2,2)*o2 - 2*c0*d0*f0*pow(k2,2)*o2 + 
   pow(c0,2)*g0*pow(k2,2)*o2 - 2*pow(d0,2)*e0*j2*k2*p2 + 
   4*c0*d0*f0*j2*k2*p2 - 2*pow(c0,2)*g0*j2*k2*p2 + 
   pow(d0,2)*e0*i2*pow(p2,2) - 2*c0*d0*f0*i2*pow(p2,2) + 
   pow(c0,2)*g0*i2*pow(p2,2) - 2*d0*e0*l2*n0*pow(p2,2) + 
   2*c0*f0*l2*n0*pow(p2,2) + 2*c0*d0*m2*n0*pow(p2,2) - 
   2*c0*c2*pow(n0,2)*pow(p2,2) - 2*pow(c0,2)*n0*n2*pow(p2,2) + 
   2*d0*e0*k2*n0*p2*q2 - 2*c0*f0*k2*n0*p2*q2 - 2*d0*f0*k2*l2*p2*r0 + 
   2*c0*g0*k2*l2*p2*r0 + 2*pow(d0,2)*k2*m2*p2*r0 - 2*c2*d0*k2*n0*p2*r0 - 
   2*c0*d2*k2*n0*p2*r0 - 2*c0*d0*k2*n2*p2*r0 + 2*d0*f0*pow(k2,2)*q2*r0 - 
   2*c0*g0*pow(k2,2)*q2*r0 - 2*d0*d2*pow(k2,2)*pow(r0,2) - 
   2*c0*d0*k2*n0*p2*r2 - 2*pow(d0,2)*pow(k2,2)*r0*r2 + 
   2*d0*e0*k2*l2*p2*s0 - 2*c0*f0*k2*l2*p2*s0 - 2*c0*d0*k2*m2*p2*s0 + 
   4*c0*c2*k2*n0*p2*s0 + 2*pow(c0,2)*k2*n2*p2*s0 - 
   2*d0*e0*pow(k2,2)*q2*s0 + 2*c0*f0*pow(k2,2)*q2*s0 + 
   2*c2*d0*pow(k2,2)*r0*s0 + 2*c0*d2*pow(k2,2)*r0*s0 + 
   2*c0*d0*pow(k2,2)*r2*s0 - 2*c0*c2*pow(k2,2)*pow(s0,2) + 
   2*pow(c0,2)*k2*n0*p2*s2 + 2*c0*d0*pow(k2,2)*r0*s2 - 
   2*pow(c0,2)*pow(k2,2)*s0*s2 + pow(d0,2)*e0*pow(j2,2)*t2 - 
   2*c0*d0*f0*pow(j2,2)*t2 + pow(c0,2)*g0*pow(j2,2)*t2 - 
   pow(d0,2)*e0*i2*o2*t2 + 2*c0*d0*f0*i2*o2*t2 - pow(c0,2)*g0*i2*o2*t2 + 
   2*d0*e0*l2*n0*o2*t2 - 2*c0*f0*l2*n0*o2*t2 - 2*c0*d0*m2*n0*o2*t2 + 
   2*c0*c2*pow(n0,2)*o2*t2 + 2*pow(c0,2)*n0*n2*o2*t2 - 
   2*d0*e0*j2*n0*q2*t2 + 2*c0*f0*j2*n0*q2*t2 + 
   e0*pow(n0,2)*pow(q2,2)*t2 + 2*d0*f0*j2*l2*r0*t2 - 
   2*c0*g0*j2*l2*r0*t2 - 2*pow(d0,2)*j2*m2*r0*t2 + 2*c2*d0*j2*n0*r0*t2 + 
   2*c0*d2*j2*n0*r0*t2 + 2*c0*d0*j2*n2*r0*t2 - 2*d0*f0*i2*q2*r0*t2 + 
   2*c0*g0*i2*q2*r0*t2 + 2*f0*l2*n0*q2*r0*t2 + 2*d0*m2*n0*q2*r0*t2 - 
   2*c2*pow(n0,2)*q2*r0*t2 - 4*c0*n0*n2*q2*r0*t2 + 
   2*d0*d2*i2*pow(r0,2)*t2 + g0*pow(l2,2)*pow(r0,2)*t2 - 
   2*d2*l2*n0*pow(r0,2)*t2 - 2*d0*l2*n2*pow(r0,2)*t2 + 
   2*c0*d0*j2*n0*r2*t2 - 2*c0*pow(n0,2)*q2*r2*t2 + 
   2*pow(d0,2)*i2*r0*r2*t2 - 4*d0*l2*n0*r0*r2*t2 - 2*d0*e0*j2*l2*s0*t2 + 
   2*c0*f0*j2*l2*s0*t2 + 2*c0*d0*j2*m2*s0*t2 - 4*c0*c2*j2*n0*s0*t2 - 
   2*pow(c0,2)*j2*n2*s0*t2 + 2*d0*e0*i2*q2*s0*t2 - 2*c0*f0*i2*q2*s0*t2 - 
   2*e0*l2*n0*q2*s0*t2 + 2*c0*m2*n0*q2*s0*t2 - 2*c2*d0*i2*r0*s0*t2 - 
   2*c0*d2*i2*r0*s0*t2 - 2*f0*pow(l2,2)*r0*s0*t2 + 2*d0*l2*m2*r0*s0*t2 + 
   2*c2*l2*n0*r0*s0*t2 + 2*c0*l2*n2*r0*s0*t2 - 2*c0*d0*i2*r2*s0*t2 + 
   2*c0*l2*n0*r2*s0*t2 + 2*c0*c2*i2*pow(s0,2)*t2 + 
   e0*pow(l2,2)*pow(s0,2)*t2 - 2*c0*l2*m2*pow(s0,2)*t2 - 
   2*pow(c0,2)*j2*n0*s2*t2 - 2*c0*d0*i2*r0*s2*t2 + 2*c0*l2*n0*r0*s2*t2 + 
   2*pow(c0,2)*i2*s0*s2*t2 - 2*pow(f0,2)*k2*l2*o2*u0 + 
   2*e0*g0*k2*l2*o2*u0 + 2*d0*f0*k2*m2*o2*u0 - 2*c0*g0*k2*m2*o2*u0 - 
   2*d2*e0*k2*n0*o2*u0 - 2*d0*e2*k2*n0*o2*u0 + 2*c2*f0*k2*n0*o2*u0 + 
   2*c0*f2*k2*n0*o2*u0 - 2*d0*e0*k2*n2*o2*u0 + 2*c0*f0*k2*n2*o2*u0 + 
   2*pow(f0,2)*j2*l2*p2*u0 - 2*e0*g0*j2*l2*p2*u0 - 2*d0*f0*j2*m2*p2*u0 + 
   2*c0*g0*j2*m2*p2*u0 + 2*d2*e0*j2*n0*p2*u0 + 2*d0*e2*j2*n0*p2*u0 - 
   2*c2*f0*j2*n0*p2*u0 - 2*c0*f2*j2*n0*p2*u0 + 2*d0*e0*j2*n2*p2*u0 - 
   2*c0*f0*j2*n2*p2*u0 + 2*pow(f0,2)*j2*k2*q2*u0 - 2*e0*g0*j2*k2*q2*u0 - 
   2*pow(f0,2)*i2*p2*q2*u0 + 2*e0*g0*i2*p2*q2*u0 + 4*f0*m2*n0*p2*q2*u0 - 
   2*e2*pow(n0,2)*p2*q2*u0 - 4*e0*n0*n2*p2*q2*u0 - 2*d2*f0*j2*k2*r0*u0 - 
   2*d0*f2*j2*k2*r0*u0 + 2*c2*g0*j2*k2*r0*u0 + 2*c0*g2*j2*k2*r0*u0 + 
   2*d2*f0*i2*p2*r0*u0 + 2*d0*f2*i2*p2*r0*u0 - 2*c2*g0*i2*p2*r0*u0 - 
   2*c0*g2*i2*p2*r0*u0 + 2*g0*l2*m2*p2*r0*u0 - 2*f2*l2*n0*p2*r0*u0 - 
   2*d2*m2*n0*p2*r0*u0 - 2*f0*l2*n2*p2*r0*u0 - 2*d0*m2*n2*p2*r0*u0 + 
   4*c2*n0*n2*p2*r0*u0 + 2*c0*pow(n2,2)*p2*r0*u0 + 2*g0*k2*m2*q2*r0*u0 - 
   2*f2*k2*n0*q2*r0*u0 - 2*f0*k2*n2*q2*r0*u0 - 2*g2*k2*l2*pow(r0,2)*u0 + 
   2*d2*k2*n2*pow(r0,2)*u0 - 2*d0*f0*j2*k2*r2*u0 + 2*c0*g0*j2*k2*r2*u0 + 
   2*d0*f0*i2*p2*r2*u0 - 2*c0*g0*i2*p2*r2*u0 - 2*f0*l2*n0*p2*r2*u0 - 
   2*d0*m2*n0*p2*r2*u0 + 2*c2*pow(n0,2)*p2*r2*u0 + 4*c0*n0*n2*p2*r2*u0 - 
   2*f0*k2*n0*q2*r2*u0 - 4*g0*k2*l2*r0*r2*u0 + 4*d2*k2*n0*r0*r2*u0 + 
   4*d0*k2*n2*r0*r2*u0 + 2*d0*k2*n0*pow(r2,2)*u0 + 2*d2*e0*j2*k2*s0*u0 + 
   2*d0*e2*j2*k2*s0*u0 - 2*c2*f0*j2*k2*s0*u0 - 2*c0*f2*j2*k2*s0*u0 - 
   2*d2*e0*i2*p2*s0*u0 - 2*d0*e2*i2*p2*s0*u0 + 2*c2*f0*i2*p2*s0*u0 + 
   2*c0*f2*i2*p2*s0*u0 - 2*f0*l2*m2*p2*s0*u0 + 2*d0*pow(m2,2)*p2*s0*u0 + 
   2*e2*l2*n0*p2*s0*u0 - 2*c2*m2*n0*p2*s0*u0 + 2*e0*l2*n2*p2*s0*u0 - 
   2*c0*m2*n2*p2*s0*u0 - 2*f0*k2*m2*q2*s0*u0 + 2*e2*k2*n0*q2*s0*u0 + 
   2*e0*k2*n2*q2*s0*u0 + 4*f2*k2*l2*r0*s0*u0 - 2*d2*k2*m2*r0*s0*u0 - 
   2*c2*k2*n2*r0*s0*u0 + 4*f0*k2*l2*r2*s0*u0 - 2*d0*k2*m2*r2*s0*u0 - 
   2*c2*k2*n0*r2*s0*u0 - 2*c0*k2*n2*r2*s0*u0 - 2*e2*k2*l2*pow(s0,2)*u0 + 
   2*c2*k2*m2*pow(s0,2)*u0 + 2*d0*e0*j2*k2*s2*u0 - 2*c0*f0*j2*k2*s2*u0 - 
   2*d0*e0*i2*p2*s2*u0 + 2*c0*f0*i2*p2*s2*u0 + 2*e0*l2*n0*p2*s2*u0 - 
   2*c0*m2*n0*p2*s2*u0 + 2*e0*k2*n0*q2*s2*u0 + 4*f0*k2*l2*r0*s2*u0 - 
   2*d0*k2*m2*r0*s2*u0 - 2*c2*k2*n0*r0*s2*u0 - 2*c0*k2*n2*r0*s2*u0 - 
   2*c0*k2*n0*r2*s2*u0 - 4*e0*k2*l2*s0*s2*u0 + 4*c0*k2*m2*s0*s2*u0 - 
   2*f0*f2*pow(j2,2)*pow(u0,2) + e2*g0*pow(j2,2)*pow(u0,2) + 
   e0*g2*pow(j2,2)*pow(u0,2) + 2*f0*f2*i2*o2*pow(u0,2) - 
   e2*g0*i2*o2*pow(u0,2) - e0*g2*i2*o2*pow(u0,2) + 
   g0*pow(m2,2)*o2*pow(u0,2) - 2*f2*m2*n0*o2*pow(u0,2) - 
   2*f0*m2*n2*o2*pow(u0,2) + 2*e2*n0*n2*o2*pow(u0,2) + 
   e0*pow(n2,2)*o2*pow(u0,2) - 2*g2*j2*m2*r0*pow(u0,2) + 
   2*f2*j2*n2*r0*pow(u0,2) - 2*g0*j2*m2*r2*pow(u0,2) + 
   2*f2*j2*n0*r2*pow(u0,2) + 2*f0*j2*n2*r2*pow(u0,2) + 
   2*g2*i2*r0*r2*pow(u0,2) - 2*pow(n2,2)*r0*r2*pow(u0,2) + 
   g0*i2*pow(r2,2)*pow(u0,2) - 2*n0*n2*pow(r2,2)*pow(u0,2) + 
   2*f2*j2*m2*s0*pow(u0,2) - 2*e2*j2*n2*s0*pow(u0,2) - 
   2*f2*i2*r2*s0*pow(u0,2) + 2*m2*n2*r2*s0*pow(u0,2) + 
   2*f0*j2*m2*s2*pow(u0,2) - 2*e2*j2*n0*s2*pow(u0,2) - 
   2*e0*j2*n2*s2*pow(u0,2) - 2*f2*i2*r0*s2*pow(u0,2) + 
   2*m2*n2*r0*s2*pow(u0,2) - 2*f0*i2*r2*s2*pow(u0,2) + 
   2*m2*n0*r2*s2*pow(u0,2) + 2*e2*i2*s0*s2*pow(u0,2) - 
   2*pow(m2,2)*s0*s2*pow(u0,2) + e0*i2*pow(s2,2)*pow(u0,2) - 
   2*d0*e0*k2*n0*o2*u2 + 2*c0*f0*k2*n0*o2*u2 + 2*d0*e0*j2*n0*p2*u2 - 
   2*c0*f0*j2*n0*p2*u2 - 2*e0*pow(n0,2)*p2*q2*u2 - 2*d0*f0*j2*k2*r0*u2 + 
   2*c0*g0*j2*k2*r0*u2 + 2*d0*f0*i2*p2*r0*u2 - 2*c0*g0*i2*p2*r0*u2 - 
   2*f0*l2*n0*p2*r0*u2 - 2*d0*m2*n0*p2*r0*u2 + 2*c2*pow(n0,2)*p2*r0*u2 + 
   4*c0*n0*n2*p2*r0*u2 - 2*f0*k2*n0*q2*r0*u2 - 2*g0*k2*l2*pow(r0,2)*u2 + 
   2*d2*k2*n0*pow(r0,2)*u2 + 2*d0*k2*n2*pow(r0,2)*u2 + 
   2*c0*pow(n0,2)*p2*r2*u2 + 4*d0*k2*n0*r0*r2*u2 + 2*d0*e0*j2*k2*s0*u2 - 
   2*c0*f0*j2*k2*s0*u2 - 2*d0*e0*i2*p2*s0*u2 + 2*c0*f0*i2*p2*s0*u2 + 
   2*e0*l2*n0*p2*s0*u2 - 2*c0*m2*n0*p2*s0*u2 + 2*e0*k2*n0*q2*s0*u2 + 
   4*f0*k2*l2*r0*s0*u2 - 2*d0*k2*m2*r0*s0*u2 - 2*c2*k2*n0*r0*s0*u2 - 
   2*c0*k2*n2*r0*s0*u2 - 2*c0*k2*n0*r2*s0*u2 - 2*e0*k2*l2*pow(s0,2)*u2 + 
   2*c0*k2*m2*pow(s0,2)*u2 - 2*c0*k2*n0*r0*s2*u2 - 
   2*pow(f0,2)*pow(j2,2)*u0*u2 + 2*e0*g0*pow(j2,2)*u0*u2 + 
   2*pow(f0,2)*i2*o2*u0*u2 - 2*e0*g0*i2*o2*u0*u2 - 4*f0*m2*n0*o2*u0*u2 + 
   2*e2*pow(n0,2)*o2*u0*u2 + 4*e0*n0*n2*o2*u0*u2 - 4*g0*j2*m2*r0*u0*u2 + 
   4*f2*j2*n0*r0*u0*u2 + 4*f0*j2*n2*r0*u0*u2 + 2*g2*i2*pow(r0,2)*u0*u2 - 
   2*pow(n2,2)*pow(r0,2)*u0*u2 + 4*f0*j2*n0*r2*u0*u2 + 
   4*g0*i2*r0*r2*u0*u2 - 8*n0*n2*r0*r2*u0*u2 - 
   2*pow(n0,2)*pow(r2,2)*u0*u2 + 4*f0*j2*m2*s0*u0*u2 - 
   4*e2*j2*n0*s0*u0*u2 - 4*e0*j2*n2*s0*u0*u2 - 4*f2*i2*r0*s0*u0*u2 + 
   4*m2*n2*r0*s0*u0*u2 - 4*f0*i2*r2*s0*u0*u2 + 4*m2*n0*r2*s0*u0*u2 + 
   2*e2*i2*pow(s0,2)*u0*u2 - 2*pow(m2,2)*pow(s0,2)*u0*u2 - 
   4*e0*j2*n0*s2*u0*u2 - 4*f0*i2*r0*s2*u0*u2 + 4*m2*n0*r0*s2*u0*u2 + 
   4*e0*i2*s0*s2*u0*u2 + e0*pow(n0,2)*o2*pow(u2,2) + 
   2*f0*j2*n0*r0*pow(u2,2) + g0*i2*pow(r0,2)*pow(u2,2) - 
   2*n0*n2*pow(r0,2)*pow(u2,2) - 2*pow(n0,2)*r0*r2*pow(u2,2) - 
   2*e0*j2*n0*s0*pow(u2,2) - 2*f0*i2*r0*s0*pow(u2,2) + 
   2*m2*n0*r0*s0*pow(u2,2) + e0*i2*pow(s0,2)*pow(u2,2) + 
   2*d0*f0*k2*l2*o2*v0 - 2*c0*g0*k2*l2*o2*v0 - 2*pow(d0,2)*k2*m2*o2*v0 + 
   2*c2*d0*k2*n0*o2*v0 + 2*c0*d2*k2*n0*o2*v0 + 2*c0*d0*k2*n2*o2*v0 - 
   2*d0*f0*j2*l2*p2*v0 + 2*c0*g0*j2*l2*p2*v0 + 2*pow(d0,2)*j2*m2*p2*v0 - 
   2*c2*d0*j2*n0*p2*v0 - 2*c0*d2*j2*n0*p2*v0 - 2*c0*d0*j2*n2*p2*v0 - 
   2*d0*f0*j2*k2*q2*v0 + 2*c0*g0*j2*k2*q2*v0 + 2*d0*f0*i2*p2*q2*v0 - 
   2*c0*g0*i2*p2*q2*v0 - 2*f0*l2*n0*p2*q2*v0 - 2*d0*m2*n0*p2*q2*v0 + 
   2*c2*pow(n0,2)*p2*q2*v0 + 4*c0*n0*n2*p2*q2*v0 + 
   2*f0*k2*n0*pow(q2,2)*v0 + 4*d0*d2*j2*k2*r0*v0 - 4*d0*d2*i2*p2*r0*v0 - 
   2*g0*pow(l2,2)*p2*r0*v0 + 4*d2*l2*n0*p2*r0*v0 + 4*d0*l2*n2*p2*r0*v0 + 
   2*g0*k2*l2*q2*r0*v0 - 2*d2*k2*n0*q2*r0*v0 - 2*d0*k2*n2*q2*r0*v0 + 
   2*pow(d0,2)*j2*k2*r2*v0 - 2*pow(d0,2)*i2*p2*r2*v0 + 
   4*d0*l2*n0*p2*r2*v0 - 2*d0*k2*n0*q2*r2*v0 - 2*c2*d0*j2*k2*s0*v0 - 
   2*c0*d2*j2*k2*s0*v0 + 2*c2*d0*i2*p2*s0*v0 + 2*c0*d2*i2*p2*s0*v0 + 
   2*f0*pow(l2,2)*p2*s0*v0 - 2*d0*l2*m2*p2*s0*v0 - 2*c2*l2*n0*p2*s0*v0 - 
   2*c0*l2*n2*p2*s0*v0 - 2*f0*k2*l2*q2*s0*v0 + 4*d0*k2*m2*q2*s0*v0 - 
   2*c2*k2*n0*q2*s0*v0 - 2*c0*k2*n2*q2*s0*v0 - 2*d2*k2*l2*r0*s0*v0 - 
   2*d0*k2*l2*r2*s0*v0 + 2*c2*k2*l2*pow(s0,2)*v0 - 2*c0*d0*j2*k2*s2*v0 + 
   2*c0*d0*i2*p2*s2*v0 - 2*c0*l2*n0*p2*s2*v0 - 2*c0*k2*n0*q2*s2*v0 - 
   2*d0*k2*l2*r0*s2*v0 + 4*c0*k2*l2*s0*s2*v0 + 2*d2*f0*pow(j2,2)*u0*v0 + 
   2*d0*f2*pow(j2,2)*u0*v0 - 2*c2*g0*pow(j2,2)*u0*v0 - 
   2*c0*g2*pow(j2,2)*u0*v0 - 2*d2*f0*i2*o2*u0*v0 - 2*d0*f2*i2*o2*u0*v0 + 
   2*c2*g0*i2*o2*u0*v0 + 2*c0*g2*i2*o2*u0*v0 - 2*g0*l2*m2*o2*u0*v0 + 
   2*f2*l2*n0*o2*u0*v0 + 2*d2*m2*n0*o2*u0*v0 + 2*f0*l2*n2*o2*u0*v0 + 
   2*d0*m2*n2*o2*u0*v0 - 4*c2*n0*n2*o2*u0*v0 - 2*c0*pow(n2,2)*o2*u0*v0 + 
   2*g0*j2*m2*q2*u0*v0 - 2*f2*j2*n0*q2*u0*v0 - 2*f0*j2*n2*q2*u0*v0 + 
   2*g2*j2*l2*r0*u0*v0 - 2*d2*j2*n2*r0*u0*v0 - 2*g2*i2*q2*r0*u0*v0 + 
   2*pow(n2,2)*q2*r0*u0*v0 + 2*g0*j2*l2*r2*u0*v0 - 2*d2*j2*n0*r2*u0*v0 - 
   2*d0*j2*n2*r2*u0*v0 - 2*g0*i2*q2*r2*u0*v0 + 4*n0*n2*q2*r2*u0*v0 - 
   2*f2*j2*l2*s0*u0*v0 - 2*d2*j2*m2*s0*u0*v0 + 4*c2*j2*n2*s0*u0*v0 + 
   2*f2*i2*q2*s0*u0*v0 - 2*m2*n2*q2*s0*u0*v0 + 2*d2*i2*r2*s0*u0*v0 - 
   2*l2*n2*r2*s0*u0*v0 - 2*f0*j2*l2*s2*u0*v0 - 2*d0*j2*m2*s2*u0*v0 + 
   4*c2*j2*n0*s2*u0*v0 + 4*c0*j2*n2*s2*u0*v0 + 2*f0*i2*q2*s2*u0*v0 - 
   2*m2*n0*q2*s2*u0*v0 + 2*d2*i2*r0*s2*u0*v0 - 2*l2*n2*r0*s2*u0*v0 + 
   2*d0*i2*r2*s2*u0*v0 - 2*l2*n0*r2*s2*u0*v0 - 4*c2*i2*s0*s2*u0*v0 + 
   4*l2*m2*s0*s2*u0*v0 - 2*c0*i2*pow(s2,2)*u0*v0 + 
   2*d0*f0*pow(j2,2)*u2*v0 - 2*c0*g0*pow(j2,2)*u2*v0 - 
   2*d0*f0*i2*o2*u2*v0 + 2*c0*g0*i2*o2*u2*v0 + 2*f0*l2*n0*o2*u2*v0 + 
   2*d0*m2*n0*o2*u2*v0 - 2*c2*pow(n0,2)*o2*u2*v0 - 4*c0*n0*n2*o2*u2*v0 - 
   2*f0*j2*n0*q2*u2*v0 + 2*g0*j2*l2*r0*u2*v0 - 2*d2*j2*n0*r0*u2*v0 - 
   2*d0*j2*n2*r0*u2*v0 - 2*g0*i2*q2*r0*u2*v0 + 4*n0*n2*q2*r0*u2*v0 - 
   2*d0*j2*n0*r2*u2*v0 + 2*pow(n0,2)*q2*r2*u2*v0 - 2*f0*j2*l2*s0*u2*v0 - 
   2*d0*j2*m2*s0*u2*v0 + 4*c2*j2*n0*s0*u2*v0 + 4*c0*j2*n2*s0*u2*v0 + 
   2*f0*i2*q2*s0*u2*v0 - 2*m2*n0*q2*s0*u2*v0 + 2*d2*i2*r0*s0*u2*v0 - 
   2*l2*n2*r0*s0*u2*v0 + 2*d0*i2*r2*s0*u2*v0 - 2*l2*n0*r2*s0*u2*v0 - 
   2*c2*i2*pow(s0,2)*u2*v0 + 2*l2*m2*pow(s0,2)*u2*v0 + 
   4*c0*j2*n0*s2*u2*v0 + 2*d0*i2*r0*s2*u2*v0 - 2*l2*n0*r0*s2*u2*v0 - 
   4*c0*i2*s0*s2*u2*v0 - 2*d0*d2*pow(j2,2)*pow(v0,2) + 
   2*d0*d2*i2*o2*pow(v0,2) + g0*pow(l2,2)*o2*pow(v0,2) - 
   2*d2*l2*n0*o2*pow(v0,2) - 2*d0*l2*n2*o2*pow(v0,2) - 
   2*g0*j2*l2*q2*pow(v0,2) + 2*d2*j2*n0*q2*pow(v0,2) + 
   2*d0*j2*n2*q2*pow(v0,2) + g0*i2*pow(q2,2)*pow(v0,2) - 
   2*n0*n2*pow(q2,2)*pow(v0,2) + 2*d2*j2*l2*s0*pow(v0,2) - 
   2*d2*i2*q2*s0*pow(v0,2) + 2*l2*n2*q2*s0*pow(v0,2) + 
   2*d0*j2*l2*s2*pow(v0,2) - 2*d0*i2*q2*s2*pow(v0,2) + 
   2*l2*n0*q2*s2*pow(v0,2) - 2*pow(l2,2)*s0*s2*pow(v0,2) + 
   2*c0*d0*k2*n0*o2*v2 - 2*c0*d0*j2*n0*p2*v2 + 2*c0*pow(n0,2)*p2*q2*v2 + 
   2*pow(d0,2)*j2*k2*r0*v2 - 2*pow(d0,2)*i2*p2*r0*v2 + 
   4*d0*l2*n0*p2*r0*v2 - 2*d0*k2*n0*q2*r0*v2 - 2*c0*d0*j2*k2*s0*v2 + 
   2*c0*d0*i2*p2*s0*v2 - 2*c0*l2*n0*p2*s0*v2 - 2*c0*k2*n0*q2*s0*v2 - 
   2*d0*k2*l2*r0*s0*v2 + 2*c0*k2*l2*pow(s0,2)*v2 + 
   2*d0*f0*pow(j2,2)*u0*v2 - 2*c0*g0*pow(j2,2)*u0*v2 - 
   2*d0*f0*i2*o2*u0*v2 + 2*c0*g0*i2*o2*u0*v2 + 2*f0*l2*n0*o2*u0*v2 + 
   2*d0*m2*n0*o2*u0*v2 - 2*c2*pow(n0,2)*o2*u0*v2 - 4*c0*n0*n2*o2*u0*v2 - 
   2*f0*j2*n0*q2*u0*v2 + 2*g0*j2*l2*r0*u0*v2 - 2*d2*j2*n0*r0*u0*v2 - 
   2*d0*j2*n2*r0*u0*v2 - 2*g0*i2*q2*r0*u0*v2 + 4*n0*n2*q2*r0*u0*v2 - 
   2*d0*j2*n0*r2*u0*v2 + 2*pow(n0,2)*q2*r2*u0*v2 - 2*f0*j2*l2*s0*u0*v2 - 
   2*d0*j2*m2*s0*u0*v2 + 4*c2*j2*n0*s0*u0*v2 + 4*c0*j2*n2*s0*u0*v2 + 
   2*f0*i2*q2*s0*u0*v2 - 2*m2*n0*q2*s0*u0*v2 + 2*d2*i2*r0*s0*u0*v2 - 
   2*l2*n2*r0*s0*u0*v2 + 2*d0*i2*r2*s0*u0*v2 - 2*l2*n0*r2*s0*u0*v2 - 
   2*c2*i2*pow(s0,2)*u0*v2 + 2*l2*m2*pow(s0,2)*u0*v2 + 
   4*c0*j2*n0*s2*u0*v2 + 2*d0*i2*r0*s2*u0*v2 - 2*l2*n0*r0*s2*u0*v2 - 
   4*c0*i2*s0*s2*u0*v2 - 2*c0*pow(n0,2)*o2*u2*v2 - 2*d0*j2*n0*r0*u2*v2 + 
   2*pow(n0,2)*q2*r0*u2*v2 + 4*c0*j2*n0*s0*u2*v2 + 2*d0*i2*r0*s0*u2*v2 - 
   2*l2*n0*r0*s0*u2*v2 - 2*c0*i2*pow(s0,2)*u2*v2 - 
   2*pow(d0,2)*pow(j2,2)*v0*v2 + 2*pow(d0,2)*i2*o2*v0*v2 - 
   4*d0*l2*n0*o2*v0*v2 + 4*d0*j2*n0*q2*v0*v2 - 
   2*pow(n0,2)*pow(q2,2)*v0*v2 + 4*d0*j2*l2*s0*v0*v2 - 
   4*d0*i2*q2*s0*v0*v2 + 4*l2*n0*q2*s0*v0*v2 - 
   2*pow(l2,2)*pow(s0,2)*v0*v2 - 2*d0*e0*k2*l2*o2*w0 + 
   2*c0*f0*k2*l2*o2*w0 + 2*c0*d0*k2*m2*o2*w0 - 4*c0*c2*k2*n0*o2*w0 - 
   2*pow(c0,2)*k2*n2*o2*w0 + 2*d0*e0*j2*l2*p2*w0 - 2*c0*f0*j2*l2*p2*w0 - 
   2*c0*d0*j2*m2*p2*w0 + 4*c0*c2*j2*n0*p2*w0 + 2*pow(c0,2)*j2*n2*p2*w0 + 
   2*d0*e0*j2*k2*q2*w0 - 2*c0*f0*j2*k2*q2*w0 - 2*d0*e0*i2*p2*q2*w0 + 
   2*c0*f0*i2*p2*q2*w0 + 2*e0*l2*n0*p2*q2*w0 - 2*c0*m2*n0*p2*q2*w0 - 
   2*e0*k2*n0*pow(q2,2)*w0 - 2*c2*d0*j2*k2*r0*w0 - 2*c0*d2*j2*k2*r0*w0 + 
   2*c2*d0*i2*p2*r0*w0 + 2*c0*d2*i2*p2*r0*w0 + 2*f0*pow(l2,2)*p2*r0*w0 - 
   2*d0*l2*m2*p2*r0*w0 - 2*c2*l2*n0*p2*r0*w0 - 2*c0*l2*n2*p2*r0*w0 - 
   2*f0*k2*l2*q2*r0*w0 - 2*d0*k2*m2*q2*r0*w0 + 4*c2*k2*n0*q2*r0*w0 + 
   4*c0*k2*n2*q2*r0*w0 + 2*d2*k2*l2*pow(r0,2)*w0 - 2*c0*d0*j2*k2*r2*w0 + 
   2*c0*d0*i2*p2*r2*w0 - 2*c0*l2*n0*p2*r2*w0 + 4*c0*k2*n0*q2*r2*w0 + 
   4*d0*k2*l2*r0*r2*w0 + 4*c0*c2*j2*k2*s0*w0 - 4*c0*c2*i2*p2*s0*w0 - 
   2*e0*pow(l2,2)*p2*s0*w0 + 4*c0*l2*m2*p2*s0*w0 + 2*e0*k2*l2*q2*s0*w0 - 
   2*c0*k2*m2*q2*s0*w0 - 2*c2*k2*l2*r0*s0*w0 - 2*c0*k2*l2*r2*s0*w0 + 
   2*pow(c0,2)*j2*k2*s2*w0 - 2*pow(c0,2)*i2*p2*s2*w0 - 
   2*c0*k2*l2*r0*s2*w0 - 2*d2*e0*pow(j2,2)*u0*w0 - 
   2*d0*e2*pow(j2,2)*u0*w0 + 2*c2*f0*pow(j2,2)*u0*w0 + 
   2*c0*f2*pow(j2,2)*u0*w0 + 2*d2*e0*i2*o2*u0*w0 + 2*d0*e2*i2*o2*u0*w0 - 
   2*c2*f0*i2*o2*u0*w0 - 2*c0*f2*i2*o2*u0*w0 + 2*f0*l2*m2*o2*u0*w0 - 
   2*d0*pow(m2,2)*o2*u0*w0 - 2*e2*l2*n0*o2*u0*w0 + 2*c2*m2*n0*o2*u0*w0 - 
   2*e0*l2*n2*o2*u0*w0 + 2*c0*m2*n2*o2*u0*w0 - 2*f0*j2*m2*q2*u0*w0 + 
   2*e2*j2*n0*q2*u0*w0 + 2*e0*j2*n2*q2*u0*w0 - 2*f2*j2*l2*r0*u0*w0 + 
   4*d2*j2*m2*r0*u0*w0 - 2*c2*j2*n2*r0*u0*w0 + 2*f2*i2*q2*r0*u0*w0 - 
   2*m2*n2*q2*r0*u0*w0 - 2*f0*j2*l2*r2*u0*w0 + 4*d0*j2*m2*r2*u0*w0 - 
   2*c2*j2*n0*r2*u0*w0 - 2*c0*j2*n2*r2*u0*w0 + 2*f0*i2*q2*r2*u0*w0 - 
   2*m2*n0*q2*r2*u0*w0 - 4*d2*i2*r0*r2*u0*w0 + 4*l2*n2*r0*r2*u0*w0 - 
   2*d0*i2*pow(r2,2)*u0*w0 + 2*l2*n0*pow(r2,2)*u0*w0 + 
   2*e2*j2*l2*s0*u0*w0 - 2*c2*j2*m2*s0*u0*w0 - 2*e2*i2*q2*s0*u0*w0 + 
   2*pow(m2,2)*q2*s0*u0*w0 + 2*c2*i2*r2*s0*u0*w0 - 2*l2*m2*r2*s0*u0*w0 + 
   2*e0*j2*l2*s2*u0*w0 - 2*c0*j2*m2*s2*u0*w0 - 2*e0*i2*q2*s2*u0*w0 + 
   2*c2*i2*r0*s2*u0*w0 - 2*l2*m2*r0*s2*u0*w0 + 2*c0*i2*r2*s2*u0*w0 - 
   2*d0*e0*pow(j2,2)*u2*w0 + 2*c0*f0*pow(j2,2)*u2*w0 + 
   2*d0*e0*i2*o2*u2*w0 - 2*c0*f0*i2*o2*u2*w0 - 2*e0*l2*n0*o2*u2*w0 + 
   2*c0*m2*n0*o2*u2*w0 + 2*e0*j2*n0*q2*u2*w0 - 2*f0*j2*l2*r0*u2*w0 + 
   4*d0*j2*m2*r0*u2*w0 - 2*c2*j2*n0*r0*u2*w0 - 2*c0*j2*n2*r0*u2*w0 + 
   2*f0*i2*q2*r0*u2*w0 - 2*m2*n0*q2*r0*u2*w0 - 2*d2*i2*pow(r0,2)*u2*w0 + 
   2*l2*n2*pow(r0,2)*u2*w0 - 2*c0*j2*n0*r2*u2*w0 - 4*d0*i2*r0*r2*u2*w0 + 
   4*l2*n0*r0*r2*u2*w0 + 2*e0*j2*l2*s0*u2*w0 - 2*c0*j2*m2*s0*u2*w0 - 
   2*e0*i2*q2*s0*u2*w0 + 2*c2*i2*r0*s0*u2*w0 - 2*l2*m2*r0*s0*u2*w0 + 
   2*c0*i2*r2*s0*u2*w0 + 2*c0*i2*r0*s2*u2*w0 + 2*c2*d0*pow(j2,2)*v0*w0 + 
   2*c0*d2*pow(j2,2)*v0*w0 - 2*c2*d0*i2*o2*v0*w0 - 2*c0*d2*i2*o2*v0*w0 - 
   2*f0*pow(l2,2)*o2*v0*w0 + 2*d0*l2*m2*o2*v0*w0 + 2*c2*l2*n0*o2*v0*w0 + 
   2*c0*l2*n2*o2*v0*w0 + 4*f0*j2*l2*q2*v0*w0 - 2*d0*j2*m2*q2*v0*w0 - 
   2*c2*j2*n0*q2*v0*w0 - 2*c0*j2*n2*q2*v0*w0 - 2*f0*i2*pow(q2,2)*v0*w0 + 
   2*m2*n0*pow(q2,2)*v0*w0 - 2*d2*j2*l2*r0*v0*w0 + 2*d2*i2*q2*r0*v0*w0 - 
   2*l2*n2*q2*r0*v0*w0 - 2*d0*j2*l2*r2*v0*w0 + 2*d0*i2*q2*r2*v0*w0 - 
   2*l2*n0*q2*r2*v0*w0 - 2*c2*j2*l2*s0*v0*w0 + 2*c2*i2*q2*s0*v0*w0 - 
   2*l2*m2*q2*s0*v0*w0 + 2*pow(l2,2)*r2*s0*v0*w0 - 2*c0*j2*l2*s2*v0*w0 + 
   2*c0*i2*q2*s2*v0*w0 + 2*pow(l2,2)*r0*s2*v0*w0 + 
   2*c0*d0*pow(j2,2)*v2*w0 - 2*c0*d0*i2*o2*v2*w0 + 2*c0*l2*n0*o2*v2*w0 - 
   2*c0*j2*n0*q2*v2*w0 - 2*d0*j2*l2*r0*v2*w0 + 2*d0*i2*q2*r0*v2*w0 - 
   2*l2*n0*q2*r0*v2*w0 - 2*c0*j2*l2*s0*v2*w0 + 2*c0*i2*q2*s0*v2*w0 + 
   2*pow(l2,2)*r0*s0*v2*w0 - 2*c0*c2*pow(j2,2)*pow(w0,2) + 
   2*c0*c2*i2*o2*pow(w0,2) + e0*pow(l2,2)*o2*pow(w0,2) - 
   2*c0*l2*m2*o2*pow(w0,2) - 2*e0*j2*l2*q2*pow(w0,2) + 
   2*c0*j2*m2*q2*pow(w0,2) + e0*i2*pow(q2,2)*pow(w0,2) + 
   2*c2*j2*l2*r0*pow(w0,2) - 2*c2*i2*q2*r0*pow(w0,2) + 
   2*l2*m2*q2*r0*pow(w0,2) + 2*c0*j2*l2*r2*pow(w0,2) - 
   2*c0*i2*q2*r2*pow(w0,2) - 2*pow(l2,2)*r0*r2*pow(w0,2) - 
   2*pow(c0,2)*k2*n0*o2*w2 + 2*pow(c0,2)*j2*n0*p2*w2 - 
   2*c0*d0*j2*k2*r0*w2 + 2*c0*d0*i2*p2*r0*w2 - 2*c0*l2*n0*p2*r0*w2 + 
   4*c0*k2*n0*q2*r0*w2 + 2*d0*k2*l2*pow(r0,2)*w2 + 
   2*pow(c0,2)*j2*k2*s0*w2 - 2*pow(c0,2)*i2*p2*s0*w2 - 
   2*c0*k2*l2*r0*s0*w2 - 2*d0*e0*pow(j2,2)*u0*w2 + 
   2*c0*f0*pow(j2,2)*u0*w2 + 2*d0*e0*i2*o2*u0*w2 - 2*c0*f0*i2*o2*u0*w2 - 
   2*e0*l2*n0*o2*u0*w2 + 2*c0*m2*n0*o2*u0*w2 + 2*e0*j2*n0*q2*u0*w2 - 
   2*f0*j2*l2*r0*u0*w2 + 4*d0*j2*m2*r0*u0*w2 - 2*c2*j2*n0*r0*u0*w2 - 
   2*c0*j2*n2*r0*u0*w2 + 2*f0*i2*q2*r0*u0*w2 - 2*m2*n0*q2*r0*u0*w2 - 
   2*d2*i2*pow(r0,2)*u0*w2 + 2*l2*n2*pow(r0,2)*u0*w2 - 
   2*c0*j2*n0*r2*u0*w2 - 4*d0*i2*r0*r2*u0*w2 + 4*l2*n0*r0*r2*u0*w2 + 
   2*e0*j2*l2*s0*u0*w2 - 2*c0*j2*m2*s0*u0*w2 - 2*e0*i2*q2*s0*u0*w2 + 
   2*c2*i2*r0*s0*u0*w2 - 2*l2*m2*r0*s0*u0*w2 + 2*c0*i2*r2*s0*u0*w2 + 
   2*c0*i2*r0*s2*u0*w2 - 2*c0*j2*n0*r0*u2*w2 - 2*d0*i2*pow(r0,2)*u2*w2 + 
   2*l2*n0*pow(r0,2)*u2*w2 + 2*c0*i2*r0*s0*u2*w2 + 
   2*c0*d0*pow(j2,2)*v0*w2 - 2*c0*d0*i2*o2*v0*w2 + 2*c0*l2*n0*o2*v0*w2 - 
   2*c0*j2*n0*q2*v0*w2 - 2*d0*j2*l2*r0*v0*w2 + 2*d0*i2*q2*r0*v0*w2 - 
   2*l2*n0*q2*r0*v0*w2 - 2*c0*j2*l2*s0*v0*w2 + 2*c0*i2*q2*s0*v0*w2 + 
   2*pow(l2,2)*r0*s0*v0*w2 - 2*pow(c0,2)*pow(j2,2)*w0*w2 + 
   2*pow(c0,2)*i2*o2*w0*w2 + 4*c0*j2*l2*r0*w0*w2 - 4*c0*i2*q2*r0*w0*w2 - 
   2*pow(l2,2)*pow(r0,2)*w0*w2 + pow(f0,2)*pow(k2,2)*o2*z0 - 
   e0*g0*pow(k2,2)*o2*z0 - 2*pow(f0,2)*j2*k2*p2*z0 + 
   2*e0*g0*j2*k2*p2*z0 + pow(f0,2)*i2*pow(p2,2)*z0 - 
   e0*g0*i2*pow(p2,2)*z0 - 2*f0*m2*n0*pow(p2,2)*z0 + 
   e2*pow(n0,2)*pow(p2,2)*z0 + 2*e0*n0*n2*pow(p2,2)*z0 - 
   2*g0*k2*m2*p2*r0*z0 + 2*f2*k2*n0*p2*r0*z0 + 2*f0*k2*n2*p2*r0*z0 + 
   g2*pow(k2,2)*pow(r0,2)*z0 + 2*f0*k2*n0*p2*r2*z0 + 
   2*g0*pow(k2,2)*r0*r2*z0 + 2*f0*k2*m2*p2*s0*z0 - 2*e2*k2*n0*p2*s0*z0 - 
   2*e0*k2*n2*p2*s0*z0 - 2*f2*pow(k2,2)*r0*s0*z0 - 
   2*f0*pow(k2,2)*r2*s0*z0 + e2*pow(k2,2)*pow(s0,2)*z0 - 
   2*e0*k2*n0*p2*s2*z0 - 2*f0*pow(k2,2)*r0*s2*z0 + 
   2*e0*pow(k2,2)*s0*s2*z0 + pow(f0,2)*pow(j2,2)*t2*z0 - 
   e0*g0*pow(j2,2)*t2*z0 - pow(f0,2)*i2*o2*t2*z0 + e0*g0*i2*o2*t2*z0 + 
   2*f0*m2*n0*o2*t2*z0 - e2*pow(n0,2)*o2*t2*z0 - 2*e0*n0*n2*o2*t2*z0 + 
   2*g0*j2*m2*r0*t2*z0 - 2*f2*j2*n0*r0*t2*z0 - 2*f0*j2*n2*r0*t2*z0 - 
   g2*i2*pow(r0,2)*t2*z0 + pow(n2,2)*pow(r0,2)*t2*z0 - 
   2*f0*j2*n0*r2*t2*z0 - 2*g0*i2*r0*r2*t2*z0 + 4*n0*n2*r0*r2*t2*z0 + 
   pow(n0,2)*pow(r2,2)*t2*z0 - 2*f0*j2*m2*s0*t2*z0 + 
   2*e2*j2*n0*s0*t2*z0 + 2*e0*j2*n2*s0*t2*z0 + 2*f2*i2*r0*s0*t2*z0 - 
   2*m2*n2*r0*s0*t2*z0 + 2*f0*i2*r2*s0*t2*z0 - 2*m2*n0*r2*s0*t2*z0 - 
   e2*i2*pow(s0,2)*t2*z0 + pow(m2,2)*pow(s0,2)*t2*z0 + 
   2*e0*j2*n0*s2*t2*z0 + 2*f0*i2*r0*s2*t2*z0 - 2*m2*n0*r0*s2*t2*z0 - 
   2*e0*i2*s0*s2*t2*z0 + 2*g0*k2*m2*o2*v0*z0 - 2*f2*k2*n0*o2*v0*z0 - 
   2*f0*k2*n2*o2*v0*z0 - 2*g0*j2*m2*p2*v0*z0 + 2*f2*j2*n0*p2*v0*z0 + 
   2*f0*j2*n2*p2*v0*z0 - 2*g2*j2*k2*r0*v0*z0 + 2*g2*i2*p2*r0*v0*z0 - 
   2*pow(n2,2)*p2*r0*v0*z0 - 2*g0*j2*k2*r2*v0*z0 + 2*g0*i2*p2*r2*v0*z0 - 
   4*n0*n2*p2*r2*v0*z0 + 2*f2*j2*k2*s0*v0*z0 - 2*f2*i2*p2*s0*v0*z0 + 
   2*m2*n2*p2*s0*v0*z0 + 2*k2*n2*r2*s0*v0*z0 + 2*f0*j2*k2*s2*v0*z0 - 
   2*f0*i2*p2*s2*v0*z0 + 2*m2*n0*p2*s2*v0*z0 + 2*k2*n2*r0*s2*v0*z0 + 
   2*k2*n0*r2*s2*v0*z0 - 4*k2*m2*s0*s2*v0*z0 + 
   g2*pow(j2,2)*pow(v0,2)*z0 - g2*i2*o2*pow(v0,2)*z0 + 
   pow(n2,2)*o2*pow(v0,2)*z0 - 2*j2*n2*s2*pow(v0,2)*z0 + 
   i2*pow(s2,2)*pow(v0,2)*z0 - 2*f0*k2*n0*o2*v2*z0 + 
   2*f0*j2*n0*p2*v2*z0 - 2*g0*j2*k2*r0*v2*z0 + 2*g0*i2*p2*r0*v2*z0 - 
   4*n0*n2*p2*r0*v2*z0 - 2*pow(n0,2)*p2*r2*v2*z0 + 2*f0*j2*k2*s0*v2*z0 - 
   2*f0*i2*p2*s0*v2*z0 + 2*m2*n0*p2*s0*v2*z0 + 2*k2*n2*r0*s0*v2*z0 + 
   2*k2*n0*r2*s0*v2*z0 - 2*k2*m2*pow(s0,2)*v2*z0 + 2*k2*n0*r0*s2*v2*z0 + 
   2*g0*pow(j2,2)*v0*v2*z0 - 2*g0*i2*o2*v0*v2*z0 + 4*n0*n2*o2*v0*v2*z0 - 
   4*j2*n2*s0*v0*v2*z0 - 4*j2*n0*s2*v0*v2*z0 + 4*i2*s0*s2*v0*v2*z0 + 
   pow(n0,2)*o2*pow(v2,2)*z0 - 2*j2*n0*s0*pow(v2,2)*z0 + 
   i2*pow(s0,2)*pow(v2,2)*z0 - 2*f0*k2*m2*o2*w0*z0 + 
   2*e2*k2*n0*o2*w0*z0 + 2*e0*k2*n2*o2*w0*z0 + 2*f0*j2*m2*p2*w0*z0 - 
   2*e2*j2*n0*p2*w0*z0 - 2*e0*j2*n2*p2*w0*z0 + 2*f2*j2*k2*r0*w0*z0 - 
   2*f2*i2*p2*r0*w0*z0 + 2*m2*n2*p2*r0*w0*z0 + 2*f0*j2*k2*r2*w0*z0 - 
   2*f0*i2*p2*r2*w0*z0 + 2*m2*n0*p2*r2*w0*z0 - 4*k2*n2*r0*r2*w0*z0 - 
   2*k2*n0*pow(r2,2)*w0*z0 - 2*e2*j2*k2*s0*w0*z0 + 2*e2*i2*p2*s0*w0*z0 - 
   2*pow(m2,2)*p2*s0*w0*z0 + 2*k2*m2*r2*s0*w0*z0 - 2*e0*j2*k2*s2*w0*z0 + 
   2*e0*i2*p2*s2*w0*z0 + 2*k2*m2*r0*s2*w0*z0 - 2*f2*pow(j2,2)*v0*w0*z0 + 
   2*f2*i2*o2*v0*w0*z0 - 2*m2*n2*o2*v0*w0*z0 + 2*j2*n2*r2*v0*w0*z0 + 
   2*j2*m2*s2*v0*w0*z0 - 2*i2*r2*s2*v0*w0*z0 - 2*f0*pow(j2,2)*v2*w0*z0 + 
   2*f0*i2*o2*v2*w0*z0 - 2*m2*n0*o2*v2*w0*z0 + 2*j2*n2*r0*v2*w0*z0 + 
   2*j2*n0*r2*v2*w0*z0 + 2*j2*m2*s0*v2*w0*z0 - 2*i2*r2*s0*v2*w0*z0 - 
   2*i2*r0*s2*v2*w0*z0 + e2*pow(j2,2)*pow(w0,2)*z0 - 
   e2*i2*o2*pow(w0,2)*z0 + pow(m2,2)*o2*pow(w0,2)*z0 - 
   2*j2*m2*r2*pow(w0,2)*z0 + i2*pow(r2,2)*pow(w0,2)*z0 + 
   2*e0*k2*n0*o2*w2*z0 - 2*e0*j2*n0*p2*w2*z0 + 2*f0*j2*k2*r0*w2*z0 - 
   2*f0*i2*p2*r0*w2*z0 + 2*m2*n0*p2*r0*w2*z0 - 2*k2*n2*pow(r0,2)*w2*z0 - 
   4*k2*n0*r0*r2*w2*z0 - 2*e0*j2*k2*s0*w2*z0 + 2*e0*i2*p2*s0*w2*z0 + 
   2*k2*m2*r0*s0*w2*z0 - 2*f0*pow(j2,2)*v0*w2*z0 + 2*f0*i2*o2*v0*w2*z0 - 
   2*m2*n0*o2*v0*w2*z0 + 2*j2*n2*r0*v0*w2*z0 + 2*j2*n0*r2*v0*w2*z0 + 
   2*j2*m2*s0*v0*w2*z0 - 2*i2*r2*s0*v0*w2*z0 - 2*i2*r0*s2*v0*w2*z0 + 
   2*j2*n0*r0*v2*w2*z0 - 2*i2*r0*s0*v2*w2*z0 + 2*e0*pow(j2,2)*w0*w2*z0 - 
   2*e0*i2*o2*w0*w2*z0 - 4*j2*m2*r0*w0*w2*z0 + 4*i2*r0*r2*w0*w2*z0 + 
   i2*pow(r0,2)*pow(w2,2)*z0 + e0*pow(n0,2)*pow(p2,2)*z2 + 
   2*f0*k2*n0*p2*r0*z2 + g0*pow(k2,2)*pow(r0,2)*z2 - 
   2*e0*k2*n0*p2*s0*z2 - 2*f0*pow(k2,2)*r0*s0*z2 + 
   e0*pow(k2,2)*pow(s0,2)*z2 - e0*pow(n0,2)*o2*t2*z2 - 
   2*f0*j2*n0*r0*t2*z2 - g0*i2*pow(r0,2)*t2*z2 + 
   2*n0*n2*pow(r0,2)*t2*z2 + 2*pow(n0,2)*r0*r2*t2*z2 + 
   2*e0*j2*n0*s0*t2*z2 + 2*f0*i2*r0*s0*t2*z2 - 2*m2*n0*r0*s0*t2*z2 - 
   e0*i2*pow(s0,2)*t2*z2 - 2*f0*k2*n0*o2*v0*z2 + 2*f0*j2*n0*p2*v0*z2 - 
   2*g0*j2*k2*r0*v0*z2 + 2*g0*i2*p2*r0*v0*z2 - 4*n0*n2*p2*r0*v0*z2 - 
   2*pow(n0,2)*p2*r2*v0*z2 + 2*f0*j2*k2*s0*v0*z2 - 2*f0*i2*p2*s0*v0*z2 + 
   2*m2*n0*p2*s0*v0*z2 + 2*k2*n2*r0*s0*v0*z2 + 2*k2*n0*r2*s0*v0*z2 - 
   2*k2*m2*pow(s0,2)*v0*z2 + 2*k2*n0*r0*s2*v0*z2 + 
   g0*pow(j2,2)*pow(v0,2)*z2 - g0*i2*o2*pow(v0,2)*z2 + 
   2*n0*n2*o2*pow(v0,2)*z2 - 2*j2*n2*s0*pow(v0,2)*z2 - 
   2*j2*n0*s2*pow(v0,2)*z2 + 2*i2*s0*s2*pow(v0,2)*z2 - 
   2*pow(n0,2)*p2*r0*v2*z2 + 2*k2*n0*r0*s0*v2*z2 + 
   2*pow(n0,2)*o2*v0*v2*z2 - 4*j2*n0*s0*v0*v2*z2 + 
   2*i2*pow(s0,2)*v0*v2*z2 + 2*e0*k2*n0*o2*w0*z2 - 2*e0*j2*n0*p2*w0*z2 + 
   2*f0*j2*k2*r0*w0*z2 - 2*f0*i2*p2*r0*w0*z2 + 2*m2*n0*p2*r0*w0*z2 - 
   2*k2*n2*pow(r0,2)*w0*z2 - 4*k2*n0*r0*r2*w0*z2 - 2*e0*j2*k2*s0*w0*z2 + 
   2*e0*i2*p2*s0*w0*z2 + 2*k2*m2*r0*s0*w0*z2 - 2*f0*pow(j2,2)*v0*w0*z2 + 
   2*f0*i2*o2*v0*w0*z2 - 2*m2*n0*o2*v0*w0*z2 + 2*j2*n2*r0*v0*w0*z2 + 
   2*j2*n0*r2*v0*w0*z2 + 2*j2*m2*s0*v0*w0*z2 - 2*i2*r2*s0*v0*w0*z2 - 
   2*i2*r0*s2*v0*w0*z2 + 2*j2*n0*r0*v2*w0*z2 - 2*i2*r0*s0*v2*w0*z2 + 
   e0*pow(j2,2)*pow(w0,2)*z2 - e0*i2*o2*pow(w0,2)*z2 - 
   2*j2*m2*r0*pow(w0,2)*z2 + 2*i2*r0*r2*pow(w0,2)*z2 - 
   2*k2*n0*pow(r0,2)*w2*z2 + 2*j2*n0*r0*v0*w2*z2 - 2*i2*r0*s0*v0*w2*z2 + 
   2*i2*pow(r0,2)*w0*w2*z2
;

	k[7]  = 2*d0*d1*e0*pow(k2,2)*o2 + pow(d0,2)*e1*pow(k2,2)*o2 - 
   2*c1*d0*f0*pow(k2,2)*o2 - 2*c0*d1*f0*pow(k2,2)*o2 - 
   2*c0*d0*f1*pow(k2,2)*o2 + 2*c0*c1*g0*pow(k2,2)*o2 + 
   pow(c0,2)*g1*pow(k2,2)*o2 - 4*d0*d1*e0*j2*k2*p2 - 
   2*pow(d0,2)*e1*j2*k2*p2 + 4*c1*d0*f0*j2*k2*p2 + 4*c0*d1*f0*j2*k2*p2 + 
   4*c0*d0*f1*j2*k2*p2 - 4*c0*c1*g0*j2*k2*p2 - 2*pow(c0,2)*g1*j2*k2*p2 + 
   2*d0*d1*e0*i2*pow(p2,2) + pow(d0,2)*e1*i2*pow(p2,2) - 
   2*c1*d0*f0*i2*pow(p2,2) - 2*c0*d1*f0*i2*pow(p2,2) - 
   2*c0*d0*f1*i2*pow(p2,2) + 2*c0*c1*g0*i2*pow(p2,2) + 
   pow(c0,2)*g1*i2*pow(p2,2) - 2*d1*e0*l2*n0*pow(p2,2) - 
   2*d0*e1*l2*n0*pow(p2,2) + 2*c1*f0*l2*n0*pow(p2,2) + 
   2*c0*f1*l2*n0*pow(p2,2) + 2*c1*d0*m2*n0*pow(p2,2) + 
   2*c0*d1*m2*n0*pow(p2,2) - 2*c1*c2*pow(n0,2)*pow(p2,2) - 
   2*d0*e0*l2*n1*pow(p2,2) + 2*c0*f0*l2*n1*pow(p2,2) + 
   2*c0*d0*m2*n1*pow(p2,2) - 4*c0*c2*n0*n1*pow(p2,2) - 
   4*c0*c1*n0*n2*pow(p2,2) - 2*pow(c0,2)*n1*n2*pow(p2,2) + 
   2*d1*e0*k2*n0*p2*q2 + 2*d0*e1*k2*n0*p2*q2 - 2*c1*f0*k2*n0*p2*q2 - 
   2*c0*f1*k2*n0*p2*q2 + 2*d0*e0*k2*n1*p2*q2 - 2*c0*f0*k2*n1*p2*q2 - 
   2*d1*f0*k2*l2*p2*r0 - 2*d0*f1*k2*l2*p2*r0 + 2*c1*g0*k2*l2*p2*r0 + 
   2*c0*g1*k2*l2*p2*r0 + 4*d0*d1*k2*m2*p2*r0 - 2*c2*d1*k2*n0*p2*r0 - 
   2*c1*d2*k2*n0*p2*r0 - 2*c2*d0*k2*n1*p2*r0 - 2*c0*d2*k2*n1*p2*r0 - 
   2*c1*d0*k2*n2*p2*r0 - 2*c0*d1*k2*n2*p2*r0 + 2*d1*f0*pow(k2,2)*q2*r0 + 
   2*d0*f1*pow(k2,2)*q2*r0 - 2*c1*g0*pow(k2,2)*q2*r0 - 
   2*c0*g1*pow(k2,2)*q2*r0 - 2*d1*d2*pow(k2,2)*pow(r0,2) - 
   2*d0*f0*k2*l2*p2*r1 + 2*c0*g0*k2*l2*p2*r1 + 2*pow(d0,2)*k2*m2*p2*r1 - 
   2*c2*d0*k2*n0*p2*r1 - 2*c0*d2*k2*n0*p2*r1 - 2*c0*d0*k2*n2*p2*r1 + 
   2*d0*f0*pow(k2,2)*q2*r1 - 2*c0*g0*pow(k2,2)*q2*r1 - 
   4*d0*d2*pow(k2,2)*r0*r1 - 2*c1*d0*k2*n0*p2*r2 - 2*c0*d1*k2*n0*p2*r2 - 
   2*c0*d0*k2*n1*p2*r2 - 4*d0*d1*pow(k2,2)*r0*r2 - 
   2*pow(d0,2)*pow(k2,2)*r1*r2 + 2*d1*e0*k2*l2*p2*s0 + 
   2*d0*e1*k2*l2*p2*s0 - 2*c1*f0*k2*l2*p2*s0 - 2*c0*f1*k2*l2*p2*s0 - 
   2*c1*d0*k2*m2*p2*s0 - 2*c0*d1*k2*m2*p2*s0 + 4*c1*c2*k2*n0*p2*s0 + 
   4*c0*c2*k2*n1*p2*s0 + 4*c0*c1*k2*n2*p2*s0 - 2*d1*e0*pow(k2,2)*q2*s0 - 
   2*d0*e1*pow(k2,2)*q2*s0 + 2*c1*f0*pow(k2,2)*q2*s0 + 
   2*c0*f1*pow(k2,2)*q2*s0 + 2*c2*d1*pow(k2,2)*r0*s0 + 
   2*c1*d2*pow(k2,2)*r0*s0 + 2*c2*d0*pow(k2,2)*r1*s0 + 
   2*c0*d2*pow(k2,2)*r1*s0 + 2*c1*d0*pow(k2,2)*r2*s0 + 
   2*c0*d1*pow(k2,2)*r2*s0 - 2*c1*c2*pow(k2,2)*pow(s0,2) + 
   2*d0*e0*k2*l2*p2*s1 - 2*c0*f0*k2*l2*p2*s1 - 2*c0*d0*k2*m2*p2*s1 + 
   4*c0*c2*k2*n0*p2*s1 + 2*pow(c0,2)*k2*n2*p2*s1 - 
   2*d0*e0*pow(k2,2)*q2*s1 + 2*c0*f0*pow(k2,2)*q2*s1 + 
   2*c2*d0*pow(k2,2)*r0*s1 + 2*c0*d2*pow(k2,2)*r0*s1 + 
   2*c0*d0*pow(k2,2)*r2*s1 - 4*c0*c2*pow(k2,2)*s0*s1 + 
   4*c0*c1*k2*n0*p2*s2 + 2*pow(c0,2)*k2*n1*p2*s2 + 
   2*c1*d0*pow(k2,2)*r0*s2 + 2*c0*d1*pow(k2,2)*r0*s2 + 
   2*c0*d0*pow(k2,2)*r1*s2 - 4*c0*c1*pow(k2,2)*s0*s2 - 
   2*pow(c0,2)*pow(k2,2)*s1*s2 + 2*d0*d1*e0*pow(j2,2)*t2 + 
   pow(d0,2)*e1*pow(j2,2)*t2 - 2*c1*d0*f0*pow(j2,2)*t2 - 
   2*c0*d1*f0*pow(j2,2)*t2 - 2*c0*d0*f1*pow(j2,2)*t2 + 
   2*c0*c1*g0*pow(j2,2)*t2 + pow(c0,2)*g1*pow(j2,2)*t2 - 
   2*d0*d1*e0*i2*o2*t2 - pow(d0,2)*e1*i2*o2*t2 + 2*c1*d0*f0*i2*o2*t2 + 
   2*c0*d1*f0*i2*o2*t2 + 2*c0*d0*f1*i2*o2*t2 - 2*c0*c1*g0*i2*o2*t2 - 
   pow(c0,2)*g1*i2*o2*t2 + 2*d1*e0*l2*n0*o2*t2 + 2*d0*e1*l2*n0*o2*t2 - 
   2*c1*f0*l2*n0*o2*t2 - 2*c0*f1*l2*n0*o2*t2 - 2*c1*d0*m2*n0*o2*t2 - 
   2*c0*d1*m2*n0*o2*t2 + 2*c1*c2*pow(n0,2)*o2*t2 + 2*d0*e0*l2*n1*o2*t2 - 
   2*c0*f0*l2*n1*o2*t2 - 2*c0*d0*m2*n1*o2*t2 + 4*c0*c2*n0*n1*o2*t2 + 
   4*c0*c1*n0*n2*o2*t2 + 2*pow(c0,2)*n1*n2*o2*t2 - 2*d1*e0*j2*n0*q2*t2 - 
   2*d0*e1*j2*n0*q2*t2 + 2*c1*f0*j2*n0*q2*t2 + 2*c0*f1*j2*n0*q2*t2 - 
   2*d0*e0*j2*n1*q2*t2 + 2*c0*f0*j2*n1*q2*t2 + 
   e1*pow(n0,2)*pow(q2,2)*t2 + 2*e0*n0*n1*pow(q2,2)*t2 + 
   2*d1*f0*j2*l2*r0*t2 + 2*d0*f1*j2*l2*r0*t2 - 2*c1*g0*j2*l2*r0*t2 - 
   2*c0*g1*j2*l2*r0*t2 - 4*d0*d1*j2*m2*r0*t2 + 2*c2*d1*j2*n0*r0*t2 + 
   2*c1*d2*j2*n0*r0*t2 + 2*c2*d0*j2*n1*r0*t2 + 2*c0*d2*j2*n1*r0*t2 + 
   2*c1*d0*j2*n2*r0*t2 + 2*c0*d1*j2*n2*r0*t2 - 2*d1*f0*i2*q2*r0*t2 - 
   2*d0*f1*i2*q2*r0*t2 + 2*c1*g0*i2*q2*r0*t2 + 2*c0*g1*i2*q2*r0*t2 + 
   2*f1*l2*n0*q2*r0*t2 + 2*d1*m2*n0*q2*r0*t2 + 2*f0*l2*n1*q2*r0*t2 + 
   2*d0*m2*n1*q2*r0*t2 - 4*c2*n0*n1*q2*r0*t2 - 4*c1*n0*n2*q2*r0*t2 - 
   4*c0*n1*n2*q2*r0*t2 + 2*d1*d2*i2*pow(r0,2)*t2 + 
   g1*pow(l2,2)*pow(r0,2)*t2 - 2*d2*l2*n1*pow(r0,2)*t2 - 
   2*d1*l2*n2*pow(r0,2)*t2 + 2*d0*f0*j2*l2*r1*t2 - 2*c0*g0*j2*l2*r1*t2 - 
   2*pow(d0,2)*j2*m2*r1*t2 + 2*c2*d0*j2*n0*r1*t2 + 2*c0*d2*j2*n0*r1*t2 + 
   2*c0*d0*j2*n2*r1*t2 - 2*d0*f0*i2*q2*r1*t2 + 2*c0*g0*i2*q2*r1*t2 + 
   2*f0*l2*n0*q2*r1*t2 + 2*d0*m2*n0*q2*r1*t2 - 2*c2*pow(n0,2)*q2*r1*t2 - 
   4*c0*n0*n2*q2*r1*t2 + 4*d0*d2*i2*r0*r1*t2 + 2*g0*pow(l2,2)*r0*r1*t2 - 
   4*d2*l2*n0*r0*r1*t2 - 4*d0*l2*n2*r0*r1*t2 + 2*c1*d0*j2*n0*r2*t2 + 
   2*c0*d1*j2*n0*r2*t2 + 2*c0*d0*j2*n1*r2*t2 - 2*c1*pow(n0,2)*q2*r2*t2 - 
   4*c0*n0*n1*q2*r2*t2 + 4*d0*d1*i2*r0*r2*t2 - 4*d1*l2*n0*r0*r2*t2 - 
   4*d0*l2*n1*r0*r2*t2 + 2*pow(d0,2)*i2*r1*r2*t2 - 4*d0*l2*n0*r1*r2*t2 - 
   2*d1*e0*j2*l2*s0*t2 - 2*d0*e1*j2*l2*s0*t2 + 2*c1*f0*j2*l2*s0*t2 + 
   2*c0*f1*j2*l2*s0*t2 + 2*c1*d0*j2*m2*s0*t2 + 2*c0*d1*j2*m2*s0*t2 - 
   4*c1*c2*j2*n0*s0*t2 - 4*c0*c2*j2*n1*s0*t2 - 4*c0*c1*j2*n2*s0*t2 + 
   2*d1*e0*i2*q2*s0*t2 + 2*d0*e1*i2*q2*s0*t2 - 2*c1*f0*i2*q2*s0*t2 - 
   2*c0*f1*i2*q2*s0*t2 - 2*e1*l2*n0*q2*s0*t2 + 2*c1*m2*n0*q2*s0*t2 - 
   2*e0*l2*n1*q2*s0*t2 + 2*c0*m2*n1*q2*s0*t2 - 2*c2*d1*i2*r0*s0*t2 - 
   2*c1*d2*i2*r0*s0*t2 - 2*f1*pow(l2,2)*r0*s0*t2 + 2*d1*l2*m2*r0*s0*t2 + 
   2*c2*l2*n1*r0*s0*t2 + 2*c1*l2*n2*r0*s0*t2 - 2*c2*d0*i2*r1*s0*t2 - 
   2*c0*d2*i2*r1*s0*t2 - 2*f0*pow(l2,2)*r1*s0*t2 + 2*d0*l2*m2*r1*s0*t2 + 
   2*c2*l2*n0*r1*s0*t2 + 2*c0*l2*n2*r1*s0*t2 - 2*c1*d0*i2*r2*s0*t2 - 
   2*c0*d1*i2*r2*s0*t2 + 2*c1*l2*n0*r2*s0*t2 + 2*c0*l2*n1*r2*s0*t2 + 
   2*c1*c2*i2*pow(s0,2)*t2 + e1*pow(l2,2)*pow(s0,2)*t2 - 
   2*c1*l2*m2*pow(s0,2)*t2 - 2*d0*e0*j2*l2*s1*t2 + 2*c0*f0*j2*l2*s1*t2 + 
   2*c0*d0*j2*m2*s1*t2 - 4*c0*c2*j2*n0*s1*t2 - 2*pow(c0,2)*j2*n2*s1*t2 + 
   2*d0*e0*i2*q2*s1*t2 - 2*c0*f0*i2*q2*s1*t2 - 2*e0*l2*n0*q2*s1*t2 + 
   2*c0*m2*n0*q2*s1*t2 - 2*c2*d0*i2*r0*s1*t2 - 2*c0*d2*i2*r0*s1*t2 - 
   2*f0*pow(l2,2)*r0*s1*t2 + 2*d0*l2*m2*r0*s1*t2 + 2*c2*l2*n0*r0*s1*t2 + 
   2*c0*l2*n2*r0*s1*t2 - 2*c0*d0*i2*r2*s1*t2 + 2*c0*l2*n0*r2*s1*t2 + 
   4*c0*c2*i2*s0*s1*t2 + 2*e0*pow(l2,2)*s0*s1*t2 - 4*c0*l2*m2*s0*s1*t2 - 
   4*c0*c1*j2*n0*s2*t2 - 2*pow(c0,2)*j2*n1*s2*t2 - 2*c1*d0*i2*r0*s2*t2 - 
   2*c0*d1*i2*r0*s2*t2 + 2*c1*l2*n0*r0*s2*t2 + 2*c0*l2*n1*r0*s2*t2 - 
   2*c0*d0*i2*r1*s2*t2 + 2*c0*l2*n0*r1*s2*t2 + 4*c0*c1*i2*s0*s2*t2 + 
   2*pow(c0,2)*i2*s1*s2*t2 - 4*f0*f1*k2*l2*o2*u0 + 2*e1*g0*k2*l2*o2*u0 + 
   2*e0*g1*k2*l2*o2*u0 + 2*d1*f0*k2*m2*o2*u0 + 2*d0*f1*k2*m2*o2*u0 - 
   2*c1*g0*k2*m2*o2*u0 - 2*c0*g1*k2*m2*o2*u0 - 2*d2*e1*k2*n0*o2*u0 - 
   2*d1*e2*k2*n0*o2*u0 + 2*c2*f1*k2*n0*o2*u0 + 2*c1*f2*k2*n0*o2*u0 - 
   2*d2*e0*k2*n1*o2*u0 - 2*d0*e2*k2*n1*o2*u0 + 2*c2*f0*k2*n1*o2*u0 + 
   2*c0*f2*k2*n1*o2*u0 - 2*d1*e0*k2*n2*o2*u0 - 2*d0*e1*k2*n2*o2*u0 + 
   2*c1*f0*k2*n2*o2*u0 + 2*c0*f1*k2*n2*o2*u0 + 4*f0*f1*j2*l2*p2*u0 - 
   2*e1*g0*j2*l2*p2*u0 - 2*e0*g1*j2*l2*p2*u0 - 2*d1*f0*j2*m2*p2*u0 - 
   2*d0*f1*j2*m2*p2*u0 + 2*c1*g0*j2*m2*p2*u0 + 2*c0*g1*j2*m2*p2*u0 + 
   2*d2*e1*j2*n0*p2*u0 + 2*d1*e2*j2*n0*p2*u0 - 2*c2*f1*j2*n0*p2*u0 - 
   2*c1*f2*j2*n0*p2*u0 + 2*d2*e0*j2*n1*p2*u0 + 2*d0*e2*j2*n1*p2*u0 - 
   2*c2*f0*j2*n1*p2*u0 - 2*c0*f2*j2*n1*p2*u0 + 2*d1*e0*j2*n2*p2*u0 + 
   2*d0*e1*j2*n2*p2*u0 - 2*c1*f0*j2*n2*p2*u0 - 2*c0*f1*j2*n2*p2*u0 + 
   4*f0*f1*j2*k2*q2*u0 - 2*e1*g0*j2*k2*q2*u0 - 2*e0*g1*j2*k2*q2*u0 - 
   4*f0*f1*i2*p2*q2*u0 + 2*e1*g0*i2*p2*q2*u0 + 2*e0*g1*i2*p2*q2*u0 + 
   4*f1*m2*n0*p2*q2*u0 + 4*f0*m2*n1*p2*q2*u0 - 4*e2*n0*n1*p2*q2*u0 - 
   4*e1*n0*n2*p2*q2*u0 - 4*e0*n1*n2*p2*q2*u0 - 2*d2*f1*j2*k2*r0*u0 - 
   2*d1*f2*j2*k2*r0*u0 + 2*c2*g1*j2*k2*r0*u0 + 2*c1*g2*j2*k2*r0*u0 + 
   2*d2*f1*i2*p2*r0*u0 + 2*d1*f2*i2*p2*r0*u0 - 2*c2*g1*i2*p2*r0*u0 - 
   2*c1*g2*i2*p2*r0*u0 + 2*g1*l2*m2*p2*r0*u0 - 2*f2*l2*n1*p2*r0*u0 - 
   2*d2*m2*n1*p2*r0*u0 - 2*f1*l2*n2*p2*r0*u0 - 2*d1*m2*n2*p2*r0*u0 + 
   4*c2*n1*n2*p2*r0*u0 + 2*c1*pow(n2,2)*p2*r0*u0 + 2*g1*k2*m2*q2*r0*u0 - 
   2*f2*k2*n1*q2*r0*u0 - 2*f1*k2*n2*q2*r0*u0 - 2*d2*f0*j2*k2*r1*u0 - 
   2*d0*f2*j2*k2*r1*u0 + 2*c2*g0*j2*k2*r1*u0 + 2*c0*g2*j2*k2*r1*u0 + 
   2*d2*f0*i2*p2*r1*u0 + 2*d0*f2*i2*p2*r1*u0 - 2*c2*g0*i2*p2*r1*u0 - 
   2*c0*g2*i2*p2*r1*u0 + 2*g0*l2*m2*p2*r1*u0 - 2*f2*l2*n0*p2*r1*u0 - 
   2*d2*m2*n0*p2*r1*u0 - 2*f0*l2*n2*p2*r1*u0 - 2*d0*m2*n2*p2*r1*u0 + 
   4*c2*n0*n2*p2*r1*u0 + 2*c0*pow(n2,2)*p2*r1*u0 + 2*g0*k2*m2*q2*r1*u0 - 
   2*f2*k2*n0*q2*r1*u0 - 2*f0*k2*n2*q2*r1*u0 - 4*g2*k2*l2*r0*r1*u0 + 
   4*d2*k2*n2*r0*r1*u0 - 2*d1*f0*j2*k2*r2*u0 - 2*d0*f1*j2*k2*r2*u0 + 
   2*c1*g0*j2*k2*r2*u0 + 2*c0*g1*j2*k2*r2*u0 + 2*d1*f0*i2*p2*r2*u0 + 
   2*d0*f1*i2*p2*r2*u0 - 2*c1*g0*i2*p2*r2*u0 - 2*c0*g1*i2*p2*r2*u0 - 
   2*f1*l2*n0*p2*r2*u0 - 2*d1*m2*n0*p2*r2*u0 - 2*f0*l2*n1*p2*r2*u0 - 
   2*d0*m2*n1*p2*r2*u0 + 4*c2*n0*n1*p2*r2*u0 + 4*c1*n0*n2*p2*r2*u0 + 
   4*c0*n1*n2*p2*r2*u0 - 2*f1*k2*n0*q2*r2*u0 - 2*f0*k2*n1*q2*r2*u0 - 
   4*g1*k2*l2*r0*r2*u0 + 4*d2*k2*n1*r0*r2*u0 + 4*d1*k2*n2*r0*r2*u0 - 
   4*g0*k2*l2*r1*r2*u0 + 4*d2*k2*n0*r1*r2*u0 + 4*d0*k2*n2*r1*r2*u0 + 
   2*d1*k2*n0*pow(r2,2)*u0 + 2*d0*k2*n1*pow(r2,2)*u0 + 
   2*d2*e1*j2*k2*s0*u0 + 2*d1*e2*j2*k2*s0*u0 - 2*c2*f1*j2*k2*s0*u0 - 
   2*c1*f2*j2*k2*s0*u0 - 2*d2*e1*i2*p2*s0*u0 - 2*d1*e2*i2*p2*s0*u0 + 
   2*c2*f1*i2*p2*s0*u0 + 2*c1*f2*i2*p2*s0*u0 - 2*f1*l2*m2*p2*s0*u0 + 
   2*d1*pow(m2,2)*p2*s0*u0 + 2*e2*l2*n1*p2*s0*u0 - 2*c2*m2*n1*p2*s0*u0 + 
   2*e1*l2*n2*p2*s0*u0 - 2*c1*m2*n2*p2*s0*u0 - 2*f1*k2*m2*q2*s0*u0 + 
   2*e2*k2*n1*q2*s0*u0 + 2*e1*k2*n2*q2*s0*u0 + 4*f2*k2*l2*r1*s0*u0 - 
   2*d2*k2*m2*r1*s0*u0 - 2*c2*k2*n2*r1*s0*u0 + 4*f1*k2*l2*r2*s0*u0 - 
   2*d1*k2*m2*r2*s0*u0 - 2*c2*k2*n1*r2*s0*u0 - 2*c1*k2*n2*r2*s0*u0 + 
   2*d2*e0*j2*k2*s1*u0 + 2*d0*e2*j2*k2*s1*u0 - 2*c2*f0*j2*k2*s1*u0 - 
   2*c0*f2*j2*k2*s1*u0 - 2*d2*e0*i2*p2*s1*u0 - 2*d0*e2*i2*p2*s1*u0 + 
   2*c2*f0*i2*p2*s1*u0 + 2*c0*f2*i2*p2*s1*u0 - 2*f0*l2*m2*p2*s1*u0 + 
   2*d0*pow(m2,2)*p2*s1*u0 + 2*e2*l2*n0*p2*s1*u0 - 2*c2*m2*n0*p2*s1*u0 + 
   2*e0*l2*n2*p2*s1*u0 - 2*c0*m2*n2*p2*s1*u0 - 2*f0*k2*m2*q2*s1*u0 + 
   2*e2*k2*n0*q2*s1*u0 + 2*e0*k2*n2*q2*s1*u0 + 4*f2*k2*l2*r0*s1*u0 - 
   2*d2*k2*m2*r0*s1*u0 - 2*c2*k2*n2*r0*s1*u0 + 4*f0*k2*l2*r2*s1*u0 - 
   2*d0*k2*m2*r2*s1*u0 - 2*c2*k2*n0*r2*s1*u0 - 2*c0*k2*n2*r2*s1*u0 - 
   4*e2*k2*l2*s0*s1*u0 + 4*c2*k2*m2*s0*s1*u0 + 2*d1*e0*j2*k2*s2*u0 + 
   2*d0*e1*j2*k2*s2*u0 - 2*c1*f0*j2*k2*s2*u0 - 2*c0*f1*j2*k2*s2*u0 - 
   2*d1*e0*i2*p2*s2*u0 - 2*d0*e1*i2*p2*s2*u0 + 2*c1*f0*i2*p2*s2*u0 + 
   2*c0*f1*i2*p2*s2*u0 + 2*e1*l2*n0*p2*s2*u0 - 2*c1*m2*n0*p2*s2*u0 + 
   2*e0*l2*n1*p2*s2*u0 - 2*c0*m2*n1*p2*s2*u0 + 2*e1*k2*n0*q2*s2*u0 + 
   2*e0*k2*n1*q2*s2*u0 + 4*f1*k2*l2*r0*s2*u0 - 2*d1*k2*m2*r0*s2*u0 - 
   2*c2*k2*n1*r0*s2*u0 - 2*c1*k2*n2*r0*s2*u0 + 4*f0*k2*l2*r1*s2*u0 - 
   2*d0*k2*m2*r1*s2*u0 - 2*c2*k2*n0*r1*s2*u0 - 2*c0*k2*n2*r1*s2*u0 - 
   2*c1*k2*n0*r2*s2*u0 - 2*c0*k2*n1*r2*s2*u0 - 4*e1*k2*l2*s0*s2*u0 + 
   4*c1*k2*m2*s0*s2*u0 - 4*e0*k2*l2*s1*s2*u0 + 4*c0*k2*m2*s1*s2*u0 - 
   2*f1*f2*pow(j2,2)*pow(u0,2) + e2*g1*pow(j2,2)*pow(u0,2) + 
   e1*g2*pow(j2,2)*pow(u0,2) + 2*f1*f2*i2*o2*pow(u0,2) - 
   e2*g1*i2*o2*pow(u0,2) - e1*g2*i2*o2*pow(u0,2) + 
   g1*pow(m2,2)*o2*pow(u0,2) - 2*f2*m2*n1*o2*pow(u0,2) - 
   2*f1*m2*n2*o2*pow(u0,2) + 2*e2*n1*n2*o2*pow(u0,2) + 
   e1*pow(n2,2)*o2*pow(u0,2) - 2*g2*j2*m2*r1*pow(u0,2) + 
   2*f2*j2*n2*r1*pow(u0,2) - 2*g1*j2*m2*r2*pow(u0,2) + 
   2*f2*j2*n1*r2*pow(u0,2) + 2*f1*j2*n2*r2*pow(u0,2) + 
   2*g2*i2*r1*r2*pow(u0,2) - 2*pow(n2,2)*r1*r2*pow(u0,2) + 
   g1*i2*pow(r2,2)*pow(u0,2) - 2*n1*n2*pow(r2,2)*pow(u0,2) + 
   2*f2*j2*m2*s1*pow(u0,2) - 2*e2*j2*n2*s1*pow(u0,2) - 
   2*f2*i2*r2*s1*pow(u0,2) + 2*m2*n2*r2*s1*pow(u0,2) + 
   2*f1*j2*m2*s2*pow(u0,2) - 2*e2*j2*n1*s2*pow(u0,2) - 
   2*e1*j2*n2*s2*pow(u0,2) - 2*f2*i2*r1*s2*pow(u0,2) + 
   2*m2*n2*r1*s2*pow(u0,2) - 2*f1*i2*r2*s2*pow(u0,2) + 
   2*m2*n1*r2*s2*pow(u0,2) + 2*e2*i2*s1*s2*pow(u0,2) - 
   2*pow(m2,2)*s1*s2*pow(u0,2) + e1*i2*pow(s2,2)*pow(u0,2) - 
   2*pow(f0,2)*k2*l2*o2*u1 + 2*e0*g0*k2*l2*o2*u1 + 2*d0*f0*k2*m2*o2*u1 - 
   2*c0*g0*k2*m2*o2*u1 - 2*d2*e0*k2*n0*o2*u1 - 2*d0*e2*k2*n0*o2*u1 + 
   2*c2*f0*k2*n0*o2*u1 + 2*c0*f2*k2*n0*o2*u1 - 2*d0*e0*k2*n2*o2*u1 + 
   2*c0*f0*k2*n2*o2*u1 + 2*pow(f0,2)*j2*l2*p2*u1 - 2*e0*g0*j2*l2*p2*u1 - 
   2*d0*f0*j2*m2*p2*u1 + 2*c0*g0*j2*m2*p2*u1 + 2*d2*e0*j2*n0*p2*u1 + 
   2*d0*e2*j2*n0*p2*u1 - 2*c2*f0*j2*n0*p2*u1 - 2*c0*f2*j2*n0*p2*u1 + 
   2*d0*e0*j2*n2*p2*u1 - 2*c0*f0*j2*n2*p2*u1 + 2*pow(f0,2)*j2*k2*q2*u1 - 
   2*e0*g0*j2*k2*q2*u1 - 2*pow(f0,2)*i2*p2*q2*u1 + 2*e0*g0*i2*p2*q2*u1 + 
   4*f0*m2*n0*p2*q2*u1 - 2*e2*pow(n0,2)*p2*q2*u1 - 4*e0*n0*n2*p2*q2*u1 - 
   2*d2*f0*j2*k2*r0*u1 - 2*d0*f2*j2*k2*r0*u1 + 2*c2*g0*j2*k2*r0*u1 + 
   2*c0*g2*j2*k2*r0*u1 + 2*d2*f0*i2*p2*r0*u1 + 2*d0*f2*i2*p2*r0*u1 - 
   2*c2*g0*i2*p2*r0*u1 - 2*c0*g2*i2*p2*r0*u1 + 2*g0*l2*m2*p2*r0*u1 - 
   2*f2*l2*n0*p2*r0*u1 - 2*d2*m2*n0*p2*r0*u1 - 2*f0*l2*n2*p2*r0*u1 - 
   2*d0*m2*n2*p2*r0*u1 + 4*c2*n0*n2*p2*r0*u1 + 2*c0*pow(n2,2)*p2*r0*u1 + 
   2*g0*k2*m2*q2*r0*u1 - 2*f2*k2*n0*q2*r0*u1 - 2*f0*k2*n2*q2*r0*u1 - 
   2*g2*k2*l2*pow(r0,2)*u1 + 2*d2*k2*n2*pow(r0,2)*u1 - 
   2*d0*f0*j2*k2*r2*u1 + 2*c0*g0*j2*k2*r2*u1 + 2*d0*f0*i2*p2*r2*u1 - 
   2*c0*g0*i2*p2*r2*u1 - 2*f0*l2*n0*p2*r2*u1 - 2*d0*m2*n0*p2*r2*u1 + 
   2*c2*pow(n0,2)*p2*r2*u1 + 4*c0*n0*n2*p2*r2*u1 - 2*f0*k2*n0*q2*r2*u1 - 
   4*g0*k2*l2*r0*r2*u1 + 4*d2*k2*n0*r0*r2*u1 + 4*d0*k2*n2*r0*r2*u1 + 
   2*d0*k2*n0*pow(r2,2)*u1 + 2*d2*e0*j2*k2*s0*u1 + 2*d0*e2*j2*k2*s0*u1 - 
   2*c2*f0*j2*k2*s0*u1 - 2*c0*f2*j2*k2*s0*u1 - 2*d2*e0*i2*p2*s0*u1 - 
   2*d0*e2*i2*p2*s0*u1 + 2*c2*f0*i2*p2*s0*u1 + 2*c0*f2*i2*p2*s0*u1 - 
   2*f0*l2*m2*p2*s0*u1 + 2*d0*pow(m2,2)*p2*s0*u1 + 2*e2*l2*n0*p2*s0*u1 - 
   2*c2*m2*n0*p2*s0*u1 + 2*e0*l2*n2*p2*s0*u1 - 2*c0*m2*n2*p2*s0*u1 - 
   2*f0*k2*m2*q2*s0*u1 + 2*e2*k2*n0*q2*s0*u1 + 2*e0*k2*n2*q2*s0*u1 + 
   4*f2*k2*l2*r0*s0*u1 - 2*d2*k2*m2*r0*s0*u1 - 2*c2*k2*n2*r0*s0*u1 + 
   4*f0*k2*l2*r2*s0*u1 - 2*d0*k2*m2*r2*s0*u1 - 2*c2*k2*n0*r2*s0*u1 - 
   2*c0*k2*n2*r2*s0*u1 - 2*e2*k2*l2*pow(s0,2)*u1 + 
   2*c2*k2*m2*pow(s0,2)*u1 + 2*d0*e0*j2*k2*s2*u1 - 2*c0*f0*j2*k2*s2*u1 - 
   2*d0*e0*i2*p2*s2*u1 + 2*c0*f0*i2*p2*s2*u1 + 2*e0*l2*n0*p2*s2*u1 - 
   2*c0*m2*n0*p2*s2*u1 + 2*e0*k2*n0*q2*s2*u1 + 4*f0*k2*l2*r0*s2*u1 - 
   2*d0*k2*m2*r0*s2*u1 - 2*c2*k2*n0*r0*s2*u1 - 2*c0*k2*n2*r0*s2*u1 - 
   2*c0*k2*n0*r2*s2*u1 - 4*e0*k2*l2*s0*s2*u1 + 4*c0*k2*m2*s0*s2*u1 - 
   4*f0*f2*pow(j2,2)*u0*u1 + 2*e2*g0*pow(j2,2)*u0*u1 + 
   2*e0*g2*pow(j2,2)*u0*u1 + 4*f0*f2*i2*o2*u0*u1 - 2*e2*g0*i2*o2*u0*u1 - 
   2*e0*g2*i2*o2*u0*u1 + 2*g0*pow(m2,2)*o2*u0*u1 - 4*f2*m2*n0*o2*u0*u1 - 
   4*f0*m2*n2*o2*u0*u1 + 4*e2*n0*n2*o2*u0*u1 + 2*e0*pow(n2,2)*o2*u0*u1 - 
   4*g2*j2*m2*r0*u0*u1 + 4*f2*j2*n2*r0*u0*u1 - 4*g0*j2*m2*r2*u0*u1 + 
   4*f2*j2*n0*r2*u0*u1 + 4*f0*j2*n2*r2*u0*u1 + 4*g2*i2*r0*r2*u0*u1 - 
   4*pow(n2,2)*r0*r2*u0*u1 + 2*g0*i2*pow(r2,2)*u0*u1 - 
   4*n0*n2*pow(r2,2)*u0*u1 + 4*f2*j2*m2*s0*u0*u1 - 4*e2*j2*n2*s0*u0*u1 - 
   4*f2*i2*r2*s0*u0*u1 + 4*m2*n2*r2*s0*u0*u1 + 4*f0*j2*m2*s2*u0*u1 - 
   4*e2*j2*n0*s2*u0*u1 - 4*e0*j2*n2*s2*u0*u1 - 4*f2*i2*r0*s2*u0*u1 + 
   4*m2*n2*r0*s2*u0*u1 - 4*f0*i2*r2*s2*u0*u1 + 4*m2*n0*r2*s2*u0*u1 + 
   4*e2*i2*s0*s2*u0*u1 - 4*pow(m2,2)*s0*s2*u0*u1 + 
   2*e0*i2*pow(s2,2)*u0*u1 - 2*d1*e0*k2*n0*o2*u2 - 2*d0*e1*k2*n0*o2*u2 + 
   2*c1*f0*k2*n0*o2*u2 + 2*c0*f1*k2*n0*o2*u2 - 2*d0*e0*k2*n1*o2*u2 + 
   2*c0*f0*k2*n1*o2*u2 + 2*d1*e0*j2*n0*p2*u2 + 2*d0*e1*j2*n0*p2*u2 - 
   2*c1*f0*j2*n0*p2*u2 - 2*c0*f1*j2*n0*p2*u2 + 2*d0*e0*j2*n1*p2*u2 - 
   2*c0*f0*j2*n1*p2*u2 - 2*e1*pow(n0,2)*p2*q2*u2 - 4*e0*n0*n1*p2*q2*u2 - 
   2*d1*f0*j2*k2*r0*u2 - 2*d0*f1*j2*k2*r0*u2 + 2*c1*g0*j2*k2*r0*u2 + 
   2*c0*g1*j2*k2*r0*u2 + 2*d1*f0*i2*p2*r0*u2 + 2*d0*f1*i2*p2*r0*u2 - 
   2*c1*g0*i2*p2*r0*u2 - 2*c0*g1*i2*p2*r0*u2 - 2*f1*l2*n0*p2*r0*u2 - 
   2*d1*m2*n0*p2*r0*u2 - 2*f0*l2*n1*p2*r0*u2 - 2*d0*m2*n1*p2*r0*u2 + 
   4*c2*n0*n1*p2*r0*u2 + 4*c1*n0*n2*p2*r0*u2 + 4*c0*n1*n2*p2*r0*u2 - 
   2*f1*k2*n0*q2*r0*u2 - 2*f0*k2*n1*q2*r0*u2 - 2*g1*k2*l2*pow(r0,2)*u2 + 
   2*d2*k2*n1*pow(r0,2)*u2 + 2*d1*k2*n2*pow(r0,2)*u2 - 
   2*d0*f0*j2*k2*r1*u2 + 2*c0*g0*j2*k2*r1*u2 + 2*d0*f0*i2*p2*r1*u2 - 
   2*c0*g0*i2*p2*r1*u2 - 2*f0*l2*n0*p2*r1*u2 - 2*d0*m2*n0*p2*r1*u2 + 
   2*c2*pow(n0,2)*p2*r1*u2 + 4*c0*n0*n2*p2*r1*u2 - 2*f0*k2*n0*q2*r1*u2 - 
   4*g0*k2*l2*r0*r1*u2 + 4*d2*k2*n0*r0*r1*u2 + 4*d0*k2*n2*r0*r1*u2 + 
   2*c1*pow(n0,2)*p2*r2*u2 + 4*c0*n0*n1*p2*r2*u2 + 4*d1*k2*n0*r0*r2*u2 + 
   4*d0*k2*n1*r0*r2*u2 + 4*d0*k2*n0*r1*r2*u2 + 2*d1*e0*j2*k2*s0*u2 + 
   2*d0*e1*j2*k2*s0*u2 - 2*c1*f0*j2*k2*s0*u2 - 2*c0*f1*j2*k2*s0*u2 - 
   2*d1*e0*i2*p2*s0*u2 - 2*d0*e1*i2*p2*s0*u2 + 2*c1*f0*i2*p2*s0*u2 + 
   2*c0*f1*i2*p2*s0*u2 + 2*e1*l2*n0*p2*s0*u2 - 2*c1*m2*n0*p2*s0*u2 + 
   2*e0*l2*n1*p2*s0*u2 - 2*c0*m2*n1*p2*s0*u2 + 2*e1*k2*n0*q2*s0*u2 + 
   2*e0*k2*n1*q2*s0*u2 + 4*f1*k2*l2*r0*s0*u2 - 2*d1*k2*m2*r0*s0*u2 - 
   2*c2*k2*n1*r0*s0*u2 - 2*c1*k2*n2*r0*s0*u2 + 4*f0*k2*l2*r1*s0*u2 - 
   2*d0*k2*m2*r1*s0*u2 - 2*c2*k2*n0*r1*s0*u2 - 2*c0*k2*n2*r1*s0*u2 - 
   2*c1*k2*n0*r2*s0*u2 - 2*c0*k2*n1*r2*s0*u2 - 2*e1*k2*l2*pow(s0,2)*u2 + 
   2*c1*k2*m2*pow(s0,2)*u2 + 2*d0*e0*j2*k2*s1*u2 - 2*c0*f0*j2*k2*s1*u2 - 
   2*d0*e0*i2*p2*s1*u2 + 2*c0*f0*i2*p2*s1*u2 + 2*e0*l2*n0*p2*s1*u2 - 
   2*c0*m2*n0*p2*s1*u2 + 2*e0*k2*n0*q2*s1*u2 + 4*f0*k2*l2*r0*s1*u2 - 
   2*d0*k2*m2*r0*s1*u2 - 2*c2*k2*n0*r0*s1*u2 - 2*c0*k2*n2*r0*s1*u2 - 
   2*c0*k2*n0*r2*s1*u2 - 4*e0*k2*l2*s0*s1*u2 + 4*c0*k2*m2*s0*s1*u2 - 
   2*c1*k2*n0*r0*s2*u2 - 2*c0*k2*n1*r0*s2*u2 - 2*c0*k2*n0*r1*s2*u2 - 
   4*f0*f1*pow(j2,2)*u0*u2 + 2*e1*g0*pow(j2,2)*u0*u2 + 
   2*e0*g1*pow(j2,2)*u0*u2 + 4*f0*f1*i2*o2*u0*u2 - 2*e1*g0*i2*o2*u0*u2 - 
   2*e0*g1*i2*o2*u0*u2 - 4*f1*m2*n0*o2*u0*u2 - 4*f0*m2*n1*o2*u0*u2 + 
   4*e2*n0*n1*o2*u0*u2 + 4*e1*n0*n2*o2*u0*u2 + 4*e0*n1*n2*o2*u0*u2 - 
   4*g1*j2*m2*r0*u0*u2 + 4*f2*j2*n1*r0*u0*u2 + 4*f1*j2*n2*r0*u0*u2 - 
   4*g0*j2*m2*r1*u0*u2 + 4*f2*j2*n0*r1*u0*u2 + 4*f0*j2*n2*r1*u0*u2 + 
   4*g2*i2*r0*r1*u0*u2 - 4*pow(n2,2)*r0*r1*u0*u2 + 4*f1*j2*n0*r2*u0*u2 + 
   4*f0*j2*n1*r2*u0*u2 + 4*g1*i2*r0*r2*u0*u2 - 8*n1*n2*r0*r2*u0*u2 + 
   4*g0*i2*r1*r2*u0*u2 - 8*n0*n2*r1*r2*u0*u2 - 4*n0*n1*pow(r2,2)*u0*u2 + 
   4*f1*j2*m2*s0*u0*u2 - 4*e2*j2*n1*s0*u0*u2 - 4*e1*j2*n2*s0*u0*u2 - 
   4*f2*i2*r1*s0*u0*u2 + 4*m2*n2*r1*s0*u0*u2 - 4*f1*i2*r2*s0*u0*u2 + 
   4*m2*n1*r2*s0*u0*u2 + 4*f0*j2*m2*s1*u0*u2 - 4*e2*j2*n0*s1*u0*u2 - 
   4*e0*j2*n2*s1*u0*u2 - 4*f2*i2*r0*s1*u0*u2 + 4*m2*n2*r0*s1*u0*u2 - 
   4*f0*i2*r2*s1*u0*u2 + 4*m2*n0*r2*s1*u0*u2 + 4*e2*i2*s0*s1*u0*u2 - 
   4*pow(m2,2)*s0*s1*u0*u2 - 4*e1*j2*n0*s2*u0*u2 - 4*e0*j2*n1*s2*u0*u2 - 
   4*f1*i2*r0*s2*u0*u2 + 4*m2*n1*r0*s2*u0*u2 - 4*f0*i2*r1*s2*u0*u2 + 
   4*m2*n0*r1*s2*u0*u2 + 4*e1*i2*s0*s2*u0*u2 + 4*e0*i2*s1*s2*u0*u2 - 
   2*pow(f0,2)*pow(j2,2)*u1*u2 + 2*e0*g0*pow(j2,2)*u1*u2 + 
   2*pow(f0,2)*i2*o2*u1*u2 - 2*e0*g0*i2*o2*u1*u2 - 4*f0*m2*n0*o2*u1*u2 + 
   2*e2*pow(n0,2)*o2*u1*u2 + 4*e0*n0*n2*o2*u1*u2 - 4*g0*j2*m2*r0*u1*u2 + 
   4*f2*j2*n0*r0*u1*u2 + 4*f0*j2*n2*r0*u1*u2 + 2*g2*i2*pow(r0,2)*u1*u2 - 
   2*pow(n2,2)*pow(r0,2)*u1*u2 + 4*f0*j2*n0*r2*u1*u2 + 
   4*g0*i2*r0*r2*u1*u2 - 8*n0*n2*r0*r2*u1*u2 - 
   2*pow(n0,2)*pow(r2,2)*u1*u2 + 4*f0*j2*m2*s0*u1*u2 - 
   4*e2*j2*n0*s0*u1*u2 - 4*e0*j2*n2*s0*u1*u2 - 4*f2*i2*r0*s0*u1*u2 + 
   4*m2*n2*r0*s0*u1*u2 - 4*f0*i2*r2*s0*u1*u2 + 4*m2*n0*r2*s0*u1*u2 + 
   2*e2*i2*pow(s0,2)*u1*u2 - 2*pow(m2,2)*pow(s0,2)*u1*u2 - 
   4*e0*j2*n0*s2*u1*u2 - 4*f0*i2*r0*s2*u1*u2 + 4*m2*n0*r0*s2*u1*u2 + 
   4*e0*i2*s0*s2*u1*u2 + e1*pow(n0,2)*o2*pow(u2,2) + 
   2*e0*n0*n1*o2*pow(u2,2) + 2*f1*j2*n0*r0*pow(u2,2) + 
   2*f0*j2*n1*r0*pow(u2,2) + g1*i2*pow(r0,2)*pow(u2,2) - 
   2*n1*n2*pow(r0,2)*pow(u2,2) + 2*f0*j2*n0*r1*pow(u2,2) + 
   2*g0*i2*r0*r1*pow(u2,2) - 4*n0*n2*r0*r1*pow(u2,2) - 
   4*n0*n1*r0*r2*pow(u2,2) - 2*pow(n0,2)*r1*r2*pow(u2,2) - 
   2*e1*j2*n0*s0*pow(u2,2) - 2*e0*j2*n1*s0*pow(u2,2) - 
   2*f1*i2*r0*s0*pow(u2,2) + 2*m2*n1*r0*s0*pow(u2,2) - 
   2*f0*i2*r1*s0*pow(u2,2) + 2*m2*n0*r1*s0*pow(u2,2) + 
   e1*i2*pow(s0,2)*pow(u2,2) - 2*e0*j2*n0*s1*pow(u2,2) - 
   2*f0*i2*r0*s1*pow(u2,2) + 2*m2*n0*r0*s1*pow(u2,2) + 
   2*e0*i2*s0*s1*pow(u2,2) + 2*d1*f0*k2*l2*o2*v0 + 2*d0*f1*k2*l2*o2*v0 - 
   2*c1*g0*k2*l2*o2*v0 - 2*c0*g1*k2*l2*o2*v0 - 4*d0*d1*k2*m2*o2*v0 + 
   2*c2*d1*k2*n0*o2*v0 + 2*c1*d2*k2*n0*o2*v0 + 2*c2*d0*k2*n1*o2*v0 + 
   2*c0*d2*k2*n1*o2*v0 + 2*c1*d0*k2*n2*o2*v0 + 2*c0*d1*k2*n2*o2*v0 - 
   2*d1*f0*j2*l2*p2*v0 - 2*d0*f1*j2*l2*p2*v0 + 2*c1*g0*j2*l2*p2*v0 + 
   2*c0*g1*j2*l2*p2*v0 + 4*d0*d1*j2*m2*p2*v0 - 2*c2*d1*j2*n0*p2*v0 - 
   2*c1*d2*j2*n0*p2*v0 - 2*c2*d0*j2*n1*p2*v0 - 2*c0*d2*j2*n1*p2*v0 - 
   2*c1*d0*j2*n2*p2*v0 - 2*c0*d1*j2*n2*p2*v0 - 2*d1*f0*j2*k2*q2*v0 - 
   2*d0*f1*j2*k2*q2*v0 + 2*c1*g0*j2*k2*q2*v0 + 2*c0*g1*j2*k2*q2*v0 + 
   2*d1*f0*i2*p2*q2*v0 + 2*d0*f1*i2*p2*q2*v0 - 2*c1*g0*i2*p2*q2*v0 - 
   2*c0*g1*i2*p2*q2*v0 - 2*f1*l2*n0*p2*q2*v0 - 2*d1*m2*n0*p2*q2*v0 - 
   2*f0*l2*n1*p2*q2*v0 - 2*d0*m2*n1*p2*q2*v0 + 4*c2*n0*n1*p2*q2*v0 + 
   4*c1*n0*n2*p2*q2*v0 + 4*c0*n1*n2*p2*q2*v0 + 2*f1*k2*n0*pow(q2,2)*v0 + 
   2*f0*k2*n1*pow(q2,2)*v0 + 4*d1*d2*j2*k2*r0*v0 - 4*d1*d2*i2*p2*r0*v0 - 
   2*g1*pow(l2,2)*p2*r0*v0 + 4*d2*l2*n1*p2*r0*v0 + 4*d1*l2*n2*p2*r0*v0 + 
   2*g1*k2*l2*q2*r0*v0 - 2*d2*k2*n1*q2*r0*v0 - 2*d1*k2*n2*q2*r0*v0 + 
   4*d0*d2*j2*k2*r1*v0 - 4*d0*d2*i2*p2*r1*v0 - 2*g0*pow(l2,2)*p2*r1*v0 + 
   4*d2*l2*n0*p2*r1*v0 + 4*d0*l2*n2*p2*r1*v0 + 2*g0*k2*l2*q2*r1*v0 - 
   2*d2*k2*n0*q2*r1*v0 - 2*d0*k2*n2*q2*r1*v0 + 4*d0*d1*j2*k2*r2*v0 - 
   4*d0*d1*i2*p2*r2*v0 + 4*d1*l2*n0*p2*r2*v0 + 4*d0*l2*n1*p2*r2*v0 - 
   2*d1*k2*n0*q2*r2*v0 - 2*d0*k2*n1*q2*r2*v0 - 2*c2*d1*j2*k2*s0*v0 - 
   2*c1*d2*j2*k2*s0*v0 + 2*c2*d1*i2*p2*s0*v0 + 2*c1*d2*i2*p2*s0*v0 + 
   2*f1*pow(l2,2)*p2*s0*v0 - 2*d1*l2*m2*p2*s0*v0 - 2*c2*l2*n1*p2*s0*v0 - 
   2*c1*l2*n2*p2*s0*v0 - 2*f1*k2*l2*q2*s0*v0 + 4*d1*k2*m2*q2*s0*v0 - 
   2*c2*k2*n1*q2*s0*v0 - 2*c1*k2*n2*q2*s0*v0 - 2*d2*k2*l2*r1*s0*v0 - 
   2*d1*k2*l2*r2*s0*v0 - 2*c2*d0*j2*k2*s1*v0 - 2*c0*d2*j2*k2*s1*v0 + 
   2*c2*d0*i2*p2*s1*v0 + 2*c0*d2*i2*p2*s1*v0 + 2*f0*pow(l2,2)*p2*s1*v0 - 
   2*d0*l2*m2*p2*s1*v0 - 2*c2*l2*n0*p2*s1*v0 - 2*c0*l2*n2*p2*s1*v0 - 
   2*f0*k2*l2*q2*s1*v0 + 4*d0*k2*m2*q2*s1*v0 - 2*c2*k2*n0*q2*s1*v0 - 
   2*c0*k2*n2*q2*s1*v0 - 2*d2*k2*l2*r0*s1*v0 - 2*d0*k2*l2*r2*s1*v0 + 
   4*c2*k2*l2*s0*s1*v0 - 2*c1*d0*j2*k2*s2*v0 - 2*c0*d1*j2*k2*s2*v0 + 
   2*c1*d0*i2*p2*s2*v0 + 2*c0*d1*i2*p2*s2*v0 - 2*c1*l2*n0*p2*s2*v0 - 
   2*c0*l2*n1*p2*s2*v0 - 2*c1*k2*n0*q2*s2*v0 - 2*c0*k2*n1*q2*s2*v0 - 
   2*d1*k2*l2*r0*s2*v0 - 2*d0*k2*l2*r1*s2*v0 + 4*c1*k2*l2*s0*s2*v0 + 
   4*c0*k2*l2*s1*s2*v0 + 2*d2*f1*pow(j2,2)*u0*v0 + 
   2*d1*f2*pow(j2,2)*u0*v0 - 2*c2*g1*pow(j2,2)*u0*v0 - 
   2*c1*g2*pow(j2,2)*u0*v0 - 2*d2*f1*i2*o2*u0*v0 - 2*d1*f2*i2*o2*u0*v0 + 
   2*c2*g1*i2*o2*u0*v0 + 2*c1*g2*i2*o2*u0*v0 - 2*g1*l2*m2*o2*u0*v0 + 
   2*f2*l2*n1*o2*u0*v0 + 2*d2*m2*n1*o2*u0*v0 + 2*f1*l2*n2*o2*u0*v0 + 
   2*d1*m2*n2*o2*u0*v0 - 4*c2*n1*n2*o2*u0*v0 - 2*c1*pow(n2,2)*o2*u0*v0 + 
   2*g1*j2*m2*q2*u0*v0 - 2*f2*j2*n1*q2*u0*v0 - 2*f1*j2*n2*q2*u0*v0 + 
   2*g2*j2*l2*r1*u0*v0 - 2*d2*j2*n2*r1*u0*v0 - 2*g2*i2*q2*r1*u0*v0 + 
   2*pow(n2,2)*q2*r1*u0*v0 + 2*g1*j2*l2*r2*u0*v0 - 2*d2*j2*n1*r2*u0*v0 - 
   2*d1*j2*n2*r2*u0*v0 - 2*g1*i2*q2*r2*u0*v0 + 4*n1*n2*q2*r2*u0*v0 - 
   2*f2*j2*l2*s1*u0*v0 - 2*d2*j2*m2*s1*u0*v0 + 4*c2*j2*n2*s1*u0*v0 + 
   2*f2*i2*q2*s1*u0*v0 - 2*m2*n2*q2*s1*u0*v0 + 2*d2*i2*r2*s1*u0*v0 - 
   2*l2*n2*r2*s1*u0*v0 - 2*f1*j2*l2*s2*u0*v0 - 2*d1*j2*m2*s2*u0*v0 + 
   4*c2*j2*n1*s2*u0*v0 + 4*c1*j2*n2*s2*u0*v0 + 2*f1*i2*q2*s2*u0*v0 - 
   2*m2*n1*q2*s2*u0*v0 + 2*d2*i2*r1*s2*u0*v0 - 2*l2*n2*r1*s2*u0*v0 + 
   2*d1*i2*r2*s2*u0*v0 - 2*l2*n1*r2*s2*u0*v0 - 4*c2*i2*s1*s2*u0*v0 + 
   4*l2*m2*s1*s2*u0*v0 - 2*c1*i2*pow(s2,2)*u0*v0 + 
   2*d2*f0*pow(j2,2)*u1*v0 + 2*d0*f2*pow(j2,2)*u1*v0 - 
   2*c2*g0*pow(j2,2)*u1*v0 - 2*c0*g2*pow(j2,2)*u1*v0 - 
   2*d2*f0*i2*o2*u1*v0 - 2*d0*f2*i2*o2*u1*v0 + 2*c2*g0*i2*o2*u1*v0 + 
   2*c0*g2*i2*o2*u1*v0 - 2*g0*l2*m2*o2*u1*v0 + 2*f2*l2*n0*o2*u1*v0 + 
   2*d2*m2*n0*o2*u1*v0 + 2*f0*l2*n2*o2*u1*v0 + 2*d0*m2*n2*o2*u1*v0 - 
   4*c2*n0*n2*o2*u1*v0 - 2*c0*pow(n2,2)*o2*u1*v0 + 2*g0*j2*m2*q2*u1*v0 - 
   2*f2*j2*n0*q2*u1*v0 - 2*f0*j2*n2*q2*u1*v0 + 2*g2*j2*l2*r0*u1*v0 - 
   2*d2*j2*n2*r0*u1*v0 - 2*g2*i2*q2*r0*u1*v0 + 2*pow(n2,2)*q2*r0*u1*v0 + 
   2*g0*j2*l2*r2*u1*v0 - 2*d2*j2*n0*r2*u1*v0 - 2*d0*j2*n2*r2*u1*v0 - 
   2*g0*i2*q2*r2*u1*v0 + 4*n0*n2*q2*r2*u1*v0 - 2*f2*j2*l2*s0*u1*v0 - 
   2*d2*j2*m2*s0*u1*v0 + 4*c2*j2*n2*s0*u1*v0 + 2*f2*i2*q2*s0*u1*v0 - 
   2*m2*n2*q2*s0*u1*v0 + 2*d2*i2*r2*s0*u1*v0 - 2*l2*n2*r2*s0*u1*v0 - 
   2*f0*j2*l2*s2*u1*v0 - 2*d0*j2*m2*s2*u1*v0 + 4*c2*j2*n0*s2*u1*v0 + 
   4*c0*j2*n2*s2*u1*v0 + 2*f0*i2*q2*s2*u1*v0 - 2*m2*n0*q2*s2*u1*v0 + 
   2*d2*i2*r0*s2*u1*v0 - 2*l2*n2*r0*s2*u1*v0 + 2*d0*i2*r2*s2*u1*v0 - 
   2*l2*n0*r2*s2*u1*v0 - 4*c2*i2*s0*s2*u1*v0 + 4*l2*m2*s0*s2*u1*v0 - 
   2*c0*i2*pow(s2,2)*u1*v0 + 2*d1*f0*pow(j2,2)*u2*v0 + 
   2*d0*f1*pow(j2,2)*u2*v0 - 2*c1*g0*pow(j2,2)*u2*v0 - 
   2*c0*g1*pow(j2,2)*u2*v0 - 2*d1*f0*i2*o2*u2*v0 - 2*d0*f1*i2*o2*u2*v0 + 
   2*c1*g0*i2*o2*u2*v0 + 2*c0*g1*i2*o2*u2*v0 + 2*f1*l2*n0*o2*u2*v0 + 
   2*d1*m2*n0*o2*u2*v0 + 2*f0*l2*n1*o2*u2*v0 + 2*d0*m2*n1*o2*u2*v0 - 
   4*c2*n0*n1*o2*u2*v0 - 4*c1*n0*n2*o2*u2*v0 - 4*c0*n1*n2*o2*u2*v0 - 
   2*f1*j2*n0*q2*u2*v0 - 2*f0*j2*n1*q2*u2*v0 + 2*g1*j2*l2*r0*u2*v0 - 
   2*d2*j2*n1*r0*u2*v0 - 2*d1*j2*n2*r0*u2*v0 - 2*g1*i2*q2*r0*u2*v0 + 
   4*n1*n2*q2*r0*u2*v0 + 2*g0*j2*l2*r1*u2*v0 - 2*d2*j2*n0*r1*u2*v0 - 
   2*d0*j2*n2*r1*u2*v0 - 2*g0*i2*q2*r1*u2*v0 + 4*n0*n2*q2*r1*u2*v0 - 
   2*d1*j2*n0*r2*u2*v0 - 2*d0*j2*n1*r2*u2*v0 + 4*n0*n1*q2*r2*u2*v0 - 
   2*f1*j2*l2*s0*u2*v0 - 2*d1*j2*m2*s0*u2*v0 + 4*c2*j2*n1*s0*u2*v0 + 
   4*c1*j2*n2*s0*u2*v0 + 2*f1*i2*q2*s0*u2*v0 - 2*m2*n1*q2*s0*u2*v0 + 
   2*d2*i2*r1*s0*u2*v0 - 2*l2*n2*r1*s0*u2*v0 + 2*d1*i2*r2*s0*u2*v0 - 
   2*l2*n1*r2*s0*u2*v0 - 2*f0*j2*l2*s1*u2*v0 - 2*d0*j2*m2*s1*u2*v0 + 
   4*c2*j2*n0*s1*u2*v0 + 4*c0*j2*n2*s1*u2*v0 + 2*f0*i2*q2*s1*u2*v0 - 
   2*m2*n0*q2*s1*u2*v0 + 2*d2*i2*r0*s1*u2*v0 - 2*l2*n2*r0*s1*u2*v0 + 
   2*d0*i2*r2*s1*u2*v0 - 2*l2*n0*r2*s1*u2*v0 - 4*c2*i2*s0*s1*u2*v0 + 
   4*l2*m2*s0*s1*u2*v0 + 4*c1*j2*n0*s2*u2*v0 + 4*c0*j2*n1*s2*u2*v0 + 
   2*d1*i2*r0*s2*u2*v0 - 2*l2*n1*r0*s2*u2*v0 + 2*d0*i2*r1*s2*u2*v0 - 
   2*l2*n0*r1*s2*u2*v0 - 4*c1*i2*s0*s2*u2*v0 - 4*c0*i2*s1*s2*u2*v0 - 
   2*d1*d2*pow(j2,2)*pow(v0,2) + 2*d1*d2*i2*o2*pow(v0,2) + 
   g1*pow(l2,2)*o2*pow(v0,2) - 2*d2*l2*n1*o2*pow(v0,2) - 
   2*d1*l2*n2*o2*pow(v0,2) - 2*g1*j2*l2*q2*pow(v0,2) + 
   2*d2*j2*n1*q2*pow(v0,2) + 2*d1*j2*n2*q2*pow(v0,2) + 
   g1*i2*pow(q2,2)*pow(v0,2) - 2*n1*n2*pow(q2,2)*pow(v0,2) + 
   2*d2*j2*l2*s1*pow(v0,2) - 2*d2*i2*q2*s1*pow(v0,2) + 
   2*l2*n2*q2*s1*pow(v0,2) + 2*d1*j2*l2*s2*pow(v0,2) - 
   2*d1*i2*q2*s2*pow(v0,2) + 2*l2*n1*q2*s2*pow(v0,2) - 
   2*pow(l2,2)*s1*s2*pow(v0,2) + 2*d0*f0*k2*l2*o2*v1 - 
   2*c0*g0*k2*l2*o2*v1 - 2*pow(d0,2)*k2*m2*o2*v1 + 2*c2*d0*k2*n0*o2*v1 + 
   2*c0*d2*k2*n0*o2*v1 + 2*c0*d0*k2*n2*o2*v1 - 2*d0*f0*j2*l2*p2*v1 + 
   2*c0*g0*j2*l2*p2*v1 + 2*pow(d0,2)*j2*m2*p2*v1 - 2*c2*d0*j2*n0*p2*v1 - 
   2*c0*d2*j2*n0*p2*v1 - 2*c0*d0*j2*n2*p2*v1 - 2*d0*f0*j2*k2*q2*v1 + 
   2*c0*g0*j2*k2*q2*v1 + 2*d0*f0*i2*p2*q2*v1 - 2*c0*g0*i2*p2*q2*v1 - 
   2*f0*l2*n0*p2*q2*v1 - 2*d0*m2*n0*p2*q2*v1 + 2*c2*pow(n0,2)*p2*q2*v1 + 
   4*c0*n0*n2*p2*q2*v1 + 2*f0*k2*n0*pow(q2,2)*v1 + 4*d0*d2*j2*k2*r0*v1 - 
   4*d0*d2*i2*p2*r0*v1 - 2*g0*pow(l2,2)*p2*r0*v1 + 4*d2*l2*n0*p2*r0*v1 + 
   4*d0*l2*n2*p2*r0*v1 + 2*g0*k2*l2*q2*r0*v1 - 2*d2*k2*n0*q2*r0*v1 - 
   2*d0*k2*n2*q2*r0*v1 + 2*pow(d0,2)*j2*k2*r2*v1 - 
   2*pow(d0,2)*i2*p2*r2*v1 + 4*d0*l2*n0*p2*r2*v1 - 2*d0*k2*n0*q2*r2*v1 - 
   2*c2*d0*j2*k2*s0*v1 - 2*c0*d2*j2*k2*s0*v1 + 2*c2*d0*i2*p2*s0*v1 + 
   2*c0*d2*i2*p2*s0*v1 + 2*f0*pow(l2,2)*p2*s0*v1 - 2*d0*l2*m2*p2*s0*v1 - 
   2*c2*l2*n0*p2*s0*v1 - 2*c0*l2*n2*p2*s0*v1 - 2*f0*k2*l2*q2*s0*v1 + 
   4*d0*k2*m2*q2*s0*v1 - 2*c2*k2*n0*q2*s0*v1 - 2*c0*k2*n2*q2*s0*v1 - 
   2*d2*k2*l2*r0*s0*v1 - 2*d0*k2*l2*r2*s0*v1 + 2*c2*k2*l2*pow(s0,2)*v1 - 
   2*c0*d0*j2*k2*s2*v1 + 2*c0*d0*i2*p2*s2*v1 - 2*c0*l2*n0*p2*s2*v1 - 
   2*c0*k2*n0*q2*s2*v1 - 2*d0*k2*l2*r0*s2*v1 + 4*c0*k2*l2*s0*s2*v1 + 
   2*d2*f0*pow(j2,2)*u0*v1 + 2*d0*f2*pow(j2,2)*u0*v1 - 
   2*c2*g0*pow(j2,2)*u0*v1 - 2*c0*g2*pow(j2,2)*u0*v1 - 
   2*d2*f0*i2*o2*u0*v1 - 2*d0*f2*i2*o2*u0*v1 + 2*c2*g0*i2*o2*u0*v1 + 
   2*c0*g2*i2*o2*u0*v1 - 2*g0*l2*m2*o2*u0*v1 + 2*f2*l2*n0*o2*u0*v1 + 
   2*d2*m2*n0*o2*u0*v1 + 2*f0*l2*n2*o2*u0*v1 + 2*d0*m2*n2*o2*u0*v1 - 
   4*c2*n0*n2*o2*u0*v1 - 2*c0*pow(n2,2)*o2*u0*v1 + 2*g0*j2*m2*q2*u0*v1 - 
   2*f2*j2*n0*q2*u0*v1 - 2*f0*j2*n2*q2*u0*v1 + 2*g2*j2*l2*r0*u0*v1 - 
   2*d2*j2*n2*r0*u0*v1 - 2*g2*i2*q2*r0*u0*v1 + 2*pow(n2,2)*q2*r0*u0*v1 + 
   2*g0*j2*l2*r2*u0*v1 - 2*d2*j2*n0*r2*u0*v1 - 2*d0*j2*n2*r2*u0*v1 - 
   2*g0*i2*q2*r2*u0*v1 + 4*n0*n2*q2*r2*u0*v1 - 2*f2*j2*l2*s0*u0*v1 - 
   2*d2*j2*m2*s0*u0*v1 + 4*c2*j2*n2*s0*u0*v1 + 2*f2*i2*q2*s0*u0*v1 - 
   2*m2*n2*q2*s0*u0*v1 + 2*d2*i2*r2*s0*u0*v1 - 2*l2*n2*r2*s0*u0*v1 - 
   2*f0*j2*l2*s2*u0*v1 - 2*d0*j2*m2*s2*u0*v1 + 4*c2*j2*n0*s2*u0*v1 + 
   4*c0*j2*n2*s2*u0*v1 + 2*f0*i2*q2*s2*u0*v1 - 2*m2*n0*q2*s2*u0*v1 + 
   2*d2*i2*r0*s2*u0*v1 - 2*l2*n2*r0*s2*u0*v1 + 2*d0*i2*r2*s2*u0*v1 - 
   2*l2*n0*r2*s2*u0*v1 - 4*c2*i2*s0*s2*u0*v1 + 4*l2*m2*s0*s2*u0*v1 - 
   2*c0*i2*pow(s2,2)*u0*v1 + 2*d0*f0*pow(j2,2)*u2*v1 - 
   2*c0*g0*pow(j2,2)*u2*v1 - 2*d0*f0*i2*o2*u2*v1 + 2*c0*g0*i2*o2*u2*v1 + 
   2*f0*l2*n0*o2*u2*v1 + 2*d0*m2*n0*o2*u2*v1 - 2*c2*pow(n0,2)*o2*u2*v1 - 
   4*c0*n0*n2*o2*u2*v1 - 2*f0*j2*n0*q2*u2*v1 + 2*g0*j2*l2*r0*u2*v1 - 
   2*d2*j2*n0*r0*u2*v1 - 2*d0*j2*n2*r0*u2*v1 - 2*g0*i2*q2*r0*u2*v1 + 
   4*n0*n2*q2*r0*u2*v1 - 2*d0*j2*n0*r2*u2*v1 + 2*pow(n0,2)*q2*r2*u2*v1 - 
   2*f0*j2*l2*s0*u2*v1 - 2*d0*j2*m2*s0*u2*v1 + 4*c2*j2*n0*s0*u2*v1 + 
   4*c0*j2*n2*s0*u2*v1 + 2*f0*i2*q2*s0*u2*v1 - 2*m2*n0*q2*s0*u2*v1 + 
   2*d2*i2*r0*s0*u2*v1 - 2*l2*n2*r0*s0*u2*v1 + 2*d0*i2*r2*s0*u2*v1 - 
   2*l2*n0*r2*s0*u2*v1 - 2*c2*i2*pow(s0,2)*u2*v1 + 
   2*l2*m2*pow(s0,2)*u2*v1 + 4*c0*j2*n0*s2*u2*v1 + 2*d0*i2*r0*s2*u2*v1 - 
   2*l2*n0*r0*s2*u2*v1 - 4*c0*i2*s0*s2*u2*v1 - 4*d0*d2*pow(j2,2)*v0*v1 + 
   4*d0*d2*i2*o2*v0*v1 + 2*g0*pow(l2,2)*o2*v0*v1 - 4*d2*l2*n0*o2*v0*v1 - 
   4*d0*l2*n2*o2*v0*v1 - 4*g0*j2*l2*q2*v0*v1 + 4*d2*j2*n0*q2*v0*v1 + 
   4*d0*j2*n2*q2*v0*v1 + 2*g0*i2*pow(q2,2)*v0*v1 - 
   4*n0*n2*pow(q2,2)*v0*v1 + 4*d2*j2*l2*s0*v0*v1 - 4*d2*i2*q2*s0*v0*v1 + 
   4*l2*n2*q2*s0*v0*v1 + 4*d0*j2*l2*s2*v0*v1 - 4*d0*i2*q2*s2*v0*v1 + 
   4*l2*n0*q2*s2*v0*v1 - 4*pow(l2,2)*s0*s2*v0*v1 + 2*c1*d0*k2*n0*o2*v2 + 
   2*c0*d1*k2*n0*o2*v2 + 2*c0*d0*k2*n1*o2*v2 - 2*c1*d0*j2*n0*p2*v2 - 
   2*c0*d1*j2*n0*p2*v2 - 2*c0*d0*j2*n1*p2*v2 + 2*c1*pow(n0,2)*p2*q2*v2 + 
   4*c0*n0*n1*p2*q2*v2 + 4*d0*d1*j2*k2*r0*v2 - 4*d0*d1*i2*p2*r0*v2 + 
   4*d1*l2*n0*p2*r0*v2 + 4*d0*l2*n1*p2*r0*v2 - 2*d1*k2*n0*q2*r0*v2 - 
   2*d0*k2*n1*q2*r0*v2 + 2*pow(d0,2)*j2*k2*r1*v2 - 
   2*pow(d0,2)*i2*p2*r1*v2 + 4*d0*l2*n0*p2*r1*v2 - 2*d0*k2*n0*q2*r1*v2 - 
   2*c1*d0*j2*k2*s0*v2 - 2*c0*d1*j2*k2*s0*v2 + 2*c1*d0*i2*p2*s0*v2 + 
   2*c0*d1*i2*p2*s0*v2 - 2*c1*l2*n0*p2*s0*v2 - 2*c0*l2*n1*p2*s0*v2 - 
   2*c1*k2*n0*q2*s0*v2 - 2*c0*k2*n1*q2*s0*v2 - 2*d1*k2*l2*r0*s0*v2 - 
   2*d0*k2*l2*r1*s0*v2 + 2*c1*k2*l2*pow(s0,2)*v2 - 2*c0*d0*j2*k2*s1*v2 + 
   2*c0*d0*i2*p2*s1*v2 - 2*c0*l2*n0*p2*s1*v2 - 2*c0*k2*n0*q2*s1*v2 - 
   2*d0*k2*l2*r0*s1*v2 + 4*c0*k2*l2*s0*s1*v2 + 2*d1*f0*pow(j2,2)*u0*v2 + 
   2*d0*f1*pow(j2,2)*u0*v2 - 2*c1*g0*pow(j2,2)*u0*v2 - 
   2*c0*g1*pow(j2,2)*u0*v2 - 2*d1*f0*i2*o2*u0*v2 - 2*d0*f1*i2*o2*u0*v2 + 
   2*c1*g0*i2*o2*u0*v2 + 2*c0*g1*i2*o2*u0*v2 + 2*f1*l2*n0*o2*u0*v2 + 
   2*d1*m2*n0*o2*u0*v2 + 2*f0*l2*n1*o2*u0*v2 + 2*d0*m2*n1*o2*u0*v2 - 
   4*c2*n0*n1*o2*u0*v2 - 4*c1*n0*n2*o2*u0*v2 - 4*c0*n1*n2*o2*u0*v2 - 
   2*f1*j2*n0*q2*u0*v2 - 2*f0*j2*n1*q2*u0*v2 + 2*g1*j2*l2*r0*u0*v2 - 
   2*d2*j2*n1*r0*u0*v2 - 2*d1*j2*n2*r0*u0*v2 - 2*g1*i2*q2*r0*u0*v2 + 
   4*n1*n2*q2*r0*u0*v2 + 2*g0*j2*l2*r1*u0*v2 - 2*d2*j2*n0*r1*u0*v2 - 
   2*d0*j2*n2*r1*u0*v2 - 2*g0*i2*q2*r1*u0*v2 + 4*n0*n2*q2*r1*u0*v2 - 
   2*d1*j2*n0*r2*u0*v2 - 2*d0*j2*n1*r2*u0*v2 + 4*n0*n1*q2*r2*u0*v2 - 
   2*f1*j2*l2*s0*u0*v2 - 2*d1*j2*m2*s0*u0*v2 + 4*c2*j2*n1*s0*u0*v2 + 
   4*c1*j2*n2*s0*u0*v2 + 2*f1*i2*q2*s0*u0*v2 - 2*m2*n1*q2*s0*u0*v2 + 
   2*d2*i2*r1*s0*u0*v2 - 2*l2*n2*r1*s0*u0*v2 + 2*d1*i2*r2*s0*u0*v2 - 
   2*l2*n1*r2*s0*u0*v2 - 2*f0*j2*l2*s1*u0*v2 - 2*d0*j2*m2*s1*u0*v2 + 
   4*c2*j2*n0*s1*u0*v2 + 4*c0*j2*n2*s1*u0*v2 + 2*f0*i2*q2*s1*u0*v2 - 
   2*m2*n0*q2*s1*u0*v2 + 2*d2*i2*r0*s1*u0*v2 - 2*l2*n2*r0*s1*u0*v2 + 
   2*d0*i2*r2*s1*u0*v2 - 2*l2*n0*r2*s1*u0*v2 - 4*c2*i2*s0*s1*u0*v2 + 
   4*l2*m2*s0*s1*u0*v2 + 4*c1*j2*n0*s2*u0*v2 + 4*c0*j2*n1*s2*u0*v2 + 
   2*d1*i2*r0*s2*u0*v2 - 2*l2*n1*r0*s2*u0*v2 + 2*d0*i2*r1*s2*u0*v2 - 
   2*l2*n0*r1*s2*u0*v2 - 4*c1*i2*s0*s2*u0*v2 - 4*c0*i2*s1*s2*u0*v2 + 
   2*d0*f0*pow(j2,2)*u1*v2 - 2*c0*g0*pow(j2,2)*u1*v2 - 
   2*d0*f0*i2*o2*u1*v2 + 2*c0*g0*i2*o2*u1*v2 + 2*f0*l2*n0*o2*u1*v2 + 
   2*d0*m2*n0*o2*u1*v2 - 2*c2*pow(n0,2)*o2*u1*v2 - 4*c0*n0*n2*o2*u1*v2 - 
   2*f0*j2*n0*q2*u1*v2 + 2*g0*j2*l2*r0*u1*v2 - 2*d2*j2*n0*r0*u1*v2 - 
   2*d0*j2*n2*r0*u1*v2 - 2*g0*i2*q2*r0*u1*v2 + 4*n0*n2*q2*r0*u1*v2 - 
   2*d0*j2*n0*r2*u1*v2 + 2*pow(n0,2)*q2*r2*u1*v2 - 2*f0*j2*l2*s0*u1*v2 - 
   2*d0*j2*m2*s0*u1*v2 + 4*c2*j2*n0*s0*u1*v2 + 4*c0*j2*n2*s0*u1*v2 + 
   2*f0*i2*q2*s0*u1*v2 - 2*m2*n0*q2*s0*u1*v2 + 2*d2*i2*r0*s0*u1*v2 - 
   2*l2*n2*r0*s0*u1*v2 + 2*d0*i2*r2*s0*u1*v2 - 2*l2*n0*r2*s0*u1*v2 - 
   2*c2*i2*pow(s0,2)*u1*v2 + 2*l2*m2*pow(s0,2)*u1*v2 + 
   4*c0*j2*n0*s2*u1*v2 + 2*d0*i2*r0*s2*u1*v2 - 2*l2*n0*r0*s2*u1*v2 - 
   4*c0*i2*s0*s2*u1*v2 - 2*c1*pow(n0,2)*o2*u2*v2 - 4*c0*n0*n1*o2*u2*v2 - 
   2*d1*j2*n0*r0*u2*v2 - 2*d0*j2*n1*r0*u2*v2 + 4*n0*n1*q2*r0*u2*v2 - 
   2*d0*j2*n0*r1*u2*v2 + 2*pow(n0,2)*q2*r1*u2*v2 + 4*c1*j2*n0*s0*u2*v2 + 
   4*c0*j2*n1*s0*u2*v2 + 2*d1*i2*r0*s0*u2*v2 - 2*l2*n1*r0*s0*u2*v2 + 
   2*d0*i2*r1*s0*u2*v2 - 2*l2*n0*r1*s0*u2*v2 - 2*c1*i2*pow(s0,2)*u2*v2 + 
   4*c0*j2*n0*s1*u2*v2 + 2*d0*i2*r0*s1*u2*v2 - 2*l2*n0*r0*s1*u2*v2 - 
   4*c0*i2*s0*s1*u2*v2 - 4*d0*d1*pow(j2,2)*v0*v2 + 4*d0*d1*i2*o2*v0*v2 - 
   4*d1*l2*n0*o2*v0*v2 - 4*d0*l2*n1*o2*v0*v2 + 4*d1*j2*n0*q2*v0*v2 + 
   4*d0*j2*n1*q2*v0*v2 - 4*n0*n1*pow(q2,2)*v0*v2 + 4*d1*j2*l2*s0*v0*v2 - 
   4*d1*i2*q2*s0*v0*v2 + 4*l2*n1*q2*s0*v0*v2 + 4*d0*j2*l2*s1*v0*v2 - 
   4*d0*i2*q2*s1*v0*v2 + 4*l2*n0*q2*s1*v0*v2 - 4*pow(l2,2)*s0*s1*v0*v2 - 
   2*pow(d0,2)*pow(j2,2)*v1*v2 + 2*pow(d0,2)*i2*o2*v1*v2 - 
   4*d0*l2*n0*o2*v1*v2 + 4*d0*j2*n0*q2*v1*v2 - 
   2*pow(n0,2)*pow(q2,2)*v1*v2 + 4*d0*j2*l2*s0*v1*v2 - 
   4*d0*i2*q2*s0*v1*v2 + 4*l2*n0*q2*s0*v1*v2 - 
   2*pow(l2,2)*pow(s0,2)*v1*v2 - 2*d1*e0*k2*l2*o2*w0 - 
   2*d0*e1*k2*l2*o2*w0 + 2*c1*f0*k2*l2*o2*w0 + 2*c0*f1*k2*l2*o2*w0 + 
   2*c1*d0*k2*m2*o2*w0 + 2*c0*d1*k2*m2*o2*w0 - 4*c1*c2*k2*n0*o2*w0 - 
   4*c0*c2*k2*n1*o2*w0 - 4*c0*c1*k2*n2*o2*w0 + 2*d1*e0*j2*l2*p2*w0 + 
   2*d0*e1*j2*l2*p2*w0 - 2*c1*f0*j2*l2*p2*w0 - 2*c0*f1*j2*l2*p2*w0 - 
   2*c1*d0*j2*m2*p2*w0 - 2*c0*d1*j2*m2*p2*w0 + 4*c1*c2*j2*n0*p2*w0 + 
   4*c0*c2*j2*n1*p2*w0 + 4*c0*c1*j2*n2*p2*w0 + 2*d1*e0*j2*k2*q2*w0 + 
   2*d0*e1*j2*k2*q2*w0 - 2*c1*f0*j2*k2*q2*w0 - 2*c0*f1*j2*k2*q2*w0 - 
   2*d1*e0*i2*p2*q2*w0 - 2*d0*e1*i2*p2*q2*w0 + 2*c1*f0*i2*p2*q2*w0 + 
   2*c0*f1*i2*p2*q2*w0 + 2*e1*l2*n0*p2*q2*w0 - 2*c1*m2*n0*p2*q2*w0 + 
   2*e0*l2*n1*p2*q2*w0 - 2*c0*m2*n1*p2*q2*w0 - 2*e1*k2*n0*pow(q2,2)*w0 - 
   2*e0*k2*n1*pow(q2,2)*w0 - 2*c2*d1*j2*k2*r0*w0 - 2*c1*d2*j2*k2*r0*w0 + 
   2*c2*d1*i2*p2*r0*w0 + 2*c1*d2*i2*p2*r0*w0 + 2*f1*pow(l2,2)*p2*r0*w0 - 
   2*d1*l2*m2*p2*r0*w0 - 2*c2*l2*n1*p2*r0*w0 - 2*c1*l2*n2*p2*r0*w0 - 
   2*f1*k2*l2*q2*r0*w0 - 2*d1*k2*m2*q2*r0*w0 + 4*c2*k2*n1*q2*r0*w0 + 
   4*c1*k2*n2*q2*r0*w0 - 2*c2*d0*j2*k2*r1*w0 - 2*c0*d2*j2*k2*r1*w0 + 
   2*c2*d0*i2*p2*r1*w0 + 2*c0*d2*i2*p2*r1*w0 + 2*f0*pow(l2,2)*p2*r1*w0 - 
   2*d0*l2*m2*p2*r1*w0 - 2*c2*l2*n0*p2*r1*w0 - 2*c0*l2*n2*p2*r1*w0 - 
   2*f0*k2*l2*q2*r1*w0 - 2*d0*k2*m2*q2*r1*w0 + 4*c2*k2*n0*q2*r1*w0 + 
   4*c0*k2*n2*q2*r1*w0 + 4*d2*k2*l2*r0*r1*w0 - 2*c1*d0*j2*k2*r2*w0 - 
   2*c0*d1*j2*k2*r2*w0 + 2*c1*d0*i2*p2*r2*w0 + 2*c0*d1*i2*p2*r2*w0 - 
   2*c1*l2*n0*p2*r2*w0 - 2*c0*l2*n1*p2*r2*w0 + 4*c1*k2*n0*q2*r2*w0 + 
   4*c0*k2*n1*q2*r2*w0 + 4*d1*k2*l2*r0*r2*w0 + 4*d0*k2*l2*r1*r2*w0 + 
   4*c1*c2*j2*k2*s0*w0 - 4*c1*c2*i2*p2*s0*w0 - 2*e1*pow(l2,2)*p2*s0*w0 + 
   4*c1*l2*m2*p2*s0*w0 + 2*e1*k2*l2*q2*s0*w0 - 2*c1*k2*m2*q2*s0*w0 - 
   2*c2*k2*l2*r1*s0*w0 - 2*c1*k2*l2*r2*s0*w0 + 4*c0*c2*j2*k2*s1*w0 - 
   4*c0*c2*i2*p2*s1*w0 - 2*e0*pow(l2,2)*p2*s1*w0 + 4*c0*l2*m2*p2*s1*w0 + 
   2*e0*k2*l2*q2*s1*w0 - 2*c0*k2*m2*q2*s1*w0 - 2*c2*k2*l2*r0*s1*w0 - 
   2*c0*k2*l2*r2*s1*w0 + 4*c0*c1*j2*k2*s2*w0 - 4*c0*c1*i2*p2*s2*w0 - 
   2*c1*k2*l2*r0*s2*w0 - 2*c0*k2*l2*r1*s2*w0 - 2*d2*e1*pow(j2,2)*u0*w0 - 
   2*d1*e2*pow(j2,2)*u0*w0 + 2*c2*f1*pow(j2,2)*u0*w0 + 
   2*c1*f2*pow(j2,2)*u0*w0 + 2*d2*e1*i2*o2*u0*w0 + 2*d1*e2*i2*o2*u0*w0 - 
   2*c2*f1*i2*o2*u0*w0 - 2*c1*f2*i2*o2*u0*w0 + 2*f1*l2*m2*o2*u0*w0 - 
   2*d1*pow(m2,2)*o2*u0*w0 - 2*e2*l2*n1*o2*u0*w0 + 2*c2*m2*n1*o2*u0*w0 - 
   2*e1*l2*n2*o2*u0*w0 + 2*c1*m2*n2*o2*u0*w0 - 2*f1*j2*m2*q2*u0*w0 + 
   2*e2*j2*n1*q2*u0*w0 + 2*e1*j2*n2*q2*u0*w0 - 2*f2*j2*l2*r1*u0*w0 + 
   4*d2*j2*m2*r1*u0*w0 - 2*c2*j2*n2*r1*u0*w0 + 2*f2*i2*q2*r1*u0*w0 - 
   2*m2*n2*q2*r1*u0*w0 - 2*f1*j2*l2*r2*u0*w0 + 4*d1*j2*m2*r2*u0*w0 - 
   2*c2*j2*n1*r2*u0*w0 - 2*c1*j2*n2*r2*u0*w0 + 2*f1*i2*q2*r2*u0*w0 - 
   2*m2*n1*q2*r2*u0*w0 - 4*d2*i2*r1*r2*u0*w0 + 4*l2*n2*r1*r2*u0*w0 - 
   2*d1*i2*pow(r2,2)*u0*w0 + 2*l2*n1*pow(r2,2)*u0*w0 + 
   2*e2*j2*l2*s1*u0*w0 - 2*c2*j2*m2*s1*u0*w0 - 2*e2*i2*q2*s1*u0*w0 + 
   2*pow(m2,2)*q2*s1*u0*w0 + 2*c2*i2*r2*s1*u0*w0 - 2*l2*m2*r2*s1*u0*w0 + 
   2*e1*j2*l2*s2*u0*w0 - 2*c1*j2*m2*s2*u0*w0 - 2*e1*i2*q2*s2*u0*w0 + 
   2*c2*i2*r1*s2*u0*w0 - 2*l2*m2*r1*s2*u0*w0 + 2*c1*i2*r2*s2*u0*w0 - 
   2*d2*e0*pow(j2,2)*u1*w0 - 2*d0*e2*pow(j2,2)*u1*w0 + 
   2*c2*f0*pow(j2,2)*u1*w0 + 2*c0*f2*pow(j2,2)*u1*w0 + 
   2*d2*e0*i2*o2*u1*w0 + 2*d0*e2*i2*o2*u1*w0 - 2*c2*f0*i2*o2*u1*w0 - 
   2*c0*f2*i2*o2*u1*w0 + 2*f0*l2*m2*o2*u1*w0 - 2*d0*pow(m2,2)*o2*u1*w0 - 
   2*e2*l2*n0*o2*u1*w0 + 2*c2*m2*n0*o2*u1*w0 - 2*e0*l2*n2*o2*u1*w0 + 
   2*c0*m2*n2*o2*u1*w0 - 2*f0*j2*m2*q2*u1*w0 + 2*e2*j2*n0*q2*u1*w0 + 
   2*e0*j2*n2*q2*u1*w0 - 2*f2*j2*l2*r0*u1*w0 + 4*d2*j2*m2*r0*u1*w0 - 
   2*c2*j2*n2*r0*u1*w0 + 2*f2*i2*q2*r0*u1*w0 - 2*m2*n2*q2*r0*u1*w0 - 
   2*f0*j2*l2*r2*u1*w0 + 4*d0*j2*m2*r2*u1*w0 - 2*c2*j2*n0*r2*u1*w0 - 
   2*c0*j2*n2*r2*u1*w0 + 2*f0*i2*q2*r2*u1*w0 - 2*m2*n0*q2*r2*u1*w0 - 
   4*d2*i2*r0*r2*u1*w0 + 4*l2*n2*r0*r2*u1*w0 - 2*d0*i2*pow(r2,2)*u1*w0 + 
   2*l2*n0*pow(r2,2)*u1*w0 + 2*e2*j2*l2*s0*u1*w0 - 2*c2*j2*m2*s0*u1*w0 - 
   2*e2*i2*q2*s0*u1*w0 + 2*pow(m2,2)*q2*s0*u1*w0 + 2*c2*i2*r2*s0*u1*w0 - 
   2*l2*m2*r2*s0*u1*w0 + 2*e0*j2*l2*s2*u1*w0 - 2*c0*j2*m2*s2*u1*w0 - 
   2*e0*i2*q2*s2*u1*w0 + 2*c2*i2*r0*s2*u1*w0 - 2*l2*m2*r0*s2*u1*w0 + 
   2*c0*i2*r2*s2*u1*w0 - 2*d1*e0*pow(j2,2)*u2*w0 - 
   2*d0*e1*pow(j2,2)*u2*w0 + 2*c1*f0*pow(j2,2)*u2*w0 + 
   2*c0*f1*pow(j2,2)*u2*w0 + 2*d1*e0*i2*o2*u2*w0 + 2*d0*e1*i2*o2*u2*w0 - 
   2*c1*f0*i2*o2*u2*w0 - 2*c0*f1*i2*o2*u2*w0 - 2*e1*l2*n0*o2*u2*w0 + 
   2*c1*m2*n0*o2*u2*w0 - 2*e0*l2*n1*o2*u2*w0 + 2*c0*m2*n1*o2*u2*w0 + 
   2*e1*j2*n0*q2*u2*w0 + 2*e0*j2*n1*q2*u2*w0 - 2*f1*j2*l2*r0*u2*w0 + 
   4*d1*j2*m2*r0*u2*w0 - 2*c2*j2*n1*r0*u2*w0 - 2*c1*j2*n2*r0*u2*w0 + 
   2*f1*i2*q2*r0*u2*w0 - 2*m2*n1*q2*r0*u2*w0 - 2*f0*j2*l2*r1*u2*w0 + 
   4*d0*j2*m2*r1*u2*w0 - 2*c2*j2*n0*r1*u2*w0 - 2*c0*j2*n2*r1*u2*w0 + 
   2*f0*i2*q2*r1*u2*w0 - 2*m2*n0*q2*r1*u2*w0 - 4*d2*i2*r0*r1*u2*w0 + 
   4*l2*n2*r0*r1*u2*w0 - 2*c1*j2*n0*r2*u2*w0 - 2*c0*j2*n1*r2*u2*w0 - 
   4*d1*i2*r0*r2*u2*w0 + 4*l2*n1*r0*r2*u2*w0 - 4*d0*i2*r1*r2*u2*w0 + 
   4*l2*n0*r1*r2*u2*w0 + 2*e1*j2*l2*s0*u2*w0 - 2*c1*j2*m2*s0*u2*w0 - 
   2*e1*i2*q2*s0*u2*w0 + 2*c2*i2*r1*s0*u2*w0 - 2*l2*m2*r1*s0*u2*w0 + 
   2*c1*i2*r2*s0*u2*w0 + 2*e0*j2*l2*s1*u2*w0 - 2*c0*j2*m2*s1*u2*w0 - 
   2*e0*i2*q2*s1*u2*w0 + 2*c2*i2*r0*s1*u2*w0 - 2*l2*m2*r0*s1*u2*w0 + 
   2*c0*i2*r2*s1*u2*w0 + 2*c1*i2*r0*s2*u2*w0 + 2*c0*i2*r1*s2*u2*w0 + 
   2*c2*d1*pow(j2,2)*v0*w0 + 2*c1*d2*pow(j2,2)*v0*w0 - 
   2*c2*d1*i2*o2*v0*w0 - 2*c1*d2*i2*o2*v0*w0 - 2*f1*pow(l2,2)*o2*v0*w0 + 
   2*d1*l2*m2*o2*v0*w0 + 2*c2*l2*n1*o2*v0*w0 + 2*c1*l2*n2*o2*v0*w0 + 
   4*f1*j2*l2*q2*v0*w0 - 2*d1*j2*m2*q2*v0*w0 - 2*c2*j2*n1*q2*v0*w0 - 
   2*c1*j2*n2*q2*v0*w0 - 2*f1*i2*pow(q2,2)*v0*w0 + 
   2*m2*n1*pow(q2,2)*v0*w0 - 2*d2*j2*l2*r1*v0*w0 + 2*d2*i2*q2*r1*v0*w0 - 
   2*l2*n2*q2*r1*v0*w0 - 2*d1*j2*l2*r2*v0*w0 + 2*d1*i2*q2*r2*v0*w0 - 
   2*l2*n1*q2*r2*v0*w0 - 2*c2*j2*l2*s1*v0*w0 + 2*c2*i2*q2*s1*v0*w0 - 
   2*l2*m2*q2*s1*v0*w0 + 2*pow(l2,2)*r2*s1*v0*w0 - 2*c1*j2*l2*s2*v0*w0 + 
   2*c1*i2*q2*s2*v0*w0 + 2*pow(l2,2)*r1*s2*v0*w0 + 
   2*c2*d0*pow(j2,2)*v1*w0 + 2*c0*d2*pow(j2,2)*v1*w0 - 
   2*c2*d0*i2*o2*v1*w0 - 2*c0*d2*i2*o2*v1*w0 - 2*f0*pow(l2,2)*o2*v1*w0 + 
   2*d0*l2*m2*o2*v1*w0 + 2*c2*l2*n0*o2*v1*w0 + 2*c0*l2*n2*o2*v1*w0 + 
   4*f0*j2*l2*q2*v1*w0 - 2*d0*j2*m2*q2*v1*w0 - 2*c2*j2*n0*q2*v1*w0 - 
   2*c0*j2*n2*q2*v1*w0 - 2*f0*i2*pow(q2,2)*v1*w0 + 
   2*m2*n0*pow(q2,2)*v1*w0 - 2*d2*j2*l2*r0*v1*w0 + 2*d2*i2*q2*r0*v1*w0 - 
   2*l2*n2*q2*r0*v1*w0 - 2*d0*j2*l2*r2*v1*w0 + 2*d0*i2*q2*r2*v1*w0 - 
   2*l2*n0*q2*r2*v1*w0 - 2*c2*j2*l2*s0*v1*w0 + 2*c2*i2*q2*s0*v1*w0 - 
   2*l2*m2*q2*s0*v1*w0 + 2*pow(l2,2)*r2*s0*v1*w0 - 2*c0*j2*l2*s2*v1*w0 + 
   2*c0*i2*q2*s2*v1*w0 + 2*pow(l2,2)*r0*s2*v1*w0 + 
   2*c1*d0*pow(j2,2)*v2*w0 + 2*c0*d1*pow(j2,2)*v2*w0 - 
   2*c1*d0*i2*o2*v2*w0 - 2*c0*d1*i2*o2*v2*w0 + 2*c1*l2*n0*o2*v2*w0 + 
   2*c0*l2*n1*o2*v2*w0 - 2*c1*j2*n0*q2*v2*w0 - 2*c0*j2*n1*q2*v2*w0 - 
   2*d1*j2*l2*r0*v2*w0 + 2*d1*i2*q2*r0*v2*w0 - 2*l2*n1*q2*r0*v2*w0 - 
   2*d0*j2*l2*r1*v2*w0 + 2*d0*i2*q2*r1*v2*w0 - 2*l2*n0*q2*r1*v2*w0 - 
   2*c1*j2*l2*s0*v2*w0 + 2*c1*i2*q2*s0*v2*w0 + 2*pow(l2,2)*r1*s0*v2*w0 - 
   2*c0*j2*l2*s1*v2*w0 + 2*c0*i2*q2*s1*v2*w0 + 2*pow(l2,2)*r0*s1*v2*w0 - 
   2*c1*c2*pow(j2,2)*pow(w0,2) + 2*c1*c2*i2*o2*pow(w0,2) + 
   e1*pow(l2,2)*o2*pow(w0,2) - 2*c1*l2*m2*o2*pow(w0,2) - 
   2*e1*j2*l2*q2*pow(w0,2) + 2*c1*j2*m2*q2*pow(w0,2) + 
   e1*i2*pow(q2,2)*pow(w0,2) + 2*c2*j2*l2*r1*pow(w0,2) - 
   2*c2*i2*q2*r1*pow(w0,2) + 2*l2*m2*q2*r1*pow(w0,2) + 
   2*c1*j2*l2*r2*pow(w0,2) - 2*c1*i2*q2*r2*pow(w0,2) - 
   2*pow(l2,2)*r1*r2*pow(w0,2) - 2*d0*e0*k2*l2*o2*w1 + 
   2*c0*f0*k2*l2*o2*w1 + 2*c0*d0*k2*m2*o2*w1 - 4*c0*c2*k2*n0*o2*w1 - 
   2*pow(c0,2)*k2*n2*o2*w1 + 2*d0*e0*j2*l2*p2*w1 - 2*c0*f0*j2*l2*p2*w1 - 
   2*c0*d0*j2*m2*p2*w1 + 4*c0*c2*j2*n0*p2*w1 + 2*pow(c0,2)*j2*n2*p2*w1 + 
   2*d0*e0*j2*k2*q2*w1 - 2*c0*f0*j2*k2*q2*w1 - 2*d0*e0*i2*p2*q2*w1 + 
   2*c0*f0*i2*p2*q2*w1 + 2*e0*l2*n0*p2*q2*w1 - 2*c0*m2*n0*p2*q2*w1 - 
   2*e0*k2*n0*pow(q2,2)*w1 - 2*c2*d0*j2*k2*r0*w1 - 2*c0*d2*j2*k2*r0*w1 + 
   2*c2*d0*i2*p2*r0*w1 + 2*c0*d2*i2*p2*r0*w1 + 2*f0*pow(l2,2)*p2*r0*w1 - 
   2*d0*l2*m2*p2*r0*w1 - 2*c2*l2*n0*p2*r0*w1 - 2*c0*l2*n2*p2*r0*w1 - 
   2*f0*k2*l2*q2*r0*w1 - 2*d0*k2*m2*q2*r0*w1 + 4*c2*k2*n0*q2*r0*w1 + 
   4*c0*k2*n2*q2*r0*w1 + 2*d2*k2*l2*pow(r0,2)*w1 - 2*c0*d0*j2*k2*r2*w1 + 
   2*c0*d0*i2*p2*r2*w1 - 2*c0*l2*n0*p2*r2*w1 + 4*c0*k2*n0*q2*r2*w1 + 
   4*d0*k2*l2*r0*r2*w1 + 4*c0*c2*j2*k2*s0*w1 - 4*c0*c2*i2*p2*s0*w1 - 
   2*e0*pow(l2,2)*p2*s0*w1 + 4*c0*l2*m2*p2*s0*w1 + 2*e0*k2*l2*q2*s0*w1 - 
   2*c0*k2*m2*q2*s0*w1 - 2*c2*k2*l2*r0*s0*w1 - 2*c0*k2*l2*r2*s0*w1 + 
   2*pow(c0,2)*j2*k2*s2*w1 - 2*pow(c0,2)*i2*p2*s2*w1 - 
   2*c0*k2*l2*r0*s2*w1 - 2*d2*e0*pow(j2,2)*u0*w1 - 
   2*d0*e2*pow(j2,2)*u0*w1 + 2*c2*f0*pow(j2,2)*u0*w1 + 
   2*c0*f2*pow(j2,2)*u0*w1 + 2*d2*e0*i2*o2*u0*w1 + 2*d0*e2*i2*o2*u0*w1 - 
   2*c2*f0*i2*o2*u0*w1 - 2*c0*f2*i2*o2*u0*w1 + 2*f0*l2*m2*o2*u0*w1 - 
   2*d0*pow(m2,2)*o2*u0*w1 - 2*e2*l2*n0*o2*u0*w1 + 2*c2*m2*n0*o2*u0*w1 - 
   2*e0*l2*n2*o2*u0*w1 + 2*c0*m2*n2*o2*u0*w1 - 2*f0*j2*m2*q2*u0*w1 + 
   2*e2*j2*n0*q2*u0*w1 + 2*e0*j2*n2*q2*u0*w1 - 2*f2*j2*l2*r0*u0*w1 + 
   4*d2*j2*m2*r0*u0*w1 - 2*c2*j2*n2*r0*u0*w1 + 2*f2*i2*q2*r0*u0*w1 - 
   2*m2*n2*q2*r0*u0*w1 - 2*f0*j2*l2*r2*u0*w1 + 4*d0*j2*m2*r2*u0*w1 - 
   2*c2*j2*n0*r2*u0*w1 - 2*c0*j2*n2*r2*u0*w1 + 2*f0*i2*q2*r2*u0*w1 - 
   2*m2*n0*q2*r2*u0*w1 - 4*d2*i2*r0*r2*u0*w1 + 4*l2*n2*r0*r2*u0*w1 - 
   2*d0*i2*pow(r2,2)*u0*w1 + 2*l2*n0*pow(r2,2)*u0*w1 + 
   2*e2*j2*l2*s0*u0*w1 - 2*c2*j2*m2*s0*u0*w1 - 2*e2*i2*q2*s0*u0*w1 + 
   2*pow(m2,2)*q2*s0*u0*w1 + 2*c2*i2*r2*s0*u0*w1 - 2*l2*m2*r2*s0*u0*w1 + 
   2*e0*j2*l2*s2*u0*w1 - 2*c0*j2*m2*s2*u0*w1 - 2*e0*i2*q2*s2*u0*w1 + 
   2*c2*i2*r0*s2*u0*w1 - 2*l2*m2*r0*s2*u0*w1 + 2*c0*i2*r2*s2*u0*w1 - 
   2*d0*e0*pow(j2,2)*u2*w1 + 2*c0*f0*pow(j2,2)*u2*w1 + 
   2*d0*e0*i2*o2*u2*w1 - 2*c0*f0*i2*o2*u2*w1 - 2*e0*l2*n0*o2*u2*w1 + 
   2*c0*m2*n0*o2*u2*w1 + 2*e0*j2*n0*q2*u2*w1 - 2*f0*j2*l2*r0*u2*w1 + 
   4*d0*j2*m2*r0*u2*w1 - 2*c2*j2*n0*r0*u2*w1 - 2*c0*j2*n2*r0*u2*w1 + 
   2*f0*i2*q2*r0*u2*w1 - 2*m2*n0*q2*r0*u2*w1 - 2*d2*i2*pow(r0,2)*u2*w1 + 
   2*l2*n2*pow(r0,2)*u2*w1 - 2*c0*j2*n0*r2*u2*w1 - 4*d0*i2*r0*r2*u2*w1 + 
   4*l2*n0*r0*r2*u2*w1 + 2*e0*j2*l2*s0*u2*w1 - 2*c0*j2*m2*s0*u2*w1 - 
   2*e0*i2*q2*s0*u2*w1 + 2*c2*i2*r0*s0*u2*w1 - 2*l2*m2*r0*s0*u2*w1 + 
   2*c0*i2*r2*s0*u2*w1 + 2*c0*i2*r0*s2*u2*w1 + 2*c2*d0*pow(j2,2)*v0*w1 + 
   2*c0*d2*pow(j2,2)*v0*w1 - 2*c2*d0*i2*o2*v0*w1 - 2*c0*d2*i2*o2*v0*w1 - 
   2*f0*pow(l2,2)*o2*v0*w1 + 2*d0*l2*m2*o2*v0*w1 + 2*c2*l2*n0*o2*v0*w1 + 
   2*c0*l2*n2*o2*v0*w1 + 4*f0*j2*l2*q2*v0*w1 - 2*d0*j2*m2*q2*v0*w1 - 
   2*c2*j2*n0*q2*v0*w1 - 2*c0*j2*n2*q2*v0*w1 - 2*f0*i2*pow(q2,2)*v0*w1 + 
   2*m2*n0*pow(q2,2)*v0*w1 - 2*d2*j2*l2*r0*v0*w1 + 2*d2*i2*q2*r0*v0*w1 - 
   2*l2*n2*q2*r0*v0*w1 - 2*d0*j2*l2*r2*v0*w1 + 2*d0*i2*q2*r2*v0*w1 - 
   2*l2*n0*q2*r2*v0*w1 - 2*c2*j2*l2*s0*v0*w1 + 2*c2*i2*q2*s0*v0*w1 - 
   2*l2*m2*q2*s0*v0*w1 + 2*pow(l2,2)*r2*s0*v0*w1 - 2*c0*j2*l2*s2*v0*w1 + 
   2*c0*i2*q2*s2*v0*w1 + 2*pow(l2,2)*r0*s2*v0*w1 + 
   2*c0*d0*pow(j2,2)*v2*w1 - 2*c0*d0*i2*o2*v2*w1 + 2*c0*l2*n0*o2*v2*w1 - 
   2*c0*j2*n0*q2*v2*w1 - 2*d0*j2*l2*r0*v2*w1 + 2*d0*i2*q2*r0*v2*w1 - 
   2*l2*n0*q2*r0*v2*w1 - 2*c0*j2*l2*s0*v2*w1 + 2*c0*i2*q2*s0*v2*w1 + 
   2*pow(l2,2)*r0*s0*v2*w1 - 4*c0*c2*pow(j2,2)*w0*w1 + 
   4*c0*c2*i2*o2*w0*w1 + 2*e0*pow(l2,2)*o2*w0*w1 - 4*c0*l2*m2*o2*w0*w1 - 
   4*e0*j2*l2*q2*w0*w1 + 4*c0*j2*m2*q2*w0*w1 + 2*e0*i2*pow(q2,2)*w0*w1 + 
   4*c2*j2*l2*r0*w0*w1 - 4*c2*i2*q2*r0*w0*w1 + 4*l2*m2*q2*r0*w0*w1 + 
   4*c0*j2*l2*r2*w0*w1 - 4*c0*i2*q2*r2*w0*w1 - 4*pow(l2,2)*r0*r2*w0*w1 - 
   4*c0*c1*k2*n0*o2*w2 - 2*pow(c0,2)*k2*n1*o2*w2 + 4*c0*c1*j2*n0*p2*w2 + 
   2*pow(c0,2)*j2*n1*p2*w2 - 2*c1*d0*j2*k2*r0*w2 - 2*c0*d1*j2*k2*r0*w2 + 
   2*c1*d0*i2*p2*r0*w2 + 2*c0*d1*i2*p2*r0*w2 - 2*c1*l2*n0*p2*r0*w2 - 
   2*c0*l2*n1*p2*r0*w2 + 4*c1*k2*n0*q2*r0*w2 + 4*c0*k2*n1*q2*r0*w2 + 
   2*d1*k2*l2*pow(r0,2)*w2 - 2*c0*d0*j2*k2*r1*w2 + 2*c0*d0*i2*p2*r1*w2 - 
   2*c0*l2*n0*p2*r1*w2 + 4*c0*k2*n0*q2*r1*w2 + 4*d0*k2*l2*r0*r1*w2 + 
   4*c0*c1*j2*k2*s0*w2 - 4*c0*c1*i2*p2*s0*w2 - 2*c1*k2*l2*r0*s0*w2 - 
   2*c0*k2*l2*r1*s0*w2 + 2*pow(c0,2)*j2*k2*s1*w2 - 
   2*pow(c0,2)*i2*p2*s1*w2 - 2*c0*k2*l2*r0*s1*w2 - 
   2*d1*e0*pow(j2,2)*u0*w2 - 2*d0*e1*pow(j2,2)*u0*w2 + 
   2*c1*f0*pow(j2,2)*u0*w2 + 2*c0*f1*pow(j2,2)*u0*w2 + 
   2*d1*e0*i2*o2*u0*w2 + 2*d0*e1*i2*o2*u0*w2 - 2*c1*f0*i2*o2*u0*w2 - 
   2*c0*f1*i2*o2*u0*w2 - 2*e1*l2*n0*o2*u0*w2 + 2*c1*m2*n0*o2*u0*w2 - 
   2*e0*l2*n1*o2*u0*w2 + 2*c0*m2*n1*o2*u0*w2 + 2*e1*j2*n0*q2*u0*w2 + 
   2*e0*j2*n1*q2*u0*w2 - 2*f1*j2*l2*r0*u0*w2 + 4*d1*j2*m2*r0*u0*w2 - 
   2*c2*j2*n1*r0*u0*w2 - 2*c1*j2*n2*r0*u0*w2 + 2*f1*i2*q2*r0*u0*w2 - 
   2*m2*n1*q2*r0*u0*w2 - 2*f0*j2*l2*r1*u0*w2 + 4*d0*j2*m2*r1*u0*w2 - 
   2*c2*j2*n0*r1*u0*w2 - 2*c0*j2*n2*r1*u0*w2 + 2*f0*i2*q2*r1*u0*w2 - 
   2*m2*n0*q2*r1*u0*w2 - 4*d2*i2*r0*r1*u0*w2 + 4*l2*n2*r0*r1*u0*w2 - 
   2*c1*j2*n0*r2*u0*w2 - 2*c0*j2*n1*r2*u0*w2 - 4*d1*i2*r0*r2*u0*w2 + 
   4*l2*n1*r0*r2*u0*w2 - 4*d0*i2*r1*r2*u0*w2 + 4*l2*n0*r1*r2*u0*w2 + 
   2*e1*j2*l2*s0*u0*w2 - 2*c1*j2*m2*s0*u0*w2 - 2*e1*i2*q2*s0*u0*w2 + 
   2*c2*i2*r1*s0*u0*w2 - 2*l2*m2*r1*s0*u0*w2 + 2*c1*i2*r2*s0*u0*w2 + 
   2*e0*j2*l2*s1*u0*w2 - 2*c0*j2*m2*s1*u0*w2 - 2*e0*i2*q2*s1*u0*w2 + 
   2*c2*i2*r0*s1*u0*w2 - 2*l2*m2*r0*s1*u0*w2 + 2*c0*i2*r2*s1*u0*w2 + 
   2*c1*i2*r0*s2*u0*w2 + 2*c0*i2*r1*s2*u0*w2 - 2*d0*e0*pow(j2,2)*u1*w2 + 
   2*c0*f0*pow(j2,2)*u1*w2 + 2*d0*e0*i2*o2*u1*w2 - 2*c0*f0*i2*o2*u1*w2 - 
   2*e0*l2*n0*o2*u1*w2 + 2*c0*m2*n0*o2*u1*w2 + 2*e0*j2*n0*q2*u1*w2 - 
   2*f0*j2*l2*r0*u1*w2 + 4*d0*j2*m2*r0*u1*w2 - 2*c2*j2*n0*r0*u1*w2 - 
   2*c0*j2*n2*r0*u1*w2 + 2*f0*i2*q2*r0*u1*w2 - 2*m2*n0*q2*r0*u1*w2 - 
   2*d2*i2*pow(r0,2)*u1*w2 + 2*l2*n2*pow(r0,2)*u1*w2 - 
   2*c0*j2*n0*r2*u1*w2 - 4*d0*i2*r0*r2*u1*w2 + 4*l2*n0*r0*r2*u1*w2 + 
   2*e0*j2*l2*s0*u1*w2 - 2*c0*j2*m2*s0*u1*w2 - 2*e0*i2*q2*s0*u1*w2 + 
   2*c2*i2*r0*s0*u1*w2 - 2*l2*m2*r0*s0*u1*w2 + 2*c0*i2*r2*s0*u1*w2 + 
   2*c0*i2*r0*s2*u1*w2 - 2*c1*j2*n0*r0*u2*w2 - 2*c0*j2*n1*r0*u2*w2 - 
   2*d1*i2*pow(r0,2)*u2*w2 + 2*l2*n1*pow(r0,2)*u2*w2 - 
   2*c0*j2*n0*r1*u2*w2 - 4*d0*i2*r0*r1*u2*w2 + 4*l2*n0*r0*r1*u2*w2 + 
   2*c1*i2*r0*s0*u2*w2 + 2*c0*i2*r1*s0*u2*w2 + 2*c0*i2*r0*s1*u2*w2 + 
   2*c1*d0*pow(j2,2)*v0*w2 + 2*c0*d1*pow(j2,2)*v0*w2 - 
   2*c1*d0*i2*o2*v0*w2 - 2*c0*d1*i2*o2*v0*w2 + 2*c1*l2*n0*o2*v0*w2 + 
   2*c0*l2*n1*o2*v0*w2 - 2*c1*j2*n0*q2*v0*w2 - 2*c0*j2*n1*q2*v0*w2 - 
   2*d1*j2*l2*r0*v0*w2 + 2*d1*i2*q2*r0*v0*w2 - 2*l2*n1*q2*r0*v0*w2 - 
   2*d0*j2*l2*r1*v0*w2 + 2*d0*i2*q2*r1*v0*w2 - 2*l2*n0*q2*r1*v0*w2 - 
   2*c1*j2*l2*s0*v0*w2 + 2*c1*i2*q2*s0*v0*w2 + 2*pow(l2,2)*r1*s0*v0*w2 - 
   2*c0*j2*l2*s1*v0*w2 + 2*c0*i2*q2*s1*v0*w2 + 2*pow(l2,2)*r0*s1*v0*w2 + 
   2*c0*d0*pow(j2,2)*v1*w2 - 2*c0*d0*i2*o2*v1*w2 + 2*c0*l2*n0*o2*v1*w2 - 
   2*c0*j2*n0*q2*v1*w2 - 2*d0*j2*l2*r0*v1*w2 + 2*d0*i2*q2*r0*v1*w2 - 
   2*l2*n0*q2*r0*v1*w2 - 2*c0*j2*l2*s0*v1*w2 + 2*c0*i2*q2*s0*v1*w2 + 
   2*pow(l2,2)*r0*s0*v1*w2 - 4*c0*c1*pow(j2,2)*w0*w2 + 
   4*c0*c1*i2*o2*w0*w2 + 4*c1*j2*l2*r0*w0*w2 - 4*c1*i2*q2*r0*w0*w2 + 
   4*c0*j2*l2*r1*w0*w2 - 4*c0*i2*q2*r1*w0*w2 - 4*pow(l2,2)*r0*r1*w0*w2 - 
   2*pow(c0,2)*pow(j2,2)*w1*w2 + 2*pow(c0,2)*i2*o2*w1*w2 + 
   4*c0*j2*l2*r0*w1*w2 - 4*c0*i2*q2*r0*w1*w2 - 
   2*pow(l2,2)*pow(r0,2)*w1*w2 + 2*f0*f1*pow(k2,2)*o2*z0 - 
   e1*g0*pow(k2,2)*o2*z0 - e0*g1*pow(k2,2)*o2*z0 - 4*f0*f1*j2*k2*p2*z0 + 
   2*e1*g0*j2*k2*p2*z0 + 2*e0*g1*j2*k2*p2*z0 + 2*f0*f1*i2*pow(p2,2)*z0 - 
   e1*g0*i2*pow(p2,2)*z0 - e0*g1*i2*pow(p2,2)*z0 - 
   2*f1*m2*n0*pow(p2,2)*z0 - 2*f0*m2*n1*pow(p2,2)*z0 + 
   2*e2*n0*n1*pow(p2,2)*z0 + 2*e1*n0*n2*pow(p2,2)*z0 + 
   2*e0*n1*n2*pow(p2,2)*z0 - 2*g1*k2*m2*p2*r0*z0 + 2*f2*k2*n1*p2*r0*z0 + 
   2*f1*k2*n2*p2*r0*z0 - 2*g0*k2*m2*p2*r1*z0 + 2*f2*k2*n0*p2*r1*z0 + 
   2*f0*k2*n2*p2*r1*z0 + 2*g2*pow(k2,2)*r0*r1*z0 + 2*f1*k2*n0*p2*r2*z0 + 
   2*f0*k2*n1*p2*r2*z0 + 2*g1*pow(k2,2)*r0*r2*z0 + 
   2*g0*pow(k2,2)*r1*r2*z0 + 2*f1*k2*m2*p2*s0*z0 - 2*e2*k2*n1*p2*s0*z0 - 
   2*e1*k2*n2*p2*s0*z0 - 2*f2*pow(k2,2)*r1*s0*z0 - 
   2*f1*pow(k2,2)*r2*s0*z0 + 2*f0*k2*m2*p2*s1*z0 - 2*e2*k2*n0*p2*s1*z0 - 
   2*e0*k2*n2*p2*s1*z0 - 2*f2*pow(k2,2)*r0*s1*z0 - 
   2*f0*pow(k2,2)*r2*s1*z0 + 2*e2*pow(k2,2)*s0*s1*z0 - 
   2*e1*k2*n0*p2*s2*z0 - 2*e0*k2*n1*p2*s2*z0 - 2*f1*pow(k2,2)*r0*s2*z0 - 
   2*f0*pow(k2,2)*r1*s2*z0 + 2*e1*pow(k2,2)*s0*s2*z0 + 
   2*e0*pow(k2,2)*s1*s2*z0 + 2*f0*f1*pow(j2,2)*t2*z0 - 
   e1*g0*pow(j2,2)*t2*z0 - e0*g1*pow(j2,2)*t2*z0 - 2*f0*f1*i2*o2*t2*z0 + 
   e1*g0*i2*o2*t2*z0 + e0*g1*i2*o2*t2*z0 + 2*f1*m2*n0*o2*t2*z0 + 
   2*f0*m2*n1*o2*t2*z0 - 2*e2*n0*n1*o2*t2*z0 - 2*e1*n0*n2*o2*t2*z0 - 
   2*e0*n1*n2*o2*t2*z0 + 2*g1*j2*m2*r0*t2*z0 - 2*f2*j2*n1*r0*t2*z0 - 
   2*f1*j2*n2*r0*t2*z0 + 2*g0*j2*m2*r1*t2*z0 - 2*f2*j2*n0*r1*t2*z0 - 
   2*f0*j2*n2*r1*t2*z0 - 2*g2*i2*r0*r1*t2*z0 + 2*pow(n2,2)*r0*r1*t2*z0 - 
   2*f1*j2*n0*r2*t2*z0 - 2*f0*j2*n1*r2*t2*z0 - 2*g1*i2*r0*r2*t2*z0 + 
   4*n1*n2*r0*r2*t2*z0 - 2*g0*i2*r1*r2*t2*z0 + 4*n0*n2*r1*r2*t2*z0 + 
   2*n0*n1*pow(r2,2)*t2*z0 - 2*f1*j2*m2*s0*t2*z0 + 2*e2*j2*n1*s0*t2*z0 + 
   2*e1*j2*n2*s0*t2*z0 + 2*f2*i2*r1*s0*t2*z0 - 2*m2*n2*r1*s0*t2*z0 + 
   2*f1*i2*r2*s0*t2*z0 - 2*m2*n1*r2*s0*t2*z0 - 2*f0*j2*m2*s1*t2*z0 + 
   2*e2*j2*n0*s1*t2*z0 + 2*e0*j2*n2*s1*t2*z0 + 2*f2*i2*r0*s1*t2*z0 - 
   2*m2*n2*r0*s1*t2*z0 + 2*f0*i2*r2*s1*t2*z0 - 2*m2*n0*r2*s1*t2*z0 - 
   2*e2*i2*s0*s1*t2*z0 + 2*pow(m2,2)*s0*s1*t2*z0 + 2*e1*j2*n0*s2*t2*z0 + 
   2*e0*j2*n1*s2*t2*z0 + 2*f1*i2*r0*s2*t2*z0 - 2*m2*n1*r0*s2*t2*z0 + 
   2*f0*i2*r1*s2*t2*z0 - 2*m2*n0*r1*s2*t2*z0 - 2*e1*i2*s0*s2*t2*z0 - 
   2*e0*i2*s1*s2*t2*z0 + 2*g1*k2*m2*o2*v0*z0 - 2*f2*k2*n1*o2*v0*z0 - 
   2*f1*k2*n2*o2*v0*z0 - 2*g1*j2*m2*p2*v0*z0 + 2*f2*j2*n1*p2*v0*z0 + 
   2*f1*j2*n2*p2*v0*z0 - 2*g2*j2*k2*r1*v0*z0 + 2*g2*i2*p2*r1*v0*z0 - 
   2*pow(n2,2)*p2*r1*v0*z0 - 2*g1*j2*k2*r2*v0*z0 + 2*g1*i2*p2*r2*v0*z0 - 
   4*n1*n2*p2*r2*v0*z0 + 2*f2*j2*k2*s1*v0*z0 - 2*f2*i2*p2*s1*v0*z0 + 
   2*m2*n2*p2*s1*v0*z0 + 2*k2*n2*r2*s1*v0*z0 + 2*f1*j2*k2*s2*v0*z0 - 
   2*f1*i2*p2*s2*v0*z0 + 2*m2*n1*p2*s2*v0*z0 + 2*k2*n2*r1*s2*v0*z0 + 
   2*k2*n1*r2*s2*v0*z0 - 4*k2*m2*s1*s2*v0*z0 + 2*g0*k2*m2*o2*v1*z0 - 
   2*f2*k2*n0*o2*v1*z0 - 2*f0*k2*n2*o2*v1*z0 - 2*g0*j2*m2*p2*v1*z0 + 
   2*f2*j2*n0*p2*v1*z0 + 2*f0*j2*n2*p2*v1*z0 - 2*g2*j2*k2*r0*v1*z0 + 
   2*g2*i2*p2*r0*v1*z0 - 2*pow(n2,2)*p2*r0*v1*z0 - 2*g0*j2*k2*r2*v1*z0 + 
   2*g0*i2*p2*r2*v1*z0 - 4*n0*n2*p2*r2*v1*z0 + 2*f2*j2*k2*s0*v1*z0 - 
   2*f2*i2*p2*s0*v1*z0 + 2*m2*n2*p2*s0*v1*z0 + 2*k2*n2*r2*s0*v1*z0 + 
   2*f0*j2*k2*s2*v1*z0 - 2*f0*i2*p2*s2*v1*z0 + 2*m2*n0*p2*s2*v1*z0 + 
   2*k2*n2*r0*s2*v1*z0 + 2*k2*n0*r2*s2*v1*z0 - 4*k2*m2*s0*s2*v1*z0 + 
   2*g2*pow(j2,2)*v0*v1*z0 - 2*g2*i2*o2*v0*v1*z0 + 
   2*pow(n2,2)*o2*v0*v1*z0 - 4*j2*n2*s2*v0*v1*z0 + 
   2*i2*pow(s2,2)*v0*v1*z0 - 2*f1*k2*n0*o2*v2*z0 - 2*f0*k2*n1*o2*v2*z0 + 
   2*f1*j2*n0*p2*v2*z0 + 2*f0*j2*n1*p2*v2*z0 - 2*g1*j2*k2*r0*v2*z0 + 
   2*g1*i2*p2*r0*v2*z0 - 4*n1*n2*p2*r0*v2*z0 - 2*g0*j2*k2*r1*v2*z0 + 
   2*g0*i2*p2*r1*v2*z0 - 4*n0*n2*p2*r1*v2*z0 - 4*n0*n1*p2*r2*v2*z0 + 
   2*f1*j2*k2*s0*v2*z0 - 2*f1*i2*p2*s0*v2*z0 + 2*m2*n1*p2*s0*v2*z0 + 
   2*k2*n2*r1*s0*v2*z0 + 2*k2*n1*r2*s0*v2*z0 + 2*f0*j2*k2*s1*v2*z0 - 
   2*f0*i2*p2*s1*v2*z0 + 2*m2*n0*p2*s1*v2*z0 + 2*k2*n2*r0*s1*v2*z0 + 
   2*k2*n0*r2*s1*v2*z0 - 4*k2*m2*s0*s1*v2*z0 + 2*k2*n1*r0*s2*v2*z0 + 
   2*k2*n0*r1*s2*v2*z0 + 2*g1*pow(j2,2)*v0*v2*z0 - 2*g1*i2*o2*v0*v2*z0 + 
   4*n1*n2*o2*v0*v2*z0 - 4*j2*n2*s1*v0*v2*z0 - 4*j2*n1*s2*v0*v2*z0 + 
   4*i2*s1*s2*v0*v2*z0 + 2*g0*pow(j2,2)*v1*v2*z0 - 2*g0*i2*o2*v1*v2*z0 + 
   4*n0*n2*o2*v1*v2*z0 - 4*j2*n2*s0*v1*v2*z0 - 4*j2*n0*s2*v1*v2*z0 + 
   4*i2*s0*s2*v1*v2*z0 + 2*n0*n1*o2*pow(v2,2)*z0 - 
   2*j2*n1*s0*pow(v2,2)*z0 - 2*j2*n0*s1*pow(v2,2)*z0 + 
   2*i2*s0*s1*pow(v2,2)*z0 - 2*f1*k2*m2*o2*w0*z0 + 2*e2*k2*n1*o2*w0*z0 + 
   2*e1*k2*n2*o2*w0*z0 + 2*f1*j2*m2*p2*w0*z0 - 2*e2*j2*n1*p2*w0*z0 - 
   2*e1*j2*n2*p2*w0*z0 + 2*f2*j2*k2*r1*w0*z0 - 2*f2*i2*p2*r1*w0*z0 + 
   2*m2*n2*p2*r1*w0*z0 + 2*f1*j2*k2*r2*w0*z0 - 2*f1*i2*p2*r2*w0*z0 + 
   2*m2*n1*p2*r2*w0*z0 - 4*k2*n2*r1*r2*w0*z0 - 2*k2*n1*pow(r2,2)*w0*z0 - 
   2*e2*j2*k2*s1*w0*z0 + 2*e2*i2*p2*s1*w0*z0 - 2*pow(m2,2)*p2*s1*w0*z0 + 
   2*k2*m2*r2*s1*w0*z0 - 2*e1*j2*k2*s2*w0*z0 + 2*e1*i2*p2*s2*w0*z0 + 
   2*k2*m2*r1*s2*w0*z0 - 2*f2*pow(j2,2)*v1*w0*z0 + 2*f2*i2*o2*v1*w0*z0 - 
   2*m2*n2*o2*v1*w0*z0 + 2*j2*n2*r2*v1*w0*z0 + 2*j2*m2*s2*v1*w0*z0 - 
   2*i2*r2*s2*v1*w0*z0 - 2*f1*pow(j2,2)*v2*w0*z0 + 2*f1*i2*o2*v2*w0*z0 - 
   2*m2*n1*o2*v2*w0*z0 + 2*j2*n2*r1*v2*w0*z0 + 2*j2*n1*r2*v2*w0*z0 + 
   2*j2*m2*s1*v2*w0*z0 - 2*i2*r2*s1*v2*w0*z0 - 2*i2*r1*s2*v2*w0*z0 - 
   2*f0*k2*m2*o2*w1*z0 + 2*e2*k2*n0*o2*w1*z0 + 2*e0*k2*n2*o2*w1*z0 + 
   2*f0*j2*m2*p2*w1*z0 - 2*e2*j2*n0*p2*w1*z0 - 2*e0*j2*n2*p2*w1*z0 + 
   2*f2*j2*k2*r0*w1*z0 - 2*f2*i2*p2*r0*w1*z0 + 2*m2*n2*p2*r0*w1*z0 + 
   2*f0*j2*k2*r2*w1*z0 - 2*f0*i2*p2*r2*w1*z0 + 2*m2*n0*p2*r2*w1*z0 - 
   4*k2*n2*r0*r2*w1*z0 - 2*k2*n0*pow(r2,2)*w1*z0 - 2*e2*j2*k2*s0*w1*z0 + 
   2*e2*i2*p2*s0*w1*z0 - 2*pow(m2,2)*p2*s0*w1*z0 + 2*k2*m2*r2*s0*w1*z0 - 
   2*e0*j2*k2*s2*w1*z0 + 2*e0*i2*p2*s2*w1*z0 + 2*k2*m2*r0*s2*w1*z0 - 
   2*f2*pow(j2,2)*v0*w1*z0 + 2*f2*i2*o2*v0*w1*z0 - 2*m2*n2*o2*v0*w1*z0 + 
   2*j2*n2*r2*v0*w1*z0 + 2*j2*m2*s2*v0*w1*z0 - 2*i2*r2*s2*v0*w1*z0 - 
   2*f0*pow(j2,2)*v2*w1*z0 + 2*f0*i2*o2*v2*w1*z0 - 2*m2*n0*o2*v2*w1*z0 + 
   2*j2*n2*r0*v2*w1*z0 + 2*j2*n0*r2*v2*w1*z0 + 2*j2*m2*s0*v2*w1*z0 - 
   2*i2*r2*s0*v2*w1*z0 - 2*i2*r0*s2*v2*w1*z0 + 2*e2*pow(j2,2)*w0*w1*z0 - 
   2*e2*i2*o2*w0*w1*z0 + 2*pow(m2,2)*o2*w0*w1*z0 - 4*j2*m2*r2*w0*w1*z0 + 
   2*i2*pow(r2,2)*w0*w1*z0 + 2*e1*k2*n0*o2*w2*z0 + 2*e0*k2*n1*o2*w2*z0 - 
   2*e1*j2*n0*p2*w2*z0 - 2*e0*j2*n1*p2*w2*z0 + 2*f1*j2*k2*r0*w2*z0 - 
   2*f1*i2*p2*r0*w2*z0 + 2*m2*n1*p2*r0*w2*z0 + 2*f0*j2*k2*r1*w2*z0 - 
   2*f0*i2*p2*r1*w2*z0 + 2*m2*n0*p2*r1*w2*z0 - 4*k2*n2*r0*r1*w2*z0 - 
   4*k2*n1*r0*r2*w2*z0 - 4*k2*n0*r1*r2*w2*z0 - 2*e1*j2*k2*s0*w2*z0 + 
   2*e1*i2*p2*s0*w2*z0 + 2*k2*m2*r1*s0*w2*z0 - 2*e0*j2*k2*s1*w2*z0 + 
   2*e0*i2*p2*s1*w2*z0 + 2*k2*m2*r0*s1*w2*z0 - 2*f1*pow(j2,2)*v0*w2*z0 + 
   2*f1*i2*o2*v0*w2*z0 - 2*m2*n1*o2*v0*w2*z0 + 2*j2*n2*r1*v0*w2*z0 + 
   2*j2*n1*r2*v0*w2*z0 + 2*j2*m2*s1*v0*w2*z0 - 2*i2*r2*s1*v0*w2*z0 - 
   2*i2*r1*s2*v0*w2*z0 - 2*f0*pow(j2,2)*v1*w2*z0 + 2*f0*i2*o2*v1*w2*z0 - 
   2*m2*n0*o2*v1*w2*z0 + 2*j2*n2*r0*v1*w2*z0 + 2*j2*n0*r2*v1*w2*z0 + 
   2*j2*m2*s0*v1*w2*z0 - 2*i2*r2*s0*v1*w2*z0 - 2*i2*r0*s2*v1*w2*z0 + 
   2*j2*n1*r0*v2*w2*z0 + 2*j2*n0*r1*v2*w2*z0 - 2*i2*r1*s0*v2*w2*z0 - 
   2*i2*r0*s1*v2*w2*z0 + 2*e1*pow(j2,2)*w0*w2*z0 - 2*e1*i2*o2*w0*w2*z0 - 
   4*j2*m2*r1*w0*w2*z0 + 4*i2*r1*r2*w0*w2*z0 + 2*e0*pow(j2,2)*w1*w2*z0 - 
   2*e0*i2*o2*w1*w2*z0 - 4*j2*m2*r0*w1*w2*z0 + 4*i2*r0*r2*w1*w2*z0 + 
   2*i2*r0*r1*pow(w2,2)*z0 + pow(f0,2)*pow(k2,2)*o2*z1 - 
   e0*g0*pow(k2,2)*o2*z1 - 2*pow(f0,2)*j2*k2*p2*z1 + 
   2*e0*g0*j2*k2*p2*z1 + pow(f0,2)*i2*pow(p2,2)*z1 - 
   e0*g0*i2*pow(p2,2)*z1 - 2*f0*m2*n0*pow(p2,2)*z1 + 
   e2*pow(n0,2)*pow(p2,2)*z1 + 2*e0*n0*n2*pow(p2,2)*z1 - 
   2*g0*k2*m2*p2*r0*z1 + 2*f2*k2*n0*p2*r0*z1 + 2*f0*k2*n2*p2*r0*z1 + 
   g2*pow(k2,2)*pow(r0,2)*z1 + 2*f0*k2*n0*p2*r2*z1 + 
   2*g0*pow(k2,2)*r0*r2*z1 + 2*f0*k2*m2*p2*s0*z1 - 2*e2*k2*n0*p2*s0*z1 - 
   2*e0*k2*n2*p2*s0*z1 - 2*f2*pow(k2,2)*r0*s0*z1 - 
   2*f0*pow(k2,2)*r2*s0*z1 + e2*pow(k2,2)*pow(s0,2)*z1 - 
   2*e0*k2*n0*p2*s2*z1 - 2*f0*pow(k2,2)*r0*s2*z1 + 
   2*e0*pow(k2,2)*s0*s2*z1 + pow(f0,2)*pow(j2,2)*t2*z1 - 
   e0*g0*pow(j2,2)*t2*z1 - pow(f0,2)*i2*o2*t2*z1 + e0*g0*i2*o2*t2*z1 + 
   2*f0*m2*n0*o2*t2*z1 - e2*pow(n0,2)*o2*t2*z1 - 2*e0*n0*n2*o2*t2*z1 + 
   2*g0*j2*m2*r0*t2*z1 - 2*f2*j2*n0*r0*t2*z1 - 2*f0*j2*n2*r0*t2*z1 - 
   g2*i2*pow(r0,2)*t2*z1 + pow(n2,2)*pow(r0,2)*t2*z1 - 
   2*f0*j2*n0*r2*t2*z1 - 2*g0*i2*r0*r2*t2*z1 + 4*n0*n2*r0*r2*t2*z1 + 
   pow(n0,2)*pow(r2,2)*t2*z1 - 2*f0*j2*m2*s0*t2*z1 + 
   2*e2*j2*n0*s0*t2*z1 + 2*e0*j2*n2*s0*t2*z1 + 2*f2*i2*r0*s0*t2*z1 - 
   2*m2*n2*r0*s0*t2*z1 + 2*f0*i2*r2*s0*t2*z1 - 2*m2*n0*r2*s0*t2*z1 - 
   e2*i2*pow(s0,2)*t2*z1 + pow(m2,2)*pow(s0,2)*t2*z1 + 
   2*e0*j2*n0*s2*t2*z1 + 2*f0*i2*r0*s2*t2*z1 - 2*m2*n0*r0*s2*t2*z1 - 
   2*e0*i2*s0*s2*t2*z1 + 2*g0*k2*m2*o2*v0*z1 - 2*f2*k2*n0*o2*v0*z1 - 
   2*f0*k2*n2*o2*v0*z1 - 2*g0*j2*m2*p2*v0*z1 + 2*f2*j2*n0*p2*v0*z1 + 
   2*f0*j2*n2*p2*v0*z1 - 2*g2*j2*k2*r0*v0*z1 + 2*g2*i2*p2*r0*v0*z1 - 
   2*pow(n2,2)*p2*r0*v0*z1 - 2*g0*j2*k2*r2*v0*z1 + 2*g0*i2*p2*r2*v0*z1 - 
   4*n0*n2*p2*r2*v0*z1 + 2*f2*j2*k2*s0*v0*z1 - 2*f2*i2*p2*s0*v0*z1 + 
   2*m2*n2*p2*s0*v0*z1 + 2*k2*n2*r2*s0*v0*z1 + 2*f0*j2*k2*s2*v0*z1 - 
   2*f0*i2*p2*s2*v0*z1 + 2*m2*n0*p2*s2*v0*z1 + 2*k2*n2*r0*s2*v0*z1 + 
   2*k2*n0*r2*s2*v0*z1 - 4*k2*m2*s0*s2*v0*z1 + 
   g2*pow(j2,2)*pow(v0,2)*z1 - g2*i2*o2*pow(v0,2)*z1 + 
   pow(n2,2)*o2*pow(v0,2)*z1 - 2*j2*n2*s2*pow(v0,2)*z1 + 
   i2*pow(s2,2)*pow(v0,2)*z1 - 2*f0*k2*n0*o2*v2*z1 + 
   2*f0*j2*n0*p2*v2*z1 - 2*g0*j2*k2*r0*v2*z1 + 2*g0*i2*p2*r0*v2*z1 - 
   4*n0*n2*p2*r0*v2*z1 - 2*pow(n0,2)*p2*r2*v2*z1 + 2*f0*j2*k2*s0*v2*z1 - 
   2*f0*i2*p2*s0*v2*z1 + 2*m2*n0*p2*s0*v2*z1 + 2*k2*n2*r0*s0*v2*z1 + 
   2*k2*n0*r2*s0*v2*z1 - 2*k2*m2*pow(s0,2)*v2*z1 + 2*k2*n0*r0*s2*v2*z1 + 
   2*g0*pow(j2,2)*v0*v2*z1 - 2*g0*i2*o2*v0*v2*z1 + 4*n0*n2*o2*v0*v2*z1 - 
   4*j2*n2*s0*v0*v2*z1 - 4*j2*n0*s2*v0*v2*z1 + 4*i2*s0*s2*v0*v2*z1 + 
   pow(n0,2)*o2*pow(v2,2)*z1 - 2*j2*n0*s0*pow(v2,2)*z1 + 
   i2*pow(s0,2)*pow(v2,2)*z1 - 2*f0*k2*m2*o2*w0*z1 + 
   2*e2*k2*n0*o2*w0*z1 + 2*e0*k2*n2*o2*w0*z1 + 2*f0*j2*m2*p2*w0*z1 - 
   2*e2*j2*n0*p2*w0*z1 - 2*e0*j2*n2*p2*w0*z1 + 2*f2*j2*k2*r0*w0*z1 - 
   2*f2*i2*p2*r0*w0*z1 + 2*m2*n2*p2*r0*w0*z1 + 2*f0*j2*k2*r2*w0*z1 - 
   2*f0*i2*p2*r2*w0*z1 + 2*m2*n0*p2*r2*w0*z1 - 4*k2*n2*r0*r2*w0*z1 - 
   2*k2*n0*pow(r2,2)*w0*z1 - 2*e2*j2*k2*s0*w0*z1 + 2*e2*i2*p2*s0*w0*z1 - 
   2*pow(m2,2)*p2*s0*w0*z1 + 2*k2*m2*r2*s0*w0*z1 - 2*e0*j2*k2*s2*w0*z1 + 
   2*e0*i2*p2*s2*w0*z1 + 2*k2*m2*r0*s2*w0*z1 - 2*f2*pow(j2,2)*v0*w0*z1 + 
   2*f2*i2*o2*v0*w0*z1 - 2*m2*n2*o2*v0*w0*z1 + 2*j2*n2*r2*v0*w0*z1 + 
   2*j2*m2*s2*v0*w0*z1 - 2*i2*r2*s2*v0*w0*z1 - 2*f0*pow(j2,2)*v2*w0*z1 + 
   2*f0*i2*o2*v2*w0*z1 - 2*m2*n0*o2*v2*w0*z1 + 2*j2*n2*r0*v2*w0*z1 + 
   2*j2*n0*r2*v2*w0*z1 + 2*j2*m2*s0*v2*w0*z1 - 2*i2*r2*s0*v2*w0*z1 - 
   2*i2*r0*s2*v2*w0*z1 + e2*pow(j2,2)*pow(w0,2)*z1 - 
   e2*i2*o2*pow(w0,2)*z1 + pow(m2,2)*o2*pow(w0,2)*z1 - 
   2*j2*m2*r2*pow(w0,2)*z1 + i2*pow(r2,2)*pow(w0,2)*z1 + 
   2*e0*k2*n0*o2*w2*z1 - 2*e0*j2*n0*p2*w2*z1 + 2*f0*j2*k2*r0*w2*z1 - 
   2*f0*i2*p2*r0*w2*z1 + 2*m2*n0*p2*r0*w2*z1 - 2*k2*n2*pow(r0,2)*w2*z1 - 
   4*k2*n0*r0*r2*w2*z1 - 2*e0*j2*k2*s0*w2*z1 + 2*e0*i2*p2*s0*w2*z1 + 
   2*k2*m2*r0*s0*w2*z1 - 2*f0*pow(j2,2)*v0*w2*z1 + 2*f0*i2*o2*v0*w2*z1 - 
   2*m2*n0*o2*v0*w2*z1 + 2*j2*n2*r0*v0*w2*z1 + 2*j2*n0*r2*v0*w2*z1 + 
   2*j2*m2*s0*v0*w2*z1 - 2*i2*r2*s0*v0*w2*z1 - 2*i2*r0*s2*v0*w2*z1 + 
   2*j2*n0*r0*v2*w2*z1 - 2*i2*r0*s0*v2*w2*z1 + 2*e0*pow(j2,2)*w0*w2*z1 - 
   2*e0*i2*o2*w0*w2*z1 - 4*j2*m2*r0*w0*w2*z1 + 4*i2*r0*r2*w0*w2*z1 + 
   i2*pow(r0,2)*pow(w2,2)*z1 + e1*pow(n0,2)*pow(p2,2)*z2 + 
   2*e0*n0*n1*pow(p2,2)*z2 + 2*f1*k2*n0*p2*r0*z2 + 2*f0*k2*n1*p2*r0*z2 + 
   g1*pow(k2,2)*pow(r0,2)*z2 + 2*f0*k2*n0*p2*r1*z2 + 
   2*g0*pow(k2,2)*r0*r1*z2 - 2*e1*k2*n0*p2*s0*z2 - 2*e0*k2*n1*p2*s0*z2 - 
   2*f1*pow(k2,2)*r0*s0*z2 - 2*f0*pow(k2,2)*r1*s0*z2 + 
   e1*pow(k2,2)*pow(s0,2)*z2 - 2*e0*k2*n0*p2*s1*z2 - 
   2*f0*pow(k2,2)*r0*s1*z2 + 2*e0*pow(k2,2)*s0*s1*z2 - 
   e1*pow(n0,2)*o2*t2*z2 - 2*e0*n0*n1*o2*t2*z2 - 2*f1*j2*n0*r0*t2*z2 - 
   2*f0*j2*n1*r0*t2*z2 - g1*i2*pow(r0,2)*t2*z2 + 
   2*n1*n2*pow(r0,2)*t2*z2 - 2*f0*j2*n0*r1*t2*z2 - 2*g0*i2*r0*r1*t2*z2 + 
   4*n0*n2*r0*r1*t2*z2 + 4*n0*n1*r0*r2*t2*z2 + 2*pow(n0,2)*r1*r2*t2*z2 + 
   2*e1*j2*n0*s0*t2*z2 + 2*e0*j2*n1*s0*t2*z2 + 2*f1*i2*r0*s0*t2*z2 - 
   2*m2*n1*r0*s0*t2*z2 + 2*f0*i2*r1*s0*t2*z2 - 2*m2*n0*r1*s0*t2*z2 - 
   e1*i2*pow(s0,2)*t2*z2 + 2*e0*j2*n0*s1*t2*z2 + 2*f0*i2*r0*s1*t2*z2 - 
   2*m2*n0*r0*s1*t2*z2 - 2*e0*i2*s0*s1*t2*z2 - 2*f1*k2*n0*o2*v0*z2 - 
   2*f0*k2*n1*o2*v0*z2 + 2*f1*j2*n0*p2*v0*z2 + 2*f0*j2*n1*p2*v0*z2 - 
   2*g1*j2*k2*r0*v0*z2 + 2*g1*i2*p2*r0*v0*z2 - 4*n1*n2*p2*r0*v0*z2 - 
   2*g0*j2*k2*r1*v0*z2 + 2*g0*i2*p2*r1*v0*z2 - 4*n0*n2*p2*r1*v0*z2 - 
   4*n0*n1*p2*r2*v0*z2 + 2*f1*j2*k2*s0*v0*z2 - 2*f1*i2*p2*s0*v0*z2 + 
   2*m2*n1*p2*s0*v0*z2 + 2*k2*n2*r1*s0*v0*z2 + 2*k2*n1*r2*s0*v0*z2 + 
   2*f0*j2*k2*s1*v0*z2 - 2*f0*i2*p2*s1*v0*z2 + 2*m2*n0*p2*s1*v0*z2 + 
   2*k2*n2*r0*s1*v0*z2 + 2*k2*n0*r2*s1*v0*z2 - 4*k2*m2*s0*s1*v0*z2 + 
   2*k2*n1*r0*s2*v0*z2 + 2*k2*n0*r1*s2*v0*z2 + 
   g1*pow(j2,2)*pow(v0,2)*z2 - g1*i2*o2*pow(v0,2)*z2 + 
   2*n1*n2*o2*pow(v0,2)*z2 - 2*j2*n2*s1*pow(v0,2)*z2 - 
   2*j2*n1*s2*pow(v0,2)*z2 + 2*i2*s1*s2*pow(v0,2)*z2 - 
   2*f0*k2*n0*o2*v1*z2 + 2*f0*j2*n0*p2*v1*z2 - 2*g0*j2*k2*r0*v1*z2 + 
   2*g0*i2*p2*r0*v1*z2 - 4*n0*n2*p2*r0*v1*z2 - 2*pow(n0,2)*p2*r2*v1*z2 + 
   2*f0*j2*k2*s0*v1*z2 - 2*f0*i2*p2*s0*v1*z2 + 2*m2*n0*p2*s0*v1*z2 + 
   2*k2*n2*r0*s0*v1*z2 + 2*k2*n0*r2*s0*v1*z2 - 2*k2*m2*pow(s0,2)*v1*z2 + 
   2*k2*n0*r0*s2*v1*z2 + 2*g0*pow(j2,2)*v0*v1*z2 - 2*g0*i2*o2*v0*v1*z2 + 
   4*n0*n2*o2*v0*v1*z2 - 4*j2*n2*s0*v0*v1*z2 - 4*j2*n0*s2*v0*v1*z2 + 
   4*i2*s0*s2*v0*v1*z2 - 4*n0*n1*p2*r0*v2*z2 - 2*pow(n0,2)*p2*r1*v2*z2 + 
   2*k2*n1*r0*s0*v2*z2 + 2*k2*n0*r1*s0*v2*z2 + 2*k2*n0*r0*s1*v2*z2 + 
   4*n0*n1*o2*v0*v2*z2 - 4*j2*n1*s0*v0*v2*z2 - 4*j2*n0*s1*v0*v2*z2 + 
   4*i2*s0*s1*v0*v2*z2 + 2*pow(n0,2)*o2*v1*v2*z2 - 4*j2*n0*s0*v1*v2*z2 + 
   2*i2*pow(s0,2)*v1*v2*z2 + 2*e1*k2*n0*o2*w0*z2 + 2*e0*k2*n1*o2*w0*z2 - 
   2*e1*j2*n0*p2*w0*z2 - 2*e0*j2*n1*p2*w0*z2 + 2*f1*j2*k2*r0*w0*z2 - 
   2*f1*i2*p2*r0*w0*z2 + 2*m2*n1*p2*r0*w0*z2 + 2*f0*j2*k2*r1*w0*z2 - 
   2*f0*i2*p2*r1*w0*z2 + 2*m2*n0*p2*r1*w0*z2 - 4*k2*n2*r0*r1*w0*z2 - 
   4*k2*n1*r0*r2*w0*z2 - 4*k2*n0*r1*r2*w0*z2 - 2*e1*j2*k2*s0*w0*z2 + 
   2*e1*i2*p2*s0*w0*z2 + 2*k2*m2*r1*s0*w0*z2 - 2*e0*j2*k2*s1*w0*z2 + 
   2*e0*i2*p2*s1*w0*z2 + 2*k2*m2*r0*s1*w0*z2 - 2*f1*pow(j2,2)*v0*w0*z2 + 
   2*f1*i2*o2*v0*w0*z2 - 2*m2*n1*o2*v0*w0*z2 + 2*j2*n2*r1*v0*w0*z2 + 
   2*j2*n1*r2*v0*w0*z2 + 2*j2*m2*s1*v0*w0*z2 - 2*i2*r2*s1*v0*w0*z2 - 
   2*i2*r1*s2*v0*w0*z2 - 2*f0*pow(j2,2)*v1*w0*z2 + 2*f0*i2*o2*v1*w0*z2 - 
   2*m2*n0*o2*v1*w0*z2 + 2*j2*n2*r0*v1*w0*z2 + 2*j2*n0*r2*v1*w0*z2 + 
   2*j2*m2*s0*v1*w0*z2 - 2*i2*r2*s0*v1*w0*z2 - 2*i2*r0*s2*v1*w0*z2 + 
   2*j2*n1*r0*v2*w0*z2 + 2*j2*n0*r1*v2*w0*z2 - 2*i2*r1*s0*v2*w0*z2 - 
   2*i2*r0*s1*v2*w0*z2 + e1*pow(j2,2)*pow(w0,2)*z2 - 
   e1*i2*o2*pow(w0,2)*z2 - 2*j2*m2*r1*pow(w0,2)*z2 + 
   2*i2*r1*r2*pow(w0,2)*z2 + 2*e0*k2*n0*o2*w1*z2 - 2*e0*j2*n0*p2*w1*z2 + 
   2*f0*j2*k2*r0*w1*z2 - 2*f0*i2*p2*r0*w1*z2 + 2*m2*n0*p2*r0*w1*z2 - 
   2*k2*n2*pow(r0,2)*w1*z2 - 4*k2*n0*r0*r2*w1*z2 - 2*e0*j2*k2*s0*w1*z2 + 
   2*e0*i2*p2*s0*w1*z2 + 2*k2*m2*r0*s0*w1*z2 - 2*f0*pow(j2,2)*v0*w1*z2 + 
   2*f0*i2*o2*v0*w1*z2 - 2*m2*n0*o2*v0*w1*z2 + 2*j2*n2*r0*v0*w1*z2 + 
   2*j2*n0*r2*v0*w1*z2 + 2*j2*m2*s0*v0*w1*z2 - 2*i2*r2*s0*v0*w1*z2 - 
   2*i2*r0*s2*v0*w1*z2 + 2*j2*n0*r0*v2*w1*z2 - 2*i2*r0*s0*v2*w1*z2 + 
   2*e0*pow(j2,2)*w0*w1*z2 - 2*e0*i2*o2*w0*w1*z2 - 4*j2*m2*r0*w0*w1*z2 + 
   4*i2*r0*r2*w0*w1*z2 - 2*k2*n1*pow(r0,2)*w2*z2 - 4*k2*n0*r0*r1*w2*z2 + 
   2*j2*n1*r0*v0*w2*z2 + 2*j2*n0*r1*v0*w2*z2 - 2*i2*r1*s0*v0*w2*z2 - 
   2*i2*r0*s1*v0*w2*z2 + 2*j2*n0*r0*v1*w2*z2 - 2*i2*r0*s0*v1*w2*z2 + 
   4*i2*r0*r1*w0*w2*z2 + 2*i2*pow(r0,2)*w1*w2*z2
;

	k[8]  = pow(d1,2)*e0*pow(k2,2)*o2 + 2*d0*d1*e1*pow(k2,2)*o2 - 
   2*c1*d1*f0*pow(k2,2)*o2 - 2*c1*d0*f1*pow(k2,2)*o2 - 
   2*c0*d1*f1*pow(k2,2)*o2 + pow(c1,2)*g0*pow(k2,2)*o2 + 
   2*c0*c1*g1*pow(k2,2)*o2 - 2*pow(d1,2)*e0*j2*k2*p2 - 
   4*d0*d1*e1*j2*k2*p2 + 4*c1*d1*f0*j2*k2*p2 + 4*c1*d0*f1*j2*k2*p2 + 
   4*c0*d1*f1*j2*k2*p2 - 2*pow(c1,2)*g0*j2*k2*p2 - 4*c0*c1*g1*j2*k2*p2 + 
   pow(d1,2)*e0*i2*pow(p2,2) + 2*d0*d1*e1*i2*pow(p2,2) - 
   2*c1*d1*f0*i2*pow(p2,2) - 2*c1*d0*f1*i2*pow(p2,2) - 
   2*c0*d1*f1*i2*pow(p2,2) + pow(c1,2)*g0*i2*pow(p2,2) + 
   2*c0*c1*g1*i2*pow(p2,2) - 2*d1*e1*l2*n0*pow(p2,2) + 
   2*c1*f1*l2*n0*pow(p2,2) + 2*c1*d1*m2*n0*pow(p2,2) - 
   2*d1*e0*l2*n1*pow(p2,2) - 2*d0*e1*l2*n1*pow(p2,2) + 
   2*c1*f0*l2*n1*pow(p2,2) + 2*c0*f1*l2*n1*pow(p2,2) + 
   2*c1*d0*m2*n1*pow(p2,2) + 2*c0*d1*m2*n1*pow(p2,2) - 
   4*c1*c2*n0*n1*pow(p2,2) - 2*c0*c2*pow(n1,2)*pow(p2,2) - 
   2*pow(c1,2)*n0*n2*pow(p2,2) - 4*c0*c1*n1*n2*pow(p2,2) + 
   2*d1*e1*k2*n0*p2*q2 - 2*c1*f1*k2*n0*p2*q2 + 2*d1*e0*k2*n1*p2*q2 + 
   2*d0*e1*k2*n1*p2*q2 - 2*c1*f0*k2*n1*p2*q2 - 2*c0*f1*k2*n1*p2*q2 - 
   2*d1*f1*k2*l2*p2*r0 + 2*c1*g1*k2*l2*p2*r0 + 2*pow(d1,2)*k2*m2*p2*r0 - 
   2*c2*d1*k2*n1*p2*r0 - 2*c1*d2*k2*n1*p2*r0 - 2*c1*d1*k2*n2*p2*r0 + 
   2*d1*f1*pow(k2,2)*q2*r0 - 2*c1*g1*pow(k2,2)*q2*r0 - 
   2*d1*f0*k2*l2*p2*r1 - 2*d0*f1*k2*l2*p2*r1 + 2*c1*g0*k2*l2*p2*r1 + 
   2*c0*g1*k2*l2*p2*r1 + 4*d0*d1*k2*m2*p2*r1 - 2*c2*d1*k2*n0*p2*r1 - 
   2*c1*d2*k2*n0*p2*r1 - 2*c2*d0*k2*n1*p2*r1 - 2*c0*d2*k2*n1*p2*r1 - 
   2*c1*d0*k2*n2*p2*r1 - 2*c0*d1*k2*n2*p2*r1 + 2*d1*f0*pow(k2,2)*q2*r1 + 
   2*d0*f1*pow(k2,2)*q2*r1 - 2*c1*g0*pow(k2,2)*q2*r1 - 
   2*c0*g1*pow(k2,2)*q2*r1 - 4*d1*d2*pow(k2,2)*r0*r1 - 
   2*d0*d2*pow(k2,2)*pow(r1,2) - 2*c1*d1*k2*n0*p2*r2 - 
   2*c1*d0*k2*n1*p2*r2 - 2*c0*d1*k2*n1*p2*r2 - 
   2*pow(d1,2)*pow(k2,2)*r0*r2 - 4*d0*d1*pow(k2,2)*r1*r2 + 
   2*d1*e1*k2*l2*p2*s0 - 2*c1*f1*k2*l2*p2*s0 - 2*c1*d1*k2*m2*p2*s0 + 
   4*c1*c2*k2*n1*p2*s0 + 2*pow(c1,2)*k2*n2*p2*s0 - 
   2*d1*e1*pow(k2,2)*q2*s0 + 2*c1*f1*pow(k2,2)*q2*s0 + 
   2*c2*d1*pow(k2,2)*r1*s0 + 2*c1*d2*pow(k2,2)*r1*s0 + 
   2*c1*d1*pow(k2,2)*r2*s0 + 2*d1*e0*k2*l2*p2*s1 + 2*d0*e1*k2*l2*p2*s1 - 
   2*c1*f0*k2*l2*p2*s1 - 2*c0*f1*k2*l2*p2*s1 - 2*c1*d0*k2*m2*p2*s1 - 
   2*c0*d1*k2*m2*p2*s1 + 4*c1*c2*k2*n0*p2*s1 + 4*c0*c2*k2*n1*p2*s1 + 
   4*c0*c1*k2*n2*p2*s1 - 2*d1*e0*pow(k2,2)*q2*s1 - 
   2*d0*e1*pow(k2,2)*q2*s1 + 2*c1*f0*pow(k2,2)*q2*s1 + 
   2*c0*f1*pow(k2,2)*q2*s1 + 2*c2*d1*pow(k2,2)*r0*s1 + 
   2*c1*d2*pow(k2,2)*r0*s1 + 2*c2*d0*pow(k2,2)*r1*s1 + 
   2*c0*d2*pow(k2,2)*r1*s1 + 2*c1*d0*pow(k2,2)*r2*s1 + 
   2*c0*d1*pow(k2,2)*r2*s1 - 4*c1*c2*pow(k2,2)*s0*s1 - 
   2*c0*c2*pow(k2,2)*pow(s1,2) + 2*pow(c1,2)*k2*n0*p2*s2 + 
   4*c0*c1*k2*n1*p2*s2 + 2*c1*d1*pow(k2,2)*r0*s2 + 
   2*c1*d0*pow(k2,2)*r1*s2 + 2*c0*d1*pow(k2,2)*r1*s2 - 
   2*pow(c1,2)*pow(k2,2)*s0*s2 - 4*c0*c1*pow(k2,2)*s1*s2 + 
   pow(d1,2)*e0*pow(j2,2)*t2 + 2*d0*d1*e1*pow(j2,2)*t2 - 
   2*c1*d1*f0*pow(j2,2)*t2 - 2*c1*d0*f1*pow(j2,2)*t2 - 
   2*c0*d1*f1*pow(j2,2)*t2 + pow(c1,2)*g0*pow(j2,2)*t2 + 
   2*c0*c1*g1*pow(j2,2)*t2 - pow(d1,2)*e0*i2*o2*t2 - 
   2*d0*d1*e1*i2*o2*t2 + 2*c1*d1*f0*i2*o2*t2 + 2*c1*d0*f1*i2*o2*t2 + 
   2*c0*d1*f1*i2*o2*t2 - pow(c1,2)*g0*i2*o2*t2 - 2*c0*c1*g1*i2*o2*t2 + 
   2*d1*e1*l2*n0*o2*t2 - 2*c1*f1*l2*n0*o2*t2 - 2*c1*d1*m2*n0*o2*t2 + 
   2*d1*e0*l2*n1*o2*t2 + 2*d0*e1*l2*n1*o2*t2 - 2*c1*f0*l2*n1*o2*t2 - 
   2*c0*f1*l2*n1*o2*t2 - 2*c1*d0*m2*n1*o2*t2 - 2*c0*d1*m2*n1*o2*t2 + 
   4*c1*c2*n0*n1*o2*t2 + 2*c0*c2*pow(n1,2)*o2*t2 + 
   2*pow(c1,2)*n0*n2*o2*t2 + 4*c0*c1*n1*n2*o2*t2 - 2*d1*e1*j2*n0*q2*t2 + 
   2*c1*f1*j2*n0*q2*t2 - 2*d1*e0*j2*n1*q2*t2 - 2*d0*e1*j2*n1*q2*t2 + 
   2*c1*f0*j2*n1*q2*t2 + 2*c0*f1*j2*n1*q2*t2 + 2*e1*n0*n1*pow(q2,2)*t2 + 
   e0*pow(n1,2)*pow(q2,2)*t2 + 2*d1*f1*j2*l2*r0*t2 - 
   2*c1*g1*j2*l2*r0*t2 - 2*pow(d1,2)*j2*m2*r0*t2 + 2*c2*d1*j2*n1*r0*t2 + 
   2*c1*d2*j2*n1*r0*t2 + 2*c1*d1*j2*n2*r0*t2 - 2*d1*f1*i2*q2*r0*t2 + 
   2*c1*g1*i2*q2*r0*t2 + 2*f1*l2*n1*q2*r0*t2 + 2*d1*m2*n1*q2*r0*t2 - 
   2*c2*pow(n1,2)*q2*r0*t2 - 4*c1*n1*n2*q2*r0*t2 + 2*d1*f0*j2*l2*r1*t2 + 
   2*d0*f1*j2*l2*r1*t2 - 2*c1*g0*j2*l2*r1*t2 - 2*c0*g1*j2*l2*r1*t2 - 
   4*d0*d1*j2*m2*r1*t2 + 2*c2*d1*j2*n0*r1*t2 + 2*c1*d2*j2*n0*r1*t2 + 
   2*c2*d0*j2*n1*r1*t2 + 2*c0*d2*j2*n1*r1*t2 + 2*c1*d0*j2*n2*r1*t2 + 
   2*c0*d1*j2*n2*r1*t2 - 2*d1*f0*i2*q2*r1*t2 - 2*d0*f1*i2*q2*r1*t2 + 
   2*c1*g0*i2*q2*r1*t2 + 2*c0*g1*i2*q2*r1*t2 + 2*f1*l2*n0*q2*r1*t2 + 
   2*d1*m2*n0*q2*r1*t2 + 2*f0*l2*n1*q2*r1*t2 + 2*d0*m2*n1*q2*r1*t2 - 
   4*c2*n0*n1*q2*r1*t2 - 4*c1*n0*n2*q2*r1*t2 - 4*c0*n1*n2*q2*r1*t2 + 
   4*d1*d2*i2*r0*r1*t2 + 2*g1*pow(l2,2)*r0*r1*t2 - 4*d2*l2*n1*r0*r1*t2 - 
   4*d1*l2*n2*r0*r1*t2 + 2*d0*d2*i2*pow(r1,2)*t2 + 
   g0*pow(l2,2)*pow(r1,2)*t2 - 2*d2*l2*n0*pow(r1,2)*t2 - 
   2*d0*l2*n2*pow(r1,2)*t2 + 2*c1*d1*j2*n0*r2*t2 + 2*c1*d0*j2*n1*r2*t2 + 
   2*c0*d1*j2*n1*r2*t2 - 4*c1*n0*n1*q2*r2*t2 - 2*c0*pow(n1,2)*q2*r2*t2 + 
   2*pow(d1,2)*i2*r0*r2*t2 - 4*d1*l2*n1*r0*r2*t2 + 4*d0*d1*i2*r1*r2*t2 - 
   4*d1*l2*n0*r1*r2*t2 - 4*d0*l2*n1*r1*r2*t2 - 2*d1*e1*j2*l2*s0*t2 + 
   2*c1*f1*j2*l2*s0*t2 + 2*c1*d1*j2*m2*s0*t2 - 4*c1*c2*j2*n1*s0*t2 - 
   2*pow(c1,2)*j2*n2*s0*t2 + 2*d1*e1*i2*q2*s0*t2 - 2*c1*f1*i2*q2*s0*t2 - 
   2*e1*l2*n1*q2*s0*t2 + 2*c1*m2*n1*q2*s0*t2 - 2*c2*d1*i2*r1*s0*t2 - 
   2*c1*d2*i2*r1*s0*t2 - 2*f1*pow(l2,2)*r1*s0*t2 + 2*d1*l2*m2*r1*s0*t2 + 
   2*c2*l2*n1*r1*s0*t2 + 2*c1*l2*n2*r1*s0*t2 - 2*c1*d1*i2*r2*s0*t2 + 
   2*c1*l2*n1*r2*s0*t2 - 2*d1*e0*j2*l2*s1*t2 - 2*d0*e1*j2*l2*s1*t2 + 
   2*c1*f0*j2*l2*s1*t2 + 2*c0*f1*j2*l2*s1*t2 + 2*c1*d0*j2*m2*s1*t2 + 
   2*c0*d1*j2*m2*s1*t2 - 4*c1*c2*j2*n0*s1*t2 - 4*c0*c2*j2*n1*s1*t2 - 
   4*c0*c1*j2*n2*s1*t2 + 2*d1*e0*i2*q2*s1*t2 + 2*d0*e1*i2*q2*s1*t2 - 
   2*c1*f0*i2*q2*s1*t2 - 2*c0*f1*i2*q2*s1*t2 - 2*e1*l2*n0*q2*s1*t2 + 
   2*c1*m2*n0*q2*s1*t2 - 2*e0*l2*n1*q2*s1*t2 + 2*c0*m2*n1*q2*s1*t2 - 
   2*c2*d1*i2*r0*s1*t2 - 2*c1*d2*i2*r0*s1*t2 - 2*f1*pow(l2,2)*r0*s1*t2 + 
   2*d1*l2*m2*r0*s1*t2 + 2*c2*l2*n1*r0*s1*t2 + 2*c1*l2*n2*r0*s1*t2 - 
   2*c2*d0*i2*r1*s1*t2 - 2*c0*d2*i2*r1*s1*t2 - 2*f0*pow(l2,2)*r1*s1*t2 + 
   2*d0*l2*m2*r1*s1*t2 + 2*c2*l2*n0*r1*s1*t2 + 2*c0*l2*n2*r1*s1*t2 - 
   2*c1*d0*i2*r2*s1*t2 - 2*c0*d1*i2*r2*s1*t2 + 2*c1*l2*n0*r2*s1*t2 + 
   2*c0*l2*n1*r2*s1*t2 + 4*c1*c2*i2*s0*s1*t2 + 2*e1*pow(l2,2)*s0*s1*t2 - 
   4*c1*l2*m2*s0*s1*t2 + 2*c0*c2*i2*pow(s1,2)*t2 + 
   e0*pow(l2,2)*pow(s1,2)*t2 - 2*c0*l2*m2*pow(s1,2)*t2 - 
   2*pow(c1,2)*j2*n0*s2*t2 - 4*c0*c1*j2*n1*s2*t2 - 2*c1*d1*i2*r0*s2*t2 + 
   2*c1*l2*n1*r0*s2*t2 - 2*c1*d0*i2*r1*s2*t2 - 2*c0*d1*i2*r1*s2*t2 + 
   2*c1*l2*n0*r1*s2*t2 + 2*c0*l2*n1*r1*s2*t2 + 2*pow(c1,2)*i2*s0*s2*t2 + 
   4*c0*c1*i2*s1*s2*t2 - 2*pow(f1,2)*k2*l2*o2*u0 + 2*e1*g1*k2*l2*o2*u0 + 
   2*d1*f1*k2*m2*o2*u0 - 2*c1*g1*k2*m2*o2*u0 - 2*d2*e1*k2*n1*o2*u0 - 
   2*d1*e2*k2*n1*o2*u0 + 2*c2*f1*k2*n1*o2*u0 + 2*c1*f2*k2*n1*o2*u0 - 
   2*d1*e1*k2*n2*o2*u0 + 2*c1*f1*k2*n2*o2*u0 + 2*pow(f1,2)*j2*l2*p2*u0 - 
   2*e1*g1*j2*l2*p2*u0 - 2*d1*f1*j2*m2*p2*u0 + 2*c1*g1*j2*m2*p2*u0 + 
   2*d2*e1*j2*n1*p2*u0 + 2*d1*e2*j2*n1*p2*u0 - 2*c2*f1*j2*n1*p2*u0 - 
   2*c1*f2*j2*n1*p2*u0 + 2*d1*e1*j2*n2*p2*u0 - 2*c1*f1*j2*n2*p2*u0 + 
   2*pow(f1,2)*j2*k2*q2*u0 - 2*e1*g1*j2*k2*q2*u0 - 
   2*pow(f1,2)*i2*p2*q2*u0 + 2*e1*g1*i2*p2*q2*u0 + 4*f1*m2*n1*p2*q2*u0 - 
   2*e2*pow(n1,2)*p2*q2*u0 - 4*e1*n1*n2*p2*q2*u0 - 2*d2*f1*j2*k2*r1*u0 - 
   2*d1*f2*j2*k2*r1*u0 + 2*c2*g1*j2*k2*r1*u0 + 2*c1*g2*j2*k2*r1*u0 + 
   2*d2*f1*i2*p2*r1*u0 + 2*d1*f2*i2*p2*r1*u0 - 2*c2*g1*i2*p2*r1*u0 - 
   2*c1*g2*i2*p2*r1*u0 + 2*g1*l2*m2*p2*r1*u0 - 2*f2*l2*n1*p2*r1*u0 - 
   2*d2*m2*n1*p2*r1*u0 - 2*f1*l2*n2*p2*r1*u0 - 2*d1*m2*n2*p2*r1*u0 + 
   4*c2*n1*n2*p2*r1*u0 + 2*c1*pow(n2,2)*p2*r1*u0 + 2*g1*k2*m2*q2*r1*u0 - 
   2*f2*k2*n1*q2*r1*u0 - 2*f1*k2*n2*q2*r1*u0 - 2*g2*k2*l2*pow(r1,2)*u0 + 
   2*d2*k2*n2*pow(r1,2)*u0 - 2*d1*f1*j2*k2*r2*u0 + 2*c1*g1*j2*k2*r2*u0 + 
   2*d1*f1*i2*p2*r2*u0 - 2*c1*g1*i2*p2*r2*u0 - 2*f1*l2*n1*p2*r2*u0 - 
   2*d1*m2*n1*p2*r2*u0 + 2*c2*pow(n1,2)*p2*r2*u0 + 4*c1*n1*n2*p2*r2*u0 - 
   2*f1*k2*n1*q2*r2*u0 - 4*g1*k2*l2*r1*r2*u0 + 4*d2*k2*n1*r1*r2*u0 + 
   4*d1*k2*n2*r1*r2*u0 + 2*d1*k2*n1*pow(r2,2)*u0 + 2*d2*e1*j2*k2*s1*u0 + 
   2*d1*e2*j2*k2*s1*u0 - 2*c2*f1*j2*k2*s1*u0 - 2*c1*f2*j2*k2*s1*u0 - 
   2*d2*e1*i2*p2*s1*u0 - 2*d1*e2*i2*p2*s1*u0 + 2*c2*f1*i2*p2*s1*u0 + 
   2*c1*f2*i2*p2*s1*u0 - 2*f1*l2*m2*p2*s1*u0 + 2*d1*pow(m2,2)*p2*s1*u0 + 
   2*e2*l2*n1*p2*s1*u0 - 2*c2*m2*n1*p2*s1*u0 + 2*e1*l2*n2*p2*s1*u0 - 
   2*c1*m2*n2*p2*s1*u0 - 2*f1*k2*m2*q2*s1*u0 + 2*e2*k2*n1*q2*s1*u0 + 
   2*e1*k2*n2*q2*s1*u0 + 4*f2*k2*l2*r1*s1*u0 - 2*d2*k2*m2*r1*s1*u0 - 
   2*c2*k2*n2*r1*s1*u0 + 4*f1*k2*l2*r2*s1*u0 - 2*d1*k2*m2*r2*s1*u0 - 
   2*c2*k2*n1*r2*s1*u0 - 2*c1*k2*n2*r2*s1*u0 - 2*e2*k2*l2*pow(s1,2)*u0 + 
   2*c2*k2*m2*pow(s1,2)*u0 + 2*d1*e1*j2*k2*s2*u0 - 2*c1*f1*j2*k2*s2*u0 - 
   2*d1*e1*i2*p2*s2*u0 + 2*c1*f1*i2*p2*s2*u0 + 2*e1*l2*n1*p2*s2*u0 - 
   2*c1*m2*n1*p2*s2*u0 + 2*e1*k2*n1*q2*s2*u0 + 4*f1*k2*l2*r1*s2*u0 - 
   2*d1*k2*m2*r1*s2*u0 - 2*c2*k2*n1*r1*s2*u0 - 2*c1*k2*n2*r1*s2*u0 - 
   2*c1*k2*n1*r2*s2*u0 - 4*e1*k2*l2*s1*s2*u0 + 4*c1*k2*m2*s1*s2*u0 - 
   4*f0*f1*k2*l2*o2*u1 + 2*e1*g0*k2*l2*o2*u1 + 2*e0*g1*k2*l2*o2*u1 + 
   2*d1*f0*k2*m2*o2*u1 + 2*d0*f1*k2*m2*o2*u1 - 2*c1*g0*k2*m2*o2*u1 - 
   2*c0*g1*k2*m2*o2*u1 - 2*d2*e1*k2*n0*o2*u1 - 2*d1*e2*k2*n0*o2*u1 + 
   2*c2*f1*k2*n0*o2*u1 + 2*c1*f2*k2*n0*o2*u1 - 2*d2*e0*k2*n1*o2*u1 - 
   2*d0*e2*k2*n1*o2*u1 + 2*c2*f0*k2*n1*o2*u1 + 2*c0*f2*k2*n1*o2*u1 - 
   2*d1*e0*k2*n2*o2*u1 - 2*d0*e1*k2*n2*o2*u1 + 2*c1*f0*k2*n2*o2*u1 + 
   2*c0*f1*k2*n2*o2*u1 + 4*f0*f1*j2*l2*p2*u1 - 2*e1*g0*j2*l2*p2*u1 - 
   2*e0*g1*j2*l2*p2*u1 - 2*d1*f0*j2*m2*p2*u1 - 2*d0*f1*j2*m2*p2*u1 + 
   2*c1*g0*j2*m2*p2*u1 + 2*c0*g1*j2*m2*p2*u1 + 2*d2*e1*j2*n0*p2*u1 + 
   2*d1*e2*j2*n0*p2*u1 - 2*c2*f1*j2*n0*p2*u1 - 2*c1*f2*j2*n0*p2*u1 + 
   2*d2*e0*j2*n1*p2*u1 + 2*d0*e2*j2*n1*p2*u1 - 2*c2*f0*j2*n1*p2*u1 - 
   2*c0*f2*j2*n1*p2*u1 + 2*d1*e0*j2*n2*p2*u1 + 2*d0*e1*j2*n2*p2*u1 - 
   2*c1*f0*j2*n2*p2*u1 - 2*c0*f1*j2*n2*p2*u1 + 4*f0*f1*j2*k2*q2*u1 - 
   2*e1*g0*j2*k2*q2*u1 - 2*e0*g1*j2*k2*q2*u1 - 4*f0*f1*i2*p2*q2*u1 + 
   2*e1*g0*i2*p2*q2*u1 + 2*e0*g1*i2*p2*q2*u1 + 4*f1*m2*n0*p2*q2*u1 + 
   4*f0*m2*n1*p2*q2*u1 - 4*e2*n0*n1*p2*q2*u1 - 4*e1*n0*n2*p2*q2*u1 - 
   4*e0*n1*n2*p2*q2*u1 - 2*d2*f1*j2*k2*r0*u1 - 2*d1*f2*j2*k2*r0*u1 + 
   2*c2*g1*j2*k2*r0*u1 + 2*c1*g2*j2*k2*r0*u1 + 2*d2*f1*i2*p2*r0*u1 + 
   2*d1*f2*i2*p2*r0*u1 - 2*c2*g1*i2*p2*r0*u1 - 2*c1*g2*i2*p2*r0*u1 + 
   2*g1*l2*m2*p2*r0*u1 - 2*f2*l2*n1*p2*r0*u1 - 2*d2*m2*n1*p2*r0*u1 - 
   2*f1*l2*n2*p2*r0*u1 - 2*d1*m2*n2*p2*r0*u1 + 4*c2*n1*n2*p2*r0*u1 + 
   2*c1*pow(n2,2)*p2*r0*u1 + 2*g1*k2*m2*q2*r0*u1 - 2*f2*k2*n1*q2*r0*u1 - 
   2*f1*k2*n2*q2*r0*u1 - 2*d2*f0*j2*k2*r1*u1 - 2*d0*f2*j2*k2*r1*u1 + 
   2*c2*g0*j2*k2*r1*u1 + 2*c0*g2*j2*k2*r1*u1 + 2*d2*f0*i2*p2*r1*u1 + 
   2*d0*f2*i2*p2*r1*u1 - 2*c2*g0*i2*p2*r1*u1 - 2*c0*g2*i2*p2*r1*u1 + 
   2*g0*l2*m2*p2*r1*u1 - 2*f2*l2*n0*p2*r1*u1 - 2*d2*m2*n0*p2*r1*u1 - 
   2*f0*l2*n2*p2*r1*u1 - 2*d0*m2*n2*p2*r1*u1 + 4*c2*n0*n2*p2*r1*u1 + 
   2*c0*pow(n2,2)*p2*r1*u1 + 2*g0*k2*m2*q2*r1*u1 - 2*f2*k2*n0*q2*r1*u1 - 
   2*f0*k2*n2*q2*r1*u1 - 4*g2*k2*l2*r0*r1*u1 + 4*d2*k2*n2*r0*r1*u1 - 
   2*d1*f0*j2*k2*r2*u1 - 2*d0*f1*j2*k2*r2*u1 + 2*c1*g0*j2*k2*r2*u1 + 
   2*c0*g1*j2*k2*r2*u1 + 2*d1*f0*i2*p2*r2*u1 + 2*d0*f1*i2*p2*r2*u1 - 
   2*c1*g0*i2*p2*r2*u1 - 2*c0*g1*i2*p2*r2*u1 - 2*f1*l2*n0*p2*r2*u1 - 
   2*d1*m2*n0*p2*r2*u1 - 2*f0*l2*n1*p2*r2*u1 - 2*d0*m2*n1*p2*r2*u1 + 
   4*c2*n0*n1*p2*r2*u1 + 4*c1*n0*n2*p2*r2*u1 + 4*c0*n1*n2*p2*r2*u1 - 
   2*f1*k2*n0*q2*r2*u1 - 2*f0*k2*n1*q2*r2*u1 - 4*g1*k2*l2*r0*r2*u1 + 
   4*d2*k2*n1*r0*r2*u1 + 4*d1*k2*n2*r0*r2*u1 - 4*g0*k2*l2*r1*r2*u1 + 
   4*d2*k2*n0*r1*r2*u1 + 4*d0*k2*n2*r1*r2*u1 + 2*d1*k2*n0*pow(r2,2)*u1 + 
   2*d0*k2*n1*pow(r2,2)*u1 + 2*d2*e1*j2*k2*s0*u1 + 2*d1*e2*j2*k2*s0*u1 - 
   2*c2*f1*j2*k2*s0*u1 - 2*c1*f2*j2*k2*s0*u1 - 2*d2*e1*i2*p2*s0*u1 - 
   2*d1*e2*i2*p2*s0*u1 + 2*c2*f1*i2*p2*s0*u1 + 2*c1*f2*i2*p2*s0*u1 - 
   2*f1*l2*m2*p2*s0*u1 + 2*d1*pow(m2,2)*p2*s0*u1 + 2*e2*l2*n1*p2*s0*u1 - 
   2*c2*m2*n1*p2*s0*u1 + 2*e1*l2*n2*p2*s0*u1 - 2*c1*m2*n2*p2*s0*u1 - 
   2*f1*k2*m2*q2*s0*u1 + 2*e2*k2*n1*q2*s0*u1 + 2*e1*k2*n2*q2*s0*u1 + 
   4*f2*k2*l2*r1*s0*u1 - 2*d2*k2*m2*r1*s0*u1 - 2*c2*k2*n2*r1*s0*u1 + 
   4*f1*k2*l2*r2*s0*u1 - 2*d1*k2*m2*r2*s0*u1 - 2*c2*k2*n1*r2*s0*u1 - 
   2*c1*k2*n2*r2*s0*u1 + 2*d2*e0*j2*k2*s1*u1 + 2*d0*e2*j2*k2*s1*u1 - 
   2*c2*f0*j2*k2*s1*u1 - 2*c0*f2*j2*k2*s1*u1 - 2*d2*e0*i2*p2*s1*u1 - 
   2*d0*e2*i2*p2*s1*u1 + 2*c2*f0*i2*p2*s1*u1 + 2*c0*f2*i2*p2*s1*u1 - 
   2*f0*l2*m2*p2*s1*u1 + 2*d0*pow(m2,2)*p2*s1*u1 + 2*e2*l2*n0*p2*s1*u1 - 
   2*c2*m2*n0*p2*s1*u1 + 2*e0*l2*n2*p2*s1*u1 - 2*c0*m2*n2*p2*s1*u1 - 
   2*f0*k2*m2*q2*s1*u1 + 2*e2*k2*n0*q2*s1*u1 + 2*e0*k2*n2*q2*s1*u1 + 
   4*f2*k2*l2*r0*s1*u1 - 2*d2*k2*m2*r0*s1*u1 - 2*c2*k2*n2*r0*s1*u1 + 
   4*f0*k2*l2*r2*s1*u1 - 2*d0*k2*m2*r2*s1*u1 - 2*c2*k2*n0*r2*s1*u1 - 
   2*c0*k2*n2*r2*s1*u1 - 4*e2*k2*l2*s0*s1*u1 + 4*c2*k2*m2*s0*s1*u1 + 
   2*d1*e0*j2*k2*s2*u1 + 2*d0*e1*j2*k2*s2*u1 - 2*c1*f0*j2*k2*s2*u1 - 
   2*c0*f1*j2*k2*s2*u1 - 2*d1*e0*i2*p2*s2*u1 - 2*d0*e1*i2*p2*s2*u1 + 
   2*c1*f0*i2*p2*s2*u1 + 2*c0*f1*i2*p2*s2*u1 + 2*e1*l2*n0*p2*s2*u1 - 
   2*c1*m2*n0*p2*s2*u1 + 2*e0*l2*n1*p2*s2*u1 - 2*c0*m2*n1*p2*s2*u1 + 
   2*e1*k2*n0*q2*s2*u1 + 2*e0*k2*n1*q2*s2*u1 + 4*f1*k2*l2*r0*s2*u1 - 
   2*d1*k2*m2*r0*s2*u1 - 2*c2*k2*n1*r0*s2*u1 - 2*c1*k2*n2*r0*s2*u1 + 
   4*f0*k2*l2*r1*s2*u1 - 2*d0*k2*m2*r1*s2*u1 - 2*c2*k2*n0*r1*s2*u1 - 
   2*c0*k2*n2*r1*s2*u1 - 2*c1*k2*n0*r2*s2*u1 - 2*c0*k2*n1*r2*s2*u1 - 
   4*e1*k2*l2*s0*s2*u1 + 4*c1*k2*m2*s0*s2*u1 - 4*e0*k2*l2*s1*s2*u1 + 
   4*c0*k2*m2*s1*s2*u1 - 4*f1*f2*pow(j2,2)*u0*u1 + 
   2*e2*g1*pow(j2,2)*u0*u1 + 2*e1*g2*pow(j2,2)*u0*u1 + 
   4*f1*f2*i2*o2*u0*u1 - 2*e2*g1*i2*o2*u0*u1 - 2*e1*g2*i2*o2*u0*u1 + 
   2*g1*pow(m2,2)*o2*u0*u1 - 4*f2*m2*n1*o2*u0*u1 - 4*f1*m2*n2*o2*u0*u1 + 
   4*e2*n1*n2*o2*u0*u1 + 2*e1*pow(n2,2)*o2*u0*u1 - 4*g2*j2*m2*r1*u0*u1 + 
   4*f2*j2*n2*r1*u0*u1 - 4*g1*j2*m2*r2*u0*u1 + 4*f2*j2*n1*r2*u0*u1 + 
   4*f1*j2*n2*r2*u0*u1 + 4*g2*i2*r1*r2*u0*u1 - 4*pow(n2,2)*r1*r2*u0*u1 + 
   2*g1*i2*pow(r2,2)*u0*u1 - 4*n1*n2*pow(r2,2)*u0*u1 + 
   4*f2*j2*m2*s1*u0*u1 - 4*e2*j2*n2*s1*u0*u1 - 4*f2*i2*r2*s1*u0*u1 + 
   4*m2*n2*r2*s1*u0*u1 + 4*f1*j2*m2*s2*u0*u1 - 4*e2*j2*n1*s2*u0*u1 - 
   4*e1*j2*n2*s2*u0*u1 - 4*f2*i2*r1*s2*u0*u1 + 4*m2*n2*r1*s2*u0*u1 - 
   4*f1*i2*r2*s2*u0*u1 + 4*m2*n1*r2*s2*u0*u1 + 4*e2*i2*s1*s2*u0*u1 - 
   4*pow(m2,2)*s1*s2*u0*u1 + 2*e1*i2*pow(s2,2)*u0*u1 - 
   2*f0*f2*pow(j2,2)*pow(u1,2) + e2*g0*pow(j2,2)*pow(u1,2) + 
   e0*g2*pow(j2,2)*pow(u1,2) + 2*f0*f2*i2*o2*pow(u1,2) - 
   e2*g0*i2*o2*pow(u1,2) - e0*g2*i2*o2*pow(u1,2) + 
   g0*pow(m2,2)*o2*pow(u1,2) - 2*f2*m2*n0*o2*pow(u1,2) - 
   2*f0*m2*n2*o2*pow(u1,2) + 2*e2*n0*n2*o2*pow(u1,2) + 
   e0*pow(n2,2)*o2*pow(u1,2) - 2*g2*j2*m2*r0*pow(u1,2) + 
   2*f2*j2*n2*r0*pow(u1,2) - 2*g0*j2*m2*r2*pow(u1,2) + 
   2*f2*j2*n0*r2*pow(u1,2) + 2*f0*j2*n2*r2*pow(u1,2) + 
   2*g2*i2*r0*r2*pow(u1,2) - 2*pow(n2,2)*r0*r2*pow(u1,2) + 
   g0*i2*pow(r2,2)*pow(u1,2) - 2*n0*n2*pow(r2,2)*pow(u1,2) + 
   2*f2*j2*m2*s0*pow(u1,2) - 2*e2*j2*n2*s0*pow(u1,2) - 
   2*f2*i2*r2*s0*pow(u1,2) + 2*m2*n2*r2*s0*pow(u1,2) + 
   2*f0*j2*m2*s2*pow(u1,2) - 2*e2*j2*n0*s2*pow(u1,2) - 
   2*e0*j2*n2*s2*pow(u1,2) - 2*f2*i2*r0*s2*pow(u1,2) + 
   2*m2*n2*r0*s2*pow(u1,2) - 2*f0*i2*r2*s2*pow(u1,2) + 
   2*m2*n0*r2*s2*pow(u1,2) + 2*e2*i2*s0*s2*pow(u1,2) - 
   2*pow(m2,2)*s0*s2*pow(u1,2) + e0*i2*pow(s2,2)*pow(u1,2) - 
   2*d1*e1*k2*n0*o2*u2 + 2*c1*f1*k2*n0*o2*u2 - 2*d1*e0*k2*n1*o2*u2 - 
   2*d0*e1*k2*n1*o2*u2 + 2*c1*f0*k2*n1*o2*u2 + 2*c0*f1*k2*n1*o2*u2 + 
   2*d1*e1*j2*n0*p2*u2 - 2*c1*f1*j2*n0*p2*u2 + 2*d1*e0*j2*n1*p2*u2 + 
   2*d0*e1*j2*n1*p2*u2 - 2*c1*f0*j2*n1*p2*u2 - 2*c0*f1*j2*n1*p2*u2 - 
   4*e1*n0*n1*p2*q2*u2 - 2*e0*pow(n1,2)*p2*q2*u2 - 2*d1*f1*j2*k2*r0*u2 + 
   2*c1*g1*j2*k2*r0*u2 + 2*d1*f1*i2*p2*r0*u2 - 2*c1*g1*i2*p2*r0*u2 - 
   2*f1*l2*n1*p2*r0*u2 - 2*d1*m2*n1*p2*r0*u2 + 2*c2*pow(n1,2)*p2*r0*u2 + 
   4*c1*n1*n2*p2*r0*u2 - 2*f1*k2*n1*q2*r0*u2 - 2*d1*f0*j2*k2*r1*u2 - 
   2*d0*f1*j2*k2*r1*u2 + 2*c1*g0*j2*k2*r1*u2 + 2*c0*g1*j2*k2*r1*u2 + 
   2*d1*f0*i2*p2*r1*u2 + 2*d0*f1*i2*p2*r1*u2 - 2*c1*g0*i2*p2*r1*u2 - 
   2*c0*g1*i2*p2*r1*u2 - 2*f1*l2*n0*p2*r1*u2 - 2*d1*m2*n0*p2*r1*u2 - 
   2*f0*l2*n1*p2*r1*u2 - 2*d0*m2*n1*p2*r1*u2 + 4*c2*n0*n1*p2*r1*u2 + 
   4*c1*n0*n2*p2*r1*u2 + 4*c0*n1*n2*p2*r1*u2 - 2*f1*k2*n0*q2*r1*u2 - 
   2*f0*k2*n1*q2*r1*u2 - 4*g1*k2*l2*r0*r1*u2 + 4*d2*k2*n1*r0*r1*u2 + 
   4*d1*k2*n2*r0*r1*u2 - 2*g0*k2*l2*pow(r1,2)*u2 + 
   2*d2*k2*n0*pow(r1,2)*u2 + 2*d0*k2*n2*pow(r1,2)*u2 + 
   4*c1*n0*n1*p2*r2*u2 + 2*c0*pow(n1,2)*p2*r2*u2 + 4*d1*k2*n1*r0*r2*u2 + 
   4*d1*k2*n0*r1*r2*u2 + 4*d0*k2*n1*r1*r2*u2 + 2*d1*e1*j2*k2*s0*u2 - 
   2*c1*f1*j2*k2*s0*u2 - 2*d1*e1*i2*p2*s0*u2 + 2*c1*f1*i2*p2*s0*u2 + 
   2*e1*l2*n1*p2*s0*u2 - 2*c1*m2*n1*p2*s0*u2 + 2*e1*k2*n1*q2*s0*u2 + 
   4*f1*k2*l2*r1*s0*u2 - 2*d1*k2*m2*r1*s0*u2 - 2*c2*k2*n1*r1*s0*u2 - 
   2*c1*k2*n2*r1*s0*u2 - 2*c1*k2*n1*r2*s0*u2 + 2*d1*e0*j2*k2*s1*u2 + 
   2*d0*e1*j2*k2*s1*u2 - 2*c1*f0*j2*k2*s1*u2 - 2*c0*f1*j2*k2*s1*u2 - 
   2*d1*e0*i2*p2*s1*u2 - 2*d0*e1*i2*p2*s1*u2 + 2*c1*f0*i2*p2*s1*u2 + 
   2*c0*f1*i2*p2*s1*u2 + 2*e1*l2*n0*p2*s1*u2 - 2*c1*m2*n0*p2*s1*u2 + 
   2*e0*l2*n1*p2*s1*u2 - 2*c0*m2*n1*p2*s1*u2 + 2*e1*k2*n0*q2*s1*u2 + 
   2*e0*k2*n1*q2*s1*u2 + 4*f1*k2*l2*r0*s1*u2 - 2*d1*k2*m2*r0*s1*u2 - 
   2*c2*k2*n1*r0*s1*u2 - 2*c1*k2*n2*r0*s1*u2 + 4*f0*k2*l2*r1*s1*u2 - 
   2*d0*k2*m2*r1*s1*u2 - 2*c2*k2*n0*r1*s1*u2 - 2*c0*k2*n2*r1*s1*u2 - 
   2*c1*k2*n0*r2*s1*u2 - 2*c0*k2*n1*r2*s1*u2 - 4*e1*k2*l2*s0*s1*u2 + 
   4*c1*k2*m2*s0*s1*u2 - 2*e0*k2*l2*pow(s1,2)*u2 + 
   2*c0*k2*m2*pow(s1,2)*u2 - 2*c1*k2*n1*r0*s2*u2 - 2*c1*k2*n0*r1*s2*u2 - 
   2*c0*k2*n1*r1*s2*u2 - 2*pow(f1,2)*pow(j2,2)*u0*u2 + 
   2*e1*g1*pow(j2,2)*u0*u2 + 2*pow(f1,2)*i2*o2*u0*u2 - 
   2*e1*g1*i2*o2*u0*u2 - 4*f1*m2*n1*o2*u0*u2 + 2*e2*pow(n1,2)*o2*u0*u2 + 
   4*e1*n1*n2*o2*u0*u2 - 4*g1*j2*m2*r1*u0*u2 + 4*f2*j2*n1*r1*u0*u2 + 
   4*f1*j2*n2*r1*u0*u2 + 2*g2*i2*pow(r1,2)*u0*u2 - 
   2*pow(n2,2)*pow(r1,2)*u0*u2 + 4*f1*j2*n1*r2*u0*u2 + 
   4*g1*i2*r1*r2*u0*u2 - 8*n1*n2*r1*r2*u0*u2 - 
   2*pow(n1,2)*pow(r2,2)*u0*u2 + 4*f1*j2*m2*s1*u0*u2 - 
   4*e2*j2*n1*s1*u0*u2 - 4*e1*j2*n2*s1*u0*u2 - 4*f2*i2*r1*s1*u0*u2 + 
   4*m2*n2*r1*s1*u0*u2 - 4*f1*i2*r2*s1*u0*u2 + 4*m2*n1*r2*s1*u0*u2 + 
   2*e2*i2*pow(s1,2)*u0*u2 - 2*pow(m2,2)*pow(s1,2)*u0*u2 - 
   4*e1*j2*n1*s2*u0*u2 - 4*f1*i2*r1*s2*u0*u2 + 4*m2*n1*r1*s2*u0*u2 + 
   4*e1*i2*s1*s2*u0*u2 - 4*f0*f1*pow(j2,2)*u1*u2 + 
   2*e1*g0*pow(j2,2)*u1*u2 + 2*e0*g1*pow(j2,2)*u1*u2 + 
   4*f0*f1*i2*o2*u1*u2 - 2*e1*g0*i2*o2*u1*u2 - 2*e0*g1*i2*o2*u1*u2 - 
   4*f1*m2*n0*o2*u1*u2 - 4*f0*m2*n1*o2*u1*u2 + 4*e2*n0*n1*o2*u1*u2 + 
   4*e1*n0*n2*o2*u1*u2 + 4*e0*n1*n2*o2*u1*u2 - 4*g1*j2*m2*r0*u1*u2 + 
   4*f2*j2*n1*r0*u1*u2 + 4*f1*j2*n2*r0*u1*u2 - 4*g0*j2*m2*r1*u1*u2 + 
   4*f2*j2*n0*r1*u1*u2 + 4*f0*j2*n2*r1*u1*u2 + 4*g2*i2*r0*r1*u1*u2 - 
   4*pow(n2,2)*r0*r1*u1*u2 + 4*f1*j2*n0*r2*u1*u2 + 4*f0*j2*n1*r2*u1*u2 + 
   4*g1*i2*r0*r2*u1*u2 - 8*n1*n2*r0*r2*u1*u2 + 4*g0*i2*r1*r2*u1*u2 - 
   8*n0*n2*r1*r2*u1*u2 - 4*n0*n1*pow(r2,2)*u1*u2 + 4*f1*j2*m2*s0*u1*u2 - 
   4*e2*j2*n1*s0*u1*u2 - 4*e1*j2*n2*s0*u1*u2 - 4*f2*i2*r1*s0*u1*u2 + 
   4*m2*n2*r1*s0*u1*u2 - 4*f1*i2*r2*s0*u1*u2 + 4*m2*n1*r2*s0*u1*u2 + 
   4*f0*j2*m2*s1*u1*u2 - 4*e2*j2*n0*s1*u1*u2 - 4*e0*j2*n2*s1*u1*u2 - 
   4*f2*i2*r0*s1*u1*u2 + 4*m2*n2*r0*s1*u1*u2 - 4*f0*i2*r2*s1*u1*u2 + 
   4*m2*n0*r2*s1*u1*u2 + 4*e2*i2*s0*s1*u1*u2 - 4*pow(m2,2)*s0*s1*u1*u2 - 
   4*e1*j2*n0*s2*u1*u2 - 4*e0*j2*n1*s2*u1*u2 - 4*f1*i2*r0*s2*u1*u2 + 
   4*m2*n1*r0*s2*u1*u2 - 4*f0*i2*r1*s2*u1*u2 + 4*m2*n0*r1*s2*u1*u2 + 
   4*e1*i2*s0*s2*u1*u2 + 4*e0*i2*s1*s2*u1*u2 + 2*e1*n0*n1*o2*pow(u2,2) + 
   e0*pow(n1,2)*o2*pow(u2,2) + 2*f1*j2*n1*r0*pow(u2,2) + 
   2*f1*j2*n0*r1*pow(u2,2) + 2*f0*j2*n1*r1*pow(u2,2) + 
   2*g1*i2*r0*r1*pow(u2,2) - 4*n1*n2*r0*r1*pow(u2,2) + 
   g0*i2*pow(r1,2)*pow(u2,2) - 2*n0*n2*pow(r1,2)*pow(u2,2) - 
   2*pow(n1,2)*r0*r2*pow(u2,2) - 4*n0*n1*r1*r2*pow(u2,2) - 
   2*e1*j2*n1*s0*pow(u2,2) - 2*f1*i2*r1*s0*pow(u2,2) + 
   2*m2*n1*r1*s0*pow(u2,2) - 2*e1*j2*n0*s1*pow(u2,2) - 
   2*e0*j2*n1*s1*pow(u2,2) - 2*f1*i2*r0*s1*pow(u2,2) + 
   2*m2*n1*r0*s1*pow(u2,2) - 2*f0*i2*r1*s1*pow(u2,2) + 
   2*m2*n0*r1*s1*pow(u2,2) + 2*e1*i2*s0*s1*pow(u2,2) + 
   e0*i2*pow(s1,2)*pow(u2,2) + 2*d1*f1*k2*l2*o2*v0 - 
   2*c1*g1*k2*l2*o2*v0 - 2*pow(d1,2)*k2*m2*o2*v0 + 2*c2*d1*k2*n1*o2*v0 + 
   2*c1*d2*k2*n1*o2*v0 + 2*c1*d1*k2*n2*o2*v0 - 2*d1*f1*j2*l2*p2*v0 + 
   2*c1*g1*j2*l2*p2*v0 + 2*pow(d1,2)*j2*m2*p2*v0 - 2*c2*d1*j2*n1*p2*v0 - 
   2*c1*d2*j2*n1*p2*v0 - 2*c1*d1*j2*n2*p2*v0 - 2*d1*f1*j2*k2*q2*v0 + 
   2*c1*g1*j2*k2*q2*v0 + 2*d1*f1*i2*p2*q2*v0 - 2*c1*g1*i2*p2*q2*v0 - 
   2*f1*l2*n1*p2*q2*v0 - 2*d1*m2*n1*p2*q2*v0 + 2*c2*pow(n1,2)*p2*q2*v0 + 
   4*c1*n1*n2*p2*q2*v0 + 2*f1*k2*n1*pow(q2,2)*v0 + 4*d1*d2*j2*k2*r1*v0 - 
   4*d1*d2*i2*p2*r1*v0 - 2*g1*pow(l2,2)*p2*r1*v0 + 4*d2*l2*n1*p2*r1*v0 + 
   4*d1*l2*n2*p2*r1*v0 + 2*g1*k2*l2*q2*r1*v0 - 2*d2*k2*n1*q2*r1*v0 - 
   2*d1*k2*n2*q2*r1*v0 + 2*pow(d1,2)*j2*k2*r2*v0 - 
   2*pow(d1,2)*i2*p2*r2*v0 + 4*d1*l2*n1*p2*r2*v0 - 2*d1*k2*n1*q2*r2*v0 - 
   2*c2*d1*j2*k2*s1*v0 - 2*c1*d2*j2*k2*s1*v0 + 2*c2*d1*i2*p2*s1*v0 + 
   2*c1*d2*i2*p2*s1*v0 + 2*f1*pow(l2,2)*p2*s1*v0 - 2*d1*l2*m2*p2*s1*v0 - 
   2*c2*l2*n1*p2*s1*v0 - 2*c1*l2*n2*p2*s1*v0 - 2*f1*k2*l2*q2*s1*v0 + 
   4*d1*k2*m2*q2*s1*v0 - 2*c2*k2*n1*q2*s1*v0 - 2*c1*k2*n2*q2*s1*v0 - 
   2*d2*k2*l2*r1*s1*v0 - 2*d1*k2*l2*r2*s1*v0 + 2*c2*k2*l2*pow(s1,2)*v0 - 
   2*c1*d1*j2*k2*s2*v0 + 2*c1*d1*i2*p2*s2*v0 - 2*c1*l2*n1*p2*s2*v0 - 
   2*c1*k2*n1*q2*s2*v0 - 2*d1*k2*l2*r1*s2*v0 + 4*c1*k2*l2*s1*s2*v0 + 
   2*d2*f1*pow(j2,2)*u1*v0 + 2*d1*f2*pow(j2,2)*u1*v0 - 
   2*c2*g1*pow(j2,2)*u1*v0 - 2*c1*g2*pow(j2,2)*u1*v0 - 
   2*d2*f1*i2*o2*u1*v0 - 2*d1*f2*i2*o2*u1*v0 + 2*c2*g1*i2*o2*u1*v0 + 
   2*c1*g2*i2*o2*u1*v0 - 2*g1*l2*m2*o2*u1*v0 + 2*f2*l2*n1*o2*u1*v0 + 
   2*d2*m2*n1*o2*u1*v0 + 2*f1*l2*n2*o2*u1*v0 + 2*d1*m2*n2*o2*u1*v0 - 
   4*c2*n1*n2*o2*u1*v0 - 2*c1*pow(n2,2)*o2*u1*v0 + 2*g1*j2*m2*q2*u1*v0 - 
   2*f2*j2*n1*q2*u1*v0 - 2*f1*j2*n2*q2*u1*v0 + 2*g2*j2*l2*r1*u1*v0 - 
   2*d2*j2*n2*r1*u1*v0 - 2*g2*i2*q2*r1*u1*v0 + 2*pow(n2,2)*q2*r1*u1*v0 + 
   2*g1*j2*l2*r2*u1*v0 - 2*d2*j2*n1*r2*u1*v0 - 2*d1*j2*n2*r2*u1*v0 - 
   2*g1*i2*q2*r2*u1*v0 + 4*n1*n2*q2*r2*u1*v0 - 2*f2*j2*l2*s1*u1*v0 - 
   2*d2*j2*m2*s1*u1*v0 + 4*c2*j2*n2*s1*u1*v0 + 2*f2*i2*q2*s1*u1*v0 - 
   2*m2*n2*q2*s1*u1*v0 + 2*d2*i2*r2*s1*u1*v0 - 2*l2*n2*r2*s1*u1*v0 - 
   2*f1*j2*l2*s2*u1*v0 - 2*d1*j2*m2*s2*u1*v0 + 4*c2*j2*n1*s2*u1*v0 + 
   4*c1*j2*n2*s2*u1*v0 + 2*f1*i2*q2*s2*u1*v0 - 2*m2*n1*q2*s2*u1*v0 + 
   2*d2*i2*r1*s2*u1*v0 - 2*l2*n2*r1*s2*u1*v0 + 2*d1*i2*r2*s2*u1*v0 - 
   2*l2*n1*r2*s2*u1*v0 - 4*c2*i2*s1*s2*u1*v0 + 4*l2*m2*s1*s2*u1*v0 - 
   2*c1*i2*pow(s2,2)*u1*v0 + 2*d1*f1*pow(j2,2)*u2*v0 - 
   2*c1*g1*pow(j2,2)*u2*v0 - 2*d1*f1*i2*o2*u2*v0 + 2*c1*g1*i2*o2*u2*v0 + 
   2*f1*l2*n1*o2*u2*v0 + 2*d1*m2*n1*o2*u2*v0 - 2*c2*pow(n1,2)*o2*u2*v0 - 
   4*c1*n1*n2*o2*u2*v0 - 2*f1*j2*n1*q2*u2*v0 + 2*g1*j2*l2*r1*u2*v0 - 
   2*d2*j2*n1*r1*u2*v0 - 2*d1*j2*n2*r1*u2*v0 - 2*g1*i2*q2*r1*u2*v0 + 
   4*n1*n2*q2*r1*u2*v0 - 2*d1*j2*n1*r2*u2*v0 + 2*pow(n1,2)*q2*r2*u2*v0 - 
   2*f1*j2*l2*s1*u2*v0 - 2*d1*j2*m2*s1*u2*v0 + 4*c2*j2*n1*s1*u2*v0 + 
   4*c1*j2*n2*s1*u2*v0 + 2*f1*i2*q2*s1*u2*v0 - 2*m2*n1*q2*s1*u2*v0 + 
   2*d2*i2*r1*s1*u2*v0 - 2*l2*n2*r1*s1*u2*v0 + 2*d1*i2*r2*s1*u2*v0 - 
   2*l2*n1*r2*s1*u2*v0 - 2*c2*i2*pow(s1,2)*u2*v0 + 
   2*l2*m2*pow(s1,2)*u2*v0 + 4*c1*j2*n1*s2*u2*v0 + 2*d1*i2*r1*s2*u2*v0 - 
   2*l2*n1*r1*s2*u2*v0 - 4*c1*i2*s1*s2*u2*v0 + 2*d1*f0*k2*l2*o2*v1 + 
   2*d0*f1*k2*l2*o2*v1 - 2*c1*g0*k2*l2*o2*v1 - 2*c0*g1*k2*l2*o2*v1 - 
   4*d0*d1*k2*m2*o2*v1 + 2*c2*d1*k2*n0*o2*v1 + 2*c1*d2*k2*n0*o2*v1 + 
   2*c2*d0*k2*n1*o2*v1 + 2*c0*d2*k2*n1*o2*v1 + 2*c1*d0*k2*n2*o2*v1 + 
   2*c0*d1*k2*n2*o2*v1 - 2*d1*f0*j2*l2*p2*v1 - 2*d0*f1*j2*l2*p2*v1 + 
   2*c1*g0*j2*l2*p2*v1 + 2*c0*g1*j2*l2*p2*v1 + 4*d0*d1*j2*m2*p2*v1 - 
   2*c2*d1*j2*n0*p2*v1 - 2*c1*d2*j2*n0*p2*v1 - 2*c2*d0*j2*n1*p2*v1 - 
   2*c0*d2*j2*n1*p2*v1 - 2*c1*d0*j2*n2*p2*v1 - 2*c0*d1*j2*n2*p2*v1 - 
   2*d1*f0*j2*k2*q2*v1 - 2*d0*f1*j2*k2*q2*v1 + 2*c1*g0*j2*k2*q2*v1 + 
   2*c0*g1*j2*k2*q2*v1 + 2*d1*f0*i2*p2*q2*v1 + 2*d0*f1*i2*p2*q2*v1 - 
   2*c1*g0*i2*p2*q2*v1 - 2*c0*g1*i2*p2*q2*v1 - 2*f1*l2*n0*p2*q2*v1 - 
   2*d1*m2*n0*p2*q2*v1 - 2*f0*l2*n1*p2*q2*v1 - 2*d0*m2*n1*p2*q2*v1 + 
   4*c2*n0*n1*p2*q2*v1 + 4*c1*n0*n2*p2*q2*v1 + 4*c0*n1*n2*p2*q2*v1 + 
   2*f1*k2*n0*pow(q2,2)*v1 + 2*f0*k2*n1*pow(q2,2)*v1 + 
   4*d1*d2*j2*k2*r0*v1 - 4*d1*d2*i2*p2*r0*v1 - 2*g1*pow(l2,2)*p2*r0*v1 + 
   4*d2*l2*n1*p2*r0*v1 + 4*d1*l2*n2*p2*r0*v1 + 2*g1*k2*l2*q2*r0*v1 - 
   2*d2*k2*n1*q2*r0*v1 - 2*d1*k2*n2*q2*r0*v1 + 4*d0*d2*j2*k2*r1*v1 - 
   4*d0*d2*i2*p2*r1*v1 - 2*g0*pow(l2,2)*p2*r1*v1 + 4*d2*l2*n0*p2*r1*v1 + 
   4*d0*l2*n2*p2*r1*v1 + 2*g0*k2*l2*q2*r1*v1 - 2*d2*k2*n0*q2*r1*v1 - 
   2*d0*k2*n2*q2*r1*v1 + 4*d0*d1*j2*k2*r2*v1 - 4*d0*d1*i2*p2*r2*v1 + 
   4*d1*l2*n0*p2*r2*v1 + 4*d0*l2*n1*p2*r2*v1 - 2*d1*k2*n0*q2*r2*v1 - 
   2*d0*k2*n1*q2*r2*v1 - 2*c2*d1*j2*k2*s0*v1 - 2*c1*d2*j2*k2*s0*v1 + 
   2*c2*d1*i2*p2*s0*v1 + 2*c1*d2*i2*p2*s0*v1 + 2*f1*pow(l2,2)*p2*s0*v1 - 
   2*d1*l2*m2*p2*s0*v1 - 2*c2*l2*n1*p2*s0*v1 - 2*c1*l2*n2*p2*s0*v1 - 
   2*f1*k2*l2*q2*s0*v1 + 4*d1*k2*m2*q2*s0*v1 - 2*c2*k2*n1*q2*s0*v1 - 
   2*c1*k2*n2*q2*s0*v1 - 2*d2*k2*l2*r1*s0*v1 - 2*d1*k2*l2*r2*s0*v1 - 
   2*c2*d0*j2*k2*s1*v1 - 2*c0*d2*j2*k2*s1*v1 + 2*c2*d0*i2*p2*s1*v1 + 
   2*c0*d2*i2*p2*s1*v1 + 2*f0*pow(l2,2)*p2*s1*v1 - 2*d0*l2*m2*p2*s1*v1 - 
   2*c2*l2*n0*p2*s1*v1 - 2*c0*l2*n2*p2*s1*v1 - 2*f0*k2*l2*q2*s1*v1 + 
   4*d0*k2*m2*q2*s1*v1 - 2*c2*k2*n0*q2*s1*v1 - 2*c0*k2*n2*q2*s1*v1 - 
   2*d2*k2*l2*r0*s1*v1 - 2*d0*k2*l2*r2*s1*v1 + 4*c2*k2*l2*s0*s1*v1 - 
   2*c1*d0*j2*k2*s2*v1 - 2*c0*d1*j2*k2*s2*v1 + 2*c1*d0*i2*p2*s2*v1 + 
   2*c0*d1*i2*p2*s2*v1 - 2*c1*l2*n0*p2*s2*v1 - 2*c0*l2*n1*p2*s2*v1 - 
   2*c1*k2*n0*q2*s2*v1 - 2*c0*k2*n1*q2*s2*v1 - 2*d1*k2*l2*r0*s2*v1 - 
   2*d0*k2*l2*r1*s2*v1 + 4*c1*k2*l2*s0*s2*v1 + 4*c0*k2*l2*s1*s2*v1 + 
   2*d2*f1*pow(j2,2)*u0*v1 + 2*d1*f2*pow(j2,2)*u0*v1 - 
   2*c2*g1*pow(j2,2)*u0*v1 - 2*c1*g2*pow(j2,2)*u0*v1 - 
   2*d2*f1*i2*o2*u0*v1 - 2*d1*f2*i2*o2*u0*v1 + 2*c2*g1*i2*o2*u0*v1 + 
   2*c1*g2*i2*o2*u0*v1 - 2*g1*l2*m2*o2*u0*v1 + 2*f2*l2*n1*o2*u0*v1 + 
   2*d2*m2*n1*o2*u0*v1 + 2*f1*l2*n2*o2*u0*v1 + 2*d1*m2*n2*o2*u0*v1 - 
   4*c2*n1*n2*o2*u0*v1 - 2*c1*pow(n2,2)*o2*u0*v1 + 2*g1*j2*m2*q2*u0*v1 - 
   2*f2*j2*n1*q2*u0*v1 - 2*f1*j2*n2*q2*u0*v1 + 2*g2*j2*l2*r1*u0*v1 - 
   2*d2*j2*n2*r1*u0*v1 - 2*g2*i2*q2*r1*u0*v1 + 2*pow(n2,2)*q2*r1*u0*v1 + 
   2*g1*j2*l2*r2*u0*v1 - 2*d2*j2*n1*r2*u0*v1 - 2*d1*j2*n2*r2*u0*v1 - 
   2*g1*i2*q2*r2*u0*v1 + 4*n1*n2*q2*r2*u0*v1 - 2*f2*j2*l2*s1*u0*v1 - 
   2*d2*j2*m2*s1*u0*v1 + 4*c2*j2*n2*s1*u0*v1 + 2*f2*i2*q2*s1*u0*v1 - 
   2*m2*n2*q2*s1*u0*v1 + 2*d2*i2*r2*s1*u0*v1 - 2*l2*n2*r2*s1*u0*v1 - 
   2*f1*j2*l2*s2*u0*v1 - 2*d1*j2*m2*s2*u0*v1 + 4*c2*j2*n1*s2*u0*v1 + 
   4*c1*j2*n2*s2*u0*v1 + 2*f1*i2*q2*s2*u0*v1 - 2*m2*n1*q2*s2*u0*v1 + 
   2*d2*i2*r1*s2*u0*v1 - 2*l2*n2*r1*s2*u0*v1 + 2*d1*i2*r2*s2*u0*v1 - 
   2*l2*n1*r2*s2*u0*v1 - 4*c2*i2*s1*s2*u0*v1 + 4*l2*m2*s1*s2*u0*v1 - 
   2*c1*i2*pow(s2,2)*u0*v1 + 2*d2*f0*pow(j2,2)*u1*v1 + 
   2*d0*f2*pow(j2,2)*u1*v1 - 2*c2*g0*pow(j2,2)*u1*v1 - 
   2*c0*g2*pow(j2,2)*u1*v1 - 2*d2*f0*i2*o2*u1*v1 - 2*d0*f2*i2*o2*u1*v1 + 
   2*c2*g0*i2*o2*u1*v1 + 2*c0*g2*i2*o2*u1*v1 - 2*g0*l2*m2*o2*u1*v1 + 
   2*f2*l2*n0*o2*u1*v1 + 2*d2*m2*n0*o2*u1*v1 + 2*f0*l2*n2*o2*u1*v1 + 
   2*d0*m2*n2*o2*u1*v1 - 4*c2*n0*n2*o2*u1*v1 - 2*c0*pow(n2,2)*o2*u1*v1 + 
   2*g0*j2*m2*q2*u1*v1 - 2*f2*j2*n0*q2*u1*v1 - 2*f0*j2*n2*q2*u1*v1 + 
   2*g2*j2*l2*r0*u1*v1 - 2*d2*j2*n2*r0*u1*v1 - 2*g2*i2*q2*r0*u1*v1 + 
   2*pow(n2,2)*q2*r0*u1*v1 + 2*g0*j2*l2*r2*u1*v1 - 2*d2*j2*n0*r2*u1*v1 - 
   2*d0*j2*n2*r2*u1*v1 - 2*g0*i2*q2*r2*u1*v1 + 4*n0*n2*q2*r2*u1*v1 - 
   2*f2*j2*l2*s0*u1*v1 - 2*d2*j2*m2*s0*u1*v1 + 4*c2*j2*n2*s0*u1*v1 + 
   2*f2*i2*q2*s0*u1*v1 - 2*m2*n2*q2*s0*u1*v1 + 2*d2*i2*r2*s0*u1*v1 - 
   2*l2*n2*r2*s0*u1*v1 - 2*f0*j2*l2*s2*u1*v1 - 2*d0*j2*m2*s2*u1*v1 + 
   4*c2*j2*n0*s2*u1*v1 + 4*c0*j2*n2*s2*u1*v1 + 2*f0*i2*q2*s2*u1*v1 - 
   2*m2*n0*q2*s2*u1*v1 + 2*d2*i2*r0*s2*u1*v1 - 2*l2*n2*r0*s2*u1*v1 + 
   2*d0*i2*r2*s2*u1*v1 - 2*l2*n0*r2*s2*u1*v1 - 4*c2*i2*s0*s2*u1*v1 + 
   4*l2*m2*s0*s2*u1*v1 - 2*c0*i2*pow(s2,2)*u1*v1 + 
   2*d1*f0*pow(j2,2)*u2*v1 + 2*d0*f1*pow(j2,2)*u2*v1 - 
   2*c1*g0*pow(j2,2)*u2*v1 - 2*c0*g1*pow(j2,2)*u2*v1 - 
   2*d1*f0*i2*o2*u2*v1 - 2*d0*f1*i2*o2*u2*v1 + 2*c1*g0*i2*o2*u2*v1 + 
   2*c0*g1*i2*o2*u2*v1 + 2*f1*l2*n0*o2*u2*v1 + 2*d1*m2*n0*o2*u2*v1 + 
   2*f0*l2*n1*o2*u2*v1 + 2*d0*m2*n1*o2*u2*v1 - 4*c2*n0*n1*o2*u2*v1 - 
   4*c1*n0*n2*o2*u2*v1 - 4*c0*n1*n2*o2*u2*v1 - 2*f1*j2*n0*q2*u2*v1 - 
   2*f0*j2*n1*q2*u2*v1 + 2*g1*j2*l2*r0*u2*v1 - 2*d2*j2*n1*r0*u2*v1 - 
   2*d1*j2*n2*r0*u2*v1 - 2*g1*i2*q2*r0*u2*v1 + 4*n1*n2*q2*r0*u2*v1 + 
   2*g0*j2*l2*r1*u2*v1 - 2*d2*j2*n0*r1*u2*v1 - 2*d0*j2*n2*r1*u2*v1 - 
   2*g0*i2*q2*r1*u2*v1 + 4*n0*n2*q2*r1*u2*v1 - 2*d1*j2*n0*r2*u2*v1 - 
   2*d0*j2*n1*r2*u2*v1 + 4*n0*n1*q2*r2*u2*v1 - 2*f1*j2*l2*s0*u2*v1 - 
   2*d1*j2*m2*s0*u2*v1 + 4*c2*j2*n1*s0*u2*v1 + 4*c1*j2*n2*s0*u2*v1 + 
   2*f1*i2*q2*s0*u2*v1 - 2*m2*n1*q2*s0*u2*v1 + 2*d2*i2*r1*s0*u2*v1 - 
   2*l2*n2*r1*s0*u2*v1 + 2*d1*i2*r2*s0*u2*v1 - 2*l2*n1*r2*s0*u2*v1 - 
   2*f0*j2*l2*s1*u2*v1 - 2*d0*j2*m2*s1*u2*v1 + 4*c2*j2*n0*s1*u2*v1 + 
   4*c0*j2*n2*s1*u2*v1 + 2*f0*i2*q2*s1*u2*v1 - 2*m2*n0*q2*s1*u2*v1 + 
   2*d2*i2*r0*s1*u2*v1 - 2*l2*n2*r0*s1*u2*v1 + 2*d0*i2*r2*s1*u2*v1 - 
   2*l2*n0*r2*s1*u2*v1 - 4*c2*i2*s0*s1*u2*v1 + 4*l2*m2*s0*s1*u2*v1 + 
   4*c1*j2*n0*s2*u2*v1 + 4*c0*j2*n1*s2*u2*v1 + 2*d1*i2*r0*s2*u2*v1 - 
   2*l2*n1*r0*s2*u2*v1 + 2*d0*i2*r1*s2*u2*v1 - 2*l2*n0*r1*s2*u2*v1 - 
   4*c1*i2*s0*s2*u2*v1 - 4*c0*i2*s1*s2*u2*v1 - 4*d1*d2*pow(j2,2)*v0*v1 + 
   4*d1*d2*i2*o2*v0*v1 + 2*g1*pow(l2,2)*o2*v0*v1 - 4*d2*l2*n1*o2*v0*v1 - 
   4*d1*l2*n2*o2*v0*v1 - 4*g1*j2*l2*q2*v0*v1 + 4*d2*j2*n1*q2*v0*v1 + 
   4*d1*j2*n2*q2*v0*v1 + 2*g1*i2*pow(q2,2)*v0*v1 - 
   4*n1*n2*pow(q2,2)*v0*v1 + 4*d2*j2*l2*s1*v0*v1 - 4*d2*i2*q2*s1*v0*v1 + 
   4*l2*n2*q2*s1*v0*v1 + 4*d1*j2*l2*s2*v0*v1 - 4*d1*i2*q2*s2*v0*v1 + 
   4*l2*n1*q2*s2*v0*v1 - 4*pow(l2,2)*s1*s2*v0*v1 - 
   2*d0*d2*pow(j2,2)*pow(v1,2) + 2*d0*d2*i2*o2*pow(v1,2) + 
   g0*pow(l2,2)*o2*pow(v1,2) - 2*d2*l2*n0*o2*pow(v1,2) - 
   2*d0*l2*n2*o2*pow(v1,2) - 2*g0*j2*l2*q2*pow(v1,2) + 
   2*d2*j2*n0*q2*pow(v1,2) + 2*d0*j2*n2*q2*pow(v1,2) + 
   g0*i2*pow(q2,2)*pow(v1,2) - 2*n0*n2*pow(q2,2)*pow(v1,2) + 
   2*d2*j2*l2*s0*pow(v1,2) - 2*d2*i2*q2*s0*pow(v1,2) + 
   2*l2*n2*q2*s0*pow(v1,2) + 2*d0*j2*l2*s2*pow(v1,2) - 
   2*d0*i2*q2*s2*pow(v1,2) + 2*l2*n0*q2*s2*pow(v1,2) - 
   2*pow(l2,2)*s0*s2*pow(v1,2) + 2*c1*d1*k2*n0*o2*v2 + 
   2*c1*d0*k2*n1*o2*v2 + 2*c0*d1*k2*n1*o2*v2 - 2*c1*d1*j2*n0*p2*v2 - 
   2*c1*d0*j2*n1*p2*v2 - 2*c0*d1*j2*n1*p2*v2 + 4*c1*n0*n1*p2*q2*v2 + 
   2*c0*pow(n1,2)*p2*q2*v2 + 2*pow(d1,2)*j2*k2*r0*v2 - 
   2*pow(d1,2)*i2*p2*r0*v2 + 4*d1*l2*n1*p2*r0*v2 - 2*d1*k2*n1*q2*r0*v2 + 
   4*d0*d1*j2*k2*r1*v2 - 4*d0*d1*i2*p2*r1*v2 + 4*d1*l2*n0*p2*r1*v2 + 
   4*d0*l2*n1*p2*r1*v2 - 2*d1*k2*n0*q2*r1*v2 - 2*d0*k2*n1*q2*r1*v2 - 
   2*c1*d1*j2*k2*s0*v2 + 2*c1*d1*i2*p2*s0*v2 - 2*c1*l2*n1*p2*s0*v2 - 
   2*c1*k2*n1*q2*s0*v2 - 2*d1*k2*l2*r1*s0*v2 - 2*c1*d0*j2*k2*s1*v2 - 
   2*c0*d1*j2*k2*s1*v2 + 2*c1*d0*i2*p2*s1*v2 + 2*c0*d1*i2*p2*s1*v2 - 
   2*c1*l2*n0*p2*s1*v2 - 2*c0*l2*n1*p2*s1*v2 - 2*c1*k2*n0*q2*s1*v2 - 
   2*c0*k2*n1*q2*s1*v2 - 2*d1*k2*l2*r0*s1*v2 - 2*d0*k2*l2*r1*s1*v2 + 
   4*c1*k2*l2*s0*s1*v2 + 2*c0*k2*l2*pow(s1,2)*v2 + 
   2*d1*f1*pow(j2,2)*u0*v2 - 2*c1*g1*pow(j2,2)*u0*v2 - 
   2*d1*f1*i2*o2*u0*v2 + 2*c1*g1*i2*o2*u0*v2 + 2*f1*l2*n1*o2*u0*v2 + 
   2*d1*m2*n1*o2*u0*v2 - 2*c2*pow(n1,2)*o2*u0*v2 - 4*c1*n1*n2*o2*u0*v2 - 
   2*f1*j2*n1*q2*u0*v2 + 2*g1*j2*l2*r1*u0*v2 - 2*d2*j2*n1*r1*u0*v2 - 
   2*d1*j2*n2*r1*u0*v2 - 2*g1*i2*q2*r1*u0*v2 + 4*n1*n2*q2*r1*u0*v2 - 
   2*d1*j2*n1*r2*u0*v2 + 2*pow(n1,2)*q2*r2*u0*v2 - 2*f1*j2*l2*s1*u0*v2 - 
   2*d1*j2*m2*s1*u0*v2 + 4*c2*j2*n1*s1*u0*v2 + 4*c1*j2*n2*s1*u0*v2 + 
   2*f1*i2*q2*s1*u0*v2 - 2*m2*n1*q2*s1*u0*v2 + 2*d2*i2*r1*s1*u0*v2 - 
   2*l2*n2*r1*s1*u0*v2 + 2*d1*i2*r2*s1*u0*v2 - 2*l2*n1*r2*s1*u0*v2 - 
   2*c2*i2*pow(s1,2)*u0*v2 + 2*l2*m2*pow(s1,2)*u0*v2 + 
   4*c1*j2*n1*s2*u0*v2 + 2*d1*i2*r1*s2*u0*v2 - 2*l2*n1*r1*s2*u0*v2 - 
   4*c1*i2*s1*s2*u0*v2 + 2*d1*f0*pow(j2,2)*u1*v2 + 
   2*d0*f1*pow(j2,2)*u1*v2 - 2*c1*g0*pow(j2,2)*u1*v2 - 
   2*c0*g1*pow(j2,2)*u1*v2 - 2*d1*f0*i2*o2*u1*v2 - 2*d0*f1*i2*o2*u1*v2 + 
   2*c1*g0*i2*o2*u1*v2 + 2*c0*g1*i2*o2*u1*v2 + 2*f1*l2*n0*o2*u1*v2 + 
   2*d1*m2*n0*o2*u1*v2 + 2*f0*l2*n1*o2*u1*v2 + 2*d0*m2*n1*o2*u1*v2 - 
   4*c2*n0*n1*o2*u1*v2 - 4*c1*n0*n2*o2*u1*v2 - 4*c0*n1*n2*o2*u1*v2 - 
   2*f1*j2*n0*q2*u1*v2 - 2*f0*j2*n1*q2*u1*v2 + 2*g1*j2*l2*r0*u1*v2 - 
   2*d2*j2*n1*r0*u1*v2 - 2*d1*j2*n2*r0*u1*v2 - 2*g1*i2*q2*r0*u1*v2 + 
   4*n1*n2*q2*r0*u1*v2 + 2*g0*j2*l2*r1*u1*v2 - 2*d2*j2*n0*r1*u1*v2 - 
   2*d0*j2*n2*r1*u1*v2 - 2*g0*i2*q2*r1*u1*v2 + 4*n0*n2*q2*r1*u1*v2 - 
   2*d1*j2*n0*r2*u1*v2 - 2*d0*j2*n1*r2*u1*v2 + 4*n0*n1*q2*r2*u1*v2 - 
   2*f1*j2*l2*s0*u1*v2 - 2*d1*j2*m2*s0*u1*v2 + 4*c2*j2*n1*s0*u1*v2 + 
   4*c1*j2*n2*s0*u1*v2 + 2*f1*i2*q2*s0*u1*v2 - 2*m2*n1*q2*s0*u1*v2 + 
   2*d2*i2*r1*s0*u1*v2 - 2*l2*n2*r1*s0*u1*v2 + 2*d1*i2*r2*s0*u1*v2 - 
   2*l2*n1*r2*s0*u1*v2 - 2*f0*j2*l2*s1*u1*v2 - 2*d0*j2*m2*s1*u1*v2 + 
   4*c2*j2*n0*s1*u1*v2 + 4*c0*j2*n2*s1*u1*v2 + 2*f0*i2*q2*s1*u1*v2 - 
   2*m2*n0*q2*s1*u1*v2 + 2*d2*i2*r0*s1*u1*v2 - 2*l2*n2*r0*s1*u1*v2 + 
   2*d0*i2*r2*s1*u1*v2 - 2*l2*n0*r2*s1*u1*v2 - 4*c2*i2*s0*s1*u1*v2 + 
   4*l2*m2*s0*s1*u1*v2 + 4*c1*j2*n0*s2*u1*v2 + 4*c0*j2*n1*s2*u1*v2 + 
   2*d1*i2*r0*s2*u1*v2 - 2*l2*n1*r0*s2*u1*v2 + 2*d0*i2*r1*s2*u1*v2 - 
   2*l2*n0*r1*s2*u1*v2 - 4*c1*i2*s0*s2*u1*v2 - 4*c0*i2*s1*s2*u1*v2 - 
   4*c1*n0*n1*o2*u2*v2 - 2*c0*pow(n1,2)*o2*u2*v2 - 2*d1*j2*n1*r0*u2*v2 + 
   2*pow(n1,2)*q2*r0*u2*v2 - 2*d1*j2*n0*r1*u2*v2 - 2*d0*j2*n1*r1*u2*v2 + 
   4*n0*n1*q2*r1*u2*v2 + 4*c1*j2*n1*s0*u2*v2 + 2*d1*i2*r1*s0*u2*v2 - 
   2*l2*n1*r1*s0*u2*v2 + 4*c1*j2*n0*s1*u2*v2 + 4*c0*j2*n1*s1*u2*v2 + 
   2*d1*i2*r0*s1*u2*v2 - 2*l2*n1*r0*s1*u2*v2 + 2*d0*i2*r1*s1*u2*v2 - 
   2*l2*n0*r1*s1*u2*v2 - 4*c1*i2*s0*s1*u2*v2 - 2*c0*i2*pow(s1,2)*u2*v2 - 
   2*pow(d1,2)*pow(j2,2)*v0*v2 + 2*pow(d1,2)*i2*o2*v0*v2 - 
   4*d1*l2*n1*o2*v0*v2 + 4*d1*j2*n1*q2*v0*v2 - 
   2*pow(n1,2)*pow(q2,2)*v0*v2 + 4*d1*j2*l2*s1*v0*v2 - 
   4*d1*i2*q2*s1*v0*v2 + 4*l2*n1*q2*s1*v0*v2 - 
   2*pow(l2,2)*pow(s1,2)*v0*v2 - 4*d0*d1*pow(j2,2)*v1*v2 + 
   4*d0*d1*i2*o2*v1*v2 - 4*d1*l2*n0*o2*v1*v2 - 4*d0*l2*n1*o2*v1*v2 + 
   4*d1*j2*n0*q2*v1*v2 + 4*d0*j2*n1*q2*v1*v2 - 4*n0*n1*pow(q2,2)*v1*v2 + 
   4*d1*j2*l2*s0*v1*v2 - 4*d1*i2*q2*s0*v1*v2 + 4*l2*n1*q2*s0*v1*v2 + 
   4*d0*j2*l2*s1*v1*v2 - 4*d0*i2*q2*s1*v1*v2 + 4*l2*n0*q2*s1*v1*v2 - 
   4*pow(l2,2)*s0*s1*v1*v2 - 2*d1*e1*k2*l2*o2*w0 + 2*c1*f1*k2*l2*o2*w0 + 
   2*c1*d1*k2*m2*o2*w0 - 4*c1*c2*k2*n1*o2*w0 - 2*pow(c1,2)*k2*n2*o2*w0 + 
   2*d1*e1*j2*l2*p2*w0 - 2*c1*f1*j2*l2*p2*w0 - 2*c1*d1*j2*m2*p2*w0 + 
   4*c1*c2*j2*n1*p2*w0 + 2*pow(c1,2)*j2*n2*p2*w0 + 2*d1*e1*j2*k2*q2*w0 - 
   2*c1*f1*j2*k2*q2*w0 - 2*d1*e1*i2*p2*q2*w0 + 2*c1*f1*i2*p2*q2*w0 + 
   2*e1*l2*n1*p2*q2*w0 - 2*c1*m2*n1*p2*q2*w0 - 2*e1*k2*n1*pow(q2,2)*w0 - 
   2*c2*d1*j2*k2*r1*w0 - 2*c1*d2*j2*k2*r1*w0 + 2*c2*d1*i2*p2*r1*w0 + 
   2*c1*d2*i2*p2*r1*w0 + 2*f1*pow(l2,2)*p2*r1*w0 - 2*d1*l2*m2*p2*r1*w0 - 
   2*c2*l2*n1*p2*r1*w0 - 2*c1*l2*n2*p2*r1*w0 - 2*f1*k2*l2*q2*r1*w0 - 
   2*d1*k2*m2*q2*r1*w0 + 4*c2*k2*n1*q2*r1*w0 + 4*c1*k2*n2*q2*r1*w0 + 
   2*d2*k2*l2*pow(r1,2)*w0 - 2*c1*d1*j2*k2*r2*w0 + 2*c1*d1*i2*p2*r2*w0 - 
   2*c1*l2*n1*p2*r2*w0 + 4*c1*k2*n1*q2*r2*w0 + 4*d1*k2*l2*r1*r2*w0 + 
   4*c1*c2*j2*k2*s1*w0 - 4*c1*c2*i2*p2*s1*w0 - 2*e1*pow(l2,2)*p2*s1*w0 + 
   4*c1*l2*m2*p2*s1*w0 + 2*e1*k2*l2*q2*s1*w0 - 2*c1*k2*m2*q2*s1*w0 - 
   2*c2*k2*l2*r1*s1*w0 - 2*c1*k2*l2*r2*s1*w0 + 2*pow(c1,2)*j2*k2*s2*w0 - 
   2*pow(c1,2)*i2*p2*s2*w0 - 2*c1*k2*l2*r1*s2*w0 - 
   2*d2*e1*pow(j2,2)*u1*w0 - 2*d1*e2*pow(j2,2)*u1*w0 + 
   2*c2*f1*pow(j2,2)*u1*w0 + 2*c1*f2*pow(j2,2)*u1*w0 + 
   2*d2*e1*i2*o2*u1*w0 + 2*d1*e2*i2*o2*u1*w0 - 2*c2*f1*i2*o2*u1*w0 - 
   2*c1*f2*i2*o2*u1*w0 + 2*f1*l2*m2*o2*u1*w0 - 2*d1*pow(m2,2)*o2*u1*w0 - 
   2*e2*l2*n1*o2*u1*w0 + 2*c2*m2*n1*o2*u1*w0 - 2*e1*l2*n2*o2*u1*w0 + 
   2*c1*m2*n2*o2*u1*w0 - 2*f1*j2*m2*q2*u1*w0 + 2*e2*j2*n1*q2*u1*w0 + 
   2*e1*j2*n2*q2*u1*w0 - 2*f2*j2*l2*r1*u1*w0 + 4*d2*j2*m2*r1*u1*w0 - 
   2*c2*j2*n2*r1*u1*w0 + 2*f2*i2*q2*r1*u1*w0 - 2*m2*n2*q2*r1*u1*w0 - 
   2*f1*j2*l2*r2*u1*w0 + 4*d1*j2*m2*r2*u1*w0 - 2*c2*j2*n1*r2*u1*w0 - 
   2*c1*j2*n2*r2*u1*w0 + 2*f1*i2*q2*r2*u1*w0 - 2*m2*n1*q2*r2*u1*w0 - 
   4*d2*i2*r1*r2*u1*w0 + 4*l2*n2*r1*r2*u1*w0 - 2*d1*i2*pow(r2,2)*u1*w0 + 
   2*l2*n1*pow(r2,2)*u1*w0 + 2*e2*j2*l2*s1*u1*w0 - 2*c2*j2*m2*s1*u1*w0 - 
   2*e2*i2*q2*s1*u1*w0 + 2*pow(m2,2)*q2*s1*u1*w0 + 2*c2*i2*r2*s1*u1*w0 - 
   2*l2*m2*r2*s1*u1*w0 + 2*e1*j2*l2*s2*u1*w0 - 2*c1*j2*m2*s2*u1*w0 - 
   2*e1*i2*q2*s2*u1*w0 + 2*c2*i2*r1*s2*u1*w0 - 2*l2*m2*r1*s2*u1*w0 + 
   2*c1*i2*r2*s2*u1*w0 - 2*d1*e1*pow(j2,2)*u2*w0 + 
   2*c1*f1*pow(j2,2)*u2*w0 + 2*d1*e1*i2*o2*u2*w0 - 2*c1*f1*i2*o2*u2*w0 - 
   2*e1*l2*n1*o2*u2*w0 + 2*c1*m2*n1*o2*u2*w0 + 2*e1*j2*n1*q2*u2*w0 - 
   2*f1*j2*l2*r1*u2*w0 + 4*d1*j2*m2*r1*u2*w0 - 2*c2*j2*n1*r1*u2*w0 - 
   2*c1*j2*n2*r1*u2*w0 + 2*f1*i2*q2*r1*u2*w0 - 2*m2*n1*q2*r1*u2*w0 - 
   2*d2*i2*pow(r1,2)*u2*w0 + 2*l2*n2*pow(r1,2)*u2*w0 - 
   2*c1*j2*n1*r2*u2*w0 - 4*d1*i2*r1*r2*u2*w0 + 4*l2*n1*r1*r2*u2*w0 + 
   2*e1*j2*l2*s1*u2*w0 - 2*c1*j2*m2*s1*u2*w0 - 2*e1*i2*q2*s1*u2*w0 + 
   2*c2*i2*r1*s1*u2*w0 - 2*l2*m2*r1*s1*u2*w0 + 2*c1*i2*r2*s1*u2*w0 + 
   2*c1*i2*r1*s2*u2*w0 + 2*c2*d1*pow(j2,2)*v1*w0 + 
   2*c1*d2*pow(j2,2)*v1*w0 - 2*c2*d1*i2*o2*v1*w0 - 2*c1*d2*i2*o2*v1*w0 - 
   2*f1*pow(l2,2)*o2*v1*w0 + 2*d1*l2*m2*o2*v1*w0 + 2*c2*l2*n1*o2*v1*w0 + 
   2*c1*l2*n2*o2*v1*w0 + 4*f1*j2*l2*q2*v1*w0 - 2*d1*j2*m2*q2*v1*w0 - 
   2*c2*j2*n1*q2*v1*w0 - 2*c1*j2*n2*q2*v1*w0 - 2*f1*i2*pow(q2,2)*v1*w0 + 
   2*m2*n1*pow(q2,2)*v1*w0 - 2*d2*j2*l2*r1*v1*w0 + 2*d2*i2*q2*r1*v1*w0 - 
   2*l2*n2*q2*r1*v1*w0 - 2*d1*j2*l2*r2*v1*w0 + 2*d1*i2*q2*r2*v1*w0 - 
   2*l2*n1*q2*r2*v1*w0 - 2*c2*j2*l2*s1*v1*w0 + 2*c2*i2*q2*s1*v1*w0 - 
   2*l2*m2*q2*s1*v1*w0 + 2*pow(l2,2)*r2*s1*v1*w0 - 2*c1*j2*l2*s2*v1*w0 + 
   2*c1*i2*q2*s2*v1*w0 + 2*pow(l2,2)*r1*s2*v1*w0 + 
   2*c1*d1*pow(j2,2)*v2*w0 - 2*c1*d1*i2*o2*v2*w0 + 2*c1*l2*n1*o2*v2*w0 - 
   2*c1*j2*n1*q2*v2*w0 - 2*d1*j2*l2*r1*v2*w0 + 2*d1*i2*q2*r1*v2*w0 - 
   2*l2*n1*q2*r1*v2*w0 - 2*c1*j2*l2*s1*v2*w0 + 2*c1*i2*q2*s1*v2*w0 + 
   2*pow(l2,2)*r1*s1*v2*w0 - 2*d1*e0*k2*l2*o2*w1 - 2*d0*e1*k2*l2*o2*w1 + 
   2*c1*f0*k2*l2*o2*w1 + 2*c0*f1*k2*l2*o2*w1 + 2*c1*d0*k2*m2*o2*w1 + 
   2*c0*d1*k2*m2*o2*w1 - 4*c1*c2*k2*n0*o2*w1 - 4*c0*c2*k2*n1*o2*w1 - 
   4*c0*c1*k2*n2*o2*w1 + 2*d1*e0*j2*l2*p2*w1 + 2*d0*e1*j2*l2*p2*w1 - 
   2*c1*f0*j2*l2*p2*w1 - 2*c0*f1*j2*l2*p2*w1 - 2*c1*d0*j2*m2*p2*w1 - 
   2*c0*d1*j2*m2*p2*w1 + 4*c1*c2*j2*n0*p2*w1 + 4*c0*c2*j2*n1*p2*w1 + 
   4*c0*c1*j2*n2*p2*w1 + 2*d1*e0*j2*k2*q2*w1 + 2*d0*e1*j2*k2*q2*w1 - 
   2*c1*f0*j2*k2*q2*w1 - 2*c0*f1*j2*k2*q2*w1 - 2*d1*e0*i2*p2*q2*w1 - 
   2*d0*e1*i2*p2*q2*w1 + 2*c1*f0*i2*p2*q2*w1 + 2*c0*f1*i2*p2*q2*w1 + 
   2*e1*l2*n0*p2*q2*w1 - 2*c1*m2*n0*p2*q2*w1 + 2*e0*l2*n1*p2*q2*w1 - 
   2*c0*m2*n1*p2*q2*w1 - 2*e1*k2*n0*pow(q2,2)*w1 - 
   2*e0*k2*n1*pow(q2,2)*w1 - 2*c2*d1*j2*k2*r0*w1 - 2*c1*d2*j2*k2*r0*w1 + 
   2*c2*d1*i2*p2*r0*w1 + 2*c1*d2*i2*p2*r0*w1 + 2*f1*pow(l2,2)*p2*r0*w1 - 
   2*d1*l2*m2*p2*r0*w1 - 2*c2*l2*n1*p2*r0*w1 - 2*c1*l2*n2*p2*r0*w1 - 
   2*f1*k2*l2*q2*r0*w1 - 2*d1*k2*m2*q2*r0*w1 + 4*c2*k2*n1*q2*r0*w1 + 
   4*c1*k2*n2*q2*r0*w1 - 2*c2*d0*j2*k2*r1*w1 - 2*c0*d2*j2*k2*r1*w1 + 
   2*c2*d0*i2*p2*r1*w1 + 2*c0*d2*i2*p2*r1*w1 + 2*f0*pow(l2,2)*p2*r1*w1 - 
   2*d0*l2*m2*p2*r1*w1 - 2*c2*l2*n0*p2*r1*w1 - 2*c0*l2*n2*p2*r1*w1 - 
   2*f0*k2*l2*q2*r1*w1 - 2*d0*k2*m2*q2*r1*w1 + 4*c2*k2*n0*q2*r1*w1 + 
   4*c0*k2*n2*q2*r1*w1 + 4*d2*k2*l2*r0*r1*w1 - 2*c1*d0*j2*k2*r2*w1 - 
   2*c0*d1*j2*k2*r2*w1 + 2*c1*d0*i2*p2*r2*w1 + 2*c0*d1*i2*p2*r2*w1 - 
   2*c1*l2*n0*p2*r2*w1 - 2*c0*l2*n1*p2*r2*w1 + 4*c1*k2*n0*q2*r2*w1 + 
   4*c0*k2*n1*q2*r2*w1 + 4*d1*k2*l2*r0*r2*w1 + 4*d0*k2*l2*r1*r2*w1 + 
   4*c1*c2*j2*k2*s0*w1 - 4*c1*c2*i2*p2*s0*w1 - 2*e1*pow(l2,2)*p2*s0*w1 + 
   4*c1*l2*m2*p2*s0*w1 + 2*e1*k2*l2*q2*s0*w1 - 2*c1*k2*m2*q2*s0*w1 - 
   2*c2*k2*l2*r1*s0*w1 - 2*c1*k2*l2*r2*s0*w1 + 4*c0*c2*j2*k2*s1*w1 - 
   4*c0*c2*i2*p2*s1*w1 - 2*e0*pow(l2,2)*p2*s1*w1 + 4*c0*l2*m2*p2*s1*w1 + 
   2*e0*k2*l2*q2*s1*w1 - 2*c0*k2*m2*q2*s1*w1 - 2*c2*k2*l2*r0*s1*w1 - 
   2*c0*k2*l2*r2*s1*w1 + 4*c0*c1*j2*k2*s2*w1 - 4*c0*c1*i2*p2*s2*w1 - 
   2*c1*k2*l2*r0*s2*w1 - 2*c0*k2*l2*r1*s2*w1 - 2*d2*e1*pow(j2,2)*u0*w1 - 
   2*d1*e2*pow(j2,2)*u0*w1 + 2*c2*f1*pow(j2,2)*u0*w1 + 
   2*c1*f2*pow(j2,2)*u0*w1 + 2*d2*e1*i2*o2*u0*w1 + 2*d1*e2*i2*o2*u0*w1 - 
   2*c2*f1*i2*o2*u0*w1 - 2*c1*f2*i2*o2*u0*w1 + 2*f1*l2*m2*o2*u0*w1 - 
   2*d1*pow(m2,2)*o2*u0*w1 - 2*e2*l2*n1*o2*u0*w1 + 2*c2*m2*n1*o2*u0*w1 - 
   2*e1*l2*n2*o2*u0*w1 + 2*c1*m2*n2*o2*u0*w1 - 2*f1*j2*m2*q2*u0*w1 + 
   2*e2*j2*n1*q2*u0*w1 + 2*e1*j2*n2*q2*u0*w1 - 2*f2*j2*l2*r1*u0*w1 + 
   4*d2*j2*m2*r1*u0*w1 - 2*c2*j2*n2*r1*u0*w1 + 2*f2*i2*q2*r1*u0*w1 - 
   2*m2*n2*q2*r1*u0*w1 - 2*f1*j2*l2*r2*u0*w1 + 4*d1*j2*m2*r2*u0*w1 - 
   2*c2*j2*n1*r2*u0*w1 - 2*c1*j2*n2*r2*u0*w1 + 2*f1*i2*q2*r2*u0*w1 - 
   2*m2*n1*q2*r2*u0*w1 - 4*d2*i2*r1*r2*u0*w1 + 4*l2*n2*r1*r2*u0*w1 - 
   2*d1*i2*pow(r2,2)*u0*w1 + 2*l2*n1*pow(r2,2)*u0*w1 + 
   2*e2*j2*l2*s1*u0*w1 - 2*c2*j2*m2*s1*u0*w1 - 2*e2*i2*q2*s1*u0*w1 + 
   2*pow(m2,2)*q2*s1*u0*w1 + 2*c2*i2*r2*s1*u0*w1 - 2*l2*m2*r2*s1*u0*w1 + 
   2*e1*j2*l2*s2*u0*w1 - 2*c1*j2*m2*s2*u0*w1 - 2*e1*i2*q2*s2*u0*w1 + 
   2*c2*i2*r1*s2*u0*w1 - 2*l2*m2*r1*s2*u0*w1 + 2*c1*i2*r2*s2*u0*w1 - 
   2*d2*e0*pow(j2,2)*u1*w1 - 2*d0*e2*pow(j2,2)*u1*w1 + 
   2*c2*f0*pow(j2,2)*u1*w1 + 2*c0*f2*pow(j2,2)*u1*w1 + 
   2*d2*e0*i2*o2*u1*w1 + 2*d0*e2*i2*o2*u1*w1 - 2*c2*f0*i2*o2*u1*w1 - 
   2*c0*f2*i2*o2*u1*w1 + 2*f0*l2*m2*o2*u1*w1 - 2*d0*pow(m2,2)*o2*u1*w1 - 
   2*e2*l2*n0*o2*u1*w1 + 2*c2*m2*n0*o2*u1*w1 - 2*e0*l2*n2*o2*u1*w1 + 
   2*c0*m2*n2*o2*u1*w1 - 2*f0*j2*m2*q2*u1*w1 + 2*e2*j2*n0*q2*u1*w1 + 
   2*e0*j2*n2*q2*u1*w1 - 2*f2*j2*l2*r0*u1*w1 + 4*d2*j2*m2*r0*u1*w1 - 
   2*c2*j2*n2*r0*u1*w1 + 2*f2*i2*q2*r0*u1*w1 - 2*m2*n2*q2*r0*u1*w1 - 
   2*f0*j2*l2*r2*u1*w1 + 4*d0*j2*m2*r2*u1*w1 - 2*c2*j2*n0*r2*u1*w1 - 
   2*c0*j2*n2*r2*u1*w1 + 2*f0*i2*q2*r2*u1*w1 - 2*m2*n0*q2*r2*u1*w1 - 
   4*d2*i2*r0*r2*u1*w1 + 4*l2*n2*r0*r2*u1*w1 - 2*d0*i2*pow(r2,2)*u1*w1 + 
   2*l2*n0*pow(r2,2)*u1*w1 + 2*e2*j2*l2*s0*u1*w1 - 2*c2*j2*m2*s0*u1*w1 - 
   2*e2*i2*q2*s0*u1*w1 + 2*pow(m2,2)*q2*s0*u1*w1 + 2*c2*i2*r2*s0*u1*w1 - 
   2*l2*m2*r2*s0*u1*w1 + 2*e0*j2*l2*s2*u1*w1 - 2*c0*j2*m2*s2*u1*w1 - 
   2*e0*i2*q2*s2*u1*w1 + 2*c2*i2*r0*s2*u1*w1 - 2*l2*m2*r0*s2*u1*w1 + 
   2*c0*i2*r2*s2*u1*w1 - 2*d1*e0*pow(j2,2)*u2*w1 - 
   2*d0*e1*pow(j2,2)*u2*w1 + 2*c1*f0*pow(j2,2)*u2*w1 + 
   2*c0*f1*pow(j2,2)*u2*w1 + 2*d1*e0*i2*o2*u2*w1 + 2*d0*e1*i2*o2*u2*w1 - 
   2*c1*f0*i2*o2*u2*w1 - 2*c0*f1*i2*o2*u2*w1 - 2*e1*l2*n0*o2*u2*w1 + 
   2*c1*m2*n0*o2*u2*w1 - 2*e0*l2*n1*o2*u2*w1 + 2*c0*m2*n1*o2*u2*w1 + 
   2*e1*j2*n0*q2*u2*w1 + 2*e0*j2*n1*q2*u2*w1 - 2*f1*j2*l2*r0*u2*w1 + 
   4*d1*j2*m2*r0*u2*w1 - 2*c2*j2*n1*r0*u2*w1 - 2*c1*j2*n2*r0*u2*w1 + 
   2*f1*i2*q2*r0*u2*w1 - 2*m2*n1*q2*r0*u2*w1 - 2*f0*j2*l2*r1*u2*w1 + 
   4*d0*j2*m2*r1*u2*w1 - 2*c2*j2*n0*r1*u2*w1 - 2*c0*j2*n2*r1*u2*w1 + 
   2*f0*i2*q2*r1*u2*w1 - 2*m2*n0*q2*r1*u2*w1 - 4*d2*i2*r0*r1*u2*w1 + 
   4*l2*n2*r0*r1*u2*w1 - 2*c1*j2*n0*r2*u2*w1 - 2*c0*j2*n1*r2*u2*w1 - 
   4*d1*i2*r0*r2*u2*w1 + 4*l2*n1*r0*r2*u2*w1 - 4*d0*i2*r1*r2*u2*w1 + 
   4*l2*n0*r1*r2*u2*w1 + 2*e1*j2*l2*s0*u2*w1 - 2*c1*j2*m2*s0*u2*w1 - 
   2*e1*i2*q2*s0*u2*w1 + 2*c2*i2*r1*s0*u2*w1 - 2*l2*m2*r1*s0*u2*w1 + 
   2*c1*i2*r2*s0*u2*w1 + 2*e0*j2*l2*s1*u2*w1 - 2*c0*j2*m2*s1*u2*w1 - 
   2*e0*i2*q2*s1*u2*w1 + 2*c2*i2*r0*s1*u2*w1 - 2*l2*m2*r0*s1*u2*w1 + 
   2*c0*i2*r2*s1*u2*w1 + 2*c1*i2*r0*s2*u2*w1 + 2*c0*i2*r1*s2*u2*w1 + 
   2*c2*d1*pow(j2,2)*v0*w1 + 2*c1*d2*pow(j2,2)*v0*w1 - 
   2*c2*d1*i2*o2*v0*w1 - 2*c1*d2*i2*o2*v0*w1 - 2*f1*pow(l2,2)*o2*v0*w1 + 
   2*d1*l2*m2*o2*v0*w1 + 2*c2*l2*n1*o2*v0*w1 + 2*c1*l2*n2*o2*v0*w1 + 
   4*f1*j2*l2*q2*v0*w1 - 2*d1*j2*m2*q2*v0*w1 - 2*c2*j2*n1*q2*v0*w1 - 
   2*c1*j2*n2*q2*v0*w1 - 2*f1*i2*pow(q2,2)*v0*w1 + 
   2*m2*n1*pow(q2,2)*v0*w1 - 2*d2*j2*l2*r1*v0*w1 + 2*d2*i2*q2*r1*v0*w1 - 
   2*l2*n2*q2*r1*v0*w1 - 2*d1*j2*l2*r2*v0*w1 + 2*d1*i2*q2*r2*v0*w1 - 
   2*l2*n1*q2*r2*v0*w1 - 2*c2*j2*l2*s1*v0*w1 + 2*c2*i2*q2*s1*v0*w1 - 
   2*l2*m2*q2*s1*v0*w1 + 2*pow(l2,2)*r2*s1*v0*w1 - 2*c1*j2*l2*s2*v0*w1 + 
   2*c1*i2*q2*s2*v0*w1 + 2*pow(l2,2)*r1*s2*v0*w1 + 
   2*c2*d0*pow(j2,2)*v1*w1 + 2*c0*d2*pow(j2,2)*v1*w1 - 
   2*c2*d0*i2*o2*v1*w1 - 2*c0*d2*i2*o2*v1*w1 - 2*f0*pow(l2,2)*o2*v1*w1 + 
   2*d0*l2*m2*o2*v1*w1 + 2*c2*l2*n0*o2*v1*w1 + 2*c0*l2*n2*o2*v1*w1 + 
   4*f0*j2*l2*q2*v1*w1 - 2*d0*j2*m2*q2*v1*w1 - 2*c2*j2*n0*q2*v1*w1 - 
   2*c0*j2*n2*q2*v1*w1 - 2*f0*i2*pow(q2,2)*v1*w1 + 
   2*m2*n0*pow(q2,2)*v1*w1 - 2*d2*j2*l2*r0*v1*w1 + 2*d2*i2*q2*r0*v1*w1 - 
   2*l2*n2*q2*r0*v1*w1 - 2*d0*j2*l2*r2*v1*w1 + 2*d0*i2*q2*r2*v1*w1 - 
   2*l2*n0*q2*r2*v1*w1 - 2*c2*j2*l2*s0*v1*w1 + 2*c2*i2*q2*s0*v1*w1 - 
   2*l2*m2*q2*s0*v1*w1 + 2*pow(l2,2)*r2*s0*v1*w1 - 2*c0*j2*l2*s2*v1*w1 + 
   2*c0*i2*q2*s2*v1*w1 + 2*pow(l2,2)*r0*s2*v1*w1 + 
   2*c1*d0*pow(j2,2)*v2*w1 + 2*c0*d1*pow(j2,2)*v2*w1 - 
   2*c1*d0*i2*o2*v2*w1 - 2*c0*d1*i2*o2*v2*w1 + 2*c1*l2*n0*o2*v2*w1 + 
   2*c0*l2*n1*o2*v2*w1 - 2*c1*j2*n0*q2*v2*w1 - 2*c0*j2*n1*q2*v2*w1 - 
   2*d1*j2*l2*r0*v2*w1 + 2*d1*i2*q2*r0*v2*w1 - 2*l2*n1*q2*r0*v2*w1 - 
   2*d0*j2*l2*r1*v2*w1 + 2*d0*i2*q2*r1*v2*w1 - 2*l2*n0*q2*r1*v2*w1 - 
   2*c1*j2*l2*s0*v2*w1 + 2*c1*i2*q2*s0*v2*w1 + 2*pow(l2,2)*r1*s0*v2*w1 - 
   2*c0*j2*l2*s1*v2*w1 + 2*c0*i2*q2*s1*v2*w1 + 2*pow(l2,2)*r0*s1*v2*w1 - 
   4*c1*c2*pow(j2,2)*w0*w1 + 4*c1*c2*i2*o2*w0*w1 + 
   2*e1*pow(l2,2)*o2*w0*w1 - 4*c1*l2*m2*o2*w0*w1 - 4*e1*j2*l2*q2*w0*w1 + 
   4*c1*j2*m2*q2*w0*w1 + 2*e1*i2*pow(q2,2)*w0*w1 + 4*c2*j2*l2*r1*w0*w1 - 
   4*c2*i2*q2*r1*w0*w1 + 4*l2*m2*q2*r1*w0*w1 + 4*c1*j2*l2*r2*w0*w1 - 
   4*c1*i2*q2*r2*w0*w1 - 4*pow(l2,2)*r1*r2*w0*w1 - 
   2*c0*c2*pow(j2,2)*pow(w1,2) + 2*c0*c2*i2*o2*pow(w1,2) + 
   e0*pow(l2,2)*o2*pow(w1,2) - 2*c0*l2*m2*o2*pow(w1,2) - 
   2*e0*j2*l2*q2*pow(w1,2) + 2*c0*j2*m2*q2*pow(w1,2) + 
   e0*i2*pow(q2,2)*pow(w1,2) + 2*c2*j2*l2*r0*pow(w1,2) - 
   2*c2*i2*q2*r0*pow(w1,2) + 2*l2*m2*q2*r0*pow(w1,2) + 
   2*c0*j2*l2*r2*pow(w1,2) - 2*c0*i2*q2*r2*pow(w1,2) - 
   2*pow(l2,2)*r0*r2*pow(w1,2) - 2*pow(c1,2)*k2*n0*o2*w2 - 
   4*c0*c1*k2*n1*o2*w2 + 2*pow(c1,2)*j2*n0*p2*w2 + 4*c0*c1*j2*n1*p2*w2 - 
   2*c1*d1*j2*k2*r0*w2 + 2*c1*d1*i2*p2*r0*w2 - 2*c1*l2*n1*p2*r0*w2 + 
   4*c1*k2*n1*q2*r0*w2 - 2*c1*d0*j2*k2*r1*w2 - 2*c0*d1*j2*k2*r1*w2 + 
   2*c1*d0*i2*p2*r1*w2 + 2*c0*d1*i2*p2*r1*w2 - 2*c1*l2*n0*p2*r1*w2 - 
   2*c0*l2*n1*p2*r1*w2 + 4*c1*k2*n0*q2*r1*w2 + 4*c0*k2*n1*q2*r1*w2 + 
   4*d1*k2*l2*r0*r1*w2 + 2*d0*k2*l2*pow(r1,2)*w2 + 
   2*pow(c1,2)*j2*k2*s0*w2 - 2*pow(c1,2)*i2*p2*s0*w2 - 
   2*c1*k2*l2*r1*s0*w2 + 4*c0*c1*j2*k2*s1*w2 - 4*c0*c1*i2*p2*s1*w2 - 
   2*c1*k2*l2*r0*s1*w2 - 2*c0*k2*l2*r1*s1*w2 - 2*d1*e1*pow(j2,2)*u0*w2 + 
   2*c1*f1*pow(j2,2)*u0*w2 + 2*d1*e1*i2*o2*u0*w2 - 2*c1*f1*i2*o2*u0*w2 - 
   2*e1*l2*n1*o2*u0*w2 + 2*c1*m2*n1*o2*u0*w2 + 2*e1*j2*n1*q2*u0*w2 - 
   2*f1*j2*l2*r1*u0*w2 + 4*d1*j2*m2*r1*u0*w2 - 2*c2*j2*n1*r1*u0*w2 - 
   2*c1*j2*n2*r1*u0*w2 + 2*f1*i2*q2*r1*u0*w2 - 2*m2*n1*q2*r1*u0*w2 - 
   2*d2*i2*pow(r1,2)*u0*w2 + 2*l2*n2*pow(r1,2)*u0*w2 - 
   2*c1*j2*n1*r2*u0*w2 - 4*d1*i2*r1*r2*u0*w2 + 4*l2*n1*r1*r2*u0*w2 + 
   2*e1*j2*l2*s1*u0*w2 - 2*c1*j2*m2*s1*u0*w2 - 2*e1*i2*q2*s1*u0*w2 + 
   2*c2*i2*r1*s1*u0*w2 - 2*l2*m2*r1*s1*u0*w2 + 2*c1*i2*r2*s1*u0*w2 + 
   2*c1*i2*r1*s2*u0*w2 - 2*d1*e0*pow(j2,2)*u1*w2 - 
   2*d0*e1*pow(j2,2)*u1*w2 + 2*c1*f0*pow(j2,2)*u1*w2 + 
   2*c0*f1*pow(j2,2)*u1*w2 + 2*d1*e0*i2*o2*u1*w2 + 2*d0*e1*i2*o2*u1*w2 - 
   2*c1*f0*i2*o2*u1*w2 - 2*c0*f1*i2*o2*u1*w2 - 2*e1*l2*n0*o2*u1*w2 + 
   2*c1*m2*n0*o2*u1*w2 - 2*e0*l2*n1*o2*u1*w2 + 2*c0*m2*n1*o2*u1*w2 + 
   2*e1*j2*n0*q2*u1*w2 + 2*e0*j2*n1*q2*u1*w2 - 2*f1*j2*l2*r0*u1*w2 + 
   4*d1*j2*m2*r0*u1*w2 - 2*c2*j2*n1*r0*u1*w2 - 2*c1*j2*n2*r0*u1*w2 + 
   2*f1*i2*q2*r0*u1*w2 - 2*m2*n1*q2*r0*u1*w2 - 2*f0*j2*l2*r1*u1*w2 + 
   4*d0*j2*m2*r1*u1*w2 - 2*c2*j2*n0*r1*u1*w2 - 2*c0*j2*n2*r1*u1*w2 + 
   2*f0*i2*q2*r1*u1*w2 - 2*m2*n0*q2*r1*u1*w2 - 4*d2*i2*r0*r1*u1*w2 + 
   4*l2*n2*r0*r1*u1*w2 - 2*c1*j2*n0*r2*u1*w2 - 2*c0*j2*n1*r2*u1*w2 - 
   4*d1*i2*r0*r2*u1*w2 + 4*l2*n1*r0*r2*u1*w2 - 4*d0*i2*r1*r2*u1*w2 + 
   4*l2*n0*r1*r2*u1*w2 + 2*e1*j2*l2*s0*u1*w2 - 2*c1*j2*m2*s0*u1*w2 - 
   2*e1*i2*q2*s0*u1*w2 + 2*c2*i2*r1*s0*u1*w2 - 2*l2*m2*r1*s0*u1*w2 + 
   2*c1*i2*r2*s0*u1*w2 + 2*e0*j2*l2*s1*u1*w2 - 2*c0*j2*m2*s1*u1*w2 - 
   2*e0*i2*q2*s1*u1*w2 + 2*c2*i2*r0*s1*u1*w2 - 2*l2*m2*r0*s1*u1*w2 + 
   2*c0*i2*r2*s1*u1*w2 + 2*c1*i2*r0*s2*u1*w2 + 2*c0*i2*r1*s2*u1*w2 - 
   2*c1*j2*n1*r0*u2*w2 - 2*c1*j2*n0*r1*u2*w2 - 2*c0*j2*n1*r1*u2*w2 - 
   4*d1*i2*r0*r1*u2*w2 + 4*l2*n1*r0*r1*u2*w2 - 2*d0*i2*pow(r1,2)*u2*w2 + 
   2*l2*n0*pow(r1,2)*u2*w2 + 2*c1*i2*r1*s0*u2*w2 + 2*c1*i2*r0*s1*u2*w2 + 
   2*c0*i2*r1*s1*u2*w2 + 2*c1*d1*pow(j2,2)*v0*w2 - 2*c1*d1*i2*o2*v0*w2 + 
   2*c1*l2*n1*o2*v0*w2 - 2*c1*j2*n1*q2*v0*w2 - 2*d1*j2*l2*r1*v0*w2 + 
   2*d1*i2*q2*r1*v0*w2 - 2*l2*n1*q2*r1*v0*w2 - 2*c1*j2*l2*s1*v0*w2 + 
   2*c1*i2*q2*s1*v0*w2 + 2*pow(l2,2)*r1*s1*v0*w2 + 
   2*c1*d0*pow(j2,2)*v1*w2 + 2*c0*d1*pow(j2,2)*v1*w2 - 
   2*c1*d0*i2*o2*v1*w2 - 2*c0*d1*i2*o2*v1*w2 + 2*c1*l2*n0*o2*v1*w2 + 
   2*c0*l2*n1*o2*v1*w2 - 2*c1*j2*n0*q2*v1*w2 - 2*c0*j2*n1*q2*v1*w2 - 
   2*d1*j2*l2*r0*v1*w2 + 2*d1*i2*q2*r0*v1*w2 - 2*l2*n1*q2*r0*v1*w2 - 
   2*d0*j2*l2*r1*v1*w2 + 2*d0*i2*q2*r1*v1*w2 - 2*l2*n0*q2*r1*v1*w2 - 
   2*c1*j2*l2*s0*v1*w2 + 2*c1*i2*q2*s0*v1*w2 + 2*pow(l2,2)*r1*s0*v1*w2 - 
   2*c0*j2*l2*s1*v1*w2 + 2*c0*i2*q2*s1*v1*w2 + 2*pow(l2,2)*r0*s1*v1*w2 - 
   2*pow(c1,2)*pow(j2,2)*w0*w2 + 2*pow(c1,2)*i2*o2*w0*w2 + 
   4*c1*j2*l2*r1*w0*w2 - 4*c1*i2*q2*r1*w0*w2 - 
   2*pow(l2,2)*pow(r1,2)*w0*w2 - 4*c0*c1*pow(j2,2)*w1*w2 + 
   4*c0*c1*i2*o2*w1*w2 + 4*c1*j2*l2*r0*w1*w2 - 4*c1*i2*q2*r0*w1*w2 + 
   4*c0*j2*l2*r1*w1*w2 - 4*c0*i2*q2*r1*w1*w2 - 4*pow(l2,2)*r0*r1*w1*w2 + 
   pow(f1,2)*pow(k2,2)*o2*z0 - e1*g1*pow(k2,2)*o2*z0 - 
   2*pow(f1,2)*j2*k2*p2*z0 + 2*e1*g1*j2*k2*p2*z0 + 
   pow(f1,2)*i2*pow(p2,2)*z0 - e1*g1*i2*pow(p2,2)*z0 - 
   2*f1*m2*n1*pow(p2,2)*z0 + e2*pow(n1,2)*pow(p2,2)*z0 + 
   2*e1*n1*n2*pow(p2,2)*z0 - 2*g1*k2*m2*p2*r1*z0 + 2*f2*k2*n1*p2*r1*z0 + 
   2*f1*k2*n2*p2*r1*z0 + g2*pow(k2,2)*pow(r1,2)*z0 + 
   2*f1*k2*n1*p2*r2*z0 + 2*g1*pow(k2,2)*r1*r2*z0 + 2*f1*k2*m2*p2*s1*z0 - 
   2*e2*k2*n1*p2*s1*z0 - 2*e1*k2*n2*p2*s1*z0 - 2*f2*pow(k2,2)*r1*s1*z0 - 
   2*f1*pow(k2,2)*r2*s1*z0 + e2*pow(k2,2)*pow(s1,2)*z0 - 
   2*e1*k2*n1*p2*s2*z0 - 2*f1*pow(k2,2)*r1*s2*z0 + 
   2*e1*pow(k2,2)*s1*s2*z0 + pow(f1,2)*pow(j2,2)*t2*z0 - 
   e1*g1*pow(j2,2)*t2*z0 - pow(f1,2)*i2*o2*t2*z0 + e1*g1*i2*o2*t2*z0 + 
   2*f1*m2*n1*o2*t2*z0 - e2*pow(n1,2)*o2*t2*z0 - 2*e1*n1*n2*o2*t2*z0 + 
   2*g1*j2*m2*r1*t2*z0 - 2*f2*j2*n1*r1*t2*z0 - 2*f1*j2*n2*r1*t2*z0 - 
   g2*i2*pow(r1,2)*t2*z0 + pow(n2,2)*pow(r1,2)*t2*z0 - 
   2*f1*j2*n1*r2*t2*z0 - 2*g1*i2*r1*r2*t2*z0 + 4*n1*n2*r1*r2*t2*z0 + 
   pow(n1,2)*pow(r2,2)*t2*z0 - 2*f1*j2*m2*s1*t2*z0 + 
   2*e2*j2*n1*s1*t2*z0 + 2*e1*j2*n2*s1*t2*z0 + 2*f2*i2*r1*s1*t2*z0 - 
   2*m2*n2*r1*s1*t2*z0 + 2*f1*i2*r2*s1*t2*z0 - 2*m2*n1*r2*s1*t2*z0 - 
   e2*i2*pow(s1,2)*t2*z0 + pow(m2,2)*pow(s1,2)*t2*z0 + 
   2*e1*j2*n1*s2*t2*z0 + 2*f1*i2*r1*s2*t2*z0 - 2*m2*n1*r1*s2*t2*z0 - 
   2*e1*i2*s1*s2*t2*z0 + 2*g1*k2*m2*o2*v1*z0 - 2*f2*k2*n1*o2*v1*z0 - 
   2*f1*k2*n2*o2*v1*z0 - 2*g1*j2*m2*p2*v1*z0 + 2*f2*j2*n1*p2*v1*z0 + 
   2*f1*j2*n2*p2*v1*z0 - 2*g2*j2*k2*r1*v1*z0 + 2*g2*i2*p2*r1*v1*z0 - 
   2*pow(n2,2)*p2*r1*v1*z0 - 2*g1*j2*k2*r2*v1*z0 + 2*g1*i2*p2*r2*v1*z0 - 
   4*n1*n2*p2*r2*v1*z0 + 2*f2*j2*k2*s1*v1*z0 - 2*f2*i2*p2*s1*v1*z0 + 
   2*m2*n2*p2*s1*v1*z0 + 2*k2*n2*r2*s1*v1*z0 + 2*f1*j2*k2*s2*v1*z0 - 
   2*f1*i2*p2*s2*v1*z0 + 2*m2*n1*p2*s2*v1*z0 + 2*k2*n2*r1*s2*v1*z0 + 
   2*k2*n1*r2*s2*v1*z0 - 4*k2*m2*s1*s2*v1*z0 + 
   g2*pow(j2,2)*pow(v1,2)*z0 - g2*i2*o2*pow(v1,2)*z0 + 
   pow(n2,2)*o2*pow(v1,2)*z0 - 2*j2*n2*s2*pow(v1,2)*z0 + 
   i2*pow(s2,2)*pow(v1,2)*z0 - 2*f1*k2*n1*o2*v2*z0 + 
   2*f1*j2*n1*p2*v2*z0 - 2*g1*j2*k2*r1*v2*z0 + 2*g1*i2*p2*r1*v2*z0 - 
   4*n1*n2*p2*r1*v2*z0 - 2*pow(n1,2)*p2*r2*v2*z0 + 2*f1*j2*k2*s1*v2*z0 - 
   2*f1*i2*p2*s1*v2*z0 + 2*m2*n1*p2*s1*v2*z0 + 2*k2*n2*r1*s1*v2*z0 + 
   2*k2*n1*r2*s1*v2*z0 - 2*k2*m2*pow(s1,2)*v2*z0 + 2*k2*n1*r1*s2*v2*z0 + 
   2*g1*pow(j2,2)*v1*v2*z0 - 2*g1*i2*o2*v1*v2*z0 + 4*n1*n2*o2*v1*v2*z0 - 
   4*j2*n2*s1*v1*v2*z0 - 4*j2*n1*s2*v1*v2*z0 + 4*i2*s1*s2*v1*v2*z0 + 
   pow(n1,2)*o2*pow(v2,2)*z0 - 2*j2*n1*s1*pow(v2,2)*z0 + 
   i2*pow(s1,2)*pow(v2,2)*z0 - 2*f1*k2*m2*o2*w1*z0 + 
   2*e2*k2*n1*o2*w1*z0 + 2*e1*k2*n2*o2*w1*z0 + 2*f1*j2*m2*p2*w1*z0 - 
   2*e2*j2*n1*p2*w1*z0 - 2*e1*j2*n2*p2*w1*z0 + 2*f2*j2*k2*r1*w1*z0 - 
   2*f2*i2*p2*r1*w1*z0 + 2*m2*n2*p2*r1*w1*z0 + 2*f1*j2*k2*r2*w1*z0 - 
   2*f1*i2*p2*r2*w1*z0 + 2*m2*n1*p2*r2*w1*z0 - 4*k2*n2*r1*r2*w1*z0 - 
   2*k2*n1*pow(r2,2)*w1*z0 - 2*e2*j2*k2*s1*w1*z0 + 2*e2*i2*p2*s1*w1*z0 - 
   2*pow(m2,2)*p2*s1*w1*z0 + 2*k2*m2*r2*s1*w1*z0 - 2*e1*j2*k2*s2*w1*z0 + 
   2*e1*i2*p2*s2*w1*z0 + 2*k2*m2*r1*s2*w1*z0 - 2*f2*pow(j2,2)*v1*w1*z0 + 
   2*f2*i2*o2*v1*w1*z0 - 2*m2*n2*o2*v1*w1*z0 + 2*j2*n2*r2*v1*w1*z0 + 
   2*j2*m2*s2*v1*w1*z0 - 2*i2*r2*s2*v1*w1*z0 - 2*f1*pow(j2,2)*v2*w1*z0 + 
   2*f1*i2*o2*v2*w1*z0 - 2*m2*n1*o2*v2*w1*z0 + 2*j2*n2*r1*v2*w1*z0 + 
   2*j2*n1*r2*v2*w1*z0 + 2*j2*m2*s1*v2*w1*z0 - 2*i2*r2*s1*v2*w1*z0 - 
   2*i2*r1*s2*v2*w1*z0 + e2*pow(j2,2)*pow(w1,2)*z0 - 
   e2*i2*o2*pow(w1,2)*z0 + pow(m2,2)*o2*pow(w1,2)*z0 - 
   2*j2*m2*r2*pow(w1,2)*z0 + i2*pow(r2,2)*pow(w1,2)*z0 + 
   2*e1*k2*n1*o2*w2*z0 - 2*e1*j2*n1*p2*w2*z0 + 2*f1*j2*k2*r1*w2*z0 - 
   2*f1*i2*p2*r1*w2*z0 + 2*m2*n1*p2*r1*w2*z0 - 2*k2*n2*pow(r1,2)*w2*z0 - 
   4*k2*n1*r1*r2*w2*z0 - 2*e1*j2*k2*s1*w2*z0 + 2*e1*i2*p2*s1*w2*z0 + 
   2*k2*m2*r1*s1*w2*z0 - 2*f1*pow(j2,2)*v1*w2*z0 + 2*f1*i2*o2*v1*w2*z0 - 
   2*m2*n1*o2*v1*w2*z0 + 2*j2*n2*r1*v1*w2*z0 + 2*j2*n1*r2*v1*w2*z0 + 
   2*j2*m2*s1*v1*w2*z0 - 2*i2*r2*s1*v1*w2*z0 - 2*i2*r1*s2*v1*w2*z0 + 
   2*j2*n1*r1*v2*w2*z0 - 2*i2*r1*s1*v2*w2*z0 + 2*e1*pow(j2,2)*w1*w2*z0 - 
   2*e1*i2*o2*w1*w2*z0 - 4*j2*m2*r1*w1*w2*z0 + 4*i2*r1*r2*w1*w2*z0 + 
   i2*pow(r1,2)*pow(w2,2)*z0 + 2*f0*f1*pow(k2,2)*o2*z1 - 
   e1*g0*pow(k2,2)*o2*z1 - e0*g1*pow(k2,2)*o2*z1 - 4*f0*f1*j2*k2*p2*z1 + 
   2*e1*g0*j2*k2*p2*z1 + 2*e0*g1*j2*k2*p2*z1 + 2*f0*f1*i2*pow(p2,2)*z1 - 
   e1*g0*i2*pow(p2,2)*z1 - e0*g1*i2*pow(p2,2)*z1 - 
   2*f1*m2*n0*pow(p2,2)*z1 - 2*f0*m2*n1*pow(p2,2)*z1 + 
   2*e2*n0*n1*pow(p2,2)*z1 + 2*e1*n0*n2*pow(p2,2)*z1 + 
   2*e0*n1*n2*pow(p2,2)*z1 - 2*g1*k2*m2*p2*r0*z1 + 2*f2*k2*n1*p2*r0*z1 + 
   2*f1*k2*n2*p2*r0*z1 - 2*g0*k2*m2*p2*r1*z1 + 2*f2*k2*n0*p2*r1*z1 + 
   2*f0*k2*n2*p2*r1*z1 + 2*g2*pow(k2,2)*r0*r1*z1 + 2*f1*k2*n0*p2*r2*z1 + 
   2*f0*k2*n1*p2*r2*z1 + 2*g1*pow(k2,2)*r0*r2*z1 + 
   2*g0*pow(k2,2)*r1*r2*z1 + 2*f1*k2*m2*p2*s0*z1 - 2*e2*k2*n1*p2*s0*z1 - 
   2*e1*k2*n2*p2*s0*z1 - 2*f2*pow(k2,2)*r1*s0*z1 - 
   2*f1*pow(k2,2)*r2*s0*z1 + 2*f0*k2*m2*p2*s1*z1 - 2*e2*k2*n0*p2*s1*z1 - 
   2*e0*k2*n2*p2*s1*z1 - 2*f2*pow(k2,2)*r0*s1*z1 - 
   2*f0*pow(k2,2)*r2*s1*z1 + 2*e2*pow(k2,2)*s0*s1*z1 - 
   2*e1*k2*n0*p2*s2*z1 - 2*e0*k2*n1*p2*s2*z1 - 2*f1*pow(k2,2)*r0*s2*z1 - 
   2*f0*pow(k2,2)*r1*s2*z1 + 2*e1*pow(k2,2)*s0*s2*z1 + 
   2*e0*pow(k2,2)*s1*s2*z1 + 2*f0*f1*pow(j2,2)*t2*z1 - 
   e1*g0*pow(j2,2)*t2*z1 - e0*g1*pow(j2,2)*t2*z1 - 2*f0*f1*i2*o2*t2*z1 + 
   e1*g0*i2*o2*t2*z1 + e0*g1*i2*o2*t2*z1 + 2*f1*m2*n0*o2*t2*z1 + 
   2*f0*m2*n1*o2*t2*z1 - 2*e2*n0*n1*o2*t2*z1 - 2*e1*n0*n2*o2*t2*z1 - 
   2*e0*n1*n2*o2*t2*z1 + 2*g1*j2*m2*r0*t2*z1 - 2*f2*j2*n1*r0*t2*z1 - 
   2*f1*j2*n2*r0*t2*z1 + 2*g0*j2*m2*r1*t2*z1 - 2*f2*j2*n0*r1*t2*z1 - 
   2*f0*j2*n2*r1*t2*z1 - 2*g2*i2*r0*r1*t2*z1 + 2*pow(n2,2)*r0*r1*t2*z1 - 
   2*f1*j2*n0*r2*t2*z1 - 2*f0*j2*n1*r2*t2*z1 - 2*g1*i2*r0*r2*t2*z1 + 
   4*n1*n2*r0*r2*t2*z1 - 2*g0*i2*r1*r2*t2*z1 + 4*n0*n2*r1*r2*t2*z1 + 
   2*n0*n1*pow(r2,2)*t2*z1 - 2*f1*j2*m2*s0*t2*z1 + 2*e2*j2*n1*s0*t2*z1 + 
   2*e1*j2*n2*s0*t2*z1 + 2*f2*i2*r1*s0*t2*z1 - 2*m2*n2*r1*s0*t2*z1 + 
   2*f1*i2*r2*s0*t2*z1 - 2*m2*n1*r2*s0*t2*z1 - 2*f0*j2*m2*s1*t2*z1 + 
   2*e2*j2*n0*s1*t2*z1 + 2*e0*j2*n2*s1*t2*z1 + 2*f2*i2*r0*s1*t2*z1 - 
   2*m2*n2*r0*s1*t2*z1 + 2*f0*i2*r2*s1*t2*z1 - 2*m2*n0*r2*s1*t2*z1 - 
   2*e2*i2*s0*s1*t2*z1 + 2*pow(m2,2)*s0*s1*t2*z1 + 2*e1*j2*n0*s2*t2*z1 + 
   2*e0*j2*n1*s2*t2*z1 + 2*f1*i2*r0*s2*t2*z1 - 2*m2*n1*r0*s2*t2*z1 + 
   2*f0*i2*r1*s2*t2*z1 - 2*m2*n0*r1*s2*t2*z1 - 2*e1*i2*s0*s2*t2*z1 - 
   2*e0*i2*s1*s2*t2*z1 + 2*g1*k2*m2*o2*v0*z1 - 2*f2*k2*n1*o2*v0*z1 - 
   2*f1*k2*n2*o2*v0*z1 - 2*g1*j2*m2*p2*v0*z1 + 2*f2*j2*n1*p2*v0*z1 + 
   2*f1*j2*n2*p2*v0*z1 - 2*g2*j2*k2*r1*v0*z1 + 2*g2*i2*p2*r1*v0*z1 - 
   2*pow(n2,2)*p2*r1*v0*z1 - 2*g1*j2*k2*r2*v0*z1 + 2*g1*i2*p2*r2*v0*z1 - 
   4*n1*n2*p2*r2*v0*z1 + 2*f2*j2*k2*s1*v0*z1 - 2*f2*i2*p2*s1*v0*z1 + 
   2*m2*n2*p2*s1*v0*z1 + 2*k2*n2*r2*s1*v0*z1 + 2*f1*j2*k2*s2*v0*z1 - 
   2*f1*i2*p2*s2*v0*z1 + 2*m2*n1*p2*s2*v0*z1 + 2*k2*n2*r1*s2*v0*z1 + 
   2*k2*n1*r2*s2*v0*z1 - 4*k2*m2*s1*s2*v0*z1 + 2*g0*k2*m2*o2*v1*z1 - 
   2*f2*k2*n0*o2*v1*z1 - 2*f0*k2*n2*o2*v1*z1 - 2*g0*j2*m2*p2*v1*z1 + 
   2*f2*j2*n0*p2*v1*z1 + 2*f0*j2*n2*p2*v1*z1 - 2*g2*j2*k2*r0*v1*z1 + 
   2*g2*i2*p2*r0*v1*z1 - 2*pow(n2,2)*p2*r0*v1*z1 - 2*g0*j2*k2*r2*v1*z1 + 
   2*g0*i2*p2*r2*v1*z1 - 4*n0*n2*p2*r2*v1*z1 + 2*f2*j2*k2*s0*v1*z1 - 
   2*f2*i2*p2*s0*v1*z1 + 2*m2*n2*p2*s0*v1*z1 + 2*k2*n2*r2*s0*v1*z1 + 
   2*f0*j2*k2*s2*v1*z1 - 2*f0*i2*p2*s2*v1*z1 + 2*m2*n0*p2*s2*v1*z1 + 
   2*k2*n2*r0*s2*v1*z1 + 2*k2*n0*r2*s2*v1*z1 - 4*k2*m2*s0*s2*v1*z1 + 
   2*g2*pow(j2,2)*v0*v1*z1 - 2*g2*i2*o2*v0*v1*z1 + 
   2*pow(n2,2)*o2*v0*v1*z1 - 4*j2*n2*s2*v0*v1*z1 + 
   2*i2*pow(s2,2)*v0*v1*z1 - 2*f1*k2*n0*o2*v2*z1 - 2*f0*k2*n1*o2*v2*z1 + 
   2*f1*j2*n0*p2*v2*z1 + 2*f0*j2*n1*p2*v2*z1 - 2*g1*j2*k2*r0*v2*z1 + 
   2*g1*i2*p2*r0*v2*z1 - 4*n1*n2*p2*r0*v2*z1 - 2*g0*j2*k2*r1*v2*z1 + 
   2*g0*i2*p2*r1*v2*z1 - 4*n0*n2*p2*r1*v2*z1 - 4*n0*n1*p2*r2*v2*z1 + 
   2*f1*j2*k2*s0*v2*z1 - 2*f1*i2*p2*s0*v2*z1 + 2*m2*n1*p2*s0*v2*z1 + 
   2*k2*n2*r1*s0*v2*z1 + 2*k2*n1*r2*s0*v2*z1 + 2*f0*j2*k2*s1*v2*z1 - 
   2*f0*i2*p2*s1*v2*z1 + 2*m2*n0*p2*s1*v2*z1 + 2*k2*n2*r0*s1*v2*z1 + 
   2*k2*n0*r2*s1*v2*z1 - 4*k2*m2*s0*s1*v2*z1 + 2*k2*n1*r0*s2*v2*z1 + 
   2*k2*n0*r1*s2*v2*z1 + 2*g1*pow(j2,2)*v0*v2*z1 - 2*g1*i2*o2*v0*v2*z1 + 
   4*n1*n2*o2*v0*v2*z1 - 4*j2*n2*s1*v0*v2*z1 - 4*j2*n1*s2*v0*v2*z1 + 
   4*i2*s1*s2*v0*v2*z1 + 2*g0*pow(j2,2)*v1*v2*z1 - 2*g0*i2*o2*v1*v2*z1 + 
   4*n0*n2*o2*v1*v2*z1 - 4*j2*n2*s0*v1*v2*z1 - 4*j2*n0*s2*v1*v2*z1 + 
   4*i2*s0*s2*v1*v2*z1 + 2*n0*n1*o2*pow(v2,2)*z1 - 
   2*j2*n1*s0*pow(v2,2)*z1 - 2*j2*n0*s1*pow(v2,2)*z1 + 
   2*i2*s0*s1*pow(v2,2)*z1 - 2*f1*k2*m2*o2*w0*z1 + 2*e2*k2*n1*o2*w0*z1 + 
   2*e1*k2*n2*o2*w0*z1 + 2*f1*j2*m2*p2*w0*z1 - 2*e2*j2*n1*p2*w0*z1 - 
   2*e1*j2*n2*p2*w0*z1 + 2*f2*j2*k2*r1*w0*z1 - 2*f2*i2*p2*r1*w0*z1 + 
   2*m2*n2*p2*r1*w0*z1 + 2*f1*j2*k2*r2*w0*z1 - 2*f1*i2*p2*r2*w0*z1 + 
   2*m2*n1*p2*r2*w0*z1 - 4*k2*n2*r1*r2*w0*z1 - 2*k2*n1*pow(r2,2)*w0*z1 - 
   2*e2*j2*k2*s1*w0*z1 + 2*e2*i2*p2*s1*w0*z1 - 2*pow(m2,2)*p2*s1*w0*z1 + 
   2*k2*m2*r2*s1*w0*z1 - 2*e1*j2*k2*s2*w0*z1 + 2*e1*i2*p2*s2*w0*z1 + 
   2*k2*m2*r1*s2*w0*z1 - 2*f2*pow(j2,2)*v1*w0*z1 + 2*f2*i2*o2*v1*w0*z1 - 
   2*m2*n2*o2*v1*w0*z1 + 2*j2*n2*r2*v1*w0*z1 + 2*j2*m2*s2*v1*w0*z1 - 
   2*i2*r2*s2*v1*w0*z1 - 2*f1*pow(j2,2)*v2*w0*z1 + 2*f1*i2*o2*v2*w0*z1 - 
   2*m2*n1*o2*v2*w0*z1 + 2*j2*n2*r1*v2*w0*z1 + 2*j2*n1*r2*v2*w0*z1 + 
   2*j2*m2*s1*v2*w0*z1 - 2*i2*r2*s1*v2*w0*z1 - 2*i2*r1*s2*v2*w0*z1 - 
   2*f0*k2*m2*o2*w1*z1 + 2*e2*k2*n0*o2*w1*z1 + 2*e0*k2*n2*o2*w1*z1 + 
   2*f0*j2*m2*p2*w1*z1 - 2*e2*j2*n0*p2*w1*z1 - 2*e0*j2*n2*p2*w1*z1 + 
   2*f2*j2*k2*r0*w1*z1 - 2*f2*i2*p2*r0*w1*z1 + 2*m2*n2*p2*r0*w1*z1 + 
   2*f0*j2*k2*r2*w1*z1 - 2*f0*i2*p2*r2*w1*z1 + 2*m2*n0*p2*r2*w1*z1 - 
   4*k2*n2*r0*r2*w1*z1 - 2*k2*n0*pow(r2,2)*w1*z1 - 2*e2*j2*k2*s0*w1*z1 + 
   2*e2*i2*p2*s0*w1*z1 - 2*pow(m2,2)*p2*s0*w1*z1 + 2*k2*m2*r2*s0*w1*z1 - 
   2*e0*j2*k2*s2*w1*z1 + 2*e0*i2*p2*s2*w1*z1 + 2*k2*m2*r0*s2*w1*z1 - 
   2*f2*pow(j2,2)*v0*w1*z1 + 2*f2*i2*o2*v0*w1*z1 - 2*m2*n2*o2*v0*w1*z1 + 
   2*j2*n2*r2*v0*w1*z1 + 2*j2*m2*s2*v0*w1*z1 - 2*i2*r2*s2*v0*w1*z1 - 
   2*f0*pow(j2,2)*v2*w1*z1 + 2*f0*i2*o2*v2*w1*z1 - 2*m2*n0*o2*v2*w1*z1 + 
   2*j2*n2*r0*v2*w1*z1 + 2*j2*n0*r2*v2*w1*z1 + 2*j2*m2*s0*v2*w1*z1 - 
   2*i2*r2*s0*v2*w1*z1 - 2*i2*r0*s2*v2*w1*z1 + 2*e2*pow(j2,2)*w0*w1*z1 - 
   2*e2*i2*o2*w0*w1*z1 + 2*pow(m2,2)*o2*w0*w1*z1 - 4*j2*m2*r2*w0*w1*z1 + 
   2*i2*pow(r2,2)*w0*w1*z1 + 2*e1*k2*n0*o2*w2*z1 + 2*e0*k2*n1*o2*w2*z1 - 
   2*e1*j2*n0*p2*w2*z1 - 2*e0*j2*n1*p2*w2*z1 + 2*f1*j2*k2*r0*w2*z1 - 
   2*f1*i2*p2*r0*w2*z1 + 2*m2*n1*p2*r0*w2*z1 + 2*f0*j2*k2*r1*w2*z1 - 
   2*f0*i2*p2*r1*w2*z1 + 2*m2*n0*p2*r1*w2*z1 - 4*k2*n2*r0*r1*w2*z1 - 
   4*k2*n1*r0*r2*w2*z1 - 4*k2*n0*r1*r2*w2*z1 - 2*e1*j2*k2*s0*w2*z1 + 
   2*e1*i2*p2*s0*w2*z1 + 2*k2*m2*r1*s0*w2*z1 - 2*e0*j2*k2*s1*w2*z1 + 
   2*e0*i2*p2*s1*w2*z1 + 2*k2*m2*r0*s1*w2*z1 - 2*f1*pow(j2,2)*v0*w2*z1 + 
   2*f1*i2*o2*v0*w2*z1 - 2*m2*n1*o2*v0*w2*z1 + 2*j2*n2*r1*v0*w2*z1 + 
   2*j2*n1*r2*v0*w2*z1 + 2*j2*m2*s1*v0*w2*z1 - 2*i2*r2*s1*v0*w2*z1 - 
   2*i2*r1*s2*v0*w2*z1 - 2*f0*pow(j2,2)*v1*w2*z1 + 2*f0*i2*o2*v1*w2*z1 - 
   2*m2*n0*o2*v1*w2*z1 + 2*j2*n2*r0*v1*w2*z1 + 2*j2*n0*r2*v1*w2*z1 + 
   2*j2*m2*s0*v1*w2*z1 - 2*i2*r2*s0*v1*w2*z1 - 2*i2*r0*s2*v1*w2*z1 + 
   2*j2*n1*r0*v2*w2*z1 + 2*j2*n0*r1*v2*w2*z1 - 2*i2*r1*s0*v2*w2*z1 - 
   2*i2*r0*s1*v2*w2*z1 + 2*e1*pow(j2,2)*w0*w2*z1 - 2*e1*i2*o2*w0*w2*z1 - 
   4*j2*m2*r1*w0*w2*z1 + 4*i2*r1*r2*w0*w2*z1 + 2*e0*pow(j2,2)*w1*w2*z1 - 
   2*e0*i2*o2*w1*w2*z1 - 4*j2*m2*r0*w1*w2*z1 + 4*i2*r0*r2*w1*w2*z1 + 
   2*i2*r0*r1*pow(w2,2)*z1 + 2*e1*n0*n1*pow(p2,2)*z2 + 
   e0*pow(n1,2)*pow(p2,2)*z2 + 2*f1*k2*n1*p2*r0*z2 + 
   2*f1*k2*n0*p2*r1*z2 + 2*f0*k2*n1*p2*r1*z2 + 2*g1*pow(k2,2)*r0*r1*z2 + 
   g0*pow(k2,2)*pow(r1,2)*z2 - 2*e1*k2*n1*p2*s0*z2 - 
   2*f1*pow(k2,2)*r1*s0*z2 - 2*e1*k2*n0*p2*s1*z2 - 2*e0*k2*n1*p2*s1*z2 - 
   2*f1*pow(k2,2)*r0*s1*z2 - 2*f0*pow(k2,2)*r1*s1*z2 + 
   2*e1*pow(k2,2)*s0*s1*z2 + e0*pow(k2,2)*pow(s1,2)*z2 - 
   2*e1*n0*n1*o2*t2*z2 - e0*pow(n1,2)*o2*t2*z2 - 2*f1*j2*n1*r0*t2*z2 - 
   2*f1*j2*n0*r1*t2*z2 - 2*f0*j2*n1*r1*t2*z2 - 2*g1*i2*r0*r1*t2*z2 + 
   4*n1*n2*r0*r1*t2*z2 - g0*i2*pow(r1,2)*t2*z2 + 
   2*n0*n2*pow(r1,2)*t2*z2 + 2*pow(n1,2)*r0*r2*t2*z2 + 
   4*n0*n1*r1*r2*t2*z2 + 2*e1*j2*n1*s0*t2*z2 + 2*f1*i2*r1*s0*t2*z2 - 
   2*m2*n1*r1*s0*t2*z2 + 2*e1*j2*n0*s1*t2*z2 + 2*e0*j2*n1*s1*t2*z2 + 
   2*f1*i2*r0*s1*t2*z2 - 2*m2*n1*r0*s1*t2*z2 + 2*f0*i2*r1*s1*t2*z2 - 
   2*m2*n0*r1*s1*t2*z2 - 2*e1*i2*s0*s1*t2*z2 - e0*i2*pow(s1,2)*t2*z2 - 
   2*f1*k2*n1*o2*v0*z2 + 2*f1*j2*n1*p2*v0*z2 - 2*g1*j2*k2*r1*v0*z2 + 
   2*g1*i2*p2*r1*v0*z2 - 4*n1*n2*p2*r1*v0*z2 - 2*pow(n1,2)*p2*r2*v0*z2 + 
   2*f1*j2*k2*s1*v0*z2 - 2*f1*i2*p2*s1*v0*z2 + 2*m2*n1*p2*s1*v0*z2 + 
   2*k2*n2*r1*s1*v0*z2 + 2*k2*n1*r2*s1*v0*z2 - 2*k2*m2*pow(s1,2)*v0*z2 + 
   2*k2*n1*r1*s2*v0*z2 - 2*f1*k2*n0*o2*v1*z2 - 2*f0*k2*n1*o2*v1*z2 + 
   2*f1*j2*n0*p2*v1*z2 + 2*f0*j2*n1*p2*v1*z2 - 2*g1*j2*k2*r0*v1*z2 + 
   2*g1*i2*p2*r0*v1*z2 - 4*n1*n2*p2*r0*v1*z2 - 2*g0*j2*k2*r1*v1*z2 + 
   2*g0*i2*p2*r1*v1*z2 - 4*n0*n2*p2*r1*v1*z2 - 4*n0*n1*p2*r2*v1*z2 + 
   2*f1*j2*k2*s0*v1*z2 - 2*f1*i2*p2*s0*v1*z2 + 2*m2*n1*p2*s0*v1*z2 + 
   2*k2*n2*r1*s0*v1*z2 + 2*k2*n1*r2*s0*v1*z2 + 2*f0*j2*k2*s1*v1*z2 - 
   2*f0*i2*p2*s1*v1*z2 + 2*m2*n0*p2*s1*v1*z2 + 2*k2*n2*r0*s1*v1*z2 + 
   2*k2*n0*r2*s1*v1*z2 - 4*k2*m2*s0*s1*v1*z2 + 2*k2*n1*r0*s2*v1*z2 + 
   2*k2*n0*r1*s2*v1*z2 + 2*g1*pow(j2,2)*v0*v1*z2 - 2*g1*i2*o2*v0*v1*z2 + 
   4*n1*n2*o2*v0*v1*z2 - 4*j2*n2*s1*v0*v1*z2 - 4*j2*n1*s2*v0*v1*z2 + 
   4*i2*s1*s2*v0*v1*z2 + g0*pow(j2,2)*pow(v1,2)*z2 - 
   g0*i2*o2*pow(v1,2)*z2 + 2*n0*n2*o2*pow(v1,2)*z2 - 
   2*j2*n2*s0*pow(v1,2)*z2 - 2*j2*n0*s2*pow(v1,2)*z2 + 
   2*i2*s0*s2*pow(v1,2)*z2 - 2*pow(n1,2)*p2*r0*v2*z2 - 
   4*n0*n1*p2*r1*v2*z2 + 2*k2*n1*r1*s0*v2*z2 + 2*k2*n1*r0*s1*v2*z2 + 
   2*k2*n0*r1*s1*v2*z2 + 2*pow(n1,2)*o2*v0*v2*z2 - 4*j2*n1*s1*v0*v2*z2 + 
   2*i2*pow(s1,2)*v0*v2*z2 + 4*n0*n1*o2*v1*v2*z2 - 4*j2*n1*s0*v1*v2*z2 - 
   4*j2*n0*s1*v1*v2*z2 + 4*i2*s0*s1*v1*v2*z2 + 2*e1*k2*n1*o2*w0*z2 - 
   2*e1*j2*n1*p2*w0*z2 + 2*f1*j2*k2*r1*w0*z2 - 2*f1*i2*p2*r1*w0*z2 + 
   2*m2*n1*p2*r1*w0*z2 - 2*k2*n2*pow(r1,2)*w0*z2 - 4*k2*n1*r1*r2*w0*z2 - 
   2*e1*j2*k2*s1*w0*z2 + 2*e1*i2*p2*s1*w0*z2 + 2*k2*m2*r1*s1*w0*z2 - 
   2*f1*pow(j2,2)*v1*w0*z2 + 2*f1*i2*o2*v1*w0*z2 - 2*m2*n1*o2*v1*w0*z2 + 
   2*j2*n2*r1*v1*w0*z2 + 2*j2*n1*r2*v1*w0*z2 + 2*j2*m2*s1*v1*w0*z2 - 
   2*i2*r2*s1*v1*w0*z2 - 2*i2*r1*s2*v1*w0*z2 + 2*j2*n1*r1*v2*w0*z2 - 
   2*i2*r1*s1*v2*w0*z2 + 2*e1*k2*n0*o2*w1*z2 + 2*e0*k2*n1*o2*w1*z2 - 
   2*e1*j2*n0*p2*w1*z2 - 2*e0*j2*n1*p2*w1*z2 + 2*f1*j2*k2*r0*w1*z2 - 
   2*f1*i2*p2*r0*w1*z2 + 2*m2*n1*p2*r0*w1*z2 + 2*f0*j2*k2*r1*w1*z2 - 
   2*f0*i2*p2*r1*w1*z2 + 2*m2*n0*p2*r1*w1*z2 - 4*k2*n2*r0*r1*w1*z2 - 
   4*k2*n1*r0*r2*w1*z2 - 4*k2*n0*r1*r2*w1*z2 - 2*e1*j2*k2*s0*w1*z2 + 
   2*e1*i2*p2*s0*w1*z2 + 2*k2*m2*r1*s0*w1*z2 - 2*e0*j2*k2*s1*w1*z2 + 
   2*e0*i2*p2*s1*w1*z2 + 2*k2*m2*r0*s1*w1*z2 - 2*f1*pow(j2,2)*v0*w1*z2 + 
   2*f1*i2*o2*v0*w1*z2 - 2*m2*n1*o2*v0*w1*z2 + 2*j2*n2*r1*v0*w1*z2 + 
   2*j2*n1*r2*v0*w1*z2 + 2*j2*m2*s1*v0*w1*z2 - 2*i2*r2*s1*v0*w1*z2 - 
   2*i2*r1*s2*v0*w1*z2 - 2*f0*pow(j2,2)*v1*w1*z2 + 2*f0*i2*o2*v1*w1*z2 - 
   2*m2*n0*o2*v1*w1*z2 + 2*j2*n2*r0*v1*w1*z2 + 2*j2*n0*r2*v1*w1*z2 + 
   2*j2*m2*s0*v1*w1*z2 - 2*i2*r2*s0*v1*w1*z2 - 2*i2*r0*s2*v1*w1*z2 + 
   2*j2*n1*r0*v2*w1*z2 + 2*j2*n0*r1*v2*w1*z2 - 2*i2*r1*s0*v2*w1*z2 - 
   2*i2*r0*s1*v2*w1*z2 + 2*e1*pow(j2,2)*w0*w1*z2 - 2*e1*i2*o2*w0*w1*z2 - 
   4*j2*m2*r1*w0*w1*z2 + 4*i2*r1*r2*w0*w1*z2 + 
   e0*pow(j2,2)*pow(w1,2)*z2 - e0*i2*o2*pow(w1,2)*z2 - 
   2*j2*m2*r0*pow(w1,2)*z2 + 2*i2*r0*r2*pow(w1,2)*z2 - 
   4*k2*n1*r0*r1*w2*z2 - 2*k2*n0*pow(r1,2)*w2*z2 + 2*j2*n1*r1*v0*w2*z2 - 
   2*i2*r1*s1*v0*w2*z2 + 2*j2*n1*r0*v1*w2*z2 + 2*j2*n0*r1*v1*w2*z2 - 
   2*i2*r1*s0*v1*w2*z2 - 2*i2*r0*s1*v1*w2*z2 + 2*i2*pow(r1,2)*w0*w2*z2 + 
   4*i2*r0*r1*w1*w2*z2
;

	k[9]  = pow(d1,2)*e1*pow(k2,2)*o2 - 2*c1*d1*f1*pow(k2,2)*o2 + 
   pow(c1,2)*g1*pow(k2,2)*o2 - 2*pow(d1,2)*e1*j2*k2*p2 + 
   4*c1*d1*f1*j2*k2*p2 - 2*pow(c1,2)*g1*j2*k2*p2 + 
   pow(d1,2)*e1*i2*pow(p2,2) - 2*c1*d1*f1*i2*pow(p2,2) + 
   pow(c1,2)*g1*i2*pow(p2,2) - 2*d1*e1*l2*n1*pow(p2,2) + 
   2*c1*f1*l2*n1*pow(p2,2) + 2*c1*d1*m2*n1*pow(p2,2) - 
   2*c1*c2*pow(n1,2)*pow(p2,2) - 2*pow(c1,2)*n1*n2*pow(p2,2) + 
   2*d1*e1*k2*n1*p2*q2 - 2*c1*f1*k2*n1*p2*q2 - 2*d1*f1*k2*l2*p2*r1 + 
   2*c1*g1*k2*l2*p2*r1 + 2*pow(d1,2)*k2*m2*p2*r1 - 2*c2*d1*k2*n1*p2*r1 - 
   2*c1*d2*k2*n1*p2*r1 - 2*c1*d1*k2*n2*p2*r1 + 2*d1*f1*pow(k2,2)*q2*r1 - 
   2*c1*g1*pow(k2,2)*q2*r1 - 2*d1*d2*pow(k2,2)*pow(r1,2) - 
   2*c1*d1*k2*n1*p2*r2 - 2*pow(d1,2)*pow(k2,2)*r1*r2 + 
   2*d1*e1*k2*l2*p2*s1 - 2*c1*f1*k2*l2*p2*s1 - 2*c1*d1*k2*m2*p2*s1 + 
   4*c1*c2*k2*n1*p2*s1 + 2*pow(c1,2)*k2*n2*p2*s1 - 
   2*d1*e1*pow(k2,2)*q2*s1 + 2*c1*f1*pow(k2,2)*q2*s1 + 
   2*c2*d1*pow(k2,2)*r1*s1 + 2*c1*d2*pow(k2,2)*r1*s1 + 
   2*c1*d1*pow(k2,2)*r2*s1 - 2*c1*c2*pow(k2,2)*pow(s1,2) + 
   2*pow(c1,2)*k2*n1*p2*s2 + 2*c1*d1*pow(k2,2)*r1*s2 - 
   2*pow(c1,2)*pow(k2,2)*s1*s2 + pow(d1,2)*e1*pow(j2,2)*t2 - 
   2*c1*d1*f1*pow(j2,2)*t2 + pow(c1,2)*g1*pow(j2,2)*t2 - 
   pow(d1,2)*e1*i2*o2*t2 + 2*c1*d1*f1*i2*o2*t2 - pow(c1,2)*g1*i2*o2*t2 + 
   2*d1*e1*l2*n1*o2*t2 - 2*c1*f1*l2*n1*o2*t2 - 2*c1*d1*m2*n1*o2*t2 + 
   2*c1*c2*pow(n1,2)*o2*t2 + 2*pow(c1,2)*n1*n2*o2*t2 - 
   2*d1*e1*j2*n1*q2*t2 + 2*c1*f1*j2*n1*q2*t2 + 
   e1*pow(n1,2)*pow(q2,2)*t2 + 2*d1*f1*j2*l2*r1*t2 - 
   2*c1*g1*j2*l2*r1*t2 - 2*pow(d1,2)*j2*m2*r1*t2 + 2*c2*d1*j2*n1*r1*t2 + 
   2*c1*d2*j2*n1*r1*t2 + 2*c1*d1*j2*n2*r1*t2 - 2*d1*f1*i2*q2*r1*t2 + 
   2*c1*g1*i2*q2*r1*t2 + 2*f1*l2*n1*q2*r1*t2 + 2*d1*m2*n1*q2*r1*t2 - 
   2*c2*pow(n1,2)*q2*r1*t2 - 4*c1*n1*n2*q2*r1*t2 + 
   2*d1*d2*i2*pow(r1,2)*t2 + g1*pow(l2,2)*pow(r1,2)*t2 - 
   2*d2*l2*n1*pow(r1,2)*t2 - 2*d1*l2*n2*pow(r1,2)*t2 + 
   2*c1*d1*j2*n1*r2*t2 - 2*c1*pow(n1,2)*q2*r2*t2 + 
   2*pow(d1,2)*i2*r1*r2*t2 - 4*d1*l2*n1*r1*r2*t2 - 2*d1*e1*j2*l2*s1*t2 + 
   2*c1*f1*j2*l2*s1*t2 + 2*c1*d1*j2*m2*s1*t2 - 4*c1*c2*j2*n1*s1*t2 - 
   2*pow(c1,2)*j2*n2*s1*t2 + 2*d1*e1*i2*q2*s1*t2 - 2*c1*f1*i2*q2*s1*t2 - 
   2*e1*l2*n1*q2*s1*t2 + 2*c1*m2*n1*q2*s1*t2 - 2*c2*d1*i2*r1*s1*t2 - 
   2*c1*d2*i2*r1*s1*t2 - 2*f1*pow(l2,2)*r1*s1*t2 + 2*d1*l2*m2*r1*s1*t2 + 
   2*c2*l2*n1*r1*s1*t2 + 2*c1*l2*n2*r1*s1*t2 - 2*c1*d1*i2*r2*s1*t2 + 
   2*c1*l2*n1*r2*s1*t2 + 2*c1*c2*i2*pow(s1,2)*t2 + 
   e1*pow(l2,2)*pow(s1,2)*t2 - 2*c1*l2*m2*pow(s1,2)*t2 - 
   2*pow(c1,2)*j2*n1*s2*t2 - 2*c1*d1*i2*r1*s2*t2 + 2*c1*l2*n1*r1*s2*t2 + 
   2*pow(c1,2)*i2*s1*s2*t2 - 2*pow(f1,2)*k2*l2*o2*u1 + 
   2*e1*g1*k2*l2*o2*u1 + 2*d1*f1*k2*m2*o2*u1 - 2*c1*g1*k2*m2*o2*u1 - 
   2*d2*e1*k2*n1*o2*u1 - 2*d1*e2*k2*n1*o2*u1 + 2*c2*f1*k2*n1*o2*u1 + 
   2*c1*f2*k2*n1*o2*u1 - 2*d1*e1*k2*n2*o2*u1 + 2*c1*f1*k2*n2*o2*u1 + 
   2*pow(f1,2)*j2*l2*p2*u1 - 2*e1*g1*j2*l2*p2*u1 - 2*d1*f1*j2*m2*p2*u1 + 
   2*c1*g1*j2*m2*p2*u1 + 2*d2*e1*j2*n1*p2*u1 + 2*d1*e2*j2*n1*p2*u1 - 
   2*c2*f1*j2*n1*p2*u1 - 2*c1*f2*j2*n1*p2*u1 + 2*d1*e1*j2*n2*p2*u1 - 
   2*c1*f1*j2*n2*p2*u1 + 2*pow(f1,2)*j2*k2*q2*u1 - 2*e1*g1*j2*k2*q2*u1 - 
   2*pow(f1,2)*i2*p2*q2*u1 + 2*e1*g1*i2*p2*q2*u1 + 4*f1*m2*n1*p2*q2*u1 - 
   2*e2*pow(n1,2)*p2*q2*u1 - 4*e1*n1*n2*p2*q2*u1 - 2*d2*f1*j2*k2*r1*u1 - 
   2*d1*f2*j2*k2*r1*u1 + 2*c2*g1*j2*k2*r1*u1 + 2*c1*g2*j2*k2*r1*u1 + 
   2*d2*f1*i2*p2*r1*u1 + 2*d1*f2*i2*p2*r1*u1 - 2*c2*g1*i2*p2*r1*u1 - 
   2*c1*g2*i2*p2*r1*u1 + 2*g1*l2*m2*p2*r1*u1 - 2*f2*l2*n1*p2*r1*u1 - 
   2*d2*m2*n1*p2*r1*u1 - 2*f1*l2*n2*p2*r1*u1 - 2*d1*m2*n2*p2*r1*u1 + 
   4*c2*n1*n2*p2*r1*u1 + 2*c1*pow(n2,2)*p2*r1*u1 + 2*g1*k2*m2*q2*r1*u1 - 
   2*f2*k2*n1*q2*r1*u1 - 2*f1*k2*n2*q2*r1*u1 - 2*g2*k2*l2*pow(r1,2)*u1 + 
   2*d2*k2*n2*pow(r1,2)*u1 - 2*d1*f1*j2*k2*r2*u1 + 2*c1*g1*j2*k2*r2*u1 + 
   2*d1*f1*i2*p2*r2*u1 - 2*c1*g1*i2*p2*r2*u1 - 2*f1*l2*n1*p2*r2*u1 - 
   2*d1*m2*n1*p2*r2*u1 + 2*c2*pow(n1,2)*p2*r2*u1 + 4*c1*n1*n2*p2*r2*u1 - 
   2*f1*k2*n1*q2*r2*u1 - 4*g1*k2*l2*r1*r2*u1 + 4*d2*k2*n1*r1*r2*u1 + 
   4*d1*k2*n2*r1*r2*u1 + 2*d1*k2*n1*pow(r2,2)*u1 + 2*d2*e1*j2*k2*s1*u1 + 
   2*d1*e2*j2*k2*s1*u1 - 2*c2*f1*j2*k2*s1*u1 - 2*c1*f2*j2*k2*s1*u1 - 
   2*d2*e1*i2*p2*s1*u1 - 2*d1*e2*i2*p2*s1*u1 + 2*c2*f1*i2*p2*s1*u1 + 
   2*c1*f2*i2*p2*s1*u1 - 2*f1*l2*m2*p2*s1*u1 + 2*d1*pow(m2,2)*p2*s1*u1 + 
   2*e2*l2*n1*p2*s1*u1 - 2*c2*m2*n1*p2*s1*u1 + 2*e1*l2*n2*p2*s1*u1 - 
   2*c1*m2*n2*p2*s1*u1 - 2*f1*k2*m2*q2*s1*u1 + 2*e2*k2*n1*q2*s1*u1 + 
   2*e1*k2*n2*q2*s1*u1 + 4*f2*k2*l2*r1*s1*u1 - 2*d2*k2*m2*r1*s1*u1 - 
   2*c2*k2*n2*r1*s1*u1 + 4*f1*k2*l2*r2*s1*u1 - 2*d1*k2*m2*r2*s1*u1 - 
   2*c2*k2*n1*r2*s1*u1 - 2*c1*k2*n2*r2*s1*u1 - 2*e2*k2*l2*pow(s1,2)*u1 + 
   2*c2*k2*m2*pow(s1,2)*u1 + 2*d1*e1*j2*k2*s2*u1 - 2*c1*f1*j2*k2*s2*u1 - 
   2*d1*e1*i2*p2*s2*u1 + 2*c1*f1*i2*p2*s2*u1 + 2*e1*l2*n1*p2*s2*u1 - 
   2*c1*m2*n1*p2*s2*u1 + 2*e1*k2*n1*q2*s2*u1 + 4*f1*k2*l2*r1*s2*u1 - 
   2*d1*k2*m2*r1*s2*u1 - 2*c2*k2*n1*r1*s2*u1 - 2*c1*k2*n2*r1*s2*u1 - 
   2*c1*k2*n1*r2*s2*u1 - 4*e1*k2*l2*s1*s2*u1 + 4*c1*k2*m2*s1*s2*u1 - 
   2*f1*f2*pow(j2,2)*pow(u1,2) + e2*g1*pow(j2,2)*pow(u1,2) + 
   e1*g2*pow(j2,2)*pow(u1,2) + 2*f1*f2*i2*o2*pow(u1,2) - 
   e2*g1*i2*o2*pow(u1,2) - e1*g2*i2*o2*pow(u1,2) + 
   g1*pow(m2,2)*o2*pow(u1,2) - 2*f2*m2*n1*o2*pow(u1,2) - 
   2*f1*m2*n2*o2*pow(u1,2) + 2*e2*n1*n2*o2*pow(u1,2) + 
   e1*pow(n2,2)*o2*pow(u1,2) - 2*g2*j2*m2*r1*pow(u1,2) + 
   2*f2*j2*n2*r1*pow(u1,2) - 2*g1*j2*m2*r2*pow(u1,2) + 
   2*f2*j2*n1*r2*pow(u1,2) + 2*f1*j2*n2*r2*pow(u1,2) + 
   2*g2*i2*r1*r2*pow(u1,2) - 2*pow(n2,2)*r1*r2*pow(u1,2) + 
   g1*i2*pow(r2,2)*pow(u1,2) - 2*n1*n2*pow(r2,2)*pow(u1,2) + 
   2*f2*j2*m2*s1*pow(u1,2) - 2*e2*j2*n2*s1*pow(u1,2) - 
   2*f2*i2*r2*s1*pow(u1,2) + 2*m2*n2*r2*s1*pow(u1,2) + 
   2*f1*j2*m2*s2*pow(u1,2) - 2*e2*j2*n1*s2*pow(u1,2) - 
   2*e1*j2*n2*s2*pow(u1,2) - 2*f2*i2*r1*s2*pow(u1,2) + 
   2*m2*n2*r1*s2*pow(u1,2) - 2*f1*i2*r2*s2*pow(u1,2) + 
   2*m2*n1*r2*s2*pow(u1,2) + 2*e2*i2*s1*s2*pow(u1,2) - 
   2*pow(m2,2)*s1*s2*pow(u1,2) + e1*i2*pow(s2,2)*pow(u1,2) - 
   2*d1*e1*k2*n1*o2*u2 + 2*c1*f1*k2*n1*o2*u2 + 2*d1*e1*j2*n1*p2*u2 - 
   2*c1*f1*j2*n1*p2*u2 - 2*e1*pow(n1,2)*p2*q2*u2 - 2*d1*f1*j2*k2*r1*u2 + 
   2*c1*g1*j2*k2*r1*u2 + 2*d1*f1*i2*p2*r1*u2 - 2*c1*g1*i2*p2*r1*u2 - 
   2*f1*l2*n1*p2*r1*u2 - 2*d1*m2*n1*p2*r1*u2 + 2*c2*pow(n1,2)*p2*r1*u2 + 
   4*c1*n1*n2*p2*r1*u2 - 2*f1*k2*n1*q2*r1*u2 - 2*g1*k2*l2*pow(r1,2)*u2 + 
   2*d2*k2*n1*pow(r1,2)*u2 + 2*d1*k2*n2*pow(r1,2)*u2 + 
   2*c1*pow(n1,2)*p2*r2*u2 + 4*d1*k2*n1*r1*r2*u2 + 2*d1*e1*j2*k2*s1*u2 - 
   2*c1*f1*j2*k2*s1*u2 - 2*d1*e1*i2*p2*s1*u2 + 2*c1*f1*i2*p2*s1*u2 + 
   2*e1*l2*n1*p2*s1*u2 - 2*c1*m2*n1*p2*s1*u2 + 2*e1*k2*n1*q2*s1*u2 + 
   4*f1*k2*l2*r1*s1*u2 - 2*d1*k2*m2*r1*s1*u2 - 2*c2*k2*n1*r1*s1*u2 - 
   2*c1*k2*n2*r1*s1*u2 - 2*c1*k2*n1*r2*s1*u2 - 2*e1*k2*l2*pow(s1,2)*u2 + 
   2*c1*k2*m2*pow(s1,2)*u2 - 2*c1*k2*n1*r1*s2*u2 - 
   2*pow(f1,2)*pow(j2,2)*u1*u2 + 2*e1*g1*pow(j2,2)*u1*u2 + 
   2*pow(f1,2)*i2*o2*u1*u2 - 2*e1*g1*i2*o2*u1*u2 - 4*f1*m2*n1*o2*u1*u2 + 
   2*e2*pow(n1,2)*o2*u1*u2 + 4*e1*n1*n2*o2*u1*u2 - 4*g1*j2*m2*r1*u1*u2 + 
   4*f2*j2*n1*r1*u1*u2 + 4*f1*j2*n2*r1*u1*u2 + 2*g2*i2*pow(r1,2)*u1*u2 - 
   2*pow(n2,2)*pow(r1,2)*u1*u2 + 4*f1*j2*n1*r2*u1*u2 + 
   4*g1*i2*r1*r2*u1*u2 - 8*n1*n2*r1*r2*u1*u2 - 
   2*pow(n1,2)*pow(r2,2)*u1*u2 + 4*f1*j2*m2*s1*u1*u2 - 
   4*e2*j2*n1*s1*u1*u2 - 4*e1*j2*n2*s1*u1*u2 - 4*f2*i2*r1*s1*u1*u2 + 
   4*m2*n2*r1*s1*u1*u2 - 4*f1*i2*r2*s1*u1*u2 + 4*m2*n1*r2*s1*u1*u2 + 
   2*e2*i2*pow(s1,2)*u1*u2 - 2*pow(m2,2)*pow(s1,2)*u1*u2 - 
   4*e1*j2*n1*s2*u1*u2 - 4*f1*i2*r1*s2*u1*u2 + 4*m2*n1*r1*s2*u1*u2 + 
   4*e1*i2*s1*s2*u1*u2 + e1*pow(n1,2)*o2*pow(u2,2) + 
   2*f1*j2*n1*r1*pow(u2,2) + g1*i2*pow(r1,2)*pow(u2,2) - 
   2*n1*n2*pow(r1,2)*pow(u2,2) - 2*pow(n1,2)*r1*r2*pow(u2,2) - 
   2*e1*j2*n1*s1*pow(u2,2) - 2*f1*i2*r1*s1*pow(u2,2) + 
   2*m2*n1*r1*s1*pow(u2,2) + e1*i2*pow(s1,2)*pow(u2,2) + 
   2*d1*f1*k2*l2*o2*v1 - 2*c1*g1*k2*l2*o2*v1 - 2*pow(d1,2)*k2*m2*o2*v1 + 
   2*c2*d1*k2*n1*o2*v1 + 2*c1*d2*k2*n1*o2*v1 + 2*c1*d1*k2*n2*o2*v1 - 
   2*d1*f1*j2*l2*p2*v1 + 2*c1*g1*j2*l2*p2*v1 + 2*pow(d1,2)*j2*m2*p2*v1 - 
   2*c2*d1*j2*n1*p2*v1 - 2*c1*d2*j2*n1*p2*v1 - 2*c1*d1*j2*n2*p2*v1 - 
   2*d1*f1*j2*k2*q2*v1 + 2*c1*g1*j2*k2*q2*v1 + 2*d1*f1*i2*p2*q2*v1 - 
   2*c1*g1*i2*p2*q2*v1 - 2*f1*l2*n1*p2*q2*v1 - 2*d1*m2*n1*p2*q2*v1 + 
   2*c2*pow(n1,2)*p2*q2*v1 + 4*c1*n1*n2*p2*q2*v1 + 
   2*f1*k2*n1*pow(q2,2)*v1 + 4*d1*d2*j2*k2*r1*v1 - 4*d1*d2*i2*p2*r1*v1 - 
   2*g1*pow(l2,2)*p2*r1*v1 + 4*d2*l2*n1*p2*r1*v1 + 4*d1*l2*n2*p2*r1*v1 + 
   2*g1*k2*l2*q2*r1*v1 - 2*d2*k2*n1*q2*r1*v1 - 2*d1*k2*n2*q2*r1*v1 + 
   2*pow(d1,2)*j2*k2*r2*v1 - 2*pow(d1,2)*i2*p2*r2*v1 + 
   4*d1*l2*n1*p2*r2*v1 - 2*d1*k2*n1*q2*r2*v1 - 2*c2*d1*j2*k2*s1*v1 - 
   2*c1*d2*j2*k2*s1*v1 + 2*c2*d1*i2*p2*s1*v1 + 2*c1*d2*i2*p2*s1*v1 + 
   2*f1*pow(l2,2)*p2*s1*v1 - 2*d1*l2*m2*p2*s1*v1 - 2*c2*l2*n1*p2*s1*v1 - 
   2*c1*l2*n2*p2*s1*v1 - 2*f1*k2*l2*q2*s1*v1 + 4*d1*k2*m2*q2*s1*v1 - 
   2*c2*k2*n1*q2*s1*v1 - 2*c1*k2*n2*q2*s1*v1 - 2*d2*k2*l2*r1*s1*v1 - 
   2*d1*k2*l2*r2*s1*v1 + 2*c2*k2*l2*pow(s1,2)*v1 - 2*c1*d1*j2*k2*s2*v1 + 
   2*c1*d1*i2*p2*s2*v1 - 2*c1*l2*n1*p2*s2*v1 - 2*c1*k2*n1*q2*s2*v1 - 
   2*d1*k2*l2*r1*s2*v1 + 4*c1*k2*l2*s1*s2*v1 + 2*d2*f1*pow(j2,2)*u1*v1 + 
   2*d1*f2*pow(j2,2)*u1*v1 - 2*c2*g1*pow(j2,2)*u1*v1 - 
   2*c1*g2*pow(j2,2)*u1*v1 - 2*d2*f1*i2*o2*u1*v1 - 2*d1*f2*i2*o2*u1*v1 + 
   2*c2*g1*i2*o2*u1*v1 + 2*c1*g2*i2*o2*u1*v1 - 2*g1*l2*m2*o2*u1*v1 + 
   2*f2*l2*n1*o2*u1*v1 + 2*d2*m2*n1*o2*u1*v1 + 2*f1*l2*n2*o2*u1*v1 + 
   2*d1*m2*n2*o2*u1*v1 - 4*c2*n1*n2*o2*u1*v1 - 2*c1*pow(n2,2)*o2*u1*v1 + 
   2*g1*j2*m2*q2*u1*v1 - 2*f2*j2*n1*q2*u1*v1 - 2*f1*j2*n2*q2*u1*v1 + 
   2*g2*j2*l2*r1*u1*v1 - 2*d2*j2*n2*r1*u1*v1 - 2*g2*i2*q2*r1*u1*v1 + 
   2*pow(n2,2)*q2*r1*u1*v1 + 2*g1*j2*l2*r2*u1*v1 - 2*d2*j2*n1*r2*u1*v1 - 
   2*d1*j2*n2*r2*u1*v1 - 2*g1*i2*q2*r2*u1*v1 + 4*n1*n2*q2*r2*u1*v1 - 
   2*f2*j2*l2*s1*u1*v1 - 2*d2*j2*m2*s1*u1*v1 + 4*c2*j2*n2*s1*u1*v1 + 
   2*f2*i2*q2*s1*u1*v1 - 2*m2*n2*q2*s1*u1*v1 + 2*d2*i2*r2*s1*u1*v1 - 
   2*l2*n2*r2*s1*u1*v1 - 2*f1*j2*l2*s2*u1*v1 - 2*d1*j2*m2*s2*u1*v1 + 
   4*c2*j2*n1*s2*u1*v1 + 4*c1*j2*n2*s2*u1*v1 + 2*f1*i2*q2*s2*u1*v1 - 
   2*m2*n1*q2*s2*u1*v1 + 2*d2*i2*r1*s2*u1*v1 - 2*l2*n2*r1*s2*u1*v1 + 
   2*d1*i2*r2*s2*u1*v1 - 2*l2*n1*r2*s2*u1*v1 - 4*c2*i2*s1*s2*u1*v1 + 
   4*l2*m2*s1*s2*u1*v1 - 2*c1*i2*pow(s2,2)*u1*v1 + 
   2*d1*f1*pow(j2,2)*u2*v1 - 2*c1*g1*pow(j2,2)*u2*v1 - 
   2*d1*f1*i2*o2*u2*v1 + 2*c1*g1*i2*o2*u2*v1 + 2*f1*l2*n1*o2*u2*v1 + 
   2*d1*m2*n1*o2*u2*v1 - 2*c2*pow(n1,2)*o2*u2*v1 - 4*c1*n1*n2*o2*u2*v1 - 
   2*f1*j2*n1*q2*u2*v1 + 2*g1*j2*l2*r1*u2*v1 - 2*d2*j2*n1*r1*u2*v1 - 
   2*d1*j2*n2*r1*u2*v1 - 2*g1*i2*q2*r1*u2*v1 + 4*n1*n2*q2*r1*u2*v1 - 
   2*d1*j2*n1*r2*u2*v1 + 2*pow(n1,2)*q2*r2*u2*v1 - 2*f1*j2*l2*s1*u2*v1 - 
   2*d1*j2*m2*s1*u2*v1 + 4*c2*j2*n1*s1*u2*v1 + 4*c1*j2*n2*s1*u2*v1 + 
   2*f1*i2*q2*s1*u2*v1 - 2*m2*n1*q2*s1*u2*v1 + 2*d2*i2*r1*s1*u2*v1 - 
   2*l2*n2*r1*s1*u2*v1 + 2*d1*i2*r2*s1*u2*v1 - 2*l2*n1*r2*s1*u2*v1 - 
   2*c2*i2*pow(s1,2)*u2*v1 + 2*l2*m2*pow(s1,2)*u2*v1 + 
   4*c1*j2*n1*s2*u2*v1 + 2*d1*i2*r1*s2*u2*v1 - 2*l2*n1*r1*s2*u2*v1 - 
   4*c1*i2*s1*s2*u2*v1 - 2*d1*d2*pow(j2,2)*pow(v1,2) + 
   2*d1*d2*i2*o2*pow(v1,2) + g1*pow(l2,2)*o2*pow(v1,2) - 
   2*d2*l2*n1*o2*pow(v1,2) - 2*d1*l2*n2*o2*pow(v1,2) - 
   2*g1*j2*l2*q2*pow(v1,2) + 2*d2*j2*n1*q2*pow(v1,2) + 
   2*d1*j2*n2*q2*pow(v1,2) + g1*i2*pow(q2,2)*pow(v1,2) - 
   2*n1*n2*pow(q2,2)*pow(v1,2) + 2*d2*j2*l2*s1*pow(v1,2) - 
   2*d2*i2*q2*s1*pow(v1,2) + 2*l2*n2*q2*s1*pow(v1,2) + 
   2*d1*j2*l2*s2*pow(v1,2) - 2*d1*i2*q2*s2*pow(v1,2) + 
   2*l2*n1*q2*s2*pow(v1,2) - 2*pow(l2,2)*s1*s2*pow(v1,2) + 
   2*c1*d1*k2*n1*o2*v2 - 2*c1*d1*j2*n1*p2*v2 + 2*c1*pow(n1,2)*p2*q2*v2 + 
   2*pow(d1,2)*j2*k2*r1*v2 - 2*pow(d1,2)*i2*p2*r1*v2 + 
   4*d1*l2*n1*p2*r1*v2 - 2*d1*k2*n1*q2*r1*v2 - 2*c1*d1*j2*k2*s1*v2 + 
   2*c1*d1*i2*p2*s1*v2 - 2*c1*l2*n1*p2*s1*v2 - 2*c1*k2*n1*q2*s1*v2 - 
   2*d1*k2*l2*r1*s1*v2 + 2*c1*k2*l2*pow(s1,2)*v2 + 
   2*d1*f1*pow(j2,2)*u1*v2 - 2*c1*g1*pow(j2,2)*u1*v2 - 
   2*d1*f1*i2*o2*u1*v2 + 2*c1*g1*i2*o2*u1*v2 + 2*f1*l2*n1*o2*u1*v2 + 
   2*d1*m2*n1*o2*u1*v2 - 2*c2*pow(n1,2)*o2*u1*v2 - 4*c1*n1*n2*o2*u1*v2 - 
   2*f1*j2*n1*q2*u1*v2 + 2*g1*j2*l2*r1*u1*v2 - 2*d2*j2*n1*r1*u1*v2 - 
   2*d1*j2*n2*r1*u1*v2 - 2*g1*i2*q2*r1*u1*v2 + 4*n1*n2*q2*r1*u1*v2 - 
   2*d1*j2*n1*r2*u1*v2 + 2*pow(n1,2)*q2*r2*u1*v2 - 2*f1*j2*l2*s1*u1*v2 - 
   2*d1*j2*m2*s1*u1*v2 + 4*c2*j2*n1*s1*u1*v2 + 4*c1*j2*n2*s1*u1*v2 + 
   2*f1*i2*q2*s1*u1*v2 - 2*m2*n1*q2*s1*u1*v2 + 2*d2*i2*r1*s1*u1*v2 - 
   2*l2*n2*r1*s1*u1*v2 + 2*d1*i2*r2*s1*u1*v2 - 2*l2*n1*r2*s1*u1*v2 - 
   2*c2*i2*pow(s1,2)*u1*v2 + 2*l2*m2*pow(s1,2)*u1*v2 + 
   4*c1*j2*n1*s2*u1*v2 + 2*d1*i2*r1*s2*u1*v2 - 2*l2*n1*r1*s2*u1*v2 - 
   4*c1*i2*s1*s2*u1*v2 - 2*c1*pow(n1,2)*o2*u2*v2 - 2*d1*j2*n1*r1*u2*v2 + 
   2*pow(n1,2)*q2*r1*u2*v2 + 4*c1*j2*n1*s1*u2*v2 + 2*d1*i2*r1*s1*u2*v2 - 
   2*l2*n1*r1*s1*u2*v2 - 2*c1*i2*pow(s1,2)*u2*v2 - 
   2*pow(d1,2)*pow(j2,2)*v1*v2 + 2*pow(d1,2)*i2*o2*v1*v2 - 
   4*d1*l2*n1*o2*v1*v2 + 4*d1*j2*n1*q2*v1*v2 - 
   2*pow(n1,2)*pow(q2,2)*v1*v2 + 4*d1*j2*l2*s1*v1*v2 - 
   4*d1*i2*q2*s1*v1*v2 + 4*l2*n1*q2*s1*v1*v2 - 
   2*pow(l2,2)*pow(s1,2)*v1*v2 - 2*d1*e1*k2*l2*o2*w1 + 
   2*c1*f1*k2*l2*o2*w1 + 2*c1*d1*k2*m2*o2*w1 - 4*c1*c2*k2*n1*o2*w1 - 
   2*pow(c1,2)*k2*n2*o2*w1 + 2*d1*e1*j2*l2*p2*w1 - 2*c1*f1*j2*l2*p2*w1 - 
   2*c1*d1*j2*m2*p2*w1 + 4*c1*c2*j2*n1*p2*w1 + 2*pow(c1,2)*j2*n2*p2*w1 + 
   2*d1*e1*j2*k2*q2*w1 - 2*c1*f1*j2*k2*q2*w1 - 2*d1*e1*i2*p2*q2*w1 + 
   2*c1*f1*i2*p2*q2*w1 + 2*e1*l2*n1*p2*q2*w1 - 2*c1*m2*n1*p2*q2*w1 - 
   2*e1*k2*n1*pow(q2,2)*w1 - 2*c2*d1*j2*k2*r1*w1 - 2*c1*d2*j2*k2*r1*w1 + 
   2*c2*d1*i2*p2*r1*w1 + 2*c1*d2*i2*p2*r1*w1 + 2*f1*pow(l2,2)*p2*r1*w1 - 
   2*d1*l2*m2*p2*r1*w1 - 2*c2*l2*n1*p2*r1*w1 - 2*c1*l2*n2*p2*r1*w1 - 
   2*f1*k2*l2*q2*r1*w1 - 2*d1*k2*m2*q2*r1*w1 + 4*c2*k2*n1*q2*r1*w1 + 
   4*c1*k2*n2*q2*r1*w1 + 2*d2*k2*l2*pow(r1,2)*w1 - 2*c1*d1*j2*k2*r2*w1 + 
   2*c1*d1*i2*p2*r2*w1 - 2*c1*l2*n1*p2*r2*w1 + 4*c1*k2*n1*q2*r2*w1 + 
   4*d1*k2*l2*r1*r2*w1 + 4*c1*c2*j2*k2*s1*w1 - 4*c1*c2*i2*p2*s1*w1 - 
   2*e1*pow(l2,2)*p2*s1*w1 + 4*c1*l2*m2*p2*s1*w1 + 2*e1*k2*l2*q2*s1*w1 - 
   2*c1*k2*m2*q2*s1*w1 - 2*c2*k2*l2*r1*s1*w1 - 2*c1*k2*l2*r2*s1*w1 + 
   2*pow(c1,2)*j2*k2*s2*w1 - 2*pow(c1,2)*i2*p2*s2*w1 - 
   2*c1*k2*l2*r1*s2*w1 - 2*d2*e1*pow(j2,2)*u1*w1 - 
   2*d1*e2*pow(j2,2)*u1*w1 + 2*c2*f1*pow(j2,2)*u1*w1 + 
   2*c1*f2*pow(j2,2)*u1*w1 + 2*d2*e1*i2*o2*u1*w1 + 2*d1*e2*i2*o2*u1*w1 - 
   2*c2*f1*i2*o2*u1*w1 - 2*c1*f2*i2*o2*u1*w1 + 2*f1*l2*m2*o2*u1*w1 - 
   2*d1*pow(m2,2)*o2*u1*w1 - 2*e2*l2*n1*o2*u1*w1 + 2*c2*m2*n1*o2*u1*w1 - 
   2*e1*l2*n2*o2*u1*w1 + 2*c1*m2*n2*o2*u1*w1 - 2*f1*j2*m2*q2*u1*w1 + 
   2*e2*j2*n1*q2*u1*w1 + 2*e1*j2*n2*q2*u1*w1 - 2*f2*j2*l2*r1*u1*w1 + 
   4*d2*j2*m2*r1*u1*w1 - 2*c2*j2*n2*r1*u1*w1 + 2*f2*i2*q2*r1*u1*w1 - 
   2*m2*n2*q2*r1*u1*w1 - 2*f1*j2*l2*r2*u1*w1 + 4*d1*j2*m2*r2*u1*w1 - 
   2*c2*j2*n1*r2*u1*w1 - 2*c1*j2*n2*r2*u1*w1 + 2*f1*i2*q2*r2*u1*w1 - 
   2*m2*n1*q2*r2*u1*w1 - 4*d2*i2*r1*r2*u1*w1 + 4*l2*n2*r1*r2*u1*w1 - 
   2*d1*i2*pow(r2,2)*u1*w1 + 2*l2*n1*pow(r2,2)*u1*w1 + 
   2*e2*j2*l2*s1*u1*w1 - 2*c2*j2*m2*s1*u1*w1 - 2*e2*i2*q2*s1*u1*w1 + 
   2*pow(m2,2)*q2*s1*u1*w1 + 2*c2*i2*r2*s1*u1*w1 - 2*l2*m2*r2*s1*u1*w1 + 
   2*e1*j2*l2*s2*u1*w1 - 2*c1*j2*m2*s2*u1*w1 - 2*e1*i2*q2*s2*u1*w1 + 
   2*c2*i2*r1*s2*u1*w1 - 2*l2*m2*r1*s2*u1*w1 + 2*c1*i2*r2*s2*u1*w1 - 
   2*d1*e1*pow(j2,2)*u2*w1 + 2*c1*f1*pow(j2,2)*u2*w1 + 
   2*d1*e1*i2*o2*u2*w1 - 2*c1*f1*i2*o2*u2*w1 - 2*e1*l2*n1*o2*u2*w1 + 
   2*c1*m2*n1*o2*u2*w1 + 2*e1*j2*n1*q2*u2*w1 - 2*f1*j2*l2*r1*u2*w1 + 
   4*d1*j2*m2*r1*u2*w1 - 2*c2*j2*n1*r1*u2*w1 - 2*c1*j2*n2*r1*u2*w1 + 
   2*f1*i2*q2*r1*u2*w1 - 2*m2*n1*q2*r1*u2*w1 - 2*d2*i2*pow(r1,2)*u2*w1 + 
   2*l2*n2*pow(r1,2)*u2*w1 - 2*c1*j2*n1*r2*u2*w1 - 4*d1*i2*r1*r2*u2*w1 + 
   4*l2*n1*r1*r2*u2*w1 + 2*e1*j2*l2*s1*u2*w1 - 2*c1*j2*m2*s1*u2*w1 - 
   2*e1*i2*q2*s1*u2*w1 + 2*c2*i2*r1*s1*u2*w1 - 2*l2*m2*r1*s1*u2*w1 + 
   2*c1*i2*r2*s1*u2*w1 + 2*c1*i2*r1*s2*u2*w1 + 2*c2*d1*pow(j2,2)*v1*w1 + 
   2*c1*d2*pow(j2,2)*v1*w1 - 2*c2*d1*i2*o2*v1*w1 - 2*c1*d2*i2*o2*v1*w1 - 
   2*f1*pow(l2,2)*o2*v1*w1 + 2*d1*l2*m2*o2*v1*w1 + 2*c2*l2*n1*o2*v1*w1 + 
   2*c1*l2*n2*o2*v1*w1 + 4*f1*j2*l2*q2*v1*w1 - 2*d1*j2*m2*q2*v1*w1 - 
   2*c2*j2*n1*q2*v1*w1 - 2*c1*j2*n2*q2*v1*w1 - 2*f1*i2*pow(q2,2)*v1*w1 + 
   2*m2*n1*pow(q2,2)*v1*w1 - 2*d2*j2*l2*r1*v1*w1 + 2*d2*i2*q2*r1*v1*w1 - 
   2*l2*n2*q2*r1*v1*w1 - 2*d1*j2*l2*r2*v1*w1 + 2*d1*i2*q2*r2*v1*w1 - 
   2*l2*n1*q2*r2*v1*w1 - 2*c2*j2*l2*s1*v1*w1 + 2*c2*i2*q2*s1*v1*w1 - 
   2*l2*m2*q2*s1*v1*w1 + 2*pow(l2,2)*r2*s1*v1*w1 - 2*c1*j2*l2*s2*v1*w1 + 
   2*c1*i2*q2*s2*v1*w1 + 2*pow(l2,2)*r1*s2*v1*w1 + 
   2*c1*d1*pow(j2,2)*v2*w1 - 2*c1*d1*i2*o2*v2*w1 + 2*c1*l2*n1*o2*v2*w1 - 
   2*c1*j2*n1*q2*v2*w1 - 2*d1*j2*l2*r1*v2*w1 + 2*d1*i2*q2*r1*v2*w1 - 
   2*l2*n1*q2*r1*v2*w1 - 2*c1*j2*l2*s1*v2*w1 + 2*c1*i2*q2*s1*v2*w1 + 
   2*pow(l2,2)*r1*s1*v2*w1 - 2*c1*c2*pow(j2,2)*pow(w1,2) + 
   2*c1*c2*i2*o2*pow(w1,2) + e1*pow(l2,2)*o2*pow(w1,2) - 
   2*c1*l2*m2*o2*pow(w1,2) - 2*e1*j2*l2*q2*pow(w1,2) + 
   2*c1*j2*m2*q2*pow(w1,2) + e1*i2*pow(q2,2)*pow(w1,2) + 
   2*c2*j2*l2*r1*pow(w1,2) - 2*c2*i2*q2*r1*pow(w1,2) + 
   2*l2*m2*q2*r1*pow(w1,2) + 2*c1*j2*l2*r2*pow(w1,2) - 
   2*c1*i2*q2*r2*pow(w1,2) - 2*pow(l2,2)*r1*r2*pow(w1,2) - 
   2*pow(c1,2)*k2*n1*o2*w2 + 2*pow(c1,2)*j2*n1*p2*w2 - 
   2*c1*d1*j2*k2*r1*w2 + 2*c1*d1*i2*p2*r1*w2 - 2*c1*l2*n1*p2*r1*w2 + 
   4*c1*k2*n1*q2*r1*w2 + 2*d1*k2*l2*pow(r1,2)*w2 + 
   2*pow(c1,2)*j2*k2*s1*w2 - 2*pow(c1,2)*i2*p2*s1*w2 - 
   2*c1*k2*l2*r1*s1*w2 - 2*d1*e1*pow(j2,2)*u1*w2 + 
   2*c1*f1*pow(j2,2)*u1*w2 + 2*d1*e1*i2*o2*u1*w2 - 2*c1*f1*i2*o2*u1*w2 - 
   2*e1*l2*n1*o2*u1*w2 + 2*c1*m2*n1*o2*u1*w2 + 2*e1*j2*n1*q2*u1*w2 - 
   2*f1*j2*l2*r1*u1*w2 + 4*d1*j2*m2*r1*u1*w2 - 2*c2*j2*n1*r1*u1*w2 - 
   2*c1*j2*n2*r1*u1*w2 + 2*f1*i2*q2*r1*u1*w2 - 2*m2*n1*q2*r1*u1*w2 - 
   2*d2*i2*pow(r1,2)*u1*w2 + 2*l2*n2*pow(r1,2)*u1*w2 - 
   2*c1*j2*n1*r2*u1*w2 - 4*d1*i2*r1*r2*u1*w2 + 4*l2*n1*r1*r2*u1*w2 + 
   2*e1*j2*l2*s1*u1*w2 - 2*c1*j2*m2*s1*u1*w2 - 2*e1*i2*q2*s1*u1*w2 + 
   2*c2*i2*r1*s1*u1*w2 - 2*l2*m2*r1*s1*u1*w2 + 2*c1*i2*r2*s1*u1*w2 + 
   2*c1*i2*r1*s2*u1*w2 - 2*c1*j2*n1*r1*u2*w2 - 2*d1*i2*pow(r1,2)*u2*w2 + 
   2*l2*n1*pow(r1,2)*u2*w2 + 2*c1*i2*r1*s1*u2*w2 + 
   2*c1*d1*pow(j2,2)*v1*w2 - 2*c1*d1*i2*o2*v1*w2 + 2*c1*l2*n1*o2*v1*w2 - 
   2*c1*j2*n1*q2*v1*w2 - 2*d1*j2*l2*r1*v1*w2 + 2*d1*i2*q2*r1*v1*w2 - 
   2*l2*n1*q2*r1*v1*w2 - 2*c1*j2*l2*s1*v1*w2 + 2*c1*i2*q2*s1*v1*w2 + 
   2*pow(l2,2)*r1*s1*v1*w2 - 2*pow(c1,2)*pow(j2,2)*w1*w2 + 
   2*pow(c1,2)*i2*o2*w1*w2 + 4*c1*j2*l2*r1*w1*w2 - 4*c1*i2*q2*r1*w1*w2 - 
   2*pow(l2,2)*pow(r1,2)*w1*w2 + pow(f1,2)*pow(k2,2)*o2*z1 - 
   e1*g1*pow(k2,2)*o2*z1 - 2*pow(f1,2)*j2*k2*p2*z1 + 
   2*e1*g1*j2*k2*p2*z1 + pow(f1,2)*i2*pow(p2,2)*z1 - 
   e1*g1*i2*pow(p2,2)*z1 - 2*f1*m2*n1*pow(p2,2)*z1 + 
   e2*pow(n1,2)*pow(p2,2)*z1 + 2*e1*n1*n2*pow(p2,2)*z1 - 
   2*g1*k2*m2*p2*r1*z1 + 2*f2*k2*n1*p2*r1*z1 + 2*f1*k2*n2*p2*r1*z1 + 
   g2*pow(k2,2)*pow(r1,2)*z1 + 2*f1*k2*n1*p2*r2*z1 + 
   2*g1*pow(k2,2)*r1*r2*z1 + 2*f1*k2*m2*p2*s1*z1 - 2*e2*k2*n1*p2*s1*z1 - 
   2*e1*k2*n2*p2*s1*z1 - 2*f2*pow(k2,2)*r1*s1*z1 - 
   2*f1*pow(k2,2)*r2*s1*z1 + e2*pow(k2,2)*pow(s1,2)*z1 - 
   2*e1*k2*n1*p2*s2*z1 - 2*f1*pow(k2,2)*r1*s2*z1 + 
   2*e1*pow(k2,2)*s1*s2*z1 + pow(f1,2)*pow(j2,2)*t2*z1 - 
   e1*g1*pow(j2,2)*t2*z1 - pow(f1,2)*i2*o2*t2*z1 + e1*g1*i2*o2*t2*z1 + 
   2*f1*m2*n1*o2*t2*z1 - e2*pow(n1,2)*o2*t2*z1 - 2*e1*n1*n2*o2*t2*z1 + 
   2*g1*j2*m2*r1*t2*z1 - 2*f2*j2*n1*r1*t2*z1 - 2*f1*j2*n2*r1*t2*z1 - 
   g2*i2*pow(r1,2)*t2*z1 + pow(n2,2)*pow(r1,2)*t2*z1 - 
   2*f1*j2*n1*r2*t2*z1 - 2*g1*i2*r1*r2*t2*z1 + 4*n1*n2*r1*r2*t2*z1 + 
   pow(n1,2)*pow(r2,2)*t2*z1 - 2*f1*j2*m2*s1*t2*z1 + 
   2*e2*j2*n1*s1*t2*z1 + 2*e1*j2*n2*s1*t2*z1 + 2*f2*i2*r1*s1*t2*z1 - 
   2*m2*n2*r1*s1*t2*z1 + 2*f1*i2*r2*s1*t2*z1 - 2*m2*n1*r2*s1*t2*z1 - 
   e2*i2*pow(s1,2)*t2*z1 + pow(m2,2)*pow(s1,2)*t2*z1 + 
   2*e1*j2*n1*s2*t2*z1 + 2*f1*i2*r1*s2*t2*z1 - 2*m2*n1*r1*s2*t2*z1 - 
   2*e1*i2*s1*s2*t2*z1 + 2*g1*k2*m2*o2*v1*z1 - 2*f2*k2*n1*o2*v1*z1 - 
   2*f1*k2*n2*o2*v1*z1 - 2*g1*j2*m2*p2*v1*z1 + 2*f2*j2*n1*p2*v1*z1 + 
   2*f1*j2*n2*p2*v1*z1 - 2*g2*j2*k2*r1*v1*z1 + 2*g2*i2*p2*r1*v1*z1 - 
   2*pow(n2,2)*p2*r1*v1*z1 - 2*g1*j2*k2*r2*v1*z1 + 2*g1*i2*p2*r2*v1*z1 - 
   4*n1*n2*p2*r2*v1*z1 + 2*f2*j2*k2*s1*v1*z1 - 2*f2*i2*p2*s1*v1*z1 + 
   2*m2*n2*p2*s1*v1*z1 + 2*k2*n2*r2*s1*v1*z1 + 2*f1*j2*k2*s2*v1*z1 - 
   2*f1*i2*p2*s2*v1*z1 + 2*m2*n1*p2*s2*v1*z1 + 2*k2*n2*r1*s2*v1*z1 + 
   2*k2*n1*r2*s2*v1*z1 - 4*k2*m2*s1*s2*v1*z1 + 
   g2*pow(j2,2)*pow(v1,2)*z1 - g2*i2*o2*pow(v1,2)*z1 + 
   pow(n2,2)*o2*pow(v1,2)*z1 - 2*j2*n2*s2*pow(v1,2)*z1 + 
   i2*pow(s2,2)*pow(v1,2)*z1 - 2*f1*k2*n1*o2*v2*z1 + 
   2*f1*j2*n1*p2*v2*z1 - 2*g1*j2*k2*r1*v2*z1 + 2*g1*i2*p2*r1*v2*z1 - 
   4*n1*n2*p2*r1*v2*z1 - 2*pow(n1,2)*p2*r2*v2*z1 + 2*f1*j2*k2*s1*v2*z1 - 
   2*f1*i2*p2*s1*v2*z1 + 2*m2*n1*p2*s1*v2*z1 + 2*k2*n2*r1*s1*v2*z1 + 
   2*k2*n1*r2*s1*v2*z1 - 2*k2*m2*pow(s1,2)*v2*z1 + 2*k2*n1*r1*s2*v2*z1 + 
   2*g1*pow(j2,2)*v1*v2*z1 - 2*g1*i2*o2*v1*v2*z1 + 4*n1*n2*o2*v1*v2*z1 - 
   4*j2*n2*s1*v1*v2*z1 - 4*j2*n1*s2*v1*v2*z1 + 4*i2*s1*s2*v1*v2*z1 + 
   pow(n1,2)*o2*pow(v2,2)*z1 - 2*j2*n1*s1*pow(v2,2)*z1 + 
   i2*pow(s1,2)*pow(v2,2)*z1 - 2*f1*k2*m2*o2*w1*z1 + 
   2*e2*k2*n1*o2*w1*z1 + 2*e1*k2*n2*o2*w1*z1 + 2*f1*j2*m2*p2*w1*z1 - 
   2*e2*j2*n1*p2*w1*z1 - 2*e1*j2*n2*p2*w1*z1 + 2*f2*j2*k2*r1*w1*z1 - 
   2*f2*i2*p2*r1*w1*z1 + 2*m2*n2*p2*r1*w1*z1 + 2*f1*j2*k2*r2*w1*z1 - 
   2*f1*i2*p2*r2*w1*z1 + 2*m2*n1*p2*r2*w1*z1 - 4*k2*n2*r1*r2*w1*z1 - 
   2*k2*n1*pow(r2,2)*w1*z1 - 2*e2*j2*k2*s1*w1*z1 + 2*e2*i2*p2*s1*w1*z1 - 
   2*pow(m2,2)*p2*s1*w1*z1 + 2*k2*m2*r2*s1*w1*z1 - 2*e1*j2*k2*s2*w1*z1 + 
   2*e1*i2*p2*s2*w1*z1 + 2*k2*m2*r1*s2*w1*z1 - 2*f2*pow(j2,2)*v1*w1*z1 + 
   2*f2*i2*o2*v1*w1*z1 - 2*m2*n2*o2*v1*w1*z1 + 2*j2*n2*r2*v1*w1*z1 + 
   2*j2*m2*s2*v1*w1*z1 - 2*i2*r2*s2*v1*w1*z1 - 2*f1*pow(j2,2)*v2*w1*z1 + 
   2*f1*i2*o2*v2*w1*z1 - 2*m2*n1*o2*v2*w1*z1 + 2*j2*n2*r1*v2*w1*z1 + 
   2*j2*n1*r2*v2*w1*z1 + 2*j2*m2*s1*v2*w1*z1 - 2*i2*r2*s1*v2*w1*z1 - 
   2*i2*r1*s2*v2*w1*z1 + e2*pow(j2,2)*pow(w1,2)*z1 - 
   e2*i2*o2*pow(w1,2)*z1 + pow(m2,2)*o2*pow(w1,2)*z1 - 
   2*j2*m2*r2*pow(w1,2)*z1 + i2*pow(r2,2)*pow(w1,2)*z1 + 
   2*e1*k2*n1*o2*w2*z1 - 2*e1*j2*n1*p2*w2*z1 + 2*f1*j2*k2*r1*w2*z1 - 
   2*f1*i2*p2*r1*w2*z1 + 2*m2*n1*p2*r1*w2*z1 - 2*k2*n2*pow(r1,2)*w2*z1 - 
   4*k2*n1*r1*r2*w2*z1 - 2*e1*j2*k2*s1*w2*z1 + 2*e1*i2*p2*s1*w2*z1 + 
   2*k2*m2*r1*s1*w2*z1 - 2*f1*pow(j2,2)*v1*w2*z1 + 2*f1*i2*o2*v1*w2*z1 - 
   2*m2*n1*o2*v1*w2*z1 + 2*j2*n2*r1*v1*w2*z1 + 2*j2*n1*r2*v1*w2*z1 + 
   2*j2*m2*s1*v1*w2*z1 - 2*i2*r2*s1*v1*w2*z1 - 2*i2*r1*s2*v1*w2*z1 + 
   2*j2*n1*r1*v2*w2*z1 - 2*i2*r1*s1*v2*w2*z1 + 2*e1*pow(j2,2)*w1*w2*z1 - 
   2*e1*i2*o2*w1*w2*z1 - 4*j2*m2*r1*w1*w2*z1 + 4*i2*r1*r2*w1*w2*z1 + 
   i2*pow(r1,2)*pow(w2,2)*z1 + e1*pow(n1,2)*pow(p2,2)*z2 + 
   2*f1*k2*n1*p2*r1*z2 + g1*pow(k2,2)*pow(r1,2)*z2 - 
   2*e1*k2*n1*p2*s1*z2 - 2*f1*pow(k2,2)*r1*s1*z2 + 
   e1*pow(k2,2)*pow(s1,2)*z2 - e1*pow(n1,2)*o2*t2*z2 - 
   2*f1*j2*n1*r1*t2*z2 - g1*i2*pow(r1,2)*t2*z2 + 
   2*n1*n2*pow(r1,2)*t2*z2 + 2*pow(n1,2)*r1*r2*t2*z2 + 
   2*e1*j2*n1*s1*t2*z2 + 2*f1*i2*r1*s1*t2*z2 - 2*m2*n1*r1*s1*t2*z2 - 
   e1*i2*pow(s1,2)*t2*z2 - 2*f1*k2*n1*o2*v1*z2 + 2*f1*j2*n1*p2*v1*z2 - 
   2*g1*j2*k2*r1*v1*z2 + 2*g1*i2*p2*r1*v1*z2 - 4*n1*n2*p2*r1*v1*z2 - 
   2*pow(n1,2)*p2*r2*v1*z2 + 2*f1*j2*k2*s1*v1*z2 - 2*f1*i2*p2*s1*v1*z2 + 
   2*m2*n1*p2*s1*v1*z2 + 2*k2*n2*r1*s1*v1*z2 + 2*k2*n1*r2*s1*v1*z2 - 
   2*k2*m2*pow(s1,2)*v1*z2 + 2*k2*n1*r1*s2*v1*z2 + 
   g1*pow(j2,2)*pow(v1,2)*z2 - g1*i2*o2*pow(v1,2)*z2 + 
   2*n1*n2*o2*pow(v1,2)*z2 - 2*j2*n2*s1*pow(v1,2)*z2 - 
   2*j2*n1*s2*pow(v1,2)*z2 + 2*i2*s1*s2*pow(v1,2)*z2 - 
   2*pow(n1,2)*p2*r1*v2*z2 + 2*k2*n1*r1*s1*v2*z2 + 
   2*pow(n1,2)*o2*v1*v2*z2 - 4*j2*n1*s1*v1*v2*z2 + 
   2*i2*pow(s1,2)*v1*v2*z2 + 2*e1*k2*n1*o2*w1*z2 - 2*e1*j2*n1*p2*w1*z2 + 
   2*f1*j2*k2*r1*w1*z2 - 2*f1*i2*p2*r1*w1*z2 + 2*m2*n1*p2*r1*w1*z2 - 
   2*k2*n2*pow(r1,2)*w1*z2 - 4*k2*n1*r1*r2*w1*z2 - 2*e1*j2*k2*s1*w1*z2 + 
   2*e1*i2*p2*s1*w1*z2 + 2*k2*m2*r1*s1*w1*z2 - 2*f1*pow(j2,2)*v1*w1*z2 + 
   2*f1*i2*o2*v1*w1*z2 - 2*m2*n1*o2*v1*w1*z2 + 2*j2*n2*r1*v1*w1*z2 + 
   2*j2*n1*r2*v1*w1*z2 + 2*j2*m2*s1*v1*w1*z2 - 2*i2*r2*s1*v1*w1*z2 - 
   2*i2*r1*s2*v1*w1*z2 + 2*j2*n1*r1*v2*w1*z2 - 2*i2*r1*s1*v2*w1*z2 + 
   e1*pow(j2,2)*pow(w1,2)*z2 - e1*i2*o2*pow(w1,2)*z2 - 
   2*j2*m2*r1*pow(w1,2)*z2 + 2*i2*r1*r2*pow(w1,2)*z2 - 
   2*k2*n1*pow(r1,2)*w2*z2 + 2*j2*n1*r1*v1*w2*z2 - 2*i2*r1*s1*v1*w2*z2 + 
   2*i2*pow(r1,2)*w1*w2*z2
;

	k[10]  = -(pow(c0,2)*pow(n0,2)*pow(p2,2)) - 2*c0*d0*k2*n0*p2*r0 - 
   pow(d0,2)*pow(k2,2)*pow(r0,2) + 2*pow(c0,2)*k2*n0*p2*s0 + 
   2*c0*d0*pow(k2,2)*r0*s0 - pow(c0,2)*pow(k2,2)*pow(s0,2) + 
   pow(c0,2)*pow(n0,2)*o2*t2 + 2*c0*d0*j2*n0*r0*t2 - 
   2*c0*pow(n0,2)*q2*r0*t2 + pow(d0,2)*i2*pow(r0,2)*t2 - 
   2*d0*l2*n0*pow(r0,2)*t2 - 2*pow(c0,2)*j2*n0*s0*t2 - 
   2*c0*d0*i2*r0*s0*t2 + 2*c0*l2*n0*r0*s0*t2 + 
   pow(c0,2)*i2*pow(s0,2)*t2 - 2*d0*e0*k2*n0*o2*u0 + 
   2*c0*f0*k2*n0*o2*u0 + 2*d0*e0*j2*n0*p2*u0 - 2*c0*f0*j2*n0*p2*u0 - 
   2*e0*pow(n0,2)*p2*q2*u0 - 2*d0*f0*j2*k2*r0*u0 + 2*c0*g0*j2*k2*r0*u0 + 
   2*d0*f0*i2*p2*r0*u0 - 2*c0*g0*i2*p2*r0*u0 - 2*f0*l2*n0*p2*r0*u0 - 
   2*d0*m2*n0*p2*r0*u0 + 2*c2*pow(n0,2)*p2*r0*u0 + 4*c0*n0*n2*p2*r0*u0 - 
   2*f0*k2*n0*q2*r0*u0 - 2*g0*k2*l2*pow(r0,2)*u0 + 
   2*d2*k2*n0*pow(r0,2)*u0 + 2*d0*k2*n2*pow(r0,2)*u0 + 
   2*c0*pow(n0,2)*p2*r2*u0 + 4*d0*k2*n0*r0*r2*u0 + 2*d0*e0*j2*k2*s0*u0 - 
   2*c0*f0*j2*k2*s0*u0 - 2*d0*e0*i2*p2*s0*u0 + 2*c0*f0*i2*p2*s0*u0 + 
   2*e0*l2*n0*p2*s0*u0 - 2*c0*m2*n0*p2*s0*u0 + 2*e0*k2*n0*q2*s0*u0 + 
   4*f0*k2*l2*r0*s0*u0 - 2*d0*k2*m2*r0*s0*u0 - 2*c2*k2*n0*r0*s0*u0 - 
   2*c0*k2*n2*r0*s0*u0 - 2*c0*k2*n0*r2*s0*u0 - 2*e0*k2*l2*pow(s0,2)*u0 + 
   2*c0*k2*m2*pow(s0,2)*u0 - 2*c0*k2*n0*r0*s2*u0 - 
   pow(f0,2)*pow(j2,2)*pow(u0,2) + e0*g0*pow(j2,2)*pow(u0,2) + 
   pow(f0,2)*i2*o2*pow(u0,2) - e0*g0*i2*o2*pow(u0,2) - 
   2*f0*m2*n0*o2*pow(u0,2) + e2*pow(n0,2)*o2*pow(u0,2) + 
   2*e0*n0*n2*o2*pow(u0,2) - 2*g0*j2*m2*r0*pow(u0,2) + 
   2*f2*j2*n0*r0*pow(u0,2) + 2*f0*j2*n2*r0*pow(u0,2) + 
   g2*i2*pow(r0,2)*pow(u0,2) - pow(n2,2)*pow(r0,2)*pow(u0,2) + 
   2*f0*j2*n0*r2*pow(u0,2) + 2*g0*i2*r0*r2*pow(u0,2) - 
   4*n0*n2*r0*r2*pow(u0,2) - pow(n0,2)*pow(r2,2)*pow(u0,2) + 
   2*f0*j2*m2*s0*pow(u0,2) - 2*e2*j2*n0*s0*pow(u0,2) - 
   2*e0*j2*n2*s0*pow(u0,2) - 2*f2*i2*r0*s0*pow(u0,2) + 
   2*m2*n2*r0*s0*pow(u0,2) - 2*f0*i2*r2*s0*pow(u0,2) + 
   2*m2*n0*r2*s0*pow(u0,2) + e2*i2*pow(s0,2)*pow(u0,2) - 
   pow(m2,2)*pow(s0,2)*pow(u0,2) - 2*e0*j2*n0*s2*pow(u0,2) - 
   2*f0*i2*r0*s2*pow(u0,2) + 2*m2*n0*r0*s2*pow(u0,2) + 
   2*e0*i2*s0*s2*pow(u0,2) + 2*c0*pow(n0,2)*p2*r0*u2 + 
   2*d0*k2*n0*pow(r0,2)*u2 - 2*c0*k2*n0*r0*s0*u2 + 
   2*e0*pow(n0,2)*o2*u0*u2 + 4*f0*j2*n0*r0*u0*u2 + 
   2*g0*i2*pow(r0,2)*u0*u2 - 4*n0*n2*pow(r0,2)*u0*u2 - 
   4*pow(n0,2)*r0*r2*u0*u2 - 4*e0*j2*n0*s0*u0*u2 - 4*f0*i2*r0*s0*u0*u2 + 
   4*m2*n0*r0*s0*u0*u2 + 2*e0*i2*pow(s0,2)*u0*u2 - 
   pow(n0,2)*pow(r0,2)*pow(u2,2) + 2*c0*d0*k2*n0*o2*v0 - 
   2*c0*d0*j2*n0*p2*v0 + 2*c0*pow(n0,2)*p2*q2*v0 + 
   2*pow(d0,2)*j2*k2*r0*v0 - 2*pow(d0,2)*i2*p2*r0*v0 + 
   4*d0*l2*n0*p2*r0*v0 - 2*d0*k2*n0*q2*r0*v0 - 2*c0*d0*j2*k2*s0*v0 + 
   2*c0*d0*i2*p2*s0*v0 - 2*c0*l2*n0*p2*s0*v0 - 2*c0*k2*n0*q2*s0*v0 - 
   2*d0*k2*l2*r0*s0*v0 + 2*c0*k2*l2*pow(s0,2)*v0 + 
   2*d0*f0*pow(j2,2)*u0*v0 - 2*c0*g0*pow(j2,2)*u0*v0 - 
   2*d0*f0*i2*o2*u0*v0 + 2*c0*g0*i2*o2*u0*v0 + 2*f0*l2*n0*o2*u0*v0 + 
   2*d0*m2*n0*o2*u0*v0 - 2*c2*pow(n0,2)*o2*u0*v0 - 4*c0*n0*n2*o2*u0*v0 - 
   2*f0*j2*n0*q2*u0*v0 + 2*g0*j2*l2*r0*u0*v0 - 2*d2*j2*n0*r0*u0*v0 - 
   2*d0*j2*n2*r0*u0*v0 - 2*g0*i2*q2*r0*u0*v0 + 4*n0*n2*q2*r0*u0*v0 - 
   2*d0*j2*n0*r2*u0*v0 + 2*pow(n0,2)*q2*r2*u0*v0 - 2*f0*j2*l2*s0*u0*v0 - 
   2*d0*j2*m2*s0*u0*v0 + 4*c2*j2*n0*s0*u0*v0 + 4*c0*j2*n2*s0*u0*v0 + 
   2*f0*i2*q2*s0*u0*v0 - 2*m2*n0*q2*s0*u0*v0 + 2*d2*i2*r0*s0*u0*v0 - 
   2*l2*n2*r0*s0*u0*v0 + 2*d0*i2*r2*s0*u0*v0 - 2*l2*n0*r2*s0*u0*v0 - 
   2*c2*i2*pow(s0,2)*u0*v0 + 2*l2*m2*pow(s0,2)*u0*v0 + 
   4*c0*j2*n0*s2*u0*v0 + 2*d0*i2*r0*s2*u0*v0 - 2*l2*n0*r0*s2*u0*v0 - 
   4*c0*i2*s0*s2*u0*v0 - 2*c0*pow(n0,2)*o2*u2*v0 - 2*d0*j2*n0*r0*u2*v0 + 
   2*pow(n0,2)*q2*r0*u2*v0 + 4*c0*j2*n0*s0*u2*v0 + 2*d0*i2*r0*s0*u2*v0 - 
   2*l2*n0*r0*s0*u2*v0 - 2*c0*i2*pow(s0,2)*u2*v0 - 
   pow(d0,2)*pow(j2,2)*pow(v0,2) + pow(d0,2)*i2*o2*pow(v0,2) - 
   2*d0*l2*n0*o2*pow(v0,2) + 2*d0*j2*n0*q2*pow(v0,2) - 
   pow(n0,2)*pow(q2,2)*pow(v0,2) + 2*d0*j2*l2*s0*pow(v0,2) - 
   2*d0*i2*q2*s0*pow(v0,2) + 2*l2*n0*q2*s0*pow(v0,2) - 
   pow(l2,2)*pow(s0,2)*pow(v0,2) - 2*c0*pow(n0,2)*o2*u0*v2 - 
   2*d0*j2*n0*r0*u0*v2 + 2*pow(n0,2)*q2*r0*u0*v2 + 4*c0*j2*n0*s0*u0*v2 + 
   2*d0*i2*r0*s0*u0*v2 - 2*l2*n0*r0*s0*u0*v2 - 2*c0*i2*pow(s0,2)*u0*v2 - 
   2*pow(c0,2)*k2*n0*o2*w0 + 2*pow(c0,2)*j2*n0*p2*w0 - 
   2*c0*d0*j2*k2*r0*w0 + 2*c0*d0*i2*p2*r0*w0 - 2*c0*l2*n0*p2*r0*w0 + 
   4*c0*k2*n0*q2*r0*w0 + 2*d0*k2*l2*pow(r0,2)*w0 + 
   2*pow(c0,2)*j2*k2*s0*w0 - 2*pow(c0,2)*i2*p2*s0*w0 - 
   2*c0*k2*l2*r0*s0*w0 - 2*d0*e0*pow(j2,2)*u0*w0 + 
   2*c0*f0*pow(j2,2)*u0*w0 + 2*d0*e0*i2*o2*u0*w0 - 2*c0*f0*i2*o2*u0*w0 - 
   2*e0*l2*n0*o2*u0*w0 + 2*c0*m2*n0*o2*u0*w0 + 2*e0*j2*n0*q2*u0*w0 - 
   2*f0*j2*l2*r0*u0*w0 + 4*d0*j2*m2*r0*u0*w0 - 2*c2*j2*n0*r0*u0*w0 - 
   2*c0*j2*n2*r0*u0*w0 + 2*f0*i2*q2*r0*u0*w0 - 2*m2*n0*q2*r0*u0*w0 - 
   2*d2*i2*pow(r0,2)*u0*w0 + 2*l2*n2*pow(r0,2)*u0*w0 - 
   2*c0*j2*n0*r2*u0*w0 - 4*d0*i2*r0*r2*u0*w0 + 4*l2*n0*r0*r2*u0*w0 + 
   2*e0*j2*l2*s0*u0*w0 - 2*c0*j2*m2*s0*u0*w0 - 2*e0*i2*q2*s0*u0*w0 + 
   2*c2*i2*r0*s0*u0*w0 - 2*l2*m2*r0*s0*u0*w0 + 2*c0*i2*r2*s0*u0*w0 + 
   2*c0*i2*r0*s2*u0*w0 - 2*c0*j2*n0*r0*u2*w0 - 2*d0*i2*pow(r0,2)*u2*w0 + 
   2*l2*n0*pow(r0,2)*u2*w0 + 2*c0*i2*r0*s0*u2*w0 + 
   2*c0*d0*pow(j2,2)*v0*w0 - 2*c0*d0*i2*o2*v0*w0 + 2*c0*l2*n0*o2*v0*w0 - 
   2*c0*j2*n0*q2*v0*w0 - 2*d0*j2*l2*r0*v0*w0 + 2*d0*i2*q2*r0*v0*w0 - 
   2*l2*n0*q2*r0*v0*w0 - 2*c0*j2*l2*s0*v0*w0 + 2*c0*i2*q2*s0*v0*w0 + 
   2*pow(l2,2)*r0*s0*v0*w0 - pow(c0,2)*pow(j2,2)*pow(w0,2) + 
   pow(c0,2)*i2*o2*pow(w0,2) + 2*c0*j2*l2*r0*pow(w0,2) - 
   2*c0*i2*q2*r0*pow(w0,2) - pow(l2,2)*pow(r0,2)*pow(w0,2) - 
   2*c0*j2*n0*r0*u0*w2 - 2*d0*i2*pow(r0,2)*u0*w2 + 
   2*l2*n0*pow(r0,2)*u0*w2 + 2*c0*i2*r0*s0*u0*w2 + 
   e0*pow(n0,2)*pow(p2,2)*z0 + 2*f0*k2*n0*p2*r0*z0 + 
   g0*pow(k2,2)*pow(r0,2)*z0 - 2*e0*k2*n0*p2*s0*z0 - 
   2*f0*pow(k2,2)*r0*s0*z0 + e0*pow(k2,2)*pow(s0,2)*z0 - 
   e0*pow(n0,2)*o2*t2*z0 - 2*f0*j2*n0*r0*t2*z0 - g0*i2*pow(r0,2)*t2*z0 + 
   2*n0*n2*pow(r0,2)*t2*z0 + 2*pow(n0,2)*r0*r2*t2*z0 + 
   2*e0*j2*n0*s0*t2*z0 + 2*f0*i2*r0*s0*t2*z0 - 2*m2*n0*r0*s0*t2*z0 - 
   e0*i2*pow(s0,2)*t2*z0 - 2*f0*k2*n0*o2*v0*z0 + 2*f0*j2*n0*p2*v0*z0 - 
   2*g0*j2*k2*r0*v0*z0 + 2*g0*i2*p2*r0*v0*z0 - 4*n0*n2*p2*r0*v0*z0 - 
   2*pow(n0,2)*p2*r2*v0*z0 + 2*f0*j2*k2*s0*v0*z0 - 2*f0*i2*p2*s0*v0*z0 + 
   2*m2*n0*p2*s0*v0*z0 + 2*k2*n2*r0*s0*v0*z0 + 2*k2*n0*r2*s0*v0*z0 - 
   2*k2*m2*pow(s0,2)*v0*z0 + 2*k2*n0*r0*s2*v0*z0 + 
   g0*pow(j2,2)*pow(v0,2)*z0 - g0*i2*o2*pow(v0,2)*z0 + 
   2*n0*n2*o2*pow(v0,2)*z0 - 2*j2*n2*s0*pow(v0,2)*z0 - 
   2*j2*n0*s2*pow(v0,2)*z0 + 2*i2*s0*s2*pow(v0,2)*z0 - 
   2*pow(n0,2)*p2*r0*v2*z0 + 2*k2*n0*r0*s0*v2*z0 + 
   2*pow(n0,2)*o2*v0*v2*z0 - 4*j2*n0*s0*v0*v2*z0 + 
   2*i2*pow(s0,2)*v0*v2*z0 + 2*e0*k2*n0*o2*w0*z0 - 2*e0*j2*n0*p2*w0*z0 + 
   2*f0*j2*k2*r0*w0*z0 - 2*f0*i2*p2*r0*w0*z0 + 2*m2*n0*p2*r0*w0*z0 - 
   2*k2*n2*pow(r0,2)*w0*z0 - 4*k2*n0*r0*r2*w0*z0 - 2*e0*j2*k2*s0*w0*z0 + 
   2*e0*i2*p2*s0*w0*z0 + 2*k2*m2*r0*s0*w0*z0 - 2*f0*pow(j2,2)*v0*w0*z0 + 
   2*f0*i2*o2*v0*w0*z0 - 2*m2*n0*o2*v0*w0*z0 + 2*j2*n2*r0*v0*w0*z0 + 
   2*j2*n0*r2*v0*w0*z0 + 2*j2*m2*s0*v0*w0*z0 - 2*i2*r2*s0*v0*w0*z0 - 
   2*i2*r0*s2*v0*w0*z0 + 2*j2*n0*r0*v2*w0*z0 - 2*i2*r0*s0*v2*w0*z0 + 
   e0*pow(j2,2)*pow(w0,2)*z0 - e0*i2*o2*pow(w0,2)*z0 - 
   2*j2*m2*r0*pow(w0,2)*z0 + 2*i2*r0*r2*pow(w0,2)*z0 - 
   2*k2*n0*pow(r0,2)*w2*z0 + 2*j2*n0*r0*v0*w2*z0 - 2*i2*r0*s0*v0*w2*z0 + 
   2*i2*pow(r0,2)*w0*w2*z0 + pow(n0,2)*pow(r0,2)*t2*z2 - 
   2*pow(n0,2)*p2*r0*v0*z2 + 2*k2*n0*r0*s0*v0*z2 + 
   pow(n0,2)*o2*pow(v0,2)*z2 - 2*j2*n0*s0*pow(v0,2)*z2 + 
   i2*pow(s0,2)*pow(v0,2)*z2 - 2*k2*n0*pow(r0,2)*w0*z2 + 
   2*j2*n0*r0*v0*w0*z2 - 2*i2*r0*s0*v0*w0*z2 + i2*pow(r0,2)*pow(w0,2)*z2
;

	k[11]  = -2*c0*c1*pow(n0,2)*pow(p2,2) - 2*pow(c0,2)*n0*n1*pow(p2,2) - 
   2*c1*d0*k2*n0*p2*r0 - 2*c0*d1*k2*n0*p2*r0 - 2*c0*d0*k2*n1*p2*r0 - 
   2*d0*d1*pow(k2,2)*pow(r0,2) - 2*c0*d0*k2*n0*p2*r1 - 
   2*pow(d0,2)*pow(k2,2)*r0*r1 + 4*c0*c1*k2*n0*p2*s0 + 
   2*pow(c0,2)*k2*n1*p2*s0 + 2*c1*d0*pow(k2,2)*r0*s0 + 
   2*c0*d1*pow(k2,2)*r0*s0 + 2*c0*d0*pow(k2,2)*r1*s0 - 
   2*c0*c1*pow(k2,2)*pow(s0,2) + 2*pow(c0,2)*k2*n0*p2*s1 + 
   2*c0*d0*pow(k2,2)*r0*s1 - 2*pow(c0,2)*pow(k2,2)*s0*s1 + 
   2*c0*c1*pow(n0,2)*o2*t2 + 2*pow(c0,2)*n0*n1*o2*t2 + 
   2*c1*d0*j2*n0*r0*t2 + 2*c0*d1*j2*n0*r0*t2 + 2*c0*d0*j2*n1*r0*t2 - 
   2*c1*pow(n0,2)*q2*r0*t2 - 4*c0*n0*n1*q2*r0*t2 + 
   2*d0*d1*i2*pow(r0,2)*t2 - 2*d1*l2*n0*pow(r0,2)*t2 - 
   2*d0*l2*n1*pow(r0,2)*t2 + 2*c0*d0*j2*n0*r1*t2 - 
   2*c0*pow(n0,2)*q2*r1*t2 + 2*pow(d0,2)*i2*r0*r1*t2 - 
   4*d0*l2*n0*r0*r1*t2 - 4*c0*c1*j2*n0*s0*t2 - 2*pow(c0,2)*j2*n1*s0*t2 - 
   2*c1*d0*i2*r0*s0*t2 - 2*c0*d1*i2*r0*s0*t2 + 2*c1*l2*n0*r0*s0*t2 + 
   2*c0*l2*n1*r0*s0*t2 - 2*c0*d0*i2*r1*s0*t2 + 2*c0*l2*n0*r1*s0*t2 + 
   2*c0*c1*i2*pow(s0,2)*t2 - 2*pow(c0,2)*j2*n0*s1*t2 - 
   2*c0*d0*i2*r0*s1*t2 + 2*c0*l2*n0*r0*s1*t2 + 2*pow(c0,2)*i2*s0*s1*t2 - 
   2*d1*e0*k2*n0*o2*u0 - 2*d0*e1*k2*n0*o2*u0 + 2*c1*f0*k2*n0*o2*u0 + 
   2*c0*f1*k2*n0*o2*u0 - 2*d0*e0*k2*n1*o2*u0 + 2*c0*f0*k2*n1*o2*u0 + 
   2*d1*e0*j2*n0*p2*u0 + 2*d0*e1*j2*n0*p2*u0 - 2*c1*f0*j2*n0*p2*u0 - 
   2*c0*f1*j2*n0*p2*u0 + 2*d0*e0*j2*n1*p2*u0 - 2*c0*f0*j2*n1*p2*u0 - 
   2*e1*pow(n0,2)*p2*q2*u0 - 4*e0*n0*n1*p2*q2*u0 - 2*d1*f0*j2*k2*r0*u0 - 
   2*d0*f1*j2*k2*r0*u0 + 2*c1*g0*j2*k2*r0*u0 + 2*c0*g1*j2*k2*r0*u0 + 
   2*d1*f0*i2*p2*r0*u0 + 2*d0*f1*i2*p2*r0*u0 - 2*c1*g0*i2*p2*r0*u0 - 
   2*c0*g1*i2*p2*r0*u0 - 2*f1*l2*n0*p2*r0*u0 - 2*d1*m2*n0*p2*r0*u0 - 
   2*f0*l2*n1*p2*r0*u0 - 2*d0*m2*n1*p2*r0*u0 + 4*c2*n0*n1*p2*r0*u0 + 
   4*c1*n0*n2*p2*r0*u0 + 4*c0*n1*n2*p2*r0*u0 - 2*f1*k2*n0*q2*r0*u0 - 
   2*f0*k2*n1*q2*r0*u0 - 2*g1*k2*l2*pow(r0,2)*u0 + 
   2*d2*k2*n1*pow(r0,2)*u0 + 2*d1*k2*n2*pow(r0,2)*u0 - 
   2*d0*f0*j2*k2*r1*u0 + 2*c0*g0*j2*k2*r1*u0 + 2*d0*f0*i2*p2*r1*u0 - 
   2*c0*g0*i2*p2*r1*u0 - 2*f0*l2*n0*p2*r1*u0 - 2*d0*m2*n0*p2*r1*u0 + 
   2*c2*pow(n0,2)*p2*r1*u0 + 4*c0*n0*n2*p2*r1*u0 - 2*f0*k2*n0*q2*r1*u0 - 
   4*g0*k2*l2*r0*r1*u0 + 4*d2*k2*n0*r0*r1*u0 + 4*d0*k2*n2*r0*r1*u0 + 
   2*c1*pow(n0,2)*p2*r2*u0 + 4*c0*n0*n1*p2*r2*u0 + 4*d1*k2*n0*r0*r2*u0 + 
   4*d0*k2*n1*r0*r2*u0 + 4*d0*k2*n0*r1*r2*u0 + 2*d1*e0*j2*k2*s0*u0 + 
   2*d0*e1*j2*k2*s0*u0 - 2*c1*f0*j2*k2*s0*u0 - 2*c0*f1*j2*k2*s0*u0 - 
   2*d1*e0*i2*p2*s0*u0 - 2*d0*e1*i2*p2*s0*u0 + 2*c1*f0*i2*p2*s0*u0 + 
   2*c0*f1*i2*p2*s0*u0 + 2*e1*l2*n0*p2*s0*u0 - 2*c1*m2*n0*p2*s0*u0 + 
   2*e0*l2*n1*p2*s0*u0 - 2*c0*m2*n1*p2*s0*u0 + 2*e1*k2*n0*q2*s0*u0 + 
   2*e0*k2*n1*q2*s0*u0 + 4*f1*k2*l2*r0*s0*u0 - 2*d1*k2*m2*r0*s0*u0 - 
   2*c2*k2*n1*r0*s0*u0 - 2*c1*k2*n2*r0*s0*u0 + 4*f0*k2*l2*r1*s0*u0 - 
   2*d0*k2*m2*r1*s0*u0 - 2*c2*k2*n0*r1*s0*u0 - 2*c0*k2*n2*r1*s0*u0 - 
   2*c1*k2*n0*r2*s0*u0 - 2*c0*k2*n1*r2*s0*u0 - 2*e1*k2*l2*pow(s0,2)*u0 + 
   2*c1*k2*m2*pow(s0,2)*u0 + 2*d0*e0*j2*k2*s1*u0 - 2*c0*f0*j2*k2*s1*u0 - 
   2*d0*e0*i2*p2*s1*u0 + 2*c0*f0*i2*p2*s1*u0 + 2*e0*l2*n0*p2*s1*u0 - 
   2*c0*m2*n0*p2*s1*u0 + 2*e0*k2*n0*q2*s1*u0 + 4*f0*k2*l2*r0*s1*u0 - 
   2*d0*k2*m2*r0*s1*u0 - 2*c2*k2*n0*r0*s1*u0 - 2*c0*k2*n2*r0*s1*u0 - 
   2*c0*k2*n0*r2*s1*u0 - 4*e0*k2*l2*s0*s1*u0 + 4*c0*k2*m2*s0*s1*u0 - 
   2*c1*k2*n0*r0*s2*u0 - 2*c0*k2*n1*r0*s2*u0 - 2*c0*k2*n0*r1*s2*u0 - 
   2*f0*f1*pow(j2,2)*pow(u0,2) + e1*g0*pow(j2,2)*pow(u0,2) + 
   e0*g1*pow(j2,2)*pow(u0,2) + 2*f0*f1*i2*o2*pow(u0,2) - 
   e1*g0*i2*o2*pow(u0,2) - e0*g1*i2*o2*pow(u0,2) - 
   2*f1*m2*n0*o2*pow(u0,2) - 2*f0*m2*n1*o2*pow(u0,2) + 
   2*e2*n0*n1*o2*pow(u0,2) + 2*e1*n0*n2*o2*pow(u0,2) + 
   2*e0*n1*n2*o2*pow(u0,2) - 2*g1*j2*m2*r0*pow(u0,2) + 
   2*f2*j2*n1*r0*pow(u0,2) + 2*f1*j2*n2*r0*pow(u0,2) - 
   2*g0*j2*m2*r1*pow(u0,2) + 2*f2*j2*n0*r1*pow(u0,2) + 
   2*f0*j2*n2*r1*pow(u0,2) + 2*g2*i2*r0*r1*pow(u0,2) - 
   2*pow(n2,2)*r0*r1*pow(u0,2) + 2*f1*j2*n0*r2*pow(u0,2) + 
   2*f0*j2*n1*r2*pow(u0,2) + 2*g1*i2*r0*r2*pow(u0,2) - 
   4*n1*n2*r0*r2*pow(u0,2) + 2*g0*i2*r1*r2*pow(u0,2) - 
   4*n0*n2*r1*r2*pow(u0,2) - 2*n0*n1*pow(r2,2)*pow(u0,2) + 
   2*f1*j2*m2*s0*pow(u0,2) - 2*e2*j2*n1*s0*pow(u0,2) - 
   2*e1*j2*n2*s0*pow(u0,2) - 2*f2*i2*r1*s0*pow(u0,2) + 
   2*m2*n2*r1*s0*pow(u0,2) - 2*f1*i2*r2*s0*pow(u0,2) + 
   2*m2*n1*r2*s0*pow(u0,2) + 2*f0*j2*m2*s1*pow(u0,2) - 
   2*e2*j2*n0*s1*pow(u0,2) - 2*e0*j2*n2*s1*pow(u0,2) - 
   2*f2*i2*r0*s1*pow(u0,2) + 2*m2*n2*r0*s1*pow(u0,2) - 
   2*f0*i2*r2*s1*pow(u0,2) + 2*m2*n0*r2*s1*pow(u0,2) + 
   2*e2*i2*s0*s1*pow(u0,2) - 2*pow(m2,2)*s0*s1*pow(u0,2) - 
   2*e1*j2*n0*s2*pow(u0,2) - 2*e0*j2*n1*s2*pow(u0,2) - 
   2*f1*i2*r0*s2*pow(u0,2) + 2*m2*n1*r0*s2*pow(u0,2) - 
   2*f0*i2*r1*s2*pow(u0,2) + 2*m2*n0*r1*s2*pow(u0,2) + 
   2*e1*i2*s0*s2*pow(u0,2) + 2*e0*i2*s1*s2*pow(u0,2) - 
   2*d0*e0*k2*n0*o2*u1 + 2*c0*f0*k2*n0*o2*u1 + 2*d0*e0*j2*n0*p2*u1 - 
   2*c0*f0*j2*n0*p2*u1 - 2*e0*pow(n0,2)*p2*q2*u1 - 2*d0*f0*j2*k2*r0*u1 + 
   2*c0*g0*j2*k2*r0*u1 + 2*d0*f0*i2*p2*r0*u1 - 2*c0*g0*i2*p2*r0*u1 - 
   2*f0*l2*n0*p2*r0*u1 - 2*d0*m2*n0*p2*r0*u1 + 2*c2*pow(n0,2)*p2*r0*u1 + 
   4*c0*n0*n2*p2*r0*u1 - 2*f0*k2*n0*q2*r0*u1 - 2*g0*k2*l2*pow(r0,2)*u1 + 
   2*d2*k2*n0*pow(r0,2)*u1 + 2*d0*k2*n2*pow(r0,2)*u1 + 
   2*c0*pow(n0,2)*p2*r2*u1 + 4*d0*k2*n0*r0*r2*u1 + 2*d0*e0*j2*k2*s0*u1 - 
   2*c0*f0*j2*k2*s0*u1 - 2*d0*e0*i2*p2*s0*u1 + 2*c0*f0*i2*p2*s0*u1 + 
   2*e0*l2*n0*p2*s0*u1 - 2*c0*m2*n0*p2*s0*u1 + 2*e0*k2*n0*q2*s0*u1 + 
   4*f0*k2*l2*r0*s0*u1 - 2*d0*k2*m2*r0*s0*u1 - 2*c2*k2*n0*r0*s0*u1 - 
   2*c0*k2*n2*r0*s0*u1 - 2*c0*k2*n0*r2*s0*u1 - 2*e0*k2*l2*pow(s0,2)*u1 + 
   2*c0*k2*m2*pow(s0,2)*u1 - 2*c0*k2*n0*r0*s2*u1 - 
   2*pow(f0,2)*pow(j2,2)*u0*u1 + 2*e0*g0*pow(j2,2)*u0*u1 + 
   2*pow(f0,2)*i2*o2*u0*u1 - 2*e0*g0*i2*o2*u0*u1 - 4*f0*m2*n0*o2*u0*u1 + 
   2*e2*pow(n0,2)*o2*u0*u1 + 4*e0*n0*n2*o2*u0*u1 - 4*g0*j2*m2*r0*u0*u1 + 
   4*f2*j2*n0*r0*u0*u1 + 4*f0*j2*n2*r0*u0*u1 + 2*g2*i2*pow(r0,2)*u0*u1 - 
   2*pow(n2,2)*pow(r0,2)*u0*u1 + 4*f0*j2*n0*r2*u0*u1 + 
   4*g0*i2*r0*r2*u0*u1 - 8*n0*n2*r0*r2*u0*u1 - 
   2*pow(n0,2)*pow(r2,2)*u0*u1 + 4*f0*j2*m2*s0*u0*u1 - 
   4*e2*j2*n0*s0*u0*u1 - 4*e0*j2*n2*s0*u0*u1 - 4*f2*i2*r0*s0*u0*u1 + 
   4*m2*n2*r0*s0*u0*u1 - 4*f0*i2*r2*s0*u0*u1 + 4*m2*n0*r2*s0*u0*u1 + 
   2*e2*i2*pow(s0,2)*u0*u1 - 2*pow(m2,2)*pow(s0,2)*u0*u1 - 
   4*e0*j2*n0*s2*u0*u1 - 4*f0*i2*r0*s2*u0*u1 + 4*m2*n0*r0*s2*u0*u1 + 
   4*e0*i2*s0*s2*u0*u1 + 2*c1*pow(n0,2)*p2*r0*u2 + 4*c0*n0*n1*p2*r0*u2 + 
   2*d1*k2*n0*pow(r0,2)*u2 + 2*d0*k2*n1*pow(r0,2)*u2 + 
   2*c0*pow(n0,2)*p2*r1*u2 + 4*d0*k2*n0*r0*r1*u2 - 2*c1*k2*n0*r0*s0*u2 - 
   2*c0*k2*n1*r0*s0*u2 - 2*c0*k2*n0*r1*s0*u2 - 2*c0*k2*n0*r0*s1*u2 + 
   2*e1*pow(n0,2)*o2*u0*u2 + 4*e0*n0*n1*o2*u0*u2 + 4*f1*j2*n0*r0*u0*u2 + 
   4*f0*j2*n1*r0*u0*u2 + 2*g1*i2*pow(r0,2)*u0*u2 - 
   4*n1*n2*pow(r0,2)*u0*u2 + 4*f0*j2*n0*r1*u0*u2 + 4*g0*i2*r0*r1*u0*u2 - 
   8*n0*n2*r0*r1*u0*u2 - 8*n0*n1*r0*r2*u0*u2 - 4*pow(n0,2)*r1*r2*u0*u2 - 
   4*e1*j2*n0*s0*u0*u2 - 4*e0*j2*n1*s0*u0*u2 - 4*f1*i2*r0*s0*u0*u2 + 
   4*m2*n1*r0*s0*u0*u2 - 4*f0*i2*r1*s0*u0*u2 + 4*m2*n0*r1*s0*u0*u2 + 
   2*e1*i2*pow(s0,2)*u0*u2 - 4*e0*j2*n0*s1*u0*u2 - 4*f0*i2*r0*s1*u0*u2 + 
   4*m2*n0*r0*s1*u0*u2 + 4*e0*i2*s0*s1*u0*u2 + 2*e0*pow(n0,2)*o2*u1*u2 + 
   4*f0*j2*n0*r0*u1*u2 + 2*g0*i2*pow(r0,2)*u1*u2 - 
   4*n0*n2*pow(r0,2)*u1*u2 - 4*pow(n0,2)*r0*r2*u1*u2 - 
   4*e0*j2*n0*s0*u1*u2 - 4*f0*i2*r0*s0*u1*u2 + 4*m2*n0*r0*s0*u1*u2 + 
   2*e0*i2*pow(s0,2)*u1*u2 - 2*n0*n1*pow(r0,2)*pow(u2,2) - 
   2*pow(n0,2)*r0*r1*pow(u2,2) + 2*c1*d0*k2*n0*o2*v0 + 
   2*c0*d1*k2*n0*o2*v0 + 2*c0*d0*k2*n1*o2*v0 - 2*c1*d0*j2*n0*p2*v0 - 
   2*c0*d1*j2*n0*p2*v0 - 2*c0*d0*j2*n1*p2*v0 + 2*c1*pow(n0,2)*p2*q2*v0 + 
   4*c0*n0*n1*p2*q2*v0 + 4*d0*d1*j2*k2*r0*v0 - 4*d0*d1*i2*p2*r0*v0 + 
   4*d1*l2*n0*p2*r0*v0 + 4*d0*l2*n1*p2*r0*v0 - 2*d1*k2*n0*q2*r0*v0 - 
   2*d0*k2*n1*q2*r0*v0 + 2*pow(d0,2)*j2*k2*r1*v0 - 
   2*pow(d0,2)*i2*p2*r1*v0 + 4*d0*l2*n0*p2*r1*v0 - 2*d0*k2*n0*q2*r1*v0 - 
   2*c1*d0*j2*k2*s0*v0 - 2*c0*d1*j2*k2*s0*v0 + 2*c1*d0*i2*p2*s0*v0 + 
   2*c0*d1*i2*p2*s0*v0 - 2*c1*l2*n0*p2*s0*v0 - 2*c0*l2*n1*p2*s0*v0 - 
   2*c1*k2*n0*q2*s0*v0 - 2*c0*k2*n1*q2*s0*v0 - 2*d1*k2*l2*r0*s0*v0 - 
   2*d0*k2*l2*r1*s0*v0 + 2*c1*k2*l2*pow(s0,2)*v0 - 2*c0*d0*j2*k2*s1*v0 + 
   2*c0*d0*i2*p2*s1*v0 - 2*c0*l2*n0*p2*s1*v0 - 2*c0*k2*n0*q2*s1*v0 - 
   2*d0*k2*l2*r0*s1*v0 + 4*c0*k2*l2*s0*s1*v0 + 2*d1*f0*pow(j2,2)*u0*v0 + 
   2*d0*f1*pow(j2,2)*u0*v0 - 2*c1*g0*pow(j2,2)*u0*v0 - 
   2*c0*g1*pow(j2,2)*u0*v0 - 2*d1*f0*i2*o2*u0*v0 - 2*d0*f1*i2*o2*u0*v0 + 
   2*c1*g0*i2*o2*u0*v0 + 2*c0*g1*i2*o2*u0*v0 + 2*f1*l2*n0*o2*u0*v0 + 
   2*d1*m2*n0*o2*u0*v0 + 2*f0*l2*n1*o2*u0*v0 + 2*d0*m2*n1*o2*u0*v0 - 
   4*c2*n0*n1*o2*u0*v0 - 4*c1*n0*n2*o2*u0*v0 - 4*c0*n1*n2*o2*u0*v0 - 
   2*f1*j2*n0*q2*u0*v0 - 2*f0*j2*n1*q2*u0*v0 + 2*g1*j2*l2*r0*u0*v0 - 
   2*d2*j2*n1*r0*u0*v0 - 2*d1*j2*n2*r0*u0*v0 - 2*g1*i2*q2*r0*u0*v0 + 
   4*n1*n2*q2*r0*u0*v0 + 2*g0*j2*l2*r1*u0*v0 - 2*d2*j2*n0*r1*u0*v0 - 
   2*d0*j2*n2*r1*u0*v0 - 2*g0*i2*q2*r1*u0*v0 + 4*n0*n2*q2*r1*u0*v0 - 
   2*d1*j2*n0*r2*u0*v0 - 2*d0*j2*n1*r2*u0*v0 + 4*n0*n1*q2*r2*u0*v0 - 
   2*f1*j2*l2*s0*u0*v0 - 2*d1*j2*m2*s0*u0*v0 + 4*c2*j2*n1*s0*u0*v0 + 
   4*c1*j2*n2*s0*u0*v0 + 2*f1*i2*q2*s0*u0*v0 - 2*m2*n1*q2*s0*u0*v0 + 
   2*d2*i2*r1*s0*u0*v0 - 2*l2*n2*r1*s0*u0*v0 + 2*d1*i2*r2*s0*u0*v0 - 
   2*l2*n1*r2*s0*u0*v0 - 2*f0*j2*l2*s1*u0*v0 - 2*d0*j2*m2*s1*u0*v0 + 
   4*c2*j2*n0*s1*u0*v0 + 4*c0*j2*n2*s1*u0*v0 + 2*f0*i2*q2*s1*u0*v0 - 
   2*m2*n0*q2*s1*u0*v0 + 2*d2*i2*r0*s1*u0*v0 - 2*l2*n2*r0*s1*u0*v0 + 
   2*d0*i2*r2*s1*u0*v0 - 2*l2*n0*r2*s1*u0*v0 - 4*c2*i2*s0*s1*u0*v0 + 
   4*l2*m2*s0*s1*u0*v0 + 4*c1*j2*n0*s2*u0*v0 + 4*c0*j2*n1*s2*u0*v0 + 
   2*d1*i2*r0*s2*u0*v0 - 2*l2*n1*r0*s2*u0*v0 + 2*d0*i2*r1*s2*u0*v0 - 
   2*l2*n0*r1*s2*u0*v0 - 4*c1*i2*s0*s2*u0*v0 - 4*c0*i2*s1*s2*u0*v0 + 
   2*d0*f0*pow(j2,2)*u1*v0 - 2*c0*g0*pow(j2,2)*u1*v0 - 
   2*d0*f0*i2*o2*u1*v0 + 2*c0*g0*i2*o2*u1*v0 + 2*f0*l2*n0*o2*u1*v0 + 
   2*d0*m2*n0*o2*u1*v0 - 2*c2*pow(n0,2)*o2*u1*v0 - 4*c0*n0*n2*o2*u1*v0 - 
   2*f0*j2*n0*q2*u1*v0 + 2*g0*j2*l2*r0*u1*v0 - 2*d2*j2*n0*r0*u1*v0 - 
   2*d0*j2*n2*r0*u1*v0 - 2*g0*i2*q2*r0*u1*v0 + 4*n0*n2*q2*r0*u1*v0 - 
   2*d0*j2*n0*r2*u1*v0 + 2*pow(n0,2)*q2*r2*u1*v0 - 2*f0*j2*l2*s0*u1*v0 - 
   2*d0*j2*m2*s0*u1*v0 + 4*c2*j2*n0*s0*u1*v0 + 4*c0*j2*n2*s0*u1*v0 + 
   2*f0*i2*q2*s0*u1*v0 - 2*m2*n0*q2*s0*u1*v0 + 2*d2*i2*r0*s0*u1*v0 - 
   2*l2*n2*r0*s0*u1*v0 + 2*d0*i2*r2*s0*u1*v0 - 2*l2*n0*r2*s0*u1*v0 - 
   2*c2*i2*pow(s0,2)*u1*v0 + 2*l2*m2*pow(s0,2)*u1*v0 + 
   4*c0*j2*n0*s2*u1*v0 + 2*d0*i2*r0*s2*u1*v0 - 2*l2*n0*r0*s2*u1*v0 - 
   4*c0*i2*s0*s2*u1*v0 - 2*c1*pow(n0,2)*o2*u2*v0 - 4*c0*n0*n1*o2*u2*v0 - 
   2*d1*j2*n0*r0*u2*v0 - 2*d0*j2*n1*r0*u2*v0 + 4*n0*n1*q2*r0*u2*v0 - 
   2*d0*j2*n0*r1*u2*v0 + 2*pow(n0,2)*q2*r1*u2*v0 + 4*c1*j2*n0*s0*u2*v0 + 
   4*c0*j2*n1*s0*u2*v0 + 2*d1*i2*r0*s0*u2*v0 - 2*l2*n1*r0*s0*u2*v0 + 
   2*d0*i2*r1*s0*u2*v0 - 2*l2*n0*r1*s0*u2*v0 - 2*c1*i2*pow(s0,2)*u2*v0 + 
   4*c0*j2*n0*s1*u2*v0 + 2*d0*i2*r0*s1*u2*v0 - 2*l2*n0*r0*s1*u2*v0 - 
   4*c0*i2*s0*s1*u2*v0 - 2*d0*d1*pow(j2,2)*pow(v0,2) + 
   2*d0*d1*i2*o2*pow(v0,2) - 2*d1*l2*n0*o2*pow(v0,2) - 
   2*d0*l2*n1*o2*pow(v0,2) + 2*d1*j2*n0*q2*pow(v0,2) + 
   2*d0*j2*n1*q2*pow(v0,2) - 2*n0*n1*pow(q2,2)*pow(v0,2) + 
   2*d1*j2*l2*s0*pow(v0,2) - 2*d1*i2*q2*s0*pow(v0,2) + 
   2*l2*n1*q2*s0*pow(v0,2) + 2*d0*j2*l2*s1*pow(v0,2) - 
   2*d0*i2*q2*s1*pow(v0,2) + 2*l2*n0*q2*s1*pow(v0,2) - 
   2*pow(l2,2)*s0*s1*pow(v0,2) + 2*c0*d0*k2*n0*o2*v1 - 
   2*c0*d0*j2*n0*p2*v1 + 2*c0*pow(n0,2)*p2*q2*v1 + 
   2*pow(d0,2)*j2*k2*r0*v1 - 2*pow(d0,2)*i2*p2*r0*v1 + 
   4*d0*l2*n0*p2*r0*v1 - 2*d0*k2*n0*q2*r0*v1 - 2*c0*d0*j2*k2*s0*v1 + 
   2*c0*d0*i2*p2*s0*v1 - 2*c0*l2*n0*p2*s0*v1 - 2*c0*k2*n0*q2*s0*v1 - 
   2*d0*k2*l2*r0*s0*v1 + 2*c0*k2*l2*pow(s0,2)*v1 + 
   2*d0*f0*pow(j2,2)*u0*v1 - 2*c0*g0*pow(j2,2)*u0*v1 - 
   2*d0*f0*i2*o2*u0*v1 + 2*c0*g0*i2*o2*u0*v1 + 2*f0*l2*n0*o2*u0*v1 + 
   2*d0*m2*n0*o2*u0*v1 - 2*c2*pow(n0,2)*o2*u0*v1 - 4*c0*n0*n2*o2*u0*v1 - 
   2*f0*j2*n0*q2*u0*v1 + 2*g0*j2*l2*r0*u0*v1 - 2*d2*j2*n0*r0*u0*v1 - 
   2*d0*j2*n2*r0*u0*v1 - 2*g0*i2*q2*r0*u0*v1 + 4*n0*n2*q2*r0*u0*v1 - 
   2*d0*j2*n0*r2*u0*v1 + 2*pow(n0,2)*q2*r2*u0*v1 - 2*f0*j2*l2*s0*u0*v1 - 
   2*d0*j2*m2*s0*u0*v1 + 4*c2*j2*n0*s0*u0*v1 + 4*c0*j2*n2*s0*u0*v1 + 
   2*f0*i2*q2*s0*u0*v1 - 2*m2*n0*q2*s0*u0*v1 + 2*d2*i2*r0*s0*u0*v1 - 
   2*l2*n2*r0*s0*u0*v1 + 2*d0*i2*r2*s0*u0*v1 - 2*l2*n0*r2*s0*u0*v1 - 
   2*c2*i2*pow(s0,2)*u0*v1 + 2*l2*m2*pow(s0,2)*u0*v1 + 
   4*c0*j2*n0*s2*u0*v1 + 2*d0*i2*r0*s2*u0*v1 - 2*l2*n0*r0*s2*u0*v1 - 
   4*c0*i2*s0*s2*u0*v1 - 2*c0*pow(n0,2)*o2*u2*v1 - 2*d0*j2*n0*r0*u2*v1 + 
   2*pow(n0,2)*q2*r0*u2*v1 + 4*c0*j2*n0*s0*u2*v1 + 2*d0*i2*r0*s0*u2*v1 - 
   2*l2*n0*r0*s0*u2*v1 - 2*c0*i2*pow(s0,2)*u2*v1 - 
   2*pow(d0,2)*pow(j2,2)*v0*v1 + 2*pow(d0,2)*i2*o2*v0*v1 - 
   4*d0*l2*n0*o2*v0*v1 + 4*d0*j2*n0*q2*v0*v1 - 
   2*pow(n0,2)*pow(q2,2)*v0*v1 + 4*d0*j2*l2*s0*v0*v1 - 
   4*d0*i2*q2*s0*v0*v1 + 4*l2*n0*q2*s0*v0*v1 - 
   2*pow(l2,2)*pow(s0,2)*v0*v1 - 2*c1*pow(n0,2)*o2*u0*v2 - 
   4*c0*n0*n1*o2*u0*v2 - 2*d1*j2*n0*r0*u0*v2 - 2*d0*j2*n1*r0*u0*v2 + 
   4*n0*n1*q2*r0*u0*v2 - 2*d0*j2*n0*r1*u0*v2 + 2*pow(n0,2)*q2*r1*u0*v2 + 
   4*c1*j2*n0*s0*u0*v2 + 4*c0*j2*n1*s0*u0*v2 + 2*d1*i2*r0*s0*u0*v2 - 
   2*l2*n1*r0*s0*u0*v2 + 2*d0*i2*r1*s0*u0*v2 - 2*l2*n0*r1*s0*u0*v2 - 
   2*c1*i2*pow(s0,2)*u0*v2 + 4*c0*j2*n0*s1*u0*v2 + 2*d0*i2*r0*s1*u0*v2 - 
   2*l2*n0*r0*s1*u0*v2 - 4*c0*i2*s0*s1*u0*v2 - 2*c0*pow(n0,2)*o2*u1*v2 - 
   2*d0*j2*n0*r0*u1*v2 + 2*pow(n0,2)*q2*r0*u1*v2 + 4*c0*j2*n0*s0*u1*v2 + 
   2*d0*i2*r0*s0*u1*v2 - 2*l2*n0*r0*s0*u1*v2 - 2*c0*i2*pow(s0,2)*u1*v2 - 
   4*c0*c1*k2*n0*o2*w0 - 2*pow(c0,2)*k2*n1*o2*w0 + 4*c0*c1*j2*n0*p2*w0 + 
   2*pow(c0,2)*j2*n1*p2*w0 - 2*c1*d0*j2*k2*r0*w0 - 2*c0*d1*j2*k2*r0*w0 + 
   2*c1*d0*i2*p2*r0*w0 + 2*c0*d1*i2*p2*r0*w0 - 2*c1*l2*n0*p2*r0*w0 - 
   2*c0*l2*n1*p2*r0*w0 + 4*c1*k2*n0*q2*r0*w0 + 4*c0*k2*n1*q2*r0*w0 + 
   2*d1*k2*l2*pow(r0,2)*w0 - 2*c0*d0*j2*k2*r1*w0 + 2*c0*d0*i2*p2*r1*w0 - 
   2*c0*l2*n0*p2*r1*w0 + 4*c0*k2*n0*q2*r1*w0 + 4*d0*k2*l2*r0*r1*w0 + 
   4*c0*c1*j2*k2*s0*w0 - 4*c0*c1*i2*p2*s0*w0 - 2*c1*k2*l2*r0*s0*w0 - 
   2*c0*k2*l2*r1*s0*w0 + 2*pow(c0,2)*j2*k2*s1*w0 - 
   2*pow(c0,2)*i2*p2*s1*w0 - 2*c0*k2*l2*r0*s1*w0 - 
   2*d1*e0*pow(j2,2)*u0*w0 - 2*d0*e1*pow(j2,2)*u0*w0 + 
   2*c1*f0*pow(j2,2)*u0*w0 + 2*c0*f1*pow(j2,2)*u0*w0 + 
   2*d1*e0*i2*o2*u0*w0 + 2*d0*e1*i2*o2*u0*w0 - 2*c1*f0*i2*o2*u0*w0 - 
   2*c0*f1*i2*o2*u0*w0 - 2*e1*l2*n0*o2*u0*w0 + 2*c1*m2*n0*o2*u0*w0 - 
   2*e0*l2*n1*o2*u0*w0 + 2*c0*m2*n1*o2*u0*w0 + 2*e1*j2*n0*q2*u0*w0 + 
   2*e0*j2*n1*q2*u0*w0 - 2*f1*j2*l2*r0*u0*w0 + 4*d1*j2*m2*r0*u0*w0 - 
   2*c2*j2*n1*r0*u0*w0 - 2*c1*j2*n2*r0*u0*w0 + 2*f1*i2*q2*r0*u0*w0 - 
   2*m2*n1*q2*r0*u0*w0 - 2*f0*j2*l2*r1*u0*w0 + 4*d0*j2*m2*r1*u0*w0 - 
   2*c2*j2*n0*r1*u0*w0 - 2*c0*j2*n2*r1*u0*w0 + 2*f0*i2*q2*r1*u0*w0 - 
   2*m2*n0*q2*r1*u0*w0 - 4*d2*i2*r0*r1*u0*w0 + 4*l2*n2*r0*r1*u0*w0 - 
   2*c1*j2*n0*r2*u0*w0 - 2*c0*j2*n1*r2*u0*w0 - 4*d1*i2*r0*r2*u0*w0 + 
   4*l2*n1*r0*r2*u0*w0 - 4*d0*i2*r1*r2*u0*w0 + 4*l2*n0*r1*r2*u0*w0 + 
   2*e1*j2*l2*s0*u0*w0 - 2*c1*j2*m2*s0*u0*w0 - 2*e1*i2*q2*s0*u0*w0 + 
   2*c2*i2*r1*s0*u0*w0 - 2*l2*m2*r1*s0*u0*w0 + 2*c1*i2*r2*s0*u0*w0 + 
   2*e0*j2*l2*s1*u0*w0 - 2*c0*j2*m2*s1*u0*w0 - 2*e0*i2*q2*s1*u0*w0 + 
   2*c2*i2*r0*s1*u0*w0 - 2*l2*m2*r0*s1*u0*w0 + 2*c0*i2*r2*s1*u0*w0 + 
   2*c1*i2*r0*s2*u0*w0 + 2*c0*i2*r1*s2*u0*w0 - 2*d0*e0*pow(j2,2)*u1*w0 + 
   2*c0*f0*pow(j2,2)*u1*w0 + 2*d0*e0*i2*o2*u1*w0 - 2*c0*f0*i2*o2*u1*w0 - 
   2*e0*l2*n0*o2*u1*w0 + 2*c0*m2*n0*o2*u1*w0 + 2*e0*j2*n0*q2*u1*w0 - 
   2*f0*j2*l2*r0*u1*w0 + 4*d0*j2*m2*r0*u1*w0 - 2*c2*j2*n0*r0*u1*w0 - 
   2*c0*j2*n2*r0*u1*w0 + 2*f0*i2*q2*r0*u1*w0 - 2*m2*n0*q2*r0*u1*w0 - 
   2*d2*i2*pow(r0,2)*u1*w0 + 2*l2*n2*pow(r0,2)*u1*w0 - 
   2*c0*j2*n0*r2*u1*w0 - 4*d0*i2*r0*r2*u1*w0 + 4*l2*n0*r0*r2*u1*w0 + 
   2*e0*j2*l2*s0*u1*w0 - 2*c0*j2*m2*s0*u1*w0 - 2*e0*i2*q2*s0*u1*w0 + 
   2*c2*i2*r0*s0*u1*w0 - 2*l2*m2*r0*s0*u1*w0 + 2*c0*i2*r2*s0*u1*w0 + 
   2*c0*i2*r0*s2*u1*w0 - 2*c1*j2*n0*r0*u2*w0 - 2*c0*j2*n1*r0*u2*w0 - 
   2*d1*i2*pow(r0,2)*u2*w0 + 2*l2*n1*pow(r0,2)*u2*w0 - 
   2*c0*j2*n0*r1*u2*w0 - 4*d0*i2*r0*r1*u2*w0 + 4*l2*n0*r0*r1*u2*w0 + 
   2*c1*i2*r0*s0*u2*w0 + 2*c0*i2*r1*s0*u2*w0 + 2*c0*i2*r0*s1*u2*w0 + 
   2*c1*d0*pow(j2,2)*v0*w0 + 2*c0*d1*pow(j2,2)*v0*w0 - 
   2*c1*d0*i2*o2*v0*w0 - 2*c0*d1*i2*o2*v0*w0 + 2*c1*l2*n0*o2*v0*w0 + 
   2*c0*l2*n1*o2*v0*w0 - 2*c1*j2*n0*q2*v0*w0 - 2*c0*j2*n1*q2*v0*w0 - 
   2*d1*j2*l2*r0*v0*w0 + 2*d1*i2*q2*r0*v0*w0 - 2*l2*n1*q2*r0*v0*w0 - 
   2*d0*j2*l2*r1*v0*w0 + 2*d0*i2*q2*r1*v0*w0 - 2*l2*n0*q2*r1*v0*w0 - 
   2*c1*j2*l2*s0*v0*w0 + 2*c1*i2*q2*s0*v0*w0 + 2*pow(l2,2)*r1*s0*v0*w0 - 
   2*c0*j2*l2*s1*v0*w0 + 2*c0*i2*q2*s1*v0*w0 + 2*pow(l2,2)*r0*s1*v0*w0 + 
   2*c0*d0*pow(j2,2)*v1*w0 - 2*c0*d0*i2*o2*v1*w0 + 2*c0*l2*n0*o2*v1*w0 - 
   2*c0*j2*n0*q2*v1*w0 - 2*d0*j2*l2*r0*v1*w0 + 2*d0*i2*q2*r0*v1*w0 - 
   2*l2*n0*q2*r0*v1*w0 - 2*c0*j2*l2*s0*v1*w0 + 2*c0*i2*q2*s0*v1*w0 + 
   2*pow(l2,2)*r0*s0*v1*w0 - 2*c0*c1*pow(j2,2)*pow(w0,2) + 
   2*c0*c1*i2*o2*pow(w0,2) + 2*c1*j2*l2*r0*pow(w0,2) - 
   2*c1*i2*q2*r0*pow(w0,2) + 2*c0*j2*l2*r1*pow(w0,2) - 
   2*c0*i2*q2*r1*pow(w0,2) - 2*pow(l2,2)*r0*r1*pow(w0,2) - 
   2*pow(c0,2)*k2*n0*o2*w1 + 2*pow(c0,2)*j2*n0*p2*w1 - 
   2*c0*d0*j2*k2*r0*w1 + 2*c0*d0*i2*p2*r0*w1 - 2*c0*l2*n0*p2*r0*w1 + 
   4*c0*k2*n0*q2*r0*w1 + 2*d0*k2*l2*pow(r0,2)*w1 + 
   2*pow(c0,2)*j2*k2*s0*w1 - 2*pow(c0,2)*i2*p2*s0*w1 - 
   2*c0*k2*l2*r0*s0*w1 - 2*d0*e0*pow(j2,2)*u0*w1 + 
   2*c0*f0*pow(j2,2)*u0*w1 + 2*d0*e0*i2*o2*u0*w1 - 2*c0*f0*i2*o2*u0*w1 - 
   2*e0*l2*n0*o2*u0*w1 + 2*c0*m2*n0*o2*u0*w1 + 2*e0*j2*n0*q2*u0*w1 - 
   2*f0*j2*l2*r0*u0*w1 + 4*d0*j2*m2*r0*u0*w1 - 2*c2*j2*n0*r0*u0*w1 - 
   2*c0*j2*n2*r0*u0*w1 + 2*f0*i2*q2*r0*u0*w1 - 2*m2*n0*q2*r0*u0*w1 - 
   2*d2*i2*pow(r0,2)*u0*w1 + 2*l2*n2*pow(r0,2)*u0*w1 - 
   2*c0*j2*n0*r2*u0*w1 - 4*d0*i2*r0*r2*u0*w1 + 4*l2*n0*r0*r2*u0*w1 + 
   2*e0*j2*l2*s0*u0*w1 - 2*c0*j2*m2*s0*u0*w1 - 2*e0*i2*q2*s0*u0*w1 + 
   2*c2*i2*r0*s0*u0*w1 - 2*l2*m2*r0*s0*u0*w1 + 2*c0*i2*r2*s0*u0*w1 + 
   2*c0*i2*r0*s2*u0*w1 - 2*c0*j2*n0*r0*u2*w1 - 2*d0*i2*pow(r0,2)*u2*w1 + 
   2*l2*n0*pow(r0,2)*u2*w1 + 2*c0*i2*r0*s0*u2*w1 + 
   2*c0*d0*pow(j2,2)*v0*w1 - 2*c0*d0*i2*o2*v0*w1 + 2*c0*l2*n0*o2*v0*w1 - 
   2*c0*j2*n0*q2*v0*w1 - 2*d0*j2*l2*r0*v0*w1 + 2*d0*i2*q2*r0*v0*w1 - 
   2*l2*n0*q2*r0*v0*w1 - 2*c0*j2*l2*s0*v0*w1 + 2*c0*i2*q2*s0*v0*w1 + 
   2*pow(l2,2)*r0*s0*v0*w1 - 2*pow(c0,2)*pow(j2,2)*w0*w1 + 
   2*pow(c0,2)*i2*o2*w0*w1 + 4*c0*j2*l2*r0*w0*w1 - 4*c0*i2*q2*r0*w0*w1 - 
   2*pow(l2,2)*pow(r0,2)*w0*w1 - 2*c1*j2*n0*r0*u0*w2 - 
   2*c0*j2*n1*r0*u0*w2 - 2*d1*i2*pow(r0,2)*u0*w2 + 
   2*l2*n1*pow(r0,2)*u0*w2 - 2*c0*j2*n0*r1*u0*w2 - 4*d0*i2*r0*r1*u0*w2 + 
   4*l2*n0*r0*r1*u0*w2 + 2*c1*i2*r0*s0*u0*w2 + 2*c0*i2*r1*s0*u0*w2 + 
   2*c0*i2*r0*s1*u0*w2 - 2*c0*j2*n0*r0*u1*w2 - 2*d0*i2*pow(r0,2)*u1*w2 + 
   2*l2*n0*pow(r0,2)*u1*w2 + 2*c0*i2*r0*s0*u1*w2 + 
   e1*pow(n0,2)*pow(p2,2)*z0 + 2*e0*n0*n1*pow(p2,2)*z0 + 
   2*f1*k2*n0*p2*r0*z0 + 2*f0*k2*n1*p2*r0*z0 + 
   g1*pow(k2,2)*pow(r0,2)*z0 + 2*f0*k2*n0*p2*r1*z0 + 
   2*g0*pow(k2,2)*r0*r1*z0 - 2*e1*k2*n0*p2*s0*z0 - 2*e0*k2*n1*p2*s0*z0 - 
   2*f1*pow(k2,2)*r0*s0*z0 - 2*f0*pow(k2,2)*r1*s0*z0 + 
   e1*pow(k2,2)*pow(s0,2)*z0 - 2*e0*k2*n0*p2*s1*z0 - 
   2*f0*pow(k2,2)*r0*s1*z0 + 2*e0*pow(k2,2)*s0*s1*z0 - 
   e1*pow(n0,2)*o2*t2*z0 - 2*e0*n0*n1*o2*t2*z0 - 2*f1*j2*n0*r0*t2*z0 - 
   2*f0*j2*n1*r0*t2*z0 - g1*i2*pow(r0,2)*t2*z0 + 
   2*n1*n2*pow(r0,2)*t2*z0 - 2*f0*j2*n0*r1*t2*z0 - 2*g0*i2*r0*r1*t2*z0 + 
   4*n0*n2*r0*r1*t2*z0 + 4*n0*n1*r0*r2*t2*z0 + 2*pow(n0,2)*r1*r2*t2*z0 + 
   2*e1*j2*n0*s0*t2*z0 + 2*e0*j2*n1*s0*t2*z0 + 2*f1*i2*r0*s0*t2*z0 - 
   2*m2*n1*r0*s0*t2*z0 + 2*f0*i2*r1*s0*t2*z0 - 2*m2*n0*r1*s0*t2*z0 - 
   e1*i2*pow(s0,2)*t2*z0 + 2*e0*j2*n0*s1*t2*z0 + 2*f0*i2*r0*s1*t2*z0 - 
   2*m2*n0*r0*s1*t2*z0 - 2*e0*i2*s0*s1*t2*z0 - 2*f1*k2*n0*o2*v0*z0 - 
   2*f0*k2*n1*o2*v0*z0 + 2*f1*j2*n0*p2*v0*z0 + 2*f0*j2*n1*p2*v0*z0 - 
   2*g1*j2*k2*r0*v0*z0 + 2*g1*i2*p2*r0*v0*z0 - 4*n1*n2*p2*r0*v0*z0 - 
   2*g0*j2*k2*r1*v0*z0 + 2*g0*i2*p2*r1*v0*z0 - 4*n0*n2*p2*r1*v0*z0 - 
   4*n0*n1*p2*r2*v0*z0 + 2*f1*j2*k2*s0*v0*z0 - 2*f1*i2*p2*s0*v0*z0 + 
   2*m2*n1*p2*s0*v0*z0 + 2*k2*n2*r1*s0*v0*z0 + 2*k2*n1*r2*s0*v0*z0 + 
   2*f0*j2*k2*s1*v0*z0 - 2*f0*i2*p2*s1*v0*z0 + 2*m2*n0*p2*s1*v0*z0 + 
   2*k2*n2*r0*s1*v0*z0 + 2*k2*n0*r2*s1*v0*z0 - 4*k2*m2*s0*s1*v0*z0 + 
   2*k2*n1*r0*s2*v0*z0 + 2*k2*n0*r1*s2*v0*z0 + 
   g1*pow(j2,2)*pow(v0,2)*z0 - g1*i2*o2*pow(v0,2)*z0 + 
   2*n1*n2*o2*pow(v0,2)*z0 - 2*j2*n2*s1*pow(v0,2)*z0 - 
   2*j2*n1*s2*pow(v0,2)*z0 + 2*i2*s1*s2*pow(v0,2)*z0 - 
   2*f0*k2*n0*o2*v1*z0 + 2*f0*j2*n0*p2*v1*z0 - 2*g0*j2*k2*r0*v1*z0 + 
   2*g0*i2*p2*r0*v1*z0 - 4*n0*n2*p2*r0*v1*z0 - 2*pow(n0,2)*p2*r2*v1*z0 + 
   2*f0*j2*k2*s0*v1*z0 - 2*f0*i2*p2*s0*v1*z0 + 2*m2*n0*p2*s0*v1*z0 + 
   2*k2*n2*r0*s0*v1*z0 + 2*k2*n0*r2*s0*v1*z0 - 2*k2*m2*pow(s0,2)*v1*z0 + 
   2*k2*n0*r0*s2*v1*z0 + 2*g0*pow(j2,2)*v0*v1*z0 - 2*g0*i2*o2*v0*v1*z0 + 
   4*n0*n2*o2*v0*v1*z0 - 4*j2*n2*s0*v0*v1*z0 - 4*j2*n0*s2*v0*v1*z0 + 
   4*i2*s0*s2*v0*v1*z0 - 4*n0*n1*p2*r0*v2*z0 - 2*pow(n0,2)*p2*r1*v2*z0 + 
   2*k2*n1*r0*s0*v2*z0 + 2*k2*n0*r1*s0*v2*z0 + 2*k2*n0*r0*s1*v2*z0 + 
   4*n0*n1*o2*v0*v2*z0 - 4*j2*n1*s0*v0*v2*z0 - 4*j2*n0*s1*v0*v2*z0 + 
   4*i2*s0*s1*v0*v2*z0 + 2*pow(n0,2)*o2*v1*v2*z0 - 4*j2*n0*s0*v1*v2*z0 + 
   2*i2*pow(s0,2)*v1*v2*z0 + 2*e1*k2*n0*o2*w0*z0 + 2*e0*k2*n1*o2*w0*z0 - 
   2*e1*j2*n0*p2*w0*z0 - 2*e0*j2*n1*p2*w0*z0 + 2*f1*j2*k2*r0*w0*z0 - 
   2*f1*i2*p2*r0*w0*z0 + 2*m2*n1*p2*r0*w0*z0 + 2*f0*j2*k2*r1*w0*z0 - 
   2*f0*i2*p2*r1*w0*z0 + 2*m2*n0*p2*r1*w0*z0 - 4*k2*n2*r0*r1*w0*z0 - 
   4*k2*n1*r0*r2*w0*z0 - 4*k2*n0*r1*r2*w0*z0 - 2*e1*j2*k2*s0*w0*z0 + 
   2*e1*i2*p2*s0*w0*z0 + 2*k2*m2*r1*s0*w0*z0 - 2*e0*j2*k2*s1*w0*z0 + 
   2*e0*i2*p2*s1*w0*z0 + 2*k2*m2*r0*s1*w0*z0 - 2*f1*pow(j2,2)*v0*w0*z0 + 
   2*f1*i2*o2*v0*w0*z0 - 2*m2*n1*o2*v0*w0*z0 + 2*j2*n2*r1*v0*w0*z0 + 
   2*j2*n1*r2*v0*w0*z0 + 2*j2*m2*s1*v0*w0*z0 - 2*i2*r2*s1*v0*w0*z0 - 
   2*i2*r1*s2*v0*w0*z0 - 2*f0*pow(j2,2)*v1*w0*z0 + 2*f0*i2*o2*v1*w0*z0 - 
   2*m2*n0*o2*v1*w0*z0 + 2*j2*n2*r0*v1*w0*z0 + 2*j2*n0*r2*v1*w0*z0 + 
   2*j2*m2*s0*v1*w0*z0 - 2*i2*r2*s0*v1*w0*z0 - 2*i2*r0*s2*v1*w0*z0 + 
   2*j2*n1*r0*v2*w0*z0 + 2*j2*n0*r1*v2*w0*z0 - 2*i2*r1*s0*v2*w0*z0 - 
   2*i2*r0*s1*v2*w0*z0 + e1*pow(j2,2)*pow(w0,2)*z0 - 
   e1*i2*o2*pow(w0,2)*z0 - 2*j2*m2*r1*pow(w0,2)*z0 + 
   2*i2*r1*r2*pow(w0,2)*z0 + 2*e0*k2*n0*o2*w1*z0 - 2*e0*j2*n0*p2*w1*z0 + 
   2*f0*j2*k2*r0*w1*z0 - 2*f0*i2*p2*r0*w1*z0 + 2*m2*n0*p2*r0*w1*z0 - 
   2*k2*n2*pow(r0,2)*w1*z0 - 4*k2*n0*r0*r2*w1*z0 - 2*e0*j2*k2*s0*w1*z0 + 
   2*e0*i2*p2*s0*w1*z0 + 2*k2*m2*r0*s0*w1*z0 - 2*f0*pow(j2,2)*v0*w1*z0 + 
   2*f0*i2*o2*v0*w1*z0 - 2*m2*n0*o2*v0*w1*z0 + 2*j2*n2*r0*v0*w1*z0 + 
   2*j2*n0*r2*v0*w1*z0 + 2*j2*m2*s0*v0*w1*z0 - 2*i2*r2*s0*v0*w1*z0 - 
   2*i2*r0*s2*v0*w1*z0 + 2*j2*n0*r0*v2*w1*z0 - 2*i2*r0*s0*v2*w1*z0 + 
   2*e0*pow(j2,2)*w0*w1*z0 - 2*e0*i2*o2*w0*w1*z0 - 4*j2*m2*r0*w0*w1*z0 + 
   4*i2*r0*r2*w0*w1*z0 - 2*k2*n1*pow(r0,2)*w2*z0 - 4*k2*n0*r0*r1*w2*z0 + 
   2*j2*n1*r0*v0*w2*z0 + 2*j2*n0*r1*v0*w2*z0 - 2*i2*r1*s0*v0*w2*z0 - 
   2*i2*r0*s1*v0*w2*z0 + 2*j2*n0*r0*v1*w2*z0 - 2*i2*r0*s0*v1*w2*z0 + 
   4*i2*r0*r1*w0*w2*z0 + 2*i2*pow(r0,2)*w1*w2*z0 + 
   e0*pow(n0,2)*pow(p2,2)*z1 + 2*f0*k2*n0*p2*r0*z1 + 
   g0*pow(k2,2)*pow(r0,2)*z1 - 2*e0*k2*n0*p2*s0*z1 - 
   2*f0*pow(k2,2)*r0*s0*z1 + e0*pow(k2,2)*pow(s0,2)*z1 - 
   e0*pow(n0,2)*o2*t2*z1 - 2*f0*j2*n0*r0*t2*z1 - g0*i2*pow(r0,2)*t2*z1 + 
   2*n0*n2*pow(r0,2)*t2*z1 + 2*pow(n0,2)*r0*r2*t2*z1 + 
   2*e0*j2*n0*s0*t2*z1 + 2*f0*i2*r0*s0*t2*z1 - 2*m2*n0*r0*s0*t2*z1 - 
   e0*i2*pow(s0,2)*t2*z1 - 2*f0*k2*n0*o2*v0*z1 + 2*f0*j2*n0*p2*v0*z1 - 
   2*g0*j2*k2*r0*v0*z1 + 2*g0*i2*p2*r0*v0*z1 - 4*n0*n2*p2*r0*v0*z1 - 
   2*pow(n0,2)*p2*r2*v0*z1 + 2*f0*j2*k2*s0*v0*z1 - 2*f0*i2*p2*s0*v0*z1 + 
   2*m2*n0*p2*s0*v0*z1 + 2*k2*n2*r0*s0*v0*z1 + 2*k2*n0*r2*s0*v0*z1 - 
   2*k2*m2*pow(s0,2)*v0*z1 + 2*k2*n0*r0*s2*v0*z1 + 
   g0*pow(j2,2)*pow(v0,2)*z1 - g0*i2*o2*pow(v0,2)*z1 + 
   2*n0*n2*o2*pow(v0,2)*z1 - 2*j2*n2*s0*pow(v0,2)*z1 - 
   2*j2*n0*s2*pow(v0,2)*z1 + 2*i2*s0*s2*pow(v0,2)*z1 - 
   2*pow(n0,2)*p2*r0*v2*z1 + 2*k2*n0*r0*s0*v2*z1 + 
   2*pow(n0,2)*o2*v0*v2*z1 - 4*j2*n0*s0*v0*v2*z1 + 
   2*i2*pow(s0,2)*v0*v2*z1 + 2*e0*k2*n0*o2*w0*z1 - 2*e0*j2*n0*p2*w0*z1 + 
   2*f0*j2*k2*r0*w0*z1 - 2*f0*i2*p2*r0*w0*z1 + 2*m2*n0*p2*r0*w0*z1 - 
   2*k2*n2*pow(r0,2)*w0*z1 - 4*k2*n0*r0*r2*w0*z1 - 2*e0*j2*k2*s0*w0*z1 + 
   2*e0*i2*p2*s0*w0*z1 + 2*k2*m2*r0*s0*w0*z1 - 2*f0*pow(j2,2)*v0*w0*z1 + 
   2*f0*i2*o2*v0*w0*z1 - 2*m2*n0*o2*v0*w0*z1 + 2*j2*n2*r0*v0*w0*z1 + 
   2*j2*n0*r2*v0*w0*z1 + 2*j2*m2*s0*v0*w0*z1 - 2*i2*r2*s0*v0*w0*z1 - 
   2*i2*r0*s2*v0*w0*z1 + 2*j2*n0*r0*v2*w0*z1 - 2*i2*r0*s0*v2*w0*z1 + 
   e0*pow(j2,2)*pow(w0,2)*z1 - e0*i2*o2*pow(w0,2)*z1 - 
   2*j2*m2*r0*pow(w0,2)*z1 + 2*i2*r0*r2*pow(w0,2)*z1 - 
   2*k2*n0*pow(r0,2)*w2*z1 + 2*j2*n0*r0*v0*w2*z1 - 2*i2*r0*s0*v0*w2*z1 + 
   2*i2*pow(r0,2)*w0*w2*z1 + 2*n0*n1*pow(r0,2)*t2*z2 + 
   2*pow(n0,2)*r0*r1*t2*z2 - 4*n0*n1*p2*r0*v0*z2 - 
   2*pow(n0,2)*p2*r1*v0*z2 + 2*k2*n1*r0*s0*v0*z2 + 2*k2*n0*r1*s0*v0*z2 + 
   2*k2*n0*r0*s1*v0*z2 + 2*n0*n1*o2*pow(v0,2)*z2 - 
   2*j2*n1*s0*pow(v0,2)*z2 - 2*j2*n0*s1*pow(v0,2)*z2 + 
   2*i2*s0*s1*pow(v0,2)*z2 - 2*pow(n0,2)*p2*r0*v1*z2 + 
   2*k2*n0*r0*s0*v1*z2 + 2*pow(n0,2)*o2*v0*v1*z2 - 4*j2*n0*s0*v0*v1*z2 + 
   2*i2*pow(s0,2)*v0*v1*z2 - 2*k2*n1*pow(r0,2)*w0*z2 - 
   4*k2*n0*r0*r1*w0*z2 + 2*j2*n1*r0*v0*w0*z2 + 2*j2*n0*r1*v0*w0*z2 - 
   2*i2*r1*s0*v0*w0*z2 - 2*i2*r0*s1*v0*w0*z2 + 2*j2*n0*r0*v1*w0*z2 - 
   2*i2*r0*s0*v1*w0*z2 + 2*i2*r0*r1*pow(w0,2)*z2 - 
   2*k2*n0*pow(r0,2)*w1*z2 + 2*j2*n0*r0*v0*w1*z2 - 2*i2*r0*s0*v0*w1*z2 + 
   2*i2*pow(r0,2)*w0*w1*z2
;

	k[12]  = -(pow(c1,2)*pow(n0,2)*pow(p2,2)) - 4*c0*c1*n0*n1*pow(p2,2) - 
   pow(c0,2)*pow(n1,2)*pow(p2,2) - 2*c1*d1*k2*n0*p2*r0 - 
   2*c1*d0*k2*n1*p2*r0 - 2*c0*d1*k2*n1*p2*r0 - 
   pow(d1,2)*pow(k2,2)*pow(r0,2) - 2*c1*d0*k2*n0*p2*r1 - 
   2*c0*d1*k2*n0*p2*r1 - 2*c0*d0*k2*n1*p2*r1 - 4*d0*d1*pow(k2,2)*r0*r1 - 
   pow(d0,2)*pow(k2,2)*pow(r1,2) + 2*pow(c1,2)*k2*n0*p2*s0 + 
   4*c0*c1*k2*n1*p2*s0 + 2*c1*d1*pow(k2,2)*r0*s0 + 
   2*c1*d0*pow(k2,2)*r1*s0 + 2*c0*d1*pow(k2,2)*r1*s0 - 
   pow(c1,2)*pow(k2,2)*pow(s0,2) + 4*c0*c1*k2*n0*p2*s1 + 
   2*pow(c0,2)*k2*n1*p2*s1 + 2*c1*d0*pow(k2,2)*r0*s1 + 
   2*c0*d1*pow(k2,2)*r0*s1 + 2*c0*d0*pow(k2,2)*r1*s1 - 
   4*c0*c1*pow(k2,2)*s0*s1 - pow(c0,2)*pow(k2,2)*pow(s1,2) + 
   pow(c1,2)*pow(n0,2)*o2*t2 + 4*c0*c1*n0*n1*o2*t2 + 
   pow(c0,2)*pow(n1,2)*o2*t2 + 2*c1*d1*j2*n0*r0*t2 + 
   2*c1*d0*j2*n1*r0*t2 + 2*c0*d1*j2*n1*r0*t2 - 4*c1*n0*n1*q2*r0*t2 - 
   2*c0*pow(n1,2)*q2*r0*t2 + pow(d1,2)*i2*pow(r0,2)*t2 - 
   2*d1*l2*n1*pow(r0,2)*t2 + 2*c1*d0*j2*n0*r1*t2 + 2*c0*d1*j2*n0*r1*t2 + 
   2*c0*d0*j2*n1*r1*t2 - 2*c1*pow(n0,2)*q2*r1*t2 - 4*c0*n0*n1*q2*r1*t2 + 
   4*d0*d1*i2*r0*r1*t2 - 4*d1*l2*n0*r0*r1*t2 - 4*d0*l2*n1*r0*r1*t2 + 
   pow(d0,2)*i2*pow(r1,2)*t2 - 2*d0*l2*n0*pow(r1,2)*t2 - 
   2*pow(c1,2)*j2*n0*s0*t2 - 4*c0*c1*j2*n1*s0*t2 - 2*c1*d1*i2*r0*s0*t2 + 
   2*c1*l2*n1*r0*s0*t2 - 2*c1*d0*i2*r1*s0*t2 - 2*c0*d1*i2*r1*s0*t2 + 
   2*c1*l2*n0*r1*s0*t2 + 2*c0*l2*n1*r1*s0*t2 + 
   pow(c1,2)*i2*pow(s0,2)*t2 - 4*c0*c1*j2*n0*s1*t2 - 
   2*pow(c0,2)*j2*n1*s1*t2 - 2*c1*d0*i2*r0*s1*t2 - 2*c0*d1*i2*r0*s1*t2 + 
   2*c1*l2*n0*r0*s1*t2 + 2*c0*l2*n1*r0*s1*t2 - 2*c0*d0*i2*r1*s1*t2 + 
   2*c0*l2*n0*r1*s1*t2 + 4*c0*c1*i2*s0*s1*t2 + 
   pow(c0,2)*i2*pow(s1,2)*t2 - 2*d1*e1*k2*n0*o2*u0 + 
   2*c1*f1*k2*n0*o2*u0 - 2*d1*e0*k2*n1*o2*u0 - 2*d0*e1*k2*n1*o2*u0 + 
   2*c1*f0*k2*n1*o2*u0 + 2*c0*f1*k2*n1*o2*u0 + 2*d1*e1*j2*n0*p2*u0 - 
   2*c1*f1*j2*n0*p2*u0 + 2*d1*e0*j2*n1*p2*u0 + 2*d0*e1*j2*n1*p2*u0 - 
   2*c1*f0*j2*n1*p2*u0 - 2*c0*f1*j2*n1*p2*u0 - 4*e1*n0*n1*p2*q2*u0 - 
   2*e0*pow(n1,2)*p2*q2*u0 - 2*d1*f1*j2*k2*r0*u0 + 2*c1*g1*j2*k2*r0*u0 + 
   2*d1*f1*i2*p2*r0*u0 - 2*c1*g1*i2*p2*r0*u0 - 2*f1*l2*n1*p2*r0*u0 - 
   2*d1*m2*n1*p2*r0*u0 + 2*c2*pow(n1,2)*p2*r0*u0 + 4*c1*n1*n2*p2*r0*u0 - 
   2*f1*k2*n1*q2*r0*u0 - 2*d1*f0*j2*k2*r1*u0 - 2*d0*f1*j2*k2*r1*u0 + 
   2*c1*g0*j2*k2*r1*u0 + 2*c0*g1*j2*k2*r1*u0 + 2*d1*f0*i2*p2*r1*u0 + 
   2*d0*f1*i2*p2*r1*u0 - 2*c1*g0*i2*p2*r1*u0 - 2*c0*g1*i2*p2*r1*u0 - 
   2*f1*l2*n0*p2*r1*u0 - 2*d1*m2*n0*p2*r1*u0 - 2*f0*l2*n1*p2*r1*u0 - 
   2*d0*m2*n1*p2*r1*u0 + 4*c2*n0*n1*p2*r1*u0 + 4*c1*n0*n2*p2*r1*u0 + 
   4*c0*n1*n2*p2*r1*u0 - 2*f1*k2*n0*q2*r1*u0 - 2*f0*k2*n1*q2*r1*u0 - 
   4*g1*k2*l2*r0*r1*u0 + 4*d2*k2*n1*r0*r1*u0 + 4*d1*k2*n2*r0*r1*u0 - 
   2*g0*k2*l2*pow(r1,2)*u0 + 2*d2*k2*n0*pow(r1,2)*u0 + 
   2*d0*k2*n2*pow(r1,2)*u0 + 4*c1*n0*n1*p2*r2*u0 + 
   2*c0*pow(n1,2)*p2*r2*u0 + 4*d1*k2*n1*r0*r2*u0 + 4*d1*k2*n0*r1*r2*u0 + 
   4*d0*k2*n1*r1*r2*u0 + 2*d1*e1*j2*k2*s0*u0 - 2*c1*f1*j2*k2*s0*u0 - 
   2*d1*e1*i2*p2*s0*u0 + 2*c1*f1*i2*p2*s0*u0 + 2*e1*l2*n1*p2*s0*u0 - 
   2*c1*m2*n1*p2*s0*u0 + 2*e1*k2*n1*q2*s0*u0 + 4*f1*k2*l2*r1*s0*u0 - 
   2*d1*k2*m2*r1*s0*u0 - 2*c2*k2*n1*r1*s0*u0 - 2*c1*k2*n2*r1*s0*u0 - 
   2*c1*k2*n1*r2*s0*u0 + 2*d1*e0*j2*k2*s1*u0 + 2*d0*e1*j2*k2*s1*u0 - 
   2*c1*f0*j2*k2*s1*u0 - 2*c0*f1*j2*k2*s1*u0 - 2*d1*e0*i2*p2*s1*u0 - 
   2*d0*e1*i2*p2*s1*u0 + 2*c1*f0*i2*p2*s1*u0 + 2*c0*f1*i2*p2*s1*u0 + 
   2*e1*l2*n0*p2*s1*u0 - 2*c1*m2*n0*p2*s1*u0 + 2*e0*l2*n1*p2*s1*u0 - 
   2*c0*m2*n1*p2*s1*u0 + 2*e1*k2*n0*q2*s1*u0 + 2*e0*k2*n1*q2*s1*u0 + 
   4*f1*k2*l2*r0*s1*u0 - 2*d1*k2*m2*r0*s1*u0 - 2*c2*k2*n1*r0*s1*u0 - 
   2*c1*k2*n2*r0*s1*u0 + 4*f0*k2*l2*r1*s1*u0 - 2*d0*k2*m2*r1*s1*u0 - 
   2*c2*k2*n0*r1*s1*u0 - 2*c0*k2*n2*r1*s1*u0 - 2*c1*k2*n0*r2*s1*u0 - 
   2*c0*k2*n1*r2*s1*u0 - 4*e1*k2*l2*s0*s1*u0 + 4*c1*k2*m2*s0*s1*u0 - 
   2*e0*k2*l2*pow(s1,2)*u0 + 2*c0*k2*m2*pow(s1,2)*u0 - 
   2*c1*k2*n1*r0*s2*u0 - 2*c1*k2*n0*r1*s2*u0 - 2*c0*k2*n1*r1*s2*u0 - 
   pow(f1,2)*pow(j2,2)*pow(u0,2) + e1*g1*pow(j2,2)*pow(u0,2) + 
   pow(f1,2)*i2*o2*pow(u0,2) - e1*g1*i2*o2*pow(u0,2) - 
   2*f1*m2*n1*o2*pow(u0,2) + e2*pow(n1,2)*o2*pow(u0,2) + 
   2*e1*n1*n2*o2*pow(u0,2) - 2*g1*j2*m2*r1*pow(u0,2) + 
   2*f2*j2*n1*r1*pow(u0,2) + 2*f1*j2*n2*r1*pow(u0,2) + 
   g2*i2*pow(r1,2)*pow(u0,2) - pow(n2,2)*pow(r1,2)*pow(u0,2) + 
   2*f1*j2*n1*r2*pow(u0,2) + 2*g1*i2*r1*r2*pow(u0,2) - 
   4*n1*n2*r1*r2*pow(u0,2) - pow(n1,2)*pow(r2,2)*pow(u0,2) + 
   2*f1*j2*m2*s1*pow(u0,2) - 2*e2*j2*n1*s1*pow(u0,2) - 
   2*e1*j2*n2*s1*pow(u0,2) - 2*f2*i2*r1*s1*pow(u0,2) + 
   2*m2*n2*r1*s1*pow(u0,2) - 2*f1*i2*r2*s1*pow(u0,2) + 
   2*m2*n1*r2*s1*pow(u0,2) + e2*i2*pow(s1,2)*pow(u0,2) - 
   pow(m2,2)*pow(s1,2)*pow(u0,2) - 2*e1*j2*n1*s2*pow(u0,2) - 
   2*f1*i2*r1*s2*pow(u0,2) + 2*m2*n1*r1*s2*pow(u0,2) + 
   2*e1*i2*s1*s2*pow(u0,2) - 2*d1*e0*k2*n0*o2*u1 - 2*d0*e1*k2*n0*o2*u1 + 
   2*c1*f0*k2*n0*o2*u1 + 2*c0*f1*k2*n0*o2*u1 - 2*d0*e0*k2*n1*o2*u1 + 
   2*c0*f0*k2*n1*o2*u1 + 2*d1*e0*j2*n0*p2*u1 + 2*d0*e1*j2*n0*p2*u1 - 
   2*c1*f0*j2*n0*p2*u1 - 2*c0*f1*j2*n0*p2*u1 + 2*d0*e0*j2*n1*p2*u1 - 
   2*c0*f0*j2*n1*p2*u1 - 2*e1*pow(n0,2)*p2*q2*u1 - 4*e0*n0*n1*p2*q2*u1 - 
   2*d1*f0*j2*k2*r0*u1 - 2*d0*f1*j2*k2*r0*u1 + 2*c1*g0*j2*k2*r0*u1 + 
   2*c0*g1*j2*k2*r0*u1 + 2*d1*f0*i2*p2*r0*u1 + 2*d0*f1*i2*p2*r0*u1 - 
   2*c1*g0*i2*p2*r0*u1 - 2*c0*g1*i2*p2*r0*u1 - 2*f1*l2*n0*p2*r0*u1 - 
   2*d1*m2*n0*p2*r0*u1 - 2*f0*l2*n1*p2*r0*u1 - 2*d0*m2*n1*p2*r0*u1 + 
   4*c2*n0*n1*p2*r0*u1 + 4*c1*n0*n2*p2*r0*u1 + 4*c0*n1*n2*p2*r0*u1 - 
   2*f1*k2*n0*q2*r0*u1 - 2*f0*k2*n1*q2*r0*u1 - 2*g1*k2*l2*pow(r0,2)*u1 + 
   2*d2*k2*n1*pow(r0,2)*u1 + 2*d1*k2*n2*pow(r0,2)*u1 - 
   2*d0*f0*j2*k2*r1*u1 + 2*c0*g0*j2*k2*r1*u1 + 2*d0*f0*i2*p2*r1*u1 - 
   2*c0*g0*i2*p2*r1*u1 - 2*f0*l2*n0*p2*r1*u1 - 2*d0*m2*n0*p2*r1*u1 + 
   2*c2*pow(n0,2)*p2*r1*u1 + 4*c0*n0*n2*p2*r1*u1 - 2*f0*k2*n0*q2*r1*u1 - 
   4*g0*k2*l2*r0*r1*u1 + 4*d2*k2*n0*r0*r1*u1 + 4*d0*k2*n2*r0*r1*u1 + 
   2*c1*pow(n0,2)*p2*r2*u1 + 4*c0*n0*n1*p2*r2*u1 + 4*d1*k2*n0*r0*r2*u1 + 
   4*d0*k2*n1*r0*r2*u1 + 4*d0*k2*n0*r1*r2*u1 + 2*d1*e0*j2*k2*s0*u1 + 
   2*d0*e1*j2*k2*s0*u1 - 2*c1*f0*j2*k2*s0*u1 - 2*c0*f1*j2*k2*s0*u1 - 
   2*d1*e0*i2*p2*s0*u1 - 2*d0*e1*i2*p2*s0*u1 + 2*c1*f0*i2*p2*s0*u1 + 
   2*c0*f1*i2*p2*s0*u1 + 2*e1*l2*n0*p2*s0*u1 - 2*c1*m2*n0*p2*s0*u1 + 
   2*e0*l2*n1*p2*s0*u1 - 2*c0*m2*n1*p2*s0*u1 + 2*e1*k2*n0*q2*s0*u1 + 
   2*e0*k2*n1*q2*s0*u1 + 4*f1*k2*l2*r0*s0*u1 - 2*d1*k2*m2*r0*s0*u1 - 
   2*c2*k2*n1*r0*s0*u1 - 2*c1*k2*n2*r0*s0*u1 + 4*f0*k2*l2*r1*s0*u1 - 
   2*d0*k2*m2*r1*s0*u1 - 2*c2*k2*n0*r1*s0*u1 - 2*c0*k2*n2*r1*s0*u1 - 
   2*c1*k2*n0*r2*s0*u1 - 2*c0*k2*n1*r2*s0*u1 - 2*e1*k2*l2*pow(s0,2)*u1 + 
   2*c1*k2*m2*pow(s0,2)*u1 + 2*d0*e0*j2*k2*s1*u1 - 2*c0*f0*j2*k2*s1*u1 - 
   2*d0*e0*i2*p2*s1*u1 + 2*c0*f0*i2*p2*s1*u1 + 2*e0*l2*n0*p2*s1*u1 - 
   2*c0*m2*n0*p2*s1*u1 + 2*e0*k2*n0*q2*s1*u1 + 4*f0*k2*l2*r0*s1*u1 - 
   2*d0*k2*m2*r0*s1*u1 - 2*c2*k2*n0*r0*s1*u1 - 2*c0*k2*n2*r0*s1*u1 - 
   2*c0*k2*n0*r2*s1*u1 - 4*e0*k2*l2*s0*s1*u1 + 4*c0*k2*m2*s0*s1*u1 - 
   2*c1*k2*n0*r0*s2*u1 - 2*c0*k2*n1*r0*s2*u1 - 2*c0*k2*n0*r1*s2*u1 - 
   4*f0*f1*pow(j2,2)*u0*u1 + 2*e1*g0*pow(j2,2)*u0*u1 + 
   2*e0*g1*pow(j2,2)*u0*u1 + 4*f0*f1*i2*o2*u0*u1 - 2*e1*g0*i2*o2*u0*u1 - 
   2*e0*g1*i2*o2*u0*u1 - 4*f1*m2*n0*o2*u0*u1 - 4*f0*m2*n1*o2*u0*u1 + 
   4*e2*n0*n1*o2*u0*u1 + 4*e1*n0*n2*o2*u0*u1 + 4*e0*n1*n2*o2*u0*u1 - 
   4*g1*j2*m2*r0*u0*u1 + 4*f2*j2*n1*r0*u0*u1 + 4*f1*j2*n2*r0*u0*u1 - 
   4*g0*j2*m2*r1*u0*u1 + 4*f2*j2*n0*r1*u0*u1 + 4*f0*j2*n2*r1*u0*u1 + 
   4*g2*i2*r0*r1*u0*u1 - 4*pow(n2,2)*r0*r1*u0*u1 + 4*f1*j2*n0*r2*u0*u1 + 
   4*f0*j2*n1*r2*u0*u1 + 4*g1*i2*r0*r2*u0*u1 - 8*n1*n2*r0*r2*u0*u1 + 
   4*g0*i2*r1*r2*u0*u1 - 8*n0*n2*r1*r2*u0*u1 - 4*n0*n1*pow(r2,2)*u0*u1 + 
   4*f1*j2*m2*s0*u0*u1 - 4*e2*j2*n1*s0*u0*u1 - 4*e1*j2*n2*s0*u0*u1 - 
   4*f2*i2*r1*s0*u0*u1 + 4*m2*n2*r1*s0*u0*u1 - 4*f1*i2*r2*s0*u0*u1 + 
   4*m2*n1*r2*s0*u0*u1 + 4*f0*j2*m2*s1*u0*u1 - 4*e2*j2*n0*s1*u0*u1 - 
   4*e0*j2*n2*s1*u0*u1 - 4*f2*i2*r0*s1*u0*u1 + 4*m2*n2*r0*s1*u0*u1 - 
   4*f0*i2*r2*s1*u0*u1 + 4*m2*n0*r2*s1*u0*u1 + 4*e2*i2*s0*s1*u0*u1 - 
   4*pow(m2,2)*s0*s1*u0*u1 - 4*e1*j2*n0*s2*u0*u1 - 4*e0*j2*n1*s2*u0*u1 - 
   4*f1*i2*r0*s2*u0*u1 + 4*m2*n1*r0*s2*u0*u1 - 4*f0*i2*r1*s2*u0*u1 + 
   4*m2*n0*r1*s2*u0*u1 + 4*e1*i2*s0*s2*u0*u1 + 4*e0*i2*s1*s2*u0*u1 - 
   pow(f0,2)*pow(j2,2)*pow(u1,2) + e0*g0*pow(j2,2)*pow(u1,2) + 
   pow(f0,2)*i2*o2*pow(u1,2) - e0*g0*i2*o2*pow(u1,2) - 
   2*f0*m2*n0*o2*pow(u1,2) + e2*pow(n0,2)*o2*pow(u1,2) + 
   2*e0*n0*n2*o2*pow(u1,2) - 2*g0*j2*m2*r0*pow(u1,2) + 
   2*f2*j2*n0*r0*pow(u1,2) + 2*f0*j2*n2*r0*pow(u1,2) + 
   g2*i2*pow(r0,2)*pow(u1,2) - pow(n2,2)*pow(r0,2)*pow(u1,2) + 
   2*f0*j2*n0*r2*pow(u1,2) + 2*g0*i2*r0*r2*pow(u1,2) - 
   4*n0*n2*r0*r2*pow(u1,2) - pow(n0,2)*pow(r2,2)*pow(u1,2) + 
   2*f0*j2*m2*s0*pow(u1,2) - 2*e2*j2*n0*s0*pow(u1,2) - 
   2*e0*j2*n2*s0*pow(u1,2) - 2*f2*i2*r0*s0*pow(u1,2) + 
   2*m2*n2*r0*s0*pow(u1,2) - 2*f0*i2*r2*s0*pow(u1,2) + 
   2*m2*n0*r2*s0*pow(u1,2) + e2*i2*pow(s0,2)*pow(u1,2) - 
   pow(m2,2)*pow(s0,2)*pow(u1,2) - 2*e0*j2*n0*s2*pow(u1,2) - 
   2*f0*i2*r0*s2*pow(u1,2) + 2*m2*n0*r0*s2*pow(u1,2) + 
   2*e0*i2*s0*s2*pow(u1,2) + 4*c1*n0*n1*p2*r0*u2 + 
   2*c0*pow(n1,2)*p2*r0*u2 + 2*d1*k2*n1*pow(r0,2)*u2 + 
   2*c1*pow(n0,2)*p2*r1*u2 + 4*c0*n0*n1*p2*r1*u2 + 4*d1*k2*n0*r0*r1*u2 + 
   4*d0*k2*n1*r0*r1*u2 + 2*d0*k2*n0*pow(r1,2)*u2 - 2*c1*k2*n1*r0*s0*u2 - 
   2*c1*k2*n0*r1*s0*u2 - 2*c0*k2*n1*r1*s0*u2 - 2*c1*k2*n0*r0*s1*u2 - 
   2*c0*k2*n1*r0*s1*u2 - 2*c0*k2*n0*r1*s1*u2 + 4*e1*n0*n1*o2*u0*u2 + 
   2*e0*pow(n1,2)*o2*u0*u2 + 4*f1*j2*n1*r0*u0*u2 + 4*f1*j2*n0*r1*u0*u2 + 
   4*f0*j2*n1*r1*u0*u2 + 4*g1*i2*r0*r1*u0*u2 - 8*n1*n2*r0*r1*u0*u2 + 
   2*g0*i2*pow(r1,2)*u0*u2 - 4*n0*n2*pow(r1,2)*u0*u2 - 
   4*pow(n1,2)*r0*r2*u0*u2 - 8*n0*n1*r1*r2*u0*u2 - 4*e1*j2*n1*s0*u0*u2 - 
   4*f1*i2*r1*s0*u0*u2 + 4*m2*n1*r1*s0*u0*u2 - 4*e1*j2*n0*s1*u0*u2 - 
   4*e0*j2*n1*s1*u0*u2 - 4*f1*i2*r0*s1*u0*u2 + 4*m2*n1*r0*s1*u0*u2 - 
   4*f0*i2*r1*s1*u0*u2 + 4*m2*n0*r1*s1*u0*u2 + 4*e1*i2*s0*s1*u0*u2 + 
   2*e0*i2*pow(s1,2)*u0*u2 + 2*e1*pow(n0,2)*o2*u1*u2 + 
   4*e0*n0*n1*o2*u1*u2 + 4*f1*j2*n0*r0*u1*u2 + 4*f0*j2*n1*r0*u1*u2 + 
   2*g1*i2*pow(r0,2)*u1*u2 - 4*n1*n2*pow(r0,2)*u1*u2 + 
   4*f0*j2*n0*r1*u1*u2 + 4*g0*i2*r0*r1*u1*u2 - 8*n0*n2*r0*r1*u1*u2 - 
   8*n0*n1*r0*r2*u1*u2 - 4*pow(n0,2)*r1*r2*u1*u2 - 4*e1*j2*n0*s0*u1*u2 - 
   4*e0*j2*n1*s0*u1*u2 - 4*f1*i2*r0*s0*u1*u2 + 4*m2*n1*r0*s0*u1*u2 - 
   4*f0*i2*r1*s0*u1*u2 + 4*m2*n0*r1*s0*u1*u2 + 2*e1*i2*pow(s0,2)*u1*u2 - 
   4*e0*j2*n0*s1*u1*u2 - 4*f0*i2*r0*s1*u1*u2 + 4*m2*n0*r0*s1*u1*u2 + 
   4*e0*i2*s0*s1*u1*u2 - pow(n1,2)*pow(r0,2)*pow(u2,2) - 
   4*n0*n1*r0*r1*pow(u2,2) - pow(n0,2)*pow(r1,2)*pow(u2,2) + 
   2*c1*d1*k2*n0*o2*v0 + 2*c1*d0*k2*n1*o2*v0 + 2*c0*d1*k2*n1*o2*v0 - 
   2*c1*d1*j2*n0*p2*v0 - 2*c1*d0*j2*n1*p2*v0 - 2*c0*d1*j2*n1*p2*v0 + 
   4*c1*n0*n1*p2*q2*v0 + 2*c0*pow(n1,2)*p2*q2*v0 + 
   2*pow(d1,2)*j2*k2*r0*v0 - 2*pow(d1,2)*i2*p2*r0*v0 + 
   4*d1*l2*n1*p2*r0*v0 - 2*d1*k2*n1*q2*r0*v0 + 4*d0*d1*j2*k2*r1*v0 - 
   4*d0*d1*i2*p2*r1*v0 + 4*d1*l2*n0*p2*r1*v0 + 4*d0*l2*n1*p2*r1*v0 - 
   2*d1*k2*n0*q2*r1*v0 - 2*d0*k2*n1*q2*r1*v0 - 2*c1*d1*j2*k2*s0*v0 + 
   2*c1*d1*i2*p2*s0*v0 - 2*c1*l2*n1*p2*s0*v0 - 2*c1*k2*n1*q2*s0*v0 - 
   2*d1*k2*l2*r1*s0*v0 - 2*c1*d0*j2*k2*s1*v0 - 2*c0*d1*j2*k2*s1*v0 + 
   2*c1*d0*i2*p2*s1*v0 + 2*c0*d1*i2*p2*s1*v0 - 2*c1*l2*n0*p2*s1*v0 - 
   2*c0*l2*n1*p2*s1*v0 - 2*c1*k2*n0*q2*s1*v0 - 2*c0*k2*n1*q2*s1*v0 - 
   2*d1*k2*l2*r0*s1*v0 - 2*d0*k2*l2*r1*s1*v0 + 4*c1*k2*l2*s0*s1*v0 + 
   2*c0*k2*l2*pow(s1,2)*v0 + 2*d1*f1*pow(j2,2)*u0*v0 - 
   2*c1*g1*pow(j2,2)*u0*v0 - 2*d1*f1*i2*o2*u0*v0 + 2*c1*g1*i2*o2*u0*v0 + 
   2*f1*l2*n1*o2*u0*v0 + 2*d1*m2*n1*o2*u0*v0 - 2*c2*pow(n1,2)*o2*u0*v0 - 
   4*c1*n1*n2*o2*u0*v0 - 2*f1*j2*n1*q2*u0*v0 + 2*g1*j2*l2*r1*u0*v0 - 
   2*d2*j2*n1*r1*u0*v0 - 2*d1*j2*n2*r1*u0*v0 - 2*g1*i2*q2*r1*u0*v0 + 
   4*n1*n2*q2*r1*u0*v0 - 2*d1*j2*n1*r2*u0*v0 + 2*pow(n1,2)*q2*r2*u0*v0 - 
   2*f1*j2*l2*s1*u0*v0 - 2*d1*j2*m2*s1*u0*v0 + 4*c2*j2*n1*s1*u0*v0 + 
   4*c1*j2*n2*s1*u0*v0 + 2*f1*i2*q2*s1*u0*v0 - 2*m2*n1*q2*s1*u0*v0 + 
   2*d2*i2*r1*s1*u0*v0 - 2*l2*n2*r1*s1*u0*v0 + 2*d1*i2*r2*s1*u0*v0 - 
   2*l2*n1*r2*s1*u0*v0 - 2*c2*i2*pow(s1,2)*u0*v0 + 
   2*l2*m2*pow(s1,2)*u0*v0 + 4*c1*j2*n1*s2*u0*v0 + 2*d1*i2*r1*s2*u0*v0 - 
   2*l2*n1*r1*s2*u0*v0 - 4*c1*i2*s1*s2*u0*v0 + 2*d1*f0*pow(j2,2)*u1*v0 + 
   2*d0*f1*pow(j2,2)*u1*v0 - 2*c1*g0*pow(j2,2)*u1*v0 - 
   2*c0*g1*pow(j2,2)*u1*v0 - 2*d1*f0*i2*o2*u1*v0 - 2*d0*f1*i2*o2*u1*v0 + 
   2*c1*g0*i2*o2*u1*v0 + 2*c0*g1*i2*o2*u1*v0 + 2*f1*l2*n0*o2*u1*v0 + 
   2*d1*m2*n0*o2*u1*v0 + 2*f0*l2*n1*o2*u1*v0 + 2*d0*m2*n1*o2*u1*v0 - 
   4*c2*n0*n1*o2*u1*v0 - 4*c1*n0*n2*o2*u1*v0 - 4*c0*n1*n2*o2*u1*v0 - 
   2*f1*j2*n0*q2*u1*v0 - 2*f0*j2*n1*q2*u1*v0 + 2*g1*j2*l2*r0*u1*v0 - 
   2*d2*j2*n1*r0*u1*v0 - 2*d1*j2*n2*r0*u1*v0 - 2*g1*i2*q2*r0*u1*v0 + 
   4*n1*n2*q2*r0*u1*v0 + 2*g0*j2*l2*r1*u1*v0 - 2*d2*j2*n0*r1*u1*v0 - 
   2*d0*j2*n2*r1*u1*v0 - 2*g0*i2*q2*r1*u1*v0 + 4*n0*n2*q2*r1*u1*v0 - 
   2*d1*j2*n0*r2*u1*v0 - 2*d0*j2*n1*r2*u1*v0 + 4*n0*n1*q2*r2*u1*v0 - 
   2*f1*j2*l2*s0*u1*v0 - 2*d1*j2*m2*s0*u1*v0 + 4*c2*j2*n1*s0*u1*v0 + 
   4*c1*j2*n2*s0*u1*v0 + 2*f1*i2*q2*s0*u1*v0 - 2*m2*n1*q2*s0*u1*v0 + 
   2*d2*i2*r1*s0*u1*v0 - 2*l2*n2*r1*s0*u1*v0 + 2*d1*i2*r2*s0*u1*v0 - 
   2*l2*n1*r2*s0*u1*v0 - 2*f0*j2*l2*s1*u1*v0 - 2*d0*j2*m2*s1*u1*v0 + 
   4*c2*j2*n0*s1*u1*v0 + 4*c0*j2*n2*s1*u1*v0 + 2*f0*i2*q2*s1*u1*v0 - 
   2*m2*n0*q2*s1*u1*v0 + 2*d2*i2*r0*s1*u1*v0 - 2*l2*n2*r0*s1*u1*v0 + 
   2*d0*i2*r2*s1*u1*v0 - 2*l2*n0*r2*s1*u1*v0 - 4*c2*i2*s0*s1*u1*v0 + 
   4*l2*m2*s0*s1*u1*v0 + 4*c1*j2*n0*s2*u1*v0 + 4*c0*j2*n1*s2*u1*v0 + 
   2*d1*i2*r0*s2*u1*v0 - 2*l2*n1*r0*s2*u1*v0 + 2*d0*i2*r1*s2*u1*v0 - 
   2*l2*n0*r1*s2*u1*v0 - 4*c1*i2*s0*s2*u1*v0 - 4*c0*i2*s1*s2*u1*v0 - 
   4*c1*n0*n1*o2*u2*v0 - 2*c0*pow(n1,2)*o2*u2*v0 - 2*d1*j2*n1*r0*u2*v0 + 
   2*pow(n1,2)*q2*r0*u2*v0 - 2*d1*j2*n0*r1*u2*v0 - 2*d0*j2*n1*r1*u2*v0 + 
   4*n0*n1*q2*r1*u2*v0 + 4*c1*j2*n1*s0*u2*v0 + 2*d1*i2*r1*s0*u2*v0 - 
   2*l2*n1*r1*s0*u2*v0 + 4*c1*j2*n0*s1*u2*v0 + 4*c0*j2*n1*s1*u2*v0 + 
   2*d1*i2*r0*s1*u2*v0 - 2*l2*n1*r0*s1*u2*v0 + 2*d0*i2*r1*s1*u2*v0 - 
   2*l2*n0*r1*s1*u2*v0 - 4*c1*i2*s0*s1*u2*v0 - 2*c0*i2*pow(s1,2)*u2*v0 - 
   pow(d1,2)*pow(j2,2)*pow(v0,2) + pow(d1,2)*i2*o2*pow(v0,2) - 
   2*d1*l2*n1*o2*pow(v0,2) + 2*d1*j2*n1*q2*pow(v0,2) - 
   pow(n1,2)*pow(q2,2)*pow(v0,2) + 2*d1*j2*l2*s1*pow(v0,2) - 
   2*d1*i2*q2*s1*pow(v0,2) + 2*l2*n1*q2*s1*pow(v0,2) - 
   pow(l2,2)*pow(s1,2)*pow(v0,2) + 2*c1*d0*k2*n0*o2*v1 + 
   2*c0*d1*k2*n0*o2*v1 + 2*c0*d0*k2*n1*o2*v1 - 2*c1*d0*j2*n0*p2*v1 - 
   2*c0*d1*j2*n0*p2*v1 - 2*c0*d0*j2*n1*p2*v1 + 2*c1*pow(n0,2)*p2*q2*v1 + 
   4*c0*n0*n1*p2*q2*v1 + 4*d0*d1*j2*k2*r0*v1 - 4*d0*d1*i2*p2*r0*v1 + 
   4*d1*l2*n0*p2*r0*v1 + 4*d0*l2*n1*p2*r0*v1 - 2*d1*k2*n0*q2*r0*v1 - 
   2*d0*k2*n1*q2*r0*v1 + 2*pow(d0,2)*j2*k2*r1*v1 - 
   2*pow(d0,2)*i2*p2*r1*v1 + 4*d0*l2*n0*p2*r1*v1 - 2*d0*k2*n0*q2*r1*v1 - 
   2*c1*d0*j2*k2*s0*v1 - 2*c0*d1*j2*k2*s0*v1 + 2*c1*d0*i2*p2*s0*v1 + 
   2*c0*d1*i2*p2*s0*v1 - 2*c1*l2*n0*p2*s0*v1 - 2*c0*l2*n1*p2*s0*v1 - 
   2*c1*k2*n0*q2*s0*v1 - 2*c0*k2*n1*q2*s0*v1 - 2*d1*k2*l2*r0*s0*v1 - 
   2*d0*k2*l2*r1*s0*v1 + 2*c1*k2*l2*pow(s0,2)*v1 - 2*c0*d0*j2*k2*s1*v1 + 
   2*c0*d0*i2*p2*s1*v1 - 2*c0*l2*n0*p2*s1*v1 - 2*c0*k2*n0*q2*s1*v1 - 
   2*d0*k2*l2*r0*s1*v1 + 4*c0*k2*l2*s0*s1*v1 + 2*d1*f0*pow(j2,2)*u0*v1 + 
   2*d0*f1*pow(j2,2)*u0*v1 - 2*c1*g0*pow(j2,2)*u0*v1 - 
   2*c0*g1*pow(j2,2)*u0*v1 - 2*d1*f0*i2*o2*u0*v1 - 2*d0*f1*i2*o2*u0*v1 + 
   2*c1*g0*i2*o2*u0*v1 + 2*c0*g1*i2*o2*u0*v1 + 2*f1*l2*n0*o2*u0*v1 + 
   2*d1*m2*n0*o2*u0*v1 + 2*f0*l2*n1*o2*u0*v1 + 2*d0*m2*n1*o2*u0*v1 - 
   4*c2*n0*n1*o2*u0*v1 - 4*c1*n0*n2*o2*u0*v1 - 4*c0*n1*n2*o2*u0*v1 - 
   2*f1*j2*n0*q2*u0*v1 - 2*f0*j2*n1*q2*u0*v1 + 2*g1*j2*l2*r0*u0*v1 - 
   2*d2*j2*n1*r0*u0*v1 - 2*d1*j2*n2*r0*u0*v1 - 2*g1*i2*q2*r0*u0*v1 + 
   4*n1*n2*q2*r0*u0*v1 + 2*g0*j2*l2*r1*u0*v1 - 2*d2*j2*n0*r1*u0*v1 - 
   2*d0*j2*n2*r1*u0*v1 - 2*g0*i2*q2*r1*u0*v1 + 4*n0*n2*q2*r1*u0*v1 - 
   2*d1*j2*n0*r2*u0*v1 - 2*d0*j2*n1*r2*u0*v1 + 4*n0*n1*q2*r2*u0*v1 - 
   2*f1*j2*l2*s0*u0*v1 - 2*d1*j2*m2*s0*u0*v1 + 4*c2*j2*n1*s0*u0*v1 + 
   4*c1*j2*n2*s0*u0*v1 + 2*f1*i2*q2*s0*u0*v1 - 2*m2*n1*q2*s0*u0*v1 + 
   2*d2*i2*r1*s0*u0*v1 - 2*l2*n2*r1*s0*u0*v1 + 2*d1*i2*r2*s0*u0*v1 - 
   2*l2*n1*r2*s0*u0*v1 - 2*f0*j2*l2*s1*u0*v1 - 2*d0*j2*m2*s1*u0*v1 + 
   4*c2*j2*n0*s1*u0*v1 + 4*c0*j2*n2*s1*u0*v1 + 2*f0*i2*q2*s1*u0*v1 - 
   2*m2*n0*q2*s1*u0*v1 + 2*d2*i2*r0*s1*u0*v1 - 2*l2*n2*r0*s1*u0*v1 + 
   2*d0*i2*r2*s1*u0*v1 - 2*l2*n0*r2*s1*u0*v1 - 4*c2*i2*s0*s1*u0*v1 + 
   4*l2*m2*s0*s1*u0*v1 + 4*c1*j2*n0*s2*u0*v1 + 4*c0*j2*n1*s2*u0*v1 + 
   2*d1*i2*r0*s2*u0*v1 - 2*l2*n1*r0*s2*u0*v1 + 2*d0*i2*r1*s2*u0*v1 - 
   2*l2*n0*r1*s2*u0*v1 - 4*c1*i2*s0*s2*u0*v1 - 4*c0*i2*s1*s2*u0*v1 + 
   2*d0*f0*pow(j2,2)*u1*v1 - 2*c0*g0*pow(j2,2)*u1*v1 - 
   2*d0*f0*i2*o2*u1*v1 + 2*c0*g0*i2*o2*u1*v1 + 2*f0*l2*n0*o2*u1*v1 + 
   2*d0*m2*n0*o2*u1*v1 - 2*c2*pow(n0,2)*o2*u1*v1 - 4*c0*n0*n2*o2*u1*v1 - 
   2*f0*j2*n0*q2*u1*v1 + 2*g0*j2*l2*r0*u1*v1 - 2*d2*j2*n0*r0*u1*v1 - 
   2*d0*j2*n2*r0*u1*v1 - 2*g0*i2*q2*r0*u1*v1 + 4*n0*n2*q2*r0*u1*v1 - 
   2*d0*j2*n0*r2*u1*v1 + 2*pow(n0,2)*q2*r2*u1*v1 - 2*f0*j2*l2*s0*u1*v1 - 
   2*d0*j2*m2*s0*u1*v1 + 4*c2*j2*n0*s0*u1*v1 + 4*c0*j2*n2*s0*u1*v1 + 
   2*f0*i2*q2*s0*u1*v1 - 2*m2*n0*q2*s0*u1*v1 + 2*d2*i2*r0*s0*u1*v1 - 
   2*l2*n2*r0*s0*u1*v1 + 2*d0*i2*r2*s0*u1*v1 - 2*l2*n0*r2*s0*u1*v1 - 
   2*c2*i2*pow(s0,2)*u1*v1 + 2*l2*m2*pow(s0,2)*u1*v1 + 
   4*c0*j2*n0*s2*u1*v1 + 2*d0*i2*r0*s2*u1*v1 - 2*l2*n0*r0*s2*u1*v1 - 
   4*c0*i2*s0*s2*u1*v1 - 2*c1*pow(n0,2)*o2*u2*v1 - 4*c0*n0*n1*o2*u2*v1 - 
   2*d1*j2*n0*r0*u2*v1 - 2*d0*j2*n1*r0*u2*v1 + 4*n0*n1*q2*r0*u2*v1 - 
   2*d0*j2*n0*r1*u2*v1 + 2*pow(n0,2)*q2*r1*u2*v1 + 4*c1*j2*n0*s0*u2*v1 + 
   4*c0*j2*n1*s0*u2*v1 + 2*d1*i2*r0*s0*u2*v1 - 2*l2*n1*r0*s0*u2*v1 + 
   2*d0*i2*r1*s0*u2*v1 - 2*l2*n0*r1*s0*u2*v1 - 2*c1*i2*pow(s0,2)*u2*v1 + 
   4*c0*j2*n0*s1*u2*v1 + 2*d0*i2*r0*s1*u2*v1 - 2*l2*n0*r0*s1*u2*v1 - 
   4*c0*i2*s0*s1*u2*v1 - 4*d0*d1*pow(j2,2)*v0*v1 + 4*d0*d1*i2*o2*v0*v1 - 
   4*d1*l2*n0*o2*v0*v1 - 4*d0*l2*n1*o2*v0*v1 + 4*d1*j2*n0*q2*v0*v1 + 
   4*d0*j2*n1*q2*v0*v1 - 4*n0*n1*pow(q2,2)*v0*v1 + 4*d1*j2*l2*s0*v0*v1 - 
   4*d1*i2*q2*s0*v0*v1 + 4*l2*n1*q2*s0*v0*v1 + 4*d0*j2*l2*s1*v0*v1 - 
   4*d0*i2*q2*s1*v0*v1 + 4*l2*n0*q2*s1*v0*v1 - 4*pow(l2,2)*s0*s1*v0*v1 - 
   pow(d0,2)*pow(j2,2)*pow(v1,2) + pow(d0,2)*i2*o2*pow(v1,2) - 
   2*d0*l2*n0*o2*pow(v1,2) + 2*d0*j2*n0*q2*pow(v1,2) - 
   pow(n0,2)*pow(q2,2)*pow(v1,2) + 2*d0*j2*l2*s0*pow(v1,2) - 
   2*d0*i2*q2*s0*pow(v1,2) + 2*l2*n0*q2*s0*pow(v1,2) - 
   pow(l2,2)*pow(s0,2)*pow(v1,2) - 4*c1*n0*n1*o2*u0*v2 - 
   2*c0*pow(n1,2)*o2*u0*v2 - 2*d1*j2*n1*r0*u0*v2 + 
   2*pow(n1,2)*q2*r0*u0*v2 - 2*d1*j2*n0*r1*u0*v2 - 2*d0*j2*n1*r1*u0*v2 + 
   4*n0*n1*q2*r1*u0*v2 + 4*c1*j2*n1*s0*u0*v2 + 2*d1*i2*r1*s0*u0*v2 - 
   2*l2*n1*r1*s0*u0*v2 + 4*c1*j2*n0*s1*u0*v2 + 4*c0*j2*n1*s1*u0*v2 + 
   2*d1*i2*r0*s1*u0*v2 - 2*l2*n1*r0*s1*u0*v2 + 2*d0*i2*r1*s1*u0*v2 - 
   2*l2*n0*r1*s1*u0*v2 - 4*c1*i2*s0*s1*u0*v2 - 2*c0*i2*pow(s1,2)*u0*v2 - 
   2*c1*pow(n0,2)*o2*u1*v2 - 4*c0*n0*n1*o2*u1*v2 - 2*d1*j2*n0*r0*u1*v2 - 
   2*d0*j2*n1*r0*u1*v2 + 4*n0*n1*q2*r0*u1*v2 - 2*d0*j2*n0*r1*u1*v2 + 
   2*pow(n0,2)*q2*r1*u1*v2 + 4*c1*j2*n0*s0*u1*v2 + 4*c0*j2*n1*s0*u1*v2 + 
   2*d1*i2*r0*s0*u1*v2 - 2*l2*n1*r0*s0*u1*v2 + 2*d0*i2*r1*s0*u1*v2 - 
   2*l2*n0*r1*s0*u1*v2 - 2*c1*i2*pow(s0,2)*u1*v2 + 4*c0*j2*n0*s1*u1*v2 + 
   2*d0*i2*r0*s1*u1*v2 - 2*l2*n0*r0*s1*u1*v2 - 4*c0*i2*s0*s1*u1*v2 - 
   2*pow(c1,2)*k2*n0*o2*w0 - 4*c0*c1*k2*n1*o2*w0 + 
   2*pow(c1,2)*j2*n0*p2*w0 + 4*c0*c1*j2*n1*p2*w0 - 2*c1*d1*j2*k2*r0*w0 + 
   2*c1*d1*i2*p2*r0*w0 - 2*c1*l2*n1*p2*r0*w0 + 4*c1*k2*n1*q2*r0*w0 - 
   2*c1*d0*j2*k2*r1*w0 - 2*c0*d1*j2*k2*r1*w0 + 2*c1*d0*i2*p2*r1*w0 + 
   2*c0*d1*i2*p2*r1*w0 - 2*c1*l2*n0*p2*r1*w0 - 2*c0*l2*n1*p2*r1*w0 + 
   4*c1*k2*n0*q2*r1*w0 + 4*c0*k2*n1*q2*r1*w0 + 4*d1*k2*l2*r0*r1*w0 + 
   2*d0*k2*l2*pow(r1,2)*w0 + 2*pow(c1,2)*j2*k2*s0*w0 - 
   2*pow(c1,2)*i2*p2*s0*w0 - 2*c1*k2*l2*r1*s0*w0 + 4*c0*c1*j2*k2*s1*w0 - 
   4*c0*c1*i2*p2*s1*w0 - 2*c1*k2*l2*r0*s1*w0 - 2*c0*k2*l2*r1*s1*w0 - 
   2*d1*e1*pow(j2,2)*u0*w0 + 2*c1*f1*pow(j2,2)*u0*w0 + 
   2*d1*e1*i2*o2*u0*w0 - 2*c1*f1*i2*o2*u0*w0 - 2*e1*l2*n1*o2*u0*w0 + 
   2*c1*m2*n1*o2*u0*w0 + 2*e1*j2*n1*q2*u0*w0 - 2*f1*j2*l2*r1*u0*w0 + 
   4*d1*j2*m2*r1*u0*w0 - 2*c2*j2*n1*r1*u0*w0 - 2*c1*j2*n2*r1*u0*w0 + 
   2*f1*i2*q2*r1*u0*w0 - 2*m2*n1*q2*r1*u0*w0 - 2*d2*i2*pow(r1,2)*u0*w0 + 
   2*l2*n2*pow(r1,2)*u0*w0 - 2*c1*j2*n1*r2*u0*w0 - 4*d1*i2*r1*r2*u0*w0 + 
   4*l2*n1*r1*r2*u0*w0 + 2*e1*j2*l2*s1*u0*w0 - 2*c1*j2*m2*s1*u0*w0 - 
   2*e1*i2*q2*s1*u0*w0 + 2*c2*i2*r1*s1*u0*w0 - 2*l2*m2*r1*s1*u0*w0 + 
   2*c1*i2*r2*s1*u0*w0 + 2*c1*i2*r1*s2*u0*w0 - 2*d1*e0*pow(j2,2)*u1*w0 - 
   2*d0*e1*pow(j2,2)*u1*w0 + 2*c1*f0*pow(j2,2)*u1*w0 + 
   2*c0*f1*pow(j2,2)*u1*w0 + 2*d1*e0*i2*o2*u1*w0 + 2*d0*e1*i2*o2*u1*w0 - 
   2*c1*f0*i2*o2*u1*w0 - 2*c0*f1*i2*o2*u1*w0 - 2*e1*l2*n0*o2*u1*w0 + 
   2*c1*m2*n0*o2*u1*w0 - 2*e0*l2*n1*o2*u1*w0 + 2*c0*m2*n1*o2*u1*w0 + 
   2*e1*j2*n0*q2*u1*w0 + 2*e0*j2*n1*q2*u1*w0 - 2*f1*j2*l2*r0*u1*w0 + 
   4*d1*j2*m2*r0*u1*w0 - 2*c2*j2*n1*r0*u1*w0 - 2*c1*j2*n2*r0*u1*w0 + 
   2*f1*i2*q2*r0*u1*w0 - 2*m2*n1*q2*r0*u1*w0 - 2*f0*j2*l2*r1*u1*w0 + 
   4*d0*j2*m2*r1*u1*w0 - 2*c2*j2*n0*r1*u1*w0 - 2*c0*j2*n2*r1*u1*w0 + 
   2*f0*i2*q2*r1*u1*w0 - 2*m2*n0*q2*r1*u1*w0 - 4*d2*i2*r0*r1*u1*w0 + 
   4*l2*n2*r0*r1*u1*w0 - 2*c1*j2*n0*r2*u1*w0 - 2*c0*j2*n1*r2*u1*w0 - 
   4*d1*i2*r0*r2*u1*w0 + 4*l2*n1*r0*r2*u1*w0 - 4*d0*i2*r1*r2*u1*w0 + 
   4*l2*n0*r1*r2*u1*w0 + 2*e1*j2*l2*s0*u1*w0 - 2*c1*j2*m2*s0*u1*w0 - 
   2*e1*i2*q2*s0*u1*w0 + 2*c2*i2*r1*s0*u1*w0 - 2*l2*m2*r1*s0*u1*w0 + 
   2*c1*i2*r2*s0*u1*w0 + 2*e0*j2*l2*s1*u1*w0 - 2*c0*j2*m2*s1*u1*w0 - 
   2*e0*i2*q2*s1*u1*w0 + 2*c2*i2*r0*s1*u1*w0 - 2*l2*m2*r0*s1*u1*w0 + 
   2*c0*i2*r2*s1*u1*w0 + 2*c1*i2*r0*s2*u1*w0 + 2*c0*i2*r1*s2*u1*w0 - 
   2*c1*j2*n1*r0*u2*w0 - 2*c1*j2*n0*r1*u2*w0 - 2*c0*j2*n1*r1*u2*w0 - 
   4*d1*i2*r0*r1*u2*w0 + 4*l2*n1*r0*r1*u2*w0 - 2*d0*i2*pow(r1,2)*u2*w0 + 
   2*l2*n0*pow(r1,2)*u2*w0 + 2*c1*i2*r1*s0*u2*w0 + 2*c1*i2*r0*s1*u2*w0 + 
   2*c0*i2*r1*s1*u2*w0 + 2*c1*d1*pow(j2,2)*v0*w0 - 2*c1*d1*i2*o2*v0*w0 + 
   2*c1*l2*n1*o2*v0*w0 - 2*c1*j2*n1*q2*v0*w0 - 2*d1*j2*l2*r1*v0*w0 + 
   2*d1*i2*q2*r1*v0*w0 - 2*l2*n1*q2*r1*v0*w0 - 2*c1*j2*l2*s1*v0*w0 + 
   2*c1*i2*q2*s1*v0*w0 + 2*pow(l2,2)*r1*s1*v0*w0 + 
   2*c1*d0*pow(j2,2)*v1*w0 + 2*c0*d1*pow(j2,2)*v1*w0 - 
   2*c1*d0*i2*o2*v1*w0 - 2*c0*d1*i2*o2*v1*w0 + 2*c1*l2*n0*o2*v1*w0 + 
   2*c0*l2*n1*o2*v1*w0 - 2*c1*j2*n0*q2*v1*w0 - 2*c0*j2*n1*q2*v1*w0 - 
   2*d1*j2*l2*r0*v1*w0 + 2*d1*i2*q2*r0*v1*w0 - 2*l2*n1*q2*r0*v1*w0 - 
   2*d0*j2*l2*r1*v1*w0 + 2*d0*i2*q2*r1*v1*w0 - 2*l2*n0*q2*r1*v1*w0 - 
   2*c1*j2*l2*s0*v1*w0 + 2*c1*i2*q2*s0*v1*w0 + 2*pow(l2,2)*r1*s0*v1*w0 - 
   2*c0*j2*l2*s1*v1*w0 + 2*c0*i2*q2*s1*v1*w0 + 2*pow(l2,2)*r0*s1*v1*w0 - 
   pow(c1,2)*pow(j2,2)*pow(w0,2) + pow(c1,2)*i2*o2*pow(w0,2) + 
   2*c1*j2*l2*r1*pow(w0,2) - 2*c1*i2*q2*r1*pow(w0,2) - 
   pow(l2,2)*pow(r1,2)*pow(w0,2) - 4*c0*c1*k2*n0*o2*w1 - 
   2*pow(c0,2)*k2*n1*o2*w1 + 4*c0*c1*j2*n0*p2*w1 + 
   2*pow(c0,2)*j2*n1*p2*w1 - 2*c1*d0*j2*k2*r0*w1 - 2*c0*d1*j2*k2*r0*w1 + 
   2*c1*d0*i2*p2*r0*w1 + 2*c0*d1*i2*p2*r0*w1 - 2*c1*l2*n0*p2*r0*w1 - 
   2*c0*l2*n1*p2*r0*w1 + 4*c1*k2*n0*q2*r0*w1 + 4*c0*k2*n1*q2*r0*w1 + 
   2*d1*k2*l2*pow(r0,2)*w1 - 2*c0*d0*j2*k2*r1*w1 + 2*c0*d0*i2*p2*r1*w1 - 
   2*c0*l2*n0*p2*r1*w1 + 4*c0*k2*n0*q2*r1*w1 + 4*d0*k2*l2*r0*r1*w1 + 
   4*c0*c1*j2*k2*s0*w1 - 4*c0*c1*i2*p2*s0*w1 - 2*c1*k2*l2*r0*s0*w1 - 
   2*c0*k2*l2*r1*s0*w1 + 2*pow(c0,2)*j2*k2*s1*w1 - 
   2*pow(c0,2)*i2*p2*s1*w1 - 2*c0*k2*l2*r0*s1*w1 - 
   2*d1*e0*pow(j2,2)*u0*w1 - 2*d0*e1*pow(j2,2)*u0*w1 + 
   2*c1*f0*pow(j2,2)*u0*w1 + 2*c0*f1*pow(j2,2)*u0*w1 + 
   2*d1*e0*i2*o2*u0*w1 + 2*d0*e1*i2*o2*u0*w1 - 2*c1*f0*i2*o2*u0*w1 - 
   2*c0*f1*i2*o2*u0*w1 - 2*e1*l2*n0*o2*u0*w1 + 2*c1*m2*n0*o2*u0*w1 - 
   2*e0*l2*n1*o2*u0*w1 + 2*c0*m2*n1*o2*u0*w1 + 2*e1*j2*n0*q2*u0*w1 + 
   2*e0*j2*n1*q2*u0*w1 - 2*f1*j2*l2*r0*u0*w1 + 4*d1*j2*m2*r0*u0*w1 - 
   2*c2*j2*n1*r0*u0*w1 - 2*c1*j2*n2*r0*u0*w1 + 2*f1*i2*q2*r0*u0*w1 - 
   2*m2*n1*q2*r0*u0*w1 - 2*f0*j2*l2*r1*u0*w1 + 4*d0*j2*m2*r1*u0*w1 - 
   2*c2*j2*n0*r1*u0*w1 - 2*c0*j2*n2*r1*u0*w1 + 2*f0*i2*q2*r1*u0*w1 - 
   2*m2*n0*q2*r1*u0*w1 - 4*d2*i2*r0*r1*u0*w1 + 4*l2*n2*r0*r1*u0*w1 - 
   2*c1*j2*n0*r2*u0*w1 - 2*c0*j2*n1*r2*u0*w1 - 4*d1*i2*r0*r2*u0*w1 + 
   4*l2*n1*r0*r2*u0*w1 - 4*d0*i2*r1*r2*u0*w1 + 4*l2*n0*r1*r2*u0*w1 + 
   2*e1*j2*l2*s0*u0*w1 - 2*c1*j2*m2*s0*u0*w1 - 2*e1*i2*q2*s0*u0*w1 + 
   2*c2*i2*r1*s0*u0*w1 - 2*l2*m2*r1*s0*u0*w1 + 2*c1*i2*r2*s0*u0*w1 + 
   2*e0*j2*l2*s1*u0*w1 - 2*c0*j2*m2*s1*u0*w1 - 2*e0*i2*q2*s1*u0*w1 + 
   2*c2*i2*r0*s1*u0*w1 - 2*l2*m2*r0*s1*u0*w1 + 2*c0*i2*r2*s1*u0*w1 + 
   2*c1*i2*r0*s2*u0*w1 + 2*c0*i2*r1*s2*u0*w1 - 2*d0*e0*pow(j2,2)*u1*w1 + 
   2*c0*f0*pow(j2,2)*u1*w1 + 2*d0*e0*i2*o2*u1*w1 - 2*c0*f0*i2*o2*u1*w1 - 
   2*e0*l2*n0*o2*u1*w1 + 2*c0*m2*n0*o2*u1*w1 + 2*e0*j2*n0*q2*u1*w1 - 
   2*f0*j2*l2*r0*u1*w1 + 4*d0*j2*m2*r0*u1*w1 - 2*c2*j2*n0*r0*u1*w1 - 
   2*c0*j2*n2*r0*u1*w1 + 2*f0*i2*q2*r0*u1*w1 - 2*m2*n0*q2*r0*u1*w1 - 
   2*d2*i2*pow(r0,2)*u1*w1 + 2*l2*n2*pow(r0,2)*u1*w1 - 
   2*c0*j2*n0*r2*u1*w1 - 4*d0*i2*r0*r2*u1*w1 + 4*l2*n0*r0*r2*u1*w1 + 
   2*e0*j2*l2*s0*u1*w1 - 2*c0*j2*m2*s0*u1*w1 - 2*e0*i2*q2*s0*u1*w1 + 
   2*c2*i2*r0*s0*u1*w1 - 2*l2*m2*r0*s0*u1*w1 + 2*c0*i2*r2*s0*u1*w1 + 
   2*c0*i2*r0*s2*u1*w1 - 2*c1*j2*n0*r0*u2*w1 - 2*c0*j2*n1*r0*u2*w1 - 
   2*d1*i2*pow(r0,2)*u2*w1 + 2*l2*n1*pow(r0,2)*u2*w1 - 
   2*c0*j2*n0*r1*u2*w1 - 4*d0*i2*r0*r1*u2*w1 + 4*l2*n0*r0*r1*u2*w1 + 
   2*c1*i2*r0*s0*u2*w1 + 2*c0*i2*r1*s0*u2*w1 + 2*c0*i2*r0*s1*u2*w1 + 
   2*c1*d0*pow(j2,2)*v0*w1 + 2*c0*d1*pow(j2,2)*v0*w1 - 
   2*c1*d0*i2*o2*v0*w1 - 2*c0*d1*i2*o2*v0*w1 + 2*c1*l2*n0*o2*v0*w1 + 
   2*c0*l2*n1*o2*v0*w1 - 2*c1*j2*n0*q2*v0*w1 - 2*c0*j2*n1*q2*v0*w1 - 
   2*d1*j2*l2*r0*v0*w1 + 2*d1*i2*q2*r0*v0*w1 - 2*l2*n1*q2*r0*v0*w1 - 
   2*d0*j2*l2*r1*v0*w1 + 2*d0*i2*q2*r1*v0*w1 - 2*l2*n0*q2*r1*v0*w1 - 
   2*c1*j2*l2*s0*v0*w1 + 2*c1*i2*q2*s0*v0*w1 + 2*pow(l2,2)*r1*s0*v0*w1 - 
   2*c0*j2*l2*s1*v0*w1 + 2*c0*i2*q2*s1*v0*w1 + 2*pow(l2,2)*r0*s1*v0*w1 + 
   2*c0*d0*pow(j2,2)*v1*w1 - 2*c0*d0*i2*o2*v1*w1 + 2*c0*l2*n0*o2*v1*w1 - 
   2*c0*j2*n0*q2*v1*w1 - 2*d0*j2*l2*r0*v1*w1 + 2*d0*i2*q2*r0*v1*w1 - 
   2*l2*n0*q2*r0*v1*w1 - 2*c0*j2*l2*s0*v1*w1 + 2*c0*i2*q2*s0*v1*w1 + 
   2*pow(l2,2)*r0*s0*v1*w1 - 4*c0*c1*pow(j2,2)*w0*w1 + 
   4*c0*c1*i2*o2*w0*w1 + 4*c1*j2*l2*r0*w0*w1 - 4*c1*i2*q2*r0*w0*w1 + 
   4*c0*j2*l2*r1*w0*w1 - 4*c0*i2*q2*r1*w0*w1 - 4*pow(l2,2)*r0*r1*w0*w1 - 
   pow(c0,2)*pow(j2,2)*pow(w1,2) + pow(c0,2)*i2*o2*pow(w1,2) + 
   2*c0*j2*l2*r0*pow(w1,2) - 2*c0*i2*q2*r0*pow(w1,2) - 
   pow(l2,2)*pow(r0,2)*pow(w1,2) - 2*c1*j2*n1*r0*u0*w2 - 
   2*c1*j2*n0*r1*u0*w2 - 2*c0*j2*n1*r1*u0*w2 - 4*d1*i2*r0*r1*u0*w2 + 
   4*l2*n1*r0*r1*u0*w2 - 2*d0*i2*pow(r1,2)*u0*w2 + 
   2*l2*n0*pow(r1,2)*u0*w2 + 2*c1*i2*r1*s0*u0*w2 + 2*c1*i2*r0*s1*u0*w2 + 
   2*c0*i2*r1*s1*u0*w2 - 2*c1*j2*n0*r0*u1*w2 - 2*c0*j2*n1*r0*u1*w2 - 
   2*d1*i2*pow(r0,2)*u1*w2 + 2*l2*n1*pow(r0,2)*u1*w2 - 
   2*c0*j2*n0*r1*u1*w2 - 4*d0*i2*r0*r1*u1*w2 + 4*l2*n0*r0*r1*u1*w2 + 
   2*c1*i2*r0*s0*u1*w2 + 2*c0*i2*r1*s0*u1*w2 + 2*c0*i2*r0*s1*u1*w2 + 
   2*e1*n0*n1*pow(p2,2)*z0 + e0*pow(n1,2)*pow(p2,2)*z0 + 
   2*f1*k2*n1*p2*r0*z0 + 2*f1*k2*n0*p2*r1*z0 + 2*f0*k2*n1*p2*r1*z0 + 
   2*g1*pow(k2,2)*r0*r1*z0 + g0*pow(k2,2)*pow(r1,2)*z0 - 
   2*e1*k2*n1*p2*s0*z0 - 2*f1*pow(k2,2)*r1*s0*z0 - 2*e1*k2*n0*p2*s1*z0 - 
   2*e0*k2*n1*p2*s1*z0 - 2*f1*pow(k2,2)*r0*s1*z0 - 
   2*f0*pow(k2,2)*r1*s1*z0 + 2*e1*pow(k2,2)*s0*s1*z0 + 
   e0*pow(k2,2)*pow(s1,2)*z0 - 2*e1*n0*n1*o2*t2*z0 - 
   e0*pow(n1,2)*o2*t2*z0 - 2*f1*j2*n1*r0*t2*z0 - 2*f1*j2*n0*r1*t2*z0 - 
   2*f0*j2*n1*r1*t2*z0 - 2*g1*i2*r0*r1*t2*z0 + 4*n1*n2*r0*r1*t2*z0 - 
   g0*i2*pow(r1,2)*t2*z0 + 2*n0*n2*pow(r1,2)*t2*z0 + 
   2*pow(n1,2)*r0*r2*t2*z0 + 4*n0*n1*r1*r2*t2*z0 + 2*e1*j2*n1*s0*t2*z0 + 
   2*f1*i2*r1*s0*t2*z0 - 2*m2*n1*r1*s0*t2*z0 + 2*e1*j2*n0*s1*t2*z0 + 
   2*e0*j2*n1*s1*t2*z0 + 2*f1*i2*r0*s1*t2*z0 - 2*m2*n1*r0*s1*t2*z0 + 
   2*f0*i2*r1*s1*t2*z0 - 2*m2*n0*r1*s1*t2*z0 - 2*e1*i2*s0*s1*t2*z0 - 
   e0*i2*pow(s1,2)*t2*z0 - 2*f1*k2*n1*o2*v0*z0 + 2*f1*j2*n1*p2*v0*z0 - 
   2*g1*j2*k2*r1*v0*z0 + 2*g1*i2*p2*r1*v0*z0 - 4*n1*n2*p2*r1*v0*z0 - 
   2*pow(n1,2)*p2*r2*v0*z0 + 2*f1*j2*k2*s1*v0*z0 - 2*f1*i2*p2*s1*v0*z0 + 
   2*m2*n1*p2*s1*v0*z0 + 2*k2*n2*r1*s1*v0*z0 + 2*k2*n1*r2*s1*v0*z0 - 
   2*k2*m2*pow(s1,2)*v0*z0 + 2*k2*n1*r1*s2*v0*z0 - 2*f1*k2*n0*o2*v1*z0 - 
   2*f0*k2*n1*o2*v1*z0 + 2*f1*j2*n0*p2*v1*z0 + 2*f0*j2*n1*p2*v1*z0 - 
   2*g1*j2*k2*r0*v1*z0 + 2*g1*i2*p2*r0*v1*z0 - 4*n1*n2*p2*r0*v1*z0 - 
   2*g0*j2*k2*r1*v1*z0 + 2*g0*i2*p2*r1*v1*z0 - 4*n0*n2*p2*r1*v1*z0 - 
   4*n0*n1*p2*r2*v1*z0 + 2*f1*j2*k2*s0*v1*z0 - 2*f1*i2*p2*s0*v1*z0 + 
   2*m2*n1*p2*s0*v1*z0 + 2*k2*n2*r1*s0*v1*z0 + 2*k2*n1*r2*s0*v1*z0 + 
   2*f0*j2*k2*s1*v1*z0 - 2*f0*i2*p2*s1*v1*z0 + 2*m2*n0*p2*s1*v1*z0 + 
   2*k2*n2*r0*s1*v1*z0 + 2*k2*n0*r2*s1*v1*z0 - 4*k2*m2*s0*s1*v1*z0 + 
   2*k2*n1*r0*s2*v1*z0 + 2*k2*n0*r1*s2*v1*z0 + 2*g1*pow(j2,2)*v0*v1*z0 - 
   2*g1*i2*o2*v0*v1*z0 + 4*n1*n2*o2*v0*v1*z0 - 4*j2*n2*s1*v0*v1*z0 - 
   4*j2*n1*s2*v0*v1*z0 + 4*i2*s1*s2*v0*v1*z0 + 
   g0*pow(j2,2)*pow(v1,2)*z0 - g0*i2*o2*pow(v1,2)*z0 + 
   2*n0*n2*o2*pow(v1,2)*z0 - 2*j2*n2*s0*pow(v1,2)*z0 - 
   2*j2*n0*s2*pow(v1,2)*z0 + 2*i2*s0*s2*pow(v1,2)*z0 - 
   2*pow(n1,2)*p2*r0*v2*z0 - 4*n0*n1*p2*r1*v2*z0 + 2*k2*n1*r1*s0*v2*z0 + 
   2*k2*n1*r0*s1*v2*z0 + 2*k2*n0*r1*s1*v2*z0 + 2*pow(n1,2)*o2*v0*v2*z0 - 
   4*j2*n1*s1*v0*v2*z0 + 2*i2*pow(s1,2)*v0*v2*z0 + 4*n0*n1*o2*v1*v2*z0 - 
   4*j2*n1*s0*v1*v2*z0 - 4*j2*n0*s1*v1*v2*z0 + 4*i2*s0*s1*v1*v2*z0 + 
   2*e1*k2*n1*o2*w0*z0 - 2*e1*j2*n1*p2*w0*z0 + 2*f1*j2*k2*r1*w0*z0 - 
   2*f1*i2*p2*r1*w0*z0 + 2*m2*n1*p2*r1*w0*z0 - 2*k2*n2*pow(r1,2)*w0*z0 - 
   4*k2*n1*r1*r2*w0*z0 - 2*e1*j2*k2*s1*w0*z0 + 2*e1*i2*p2*s1*w0*z0 + 
   2*k2*m2*r1*s1*w0*z0 - 2*f1*pow(j2,2)*v1*w0*z0 + 2*f1*i2*o2*v1*w0*z0 - 
   2*m2*n1*o2*v1*w0*z0 + 2*j2*n2*r1*v1*w0*z0 + 2*j2*n1*r2*v1*w0*z0 + 
   2*j2*m2*s1*v1*w0*z0 - 2*i2*r2*s1*v1*w0*z0 - 2*i2*r1*s2*v1*w0*z0 + 
   2*j2*n1*r1*v2*w0*z0 - 2*i2*r1*s1*v2*w0*z0 + 2*e1*k2*n0*o2*w1*z0 + 
   2*e0*k2*n1*o2*w1*z0 - 2*e1*j2*n0*p2*w1*z0 - 2*e0*j2*n1*p2*w1*z0 + 
   2*f1*j2*k2*r0*w1*z0 - 2*f1*i2*p2*r0*w1*z0 + 2*m2*n1*p2*r0*w1*z0 + 
   2*f0*j2*k2*r1*w1*z0 - 2*f0*i2*p2*r1*w1*z0 + 2*m2*n0*p2*r1*w1*z0 - 
   4*k2*n2*r0*r1*w1*z0 - 4*k2*n1*r0*r2*w1*z0 - 4*k2*n0*r1*r2*w1*z0 - 
   2*e1*j2*k2*s0*w1*z0 + 2*e1*i2*p2*s0*w1*z0 + 2*k2*m2*r1*s0*w1*z0 - 
   2*e0*j2*k2*s1*w1*z0 + 2*e0*i2*p2*s1*w1*z0 + 2*k2*m2*r0*s1*w1*z0 - 
   2*f1*pow(j2,2)*v0*w1*z0 + 2*f1*i2*o2*v0*w1*z0 - 2*m2*n1*o2*v0*w1*z0 + 
   2*j2*n2*r1*v0*w1*z0 + 2*j2*n1*r2*v0*w1*z0 + 2*j2*m2*s1*v0*w1*z0 - 
   2*i2*r2*s1*v0*w1*z0 - 2*i2*r1*s2*v0*w1*z0 - 2*f0*pow(j2,2)*v1*w1*z0 + 
   2*f0*i2*o2*v1*w1*z0 - 2*m2*n0*o2*v1*w1*z0 + 2*j2*n2*r0*v1*w1*z0 + 
   2*j2*n0*r2*v1*w1*z0 + 2*j2*m2*s0*v1*w1*z0 - 2*i2*r2*s0*v1*w1*z0 - 
   2*i2*r0*s2*v1*w1*z0 + 2*j2*n1*r0*v2*w1*z0 + 2*j2*n0*r1*v2*w1*z0 - 
   2*i2*r1*s0*v2*w1*z0 - 2*i2*r0*s1*v2*w1*z0 + 2*e1*pow(j2,2)*w0*w1*z0 - 
   2*e1*i2*o2*w0*w1*z0 - 4*j2*m2*r1*w0*w1*z0 + 4*i2*r1*r2*w0*w1*z0 + 
   e0*pow(j2,2)*pow(w1,2)*z0 - e0*i2*o2*pow(w1,2)*z0 - 
   2*j2*m2*r0*pow(w1,2)*z0 + 2*i2*r0*r2*pow(w1,2)*z0 - 
   4*k2*n1*r0*r1*w2*z0 - 2*k2*n0*pow(r1,2)*w2*z0 + 2*j2*n1*r1*v0*w2*z0 - 
   2*i2*r1*s1*v0*w2*z0 + 2*j2*n1*r0*v1*w2*z0 + 2*j2*n0*r1*v1*w2*z0 - 
   2*i2*r1*s0*v1*w2*z0 - 2*i2*r0*s1*v1*w2*z0 + 2*i2*pow(r1,2)*w0*w2*z0 + 
   4*i2*r0*r1*w1*w2*z0 + e1*pow(n0,2)*pow(p2,2)*z1 + 
   2*e0*n0*n1*pow(p2,2)*z1 + 2*f1*k2*n0*p2*r0*z1 + 2*f0*k2*n1*p2*r0*z1 + 
   g1*pow(k2,2)*pow(r0,2)*z1 + 2*f0*k2*n0*p2*r1*z1 + 
   2*g0*pow(k2,2)*r0*r1*z1 - 2*e1*k2*n0*p2*s0*z1 - 2*e0*k2*n1*p2*s0*z1 - 
   2*f1*pow(k2,2)*r0*s0*z1 - 2*f0*pow(k2,2)*r1*s0*z1 + 
   e1*pow(k2,2)*pow(s0,2)*z1 - 2*e0*k2*n0*p2*s1*z1 - 
   2*f0*pow(k2,2)*r0*s1*z1 + 2*e0*pow(k2,2)*s0*s1*z1 - 
   e1*pow(n0,2)*o2*t2*z1 - 2*e0*n0*n1*o2*t2*z1 - 2*f1*j2*n0*r0*t2*z1 - 
   2*f0*j2*n1*r0*t2*z1 - g1*i2*pow(r0,2)*t2*z1 + 
   2*n1*n2*pow(r0,2)*t2*z1 - 2*f0*j2*n0*r1*t2*z1 - 2*g0*i2*r0*r1*t2*z1 + 
   4*n0*n2*r0*r1*t2*z1 + 4*n0*n1*r0*r2*t2*z1 + 2*pow(n0,2)*r1*r2*t2*z1 + 
   2*e1*j2*n0*s0*t2*z1 + 2*e0*j2*n1*s0*t2*z1 + 2*f1*i2*r0*s0*t2*z1 - 
   2*m2*n1*r0*s0*t2*z1 + 2*f0*i2*r1*s0*t2*z1 - 2*m2*n0*r1*s0*t2*z1 - 
   e1*i2*pow(s0,2)*t2*z1 + 2*e0*j2*n0*s1*t2*z1 + 2*f0*i2*r0*s1*t2*z1 - 
   2*m2*n0*r0*s1*t2*z1 - 2*e0*i2*s0*s1*t2*z1 - 2*f1*k2*n0*o2*v0*z1 - 
   2*f0*k2*n1*o2*v0*z1 + 2*f1*j2*n0*p2*v0*z1 + 2*f0*j2*n1*p2*v0*z1 - 
   2*g1*j2*k2*r0*v0*z1 + 2*g1*i2*p2*r0*v0*z1 - 4*n1*n2*p2*r0*v0*z1 - 
   2*g0*j2*k2*r1*v0*z1 + 2*g0*i2*p2*r1*v0*z1 - 4*n0*n2*p2*r1*v0*z1 - 
   4*n0*n1*p2*r2*v0*z1 + 2*f1*j2*k2*s0*v0*z1 - 2*f1*i2*p2*s0*v0*z1 + 
   2*m2*n1*p2*s0*v0*z1 + 2*k2*n2*r1*s0*v0*z1 + 2*k2*n1*r2*s0*v0*z1 + 
   2*f0*j2*k2*s1*v0*z1 - 2*f0*i2*p2*s1*v0*z1 + 2*m2*n0*p2*s1*v0*z1 + 
   2*k2*n2*r0*s1*v0*z1 + 2*k2*n0*r2*s1*v0*z1 - 4*k2*m2*s0*s1*v0*z1 + 
   2*k2*n1*r0*s2*v0*z1 + 2*k2*n0*r1*s2*v0*z1 + 
   g1*pow(j2,2)*pow(v0,2)*z1 - g1*i2*o2*pow(v0,2)*z1 + 
   2*n1*n2*o2*pow(v0,2)*z1 - 2*j2*n2*s1*pow(v0,2)*z1 - 
   2*j2*n1*s2*pow(v0,2)*z1 + 2*i2*s1*s2*pow(v0,2)*z1 - 
   2*f0*k2*n0*o2*v1*z1 + 2*f0*j2*n0*p2*v1*z1 - 2*g0*j2*k2*r0*v1*z1 + 
   2*g0*i2*p2*r0*v1*z1 - 4*n0*n2*p2*r0*v1*z1 - 2*pow(n0,2)*p2*r2*v1*z1 + 
   2*f0*j2*k2*s0*v1*z1 - 2*f0*i2*p2*s0*v1*z1 + 2*m2*n0*p2*s0*v1*z1 + 
   2*k2*n2*r0*s0*v1*z1 + 2*k2*n0*r2*s0*v1*z1 - 2*k2*m2*pow(s0,2)*v1*z1 + 
   2*k2*n0*r0*s2*v1*z1 + 2*g0*pow(j2,2)*v0*v1*z1 - 2*g0*i2*o2*v0*v1*z1 + 
   4*n0*n2*o2*v0*v1*z1 - 4*j2*n2*s0*v0*v1*z1 - 4*j2*n0*s2*v0*v1*z1 + 
   4*i2*s0*s2*v0*v1*z1 - 4*n0*n1*p2*r0*v2*z1 - 2*pow(n0,2)*p2*r1*v2*z1 + 
   2*k2*n1*r0*s0*v2*z1 + 2*k2*n0*r1*s0*v2*z1 + 2*k2*n0*r0*s1*v2*z1 + 
   4*n0*n1*o2*v0*v2*z1 - 4*j2*n1*s0*v0*v2*z1 - 4*j2*n0*s1*v0*v2*z1 + 
   4*i2*s0*s1*v0*v2*z1 + 2*pow(n0,2)*o2*v1*v2*z1 - 4*j2*n0*s0*v1*v2*z1 + 
   2*i2*pow(s0,2)*v1*v2*z1 + 2*e1*k2*n0*o2*w0*z1 + 2*e0*k2*n1*o2*w0*z1 - 
   2*e1*j2*n0*p2*w0*z1 - 2*e0*j2*n1*p2*w0*z1 + 2*f1*j2*k2*r0*w0*z1 - 
   2*f1*i2*p2*r0*w0*z1 + 2*m2*n1*p2*r0*w0*z1 + 2*f0*j2*k2*r1*w0*z1 - 
   2*f0*i2*p2*r1*w0*z1 + 2*m2*n0*p2*r1*w0*z1 - 4*k2*n2*r0*r1*w0*z1 - 
   4*k2*n1*r0*r2*w0*z1 - 4*k2*n0*r1*r2*w0*z1 - 2*e1*j2*k2*s0*w0*z1 + 
   2*e1*i2*p2*s0*w0*z1 + 2*k2*m2*r1*s0*w0*z1 - 2*e0*j2*k2*s1*w0*z1 + 
   2*e0*i2*p2*s1*w0*z1 + 2*k2*m2*r0*s1*w0*z1 - 2*f1*pow(j2,2)*v0*w0*z1 + 
   2*f1*i2*o2*v0*w0*z1 - 2*m2*n1*o2*v0*w0*z1 + 2*j2*n2*r1*v0*w0*z1 + 
   2*j2*n1*r2*v0*w0*z1 + 2*j2*m2*s1*v0*w0*z1 - 2*i2*r2*s1*v0*w0*z1 - 
   2*i2*r1*s2*v0*w0*z1 - 2*f0*pow(j2,2)*v1*w0*z1 + 2*f0*i2*o2*v1*w0*z1 - 
   2*m2*n0*o2*v1*w0*z1 + 2*j2*n2*r0*v1*w0*z1 + 2*j2*n0*r2*v1*w0*z1 + 
   2*j2*m2*s0*v1*w0*z1 - 2*i2*r2*s0*v1*w0*z1 - 2*i2*r0*s2*v1*w0*z1 + 
   2*j2*n1*r0*v2*w0*z1 + 2*j2*n0*r1*v2*w0*z1 - 2*i2*r1*s0*v2*w0*z1 - 
   2*i2*r0*s1*v2*w0*z1 + e1*pow(j2,2)*pow(w0,2)*z1 - 
   e1*i2*o2*pow(w0,2)*z1 - 2*j2*m2*r1*pow(w0,2)*z1 + 
   2*i2*r1*r2*pow(w0,2)*z1 + 2*e0*k2*n0*o2*w1*z1 - 2*e0*j2*n0*p2*w1*z1 + 
   2*f0*j2*k2*r0*w1*z1 - 2*f0*i2*p2*r0*w1*z1 + 2*m2*n0*p2*r0*w1*z1 - 
   2*k2*n2*pow(r0,2)*w1*z1 - 4*k2*n0*r0*r2*w1*z1 - 2*e0*j2*k2*s0*w1*z1 + 
   2*e0*i2*p2*s0*w1*z1 + 2*k2*m2*r0*s0*w1*z1 - 2*f0*pow(j2,2)*v0*w1*z1 + 
   2*f0*i2*o2*v0*w1*z1 - 2*m2*n0*o2*v0*w1*z1 + 2*j2*n2*r0*v0*w1*z1 + 
   2*j2*n0*r2*v0*w1*z1 + 2*j2*m2*s0*v0*w1*z1 - 2*i2*r2*s0*v0*w1*z1 - 
   2*i2*r0*s2*v0*w1*z1 + 2*j2*n0*r0*v2*w1*z1 - 2*i2*r0*s0*v2*w1*z1 + 
   2*e0*pow(j2,2)*w0*w1*z1 - 2*e0*i2*o2*w0*w1*z1 - 4*j2*m2*r0*w0*w1*z1 + 
   4*i2*r0*r2*w0*w1*z1 - 2*k2*n1*pow(r0,2)*w2*z1 - 4*k2*n0*r0*r1*w2*z1 + 
   2*j2*n1*r0*v0*w2*z1 + 2*j2*n0*r1*v0*w2*z1 - 2*i2*r1*s0*v0*w2*z1 - 
   2*i2*r0*s1*v0*w2*z1 + 2*j2*n0*r0*v1*w2*z1 - 2*i2*r0*s0*v1*w2*z1 + 
   4*i2*r0*r1*w0*w2*z1 + 2*i2*pow(r0,2)*w1*w2*z1 + 
   pow(n1,2)*pow(r0,2)*t2*z2 + 4*n0*n1*r0*r1*t2*z2 + 
   pow(n0,2)*pow(r1,2)*t2*z2 - 2*pow(n1,2)*p2*r0*v0*z2 - 
   4*n0*n1*p2*r1*v0*z2 + 2*k2*n1*r1*s0*v0*z2 + 2*k2*n1*r0*s1*v0*z2 + 
   2*k2*n0*r1*s1*v0*z2 + pow(n1,2)*o2*pow(v0,2)*z2 - 
   2*j2*n1*s1*pow(v0,2)*z2 + i2*pow(s1,2)*pow(v0,2)*z2 - 
   4*n0*n1*p2*r0*v1*z2 - 2*pow(n0,2)*p2*r1*v1*z2 + 2*k2*n1*r0*s0*v1*z2 + 
   2*k2*n0*r1*s0*v1*z2 + 2*k2*n0*r0*s1*v1*z2 + 4*n0*n1*o2*v0*v1*z2 - 
   4*j2*n1*s0*v0*v1*z2 - 4*j2*n0*s1*v0*v1*z2 + 4*i2*s0*s1*v0*v1*z2 + 
   pow(n0,2)*o2*pow(v1,2)*z2 - 2*j2*n0*s0*pow(v1,2)*z2 + 
   i2*pow(s0,2)*pow(v1,2)*z2 - 4*k2*n1*r0*r1*w0*z2 - 
   2*k2*n0*pow(r1,2)*w0*z2 + 2*j2*n1*r1*v0*w0*z2 - 2*i2*r1*s1*v0*w0*z2 + 
   2*j2*n1*r0*v1*w0*z2 + 2*j2*n0*r1*v1*w0*z2 - 2*i2*r1*s0*v1*w0*z2 - 
   2*i2*r0*s1*v1*w0*z2 + i2*pow(r1,2)*pow(w0,2)*z2 - 
   2*k2*n1*pow(r0,2)*w1*z2 - 4*k2*n0*r0*r1*w1*z2 + 2*j2*n1*r0*v0*w1*z2 + 
   2*j2*n0*r1*v0*w1*z2 - 2*i2*r1*s0*v0*w1*z2 - 2*i2*r0*s1*v0*w1*z2 + 
   2*j2*n0*r0*v1*w1*z2 - 2*i2*r0*s0*v1*w1*z2 + 4*i2*r0*r1*w0*w1*z2 + 
   i2*pow(r0,2)*pow(w1,2)*z2
;

	k[13]  = -2*pow(c1,2)*n0*n1*pow(p2,2) - 2*c0*c1*pow(n1,2)*pow(p2,2) - 
   2*c1*d1*k2*n1*p2*r0 - 2*c1*d1*k2*n0*p2*r1 - 2*c1*d0*k2*n1*p2*r1 - 
   2*c0*d1*k2*n1*p2*r1 - 2*pow(d1,2)*pow(k2,2)*r0*r1 - 
   2*d0*d1*pow(k2,2)*pow(r1,2) + 2*pow(c1,2)*k2*n1*p2*s0 + 
   2*c1*d1*pow(k2,2)*r1*s0 + 2*pow(c1,2)*k2*n0*p2*s1 + 
   4*c0*c1*k2*n1*p2*s1 + 2*c1*d1*pow(k2,2)*r0*s1 + 
   2*c1*d0*pow(k2,2)*r1*s1 + 2*c0*d1*pow(k2,2)*r1*s1 - 
   2*pow(c1,2)*pow(k2,2)*s0*s1 - 2*c0*c1*pow(k2,2)*pow(s1,2) + 
   2*pow(c1,2)*n0*n1*o2*t2 + 2*c0*c1*pow(n1,2)*o2*t2 + 
   2*c1*d1*j2*n1*r0*t2 - 2*c1*pow(n1,2)*q2*r0*t2 + 2*c1*d1*j2*n0*r1*t2 + 
   2*c1*d0*j2*n1*r1*t2 + 2*c0*d1*j2*n1*r1*t2 - 4*c1*n0*n1*q2*r1*t2 - 
   2*c0*pow(n1,2)*q2*r1*t2 + 2*pow(d1,2)*i2*r0*r1*t2 - 
   4*d1*l2*n1*r0*r1*t2 + 2*d0*d1*i2*pow(r1,2)*t2 - 
   2*d1*l2*n0*pow(r1,2)*t2 - 2*d0*l2*n1*pow(r1,2)*t2 - 
   2*pow(c1,2)*j2*n1*s0*t2 - 2*c1*d1*i2*r1*s0*t2 + 2*c1*l2*n1*r1*s0*t2 - 
   2*pow(c1,2)*j2*n0*s1*t2 - 4*c0*c1*j2*n1*s1*t2 - 2*c1*d1*i2*r0*s1*t2 + 
   2*c1*l2*n1*r0*s1*t2 - 2*c1*d0*i2*r1*s1*t2 - 2*c0*d1*i2*r1*s1*t2 + 
   2*c1*l2*n0*r1*s1*t2 + 2*c0*l2*n1*r1*s1*t2 + 2*pow(c1,2)*i2*s0*s1*t2 + 
   2*c0*c1*i2*pow(s1,2)*t2 - 2*d1*e1*k2*n1*o2*u0 + 2*c1*f1*k2*n1*o2*u0 + 
   2*d1*e1*j2*n1*p2*u0 - 2*c1*f1*j2*n1*p2*u0 - 2*e1*pow(n1,2)*p2*q2*u0 - 
   2*d1*f1*j2*k2*r1*u0 + 2*c1*g1*j2*k2*r1*u0 + 2*d1*f1*i2*p2*r1*u0 - 
   2*c1*g1*i2*p2*r1*u0 - 2*f1*l2*n1*p2*r1*u0 - 2*d1*m2*n1*p2*r1*u0 + 
   2*c2*pow(n1,2)*p2*r1*u0 + 4*c1*n1*n2*p2*r1*u0 - 2*f1*k2*n1*q2*r1*u0 - 
   2*g1*k2*l2*pow(r1,2)*u0 + 2*d2*k2*n1*pow(r1,2)*u0 + 
   2*d1*k2*n2*pow(r1,2)*u0 + 2*c1*pow(n1,2)*p2*r2*u0 + 
   4*d1*k2*n1*r1*r2*u0 + 2*d1*e1*j2*k2*s1*u0 - 2*c1*f1*j2*k2*s1*u0 - 
   2*d1*e1*i2*p2*s1*u0 + 2*c1*f1*i2*p2*s1*u0 + 2*e1*l2*n1*p2*s1*u0 - 
   2*c1*m2*n1*p2*s1*u0 + 2*e1*k2*n1*q2*s1*u0 + 4*f1*k2*l2*r1*s1*u0 - 
   2*d1*k2*m2*r1*s1*u0 - 2*c2*k2*n1*r1*s1*u0 - 2*c1*k2*n2*r1*s1*u0 - 
   2*c1*k2*n1*r2*s1*u0 - 2*e1*k2*l2*pow(s1,2)*u0 + 
   2*c1*k2*m2*pow(s1,2)*u0 - 2*c1*k2*n1*r1*s2*u0 - 2*d1*e1*k2*n0*o2*u1 + 
   2*c1*f1*k2*n0*o2*u1 - 2*d1*e0*k2*n1*o2*u1 - 2*d0*e1*k2*n1*o2*u1 + 
   2*c1*f0*k2*n1*o2*u1 + 2*c0*f1*k2*n1*o2*u1 + 2*d1*e1*j2*n0*p2*u1 - 
   2*c1*f1*j2*n0*p2*u1 + 2*d1*e0*j2*n1*p2*u1 + 2*d0*e1*j2*n1*p2*u1 - 
   2*c1*f0*j2*n1*p2*u1 - 2*c0*f1*j2*n1*p2*u1 - 4*e1*n0*n1*p2*q2*u1 - 
   2*e0*pow(n1,2)*p2*q2*u1 - 2*d1*f1*j2*k2*r0*u1 + 2*c1*g1*j2*k2*r0*u1 + 
   2*d1*f1*i2*p2*r0*u1 - 2*c1*g1*i2*p2*r0*u1 - 2*f1*l2*n1*p2*r0*u1 - 
   2*d1*m2*n1*p2*r0*u1 + 2*c2*pow(n1,2)*p2*r0*u1 + 4*c1*n1*n2*p2*r0*u1 - 
   2*f1*k2*n1*q2*r0*u1 - 2*d1*f0*j2*k2*r1*u1 - 2*d0*f1*j2*k2*r1*u1 + 
   2*c1*g0*j2*k2*r1*u1 + 2*c0*g1*j2*k2*r1*u1 + 2*d1*f0*i2*p2*r1*u1 + 
   2*d0*f1*i2*p2*r1*u1 - 2*c1*g0*i2*p2*r1*u1 - 2*c0*g1*i2*p2*r1*u1 - 
   2*f1*l2*n0*p2*r1*u1 - 2*d1*m2*n0*p2*r1*u1 - 2*f0*l2*n1*p2*r1*u1 - 
   2*d0*m2*n1*p2*r1*u1 + 4*c2*n0*n1*p2*r1*u1 + 4*c1*n0*n2*p2*r1*u1 + 
   4*c0*n1*n2*p2*r1*u1 - 2*f1*k2*n0*q2*r1*u1 - 2*f0*k2*n1*q2*r1*u1 - 
   4*g1*k2*l2*r0*r1*u1 + 4*d2*k2*n1*r0*r1*u1 + 4*d1*k2*n2*r0*r1*u1 - 
   2*g0*k2*l2*pow(r1,2)*u1 + 2*d2*k2*n0*pow(r1,2)*u1 + 
   2*d0*k2*n2*pow(r1,2)*u1 + 4*c1*n0*n1*p2*r2*u1 + 
   2*c0*pow(n1,2)*p2*r2*u1 + 4*d1*k2*n1*r0*r2*u1 + 4*d1*k2*n0*r1*r2*u1 + 
   4*d0*k2*n1*r1*r2*u1 + 2*d1*e1*j2*k2*s0*u1 - 2*c1*f1*j2*k2*s0*u1 - 
   2*d1*e1*i2*p2*s0*u1 + 2*c1*f1*i2*p2*s0*u1 + 2*e1*l2*n1*p2*s0*u1 - 
   2*c1*m2*n1*p2*s0*u1 + 2*e1*k2*n1*q2*s0*u1 + 4*f1*k2*l2*r1*s0*u1 - 
   2*d1*k2*m2*r1*s0*u1 - 2*c2*k2*n1*r1*s0*u1 - 2*c1*k2*n2*r1*s0*u1 - 
   2*c1*k2*n1*r2*s0*u1 + 2*d1*e0*j2*k2*s1*u1 + 2*d0*e1*j2*k2*s1*u1 - 
   2*c1*f0*j2*k2*s1*u1 - 2*c0*f1*j2*k2*s1*u1 - 2*d1*e0*i2*p2*s1*u1 - 
   2*d0*e1*i2*p2*s1*u1 + 2*c1*f0*i2*p2*s1*u1 + 2*c0*f1*i2*p2*s1*u1 + 
   2*e1*l2*n0*p2*s1*u1 - 2*c1*m2*n0*p2*s1*u1 + 2*e0*l2*n1*p2*s1*u1 - 
   2*c0*m2*n1*p2*s1*u1 + 2*e1*k2*n0*q2*s1*u1 + 2*e0*k2*n1*q2*s1*u1 + 
   4*f1*k2*l2*r0*s1*u1 - 2*d1*k2*m2*r0*s1*u1 - 2*c2*k2*n1*r0*s1*u1 - 
   2*c1*k2*n2*r0*s1*u1 + 4*f0*k2*l2*r1*s1*u1 - 2*d0*k2*m2*r1*s1*u1 - 
   2*c2*k2*n0*r1*s1*u1 - 2*c0*k2*n2*r1*s1*u1 - 2*c1*k2*n0*r2*s1*u1 - 
   2*c0*k2*n1*r2*s1*u1 - 4*e1*k2*l2*s0*s1*u1 + 4*c1*k2*m2*s0*s1*u1 - 
   2*e0*k2*l2*pow(s1,2)*u1 + 2*c0*k2*m2*pow(s1,2)*u1 - 
   2*c1*k2*n1*r0*s2*u1 - 2*c1*k2*n0*r1*s2*u1 - 2*c0*k2*n1*r1*s2*u1 - 
   2*pow(f1,2)*pow(j2,2)*u0*u1 + 2*e1*g1*pow(j2,2)*u0*u1 + 
   2*pow(f1,2)*i2*o2*u0*u1 - 2*e1*g1*i2*o2*u0*u1 - 4*f1*m2*n1*o2*u0*u1 + 
   2*e2*pow(n1,2)*o2*u0*u1 + 4*e1*n1*n2*o2*u0*u1 - 4*g1*j2*m2*r1*u0*u1 + 
   4*f2*j2*n1*r1*u0*u1 + 4*f1*j2*n2*r1*u0*u1 + 2*g2*i2*pow(r1,2)*u0*u1 - 
   2*pow(n2,2)*pow(r1,2)*u0*u1 + 4*f1*j2*n1*r2*u0*u1 + 
   4*g1*i2*r1*r2*u0*u1 - 8*n1*n2*r1*r2*u0*u1 - 
   2*pow(n1,2)*pow(r2,2)*u0*u1 + 4*f1*j2*m2*s1*u0*u1 - 
   4*e2*j2*n1*s1*u0*u1 - 4*e1*j2*n2*s1*u0*u1 - 4*f2*i2*r1*s1*u0*u1 + 
   4*m2*n2*r1*s1*u0*u1 - 4*f1*i2*r2*s1*u0*u1 + 4*m2*n1*r2*s1*u0*u1 + 
   2*e2*i2*pow(s1,2)*u0*u1 - 2*pow(m2,2)*pow(s1,2)*u0*u1 - 
   4*e1*j2*n1*s2*u0*u1 - 4*f1*i2*r1*s2*u0*u1 + 4*m2*n1*r1*s2*u0*u1 + 
   4*e1*i2*s1*s2*u0*u1 - 2*f0*f1*pow(j2,2)*pow(u1,2) + 
   e1*g0*pow(j2,2)*pow(u1,2) + e0*g1*pow(j2,2)*pow(u1,2) + 
   2*f0*f1*i2*o2*pow(u1,2) - e1*g0*i2*o2*pow(u1,2) - 
   e0*g1*i2*o2*pow(u1,2) - 2*f1*m2*n0*o2*pow(u1,2) - 
   2*f0*m2*n1*o2*pow(u1,2) + 2*e2*n0*n1*o2*pow(u1,2) + 
   2*e1*n0*n2*o2*pow(u1,2) + 2*e0*n1*n2*o2*pow(u1,2) - 
   2*g1*j2*m2*r0*pow(u1,2) + 2*f2*j2*n1*r0*pow(u1,2) + 
   2*f1*j2*n2*r0*pow(u1,2) - 2*g0*j2*m2*r1*pow(u1,2) + 
   2*f2*j2*n0*r1*pow(u1,2) + 2*f0*j2*n2*r1*pow(u1,2) + 
   2*g2*i2*r0*r1*pow(u1,2) - 2*pow(n2,2)*r0*r1*pow(u1,2) + 
   2*f1*j2*n0*r2*pow(u1,2) + 2*f0*j2*n1*r2*pow(u1,2) + 
   2*g1*i2*r0*r2*pow(u1,2) - 4*n1*n2*r0*r2*pow(u1,2) + 
   2*g0*i2*r1*r2*pow(u1,2) - 4*n0*n2*r1*r2*pow(u1,2) - 
   2*n0*n1*pow(r2,2)*pow(u1,2) + 2*f1*j2*m2*s0*pow(u1,2) - 
   2*e2*j2*n1*s0*pow(u1,2) - 2*e1*j2*n2*s0*pow(u1,2) - 
   2*f2*i2*r1*s0*pow(u1,2) + 2*m2*n2*r1*s0*pow(u1,2) - 
   2*f1*i2*r2*s0*pow(u1,2) + 2*m2*n1*r2*s0*pow(u1,2) + 
   2*f0*j2*m2*s1*pow(u1,2) - 2*e2*j2*n0*s1*pow(u1,2) - 
   2*e0*j2*n2*s1*pow(u1,2) - 2*f2*i2*r0*s1*pow(u1,2) + 
   2*m2*n2*r0*s1*pow(u1,2) - 2*f0*i2*r2*s1*pow(u1,2) + 
   2*m2*n0*r2*s1*pow(u1,2) + 2*e2*i2*s0*s1*pow(u1,2) - 
   2*pow(m2,2)*s0*s1*pow(u1,2) - 2*e1*j2*n0*s2*pow(u1,2) - 
   2*e0*j2*n1*s2*pow(u1,2) - 2*f1*i2*r0*s2*pow(u1,2) + 
   2*m2*n1*r0*s2*pow(u1,2) - 2*f0*i2*r1*s2*pow(u1,2) + 
   2*m2*n0*r1*s2*pow(u1,2) + 2*e1*i2*s0*s2*pow(u1,2) + 
   2*e0*i2*s1*s2*pow(u1,2) + 2*c1*pow(n1,2)*p2*r0*u2 + 
   4*c1*n0*n1*p2*r1*u2 + 2*c0*pow(n1,2)*p2*r1*u2 + 4*d1*k2*n1*r0*r1*u2 + 
   2*d1*k2*n0*pow(r1,2)*u2 + 2*d0*k2*n1*pow(r1,2)*u2 - 
   2*c1*k2*n1*r1*s0*u2 - 2*c1*k2*n1*r0*s1*u2 - 2*c1*k2*n0*r1*s1*u2 - 
   2*c0*k2*n1*r1*s1*u2 + 2*e1*pow(n1,2)*o2*u0*u2 + 4*f1*j2*n1*r1*u0*u2 + 
   2*g1*i2*pow(r1,2)*u0*u2 - 4*n1*n2*pow(r1,2)*u0*u2 - 
   4*pow(n1,2)*r1*r2*u0*u2 - 4*e1*j2*n1*s1*u0*u2 - 4*f1*i2*r1*s1*u0*u2 + 
   4*m2*n1*r1*s1*u0*u2 + 2*e1*i2*pow(s1,2)*u0*u2 + 4*e1*n0*n1*o2*u1*u2 + 
   2*e0*pow(n1,2)*o2*u1*u2 + 4*f1*j2*n1*r0*u1*u2 + 4*f1*j2*n0*r1*u1*u2 + 
   4*f0*j2*n1*r1*u1*u2 + 4*g1*i2*r0*r1*u1*u2 - 8*n1*n2*r0*r1*u1*u2 + 
   2*g0*i2*pow(r1,2)*u1*u2 - 4*n0*n2*pow(r1,2)*u1*u2 - 
   4*pow(n1,2)*r0*r2*u1*u2 - 8*n0*n1*r1*r2*u1*u2 - 4*e1*j2*n1*s0*u1*u2 - 
   4*f1*i2*r1*s0*u1*u2 + 4*m2*n1*r1*s0*u1*u2 - 4*e1*j2*n0*s1*u1*u2 - 
   4*e0*j2*n1*s1*u1*u2 - 4*f1*i2*r0*s1*u1*u2 + 4*m2*n1*r0*s1*u1*u2 - 
   4*f0*i2*r1*s1*u1*u2 + 4*m2*n0*r1*s1*u1*u2 + 4*e1*i2*s0*s1*u1*u2 + 
   2*e0*i2*pow(s1,2)*u1*u2 - 2*pow(n1,2)*r0*r1*pow(u2,2) - 
   2*n0*n1*pow(r1,2)*pow(u2,2) + 2*c1*d1*k2*n1*o2*v0 - 
   2*c1*d1*j2*n1*p2*v0 + 2*c1*pow(n1,2)*p2*q2*v0 + 
   2*pow(d1,2)*j2*k2*r1*v0 - 2*pow(d1,2)*i2*p2*r1*v0 + 
   4*d1*l2*n1*p2*r1*v0 - 2*d1*k2*n1*q2*r1*v0 - 2*c1*d1*j2*k2*s1*v0 + 
   2*c1*d1*i2*p2*s1*v0 - 2*c1*l2*n1*p2*s1*v0 - 2*c1*k2*n1*q2*s1*v0 - 
   2*d1*k2*l2*r1*s1*v0 + 2*c1*k2*l2*pow(s1,2)*v0 + 
   2*d1*f1*pow(j2,2)*u1*v0 - 2*c1*g1*pow(j2,2)*u1*v0 - 
   2*d1*f1*i2*o2*u1*v0 + 2*c1*g1*i2*o2*u1*v0 + 2*f1*l2*n1*o2*u1*v0 + 
   2*d1*m2*n1*o2*u1*v0 - 2*c2*pow(n1,2)*o2*u1*v0 - 4*c1*n1*n2*o2*u1*v0 - 
   2*f1*j2*n1*q2*u1*v0 + 2*g1*j2*l2*r1*u1*v0 - 2*d2*j2*n1*r1*u1*v0 - 
   2*d1*j2*n2*r1*u1*v0 - 2*g1*i2*q2*r1*u1*v0 + 4*n1*n2*q2*r1*u1*v0 - 
   2*d1*j2*n1*r2*u1*v0 + 2*pow(n1,2)*q2*r2*u1*v0 - 2*f1*j2*l2*s1*u1*v0 - 
   2*d1*j2*m2*s1*u1*v0 + 4*c2*j2*n1*s1*u1*v0 + 4*c1*j2*n2*s1*u1*v0 + 
   2*f1*i2*q2*s1*u1*v0 - 2*m2*n1*q2*s1*u1*v0 + 2*d2*i2*r1*s1*u1*v0 - 
   2*l2*n2*r1*s1*u1*v0 + 2*d1*i2*r2*s1*u1*v0 - 2*l2*n1*r2*s1*u1*v0 - 
   2*c2*i2*pow(s1,2)*u1*v0 + 2*l2*m2*pow(s1,2)*u1*v0 + 
   4*c1*j2*n1*s2*u1*v0 + 2*d1*i2*r1*s2*u1*v0 - 2*l2*n1*r1*s2*u1*v0 - 
   4*c1*i2*s1*s2*u1*v0 - 2*c1*pow(n1,2)*o2*u2*v0 - 2*d1*j2*n1*r1*u2*v0 + 
   2*pow(n1,2)*q2*r1*u2*v0 + 4*c1*j2*n1*s1*u2*v0 + 2*d1*i2*r1*s1*u2*v0 - 
   2*l2*n1*r1*s1*u2*v0 - 2*c1*i2*pow(s1,2)*u2*v0 + 2*c1*d1*k2*n0*o2*v1 + 
   2*c1*d0*k2*n1*o2*v1 + 2*c0*d1*k2*n1*o2*v1 - 2*c1*d1*j2*n0*p2*v1 - 
   2*c1*d0*j2*n1*p2*v1 - 2*c0*d1*j2*n1*p2*v1 + 4*c1*n0*n1*p2*q2*v1 + 
   2*c0*pow(n1,2)*p2*q2*v1 + 2*pow(d1,2)*j2*k2*r0*v1 - 
   2*pow(d1,2)*i2*p2*r0*v1 + 4*d1*l2*n1*p2*r0*v1 - 2*d1*k2*n1*q2*r0*v1 + 
   4*d0*d1*j2*k2*r1*v1 - 4*d0*d1*i2*p2*r1*v1 + 4*d1*l2*n0*p2*r1*v1 + 
   4*d0*l2*n1*p2*r1*v1 - 2*d1*k2*n0*q2*r1*v1 - 2*d0*k2*n1*q2*r1*v1 - 
   2*c1*d1*j2*k2*s0*v1 + 2*c1*d1*i2*p2*s0*v1 - 2*c1*l2*n1*p2*s0*v1 - 
   2*c1*k2*n1*q2*s0*v1 - 2*d1*k2*l2*r1*s0*v1 - 2*c1*d0*j2*k2*s1*v1 - 
   2*c0*d1*j2*k2*s1*v1 + 2*c1*d0*i2*p2*s1*v1 + 2*c0*d1*i2*p2*s1*v1 - 
   2*c1*l2*n0*p2*s1*v1 - 2*c0*l2*n1*p2*s1*v1 - 2*c1*k2*n0*q2*s1*v1 - 
   2*c0*k2*n1*q2*s1*v1 - 2*d1*k2*l2*r0*s1*v1 - 2*d0*k2*l2*r1*s1*v1 + 
   4*c1*k2*l2*s0*s1*v1 + 2*c0*k2*l2*pow(s1,2)*v1 + 
   2*d1*f1*pow(j2,2)*u0*v1 - 2*c1*g1*pow(j2,2)*u0*v1 - 
   2*d1*f1*i2*o2*u0*v1 + 2*c1*g1*i2*o2*u0*v1 + 2*f1*l2*n1*o2*u0*v1 + 
   2*d1*m2*n1*o2*u0*v1 - 2*c2*pow(n1,2)*o2*u0*v1 - 4*c1*n1*n2*o2*u0*v1 - 
   2*f1*j2*n1*q2*u0*v1 + 2*g1*j2*l2*r1*u0*v1 - 2*d2*j2*n1*r1*u0*v1 - 
   2*d1*j2*n2*r1*u0*v1 - 2*g1*i2*q2*r1*u0*v1 + 4*n1*n2*q2*r1*u0*v1 - 
   2*d1*j2*n1*r2*u0*v1 + 2*pow(n1,2)*q2*r2*u0*v1 - 2*f1*j2*l2*s1*u0*v1 - 
   2*d1*j2*m2*s1*u0*v1 + 4*c2*j2*n1*s1*u0*v1 + 4*c1*j2*n2*s1*u0*v1 + 
   2*f1*i2*q2*s1*u0*v1 - 2*m2*n1*q2*s1*u0*v1 + 2*d2*i2*r1*s1*u0*v1 - 
   2*l2*n2*r1*s1*u0*v1 + 2*d1*i2*r2*s1*u0*v1 - 2*l2*n1*r2*s1*u0*v1 - 
   2*c2*i2*pow(s1,2)*u0*v1 + 2*l2*m2*pow(s1,2)*u0*v1 + 
   4*c1*j2*n1*s2*u0*v1 + 2*d1*i2*r1*s2*u0*v1 - 2*l2*n1*r1*s2*u0*v1 - 
   4*c1*i2*s1*s2*u0*v1 + 2*d1*f0*pow(j2,2)*u1*v1 + 
   2*d0*f1*pow(j2,2)*u1*v1 - 2*c1*g0*pow(j2,2)*u1*v1 - 
   2*c0*g1*pow(j2,2)*u1*v1 - 2*d1*f0*i2*o2*u1*v1 - 2*d0*f1*i2*o2*u1*v1 + 
   2*c1*g0*i2*o2*u1*v1 + 2*c0*g1*i2*o2*u1*v1 + 2*f1*l2*n0*o2*u1*v1 + 
   2*d1*m2*n0*o2*u1*v1 + 2*f0*l2*n1*o2*u1*v1 + 2*d0*m2*n1*o2*u1*v1 - 
   4*c2*n0*n1*o2*u1*v1 - 4*c1*n0*n2*o2*u1*v1 - 4*c0*n1*n2*o2*u1*v1 - 
   2*f1*j2*n0*q2*u1*v1 - 2*f0*j2*n1*q2*u1*v1 + 2*g1*j2*l2*r0*u1*v1 - 
   2*d2*j2*n1*r0*u1*v1 - 2*d1*j2*n2*r0*u1*v1 - 2*g1*i2*q2*r0*u1*v1 + 
   4*n1*n2*q2*r0*u1*v1 + 2*g0*j2*l2*r1*u1*v1 - 2*d2*j2*n0*r1*u1*v1 - 
   2*d0*j2*n2*r1*u1*v1 - 2*g0*i2*q2*r1*u1*v1 + 4*n0*n2*q2*r1*u1*v1 - 
   2*d1*j2*n0*r2*u1*v1 - 2*d0*j2*n1*r2*u1*v1 + 4*n0*n1*q2*r2*u1*v1 - 
   2*f1*j2*l2*s0*u1*v1 - 2*d1*j2*m2*s0*u1*v1 + 4*c2*j2*n1*s0*u1*v1 + 
   4*c1*j2*n2*s0*u1*v1 + 2*f1*i2*q2*s0*u1*v1 - 2*m2*n1*q2*s0*u1*v1 + 
   2*d2*i2*r1*s0*u1*v1 - 2*l2*n2*r1*s0*u1*v1 + 2*d1*i2*r2*s0*u1*v1 - 
   2*l2*n1*r2*s0*u1*v1 - 2*f0*j2*l2*s1*u1*v1 - 2*d0*j2*m2*s1*u1*v1 + 
   4*c2*j2*n0*s1*u1*v1 + 4*c0*j2*n2*s1*u1*v1 + 2*f0*i2*q2*s1*u1*v1 - 
   2*m2*n0*q2*s1*u1*v1 + 2*d2*i2*r0*s1*u1*v1 - 2*l2*n2*r0*s1*u1*v1 + 
   2*d0*i2*r2*s1*u1*v1 - 2*l2*n0*r2*s1*u1*v1 - 4*c2*i2*s0*s1*u1*v1 + 
   4*l2*m2*s0*s1*u1*v1 + 4*c1*j2*n0*s2*u1*v1 + 4*c0*j2*n1*s2*u1*v1 + 
   2*d1*i2*r0*s2*u1*v1 - 2*l2*n1*r0*s2*u1*v1 + 2*d0*i2*r1*s2*u1*v1 - 
   2*l2*n0*r1*s2*u1*v1 - 4*c1*i2*s0*s2*u1*v1 - 4*c0*i2*s1*s2*u1*v1 - 
   4*c1*n0*n1*o2*u2*v1 - 2*c0*pow(n1,2)*o2*u2*v1 - 2*d1*j2*n1*r0*u2*v1 + 
   2*pow(n1,2)*q2*r0*u2*v1 - 2*d1*j2*n0*r1*u2*v1 - 2*d0*j2*n1*r1*u2*v1 + 
   4*n0*n1*q2*r1*u2*v1 + 4*c1*j2*n1*s0*u2*v1 + 2*d1*i2*r1*s0*u2*v1 - 
   2*l2*n1*r1*s0*u2*v1 + 4*c1*j2*n0*s1*u2*v1 + 4*c0*j2*n1*s1*u2*v1 + 
   2*d1*i2*r0*s1*u2*v1 - 2*l2*n1*r0*s1*u2*v1 + 2*d0*i2*r1*s1*u2*v1 - 
   2*l2*n0*r1*s1*u2*v1 - 4*c1*i2*s0*s1*u2*v1 - 2*c0*i2*pow(s1,2)*u2*v1 - 
   2*pow(d1,2)*pow(j2,2)*v0*v1 + 2*pow(d1,2)*i2*o2*v0*v1 - 
   4*d1*l2*n1*o2*v0*v1 + 4*d1*j2*n1*q2*v0*v1 - 
   2*pow(n1,2)*pow(q2,2)*v0*v1 + 4*d1*j2*l2*s1*v0*v1 - 
   4*d1*i2*q2*s1*v0*v1 + 4*l2*n1*q2*s1*v0*v1 - 
   2*pow(l2,2)*pow(s1,2)*v0*v1 - 2*d0*d1*pow(j2,2)*pow(v1,2) + 
   2*d0*d1*i2*o2*pow(v1,2) - 2*d1*l2*n0*o2*pow(v1,2) - 
   2*d0*l2*n1*o2*pow(v1,2) + 2*d1*j2*n0*q2*pow(v1,2) + 
   2*d0*j2*n1*q2*pow(v1,2) - 2*n0*n1*pow(q2,2)*pow(v1,2) + 
   2*d1*j2*l2*s0*pow(v1,2) - 2*d1*i2*q2*s0*pow(v1,2) + 
   2*l2*n1*q2*s0*pow(v1,2) + 2*d0*j2*l2*s1*pow(v1,2) - 
   2*d0*i2*q2*s1*pow(v1,2) + 2*l2*n0*q2*s1*pow(v1,2) - 
   2*pow(l2,2)*s0*s1*pow(v1,2) - 2*c1*pow(n1,2)*o2*u0*v2 - 
   2*d1*j2*n1*r1*u0*v2 + 2*pow(n1,2)*q2*r1*u0*v2 + 4*c1*j2*n1*s1*u0*v2 + 
   2*d1*i2*r1*s1*u0*v2 - 2*l2*n1*r1*s1*u0*v2 - 2*c1*i2*pow(s1,2)*u0*v2 - 
   4*c1*n0*n1*o2*u1*v2 - 2*c0*pow(n1,2)*o2*u1*v2 - 2*d1*j2*n1*r0*u1*v2 + 
   2*pow(n1,2)*q2*r0*u1*v2 - 2*d1*j2*n0*r1*u1*v2 - 2*d0*j2*n1*r1*u1*v2 + 
   4*n0*n1*q2*r1*u1*v2 + 4*c1*j2*n1*s0*u1*v2 + 2*d1*i2*r1*s0*u1*v2 - 
   2*l2*n1*r1*s0*u1*v2 + 4*c1*j2*n0*s1*u1*v2 + 4*c0*j2*n1*s1*u1*v2 + 
   2*d1*i2*r0*s1*u1*v2 - 2*l2*n1*r0*s1*u1*v2 + 2*d0*i2*r1*s1*u1*v2 - 
   2*l2*n0*r1*s1*u1*v2 - 4*c1*i2*s0*s1*u1*v2 - 2*c0*i2*pow(s1,2)*u1*v2 - 
   2*pow(c1,2)*k2*n1*o2*w0 + 2*pow(c1,2)*j2*n1*p2*w0 - 
   2*c1*d1*j2*k2*r1*w0 + 2*c1*d1*i2*p2*r1*w0 - 2*c1*l2*n1*p2*r1*w0 + 
   4*c1*k2*n1*q2*r1*w0 + 2*d1*k2*l2*pow(r1,2)*w0 + 
   2*pow(c1,2)*j2*k2*s1*w0 - 2*pow(c1,2)*i2*p2*s1*w0 - 
   2*c1*k2*l2*r1*s1*w0 - 2*d1*e1*pow(j2,2)*u1*w0 + 
   2*c1*f1*pow(j2,2)*u1*w0 + 2*d1*e1*i2*o2*u1*w0 - 2*c1*f1*i2*o2*u1*w0 - 
   2*e1*l2*n1*o2*u1*w0 + 2*c1*m2*n1*o2*u1*w0 + 2*e1*j2*n1*q2*u1*w0 - 
   2*f1*j2*l2*r1*u1*w0 + 4*d1*j2*m2*r1*u1*w0 - 2*c2*j2*n1*r1*u1*w0 - 
   2*c1*j2*n2*r1*u1*w0 + 2*f1*i2*q2*r1*u1*w0 - 2*m2*n1*q2*r1*u1*w0 - 
   2*d2*i2*pow(r1,2)*u1*w0 + 2*l2*n2*pow(r1,2)*u1*w0 - 
   2*c1*j2*n1*r2*u1*w0 - 4*d1*i2*r1*r2*u1*w0 + 4*l2*n1*r1*r2*u1*w0 + 
   2*e1*j2*l2*s1*u1*w0 - 2*c1*j2*m2*s1*u1*w0 - 2*e1*i2*q2*s1*u1*w0 + 
   2*c2*i2*r1*s1*u1*w0 - 2*l2*m2*r1*s1*u1*w0 + 2*c1*i2*r2*s1*u1*w0 + 
   2*c1*i2*r1*s2*u1*w0 - 2*c1*j2*n1*r1*u2*w0 - 2*d1*i2*pow(r1,2)*u2*w0 + 
   2*l2*n1*pow(r1,2)*u2*w0 + 2*c1*i2*r1*s1*u2*w0 + 
   2*c1*d1*pow(j2,2)*v1*w0 - 2*c1*d1*i2*o2*v1*w0 + 2*c1*l2*n1*o2*v1*w0 - 
   2*c1*j2*n1*q2*v1*w0 - 2*d1*j2*l2*r1*v1*w0 + 2*d1*i2*q2*r1*v1*w0 - 
   2*l2*n1*q2*r1*v1*w0 - 2*c1*j2*l2*s1*v1*w0 + 2*c1*i2*q2*s1*v1*w0 + 
   2*pow(l2,2)*r1*s1*v1*w0 - 2*pow(c1,2)*k2*n0*o2*w1 - 
   4*c0*c1*k2*n1*o2*w1 + 2*pow(c1,2)*j2*n0*p2*w1 + 4*c0*c1*j2*n1*p2*w1 - 
   2*c1*d1*j2*k2*r0*w1 + 2*c1*d1*i2*p2*r0*w1 - 2*c1*l2*n1*p2*r0*w1 + 
   4*c1*k2*n1*q2*r0*w1 - 2*c1*d0*j2*k2*r1*w1 - 2*c0*d1*j2*k2*r1*w1 + 
   2*c1*d0*i2*p2*r1*w1 + 2*c0*d1*i2*p2*r1*w1 - 2*c1*l2*n0*p2*r1*w1 - 
   2*c0*l2*n1*p2*r1*w1 + 4*c1*k2*n0*q2*r1*w1 + 4*c0*k2*n1*q2*r1*w1 + 
   4*d1*k2*l2*r0*r1*w1 + 2*d0*k2*l2*pow(r1,2)*w1 + 
   2*pow(c1,2)*j2*k2*s0*w1 - 2*pow(c1,2)*i2*p2*s0*w1 - 
   2*c1*k2*l2*r1*s0*w1 + 4*c0*c1*j2*k2*s1*w1 - 4*c0*c1*i2*p2*s1*w1 - 
   2*c1*k2*l2*r0*s1*w1 - 2*c0*k2*l2*r1*s1*w1 - 2*d1*e1*pow(j2,2)*u0*w1 + 
   2*c1*f1*pow(j2,2)*u0*w1 + 2*d1*e1*i2*o2*u0*w1 - 2*c1*f1*i2*o2*u0*w1 - 
   2*e1*l2*n1*o2*u0*w1 + 2*c1*m2*n1*o2*u0*w1 + 2*e1*j2*n1*q2*u0*w1 - 
   2*f1*j2*l2*r1*u0*w1 + 4*d1*j2*m2*r1*u0*w1 - 2*c2*j2*n1*r1*u0*w1 - 
   2*c1*j2*n2*r1*u0*w1 + 2*f1*i2*q2*r1*u0*w1 - 2*m2*n1*q2*r1*u0*w1 - 
   2*d2*i2*pow(r1,2)*u0*w1 + 2*l2*n2*pow(r1,2)*u0*w1 - 
   2*c1*j2*n1*r2*u0*w1 - 4*d1*i2*r1*r2*u0*w1 + 4*l2*n1*r1*r2*u0*w1 + 
   2*e1*j2*l2*s1*u0*w1 - 2*c1*j2*m2*s1*u0*w1 - 2*e1*i2*q2*s1*u0*w1 + 
   2*c2*i2*r1*s1*u0*w1 - 2*l2*m2*r1*s1*u0*w1 + 2*c1*i2*r2*s1*u0*w1 + 
   2*c1*i2*r1*s2*u0*w1 - 2*d1*e0*pow(j2,2)*u1*w1 - 
   2*d0*e1*pow(j2,2)*u1*w1 + 2*c1*f0*pow(j2,2)*u1*w1 + 
   2*c0*f1*pow(j2,2)*u1*w1 + 2*d1*e0*i2*o2*u1*w1 + 2*d0*e1*i2*o2*u1*w1 - 
   2*c1*f0*i2*o2*u1*w1 - 2*c0*f1*i2*o2*u1*w1 - 2*e1*l2*n0*o2*u1*w1 + 
   2*c1*m2*n0*o2*u1*w1 - 2*e0*l2*n1*o2*u1*w1 + 2*c0*m2*n1*o2*u1*w1 + 
   2*e1*j2*n0*q2*u1*w1 + 2*e0*j2*n1*q2*u1*w1 - 2*f1*j2*l2*r0*u1*w1 + 
   4*d1*j2*m2*r0*u1*w1 - 2*c2*j2*n1*r0*u1*w1 - 2*c1*j2*n2*r0*u1*w1 + 
   2*f1*i2*q2*r0*u1*w1 - 2*m2*n1*q2*r0*u1*w1 - 2*f0*j2*l2*r1*u1*w1 + 
   4*d0*j2*m2*r1*u1*w1 - 2*c2*j2*n0*r1*u1*w1 - 2*c0*j2*n2*r1*u1*w1 + 
   2*f0*i2*q2*r1*u1*w1 - 2*m2*n0*q2*r1*u1*w1 - 4*d2*i2*r0*r1*u1*w1 + 
   4*l2*n2*r0*r1*u1*w1 - 2*c1*j2*n0*r2*u1*w1 - 2*c0*j2*n1*r2*u1*w1 - 
   4*d1*i2*r0*r2*u1*w1 + 4*l2*n1*r0*r2*u1*w1 - 4*d0*i2*r1*r2*u1*w1 + 
   4*l2*n0*r1*r2*u1*w1 + 2*e1*j2*l2*s0*u1*w1 - 2*c1*j2*m2*s0*u1*w1 - 
   2*e1*i2*q2*s0*u1*w1 + 2*c2*i2*r1*s0*u1*w1 - 2*l2*m2*r1*s0*u1*w1 + 
   2*c1*i2*r2*s0*u1*w1 + 2*e0*j2*l2*s1*u1*w1 - 2*c0*j2*m2*s1*u1*w1 - 
   2*e0*i2*q2*s1*u1*w1 + 2*c2*i2*r0*s1*u1*w1 - 2*l2*m2*r0*s1*u1*w1 + 
   2*c0*i2*r2*s1*u1*w1 + 2*c1*i2*r0*s2*u1*w1 + 2*c0*i2*r1*s2*u1*w1 - 
   2*c1*j2*n1*r0*u2*w1 - 2*c1*j2*n0*r1*u2*w1 - 2*c0*j2*n1*r1*u2*w1 - 
   4*d1*i2*r0*r1*u2*w1 + 4*l2*n1*r0*r1*u2*w1 - 2*d0*i2*pow(r1,2)*u2*w1 + 
   2*l2*n0*pow(r1,2)*u2*w1 + 2*c1*i2*r1*s0*u2*w1 + 2*c1*i2*r0*s1*u2*w1 + 
   2*c0*i2*r1*s1*u2*w1 + 2*c1*d1*pow(j2,2)*v0*w1 - 2*c1*d1*i2*o2*v0*w1 + 
   2*c1*l2*n1*o2*v0*w1 - 2*c1*j2*n1*q2*v0*w1 - 2*d1*j2*l2*r1*v0*w1 + 
   2*d1*i2*q2*r1*v0*w1 - 2*l2*n1*q2*r1*v0*w1 - 2*c1*j2*l2*s1*v0*w1 + 
   2*c1*i2*q2*s1*v0*w1 + 2*pow(l2,2)*r1*s1*v0*w1 + 
   2*c1*d0*pow(j2,2)*v1*w1 + 2*c0*d1*pow(j2,2)*v1*w1 - 
   2*c1*d0*i2*o2*v1*w1 - 2*c0*d1*i2*o2*v1*w1 + 2*c1*l2*n0*o2*v1*w1 + 
   2*c0*l2*n1*o2*v1*w1 - 2*c1*j2*n0*q2*v1*w1 - 2*c0*j2*n1*q2*v1*w1 - 
   2*d1*j2*l2*r0*v1*w1 + 2*d1*i2*q2*r0*v1*w1 - 2*l2*n1*q2*r0*v1*w1 - 
   2*d0*j2*l2*r1*v1*w1 + 2*d0*i2*q2*r1*v1*w1 - 2*l2*n0*q2*r1*v1*w1 - 
   2*c1*j2*l2*s0*v1*w1 + 2*c1*i2*q2*s0*v1*w1 + 2*pow(l2,2)*r1*s0*v1*w1 - 
   2*c0*j2*l2*s1*v1*w1 + 2*c0*i2*q2*s1*v1*w1 + 2*pow(l2,2)*r0*s1*v1*w1 - 
   2*pow(c1,2)*pow(j2,2)*w0*w1 + 2*pow(c1,2)*i2*o2*w0*w1 + 
   4*c1*j2*l2*r1*w0*w1 - 4*c1*i2*q2*r1*w0*w1 - 
   2*pow(l2,2)*pow(r1,2)*w0*w1 - 2*c0*c1*pow(j2,2)*pow(w1,2) + 
   2*c0*c1*i2*o2*pow(w1,2) + 2*c1*j2*l2*r0*pow(w1,2) - 
   2*c1*i2*q2*r0*pow(w1,2) + 2*c0*j2*l2*r1*pow(w1,2) - 
   2*c0*i2*q2*r1*pow(w1,2) - 2*pow(l2,2)*r0*r1*pow(w1,2) - 
   2*c1*j2*n1*r1*u0*w2 - 2*d1*i2*pow(r1,2)*u0*w2 + 
   2*l2*n1*pow(r1,2)*u0*w2 + 2*c1*i2*r1*s1*u0*w2 - 2*c1*j2*n1*r0*u1*w2 - 
   2*c1*j2*n0*r1*u1*w2 - 2*c0*j2*n1*r1*u1*w2 - 4*d1*i2*r0*r1*u1*w2 + 
   4*l2*n1*r0*r1*u1*w2 - 2*d0*i2*pow(r1,2)*u1*w2 + 
   2*l2*n0*pow(r1,2)*u1*w2 + 2*c1*i2*r1*s0*u1*w2 + 2*c1*i2*r0*s1*u1*w2 + 
   2*c0*i2*r1*s1*u1*w2 + e1*pow(n1,2)*pow(p2,2)*z0 + 
   2*f1*k2*n1*p2*r1*z0 + g1*pow(k2,2)*pow(r1,2)*z0 - 
   2*e1*k2*n1*p2*s1*z0 - 2*f1*pow(k2,2)*r1*s1*z0 + 
   e1*pow(k2,2)*pow(s1,2)*z0 - e1*pow(n1,2)*o2*t2*z0 - 
   2*f1*j2*n1*r1*t2*z0 - g1*i2*pow(r1,2)*t2*z0 + 
   2*n1*n2*pow(r1,2)*t2*z0 + 2*pow(n1,2)*r1*r2*t2*z0 + 
   2*e1*j2*n1*s1*t2*z0 + 2*f1*i2*r1*s1*t2*z0 - 2*m2*n1*r1*s1*t2*z0 - 
   e1*i2*pow(s1,2)*t2*z0 - 2*f1*k2*n1*o2*v1*z0 + 2*f1*j2*n1*p2*v1*z0 - 
   2*g1*j2*k2*r1*v1*z0 + 2*g1*i2*p2*r1*v1*z0 - 4*n1*n2*p2*r1*v1*z0 - 
   2*pow(n1,2)*p2*r2*v1*z0 + 2*f1*j2*k2*s1*v1*z0 - 2*f1*i2*p2*s1*v1*z0 + 
   2*m2*n1*p2*s1*v1*z0 + 2*k2*n2*r1*s1*v1*z0 + 2*k2*n1*r2*s1*v1*z0 - 
   2*k2*m2*pow(s1,2)*v1*z0 + 2*k2*n1*r1*s2*v1*z0 + 
   g1*pow(j2,2)*pow(v1,2)*z0 - g1*i2*o2*pow(v1,2)*z0 + 
   2*n1*n2*o2*pow(v1,2)*z0 - 2*j2*n2*s1*pow(v1,2)*z0 - 
   2*j2*n1*s2*pow(v1,2)*z0 + 2*i2*s1*s2*pow(v1,2)*z0 - 
   2*pow(n1,2)*p2*r1*v2*z0 + 2*k2*n1*r1*s1*v2*z0 + 
   2*pow(n1,2)*o2*v1*v2*z0 - 4*j2*n1*s1*v1*v2*z0 + 
   2*i2*pow(s1,2)*v1*v2*z0 + 2*e1*k2*n1*o2*w1*z0 - 2*e1*j2*n1*p2*w1*z0 + 
   2*f1*j2*k2*r1*w1*z0 - 2*f1*i2*p2*r1*w1*z0 + 2*m2*n1*p2*r1*w1*z0 - 
   2*k2*n2*pow(r1,2)*w1*z0 - 4*k2*n1*r1*r2*w1*z0 - 2*e1*j2*k2*s1*w1*z0 + 
   2*e1*i2*p2*s1*w1*z0 + 2*k2*m2*r1*s1*w1*z0 - 2*f1*pow(j2,2)*v1*w1*z0 + 
   2*f1*i2*o2*v1*w1*z0 - 2*m2*n1*o2*v1*w1*z0 + 2*j2*n2*r1*v1*w1*z0 + 
   2*j2*n1*r2*v1*w1*z0 + 2*j2*m2*s1*v1*w1*z0 - 2*i2*r2*s1*v1*w1*z0 - 
   2*i2*r1*s2*v1*w1*z0 + 2*j2*n1*r1*v2*w1*z0 - 2*i2*r1*s1*v2*w1*z0 + 
   e1*pow(j2,2)*pow(w1,2)*z0 - e1*i2*o2*pow(w1,2)*z0 - 
   2*j2*m2*r1*pow(w1,2)*z0 + 2*i2*r1*r2*pow(w1,2)*z0 - 
   2*k2*n1*pow(r1,2)*w2*z0 + 2*j2*n1*r1*v1*w2*z0 - 2*i2*r1*s1*v1*w2*z0 + 
   2*i2*pow(r1,2)*w1*w2*z0 + 2*e1*n0*n1*pow(p2,2)*z1 + 
   e0*pow(n1,2)*pow(p2,2)*z1 + 2*f1*k2*n1*p2*r0*z1 + 
   2*f1*k2*n0*p2*r1*z1 + 2*f0*k2*n1*p2*r1*z1 + 2*g1*pow(k2,2)*r0*r1*z1 + 
   g0*pow(k2,2)*pow(r1,2)*z1 - 2*e1*k2*n1*p2*s0*z1 - 
   2*f1*pow(k2,2)*r1*s0*z1 - 2*e1*k2*n0*p2*s1*z1 - 2*e0*k2*n1*p2*s1*z1 - 
   2*f1*pow(k2,2)*r0*s1*z1 - 2*f0*pow(k2,2)*r1*s1*z1 + 
   2*e1*pow(k2,2)*s0*s1*z1 + e0*pow(k2,2)*pow(s1,2)*z1 - 
   2*e1*n0*n1*o2*t2*z1 - e0*pow(n1,2)*o2*t2*z1 - 2*f1*j2*n1*r0*t2*z1 - 
   2*f1*j2*n0*r1*t2*z1 - 2*f0*j2*n1*r1*t2*z1 - 2*g1*i2*r0*r1*t2*z1 + 
   4*n1*n2*r0*r1*t2*z1 - g0*i2*pow(r1,2)*t2*z1 + 
   2*n0*n2*pow(r1,2)*t2*z1 + 2*pow(n1,2)*r0*r2*t2*z1 + 
   4*n0*n1*r1*r2*t2*z1 + 2*e1*j2*n1*s0*t2*z1 + 2*f1*i2*r1*s0*t2*z1 - 
   2*m2*n1*r1*s0*t2*z1 + 2*e1*j2*n0*s1*t2*z1 + 2*e0*j2*n1*s1*t2*z1 + 
   2*f1*i2*r0*s1*t2*z1 - 2*m2*n1*r0*s1*t2*z1 + 2*f0*i2*r1*s1*t2*z1 - 
   2*m2*n0*r1*s1*t2*z1 - 2*e1*i2*s0*s1*t2*z1 - e0*i2*pow(s1,2)*t2*z1 - 
   2*f1*k2*n1*o2*v0*z1 + 2*f1*j2*n1*p2*v0*z1 - 2*g1*j2*k2*r1*v0*z1 + 
   2*g1*i2*p2*r1*v0*z1 - 4*n1*n2*p2*r1*v0*z1 - 2*pow(n1,2)*p2*r2*v0*z1 + 
   2*f1*j2*k2*s1*v0*z1 - 2*f1*i2*p2*s1*v0*z1 + 2*m2*n1*p2*s1*v0*z1 + 
   2*k2*n2*r1*s1*v0*z1 + 2*k2*n1*r2*s1*v0*z1 - 2*k2*m2*pow(s1,2)*v0*z1 + 
   2*k2*n1*r1*s2*v0*z1 - 2*f1*k2*n0*o2*v1*z1 - 2*f0*k2*n1*o2*v1*z1 + 
   2*f1*j2*n0*p2*v1*z1 + 2*f0*j2*n1*p2*v1*z1 - 2*g1*j2*k2*r0*v1*z1 + 
   2*g1*i2*p2*r0*v1*z1 - 4*n1*n2*p2*r0*v1*z1 - 2*g0*j2*k2*r1*v1*z1 + 
   2*g0*i2*p2*r1*v1*z1 - 4*n0*n2*p2*r1*v1*z1 - 4*n0*n1*p2*r2*v1*z1 + 
   2*f1*j2*k2*s0*v1*z1 - 2*f1*i2*p2*s0*v1*z1 + 2*m2*n1*p2*s0*v1*z1 + 
   2*k2*n2*r1*s0*v1*z1 + 2*k2*n1*r2*s0*v1*z1 + 2*f0*j2*k2*s1*v1*z1 - 
   2*f0*i2*p2*s1*v1*z1 + 2*m2*n0*p2*s1*v1*z1 + 2*k2*n2*r0*s1*v1*z1 + 
   2*k2*n0*r2*s1*v1*z1 - 4*k2*m2*s0*s1*v1*z1 + 2*k2*n1*r0*s2*v1*z1 + 
   2*k2*n0*r1*s2*v1*z1 + 2*g1*pow(j2,2)*v0*v1*z1 - 2*g1*i2*o2*v0*v1*z1 + 
   4*n1*n2*o2*v0*v1*z1 - 4*j2*n2*s1*v0*v1*z1 - 4*j2*n1*s2*v0*v1*z1 + 
   4*i2*s1*s2*v0*v1*z1 + g0*pow(j2,2)*pow(v1,2)*z1 - 
   g0*i2*o2*pow(v1,2)*z1 + 2*n0*n2*o2*pow(v1,2)*z1 - 
   2*j2*n2*s0*pow(v1,2)*z1 - 2*j2*n0*s2*pow(v1,2)*z1 + 
   2*i2*s0*s2*pow(v1,2)*z1 - 2*pow(n1,2)*p2*r0*v2*z1 - 
   4*n0*n1*p2*r1*v2*z1 + 2*k2*n1*r1*s0*v2*z1 + 2*k2*n1*r0*s1*v2*z1 + 
   2*k2*n0*r1*s1*v2*z1 + 2*pow(n1,2)*o2*v0*v2*z1 - 4*j2*n1*s1*v0*v2*z1 + 
   2*i2*pow(s1,2)*v0*v2*z1 + 4*n0*n1*o2*v1*v2*z1 - 4*j2*n1*s0*v1*v2*z1 - 
   4*j2*n0*s1*v1*v2*z1 + 4*i2*s0*s1*v1*v2*z1 + 2*e1*k2*n1*o2*w0*z1 - 
   2*e1*j2*n1*p2*w0*z1 + 2*f1*j2*k2*r1*w0*z1 - 2*f1*i2*p2*r1*w0*z1 + 
   2*m2*n1*p2*r1*w0*z1 - 2*k2*n2*pow(r1,2)*w0*z1 - 4*k2*n1*r1*r2*w0*z1 - 
   2*e1*j2*k2*s1*w0*z1 + 2*e1*i2*p2*s1*w0*z1 + 2*k2*m2*r1*s1*w0*z1 - 
   2*f1*pow(j2,2)*v1*w0*z1 + 2*f1*i2*o2*v1*w0*z1 - 2*m2*n1*o2*v1*w0*z1 + 
   2*j2*n2*r1*v1*w0*z1 + 2*j2*n1*r2*v1*w0*z1 + 2*j2*m2*s1*v1*w0*z1 - 
   2*i2*r2*s1*v1*w0*z1 - 2*i2*r1*s2*v1*w0*z1 + 2*j2*n1*r1*v2*w0*z1 - 
   2*i2*r1*s1*v2*w0*z1 + 2*e1*k2*n0*o2*w1*z1 + 2*e0*k2*n1*o2*w1*z1 - 
   2*e1*j2*n0*p2*w1*z1 - 2*e0*j2*n1*p2*w1*z1 + 2*f1*j2*k2*r0*w1*z1 - 
   2*f1*i2*p2*r0*w1*z1 + 2*m2*n1*p2*r0*w1*z1 + 2*f0*j2*k2*r1*w1*z1 - 
   2*f0*i2*p2*r1*w1*z1 + 2*m2*n0*p2*r1*w1*z1 - 4*k2*n2*r0*r1*w1*z1 - 
   4*k2*n1*r0*r2*w1*z1 - 4*k2*n0*r1*r2*w1*z1 - 2*e1*j2*k2*s0*w1*z1 + 
   2*e1*i2*p2*s0*w1*z1 + 2*k2*m2*r1*s0*w1*z1 - 2*e0*j2*k2*s1*w1*z1 + 
   2*e0*i2*p2*s1*w1*z1 + 2*k2*m2*r0*s1*w1*z1 - 2*f1*pow(j2,2)*v0*w1*z1 + 
   2*f1*i2*o2*v0*w1*z1 - 2*m2*n1*o2*v0*w1*z1 + 2*j2*n2*r1*v0*w1*z1 + 
   2*j2*n1*r2*v0*w1*z1 + 2*j2*m2*s1*v0*w1*z1 - 2*i2*r2*s1*v0*w1*z1 - 
   2*i2*r1*s2*v0*w1*z1 - 2*f0*pow(j2,2)*v1*w1*z1 + 2*f0*i2*o2*v1*w1*z1 - 
   2*m2*n0*o2*v1*w1*z1 + 2*j2*n2*r0*v1*w1*z1 + 2*j2*n0*r2*v1*w1*z1 + 
   2*j2*m2*s0*v1*w1*z1 - 2*i2*r2*s0*v1*w1*z1 - 2*i2*r0*s2*v1*w1*z1 + 
   2*j2*n1*r0*v2*w1*z1 + 2*j2*n0*r1*v2*w1*z1 - 2*i2*r1*s0*v2*w1*z1 - 
   2*i2*r0*s1*v2*w1*z1 + 2*e1*pow(j2,2)*w0*w1*z1 - 2*e1*i2*o2*w0*w1*z1 - 
   4*j2*m2*r1*w0*w1*z1 + 4*i2*r1*r2*w0*w1*z1 + 
   e0*pow(j2,2)*pow(w1,2)*z1 - e0*i2*o2*pow(w1,2)*z1 - 
   2*j2*m2*r0*pow(w1,2)*z1 + 2*i2*r0*r2*pow(w1,2)*z1 - 
   4*k2*n1*r0*r1*w2*z1 - 2*k2*n0*pow(r1,2)*w2*z1 + 2*j2*n1*r1*v0*w2*z1 - 
   2*i2*r1*s1*v0*w2*z1 + 2*j2*n1*r0*v1*w2*z1 + 2*j2*n0*r1*v1*w2*z1 - 
   2*i2*r1*s0*v1*w2*z1 - 2*i2*r0*s1*v1*w2*z1 + 2*i2*pow(r1,2)*w0*w2*z1 + 
   4*i2*r0*r1*w1*w2*z1 + 2*pow(n1,2)*r0*r1*t2*z2 + 
   2*n0*n1*pow(r1,2)*t2*z2 - 2*pow(n1,2)*p2*r1*v0*z2 + 
   2*k2*n1*r1*s1*v0*z2 - 2*pow(n1,2)*p2*r0*v1*z2 - 4*n0*n1*p2*r1*v1*z2 + 
   2*k2*n1*r1*s0*v1*z2 + 2*k2*n1*r0*s1*v1*z2 + 2*k2*n0*r1*s1*v1*z2 + 
   2*pow(n1,2)*o2*v0*v1*z2 - 4*j2*n1*s1*v0*v1*z2 + 
   2*i2*pow(s1,2)*v0*v1*z2 + 2*n0*n1*o2*pow(v1,2)*z2 - 
   2*j2*n1*s0*pow(v1,2)*z2 - 2*j2*n0*s1*pow(v1,2)*z2 + 
   2*i2*s0*s1*pow(v1,2)*z2 - 2*k2*n1*pow(r1,2)*w0*z2 + 
   2*j2*n1*r1*v1*w0*z2 - 2*i2*r1*s1*v1*w0*z2 - 4*k2*n1*r0*r1*w1*z2 - 
   2*k2*n0*pow(r1,2)*w1*z2 + 2*j2*n1*r1*v0*w1*z2 - 2*i2*r1*s1*v0*w1*z2 + 
   2*j2*n1*r0*v1*w1*z2 + 2*j2*n0*r1*v1*w1*z2 - 2*i2*r1*s0*v1*w1*z2 - 
   2*i2*r0*s1*v1*w1*z2 + 2*i2*pow(r1,2)*w0*w1*z2 + 2*i2*r0*r1*pow(w1,2)*z2
;

	k[14]  = -(pow(c1,2)*pow(n1,2)*pow(p2,2)) - 2*c1*d1*k2*n1*p2*r1 - 
   pow(d1,2)*pow(k2,2)*pow(r1,2) + 2*pow(c1,2)*k2*n1*p2*s1 + 
   2*c1*d1*pow(k2,2)*r1*s1 - pow(c1,2)*pow(k2,2)*pow(s1,2) + 
   pow(c1,2)*pow(n1,2)*o2*t2 + 2*c1*d1*j2*n1*r1*t2 - 
   2*c1*pow(n1,2)*q2*r1*t2 + pow(d1,2)*i2*pow(r1,2)*t2 - 
   2*d1*l2*n1*pow(r1,2)*t2 - 2*pow(c1,2)*j2*n1*s1*t2 - 
   2*c1*d1*i2*r1*s1*t2 + 2*c1*l2*n1*r1*s1*t2 + 
   pow(c1,2)*i2*pow(s1,2)*t2 - 2*d1*e1*k2*n1*o2*u1 + 
   2*c1*f1*k2*n1*o2*u1 + 2*d1*e1*j2*n1*p2*u1 - 2*c1*f1*j2*n1*p2*u1 - 
   2*e1*pow(n1,2)*p2*q2*u1 - 2*d1*f1*j2*k2*r1*u1 + 2*c1*g1*j2*k2*r1*u1 + 
   2*d1*f1*i2*p2*r1*u1 - 2*c1*g1*i2*p2*r1*u1 - 2*f1*l2*n1*p2*r1*u1 - 
   2*d1*m2*n1*p2*r1*u1 + 2*c2*pow(n1,2)*p2*r1*u1 + 4*c1*n1*n2*p2*r1*u1 - 
   2*f1*k2*n1*q2*r1*u1 - 2*g1*k2*l2*pow(r1,2)*u1 + 
   2*d2*k2*n1*pow(r1,2)*u1 + 2*d1*k2*n2*pow(r1,2)*u1 + 
   2*c1*pow(n1,2)*p2*r2*u1 + 4*d1*k2*n1*r1*r2*u1 + 2*d1*e1*j2*k2*s1*u1 - 
   2*c1*f1*j2*k2*s1*u1 - 2*d1*e1*i2*p2*s1*u1 + 2*c1*f1*i2*p2*s1*u1 + 
   2*e1*l2*n1*p2*s1*u1 - 2*c1*m2*n1*p2*s1*u1 + 2*e1*k2*n1*q2*s1*u1 + 
   4*f1*k2*l2*r1*s1*u1 - 2*d1*k2*m2*r1*s1*u1 - 2*c2*k2*n1*r1*s1*u1 - 
   2*c1*k2*n2*r1*s1*u1 - 2*c1*k2*n1*r2*s1*u1 - 2*e1*k2*l2*pow(s1,2)*u1 + 
   2*c1*k2*m2*pow(s1,2)*u1 - 2*c1*k2*n1*r1*s2*u1 - 
   pow(f1,2)*pow(j2,2)*pow(u1,2) + e1*g1*pow(j2,2)*pow(u1,2) + 
   pow(f1,2)*i2*o2*pow(u1,2) - e1*g1*i2*o2*pow(u1,2) - 
   2*f1*m2*n1*o2*pow(u1,2) + e2*pow(n1,2)*o2*pow(u1,2) + 
   2*e1*n1*n2*o2*pow(u1,2) - 2*g1*j2*m2*r1*pow(u1,2) + 
   2*f2*j2*n1*r1*pow(u1,2) + 2*f1*j2*n2*r1*pow(u1,2) + 
   g2*i2*pow(r1,2)*pow(u1,2) - pow(n2,2)*pow(r1,2)*pow(u1,2) + 
   2*f1*j2*n1*r2*pow(u1,2) + 2*g1*i2*r1*r2*pow(u1,2) - 
   4*n1*n2*r1*r2*pow(u1,2) - pow(n1,2)*pow(r2,2)*pow(u1,2) + 
   2*f1*j2*m2*s1*pow(u1,2) - 2*e2*j2*n1*s1*pow(u1,2) - 
   2*e1*j2*n2*s1*pow(u1,2) - 2*f2*i2*r1*s1*pow(u1,2) + 
   2*m2*n2*r1*s1*pow(u1,2) - 2*f1*i2*r2*s1*pow(u1,2) + 
   2*m2*n1*r2*s1*pow(u1,2) + e2*i2*pow(s1,2)*pow(u1,2) - 
   pow(m2,2)*pow(s1,2)*pow(u1,2) - 2*e1*j2*n1*s2*pow(u1,2) - 
   2*f1*i2*r1*s2*pow(u1,2) + 2*m2*n1*r1*s2*pow(u1,2) + 
   2*e1*i2*s1*s2*pow(u1,2) + 2*c1*pow(n1,2)*p2*r1*u2 + 
   2*d1*k2*n1*pow(r1,2)*u2 - 2*c1*k2*n1*r1*s1*u2 + 
   2*e1*pow(n1,2)*o2*u1*u2 + 4*f1*j2*n1*r1*u1*u2 + 
   2*g1*i2*pow(r1,2)*u1*u2 - 4*n1*n2*pow(r1,2)*u1*u2 - 
   4*pow(n1,2)*r1*r2*u1*u2 - 4*e1*j2*n1*s1*u1*u2 - 4*f1*i2*r1*s1*u1*u2 + 
   4*m2*n1*r1*s1*u1*u2 + 2*e1*i2*pow(s1,2)*u1*u2 - 
   pow(n1,2)*pow(r1,2)*pow(u2,2) + 2*c1*d1*k2*n1*o2*v1 - 
   2*c1*d1*j2*n1*p2*v1 + 2*c1*pow(n1,2)*p2*q2*v1 + 
   2*pow(d1,2)*j2*k2*r1*v1 - 2*pow(d1,2)*i2*p2*r1*v1 + 
   4*d1*l2*n1*p2*r1*v1 - 2*d1*k2*n1*q2*r1*v1 - 2*c1*d1*j2*k2*s1*v1 + 
   2*c1*d1*i2*p2*s1*v1 - 2*c1*l2*n1*p2*s1*v1 - 2*c1*k2*n1*q2*s1*v1 - 
   2*d1*k2*l2*r1*s1*v1 + 2*c1*k2*l2*pow(s1,2)*v1 + 
   2*d1*f1*pow(j2,2)*u1*v1 - 2*c1*g1*pow(j2,2)*u1*v1 - 
   2*d1*f1*i2*o2*u1*v1 + 2*c1*g1*i2*o2*u1*v1 + 2*f1*l2*n1*o2*u1*v1 + 
   2*d1*m2*n1*o2*u1*v1 - 2*c2*pow(n1,2)*o2*u1*v1 - 4*c1*n1*n2*o2*u1*v1 - 
   2*f1*j2*n1*q2*u1*v1 + 2*g1*j2*l2*r1*u1*v1 - 2*d2*j2*n1*r1*u1*v1 - 
   2*d1*j2*n2*r1*u1*v1 - 2*g1*i2*q2*r1*u1*v1 + 4*n1*n2*q2*r1*u1*v1 - 
   2*d1*j2*n1*r2*u1*v1 + 2*pow(n1,2)*q2*r2*u1*v1 - 2*f1*j2*l2*s1*u1*v1 - 
   2*d1*j2*m2*s1*u1*v1 + 4*c2*j2*n1*s1*u1*v1 + 4*c1*j2*n2*s1*u1*v1 + 
   2*f1*i2*q2*s1*u1*v1 - 2*m2*n1*q2*s1*u1*v1 + 2*d2*i2*r1*s1*u1*v1 - 
   2*l2*n2*r1*s1*u1*v1 + 2*d1*i2*r2*s1*u1*v1 - 2*l2*n1*r2*s1*u1*v1 - 
   2*c2*i2*pow(s1,2)*u1*v1 + 2*l2*m2*pow(s1,2)*u1*v1 + 
   4*c1*j2*n1*s2*u1*v1 + 2*d1*i2*r1*s2*u1*v1 - 2*l2*n1*r1*s2*u1*v1 - 
   4*c1*i2*s1*s2*u1*v1 - 2*c1*pow(n1,2)*o2*u2*v1 - 2*d1*j2*n1*r1*u2*v1 + 
   2*pow(n1,2)*q2*r1*u2*v1 + 4*c1*j2*n1*s1*u2*v1 + 2*d1*i2*r1*s1*u2*v1 - 
   2*l2*n1*r1*s1*u2*v1 - 2*c1*i2*pow(s1,2)*u2*v1 - 
   pow(d1,2)*pow(j2,2)*pow(v1,2) + pow(d1,2)*i2*o2*pow(v1,2) - 
   2*d1*l2*n1*o2*pow(v1,2) + 2*d1*j2*n1*q2*pow(v1,2) - 
   pow(n1,2)*pow(q2,2)*pow(v1,2) + 2*d1*j2*l2*s1*pow(v1,2) - 
   2*d1*i2*q2*s1*pow(v1,2) + 2*l2*n1*q2*s1*pow(v1,2) - 
   pow(l2,2)*pow(s1,2)*pow(v1,2) - 2*c1*pow(n1,2)*o2*u1*v2 - 
   2*d1*j2*n1*r1*u1*v2 + 2*pow(n1,2)*q2*r1*u1*v2 + 4*c1*j2*n1*s1*u1*v2 + 
   2*d1*i2*r1*s1*u1*v2 - 2*l2*n1*r1*s1*u1*v2 - 2*c1*i2*pow(s1,2)*u1*v2 - 
   2*pow(c1,2)*k2*n1*o2*w1 + 2*pow(c1,2)*j2*n1*p2*w1 - 
   2*c1*d1*j2*k2*r1*w1 + 2*c1*d1*i2*p2*r1*w1 - 2*c1*l2*n1*p2*r1*w1 + 
   4*c1*k2*n1*q2*r1*w1 + 2*d1*k2*l2*pow(r1,2)*w1 + 
   2*pow(c1,2)*j2*k2*s1*w1 - 2*pow(c1,2)*i2*p2*s1*w1 - 
   2*c1*k2*l2*r1*s1*w1 - 2*d1*e1*pow(j2,2)*u1*w1 + 
   2*c1*f1*pow(j2,2)*u1*w1 + 2*d1*e1*i2*o2*u1*w1 - 2*c1*f1*i2*o2*u1*w1 - 
   2*e1*l2*n1*o2*u1*w1 + 2*c1*m2*n1*o2*u1*w1 + 2*e1*j2*n1*q2*u1*w1 - 
   2*f1*j2*l2*r1*u1*w1 + 4*d1*j2*m2*r1*u1*w1 - 2*c2*j2*n1*r1*u1*w1 - 
   2*c1*j2*n2*r1*u1*w1 + 2*f1*i2*q2*r1*u1*w1 - 2*m2*n1*q2*r1*u1*w1 - 
   2*d2*i2*pow(r1,2)*u1*w1 + 2*l2*n2*pow(r1,2)*u1*w1 - 
   2*c1*j2*n1*r2*u1*w1 - 4*d1*i2*r1*r2*u1*w1 + 4*l2*n1*r1*r2*u1*w1 + 
   2*e1*j2*l2*s1*u1*w1 - 2*c1*j2*m2*s1*u1*w1 - 2*e1*i2*q2*s1*u1*w1 + 
   2*c2*i2*r1*s1*u1*w1 - 2*l2*m2*r1*s1*u1*w1 + 2*c1*i2*r2*s1*u1*w1 + 
   2*c1*i2*r1*s2*u1*w1 - 2*c1*j2*n1*r1*u2*w1 - 2*d1*i2*pow(r1,2)*u2*w1 + 
   2*l2*n1*pow(r1,2)*u2*w1 + 2*c1*i2*r1*s1*u2*w1 + 
   2*c1*d1*pow(j2,2)*v1*w1 - 2*c1*d1*i2*o2*v1*w1 + 2*c1*l2*n1*o2*v1*w1 - 
   2*c1*j2*n1*q2*v1*w1 - 2*d1*j2*l2*r1*v1*w1 + 2*d1*i2*q2*r1*v1*w1 - 
   2*l2*n1*q2*r1*v1*w1 - 2*c1*j2*l2*s1*v1*w1 + 2*c1*i2*q2*s1*v1*w1 + 
   2*pow(l2,2)*r1*s1*v1*w1 - pow(c1,2)*pow(j2,2)*pow(w1,2) + 
   pow(c1,2)*i2*o2*pow(w1,2) + 2*c1*j2*l2*r1*pow(w1,2) - 
   2*c1*i2*q2*r1*pow(w1,2) - pow(l2,2)*pow(r1,2)*pow(w1,2) - 
   2*c1*j2*n1*r1*u1*w2 - 2*d1*i2*pow(r1,2)*u1*w2 + 
   2*l2*n1*pow(r1,2)*u1*w2 + 2*c1*i2*r1*s1*u1*w2 + 
   e1*pow(n1,2)*pow(p2,2)*z1 + 2*f1*k2*n1*p2*r1*z1 + 
   g1*pow(k2,2)*pow(r1,2)*z1 - 2*e1*k2*n1*p2*s1*z1 - 
   2*f1*pow(k2,2)*r1*s1*z1 + e1*pow(k2,2)*pow(s1,2)*z1 - 
   e1*pow(n1,2)*o2*t2*z1 - 2*f1*j2*n1*r1*t2*z1 - g1*i2*pow(r1,2)*t2*z1 + 
   2*n1*n2*pow(r1,2)*t2*z1 + 2*pow(n1,2)*r1*r2*t2*z1 + 
   2*e1*j2*n1*s1*t2*z1 + 2*f1*i2*r1*s1*t2*z1 - 2*m2*n1*r1*s1*t2*z1 - 
   e1*i2*pow(s1,2)*t2*z1 - 2*f1*k2*n1*o2*v1*z1 + 2*f1*j2*n1*p2*v1*z1 - 
   2*g1*j2*k2*r1*v1*z1 + 2*g1*i2*p2*r1*v1*z1 - 4*n1*n2*p2*r1*v1*z1 - 
   2*pow(n1,2)*p2*r2*v1*z1 + 2*f1*j2*k2*s1*v1*z1 - 2*f1*i2*p2*s1*v1*z1 + 
   2*m2*n1*p2*s1*v1*z1 + 2*k2*n2*r1*s1*v1*z1 + 2*k2*n1*r2*s1*v1*z1 - 
   2*k2*m2*pow(s1,2)*v1*z1 + 2*k2*n1*r1*s2*v1*z1 + 
   g1*pow(j2,2)*pow(v1,2)*z1 - g1*i2*o2*pow(v1,2)*z1 + 
   2*n1*n2*o2*pow(v1,2)*z1 - 2*j2*n2*s1*pow(v1,2)*z1 - 
   2*j2*n1*s2*pow(v1,2)*z1 + 2*i2*s1*s2*pow(v1,2)*z1 - 
   2*pow(n1,2)*p2*r1*v2*z1 + 2*k2*n1*r1*s1*v2*z1 + 
   2*pow(n1,2)*o2*v1*v2*z1 - 4*j2*n1*s1*v1*v2*z1 + 
   2*i2*pow(s1,2)*v1*v2*z1 + 2*e1*k2*n1*o2*w1*z1 - 2*e1*j2*n1*p2*w1*z1 + 
   2*f1*j2*k2*r1*w1*z1 - 2*f1*i2*p2*r1*w1*z1 + 2*m2*n1*p2*r1*w1*z1 - 
   2*k2*n2*pow(r1,2)*w1*z1 - 4*k2*n1*r1*r2*w1*z1 - 2*e1*j2*k2*s1*w1*z1 + 
   2*e1*i2*p2*s1*w1*z1 + 2*k2*m2*r1*s1*w1*z1 - 2*f1*pow(j2,2)*v1*w1*z1 + 
   2*f1*i2*o2*v1*w1*z1 - 2*m2*n1*o2*v1*w1*z1 + 2*j2*n2*r1*v1*w1*z1 + 
   2*j2*n1*r2*v1*w1*z1 + 2*j2*m2*s1*v1*w1*z1 - 2*i2*r2*s1*v1*w1*z1 - 
   2*i2*r1*s2*v1*w1*z1 + 2*j2*n1*r1*v2*w1*z1 - 2*i2*r1*s1*v2*w1*z1 + 
   e1*pow(j2,2)*pow(w1,2)*z1 - e1*i2*o2*pow(w1,2)*z1 - 
   2*j2*m2*r1*pow(w1,2)*z1 + 2*i2*r1*r2*pow(w1,2)*z1 - 
   2*k2*n1*pow(r1,2)*w2*z1 + 2*j2*n1*r1*v1*w2*z1 - 2*i2*r1*s1*v1*w2*z1 + 
   2*i2*pow(r1,2)*w1*w2*z1 + pow(n1,2)*pow(r1,2)*t2*z2 - 
   2*pow(n1,2)*p2*r1*v1*z2 + 2*k2*n1*r1*s1*v1*z2 + 
   pow(n1,2)*o2*pow(v1,2)*z2 - 2*j2*n1*s1*pow(v1,2)*z2 + 
   i2*pow(s1,2)*pow(v1,2)*z2 - 2*k2*n1*pow(r1,2)*w1*z2 + 
   2*j2*n1*r1*v1*w1*z2 - 2*i2*r1*s1*v1*w1*z2 + i2*pow(r1,2)*pow(w1,2)*z2
;

	k[15]  = 2*c0*pow(n0,2)*p2*r0*u0 + 2*d0*k2*n0*pow(r0,2)*u0 - 2*c0*k2*n0*r0*s0*u0 + 
   e0*pow(n0,2)*o2*pow(u0,2) + 2*f0*j2*n0*r0*pow(u0,2) + 
   g0*i2*pow(r0,2)*pow(u0,2) - 2*n0*n2*pow(r0,2)*pow(u0,2) - 
   2*pow(n0,2)*r0*r2*pow(u0,2) - 2*e0*j2*n0*s0*pow(u0,2) - 
   2*f0*i2*r0*s0*pow(u0,2) + 2*m2*n0*r0*s0*pow(u0,2) + 
   e0*i2*pow(s0,2)*pow(u0,2) - 2*pow(n0,2)*pow(r0,2)*u0*u2 - 
   2*c0*pow(n0,2)*o2*u0*v0 - 2*d0*j2*n0*r0*u0*v0 + 
   2*pow(n0,2)*q2*r0*u0*v0 + 4*c0*j2*n0*s0*u0*v0 + 2*d0*i2*r0*s0*u0*v0 - 
   2*l2*n0*r0*s0*u0*v0 - 2*c0*i2*pow(s0,2)*u0*v0 - 2*c0*j2*n0*r0*u0*w0 - 
   2*d0*i2*pow(r0,2)*u0*w0 + 2*l2*n0*pow(r0,2)*u0*w0 + 
   2*c0*i2*r0*s0*u0*w0 + pow(n0,2)*pow(r0,2)*t2*z0 - 
   2*pow(n0,2)*p2*r0*v0*z0 + 2*k2*n0*r0*s0*v0*z0 + 
   pow(n0,2)*o2*pow(v0,2)*z0 - 2*j2*n0*s0*pow(v0,2)*z0 + 
   i2*pow(s0,2)*pow(v0,2)*z0 - 2*k2*n0*pow(r0,2)*w0*z0 + 
   2*j2*n0*r0*v0*w0*z0 - 2*i2*r0*s0*v0*w0*z0 + i2*pow(r0,2)*pow(w0,2)*z0
;

	k[16]  = 2*c1*pow(n0,2)*p2*r0*u0 + 4*c0*n0*n1*p2*r0*u0 + 2*d1*k2*n0*pow(r0,2)*u0 + 
   2*d0*k2*n1*pow(r0,2)*u0 + 2*c0*pow(n0,2)*p2*r1*u0 + 
   4*d0*k2*n0*r0*r1*u0 - 2*c1*k2*n0*r0*s0*u0 - 2*c0*k2*n1*r0*s0*u0 - 
   2*c0*k2*n0*r1*s0*u0 - 2*c0*k2*n0*r0*s1*u0 + 
   e1*pow(n0,2)*o2*pow(u0,2) + 2*e0*n0*n1*o2*pow(u0,2) + 
   2*f1*j2*n0*r0*pow(u0,2) + 2*f0*j2*n1*r0*pow(u0,2) + 
   g1*i2*pow(r0,2)*pow(u0,2) - 2*n1*n2*pow(r0,2)*pow(u0,2) + 
   2*f0*j2*n0*r1*pow(u0,2) + 2*g0*i2*r0*r1*pow(u0,2) - 
   4*n0*n2*r0*r1*pow(u0,2) - 4*n0*n1*r0*r2*pow(u0,2) - 
   2*pow(n0,2)*r1*r2*pow(u0,2) - 2*e1*j2*n0*s0*pow(u0,2) - 
   2*e0*j2*n1*s0*pow(u0,2) - 2*f1*i2*r0*s0*pow(u0,2) + 
   2*m2*n1*r0*s0*pow(u0,2) - 2*f0*i2*r1*s0*pow(u0,2) + 
   2*m2*n0*r1*s0*pow(u0,2) + e1*i2*pow(s0,2)*pow(u0,2) - 
   2*e0*j2*n0*s1*pow(u0,2) - 2*f0*i2*r0*s1*pow(u0,2) + 
   2*m2*n0*r0*s1*pow(u0,2) + 2*e0*i2*s0*s1*pow(u0,2) + 
   2*c0*pow(n0,2)*p2*r0*u1 + 2*d0*k2*n0*pow(r0,2)*u1 - 
   2*c0*k2*n0*r0*s0*u1 + 2*e0*pow(n0,2)*o2*u0*u1 + 4*f0*j2*n0*r0*u0*u1 + 
   2*g0*i2*pow(r0,2)*u0*u1 - 4*n0*n2*pow(r0,2)*u0*u1 - 
   4*pow(n0,2)*r0*r2*u0*u1 - 4*e0*j2*n0*s0*u0*u1 - 4*f0*i2*r0*s0*u0*u1 + 
   4*m2*n0*r0*s0*u0*u1 + 2*e0*i2*pow(s0,2)*u0*u1 - 
   4*n0*n1*pow(r0,2)*u0*u2 - 4*pow(n0,2)*r0*r1*u0*u2 - 
   2*pow(n0,2)*pow(r0,2)*u1*u2 - 2*c1*pow(n0,2)*o2*u0*v0 - 
   4*c0*n0*n1*o2*u0*v0 - 2*d1*j2*n0*r0*u0*v0 - 2*d0*j2*n1*r0*u0*v0 + 
   4*n0*n1*q2*r0*u0*v0 - 2*d0*j2*n0*r1*u0*v0 + 2*pow(n0,2)*q2*r1*u0*v0 + 
   4*c1*j2*n0*s0*u0*v0 + 4*c0*j2*n1*s0*u0*v0 + 2*d1*i2*r0*s0*u0*v0 - 
   2*l2*n1*r0*s0*u0*v0 + 2*d0*i2*r1*s0*u0*v0 - 2*l2*n0*r1*s0*u0*v0 - 
   2*c1*i2*pow(s0,2)*u0*v0 + 4*c0*j2*n0*s1*u0*v0 + 2*d0*i2*r0*s1*u0*v0 - 
   2*l2*n0*r0*s1*u0*v0 - 4*c0*i2*s0*s1*u0*v0 - 2*c0*pow(n0,2)*o2*u1*v0 - 
   2*d0*j2*n0*r0*u1*v0 + 2*pow(n0,2)*q2*r0*u1*v0 + 4*c0*j2*n0*s0*u1*v0 + 
   2*d0*i2*r0*s0*u1*v0 - 2*l2*n0*r0*s0*u1*v0 - 2*c0*i2*pow(s0,2)*u1*v0 - 
   2*c0*pow(n0,2)*o2*u0*v1 - 2*d0*j2*n0*r0*u0*v1 + 
   2*pow(n0,2)*q2*r0*u0*v1 + 4*c0*j2*n0*s0*u0*v1 + 2*d0*i2*r0*s0*u0*v1 - 
   2*l2*n0*r0*s0*u0*v1 - 2*c0*i2*pow(s0,2)*u0*v1 - 2*c1*j2*n0*r0*u0*w0 - 
   2*c0*j2*n1*r0*u0*w0 - 2*d1*i2*pow(r0,2)*u0*w0 + 
   2*l2*n1*pow(r0,2)*u0*w0 - 2*c0*j2*n0*r1*u0*w0 - 4*d0*i2*r0*r1*u0*w0 + 
   4*l2*n0*r0*r1*u0*w0 + 2*c1*i2*r0*s0*u0*w0 + 2*c0*i2*r1*s0*u0*w0 + 
   2*c0*i2*r0*s1*u0*w0 - 2*c0*j2*n0*r0*u1*w0 - 2*d0*i2*pow(r0,2)*u1*w0 + 
   2*l2*n0*pow(r0,2)*u1*w0 + 2*c0*i2*r0*s0*u1*w0 - 2*c0*j2*n0*r0*u0*w1 - 
   2*d0*i2*pow(r0,2)*u0*w1 + 2*l2*n0*pow(r0,2)*u0*w1 + 
   2*c0*i2*r0*s0*u0*w1 + 2*n0*n1*pow(r0,2)*t2*z0 + 
   2*pow(n0,2)*r0*r1*t2*z0 - 4*n0*n1*p2*r0*v0*z0 - 
   2*pow(n0,2)*p2*r1*v0*z0 + 2*k2*n1*r0*s0*v0*z0 + 2*k2*n0*r1*s0*v0*z0 + 
   2*k2*n0*r0*s1*v0*z0 + 2*n0*n1*o2*pow(v0,2)*z0 - 
   2*j2*n1*s0*pow(v0,2)*z0 - 2*j2*n0*s1*pow(v0,2)*z0 + 
   2*i2*s0*s1*pow(v0,2)*z0 - 2*pow(n0,2)*p2*r0*v1*z0 + 
   2*k2*n0*r0*s0*v1*z0 + 2*pow(n0,2)*o2*v0*v1*z0 - 4*j2*n0*s0*v0*v1*z0 + 
   2*i2*pow(s0,2)*v0*v1*z0 - 2*k2*n1*pow(r0,2)*w0*z0 - 
   4*k2*n0*r0*r1*w0*z0 + 2*j2*n1*r0*v0*w0*z0 + 2*j2*n0*r1*v0*w0*z0 - 
   2*i2*r1*s0*v0*w0*z0 - 2*i2*r0*s1*v0*w0*z0 + 2*j2*n0*r0*v1*w0*z0 - 
   2*i2*r0*s0*v1*w0*z0 + 2*i2*r0*r1*pow(w0,2)*z0 - 
   2*k2*n0*pow(r0,2)*w1*z0 + 2*j2*n0*r0*v0*w1*z0 - 2*i2*r0*s0*v0*w1*z0 + 
   2*i2*pow(r0,2)*w0*w1*z0 + pow(n0,2)*pow(r0,2)*t2*z1 - 
   2*pow(n0,2)*p2*r0*v0*z1 + 2*k2*n0*r0*s0*v0*z1 + 
   pow(n0,2)*o2*pow(v0,2)*z1 - 2*j2*n0*s0*pow(v0,2)*z1 + 
   i2*pow(s0,2)*pow(v0,2)*z1 - 2*k2*n0*pow(r0,2)*w0*z1 + 
   2*j2*n0*r0*v0*w0*z1 - 2*i2*r0*s0*v0*w0*z1 + i2*pow(r0,2)*pow(w0,2)*z1
;

	k[17]  = 4*c1*n0*n1*p2*r0*u0 + 2*c0*pow(n1,2)*p2*r0*u0 + 2*d1*k2*n1*pow(r0,2)*u0 + 
   2*c1*pow(n0,2)*p2*r1*u0 + 4*c0*n0*n1*p2*r1*u0 + 4*d1*k2*n0*r0*r1*u0 + 
   4*d0*k2*n1*r0*r1*u0 + 2*d0*k2*n0*pow(r1,2)*u0 - 2*c1*k2*n1*r0*s0*u0 - 
   2*c1*k2*n0*r1*s0*u0 - 2*c0*k2*n1*r1*s0*u0 - 2*c1*k2*n0*r0*s1*u0 - 
   2*c0*k2*n1*r0*s1*u0 - 2*c0*k2*n0*r1*s1*u0 + 2*e1*n0*n1*o2*pow(u0,2) + 
   e0*pow(n1,2)*o2*pow(u0,2) + 2*f1*j2*n1*r0*pow(u0,2) + 
   2*f1*j2*n0*r1*pow(u0,2) + 2*f0*j2*n1*r1*pow(u0,2) + 
   2*g1*i2*r0*r1*pow(u0,2) - 4*n1*n2*r0*r1*pow(u0,2) + 
   g0*i2*pow(r1,2)*pow(u0,2) - 2*n0*n2*pow(r1,2)*pow(u0,2) - 
   2*pow(n1,2)*r0*r2*pow(u0,2) - 4*n0*n1*r1*r2*pow(u0,2) - 
   2*e1*j2*n1*s0*pow(u0,2) - 2*f1*i2*r1*s0*pow(u0,2) + 
   2*m2*n1*r1*s0*pow(u0,2) - 2*e1*j2*n0*s1*pow(u0,2) - 
   2*e0*j2*n1*s1*pow(u0,2) - 2*f1*i2*r0*s1*pow(u0,2) + 
   2*m2*n1*r0*s1*pow(u0,2) - 2*f0*i2*r1*s1*pow(u0,2) + 
   2*m2*n0*r1*s1*pow(u0,2) + 2*e1*i2*s0*s1*pow(u0,2) + 
   e0*i2*pow(s1,2)*pow(u0,2) + 2*c1*pow(n0,2)*p2*r0*u1 + 
   4*c0*n0*n1*p2*r0*u1 + 2*d1*k2*n0*pow(r0,2)*u1 + 
   2*d0*k2*n1*pow(r0,2)*u1 + 2*c0*pow(n0,2)*p2*r1*u1 + 
   4*d0*k2*n0*r0*r1*u1 - 2*c1*k2*n0*r0*s0*u1 - 2*c0*k2*n1*r0*s0*u1 - 
   2*c0*k2*n0*r1*s0*u1 - 2*c0*k2*n0*r0*s1*u1 + 2*e1*pow(n0,2)*o2*u0*u1 + 
   4*e0*n0*n1*o2*u0*u1 + 4*f1*j2*n0*r0*u0*u1 + 4*f0*j2*n1*r0*u0*u1 + 
   2*g1*i2*pow(r0,2)*u0*u1 - 4*n1*n2*pow(r0,2)*u0*u1 + 
   4*f0*j2*n0*r1*u0*u1 + 4*g0*i2*r0*r1*u0*u1 - 8*n0*n2*r0*r1*u0*u1 - 
   8*n0*n1*r0*r2*u0*u1 - 4*pow(n0,2)*r1*r2*u0*u1 - 4*e1*j2*n0*s0*u0*u1 - 
   4*e0*j2*n1*s0*u0*u1 - 4*f1*i2*r0*s0*u0*u1 + 4*m2*n1*r0*s0*u0*u1 - 
   4*f0*i2*r1*s0*u0*u1 + 4*m2*n0*r1*s0*u0*u1 + 2*e1*i2*pow(s0,2)*u0*u1 - 
   4*e0*j2*n0*s1*u0*u1 - 4*f0*i2*r0*s1*u0*u1 + 4*m2*n0*r0*s1*u0*u1 + 
   4*e0*i2*s0*s1*u0*u1 + e0*pow(n0,2)*o2*pow(u1,2) + 
   2*f0*j2*n0*r0*pow(u1,2) + g0*i2*pow(r0,2)*pow(u1,2) - 
   2*n0*n2*pow(r0,2)*pow(u1,2) - 2*pow(n0,2)*r0*r2*pow(u1,2) - 
   2*e0*j2*n0*s0*pow(u1,2) - 2*f0*i2*r0*s0*pow(u1,2) + 
   2*m2*n0*r0*s0*pow(u1,2) + e0*i2*pow(s0,2)*pow(u1,2) - 
   2*pow(n1,2)*pow(r0,2)*u0*u2 - 8*n0*n1*r0*r1*u0*u2 - 
   2*pow(n0,2)*pow(r1,2)*u0*u2 - 4*n0*n1*pow(r0,2)*u1*u2 - 
   4*pow(n0,2)*r0*r1*u1*u2 - 4*c1*n0*n1*o2*u0*v0 - 
   2*c0*pow(n1,2)*o2*u0*v0 - 2*d1*j2*n1*r0*u0*v0 + 
   2*pow(n1,2)*q2*r0*u0*v0 - 2*d1*j2*n0*r1*u0*v0 - 2*d0*j2*n1*r1*u0*v0 + 
   4*n0*n1*q2*r1*u0*v0 + 4*c1*j2*n1*s0*u0*v0 + 2*d1*i2*r1*s0*u0*v0 - 
   2*l2*n1*r1*s0*u0*v0 + 4*c1*j2*n0*s1*u0*v0 + 4*c0*j2*n1*s1*u0*v0 + 
   2*d1*i2*r0*s1*u0*v0 - 2*l2*n1*r0*s1*u0*v0 + 2*d0*i2*r1*s1*u0*v0 - 
   2*l2*n0*r1*s1*u0*v0 - 4*c1*i2*s0*s1*u0*v0 - 2*c0*i2*pow(s1,2)*u0*v0 - 
   2*c1*pow(n0,2)*o2*u1*v0 - 4*c0*n0*n1*o2*u1*v0 - 2*d1*j2*n0*r0*u1*v0 - 
   2*d0*j2*n1*r0*u1*v0 + 4*n0*n1*q2*r0*u1*v0 - 2*d0*j2*n0*r1*u1*v0 + 
   2*pow(n0,2)*q2*r1*u1*v0 + 4*c1*j2*n0*s0*u1*v0 + 4*c0*j2*n1*s0*u1*v0 + 
   2*d1*i2*r0*s0*u1*v0 - 2*l2*n1*r0*s0*u1*v0 + 2*d0*i2*r1*s0*u1*v0 - 
   2*l2*n0*r1*s0*u1*v0 - 2*c1*i2*pow(s0,2)*u1*v0 + 4*c0*j2*n0*s1*u1*v0 + 
   2*d0*i2*r0*s1*u1*v0 - 2*l2*n0*r0*s1*u1*v0 - 4*c0*i2*s0*s1*u1*v0 - 
   2*c1*pow(n0,2)*o2*u0*v1 - 4*c0*n0*n1*o2*u0*v1 - 2*d1*j2*n0*r0*u0*v1 - 
   2*d0*j2*n1*r0*u0*v1 + 4*n0*n1*q2*r0*u0*v1 - 2*d0*j2*n0*r1*u0*v1 + 
   2*pow(n0,2)*q2*r1*u0*v1 + 4*c1*j2*n0*s0*u0*v1 + 4*c0*j2*n1*s0*u0*v1 + 
   2*d1*i2*r0*s0*u0*v1 - 2*l2*n1*r0*s0*u0*v1 + 2*d0*i2*r1*s0*u0*v1 - 
   2*l2*n0*r1*s0*u0*v1 - 2*c1*i2*pow(s0,2)*u0*v1 + 4*c0*j2*n0*s1*u0*v1 + 
   2*d0*i2*r0*s1*u0*v1 - 2*l2*n0*r0*s1*u0*v1 - 4*c0*i2*s0*s1*u0*v1 - 
   2*c0*pow(n0,2)*o2*u1*v1 - 2*d0*j2*n0*r0*u1*v1 + 
   2*pow(n0,2)*q2*r0*u1*v1 + 4*c0*j2*n0*s0*u1*v1 + 2*d0*i2*r0*s0*u1*v1 - 
   2*l2*n0*r0*s0*u1*v1 - 2*c0*i2*pow(s0,2)*u1*v1 - 2*c1*j2*n1*r0*u0*w0 - 
   2*c1*j2*n0*r1*u0*w0 - 2*c0*j2*n1*r1*u0*w0 - 4*d1*i2*r0*r1*u0*w0 + 
   4*l2*n1*r0*r1*u0*w0 - 2*d0*i2*pow(r1,2)*u0*w0 + 
   2*l2*n0*pow(r1,2)*u0*w0 + 2*c1*i2*r1*s0*u0*w0 + 2*c1*i2*r0*s1*u0*w0 + 
   2*c0*i2*r1*s1*u0*w0 - 2*c1*j2*n0*r0*u1*w0 - 2*c0*j2*n1*r0*u1*w0 - 
   2*d1*i2*pow(r0,2)*u1*w0 + 2*l2*n1*pow(r0,2)*u1*w0 - 
   2*c0*j2*n0*r1*u1*w0 - 4*d0*i2*r0*r1*u1*w0 + 4*l2*n0*r0*r1*u1*w0 + 
   2*c1*i2*r0*s0*u1*w0 + 2*c0*i2*r1*s0*u1*w0 + 2*c0*i2*r0*s1*u1*w0 - 
   2*c1*j2*n0*r0*u0*w1 - 2*c0*j2*n1*r0*u0*w1 - 2*d1*i2*pow(r0,2)*u0*w1 + 
   2*l2*n1*pow(r0,2)*u0*w1 - 2*c0*j2*n0*r1*u0*w1 - 4*d0*i2*r0*r1*u0*w1 + 
   4*l2*n0*r0*r1*u0*w1 + 2*c1*i2*r0*s0*u0*w1 + 2*c0*i2*r1*s0*u0*w1 + 
   2*c0*i2*r0*s1*u0*w1 - 2*c0*j2*n0*r0*u1*w1 - 2*d0*i2*pow(r0,2)*u1*w1 + 
   2*l2*n0*pow(r0,2)*u1*w1 + 2*c0*i2*r0*s0*u1*w1 + 
   pow(n1,2)*pow(r0,2)*t2*z0 + 4*n0*n1*r0*r1*t2*z0 + 
   pow(n0,2)*pow(r1,2)*t2*z0 - 2*pow(n1,2)*p2*r0*v0*z0 - 
   4*n0*n1*p2*r1*v0*z0 + 2*k2*n1*r1*s0*v0*z0 + 2*k2*n1*r0*s1*v0*z0 + 
   2*k2*n0*r1*s1*v0*z0 + pow(n1,2)*o2*pow(v0,2)*z0 - 
   2*j2*n1*s1*pow(v0,2)*z0 + i2*pow(s1,2)*pow(v0,2)*z0 - 
   4*n0*n1*p2*r0*v1*z0 - 2*pow(n0,2)*p2*r1*v1*z0 + 2*k2*n1*r0*s0*v1*z0 + 
   2*k2*n0*r1*s0*v1*z0 + 2*k2*n0*r0*s1*v1*z0 + 4*n0*n1*o2*v0*v1*z0 - 
   4*j2*n1*s0*v0*v1*z0 - 4*j2*n0*s1*v0*v1*z0 + 4*i2*s0*s1*v0*v1*z0 + 
   pow(n0,2)*o2*pow(v1,2)*z0 - 2*j2*n0*s0*pow(v1,2)*z0 + 
   i2*pow(s0,2)*pow(v1,2)*z0 - 4*k2*n1*r0*r1*w0*z0 - 
   2*k2*n0*pow(r1,2)*w0*z0 + 2*j2*n1*r1*v0*w0*z0 - 2*i2*r1*s1*v0*w0*z0 + 
   2*j2*n1*r0*v1*w0*z0 + 2*j2*n0*r1*v1*w0*z0 - 2*i2*r1*s0*v1*w0*z0 - 
   2*i2*r0*s1*v1*w0*z0 + i2*pow(r1,2)*pow(w0,2)*z0 - 
   2*k2*n1*pow(r0,2)*w1*z0 - 4*k2*n0*r0*r1*w1*z0 + 2*j2*n1*r0*v0*w1*z0 + 
   2*j2*n0*r1*v0*w1*z0 - 2*i2*r1*s0*v0*w1*z0 - 2*i2*r0*s1*v0*w1*z0 + 
   2*j2*n0*r0*v1*w1*z0 - 2*i2*r0*s0*v1*w1*z0 + 4*i2*r0*r1*w0*w1*z0 + 
   i2*pow(r0,2)*pow(w1,2)*z0 + 2*n0*n1*pow(r0,2)*t2*z1 + 
   2*pow(n0,2)*r0*r1*t2*z1 - 4*n0*n1*p2*r0*v0*z1 - 
   2*pow(n0,2)*p2*r1*v0*z1 + 2*k2*n1*r0*s0*v0*z1 + 2*k2*n0*r1*s0*v0*z1 + 
   2*k2*n0*r0*s1*v0*z1 + 2*n0*n1*o2*pow(v0,2)*z1 - 
   2*j2*n1*s0*pow(v0,2)*z1 - 2*j2*n0*s1*pow(v0,2)*z1 + 
   2*i2*s0*s1*pow(v0,2)*z1 - 2*pow(n0,2)*p2*r0*v1*z1 + 
   2*k2*n0*r0*s0*v1*z1 + 2*pow(n0,2)*o2*v0*v1*z1 - 4*j2*n0*s0*v0*v1*z1 + 
   2*i2*pow(s0,2)*v0*v1*z1 - 2*k2*n1*pow(r0,2)*w0*z1 - 
   4*k2*n0*r0*r1*w0*z1 + 2*j2*n1*r0*v0*w0*z1 + 2*j2*n0*r1*v0*w0*z1 - 
   2*i2*r1*s0*v0*w0*z1 - 2*i2*r0*s1*v0*w0*z1 + 2*j2*n0*r0*v1*w0*z1 - 
   2*i2*r0*s0*v1*w0*z1 + 2*i2*r0*r1*pow(w0,2)*z1 - 
   2*k2*n0*pow(r0,2)*w1*z1 + 2*j2*n0*r0*v0*w1*z1 - 2*i2*r0*s0*v0*w1*z1 + 
   2*i2*pow(r0,2)*w0*w1*z1
;

	k[18]  = 2*c1*pow(n1,2)*p2*r0*u0 + 4*c1*n0*n1*p2*r1*u0 + 2*c0*pow(n1,2)*p2*r1*u0 + 
   4*d1*k2*n1*r0*r1*u0 + 2*d1*k2*n0*pow(r1,2)*u0 + 
   2*d0*k2*n1*pow(r1,2)*u0 - 2*c1*k2*n1*r1*s0*u0 - 2*c1*k2*n1*r0*s1*u0 - 
   2*c1*k2*n0*r1*s1*u0 - 2*c0*k2*n1*r1*s1*u0 + 
   e1*pow(n1,2)*o2*pow(u0,2) + 2*f1*j2*n1*r1*pow(u0,2) + 
   g1*i2*pow(r1,2)*pow(u0,2) - 2*n1*n2*pow(r1,2)*pow(u0,2) - 
   2*pow(n1,2)*r1*r2*pow(u0,2) - 2*e1*j2*n1*s1*pow(u0,2) - 
   2*f1*i2*r1*s1*pow(u0,2) + 2*m2*n1*r1*s1*pow(u0,2) + 
   e1*i2*pow(s1,2)*pow(u0,2) + 4*c1*n0*n1*p2*r0*u1 + 
   2*c0*pow(n1,2)*p2*r0*u1 + 2*d1*k2*n1*pow(r0,2)*u1 + 
   2*c1*pow(n0,2)*p2*r1*u1 + 4*c0*n0*n1*p2*r1*u1 + 4*d1*k2*n0*r0*r1*u1 + 
   4*d0*k2*n1*r0*r1*u1 + 2*d0*k2*n0*pow(r1,2)*u1 - 2*c1*k2*n1*r0*s0*u1 - 
   2*c1*k2*n0*r1*s0*u1 - 2*c0*k2*n1*r1*s0*u1 - 2*c1*k2*n0*r0*s1*u1 - 
   2*c0*k2*n1*r0*s1*u1 - 2*c0*k2*n0*r1*s1*u1 + 4*e1*n0*n1*o2*u0*u1 + 
   2*e0*pow(n1,2)*o2*u0*u1 + 4*f1*j2*n1*r0*u0*u1 + 4*f1*j2*n0*r1*u0*u1 + 
   4*f0*j2*n1*r1*u0*u1 + 4*g1*i2*r0*r1*u0*u1 - 8*n1*n2*r0*r1*u0*u1 + 
   2*g0*i2*pow(r1,2)*u0*u1 - 4*n0*n2*pow(r1,2)*u0*u1 - 
   4*pow(n1,2)*r0*r2*u0*u1 - 8*n0*n1*r1*r2*u0*u1 - 4*e1*j2*n1*s0*u0*u1 - 
   4*f1*i2*r1*s0*u0*u1 + 4*m2*n1*r1*s0*u0*u1 - 4*e1*j2*n0*s1*u0*u1 - 
   4*e0*j2*n1*s1*u0*u1 - 4*f1*i2*r0*s1*u0*u1 + 4*m2*n1*r0*s1*u0*u1 - 
   4*f0*i2*r1*s1*u0*u1 + 4*m2*n0*r1*s1*u0*u1 + 4*e1*i2*s0*s1*u0*u1 + 
   2*e0*i2*pow(s1,2)*u0*u1 + e1*pow(n0,2)*o2*pow(u1,2) + 
   2*e0*n0*n1*o2*pow(u1,2) + 2*f1*j2*n0*r0*pow(u1,2) + 
   2*f0*j2*n1*r0*pow(u1,2) + g1*i2*pow(r0,2)*pow(u1,2) - 
   2*n1*n2*pow(r0,2)*pow(u1,2) + 2*f0*j2*n0*r1*pow(u1,2) + 
   2*g0*i2*r0*r1*pow(u1,2) - 4*n0*n2*r0*r1*pow(u1,2) - 
   4*n0*n1*r0*r2*pow(u1,2) - 2*pow(n0,2)*r1*r2*pow(u1,2) - 
   2*e1*j2*n0*s0*pow(u1,2) - 2*e0*j2*n1*s0*pow(u1,2) - 
   2*f1*i2*r0*s0*pow(u1,2) + 2*m2*n1*r0*s0*pow(u1,2) - 
   2*f0*i2*r1*s0*pow(u1,2) + 2*m2*n0*r1*s0*pow(u1,2) + 
   e1*i2*pow(s0,2)*pow(u1,2) - 2*e0*j2*n0*s1*pow(u1,2) - 
   2*f0*i2*r0*s1*pow(u1,2) + 2*m2*n0*r0*s1*pow(u1,2) + 
   2*e0*i2*s0*s1*pow(u1,2) - 4*pow(n1,2)*r0*r1*u0*u2 - 
   4*n0*n1*pow(r1,2)*u0*u2 - 2*pow(n1,2)*pow(r0,2)*u1*u2 - 
   8*n0*n1*r0*r1*u1*u2 - 2*pow(n0,2)*pow(r1,2)*u1*u2 - 
   2*c1*pow(n1,2)*o2*u0*v0 - 2*d1*j2*n1*r1*u0*v0 + 
   2*pow(n1,2)*q2*r1*u0*v0 + 4*c1*j2*n1*s1*u0*v0 + 2*d1*i2*r1*s1*u0*v0 - 
   2*l2*n1*r1*s1*u0*v0 - 2*c1*i2*pow(s1,2)*u0*v0 - 4*c1*n0*n1*o2*u1*v0 - 
   2*c0*pow(n1,2)*o2*u1*v0 - 2*d1*j2*n1*r0*u1*v0 + 
   2*pow(n1,2)*q2*r0*u1*v0 - 2*d1*j2*n0*r1*u1*v0 - 2*d0*j2*n1*r1*u1*v0 + 
   4*n0*n1*q2*r1*u1*v0 + 4*c1*j2*n1*s0*u1*v0 + 2*d1*i2*r1*s0*u1*v0 - 
   2*l2*n1*r1*s0*u1*v0 + 4*c1*j2*n0*s1*u1*v0 + 4*c0*j2*n1*s1*u1*v0 + 
   2*d1*i2*r0*s1*u1*v0 - 2*l2*n1*r0*s1*u1*v0 + 2*d0*i2*r1*s1*u1*v0 - 
   2*l2*n0*r1*s1*u1*v0 - 4*c1*i2*s0*s1*u1*v0 - 2*c0*i2*pow(s1,2)*u1*v0 - 
   4*c1*n0*n1*o2*u0*v1 - 2*c0*pow(n1,2)*o2*u0*v1 - 2*d1*j2*n1*r0*u0*v1 + 
   2*pow(n1,2)*q2*r0*u0*v1 - 2*d1*j2*n0*r1*u0*v1 - 2*d0*j2*n1*r1*u0*v1 + 
   4*n0*n1*q2*r1*u0*v1 + 4*c1*j2*n1*s0*u0*v1 + 2*d1*i2*r1*s0*u0*v1 - 
   2*l2*n1*r1*s0*u0*v1 + 4*c1*j2*n0*s1*u0*v1 + 4*c0*j2*n1*s1*u0*v1 + 
   2*d1*i2*r0*s1*u0*v1 - 2*l2*n1*r0*s1*u0*v1 + 2*d0*i2*r1*s1*u0*v1 - 
   2*l2*n0*r1*s1*u0*v1 - 4*c1*i2*s0*s1*u0*v1 - 2*c0*i2*pow(s1,2)*u0*v1 - 
   2*c1*pow(n0,2)*o2*u1*v1 - 4*c0*n0*n1*o2*u1*v1 - 2*d1*j2*n0*r0*u1*v1 - 
   2*d0*j2*n1*r0*u1*v1 + 4*n0*n1*q2*r0*u1*v1 - 2*d0*j2*n0*r1*u1*v1 + 
   2*pow(n0,2)*q2*r1*u1*v1 + 4*c1*j2*n0*s0*u1*v1 + 4*c0*j2*n1*s0*u1*v1 + 
   2*d1*i2*r0*s0*u1*v1 - 2*l2*n1*r0*s0*u1*v1 + 2*d0*i2*r1*s0*u1*v1 - 
   2*l2*n0*r1*s0*u1*v1 - 2*c1*i2*pow(s0,2)*u1*v1 + 4*c0*j2*n0*s1*u1*v1 + 
   2*d0*i2*r0*s1*u1*v1 - 2*l2*n0*r0*s1*u1*v1 - 4*c0*i2*s0*s1*u1*v1 - 
   2*c1*j2*n1*r1*u0*w0 - 2*d1*i2*pow(r1,2)*u0*w0 + 
   2*l2*n1*pow(r1,2)*u0*w0 + 2*c1*i2*r1*s1*u0*w0 - 2*c1*j2*n1*r0*u1*w0 - 
   2*c1*j2*n0*r1*u1*w0 - 2*c0*j2*n1*r1*u1*w0 - 4*d1*i2*r0*r1*u1*w0 + 
   4*l2*n1*r0*r1*u1*w0 - 2*d0*i2*pow(r1,2)*u1*w0 + 
   2*l2*n0*pow(r1,2)*u1*w0 + 2*c1*i2*r1*s0*u1*w0 + 2*c1*i2*r0*s1*u1*w0 + 
   2*c0*i2*r1*s1*u1*w0 - 2*c1*j2*n1*r0*u0*w1 - 2*c1*j2*n0*r1*u0*w1 - 
   2*c0*j2*n1*r1*u0*w1 - 4*d1*i2*r0*r1*u0*w1 + 4*l2*n1*r0*r1*u0*w1 - 
   2*d0*i2*pow(r1,2)*u0*w1 + 2*l2*n0*pow(r1,2)*u0*w1 + 
   2*c1*i2*r1*s0*u0*w1 + 2*c1*i2*r0*s1*u0*w1 + 2*c0*i2*r1*s1*u0*w1 - 
   2*c1*j2*n0*r0*u1*w1 - 2*c0*j2*n1*r0*u1*w1 - 2*d1*i2*pow(r0,2)*u1*w1 + 
   2*l2*n1*pow(r0,2)*u1*w1 - 2*c0*j2*n0*r1*u1*w1 - 4*d0*i2*r0*r1*u1*w1 + 
   4*l2*n0*r0*r1*u1*w1 + 2*c1*i2*r0*s0*u1*w1 + 2*c0*i2*r1*s0*u1*w1 + 
   2*c0*i2*r0*s1*u1*w1 + 2*pow(n1,2)*r0*r1*t2*z0 + 
   2*n0*n1*pow(r1,2)*t2*z0 - 2*pow(n1,2)*p2*r1*v0*z0 + 
   2*k2*n1*r1*s1*v0*z0 - 2*pow(n1,2)*p2*r0*v1*z0 - 4*n0*n1*p2*r1*v1*z0 + 
   2*k2*n1*r1*s0*v1*z0 + 2*k2*n1*r0*s1*v1*z0 + 2*k2*n0*r1*s1*v1*z0 + 
   2*pow(n1,2)*o2*v0*v1*z0 - 4*j2*n1*s1*v0*v1*z0 + 
   2*i2*pow(s1,2)*v0*v1*z0 + 2*n0*n1*o2*pow(v1,2)*z0 - 
   2*j2*n1*s0*pow(v1,2)*z0 - 2*j2*n0*s1*pow(v1,2)*z0 + 
   2*i2*s0*s1*pow(v1,2)*z0 - 2*k2*n1*pow(r1,2)*w0*z0 + 
   2*j2*n1*r1*v1*w0*z0 - 2*i2*r1*s1*v1*w0*z0 - 4*k2*n1*r0*r1*w1*z0 - 
   2*k2*n0*pow(r1,2)*w1*z0 + 2*j2*n1*r1*v0*w1*z0 - 2*i2*r1*s1*v0*w1*z0 + 
   2*j2*n1*r0*v1*w1*z0 + 2*j2*n0*r1*v1*w1*z0 - 2*i2*r1*s0*v1*w1*z0 - 
   2*i2*r0*s1*v1*w1*z0 + 2*i2*pow(r1,2)*w0*w1*z0 + 
   2*i2*r0*r1*pow(w1,2)*z0 + pow(n1,2)*pow(r0,2)*t2*z1 + 
   4*n0*n1*r0*r1*t2*z1 + pow(n0,2)*pow(r1,2)*t2*z1 - 
   2*pow(n1,2)*p2*r0*v0*z1 - 4*n0*n1*p2*r1*v0*z1 + 2*k2*n1*r1*s0*v0*z1 + 
   2*k2*n1*r0*s1*v0*z1 + 2*k2*n0*r1*s1*v0*z1 + 
   pow(n1,2)*o2*pow(v0,2)*z1 - 2*j2*n1*s1*pow(v0,2)*z1 + 
   i2*pow(s1,2)*pow(v0,2)*z1 - 4*n0*n1*p2*r0*v1*z1 - 
   2*pow(n0,2)*p2*r1*v1*z1 + 2*k2*n1*r0*s0*v1*z1 + 2*k2*n0*r1*s0*v1*z1 + 
   2*k2*n0*r0*s1*v1*z1 + 4*n0*n1*o2*v0*v1*z1 - 4*j2*n1*s0*v0*v1*z1 - 
   4*j2*n0*s1*v0*v1*z1 + 4*i2*s0*s1*v0*v1*z1 + 
   pow(n0,2)*o2*pow(v1,2)*z1 - 2*j2*n0*s0*pow(v1,2)*z1 + 
   i2*pow(s0,2)*pow(v1,2)*z1 - 4*k2*n1*r0*r1*w0*z1 - 
   2*k2*n0*pow(r1,2)*w0*z1 + 2*j2*n1*r1*v0*w0*z1 - 2*i2*r1*s1*v0*w0*z1 + 
   2*j2*n1*r0*v1*w0*z1 + 2*j2*n0*r1*v1*w0*z1 - 2*i2*r1*s0*v1*w0*z1 - 
   2*i2*r0*s1*v1*w0*z1 + i2*pow(r1,2)*pow(w0,2)*z1 - 
   2*k2*n1*pow(r0,2)*w1*z1 - 4*k2*n0*r0*r1*w1*z1 + 2*j2*n1*r0*v0*w1*z1 + 
   2*j2*n0*r1*v0*w1*z1 - 2*i2*r1*s0*v0*w1*z1 - 2*i2*r0*s1*v0*w1*z1 + 
   2*j2*n0*r0*v1*w1*z1 - 2*i2*r0*s0*v1*w1*z1 + 4*i2*r0*r1*w0*w1*z1 + 
   i2*pow(r0,2)*pow(w1,2)*z1
;

	k[19]  = 2*c1*pow(n1,2)*p2*r1*u0 + 2*d1*k2*n1*pow(r1,2)*u0 - 2*c1*k2*n1*r1*s1*u0 + 
   2*c1*pow(n1,2)*p2*r0*u1 + 4*c1*n0*n1*p2*r1*u1 + 
   2*c0*pow(n1,2)*p2*r1*u1 + 4*d1*k2*n1*r0*r1*u1 + 
   2*d1*k2*n0*pow(r1,2)*u1 + 2*d0*k2*n1*pow(r1,2)*u1 - 
   2*c1*k2*n1*r1*s0*u1 - 2*c1*k2*n1*r0*s1*u1 - 2*c1*k2*n0*r1*s1*u1 - 
   2*c0*k2*n1*r1*s1*u1 + 2*e1*pow(n1,2)*o2*u0*u1 + 4*f1*j2*n1*r1*u0*u1 + 
   2*g1*i2*pow(r1,2)*u0*u1 - 4*n1*n2*pow(r1,2)*u0*u1 - 
   4*pow(n1,2)*r1*r2*u0*u1 - 4*e1*j2*n1*s1*u0*u1 - 4*f1*i2*r1*s1*u0*u1 + 
   4*m2*n1*r1*s1*u0*u1 + 2*e1*i2*pow(s1,2)*u0*u1 + 
   2*e1*n0*n1*o2*pow(u1,2) + e0*pow(n1,2)*o2*pow(u1,2) + 
   2*f1*j2*n1*r0*pow(u1,2) + 2*f1*j2*n0*r1*pow(u1,2) + 
   2*f0*j2*n1*r1*pow(u1,2) + 2*g1*i2*r0*r1*pow(u1,2) - 
   4*n1*n2*r0*r1*pow(u1,2) + g0*i2*pow(r1,2)*pow(u1,2) - 
   2*n0*n2*pow(r1,2)*pow(u1,2) - 2*pow(n1,2)*r0*r2*pow(u1,2) - 
   4*n0*n1*r1*r2*pow(u1,2) - 2*e1*j2*n1*s0*pow(u1,2) - 
   2*f1*i2*r1*s0*pow(u1,2) + 2*m2*n1*r1*s0*pow(u1,2) - 
   2*e1*j2*n0*s1*pow(u1,2) - 2*e0*j2*n1*s1*pow(u1,2) - 
   2*f1*i2*r0*s1*pow(u1,2) + 2*m2*n1*r0*s1*pow(u1,2) - 
   2*f0*i2*r1*s1*pow(u1,2) + 2*m2*n0*r1*s1*pow(u1,2) + 
   2*e1*i2*s0*s1*pow(u1,2) + e0*i2*pow(s1,2)*pow(u1,2) - 
   2*pow(n1,2)*pow(r1,2)*u0*u2 - 4*pow(n1,2)*r0*r1*u1*u2 - 
   4*n0*n1*pow(r1,2)*u1*u2 - 2*c1*pow(n1,2)*o2*u1*v0 - 
   2*d1*j2*n1*r1*u1*v0 + 2*pow(n1,2)*q2*r1*u1*v0 + 4*c1*j2*n1*s1*u1*v0 + 
   2*d1*i2*r1*s1*u1*v0 - 2*l2*n1*r1*s1*u1*v0 - 2*c1*i2*pow(s1,2)*u1*v0 - 
   2*c1*pow(n1,2)*o2*u0*v1 - 2*d1*j2*n1*r1*u0*v1 + 
   2*pow(n1,2)*q2*r1*u0*v1 + 4*c1*j2*n1*s1*u0*v1 + 2*d1*i2*r1*s1*u0*v1 - 
   2*l2*n1*r1*s1*u0*v1 - 2*c1*i2*pow(s1,2)*u0*v1 - 4*c1*n0*n1*o2*u1*v1 - 
   2*c0*pow(n1,2)*o2*u1*v1 - 2*d1*j2*n1*r0*u1*v1 + 
   2*pow(n1,2)*q2*r0*u1*v1 - 2*d1*j2*n0*r1*u1*v1 - 2*d0*j2*n1*r1*u1*v1 + 
   4*n0*n1*q2*r1*u1*v1 + 4*c1*j2*n1*s0*u1*v1 + 2*d1*i2*r1*s0*u1*v1 - 
   2*l2*n1*r1*s0*u1*v1 + 4*c1*j2*n0*s1*u1*v1 + 4*c0*j2*n1*s1*u1*v1 + 
   2*d1*i2*r0*s1*u1*v1 - 2*l2*n1*r0*s1*u1*v1 + 2*d0*i2*r1*s1*u1*v1 - 
   2*l2*n0*r1*s1*u1*v1 - 4*c1*i2*s0*s1*u1*v1 - 2*c0*i2*pow(s1,2)*u1*v1 - 
   2*c1*j2*n1*r1*u1*w0 - 2*d1*i2*pow(r1,2)*u1*w0 + 
   2*l2*n1*pow(r1,2)*u1*w0 + 2*c1*i2*r1*s1*u1*w0 - 2*c1*j2*n1*r1*u0*w1 - 
   2*d1*i2*pow(r1,2)*u0*w1 + 2*l2*n1*pow(r1,2)*u0*w1 + 
   2*c1*i2*r1*s1*u0*w1 - 2*c1*j2*n1*r0*u1*w1 - 2*c1*j2*n0*r1*u1*w1 - 
   2*c0*j2*n1*r1*u1*w1 - 4*d1*i2*r0*r1*u1*w1 + 4*l2*n1*r0*r1*u1*w1 - 
   2*d0*i2*pow(r1,2)*u1*w1 + 2*l2*n0*pow(r1,2)*u1*w1 + 
   2*c1*i2*r1*s0*u1*w1 + 2*c1*i2*r0*s1*u1*w1 + 2*c0*i2*r1*s1*u1*w1 + 
   pow(n1,2)*pow(r1,2)*t2*z0 - 2*pow(n1,2)*p2*r1*v1*z0 + 
   2*k2*n1*r1*s1*v1*z0 + pow(n1,2)*o2*pow(v1,2)*z0 - 
   2*j2*n1*s1*pow(v1,2)*z0 + i2*pow(s1,2)*pow(v1,2)*z0 - 
   2*k2*n1*pow(r1,2)*w1*z0 + 2*j2*n1*r1*v1*w1*z0 - 2*i2*r1*s1*v1*w1*z0 + 
   i2*pow(r1,2)*pow(w1,2)*z0 + 2*pow(n1,2)*r0*r1*t2*z1 + 
   2*n0*n1*pow(r1,2)*t2*z1 - 2*pow(n1,2)*p2*r1*v0*z1 + 
   2*k2*n1*r1*s1*v0*z1 - 2*pow(n1,2)*p2*r0*v1*z1 - 4*n0*n1*p2*r1*v1*z1 + 
   2*k2*n1*r1*s0*v1*z1 + 2*k2*n1*r0*s1*v1*z1 + 2*k2*n0*r1*s1*v1*z1 + 
   2*pow(n1,2)*o2*v0*v1*z1 - 4*j2*n1*s1*v0*v1*z1 + 
   2*i2*pow(s1,2)*v0*v1*z1 + 2*n0*n1*o2*pow(v1,2)*z1 - 
   2*j2*n1*s0*pow(v1,2)*z1 - 2*j2*n0*s1*pow(v1,2)*z1 + 
   2*i2*s0*s1*pow(v1,2)*z1 - 2*k2*n1*pow(r1,2)*w0*z1 + 
   2*j2*n1*r1*v1*w0*z1 - 2*i2*r1*s1*v1*w0*z1 - 4*k2*n1*r0*r1*w1*z1 - 
   2*k2*n0*pow(r1,2)*w1*z1 + 2*j2*n1*r1*v0*w1*z1 - 2*i2*r1*s1*v0*w1*z1 + 
   2*j2*n1*r0*v1*w1*z1 + 2*j2*n0*r1*v1*w1*z1 - 2*i2*r1*s0*v1*w1*z1 - 
   2*i2*r0*s1*v1*w1*z1 + 2*i2*pow(r1,2)*w0*w1*z1 + 2*i2*r0*r1*pow(w1,2)*z1
;

	k[20]  = 2*c1*pow(n1,2)*p2*r1*u1 + 2*d1*k2*n1*pow(r1,2)*u1 - 2*c1*k2*n1*r1*s1*u1 + 
   e1*pow(n1,2)*o2*pow(u1,2) + 2*f1*j2*n1*r1*pow(u1,2) + 
   g1*i2*pow(r1,2)*pow(u1,2) - 2*n1*n2*pow(r1,2)*pow(u1,2) - 
   2*pow(n1,2)*r1*r2*pow(u1,2) - 2*e1*j2*n1*s1*pow(u1,2) - 
   2*f1*i2*r1*s1*pow(u1,2) + 2*m2*n1*r1*s1*pow(u1,2) + 
   e1*i2*pow(s1,2)*pow(u1,2) - 2*pow(n1,2)*pow(r1,2)*u1*u2 - 
   2*c1*pow(n1,2)*o2*u1*v1 - 2*d1*j2*n1*r1*u1*v1 + 
   2*pow(n1,2)*q2*r1*u1*v1 + 4*c1*j2*n1*s1*u1*v1 + 2*d1*i2*r1*s1*u1*v1 - 
   2*l2*n1*r1*s1*u1*v1 - 2*c1*i2*pow(s1,2)*u1*v1 - 2*c1*j2*n1*r1*u1*w1 - 
   2*d1*i2*pow(r1,2)*u1*w1 + 2*l2*n1*pow(r1,2)*u1*w1 + 
   2*c1*i2*r1*s1*u1*w1 + pow(n1,2)*pow(r1,2)*t2*z1 - 
   2*pow(n1,2)*p2*r1*v1*z1 + 2*k2*n1*r1*s1*v1*z1 + 
   pow(n1,2)*o2*pow(v1,2)*z1 - 2*j2*n1*s1*pow(v1,2)*z1 + 
   i2*pow(s1,2)*pow(v1,2)*z1 - 2*k2*n1*pow(r1,2)*w1*z1 + 
   2*j2*n1*r1*v1*w1*z1 - 2*i2*r1*s1*v1*w1*z1 + i2*pow(r1,2)*pow(w1,2)*z1
;

	k[21]  = -(pow(n0,2)*pow(r0,2)*pow(u0,2))
;

	k[22]  = -2*n0*n1*pow(r0,2)*pow(u0,2) - 2*pow(n0,2)*r0*r1*pow(u0,2) - 
   2*pow(n0,2)*pow(r0,2)*u0*u1
;

	k[23]  = -(pow(n1,2)*pow(r0,2)*pow(u0,2)) - 4*n0*n1*r0*r1*pow(u0,2) - 
   pow(n0,2)*pow(r1,2)*pow(u0,2) - 4*n0*n1*pow(r0,2)*u0*u1 - 
   4*pow(n0,2)*r0*r1*u0*u1 - pow(n0,2)*pow(r0,2)*pow(u1,2)
;

	k[24]  = -2*pow(n1,2)*r0*r1*pow(u0,2) - 2*n0*n1*pow(r1,2)*pow(u0,2) - 
   2*pow(n1,2)*pow(r0,2)*u0*u1 - 8*n0*n1*r0*r1*u0*u1 - 
   2*pow(n0,2)*pow(r1,2)*u0*u1 - 2*n0*n1*pow(r0,2)*pow(u1,2) - 
   2*pow(n0,2)*r0*r1*pow(u1,2)
;

	k[25]  = -(pow(n1,2)*pow(r1,2)*pow(u0,2)) - 4*pow(n1,2)*r0*r1*u0*u1 - 
   4*n0*n1*pow(r1,2)*u0*u1 - pow(n1,2)*pow(r0,2)*pow(u1,2) - 
   4*n0*n1*r0*r1*pow(u1,2) - pow(n0,2)*pow(r1,2)*pow(u1,2)
;
	
	k[26]  = -2*pow(n1,2)*pow(r1,2)*u0*u1 - 2*pow(n1,2)*r0*r1*pow(u1,2) - 
   2*n0*n1*pow(r1,2)*pow(u1,2)
;
	
	k[27]  = -(pow(n1,2)*pow(r1,2)*pow(u1,2))
;

	return 	rg_ImplicitEquation(6, k);
}

rg_REAL rg_IntsecImplicit6::inversion(const rg_BzCurve2D &curve, const rg_Point2D &point)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL x3 = curve.getCtrlPt(3).getX();
	rg_REAL x4 = curve.getCtrlPt(4).getX();
    rg_REAL x5 = curve.getCtrlPt(5).getX();
	rg_REAL x6 = curve.getCtrlPt(6).getX();

	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL y3 = curve.getCtrlPt(3).getY();
	rg_REAL y4 = curve.getCtrlPt(4).getY();
    rg_REAL y5 = curve.getCtrlPt(5).getY();
	rg_REAL y6 = curve.getCtrlPt(6).getY();

//  If
//		    
//	x(t) = a6 t^6 + a5 t^5 + a4 t^4 + a3 t^3 + a2 t^2 + a1 t + a0
//
//	y(t) = b6 t^6 + b5 t^5 + b4 t^4 + b3 t^3 + b2 t^2 + b1 t + b0
//  
//  then the resultant of p(x, t) and q(y, t) is like the following.

	rg_REAL a6 =    x0 - 6*x1 + 15*x2 - 20*x3 + 15*x4 - 6*x5 + x6;
    rg_REAL a5 = -6*x0 + 30*x1 - 60*x2 + 60*x3 - 30*x4 + 6*x5;
	rg_REAL a4 = 15*x0 - 60*x1 + 90*x2 - 60*x3 + 15*x4;
	rg_REAL a3 =-20*x0 + 60*x1 - 60*x2 + 20*x3;
	rg_REAL a2 = 15*x0 - 30*x1 + 15*x2;
	rg_REAL a1 = -6*x0 + 6*x1;
	rg_REAL a0 =    x0;

	rg_REAL b6 = y0 - 6*y1 + 15*y2 - 20*y3 + 15*y4 - 6*y5 + y6;
    rg_REAL b5 = -6*y0 + 30*y1 - 60*y2 + 60*y3 - 30*y4 + 6*y5;
	rg_REAL b4 = 15*y0 - 60*y1 + 90*y2 - 60*y3 + 15*y4;
	rg_REAL b3 = -20*y0 + 60*y1 - 60*y2 + 20*y3;
	rg_REAL b2 = 15*y0 - 30*y1 + 15*y2;
	rg_REAL b1 = -6*y0 + 6*y1;
	rg_REAL b0 =    y0;

  //			|														                                                       |
  //			|	i0 x + i1 y + i2  j0 x + j1 y + j2	k0 x + k1 y + k2  l0 x + l1 y + l2  m0 x + m1 y + m2  n0 x + n1 y + n2 |
  //			|														                                                       |
  // R(P, Q) =	|	j0 x + j1 y + j2  o0 x + o1 y + o2	p0 x + p1 y + p2  q0 x + q1 y + q2  r0 x + r1 y + r2  s0 x + s1 y + s2 |
  //			|														                                                       |
  //			|	k0 x + k1 y + k2  p0 x + p1 y + p2	t0 x + t1 y + t2  u0 x + u1 y + u2  v0 x + v1 y + v2  w0 x + w1 y + w2 |
  //			|														                                                       |
  //            |   l0 x + l1 y + l2  q0 x + q1 y + q2  u0 x + u1 y + u2  z0 x + z1 y + z2  c0 x + c1 y + c2  d0 x + d1 y + d2 |
  //            |                                                                                                              |
  //            |   m0 x + m1 y + m2  r0 x + r1 y + r2  v0 x + v1 y + v2  c0 x + c1 y + c2  e0 x + e1 y + e2  f0 x + f1 y + f2 |
  //            |                                                                                                              |
  //            |   n0 x + n1 y + n2  s0 x + s1 y + s2  w0 x + w1 y + w2  d0 x + d1 y + d2  f0 x + f1 y + f2  g0 x + g1 y + g2 |
  //            |                                                                                                              |

	rg_REAL i0 =   0.0;
	rg_REAL i1 =   0.0;
	rg_REAL i2 =   a6*b5 - a5*b6;

	rg_REAL j0 =   0.0;
	rg_REAL j1 =   0.0;
	rg_REAL j2 =   a6*b4 - a4*b6;

	rg_REAL k0 =   0.0;
	rg_REAL k1 =   0.0;
	rg_REAL k2 =   a6*b3 - a3*b6;

	rg_REAL l0 =   0.0;
	rg_REAL l1 =   0.0;
	rg_REAL l2 =   a6*b2 - a2*b6;

	rg_REAL m0 =   0.0;
	rg_REAL m1 =   0.0;
	rg_REAL m2 =   a6*b1 - a1*b6;

	rg_REAL n0 =    b6;
	rg_REAL n1 =   -a6;
	rg_REAL n2 =   a6*b0 - a0*b6;

	rg_REAL o0 =   0.0;
	rg_REAL o1 =   0.0;
	rg_REAL o2 =   a5*b4 - a4*b5 + a6*b3 - a3*b6;

    rg_REAL p0 =   0.0;
	rg_REAL p1 =   0.0;
	rg_REAL p2 =   a5*b3 - a3*b5 + a6*b2 - a2*b6;

  	rg_REAL q0 =   0.0;
	rg_REAL q1 =   0.0;
	rg_REAL q2 =   a5*b2 - a2*b5 + a6*b1 - a1*b6;

  	rg_REAL r0 =   b6;
	rg_REAL r1 = - a6;
	rg_REAL r2 =   a5*b1 - a1*b5 + a6*b0 - a0*b6;

  	rg_REAL s0 =   b5;
	rg_REAL s1 = - a5;
	rg_REAL s2 =   a5*b0 - a0*b5;

   	rg_REAL t0 =   0.0;
	rg_REAL t1 =   0.0;
	rg_REAL t2 =   a4*b3 - a3*b4 + a5*b2 - a2*b5 + a6*b1 - a1*b6;

   	rg_REAL u0 =   b6;
	rg_REAL u1 = - a6;
	rg_REAL u2 =   a4*b2 - a2*b4 + a5*b1 - a1*b5 + a6*b0 - a0*b6;

	rg_REAL v0 =   b5;
	rg_REAL v1 = - a5;
	rg_REAL v2 =   a4*b1 - a1*b4 + a5*b0 - a0*b5;

   	rg_REAL w0 =   b4;
	rg_REAL w1 = - a4;
	rg_REAL w2 =   a4*b0 - a0*b4;

   	rg_REAL z0 =   b5;
	rg_REAL z1 = - a5;
	rg_REAL z2 =   a3*b2 - a2*b3 + a4*b1 - a1*b4 + a5*b0 - a0*b5;

   	rg_REAL c0 =   b4;
	rg_REAL c1 = - a4;
	rg_REAL c2 =   a3*b1 - a1*b3 + a4*b0 - a0*b4;

   	rg_REAL d0 =   b3;
	rg_REAL d1 = - a3;
	rg_REAL d2 =   a3*b0 - a0*b3;

   	rg_REAL e0 =   b3;
	rg_REAL e1 = - a3;
	rg_REAL e2 =   a2*b1 - a1*b2 + a3*b0 - a0*b3;

   	rg_REAL f0 =   b2;
	rg_REAL f1 = - a2;
	rg_REAL f2 =   a2*b0 - a0*b2;

   	rg_REAL g0 =   b1;
	rg_REAL g1 = - a1;
	rg_REAL g2 =   a1*b0 - a0*b1;

// 고쳐야할 부분
	rg_REAL x  =  point.getX();
	rg_REAL y  =  point.getY();

	rg_REAL t_4 = d2*e2*k2*l2*p2 - c2*f2*k2*l2*p2 - c2*d2*k2*m2*p2 + pow(c2,2)*k2*n2*p2 - 
   d2*e2*pow(k2,2)*q2 + c2*f2*pow(k2,2)*q2 + c2*d2*pow(k2,2)*r2 - 
   pow(c2,2)*pow(k2,2)*s2 - d2*e2*j2*l2*t2 + c2*f2*j2*l2*t2 + 
   c2*d2*j2*m2*t2 - pow(c2,2)*j2*n2*t2 + d2*e2*i2*q2*t2 - c2*f2*i2*q2*t2 + 
   f2*l2*m2*q2*t2 - d2*pow(m2,2)*q2*t2 - e2*l2*n2*q2*t2 + c2*m2*n2*q2*t2 - 
   c2*d2*i2*r2*t2 - f2*pow(l2,2)*r2*t2 + d2*l2*m2*r2*t2 + c2*l2*n2*r2*t2 + 
   pow(c2,2)*i2*s2*t2 + e2*pow(l2,2)*s2*t2 - 2*c2*l2*m2*s2*t2 + 
   d2*e2*j2*k2*u2 - c2*f2*j2*k2*u2 - d2*e2*i2*p2*u2 + c2*f2*i2*p2*u2 - 
   f2*l2*m2*p2*u2 + d2*pow(m2,2)*p2*u2 + e2*l2*n2*p2*u2 - c2*m2*n2*p2*u2 - 
   f2*k2*m2*q2*u2 + e2*k2*n2*q2*u2 + 2*f2*k2*l2*r2*u2 - d2*k2*m2*r2*u2 - 
   c2*k2*n2*r2*u2 - 2*e2*k2*l2*s2*u2 + 2*c2*k2*m2*s2*u2 + 
   f2*j2*m2*pow(u2,2) - e2*j2*n2*pow(u2,2) - f2*i2*r2*pow(u2,2) + 
   m2*n2*r2*pow(u2,2) + e2*i2*s2*pow(u2,2) - pow(m2,2)*s2*pow(u2,2) - 
   c2*d2*j2*k2*v2 + c2*d2*i2*p2*v2 + f2*pow(l2,2)*p2*v2 - d2*l2*m2*p2*v2 - 
   c2*l2*n2*p2*v2 - f2*k2*l2*q2*v2 + 2*d2*k2*m2*q2*v2 - c2*k2*n2*q2*v2 - 
   d2*k2*l2*r2*v2 + 2*c2*k2*l2*s2*v2 - f2*j2*l2*u2*v2 - d2*j2*m2*u2*v2 + 
   2*c2*j2*n2*u2*v2 + f2*i2*q2*u2*v2 - m2*n2*q2*u2*v2 + d2*i2*r2*u2*v2 - 
   l2*n2*r2*u2*v2 - 2*c2*i2*s2*u2*v2 + 2*l2*m2*s2*u2*v2 + 
   d2*j2*l2*pow(v2,2) - d2*i2*q2*pow(v2,2) + l2*n2*q2*pow(v2,2) - 
   pow(l2,2)*s2*pow(v2,2) + pow(c2,2)*j2*k2*w2 - pow(c2,2)*i2*p2*w2 - 
   e2*pow(l2,2)*p2*w2 + 2*c2*l2*m2*p2*w2 + e2*k2*l2*q2*w2 - 
   c2*k2*m2*q2*w2 - c2*k2*l2*r2*w2 + e2*j2*l2*u2*w2 - c2*j2*m2*u2*w2 - 
   e2*i2*q2*u2*w2 + pow(m2,2)*q2*u2*w2 + c2*i2*r2*u2*w2 - l2*m2*r2*u2*w2 - 
   c2*j2*l2*v2*w2 + c2*i2*q2*v2*w2 - l2*m2*q2*v2*w2 + pow(l2,2)*r2*v2*w2 + 
   d2*e0*k2*l2*p2*x + d0*e2*k2*l2*p2*x - c2*f0*k2*l2*p2*x - 
   c0*f2*k2*l2*p2*x - c2*d0*k2*m2*p2*x - c0*d2*k2*m2*p2*x + 
   pow(c2,2)*k2*n0*p2*x + 2*c0*c2*k2*n2*p2*x - d2*e0*pow(k2,2)*q2*x - 
   d0*e2*pow(k2,2)*q2*x + c2*f0*pow(k2,2)*q2*x + c0*f2*pow(k2,2)*q2*x + 
   c2*d2*pow(k2,2)*r0*x + c2*d0*pow(k2,2)*r2*x + c0*d2*pow(k2,2)*r2*x - 
   pow(c2,2)*pow(k2,2)*s0*x - 2*c0*c2*pow(k2,2)*s2*x - 
   d2*e0*j2*l2*t2*x - d0*e2*j2*l2*t2*x + c2*f0*j2*l2*t2*x + 
   c0*f2*j2*l2*t2*x + c2*d0*j2*m2*t2*x + c0*d2*j2*m2*t2*x - 
   pow(c2,2)*j2*n0*t2*x - 2*c0*c2*j2*n2*t2*x + d2*e0*i2*q2*t2*x + 
   d0*e2*i2*q2*t2*x - c2*f0*i2*q2*t2*x - c0*f2*i2*q2*t2*x + 
   f0*l2*m2*q2*t2*x - d0*pow(m2,2)*q2*t2*x - e2*l2*n0*q2*t2*x + 
   c2*m2*n0*q2*t2*x - e0*l2*n2*q2*t2*x + c0*m2*n2*q2*t2*x - 
   c2*d2*i2*r0*t2*x - f2*pow(l2,2)*r0*t2*x + d2*l2*m2*r0*t2*x + 
   c2*l2*n2*r0*t2*x - c2*d0*i2*r2*t2*x - c0*d2*i2*r2*t2*x - 
   f0*pow(l2,2)*r2*t2*x + d0*l2*m2*r2*t2*x + c2*l2*n0*r2*t2*x + 
   c0*l2*n2*r2*t2*x + pow(c2,2)*i2*s0*t2*x + e2*pow(l2,2)*s0*t2*x - 
   2*c2*l2*m2*s0*t2*x + 2*c0*c2*i2*s2*t2*x + e0*pow(l2,2)*s2*t2*x - 
   2*c0*l2*m2*s2*t2*x + d2*e2*j2*k2*u0*x - c2*f2*j2*k2*u0*x - 
   d2*e2*i2*p2*u0*x + c2*f2*i2*p2*u0*x - f2*l2*m2*p2*u0*x + 
   d2*pow(m2,2)*p2*u0*x + e2*l2*n2*p2*u0*x - c2*m2*n2*p2*u0*x - 
   f2*k2*m2*q2*u0*x + e2*k2*n2*q2*u0*x + 2*f2*k2*l2*r2*u0*x - 
   d2*k2*m2*r2*u0*x - c2*k2*n2*r2*u0*x - 2*e2*k2*l2*s2*u0*x + 
   2*c2*k2*m2*s2*u0*x + d2*e0*j2*k2*u2*x + d0*e2*j2*k2*u2*x - 
   c2*f0*j2*k2*u2*x - c0*f2*j2*k2*u2*x - d2*e0*i2*p2*u2*x - 
   d0*e2*i2*p2*u2*x + c2*f0*i2*p2*u2*x + c0*f2*i2*p2*u2*x - 
   f0*l2*m2*p2*u2*x + d0*pow(m2,2)*p2*u2*x + e2*l2*n0*p2*u2*x - 
   c2*m2*n0*p2*u2*x + e0*l2*n2*p2*u2*x - c0*m2*n2*p2*u2*x - 
   f0*k2*m2*q2*u2*x + e2*k2*n0*q2*u2*x + e0*k2*n2*q2*u2*x + 
   2*f2*k2*l2*r0*u2*x - d2*k2*m2*r0*u2*x - c2*k2*n2*r0*u2*x + 
   2*f0*k2*l2*r2*u2*x - d0*k2*m2*r2*u2*x - c2*k2*n0*r2*u2*x - 
   c0*k2*n2*r2*u2*x - 2*e2*k2*l2*s0*u2*x + 2*c2*k2*m2*s0*u2*x - 
   2*e0*k2*l2*s2*u2*x + 2*c0*k2*m2*s2*u2*x + 2*f2*j2*m2*u0*u2*x - 
   2*e2*j2*n2*u0*u2*x - 2*f2*i2*r2*u0*u2*x + 2*m2*n2*r2*u0*u2*x + 
   2*e2*i2*s2*u0*u2*x - 2*pow(m2,2)*s2*u0*u2*x + f0*j2*m2*pow(u2,2)*x - 
   e2*j2*n0*pow(u2,2)*x - e0*j2*n2*pow(u2,2)*x - f2*i2*r0*pow(u2,2)*x + 
   m2*n2*r0*pow(u2,2)*x - f0*i2*r2*pow(u2,2)*x + m2*n0*r2*pow(u2,2)*x + 
   e2*i2*s0*pow(u2,2)*x - pow(m2,2)*s0*pow(u2,2)*x + 
   e0*i2*s2*pow(u2,2)*x - c2*d2*j2*k2*v0*x + c2*d2*i2*p2*v0*x + 
   f2*pow(l2,2)*p2*v0*x - d2*l2*m2*p2*v0*x - c2*l2*n2*p2*v0*x - 
   f2*k2*l2*q2*v0*x + 2*d2*k2*m2*q2*v0*x - c2*k2*n2*q2*v0*x - 
   d2*k2*l2*r2*v0*x + 2*c2*k2*l2*s2*v0*x - f2*j2*l2*u2*v0*x - 
   d2*j2*m2*u2*v0*x + 2*c2*j2*n2*u2*v0*x + f2*i2*q2*u2*v0*x - 
   m2*n2*q2*u2*v0*x + d2*i2*r2*u2*v0*x - l2*n2*r2*u2*v0*x - 
   2*c2*i2*s2*u2*v0*x + 2*l2*m2*s2*u2*v0*x - c2*d0*j2*k2*v2*x - 
   c0*d2*j2*k2*v2*x + c2*d0*i2*p2*v2*x + c0*d2*i2*p2*v2*x + 
   f0*pow(l2,2)*p2*v2*x - d0*l2*m2*p2*v2*x - c2*l2*n0*p2*v2*x - 
   c0*l2*n2*p2*v2*x - f0*k2*l2*q2*v2*x + 2*d0*k2*m2*q2*v2*x - 
   c2*k2*n0*q2*v2*x - c0*k2*n2*q2*v2*x - d2*k2*l2*r0*v2*x - 
   d0*k2*l2*r2*v2*x + 2*c2*k2*l2*s0*v2*x + 2*c0*k2*l2*s2*v2*x - 
   f2*j2*l2*u0*v2*x - d2*j2*m2*u0*v2*x + 2*c2*j2*n2*u0*v2*x + 
   f2*i2*q2*u0*v2*x - m2*n2*q2*u0*v2*x + d2*i2*r2*u0*v2*x - 
   l2*n2*r2*u0*v2*x - 2*c2*i2*s2*u0*v2*x + 2*l2*m2*s2*u0*v2*x - 
   f0*j2*l2*u2*v2*x - d0*j2*m2*u2*v2*x + 2*c2*j2*n0*u2*v2*x + 
   2*c0*j2*n2*u2*v2*x + f0*i2*q2*u2*v2*x - m2*n0*q2*u2*v2*x + 
   d2*i2*r0*u2*v2*x - l2*n2*r0*u2*v2*x + d0*i2*r2*u2*v2*x - 
   l2*n0*r2*u2*v2*x - 2*c2*i2*s0*u2*v2*x + 2*l2*m2*s0*u2*v2*x - 
   2*c0*i2*s2*u2*v2*x + 2*d2*j2*l2*v0*v2*x - 2*d2*i2*q2*v0*v2*x + 
   2*l2*n2*q2*v0*v2*x - 2*pow(l2,2)*s2*v0*v2*x + d0*j2*l2*pow(v2,2)*x - 
   d0*i2*q2*pow(v2,2)*x + l2*n0*q2*pow(v2,2)*x - 
   pow(l2,2)*s0*pow(v2,2)*x + pow(c2,2)*j2*k2*w0*x - 
   pow(c2,2)*i2*p2*w0*x - e2*pow(l2,2)*p2*w0*x + 2*c2*l2*m2*p2*w0*x + 
   e2*k2*l2*q2*w0*x - c2*k2*m2*q2*w0*x - c2*k2*l2*r2*w0*x + 
   e2*j2*l2*u2*w0*x - c2*j2*m2*u2*w0*x - e2*i2*q2*u2*w0*x + 
   pow(m2,2)*q2*u2*w0*x + c2*i2*r2*u2*w0*x - l2*m2*r2*u2*w0*x - 
   c2*j2*l2*v2*w0*x + c2*i2*q2*v2*w0*x - l2*m2*q2*v2*w0*x + 
   pow(l2,2)*r2*v2*w0*x + 2*c0*c2*j2*k2*w2*x - 2*c0*c2*i2*p2*w2*x - 
   e0*pow(l2,2)*p2*w2*x + 2*c0*l2*m2*p2*w2*x + e0*k2*l2*q2*w2*x - 
   c0*k2*m2*q2*w2*x - c2*k2*l2*r0*w2*x - c0*k2*l2*r2*w2*x + 
   e2*j2*l2*u0*w2*x - c2*j2*m2*u0*w2*x - e2*i2*q2*u0*w2*x + 
   pow(m2,2)*q2*u0*w2*x + c2*i2*r2*u0*w2*x - l2*m2*r2*u0*w2*x + 
   e0*j2*l2*u2*w2*x - c0*j2*m2*u2*w2*x - e0*i2*q2*u2*w2*x + 
   c2*i2*r0*u2*w2*x - l2*m2*r0*u2*w2*x + c0*i2*r2*u2*w2*x - 
   c2*j2*l2*v0*w2*x + c2*i2*q2*v0*w2*x - l2*m2*q2*v0*w2*x + 
   pow(l2,2)*r2*v0*w2*x - c0*j2*l2*v2*w2*x + c0*i2*q2*v2*w2*x + 
   pow(l2,2)*r0*v2*w2*x + d0*e0*k2*l2*p2*pow(x,2) - 
   c0*f0*k2*l2*p2*pow(x,2) - c0*d0*k2*m2*p2*pow(x,2) + 
   2*c0*c2*k2*n0*p2*pow(x,2) + pow(c0,2)*k2*n2*p2*pow(x,2) - 
   d0*e0*pow(k2,2)*q2*pow(x,2) + c0*f0*pow(k2,2)*q2*pow(x,2) + 
   c2*d0*pow(k2,2)*r0*pow(x,2) + c0*d2*pow(k2,2)*r0*pow(x,2) + 
   c0*d0*pow(k2,2)*r2*pow(x,2) - 2*c0*c2*pow(k2,2)*s0*pow(x,2) - 
   pow(c0,2)*pow(k2,2)*s2*pow(x,2) - d0*e0*j2*l2*t2*pow(x,2) + 
   c0*f0*j2*l2*t2*pow(x,2) + c0*d0*j2*m2*t2*pow(x,2) - 
   2*c0*c2*j2*n0*t2*pow(x,2) - pow(c0,2)*j2*n2*t2*pow(x,2) + 
   d0*e0*i2*q2*t2*pow(x,2) - c0*f0*i2*q2*t2*pow(x,2) - 
   e0*l2*n0*q2*t2*pow(x,2) + c0*m2*n0*q2*t2*pow(x,2) - 
   c2*d0*i2*r0*t2*pow(x,2) - c0*d2*i2*r0*t2*pow(x,2) - 
   f0*pow(l2,2)*r0*t2*pow(x,2) + d0*l2*m2*r0*t2*pow(x,2) + 
   c2*l2*n0*r0*t2*pow(x,2) + c0*l2*n2*r0*t2*pow(x,2) - 
   c0*d0*i2*r2*t2*pow(x,2) + c0*l2*n0*r2*t2*pow(x,2) + 
   2*c0*c2*i2*s0*t2*pow(x,2) + e0*pow(l2,2)*s0*t2*pow(x,2) - 
   2*c0*l2*m2*s0*t2*pow(x,2) + pow(c0,2)*i2*s2*t2*pow(x,2) + 
   d2*e0*j2*k2*u0*pow(x,2) + d0*e2*j2*k2*u0*pow(x,2) - 
   c2*f0*j2*k2*u0*pow(x,2) - c0*f2*j2*k2*u0*pow(x,2) - 
   d2*e0*i2*p2*u0*pow(x,2) - d0*e2*i2*p2*u0*pow(x,2) + 
   c2*f0*i2*p2*u0*pow(x,2) + c0*f2*i2*p2*u0*pow(x,2) - 
   f0*l2*m2*p2*u0*pow(x,2) + d0*pow(m2,2)*p2*u0*pow(x,2) + 
   e2*l2*n0*p2*u0*pow(x,2) - c2*m2*n0*p2*u0*pow(x,2) + 
   e0*l2*n2*p2*u0*pow(x,2) - c0*m2*n2*p2*u0*pow(x,2) - 
   f0*k2*m2*q2*u0*pow(x,2) + e2*k2*n0*q2*u0*pow(x,2) + 
   e0*k2*n2*q2*u0*pow(x,2) + 2*f2*k2*l2*r0*u0*pow(x,2) - 
   d2*k2*m2*r0*u0*pow(x,2) - c2*k2*n2*r0*u0*pow(x,2) + 
   2*f0*k2*l2*r2*u0*pow(x,2) - d0*k2*m2*r2*u0*pow(x,2) - 
   c2*k2*n0*r2*u0*pow(x,2) - c0*k2*n2*r2*u0*pow(x,2) - 
   2*e2*k2*l2*s0*u0*pow(x,2) + 2*c2*k2*m2*s0*u0*pow(x,2) - 
   2*e0*k2*l2*s2*u0*pow(x,2) + 2*c0*k2*m2*s2*u0*pow(x,2) + 
   f2*j2*m2*pow(u0,2)*pow(x,2) - e2*j2*n2*pow(u0,2)*pow(x,2) - 
   f2*i2*r2*pow(u0,2)*pow(x,2) + m2*n2*r2*pow(u0,2)*pow(x,2) + 
   e2*i2*s2*pow(u0,2)*pow(x,2) - pow(m2,2)*s2*pow(u0,2)*pow(x,2) + 
   d0*e0*j2*k2*u2*pow(x,2) - c0*f0*j2*k2*u2*pow(x,2) - 
   d0*e0*i2*p2*u2*pow(x,2) + c0*f0*i2*p2*u2*pow(x,2) + 
   e0*l2*n0*p2*u2*pow(x,2) - c0*m2*n0*p2*u2*pow(x,2) + 
   e0*k2*n0*q2*u2*pow(x,2) + 2*f0*k2*l2*r0*u2*pow(x,2) - 
   d0*k2*m2*r0*u2*pow(x,2) - c2*k2*n0*r0*u2*pow(x,2) - 
   c0*k2*n2*r0*u2*pow(x,2) - c0*k2*n0*r2*u2*pow(x,2) - 
   2*e0*k2*l2*s0*u2*pow(x,2) + 2*c0*k2*m2*s0*u2*pow(x,2) + 
   2*f0*j2*m2*u0*u2*pow(x,2) - 2*e2*j2*n0*u0*u2*pow(x,2) - 
   2*e0*j2*n2*u0*u2*pow(x,2) - 2*f2*i2*r0*u0*u2*pow(x,2) + 
   2*m2*n2*r0*u0*u2*pow(x,2) - 2*f0*i2*r2*u0*u2*pow(x,2) + 
   2*m2*n0*r2*u0*u2*pow(x,2) + 2*e2*i2*s0*u0*u2*pow(x,2) - 
   2*pow(m2,2)*s0*u0*u2*pow(x,2) + 2*e0*i2*s2*u0*u2*pow(x,2) - 
   e0*j2*n0*pow(u2,2)*pow(x,2) - f0*i2*r0*pow(u2,2)*pow(x,2) + 
   m2*n0*r0*pow(u2,2)*pow(x,2) + e0*i2*s0*pow(u2,2)*pow(x,2) - 
   c2*d0*j2*k2*v0*pow(x,2) - c0*d2*j2*k2*v0*pow(x,2) + 
   c2*d0*i2*p2*v0*pow(x,2) + c0*d2*i2*p2*v0*pow(x,2) + 
   f0*pow(l2,2)*p2*v0*pow(x,2) - d0*l2*m2*p2*v0*pow(x,2) - 
   c2*l2*n0*p2*v0*pow(x,2) - c0*l2*n2*p2*v0*pow(x,2) - 
   f0*k2*l2*q2*v0*pow(x,2) + 2*d0*k2*m2*q2*v0*pow(x,2) - 
   c2*k2*n0*q2*v0*pow(x,2) - c0*k2*n2*q2*v0*pow(x,2) - 
   d2*k2*l2*r0*v0*pow(x,2) - d0*k2*l2*r2*v0*pow(x,2) + 
   2*c2*k2*l2*s0*v0*pow(x,2) + 2*c0*k2*l2*s2*v0*pow(x,2) - 
   f2*j2*l2*u0*v0*pow(x,2) - d2*j2*m2*u0*v0*pow(x,2) + 
   2*c2*j2*n2*u0*v0*pow(x,2) + f2*i2*q2*u0*v0*pow(x,2) - 
   m2*n2*q2*u0*v0*pow(x,2) + d2*i2*r2*u0*v0*pow(x,2) - 
   l2*n2*r2*u0*v0*pow(x,2) - 2*c2*i2*s2*u0*v0*pow(x,2) + 
   2*l2*m2*s2*u0*v0*pow(x,2) - f0*j2*l2*u2*v0*pow(x,2) - 
   d0*j2*m2*u2*v0*pow(x,2) + 2*c2*j2*n0*u2*v0*pow(x,2) + 
   2*c0*j2*n2*u2*v0*pow(x,2) + f0*i2*q2*u2*v0*pow(x,2) - 
   m2*n0*q2*u2*v0*pow(x,2) + d2*i2*r0*u2*v0*pow(x,2) - 
   l2*n2*r0*u2*v0*pow(x,2) + d0*i2*r2*u2*v0*pow(x,2) - 
   l2*n0*r2*u2*v0*pow(x,2) - 2*c2*i2*s0*u2*v0*pow(x,2) + 
   2*l2*m2*s0*u2*v0*pow(x,2) - 2*c0*i2*s2*u2*v0*pow(x,2) + 
   d2*j2*l2*pow(v0,2)*pow(x,2) - d2*i2*q2*pow(v0,2)*pow(x,2) + 
   l2*n2*q2*pow(v0,2)*pow(x,2) - pow(l2,2)*s2*pow(v0,2)*pow(x,2) - 
   c0*d0*j2*k2*v2*pow(x,2) + c0*d0*i2*p2*v2*pow(x,2) - 
   c0*l2*n0*p2*v2*pow(x,2) - c0*k2*n0*q2*v2*pow(x,2) - 
   d0*k2*l2*r0*v2*pow(x,2) + 2*c0*k2*l2*s0*v2*pow(x,2) - 
   f0*j2*l2*u0*v2*pow(x,2) - d0*j2*m2*u0*v2*pow(x,2) + 
   2*c2*j2*n0*u0*v2*pow(x,2) + 2*c0*j2*n2*u0*v2*pow(x,2) + 
   f0*i2*q2*u0*v2*pow(x,2) - m2*n0*q2*u0*v2*pow(x,2) + 
   d2*i2*r0*u0*v2*pow(x,2) - l2*n2*r0*u0*v2*pow(x,2) + 
   d0*i2*r2*u0*v2*pow(x,2) - l2*n0*r2*u0*v2*pow(x,2) - 
   2*c2*i2*s0*u0*v2*pow(x,2) + 2*l2*m2*s0*u0*v2*pow(x,2) - 
   2*c0*i2*s2*u0*v2*pow(x,2) + 2*c0*j2*n0*u2*v2*pow(x,2) + 
   d0*i2*r0*u2*v2*pow(x,2) - l2*n0*r0*u2*v2*pow(x,2) - 
   2*c0*i2*s0*u2*v2*pow(x,2) + 2*d0*j2*l2*v0*v2*pow(x,2) - 
   2*d0*i2*q2*v0*v2*pow(x,2) + 2*l2*n0*q2*v0*v2*pow(x,2) - 
   2*pow(l2,2)*s0*v0*v2*pow(x,2) + 2*c0*c2*j2*k2*w0*pow(x,2) - 
   2*c0*c2*i2*p2*w0*pow(x,2) - e0*pow(l2,2)*p2*w0*pow(x,2) + 
   2*c0*l2*m2*p2*w0*pow(x,2) + e0*k2*l2*q2*w0*pow(x,2) - 
   c0*k2*m2*q2*w0*pow(x,2) - c2*k2*l2*r0*w0*pow(x,2) - 
   c0*k2*l2*r2*w0*pow(x,2) + e2*j2*l2*u0*w0*pow(x,2) - 
   c2*j2*m2*u0*w0*pow(x,2) - e2*i2*q2*u0*w0*pow(x,2) + 
   pow(m2,2)*q2*u0*w0*pow(x,2) + c2*i2*r2*u0*w0*pow(x,2) - 
   l2*m2*r2*u0*w0*pow(x,2) + e0*j2*l2*u2*w0*pow(x,2) - 
   c0*j2*m2*u2*w0*pow(x,2) - e0*i2*q2*u2*w0*pow(x,2) + 
   c2*i2*r0*u2*w0*pow(x,2) - l2*m2*r0*u2*w0*pow(x,2) + 
   c0*i2*r2*u2*w0*pow(x,2) - c2*j2*l2*v0*w0*pow(x,2) + 
   c2*i2*q2*v0*w0*pow(x,2) - l2*m2*q2*v0*w0*pow(x,2) + 
   pow(l2,2)*r2*v0*w0*pow(x,2) - c0*j2*l2*v2*w0*pow(x,2) + 
   c0*i2*q2*v2*w0*pow(x,2) + pow(l2,2)*r0*v2*w0*pow(x,2) + 
   pow(c0,2)*j2*k2*w2*pow(x,2) - pow(c0,2)*i2*p2*w2*pow(x,2) - 
   c0*k2*l2*r0*w2*pow(x,2) + e0*j2*l2*u0*w2*pow(x,2) - 
   c0*j2*m2*u0*w2*pow(x,2) - e0*i2*q2*u0*w2*pow(x,2) + 
   c2*i2*r0*u0*w2*pow(x,2) - l2*m2*r0*u0*w2*pow(x,2) + 
   c0*i2*r2*u0*w2*pow(x,2) + c0*i2*r0*u2*w2*pow(x,2) - 
   c0*j2*l2*v0*w2*pow(x,2) + c0*i2*q2*v0*w2*pow(x,2) + 
   pow(l2,2)*r0*v0*w2*pow(x,2) + pow(c0,2)*k2*n0*p2*pow(x,3) + 
   c0*d0*pow(k2,2)*r0*pow(x,3) - pow(c0,2)*pow(k2,2)*s0*pow(x,3) - 
   pow(c0,2)*j2*n0*t2*pow(x,3) - c0*d0*i2*r0*t2*pow(x,3) + 
   c0*l2*n0*r0*t2*pow(x,3) + pow(c0,2)*i2*s0*t2*pow(x,3) + 
   d0*e0*j2*k2*u0*pow(x,3) - c0*f0*j2*k2*u0*pow(x,3) - 
   d0*e0*i2*p2*u0*pow(x,3) + c0*f0*i2*p2*u0*pow(x,3) + 
   e0*l2*n0*p2*u0*pow(x,3) - c0*m2*n0*p2*u0*pow(x,3) + 
   e0*k2*n0*q2*u0*pow(x,3) + 2*f0*k2*l2*r0*u0*pow(x,3) - 
   d0*k2*m2*r0*u0*pow(x,3) - c2*k2*n0*r0*u0*pow(x,3) - 
   c0*k2*n2*r0*u0*pow(x,3) - c0*k2*n0*r2*u0*pow(x,3) - 
   2*e0*k2*l2*s0*u0*pow(x,3) + 2*c0*k2*m2*s0*u0*pow(x,3) + 
   f0*j2*m2*pow(u0,2)*pow(x,3) - e2*j2*n0*pow(u0,2)*pow(x,3) - 
   e0*j2*n2*pow(u0,2)*pow(x,3) - f2*i2*r0*pow(u0,2)*pow(x,3) + 
   m2*n2*r0*pow(u0,2)*pow(x,3) - f0*i2*r2*pow(u0,2)*pow(x,3) + 
   m2*n0*r2*pow(u0,2)*pow(x,3) + e2*i2*s0*pow(u0,2)*pow(x,3) - 
   pow(m2,2)*s0*pow(u0,2)*pow(x,3) + e0*i2*s2*pow(u0,2)*pow(x,3) - 
   c0*k2*n0*r0*u2*pow(x,3) - 2*e0*j2*n0*u0*u2*pow(x,3) - 
   2*f0*i2*r0*u0*u2*pow(x,3) + 2*m2*n0*r0*u0*u2*pow(x,3) + 
   2*e0*i2*s0*u0*u2*pow(x,3) - c0*d0*j2*k2*v0*pow(x,3) + 
   c0*d0*i2*p2*v0*pow(x,3) - c0*l2*n0*p2*v0*pow(x,3) - 
   c0*k2*n0*q2*v0*pow(x,3) - d0*k2*l2*r0*v0*pow(x,3) + 
   2*c0*k2*l2*s0*v0*pow(x,3) - f0*j2*l2*u0*v0*pow(x,3) - 
   d0*j2*m2*u0*v0*pow(x,3) + 2*c2*j2*n0*u0*v0*pow(x,3) + 
   2*c0*j2*n2*u0*v0*pow(x,3) + f0*i2*q2*u0*v0*pow(x,3) - 
   m2*n0*q2*u0*v0*pow(x,3) + d2*i2*r0*u0*v0*pow(x,3) - 
   l2*n2*r0*u0*v0*pow(x,3) + d0*i2*r2*u0*v0*pow(x,3) - 
   l2*n0*r2*u0*v0*pow(x,3) - 2*c2*i2*s0*u0*v0*pow(x,3) + 
   2*l2*m2*s0*u0*v0*pow(x,3) - 2*c0*i2*s2*u0*v0*pow(x,3) + 
   2*c0*j2*n0*u2*v0*pow(x,3) + d0*i2*r0*u2*v0*pow(x,3) - 
   l2*n0*r0*u2*v0*pow(x,3) - 2*c0*i2*s0*u2*v0*pow(x,3) + 
   d0*j2*l2*pow(v0,2)*pow(x,3) - d0*i2*q2*pow(v0,2)*pow(x,3) + 
   l2*n0*q2*pow(v0,2)*pow(x,3) - pow(l2,2)*s0*pow(v0,2)*pow(x,3) + 
   2*c0*j2*n0*u0*v2*pow(x,3) + d0*i2*r0*u0*v2*pow(x,3) - 
   l2*n0*r0*u0*v2*pow(x,3) - 2*c0*i2*s0*u0*v2*pow(x,3) + 
   pow(c0,2)*j2*k2*w0*pow(x,3) - pow(c0,2)*i2*p2*w0*pow(x,3) - 
   c0*k2*l2*r0*w0*pow(x,3) + e0*j2*l2*u0*w0*pow(x,3) - 
   c0*j2*m2*u0*w0*pow(x,3) - e0*i2*q2*u0*w0*pow(x,3) + 
   c2*i2*r0*u0*w0*pow(x,3) - l2*m2*r0*u0*w0*pow(x,3) + 
   c0*i2*r2*u0*w0*pow(x,3) + c0*i2*r0*u2*w0*pow(x,3) - 
   c0*j2*l2*v0*w0*pow(x,3) + c0*i2*q2*v0*w0*pow(x,3) + 
   pow(l2,2)*r0*v0*w0*pow(x,3) + c0*i2*r0*u0*w2*pow(x,3) - 
   c0*k2*n0*r0*u0*pow(x,4) - e0*j2*n0*pow(u0,2)*pow(x,4) - 
   f0*i2*r0*pow(u0,2)*pow(x,4) + m2*n0*r0*pow(u0,2)*pow(x,4) + 
   e0*i2*s0*pow(u0,2)*pow(x,4) + 2*c0*j2*n0*u0*v0*pow(x,4) + 
   d0*i2*r0*u0*v0*pow(x,4) - l2*n0*r0*u0*v0*pow(x,4) - 
   2*c0*i2*s0*u0*v0*pow(x,4) + c0*i2*r0*u0*w0*pow(x,4) + 
   d2*e1*k2*l2*p2*y + d1*e2*k2*l2*p2*y - c2*f1*k2*l2*p2*y - 
   c1*f2*k2*l2*p2*y - c2*d1*k2*m2*p2*y - c1*d2*k2*m2*p2*y + 
   pow(c2,2)*k2*n1*p2*y + 2*c1*c2*k2*n2*p2*y - d2*e1*pow(k2,2)*q2*y - 
   d1*e2*pow(k2,2)*q2*y + c2*f1*pow(k2,2)*q2*y + c1*f2*pow(k2,2)*q2*y + 
   c2*d2*pow(k2,2)*r1*y + c2*d1*pow(k2,2)*r2*y + c1*d2*pow(k2,2)*r2*y - 
   pow(c2,2)*pow(k2,2)*s1*y - 2*c1*c2*pow(k2,2)*s2*y - 
   d2*e1*j2*l2*t2*y - d1*e2*j2*l2*t2*y + c2*f1*j2*l2*t2*y + 
   c1*f2*j2*l2*t2*y + c2*d1*j2*m2*t2*y + c1*d2*j2*m2*t2*y - 
   pow(c2,2)*j2*n1*t2*y - 2*c1*c2*j2*n2*t2*y + d2*e1*i2*q2*t2*y + 
   d1*e2*i2*q2*t2*y - c2*f1*i2*q2*t2*y - c1*f2*i2*q2*t2*y + 
   f1*l2*m2*q2*t2*y - d1*pow(m2,2)*q2*t2*y - e2*l2*n1*q2*t2*y + 
   c2*m2*n1*q2*t2*y - e1*l2*n2*q2*t2*y + c1*m2*n2*q2*t2*y - 
   c2*d2*i2*r1*t2*y - f2*pow(l2,2)*r1*t2*y + d2*l2*m2*r1*t2*y + 
   c2*l2*n2*r1*t2*y - c2*d1*i2*r2*t2*y - c1*d2*i2*r2*t2*y - 
   f1*pow(l2,2)*r2*t2*y + d1*l2*m2*r2*t2*y + c2*l2*n1*r2*t2*y + 
   c1*l2*n2*r2*t2*y + pow(c2,2)*i2*s1*t2*y + e2*pow(l2,2)*s1*t2*y - 
   2*c2*l2*m2*s1*t2*y + 2*c1*c2*i2*s2*t2*y + e1*pow(l2,2)*s2*t2*y - 
   2*c1*l2*m2*s2*t2*y + d2*e2*j2*k2*u1*y - c2*f2*j2*k2*u1*y - 
   d2*e2*i2*p2*u1*y + c2*f2*i2*p2*u1*y - f2*l2*m2*p2*u1*y + 
   d2*pow(m2,2)*p2*u1*y + e2*l2*n2*p2*u1*y - c2*m2*n2*p2*u1*y - 
   f2*k2*m2*q2*u1*y + e2*k2*n2*q2*u1*y + 2*f2*k2*l2*r2*u1*y - 
   d2*k2*m2*r2*u1*y - c2*k2*n2*r2*u1*y - 2*e2*k2*l2*s2*u1*y + 
   2*c2*k2*m2*s2*u1*y + d2*e1*j2*k2*u2*y + d1*e2*j2*k2*u2*y - 
   c2*f1*j2*k2*u2*y - c1*f2*j2*k2*u2*y - d2*e1*i2*p2*u2*y - 
   d1*e2*i2*p2*u2*y + c2*f1*i2*p2*u2*y + c1*f2*i2*p2*u2*y - 
   f1*l2*m2*p2*u2*y + d1*pow(m2,2)*p2*u2*y + e2*l2*n1*p2*u2*y - 
   c2*m2*n1*p2*u2*y + e1*l2*n2*p2*u2*y - c1*m2*n2*p2*u2*y - 
   f1*k2*m2*q2*u2*y + e2*k2*n1*q2*u2*y + e1*k2*n2*q2*u2*y + 
   2*f2*k2*l2*r1*u2*y - d2*k2*m2*r1*u2*y - c2*k2*n2*r1*u2*y + 
   2*f1*k2*l2*r2*u2*y - d1*k2*m2*r2*u2*y - c2*k2*n1*r2*u2*y - 
   c1*k2*n2*r2*u2*y - 2*e2*k2*l2*s1*u2*y + 2*c2*k2*m2*s1*u2*y - 
   2*e1*k2*l2*s2*u2*y + 2*c1*k2*m2*s2*u2*y + 2*f2*j2*m2*u1*u2*y - 
   2*e2*j2*n2*u1*u2*y - 2*f2*i2*r2*u1*u2*y + 2*m2*n2*r2*u1*u2*y + 
   2*e2*i2*s2*u1*u2*y - 2*pow(m2,2)*s2*u1*u2*y + f1*j2*m2*pow(u2,2)*y - 
   e2*j2*n1*pow(u2,2)*y - e1*j2*n2*pow(u2,2)*y - f2*i2*r1*pow(u2,2)*y + 
   m2*n2*r1*pow(u2,2)*y - f1*i2*r2*pow(u2,2)*y + m2*n1*r2*pow(u2,2)*y + 
   e2*i2*s1*pow(u2,2)*y - pow(m2,2)*s1*pow(u2,2)*y + 
   e1*i2*s2*pow(u2,2)*y - c2*d2*j2*k2*v1*y + c2*d2*i2*p2*v1*y + 
   f2*pow(l2,2)*p2*v1*y - d2*l2*m2*p2*v1*y - c2*l2*n2*p2*v1*y - 
   f2*k2*l2*q2*v1*y + 2*d2*k2*m2*q2*v1*y - c2*k2*n2*q2*v1*y - 
   d2*k2*l2*r2*v1*y + 2*c2*k2*l2*s2*v1*y - f2*j2*l2*u2*v1*y - 
   d2*j2*m2*u2*v1*y + 2*c2*j2*n2*u2*v1*y + f2*i2*q2*u2*v1*y - 
   m2*n2*q2*u2*v1*y + d2*i2*r2*u2*v1*y - l2*n2*r2*u2*v1*y - 
   2*c2*i2*s2*u2*v1*y + 2*l2*m2*s2*u2*v1*y - c2*d1*j2*k2*v2*y - 
   c1*d2*j2*k2*v2*y + c2*d1*i2*p2*v2*y + c1*d2*i2*p2*v2*y + 
   f1*pow(l2,2)*p2*v2*y - d1*l2*m2*p2*v2*y - c2*l2*n1*p2*v2*y - 
   c1*l2*n2*p2*v2*y - f1*k2*l2*q2*v2*y + 2*d1*k2*m2*q2*v2*y - 
   c2*k2*n1*q2*v2*y - c1*k2*n2*q2*v2*y - d2*k2*l2*r1*v2*y - 
   d1*k2*l2*r2*v2*y + 2*c2*k2*l2*s1*v2*y + 2*c1*k2*l2*s2*v2*y - 
   f2*j2*l2*u1*v2*y - d2*j2*m2*u1*v2*y + 2*c2*j2*n2*u1*v2*y + 
   f2*i2*q2*u1*v2*y - m2*n2*q2*u1*v2*y + d2*i2*r2*u1*v2*y - 
   l2*n2*r2*u1*v2*y - 2*c2*i2*s2*u1*v2*y + 2*l2*m2*s2*u1*v2*y - 
   f1*j2*l2*u2*v2*y - d1*j2*m2*u2*v2*y + 2*c2*j2*n1*u2*v2*y + 
   2*c1*j2*n2*u2*v2*y + f1*i2*q2*u2*v2*y - m2*n1*q2*u2*v2*y + 
   d2*i2*r1*u2*v2*y - l2*n2*r1*u2*v2*y + d1*i2*r2*u2*v2*y - 
   l2*n1*r2*u2*v2*y - 2*c2*i2*s1*u2*v2*y + 2*l2*m2*s1*u2*v2*y - 
   2*c1*i2*s2*u2*v2*y + 2*d2*j2*l2*v1*v2*y - 2*d2*i2*q2*v1*v2*y + 
   2*l2*n2*q2*v1*v2*y - 2*pow(l2,2)*s2*v1*v2*y + d1*j2*l2*pow(v2,2)*y - 
   d1*i2*q2*pow(v2,2)*y + l2*n1*q2*pow(v2,2)*y - 
   pow(l2,2)*s1*pow(v2,2)*y + pow(c2,2)*j2*k2*w1*y - 
   pow(c2,2)*i2*p2*w1*y - e2*pow(l2,2)*p2*w1*y + 2*c2*l2*m2*p2*w1*y + 
   e2*k2*l2*q2*w1*y - c2*k2*m2*q2*w1*y - c2*k2*l2*r2*w1*y + 
   e2*j2*l2*u2*w1*y - c2*j2*m2*u2*w1*y - e2*i2*q2*u2*w1*y + 
   pow(m2,2)*q2*u2*w1*y + c2*i2*r2*u2*w1*y - l2*m2*r2*u2*w1*y - 
   c2*j2*l2*v2*w1*y + c2*i2*q2*v2*w1*y - l2*m2*q2*v2*w1*y + 
   pow(l2,2)*r2*v2*w1*y + 2*c1*c2*j2*k2*w2*y - 2*c1*c2*i2*p2*w2*y - 
   e1*pow(l2,2)*p2*w2*y + 2*c1*l2*m2*p2*w2*y + e1*k2*l2*q2*w2*y - 
   c1*k2*m2*q2*w2*y - c2*k2*l2*r1*w2*y - c1*k2*l2*r2*w2*y + 
   e2*j2*l2*u1*w2*y - c2*j2*m2*u1*w2*y - e2*i2*q2*u1*w2*y + 
   pow(m2,2)*q2*u1*w2*y + c2*i2*r2*u1*w2*y - l2*m2*r2*u1*w2*y + 
   e1*j2*l2*u2*w2*y - c1*j2*m2*u2*w2*y - e1*i2*q2*u2*w2*y + 
   c2*i2*r1*u2*w2*y - l2*m2*r1*u2*w2*y + c1*i2*r2*u2*w2*y - 
   c2*j2*l2*v1*w2*y + c2*i2*q2*v1*w2*y - l2*m2*q2*v1*w2*y + 
   pow(l2,2)*r2*v1*w2*y - c1*j2*l2*v2*w2*y + c1*i2*q2*v2*w2*y + 
   pow(l2,2)*r1*v2*w2*y + d1*e0*k2*l2*p2*x*y + d0*e1*k2*l2*p2*x*y - 
   c1*f0*k2*l2*p2*x*y - c0*f1*k2*l2*p2*x*y - c1*d0*k2*m2*p2*x*y - 
   c0*d1*k2*m2*p2*x*y + 2*c1*c2*k2*n0*p2*x*y + 2*c0*c2*k2*n1*p2*x*y + 
   2*c0*c1*k2*n2*p2*x*y - d1*e0*pow(k2,2)*q2*x*y - 
   d0*e1*pow(k2,2)*q2*x*y + c1*f0*pow(k2,2)*q2*x*y + 
   c0*f1*pow(k2,2)*q2*x*y + c2*d1*pow(k2,2)*r0*x*y + 
   c1*d2*pow(k2,2)*r0*x*y + c2*d0*pow(k2,2)*r1*x*y + 
   c0*d2*pow(k2,2)*r1*x*y + c1*d0*pow(k2,2)*r2*x*y + 
   c0*d1*pow(k2,2)*r2*x*y - 2*c1*c2*pow(k2,2)*s0*x*y - 
   2*c0*c2*pow(k2,2)*s1*x*y - 2*c0*c1*pow(k2,2)*s2*x*y - 
   d1*e0*j2*l2*t2*x*y - d0*e1*j2*l2*t2*x*y + c1*f0*j2*l2*t2*x*y + 
   c0*f1*j2*l2*t2*x*y + c1*d0*j2*m2*t2*x*y + c0*d1*j2*m2*t2*x*y - 
   2*c1*c2*j2*n0*t2*x*y - 2*c0*c2*j2*n1*t2*x*y - 2*c0*c1*j2*n2*t2*x*y + 
   d1*e0*i2*q2*t2*x*y + d0*e1*i2*q2*t2*x*y - c1*f0*i2*q2*t2*x*y - 
   c0*f1*i2*q2*t2*x*y - e1*l2*n0*q2*t2*x*y + c1*m2*n0*q2*t2*x*y - 
   e0*l2*n1*q2*t2*x*y + c0*m2*n1*q2*t2*x*y - c2*d1*i2*r0*t2*x*y - 
   c1*d2*i2*r0*t2*x*y - f1*pow(l2,2)*r0*t2*x*y + d1*l2*m2*r0*t2*x*y + 
   c2*l2*n1*r0*t2*x*y + c1*l2*n2*r0*t2*x*y - c2*d0*i2*r1*t2*x*y - 
   c0*d2*i2*r1*t2*x*y - f0*pow(l2,2)*r1*t2*x*y + d0*l2*m2*r1*t2*x*y + 
   c2*l2*n0*r1*t2*x*y + c0*l2*n2*r1*t2*x*y - c1*d0*i2*r2*t2*x*y - 
   c0*d1*i2*r2*t2*x*y + c1*l2*n0*r2*t2*x*y + c0*l2*n1*r2*t2*x*y + 
   2*c1*c2*i2*s0*t2*x*y + e1*pow(l2,2)*s0*t2*x*y - 2*c1*l2*m2*s0*t2*x*y + 
   2*c0*c2*i2*s1*t2*x*y + e0*pow(l2,2)*s1*t2*x*y - 2*c0*l2*m2*s1*t2*x*y + 
   2*c0*c1*i2*s2*t2*x*y + d2*e1*j2*k2*u0*x*y + d1*e2*j2*k2*u0*x*y - 
   c2*f1*j2*k2*u0*x*y - c1*f2*j2*k2*u0*x*y - d2*e1*i2*p2*u0*x*y - 
   d1*e2*i2*p2*u0*x*y + c2*f1*i2*p2*u0*x*y + c1*f2*i2*p2*u0*x*y - 
   f1*l2*m2*p2*u0*x*y + d1*pow(m2,2)*p2*u0*x*y + e2*l2*n1*p2*u0*x*y - 
   c2*m2*n1*p2*u0*x*y + e1*l2*n2*p2*u0*x*y - c1*m2*n2*p2*u0*x*y - 
   f1*k2*m2*q2*u0*x*y + e2*k2*n1*q2*u0*x*y + e1*k2*n2*q2*u0*x*y + 
   2*f2*k2*l2*r1*u0*x*y - d2*k2*m2*r1*u0*x*y - c2*k2*n2*r1*u0*x*y + 
   2*f1*k2*l2*r2*u0*x*y - d1*k2*m2*r2*u0*x*y - c2*k2*n1*r2*u0*x*y - 
   c1*k2*n2*r2*u0*x*y - 2*e2*k2*l2*s1*u0*x*y + 2*c2*k2*m2*s1*u0*x*y - 
   2*e1*k2*l2*s2*u0*x*y + 2*c1*k2*m2*s2*u0*x*y + d2*e0*j2*k2*u1*x*y + 
   d0*e2*j2*k2*u1*x*y - c2*f0*j2*k2*u1*x*y - c0*f2*j2*k2*u1*x*y - 
   d2*e0*i2*p2*u1*x*y - d0*e2*i2*p2*u1*x*y + c2*f0*i2*p2*u1*x*y + 
   c0*f2*i2*p2*u1*x*y - f0*l2*m2*p2*u1*x*y + d0*pow(m2,2)*p2*u1*x*y + 
   e2*l2*n0*p2*u1*x*y - c2*m2*n0*p2*u1*x*y + e0*l2*n2*p2*u1*x*y - 
   c0*m2*n2*p2*u1*x*y - f0*k2*m2*q2*u1*x*y + e2*k2*n0*q2*u1*x*y + 
   e0*k2*n2*q2*u1*x*y + 2*f2*k2*l2*r0*u1*x*y - d2*k2*m2*r0*u1*x*y - 
   c2*k2*n2*r0*u1*x*y + 2*f0*k2*l2*r2*u1*x*y - d0*k2*m2*r2*u1*x*y - 
   c2*k2*n0*r2*u1*x*y - c0*k2*n2*r2*u1*x*y - 2*e2*k2*l2*s0*u1*x*y + 
   2*c2*k2*m2*s0*u1*x*y - 2*e0*k2*l2*s2*u1*x*y + 2*c0*k2*m2*s2*u1*x*y + 
   2*f2*j2*m2*u0*u1*x*y - 2*e2*j2*n2*u0*u1*x*y - 2*f2*i2*r2*u0*u1*x*y + 
   2*m2*n2*r2*u0*u1*x*y + 2*e2*i2*s2*u0*u1*x*y - 2*pow(m2,2)*s2*u0*u1*x*y + 
   d1*e0*j2*k2*u2*x*y + d0*e1*j2*k2*u2*x*y - c1*f0*j2*k2*u2*x*y - 
   c0*f1*j2*k2*u2*x*y - d1*e0*i2*p2*u2*x*y - d0*e1*i2*p2*u2*x*y + 
   c1*f0*i2*p2*u2*x*y + c0*f1*i2*p2*u2*x*y + e1*l2*n0*p2*u2*x*y - 
   c1*m2*n0*p2*u2*x*y + e0*l2*n1*p2*u2*x*y - c0*m2*n1*p2*u2*x*y + 
   e1*k2*n0*q2*u2*x*y + e0*k2*n1*q2*u2*x*y + 2*f1*k2*l2*r0*u2*x*y - 
   d1*k2*m2*r0*u2*x*y - c2*k2*n1*r0*u2*x*y - c1*k2*n2*r0*u2*x*y + 
   2*f0*k2*l2*r1*u2*x*y - d0*k2*m2*r1*u2*x*y - c2*k2*n0*r1*u2*x*y - 
   c0*k2*n2*r1*u2*x*y - c1*k2*n0*r2*u2*x*y - c0*k2*n1*r2*u2*x*y - 
   2*e1*k2*l2*s0*u2*x*y + 2*c1*k2*m2*s0*u2*x*y - 2*e0*k2*l2*s1*u2*x*y + 
   2*c0*k2*m2*s1*u2*x*y + 2*f1*j2*m2*u0*u2*x*y - 2*e2*j2*n1*u0*u2*x*y - 
   2*e1*j2*n2*u0*u2*x*y - 2*f2*i2*r1*u0*u2*x*y + 2*m2*n2*r1*u0*u2*x*y - 
   2*f1*i2*r2*u0*u2*x*y + 2*m2*n1*r2*u0*u2*x*y + 2*e2*i2*s1*u0*u2*x*y - 
   2*pow(m2,2)*s1*u0*u2*x*y + 2*e1*i2*s2*u0*u2*x*y + 2*f0*j2*m2*u1*u2*x*y - 
   2*e2*j2*n0*u1*u2*x*y - 2*e0*j2*n2*u1*u2*x*y - 2*f2*i2*r0*u1*u2*x*y + 
   2*m2*n2*r0*u1*u2*x*y - 2*f0*i2*r2*u1*u2*x*y + 2*m2*n0*r2*u1*u2*x*y + 
   2*e2*i2*s0*u1*u2*x*y - 2*pow(m2,2)*s0*u1*u2*x*y + 2*e0*i2*s2*u1*u2*x*y - 
   e1*j2*n0*pow(u2,2)*x*y - e0*j2*n1*pow(u2,2)*x*y - 
   f1*i2*r0*pow(u2,2)*x*y + m2*n1*r0*pow(u2,2)*x*y - 
   f0*i2*r1*pow(u2,2)*x*y + m2*n0*r1*pow(u2,2)*x*y + 
   e1*i2*s0*pow(u2,2)*x*y + e0*i2*s1*pow(u2,2)*x*y - c2*d1*j2*k2*v0*x*y - 
   c1*d2*j2*k2*v0*x*y + c2*d1*i2*p2*v0*x*y + c1*d2*i2*p2*v0*x*y + 
   f1*pow(l2,2)*p2*v0*x*y - d1*l2*m2*p2*v0*x*y - c2*l2*n1*p2*v0*x*y - 
   c1*l2*n2*p2*v0*x*y - f1*k2*l2*q2*v0*x*y + 2*d1*k2*m2*q2*v0*x*y - 
   c2*k2*n1*q2*v0*x*y - c1*k2*n2*q2*v0*x*y - d2*k2*l2*r1*v0*x*y - 
   d1*k2*l2*r2*v0*x*y + 2*c2*k2*l2*s1*v0*x*y + 2*c1*k2*l2*s2*v0*x*y - 
   f2*j2*l2*u1*v0*x*y - d2*j2*m2*u1*v0*x*y + 2*c2*j2*n2*u1*v0*x*y + 
   f2*i2*q2*u1*v0*x*y - m2*n2*q2*u1*v0*x*y + d2*i2*r2*u1*v0*x*y - 
   l2*n2*r2*u1*v0*x*y - 2*c2*i2*s2*u1*v0*x*y + 2*l2*m2*s2*u1*v0*x*y - 
   f1*j2*l2*u2*v0*x*y - d1*j2*m2*u2*v0*x*y + 2*c2*j2*n1*u2*v0*x*y + 
   2*c1*j2*n2*u2*v0*x*y + f1*i2*q2*u2*v0*x*y - m2*n1*q2*u2*v0*x*y + 
   d2*i2*r1*u2*v0*x*y - l2*n2*r1*u2*v0*x*y + d1*i2*r2*u2*v0*x*y - 
   l2*n1*r2*u2*v0*x*y - 2*c2*i2*s1*u2*v0*x*y + 2*l2*m2*s1*u2*v0*x*y - 
   2*c1*i2*s2*u2*v0*x*y - c2*d0*j2*k2*v1*x*y - c0*d2*j2*k2*v1*x*y + 
   c2*d0*i2*p2*v1*x*y + c0*d2*i2*p2*v1*x*y + f0*pow(l2,2)*p2*v1*x*y - 
   d0*l2*m2*p2*v1*x*y - c2*l2*n0*p2*v1*x*y - c0*l2*n2*p2*v1*x*y - 
   f0*k2*l2*q2*v1*x*y + 2*d0*k2*m2*q2*v1*x*y - c2*k2*n0*q2*v1*x*y - 
   c0*k2*n2*q2*v1*x*y - d2*k2*l2*r0*v1*x*y - d0*k2*l2*r2*v1*x*y + 
   2*c2*k2*l2*s0*v1*x*y + 2*c0*k2*l2*s2*v1*x*y - f2*j2*l2*u0*v1*x*y - 
   d2*j2*m2*u0*v1*x*y + 2*c2*j2*n2*u0*v1*x*y + f2*i2*q2*u0*v1*x*y - 
   m2*n2*q2*u0*v1*x*y + d2*i2*r2*u0*v1*x*y - l2*n2*r2*u0*v1*x*y - 
   2*c2*i2*s2*u0*v1*x*y + 2*l2*m2*s2*u0*v1*x*y - f0*j2*l2*u2*v1*x*y - 
   d0*j2*m2*u2*v1*x*y + 2*c2*j2*n0*u2*v1*x*y + 2*c0*j2*n2*u2*v1*x*y + 
   f0*i2*q2*u2*v1*x*y - m2*n0*q2*u2*v1*x*y + d2*i2*r0*u2*v1*x*y - 
   l2*n2*r0*u2*v1*x*y + d0*i2*r2*u2*v1*x*y - l2*n0*r2*u2*v1*x*y - 
   2*c2*i2*s0*u2*v1*x*y + 2*l2*m2*s0*u2*v1*x*y - 2*c0*i2*s2*u2*v1*x*y + 
   2*d2*j2*l2*v0*v1*x*y - 2*d2*i2*q2*v0*v1*x*y + 2*l2*n2*q2*v0*v1*x*y - 
   2*pow(l2,2)*s2*v0*v1*x*y - c1*d0*j2*k2*v2*x*y - c0*d1*j2*k2*v2*x*y + 
   c1*d0*i2*p2*v2*x*y + c0*d1*i2*p2*v2*x*y - c1*l2*n0*p2*v2*x*y - 
   c0*l2*n1*p2*v2*x*y - c1*k2*n0*q2*v2*x*y - c0*k2*n1*q2*v2*x*y - 
   d1*k2*l2*r0*v2*x*y - d0*k2*l2*r1*v2*x*y + 2*c1*k2*l2*s0*v2*x*y + 
   2*c0*k2*l2*s1*v2*x*y - f1*j2*l2*u0*v2*x*y - d1*j2*m2*u0*v2*x*y + 
   2*c2*j2*n1*u0*v2*x*y + 2*c1*j2*n2*u0*v2*x*y + f1*i2*q2*u0*v2*x*y - 
   m2*n1*q2*u0*v2*x*y + d2*i2*r1*u0*v2*x*y - l2*n2*r1*u0*v2*x*y + 
   d1*i2*r2*u0*v2*x*y - l2*n1*r2*u0*v2*x*y - 2*c2*i2*s1*u0*v2*x*y + 
   2*l2*m2*s1*u0*v2*x*y - 2*c1*i2*s2*u0*v2*x*y - f0*j2*l2*u1*v2*x*y - 
   d0*j2*m2*u1*v2*x*y + 2*c2*j2*n0*u1*v2*x*y + 2*c0*j2*n2*u1*v2*x*y + 
   f0*i2*q2*u1*v2*x*y - m2*n0*q2*u1*v2*x*y + d2*i2*r0*u1*v2*x*y - 
   l2*n2*r0*u1*v2*x*y + d0*i2*r2*u1*v2*x*y - l2*n0*r2*u1*v2*x*y - 
   2*c2*i2*s0*u1*v2*x*y + 2*l2*m2*s0*u1*v2*x*y - 2*c0*i2*s2*u1*v2*x*y + 
   2*c1*j2*n0*u2*v2*x*y + 2*c0*j2*n1*u2*v2*x*y + d1*i2*r0*u2*v2*x*y - 
   l2*n1*r0*u2*v2*x*y + d0*i2*r1*u2*v2*x*y - l2*n0*r1*u2*v2*x*y - 
   2*c1*i2*s0*u2*v2*x*y - 2*c0*i2*s1*u2*v2*x*y + 2*d1*j2*l2*v0*v2*x*y - 
   2*d1*i2*q2*v0*v2*x*y + 2*l2*n1*q2*v0*v2*x*y - 2*pow(l2,2)*s1*v0*v2*x*y + 
   2*d0*j2*l2*v1*v2*x*y - 2*d0*i2*q2*v1*v2*x*y + 2*l2*n0*q2*v1*v2*x*y - 
   2*pow(l2,2)*s0*v1*v2*x*y + 2*c1*c2*j2*k2*w0*x*y - 2*c1*c2*i2*p2*w0*x*y - 
   e1*pow(l2,2)*p2*w0*x*y + 2*c1*l2*m2*p2*w0*x*y + e1*k2*l2*q2*w0*x*y - 
   c1*k2*m2*q2*w0*x*y - c2*k2*l2*r1*w0*x*y - c1*k2*l2*r2*w0*x*y + 
   e2*j2*l2*u1*w0*x*y - c2*j2*m2*u1*w0*x*y - e2*i2*q2*u1*w0*x*y + 
   pow(m2,2)*q2*u1*w0*x*y + c2*i2*r2*u1*w0*x*y - l2*m2*r2*u1*w0*x*y + 
   e1*j2*l2*u2*w0*x*y - c1*j2*m2*u2*w0*x*y - e1*i2*q2*u2*w0*x*y + 
   c2*i2*r1*u2*w0*x*y - l2*m2*r1*u2*w0*x*y + c1*i2*r2*u2*w0*x*y - 
   c2*j2*l2*v1*w0*x*y + c2*i2*q2*v1*w0*x*y - l2*m2*q2*v1*w0*x*y + 
   pow(l2,2)*r2*v1*w0*x*y - c1*j2*l2*v2*w0*x*y + c1*i2*q2*v2*w0*x*y + 
   pow(l2,2)*r1*v2*w0*x*y + 2*c0*c2*j2*k2*w1*x*y - 2*c0*c2*i2*p2*w1*x*y - 
   e0*pow(l2,2)*p2*w1*x*y + 2*c0*l2*m2*p2*w1*x*y + e0*k2*l2*q2*w1*x*y - 
   c0*k2*m2*q2*w1*x*y - c2*k2*l2*r0*w1*x*y - c0*k2*l2*r2*w1*x*y + 
   e2*j2*l2*u0*w1*x*y - c2*j2*m2*u0*w1*x*y - e2*i2*q2*u0*w1*x*y + 
   pow(m2,2)*q2*u0*w1*x*y + c2*i2*r2*u0*w1*x*y - l2*m2*r2*u0*w1*x*y + 
   e0*j2*l2*u2*w1*x*y - c0*j2*m2*u2*w1*x*y - e0*i2*q2*u2*w1*x*y + 
   c2*i2*r0*u2*w1*x*y - l2*m2*r0*u2*w1*x*y + c0*i2*r2*u2*w1*x*y - 
   c2*j2*l2*v0*w1*x*y + c2*i2*q2*v0*w1*x*y - l2*m2*q2*v0*w1*x*y + 
   pow(l2,2)*r2*v0*w1*x*y - c0*j2*l2*v2*w1*x*y + c0*i2*q2*v2*w1*x*y + 
   pow(l2,2)*r0*v2*w1*x*y + 2*c0*c1*j2*k2*w2*x*y - 2*c0*c1*i2*p2*w2*x*y - 
   c1*k2*l2*r0*w2*x*y - c0*k2*l2*r1*w2*x*y + e1*j2*l2*u0*w2*x*y - 
   c1*j2*m2*u0*w2*x*y - e1*i2*q2*u0*w2*x*y + c2*i2*r1*u0*w2*x*y - 
   l2*m2*r1*u0*w2*x*y + c1*i2*r2*u0*w2*x*y + e0*j2*l2*u1*w2*x*y - 
   c0*j2*m2*u1*w2*x*y - e0*i2*q2*u1*w2*x*y + c2*i2*r0*u1*w2*x*y - 
   l2*m2*r0*u1*w2*x*y + c0*i2*r2*u1*w2*x*y + c1*i2*r0*u2*w2*x*y + 
   c0*i2*r1*u2*w2*x*y - c1*j2*l2*v0*w2*x*y + c1*i2*q2*v0*w2*x*y + 
   pow(l2,2)*r1*v0*w2*x*y - c0*j2*l2*v1*w2*x*y + c0*i2*q2*v1*w2*x*y + 
   pow(l2,2)*r0*v1*w2*x*y + 2*c0*c1*k2*n0*p2*pow(x,2)*y + 
   pow(c0,2)*k2*n1*p2*pow(x,2)*y + c1*d0*pow(k2,2)*r0*pow(x,2)*y + 
   c0*d1*pow(k2,2)*r0*pow(x,2)*y + c0*d0*pow(k2,2)*r1*pow(x,2)*y - 
   2*c0*c1*pow(k2,2)*s0*pow(x,2)*y - 
   pow(c0,2)*pow(k2,2)*s1*pow(x,2)*y - 2*c0*c1*j2*n0*t2*pow(x,2)*y - 
   pow(c0,2)*j2*n1*t2*pow(x,2)*y - c1*d0*i2*r0*t2*pow(x,2)*y - 
   c0*d1*i2*r0*t2*pow(x,2)*y + c1*l2*n0*r0*t2*pow(x,2)*y + 
   c0*l2*n1*r0*t2*pow(x,2)*y - c0*d0*i2*r1*t2*pow(x,2)*y + 
   c0*l2*n0*r1*t2*pow(x,2)*y + 2*c0*c1*i2*s0*t2*pow(x,2)*y + 
   pow(c0,2)*i2*s1*t2*pow(x,2)*y + d1*e0*j2*k2*u0*pow(x,2)*y + 
   d0*e1*j2*k2*u0*pow(x,2)*y - c1*f0*j2*k2*u0*pow(x,2)*y - 
   c0*f1*j2*k2*u0*pow(x,2)*y - d1*e0*i2*p2*u0*pow(x,2)*y - 
   d0*e1*i2*p2*u0*pow(x,2)*y + c1*f0*i2*p2*u0*pow(x,2)*y + 
   c0*f1*i2*p2*u0*pow(x,2)*y + e1*l2*n0*p2*u0*pow(x,2)*y - 
   c1*m2*n0*p2*u0*pow(x,2)*y + e0*l2*n1*p2*u0*pow(x,2)*y - 
   c0*m2*n1*p2*u0*pow(x,2)*y + e1*k2*n0*q2*u0*pow(x,2)*y + 
   e0*k2*n1*q2*u0*pow(x,2)*y + 2*f1*k2*l2*r0*u0*pow(x,2)*y - 
   d1*k2*m2*r0*u0*pow(x,2)*y - c2*k2*n1*r0*u0*pow(x,2)*y - 
   c1*k2*n2*r0*u0*pow(x,2)*y + 2*f0*k2*l2*r1*u0*pow(x,2)*y - 
   d0*k2*m2*r1*u0*pow(x,2)*y - c2*k2*n0*r1*u0*pow(x,2)*y - 
   c0*k2*n2*r1*u0*pow(x,2)*y - c1*k2*n0*r2*u0*pow(x,2)*y - 
   c0*k2*n1*r2*u0*pow(x,2)*y - 2*e1*k2*l2*s0*u0*pow(x,2)*y + 
   2*c1*k2*m2*s0*u0*pow(x,2)*y - 2*e0*k2*l2*s1*u0*pow(x,2)*y + 
   2*c0*k2*m2*s1*u0*pow(x,2)*y + f1*j2*m2*pow(u0,2)*pow(x,2)*y - 
   e2*j2*n1*pow(u0,2)*pow(x,2)*y - e1*j2*n2*pow(u0,2)*pow(x,2)*y - 
   f2*i2*r1*pow(u0,2)*pow(x,2)*y + m2*n2*r1*pow(u0,2)*pow(x,2)*y - 
   f1*i2*r2*pow(u0,2)*pow(x,2)*y + m2*n1*r2*pow(u0,2)*pow(x,2)*y + 
   e2*i2*s1*pow(u0,2)*pow(x,2)*y - 
   pow(m2,2)*s1*pow(u0,2)*pow(x,2)*y + 
   e1*i2*s2*pow(u0,2)*pow(x,2)*y + d0*e0*j2*k2*u1*pow(x,2)*y - 
   c0*f0*j2*k2*u1*pow(x,2)*y - d0*e0*i2*p2*u1*pow(x,2)*y + 
   c0*f0*i2*p2*u1*pow(x,2)*y + e0*l2*n0*p2*u1*pow(x,2)*y - 
   c0*m2*n0*p2*u1*pow(x,2)*y + e0*k2*n0*q2*u1*pow(x,2)*y + 
   2*f0*k2*l2*r0*u1*pow(x,2)*y - d0*k2*m2*r0*u1*pow(x,2)*y - 
   c2*k2*n0*r0*u1*pow(x,2)*y - c0*k2*n2*r0*u1*pow(x,2)*y - 
   c0*k2*n0*r2*u1*pow(x,2)*y - 2*e0*k2*l2*s0*u1*pow(x,2)*y + 
   2*c0*k2*m2*s0*u1*pow(x,2)*y + 2*f0*j2*m2*u0*u1*pow(x,2)*y - 
   2*e2*j2*n0*u0*u1*pow(x,2)*y - 2*e0*j2*n2*u0*u1*pow(x,2)*y - 
   2*f2*i2*r0*u0*u1*pow(x,2)*y + 2*m2*n2*r0*u0*u1*pow(x,2)*y - 
   2*f0*i2*r2*u0*u1*pow(x,2)*y + 2*m2*n0*r2*u0*u1*pow(x,2)*y + 
   2*e2*i2*s0*u0*u1*pow(x,2)*y - 2*pow(m2,2)*s0*u0*u1*pow(x,2)*y + 
   2*e0*i2*s2*u0*u1*pow(x,2)*y - c1*k2*n0*r0*u2*pow(x,2)*y - 
   c0*k2*n1*r0*u2*pow(x,2)*y - c0*k2*n0*r1*u2*pow(x,2)*y - 
   2*e1*j2*n0*u0*u2*pow(x,2)*y - 2*e0*j2*n1*u0*u2*pow(x,2)*y - 
   2*f1*i2*r0*u0*u2*pow(x,2)*y + 2*m2*n1*r0*u0*u2*pow(x,2)*y - 
   2*f0*i2*r1*u0*u2*pow(x,2)*y + 2*m2*n0*r1*u0*u2*pow(x,2)*y + 
   2*e1*i2*s0*u0*u2*pow(x,2)*y + 2*e0*i2*s1*u0*u2*pow(x,2)*y - 
   2*e0*j2*n0*u1*u2*pow(x,2)*y - 2*f0*i2*r0*u1*u2*pow(x,2)*y + 
   2*m2*n0*r0*u1*u2*pow(x,2)*y + 2*e0*i2*s0*u1*u2*pow(x,2)*y - 
   c1*d0*j2*k2*v0*pow(x,2)*y - c0*d1*j2*k2*v0*pow(x,2)*y + 
   c1*d0*i2*p2*v0*pow(x,2)*y + c0*d1*i2*p2*v0*pow(x,2)*y - 
   c1*l2*n0*p2*v0*pow(x,2)*y - c0*l2*n1*p2*v0*pow(x,2)*y - 
   c1*k2*n0*q2*v0*pow(x,2)*y - c0*k2*n1*q2*v0*pow(x,2)*y - 
   d1*k2*l2*r0*v0*pow(x,2)*y - d0*k2*l2*r1*v0*pow(x,2)*y + 
   2*c1*k2*l2*s0*v0*pow(x,2)*y + 2*c0*k2*l2*s1*v0*pow(x,2)*y - 
   f1*j2*l2*u0*v0*pow(x,2)*y - d1*j2*m2*u0*v0*pow(x,2)*y + 
   2*c2*j2*n1*u0*v0*pow(x,2)*y + 2*c1*j2*n2*u0*v0*pow(x,2)*y + 
   f1*i2*q2*u0*v0*pow(x,2)*y - m2*n1*q2*u0*v0*pow(x,2)*y + 
   d2*i2*r1*u0*v0*pow(x,2)*y - l2*n2*r1*u0*v0*pow(x,2)*y + 
   d1*i2*r2*u0*v0*pow(x,2)*y - l2*n1*r2*u0*v0*pow(x,2)*y - 
   2*c2*i2*s1*u0*v0*pow(x,2)*y + 2*l2*m2*s1*u0*v0*pow(x,2)*y - 
   2*c1*i2*s2*u0*v0*pow(x,2)*y - f0*j2*l2*u1*v0*pow(x,2)*y - 
   d0*j2*m2*u1*v0*pow(x,2)*y + 2*c2*j2*n0*u1*v0*pow(x,2)*y + 
   2*c0*j2*n2*u1*v0*pow(x,2)*y + f0*i2*q2*u1*v0*pow(x,2)*y - 
   m2*n0*q2*u1*v0*pow(x,2)*y + d2*i2*r0*u1*v0*pow(x,2)*y - 
   l2*n2*r0*u1*v0*pow(x,2)*y + d0*i2*r2*u1*v0*pow(x,2)*y - 
   l2*n0*r2*u1*v0*pow(x,2)*y - 2*c2*i2*s0*u1*v0*pow(x,2)*y + 
   2*l2*m2*s0*u1*v0*pow(x,2)*y - 2*c0*i2*s2*u1*v0*pow(x,2)*y + 
   2*c1*j2*n0*u2*v0*pow(x,2)*y + 2*c0*j2*n1*u2*v0*pow(x,2)*y + 
   d1*i2*r0*u2*v0*pow(x,2)*y - l2*n1*r0*u2*v0*pow(x,2)*y + 
   d0*i2*r1*u2*v0*pow(x,2)*y - l2*n0*r1*u2*v0*pow(x,2)*y - 
   2*c1*i2*s0*u2*v0*pow(x,2)*y - 2*c0*i2*s1*u2*v0*pow(x,2)*y + 
   d1*j2*l2*pow(v0,2)*pow(x,2)*y - d1*i2*q2*pow(v0,2)*pow(x,2)*y + 
   l2*n1*q2*pow(v0,2)*pow(x,2)*y - 
   pow(l2,2)*s1*pow(v0,2)*pow(x,2)*y - c0*d0*j2*k2*v1*pow(x,2)*y + 
   c0*d0*i2*p2*v1*pow(x,2)*y - c0*l2*n0*p2*v1*pow(x,2)*y - 
   c0*k2*n0*q2*v1*pow(x,2)*y - d0*k2*l2*r0*v1*pow(x,2)*y + 
   2*c0*k2*l2*s0*v1*pow(x,2)*y - f0*j2*l2*u0*v1*pow(x,2)*y - 
   d0*j2*m2*u0*v1*pow(x,2)*y + 2*c2*j2*n0*u0*v1*pow(x,2)*y + 
   2*c0*j2*n2*u0*v1*pow(x,2)*y + f0*i2*q2*u0*v1*pow(x,2)*y - 
   m2*n0*q2*u0*v1*pow(x,2)*y + d2*i2*r0*u0*v1*pow(x,2)*y - 
   l2*n2*r0*u0*v1*pow(x,2)*y + d0*i2*r2*u0*v1*pow(x,2)*y - 
   l2*n0*r2*u0*v1*pow(x,2)*y - 2*c2*i2*s0*u0*v1*pow(x,2)*y + 
   2*l2*m2*s0*u0*v1*pow(x,2)*y - 2*c0*i2*s2*u0*v1*pow(x,2)*y + 
   2*c0*j2*n0*u2*v1*pow(x,2)*y + d0*i2*r0*u2*v1*pow(x,2)*y - 
   l2*n0*r0*u2*v1*pow(x,2)*y - 2*c0*i2*s0*u2*v1*pow(x,2)*y + 
   2*d0*j2*l2*v0*v1*pow(x,2)*y - 2*d0*i2*q2*v0*v1*pow(x,2)*y + 
   2*l2*n0*q2*v0*v1*pow(x,2)*y - 2*pow(l2,2)*s0*v0*v1*pow(x,2)*y + 
   2*c1*j2*n0*u0*v2*pow(x,2)*y + 2*c0*j2*n1*u0*v2*pow(x,2)*y + 
   d1*i2*r0*u0*v2*pow(x,2)*y - l2*n1*r0*u0*v2*pow(x,2)*y + 
   d0*i2*r1*u0*v2*pow(x,2)*y - l2*n0*r1*u0*v2*pow(x,2)*y - 
   2*c1*i2*s0*u0*v2*pow(x,2)*y - 2*c0*i2*s1*u0*v2*pow(x,2)*y + 
   2*c0*j2*n0*u1*v2*pow(x,2)*y + d0*i2*r0*u1*v2*pow(x,2)*y - 
   l2*n0*r0*u1*v2*pow(x,2)*y - 2*c0*i2*s0*u1*v2*pow(x,2)*y + 
   2*c0*c1*j2*k2*w0*pow(x,2)*y - 2*c0*c1*i2*p2*w0*pow(x,2)*y - 
   c1*k2*l2*r0*w0*pow(x,2)*y - c0*k2*l2*r1*w0*pow(x,2)*y + 
   e1*j2*l2*u0*w0*pow(x,2)*y - c1*j2*m2*u0*w0*pow(x,2)*y - 
   e1*i2*q2*u0*w0*pow(x,2)*y + c2*i2*r1*u0*w0*pow(x,2)*y - 
   l2*m2*r1*u0*w0*pow(x,2)*y + c1*i2*r2*u0*w0*pow(x,2)*y + 
   e0*j2*l2*u1*w0*pow(x,2)*y - c0*j2*m2*u1*w0*pow(x,2)*y - 
   e0*i2*q2*u1*w0*pow(x,2)*y + c2*i2*r0*u1*w0*pow(x,2)*y - 
   l2*m2*r0*u1*w0*pow(x,2)*y + c0*i2*r2*u1*w0*pow(x,2)*y + 
   c1*i2*r0*u2*w0*pow(x,2)*y + c0*i2*r1*u2*w0*pow(x,2)*y - 
   c1*j2*l2*v0*w0*pow(x,2)*y + c1*i2*q2*v0*w0*pow(x,2)*y + 
   pow(l2,2)*r1*v0*w0*pow(x,2)*y - c0*j2*l2*v1*w0*pow(x,2)*y + 
   c0*i2*q2*v1*w0*pow(x,2)*y + pow(l2,2)*r0*v1*w0*pow(x,2)*y + 
   pow(c0,2)*j2*k2*w1*pow(x,2)*y - pow(c0,2)*i2*p2*w1*pow(x,2)*y - 
   c0*k2*l2*r0*w1*pow(x,2)*y + e0*j2*l2*u0*w1*pow(x,2)*y - 
   c0*j2*m2*u0*w1*pow(x,2)*y - e0*i2*q2*u0*w1*pow(x,2)*y + 
   c2*i2*r0*u0*w1*pow(x,2)*y - l2*m2*r0*u0*w1*pow(x,2)*y + 
   c0*i2*r2*u0*w1*pow(x,2)*y + c0*i2*r0*u2*w1*pow(x,2)*y - 
   c0*j2*l2*v0*w1*pow(x,2)*y + c0*i2*q2*v0*w1*pow(x,2)*y + 
   pow(l2,2)*r0*v0*w1*pow(x,2)*y + c1*i2*r0*u0*w2*pow(x,2)*y + 
   c0*i2*r1*u0*w2*pow(x,2)*y + c0*i2*r0*u1*w2*pow(x,2)*y - 
   c1*k2*n0*r0*u0*pow(x,3)*y - c0*k2*n1*r0*u0*pow(x,3)*y - 
   c0*k2*n0*r1*u0*pow(x,3)*y - e1*j2*n0*pow(u0,2)*pow(x,3)*y - 
   e0*j2*n1*pow(u0,2)*pow(x,3)*y - f1*i2*r0*pow(u0,2)*pow(x,3)*y + 
   m2*n1*r0*pow(u0,2)*pow(x,3)*y - f0*i2*r1*pow(u0,2)*pow(x,3)*y + 
   m2*n0*r1*pow(u0,2)*pow(x,3)*y + e1*i2*s0*pow(u0,2)*pow(x,3)*y + 
   e0*i2*s1*pow(u0,2)*pow(x,3)*y - c0*k2*n0*r0*u1*pow(x,3)*y - 
   2*e0*j2*n0*u0*u1*pow(x,3)*y - 2*f0*i2*r0*u0*u1*pow(x,3)*y + 
   2*m2*n0*r0*u0*u1*pow(x,3)*y + 2*e0*i2*s0*u0*u1*pow(x,3)*y + 
   2*c1*j2*n0*u0*v0*pow(x,3)*y + 2*c0*j2*n1*u0*v0*pow(x,3)*y + 
   d1*i2*r0*u0*v0*pow(x,3)*y - l2*n1*r0*u0*v0*pow(x,3)*y + 
   d0*i2*r1*u0*v0*pow(x,3)*y - l2*n0*r1*u0*v0*pow(x,3)*y - 
   2*c1*i2*s0*u0*v0*pow(x,3)*y - 2*c0*i2*s1*u0*v0*pow(x,3)*y + 
   2*c0*j2*n0*u1*v0*pow(x,3)*y + d0*i2*r0*u1*v0*pow(x,3)*y - 
   l2*n0*r0*u1*v0*pow(x,3)*y - 2*c0*i2*s0*u1*v0*pow(x,3)*y + 
   2*c0*j2*n0*u0*v1*pow(x,3)*y + d0*i2*r0*u0*v1*pow(x,3)*y - 
   l2*n0*r0*u0*v1*pow(x,3)*y - 2*c0*i2*s0*u0*v1*pow(x,3)*y + 
   c1*i2*r0*u0*w0*pow(x,3)*y + c0*i2*r1*u0*w0*pow(x,3)*y + 
   c0*i2*r0*u1*w0*pow(x,3)*y + c0*i2*r0*u0*w1*pow(x,3)*y + 
   d1*e1*k2*l2*p2*pow(y,2) - c1*f1*k2*l2*p2*pow(y,2) - 
   c1*d1*k2*m2*p2*pow(y,2) + 2*c1*c2*k2*n1*p2*pow(y,2) + 
   pow(c1,2)*k2*n2*p2*pow(y,2) - d1*e1*pow(k2,2)*q2*pow(y,2) + 
   c1*f1*pow(k2,2)*q2*pow(y,2) + c2*d1*pow(k2,2)*r1*pow(y,2) + 
   c1*d2*pow(k2,2)*r1*pow(y,2) + c1*d1*pow(k2,2)*r2*pow(y,2) - 
   2*c1*c2*pow(k2,2)*s1*pow(y,2) - 
   pow(c1,2)*pow(k2,2)*s2*pow(y,2) - d1*e1*j2*l2*t2*pow(y,2) + 
   c1*f1*j2*l2*t2*pow(y,2) + c1*d1*j2*m2*t2*pow(y,2) - 
   2*c1*c2*j2*n1*t2*pow(y,2) - pow(c1,2)*j2*n2*t2*pow(y,2) + 
   d1*e1*i2*q2*t2*pow(y,2) - c1*f1*i2*q2*t2*pow(y,2) - 
   e1*l2*n1*q2*t2*pow(y,2) + c1*m2*n1*q2*t2*pow(y,2) - 
   c2*d1*i2*r1*t2*pow(y,2) - c1*d2*i2*r1*t2*pow(y,2) - 
   f1*pow(l2,2)*r1*t2*pow(y,2) + d1*l2*m2*r1*t2*pow(y,2) + 
   c2*l2*n1*r1*t2*pow(y,2) + c1*l2*n2*r1*t2*pow(y,2) - 
   c1*d1*i2*r2*t2*pow(y,2) + c1*l2*n1*r2*t2*pow(y,2) + 
   2*c1*c2*i2*s1*t2*pow(y,2) + e1*pow(l2,2)*s1*t2*pow(y,2) - 
   2*c1*l2*m2*s1*t2*pow(y,2) + pow(c1,2)*i2*s2*t2*pow(y,2) + 
   d2*e1*j2*k2*u1*pow(y,2) + d1*e2*j2*k2*u1*pow(y,2) - 
   c2*f1*j2*k2*u1*pow(y,2) - c1*f2*j2*k2*u1*pow(y,2) - 
   d2*e1*i2*p2*u1*pow(y,2) - d1*e2*i2*p2*u1*pow(y,2) + 
   c2*f1*i2*p2*u1*pow(y,2) + c1*f2*i2*p2*u1*pow(y,2) - 
   f1*l2*m2*p2*u1*pow(y,2) + d1*pow(m2,2)*p2*u1*pow(y,2) + 
   e2*l2*n1*p2*u1*pow(y,2) - c2*m2*n1*p2*u1*pow(y,2) + 
   e1*l2*n2*p2*u1*pow(y,2) - c1*m2*n2*p2*u1*pow(y,2) - 
   f1*k2*m2*q2*u1*pow(y,2) + e2*k2*n1*q2*u1*pow(y,2) + 
   e1*k2*n2*q2*u1*pow(y,2) + 2*f2*k2*l2*r1*u1*pow(y,2) - 
   d2*k2*m2*r1*u1*pow(y,2) - c2*k2*n2*r1*u1*pow(y,2) + 
   2*f1*k2*l2*r2*u1*pow(y,2) - d1*k2*m2*r2*u1*pow(y,2) - 
   c2*k2*n1*r2*u1*pow(y,2) - c1*k2*n2*r2*u1*pow(y,2) - 
   2*e2*k2*l2*s1*u1*pow(y,2) + 2*c2*k2*m2*s1*u1*pow(y,2) - 
   2*e1*k2*l2*s2*u1*pow(y,2) + 2*c1*k2*m2*s2*u1*pow(y,2) + 
   f2*j2*m2*pow(u1,2)*pow(y,2) - e2*j2*n2*pow(u1,2)*pow(y,2) - 
   f2*i2*r2*pow(u1,2)*pow(y,2) + m2*n2*r2*pow(u1,2)*pow(y,2) + 
   e2*i2*s2*pow(u1,2)*pow(y,2) - pow(m2,2)*s2*pow(u1,2)*pow(y,2) + 
   d1*e1*j2*k2*u2*pow(y,2) - c1*f1*j2*k2*u2*pow(y,2) - 
   d1*e1*i2*p2*u2*pow(y,2) + c1*f1*i2*p2*u2*pow(y,2) + 
   e1*l2*n1*p2*u2*pow(y,2) - c1*m2*n1*p2*u2*pow(y,2) + 
   e1*k2*n1*q2*u2*pow(y,2) + 2*f1*k2*l2*r1*u2*pow(y,2) - 
   d1*k2*m2*r1*u2*pow(y,2) - c2*k2*n1*r1*u2*pow(y,2) - 
   c1*k2*n2*r1*u2*pow(y,2) - c1*k2*n1*r2*u2*pow(y,2) - 
   2*e1*k2*l2*s1*u2*pow(y,2) + 2*c1*k2*m2*s1*u2*pow(y,2) + 
   2*f1*j2*m2*u1*u2*pow(y,2) - 2*e2*j2*n1*u1*u2*pow(y,2) - 
   2*e1*j2*n2*u1*u2*pow(y,2) - 2*f2*i2*r1*u1*u2*pow(y,2) + 
   2*m2*n2*r1*u1*u2*pow(y,2) - 2*f1*i2*r2*u1*u2*pow(y,2) + 
   2*m2*n1*r2*u1*u2*pow(y,2) + 2*e2*i2*s1*u1*u2*pow(y,2) - 
   2*pow(m2,2)*s1*u1*u2*pow(y,2) + 2*e1*i2*s2*u1*u2*pow(y,2) - 
   e1*j2*n1*pow(u2,2)*pow(y,2) - f1*i2*r1*pow(u2,2)*pow(y,2) + 
   m2*n1*r1*pow(u2,2)*pow(y,2) + e1*i2*s1*pow(u2,2)*pow(y,2) - 
   c2*d1*j2*k2*v1*pow(y,2) - c1*d2*j2*k2*v1*pow(y,2) + 
   c2*d1*i2*p2*v1*pow(y,2) + c1*d2*i2*p2*v1*pow(y,2) + 
   f1*pow(l2,2)*p2*v1*pow(y,2) - d1*l2*m2*p2*v1*pow(y,2) - 
   c2*l2*n1*p2*v1*pow(y,2) - c1*l2*n2*p2*v1*pow(y,2) - 
   f1*k2*l2*q2*v1*pow(y,2) + 2*d1*k2*m2*q2*v1*pow(y,2) - 
   c2*k2*n1*q2*v1*pow(y,2) - c1*k2*n2*q2*v1*pow(y,2) - 
   d2*k2*l2*r1*v1*pow(y,2) - d1*k2*l2*r2*v1*pow(y,2) + 
   2*c2*k2*l2*s1*v1*pow(y,2) + 2*c1*k2*l2*s2*v1*pow(y,2) - 
   f2*j2*l2*u1*v1*pow(y,2) - d2*j2*m2*u1*v1*pow(y,2) + 
   2*c2*j2*n2*u1*v1*pow(y,2) + f2*i2*q2*u1*v1*pow(y,2) - 
   m2*n2*q2*u1*v1*pow(y,2) + d2*i2*r2*u1*v1*pow(y,2) - 
   l2*n2*r2*u1*v1*pow(y,2) - 2*c2*i2*s2*u1*v1*pow(y,2) + 
   2*l2*m2*s2*u1*v1*pow(y,2) - f1*j2*l2*u2*v1*pow(y,2) - 
   d1*j2*m2*u2*v1*pow(y,2) + 2*c2*j2*n1*u2*v1*pow(y,2) + 
   2*c1*j2*n2*u2*v1*pow(y,2) + f1*i2*q2*u2*v1*pow(y,2) - 
   m2*n1*q2*u2*v1*pow(y,2) + d2*i2*r1*u2*v1*pow(y,2) - 
   l2*n2*r1*u2*v1*pow(y,2) + d1*i2*r2*u2*v1*pow(y,2) - 
   l2*n1*r2*u2*v1*pow(y,2) - 2*c2*i2*s1*u2*v1*pow(y,2) + 
   2*l2*m2*s1*u2*v1*pow(y,2) - 2*c1*i2*s2*u2*v1*pow(y,2) + 
   d2*j2*l2*pow(v1,2)*pow(y,2) - d2*i2*q2*pow(v1,2)*pow(y,2) + 
   l2*n2*q2*pow(v1,2)*pow(y,2) - pow(l2,2)*s2*pow(v1,2)*pow(y,2) - 
   c1*d1*j2*k2*v2*pow(y,2) + c1*d1*i2*p2*v2*pow(y,2) - 
   c1*l2*n1*p2*v2*pow(y,2) - c1*k2*n1*q2*v2*pow(y,2) - 
   d1*k2*l2*r1*v2*pow(y,2) + 2*c1*k2*l2*s1*v2*pow(y,2) - 
   f1*j2*l2*u1*v2*pow(y,2) - d1*j2*m2*u1*v2*pow(y,2) + 
   2*c2*j2*n1*u1*v2*pow(y,2) + 2*c1*j2*n2*u1*v2*pow(y,2) + 
   f1*i2*q2*u1*v2*pow(y,2) - m2*n1*q2*u1*v2*pow(y,2) + 
   d2*i2*r1*u1*v2*pow(y,2) - l2*n2*r1*u1*v2*pow(y,2) + 
   d1*i2*r2*u1*v2*pow(y,2) - l2*n1*r2*u1*v2*pow(y,2) - 
   2*c2*i2*s1*u1*v2*pow(y,2) + 2*l2*m2*s1*u1*v2*pow(y,2) - 
   2*c1*i2*s2*u1*v2*pow(y,2) + 2*c1*j2*n1*u2*v2*pow(y,2) + 
   d1*i2*r1*u2*v2*pow(y,2) - l2*n1*r1*u2*v2*pow(y,2) - 
   2*c1*i2*s1*u2*v2*pow(y,2) + 2*d1*j2*l2*v1*v2*pow(y,2) - 
   2*d1*i2*q2*v1*v2*pow(y,2) + 2*l2*n1*q2*v1*v2*pow(y,2) - 
   2*pow(l2,2)*s1*v1*v2*pow(y,2) + 2*c1*c2*j2*k2*w1*pow(y,2) - 
   2*c1*c2*i2*p2*w1*pow(y,2) - e1*pow(l2,2)*p2*w1*pow(y,2) + 
   2*c1*l2*m2*p2*w1*pow(y,2) + e1*k2*l2*q2*w1*pow(y,2) - 
   c1*k2*m2*q2*w1*pow(y,2) - c2*k2*l2*r1*w1*pow(y,2) - 
   c1*k2*l2*r2*w1*pow(y,2) + e2*j2*l2*u1*w1*pow(y,2) - 
   c2*j2*m2*u1*w1*pow(y,2) - e2*i2*q2*u1*w1*pow(y,2) + 
   pow(m2,2)*q2*u1*w1*pow(y,2) + c2*i2*r2*u1*w1*pow(y,2) - 
   l2*m2*r2*u1*w1*pow(y,2) + e1*j2*l2*u2*w1*pow(y,2) - 
   c1*j2*m2*u2*w1*pow(y,2) - e1*i2*q2*u2*w1*pow(y,2) + 
   c2*i2*r1*u2*w1*pow(y,2) - l2*m2*r1*u2*w1*pow(y,2) + 
   c1*i2*r2*u2*w1*pow(y,2) - c2*j2*l2*v1*w1*pow(y,2) + 
   c2*i2*q2*v1*w1*pow(y,2) - l2*m2*q2*v1*w1*pow(y,2) + 
   pow(l2,2)*r2*v1*w1*pow(y,2) - c1*j2*l2*v2*w1*pow(y,2) + 
   c1*i2*q2*v2*w1*pow(y,2) + pow(l2,2)*r1*v2*w1*pow(y,2) + 
   pow(c1,2)*j2*k2*w2*pow(y,2) - pow(c1,2)*i2*p2*w2*pow(y,2) - 
   c1*k2*l2*r1*w2*pow(y,2) + e1*j2*l2*u1*w2*pow(y,2) - 
   c1*j2*m2*u1*w2*pow(y,2) - e1*i2*q2*u1*w2*pow(y,2) + 
   c2*i2*r1*u1*w2*pow(y,2) - l2*m2*r1*u1*w2*pow(y,2) + 
   c1*i2*r2*u1*w2*pow(y,2) + c1*i2*r1*u2*w2*pow(y,2) - 
   c1*j2*l2*v1*w2*pow(y,2) + c1*i2*q2*v1*w2*pow(y,2) + 
   pow(l2,2)*r1*v1*w2*pow(y,2) + pow(c1,2)*k2*n0*p2*x*pow(y,2) + 
   2*c0*c1*k2*n1*p2*x*pow(y,2) + c1*d1*pow(k2,2)*r0*x*pow(y,2) + 
   c1*d0*pow(k2,2)*r1*x*pow(y,2) + c0*d1*pow(k2,2)*r1*x*pow(y,2) - 
   pow(c1,2)*pow(k2,2)*s0*x*pow(y,2) - 
   2*c0*c1*pow(k2,2)*s1*x*pow(y,2) - pow(c1,2)*j2*n0*t2*x*pow(y,2) - 
   2*c0*c1*j2*n1*t2*x*pow(y,2) - c1*d1*i2*r0*t2*x*pow(y,2) + 
   c1*l2*n1*r0*t2*x*pow(y,2) - c1*d0*i2*r1*t2*x*pow(y,2) - 
   c0*d1*i2*r1*t2*x*pow(y,2) + c1*l2*n0*r1*t2*x*pow(y,2) + 
   c0*l2*n1*r1*t2*x*pow(y,2) + pow(c1,2)*i2*s0*t2*x*pow(y,2) + 
   2*c0*c1*i2*s1*t2*x*pow(y,2) + d1*e1*j2*k2*u0*x*pow(y,2) - 
   c1*f1*j2*k2*u0*x*pow(y,2) - d1*e1*i2*p2*u0*x*pow(y,2) + 
   c1*f1*i2*p2*u0*x*pow(y,2) + e1*l2*n1*p2*u0*x*pow(y,2) - 
   c1*m2*n1*p2*u0*x*pow(y,2) + e1*k2*n1*q2*u0*x*pow(y,2) + 
   2*f1*k2*l2*r1*u0*x*pow(y,2) - d1*k2*m2*r1*u0*x*pow(y,2) - 
   c2*k2*n1*r1*u0*x*pow(y,2) - c1*k2*n2*r1*u0*x*pow(y,2) - 
   c1*k2*n1*r2*u0*x*pow(y,2) - 2*e1*k2*l2*s1*u0*x*pow(y,2) + 
   2*c1*k2*m2*s1*u0*x*pow(y,2) + d1*e0*j2*k2*u1*x*pow(y,2) + 
   d0*e1*j2*k2*u1*x*pow(y,2) - c1*f0*j2*k2*u1*x*pow(y,2) - 
   c0*f1*j2*k2*u1*x*pow(y,2) - d1*e0*i2*p2*u1*x*pow(y,2) - 
   d0*e1*i2*p2*u1*x*pow(y,2) + c1*f0*i2*p2*u1*x*pow(y,2) + 
   c0*f1*i2*p2*u1*x*pow(y,2) + e1*l2*n0*p2*u1*x*pow(y,2) - 
   c1*m2*n0*p2*u1*x*pow(y,2) + e0*l2*n1*p2*u1*x*pow(y,2) - 
   c0*m2*n1*p2*u1*x*pow(y,2) + e1*k2*n0*q2*u1*x*pow(y,2) + 
   e0*k2*n1*q2*u1*x*pow(y,2) + 2*f1*k2*l2*r0*u1*x*pow(y,2) - 
   d1*k2*m2*r0*u1*x*pow(y,2) - c2*k2*n1*r0*u1*x*pow(y,2) - 
   c1*k2*n2*r0*u1*x*pow(y,2) + 2*f0*k2*l2*r1*u1*x*pow(y,2) - 
   d0*k2*m2*r1*u1*x*pow(y,2) - c2*k2*n0*r1*u1*x*pow(y,2) - 
   c0*k2*n2*r1*u1*x*pow(y,2) - c1*k2*n0*r2*u1*x*pow(y,2) - 
   c0*k2*n1*r2*u1*x*pow(y,2) - 2*e1*k2*l2*s0*u1*x*pow(y,2) + 
   2*c1*k2*m2*s0*u1*x*pow(y,2) - 2*e0*k2*l2*s1*u1*x*pow(y,2) + 
   2*c0*k2*m2*s1*u1*x*pow(y,2) + 2*f1*j2*m2*u0*u1*x*pow(y,2) - 
   2*e2*j2*n1*u0*u1*x*pow(y,2) - 2*e1*j2*n2*u0*u1*x*pow(y,2) - 
   2*f2*i2*r1*u0*u1*x*pow(y,2) + 2*m2*n2*r1*u0*u1*x*pow(y,2) - 
   2*f1*i2*r2*u0*u1*x*pow(y,2) + 2*m2*n1*r2*u0*u1*x*pow(y,2) + 
   2*e2*i2*s1*u0*u1*x*pow(y,2) - 2*pow(m2,2)*s1*u0*u1*x*pow(y,2) + 
   2*e1*i2*s2*u0*u1*x*pow(y,2) + f0*j2*m2*pow(u1,2)*x*pow(y,2) - 
   e2*j2*n0*pow(u1,2)*x*pow(y,2) - e0*j2*n2*pow(u1,2)*x*pow(y,2) - 
   f2*i2*r0*pow(u1,2)*x*pow(y,2) + m2*n2*r0*pow(u1,2)*x*pow(y,2) - 
   f0*i2*r2*pow(u1,2)*x*pow(y,2) + m2*n0*r2*pow(u1,2)*x*pow(y,2) + 
   e2*i2*s0*pow(u1,2)*x*pow(y,2) - 
   pow(m2,2)*s0*pow(u1,2)*x*pow(y,2) + 
   e0*i2*s2*pow(u1,2)*x*pow(y,2) - c1*k2*n1*r0*u2*x*pow(y,2) - 
   c1*k2*n0*r1*u2*x*pow(y,2) - c0*k2*n1*r1*u2*x*pow(y,2) - 
   2*e1*j2*n1*u0*u2*x*pow(y,2) - 2*f1*i2*r1*u0*u2*x*pow(y,2) + 
   2*m2*n1*r1*u0*u2*x*pow(y,2) + 2*e1*i2*s1*u0*u2*x*pow(y,2) - 
   2*e1*j2*n0*u1*u2*x*pow(y,2) - 2*e0*j2*n1*u1*u2*x*pow(y,2) - 
   2*f1*i2*r0*u1*u2*x*pow(y,2) + 2*m2*n1*r0*u1*u2*x*pow(y,2) - 
   2*f0*i2*r1*u1*u2*x*pow(y,2) + 2*m2*n0*r1*u1*u2*x*pow(y,2) + 
   2*e1*i2*s0*u1*u2*x*pow(y,2) + 2*e0*i2*s1*u1*u2*x*pow(y,2) - 
   c1*d1*j2*k2*v0*x*pow(y,2) + c1*d1*i2*p2*v0*x*pow(y,2) - 
   c1*l2*n1*p2*v0*x*pow(y,2) - c1*k2*n1*q2*v0*x*pow(y,2) - 
   d1*k2*l2*r1*v0*x*pow(y,2) + 2*c1*k2*l2*s1*v0*x*pow(y,2) - 
   f1*j2*l2*u1*v0*x*pow(y,2) - d1*j2*m2*u1*v0*x*pow(y,2) + 
   2*c2*j2*n1*u1*v0*x*pow(y,2) + 2*c1*j2*n2*u1*v0*x*pow(y,2) + 
   f1*i2*q2*u1*v0*x*pow(y,2) - m2*n1*q2*u1*v0*x*pow(y,2) + 
   d2*i2*r1*u1*v0*x*pow(y,2) - l2*n2*r1*u1*v0*x*pow(y,2) + 
   d1*i2*r2*u1*v0*x*pow(y,2) - l2*n1*r2*u1*v0*x*pow(y,2) - 
   2*c2*i2*s1*u1*v0*x*pow(y,2) + 2*l2*m2*s1*u1*v0*x*pow(y,2) - 
   2*c1*i2*s2*u1*v0*x*pow(y,2) + 2*c1*j2*n1*u2*v0*x*pow(y,2) + 
   d1*i2*r1*u2*v0*x*pow(y,2) - l2*n1*r1*u2*v0*x*pow(y,2) - 
   2*c1*i2*s1*u2*v0*x*pow(y,2) - c1*d0*j2*k2*v1*x*pow(y,2) - 
   c0*d1*j2*k2*v1*x*pow(y,2) + c1*d0*i2*p2*v1*x*pow(y,2) + 
   c0*d1*i2*p2*v1*x*pow(y,2) - c1*l2*n0*p2*v1*x*pow(y,2) - 
   c0*l2*n1*p2*v1*x*pow(y,2) - c1*k2*n0*q2*v1*x*pow(y,2) - 
   c0*k2*n1*q2*v1*x*pow(y,2) - d1*k2*l2*r0*v1*x*pow(y,2) - 
   d0*k2*l2*r1*v1*x*pow(y,2) + 2*c1*k2*l2*s0*v1*x*pow(y,2) + 
   2*c0*k2*l2*s1*v1*x*pow(y,2) - f1*j2*l2*u0*v1*x*pow(y,2) - 
   d1*j2*m2*u0*v1*x*pow(y,2) + 2*c2*j2*n1*u0*v1*x*pow(y,2) + 
   2*c1*j2*n2*u0*v1*x*pow(y,2) + f1*i2*q2*u0*v1*x*pow(y,2) - 
   m2*n1*q2*u0*v1*x*pow(y,2) + d2*i2*r1*u0*v1*x*pow(y,2) - 
   l2*n2*r1*u0*v1*x*pow(y,2) + d1*i2*r2*u0*v1*x*pow(y,2) - 
   l2*n1*r2*u0*v1*x*pow(y,2) - 2*c2*i2*s1*u0*v1*x*pow(y,2) + 
   2*l2*m2*s1*u0*v1*x*pow(y,2) - 2*c1*i2*s2*u0*v1*x*pow(y,2) - 
   f0*j2*l2*u1*v1*x*pow(y,2) - d0*j2*m2*u1*v1*x*pow(y,2) + 
   2*c2*j2*n0*u1*v1*x*pow(y,2) + 2*c0*j2*n2*u1*v1*x*pow(y,2) + 
   f0*i2*q2*u1*v1*x*pow(y,2) - m2*n0*q2*u1*v1*x*pow(y,2) + 
   d2*i2*r0*u1*v1*x*pow(y,2) - l2*n2*r0*u1*v1*x*pow(y,2) + 
   d0*i2*r2*u1*v1*x*pow(y,2) - l2*n0*r2*u1*v1*x*pow(y,2) - 
   2*c2*i2*s0*u1*v1*x*pow(y,2) + 2*l2*m2*s0*u1*v1*x*pow(y,2) - 
   2*c0*i2*s2*u1*v1*x*pow(y,2) + 2*c1*j2*n0*u2*v1*x*pow(y,2) + 
   2*c0*j2*n1*u2*v1*x*pow(y,2) + d1*i2*r0*u2*v1*x*pow(y,2) - 
   l2*n1*r0*u2*v1*x*pow(y,2) + d0*i2*r1*u2*v1*x*pow(y,2) - 
   l2*n0*r1*u2*v1*x*pow(y,2) - 2*c1*i2*s0*u2*v1*x*pow(y,2) - 
   2*c0*i2*s1*u2*v1*x*pow(y,2) + 2*d1*j2*l2*v0*v1*x*pow(y,2) - 
   2*d1*i2*q2*v0*v1*x*pow(y,2) + 2*l2*n1*q2*v0*v1*x*pow(y,2) - 
   2*pow(l2,2)*s1*v0*v1*x*pow(y,2) + d0*j2*l2*pow(v1,2)*x*pow(y,2) - 
   d0*i2*q2*pow(v1,2)*x*pow(y,2) + l2*n0*q2*pow(v1,2)*x*pow(y,2) - 
   pow(l2,2)*s0*pow(v1,2)*x*pow(y,2) + 2*c1*j2*n1*u0*v2*x*pow(y,2) + 
   d1*i2*r1*u0*v2*x*pow(y,2) - l2*n1*r1*u0*v2*x*pow(y,2) - 
   2*c1*i2*s1*u0*v2*x*pow(y,2) + 2*c1*j2*n0*u1*v2*x*pow(y,2) + 
   2*c0*j2*n1*u1*v2*x*pow(y,2) + d1*i2*r0*u1*v2*x*pow(y,2) - 
   l2*n1*r0*u1*v2*x*pow(y,2) + d0*i2*r1*u1*v2*x*pow(y,2) - 
   l2*n0*r1*u1*v2*x*pow(y,2) - 2*c1*i2*s0*u1*v2*x*pow(y,2) - 
   2*c0*i2*s1*u1*v2*x*pow(y,2) + pow(c1,2)*j2*k2*w0*x*pow(y,2) - 
   pow(c1,2)*i2*p2*w0*x*pow(y,2) - c1*k2*l2*r1*w0*x*pow(y,2) + 
   e1*j2*l2*u1*w0*x*pow(y,2) - c1*j2*m2*u1*w0*x*pow(y,2) - 
   e1*i2*q2*u1*w0*x*pow(y,2) + c2*i2*r1*u1*w0*x*pow(y,2) - 
   l2*m2*r1*u1*w0*x*pow(y,2) + c1*i2*r2*u1*w0*x*pow(y,2) + 
   c1*i2*r1*u2*w0*x*pow(y,2) - c1*j2*l2*v1*w0*x*pow(y,2) + 
   c1*i2*q2*v1*w0*x*pow(y,2) + pow(l2,2)*r1*v1*w0*x*pow(y,2) + 
   2*c0*c1*j2*k2*w1*x*pow(y,2) - 2*c0*c1*i2*p2*w1*x*pow(y,2) - 
   c1*k2*l2*r0*w1*x*pow(y,2) - c0*k2*l2*r1*w1*x*pow(y,2) + 
   e1*j2*l2*u0*w1*x*pow(y,2) - c1*j2*m2*u0*w1*x*pow(y,2) - 
   e1*i2*q2*u0*w1*x*pow(y,2) + c2*i2*r1*u0*w1*x*pow(y,2) - 
   l2*m2*r1*u0*w1*x*pow(y,2) + c1*i2*r2*u0*w1*x*pow(y,2) + 
   e0*j2*l2*u1*w1*x*pow(y,2) - c0*j2*m2*u1*w1*x*pow(y,2) - 
   e0*i2*q2*u1*w1*x*pow(y,2) + c2*i2*r0*u1*w1*x*pow(y,2) - 
   l2*m2*r0*u1*w1*x*pow(y,2) + c0*i2*r2*u1*w1*x*pow(y,2) + 
   c1*i2*r0*u2*w1*x*pow(y,2) + c0*i2*r1*u2*w1*x*pow(y,2) - 
   c1*j2*l2*v0*w1*x*pow(y,2) + c1*i2*q2*v0*w1*x*pow(y,2) + 
   pow(l2,2)*r1*v0*w1*x*pow(y,2) - c0*j2*l2*v1*w1*x*pow(y,2) + 
   c0*i2*q2*v1*w1*x*pow(y,2) + pow(l2,2)*r0*v1*w1*x*pow(y,2) + 
   c1*i2*r1*u0*w2*x*pow(y,2) + c1*i2*r0*u1*w2*x*pow(y,2) + 
   c0*i2*r1*u1*w2*x*pow(y,2) - c1*k2*n1*r0*u0*pow(x,2)*pow(y,2) - 
   c1*k2*n0*r1*u0*pow(x,2)*pow(y,2) - 
   c0*k2*n1*r1*u0*pow(x,2)*pow(y,2) - 
   e1*j2*n1*pow(u0,2)*pow(x,2)*pow(y,2) - 
   f1*i2*r1*pow(u0,2)*pow(x,2)*pow(y,2) + 
   m2*n1*r1*pow(u0,2)*pow(x,2)*pow(y,2) + 
   e1*i2*s1*pow(u0,2)*pow(x,2)*pow(y,2) - 
   c1*k2*n0*r0*u1*pow(x,2)*pow(y,2) - 
   c0*k2*n1*r0*u1*pow(x,2)*pow(y,2) - 
   c0*k2*n0*r1*u1*pow(x,2)*pow(y,2) - 
   2*e1*j2*n0*u0*u1*pow(x,2)*pow(y,2) - 
   2*e0*j2*n1*u0*u1*pow(x,2)*pow(y,2) - 
   2*f1*i2*r0*u0*u1*pow(x,2)*pow(y,2) + 
   2*m2*n1*r0*u0*u1*pow(x,2)*pow(y,2) - 
   2*f0*i2*r1*u0*u1*pow(x,2)*pow(y,2) + 
   2*m2*n0*r1*u0*u1*pow(x,2)*pow(y,2) + 
   2*e1*i2*s0*u0*u1*pow(x,2)*pow(y,2) + 
   2*e0*i2*s1*u0*u1*pow(x,2)*pow(y,2) - 
   e0*j2*n0*pow(u1,2)*pow(x,2)*pow(y,2) - 
   f0*i2*r0*pow(u1,2)*pow(x,2)*pow(y,2) + 
   m2*n0*r0*pow(u1,2)*pow(x,2)*pow(y,2) + 
   e0*i2*s0*pow(u1,2)*pow(x,2)*pow(y,2) + 
   2*c1*j2*n1*u0*v0*pow(x,2)*pow(y,2) + 
   d1*i2*r1*u0*v0*pow(x,2)*pow(y,2) - 
   l2*n1*r1*u0*v0*pow(x,2)*pow(y,2) - 
   2*c1*i2*s1*u0*v0*pow(x,2)*pow(y,2) + 
   2*c1*j2*n0*u1*v0*pow(x,2)*pow(y,2) + 
   2*c0*j2*n1*u1*v0*pow(x,2)*pow(y,2) + 
   d1*i2*r0*u1*v0*pow(x,2)*pow(y,2) - 
   l2*n1*r0*u1*v0*pow(x,2)*pow(y,2) + 
   d0*i2*r1*u1*v0*pow(x,2)*pow(y,2) - 
   l2*n0*r1*u1*v0*pow(x,2)*pow(y,2) - 
   2*c1*i2*s0*u1*v0*pow(x,2)*pow(y,2) - 
   2*c0*i2*s1*u1*v0*pow(x,2)*pow(y,2) + 
   2*c1*j2*n0*u0*v1*pow(x,2)*pow(y,2) + 
   2*c0*j2*n1*u0*v1*pow(x,2)*pow(y,2) + 
   d1*i2*r0*u0*v1*pow(x,2)*pow(y,2) - 
   l2*n1*r0*u0*v1*pow(x,2)*pow(y,2) + 
   d0*i2*r1*u0*v1*pow(x,2)*pow(y,2) - 
   l2*n0*r1*u0*v1*pow(x,2)*pow(y,2) - 
   2*c1*i2*s0*u0*v1*pow(x,2)*pow(y,2) - 
   2*c0*i2*s1*u0*v1*pow(x,2)*pow(y,2) + 
   2*c0*j2*n0*u1*v1*pow(x,2)*pow(y,2) + 
   d0*i2*r0*u1*v1*pow(x,2)*pow(y,2) - 
   l2*n0*r0*u1*v1*pow(x,2)*pow(y,2) - 
   2*c0*i2*s0*u1*v1*pow(x,2)*pow(y,2) + 
   c1*i2*r1*u0*w0*pow(x,2)*pow(y,2) + 
   c1*i2*r0*u1*w0*pow(x,2)*pow(y,2) + 
   c0*i2*r1*u1*w0*pow(x,2)*pow(y,2) + 
   c1*i2*r0*u0*w1*pow(x,2)*pow(y,2) + 
   c0*i2*r1*u0*w1*pow(x,2)*pow(y,2) + 
   c0*i2*r0*u1*w1*pow(x,2)*pow(y,2) + pow(c1,2)*k2*n1*p2*pow(y,3) + 
   c1*d1*pow(k2,2)*r1*pow(y,3) - pow(c1,2)*pow(k2,2)*s1*pow(y,3) - 
   pow(c1,2)*j2*n1*t2*pow(y,3) - c1*d1*i2*r1*t2*pow(y,3) + 
   c1*l2*n1*r1*t2*pow(y,3) + pow(c1,2)*i2*s1*t2*pow(y,3) + 
   d1*e1*j2*k2*u1*pow(y,3) - c1*f1*j2*k2*u1*pow(y,3) - 
   d1*e1*i2*p2*u1*pow(y,3) + c1*f1*i2*p2*u1*pow(y,3) + 
   e1*l2*n1*p2*u1*pow(y,3) - c1*m2*n1*p2*u1*pow(y,3) + 
   e1*k2*n1*q2*u1*pow(y,3) + 2*f1*k2*l2*r1*u1*pow(y,3) - 
   d1*k2*m2*r1*u1*pow(y,3) - c2*k2*n1*r1*u1*pow(y,3) - 
   c1*k2*n2*r1*u1*pow(y,3) - c1*k2*n1*r2*u1*pow(y,3) - 
   2*e1*k2*l2*s1*u1*pow(y,3) + 2*c1*k2*m2*s1*u1*pow(y,3) + 
   f1*j2*m2*pow(u1,2)*pow(y,3) - e2*j2*n1*pow(u1,2)*pow(y,3) - 
   e1*j2*n2*pow(u1,2)*pow(y,3) - f2*i2*r1*pow(u1,2)*pow(y,3) + 
   m2*n2*r1*pow(u1,2)*pow(y,3) - f1*i2*r2*pow(u1,2)*pow(y,3) + 
   m2*n1*r2*pow(u1,2)*pow(y,3) + e2*i2*s1*pow(u1,2)*pow(y,3) - 
   pow(m2,2)*s1*pow(u1,2)*pow(y,3) + e1*i2*s2*pow(u1,2)*pow(y,3) - 
   c1*k2*n1*r1*u2*pow(y,3) - 2*e1*j2*n1*u1*u2*pow(y,3) - 
   2*f1*i2*r1*u1*u2*pow(y,3) + 2*m2*n1*r1*u1*u2*pow(y,3) + 
   2*e1*i2*s1*u1*u2*pow(y,3) - c1*d1*j2*k2*v1*pow(y,3) + 
   c1*d1*i2*p2*v1*pow(y,3) - c1*l2*n1*p2*v1*pow(y,3) - 
   c1*k2*n1*q2*v1*pow(y,3) - d1*k2*l2*r1*v1*pow(y,3) + 
   2*c1*k2*l2*s1*v1*pow(y,3) - f1*j2*l2*u1*v1*pow(y,3) - 
   d1*j2*m2*u1*v1*pow(y,3) + 2*c2*j2*n1*u1*v1*pow(y,3) + 
   2*c1*j2*n2*u1*v1*pow(y,3) + f1*i2*q2*u1*v1*pow(y,3) - 
   m2*n1*q2*u1*v1*pow(y,3) + d2*i2*r1*u1*v1*pow(y,3) - 
   l2*n2*r1*u1*v1*pow(y,3) + d1*i2*r2*u1*v1*pow(y,3) - 
   l2*n1*r2*u1*v1*pow(y,3) - 2*c2*i2*s1*u1*v1*pow(y,3) + 
   2*l2*m2*s1*u1*v1*pow(y,3) - 2*c1*i2*s2*u1*v1*pow(y,3) + 
   2*c1*j2*n1*u2*v1*pow(y,3) + d1*i2*r1*u2*v1*pow(y,3) - 
   l2*n1*r1*u2*v1*pow(y,3) - 2*c1*i2*s1*u2*v1*pow(y,3) + 
   d1*j2*l2*pow(v1,2)*pow(y,3) - d1*i2*q2*pow(v1,2)*pow(y,3) + 
   l2*n1*q2*pow(v1,2)*pow(y,3) - pow(l2,2)*s1*pow(v1,2)*pow(y,3) + 
   2*c1*j2*n1*u1*v2*pow(y,3) + d1*i2*r1*u1*v2*pow(y,3) - 
   l2*n1*r1*u1*v2*pow(y,3) - 2*c1*i2*s1*u1*v2*pow(y,3) + 
   pow(c1,2)*j2*k2*w1*pow(y,3) - pow(c1,2)*i2*p2*w1*pow(y,3) - 
   c1*k2*l2*r1*w1*pow(y,3) + e1*j2*l2*u1*w1*pow(y,3) - 
   c1*j2*m2*u1*w1*pow(y,3) - e1*i2*q2*u1*w1*pow(y,3) + 
   c2*i2*r1*u1*w1*pow(y,3) - l2*m2*r1*u1*w1*pow(y,3) + 
   c1*i2*r2*u1*w1*pow(y,3) + c1*i2*r1*u2*w1*pow(y,3) - 
   c1*j2*l2*v1*w1*pow(y,3) + c1*i2*q2*v1*w1*pow(y,3) + 
   pow(l2,2)*r1*v1*w1*pow(y,3) + c1*i2*r1*u1*w2*pow(y,3) - 
   c1*k2*n1*r1*u0*x*pow(y,3) - c1*k2*n1*r0*u1*x*pow(y,3) - 
   c1*k2*n0*r1*u1*x*pow(y,3) - c0*k2*n1*r1*u1*x*pow(y,3) - 
   2*e1*j2*n1*u0*u1*x*pow(y,3) - 2*f1*i2*r1*u0*u1*x*pow(y,3) + 
   2*m2*n1*r1*u0*u1*x*pow(y,3) + 2*e1*i2*s1*u0*u1*x*pow(y,3) - 
   e1*j2*n0*pow(u1,2)*x*pow(y,3) - e0*j2*n1*pow(u1,2)*x*pow(y,3) - 
   f1*i2*r0*pow(u1,2)*x*pow(y,3) + m2*n1*r0*pow(u1,2)*x*pow(y,3) - 
   f0*i2*r1*pow(u1,2)*x*pow(y,3) + m2*n0*r1*pow(u1,2)*x*pow(y,3) + 
   e1*i2*s0*pow(u1,2)*x*pow(y,3) + e0*i2*s1*pow(u1,2)*x*pow(y,3) + 
   2*c1*j2*n1*u1*v0*x*pow(y,3) + d1*i2*r1*u1*v0*x*pow(y,3) - 
   l2*n1*r1*u1*v0*x*pow(y,3) - 2*c1*i2*s1*u1*v0*x*pow(y,3) + 
   2*c1*j2*n1*u0*v1*x*pow(y,3) + d1*i2*r1*u0*v1*x*pow(y,3) - 
   l2*n1*r1*u0*v1*x*pow(y,3) - 2*c1*i2*s1*u0*v1*x*pow(y,3) + 
   2*c1*j2*n0*u1*v1*x*pow(y,3) + 2*c0*j2*n1*u1*v1*x*pow(y,3) + 
   d1*i2*r0*u1*v1*x*pow(y,3) - l2*n1*r0*u1*v1*x*pow(y,3) + 
   d0*i2*r1*u1*v1*x*pow(y,3) - l2*n0*r1*u1*v1*x*pow(y,3) - 
   2*c1*i2*s0*u1*v1*x*pow(y,3) - 2*c0*i2*s1*u1*v1*x*pow(y,3) + 
   c1*i2*r1*u1*w0*x*pow(y,3) + c1*i2*r1*u0*w1*x*pow(y,3) + 
   c1*i2*r0*u1*w1*x*pow(y,3) + c0*i2*r1*u1*w1*x*pow(y,3) - 
   c1*k2*n1*r1*u1*pow(y,4) - e1*j2*n1*pow(u1,2)*pow(y,4) - 
   f1*i2*r1*pow(u1,2)*pow(y,4) + m2*n1*r1*pow(u1,2)*pow(y,4) + 
   e1*i2*s1*pow(u1,2)*pow(y,4) + 2*c1*j2*n1*u1*v1*pow(y,4) + 
   d1*i2*r1*u1*v1*pow(y,4) - l2*n1*r1*u1*v1*pow(y,4) - 
   2*c1*i2*s1*u1*v1*pow(y,4) + c1*i2*r1*u1*w1*pow(y,4) + 
   f2*k2*m2*p2*x*z0 - e2*k2*n2*p2*x*z0 - f2*pow(k2,2)*r2*x*z0 + 
   e2*pow(k2,2)*s2*x*z0 - f2*j2*m2*t2*x*z0 + e2*j2*n2*t2*x*z0 + 
   f2*i2*r2*t2*x*z0 - m2*n2*r2*t2*x*z0 - e2*i2*s2*t2*x*z0 + 
   pow(m2,2)*s2*t2*x*z0 + f2*j2*k2*v2*x*z0 - f2*i2*p2*v2*x*z0 + 
   m2*n2*p2*v2*x*z0 + k2*n2*r2*v2*x*z0 - 2*k2*m2*s2*v2*x*z0 - 
   j2*n2*pow(v2,2)*x*z0 + i2*s2*pow(v2,2)*x*z0 - e2*j2*k2*w2*x*z0 + 
   e2*i2*p2*w2*x*z0 - pow(m2,2)*p2*w2*x*z0 + k2*m2*r2*w2*x*z0 + 
   j2*m2*v2*w2*x*z0 - i2*r2*v2*w2*x*z0 + f0*k2*m2*p2*pow(x,2)*z0 - 
   e2*k2*n0*p2*pow(x,2)*z0 - e0*k2*n2*p2*pow(x,2)*z0 - 
   f2*pow(k2,2)*r0*pow(x,2)*z0 - f0*pow(k2,2)*r2*pow(x,2)*z0 + 
   e2*pow(k2,2)*s0*pow(x,2)*z0 + e0*pow(k2,2)*s2*pow(x,2)*z0 - 
   f0*j2*m2*t2*pow(x,2)*z0 + e2*j2*n0*t2*pow(x,2)*z0 + 
   e0*j2*n2*t2*pow(x,2)*z0 + f2*i2*r0*t2*pow(x,2)*z0 - 
   m2*n2*r0*t2*pow(x,2)*z0 + f0*i2*r2*t2*pow(x,2)*z0 - 
   m2*n0*r2*t2*pow(x,2)*z0 - e2*i2*s0*t2*pow(x,2)*z0 + 
   pow(m2,2)*s0*t2*pow(x,2)*z0 - e0*i2*s2*t2*pow(x,2)*z0 + 
   f2*j2*k2*v0*pow(x,2)*z0 - f2*i2*p2*v0*pow(x,2)*z0 + 
   m2*n2*p2*v0*pow(x,2)*z0 + k2*n2*r2*v0*pow(x,2)*z0 - 
   2*k2*m2*s2*v0*pow(x,2)*z0 + f0*j2*k2*v2*pow(x,2)*z0 - 
   f0*i2*p2*v2*pow(x,2)*z0 + m2*n0*p2*v2*pow(x,2)*z0 + 
   k2*n2*r0*v2*pow(x,2)*z0 + k2*n0*r2*v2*pow(x,2)*z0 - 
   2*k2*m2*s0*v2*pow(x,2)*z0 - 2*j2*n2*v0*v2*pow(x,2)*z0 + 
   2*i2*s2*v0*v2*pow(x,2)*z0 - j2*n0*pow(v2,2)*pow(x,2)*z0 + 
   i2*s0*pow(v2,2)*pow(x,2)*z0 - e2*j2*k2*w0*pow(x,2)*z0 + 
   e2*i2*p2*w0*pow(x,2)*z0 - pow(m2,2)*p2*w0*pow(x,2)*z0 + 
   k2*m2*r2*w0*pow(x,2)*z0 + j2*m2*v2*w0*pow(x,2)*z0 - 
   i2*r2*v2*w0*pow(x,2)*z0 - e0*j2*k2*w2*pow(x,2)*z0 + 
   e0*i2*p2*w2*pow(x,2)*z0 + k2*m2*r0*w2*pow(x,2)*z0 + 
   j2*m2*v0*w2*pow(x,2)*z0 - i2*r2*v0*w2*pow(x,2)*z0 - 
   i2*r0*v2*w2*pow(x,2)*z0 - e0*k2*n0*p2*pow(x,3)*z0 - 
   f0*pow(k2,2)*r0*pow(x,3)*z0 + e0*pow(k2,2)*s0*pow(x,3)*z0 + 
   e0*j2*n0*t2*pow(x,3)*z0 + f0*i2*r0*t2*pow(x,3)*z0 - 
   m2*n0*r0*t2*pow(x,3)*z0 - e0*i2*s0*t2*pow(x,3)*z0 + 
   f0*j2*k2*v0*pow(x,3)*z0 - f0*i2*p2*v0*pow(x,3)*z0 + 
   m2*n0*p2*v0*pow(x,3)*z0 + k2*n2*r0*v0*pow(x,3)*z0 + 
   k2*n0*r2*v0*pow(x,3)*z0 - 2*k2*m2*s0*v0*pow(x,3)*z0 - 
   j2*n2*pow(v0,2)*pow(x,3)*z0 + i2*s2*pow(v0,2)*pow(x,3)*z0 + 
   k2*n0*r0*v2*pow(x,3)*z0 - 2*j2*n0*v0*v2*pow(x,3)*z0 + 
   2*i2*s0*v0*v2*pow(x,3)*z0 - e0*j2*k2*w0*pow(x,3)*z0 + 
   e0*i2*p2*w0*pow(x,3)*z0 + k2*m2*r0*w0*pow(x,3)*z0 + 
   j2*m2*v0*w0*pow(x,3)*z0 - i2*r2*v0*w0*pow(x,3)*z0 - 
   i2*r0*v2*w0*pow(x,3)*z0 - i2*r0*v0*w2*pow(x,3)*z0 + 
   k2*n0*r0*v0*pow(x,4)*z0 - j2*n0*pow(v0,2)*pow(x,4)*z0 + 
   i2*s0*pow(v0,2)*pow(x,4)*z0 - i2*r0*v0*w0*pow(x,4)*z0 + 
   f1*k2*m2*p2*x*y*z0 - e2*k2*n1*p2*x*y*z0 - e1*k2*n2*p2*x*y*z0 - 
   f2*pow(k2,2)*r1*x*y*z0 - f1*pow(k2,2)*r2*x*y*z0 + 
   e2*pow(k2,2)*s1*x*y*z0 + e1*pow(k2,2)*s2*x*y*z0 - f1*j2*m2*t2*x*y*z0 + 
   e2*j2*n1*t2*x*y*z0 + e1*j2*n2*t2*x*y*z0 + f2*i2*r1*t2*x*y*z0 - 
   m2*n2*r1*t2*x*y*z0 + f1*i2*r2*t2*x*y*z0 - m2*n1*r2*t2*x*y*z0 - 
   e2*i2*s1*t2*x*y*z0 + pow(m2,2)*s1*t2*x*y*z0 - e1*i2*s2*t2*x*y*z0 + 
   f2*j2*k2*v1*x*y*z0 - f2*i2*p2*v1*x*y*z0 + m2*n2*p2*v1*x*y*z0 + 
   k2*n2*r2*v1*x*y*z0 - 2*k2*m2*s2*v1*x*y*z0 + f1*j2*k2*v2*x*y*z0 - 
   f1*i2*p2*v2*x*y*z0 + m2*n1*p2*v2*x*y*z0 + k2*n2*r1*v2*x*y*z0 + 
   k2*n1*r2*v2*x*y*z0 - 2*k2*m2*s1*v2*x*y*z0 - 2*j2*n2*v1*v2*x*y*z0 + 
   2*i2*s2*v1*v2*x*y*z0 - j2*n1*pow(v2,2)*x*y*z0 + 
   i2*s1*pow(v2,2)*x*y*z0 - e2*j2*k2*w1*x*y*z0 + e2*i2*p2*w1*x*y*z0 - 
   pow(m2,2)*p2*w1*x*y*z0 + k2*m2*r2*w1*x*y*z0 + j2*m2*v2*w1*x*y*z0 - 
   i2*r2*v2*w1*x*y*z0 - e1*j2*k2*w2*x*y*z0 + e1*i2*p2*w2*x*y*z0 + 
   k2*m2*r1*w2*x*y*z0 + j2*m2*v1*w2*x*y*z0 - i2*r2*v1*w2*x*y*z0 - 
   i2*r1*v2*w2*x*y*z0 - e1*k2*n0*p2*pow(x,2)*y*z0 - 
   e0*k2*n1*p2*pow(x,2)*y*z0 - f1*pow(k2,2)*r0*pow(x,2)*y*z0 - 
   f0*pow(k2,2)*r1*pow(x,2)*y*z0 + e1*pow(k2,2)*s0*pow(x,2)*y*z0 + 
   e0*pow(k2,2)*s1*pow(x,2)*y*z0 + e1*j2*n0*t2*pow(x,2)*y*z0 + 
   e0*j2*n1*t2*pow(x,2)*y*z0 + f1*i2*r0*t2*pow(x,2)*y*z0 - 
   m2*n1*r0*t2*pow(x,2)*y*z0 + f0*i2*r1*t2*pow(x,2)*y*z0 - 
   m2*n0*r1*t2*pow(x,2)*y*z0 - e1*i2*s0*t2*pow(x,2)*y*z0 - 
   e0*i2*s1*t2*pow(x,2)*y*z0 + f1*j2*k2*v0*pow(x,2)*y*z0 - 
   f1*i2*p2*v0*pow(x,2)*y*z0 + m2*n1*p2*v0*pow(x,2)*y*z0 + 
   k2*n2*r1*v0*pow(x,2)*y*z0 + k2*n1*r2*v0*pow(x,2)*y*z0 - 
   2*k2*m2*s1*v0*pow(x,2)*y*z0 + f0*j2*k2*v1*pow(x,2)*y*z0 - 
   f0*i2*p2*v1*pow(x,2)*y*z0 + m2*n0*p2*v1*pow(x,2)*y*z0 + 
   k2*n2*r0*v1*pow(x,2)*y*z0 + k2*n0*r2*v1*pow(x,2)*y*z0 - 
   2*k2*m2*s0*v1*pow(x,2)*y*z0 - 2*j2*n2*v0*v1*pow(x,2)*y*z0 + 
   2*i2*s2*v0*v1*pow(x,2)*y*z0 + k2*n1*r0*v2*pow(x,2)*y*z0 + 
   k2*n0*r1*v2*pow(x,2)*y*z0 - 2*j2*n1*v0*v2*pow(x,2)*y*z0 + 
   2*i2*s1*v0*v2*pow(x,2)*y*z0 - 2*j2*n0*v1*v2*pow(x,2)*y*z0 + 
   2*i2*s0*v1*v2*pow(x,2)*y*z0 - e1*j2*k2*w0*pow(x,2)*y*z0 + 
   e1*i2*p2*w0*pow(x,2)*y*z0 + k2*m2*r1*w0*pow(x,2)*y*z0 + 
   j2*m2*v1*w0*pow(x,2)*y*z0 - i2*r2*v1*w0*pow(x,2)*y*z0 - 
   i2*r1*v2*w0*pow(x,2)*y*z0 - e0*j2*k2*w1*pow(x,2)*y*z0 + 
   e0*i2*p2*w1*pow(x,2)*y*z0 + k2*m2*r0*w1*pow(x,2)*y*z0 + 
   j2*m2*v0*w1*pow(x,2)*y*z0 - i2*r2*v0*w1*pow(x,2)*y*z0 - 
   i2*r0*v2*w1*pow(x,2)*y*z0 - i2*r1*v0*w2*pow(x,2)*y*z0 - 
   i2*r0*v1*w2*pow(x,2)*y*z0 + k2*n1*r0*v0*pow(x,3)*y*z0 + 
   k2*n0*r1*v0*pow(x,3)*y*z0 - j2*n1*pow(v0,2)*pow(x,3)*y*z0 + 
   i2*s1*pow(v0,2)*pow(x,3)*y*z0 + k2*n0*r0*v1*pow(x,3)*y*z0 - 
   2*j2*n0*v0*v1*pow(x,3)*y*z0 + 2*i2*s0*v0*v1*pow(x,3)*y*z0 - 
   i2*r1*v0*w0*pow(x,3)*y*z0 - i2*r0*v1*w0*pow(x,3)*y*z0 - 
   i2*r0*v0*w1*pow(x,3)*y*z0 - e1*k2*n1*p2*x*pow(y,2)*z0 - 
   f1*pow(k2,2)*r1*x*pow(y,2)*z0 + e1*pow(k2,2)*s1*x*pow(y,2)*z0 + 
   e1*j2*n1*t2*x*pow(y,2)*z0 + f1*i2*r1*t2*x*pow(y,2)*z0 - 
   m2*n1*r1*t2*x*pow(y,2)*z0 - e1*i2*s1*t2*x*pow(y,2)*z0 + 
   f1*j2*k2*v1*x*pow(y,2)*z0 - f1*i2*p2*v1*x*pow(y,2)*z0 + 
   m2*n1*p2*v1*x*pow(y,2)*z0 + k2*n2*r1*v1*x*pow(y,2)*z0 + 
   k2*n1*r2*v1*x*pow(y,2)*z0 - 2*k2*m2*s1*v1*x*pow(y,2)*z0 - 
   j2*n2*pow(v1,2)*x*pow(y,2)*z0 + i2*s2*pow(v1,2)*x*pow(y,2)*z0 + 
   k2*n1*r1*v2*x*pow(y,2)*z0 - 2*j2*n1*v1*v2*x*pow(y,2)*z0 + 
   2*i2*s1*v1*v2*x*pow(y,2)*z0 - e1*j2*k2*w1*x*pow(y,2)*z0 + 
   e1*i2*p2*w1*x*pow(y,2)*z0 + k2*m2*r1*w1*x*pow(y,2)*z0 + 
   j2*m2*v1*w1*x*pow(y,2)*z0 - i2*r2*v1*w1*x*pow(y,2)*z0 - 
   i2*r1*v2*w1*x*pow(y,2)*z0 - i2*r1*v1*w2*x*pow(y,2)*z0 + 
   k2*n1*r1*v0*pow(x,2)*pow(y,2)*z0 + 
   k2*n1*r0*v1*pow(x,2)*pow(y,2)*z0 + 
   k2*n0*r1*v1*pow(x,2)*pow(y,2)*z0 - 
   2*j2*n1*v0*v1*pow(x,2)*pow(y,2)*z0 + 
   2*i2*s1*v0*v1*pow(x,2)*pow(y,2)*z0 - 
   j2*n0*pow(v1,2)*pow(x,2)*pow(y,2)*z0 + 
   i2*s0*pow(v1,2)*pow(x,2)*pow(y,2)*z0 - 
   i2*r1*v1*w0*pow(x,2)*pow(y,2)*z0 - 
   i2*r1*v0*w1*pow(x,2)*pow(y,2)*z0 - 
   i2*r0*v1*w1*pow(x,2)*pow(y,2)*z0 + k2*n1*r1*v1*x*pow(y,3)*z0 - 
   j2*n1*pow(v1,2)*x*pow(y,3)*z0 + i2*s1*pow(v1,2)*x*pow(y,3)*z0 - 
   i2*r1*v1*w1*x*pow(y,3)*z0 + f2*k2*m2*p2*y*z1 - e2*k2*n2*p2*y*z1 - 
   f2*pow(k2,2)*r2*y*z1 + e2*pow(k2,2)*s2*y*z1 - f2*j2*m2*t2*y*z1 + 
   e2*j2*n2*t2*y*z1 + f2*i2*r2*t2*y*z1 - m2*n2*r2*t2*y*z1 - 
   e2*i2*s2*t2*y*z1 + pow(m2,2)*s2*t2*y*z1 + f2*j2*k2*v2*y*z1 - 
   f2*i2*p2*v2*y*z1 + m2*n2*p2*v2*y*z1 + k2*n2*r2*v2*y*z1 - 
   2*k2*m2*s2*v2*y*z1 - j2*n2*pow(v2,2)*y*z1 + i2*s2*pow(v2,2)*y*z1 - 
   e2*j2*k2*w2*y*z1 + e2*i2*p2*w2*y*z1 - pow(m2,2)*p2*w2*y*z1 + 
   k2*m2*r2*w2*y*z1 + j2*m2*v2*w2*y*z1 - i2*r2*v2*w2*y*z1 + 
   f0*k2*m2*p2*x*y*z1 - e2*k2*n0*p2*x*y*z1 - e0*k2*n2*p2*x*y*z1 - 
   f2*pow(k2,2)*r0*x*y*z1 - f0*pow(k2,2)*r2*x*y*z1 + 
   e2*pow(k2,2)*s0*x*y*z1 + e0*pow(k2,2)*s2*x*y*z1 - f0*j2*m2*t2*x*y*z1 + 
   e2*j2*n0*t2*x*y*z1 + e0*j2*n2*t2*x*y*z1 + f2*i2*r0*t2*x*y*z1 - 
   m2*n2*r0*t2*x*y*z1 + f0*i2*r2*t2*x*y*z1 - m2*n0*r2*t2*x*y*z1 - 
   e2*i2*s0*t2*x*y*z1 + pow(m2,2)*s0*t2*x*y*z1 - e0*i2*s2*t2*x*y*z1 + 
   f2*j2*k2*v0*x*y*z1 - f2*i2*p2*v0*x*y*z1 + m2*n2*p2*v0*x*y*z1 + 
   k2*n2*r2*v0*x*y*z1 - 2*k2*m2*s2*v0*x*y*z1 + f0*j2*k2*v2*x*y*z1 - 
   f0*i2*p2*v2*x*y*z1 + m2*n0*p2*v2*x*y*z1 + k2*n2*r0*v2*x*y*z1 + 
   k2*n0*r2*v2*x*y*z1 - 2*k2*m2*s0*v2*x*y*z1 - 2*j2*n2*v0*v2*x*y*z1 + 
   2*i2*s2*v0*v2*x*y*z1 - j2*n0*pow(v2,2)*x*y*z1 + 
   i2*s0*pow(v2,2)*x*y*z1 - e2*j2*k2*w0*x*y*z1 + e2*i2*p2*w0*x*y*z1 - 
   pow(m2,2)*p2*w0*x*y*z1 + k2*m2*r2*w0*x*y*z1 + j2*m2*v2*w0*x*y*z1 - 
   i2*r2*v2*w0*x*y*z1 - e0*j2*k2*w2*x*y*z1 + e0*i2*p2*w2*x*y*z1 + 
   k2*m2*r0*w2*x*y*z1 + j2*m2*v0*w2*x*y*z1 - i2*r2*v0*w2*x*y*z1 - 
   i2*r0*v2*w2*x*y*z1 - e0*k2*n0*p2*pow(x,2)*y*z1 - 
   f0*pow(k2,2)*r0*pow(x,2)*y*z1 + e0*pow(k2,2)*s0*pow(x,2)*y*z1 + 
   e0*j2*n0*t2*pow(x,2)*y*z1 + f0*i2*r0*t2*pow(x,2)*y*z1 - 
   m2*n0*r0*t2*pow(x,2)*y*z1 - e0*i2*s0*t2*pow(x,2)*y*z1 + 
   f0*j2*k2*v0*pow(x,2)*y*z1 - f0*i2*p2*v0*pow(x,2)*y*z1 + 
   m2*n0*p2*v0*pow(x,2)*y*z1 + k2*n2*r0*v0*pow(x,2)*y*z1 + 
   k2*n0*r2*v0*pow(x,2)*y*z1 - 2*k2*m2*s0*v0*pow(x,2)*y*z1 - 
   j2*n2*pow(v0,2)*pow(x,2)*y*z1 + i2*s2*pow(v0,2)*pow(x,2)*y*z1 + 
   k2*n0*r0*v2*pow(x,2)*y*z1 - 2*j2*n0*v0*v2*pow(x,2)*y*z1 + 
   2*i2*s0*v0*v2*pow(x,2)*y*z1 - e0*j2*k2*w0*pow(x,2)*y*z1 + 
   e0*i2*p2*w0*pow(x,2)*y*z1 + k2*m2*r0*w0*pow(x,2)*y*z1 + 
   j2*m2*v0*w0*pow(x,2)*y*z1 - i2*r2*v0*w0*pow(x,2)*y*z1 - 
   i2*r0*v2*w0*pow(x,2)*y*z1 - i2*r0*v0*w2*pow(x,2)*y*z1 + 
   k2*n0*r0*v0*pow(x,3)*y*z1 - j2*n0*pow(v0,2)*pow(x,3)*y*z1 + 
   i2*s0*pow(v0,2)*pow(x,3)*y*z1 - i2*r0*v0*w0*pow(x,3)*y*z1 + 
   f1*k2*m2*p2*pow(y,2)*z1 - e2*k2*n1*p2*pow(y,2)*z1 - 
   e1*k2*n2*p2*pow(y,2)*z1 - f2*pow(k2,2)*r1*pow(y,2)*z1 - 
   f1*pow(k2,2)*r2*pow(y,2)*z1 + e2*pow(k2,2)*s1*pow(y,2)*z1 + 
   e1*pow(k2,2)*s2*pow(y,2)*z1 - f1*j2*m2*t2*pow(y,2)*z1 + 
   e2*j2*n1*t2*pow(y,2)*z1 + e1*j2*n2*t2*pow(y,2)*z1 + 
   f2*i2*r1*t2*pow(y,2)*z1 - m2*n2*r1*t2*pow(y,2)*z1 + 
   f1*i2*r2*t2*pow(y,2)*z1 - m2*n1*r2*t2*pow(y,2)*z1 - 
   e2*i2*s1*t2*pow(y,2)*z1 + pow(m2,2)*s1*t2*pow(y,2)*z1 - 
   e1*i2*s2*t2*pow(y,2)*z1 + f2*j2*k2*v1*pow(y,2)*z1 - 
   f2*i2*p2*v1*pow(y,2)*z1 + m2*n2*p2*v1*pow(y,2)*z1 + 
   k2*n2*r2*v1*pow(y,2)*z1 - 2*k2*m2*s2*v1*pow(y,2)*z1 + 
   f1*j2*k2*v2*pow(y,2)*z1 - f1*i2*p2*v2*pow(y,2)*z1 + 
   m2*n1*p2*v2*pow(y,2)*z1 + k2*n2*r1*v2*pow(y,2)*z1 + 
   k2*n1*r2*v2*pow(y,2)*z1 - 2*k2*m2*s1*v2*pow(y,2)*z1 - 
   2*j2*n2*v1*v2*pow(y,2)*z1 + 2*i2*s2*v1*v2*pow(y,2)*z1 - 
   j2*n1*pow(v2,2)*pow(y,2)*z1 + i2*s1*pow(v2,2)*pow(y,2)*z1 - 
   e2*j2*k2*w1*pow(y,2)*z1 + e2*i2*p2*w1*pow(y,2)*z1 - 
   pow(m2,2)*p2*w1*pow(y,2)*z1 + k2*m2*r2*w1*pow(y,2)*z1 + 
   j2*m2*v2*w1*pow(y,2)*z1 - i2*r2*v2*w1*pow(y,2)*z1 - 
   e1*j2*k2*w2*pow(y,2)*z1 + e1*i2*p2*w2*pow(y,2)*z1 + 
   k2*m2*r1*w2*pow(y,2)*z1 + j2*m2*v1*w2*pow(y,2)*z1 - 
   i2*r2*v1*w2*pow(y,2)*z1 - i2*r1*v2*w2*pow(y,2)*z1 - 
   e1*k2*n0*p2*x*pow(y,2)*z1 - e0*k2*n1*p2*x*pow(y,2)*z1 - 
   f1*pow(k2,2)*r0*x*pow(y,2)*z1 - f0*pow(k2,2)*r1*x*pow(y,2)*z1 + 
   e1*pow(k2,2)*s0*x*pow(y,2)*z1 + e0*pow(k2,2)*s1*x*pow(y,2)*z1 + 
   e1*j2*n0*t2*x*pow(y,2)*z1 + e0*j2*n1*t2*x*pow(y,2)*z1 + 
   f1*i2*r0*t2*x*pow(y,2)*z1 - m2*n1*r0*t2*x*pow(y,2)*z1 + 
   f0*i2*r1*t2*x*pow(y,2)*z1 - m2*n0*r1*t2*x*pow(y,2)*z1 - 
   e1*i2*s0*t2*x*pow(y,2)*z1 - e0*i2*s1*t2*x*pow(y,2)*z1 + 
   f1*j2*k2*v0*x*pow(y,2)*z1 - f1*i2*p2*v0*x*pow(y,2)*z1 + 
   m2*n1*p2*v0*x*pow(y,2)*z1 + k2*n2*r1*v0*x*pow(y,2)*z1 + 
   k2*n1*r2*v0*x*pow(y,2)*z1 - 2*k2*m2*s1*v0*x*pow(y,2)*z1 + 
   f0*j2*k2*v1*x*pow(y,2)*z1 - f0*i2*p2*v1*x*pow(y,2)*z1 + 
   m2*n0*p2*v1*x*pow(y,2)*z1 + k2*n2*r0*v1*x*pow(y,2)*z1 + 
   k2*n0*r2*v1*x*pow(y,2)*z1 - 2*k2*m2*s0*v1*x*pow(y,2)*z1 - 
   2*j2*n2*v0*v1*x*pow(y,2)*z1 + 2*i2*s2*v0*v1*x*pow(y,2)*z1 + 
   k2*n1*r0*v2*x*pow(y,2)*z1 + k2*n0*r1*v2*x*pow(y,2)*z1 - 
   2*j2*n1*v0*v2*x*pow(y,2)*z1 + 2*i2*s1*v0*v2*x*pow(y,2)*z1 - 
   2*j2*n0*v1*v2*x*pow(y,2)*z1 + 2*i2*s0*v1*v2*x*pow(y,2)*z1 - 
   e1*j2*k2*w0*x*pow(y,2)*z1 + e1*i2*p2*w0*x*pow(y,2)*z1 + 
   k2*m2*r1*w0*x*pow(y,2)*z1 + j2*m2*v1*w0*x*pow(y,2)*z1 - 
   i2*r2*v1*w0*x*pow(y,2)*z1 - i2*r1*v2*w0*x*pow(y,2)*z1 - 
   e0*j2*k2*w1*x*pow(y,2)*z1 + e0*i2*p2*w1*x*pow(y,2)*z1 + 
   k2*m2*r0*w1*x*pow(y,2)*z1 + j2*m2*v0*w1*x*pow(y,2)*z1 - 
   i2*r2*v0*w1*x*pow(y,2)*z1 - i2*r0*v2*w1*x*pow(y,2)*z1 - 
   i2*r1*v0*w2*x*pow(y,2)*z1 - i2*r0*v1*w2*x*pow(y,2)*z1 + 
   k2*n1*r0*v0*pow(x,2)*pow(y,2)*z1 + 
   k2*n0*r1*v0*pow(x,2)*pow(y,2)*z1 - 
   j2*n1*pow(v0,2)*pow(x,2)*pow(y,2)*z1 + 
   i2*s1*pow(v0,2)*pow(x,2)*pow(y,2)*z1 + 
   k2*n0*r0*v1*pow(x,2)*pow(y,2)*z1 - 
   2*j2*n0*v0*v1*pow(x,2)*pow(y,2)*z1 + 
   2*i2*s0*v0*v1*pow(x,2)*pow(y,2)*z1 - 
   i2*r1*v0*w0*pow(x,2)*pow(y,2)*z1 - 
   i2*r0*v1*w0*pow(x,2)*pow(y,2)*z1 - 
   i2*r0*v0*w1*pow(x,2)*pow(y,2)*z1 - e1*k2*n1*p2*pow(y,3)*z1 - 
   f1*pow(k2,2)*r1*pow(y,3)*z1 + e1*pow(k2,2)*s1*pow(y,3)*z1 + 
   e1*j2*n1*t2*pow(y,3)*z1 + f1*i2*r1*t2*pow(y,3)*z1 - 
   m2*n1*r1*t2*pow(y,3)*z1 - e1*i2*s1*t2*pow(y,3)*z1 + 
   f1*j2*k2*v1*pow(y,3)*z1 - f1*i2*p2*v1*pow(y,3)*z1 + 
   m2*n1*p2*v1*pow(y,3)*z1 + k2*n2*r1*v1*pow(y,3)*z1 + 
   k2*n1*r2*v1*pow(y,3)*z1 - 2*k2*m2*s1*v1*pow(y,3)*z1 - 
   j2*n2*pow(v1,2)*pow(y,3)*z1 + i2*s2*pow(v1,2)*pow(y,3)*z1 + 
   k2*n1*r1*v2*pow(y,3)*z1 - 2*j2*n1*v1*v2*pow(y,3)*z1 + 
   2*i2*s1*v1*v2*pow(y,3)*z1 - e1*j2*k2*w1*pow(y,3)*z1 + 
   e1*i2*p2*w1*pow(y,3)*z1 + k2*m2*r1*w1*pow(y,3)*z1 + 
   j2*m2*v1*w1*pow(y,3)*z1 - i2*r2*v1*w1*pow(y,3)*z1 - 
   i2*r1*v2*w1*pow(y,3)*z1 - i2*r1*v1*w2*pow(y,3)*z1 + 
   k2*n1*r1*v0*x*pow(y,3)*z1 + k2*n1*r0*v1*x*pow(y,3)*z1 + 
   k2*n0*r1*v1*x*pow(y,3)*z1 - 2*j2*n1*v0*v1*x*pow(y,3)*z1 + 
   2*i2*s1*v0*v1*x*pow(y,3)*z1 - j2*n0*pow(v1,2)*x*pow(y,3)*z1 + 
   i2*s0*pow(v1,2)*x*pow(y,3)*z1 - i2*r1*v1*w0*x*pow(y,3)*z1 - 
   i2*r1*v0*w1*x*pow(y,3)*z1 - i2*r0*v1*w1*x*pow(y,3)*z1 + 
   k2*n1*r1*v1*pow(y,4)*z1 - j2*n1*pow(v1,2)*pow(y,4)*z1 + 
   i2*s1*pow(v1,2)*pow(y,4)*z1 - i2*r1*v1*w1*pow(y,4)*z1 + 
   f2*k2*m2*p2*z2 - e2*k2*n2*p2*z2 - f2*pow(k2,2)*r2*z2 + 
   e2*pow(k2,2)*s2*z2 - f2*j2*m2*t2*z2 + e2*j2*n2*t2*z2 + f2*i2*r2*t2*z2 - 
   m2*n2*r2*t2*z2 - e2*i2*s2*t2*z2 + pow(m2,2)*s2*t2*z2 + f2*j2*k2*v2*z2 - 
   f2*i2*p2*v2*z2 + m2*n2*p2*v2*z2 + k2*n2*r2*v2*z2 - 2*k2*m2*s2*v2*z2 - 
   j2*n2*pow(v2,2)*z2 + i2*s2*pow(v2,2)*z2 - e2*j2*k2*w2*z2 + 
   e2*i2*p2*w2*z2 - pow(m2,2)*p2*w2*z2 + k2*m2*r2*w2*z2 + j2*m2*v2*w2*z2 - 
   i2*r2*v2*w2*z2 + f0*k2*m2*p2*x*z2 - e2*k2*n0*p2*x*z2 - e0*k2*n2*p2*x*z2 - 
   f2*pow(k2,2)*r0*x*z2 - f0*pow(k2,2)*r2*x*z2 + e2*pow(k2,2)*s0*x*z2 + 
   e0*pow(k2,2)*s2*x*z2 - f0*j2*m2*t2*x*z2 + e2*j2*n0*t2*x*z2 + 
   e0*j2*n2*t2*x*z2 + f2*i2*r0*t2*x*z2 - m2*n2*r0*t2*x*z2 + 
   f0*i2*r2*t2*x*z2 - m2*n0*r2*t2*x*z2 - e2*i2*s0*t2*x*z2 + 
   pow(m2,2)*s0*t2*x*z2 - e0*i2*s2*t2*x*z2 + f2*j2*k2*v0*x*z2 - 
   f2*i2*p2*v0*x*z2 + m2*n2*p2*v0*x*z2 + k2*n2*r2*v0*x*z2 - 
   2*k2*m2*s2*v0*x*z2 + f0*j2*k2*v2*x*z2 - f0*i2*p2*v2*x*z2 + 
   m2*n0*p2*v2*x*z2 + k2*n2*r0*v2*x*z2 + k2*n0*r2*v2*x*z2 - 
   2*k2*m2*s0*v2*x*z2 - 2*j2*n2*v0*v2*x*z2 + 2*i2*s2*v0*v2*x*z2 - 
   j2*n0*pow(v2,2)*x*z2 + i2*s0*pow(v2,2)*x*z2 - e2*j2*k2*w0*x*z2 + 
   e2*i2*p2*w0*x*z2 - pow(m2,2)*p2*w0*x*z2 + k2*m2*r2*w0*x*z2 + 
   j2*m2*v2*w0*x*z2 - i2*r2*v2*w0*x*z2 - e0*j2*k2*w2*x*z2 + 
   e0*i2*p2*w2*x*z2 + k2*m2*r0*w2*x*z2 + j2*m2*v0*w2*x*z2 - 
   i2*r2*v0*w2*x*z2 - i2*r0*v2*w2*x*z2 - e0*k2*n0*p2*pow(x,2)*z2 - 
   f0*pow(k2,2)*r0*pow(x,2)*z2 + e0*pow(k2,2)*s0*pow(x,2)*z2 + 
   e0*j2*n0*t2*pow(x,2)*z2 + f0*i2*r0*t2*pow(x,2)*z2 - 
   m2*n0*r0*t2*pow(x,2)*z2 - e0*i2*s0*t2*pow(x,2)*z2 + 
   f0*j2*k2*v0*pow(x,2)*z2 - f0*i2*p2*v0*pow(x,2)*z2 + 
   m2*n0*p2*v0*pow(x,2)*z2 + k2*n2*r0*v0*pow(x,2)*z2 + 
   k2*n0*r2*v0*pow(x,2)*z2 - 2*k2*m2*s0*v0*pow(x,2)*z2 - 
   j2*n2*pow(v0,2)*pow(x,2)*z2 + i2*s2*pow(v0,2)*pow(x,2)*z2 + 
   k2*n0*r0*v2*pow(x,2)*z2 - 2*j2*n0*v0*v2*pow(x,2)*z2 + 
   2*i2*s0*v0*v2*pow(x,2)*z2 - e0*j2*k2*w0*pow(x,2)*z2 + 
   e0*i2*p2*w0*pow(x,2)*z2 + k2*m2*r0*w0*pow(x,2)*z2 + 
   j2*m2*v0*w0*pow(x,2)*z2 - i2*r2*v0*w0*pow(x,2)*z2 - 
   i2*r0*v2*w0*pow(x,2)*z2 - i2*r0*v0*w2*pow(x,2)*z2 + 
   k2*n0*r0*v0*pow(x,3)*z2 - j2*n0*pow(v0,2)*pow(x,3)*z2 + 
   i2*s0*pow(v0,2)*pow(x,3)*z2 - i2*r0*v0*w0*pow(x,3)*z2 + 
   f1*k2*m2*p2*y*z2 - e2*k2*n1*p2*y*z2 - e1*k2*n2*p2*y*z2 - 
   f2*pow(k2,2)*r1*y*z2 - f1*pow(k2,2)*r2*y*z2 + e2*pow(k2,2)*s1*y*z2 + 
   e1*pow(k2,2)*s2*y*z2 - f1*j2*m2*t2*y*z2 + e2*j2*n1*t2*y*z2 + 
   e1*j2*n2*t2*y*z2 + f2*i2*r1*t2*y*z2 - m2*n2*r1*t2*y*z2 + 
   f1*i2*r2*t2*y*z2 - m2*n1*r2*t2*y*z2 - e2*i2*s1*t2*y*z2 + 
   pow(m2,2)*s1*t2*y*z2 - e1*i2*s2*t2*y*z2 + f2*j2*k2*v1*y*z2 - 
   f2*i2*p2*v1*y*z2 + m2*n2*p2*v1*y*z2 + k2*n2*r2*v1*y*z2 - 
   2*k2*m2*s2*v1*y*z2 + f1*j2*k2*v2*y*z2 - f1*i2*p2*v2*y*z2 + 
   m2*n1*p2*v2*y*z2 + k2*n2*r1*v2*y*z2 + k2*n1*r2*v2*y*z2 - 
   2*k2*m2*s1*v2*y*z2 - 2*j2*n2*v1*v2*y*z2 + 2*i2*s2*v1*v2*y*z2 - 
   j2*n1*pow(v2,2)*y*z2 + i2*s1*pow(v2,2)*y*z2 - e2*j2*k2*w1*y*z2 + 
   e2*i2*p2*w1*y*z2 - pow(m2,2)*p2*w1*y*z2 + k2*m2*r2*w1*y*z2 + 
   j2*m2*v2*w1*y*z2 - i2*r2*v2*w1*y*z2 - e1*j2*k2*w2*y*z2 + 
   e1*i2*p2*w2*y*z2 + k2*m2*r1*w2*y*z2 + j2*m2*v1*w2*y*z2 - 
   i2*r2*v1*w2*y*z2 - i2*r1*v2*w2*y*z2 - e1*k2*n0*p2*x*y*z2 - 
   e0*k2*n1*p2*x*y*z2 - f1*pow(k2,2)*r0*x*y*z2 - f0*pow(k2,2)*r1*x*y*z2 + 
   e1*pow(k2,2)*s0*x*y*z2 + e0*pow(k2,2)*s1*x*y*z2 + e1*j2*n0*t2*x*y*z2 + 
   e0*j2*n1*t2*x*y*z2 + f1*i2*r0*t2*x*y*z2 - m2*n1*r0*t2*x*y*z2 + 
   f0*i2*r1*t2*x*y*z2 - m2*n0*r1*t2*x*y*z2 - e1*i2*s0*t2*x*y*z2 - 
   e0*i2*s1*t2*x*y*z2 + f1*j2*k2*v0*x*y*z2 - f1*i2*p2*v0*x*y*z2 + 
   m2*n1*p2*v0*x*y*z2 + k2*n2*r1*v0*x*y*z2 + k2*n1*r2*v0*x*y*z2 - 
   2*k2*m2*s1*v0*x*y*z2 + f0*j2*k2*v1*x*y*z2 - f0*i2*p2*v1*x*y*z2 + 
   m2*n0*p2*v1*x*y*z2 + k2*n2*r0*v1*x*y*z2 + k2*n0*r2*v1*x*y*z2 - 
   2*k2*m2*s0*v1*x*y*z2 - 2*j2*n2*v0*v1*x*y*z2 + 2*i2*s2*v0*v1*x*y*z2 + 
   k2*n1*r0*v2*x*y*z2 + k2*n0*r1*v2*x*y*z2 - 2*j2*n1*v0*v2*x*y*z2 + 
   2*i2*s1*v0*v2*x*y*z2 - 2*j2*n0*v1*v2*x*y*z2 + 2*i2*s0*v1*v2*x*y*z2 - 
   e1*j2*k2*w0*x*y*z2 + e1*i2*p2*w0*x*y*z2 + k2*m2*r1*w0*x*y*z2 + 
   j2*m2*v1*w0*x*y*z2 - i2*r2*v1*w0*x*y*z2 - i2*r1*v2*w0*x*y*z2 - 
   e0*j2*k2*w1*x*y*z2 + e0*i2*p2*w1*x*y*z2 + k2*m2*r0*w1*x*y*z2 + 
   j2*m2*v0*w1*x*y*z2 - i2*r2*v0*w1*x*y*z2 - i2*r0*v2*w1*x*y*z2 - 
   i2*r1*v0*w2*x*y*z2 - i2*r0*v1*w2*x*y*z2 + k2*n1*r0*v0*pow(x,2)*y*z2 + 
   k2*n0*r1*v0*pow(x,2)*y*z2 - j2*n1*pow(v0,2)*pow(x,2)*y*z2 + 
   i2*s1*pow(v0,2)*pow(x,2)*y*z2 + k2*n0*r0*v1*pow(x,2)*y*z2 - 
   2*j2*n0*v0*v1*pow(x,2)*y*z2 + 2*i2*s0*v0*v1*pow(x,2)*y*z2 - 
   i2*r1*v0*w0*pow(x,2)*y*z2 - i2*r0*v1*w0*pow(x,2)*y*z2 - 
   i2*r0*v0*w1*pow(x,2)*y*z2 - e1*k2*n1*p2*pow(y,2)*z2 - 
   f1*pow(k2,2)*r1*pow(y,2)*z2 + e1*pow(k2,2)*s1*pow(y,2)*z2 + 
   e1*j2*n1*t2*pow(y,2)*z2 + f1*i2*r1*t2*pow(y,2)*z2 - 
   m2*n1*r1*t2*pow(y,2)*z2 - e1*i2*s1*t2*pow(y,2)*z2 + 
   f1*j2*k2*v1*pow(y,2)*z2 - f1*i2*p2*v1*pow(y,2)*z2 + 
   m2*n1*p2*v1*pow(y,2)*z2 + k2*n2*r1*v1*pow(y,2)*z2 + 
   k2*n1*r2*v1*pow(y,2)*z2 - 2*k2*m2*s1*v1*pow(y,2)*z2 - 
   j2*n2*pow(v1,2)*pow(y,2)*z2 + i2*s2*pow(v1,2)*pow(y,2)*z2 + 
   k2*n1*r1*v2*pow(y,2)*z2 - 2*j2*n1*v1*v2*pow(y,2)*z2 + 
   2*i2*s1*v1*v2*pow(y,2)*z2 - e1*j2*k2*w1*pow(y,2)*z2 + 
   e1*i2*p2*w1*pow(y,2)*z2 + k2*m2*r1*w1*pow(y,2)*z2 + 
   j2*m2*v1*w1*pow(y,2)*z2 - i2*r2*v1*w1*pow(y,2)*z2 - 
   i2*r1*v2*w1*pow(y,2)*z2 - i2*r1*v1*w2*pow(y,2)*z2 + 
   k2*n1*r1*v0*x*pow(y,2)*z2 + k2*n1*r0*v1*x*pow(y,2)*z2 + 
   k2*n0*r1*v1*x*pow(y,2)*z2 - 2*j2*n1*v0*v1*x*pow(y,2)*z2 + 
   2*i2*s1*v0*v1*x*pow(y,2)*z2 - j2*n0*pow(v1,2)*x*pow(y,2)*z2 + 
   i2*s0*pow(v1,2)*x*pow(y,2)*z2 - i2*r1*v1*w0*x*pow(y,2)*z2 - 
   i2*r1*v0*w1*x*pow(y,2)*z2 - i2*r0*v1*w1*x*pow(y,2)*z2 + 
   k2*n1*r1*v1*pow(y,3)*z2 - j2*n1*pow(v1,2)*pow(y,3)*z2 + 
   i2*s1*pow(v1,2)*pow(y,3)*z2 - i2*r1*v1*w1*pow(y,3)*z2
;

	rg_REAL t_5 = d2*e2*l2*pow(p2,2) - c2*f2*l2*pow(p2,2) - c2*d2*m2*pow(p2,2) + 
   pow(c2,2)*n2*pow(p2,2) - d2*e2*k2*p2*q2 + c2*f2*k2*p2*q2 + 
   c2*d2*k2*p2*r2 - pow(c2,2)*k2*p2*s2 - d2*e2*l2*o2*t2 + c2*f2*l2*o2*t2 + 
   c2*d2*m2*o2*t2 - pow(c2,2)*n2*o2*t2 + d2*e2*j2*q2*t2 - c2*f2*j2*q2*t2 + 
   f2*m2*pow(q2,2)*t2 - e2*n2*pow(q2,2)*t2 - c2*d2*j2*r2*t2 - 
   f2*l2*q2*r2*t2 - d2*m2*q2*r2*t2 + 2*c2*n2*q2*r2*t2 + 
   d2*l2*pow(r2,2)*t2 + pow(c2,2)*j2*s2*t2 + e2*l2*q2*s2*t2 - 
   c2*m2*q2*s2*t2 - c2*l2*r2*s2*t2 + d2*e2*k2*o2*u2 - c2*f2*k2*o2*u2 - 
   d2*e2*j2*p2*u2 + c2*f2*j2*p2*u2 - 2*f2*m2*p2*q2*u2 + 2*e2*n2*p2*q2*u2 + 
   f2*l2*p2*r2*u2 + d2*m2*p2*r2*u2 - 2*c2*n2*p2*r2*u2 + f2*k2*q2*r2*u2 - 
   d2*k2*pow(r2,2)*u2 - e2*l2*p2*s2*u2 + c2*m2*p2*s2*u2 - e2*k2*q2*s2*u2 + 
   c2*k2*r2*s2*u2 + f2*m2*o2*pow(u2,2) - e2*n2*o2*pow(u2,2) - 
   f2*j2*r2*pow(u2,2) + n2*pow(r2,2)*pow(u2,2) + e2*j2*s2*pow(u2,2) - 
   m2*r2*s2*pow(u2,2) - c2*d2*k2*o2*v2 + c2*d2*j2*p2*v2 + f2*l2*p2*q2*v2 + 
   d2*m2*p2*q2*v2 - 2*c2*n2*p2*q2*v2 - f2*k2*pow(q2,2)*v2 - 
   2*d2*l2*p2*r2*v2 + d2*k2*q2*r2*v2 + c2*l2*p2*s2*v2 + c2*k2*q2*s2*v2 - 
   f2*l2*o2*u2*v2 - d2*m2*o2*u2*v2 + 2*c2*n2*o2*u2*v2 + f2*j2*q2*u2*v2 + 
   d2*j2*r2*u2*v2 - 2*n2*q2*r2*u2*v2 - 2*c2*j2*s2*u2*v2 + m2*q2*s2*u2*v2 + 
   l2*r2*s2*u2*v2 + d2*l2*o2*pow(v2,2) - d2*j2*q2*pow(v2,2) + 
   n2*pow(q2,2)*pow(v2,2) - l2*q2*s2*pow(v2,2) + pow(c2,2)*k2*o2*w2 - 
   pow(c2,2)*j2*p2*w2 - e2*l2*p2*q2*w2 + c2*m2*p2*q2*w2 + 
   e2*k2*pow(q2,2)*w2 + c2*l2*p2*r2*w2 - 2*c2*k2*q2*r2*w2 + 
   e2*l2*o2*u2*w2 - c2*m2*o2*u2*w2 - e2*j2*q2*u2*w2 + c2*j2*r2*u2*w2 + 
   m2*q2*r2*u2*w2 - l2*pow(r2,2)*u2*w2 - c2*l2*o2*v2*w2 + c2*j2*q2*v2*w2 - 
   m2*pow(q2,2)*v2*w2 + l2*q2*r2*v2*w2 + d2*e0*l2*pow(p2,2)*x + 
   d0*e2*l2*pow(p2,2)*x - c2*f0*l2*pow(p2,2)*x - c0*f2*l2*pow(p2,2)*x - 
   c2*d0*m2*pow(p2,2)*x - c0*d2*m2*pow(p2,2)*x + 
   pow(c2,2)*n0*pow(p2,2)*x + 2*c0*c2*n2*pow(p2,2)*x - 
   d2*e0*k2*p2*q2*x - d0*e2*k2*p2*q2*x + c2*f0*k2*p2*q2*x + 
   c0*f2*k2*p2*q2*x + c2*d2*k2*p2*r0*x + c2*d0*k2*p2*r2*x + 
   c0*d2*k2*p2*r2*x - pow(c2,2)*k2*p2*s0*x - 2*c0*c2*k2*p2*s2*x - 
   d2*e0*l2*o2*t2*x - d0*e2*l2*o2*t2*x + c2*f0*l2*o2*t2*x + 
   c0*f2*l2*o2*t2*x + c2*d0*m2*o2*t2*x + c0*d2*m2*o2*t2*x - 
   pow(c2,2)*n0*o2*t2*x - 2*c0*c2*n2*o2*t2*x + d2*e0*j2*q2*t2*x + 
   d0*e2*j2*q2*t2*x - c2*f0*j2*q2*t2*x - c0*f2*j2*q2*t2*x + 
   f0*m2*pow(q2,2)*t2*x - e2*n0*pow(q2,2)*t2*x - e0*n2*pow(q2,2)*t2*x - 
   c2*d2*j2*r0*t2*x - f2*l2*q2*r0*t2*x - d2*m2*q2*r0*t2*x + 
   2*c2*n2*q2*r0*t2*x - c2*d0*j2*r2*t2*x - c0*d2*j2*r2*t2*x - 
   f0*l2*q2*r2*t2*x - d0*m2*q2*r2*t2*x + 2*c2*n0*q2*r2*t2*x + 
   2*c0*n2*q2*r2*t2*x + 2*d2*l2*r0*r2*t2*x + d0*l2*pow(r2,2)*t2*x + 
   pow(c2,2)*j2*s0*t2*x + e2*l2*q2*s0*t2*x - c2*m2*q2*s0*t2*x - 
   c2*l2*r2*s0*t2*x + 2*c0*c2*j2*s2*t2*x + e0*l2*q2*s2*t2*x - 
   c0*m2*q2*s2*t2*x - c2*l2*r0*s2*t2*x - c0*l2*r2*s2*t2*x + 
   d2*e2*k2*o2*u0*x - c2*f2*k2*o2*u0*x - d2*e2*j2*p2*u0*x + 
   c2*f2*j2*p2*u0*x - 2*f2*m2*p2*q2*u0*x + 2*e2*n2*p2*q2*u0*x + 
   f2*l2*p2*r2*u0*x + d2*m2*p2*r2*u0*x - 2*c2*n2*p2*r2*u0*x + 
   f2*k2*q2*r2*u0*x - d2*k2*pow(r2,2)*u0*x - e2*l2*p2*s2*u0*x + 
   c2*m2*p2*s2*u0*x - e2*k2*q2*s2*u0*x + c2*k2*r2*s2*u0*x + 
   d2*e0*k2*o2*u2*x + d0*e2*k2*o2*u2*x - c2*f0*k2*o2*u2*x - 
   c0*f2*k2*o2*u2*x - d2*e0*j2*p2*u2*x - d0*e2*j2*p2*u2*x + 
   c2*f0*j2*p2*u2*x + c0*f2*j2*p2*u2*x - 2*f0*m2*p2*q2*u2*x + 
   2*e2*n0*p2*q2*u2*x + 2*e0*n2*p2*q2*u2*x + f2*l2*p2*r0*u2*x + 
   d2*m2*p2*r0*u2*x - 2*c2*n2*p2*r0*u2*x + f2*k2*q2*r0*u2*x + 
   f0*l2*p2*r2*u2*x + d0*m2*p2*r2*u2*x - 2*c2*n0*p2*r2*u2*x - 
   2*c0*n2*p2*r2*u2*x + f0*k2*q2*r2*u2*x - 2*d2*k2*r0*r2*u2*x - 
   d0*k2*pow(r2,2)*u2*x - e2*l2*p2*s0*u2*x + c2*m2*p2*s0*u2*x - 
   e2*k2*q2*s0*u2*x + c2*k2*r2*s0*u2*x - e0*l2*p2*s2*u2*x + 
   c0*m2*p2*s2*u2*x - e0*k2*q2*s2*u2*x + c2*k2*r0*s2*u2*x + 
   c0*k2*r2*s2*u2*x + 2*f2*m2*o2*u0*u2*x - 2*e2*n2*o2*u0*u2*x - 
   2*f2*j2*r2*u0*u2*x + 2*n2*pow(r2,2)*u0*u2*x + 2*e2*j2*s2*u0*u2*x - 
   2*m2*r2*s2*u0*u2*x + f0*m2*o2*pow(u2,2)*x - e2*n0*o2*pow(u2,2)*x - 
   e0*n2*o2*pow(u2,2)*x - f2*j2*r0*pow(u2,2)*x - f0*j2*r2*pow(u2,2)*x + 
   2*n2*r0*r2*pow(u2,2)*x + n0*pow(r2,2)*pow(u2,2)*x + 
   e2*j2*s0*pow(u2,2)*x - m2*r2*s0*pow(u2,2)*x + e0*j2*s2*pow(u2,2)*x - 
   m2*r0*s2*pow(u2,2)*x - c2*d2*k2*o2*v0*x + c2*d2*j2*p2*v0*x + 
   f2*l2*p2*q2*v0*x + d2*m2*p2*q2*v0*x - 2*c2*n2*p2*q2*v0*x - 
   f2*k2*pow(q2,2)*v0*x - 2*d2*l2*p2*r2*v0*x + d2*k2*q2*r2*v0*x + 
   c2*l2*p2*s2*v0*x + c2*k2*q2*s2*v0*x - f2*l2*o2*u2*v0*x - 
   d2*m2*o2*u2*v0*x + 2*c2*n2*o2*u2*v0*x + f2*j2*q2*u2*v0*x + 
   d2*j2*r2*u2*v0*x - 2*n2*q2*r2*u2*v0*x - 2*c2*j2*s2*u2*v0*x + 
   m2*q2*s2*u2*v0*x + l2*r2*s2*u2*v0*x - c2*d0*k2*o2*v2*x - 
   c0*d2*k2*o2*v2*x + c2*d0*j2*p2*v2*x + c0*d2*j2*p2*v2*x + 
   f0*l2*p2*q2*v2*x + d0*m2*p2*q2*v2*x - 2*c2*n0*p2*q2*v2*x - 
   2*c0*n2*p2*q2*v2*x - f0*k2*pow(q2,2)*v2*x - 2*d2*l2*p2*r0*v2*x + 
   d2*k2*q2*r0*v2*x - 2*d0*l2*p2*r2*v2*x + d0*k2*q2*r2*v2*x + 
   c2*l2*p2*s0*v2*x + c2*k2*q2*s0*v2*x + c0*l2*p2*s2*v2*x + 
   c0*k2*q2*s2*v2*x - f2*l2*o2*u0*v2*x - d2*m2*o2*u0*v2*x + 
   2*c2*n2*o2*u0*v2*x + f2*j2*q2*u0*v2*x + d2*j2*r2*u0*v2*x - 
   2*n2*q2*r2*u0*v2*x - 2*c2*j2*s2*u0*v2*x + m2*q2*s2*u0*v2*x + 
   l2*r2*s2*u0*v2*x - f0*l2*o2*u2*v2*x - d0*m2*o2*u2*v2*x + 
   2*c2*n0*o2*u2*v2*x + 2*c0*n2*o2*u2*v2*x + f0*j2*q2*u2*v2*x + 
   d2*j2*r0*u2*v2*x - 2*n2*q2*r0*u2*v2*x + d0*j2*r2*u2*v2*x - 
   2*n0*q2*r2*u2*v2*x - 2*c2*j2*s0*u2*v2*x + m2*q2*s0*u2*v2*x + 
   l2*r2*s0*u2*v2*x - 2*c0*j2*s2*u2*v2*x + l2*r0*s2*u2*v2*x + 
   2*d2*l2*o2*v0*v2*x - 2*d2*j2*q2*v0*v2*x + 2*n2*pow(q2,2)*v0*v2*x - 
   2*l2*q2*s2*v0*v2*x + d0*l2*o2*pow(v2,2)*x - d0*j2*q2*pow(v2,2)*x + 
   n0*pow(q2,2)*pow(v2,2)*x - l2*q2*s0*pow(v2,2)*x + 
   pow(c2,2)*k2*o2*w0*x - pow(c2,2)*j2*p2*w0*x - e2*l2*p2*q2*w0*x + 
   c2*m2*p2*q2*w0*x + e2*k2*pow(q2,2)*w0*x + c2*l2*p2*r2*w0*x - 
   2*c2*k2*q2*r2*w0*x + e2*l2*o2*u2*w0*x - c2*m2*o2*u2*w0*x - 
   e2*j2*q2*u2*w0*x + c2*j2*r2*u2*w0*x + m2*q2*r2*u2*w0*x - 
   l2*pow(r2,2)*u2*w0*x - c2*l2*o2*v2*w0*x + c2*j2*q2*v2*w0*x - 
   m2*pow(q2,2)*v2*w0*x + l2*q2*r2*v2*w0*x + 2*c0*c2*k2*o2*w2*x - 
   2*c0*c2*j2*p2*w2*x - e0*l2*p2*q2*w2*x + c0*m2*p2*q2*w2*x + 
   e0*k2*pow(q2,2)*w2*x + c2*l2*p2*r0*w2*x - 2*c2*k2*q2*r0*w2*x + 
   c0*l2*p2*r2*w2*x - 2*c0*k2*q2*r2*w2*x + e2*l2*o2*u0*w2*x - 
   c2*m2*o2*u0*w2*x - e2*j2*q2*u0*w2*x + c2*j2*r2*u0*w2*x + 
   m2*q2*r2*u0*w2*x - l2*pow(r2,2)*u0*w2*x + e0*l2*o2*u2*w2*x - 
   c0*m2*o2*u2*w2*x - e0*j2*q2*u2*w2*x + c2*j2*r0*u2*w2*x + 
   m2*q2*r0*u2*w2*x + c0*j2*r2*u2*w2*x - 2*l2*r0*r2*u2*w2*x - 
   c2*l2*o2*v0*w2*x + c2*j2*q2*v0*w2*x - m2*pow(q2,2)*v0*w2*x + 
   l2*q2*r2*v0*w2*x - c0*l2*o2*v2*w2*x + c0*j2*q2*v2*w2*x + 
   l2*q2*r0*v2*w2*x + d0*e0*l2*pow(p2,2)*pow(x,2) - 
   c0*f0*l2*pow(p2,2)*pow(x,2) - c0*d0*m2*pow(p2,2)*pow(x,2) + 
   2*c0*c2*n0*pow(p2,2)*pow(x,2) + 
   pow(c0,2)*n2*pow(p2,2)*pow(x,2) - d0*e0*k2*p2*q2*pow(x,2) + 
   c0*f0*k2*p2*q2*pow(x,2) + c2*d0*k2*p2*r0*pow(x,2) + 
   c0*d2*k2*p2*r0*pow(x,2) + c0*d0*k2*p2*r2*pow(x,2) - 
   2*c0*c2*k2*p2*s0*pow(x,2) - pow(c0,2)*k2*p2*s2*pow(x,2) - 
   d0*e0*l2*o2*t2*pow(x,2) + c0*f0*l2*o2*t2*pow(x,2) + 
   c0*d0*m2*o2*t2*pow(x,2) - 2*c0*c2*n0*o2*t2*pow(x,2) - 
   pow(c0,2)*n2*o2*t2*pow(x,2) + d0*e0*j2*q2*t2*pow(x,2) - 
   c0*f0*j2*q2*t2*pow(x,2) - e0*n0*pow(q2,2)*t2*pow(x,2) - 
   c2*d0*j2*r0*t2*pow(x,2) - c0*d2*j2*r0*t2*pow(x,2) - 
   f0*l2*q2*r0*t2*pow(x,2) - d0*m2*q2*r0*t2*pow(x,2) + 
   2*c2*n0*q2*r0*t2*pow(x,2) + 2*c0*n2*q2*r0*t2*pow(x,2) + 
   d2*l2*pow(r0,2)*t2*pow(x,2) - c0*d0*j2*r2*t2*pow(x,2) + 
   2*c0*n0*q2*r2*t2*pow(x,2) + 2*d0*l2*r0*r2*t2*pow(x,2) + 
   2*c0*c2*j2*s0*t2*pow(x,2) + e0*l2*q2*s0*t2*pow(x,2) - 
   c0*m2*q2*s0*t2*pow(x,2) - c2*l2*r0*s0*t2*pow(x,2) - 
   c0*l2*r2*s0*t2*pow(x,2) + pow(c0,2)*j2*s2*t2*pow(x,2) - 
   c0*l2*r0*s2*t2*pow(x,2) + d2*e0*k2*o2*u0*pow(x,2) + 
   d0*e2*k2*o2*u0*pow(x,2) - c2*f0*k2*o2*u0*pow(x,2) - 
   c0*f2*k2*o2*u0*pow(x,2) - d2*e0*j2*p2*u0*pow(x,2) - 
   d0*e2*j2*p2*u0*pow(x,2) + c2*f0*j2*p2*u0*pow(x,2) + 
   c0*f2*j2*p2*u0*pow(x,2) - 2*f0*m2*p2*q2*u0*pow(x,2) + 
   2*e2*n0*p2*q2*u0*pow(x,2) + 2*e0*n2*p2*q2*u0*pow(x,2) + 
   f2*l2*p2*r0*u0*pow(x,2) + d2*m2*p2*r0*u0*pow(x,2) - 
   2*c2*n2*p2*r0*u0*pow(x,2) + f2*k2*q2*r0*u0*pow(x,2) + 
   f0*l2*p2*r2*u0*pow(x,2) + d0*m2*p2*r2*u0*pow(x,2) - 
   2*c2*n0*p2*r2*u0*pow(x,2) - 2*c0*n2*p2*r2*u0*pow(x,2) + 
   f0*k2*q2*r2*u0*pow(x,2) - 2*d2*k2*r0*r2*u0*pow(x,2) - 
   d0*k2*pow(r2,2)*u0*pow(x,2) - e2*l2*p2*s0*u0*pow(x,2) + 
   c2*m2*p2*s0*u0*pow(x,2) - e2*k2*q2*s0*u0*pow(x,2) + 
   c2*k2*r2*s0*u0*pow(x,2) - e0*l2*p2*s2*u0*pow(x,2) + 
   c0*m2*p2*s2*u0*pow(x,2) - e0*k2*q2*s2*u0*pow(x,2) + 
   c2*k2*r0*s2*u0*pow(x,2) + c0*k2*r2*s2*u0*pow(x,2) + 
   f2*m2*o2*pow(u0,2)*pow(x,2) - e2*n2*o2*pow(u0,2)*pow(x,2) - 
   f2*j2*r2*pow(u0,2)*pow(x,2) + n2*pow(r2,2)*pow(u0,2)*pow(x,2) + 
   e2*j2*s2*pow(u0,2)*pow(x,2) - m2*r2*s2*pow(u0,2)*pow(x,2) + 
   d0*e0*k2*o2*u2*pow(x,2) - c0*f0*k2*o2*u2*pow(x,2) - 
   d0*e0*j2*p2*u2*pow(x,2) + c0*f0*j2*p2*u2*pow(x,2) + 
   2*e0*n0*p2*q2*u2*pow(x,2) + f0*l2*p2*r0*u2*pow(x,2) + 
   d0*m2*p2*r0*u2*pow(x,2) - 2*c2*n0*p2*r0*u2*pow(x,2) - 
   2*c0*n2*p2*r0*u2*pow(x,2) + f0*k2*q2*r0*u2*pow(x,2) - 
   d2*k2*pow(r0,2)*u2*pow(x,2) - 2*c0*n0*p2*r2*u2*pow(x,2) - 
   2*d0*k2*r0*r2*u2*pow(x,2) - e0*l2*p2*s0*u2*pow(x,2) + 
   c0*m2*p2*s0*u2*pow(x,2) - e0*k2*q2*s0*u2*pow(x,2) + 
   c2*k2*r0*s0*u2*pow(x,2) + c0*k2*r2*s0*u2*pow(x,2) + 
   c0*k2*r0*s2*u2*pow(x,2) + 2*f0*m2*o2*u0*u2*pow(x,2) - 
   2*e2*n0*o2*u0*u2*pow(x,2) - 2*e0*n2*o2*u0*u2*pow(x,2) - 
   2*f2*j2*r0*u0*u2*pow(x,2) - 2*f0*j2*r2*u0*u2*pow(x,2) + 
   4*n2*r0*r2*u0*u2*pow(x,2) + 2*n0*pow(r2,2)*u0*u2*pow(x,2) + 
   2*e2*j2*s0*u0*u2*pow(x,2) - 2*m2*r2*s0*u0*u2*pow(x,2) + 
   2*e0*j2*s2*u0*u2*pow(x,2) - 2*m2*r0*s2*u0*u2*pow(x,2) - 
   e0*n0*o2*pow(u2,2)*pow(x,2) - f0*j2*r0*pow(u2,2)*pow(x,2) + 
   n2*pow(r0,2)*pow(u2,2)*pow(x,2) + 
   2*n0*r0*r2*pow(u2,2)*pow(x,2) + e0*j2*s0*pow(u2,2)*pow(x,2) - 
   m2*r0*s0*pow(u2,2)*pow(x,2) - c2*d0*k2*o2*v0*pow(x,2) - 
   c0*d2*k2*o2*v0*pow(x,2) + c2*d0*j2*p2*v0*pow(x,2) + 
   c0*d2*j2*p2*v0*pow(x,2) + f0*l2*p2*q2*v0*pow(x,2) + 
   d0*m2*p2*q2*v0*pow(x,2) - 2*c2*n0*p2*q2*v0*pow(x,2) - 
   2*c0*n2*p2*q2*v0*pow(x,2) - f0*k2*pow(q2,2)*v0*pow(x,2) - 
   2*d2*l2*p2*r0*v0*pow(x,2) + d2*k2*q2*r0*v0*pow(x,2) - 
   2*d0*l2*p2*r2*v0*pow(x,2) + d0*k2*q2*r2*v0*pow(x,2) + 
   c2*l2*p2*s0*v0*pow(x,2) + c2*k2*q2*s0*v0*pow(x,2) + 
   c0*l2*p2*s2*v0*pow(x,2) + c0*k2*q2*s2*v0*pow(x,2) - 
   f2*l2*o2*u0*v0*pow(x,2) - d2*m2*o2*u0*v0*pow(x,2) + 
   2*c2*n2*o2*u0*v0*pow(x,2) + f2*j2*q2*u0*v0*pow(x,2) + 
   d2*j2*r2*u0*v0*pow(x,2) - 2*n2*q2*r2*u0*v0*pow(x,2) - 
   2*c2*j2*s2*u0*v0*pow(x,2) + m2*q2*s2*u0*v0*pow(x,2) + 
   l2*r2*s2*u0*v0*pow(x,2) - f0*l2*o2*u2*v0*pow(x,2) - 
   d0*m2*o2*u2*v0*pow(x,2) + 2*c2*n0*o2*u2*v0*pow(x,2) + 
   2*c0*n2*o2*u2*v0*pow(x,2) + f0*j2*q2*u2*v0*pow(x,2) + 
   d2*j2*r0*u2*v0*pow(x,2) - 2*n2*q2*r0*u2*v0*pow(x,2) + 
   d0*j2*r2*u2*v0*pow(x,2) - 2*n0*q2*r2*u2*v0*pow(x,2) - 
   2*c2*j2*s0*u2*v0*pow(x,2) + m2*q2*s0*u2*v0*pow(x,2) + 
   l2*r2*s0*u2*v0*pow(x,2) - 2*c0*j2*s2*u2*v0*pow(x,2) + 
   l2*r0*s2*u2*v0*pow(x,2) + d2*l2*o2*pow(v0,2)*pow(x,2) - 
   d2*j2*q2*pow(v0,2)*pow(x,2) + n2*pow(q2,2)*pow(v0,2)*pow(x,2) - 
   l2*q2*s2*pow(v0,2)*pow(x,2) - c0*d0*k2*o2*v2*pow(x,2) + 
   c0*d0*j2*p2*v2*pow(x,2) - 2*c0*n0*p2*q2*v2*pow(x,2) - 
   2*d0*l2*p2*r0*v2*pow(x,2) + d0*k2*q2*r0*v2*pow(x,2) + 
   c0*l2*p2*s0*v2*pow(x,2) + c0*k2*q2*s0*v2*pow(x,2) - 
   f0*l2*o2*u0*v2*pow(x,2) - d0*m2*o2*u0*v2*pow(x,2) + 
   2*c2*n0*o2*u0*v2*pow(x,2) + 2*c0*n2*o2*u0*v2*pow(x,2) + 
   f0*j2*q2*u0*v2*pow(x,2) + d2*j2*r0*u0*v2*pow(x,2) - 
   2*n2*q2*r0*u0*v2*pow(x,2) + d0*j2*r2*u0*v2*pow(x,2) - 
   2*n0*q2*r2*u0*v2*pow(x,2) - 2*c2*j2*s0*u0*v2*pow(x,2) + 
   m2*q2*s0*u0*v2*pow(x,2) + l2*r2*s0*u0*v2*pow(x,2) - 
   2*c0*j2*s2*u0*v2*pow(x,2) + l2*r0*s2*u0*v2*pow(x,2) + 
   2*c0*n0*o2*u2*v2*pow(x,2) + d0*j2*r0*u2*v2*pow(x,2) - 
   2*n0*q2*r0*u2*v2*pow(x,2) - 2*c0*j2*s0*u2*v2*pow(x,2) + 
   l2*r0*s0*u2*v2*pow(x,2) + 2*d0*l2*o2*v0*v2*pow(x,2) - 
   2*d0*j2*q2*v0*v2*pow(x,2) + 2*n0*pow(q2,2)*v0*v2*pow(x,2) - 
   2*l2*q2*s0*v0*v2*pow(x,2) + 2*c0*c2*k2*o2*w0*pow(x,2) - 
   2*c0*c2*j2*p2*w0*pow(x,2) - e0*l2*p2*q2*w0*pow(x,2) + 
   c0*m2*p2*q2*w0*pow(x,2) + e0*k2*pow(q2,2)*w0*pow(x,2) + 
   c2*l2*p2*r0*w0*pow(x,2) - 2*c2*k2*q2*r0*w0*pow(x,2) + 
   c0*l2*p2*r2*w0*pow(x,2) - 2*c0*k2*q2*r2*w0*pow(x,2) + 
   e2*l2*o2*u0*w0*pow(x,2) - c2*m2*o2*u0*w0*pow(x,2) - 
   e2*j2*q2*u0*w0*pow(x,2) + c2*j2*r2*u0*w0*pow(x,2) + 
   m2*q2*r2*u0*w0*pow(x,2) - l2*pow(r2,2)*u0*w0*pow(x,2) + 
   e0*l2*o2*u2*w0*pow(x,2) - c0*m2*o2*u2*w0*pow(x,2) - 
   e0*j2*q2*u2*w0*pow(x,2) + c2*j2*r0*u2*w0*pow(x,2) + 
   m2*q2*r0*u2*w0*pow(x,2) + c0*j2*r2*u2*w0*pow(x,2) - 
   2*l2*r0*r2*u2*w0*pow(x,2) - c2*l2*o2*v0*w0*pow(x,2) + 
   c2*j2*q2*v0*w0*pow(x,2) - m2*pow(q2,2)*v0*w0*pow(x,2) + 
   l2*q2*r2*v0*w0*pow(x,2) - c0*l2*o2*v2*w0*pow(x,2) + 
   c0*j2*q2*v2*w0*pow(x,2) + l2*q2*r0*v2*w0*pow(x,2) + 
   pow(c0,2)*k2*o2*w2*pow(x,2) - pow(c0,2)*j2*p2*w2*pow(x,2) + 
   c0*l2*p2*r0*w2*pow(x,2) - 2*c0*k2*q2*r0*w2*pow(x,2) + 
   e0*l2*o2*u0*w2*pow(x,2) - c0*m2*o2*u0*w2*pow(x,2) - 
   e0*j2*q2*u0*w2*pow(x,2) + c2*j2*r0*u0*w2*pow(x,2) + 
   m2*q2*r0*u0*w2*pow(x,2) + c0*j2*r2*u0*w2*pow(x,2) - 
   2*l2*r0*r2*u0*w2*pow(x,2) + c0*j2*r0*u2*w2*pow(x,2) - 
   l2*pow(r0,2)*u2*w2*pow(x,2) - c0*l2*o2*v0*w2*pow(x,2) + 
   c0*j2*q2*v0*w2*pow(x,2) + l2*q2*r0*v0*w2*pow(x,2) + 
   pow(c0,2)*n0*pow(p2,2)*pow(x,3) + c0*d0*k2*p2*r0*pow(x,3) - 
   pow(c0,2)*k2*p2*s0*pow(x,3) - pow(c0,2)*n0*o2*t2*pow(x,3) - 
   c0*d0*j2*r0*t2*pow(x,3) + 2*c0*n0*q2*r0*t2*pow(x,3) + 
   d0*l2*pow(r0,2)*t2*pow(x,3) + pow(c0,2)*j2*s0*t2*pow(x,3) - 
   c0*l2*r0*s0*t2*pow(x,3) + d0*e0*k2*o2*u0*pow(x,3) - 
   c0*f0*k2*o2*u0*pow(x,3) - d0*e0*j2*p2*u0*pow(x,3) + 
   c0*f0*j2*p2*u0*pow(x,3) + 2*e0*n0*p2*q2*u0*pow(x,3) + 
   f0*l2*p2*r0*u0*pow(x,3) + d0*m2*p2*r0*u0*pow(x,3) - 
   2*c2*n0*p2*r0*u0*pow(x,3) - 2*c0*n2*p2*r0*u0*pow(x,3) + 
   f0*k2*q2*r0*u0*pow(x,3) - d2*k2*pow(r0,2)*u0*pow(x,3) - 
   2*c0*n0*p2*r2*u0*pow(x,3) - 2*d0*k2*r0*r2*u0*pow(x,3) - 
   e0*l2*p2*s0*u0*pow(x,3) + c0*m2*p2*s0*u0*pow(x,3) - 
   e0*k2*q2*s0*u0*pow(x,3) + c2*k2*r0*s0*u0*pow(x,3) + 
   c0*k2*r2*s0*u0*pow(x,3) + c0*k2*r0*s2*u0*pow(x,3) + 
   f0*m2*o2*pow(u0,2)*pow(x,3) - e2*n0*o2*pow(u0,2)*pow(x,3) - 
   e0*n2*o2*pow(u0,2)*pow(x,3) - f2*j2*r0*pow(u0,2)*pow(x,3) - 
   f0*j2*r2*pow(u0,2)*pow(x,3) + 2*n2*r0*r2*pow(u0,2)*pow(x,3) + 
   n0*pow(r2,2)*pow(u0,2)*pow(x,3) + e2*j2*s0*pow(u0,2)*pow(x,3) - 
   m2*r2*s0*pow(u0,2)*pow(x,3) + e0*j2*s2*pow(u0,2)*pow(x,3) - 
   m2*r0*s2*pow(u0,2)*pow(x,3) - 2*c0*n0*p2*r0*u2*pow(x,3) - 
   d0*k2*pow(r0,2)*u2*pow(x,3) + c0*k2*r0*s0*u2*pow(x,3) - 
   2*e0*n0*o2*u0*u2*pow(x,3) - 2*f0*j2*r0*u0*u2*pow(x,3) + 
   2*n2*pow(r0,2)*u0*u2*pow(x,3) + 4*n0*r0*r2*u0*u2*pow(x,3) + 
   2*e0*j2*s0*u0*u2*pow(x,3) - 2*m2*r0*s0*u0*u2*pow(x,3) + 
   n0*pow(r0,2)*pow(u2,2)*pow(x,3) - c0*d0*k2*o2*v0*pow(x,3) + 
   c0*d0*j2*p2*v0*pow(x,3) - 2*c0*n0*p2*q2*v0*pow(x,3) - 
   2*d0*l2*p2*r0*v0*pow(x,3) + d0*k2*q2*r0*v0*pow(x,3) + 
   c0*l2*p2*s0*v0*pow(x,3) + c0*k2*q2*s0*v0*pow(x,3) - 
   f0*l2*o2*u0*v0*pow(x,3) - d0*m2*o2*u0*v0*pow(x,3) + 
   2*c2*n0*o2*u0*v0*pow(x,3) + 2*c0*n2*o2*u0*v0*pow(x,3) + 
   f0*j2*q2*u0*v0*pow(x,3) + d2*j2*r0*u0*v0*pow(x,3) - 
   2*n2*q2*r0*u0*v0*pow(x,3) + d0*j2*r2*u0*v0*pow(x,3) - 
   2*n0*q2*r2*u0*v0*pow(x,3) - 2*c2*j2*s0*u0*v0*pow(x,3) + 
   m2*q2*s0*u0*v0*pow(x,3) + l2*r2*s0*u0*v0*pow(x,3) - 
   2*c0*j2*s2*u0*v0*pow(x,3) + l2*r0*s2*u0*v0*pow(x,3) + 
   2*c0*n0*o2*u2*v0*pow(x,3) + d0*j2*r0*u2*v0*pow(x,3) - 
   2*n0*q2*r0*u2*v0*pow(x,3) - 2*c0*j2*s0*u2*v0*pow(x,3) + 
   l2*r0*s0*u2*v0*pow(x,3) + d0*l2*o2*pow(v0,2)*pow(x,3) - 
   d0*j2*q2*pow(v0,2)*pow(x,3) + n0*pow(q2,2)*pow(v0,2)*pow(x,3) - 
   l2*q2*s0*pow(v0,2)*pow(x,3) + 2*c0*n0*o2*u0*v2*pow(x,3) + 
   d0*j2*r0*u0*v2*pow(x,3) - 2*n0*q2*r0*u0*v2*pow(x,3) - 
   2*c0*j2*s0*u0*v2*pow(x,3) + l2*r0*s0*u0*v2*pow(x,3) + 
   pow(c0,2)*k2*o2*w0*pow(x,3) - pow(c0,2)*j2*p2*w0*pow(x,3) + 
   c0*l2*p2*r0*w0*pow(x,3) - 2*c0*k2*q2*r0*w0*pow(x,3) + 
   e0*l2*o2*u0*w0*pow(x,3) - c0*m2*o2*u0*w0*pow(x,3) - 
   e0*j2*q2*u0*w0*pow(x,3) + c2*j2*r0*u0*w0*pow(x,3) + 
   m2*q2*r0*u0*w0*pow(x,3) + c0*j2*r2*u0*w0*pow(x,3) - 
   2*l2*r0*r2*u0*w0*pow(x,3) + c0*j2*r0*u2*w0*pow(x,3) - 
   l2*pow(r0,2)*u2*w0*pow(x,3) - c0*l2*o2*v0*w0*pow(x,3) + 
   c0*j2*q2*v0*w0*pow(x,3) + l2*q2*r0*v0*w0*pow(x,3) + 
   c0*j2*r0*u0*w2*pow(x,3) - l2*pow(r0,2)*u0*w2*pow(x,3) - 
   2*c0*n0*p2*r0*u0*pow(x,4) - d0*k2*pow(r0,2)*u0*pow(x,4) + 
   c0*k2*r0*s0*u0*pow(x,4) - e0*n0*o2*pow(u0,2)*pow(x,4) - 
   f0*j2*r0*pow(u0,2)*pow(x,4) + n2*pow(r0,2)*pow(u0,2)*pow(x,4) + 
   2*n0*r0*r2*pow(u0,2)*pow(x,4) + e0*j2*s0*pow(u0,2)*pow(x,4) - 
   m2*r0*s0*pow(u0,2)*pow(x,4) + 2*n0*pow(r0,2)*u0*u2*pow(x,4) + 
   2*c0*n0*o2*u0*v0*pow(x,4) + d0*j2*r0*u0*v0*pow(x,4) - 
   2*n0*q2*r0*u0*v0*pow(x,4) - 2*c0*j2*s0*u0*v0*pow(x,4) + 
   l2*r0*s0*u0*v0*pow(x,4) + c0*j2*r0*u0*w0*pow(x,4) - 
   l2*pow(r0,2)*u0*w0*pow(x,4) + n0*pow(r0,2)*pow(u0,2)*pow(x,5) + 
   d2*e1*l2*pow(p2,2)*y + d1*e2*l2*pow(p2,2)*y - c2*f1*l2*pow(p2,2)*y - 
   c1*f2*l2*pow(p2,2)*y - c2*d1*m2*pow(p2,2)*y - c1*d2*m2*pow(p2,2)*y + 
   pow(c2,2)*n1*pow(p2,2)*y + 2*c1*c2*n2*pow(p2,2)*y - 
   d2*e1*k2*p2*q2*y - d1*e2*k2*p2*q2*y + c2*f1*k2*p2*q2*y + 
   c1*f2*k2*p2*q2*y + c2*d2*k2*p2*r1*y + c2*d1*k2*p2*r2*y + 
   c1*d2*k2*p2*r2*y - pow(c2,2)*k2*p2*s1*y - 2*c1*c2*k2*p2*s2*y - 
   d2*e1*l2*o2*t2*y - d1*e2*l2*o2*t2*y + c2*f1*l2*o2*t2*y + 
   c1*f2*l2*o2*t2*y + c2*d1*m2*o2*t2*y + c1*d2*m2*o2*t2*y - 
   pow(c2,2)*n1*o2*t2*y - 2*c1*c2*n2*o2*t2*y + d2*e1*j2*q2*t2*y + 
   d1*e2*j2*q2*t2*y - c2*f1*j2*q2*t2*y - c1*f2*j2*q2*t2*y + 
   f1*m2*pow(q2,2)*t2*y - e2*n1*pow(q2,2)*t2*y - e1*n2*pow(q2,2)*t2*y - 
   c2*d2*j2*r1*t2*y - f2*l2*q2*r1*t2*y - d2*m2*q2*r1*t2*y + 
   2*c2*n2*q2*r1*t2*y - c2*d1*j2*r2*t2*y - c1*d2*j2*r2*t2*y - 
   f1*l2*q2*r2*t2*y - d1*m2*q2*r2*t2*y + 2*c2*n1*q2*r2*t2*y + 
   2*c1*n2*q2*r2*t2*y + 2*d2*l2*r1*r2*t2*y + d1*l2*pow(r2,2)*t2*y + 
   pow(c2,2)*j2*s1*t2*y + e2*l2*q2*s1*t2*y - c2*m2*q2*s1*t2*y - 
   c2*l2*r2*s1*t2*y + 2*c1*c2*j2*s2*t2*y + e1*l2*q2*s2*t2*y - 
   c1*m2*q2*s2*t2*y - c2*l2*r1*s2*t2*y - c1*l2*r2*s2*t2*y + 
   d2*e2*k2*o2*u1*y - c2*f2*k2*o2*u1*y - d2*e2*j2*p2*u1*y + 
   c2*f2*j2*p2*u1*y - 2*f2*m2*p2*q2*u1*y + 2*e2*n2*p2*q2*u1*y + 
   f2*l2*p2*r2*u1*y + d2*m2*p2*r2*u1*y - 2*c2*n2*p2*r2*u1*y + 
   f2*k2*q2*r2*u1*y - d2*k2*pow(r2,2)*u1*y - e2*l2*p2*s2*u1*y + 
   c2*m2*p2*s2*u1*y - e2*k2*q2*s2*u1*y + c2*k2*r2*s2*u1*y + 
   d2*e1*k2*o2*u2*y + d1*e2*k2*o2*u2*y - c2*f1*k2*o2*u2*y - 
   c1*f2*k2*o2*u2*y - d2*e1*j2*p2*u2*y - d1*e2*j2*p2*u2*y + 
   c2*f1*j2*p2*u2*y + c1*f2*j2*p2*u2*y - 2*f1*m2*p2*q2*u2*y + 
   2*e2*n1*p2*q2*u2*y + 2*e1*n2*p2*q2*u2*y + f2*l2*p2*r1*u2*y + 
   d2*m2*p2*r1*u2*y - 2*c2*n2*p2*r1*u2*y + f2*k2*q2*r1*u2*y + 
   f1*l2*p2*r2*u2*y + d1*m2*p2*r2*u2*y - 2*c2*n1*p2*r2*u2*y - 
   2*c1*n2*p2*r2*u2*y + f1*k2*q2*r2*u2*y - 2*d2*k2*r1*r2*u2*y - 
   d1*k2*pow(r2,2)*u2*y - e2*l2*p2*s1*u2*y + c2*m2*p2*s1*u2*y - 
   e2*k2*q2*s1*u2*y + c2*k2*r2*s1*u2*y - e1*l2*p2*s2*u2*y + 
   c1*m2*p2*s2*u2*y - e1*k2*q2*s2*u2*y + c2*k2*r1*s2*u2*y + 
   c1*k2*r2*s2*u2*y + 2*f2*m2*o2*u1*u2*y - 2*e2*n2*o2*u1*u2*y - 
   2*f2*j2*r2*u1*u2*y + 2*n2*pow(r2,2)*u1*u2*y + 2*e2*j2*s2*u1*u2*y - 
   2*m2*r2*s2*u1*u2*y + f1*m2*o2*pow(u2,2)*y - e2*n1*o2*pow(u2,2)*y - 
   e1*n2*o2*pow(u2,2)*y - f2*j2*r1*pow(u2,2)*y - f1*j2*r2*pow(u2,2)*y + 
   2*n2*r1*r2*pow(u2,2)*y + n1*pow(r2,2)*pow(u2,2)*y + 
   e2*j2*s1*pow(u2,2)*y - m2*r2*s1*pow(u2,2)*y + e1*j2*s2*pow(u2,2)*y - 
   m2*r1*s2*pow(u2,2)*y - c2*d2*k2*o2*v1*y + c2*d2*j2*p2*v1*y + 
   f2*l2*p2*q2*v1*y + d2*m2*p2*q2*v1*y - 2*c2*n2*p2*q2*v1*y - 
   f2*k2*pow(q2,2)*v1*y - 2*d2*l2*p2*r2*v1*y + d2*k2*q2*r2*v1*y + 
   c2*l2*p2*s2*v1*y + c2*k2*q2*s2*v1*y - f2*l2*o2*u2*v1*y - 
   d2*m2*o2*u2*v1*y + 2*c2*n2*o2*u2*v1*y + f2*j2*q2*u2*v1*y + 
   d2*j2*r2*u2*v1*y - 2*n2*q2*r2*u2*v1*y - 2*c2*j2*s2*u2*v1*y + 
   m2*q2*s2*u2*v1*y + l2*r2*s2*u2*v1*y - c2*d1*k2*o2*v2*y - 
   c1*d2*k2*o2*v2*y + c2*d1*j2*p2*v2*y + c1*d2*j2*p2*v2*y + 
   f1*l2*p2*q2*v2*y + d1*m2*p2*q2*v2*y - 2*c2*n1*p2*q2*v2*y - 
   2*c1*n2*p2*q2*v2*y - f1*k2*pow(q2,2)*v2*y - 2*d2*l2*p2*r1*v2*y + 
   d2*k2*q2*r1*v2*y - 2*d1*l2*p2*r2*v2*y + d1*k2*q2*r2*v2*y + 
   c2*l2*p2*s1*v2*y + c2*k2*q2*s1*v2*y + c1*l2*p2*s2*v2*y + 
   c1*k2*q2*s2*v2*y - f2*l2*o2*u1*v2*y - d2*m2*o2*u1*v2*y + 
   2*c2*n2*o2*u1*v2*y + f2*j2*q2*u1*v2*y + d2*j2*r2*u1*v2*y - 
   2*n2*q2*r2*u1*v2*y - 2*c2*j2*s2*u1*v2*y + m2*q2*s2*u1*v2*y + 
   l2*r2*s2*u1*v2*y - f1*l2*o2*u2*v2*y - d1*m2*o2*u2*v2*y + 
   2*c2*n1*o2*u2*v2*y + 2*c1*n2*o2*u2*v2*y + f1*j2*q2*u2*v2*y + 
   d2*j2*r1*u2*v2*y - 2*n2*q2*r1*u2*v2*y + d1*j2*r2*u2*v2*y - 
   2*n1*q2*r2*u2*v2*y - 2*c2*j2*s1*u2*v2*y + m2*q2*s1*u2*v2*y + 
   l2*r2*s1*u2*v2*y - 2*c1*j2*s2*u2*v2*y + l2*r1*s2*u2*v2*y + 
   2*d2*l2*o2*v1*v2*y - 2*d2*j2*q2*v1*v2*y + 2*n2*pow(q2,2)*v1*v2*y - 
   2*l2*q2*s2*v1*v2*y + d1*l2*o2*pow(v2,2)*y - d1*j2*q2*pow(v2,2)*y + 
   n1*pow(q2,2)*pow(v2,2)*y - l2*q2*s1*pow(v2,2)*y + 
   pow(c2,2)*k2*o2*w1*y - pow(c2,2)*j2*p2*w1*y - e2*l2*p2*q2*w1*y + 
   c2*m2*p2*q2*w1*y + e2*k2*pow(q2,2)*w1*y + c2*l2*p2*r2*w1*y - 
   2*c2*k2*q2*r2*w1*y + e2*l2*o2*u2*w1*y - c2*m2*o2*u2*w1*y - 
   e2*j2*q2*u2*w1*y + c2*j2*r2*u2*w1*y + m2*q2*r2*u2*w1*y - 
   l2*pow(r2,2)*u2*w1*y - c2*l2*o2*v2*w1*y + c2*j2*q2*v2*w1*y - 
   m2*pow(q2,2)*v2*w1*y + l2*q2*r2*v2*w1*y + 2*c1*c2*k2*o2*w2*y - 
   2*c1*c2*j2*p2*w2*y - e1*l2*p2*q2*w2*y + c1*m2*p2*q2*w2*y + 
   e1*k2*pow(q2,2)*w2*y + c2*l2*p2*r1*w2*y - 2*c2*k2*q2*r1*w2*y + 
   c1*l2*p2*r2*w2*y - 2*c1*k2*q2*r2*w2*y + e2*l2*o2*u1*w2*y - 
   c2*m2*o2*u1*w2*y - e2*j2*q2*u1*w2*y + c2*j2*r2*u1*w2*y + 
   m2*q2*r2*u1*w2*y - l2*pow(r2,2)*u1*w2*y + e1*l2*o2*u2*w2*y - 
   c1*m2*o2*u2*w2*y - e1*j2*q2*u2*w2*y + c2*j2*r1*u2*w2*y + 
   m2*q2*r1*u2*w2*y + c1*j2*r2*u2*w2*y - 2*l2*r1*r2*u2*w2*y - 
   c2*l2*o2*v1*w2*y + c2*j2*q2*v1*w2*y - m2*pow(q2,2)*v1*w2*y + 
   l2*q2*r2*v1*w2*y - c1*l2*o2*v2*w2*y + c1*j2*q2*v2*w2*y + 
   l2*q2*r1*v2*w2*y + d1*e0*l2*pow(p2,2)*x*y + d0*e1*l2*pow(p2,2)*x*y - 
   c1*f0*l2*pow(p2,2)*x*y - c0*f1*l2*pow(p2,2)*x*y - 
   c1*d0*m2*pow(p2,2)*x*y - c0*d1*m2*pow(p2,2)*x*y + 
   2*c1*c2*n0*pow(p2,2)*x*y + 2*c0*c2*n1*pow(p2,2)*x*y + 
   2*c0*c1*n2*pow(p2,2)*x*y - d1*e0*k2*p2*q2*x*y - d0*e1*k2*p2*q2*x*y + 
   c1*f0*k2*p2*q2*x*y + c0*f1*k2*p2*q2*x*y + c2*d1*k2*p2*r0*x*y + 
   c1*d2*k2*p2*r0*x*y + c2*d0*k2*p2*r1*x*y + c0*d2*k2*p2*r1*x*y + 
   c1*d0*k2*p2*r2*x*y + c0*d1*k2*p2*r2*x*y - 2*c1*c2*k2*p2*s0*x*y - 
   2*c0*c2*k2*p2*s1*x*y - 2*c0*c1*k2*p2*s2*x*y - d1*e0*l2*o2*t2*x*y - 
   d0*e1*l2*o2*t2*x*y + c1*f0*l2*o2*t2*x*y + c0*f1*l2*o2*t2*x*y + 
   c1*d0*m2*o2*t2*x*y + c0*d1*m2*o2*t2*x*y - 2*c1*c2*n0*o2*t2*x*y - 
   2*c0*c2*n1*o2*t2*x*y - 2*c0*c1*n2*o2*t2*x*y + d1*e0*j2*q2*t2*x*y + 
   d0*e1*j2*q2*t2*x*y - c1*f0*j2*q2*t2*x*y - c0*f1*j2*q2*t2*x*y - 
   e1*n0*pow(q2,2)*t2*x*y - e0*n1*pow(q2,2)*t2*x*y - c2*d1*j2*r0*t2*x*y - 
   c1*d2*j2*r0*t2*x*y - f1*l2*q2*r0*t2*x*y - d1*m2*q2*r0*t2*x*y + 
   2*c2*n1*q2*r0*t2*x*y + 2*c1*n2*q2*r0*t2*x*y - c2*d0*j2*r1*t2*x*y - 
   c0*d2*j2*r1*t2*x*y - f0*l2*q2*r1*t2*x*y - d0*m2*q2*r1*t2*x*y + 
   2*c2*n0*q2*r1*t2*x*y + 2*c0*n2*q2*r1*t2*x*y + 2*d2*l2*r0*r1*t2*x*y - 
   c1*d0*j2*r2*t2*x*y - c0*d1*j2*r2*t2*x*y + 2*c1*n0*q2*r2*t2*x*y + 
   2*c0*n1*q2*r2*t2*x*y + 2*d1*l2*r0*r2*t2*x*y + 2*d0*l2*r1*r2*t2*x*y + 
   2*c1*c2*j2*s0*t2*x*y + e1*l2*q2*s0*t2*x*y - c1*m2*q2*s0*t2*x*y - 
   c2*l2*r1*s0*t2*x*y - c1*l2*r2*s0*t2*x*y + 2*c0*c2*j2*s1*t2*x*y + 
   e0*l2*q2*s1*t2*x*y - c0*m2*q2*s1*t2*x*y - c2*l2*r0*s1*t2*x*y - 
   c0*l2*r2*s1*t2*x*y + 2*c0*c1*j2*s2*t2*x*y - c1*l2*r0*s2*t2*x*y - 
   c0*l2*r1*s2*t2*x*y + d2*e1*k2*o2*u0*x*y + d1*e2*k2*o2*u0*x*y - 
   c2*f1*k2*o2*u0*x*y - c1*f2*k2*o2*u0*x*y - d2*e1*j2*p2*u0*x*y - 
   d1*e2*j2*p2*u0*x*y + c2*f1*j2*p2*u0*x*y + c1*f2*j2*p2*u0*x*y - 
   2*f1*m2*p2*q2*u0*x*y + 2*e2*n1*p2*q2*u0*x*y + 2*e1*n2*p2*q2*u0*x*y + 
   f2*l2*p2*r1*u0*x*y + d2*m2*p2*r1*u0*x*y - 2*c2*n2*p2*r1*u0*x*y + 
   f2*k2*q2*r1*u0*x*y + f1*l2*p2*r2*u0*x*y + d1*m2*p2*r2*u0*x*y - 
   2*c2*n1*p2*r2*u0*x*y - 2*c1*n2*p2*r2*u0*x*y + f1*k2*q2*r2*u0*x*y - 
   2*d2*k2*r1*r2*u0*x*y - d1*k2*pow(r2,2)*u0*x*y - e2*l2*p2*s1*u0*x*y + 
   c2*m2*p2*s1*u0*x*y - e2*k2*q2*s1*u0*x*y + c2*k2*r2*s1*u0*x*y - 
   e1*l2*p2*s2*u0*x*y + c1*m2*p2*s2*u0*x*y - e1*k2*q2*s2*u0*x*y + 
   c2*k2*r1*s2*u0*x*y + c1*k2*r2*s2*u0*x*y + d2*e0*k2*o2*u1*x*y + 
   d0*e2*k2*o2*u1*x*y - c2*f0*k2*o2*u1*x*y - c0*f2*k2*o2*u1*x*y - 
   d2*e0*j2*p2*u1*x*y - d0*e2*j2*p2*u1*x*y + c2*f0*j2*p2*u1*x*y + 
   c0*f2*j2*p2*u1*x*y - 2*f0*m2*p2*q2*u1*x*y + 2*e2*n0*p2*q2*u1*x*y + 
   2*e0*n2*p2*q2*u1*x*y + f2*l2*p2*r0*u1*x*y + d2*m2*p2*r0*u1*x*y - 
   2*c2*n2*p2*r0*u1*x*y + f2*k2*q2*r0*u1*x*y + f0*l2*p2*r2*u1*x*y + 
   d0*m2*p2*r2*u1*x*y - 2*c2*n0*p2*r2*u1*x*y - 2*c0*n2*p2*r2*u1*x*y + 
   f0*k2*q2*r2*u1*x*y - 2*d2*k2*r0*r2*u1*x*y - d0*k2*pow(r2,2)*u1*x*y - 
   e2*l2*p2*s0*u1*x*y + c2*m2*p2*s0*u1*x*y - e2*k2*q2*s0*u1*x*y + 
   c2*k2*r2*s0*u1*x*y - e0*l2*p2*s2*u1*x*y + c0*m2*p2*s2*u1*x*y - 
   e0*k2*q2*s2*u1*x*y + c2*k2*r0*s2*u1*x*y + c0*k2*r2*s2*u1*x*y + 
   2*f2*m2*o2*u0*u1*x*y - 2*e2*n2*o2*u0*u1*x*y - 2*f2*j2*r2*u0*u1*x*y + 
   2*n2*pow(r2,2)*u0*u1*x*y + 2*e2*j2*s2*u0*u1*x*y - 2*m2*r2*s2*u0*u1*x*y + 
   d1*e0*k2*o2*u2*x*y + d0*e1*k2*o2*u2*x*y - c1*f0*k2*o2*u2*x*y - 
   c0*f1*k2*o2*u2*x*y - d1*e0*j2*p2*u2*x*y - d0*e1*j2*p2*u2*x*y + 
   c1*f0*j2*p2*u2*x*y + c0*f1*j2*p2*u2*x*y + 2*e1*n0*p2*q2*u2*x*y + 
   2*e0*n1*p2*q2*u2*x*y + f1*l2*p2*r0*u2*x*y + d1*m2*p2*r0*u2*x*y - 
   2*c2*n1*p2*r0*u2*x*y - 2*c1*n2*p2*r0*u2*x*y + f1*k2*q2*r0*u2*x*y + 
   f0*l2*p2*r1*u2*x*y + d0*m2*p2*r1*u2*x*y - 2*c2*n0*p2*r1*u2*x*y - 
   2*c0*n2*p2*r1*u2*x*y + f0*k2*q2*r1*u2*x*y - 2*d2*k2*r0*r1*u2*x*y - 
   2*c1*n0*p2*r2*u2*x*y - 2*c0*n1*p2*r2*u2*x*y - 2*d1*k2*r0*r2*u2*x*y - 
   2*d0*k2*r1*r2*u2*x*y - e1*l2*p2*s0*u2*x*y + c1*m2*p2*s0*u2*x*y - 
   e1*k2*q2*s0*u2*x*y + c2*k2*r1*s0*u2*x*y + c1*k2*r2*s0*u2*x*y - 
   e0*l2*p2*s1*u2*x*y + c0*m2*p2*s1*u2*x*y - e0*k2*q2*s1*u2*x*y + 
   c2*k2*r0*s1*u2*x*y + c0*k2*r2*s1*u2*x*y + c1*k2*r0*s2*u2*x*y + 
   c0*k2*r1*s2*u2*x*y + 2*f1*m2*o2*u0*u2*x*y - 2*e2*n1*o2*u0*u2*x*y - 
   2*e1*n2*o2*u0*u2*x*y - 2*f2*j2*r1*u0*u2*x*y - 2*f1*j2*r2*u0*u2*x*y + 
   4*n2*r1*r2*u0*u2*x*y + 2*n1*pow(r2,2)*u0*u2*x*y + 2*e2*j2*s1*u0*u2*x*y - 
   2*m2*r2*s1*u0*u2*x*y + 2*e1*j2*s2*u0*u2*x*y - 2*m2*r1*s2*u0*u2*x*y + 
   2*f0*m2*o2*u1*u2*x*y - 2*e2*n0*o2*u1*u2*x*y - 2*e0*n2*o2*u1*u2*x*y - 
   2*f2*j2*r0*u1*u2*x*y - 2*f0*j2*r2*u1*u2*x*y + 4*n2*r0*r2*u1*u2*x*y + 
   2*n0*pow(r2,2)*u1*u2*x*y + 2*e2*j2*s0*u1*u2*x*y - 2*m2*r2*s0*u1*u2*x*y + 
   2*e0*j2*s2*u1*u2*x*y - 2*m2*r0*s2*u1*u2*x*y - e1*n0*o2*pow(u2,2)*x*y - 
   e0*n1*o2*pow(u2,2)*x*y - f1*j2*r0*pow(u2,2)*x*y - 
   f0*j2*r1*pow(u2,2)*x*y + 2*n2*r0*r1*pow(u2,2)*x*y + 
   2*n1*r0*r2*pow(u2,2)*x*y + 2*n0*r1*r2*pow(u2,2)*x*y + 
   e1*j2*s0*pow(u2,2)*x*y - m2*r1*s0*pow(u2,2)*x*y + 
   e0*j2*s1*pow(u2,2)*x*y - m2*r0*s1*pow(u2,2)*x*y - c2*d1*k2*o2*v0*x*y - 
   c1*d2*k2*o2*v0*x*y + c2*d1*j2*p2*v0*x*y + c1*d2*j2*p2*v0*x*y + 
   f1*l2*p2*q2*v0*x*y + d1*m2*p2*q2*v0*x*y - 2*c2*n1*p2*q2*v0*x*y - 
   2*c1*n2*p2*q2*v0*x*y - f1*k2*pow(q2,2)*v0*x*y - 2*d2*l2*p2*r1*v0*x*y + 
   d2*k2*q2*r1*v0*x*y - 2*d1*l2*p2*r2*v0*x*y + d1*k2*q2*r2*v0*x*y + 
   c2*l2*p2*s1*v0*x*y + c2*k2*q2*s1*v0*x*y + c1*l2*p2*s2*v0*x*y + 
   c1*k2*q2*s2*v0*x*y - f2*l2*o2*u1*v0*x*y - d2*m2*o2*u1*v0*x*y + 
   2*c2*n2*o2*u1*v0*x*y + f2*j2*q2*u1*v0*x*y + d2*j2*r2*u1*v0*x*y - 
   2*n2*q2*r2*u1*v0*x*y - 2*c2*j2*s2*u1*v0*x*y + m2*q2*s2*u1*v0*x*y + 
   l2*r2*s2*u1*v0*x*y - f1*l2*o2*u2*v0*x*y - d1*m2*o2*u2*v0*x*y + 
   2*c2*n1*o2*u2*v0*x*y + 2*c1*n2*o2*u2*v0*x*y + f1*j2*q2*u2*v0*x*y + 
   d2*j2*r1*u2*v0*x*y - 2*n2*q2*r1*u2*v0*x*y + d1*j2*r2*u2*v0*x*y - 
   2*n1*q2*r2*u2*v0*x*y - 2*c2*j2*s1*u2*v0*x*y + m2*q2*s1*u2*v0*x*y + 
   l2*r2*s1*u2*v0*x*y - 2*c1*j2*s2*u2*v0*x*y + l2*r1*s2*u2*v0*x*y - 
   c2*d0*k2*o2*v1*x*y - c0*d2*k2*o2*v1*x*y + c2*d0*j2*p2*v1*x*y + 
   c0*d2*j2*p2*v1*x*y + f0*l2*p2*q2*v1*x*y + d0*m2*p2*q2*v1*x*y - 
   2*c2*n0*p2*q2*v1*x*y - 2*c0*n2*p2*q2*v1*x*y - f0*k2*pow(q2,2)*v1*x*y - 
   2*d2*l2*p2*r0*v1*x*y + d2*k2*q2*r0*v1*x*y - 2*d0*l2*p2*r2*v1*x*y + 
   d0*k2*q2*r2*v1*x*y + c2*l2*p2*s0*v1*x*y + c2*k2*q2*s0*v1*x*y + 
   c0*l2*p2*s2*v1*x*y + c0*k2*q2*s2*v1*x*y - f2*l2*o2*u0*v1*x*y - 
   d2*m2*o2*u0*v1*x*y + 2*c2*n2*o2*u0*v1*x*y + f2*j2*q2*u0*v1*x*y + 
   d2*j2*r2*u0*v1*x*y - 2*n2*q2*r2*u0*v1*x*y - 2*c2*j2*s2*u0*v1*x*y + 
   m2*q2*s2*u0*v1*x*y + l2*r2*s2*u0*v1*x*y - f0*l2*o2*u2*v1*x*y - 
   d0*m2*o2*u2*v1*x*y + 2*c2*n0*o2*u2*v1*x*y + 2*c0*n2*o2*u2*v1*x*y + 
   f0*j2*q2*u2*v1*x*y + d2*j2*r0*u2*v1*x*y - 2*n2*q2*r0*u2*v1*x*y + 
   d0*j2*r2*u2*v1*x*y - 2*n0*q2*r2*u2*v1*x*y - 2*c2*j2*s0*u2*v1*x*y + 
   m2*q2*s0*u2*v1*x*y + l2*r2*s0*u2*v1*x*y - 2*c0*j2*s2*u2*v1*x*y + 
   l2*r0*s2*u2*v1*x*y + 2*d2*l2*o2*v0*v1*x*y - 2*d2*j2*q2*v0*v1*x*y + 
   2*n2*pow(q2,2)*v0*v1*x*y - 2*l2*q2*s2*v0*v1*x*y - c1*d0*k2*o2*v2*x*y - 
   c0*d1*k2*o2*v2*x*y + c1*d0*j2*p2*v2*x*y + c0*d1*j2*p2*v2*x*y - 
   2*c1*n0*p2*q2*v2*x*y - 2*c0*n1*p2*q2*v2*x*y - 2*d1*l2*p2*r0*v2*x*y + 
   d1*k2*q2*r0*v2*x*y - 2*d0*l2*p2*r1*v2*x*y + d0*k2*q2*r1*v2*x*y + 
   c1*l2*p2*s0*v2*x*y + c1*k2*q2*s0*v2*x*y + c0*l2*p2*s1*v2*x*y + 
   c0*k2*q2*s1*v2*x*y - f1*l2*o2*u0*v2*x*y - d1*m2*o2*u0*v2*x*y + 
   2*c2*n1*o2*u0*v2*x*y + 2*c1*n2*o2*u0*v2*x*y + f1*j2*q2*u0*v2*x*y + 
   d2*j2*r1*u0*v2*x*y - 2*n2*q2*r1*u0*v2*x*y + d1*j2*r2*u0*v2*x*y - 
   2*n1*q2*r2*u0*v2*x*y - 2*c2*j2*s1*u0*v2*x*y + m2*q2*s1*u0*v2*x*y + 
   l2*r2*s1*u0*v2*x*y - 2*c1*j2*s2*u0*v2*x*y + l2*r1*s2*u0*v2*x*y - 
   f0*l2*o2*u1*v2*x*y - d0*m2*o2*u1*v2*x*y + 2*c2*n0*o2*u1*v2*x*y + 
   2*c0*n2*o2*u1*v2*x*y + f0*j2*q2*u1*v2*x*y + d2*j2*r0*u1*v2*x*y - 
   2*n2*q2*r0*u1*v2*x*y + d0*j2*r2*u1*v2*x*y - 2*n0*q2*r2*u1*v2*x*y - 
   2*c2*j2*s0*u1*v2*x*y + m2*q2*s0*u1*v2*x*y + l2*r2*s0*u1*v2*x*y - 
   2*c0*j2*s2*u1*v2*x*y + l2*r0*s2*u1*v2*x*y + 2*c1*n0*o2*u2*v2*x*y + 
   2*c0*n1*o2*u2*v2*x*y + d1*j2*r0*u2*v2*x*y - 2*n1*q2*r0*u2*v2*x*y + 
   d0*j2*r1*u2*v2*x*y - 2*n0*q2*r1*u2*v2*x*y - 2*c1*j2*s0*u2*v2*x*y + 
   l2*r1*s0*u2*v2*x*y - 2*c0*j2*s1*u2*v2*x*y + l2*r0*s1*u2*v2*x*y + 
   2*d1*l2*o2*v0*v2*x*y - 2*d1*j2*q2*v0*v2*x*y + 2*n1*pow(q2,2)*v0*v2*x*y - 
   2*l2*q2*s1*v0*v2*x*y + 2*d0*l2*o2*v1*v2*x*y - 2*d0*j2*q2*v1*v2*x*y + 
   2*n0*pow(q2,2)*v1*v2*x*y - 2*l2*q2*s0*v1*v2*x*y + 2*c1*c2*k2*o2*w0*x*y - 
   2*c1*c2*j2*p2*w0*x*y - e1*l2*p2*q2*w0*x*y + c1*m2*p2*q2*w0*x*y + 
   e1*k2*pow(q2,2)*w0*x*y + c2*l2*p2*r1*w0*x*y - 2*c2*k2*q2*r1*w0*x*y + 
   c1*l2*p2*r2*w0*x*y - 2*c1*k2*q2*r2*w0*x*y + e2*l2*o2*u1*w0*x*y - 
   c2*m2*o2*u1*w0*x*y - e2*j2*q2*u1*w0*x*y + c2*j2*r2*u1*w0*x*y + 
   m2*q2*r2*u1*w0*x*y - l2*pow(r2,2)*u1*w0*x*y + e1*l2*o2*u2*w0*x*y - 
   c1*m2*o2*u2*w0*x*y - e1*j2*q2*u2*w0*x*y + c2*j2*r1*u2*w0*x*y + 
   m2*q2*r1*u2*w0*x*y + c1*j2*r2*u2*w0*x*y - 2*l2*r1*r2*u2*w0*x*y - 
   c2*l2*o2*v1*w0*x*y + c2*j2*q2*v1*w0*x*y - m2*pow(q2,2)*v1*w0*x*y + 
   l2*q2*r2*v1*w0*x*y - c1*l2*o2*v2*w0*x*y + c1*j2*q2*v2*w0*x*y + 
   l2*q2*r1*v2*w0*x*y + 2*c0*c2*k2*o2*w1*x*y - 2*c0*c2*j2*p2*w1*x*y - 
   e0*l2*p2*q2*w1*x*y + c0*m2*p2*q2*w1*x*y + e0*k2*pow(q2,2)*w1*x*y + 
   c2*l2*p2*r0*w1*x*y - 2*c2*k2*q2*r0*w1*x*y + c0*l2*p2*r2*w1*x*y - 
   2*c0*k2*q2*r2*w1*x*y + e2*l2*o2*u0*w1*x*y - c2*m2*o2*u0*w1*x*y - 
   e2*j2*q2*u0*w1*x*y + c2*j2*r2*u0*w1*x*y + m2*q2*r2*u0*w1*x*y - 
   l2*pow(r2,2)*u0*w1*x*y + e0*l2*o2*u2*w1*x*y - c0*m2*o2*u2*w1*x*y - 
   e0*j2*q2*u2*w1*x*y + c2*j2*r0*u2*w1*x*y + m2*q2*r0*u2*w1*x*y + 
   c0*j2*r2*u2*w1*x*y - 2*l2*r0*r2*u2*w1*x*y - c2*l2*o2*v0*w1*x*y + 
   c2*j2*q2*v0*w1*x*y - m2*pow(q2,2)*v0*w1*x*y + l2*q2*r2*v0*w1*x*y - 
   c0*l2*o2*v2*w1*x*y + c0*j2*q2*v2*w1*x*y + l2*q2*r0*v2*w1*x*y + 
   2*c0*c1*k2*o2*w2*x*y - 2*c0*c1*j2*p2*w2*x*y + c1*l2*p2*r0*w2*x*y - 
   2*c1*k2*q2*r0*w2*x*y + c0*l2*p2*r1*w2*x*y - 2*c0*k2*q2*r1*w2*x*y + 
   e1*l2*o2*u0*w2*x*y - c1*m2*o2*u0*w2*x*y - e1*j2*q2*u0*w2*x*y + 
   c2*j2*r1*u0*w2*x*y + m2*q2*r1*u0*w2*x*y + c1*j2*r2*u0*w2*x*y - 
   2*l2*r1*r2*u0*w2*x*y + e0*l2*o2*u1*w2*x*y - c0*m2*o2*u1*w2*x*y - 
   e0*j2*q2*u1*w2*x*y + c2*j2*r0*u1*w2*x*y + m2*q2*r0*u1*w2*x*y + 
   c0*j2*r2*u1*w2*x*y - 2*l2*r0*r2*u1*w2*x*y + c1*j2*r0*u2*w2*x*y + 
   c0*j2*r1*u2*w2*x*y - 2*l2*r0*r1*u2*w2*x*y - c1*l2*o2*v0*w2*x*y + 
   c1*j2*q2*v0*w2*x*y + l2*q2*r1*v0*w2*x*y - c0*l2*o2*v1*w2*x*y + 
   c0*j2*q2*v1*w2*x*y + l2*q2*r0*v1*w2*x*y + 
   2*c0*c1*n0*pow(p2,2)*pow(x,2)*y + 
   pow(c0,2)*n1*pow(p2,2)*pow(x,2)*y + c1*d0*k2*p2*r0*pow(x,2)*y + 
   c0*d1*k2*p2*r0*pow(x,2)*y + c0*d0*k2*p2*r1*pow(x,2)*y - 
   2*c0*c1*k2*p2*s0*pow(x,2)*y - pow(c0,2)*k2*p2*s1*pow(x,2)*y - 
   2*c0*c1*n0*o2*t2*pow(x,2)*y - pow(c0,2)*n1*o2*t2*pow(x,2)*y - 
   c1*d0*j2*r0*t2*pow(x,2)*y - c0*d1*j2*r0*t2*pow(x,2)*y + 
   2*c1*n0*q2*r0*t2*pow(x,2)*y + 2*c0*n1*q2*r0*t2*pow(x,2)*y + 
   d1*l2*pow(r0,2)*t2*pow(x,2)*y - c0*d0*j2*r1*t2*pow(x,2)*y + 
   2*c0*n0*q2*r1*t2*pow(x,2)*y + 2*d0*l2*r0*r1*t2*pow(x,2)*y + 
   2*c0*c1*j2*s0*t2*pow(x,2)*y - c1*l2*r0*s0*t2*pow(x,2)*y - 
   c0*l2*r1*s0*t2*pow(x,2)*y + pow(c0,2)*j2*s1*t2*pow(x,2)*y - 
   c0*l2*r0*s1*t2*pow(x,2)*y + d1*e0*k2*o2*u0*pow(x,2)*y + 
   d0*e1*k2*o2*u0*pow(x,2)*y - c1*f0*k2*o2*u0*pow(x,2)*y - 
   c0*f1*k2*o2*u0*pow(x,2)*y - d1*e0*j2*p2*u0*pow(x,2)*y - 
   d0*e1*j2*p2*u0*pow(x,2)*y + c1*f0*j2*p2*u0*pow(x,2)*y + 
   c0*f1*j2*p2*u0*pow(x,2)*y + 2*e1*n0*p2*q2*u0*pow(x,2)*y + 
   2*e0*n1*p2*q2*u0*pow(x,2)*y + f1*l2*p2*r0*u0*pow(x,2)*y + 
   d1*m2*p2*r0*u0*pow(x,2)*y - 2*c2*n1*p2*r0*u0*pow(x,2)*y - 
   2*c1*n2*p2*r0*u0*pow(x,2)*y + f1*k2*q2*r0*u0*pow(x,2)*y + 
   f0*l2*p2*r1*u0*pow(x,2)*y + d0*m2*p2*r1*u0*pow(x,2)*y - 
   2*c2*n0*p2*r1*u0*pow(x,2)*y - 2*c0*n2*p2*r1*u0*pow(x,2)*y + 
   f0*k2*q2*r1*u0*pow(x,2)*y - 2*d2*k2*r0*r1*u0*pow(x,2)*y - 
   2*c1*n0*p2*r2*u0*pow(x,2)*y - 2*c0*n1*p2*r2*u0*pow(x,2)*y - 
   2*d1*k2*r0*r2*u0*pow(x,2)*y - 2*d0*k2*r1*r2*u0*pow(x,2)*y - 
   e1*l2*p2*s0*u0*pow(x,2)*y + c1*m2*p2*s0*u0*pow(x,2)*y - 
   e1*k2*q2*s0*u0*pow(x,2)*y + c2*k2*r1*s0*u0*pow(x,2)*y + 
   c1*k2*r2*s0*u0*pow(x,2)*y - e0*l2*p2*s1*u0*pow(x,2)*y + 
   c0*m2*p2*s1*u0*pow(x,2)*y - e0*k2*q2*s1*u0*pow(x,2)*y + 
   c2*k2*r0*s1*u0*pow(x,2)*y + c0*k2*r2*s1*u0*pow(x,2)*y + 
   c1*k2*r0*s2*u0*pow(x,2)*y + c0*k2*r1*s2*u0*pow(x,2)*y + 
   f1*m2*o2*pow(u0,2)*pow(x,2)*y - e2*n1*o2*pow(u0,2)*pow(x,2)*y - 
   e1*n2*o2*pow(u0,2)*pow(x,2)*y - f2*j2*r1*pow(u0,2)*pow(x,2)*y - 
   f1*j2*r2*pow(u0,2)*pow(x,2)*y + 2*n2*r1*r2*pow(u0,2)*pow(x,2)*y + 
   n1*pow(r2,2)*pow(u0,2)*pow(x,2)*y + 
   e2*j2*s1*pow(u0,2)*pow(x,2)*y - m2*r2*s1*pow(u0,2)*pow(x,2)*y + 
   e1*j2*s2*pow(u0,2)*pow(x,2)*y - m2*r1*s2*pow(u0,2)*pow(x,2)*y + 
   d0*e0*k2*o2*u1*pow(x,2)*y - c0*f0*k2*o2*u1*pow(x,2)*y - 
   d0*e0*j2*p2*u1*pow(x,2)*y + c0*f0*j2*p2*u1*pow(x,2)*y + 
   2*e0*n0*p2*q2*u1*pow(x,2)*y + f0*l2*p2*r0*u1*pow(x,2)*y + 
   d0*m2*p2*r0*u1*pow(x,2)*y - 2*c2*n0*p2*r0*u1*pow(x,2)*y - 
   2*c0*n2*p2*r0*u1*pow(x,2)*y + f0*k2*q2*r0*u1*pow(x,2)*y - 
   d2*k2*pow(r0,2)*u1*pow(x,2)*y - 2*c0*n0*p2*r2*u1*pow(x,2)*y - 
   2*d0*k2*r0*r2*u1*pow(x,2)*y - e0*l2*p2*s0*u1*pow(x,2)*y + 
   c0*m2*p2*s0*u1*pow(x,2)*y - e0*k2*q2*s0*u1*pow(x,2)*y + 
   c2*k2*r0*s0*u1*pow(x,2)*y + c0*k2*r2*s0*u1*pow(x,2)*y + 
   c0*k2*r0*s2*u1*pow(x,2)*y + 2*f0*m2*o2*u0*u1*pow(x,2)*y - 
   2*e2*n0*o2*u0*u1*pow(x,2)*y - 2*e0*n2*o2*u0*u1*pow(x,2)*y - 
   2*f2*j2*r0*u0*u1*pow(x,2)*y - 2*f0*j2*r2*u0*u1*pow(x,2)*y + 
   4*n2*r0*r2*u0*u1*pow(x,2)*y + 2*n0*pow(r2,2)*u0*u1*pow(x,2)*y + 
   2*e2*j2*s0*u0*u1*pow(x,2)*y - 2*m2*r2*s0*u0*u1*pow(x,2)*y + 
   2*e0*j2*s2*u0*u1*pow(x,2)*y - 2*m2*r0*s2*u0*u1*pow(x,2)*y - 
   2*c1*n0*p2*r0*u2*pow(x,2)*y - 2*c0*n1*p2*r0*u2*pow(x,2)*y - 
   d1*k2*pow(r0,2)*u2*pow(x,2)*y - 2*c0*n0*p2*r1*u2*pow(x,2)*y - 
   2*d0*k2*r0*r1*u2*pow(x,2)*y + c1*k2*r0*s0*u2*pow(x,2)*y + 
   c0*k2*r1*s0*u2*pow(x,2)*y + c0*k2*r0*s1*u2*pow(x,2)*y - 
   2*e1*n0*o2*u0*u2*pow(x,2)*y - 2*e0*n1*o2*u0*u2*pow(x,2)*y - 
   2*f1*j2*r0*u0*u2*pow(x,2)*y - 2*f0*j2*r1*u0*u2*pow(x,2)*y + 
   4*n2*r0*r1*u0*u2*pow(x,2)*y + 4*n1*r0*r2*u0*u2*pow(x,2)*y + 
   4*n0*r1*r2*u0*u2*pow(x,2)*y + 2*e1*j2*s0*u0*u2*pow(x,2)*y - 
   2*m2*r1*s0*u0*u2*pow(x,2)*y + 2*e0*j2*s1*u0*u2*pow(x,2)*y - 
   2*m2*r0*s1*u0*u2*pow(x,2)*y - 2*e0*n0*o2*u1*u2*pow(x,2)*y - 
   2*f0*j2*r0*u1*u2*pow(x,2)*y + 2*n2*pow(r0,2)*u1*u2*pow(x,2)*y + 
   4*n0*r0*r2*u1*u2*pow(x,2)*y + 2*e0*j2*s0*u1*u2*pow(x,2)*y - 
   2*m2*r0*s0*u1*u2*pow(x,2)*y + n1*pow(r0,2)*pow(u2,2)*pow(x,2)*y + 
   2*n0*r0*r1*pow(u2,2)*pow(x,2)*y - c1*d0*k2*o2*v0*pow(x,2)*y - 
   c0*d1*k2*o2*v0*pow(x,2)*y + c1*d0*j2*p2*v0*pow(x,2)*y + 
   c0*d1*j2*p2*v0*pow(x,2)*y - 2*c1*n0*p2*q2*v0*pow(x,2)*y - 
   2*c0*n1*p2*q2*v0*pow(x,2)*y - 2*d1*l2*p2*r0*v0*pow(x,2)*y + 
   d1*k2*q2*r0*v0*pow(x,2)*y - 2*d0*l2*p2*r1*v0*pow(x,2)*y + 
   d0*k2*q2*r1*v0*pow(x,2)*y + c1*l2*p2*s0*v0*pow(x,2)*y + 
   c1*k2*q2*s0*v0*pow(x,2)*y + c0*l2*p2*s1*v0*pow(x,2)*y + 
   c0*k2*q2*s1*v0*pow(x,2)*y - f1*l2*o2*u0*v0*pow(x,2)*y - 
   d1*m2*o2*u0*v0*pow(x,2)*y + 2*c2*n1*o2*u0*v0*pow(x,2)*y + 
   2*c1*n2*o2*u0*v0*pow(x,2)*y + f1*j2*q2*u0*v0*pow(x,2)*y + 
   d2*j2*r1*u0*v0*pow(x,2)*y - 2*n2*q2*r1*u0*v0*pow(x,2)*y + 
   d1*j2*r2*u0*v0*pow(x,2)*y - 2*n1*q2*r2*u0*v0*pow(x,2)*y - 
   2*c2*j2*s1*u0*v0*pow(x,2)*y + m2*q2*s1*u0*v0*pow(x,2)*y + 
   l2*r2*s1*u0*v0*pow(x,2)*y - 2*c1*j2*s2*u0*v0*pow(x,2)*y + 
   l2*r1*s2*u0*v0*pow(x,2)*y - f0*l2*o2*u1*v0*pow(x,2)*y - 
   d0*m2*o2*u1*v0*pow(x,2)*y + 2*c2*n0*o2*u1*v0*pow(x,2)*y + 
   2*c0*n2*o2*u1*v0*pow(x,2)*y + f0*j2*q2*u1*v0*pow(x,2)*y + 
   d2*j2*r0*u1*v0*pow(x,2)*y - 2*n2*q2*r0*u1*v0*pow(x,2)*y + 
   d0*j2*r2*u1*v0*pow(x,2)*y - 2*n0*q2*r2*u1*v0*pow(x,2)*y - 
   2*c2*j2*s0*u1*v0*pow(x,2)*y + m2*q2*s0*u1*v0*pow(x,2)*y + 
   l2*r2*s0*u1*v0*pow(x,2)*y - 2*c0*j2*s2*u1*v0*pow(x,2)*y + 
   l2*r0*s2*u1*v0*pow(x,2)*y + 2*c1*n0*o2*u2*v0*pow(x,2)*y + 
   2*c0*n1*o2*u2*v0*pow(x,2)*y + d1*j2*r0*u2*v0*pow(x,2)*y - 
   2*n1*q2*r0*u2*v0*pow(x,2)*y + d0*j2*r1*u2*v0*pow(x,2)*y - 
   2*n0*q2*r1*u2*v0*pow(x,2)*y - 2*c1*j2*s0*u2*v0*pow(x,2)*y + 
   l2*r1*s0*u2*v0*pow(x,2)*y - 2*c0*j2*s1*u2*v0*pow(x,2)*y + 
   l2*r0*s1*u2*v0*pow(x,2)*y + d1*l2*o2*pow(v0,2)*pow(x,2)*y - 
   d1*j2*q2*pow(v0,2)*pow(x,2)*y + 
   n1*pow(q2,2)*pow(v0,2)*pow(x,2)*y - 
   l2*q2*s1*pow(v0,2)*pow(x,2)*y - c0*d0*k2*o2*v1*pow(x,2)*y + 
   c0*d0*j2*p2*v1*pow(x,2)*y - 2*c0*n0*p2*q2*v1*pow(x,2)*y - 
   2*d0*l2*p2*r0*v1*pow(x,2)*y + d0*k2*q2*r0*v1*pow(x,2)*y + 
   c0*l2*p2*s0*v1*pow(x,2)*y + c0*k2*q2*s0*v1*pow(x,2)*y - 
   f0*l2*o2*u0*v1*pow(x,2)*y - d0*m2*o2*u0*v1*pow(x,2)*y + 
   2*c2*n0*o2*u0*v1*pow(x,2)*y + 2*c0*n2*o2*u0*v1*pow(x,2)*y + 
   f0*j2*q2*u0*v1*pow(x,2)*y + d2*j2*r0*u0*v1*pow(x,2)*y - 
   2*n2*q2*r0*u0*v1*pow(x,2)*y + d0*j2*r2*u0*v1*pow(x,2)*y - 
   2*n0*q2*r2*u0*v1*pow(x,2)*y - 2*c2*j2*s0*u0*v1*pow(x,2)*y + 
   m2*q2*s0*u0*v1*pow(x,2)*y + l2*r2*s0*u0*v1*pow(x,2)*y - 
   2*c0*j2*s2*u0*v1*pow(x,2)*y + l2*r0*s2*u0*v1*pow(x,2)*y + 
   2*c0*n0*o2*u2*v1*pow(x,2)*y + d0*j2*r0*u2*v1*pow(x,2)*y - 
   2*n0*q2*r0*u2*v1*pow(x,2)*y - 2*c0*j2*s0*u2*v1*pow(x,2)*y + 
   l2*r0*s0*u2*v1*pow(x,2)*y + 2*d0*l2*o2*v0*v1*pow(x,2)*y - 
   2*d0*j2*q2*v0*v1*pow(x,2)*y + 2*n0*pow(q2,2)*v0*v1*pow(x,2)*y - 
   2*l2*q2*s0*v0*v1*pow(x,2)*y + 2*c1*n0*o2*u0*v2*pow(x,2)*y + 
   2*c0*n1*o2*u0*v2*pow(x,2)*y + d1*j2*r0*u0*v2*pow(x,2)*y - 
   2*n1*q2*r0*u0*v2*pow(x,2)*y + d0*j2*r1*u0*v2*pow(x,2)*y - 
   2*n0*q2*r1*u0*v2*pow(x,2)*y - 2*c1*j2*s0*u0*v2*pow(x,2)*y + 
   l2*r1*s0*u0*v2*pow(x,2)*y - 2*c0*j2*s1*u0*v2*pow(x,2)*y + 
   l2*r0*s1*u0*v2*pow(x,2)*y + 2*c0*n0*o2*u1*v2*pow(x,2)*y + 
   d0*j2*r0*u1*v2*pow(x,2)*y - 2*n0*q2*r0*u1*v2*pow(x,2)*y - 
   2*c0*j2*s0*u1*v2*pow(x,2)*y + l2*r0*s0*u1*v2*pow(x,2)*y + 
   2*c0*c1*k2*o2*w0*pow(x,2)*y - 2*c0*c1*j2*p2*w0*pow(x,2)*y + 
   c1*l2*p2*r0*w0*pow(x,2)*y - 2*c1*k2*q2*r0*w0*pow(x,2)*y + 
   c0*l2*p2*r1*w0*pow(x,2)*y - 2*c0*k2*q2*r1*w0*pow(x,2)*y + 
   e1*l2*o2*u0*w0*pow(x,2)*y - c1*m2*o2*u0*w0*pow(x,2)*y - 
   e1*j2*q2*u0*w0*pow(x,2)*y + c2*j2*r1*u0*w0*pow(x,2)*y + 
   m2*q2*r1*u0*w0*pow(x,2)*y + c1*j2*r2*u0*w0*pow(x,2)*y - 
   2*l2*r1*r2*u0*w0*pow(x,2)*y + e0*l2*o2*u1*w0*pow(x,2)*y - 
   c0*m2*o2*u1*w0*pow(x,2)*y - e0*j2*q2*u1*w0*pow(x,2)*y + 
   c2*j2*r0*u1*w0*pow(x,2)*y + m2*q2*r0*u1*w0*pow(x,2)*y + 
   c0*j2*r2*u1*w0*pow(x,2)*y - 2*l2*r0*r2*u1*w0*pow(x,2)*y + 
   c1*j2*r0*u2*w0*pow(x,2)*y + c0*j2*r1*u2*w0*pow(x,2)*y - 
   2*l2*r0*r1*u2*w0*pow(x,2)*y - c1*l2*o2*v0*w0*pow(x,2)*y + 
   c1*j2*q2*v0*w0*pow(x,2)*y + l2*q2*r1*v0*w0*pow(x,2)*y - 
   c0*l2*o2*v1*w0*pow(x,2)*y + c0*j2*q2*v1*w0*pow(x,2)*y + 
   l2*q2*r0*v1*w0*pow(x,2)*y + pow(c0,2)*k2*o2*w1*pow(x,2)*y - 
   pow(c0,2)*j2*p2*w1*pow(x,2)*y + c0*l2*p2*r0*w1*pow(x,2)*y - 
   2*c0*k2*q2*r0*w1*pow(x,2)*y + e0*l2*o2*u0*w1*pow(x,2)*y - 
   c0*m2*o2*u0*w1*pow(x,2)*y - e0*j2*q2*u0*w1*pow(x,2)*y + 
   c2*j2*r0*u0*w1*pow(x,2)*y + m2*q2*r0*u0*w1*pow(x,2)*y + 
   c0*j2*r2*u0*w1*pow(x,2)*y - 2*l2*r0*r2*u0*w1*pow(x,2)*y + 
   c0*j2*r0*u2*w1*pow(x,2)*y - l2*pow(r0,2)*u2*w1*pow(x,2)*y - 
   c0*l2*o2*v0*w1*pow(x,2)*y + c0*j2*q2*v0*w1*pow(x,2)*y + 
   l2*q2*r0*v0*w1*pow(x,2)*y + c1*j2*r0*u0*w2*pow(x,2)*y + 
   c0*j2*r1*u0*w2*pow(x,2)*y - 2*l2*r0*r1*u0*w2*pow(x,2)*y + 
   c0*j2*r0*u1*w2*pow(x,2)*y - l2*pow(r0,2)*u1*w2*pow(x,2)*y - 
   2*c1*n0*p2*r0*u0*pow(x,3)*y - 2*c0*n1*p2*r0*u0*pow(x,3)*y - 
   d1*k2*pow(r0,2)*u0*pow(x,3)*y - 2*c0*n0*p2*r1*u0*pow(x,3)*y - 
   2*d0*k2*r0*r1*u0*pow(x,3)*y + c1*k2*r0*s0*u0*pow(x,3)*y + 
   c0*k2*r1*s0*u0*pow(x,3)*y + c0*k2*r0*s1*u0*pow(x,3)*y - 
   e1*n0*o2*pow(u0,2)*pow(x,3)*y - e0*n1*o2*pow(u0,2)*pow(x,3)*y - 
   f1*j2*r0*pow(u0,2)*pow(x,3)*y - f0*j2*r1*pow(u0,2)*pow(x,3)*y + 
   2*n2*r0*r1*pow(u0,2)*pow(x,3)*y + 
   2*n1*r0*r2*pow(u0,2)*pow(x,3)*y + 
   2*n0*r1*r2*pow(u0,2)*pow(x,3)*y + e1*j2*s0*pow(u0,2)*pow(x,3)*y - 
   m2*r1*s0*pow(u0,2)*pow(x,3)*y + e0*j2*s1*pow(u0,2)*pow(x,3)*y - 
   m2*r0*s1*pow(u0,2)*pow(x,3)*y - 2*c0*n0*p2*r0*u1*pow(x,3)*y - 
   d0*k2*pow(r0,2)*u1*pow(x,3)*y + c0*k2*r0*s0*u1*pow(x,3)*y - 
   2*e0*n0*o2*u0*u1*pow(x,3)*y - 2*f0*j2*r0*u0*u1*pow(x,3)*y + 
   2*n2*pow(r0,2)*u0*u1*pow(x,3)*y + 4*n0*r0*r2*u0*u1*pow(x,3)*y + 
   2*e0*j2*s0*u0*u1*pow(x,3)*y - 2*m2*r0*s0*u0*u1*pow(x,3)*y + 
   2*n1*pow(r0,2)*u0*u2*pow(x,3)*y + 4*n0*r0*r1*u0*u2*pow(x,3)*y + 
   2*n0*pow(r0,2)*u1*u2*pow(x,3)*y + 2*c1*n0*o2*u0*v0*pow(x,3)*y + 
   2*c0*n1*o2*u0*v0*pow(x,3)*y + d1*j2*r0*u0*v0*pow(x,3)*y - 
   2*n1*q2*r0*u0*v0*pow(x,3)*y + d0*j2*r1*u0*v0*pow(x,3)*y - 
   2*n0*q2*r1*u0*v0*pow(x,3)*y - 2*c1*j2*s0*u0*v0*pow(x,3)*y + 
   l2*r1*s0*u0*v0*pow(x,3)*y - 2*c0*j2*s1*u0*v0*pow(x,3)*y + 
   l2*r0*s1*u0*v0*pow(x,3)*y + 2*c0*n0*o2*u1*v0*pow(x,3)*y + 
   d0*j2*r0*u1*v0*pow(x,3)*y - 2*n0*q2*r0*u1*v0*pow(x,3)*y - 
   2*c0*j2*s0*u1*v0*pow(x,3)*y + l2*r0*s0*u1*v0*pow(x,3)*y + 
   2*c0*n0*o2*u0*v1*pow(x,3)*y + d0*j2*r0*u0*v1*pow(x,3)*y - 
   2*n0*q2*r0*u0*v1*pow(x,3)*y - 2*c0*j2*s0*u0*v1*pow(x,3)*y + 
   l2*r0*s0*u0*v1*pow(x,3)*y + c1*j2*r0*u0*w0*pow(x,3)*y + 
   c0*j2*r1*u0*w0*pow(x,3)*y - 2*l2*r0*r1*u0*w0*pow(x,3)*y + 
   c0*j2*r0*u1*w0*pow(x,3)*y - l2*pow(r0,2)*u1*w0*pow(x,3)*y + 
   c0*j2*r0*u0*w1*pow(x,3)*y - l2*pow(r0,2)*u0*w1*pow(x,3)*y + 
   n1*pow(r0,2)*pow(u0,2)*pow(x,4)*y + 
   2*n0*r0*r1*pow(u0,2)*pow(x,4)*y + 
   2*n0*pow(r0,2)*u0*u1*pow(x,4)*y + d1*e1*l2*pow(p2,2)*pow(y,2) - 
   c1*f1*l2*pow(p2,2)*pow(y,2) - c1*d1*m2*pow(p2,2)*pow(y,2) + 
   2*c1*c2*n1*pow(p2,2)*pow(y,2) + 
   pow(c1,2)*n2*pow(p2,2)*pow(y,2) - d1*e1*k2*p2*q2*pow(y,2) + 
   c1*f1*k2*p2*q2*pow(y,2) + c2*d1*k2*p2*r1*pow(y,2) + 
   c1*d2*k2*p2*r1*pow(y,2) + c1*d1*k2*p2*r2*pow(y,2) - 
   2*c1*c2*k2*p2*s1*pow(y,2) - pow(c1,2)*k2*p2*s2*pow(y,2) - 
   d1*e1*l2*o2*t2*pow(y,2) + c1*f1*l2*o2*t2*pow(y,2) + 
   c1*d1*m2*o2*t2*pow(y,2) - 2*c1*c2*n1*o2*t2*pow(y,2) - 
   pow(c1,2)*n2*o2*t2*pow(y,2) + d1*e1*j2*q2*t2*pow(y,2) - 
   c1*f1*j2*q2*t2*pow(y,2) - e1*n1*pow(q2,2)*t2*pow(y,2) - 
   c2*d1*j2*r1*t2*pow(y,2) - c1*d2*j2*r1*t2*pow(y,2) - 
   f1*l2*q2*r1*t2*pow(y,2) - d1*m2*q2*r1*t2*pow(y,2) + 
   2*c2*n1*q2*r1*t2*pow(y,2) + 2*c1*n2*q2*r1*t2*pow(y,2) + 
   d2*l2*pow(r1,2)*t2*pow(y,2) - c1*d1*j2*r2*t2*pow(y,2) + 
   2*c1*n1*q2*r2*t2*pow(y,2) + 2*d1*l2*r1*r2*t2*pow(y,2) + 
   2*c1*c2*j2*s1*t2*pow(y,2) + e1*l2*q2*s1*t2*pow(y,2) - 
   c1*m2*q2*s1*t2*pow(y,2) - c2*l2*r1*s1*t2*pow(y,2) - 
   c1*l2*r2*s1*t2*pow(y,2) + pow(c1,2)*j2*s2*t2*pow(y,2) - 
   c1*l2*r1*s2*t2*pow(y,2) + d2*e1*k2*o2*u1*pow(y,2) + 
   d1*e2*k2*o2*u1*pow(y,2) - c2*f1*k2*o2*u1*pow(y,2) - 
   c1*f2*k2*o2*u1*pow(y,2) - d2*e1*j2*p2*u1*pow(y,2) - 
   d1*e2*j2*p2*u1*pow(y,2) + c2*f1*j2*p2*u1*pow(y,2) + 
   c1*f2*j2*p2*u1*pow(y,2) - 2*f1*m2*p2*q2*u1*pow(y,2) + 
   2*e2*n1*p2*q2*u1*pow(y,2) + 2*e1*n2*p2*q2*u1*pow(y,2) + 
   f2*l2*p2*r1*u1*pow(y,2) + d2*m2*p2*r1*u1*pow(y,2) - 
   2*c2*n2*p2*r1*u1*pow(y,2) + f2*k2*q2*r1*u1*pow(y,2) + 
   f1*l2*p2*r2*u1*pow(y,2) + d1*m2*p2*r2*u1*pow(y,2) - 
   2*c2*n1*p2*r2*u1*pow(y,2) - 2*c1*n2*p2*r2*u1*pow(y,2) + 
   f1*k2*q2*r2*u1*pow(y,2) - 2*d2*k2*r1*r2*u1*pow(y,2) - 
   d1*k2*pow(r2,2)*u1*pow(y,2) - e2*l2*p2*s1*u1*pow(y,2) + 
   c2*m2*p2*s1*u1*pow(y,2) - e2*k2*q2*s1*u1*pow(y,2) + 
   c2*k2*r2*s1*u1*pow(y,2) - e1*l2*p2*s2*u1*pow(y,2) + 
   c1*m2*p2*s2*u1*pow(y,2) - e1*k2*q2*s2*u1*pow(y,2) + 
   c2*k2*r1*s2*u1*pow(y,2) + c1*k2*r2*s2*u1*pow(y,2) + 
   f2*m2*o2*pow(u1,2)*pow(y,2) - e2*n2*o2*pow(u1,2)*pow(y,2) - 
   f2*j2*r2*pow(u1,2)*pow(y,2) + n2*pow(r2,2)*pow(u1,2)*pow(y,2) + 
   e2*j2*s2*pow(u1,2)*pow(y,2) - m2*r2*s2*pow(u1,2)*pow(y,2) + 
   d1*e1*k2*o2*u2*pow(y,2) - c1*f1*k2*o2*u2*pow(y,2) - 
   d1*e1*j2*p2*u2*pow(y,2) + c1*f1*j2*p2*u2*pow(y,2) + 
   2*e1*n1*p2*q2*u2*pow(y,2) + f1*l2*p2*r1*u2*pow(y,2) + 
   d1*m2*p2*r1*u2*pow(y,2) - 2*c2*n1*p2*r1*u2*pow(y,2) - 
   2*c1*n2*p2*r1*u2*pow(y,2) + f1*k2*q2*r1*u2*pow(y,2) - 
   d2*k2*pow(r1,2)*u2*pow(y,2) - 2*c1*n1*p2*r2*u2*pow(y,2) - 
   2*d1*k2*r1*r2*u2*pow(y,2) - e1*l2*p2*s1*u2*pow(y,2) + 
   c1*m2*p2*s1*u2*pow(y,2) - e1*k2*q2*s1*u2*pow(y,2) + 
   c2*k2*r1*s1*u2*pow(y,2) + c1*k2*r2*s1*u2*pow(y,2) + 
   c1*k2*r1*s2*u2*pow(y,2) + 2*f1*m2*o2*u1*u2*pow(y,2) - 
   2*e2*n1*o2*u1*u2*pow(y,2) - 2*e1*n2*o2*u1*u2*pow(y,2) - 
   2*f2*j2*r1*u1*u2*pow(y,2) - 2*f1*j2*r2*u1*u2*pow(y,2) + 
   4*n2*r1*r2*u1*u2*pow(y,2) + 2*n1*pow(r2,2)*u1*u2*pow(y,2) + 
   2*e2*j2*s1*u1*u2*pow(y,2) - 2*m2*r2*s1*u1*u2*pow(y,2) + 
   2*e1*j2*s2*u1*u2*pow(y,2) - 2*m2*r1*s2*u1*u2*pow(y,2) - 
   e1*n1*o2*pow(u2,2)*pow(y,2) - f1*j2*r1*pow(u2,2)*pow(y,2) + 
   n2*pow(r1,2)*pow(u2,2)*pow(y,2) + 
   2*n1*r1*r2*pow(u2,2)*pow(y,2) + e1*j2*s1*pow(u2,2)*pow(y,2) - 
   m2*r1*s1*pow(u2,2)*pow(y,2) - c2*d1*k2*o2*v1*pow(y,2) - 
   c1*d2*k2*o2*v1*pow(y,2) + c2*d1*j2*p2*v1*pow(y,2) + 
   c1*d2*j2*p2*v1*pow(y,2) + f1*l2*p2*q2*v1*pow(y,2) + 
   d1*m2*p2*q2*v1*pow(y,2) - 2*c2*n1*p2*q2*v1*pow(y,2) - 
   2*c1*n2*p2*q2*v1*pow(y,2) - f1*k2*pow(q2,2)*v1*pow(y,2) - 
   2*d2*l2*p2*r1*v1*pow(y,2) + d2*k2*q2*r1*v1*pow(y,2) - 
   2*d1*l2*p2*r2*v1*pow(y,2) + d1*k2*q2*r2*v1*pow(y,2) + 
   c2*l2*p2*s1*v1*pow(y,2) + c2*k2*q2*s1*v1*pow(y,2) + 
   c1*l2*p2*s2*v1*pow(y,2) + c1*k2*q2*s2*v1*pow(y,2) - 
   f2*l2*o2*u1*v1*pow(y,2) - d2*m2*o2*u1*v1*pow(y,2) + 
   2*c2*n2*o2*u1*v1*pow(y,2) + f2*j2*q2*u1*v1*pow(y,2) + 
   d2*j2*r2*u1*v1*pow(y,2) - 2*n2*q2*r2*u1*v1*pow(y,2) - 
   2*c2*j2*s2*u1*v1*pow(y,2) + m2*q2*s2*u1*v1*pow(y,2) + 
   l2*r2*s2*u1*v1*pow(y,2) - f1*l2*o2*u2*v1*pow(y,2) - 
   d1*m2*o2*u2*v1*pow(y,2) + 2*c2*n1*o2*u2*v1*pow(y,2) + 
   2*c1*n2*o2*u2*v1*pow(y,2) + f1*j2*q2*u2*v1*pow(y,2) + 
   d2*j2*r1*u2*v1*pow(y,2) - 2*n2*q2*r1*u2*v1*pow(y,2) + 
   d1*j2*r2*u2*v1*pow(y,2) - 2*n1*q2*r2*u2*v1*pow(y,2) - 
   2*c2*j2*s1*u2*v1*pow(y,2) + m2*q2*s1*u2*v1*pow(y,2) + 
   l2*r2*s1*u2*v1*pow(y,2) - 2*c1*j2*s2*u2*v1*pow(y,2) + 
   l2*r1*s2*u2*v1*pow(y,2) + d2*l2*o2*pow(v1,2)*pow(y,2) - 
   d2*j2*q2*pow(v1,2)*pow(y,2) + n2*pow(q2,2)*pow(v1,2)*pow(y,2) - 
   l2*q2*s2*pow(v1,2)*pow(y,2) - c1*d1*k2*o2*v2*pow(y,2) + 
   c1*d1*j2*p2*v2*pow(y,2) - 2*c1*n1*p2*q2*v2*pow(y,2) - 
   2*d1*l2*p2*r1*v2*pow(y,2) + d1*k2*q2*r1*v2*pow(y,2) + 
   c1*l2*p2*s1*v2*pow(y,2) + c1*k2*q2*s1*v2*pow(y,2) - 
   f1*l2*o2*u1*v2*pow(y,2) - d1*m2*o2*u1*v2*pow(y,2) + 
   2*c2*n1*o2*u1*v2*pow(y,2) + 2*c1*n2*o2*u1*v2*pow(y,2) + 
   f1*j2*q2*u1*v2*pow(y,2) + d2*j2*r1*u1*v2*pow(y,2) - 
   2*n2*q2*r1*u1*v2*pow(y,2) + d1*j2*r2*u1*v2*pow(y,2) - 
   2*n1*q2*r2*u1*v2*pow(y,2) - 2*c2*j2*s1*u1*v2*pow(y,2) + 
   m2*q2*s1*u1*v2*pow(y,2) + l2*r2*s1*u1*v2*pow(y,2) - 
   2*c1*j2*s2*u1*v2*pow(y,2) + l2*r1*s2*u1*v2*pow(y,2) + 
   2*c1*n1*o2*u2*v2*pow(y,2) + d1*j2*r1*u2*v2*pow(y,2) - 
   2*n1*q2*r1*u2*v2*pow(y,2) - 2*c1*j2*s1*u2*v2*pow(y,2) + 
   l2*r1*s1*u2*v2*pow(y,2) + 2*d1*l2*o2*v1*v2*pow(y,2) - 
   2*d1*j2*q2*v1*v2*pow(y,2) + 2*n1*pow(q2,2)*v1*v2*pow(y,2) - 
   2*l2*q2*s1*v1*v2*pow(y,2) + 2*c1*c2*k2*o2*w1*pow(y,2) - 
   2*c1*c2*j2*p2*w1*pow(y,2) - e1*l2*p2*q2*w1*pow(y,2) + 
   c1*m2*p2*q2*w1*pow(y,2) + e1*k2*pow(q2,2)*w1*pow(y,2) + 
   c2*l2*p2*r1*w1*pow(y,2) - 2*c2*k2*q2*r1*w1*pow(y,2) + 
   c1*l2*p2*r2*w1*pow(y,2) - 2*c1*k2*q2*r2*w1*pow(y,2) + 
   e2*l2*o2*u1*w1*pow(y,2) - c2*m2*o2*u1*w1*pow(y,2) - 
   e2*j2*q2*u1*w1*pow(y,2) + c2*j2*r2*u1*w1*pow(y,2) + 
   m2*q2*r2*u1*w1*pow(y,2) - l2*pow(r2,2)*u1*w1*pow(y,2) + 
   e1*l2*o2*u2*w1*pow(y,2) - c1*m2*o2*u2*w1*pow(y,2) - 
   e1*j2*q2*u2*w1*pow(y,2) + c2*j2*r1*u2*w1*pow(y,2) + 
   m2*q2*r1*u2*w1*pow(y,2) + c1*j2*r2*u2*w1*pow(y,2) - 
   2*l2*r1*r2*u2*w1*pow(y,2) - c2*l2*o2*v1*w1*pow(y,2) + 
   c2*j2*q2*v1*w1*pow(y,2) - m2*pow(q2,2)*v1*w1*pow(y,2) + 
   l2*q2*r2*v1*w1*pow(y,2) - c1*l2*o2*v2*w1*pow(y,2) + 
   c1*j2*q2*v2*w1*pow(y,2) + l2*q2*r1*v2*w1*pow(y,2) + 
   pow(c1,2)*k2*o2*w2*pow(y,2) - pow(c1,2)*j2*p2*w2*pow(y,2) + 
   c1*l2*p2*r1*w2*pow(y,2) - 2*c1*k2*q2*r1*w2*pow(y,2) + 
   e1*l2*o2*u1*w2*pow(y,2) - c1*m2*o2*u1*w2*pow(y,2) - 
   e1*j2*q2*u1*w2*pow(y,2) + c2*j2*r1*u1*w2*pow(y,2) + 
   m2*q2*r1*u1*w2*pow(y,2) + c1*j2*r2*u1*w2*pow(y,2) - 
   2*l2*r1*r2*u1*w2*pow(y,2) + c1*j2*r1*u2*w2*pow(y,2) - 
   l2*pow(r1,2)*u2*w2*pow(y,2) - c1*l2*o2*v1*w2*pow(y,2) + 
   c1*j2*q2*v1*w2*pow(y,2) + l2*q2*r1*v1*w2*pow(y,2) + 
   pow(c1,2)*n0*pow(p2,2)*x*pow(y,2) + 
   2*c0*c1*n1*pow(p2,2)*x*pow(y,2) + c1*d1*k2*p2*r0*x*pow(y,2) + 
   c1*d0*k2*p2*r1*x*pow(y,2) + c0*d1*k2*p2*r1*x*pow(y,2) - 
   pow(c1,2)*k2*p2*s0*x*pow(y,2) - 2*c0*c1*k2*p2*s1*x*pow(y,2) - 
   pow(c1,2)*n0*o2*t2*x*pow(y,2) - 2*c0*c1*n1*o2*t2*x*pow(y,2) - 
   c1*d1*j2*r0*t2*x*pow(y,2) + 2*c1*n1*q2*r0*t2*x*pow(y,2) - 
   c1*d0*j2*r1*t2*x*pow(y,2) - c0*d1*j2*r1*t2*x*pow(y,2) + 
   2*c1*n0*q2*r1*t2*x*pow(y,2) + 2*c0*n1*q2*r1*t2*x*pow(y,2) + 
   2*d1*l2*r0*r1*t2*x*pow(y,2) + d0*l2*pow(r1,2)*t2*x*pow(y,2) + 
   pow(c1,2)*j2*s0*t2*x*pow(y,2) - c1*l2*r1*s0*t2*x*pow(y,2) + 
   2*c0*c1*j2*s1*t2*x*pow(y,2) - c1*l2*r0*s1*t2*x*pow(y,2) - 
   c0*l2*r1*s1*t2*x*pow(y,2) + d1*e1*k2*o2*u0*x*pow(y,2) - 
   c1*f1*k2*o2*u0*x*pow(y,2) - d1*e1*j2*p2*u0*x*pow(y,2) + 
   c1*f1*j2*p2*u0*x*pow(y,2) + 2*e1*n1*p2*q2*u0*x*pow(y,2) + 
   f1*l2*p2*r1*u0*x*pow(y,2) + d1*m2*p2*r1*u0*x*pow(y,2) - 
   2*c2*n1*p2*r1*u0*x*pow(y,2) - 2*c1*n2*p2*r1*u0*x*pow(y,2) + 
   f1*k2*q2*r1*u0*x*pow(y,2) - d2*k2*pow(r1,2)*u0*x*pow(y,2) - 
   2*c1*n1*p2*r2*u0*x*pow(y,2) - 2*d1*k2*r1*r2*u0*x*pow(y,2) - 
   e1*l2*p2*s1*u0*x*pow(y,2) + c1*m2*p2*s1*u0*x*pow(y,2) - 
   e1*k2*q2*s1*u0*x*pow(y,2) + c2*k2*r1*s1*u0*x*pow(y,2) + 
   c1*k2*r2*s1*u0*x*pow(y,2) + c1*k2*r1*s2*u0*x*pow(y,2) + 
   d1*e0*k2*o2*u1*x*pow(y,2) + d0*e1*k2*o2*u1*x*pow(y,2) - 
   c1*f0*k2*o2*u1*x*pow(y,2) - c0*f1*k2*o2*u1*x*pow(y,2) - 
   d1*e0*j2*p2*u1*x*pow(y,2) - d0*e1*j2*p2*u1*x*pow(y,2) + 
   c1*f0*j2*p2*u1*x*pow(y,2) + c0*f1*j2*p2*u1*x*pow(y,2) + 
   2*e1*n0*p2*q2*u1*x*pow(y,2) + 2*e0*n1*p2*q2*u1*x*pow(y,2) + 
   f1*l2*p2*r0*u1*x*pow(y,2) + d1*m2*p2*r0*u1*x*pow(y,2) - 
   2*c2*n1*p2*r0*u1*x*pow(y,2) - 2*c1*n2*p2*r0*u1*x*pow(y,2) + 
   f1*k2*q2*r0*u1*x*pow(y,2) + f0*l2*p2*r1*u1*x*pow(y,2) + 
   d0*m2*p2*r1*u1*x*pow(y,2) - 2*c2*n0*p2*r1*u1*x*pow(y,2) - 
   2*c0*n2*p2*r1*u1*x*pow(y,2) + f0*k2*q2*r1*u1*x*pow(y,2) - 
   2*d2*k2*r0*r1*u1*x*pow(y,2) - 2*c1*n0*p2*r2*u1*x*pow(y,2) - 
   2*c0*n1*p2*r2*u1*x*pow(y,2) - 2*d1*k2*r0*r2*u1*x*pow(y,2) - 
   2*d0*k2*r1*r2*u1*x*pow(y,2) - e1*l2*p2*s0*u1*x*pow(y,2) + 
   c1*m2*p2*s0*u1*x*pow(y,2) - e1*k2*q2*s0*u1*x*pow(y,2) + 
   c2*k2*r1*s0*u1*x*pow(y,2) + c1*k2*r2*s0*u1*x*pow(y,2) - 
   e0*l2*p2*s1*u1*x*pow(y,2) + c0*m2*p2*s1*u1*x*pow(y,2) - 
   e0*k2*q2*s1*u1*x*pow(y,2) + c2*k2*r0*s1*u1*x*pow(y,2) + 
   c0*k2*r2*s1*u1*x*pow(y,2) + c1*k2*r0*s2*u1*x*pow(y,2) + 
   c0*k2*r1*s2*u1*x*pow(y,2) + 2*f1*m2*o2*u0*u1*x*pow(y,2) - 
   2*e2*n1*o2*u0*u1*x*pow(y,2) - 2*e1*n2*o2*u0*u1*x*pow(y,2) - 
   2*f2*j2*r1*u0*u1*x*pow(y,2) - 2*f1*j2*r2*u0*u1*x*pow(y,2) + 
   4*n2*r1*r2*u0*u1*x*pow(y,2) + 2*n1*pow(r2,2)*u0*u1*x*pow(y,2) + 
   2*e2*j2*s1*u0*u1*x*pow(y,2) - 2*m2*r2*s1*u0*u1*x*pow(y,2) + 
   2*e1*j2*s2*u0*u1*x*pow(y,2) - 2*m2*r1*s2*u0*u1*x*pow(y,2) + 
   f0*m2*o2*pow(u1,2)*x*pow(y,2) - e2*n0*o2*pow(u1,2)*x*pow(y,2) - 
   e0*n2*o2*pow(u1,2)*x*pow(y,2) - f2*j2*r0*pow(u1,2)*x*pow(y,2) - 
   f0*j2*r2*pow(u1,2)*x*pow(y,2) + 2*n2*r0*r2*pow(u1,2)*x*pow(y,2) + 
   n0*pow(r2,2)*pow(u1,2)*x*pow(y,2) + 
   e2*j2*s0*pow(u1,2)*x*pow(y,2) - m2*r2*s0*pow(u1,2)*x*pow(y,2) + 
   e0*j2*s2*pow(u1,2)*x*pow(y,2) - m2*r0*s2*pow(u1,2)*x*pow(y,2) - 
   2*c1*n1*p2*r0*u2*x*pow(y,2) - 2*c1*n0*p2*r1*u2*x*pow(y,2) - 
   2*c0*n1*p2*r1*u2*x*pow(y,2) - 2*d1*k2*r0*r1*u2*x*pow(y,2) - 
   d0*k2*pow(r1,2)*u2*x*pow(y,2) + c1*k2*r1*s0*u2*x*pow(y,2) + 
   c1*k2*r0*s1*u2*x*pow(y,2) + c0*k2*r1*s1*u2*x*pow(y,2) - 
   2*e1*n1*o2*u0*u2*x*pow(y,2) - 2*f1*j2*r1*u0*u2*x*pow(y,2) + 
   2*n2*pow(r1,2)*u0*u2*x*pow(y,2) + 4*n1*r1*r2*u0*u2*x*pow(y,2) + 
   2*e1*j2*s1*u0*u2*x*pow(y,2) - 2*m2*r1*s1*u0*u2*x*pow(y,2) - 
   2*e1*n0*o2*u1*u2*x*pow(y,2) - 2*e0*n1*o2*u1*u2*x*pow(y,2) - 
   2*f1*j2*r0*u1*u2*x*pow(y,2) - 2*f0*j2*r1*u1*u2*x*pow(y,2) + 
   4*n2*r0*r1*u1*u2*x*pow(y,2) + 4*n1*r0*r2*u1*u2*x*pow(y,2) + 
   4*n0*r1*r2*u1*u2*x*pow(y,2) + 2*e1*j2*s0*u1*u2*x*pow(y,2) - 
   2*m2*r1*s0*u1*u2*x*pow(y,2) + 2*e0*j2*s1*u1*u2*x*pow(y,2) - 
   2*m2*r0*s1*u1*u2*x*pow(y,2) + 2*n1*r0*r1*pow(u2,2)*x*pow(y,2) + 
   n0*pow(r1,2)*pow(u2,2)*x*pow(y,2) - c1*d1*k2*o2*v0*x*pow(y,2) + 
   c1*d1*j2*p2*v0*x*pow(y,2) - 2*c1*n1*p2*q2*v0*x*pow(y,2) - 
   2*d1*l2*p2*r1*v0*x*pow(y,2) + d1*k2*q2*r1*v0*x*pow(y,2) + 
   c1*l2*p2*s1*v0*x*pow(y,2) + c1*k2*q2*s1*v0*x*pow(y,2) - 
   f1*l2*o2*u1*v0*x*pow(y,2) - d1*m2*o2*u1*v0*x*pow(y,2) + 
   2*c2*n1*o2*u1*v0*x*pow(y,2) + 2*c1*n2*o2*u1*v0*x*pow(y,2) + 
   f1*j2*q2*u1*v0*x*pow(y,2) + d2*j2*r1*u1*v0*x*pow(y,2) - 
   2*n2*q2*r1*u1*v0*x*pow(y,2) + d1*j2*r2*u1*v0*x*pow(y,2) - 
   2*n1*q2*r2*u1*v0*x*pow(y,2) - 2*c2*j2*s1*u1*v0*x*pow(y,2) + 
   m2*q2*s1*u1*v0*x*pow(y,2) + l2*r2*s1*u1*v0*x*pow(y,2) - 
   2*c1*j2*s2*u1*v0*x*pow(y,2) + l2*r1*s2*u1*v0*x*pow(y,2) + 
   2*c1*n1*o2*u2*v0*x*pow(y,2) + d1*j2*r1*u2*v0*x*pow(y,2) - 
   2*n1*q2*r1*u2*v0*x*pow(y,2) - 2*c1*j2*s1*u2*v0*x*pow(y,2) + 
   l2*r1*s1*u2*v0*x*pow(y,2) - c1*d0*k2*o2*v1*x*pow(y,2) - 
   c0*d1*k2*o2*v1*x*pow(y,2) + c1*d0*j2*p2*v1*x*pow(y,2) + 
   c0*d1*j2*p2*v1*x*pow(y,2) - 2*c1*n0*p2*q2*v1*x*pow(y,2) - 
   2*c0*n1*p2*q2*v1*x*pow(y,2) - 2*d1*l2*p2*r0*v1*x*pow(y,2) + 
   d1*k2*q2*r0*v1*x*pow(y,2) - 2*d0*l2*p2*r1*v1*x*pow(y,2) + 
   d0*k2*q2*r1*v1*x*pow(y,2) + c1*l2*p2*s0*v1*x*pow(y,2) + 
   c1*k2*q2*s0*v1*x*pow(y,2) + c0*l2*p2*s1*v1*x*pow(y,2) + 
   c0*k2*q2*s1*v1*x*pow(y,2) - f1*l2*o2*u0*v1*x*pow(y,2) - 
   d1*m2*o2*u0*v1*x*pow(y,2) + 2*c2*n1*o2*u0*v1*x*pow(y,2) + 
   2*c1*n2*o2*u0*v1*x*pow(y,2) + f1*j2*q2*u0*v1*x*pow(y,2) + 
   d2*j2*r1*u0*v1*x*pow(y,2) - 2*n2*q2*r1*u0*v1*x*pow(y,2) + 
   d1*j2*r2*u0*v1*x*pow(y,2) - 2*n1*q2*r2*u0*v1*x*pow(y,2) - 
   2*c2*j2*s1*u0*v1*x*pow(y,2) + m2*q2*s1*u0*v1*x*pow(y,2) + 
   l2*r2*s1*u0*v1*x*pow(y,2) - 2*c1*j2*s2*u0*v1*x*pow(y,2) + 
   l2*r1*s2*u0*v1*x*pow(y,2) - f0*l2*o2*u1*v1*x*pow(y,2) - 
   d0*m2*o2*u1*v1*x*pow(y,2) + 2*c2*n0*o2*u1*v1*x*pow(y,2) + 
   2*c0*n2*o2*u1*v1*x*pow(y,2) + f0*j2*q2*u1*v1*x*pow(y,2) + 
   d2*j2*r0*u1*v1*x*pow(y,2) - 2*n2*q2*r0*u1*v1*x*pow(y,2) + 
   d0*j2*r2*u1*v1*x*pow(y,2) - 2*n0*q2*r2*u1*v1*x*pow(y,2) - 
   2*c2*j2*s0*u1*v1*x*pow(y,2) + m2*q2*s0*u1*v1*x*pow(y,2) + 
   l2*r2*s0*u1*v1*x*pow(y,2) - 2*c0*j2*s2*u1*v1*x*pow(y,2) + 
   l2*r0*s2*u1*v1*x*pow(y,2) + 2*c1*n0*o2*u2*v1*x*pow(y,2) + 
   2*c0*n1*o2*u2*v1*x*pow(y,2) + d1*j2*r0*u2*v1*x*pow(y,2) - 
   2*n1*q2*r0*u2*v1*x*pow(y,2) + d0*j2*r1*u2*v1*x*pow(y,2) - 
   2*n0*q2*r1*u2*v1*x*pow(y,2) - 2*c1*j2*s0*u2*v1*x*pow(y,2) + 
   l2*r1*s0*u2*v1*x*pow(y,2) - 2*c0*j2*s1*u2*v1*x*pow(y,2) + 
   l2*r0*s1*u2*v1*x*pow(y,2) + 2*d1*l2*o2*v0*v1*x*pow(y,2) - 
   2*d1*j2*q2*v0*v1*x*pow(y,2) + 2*n1*pow(q2,2)*v0*v1*x*pow(y,2) - 
   2*l2*q2*s1*v0*v1*x*pow(y,2) + d0*l2*o2*pow(v1,2)*x*pow(y,2) - 
   d0*j2*q2*pow(v1,2)*x*pow(y,2) + 
   n0*pow(q2,2)*pow(v1,2)*x*pow(y,2) - 
   l2*q2*s0*pow(v1,2)*x*pow(y,2) + 2*c1*n1*o2*u0*v2*x*pow(y,2) + 
   d1*j2*r1*u0*v2*x*pow(y,2) - 2*n1*q2*r1*u0*v2*x*pow(y,2) - 
   2*c1*j2*s1*u0*v2*x*pow(y,2) + l2*r1*s1*u0*v2*x*pow(y,2) + 
   2*c1*n0*o2*u1*v2*x*pow(y,2) + 2*c0*n1*o2*u1*v2*x*pow(y,2) + 
   d1*j2*r0*u1*v2*x*pow(y,2) - 2*n1*q2*r0*u1*v2*x*pow(y,2) + 
   d0*j2*r1*u1*v2*x*pow(y,2) - 2*n0*q2*r1*u1*v2*x*pow(y,2) - 
   2*c1*j2*s0*u1*v2*x*pow(y,2) + l2*r1*s0*u1*v2*x*pow(y,2) - 
   2*c0*j2*s1*u1*v2*x*pow(y,2) + l2*r0*s1*u1*v2*x*pow(y,2) + 
   pow(c1,2)*k2*o2*w0*x*pow(y,2) - pow(c1,2)*j2*p2*w0*x*pow(y,2) + 
   c1*l2*p2*r1*w0*x*pow(y,2) - 2*c1*k2*q2*r1*w0*x*pow(y,2) + 
   e1*l2*o2*u1*w0*x*pow(y,2) - c1*m2*o2*u1*w0*x*pow(y,2) - 
   e1*j2*q2*u1*w0*x*pow(y,2) + c2*j2*r1*u1*w0*x*pow(y,2) + 
   m2*q2*r1*u1*w0*x*pow(y,2) + c1*j2*r2*u1*w0*x*pow(y,2) - 
   2*l2*r1*r2*u1*w0*x*pow(y,2) + c1*j2*r1*u2*w0*x*pow(y,2) - 
   l2*pow(r1,2)*u2*w0*x*pow(y,2) - c1*l2*o2*v1*w0*x*pow(y,2) + 
   c1*j2*q2*v1*w0*x*pow(y,2) + l2*q2*r1*v1*w0*x*pow(y,2) + 
   2*c0*c1*k2*o2*w1*x*pow(y,2) - 2*c0*c1*j2*p2*w1*x*pow(y,2) + 
   c1*l2*p2*r0*w1*x*pow(y,2) - 2*c1*k2*q2*r0*w1*x*pow(y,2) + 
   c0*l2*p2*r1*w1*x*pow(y,2) - 2*c0*k2*q2*r1*w1*x*pow(y,2) + 
   e1*l2*o2*u0*w1*x*pow(y,2) - c1*m2*o2*u0*w1*x*pow(y,2) - 
   e1*j2*q2*u0*w1*x*pow(y,2) + c2*j2*r1*u0*w1*x*pow(y,2) + 
   m2*q2*r1*u0*w1*x*pow(y,2) + c1*j2*r2*u0*w1*x*pow(y,2) - 
   2*l2*r1*r2*u0*w1*x*pow(y,2) + e0*l2*o2*u1*w1*x*pow(y,2) - 
   c0*m2*o2*u1*w1*x*pow(y,2) - e0*j2*q2*u1*w1*x*pow(y,2) + 
   c2*j2*r0*u1*w1*x*pow(y,2) + m2*q2*r0*u1*w1*x*pow(y,2) + 
   c0*j2*r2*u1*w1*x*pow(y,2) - 2*l2*r0*r2*u1*w1*x*pow(y,2) + 
   c1*j2*r0*u2*w1*x*pow(y,2) + c0*j2*r1*u2*w1*x*pow(y,2) - 
   2*l2*r0*r1*u2*w1*x*pow(y,2) - c1*l2*o2*v0*w1*x*pow(y,2) + 
   c1*j2*q2*v0*w1*x*pow(y,2) + l2*q2*r1*v0*w1*x*pow(y,2) - 
   c0*l2*o2*v1*w1*x*pow(y,2) + c0*j2*q2*v1*w1*x*pow(y,2) + 
   l2*q2*r0*v1*w1*x*pow(y,2) + c1*j2*r1*u0*w2*x*pow(y,2) - 
   l2*pow(r1,2)*u0*w2*x*pow(y,2) + c1*j2*r0*u1*w2*x*pow(y,2) + 
   c0*j2*r1*u1*w2*x*pow(y,2) - 2*l2*r0*r1*u1*w2*x*pow(y,2) - 
   2*c1*n1*p2*r0*u0*pow(x,2)*pow(y,2) - 
   2*c1*n0*p2*r1*u0*pow(x,2)*pow(y,2) - 
   2*c0*n1*p2*r1*u0*pow(x,2)*pow(y,2) - 
   2*d1*k2*r0*r1*u0*pow(x,2)*pow(y,2) - 
   d0*k2*pow(r1,2)*u0*pow(x,2)*pow(y,2) + 
   c1*k2*r1*s0*u0*pow(x,2)*pow(y,2) + 
   c1*k2*r0*s1*u0*pow(x,2)*pow(y,2) + 
   c0*k2*r1*s1*u0*pow(x,2)*pow(y,2) - 
   e1*n1*o2*pow(u0,2)*pow(x,2)*pow(y,2) - 
   f1*j2*r1*pow(u0,2)*pow(x,2)*pow(y,2) + 
   n2*pow(r1,2)*pow(u0,2)*pow(x,2)*pow(y,2) + 
   2*n1*r1*r2*pow(u0,2)*pow(x,2)*pow(y,2) + 
   e1*j2*s1*pow(u0,2)*pow(x,2)*pow(y,2) - 
   m2*r1*s1*pow(u0,2)*pow(x,2)*pow(y,2) - 
   2*c1*n0*p2*r0*u1*pow(x,2)*pow(y,2) - 
   2*c0*n1*p2*r0*u1*pow(x,2)*pow(y,2) - 
   d1*k2*pow(r0,2)*u1*pow(x,2)*pow(y,2) - 
   2*c0*n0*p2*r1*u1*pow(x,2)*pow(y,2) - 
   2*d0*k2*r0*r1*u1*pow(x,2)*pow(y,2) + 
   c1*k2*r0*s0*u1*pow(x,2)*pow(y,2) + 
   c0*k2*r1*s0*u1*pow(x,2)*pow(y,2) + 
   c0*k2*r0*s1*u1*pow(x,2)*pow(y,2) - 
   2*e1*n0*o2*u0*u1*pow(x,2)*pow(y,2) - 
   2*e0*n1*o2*u0*u1*pow(x,2)*pow(y,2) - 
   2*f1*j2*r0*u0*u1*pow(x,2)*pow(y,2) - 
   2*f0*j2*r1*u0*u1*pow(x,2)*pow(y,2) + 
   4*n2*r0*r1*u0*u1*pow(x,2)*pow(y,2) + 
   4*n1*r0*r2*u0*u1*pow(x,2)*pow(y,2) + 
   4*n0*r1*r2*u0*u1*pow(x,2)*pow(y,2) + 
   2*e1*j2*s0*u0*u1*pow(x,2)*pow(y,2) - 
   2*m2*r1*s0*u0*u1*pow(x,2)*pow(y,2) + 
   2*e0*j2*s1*u0*u1*pow(x,2)*pow(y,2) - 
   2*m2*r0*s1*u0*u1*pow(x,2)*pow(y,2) - 
   e0*n0*o2*pow(u1,2)*pow(x,2)*pow(y,2) - 
   f0*j2*r0*pow(u1,2)*pow(x,2)*pow(y,2) + 
   n2*pow(r0,2)*pow(u1,2)*pow(x,2)*pow(y,2) + 
   2*n0*r0*r2*pow(u1,2)*pow(x,2)*pow(y,2) + 
   e0*j2*s0*pow(u1,2)*pow(x,2)*pow(y,2) - 
   m2*r0*s0*pow(u1,2)*pow(x,2)*pow(y,2) + 
   4*n1*r0*r1*u0*u2*pow(x,2)*pow(y,2) + 
   2*n0*pow(r1,2)*u0*u2*pow(x,2)*pow(y,2) + 
   2*n1*pow(r0,2)*u1*u2*pow(x,2)*pow(y,2) + 
   4*n0*r0*r1*u1*u2*pow(x,2)*pow(y,2) + 
   2*c1*n1*o2*u0*v0*pow(x,2)*pow(y,2) + 
   d1*j2*r1*u0*v0*pow(x,2)*pow(y,2) - 
   2*n1*q2*r1*u0*v0*pow(x,2)*pow(y,2) - 
   2*c1*j2*s1*u0*v0*pow(x,2)*pow(y,2) + 
   l2*r1*s1*u0*v0*pow(x,2)*pow(y,2) + 
   2*c1*n0*o2*u1*v0*pow(x,2)*pow(y,2) + 
   2*c0*n1*o2*u1*v0*pow(x,2)*pow(y,2) + 
   d1*j2*r0*u1*v0*pow(x,2)*pow(y,2) - 
   2*n1*q2*r0*u1*v0*pow(x,2)*pow(y,2) + 
   d0*j2*r1*u1*v0*pow(x,2)*pow(y,2) - 
   2*n0*q2*r1*u1*v0*pow(x,2)*pow(y,2) - 
   2*c1*j2*s0*u1*v0*pow(x,2)*pow(y,2) + 
   l2*r1*s0*u1*v0*pow(x,2)*pow(y,2) - 
   2*c0*j2*s1*u1*v0*pow(x,2)*pow(y,2) + 
   l2*r0*s1*u1*v0*pow(x,2)*pow(y,2) + 
   2*c1*n0*o2*u0*v1*pow(x,2)*pow(y,2) + 
   2*c0*n1*o2*u0*v1*pow(x,2)*pow(y,2) + 
   d1*j2*r0*u0*v1*pow(x,2)*pow(y,2) - 
   2*n1*q2*r0*u0*v1*pow(x,2)*pow(y,2) + 
   d0*j2*r1*u0*v1*pow(x,2)*pow(y,2) - 
   2*n0*q2*r1*u0*v1*pow(x,2)*pow(y,2) - 
   2*c1*j2*s0*u0*v1*pow(x,2)*pow(y,2) + 
   l2*r1*s0*u0*v1*pow(x,2)*pow(y,2) - 
   2*c0*j2*s1*u0*v1*pow(x,2)*pow(y,2) + 
   l2*r0*s1*u0*v1*pow(x,2)*pow(y,2) + 
   2*c0*n0*o2*u1*v1*pow(x,2)*pow(y,2) + 
   d0*j2*r0*u1*v1*pow(x,2)*pow(y,2) - 
   2*n0*q2*r0*u1*v1*pow(x,2)*pow(y,2) - 
   2*c0*j2*s0*u1*v1*pow(x,2)*pow(y,2) + 
   l2*r0*s0*u1*v1*pow(x,2)*pow(y,2) + 
   c1*j2*r1*u0*w0*pow(x,2)*pow(y,2) - 
   l2*pow(r1,2)*u0*w0*pow(x,2)*pow(y,2) + 
   c1*j2*r0*u1*w0*pow(x,2)*pow(y,2) + 
   c0*j2*r1*u1*w0*pow(x,2)*pow(y,2) - 
   2*l2*r0*r1*u1*w0*pow(x,2)*pow(y,2) + 
   c1*j2*r0*u0*w1*pow(x,2)*pow(y,2) + 
   c0*j2*r1*u0*w1*pow(x,2)*pow(y,2) - 
   2*l2*r0*r1*u0*w1*pow(x,2)*pow(y,2) + 
   c0*j2*r0*u1*w1*pow(x,2)*pow(y,2) - 
   l2*pow(r0,2)*u1*w1*pow(x,2)*pow(y,2) + 
   2*n1*r0*r1*pow(u0,2)*pow(x,3)*pow(y,2) + 
   n0*pow(r1,2)*pow(u0,2)*pow(x,3)*pow(y,2) + 
   2*n1*pow(r0,2)*u0*u1*pow(x,3)*pow(y,2) + 
   4*n0*r0*r1*u0*u1*pow(x,3)*pow(y,2) + 
   n0*pow(r0,2)*pow(u1,2)*pow(x,3)*pow(y,2) + 
   pow(c1,2)*n1*pow(p2,2)*pow(y,3) + c1*d1*k2*p2*r1*pow(y,3) - 
   pow(c1,2)*k2*p2*s1*pow(y,3) - pow(c1,2)*n1*o2*t2*pow(y,3) - 
   c1*d1*j2*r1*t2*pow(y,3) + 2*c1*n1*q2*r1*t2*pow(y,3) + 
   d1*l2*pow(r1,2)*t2*pow(y,3) + pow(c1,2)*j2*s1*t2*pow(y,3) - 
   c1*l2*r1*s1*t2*pow(y,3) + d1*e1*k2*o2*u1*pow(y,3) - 
   c1*f1*k2*o2*u1*pow(y,3) - d1*e1*j2*p2*u1*pow(y,3) + 
   c1*f1*j2*p2*u1*pow(y,3) + 2*e1*n1*p2*q2*u1*pow(y,3) + 
   f1*l2*p2*r1*u1*pow(y,3) + d1*m2*p2*r1*u1*pow(y,3) - 
   2*c2*n1*p2*r1*u1*pow(y,3) - 2*c1*n2*p2*r1*u1*pow(y,3) + 
   f1*k2*q2*r1*u1*pow(y,3) - d2*k2*pow(r1,2)*u1*pow(y,3) - 
   2*c1*n1*p2*r2*u1*pow(y,3) - 2*d1*k2*r1*r2*u1*pow(y,3) - 
   e1*l2*p2*s1*u1*pow(y,3) + c1*m2*p2*s1*u1*pow(y,3) - 
   e1*k2*q2*s1*u1*pow(y,3) + c2*k2*r1*s1*u1*pow(y,3) + 
   c1*k2*r2*s1*u1*pow(y,3) + c1*k2*r1*s2*u1*pow(y,3) + 
   f1*m2*o2*pow(u1,2)*pow(y,3) - e2*n1*o2*pow(u1,2)*pow(y,3) - 
   e1*n2*o2*pow(u1,2)*pow(y,3) - f2*j2*r1*pow(u1,2)*pow(y,3) - 
   f1*j2*r2*pow(u1,2)*pow(y,3) + 2*n2*r1*r2*pow(u1,2)*pow(y,3) + 
   n1*pow(r2,2)*pow(u1,2)*pow(y,3) + e2*j2*s1*pow(u1,2)*pow(y,3) - 
   m2*r2*s1*pow(u1,2)*pow(y,3) + e1*j2*s2*pow(u1,2)*pow(y,3) - 
   m2*r1*s2*pow(u1,2)*pow(y,3) - 2*c1*n1*p2*r1*u2*pow(y,3) - 
   d1*k2*pow(r1,2)*u2*pow(y,3) + c1*k2*r1*s1*u2*pow(y,3) - 
   2*e1*n1*o2*u1*u2*pow(y,3) - 2*f1*j2*r1*u1*u2*pow(y,3) + 
   2*n2*pow(r1,2)*u1*u2*pow(y,3) + 4*n1*r1*r2*u1*u2*pow(y,3) + 
   2*e1*j2*s1*u1*u2*pow(y,3) - 2*m2*r1*s1*u1*u2*pow(y,3) + 
   n1*pow(r1,2)*pow(u2,2)*pow(y,3) - c1*d1*k2*o2*v1*pow(y,3) + 
   c1*d1*j2*p2*v1*pow(y,3) - 2*c1*n1*p2*q2*v1*pow(y,3) - 
   2*d1*l2*p2*r1*v1*pow(y,3) + d1*k2*q2*r1*v1*pow(y,3) + 
   c1*l2*p2*s1*v1*pow(y,3) + c1*k2*q2*s1*v1*pow(y,3) - 
   f1*l2*o2*u1*v1*pow(y,3) - d1*m2*o2*u1*v1*pow(y,3) + 
   2*c2*n1*o2*u1*v1*pow(y,3) + 2*c1*n2*o2*u1*v1*pow(y,3) + 
   f1*j2*q2*u1*v1*pow(y,3) + d2*j2*r1*u1*v1*pow(y,3) - 
   2*n2*q2*r1*u1*v1*pow(y,3) + d1*j2*r2*u1*v1*pow(y,3) - 
   2*n1*q2*r2*u1*v1*pow(y,3) - 2*c2*j2*s1*u1*v1*pow(y,3) + 
   m2*q2*s1*u1*v1*pow(y,3) + l2*r2*s1*u1*v1*pow(y,3) - 
   2*c1*j2*s2*u1*v1*pow(y,3) + l2*r1*s2*u1*v1*pow(y,3) + 
   2*c1*n1*o2*u2*v1*pow(y,3) + d1*j2*r1*u2*v1*pow(y,3) - 
   2*n1*q2*r1*u2*v1*pow(y,3) - 2*c1*j2*s1*u2*v1*pow(y,3) + 
   l2*r1*s1*u2*v1*pow(y,3) + d1*l2*o2*pow(v1,2)*pow(y,3) - 
   d1*j2*q2*pow(v1,2)*pow(y,3) + n1*pow(q2,2)*pow(v1,2)*pow(y,3) - 
   l2*q2*s1*pow(v1,2)*pow(y,3) + 2*c1*n1*o2*u1*v2*pow(y,3) + 
   d1*j2*r1*u1*v2*pow(y,3) - 2*n1*q2*r1*u1*v2*pow(y,3) - 
   2*c1*j2*s1*u1*v2*pow(y,3) + l2*r1*s1*u1*v2*pow(y,3) + 
   pow(c1,2)*k2*o2*w1*pow(y,3) - pow(c1,2)*j2*p2*w1*pow(y,3) + 
   c1*l2*p2*r1*w1*pow(y,3) - 2*c1*k2*q2*r1*w1*pow(y,3) + 
   e1*l2*o2*u1*w1*pow(y,3) - c1*m2*o2*u1*w1*pow(y,3) - 
   e1*j2*q2*u1*w1*pow(y,3) + c2*j2*r1*u1*w1*pow(y,3) + 
   m2*q2*r1*u1*w1*pow(y,3) + c1*j2*r2*u1*w1*pow(y,3) - 
   2*l2*r1*r2*u1*w1*pow(y,3) + c1*j2*r1*u2*w1*pow(y,3) - 
   l2*pow(r1,2)*u2*w1*pow(y,3) - c1*l2*o2*v1*w1*pow(y,3) + 
   c1*j2*q2*v1*w1*pow(y,3) + l2*q2*r1*v1*w1*pow(y,3) + 
   c1*j2*r1*u1*w2*pow(y,3) - l2*pow(r1,2)*u1*w2*pow(y,3) - 
   2*c1*n1*p2*r1*u0*x*pow(y,3) - d1*k2*pow(r1,2)*u0*x*pow(y,3) + 
   c1*k2*r1*s1*u0*x*pow(y,3) - 2*c1*n1*p2*r0*u1*x*pow(y,3) - 
   2*c1*n0*p2*r1*u1*x*pow(y,3) - 2*c0*n1*p2*r1*u1*x*pow(y,3) - 
   2*d1*k2*r0*r1*u1*x*pow(y,3) - d0*k2*pow(r1,2)*u1*x*pow(y,3) + 
   c1*k2*r1*s0*u1*x*pow(y,3) + c1*k2*r0*s1*u1*x*pow(y,3) + 
   c0*k2*r1*s1*u1*x*pow(y,3) - 2*e1*n1*o2*u0*u1*x*pow(y,3) - 
   2*f1*j2*r1*u0*u1*x*pow(y,3) + 2*n2*pow(r1,2)*u0*u1*x*pow(y,3) + 
   4*n1*r1*r2*u0*u1*x*pow(y,3) + 2*e1*j2*s1*u0*u1*x*pow(y,3) - 
   2*m2*r1*s1*u0*u1*x*pow(y,3) - e1*n0*o2*pow(u1,2)*x*pow(y,3) - 
   e0*n1*o2*pow(u1,2)*x*pow(y,3) - f1*j2*r0*pow(u1,2)*x*pow(y,3) - 
   f0*j2*r1*pow(u1,2)*x*pow(y,3) + 2*n2*r0*r1*pow(u1,2)*x*pow(y,3) + 
   2*n1*r0*r2*pow(u1,2)*x*pow(y,3) + 
   2*n0*r1*r2*pow(u1,2)*x*pow(y,3) + e1*j2*s0*pow(u1,2)*x*pow(y,3) - 
   m2*r1*s0*pow(u1,2)*x*pow(y,3) + e0*j2*s1*pow(u1,2)*x*pow(y,3) - 
   m2*r0*s1*pow(u1,2)*x*pow(y,3) + 2*n1*pow(r1,2)*u0*u2*x*pow(y,3) + 
   4*n1*r0*r1*u1*u2*x*pow(y,3) + 2*n0*pow(r1,2)*u1*u2*x*pow(y,3) + 
   2*c1*n1*o2*u1*v0*x*pow(y,3) + d1*j2*r1*u1*v0*x*pow(y,3) - 
   2*n1*q2*r1*u1*v0*x*pow(y,3) - 2*c1*j2*s1*u1*v0*x*pow(y,3) + 
   l2*r1*s1*u1*v0*x*pow(y,3) + 2*c1*n1*o2*u0*v1*x*pow(y,3) + 
   d1*j2*r1*u0*v1*x*pow(y,3) - 2*n1*q2*r1*u0*v1*x*pow(y,3) - 
   2*c1*j2*s1*u0*v1*x*pow(y,3) + l2*r1*s1*u0*v1*x*pow(y,3) + 
   2*c1*n0*o2*u1*v1*x*pow(y,3) + 2*c0*n1*o2*u1*v1*x*pow(y,3) + 
   d1*j2*r0*u1*v1*x*pow(y,3) - 2*n1*q2*r0*u1*v1*x*pow(y,3) + 
   d0*j2*r1*u1*v1*x*pow(y,3) - 2*n0*q2*r1*u1*v1*x*pow(y,3) - 
   2*c1*j2*s0*u1*v1*x*pow(y,3) + l2*r1*s0*u1*v1*x*pow(y,3) - 
   2*c0*j2*s1*u1*v1*x*pow(y,3) + l2*r0*s1*u1*v1*x*pow(y,3) + 
   c1*j2*r1*u1*w0*x*pow(y,3) - l2*pow(r1,2)*u1*w0*x*pow(y,3) + 
   c1*j2*r1*u0*w1*x*pow(y,3) - l2*pow(r1,2)*u0*w1*x*pow(y,3) + 
   c1*j2*r0*u1*w1*x*pow(y,3) + c0*j2*r1*u1*w1*x*pow(y,3) - 
   2*l2*r0*r1*u1*w1*x*pow(y,3) + 
   n1*pow(r1,2)*pow(u0,2)*pow(x,2)*pow(y,3) + 
   4*n1*r0*r1*u0*u1*pow(x,2)*pow(y,3) + 
   2*n0*pow(r1,2)*u0*u1*pow(x,2)*pow(y,3) + 
   n1*pow(r0,2)*pow(u1,2)*pow(x,2)*pow(y,3) + 
   2*n0*r0*r1*pow(u1,2)*pow(x,2)*pow(y,3) - 
   2*c1*n1*p2*r1*u1*pow(y,4) - d1*k2*pow(r1,2)*u1*pow(y,4) + 
   c1*k2*r1*s1*u1*pow(y,4) - e1*n1*o2*pow(u1,2)*pow(y,4) - 
   f1*j2*r1*pow(u1,2)*pow(y,4) + n2*pow(r1,2)*pow(u1,2)*pow(y,4) + 
   2*n1*r1*r2*pow(u1,2)*pow(y,4) + e1*j2*s1*pow(u1,2)*pow(y,4) - 
   m2*r1*s1*pow(u1,2)*pow(y,4) + 2*n1*pow(r1,2)*u1*u2*pow(y,4) + 
   2*c1*n1*o2*u1*v1*pow(y,4) + d1*j2*r1*u1*v1*pow(y,4) - 
   2*n1*q2*r1*u1*v1*pow(y,4) - 2*c1*j2*s1*u1*v1*pow(y,4) + 
   l2*r1*s1*u1*v1*pow(y,4) + c1*j2*r1*u1*w1*pow(y,4) - 
   l2*pow(r1,2)*u1*w1*pow(y,4) + 2*n1*pow(r1,2)*u0*u1*x*pow(y,4) + 
   2*n1*r0*r1*pow(u1,2)*x*pow(y,4) + 
   n0*pow(r1,2)*pow(u1,2)*x*pow(y,4) + 
   n1*pow(r1,2)*pow(u1,2)*pow(y,5) + f2*m2*pow(p2,2)*x*z0 - 
   e2*n2*pow(p2,2)*x*z0 - f2*k2*p2*r2*x*z0 + e2*k2*p2*s2*x*z0 - 
   f2*m2*o2*t2*x*z0 + e2*n2*o2*t2*x*z0 + f2*j2*r2*t2*x*z0 - 
   n2*pow(r2,2)*t2*x*z0 - e2*j2*s2*t2*x*z0 + m2*r2*s2*t2*x*z0 + 
   f2*k2*o2*v2*x*z0 - f2*j2*p2*v2*x*z0 + 2*n2*p2*r2*v2*x*z0 - 
   m2*p2*s2*v2*x*z0 - k2*r2*s2*v2*x*z0 - n2*o2*pow(v2,2)*x*z0 + 
   j2*s2*pow(v2,2)*x*z0 - e2*k2*o2*w2*x*z0 + e2*j2*p2*w2*x*z0 - 
   m2*p2*r2*w2*x*z0 + k2*pow(r2,2)*w2*x*z0 + m2*o2*v2*w2*x*z0 - 
   j2*r2*v2*w2*x*z0 + f0*m2*pow(p2,2)*pow(x,2)*z0 - 
   e2*n0*pow(p2,2)*pow(x,2)*z0 - e0*n2*pow(p2,2)*pow(x,2)*z0 - 
   f2*k2*p2*r0*pow(x,2)*z0 - f0*k2*p2*r2*pow(x,2)*z0 + 
   e2*k2*p2*s0*pow(x,2)*z0 + e0*k2*p2*s2*pow(x,2)*z0 - 
   f0*m2*o2*t2*pow(x,2)*z0 + e2*n0*o2*t2*pow(x,2)*z0 + 
   e0*n2*o2*t2*pow(x,2)*z0 + f2*j2*r0*t2*pow(x,2)*z0 + 
   f0*j2*r2*t2*pow(x,2)*z0 - 2*n2*r0*r2*t2*pow(x,2)*z0 - 
   n0*pow(r2,2)*t2*pow(x,2)*z0 - e2*j2*s0*t2*pow(x,2)*z0 + 
   m2*r2*s0*t2*pow(x,2)*z0 - e0*j2*s2*t2*pow(x,2)*z0 + 
   m2*r0*s2*t2*pow(x,2)*z0 + f2*k2*o2*v0*pow(x,2)*z0 - 
   f2*j2*p2*v0*pow(x,2)*z0 + 2*n2*p2*r2*v0*pow(x,2)*z0 - 
   m2*p2*s2*v0*pow(x,2)*z0 - k2*r2*s2*v0*pow(x,2)*z0 + 
   f0*k2*o2*v2*pow(x,2)*z0 - f0*j2*p2*v2*pow(x,2)*z0 + 
   2*n2*p2*r0*v2*pow(x,2)*z0 + 2*n0*p2*r2*v2*pow(x,2)*z0 - 
   m2*p2*s0*v2*pow(x,2)*z0 - k2*r2*s0*v2*pow(x,2)*z0 - 
   k2*r0*s2*v2*pow(x,2)*z0 - 2*n2*o2*v0*v2*pow(x,2)*z0 + 
   2*j2*s2*v0*v2*pow(x,2)*z0 - n0*o2*pow(v2,2)*pow(x,2)*z0 + 
   j2*s0*pow(v2,2)*pow(x,2)*z0 - e2*k2*o2*w0*pow(x,2)*z0 + 
   e2*j2*p2*w0*pow(x,2)*z0 - m2*p2*r2*w0*pow(x,2)*z0 + 
   k2*pow(r2,2)*w0*pow(x,2)*z0 + m2*o2*v2*w0*pow(x,2)*z0 - 
   j2*r2*v2*w0*pow(x,2)*z0 - e0*k2*o2*w2*pow(x,2)*z0 + 
   e0*j2*p2*w2*pow(x,2)*z0 - m2*p2*r0*w2*pow(x,2)*z0 + 
   2*k2*r0*r2*w2*pow(x,2)*z0 + m2*o2*v0*w2*pow(x,2)*z0 - 
   j2*r2*v0*w2*pow(x,2)*z0 - j2*r0*v2*w2*pow(x,2)*z0 - 
   e0*n0*pow(p2,2)*pow(x,3)*z0 - f0*k2*p2*r0*pow(x,3)*z0 + 
   e0*k2*p2*s0*pow(x,3)*z0 + e0*n0*o2*t2*pow(x,3)*z0 + 
   f0*j2*r0*t2*pow(x,3)*z0 - n2*pow(r0,2)*t2*pow(x,3)*z0 - 
   2*n0*r0*r2*t2*pow(x,3)*z0 - e0*j2*s0*t2*pow(x,3)*z0 + 
   m2*r0*s0*t2*pow(x,3)*z0 + f0*k2*o2*v0*pow(x,3)*z0 - 
   f0*j2*p2*v0*pow(x,3)*z0 + 2*n2*p2*r0*v0*pow(x,3)*z0 + 
   2*n0*p2*r2*v0*pow(x,3)*z0 - m2*p2*s0*v0*pow(x,3)*z0 - 
   k2*r2*s0*v0*pow(x,3)*z0 - k2*r0*s2*v0*pow(x,3)*z0 - 
   n2*o2*pow(v0,2)*pow(x,3)*z0 + j2*s2*pow(v0,2)*pow(x,3)*z0 + 
   2*n0*p2*r0*v2*pow(x,3)*z0 - k2*r0*s0*v2*pow(x,3)*z0 - 
   2*n0*o2*v0*v2*pow(x,3)*z0 + 2*j2*s0*v0*v2*pow(x,3)*z0 - 
   e0*k2*o2*w0*pow(x,3)*z0 + e0*j2*p2*w0*pow(x,3)*z0 - 
   m2*p2*r0*w0*pow(x,3)*z0 + 2*k2*r0*r2*w0*pow(x,3)*z0 + 
   m2*o2*v0*w0*pow(x,3)*z0 - j2*r2*v0*w0*pow(x,3)*z0 - 
   j2*r0*v2*w0*pow(x,3)*z0 + k2*pow(r0,2)*w2*pow(x,3)*z0 - 
   j2*r0*v0*w2*pow(x,3)*z0 - n0*pow(r0,2)*t2*pow(x,4)*z0 + 
   2*n0*p2*r0*v0*pow(x,4)*z0 - k2*r0*s0*v0*pow(x,4)*z0 - 
   n0*o2*pow(v0,2)*pow(x,4)*z0 + j2*s0*pow(v0,2)*pow(x,4)*z0 + 
   k2*pow(r0,2)*w0*pow(x,4)*z0 - j2*r0*v0*w0*pow(x,4)*z0 + 
   f1*m2*pow(p2,2)*x*y*z0 - e2*n1*pow(p2,2)*x*y*z0 - 
   e1*n2*pow(p2,2)*x*y*z0 - f2*k2*p2*r1*x*y*z0 - f1*k2*p2*r2*x*y*z0 + 
   e2*k2*p2*s1*x*y*z0 + e1*k2*p2*s2*x*y*z0 - f1*m2*o2*t2*x*y*z0 + 
   e2*n1*o2*t2*x*y*z0 + e1*n2*o2*t2*x*y*z0 + f2*j2*r1*t2*x*y*z0 + 
   f1*j2*r2*t2*x*y*z0 - 2*n2*r1*r2*t2*x*y*z0 - n1*pow(r2,2)*t2*x*y*z0 - 
   e2*j2*s1*t2*x*y*z0 + m2*r2*s1*t2*x*y*z0 - e1*j2*s2*t2*x*y*z0 + 
   m2*r1*s2*t2*x*y*z0 + f2*k2*o2*v1*x*y*z0 - f2*j2*p2*v1*x*y*z0 + 
   2*n2*p2*r2*v1*x*y*z0 - m2*p2*s2*v1*x*y*z0 - k2*r2*s2*v1*x*y*z0 + 
   f1*k2*o2*v2*x*y*z0 - f1*j2*p2*v2*x*y*z0 + 2*n2*p2*r1*v2*x*y*z0 + 
   2*n1*p2*r2*v2*x*y*z0 - m2*p2*s1*v2*x*y*z0 - k2*r2*s1*v2*x*y*z0 - 
   k2*r1*s2*v2*x*y*z0 - 2*n2*o2*v1*v2*x*y*z0 + 2*j2*s2*v1*v2*x*y*z0 - 
   n1*o2*pow(v2,2)*x*y*z0 + j2*s1*pow(v2,2)*x*y*z0 - e2*k2*o2*w1*x*y*z0 + 
   e2*j2*p2*w1*x*y*z0 - m2*p2*r2*w1*x*y*z0 + k2*pow(r2,2)*w1*x*y*z0 + 
   m2*o2*v2*w1*x*y*z0 - j2*r2*v2*w1*x*y*z0 - e1*k2*o2*w2*x*y*z0 + 
   e1*j2*p2*w2*x*y*z0 - m2*p2*r1*w2*x*y*z0 + 2*k2*r1*r2*w2*x*y*z0 + 
   m2*o2*v1*w2*x*y*z0 - j2*r2*v1*w2*x*y*z0 - j2*r1*v2*w2*x*y*z0 - 
   e1*n0*pow(p2,2)*pow(x,2)*y*z0 - e0*n1*pow(p2,2)*pow(x,2)*y*z0 - 
   f1*k2*p2*r0*pow(x,2)*y*z0 - f0*k2*p2*r1*pow(x,2)*y*z0 + 
   e1*k2*p2*s0*pow(x,2)*y*z0 + e0*k2*p2*s1*pow(x,2)*y*z0 + 
   e1*n0*o2*t2*pow(x,2)*y*z0 + e0*n1*o2*t2*pow(x,2)*y*z0 + 
   f1*j2*r0*t2*pow(x,2)*y*z0 + f0*j2*r1*t2*pow(x,2)*y*z0 - 
   2*n2*r0*r1*t2*pow(x,2)*y*z0 - 2*n1*r0*r2*t2*pow(x,2)*y*z0 - 
   2*n0*r1*r2*t2*pow(x,2)*y*z0 - e1*j2*s0*t2*pow(x,2)*y*z0 + 
   m2*r1*s0*t2*pow(x,2)*y*z0 - e0*j2*s1*t2*pow(x,2)*y*z0 + 
   m2*r0*s1*t2*pow(x,2)*y*z0 + f1*k2*o2*v0*pow(x,2)*y*z0 - 
   f1*j2*p2*v0*pow(x,2)*y*z0 + 2*n2*p2*r1*v0*pow(x,2)*y*z0 + 
   2*n1*p2*r2*v0*pow(x,2)*y*z0 - m2*p2*s1*v0*pow(x,2)*y*z0 - 
   k2*r2*s1*v0*pow(x,2)*y*z0 - k2*r1*s2*v0*pow(x,2)*y*z0 + 
   f0*k2*o2*v1*pow(x,2)*y*z0 - f0*j2*p2*v1*pow(x,2)*y*z0 + 
   2*n2*p2*r0*v1*pow(x,2)*y*z0 + 2*n0*p2*r2*v1*pow(x,2)*y*z0 - 
   m2*p2*s0*v1*pow(x,2)*y*z0 - k2*r2*s0*v1*pow(x,2)*y*z0 - 
   k2*r0*s2*v1*pow(x,2)*y*z0 - 2*n2*o2*v0*v1*pow(x,2)*y*z0 + 
   2*j2*s2*v0*v1*pow(x,2)*y*z0 + 2*n1*p2*r0*v2*pow(x,2)*y*z0 + 
   2*n0*p2*r1*v2*pow(x,2)*y*z0 - k2*r1*s0*v2*pow(x,2)*y*z0 - 
   k2*r0*s1*v2*pow(x,2)*y*z0 - 2*n1*o2*v0*v2*pow(x,2)*y*z0 + 
   2*j2*s1*v0*v2*pow(x,2)*y*z0 - 2*n0*o2*v1*v2*pow(x,2)*y*z0 + 
   2*j2*s0*v1*v2*pow(x,2)*y*z0 - e1*k2*o2*w0*pow(x,2)*y*z0 + 
   e1*j2*p2*w0*pow(x,2)*y*z0 - m2*p2*r1*w0*pow(x,2)*y*z0 + 
   2*k2*r1*r2*w0*pow(x,2)*y*z0 + m2*o2*v1*w0*pow(x,2)*y*z0 - 
   j2*r2*v1*w0*pow(x,2)*y*z0 - j2*r1*v2*w0*pow(x,2)*y*z0 - 
   e0*k2*o2*w1*pow(x,2)*y*z0 + e0*j2*p2*w1*pow(x,2)*y*z0 - 
   m2*p2*r0*w1*pow(x,2)*y*z0 + 2*k2*r0*r2*w1*pow(x,2)*y*z0 + 
   m2*o2*v0*w1*pow(x,2)*y*z0 - j2*r2*v0*w1*pow(x,2)*y*z0 - 
   j2*r0*v2*w1*pow(x,2)*y*z0 + 2*k2*r0*r1*w2*pow(x,2)*y*z0 - 
   j2*r1*v0*w2*pow(x,2)*y*z0 - j2*r0*v1*w2*pow(x,2)*y*z0 - 
   n1*pow(r0,2)*t2*pow(x,3)*y*z0 - 2*n0*r0*r1*t2*pow(x,3)*y*z0 + 
   2*n1*p2*r0*v0*pow(x,3)*y*z0 + 2*n0*p2*r1*v0*pow(x,3)*y*z0 - 
   k2*r1*s0*v0*pow(x,3)*y*z0 - k2*r0*s1*v0*pow(x,3)*y*z0 - 
   n1*o2*pow(v0,2)*pow(x,3)*y*z0 + j2*s1*pow(v0,2)*pow(x,3)*y*z0 + 
   2*n0*p2*r0*v1*pow(x,3)*y*z0 - k2*r0*s0*v1*pow(x,3)*y*z0 - 
   2*n0*o2*v0*v1*pow(x,3)*y*z0 + 2*j2*s0*v0*v1*pow(x,3)*y*z0 + 
   2*k2*r0*r1*w0*pow(x,3)*y*z0 - j2*r1*v0*w0*pow(x,3)*y*z0 - 
   j2*r0*v1*w0*pow(x,3)*y*z0 + k2*pow(r0,2)*w1*pow(x,3)*y*z0 - 
   j2*r0*v0*w1*pow(x,3)*y*z0 - e1*n1*pow(p2,2)*x*pow(y,2)*z0 - 
   f1*k2*p2*r1*x*pow(y,2)*z0 + e1*k2*p2*s1*x*pow(y,2)*z0 + 
   e1*n1*o2*t2*x*pow(y,2)*z0 + f1*j2*r1*t2*x*pow(y,2)*z0 - 
   n2*pow(r1,2)*t2*x*pow(y,2)*z0 - 2*n1*r1*r2*t2*x*pow(y,2)*z0 - 
   e1*j2*s1*t2*x*pow(y,2)*z0 + m2*r1*s1*t2*x*pow(y,2)*z0 + 
   f1*k2*o2*v1*x*pow(y,2)*z0 - f1*j2*p2*v1*x*pow(y,2)*z0 + 
   2*n2*p2*r1*v1*x*pow(y,2)*z0 + 2*n1*p2*r2*v1*x*pow(y,2)*z0 - 
   m2*p2*s1*v1*x*pow(y,2)*z0 - k2*r2*s1*v1*x*pow(y,2)*z0 - 
   k2*r1*s2*v1*x*pow(y,2)*z0 - n2*o2*pow(v1,2)*x*pow(y,2)*z0 + 
   j2*s2*pow(v1,2)*x*pow(y,2)*z0 + 2*n1*p2*r1*v2*x*pow(y,2)*z0 - 
   k2*r1*s1*v2*x*pow(y,2)*z0 - 2*n1*o2*v1*v2*x*pow(y,2)*z0 + 
   2*j2*s1*v1*v2*x*pow(y,2)*z0 - e1*k2*o2*w1*x*pow(y,2)*z0 + 
   e1*j2*p2*w1*x*pow(y,2)*z0 - m2*p2*r1*w1*x*pow(y,2)*z0 + 
   2*k2*r1*r2*w1*x*pow(y,2)*z0 + m2*o2*v1*w1*x*pow(y,2)*z0 - 
   j2*r2*v1*w1*x*pow(y,2)*z0 - j2*r1*v2*w1*x*pow(y,2)*z0 + 
   k2*pow(r1,2)*w2*x*pow(y,2)*z0 - j2*r1*v1*w2*x*pow(y,2)*z0 - 
   2*n1*r0*r1*t2*pow(x,2)*pow(y,2)*z0 - 
   n0*pow(r1,2)*t2*pow(x,2)*pow(y,2)*z0 + 
   2*n1*p2*r1*v0*pow(x,2)*pow(y,2)*z0 - 
   k2*r1*s1*v0*pow(x,2)*pow(y,2)*z0 + 
   2*n1*p2*r0*v1*pow(x,2)*pow(y,2)*z0 + 
   2*n0*p2*r1*v1*pow(x,2)*pow(y,2)*z0 - 
   k2*r1*s0*v1*pow(x,2)*pow(y,2)*z0 - 
   k2*r0*s1*v1*pow(x,2)*pow(y,2)*z0 - 
   2*n1*o2*v0*v1*pow(x,2)*pow(y,2)*z0 + 
   2*j2*s1*v0*v1*pow(x,2)*pow(y,2)*z0 - 
   n0*o2*pow(v1,2)*pow(x,2)*pow(y,2)*z0 + 
   j2*s0*pow(v1,2)*pow(x,2)*pow(y,2)*z0 + 
   k2*pow(r1,2)*w0*pow(x,2)*pow(y,2)*z0 - 
   j2*r1*v1*w0*pow(x,2)*pow(y,2)*z0 + 
   2*k2*r0*r1*w1*pow(x,2)*pow(y,2)*z0 - 
   j2*r1*v0*w1*pow(x,2)*pow(y,2)*z0 - 
   j2*r0*v1*w1*pow(x,2)*pow(y,2)*z0 - n1*pow(r1,2)*t2*x*pow(y,3)*z0 + 
   2*n1*p2*r1*v1*x*pow(y,3)*z0 - k2*r1*s1*v1*x*pow(y,3)*z0 - 
   n1*o2*pow(v1,2)*x*pow(y,3)*z0 + j2*s1*pow(v1,2)*x*pow(y,3)*z0 + 
   k2*pow(r1,2)*w1*x*pow(y,3)*z0 - j2*r1*v1*w1*x*pow(y,3)*z0 + 
   f2*m2*pow(p2,2)*y*z1 - e2*n2*pow(p2,2)*y*z1 - f2*k2*p2*r2*y*z1 + 
   e2*k2*p2*s2*y*z1 - f2*m2*o2*t2*y*z1 + e2*n2*o2*t2*y*z1 + 
   f2*j2*r2*t2*y*z1 - n2*pow(r2,2)*t2*y*z1 - e2*j2*s2*t2*y*z1 + 
   m2*r2*s2*t2*y*z1 + f2*k2*o2*v2*y*z1 - f2*j2*p2*v2*y*z1 + 
   2*n2*p2*r2*v2*y*z1 - m2*p2*s2*v2*y*z1 - k2*r2*s2*v2*y*z1 - 
   n2*o2*pow(v2,2)*y*z1 + j2*s2*pow(v2,2)*y*z1 - e2*k2*o2*w2*y*z1 + 
   e2*j2*p2*w2*y*z1 - m2*p2*r2*w2*y*z1 + k2*pow(r2,2)*w2*y*z1 + 
   m2*o2*v2*w2*y*z1 - j2*r2*v2*w2*y*z1 + f0*m2*pow(p2,2)*x*y*z1 - 
   e2*n0*pow(p2,2)*x*y*z1 - e0*n2*pow(p2,2)*x*y*z1 - f2*k2*p2*r0*x*y*z1 - 
   f0*k2*p2*r2*x*y*z1 + e2*k2*p2*s0*x*y*z1 + e0*k2*p2*s2*x*y*z1 - 
   f0*m2*o2*t2*x*y*z1 + e2*n0*o2*t2*x*y*z1 + e0*n2*o2*t2*x*y*z1 + 
   f2*j2*r0*t2*x*y*z1 + f0*j2*r2*t2*x*y*z1 - 2*n2*r0*r2*t2*x*y*z1 - 
   n0*pow(r2,2)*t2*x*y*z1 - e2*j2*s0*t2*x*y*z1 + m2*r2*s0*t2*x*y*z1 - 
   e0*j2*s2*t2*x*y*z1 + m2*r0*s2*t2*x*y*z1 + f2*k2*o2*v0*x*y*z1 - 
   f2*j2*p2*v0*x*y*z1 + 2*n2*p2*r2*v0*x*y*z1 - m2*p2*s2*v0*x*y*z1 - 
   k2*r2*s2*v0*x*y*z1 + f0*k2*o2*v2*x*y*z1 - f0*j2*p2*v2*x*y*z1 + 
   2*n2*p2*r0*v2*x*y*z1 + 2*n0*p2*r2*v2*x*y*z1 - m2*p2*s0*v2*x*y*z1 - 
   k2*r2*s0*v2*x*y*z1 - k2*r0*s2*v2*x*y*z1 - 2*n2*o2*v0*v2*x*y*z1 + 
   2*j2*s2*v0*v2*x*y*z1 - n0*o2*pow(v2,2)*x*y*z1 + 
   j2*s0*pow(v2,2)*x*y*z1 - e2*k2*o2*w0*x*y*z1 + e2*j2*p2*w0*x*y*z1 - 
   m2*p2*r2*w0*x*y*z1 + k2*pow(r2,2)*w0*x*y*z1 + m2*o2*v2*w0*x*y*z1 - 
   j2*r2*v2*w0*x*y*z1 - e0*k2*o2*w2*x*y*z1 + e0*j2*p2*w2*x*y*z1 - 
   m2*p2*r0*w2*x*y*z1 + 2*k2*r0*r2*w2*x*y*z1 + m2*o2*v0*w2*x*y*z1 - 
   j2*r2*v0*w2*x*y*z1 - j2*r0*v2*w2*x*y*z1 - 
   e0*n0*pow(p2,2)*pow(x,2)*y*z1 - f0*k2*p2*r0*pow(x,2)*y*z1 + 
   e0*k2*p2*s0*pow(x,2)*y*z1 + e0*n0*o2*t2*pow(x,2)*y*z1 + 
   f0*j2*r0*t2*pow(x,2)*y*z1 - n2*pow(r0,2)*t2*pow(x,2)*y*z1 - 
   2*n0*r0*r2*t2*pow(x,2)*y*z1 - e0*j2*s0*t2*pow(x,2)*y*z1 + 
   m2*r0*s0*t2*pow(x,2)*y*z1 + f0*k2*o2*v0*pow(x,2)*y*z1 - 
   f0*j2*p2*v0*pow(x,2)*y*z1 + 2*n2*p2*r0*v0*pow(x,2)*y*z1 + 
   2*n0*p2*r2*v0*pow(x,2)*y*z1 - m2*p2*s0*v0*pow(x,2)*y*z1 - 
   k2*r2*s0*v0*pow(x,2)*y*z1 - k2*r0*s2*v0*pow(x,2)*y*z1 - 
   n2*o2*pow(v0,2)*pow(x,2)*y*z1 + j2*s2*pow(v0,2)*pow(x,2)*y*z1 + 
   2*n0*p2*r0*v2*pow(x,2)*y*z1 - k2*r0*s0*v2*pow(x,2)*y*z1 - 
   2*n0*o2*v0*v2*pow(x,2)*y*z1 + 2*j2*s0*v0*v2*pow(x,2)*y*z1 - 
   e0*k2*o2*w0*pow(x,2)*y*z1 + e0*j2*p2*w0*pow(x,2)*y*z1 - 
   m2*p2*r0*w0*pow(x,2)*y*z1 + 2*k2*r0*r2*w0*pow(x,2)*y*z1 + 
   m2*o2*v0*w0*pow(x,2)*y*z1 - j2*r2*v0*w0*pow(x,2)*y*z1 - 
   j2*r0*v2*w0*pow(x,2)*y*z1 + k2*pow(r0,2)*w2*pow(x,2)*y*z1 - 
   j2*r0*v0*w2*pow(x,2)*y*z1 - n0*pow(r0,2)*t2*pow(x,3)*y*z1 + 
   2*n0*p2*r0*v0*pow(x,3)*y*z1 - k2*r0*s0*v0*pow(x,3)*y*z1 - 
   n0*o2*pow(v0,2)*pow(x,3)*y*z1 + j2*s0*pow(v0,2)*pow(x,3)*y*z1 + 
   k2*pow(r0,2)*w0*pow(x,3)*y*z1 - j2*r0*v0*w0*pow(x,3)*y*z1 + 
   f1*m2*pow(p2,2)*pow(y,2)*z1 - e2*n1*pow(p2,2)*pow(y,2)*z1 - 
   e1*n2*pow(p2,2)*pow(y,2)*z1 - f2*k2*p2*r1*pow(y,2)*z1 - 
   f1*k2*p2*r2*pow(y,2)*z1 + e2*k2*p2*s1*pow(y,2)*z1 + 
   e1*k2*p2*s2*pow(y,2)*z1 - f1*m2*o2*t2*pow(y,2)*z1 + 
   e2*n1*o2*t2*pow(y,2)*z1 + e1*n2*o2*t2*pow(y,2)*z1 + 
   f2*j2*r1*t2*pow(y,2)*z1 + f1*j2*r2*t2*pow(y,2)*z1 - 
   2*n2*r1*r2*t2*pow(y,2)*z1 - n1*pow(r2,2)*t2*pow(y,2)*z1 - 
   e2*j2*s1*t2*pow(y,2)*z1 + m2*r2*s1*t2*pow(y,2)*z1 - 
   e1*j2*s2*t2*pow(y,2)*z1 + m2*r1*s2*t2*pow(y,2)*z1 + 
   f2*k2*o2*v1*pow(y,2)*z1 - f2*j2*p2*v1*pow(y,2)*z1 + 
   2*n2*p2*r2*v1*pow(y,2)*z1 - m2*p2*s2*v1*pow(y,2)*z1 - 
   k2*r2*s2*v1*pow(y,2)*z1 + f1*k2*o2*v2*pow(y,2)*z1 - 
   f1*j2*p2*v2*pow(y,2)*z1 + 2*n2*p2*r1*v2*pow(y,2)*z1 + 
   2*n1*p2*r2*v2*pow(y,2)*z1 - m2*p2*s1*v2*pow(y,2)*z1 - 
   k2*r2*s1*v2*pow(y,2)*z1 - k2*r1*s2*v2*pow(y,2)*z1 - 
   2*n2*o2*v1*v2*pow(y,2)*z1 + 2*j2*s2*v1*v2*pow(y,2)*z1 - 
   n1*o2*pow(v2,2)*pow(y,2)*z1 + j2*s1*pow(v2,2)*pow(y,2)*z1 - 
   e2*k2*o2*w1*pow(y,2)*z1 + e2*j2*p2*w1*pow(y,2)*z1 - 
   m2*p2*r2*w1*pow(y,2)*z1 + k2*pow(r2,2)*w1*pow(y,2)*z1 + 
   m2*o2*v2*w1*pow(y,2)*z1 - j2*r2*v2*w1*pow(y,2)*z1 - 
   e1*k2*o2*w2*pow(y,2)*z1 + e1*j2*p2*w2*pow(y,2)*z1 - 
   m2*p2*r1*w2*pow(y,2)*z1 + 2*k2*r1*r2*w2*pow(y,2)*z1 + 
   m2*o2*v1*w2*pow(y,2)*z1 - j2*r2*v1*w2*pow(y,2)*z1 - 
   j2*r1*v2*w2*pow(y,2)*z1 - e1*n0*pow(p2,2)*x*pow(y,2)*z1 - 
   e0*n1*pow(p2,2)*x*pow(y,2)*z1 - f1*k2*p2*r0*x*pow(y,2)*z1 - 
   f0*k2*p2*r1*x*pow(y,2)*z1 + e1*k2*p2*s0*x*pow(y,2)*z1 + 
   e0*k2*p2*s1*x*pow(y,2)*z1 + e1*n0*o2*t2*x*pow(y,2)*z1 + 
   e0*n1*o2*t2*x*pow(y,2)*z1 + f1*j2*r0*t2*x*pow(y,2)*z1 + 
   f0*j2*r1*t2*x*pow(y,2)*z1 - 2*n2*r0*r1*t2*x*pow(y,2)*z1 - 
   2*n1*r0*r2*t2*x*pow(y,2)*z1 - 2*n0*r1*r2*t2*x*pow(y,2)*z1 - 
   e1*j2*s0*t2*x*pow(y,2)*z1 + m2*r1*s0*t2*x*pow(y,2)*z1 - 
   e0*j2*s1*t2*x*pow(y,2)*z1 + m2*r0*s1*t2*x*pow(y,2)*z1 + 
   f1*k2*o2*v0*x*pow(y,2)*z1 - f1*j2*p2*v0*x*pow(y,2)*z1 + 
   2*n2*p2*r1*v0*x*pow(y,2)*z1 + 2*n1*p2*r2*v0*x*pow(y,2)*z1 - 
   m2*p2*s1*v0*x*pow(y,2)*z1 - k2*r2*s1*v0*x*pow(y,2)*z1 - 
   k2*r1*s2*v0*x*pow(y,2)*z1 + f0*k2*o2*v1*x*pow(y,2)*z1 - 
   f0*j2*p2*v1*x*pow(y,2)*z1 + 2*n2*p2*r0*v1*x*pow(y,2)*z1 + 
   2*n0*p2*r2*v1*x*pow(y,2)*z1 - m2*p2*s0*v1*x*pow(y,2)*z1 - 
   k2*r2*s0*v1*x*pow(y,2)*z1 - k2*r0*s2*v1*x*pow(y,2)*z1 - 
   2*n2*o2*v0*v1*x*pow(y,2)*z1 + 2*j2*s2*v0*v1*x*pow(y,2)*z1 + 
   2*n1*p2*r0*v2*x*pow(y,2)*z1 + 2*n0*p2*r1*v2*x*pow(y,2)*z1 - 
   k2*r1*s0*v2*x*pow(y,2)*z1 - k2*r0*s1*v2*x*pow(y,2)*z1 - 
   2*n1*o2*v0*v2*x*pow(y,2)*z1 + 2*j2*s1*v0*v2*x*pow(y,2)*z1 - 
   2*n0*o2*v1*v2*x*pow(y,2)*z1 + 2*j2*s0*v1*v2*x*pow(y,2)*z1 - 
   e1*k2*o2*w0*x*pow(y,2)*z1 + e1*j2*p2*w0*x*pow(y,2)*z1 - 
   m2*p2*r1*w0*x*pow(y,2)*z1 + 2*k2*r1*r2*w0*x*pow(y,2)*z1 + 
   m2*o2*v1*w0*x*pow(y,2)*z1 - j2*r2*v1*w0*x*pow(y,2)*z1 - 
   j2*r1*v2*w0*x*pow(y,2)*z1 - e0*k2*o2*w1*x*pow(y,2)*z1 + 
   e0*j2*p2*w1*x*pow(y,2)*z1 - m2*p2*r0*w1*x*pow(y,2)*z1 + 
   2*k2*r0*r2*w1*x*pow(y,2)*z1 + m2*o2*v0*w1*x*pow(y,2)*z1 - 
   j2*r2*v0*w1*x*pow(y,2)*z1 - j2*r0*v2*w1*x*pow(y,2)*z1 + 
   2*k2*r0*r1*w2*x*pow(y,2)*z1 - j2*r1*v0*w2*x*pow(y,2)*z1 - 
   j2*r0*v1*w2*x*pow(y,2)*z1 - n1*pow(r0,2)*t2*pow(x,2)*pow(y,2)*z1 - 
   2*n0*r0*r1*t2*pow(x,2)*pow(y,2)*z1 + 
   2*n1*p2*r0*v0*pow(x,2)*pow(y,2)*z1 + 
   2*n0*p2*r1*v0*pow(x,2)*pow(y,2)*z1 - 
   k2*r1*s0*v0*pow(x,2)*pow(y,2)*z1 - 
   k2*r0*s1*v0*pow(x,2)*pow(y,2)*z1 - 
   n1*o2*pow(v0,2)*pow(x,2)*pow(y,2)*z1 + 
   j2*s1*pow(v0,2)*pow(x,2)*pow(y,2)*z1 + 
   2*n0*p2*r0*v1*pow(x,2)*pow(y,2)*z1 - 
   k2*r0*s0*v1*pow(x,2)*pow(y,2)*z1 - 
   2*n0*o2*v0*v1*pow(x,2)*pow(y,2)*z1 + 
   2*j2*s0*v0*v1*pow(x,2)*pow(y,2)*z1 + 
   2*k2*r0*r1*w0*pow(x,2)*pow(y,2)*z1 - 
   j2*r1*v0*w0*pow(x,2)*pow(y,2)*z1 - 
   j2*r0*v1*w0*pow(x,2)*pow(y,2)*z1 + 
   k2*pow(r0,2)*w1*pow(x,2)*pow(y,2)*z1 - 
   j2*r0*v0*w1*pow(x,2)*pow(y,2)*z1 - e1*n1*pow(p2,2)*pow(y,3)*z1 - 
   f1*k2*p2*r1*pow(y,3)*z1 + e1*k2*p2*s1*pow(y,3)*z1 + 
   e1*n1*o2*t2*pow(y,3)*z1 + f1*j2*r1*t2*pow(y,3)*z1 - 
   n2*pow(r1,2)*t2*pow(y,3)*z1 - 2*n1*r1*r2*t2*pow(y,3)*z1 - 
   e1*j2*s1*t2*pow(y,3)*z1 + m2*r1*s1*t2*pow(y,3)*z1 + 
   f1*k2*o2*v1*pow(y,3)*z1 - f1*j2*p2*v1*pow(y,3)*z1 + 
   2*n2*p2*r1*v1*pow(y,3)*z1 + 2*n1*p2*r2*v1*pow(y,3)*z1 - 
   m2*p2*s1*v1*pow(y,3)*z1 - k2*r2*s1*v1*pow(y,3)*z1 - 
   k2*r1*s2*v1*pow(y,3)*z1 - n2*o2*pow(v1,2)*pow(y,3)*z1 + 
   j2*s2*pow(v1,2)*pow(y,3)*z1 + 2*n1*p2*r1*v2*pow(y,3)*z1 - 
   k2*r1*s1*v2*pow(y,3)*z1 - 2*n1*o2*v1*v2*pow(y,3)*z1 + 
   2*j2*s1*v1*v2*pow(y,3)*z1 - e1*k2*o2*w1*pow(y,3)*z1 + 
   e1*j2*p2*w1*pow(y,3)*z1 - m2*p2*r1*w1*pow(y,3)*z1 + 
   2*k2*r1*r2*w1*pow(y,3)*z1 + m2*o2*v1*w1*pow(y,3)*z1 - 
   j2*r2*v1*w1*pow(y,3)*z1 - j2*r1*v2*w1*pow(y,3)*z1 + 
   k2*pow(r1,2)*w2*pow(y,3)*z1 - j2*r1*v1*w2*pow(y,3)*z1 - 
   2*n1*r0*r1*t2*x*pow(y,3)*z1 - n0*pow(r1,2)*t2*x*pow(y,3)*z1 + 
   2*n1*p2*r1*v0*x*pow(y,3)*z1 - k2*r1*s1*v0*x*pow(y,3)*z1 + 
   2*n1*p2*r0*v1*x*pow(y,3)*z1 + 2*n0*p2*r1*v1*x*pow(y,3)*z1 - 
   k2*r1*s0*v1*x*pow(y,3)*z1 - k2*r0*s1*v1*x*pow(y,3)*z1 - 
   2*n1*o2*v0*v1*x*pow(y,3)*z1 + 2*j2*s1*v0*v1*x*pow(y,3)*z1 - 
   n0*o2*pow(v1,2)*x*pow(y,3)*z1 + j2*s0*pow(v1,2)*x*pow(y,3)*z1 + 
   k2*pow(r1,2)*w0*x*pow(y,3)*z1 - j2*r1*v1*w0*x*pow(y,3)*z1 + 
   2*k2*r0*r1*w1*x*pow(y,3)*z1 - j2*r1*v0*w1*x*pow(y,3)*z1 - 
   j2*r0*v1*w1*x*pow(y,3)*z1 - n1*pow(r1,2)*t2*pow(y,4)*z1 + 
   2*n1*p2*r1*v1*pow(y,4)*z1 - k2*r1*s1*v1*pow(y,4)*z1 - 
   n1*o2*pow(v1,2)*pow(y,4)*z1 + j2*s1*pow(v1,2)*pow(y,4)*z1 + 
   k2*pow(r1,2)*w1*pow(y,4)*z1 - j2*r1*v1*w1*pow(y,4)*z1 + 
   f2*m2*pow(p2,2)*z2 - e2*n2*pow(p2,2)*z2 - f2*k2*p2*r2*z2 + 
   e2*k2*p2*s2*z2 - f2*m2*o2*t2*z2 + e2*n2*o2*t2*z2 + f2*j2*r2*t2*z2 - 
   n2*pow(r2,2)*t2*z2 - e2*j2*s2*t2*z2 + m2*r2*s2*t2*z2 + f2*k2*o2*v2*z2 - 
   f2*j2*p2*v2*z2 + 2*n2*p2*r2*v2*z2 - m2*p2*s2*v2*z2 - k2*r2*s2*v2*z2 - 
   n2*o2*pow(v2,2)*z2 + j2*s2*pow(v2,2)*z2 - e2*k2*o2*w2*z2 + 
   e2*j2*p2*w2*z2 - m2*p2*r2*w2*z2 + k2*pow(r2,2)*w2*z2 + m2*o2*v2*w2*z2 - 
   j2*r2*v2*w2*z2 + f0*m2*pow(p2,2)*x*z2 - e2*n0*pow(p2,2)*x*z2 - 
   e0*n2*pow(p2,2)*x*z2 - f2*k2*p2*r0*x*z2 - f0*k2*p2*r2*x*z2 + 
   e2*k2*p2*s0*x*z2 + e0*k2*p2*s2*x*z2 - f0*m2*o2*t2*x*z2 + 
   e2*n0*o2*t2*x*z2 + e0*n2*o2*t2*x*z2 + f2*j2*r0*t2*x*z2 + 
   f0*j2*r2*t2*x*z2 - 2*n2*r0*r2*t2*x*z2 - n0*pow(r2,2)*t2*x*z2 - 
   e2*j2*s0*t2*x*z2 + m2*r2*s0*t2*x*z2 - e0*j2*s2*t2*x*z2 + 
   m2*r0*s2*t2*x*z2 + f2*k2*o2*v0*x*z2 - f2*j2*p2*v0*x*z2 + 
   2*n2*p2*r2*v0*x*z2 - m2*p2*s2*v0*x*z2 - k2*r2*s2*v0*x*z2 + 
   f0*k2*o2*v2*x*z2 - f0*j2*p2*v2*x*z2 + 2*n2*p2*r0*v2*x*z2 + 
   2*n0*p2*r2*v2*x*z2 - m2*p2*s0*v2*x*z2 - k2*r2*s0*v2*x*z2 - 
   k2*r0*s2*v2*x*z2 - 2*n2*o2*v0*v2*x*z2 + 2*j2*s2*v0*v2*x*z2 - 
   n0*o2*pow(v2,2)*x*z2 + j2*s0*pow(v2,2)*x*z2 - e2*k2*o2*w0*x*z2 + 
   e2*j2*p2*w0*x*z2 - m2*p2*r2*w0*x*z2 + k2*pow(r2,2)*w0*x*z2 + 
   m2*o2*v2*w0*x*z2 - j2*r2*v2*w0*x*z2 - e0*k2*o2*w2*x*z2 + 
   e0*j2*p2*w2*x*z2 - m2*p2*r0*w2*x*z2 + 2*k2*r0*r2*w2*x*z2 + 
   m2*o2*v0*w2*x*z2 - j2*r2*v0*w2*x*z2 - j2*r0*v2*w2*x*z2 - 
   e0*n0*pow(p2,2)*pow(x,2)*z2 - f0*k2*p2*r0*pow(x,2)*z2 + 
   e0*k2*p2*s0*pow(x,2)*z2 + e0*n0*o2*t2*pow(x,2)*z2 + 
   f0*j2*r0*t2*pow(x,2)*z2 - n2*pow(r0,2)*t2*pow(x,2)*z2 - 
   2*n0*r0*r2*t2*pow(x,2)*z2 - e0*j2*s0*t2*pow(x,2)*z2 + 
   m2*r0*s0*t2*pow(x,2)*z2 + f0*k2*o2*v0*pow(x,2)*z2 - 
   f0*j2*p2*v0*pow(x,2)*z2 + 2*n2*p2*r0*v0*pow(x,2)*z2 + 
   2*n0*p2*r2*v0*pow(x,2)*z2 - m2*p2*s0*v0*pow(x,2)*z2 - 
   k2*r2*s0*v0*pow(x,2)*z2 - k2*r0*s2*v0*pow(x,2)*z2 - 
   n2*o2*pow(v0,2)*pow(x,2)*z2 + j2*s2*pow(v0,2)*pow(x,2)*z2 + 
   2*n0*p2*r0*v2*pow(x,2)*z2 - k2*r0*s0*v2*pow(x,2)*z2 - 
   2*n0*o2*v0*v2*pow(x,2)*z2 + 2*j2*s0*v0*v2*pow(x,2)*z2 - 
   e0*k2*o2*w0*pow(x,2)*z2 + e0*j2*p2*w0*pow(x,2)*z2 - 
   m2*p2*r0*w0*pow(x,2)*z2 + 2*k2*r0*r2*w0*pow(x,2)*z2 + 
   m2*o2*v0*w0*pow(x,2)*z2 - j2*r2*v0*w0*pow(x,2)*z2 - 
   j2*r0*v2*w0*pow(x,2)*z2 + k2*pow(r0,2)*w2*pow(x,2)*z2 - 
   j2*r0*v0*w2*pow(x,2)*z2 - n0*pow(r0,2)*t2*pow(x,3)*z2 + 
   2*n0*p2*r0*v0*pow(x,3)*z2 - k2*r0*s0*v0*pow(x,3)*z2 - 
   n0*o2*pow(v0,2)*pow(x,3)*z2 + j2*s0*pow(v0,2)*pow(x,3)*z2 + 
   k2*pow(r0,2)*w0*pow(x,3)*z2 - j2*r0*v0*w0*pow(x,3)*z2 + 
   f1*m2*pow(p2,2)*y*z2 - e2*n1*pow(p2,2)*y*z2 - e1*n2*pow(p2,2)*y*z2 - 
   f2*k2*p2*r1*y*z2 - f1*k2*p2*r2*y*z2 + e2*k2*p2*s1*y*z2 + 
   e1*k2*p2*s2*y*z2 - f1*m2*o2*t2*y*z2 + e2*n1*o2*t2*y*z2 + 
   e1*n2*o2*t2*y*z2 + f2*j2*r1*t2*y*z2 + f1*j2*r2*t2*y*z2 - 
   2*n2*r1*r2*t2*y*z2 - n1*pow(r2,2)*t2*y*z2 - e2*j2*s1*t2*y*z2 + 
   m2*r2*s1*t2*y*z2 - e1*j2*s2*t2*y*z2 + m2*r1*s2*t2*y*z2 + 
   f2*k2*o2*v1*y*z2 - f2*j2*p2*v1*y*z2 + 2*n2*p2*r2*v1*y*z2 - 
   m2*p2*s2*v1*y*z2 - k2*r2*s2*v1*y*z2 + f1*k2*o2*v2*y*z2 - 
   f1*j2*p2*v2*y*z2 + 2*n2*p2*r1*v2*y*z2 + 2*n1*p2*r2*v2*y*z2 - 
   m2*p2*s1*v2*y*z2 - k2*r2*s1*v2*y*z2 - k2*r1*s2*v2*y*z2 - 
   2*n2*o2*v1*v2*y*z2 + 2*j2*s2*v1*v2*y*z2 - n1*o2*pow(v2,2)*y*z2 + 
   j2*s1*pow(v2,2)*y*z2 - e2*k2*o2*w1*y*z2 + e2*j2*p2*w1*y*z2 - 
   m2*p2*r2*w1*y*z2 + k2*pow(r2,2)*w1*y*z2 + m2*o2*v2*w1*y*z2 - 
   j2*r2*v2*w1*y*z2 - e1*k2*o2*w2*y*z2 + e1*j2*p2*w2*y*z2 - 
   m2*p2*r1*w2*y*z2 + 2*k2*r1*r2*w2*y*z2 + m2*o2*v1*w2*y*z2 - 
   j2*r2*v1*w2*y*z2 - j2*r1*v2*w2*y*z2 - e1*n0*pow(p2,2)*x*y*z2 - 
   e0*n1*pow(p2,2)*x*y*z2 - f1*k2*p2*r0*x*y*z2 - f0*k2*p2*r1*x*y*z2 + 
   e1*k2*p2*s0*x*y*z2 + e0*k2*p2*s1*x*y*z2 + e1*n0*o2*t2*x*y*z2 + 
   e0*n1*o2*t2*x*y*z2 + f1*j2*r0*t2*x*y*z2 + f0*j2*r1*t2*x*y*z2 - 
   2*n2*r0*r1*t2*x*y*z2 - 2*n1*r0*r2*t2*x*y*z2 - 2*n0*r1*r2*t2*x*y*z2 - 
   e1*j2*s0*t2*x*y*z2 + m2*r1*s0*t2*x*y*z2 - e0*j2*s1*t2*x*y*z2 + 
   m2*r0*s1*t2*x*y*z2 + f1*k2*o2*v0*x*y*z2 - f1*j2*p2*v0*x*y*z2 + 
   2*n2*p2*r1*v0*x*y*z2 + 2*n1*p2*r2*v0*x*y*z2 - m2*p2*s1*v0*x*y*z2 - 
   k2*r2*s1*v0*x*y*z2 - k2*r1*s2*v0*x*y*z2 + f0*k2*o2*v1*x*y*z2 - 
   f0*j2*p2*v1*x*y*z2 + 2*n2*p2*r0*v1*x*y*z2 + 2*n0*p2*r2*v1*x*y*z2 - 
   m2*p2*s0*v1*x*y*z2 - k2*r2*s0*v1*x*y*z2 - k2*r0*s2*v1*x*y*z2 - 
   2*n2*o2*v0*v1*x*y*z2 + 2*j2*s2*v0*v1*x*y*z2 + 2*n1*p2*r0*v2*x*y*z2 + 
   2*n0*p2*r1*v2*x*y*z2 - k2*r1*s0*v2*x*y*z2 - k2*r0*s1*v2*x*y*z2 - 
   2*n1*o2*v0*v2*x*y*z2 + 2*j2*s1*v0*v2*x*y*z2 - 2*n0*o2*v1*v2*x*y*z2 + 
   2*j2*s0*v1*v2*x*y*z2 - e1*k2*o2*w0*x*y*z2 + e1*j2*p2*w0*x*y*z2 - 
   m2*p2*r1*w0*x*y*z2 + 2*k2*r1*r2*w0*x*y*z2 + m2*o2*v1*w0*x*y*z2 - 
   j2*r2*v1*w0*x*y*z2 - j2*r1*v2*w0*x*y*z2 - e0*k2*o2*w1*x*y*z2 + 
   e0*j2*p2*w1*x*y*z2 - m2*p2*r0*w1*x*y*z2 + 2*k2*r0*r2*w1*x*y*z2 + 
   m2*o2*v0*w1*x*y*z2 - j2*r2*v0*w1*x*y*z2 - j2*r0*v2*w1*x*y*z2 + 
   2*k2*r0*r1*w2*x*y*z2 - j2*r1*v0*w2*x*y*z2 - j2*r0*v1*w2*x*y*z2 - 
   n1*pow(r0,2)*t2*pow(x,2)*y*z2 - 2*n0*r0*r1*t2*pow(x,2)*y*z2 + 
   2*n1*p2*r0*v0*pow(x,2)*y*z2 + 2*n0*p2*r1*v0*pow(x,2)*y*z2 - 
   k2*r1*s0*v0*pow(x,2)*y*z2 - k2*r0*s1*v0*pow(x,2)*y*z2 - 
   n1*o2*pow(v0,2)*pow(x,2)*y*z2 + j2*s1*pow(v0,2)*pow(x,2)*y*z2 + 
   2*n0*p2*r0*v1*pow(x,2)*y*z2 - k2*r0*s0*v1*pow(x,2)*y*z2 - 
   2*n0*o2*v0*v1*pow(x,2)*y*z2 + 2*j2*s0*v0*v1*pow(x,2)*y*z2 + 
   2*k2*r0*r1*w0*pow(x,2)*y*z2 - j2*r1*v0*w0*pow(x,2)*y*z2 - 
   j2*r0*v1*w0*pow(x,2)*y*z2 + k2*pow(r0,2)*w1*pow(x,2)*y*z2 - 
   j2*r0*v0*w1*pow(x,2)*y*z2 - e1*n1*pow(p2,2)*pow(y,2)*z2 - 
   f1*k2*p2*r1*pow(y,2)*z2 + e1*k2*p2*s1*pow(y,2)*z2 + 
   e1*n1*o2*t2*pow(y,2)*z2 + f1*j2*r1*t2*pow(y,2)*z2 - 
   n2*pow(r1,2)*t2*pow(y,2)*z2 - 2*n1*r1*r2*t2*pow(y,2)*z2 - 
   e1*j2*s1*t2*pow(y,2)*z2 + m2*r1*s1*t2*pow(y,2)*z2 + 
   f1*k2*o2*v1*pow(y,2)*z2 - f1*j2*p2*v1*pow(y,2)*z2 + 
   2*n2*p2*r1*v1*pow(y,2)*z2 + 2*n1*p2*r2*v1*pow(y,2)*z2 - 
   m2*p2*s1*v1*pow(y,2)*z2 - k2*r2*s1*v1*pow(y,2)*z2 - 
   k2*r1*s2*v1*pow(y,2)*z2 - n2*o2*pow(v1,2)*pow(y,2)*z2 + 
   j2*s2*pow(v1,2)*pow(y,2)*z2 + 2*n1*p2*r1*v2*pow(y,2)*z2 - 
   k2*r1*s1*v2*pow(y,2)*z2 - 2*n1*o2*v1*v2*pow(y,2)*z2 + 
   2*j2*s1*v1*v2*pow(y,2)*z2 - e1*k2*o2*w1*pow(y,2)*z2 + 
   e1*j2*p2*w1*pow(y,2)*z2 - m2*p2*r1*w1*pow(y,2)*z2 + 
   2*k2*r1*r2*w1*pow(y,2)*z2 + m2*o2*v1*w1*pow(y,2)*z2 - 
   j2*r2*v1*w1*pow(y,2)*z2 - j2*r1*v2*w1*pow(y,2)*z2 + 
   k2*pow(r1,2)*w2*pow(y,2)*z2 - j2*r1*v1*w2*pow(y,2)*z2 - 
   2*n1*r0*r1*t2*x*pow(y,2)*z2 - n0*pow(r1,2)*t2*x*pow(y,2)*z2 + 
   2*n1*p2*r1*v0*x*pow(y,2)*z2 - k2*r1*s1*v0*x*pow(y,2)*z2 + 
   2*n1*p2*r0*v1*x*pow(y,2)*z2 + 2*n0*p2*r1*v1*x*pow(y,2)*z2 - 
   k2*r1*s0*v1*x*pow(y,2)*z2 - k2*r0*s1*v1*x*pow(y,2)*z2 - 
   2*n1*o2*v0*v1*x*pow(y,2)*z2 + 2*j2*s1*v0*v1*x*pow(y,2)*z2 - 
   n0*o2*pow(v1,2)*x*pow(y,2)*z2 + j2*s0*pow(v1,2)*x*pow(y,2)*z2 + 
   k2*pow(r1,2)*w0*x*pow(y,2)*z2 - j2*r1*v1*w0*x*pow(y,2)*z2 + 
   2*k2*r0*r1*w1*x*pow(y,2)*z2 - j2*r1*v0*w1*x*pow(y,2)*z2 - 
   j2*r0*v1*w1*x*pow(y,2)*z2 - n1*pow(r1,2)*t2*pow(y,3)*z2 + 
   2*n1*p2*r1*v1*pow(y,3)*z2 - k2*r1*s1*v1*pow(y,3)*z2 - 
   n1*o2*pow(v1,2)*pow(y,3)*z2 + j2*s1*pow(v1,2)*pow(y,3)*z2 + 
   k2*pow(r1,2)*w1*pow(y,3)*z2 - j2*r1*v1*w1*pow(y,3)*z2
;

	return fabs(t_5/t_4);
}


