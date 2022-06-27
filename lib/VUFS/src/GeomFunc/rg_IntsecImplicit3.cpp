#include "rg_IntsecImplicit3.h"
#include "rg_RelativeOp.h"

// inserted by Jung-Hyun Ryu
#include "rg_IntersectFunc.h"

//#include <fstream>
#include <time.h>

rg_IntsecImplicit3::rg_IntsecImplicit3()
{

}

rg_IntsecImplicit3::~rg_IntsecImplicit3()
{

}

rg_dList<rg_Point2D> rg_IntsecImplicit3::intersectBzCurveVsBzCurve(const rg_BzCurve2D &curve_s, const rg_BzCurve2D &curve_t, rg_REAL &time)
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
	//////////////////////////////////////////  test

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
//	ofstream fout("timeOut.dat", ios::app);
//	fout << time << endl;
	
 	return intersectPointList;
}  

rg_ImplicitEquation rg_IntsecImplicit3::implicitize(const rg_BzCurve2D &curve)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL x3 = curve.getCtrlPt(3).getX();
	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL y3 = curve.getCtrlPt(3).getY();

//  If
//		    
//	x(t) = a3 t^3 + a2 t^2 + a1 t + a0
//
//			
//	y(t) = b3 t^3 + b2 t^2 + b1 t + b0
//  
//
//  then the resultant of p(x, t) and q(y, t) is like the following.

	rg_REAL a3 =   -x0 + 3*x1 - 3*x2 + x3;
	rg_REAL a2 =  3*x0 - 6*x1 + 3*x2;
	rg_REAL a1 = -3*x0 + 3*x1;
	rg_REAL a0 =    x0;

	rg_REAL b3 =   -y0 + 3*y1 - 3*y2 + y3;
	rg_REAL b2 =  3*y0 - 6*y1 + 3*y2;
	rg_REAL b1 = -3*y0 + 3*y1;
	rg_REAL b0 =    y0;

  //			|															  |
  //			|	m0 x + m1 y + m2	n0 x + n1 y + n2	o0 x + o1 y + o2  |
  //			|															  |
  // R(P, Q) =	|	n0 x + n1 y + n2	p0 x + p1 y + p2	q0 x + q1 y + q2  |
  //			|															  |
  //			|	o0 x + o1 y + o2	q0 x + q1 y + q2	r0 x + r1 y + r2  |
  //			|															  |

	rg_REAL m0 =   0.0;
	rg_REAL m1 =   0.0;
	rg_REAL m2 =   a3*b2 - a2*b3;

  	rg_REAL n0 =   0.0;
	rg_REAL n1 =   0.0;
	rg_REAL n2 =   a3*b1 - a1*b3;

  	rg_REAL o0 =   b3;
	rg_REAL o1 = - a3;
	rg_REAL o2 =   a3*b0 - a0*b3;

   	rg_REAL p0 =   b3;
	rg_REAL p1 = - a3;
	rg_REAL p2 =   a3*b0 + a2*b1 - a1*b2 - a0*b3;

   	rg_REAL q0 =   b2;
	rg_REAL q1 = - a2;
	rg_REAL q2 =   a2*b0 - a0*b2;

	rg_REAL r0 =   b1;
	rg_REAL r1 = - a1;
	rg_REAL r2 =   a1*b0 - a0*b1;

// If the form of the implicitization is 
//        k0     
//      + k1 x   + k2 y
//      + k3 x^2 + k4 xy + k5 y^2
//      + k6 x^3 + k7 x^2 y + k8 x y^2 + k9 y^3
    
    rg_REAL *k = new rg_REAL[10];

	k[0] = -o2*o2*p2 + 2*n2*o2*q2 - m2*q2*q2 - n2*n2*r2 + m2*p2*r2;

    k[1] = - o2*o2*p0 - 2*o0*o2*p2 + 2*n2*o2*q0 + 2*n2*o0*q2 
		   + 2*n0*o2*q2 - 2*m2*q0*q2 - m0*q2*q2 - n2*n2*r0 
		   + m2*p2*r0 - 2*n0*n2*r2 + m2*p0*r2 + m0*p2*r2;
	
    k[2] = - o2*o2*p1 - 2*o1*o2*p2 + 2*n2*o2*q1 + 2*n2*o1*q2 
		   + 2*n1*o2*q2 - 2*m2*q1*q2 - m1*q2*q2 - n2*n2*r1 
		   + m2*p2*r1 - 2*n1*n2*r2	+ m2*p1*r2 + m1*p2*r2;
	
    k[3] = - 2*o0*o2*p0 - o0*o0*p2 + 2*n2*o0*q0 + 2*n0*o2*q0 
		   - m2*q0*q0 + 2*n0*o0*q2 - 2*m0*q0*q2 - 2*n0*n2*r0 
		   + m2*p0*r0 + m0*p2*r0 - n0*n0*r2 + m0*p0*r2;
	
    k[4] = - 2*o1*o2*p0 - 2*o0*o2*p1 - 2*o0*o1*p2 + 2*n2*o1*q0 
		   + 2*n1*o2*q0 + 2*n2*o0*q1 + 2*n0*o2*q1 - 2*m2*q0*q1 
		   + 2*n1*o0*q2 + 2*n0*o1*q2 - 2*m1*q0*q2 - 2*m0*q1*q2 
		   - 2*n1*n2*r0 + m2*p1*r0 + m1*p2*r0 - 2*n0*n2*r1 
		   + m2*p0*r1 + m0*p2*r1 - 2*n0*n1*r2 + m1*p0*r2 + m0*p1*r2;

    k[5] = - 2*o1*o2*p1 - o1*o1*p2 + 2*n2*o1*q1 + 2*n1*o2*q1 
		   - m2*q1*q1 + 2*n1*o1*q2 - 2*m1*q1*q2 - 2*n1*n2*r1 
		   + m2*p1*r1 + m1*p2*r1 - n1*n1*r2 + m1*p1*r2;
	
    k[6] = - o0*o0*p0 + 2*n0*o0*q0 - m0*q0*q0 - n0*n0*r0 
		   + m0*p0*r0;

    k[7] = - 2*o0*o1*p0 - o0*o0*p1 + 2*n1*o0*q0 + 2*n0*o1*q0 
		   - m1*q0*q0 + 2*n0*o0*q1 - 2*m0*q0*q1 - 2*n0*n1*r0 
		   + m1*p0*r0 + m0*p1*r0 - n0*n0*r1 + m0*p0*r1;
	
    k[8] = - o1*o1*p0 - 2*o0*o1*p1 + 2*n1*o1*q0 + 2*n1*o0*q1 
		   + 2*n0*o1*q1 - 2*m1*q0*q1 - m0*q1*q1 - n1*n1*r0 
		   + m1*p1*r0 - 2*n0*n1*r1 + m1*p0*r1 + m0*p1*r1;
	
    k[9] = - o1*o1*p1 + 2*n1*o1*q1 - m1*q1*q1 - n1*n1*r1 + m1*p1*r1;

	return 	rg_ImplicitEquation(3, k);
}

rg_REAL rg_IntsecImplicit3::inversion(const rg_BzCurve2D &curve, const rg_Point2D &point)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL x3 = curve.getCtrlPt(3).getX();
	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL y3 = curve.getCtrlPt(3).getY();

	rg_REAL a3 =   -x0 + 3*x1 - 3*x2 + x3;
	rg_REAL a2 =  3*x0 - 6*x1 + 3*x2;
	rg_REAL a1 = -3*x0 + 3*x1;
	rg_REAL a0 =    x0;

	rg_REAL b3 =   -y0 + 3*y1 - 3*y2 + y3;
	rg_REAL b2 =  3*y0 - 6*y1 + 3*y2;
	rg_REAL b1 = -3*y0 + 3*y1;
	rg_REAL b0 =    y0;

  //			|															  |
  //			|	m0 x + m1 y + m2	n0 x + n1 y + n2	o0 x + o1 y + o2  |
  //			|															  |
  // R(P, Q) =	|	n0 x + n1 y + n2	p0 x + p1 y + p2	q0 x + q1 y + q2  |
  //			|															  |
  //			|	o0 x + o1 y + o2	q0 x + q1 y + q2	r0 x + r1 y + r2  |
  //			|															  |

  	rg_REAL m0 =   0.0;
	rg_REAL m1 =   0.0;
	rg_REAL m2 =   a3*b2 - a2*b3;

  	rg_REAL n0 =   0.0;
	rg_REAL n1 =   0.0;
	rg_REAL n2 =   a3*b1 - a1*b3;

  	rg_REAL o0 =   b3;
	rg_REAL o1 = - a3;
	rg_REAL o2 =   a3*b0 - a0*b3;

   	rg_REAL p0 =   b3;
	rg_REAL p1 = - a3;
	rg_REAL p2 =   a3*b0 + a2*b1 - a1*b2 - a0*b3;

   	rg_REAL q0 =   b2;
	rg_REAL q1 = - a2;
	rg_REAL q2 =   a2*b0 - a0*b2;

	rg_REAL r0 =   b1;
	rg_REAL r1 = - a1;
	rg_REAL r2 =   a1*b0 - a0*b1;

	rg_REAL x  =  point.getX();
	rg_REAL y  =  point.getY();

	rg_REAL t2 = -(o2*p2) + n2*q2 
				- o2*p0*x - o0*p2*x + n2*q0*x + n0*q2*x 
				- o0*p0*x*x + n0*q0*x*x 
				- o2*p1*y - o1*p2*y + n2*q1*y + n1*q2*y 
				- o1*p0*x*y - o0*p1*x*y + n1*q0*x*y + n0*q1*x*y 
				- o1*p1*y*y + n1*q1*y*y;

	rg_REAL t  = -(n2*o2) + m2*q2 
				- n2*o0*x - n0*o2*x + m2*q0*x + m0*q2*x 
				- n0*o0*x*x + m0*q0*x*x 
				- n2*o1*y - n1*o2*y + m2*q1*y + m1*q2*y 
				- n1*o0*x*y - n0*o1*x*y + m1*q0*x*y + m0*q1*x*y 
				- n1*o1*y*y + m1*q1*y*y;

	return fabs(t2/t);
}


