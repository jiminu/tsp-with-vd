#include "rg_IntsecImplicitNotUsingMathematica.h"
#include "rg_RelativeOp.h"

// inserted by Jung-Hyun Ryu
#include "rg_IntersectFunc.h"
#include "rg_Matrix.h"

//#include <fstream>
#include <time.h>
#include "rg_ImplicitEquation.h"
#include "rg_Polynomial.h"


rg_IntsecImplicitNotUsingMathematica::rg_IntsecImplicitNotUsingMathematica()
{

}

rg_IntsecImplicitNotUsingMathematica::~rg_IntsecImplicitNotUsingMathematica()
{

}

rg_dList<rg_Point2D> rg_IntsecImplicitNotUsingMathematica::intersectBzCurveVsBzCurve1(const rg_BzCurve2D &curve_s, const rg_BzCurve2D &curve_t, rg_REAL &time)
{
	rg_dList<rg_Point2D> intersectPointList;
	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	//clock_t StartTime = clock(); 

	//for(rg_INT i=0; i<7; i++)
	//{

		intersectPointList.removeAll();
		clock_t StartTime = clock(); 
		//for(rg_INT o = 0;o < 19;o++) implicitize1(curve_t);
		rg_ImplicitEquation  impEq   = implicitize1(curve_t);
		clock_t EndTime = clock();
		time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
//		std::ofstream fout("timeOutBz2.dat", std::ios::app);
//		fout << time << '\n';

		rg_Polynomial*       cs      = curve_s.convertBzCurve2Polynomial();
		rg_Polynomial        poly    = impEq.substitute(cs[0], cs[1]);
		delete[] cs;

		rg_ComplexNumber*    param_s = poly.solve();

		for(rg_INT i = 0; i < poly.getDegree(); i++)
		{
			if(param_s[i].isPureRealNumber() && rg_BTOR(0.0, param_s[i].getRealNumber(), 1.0) )
			{
				rg_Point2D point_s = curve_s.evaluatePt(param_s[i].getRealNumber());
				//rg_REAL  param_t = inversion(curve_t, point_s); 
				rg_REAL  param_t = inversionUsingLinearEquation(curve_t, point_s);
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


	//clock_t EndTime = clock();
	//time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
//	ofstream fout("timeOut.dat", ios::app);
//	fout << time << endl;
	
 	return intersectPointList;
}


rg_dList<rg_Point2D> rg_IntsecImplicitNotUsingMathematica::intersectBzCurveVsBzCurve2(const rg_BzCurve2D &curve_s, const rg_BzCurve2D &curve_t, rg_REAL &time)
{
	rg_dList<rg_Point2D> intersectPointList;
	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime = clock(); 

	//for(rg_INT i=0; i<3; i++)
	//{

		intersectPointList.removeAll();
		rg_ImplicitEquation  impEq    = implicitize2(curve_t);

		rg_Polynomial*       cs      = curve_s.convertBzCurve2Polynomial();
		rg_Polynomial        poly    = impEq.substitute(cs[0], cs[1]);
		delete[] cs;

		rg_ComplexNumber*    param_s = poly.solve();

		for(rg_INT i = 0; i < poly.getDegree(); i++)
		{
			if(param_s[i].isPureRealNumber() && rg_BTOR(0.0, param_s[i].getRealNumber(), 1.0) )
			{
				rg_Point2D point_s = curve_s.evaluatePt(param_s[i].getRealNumber());
				//rg_REAL  param_t = inversion(curve_t, point_s); 
				rg_REAL  param_t = inversionUsingLinearEquation(curve_t, point_s);
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


rg_ImplicitEquation rg_IntsecImplicitNotUsingMathematica::implicitize1(const rg_BzCurve2D &curve)
{
	return rg_IntersectFunc::implicitizeNotUsingMathematica1(curve);
}

rg_ImplicitEquation rg_IntsecImplicitNotUsingMathematica::implicitize2(const rg_BzCurve2D &curve)
{
	return rg_IntersectFunc::implicitizeNotUsingMathematica2(curve);
}

rg_REAL rg_IntsecImplicitNotUsingMathematica::inversion(const rg_BzCurve2D &curve, const rg_Point2D &point)
{

// Using determinant or gaussian elimination 
/*
	rg_DEGREE degree = curve.getDegree();
	rg_REAL resultantValue = 0.0;
	rg_Matrix numerator(degree - 1, degree - 1);
	rg_Matrix denominator(degree - 1, degree - 1);
*/

// Using the GCD(Greatest Common Divisor)

	rg_DEGREE degree = curve.getDegree();
	rg_Polynomial* poly = curve.convertBzCurve2Polynomial();

	rg_Polynomial temp_x(0), temp_y(0);

	temp_x.setCoefficient(0, point.getX());
	rg_Polynomial x_t = temp_x - poly[0];

	temp_y.setCoefficient(0, point.getY());
	rg_Polynomial y_t = temp_y - poly[1];

	delete[] poly;

	// initialization
	rg_Polynomial dividend = (x_t.getCoefficient(degree) > y_t.getCoefficient(degree) ? x_t : y_t);
	rg_Polynomial divisor = (dividend == x_t) ? y_t : x_t;
	rg_Polynomial times = dividend / divisor;
	rg_Polynomial preRemainder, nextRemainder;
	preRemainder = nextRemainder = dividend - divisor * times;

	while(!( (nextRemainder.getDegree() == 0) ) ) // && ( rg_BTOR(0.0, nextRemainder.getCoefficient(0), 1.0) ) ) )
	{
		preRemainder = nextRemainder;
		dividend = divisor;
		divisor = preRemainder;

		times = dividend / divisor;
		nextRemainder = dividend - divisor * times;
		nextRemainder.reconstruct();
	}
	
	rg_ComplexNumber* param_t = preRemainder.solve();

	rg_REAL parameter_t;

	for(rg_INDEX i = 0 ; i < preRemainder.getDegree() ; i++)
	{
		if(param_t[i].isPureRealNumber() && rg_BTOR(0.0, param_t[i].getRealNumber(), 1.0))
		{
			parameter_t = param_t[i].getRealNumber();
		}
	}

	return parameter_t;
}

rg_REAL rg_IntsecImplicitNotUsingMathematica::inversionUsingLinearEquation(const rg_BzCurve2D &curve, const rg_Point2D &point)
{
	rg_DEGREE degree = curve.getDegree();
	rg_Polynomial* polyForm = curve.convertBzCurve2Polynomial();
	rg_Polynomial x(0), y(0);
	x.setCoefficient(0, point.getX());
	y.setCoefficient(0, point.getY());

	polyForm[ 0 ] = polyForm[ 0 ] - x;
	polyForm[ 1 ] = polyForm[ 1 ] - y;
	
	rg_Matrix resultantMatrix(degree, degree);

	// fills elements in coefficient matrix necessary for calculating resultant

	for(rg_INT i = 0;i < degree;i++)
	{
		for(rg_INT j = 0;j < degree;j++)
		{
			rg_INT p = ((degree - i > degree - j) ? degree - i : degree - j);
			rg_INT q = 2*degree - i - j - 1 - p;

			while(p <= degree && q >= 0)
			{
				rg_Matrix temp(2, 2);
				temp[0][0] = polyForm[ 0 ].getCoefficient(p);
				temp[0][1] = polyForm[ 0 ].getCoefficient(q);
				temp[1][0] = polyForm[ 1 ].getCoefficient(p);
				temp[1][1] = polyForm[ 1 ].getCoefficient(q);
				resultantMatrix[ i ][ j ] = resultantMatrix[ i ][ j ] + temp.determinant();
				p++;
			    q--;
			}
		}
	}

	delete[] polyForm;

	// delete the last row
	rg_Matrix lastRowDeletedresultantMatrix = resultantMatrix.getOneRowDeletedSubmatrix(degree - 1);

	rg_Matrix lastRowAndOneColDeletedresultantMatrix1 = lastRowDeletedresultantMatrix.getOneColDeletedSubmatrix( 0 );
	rg_Matrix lastRowAndOneColDeletedresultantMatrix2 = lastRowDeletedresultantMatrix.getOneColDeletedSubmatrix( 1 );

	rg_REAL t_degree_1 = lastRowAndOneColDeletedresultantMatrix1.determinant();
	rg_REAL t_degree_2 = lastRowAndOneColDeletedresultantMatrix2.determinant();

	return fabs(t_degree_1 / t_degree_2);
}


