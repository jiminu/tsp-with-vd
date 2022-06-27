/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : intersectFunc.cpp
//	  
//    DESCRIPTION : collections of function related with intersection
//                          
//	  CLASS NAME  : 
//
//    BASE CLASS  : 
//      
//    AUTHOR      : Lee, Soon-Woong
//
//    START DATE  :     
//    
//    HISTORY :     1. Ryu, Jung-Hyun inserted the following functions on Feb. 25. 1998
//
//                       rg_ImplicitEquation rg_IntersectFunc::getResultantByBezout(rg_Polynomial x_t, rg_Polynomial y_t);
//                       rg_ImplicitEquation rg_IntersectFunc::getResultantByBezout(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t);
//					     rg_ImplicitEquation implicitizeNotUsingMathematica(const rg_RQBzCurve2D &curve);
//                       rg_ImplicitEquation implicitizeNotUsingMathematica(const rg_BzCurve2D& curve);
//
//                  2. Ryu, Jung-Hyun inserted the following functions on Mar. 5. 1998
//
//                       rg_ImplicitEquation rg_IntersectFunc::getResultantBySylvester(rg_Polynomial x_t, rg_Polynomial y_t);
//                       rg_ImplicitEquation rg_IntersectFunc::getResultantBySylvester(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t);
//					 
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include "rg_RelativeOp.h"
#include "rg_IntersectFunc.h"
#include "rg_MathFunc.h"
#include "rg_QuarticPolynomial.h"

// Inserted by Jung-Hyun Ryu
#include "rg_Polynomial.h"
#include "rg_ImplicitEqMatrix.h"
#include "sortFunc.h"

rg_ImplicitEquation rg_IntersectFunc::implicitize(const rg_RQBzCurve2D &curve)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL w0 = curve.getWeight(0);
	rg_REAL w1 = curve.getWeight(1);
	rg_REAL w2 = curve.getWeight(2);

//  If
//		    a2 t^2 + a1 t^1 + a0
//	x(t) = ----------------------
//          d2 t^2 + d1 t^1 + d0
//
//
//			b2 t^2 + b1 t^1 + b0
//	y(t) = -----------------------------
//	        d2 t^3 + d1 t^1 + d0
//  
//
//  then the resultant of p(x, t) and q(y, t) is like the following.	
    
    rg_REAL a2 =    w0*x0 - 2*w1*x1 + w2*x2;
	rg_REAL a1 = -2*w0*x0 + 2*w1*x1;
	rg_REAL a0 =    w0*x0;

	rg_REAL b2 =    w0*y0 - 2*w1*y1 + w2*y2;
	rg_REAL b1 = -2*w0*y0 + 2*w1*y1;
	rg_REAL b0 =    w0*y0;

	rg_REAL d2 =    w0 - 2*w1 + w2;
	rg_REAL d1 = -2*w0 + 2*w1;
	rg_REAL d0 =    w0;

  //			|										  |
  //			|	m0 x + m1 y + m2	n0 x + n1 y + n2  |
  // R(P, Q) =	|										  |
  //       	    |	n0 x + n1 y + n2	p0 x + p1 y + p2  |
  //			|										  |

	rg_REAL m0 =  b2*d1 - b1*d2;
	rg_REAL m1 = -a2*d1 + a1*d2;
	rg_REAL m2 =  a2*b1 - a1*b2;

	rg_REAL n0 =  b2*d0 - b0*d2;
	rg_REAL n1 = -a2*d0 + a0*d2;
	rg_REAL n2 =  a2*b0 - a0*b2;

	rg_REAL p0 =  b1*d0 - b0*d1;
	rg_REAL p1 = -a1*d0 + a0*d1;
	rg_REAL p2 =  a1*b0 - a0*b1;

    // k0 + k1x + k2y + k3x^2 + k4xy + k5y^2
    rg_REAL *k = new rg_REAL[6];

    k[0] = -  n2*n2 + m2*p2;
    k[1] = -2*n0*n2 + m2*p0 + m0*p2;
    k[2] = -2*n1*n2 + m2*p1 + m1*p2; 
    k[3] = -  n0*n0 + m0*p0;
    k[4] = -2*n0*n1 + m1*p0 + m0*p1;
    k[5] = -  n1*n1 + m1*p1;

    return rg_ImplicitEquation(2, k); 
}

//rg_ImplicitEquation getResultant(const rg_BzCurve2D& curve)
rg_ImplicitEquation rg_IntersectFunc::getResultantByBezout(rg_Polynomial x_t, rg_Polynomial y_t) //test
{
	//rg_DEGREE degree = curve.getDegree();
	rg_DEGREE degree = (x_t.getDegree() > y_t.getDegree()) ? x_t.getDegree() : y_t.getDegree();
	//rg_Polynomial* polyForm = curve.convertBzCurve2Polynomial();
	//rg_Polynomial polyForm[2];
	//polyForm[0] = x_t;
	//polyForm[1] = y_t;

	rg_ImplicitEqMatrix tmatrix(degree, degree);
	rg_ImplicitEquation resultant;

	// definition of resultant

	for(rg_INT i = 0;i < degree;i++)
	{
		for(rg_INT j = 0;j < degree;j++)
		{
			tmatrix[ i ][ j ].setDegree(1); // initialization
			rg_INT p = ((degree - i > degree - j) ? degree - i : degree - j);
			//rg_INT p = ((degree - 1 - i > degree - 1 - j) ? degree - 1 - i : degree - 1 - j);
			rg_INT q = 2*degree - i - j - 1 - p;
			//rg_INT q = 2*(degree - 1) - i - j - 1 - p;

			while(p <= degree && q >= 0)
			{
				rg_ImplicitEqMatrix temp(2, 2);
				rg_ImplicitEquation tempX_t(1);
				rg_ImplicitEquation tempY_t(1);

				// the case that both p and q are equal to zero does not occur.
				if(p == 0)
				{
					tempX_t[0][0] = -x_t.getCoefficient(p);
					tempX_t[1][0] = 1.0;
					temp[0][0] = tempX_t;

					tempX_t[0][0] = -x_t.getCoefficient(q);
					temp[0][1] = tempX_t;

					tempY_t[0][0] = -y_t.getCoefficient(p);
					tempY_t[0][1] = 1.0;
					temp[1][0] = tempY_t;

					tempY_t[0][0] = -y_t.getCoefficient(q);										
					temp[1][1] = tempY_t;
				}
				else if(q == 0)
				{
					tempX_t[0][0] = -x_t.getCoefficient(p);
					temp[0][0] = tempX_t;
					//temp[0][0].show();

					tempX_t[0][0] = -x_t.getCoefficient(q);
					tempX_t[1][0] = 1.0;
					temp[0][1] = tempX_t;
					//temp[0][1].show();

					tempY_t[0][0] = -y_t.getCoefficient(p);
					temp[1][0] = tempY_t;
					//temp[1][0].show();

					tempY_t[0][0] = -y_t.getCoefficient(q);
					tempY_t[0][1] = 1.0;
					temp[1][1] = tempY_t;
					//temp[1][1].show();
				}
				else
				{
					tempX_t[0][0] = -x_t.getCoefficient(p);
					temp[0][0] = tempX_t;
					tempX_t[0][0] = -x_t.getCoefficient(q);
					temp[0][1] = tempX_t;
					tempY_t[0][0] = -y_t.getCoefficient(p);
					temp[1][0] = tempY_t;
					tempY_t[0][0] = -y_t.getCoefficient(q);
					temp[1][1] = tempY_t;
				}

				
				tmatrix[ i ][ j ] += temp.getDeterminant();
				//tmatrix[ i ][ j ].show();//test
				p++;
			    q--;
			}
		}
	}
	
	resultant = tmatrix.getDeterminant();

	//resultant.reconstruct();//test

	return resultant;
}

rg_ImplicitEquation rg_IntersectFunc::getResultantBySylvesterOld(rg_Polynomial x_t, rg_Polynomial y_t)//antique
{
	rg_DEGREE largerDegree, smallerDegree, degree;
	rg_Polynomial larger, smaller;
	rg_FLAG isXlarger = rg_TRUE;// X(t) has the larger degree than Y(t)

	if(x_t.getDegree() > y_t.getDegree())
	{
		largerDegree = x_t.getDegree();
		larger = x_t;
		smallerDegree = y_t.getDegree();
		smaller = y_t;
	}
	else
	{
		largerDegree = y_t.getDegree();
		larger = y_t;
		smallerDegree = x_t.getDegree();
		smaller = x_t;
		isXlarger = rg_FALSE;// Y(t) has the larger degree than Y(t)
	}

	degree = smallerDegree + largerDegree;

	rg_ImplicitEqMatrix tmatrix(degree, degree);

	rg_ImplicitEquation resultant;


	rg_INT i = 0;
	for(i = 0;i < smallerDegree;i++)
	{
		rg_INT index = largerDegree;
		for(rg_INT j = 0;j < degree;j++)
		{
			rg_ImplicitEquation impEqZeroDeg(0);
			rg_ImplicitEquation impEqOneDeg(1);

			if(i <= j && j <= (i + largerDegree))
			{
				if(j == i + largerDegree)
				{
					if(isXlarger)
					{
						impEqOneDeg[ 0 ][ 0 ] = - larger.getCoefficient(index);
						impEqOneDeg[ 1 ][ 0 ] = 1;
					}
					else
					{
						impEqOneDeg[ 0 ][ 0 ] = - larger.getCoefficient(index);
						impEqOneDeg[ 0 ][ 1 ] = 1;
					}
					tmatrix[ i ][ j ] = impEqOneDeg;
					//tmatrix[ i ][ j ].show();//test
				}
				else
				{
					impEqZeroDeg[ 0 ][ 0 ] = - larger.getCoefficient(index);
					tmatrix[ i ][ j ] = impEqZeroDeg;
					//tmatrix[ i ][ j ].show();//test
					index--;
				}
			}			
			else
			{
				impEqZeroDeg[ 0 ][ 0 ] = 0.0;
				tmatrix[ i ][ j ] = impEqZeroDeg;
				//tmatrix[ i ][ j ].show();//test
			}
		}
	}

	for(i = smallerDegree;i < degree;i++)
	{
		rg_INT index = smallerDegree;
		for(rg_INT j = 0;j < degree;j++)
		{
			rg_ImplicitEquation impEqZeroDeg(0);
			rg_ImplicitEquation impEqOneDeg(1);

			if((degree - 1  - i) <= j && j <= (degree - 1  - i + smallerDegree))
			//if(i - smallerDegree <= j && j <= (i + smallerDegree + smallerDegree)) //test
			{
				if(j == degree - 1  - i + smallerDegree)
				//if(j == i)//test
				{
					if(isXlarger)
					{
						impEqOneDeg[ 0 ][ 0 ] = - smaller.getCoefficient(index);
						impEqOneDeg[ 0 ][ 1 ] = 1;
						//impEqOneDeg[ 0 ][ 0 ] = - smaller.getCoefficient(index);
						//impEqOneDeg[ 1 ][ 0 ] = 1;
					}
					else
					{
						//impEqOneDeg[ 0 ][ 0 ] = - smaller.getCoefficient(index);
						//impEqOneDeg[ 0 ][ 1 ] = 1;
						impEqOneDeg[ 0 ][ 0 ] = - smaller.getCoefficient(index);
						impEqOneDeg[ 1 ][ 0 ] = 1;
					}
					tmatrix[ i ][ j ] = impEqOneDeg;
					//tmatrix[ i ][ j ].show();//test
				}
				else
				{
					impEqZeroDeg[ 0 ][ 0 ] = - smaller.getCoefficient(index);
					tmatrix[ i ][ j ] = impEqZeroDeg;
					//tmatrix[ i ][ j ].show();//test
					index--;
				}
			}
			else
			{
				impEqZeroDeg[ 0 ][ 0 ] = 0.0;
				tmatrix[ i ][ j ] = impEqZeroDeg;
				//tmatrix[ i ][ j ].show();//test
			}
		}
	}

	resultant = tmatrix.getDeterminant();
	resultant.reconstruct();//test

	return resultant;	
}

rg_ImplicitEquation rg_IntersectFunc::getResultantBySylvester(rg_Polynomial x_t, rg_Polynomial y_t)
{
	rg_DEGREE largerDegree, smallerDegree, degree;
	rg_Polynomial larger, smaller;
	rg_FLAG isXlarger = rg_TRUE;// X(t) has the larger degree than Y(t)

	if(x_t.getDegree() > y_t.getDegree())
	{
		largerDegree = x_t.getDegree();
		larger = x_t;
		smallerDegree = y_t.getDegree();
		smaller = y_t;
	}
	else
	{
		largerDegree = y_t.getDegree();
		larger = y_t;
		smallerDegree = x_t.getDegree();
		smaller = x_t;
		isXlarger = rg_FALSE;// Y(t) has the larger degree than X(t)
	}

	degree = smallerDegree + largerDegree;

	rg_ImplicitEqMatrix tmatrix(degree, degree);

	rg_ImplicitEquation resultant;

	rg_INT i = 0;
	for(i = 0;i < smallerDegree;i++)
	{
		rg_ImplicitEquation impEqOneDeg(1);
		rg_ImplicitEquation impEqZeroDeg(0);
		rg_Polynomial temp(1);
		temp[ 1 ] = 1.0;

		if(isXlarger)
		{
			impEqOneDeg[ 1 ][ 0 ] = 1.0;			
		}
		else
		{
			impEqOneDeg[ 0 ][ 1 ] = 1.0;			
		}

		//impEqOneDeg.show();//test

		temp = larger * (temp^(smallerDegree - 1 - i));
		
		//temp.reconstruct();

		//temp.show();//test

		rg_INT index = largerDegree + smallerDegree - 1;
		for(rg_INT j = 0;j < degree;j++)
		{
			if(j == (i + largerDegree))
			{
				impEqOneDeg[ 0 ][ 0 ] = - temp.getCoefficient( index );
				tmatrix[ i ][ j ] = impEqOneDeg;
				index--;
			}
			else
			{
				impEqZeroDeg[ 0 ][ 0 ] = - temp.getCoefficient( index );
				tmatrix[ i ][ j ] = impEqZeroDeg;
				index--;
			}
			//tmatrix[ i ][ j ].show();//test
		}
	}

	for(i = smallerDegree;i < degree;i++)
	{
		rg_ImplicitEquation impEqOneDeg(1);
		rg_ImplicitEquation impEqZeroDeg(0);
		rg_Polynomial temp(1);
		temp[ 1 ] = 1.0;

		if(isXlarger)
		{
			impEqOneDeg[ 0 ][ 1 ] = 1.0;
		}
		else
		{
			impEqOneDeg[ 1 ][ 0 ] = 1.0;
		}
		
		//impEqOneDeg.show();//test

		temp = smaller * (temp^(largerDegree + smallerDegree - 1 - i));

		//temp.reconstruct();

		//temp.show();//test

		rg_INT index = largerDegree + smallerDegree - 1;
		for(rg_INT j = 0;j < degree;j++)
		{
			if(j == ((i - smallerDegree) + smallerDegree))
			{
				impEqOneDeg[ 0 ][ 0 ] = - temp.getCoefficient( index );
				tmatrix[ i ][ j ] = impEqOneDeg;
				index--;
			}
			else
			{
				impEqZeroDeg[ 0 ][ 0 ] = - temp.getCoefficient( index );
				tmatrix[ i ][ j ] = impEqZeroDeg;
				index--;
			}
			//tmatrix[ i ][ j ].show();//test
		}
	}

	resultant = tmatrix.getDeterminant();
	//resultant.reconstruct();//test

	return resultant;
}

rg_ImplicitEquation rg_IntersectFunc::getResultantByBezout(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t)
{
	rg_DEGREE degree = (x_t.getDegree() > y_t.getDegree()) ? x_t.getDegree() : y_t.getDegree();
	degree = (degree > w_t.getDegree()) ? degree : w_t.getDegree();

	rg_ImplicitEqMatrix tmatrix(degree, degree);
	rg_ImplicitEqMatrix temp(2, 2);
	rg_ImplicitEquation tempX_t(1);
	rg_ImplicitEquation tempY_t(1);
	rg_ImplicitEquation resultant;
	
	for(rg_INT i = 0;i < degree;i++)
	{
		for(rg_INT j = 0;j < degree;j++)
		{
			tmatrix[ i ][ j ].setDegree(1); // initialization
			rg_INT p = ((degree - i > degree - j) ? degree - i : degree - j);
			rg_INT q = 2*degree - i - j - 1 - p;
			while(p <= degree && q >= 0)
			{
				tempX_t[0][0] = -x_t.getCoefficient(p);
				tempX_t[1][0] = w_t.getCoefficient(p);
				//w_t.show();//w_t.show();//test
				temp[0][0] = tempX_t;
				//temp[0][0].show();//test

				tempX_t[0][0] = -x_t.getCoefficient(q);
				tempX_t[1][0] = w_t.getCoefficient(q);
				temp[0][1] = tempX_t;
				//temp[0][1].show();//test

				tempY_t[0][0] = -y_t.getCoefficient(p);
				tempY_t[0][1] = w_t.getCoefficient(p);
				temp[1][0] = tempY_t;
				//temp[1][0].show();//test

				tempY_t[0][0] = -y_t.getCoefficient(q);
				tempY_t[0][1] = w_t.getCoefficient(q);
				temp[1][1] = tempY_t;
				//temp[1][1].show();//test
				
				tmatrix[ i ][ j ] += temp.getDeterminant();
				//tmatrix[ i ][ j ].show();//test
				p++;
			    q--;
			}
		}
	}
	
	resultant = tmatrix.getDeterminant();
	//resultant.reconstruct();//test
	
	return resultant;
}

rg_ImplicitEquation rg_IntersectFunc::getResultantBySylvesterOld(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t)//antique
{
	rg_DEGREE degreeArray[ 3 ];
	degreeArray[ 0 ] = x_t.getDegree();
	degreeArray[ 1 ] = y_t.getDegree();
	degreeArray[ 2 ] = w_t.getDegree();
	QuickSort(degreeArray, 0, 2);
	
	rg_DEGREE largerDegree, smallerDegree, degree;
	rg_Polynomial larger, smaller;
	rg_FLAG isXlarger = rg_TRUE;// X(t) has the larger degree than Y(t)

	if(w_t.getDegree() == degreeArray[ 2 ]) // if w_t has the largest degree, ...
	{
		largerDegree = smallerDegree = degreeArray[ 2 ];
		larger = x_t;
		smaller = y_t;
	}
	else // if x_t or y_t has the largest degree, ...
	{
		if(x_t.getDegree() > y_t.getDegree())
		{
			largerDegree = x_t.getDegree();
			larger = x_t;
			smallerDegree = y_t.getDegree();
			smaller = y_t;
		}
		else
		{
			largerDegree = y_t.getDegree();
			larger = y_t;
			smallerDegree = x_t.getDegree();
			smaller = x_t;
			isXlarger = rg_FALSE;// Y(t) has the larger degree than Y(t)
		}
	}

	degree = smallerDegree + largerDegree;

	rg_ImplicitEqMatrix tmatrix(degree, degree);

	rg_ImplicitEquation resultant;


	rg_INT i = 0;
	for(i = 0;i < smallerDegree;i++)
	{
		rg_INT index = largerDegree;
		for(rg_INT j = 0;j < degree;j++)
		{
			rg_ImplicitEquation impEqZeroDeg(0);
			rg_ImplicitEquation impEqOneDeg(1);

			if(i <= j && j <= (i + largerDegree))
			{
    			if(isXlarger)
				{
					impEqOneDeg[ 0 ][ 0 ] = - larger.getCoefficient(index);
					impEqOneDeg[ 1 ][ 0 ] = w_t.getCoefficient(index);
				}
				else
				{
					impEqOneDeg[ 0 ][ 0 ] = - larger.getCoefficient(index);
					impEqOneDeg[ 0 ][ 1 ] = w_t.getCoefficient(index);
				}
				tmatrix[ i ][ j ] = impEqOneDeg;
				//tmatrix[ i ][ j ].show();//test
				index--;
			}			
			else
			{
				impEqZeroDeg[ 0 ][ 0 ] = 0.0;
				tmatrix[ i ][ j ] = impEqZeroDeg;
				//tmatrix[ i ][ j ].show();//test
			}
		}
	}

	for(i = smallerDegree;i < degree;i++)
	{
		rg_INT index = smallerDegree;
		for(rg_INT j = 0;j < degree;j++)
		{
			rg_ImplicitEquation impEqZeroDeg(0);
			rg_ImplicitEquation impEqOneDeg(1);

			if((degree - 1  - i) <= j && j <= (degree - 1  - i + smallerDegree))
			{
				if(isXlarger)
				{
					impEqOneDeg[ 0 ][ 0 ] = - smaller.getCoefficient(index);
					impEqOneDeg[ 0 ][ 1 ] = w_t.getCoefficient(index);
					//impEqOneDeg[ 0 ][ 0 ] = - smaller.getCoefficient(index);
					//impEqOneDeg[ 1 ][ 0 ] = 1;
				}
				else
				{
					//impEqOneDeg[ 0 ][ 0 ] = - smaller.getCoefficient(index);
					//impEqOneDeg[ 0 ][ 1 ] = 1;
    				impEqOneDeg[ 0 ][ 0 ] = - smaller.getCoefficient(index);
					impEqOneDeg[ 1 ][ 0 ] = w_t.getCoefficient(index);
				}
				tmatrix[ i ][ j ] = impEqOneDeg;
				//tmatrix[ i ][ j ].show();//test
				index--;
			}
			else
			{
				impEqZeroDeg[ 0 ][ 0 ] = 0.0;
				tmatrix[ i ][ j ] = impEqZeroDeg;
				//tmatrix[ i ][ j ].show();//test
			}
		}
	}

	resultant = tmatrix.getDeterminant();
	resultant.reconstruct();//test

	return resultant;	
}

rg_ImplicitEquation rg_IntersectFunc::getResultantBySylvester(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t)
{
	rg_DEGREE degreeArray[ 3 ];
	degreeArray[ 0 ] = x_t.getDegree();
	degreeArray[ 1 ] = y_t.getDegree();
	degreeArray[ 2 ] = w_t.getDegree();
	QuickSort(degreeArray, 0, 2);
	
	rg_DEGREE largerDegree, smallerDegree, degree;
	rg_Polynomial larger, smaller;
	rg_FLAG isXlarger = rg_TRUE;// X(t) has the larger degree than Y(t)
	//rg_FLAG isWlargest = rg_FALSE;

	if(w_t.getDegree() == degreeArray[ 2 ]) // if w_t has the largest degree, ...
	{
		largerDegree = smallerDegree = degreeArray[ 2 ];
		larger = x_t;
		smaller = y_t;
		//isWlargest = rg_TRUE;
	}
	else // if x_t or y_t has the largest degree, ...
	{
		if(x_t.getDegree() > y_t.getDegree())
		{
			largerDegree = x_t.getDegree();
			larger = x_t;
			smallerDegree = y_t.getDegree();
			smaller = y_t;
		}
		else
		{
			largerDegree = y_t.getDegree();
			larger = y_t;
			smallerDegree = x_t.getDegree();
			smaller = x_t;
			isXlarger = rg_FALSE;// Y(t) has the larger degree than Y(t)
		}
	}

	degree = smallerDegree + largerDegree;

	rg_ImplicitEqMatrix tmatrix(degree, degree);

	rg_ImplicitEquation resultant;

	rg_INT i = 0;
	for(i = 0;i < smallerDegree;i++)
	{
		rg_ImplicitEquation impEqOneDeg(1), impEqZeroDeg(0);
		rg_Polynomial temp(1), temp2(1);
		temp[ 1 ] = 1.0;
		temp2[ 1 ] = 1.0;

		//if(isWlargest)
		//{
			temp2 = w_t * (temp2^(smallerDegree - 1 - i));
		//}
		temp = larger * (temp^(smallerDegree - 1 - i));
		//temp.show();//test
		//temp2.show();//test

		rg_INT index = largerDegree + smallerDegree - 1;
		for(rg_INT j = 0;j < degree;j++)
		{
			/*
			if(j == (i + largerDegree))
			{
				impEqZeroDeg[ 0 ][ 0 ] = temp2.getCoefficient( index ) - temp.getCoefficient( index );
				tmatrix[ i ][ j ] = impEqZeroDeg;
				index--;
			}
			else
			{*/
				impEqOneDeg[ 0 ][ 0 ] = - temp.getCoefficient( index );
				if(isXlarger)
				{
					impEqOneDeg[ 1 ][ 0 ] = temp2.getCoefficient( index );
				}
				else
				{
					impEqOneDeg[ 0 ][ 1 ] = temp2.getCoefficient( index );
				}
				tmatrix[ i ][ j ] = impEqOneDeg;
				index--;
			//}
			//tmatrix[ i ][ j ].show();//test
		}
	}

	for(i = smallerDegree;i < degree;i++)
	{
		rg_ImplicitEquation impEqOneDeg(1), impEqZeroDeg(0);
		rg_Polynomial temp(1), temp2(1);
		temp[ 1 ] = 1.0;
		temp2[ 1 ] = 1.0;

		//if(isWlargest)
		//{
			temp2 = w_t * (temp^(largerDegree + smallerDegree - 1 - i));
		//}
		temp = smaller * (temp^(largerDegree + smallerDegree - 1 - i));
		//temp.show();//test
		//temp2.show();//test

		rg_INT index = largerDegree + smallerDegree - 1;
		for(rg_INT j = 0;j < degree;j++)
		{
			/*
			if(j == ((i - smallerDegree) + smallerDegree))
			{
				impEqZeroDeg[ 0 ][ 0 ] = temp2.getCoefficient( index ) - temp.getCoefficient( index );
				tmatrix[ i ][ j ] = impEqZeroDeg;
				index--;
			}
			else
			{*/
				impEqOneDeg[ 0 ][ 0 ] = - temp.getCoefficient( index );
				if(isXlarger)
				{
					impEqOneDeg[ 0 ][ 1 ] = temp2.getCoefficient( index );
				}
				else
				{
					impEqOneDeg[ 1 ][ 0 ] = temp2.getCoefficient( index );
				}
				tmatrix[ i ][ j ] = impEqOneDeg;
				index--;
			//}
			//tmatrix[ i ][ j ].show();//test
		}
	}

	resultant = tmatrix.getDeterminant();
	resultant.reconstruct();//test

	return resultant;
}

rg_ImplicitEquation rg_IntersectFunc::implicitizeNotUsingMathematica1(const rg_RQBzCurve2D &curve)
{
	rg_Polynomial* numerator = curve.convertNumerator2Polynomial();
	rg_Polynomial x_t = numerator[0];
	rg_Polynomial y_t = numerator[1];
	rg_Polynomial w_t = curve.convertDenominator2Polynomial();

	delete[] numerator;

	return rg_IntersectFunc::getResultantByBezout(x_t, y_t, w_t);	
}

rg_ImplicitEquation rg_IntersectFunc::implicitizeNotUsingMathematica2(const rg_RQBzCurve2D &curve)
{
	rg_Polynomial* numerator = curve.convertNumerator2Polynomial();
	rg_Polynomial x_t = numerator[0];
	rg_Polynomial y_t = numerator[1];
	rg_Polynomial w_t = curve.convertDenominator2Polynomial();

	delete[] numerator;

	//return rg_IntersectFunc::getResultantBySylvesterOld(x_t, y_t, w_t); //test
	return rg_IntersectFunc::getResultantBySylvester(x_t, y_t, w_t); //test
}

rg_ImplicitEquation rg_IntersectFunc::implicitizeNotUsingMathematica1(const rg_BzCurve2D& curve)
{
	rg_Polynomial* polyForm = curve.convertBzCurve2Polynomial();
	rg_Polynomial x_t = polyForm[0];
	rg_Polynomial y_t = polyForm[1];
	
	delete[] polyForm;

	return rg_IntersectFunc::getResultantByBezout(x_t, y_t);
}

rg_ImplicitEquation rg_IntersectFunc::implicitizeNotUsingMathematica2(const rg_BzCurve2D& curve)
{
	rg_Polynomial* polyForm = curve.convertBzCurve2Polynomial();
	rg_Polynomial x_t = polyForm[0];
	rg_Polynomial y_t = polyForm[1];
	
	delete[] polyForm;

    //return rg_IntersectFunc::getResultantBySylvesterOld(x_t, y_t);//test
	return rg_IntersectFunc::getResultantBySylvester(x_t, y_t);//test
}


/////////////////////////////////////////////////////////////////////////////////////////

// Evaluating resultant by eigen decomposition

// Given polynomial f(x), g(x),
//								 f(s)g(t) - g(t)f(s)
// extract numerica value from  ---------------------
//								        s - t
// This function is called during evaluating resultant by eigen decomposition


rg_Matrix rg_IntersectFunc::getNumericValueFromOperationBtnTwoDifferentPolynomial(const rg_Polynomial& f, const rg_Polynomial& g)
{
	rg_DEGREE m = (f.getDegree() > g.getDegree()) ? f.getDegree() : g.getDegree();

	rg_Matrix numericValuedMatrix(m, m);

	for(rg_INT i = 0;i <= m;i++)
	{
		rg_INT j = 0;
		for(j = 0;j <= i;j++)
		{
			rg_REAL numericValue = f.getCoefficient( i ) * g.getCoefficient( j );
			rg_INT bound = i - j - 1;
			
			for(rg_INT k = 0;k <= bound;k++)
			{
				numericValuedMatrix[i - 1 - k][k + j] += numericValue;
			}
		}

		for(j = i + 1;j <= m;j++)
		{
			rg_REAL numericValue = - f.getCoefficient( i ) * g.getCoefficient( j );
			rg_INT bound = j - i - 1;
			
			for(rg_INT k = 0;k <= bound;k++)
			{
				numericValuedMatrix[j - 1 - k][k + i] += numericValue;
			}
		}
	}

	return numericValuedMatrix;
}

rg_Matrix* rg_IntersectFunc::getCoefficientMatrix(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t)
{
	rg_Matrix* coefficientMatrix = new rg_Matrix[ 3 ];

	coefficientMatrix[ 0 ] = rg_IntersectFunc::getNumericValueFromOperationBtnTwoDifferentPolynomial(w_t, y_t);
	coefficientMatrix[ 1 ] = rg_IntersectFunc::getNumericValueFromOperationBtnTwoDifferentPolynomial(x_t, w_t);
	coefficientMatrix[ 2 ] = rg_IntersectFunc::getNumericValueFromOperationBtnTwoDifferentPolynomial(x_t, y_t);

	return coefficientMatrix;
}

rg_Matrix* rg_IntersectFunc::getCoefficientMatrix(rg_Polynomial x_t, rg_Polynomial y_t)
{
	rg_Matrix* coefficientMatrix = new rg_Matrix[ 3 ];

	return coefficientMatrix;
}

/////////////////////////////////////////////////////////////////////////////////////////



rg_QuarticPolynomial rg_IntersectFunc::substituteParametricIntoImplicit(const rg_RQBzCurve2D &curve, const rg_ImplicitEquation &implicit)
{
	rg_REAL x0 = curve.getCtrlPt(0).getX();
	rg_REAL x1 = curve.getCtrlPt(1).getX();
	rg_REAL x2 = curve.getCtrlPt(2).getX();
	rg_REAL y0 = curve.getCtrlPt(0).getY();
	rg_REAL y1 = curve.getCtrlPt(1).getY();
	rg_REAL y2 = curve.getCtrlPt(2).getY();
	rg_REAL w0 = curve.getWeight(0);
	rg_REAL w1 = curve.getWeight(1);
	rg_REAL w2 = curve.getWeight(2);

	rg_REAL a2 =    w0*x0 - 2*w1*x1 + w2*x2;
	rg_REAL a1 = -2*w0*x0 + 2*w1*x1;
	rg_REAL a0 =    w0*x0;

	rg_REAL b2 =    w0*y0 - 2*w1*y1 + w2*y2;
	rg_REAL b1 = -2*w0*y0 + 2*w1*y1;
	rg_REAL b0 =    w0*y0;

	rg_REAL d2 =    w0 - 2*w1 + w2;
	rg_REAL d1 = -2*w0 + 2*w1;
	rg_REAL d0 =    w0;


	rg_REAL *k = implicit.getAllCoeff();
  
    rg_REAL *coeff = new rg_REAL[5];
    coeff[0] =    d0*d0*k[0] +   a0*d0*k[1] +   b0*d0*k[2] +   a0*a0*k[3] +   a0*b0*k[4] +   b0*b0*k[5];

    coeff[1] =  2*d0*d1*k[0] +   a1*d0*k[1] +   a0*d1*k[1] +   b1*d0*k[2] +   b0*d1*k[2] + 2*a0*a1*k[3]
               +  a1*b0*k[4] +   a0*b1*k[4] + 2*b0*b1*k[5];

    coeff[2] =    d1*d1*k[0] + 2*d0*d2*k[0] +   a2*d0*k[1] +   a1*d1*k[1] +   a0*d2*k[1]
               +  b2*d0*k[2] +   b1*d1*k[2] +   b0*d2*k[2] +   a1*a1*k[3] + 2*a0*a2*k[3]
               +  a2*b0*k[4] +   a1*b1*k[4] +   a0*b2*k[4] +   b1*b1*k[5] + 2*b0*b2*k[5];

    coeff[3] =  2*d1*d2*k[0] +   a2*d1*k[1] +   a1*d2*k[1] +   b2*d1*k[2] +   b1*d2*k[2]
               +2*a1*a2*k[3] +   a2*b1*k[4] +   a1*b2*k[4] + 2*b1*b2*k[5];
    
    coeff[4] =    d2*d2*k[0] +   a2*d2*k[1] +   b2*d2*k[2] +   a2*a2*k[3] +   a2*b2*k[4] +   b2*b2*k[5];

    return rg_QuarticPolynomial(coeff);
}

rg_Point2D rg_IntersectFunc::intersectLineVsLine(const rg_Line<rg_Point2D> &pLine, const rg_Line<rg_Point2D> &nLine)
{
	// pLine equation : ax + by + c = 0
	// nLine equation : dx + ey + f = 0
		
	rg_REAL x;
	rg_REAL y;
	
	rg_REAL a;
	rg_REAL b;
	rg_REAL c;
	
	rg_REAL d;
	rg_REAL e;
	rg_REAL f;        

	// define the line equation of one line 
	if( rg_ZERO(pLine.getEP().getX() - pLine.getSP().getX(), resNeg10) ) 
	{
		a = 1;
		b = 0;
		c = -pLine.getSP().getX();
	}	
	else if( rg_ZERO(pLine.getEP().getY() - pLine.getSP().getY(), resNeg10) )
	{
		a = 0;
		b = -1;
		c = pLine.getSP().getY();
	}	
	else
	{
		a = (pLine.getEP().getY()-pLine.getSP().getY())
		   /(pLine.getEP().getX()-pLine.getSP().getX());
		b = -1;                  
		c = -a*pLine.getSP().getX() + pLine.getSP().getY();
	}

	// define the line equation of the other
	if( rg_ZERO(nLine.getEP().getX() - nLine.getSP().getX(), resNeg10) )
	{
		d = 1;
		e = 0;            
		f = -nLine.getSP().getX();
	}	   
	else if ( rg_ZERO(nLine.getEP().getY() - nLine.getSP().getY(), resNeg10) ) 
	{
		d = 0;
		e = -1;
		f = nLine.getSP().getY();
	}	
	else 
	{
		d = (nLine.getEP().getY()-nLine.getSP().getY())
		   /(nLine.getEP().getX()-nLine.getSP().getX());
		e = -1;
		f = -d*nLine.getSP().getX() + nLine.getSP().getY();
	}

    
    rg_Point2D intersection;
	rg_REAL D = a*e - b*d;
    if ( D != 0.0 ) {
        x = (b*f - c*e)/D;
        y = (c*d - a*f)/D;
        intersection.setPoint(x, y);
    }

    return intersection;
    

    //if ( rg_NE(a*e, b*d) ) 
    //{	                
    //    rg_REAL D = a*e - b*d;
    //    x = (b*f - c*e)/D;
    //    y = (c*d - a*f)/D;
    //    return rg_Point2D(x, y);
    //}       

    //return rg_Point2D(0.0, 0.0);
}

rg_Point2D *rg_IntersectFunc::intersectQBezierVsLine( rg_QBzCurve2D &curve, const rg_Line<rg_Point2D> &line)
{
	rg_REAL k    = (line.getEP().getY()-line.getSP().getY())
		         /(line.getEP().getX()-line.getSP().getX());
	rg_Point2D *pt = curve.getCtrlPts();

	rg_Point2D a(   pt[0] - 2*pt[1] + pt[2] );
	rg_Point2D b(2*(pt[1] -   pt[0])        );
	rg_Point2D c(   pt[0]                   );

	rg_REAL  p = k*a.getX() - a.getY();
	rg_REAL  q = k*b.getX() - b.getY();
	rg_REAL  r = k*c.getX() - c.getY();

	rg_REAL *t = rg_MathFunc::solveQuadraticEq(p, q, r);
				                
	rg_Point2D *intersection = new rg_Point2D[2];
	for(rg_INT i=0; i<2; i++)
		intersection[i] = a*t[i]*t[i] + b*t[i] + c;
	
	return intersection;
}

rg_Point3D* rg_IntersectFunc::intersectUnboundedLine3Ds(const rg_Line3D& line1, rg_Line3D& line2)
{
    rg_Point3D sp[2];
    sp[0]= line1.getSP();
    sp[1]= line2.getSP();

    rg_Point3D ep[2];
    ep[0]= line1.getEP();
    ep[1]= line2.getEP();

    rg_REAL* param = solveSimultaneousEq(
                      (ep[0] - sp[0]),
                      (sp[1] - ep[1]),
                      (sp[1] - sp[0]) );

    rg_Point3D* intersectPt = rg_NULL;
    if ( param != rg_NULL )
    {
        intersectPt = new rg_Point3D(( (1-param[0]) * sp[0] ) + (param[0] * ep[0]) );
        delete [] param;

        return intersectPt;
    }
    else
    {
        if ( line1.isColinearWith(line2) == rg_TRUE )
        {
            intersectPt = new rg_Point3D(sp[0]);
        }

        return intersectPt;
    }
}




void rg_IntersectFunc::intersectRQBzCurveVsRQBzCurve(const rg_RQBzCurve2D &curve1, const rg_RQBzCurve2D &curve2, rg_dList <rg_Point2D> &intersectList)
{
	//rg_QuarticPolynomial eq(rg_IntersectFunc::substituteParametricIntoImplicit(curve1, rg_IntersectFunc::implicitize(curve2)));
	rg_QuarticPolynomial eq(rg_IntersectFunc::substituteParametricIntoImplicit(curve1, rg_IntersectFunc::implicitizeNotUsingMathematica2(curve2))); //test
	
	// if roots are not rg_NULL, insert in the linked list
	rg_ComplexNumber *intersect_param = eq.solve();
	if(intersect_param)
	{
		for(rg_INT k = 0; k < 4; k++)
		{
			if(intersect_param[k].isPureRealNumber() && rg_BTORexclusive(0.0, intersect_param[k].getRealNumber(), 1.0))
			{
				rg_Point2D intersect_point = curve1.evaluatePt(intersect_param[k].getRealNumber());
				rg_REAL param = curve2.getParameter4Point(intersect_point);
				if( rg_BTORexclusive(0.0, param, 1.0) )
					intersectList.add(intersect_point);
			}
		}
	}
}


rg_INT rg_IntersectFunc::intersectRQBzCurveVsLine(const rg_RQBzCurve2D & rqbzCurve, const rg_Line2D & line, list<rg_Point2D>& intersetionPoints)
{
    rg_Polynomial* numerator = rg_NULL;
    numerator = rqbzCurve.convertNumerator2Polynomial();
    rg_Polynomial denominator = rqbzCurve.convertDenominator2Polynomial();

    rg_REAL coefficientOfX, coefficientOfY, coefficientOfConstant;
    line.get_coefficients_of_implicit_form_of_line_equation(coefficientOfX, coefficientOfY, coefficientOfConstant);

    rg_Polynomial equationForIntersections = coefficientOfX * numerator[0] + coefficientOfY * numerator[1] + coefficientOfConstant * denominator;
    rg_ComplexNumber* parametersOfIntersections = rg_NULL;
    parametersOfIntersections = equationForIntersections.solve();

    if (parametersOfIntersections)
    {
        for (int i = 0; i < 2; ++i)
        {
            if (parametersOfIntersections[i].isPureRealNumber() && rg_BTORexclusive(0.0, parametersOfIntersections[i].getRealNumber(), 1.0))
            {
                rg_Point2D intersectPt = rqbzCurve.evaluatePt(parametersOfIntersections[i].getRealNumber());
                rg_REAL parameterOfLineSeg = line.get_parameter_of_point_on_line_segment(intersectPt);
                if (rg_BTORexclusive(0.0, parameterOfLineSeg, 1.0))
                {
                    intersetionPoints.push_back(intersectPt);
                }
            }
        }
    }

    if (numerator != rg_NULL)
        delete[] numerator;

    if (parametersOfIntersections != rg_NULL)
        delete[] parametersOfIntersections;

    return intersetionPoints.size();
}



void rg_IntersectFunc::intersectRQBzCurveVsRQBzCurve(const rg_RQBzCurve2D &curve1, const rg_REAL &s0, const rg_REAL &s1,
								   const rg_RQBzCurve2D &curve2, const rg_REAL &t0, const rg_REAL &t1, 
								   rg_dList <rg_Point2D> &intersectList, 
								   rg_dList <rg_REAL*> &cubicParam)
{
	rg_QuarticPolynomial eq(rg_IntersectFunc::substituteParametricIntoImplicit(curve1, rg_IntersectFunc::implicitize(curve2)));
	//rg_QuarticPolynomial eq(rg_IntersectFunc::substituteParametricIntoImplicit(curve1, rg_IntersectFunc::implicitizeNotUsingMathematica2(curve2))); // test
	//rg_QuarticPolynomial eq(rg_IntersectFunc::substituteParametricIntoImplicit(curve1, rg_IntersectFunc::implicitizeNotUsingMathematica1(curve2))); // test

	// if roots are not rg_NULL, insert in the linked list
	rg_ComplexNumber *intersect_param_s = eq.solve();
	if(intersect_param_s)
	{
		for(rg_INT k = 0; k < 4; k++)
		{
			if(intersect_param_s[k].isPureRealNumber() && rg_BTORexclusive(0.0, intersect_param_s[k].getRealNumber(), 1.0))
			{
				rg_Point2D intersect_point   = curve1.evaluatePt(intersect_param_s[k].getRealNumber());
				rg_REAL  intersect_param_t = curve2.getParameter4Point(intersect_point);
				if( rg_BTORexclusive(0.0, intersect_param_t, 1.0) )
				{
					rg_REAL rqParam_s = intersect_param_s[k].getRealNumber();
					rg_REAL rqParam_t = intersect_param_t;

					rg_REAL *tmParam = new rg_REAL[2];
					tmParam[0] = s0 + (s1 - s0) * rqParam_s;
					tmParam[1] = t0 + (t1 - t0) * rqParam_t;
					cubicParam.add(tmParam);

					intersectList.add(intersect_point);
				}
			}
		}
	}
}

/*
rg_dList <rg_Point2D>intersectCBzCurveVsCBzCurve(const rg_CBzCurve2D &curve_s, const rg_CBzCurve2D &curve_t)
{
	rg_dList<rg_RQBzCurve2D> rqBzCurveList_s;
	rg_dList<rg_RQBzCurve2D> rqBzCurveList_t;

	rg_dList<rg_REAL> param_s = curve_s.approximateRQBzCurves(0.0, 1.0, rqBzCurveList_s);
	rg_dList<rg_REAL> param_t = curve_t.approximateRQBzCurves(0.0, 1.0, rqBzCurveList_t);

	// for each rational quadratic curve, compute intersection
	rg_dList<rg_Point2D> intersectList;
	rqBzCurveList_s.reset();
	do
	{
		rqBzCurveList_t.reset();
		
		do
		{
			rg_IntersectFunc::intersectRQBzCurveVsRQBzCurve(rqBzCurveList_s.getEntity(), rqBzCurveList_t.getEntity(), intersectList);
			rqBzCurveList_t.setCurrentNext();
		} while(!rqBzCurveList_t.isTail());

		rqBzCurveList_s.setCurrentNext();
	} while(!rqBzCurveList_s.isTail());

	return intersectList;
}

rg_RQBzCurve2D approximateCBzCurve2RQBzCurve(rg_CBzCurve2D &curve)
{
	rg_RQBzCurve2D rqcurve;
	
	rqcurve.makeRQBezier(curve);
	return rqcurve;
}
*/


