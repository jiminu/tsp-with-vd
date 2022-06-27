//********************************************************************
//
//	  FILENAME    : CurveSurfaceFunc.cpp
//	  
//    DESCRIPTION : 
//           This is the implementatione of external functions for curve and surface
//                          
//
//    AUTHOR      : Deok-Soo Kim, Dong-Gyou Lee
//
//    START DATE  : 17 Mar. 1998    
//
//    HISTORY     :
//				
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#include "rg_CurveSurfaceFunc.h"
#include "rg_MathFunc.h"

#include "rg_Point3D.h"
#include "rg_Matrix.h"
#include "rg_IntersectFunc.h"

//*****************************************************************************
//
//    FUNCTION    : rg_CurveSurfaceFunc::bezierToPowerMatrix
//    DESCRIPTION : 
//                  This function is for conversion between b-spline and 
//                  piecewise power basis form.
//                  Return value is Mp matrix!!!!
//                  
//                  form : Cj(s) = [Bi,p(s)][Pi] = [si]Mp[Pi] = [si][ai]
//
//                         Mp(j,k) = | 0                               j=0,...,k-1
//                                   | (-1)^(j-k) * pCk * (p-k)C(j-k)  j=k,...,p
//
//    AUTHOR      : Dong-Gyou Lee
//    START DATE  : 17 Mar. 1998   
//    REFERENCE   : The NURBS Book, Les Piegl & Wayne Tiller, p.265 
//
//*****************************************************************************
rg_Matrix rg_CurveSurfaceFunc::bezierToPowerMatrix( const rg_INT& order )
{
	rg_Matrix Mp( order, order );

	for( rg_INT k=0; k < order; k++ )
	{
		for( rg_INT j=k; j < order; j++ )
		{
			if( ( j-k )%2 )
				Mp[j][k] = (rg_INT) -( rg_MathFunc::combination( order-1, k ) * rg_MathFunc::combination( (order-1)-k, ( j-k ) ) );
			else
				Mp[j][k] = (rg_INT) ( rg_MathFunc::combination( order-1, k ) * rg_MathFunc::combination( (order-1)-k, ( j-k ) ) );
		}
	}
	return Mp;
}

//*****************************************************************************
//
//    FUNCTION    : rg_CurveSurfaceFunc::powerToBezierMatrix
//    DESCRIPTION : 
//                  This function is for conversion between b-spline and 
//                  piecewise power basis form.
//                  Return value is Mi(inverse of Mp) matrix!!!!
//                  
//                  form : Cj(s) = [Bi,p(s)][Pi] = [si]Mp[Pi] 
//                               = [ui]RpMp[Pi]  = [ui][ai]    = Cj(u)
//
//                         So, [Pi] = MiRi[ai] 
//
//                         Mi(j,k) = | 0                               j=0,...,k-1
//                                   | 1/Mp(k,k)                       j=k
//                                   | -{Sum(Mp(j,i)Mi(i,k))}/Mp(j,j)  j=k+1,...,p
//
//    AUTHOR      : Dong-Gyou Lee
//    START DATE  : 17 Mar. 1998   
//    REFERENCE   : The NURBS Book, Les Piegl & Wayne Tiller, p.271 
//
//*****************************************************************************
rg_Matrix rg_CurveSurfaceFunc::powerToBezierMatrix( const rg_INT& order )
{
	rg_Matrix Mp( rg_CurveSurfaceFunc::bezierToPowerMatrix( order ) );
	rg_Matrix Mi( order, order );

	for( rg_INT k=0; k < order; k++ )
	{
		Mi[k][k] = 1./Mp[k][k];
		
		for( rg_INT j=(k+1); j < order; j++ )
		{
			rg_REAL sum = 0.;
			for( rg_INT i=k; i < j; i++ )
				sum = sum + Mp[j][i] * Mi[i][k];

			Mi[j][k] = -sum/Mp[j][j];
		}
	}

	return Mi;
}

rg_RPolynomialCurve2D rg_CurveSurfaceFunc::convertToRPolynomialCurve2D( const rg_RBzCurve2D& curve)
{
    const rg_INT degree=curve.getDegree();
    rg_REAL*     weight=curve.getWeightVector();
    rg_Point2D*   ctrlPts=curve.getCtrlPts();

    rg_REAL*    coeffOfDenominator= new rg_REAL[degree+1];
    rg_Point2D* coefficient=new rg_Point2D[degree+1];

    rg_INT i = 0;
	for( i=0 ; i < degree+1; i++ )
    {
        coefficient[i].setPoint(0.0,0.0);
        coeffOfDenominator[i]=0.0;
    }

    for ( i=0 ; i < degree+1; i++)
    {
        for( rg_INT j=0; j <= i; j++ )
        {
            if ( (i-j)%2 == 0 ) 
            {
                coefficient[i]+=ctrlPts[j]*weight[j]*( (rg_REAL)rg_MathFunc::combination(i,j) );
                coeffOfDenominator[i]+=weight[j]*( (rg_REAL)rg_MathFunc::combination(i,j) );
            }
            else
            {
                coefficient[i]-=ctrlPts[j]*weight[j]*( (rg_REAL)rg_MathFunc::combination(i,j) );
                coeffOfDenominator[i]-=weight[j]*( (rg_REAL)rg_MathFunc::combination(i,j) );

            }
        }
        coefficient[i]=coefficient[i]*(rg_REAL)rg_MathFunc::combination(degree,i);
        coeffOfDenominator[i]=coeffOfDenominator[i]*(rg_REAL)rg_MathFunc::combination(degree,i);
    }

    rg_RPolynomialCurve2D output(degree);

    output.setCoefficients(coefficient);
    output.setCoefficientsOfDenominator(coeffOfDenominator);

    delete[] weight;
    delete[] ctrlPts;
    delete[] coeffOfDenominator;
    delete[] coefficient;

    return output;
}

////  For Reparameterization

//*****************************************************************************
//
//    FUNCTION    : rg_CurveSurfaceFunc::reparameterMatrix
//    DESCRIPTION : 
//                  This function is for reparameterization between b-spline
//                  segment and bezier.
//                  The jth b-spline curve segment Cj(u) is defined on u 
//                  in [uj, uj+1] and Cj(s) is bezier form of that segment on
//                  s in [0, 1].
//                               1              uj
//                  so, s = ----------- u - ----------- = cu + d
//                           uj+1 - uj       uj+1 - uj
//
//                  and u = ( uj+1 - uj )s + uj         = cs + d
//
//                      R(i,j) = | 0                     j < i
//                               | jCi * c^i * d^(j-i)   j >= i  i,j = 0,...,p
//
//    AUTHOR      : Dong-Gyou Lee
//    START DATE  : 17 Mar. 1998   
//    REFERENCE   : The NURBS Book, Les Piegl & Wayne Tiller, p.270 
//
//*****************************************************************************
rg_Matrix rg_CurveSurfaceFunc::reparameterMatrix( const rg_INT& order, const rg_REAL& multiValue, const rg_REAL& addValue )
{
	rg_Matrix R( order, order );

	for( rg_INT i=0; i < order; i++ )
	{
		for( rg_INT j=i; j < order; j++ )
			R[i][j] = rg_MathFunc::combination( j, i ) * pow( multiValue, i ) * pow( addValue, (j-i) );
	}
	return R;
}




              





    

