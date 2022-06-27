/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_BzCurve2D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_BzCurve2D 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jun 1997    
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include <math.h>
#include "rg_BzCurve2D.h"
#include "rg_RelativeOp.h"
#include "rg_MathFunc.h"
#include "rg_dList.h"
#include "sortFunc.h"// for the usage of QuickSort() in blossomAt()

rg_BzCurve2D::rg_BzCurve2D()
{
    degree=-1;
    ctrlPts = rg_NULL;
}

rg_BzCurve2D::rg_BzCurve2D(const rg_DEGREE &n) :degree(n)
{
	ctrlPts = new rg_Point2D[degree+1];
}

rg_BzCurve2D::rg_BzCurve2D(const rg_DEGREE &n, const rg_Point2D *ctrlpt) :degree(n)
{
	ctrlPts = new rg_Point2D[degree+1];
	for(rg_INT i=0; i<=degree; i++)
		ctrlPts[i] = ctrlpt[i];
}

rg_BzCurve2D::rg_BzCurve2D(const rg_BzCurve2D &curve) :degree(curve.degree)
{
	ctrlPts = new rg_Point2D[degree+1];
	for(rg_INT i=0; i<=degree;i++)
		ctrlPts[i] = curve.ctrlPts[i];
}
/*
rg_BzCurve2D::rg_BzCurve2D(const rg_BzCurve3D &curve) :degree(curve.getDegree())
{
	ctrlPts = new rg_Point2D[degree+1];
	for(rg_INT i=0; i<=degree;i++)
		ctrlPts[i] = curve.getCtrlPt( i );
}
*/
rg_BzCurve2D::~rg_BzCurve2D()
{
	if(ctrlPts != rg_NULL) delete []ctrlPts;
}

rg_Point2D rg_BzCurve2D::evaluatePt(const rg_PARAMETER &t) const
{

	const rg_INT size=rg_BzCurve2D::getDegree()+1;
    rg_Point2D* consideringPts= getCtrlPts();

    const rg_REAL rearCoefficient=t;
	const rg_REAL preCoefficient=1.0-t;
    for( rg_INT i=0 ; i < size; i++ )
    {
        for( rg_INT j=0; j < size-1-i; j++ )
        {
            consideringPts[j]=   preCoefficient* consideringPts[j]
                              + rearCoefficient* consideringPts[j+1];
        }
    }
    const rg_Point2D output=consideringPts[0];
    delete[] consideringPts;

    return output;
/*

	rg_Point2D  pt;
	for(rg_INT i = 0; i <= degree; i++)
		pt += ctrlPts[i]*bernstein(i, degree, t); 

	return pt;

*/

}

rg_BzCurve2D rg_BzCurve2D::makeDerivative() const
{
	rg_Point2D *ctrlpt = new rg_Point2D[degree];
	for(rg_INT i=0; i<=degree-1; i++)
		ctrlpt[i] = degree*(ctrlPts[i+1]-ctrlPts[i]);
	
	return rg_BzCurve2D(degree-1, ctrlpt);
}

rg_Point2D** rg_BzCurve2D::deCasteljau(const rg_PARAMETER &t)
{
	rg_Point2D **b = new rg_Point2D*[degree+1];
	for(rg_INT j=0; j<=degree; j++)
		b[j] = new rg_Point2D[degree+1];
	
	// range of parameter value : [0, 1] 
	// if not satisfied, exit.
	for(rg_INT i=0; i<=degree; i++)
		b[0][i] = ctrlPts[i];

	for(rg_INT r=1; r<=degree; r++)
		for(rg_INT i=0; i<=degree-r; i++)
			b[r][i] = (1-t)*b[r-1][i] + t*b[r-1][i+1];

	return b;
}

// The Berstein polynomial of the following nonrecursive form is faster than that of the recursive form
//
// Nonrecursive form : Bi, n(t) = nCi t^i (1-t)^(n-i)
// Recursive form      : Bi, n(t) = (1-t) Bi, n-1(t) + t Bn-1, i-1(t)
// Farin's book pp.41 is recommended. 

rg_REAL rg_BzCurve2D::bernstein(const rg_DEGREE &i, const rg_DEGREE &n, const rg_PARAMETER &t)
{
	if( i < 0 || n-i < 0 )
    {
        return 0.0;	
    }
	else
    {
		return (rg_REAL)rg_MathFunc::combination(n, i)* pow(t, i) * pow(1-t, n-i);
	}
}

rg_Polynomial rg_BzCurve2D::bernstein(const rg_DEGREE &n, const rg_INDEX &i) const
{
	// define t
	rg_REAL    *t    = new rg_REAL[2];
	t[0] = 0.0;
	t[1] = 1.0;

	// define (1-t)
	rg_REAL    *t_1  = new rg_REAL[2];
	t_1[0] =  1.0;
    t_1[1] = -1.0;

	rg_Polynomial t_i(1, t);
	rg_Polynomial t_1_n_i(1, t_1);
    delete[] t;
    delete[] t_1;

	// define the Bernstein polynomial : (n, i) t^i (1-t)^(n-i)
	rg_Polynomial result = rg_MathFunc::combination(n, i)*t_i.power(i)*t_1_n_i.power(n-i);
    
    return result;
}

rg_Point2D rg_BzCurve2D::blossomAt(const rg_PARAMETER* parameterList, const rg_REAL& lowerBound, const rg_REAL& upperBound) const
{
	rg_Point2D** b = new rg_Point2D* [degree + 1];
	rg_INT i = 0;
	for(i = 0;i < degree + 1;i++)
		b[ i ] = new rg_Point2D[degree + 1 - i];

	for(i = 0;i < degree + 1;i++)
		b[ 0 ][ i ] = ctrlPts[ i ];
	
	for(i = 1;i < degree + 1;i++)
	{
		for(rg_INT j = 0;j < degree + 1 - i;j++)
		{
			b[ i ][ j ] = (upperBound - parameterList[i - 1])/(upperBound - lowerBound) * b[i - 1][ j ] + (parameterList[i - 1] - lowerBound)/(upperBound - lowerBound)* b[i - 1][j + 1];
		}
	}

	for(i = 0;i < degree;i++)
		delete[] b[ i ];

	return b[degree][ 0 ];
}

/*
rg_REAL rg_BzCurve2D::rg_MathFunc::combination(const rg_REAL& n, const rg_REAL& i) const
{
	if( rg_GE(i, 0.0) && rg_LE(i, n) )
	{
		return rg_MathFunc::factorial(n)/rg_MathFunc::factorial(i)/rg_MathFunc::factorial(n-i);
	}

	else 
	{
		return 0.0;
	}
}

rg_REAL rg_BzCurve2D::rg_MathFunc::factorial(const rg_REAL &n) const
{
	if( n ) return n * rg_MathFunc::factorial(n-1);
	else return 1;
}
*/

void rg_BzCurve2D::setOrder(const rg_DEGREE &n)
{

	degree = n-1;
    if ( ctrlPts != rg_NULL)
    {
        delete[] ctrlPts;
    }

	if ( degree < 0 ) 
	{
		ctrlPts=rg_NULL;
		return;
	}

    ctrlPts= new rg_Point2D[degree+1];
    for( rg_INT i=0; i < degree+1; i++ )
    {
        ctrlPts[i]=rg_Point2D(0.0,0.0);
    }

}

void rg_BzCurve2D::setDegree(const rg_DEGREE &n)
{
	degree = n;
    if ( ctrlPts != rg_NULL)
    {
        delete[] ctrlPts;
    }
	if ( degree < 0 ) 
	{
		ctrlPts=rg_NULL;
		return;
	}

    ctrlPts= new rg_Point2D[degree+1];
    for( rg_INT i=0; i < degree+1; i++ )
    {
        ctrlPts[i]=rg_Point2D(0.0,0.0);
    }

}


rg_Point2D* rg_BzCurve2D::getCtrlPts() const
{
    rg_Point2D* output=new rg_Point2D[degree+1];
    for( rg_INT i=0; i < degree+1; i++ )
    {
        output[i]=ctrlPts[i];
    }

	return output;
}

void rg_BzCurve2D::setCtrlPt(const rg_INDEX& i, const rg_Point2D& point)
{
	if( (ctrlPts != rg_NULL) 
		&& (degree >= i) )
	{
		ctrlPts[i] = point;
	}
}

void rg_BzCurve2D::setCtrlPts(const rg_Point2D *point)
{
	if(ctrlPts != rg_NULL) delete[] ctrlPts;

	ctrlPts = new rg_Point2D[degree+1];
	for(rg_INT i=0; i <= degree;i++)
		ctrlPts[i] = point[i];
}

void rg_BzCurve2D::setCurve(const rg_BzCurve2D &curve)
{
	if(ctrlPts != rg_NULL) delete[] ctrlPts;

    degree=curve.degree;
	ctrlPts = new rg_Point2D[degree+1];
	for(rg_INT i=0; i<=degree;i++)
		ctrlPts[i] = curve.ctrlPts[i];
}

rg_BzCurve2D& rg_BzCurve2D::operator=(const rg_BzCurve2D& curve)
{
	setCurve(curve);
	
	return *this;
}

rg_Polynomial* rg_BzCurve2D::convertBzCurve2Polynomial() const
{
    rg_REAL  *x = new rg_REAL[degree+1];
    rg_REAL  *y = new rg_REAL[degree+1];

    rg_INT i = 0;
	for(i = 0; i < degree+1; i++)
	{
		x[i] = ctrlPts[i].getX();
		y[i] = ctrlPts[i].getY();
	}

    rg_Polynomial *ct = new rg_Polynomial[2];	
    ct[0] = x[0]*bernstein(degree, 0);     // x(t)
    ct[1] = y[0]*bernstein(degree, 0);     // y(t)

    for(i = 1; i < degree + 1; i++)
    {
        ct[0] = ct[0] + x[i]*bernstein(degree, i);
        ct[1] = ct[1] + y[i]*bernstein(degree, i);
    }

    delete[] x;
    delete[] y;

    return ct;
}

void rg_BzCurve2D::raiseDegree(const rg_DEGREE& raisingTimes)
{
    rg_DEGREE   oldDegree=degree;
    degree=oldDegree+raisingTimes;
    rg_Point2D* oldCtrlPts=ctrlPts;
    ctrlPts= new rg_Point2D[degree+1];

    for( rg_INT i=0; i < degree+1; i++ )
    {
        ctrlPts[i]=rg_Point2D( 0.0, 0.0);

        rg_INT minJ=( i-raisingTimes < 0 ) ? 0:i-raisingTimes;
        rg_INT maxJ=( i > oldDegree ) ? oldDegree:i;
        for( rg_INT j= minJ ; j <= maxJ; j++ )
        {
            ctrlPts[i]+= oldCtrlPts[j]
                          *rg_MathFunc::combination(oldDegree,j)
                          *rg_MathFunc::combination(raisingTimes,i-j)
                          /rg_MathFunc::combination(degree,i);
        }
    }
    delete []oldCtrlPts;
}

rg_REAL* rg_BzCurve2D::inflectionPointByHodograph(rg_INT& numIPts) const
{
	rg_BzCurve2D first_derivative  = makeDerivative();
	rg_BzCurve2D second_derivative = first_derivative.makeDerivative();

	rg_Polynomial* first  = first_derivative.convertBzCurve2Polynomial();
	rg_Polynomial* second = second_derivative.convertBzCurve2Polynomial();

	rg_Polynomial target = first[0] * second[1] - first[1] * second[0];

	delete[] first;
	delete[] second;

	target.updateDegree();

	rg_ComplexNumber* solution = target.solve();

	rg_dList<rg_REAL> solList;
	for( rg_INT i = 0; i<target.getDegree(); i++ )
	{
		if( rg_ZERO( solution[i].getImaginaryNumber() )
			&& rg_GE( solution[i].getRealNumber(), 0.0 )
			&& rg_LE( solution[i].getRealNumber(), 1.0 ) )
		{
			solList.add( solution[i].getRealNumber() );
		}
	}

	delete[] solution;

	numIPts = solList.getSize();

	return solList.getArray();
}


