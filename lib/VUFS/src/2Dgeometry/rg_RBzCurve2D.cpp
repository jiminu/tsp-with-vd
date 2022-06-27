/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_RBzCurve2D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_RBzCurve2D 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jun 1997   
//    HISTORY     :
//         1. BY Hyung-Joo Lee 14 Jan. 1998
//             make   :  rg_BzCurve2D rg_RBzCurve2D::scaledHodo()
//         2. By Hyung-Joo Lee 15 Jan. 1998
//             make   :  rg_Polynomial* rg_RBzCurve2D::numeratorDerivative()
//         3. By Hyung-Joo Lee 15 Jan. 1998
//             make   :  rg_Polynomial rg_RBzCurve2D::denominatorDerivative()
//         4. By Hyung-Joo Lee 16 Jan. 1998
//             make   :  rg_Point2D* rg_RBzCurve2D::inflectionPointByHodograph()
//             make   :  rg_INT rg_RBzCurve2D::countInflectionPointByHodograph()
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include "rg_RBzCurve2D.h"
// Due to scaledHodo() 
#include "rg_RelativeOp.h"
#include "rg_MathFunc.h"
#include "rg_Polynomial.h"
#include "rg_dList.h"

rg_RBzCurve2D::rg_RBzCurve2D() :rg_BzCurve2D()
{
	weight = rg_NULL;
}

rg_RBzCurve2D::rg_RBzCurve2D(const rg_DEGREE &n) :rg_BzCurve2D(n)
{
	weight = new rg_REAL[degree+1];
	for(rg_INT i = 0; i <= degree; i++)
		weight[i] = 1.0;
}

rg_RBzCurve2D::rg_RBzCurve2D(const rg_DEGREE &n, const rg_Point2D *ctrlpt) :rg_BzCurve2D(n, ctrlpt)
{
	weight = new rg_REAL[degree+1];
	for(rg_INT i = 0; i <= degree; i++)
		weight[i] = 1.0;
}

rg_RBzCurve2D::rg_RBzCurve2D(const rg_RBzCurve2D &curve) :rg_BzCurve2D(curve) 
{
	weight = new rg_REAL[degree+1];
	for(rg_INT i = 0; i <= degree; i++)
		weight[i] = curve.weight[i];
}

rg_RBzCurve2D::~rg_RBzCurve2D()
{
	if(weight != rg_NULL) delete []weight;
}

rg_Point2D rg_RBzCurve2D::evaluatePt(const rg_PARAMETER &t) const
{
/*
	rg_Point2D  p_t;
	rg_REAL     w_t = 0;

	for(rg_INT i = 0; i <= degree; i++) {
		p_t += weight[i]*ctrlPts[i]*bernstein(i, degree, t); 
		w_t += weight[i]*bernstein(i, degree, t);
	}

	return rg_Point2D( p_t/w_t );
*/

	const rg_INT size=rg_BzCurve2D::getDegree()+1;
    rg_REAL* consideringWeights= rg_RBzCurve2D::getWeightVector();
    rg_Point2D* consideringPts= rg_BzCurve2D::getCtrlPts();

    rg_INT i = 0;
    for( i=0; i < size; i++ )
    {
        consideringPts[i]=consideringWeights[i]*consideringPts[i];
    }

	const rg_REAL rearCoefficient=t;
	const rg_REAL preCoefficient=1.0-t;

    for( i=0 ; i < size; i++ )
    {
        for( rg_INT j=0; j < size-1-i; j++ )
        {
            consideringWeights[j]= preCoefficient*consideringWeights[j]
                                  + rearCoefficient* consideringWeights[j+1];
            consideringPts[j]= preCoefficient*consideringPts[j]
                              + rearCoefficient* consideringPts[j+1];
        }
    }
    rg_Point2D output=consideringPts[0]/consideringWeights[0];
    delete[] consideringWeights;
    delete[] consideringPts;
    return output;

}

rg_Point2D rg_RBzCurve2D::evaluateDerivative(const rg_PARAMETER &t) const
{
	rg_Point2D  p_t;
	rg_REAL     w_t = 0;
	rg_Point2D  p_prime_t;
	rg_REAL     w_prime_t = 0;

	rg_INT i =0;
	for(i = 0; i <= degree; i++) {
		p_t += weight[i]*ctrlPts[i]*bernstein(i, degree, t); 
		w_t += weight[i]*bernstein(i, degree, t);
	}

	for(i = 0; i <= degree-1; i++) {
		p_prime_t += degree*(weight[i+1]*ctrlPts[i+1]-weight[i]*ctrlPts[i])
			               *bernstein(i, degree-1, t);
		w_prime_t += degree*(weight[i+1]-weight[i])*bernstein(i, degree-1, t);;
	}

	return rg_Point2D((p_prime_t*w_t - p_t*w_prime_t)/(w_t*w_t));
}

rg_REAL*   rg_RBzCurve2D::getWeightVector() const
{
	rg_REAL* output=new rg_REAL[degree+1];

	for( rg_INT i=0 ; i <= degree; i++ )
	{
		output[i]=weight[i];
	}
	return output;
}
void rg_RBzCurve2D::setDegree(const rg_DEGREE& tDegree)
{
    rg_BzCurve2D::setDegree(tDegree);
    if ( weight != rg_NULL )
    {
        delete[] weight;
    }

	weight = new rg_REAL[tDegree+1];
	for(rg_INT i=0; i<= tDegree; i++)
		weight[i] = 1.0;
}


void rg_RBzCurve2D::setWeight(const rg_REAL *w)
{
	if(weight != rg_NULL) delete[] weight;

	weight = new rg_REAL[degree+1];
	for(rg_INT i=0; i<= degree; i++)
		weight[i] = w[i];
}

// insert in History 1
//    due to setting only one weight.
void rg_RBzCurve2D::setWeight(const rg_INDEX& i,
                         const rg_REAL&  newWeight)
{
    if ( i >= 0 && i <= degree )
    {
        weight[i]=newWeight;
    }
}

void rg_RBzCurve2D::setCurve(const rg_RBzCurve2D &curve)
{
    if ( this == &curve )
    {
        return;
    }

    if ( ctrlPts != rg_NULL ) delete[] ctrlPts;
    if ( weight != rg_NULL ) delete[] weight;


    degree=curve.degree;
    ctrlPts=new rg_Point2D[degree+1];
    weight=new rg_REAL[degree+1];
    for ( rg_INT i=0; i < degree+1; i++ )
    {
	    ctrlPts[i]=curve.ctrlPts[i];
        weight[i]=curve.weight[i];
    }

}

rg_RBzCurve2D& rg_RBzCurve2D::operator=(const rg_RBzCurve2D& curve)
{
    if ( this != & curve )
    {
	    rg_RBzCurve2D::setCurve(curve);
    }

	return *this;
}
// inserted in histrory by Hyung-Joo Lee
// function to find scaled hodograph of rational bezier curve
// The code for scaledHodo() start here

rg_BzCurve2D rg_RBzCurve2D::scaledHodo() const
{
	rg_INT newDegree = 2 * getDegree() - 2;

    rg_Point2D* ptsOfHodo = new rg_Point2D[newDegree + 1];
	for( rg_INT k = 0; k < newDegree+1 ; k++ )
	{
		ptsOfHodo[k]=rg_Point2D( 0.0, 0.0);
		rg_Point2D dPoint;
		for( rg_INT i = (rg_INT)rg_Max( 0, k - getDegree() + 1 ); i <= k/2; i++ )
		{
			dPoint = getWeight(i) * getWeight(k-i+1) * ( getCtrlPt(k-i+1) - getCtrlPt(i) );
			rg_REAL coeff = ((k-2*i+1)*rg_MathFunc::combination(getDegree(), i)*rg_MathFunc::combination(getDegree(), k-i+1)
				/ rg_MathFunc::combination(2*getDegree()-2, k));
			ptsOfHodo[k] = ptsOfHodo[k] + coeff * dPoint;
		}
	}

	rg_BzCurve2D output( newDegree, ptsOfHodo );
	delete[] ptsOfHodo;

	return output;
}
// The code for scaledHodo() end here

// inserted in histrory by Hyung-Joo Lee
// function to find scaled hodograph of rational bezier curve
// The code for numeratorDerivative() start here

rg_Polynomial* rg_RBzCurve2D::numeratorDerivative() const
{
	rg_BzCurve2D derivative;
	derivative.setDegree( degree - 1 );
	rg_Point2D* temp = new rg_Point2D[derivative.getDegree() + 1];
	
	for( rg_INT i = 0; i <= derivative.getDegree(); i++ )
	{
		temp[i] = degree
			      * ( weight[i + 1] * ctrlPts[i + 1] - weight[i] * ctrlPts[i] );
	}

	derivative.setCtrlPts( temp );
	delete[] temp;

	rg_Polynomial* polyDerivative = new rg_Polynomial[2];

	polyDerivative = derivative.convertBzCurve2Polynomial();

	return polyDerivative;
}
// The code for numeratorDerivative() end here

// inserted in histrory by Hyung-Joo Lee
// function to find scaled hodograph of rational bezier curve
// The code for denumeratorDerivative() start here

rg_Polynomial rg_RBzCurve2D::denominatorDerivative() const
{
	rg_REAL* temp = new rg_REAL[degree];
	rg_INT i = 0;
	for( i = 0; i < degree; i++ )
	{
		temp[i] = degree * ( weight[i + 1] - weight[i] );
	}

    rg_Polynomial poly = temp[0] * bernstein(degree, 0);

    for( i = 1; i < degree + 1; i++)
    {
        poly = poly + temp[i] * bernstein(degree, i);
    }

	rg_Polynomial derivative = poly.makeDerivative();

	return derivative;
}
// The code for denumeratorDerivative() end here

// inserted in histrory by Hyung-Joo Lee
// function to find scaled hodograph of rational bezier curve
// The code for inflectionPointByHodograph() start here
rg_Point2D* rg_RBzCurve2D::inflectionPointByHodograph() const
{
	rg_BzCurve2D scaledHodo = this->scaledHodo();
	rg_BzCurve2D first_scaledHodo = scaledHodo.makeDerivative();

	rg_Polynomial* first = scaledHodo.convertBzCurve2Polynomial();
	rg_Polynomial* second = first_scaledHodo.convertBzCurve2Polynomial();

	rg_Polynomial target = first[0] * second[1] - first[1] * second[0];

	delete[] first;
	delete[] second;

	target.updateDegree();

	rg_ComplexNumber* solution = target.solve();

	rg_dList<rg_Point2D> solList;
	rg_INT dg = target.getDegree();
	for( rg_INT j = 0; j < dg; j++ )
	{
		if( rg_ZERO( solution[j].getImaginaryNumber(), 0.0000001 )
			&& rg_GE( solution[j].getRealNumber(), 0.0 )
			&& rg_LE( solution[j].getRealNumber(), 1.0 ) )
		{
			solList.add( evaluatePt( solution[j].getRealNumber() ) );
		}
	}

	delete[] solution;

	return solList.getArray();
}
// The code for inflectionPointByHodograph() end here

// inserted in histrory by Hyung-Joo Lee
// function to find scaled hodograph of rational bezier curve
// The code for countInflectionPointByHodograph() start here

rg_INT rg_RBzCurve2D::countInflectionPointByHodograph() const
{
	rg_BzCurve2D scaledHodo = this->scaledHodo();
	rg_BzCurve2D first_scaledHodo = scaledHodo.makeDerivative();

	rg_Polynomial* first = scaledHodo.convertBzCurve2Polynomial();
	rg_Polynomial* second = first_scaledHodo.convertBzCurve2Polynomial();

	rg_Polynomial target = first[0] * second[1] - first[1] * second[0];

	delete[] first;
	delete[] second;

	target.updateDegree();

	rg_ComplexNumber* solution = target.solve();

	rg_INT count = 0;
	for( rg_INT j = 0; j < target.getDegree(); j++ )
	{
		if( rg_ZERO( solution[j].getImaginaryNumber(), 0.0000001 )
			&& rg_GE( solution[j].getRealNumber(), 0.0 )
			&& rg_LE( solution[j].getRealNumber(), 1.0 ) )
			count++;
	}

	delete[] solution;

	return count;
}
// The code for countInflectionPointByHodograph() end here

// inserted in histrory by Hyung-Joo Lee
// function to find scaled hodograph of rational bezier curve
// The code for inflectionPointByLeeMethod() start here
rg_Point2D* rg_RBzCurve2D::inflectionPointByLeeMethod() const
{
	rg_Polynomial* n = convertNumerator2Polynomial();
	rg_Polynomial d = convertDenominator2Polynomial();

	d.updateDegree();

	rg_Polynomial* first_n = new rg_Polynomial[2];
	rg_Polynomial* second_n = new rg_Polynomial[2];

	rg_Polynomial first_d;
	rg_Polynomial second_d;

	for( rg_INT i = 0; i < 2; i++ )
	{
		first_n[i] = n[i].makeDerivative();
		second_n[i] = first_n[i].makeDerivative();
		first_n[i].updateDegree();
		second_n[i].updateDegree();
	}

	first_d = d.makeDerivative();
	first_d.updateDegree();
	second_d = first_d.makeDerivative();
	second_d.updateDegree();

	rg_Polynomial target = d * ( first_n[0] * second_n[1] - first_n[1] * second_n[0] )
						+ first_d * ( second_n[0] * n[1] - second_n[1] * n[0] )
						+ second_d * ( n[0] * first_n[1] - n[1] * first_n[0] );
	target.updateDegree();
	rg_ComplexNumber* solution = target.solve();

	rg_dList<rg_Point2D> solList;
	rg_INT dg = target.getDegree();
	for( rg_INT j = 0; j < dg; j++ )
	{
		if( rg_ZERO( solution[j].getImaginaryNumber(), 0.0000001 )
			&& rg_GE( solution[j].getRealNumber(), 0.0 )
			&& rg_LE( solution[j].getRealNumber(), 1.0 ) )
		{
			solList.add( evaluatePt( solution[j].getRealNumber() ) );
		}
	}

	delete[] n;
	delete[] first_n;
	delete[] second_n;
	delete[] solution;

	return solList.getArray();
}
// The code for inflectionPointByLeeMethod() end here

// inserted in histrory by Hyung-Joo Lee
// function to find scaled hodograph of rational bezier curve
// The code for countInflectionPointByLeeMethod() start here
rg_INT rg_RBzCurve2D::countInflectionPointByLeeMethod() const
{
	rg_Polynomial* n = convertNumerator2Polynomial();
	rg_Polynomial d = convertDenominator2Polynomial();

	d.updateDegree();

	rg_Polynomial* first_n = new rg_Polynomial[2];
	rg_Polynomial* second_n = new rg_Polynomial[2];

	rg_Polynomial first_d;
	rg_Polynomial second_d;

	for( rg_INT i = 0; i < 2; i++ )
	{
		first_n[i] = n[i].makeDerivative();
		second_n[i] = first_n[i].makeDerivative();
		first_n[i].updateDegree();
		second_n[i].updateDegree();
	}

	first_d = d.makeDerivative();
	first_d.updateDegree();
	second_d = first_d.makeDerivative();
	second_d.updateDegree();

	rg_Polynomial target = d * ( first_n[0] * second_n[1] - first_n[1] * second_n[0] )
						+ first_d * ( second_n[0] * n[1] - second_n[1] * n[0] )
						+ second_d * ( n[0] * first_n[1] - n[1] * first_n[0] );

	target.updateDegree();
	rg_ComplexNumber* solution = target.solve();

	rg_INT count = 0;
	rg_INT dg = target.getDegree();
	for( rg_INT j = 0; j < dg; j++ )
	{
		if( rg_ZERO( solution[j].getImaginaryNumber(), 0.0000001 )
			&& rg_GE( solution[j].getRealNumber(), 0.0 )
			&& rg_LE( solution[j].getRealNumber(), 1.0 ) )
		{
			count++;
		}
	}

	delete[] n;
	delete[] first_n;
	delete[] second_n;
	delete[] solution;

	return count;
}
// The code for countInflectionPointByLeeMethod() end here
rg_Point2D* rg_RBzCurve2D::inflectionPointByNormalMethod() const
{
	rg_RBzCurve2D firstDerivative = makeDerivative();
	rg_RBzCurve2D secondDerivative  = firstDerivative.makeDerivative();

	rg_Polynomial* first_n = firstDerivative.convertNumerator2Polynomial();
	rg_Polynomial* second_n = secondDerivative.convertNumerator2Polynomial();

	rg_Polynomial target = first_n[0] * second_n[1] - first_n[1] * second_n[0];
	target.updateDegree();

	delete[] first_n;
	delete[] second_n;

	rg_ComplexNumber* solution = target.solve();

	rg_dList<rg_Point2D> solList;
	rg_INT dg = target.getDegree();
	for( rg_INT j = 0; j < dg; j++ )
	{
		if( rg_ZERO( solution[j].getImaginaryNumber(), 0.0000001 )
			&& rg_GE( solution[j].getRealNumber(), 0.0 )
			&& rg_LE( solution[j].getRealNumber(), 1.0 ) )
		{
			solList.add( evaluatePt( solution[j].getRealNumber() ) );
		}
	}

	delete[] solution;

	return solList.getArray();
}

rg_INT rg_RBzCurve2D::countInflectionPointByNormalMethod() const
{
	rg_RBzCurve2D firstDerivative = makeDerivative();
	rg_RBzCurve2D secondDerivative  = firstDerivative.makeDerivative();

	rg_Polynomial* first_n = firstDerivative.convertNumerator2Polynomial();
	rg_Polynomial* second_n = secondDerivative.convertNumerator2Polynomial();

	rg_Polynomial target = first_n[0] * second_n[1] - first_n[1] * second_n[0];
	target.updateDegree();

	rg_ComplexNumber* solution = target.solve();

	delete[] first_n;
	delete[] second_n;

	rg_INT count = 0;
	rg_INT dg = target.getDegree();
	for( rg_INT j = 0; j < dg; j++ )
	{
		if( rg_ZERO( solution[j].getImaginaryNumber(), 0.0000001 )
			&& rg_GE( solution[j].getRealNumber(), 0.0 )
			&& rg_LE( solution[j].getRealNumber(), 1.0 ) )
		{
			count++;
		}
	}

	delete[] solution;

	return count;
}

rg_REAL* rg_RBzCurve2D::inflectionPointByDM( rg_INT& num ) const
{
	rg_Polynomial* n = convertNumerator2Polynomial();
	rg_Polynomial d = convertDenominator2Polynomial();

	d.updateDegree();

	rg_Polynomial* first_n = new rg_Polynomial[2];
	rg_Polynomial* second_n = new rg_Polynomial[2];

	rg_Polynomial first_d;
	rg_Polynomial second_d;

	for( rg_INT i = 0; i < 2; i++ )
	{
		first_n[i] = n[i].makeDerivative();
		second_n[i] = first_n[i].makeDerivative();
		first_n[i].updateDegree();
		second_n[i].updateDegree();
	}

	first_d = d.makeDerivative();
	first_d.updateDegree();
	second_d = first_d.makeDerivative();
	second_d.updateDegree();

	rg_Polynomial target = d * ( first_n[0] * second_n[1] - first_n[1] * second_n[0] )
						+ first_d * ( second_n[0] * n[1] - second_n[1] * n[0] )
						+ second_d * ( n[0] * first_n[1] - n[1] * first_n[0] );

	target.updateDegree();
	rg_ComplexNumber* solution = target.solve();

	rg_INT dg = target.getDegree();
	rg_dList<rg_REAL> solList;
	for( rg_INT j = 0; j < dg; j++ )
	{
		if( rg_ZERO( solution[j].getImaginaryNumber(), 0.0000001 )
			&& rg_GE( solution[j].getRealNumber(), 0.0 )
			&& rg_LE( solution[j].getRealNumber(), 1.0 ) )
		{
			solList.add( solution[j].getRealNumber() );
			num++;
		}
	}

	delete[] solution;

	return solList.getArray();
}


// inserted in histrory by Hyung-Joo Lee
// function to find scaled hodograph of rational bezier curve
// The code for convertNumerator2Polynomial() start here
rg_Polynomial* rg_RBzCurve2D::convertNumerator2Polynomial() const
{
    rg_REAL  *x = new rg_REAL[degree+1];
    rg_REAL  *y = new rg_REAL[degree+1];

    rg_INT i = 0;
	for(i = 0; i < degree+1; i++)
	{
		x[i] = weight[i] * ctrlPts[i].getX();
		y[i] = weight[i] * ctrlPts[i].getY();
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
// The code for convertNumerator2Polynomial() end here

// inserted in histrory by Hyung-Joo Lee
// function to find scaled hodograph of rational bezier curve
// The code for convertDenumeratorPolynomial() start here
rg_Polynomial rg_RBzCurve2D::convertDenominator2Polynomial() const
{
	rg_Polynomial ct;

	ct = weight[0]*bernstein(degree, 0);
    for( rg_INT i = 1; i < degree + 1; i++)
    {
        ct = ct + weight[i]*bernstein(degree, i);
    }

	return ct;
}
// The code for convertDenumeratorPolynomial() end here

rg_RBzCurve2D rg_RBzCurve2D::makeDerivative() const
{
    const rg_INT oldDegree=rg_BzCurve2D::getDegree();
    const rg_INT newDegree=2*oldDegree;
    rg_Point2D* ctrlPts=rg_BzCurve2D::getCtrlPts();
    rg_RBzCurve2D output(newDegree);

    for (rg_INT i=0; i <= newDegree; i++ )
    {
        rg_Point2D ctrlPt( 0.0, 0.0);
        rg_REAL  Weight=0.0;

        for( rg_INT j= rg_Max( 0, i-oldDegree); j <= rg_Min(oldDegree,i); j++ )
        {
            rg_Point2D temp1(0.0,0.0);
            rg_Point2D temp2(0.0,0.0);
            
            if ( j != 0 )
            {
                temp1 =   weight[i-j]
                        * ( weight[j]* weight[j]  - weight[j-1]* weight[j-1] );
                temp1 =  temp1
                       -  weight[i-j]*ctrlPts[i-j]
                        * (  weight[j] - weight[j-1] ) ;
                temp1 = j* temp1;

            }

            if ( oldDegree-j != 0 )
            {
                temp2 = weight[i-j] * ( weight[j+1]* ctrlPts[j+1]  - weight[j]* ctrlPts[j] );
                temp2 =  temp2
                       -  weight[i-j]*ctrlPts[i-j]
                       * (  weight[j+1] - weight[j]) ;
                temp2 = (oldDegree-j)* temp2;
            }

            rg_REAL coefficient =  (rg_REAL)rg_MathFunc::combination( oldDegree, i-j)
                              * (rg_REAL)rg_MathFunc::combination( oldDegree,   j)
                              / (rg_REAL)rg_MathFunc::combination( newDegree,   i);


            ctrlPt+=  (( temp1+ temp2 )*coefficient);
            Weight+=  weight[i-j]
                    * weight[j]
                    * coefficient;

        }
        output.setWeight(i,Weight);
        output.setCtrlPt(i,ctrlPt/Weight);
    }
    delete[] ctrlPts;
    return output;
/*
    rg_DEGREE degreeOfHodo=2*degree;
    rg_REAL*  weightsOfHodo=new rg_REAL[degreeOfHodo+1];
    
    // find out weights of derivative
    for( rg_INT i=0; i < degreeOfHodo+1; i++ )
    {
        weightsOfHodo[i]=0.0;
        rg_INT max= ( 0 >= i-degree ) ? 0 : i- degree;
        for( rg_INT j= max; j <= (i-1)/2 && i != 0 ; j++ )
        {
            weightsOfHodo[i]+= 2*weight[j]*weight[i-j]
                              *rg_MathFunc::combination(degree,j)*rg_MathFunc::combination(degree,i-j)
                              /rg_MathFunc::combination(degreeOfHodo,i);
        }

        if ( i%2 == 0 )
        {
            weightsOfHodo[i]+=  weight[i/2]*pow( rg_MathFunc::combination(degree,i/2), 2 )
                              / rg_MathFunc::combination( degreeOfHodo, i);
        }
    }

    // find out ctrl points of derivative
    
    rg_BzCurve2D raisedScaledHodo( scaledHodo() );
    raisedScaledHodo.raiseDegree(2);
    rg_Point2D* ctrlPts=new rg_Point2D[degreeOfHodo+1];

    for( i=0 ; i < degreeOfHodo+1; i++ )
    {
        ctrlPts[i]=raisedScaledHodo.getCtrlPoint(i)/weightsOfHodo[i];
    }

    // make derivative
    rg_RBzCurve2D output(degreeOfHodo,ctrlPts);
    output.setWeight(weightsOfHodo);
    delete []ctrlPts;
    delete []weightsOfHodo;
    return output;
*/
}


