//********************************************************************
//
//	  FILENAME    : rg_RBzCurve3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_BzCurve3D
//           which define a rational Bezier rg_Curve in 3-D and its property. 
//                          
//	  CLASS NAME  : rg_RBzCurve3D
//
//    BASE CLASS  : rg_BzCurve3D
//      
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 11 Jul. 1997    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#include "rg_RBzCurve3D.h"
#include "rg_CurveSurfaceFunc.h"
#include "rg_dList.h"
#include "rg_MathFunc.h"
#include "rg_IntersectFunc.h"
rg_RBzCurve3D::rg_RBzCurve3D()
{
    weight = rg_NULL;
}

rg_RBzCurve3D::rg_RBzCurve3D(const rg_DEGREE &dgr)
: rg_BzCurve3D(dgr)
{
    weight = new rg_REAL[rg_BzCurve3D::getDegree() + 1];
}

rg_RBzCurve3D::rg_RBzCurve3D(const rg_DEGREE &dgr, const rg_Point3D* ctrlPts)
: rg_BzCurve3D(dgr, ctrlPts)
{
    weight = new rg_REAL[dgr + 1];
    for (rg_INDEX i=0; i<=dgr; i++)
        weight[i] = 1.0;
}

rg_RBzCurve3D::rg_RBzCurve3D(const rg_DEGREE &dgr, const rg_Point3D* ctrlPts, const rg_REAL* wght)
: rg_BzCurve3D(dgr, ctrlPts)
{
    weight = new rg_REAL [dgr + 1];
    for (rg_INDEX i=0; i<=dgr; i++)
        weight[i] = wght[i];
}

rg_RBzCurve3D::rg_RBzCurve3D(const rg_RBzCurve3D &curve)
: rg_BzCurve3D( curve )
{
    rg_DEGREE dgr = rg_BzCurve3D::getDegree();

    weight = new rg_REAL [dgr + 1];
    for (rg_INDEX i=0; i<=dgr; i++)
        weight[i] = curve.weight[i];
}

rg_RBzCurve3D::~rg_RBzCurve3D()
{
    if (weight != rg_NULL)
        delete [] weight;
}

//Operations-------------------------------------------------------------------

rg_Point3D rg_RBzCurve3D::evaluatePt(const rg_PARAMETER &t)
{
    rg_DEGREE dgr = rg_BzCurve3D::getDegree();

    rg_Point3D   ptOnCurve(0.0, 0.0, 0.0);
    rg_REAL      sumOfWeight = 0.0;

    for (rg_INDEX i = 0; i <= dgr; i++) {
        double bernsteinValue = bernstein(i, dgr, t);
        ptOnCurve += weight[i] * getCtrlPt(i)*bernsteinValue;
        sumOfWeight += weight[i] * bernsteinValue;
    }

    return rg_Point3D( ptOnCurve/sumOfWeight );
}

rg_Point3D rg_RBzCurve3D::evaluateDerivative(const rg_PARAMETER &t)
{
    rg_DEGREE dgr = rg_BzCurve3D::getDegree();

	rg_Point3D  p_t(0, 0, 0);
	rg_REAL   w_t = 0;
	rg_Point3D  p_prime_t;
	rg_REAL   w_prime_t = 0;
	
	rg_INT i = 0;
	for(i = 0; i <= dgr; i++) 
    {
		p_t += weight[i]*getCtrlPt(i)*bernstein(i, dgr, t); 
		w_t += weight[i]*bernstein(i, dgr, t);
	}

	for(i = 0; i <= dgr-1; i++) 
    {
		p_prime_t += dgr*(weight[i+1]*getCtrlPt(i+1)-weight[i]*getCtrlPt(i))
			               *bernstein(i, dgr-1, t);
		w_prime_t += dgr*(weight[i+1]-weight[i])*bernstein(i, dgr-1, t);
	}

	return rg_Point3D((p_prime_t*w_t - p_t*w_prime_t)/(w_t*w_t));

}


//Access elements--------------------------------------------------------------

rg_REAL rg_RBzCurve3D::getWeight(const rg_INDEX &i) const
{
    return weight[i];    
}

rg_REAL* rg_RBzCurve3D::getWeight() const
{
    return weight;
}

void rg_RBzCurve3D::setWeight(const rg_INDEX &i, const rg_REAL &w)
{
    weight[i] = w;
}

void rg_RBzCurve3D::setWeight(const rg_REAL *w)
{
    if (weight != rg_NULL)
        delete [] weight;

    rg_DEGREE dgr = rg_BzCurve3D::getDegree();

    weight = new rg_REAL[dgr + 1];
    for (rg_INDEX i=0; i<=dgr; i++)
        weight[i] = w[i];
}
void    rg_RBzCurve3D::setDegree(const rg_DEGREE& tDegree)
{
    rg_BzCurve3D::setDegree(tDegree);
    weight = new rg_REAL[tDegree + 1];
    for (rg_INT i=0; i< tDegree + 1; i++)
    {
        weight[i] = 1.0;
    }

}
void rg_RBzCurve3D::setCurve(const rg_RBzCurve3D &curve)
{
    rg_BzCurve3D::setCtrlPts(curve. getDegree()+1, curve.getCtrlPts());
    setWeight(curve.weight);
}

//  Intersection---------------------------------------------------------------

//  numOfIntersect : number of intersection point of Rational Quadratic 
//                   Bezier curve and plane.
//                   It is determined in this function. 
rg_Point3D* rg_RBzCurve3D::intersectRationalQuadraticBezierAndPlane(
                       const rg_Plane3D &plane, rg_INT &numOfIntersect) 
{
    if ( getDegree() != rg_QUADRATIC )
        return rg_NULL;
    
    //  Intersection of Rational Quadratic Bezier rg_Curve and rg_Point3D
    //
    //    expand rational quadratic Bezier curve about parameter, t
    //
    //              (w2P2 - 2w1P1 + w0P0)t^2 + (2w1P1 - 2w0P0)t + w0P0 
    //      B(t) = ----------------------------------------------------
    //                   (w2 - 2w1 + w0)t^2 + (2w1 - 2w0)t + w0
    //
    //      plane : Pl(u, v) = A + uB + vC
    //
    //      B(t) = Pl(u, v)
    //
    //      inner product (B * C)
    //
    //      C1t^2 + C2t + C3 = A'(D1t^2 + D2t + D3)
    //                  C1 = (w2P2 - 2w1P1 + w0P0)%(B * C)
    //                  C2 = (2w1P1 - 2w0P0)%(B * C)
    //                  C3 = w0P0%(B * C)
    //                  A' = A%(B * C)
    //                  D1 = (w2 - 2w1 + w0)
    //                  D2 = (2w1 - 2w0)
    //                  D3 = w0
    //
    //      -> (C1 - A'D1)t^2 + (C2 - A'D2)t + (C3 - A'D3) = 0
    //
    //---------------------------------------------------------------------
    
    rg_Point3D ctrlPt[rg_QUADRATIC+1];
    for (rg_INDEX i=0; i<=rg_QUADRATIC; i++)
        ctrlPt[i] = getCtrlPt(i);

    rg_Point3D coeffOfImplicitForm[3]
            = { (weight[2]*ctrlPt[2] - 2*weight[1]*ctrlPt[1] + weight[0]*ctrlPt[0]),
                (2*weight[1]*ctrlPt[1] - 2*weight[0]*ctrlPt[0]),
                (weight[0]*ctrlPt[0]) };

    rg_Point3D crsProduct = plane.getUVector()*plane.getVVector();

    rg_REAL coeff[3] = { (coeffOfImplicitForm[0]%crsProduct),
                      (coeffOfImplicitForm[1]%crsProduct),
                      (coeffOfImplicitForm[2]%crsProduct) };

    rg_REAL rightPart = plane.getPosVector()%crsProduct;

    coeff[0] -= rightPart*(weight[2] - 2*weight[1] + weight[0]);
    coeff[1] -= rightPart*(2*weight[1] - 2*weight[0]);
    coeff[2] -= rightPart*weight[0];

    rg_REAL innerRoot = sqrt( coeff[1]*coeff[1] - 4*coeff[0]*coeff[2] );
    if ( rg_LT( innerRoot, 0.0) )
    {
        numOfIntersect = 0;
        return rg_NULL;
    }
    else if ( rg_ZERO( innerRoot ) )
    {
        numOfIntersect = 1;
        rg_Point3D* intersectPt = new rg_Point3D;
        *intersectPt = evaluatePt( (-coeff[1])/(2*coeff[0]) );
        return intersectPt;
    }
    else
    {
        numOfIntersect = 2;
        rg_Point3D* intersectPt = new rg_Point3D[2];
        intersectPt[0] = evaluatePt( (-coeff[1] + innerRoot)/(2*coeff[0]) );
        intersectPt[1] = evaluatePt( (-coeff[1] - innerRoot)/(2*coeff[0]) );
        return intersectPt;
    }
}

// temporary member function for algrithm for computing between planar curves
// to be deleted !!
rg_Point3D* rg_RBzCurve3D::inflectionPointByLeeMethod(rg_ComplexNumber*& solution, rg_INT& count)
{
	//rg_Polynomial* n = convertNumerator2Polynomial();
	//rg_Polynomial d = convertDenominator2Polynomial();

	rg_DEGREE p = getDegree();
	rg_Matrix P(p + 1, 4);
	rg_Matrix M_p = rg_CurveSurfaceFunc::bezierToPowerMatrix(p + 1);

	rg_Polynomial* n = new rg_Polynomial[ 3 ];

	rg_Polynomial d( p );

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		P[ i ][ 0 ] = weight[ i ] * getCtrlPt( i ).getX(); // wX : X coordinate of Homogeneous Coordinate
		P[ i ][ 1 ] = weight[ i ] * getCtrlPt( i ).getY(); // wY : Y coordinate of Homogeneous Coordinate
		P[ i ][ 2 ] = weight[ i ] * getCtrlPt( i ).getZ(); // wZ : Z coordinate of Homogeneous Coordinate
		P[ i ][ 3 ] = weight[ i ]                           ; // w  : W weight
	}

	rg_Matrix coefficient = M_p * P;
	
	rg_INT j = 0;
	for(j = 0;j < 3;j++)
	{
		n[ j ].setDegree( p );
		for(rg_INT k = 0;k < p + 1;k++)
		{
			n[ j ].setCoefficient(k, coefficient[ k ][ j ]);
			d.setCoefficient(k, coefficient[ k ][ 3 ]);
		}
	}

	d.reconstruct();

	rg_Polynomial* first_n = new rg_Polynomial[2];
	rg_Polynomial* second_n = new rg_Polynomial[2];

	rg_Polynomial first_d;
	rg_Polynomial second_d;

	for( i = 0; i < 2; i++ )
	//for(rg_INT i = 0; i < 2; i++ )
	{
		first_n[i] = n[i].makeDerivative();
		second_n[i] = first_n[i].makeDerivative();
		first_n[i].reconstruct();
		second_n[i].reconstruct();
	}

	first_d = d.makeDerivative();
	first_d.reconstruct();
	second_d = first_d.makeDerivative();
	second_d.reconstruct();

	rg_Polynomial target = d * ( first_n[0] * second_n[1] - first_n[1] * second_n[0] )
						+ first_d * ( second_n[0] * n[1] - second_n[1] * n[0] )
						+ second_d * ( n[0] * first_n[1] - n[1] * first_n[0] );
	target.reconstruct();
	//rg_ComplexNumber* solution = target.solve();
	rg_ComplexNumber* sol = target.solve();

	//rg_INT count = 0;
	rg_dList<rg_Point3D> solList;
	rg_INT dg = target.getDegree();
	
	solution = new rg_ComplexNumber[dg];
	count = 0;
	for( j = 0; j < dg; j++ )
	//for(rg_INT j = 0; j < dg; j++ )
	{
		if( rg_ZERO( sol[j].getImaginaryNumber(), 0.0000001 )
			&& rg_GE( sol[j].getRealNumber(), 0.0 )
			&& rg_LE( sol[j].getRealNumber(), 1.0 ) )
		{
			solList.add( evaluatePt( sol[j].getRealNumber() ) );
			solution[count] = sol[ j ];
			count++;
		}
	}


	delete[] n;
	delete[] first_n;
	delete[] second_n;

	return solList.getArray();
}

rg_Polynomial* rg_RBzCurve3D::convertNumerator2Polynomial() const
{
	rg_DEGREE degree = getDegree();
    rg_REAL  *x = new rg_REAL[degree+1];
    rg_REAL  *y = new rg_REAL[degree+1];	
	rg_REAL  *z = new rg_REAL[degree+1];

	rg_INT i = 0;
	for(i = 0; i < degree+1; i++)
	{
		x[i] = weight[i] * getCtrlPt(i).getX();
		y[i] = weight[i] * getCtrlPt(i).getY();
		z[i] = weight[i] * getCtrlPt(i).getZ();
	}

    rg_Polynomial *ct = new rg_Polynomial[3];
    ct[0] = x[0]*bernstein(degree, 0);     // x(t)
    ct[1] = y[0]*bernstein(degree, 0);     // y(t)
	ct[2] = z[0]*bernstein(degree, 0);     // z(t)

    for(i = 1; i < degree + 1; i++)
    {
        ct[0] = ct[0] + x[i]*bernstein(degree, i);
        ct[1] = ct[1] + y[i]*bernstein(degree, i);
		ct[2] = ct[2] + z[i]*bernstein(degree, i);
    }

    delete[] x;
    delete[] y;
	delete[] z;

    return ct;
}

rg_Polynomial  rg_RBzCurve3D::convertDenominator2Polynomial() const
{
	rg_Polynomial ct;
	rg_DEGREE degree = getDegree();

	ct = weight[0]*bernstein(degree, 0);
    for( rg_INT i = 1; i < degree + 1; i++)
    {
        ct = ct + weight[i]*bernstein(degree, i);
    }

	return ct;
}

rg_RBzCurve3D rg_RBzCurve3D::getDerivative() const
{
    const rg_INT oldDegree=getDegree();
    const rg_INT newDegree=2*oldDegree;
    rg_RBzCurve3D output(newDegree);

    for (rg_INT i=0; i <= newDegree; i++ )
    {
        rg_Point3D ctrlPts( 0.0, 0.0,0.0);
        rg_REAL  weight=0.0;

        for( rg_INT j= rg_Max( 0, i-oldDegree); j <= rg_Min(oldDegree,i); j++ )
        {
            rg_Point3D temp1(0.0,0.0,0.0);
            rg_Point3D temp2(0.0,0.0,0.0);
            
            if ( j != 0 )
            {
                temp1 =   getWeight(i-j)
                        * ( getWeight(j)* getCtrlPt(j)  - getWeight(j-1)* getCtrlPt(j-1) );
                temp1 =  temp1
                       -  getWeight(i-j)*getCtrlPt(i-j)
                        * (  getWeight(j) - getWeight(j-1)) ;
                temp1 = j* temp1;

            }

            if ( oldDegree-j != 0 )
            {
                temp2 = getWeight(i-j) * ( getWeight(j+1)* getCtrlPt(j+1)  - getWeight(j)* getCtrlPt(j) );
                temp2 =  temp2
                       -  getWeight(i-j)*getCtrlPt(i-j)
                       * (  getWeight(j+1) - getWeight(j)) ;
                temp2 = (oldDegree-j)* temp2;
            }

            rg_REAL coefficient =  rg_MathFunc::combination( oldDegree, i-j)
                              * rg_MathFunc::combination( oldDegree,   j)
                              / rg_MathFunc::combination( newDegree,   i);


            ctrlPts+=  (( temp1+ temp2 )*coefficient);
            weight+=  getWeight(i-j)
                    * getWeight(j)
                    * coefficient;

        }
        output.setWeight(i,weight);
        output.setCtrlPt(i,ctrlPts/weight);
    }

    return output;
}

rg_RBzCurve3D rg_RBzCurve3D::makeHodograph() const
{
    const rg_INT oldDegree=getDegree();
    const rg_INT newDegree=2*oldDegree;
    rg_RBzCurve3D output(newDegree);

    for (rg_INT i=0; i <= newDegree; i++ )
    {
        rg_Point3D ctrlPts( 0.0, 0.0,0.0);
        rg_REAL  weight=0.0;

        for( rg_INT j= rg_Max( 0, i-oldDegree); j <= rg_Min(oldDegree,i); j++ )
        {
            rg_Point3D temp1(0.0,0.0,0.0);
            rg_Point3D temp2(0.0,0.0,0.0);
            
            if ( j != 0 )
            {
                temp1 =   getWeight(i-j)
                        * ( getWeight(j)* getCtrlPt(j)  - getWeight(j-1)* getCtrlPt(j-1) );
                temp1 =  temp1
                       -  getWeight(i-j)*getCtrlPt(i-j)
                        * (  getWeight(j) - getWeight(j-1)) ;
                temp1 = j* temp1;

            }

            if ( oldDegree-j != 0 )
            {
                temp2 = getWeight(i-j) * ( getWeight(j+1)* getCtrlPt(j+1)  - getWeight(j)* getCtrlPt(j) );
                temp2 =  temp2
                       -  getWeight(i-j)*getCtrlPt(i-j)
                       * (  getWeight(j+1) - getWeight(j)) ;
                temp2 = (oldDegree-j)* temp2;
            }

            rg_REAL coefficient =  (rg_REAL)rg_MathFunc::combination( oldDegree, i-j)
                              * rg_MathFunc::combination( oldDegree,   j)
                              / rg_MathFunc::combination( newDegree,   i);


            ctrlPts+=  (( temp1+ temp2 )*coefficient);
            weight+=  getWeight(i-j)
                    * getWeight(j)
                    * coefficient;

        }
        output.setWeight(i,weight);
        output.setCtrlPt(i,ctrlPts/weight);
    }

    return output;
}

rg_BzCurve3D rg_RBzCurve3D::makeScaledHodograph() const
{
	rg_INT newDegree = 2 * getDegree() - 2;

    rg_Point3D* ptsOfHodo = new rg_Point3D[newDegree + 1];
	for( rg_INT k = 0; k < newDegree+1 ; k++ )
	{
		ptsOfHodo[k]=rg_Point3D( 0.0, 0.0,0.0);
		rg_Point3D dPoint;
		for( rg_INT i = (rg_INT)rg_Max( 0, k - getDegree() + 1 ); i <= k/2; i++ )
		{
			dPoint = getWeight(i) * getWeight(k-i+1) * ( getCtrlPt(k-i+1) - getCtrlPt(i) );
			rg_REAL coeff = ((k-2*i+1)*rg_MathFunc::combination(getDegree(), i)*rg_MathFunc::combination(getDegree(), k-i+1)
				/ rg_MathFunc::combination(2*getDegree()-2, k));
			ptsOfHodo[k] = ptsOfHodo[k] + coeff * dPoint;
		}
	}

	rg_BzCurve3D output( newDegree, ptsOfHodo );
	delete[] ptsOfHodo;

	return output;
}
void rg_RBzCurve3D::makeConic( const rg_Point3D &startPt,
			             	   const rg_Point3D &startTangent,
			  			       const rg_Point3D &endPt,
						       const rg_Point3D &endTangent,
						       const rg_Point3D &passingPt )
{
    rg_Point3D tCtrlPts[3];
    tCtrlPts[0]=startPt;
    tCtrlPts[2]=endPt;

    rg_REAL    tWeights[3]={1.0,1.0,1.0};
    tWeights[0]=1.0;
    tWeights[1]=1.0;
    tWeights[2]=1.0;
    rg_Line3D line1(startPt,startPt+startTangent);
    rg_Line3D line2(endPt,endPt-endTangent);
    rg_Line3D baseLine(startPt,endPt);

    // find the 2nd control point
    rg_Point3D* temp=rg_IntersectFunc::intersectUnboundedLine3Ds(line1,line2);
    tCtrlPts[1]=*temp;
    delete temp;

    // find the parameter of passing point
    rg_Line3D   line3(tCtrlPts[1],passingPt);
    temp=rg_IntersectFunc::intersectUnboundedLine3Ds(line3,baseLine);
    rg_Point3D basePt=*temp;
    delete temp;

    rg_REAL a=sqrt( (tCtrlPts[0]-basePt).magnitude()
                /(tCtrlPts[2]-basePt).magnitude()   );
    rg_REAL pOfPassingPt=a/(1.0+a);

    // find the 2nd weight
    tWeights[1] =    pow(1.0-pOfPassingPt,2.0)*( (passingPt-tCtrlPts[0])%(tCtrlPts[1]-passingPt) )
                   + pow(pOfPassingPt,2.0)* ( (passingPt-tCtrlPts[2])%(tCtrlPts[1]-passingPt) );

    tWeights[1]=  tWeights[1]
              / (  2.0 * pOfPassingPt
                  *(1.0-pOfPassingPt)
                  *pow( (tCtrlPts[1]-passingPt).magnitude(), 2.0 ) ) ;
    const rg_INT tDgree=2;
    (*this).setDegree(tDgree);
    
    (*this).setCtrlPts(tDgree+1, tCtrlPts);
    (*this).setWeight(tWeights);
}

rg_RBzCurve3D& rg_RBzCurve3D::operator=(const rg_RBzCurve3D& temp)
{
	if ( this == &temp )
	{
		return *this;
	}
	rg_RBzCurve3D::setCurve(temp);

	return *this;
}


    
rg_REAL rg_RBzCurve3D::computeLength(const rg_INT& numSamplingPoint)
{
    rg_Point3D* ptOnCurve = new rg_Point3D[numSamplingPoint+1];

    int i=0;
	for ( i=0; i<=numSamplingPoint; ++i )  { 
		ptOnCurve[i] = evaluatePt( (rg_REAL)i/numSamplingPoint );
	}

    double length = 0.0;
    for ( i=0; i<numSamplingPoint; ++i ) {
        length += ptOnCurve[i].distance( ptOnCurve[i+1] );
    }

    delete [] ptOnCurve;

    return length;
}
