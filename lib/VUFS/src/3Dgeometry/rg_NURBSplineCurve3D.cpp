//********************************************************************
//
//	  FILENAME    : rg_NURBSplineCurve3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_NURBSplineCurve3D 
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun 1996    
//    History     : 
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************


#include "rg_NUBSplineCurve3D.h"
#include "rg_NURBSplineCurve3D.h"

#include "rg_BandedMatrix.h"
#include "rg_RelativeOp.h"
#include "rg_CurveSurfaceFunc.h"
#include "rg_TMatrix3D.h"
#include "rg_NURBSCurveIntersector.h"
#include "rg_MathFunc.h"
#include "rg_GeoFunc.h"
#include "rg_BoundingBox3D.h"
rg_NURBSplineCurve3D::rg_NURBSplineCurve3D()
: rg_NUBSplineCurve3D()
{
	weights = rg_NULL;
}

rg_NURBSplineCurve3D::rg_NURBSplineCurve3D( const rg_ORDER &newOrder )
: rg_NUBSplineCurve3D(newOrder)
{
	weights = rg_NULL;
}

//// Constructor : March 13 1997
rg_NURBSplineCurve3D::rg_NURBSplineCurve3D( const unsigned rg_INT &newID, 
                                      const rg_Planarity    &newPlanarity,
                                      const rg_ORDER        &newOrder, 
                                      const rg_INT          &num, 
                                      rg_Point3D*              newControlP,
                                      rg_REAL*               newKnotVector,
                                      rg_REAL*               newWeightVector )
: rg_NUBSplineCurve3D(newID,
                   newPlanarity,
                   newOrder,
                   num,
                   newControlP,
                   newKnotVector)
{
    rg_INT n = rg_BSplineCurve3D::getNumOfCtrlPts();

	if ( n > 0 )
	{
		weights = new rg_REAL [n];
	
		for (rg_INT i=0; i<n; i++)
			weights[i] = newWeightVector[i];
	}
	else
	{
		weights=rg_NULL;
	}

}

rg_NURBSplineCurve3D::rg_NURBSplineCurve3D( const rg_BSplineCurve3D &curve )
: rg_NUBSplineCurve3D( curve )
{
    rg_INT n = curve.getNumOfCtrlPts();

	if ( n > 0 )
	{
		weights = new rg_REAL [n];
	
		for (rg_INT i=0; i<n; i++)
			weights[i] = 1.0;
	}
	else
	{
		weights =rg_NULL;
	}
}

rg_NURBSplineCurve3D::rg_NURBSplineCurve3D( const rg_NUBSplineCurve3D &curve )
: rg_NUBSplineCurve3D( curve )
{
    rg_INT n = curve.getNumOfCtrlPts();

	if ( n > 0 )
	{
		weights = new rg_REAL [n];
	    for(rg_INT i=0; i<n; i++)
		    weights[i] = 1.0;
	}
	else
	{
		weights =rg_NULL;
	}

}
    
//// Copy Constructor : March 13 1997
rg_NURBSplineCurve3D::rg_NURBSplineCurve3D( const rg_NURBSplineCurve3D &curve )
: rg_NUBSplineCurve3D(curve.getID(),
                   curve.getPlanarity(),
                   curve.getOrder(),
                   curve.getNumOfCtrlPts(),
                   curve.getCtrlPts(),
                   curve.getKnotVector() )
{
    rg_INT n = curve.getNumOfCtrlPts();

	if ( n > 0 )
	{
		weights = new rg_REAL [n];
	    for(rg_INT i=0; i<n; i++)
		    weights[i] = curve.weights[i];
	}
	else
	{
		weights =rg_NULL;
	}

}
            
rg_NURBSplineCurve3D::~rg_NURBSplineCurve3D()
{
    if ( weights != rg_NULL )
        delete [] weights;
}

////    Get Functions. : March 13 1997

rg_REAL rg_NURBSplineCurve3D::getWeight(const rg_INT &i) const
{
    return weights[i];
}

rg_REAL* rg_NURBSplineCurve3D::getWeightVector() const
{
    return weights;
}

rg_RationalPolynomial** rg_NURBSplineCurve3D::makePolynomialFormCurveUsingCurveDecomp() const
{
	/*
	// In these codes, reparameterization of each Bezier curve segment is not included
	rg_INT NumOfNonZeroLengthKnotSpan = getNumOfNonZeroLengthKnotSpan();
	rg_DEGREE p = getOrder() - 1;

	rg_BzCurve3D* BezierCurveList = decomposeCurveIntoBezierSegment();

	rg_PolynomialWithBound** polynomialFormCurve = new rg_PolynomialWithBound* [ 3 ];

	rg_INT i = 0;
	for(i = 0;i < 3;i++)
		polynomialFormCurve[ i ] = new rg_PolynomialWithBound[NumOfNonZeroLengthKnotSpan];

	rg_PolynomialWithBound tempX( p ), tempY( p ), tempZ( p );
	rg_REAL* distinctKnotValue = getDistinctKnotValues();

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		rg_Polynomial* PolyFormBezierCurve = BezierCurveList[ i ].convertBzCurve2Polynomial();

		tempX = PolyFormBezierCurve[ 0 ];
		tempX.setBound(distinctKnotValue[ i ], distinctKnotValue[i + 1]);
		tempY = PolyFormBezierCurve[ 1 ];
		tempY.setBound(distinctKnotValue[ i ], distinctKnotValue[i + 1]);
		tempZ = PolyFormBezierCurve[ 2 ];
		tempZ.setBound(distinctKnotValue[ i ], distinctKnotValue[i + 1]);

		polynomialFormCurve[ 0 ][ i ] = tempX;
		polynomialFormCurve[ 1 ][ i ] = tempY;
		polynomialFormCurve[ 2 ][ i ] = tempZ;
	}

	delete[] BezierCurveList;

	return polynomialFormCurve;
	*/

	rg_INT NumOfNonZeroLengthKnotSpan = getNumOfNonZeroLengthKnotSpan();
	rg_DEGREE p = getOrder() - 1;

	rg_RBzCurve3D* RBezierCurveList = decomposeCurveIntoBezierSegment();

	rg_RationalPolynomial** polynomialFormCurve = new rg_RationalPolynomial* [ 3 ];

	rg_INT i = 0;
	for(i = 0;i < 3;i++)
		polynomialFormCurve[ i ] = new rg_RationalPolynomial[NumOfNonZeroLengthKnotSpan];

	rg_REAL* distinctKnotValue = getDistinctKnotValues();
	rg_Matrix P(p + 1, 4);
	rg_Matrix M_p = rg_CurveSurfaceFunc::bezierToPowerMatrix(p + 1);

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
		                               1/(distinctKnotValue[i + 1] - distinctKnotValue[ i ]),
			                           - distinctKnotValue[ i ] /(distinctKnotValue[i + 1] - distinctKnotValue[ i ]));
		rg_INT j = 0;
		for(j = 0;j < p + 1;j++)
		{
			P[ j ][ 0 ] = RBezierCurveList[ i ].getWeight( j ) * (RBezierCurveList[ i ].getCtrlPt( j )).getX(); // wX : X coordinate of ctrl pt
			P[ j ][ 1 ] = RBezierCurveList[ i ].getWeight( j ) * (RBezierCurveList[ i ].getCtrlPt( j )).getY(); // wY : Y coordinate of ctrl pt
			P[ j ][ 2 ] = RBezierCurveList[ i ].getWeight( j ) * (RBezierCurveList[ i ].getCtrlPt( j )).getZ(); // wZ : Z coordinate of ctrl pt
			P[ j ][ 3 ] = RBezierCurveList[ i ].getWeight( j )                                                   ; // w  : W weight
		}

		rg_Matrix coefficient = R_p * M_p * P; // Coefficient of power basis
										    // (p + 1) by 4 matrix
		
		rg_RationalPolynomial* temp = new rg_RationalPolynomial[ 3 ];

		for(j = 0;j < 3;j++)
		{
			temp[ j ].setDegree(p, p);
			for(rg_INT k = 0;k < p + 1;k++)
			{
				temp[ j ].setNumerCoefficient(k, coefficient[ k ][ j ]);
				temp[ j ].setDenomiCoefficient(k, coefficient[ k ][ 3 ]);
			}
			polynomialFormCurve[ j ][ i ] = temp[ j ];			
		}

		delete[] temp;
	}

	delete[] RBezierCurveList;

	return polynomialFormCurve;
}


void rg_NURBSplineCurve3D::getEntirePolynomialCurve(rg_REAL** & coeffOfXCoordinate,
												 rg_REAL** & coeffOfYCoordinate,
												 rg_REAL** & coeffOfZCoordinate,
												 rg_REAL** & coeffOfWeight      ) const
{
	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT*** allPossiblePath = new rg_INT** [p + 1]; // all possible paths from (p + 1) basic graphs
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;

	// find all the possible paths from (p + 1) basic graphs
	rg_INT* noOfAllPossiblePathsInEachGraph = new rg_INT[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		allPossiblePath[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, p - i);
		noOfAllPossiblePathsInEachGraph[ i ] = rg_MathFunc::combination(p, i);
	}
	//------------------------------------------

	// initialization of coefficient of piecewise polyomial curve
	coeffOfXCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfYCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfZCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfWeight      = new rg_REAL* [NoOfNonZeroLengthKnotSpan]; 

	rg_REAL* tempCoeff;

	for(i = 0;i < NoOfNonZeroLengthKnotSpan;i++)
	{
		coeffOfXCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfYCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfZCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfWeight[ i ]      = new rg_REAL[p + 1]; 
		for(rg_INT j = 0;j < p + 1;j++)
		{
			coeffOfXCoordinate[ i ][ j ] = 0.0;
			coeffOfYCoordinate[ i ][ j ] = 0.0;
			coeffOfZCoordinate[ i ][ j ] = 0.0;
			coeffOfWeight[ i ][ j ]      = 0.0;
		}
	}

	// generating the piecewise polynomial curve

	rg_Point3D* ctrlPt = getCtrlPts();
	rg_INT indexOfPolyCurve = 0;//, indexOfPolyCurveCoeff = 0;

	for(i = p;i <= m - p - 1 && indexOfPolyCurve < NoOfNonZeroLengthKnotSpan;i++)
	// loop for interior nonzero length knot spans
	{
		for(rg_INT j = i - p;j <= i;j++)
		// 
		{
			tempCoeff = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);

			for(rg_INT indexOfPolyCurveCoeff = 0;indexOfPolyCurveCoeff < p + 1;indexOfPolyCurveCoeff++)
			{
				coeffOfXCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (weights[ j ] * ctrlPt[ j ].getX() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfYCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (weights[ j ] * ctrlPt[ j ].getY() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfZCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (weights[ j ] * ctrlPt[ j ].getZ() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfWeight[indexOfPolyCurve][indexOfPolyCurveCoeff]      += (weights[ j ] * tempCoeff[indexOfPolyCurveCoeff]);
			}
			delete[] tempCoeff;			
		}
		indexOfPolyCurve++;
	}
	//------------------------------------------

	delete[] allPossiblePath;
	delete[] noOfAllPossiblePathsInEachGraph;
}

rg_NUBSplineCurve3D rg_NURBSplineCurve3D::getNumeratorInNUBS() const
{
    // Only the point is changed
    // Ohter data such as knot, order is same as this!!
    rg_NUBSplineCurve3D output(*this);

    rg_INT numOfCtrlPts=this->getNumOfCtrlPts();
    for ( rg_INT i=0 ; i < numOfCtrlPts;i++ )
    {
        rg_Point3D newCtrlPt= (this->getCtrlPt(i))*(this->getWeight(i));
        output.setCtrlPt(i, newCtrlPt);
    }
    return output;
}


rg_NUBSplineCurve3D rg_NURBSplineCurve3D::getDenominatorInNUBS() const // Onliy x coordinate has meaning 
                                                   // others coordinates are setted to 0
{
    // As there is not NUBs function with scalar value
    //    NUB curve is used. Only X-coordinate has meaning.
    // Only the point is changed
    // Ohter data such as knot, order is same as this!!
    rg_NUBSplineCurve3D output(*this);

    rg_INT numOfCtrlPts=this->getNumOfCtrlPts();
    for ( rg_INT i=0 ; i < numOfCtrlPts;i++ )
    {
        rg_Point3D newCtrlPt( (this->getWeight(i) ), 0.0, 0.0);
        output.setCtrlPt(i, newCtrlPt);
    }
    return output;
}

rg_INT   rg_NURBSplineCurve3D::getIndexOfNearestCtrlPt(const rg_Point3D& pt) const
{
    rg_INT numOfCtrlPts=this->getNumOfCtrlPts();

    if ( numOfCtrlPts < 1 )
    {
        return -1;
    }
    

    rg_REAL  min=(pt-(this->getCtrlPt(0))).magnitude();
    rg_INT   indexOfMin=0;
    for( rg_INT i=1; i < numOfCtrlPts; i++)
    {
        rg_REAL candidate=pt.distance(this->getCtrlPt(i));
        if ( rg_GT(min, candidate) )
        {
            min=candidate;
            indexOfMin=i;
        }
    }
    return indexOfMin;
}
rg_REAL  rg_NURBSplineCurve3D::getParameterOfNearestPt(const rg_Point3D& pt,
                                                 const rg_REAL&    distanceTol,
                                                 const rg_REAL&    cosAngleTol) const
{
    rg_NURBSplineCurve3D curve(*this);

    rg_NURBSCurveIntersector intersector;
    rg_dList<rg_REAL> characteristics=intersector.findCharParam(curve);

    // for case that there is no characteristic points
    if ( characteristics.getSize() == 0 ) 
    {
        rg_INT  numOfCtrlPts=this->getNumOfCtrlPts();
        rg_INT  indexOfNearestCtrlPt=this->getIndexOfNearestCtrlPt(pt);

        rg_REAL  seed=  ( (rg_REAL) indexOfNearestCtrlPt )
                    /( (rg_REAL) numOfCtrlPts-1.0     );

        rg_REAL output=getParameterOfNearestPt(pt,distanceTol,cosAngleTol,seed);
        
        return output;
    }

    rg_REAL startParameter=curve.getStartParameter();
    rg_REAL endParameter=curve.getEndParameter();
    //characteristics.addW(endParameter);
    
    rg_REAL output=0.5;
    rg_Point3D currentPt=curve.evaluatePt(output);
    

    characteristics.reset4Loop();
    characteristics.setNext4Loop();// first parameter is always start parater

    rg_REAL u=startParameter;
    rg_REAL prevU=0.0;

    while( characteristics.setNext4Loop() )
    {
        prevU=u;
        u=characteristics.getEntity();

        rg_REAL seed=(prevU+u)*0.5;

        rg_REAL sol=getParameterOfNearestPt(pt,distanceTol,cosAngleTol,seed);

        rg_Point3D newPt=curve.evaluatePt(sol);
        
        if ( rg_LT( newPt.distance(pt) , currentPt.distance(pt) ) )
        {
            output=sol;
            currentPt=newPt;
        }
    }

    return output;

}


rg_REAL  rg_NURBSplineCurve3D::getParameterOfNearestPt(const rg_Point3D& pt,
                                                 const rg_REAL&    distanceTol,
                                                 const rg_REAL&    cosAngleTol,
                                                 const rg_REAL&    seed) const
{
    // This problem is corresponding to the 
    //            the problem of finding  the u 
    //                    such as f(u)=C(u)'(C(u)-pt) = 0
    //                     if C(u) ~ this curve

    // find the seed parameter

    rg_INT  numOfCtrlPts=this->getNumOfCtrlPts();
    
    rg_INT  indexOfNearestCtrlPt=this->getIndexOfNearestCtrlPt(pt);
/*
    rg_REAL  seed=  ( (rg_REAL) indexOfNearestCtrlPt )
                /( (rg_REAL) numOfCtrlPts-1.0     );
*/  
    rg_Point3D seedPt=this->evaluatePt(seed);
    rg_Point3D seedTangent=(this->makeDerivative(seed)).getUnitVector();
    rg_Point3D projectDir=seedPt-pt;
    rg_Point3D unitProjectDir=projectDir.getUnitVector();
    
    // For case that seed is nearest point!!
    if (    rg_ZERO( projectDir.distance(), distanceTol ) 
         && rg_ZERO( unitProjectDir%seedTangent, cosAngleTol ) )
    {
        return seed;
    }

    rg_REAL  u=seed;// u is iterated 
    rg_REAL  prevU=seed;

    // [0]-> numerator
    // [1]-> 1st derivative of numerator
    // [2]-> 2nd derivative of numerator
    rg_NUBSplineCurve3D numerator[3];
    numerator[0]=this->getNumeratorInNUBS();
    numerator[1]=numerator[0].makeDerivative();
    numerator[2]=numerator[1].makeDerivative();

    // [0]-> denominator
    // [1]-> 1st derivative of denominator
    // [2]-> 2nd derivative of denominator
    rg_NUBSplineCurve3D denominator[3];
    denominator[0]=this->getDenominatorInNUBS();
    denominator[1]=denominator[0].makeDerivative();
    denominator[2]=denominator[1].makeDerivative();


    // compute C(u), C'(u), C''(u)
    // c[0] -> c(u)
    // c[1] -> c'(u)
    // c[2] -> c''(u)
    rg_REAL startParameter=rg_NURBSplineCurve3D::getStartParameter();
    rg_REAL endParameter=rg_NURBSplineCurve3D::getEndParameter();
    rg_FLAG doEndIteration=rg_FALSE;
    do
    {
        prevU=u;

        rg_Point3D A[3];
        A[0]=numerator[0].evaluatePt(u);
        A[1]=numerator[1].evaluatePt(u);
        A[2]=numerator[2].evaluatePt(u);

        rg_REAL w[3];
        w[0]=(denominator[0].evaluatePt(u)).getX();
        w[1]=(denominator[1].evaluatePt(u)).getX();
        w[2]=(denominator[2].evaluatePt(u)).getX();

        rg_Point3D c[3];
        c[0]=A[0]/w[0];
        c[1]=       A[1]-w[1]*c[0];
        c[1]=c[1] /( w[0] );
        c[2]=      A[2] - 2.0*w[1]*c[1]-w[2]*c[0];
        c[2]=c[2]/( w[0] );

        rg_Point3D projectDir=c[0]-pt;
        rg_Point3D unitProjectDir=projectDir.getUnitVector();


        u=prevU - ( c[1]%projectDir ) / ( c[2]%projectDir + pow(c[1].distance(),2) );

        u = ( rg_LT( u, startParameter) ) ? startParameter:u;
        u = ( rg_GT( u, endParameter) ) ? endParameter:u;

        rg_FLAG isCoincident    = rg_ZERO( projectDir.magnitude(), distanceTol);
        rg_FLAG isPerpendicular = rg_ZERO( unitProjectDir%(c[1].getUnitVector()), cosAngleTol);
        rg_FLAG isNoMoreChanged = rg_ZERO( (u-prevU)*c[1].distance() , distanceTol);
        doEndIteration =   ( isCoincident && isPerpendicular ) 
                         ||  isNoMoreChanged;

    }while ( doEndIteration == rg_FALSE );
   
    return u;
}
////    Set Functions. : March 13 1997
void rg_NURBSplineCurve3D::setWeight(const rg_INT &i, const rg_REAL &weight)
{
    weights[i] = weight;
}

void rg_NURBSplineCurve3D::setWeightVector(rg_REAL* weightVector)
{
    if (weights != rg_NULL)
        delete[] weights;

//  commented by rockhead!( 97.10.7 )
//  weights = weightVector;   
	
	rg_INT num = rg_BSplineCurve3D::getNumOfCtrlPts();
    
	weights = new rg_REAL [num];
    for (rg_INT i=0; i<num; i++)
        weights[i] = weightVector[i];

}

void rg_NURBSplineCurve3D::setNURBSCurve(const unsigned rg_INT &newID, 
                                      const rg_Planarity    &newPlanarity,
                                      const unsigned rg_INT &newOrder, 
                                      const rg_INT          &num, 
                                      rg_Point3D*             newControlP,
                                      const rg_REAL* const  newKnotVector,
                                      const rg_REAL* const  newWeightVector )
{
    rg_Curve::setID(newID);
    rg_Curve::setPlanarity(newPlanarity);
    
    rg_BSplineCurve3D::setOrder(newOrder);
    rg_BSplineCurve3D::setCtrlPts(num, newControlP);

    rg_NUBSplineCurve3D::setKnotVector(num+newOrder, newKnotVector);

    if (weights != rg_NULL)
        delete[] weights;

    weights = new rg_REAL [num];
    for (rg_INT i=0; i<num; i++)
        weights[i] = newWeightVector[i];
}
void rg_NURBSplineCurve3D::setNumOfCtrlPts(const rg_INT& newNumOfCtrlPts)
{
	rg_NUBSplineCurve3D::setNumOfCtrlPts(newNumOfCtrlPts);
    if (weights != rg_NULL)
        delete[] weights;

//  commented by rockhead!( 97.10.7 )
//  weights = weightVector;   
	
	weights = new rg_REAL [newNumOfCtrlPts];
    for (rg_INT i=0; i<newNumOfCtrlPts; i++)
        weights[i] = 1.0;
}
void rg_NURBSplineCurve3D::setCurve(const rg_RBzCurve3D& curve)
{
    const rg_INT order=curve.getDegree()+1;
    rg_BSplineCurve3D::setOrder(order);
    rg_Point3D* ctrlPts=curve.getCtrlPts();
    rg_REAL* weights=curve.getWeight();
    rg_BSplineCurve3D::setCtrlPts(order,ctrlPts);
    rg_NURBSplineCurve3D::setWeightVector(weights);

    delete[] ctrlPts;
    delete[] weights;

    const rg_INT numOfKnots=2*order;
    rg_REAL* knots=new rg_REAL[numOfKnots];

    for ( rg_INT i=0; i < order; i++ )
    {
        knots[i]=0.0;
        knots[order+i]=1.0;
    }
    rg_NUBSplineCurve3D::setKnotVector(numOfKnots,knots);
    delete[] knots;
}
////    BasisFunction & rg_Point3D Evaluating. : March 13 1997
//
//  March 13 1997 : Made
//  April  3 1997 : Modified
////////////////////////////////////////////////////////////////// 

rg_Point3D rg_NURBSplineCurve3D::evaluatePt(const rg_REAL &param) const
{
/*
    if( rg_NUBSplineCurve3D::getKnotVector() == rg_NULL )
        rg_NUBSplineCurve3D::setInitialKnotVector();
*/
/*
    rg_INT n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT order = (rg_INT) rg_BSplineCurve3D::getOrder();    

    rg_Point3D ptOnCurve(0.0, 0.0, 0.0);

    rg_REAL basisValue   = 0.;
    rg_REAL rationalPart = 0.;
 
    for (rg_INT i=0; i<n; i++)
    {
        basisValue = rg_NUBSplineCurve3D::evaluateBasisFunc(i, param, order);
        if ( rg_NZERO(basisValue) )
        {
            rationalPart += basisValue * weights[i];

            ptOnCurve += rg_BSplineCurve3D::getCtrlPt(i) 
                         * basisValue * weights[i];
        }
    }

    return ptOnCurve/rationalPart;
*/
    rg_INT  n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT  order = rg_BSplineCurve3D::getOrder();    

    rg_INDEX validKnot = 0;

    rg_INT first = order-1;
    rg_INT last  = n;
    rg_INT middle;

    while ( last >= first )
    {
        middle = (first + last)/2;

        if ( rg_LT(param, rg_NUBSplineCurve3D::getKnotValue(middle)) )
            last = middle;
        else if ( rg_GT(param, rg_NUBSplineCurve3D::getKnotValue(middle + 1)) )
            first = middle + 1;
        else
        {
			if ( middle != first )
			{
	            while ( rg_EQ( getKnotValue(middle), getKnotValue(middle+1) ) )
		            middle++;
	            validKnot = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
			else
			{
	            validKnot = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
            break;        }
    }

    rg_REAL* nonZeroBasis;

//    if ( rg_EQ(param, 1.0) )
//    {
//        nonZeroBasis = new rg_REAL[order];
//        for(rg_INDEX i=0; i<order-1; i++)
//            nonZeroBasis[i] = 0.0;
//        nonZeroBasis[order-1] = 1.0;
//    }
//    else 
//    {
        nonZeroBasis = evaluateMultiBasisFunc( validKnot, param, order );
//    }

    rg_Point3D ptOnCurve;
    rg_REAL  rationalPart = 0.;

    for (rg_INDEX i=validKnot-order+1; i<=validKnot; i++)
    {
        rationalPart += nonZeroBasis[i - validKnot + order -1]
                        * weights[i];
        ptOnCurve    += nonZeroBasis[i - validKnot + order -1]
                        * weights[i]
                        * rg_BSplineCurve3D::getCtrlPt(i); 
                          
    }

    delete [] nonZeroBasis;

    return ptOnCurve/rationalPart;

}

//  March 13 1997 : Made
//  April  3 1997 : Modified
////////////////////////////////////////////////////////////////// 
rg_Point3D* rg_NURBSplineCurve3D::evaluatePtsInEvenParameter(const rg_INT &numOfPtOnCurve) const 
{
    /*
    if( rg_NUBSplineCurve3D::getKnotVector() == rg_NULL )
        rg_NUBSplineCurve3D::setInitialKnotVector();
    */

    rg_Point3D* ptOnCurve = new rg_Point3D [numOfPtOnCurve];

    rg_REAL increment = 1. / (numOfPtOnCurve-1);
    rg_REAL u         = 0.;

    rg_INT n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT order = (rg_INT) rg_BSplineCurve3D::getOrder();

    rg_INDEX validKnot    = 0;
    rg_REAL* nonZeroBasis = rg_NULL;

    ptOnCurve[0] = rg_BSplineCurve3D::getCtrlPt(0);
    ptOnCurve[numOfPtOnCurve-1] = rg_BSplineCurve3D::getCtrlPt(n-1);

    u += increment;
    for (rg_INT i=1; i<numOfPtOnCurve-1; i++)
    {
//        evaluatedPoint[i] = evaluatePt(u);

        while ( rg_LT(u, rg_NUBSplineCurve3D::getKnotValue(validKnot))
                || rg_GE(u, rg_NUBSplineCurve3D::getKnotValue(validKnot + 1)) )          
        {
            validKnot++;
        }

        nonZeroBasis = evaluateMultiBasisFunc( validKnot, u, order );
        for (rg_INDEX j=validKnot-order+1; j<=validKnot; j++)
        {
            ptOnCurve[i] += rg_BSplineCurve3D::getCtrlPt(j) 
                                 * nonZeroBasis[j - validKnot + order -1];  
        }
        delete [] nonZeroBasis;

        u += increment;
    }

    return ptOnCurve;

/*	
    rg_Point3D tempPt;
	rg_REAL  basisValue   = 0.;
    rg_REAL  rationalPart = 0.;
    for (rg_INT i=0; i<numOfPtOnCurve; i++)
    {
        tempPt = 0.;

        rationalPart = 0.0;
        for (rg_INT j=0; j<n; j++)
        {
            basisValue = rg_NUBSplineCurve3D::evaluateBasisFunc(j, u, order);
            if ( rg_NZERO(basisValue) )
            {
                rationalPart += basisValue * weights[j];

                tempPt += rg_BSplineCurve3D::getCtrlPt(j) 
                          * basisValue * weights[j];
            }
        }
        ptOnCurve[i].setPoint(tempPt/rationalPart);
        u += increment;
    }

    return ptOnCurve;
*/
}
rg_Polyline2D rg_NURBSplineCurve3D::makePolyline2DInEvenParameter(const rg_INT &numOfPts) const
{
    if ( numOfPts < 1 )
    {
        return rg_Polyline2D();
    }
    rg_Point3D* tPts=evaluatePtsInEvenParameter(numOfPts);
    rg_Point2D*  pts=new rg_Point2D[numOfPts];
    for( rg_INT i=0; i < numOfPts; i++ )
    {
        pts[i]=tPts[i].evaluatePt2D();
    }
    
    rg_Polyline2D output(numOfPts, pts);
    delete[] pts;
    delete[] tPts;

    return output;
}

rg_Polyline2D rg_NURBSplineCurve3D::makePolyline2DConsideringKnots(const rg_INT &numOfPts) const
{
    if ( isNull() || numOfPts < 1 )
    {
        return rg_Polyline2D();
    }
    rg_Point2D* pts=new rg_Point2D[numOfPts];

    rg_INT numOfSegments=getNumOfNonZeroLengthKnotSpan();
    rg_REAL* parametersOfSegments=getDistinctKnotValues();
    rg_INT   numOfPtsOfSegment=(numOfPts-1)/(numOfSegments);
    rg_INT currentIndex=0;

    for( rg_INT i=0; i < numOfSegments-1; i++ )
    {
        rg_REAL current=parametersOfSegments[i];
        rg_REAL step=  (parametersOfSegments[i+1]-parametersOfSegments[i])
                      /(rg_REAL)numOfPtsOfSegment;
        for( rg_INT j=0 ; j < numOfPtsOfSegment; j++ )
        {
            rg_Point3D pt=evaluatePt(current);
            pts[currentIndex]=pt.evaluatePt2D();
            currentIndex++;
            current+=step;
        }
    }
    rg_INT remaindedNumOfPts=numOfPts-currentIndex;
    rg_REAL current=parametersOfSegments[numOfSegments-1];
    rg_REAL step = (parametersOfSegments[numOfSegments]-current)
                   /(remaindedNumOfPts-1);
    for( rg_INT j=0; j < remaindedNumOfPts; j++)
    {
        rg_Point3D pt=evaluatePt(current);
        pts[currentIndex]=pt.evaluatePt2D();
        currentIndex++;
        current+=step;
    }

    rg_Polyline2D output(numOfPts, pts);

    delete[] parametersOfSegments;
    delete[] pts;
    return output;
}

rg_Polyline3D rg_NURBSplineCurve3D::makePolyline3DInEvenParameter(const rg_INT &numOfPts) const
{
    if ( numOfPts < 1 )
    {
        return rg_Polyline3D();
    }
    rg_Point3D* pts=evaluatePtsInEvenParameter(numOfPts);
    
    rg_Polyline3D output(numOfPts, pts);
    delete[] pts;

    return output;
}


////    Derivative rg_Curve & Curvature
rg_NURBSplineCurve3D rg_NURBSplineCurve3D::makeHodograph() const
{
	const rg_INT numOfSegments=getNumOfNonZeroLengthKnotSpan();
	rg_RBzCurve3D* segments=decomposeCurveIntoBezierSegment();
	const rg_INT numOfDistinctKnots=numOfSegments+1;
	const rg_INT numOfKnots=getNumberOfKnotValues();

	// find orginal knot vector information
    rg_REAL* distinctKnots= new rg_REAL[numOfDistinctKnots];
	rg_INT* knotMultiplicities= new rg_INT[numOfDistinctKnots];
	distinctKnots[0]=getStartParameter();
	knotMultiplicities[0]=order;
	distinctKnots[numOfDistinctKnots-1]=getEndParameter();
	knotMultiplicities[numOfDistinctKnots-1]=order;	
    rg_INT i = 1;      //  i traces the array multiplicity.
    rg_INT j = order;    //  j traces knotVector.
	while( i < numOfDistinctKnots-1 )
	{
		rg_INT s = 1;
		while(   (j < numOfKnots - order) 
			  && rg_EQ(getKnotValue(j), getKnotValue(j+1)) )
		{
			s++;
			j++;
		}
		distinctKnots[i]= getKnotValue(j);  // replicated knot value in the interior knot span
		knotMultiplicities[i] = s;                // multiplicity of replicated knot value in the interior knot span
		i++;
        j++;
	}	

	// obtain hodographs of Bezier segments
	const rg_INT orderOfHodographs=2*order-1;
	rg_RBzCurve3D* hodographsOfSegments=new rg_RBzCurve3D[numOfSegments];
	for (i=0 ; i < numOfSegments;i++ )
	{
		hodographsOfSegments[i]= segments[i].makeHodograph();
	}

	// compose hodographs of Bezier segments into NURBsplines
	rg_NURBSplineCurve3D hodograph(orderOfHodographs);
	hodograph.setNumOfCtrlPts(numOfSegments*orderOfHodographs);

	i=0; // i for tarcing ctrlPts
	for( j=0; j < numOfSegments; j++ )// j for tracing segments
	{
		for ( int k=0; k < orderOfHodographs; k++)
		{

			const rg_Point3D pt=hodographsOfSegments[j].getCtrlPt(k) ;
			const rg_REAL    weight=hodographsOfSegments[j].getWeight(k) ;
			hodograph.setCtrlPt(i,pt);
			hodograph.setWeight(i,weight);
			i++;
		}
	}

	hodograph.setInitialKnotVector();
	i=0; // i for tarcing knot vector
	for( j=0; j < numOfDistinctKnots; j++ ) // j for tracing origianl discrete knot 
	{
		for ( int k=0; k < orderOfHodographs; k++ )
		{
			hodograph.setKnotValue(i,distinctKnots[j]);
			i++;
		}
	}


	for ( i=1; i < numOfDistinctKnots-1; i++ )
	{
		for( j=0; j < order-1- knotMultiplicities[i]; j++ )
		{
			hodograph.removeRedundantKnot(distinctKnots[i]);
		}
	}




	delete[] segments;
	delete[]  hodographsOfSegments;
	delete[] distinctKnots;
	delete[] knotMultiplicities;


	return hodograph;


}

rg_NURBSplineCurve3D rg_NURBSplineCurve3D::makeRedundantHodograph() const
{
	const rg_INT numOfSegments=getNumOfNonZeroLengthKnotSpan();
	rg_RBzCurve3D* segments=decomposeCurveIntoBezierSegment();
	const rg_INT numOfDistinctKnots=numOfSegments+1;
	const rg_INT numOfKnots=getNumberOfKnotValues();

	// find orginal knot vector information
    rg_REAL* distinctKnots= new rg_REAL[numOfDistinctKnots];
	rg_INT* knotMultiplicities= new rg_INT[numOfDistinctKnots];
	distinctKnots[0]=getStartParameter();
	knotMultiplicities[0]=order;
	distinctKnots[numOfDistinctKnots-1]=getEndParameter();
	knotMultiplicities[numOfDistinctKnots-1]=order;	
    rg_INT i = 1;      //  i traces the array multiplicity.
    rg_INT j = order;    //  j traces knotVector.
	while( i < numOfDistinctKnots-1 )
	{
		rg_INT s = 1;
		while(   (j < numOfKnots - order) 
			  && rg_EQ(getKnotValue(j), getKnotValue(j+1)) )
		{
			s++;
			j++;
		}
		distinctKnots[i]= getKnotValue(j);  // replicated knot value in the interior knot span
		knotMultiplicities[i] = s;                // multiplicity of replicated knot value in the interior knot span
		i++;
        j++;
	}	

	// obtain hodographs of Bezier segments
	const rg_INT orderOfHodographs=2*order-1;
	rg_RBzCurve3D* hodographsOfSegments=new rg_RBzCurve3D[numOfSegments];
	for (i=0 ; i < numOfSegments;i++ )
	{
		hodographsOfSegments[i]= segments[i].makeHodograph();
	}

	// compose hodographs of Bezier segments into NURBsplines
	rg_NURBSplineCurve3D hodograph(orderOfHodographs);
	hodograph.setNumOfCtrlPts(numOfSegments*orderOfHodographs);

	i=0; // i for tarcing ctrlPts
	for( j=0; j < numOfSegments; j++ )// j for tracing segments
	{
		for ( int k=0; k < orderOfHodographs; k++)
		{

			const rg_Point3D pt=hodographsOfSegments[j].getCtrlPt(k) ;
			const rg_REAL    weight=hodographsOfSegments[j].getWeight(k) ;
			hodograph.setCtrlPt(i,pt);
			hodograph.setWeight(i,weight);
			i++;
		}
	}

	hodograph.setInitialKnotVector();
	i=0; // i for tarcing knot vector
	for( j=0; j < numOfDistinctKnots; j++ ) // j for tracing origianl discrete knot 
	{
		for ( int k=0; k < orderOfHodographs; k++ )
		{
			hodograph.setKnotValue(i,distinctKnots[j]);
			i++;
		}
	}


	delete[] segments;
	delete[]  hodographsOfSegments;
	delete[] distinctKnots;
	delete[] knotMultiplicities;


	return hodograph;


}	


//  April 7 1997 : made.
rg_Point3D rg_NURBSplineCurve3D::makeDerivative(const rg_REAL &u) const
{
    //  C(u) : NURBS rg_Curve.
    //
    //         A(u)              A'(u) - w'(u)C(u)
    //  C(u) = ---- ,   C'(u) = -------------------
    //         w(u)                    w(u)
    //
    //___________________________________________________

    rg_INT n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT order = (rg_INT) rg_BSplineCurve3D::getOrder();

    rg_NUBSplineCurve3D forBasisFunc(n-1);
    forBasisFunc.rg_BSplineCurve3D::setOrder(order-1);

    rg_REAL* tempKnot = new rg_REAL [n+order-2];
    rg_INT i = 0;
	for (i=0; i<n+order-2; i++)
        tempKnot[i] = rg_NUBSplineCurve3D::getKnotValue(i+1);
    
    forBasisFunc.rg_NUBSplineCurve3D::setKnotVector(n+order-2, tempKnot);                               
    delete [] tempKnot;

    rg_Point3D AprimeU;        // A'(u)
    rg_Point3D CU;             // C(u)
    rg_REAL  wprimeU = 0.0;  // w'(u)
    rg_REAL  wU      = 0.0;  // w(u)

    rg_REAL  knotDelta = 0.0;
    rg_REAL  basisFunc = 0.0;
    for (i=0; i<n-1; i++)
    {
//	commented by rockhead!!!( 97.10.7 )
/*
        knotDelta = rg_NUBSplineCurve3D::getKnotValue(i+order-1)
                    - rg_NUBSplineCurve3D::getKnotValue(i+1);
*/
        knotDelta = rg_NUBSplineCurve3D::getKnotValue(i+order)
                    - rg_NUBSplineCurve3D::getKnotValue(i+1);

        basisFunc = forBasisFunc.rg_NUBSplineCurve3D::evaluateBasisFunc(i, u, order-1);
        if ( rg_NZERO(basisFunc) )
        {
            AprimeU += basisFunc
                       * ( weights[i+1]*rg_BSplineCurve3D::getCtrlPt(i+1)
                         - weights[i]*rg_BSplineCurve3D::getCtrlPt(i) )
                       / knotDelta;

            wprimeU += basisFunc
                       * ( weights[i+1] - weights[i] )
                       / knotDelta;
        }
    }
    AprimeU = (order-1)*AprimeU;
    wprimeU = (order-1)*wprimeU;

    CU = evaluatePt(u);

    for (i=0; i<n; i++)
        wU += rg_NUBSplineCurve3D::evaluateBasisFunc(i, u, order)*weights[i];

    rg_Point3D derivativeAtU = ( AprimeU - wprimeU * CU ) / wU;

    return derivativeAtU;
}

////    represent the basic geometry with NURBS
////    represent the basic geometry with NURBS

// formArc
//   start and end is arranged in antri-clockwise
//   localZ = vector from arc to eye point

void rg_NURBSplineCurve3D::formArc( const rg_Point3D& center,
                                 const rg_Point3D& start,
                                 const rg_Point3D& end,
                                 const rg_Point3D& localZ  )
{
    const rg_Point3D axisOfArc      = (start-center)*(end-center);
    rg_REAL        cosOfAngle     = ( (start-center).getUnitVector() ) % ( (end-center).getUnitVector() ) ;
    rg_REAL        angle          = acos( cosOfAngle );
    rg_REAL        radius         = ( start-center ).magnitude();
    //const rg_REAL  phi            = 3.141592654;
	const rg_REAL phi			  = rg_PI;	 

    if ( rg_LE( axisOfArc%localZ, 0.0 ) )
    {
        angle=2*phi-angle;
    }

    // get number of partion to each arc segment's angle between min (angle,phi/4) and phi/21
    const rg_REAL ratioOfPartition=angle/phi*2;
    rg_INT        numOfPartition=(rg_INT) ratioOfPartition;
    numOfPartition = ( rg_EQ( ratioOfPartition,numOfPartition ) ) ? numOfPartition:numOfPartition+1;

    const rg_INT  newNumOfPoints = 2*numOfPartition+1;
    const rg_INT  newOrder       = 3;
    const rg_REAL stepAngle=angle/numOfPartition;
    const rg_REAL halfStepAngle=stepAngle/2;
    const rg_REAL cosOfHalfStepAngle=cos( halfStepAngle);

    rg_Point3D* newCtrlPts = new rg_Point3D[newNumOfPoints];
    rg_REAL*  newWeights = new rg_REAL[newNumOfPoints];
    rg_REAL*  newKnots   = new rg_REAL[newNumOfPoints+newOrder];

    rg_INT i = 0;
	for( i=0 ; i < newOrder; i++)
    {
        newKnots[0]=0.0;
        newKnots[newNumOfPoints+i]=1.0;
    }

    rg_TMatrix3D transform;
    transform.rotateArbitraryAxis( localZ, stepAngle );

    rg_Point3D  front(start);
    rg_Point3D  rear;
    rg_Point3D  bisectVector;
    rg_Point3D  middle;
 
    rg_INT ptsIndex=0;
    rg_INT knotIndex=1;

    for(i=0; i < numOfPartition; i++ )
    {
        rear  = transform*(front-center)+center;
        bisectVector= (front+rear)/2 - center;
        middle = radius/cosOfHalfStepAngle*bisectVector.getUnitVector()+ center; 

        newCtrlPts[ptsIndex]    = front;
        newCtrlPts[ptsIndex+1]  = middle;
        newWeights[ptsIndex]    = 1.0;
        newWeights[ptsIndex+1]  = cosOfHalfStepAngle;
        ptsIndex += 2;
        newKnots[knotIndex]   = i/(rg_REAL)numOfPartition;
        newKnots[knotIndex+1] = i/(rg_REAL)numOfPartition;
        knotIndex += 2;
        front=rear;
    }
    newCtrlPts[newNumOfPoints-1]=end;
    newWeights[newNumOfPoints-1]=1.0;

    rg_BSplineCurve3D::setOrder(newOrder);
    rg_BSplineCurve3D::setCtrlPts(newNumOfPoints,newCtrlPts);
    rg_NUBSplineCurve3D::setKnotVector(newNumOfPoints+newOrder,newKnots);
    setWeightVector(newWeights);
    delete[] newCtrlPts;
    delete[] newKnots;
    delete[] newWeights;
}

void  rg_NURBSplineCurve3D::formEllipticArc( const rg_EllipticArc3D& ellipticArc)
{
    const rg_REAL pi=4.0*atan(1.0);
    rg_REAL nStartAngle=rg_GeoFunc::getNormalizedAngle(ellipticArc.getStartAngle());
    rg_REAL nEndAngle=rg_GeoFunc::getNormalizedAngle(ellipticArc.getEndAngle());

    // make 0 <= start < end <  4*pi
    rg_REAL start=nStartAngle;
    rg_REAL end=nEndAngle;

    if ( rg_GE(nStartAngle, nEndAngle) )
    {
        end=nEndAngle+2*pi;
    }


    // find angles and numOfSegs
    rg_REAL splittedAngles[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    rg_INT  numOfSegs=0;
    
    splittedAngles[0]=start;
    rg_REAL current=ceil(start/(pi/2.0))*(pi/2.0);

    if ( rg_EQ( current , start) )
    {
        current=current+pi/2;
    }


    while( rg_LT(current,end ) )
    {
        numOfSegs++;
        splittedAngles[numOfSegs]=current;
        current=current+pi/2;
    }

    numOfSegs++;
    splittedAngles[numOfSegs]=end;
    
    //
    const rg_INT  newNumOfPoints = 2*numOfSegs+1;
    const rg_INT  newOrder       = 3;
    rg_Point3D* newCtrlPts = new rg_Point3D[newNumOfPoints];
    rg_REAL*  newWeights = new rg_REAL[newNumOfPoints];
    rg_REAL*  newKnots   = new rg_REAL[newNumOfPoints+newOrder];

    newKnots[0]=0.0;
//    newKnots[1]=0.0;
//    newKnots[2]=0.0;

    //rg_INT knotStart=3;
    rg_INT knotIndex=1;
    rg_INT ptsIndex=0;
    rg_INT i = 0;
	for( i=0; i < numOfSegs; i++ )
    {
        rg_REAL sAngleOfCurrent=splittedAngles[i];
        rg_REAL eAngleOfCurrent=splittedAngles[i+1];
        rg_REAL middleAngle=(sAngleOfCurrent+eAngleOfCurrent)*0.5;
        rg_EllipticArc3D currentArc(ellipticArc);
        currentArc.setStartAngle(sAngleOfCurrent);
        currentArc.setEndAngle(eAngleOfCurrent);
        rg_RBzCurve3D currentSeg;
        currentSeg.makeConic( currentArc.evaluatePt(sAngleOfCurrent)
                             ,currentArc.evaluateDerivative(sAngleOfCurrent)
                             ,currentArc.evaluatePt(eAngleOfCurrent)
                             ,currentArc.evaluateDerivative(eAngleOfCurrent)
                             ,currentArc.evaluatePt(middleAngle) );
/*
        newCtrlPts[2*i]  =currentSeg.getCtrlPt(0);
        newCtrlPts[2*i+1]=currentSeg.getCtrlPt(1);
        newWeights[2*i]  =currentSeg.getWeight(0);
        newWeights[2*i+1]=currentSeg.getWeight(1);

        newKnots[knotStart+i]  =(eAngleOfCurrent-start)/(end-start);
*/
        newCtrlPts[ptsIndex]  =currentSeg.getCtrlPt(0);
        newWeights[ptsIndex]  =currentSeg.getWeight(0);

        newCtrlPts[ptsIndex+1]=currentSeg.getCtrlPt(1);
        newWeights[ptsIndex+1]=currentSeg.getWeight(1);
        
        ptsIndex=ptsIndex+2;

        newKnots[knotIndex]    =(sAngleOfCurrent-start)/(end-start);
        newKnots[knotIndex+1]  =(sAngleOfCurrent-start)/(end-start);

        knotIndex=knotIndex+2;

    }
    newCtrlPts[newNumOfPoints-1]=ellipticArc.evaluatePt(ellipticArc.getEndAngle());
    newWeights[newNumOfPoints-1]=1.0;
    newKnots[newNumOfPoints+newOrder-1]=1.0;
    newKnots[newNumOfPoints+newOrder-2]=1.0;
    newKnots[newNumOfPoints+newOrder-3]=1.0;



    rg_BSplineCurve3D::setOrder(newOrder);
    rg_BSplineCurve3D::setCtrlPts(newNumOfPoints,newCtrlPts);
    rg_NUBSplineCurve3D::setKnotVector(newNumOfPoints+newOrder,newKnots);
    setWeightVector(newWeights);
    delete[] newCtrlPts;
    delete[] newKnots;
    delete[] newWeights;

}



void rg_NURBSplineCurve3D::formPolyline3D( const rg_INT      &numOfPts, 
                                                 rg_Point3D*  polygons)
{
    const rg_INT newOrder=2;
    const rg_INT newNumOfPts=numOfPts;
    const rg_INT newNumOfKnots= newOrder+newNumOfPts;

    rg_Point3D* newCtrlPts=new rg_Point3D[newNumOfPts];
    rg_REAL*  newKnots=new rg_REAL[newNumOfKnots];
    rg_REAL*  newWeights=new rg_REAL[newNumOfPts];

    newKnots[0]=0.0;
    newKnots[newNumOfKnots-1]=1.0;
    rg_INT i = 0;
	for( i=0; i < newNumOfPts; i++ )
    {
        newKnots[1+i]=(rg_REAL)i/(rg_REAL)(newNumOfPts-1);
        newCtrlPts[i]=polygons[i];
        newWeights[i]=1.0;
    }

    rg_BSplineCurve3D::setOrder(newOrder);
    rg_BSplineCurve3D::setCtrlPts(newNumOfPts,newCtrlPts);
    rg_NUBSplineCurve3D::setKnotVector(newNumOfPts+newOrder,newKnots);
    setWeightVector(newWeights);

    delete[] newCtrlPts;
    delete[] newKnots;
    delete[] newWeights;

}

void rg_NURBSplineCurve3D::formLine( const rg_Point3D& start,
                                  const rg_Point3D& end   )
{
    const rg_INT newOrder      = 2;
    const rg_INT newNumOfPts   = 2;
    const rg_INT newNumOfKnots = newOrder + newNumOfPts;

    rg_Point3D* newCtrlPts=new rg_Point3D[newNumOfPts];
    rg_REAL*  newKnots  =new rg_REAL[newNumOfKnots];
    rg_REAL*  newWeights=new rg_REAL[newNumOfPts];

    newCtrlPts[0]=start;
    newCtrlPts[1]=end;
    newKnots[0]=0.0;

    rg_INT i = 0;
	for( i=0; i < 2; i++ )
    {
        newWeights[i]   = 1.0;
        newKnots[2*i]   = (rg_REAL)i;
        newKnots[2*i+1] = (rg_REAL)i;
    }

    rg_BSplineCurve3D::setOrder(newOrder);
    rg_BSplineCurve3D::setCtrlPts(newNumOfPts,newCtrlPts);
    rg_NUBSplineCurve3D::setKnotVector(newNumOfPts+newOrder,newKnots);
    setWeightVector(newWeights );

    delete[] newCtrlPts;
    delete[] newKnots;
    delete[] newWeights;
}

rg_REAL  rg_NURBSplineCurve3D::findParameterOfNearestPt(const rg_Point3D& pt,
							 	 					    const rg_REAL&    distanceTol,
													    const rg_REAL&    cosAngleTol,
													    const rg_REAL&    start,
													    const rg_REAL&    end) const
{
    rg_NURBSplineCurve3D curve(*this);

    rg_NURBSCurveIntersector intersector;
    rg_dList<rg_REAL> characteristics=intersector.findCharParam(curve);

    // for case that there is no characteristic points
    if ( characteristics.getSize() == 0 ) 
    {
        rg_INT  numOfCtrlPts=this->getNumOfCtrlPts();
        rg_INT  indexOfNearestCtrlPt=this->findIndexOfNearestCtrlPt(pt);

        rg_REAL seed=  ( (rg_REAL) indexOfNearestCtrlPt )
                    /( (rg_REAL) numOfCtrlPts-1.0     );

        rg_REAL output=findParameterOfNearestPt(pt,distanceTol,cosAngleTol,seed,0.0,1.0);
        
        return output;
    }

    rg_REAL startParameter=curve.getStartParameter();
    rg_REAL endParameter=curve.getEndParameter();
    //characteristics.addW(endParameter);
    
    rg_REAL output=0.5;
    rg_Point3D currentPt=curve.evaluatePt(output);
    

    characteristics.reset4Loop();
    characteristics.setNext4Loop();// first parameter is always start parater

    rg_REAL u=startParameter;
    rg_REAL prevU=0.0;

    while( characteristics.setNext4Loop() )
    {
        prevU=u;
        u=characteristics.getEntity();
		rg_NURBSplineCurve3D currentCurve;
		currentCurve=curve.evaluateCurveSegment(prevU,u);
		rg_BoundingBox3D box=currentCurve.makeBoundingBox3D();
        box=box.evaluateRelativeOffset(1.5);
		if ( box.doContain(pt)  )
		{
			rg_REAL seed=(prevU+u)*0.5;
			rg_REAL sol=findParameterOfNearestPt(pt,distanceTol,cosAngleTol,seed,prevU,u);

			rg_Point3D newPt=curve.evaluatePt(sol);
       
			if ( rg_LT( newPt.distance(pt) , currentPt.distance(pt) ) )
			{
				output=sol;
				currentPt=newPt;
			}
		}
    }

    return output;

}


rg_REAL  rg_NURBSplineCurve3D::findParameterOfNearestPt(const rg_Point3D& pt,
                                                        const rg_REAL&    distanceTol,
                                                        const rg_REAL&    cosAngleTol,
                                                        const rg_REAL&    seed,
				    				 				    const rg_REAL&    start,
					    			 				    const rg_REAL&    end) const
{ 
    // This problem is corresponding to the 
    //            the problem of finding  the u 
    //                    such as f(u)=C(u)'(C(u)-pt) = 0
    //                     if C(u) ~ this curve

    // find the seed parameter

    rg_INT  numOfCtrlPts=this->getNumOfCtrlPts();
    
    rg_INT  indexOfNearestCtrlPt=findIndexOfNearestCtrlPt(pt);
/*
    REAL  seed=  ( (REAL) indexOfNearestCtrlPt )
                /( (REAL) numOfCtrlPts-1.0     );
*/  
    rg_Point3D seedPt=this->evaluatePt(seed);
    rg_Point3D seedTangent=(this->makeDerivative(seed)).getUnitVector();
    rg_Point3D projectDir=seedPt-pt;
    rg_Point3D unitProjectDir=projectDir.getUnitVector();
    
    // For case that seed is nearest point!!
    if (    rg_ZERO( projectDir.distance(), distanceTol ) 
         && rg_ZERO( unitProjectDir%seedTangent, cosAngleTol ) )
    {
        return seed;
    }

    rg_REAL  u=seed;// u is iterated 
    rg_REAL  prevU=seed;

    // [0]-> numerator
    // [1]-> 1st derivative of numerator
    // [2]-> 2nd derivative of numerator
    rg_NUBSplineCurve3D numerator[3];
    numerator[0]=this->getNumeratorInNUBS();
    numerator[1]=numerator[0].makeDerivative();
    numerator[2]=numerator[1].makeDerivative();

    // [0]-> denominator
    // [1]-> 1st derivative of denominator
    // [2]-> 2nd derivative of denominator
    rg_NUBSplineCurve3D denominator[3];
    denominator[0]=this->getDenominatorInNUBS();
    denominator[1]=denominator[0].makeDerivative();
    denominator[2]=denominator[1].makeDerivative();


    // compute C(u), C'(u), C''(u)
    // c[0] -> c(u)
    // c[1] -> c'(u)
    // c[2] -> c''(u)
    rg_REAL startParameter=rg_NUBSplineCurve3D::getStartParameter();
    rg_REAL endParameter=rg_NUBSplineCurve3D::getEndParameter();
    rg_FLAG doEndIteration=rg_FALSE;
    do
    {
        prevU=u;

        rg_Point3D A[3];
        A[0]=numerator[0].evaluatePt(u);
        A[1]=numerator[1].evaluatePt(u);
        A[2]=numerator[2].evaluatePt(u);

        rg_REAL w[3];
        w[0]=(denominator[0].evaluatePt(u)).getX();
        w[1]=(denominator[1].evaluatePt(u)).getX();
        w[2]=(denominator[2].evaluatePt(u)).getX();

        rg_Point3D c[3];
        c[0]=A[0]/w[0];
        c[1]=       A[1]-w[1]*c[0];
        c[1]=c[1] /( w[0] );
        c[2]=      A[2] - 2.0*w[1]*c[1]-w[2]*c[0];
        c[2]=c[2]/( w[0] );

        rg_Point3D projectDir=c[0]-pt;
        rg_Point3D unitProjectDir=projectDir.getUnitVector();


        u=prevU - ( c[1]%projectDir ) / ( c[2]%projectDir + pow(c[1].distance(),2) );

        u = ( rg_LT( u, start) ) ? (prevU+start)*0.5:u;
        u = ( rg_GT( u, end) ) ? (prevU+end)*0.5:u;

        rg_FLAG isCoincident    = rg_ZERO( projectDir.magnitude(), distanceTol);
        rg_FLAG isPerpendicular = rg_ZERO( unitProjectDir%(c[1].getUnitVector()), cosAngleTol);
        rg_FLAG isNoMoreChanged = rg_ZERO( (u-prevU)*c[1].distance() , distanceTol);
        doEndIteration =   ( isCoincident && isPerpendicular ) 
                         ||  isNoMoreChanged;

    }while ( doEndIteration == rg_FALSE );
   
    return u;
}
void   rg_NURBSplineCurve3D::reverseTrace()
{
    rg_INT halfNumOfCtrlPts=numOfCtrlPts/2;
    rg_INT i = 0;
	for( i=0; i < halfNumOfCtrlPts; i++ )
    {
        const rg_INT before=i;
        const rg_INT after=numOfCtrlPts-1-i;
        const rg_Point3D temp=ctrlPts[before];
        ctrlPts[before]=ctrlPts[after];
        ctrlPts[after]=temp;

        const rg_REAL tWeight=weights[before];
        weights[before]=weights[after];
        weights[after]=tWeight;
    }
    const rg_INT numOfKnots=getNumberOfKnotValues();
    rg_INT halfNumOfKnots=numOfKnots/2;
    for( i=0; i < halfNumOfKnots; i++ )
    {
        const rg_INT before=i;
        const rg_INT after=numOfCtrlPts-1-i;
        const rg_REAL temp=knotVector[before];
        knotVector[before]=knotVector[after];
        knotVector[after]=temp;
    }

}

//  March 26 1997
void rg_NURBSplineCurve3D::curveInterpolation(const rg_INT          &n, 
                                           const rg_Point3D* const  ptsPassedThrough,
                                           const rg_INT          &order,
                                           rg_REAL*               param,
                                           rg_REAL*               weightVector)
{
    if (param == rg_NULL)
        //param = rg_GeoFunc::rg_GeoFunc::chordLength(n, ptsPassedThrough);
        param = rg_GeoFunc::chordLength(n, ptsPassedThrough);
    if (weightVector == rg_NULL)
    {
        weightVector = new rg_REAL [n+2];
        for (rg_INT i=0; i<n+2; i++)
            weightVector[i] = 1.0;
    }

	rg_Point3D m0(ptsPassedThrough[1] - ptsPassedThrough[0]);
	rg_Point3D mL(ptsPassedThrough[n-1] - ptsPassedThrough[n-2]);

	rg_Point3D *data = new rg_Point3D[n];
	rg_INT i = 0;
	for (i = 0; i < n; i++) 
		data[i] = ptsPassedThrough[i];
	
	rg_Point3D *d = new rg_Point3D[n];

    rg_REAL **mat = rg_NUBSplineCurve3D::makeMatrix(data, param,  n-1, m0, mL);

    rg_INT bandSize=1;
	rg_BandedMatrix<rg_Point3D> bm2(n, bandSize);

	bm2.findVariable(mat, d, data);

    //  Set order.
    rg_BSplineCurve3D::setOrder(order);
    //  Set control polygon. 
	rg_BSplineCurve3D::setNumOfCtrlPts(n + 2);

	rg_BSplineCurve3D::setCtrlPt(0, ptsPassedThrough[0]);
	for (i=0; i<n; i++)
		rg_BSplineCurve3D::setCtrlPt(i+1, d[i]);
	rg_BSplineCurve3D::setCtrlPt(n+1, ptsPassedThrough[n-1]);
	
    //  Set knot vector.
	rg_INT numOfKnot = n + 2 + order;

    rg_REAL* knotVector = new rg_REAL [numOfKnot];
	for (i = 0; i<(order-1); i++)
	{
		knotVector[i] = 0.;
		knotVector[numOfKnot - 1 - i] = 1.;
	}

	for (i = order-1; i<(numOfKnot-(order -1)); i++)
		knotVector[i] = param[i - (order-1)];
	
    rg_NUBSplineCurve3D::setKnotVector(numOfKnot, knotVector);

    //  Set weight.
    weights = weightVector;

    delete [] knotVector;
	delete [] d;
    delete [] data;
	for(i=0; i<=2*bandSize; i++)
		delete [] mat[i];
	delete [] mat;
}

void rg_NURBSplineCurve3D::knotInsertion( const rg_REAL &insertingKnot )
{
//    if ( rg_LT(insertingKnot, 0.) || rg_GT(insertingKnot, 1.0) )
    if ( !rg_BTOR(0., insertingKnot, 1.0) )
        return;

    //if (knotVector == rg_NULL)
        //setInitialKnotVector();

    rg_INT n     = rg_BSplineCurve3D::getNumOfCtrlPts();
	rg_INT order = rg_BSplineCurve3D::getOrder();

    //  1.  Find the k-th knot-span, [U_k, U_k+1) lying insertingKnot.
    rg_INT k = 0;
	for (k=order-1; k<n; k++)
    {
        if (    rg_GE( insertingKnot, getKnotValue(k)) 
             && rg_LT( insertingKnot, getKnotValue(k+1)) )
            break;
    }	

    //  2.  Construct a new knot vector
    rg_REAL* newKnotVector = new rg_REAL[n+order+1];

    rg_INT i = 0;
	for (i=0; i<=k; i++)
        newKnotVector[i] = getKnotValue(i);
    newKnotVector[k+1] = insertingKnot;
    for (i=k+2; i<n+order+1 ; i++)
        newKnotVector[i] = getKnotValue(i-1);

    //  3.  Construct a new control polygon influenced by knot insertion. 
    //
    //  3.1 Store control points that aren't influenced by knot insertion
    //      to array for the new control polygon.
    rg_Point3D* newCtrlPlygn = new rg_Point3D[n+1];
	rg_REAL*  newWeight_vector = new rg_REAL[n+1];

    rg_INT j = 0;
	for (j=0; j<=(k-order+1); j++)
    {
        newCtrlPlygn[j] = rg_BSplineCurve3D::getCtrlPt(j);
		newWeight_vector[j] = getWeight(j);
    }
    //  3.2 Store (order-1) control points influenced by knot insertion
    //      to array. 
    for (j=(k-order+2); j<=k; j++)
    {
        rg_REAL alpha = ( insertingKnot - getKnotValue(j) )
                     / ( getKnotValue(j+order-1) - getKnotValue(j) );

        //newCtrlPlygn[j] = alpha*rg_BSplineCurve3D::getCtrlPt(j)
        //                  + (1-alpha)*rg_BSplineCurve3D::getCtrlPt(j-1);

        newCtrlPlygn[j] = alpha*getWeight(j)*rg_BSplineCurve3D::getCtrlPt(j)
                          + (1-alpha)*getWeight(j-1)*rg_BSplineCurve3D::getCtrlPt(j-1);//test

        newWeight_vector[j] = alpha*getWeight(j)
                          + (1-alpha)*getWeight(j-1);

		newCtrlPlygn[j] = newCtrlPlygn[j] / newWeight_vector[j]; //test
    }
    for (j=k+1; j<n+1; j++)
    {
        newCtrlPlygn[j] = rg_BSplineCurve3D::getCtrlPt(j-1);
		newWeight_vector[j] = getWeight(j-1);
    }


    rg_BSplineCurve3D::setCtrlPts(n+1, newCtrlPlygn);

	rg_NUBSplineCurve3D::setKnotVector(n+order+1, newKnotVector);

	setWeightVector(newWeight_vector);

}
void rg_NURBSplineCurve3D::removeRedundantKnot(const rg_REAL& knot)
{
	const rg_INT indexOfKnot=getIndexOfKnotSpan(knot);
	const rg_INT multiplicityOfKnot=getKnotMultiplicity(knot);

	if ( rg_LT(indexOfKnot , 0.0) )
	{
		return;
	}

	rg_NURBSplineCurve3D curve(order);
	curve.setNumOfCtrlPts(numOfCtrlPts-1);
	curve.setInitialKnotVector();

	const rg_INT secondStart=indexOfKnot-order+1;
	const rg_INT thirdStart=indexOfKnot-multiplicityOfKnot+1;

	rg_INT i = 0;
	for( i=0; i < secondStart; i++ )
	{
		rg_Point3D oldPt=getCtrlPt(i);
		rg_REAL    oldWeight=getWeight(i);
		curve.setCtrlPt(i, oldPt);
		curve.setWeight(i, oldWeight);
	}

	for( i=secondStart; i < thirdStart; i++ )
	{
		rg_Point3D oldHpt=getCtrlPt(i);
		rg_REAL    oldWeight=getWeight(i);
		oldHpt=oldWeight*oldHpt;

		rg_Point3D   prevHpt=curve.getCtrlPt(i-1);
		rg_REAL    prevWeight=curve.getWeight(i-1);
		prevHpt=prevWeight*prevHpt;

		const rg_REAL ti=getKnotValue(i);
		const rg_REAL ti_=getKnotValue(i+order);
		const rg_REAL ratio=(knot-ti)/(ti_-ti);

		rg_Point3D hPt=(oldHpt-(1-ratio)*prevHpt)/ratio;
		rg_REAL    weight=(oldWeight-(1-ratio)*prevWeight)/ratio;
		rg_Point3D pt=hPt/weight;

		curve.setCtrlPt(i, pt);
		curve.setWeight(i, weight);
	}

	for( i=thirdStart; i < numOfCtrlPts-1; i++ )
	{
		rg_Point3D oldPt=getCtrlPt(i+1);
		rg_REAL    oldWeight=getWeight(i+1);
		curve.setCtrlPt(i, oldPt);
		curve.setWeight(i, oldWeight);
	}

	const rg_INT numOfKnots=this->getNumberOfKnotValues();
	for( i=0 ; i < indexOfKnot; i++ )
	{
		const rg_REAL knot=getKnotValue(i);
		curve.setKnotValue(i,knot);
	}

	for( i=indexOfKnot+1 ; i <numOfKnots ; i++ )
	{
		const rg_REAL knot=getKnotValue(i);
		curve.setKnotValue(i-1,knot);
	}
	
	*this=curve;
}

rg_RBzCurve3D* rg_NURBSplineCurve3D::decomposeCurveIntoBezierSegment() const
{

    //////////////////////////////////////////////////////////////////////////////////////////
	rg_NURBSplineCurve3D duplicatedNURBSpline = (*this);

	//duplicatedNURBSpline.reparameterizationKnotVector();

	// record the index, the knot value 
	// and multiplicity whose multiplicity is more than zero
	// in the "multiplicity[][]"
	//////////////////////////////////////////////////////

	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT p = getOrder() - 1;
	rg_INT m = n + p + 1;

	rg_INT numberOfKnotSpan     = getNumOfNonZeroLengthKnotSpan(); // the number of nonzero length knot span
    rg_INT numberOfInteriorKnot = numberOfKnotSpan-1;
	
    rg_REAL** multiplicity = new rg_REAL* [numberOfInteriorKnot];

	rg_INT i = 0;
	for(i = 0;i < numberOfInteriorKnot;i++)
		multiplicity[ i ] = new rg_REAL[ 2 ];

    i     = 0;      //  i traces the array multiplicity.
    rg_INT j = p+1;    //  j traces knotVector.
	while( i < numberOfInteriorKnot )
	{
		rg_INT s = 0;
		while( (j < m - p) && rg_EQ(getKnotValue(j), getKnotValue(j+1)))
		{
			s++;
			j++;
		}
		multiplicity[i][ 0 ] = getKnotValue(j);  // replicated knot value in the interior knot span
		multiplicity[i][ 1 ] = s;                // multiplicity of replicated knot value in the interior knot span
		i++;
        j++;
	}
	//////////////////////////////////////////////////////

	
	// knot insertion is applied into the interior knot whose multiplicity is less than p(degree)
	//////////////////////////////////////////////////////

    for (i=0; i<numberOfInteriorKnot; i++)
	{
		rg_INT j = 0;
		// knot insertion is repeated until the total multiplicity of any one interior knot is p!!
		while(j < (p - multiplicity[i][1] - 1))
		{
			duplicatedNURBSpline.knotInsertion(multiplicity[i][0]);
			j++;
		}
	}

	//////////////////////////////////////////////////////

	//////////////////////////////////////////////////////
	// provide the storage for Bezier segments
	//////////////////////////////////////////////////////

    rg_RBzCurve3D* RBezierCurveList = new rg_RBzCurve3D[numberOfKnotSpan];

	//////////////////////////////////////////////////////
	// assign the respective newly-generated control points to each Bezier curve segment.
	//////////////////////////////////////////////////////

    for (i=0; i<numberOfKnotSpan; i++)
    {
		RBezierCurveList[ i ].setDegree( p );
		for(rg_INT j = 0;j < p + 1;j++)
		{
			RBezierCurveList[ i ].setCtrlPt(j, duplicatedNURBSpline.getCtrlPt(i * p + j));
			RBezierCurveList[ i ].setWeight(j, duplicatedNURBSpline.getWeight(i * p + j));
		}
	}


	//////////////////////////////////////////////////////

	for(i = 0;i < numberOfInteriorKnot;i++)
		delete [] multiplicity[ i ];

	delete[] multiplicity;

	return RBezierCurveList;
}
rg_NURBSplineCurve3D rg_NURBSplineCurve3D::decomposeCurveIntoBezierSegmentInNURBS() const
{
	rg_NURBSplineCurve3D duplicatedNURBSpline = (*this);

	//duplicatedNURBSpline.reparameterizationKnotVector();

	// record the index, the knot value 
	// and multiplicity whose multiplicity is more than zero
	// in the "multiplicity[][]"
	//////////////////////////////////////////////////////

	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT p = getOrder() - 1;
	rg_INT m = n + p + 1;

	rg_INT numberOfKnotSpan     = getNumOfNonZeroLengthKnotSpan(); // the number of nonzero length knot span
    rg_INT numberOfInteriorKnot = numberOfKnotSpan-1;
	
    rg_REAL** multiplicity = new rg_REAL* [numberOfInteriorKnot];

	rg_INT i = 0;
	for(i = 0;i < numberOfInteriorKnot;i++)
		multiplicity[ i ] = new rg_REAL[ 2 ];

    i     = 0;      //  i traces the array multiplicity.
    rg_INT j = p+1;    //  j traces knotVector.
	while( i < numberOfInteriorKnot )
	{
		rg_INT s = 0;
		while( (j < m - p) && rg_EQ(getKnotValue(j), getKnotValue(j+1)))
		{
			s++;
			j++;
		}
		multiplicity[i][ 0 ] = getKnotValue(j);  // replicated knot value in the interior knot span
		multiplicity[i][ 1 ] = s;                // multiplicity of replicated knot value in the interior knot span
		i++;
        j++;
	}
	//////////////////////////////////////////////////////

	
	// knot insertion is applied into the interior knot whose multiplicity is less than p(degree)
	//////////////////////////////////////////////////////

    for (i=0; i<numberOfInteriorKnot; i++)
	{
		rg_INT j = 0;
		// knot insertion is repeated until the total multiplicity of any one interior knot is p!!
		while(j < (p - multiplicity[i][1] - 1))
		{
			duplicatedNURBSpline.knotInsertion(multiplicity[i][0]);
			j++;
		}
	}

	return duplicatedNURBSpline;

}
////    Operator Overloading.
////    Functions to obtain the information of KNOT

rg_NURBSplineCurve3D rg_NURBSplineCurve3D::evaluateCurveSegment( const rg_REAL& start,
                                                      const rg_REAL& end ) const
{
    rg_INT order=rg_BSplineCurve3D::getOrder();
    rg_INT numOfKnots=rg_NUBSplineCurve3D::getNumberOfKnotValues();

    rg_REAL startParam=rg_NUBSplineCurve3D::getKnotValue(0);
    rg_REAL endParam=rg_NUBSplineCurve3D::getKnotValue(numOfKnots-1);

    if (   !rg_BTOR(startParam,start,endParam)
        || !rg_BTOR(startParam,end,endParam)
        || !rg_LT( start, end ) )
    {
        return rg_NURBSplineCurve3D();
    }


    // Insert knot whose mutiplicity is
    //       order      if parameter is one of end parameters
    //       order-1    if parameter is in (start,end)
    
    rg_NURBSplineCurve3D curve(*this);
    rg_INT multiplicity=0;
    rg_INT insertingTimes=0;

    // insert knot at start
    multiplicity=rg_NUBSplineCurve3D::getKnotMultiplicity(start);
    insertingTimes=order-1-multiplicity;
    if( rg_EQ( start, startParam ) )
    {
        insertingTimes++;
    }
    curve.knotInsertion(start, insertingTimes);

    // insert knot at end    
    multiplicity=rg_NUBSplineCurve3D::getKnotMultiplicity(end);
    insertingTimes=order-1-multiplicity;
    if( rg_EQ( end, endParam ) )
    {
        insertingTimes++;
    }
    curve.knotInsertion(end, insertingTimes);
    
    // compute knot vector
    rg_INT newOrder=rg_BSplineCurve3D::getOrder();

    rg_INT startKnotIndex=curve.rg_NUBSplineCurve3D::getIndexOfKnotSpan(start);
    rg_INT endKnotIndex=curve.rg_NUBSplineCurve3D::getIndexOfKnot(end);
    rg_INT newNumOfKnots=endKnotIndex-startKnotIndex+2*newOrder-1;
    rg_REAL* newKnots=new rg_REAL[newNumOfKnots];

    rg_INT knotIndex=0;
    rg_INT i = 0;
	for( i=0; i < newOrder-1; i++ )
    {
        newKnots[knotIndex]=0.0;
        knotIndex++;
    }

    rg_REAL startKnot=curve.getKnotValue(startKnotIndex);
    rg_REAL endKnot=curve.getKnotValue(endKnotIndex);
    for( i=0; i < endKnotIndex-startKnotIndex+1; i++ )
    {
        rg_REAL currentKnot=curve.getKnotValue(startKnotIndex+i);
        newKnots[knotIndex]=(currentKnot-startKnot)/(endKnot-startKnot);
        knotIndex++;
    }
    for( i=0; i < newOrder-1; i++ )
    {
        newKnots[knotIndex]=1.0;
        knotIndex++;
    }

    // compute control points and weights
    rg_INT startCtrlPtIndex=startKnotIndex-newOrder+1;
    rg_INT newNumOfCtrlPts=newNumOfKnots-newOrder;
    rg_Point3D* newCtrlPts=new rg_Point3D[newNumOfCtrlPts];
    rg_REAL*    newWeights=new rg_REAL[newNumOfCtrlPts];

    for( i=0; i < newNumOfCtrlPts; i++ )
    {
        newCtrlPts[i]=curve.getCtrlPt(startCtrlPtIndex+i);
        newWeights[i]=curve.getWeight(startCtrlPtIndex+i);
    }

    rg_NURBSplineCurve3D output;
    output.setOrder(newOrder);
    output.setCtrlPts(newNumOfCtrlPts,newCtrlPts);
    output.setKnotVector(newNumOfKnots,newKnots);
    output.setWeightVector(newWeights);

    delete[] newKnots;
    delete[] newCtrlPts;
    delete[] newWeights;

    return output;
    
}

void rg_NURBSplineCurve3D::knotInsertion( const rg_REAL& insertingKnot,
                                       const rg_INT&  insertingTimes)
{
    for( rg_INT i =0 ; i < insertingTimes; i++)
    {
        rg_NURBSplineCurve3D::knotInsertion(insertingKnot);
    }
}

//  April 7 1997 : made.

rg_NURBSplineCurve3D& rg_NURBSplineCurve3D::operator =(const rg_NURBSplineCurve3D &curve)
{
    if (this == &curve)
        return *this;

    rg_INT n     = curve.rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT order = curve.rg_BSplineCurve3D::getOrder();

    rg_Curve::setID(curve.getID());
    rg_Curve::setPlanarity(curve.getPlanarity());

    rg_BSplineCurve3D::setOrder( order );

    rg_BSplineCurve3D::setCtrlPts( n,curve.rg_BSplineCurve3D::getCtrlPts() );

    rg_NUBSplineCurve3D::setKnotVector( n+order, 
                                     curve.rg_NUBSplineCurve3D::getKnotVector() );

    if ( weights !=rg_NULL )
        delete [] weights;

    weights = new rg_REAL [n];

    for (rg_INT i=0; i<n; i++)
        weights[i] = curve.weights[i];

    return *this;
}

void  rg_NURBSplineCurve3D::removeAll()
{
    if ( weights != rg_NULL )
    {
        delete [] weights;
        weights=rg_NULL;
    }

    rg_NUBSplineCurve3D::removeAll();
}
    



