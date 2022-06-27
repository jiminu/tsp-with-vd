//******************************************************************************
//
//	  FILENAME    : rg_NUBSplineCurve3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_NUBSplineCurve3D 
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//
//    START DATE  : 21 Jun 1996    
//
//    HISTORY     :
//          BY Young-Song Cho.  13 Jul. 1997
//            make : rg_INT        getNumOfKnotSpan() 
//            make : rg_BzCurve3D* separateBSplineIntoBezier()   
//            make : rg_sListByPtr*     intersectOfCubicBSplineAndPlane(const rg_Plane3D &plane)
//
//          BY Young-Song Cho.  24 Jul. 1997
//            make : rg_BzCurve3D* pullOutCubicBezierForKnotSpan(
//                                  const rg_REAL &prevKnot, 
//                                  const rg_REAL &nextKnot ) 
//            make : rg_BzCurve3D* pullOutCubicBezierForKnotSpan(
//                                  const rg_INDEX &indexOfKnotSpan) 
//
//          BY Young-Song Cho.  14 Aug. 1997
//            make : rg_REAL* evaluateMultiBasisFunc( const rg_INDEX     &knotIndex,
//                                                 const rg_PARAMETER &u,
//                                                 const rg_ORDER     &order)
//            modify : rg_Point3D evaluatePt( const rg_REAL &param )
//
//          By Dong-Gyou Lee 18 Mar. 1998
//                	void powerSplineToNUBSplineCurve( const rg_DEGREE& dgr,
//                                                    const rg_INT& numOfSegment,
//                                                    rg_REAL** paramValues,
//                                                    const rg_Matrix& powerCoeff )
//
//                  void bzCurvesToNUBSplineCurve( const rg_INT& numOfSegment,
//                                                 const rg_BzCurve3D* bzCurves )
//
//          By Dong-Gyou Lee 26 Mar. 1998
//                  void reparameterizationKnotVector() 
//
//			By Jung-Hyun Ryu Apr.11 1998
//					void multipleKnotInsertion();
//
//           By Young-Song Cho Jul.7 1998 
//                  rg_FLAG makeCompositeCurveWithC0(const rg_NUBSplineCurve3D &curve1,
//                                                const rg_NUBSplineCurve3D &curve2);
//           By Young-Song Cho Jul.7 1998 
//                  void formLine( const rg_Point3D& start, const rg_Point3D& end, const rg_INT &order = 4 );
//
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering   
//                          Hanyang University, Seoul Korea	  	
//
//******************************************************************************


#include <math.h>

#include "rg_BandedMatrix.h"
#include "rg_SignArray.h"
#include "rg_RelativeOp.h"
#include "rg_CurveSurfaceFunc.h"

#include "rg_ComplexNumber.h"
#include "rg_TMatrix3D.h"
#include "rg_BzCurve2D.h"
#include "rg_dList.h"
#include "rg_Matrix.h"
#include "rg_NUBSplineCurve3D.h"
#include "rg_GeoFunc.h"
#include "rg_MathFunc.h"
#include "rg_NURBSCurveIntersector.h"
// for test.
//#include <fstream>

#include <time.h>

// private member functions
rg_Point3D*   rg_NUBSplineCurve3D::makeFilteredPassingPts(       rg_INT&     numOfPts,
		                                                   const rg_Point3D* const pts )
{
	rg_Point3D* newPts= new rg_Point3D[numOfPts];
	rg_INT newNumOfPts=1;
	rg_Point3D prevPt=pts[0];
	newPts[0]=prevPt;
	
	rg_INT i = 0;
	for( i=1; i < numOfPts; i++ )
	{
		if ( prevPt == pts[i] )
		{
			continue;
		}
		else
		{
			newPts[newNumOfPts]=pts[i];
			newNumOfPts++;
		}
	}
	numOfPts=newNumOfPts;
	rg_Point3D* output=new rg_Point3D[newNumOfPts];
	for( i=0; i < newNumOfPts; i++ )
	{
		output[i]=newPts[i];
	}
	delete[] newPts;

	return output;
}






////	Constructor & Destructor-----------------------------------------------
rg_NUBSplineCurve3D::rg_NUBSplineCurve3D()
: rg_BSplineCurve3D()
{
    knotVector = rg_NULL;
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                                    const rg_Planarity    &newPlanarity )
: rg_BSplineCurve3D(newID, newPlanarity)
{
    knotVector = rg_NULL;
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                                    const rg_Planarity    &newPlanarity,
                                    const rg_ORDER        &newOrder )
: rg_BSplineCurve3D(newID, newPlanarity, newOrder)
{
    knotVector = rg_NULL;
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                                    const rg_Planarity    &newPlanarity,
                                    const rg_INT          &num )
: rg_BSplineCurve3D(newID, newPlanarity, num)
{
    knotVector = rg_NULL;
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                                    const rg_Planarity    &newPlanarity,
                                    const rg_INT          &num, 
                                    rg_Point3D*              newControlP)
: rg_BSplineCurve3D(newID, newPlanarity, num, newControlP)
{
    knotVector = rg_NULL;
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const rg_ORDER &newOrder )
: rg_BSplineCurve3D(newOrder)
{
    knotVector = rg_NULL;
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const rg_INT &num )
: rg_BSplineCurve3D(num)
{
    knotVector = rg_NULL;
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const rg_INT &num, 
                                    rg_Point3D*     newControlP )
: rg_BSplineCurve3D(num, newControlP)
{
    knotVector = rg_NULL;
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const rg_ORDER &newOrder, 
                                    const rg_INT   &num, 
                                    rg_Point3D*       newControlP )
: rg_BSplineCurve3D(newOrder, num, newControlP)
{
    knotVector = rg_NULL;
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const rg_ORDER &newOrder, 
                                    const rg_INT   &num, 
                                    rg_Point3D*       newControlP,
                                    rg_REAL*        newKnotVector )
: rg_BSplineCurve3D(newOrder, num, newControlP)
{
    knotVector = new rg_REAL[ num + newOrder ];

    for (rg_INDEX i=0; i<(num+newOrder); i++)
        knotVector[i] = newKnotVector[i];
}

////  Constructor      : March 13 1997
rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                                    const rg_Planarity    &newPlanarity,
                                    const rg_ORDER        &newOrder, 
                                    const rg_INT          &num, 
                                    rg_Point3D*              newControlP,
                                    rg_REAL*               newKnotVector )
: rg_BSplineCurve3D(newID,
                 newPlanarity,
                 newOrder,
                 num,
                 newControlP)
{
    knotVector = new rg_REAL [newOrder+num];

    for (rg_INT i=0; i<(rg_INT)(newOrder+num); i++)
        knotVector[i] = newKnotVector[i];
}

rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const rg_BSplineCurve3D &curve)
: rg_BSplineCurve3D( curve )
{
    knotVector = rg_NULL;
}

////  Copy Constructor March 13 1997
rg_NUBSplineCurve3D::rg_NUBSplineCurve3D( const rg_NUBSplineCurve3D &curve)
: rg_BSplineCurve3D(curve)
{
	if ( curve.isNull() )
	{
		knotVector=NULL;
		return;
	}

    rg_INT n = rg_BSplineCurve3D::getOrder() + rg_BSplineCurve3D::getNumOfCtrlPts();

    knotVector = new rg_REAL [n];

    for (rg_INT i=0; i<n; i++)
        knotVector[i] = curve.knotVector[i];
}

rg_NUBSplineCurve3D::~rg_NUBSplineCurve3D()
{
    if ( knotVector != rg_NULL )
		delete [] knotVector;
}

////	Get Functions.---------------------------------------------------------
rg_REAL rg_NUBSplineCurve3D::getKnotValue( const rg_INT &kIndex ) const
{
    return knotVector[kIndex];
}

rg_REAL  rg_NUBSplineCurve3D::getStartParameter() const
{
    return knotVector[0];
}

rg_REAL  rg_NUBSplineCurve3D::getEndParameter() const
{
    rg_INT order=rg_BSplineCurve3D::getOrder();
    rg_INT numOfCtrlPts=rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT endIndex=order+numOfCtrlPts-1;

    return knotVector[endIndex];
}

rg_Point3D  rg_NUBSplineCurve3D::getStartPoint() const
{
    const rg_REAL startParameter=rg_NUBSplineCurve3D::getStartParameter();
    rg_Point3D output=evaluatePt(startParameter);

    return output;
}

rg_Point3D  rg_NUBSplineCurve3D::getEndPoint() const
{
    const rg_REAL endParameter=rg_NUBSplineCurve3D::getEndParameter();
    rg_Point3D output=evaluatePt(endParameter);

    return output;
}

rg_REAL* rg_NUBSplineCurve3D::getKnotVector() const
{
    rg_INT numOfCtrlPts     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT order            = rg_BSplineCurve3D::getOrder();
    rg_INT numOfKnotValue   = numOfCtrlPts + order;
    rg_REAL* output= new rg_REAL[numOfKnotValue];

    for ( rg_INT i=0 ; i < numOfKnotValue; i++ )
    {
        output[i]=knotVector[i];
    }
    return output;
}

// Assume that the multiplicity of all knots is 1.
rg_INT rg_NUBSplineCurve3D::getNumOfKnotSpan() const
{
    rg_INT n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT order = (rg_INT) rg_BSplineCurve3D::getOrder();

    return n + order - 2*( order -1 ) - 1;
}

// Regardless of the multiplicity of all knots(1 or more than 1)
// find the No. of NonZeroLengthKnotSpan
rg_INT rg_NUBSplineCurve3D::getNumOfNonZeroLengthKnotSpan() const
{
	rg_INT NumOfNonZeroLengthKnotSpan = 0;
	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT p = getOrder() - 1;
	rg_INT i = p;
	rg_INT m = n + p + 1;

	while(i <= m - p - 1)
	{
		if( rg_NE(knotVector[i], knotVector[i+1]) )
			NumOfNonZeroLengthKnotSpan++;
		i++;
	}
	return NumOfNonZeroLengthKnotSpan;
}



rg_REAL*  rg_NUBSplineCurve3D::getDistinctKnotValues() const
{
	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT p = getOrder() - 1;
	rg_INT i = p;
	rg_INT m = n + p + 1;

	rg_INT NumOfNonZeroLengthKnotSpan = getNumOfNonZeroLengthKnotSpan();
	rg_REAL* distinctKnotValue = new rg_REAL[NumOfNonZeroLengthKnotSpan + 1];

	rg_INT index = 0;
	while((index < NumOfNonZeroLengthKnotSpan + 1) && (i <= m))
	{
		if(rg_NE(knotVector[ i ], knotVector[i + 1]))
		{
			distinctKnotValue[index] = knotVector[ i ];
			index++;
		}
		i++;
	}

	return distinctKnotValue;
}

rg_INT* rg_NUBSplineCurve3D::getIndexOfNonZeroBasisInKnotSpan(const rg_INDEX& indexOfKnotValue) const
{
	rg_DEGREE p = getOrder() - 1;

	rg_INT* ctrlPointIndexOfBasisThatIsNonZeroInGivenKnotSpan = new rg_INT[p + 1 + 1];

	rg_INT indexOfCtrlPoint = 0;
	rg_INT indexOfKnotVector = indexOfCtrlPoint;
	rg_INT index = 0;

	// Write the index of control point for nonzero basis N_i_p in given knot value

	while(rg_LE(knotVector[indexOfKnotVector], knotVector[indexOfKnotValue]) && (index < p + 1))
	{
		while(indexOfKnotVector <= indexOfCtrlPoint + p + 1)
		{
			if((rg_EQ(knotVector[indexOfKnotVector], knotVector[indexOfKnotValue])) && (rg_NE(knotVector[indexOfKnotVector], knotVector[indexOfKnotVector + 1])))
			{
				ctrlPointIndexOfBasisThatIsNonZeroInGivenKnotSpan[index] = indexOfCtrlPoint;
				index++;
			}
			indexOfKnotVector++;
		}

		indexOfCtrlPoint++;
		indexOfKnotVector = indexOfCtrlPoint;
	}

	// Write the corresponding no. of nonzero basis in (p + 2)th element

	ctrlPointIndexOfBasisThatIsNonZeroInGivenKnotSpan[p + 1] = index;

	return ctrlPointIndexOfBasisThatIsNonZeroInGivenKnotSpan;
}


// Given i(indexOfCtrlPoint), p(degree) and i(indexOfKnotValue), 
// construct the polynomial form of N_i_j(u).

rg_Polynomial rg_NUBSplineCurve3D::getBasisWithPolynomialFormIn(const rg_INDEX& indexOfCtrlPoint, // index for control point(index of a basis function)
														  const rg_DEGREE& degree, // degree
														  const rg_INDEX& indexOfKnotValue /*index for knot value*/) const
{
	if(degree == 0)
	{
		rg_Polynomial basis( 0 );

		if(indexOfCtrlPoint == indexOfKnotValue)
		{
			basis.setCoefficient(0, 1.0);
			return basis;
		}
		else
		{
			basis.setCoefficient(0, 0.0);
			return basis;
		}
	}
	else
	{
		rg_Polynomial firstMultiplier;
		rg_Polynomial secondMultiplier;

		if(rg_EQ(knotVector[indexOfCtrlPoint + degree], knotVector[indexOfCtrlPoint]))
		{
			firstMultiplier.setDegree( 0 );
			firstMultiplier.setCoefficient(0, 0.0);
		}
		else
		{
			firstMultiplier.setDegree( 1 );
			firstMultiplier.setCoefficient(0, - knotVector[indexOfCtrlPoint] / (knotVector[indexOfCtrlPoint + degree] - knotVector[indexOfCtrlPoint]));
			firstMultiplier.setCoefficient(1, 1.0 / (knotVector[indexOfCtrlPoint + degree] - knotVector[indexOfCtrlPoint]));
		}

		if(rg_EQ(knotVector[indexOfCtrlPoint + degree + 1], knotVector[indexOfCtrlPoint + 1]))
		{
			secondMultiplier.setDegree( 0 );
			secondMultiplier.setCoefficient(0, 0.0);
		}
		else
		{
			secondMultiplier.setDegree( 1 );
			secondMultiplier.setCoefficient(0, knotVector[indexOfCtrlPoint + degree + 1] / (knotVector[indexOfCtrlPoint + degree + 1] - knotVector[indexOfCtrlPoint + 1]));
			secondMultiplier.setCoefficient(1, - 1.0 / (knotVector[indexOfCtrlPoint + degree + 1] - knotVector[indexOfCtrlPoint + 1]));
		}

		return ((firstMultiplier * getBasisWithPolynomialFormIn(indexOfCtrlPoint, degree - 1, indexOfKnotValue)) 
			  + (secondMultiplier * getBasisWithPolynomialFormIn(indexOfCtrlPoint + 1, degree - 1, indexOfKnotValue)));
	}
}

rg_PolynomialWithBound* rg_NUBSplineCurve3D::getBasisWithPolynomialFormIn(const rg_INDEX& indexOfCtrlPoint) const
{
	rg_INT NumOfNonZeroLengthKnotSpan = 0;
	rg_INT n = getNumOfCtrlPts() - 1;
	rg_DEGREE p = getOrder() - 1;
	rg_INT m = n + p + 1;

	rg_INDEX i = indexOfCtrlPoint; 

	while(i <= indexOfCtrlPoint + p + 1)
	{
		if( rg_NE(knotVector[i], knotVector[i+1]) )
			NumOfNonZeroLengthKnotSpan++;
		i++;
	}

	rg_PolynomialWithBound* basis = new rg_PolynomialWithBound [NumOfNonZeroLengthKnotSpan];

	i = indexOfCtrlPoint; // initialization of index for touring the index of knot value
	//rg_INT j = 0; // index of degree
	rg_INT index = 0; // index of rg_Polynomial List

	rg_PolynomialWithBound onePolynomialFormBasis;

	while((index < NumOfNonZeroLengthKnotSpan) && (i < indexOfCtrlPoint + p + 1))
	{
		if(rg_NE(knotVector[ i ], knotVector[i + 1]))
		{
			onePolynomialFormBasis = getBasisWithPolynomialFormIn(indexOfCtrlPoint, p, i);
			onePolynomialFormBasis.setBound(knotVector[i], knotVector[i + 1]);
			basis[index] = onePolynomialFormBasis;
			index++;
		}
		// One polynomial is consructed.
		i++;
	}
	// One rg_Polynomial list is constructed.

	return basis;
}


rg_PolynomialWithBound** rg_NUBSplineCurve3D::makePolynomialFormCurve() const
{
	rg_INT NumOfNonZeroLengthKnotSpan = getNumOfNonZeroLengthKnotSpan();
	rg_INT n = getNumOfCtrlPts() - 1;
	rg_DEGREE p = getOrder() - 1;
	rg_INT m = n + p + 1;

	rg_PolynomialWithBound** polynomialFormCurve = new rg_PolynomialWithBound* [ 3 ];

	for(rg_INT i = 0;i < 3;i++)
		polynomialFormCurve[ i ] = new rg_PolynomialWithBound[NumOfNonZeroLengthKnotSpan];

	rg_INT indexOfKnotVector = p;
	rg_INT index = 0;

	while((indexOfKnotVector <= m - p - 1) && (index < NumOfNonZeroLengthKnotSpan))
	{
		if(rg_NE(knotVector[indexOfKnotVector], knotVector[indexOfKnotVector + 1]))
		{
			rg_INT* indexList = getIndexOfNonZeroBasisInKnotSpan(indexOfKnotVector);

			rg_PolynomialWithBound tempX( p ), tempY( p ), tempZ( p );
			rg_INT i = 0;

			// indexList[p + 1]에 현재의 knot span에 영향을 주는 basis function의 개수를 저당
			// 최대로 p + 1개가 될 수 있다.

			while(i < indexList[p + 1])
			{
				tempX = tempX + getCtrlPt(indexList[ i ]).getX() * getBasisWithPolynomialFormIn(indexList[ i ], p, indexOfKnotVector);
				tempY = tempY + getCtrlPt(indexList[ i ]).getY() * getBasisWithPolynomialFormIn(indexList[ i ], p, indexOfKnotVector);
				tempZ = tempZ + getCtrlPt(indexList[ i ]).getZ() * getBasisWithPolynomialFormIn(indexList[ i ], p, indexOfKnotVector);
				// PolynomialwithBound간의 operation을 재정의 할 경우 수정이 필요할 지도 모른다?
				i++;
			}

			tempX.setBound(knotVector[indexOfKnotVector], knotVector[indexOfKnotVector + 1]);
			tempY.setBound(knotVector[indexOfKnotVector], knotVector[indexOfKnotVector + 1]);
			tempZ.setBound(knotVector[indexOfKnotVector], knotVector[indexOfKnotVector + 1]);

			delete[] indexList;

			polynomialFormCurve[ 0 ][index] = tempX;
			polynomialFormCurve[ 1 ][index] = tempY;
			polynomialFormCurve[ 2 ][index] = tempZ;

			index++;
		}
		indexOfKnotVector++;
	}

	return polynomialFormCurve;
}


rg_PolynomialWithBound** rg_NUBSplineCurve3D::makePolynomialFormCurveUsingCurveDecomp() const
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

	//rg_BzCurve3D* BezierCurveList = decomposeCurveIntoBezierSegment();
	rg_BzCurve3D* BezierCurveList = decomposeCurveIntoBezierSegmentUsingKnotRefinement();

	rg_PolynomialWithBound** polynomialFormCurve = new rg_PolynomialWithBound* [ 3 ];

	rg_INT i = 0;
	for(i = 0;i < 3;i++)
		polynomialFormCurve[ i ] = new rg_PolynomialWithBound[NumOfNonZeroLengthKnotSpan];

	rg_REAL* distinctKnotValue = getDistinctKnotValues();
	rg_Matrix P(p + 1, 3);
	rg_Matrix M_p = rg_CurveSurfaceFunc::bezierToPowerMatrix(p + 1);

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
		                               1/(distinctKnotValue[i + 1] - distinctKnotValue[ i ]),
			                           - distinctKnotValue[ i ] /(distinctKnotValue[i + 1] - distinctKnotValue[ i ]));
		rg_INT j = 0;
		for(j = 0;j < p + 1;j++)
		{
			P[ j ][ 0 ] = (BezierCurveList[ i ].getCtrlPt( j )).getX(); // X coordinate of ctrl pt
			P[ j ][ 1 ] = (BezierCurveList[ i ].getCtrlPt( j )).getY(); // Y coordinate of ctrl pt
			P[ j ][ 2 ] = (BezierCurveList[ i ].getCtrlPt( j )).getZ(); // Z coordinate of ctrl pt
		}

		rg_Matrix coefficient = R_p * M_p * P; // Coefficient of power basis
										    // (p + 1) by 3 matrix
		
		rg_PolynomialWithBound* temp = new rg_PolynomialWithBound[ 3 ];

		for(j = 0;j < 3;j++)
		{
			temp[ j ].setDegree( p );
			for(rg_INT k = 0;k < p + 1;k++)
				temp[ j ].setCoefficient(k, coefficient[ k ][ j ]);
			temp[ j ].setBound(distinctKnotValue[ i ], distinctKnotValue[i + 1]);
			polynomialFormCurve[ j ][ i ] = temp[ j ];			
		}

		delete[] temp;
	}

	delete[] BezierCurveList;

	return polynomialFormCurve;
}


rg_PolynomialWithBound** rg_NUBSplineCurve3D::makePolynomialFormCurveUsingCurveDecomp(rg_REAL& timeForKnotRefinement1,
	                                                                            rg_REAL& timeForKnotRefinement2) const
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

	// checking for time1======================
	clock_t StartTime, EndTime;
    StartTime = clock();

	rg_INT NumOfNonZeroLengthKnotSpan = getNumOfNonZeroLengthKnotSpan();
	rg_DEGREE p = getOrder() - 1;

	//rg_BzCurve3D* BezierCurveList = decomposeCurveIntoBezierSegment();
	rg_BzCurve3D* BezierCurveList = decomposeCurveIntoBezierSegmentUsingKnotRefinement();

	EndTime = clock();
	timeForKnotRefinement1 += (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	//========================================

	// checking for time2======================
	StartTime = clock();

	rg_PolynomialWithBound** polynomialFormCurve = new rg_PolynomialWithBound* [ 3 ];

	rg_INT i = 0;
	for(i = 0;i < 3;i++)
		polynomialFormCurve[ i ] = new rg_PolynomialWithBound[NumOfNonZeroLengthKnotSpan];

	rg_REAL* distinctKnotValue = getDistinctKnotValues();
	rg_Matrix P(p + 1, 3);
	rg_Matrix M_p = rg_CurveSurfaceFunc::bezierToPowerMatrix(p + 1);

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
		                               1/(distinctKnotValue[i + 1] - distinctKnotValue[ i ]),
			                           - distinctKnotValue[ i ] /(distinctKnotValue[i + 1] - distinctKnotValue[ i ]));
		rg_INT j = 0;
		for(j = 0;j < p + 1;j++)
		{
			P[ j ][ 0 ] = (BezierCurveList[ i ].getCtrlPt( j )).getX(); // X coordinate of ctrl pt
			P[ j ][ 1 ] = (BezierCurveList[ i ].getCtrlPt( j )).getY(); // Y coordinate of ctrl pt
			P[ j ][ 2 ] = (BezierCurveList[ i ].getCtrlPt( j )).getZ(); // Z coordinate of ctrl pt
		}

		rg_Matrix coefficient = R_p * M_p * P; // Coefficient of power basis
										    // (p + 1) by 3 matrix
		
		rg_PolynomialWithBound* temp = new rg_PolynomialWithBound[ 3 ];

		for(j = 0;j < 3;j++)
		{
			temp[ j ].setDegree( p );
			for(rg_INT k = 0;k < p + 1;k++)
				temp[ j ].setCoefficient(k, coefficient[ k ][ j ]);
			temp[ j ].setBound(distinctKnotValue[ i ], distinctKnotValue[i + 1]);
			polynomialFormCurve[ j ][ i ] = temp[ j ];			
		}

		delete[] temp;
	}

	delete[] BezierCurveList;

	EndTime = clock();
	timeForKnotRefinement2 += (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	//========================================

	return polynomialFormCurve;
}

rg_REAL*  rg_NUBSplineCurve3D::getDistributionPolynomialsInOneGraph(rg_INT contributionKnotSpan,/* rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths) const
{
	rg_INT p = getOrder() - 1;

	rg_REAL* coefficientOfPolynomial = new rg_REAL[p + 1];
	rg_REAL* tempCoeff = new rg_REAL[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
		coefficientOfPolynomial[ i ] = 0.0;

	rg_REAL*** coeffOfLinearTerms = new rg_REAL** [noOfAllPossiblePaths];
	
	// locating two coefficient of each linear polynomial
	rg_INT degree;
	rg_INT index;
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		degree = 0;
		index = contributionKnotSpan;
		coeffOfLinearTerms[ i ] = new rg_REAL* [ p ];
		for(rg_INT j = 0;j < p;j++)
		{
			coeffOfLinearTerms[ i ][ j ] = new rg_REAL [ 2 ];
			degree++;

			if(graph[ i ][ j ] == 0)
			{
				coeffOfLinearTerms[ i ][ j ][ 0 ] = - knotVector[index] / (knotVector[index + degree] - knotVector[index]);
				coeffOfLinearTerms[ i ][ j ][ 1 ] = 1.0 / (knotVector[index + degree] - knotVector[index]);
			}
			else
			{
				index--;
				coeffOfLinearTerms[ i ][ j ][ 0 ] = knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]);
				coeffOfLinearTerms[ i ][ j ][ 1 ] = - 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1]);
			}
		}
	}
	//==================================================

	// Expanding sum of multiplications of linear terms
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		tempCoeff[ 0 ] = coeffOfLinearTerms[ i ][ 0 ][ 0 ];
		tempCoeff[ 1 ] = coeffOfLinearTerms[ i ][ 0 ][ 1 ];
		for(rg_INT j = 1;j <= p - 1;j++)
		{
			tempCoeff[j + 1] = tempCoeff[ j ] * coeffOfLinearTerms[ i ][ j ][ 1 ];
			for(rg_INT k = j;k >= 1;k--)
			{
				tempCoeff[ k ] = tempCoeff[k - 1] * coeffOfLinearTerms[ i ][ j ][ 1 ] + tempCoeff[ k ] * coeffOfLinearTerms[ i ][ j ][ 0 ];
			}
			tempCoeff[ 0 ] *= coeffOfLinearTerms[ i ][ j ][ 0 ];
		}
		for(rg_INT l = 0;l < p + 1;l++)
			coefficientOfPolynomial[ l ] += tempCoeff[ l ]; 
	}
	//==================================================

	delete[] tempCoeff;

	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		for(rg_INT j = 0;j < p;j++)
		{
			delete[] coeffOfLinearTerms[ i ][ j ];
		}
		delete[] coeffOfLinearTerms[ i ];
	}
	delete[] coeffOfLinearTerms;

	return coefficientOfPolynomial;
}

rg_REAL* rg_NUBSplineCurve3D::getDistributionPolynomialsInOneGraphRevised(rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths) const
{
	rg_INT p = getOrder() - 1;

	rg_REAL* coefficientOfPolynomial = new rg_REAL[p + 1];
	rg_REAL* tempCoeff = new rg_REAL[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
		coefficientOfPolynomial[ i ] = 0.0;

	// Expanding sum of multiplications of linear terms according to basic graph
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		rg_INT degree = 0;
		rg_INT index = contributionKnotSpan;

		degree++;

		if(graph[ i ][ 0 ] == 0)
		{
			tempCoeff[ 0 ] = - knotVector[index] / (knotVector[index + degree] - knotVector[index]);
			tempCoeff[ 1 ] = 1.0 / (knotVector[index + degree] - knotVector[index]);
		}
		else
		{
			index--;
			tempCoeff[ 0 ] = knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]);
			tempCoeff[ 1 ] = - 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1]);
		}
		//tempCoeff[ 0 ] = coeffOfLinearTerms[ i ][ 0 ][ 0 ];
		//tempCoeff[ 1 ] = coeffOfLinearTerms[ i ][ 0 ][ 1 ];
		for(rg_INT j = 1;j <= p - 1;j++)
		{
			degree++;

			if(graph[ i ][ j ] == 0)
			{
				tempCoeff[j + 1] = tempCoeff[ j ] * (1.0 / (knotVector[index + degree] - knotVector[index]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempCoeff[ k ] = tempCoeff[k - 1] * (1.0 / (knotVector[index + degree] - knotVector[index])) + tempCoeff[ k ] * (- knotVector[index] / (knotVector[index + degree] - knotVector[index]));
				}
				tempCoeff[ 0 ] *= (- knotVector[index] / (knotVector[index + degree] - knotVector[index]));
			}
			else
			{
				index--;
				tempCoeff[j + 1] = tempCoeff[ j ] * (- 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempCoeff[ k ] = tempCoeff[k - 1] * (- 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1])) + tempCoeff[ k ] * (knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]));
				}
				tempCoeff[ 0 ] *= (knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]));
			}

			//tempCoeff[j + 1] = tempCoeff[ j ] * coeffOfLinearTerms[ i ][ j ][ 1 ];
			//for(rg_INT k = j;k >= 1;k--)
			//{
			//	tempCoeff[ k ] = tempCoeff[k - 1] * coeffOfLinearTerms[ i ][ j ][ 1 ] + tempCoeff[ k ] * coeffOfLinearTerms[ i ][ j ][ 0 ];
			//}
			//tempCoeff[ 0 ] *= coeffOfLinearTerms[ i ][ j ][ 0 ];
		}
		for(rg_INT l = 0;l < p + 1;l++)
			coefficientOfPolynomial[ l ] += tempCoeff[ l ]; 
	}

	return coefficientOfPolynomial;
}


rg_REAL*  rg_NUBSplineCurve3D::getDistributionPolynomialsInOneGraph(rg_INT contributionKnotSpan,/* rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths,
															  rg_REAL& time1, rg_REAL& time2) const
{
	// checking for time1======================
	clock_t StartTime, EndTime;
    StartTime = clock();

	rg_INT p = getOrder() - 1;

	rg_REAL* coefficientOfPolynomial = new rg_REAL[p + 1];
	rg_REAL* tempCoeff = new rg_REAL[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
		coefficientOfPolynomial[ i ] = 0.0;

	rg_REAL*** coeffOfLinearTerms = new rg_REAL** [noOfAllPossiblePaths];
	
	// locating two coefficient of each linear polynomial
	rg_INT degree;
	rg_INT index;
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		degree = 0;
		index = contributionKnotSpan;
		coeffOfLinearTerms[ i ] = new rg_REAL* [ p ];
		for(rg_INT j = 0;j < p;j++)
		{
			coeffOfLinearTerms[ i ][ j ] = new rg_REAL [ 2 ];
			degree++;

			if(graph[ i ][ j ] == 0)
			{
				coeffOfLinearTerms[ i ][ j ][ 0 ] = - knotVector[index] / (knotVector[index + degree] - knotVector[index]);
				coeffOfLinearTerms[ i ][ j ][ 1 ] = 1.0 / (knotVector[index + degree] - knotVector[index]);
			}
			else
			{
				index--;
				coeffOfLinearTerms[ i ][ j ][ 0 ] = knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]);
				coeffOfLinearTerms[ i ][ j ][ 1 ] = - 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1]);
			}
		}
	}
	//==================================================

	EndTime = clock();
	time1 += (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	//========================================

	// checking for time2======================
    StartTime = clock();

	// Expanding sum of multiplications of linear terms
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		tempCoeff[ 0 ] = coeffOfLinearTerms[ i ][ 0 ][ 0 ];
		tempCoeff[ 1 ] = coeffOfLinearTerms[ i ][ 0 ][ 1 ];
		for(rg_INT j = 1;j <= p - 1;j++)
		{
			tempCoeff[j + 1] = tempCoeff[ j ] * coeffOfLinearTerms[ i ][ j ][ 1 ];
			for(rg_INT k = j;k >= 1;k--)
			{
				tempCoeff[ k ] = tempCoeff[k - 1] * coeffOfLinearTerms[ i ][ j ][ 1 ] + tempCoeff[ k ] * coeffOfLinearTerms[ i ][ j ][ 0 ];
			}
			tempCoeff[ 0 ] *= coeffOfLinearTerms[ i ][ j ][ 0 ];
		}
		for(rg_INT l = 0;l < p + 1;l++)
			coefficientOfPolynomial[ l ] += tempCoeff[ l ]; 
	}
	//==================================================

	delete[] tempCoeff;

	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		for(rg_INT j = 0;j < p;j++)
		{
			delete[] coeffOfLinearTerms[ i ][ j ];
		}
		delete[] coeffOfLinearTerms[ i ];
	}
	delete[] coeffOfLinearTerms;

	EndTime = clock();
	time2 += (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	//========================================	

	return coefficientOfPolynomial;
}

rg_REAL* rg_NUBSplineCurve3D::getDistributionPolynomialsInOneGraph1(rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths,
											rg_REAL& time1) const
{
	// checking for time1======================
	clock_t StartTime, EndTime;
    StartTime = clock();

	rg_INT p = getOrder() - 1;

	rg_REAL* coefficientOfPolynomial = new rg_REAL[p + 1];
	rg_REAL* tempCoeff = new rg_REAL[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
		coefficientOfPolynomial[ i ] = 0.0;

	rg_REAL*** coeffOfLinearTerms = new rg_REAL** [noOfAllPossiblePaths];
	
	// locating two coefficient of each linear polynomial
	rg_INT degree;
	rg_INT index;
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		degree = 0;
		index = contributionKnotSpan;
		coeffOfLinearTerms[ i ] = new rg_REAL* [ p ];
		for(rg_INT j = 0;j < p;j++)
		{
			coeffOfLinearTerms[ i ][ j ] = new rg_REAL [ 2 ];
			degree++;

			if(graph[ i ][ j ] == 0)
			{
				coeffOfLinearTerms[ i ][ j ][ 0 ] = - knotVector[index] / (knotVector[index + degree] - knotVector[index]);
				coeffOfLinearTerms[ i ][ j ][ 1 ] = 1.0 / (knotVector[index + degree] - knotVector[index]);
			}
			else
			{
				index--;
				coeffOfLinearTerms[ i ][ j ][ 0 ] = knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]);
				coeffOfLinearTerms[ i ][ j ][ 1 ] = - 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1]);
			}
		}
	}
	//==================================================

	EndTime = clock();
	time1 += (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	//========================================

	// Expanding sum of multiplications of linear terms
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		tempCoeff[ 0 ] = coeffOfLinearTerms[ i ][ 0 ][ 0 ];
		tempCoeff[ 1 ] = coeffOfLinearTerms[ i ][ 0 ][ 1 ];
		for(rg_INT j = 1;j <= p - 1;j++)
		{
			tempCoeff[j + 1] = tempCoeff[ j ] * coeffOfLinearTerms[ i ][ j ][ 1 ];
			for(rg_INT k = j;k >= 1;k--)
			{
				tempCoeff[ k ] = tempCoeff[k - 1] * coeffOfLinearTerms[ i ][ j ][ 1 ] + tempCoeff[ k ] * coeffOfLinearTerms[ i ][ j ][ 0 ];
			}
			tempCoeff[ 0 ] *= coeffOfLinearTerms[ i ][ j ][ 0 ];
		}
		for(rg_INT l = 0;l < p + 1;l++)
			coefficientOfPolynomial[ l ] += tempCoeff[ l ]; 
	}
	//==================================================

	delete[] tempCoeff;

	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		for(rg_INT j = 0;j < p;j++)
		{
			delete[] coeffOfLinearTerms[ i ][ j ];
		}
		delete[] coeffOfLinearTerms[ i ];
	}
	delete[] coeffOfLinearTerms;

	return coefficientOfPolynomial;
}

rg_REAL* rg_NUBSplineCurve3D::getDistributionPolynomialsInOneGraph2(rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths,
                											  rg_REAL& time2) const
{
	clock_t StartTime, EndTime;

	rg_INT p = getOrder() - 1;

	rg_REAL* coefficientOfPolynomial = new rg_REAL[p + 1];
	rg_REAL* tempCoeff = new rg_REAL[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
		coefficientOfPolynomial[ i ] = 0.0;

	rg_REAL*** coeffOfLinearTerms = new rg_REAL** [noOfAllPossiblePaths];
	
	// locating two coefficient of each linear polynomial
	rg_INT degree;
	rg_INT index;
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		degree = 0;
		index = contributionKnotSpan;
		coeffOfLinearTerms[ i ] = new rg_REAL* [ p ];
		for(rg_INT j = 0;j < p;j++)
		{
			coeffOfLinearTerms[ i ][ j ] = new rg_REAL [ 2 ];
			degree++;

			if(graph[ i ][ j ] == 0)
			{
				coeffOfLinearTerms[ i ][ j ][ 0 ] = - knotVector[index] / (knotVector[index + degree] - knotVector[index]);
				coeffOfLinearTerms[ i ][ j ][ 1 ] = 1.0 / (knotVector[index + degree] - knotVector[index]);
			}
			else
			{
				index--;
				coeffOfLinearTerms[ i ][ j ][ 0 ] = knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]);
				coeffOfLinearTerms[ i ][ j ][ 1 ] = - 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1]);
			}
		}
	}
	//==================================================

	// checking for time2======================
    StartTime = clock();

	// Expanding sum of multiplications of linear terms
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		tempCoeff[ 0 ] = coeffOfLinearTerms[ i ][ 0 ][ 0 ];
		tempCoeff[ 1 ] = coeffOfLinearTerms[ i ][ 0 ][ 1 ];
		for(rg_INT j = 1;j <= p - 1;j++)
		{
			tempCoeff[j + 1] = tempCoeff[ j ] * coeffOfLinearTerms[ i ][ j ][ 1 ];
			for(rg_INT k = j;k >= 1;k--)
			{
				tempCoeff[ k ] = tempCoeff[k - 1] * coeffOfLinearTerms[ i ][ j ][ 1 ] + tempCoeff[ k ] * coeffOfLinearTerms[ i ][ j ][ 0 ];
			}
			tempCoeff[ 0 ] *= coeffOfLinearTerms[ i ][ j ][ 0 ];
		}
		for(rg_INT l = 0;l < p + 1;l++)
			coefficientOfPolynomial[ l ] += tempCoeff[ l ]; 
	}
	//==================================================

	delete[] tempCoeff;

	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		for(rg_INT j = 0;j < p;j++)
		{
			delete[] coeffOfLinearTerms[ i ][ j ];
		}
		delete[] coeffOfLinearTerms[ i ];
	}
	delete[] coeffOfLinearTerms;

	EndTime = clock();
	time2 += (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	//========================================	

	return coefficientOfPolynomial;
}


void rg_NUBSplineCurve3D::getEntirePolynomialCurve(rg_REAL** & coeffOfXCoordinate,
												rg_REAL** & coeffOfYCoordinate,
												rg_REAL** & coeffOfZCoordinate) const
{
	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT*** allPossiblePath = new rg_INT** [p + 1]; // all possible paths from (p + 1) basic graphs
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;

	// find all the possible paths from (p + 1) basic graphs
/*
	//1. special case
	for(rg_INT i = 0;i < degree;i++)
	{
		allPossiblePath[ 0 ][ i ] = 1;
		allPossiblePath[degree][ i ] = 0;
	}
*/	
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

	rg_REAL* tempCoeff;
	//rg_REAL* tempCoeffOfX;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfY;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfZ;// = new rg_REAL[p + 1];
	

	for(i = 0;i < NoOfNonZeroLengthKnotSpan;i++)
	{
		coeffOfXCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfYCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfZCoordinate[ i ] = new rg_REAL[p + 1];
		for(rg_INT j = 0;j < p + 1;j++)
		{
			/* tempCoeffOfX[ j ] = */coeffOfXCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfY[ j ] = */coeffOfYCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfZ[ j ] = */coeffOfZCoordinate[ i ][ j ] = 0.0;
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
			//tempCoeff = getDistributionPolynomialsInOneGraph(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);
			tempCoeff = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);

			for(rg_INT indexOfPolyCurveCoeff = 0;indexOfPolyCurveCoeff < p + 1;indexOfPolyCurveCoeff++)
			{
				coeffOfXCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getX() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfYCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getY() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfZCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getZ() * tempCoeff[indexOfPolyCurveCoeff]);
			}
			delete[] tempCoeff;			
		}
		indexOfPolyCurve++;
	}
	//------------------------------------------

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraph[ i ];j++)
		{
			delete[] allPossiblePath[ i ][ j ];
		}
		delete[] allPossiblePath[ i ];
	}
	delete[] allPossiblePath;
	delete[] noOfAllPossiblePathsInEachGraph;

	//delete[] allPossiblePath;
	//delete[] noOfAllPossiblePathsInEachGraph;
}


void rg_NUBSplineCurve3D::getEntirePolynomialCurve(rg_REAL** & coeffOfXCoordinate,
												rg_REAL** & coeffOfYCoordinate,
												rg_REAL** & coeffOfZCoordinate,
												rg_REAL   & timeForOurAlgorithm1,
												rg_REAL   & timeForOurAlgorithm2,
												rg_REAL   & timeForOurAlgorithm3) const
{
	// checking for time1======================

	clock_t StartTime, EndTime;
    StartTime = clock();

	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT*** allPossiblePath = new rg_INT** [p + 1]; // all possible paths from (p + 1) basic graphs
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;

	// find all the possible paths from (p + 1) basic graphs
/*
	//1. special case
	for(rg_INT i = 0;i < degree;i++)
	{
		allPossiblePath[ 0 ][ i ] = 1;
		allPossiblePath[degree][ i ] = 0;
	}
*/	
	rg_INT* noOfAllPossiblePathsInEachGraph = new rg_INT[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		allPossiblePath[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, p - i);
		noOfAllPossiblePathsInEachGraph[ i ] = rg_MathFunc::combination(p, i);
	}
	//------------------------------------------
	
	EndTime = clock();
	timeForOurAlgorithm1 += (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	//========================================

	// checking for time2 and time3======================

	StartTime = clock();
	rg_REAL timeForLinearTermComputation1;
	//rg_REAL timeForLinearTermComputation2 = 0.0;

	// initialization of coefficient of piecewise polyomial curve
	coeffOfXCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfYCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfZCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];

	rg_REAL* tempCoeff;
	//rg_REAL* tempCoeffOfX;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfY;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfZ;// = new rg_REAL[p + 1];
	

	for(i = 0;i < NoOfNonZeroLengthKnotSpan;i++)
	{
		coeffOfXCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfYCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfZCoordinate[ i ] = new rg_REAL[p + 1];
		for(rg_INT j = 0;j < p + 1;j++)
		{
			/* tempCoeffOfX[ j ] = */coeffOfXCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfY[ j ] = */coeffOfYCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfZ[ j ] = */coeffOfZCoordinate[ i ][ j ] = 0.0;
		}
	}

	// generating the piecewise polynomial curve

	rg_Point3D* ctrlPt = getCtrlPts();
	rg_INT indexOfPolyCurve = 0;//, indexOfPolyCurveCoeff = 0;

	timeForLinearTermComputation1 = 0.0;

	for(i = p;i <= m - p - 1 && indexOfPolyCurve < NoOfNonZeroLengthKnotSpan;i++)
	// loop for interior nonzero length knot spans
	{
		for(rg_INT j = i - p;j <= i;j++)
		// 
		{
			clock_t StartTime2, EndTime2;
			StartTime2 = clock();

			//tempCoeff = getDistributionPolynomialsInOneGraph(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);
			tempCoeff = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);

			EndTime2 = clock();
			timeForLinearTermComputation1 += (rg_REAL)(EndTime2 - StartTime2)/CLOCKS_PER_SEC;

			for(rg_INT indexOfPolyCurveCoeff = 0;indexOfPolyCurveCoeff < p + 1;indexOfPolyCurveCoeff++)
			{
				coeffOfXCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getX() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfYCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getY() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfZCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getZ() * tempCoeff[indexOfPolyCurveCoeff]);
			}
			delete[] tempCoeff;			
		}
		indexOfPolyCurve++;
	}
	//------------------------------------------

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraph[ i ];j++)
		{
			delete[] allPossiblePath[ i ][ j ];
		}
		delete[] allPossiblePath[ i ];
	}
	delete[] allPossiblePath;
	delete[] noOfAllPossiblePathsInEachGraph;

	//delete[] allPossiblePath;
	//delete[] noOfAllPossiblePathsInEachGraph;

	EndTime = clock();
	
	timeForOurAlgorithm2 += timeForLinearTermComputation1;
	//==========================================

	// checking for time3======================

	timeForOurAlgorithm3 += ( (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC - timeForLinearTermComputation1 );
	//==========================================
}

void rg_NUBSplineCurve3D::getEntirePolynomialCurve1(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
											    rg_REAL   & timeForOurAlgorithm1) const
{
	// checking for time1======================

	clock_t StartTime, EndTime;
    StartTime = clock();

	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT*** allPossiblePath = new rg_INT** [p + 1]; // all possible paths from (p + 1) basic graphs
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;

	// find all the possible paths from (p + 1) basic graphs
/*
	//1. special case
	for(rg_INT i = 0;i < degree;i++)
	{
		allPossiblePath[ 0 ][ i ] = 1;
		allPossiblePath[degree][ i ] = 0;
	}
*/	
	rg_INT* noOfAllPossiblePathsInEachGraph = new rg_INT[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		allPossiblePath[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, p - i);
		noOfAllPossiblePathsInEachGraph[ i ] = rg_MathFunc::combination(p, i);
	}
	//------------------------------------------
	
	EndTime = clock();
	timeForOurAlgorithm1 += (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	//========================================

	// initialization of coefficient of piecewise polyomial curve
	coeffOfXCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfYCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfZCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];

	rg_REAL* tempCoeff;
	//rg_REAL* tempCoeffOfX;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfY;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfZ;// = new rg_REAL[p + 1];
	

	for(i = 0;i < NoOfNonZeroLengthKnotSpan;i++)
	{
		coeffOfXCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfYCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfZCoordinate[ i ] = new rg_REAL[p + 1];
		for(rg_INT j = 0;j < p + 1;j++)
		{
			/* tempCoeffOfX[ j ] = */coeffOfXCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfY[ j ] = */coeffOfYCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfZ[ j ] = */coeffOfZCoordinate[ i ][ j ] = 0.0;
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
			//tempCoeff = getDistributionPolynomialsInOneGraph(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);
			tempCoeff = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);

			for(rg_INT indexOfPolyCurveCoeff = 0;indexOfPolyCurveCoeff < p + 1;indexOfPolyCurveCoeff++)
			{
				coeffOfXCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getX() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfYCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getY() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfZCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getZ() * tempCoeff[indexOfPolyCurveCoeff]);
			}
			delete[] tempCoeff;			
		}
		indexOfPolyCurve++;
	}
	//------------------------------------------

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraph[ i ];j++)
		{
			delete[] allPossiblePath[ i ][ j ];
		}
		delete[] allPossiblePath[ i ];
	}
	delete[] allPossiblePath;
	delete[] noOfAllPossiblePathsInEachGraph;

	//delete[] allPossiblePath;
	//delete[] noOfAllPossiblePathsInEachGraph;
}

void rg_NUBSplineCurve3D::getEntirePolynomialCurve2(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
        					   rg_REAL   & timeForOurAlgorithm2) const
{
	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT*** allPossiblePath = new rg_INT** [p + 1]; // all possible paths from (p + 1) basic graphs
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;

	// find all the possible paths from (p + 1) basic graphs
/*
	//1. special case
	for(rg_INT i = 0;i < degree;i++)
	{
		allPossiblePath[ 0 ][ i ] = 1;
		allPossiblePath[degree][ i ] = 0;
	}
*/	
	rg_INT* noOfAllPossiblePathsInEachGraph = new rg_INT[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		allPossiblePath[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, p - i);
		noOfAllPossiblePathsInEachGraph[ i ] = rg_MathFunc::combination(p, i);
	}
	//------------------------------------------

	// checking for time2======================
	clock_t StartTime, EndTime;

	// initialization of coefficient of piecewise polyomial curve
	coeffOfXCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfYCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfZCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];

	rg_REAL* tempCoeff;
	//rg_REAL* tempCoeffOfX;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfY;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfZ;// = new rg_REAL[p + 1];
	

	for(i = 0;i < NoOfNonZeroLengthKnotSpan;i++)
	{
		coeffOfXCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfYCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfZCoordinate[ i ] = new rg_REAL[p + 1];
		for(rg_INT j = 0;j < p + 1;j++)
		{
			/* tempCoeffOfX[ j ] = */coeffOfXCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfY[ j ] = */coeffOfYCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfZ[ j ] = */coeffOfZCoordinate[ i ][ j ] = 0.0;
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
			StartTime = clock();

			//tempCoeff = getDistributionPolynomialsInOneGraph(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);
			tempCoeff = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);

			EndTime = clock();
			timeForOurAlgorithm2 += (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;

			for(rg_INT indexOfPolyCurveCoeff = 0;indexOfPolyCurveCoeff < p + 1;indexOfPolyCurveCoeff++)
			{
				coeffOfXCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getX() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfYCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getY() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfZCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getZ() * tempCoeff[indexOfPolyCurveCoeff]);
			}
			delete[] tempCoeff;			
		}
		indexOfPolyCurve++;
	}
	//------------------------------------------

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraph[ i ];j++)
		{
			delete[] allPossiblePath[ i ][ j ];
		}
		delete[] allPossiblePath[ i ];
	}
	delete[] allPossiblePath;
	delete[] noOfAllPossiblePathsInEachGraph;

	//delete[] allPossiblePath;
	//delete[] noOfAllPossiblePathsInEachGraph;

	
	//==========================================
}

void rg_NUBSplineCurve3D::getEntirePolynomialCurve3(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
                               rg_REAL   & timeForOurAlgorithm3) const
{
	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT*** allPossiblePath = new rg_INT** [p + 1]; // all possible paths from (p + 1) basic graphs
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;

	// find all the possible paths from (p + 1) basic graphs
/*
	//1. special case
	for(rg_INT i = 0;i < degree;i++)
	{
		allPossiblePath[ 0 ][ i ] = 1;
		allPossiblePath[degree][ i ] = 0;
	}
*/	
	rg_INT* noOfAllPossiblePathsInEachGraph = new rg_INT[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		allPossiblePath[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, p - i);
		noOfAllPossiblePathsInEachGraph[ i ] = rg_MathFunc::combination(p, i);
	}
	//------------------------------------------

	// checking for time3======================
	clock_t StartTime, EndTime;
	StartTime = clock();
	rg_REAL timeForLinearTermComputation1;	
	//rg_REAL timeForLinearTermComputation2 = 0.0;

	// initialization of coefficient of piecewise polyomial curve
	coeffOfXCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfYCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfZCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];

	rg_REAL* tempCoeff;
	//rg_REAL* tempCoeffOfX;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfY;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfZ;// = new rg_REAL[p + 1];
	

	for(i = 0;i < NoOfNonZeroLengthKnotSpan;i++)
	{
		coeffOfXCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfYCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfZCoordinate[ i ] = new rg_REAL[p + 1];
		for(rg_INT j = 0;j < p + 1;j++)
		{
			/* tempCoeffOfX[ j ] = */coeffOfXCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfY[ j ] = */coeffOfYCoordinate[ i ][ j ] = 0.0;
			/* tempCoeffOfZ[ j ] = */coeffOfZCoordinate[ i ][ j ] = 0.0;
		}
	}

	// generating the piecewise polynomial curve

	rg_Point3D* ctrlPt = getCtrlPts();
	rg_INT indexOfPolyCurve = 0;//, indexOfPolyCurveCoeff = 0;

	timeForLinearTermComputation1 = 0.0;

	for(i = p;i <= m - p - 1 && indexOfPolyCurve < NoOfNonZeroLengthKnotSpan;i++)
	// loop for interior nonzero length knot spans
	{
		for(rg_INT j = i - p;j <= i;j++)
		// 
		{
			clock_t StartTime2, EndTime2;
			StartTime2 = clock();

			//tempCoeff = getDistributionPolynomialsInOneGraph1(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p],
			//	                                            timeForOurAlgorithm3);//timeForLinearTermComputation1);//, timeForLinearTermComputation2);
			tempCoeff = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);

			EndTime2 = clock();
			timeForLinearTermComputation1 += (rg_REAL)(EndTime2 - StartTime2)/CLOCKS_PER_SEC;

			for(rg_INT indexOfPolyCurveCoeff = 0;indexOfPolyCurveCoeff < p + 1;indexOfPolyCurveCoeff++)
			{
				coeffOfXCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getX() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfYCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getY() * tempCoeff[indexOfPolyCurveCoeff]);
				coeffOfZCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getZ() * tempCoeff[indexOfPolyCurveCoeff]);
			}
			delete[] tempCoeff;			
		}
		indexOfPolyCurve++;
	}
	//------------------------------------------

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraph[ i ];j++)
		{
			delete[] allPossiblePath[ i ][ j ];
		}
		delete[] allPossiblePath[ i ];
	}
	delete[] allPossiblePath;
	delete[] noOfAllPossiblePathsInEachGraph;

	//delete[] allPossiblePath;
	//delete[] noOfAllPossiblePathsInEachGraph;

	EndTime = clock();
	

	// checking for time3======================

	timeForOurAlgorithm3 += ( (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC - timeForLinearTermComputation1 );
	//==========================================
}

rg_REAL rg_NUBSplineCurve3D::getCoeffOfMatrixRepresentation(rg_INT multiplicity, rg_INT row, rg_INT col, rg_INT indexOfInsertedKnot, rg_INT indexOfCtrlPt)
{
	rg_INT r = multiplicity;
	rg_INT k = indexOfInsertedKnot;
	rg_INT p = getOrder() - 1;
	rg_INT i = indexOfCtrlPt;

	if(r == 0)
	{
		if((k - p + col == indexOfCtrlPt) && (row == 0))
			return 1.0;
		else
			return 0.0;
	}
	else
	{
		return  1.0 / (knotVector[i + p - r + 1] - knotVector[ i ]) * getCoeffOfMatrixRepresentation(r - 1, row - 1, col, k, i)
		        + - knotVector[ i ] / (knotVector[i + p - r + 1] - knotVector[ i ]) * getCoeffOfMatrixRepresentation(r - 1, row, col, k, i)
				+ knotVector[i + p - r + 1] / (knotVector[i + p - r + 1] - knotVector[ i ]) * getCoeffOfMatrixRepresentation(r - 1, row, col, k, i - 1)
				+ - 1.0 / (knotVector[i + p - r + 1] - knotVector[ i ]) * getCoeffOfMatrixRepresentation(r - 1, row - 1, col, k, i - 1);
	}
}

// Assume that the multiplicity of the interior knot is 1.
rg_Matrix* rg_NUBSplineCurve3D::getMatrixRepresentation()
{
	rg_INT NumOfNonZeroLengthKnotSpan = getNumOfNonZeroLengthKnotSpan();
	rg_DEGREE p = getOrder() - 1;

	rg_Matrix* coeffMatrixList = new rg_Matrix[NumOfNonZeroLengthKnotSpan];

	rg_INT i = 0;
	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
		coeffMatrixList[ i ].setSize(p + 1, p + 1);

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		for(rg_INT j = 0;j < p + 1;j++)
		{
			for(rg_INT k = 0;k < p + 1;k++)
			{
				coeffMatrixList[ i ].setElement(j, k, getCoeffOfMatrixRepresentation(p, j, k, i + p, i + p));
			}
		}
	}

	return coeffMatrixList;
}

rg_REAL rg_NUBSplineCurve3D::getCoeffOfMatrixRepresentationInUniformKnot(rg_INT row, rg_INT col)
{
	rg_REAL coefficient = 0.0;
	rg_INT p = getOrder() - 1;

	for(rg_INT i = col;i <= p;i++)
	{
		coefficient += ( pow((double)(p - i), row) * pow(-1., i - col) * rg_MathFunc::combination(p + 1, i - col) );
	}

	coefficient *= (1 / (rg_MathFunc::factorial(row) * rg_MathFunc::factorial(p - row)));

	return coefficient;
}

// By Hoschek and Yamaguchi 's

rg_Matrix* rg_NUBSplineCurve3D::getMatrixRepresentationInUniformKnot()
{
	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1

	rg_INT NumOfNonZeroLengthKnotSpan = m - 2 * p;

	rg_Matrix* coeffMatrixList = new rg_Matrix[NumOfNonZeroLengthKnotSpan];

	rg_INT i = 0;
	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
		coeffMatrixList[ i ].setSize(p + 1, p + 1);

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		for(rg_INT j = 0;j < p + 1;j++)
		{
			for(rg_INT k = 0;k < p + 1;k++)
			{
				coeffMatrixList[ i ].setElement(j, k, getCoeffOfMatrixRepresentationInUniformKnot(j, k));
					//(knotVector[i + 1] - knotVector[ i ]) * getCoeffOfMatrixRepresentationInUniformKnot(j, k) + knotVector[ i ]);

			}
		}
        rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
                                       1/(knotVector[i + p + 1] - knotVector[ i  + p ]),
	                                   - knotVector[ i  + p ] /(knotVector[i  + p + 1] - knotVector[ i  + p ]));

		coeffMatrixList[ i ] = R_p * coeffMatrixList[ i ];
	}

	return coeffMatrixList;
}


rg_Matrix* rg_NUBSplineCurve3D::getMatrixRepresentationByDE()
{
	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;
	rg_INT*** allPossiblePath = new rg_INT** [p + 1]; // all possible paths from (p + 1) basic graphs

	rg_INT* noOfAllPossiblePathsInEachGraph = new rg_INT[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		allPossiblePath[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, p - i);
		noOfAllPossiblePathsInEachGraph[ i ] = rg_MathFunc::combination(p, i);
	}

	rg_Matrix* coeffMatrixList = new rg_Matrix[NoOfNonZeroLengthKnotSpan];

	for(i = 0;i < NoOfNonZeroLengthKnotSpan;i++)
		coeffMatrixList[ i ].setSize(p + 1, p + 1);

	rg_INT indexOfPolyCurve = 0;

	for(i = p;i <= m - p - 1;i++, indexOfPolyCurve++)
	{
		rg_INT col = 0;
		rg_REAL* tempCoeff;

		for(rg_INT j = i - p;j <= i;j++, col++)
		{
			tempCoeff = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);

			for(rg_INT row = 0;row < p + 1;row++)
			{
				coeffMatrixList[indexOfPolyCurve].setElement(row, col, tempCoeff[row]);
			}
			delete[] tempCoeff;			
		}
	}

	return coeffMatrixList;
}

// Using Taylor's series

void rg_NUBSplineCurve3D::getPiecewisePolynomialCurveInaPowerForm(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate)
{
	// Calculate derivative of curve up to degree p

	rg_INT p = getOrder() - 1; // degree

	rg_NUBSplineCurve3D* derivatives = new rg_NUBSplineCurve3D [ p ];

	derivatives[ 0 ] = makeDerivative();

	rg_INT i = 0;
	for(i = 0;i < p - 1;i++)
		derivatives[i + 1] = derivatives[ i ].makeDerivative();

	// -----------------------------------------------

	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;

	rg_REAL** TaylorExpansionInMiddle = new rg_REAL* [ 3 ];

	for(i = 0;i < 3;i++)
		TaylorExpansionInMiddle[ i ] = new rg_REAL[p + 1];

	//-----------------------------------------------

	// initialization of coefficient of piecewise polyomial curve
	coeffOfXCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfYCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];
	coeffOfZCoordinate = new rg_REAL* [NoOfNonZeroLengthKnotSpan];

	for(i = 0;i < NoOfNonZeroLengthKnotSpan;i++)
	{
		coeffOfXCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfYCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfZCoordinate[ i ] = new rg_REAL[p + 1];
		/*
		for(rg_INT j = 0;j < p + 1;j++)
		{
			coeffOfXCoordinate[ i ][ j ] = 0.0;
			coeffOfYCoordinate[ i ][ j ] = 0.0;
			coeffOfZCoordinate[ i ][ j ] = 0.0;
		}
		*/
	}


	// For each knot span of nonzero length
	for(i = p;i <= m - p - 1;i++)
	{
		rg_REAL middleKnot = (knotVector[ i ] + knotVector[i + 1]) / 2;

		TaylorExpansionInMiddle[ 0 ][ 0 ] = evaluatePt(middleKnot).getX(); // X-coordinate
		TaylorExpansionInMiddle[ 1 ][ 0 ] = evaluatePt(middleKnot).getY(); // Y-coordinate
		TaylorExpansionInMiddle[ 2 ][ 0 ] = evaluatePt(middleKnot).getZ(); // Z-coordinate

		rg_INT j = 0;
		for(j = 1;j <= p;j++)
		{
			rg_Point3D derivative = derivatives[j - 1].evaluatePt(middleKnot);

			TaylorExpansionInMiddle[ 0 ][ j ] = derivative.getX() / rg_MathFunc::factorial( j );
			TaylorExpansionInMiddle[ 1 ][ j ] = derivative.getY() / rg_MathFunc::factorial( j );
			TaylorExpansionInMiddle[ 2 ][ j ] = derivative.getZ() / rg_MathFunc::factorial( j );
		}

		// Extracting the coefficients of each polynomial curve
		for(j = 0;j < p;j++)
		{
			for(rg_INT k = p - 1;k >= j;k--)
			{
				TaylorExpansionInMiddle[ 0 ][ k ] += (TaylorExpansionInMiddle[ 0 ][k + 1] * (- middleKnot));
				TaylorExpansionInMiddle[ 1 ][ k ] += (TaylorExpansionInMiddle[ 1 ][k + 1] * (- middleKnot));
				TaylorExpansionInMiddle[ 2 ][ k ] += (TaylorExpansionInMiddle[ 2 ][k + 1] * (- middleKnot));
			}
		}

		for(j = 0;j <= p;j++)
		{
			coeffOfXCoordinate[i - p][ j ] = TaylorExpansionInMiddle[ 0 ][ j ];
			coeffOfYCoordinate[i - p][ j ] = TaylorExpansionInMiddle[ 1 ][ j ];
			coeffOfZCoordinate[i - p][ j ] = TaylorExpansionInMiddle[ 2 ][ j ];
		}
	}

	for(i = 0;i < 3;i++)
		delete[] TaylorExpansionInMiddle[ i ];

	delete[] TaylorExpansionInMiddle;

	delete[] derivatives;
}

rg_Polynomial* rg_NUBSplineCurve3D::getEntirePolynomialOfX() const
{
	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT*** allPossiblePath = new rg_INT** [p + 1]; // all possible paths from (p + 1) basic graphs
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;

	// find all the possible paths from (p + 1) basic graphs
/*
	//1. special case
	for(rg_INT i = 0;i < degree;i++)
	{
		allPossiblePath[ 0 ][ i ] = 1;
		allPossiblePath[degree][ i ] = 0;
	}
*/	
	rg_INT* noOfAllPossiblePathsInEachGraph = new rg_INT[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		allPossiblePath[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, p - i);
		noOfAllPossiblePathsInEachGraph[ i ] = rg_MathFunc::combination(p, i);
	}
	//------------------------------------------

    const rg_INT dimension=3;
	// initialization of coefficient of piecewise polyomial curve
	rg_Polynomial* output = new rg_Polynomial[NoOfNonZeroLengthKnotSpan];
    for (  i=0; i < NoOfNonZeroLengthKnotSpan; i++)
    {
        output[i].setDegree(p);
    }
	rg_REAL* tempCoeff;
	//rg_REAL* tempCoeffOfX;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfY;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfZ;// = new rg_REAL[p + 1];
	
	// generating the piecewise polynomial curve

	rg_Point3D* ctrlPt = getCtrlPts();
	rg_INT indexOfPolyCurve = 0;//, indexOfPolyCurveCoeff = 0;

	for(i = p;i <= m - p - 1 && indexOfPolyCurve < NoOfNonZeroLengthKnotSpan;i++)
	// loop for interior nonzero length knot spans
	{
		for(rg_INT j = i - p;j <= i;j++)
		// 
		{
			//tempCoeff = getDistributionPolynomialsInOneGraph(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);
			tempCoeff = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);

			for(rg_INT indexOfPolyCurveCoeff = 0;indexOfPolyCurveCoeff < p + 1;indexOfPolyCurveCoeff++)
			{
				output[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getX() * tempCoeff[indexOfPolyCurveCoeff]);
			}
			delete[] tempCoeff;			
		}
		indexOfPolyCurve++;
	}
	//------------------------------------------

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraph[ i ];j++)
		{
			delete[] allPossiblePath[ i ][ j ];
		}
		delete[] allPossiblePath[ i ];
	}
	delete[] allPossiblePath;
	delete[] noOfAllPossiblePathsInEachGraph;

	//delete[] allPossiblePath;
	//delete[] noOfAllPossiblePathsInEachGraph;
    return output;
}

rg_Polynomial** rg_NUBSplineCurve3D::getEntirePolynomialCurve() const
{
	rg_INT p = getOrder() - 1;                     // degree
	rg_INT n = getNumOfCtrlPts() - 1;      // No. of control points - 1
	rg_INT m = n + p + 1;                          // no. of knotvector - 1
	rg_INT*** allPossiblePath = new rg_INT** [p + 1]; // all possible paths from (p + 1) basic graphs
	rg_INT NoOfNonZeroLengthKnotSpan = m - 2 * p;

	// find all the possible paths from (p + 1) basic graphs
/*
	//1. special case
	for(rg_INT i = 0;i < degree;i++)
	{
		allPossiblePath[ 0 ][ i ] = 1;
		allPossiblePath[degree][ i ] = 0;
	}
*/	
	rg_INT* noOfAllPossiblePathsInEachGraph = new rg_INT[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		allPossiblePath[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, p - i);
		noOfAllPossiblePathsInEachGraph[ i ] = rg_MathFunc::combination(p, i);
	}
	//------------------------------------------

    const rg_INT dimension=3;
	// initialization of coefficient of piecewise polyomial curve
	rg_Polynomial** output = new rg_Polynomial*[NoOfNonZeroLengthKnotSpan];
    for (  i=0; i < NoOfNonZeroLengthKnotSpan; i++)
    {
        output[i]= new rg_Polynomial[dimension];
        for( rg_INT j=0; j < dimension; j++ )
        {
            output[i][j].setDegree(p);
        }
    }
	rg_REAL* tempCoeff;
	//rg_REAL* tempCoeffOfX;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfY;// = new rg_REAL[p + 1];
	//rg_REAL* tempCoeffOfZ;// = new rg_REAL[p + 1];
	
	// generating the piecewise polynomial curve

	rg_Point3D* ctrlPt = getCtrlPts();
	rg_INT indexOfPolyCurve = 0;//, indexOfPolyCurveCoeff = 0;

	for(i = p;i <= m - p - 1 && indexOfPolyCurve < NoOfNonZeroLengthKnotSpan;i++)
	// loop for interior nonzero length knot spans
	{
		for(rg_INT j = i - p;j <= i;j++)
		// 
		{
			//tempCoeff = getDistributionPolynomialsInOneGraph(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);
			tempCoeff = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);

			for(rg_INT indexOfPolyCurveCoeff = 0;indexOfPolyCurveCoeff < p + 1;indexOfPolyCurveCoeff++)
			{
				output[indexOfPolyCurve][0][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getX() * tempCoeff[indexOfPolyCurveCoeff]);
                output[indexOfPolyCurve][1][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getY() * tempCoeff[indexOfPolyCurveCoeff]);
				output[indexOfPolyCurve][2][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getZ() * tempCoeff[indexOfPolyCurveCoeff]);

			}
			delete[] tempCoeff;			
		}
		indexOfPolyCurve++;
	}
	//------------------------------------------

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraph[ i ];j++)
		{
			delete[] allPossiblePath[ i ][ j ];
		}
		delete[] allPossiblePath[ i ];
	}
	delete[] allPossiblePath;
	delete[] noOfAllPossiblePathsInEachGraph;

	//delete[] allPossiblePath;
	//delete[] noOfAllPossiblePathsInEachGraph;
    return output;
}

void rg_NUBSplineCurve3D::makePiecewisePowerFormPolyUsingDE(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate, rg_REAL*** & truncatedBasisPolys)
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
	truncatedBasisPolys = new rg_REAL** [NoOfNonZeroLengthKnotSpan];

	for(i = 0;i < NoOfNonZeroLengthKnotSpan;i++)
	{
		coeffOfXCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfYCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfZCoordinate[ i ] = new rg_REAL[p + 1];
		truncatedBasisPolys[ i ] = new rg_REAL* [p + 1];
		for(rg_INT j = 0;j < p + 1;j++)
		{
			coeffOfXCoordinate[ i ][ j ] = 0.0;
			coeffOfYCoordinate[ i ][ j ] = 0.0;
			coeffOfZCoordinate[ i ][ j ] = 0.0;
		}
	}

	// generating the piecewise polynomial curve

	rg_Point3D* ctrlPt = getCtrlPts();
	rg_INT indexOfPolyCurve = 0;

	for(i = p;i <= m - p - 1 && indexOfPolyCurve < NoOfNonZeroLengthKnotSpan;i++)
	// loop for interior nonzero length knot spans
	{
		for(rg_INT j = i - p;j <= i;j++)
		// 
		{
			makeTruncatedBasisPolynomial(truncatedBasisPolys[i - p][j - (i - p)], i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);
			//truncatedBasisPolys[i - p][ j ] = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + p], noOfAllPossiblePathsInEachGraph[j - i + p]);
			for(rg_INT indexOfPolyCurveCoeff = 0;indexOfPolyCurveCoeff < p + 1;indexOfPolyCurveCoeff++)
			{
				coeffOfXCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getX() * truncatedBasisPolys[i - p][j - (i - p)][indexOfPolyCurveCoeff]);
				coeffOfYCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getY() * truncatedBasisPolys[i - p][j - (i - p)][indexOfPolyCurveCoeff]);
				coeffOfZCoordinate[indexOfPolyCurve][indexOfPolyCurveCoeff] += (ctrlPt[ j ].getZ() * truncatedBasisPolys[i - p][j - (i - p)][indexOfPolyCurveCoeff]);
			}
		}
		indexOfPolyCurve++;
	}
	//------------------------------------------

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraph[ i ];j++)
		{
			delete[] allPossiblePath[ i ][ j ];
		}
		delete[] allPossiblePath[ i ];
	}
	delete[] allPossiblePath;
	delete[] noOfAllPossiblePathsInEachGraph;

	//delete[] allPossiblePath;
	//delete[] noOfAllPossiblePathsInEachGraph;
}

void rg_NUBSplineCurve3D::makeTruncatedBasisPolynomial(rg_REAL* & truncatedBasisPolynomial, rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths)
{
	rg_INT p = getOrder() - 1;

	truncatedBasisPolynomial = new rg_REAL[p + 1];
	rg_REAL* tempCoeff = new rg_REAL[p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
		truncatedBasisPolynomial[ i ] = 0.0;

	// Expanding sum of multiplications of linear terms according to basic graph
	for(i = 0;i < noOfAllPossiblePaths;i++)
	{
		rg_INT degree = 0;
		rg_INT index = contributionKnotSpan;

		degree++;

		if(graph[ i ][ 0 ] == 0)
		{
			tempCoeff[ 0 ] = - knotVector[index] / (knotVector[index + degree] - knotVector[index]);
			tempCoeff[ 1 ] = 1.0 / (knotVector[index + degree] - knotVector[index]);
		}
		else
		{
			index--;
			tempCoeff[ 0 ] = knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]);
			tempCoeff[ 1 ] = - 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1]);
		}
		for(rg_INT j = 1;j <= p - 1;j++)
		{
			degree++;

			if(graph[ i ][ j ] == 0)
			{
				tempCoeff[j + 1] = tempCoeff[ j ] * (1.0 / (knotVector[index + degree] - knotVector[index]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempCoeff[ k ] = tempCoeff[k - 1] * (1.0 / (knotVector[index + degree] - knotVector[index])) + tempCoeff[ k ] * (- knotVector[index] / (knotVector[index + degree] - knotVector[index]));
				}
				tempCoeff[ 0 ] *= (- knotVector[index] / (knotVector[index + degree] - knotVector[index]));
			}
			else
			{
				index--;
				tempCoeff[j + 1] = tempCoeff[ j ] * (- 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempCoeff[ k ] = tempCoeff[k - 1] * (- 1.0 / (knotVector[index + degree + 1] - knotVector[index + 1])) + tempCoeff[ k ] * (knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]));
				}
				tempCoeff[ 0 ] *= (knotVector[index + degree + 1] / (knotVector[index + degree + 1] - knotVector[index + 1]));
			}
		}
		for(rg_INT l = 0;l < p + 1;l++)
			truncatedBasisPolynomial[ l ] += tempCoeff[ l ]; 
	}

	delete[] tempCoeff;
}

void rg_NUBSplineCurve3D::makeUpdatedCurveUsingDE(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate, rg_REAL*** & truncatedBasisPolys,
		                                       rg_INT* indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts)
{
	rg_INT p = getOrder() - 1;               
	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT m = n + p + 1;
	
	rg_INT i = 0;
	for(i = 0;i < numOfChangedCtrlPts;i++)
	{
		for(rg_INT j = indexOfChangedCtrlPts[ i ];j < indexOfChangedCtrlPts[ i ] + p + 1;j++)
		{
			//if(rg_NZERO(knotVector[j + 1] - knotVector[ j ]))
			if((p <= j) && (j <= m - p - 1))
			{
				rg_Point3D ctrlPt = getCtrlPt(indexOfChangedCtrlPts[ i ]);
				for(rg_INT k = 0;k < p + 1;k++)
				{
					coeffOfXCoordinate[j - p][ k ] += (ctrlPt.getX() * truncatedBasisPolys[j - p][indexOfChangedCtrlPts[ i ] - (j - p)][ k ]);
					coeffOfYCoordinate[j - p][ k ] += (ctrlPt.getY() * truncatedBasisPolys[j - p][indexOfChangedCtrlPts[ i ] - (j - p)][ k ]);
					coeffOfZCoordinate[j - p][ k ] += (ctrlPt.getZ() * truncatedBasisPolys[j - p][indexOfChangedCtrlPts[ i ] - (j - p)][ k ]);
				}
			}
		}
	}
}

void rg_NUBSplineCurve3D::makeUpdatedCurveUsingTE(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
	                                           rg_INT* indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts)
{
	for(rg_INT i = 0;i < numOfChangedCtrlPts;i++)
	{
		rg_INT p = getOrder() - 1;               
		rg_INT n = getNumOfCtrlPts() - 1;
		rg_INT m = n + p + 1;

		// Calculate derivative of curve up to degree p

		rg_NUBSplineCurve3D* derivatives = new rg_NUBSplineCurve3D [ p ];

		derivatives[ 0 ] = makeDerivative();

		rg_INT j = 0;
		for(j = 0;j < p - 1;j++)
			derivatives[j + 1] = derivatives[ j ].makeDerivative();

		// -----------------------------------------------
		rg_REAL** TaylorExpansionInMiddle = new rg_REAL* [ 3 ];

		for(j = 0;j < 3;j++)
			TaylorExpansionInMiddle[ j ] = new rg_REAL[p + 1];

		for(j = indexOfChangedCtrlPts[ i ];j < indexOfChangedCtrlPts[ i ] + p + 1;j++)
		{
			if(rg_NZERO(knotVector[j + 1] - knotVector[ j ]))
			{
				rg_REAL middleKnot = (knotVector[ j ] + knotVector[j + 1]) / 2;

				TaylorExpansionInMiddle[ 0 ][ 0 ] = evaluatePt(middleKnot).getX(); // X-coordinate
				TaylorExpansionInMiddle[ 1 ][ 0 ] = evaluatePt(middleKnot).getY(); // Y-coordinate
				TaylorExpansionInMiddle[ 2 ][ 0 ] = evaluatePt(middleKnot).getZ(); // Z-coordinate

				rg_INT k = 0;
				for(k = 1;k <= p;k++)
				{
					rg_Point3D derivative = derivatives[k - 1].evaluatePt(middleKnot);

					TaylorExpansionInMiddle[ 0 ][ k ] = derivative.getX() / rg_MathFunc::factorial( k );
					TaylorExpansionInMiddle[ 1 ][ k ] = derivative.getY() / rg_MathFunc::factorial( k );
					TaylorExpansionInMiddle[ 2 ][ k ] = derivative.getZ() / rg_MathFunc::factorial( k );
				}

				// Extracting the coefficients of each polynomial curve
				for(k = 0;k < p;k++)
				{
					for(rg_INT l = p - 1;l >= k;l--)
					{
						TaylorExpansionInMiddle[ 0 ][ l ] += (TaylorExpansionInMiddle[ 0 ][l + 1] * (- middleKnot));
						TaylorExpansionInMiddle[ 1 ][ l ] += (TaylorExpansionInMiddle[ 1 ][l + 1] * (- middleKnot));
						TaylorExpansionInMiddle[ 2 ][ l ] += (TaylorExpansionInMiddle[ 2 ][l + 1] * (- middleKnot));
					}
				}

				for(k = 0;k <= p;k++)
				{
					coeffOfXCoordinate[j - p][ k ] = TaylorExpansionInMiddle[ 0 ][ k ];
					coeffOfYCoordinate[j - p][ k ] = TaylorExpansionInMiddle[ 1 ][ k ];
					coeffOfZCoordinate[j - p][ k ] = TaylorExpansionInMiddle[ 2 ][ k ];
				}
			}
		}
		delete[] derivatives;

		for(j = 0;j < 3;j++)
			delete[] TaylorExpansionInMiddle[ j ];
		delete[] TaylorExpansionInMiddle;
	}
}

void rg_NUBSplineCurve3D::makePiecewisePowerFormPolyUsingKR(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate)
{
	rg_INT NumOfNonZeroLengthKnotSpan = getNumOfNonZeroLengthKnotSpan();
	rg_DEGREE p = getOrder() - 1;

	rg_BzCurve3D* BezierCurveList = decomposeCurveIntoBezierSegmentUsingKnotRefinement();

	rg_REAL* distinctKnotValue = getDistinctKnotValues();
	rg_Matrix P(p + 1, 3);
	rg_Matrix M_p = rg_CurveSurfaceFunc::bezierToPowerMatrix(p + 1);

	coeffOfXCoordinate = new rg_REAL* [NumOfNonZeroLengthKnotSpan];
	coeffOfYCoordinate = new rg_REAL* [NumOfNonZeroLengthKnotSpan];
	coeffOfZCoordinate = new rg_REAL* [NumOfNonZeroLengthKnotSpan];

	rg_INT i = 0;
	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		coeffOfXCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfYCoordinate[ i ] = new rg_REAL[p + 1];
		coeffOfZCoordinate[ i ] = new rg_REAL[p + 1];
	}

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
		                               1/(distinctKnotValue[i + 1] - distinctKnotValue[ i ]),
			                           - distinctKnotValue[ i ] /(distinctKnotValue[i + 1] - distinctKnotValue[ i ]));
		rg_INT j = 0;
		for(j = 0;j < p + 1;j++)
		{
			P[ j ][ 0 ] = (BezierCurveList[ i ].getCtrlPt( j )).getX(); // X coordinate of ctrl pt
			P[ j ][ 1 ] = (BezierCurveList[ i ].getCtrlPt( j )).getY(); // Y coordinate of ctrl pt
			P[ j ][ 2 ] = (BezierCurveList[ i ].getCtrlPt( j )).getZ(); // Z coordinate of ctrl pt
		}

		rg_Matrix coefficient = R_p * M_p * P; // Coefficient of power basis
										    // (p + 1) by 3 matrix

		for(j = 0;j < p + 1;j++)
		{
			coeffOfXCoordinate[ i ][ j ] = coefficient[ j ][ 0 ];
			coeffOfYCoordinate[ i ][ j ] = coefficient[ j ][ 1 ];
			coeffOfZCoordinate[ i ][ j ] = coefficient[ j ][ 2 ];
		}
	}

	delete[] BezierCurveList;
}

void rg_NUBSplineCurve3D::makeUpdatedCurveUsingKR(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
	                                           rg_INT* indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts)
{

	rg_INT p = getOrder() - 1;               
	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT m = n + p + 1;

	//rg_NUBSplineCurve3D duplicatedNUBSpline = (*this);

	rg_INT i = 0;
	for(i = 0;i < numOfChangedCtrlPts;i++)
	{
		//rg_NUBSplineCurve3D duplicatedNUBSpline = (*this);

		//rg_INT upperBnd = (rg_INT)rg_Min(indexOfChangedCtrlPts[ i ] + p + 1, n + p + 1);
		//rg_INT upperBnd = (rg_INT)rg_Min(indexOfChangedCtrlPts[ i ] + p + 1, m - p - 1);
		rg_REAL* distinctKnotValue = getDistinctKnotValues();

		for(rg_INT j = indexOfChangedCtrlPts[ i ];j < indexOfChangedCtrlPts[ i ] + p + 1;j++)
		//for(rg_INT j = indexOfChangedCtrlPts[ i ];j < indexOfChangedCtrlPts[ i ] + p;j++)
		//for(rg_INT j = indexOfChangedCtrlPts[ i ];j < upperBnd;j++)
		{
			//rg_NUBSplineCurve3D duplicatedNUBSpline = (*this);
/*
			if(rg_NZERO(knotVector[j + 1] - knotVector[ j ]))
			{
				rg_NUBSplineCurve3D duplicatedNUBSpline = (*this);

				rg_INT multiplicity = 0;

				for(rg_INT k = p;k <= m - p -1;k++)
				{
					if(rg_POS(knotVector[ k ] - knotVector[j + 1]))
					{
						break;
					}
					else
					{
						if(rg_ZERO(knotVector[j + 1] - knotVector[ k ]))
						{
							multiplicity++;
						}
					}
				}

				rg_INT numOfInsertingKnotValues = p - multiplicity;
				rg_REAL* insertingKnotValues = new rg_REAL [numOfInsertingKnotValues];

				for(k = 0;k < numOfInsertingKnotValues;k++)
					insertingKnotValues[ k ] = knotVector[j + 1];

				if(numOfInsertingKnotValues > 0)
					duplicatedNUBSpline.knotRefinement(insertingKnotValues, numOfInsertingKnotValues);

				rg_REAL* distinctKnotValue = getDistinctKnotValues();
				rg_Matrix P(p + 1, 3);
				rg_Matrix M_p = rg_CurveSurfaceFunc::bezierToPowerMatrix(p + 1);

				rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
											   1/(distinctKnotValue[j - p + 1] - distinctKnotValue[ j - p ]),
											   - distinctKnotValue[ j - p ] /(distinctKnotValue[j - p + 1] - distinctKnotValue[ j - p ]));

				for(k = 0;k < p + 1;k++)
				{
					P[ k ][ 0 ] = (duplicatedNUBSpline.getCtrlPt((j - p) * p + k)).getX(); // X coordinate of ctrl pt
					P[ k ][ 1 ] = (duplicatedNUBSpline.getCtrlPt((j - p) * p + k)).getY(); // Y coordinate of ctrl pt
					P[ k ][ 2 ] = (duplicatedNUBSpline.getCtrlPt((j - p) * p + k)).getZ(); // Z coordinate of ctrl pt
				}
				rg_Matrix coefficient = R_p * M_p * P; // Coefficient of power basis
													// (p + 1) by 3 matrix
				for(k = 0;k < p + 1;k++)
				{
					coeffOfXCoordinate[j - p][ k ] = coefficient[ k ][ 0 ];
					coeffOfYCoordinate[j - p][ k ] = coefficient[ k ][ 1 ];
					coeffOfZCoordinate[j - p][ k ] = coefficient[ k ][ 2 ];
				}
				delete[] insertingKnotValues;
			}
  */
			if(rg_NZERO(knotVector[j + 1] - knotVector[ j ]))
			{
				rg_NUBSplineCurve3D duplicatedNUBSpline = (*this);
				rg_INT multiplicity1 = 1;

				rg_INT k = j;
				while(k >= 1)
				{
					if(rg_ZERO(knotVector[ j ] - knotVector[k - 1]))
						multiplicity1++;
					else
						break;
					k--;
				}

				rg_INT multiplicity2 = 1;

				k = j + 1;
				while(k <= m - p - 2)
				{
					if(rg_ZERO(knotVector[j + 1] - knotVector[k + 1]))
						multiplicity2++;
					else
						break;
					k++;
				}

				rg_INT numOfInsertingKnotValues1 = rg_Max(0, p - multiplicity1);
				rg_INT numOfInsertingKnotValues2 = rg_Max(0, p - multiplicity2);

				rg_REAL* insertingKnotValues1 = new rg_REAL [numOfInsertingKnotValues1];
				rg_REAL* insertingKnotValues2 = new rg_REAL [numOfInsertingKnotValues2];

				for(k = 0;k < numOfInsertingKnotValues1;k++)
					insertingKnotValues1[ k ] = knotVector[ j ];
				for(k = 0;k < numOfInsertingKnotValues2;k++)
					insertingKnotValues2[ k ] = knotVector[j + 1];

				if(numOfInsertingKnotValues1 > 0)
					duplicatedNUBSpline.knotRefinement(insertingKnotValues1, numOfInsertingKnotValues1);
				if(numOfInsertingKnotValues2 > 0)
					duplicatedNUBSpline.knotRefinement(insertingKnotValues2, numOfInsertingKnotValues2);

				//rg_REAL* distinctKnotValue = getDistinctKnotValues();
				rg_Matrix P(p + 1, 3);
				rg_Matrix M_p = rg_CurveSurfaceFunc::bezierToPowerMatrix(p + 1);

				rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
											   1/(distinctKnotValue[j - p + 1] - distinctKnotValue[ j - p ]),
											   - distinctKnotValue[ j - p ] /(distinctKnotValue[j - p + 1] - distinctKnotValue[ j - p ]));

				// determine the index of starting Bezier points
				rg_INT indexOfStartCtrlPt;
				
				if(j <= p + 1)
					indexOfStartCtrlPt = (j - p) * p;
				else
					indexOfStartCtrlPt = (j - p) + (p - 1); // fucking index!!

				for(k = 0;k < p + 1;k++)
				{
					/*
					P[ k ][ 0 ] = (duplicatedNUBSpline.getCtrlPt((j - p) * p + k)).getX(); // X coordinate of ctrl pt
					P[ k ][ 1 ] = (duplicatedNUBSpline.getCtrlPt((j - p) * p + k)).getY(); // Y coordinate of ctrl pt
					P[ k ][ 2 ] = (duplicatedNUBSpline.getCtrlPt((j - p) * p + k)).getZ(); // Z coordinate of ctrl pt
					*/
					P[ k ][ 0 ] = (duplicatedNUBSpline.getCtrlPt(indexOfStartCtrlPt + k)).getX(); // X coordinate of ctrl pt
					P[ k ][ 1 ] = (duplicatedNUBSpline.getCtrlPt(indexOfStartCtrlPt + k)).getY(); // Y coordinate of ctrl pt
					P[ k ][ 2 ] = (duplicatedNUBSpline.getCtrlPt(indexOfStartCtrlPt + k)).getZ(); // Z coordinate of ctrl pt
				}
				rg_Matrix coefficient = R_p * M_p * P; // Coefficient of power basis
													// (p + 1) by 3 matrix
				for(k = 0;k < p + 1;k++)
				{
					coeffOfXCoordinate[j - p][ k ] = coefficient[ k ][ 0 ];
					coeffOfYCoordinate[j - p][ k ] = coefficient[ k ][ 1 ];
					coeffOfZCoordinate[j - p][ k ] = coefficient[ k ][ 2 ];
				}
				delete[] insertingKnotValues1;
				delete[] insertingKnotValues2;
			}
		}
		delete[] distinctKnotValue;
	}
}

void rg_NUBSplineCurve3D::makePiecewisePowerFormPolynomialOfBSplineBasisUsingDE(rg_REAL*** & piecewisePolynomialsInPowerFormOfBasis)
{
	rg_INT degree = getOrder() - 1;
	rg_INT numOfCtrlPts = getNumOfCtrlPts() - 1;
	rg_INT m = numOfCtrlPts + degree + 1;                          // no. of knotvector - 1
	rg_INT numOfNonzeroLengthKnotspan = getNumOfNonZeroLengthKnotSpan();

	piecewisePolynomialsInPowerFormOfBasis = new rg_REAL** [numOfNonzeroLengthKnotspan];
	rg_INT i = 0;
	for(i = 0;i < numOfNonzeroLengthKnotspan;i++)
		piecewisePolynomialsInPowerFormOfBasis[ i ] = new rg_REAL* [degree + 1];

	rg_INT*** allPossiblePath = new rg_INT** [degree + 1]; // all possible paths from (p + 1) basic graphs
	rg_INT* numOfAllPossiblePathsInEachGraph = new rg_INT[degree + 1];

	for(i = 0;i < degree + 1;i++)
	{
		allPossiblePath[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, degree - i);
		numOfAllPossiblePathsInEachGraph[ i ] = rg_MathFunc::combination(degree, i);
	}
			
	// generating the piecewise polynomial curve

	rg_INT indexOfKnotSpanOfNonzeroLength = 0;//, indexOfPolyCurveCoeff = 0;

	for(i = degree;i <= m - degree - 1 && indexOfKnotSpanOfNonzeroLength < numOfNonzeroLengthKnotspan;i++)
	// loop for interior nonzero length knot spans
	{
		rg_INT indexOfTruncatedBasisInEachKnotSpan = 0;
		for(rg_INT j = i - degree;j <= i;j++)
		// 
		{
			piecewisePolynomialsInPowerFormOfBasis[indexOfKnotSpanOfNonzeroLength][indexOfTruncatedBasisInEachKnotSpan] = getDistributionPolynomialsInOneGraphRevised(i,/* j,*/ allPossiblePath[j - i + degree], numOfAllPossiblePathsInEachGraph[j - i + degree]);
			indexOfTruncatedBasisInEachKnotSpan++;
		}
		indexOfKnotSpanOfNonzeroLength++;
	}
	//------------------------------------------
	delete[] allPossiblePath;
	delete[] numOfAllPossiblePathsInEachGraph;
}

rg_REAL* rg_NUBSplineCurve3D::getPowerFormPolynomialUsingBD(const rg_INDEX& index, const rg_INDEX& knotSpanIndex)
{
	rg_INT degree = getOrder() - 1;
	rg_REAL* powerFormPolynomial = new rg_REAL[degree + 1];

	rg_REAL middleKnot = (knotVector[knotSpanIndex] + knotVector[knotSpanIndex + 1]) / 2;
	powerFormPolynomial[ 0 ] = evaluateBasisFunc(index, middleKnot, degree + 1);

	rg_INT i = 0;
	for(i = 1;i <= degree;i++)
		powerFormPolynomial[ i ] = evaluateBasisFuncDerivative(index, middleKnot, i) / rg_MathFunc::factorial( i );
	
	// Extract the coefficient by translating polynomial
	for(i = 0;i < degree;i++)
		for(rg_INT j = degree - 1;j >= i;j--)
			powerFormPolynomial[ j ] += (powerFormPolynomial[j + 1] * (-middleKnot));

	return powerFormPolynomial;
}

rg_REAL* rg_NUBSplineCurve3D::getPowerFormPolynomialUsingBDRecursion(const rg_INDEX& index, const rg_INDEX& knotSpanIndex)
{
	rg_INT degree = getOrder() - 1;
	rg_REAL* powerFormPolynomial = new rg_REAL[degree + 1];

	rg_REAL middleKnot = (knotVector[knotSpanIndex] + knotVector[knotSpanIndex + 1]) / 2;
	powerFormPolynomial[ 0 ] = evaluateBasisFunc(index, middleKnot, degree + 1);

	rg_INT i = 0;
	for(i = 1;i <= degree;i++)
		powerFormPolynomial[ i ] = evaluateBasisFuncDerivativeUsingRecursion(index, middleKnot, i, degree) / rg_MathFunc::factorial( i );
	
	// Extract the coefficient by translating polynomial
	for(i = 0;i < degree;i++)
		for(rg_INT j = degree - 1;j >= i;j--)
			powerFormPolynomial[ j ] += (powerFormPolynomial[j + 1] * (-middleKnot));

	return powerFormPolynomial;
}

void rg_NUBSplineCurve3D::makePiecewisePowerFormPolynomialOfBSplineBasisUsingBD(rg_REAL*** & piecewisePolynomialsInPowerFormOfBasis)
{
	rg_INT degree = getOrder() - 1;
	rg_INT numOfCtrlPts = getNumOfCtrlPts() - 1;
	rg_INT m = numOfCtrlPts + degree + 1;                          // no. of knotvector - 1
	rg_INT numOfNonzeroLengthKnotspan = getNumOfNonZeroLengthKnotSpan();

	piecewisePolynomialsInPowerFormOfBasis = new rg_REAL** [numOfNonzeroLengthKnotspan];
	rg_INT i = 0;
	for(i = 0;i < numOfNonzeroLengthKnotspan;i++)
		piecewisePolynomialsInPowerFormOfBasis[ i ] = new rg_REAL* [degree + 1];

	// generating the piecewise polynomial curve

	rg_INT indexOfKnotSpanOfNonzeroLength = 0;//, indexOfPolyCurveCoeff = 0;

	for(i = degree;i <= m - degree - 1;i++)
	// loop for interior nonzero length knot spans
	{
		rg_INT indexOfTruncatedBasisInEachKnotSpan = 0;
		for(rg_INT j = i - degree;j <= i && indexOfKnotSpanOfNonzeroLength < numOfNonzeroLengthKnotspan;j++)
		// 
		{
			piecewisePolynomialsInPowerFormOfBasis[indexOfKnotSpanOfNonzeroLength][indexOfTruncatedBasisInEachKnotSpan] = getPowerFormPolynomialUsingBD(j, i);
			indexOfTruncatedBasisInEachKnotSpan++;
		}
		indexOfKnotSpanOfNonzeroLength++;
	}
	//------------------------------------------
}

void rg_NUBSplineCurve3D::makePiecewisePowerFormPolynomialOfBSplineBasisUsingBDRecursion(rg_REAL*** & piecewisePolynomialsInPowerFormOfBasis)
{
	rg_INT degree = getOrder() - 1;
	rg_INT numOfCtrlPts = getNumOfCtrlPts() - 1;
	rg_INT m = numOfCtrlPts + degree + 1;                          // no. of knotvector - 1
	rg_INT numOfNonzeroLengthKnotspan = getNumOfNonZeroLengthKnotSpan();

	piecewisePolynomialsInPowerFormOfBasis = new rg_REAL** [numOfNonzeroLengthKnotspan];
	
	rg_INT i = 0;
	for(i = 0;i < numOfNonzeroLengthKnotspan;i++)
		piecewisePolynomialsInPowerFormOfBasis[ i ] = new rg_REAL* [degree + 1];

	// generating the piecewise polynomial curve

	rg_INT indexOfKnotSpanOfNonzeroLength = 0;//, indexOfPolyCurveCoeff = 0;

	for(i = degree;i <= m - degree - 1;i++)
	// loop for interior nonzero length knot spans
	{
		rg_INT indexOfTruncatedBasisInEachKnotSpan = 0;
		for(rg_INT j = i - degree;j <= i && indexOfKnotSpanOfNonzeroLength < numOfNonzeroLengthKnotspan;j++)
		// 
		{
			piecewisePolynomialsInPowerFormOfBasis[indexOfKnotSpanOfNonzeroLength][indexOfTruncatedBasisInEachKnotSpan] = getPowerFormPolynomialUsingBDRecursion(j, i);
			indexOfTruncatedBasisInEachKnotSpan++;
		}
		indexOfKnotSpanOfNonzeroLength++;
	}
	//------------------------------------------
}
rg_NUBSplineCurve3D rg_NUBSplineCurve3D::evaluateCurveSegment( const rg_REAL& start,
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
    
    rg_NUBSplineCurve3D curve(*this);
    rg_INT multiplicity=0;
    rg_INT insertingTimes=0;

    // insert knot at start
    multiplicity=rg_NUBSplineCurve3D::getKnotMultiplicity(start);
    insertingTimes=order-1-multiplicity;
    if( rg_EQ( start, startParam ) )
    {
        insertingTimes++;
    }
    rg_INT i = 0;
	for( i=0; i < insertingTimes; i++ )
        curve.knotInsertion(start);

    // insert knot at end    
    multiplicity=rg_NUBSplineCurve3D::getKnotMultiplicity(end);
    insertingTimes=order-1-multiplicity;
    if( rg_EQ( end, endParam ) )
    {
        insertingTimes++;
    }
    for( i=0; i < insertingTimes; i++ )
    {
        curve.knotInsertion(end);
    }
    
    // compute knot vector
    rg_INT newOrder=rg_BSplineCurve3D::getOrder();

    rg_INT startKnotIndex=curve.rg_NUBSplineCurve3D::getIndexOfKnotSpan(start);
    rg_INT endKnotIndex=curve.rg_NUBSplineCurve3D::getIndexOfKnot(end);
    rg_INT newNumOfKnots=endKnotIndex-startKnotIndex+2*newOrder-1;
    rg_REAL* newKnots=new rg_REAL[newNumOfKnots];

    rg_INT knotIndex=0;
    for(  i=0; i < newOrder-1; i++ )
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
    }

    rg_NUBSplineCurve3D output;
    output.setOrder(newOrder);
    output.setCtrlPts(newNumOfCtrlPts,newCtrlPts);
    output.setKnotVector(newNumOfKnots,newKnots);

    delete[] newKnots;
    delete[] newCtrlPts;
    delete[] newWeights;

    return output;
    
}
////	Set Functions.---------------------------------------------------------
void rg_NUBSplineCurve3D::setInitialKnotVector()
{
    rg_INT  n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT  order = (rg_INT) rg_BSplineCurve3D::getOrder();
    rg_INT  maxI  = n + order;

    knotVector = new rg_REAL[maxI];

    for (rg_INT i=0; i<maxI; i++)
    {
        if ( i < order )
            knotVector[i] = 0.;
        else if ( i < n )
            knotVector[i] = (i - order + 1.) / (n - order + 1.);
        else
            knotVector[i] = 1.;
    }
}

void rg_NUBSplineCurve3D::setKnotValue( const rg_INT    &kIndex, 
									 const rg_REAL &newKnotValue )
{
    knotVector[kIndex] = newKnotValue;
}

void rg_NUBSplineCurve3D::setKnotVector( const rg_INT &numOfKnot,
                                      const rg_REAL* const newKnotVector )
{
    rg_INT  n     = rg_BSplineCurve3D::getNumOfCtrlPts();
	rg_INT  order = (rg_INT) rg_BSplineCurve3D::getOrder();

    if ( numOfKnot != n+order )
        return;

    if ( knotVector != rg_NULL )
		delete [] knotVector;

    knotVector = new rg_REAL [numOfKnot];
    for (rg_INT i=0; i<n+order; i++)
        knotVector[i] = newKnotVector[i];
}

////	BasisFunction & rg_Point3D Evaluating.--------------------------------------
rg_REAL rg_NUBSplineCurve3D::evaluateBasisFunc( const rg_INT  &index, 
                                          const rg_REAL &param,
                                          const rg_INT  &Order ) const
{
    rg_INT  n     = rg_BSplineCurve3D::getNumOfCtrlPts() - 1;
//    if (Order == -1)
//        Order = (rg_INT) rg_BSplineCurve3D::getOrder();
   
	if (    ( index==0 && rg_EQ(param, knotVector[0]) )
		 || ( index==n && rg_EQ(param, knotVector[n + Order]) ) )
	{
        return 1.0;
	}
    else if (    rg_LT(param, knotVector[index]) 
		      || rg_GE(param, knotVector[index + Order]) )
	{
        return 0.0;
    }
	else 
	{}

    rg_REAL Uleft  = 0.;
	rg_REAL Uright = 0.;
	rg_REAL temp   = 0.;
	rg_REAL saved  = 0.;

    rg_REAL* triN = new rg_REAL [Order];

    rg_INT i = 0;
	for (i=0; i<Order; i++)
    {
        if (    rg_LT(param, knotVector[index + i + 1])
			 && rg_GE(param, knotVector[index + i]) )
		{
            triN[i] = 1.0;
		}
        else 
            triN[i] = 0.0;
    }	

    for (i=1; i<Order; i++)
    {
        if (triN[0] == 0.0) 
            saved = 0.0;
        else 
		{
            saved = ( (param - knotVector[index]) * triN[0] ) 
			        / ( knotVector[index+i] - knotVector[index] );
		}
	
        for (rg_INT j=0; j<(Order - i); j++)
        {
            Uleft  = knotVector[index + j + 1];
            Uright = knotVector[index + j + i + 1];

            if ( rg_EQ(triN[j+1], 0.0) )
            {
                triN[j] = saved;
                saved = 0.0;
            }
            else
            {
                temp    = triN[j+1] / (Uright - Uleft);
                triN[j] = saved + (Uright - param) * temp;
                saved   = (param - Uleft) * temp;
            }
        }
    }
	
    rg_REAL returnValue = triN[0];

    delete [] triN;

    return returnValue;
}

rg_REAL* rg_NUBSplineCurve3D::evaluateMultiBasisFunc( const rg_INDEX     &knotIndex,
                                                const rg_PARAMETER &u,
                                                const rg_ORDER     &order) const
{
    rg_REAL* nonZeroBasis = new rg_REAL[order];

    rg_REAL  tempBasis  = 0.0;
    rg_INDEX basisIndex = 0;

    nonZeroBasis[0] = 1.0;
    //  N            (u)
    //   knotIndex, 1   
    for (rg_INDEX basisOrder=2; basisOrder<=order; basisOrder++)
    {
        for (rg_INDEX j=basisOrder-1; j>=0; j--)
        {
            basisIndex = knotIndex-(basisOrder-j)+1;
            if ( j == (basisOrder-1) )
            {
                nonZeroBasis[j] 
//                tempBasis
                    = (u - knotVector[basisIndex]) * nonZeroBasis[j-1] 
                      / (knotVector[basisIndex + basisOrder - 1] - knotVector[basisIndex]);
            }
            else if ( j != 0 ) // else if ( j != 0) )
            {
                nonZeroBasis[j] 
//                tempBasis
                    = (u - knotVector[basisIndex]) * nonZeroBasis[j-1] 
                      / (knotVector[basisIndex +  basisOrder - 1] - knotVector[basisIndex])
                      + 
                      (knotVector[basisIndex + basisOrder] - u) * nonZeroBasis[j]
                      / (knotVector[basisIndex + basisOrder] - knotVector[basisIndex + 1]);       
            }
            else // if (j == 0)
            {
                nonZeroBasis[j]
//                tempBasis
                    = (knotVector[basisIndex + basisOrder] - u) * nonZeroBasis[0]
                      / (knotVector[basisIndex + basisOrder] - knotVector[basisIndex + 1]);
            }

//            nonZeroBasis[j] = tempBasis;
        }
    }

    return nonZeroBasis;
}

rg_REAL rg_NUBSplineCurve3D::getCoefficientForDerivativeEvaluation(const rg_INDEX& index1, const rg_INDEX& index2, const rg_INDEX& index3) const
{
	rg_REAL coefficient = 0.0;
	rg_ORDER degree = getOrder() - 1;

	if(index1 == 0 && index2 == 0)
		return 1.0;
	else if(index2 == 0)
	{
		rg_REAL denominator = (knotVector[index3 + degree - index1 + 1] - knotVector[index3]);
		if(rg_NZERO(denominator))
			return getCoefficientForDerivativeEvaluation(index1 - 1, 0, index3) / denominator;
		else
			return 0.0;
	}
	else if(index2 == index1)
	{
		rg_REAL denominator = (knotVector[index3 + degree + 1] - knotVector[index3 + index1]);
		if(rg_NZERO(denominator))
			return - getCoefficientForDerivativeEvaluation(index1 - 1, index1 - 1, index3) / denominator;
		else
			return 0.0;
	}
	else
	{
		rg_REAL denominator = (knotVector[index3 + degree + index2 - index1 + 1] - knotVector[index3 + index2]);
		if(rg_NZERO(denominator))
			return (getCoefficientForDerivativeEvaluation(index1 - 1, index2, index3) - getCoefficientForDerivativeEvaluation(index1 - 1, index2 - 1, index3)) / denominator;
		else
			return 0.0;
	}
}

rg_REAL rg_NUBSplineCurve3D::evaluateBasisFuncDerivative(const rg_INDEX& index, const rg_PARAMETER& u, const rg_ORDER& derivativeOrder) const
{
	rg_REAL derivative = 0.0;

	rg_ORDER degree = getOrder() - 1;

	for(rg_INT i = 0;i <= derivativeOrder;i++)
	{
		derivative += (getCoefficientForDerivativeEvaluation(derivativeOrder, i, index) * evaluateBasisFunc(index + i, u, degree - derivativeOrder + 1));
	}

	derivative = (rg_MathFunc::factorial(degree) / rg_MathFunc::factorial(degree - derivativeOrder)) * derivative;

	return derivative;
}

rg_REAL rg_NUBSplineCurve3D::evaluateBasisFuncDerivativeUsingRecursion(const rg_INDEX& index, const rg_PARAMETER& u, const rg_ORDER& derivaOrder, const rg_ORDER& degree) const
{
	rg_REAL derivative = 0.0;

	if(derivaOrder == 0)
		return evaluateBasisFunc(index, u, degree + 1);
	else
	{
		rg_REAL denominator1 = knotVector[index + degree] - knotVector[index];
		rg_REAL denominator2 = knotVector[index + degree + 1] - knotVector[index + 1];

		rg_REAL leftTerm = 0.0;
		rg_REAL rightTerm = 0.0;

		if(rg_NZERO(denominator1))
			leftTerm = evaluateBasisFuncDerivativeUsingRecursion(index, u, derivaOrder - 1, degree - 1) / denominator1;

		if(rg_NZERO(denominator2))
			rightTerm = evaluateBasisFuncDerivativeUsingRecursion(index + 1, u, derivaOrder - 1, degree - 1) / denominator2;

		return degree * (leftTerm - rightTerm);
	}
}

//  April  3 1997 : Modified
////////////////////////////////////////////////////////////////// 
rg_Point3D rg_NUBSplineCurve3D::evaluatePt( const rg_REAL &param )
{
	if( knotVector == rg_NULL )
		setInitialKnotVector();
/*
    rg_INT  n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT  order = (rg_INT) rg_BSplineCurve3D::getOrder();    

    rg_Point3D ptOnCurve;
    rg_REAL  basisValue = 0.;
    for (rg_INT i=0; i<n; i++)
    {
		basisValue = evaluateBasisFunc(i, param, order);
        if ( rg_NZERO(basisValue) )
            ptOnCurve += rg_BSplineCurve3D::getCtrlPt(i) * basisValue;
    }

    return ptOnCurve;
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

        if ( rg_LT(param, knotVector[middle]) )
            last = middle;
        else if ( rg_GT(param, knotVector[middle + 1]) )
            first = middle + 1;
        else
        {
			if ( middle != first )
			{
	            while ( rg_EQ(knotVector[middle], knotVector[middle+1]) )
		            middle++;
	            validKnot = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
			else
			{
	            validKnot = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
            break;
        }
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
    for (rg_INDEX i=validKnot-order+1; i<=validKnot; i++)
        ptOnCurve += rg_BSplineCurve3D::getCtrlPt(i) 
                     * nonZeroBasis[i - validKnot + order -1];  

    delete [] nonZeroBasis;

    return ptOnCurve;

}

rg_Point3D   rg_NUBSplineCurve3D::evaluatePt( const rg_REAL &param ) const
{

/*
    rg_INT  n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT  order = (rg_INT) rg_BSplineCurve3D::getOrder();    

    rg_Point3D ptOnCurve;
    rg_REAL  basisValue = 0.;
    for (rg_INT i=0; i<n; i++)
    {
		basisValue = evaluateBasisFunc(i, param, order);
        if ( rg_NZERO(basisValue) )
            ptOnCurve += rg_BSplineCurve3D::getCtrlPt(i) * basisValue;
    }

    return ptOnCurve;
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

        if ( rg_LT(param, knotVector[middle]) )
            last = middle;
        else if ( rg_GT(param, knotVector[middle + 1]) )
            first = middle + 1;
        else
        {
			if ( middle != first )
			{
	            while ( rg_EQ(knotVector[middle], knotVector[middle+1]) )
		            middle++;
	            validKnot = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
			else
			{
	            validKnot = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
            break;
        }
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
    for (rg_INDEX i=validKnot-order+1; i<=validKnot; i++)
        ptOnCurve += rg_BSplineCurve3D::getCtrlPt(i) 
                     * nonZeroBasis[i - validKnot + order -1];  

    delete [] nonZeroBasis;

    return ptOnCurve;

}

//  April  3 1997 : Modified
////////////////////////////////////////////////////////////////// 

rg_Point3D* rg_NUBSplineCurve3D::evaluatePtsInEvenParameter( const rg_INT &noOfEvaluatedPoint ) const
{
    /*
    if( rg_NUBSplineCurve3D::getKnotVector() == rg_NULL )
        rg_NUBSplineCurve3D::setInitialKnotVector();
    */

    rg_Point3D* evaluatedPoint = new rg_Point3D [noOfEvaluatedPoint];

    rg_REAL increment = 1. / (noOfEvaluatedPoint-1);
    rg_REAL u         = 0.;

    rg_INT  n             = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_ORDER order        = rg_BSplineCurve3D::getOrder();
    rg_INDEX validKnot    = 0;
    rg_REAL* nonZeroBasis = rg_NULL;

    evaluatedPoint[0] = rg_BSplineCurve3D::getCtrlPt(0);
    evaluatedPoint[noOfEvaluatedPoint-1] = rg_BSplineCurve3D::getCtrlPt(n-1);

    u += increment;
    for (rg_INT i=1; i<noOfEvaluatedPoint-1; i++)
    {
        while ( rg_LT(u, knotVector[validKnot])
                || rg_GE(u, knotVector[validKnot + 1]) )          
        {
            validKnot++;
        }

        nonZeroBasis = evaluateMultiBasisFunc( validKnot, u, order );
        for (rg_INDEX j=validKnot-order+1; j<=validKnot; j++)
        {
            evaluatedPoint[i] += rg_BSplineCurve3D::getCtrlPt(j) 
                                 * nonZeroBasis[j - validKnot + order -1];  
        }
        delete [] nonZeroBasis;

        u += increment;
    }

    return evaluatedPoint;
}

rg_Point3D* rg_NUBSplineCurve3D::evaluatePt_Plus( const rg_INT &noOfEvaluatedPoint )
{
/*
	if( knotVector == rg_NULL )
		setInitialKnotVector();

   	rg_INT    order     = (rg_INT) rg_BSplineCurve3D::getOrder();
    rg_REAL increment = 1. / noOfEvaluatedPoint;
    rg_REAL u         = 0.;

    rg_Point3D* evaluatedPoint = new rg_Point3D [noOfEvaluatedPoint];

    rg_INT    knotSpan           = order-1;
    rg_REAL lowerBndOfKnotSpan = knotVector[knotSpan];
    rg_REAL upperBndOfKnotSpan = knotVector[knotSpan + 1]; 

    for (rg_INT i=0; i<(noOfEvaluatedPoint - 1); i++)
    {
        while( rg_LT(u, lowerBndOfKnotSpan)
               || rg_GE(u, upperBndOfKnotSpan) )
        {
            knotSpan++;           
            lowerBndOfKnotSpan = knotVector[knotSpan];
            upperBndOfKnotSpan = knotVector[knotSpan + 1];
        }

        if ( rg_GE(u, lowerBndOfKnotSpan)
             && rg_LT(u, upperBndOfKnotSpan) )
        {
        	rg_REAL x = 0.;
            rg_REAL y = 0.;
            rg_REAL z = 0.;

	        rg_REAL basisValue        = 0.;
            rg_INT    firstInfluencedPt = knotSpan - (order-1);
            
            for (rg_INT j=firstInfluencedPt; j<(firstInfluencedPt+order); j++)
            {
        		basisValue = evaluateBasisFunc(j, u, order);

		        x += rg_BSplineCurve3D::getCtrlPt(j).getX() * basisValue;
        		y += rg_BSplineCurve3D::getCtrlPt(j).getY() * basisValue;
        		z += rg_BSplineCurve3D::getCtrlPt(j).getZ() * basisValue;
            }

            evaluatedPoint[i] = rg_Point3D(x, y, z);
            u += increment;
        }
    }
    
    evaluatedPoint[noOfEvaluatedPoint - 1] = evaluatePt(LASTPARAM);//1. - TOLERANCE);

    return evaluatedPoint;
*/
    return rg_NULL;
}
rg_Polyline2D rg_NUBSplineCurve3D::makePolyline2DConsideringKnots(const rg_INT &numOfPts) const
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








/*
rg_Polyline2D rg_NUBSplineCurve3D::makePolyline2DInEvenParameter(const rg_INT &numOfPts) const
{
    if ( isNull() || numOfPts < 1 )
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

rg_Polyline3D rg_NUBSplineCurve3D::makePolyline3DInEvenParameter(const rg_INT &numOfPts) const
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
*/

////	Derivative rg_Curve & Curvature-------------------------------------------

//  April 7 1997 : Modified.
//      Because of copy constructor and = operator overloading

//*-----------------------------------------------------------------------------
// rg_NUBSplineCurve3D* rg_NUBSplineCurve3D::makeDerivative()
//	
//    DESCRIPTION  
//       This function makes The derivative of non-uniform B-Spline curve 
//       which is a (current order - 1)th order non-uniform B-spline curve.
// 
//    INPUT
//       This function has no argument. But notes that this function needs a 
//       complete rg_NUBSplineCurve3D object.
//
//    OUTPUT
//       The rg_NUBSplineCurve3D object pointer to the derivative of the current 
//       non-uniform B-Spline curve.
//
//*-----------------------------------------------------------------------------
rg_NUBSplineCurve3D rg_NUBSplineCurve3D::makeDerivative() const
{
    rg_INT n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT order = (rg_INT) rg_BSplineCurve3D::getOrder();

    rg_REAL  degreeDividedByKnotDelta = 0.;
    rg_Point3D newControlPoint;

    rg_NUBSplineCurve3D derivative(n-1);

    //  Set the order Of derivative curve.
    derivative.rg_BSplineCurve3D::setOrder(order-1);

    //  Set control polygon of derivative curve.
    rg_INT i = 0;
	for (i=0; i<n-1; i++)
    {
        degreeDividedByKnotDelta = ( order-1 ) 
			                       / ( knotVector[i+order] - knotVector[i+1] );
        newControlPoint = degreeDividedByKnotDelta
			              * ( rg_BSplineCurve3D::getCtrlPt(i+1) 
					          - rg_BSplineCurve3D::getCtrlPt(i) );

        derivative.rg_BSplineCurve3D::setCtrlPt( i, newControlPoint );
    }

    //  Set knot vector of derivative curve.
    derivative.knotVector = new rg_REAL [n + order - 2];
    for (i=0; i<(n + order - 2); i++) 
        derivative.knotVector[i] = knotVector[i+1];

    return derivative;
}

//*-----------------------------------------------------------------------------
// rg_REAL getCurvatureInXY( const rg_REAL &param )
//	
//    DESCRIPTION  
//	     Given parameter, u, evaluate the curvature of non-uniform B-Spline 
//       curve in XY-plane.
// 
//    INPUT
//       param : Specifing a position on non-uniform B-Spline curve.
//               (0<= param <= 1)	
//
//    OUTPUT
//       signed curvature of non-uniform B-Spline curve in XY-plane.
//
//*-----------------------------------------------------------------------------
rg_REAL rg_NUBSplineCurve3D::getCurvatureInXY( const rg_REAL &param )
{
//	rg_BSplineCurve3D firstDeriv  = makeDerivative();
//	rg_BSplineCurve3D secondDeriv = firstDeriv.makeDerivative();
	rg_NUBSplineCurve3D firstDeriv  = makeDerivative();
	rg_NUBSplineCurve3D secondDeriv = firstDeriv.makeDerivative();

	rg_Point3D ptOnFirstDeri  = firstDeriv.evaluatePt(param);
	rg_Point3D ptOnSecondDeri = secondDeriv.evaluatePt(param);

	rg_REAL signedCurvature = 0.;
	rg_REAL xPrime          = ptOnFirstDeri.getX();
	rg_REAL yPrime          = ptOnFirstDeri.getY();
	rg_REAL xTwoPrime       = ptOnSecondDeri.getX();  
	rg_REAL yTwoPrime       = ptOnSecondDeri.getY();  

	signedCurvature = ( xTwoPrime * yPrime - yTwoPrime * xPrime )
		              / pow( (xPrime*xPrime + yPrime*yPrime), 1.5 );
	
	return signedCurvature;
}

//--------------------------------------------------------------------
//	Given parameter,u, solve the curvature in XZ-plane.
rg_REAL rg_NUBSplineCurve3D::getCurvatureInXZ( const rg_REAL &param )
{
	rg_BSplineCurve3D firstDeriv  = makeDerivative();
	rg_BSplineCurve3D secondDeriv = firstDeriv.makeDerivative();

	rg_Point3D ptOnFirstDeri  = firstDeriv.evaluatePt(param);
	rg_Point3D ptOnSecondDeri = secondDeriv.evaluatePt(param);

	rg_REAL signedCurvature = 0.;
	rg_REAL xPrime          = ptOnFirstDeri.getX();
	rg_REAL zPrime          = ptOnFirstDeri.getZ();
	rg_REAL xTwoPrime       = ptOnSecondDeri.getX();  
	rg_REAL zTwoPrime       = ptOnSecondDeri.getZ();  

	signedCurvature = ( xTwoPrime * zPrime - zTwoPrime * xPrime )
		              / pow( (xPrime*xPrime + zPrime*zPrime), 1.5 );
	
	return signedCurvature;
}
//--------------------------------------------------------------------
//	Given parameter,u, solve the curvature in YZ-plane.
rg_REAL rg_NUBSplineCurve3D::getCurvatureInYZ( const rg_REAL &param )
{
	rg_BSplineCurve3D firstDeriv  = makeDerivative();
	rg_BSplineCurve3D secondDeriv = firstDeriv.makeDerivative();

	rg_Point3D ptOnFirstDeri  = firstDeriv.evaluatePt(param);
	rg_Point3D ptOnSecondDeri = secondDeriv.evaluatePt(param);

	rg_REAL signedCurvature = 0.;
	rg_REAL yPrime          = ptOnFirstDeri.getY();
	rg_REAL zPrime          = ptOnFirstDeri.getZ();
	rg_REAL yTwoPrime       = ptOnSecondDeri.getY();  
	rg_REAL zTwoPrime       = ptOnSecondDeri.getZ();  

	signedCurvature = ( yTwoPrime * zPrime - zTwoPrime * yPrime )
		              / pow( (yPrime*yPrime + zPrime*zPrime), 1.5 );
	
	return signedCurvature;
}

////	Fundamental Geometric Algorithm.---------------------------------------


rg_INT rg_NUBSplineCurve3D::findTheCorrespondingKnotSpan(const rg_REAL& insertingKnot)
{
	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT p = getOrder() - 1;
	rg_INT m = n + p + 1;

	if(rg_EQ(insertingKnot, knotVector[n + 1]))
		return n;

	rg_INT low = p;
	rg_INT high = n + 1;
	rg_INT mid = (low + high) / 2;

	while(rg_LT(insertingKnot, knotVector[mid]) || rg_GE(insertingKnot, knotVector[mid + 1]))
	{
		if(insertingKnot < knotVector[mid])
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

//   knot insertion :  March 13 1997
void rg_NUBSplineCurve3D::knotInsertion( const rg_REAL &insertingKnotValue )
{
//    if ( rg_LT(insertingKnot, 0.) || rg_GT(insertingKnot, 1.0) )
    if ( !rg_BTOR(0., insertingKnotValue, 1.0) )
        return;

    if (knotVector == rg_NULL)
        setInitialKnotVector();

    rg_INT n     = rg_BSplineCurve3D::getNumOfCtrlPts();
	rg_INT order = rg_BSplineCurve3D::getOrder();

    //  1.  Find the k-th knot-span, [U_k, U_k+1) lying insertingKnot.
    rg_INT k = 0;
	for (k=order-1; k<n; k++)
    {
        if (    rg_GE( insertingKnotValue, knotVector[k]) 
             && rg_LT( insertingKnotValue, knotVector[k+1]) )
            break;
    }	

    //  2.  Construct a new knot vector
    rg_REAL* newKnotVector = new rg_REAL[n+order+1];

    rg_INT i = 0;
	for (i=0; i<=k; i++)
        newKnotVector[i] = knotVector[i];
    newKnotVector[k+1] = insertingKnotValue;
    for (i=k+2; i<n+order+1 ; i++)
        newKnotVector[i] = knotVector[i-1];

    //  3.  Construct a new control polygon influenced by knot insertion. 
    //
    //  3.1 Store control points that aren't influenced by knot insertion
    //      to array for the new control polygon.
    rg_Point3D* newCtrlPlygn = new rg_Point3D[n+1];
    rg_INT j = 0;
	for (j=0; j<=(k-order+1); j++)
    {
        newCtrlPlygn[j] = rg_BSplineCurve3D::getCtrlPt(j);
    }
    //  3.2 Store (order-1) control points influenced by knot insertion
    //      to array. 
    for (j=(k-order+2); j<=k; j++)
    {
        rg_REAL alpha = ( insertingKnotValue - knotVector[j] )
                     / ( knotVector[j+order-1] - knotVector[j] );

        newCtrlPlygn[j] = alpha*rg_BSplineCurve3D::getCtrlPt(j)
                          + (1-alpha)*rg_BSplineCurve3D::getCtrlPt(j-1);
    }
    for (j=k+1; j<n+1; j++)
    {
        newCtrlPlygn[j] = rg_BSplineCurve3D::getCtrlPt(j-1);
    }


    rg_BSplineCurve3D::setCtrlPts(n+1, newCtrlPlygn);

    delete [] knotVector;

    knotVector = newKnotVector;
}

/*
void rg_NUBSplineCurve3D::multipleKnotInsertion(const rg_REAL& insertingKnot, const rg_INT& NumberOfInsertingMultipleKnots)
{
    if ( !rg_BTOR(0., insertingKnot, 1.0) )
        return;

    if (knotVector == rg_NULL)
        setInitialKnotVector();

	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT p = getOrder() - 1;
	rg_INT m = n + p + 1;
	rg_INT r = NumberOfInsertingMultipleKnots;
	rg_INT s = 0; // multiplicity

    // 1.  Find the knot span that the will-be-inserted knot falls into
	//     and the multiplicity of the corresponding knot span.

    for (rg_INT i = p; i <= n;i++)
    {
        if (    rg_GE( insertingKnot, knotVector[ i ]) 
             && rg_LT( insertingKnot, knotVector[i + 1]) )
            break;
    }

	rg_INT k = i; // the corresponding knot span

	if(k <= n)
	{
		while(i == i + 1 && i <= n)
		{
			i++;
			s++;
		}
	}

	// 2. Load a new knot vector

    rg_REAL* newKnotVector = new rg_REAL[n + p + 1 + r];

    for(i = 0; i <= k;i++)
        newKnotVector[ i ] = knotVector[ i ];

	for(i = k + 1;i <= k + r;i++)
		newKnotVector[ i ] = insertingKnot;

    for(i = k + r + 1;i <= m + r;i++)
        newKnotVector[ i ] = knotVector[i - r];

	// 3. Save unaltered control points

	rg_Point3D* newControlPolygon = new rg_Point3D[n + r];

	for(i = 0;i <= k - p;i++)
		newControlPolygon[ i ] = getCtrlPt( i );

	for(i = k;i <= n;i++)
		newControlPolygon[i + r] = getCtrlPt( i );

	rg_Point3D* tempControlPolygon = new rg_Point3D[p + 1];

	for(i = 0;i <= p - s;i++)
		tempControlPolygon[ i ] = getCtrlPt(k - p + i);		

	// 4. Determine the newly-generated control points during inserting the knot r times.

	rg_INT L = 0;

	for(rg_INT j = 1;j <= r;j++)
	{
		L = k - p + j;
		for(i = 0;i < p - j - s;i++)
		{
			rg_REAL alpha = (insertingKnot - getKnotValue(L + i))/(getKnotValue(i + k + 1) - getKnotValue(L + i));
			tempControlPolygon[ i ] = alpha * tempControlPolygon[i + 1] + (1.0 - alpha) * tempControlPolygon[ i ];
		}
		newControlPolygon[ L ] = tempControlPolygon[ 0 ];
		newControlPolygon[k + r - j] = tempControlPolygon[p - j];		
	}

	// 5. Load remaining control points

	for(i = L + 1;i < k;i++)
		newControlPolygon[ i ] = tempControlPolygon[i - L];

	setControlPolygon(p - s + r - 1, newControlPolygon);

	delete[] tempControlPolygon;

    delete [] knotVector;

    knotVector = newKnotVector;
}
*/

void rg_NUBSplineCurve3D::knotRefinement( rg_REAL* &insertingKnotValues, const rg_INT& numOfinsertingKnotValues)
{
	if(numOfinsertingKnotValues == 0)
		return;

	// find indices a and b such that u_a <= x_i < u_b for all i
	// ,where x_i : each inserting knot value in X

	rg_INT a = findTheCorrespondingKnotSpan(insertingKnotValues[ 0 ]);
	rg_INT b = findTheCorrespondingKnotSpan(insertingKnotValues[numOfinsertingKnotValues - 1]) + 1;
	//b++;

    rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT p = getOrder() - 1;
	rg_INT m = n + p + 1;
	rg_Point3D* ctrlPts = getCtrlPts();

	// construct new knot vector
	// and copy unchanged knot values
	rg_REAL* newKnotVector = new rg_REAL[m + 1 + numOfinsertingKnotValues];
	rg_INT i = 0;
	for(i = 0;i <= a;i++)
		newKnotVector[ i ] = knotVector[ i ];
	for(i = b + p;i <= m;i++)
		newKnotVector[i + numOfinsertingKnotValues] = knotVector[ i ];

	// construct new control points
	// and copy unchanged control points
	rg_Point3D* newCtrlPts = new rg_Point3D[n + 1 + numOfinsertingKnotValues];
	for(i = 0;i <= a - p;i++)
		newCtrlPts[ i ] = ctrlPts[ i ];
	for(i = b - 1;i <= n;i++)
		newCtrlPts[i + numOfinsertingKnotValues] = ctrlPts[ i ];

	// computing the changed control points
	// and inserting new knot values

	i = b + p - 1;
	rg_INT k = b + p + numOfinsertingKnotValues - 1;

	for(rg_INT j = numOfinsertingKnotValues - 1;j >= 0;j--)
	{
		while(rg_LE(insertingKnotValues[ j ], knotVector[ i ]) && i > a)
		{
			newCtrlPts[k - p - 1] = ctrlPts[i - p - 1];
			newKnotVector[ k ] = knotVector[ i ];
			k--;
			i--;
		}
		newCtrlPts[k - p - 1] = newCtrlPts[k - p];
		for(rg_INT l = 1;l <= p;l++)
		{
			rg_REAL alpha = newKnotVector[k + l] - insertingKnotValues[ j ];
			if(rg_ZERO(alpha))
				newCtrlPts[k - p + l - 1] = newCtrlPts[k - p + l];
			else
			{
				alpha = alpha / (newKnotVector[k + l] - knotVector[i - p + l]);
				newCtrlPts[k - p + l - 1] = alpha * newCtrlPts[k - p + l - 1] 
					                      + (1.0 - alpha) * newCtrlPts[k - p + l];
			}
		}
		newKnotVector[ k ] = insertingKnotValues[ j ];
		k--;
	}

	// set new knot vector and control points

    setCtrlPts(n + 1 + numOfinsertingKnotValues, newCtrlPts);
    delete [] knotVector;
    knotVector = newKnotVector;

}

//  This function separates cubic non-uniform B-spline curve 
//      into cubic Bezier curves.
rg_BzCurve3D* rg_NUBSplineCurve3D::separateBSplineIntoBezier() 
{
    if ( rg_BSplineCurve3D::getOrder() != (rg_CUBIC+1) 
        || rg_BSplineCurve3D::getCtrlPts() == rg_NULL )
        return rg_NULL;

    rg_INT n = getNumOfKnotSpan();

    rg_BzCurve3D* BezierCurves = new rg_BzCurve3D[n];
    
    //  Set degree of Bezier curves.
    rg_INT i = 0;
	for (i=0; i<n; i++)
        BezierCurves[i].setDegree(rg_CUBIC);

    //  Set the first control point and the last control point 
    //      of Bezier curves.
    rg_INT orderOfBSplineCurve       = rg_BSplineCurve3D::getOrder();
    rg_INT numOfCtrlPtOfBSplineCurve = rg_BSplineCurve3D::getNumOfCtrlPts();

    BezierCurves[0].setCtrlPt(0, rg_BSplineCurve3D::getCtrlPt(0));
    for (i=0; i<(n-1); i++)
    {
        rg_Point3D ctrlPt = evaluatePt( knotVector[i + orderOfBSplineCurve] );
        BezierCurves[i].setCtrlPt(3, ctrlPt);
        BezierCurves[i+1].setCtrlPt(0, ctrlPt);
    }
    BezierCurves[n-1].setCtrlPt(
        3, 
        rg_BSplineCurve3D::getCtrlPt(numOfCtrlPtOfBSplineCurve-1) );

    //  Set the second control point and the third control point
    //      of the first Bezier curve.
    rg_Point3D sndCtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(1);
    rg_Point3D trdCtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(2);

    rg_Point3D trdCtrlPtOf1stBezierCurve;

    rg_REAL fstCoeff = (knotVector[orderOfBSplineCurve + 1] - knotVector[orderOfBSplineCurve])
                    / (knotVector[orderOfBSplineCurve + 1] - knotVector[orderOfBSplineCurve-1] );
    rg_REAL sndCoeff = ( knotVector[orderOfBSplineCurve] - knotVector[orderOfBSplineCurve-1] )
                    / ( knotVector[orderOfBSplineCurve + 1] - knotVector[orderOfBSplineCurve-1] );

    trdCtrlPtOf1stBezierCurve =   (fstCoeff * sndCtrlPtOfBSplineCurve) 
                                + (sndCoeff * trdCtrlPtOfBSplineCurve);
    BezierCurves[0].setCtrlPt( 1, sndCtrlPtOfBSplineCurve );
    BezierCurves[0].setCtrlPt( 2, trdCtrlPtOf1stBezierCurve );

    //  Set the second control point and the third control point
    //      of the last Bezier curve.
    rg_Point3D n_1CtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(
                                        numOfCtrlPtOfBSplineCurve-2);
    rg_Point3D n_2CtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(
                                        numOfCtrlPtOfBSplineCurve-3);

    rg_Point3D sndCtrlPtOfLastBezierCurve;
    
    fstCoeff = (knotVector[numOfCtrlPtOfBSplineCurve] - knotVector[numOfCtrlPtOfBSplineCurve-1])
               / (knotVector[numOfCtrlPtOfBSplineCurve] - knotVector[numOfCtrlPtOfBSplineCurve-2]); 
    sndCoeff = (knotVector[numOfCtrlPtOfBSplineCurve-1] - knotVector[numOfCtrlPtOfBSplineCurve-2])
               / (knotVector[numOfCtrlPtOfBSplineCurve] - knotVector[numOfCtrlPtOfBSplineCurve-2]); 

    sndCtrlPtOfLastBezierCurve =   (fstCoeff * n_2CtrlPtOfBSplineCurve)
                                 + (sndCoeff * n_1CtrlPtOfBSplineCurve);
    BezierCurves[n-1].setCtrlPt( 1, sndCtrlPtOfLastBezierCurve );
    BezierCurves[n-1].setCtrlPt( 2, n_1CtrlPtOfBSplineCurve );

    //  Set the second control point and the third control point
    //      of the other Bezier curves.
    for (i=1; i<(n-1); i++)
    {
        rg_Point3D curPt  = rg_BSplineCurve3D::getCtrlPt(i+2);
        rg_Point3D prevPt = rg_BSplineCurve3D::getCtrlPt(i+1);

        rg_REAL fstCoeffOf2ndCtrlPt = 0.0;
        rg_REAL sndCoeffOf2ndCtrlPt = 0.0;
        rg_REAL fstCoeffOf3rdCtrlPt = 0.0;
        rg_REAL sndCoeffOf3rdCtrlPt = 0.0;

        fstCoeffOf2ndCtrlPt 
            = (knotVector[orderOfBSplineCurve + i + 1] - knotVector[orderOfBSplineCurve + i - 1])
              / (knotVector[orderOfBSplineCurve + i + 1] - knotVector[orderOfBSplineCurve + i - 2]);
        sndCoeffOf2ndCtrlPt 
            = (knotVector[orderOfBSplineCurve + i - 1] - knotVector[orderOfBSplineCurve + i - 2])
              / (knotVector[orderOfBSplineCurve + i + 1] - knotVector[orderOfBSplineCurve + i - 2]);
        fstCoeffOf3rdCtrlPt 
            = (knotVector[orderOfBSplineCurve + i + 1] - knotVector[orderOfBSplineCurve + i])
              / (knotVector[orderOfBSplineCurve + i + 1] - knotVector[orderOfBSplineCurve + i - 2]);
        sndCoeffOf3rdCtrlPt 
            = (knotVector[orderOfBSplineCurve + i] - knotVector[orderOfBSplineCurve + i - 2])
              / (knotVector[orderOfBSplineCurve + i + 1] - knotVector[orderOfBSplineCurve + i - 2]);

        BezierCurves[i].setCtrlPt( 
            1, 
            (fstCoeffOf2ndCtrlPt*prevPt) + (sndCoeffOf2ndCtrlPt*curPt) );
        BezierCurves[i].setCtrlPt( 
            2, 
            (fstCoeffOf3rdCtrlPt*prevPt) + (sndCoeffOf3rdCtrlPt*curPt) );
    }

    return BezierCurves;
}

rg_BzCurve3D* rg_NUBSplineCurve3D::decomposeCurveIntoBezierSegment() const
{

    //////////////////////////////////////////////////////////////////////////////////////////
	rg_NUBSplineCurve3D duplicatedNUBSpline = (*this);

	//duplicatedNUBSpline.reparameterizationKnotVector();

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
		while( (j < m - p) && rg_EQ(knotVector[j], knotVector[j+1]))
		{
			s++;
			j++;
		}
		multiplicity[i][ 0 ] = knotVector[j];  // replicated knot value in the interior knot span
		multiplicity[i][ 1 ] = s;               // multiplicity of replicated knot value in the interior knot span
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
			duplicatedNUBSpline.knotInsertion(multiplicity[i][0]);
			j++;
		}
	}

	//////////////////////////////////////////////////////

	//////////////////////////////////////////////////////
	// provide the storage for Bezier segments
	//////////////////////////////////////////////////////

    rg_BzCurve3D* BezierCurveList = new rg_BzCurve3D[numberOfKnotSpan];

	//////////////////////////////////////////////////////
	// assign the respective newly-generated control points to each Bezier curve segment.
	//////////////////////////////////////////////////////

    for (i=0; i<numberOfKnotSpan; i++)
    {
		BezierCurveList[ i ].setDegree( p );
		for(rg_INT j = 0;j < p + 1;j++)
		{
			BezierCurveList[ i ].setCtrlPt(j, duplicatedNUBSpline.getCtrlPt(i * p + j));
		}
	}


	//////////////////////////////////////////////////////

	for(i = 0;i < numberOfInteriorKnot;i++)
		delete [] multiplicity[ i ];

	delete[] multiplicity;

	return BezierCurveList;
/*
    rg_NUBSplineCurve3D duplicatedNUBSpline = (*this);

	//duplicatedNUBSpline.reparameterizationKnotVector();

	// record the index, the knot value 
	// and multiplicity whose multiplicity is more than zero
	// in the "multiplicity[][]"
	//////////////////////////////////////////////////////

	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT p = getOrder() - 1;
	rg_INT m = n + p + 1;

	rg_INT numberOfKnotSpan     = getNumOfNonZeroLengthKnotSpan(); // the number of nonzero length knot span
	
    rg_REAL** multiplicity = new rg_REAL* [numberOfKnotSpan];

	for(rg_INT i = 0;i < numberOfKnotSpan;i++)
		multiplicity[ i ] = new rg_REAL[ 2 ];

	i = p + 1;
	while(i <= m - p - 1)
	{
		rg_INT s = 0;
		rg_INT j = i + 1;

		while((j <= m - p) && rg_EQ(knotVector[ i ], knotVector[ j ]))
		{
			s++;
			j++;
		}
		//multiplicity[i - p - 1][ 0 ] = i;				  // Index of replicated knot value in the interior knot span
		multiplicity[i - p - 1][ 0 ] = knotVector[ i ];  // replicated knot value in the interior knot span
		multiplicity[i - p - 1][ 1 ] = s;				  // multiplicity of replicated knot value in the interior knot span
		i = j;
	}
	//////////////////////////////////////////////////////

	
	// knot insertion is applied into the interior knot whose multiplicity is less than p(degree)
	//////////////////////////////////////////////////////
	i = p + 1;
	while(i <= m - p - 1)
	{
		rg_INT j = 0;
		// knot insertion is repeated until the total multiplicity of any one interior knot is p!!
		if(multiplicity[i - p - 1][ 1 ] < p)
		{
			while(j < (p - multiplicity[i - p - 1][ 1 ] - 1))
			{
				duplicatedNUBSpline.knotInsertion(multiplicity[i - p - 1][ 0 ]);
				j++;
			}
		}
		i++;
	}

	//////////////////////////////////////////////////////

	//////////////////////////////////////////////////////
	// provide the storage for Bezier segments
	//////////////////////////////////////////////////////

    rg_BzCurve3D* BezierCurveList = new rg_BzCurve3D[numberOfKnotSpan];

	//////////////////////////////////////////////////////
	// assign the respective newly-generated control points to each Bezier curve segment.
	//////////////////////////////////////////////////////

	i = 0;
	while(i < numberOfKnotSpan)
	{
		BezierCurveList[ i ].setDegree( p );
		for(rg_INT j = 0;j < p + 1;j++)
		{
			BezierCurveList[ i ].setCtrlPt(j, duplicatedNUBSpline.getCtrlPt(i * p + j));
		}
		i++;
	}


	//////////////////////////////////////////////////////

	for(i = 0;i < numberOfKnotSpan;i++)
		delete [] multiplicity[ i ];

	delete[] multiplicity;

	return BezierCurveList;
*/
}

rg_BzCurve3D* rg_NUBSplineCurve3D::decomposeCurveIntoBezierSegmentUsingKnotRefinement() const
{

    //////////////////////////////////////////////////////////////////////////////////////////
	rg_NUBSplineCurve3D duplicatedNUBSpline = (*this);

	//duplicatedNUBSpline.reparameterizationKnotVector();

	// record the index, the knot value 
	// and multiplicity whose multiplicity is more than zero
	// in the "multiplicity[][]"
	//////////////////////////////////////////////////////

	rg_INT n = getNumOfCtrlPts() - 1;
	rg_INT p = getOrder() - 1;
	rg_INT m = n + p + 1;

	rg_INT numberOfKnotSpan     = getNumOfNonZeroLengthKnotSpan(); // the number of nonzero length knot span
    rg_INT numberOfInteriorKnot = numberOfKnotSpan-1;

	//rg_INT numOfMultiplicity = 0;
	rg_INT numOfInsertingKnotValues = 0;
    rg_REAL** multiplicity = new rg_REAL* [numberOfInteriorKnot];

	rg_INT i = 0;
	for(i = 0;i < numberOfInteriorKnot;i++)
		multiplicity[ i ] = new rg_REAL[ 2 ];

    i     = 0;      //  i traces the array multiplicity.
    rg_INT j = p+1;    //  j traces knotVector.
	while( i < numberOfInteriorKnot )
	{
		rg_INT s = 1;
		while( (j < m - p) && rg_EQ(knotVector[j], knotVector[j+1]))
		{
			s++;
			//numOfMultiplicity++;
			j++;
		}
		multiplicity[i][ 0 ] = knotVector[j];  // replicated knot value in the interior knot span
		multiplicity[i][ 1 ] = s;               // multiplicity of replicated knot value in the interior knot span
		numOfInsertingKnotValues += (p - s);
		i++;
        j++;
	}
	//////////////////////////////////////////////////////
	//rg_INT numOfInsertingKnotValues = numOfMultiplicity + numberOfInteriorKnot;
    rg_REAL* insertingKnotValues = new rg_REAL [numOfInsertingKnotValues];

	i = 0;
	rg_INT indexOfinsertingKnotValues = 0;
	while((i < numberOfInteriorKnot) && (indexOfinsertingKnotValues < numOfInsertingKnotValues))
	{
		rg_INT s = 1;
		while(s <= p - multiplicity[ i ][ 1 ])
		{
			insertingKnotValues[indexOfinsertingKnotValues] = multiplicity[ i ][ 0 ];
			indexOfinsertingKnotValues++;
			s++;
		}
		i++;
	}
	
	// knot refinement is performed in case of non-zero interior knot
	//////////////////////////////////////////////////////

	if(numberOfInteriorKnot > 0)
		duplicatedNUBSpline.knotRefinement(insertingKnotValues, numOfInsertingKnotValues);

	//////////////////////////////////////////////////////

	//////////////////////////////////////////////////////
	// provide the storage for Bezier segments
	//////////////////////////////////////////////////////

    rg_BzCurve3D* BezierCurveList = new rg_BzCurve3D[numberOfKnotSpan];

	//////////////////////////////////////////////////////
	// assign the respective newly-generated control points to each Bezier curve segment.
	//////////////////////////////////////////////////////

    for (i=0; i<numberOfKnotSpan; i++)
    {
		BezierCurveList[ i ].setDegree( p );
		for(rg_INT j = 0;j < p + 1;j++)
		{
			BezierCurveList[ i ].setCtrlPt(j, duplicatedNUBSpline.getCtrlPt(i * p + j));
		}
	}


	//////////////////////////////////////////////////////

	for(i = 0;i < numberOfInteriorKnot;i++)
		delete [] multiplicity[ i ];

	delete[] multiplicity;

	return BezierCurveList;
}


//  This function pulls out a cubic Bezier curve with responding to 
//  a single knot span.
rg_BzCurve3D* rg_NUBSplineCurve3D::pullOutCubicBezierForKnotSpan(
                                const rg_REAL &prevKnot, 
                                const rg_REAL &nextKnot) 
{
    
    rg_INT numOfKnotSpan = getNumOfKnotSpan();
    rg_INT order         = rg_BSplineCurve3D::getOrder();
    rg_INT numOfCtrlPt   = rg_BSplineCurve3D::getNumOfCtrlPts();

    if ( order != (rg_CUBIC+1) )
        return rg_NULL;

    //  determine the index of knot span.
    rg_INT indexOfKnotSpan;
    for( rg_INDEX i=(order-1) ; i<numOfCtrlPt; i++)
    {
        if ( rg_EQ(prevKnot, knotVector[i]) )
        {
            if ( rg_EQ(nextKnot, knotVector[i+1]) )
                indexOfKnotSpan = i;
        }
    }

    rg_BzCurve3D* BezierCurve = new rg_BzCurve3D(rg_CUBIC);

    if ( indexOfKnotSpan == 0 )
    {
        BezierCurve->setCtrlPt(0, rg_BSplineCurve3D::getCtrlPt(0));
        BezierCurve->setCtrlPt(3, evaluatePt(knotVector[order]));

        rg_Point3D sndCtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(1);
        rg_Point3D trdCtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(2);

        rg_Point3D trdCtrlPtOf1stBezierCurve;

        rg_REAL fstCoeff = (knotVector[order + 1] - knotVector[order])
                        / (knotVector[order + 1] - knotVector[order - 1] );
        rg_REAL sndCoeff = ( knotVector[order] - knotVector[order - 1] )
                        / ( knotVector[order + 1] - knotVector[order - 1] );

        trdCtrlPtOf1stBezierCurve =   (fstCoeff * sndCtrlPtOfBSplineCurve) 
                                    + (sndCoeff * trdCtrlPtOfBSplineCurve);
        BezierCurve->setCtrlPt( 1, sndCtrlPtOfBSplineCurve );
        BezierCurve->setCtrlPt( 2, trdCtrlPtOf1stBezierCurve );

        return BezierCurve;
    }
    else if ( indexOfKnotSpan == (numOfKnotSpan-1) )
    {
        BezierCurve->setCtrlPt(0, evaluatePt(knotVector[order + numOfKnotSpan - 1]));
        BezierCurve->setCtrlPt(3, rg_BSplineCurve3D::getCtrlPt(numOfCtrlPt));

        rg_Point3D n_1CtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(
                                            numOfCtrlPt-2);
        rg_Point3D n_2CtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(
                                            numOfCtrlPt-3);

        rg_Point3D sndCtrlPtOfLastBezierCurve;
    
        rg_REAL fstCoeff = (knotVector[numOfCtrlPt] - knotVector[numOfCtrlPt - 1])
                        / (knotVector[numOfCtrlPt] - knotVector[numOfCtrlPt - 2]); 
        rg_REAL sndCoeff = (knotVector[numOfCtrlPt - 1] - knotVector[numOfCtrlPt - 2])
                        / (knotVector[numOfCtrlPt] - knotVector[numOfCtrlPt - 2]); 

        sndCtrlPtOfLastBezierCurve =   (fstCoeff * n_2CtrlPtOfBSplineCurve)
                                     + (sndCoeff * n_1CtrlPtOfBSplineCurve);
        BezierCurve->setCtrlPt( 1, sndCtrlPtOfLastBezierCurve );
        BezierCurve->setCtrlPt( 2, n_1CtrlPtOfBSplineCurve );

        return BezierCurve;
    }
    else
    {
        BezierCurve->setCtrlPt(0, evaluatePt(knotVector[indexOfKnotSpan + order - 1]));
        BezierCurve->setCtrlPt(3, evaluatePt(knotVector[indexOfKnotSpan + order]));

        rg_Point3D curPt  = rg_BSplineCurve3D::getCtrlPt(indexOfKnotSpan + 2);
        rg_Point3D prevPt = rg_BSplineCurve3D::getCtrlPt(indexOfKnotSpan + 1);

        rg_REAL fstCoeffOf2ndCtrlPt = 0.0;
        rg_REAL sndCoeffOf2ndCtrlPt = 0.0;
        rg_REAL fstCoeffOf3rdCtrlPt = 0.0;
        rg_REAL sndCoeffOf3rdCtrlPt = 0.0;

        fstCoeffOf2ndCtrlPt 
            = (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 1])
              / (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 2]);
        sndCoeffOf2ndCtrlPt 
            = (knotVector[order + indexOfKnotSpan - 1] - knotVector[order + indexOfKnotSpan - 2])
              / (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 2]);
        fstCoeffOf3rdCtrlPt 
            = (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan])
              / (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 2]);
        sndCoeffOf3rdCtrlPt 
            = (knotVector[order + indexOfKnotSpan] - knotVector[order + indexOfKnotSpan - 2])
              / (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 2]);

        BezierCurve->setCtrlPt( 
            1, 
            (fstCoeffOf2ndCtrlPt*prevPt) + (sndCoeffOf2ndCtrlPt*curPt) );
        BezierCurve->setCtrlPt( 
            2, 
            (fstCoeffOf3rdCtrlPt*prevPt) + (sndCoeffOf3rdCtrlPt*curPt) );

        return BezierCurve;
    }
}

rg_BzCurve3D* rg_NUBSplineCurve3D::pullOutCubicBezierForKnotSpan(const rg_INDEX &indexOfKnotSpan) 
{    
    rg_INT numOfKnotSpan = getNumOfKnotSpan();
    rg_INT order         = rg_BSplineCurve3D::getOrder();
    rg_INT numOfCtrlPt   = rg_BSplineCurve3D::getNumOfCtrlPts();

    if ( order != (rg_CUBIC+1) )
        return rg_NULL;

    rg_BzCurve3D* BezierCurve = new rg_BzCurve3D(rg_CUBIC);

    if ( indexOfKnotSpan == 0 )
    {
        BezierCurve->setCtrlPt(0, rg_BSplineCurve3D::getCtrlPt(0));
        BezierCurve->setCtrlPt(3, evaluatePt(knotVector[order]));

        rg_Point3D sndCtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(1);
        rg_Point3D trdCtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(2);

        rg_Point3D trdCtrlPtOf1stBezierCurve;

        rg_REAL fstCoeff = (knotVector[order + 1] - knotVector[order])
                        / (knotVector[order + 1] - knotVector[order - 1] );
        rg_REAL sndCoeff = ( knotVector[order] - knotVector[order - 1] )
                        / ( knotVector[order + 1] - knotVector[order - 1] );

        trdCtrlPtOf1stBezierCurve =   (fstCoeff * sndCtrlPtOfBSplineCurve) 
                                    + (sndCoeff * trdCtrlPtOfBSplineCurve);
        BezierCurve->setCtrlPt( 1, sndCtrlPtOfBSplineCurve );
        BezierCurve->setCtrlPt( 2, trdCtrlPtOf1stBezierCurve );

        return BezierCurve;
    }
    else if ( indexOfKnotSpan == (numOfKnotSpan-1) )
    {
        BezierCurve->setCtrlPt(0, evaluatePt(knotVector[order + numOfKnotSpan - 1]));
        BezierCurve->setCtrlPt(3, rg_BSplineCurve3D::getCtrlPt(numOfCtrlPt));

        rg_Point3D n_1CtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(
                                            numOfCtrlPt-2);
        rg_Point3D n_2CtrlPtOfBSplineCurve = rg_BSplineCurve3D::getCtrlPt(
                                            numOfCtrlPt-3);

        rg_Point3D sndCtrlPtOfLastBezierCurve;
    
        rg_REAL fstCoeff = (knotVector[numOfCtrlPt] - knotVector[numOfCtrlPt - 1])
                        / (knotVector[numOfCtrlPt] - knotVector[numOfCtrlPt - 2]); 
        rg_REAL sndCoeff = (knotVector[numOfCtrlPt - 1] - knotVector[numOfCtrlPt - 2])
                        / (knotVector[numOfCtrlPt] - knotVector[numOfCtrlPt - 2]); 

        sndCtrlPtOfLastBezierCurve =   (fstCoeff * n_2CtrlPtOfBSplineCurve)
                                     + (sndCoeff * n_1CtrlPtOfBSplineCurve);
        BezierCurve->setCtrlPt( 1, sndCtrlPtOfLastBezierCurve );
        BezierCurve->setCtrlPt( 2, n_1CtrlPtOfBSplineCurve );

        return BezierCurve;
    }
    else
    {
        BezierCurve->setCtrlPt(0, evaluatePt(knotVector[indexOfKnotSpan + order - 1]));
        BezierCurve->setCtrlPt(3, evaluatePt(knotVector[indexOfKnotSpan + order]));

        rg_Point3D curPt  = rg_BSplineCurve3D::getCtrlPt(indexOfKnotSpan + 2);
        rg_Point3D prevPt = rg_BSplineCurve3D::getCtrlPt(indexOfKnotSpan + 1);

        rg_REAL fstCoeffOf2ndCtrlPt = 0.0;
        rg_REAL sndCoeffOf2ndCtrlPt = 0.0;
        rg_REAL fstCoeffOf3rdCtrlPt = 0.0;
        rg_REAL sndCoeffOf3rdCtrlPt = 0.0;

        fstCoeffOf2ndCtrlPt 
            = (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 1])
              / (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 2]);
        sndCoeffOf2ndCtrlPt 
            = (knotVector[order + indexOfKnotSpan - 1] - knotVector[order + indexOfKnotSpan - 2])
              / (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 2]);
        fstCoeffOf3rdCtrlPt 
            = (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan])
              / (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 2]);
        sndCoeffOf3rdCtrlPt 
            = (knotVector[order + indexOfKnotSpan] - knotVector[order + indexOfKnotSpan - 2])
              / (knotVector[order + indexOfKnotSpan + 1] - knotVector[order + indexOfKnotSpan - 2]);

        BezierCurve->setCtrlPt( 
            1, 
            (fstCoeffOf2ndCtrlPt*prevPt) + (sndCoeffOf2ndCtrlPt*curPt) );
        BezierCurve->setCtrlPt( 
            2, 
            (fstCoeffOf3rdCtrlPt*prevPt) + (sndCoeffOf3rdCtrlPt*curPt) );

        return BezierCurve;
    }
}

rg_FLAG rg_NUBSplineCurve3D::makeCompositeCurveWithC0(const rg_NUBSplineCurve3D &curve1,
                                                const rg_NUBSplineCurve3D &curve2)
{
    rg_INT n1 = curve1.getNumOfCtrlPts();
    rg_INT n2 = curve2.getNumOfCtrlPts();

    rg_NUBSplineCurve3D firstCurve;
    rg_NUBSplineCurve3D secondCurve;
    //  check whether two curves have the same end point or not.
    if ( curve1.rg_BSplineCurve3D::getCtrlPt(n1-1) == curve2.rg_BSplineCurve3D::getCtrlPt(0) 
         && curve1.rg_BSplineCurve3D::getOrder() == curve2.rg_BSplineCurve3D::getOrder() )
    {
        firstCurve = curve1;
        secondCurve = curve2;
    }
    else if ( curve1.rg_BSplineCurve3D::getCtrlPt(0) == curve2.rg_BSplineCurve3D::getCtrlPt(n2-1) 
              && curve1.rg_BSplineCurve3D::getOrder() == curve2.rg_BSplineCurve3D::getOrder() )
    {
        firstCurve = curve2;
        secondCurve = curve1;
    }
    else
    {
        return rg_FALSE;
    }

    n1 = firstCurve.getNumOfCtrlPts();
    n2 = secondCurve.getNumOfCtrlPts();

    ///////////////////////////////////////////////////////////////// 
    rg_Point3D* ctrlPtsOfCompositeCurve = new rg_Point3D [n1 + n2 - 1];
    rg_INT i = 0;
	for (i=0; i<n1; i++)
        ctrlPtsOfCompositeCurve[i] = firstCurve.getCtrlPt(i);

    for (i=n1; i<(n1+n2-1); i++)
        ctrlPtsOfCompositeCurve[i] = secondCurve.getCtrlPt(i-n1+1);

    rg_INT order = firstCurve.rg_BSplineCurve3D::getOrder();

    //  Set the control polygon and the order of new composite curve.
    rg_BSplineCurve3D::setCtrlPts( n1+n2-1, ctrlPtsOfCompositeCurve );
    rg_BSplineCurve3D::setOrder( order );

    //////////////////////////////////////////////////////////////////
    rg_REAL* knotOfCompositeCurve = new rg_REAL[ n1+n2-1 + order ];
    for (i=0; i<n1+order-1; i++)
        knotOfCompositeCurve[i] = firstCurve.knotVector[i];
    for (i=n1+order-1; i<n1+n2-1+order; i++)
        knotOfCompositeCurve[i] = secondCurve.knotVector[i-(n1+order-1) + order] + 1.0;

    //  Normalize the knot vector of new composite curve.
    for (i=0; i<n1+n2-1+order; i++)
        knotOfCompositeCurve[i] = knotOfCompositeCurve[i]/2.0; 

    setKnotVector( n1+n2-1 + order, knotOfCompositeCurve );

    delete [] ctrlPtsOfCompositeCurve;
    delete [] knotOfCompositeCurve;

    return rg_TRUE;
}

void rg_NUBSplineCurve3D::formLine( const rg_Point3D& start,
                                 const rg_Point3D& end, const rg_INT &order )
{
    if ( order < 2 )
        return;

    rg_INT newNumOfPts   = order;
    rg_INT newNumOfKnots = order + newNumOfPts;

    rg_Point3D* newCtrlPts=new rg_Point3D[newNumOfPts];
    rg_REAL*  newKnots  =new rg_REAL[newNumOfKnots];

    newCtrlPts[0]             = start;
    newCtrlPts[newNumOfPts-1] = end;
    rg_INT i = 0;
	for (i=1; i<newNumOfPts-1; i++)
    {
        rg_REAL t = (rg_REAL)i/(order-1);
        newCtrlPts[i] = (1-t)*start + t*end;
    }

    for( i=0; i < order; i++ )
    {
        newKnots[i]       = 0.0;
        newKnots[order+i] = 1.0;
    }

    rg_BSplineCurve3D::setOrder(order);
    rg_BSplineCurve3D::setCtrlPts(newNumOfPts,newCtrlPts);
    setKnotVector(newNumOfPts+order,newKnots);

    delete[] newCtrlPts;
    delete[] newKnots;
}


rg_sListByPtr* rg_NUBSplineCurve3D::intersectOfCubicBSplineAndPlane(const rg_Plane3D &plane) const
{

    if ( rg_BSplineCurve3D::getOrder() != (rg_CUBIC+1) 
        || rg_BSplineCurve3D::getNumOfCtrlPts() == 0 )
        return rg_NULL;

    rg_sListByPtr* intersectPoint = new rg_sListByPtr;

    rg_INT numOfBzCurve            = getNumOfNonZeroLengthKnotSpan();
    rg_BzCurve3D* cubicBezierCurve = decomposeCurveIntoBezierSegment();
    for (rg_INDEX i=0; i<numOfBzCurve; i++)
    {
        rg_sListByPtr* intersectBezierAndPlane 
            = cubicBezierCurve[i].intersectOfCubicBezierAndPlane(plane);
        if ( intersectBezierAndPlane->getNumberOfEntity() < 1 )
        {
            delete intersectBezierAndPlane;
            continue;
        }

        if ( intersectPoint->getNumberOfEntity() >= 1 )
        {
            rg_Point3D* prevPt = (rg_Point3D*)intersectPoint->getLastEntity();
            rg_Point3D* nextPt = (rg_Point3D*)intersectBezierAndPlane->getFirstEntity();

            if ( *prevPt == *nextPt )
            {
                intersectBezierAndPlane->killOneNodeNentity( intersectBezierAndPlane->getFirstEntity() );
            }
        }
        intersectPoint->concatenate(intersectBezierAndPlane);

/*
        //  # of rg_MathFunc::root is 3.
        rg_ComplexNumber* rg_MathFunc::root = cubicBezierCurve[i].intersectOfCubicAndPlane(plane);

        for (rg_INDEX j=0; j<rg_CUBIC; j++)
        {
            if ( rg_ZERO( rg_MathFunc::root[j].getImaginaryNumber() ) ) 
            {
                rg_REAL realRoot = rg_MathFunc::root[j].getRealNumber();
                if ( rg_BTOR(0.0, realRoot, 1.0) )
                {
                    rg_Point3D* pt = new rg_Point3D(cubicBezierCurve[i].evaluatePt(realRoot));
                    if ( intersectPoint->getNumberOfEntity() == 0 )
                        intersectPoint->append(pt);
                    else if ( *pt != *(rg_Point3D*)intersectPoint->getLastEntity() )
                        intersectPoint->append(pt);
                    else;
                }
            }
        }
        delete [] rg_MathFunc::root;
*/
    }

    delete [] cubicBezierCurve;

    return intersectPoint;                    
}

rg_dList<rg_Point3D> rg_NUBSplineCurve3D::intersectWithPlaneForCubic(const rg_Plane3D& plane) const
{
    rg_dList<rg_Point3D> output;

    if ( rg_BSplineCurve3D::getOrder() != (rg_CUBIC+1) 
        || rg_BSplineCurve3D::getNumOfCtrlPts() == 0 )
    {
        return output;
    }
    

    rg_INT numOfBzCurve            = getNumOfNonZeroLengthKnotSpan();
    rg_BzCurve3D* cubicBezierCurve = decomposeCurveIntoBezierSegment();
    for (rg_INDEX i=0; i<numOfBzCurve; i++)
    {
        rg_dList<rg_Point3D>  intersectionsWithOnePlane
                      = cubicBezierCurve[i].intersectWithPlaneForCubic(plane);

        if (   intersectionsWithOnePlane.getSize() > 0 )
        {
            if ( output.getSize() > 0 )
            {
                rg_Point3D prevPt = output.getLastEntity();
                rg_Point3D nextPt = intersectionsWithOnePlane.getFirstEntity();

                if ( prevPt == nextPt )
                {
                    intersectionsWithOnePlane.killHead();
                }
            }
            output.appendTail(intersectionsWithOnePlane);
        }


    }
    delete [] cubicBezierCurve;

    return output;                    
}


//  3D curve에서 inflection point는 의미없음.
//  이 함수는 3D 공간상의 평면 곡선에 대한 함수임. 
rg_REAL* rg_NUBSplineCurve3D::inflectionPointByHodograph(rg_INT& numIPts) const
{
	rg_Point3D normal;
	//  1. 평면곡선인지 아닌지를 검사한다. 
	//     평면곡선일 경우 normal vector를 얻는다.
	if ( !rg_BSplineCurve3D::isPlanarCurve(normal) )  {
		numIPts = 0;
		return rg_NULL;
	}

	rg_NUBSplineCurve3D copyCurve(*this);
	//  2. normal vector를 이용하여 평면곡선을 xy평면에 놓이도록 rotation한다.
	rg_TMatrix3D rMat;
	rMat.rotate(normal, rg_Point3D(0.0, 0.0, 1.0));

	rg_INT      degree     = copyCurve.getOrder() - 1;
	rg_Point3D* ctrlPts    = copyCurve.getCtrlPts();
	rg_INT      numCtrlPts = copyCurve.getNumOfCtrlPts();
	rg_INT i = 0;
	for(i=0; i<numCtrlPts; i++)  {
		ctrlPts[i] = rMat*ctrlPts[i];
	}
	copyCurve.setCtrlPts(numCtrlPts, ctrlPts);
	delete [] ctrlPts;

	//  3. 3D Bezier curve segment로 나눈다.
	rg_INT numOf3DBzSegs = copyCurve.getNumOfNonZeroLengthKnotSpan(); 
	rg_BzCurve3D* bzSegs = copyCurve.decomposeCurveIntoBezierSegment();

	//  4. 3D Bezier curve segment들을 2D Bezier curve로 변환한다.
	rg_BzCurve2D* bzSegsInXY = new rg_BzCurve2D[ numOf3DBzSegs ];
	for (i=0; i<numOf3DBzSegs; i++)  {
		bzSegsInXY[i].setDegree( degree );

		rg_Point2D* bzCtrlPts = new rg_Point2D[degree+1];
		for(rg_INT j=0; j<=degree; j++)  {
			rg_Point3D   pt = bzSegs[i].getCtrlPt( j );
			bzCtrlPts[j] = pt.evaluatePt2D();
		}
		bzSegsInXY[i].setCtrlPts(bzCtrlPts);
		delete [] bzCtrlPts;
	}

	//  5. calculate inflection points of 2D bezier curves.
	rg_dList<rg_REAL> inflection;
	rg_REAL* distinctKnot = getDistinctKnotValues();
	for (i=0; i<numOf3DBzSegs; i++)  {
		rg_INT numInf;
		rg_REAL* infOfBezierCurve = bzSegsInXY[i].inflectionPointByHodograph( numInf );
		for (rg_INT j=0; j<numInf; j++)  {
			rg_REAL param;
			//  must be reparameterized.
			param = distinctKnot[i]
				    +infOfBezierCurve[j]*(distinctKnot[i+1]-distinctKnot[i]);
			if ( i == 0 )  {
				inflection.add( param );
			}
			else if ( i>0 && rg_NE(inflection.getLastEntity(), param) )  {
				inflection.add( param );
			}
			else{}
		}
	}
	numIPts = inflection.getSize();

	return inflection.getArray();
}
void  rg_NUBSplineCurve3D::appendWithC0(      rg_NUBSplineCurve3D next,
                                        const rg_REAL connectedKnot )
{
    rg_INT newOrder=this->getOrder();

    if ( this->isNull() == rg_TRUE )
    {
        *this=next;
        return;
    }

    if (  next.getOrder() != newOrder ) 
    {
        return;
    }
    rg_Point3D preHead=getStartPoint();
    rg_Point3D preTail=getEndPoint();
    rg_Point3D nextHead=next.getStartPoint();
    rg_Point3D nextTail=next.getEndPoint();
    if (  preTail != nextHead )
    {
        if ( preTail == nextTail )
        {
            next.reverseTrace();
        }
        else if ( preHead == nextHead )
        {
            this->reverseTrace();
        }
        else if ( preHead == nextTail )
        {
            rg_NUBSplineCurve3D temp=*this;
            *this=next;
            next=temp;
        }
        else
        {
            return;
        }
    }

    // compute control points
    const rg_INT prevNumOfCtrlPts=this->getNumOfCtrlPts();
    const rg_INT nextNumOfCtrlPts=next.getNumOfCtrlPts();
    const rg_INT newNumOfCtrlPts= prevNumOfCtrlPts+ nextNumOfCtrlPts-1;
    rg_Point3D*   newCtrlPts= new rg_Point3D[newNumOfCtrlPts];

    rg_INT ctrlPtsIndex=0;
    rg_INT i = 0;
	for( i=0; i < prevNumOfCtrlPts-1; i++)
    {
        newCtrlPts[ctrlPtsIndex]=this->getCtrlPt(i);
        ctrlPtsIndex++;
    }
    for(  i=0; i < nextNumOfCtrlPts; i++)
    {
        newCtrlPts[ctrlPtsIndex]=next.getCtrlPt(i);
        ctrlPtsIndex++;
    }

    // compute Knot vector
    const rg_REAL prevScale=connectedKnot;
    const rg_REAL nextScale=1.0-connectedKnot;
    const rg_INT  prevNumOfKnots=this->getNumberOfKnotValues();
    const rg_INT  nextNumOfKnots=next.getNumberOfKnotValues();
    const rg_INT  newNumOfKnots= newOrder+newNumOfCtrlPts;
    rg_REAL*      newKnotVector=new rg_REAL[newNumOfKnots];
    rg_INT  knotIndex=0;

    for ( i=0 ; i <prevNumOfKnots-1; i++ )
    {
        newKnotVector[knotIndex]=prevScale*(this->getKnotValue(i));
        knotIndex++;
    }

    for( i=newOrder; i < nextNumOfKnots; i++ )
    {
        newKnotVector[knotIndex]=connectedKnot+nextScale*(next.getKnotValue(i));
        knotIndex++;
    }

    // set this curve
    this->setCtrlPts(newNumOfCtrlPts,newCtrlPts);
    this->setKnotVector(newNumOfKnots, newKnotVector);

    delete[] newCtrlPts;
    delete[] newKnotVector;

}

void   rg_NUBSplineCurve3D::reverseTrace()
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
    }
    const rg_INT numOfKnots=getNumberOfKnotValues();
    rg_INT halfNumOfKnots=numOfKnots/2;
    for( i=0; i < halfNumOfKnots; i++ )
    {
        const rg_INT before=i;
        const rg_INT after=numOfCtrlPts-1-i;
        const rg_REAL temp=knotVector[before];
        knotVector[before]=1.0-knotVector[after];
        knotVector[after]=1.0-temp;
    }

}
    


////	rg_Curve Generating function----------------------------------------------
rg_REAL** rg_NUBSplineCurve3D::makeMatrix( rg_Point3D*      b, 
	    						     const rg_REAL*   u, 
									 const rg_INT     &L, 
									 const rg_Point3D &m0, 
									 const rg_Point3D &mL )
{
	///////////////////////////
	rg_SignArray<rg_REAL> delta(-1, L+2);

	delta[-1] = delta[L] = 0.0;
	rg_INT i = 0;
	for (i=0; i<=L-1; i++)	delta[i] = u[i+1] - u[i];

	////////////////////////////
	rg_Point3D *r = new rg_Point3D[L+1];
	r[0] = b[0] + delta[0] / 3.0 * m0;
	r[L] = b[L] - delta[L-1] / 3.0 * mL;

	for (i=1; i<=L-1; i++) 
		r[i] = (delta[i-1]+delta[i])*b[i];

	for (i=0; i<=L; i++)  
		b[i] = r[i];

	////////////////////////////
	rg_SignArray<rg_REAL> alpha(1, L-1);
	rg_SignArray<rg_REAL> beta (1, L-1);
	rg_SignArray<rg_REAL> gamma(1, L-1);

	for(i=1; i<=L-1; i++) 
	{
		alpha[i] = pow(delta[i], 2) / (delta[i-2] + delta[i-1] + delta[i]);
		
		beta[i]  =   delta[i] 
			         * (delta[i-2] + delta[i-1]) 
					 / (delta[i-2] + delta[i-1] + delta[i])
				   + delta[i-1] 
				     * (delta[i] + delta[i+1]) 
					 / (delta[i-1] + delta[i] + delta[i+1]);
		gamma[i] = pow(delta[i-1], 2) / (delta[i-1] + delta[i] + delta[i+1]);
	}

	////////////////////////////
	rg_INT    s=0;
    rg_INT    bandSize=1;
	rg_REAL **m = new rg_REAL* [2 * bandSize + 1];

	for (i=0; i<=(2 * bandSize); i++) 
	{
		if( i <= bandSize ) 
			s = bandSize-i;
		else 
			s = i-bandSize;

		m[i] = new rg_REAL [ L + 1 - s];
	}

///////////////////////////////////////
	for (i=0; i<=L-2; i++) 
		m[0][i] = alpha[i+1];

	m[0][L-1] = 0;

	for (i=1; i<=L-1; i++)
		m[1][i] = beta[i];

	m[1][0] = m[1][L] = 1;

	for (i=1; i<=L-1; i++) 
		m[2][i] = gamma[i];

	m[2][0] = 0;
//////////////////////////////////////
/*
	for( i=0; i < L-1; i++ )
	{
		m[0][i]=alpha[i+1];
		m[1][i+1]=beta[i];
		m[2][i+1]=gamma[i];
	}

	m[0][L-1] = 0.0;
	m[1][0] = m[1][L] = 1.0;
	m[2][0] = 0.0;
*/
	delete [] r;

	return m;
}

void rg_NUBSplineCurve3D::interpolateWithEndCondition( const rg_INT &noOfData, 
									                   rg_Point3D*     fittingData,
                                                       rg_Point3D      startTngnt,
                                                       rg_Point3D      endTngnt,
                                                       rg_REAL*      parametrization)
{
    const rg_INT order=4;
    rg_REAL* param = parametrization;
    if (param == rg_NULL)
//      param = rg_GeoFunc::rg_GeoFunc::chordLength(noOfData, fittingData);
        param = rg_GeoFunc::centripetal(noOfData, fittingData);

    if ( rg_ZERO(startTngnt.getX()) && rg_ZERO(startTngnt.getY()) && rg_ZERO(startTngnt.getZ()) )
        startTngnt = fittingData[1] - fittingData[0];
//	    rg_Point3D m0(0, (fittingData[1].getY() - fittingData[0].getY()), 0);

    if ( rg_ZERO(endTngnt.getX()) && rg_ZERO(endTngnt.getY()) && rg_ZERO(endTngnt.getZ()) )
        endTngnt = fittingData[noOfData-1] - fittingData[noOfData-2];
//	rg_Point3D mL(fittingData[noOfData-1] - fittingData[noOfData-2]);

	rg_Point3D *data = new rg_Point3D[noOfData];
	rg_INT i = 0;
	for (i = 0; i < noOfData; i++) 
		data[i] = fittingData[i];
	
	rg_Point3D *d = new rg_Point3D[noOfData];

	rg_REAL **mat = makeMatrix(data, param,  noOfData-1, startTngnt, endTngnt);

    rg_INT bandSize=1;
	rg_BandedMatrix<rg_Point3D> bm2(noOfData, bandSize);

	bm2.findVariable(mat, d, data);

    //  set order
    rg_BSplineCurve3D::setOrder(order);

    //  set control polygon
	rg_BSplineCurve3D::setNumOfCtrlPts(noOfData + 2);

	rg_BSplineCurve3D::setCtrlPt(0, fittingData[0]);
	for (i=0; i<noOfData; i++)
		rg_BSplineCurve3D::setCtrlPt(i+1, d[i]);
	rg_BSplineCurve3D::setCtrlPt(noOfData+1, fittingData[noOfData-1]);
	
    //  set knot
	rg_INT degree = (rg_INT)(rg_BSplineCurve3D::getOrder() - 1); 
	rg_INT noOfKnot = noOfData + 2*degree;
	knotVector = new rg_REAL [noOfKnot];

	for (i = 0; i<degree; i++)
	{
		knotVector[i] = 0.;
		knotVector[noOfKnot - 1 - i] = 1.;
	}

	for (i = degree; i<(noOfKnot-degree); i++)
		knotVector[i] = param[i - degree];
	
	delete [] data;
	delete [] d;

	for(i=0; i<=2*bandSize; i++)
		delete [] mat[i];
	delete [] mat;

    if ( parametrization == rg_NULL )
        delete [] param;
}

void rg_NUBSplineCurve3D::interpolate(const rg_INT& numPts,
								 	  rg_Point3D*   pts,
									  const rg_INT& order,
									  rg_REAL*      parametrization)
{
    // In the case order > numPts
    // , following procedure is done
    //    1. order = numPts 
    //    2. interpolate pts ( Always, this curve can be converted Bezier curve )
    //    3. get BzCurve which is identical to this curve
    //    4. raise degree of BzCurve by numPts-order
    //    5. convert BzCurve to this Curve in NUBSpline curve form

    if ( numPts < 2 )
    {
        removeAll();
        return;
    }

    rg_INT raisingTimes = order - numPts;
    rg_INT tOrder = ( raisingTimes > 0  ) ? numPts: order;
    //     set order
    rg_BSplineCurve3D::setOrder(tOrder);

    //     set control points.
	rg_BSplineCurve3D::setNumOfCtrlPts( numPts );


	//  1. Determine the parametrization.    
	rg_REAL* param = parametrization;
    if (param == rg_NULL)  
    {
//      param = rg_GeoFunc::rg_GeoFunc::chordLength(numPts, pts);
        param = rg_GeoFunc::centripetal(numPts, pts);
	}

	//  2. Determine the knot vector.
	if ( knotVector )  {
		delete [] knotVector;
	}
	knotVector = new rg_REAL[ numPts + tOrder ];

	rg_INT i, j;
	for (i=0; i<tOrder; i++)  {
		knotVector[i] = 0.0;
		knotVector[numPts+tOrder-1 - i] = 1.0;
	}
	//  knot vector obtained by averaging parameters.
	for (i=1; i<numPts-tOrder+1; i++)  {
		rg_REAL sumParam = 0.0;
		for(j=i; j<i+tOrder-1; j++)  {
			sumParam += param[j];
		}
		knotVector[i+tOrder-1] = (sumParam)/(tOrder-1);
	}

	//  3. Solve the system of linear equation.
	//           cofMat*sol_ = b_
	//     sol_x : x's coordinate value of control points.
	//     sol_y : y's coordinate value of control points.
	//     sol_z : z's coordinate value of control points.

	//  3.1  Construct the system of linear equation.
	rg_Matrix cofMat(numPts, numPts);

	rg_Matrix b_x(numPts, 1);
	rg_Matrix b_y(numPts, 1);
	rg_Matrix b_z(numPts, 1);

	rg_Matrix sol_x(numPts, 1);
	rg_Matrix sol_y(numPts, 1);
	rg_Matrix sol_z(numPts, 1);

	for (i=0; i<numPts; i++)  {
		for (j=0; j<numPts; j++)  {
			rg_REAL N_j_order = evaluateBasisFunc(j, param[i], tOrder);
			cofMat.setElement(i, j, N_j_order );
//			cofMat.setElement(i, j, evaluateBasisFunc(j, param[i], order) );
		}
		b_x.setElement(i, 0, pts[i].getX() );
		b_y.setElement(i, 0, pts[i].getY() );
		b_z.setElement(i, 0, pts[i].getZ() );
	}
/*
	ofstream fout("cnsMatrix.dat");
	fout << cofMat << endl;
	fout << b_x << endl;
	fout << b_y << endl;
	fout << b_z << endl;
*/
	//  3.2  Solve the system of linear equation.
	cofMat.gaussianElimination( sol_x, b_x );
	cofMat.gaussianElimination( sol_y, b_y );
	cofMat.gaussianElimination( sol_z, b_z );

//	rg_Matrix cofMat_1; 
//	cofMat_1 = cofMat.inverse();
//	sol_x = cofMat_1*b_x;
//	sol_y = cofMat_1*b_y;
//	sol_z = cofMat_1*b_z;
/*
	fout << sol_x << endl;
	fout << sol_y << endl;
	fout << sol_z << endl;
*/
	//  4. Set non-uniform B-spline curve.
	for (i=0; i<numPts; i++)  {
		rg_Point3D ctrlPt(sol_x[i][0], sol_y[i][0], sol_z[i][0]);
		rg_BSplineCurve3D::setCtrlPt(i, ctrlPt );
	}
    
	if ( parametrization == rg_NULL )  {
        delete [] param;
	}

    if ( raisingTimes > 0 )
    {
        rg_BzCurve3D* curve=this->decomposeCurveIntoBezierSegment();
        curve->raiseDegree(raisingTimes);
        this->bzCurvesToNUBSplineCurve(1,curve);
        delete[] curve;
    }
}

void     rg_NUBSplineCurve3D::interpolateWithFiltering(const rg_INT& numPts,
													   rg_Point3D*   pts,
													   const rg_INT& order ,
													   rg_REAL*      parametrization)
{
	rg_INT numOfFilteredPts=numPts;
	rg_Point3D* filteredPts=makeFilteredPassingPts(numOfFilteredPts, pts);
	interpolate(numOfFilteredPts,filteredPts, order, parametrization);
	delete[] filteredPts;
}

rg_INT  rg_NUBSplineCurve3D::getNumberOfKnotValues() const
{
    rg_INT order=this->getOrder();
    rg_INT numOfCtrlPts=this->getNumOfCtrlPts();

    return order+numOfCtrlPts;
}

rg_INT  rg_NUBSplineCurve3D::getIndexOfKnot(const rg_REAL& tKnotValue) const
{
    
    rg_INT numOfKnots=this->getNumberOfKnotValues();

    if ( !rg_BTOR(knotVector[0],tKnotValue,knotVector[numOfKnots-1]) )
    {
        return -1;
    }

    rg_INT output=-1;
    for(rg_INT i =0; i < numOfKnots; i++ )
    {
        if ( rg_EQ(knotVector[i],tKnotValue) )
        {
            output=i;
            break;
        }
    }

    return output;
}



rg_INT  rg_NUBSplineCurve3D::getKnotMultiplicity(const rg_REAL& tKnotValue) const
{
    rg_INT numOfKnots=this->getNumberOfKnotValues();
    rg_INT knotIndex=getIndexOfKnot(tKnotValue);
   
    if (    knotIndex <= -1 
         || knotIndex >= numOfKnots )
    {
        return 0;
    }
    
    rg_INT multiplicity=1;
    rg_INT i=knotIndex-1;

    // Add multiplicity corresponding knots before knot[knotIndex]
    while(    i >  -1 
           && rg_EQ(knotVector[i],tKnotValue) )
    {
        multiplicity++;
        i--;
    }

    // Add multiplicity corresponding knots after knot[knotIndex]
    i=knotIndex+1;

    while(   i < numOfKnots
          && rg_EQ( knotVector[i], tKnotValue) )
    {
        multiplicity++;
        i++;
    }

    return multiplicity;
}

// return the index of knotspan containing tParameter

rg_INT  rg_NUBSplineCurve3D::getIndexOfKnotSpan( const rg_REAL& tParameter) const
{
    rg_INT numOfKnots=this->getNumberOfKnotValues();
    rg_INT order=this->getOrder();

    rg_INT indexOfKnotSpan=0;
    // For the case but first knot <= tParameter <= kno_tVector
    if ( !rg_BTOR(knotVector[0], tParameter, knotVector[numOfKnots -1] ) )
    {
        return -1;
    }

    // For the special case that tParameter = last knot value
    if ( rg_EQ( knotVector[numOfKnots-1] , tParameter ) )
    {
        indexOfKnotSpan=numOfKnots-2;

        while( rg_LT(knotVector[indexOfKnotSpan], tParameter ) )
        {
            indexOfKnotSpan--;
        }

        return indexOfKnotSpan;

    }
    
    // For the general case
    indexOfKnotSpan=0;

    while ( !rg_LT( tParameter, knotVector[indexOfKnotSpan+1] ) )
    {
        indexOfKnotSpan++;
    }

    return indexOfKnotSpan;
}

rg_REAL  rg_NUBSplineCurve3D::findParameterOfNearestPt(const rg_Point3D& pt,
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
		rg_NUBSplineCurve3D currentCurve;
        if ( rg_EQ(u,prevU) )
        {
            continue;
        }
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


rg_REAL  rg_NUBSplineCurve3D::findParameterOfNearestPt(const rg_Point3D& pt,
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
    rg_NUBSplineCurve3D derivative[2];
    derivative[0]=this->makeDerivative();
    derivative[1]=derivative[0].makeDerivative();

    rg_Point3D seedPt=this->evaluatePt(seed);
    rg_Point3D seedTangent=(derivative[0].evaluatePt(seed)).getUnitVector();
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
    
    rg_REAL startParameter=rg_NUBSplineCurve3D::getStartParameter();
    rg_REAL endParameter=rg_NUBSplineCurve3D::getEndParameter();
    rg_FLAG doEndIteration=rg_FALSE;
    do
    {
        prevU=u;

        rg_Point3D c[3];
        c[0]=this->evaluatePt(u);
        c[1]=derivative[0].evaluatePt(u);
        c[2]=derivative[1].evaluatePt(u);


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
////    Operator Overloading.--------------------------------------------------
rg_NUBSplineCurve3D& rg_NUBSplineCurve3D::operator =(const rg_NUBSplineCurve3D &curve)
{
    if (this == &curve )
        return *this;

    rg_INT n     = curve.rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT order = curve.rg_BSplineCurve3D::getOrder();

    rg_Curve::setID(curve.getID());
    rg_Curve::setPlanarity(curve.getPlanarity());

    rg_BSplineCurve3D::setOrder( order );

	rg_Point3D* ctrlPtsOnCurve=curve.rg_BSplineCurve3D::getCtrlPts();
	rg_BSplineCurve3D::setCtrlPts( n,ctrlPtsOnCurve );
	delete[] ctrlPtsOnCurve;



    if (knotVector != rg_NULL)
        delete [] knotVector;

   	if ( curve.isNull() )
	{
		knotVector=NULL;

	}
	else
	{
		knotVector = new rg_REAL [n+order];
		for (rg_INT i=0; i<n+order; i++)
			knotVector[i] = curve.knotVector[i];
	}

    return *this;
}

////	Conversion between power spline & b-spline form.
void rg_NUBSplineCurve3D::powerSplineToNUBSplineCurve( const rg_DEGREE& dgr,
							                        const rg_INT& numOfSegment, 
                                                    rg_REAL** paramValues,
                                                    const rg_Matrix* powerCoeff )
{
	rg_BzCurve3D* curveSegment = new rg_BzCurve3D[numOfSegment];
	
	for( rg_INDEX i=0; i < numOfSegment; i++ )
	{
		rg_REAL interval[2];
		interval[0] = paramValues[i][0];
		interval[1] = paramValues[i][1];

		curveSegment[i].powerToBezierCurve( dgr, interval, powerCoeff[i] );
	}

	bzCurvesToNUBSplineCurve( numOfSegment, curveSegment );

	delete[] curveSegment;
}

void rg_NUBSplineCurve3D::bzCurvesToNUBSplineCurve( const rg_INT& numOfSegment,
												 const rg_BzCurve3D* bzCurves )
{
	if( numOfSegment <= 0 )
		return;

	if( numOfSegment >1 )
	{
		for( rg_INDEX i=0; i < (numOfSegment - 1); i++ )
		{
			if( bzCurves[i].getDegree() != bzCurves[i+1].getDegree() )
				return;
		}
	}
	
	// set order & control points.
	rg_ORDER ord = bzCurves[0].getDegree() + 1;	
	rg_INT numOfCtrlPt = ord * numOfSegment;

	rg_BSplineCurve3D::setOrder( ord );
	rg_BSplineCurve3D::setNumOfCtrlPts( numOfCtrlPt );
	
	for( rg_INDEX i=0; i < numOfSegment; i++ )
	{
		for( rg_INDEX j=0; j < ord; j++ )
			rg_BSplineCurve3D::setCtrlPt( (i*ord)+j, bzCurves[i].getCtrlPt( j ) );
	}

	// initialize knot.
	setInitialKnotVector();

	// set remaining knots.
	if( numOfSegment > 1 )
	{
		for( rg_INDEX i=0; i < (numOfSegment-1); i++ )
		{
			rg_REAL newKnot = (rg_REAL)(i+1)/(rg_REAL)numOfSegment;
		
			for( rg_INT j=0; j < ord; j++ )
				setKnotValue( (i+1)*ord+j, newKnot );
		}
	}
}

////	Reparameterization
////	This function reparameterize knots from 0 to 1.
void rg_NUBSplineCurve3D::reparameterizationKnotVector()
{
    rg_INT n     = rg_BSplineCurve3D::getNumOfCtrlPts();
    rg_INT order = rg_BSplineCurve3D::getOrder();
	
	rg_INT numOfKnots = n + order;

	if( rg_EQ( knotVector[0], 0. ) && rg_EQ( knotVector[numOfKnots-1], 1. ) )
		return;

	rg_REAL firstKnot = knotVector[0];

	rg_INT i = 0;
	for( i=0; i < numOfKnots; i++ )
		knotVector[i] = knotVector[i] - firstKnot;
	
	rg_REAL interval = knotVector[numOfKnots - 1] - knotVector[0];

	for( i=0; i < numOfKnots; i++ )
		knotVector[i] = knotVector[i] / interval;
}

/*
////  Parameterization.--------------------------------------------------------
////  functions below moved to 'rg_CurveSurfaceFunc.h & cpp' by Lee, Dong-Gyou 17 Mar 1998 

rg_REAL* rg_GeoFunc::rg_GeoFunc::chordLength(const rg_INT &n, const rg_Point3D* const ptsPassedThrough)
{
	rg_REAL d = 0.0;
	for(rg_INT i=1; i<n; i++) 
        d += ptsPassedThrough[i].distance(ptsPassedThrough[i-1]);

	///////////////////////////
	rg_REAL *u = new rg_REAL[n];

	u[0] = 0;
	u[n-1] = 1;

	for(i=1; i<n-1; i++)
        u[i] = u[i-1] 
               + ptsPassedThrough[i].distance(ptsPassedThrough[i-1])/d;

	return u;
}

rg_REAL* rg_GeoFunc::centripetal(const rg_INT &n, const rg_Point3D* const ptsPassedThrough)
{
	rg_REAL d = 0;
	rg_INT i = 0;
	for (i=1; i<n; i++) 
        d += sqrt(ptsPassedThrough[i].distance(ptsPassedThrough[i-1]));

	///////////////////////////
	rg_REAL *u = new rg_REAL[n];

	u[0] = 0;
	u[n-1] = 1;

	for (i=1; i<n-1; i++)
        u[i] = u[i-1] 
               + sqrt(ptsPassedThrough[i].distance(ptsPassedThrough[i-1]))/d;

	return u;
}
*/

void rg_NUBSplineCurve3D::removeAll()
{
    if ( knotVector != rg_NULL )
    {
		delete [] knotVector;
        knotVector=rg_NULL;
    }
    rg_BSplineCurve3D::removeAll();

}



