#include <math.h>
//#include <fstream>
#include <time.h>

#include "rg_BSplineIntersector.h"
#include "rg_BzCurve3D.h"
#include "rg_PolynomialWithBound.h"
#include "rg_IntersectFunc.h"
#include "sortFunc.h"
#include "rg_RelativeOp.h"

const rg_REAL pi = 4.0*atan(1.0);


rg_BSplineIntersector::rg_BSplineIntersector()
{
}

rg_BSplineIntersector::~rg_BSplineIntersector()
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Decompose the given NUBSpline curve into Bezier curves and 
// then apply Cocktail Algorithm to each pair of Bezier rg_Curve !!(filtering process not yet implemented)

rg_dList<rg_Point2D> rg_BSplineIntersector::intersectBSplineCurveVsBSplineCurveUsingCurveDecomposition(rg_NUBSplineCurve3D& curve_s,
																							  rg_NUBSplineCurve3D& curve_t,
																							  rg_dList<rg_REAL*> & seedParamList4TwoCurve,
																							  rg_REAL& time)
{
	// splitting into Bezier rg_Curve
	// apply Our Algorithm to each pair of Bezier rg_Curve !!(filtering process not yet implemented)

	clock_t StartTime, EndTime;
    StartTime = clock();

	//rg_BzCurve3D* BezierCurveList_s = curve_s.decomposeCurveIntoBezierSegment();
	//rg_BzCurve3D* BezierCurveList_t = curve_t.decomposeCurveIntoBezierSegment();

	rg_BzCurve3D* BezierCurveList_s = curve_s.decomposeCurveIntoBezierSegmentUsingKnotRefinement();
	rg_BzCurve3D* BezierCurveList_t = curve_t.decomposeCurveIntoBezierSegmentUsingKnotRefinement();

	rg_INT numberOfBzCurvesInCurve_s = curve_s.getNumOfNonZeroLengthKnotSpan();
	rg_INT numberOfBzCurvesInCurve_t = curve_t.getNumOfNonZeroLengthKnotSpan();

	rg_BzIntersector BzCurveIntersector;

    rg_dList<rg_Point2D> intersectPointList;
		
	for(rg_INT i = 0;i < numberOfBzCurvesInCurve_s;i++)
	{
		rg_dList<rg_REAL*> tSeedParam4TwoCurve;
		rg_dList<rg_Point2D> tIntersectPointList;

		for(rg_INT j = 0;j < numberOfBzCurvesInCurve_t;j++)
		{
			tIntersectPointList = BzCurveIntersector.intersectBzCurveVsBzCurve( BezierCurveList_s[ i ].evaluateBzCurve2D(),
																				BezierCurveList_t[ j ].evaluateBzCurve2D(),
																				tSeedParam4TwoCurve    );
			intersectPointList.append(tIntersectPointList);
			seedParamList4TwoCurve.append(tSeedParam4TwoCurve);
			tIntersectPointList.removeAll();
			tSeedParam4TwoCurve.removeAll();
		}

	}

	EndTime = clock();
	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;

	return intersectPointList;
}


rg_dList<rg_Point3D> rg_BSplineIntersector::findCharPointUsingCurveDecomposition(rg_NUBSplineCurve3D &curve)
{
	rg_dList<rg_Point3D> charPoint;

	rg_BzCurve3D* BezierCurveList = curve.decomposeCurveIntoBezierSegment();

	rg_INT numberOfBzCurvesInCurve = curve.getNumOfNonZeroLengthKnotSpan();

	for(rg_INT i = 0;i < numberOfBzCurvesInCurve;i++)
	{
		rg_BzIntersector intersector;

		rg_dList<rg_REAL> charParam = intersector.findCharParam(BezierCurveList[ i ].evaluateBzCurve2D());

		charParam.reset();

		do{
			charPoint.add(BezierCurveList[ i ].evaluatePt(charParam.getEntity()));
			charParam.setCurrentNext();
		}while(!charParam.isHead());
	}

	return charPoint;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



rg_dList<rg_Point2D> rg_BSplineIntersector::intersectBSplineCurveVsBSplineCurve(rg_NUBSplineCurve3D& curve_s, 
																	   rg_NUBSplineCurve3D& curve_t,
																	   rg_dList<rg_REAL*> & seedParam4TwoCurve, rg_REAL & time)
{
    rg_dList<rg_Point2D> intersectPointList;

	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime, EndTime;
    StartTime = clock();


	//for(rg_INT i=0; i<20; i++)
	//{

        intersectPointList.removeAll();
  
        rg_dList<rg_REAL>    param_s = findCharParam(curve_s);
        rg_dList<rg_REAL>    param_t = findCharParam(curve_t);

        rg_dList<rg_RQBzCurve2D> rqBzCurveList_s    = approximateRQBzCurves(curve_s, param_s);
        rg_dList<rg_RQBzCurve2D> rqBzCurveList_t    = approximateRQBzCurves(curve_t, param_t);

	    seedParam4TwoCurve = makeSeed(rqBzCurveList_s, param_s, 
								      rqBzCurveList_t, param_t);

//    rg_dList<rg_Point2D> intersectPointList;
        if(seedParam4TwoCurve.getSize()==0) return intersectPointList;

        seedParam4TwoCurve.reset();

	    do
	    {
		    rg_REAL *parameter = seedParam4TwoCurve.getEntity();
		    rg_Point2D intersectPt = iterationWithSeed(curve_s, parameter[0], 
			                                        curve_t, parameter[1]);
		    
		    intersectPointList.add(intersectPt);
		    seedParam4TwoCurve.setCurrentNext();
	    } while(!seedParam4TwoCurve.isHead());

	//}

	EndTime = clock();
	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
//    ofstream fout("timeOut.dat", ios::app);
//    fout << time << endl;

	return intersectPointList;
}


rg_dList<rg_Point2D> rg_BSplineIntersector::intersectBSplineCurveVsBSplineCurveUsingDirectPolynomialGeneration(rg_NUBSplineCurve3D& curve_s, 
																	                                  rg_NUBSplineCurve3D& curve_t,
																	                                  rg_dList<rg_REAL*> & seedParam4TwoCurve, rg_REAL & time)
{
    rg_dList<rg_Point2D> intersectPointList;

	// This part is to find a seed
	// record the time to find the seed using implicitization approach
	clock_t StartTime, EndTime;
    StartTime = clock();


	//for(rg_INT i=0; i<20; i++)
	//{

        intersectPointList.removeAll();
  
        rg_dList<rg_REAL>    param_s = findCharParamUsingDirectPolynomialGeneration(curve_s);
        rg_dList<rg_REAL>    param_t = findCharParamUsingDirectPolynomialGeneration(curve_t);

        rg_dList<rg_RQBzCurve2D> rqBzCurveList_s    = approximateRQBzCurves(curve_s, param_s);
        rg_dList<rg_RQBzCurve2D> rqBzCurveList_t    = approximateRQBzCurves(curve_t, param_t);

	    seedParam4TwoCurve = makeSeed(rqBzCurveList_s, param_s, 
								      rqBzCurveList_t, param_t);

//    rg_dList<rg_Point2D> intersectPointList;
        if(seedParam4TwoCurve.getSize()==0) return intersectPointList;

        seedParam4TwoCurve.reset();

	    do
	    {
		    rg_REAL *parameter = seedParam4TwoCurve.getEntity();
		    rg_Point2D intersectPt = iterationWithSeed(curve_s, parameter[0], 
			                                        curve_t, parameter[1]);
		    
		    intersectPointList.add(intersectPt);
		    seedParam4TwoCurve.setCurrentNext();
	    } while(!seedParam4TwoCurve.isHead());

	//}

	EndTime = clock();
	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
//    ofstream fout("timeOut.dat", ios::app);
//    fout << time << endl;

	return intersectPointList;
}


// Has not yet considered the case where curve has the conjugate tagent vector !!

rg_dList<rg_REAL> rg_BSplineIntersector::findCharParam(const rg_NUBSplineCurve3D &curve)
{
	rg_INT NumOfNonZeroLengthKnotSpan = curve.getNumOfNonZeroLengthKnotSpan();

	//rg_PolynomialWithBound** polyFormCurve = curve.makePolynomialFormCurve();
	rg_PolynomialWithBound** polyFormCurve = curve.makePolynomialFormCurveUsingCurveDecomp();

	rg_PolynomialWithBound** poly_ht     = new rg_PolynomialWithBound* [ 3 ];
	rg_PolynomialWithBound** poly_ht_dev = new rg_PolynomialWithBound* [ 3 ];

	rg_INT i = 0;
	for(i = 0;i < 3;i++)
	{
		poly_ht[ i ]     = new rg_PolynomialWithBound[NumOfNonZeroLengthKnotSpan];
		poly_ht_dev[ i ] = new rg_PolynomialWithBound[NumOfNonZeroLengthKnotSpan];

		for(rg_INT j = 0;j < NumOfNonZeroLengthKnotSpan;j++)
		{
			poly_ht[ i ][ j ] = polyFormCurve[ i ][ j ].makeDerivative();
			poly_ht_dev[ i ][ j ] = poly_ht[ i ][ j ].makeDerivative();
		}
	}

	//------------------------------------------------------------------------------
	//rg_BzCurve3D* BezierCurveList = curve.decomposeCurveIntoBezierSegment(); //test
	//------------------------------------------------------------------------------


    rg_dList<rg_REAL> inflectParam;

	//rg_INT indexOfDistinctKnotValue = 0;
	rg_REAL* distinctKnotValue = curve.getDistinctKnotValues();

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		rg_PolynomialWithBound polynomialFormCurve = poly_ht[ 0 ][ i ] * poly_ht_dev[ 1 ][ i ] - poly_ht[ 1 ][ i ] * poly_ht_dev[ 0 ][ i ];
		rg_ComplexNumber *root = polynomialFormCurve.solve();

/*
		rg_REAL bounds[ 8 ];
		bounds[ 0 ] = poly_ht[ 0 ][ i ].getLowerBound();
		bounds[ 1 ] = poly_ht[ 0 ][ i ].getUpperBound();
		bounds[ 2 ] = poly_ht[ 1 ][ i ].getLowerBound();
		bounds[ 3 ] = poly_ht[ 1 ][ i ].getUpperBound();
		bounds[ 4 ] = poly_ht_dev[ 0 ][ i ].getLowerBound();
		bounds[ 5 ] = poly_ht_dev[ 0 ][ i ].getUpperBound();
		bounds[ 6 ] = poly_ht_dev[ 1 ][ i ].getLowerBound();
		bounds[ 7 ] = poly_ht_dev[ 1 ][ i ].getUpperBound();

		QuickSort(bounds, 0, 7);

		rg_REAL lowerBound = bounds[ 0 ];
		rg_REAL upperBound = bounds[ 7 ];
*/		

//----------------------------------------------------------------------------------------

		rg_INT numOfTotalRoot = polynomialFormCurve.getDegree();
		rg_REAL *r = new rg_REAL[numOfTotalRoot];
		rg_INT j = 0;
		for(j = 0; j < numOfTotalRoot; j++) r[j] = 0.0;
        //이 부분을 다시 살리고 아래 // 표시가 있는 부분까지 지운다.
/*
		rg_ComplexNumber* root;
		rg_INT numOfTotalRoot;

		BezierCurveList[ i ].inflectionPointByLeeMethod(root, numOfTotalRoot);
		rg_REAL *r = new rg_REAL[numOfTotalRoot];
		for(rg_INT j = 0; j < numOfTotalRoot; j++) r[j] = 0.0;*/
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------

		j = 0;
		for(rg_INT k = 0; k < numOfTotalRoot; k++)
		{
			if(root[k].isPureRealNumber() && rg_LE(distinctKnotValue[ i ], root[k].getRealNumber()) && rg_LT(root[k].getRealNumber(), distinctKnotValue[i + 1]))
			{
				r[j++] = root[k].getRealNumber();
			}
		}
        //이 부분을 다시 살리고 아래 // 표시가 있는 부분까지 지운다.

/*		j = 0;
		for(rg_INT k = 0; k < numOfTotalRoot; k++)
		{
			if(rg_LE(distinctKnotValue[ i ], distinctKnotValue[ i ] + root[k].getRealNumber()*(distinctKnotValue[i + 1] - distinctKnotValue[ i ])) 
		    && rg_LT(distinctKnotValue[ i ] + root[k].getRealNumber()*(distinctKnotValue[i + 1] - distinctKnotValue[ i ]), distinctKnotValue[i + 1]))
			{
				r[j++] = root[k].getRealNumber();
			}
		}*/
//----------------------------------------------------------------------------------------

		delete[] root;
		QuickSort(r, 0, numOfTotalRoot-1);

		//inflectParam.add(distinctKnotValue[indexOfDistinctKnotValue]);
		//indexOfDistinctKnotValue++;
		inflectParam.add(distinctKnotValue[ i ]);
		
		for(rg_INT l = 0; l < numOfTotalRoot; l++)
		{
			if( rg_NZERO(r[l]) )
				inflectParam.add(r[l]);
		}
		delete[] r;

	}
/*
	if(indexOfDistinctKnotValue < NumOfNonZeroLengthKnotSpan + 2)
	{
		for(i = indexOfDistinctKnotValue;i < NumOfNonZeroLengthKnotSpan + 1;i++)
			inflectParam.add(distinctKnotValue[ i ]);
	}
*/
	inflectParam.add(distinctKnotValue[NumOfNonZeroLengthKnotSpan]);

	delete[] distinctKnotValue;

    inflectParam.reset();

    rg_dList<rg_REAL> charParam;
   
    rg_REAL t0 = inflectParam.getEntity();
    charParam.add(t0);
    inflectParam.setCurrentNext();

	rg_NUBSplineCurve3D ht = curve.makeDerivative();

    do
    {
        rg_REAL t1 = inflectParam.getEntity();
        
        makeSimpleParam(ht, t0, t1, charParam);
        charParam.add(t1);
        t0 = t1;
        inflectParam.setCurrentNext();
    } while(!inflectParam.isHead());


	for(i = 0;i < 3;i++)
	{
		delete[] poly_ht[ i ];
		delete[] poly_ht_dev[ i ];
	}

	delete[] poly_ht;
	delete[] poly_ht_dev;

	return charParam;

/*
	rg_NUBSplineCurve3D ht     = curve.makeDerivative();
	rg_NUBSplineCurve3D ht_dev = ht.makeDerivative();
	
	rg_PolynomialWithBound** poly_ht     =     ht.makePolynomialFormCurve();
	rg_PolynomialWithBound** poly_ht_dev = ht_dev.makePolynomialFormCurve();

	rg_INT NumOfNonZeroLengthKnotSpan = curve.getNumOfNonZeroLengthKnotSpan();

    rg_dList<rg_REAL> inflectParam;

	rg_INT indexOfDistinctKnotValue = 0;
	rg_REAL* distinctKnotValue = curve.getDistinctKnotValues();

	for(rg_INT i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		rg_Polynomial polynomialFormCurve = poly_ht[ 0 ][ i ] * poly_ht_dev[ 1 ][ i ] - poly_ht[ 1 ][ i ] * poly_ht_dev[ 0 ][ i ];
		rg_ComplexNumber *root = polynomialFormCurve.solve();
		
		rg_INT numOfTotalRoot = polynomialFormCurve.getDegree();
		rg_REAL *r = new rg_REAL[numOfTotalRoot];
		for(rg_INT j = 0; j < numOfTotalRoot; j++) r[j] = 0.0;
		
		j = 0;
		for(rg_INT k = 0; k < numOfTotalRoot; k++)
		{
			if(root[k].isPureRealNumber() && rg_BTORexclusive(distinctKnotValue[ i ], root[k].getRealNumber(), distinctKnotValue[i + 1]))
			{
				r[j++] = root[k].getRealNumber();
			}
		}
		
		delete[] root;
		QuickSort(r, 0, numOfTotalRoot-1);

		inflectParam.add(distinctKnotValue[indexOfDistinctKnotValue]);
		indexOfDistinctKnotValue++;
		
		for(rg_INT l = 0; l < numOfTotalRoot; l++)
		{
			if( rg_NZERO(r[l]) )
				inflectParam.add(r[l]);
		}
		delete[] r;

	}

	if(indexOfDistinctKnotValue < NumOfNonZeroLengthKnotSpan)
	{
		for(i = indexOfDistinctKnotValue;i <= NumOfNonZeroLengthKnotSpan;i++)
			inflectParam.add(distinctKnotValue[ i ]);
	}

	delete[] distinctKnotValue;

	return inflectParam;
*/

/*
	rg_NUBSplineCurve3D ht     = curve.makeDerivative();
	rg_NUBSplineCurve3D ht_dev = ht.makeDerivative();
	
	rg_PolynomialWithBound** poly_ht     =     ht.makePolynomialFormCurve();
	rg_PolynomialWithBound** poly_ht_dev = ht_dev.makePolynomialFormCurve();

	rg_INT NumOfNonZeroLengthKnotSpan = curve.getNumOfNonZeroLengthKnotSpan();

    rg_dList<rg_REAL> inflectParam;

	rg_INT indexOfDistinctKnotValue = 0;
	rg_REAL* distinctKnotValue = curve.getDistinctKnotValues();

	for(rg_INT i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		rg_Polynomial polynomialFormCurve = poly_ht[ 0 ][ i ] * poly_ht_dev[ 1 ][ i ] - poly_ht[ 1 ][ i ] * poly_ht_dev[ 0 ][ i ];
		rg_ComplexNumber *root = polynomialFormCurve.solve();
		
		rg_INT numOfTotalRoot = polynomialFormCurve.getDegree();
		rg_REAL *r = new rg_REAL[numOfTotalRoot];
		for(rg_INT j = 0; j < numOfTotalRoot; j++) r[j] = 0.0;
		
		j = 0;
		for(rg_INT k = 0; k < numOfTotalRoot; k++)
		{
			if(root[k].isPureRealNumber() && rg_BTORexclusive(curve.getKnotValue( i ), root[k].getRealNumber(), curve.getKnotValue(i + 1)))
			{
				r[j++] = root[k].getRealNumber();
			}
		}
		
		delete[] root;
		QuickSort(r, 0, numOfTotalRoot-1);

		inflectParam.add(distinctKnotValue[indexOfDistinctKnotValue]);
		indexOfDistinctKnotValue++;
		
		for(rg_INT l = 0; l < numOfTotalRoot; l++)
		{
			//if(rg_NE(r[l], curve.getKnotValue( l )))
				inflectParam.add(r[l]);
		}
		delete[] r;

	}

	delete[] distinctKnotValue;

	return inflectParam;
*/

}

rg_dList<rg_REAL> rg_BSplineIntersector::findCharParamUsingDirectPolynomialGeneration(const rg_NUBSplineCurve3D &curve)
{
	rg_INT NumOfNonZeroLengthKnotSpan = curve.getNumOfNonZeroLengthKnotSpan();

	rg_PolynomialWithBound** polyFormCurve = curve.makePolynomialFormCurve();
	//rg_PolynomialWithBound** polyFormCurve = curve.makePolynomialFormCurveUsingCurveDecomp();

	rg_PolynomialWithBound** poly_ht     = new rg_PolynomialWithBound* [ 3 ];
	rg_PolynomialWithBound** poly_ht_dev = new rg_PolynomialWithBound* [ 3 ];

	rg_INT i = 0;
	for(i = 0;i < 3;i++)
	{
		poly_ht[ i ]     = new rg_PolynomialWithBound[NumOfNonZeroLengthKnotSpan];
		poly_ht_dev[ i ] = new rg_PolynomialWithBound[NumOfNonZeroLengthKnotSpan];

		for(rg_INT j = 0;j < NumOfNonZeroLengthKnotSpan;j++)
		{
			poly_ht[ i ][ j ] = polyFormCurve[ i ][ j ].makeDerivative();
			poly_ht_dev[ i ][ j ] = poly_ht[ i ][ j ].makeDerivative();
		}
	}

	//------------------------------------------------------------------------------
	//rg_BzCurve3D* BezierCurveList = curve.decomposeCurveIntoBezierSegment(); //test
	//------------------------------------------------------------------------------


    rg_dList<rg_REAL> inflectParam;

	//rg_INT indexOfDistinctKnotValue = 0;
	rg_REAL* distinctKnotValue = curve.getDistinctKnotValues();

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		rg_PolynomialWithBound polynomialFormCurve = poly_ht[ 0 ][ i ] * poly_ht_dev[ 1 ][ i ] - poly_ht[ 1 ][ i ] * poly_ht_dev[ 0 ][ i ];
		rg_ComplexNumber *root = polynomialFormCurve.solve();

/*
		rg_REAL bounds[ 8 ];
		bounds[ 0 ] = poly_ht[ 0 ][ i ].getLowerBound();
		bounds[ 1 ] = poly_ht[ 0 ][ i ].getUpperBound();
		bounds[ 2 ] = poly_ht[ 1 ][ i ].getLowerBound();
		bounds[ 3 ] = poly_ht[ 1 ][ i ].getUpperBound();
		bounds[ 4 ] = poly_ht_dev[ 0 ][ i ].getLowerBound();
		bounds[ 5 ] = poly_ht_dev[ 0 ][ i ].getUpperBound();
		bounds[ 6 ] = poly_ht_dev[ 1 ][ i ].getLowerBound();
		bounds[ 7 ] = poly_ht_dev[ 1 ][ i ].getUpperBound();

		QuickSort(bounds, 0, 7);

		rg_REAL lowerBound = bounds[ 0 ];
		rg_REAL upperBound = bounds[ 7 ];
*/		

//----------------------------------------------------------------------------------------

		rg_INT numOfTotalRoot = polynomialFormCurve.getDegree();
		rg_REAL *r = new rg_REAL[numOfTotalRoot];
		
		rg_INT j = 0;
		for(j = 0; j < numOfTotalRoot; j++) r[j] = 0.0;
        //이 부분을 다시 살리고 아래 // 표시가 있는 부분까지 지운다.
/*
		rg_ComplexNumber* root;
		rg_INT numOfTotalRoot;

		BezierCurveList[ i ].inflectionPointByLeeMethod(root, numOfTotalRoot);
		rg_REAL *r = new rg_REAL[numOfTotalRoot];
		for(rg_INT j = 0; j < numOfTotalRoot; j++) r[j] = 0.0;*/
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------

		j = 0;
		for(rg_INT k = 0; k < numOfTotalRoot; k++)
		{
			if(root[k].isPureRealNumber() && rg_LE(distinctKnotValue[ i ], root[k].getRealNumber()) && rg_LT(root[k].getRealNumber(), distinctKnotValue[i + 1]))
			{
				r[j++] = root[k].getRealNumber();
			}
		}
        //이 부분을 다시 살리고 아래 // 표시가 있는 부분까지 지운다.

/*		j = 0;
		for(rg_INT k = 0; k < numOfTotalRoot; k++)
		{
			if(rg_LE(distinctKnotValue[ i ], distinctKnotValue[ i ] + root[k].getRealNumber()*(distinctKnotValue[i + 1] - distinctKnotValue[ i ])) 
		    && rg_LT(distinctKnotValue[ i ] + root[k].getRealNumber()*(distinctKnotValue[i + 1] - distinctKnotValue[ i ]), distinctKnotValue[i + 1]))
			{
				r[j++] = root[k].getRealNumber();
			}
		}*/
//----------------------------------------------------------------------------------------

		delete[] root;
		QuickSort(r, 0, numOfTotalRoot-1);

		//inflectParam.add(distinctKnotValue[indexOfDistinctKnotValue]);
		//indexOfDistinctKnotValue++;
		inflectParam.add(distinctKnotValue[ i ]);
		
		for(rg_INT l = 0; l < numOfTotalRoot; l++)
		{
			if( rg_NZERO(r[l]) )
				inflectParam.add(r[l]);
		}
		delete[] r;

	}
/*
	if(indexOfDistinctKnotValue < NumOfNonZeroLengthKnotSpan + 2)
	{
		for(i = indexOfDistinctKnotValue;i < NumOfNonZeroLengthKnotSpan + 1;i++)
			inflectParam.add(distinctKnotValue[ i ]);
	}
*/
	inflectParam.add(distinctKnotValue[NumOfNonZeroLengthKnotSpan]);

	delete[] distinctKnotValue;

    inflectParam.reset();

    rg_dList<rg_REAL> charParam;
   
    rg_REAL t0 = inflectParam.getEntity();
    charParam.add(t0);
    inflectParam.setCurrentNext();

	rg_NUBSplineCurve3D ht = curve.makeDerivative();

    do
    {
        rg_REAL t1 = inflectParam.getEntity();
        
        makeSimpleParam(ht, t0, t1, charParam);
        charParam.add(t1);
        t0 = t1;
        inflectParam.setCurrentNext();
    } while(!inflectParam.isHead());


	for(i = 0;i < 3;i++)
	{
		delete[] poly_ht[ i ];
		delete[] poly_ht_dev[ i ];
	}

	delete[] poly_ht;
	delete[] poly_ht_dev;

	return charParam;
}


rg_dList<rg_RQBzCurve2D> rg_BSplineIntersector::approximateRQBzCurves(rg_NUBSplineCurve3D &curve, rg_dList<rg_REAL> &param)
{
	rg_dList<rg_RQBzCurve2D> rqBzCurveList;

	param.reset();
	rg_REAL t0 = param.getEntity();
	param.setCurrentNext();
	
	do
	{
		rg_REAL t1 = param.getEntity();
		rqBzCurveList.add(makeOneRQBzCurve(curve, t0, t1));
		t0 = t1;
		param.setCurrentNext();
	} while(!param.isHead());

	return rqBzCurveList;
}

rg_RQBzCurve2D rg_BSplineIntersector::makeOneRQBzCurve(rg_NUBSplineCurve3D &curve, const rg_REAL &t0, const rg_REAL &t1)
{
	
	rg_NUBSplineCurve3D  curveDev = curve.makeDerivative();
	rg_Point2D pt[2];
	rg_Point2D pt_prime[2];

	pt[0]       = (curve.evaluatePt(t0)).evaluatePt2D();
	pt_prime[0] = (curveDev.evaluatePt(t0)).evaluatePt2D();

	pt[1]       = (curve.evaluatePt(t1)).evaluatePt2D();
	pt_prime[1] = (curveDev.evaluatePt(t1)).evaluatePt2D();

	rg_RQBzCurve2D rqcurve;
	rqcurve.makeRQBezier(pt[0], pt_prime[0], 
		                 pt[1], pt_prime[1],
						 (curve.evaluatePt((t0+t1)/2.0)).evaluatePt2D());

	return rqcurve;
/*
	rg_PolynomialWithBound** polyFormCurve = curve.makePolynomialFormCurve();
	rg_INT NumOfNonZeroLengthKnotSpan = curve.getNumOfNonZeroLengthKnotSpan();
	rg_PolynomialWithBound** curveDev = new rg_PolynomialWithBound* [ 3 ];

	for(rg_INT i = 0;i < 3;i++)
	{
		curveDev[ i ] = new rg_PolynomialWithBound[NumOfNonZeroLengthKnotSpan];

		for(rg_INT j = 0;j < NumOfNonZeroLengthKnotSpan;j++)
		{
			curveDev[ i ][ j ] = polyFormCurve[ i ][ j ].makeDerivative();
		}
	}	

	rg_REAL* DistinctKnotValues = curve.getDistinctKnotValues();

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		if(rg_EQ(DistinctKnotValues[ i ], t0))
			break;
	}


	rg_Point2D pt[2];
	rg_Point2D pt_prime[2];
	

	pt[0]       = curve.evaluatePt(t0);
	pt_prime[0].setX(pt[0].getX() + curveDev[ 0 ][ i ].evaluatePolynomial(t0));
	pt_prime[0].setY(pt[0].getY() + curveDev[ 0 ][ i ].evaluatePolynomial(t0));

	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		if(rg_EQ(DistinctKnotValues[ i ], t1))
			break;
	}

	delete[] DistinctKnotValues;

	pt[1]       = curve.evaluatePt(t1);
	pt_prime[1].setX(pt[1].getX() + curveDev[ 1 ][ i ].evaluatePolynomial(t1));
	pt_prime[1].setY(pt[1].getX() + curveDev[ 1 ][ i ].evaluatePolynomial(t1));


	for(i = 0;i < 3;i++)
	{
		delete[] curveDev[ i ];
		delete[] polyFormCurve[ i ];
	}

	delete[] curveDev;
	delete[] polyFormCurve;

	rg_RQBzCurve2D rqcurve;
	rqcurve.makeRQBezier(pt[0], pt_prime[0], 
		                 pt[1], pt_prime[1],
						 rg_Point2D(curve.evaluatePt((t0+t1)/2.0)));

	return rqcurve;*/
}

rg_dList<rg_REAL*> rg_BSplineIntersector::makeSeed(rg_dList<rg_RQBzCurve2D> &rqcurve_s, rg_dList<rg_REAL> subParam_s,
	                                        rg_dList<rg_RQBzCurve2D> &rqcurve_t, rg_dList<rg_REAL> subParam_t)
{
	rg_dList <rg_REAL*> seedParam4TwoCurve;
	rg_dList <rg_Point2D> intersectPointList;

	rqcurve_s.reset();
	subParam_s.reset();
	rg_REAL s0 = subParam_s.getEntity();
	subParam_s.setCurrentNext();

	do
	{
		rg_REAL s1 = subParam_s.getEntity();

		rqcurve_t.reset();
		subParam_t.reset();
		rg_REAL t0 = subParam_t.getEntity();
		subParam_t.setCurrentNext();
		
		do
		{
			rg_REAL t1 = subParam_t.getEntity();
	
			// Using the bounding box, 
            // it can be omitted to  compute the intersetions between two curves.
            rg_RQBzCurve2D   rqbz_s = rqcurve_s.getEntity();
            rg_RQBzCurve2D   rqbz_t = rqcurve_t.getEntity();
            rg_BoundingBox2D box_s  = rqbz_s.makeBoundingBox();
            rg_BoundingBox2D box_t  = rqbz_t.makeBoundingBox();
			if( box_s.isOverlapped(box_t) )
			{
				rg_IntersectFunc::intersectRQBzCurveVsRQBzCurve(rqbz_s, s0, s1, 
											  rqbz_t, t0, t1, 
											  intersectPointList, 
											  seedParam4TwoCurve);
			}
			t0 = t1;
			rqcurve_t.setCurrentNext();
			subParam_t.setCurrentNext();
		} while(!rqcurve_t.isHead());

		s0 = s1;
		rqcurve_s.setCurrentNext();
		subParam_s.setCurrentNext();
	} while(!rqcurve_s.isHead());

	return seedParam4TwoCurve;
}

rg_Point2D rg_BSplineIntersector::iterationWithSeed( rg_NUBSplineCurve3D &curve_s, const rg_REAL &param_s,
										       rg_NUBSplineCurve3D &curve_t, const rg_REAL &param_t)
{
/*
	rg_REAL newParam_s = param_s;
	rg_REAL newParam_t = param_t;

	rg_Point2D initialPtOnCurve_s = curve_s.evaluatePt(newParam_s);
	rg_Point2D initialPtOnCurve_t = curve_t.evaluatePt(newParam_t);

	rg_Point2D ptOnCurve_s = initialPtOnCurve_s;
	rg_Point2D ptOnCurve_t = initialPtOnCurve_t;

	rg_Point2D distance = ptOnCurve_s - ptOnCurve_t;

	rg_BzCurve2D devCurve_s = curve_s.makeDerivative();
	rg_BzCurve2D devCurve_t = curve_t.makeDerivative();

	rg_INT count = 0;
	while( rg_GT(distance.magnitudeSquare(), 1.0e-12) && (count < 8) )
	{
		rg_Point2D tangentVector_s = devCurve_s.evaluatePt(newParam_s);
		rg_Point2D tangentVector_t = devCurve_t.evaluatePt(newParam_t);

		rg_REAL temp = tangentVector_s*tangentVector_t;
		newParam_s += (ptOnCurve_t-ptOnCurve_s)*tangentVector_t/temp;
		newParam_t += (ptOnCurve_s-ptOnCurve_t)*tangentVector_s/(-temp);
		
		ptOnCurve_s.setPoint( curve_s.evaluatePt(newParam_s) );
		ptOnCurve_t.setPoint( curve_t.evaluatePt(newParam_t) );

		distance = ptOnCurve_s - ptOnCurve_t;
		count++;
	}

	return (ptOnCurve_s + ptOnCurve_t)/2;
*/

	rg_REAL newParam_s = param_s;
	rg_REAL newParam_t = param_t;

	rg_Point2D initialPtOnCurve_s = (curve_s.evaluatePt(newParam_s)).evaluatePt2D();
	rg_Point2D initialPtOnCurve_t = (curve_t.evaluatePt(newParam_t)).evaluatePt2D();

	rg_Point2D ptOnCurve_s = initialPtOnCurve_s;
	rg_Point2D ptOnCurve_t = initialPtOnCurve_t;

	rg_Point2D distance = ptOnCurve_s - ptOnCurve_t;

	rg_NUBSplineCurve3D devCurve_s = curve_s.makeDerivative();
	rg_NUBSplineCurve3D devCurve_t = curve_t.makeDerivative();

	rg_INT count = 0;
	while( rg_GT(distance.magnitudeSquare(), 1.0e-12) && (count < 8) )
	{
		rg_Point2D tangentVector_s = (devCurve_s.evaluatePt(newParam_s)).evaluatePt2D();
		rg_Point2D tangentVector_t = (devCurve_t.evaluatePt(newParam_t)).evaluatePt2D();

		rg_REAL temp = tangentVector_s*tangentVector_t;
		newParam_s += (ptOnCurve_t-ptOnCurve_s)*tangentVector_t/temp;
		newParam_t += (ptOnCurve_s-ptOnCurve_t)*tangentVector_s/(-temp);
		
		ptOnCurve_s.setPoint( (curve_s.evaluatePt(newParam_s)).evaluatePt2D() );
		ptOnCurve_t.setPoint( (curve_t.evaluatePt(newParam_t)).evaluatePt2D() );

		distance = ptOnCurve_s - ptOnCurve_t;
		count++;
	}

	return (ptOnCurve_s + ptOnCurve_t)/2;

}

void rg_BSplineIntersector::makeSimpleParam(rg_NUBSplineCurve3D& curve, const rg_REAL& t0, const rg_REAL& t1,
                                         rg_dList<rg_REAL>& charParam)
{
    rg_Point2D v1 = (curve.evaluatePt(t0)).evaluatePt2D();
    rg_Point2D v2 = (curve.evaluatePt(t1)).evaluatePt2D();
    rg_Point2D p  = (curve.evaluatePt((t0+t1)/2.0)).evaluatePt2D();

    rg_REAL  D  = v1.getX()*v2.getY() - v2.getX()*v1.getY();

    rg_REAL  a  = (v2.getY()*p.getX() - v2.getX()*p.getY())/D;
    rg_REAL  b  = (v1.getX()*p.getY() - v1.getY()*p.getX())/D;

    rg_REAL  theta = acos(v1%v2/v1.magnitude()/v2.magnitude());

    if(rg_POS(a) && rg_POS(b))
    {
        return;
    }

    else
    {
        theta = 2*pi - theta;
        
        rg_REAL *t = getConjugateTangentParam(curve, t0, t1);

        charParam.add(t[0]);
        charParam.add(t[1]);
            
        delete[] t;
    }
}

rg_REAL* rg_BSplineIntersector::getConjugateTangentParam(rg_NUBSplineCurve3D& curve, const rg_REAL& t0, const rg_REAL& t1)
{
    rg_Point2D v0 = (curve.evaluatePt(t0)).evaluatePt2D();
    rg_Point2D v1 = (curve.evaluatePt(t1)).evaluatePt2D();
 
    rg_REAL tangent0 = v0.getY()/v0.getX();
    rg_REAL tangent1 = v1.getY()/v1.getX();


	rg_PolynomialWithBound** curvePolynomial = curve.makePolynomialFormCurveUsingCurveDecomp();
	rg_INT NumOfNonZeroLengthKnotSpan = curve.getNumOfNonZeroLengthKnotSpan();

	rg_INT i = 0;
	for(i = 0;i < NumOfNonZeroLengthKnotSpan;i++)
	{
		if(rg_LE(curvePolynomial[ 0 ][ i ].getLowerBound(), t0) && rg_LT(t0, curvePolynomial[ 0 ][ i ].getUpperBound()))
			break;
	}

    rg_Polynomial p0 = curvePolynomial[1][i] - tangent0*curvePolynomial[0][i];
    rg_Polynomial p1 = curvePolynomial[1][i] - tangent1*curvePolynomial[0][i];

    delete[] curvePolynomial;

    rg_ComplexNumber *root0 = p0.solve();
    rg_ComplexNumber *root1 = p1.solve();

    rg_REAL *result = new rg_REAL[2];
    result[0] = t0;
    result[1] = t1;
    
    const rg_REAL errBound = 0.1;
        
    for(i=0; i<p0.getDegree(); i++)
    {
        if( root0[i].isPureRealNumber() )
        {
            rg_REAL r = root0[i].getRealNumber();
            if(rg_BTORexclusive(result[0], r, t1) ) 
            {
                result[1] = r;
            }
        }

        if( root1[i].isPureRealNumber() )
        {
            rg_REAL r = root1[i].getRealNumber();
            if(rg_BTORexclusive(t0, r, result[1]) )
            {
                result[0] = r;
            }
        }
    }

    delete[] root0;
    delete[] root1;
    return result;
}


