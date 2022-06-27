#include <math.h>
//#include <fstream>
#include <time.h>

#include "rg_NUBSplineIntersector.h"
#include "rg_BzCurve3D.h"
#include "rg_Polynomial.h"
#include "rg_IntersectFunc.h"
#include "sortFunc.h"
#include "rg_RelativeOp.h"
#include "rg_IntervalIntersector.h"
#include "rg_IntsecImplicit3.h"
#include "rg_IntsecImplicit4.h"
#include "rg_IntsecImplicit5.h"
#include "rg_IntsecImplicit6.h"
#include "rg_IntsecImplicitNotUsingMathematica.h"


rg_NUBSplineIntersector::rg_NUBSplineIntersector()
{
}

rg_NUBSplineIntersector::~rg_NUBSplineIntersector()
{
}


rg_dList<rg_Point2D> rg_NUBSplineIntersector::intersectBSplineCurveVsBSplineCurveUsingIntsecInterval(rg_NUBSplineCurve3D& curve_s,
																							 rg_NUBSplineCurve3D& curve_t,
																							 rg_REAL & time)
{
	// splitting into Bezier rg_Curve
	// apply Our Algorithm to each pair of Bezier rg_Curve !!(filtering process not yet implemented)

	clock_t StartTime, EndTime;
    StartTime = clock();

	rg_BzCurve3D* BezierCurveList_s = curve_s.decomposeCurveIntoBezierSegment();
	rg_BzCurve3D* BezierCurveList_t = curve_t.decomposeCurveIntoBezierSegment();

	rg_INT numberOfBzCurvesInCurve_s = curve_s.getNumOfNonZeroLengthKnotSpan();
	rg_INT numberOfBzCurvesInCurve_t = curve_t.getNumOfNonZeroLengthKnotSpan();

	rg_IntervalIntersector BzCurveIntersector;

    rg_dList<rg_Point2D> intersectPointList;
		
	for(rg_INT i = 0;i < numberOfBzCurvesInCurve_s;i++)
	{
		rg_dList<rg_Point2D> tIntersectPointList;

		for(rg_INT j = 0;j < numberOfBzCurvesInCurve_t;j++)
		{
			BzCurveIntersector.intersectBzCurveVsBzCurve( BezierCurveList_s[ i ].evaluateBzCurve2D(),
			                                              BezierCurveList_t[ j ].evaluateBzCurve2D(),
														  time    );

			tIntersectPointList = BzCurveIntersector.getIntersectionPointList();

			intersectPointList.append(tIntersectPointList);
			tIntersectPointList.removeAll();
		}

	}

	EndTime = clock();
	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;

	return intersectPointList;
}


rg_dList<rg_Point2D> rg_NUBSplineIntersector::intersectBSplineCurveVsBSplineCurveUsingHardCodedImplicitization(rg_NUBSplineCurve3D& curve_s,
																									   rg_NUBSplineCurve3D& curve_t,
																									   rg_REAL & time)
{
	// splitting into Bezier rg_Curve
	// apply Our Algorithm to each pair of Bezier rg_Curve !!(filtering process not yet implemented)


	rg_BzCurve3D* BezierCurveList_s = curve_s.decomposeCurveIntoBezierSegment();
	rg_BzCurve3D* BezierCurveList_t = curve_t.decomposeCurveIntoBezierSegment();

	//rg_BzCurve3D* BezierCurveList_s = curve_s.decomposeCurveIntoBezierSegmentUsingKnotRefinement();
	//rg_BzCurve3D* BezierCurveList_t = curve_t.decomposeCurveIntoBezierSegmentUsingKnotRefinement();

	rg_INT numberOfBzCurvesInCurve_s = curve_s.getNumOfNonZeroLengthKnotSpan();
	rg_INT numberOfBzCurvesInCurve_t = curve_t.getNumOfNonZeroLengthKnotSpan();

	rg_dList<rg_Point2D> intersectPointList;

	if(curve_t.getOrder() - 1 == 3)
	{
		clock_t StartTime, EndTime;
		StartTime = clock();

		rg_IntsecImplicit3 BzCurveIntersector;	
		
		for(rg_INT i = 0;i < numberOfBzCurvesInCurve_s;i++)
		{
			rg_dList<rg_Point2D> tIntersectPointList;
			
			for(rg_INT j = 0;j < numberOfBzCurvesInCurve_t;j++)
			{
				tIntersectPointList = BzCurveIntersector.intersectBzCurveVsBzCurve( BezierCurveList_s[ i ].evaluateBzCurve2D(),
																					BezierCurveList_t[ j ].evaluateBzCurve2D(),	
																					time    );
				intersectPointList.append(tIntersectPointList);
				tIntersectPointList.removeAll();
			}
		}
		
		EndTime = clock();

		time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	}

	if(curve_t.getOrder() - 1 == 4)
	{
		clock_t StartTime, EndTime;
		StartTime = clock();

		rg_IntsecImplicit4 BzCurveIntersector;	
		
		for(rg_INT i = 0;i < numberOfBzCurvesInCurve_s;i++)
		{
			rg_dList<rg_Point2D> tIntersectPointList;
			
			for(rg_INT j = 0;j < numberOfBzCurvesInCurve_t;j++)
			{
				tIntersectPointList = BzCurveIntersector.intersectBzCurveVsBzCurve( BezierCurveList_s[ i ].evaluateBzCurve2D(),
																					BezierCurveList_t[ j ].evaluateBzCurve2D(),	
																					time    );
				intersectPointList.append(tIntersectPointList);
				tIntersectPointList.removeAll();
			}
		}
		
		EndTime = clock();

		time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	}

	if(curve_t.getOrder() - 1 == 5)
	{
		clock_t StartTime, EndTime;
		StartTime = clock();

		rg_IntsecImplicit5 BzCurveIntersector;	
		
		for(rg_INT i = 0;i < numberOfBzCurvesInCurve_s;i++)
		{
			rg_dList<rg_Point2D> tIntersectPointList;
			
			for(rg_INT j = 0;j < numberOfBzCurvesInCurve_t;j++)
			{
				tIntersectPointList = BzCurveIntersector.intersectBzCurveVsBzCurve( BezierCurveList_s[ i ].evaluateBzCurve2D(),
																					BezierCurveList_t[ j ].evaluateBzCurve2D(),	
																					time    );
				intersectPointList.append(tIntersectPointList);
				tIntersectPointList.removeAll();
			}
		}
		
		EndTime = clock();

		time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	}

	
	if(curve_t.getOrder() - 1 == 6)
	{
		clock_t StartTime, EndTime;
		StartTime = clock();

		rg_IntsecImplicit6 BzCurveIntersector;	
		
		for(rg_INT i = 0;i < numberOfBzCurvesInCurve_s;i++)
		{
			rg_dList<rg_Point2D> tIntersectPointList;
			
			for(rg_INT j = 0;j < numberOfBzCurvesInCurve_t;j++)
			{
				tIntersectPointList = BzCurveIntersector.intersectBzCurveVsBzCurve( BezierCurveList_s[ i ].evaluateBzCurve2D(),
																					BezierCurveList_t[ j ].evaluateBzCurve2D(),	
																					time    );
				intersectPointList.append(tIntersectPointList);
				tIntersectPointList.removeAll();
			}
		}
		
		EndTime = clock();

		time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;
	}


	return intersectPointList;
}

rg_dList<rg_Point2D> rg_NUBSplineIntersector::intersectBSplineCurveVsBSplineCurveUsingImplicitization(rg_NUBSplineCurve3D& curve_s,
																							  rg_NUBSplineCurve3D& curve_t,
																							  rg_REAL & time)
{
	// splitting into Bezier rg_Curve
	// apply Our Algorithm to each pair of Bezier rg_Curve !!(filtering process not yet implemented)

	clock_t StartTime, EndTime;
    StartTime = clock();

	rg_BzCurve3D* BezierCurveList_s = curve_s.decomposeCurveIntoBezierSegment();
	rg_BzCurve3D* BezierCurveList_t = curve_t.decomposeCurveIntoBezierSegment();

	rg_INT numberOfBzCurvesInCurve_s = curve_s.getNumOfNonZeroLengthKnotSpan();
	rg_INT numberOfBzCurvesInCurve_t = curve_t.getNumOfNonZeroLengthKnotSpan();

	rg_IntsecImplicitNotUsingMathematica BzCurveIntersector;

    rg_dList<rg_Point2D> intersectPointList;
		
	for(rg_INT i = 0;i < numberOfBzCurvesInCurve_s;i++)
	{
		rg_dList<rg_Point2D> tIntersectPointList;

		for(rg_INT j = 0;j < numberOfBzCurvesInCurve_t;j++)
		{
			tIntersectPointList = BzCurveIntersector.intersectBzCurveVsBzCurve1( BezierCurveList_s[ i ].evaluateBzCurve2D(),
																			     BezierCurveList_t[ j ].evaluateBzCurve2D(),
																				 time    );
			intersectPointList.append(tIntersectPointList);
			tIntersectPointList.removeAll();
		}

	}

	EndTime = clock();
	time = (rg_REAL)(EndTime - StartTime)/CLOCKS_PER_SEC;

	return intersectPointList;
}


