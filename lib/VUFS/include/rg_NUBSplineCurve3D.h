//******************************************************************************
//
//	  FILENAME    : rg_NUBSplineCurve3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_NUBSplineCurve3D 
//           which define Non-uniform B-Spline rg_Curve and its property. 
//                          
//	  CLASS NAME  : rg_NUBSplineCurve3D
//
//    BASE CLASS  : rg_BSplineCurve3D
//      
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//
//    START DATE  : 21 Jun 1996    
//
//    HISTORY     :
//          BY Young-Song Cho.  13 Jul. 1997
//            make : rg_INT        getNumOfKnotSpan(); 
//            make : rg_BzCurve3D* separateBSplineIntoBezier();   
//            make : rg_sListByPtr*     intersectOfCubicBSplineAndPlane(const rg_Point3D &plane);
//
//          BY Young-Song Cho.  24 Jul. 1997
//            make : rg_BzCurve3D* pullOutCubicBezierForKnotSpan(
//                                  const rg_REAL &prevKnot, 
//                                  const rg_REAL &nextKnot );
//            make : rg_BzCurve3D* pullOutCubicBezierForKnotSpan(
//                                  const rg_INDEX &indexOfKnotSpan);
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
//                                                    const rg_Matrix* powerCoeff )
//
//                  void bzCurvesToNUBSplineCurve( const rg_INT& numOfSegment, 
//                                                 const rg_BzCurve3D* bzCurves )
//
//           By Dong-Gyou Lee 26 Mar. 1998
//                  void reparameterizationKnotVector()
//
//			 By JoongHyun Ryu Apr.11 1998
//					void multipleKnotInsertion(const rg_REAL& insertedKnot, const rg_INT& NumberOfInsertingMultipleKnots);
//
//			 By JoongHyun Ryu Apr.11 1998
//					rg_INT findTheCorrespondingKnotSpan(const rg_REAL& insertedKnot);
//
//			 By JoongHyun Ryu Apr.16 1998
//					rg_INT getNumOfNonZeroLengthKnotSpan() const;
//
//			 By JoongHyun Ryu Apr.16 1998
//					rg_BzCurve3D* decomposeCurveIntoBezierSegment();
//					
//           By Young-Song Cho Jul.7 1998 
//                  rg_FLAG makeCompositeCurveWithC0(const rg_NUBSplineCurve3D &curve1,
//                                                const rg_NUBSplineCurve3D &curve2);
//           By Young-Song Cho Jul.7 1998 
//                  void formLine( const rg_Point3D& start, const rg_Point3D& end, const rg_INT &order = 4 );
//
//			 By JoongHyun Ryu May 28 1999
//						rg_REAL** getEntirePolynomialCurve() const;
//
//           By Joonghyun Ryu Jun 23 1999
//                      void knotRefinement( const rg_REAL* &insertingKnotValues, const rg_INT& numOfinsertingKnotValues); 
//
//	         By Taeboom Jang 1999. 8.27
//                     rg_dList<rg_Point3D> intersectWithPlaneForCubic(const rg_Plane3D& plane)
//
//	         By Taeboom Jang 1999. 8.31
//                     rg_REAL getStartParameter() const;
//                     rg_REAL getEndParameter() const;
//
//           By Taeboom Jang    1999. 10.14 
//                     rg_Point3D  getStartPoint() const;
//                     rg_Point3D  getEndPoint() const;
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea   	  	
//
//******************************************************************************

#ifndef _RG_NUBSPLINECURVE3D_H
#define _RG_NUBSPLINECURVE3D_H

#include "rg_BSplineCurve3D.h"

#include "rg_Const.h"
//#include "DefConst.h"
#include "rg_Point3D.h"
#include "rg_BzCurve3D.h"
#include "rg_ListByPtr.h"
#include "rg_Polynomial.h"
#include "rg_PolynomialWithBound.h"
#include "rg_Polyline3D.h"
#include "rg_Polyline2D.h"

class rg_NUBSplineCurve3D : public rg_BSplineCurve3D
{
protected:
    rg_REAL* knotVector;

protected:
static 	rg_Point3D*   makeFilteredPassingPts(       rg_INT&     numOfPts,
	    	                                  const rg_Point3D* const pts );
		                                       

public:
////    Constructor & Destructor
    rg_NUBSplineCurve3D();
    rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                      const rg_Planarity    &newPlanarity );
    rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                      const rg_Planarity    &newPlanarity,
                      const rg_ORDER        &newOrder );
    rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                      const rg_Planarity    &newPlanarity,
                      const rg_INT          &num );
    rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                      const rg_Planarity    &newPlanarity,
                      const rg_INT          &num, 
                      rg_Point3D*              newControlP );
    rg_NUBSplineCurve3D( const rg_ORDER &newOrder );
    rg_NUBSplineCurve3D( const rg_INT &num );
    rg_NUBSplineCurve3D( const rg_INT &num, 
                      rg_Point3D*     newControlP );
    rg_NUBSplineCurve3D( const rg_ORDER &newOrder, 
                      const rg_INT   &num, 
                      rg_Point3D*       newControlP );
    rg_NUBSplineCurve3D( const rg_ORDER &newOrder, 
                      const rg_INT   &num, 
                      rg_Point3D*       newControlP,
                      rg_REAL*        newKnotVector );
    ////  Constructor      : March 13 1997
    rg_NUBSplineCurve3D( const unsigned rg_INT &newID, 
                      const rg_Planarity    &newPlanarity,
                      const rg_ORDER        &newOrder, 
                      const rg_INT          &num, 
                      rg_Point3D*              newControlP,
                      rg_REAL*               newKnotVector );
    rg_NUBSplineCurve3D( const rg_BSplineCurve3D &curve);
    ////  Copy Constructor : March 13 1997
    rg_NUBSplineCurve3D( const rg_NUBSplineCurve3D &curve);

    virtual ~rg_NUBSplineCurve3D();

////    Get Functions.
    rg_REAL  getKnotValue( const rg_INT &kIndex ) const;
    rg_REAL* getKnotVector() const;
    rg_INT   getNumOfKnotSpan() const;
    rg_INT   getNumOfNonZeroLengthKnotSpan() const;
    rg_REAL  getStartParameter() const;
    rg_REAL  getEndParameter() const;

    rg_Point3D  getStartPoint() const;
    rg_Point3D  getEndPoint() const;

	rg_REAL*  getDistinctKnotValues() const;
	rg_INT*  getIndexOfNonZeroBasisInKnotSpan(const rg_INDEX& indexOfKnotValue) const;
    rg_Polynomial getBasisWithPolynomialFormIn(const rg_INDEX& indexOfCtrlPoint, const rg_DEGREE& degree, const rg_INDEX& indexOfKnotValue) const;
    rg_PolynomialWithBound* getBasisWithPolynomialFormIn(const rg_INDEX& indexOfCtrlPoint) const;
	rg_PolynomialWithBound** makePolynomialFormCurve() const;
	rg_PolynomialWithBound** makePolynomialFormCurveUsingCurveDecomp() const;
	rg_PolynomialWithBound** makePolynomialFormCurveUsingCurveDecomp(rg_REAL& timeForKnotRefinement1,
	                                                              rg_REAL& timeForKnotRefinement2) const;

	// B-spline basis to Power basis
	rg_REAL*  getDistributionPolynomialsInOneGraph(rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths) const;
	rg_REAL*  getDistributionPolynomialsInOneGraphRevised(rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths) const;
	rg_REAL*  getDistributionPolynomialsInOneGraph(rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths,
												rg_REAL& time1, rg_REAL& time2) const;
	rg_REAL*  getDistributionPolynomialsInOneGraph1(rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths,
												rg_REAL& time1) const;
	rg_REAL*  getDistributionPolynomialsInOneGraph2(rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths,
												rg_REAL& time2) const;
	void getEntirePolynomialCurve(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate) const;
	void getEntirePolynomialCurve(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
												rg_REAL   & timeForOurAlgorithm1,
												rg_REAL   & timeForOurAlgorithm2,
												rg_REAL   & timeForOurAlgorithm3) const;
	void getEntirePolynomialCurve1(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
												rg_REAL   & timeForOurAlgorithm1) const;
	void getEntirePolynomialCurve2(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
												rg_REAL   & timeForOurAlgorithm2) const;
	void getEntirePolynomialCurve3(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
												rg_REAL   & timeForOurAlgorithm3) const;
    rg_Polynomial*  getEntirePolynomialOfX() const;
    rg_Polynomial** getEntirePolynomialCurve() const;    

    rg_REAL getCoeffOfMatrixRepresentation(rg_INT multiplicity, rg_INT row, rg_INT col, rg_INT indexOfInsertedKnot, rg_INT indexOfCtrlPt);
	rg_Matrix* getMatrixRepresentation();

	// for only non-periodic B-spline !!
	// (for experimental analysis alone)
	rg_REAL getCoeffOfMatrixRepresentationInUniformKnot(rg_INT row, rg_INT col);
	rg_Matrix* getMatrixRepresentationInUniformKnot();

	rg_Matrix* getMatrixRepresentationByDE();

	// Using Taylor's Expansion
	void getPiecewisePolynomialCurveInaPowerForm(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate);

	// For experimental analysis--------------------------------

	// Saving truncated basis polynomials
	void makePiecewisePowerFormPolyUsingDE(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate, rg_REAL*** & truncatedBasisPolys);

	void makeUpdatedCurveUsingDE(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate, rg_REAL*** & truncatedBasisPolys,
		                         rg_INT* indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts);

	void makeTruncatedBasisPolynomial(rg_REAL* & truncatedBasisPolynomial, rg_INT contributionKnotSpan, /*rg_INT IndexOfDestination,*/ rg_INT** graph, rg_INT noOfAllPossiblePaths);

	// 
	void makeUpdatedCurveUsingTE(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
		                         rg_INT* indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts);

	void makePiecewisePowerFormPolyUsingKR(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate);
    void makeUpdatedCurveUsingKR(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate,
	                             rg_INT* indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts);

	//----------------------------------------------------------

	void makePiecewisePowerFormPolynomialOfBSplineBasisUsingDE(rg_REAL*** & piecewisePolynomialsInPowerFormOfBasis);
	rg_REAL* getPowerFormPolynomialUsingBD(const rg_INDEX& index, const rg_INDEX& knotSpanIndex);
	rg_REAL* getPowerFormPolynomialUsingBDRecursion(const rg_INDEX& index, const rg_INDEX& knotSpanIndex);
	void makePiecewisePowerFormPolynomialOfBSplineBasisUsingBD(rg_REAL*** & piecewisePolynomialsInPowerFormOfBasis);

	void makePiecewisePowerFormPolynomialOfBSplineBasisUsingBDRecursion(rg_REAL*** & piecewisePolynomialsInPowerFormOfBasis);

	// -----------------------------
    rg_NUBSplineCurve3D evaluateCurveSegment( const rg_REAL& start,
                                              const rg_REAL& end ) const;

////    Set Functions.
    void    setInitialKnotVector();
    void    setKnotValue( const rg_INT  &kIndex, 
                          const rg_REAL &newKnotValue );
    void    setKnotVector( const rg_INT &numOfKnot,
                           const rg_REAL* const newKnotVector );

////    BasisFunction & rg_Point3D Evaluating.
    virtual rg_REAL evaluateBasisFunc( const rg_INT  &index, 
                                    const rg_REAL &param,
                                    const rg_INT  &Order ) const;

    rg_REAL* evaluateMultiBasisFunc( const rg_INDEX     &knotIndex,
                                  const rg_PARAMETER &u,
                                  const rg_ORDER     &order) const;

	rg_REAL getCoefficientForDerivativeEvaluation(const rg_INDEX& index1, const rg_INDEX& index2, const rg_INDEX& index3) const;

	rg_REAL evaluateBasisFuncDerivative(const rg_INDEX& index, const rg_PARAMETER& u, const rg_ORDER& derivaOrder) const;

	rg_REAL evaluateBasisFuncDerivativeUsingRecursion(const rg_INDEX& index, const rg_PARAMETER& u, const rg_ORDER& derivaOrder, const rg_ORDER& degree) const;

    virtual rg_Point3D   evaluatePt( const rg_REAL &param );
    virtual rg_Point3D   evaluatePt( const rg_REAL &param ) const;
    virtual rg_Point3D*  evaluatePtsInEvenParameter( const rg_INT &noOfEvaluatedPoint ) const;
    rg_Point3D*  evaluatePt_Plus( const rg_INT &noOfEvaluatedPoint );
//    virtual rg_Polyline2D makePolyline2DInEvenParameter(const rg_INT &numOfPts) const;
//    virtual rg_Polyline3D makePolyline3DInEvenParameter(const rg_INT &numOfPts) const;
    rg_Polyline2D makePolyline2DConsideringKnots(const rg_INT &numOfPts) const;
////    Derivative rg_Curve & Curvature
    rg_NUBSplineCurve3D makeDerivative() const;

	//  signed curvature.
    rg_REAL  getCurvatureInXY( const rg_REAL &param );
    rg_REAL  getCurvatureInXZ( const rg_REAL &param );
    rg_REAL  getCurvatureInYZ( const rg_REAL &param );

////    Fundamental Geometric Algorithm.
	rg_INT findTheCorrespondingKnotSpan(const rg_REAL& insertingKnot);

    void knotInsertion( const rg_REAL &insertingKnotValue );

	//void multipleKnotInsertion(const rg_REAL& insertingKnot, const rg_INT& NumberOfInsertingMultipleKnots);

	void knotRefinement( rg_REAL* &insertingKnotValues, const rg_INT& numOfinsertingKnotValues);

    //  This function separates cubic non-uniform B-spline curve 
    //      into cubic Bezier curves.
    rg_BzCurve3D* separateBSplineIntoBezier();

    //  This function decompose non-uniform B-spline curve into Bezier curves.
	rg_BzCurve3D* decomposeCurveIntoBezierSegment() const;

	rg_BzCurve3D* decomposeCurveIntoBezierSegmentUsingKnotRefinement() const;

    //  This function pulls out a cubic Bezier curve with responding to 
    //  a single knot span.
    rg_BzCurve3D* pullOutCubicBezierForKnotSpan(const rg_REAL &prevKnot, 
                                              const rg_REAL &nextKnot);
    rg_BzCurve3D* pullOutCubicBezierForKnotSpan(const rg_INDEX &indexOfKnotSpan);

    rg_FLAG makeCompositeCurveWithC0(const rg_NUBSplineCurve3D &curve1,
                                  const rg_NUBSplineCurve3D &curve2);
    void formLine( const rg_Point3D& start, const rg_Point3D& end, const rg_INT &order = 4 );

    rg_sListByPtr* intersectOfCubicBSplineAndPlane(const rg_Plane3D &plane) const;
    rg_dList<rg_Point3D> intersectWithPlaneForCubic(const rg_Plane3D& plane) const;
	rg_REAL* inflectionPointByHodograph(rg_INT& numIPts) const;

    void  appendWithC0(      rg_NUBSplineCurve3D   next,
                       const rg_REAL connectedKnot =0.5 );

    void  reverseTrace();
////	Conversion between power spline & b-spline form.
   	void powerSplineToNUBSplineCurve( const rg_DEGREE& dgr,
			                          const rg_INT& numOfSegment,
                                      rg_REAL** paramValues,
							          const rg_Matrix* powerCoeff );

	void bzCurvesToNUBSplineCurve( const rg_INT& numOfSegment,
		                           const rg_BzCurve3D* bzCurves );

////	Reparameterization
	void reparameterizationKnotVector();


////    rg_Curve Generating function
    rg_REAL** makeMatrix( rg_Point3D*         b, 
                         const rg_REAL*  u, 
                         const rg_INT     &L, 
                         const rg_Point3D   &m0, 
                         const rg_Point3D   &mL );

    void     interpolateWithEndCondition( const rg_INT &noOfData, 
                                          rg_Point3D*     fittingData,
                                          rg_Point3D      startTngnt,
                                          rg_Point3D      endTngnt,
                                          rg_REAL*        parametrization = rg_NULL);
	 
	void     interpolate(const rg_INT& numPts,
						 rg_Point3D*   pts,
						 const rg_INT& order = 4,
						 rg_REAL*      parametrization = rg_NULL);
	void     interpolateWithFiltering(const rg_INT& numPts,
						              rg_Point3D*   pts,
						              const rg_INT& order = 4,
						              rg_REAL*      parametrization = rg_NULL);
////    Functions to obtain the information of KNOT
    rg_INT  getNumberOfKnotValues() const;
    rg_INT  getIndexOfKnot(const rg_REAL& tKnotValue) const;
    rg_INT  getKnotMultiplicity(const rg_REAL& tKnotValue) const;
    rg_INT  getIndexOfKnotSpan( const rg_REAL& tParameter) const;
    rg_REAL findParameterOfNearestPt(const rg_Point3D& pt,
							 	 	 const rg_REAL&    distanceTol=1.0e-3,
									 const rg_REAL&    cosAngleTol=1.0e-3,
									 const rg_REAL&    start=0.0,
									 const rg_REAL&    end=1.0) const;

    rg_REAL findParameterOfNearestPt(const rg_Point3D& pt,
                                     const rg_REAL&    distanceTol,
                                     const rg_REAL&    cosAngleTol,
                                     const rg_REAL&    seed,
			    				     const rg_REAL&    start,
				    			     const rg_REAL&    end) const;

////    Operator Overloading.
    rg_NUBSplineCurve3D& operator =(const rg_NUBSplineCurve3D &curve);

    void removeAll();
	

};

/*
////  Parameterization.
////  functions below moved to 'rg_CurveSurfaceFunc.h & cpp' by Lee, Dong-Gyou 17 Mar 1998 

rg_REAL* rg_GeoFunc::rg_GeoFunc::chordLength(const rg_INT &n, const rg_Point3D* const ptsPassedThrough);
rg_REAL* rg_GeoFunc::centripetal(const rg_INT &n, const rg_Point3D* const ptsPassedThrough);
*/

#endif


