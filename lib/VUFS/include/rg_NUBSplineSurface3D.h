//********************************************************************
//
//	  FILENAME    : rg_NUBSplineSurface3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_NUBSplineSurface3D 
//           which define Non-Uniform B-Spline rg_Surface and its property. 
//                          
//	  CLASS NAME  : rg_NUBSplineSurface3D
//
//    BASE CLASS  : rg_BSplineSurface3D
//      
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun. 1996    
//
//    HISTORY     :
//          BY Young-Song Cho.  14 Aug. 1997
//            make : rg_REAL* evaluateMultiBasisFuncU( const rg_INDEX     &uKnotIndex,
//                                                  const rg_PARAMETER &u,
//                                                  const rg_ORDER     &uOrder)
//            make : rg_REAL* evaluateMultiBasisFuncV( const rg_INDEX     &vKnotIndex,
//                                                  const rg_PARAMETER &v,
//                                                  const rg_ORDER     &vOrder)
//            modify : rg_Point3D evaluatePt( const rg_PARAMETER &u, 
//                                       const rg_PARAMETER &v )
//          BY Young-Song Cho.   14 Oct. 1997
//            make : rg_REAL  getGaussianCurvature( const rg_PARAMETER &u, 
//                                               const rg_PARAMETER &v)
//
//          By Dong-Gyou Lee.	24 Mar. 1998
//                	void powerSplineToNUBSplineSurface( const rg_DEGREE& uDegree,
//													    const rg_DEGREE& vDegree,
//														const rg_INT& numOfUPatch,
//														const rg_INT& numOfVPatch,
//														rg_REAL** paramValuesOfU,
//														rg_REAL** paramValuesOfV,
//														rg_Matrix** powerCoeff )
//
//                  void bzSurfacesToNUBSplineSurface( const rg_INT& numOfUPatch,
//													   const rg_INT& numOfVPatch,
//													   rg_BzSurface3D** bzSurfaces )
//
//           By Dong-Gyou Lee 26 Mar. 1998
//                  void reparameterizationKnotVector() 
//
//			 By Joonghyun Ryu 6 Aug. 2001
//					void knotRefinement(const rg_REAL* &insertingKnotValuesInU, const rg_INT& numOfinsertingKnotValuesInU,
//		                                const rg_REAL* &insertingKnotValuesInV, const rg_INT& numOfinsertingKnotValuesInV)
//					rg_BzSurface3D** decomposeSurfaceIntoBezierPatchesUsingKnotRefinement()
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_NUBSPLINESURFACE3D_H
#define _RG_NUBSPLINESURFACE3D_H

#include "rg_BSplineSurface3D.h"

#include "rg_Const.h"
//#include "DefConst.h"
#include "rg_ListByPtr.h"
#include "rg_Point3D.h"
#include "rg_BzSurface3D.h"

#include "rg_NUBSplineCurve3D.h"

class rg_NUBSplineSurface3D : public rg_BSplineSurface3D
{
protected:
	rg_REAL* u_knot_vector;
	rg_REAL* v_knot_vector;

public:
////	Constructor & Destructor.----------------------------------------------
	rg_NUBSplineSurface3D();
	rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
		                const rg_Planarity    &newPlanarity );
	rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
		                const rg_Planarity    &newPlanarity, 
						const rg_INT          &row, 
						const rg_INT          &col ); 
	rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
		                const rg_Planarity    &newPlanarity, 
						const rg_ORDER        &uOrder, 
						const rg_ORDER        &vOrder );
	rg_NUBSplineSurface3D( const rg_INT &row, 
		                const rg_INT &col );
	rg_NUBSplineSurface3D( const rg_ORDER &uOrder, 
		                const rg_ORDER &vOrder );
	rg_NUBSplineSurface3D( const rg_INT   &row, 
		                const rg_INT   &col,
						const rg_ORDER &uOrder, 
						const rg_ORDER &vOrder );
	rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
		                const rg_Planarity    &newPlanarity, 
						const rg_INT          &row, 
						const rg_INT          &col, 
						const rg_ORDER        &uOrder, 
						const rg_ORDER        &vOrder );
    //// Constructor     : March 13 1997
    rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
		                const rg_Planarity    &newPlanarity, 
						const rg_INT          &row, 
						const rg_INT          &col, 
						const rg_ORDER        &uOrder, 
						const rg_ORDER        &vOrder,
                        rg_Point3D**             ctrlNet,
                        rg_REAL*               uKnotVector,
                        rg_REAL*               vKnotVector );
    //// Copy Constructor     : March 13 1997
    rg_NUBSplineSurface3D(const rg_NUBSplineSurface3D &surface);

    virtual ~rg_NUBSplineSurface3D();

////	Get Functions.---------------------------------------------------------
	rg_REAL  getKnotValueOfU( const rg_INDEX &kIndex ) const;	
    rg_REAL* getKnotVectorOfU() const;
	rg_REAL  getKnotValueOfV( const rg_INDEX &kIndex ) const;
    rg_REAL* getKnotVectorOfV() const;

	rg_INT getNumOfNonZeroLengthKnotSpanOfKnotVectorU() const;
	rg_INT getNumOfNonZeroLengthKnotSpanOfKnotVectorV() const;

	rg_REAL* getDistinctKnotValuesU() const;
	rg_REAL* getDistinctKnotValuesV() const;

    rg_NUBSplineCurve3D getUIsoparamCurve(const rg_PARAMETER &u);
    rg_NUBSplineCurve3D getVIsoparamCurve(const rg_PARAMETER &v);

    rg_Point3D getUnitNormalVector(const rg_PARAMETER &u, const rg_PARAMETER &v);
    rg_REAL  getGaussianCurvature(const rg_PARAMETER &u, const rg_PARAMETER &v);

////    Set Functions.---------------------------------------------------------
	void   setInitialKnotVectorOfU();
	void   setInitialKnotVectorOfV();
	void   setKnotValueOfU( const rg_INDEX &kIndex, 
		                    const rg_REAL  &knotValue );
	void   setKnotValueOfV( const rg_INDEX &kIndex, 
		                    const rg_REAL  &knotValue );
    void   setKnotVectorOfU( const rg_INT &numOfKnot,
                             const rg_REAL* const newKnotVector );
    void   setKnotVectorOfV( const rg_INT &numOfKnot, 
                             const rg_REAL* const newKnotVector );

////	Operating & Calculating.-----------------------------------------------
    virtual rg_REAL evaluateBasisFuncU( const rg_INDEX     &index, 
                                     const rg_PARAMETER &u,
                                     const rg_ORDER     &uOrder );
    virtual rg_REAL evaluateBasisFuncV( const rg_INDEX     &index, 
                                     const rg_PARAMETER &v,
                                     const rg_ORDER     &vOrder );

    rg_REAL* evaluateMultiBasisFuncU( const rg_INDEX     &uKnotIndex,
                                   const rg_PARAMETER &u,
                                   const rg_ORDER     &uOrder);
    rg_REAL* evaluateMultiBasisFuncV( const rg_INDEX     &vKnotIndex,
                                   const rg_PARAMETER &v,
                                   const rg_ORDER     &vOrder);

	virtual rg_Point3D evaluatePt( const rg_PARAMETER &u, 
                              const rg_PARAMETER &v );

////    Derivative.------------------------------------------------------------
    rg_NUBSplineSurface3D derivativeSurfaceOfU() const;
    rg_NUBSplineSurface3D derivativeSurfaceOfV() const;
    rg_NUBSplineSurface3D derivativeSurfaceOfUV() const;

////    Fundamental Geometric Algorithm----------------------------------------

	rg_INT findTheCorrespondingKnotSpanInU(const rg_REAL& insertingKnot);
	rg_INT findTheCorrespondingKnotSpanInV(const rg_REAL& insertingKnot);

	void knotRefinement(rg_REAL* &insertingKnotValuesInU, const rg_INT& numOfinsertingKnotValuesInU,
		                rg_REAL* &insertingKnotValuesInV, const rg_INT& numOfinsertingKnotValuesInV); // efficient surface-oriented knot refinement

	void knotRefinementOfU(rg_REAL* &insertingKnotValuesInU, const rg_INT& numOfinsertingKnotValuesInU);
	void knotRefinementOfV(rg_REAL* &insertingKnotValuesInV, const rg_INT& numOfinsertingKnotValuesInV);

	rg_BzSurface3D** decomposeSurfaceIntoBezierPatches() const;
	rg_BzSurface3D** decomposeSurfaceIntoBezierPatchesUsingKnotRefinement() const;


////	Conversion to power basis
	rg_REAL** getDistributionPolynomialInOneGraphPrimitive(rg_INT contributionKnotSpanInU, rg_INT contributionKnotSpanInV, rg_INT** graphInU, rg_INT** graphInV, rg_INT noOfAllPossiblePathsInU, rg_INT noOfAllPossiblePathsInV) const;
	void getPiecewiseSurfaceInPowerForm(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ) const;
	void makePiecewisePowerFormPolyUsingDE(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_REAL****** & truncatedBasisPolys) const;
	void makePiecewiseSurfaceInPowerFormUsingKR(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ) const;
	void makePiecewiseSurfaceInPowerFormUsingTE(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ);

	void makeUpdatedSurfaceInPowerBasisUsingDE(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_REAL****** & truncatedBasisPolys, rg_INT** indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts);
	void makeUpdatedSurfaceInPowerBasisUsingKR(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_INT** indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts);
	void makeUpdatedSurfaceInPowerBasisUsingTE(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_INT** indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts);


	void makeTruncatedBasisPolynomial(rg_REAL** & truncatedBasisPolynomial, rg_INT contributionKnotSpanInU, rg_INT contributionKnotSpanInV, rg_INT** graphInU, rg_INT** graphInV, rg_INT noOfAllPossiblePathsInU, rg_INT noOfAllPossiblePathsInV) const;


	
////	Skinned rg_Surface.-------------------------------------------------------

	void   makeSkinnedSurface( const rg_INT         &noOfSectionC, 
		                       rg_NUBSplineCurve3D* sectionCurves ); 
	void   setUOrdersOfSkinnedSurface( const rg_INT         &noOfSectionC,
									   rg_NUBSplineCurve3D* sectionCurves, 
									   rg_sListByPtr*            modifiedControlPolygons ); 
	void   setNewKnotVectorInSkinnedSurface( const rg_INT         &noOfSectionC, 
									         rg_NUBSplineCurve3D* sectionCurves, 
											 rg_INT               &noOfNewKnot );
	void   knotInsertionInSkinnedSurface ( rg_sListByPtr*     curKnotVector, 
		                                   rg_sListByPtr      &curControlPolygon, 
						                   const rg_REAL &insertingKnot );
	void   modifyControlPolygonsWithKnotInsertion( const rg_INT         &noOfNewKnot, 
		                                           const rg_INT         &noOfSectionC, 
									               rg_NUBSplineCurve3D* sectionCurves, 
												   rg_sListByPtr*            controlPolygons );
    void   setKnotVectorInVDirection( const rg_INT &noOfSectionC, 
                                      rg_sListByPtr*    controlPolygons );
	void   makeControlNetForSkinnedSurface( const rg_INT &noOfKnotVector, 
									        const rg_INT &noOfSectionC, 
											rg_sListByPtr*    controlPolygons );

////    Operator Overloading.--------------------------------------------------
    rg_NUBSplineSurface3D& operator =(const rg_NUBSplineSurface3D &surface);

////	File-Out Function.-----------------------------------------------------
	//void   fileOut( const char* fileName );   

////	Conversion between power spline & b-spline form.
	void powerSplineToNUBSplineSurface( const rg_DEGREE& uDegree,
									    const rg_DEGREE& vDegree,
										const rg_INT& numOfUPatch,
										const rg_INT& numOfVPatch,
										rg_REAL** paramValuesOfU,
										rg_REAL** paramValuesOfV,
										rg_Matrix** powerCoeff );

    void bzSurfacesToNUBSplineSurface( const rg_INT& numOfUPatch,
									   const rg_INT& numOfVPatch,
									   rg_BzSurface3D** bzSurfaces );

////	Reparameterization.
	void reparameterizationKnotVector(); 

};

#endif


