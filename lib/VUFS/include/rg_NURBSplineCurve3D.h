//******************************************************************************
//
//	  FILENAME    : rg_NURBSplineCurve3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_NURBSplineCurve3D 
//           which define Non-uniform Rational B-Spline rg_Curve and its property. 
//                          
//	  CLASS NAME  : rg_NURBSplineCurve3D
//
//    BASE CLASS  : rg_NUBSplineCurve3D
//      
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//
//    HISTORY     :   
//         1. By Tae-Bum Jang.    1997.11.19. 
//                  void  formArc ( const rg_Point3D& center,
//                                  const rg_Point3D& start,
//                                  const rg_Point3D& end,
//                                  const rg_Point3D& localZ  )
//
//                  void  formLine( const rg_Point3D& start,
//                                  const rg_Point3D& end   )
//
//         2. By Tae-Bum Jang.    1997. 8.20. 
//                  void setCurve(const rg_RBzCurve3D& curve);
//
//         3. By Taeboom Jang   1998. 9. 7
//                  void formPolygon( rg_Point3D* polygons, const rg_INT& numOfPts);
//
//		   4. By JoongHyun Ryu  1998. 10.
//					rg_RBzCurve3D* decomposeCurveIntoBezierSegment() const;
//
//    START DATE  : 21 Jun 1996 
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************

#ifndef _RG_NURBSPLINECURVE3D_H
#define _RG_NURBSPLINECURVE3D_H

#include "rg_NUBSplineCurve3D.h"
#include "rg_RBzCurve3D.h"
#include "rg_RationalPolynomial.h"

//#include "DefConst.h"
#include "rg_Point3D.h"
#include "rg_EllipticArc3D.h"
#include "rg_Polyline3D.h"
#include "rg_Polyline2D.h"

class rg_NURBSplineCurve3D : public rg_NUBSplineCurve3D
{
protected:
    rg_REAL* weights;

public:
////    Constructor & Destructor
    rg_NURBSplineCurve3D();
    rg_NURBSplineCurve3D( const rg_ORDER &newOrder );
    //// Constructor      : March 13 1997
    rg_NURBSplineCurve3D( const unsigned rg_INT &newID, 
                       const rg_Planarity    &newPlanarity,
                       const rg_ORDER        &newOrder, 
                       const rg_INT          &num, 
                       rg_Point3D*              newControlP,
                       rg_REAL*               newKnotVector,
                       rg_REAL*               weight_vector );
    rg_NURBSplineCurve3D( const rg_BSplineCurve3D &curve );
    rg_NURBSplineCurve3D( const rg_NUBSplineCurve3D &curve );
    //// Copy Constructor : March 13 1997
    rg_NURBSplineCurve3D( const rg_NURBSplineCurve3D &curve );
    
    virtual ~rg_NURBSplineCurve3D();

////    Get Functions.
    //  March 13 1997
    rg_REAL  getWeight(const rg_INT &i) const;
    rg_REAL* getWeightVector() const;
	rg_RationalPolynomial** makePolynomialFormCurveUsingCurveDecomp() const;

// B-spline basis to Power basis
	void getEntirePolynomialCurve(rg_REAL** & coeffOfXCoordinate, rg_REAL** & coeffOfYCoordinate, rg_REAL** & coeffOfZCoordinate, rg_REAL** & coeffOfWeight) const;

    rg_NUBSplineCurve3D getNumeratorInNUBS() const;
    rg_NUBSplineCurve3D getDenominatorInNUBS() const; // Only x coordinate has meaning 
                                                   // others coordinates are setted to 0

    rg_REAL             getParameterOfNearestPt(const rg_Point3D& pt,
                                             const rg_REAL&    distanceTol,
                                             const rg_REAL&    cosAngleTol) const;
    rg_REAL             getParameterOfNearestPt(const rg_Point3D& pt,
                                             const rg_REAL&    distanceTol,
                                             const rg_REAL&    cosAngleTol,
                                             const rg_REAL&    seed) const;

    rg_INT              getIndexOfNearestCtrlPt(const rg_Point3D& pt) const;
////    Set Functions.
    //  March 13 1997
    void setWeight(const rg_INT &i, const rg_REAL &weight);
    void setWeightVector(rg_REAL* weightVector);
    void setNURBSCurve(const unsigned rg_INT &newID, 
                       const rg_Planarity    &newPlanarity,
                       const unsigned rg_INT &newOrder, 
                       const rg_INT          &num, 
                       rg_Point3D*             newControlP,
                       const rg_REAL* const  newKnotVector,
                       const rg_REAL* const  newWeightVector );
    void setCurve(const rg_RBzCurve3D& curve);
	void setNumOfCtrlPts(const rg_INT& newNumOfCtrlpts);
////    BasisFunction & rg_Point3D Evaluating.
    //  March 13 1997
    virtual rg_Point3D  evaluatePt(const rg_REAL &param) const;
    virtual rg_Point3D* evaluatePtsInEvenParameter(const rg_INT  &numOfPtOnCurve) const;
    virtual rg_Polyline2D makePolyline2DInEvenParameter(const rg_INT &numOfPts) const;
    virtual rg_Polyline3D makePolyline3DInEvenParameter(const rg_INT &numOfPts) const;
	virtual rg_Polyline2D makePolyline2DConsideringKnots(const rg_INT &numOfPts) const;

////    Derivative rg_Curve & Curvature
    rg_Point3D  makeDerivative(const rg_REAL &u) const;
	rg_NURBSplineCurve3D makeHodograph() const;
	rg_NURBSplineCurve3D makeRedundantHodograph() const;

////    represent the basic geometry with NURBS
    void    formArc( const rg_Point3D& center,
                     const rg_Point3D& start,
                     const rg_Point3D& end,
                     const rg_Point3D& localZ  );

    void    formLine( const rg_Point3D& start,
                      const rg_Point3D& end   );

    void    formPolyline3D( const rg_INT      &numOfPts, 
                                  rg_Point3D*  polygons);
    void    formEllipticArc( const rg_EllipticArc3D& ellipticArc);
////    Fundamental Geometric Algorithm.
    void reverseTrace();
    void curveInterpolation(const rg_INT          &n, 
                            const rg_Point3D* const  ptsPassedThrough,
                            const rg_INT          &order = 4,
                            rg_REAL*               param = rg_NULL,
                            rg_REAL*               weightVector = rg_NULL);

    void knotInsertion( const rg_REAL &insertingKnot );

	rg_RBzCurve3D* decomposeCurveIntoBezierSegment() const;
	rg_NURBSplineCurve3D decomposeCurveIntoBezierSegmentInNURBS() const;

    
    rg_NURBSplineCurve3D evaluateCurveSegment( const rg_REAL& start,
                                       const rg_REAL& end ) const;
    void knotInsertion( const rg_REAL& insertingKnot,
                        const rg_INT&  insertinTimes);

	void removeRedundantKnot(const rg_REAL& knot); // remove knot which must be reduant.
	rg_REAL  findParameterOfNearestPt(const rg_Point3D& pt,
                                     const rg_REAL&    distanceTol,
                                     const rg_REAL&    cosAngleTol,
			   					     const rg_REAL&    start,
									 const rg_REAL&    end) const;
	rg_REAL  findParameterOfNearestPt(const rg_Point3D& pt,
                                      const rg_REAL&    distanceTol,
                                      const rg_REAL&    cosAngleTol,
                                      const rg_REAL&    seed,
				    				  const rg_REAL&    start,
					    			  const rg_REAL&    end) const;

////    Operator Overloading.
    rg_NURBSplineCurve3D& operator =(const rg_NURBSplineCurve3D &curve);

    void removeAll();
};


#endif


