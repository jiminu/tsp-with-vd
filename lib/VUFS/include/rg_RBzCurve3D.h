//********************************************************************
//
//	  FILENAME    : rg_RBzCurve3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_BzCurve3D
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

#ifndef _RBZCURVE3D_H
#define _RBZCURVE3D_H

#include <math.h>
#include "rg_Const.h"
#include "rg_Point3D.h"
#include "rg_Point3D.h"
#include "rg_BzCurve3D.h"

class rg_RBzCurve3D :public rg_BzCurve3D
{
protected:
	rg_REAL* weight;

public:
	rg_RBzCurve3D();
	rg_RBzCurve3D(const rg_DEGREE &dgr);
	rg_RBzCurve3D(const rg_DEGREE &dgr, const rg_Point3D* ctrlPts);
	rg_RBzCurve3D(const rg_DEGREE &dgr, const rg_Point3D* ctrlPts, const rg_REAL* wght);
	rg_RBzCurve3D(const rg_RBzCurve3D &curve);
	virtual ~rg_RBzCurve3D();

	//Operations
	rg_Point3D evaluatePt(const rg_PARAMETER &t);
	rg_Point3D evaluateDerivative(const rg_PARAMETER &t);

	//Access elements
	rg_REAL    getWeight(const rg_INDEX &i) const;
    rg_REAL*   getWeight() const;

    void    setDegree(const rg_DEGREE& tDegree);
    void    setWeight(const rg_INDEX &i, const rg_REAL &w);
    void    setWeight(const rg_REAL *w);
	void    setCurve(const rg_RBzCurve3D &curve);

    //  Intersection---------------------------------------------------------------

    //  numOfIntersect : number of intersection point of Rational Quadratic 
    //                   Bezier curve and plane.
    //                   It is determined in this function. 
    rg_Point3D* intersectRationalQuadraticBezierAndPlane(
                const rg_Plane3D &plane, rg_INT &numOfIntersect);

	rg_Polynomial* convertNumerator2Polynomial() const;
	rg_Polynomial  convertDenominator2Polynomial() const;

	// temporary member function for algrithm for computing between planar curves
	// to be deleted !!
	rg_Point3D* inflectionPointByLeeMethod(rg_ComplexNumber* &solution, rg_INT& count);

    rg_RBzCurve3D getDerivative() const;
	rg_RBzCurve3D makeHodograph() const;
	rg_BzCurve3D makeScaledHodograph() const;
    // make conic
    void  makeConic( const rg_Point3D &startPt,
                     const rg_Point3D &startTangent,
                     const rg_Point3D &endPt,
                     const rg_Point3D &endTangent,
                     const rg_Point3D &passingPt );

	rg_RBzCurve3D& operator=(const rg_RBzCurve3D& temp);

	// Conversion to power basis

    rg_REAL computeLength(const rg_INT& numSamplingPoint);
};

#endif  // rg_RBzCurve3D class


