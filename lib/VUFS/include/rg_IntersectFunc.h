/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_IntersectFunc.h
//	  
//    DESCRIPTION : collections of function related with intersection
//                          
//	  CLASS NAME  : 
//
//    BASE CLASS  : 
//      
//    AUTHOR      : Lee, Soon-Woong
//
//    START DATE  :     
//    
//    HISTORY :     1. Ryu, Jung-Hyun inserted the following functions on Feb. 25. 1998
//
//                       rg_ImplicitEquation getResultantByBezout(rg_Polynomial x_t, rg_Polynomial y_t);
//                       rg_ImplicitEquation getResultantByBezout(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t);
//					     rg_ImplicitEquation implicitizeNotUsingMathematica(const rg_RQBzCurve2D &curve);
//                       rg_ImplicitEquation implicitizeNotUsingMathematica(const rg_BzCurve2D& curve);
//
//                  2. Ryu, Jung-Hyun inserted the following functions on Mar. 5. 1998
//
//                       rg_ImplicitEquation getResultantBySylvester(rg_Polynomial x_t, rg_Polynomial y_t);
//                       rg_ImplicitEquation getResultantBySylvester(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t);
//
//
//					 
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _INTERSECT_included
#define _INTERSECT_included

#include "rg_RQBzCurve2D.h"
#include "rg_QuarticPolynomial.h"
#include "rg_Point2D.h"
#include "rg_Point3D.h"
#include "rg_Line.h"
#include "rg_Line2D.h"
#include "rg_Line3D.h"
#include "rg_QBzCurve2D.h"
#include "rg_CBzCurve2D.h"
#include "rg_dList.h"
#include "rg_ImplicitEquation.h"
#include "rg_Matrix.h"

#include "rg_Const.h"

#include <list>
using namespace std;



enum TypeLineLineIntersection { PPI_NO_INTERSECTION, PPI_ONE_INTERSECTION, PPI_ALL_INTERSECTION};

class rg_IntersectFunc
{
public:

static rg_ImplicitEquation          implicitize(const rg_RQBzCurve2D &curve);

// Made by Ryu, Jung-Hyun 

static rg_ImplicitEquation          getResultantByBezout(rg_Polynomial x_t, rg_Polynomial y_t);// test
static rg_ImplicitEquation          getResultantBySylvesterOld(rg_Polynomial x_t, rg_Polynomial y_t);//antique
static rg_ImplicitEquation          getResultantBySylvester(rg_Polynomial x_t, rg_Polynomial y_t);

//rg_ImplicitEquation          getResultant(const rg_BzCurve2D& curve);
static rg_ImplicitEquation          getResultantByBezout(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t);
static rg_ImplicitEquation          getResultantBySylvesterOld(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t);//antique
static rg_ImplicitEquation          getResultantBySylvester(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t);

static rg_ImplicitEquation          implicitizeNotUsingMathematica1(const rg_RQBzCurve2D &curve);
static rg_ImplicitEquation          implicitizeNotUsingMathematica2(const rg_RQBzCurve2D &curve);

static rg_ImplicitEquation          implicitizeNotUsingMathematica1(const rg_BzCurve2D& curve);
static rg_ImplicitEquation          implicitizeNotUsingMathematica2(const rg_BzCurve2D& curve);

static rg_Matrix getNumericValueFromOperationBtnTwoDifferentPolynomial(const rg_Polynomial& f, const rg_Polynomial& g);
static rg_Matrix* getCoefficientMatrix(rg_Polynomial x_t, rg_Polynomial y_t, rg_Polynomial w_t);
static rg_Matrix* getCoefficientMatrix(rg_Polynomial x_t, rg_Polynomial y_t);

//

static rg_QuarticPolynomial          substituteParametricIntoImplicit(const rg_RQBzCurve2D &curve, const rg_ImplicitEquation &implicit);

static rg_Point2D            *intersectQBezierVsLine( rg_QBzCurve2D &curve, const rg_Line<rg_Point2D> &line);
static rg_Point2D             intersectLineVsLine(const rg_Line<rg_Point2D> &pLine, const rg_Line<rg_Point2D> &nLine);
static void                intersectRQBzCurveVsRQBzCurve(const rg_RQBzCurve2D &curve1, const rg_RQBzCurve2D &curve2, 
                                                  rg_dList <rg_Point2D> &intersectList);
static rg_INT              intersectRQBzCurveVsLine(const rg_RQBzCurve2D &rqbzCurve, const rg_Line2D& line, list<rg_Point2D>& intersetionPoints);
static rg_Point3D*			intersectUnboundedLine3Ds(const rg_Line3D& line1, rg_Line3D& line2);


static void intersectRQBzCurveVsRQBzCurve(const rg_RQBzCurve2D &curve1, const rg_REAL &s0, const rg_REAL &s1,  
										  const rg_RQBzCurve2D &curve2, const rg_REAL &t0, const rg_REAL &t1, 
										  rg_dList <rg_Point2D> &intersectList, 
										  rg_dList <rg_REAL*> &cubicParam);
//rg_dList<rg_Point2D> intersectCBzCurveVsCBzCurve(const rg_CBzCurve2D &curve1, const rg_CBzCurve2D &curve2);
//rg_RQBzCurve2D           approximateCBzCurve2RQBzCurve(rg_CBzCurve2D &curve);
};
#endif

