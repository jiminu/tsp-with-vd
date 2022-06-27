//********************************************************************
//
//	  FILENAME    : rg_CurveSurfaceFunc.h
//	  
//    DESCRIPTION : 
//           This is the interface of external functions for curve and surface
//                          
//
//    AUTHOR      : Deok-Soo Kim, Dong-Gyou Lee Lee
//    START DATE  : 17 Mar. 1998    
//
//    HISTORY     :
//				
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _CURVESURFACEFUNC_INCLUDED
#define _CURVESURFACEFUNC_INCLUDED

#include "rg_Point3D.h"
#include "rg_Matrix.h"
#include "rg_RPolynomialCurve2D.h"
#include "rg_RBzCurve2D.h"
#include "rg_RBzCurve3D.h"

class rg_CurveSurfaceFunc
{
public:
////  For Power Basis rg_Curve & rg_Surface
static rg_Matrix bezierToPowerMatrix( const rg_INT& order );
static rg_Matrix powerToBezierMatrix( const rg_INT& order );
static rg_RPolynomialCurve2D convertToRPolynomialCurve2D( const rg_RBzCurve2D &curve);


////  For Reparameterization
static rg_Matrix reparameterMatrix( const rg_INT& order, const rg_REAL& multiValue, const rg_REAL& addValue ); 

};
////  For make curves and surfaces

#endif


