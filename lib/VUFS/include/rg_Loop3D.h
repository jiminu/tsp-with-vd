//********************************************************************
//
//	  FILENAME    : rg_Loop3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_Loop3D
//           which define a loop of curves and its property. 
//                          
//	  CLASS NAME  : rg_Loop3D
//
//    BASE CLASS  : None
//      
//    AUTHOR      : Deok-Soo Kim, Dong-Gyou Lee
//
//    HISTORY     : 	
//
//    START DATE  : 3 Apr. 1998    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_LOOP3D_H
#define _RG_LOOP3D_H

#include "rg_dList.h"
#include "rg_NURBSplineCurve3D.h"

class rg_Loop3D
{
private:
    rg_dList<rg_NURBSplineCurve3D> curvesOnDomain;
	rg_dList<rg_NURBSplineCurve3D> curvesOnObject;

public:

//	Constructor & Destructor.
	rg_Loop3D();
	rg_Loop3D( const rg_Loop3D& loop );
	rg_Loop3D( const rg_dList<rg_NURBSplineCurve3D>& domainCurves,
		  const rg_dList<rg_NURBSplineCurve3D>& objectCurves );

	virtual ~rg_Loop3D();

//	Set functions.
	void setLoop( const rg_dList<rg_NURBSplineCurve3D>& domainCurves,
		          const rg_dList<rg_NURBSplineCurve3D>& objectCurves );
	
	void setCurvesOnDomain( const rg_dList<rg_NURBSplineCurve3D>& domainCurves );
	void setCurvesOnObject( const rg_dList<rg_NURBSplineCurve3D>& objectCurves );

	void setOneCurveOnDomain( const rg_NURBSplineCurve3D& domainCurve );
	void setOneCurveOnObject( const rg_NURBSplineCurve3D& objectCurve );
	
//	Get functions.
	rg_dList<rg_NURBSplineCurve3D> getCurvesOnDomain() const;
	rg_dList<rg_NURBSplineCurve3D> getCurvesOnObject() const;

	rg_INT getNumOfDomainCurves() const;
	rg_INT getNumOfObjectCurves() const;

//	Add curves.
	void addCurvesOnDomain( const rg_dList<rg_NURBSplineCurve3D>& domainCurves );
	void addCurvesOnObject( const rg_dList<rg_NURBSplineCurve3D>& objectCurves );

	void addOneCurveOnDomain( const rg_NURBSplineCurve3D& domainCurve );
	void addOneCurveOnObject( const rg_NURBSplineCurve3D& objectCurve );

//	Operator overloading
	rg_Loop3D& operator=( const rg_Loop3D& loop ); 
};

#endif


