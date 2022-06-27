//********************************************************************
//
//	  FILENAME    : rg_Boundary3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_Boundary3D
//           which define a boundary of surface and its property. 
//                          
//	  CLASS NAME  : rg_Boundary3D
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

#ifndef _RG_BOUNDARY_H
#define _RG_BOUNDARY_H

#include "rg_dList.h"
#include "rg_Loop3D.h"

class rg_Boundary3D
{
private:
    rg_Loop3D        outerLoop;
	rg_dList<rg_Loop3D> innerLoop;

public:

//	Constructor & Destructor.
	rg_Boundary3D();
	rg_Boundary3D( const rg_Boundary3D& boundary );
	rg_Boundary3D( const rg_Loop3D& outer, const rg_dList<rg_Loop3D>& inners );

	virtual ~rg_Boundary3D();

//	Set functions.
	void setBoundary( const rg_Loop3D& outer, const rg_dList<rg_Loop3D>& inners );
	void setBoundary( const rg_Loop3D& outer, const rg_Loop3D& inner );

	void setOuterLoop( const rg_Loop3D& outer );
	
	void setInnerLoop( const rg_dList<rg_Loop3D>& inners );
	void setOneInnerLoop( const rg_Loop3D& inner );

//	Get functions.
	rg_Loop3D        getOuterLoop() const;
	rg_dList<rg_Loop3D> getInnerLoop() const;
	rg_Loop3D        getOneInnerLoop( const rg_INDEX& index ) const;

	rg_INT getNumOfInnerLoop() const;

//	Add loop.
	void addCurveOnDomainOfOuterLoop( const rg_NURBSplineCurve3D& domainCurve );
	void addCurveOnObjectOfOuterLoop( const rg_NURBSplineCurve3D& objectCurve );

	void addInnerLoop( const rg_dList<rg_Loop3D>& inners );
	void addOneInnerLoop( const rg_Loop3D& inner );

	bool addCurveOnDomainOfInnerLoop( const rg_INT& index, const rg_NURBSplineCurve3D& domainCurve );
	bool addCurveOnObjectOfInnerLoop( const rg_INT& index, const rg_NURBSplineCurve3D& objectCurve );

// Operator overloading.
	rg_Boundary3D& operator=( const rg_Boundary3D& boundary );
};

#endif


