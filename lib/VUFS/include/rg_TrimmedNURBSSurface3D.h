//********************************************************************
//
//	  FILENAME    : rg_TrimmedNURBSSurface3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_TrimmedNURBSSurface3D
//           which define a Trimmed NURBS rg_Surface and its property. 
//                          
//	  CLASS NAME  : rg_TrimmedNURBSSurface3D
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

#ifndef _RG_TRIMMEDNURBSSURFACE3D_H
#define _RG_TRIMMEDNURBSSURFACE3D_H

#include "rg_NURBSplineSurface3D.h"
#include "rg_TrimmedSurface3D.h"

class rg_TrimmedNURBSSurface3D : public rg_NURBSplineSurface3D, public rg_TrimmedSurface3D
{
private:

public:

//	Constructor & Destructor.
	rg_TrimmedNURBSSurface3D();
	rg_TrimmedNURBSSurface3D( const rg_TrimmedNURBSSurface3D& surface );
	rg_TrimmedNURBSSurface3D( const rg_NURBSplineSurface3D& surface );

	virtual ~rg_TrimmedNURBSSurface3D();

//	Set functions.
	void setBaseSurface( const rg_NURBSplineSurface3D& surface );
	
//	Get functions.
	rg_NURBSplineSurface3D getBaseSurface() const;

//	Operations
	bool selfTrimming();

//	Operator overloading.
	rg_TrimmedNURBSSurface3D& operator=( const rg_TrimmedNURBSSurface3D& surface );

};

#endif


