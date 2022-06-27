//********************************************************************
//
//	  FILENAME    : rg_TrimmedSurface3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_TrimmedSurface3D
//           which define a Trimmed rg_Surface and its property. 
//                          
//	  CLASS NAME  : rg_TrimmedSurface3D
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

#ifndef _RG_TRIMMEDSURFACE3D_H
#define _RG_TRIMMEDSURFACE3D_H

#include "rg_dList.h"
#include "rg_Boundary3D.h"

class rg_TrimmedSurface3D
{
private:
    rg_dList<rg_Boundary3D> boundaries;

public:

//	Constructor & Destructor.
	rg_TrimmedSurface3D();
	rg_TrimmedSurface3D( const rg_TrimmedSurface3D& surface );
	rg_TrimmedSurface3D( const rg_dList<rg_Boundary3D>& boundaryList );

	virtual ~rg_TrimmedSurface3D();

//	Set functions.
	void setBoundaries( const rg_dList<rg_Boundary3D>& boundaryList );
	void setOneBoundary( const rg_Boundary3D& boundary );

//	Get functions.
	rg_dList<rg_Boundary3D> getBoundaries() const;
	rg_Boundary3D        getOneBoundary( const rg_INDEX& index ) const;
	
	rg_INT getNumOfBoundaries() const;

//	Operator overloading.
	rg_TrimmedSurface3D& operator=( const rg_TrimmedSurface3D& surface );

};

#endif

