//********************************************************************
//
//	  FILENAME    : rg_TrimmedSurface3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_TrimmedSurface3D
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

#include "rg_TrimmedSurface3D.h"
#include "rg_dList.h"


//	Constructor & Destructor.
rg_TrimmedSurface3D::rg_TrimmedSurface3D()
{
}

rg_TrimmedSurface3D::rg_TrimmedSurface3D( const rg_TrimmedSurface3D& surface )
{
	boundaries = surface.boundaries;
}

rg_TrimmedSurface3D::rg_TrimmedSurface3D( const rg_dList<rg_Boundary3D>& boundaryList )
{
	if( boundaryList.getSize() )
		boundaries = boundaryList;
}

rg_TrimmedSurface3D::~rg_TrimmedSurface3D()
{
	if( boundaries.getSize() )
		boundaries.removeAll();
}

//	Set functions.
void rg_TrimmedSurface3D::setBoundaries( const rg_dList<rg_Boundary3D>& boundaryList )
{
	if( boundaryList.getSize() )
		boundaries = boundaryList;
}

void rg_TrimmedSurface3D::setOneBoundary( const rg_Boundary3D& boundary )
{
	if( boundaries.getSize() )
		boundaries.removeAll();

	boundaries.add( boundary );
}

//	Get functions.
rg_dList<rg_Boundary3D> rg_TrimmedSurface3D::getBoundaries() const
{
	return boundaries;
}

rg_INT rg_TrimmedSurface3D::getNumOfBoundaries() const
{
	return boundaries.getSize();
}

rg_Boundary3D rg_TrimmedSurface3D::getOneBoundary( const rg_INDEX& index ) const
{
	return boundaries.elementAt( index );
}


//	Operator overloading.
rg_TrimmedSurface3D& rg_TrimmedSurface3D::operator=( const rg_TrimmedSurface3D& surface )
{
	if( this == &surface )
		return *this;

	if( boundaries.getSize() )
		boundaries.removeAll();

	boundaries = surface.boundaries;

	return *this;
}



