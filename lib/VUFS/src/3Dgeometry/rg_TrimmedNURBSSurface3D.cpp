//********************************************************************
//
//	  FILENAME    : rg_TrimmedNURBSSurface3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_TrimmedNURBSSurface3D
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

#include "rg_TrimmedNURBSSurface3D.h"


//	Constructor & Destructor.
rg_TrimmedNURBSSurface3D::rg_TrimmedNURBSSurface3D()
: rg_NURBSplineSurface3D(), rg_TrimmedSurface3D()
{
}

rg_TrimmedNURBSSurface3D::rg_TrimmedNURBSSurface3D( const rg_TrimmedNURBSSurface3D& surface )
: rg_NURBSplineSurface3D( surface.getID(),
                       surface.getPlanarity(),
                       surface.getRowOfControlNet(),
                       surface.getColumnOfControlNet(),
                       surface.getOrderOfU(),
                       surface.getOrderOfV(),
                       surface.getControlNet(),
                       surface.getKnotVectorOfU(),
                       surface.getKnotVectorOfV(),
                       surface.getWeightVector() ),
   rg_TrimmedSurface3D( surface.getBoundaries() )
{
}

rg_TrimmedNURBSSurface3D::rg_TrimmedNURBSSurface3D( const rg_NURBSplineSurface3D& surface )
: rg_NURBSplineSurface3D( surface.getID(),
                       surface.getPlanarity(),
                       surface.getRowOfControlNet(),
                       surface.getColumnOfControlNet(),
                       surface.getOrderOfU(),
                       surface.getOrderOfV(),
                       surface.getControlNet(),
                       surface.getKnotVectorOfU(),
                       surface.getKnotVectorOfV(),
                       surface.getWeightVector() )
{ 
}

rg_TrimmedNURBSSurface3D::~rg_TrimmedNURBSSurface3D()
{
}

//	Set functions.
void rg_TrimmedNURBSSurface3D::setBaseSurface( const rg_NURBSplineSurface3D& surface )
{
    rg_Surface::setID( surface.rg_Surface::getID() );
    rg_Surface::setPlanarity( surface.rg_Surface::getPlanarity() );

	rg_INT row = surface.rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT col = surface.rg_BSplineSurface3D::getColumnOfControlNet();
	
	rg_INT uOrder = surface.rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = surface.rg_BSplineSurface3D::getOrderOfV();
    
	rg_BSplineSurface3D::setControlNet( row, col,
									 surface.rg_BSplineSurface3D::getControlNet() );

    rg_BSplineSurface3D::setOrderOfSurface( uOrder, vOrder );

    rg_NUBSplineSurface3D::setKnotVectorOfU( row + uOrder,
		                                  surface.rg_NUBSplineSurface3D::getKnotVectorOfU() );
    rg_NUBSplineSurface3D::setKnotVectorOfV( col + vOrder,
		                                  surface.rg_NUBSplineSurface3D::getKnotVectorOfV() );

	rg_NURBSplineSurface3D::setWeightVector( surface.rg_NURBSplineSurface3D::getWeightVector() );
}

//	Get functions.
rg_NURBSplineSurface3D rg_TrimmedNURBSSurface3D::getBaseSurface() const
{
	rg_NURBSplineSurface3D baseSurface( getID(),
                                     getPlanarity(),
                                     getRowOfControlNet(),
                                     getColumnOfControlNet(),
                                     getOrderOfU(),
                                     getOrderOfV(),
                                     getControlNet(),
                                     getKnotVectorOfU(),
                                     getKnotVectorOfV(),
                                     getWeightVector() );

	return baseSurface;
}

//	Operations
bool rg_TrimmedNURBSSurface3D::selfTrimming()
{
	if( getRowOfControlNet() <= 0 || getColumnOfControlNet() <= 0 )
		return false;

	rg_Boundary3D surfaceBoundary;
	rg_NURBSplineCurve3D boundaryCurve[4];

	boundaryCurve[0] = getVIsoparamCurve( 0. );
	boundaryCurve[1] = getUIsoparamCurve( 1. );
	boundaryCurve[2] = getVIsoparamCurve( 1. );
	boundaryCurve[3] = getUIsoparamCurve( 0. );
	
	for( rg_INDEX i=0; i < 4; i++ )
		surfaceBoundary.addCurveOnObjectOfOuterLoop( boundaryCurve[i] );

	setOneBoundary( surfaceBoundary );

	return true;
}
	
//	Operator overloading.
rg_TrimmedNURBSSurface3D& rg_TrimmedNURBSSurface3D::operator=( const rg_TrimmedNURBSSurface3D& surface )
{
	if( this == &surface )
		return *this;

    rg_Surface::setID( surface.rg_Surface::getID() );
    rg_Surface::setPlanarity( surface.rg_Surface::getPlanarity() );

	rg_INT row = surface.rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT col = surface.rg_BSplineSurface3D::getColumnOfControlNet();
	
	rg_INT uOrder = surface.rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = surface.rg_BSplineSurface3D::getOrderOfV();
    
	rg_BSplineSurface3D::setControlNet( row, col,
									 surface.rg_BSplineSurface3D::getControlNet() );

    rg_BSplineSurface3D::setOrderOfSurface( uOrder, vOrder );

    rg_NUBSplineSurface3D::setKnotVectorOfU( row + uOrder,
		                                  surface.rg_NUBSplineSurface3D::getKnotVectorOfU() );
    rg_NUBSplineSurface3D::setKnotVectorOfV( col + vOrder,
		                                  surface.rg_NUBSplineSurface3D::getKnotVectorOfV() );

	rg_NURBSplineSurface3D::setWeightVector( surface.rg_NURBSplineSurface3D::getWeightVector() );

	rg_TrimmedSurface3D::setBoundaries( surface.rg_TrimmedSurface3D::getBoundaries() );
    
	return *this;
}



