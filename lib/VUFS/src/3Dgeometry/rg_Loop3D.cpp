//********************************************************************
//
//	  FILENAME    : rg_Loop3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Loop3D
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

#include "rg_Loop3D.h"
#include "rg_dList.h"

//	Constructor & Destructor.
rg_Loop3D::rg_Loop3D()
{
}

rg_Loop3D::rg_Loop3D( const rg_Loop3D& loop )
{
	if( loop.curvesOnDomain.getSize() )
		curvesOnDomain = loop.curvesOnDomain;
	
	if( loop.curvesOnObject.getSize() )
		curvesOnObject = loop.curvesOnObject;

}

rg_Loop3D::rg_Loop3D( const rg_dList<rg_NURBSplineCurve3D>& domainCurves, 
		    const rg_dList<rg_NURBSplineCurve3D>& objectCurves )
{
	curvesOnDomain = domainCurves;
	curvesOnObject = objectCurves;
}

rg_Loop3D::~rg_Loop3D()
{
	if( curvesOnDomain.getSize() )
		curvesOnDomain.removeAll();

	if( curvesOnObject.getSize() )
		curvesOnObject.removeAll();
}

//	Set functions.
void rg_Loop3D::setLoop( const rg_dList<rg_NURBSplineCurve3D>& domainCurves,
				    const rg_dList<rg_NURBSplineCurve3D>& objectCurves )
{
	if( curvesOnDomain.getSize() )
		curvesOnDomain.removeAll();

	if( curvesOnObject.getSize() )
		curvesOnObject.removeAll();

	curvesOnDomain = domainCurves;
	curvesOnObject = objectCurves;
}

void rg_Loop3D::setCurvesOnDomain( const rg_dList<rg_NURBSplineCurve3D>& domainCurves )
{
	if( curvesOnDomain.getSize() )
		curvesOnDomain.removeAll();

	curvesOnDomain = domainCurves;
}

void rg_Loop3D::setCurvesOnObject( const rg_dList<rg_NURBSplineCurve3D>& objectCurves )
{
	if( curvesOnObject.getSize() )
		curvesOnObject.removeAll();

	curvesOnObject = objectCurves;
}

void rg_Loop3D::setOneCurveOnDomain( const rg_NURBSplineCurve3D& domainCurve )
{
	if( curvesOnDomain.getSize() )
		curvesOnDomain.removeAll();

	curvesOnDomain.add( domainCurve );
}

void rg_Loop3D::setOneCurveOnObject( const rg_NURBSplineCurve3D& objectCurve )
{
	if( curvesOnObject.getSize() )
		curvesOnObject.removeAll();

	curvesOnObject.add( objectCurve );
}

//	Get functions.
rg_dList<rg_NURBSplineCurve3D> rg_Loop3D::getCurvesOnDomain() const
{
	return curvesOnDomain;
}

rg_dList<rg_NURBSplineCurve3D> rg_Loop3D::getCurvesOnObject() const
{
	return curvesOnObject;
}

rg_INT rg_Loop3D::getNumOfDomainCurves() const
{
	return curvesOnDomain.getSize();
}

rg_INT rg_Loop3D::getNumOfObjectCurves() const
{
	return curvesOnObject.getSize();
}

//	Add curves.
void rg_Loop3D::addCurvesOnDomain( const rg_dList<rg_NURBSplineCurve3D>& domainCurves )
{
	curvesOnDomain.append( domainCurves );
}

void rg_Loop3D::addCurvesOnObject( const rg_dList<rg_NURBSplineCurve3D>& objectCurves )
{
	curvesOnObject.append( objectCurves );
}

void rg_Loop3D::addOneCurveOnDomain( const rg_NURBSplineCurve3D& domainCurve )
{
	curvesOnDomain.add( domainCurve );
}

void rg_Loop3D::addOneCurveOnObject( const rg_NURBSplineCurve3D& objectCurve )
{
	curvesOnObject.add( objectCurve );
}

//	Operator overloading
rg_Loop3D& rg_Loop3D::operator=( const rg_Loop3D& loop ) 
{
	if( this == &loop )
		return *this;

	if( curvesOnDomain.getSize() )
		curvesOnDomain.removeAll();

	if( curvesOnObject.getSize() )
		curvesOnObject.removeAll();

	if( loop.curvesOnDomain.getSize() )
		curvesOnDomain = loop.curvesOnDomain;
	
	if( loop.curvesOnObject.getSize() )
		curvesOnObject = loop.curvesOnObject;

	return *this;
}








