//********************************************************************
//
//	  FILENAME    : rg_Boundary3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Boundary3D
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

#include "rg_Boundary3D.h"
#include "rg_dList.h"

//	Constructor & Destructor.
rg_Boundary3D::rg_Boundary3D()
{
}

rg_Boundary3D::rg_Boundary3D( const rg_Boundary3D& boundary )
{
	outerLoop = boundary.outerLoop;

	if( boundary.innerLoop.getSize() )
		innerLoop = boundary.innerLoop;
}

rg_Boundary3D::rg_Boundary3D( const rg_Loop3D& outer, const rg_dList<rg_Loop3D>& inners )
{
	outerLoop = outer;

	if( inners.getSize() )
		innerLoop = inners;
}

rg_Boundary3D::~rg_Boundary3D()
{
	if( innerLoop.getSize() )
		innerLoop.removeAll();
}

//	Set functions.
void rg_Boundary3D::setBoundary( const rg_Loop3D& outer, const rg_dList<rg_Loop3D>& inners )
{
	outerLoop = outer;

	if( innerLoop.getSize() )
		innerLoop.removeAll();

	innerLoop = inners;
}

void rg_Boundary3D::setBoundary( const rg_Loop3D& outer, const rg_Loop3D& inner )
{
	outerLoop = outer;

	if( innerLoop.getSize() )
		innerLoop.removeAll();
	
	innerLoop.add( inner );
}

void rg_Boundary3D::setOuterLoop( const rg_Loop3D& outer )
{
	outerLoop = outer;
}

void rg_Boundary3D::setInnerLoop( const rg_dList<rg_Loop3D>& inners )
{
	if( innerLoop.getSize() )
		innerLoop.removeAll();
	
	innerLoop = inners;
}

void rg_Boundary3D::setOneInnerLoop( const rg_Loop3D& inner )
{
	if( innerLoop.getSize() )
		innerLoop.removeAll();
	
	innerLoop.add( inner );
}

//	Get functions.
rg_Loop3D rg_Boundary3D::getOuterLoop() const
{
	return outerLoop;
}

rg_dList<rg_Loop3D> rg_Boundary3D::getInnerLoop() const
{
	return innerLoop;
}

rg_Loop3D rg_Boundary3D::getOneInnerLoop( const rg_INDEX& index ) const
{
	return innerLoop.elementAt( index );
}

rg_INT rg_Boundary3D::getNumOfInnerLoop() const
{
	return innerLoop.getSize();
}

//	Add loop.
void rg_Boundary3D::addCurveOnDomainOfOuterLoop( const rg_NURBSplineCurve3D& domainCurve )
{
	outerLoop.addOneCurveOnDomain( domainCurve );
}

void rg_Boundary3D::addCurveOnObjectOfOuterLoop( const rg_NURBSplineCurve3D& objectCurve )
{
	outerLoop.addOneCurveOnObject( objectCurve );
}

void rg_Boundary3D::addInnerLoop( const rg_dList<rg_Loop3D>& inners )
{
	innerLoop.append( inners );
}

void rg_Boundary3D::addOneInnerLoop( const rg_Loop3D& inner )
{
	innerLoop.add( inner );
}

bool rg_Boundary3D::addCurveOnDomainOfInnerLoop( const rg_INT& index, const rg_NURBSplineCurve3D& domainCurve )
{
	if( innerLoop.getSize() < index )
		return false;

	innerLoop[index].addOneCurveOnDomain( domainCurve );

	return true;
}

bool rg_Boundary3D::addCurveOnObjectOfInnerLoop( const rg_INT& index, const rg_NURBSplineCurve3D& objectCurve )
{
	if( innerLoop.getSize() < index )
		return false;

	innerLoop[index].addOneCurveOnObject( objectCurve );

	return true;

}


// Operator overloading.
rg_Boundary3D& rg_Boundary3D::operator=( const rg_Boundary3D& boundary )
{
	if( this == &boundary )
		return *this;

	outerLoop = boundary.outerLoop;

	if( innerLoop.getSize() )
		innerLoop.removeAll();

	if( boundary.innerLoop.getSize() )
		innerLoop = boundary.innerLoop;

	return *this;
}



