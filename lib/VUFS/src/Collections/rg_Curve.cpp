//********************************************************************
//
//	  FILENAME    : rg_Curve.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Curve 
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun 1996    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#include "rg_Curve.h"

////	Constructor & Destructor.
rg_Curve::rg_Curve()
: ID(0), planarOrSpatial(rg_PLANAR)
{
}

rg_Curve::rg_Curve( const unsigned rg_INT &newID, 
			  const rg_Planarity &newPlanarity)
: ID(newID), planarOrSpatial(newPlanarity)
{
}

////  Copy Constructor March 13 1997
rg_Curve::rg_Curve( const rg_Curve &curve )
: ID(curve.ID), planarOrSpatial(curve.planarOrSpatial)
{
}

rg_Curve::~rg_Curve()
{
}

////	Get Functions.
unsigned rg_INT rg_Curve::getID() const
{
	return ID;
}

rg_Planarity rg_Curve::getPlanarity() const 
{
	return planarOrSpatial;
}

////	Set Functions.
void rg_Curve::setID(const unsigned rg_INT &newID)
{
	ID = newID;
}

void rg_Curve::setPlanarity(const rg_Planarity &newPlanarity)
{
	planarOrSpatial = newPlanarity;
}

void rg_Curve::removeAll()
{
    ID=0;
    planarOrSpatial=rg_PLANAR;
}




