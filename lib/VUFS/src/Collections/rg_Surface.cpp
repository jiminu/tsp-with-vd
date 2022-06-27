//********************************************************************
//
//	  FILENAME    : rg_Surface.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Surface 
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun 1996    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#include "rg_Surface.h"

////	Constructor & Destructor.
rg_Surface::rg_Surface()
: ID(0), planarOrSpatial(rg_SPATIAL)
{
}

rg_Surface::rg_Surface( const unsigned rg_INT &newID, 
				  const rg_Planarity &newPlanarity )
: ID(newID), planarOrSpatial(newPlanarity)
{
}

//// Copy constructor : March 13 1997
rg_Surface::rg_Surface( const rg_Surface &surface)
: ID(surface.ID), planarOrSpatial(surface.planarOrSpatial)
{
}

rg_Surface::~rg_Surface()
{
}

////	Get Functions.
unsigned rg_INT rg_Surface::getID() const
{
	return ID;
}

rg_Planarity rg_Surface::getPlanarity() const
{
	return planarOrSpatial;
}

////	Set Functions.
void rg_Surface::setID( const unsigned rg_INT &newID )
{
	ID = newID;
}

void rg_Surface::setPlanarity( const rg_Planarity &newPlanarity )
{
	planarOrSpatial = newPlanarity;
}






