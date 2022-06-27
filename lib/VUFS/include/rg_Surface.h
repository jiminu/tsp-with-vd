//******************************************************************************
//
//	  FILENAME    : rg_Surface.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_Surface 
//           which will be the base class of all kinds of surfaces. 
//                          
//	  CLASS NAME  : rg_Surface
//
//    BASE CLASS  : NONE
//      
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun 1996    
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************

#ifndef _RG_SURFACE_H
#define _RG_SURFACE_H

#include "rg_Const.h"
//#include "Defconst.h"

class rg_Surface
{
protected:
    unsigned rg_INT ID;
    rg_Planarity    planarOrSpatial;

public:
////    Constructor & Destructor.
    rg_Surface();
    rg_Surface( const unsigned rg_INT &newID, 
             const rg_Planarity    &newPlanarity );

    //// Copy constructor : March 13 1997
    rg_Surface( const rg_Surface &surface);
    virtual ~rg_Surface();

////    Get Functions.
    unsigned rg_INT getID() const;
    rg_Planarity    getPlanarity() const;

////    Set Functions.
    void         setID( const unsigned rg_INT &newID );
    void         setPlanarity( const rg_Planarity &newPlanarity );
};

#endif


