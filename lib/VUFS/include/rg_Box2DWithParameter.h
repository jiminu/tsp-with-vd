//******************************************************************************
//
//	  FILENAME    : Box2D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_PrimitiveBox3D
//           which will be the base of all kinds of rg_Box3D
//
//	  CLASS NAME  : Box2D
//
//    BASE CLASS  : NONE
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 4. 3
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************



#ifndef _RG_BOX2DWITHPARAMETER_H
#define _RG_BOX2DWITHPARAMETER_H

#include "rg_Point2D.h"

class rg_Box2DWithParameter
{
public:
	rg_Point2D SP;
	rg_Point2D EP;
	rg_REAL  ST;
	rg_REAL  ET;
public:
	rg_Box2DWithParameter() { }
	rg_Box2DWithParameter(const rg_Point2D &pt1, const rg_Point2D &pt2, 
	      const rg_REAL  &t1,  const rg_REAL  &t2):SP(pt1), EP(pt2), ST(t1), ET(t2){ }
	~rg_Box2DWithParameter() { }
};

#endif

