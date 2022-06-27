//******************************************************************************
//
//	  FILENAME    : rg_PrimitiveCylinder3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_PrimitiveCylinder3D 
//           which will be the base of all kinds of rg_Cylinder3D
//                          
//	  CLASS NAME  : rg_PrimitiveCylinder3D
//
//    BASE CLASS  : rg_Geometry
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 4.2
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************

#ifndef _RG_PRIMITIVECYLINDER3D_H
#define _RG_PRIMITIVECYLINDER3D_H

#include "rg_Geometry.h"
#include "rg_Const.h"

class rg_PrimitiveCylinder3D: virtual public rg_Geometry
{
protected:
    rg_REAL  radius;
    rg_REAL  height;
public:
    rg_PrimitiveCylinder3D();
    rg_PrimitiveCylinder3D(const GeometryType& tType,const unsigned rg_INT& tID=0);
    rg_PrimitiveCylinder3D(const rg_PrimitiveCylinder3D& temp);
    ~rg_PrimitiveCylinder3D();
    rg_REAL   getRadius() const;
    rg_REAL   getHeight() const;

    void   setRadius(const rg_REAL& tRadius);
    void   setHeight(const rg_REAL& tHeight);
    void   setPrimitiveCylinder(const rg_PrimitiveCylinder3D& temp);
    rg_PrimitiveCylinder3D& operator=(const rg_PrimitiveCylinder3D& temp);
};

#endif


