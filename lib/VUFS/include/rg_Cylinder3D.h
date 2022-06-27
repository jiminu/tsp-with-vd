//******************************************************************************
//
//	  FILENAME    : rg_Cylinder3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_Cylinder3D 
//                          
//	  CLASS NAME  : rg_PrimitiveCylinder3D
//
//    BASE CLASS  : NONE
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 4.2
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************

#ifndef _RG_CYLINDER3D_H
#define _RG_CYLINDER3D_H

#include "rg_PrimitiveCylinder3D.h"
#include "rg_Point3D.h"

class rg_Cylinder3D: virtual public rg_PrimitiveCylinder3D
{
protected:
    rg_Point3D localOrigin;
    rg_Point3D localZAxis;
public:
    rg_Cylinder3D(const unsigned rg_INT& tID=0);
    rg_Cylinder3D(const GeometryType& tType, const unsigned rg_INT& tID=0);
    rg_Cylinder3D(const rg_Cylinder3D& temp);
    ~rg_Cylinder3D();

    rg_Point3D  getLocalOrigin() const;
    rg_Point3D  getLocalZAxis() const;

    void   setLocalOrigin(const rg_Point3D& tLocalOrigin);
    void   setLocalZAxis(const rg_Point3D& tLocalZAxis);
    void   setCylinder(const rg_Cylinder3D& temp);
    rg_Cylinder3D& operator=(const rg_Cylinder3D& temp);
};
#endif


