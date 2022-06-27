//*****************************************************************************
//
//	  FILENAME    : rg_PartialCylinder3D.cpp
//	  
//    DESCRIPTION : 
//           This is the interfacef the class rg_PartialCylinder3D 
//                          
//	  CLASS NAME  : rg_PartialCylinder3D
//
//    BASE CLASS  : rg_Cylinder3D
//
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 4.2
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//*****************************************************************************
#ifndef _RG_PARTIALCYLINDER3D_H
#define _RG_PARTIALCYLINDER3D_H

#include "rg_Const.h"
#include "rg_Cylinder3D.h"
#include "rg_RelativeOp.h"

class rg_PartialCylinder3D: virtual public rg_Cylinder3D
{
protected:
    rg_Point3D localXAxis;
    rg_REAL  angle;

public:
    rg_PartialCylinder3D(const unsigned rg_INT& tID);
    rg_PartialCylinder3D(const GeometryType& tType, 
                    const unsigned rg_INT& tID);
    rg_PartialCylinder3D(const rg_Cylinder3D& tCylinder);
    rg_PartialCylinder3D(const rg_PartialCylinder3D& temp);
    ~rg_PartialCylinder3D();

    rg_REAL convertAngle(const rg_REAL& tAngle) const;
    rg_Point3D getLocalXAxis() const;

    rg_REAL getAngle() const;

    void   setLocalXAxis (const rg_Point3D& tLocalXAxis);
    void   setAngle( const rg_REAL& tAngle);
    void   setPartialCylinder(const rg_Cylinder3D& tCylinder);
    void   setPartialCylinder(const rg_PartialCylinder3D& temp);
    rg_PartialCylinder3D& operator=(const rg_Cylinder3D& tCylinder);
    rg_PartialCylinder3D& operator=(const rg_PartialCylinder3D& temp);
};

#endif


