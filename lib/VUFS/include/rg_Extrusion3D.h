//******************************************************************************
//
//	  FILENAME    : rg_Extrusion3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_Extrusion3D 
//                          
//	  CLASS NAME  : rg_Extrusion3D
//
//    BASE CLASS  : rg_PrimitiveExtrusion3D
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 8.22
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************

#ifndef _RG_EXTRUSION3D_H
#define _RG_EXTRUSION3D_H

#include "rg_PrimitiveExtrusion3D.h"
#include "rg_Const.h"

class rg_Extrusion3D: virtual public rg_PrimitiveExtrusion3D
{
protected:
    rg_Point3D localOrigin;
    rg_Point3D localZAxis;
    rg_Point3D localXAxis;
public:
    rg_Extrusion3D();
    rg_Extrusion3D(const GeometryType& tType,const unsigned rg_INT& tID=0);
    rg_Extrusion3D(const rg_PrimitiveExtrusion3D& temp);
    rg_Extrusion3D(const rg_Extrusion3D& temp);
    ~rg_Extrusion3D();

    void  inititializeLocalCoordinates();
    // get function
    rg_Point3D getLocalOrigin() const;
    rg_Point3D getLocalXAxis() const;
    rg_Point3D getLocalZAxis() const;

    // set function
    void  setLocalOrigin(const rg_Point3D& tLocalOrigin);
    void  setLocalXAxis(const rg_Point3D& tLocalXAxis);
    void  setLocalZAxis(const rg_Point3D& tLocalZAxis);
    void  setExtrusion(const rg_Extrusion3D& temp);

    rg_Extrusion3D& operator=(const rg_Extrusion3D& temp);
};

#endif


