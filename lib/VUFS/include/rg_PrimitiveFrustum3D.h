//******************************************************************************
//
//	  FILENAME    : rg_PrimitiveFrustum3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_PrimitiveFrustum3D 
//           which will be the base of all kinds of rg_Frustum3D
//                          
//	  CLASS NAME  : rg_PrimitiveFrustum3D
//
//    BASE CLASS  : rg_PrimitiveCylinder3D
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 4.2
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************

#ifndef _RG_PRIMITIVEFRUSTUM3D_H
#define _RG_PRIMITIVEFRUSTUM3D_H

#include "rg_PrimitiveCylinder3D.h"
#include "rg_Const.h"

class rg_PrimitiveFrustum3D: virtual public rg_PrimitiveCylinder3D
{
protected:
    rg_REAL  topRadius;

public:
    rg_PrimitiveFrustum3D(const unsigned rg_INT& tID=0);
    rg_PrimitiveFrustum3D(const GeometryType& tType,const unsigned rg_INT& tID=0);
    rg_PrimitiveFrustum3D(const rg_PrimitiveFrustum3D& temp);
    ~rg_PrimitiveFrustum3D();
    rg_REAL   getBottomRadius() const;
    rg_REAL   getTopRadius() const;

    void   setBottomRadius(const rg_REAL& tBottomRadius);
    void   setTopRadius(const rg_REAL& tTopRadius);

    void   setPrimitiveFrustum(const rg_PrimitiveFrustum3D& temp);
    rg_PrimitiveFrustum3D& operator=(const rg_PrimitiveFrustum3D& temp);
};

#endif


