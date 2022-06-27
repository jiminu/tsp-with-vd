//******************************************************************************
//
//    FILENAME    : rg_Frustum3D.h
//    
//    DESCRIPTION : 
//           This is the interface of the class rg_Cylinder3D 
//                          
//    CLASS NAME  : rg_Frustum3D
//
//    BASE CLASS  : rg_PrimitiveFrustum3D, rg_Cylinder3D
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 6.1
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea         
//
//******************************************************************************

#ifndef _RG_FRUSTUM3D_H
#define _RG_FRUSTUM3D_H

#include "rg_PrimitiveFrustum3D.h"
#include "rg_Cylinder3D.h"
#include "rg_Point3D.h"

class rg_Frustum3D: public rg_PrimitiveFrustum3D, public rg_Cylinder3D
{
public:
    rg_Frustum3D(const unsigned rg_INT& tID=0);
    rg_Frustum3D(const GeometryType& tType, const unsigned rg_INT& tID=0);
    rg_Frustum3D(const rg_Frustum3D& temp);
    ~rg_Frustum3D();
    void   setFrustum(const rg_Frustum3D& temp);
    rg_Frustum3D& operator=(const rg_Frustum3D& temp);
};
#endif

