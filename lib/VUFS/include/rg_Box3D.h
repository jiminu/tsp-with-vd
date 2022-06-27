//******************************************************************************
//
//	  FILENAME    : rg_Box3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_Box3D
//
//                local z axis        
//                 |
//                 + height
//                 |
//                 |    
//    local origin +---------+--- local y axis
//                /         width 
//               / 
//              + breadth 
//             /
//           local  x axis
// 
//	  CLASS NAME  : rg_Box3D
//
//    BASE CLASS  : rg_PrimitiveBox3D
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 4. 3
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************


#ifndef _rg_BOX_H
#define _rg_BOX_H

#include "rg_PrimitiveBox3D.h"
#include "rg_Point3D.h"
#include "rg_Const.h"

class rg_Box3D: public rg_PrimitiveBox3D
{
private:
    rg_Point3D localXAxis;
    rg_Point3D localZAxis;

    rg_Point3D localOrigin;
public:
    rg_Box3D(const unsigned rg_INT& tID=0);
    rg_Box3D(const GeometryType& tType,const unsigned rg_INT& tID=0);
    rg_Box3D(const rg_Box3D& temp);
    ~rg_Box3D();

    // Get function
    rg_Point3D getLocalXAxis() const;
    rg_Point3D getLocalYAxis() const;
    rg_Point3D getLocalZAxis() const;
    rg_Point3D getLocalOrigin() const;
 
    // set function
    void setLocalXAxis(const rg_Point3D& tLocalXAxis);
    void setLocalZAxis(const rg_Point3D& tLocalZAxis);
    void setLocalOrigin(const rg_Point3D& tLocalOrigin);
    void setBox(const rg_Box3D& temp);
    rg_Box3D& operator=(const rg_Box3D& temp);
};

#endif

