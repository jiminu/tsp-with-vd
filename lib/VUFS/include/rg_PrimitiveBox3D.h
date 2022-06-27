//******************************************************************************
//
//	  FILENAME    : rg_PrimitiveBox3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_PrimitiveBox3D
//           which will be the base of all kinds of rg_Box3D
//                 z         
//                 |
//                 + height
//                 |
//                 |    
//                 +---------+--- y
//                /         width 
//               / 
//              + breadth 
//             /
//             x
//	  CLASS NAME  : rg_PrimitiveBox3D
//
//    BASE CLASS  : rg_Geometry
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 4. 3
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************


#ifndef _RG_PRIMITIVEBOX3D_H
#define _RG_PRIMITIVEBOX3D_H

#include "rg_Geometry.h"
#include "rg_Const.h"

class rg_PrimitiveBox3D: virtual public rg_Geometry
{
protected:
    rg_REAL breadth;
    rg_REAL width;
    rg_REAL height;

public:
    rg_PrimitiveBox3D(const unsigned rg_INT& tID=0);
    rg_PrimitiveBox3D(const GeometryType& tType,const unsigned rg_INT& tID=0);
    rg_PrimitiveBox3D(const rg_PrimitiveBox3D& temp);
    ~rg_PrimitiveBox3D();
    rg_REAL   getBreadth() const;
    rg_REAL   getWidth() const;
    rg_REAL   getHeight() const;

    void   setBreadth( const rg_REAL& tBreadth);
    void   setWidth( const rg_REAL& tWidth);
    void   setHeight( const rg_REAL& tHeight);
    void   setPrimitiveBox( const rg_PrimitiveBox3D& temp);
    rg_PrimitiveBox3D& operator=( const rg_PrimitiveBox3D& temp );
};

#endif


