//******************************************************************************
//
//	  FILENAME    : rg_Box3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Box3D
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



#include "rg_Box3D.h"

rg_Box3D::rg_Box3D(const unsigned rg_INT& tID)
: rg_Geometry(BOX,tID),rg_PrimitiveBox3D(BOX,tID)
, localXAxis( 1.0, 0.0, 0.0 )
, localZAxis( 0.0, 0.0, 1.0 )
, localOrigin( 0.0, 0.0, 0.0 )
{
}

rg_Box3D::rg_Box3D(const GeometryType& tType,const unsigned rg_INT& tID)
: rg_Geometry(tType,tID),rg_PrimitiveBox3D(tType,tID)
, localXAxis( 1.0, 0.0, 0.0 )
, localZAxis( 0.0, 0.0, 1.0 )
, localOrigin( 0.0, 0.0, 0.0 )
{
}

rg_Box3D::rg_Box3D(const rg_Box3D& temp)
: rg_Geometry(temp),rg_PrimitiveBox3D(temp)
, localXAxis( temp.localXAxis )
, localZAxis( temp.localZAxis )
, localOrigin( temp.localOrigin )
{
}

rg_Box3D::~rg_Box3D()
{
}

    // Get function
rg_Point3D rg_Box3D::getLocalXAxis() const
{
    return localXAxis;
}

rg_Point3D rg_Box3D::getLocalYAxis() const
{
    return localZAxis*localXAxis;
}


rg_Point3D rg_Box3D::getLocalZAxis() const
{
    return localZAxis;
}

rg_Point3D rg_Box3D::getLocalOrigin() const
{
    return localOrigin;
}


// set function
void rg_Box3D::setLocalXAxis(const rg_Point3D& tLocalXAxis)
{
    localXAxis=tLocalXAxis;
}

void rg_Box3D::setLocalZAxis(const rg_Point3D& tLocalZAxis)
{
    localZAxis=tLocalZAxis;
}

void rg_Box3D::setLocalOrigin(const rg_Point3D& tLocalOrigin)
{
    localOrigin=tLocalOrigin;
}

void rg_Box3D::setBox(const rg_Box3D& temp)
{
    rg_PrimitiveBox3D::setPrimitiveBox(temp);
    localXAxis=temp.localXAxis;
    localZAxis=temp.localZAxis;
    localOrigin=temp.localOrigin;
}

rg_Box3D& rg_Box3D::operator=(const rg_Box3D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }
    
    // set this box
    rg_PrimitiveBox3D::setPrimitiveBox(temp);
    localXAxis=temp.localXAxis;
    localZAxis=temp.localZAxis;
    localOrigin=temp.localOrigin;

    return *this;
}


