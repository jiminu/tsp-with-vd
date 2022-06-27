//******************************************************************************
//
//	  FILENAME    : rg_PrimitiveBox3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_PrimitiveBox3D
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

#include "rg_PrimitiveBox3D.h"

rg_PrimitiveBox3D::rg_PrimitiveBox3D(const unsigned rg_INT& tID)
: rg_Geometry(PRIMITIVE_BOX,tID),breadth(0.0),width(0.0),height(0.0)
{
}

rg_PrimitiveBox3D::rg_PrimitiveBox3D(const GeometryType& tType,const unsigned rg_INT& tID)
: rg_Geometry(tType,tID),breadth(0.0),width(0.0),height(0.0)
{
}

rg_PrimitiveBox3D::rg_PrimitiveBox3D(const rg_PrimitiveBox3D& temp)
: rg_Geometry(temp)
, breadth( temp.breadth)
, width( temp.width)
, height( temp.height)
{
}

rg_PrimitiveBox3D::~rg_PrimitiveBox3D()
{
}

rg_REAL   rg_PrimitiveBox3D::getBreadth() const
{
    return breadth;
}

rg_REAL   rg_PrimitiveBox3D::getWidth() const
{
    return width;
}

rg_REAL   rg_PrimitiveBox3D::getHeight() const
{
    return height;
}

void   rg_PrimitiveBox3D::setBreadth(const rg_REAL& tBreadth) 
{
    breadth=tBreadth;
}

void   rg_PrimitiveBox3D::setWidth(const rg_REAL& tWidth) 
{
    width=tWidth;
}

void   rg_PrimitiveBox3D::setHeight(const rg_REAL& tHeight) 
{
    height=tHeight;
}


void   rg_PrimitiveBox3D::setPrimitiveBox( const rg_PrimitiveBox3D& temp)
{
    rg_Geometry::setGeometry(temp);
    breadth=temp.breadth;
    width=temp.width;
    height=temp.height;
}

rg_PrimitiveBox3D& rg_PrimitiveBox3D::operator=( const rg_PrimitiveBox3D& temp )
{
    if ( this == &temp )
    {
        return *this;
    }

    rg_Geometry::setGeometry(temp);
    breadth=temp.breadth;
    width=temp.width;
    height=temp.height;

    return *this;
}


