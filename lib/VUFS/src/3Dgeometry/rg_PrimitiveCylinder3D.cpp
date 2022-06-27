//******************************************************************************
//
//	  FILENAME    : rg_PrimitiveCylinder3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_PrimitiveCylinder3D 
//           which will be the base class of all kinds of Cylinders. 
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

#include "rg_PrimitiveCylinder3D.h"


rg_PrimitiveCylinder3D::rg_PrimitiveCylinder3D()
: rg_Geometry(PRIMITIVE_CYLINDER),radius(0.0), height(0.0)
{
}

rg_PrimitiveCylinder3D::rg_PrimitiveCylinder3D(const GeometryType& tType,const unsigned rg_INT& tID)
: rg_Geometry(tType,tID),radius(0.0), height(0.0)
{
}

rg_PrimitiveCylinder3D::rg_PrimitiveCylinder3D(const rg_PrimitiveCylinder3D& temp)
: rg_Geometry(temp), radius(temp.radius), height(temp.height)
{
}

rg_PrimitiveCylinder3D::~rg_PrimitiveCylinder3D()
{
}

rg_REAL   rg_PrimitiveCylinder3D::getRadius() const
{
    return radius;
}

rg_REAL   rg_PrimitiveCylinder3D::getHeight() const
{
    return height;
}

void   rg_PrimitiveCylinder3D::setRadius(const rg_REAL& tRadius) 
{
    radius=tRadius;
}

void   rg_PrimitiveCylinder3D::setHeight(const rg_REAL& tHeight) 
{
    height=tHeight;
}

void   rg_PrimitiveCylinder3D::setPrimitiveCylinder(const rg_PrimitiveCylinder3D& temp)
{
    rg_Geometry::setGeometry(temp);
    radius=temp.radius;
    height=temp.height;
}

rg_PrimitiveCylinder3D& rg_PrimitiveCylinder3D::operator=(const rg_PrimitiveCylinder3D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }

    rg_Geometry::setGeometry(temp);
    radius=temp.radius;
    height=temp.height;

    return *this;
}


