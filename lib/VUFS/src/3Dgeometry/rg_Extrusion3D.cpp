//******************************************************************************
//
//	  FILENAME    : rg_Extrusion3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Extrusion3D 
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

#include "rg_Extrusion3D.h"

void  rg_Extrusion3D::inititializeLocalCoordinates()
{
    localOrigin.setPoint(0.0,0.0,0.0);
    localXAxis.setPoint(1.0,0.0,0.0);
    localZAxis.setPoint(0.0,0.0,1.0);
}

rg_Extrusion3D::rg_Extrusion3D()
:rg_PrimitiveExtrusion3D()
{
}

rg_Extrusion3D::rg_Extrusion3D(const GeometryType& tType,const unsigned rg_INT& tID)
:rg_PrimitiveExtrusion3D(tType,tID)
{
}

rg_Extrusion3D::rg_Extrusion3D(const rg_PrimitiveExtrusion3D& temp)
:rg_PrimitiveExtrusion3D(temp)
{
    inititializeLocalCoordinates();
}

rg_Extrusion3D::rg_Extrusion3D(const rg_Extrusion3D& temp)
:rg_PrimitiveExtrusion3D(temp)
{
    localOrigin=temp.localOrigin;
    localXAxis=temp.localXAxis;
    localZAxis=temp.localZAxis;
}

rg_Extrusion3D::~rg_Extrusion3D()
{
}

// get function
rg_Point3D rg_Extrusion3D::getLocalOrigin() const
{
    return localOrigin;
}

rg_Point3D rg_Extrusion3D::getLocalXAxis() const
{
    return localXAxis;
}

rg_Point3D rg_Extrusion3D::getLocalZAxis() const
{
    return localZAxis;
}

// set function
void  rg_Extrusion3D::setLocalOrigin(const rg_Point3D& tLocalOrigin)
{
    localOrigin=tLocalOrigin;
}

void  rg_Extrusion3D::setLocalXAxis(const rg_Point3D& tLocalXAxis)
{
    localXAxis=tLocalXAxis;
}

void  rg_Extrusion3D::setLocalZAxis(const rg_Point3D& tLocalZAxis)
{
    localZAxis=tLocalZAxis;
}

void  rg_Extrusion3D::setExtrusion(const rg_Extrusion3D& temp)
{
    rg_PrimitiveExtrusion3D::setPrimitiveExtrusion(temp);
    localOrigin=temp.localOrigin;
    localXAxis=temp.localXAxis;
    localZAxis=temp.localZAxis;
}

rg_Extrusion3D& rg_Extrusion3D::operator=(const rg_Extrusion3D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }

    rg_Extrusion3D::setExtrusion(temp);
    return *this;
}


