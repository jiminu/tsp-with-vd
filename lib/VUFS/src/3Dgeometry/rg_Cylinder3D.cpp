//******************************************************************************
//
//	  FILENAME    : rg_Cylinder3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Cylinder3D 
//                          
//	  CLASS NAME  : rg_Cylinder3D
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
#include "rg_Cylinder3D.h"


rg_Cylinder3D::rg_Cylinder3D(const unsigned rg_INT& tID)
:rg_Geometry(CYLINDER,tID),rg_PrimitiveCylinder3D(CYLINDER,tID), localOrigin(0.0,0.0,0.0), localZAxis(0.0,0.0,1.0)
{
}

rg_Cylinder3D::rg_Cylinder3D(const GeometryType& tType, const unsigned rg_INT& tID)
:rg_Geometry(tType,tID),rg_PrimitiveCylinder3D(tType,tID), localOrigin(0.0,0.0,0.0), localZAxis(0.0,0.0,1.0)
{
}

rg_Cylinder3D::rg_Cylinder3D(const rg_Cylinder3D& temp)
:rg_Geometry(temp),rg_PrimitiveCylinder3D(temp),localOrigin(temp.localOrigin), localZAxis(temp.localZAxis)
{
}

rg_Cylinder3D::~rg_Cylinder3D()
{
}

rg_Point3D  rg_Cylinder3D::getLocalOrigin() const
{
    return localOrigin;
}

rg_Point3D   rg_Cylinder3D::getLocalZAxis() const
{
    return localZAxis;
}

void   rg_Cylinder3D::setLocalOrigin(const rg_Point3D& tLocalOrigin)
{
    localOrigin=tLocalOrigin;
}

void   rg_Cylinder3D::setLocalZAxis(const rg_Point3D& tLocalZAxis)
{
    localZAxis=tLocalZAxis.getUnitVector();
}

void   rg_Cylinder3D::setCylinder(const rg_Cylinder3D& temp)
{
    rg_PrimitiveCylinder3D::setPrimitiveCylinder(temp);
    localOrigin=temp.localOrigin;
    localZAxis=temp.localZAxis;
}

rg_Cylinder3D& rg_Cylinder3D::operator=(const rg_Cylinder3D& temp)
{
    if ( this == &temp)
    {
        return *this;
    }

    rg_PrimitiveCylinder3D::setPrimitiveCylinder(temp);
    localOrigin=temp.localOrigin;
    localZAxis=temp.localZAxis;
    
    return *this;
}


