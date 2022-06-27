//******************************************************************************
//
//	  FILENAME    : rg_PartialCylinder3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Cylinder3D 
//                          
//	  CLASS NAME  : rg_PartialCylinder3D
//
//    BASE CLASS  : rg_Cylinder3D
//
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 4.2
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************
#include "rg_PartialCylinder3D.h"
#include "rg_RelativeOp.h"


rg_PartialCylinder3D::rg_PartialCylinder3D(const unsigned rg_INT& tID)
: rg_Geometry(PARTIAL_CYLINDER,tID),rg_Cylinder3D(PARTIAL_CYLINDER,tID), 
  localXAxis(1.0,0.0,0.0), angle(0.0)
{
}

rg_PartialCylinder3D::rg_PartialCylinder3D(const GeometryType& tType, const unsigned rg_INT& tID)
: rg_Geometry(tType,tID),rg_Cylinder3D(tType,tID), 
  localXAxis(1.0,0.0,0.0), angle(0.0)
{
}

rg_PartialCylinder3D::rg_PartialCylinder3D(const rg_Cylinder3D& tCylinder)
: rg_Geometry(tCylinder),rg_Cylinder3D(tCylinder),
  localXAxis(1.0,0.0,0.0),angle(360.0)
{
}

rg_PartialCylinder3D::rg_PartialCylinder3D(const rg_PartialCylinder3D& temp)
: rg_Geometry(temp),rg_Cylinder3D(temp), 
  localXAxis(temp.localXAxis), angle(temp.angle)
{
}

rg_PartialCylinder3D::~rg_PartialCylinder3D()
{
}

rg_REAL rg_PartialCylinder3D::convertAngle(const rg_REAL& tAngle) const
{
    rg_REAL output=fmod(tAngle,360);
    if ( rg_LT(output,0.0) )
    {
        output=output+360.0;
    }
    return output;
}


rg_Point3D rg_PartialCylinder3D::getLocalXAxis() const
{
    return localXAxis;
}

rg_REAL rg_PartialCylinder3D::getAngle() const
{
    return angle;
}

void   rg_PartialCylinder3D::setLocalXAxis (const rg_Point3D& tLocalXAxis)
{
    localXAxis=tLocalXAxis.getUnitVector();
}
    
void   rg_PartialCylinder3D::setAngle( const rg_REAL& tAngle)
{
    angle=convertAngle(tAngle);
}


void   rg_PartialCylinder3D::setPartialCylinder(const rg_Cylinder3D& tCylinder)
{
    rg_Cylinder3D::setCylinder(tCylinder);
    localXAxis=rg_Point3D(1.0, 0.0, 0.0); 
    angle=360.0;
}

void   rg_PartialCylinder3D::setPartialCylinder(const rg_PartialCylinder3D& temp)
{
    rg_Cylinder3D::setCylinder(temp);
    localXAxis=temp.localXAxis;
    angle=temp.angle;
}

rg_PartialCylinder3D& rg_PartialCylinder3D::operator=(const rg_Cylinder3D& tCylinder)
{
    if ( this == &tCylinder)
    {
        return *this;
    }

    rg_Cylinder3D::setCylinder(tCylinder);
    localXAxis=rg_Point3D( 1.0, 0.0, 0.0);
    angle=360.0;
    return *this;
}

rg_PartialCylinder3D& rg_PartialCylinder3D::operator=(const rg_PartialCylinder3D& temp)
{
    if ( this == &temp)
    {
        return *this;
    }

    rg_Cylinder3D::setCylinder(temp);
    localXAxis=temp.localXAxis;
    angle=temp.angle;
    return *this;
}


