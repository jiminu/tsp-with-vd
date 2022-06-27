//******************************************************************************
//
//	  FILENAME    : rg_Link3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implements of the class 'rg_Link3D' 
//                          
//	  CLASS NAME  : rg_Link3D
//
//    BASE CLASS  : 
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 5.20
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************

#include "rg_Link3D.h"
rg_Link3D::rg_Link3D(const unsigned rg_INT& tID)
: rg_Geometry(LINK,tID), startPt(0.0,0.0,0.0), endPt(0.0,0.0,0.0), localZAxis(0.0,0.0,0.0)
{
    width=0.0;
    height=0.0;
}

rg_Link3D::rg_Link3D(const GeometryType& tType,const unsigned rg_INT& tID)
: rg_Geometry(tType,tID), startPt(0.0,0.0,0.0), endPt(0.0,0.0,0.0), localZAxis(0.0,0.0,0.0)
{
    width=0.0;
    height=0.0;
}

rg_Link3D::rg_Link3D(const rg_Link3D& temp)
: rg_Geometry( temp ), startPt(temp.startPt), endPt(temp.endPt), localZAxis(temp.localZAxis)
{
    width=temp.width;
    height=temp.height;
}

rg_Link3D::~rg_Link3D()
{
}

rg_REAL   rg_Link3D::getWidth() const
{
    return width;
}

rg_REAL   rg_Link3D::getHeight() const
{
    return height;
}

rg_Point3D  rg_Link3D::getStartPt() const
{
    return startPt;
}

rg_Point3D  rg_Link3D::getEndPt() const
{
    return endPt;
}

rg_Point3D  rg_Link3D::getLocalZAxis() const
{
    return localZAxis;
}

rg_Point3D  rg_Link3D::getLocalXAxis() const
{
    return (endPt-startPt).getUnitVector();
}

rg_Point3D  rg_Link3D::getLocalOrigin() const
{
    return startPt;
}

void   rg_Link3D::setWidth( const rg_REAL& tWidth)
{
    width=tWidth;
}

void   rg_Link3D::setHeight( const rg_REAL& tHeight)
{
    height=tHeight;
}

void   rg_Link3D::setLink( const rg_Link3D& temp)
{
    rg_Geometry::setGeometry(temp);
    height=temp.height;
    width=temp.width;
    startPt=temp.startPt;
    endPt=temp.endPt;
    localZAxis=temp.localZAxis;
}

void   rg_Link3D::setStartPt(const rg_Point3D& tStartPt)
{
    startPt=tStartPt;
}

void   rg_Link3D::setEndPt(const rg_Point3D& tEndPt)
{
    endPt=tEndPt;
}

void   rg_Link3D::setLocalZAxis(const rg_Point3D& tLocalZAxis)
{
    localZAxis=tLocalZAxis;
}


