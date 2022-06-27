//******************************************************************************
//
//	  FILENAME    : rg_Link3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class 'rg_Link3D' 
//           which will be the base class of 'rg_Link3D'.
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
#ifndef _RG_LINK3D_H
#define _RG_LINK3D_H
#include "rg_Point3D.h"
#include "rg_Geometry.h"

class rg_Link3D : virtual public rg_Geometry
{
protected:
	rg_Point3D startPt;
	rg_Point3D endPt;
	rg_Point3D localZAxis;

	rg_REAL  width;
	rg_REAL  height;

public:
    rg_Link3D(const unsigned rg_INT& tID=0);
    rg_Link3D(const GeometryType& tType,const unsigned rg_INT& tID=0);
    rg_Link3D(const rg_Link3D& temp);
    ~rg_Link3D();

    rg_REAL   getWidth() const;
    rg_REAL   getHeight() const;

	rg_Point3D  getStartPt() const;
	rg_Point3D  getEndPt() const;
    rg_Point3D  getLocalZAxis() const;
    rg_Point3D  getLocalXAxis() const;
    rg_Point3D  getLocalOrigin() const;

    void   setWidth( const rg_REAL& tWidth);
    void   setHeight( const rg_REAL& tHeight);
    void   setLink( const rg_Link3D& temp);

    void   setStartPt(const rg_Point3D& tStartPt);
    void   setEndPt(const rg_Point3D& tEndPt);
    void   setLocalZAxis(const rg_Point3D& tLocalZAxis);

    rg_Link3D& operator=( const rg_Link3D& temp );
};

#endif


