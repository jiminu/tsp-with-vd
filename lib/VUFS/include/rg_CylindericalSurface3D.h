//******************************************************************************
//
//	  FILENAME    : rg_CylindericalSurface3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class Toroidal surface
//                          
//	  CLASS NAME  : rg_CylindericalSurface3D
//
//    BASE CLASS  : NONE
//      
//
//    AUTHOR      : Ryu, Joonghyun
//    START DATE  : 1999. 10.28
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************

#ifndef _RG_CYLINDERICALSURFACE3D_H
#define _RG_CYLINDERICALSURFACE3D_H

#include "rg_Point3D.h"

class rg_CylindericalSurface3D
{
// Assume that parameterization is like followings
// C(u, v) = C + R(cos(u)X + (sin(u)Y) + vZ

protected:
	rg_Point3D localOrigin;
	rg_Point3D localXAxis;
	rg_Point3D localYAxis;
	rg_Point3D localZAxis;
	rg_REAL    radius;
	rg_REAL    height;

public:
	rg_CylindericalSurface3D();
	~rg_CylindericalSurface3D();

	rg_Point3D getlocalOrigin() const;
	rg_Point3D getlocalXAxis() const;
	rg_Point3D getlocalYAxis() const;
	rg_Point3D getlocalZAxis() const;
	rg_REAL    getRadius() const;
	rg_REAL    getHeight() const;

	void setLocalOrigin(const rg_Point3D& localOrg);
	void setLocalXAxis(const rg_Point3D& localX);
	void setLocalYAxis(const rg_Point3D& localY);
	void setLocalZAxis(const rg_Point3D& localZ);
	void    setRadius(const rg_REAL& r);
	void    setHeight(const rg_REAL& h);
	rg_Point3D** triangulateSurface(rg_INT& numOfRow, rg_INT& numOfCol, rg_REAL angleTol = 20);

};

#endif


