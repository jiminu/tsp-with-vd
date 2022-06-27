//******************************************************************************
//
//	  FILENAME    : rg_ToroidalSurface3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class Toroidal surface
//                          
//	  CLASS NAME  : rg_ToroidalSurface3D
//
//    BASE CLASS  : NONE
//      
//
//    AUTHOR      : Ryu, Joonghyun
//    START DATE  : 1999. 10.23
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************

#ifndef _RG_TOROIDALSURFACE3D_H
#define _RG_TOROIDALSURFACE3D_H

#include "rg_Point3D.h"

class rg_ToroidalSurface3D
{

// Assume that parameterization is like followings
// T(u, v) = C + (R + rcos(v))(cos(u)X + sin(u)Y) + rsin(v)Z

private:
	rg_REAL minorRadius;
	rg_REAL majorRadius;
	rg_Point3D localOrigin;
	rg_Point3D localXAxis;
	rg_Point3D localYAxis;
	rg_Point3D localZAxis;

public:
	rg_ToroidalSurface3D();
	rg_ToroidalSurface3D(const rg_ToroidalSurface3D & sourceObj);
	~rg_ToroidalSurface3D();

	rg_REAL getMinorRadius() const;
	rg_REAL getMajorRadius() const;
	rg_Point3D getLocalOrigin() const;
	rg_Point3D getLocalXAxis() const;
	rg_Point3D getLocalYAxis() const;
	rg_Point3D getLocalZAxis() const;

	void setMinorRadius(const rg_REAL& minorR);
	void setMajorRadius(const rg_REAL& majorR);
	void setLocalOrigin(const rg_Point3D& localOrg);
	void setLocalXAxis(const rg_Point3D& localX);
	void setLocalYAxis(const rg_Point3D& localY);
	void setLocalZAxis(const rg_Point3D& localZ);

	rg_Point3D** triangulateSurface(rg_INT& numOfRow, rg_INT& numOfCol, rg_REAL majorTol = 30, rg_REAL minorTol = 30);
	
	rg_ToroidalSurface3D& operator=(const rg_ToroidalSurface3D& sourceObj);
};

#endif


