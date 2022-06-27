//********************************************************************
//
//	  FILENAME    : rg_BoundingBox3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_BoundingBox3D
//           which define a bounding box of curve and surface in 3D space. 
//                          
//	  CLASS NAME  : rg_BoundingBox3D
//
//    BASE CLASS  : None
//      
//    AUTHOR      : Deok-Soo Kim, Youngsong Cho
//
//    HISTORY     : 
//
//    START DATE  : 6 Mar. 2000    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_BOUNDINGBOX3D_H
#define _RG_BOUNDINGBOX3D_H

#include "rg_Const.h"
#include "rg_Point3D.h"
#include "rg_BzSurface3D.h"
#include "Sphere.h"

class rg_BoundingBox3D
{
private:
	rg_Point3D minPt;
			//	front  : max of x
			//	left   : min of y
			//	top    : max of z
	rg_Point3D maxPt;
			//	back   : min of x
			//	right  : max of y
			//	bottom : min of z

public:
	rg_BoundingBox3D();
	rg_BoundingBox3D(const rg_Point3D& tMinPt, const rg_Point3D& tMaxPt);
	rg_BoundingBox3D(const rg_BoundingBox3D& temp);
	rg_BoundingBox3D( rg_dList<Sphere>& spheres );
	~rg_BoundingBox3D();

	//  get & set functions.
	rg_Point3D getMinPt() const;
	rg_Point3D getMaxPt() const;
	rg_Point3D getCenterPt() const;
	rg_REAL    getXLength() const;   
	rg_REAL    getYLength() const;   
	rg_REAL    getZLength() const;   
	rg_REAL    getLongestLength() const;

	void    setMinPt(const rg_Point3D& tlf);
	void    setMaxPt(const rg_Point3D& brb);
	void    setAll(const rg_Point3D& tMinPt, const rg_Point3D& tMaxPt);
	void    setAll(const rg_BoundingBox3D& temp);

	void    reset();

    rg_BoundingBox3D evaluateRelativeOffset(const rg_REAL& ratio) const;
    rg_BoundingBox3D evaluateAbsoluteOffset(const rg_REAL& dist) const;

	//  operator overloading
	rg_BoundingBox3D& operator =(const rg_BoundingBox3D& temp);

	//  bounding box 3d of curve and surface.
//	void calculateBoundingBox3DOfBezierSurface(const rg_BzSurface3D& aBzSurface);
	void contain(const rg_INT& numPts, rg_Point3D*   pts);
	void contain(const rg_Point3D& pt);

	rg_FLAG  doContain(const rg_Point3D& pt) const;
	rg_FLAG  isNull() const;

	// functions for constructing the bounding box of a sphere set
	void updateBoxByAddingSphere(const Sphere& sphere);
	void constructBoxByAddingSpheres(Sphere* spheres, const int& numSpheres);
	void constructBoxByAddingSpheres(::rg_dList<Sphere>& spheres);

};

#endif


