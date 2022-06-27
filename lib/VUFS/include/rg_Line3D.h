//********************************************************************
//
//	  FILENAME    : rg_Line3D.h
//	  
//    DESCRIPTION : 
//           This consists of the interface of class rg_Line3D.      
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 11 Jul. 1997    
//
//    HISTROY     :
//     1. By Young-Song Cho.   19 Jul. 1997
//         make : 
//            rg_FLAG isParallelToLine(const rg_Line3D &line) const;
//            rg_FLAG isPerpendicularToLine(const rg_Line3D &line) const;
//            rg_FLAG isThisPointOnLine(const rg_Point3D &pt, rg_REAL &u) const;
//
//     2. By Tae-bum Jang   1998. 7.13
//         make : 
//            rg_Point3D getUnitDirection() const;
//            rg_REAL getDistance(const rg_Point3D& pt) const;
//         update :
//            rg_REAL getLength(); -->  rg_REAL getLength() const;
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_LINE3D_H
#define _RG_LINE3D_H

#include "rg_Const.h"
#include "rg_Point3D.h"

enum INTERSECT_TYPE  { NON_INTERSECT, POINT_INTERSECT,
                       LINE_INTERSECT, PLANE_INTERSECT };
class rg_Line3D
{
protected:
	rg_Point3D sp;
	rg_Point3D ep;
public: 

	// Construction/Destruction
	rg_Line3D();
	rg_Line3D(const rg_Point3D &sPt, const rg_Point3D &ePt);
	rg_Line3D(const rg_Line3D &line);
	~rg_Line3D();            
              
	// Access elements	                       
	rg_Point3D  getSP() const;
	rg_Point3D  getEP() const;  

//	In 2D the following function is meaningful but in 3D is not.
//	rg_REAL  getTangent() const;
	
	void   setSP(const rg_Point3D &pt);
	void   setEP(const rg_Point3D &pt);
	void   setLine(const rg_Point3D &sPt, const rg_Point3D &ePt);
	void   setLine(const rg_Line3D &line);

	// Operations
	rg_REAL getLength() const;
    rg_REAL getDistance(const rg_Point3D& pt) const;
    rg_Point3D getUnitDirection() const;

    rg_FLAG isParallelToLine(const rg_Line3D &line) const;
    rg_FLAG isColinearWith(const rg_Line3D &line) const;

    rg_FLAG isPerpendicularToLine(const rg_Line3D &line) const;
    rg_FLAG isThisPointOnLine(const rg_Point3D &pt, rg_REAL &u) const;

    // Operator Overloading
    rg_FLAG operator==(const rg_Line3D &line) const;

    // Intersection
    rg_Point3D* intersectLine(const rg_Line3D &line, INTERSECT_TYPE &type) const;

};

//  uA + vB = C
//  find u, v 
//  The first returned values is u, the second v.
rg_REAL* solveSimultaneousEq(const rg_Point3D &A, const rg_Point3D &B, const rg_Point3D &C);

#endif


