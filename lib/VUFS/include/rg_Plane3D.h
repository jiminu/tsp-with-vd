//********************************************************************
//
//	  FILENAME    : rg_Plane3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of class rg_Point3D.      
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 8 Jul. 1997    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_PLANE3D_H
#define _RG_PLANE3D_H

#include "rg_Const.h"
#include "rg_Point3D.h"
#include "rg_Line3D.h"
#include "rg_ListByPtr.h"


//                    ------------ 
//              vVec /          /
//                  /          /
//           z |   ------------  
//             |  /       uVec
//             | / posVec 
//             |/
//     origin  ------------- 
//            /           y
//           /
//         x/
//

class rg_Plane3D
{
protected:
    rg_Point3D  posVec;
    rg_Point3D  uVec;
    rg_Point3D  vVec;

public:
    //  Constructor & Destructor.
    rg_Plane3D();
    rg_Plane3D(const rg_Point3D& pVector, const rg_Point3D &uVector, const rg_Point3D &vVector);
    rg_Plane3D(const rg_Plane3D &plane);
    virtual ~rg_Plane3D();

    //  Access Functions.
    rg_Point3D   getPosVector() const;
    rg_Point3D   getUVector() const;
    rg_Point3D   getVVector() const;

    void setPosVector(const rg_Point3D &pVector);
    void setUVector(const rg_Point3D &uVector);
    void setVVector(const rg_Point3D &vVector);
    void setPlane(const rg_Point3D &pVector, const rg_Point3D &uVector, const rg_Point3D &vVector);
    //  cornerPt consists of 4 corner points of rg_Point3D.
    void setPlane(const rg_Point3D* const cornerPt);
    void setPlane(const rg_Plane3D &plane);

    //  Evaluate rg_Point3D.
    virtual rg_Point3D   evaluatePt(const rg_PARAMETER &u, const rg_PARAMETER &V) const;
    
    //  Operations.
    rg_FLAG    isParallelTo(const rg_Line3D &line) const;
    rg_FLAG    isParallelTo(const rg_Plane3D &plane) const;
    rg_FLAG    isPerpendicularTo(const rg_Line3D &line) const;
    rg_FLAG    isPerpendicularTo(const rg_Plane3D &plane) const;
    
    rg_FLAG    isThisPointOnPlane(const rg_Point3D &pt) const;

    //  Intersect
    rg_sListByPtr*  intersectLine(const rg_Line3D &line, rg_INT &numOfintersectPt) const;
    rg_sListByPtr*  intersectPlane(const rg_Plane3D &plane, rg_INT &numOfintersectPt) const;

    //  Operator Overloading.
    rg_Plane3D& operator =(const rg_Plane3D &plane);
    rg_FLAG  operator==(const rg_Plane3D &plane) const;

};

#endif


