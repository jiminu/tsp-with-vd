//********************************************************************
//
//	  FILENAME    : rg_Point3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_Point3D which defines 
//           the property of the point and vector in a 3-dimension.
//                          
//	  CLASS NAME  : rg_Point3D
//
//    BASE CLASS  : None
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : Mar 11 1997    
//
//    HISTORY     :
//          BY Young-Song Cho.  19 Jul. 1997
//             make   :  rg_Point3D getUnitVector() const;
//             update :  rg_REAL distance(const rg_Point3D &pt) const;
//                       ->  rg_REAL distance(const rg_Point3D &pt=rg_Point3D(0., 0., 0.)) const;
//	        BY Jung-Hyun Ryu    3 Jul. 1998
//			   make   :  rg_Point3D(const rg_Point2D &pt);	
//			   make   :  rg_Point3D& operator =(const Piont2D &pt);
//
//          By Youngsong Cho    15 Jul. 2000
//             make   :  friend rg_Point3D getUnitVector(const rg_Point3D &pt);
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_POINT3D_H
#define _RG_POINT3D_H


#include "rg_Const.h"
#include "rg_Matrix.h"

class rg_Point2D;

class rg_Point3D
{
private:
    rg_REAL x;
    rg_REAL y;
    rg_REAL z;

////  member function.  ////

public:
    ////  Constructor & Destructor  ////
    rg_Point3D();
    rg_Point3D(const rg_Point3D &pt);
    rg_Point3D(const rg_Point2D &pt);
    rg_Point3D(const rg_REAL &px,
               const rg_REAL &py,
               const rg_REAL &pz);            

    ~rg_Point3D();

    ////  Access Function  ////
    inline rg_REAL    getX() const {return x;};
    inline rg_REAL    getY() const {return y;};
    inline rg_REAL    getZ() const {return z;};
	inline void       getCoordinates(rg_REAL& px, rg_REAL& py, rg_REAL& pz) { px = x; py = y; pz = z; }
    rg_Point3D getUnitVector() const;

    void setX(const rg_REAL &px);
    void setY(const rg_REAL &py);
    void setZ(const rg_REAL &pz);
    void setPoint(const rg_Point3D &pt);
    void setPoint(const rg_REAL &px,
                  const rg_REAL &py,
                  const rg_REAL &pz);                
    
    ////  Operation & Calculation  ////
    rg_REAL  distance(const rg_Point3D &pt = rg_Point3D(0., 0., 0.)) const;
    rg_REAL  squaredDistance(const rg_Point3D &pt = rg_Point3D(0., 0., 0.)) const;

	rg_REAL    squaredMagnitude() const;
    rg_REAL    magnitude() const;                     
    rg_REAL    innerProduct(const rg_Point3D &pt) const;   
    rg_Point3D crossProduct(const rg_Point3D &pt) const;
	rg_Point2D evaluatePt2D() const;

	void       normalize();
	rg_REAL    angle(const rg_Point3D& vecter );

    rg_BOOL    isEqual(const rg_Point3D &pt, const rg_REAL& tolerance = resNeg6) const;

    ////  Operator Overloading  ////
    rg_Point3D& operator =(const rg_Point3D &pt);
	rg_Point3D& operator =(const rg_Point2D &pt);

    rg_Point3D& operator =(const rg_REAL &n);
        //  n must be zero.
        //  rg_BandedMatrix which is template class uses this operator.

    rg_Point3D  operator +(const rg_Point3D &pt) const;
    rg_Point3D& operator+=(const rg_Point3D &pt);

    rg_Point3D  operator -(const rg_Point3D &pt) const;
    rg_Point3D& operator-=(const rg_Point3D &pt);
    
    rg_Point3D  operator *(const rg_REAL &n) const;
    rg_Point3D  operator *(const rg_Matrix &mat) const;   

    rg_Point3D  operator /(const rg_REAL &n) const;

    rg_FLAG   operator==(const rg_Point3D &pt) const;
    rg_FLAG   operator!=(const rg_Point3D &pt) const;

    friend rg_Point3D operator -(const rg_Point3D &pt);
    friend rg_Point3D operator *(const rg_REAL &n, const rg_Point3D &pt);
    friend rg_Point3D operator *(const rg_Point3D &pt1, 
                            const rg_Point3D &pt2) ; // Cross Product
    friend rg_REAL  operator %(const rg_Point3D &pt1,
                            const rg_Point3D &pt2) ; // Dot Product

	friend rg_Point3D getUnitVector(const rg_Point3D &pt);


    static bool areThreePointsColinear(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_REAL& res = 1.0e-6);
};

#endif







