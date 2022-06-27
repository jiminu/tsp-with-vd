//********************************************************************
//
//	  FILENAME    : rg_Point3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Point3D 
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
//             make   :  rg_Point3D rg_Point3D::getUnitVector() const
//
//	        BY Jung-Hyun Ryu    3 Jul. 1998
//			   make   :  rg_Point3D(const rg_Point2D &pt);
//			   make   :  rg_Point3D& operator =(const Piont2D &pt);
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************


#include "rg_Point3D.h"

#include <math.h>
#include "rg_RelativeOp.h"

#include "rg_Point2D.h"
#include "Plane.h"


////  Constructor & Destructor  ////
rg_Point3D::rg_Point3D()
: x(0.0), y(0.0), z(0.0)       
{
}

rg_Point3D::rg_Point3D(const rg_Point3D &pt)
: x(pt.x), y(pt.y), z(pt.z)
{
}

rg_Point3D::rg_Point3D(const rg_Point2D &pt)
: x(pt.getX()), y(pt.getY()), z(0.0)
{
}

rg_Point3D::rg_Point3D(const rg_REAL &px,
                 const rg_REAL &py,
                 const rg_REAL &pz)            
: x(px), y(py), z(pz)
{
}

rg_Point3D::~rg_Point3D()
{
}
/*
////  Access Function  ////
inline rg_REAL rg_Point3D::getX() const
{
    return x;
}

inline rg_REAL rg_Point3D::getY() const
{
    return y;
}

inline rg_REAL rg_Point3D::getZ() const
{
    return z;
}
*/

void rg_Point3D::setX(const rg_REAL &px)
{
    x = px;
}

rg_Point3D rg_Point3D::getUnitVector() const
{
    return *this/magnitude();
}

void rg_Point3D::normalize()
{
	rg_REAL mag = magnitude();

	x = x/mag;
	y = y/mag;
	z = z/mag;
}

void rg_Point3D::setY(const rg_REAL &py)
{
    y = py;
}

void rg_Point3D::setZ(const rg_REAL &pz)
{
    z = pz;
}

void rg_Point3D::setPoint(const rg_Point3D &pt)
{
    x = pt.x;
    y = pt.y;
    z = pt.z;
}

void rg_Point3D::setPoint(const rg_REAL &px,
                       const rg_REAL &py,
                       const rg_REAL &pz)                
{
    x = px;
    y = py;
    z = pz;
}

////  Operation & Calculation  ////
rg_REAL rg_Point3D::distance(const rg_Point3D &pt) const
{
    return sqrt( (x-pt.x)*(x-pt.x) +
                 (y-pt.y)*(y-pt.y) +
                 (z-pt.z)*(z-pt.z) );                 
}


rg_REAL  rg_Point3D::squaredDistance(const rg_Point3D &pt) const
{
    return ((x - pt.x)*(x - pt.x) + (y - pt.y)*(y - pt.y) + (z - pt.z)*(z - pt.z));
}


rg_REAL rg_Point3D::squaredMagnitude() const
{
    return ( x*x + y*y + z*z );
}

rg_REAL rg_Point3D::magnitude() const                     
{
    return sqrt( x*x + y*y + z*z );
}

rg_REAL rg_Point3D::innerProduct(const rg_Point3D &pt) const
{
    return x*pt.x + y*pt.y + z*pt.z;
}

rg_Point3D rg_Point3D::crossProduct(const rg_Point3D &pt) const   
//   |   i    j    k  |
//   |   x    y    z  |
//   | pt.x pt.y pt.z |
{
    return rg_Point3D( (y*pt.z - z*pt.y),
                    (z*pt.x - x*pt.z),
                    (x*pt.y - y*pt.x) );
}
rg_Point2D rg_Point3D::evaluatePt2D() const
{
	return rg_Point2D(x,y);
}

// return radian.
rg_REAL rg_Point3D::angle(const rg_Point3D& vecter )
{

	rg_REAL cosineTheta = 0.0;
	rg_REAL theta       = 0.0;
	//rg_REAL	RADIAN = 57.2957795130823;
    rg_REAL  pi=4.0*atan(1.0);
	rg_Point3D vector1 = getUnitVector();
	rg_Point3D vector2 = vecter.getUnitVector();

	cosineTheta = vector1.innerProduct(vector2);
	if ( cosineTheta > 1.0 )
	{
		theta = 0.0;
	}
	else if ( cosineTheta < -1.0 )
	{
		theta = pi;
	}
	else
	{
		theta = acos(cosineTheta);
	}

	return theta;
}



rg_BOOL rg_Point3D::isEqual(const rg_Point3D &pt, const rg_REAL& tolerance) const
{
    if ( rg_EQ(x, pt.x, tolerance) && rg_EQ(y, pt.y, tolerance) && rg_EQ(z, pt.z, tolerance) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}


////  Operator Overloading  ////
rg_Point3D& rg_Point3D::operator =(const rg_Point3D &pt)
{
    x = pt.x;
    y = pt.y;
    z = pt.z;

    return *this;
}

rg_Point3D& rg_Point3D::operator =(const rg_Point2D &pt)
{
    x = pt.getX();
    y = pt.getY();
    z = 0.0;

    return *this;
}

rg_Point3D& rg_Point3D::operator =(const rg_REAL &n)  
    //  n must be zero.
    //  bandedMatrix which is template class uses this operator.
{
    x = n;
    y = n;
    z = n;
	
    return *this;
}

rg_Point3D rg_Point3D::operator +(const rg_Point3D &pt) const
{
    return rg_Point3D(x+pt.x, y+pt.y, z+pt.z);
}

rg_Point3D& rg_Point3D::operator+=(const rg_Point3D &pt)
{
    x += pt.x;
    y += pt.y;
    z += pt.z;

    return *this;
}


rg_Point3D rg_Point3D::operator -(const rg_Point3D &pt) const
{
    return rg_Point3D(x-pt.x, y-pt.y, z-pt.z);
}

rg_Point3D& rg_Point3D::operator-=(const rg_Point3D &pt)
{
    x -= pt.x;
    y -= pt.y;
    z -= pt.z;

    return *this;
}

rg_Point3D rg_Point3D::operator *(const rg_REAL &n) const
{
    return rg_Point3D(n*x, n*y, n*z);
}

rg_Point3D rg_Point3D::operator *(const rg_Matrix &mat) const
{
    rg_REAL px = x*mat.getElement(0,0) + y*mat.getElement(1,0)
              + z*mat.getElement(2,0);        
    rg_REAL py = x*mat.getElement(0,1) + y*mat.getElement(1,1)
              + z*mat.getElement(2,1);        
    rg_REAL pz = x*mat.getElement(0,2) + y*mat.getElement(1,2)
              + z*mat.getElement(2,2);        

    return rg_Point3D(px, py, pz);
}

rg_Point3D rg_Point3D::operator /(const rg_REAL &n) const
{
    return rg_Point3D(x/n, y/n, z/n);
}

rg_FLAG rg_Point3D::operator==(const rg_Point3D &pt) const
{
    if ( rg_EQ(x, pt.x) && rg_EQ(y, pt.y) && rg_EQ(z, pt.z) )
        return rg_TRUE;
    else 
        return rg_FALSE;
}

rg_FLAG rg_Point3D::operator!=(const rg_Point3D &pt) const
{
    if ( rg_EQ(x, pt.x) && rg_EQ(y, pt.y) && rg_EQ(z, pt.z) )
        return rg_FALSE;
    else 
        return rg_TRUE;
}

////  Friend Function  ////
rg_Point3D operator -(const rg_Point3D &pt)
{
    return rg_Point3D(-pt.x, -pt.y, -pt.z);
}

rg_Point3D operator *(const rg_REAL &n, const rg_Point3D &pt)
{
    return pt*n;
}

rg_Point3D operator *(const rg_Point3D &pt1, const rg_Point3D &pt2) 
{
    return pt1.crossProduct(pt2);
}

rg_REAL operator %(const rg_Point3D &pt1, const rg_Point3D &pt2) 
{
    return pt1.innerProduct(pt2);
}

rg_Point3D getUnitVector(const rg_Point3D &pt)
{
	rg_Point3D unitPt;

	unitPt = pt.getUnitVector();

	return unitPt;
}



bool rg_Point3D::areThreePointsColinear(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_REAL& res)
{
    rg_REAL length[3] = { pt1.distance( pt2 ), pt1.distance( pt3 ), pt2.distance( pt3 ) };

    //  normal of center plane for a linear, parabolic, hyperbolic, or elliptic edge.
    if (    rg_GT(length[0]+length[1], length[2], res) 
         && rg_GT(length[0]+length[2], length[1], res) 
         && rg_GT(length[1]+length[2], length[0], res) )  {
        return false;
    }
    else {
        return true;
    }
}


