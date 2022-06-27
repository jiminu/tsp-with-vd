#ifndef _ELLIPSOID3D_H
#define _ELLIPSOID3D_H

#include <stdlib.h>


#include "Sphere.h"
#include <list>
using namespace std;

class Ellipsoid3D
{
private:
	rg_Point3D  m_center;
	rg_REAL     m_a;    // axis corresponding to X-axis (not necessarily larger than m_semiMinorAxisLength
	rg_REAL     m_b;    // axis corresponding to Y-axis
    rg_REAL     m_c;    // axis corresponding to Z-axis

public:
    Ellipsoid3D();
    ~Ellipsoid3D();
    Ellipsoid3D( const rg_Point3D& center, const rg_REAL& a, const rg_REAL& b, const rg_REAL& c );
    Ellipsoid3D( const Ellipsoid3D& ellipsoid );
    Ellipsoid3D& operator=( const Ellipsoid3D& ellipsoid );

    rg_Point3D  get_center() const;
    rg_REAL     get_a() const;
    rg_REAL     get_b() const;
    rg_REAL     get_c() const;

    void        set_center( const rg_Point3D& center );
    void        set_a( const rg_REAL& a );
    void        set_b( const rg_REAL& b );
    void        set_c( const rg_REAL& c );

	rg_Point3D  compute_normal_vector_of_point_on_ellipsoid(const rg_Point3D& givenPoint) const;

	rg_REAL     compute_perpendicular_footprint_of_point_onto_ellipsoid(const rg_Point3D & givenPoint, rg_Point3D & footOnEllipsoid) const;
	bool        contain(const rg_Point3D& point) const;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// inline functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline rg_Point3D Ellipsoid3D::get_center() const { 
    return m_center; 
}

inline rg_REAL Ellipsoid3D::get_a() const {
    return m_a;
}

inline rg_REAL Ellipsoid3D::get_b() const {
    return m_b;
}

inline rg_REAL Ellipsoid3D::get_c() const {
    return m_c;
}

inline void Ellipsoid3D::set_center( const rg_Point3D& center ) {
    m_center = center;
}

inline void Ellipsoid3D::set_a( const rg_REAL& a ) {
    m_a = a;
}

inline void Ellipsoid3D::set_b( const rg_REAL& b ) {
    m_b = b;
}

inline void Ellipsoid3D::set_c( const rg_REAL& c ) {
    m_c = c;
}

inline rg_Point3D Ellipsoid3D::compute_normal_vector_of_point_on_ellipsoid(const rg_Point3D & givenPoint) const
{
	rg_REAL xCoord = 2.0 * (givenPoint.getX() - m_center.getX());
	rg_REAL yCoord = 2.0 * (givenPoint.getY() - m_center.getY());
	rg_REAL zCoord = 2.0 * (givenPoint.getZ() - m_center.getZ());

	rg_Point3D normalVec(xCoord / (m_a*m_a), yCoord / (m_b*m_b), zCoord / (m_c*m_c));
	
	return normalVec;
}

inline bool Ellipsoid3D::contain(const rg_Point3D& point) const
{
	rg_REAL x = point.getX();
	rg_REAL y = point.getY();
	rg_REAL z = point.getZ();

	if ((x*x) / (m_a*m_a) + (y*y) / (m_b*m_b) + (z*z) / (m_c*m_c) <= 1.0)
		return true;
	else
		return false;
}

#endif