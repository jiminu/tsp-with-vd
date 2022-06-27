#ifndef _ELLIPSE2D_H
#define _ELLIPSE2D_H

//#ifdef _WIN32
//#include <afxwin.h>
//#endif
#include <stdlib.h>


#include "rg_Circle2D.h"
//#include "rg_dList.h"
#include <list>
using namespace std;

class Ellipse2D
{
private:
	rg_Point2D m_center;
	rg_REAL    m_semiMajorAxisLength; // axis corresponding to X-axis (not necessarily larger than m_semiMinorAxisLength
	rg_REAL    m_semiMinorAxisLength; // axis corresponding to Y-axis
	rg_REAL    m_angleOfMajorAxisWithXAxis;
	static const rg_INT NUM_COEFF_ELLIPSE_EQ = 6;
	rg_REAL    m_coefficient[NUM_COEFF_ELLIPSE_EQ]; // Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
	                               // m_coefficient[5] = a, m_coefficient[4] = b, m_coefficient[3] = c, 
	                               // m_coefficient[2] = d, m_coefficient[1] = e, m_coefficient[0] = f  

	void       compute_coefficients_of_ellipse_equation();
	void       initialize_coefficients_of_ellipse_equation();


public:
	Ellipse2D();
	~Ellipse2D();
	Ellipse2D(const rg_Point2D& center, const rg_REAL& semiMajorAxisLength, const rg_REAL& semiMinorAxisLength, const rg_REAL& angleOfMajorAxisWithXAxis);
	Ellipse2D(const Ellipse2D& ellipse);
	Ellipse2D& operator=(const Ellipse2D& ellipse);

	rg_Point2D get_centerPt() const;
	rg_REAL    get_angle_of_major_axis_with_X_axis() const;
	void       get_semi_axis_length(rg_REAL& semiMajorAxisLength, rg_REAL& semiMinorAxisLength) const;
	rg_REAL    get_semi_major_axis_length() const;
	rg_REAL    get_semi_minor_axis_length() const;
	void       get(rg_Point2D& center, rg_REAL& semiMajorAxisLength, rg_REAL& semiMinorAxisLength, rg_REAL& angleOfMajorAxisWithXAxis) const;
	void       set_centerPt(const rg_Point2D& center);
	void       set(const rg_Point2D& center, const rg_REAL& semiMajorAxisLength, const rg_REAL& semiMinorAxisLength, const rg_REAL& angleOfMajorAxisWithXAxis);
	void       set(const rg_Point2D& center, const rg_REAL& semiMajorAxisLength, const rg_REAL& semiMinorAxisLength, const rg_REAL& angleOfMajorAxisWithXAxis, const double coefficient[]);
	void       set_axis_length(const rg_REAL& semiMajorAxisLength, const rg_REAL& semiMinorAxisLength);
	void       evaluate_points(const rg_INT numPoints, list<rg_Point2D>& samplePoints) const;
	rg_Point2D evaluate_point(const rg_REAL& angleParameter) const;
	rg_INT     get_coefficients_of_ellipse_equation(double coefficient[]);
	rg_INT     get_num_coefficients_of_ellipse_equation() const;
	void       rotate_major_axis(const rg_REAL& angle);
	void       translate_center(const rg_Point2D& translationAmount);
	bool       compute_perpendicular_footprint_of_point_onto_ellipse(const rg_Point2D& givenPoint, rg_Point2D& footOnEllipse) const;
	rg_Point2D compute_tangent_vector_at_this_point_for_marching_direction(const rg_Point2D& givenPoint, const bool& bLeftFace, const rg_Point2D& marchingPoint);
	rg_REAL    compute_slope_at_this_point(const rg_Point2D& givenPoint);

	rg_REAL    compute_area();


    bool        does_contain(const rg_Point2D& point);
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// inline functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline rg_Point2D Ellipse2D::get_centerPt() const { return m_center; }
inline rg_REAL    Ellipse2D::get_angle_of_major_axis_with_X_axis() const { return m_angleOfMajorAxisWithXAxis; }
inline void       Ellipse2D::get_semi_axis_length(rg_REAL& semiMajorAxisLength, rg_REAL& semiMinorAxisLength) const
{
	semiMajorAxisLength = m_semiMajorAxisLength;
	semiMinorAxisLength = m_semiMinorAxisLength;
}

inline rg_REAL Ellipse2D::get_semi_major_axis_length() const { return m_semiMajorAxisLength; }
inline rg_REAL Ellipse2D::get_semi_minor_axis_length() const { return m_semiMinorAxisLength; }


inline void       Ellipse2D::get(rg_Point2D& center, rg_REAL& semiMajorAxisLength, rg_REAL& semiMinorAxisLength, rg_REAL& angleOfMajorAxisWithXAxis) const
{
	center = m_center;
	semiMajorAxisLength = m_semiMajorAxisLength;
	semiMinorAxisLength = m_semiMinorAxisLength;
	angleOfMajorAxisWithXAxis = m_angleOfMajorAxisWithXAxis;
}

inline void       Ellipse2D::set_centerPt(const rg_Point2D& center) { m_center = center; }

inline void Ellipse2D::set(const rg_Point2D& center,
	                       const rg_REAL& semiMajorAxisLength,
	                       const rg_REAL& semiMinorAxisLength,
	                       const rg_REAL& angleOfMajorAxisWithXAxis)
{
	m_center = center;
	m_semiMajorAxisLength = semiMajorAxisLength;
	m_semiMinorAxisLength = semiMinorAxisLength;
	m_angleOfMajorAxisWithXAxis = angleOfMajorAxisWithXAxis;
}

inline void Ellipse2D::set(const rg_Point2D& center, 
	                       const rg_REAL& semiMajorAxisLength, 
	                       const rg_REAL& semiMinorAxisLength, 
	                       const rg_REAL& angleOfMajorAxisWithXAxis, 
	                       const double coefficient[])
{
	set(center, semiMajorAxisLength, semiMinorAxisLength, angleOfMajorAxisWithXAxis);
	for (rg_INT i = 0; i < NUM_COEFF_ELLIPSE_EQ; i++)
	{
		m_coefficient[i] = coefficient[i];
	}
}

inline void Ellipse2D::set_axis_length(const rg_REAL& semiMajorAxisLength, const rg_REAL& semiMinorAxisLength)
{
	m_semiMajorAxisLength = semiMajorAxisLength;
	m_semiMinorAxisLength = semiMinorAxisLength;
}

inline void Ellipse2D::compute_coefficients_of_ellipse_equation()
{
	// Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
	// A = m_coefficient[5], B = m_coefficient[4], C = m_coefficient[3], 
	// D = m_coefficient[2], E = m_coefficient[1], F = m_coefficient[0]
	double semiMajorAxisLength_squared = m_semiMajorAxisLength*m_semiMajorAxisLength;
	double semiMinorAxisLength_squared = m_semiMinorAxisLength*m_semiMinorAxisLength;
	double sin_squared = sin(m_angleOfMajorAxisWithXAxis)*sin(m_angleOfMajorAxisWithXAxis);
	double cos_squared = cos(m_angleOfMajorAxisWithXAxis)*cos(m_angleOfMajorAxisWithXAxis);
	double center_x_coordinate = m_center.getX();
	double center_y_coordinate = m_center.getY();

	m_coefficient[5] = semiMajorAxisLength_squared*sin_squared + semiMinorAxisLength_squared*cos_squared;
	m_coefficient[4] = 2.0*(semiMinorAxisLength_squared - semiMajorAxisLength_squared)*sin(m_angleOfMajorAxisWithXAxis)*cos(m_angleOfMajorAxisWithXAxis);
	m_coefficient[3] = semiMajorAxisLength_squared*cos_squared + semiMinorAxisLength_squared*sin_squared;
	m_coefficient[2] = -2.0*m_coefficient[5] * center_x_coordinate - m_coefficient[4] * center_y_coordinate;
	m_coefficient[1] = -m_coefficient[4] * center_x_coordinate - 2.0*m_coefficient[3] * center_y_coordinate;
	m_coefficient[0] = m_coefficient[5] * center_x_coordinate*center_x_coordinate
		+ m_coefficient[4] * center_x_coordinate*center_y_coordinate
		+ m_coefficient[3] * center_y_coordinate*center_y_coordinate
		- semiMajorAxisLength_squared*semiMinorAxisLength_squared;
}

inline void Ellipse2D::initialize_coefficients_of_ellipse_equation()
{
	for (rg_INT i = 0; i < NUM_COEFF_ELLIPSE_EQ; i++)
	{
		m_coefficient[i] = DBL_MAX;
	}
}

inline rg_INT Ellipse2D::get_coefficients_of_ellipse_equation(double coefficient[])
{
	if (m_coefficient[0] == DBL_MAX)
		compute_coefficients_of_ellipse_equation();

	for (rg_INT i = 0; i < NUM_COEFF_ELLIPSE_EQ; i++)
	{
		coefficient[i] = m_coefficient[i];
	}
	return NUM_COEFF_ELLIPSE_EQ;
}

inline rg_INT Ellipse2D::get_num_coefficients_of_ellipse_equation() const { return NUM_COEFF_ELLIPSE_EQ; }

inline void Ellipse2D::rotate_major_axis(const rg_REAL & angle)
{
	m_angleOfMajorAxisWithXAxis = m_angleOfMajorAxisWithXAxis + angle;
}

inline void Ellipse2D::translate_center(const rg_Point2D& translationAmount)
{
	m_center = m_center + translationAmount;
}

inline rg_REAL Ellipse2D::compute_slope_at_this_point(const rg_Point2D& givenPoint) 
{ 
	if (m_coefficient[0] == DBL_MAX)
		compute_coefficients_of_ellipse_equation();

	rg_REAL xCoordinate = givenPoint.getX();
	rg_REAL yCoordinate = givenPoint.getY();

	rg_REAL denominator = (m_coefficient[4] * xCoordinate + 2.0 * m_coefficient[3] * yCoordinate + m_coefficient[1]);

	// if (rg_ZERO(denominator))
		// exit(1);
	//	AfxMessageBox(_T("Denominator near ZERO"));

	rg_REAL slopeOfGivenPoint = (-2.0 * m_coefficient[5] * xCoordinate - m_coefficient[4] * yCoordinate - m_coefficient[2]) / denominator;
	return slopeOfGivenPoint;
}

inline rg_REAL Ellipse2D::compute_area()
{
	return rg_PI * m_semiMajorAxisLength * m_semiMinorAxisLength;
}


inline bool   Ellipse2D::does_contain(const rg_Point2D& point) 
{
    if (m_coefficient[0] == DBL_MAX)
        compute_coefficients_of_ellipse_equation();

    // Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
    // A = m_coefficient[5], B = m_coefficient[4], C = m_coefficient[3], 
    // D = m_coefficient[2], E = m_coefficient[1], F = m_coefficient[0]
    rg_REAL x = point.getX();
    rg_REAL y = point.getY();
    rg_REAL discriminant = (m_coefficient[5]*x*x) + (m_coefficient[4]*x*y) + (m_coefficient[3]*y*y) + (m_coefficient[2]*x) + (m_coefficient[1]*y) + m_coefficient[0];

    return ( rg_LE(discriminant, 0.0) ) ? true : false;
}

#endif