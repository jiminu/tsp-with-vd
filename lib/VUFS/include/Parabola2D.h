#ifndef _PARABOLA2D_H
#define _PARABOLA2D_H

#include "rg_Line2D.h"

class Parabola2D
{
private:
	rg_Point2D m_focus;
	rg_Line2D  m_directrix;	
	static const rg_INT NUM_COEFF_PARABOLA_EQ = 6;
	rg_REAL    m_coefficient[NUM_COEFF_PARABOLA_EQ]; // Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
													// m_coefficient[5] = A, m_coefficient[4] = B, m_coefficient[3] = C, 
													// m_coefficient[2] = D, m_coefficient[1] = E, m_coefficient[0] = F  

	void       compute_coefficients_of_parabola_equation();
	void       initialize_coefficients_of_parabola_equation();

public:
	Parabola2D();	
	Parabola2D(const rg_Point2D& focus, const rg_Line2D& directrix);
	Parabola2D(const Parabola2D& parabola);
	Parabola2D& operator=(const Parabola2D& parabola);

	rg_Point2D get_focus() const;
	rg_Line2D  get_directrix() const;
	rg_INT     get_coefficients_of_parabola_equation(double coefficient[]);
	rg_REAL    get_rotation_angle_from_positive_X_axis_for_axis_of_symmetry() const;
	void       rotate_axis_of_symmetry(const rg_REAL& angle, const rg_Point2D& pointOfRotation);
	void       set(const rg_Point2D& focus, const rg_Line2D& directrix);
	void set_focus(const rg_Point2D& focus);
	void set_directrix(const rg_Line2D& directrix);

    rg_Point2D get_tangent_vector(const rg_Point2D& pt);
    rg_Point2D get_passing_point(const rg_Point2D& sp, const rg_Point2D& ep);
};

inline void Parabola2D::compute_coefficients_of_parabola_equation()
{
	rg_REAL a, b, c;
	m_directrix.get_coefficients_of_implicit_form_of_line_equation(a, b, c);
	rg_REAL u = m_focus.getX();
	rg_REAL v = m_focus.getY();
	
	rg_REAL a_squared = a*a;
	rg_REAL b_squared = b*b;

	m_coefficient[5] = -b_squared / (a_squared + b_squared);
	m_coefficient[4] = (2.0 * a * b) / (a_squared + b_squared);
	m_coefficient[3] = -a_squared / (a_squared + b_squared);
	m_coefficient[2] = (2.0 * a * c) / (a_squared + b_squared) + 2.0 * u;
	m_coefficient[1] = (2.0 * b * c) / (a_squared + b_squared) + 2.0 * v;
	m_coefficient[0] = (c*c) / (a_squared + b_squared) - u*u - v*v;
}

inline void Parabola2D::initialize_coefficients_of_parabola_equation()
{
	for (rg_INT i = 0; i < NUM_COEFF_PARABOLA_EQ; i++)
	{
		m_coefficient[i] = DBL_MAX;
	}
}


inline rg_Point2D Parabola2D::get_focus() const
{
	return m_focus;
}

inline rg_Line2D Parabola2D::get_directrix() const
{
	return m_directrix;
}

inline rg_INT Parabola2D::get_coefficients_of_parabola_equation(double coefficient[])
{
	if (m_coefficient[0] == DBL_MAX)
		compute_coefficients_of_parabola_equation();

	for (rg_INT i = 0; i < NUM_COEFF_PARABOLA_EQ; i++)
	{
		coefficient[i] = m_coefficient[i];
	}
	return NUM_COEFF_PARABOLA_EQ;
}

inline void Parabola2D::set_focus(const rg_Point2D & focus)
{
	m_focus = focus;
}

inline void Parabola2D::set_directrix(const rg_Line2D& directrix)
{
	m_directrix = directrix;
}


#endif