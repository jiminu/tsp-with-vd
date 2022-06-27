#include "Parabola2D.h"
#include "rg_GeoFunc.h"
#include "rg_TMatrix2D.h"

Parabola2D::Parabola2D()
{
	initialize_coefficients_of_parabola_equation();
}

Parabola2D::Parabola2D(const rg_Point2D & focus, const rg_Line2D & directrix)
{
	m_focus = focus;
	m_directrix = directrix;
	initialize_coefficients_of_parabola_equation();
}

Parabola2D::Parabola2D(const Parabola2D & parabola)
{
	set(parabola.m_focus, parabola.m_directrix);
}

Parabola2D& Parabola2D::operator=(const Parabola2D & parabola)
{
	if (this == &parabola)
		return *this;
	set(parabola.m_focus, parabola.m_directrix);
}

rg_REAL Parabola2D::get_rotation_angle_from_positive_X_axis_for_axis_of_symmetry() const
{
    rg_Point2D footprint_focus;
    m_directrix.compute_perpendicular_footprint_of_point_onto_entire_line( m_focus, footprint_focus );

    rg_Point2D vec_directrix_to_focus = m_focus - footprint_focus;
    rg_Point2D yAxisDirVec(0.0, 1.0);
    double angleOfAxisOfSymmetry = angleFromVec1toVec2(yAxisDirVec, vec_directrix_to_focus);


	//rg_REAL coeffOfX, coeffOfY, coeffOfConst;
	//m_directrix.get_coefficients_of_implicit_form_of_line_equation(coeffOfX, coeffOfY, coeffOfConst);

	//rg_REAL angleOfAxisOfSymmetry;
	//if (rg_ZERO(coeffOfX))
	//	angleOfAxisOfSymmetry = rg_PI / 2.0;
	//else if (rg_ZERO(coeffOfY))
	//	angleOfAxisOfSymmetry = 0.0;
	//else
	//{
	//	rg_Point2D lineDirVec(-coeffOfY, coeffOfX);
	//	rg_Point2D xAxisDirVec(1.0, 0.0);
	//	//angleOfAxisOfSymmetry = rg_GeoFunc::calculateAngle(lineDirVec, xAxisDirVec);
 //       angleOfAxisOfSymmetry = angleFromVec1toVec2(xAxisDirVec, lineDirVec);
	//}

	return angleOfAxisOfSymmetry;
}

void Parabola2D::rotate_axis_of_symmetry(const rg_REAL & angle, const rg_Point2D& pointOfRotation)
{
	m_directrix.rotate_this_line_about_point(angle, pointOfRotation);
	rg_TMatrix2D rotationMat;
	rotationMat.rotate(angle, pointOfRotation);
	m_focus = rotationMat * m_focus;
	initialize_coefficients_of_parabola_equation();
}

void Parabola2D::set(const rg_Point2D & focus, const rg_Line2D & directrix)
{
	set_focus(focus);
	set_directrix(directrix);
}


rg_Point2D Parabola2D::get_tangent_vector(const rg_Point2D& pt) 
{
    rg_Point2D footprint;
    m_directrix.compute_perpendicular_footprint_of_point_onto_entire_line(pt, footprint);

    rg_Point2D vec1 = (m_focus - pt).getUnitVector();
    rg_Point2D vec2 = (footprint - pt).getUnitVector();

    if( rg_ZERO( (vec1+vec2).magnitude() ) )
        return rg_Point2D( -vec1.getY(), vec1.getX() );
    else
        return vec1 + vec2;
}


rg_Point2D Parabola2D::get_passing_point(const rg_Point2D& sp, const rg_Point2D& ep) 
{
    rg_Point2D footprint_sp, footprint_ep, footprint_focus;
    m_directrix.compute_perpendicular_footprint_of_point_onto_entire_line(sp, footprint_sp);
    m_directrix.compute_perpendicular_footprint_of_point_onto_entire_line(ep, footprint_ep);
    m_directrix.compute_perpendicular_footprint_of_point_onto_entire_line(m_focus, footprint_focus);

    rg_Point2D translatePt = (m_focus + footprint_focus) / 2.0;
    rg_Point2D footprint_passingPt = (footprint_sp + footprint_ep) / 2.0;

    double xCoord = footprint_focus.distance(footprint_passingPt);
    double p      = m_directrix.getDistance(m_focus) / 2.0;
    double yCoord = xCoord * xCoord / 4.0 / p;

    rg_Point2D passingPt_before_transform(xCoord, yCoord);

    rg_TMatrix2D rotationMat;
    double angle = get_rotation_angle_from_positive_X_axis_for_axis_of_symmetry();
    rotationMat.rotate(angle);

    rg_TMatrix2D translateMat;
    translateMat.translate(translatePt);

    rg_Point2D rotatedPassingPt = rotationMat * passingPt_before_transform;
    return translateMat * rotatedPassingPt;
}