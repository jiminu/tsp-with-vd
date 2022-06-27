#include "Ellipse2D.h"
#include "rg_TMatrix2D.h"
#include "rg_QuarticPolynomial.h"
#include "rg_GeoFunc.h"

Ellipse2D::Ellipse2D() 
{
	initialize_coefficients_of_ellipse_equation();
}

Ellipse2D::~Ellipse2D() 
{}

Ellipse2D::Ellipse2D(const rg_Point2D& center, const rg_REAL& semiMajorAxisLength, const rg_REAL& semiMinorAxisLength, const rg_REAL& angleOfMajorAxisWithXAxis)
{
	set(center, semiMajorAxisLength, semiMinorAxisLength, angleOfMajorAxisWithXAxis);
	initialize_coefficients_of_ellipse_equation();
}

Ellipse2D::Ellipse2D(const Ellipse2D& ellipse)
{
	set(ellipse.m_center, ellipse.m_semiMajorAxisLength, ellipse.m_semiMinorAxisLength, ellipse.m_angleOfMajorAxisWithXAxis, ellipse.m_coefficient);
}

Ellipse2D& Ellipse2D::operator=(const Ellipse2D& ellipse)
{
	if (this != &ellipse)
	{
		set(ellipse.m_center, ellipse.m_semiMajorAxisLength, ellipse.m_semiMinorAxisLength, ellipse.m_angleOfMajorAxisWithXAxis, ellipse.m_coefficient);
	}
	return *this;
}

void Ellipse2D::evaluate_points(const rg_INT numPoints, list<rg_Point2D>& samplePoints) const
{
	rg_REAL angleIncrement = ANGLE_360_DEGREE / numPoints;

	for (rg_REAL angle = angleIncrement; angle <= ANGLE_360_DEGREE; angle+=angleIncrement)
	{
		samplePoints.push_back(evaluate_point(angle));
	}
}

rg_Point2D Ellipse2D::evaluate_point(const rg_REAL& angleParameter) const
{
	rg_TMatrix2D transformMat;
	transformMat.translate(m_center);
	transformMat.rotate(m_angleOfMajorAxisWithXAxis, m_center);	
	rg_Point2D pointOnEllipse(m_semiMajorAxisLength*cos(angleParameter), m_semiMinorAxisLength*sin(angleParameter));
	return transformMat * pointOnEllipse;
}

bool Ellipse2D::compute_perpendicular_footprint_of_point_onto_ellipse(const rg_Point2D& givenPoint, rg_Point2D& footOnEllipse) const
{
	// assume that an ellipse is placed at canonical position.
	rg_REAL a = m_semiMajorAxisLength;
	rg_REAL b = m_semiMinorAxisLength;
	// givenPoint should be translated by ellipse center
	rg_TMatrix2D transformMat;
	transformMat.translate(-m_center);
	transformMat.rotate(-m_angleOfMajorAxisWithXAxis);	
	rg_Point2D point = transformMat * givenPoint;
	rg_REAL u = point.getX();
	rg_REAL v = point.getY();
	//rg_REAL u = givenPoint.getX() - m_center.getX();
	//rg_REAL v = givenPoint.getY() - m_center.getY();

	// we can define a quartic polynomial from the condition that a vector starting from a foot (x, y) to a given givenPoint (u, v) is parallel to the normal vector at the foot
	rg_REAL* coeff = new rg_REAL[5];
	//coeff[0] = -pow(a, 0.4e1) * pow(b, 0.4e1) + b * b * v * v * pow(a, 0.4e1) + a * a * u * u * pow(b, 0.4e1);
	//coeff[1] = -0.2e1 * pow(a, 0.4e1) * b * b - 0.2e1 * a * a * pow(b, 0.4e1) + 0.2e1 * a * a * u * u * b * b + 0.2e1 * b * b * v * v * a * a;
	//coeff[2] = -pow(a, 0.4e1) - 0.4e1 * a * a * b * b + a * a * u * u - pow(b, 0.4e1) + b * b * v * v;
	//coeff[3] = -2 * a * a - 2 * b * b;
	//coeff[4] = -1.0;

	// optimized by Maple version 18
	rg_REAL t1, t2, t3, t4, t6, t7, t9, t10;
	t1 = a * a;
	t2 = t1 * t1;
	t3 = b * b;
	t4 = t3 * t3;
	t6 = v * v;
	t7 = t3 * t6;
	t9 = u * u;
	t10 = t1 * t9;
	coeff[0] = t10 * t4 - t2 * t4 + t7 * t2;
	coeff[1] = -2 * t1 * t4 + 2 * t7 * t1 + 2 * t10 * t3 - 2 * t2 * t3;
	coeff[2] = -4 * t1 * t3 + t10 - t2 - t4 + t7;
	coeff[3] = -2 * t1 - 2 * t3;
	coeff[4] = -1;

	// we should find the minimum positive real root
	// divide all coefficients by an coefficient of absolute largest value

	//rg_REAL absoluteLargestCoefficient = -DBL_MAX;
	//rg_INT  indexOfAbsoluteLargestCoefficient = -1;
	//for (rg_INT i = 0; i <= 4; i++)
	//{
	//	rg_REAL currentAbsoluteCoefficient = rg_ABS(coeff[i]);
	//	if (currentAbsoluteCoefficient > absoluteLargestCoefficient)
	//	{
	//		absoluteLargestCoefficient = currentAbsoluteCoefficient;
	//		indexOfAbsoluteLargestCoefficient = i;
	//	}
	//}

	//for (rg_INT i = 0; i <= 4; i++)
	//{
	//	if (i != indexOfAbsoluteLargestCoefficient)
	//		coeff[i] = coeff[i] / coeff[indexOfAbsoluteLargestCoefficient];
	//}
	//coeff[indexOfAbsoluteLargestCoefficient] = 1.0;

	rg_QuarticPolynomial quarticPolynomial(coeff);

	//rg_ComplexNumber* roots = quarticPolynomial.solve();
	//rg_ComplexNumber* roots = quarticPolynomial.solve_OLD();
	//rg_Polynomial quarticPolynomial(4, coeff);
	rg_ComplexNumber* roots = quarticPolynomial.solve();

	rg_REAL smallestPostiveRoot = DBL_MAX;

	for (rg_INT i = 0; i < 4; i++)
	{
		if (roots[i].isPureRealNumber())
		{
			rg_REAL currentRoot = roots[i].getRealNumber();
			if (currentRoot > 0 && currentRoot < smallestPostiveRoot)
				smallestPostiveRoot = currentRoot;
		}
	}

	if (roots != rg_NULL)
		delete[] roots;

	if (coeff != rg_NULL)
		delete[] coeff;

	//if (smallestPostiveRoot > 0.0)
	//{
	//	// foot on the ellipse in canonical position
	//	footOnEllipse.setX(a*a*u / (smallestPostiveRoot + a*a));
	//	footOnEllipse.setY(b*b*v / (smallestPostiveRoot + b*b));
	//	// transform foot back
	//	rg_TMatrix2D transformMat;
	//	transformMat.rotate_major_axis(m_angleOfMajorAxisWithXAxis);
	//	transformMat.translate(m_center);
	//	footOnEllipse = transformMat * footOnEllipse;
	//	return true;
	//}
	//else
	//	return false;

	// foot on the ellipse in canonical position
	footOnEllipse.setX(a*a*u / (smallestPostiveRoot + a*a));
	footOnEllipse.setY(b*b*v / (smallestPostiveRoot + b*b));
	// transform foot back
	rg_TMatrix2D backTransformMat;	
	backTransformMat.rotate(m_angleOfMajorAxisWithXAxis);
	backTransformMat.translate(m_center);	
	footOnEllipse = backTransformMat * footOnEllipse;

	// determine_if_footprint_is_IN_ON_OUT_of_ellipse() ON: tolerance
	// Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
	// m_coefficient[5] = A, m_coefficient[4] = B, m_coefficient[3] = C, 
	// m_coefficient[2] = D, m_coefficient[1] = E, m_coefficient[0] = F

	const_cast<Ellipse2D&>(*this).compute_coefficients_of_ellipse_equation();
	
	rg_REAL xCoordinateOfFootprint = footOnEllipse.getX();
	rg_REAL yCoordinateOfFootprint = footOnEllipse.getY();

	//rg_REAL evaluationOfFootprint = (xCoordinateOfFootprint*xCoordinateOfFootprint) / (a*a) + (yCoordinateOfFootprint*yCoordinateOfFootprint) / (b*b);
	rg_REAL evaluationOfFootprint = m_coefficient[5] * xCoordinateOfFootprint*xCoordinateOfFootprint + m_coefficient[4] * xCoordinateOfFootprint*yCoordinateOfFootprint
		                          + m_coefficient[3] * yCoordinateOfFootprint*yCoordinateOfFootprint + m_coefficient[2] * xCoordinateOfFootprint + m_coefficient[1] * yCoordinateOfFootprint + m_coefficient[0];

	if (rg_ZERO(evaluationOfFootprint))
		return true;
	else
	{
		rg_Point2D intersection[2];
		rg_Line2D oneLinePassingGivenPointAndFootprint(givenPoint, footOnEllipse);
		rg_INT numberOfIntersection = 0;
		numberOfIntersection = rg_GeoFunc::compute_intersection_between_ellipse_and_line(*this, oneLinePassingGivenPointAndFootprint, intersection);

		if (numberOfIntersection == 1)
		{
			footOnEllipse = intersection[0];
			return true;
		}
		else if (numberOfIntersection == 2)
		{
			//rg_REAL xCoordinateOfIntersection1 = intersection[0].getX();
			//rg_REAL yCoordinateOfIntersection1 = intersection[0].getY();
			//rg_REAL evaluationOfIntersection1 = (xCoordinateOfIntersection1*xCoordinateOfIntersection1) / (a*a) + (yCoordinateOfIntersection1*yCoordinateOfIntersection1) / (b*b);
			//rg_REAL xCoordinateOfIntersection2 = intersection[1].getX();
			//rg_REAL yCoordinateOfIntersection2 = intersection[1].getY();
			//rg_REAL evaluationOfIntersection2 = (xCoordinateOfIntersection2*xCoordinateOfIntersection2) / (a*a) + (yCoordinateOfIntersection2*yCoordinateOfIntersection2) / (b*b);
			rg_REAL distanceBetweenGivenPointAndIntersection1 = givenPoint.distance(intersection[0]);
			rg_REAL distanceBetweenGivenPointAndIntersection2 = givenPoint.distance(intersection[1]);
			if(distanceBetweenGivenPointAndIntersection1 < distanceBetweenGivenPointAndIntersection2)
				footOnEllipse = intersection[0];
			else
				footOnEllipse = intersection[1];
			return true;
		}
		else
			return false;
	}

	/*
	determine_if_footprint_is_IN_ON_OUT_of_ellipse() ON: tolerance
	if(ON...)
	{
	return true;
	}
	else
	{
	define oneLinePassingTwoPoints(GivenPoint, initialFootprint)
	intersect_ellipse_with_line()
	choose shorter one
	}
	*/
}


rg_Point2D Ellipse2D::compute_tangent_vector_at_this_point_for_marching_direction(const rg_Point2D& givenPoint, const bool& bLeftFace, const rg_Point2D& marchingPoint)
{
	rg_REAL slopeOfGivenPoint = compute_slope_at_this_point(givenPoint);
	rg_Point2D guideVec = marchingPoint - givenPoint;
	rg_Point2D tangentVectorOnEllipse;

	switch (bLeftFace)
	{
	case true:
	{
		if (rg_POS(guideVec.operator*(tangentVectorOnEllipse)))
			tangentVectorOnEllipse.setPoint(1.0, slopeOfGivenPoint);
		else
			tangentVectorOnEllipse.setPoint(-1.0, -slopeOfGivenPoint);
	}
	break;
	case false:
	{
		if (rg_NEG(guideVec.operator*(tangentVectorOnEllipse)))
			tangentVectorOnEllipse.setPoint(1.0, slopeOfGivenPoint);
		else
			tangentVectorOnEllipse.setPoint(-1.0, -slopeOfGivenPoint);
	}
	break;
	}

	return tangentVectorOnEllipse;
}
