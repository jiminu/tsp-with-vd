#include "rg_TMatrix2D.h"
#include <math.h>

rg_TMatrix2D::rg_TMatrix2D()
:rg_Matrix(3,3)
{
    rg_Matrix::makeIdentity();
}

rg_TMatrix2D::rg_TMatrix2D(const rg_INT &dg) :rg_Matrix(dg, dg)
{
	rg_Matrix::makeIdentity();
}

rg_TMatrix2D::rg_TMatrix2D(const rg_TMatrix2D &matrix) :rg_Matrix(matrix.row, matrix.col)
{
	for(rg_INT i = 0; i < row; i++)
		for(rg_INT j = 0; j < col; j++)
			element[i][j] = matrix.element[i][j];
}

rg_TMatrix2D::~rg_TMatrix2D()
{
}

void rg_TMatrix2D::scale(const rg_REAL &a, const rg_REAL &d)
{
    rg_TMatrix2D transform;
    
	transform.element[0][0] = a;
	transform.element[1][1] = d;

    (*this)=transform*(*this);
}

void rg_TMatrix2D::translate(const rg_REAL &m, const rg_REAL &n)
{
    rg_TMatrix2D transform;

    transform.element[0][2] = m;
	transform.element[1][2] = n;

    (*this)=transform*(*this);
}

void rg_TMatrix2D::translate(const rg_Point2D& pt)
{
    rg_TMatrix2D transform;

    transform.element[0][2] = pt.getX();
	transform.element[1][2] = pt.getY();

    (*this)=transform*(*this);
}

void rg_TMatrix2D::rotate(const rg_REAL& angle)
{
    rg_TMatrix2D transform;

	transform.element[0][0] =  cos(angle);
	transform.element[0][1] = -sin(angle);
	transform.element[1][0] =  sin(angle);
	transform.element[1][1] =  cos(angle);

    (*this)=transform*(*this);

}

void rg_TMatrix2D::rotate(const rg_REAL& angle, const rg_Point2D& center)
{
    translate(-center);
    rotate(angle);
    translate(center);
}

void rg_TMatrix2D::rotate(const rg_Point2D& fromVector, const rg_Point2D& endVector)
{
    rg_REAL magnitude=fromVector.magnitude()*endVector.magnitude();
    rg_REAL cosAngle=fromVector%endVector/magnitude;
    rg_REAL sinAngle=fromVector*endVector/magnitude;

    rg_Point2D vec(cosAngle,sinAngle);

    rotate(vec);
}

void rg_TMatrix2D::rotate(const rg_Point2D& toVector)
{
    rg_TMatrix2D transform;

    rg_Point2D unitVector=toVector.getUnitVector();
    rg_REAL cosAngle=unitVector.getX();
    rg_REAL sinAngle=unitVector.getY();

	transform.element[0][0] =  cosAngle;
	transform.element[0][1] = -sinAngle;
	transform.element[1][0] =  sinAngle;
	transform.element[1][1] =  cosAngle;

    (*this)=transform*(*this);

}

/*
void rg_TMatrix2D::rotate(const rg_REAL &m, const rg_REAL &n, const rg_REAL &theta)
{
	rg_TMatrix2D Trn(row);
	rg_TMatrix2D Rot(row);
	rg_TMatrix2D Btrn(row);
	
	Trn.translate(m,n);
	Rot.origin_rotate(theta);
	Btrn.translate(-m,-n);

	(*this) = Btrn * Rot * Trn* (*this);
}

void rg_TMatrix2D::origin_rotate(const rg_REAL &theta)
{
    rg_TMatrix2D transform;

	transform.element[0][0] =  cos(theta);
	transform.element[0][1] = -sin(theta);
	transform.element[1][0] =  sin(theta);
	transform.element[1][1] =  cos(theta);

    (*this)=transform*(*this);

}
*/
void rg_TMatrix2D::reflectToXAxis()
{
    rg_TMatrix2D transform;
    transform.element[1][1]=-1.0;
    (*this)=transform*(*this);    
}

void rg_TMatrix2D::reflectToYAxis()
{
    rg_TMatrix2D transform;
    transform.element[0][0]=-1.0;
    (*this)=transform*(*this);    

}

void rg_TMatrix2D::reflect(const rg_Line2D& line)
{
    rg_Point2D center=line.getSP();
    rg_Point2D vectorOfLine=line.evaluateVector();
    rg_Point2D xAxis(1.0,0.0);

    translate(-center);// translate center-> origin
    rotate(vectorOfLine,xAxis);// rotate vectorOfLine -> xAxis
    reflectToXAxis();
    rotate(xAxis,vectorOfLine);// rotate xAxis ->vectorOfLine
    translate(center);// translate origin-> center
}


rg_TMatrix2D& rg_TMatrix2D::operator=(const rg_TMatrix2D &matrix)
{
	rg_Matrix::removeAll();
	rg_Matrix::setSize(matrix.row, matrix.col);

	for(rg_INT i = 0; i < row; i++)
		for(rg_INT j = 0; j < col; j++)
			element[i][j] = matrix.element[i][j];

	return *(this);
}

rg_TMatrix2D rg_TMatrix2D::operator*(const rg_TMatrix2D &matrix) const
{
	rg_TMatrix2D mult(3);
	rg_INT i = 0;
	for(i = 0; i < mult.row; i++)
		for(rg_INT j = 0; j< mult.col; j++)
			mult.element[i][j] = 0.0;

	for(i = 0; i < mult.row; i++)
		for(rg_INT j = 0; j < mult.col; j++ )
			for(rg_INT k = 0; k < col ; k++)
				mult.element[i][j] += element[i][k] * matrix.element[k][j];

	return mult;
}

rg_Point2D rg_TMatrix2D::operator*(const rg_Point2D &point) const
{
	rg_REAL x = element[0][0]*point.getX()
		      +element[0][1]*point.getY()
			  +element[0][2];
	rg_REAL y = element[1][0]*point.getX()
		      +element[1][1]*point.getY()
			  +element[1][2];

	return rg_Point2D(x, y);
}


