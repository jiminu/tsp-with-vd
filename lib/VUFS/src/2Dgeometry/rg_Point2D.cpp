#include <math.h>
#include "rg_RelativeOp.h"
#include "rg_Point2D.h"

#include "rg_Line2D.h"
#include "rg_Parabola2D.h"
#include "rg_GeoFunc.h"


//////////	Constructor & Destructor	//////////
rg_Point2D::rg_Point2D()
{
	x  = 0.0;
	y  = 0.0;
}

rg_Point2D::rg_Point2D( const rg_Point2D& point )
{
	x  = point.x;
	y  = point.y;
}
/*
rg_Point2D::rg_Point2D( const rg_Point3D& point )
{
x  = point.getX();
	y  = point.getY();
}
*/
rg_Point2D::rg_Point2D( const rg_REAL &px, const rg_REAL &py ) 
{
	x = px;
	y = py;
}

rg_Point2D::~rg_Point2D()
{
}

//////////	Get functions	//////////
/*
// following 4 functions are defined in "inline function" (Nov. 23, 2017 by Youngsong)
rg_REAL rg_Point2D::getX() const
{
	return x;
}

rg_REAL rg_Point2D::getY() const
{
	return y;
}
*/

//////////	Set functions	//////////
/*
// following 4 functions are defined in "inline function" (Nov. 23, 2017 by Youngsong)
void rg_Point2D::setX( const rg_REAL &px )
{
	x = px;
}

void rg_Point2D::setY( const rg_REAL &py )
{
	y = py;
}

void rg_Point2D::setPoint(const rg_Point2D &point)
{
	x  = point.x;
	y  = point.y;
}

void rg_Point2D::setPoint(const rg_REAL &px, const rg_REAL &py)
{
	x  = px;
	y  = py;
}
*/



//////////	Operation & Calculation	//////////

// following 5 functions are defined in "inline function" (June 29, 2017 by Joonghyun)
/*
rg_REAL rg_Point2D::magnitude() const
{
	return sqrt(x*x + y*y);
}

rg_REAL rg_Point2D::magnitudeSquare() const
{
	return (x*x + y*y);
}

rg_Point2D rg_Point2D::getUnitVector() const
{
    return *this/magnitude();
}

rg_REAL rg_Point2D::distance(const rg_Point2D &point) const
{
	return sqrt( (x-point.x)*(x-point.x) +
				 (y-point.y)*(y-point.y) );
}

rg_REAL rg_Point2D::squaredDist(const rg_Point2D& point) const
{
	return (x - point.x)*(x - point.x) +
		   (y - point.y)*(y - point.y);
}
*/

//////////	Operater Overloading	//////////
rg_Point2D& rg_Point2D::operator =( const rg_Point2D &point )
{
	x = point.x;
	y = point.y;

	return *this;
}

/*
rg_Point2D& rg_Point2D::operator =( const rg_Point3D &point )
{
	x = point.getX();
	y = point.getY();

	return *this;
}
*/
rg_INT rg_Point2D::operator ==( const rg_Point2D &point ) const 
{
	// 10e-6? •?„ë§? ?•´?„ ?œ?‹¤. 10e-15ê¹Œì???Š” ?•  ?•„?š”ê°? ?—†?‹¤.
	if ( rg_EQ(x, point.x) && rg_EQ(y, point.y))
		return rg_TRUE;

	else return rg_FALSE;
}

rg_Point2D rg_Point2D::operator +( const rg_Point2D &point ) const
{
	return rg_Point2D(x+point.x, y+point.y);
}
/*
rg_Point2D rg_Point2D::operator +( const rg_Point3D &point ) const
{
	return rg_Point2D(x+point.getX(), y+point.getY());
}
*/
rg_Point2D& rg_Point2D::operator +=( const rg_Point2D &point )
{
	x += point.x;
	y += point.y;

	return *this;
}

rg_Point2D rg_Point2D::operator -( const rg_Point2D &point ) const
{
	return rg_Point2D(x-point.x, y-point.y);
}

rg_Point2D& rg_Point2D::operator -=( const rg_Point2D &point )
{
	x -= point.x;
	y -= point.y;

	return *this;
}

// The following is to define a cross product.
// symbol : *
rg_REAL rg_Point2D::operator *( const rg_Point2D &point ) const
{
	return x*point.y - y*point.x;
}

rg_Point2D rg_Point2D::operator *( const rg_REAL &n ) const
{
	return rg_Point2D(n*x, n*y);
}

rg_Point2D rg_Point2D::operator /( const rg_REAL &n ) const
{
	return rg_Point2D(x/n, y/n);
}

rg_Point2D& rg_Point2D::operator =( const rg_REAL &n )
{
	x = n;
	y = n;
	
	return *this;
}


//  Jang, Tae-Bum  inserted the following operator at 1997.7.22.
//                                          for dot product.
rg_REAL rg_Point2D::operator %( const rg_Point2D &point) const
{
    return ( x*point.getX() + y*point.getY() );
}
 


//friend function
rg_Point2D operator *( const rg_REAL n, const rg_Point2D& point )
{
	return point * n;
}

rg_Point2D operator -( const rg_Point2D &point)
{
	return rg_Point2D(-point.x, -point.y);
}

/*
rg_REAL angleFromVec1toVec2(const rg_Point2D& vector1, const rg_Point2D& vector2)
{
	rg_Point2D vec1 = vector1.getUnitVector();
	rg_Point2D vec2 = vector2.getUnitVector();
	rg_Point2D pt2pt( vec1 - vec2 );
	rg_REAL length = pt2pt.magnitude();

	rg_REAL cosine = (2. - length*length)/(2.);


    //17.10.17. cysong : added - START
    if (rg_EQ(cosine, 1.0))
    {
        cosine = 1.0;
    } 
    else if (rg_EQ(cosine, -1.0))
    {
        cosine = -1.0;
    }
    //17.10.17. cysong : added - END


	if( vec1 * vec2 > 0.0 )  //less than PI (cross product)
	{
		return acos( cosine );
	}
	else
	{
		return 2.*rg_PI - acos( cosine );
	}
}
*/



bool        rg_Point2D::calculateBisector(const rg_Point2D &pt, rg_Line2D& bisector, const rg_REAL& dist) const
{
    bool    isBisectorCalculationComplete = true;

    if (isEqual(pt)) {
        isBisectorCalculationComplete = false;
    }
    else {
        rg_Line2D lineThisToPt(*this, pt);

        double A, B, C;
        lineThisToPt.convertIntoImplicitForm(A, B, C);

        double x_pt1 = x;
        double y_pt1 = y;
        double x_pt2 = pt.x;
        double y_pt2 = pt.y;

        double x_mid = (x_pt1 + x_pt2) / 2.0;
        double y_mid = (y_pt1 + y_pt2) / 2.0;

        double x_startPtOfLn = x_mid - dist*A;
        double y_startPtOfLn = y_mid - dist*B;

        double x_endPtOfLn = x_mid + dist*A;
        double y_endPtOfLn = y_mid + dist*B;

        bisector.setSP(rg_Point2D(x_startPtOfLn, y_startPtOfLn));
        bisector.setEP(rg_Point2D(x_endPtOfLn, y_endPtOfLn));

        isBisectorCalculationComplete = true;
    }

    return isBisectorCalculationComplete;
}



bool        rg_Point2D::calculateBisector(const rg_Line2D &line, rg_Parabola2D& bisector, const rg_REAL& dist) const
{
    bool    isBisectorCalculationComplete = true;

    if (line.does_contain(*this)) {
        isBisectorCalculationComplete = false;
    }
    else {
        bisector.setFocus(*this);
        bisector.setDirectrix(line);

        isBisectorCalculationComplete = true;
    }

    return isBisectorCalculationComplete;
}



rg_REAL angleFromVec1toVec2(const rg_Point2D& vector1, const rg_Point2D& vector2)
{
    return rg_GeoFunc::angleFromVec1toVec2(vector1, vector2, resNeg10);
}

#ifdef PYVORONOI
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
namespace py = pybind11;

using Point2D = rg_Point2D;

void init_Point2D(py::module& m) {
	py::class_<Point2D, std::shared_ptr<Point2D>>(m, "Point2D")
		.def(py::init([]() { return Point2D(); }))
		.def(py::init([](const double& x, const double& y) { return Point2D(x, y); }))
		.def(py::init([](const Point2D& point) { return Point2D(point); }))
		.def_property("x", &Point2D::getX, &Point2D::setX)
		.def_property("y", &Point2D::getY, &Point2D::setY)
		.def("get_x", &Point2D::getX)
		.def("get_y", &Point2D::getY)
		.def("set_x", &Point2D::setX)
		.def("set_y", &Point2D::setY)    
		.def(py::self == py::self)
		// .def(py::self != py::self)
		;
}
#endif