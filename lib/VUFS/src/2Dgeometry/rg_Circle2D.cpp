/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_Circle2D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_Circle2D 
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Donguk Kim
//    START DATE  : 1999.  9.  8. 
//    
//            Copyright ⓒ CAD/CAM Lab.    	  	
//
/////////////////////////////////////////////////////////////////////
#include <math.h>
#include "rg_Circle2D.h"
#include "rg_RelativeOp.h"
#include "rg_TMatrix2D.h"


rg_Circle2D::rg_Circle2D()
{
	centerPt = 0.0;
	radius = 0.0;
}


rg_Circle2D::rg_Circle2D(const rg_Point2D& point0, const rg_Point2D& point1, const rg_Point2D& point2)
{
	//setCircleWithCCWThreePassingPoints(point0_CCW, point1_CCW, point2_CCW);
	setCircleWithThreePassingPoints(point0, point1, point2);
}


rg_Circle2D::rg_Circle2D(const rg_Point2D& center, const rg_REAL& r)
{
	centerPt = center;
	radius = r;
}

rg_Circle2D::rg_Circle2D(const rg_REAL& pointX, const rg_REAL& pointY, const rg_REAL& r)
{
	centerPt.setPoint(pointX, pointY);
	radius = r;
}

rg_Circle2D::rg_Circle2D(const rg_Circle2D& circle)
{
	centerPt = circle.getCenterPt();
	radius = circle.getRadius();
}

	
rg_Circle2D::~rg_Circle2D()
{
}

	
//get functions
//rg_Point2D rg_Circle2D::getCenterPt() const
//{
//	return centerPt;
//}
//
//rg_REAL rg_Circle2D::getRadius() const
//{
//	return radius;
//}
//
//rg_Circle2D rg_Circle2D::getCircle() const
//{
//	return *this;
//}

		
	
//set functions
void rg_Circle2D::setCenterPt(const rg_Point2D& center)
{
	centerPt = center;
}

void rg_Circle2D::setRadius(const rg_REAL& r)
{
	radius = r;
}

void rg_Circle2D::setCircle(const rg_Circle2D& circle)
{
	centerPt = circle.centerPt;
	radius = circle.radius;
}

void rg_Circle2D::setCircle(const rg_Point2D& center, const rg_REAL& r)
{
	centerPt = center;
	radius = r;
}


void rg_Circle2D::setCircleWithCCWThreePassingPoints(const rg_Point2D& point0_CCW, const rg_Point2D& point1_CCW, const rg_Point2D& point2_CCW)
{
    rg_REAL x0 = point0_CCW.getX();
    rg_REAL y0 = point0_CCW.getY();
    rg_REAL x1 = point1_CCW.getX();
    rg_REAL y1 = point1_CCW.getY();
    rg_REAL x2 = point2_CCW.getX();
    rg_REAL y2 = point2_CCW.getY();

    rg_REAL X1 = x1 - x0;
    rg_REAL X2 = x2 - x0;
    rg_REAL Y1 = y1 - y0;
    rg_REAL Y2 = y2 - y0;
    rg_REAL areaOfTriangle = 0.5 * (X1*Y2 - Y1*X2);
    rg_REAL L10 = point1_CCW.distance(point0_CCW);
    rg_REAL L20 = point2_CCW.distance(point0_CCW);

    centerPt.setX(x0 + (Y2*L10*L10 - Y1*L20*L20) / (4.0*areaOfTriangle));
    centerPt.setY(y0 + (X1*L20*L20 - X2*L10*L10) / (4.0*areaOfTriangle));
    radius = centerPt.distance(point0_CCW);
}


void rg_Circle2D::setCircleWithThreePassingPoints(const rg_Point2D& point1, const rg_Point2D& point2, const rg_Point2D& point3)
{
    rg_REAL x1 = point1.getX();
    rg_REAL y1 = point1.getY();
    rg_REAL x2 = point2.getX();
    rg_REAL y2 = point2.getY();
    rg_REAL x3 = point3.getX();
    rg_REAL y3 = point3.getY();

    centerPt.setX((x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y2 * y2 * y3 - y2 * y3 * y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) / 2.0);
    centerPt.setY(-(x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x2 * x2 * x3 - x3 * x3 * x2 + y1 * y1 * x2 - y3 * y3 * x2 - y1 * y1 * x3 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) / 2.0);
    rg_REAL x1_minus_centerX = (x1 - centerPt.getX());
    rg_REAL y1_minus_centery = (y1 - centerPt.getY());
    radius = sqrt(x1_minus_centerX*x1_minus_centerX - y1_minus_centery*y1_minus_centery);
}


rg_FLAG rg_Circle2D::isIntersectWith(const rg_Circle2D& circle, const rg_REAL& tolerance /*= rg_MATH_RES*/) const
{
	rg_Point2D center2center( circle.centerPt - centerPt ); 
	rg_REAL distance = center2center.magnitude();

    //두 원이 바깥에서 서로 떨어져 있는경우
	if( rg_LT(radius + circle.radius, distance, tolerance) )
	{
		return rg_FALSE;
	}
/*	//한 원이 다른 한 원을 포함하고 있는경우	
	if( rg_GT(rg_ABS( radius - circle.radius ), distance) )
	{
		return rg_FALSE;
	}
	//같은 원인 경우
	if( rg_EQ(radius, circle.radius) && rg_EQ(distance, 0.0) )
	{
		return rg_FALSE;
	}
*/
	return rg_TRUE;
}


rg_FLAG rg_Circle2D::IntersectWith_notTangentTo(const rg_Circle2D& circle, const rg_REAL& tolerance /*= rg_MATH_RES*/) const
{
	rg_Point2D center2center(circle.centerPt - centerPt);
	rg_REAL distance = center2center.magnitude();

	//두 원이 바깥에서 서로 떨어져 있거나 접하는 경우
	if (rg_LE(radius + circle.radius, distance, tolerance))
	{
		return rg_FALSE;
	}
	return rg_TRUE;
}


rg_FLAG rg_Circle2D::isTangentTo(const rg_Circle2D& circle, const rg_REAL& tolerance /*= rg_MATH_RES*/) const
{
	rg_Point2D center2center( circle.centerPt - centerPt ); 
	rg_REAL distance = center2center.magnitude();

	if( rg_EQ(radius + circle.radius, distance, tolerance) )
	{
		return rg_TRUE;
	}
	
	if( rg_EQ( rg_ABS(radius - circle.radius), distance, tolerance) )
	{
		return rg_TRUE;
	}

	return rg_FALSE;
}



rg_BOOL rg_Circle2D::isIncludedIn(const rg_Circle2D& circle, const rg_REAL& tolerance /*= rg_MATH_RES*/) const
{
	rg_Point2D center2center( circle.centerPt - centerPt ); 
	rg_REAL distance = center2center.magnitude();

	//한 원이 다른 한 원을 포함하고 있는경우	
	if( rg_GT( (circle.radius - radius), distance, tolerance) ) {
		return rg_TRUE;
	}
    else {
	    return rg_FALSE;
    }
}


rg_BOOL rg_Circle2D::doesContain(const rg_Point2D& point, const rg_REAL& tolerance /*= rg_MATH_RES*/) const
{
    rg_REAL distance = centerPt.distance(point);

    if ( rg_GT(radius, distance, tolerance) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_BOOL rg_Circle2D::contain(const rg_Circle2D& circle, const rg_REAL& tolerance) const
{
    rg_Point2D center2center(circle.centerPt - centerPt);
    rg_REAL distance = center2center.magnitude();

    if (radius >= (circle.radius + distance) - tolerance)
        return rg_TRUE;
    else
        return rg_FALSE;

}


rg_REAL rg_Circle2D::ratioOfOverlappedCircle(const rg_Circle2D& circle) const
{
    rg_REAL distBetweenCenters = centerPt.distance(circle.centerPt);

    if (rg_POS(radius + circle.radius - distBetweenCenters))
    {
        rg_REAL smallRadius = -1.0;
        if (radius > circle.radius)
            smallRadius = circle.radius;
        else
            smallRadius = radius;

        return (radius + circle.radius - distBetweenCenters) / (2.0 *smallRadius);
    }
    else
        return 0.0;
}


rg_Point2D* rg_Circle2D::getIntersectPt(const rg_Circle2D& circle) const
{
/*				 *	<-- intersect point
				/|\
			   / |  \
			  /	 |h   \
			 /_c_|______\
		center  foot   circle.center
*/	
	rg_Point2D center2center( circle.centerPt - centerPt ); 
	rg_REAL distance = center2center.magnitude();
	rg_Point2D foot;
	
	//두 원이 서로 떨어져 있는경우
	if( rg_LT(radius + circle.radius, distance) )	
	{
		return rg_NULL;
	}
	
	//한 원이 다른 한 원을 포함하고 있는경우	
	if( rg_GT(rg_ABS( radius - circle.radius ), distance) )
	{
		return rg_NULL;
	}

	//같은 원인 경우
	if( rg_EQ(radius, circle.radius) && rg_EQ(distance, 0.0) )
	{
		return rg_NULL;
	}

	rg_REAL c=0, h=0;
	c = (distance*distance + radius*radius - 
		 circle.radius*circle.radius)/(2*distance);
	
	if(radius*radius - c*c > 0)
		h = sqrt(radius*radius - c*c);
	
	rg_Point2D temp( -1*center2center.getY(), center2center.getX() );
	rg_Point2D* intersectPt= new rg_Point2D[2];  
		
	foot = centerPt + c*( center2center.getUnitVector() );
	intersectPt[0] = foot + temp.getUnitVector()*h;
	intersectPt[1] = foot - temp.getUnitVector()*h;
	
	return intersectPt;
}
	
rg_Circle2D& rg_Circle2D::operator=(const rg_Circle2D& circle)
{
	if( this == &circle ) {
		return *this;
	}
	
	centerPt = circle.centerPt;
	radius = circle.radius;
	
	return *this;
}



rg_BOOL rg_Circle2D::operator==(const rg_Circle2D& circle) const
{
    if (centerPt == circle.centerPt && radius == circle.radius) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}




//  geometric operators..
rg_REAL rg_Circle2D::distance(const rg_Point2D& point) const
{
    return (centerPt.distance(point) - radius);
}



rg_REAL rg_Circle2D::distance(const rg_Circle2D& circle) const
{	 
    rg_REAL dist = 0;
	rg_REAL distC2C = centerPt.distance( circle.centerPt );
    
    if ( rg_GT( distC2C, (radius+circle.radius) ) ) {
        dist = distC2C - (radius+circle.radius);
    }
    else {
        dist = (radius+circle.radius) - distC2C;
    }

    return dist;
}



rg_Circle2D rg_Circle2D::computeMinTangentCircle(const rg_Circle2D& circle)
{
	rg_REAL distBetweenCenters = centerPt.distance(circle.centerPt);
	rg_REAL sumOfRadii = radius + circle.radius;
	rg_REAL radiusOfMinTangentCircle = (distBetweenCenters - sumOfRadii) * 0.5;
	
	rg_Point2D vec = circle.centerPt - centerPt;
	vec = vec.getUnitVector();
	rg_Point2D centerOfMinTangentCircle = centerPt + (radius + radiusOfMinTangentCircle) * vec;

	rg_Circle2D tangentCircle(centerOfMinTangentCircle, radiusOfMinTangentCircle);
	return tangentCircle;
}



rg_Circle2D computeMinTangentCircle(const rg_Point2D& point1, const rg_Point2D& point2)
{
    rg_REAL diameter = point1.distance( point2 );
    rg_Point2D center = (point1 + point2)*0.5;

    return rg_Circle2D(center, diameter*0.5);
}



rg_INT computeCircleTangentTo3CirclesOutside( const rg_Circle2D& circle1, 
                                              const rg_Circle2D& circle2, 
                                              const rg_Circle2D& circle3, 
                                              rg_Circle2D* tangentCircle)
{
    rg_Circle2D circle[2];
    rg_REAL     minCircleX;
    rg_REAL     minCircleY;
    rg_REAL     minCircleRadius;
    if ( rg_LT(circle1.radius, circle2.radius) )  {
        if ( rg_LT(circle1.radius, circle3.radius) )  {
            circle[0] = circle2;
            circle[1] = circle3;
            minCircleX      = circle1.centerPt.getX();
            minCircleY      = circle1.centerPt.getY();
            minCircleRadius = circle1.radius;
        }
        else  {
            circle[0] = circle1;
            circle[1] = circle2;
            minCircleX      = circle3.centerPt.getX();
            minCircleY      = circle3.centerPt.getY();
            minCircleRadius = circle3.radius;            
        }
    }
    else  {
        if ( rg_LT(circle2.radius, circle3.radius) )  {
            circle[0] = circle1;
            circle[1] = circle3;
            minCircleX      = circle2.centerPt.getX();
            minCircleY      = circle2.centerPt.getY();
            minCircleRadius = circle2.radius;
        }
        else  {
            circle[0] = circle1;
            circle[1] = circle2;
            minCircleX      = circle3.centerPt.getX();
            minCircleY      = circle3.centerPt.getY();
            minCircleRadius = circle3.radius;            
        }
    }


	rg_REAL x1 = circle[0].centerPt.getX() - minCircleX;
	rg_REAL y1 = circle[0].centerPt.getY() - minCircleY;
	rg_REAL r1 = circle[0].radius          - minCircleRadius;

	rg_REAL x2 = circle[1].centerPt.getX() - minCircleX;
	rg_REAL y2 = circle[1].centerPt.getY() - minCircleY;
	rg_REAL r2 = circle[1].radius          - minCircleRadius;
    

    //system of equations
	//eq1: x1*x + y1*y + r1*r = rhs1
	//eq2: x2*x + y2*y + r2*r = rhs2
	//eq4: x*x + y*y + z*z = r*r

	rg_REAL rhs1 = (x1*x1+y1*y1-r1*r1)*0.5;
	rg_REAL rhs2 = (x2*x2+y2*y2-r2*r2)*0.5;

    rg_REAL determinant1 = x1*y2 - x2*y1;
	if( !rg_ZERO(determinant1) )  {
        //system of equations
	    //eq1: x1*x + y1*y = rhs1 - r1*r 
	    //eq2: x2*x + y2*y = rhs2 - r2*r 
        //
        //  x = a_x + b_x*r
        //  y = a_y + b_y*r
        //
        //  eq3 <- x, y
        //  a_r*r^2 + b_r*r + c_r = 0

		rg_REAL a_x = (y2*rhs1 - y1*rhs2)/determinant1;
        rg_REAL b_x = (y1*r2   - y2*r1)  /determinant1;
		rg_REAL a_y = (x1*rhs2 - x2*rhs1)/determinant1;
        rg_REAL b_y = (x2*r1   - x1*r2)  /determinant1;

		rg_REAL a_r = b_x*b_x + b_y*b_y - 1.0;
        rg_REAL b_r = 2.0*(a_x*b_x + a_y*b_y);
        rg_REAL c_r = a_x*a_x + a_y*a_y;

		rg_REAL det = b_r*b_r - 4.*a_r*c_r; //b^2 - 4ac
		if(det < 0) 
			return 0;

        rg_REAL sqrtDet = sqrt(det); 
		rg_REAL radius1 = (-b_r + sqrtDet)/(2.* a_r);
		rg_REAL radius2 = (-b_r - sqrtDet)/(2.* a_r);


        rg_REAL x;
        rg_REAL y;
        if( radius1 > 0 && radius2 > 0 )  {

			x = a_x + b_x*radius1 + minCircleX;
			y = a_y + b_y*radius1 + minCircleY;
			tangentCircle[0].centerPt.setPoint(x, y);
            tangentCircle[0].radius = radius1-minCircleRadius;

			x = a_x + b_x*radius2 + minCircleX;
			y = a_y + b_y*radius2 + minCircleY;
			tangentCircle[1].centerPt.setPoint(x, y);
            tangentCircle[1].radius = radius2-minCircleRadius;

			return 2;
		}
		else if( radius1 > 0 )  {

			x = a_x + b_x*radius1 + minCircleX;
			y = a_y + b_y*radius1 + minCircleY;
			tangentCircle[0].centerPt.setPoint(x, y);
            tangentCircle[0].radius = radius1-minCircleRadius;

			return 1;
		}
		else if( radius2 > 0 )  {

			x = a_x + b_x*radius2 + minCircleX;
			y = a_y + b_y*radius2 + minCircleY;
			tangentCircle[0].centerPt.setPoint(x, y);
            tangentCircle[0].radius = radius2-minCircleRadius;

			return 1;
		}
		else  {
			return 0;
		}    
    }


	rg_REAL determinant2 = x1*r2 - x2*r1;
    if( !rg_ZERO(determinant2) )  {
        //system of equations
	    //eq1: x1*x + r1*r = rhs1 - y1*y 
	    //eq2: x2*x + r2*r = rhs2 - y2*y 
        //
        //  x = a_x + b_x*y
        //  r = a_r + b_r*y
        //
        //  eq3 <- x, r
        //  a_y*y^2 + b_y*y + c_y = 0

		rg_REAL a_x = (r2*rhs1 - r1*rhs2)/determinant2;
        rg_REAL b_x = (r1*y2   - r2*y1)  /determinant2;
		rg_REAL a_r = (x1*rhs2 - x2*rhs1)/determinant2;
        rg_REAL b_r = (x2*y1   - x1*y2)  /determinant2;

		rg_REAL a_y = b_x*b_x + b_r*b_r - 1.0;
        rg_REAL b_y = 2.0*(a_x*b_x + a_r*b_r);
        rg_REAL c_y = a_x*a_x + a_r*a_r;

		rg_REAL det = b_y*b_y - 4.*a_y*c_y; //b^2 - 4ac
		if(det < 0) 
			return 0;

        rg_REAL sqrtDet = sqrt(det); 
		rg_REAL y1 = (-b_y + sqrtDet)/(2.* a_y);
		rg_REAL y2 = (-b_y - sqrtDet)/(2.* a_y);
        
		rg_REAL radius1 = a_r + b_r*y1;
		rg_REAL radius2 = a_r + b_r*y2;
        
        rg_REAL x;
        rg_REAL y;
        if( radius1 > 0 && radius2 > 0 )  {

			x = a_x + b_x*radius1 + minCircleX;
			y = y1                + minCircleY;
			tangentCircle[0].centerPt.setPoint(x, y);
            tangentCircle[0].radius = radius1-minCircleRadius;

			x = a_x + b_x*radius2 + minCircleX;
			y = y2                + minCircleY;
			tangentCircle[1].centerPt.setPoint(x, y);
            tangentCircle[1].radius = radius2-minCircleRadius;

			return 2;
		}
		else if( radius1 > 0 )  {

			x = a_x + b_x*radius1 + minCircleX;
			y = y1                + minCircleY;
			tangentCircle[0].centerPt.setPoint(x, y);
            tangentCircle[0].radius = radius1-minCircleRadius;

			return 1;
		}
		else if( radius2 > 0 )  {

			x = a_x + b_x*radius2 + minCircleX;
			y = y2                + minCircleY;
			tangentCircle[0].centerPt.setPoint(x, y);
            tangentCircle[0].radius = radius2-minCircleRadius;

			return 1;
		}
		else  {
			return 0;
		}        
    }

    
    rg_REAL determinant3 = y1*r2 - y2*r1;
    if( !rg_ZERO(determinant3) )  {
        //system of equations
	    //eq1: y1*y + r1*r = rhs1 - x1*x 
	    //eq2: y2*y + r2*r = rhs2 - x2*x 
        //
        //  y = a_y + b_y*x
        //  r = a_r + b_r*x
        //
        //  eq3 <- y, r
        //  a_x*x^2 + b_x*x + c_x = 0

		rg_REAL a_y = (r2*rhs1 - r1*rhs2)/determinant2;
        rg_REAL b_y = (r1*x2   - r2*x1)  /determinant2;
		rg_REAL a_r = (y1*rhs2 - y2*rhs1)/determinant2;
        rg_REAL b_r = (y2*x1   - y1*x2)  /determinant2;

		rg_REAL a_x = b_y*b_y + b_r*b_r - 1.0;
        rg_REAL b_x = 2.0*(a_y*b_y + a_r*b_r);
        rg_REAL c_x = a_y*a_y + a_r*a_r;

		rg_REAL det = b_x*b_x - 4.*a_x*c_x; //b^2 - 4ac
		if(det < 0) 
			return 0;

        rg_REAL sqrtDet = sqrt(det); 
		rg_REAL x1 = (-b_x + sqrtDet)/(2.* a_x);
		rg_REAL x2 = (-b_x - sqrtDet)/(2.* a_x);
        
		rg_REAL radius1 = a_r + b_r*x1;
		rg_REAL radius2 = a_r + b_r*x2;
        
        rg_REAL x;
        rg_REAL y;
        if( radius1 > 0 && radius2 > 0 )  {

			x = x1                + minCircleX;
			y = a_y + b_y*radius1 + minCircleY;
			tangentCircle[0].centerPt.setPoint(x, y);
            tangentCircle[0].radius = radius1-minCircleRadius;

			x = x2                + minCircleX;
			y = a_y + b_y*radius2 + minCircleY;
			tangentCircle[1].centerPt.setPoint(x, y);
            tangentCircle[1].radius = radius2-minCircleRadius;

			return 2;
		}
		else if( radius1 > 0 )  {

			x = x1                + minCircleX;
			y = a_y + b_y*radius1 + minCircleY;
			tangentCircle[0].centerPt.setPoint(x, y);
            tangentCircle[0].radius = radius1-minCircleRadius;

			return 1;
		}
		else if( radius2 > 0 )  {

			x = x2                + minCircleX;
			y = a_y + b_y*radius2 + minCircleY;
			tangentCircle[0].centerPt.setPoint(x, y);
            tangentCircle[0].radius = radius2-minCircleRadius;

			return 1;
		}
		else  {
			return 0;
		}        
    }

    return 0;
}




rg_INT computeCircleTangentTo2CirclesGivenRadius( const rg_Circle2D& circle1, 
                                                  const rg_Circle2D& circle2, 
                                                  const rg_REAL&     givenRadius,
                                                  rg_Circle2D* tangentCircle )
{
    rg_INT numTangentCircles = 0;

    rg_REAL r  = givenRadius;

    rg_REAL x1 = circle1.getCenterPt().getX();
    rg_REAL y1 = circle1.getCenterPt().getY();
    rg_REAL r1 = circle1.getRadius();
    rg_REAL x2 = circle2.getCenterPt().getX();
    rg_REAL y2 = circle2.getCenterPt().getY();
    rg_REAL r2 = circle2.getRadius();

    rg_REAL squareR1 = (r1 + r)*(r1 + r);
    rg_REAL squareR2 = (r2 + r)*(r2 + r);

    //  Given r, find the center (x, y) of tangent circles.
    //  (x - x1)^2 + (y - y1)^2 = (r1 + r)^2    -----   (1)
    //  (x - x2)^2 + (y - y2)^2 = (r2 + r)^2    -----   (2)

    //  (1) - (2)
    //  2(x2 - x1)x + 2(y2 - y1)y = (r2 + r)^2 - (r1 + r)^2 
    //                               + x2^2 - x1^2 + y2^2 - y1^2  -----   (3)
    rg_REAL a = 2.0*(x2 - x1);
    rg_REAL b = 2.0*(y2 - y1);
    rg_REAL c = squareR1 - squareR2 + x2*x2 - x1*x1 + y2*y2 - y1*y1;
    if ( rg_NZERO( a ) ) {
        //  x = ax*y + bx   -----   (4)
        rg_REAL ax = -b/a;
        rg_REAL bx = c/a;

        //  (4) -> (1)
        //  (ax^2 + 1)*y^2 + (2*ax*bx - 2*ax*x1 - 2*y1)*y
        //  + ( (bx - x1)^2 + y1^2 - (r1 + r)^2) = 0
        rg_REAL a_y = ( ax*ax + 1.0 );
        rg_REAL b_y = ( 2.0*ax*bx - 2.0*ax*x1 - 2.0*y1 );
        rg_REAL c_y = ( (bx - x1)*(bx - x1) + y1*y1 - squareR1 );
        
		rg_REAL det = b_y*b_y - 4.*a_y*c_y; //b^2 - 4ac
		if(det < 0) 
			return 0;

        rg_REAL sqrtDet   = sqrt(det); 
        rg_INT  numYRoots = rg_ZERO( sqrtDet ) ? 1: 2;
        rg_REAL ty[2] = { (-b_y + sqrtDet)/(2.* a_y), (-b_y - sqrtDet)/(2.* a_y) };

        //  ty -> (2)
        for ( rg_INT i=0; i<numYRoots; i++ ) {
            rg_REAL a_x = 1.0;
            rg_REAL b_x = -2.0*x2;
            rg_REAL c_x = x2*x2 + (ty[i] - y2)*(ty[i] - y2) - squareR2;

	        rg_REAL detX = b_x*b_x - 4.*a_x*c_x; //b^2 - 4ac
            if( detX < 0 )  {
		        continue;
            }

            rg_REAL sqrtDetX = sqrt(detX); 
            rg_INT  numXRoots = rg_ZERO( sqrtDetX ) ? 1: 2;
            rg_REAL tx[2] = { (-b_x + sqrtDetX)/(2.* a_x), (-b_x - sqrtDetX)/(2.* a_x) };

            //  tx, ty -> (1)
            //for ( rg_INT j=0; j<numYRoots; j++ ) {
			// JKKIM
			for ( rg_INT j=0; j<numXRoots; j++ ) {
				// Joonghyun June 24, 2017 ///////////////////////////////////////////////////////////////

	//          rg_REAL distC2C = sqrt( (x1-tx[j])*(x1-tx[j]) + (y1-ty[i])*(y1-ty[i]) ) - r1;
    //            if ( rg_EQ( distC2C, givenRadius ) ) {
    //                tangentCircle[numTangentCircles].setCircle( rg_Point2D(tx[j], ty[i]), givenRadius );
    //                numTangentCircles++;
    //            }
				rg_REAL distC2C1 = sqrt((x1 - tx[j])*(x1 - tx[j]) + (y1 - ty[i])*(y1 - ty[i])) - r1;
				rg_REAL distC2C2 = sqrt((x2 - tx[j])*(x2 - tx[j]) + (y2 - ty[i])*(y2 - ty[i])) - r2;
				if (rg_EQ(distC2C1, givenRadius) && rg_EQ(distC2C2, givenRadius)) {
					tangentCircle[numTangentCircles].setCircle(rg_Point2D(tx[j], ty[i]), givenRadius);
					numTangentCircles++;
				}
				//////////////////////////////////////////////////////////////////////////////////////////
            }


            //if ( rg_ZERO( sqrtDetX ) ) {
	           // rg_REAL tx1 = (-b_x)/(2.0*a_x);

            //    tangentCircle[0].setCircle( rg_Point2D(tx1, ty[i]), givenRadius );
            //    numTangentCircles = 1;
            //}
            //else {
	           // rg_REAL tx1 = (-b_x + sqrtDetX)/(2.* a_x);
	           // rg_REAL tx2 = (-b_x - sqrtDetX)/(2.* a_x);

            //    tangentCircle[0].setCircle( rg_Point2D(tx1, ty[i]), givenRadius );
            //    tangentCircle[1].setCircle( rg_Point2D(tx2, ty[i]), givenRadius );

            //    rg_REAL dist[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
            //    for ( rg_INT i=0; i<2; i++ ) {
            //        dist[i][0] = circle1.distance( tangentCircle[i] );
            //        dist[i][1] = circle2.distance( tangentCircle[i] );
            //    }

            //    numTangentCircles = 2;
            //}
            ////break;
        }
    }
    else if ( rg_NZERO( b ) ) {
        //  y = ay*x + by   -----   (5)
        rg_REAL ay = -a/b;
        rg_REAL by = c/b;

        //  (5) -> (1)
        //  (ay^2 + 1)*x^2 + (2*ay*by - 2*ay*y1 - 2*x1)*x
        //  + ( (by - y1)^2 + x1^2 - (r1 + r)^2) = 0

        rg_REAL a_x = ( ay*ay + 1.0 );
        rg_REAL b_x = ( 2.0*ay*by - 2.0*ay*y1 -2.0*x1 );
        rg_REAL c_x = ( (by - y1)*(by - y1) + x1*x1 - squareR1 );
        
		rg_REAL det = b_x*b_x - 4.*a_x*c_x; //b^2 - 4ac
		if(det < 0) 
			return 0;

        rg_REAL sqrtDet = sqrt(det); 
		rg_REAL tx1 = (-b_x + sqrtDet)/(2.* a_x);
		rg_REAL tx2 = (-b_x - sqrtDet)/(2.* a_x);

        rg_INT numXRoots = rg_ZERO( sqrtDet ) ? 1: 2;
        rg_REAL tx[2] = { tx1, tx2 };

        //  tx -> (2)
        for ( rg_INT i=0; i<numXRoots; i++ ) {
            rg_REAL a_y = 1.0;
            rg_REAL b_y = -2.0*y2;
            rg_REAL c_y = y2*y2 + (tx[i] - x2)*(tx[i] - x2) - squareR2;

	        rg_REAL detY = b_y*b_y - 4.*a_y*c_y; //b^2 - 4ac
            if( detY < 0 )  {
		        continue;
            }

            rg_REAL sqrtDetY = sqrt(detY); 
            rg_INT  numYRoots = rg_ZERO( sqrtDetY ) ? 1: 2;
            rg_REAL ty[2] = { (-b_y + sqrtDetY)/(2.* a_y), (-b_y - sqrtDetY)/(2.* a_y) };

            //  tx, ty -> (1)
            for ( rg_INT j=0; j<numYRoots; j++ ) {
				// Joonghyun June 24, 2017 ///////////////////////////////////////////////////////////////
                
				//rg_REAL distC2C = sqrt( (x1-tx[i])*(x1-tx[i]) + (y1-ty[j])*(y1-ty[j]) ) - r1;
                //if ( rg_EQ( distC2C, givenRadius ) ) {
                //    tangentCircle[numTangentCircles].setCircle( rg_Point2D(tx[i], ty[j]), givenRadius );
                //    numTangentCircles++;
                //}

				rg_REAL distC2C1 = sqrt((x1 - tx[i])*(x1 - tx[i]) + (y1 - ty[j])*(y1 - ty[j])) - r1;
				rg_REAL distC2C2 = sqrt((x2 - tx[i])*(x2 - tx[i]) + (y2 - ty[j])*(y2 - ty[j])) - r2;
				if (rg_EQ(distC2C1, givenRadius) && rg_EQ(distC2C2, givenRadius)) {
					tangentCircle[numTangentCircles].setCircle(rg_Point2D(tx[j], ty[i]), givenRadius);
					numTangentCircles++;
				}
				///////////////////////////////////////////////////////////////////////////////////////////
            }
            //if ( rg_ZERO( sqrtDetY ) ) {
	           // rg_REAL ty1 = (-b_y + sqrtDetY)/(2.* a_y);
            //    tangentCircle[0].setCircle( rg_Point2D(tx[i], ty1), givenRadius );
            //    numTangentCircles = 1;
            //}
            //else {
	           // rg_REAL ty1 = (-b_y + sqrtDetY)/(2.* a_y);
	           // rg_REAL ty2 = (-b_y - sqrtDetY)/(2.* a_y);

            //    tangentCircle[0].setCircle( rg_Point2D(tx[i], ty1), givenRadius );
            //    tangentCircle[1].setCircle( rg_Point2D(tx[i], ty2), givenRadius );

            //    rg_REAL dist[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
            //    for ( rg_INT i=0; i<2; i++ ) {
            //        dist[i][0] = circle1.distance( tangentCircle[i] );
            //        dist[i][1] = circle2.distance( tangentCircle[i] );
            //    }

            //    numTangentCircles = 2;
            //}
        }
    }
    else {
        numTangentCircles = 0;
    }


    return numTangentCircles;
}

rg_INT computeCircleTangentToOneCircleAndItsCircularContainerGivenRadius(const rg_Circle2D & circle1, const rg_Circle2D & circularContainer, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle)
{
	rg_INT numTangentCircles = 0;

	rg_REAL r = givenRadius;

	rg_REAL x1 = circle1.getCenterPt().getX();
	rg_REAL y1 = circle1.getCenterPt().getY();
	rg_REAL r1 = circle1.getRadius();
	rg_REAL x2 = circularContainer.getCenterPt().getX();
	rg_REAL y2 = circularContainer.getCenterPt().getY();
	rg_REAL r2 = circularContainer.getRadius();

	rg_REAL squareR1 = (r1 + r)*(r1 + r);
	rg_REAL squareR2 = (r2 - r)*(r2 - r);

	//  Given r, find the center (x, y) of tangent circles.
	//  (x - x1)^2 + (y - y1)^2 = (r1 + r)^2    -----   (1)
	//  (x - x2)^2 + (y - y2)^2 = (r2 + r)^2    -----   (2)

	//  (1) - (2)
	//  2(x2 - x1)x + 2(y2 - y1)y = (r1 + r)^2 - (r2 + r)^2 
	//                               + x2^2 - x1^2 + y2^2 - y1^2  -----   (3)
	rg_REAL a = 2.0*(x2 - x1);
	rg_REAL b = 2.0*(y2 - y1);
	rg_REAL c = squareR1 - squareR2 + x2*x2 - x1*x1 + y2*y2 - y1*y1;
	if (rg_NZERO(a)) {
		//  x = ax*y + bx   -----   (4)
		rg_REAL ax = -b / a;
		rg_REAL bx = c / a;

		//  (4) -> (1)
		//  (ax^2 + 1)*y^2 + (2*ax*bx - 2*ax*x1 - 2*y1)*y
		//  + ( (bx - x1)^2 + y1^2 - (r1 + r)^2) = 0
		rg_REAL a_y = (ax*ax + 1.0);
		rg_REAL b_y = (2.0*ax*bx - 2.0*ax*x1 - 2.0*y1);
		rg_REAL c_y = ((bx - x1)*(bx - x1) + y1*y1 - squareR1);

		rg_REAL det = b_y*b_y - 4.*a_y*c_y; //b^2 - 4ac
		if (det < 0)
			return 0;

		rg_REAL sqrtDet = sqrt(det);
		rg_INT  numYRoots = rg_ZERO(sqrtDet) ? 1 : 2;
		rg_REAL ty[2] = { (-b_y + sqrtDet) / (2.* a_y), (-b_y - sqrtDet) / (2.* a_y) };

		//  ty -> (2)
		for (rg_INT i = 0; i<numYRoots; i++) {
			rg_REAL a_x = 1.0;
			rg_REAL b_x = -2.0*x2;
			rg_REAL c_x = x2*x2 + (ty[i] - y2)*(ty[i] - y2) - squareR2;

			rg_REAL detX = b_x*b_x - 4.*a_x*c_x; //b^2 - 4ac
			if (detX < 0) {
				continue;
			}

			rg_REAL sqrtDetX = sqrt(detX);
			rg_INT  numXRoots = rg_ZERO(sqrtDetX) ? 1 : 2;
			rg_REAL tx[2] = { (-b_x + sqrtDetX) / (2.* a_x), (-b_x - sqrtDetX) / (2.* a_x) };

			//  tx, ty -> (1)
			//for ( rg_INT j=0; j<numYRoots; j++ ) {
			// JKKIM
			for (rg_INT j = 0; j<numXRoots; j++) {
				rg_REAL distC2C = sqrt((x1 - tx[j])*(x1 - tx[j]) + (y1 - ty[i])*(y1 - ty[i])) - r1;
				if (rg_EQ(distC2C, givenRadius)) {
					tangentCircle[numTangentCircles].setCircle(rg_Point2D(tx[j], ty[i]), givenRadius);
					numTangentCircles++;
				}
			}


			//if ( rg_ZERO( sqrtDetX ) ) {
			// rg_REAL tx1 = (-b_x)/(2.0*a_x);

			//    tangentCircle[0].setCircle( rg_Point2D(tx1, ty[i]), givenRadius );
			//    numTangentCircles = 1;
			//}
			//else {
			// rg_REAL tx1 = (-b_x + sqrtDetX)/(2.* a_x);
			// rg_REAL tx2 = (-b_x - sqrtDetX)/(2.* a_x);

			//    tangentCircle[0].setCircle( rg_Point2D(tx1, ty[i]), givenRadius );
			//    tangentCircle[1].setCircle( rg_Point2D(tx2, ty[i]), givenRadius );

			//    rg_REAL dist[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
			//    for ( rg_INT i=0; i<2; i++ ) {
			//        dist[i][0] = circle1.distance( tangentCircle[i] );
			//        dist[i][1] = circularContainer.distance( tangentCircle[i] );
			//    }

			//    numTangentCircles = 2;
			//}
			////break;
		}
	}
	else if (rg_NZERO(b)) {
		//  y = ay*x + by   -----   (5)
		rg_REAL ay = -a / b;
		rg_REAL by = c / b;

		//  (5) -> (1)
		//  (ay^2 + 1)*x^2 + (2*ay*by - 2*ay*y1 - 2*x1)*x
		//  + ( (by - y1)^2 + x1^2 - (r1 + r)^2) = 0

		rg_REAL a_x = (ay*ay + 1.0);
		rg_REAL b_x = (2.0*ay*by - 2.0*ay*y1 - 2.0*x1);
		rg_REAL c_x = ((by - y1)*(by - y1) + x1*x1 - squareR1);

		rg_REAL det = b_x*b_x - 4.*a_x*c_x; //b^2 - 4ac
		if (det < 0)
			return 0;

		rg_REAL sqrtDet = sqrt(det);
		rg_REAL tx1 = (-b_x + sqrtDet) / (2.* a_x);
		rg_REAL tx2 = (-b_x - sqrtDet) / (2.* a_x);

		rg_INT numXRoots = rg_ZERO(sqrtDet) ? 1 : 2;
		rg_REAL tx[2] = { tx1, tx2 };

		//  tx -> (2)
		for (rg_INT i = 0; i<numXRoots; i++) {
			rg_REAL a_y = 1.0;
			rg_REAL b_y = -2.0*y2;
			rg_REAL c_y = y2*y2 + (tx[i] - x2)*(tx[i] - x2) - squareR2;

			rg_REAL detY = b_y*b_y - 4.*a_y*c_y; //b^2 - 4ac
			if (detY < 0) {
				continue;
			}

			rg_REAL sqrtDetY = sqrt(detY);
			rg_INT  numYRoots = rg_ZERO(sqrtDetY) ? 1 : 2;
			rg_REAL ty1 = (-b_y + sqrtDetY) / (2.* a_y);
			rg_REAL ty2 = (-b_y - sqrtDetY) / (2.* a_y);
			rg_REAL ty[2] = { (-b_y + sqrtDetY) / (2.* a_y), (-b_y - sqrtDetY) / (2.* a_y) };

			//  tx, ty -> (1)
			for (rg_INT j = 0; j<numYRoots; j++) {
				rg_REAL distC2C = sqrt((x1 - tx[i])*(x1 - tx[i]) + (y1 - ty[j])*(y1 - ty[j])) - r1;
				if (rg_EQ(distC2C, givenRadius)) {
					tangentCircle[numTangentCircles].setCircle(rg_Point2D(tx[i], ty[j]), givenRadius);
					numTangentCircles++;
				}
			}
			//if ( rg_ZERO( sqrtDetY ) ) {
			// rg_REAL ty1 = (-b_y + sqrtDetY)/(2.* a_y);
			//    tangentCircle[0].setCircle( rg_Point2D(tx[i], ty1), givenRadius );
			//    numTangentCircles = 1;
			//}
			//else {
			// rg_REAL ty1 = (-b_y + sqrtDetY)/(2.* a_y);
			// rg_REAL ty2 = (-b_y - sqrtDetY)/(2.* a_y);

			//    tangentCircle[0].setCircle( rg_Point2D(tx[i], ty1), givenRadius );
			//    tangentCircle[1].setCircle( rg_Point2D(tx[i], ty2), givenRadius );

			//    rg_REAL dist[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
			//    for ( rg_INT i=0; i<2; i++ ) {
			//        dist[i][0] = circle1.distance( tangentCircle[i] );
			//        dist[i][1] = circularContainer.distance( tangentCircle[i] );
			//    }

			//    numTangentCircles = 2;
			//}
		}
	}
	else {
		numTangentCircles = 0;
	}


	return numTangentCircles;
}


/*  
>> Song, Chanyoung commeted out this function on April 04, 2020.
rg_INT rg_Circle2D::makeCircumcircle(const rg_Circle2D& circle1, 
                        const rg_Circle2D& circle2, 
                        const rg_Circle2D& circle3,
						rg_Circle2D& result1, 
                        rg_Circle2D& result2)
{
	rg_INT numOfCC = 0;
	
	rg_ImplicitEquation tangentLine1(1);
	rg_ImplicitEquation tangentLine2(1);
	rg_Point2D z3;
	z3 = makeExteriorTangentLinesInWPlane(circle1, circle2, circle3, tangentLine1, tangentLine2);
	rg_REAL smallest = 0;
	if( z3 == circle1.getCenterPt() )
	{
		smallest = circle1.getRadius();
	}
	else if ( z3 == circle2.getCenterPt() )
	{
		smallest = circle2.getRadius();
	}
	else
	{
		smallest = circle3.getRadius();
	}
	//하나인지 두개인지 판단 (point location problem)
	rg_REAL sign1 = tangentLine1.evaluateImpEquation(0, 0);
	rg_REAL sign2 = tangentLine2.evaluateImpEquation(0, 0);
	if( sign1 * sign2 < 0 ) //alpha region
	{
		numOfCC = 1;
		if( sign1 > 0 )
		{
			result1 = transformW2Z(tangentLine1, z3);		//+ -
			result1.setRadius( result1.getRadius() - smallest );
		}
		else
		{
			result1 = transformW2Z(tangentLine2, z3);		//- +
			result1.setRadius( result1.getRadius() - smallest );
		}
	}
	else if( sign1 >0 && sign2 >0 ) //gamma region
	{
		numOfCC = 2;
		result1 = transformW2Z(tangentLine1, z3);			//+ +
		result1.setRadius( result1.getRadius() - smallest );
		result2 = transformW2Z(tangentLine2, z3);
		result2.setRadius( result2.getRadius() - smallest );
	}
	else if( sign1 < 0 && sign2 < 0 )						//- -
	{
		numOfCC = 0;
	}
	else if( rg_EQ(sign1, 0., resNeg15) && sign2 > 0 )						//0 +
	{
		numOfCC = 1;
		result1 = transformW2Z(tangentLine2, z3);
		result1.setRadius( result1.getRadius() - smallest );
	}
	else if( sign1 > 0 && rg_EQ(sign2, 0., resNeg15) )						//+ 0
	{
		numOfCC = 1;
		result1 = transformW2Z(tangentLine1, z3);
		result1.setRadius( result1.getRadius() - smallest );
	}
	else												//0 0
	{													//0 -
		numOfCC = 0;									//- 0
	}
	return numOfCC;
}
*/


bool rg_Circle2D::compute_perpendicular_footprint_of_point_onto_circle(const rg_Point2D& givenPoint, rg_Point2D& footOnCircle) const
{
    rg_REAL distanceFromGivenPointToCenter = givenPoint.distance(centerPt);
    rg_Point2D unitDirectionVector = (centerPt - givenPoint).getUnitVector();
    footOnCircle = givenPoint + (distanceFromGivenPointToCenter - radius) * unitDirectionVector;

    return true;
}


rg_INT rg_Circle2D::makeCircumcircle(const rg_Circle2D & circle1, const rg_Circle2D & circle2, const rg_Circle2D & circle3, rg_Circle2D & result1, rg_Circle2D & result2, const rg_REAL & tolerance)
{
	rg_INT numOfCC = 0;

	rg_ImplicitEquation tangentLine1(1);
	rg_ImplicitEquation tangentLine2(1);

	rg_Point2D z3;
	z3 = makeExteriorTangentLinesInWPlane(circle1, circle2, circle3, tangentLine1, tangentLine2);

	rg_REAL smallest = 0;
	if (z3 == circle1.getCenterPt())
	{
		smallest = circle1.getRadius();
	}
	else if (z3 == circle2.getCenterPt())
	{
		smallest = circle2.getRadius();
	}
	else
	{
		smallest = circle3.getRadius();
	}

	//하나인지 두개인지 판단 (point location problem)
	rg_REAL sign1 = tangentLine1.evaluateImpEquation(0, 0);
	rg_REAL sign2 = tangentLine2.evaluateImpEquation(0, 0);

	if (sign1 * sign2 < 0) //alpha region
	{
		numOfCC = 1;
		if (sign1 > 0)
		{
			result1 = transformW2Z(tangentLine1, z3);		//+ -
			result1.setRadius(result1.getRadius() - smallest);
		}
		else
		{
			result1 = transformW2Z(tangentLine2, z3);		//- +
			result1.setRadius(result1.getRadius() - smallest);
		}
	}
	else if (sign1 >0 && sign2 >0) //gamma region
	{
		numOfCC = 2;
		result1 = transformW2Z(tangentLine1, z3);			//+ +
		result1.setRadius(result1.getRadius() - smallest);

		result2 = transformW2Z(tangentLine2, z3);
		result2.setRadius(result2.getRadius() - smallest);
	}
	else if (sign1 < 0 && sign2 < 0)						//- -
	{
		numOfCC = 0;
	}
	else if (rg_EQ(sign1, 0., tolerance) && sign2 > 0)						//0 +
	{
		numOfCC = 1;
		result1 = transformW2Z(tangentLine2, z3);
		result1.setRadius(result1.getRadius() - smallest);
	}
	else if (sign1 > 0 && rg_EQ(sign2, 0., tolerance))						//+ 0
	{
		numOfCC = 1;
		result1 = transformW2Z(tangentLine1, z3);
		result1.setRadius(result1.getRadius() - smallest);
	}
	else												//0 0
	{													//0 -
		numOfCC = 0;									//- 0
	}

	return numOfCC;
}









// Enclosing circle and two contained circles always have two tangent circles.
void rg_Circle2D::calculateTangentCircles(const rg_Circle2D& enclosingCircle,
	const rg_Circle2D& circle1,
	const rg_Circle2D& circle2,
	rg_Circle2D& result1,
	rg_Circle2D& result2)
{
	rg_Circle2D bigCircle = enclosingCircle;
	rg_Circle2D baseCircle;
	rg_Circle2D midCircle;

	if (circle1.getRadius() < circle2.getRadius())
	{
		baseCircle = circle1;
		midCircle = circle2;
	}
	else
	{
		baseCircle = circle2;
		midCircle = circle1;
	}

	//shrink and enlarge 
	midCircle.setRadius(midCircle.getRadius() - baseCircle.getRadius());
	bigCircle.setRadius(bigCircle.getRadius() + baseCircle.getRadius());

	//translation
	rg_Point2D translationVector = -(baseCircle.getCenterPt());
	midCircle.setCenterPt(midCircle.getCenterPt() + translationVector);
	bigCircle.setCenterPt(bigCircle.getCenterPt() + translationVector);

	//rotation
	rg_TMatrix2D rotateMatrix;
	rg_Point2D fromVector = midCircle.getCenterPt();
	rotateMatrix.rotate(fromVector, rg_Point2D(1, 0));

	midCircle.setCenterPt(rotateMatrix * midCircle.getCenterPt());
	bigCircle.setCenterPt(rotateMatrix * bigCircle.getCenterPt());

	// The number of tangent circles is always two.
	result1 = solution1(midCircle.getCenterPt().getX(), midCircle.getRadius(),
		bigCircle.getCenterPt().getX(), bigCircle.getCenterPt().getY(),
		bigCircle.getRadius());
	result2 = solution2(midCircle.getCenterPt().getX(), midCircle.getRadius(),
		bigCircle.getCenterPt().getX(), bigCircle.getCenterPt().getY(),
		bigCircle.getRadius());

	//rotate & translate
	rg_TMatrix2D rotateMatrix2;
	rotateMatrix2.rotate(rg_Point2D(1, 0), fromVector);
	result1.setCenterPt(rotateMatrix2 * result1.getCenterPt() - translationVector);
	result2.setCenterPt(rotateMatrix2 * result2.getCenterPt() - translationVector);

	//shrink and enlarge
	result1.setRadius(result1.getRadius() - baseCircle.getRadius());
	result2.setRadius(result2.getRadius() - baseCircle.getRadius());
}


// first tangent circles of enclosing circle and two contained circles.
rg_Circle2D rg_Circle2D::solution1(const double& x1, const double& r1, const double& x2, const double& y2, const double& r2)
{
	double x, y, r;

	if (rg_EQ(r1, 0.))
	{
		x = x1 / 2.;
		double den = (2.*(r2*r2 -
			y2*y2));
		double a = -r2*r2 +
			x2*x2 + y2*y2;
		double b = -r2*r2 +
			x1*x1 - 2 * x1*x2 +
			x2*x2 + y2*y2;
		y = (r2*r2*y2 + x1*x2*y2 -
			x2*x2*y2 -
			y2*y2*y2 -
			r2*sqrt(a*b)
			);
		y = y / den;
		r = sqrt(x*x + y*y);

		return rg_Circle2D(x, y, r);
	}

	if (rg_EQ(y2, 0.))
	{
		x = (-(pow(r1, 2)*r2) -
			r1*pow(r2, 2) +
			r2*pow(x1, 2) +
			r1*pow(x2, 2)) /
			(2 * r2*x1 + 2 * r1*x2);

		y = -sqrt(pow(r2, 2) -
			2 * pow(x2, 2) +
			pow(x2, 4) / pow(r2, 2) -
			(4 * pow(r1, 4)*pow(r2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (8 * pow(r1, 3)*
				pow(r2, 3)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (4 * pow(r1, 2)*
				pow(r2, 4)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (8 * pow(r1, 2)*
				pow(r2, 2)*pow(x1, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (8 * r1*pow(r2, 3)*
				pow(x1, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (4 * pow(r2, 2)*
				pow(x1, 4)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (4 * pow(r1, 4)*
				pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (16 * pow(r1, 3)*r2*
				pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (12 * pow(r1, 2)*
				pow(r2, 2)*pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (8 * pow(r1, 2)*
				pow(x1, 2)*pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (16 * r1*r2*pow(x1, 2)*
				pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (4 * pow(x1, 4)*
				pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (12 * pow(r1, 2)*
				pow(x2, 4)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (8 * pow(r1, 3)*
				pow(x2, 4)) /
				(r2*pow(2 * r2*x1 +
					2 * r1*x2, 2)) +
					(8 * r1*pow(x1, 2)*
						pow(x2, 4)) /
						(r2*pow(2 * r2*x1 +
							2 * r1*x2, 2)) +
							(4 * pow(r1, 2)*pow(x2, 6)) /
			(pow(r2, 2)*
				pow(2 * r2*x1 + 2 * r1*x2,
					2)) -
					(4 * pow(r1, 2)*r2*x2) /
			(2 * r2*x1 + 2 * r1*x2) -
			(4 * r1*pow(r2, 2)*x2) /
			(2 * r2*x1 + 2 * r1*x2) +
			(4 * r2*pow(x1, 2)*x2) /
			(2 * r2*x1 + 2 * r1*x2) +
			(8 * r1*pow(x2, 3)) /
			(2 * r2*x1 + 2 * r1*x2) +
			(4 * pow(r1, 2)*pow(x2, 3)) /
			(r2*(2 * r2*x1 + 2 * r1*x2)) -
			(4 * pow(x1, 2)*pow(x2, 3)) /
			(r2*(2 * r2*x1 + 2 * r1*x2)) -
			(4 * r1*pow(x2, 5)) /
			(pow(r2, 2)*
			(2 * r2*x1 + 2 * r1*x2))) / 2.;

		r = sqrt(x*x + y*y);

		return rg_Circle2D(x, y, r);
	}

	x = (-4 * pow(r1, 2)*
		pow(r2, 2)*x1 -
		4 * r1*pow(r2, 3)*x1 +
		4 * pow(r2, 2)*
		pow(x1, 3) -
		4 * pow(r1, 3)*r2*x2 -
		4 * pow(r1, 2)*pow(r2, 2)*
		x2 + 4 * r1*r2*pow(x1, 2)*
		x2 + 4 * r1*r2*x1*
		pow(x2, 2) +
		4 * pow(r1, 2)*
		pow(x2, 3) +
		4 * pow(r1, 2)*x1*
		pow(y2, 2) +
		4 * r1*r2*x1*pow(y2, 2) -
		4 * pow(x1, 3)*
		pow(y2, 2) +
		4 * pow(r1, 2)*x2*
		pow(y2, 2) +
		4 * sqrt(-(pow(r1, 6)*
			pow(r2, 2)*
			pow(y2, 2)) -
			2 * pow(r1, 5)*
			pow(r2, 3)*pow(y2, 2)
			- pow(r1, 4)*
			pow(r2, 4)*pow(y2, 2)
			+ 2 * pow(r1, 4)*
			pow(r2, 2)*
			pow(x1, 2)*pow(y2, 2)
			+ 2 * pow(r1, 3)*
			pow(r2, 3)*
			pow(x1, 2)*pow(y2, 2)
			+ pow(r1, 2)*
			pow(r2, 4)*
			pow(x1, 2)*pow(y2, 2)
			- pow(r1, 2)*
			pow(r2, 2)*
			pow(x1, 4)*pow(y2, 2)
			- 2 * pow(r1, 4)*
			pow(r2, 2)*x1*x2*
			pow(y2, 2) +
			2 * pow(r1, 2)*
			pow(r2, 2)*
			pow(x1, 3)*x2*
			pow(y2, 2) +
			pow(r1, 6)*pow(x2, 2)*
			pow(y2, 2) +
			2 * pow(r1, 5)*r2*
			pow(x2, 2)*pow(y2, 2)
			+ 2 * pow(r1, 4)*
			pow(r2, 2)*
			pow(x2, 2)*pow(y2, 2)
			- 2 * pow(r1, 4)*
			pow(x1, 2)*
			pow(x2, 2)*pow(y2, 2)
			- 2 * pow(r1, 3)*r2*
			pow(x1, 2)*
			pow(x2, 2)*pow(y2, 2)
			- 2 * pow(r1, 2)*
			pow(r2, 2)*
			pow(x1, 2)*
			pow(x2, 2)*pow(y2, 2)
			+ pow(r1, 2)*
			pow(x1, 4)*
			pow(x2, 2)*pow(y2, 2)
			+ 2 * pow(r1, 4)*x1*
			pow(x2, 3)*pow(y2, 2)
			- 2 * pow(r1, 2)*
			pow(x1, 3)*
			pow(x2, 3)*pow(y2, 2)
			- pow(r1, 4)*
			pow(x2, 4)*pow(y2, 2)
			+ pow(r1, 2)*
			pow(x1, 2)*
			pow(x2, 4)*pow(y2, 2)
			+ pow(r1, 6)*
			pow(y2, 4) +
			2 * pow(r1, 5)*r2*
			pow(y2, 4) +
			2 * pow(r1, 4)*
			pow(r2, 2)*pow(y2, 4)
			- 2 * pow(r1, 4)*
			pow(x1, 2)*pow(y2, 4)
			- 2 * pow(r1, 3)*r2*
			pow(x1, 2)*pow(y2, 4)
			- 2 * pow(r1, 2)*
			pow(r2, 2)*
			pow(x1, 2)*pow(y2, 4)
			+ pow(r1, 2)*
			pow(x1, 4)*pow(y2, 4)
			+ 2 * pow(r1, 4)*x1*x2*
			pow(y2, 4) -
			2 * pow(r1, 2)*
			pow(x1, 3)*x2*
			pow(y2, 4) -
			2 * pow(r1, 4)*
			pow(x2, 2)*pow(y2, 4)
			+ 2 * pow(r1, 2)*
			pow(x1, 2)*
			pow(x2, 2)*pow(y2, 4)
			- pow(r1, 4)*
			pow(y2, 6) +
			pow(r1, 2)*pow(x1, 2)*
			pow(y2, 6))) /
			(2.*(4 * pow(r2, 2)*
				pow(x1, 2) +
				8 * r1*r2*x1*x2 +
				4 * pow(r1, 2)*
				pow(x2, 2) +
				4 * pow(r1, 2)*
				pow(y2, 2) -
				4 * pow(x1, 2)*pow(y2, 2)
				));
	y = (-(pow(r1, 2)*r2) -
		r1*pow(r2, 2) +
		r2*pow(x1, 2) +
		r1*pow(x2, 2) +
		r1*pow(y2, 2) +
		(4 * pow(r1, 2)*pow(r2, 3)*
			pow(x1, 2)) /
			(4 * pow(r2, 2)*
				pow(x1, 2) +
				8 * r1*r2*x1*x2 +
				4 * pow(r1, 2)*
				pow(x2, 2) +
				4 * pow(r1, 2)*
				pow(y2, 2) -
				4 * pow(x1, 2)*
				pow(y2, 2)) +
				(4 * r1*pow(r2, 4)*
					pow(x1, 2)) /
					(4 * pow(r2, 2)*
						pow(x1, 2) +
						8 * r1*r2*x1*x2 +
						4 * pow(r1, 2)*
						pow(x2, 2) +
						4 * pow(r1, 2)*
						pow(y2, 2) -
						4 * pow(x1, 2)*
						pow(y2, 2)) -
						(4 * pow(r2, 3)*
							pow(x1, 4)) /
							(4 * pow(r2, 2)*
								pow(x1, 2) +
								8 * r1*r2*x1*x2 +
								4 * pow(r1, 2)*
								pow(x2, 2) +
								4 * pow(r1, 2)*
								pow(y2, 2) -
								4 * pow(x1, 2)*
								pow(y2, 2)) +
								(8 * pow(r1, 3)*pow(r2, 2)*
									x1*x2) /
									(4 * pow(r2, 2)*
										pow(x1, 2) +
										8 * r1*r2*x1*x2 +
										4 * pow(r1, 2)*
										pow(x2, 2) +
										4 * pow(r1, 2)*
										pow(y2, 2) -
										4 * pow(x1, 2)*
										pow(y2, 2)) +
										(8 * pow(r1, 2)*pow(r2, 3)*
											x1*x2) /
											(4 * pow(r2, 2)*
												pow(x1, 2) +
												8 * r1*r2*x1*x2 +
												4 * pow(r1, 2)*
												pow(x2, 2) +
												4 * pow(r1, 2)*
												pow(y2, 2) -
												4 * pow(x1, 2)*
												pow(y2, 2)) -
												(8 * r1*pow(r2, 2)*
													pow(x1, 3)*x2) /
													(4 * pow(r2, 2)*
														pow(x1, 2) +
														8 * r1*r2*x1*x2 +
														4 * pow(r1, 2)*
														pow(x2, 2) +
														4 * pow(r1, 2)*
														pow(y2, 2) -
														4 * pow(x1, 2)*
														pow(y2, 2)) +
														(4 * pow(r1, 4)*r2*
															pow(x2, 2)) /
															(4 * pow(r2, 2)*
																pow(x1, 2) +
																8 * r1*r2*x1*x2 +
																4 * pow(r1, 2)*
																pow(x2, 2) +
																4 * pow(r1, 2)*
																pow(y2, 2) -
																4 * pow(x1, 2)*
																pow(y2, 2)) +
																(4 * pow(r1, 3)*pow(r2, 2)*
																	pow(x2, 2)) /
																	(4 * pow(r2, 2)*
																		pow(x1, 2) +
																		8 * r1*r2*x1*x2 +
																		4 * pow(r1, 2)*
																		pow(x2, 2) +
																		4 * pow(r1, 2)*
																		pow(y2, 2) -
																		4 * pow(x1, 2)*
																		pow(y2, 2)) -
																		(4 * pow(r1, 2)*r2*
																			pow(x1, 2)*pow(x2, 2))
		/ (4 * pow(r2, 2)*
			pow(x1, 2) +
			8 * r1*r2*x1*x2 +
			4 * pow(r1, 2)*
			pow(x2, 2) +
			4 * pow(r1, 2)*
			pow(y2, 2) -
			4 * pow(x1, 2)*
			pow(y2, 2)) -
			(4 * r1*pow(r2, 2)*
				pow(x1, 2)*pow(x2, 2))
		/ (4 * pow(r2, 2)*
			pow(x1, 2) +
			8 * r1*r2*x1*x2 +
			4 * pow(r1, 2)*
			pow(x2, 2) +
			4 * pow(r1, 2)*
			pow(y2, 2) -
			4 * pow(x1, 2)*
			pow(y2, 2)) -
			(8 * pow(r1, 2)*r2*x1*
				pow(x2, 3)) /
				(4 * pow(r2, 2)*
					pow(x1, 2) +
					8 * r1*r2*x1*x2 +
					4 * pow(r1, 2)*
					pow(x2, 2) +
					4 * pow(r1, 2)*
					pow(y2, 2) -
					4 * pow(x1, 2)*
					pow(y2, 2)) -
					(4 * pow(r1, 3)*
						pow(x2, 4)) /
						(4 * pow(r2, 2)*
							pow(x1, 2) +
							8 * r1*r2*x1*x2 +
							4 * pow(r1, 2)*
							pow(x2, 2) +
							4 * pow(r1, 2)*
							pow(y2, 2) -
							4 * pow(x1, 2)*
							pow(y2, 2)) -
							(4 * pow(r1, 2)*r2*
								pow(x1, 2)*pow(y2, 2))
		/ (4 * pow(r2, 2)*
			pow(x1, 2) +
			8 * r1*r2*x1*x2 +
			4 * pow(r1, 2)*
			pow(x2, 2) +
			4 * pow(r1, 2)*
			pow(y2, 2) -
			4 * pow(x1, 2)*
			pow(y2, 2)) -
			(4 * r1*pow(r2, 2)*
				pow(x1, 2)*pow(y2, 2))
		/ (4 * pow(r2, 2)*
			pow(x1, 2) +
			8 * r1*r2*x1*x2 +
			4 * pow(r1, 2)*
			pow(x2, 2) +
			4 * pow(r1, 2)*
			pow(y2, 2) -
			4 * pow(x1, 2)*
			pow(y2, 2)) +
			(4 * r2*pow(x1, 4)*
				pow(y2, 2)) /
				(4 * pow(r2, 2)*
					pow(x1, 2) +
					8 * r1*r2*x1*x2 +
					4 * pow(r1, 2)*
					pow(x2, 2) +
					4 * pow(r1, 2)*
					pow(y2, 2) -
					4 * pow(x1, 2)*
					pow(y2, 2)) -
					(4 * pow(r1, 3)*x1*x2*
						pow(y2, 2)) /
						(4 * pow(r2, 2)*
							pow(x1, 2) +
							8 * r1*r2*x1*x2 +
							4 * pow(r1, 2)*
							pow(x2, 2) +
							4 * pow(r1, 2)*
							pow(y2, 2) -
							4 * pow(x1, 2)*
							pow(y2, 2)) -
							(8 * pow(r1, 2)*r2*x1*x2*
								pow(y2, 2)) /
								(4 * pow(r2, 2)*
									pow(x1, 2) +
									8 * r1*r2*x1*x2 +
									4 * pow(r1, 2)*
									pow(x2, 2) +
									4 * pow(r1, 2)*
									pow(y2, 2) -
									4 * pow(x1, 2)*
									pow(y2, 2)) +
									(4 * r1*pow(x1, 3)*x2*
										pow(y2, 2)) /
										(4 * pow(r2, 2)*
											pow(x1, 2) +
											8 * r1*r2*x1*x2 +
											4 * pow(r1, 2)*
											pow(x2, 2) +
											4 * pow(r1, 2)*
											pow(y2, 2) -
											4 * pow(x1, 2)*
											pow(y2, 2)) -
											(4 * pow(r1, 3)*pow(x2, 2)*
												pow(y2, 2)) /
												(4 * pow(r2, 2)*
													pow(x1, 2) +
													8 * r1*r2*x1*x2 +
													4 * pow(r1, 2)*
													pow(x2, 2) +
													4 * pow(r1, 2)*
													pow(y2, 2) -
													4 * pow(x1, 2)*
													pow(y2, 2)) -
													(4 * r2*x1*
														sqrt(-(pow(r1, 6)*
															pow(r2, 2)*
															pow(y2, 2)) -
															2 * pow(r1, 5)*
															pow(r2, 3)*
															pow(y2, 2) -
															pow(r1, 4)*
															pow(r2, 4)*
															pow(y2, 2) +
															2 * pow(r1, 4)*
															pow(r2, 2)*
															pow(x1, 2)*
															pow(y2, 2) +
															2 * pow(r1, 3)*
															pow(r2, 3)*
															pow(x1, 2)*
															pow(y2, 2) +
															pow(r1, 2)*
															pow(r2, 4)*
															pow(x1, 2)*
															pow(y2, 2) -
															pow(r1, 2)*
															pow(r2, 2)*
															pow(x1, 4)*
															pow(y2, 2) -
															2 * pow(r1, 4)*
															pow(r2, 2)*x1*x2*
															pow(y2, 2) +
															2 * pow(r1, 2)*
															pow(r2, 2)*
															pow(x1, 3)*x2*
															pow(y2, 2) +
															pow(r1, 6)*
															pow(x2, 2)*
															pow(y2, 2) +
															2 * pow(r1, 5)*r2*
															pow(x2, 2)*
															pow(y2, 2) +
															2 * pow(r1, 4)*
															pow(r2, 2)*
															pow(x2, 2)*
															pow(y2, 2) -
															2 * pow(r1, 4)*
															pow(x1, 2)*
															pow(x2, 2)*
															pow(y2, 2) -
															2 * pow(r1, 3)*r2*
															pow(x1, 2)*
															pow(x2, 2)*
															pow(y2, 2) -
															2 * pow(r1, 2)*
															pow(r2, 2)*
															pow(x1, 2)*
															pow(x2, 2)*
															pow(y2, 2) +
															pow(r1, 2)*
															pow(x1, 4)*
															pow(x2, 2)*
															pow(y2, 2) +
															2 * pow(r1, 4)*x1*
															pow(x2, 3)*
															pow(y2, 2) -
															2 * pow(r1, 2)*
															pow(x1, 3)*
															pow(x2, 3)*
															pow(y2, 2) -
															pow(r1, 4)*
															pow(x2, 4)*
															pow(y2, 2) +
															pow(r1, 2)*
															pow(x1, 2)*
															pow(x2, 4)*
															pow(y2, 2) +
															pow(r1, 6)*
															pow(y2, 4) +
															2 * pow(r1, 5)*r2*
															pow(y2, 4) +
															2 * pow(r1, 4)*
															pow(r2, 2)*
															pow(y2, 4) -
															2 * pow(r1, 4)*
															pow(x1, 2)*
															pow(y2, 4) -
															2 * pow(r1, 3)*r2*
															pow(x1, 2)*
															pow(y2, 4) -
															2 * pow(r1, 2)*
															pow(r2, 2)*
															pow(x1, 2)*
															pow(y2, 4) +
															pow(r1, 2)*
															pow(x1, 4)*
															pow(y2, 4) +
															2 * pow(r1, 4)*x1*x2*
															pow(y2, 4) -
															2 * pow(r1, 2)*
															pow(x1, 3)*x2*
															pow(y2, 4) -
															2 * pow(r1, 4)*
															pow(x2, 2)*
															pow(y2, 4) +
															2 * pow(r1, 2)*
															pow(x1, 2)*
															pow(x2, 2)*
															pow(y2, 4) -
															pow(r1, 4)*
															pow(y2, 6) +
															pow(r1, 2)*
															pow(x1, 2)*
															pow(y2, 6))) /
															(4 * pow(r2, 2)*
																pow(x1, 2) +
																8 * r1*r2*x1*x2 +
																4 * pow(r1, 2)*
																pow(x2, 2) +
																4 * pow(r1, 2)*
																pow(y2, 2) -
																4 * pow(x1, 2)*
																pow(y2, 2)) -
																(4 * r1*x2*
																	sqrt(-(pow(r1, 6)*
																		pow(r2, 2)*
																		pow(y2, 2)) -
																		2 * pow(r1, 5)*
																		pow(r2, 3)*
																		pow(y2, 2) -
																		pow(r1, 4)*
																		pow(r2, 4)*
																		pow(y2, 2) +
																		2 * pow(r1, 4)*
																		pow(r2, 2)*
																		pow(x1, 2)*
																		pow(y2, 2) +
																		2 * pow(r1, 3)*
																		pow(r2, 3)*
																		pow(x1, 2)*
																		pow(y2, 2) +
																		pow(r1, 2)*
																		pow(r2, 4)*
																		pow(x1, 2)*
																		pow(y2, 2) -
																		pow(r1, 2)*
																		pow(r2, 2)*
																		pow(x1, 4)*
																		pow(y2, 2) -
																		2 * pow(r1, 4)*
																		pow(r2, 2)*x1*x2*
																		pow(y2, 2) +
																		2 * pow(r1, 2)*
																		pow(r2, 2)*
																		pow(x1, 3)*x2*
																		pow(y2, 2) +
																		pow(r1, 6)*
																		pow(x2, 2)*
																		pow(y2, 2) +
																		2 * pow(r1, 5)*r2*
																		pow(x2, 2)*
																		pow(y2, 2) +
																		2 * pow(r1, 4)*
																		pow(r2, 2)*
																		pow(x2, 2)*
																		pow(y2, 2) -
																		2 * pow(r1, 4)*
																		pow(x1, 2)*
																		pow(x2, 2)*
																		pow(y2, 2) -
																		2 * pow(r1, 3)*r2*
																		pow(x1, 2)*
																		pow(x2, 2)*
																		pow(y2, 2) -
																		2 * pow(r1, 2)*
																		pow(r2, 2)*
																		pow(x1, 2)*
																		pow(x2, 2)*
																		pow(y2, 2) +
																		pow(r1, 2)*
																		pow(x1, 4)*
																		pow(x2, 2)*
																		pow(y2, 2) +
																		2 * pow(r1, 4)*x1*
																		pow(x2, 3)*
																		pow(y2, 2) -
																		2 * pow(r1, 2)*
																		pow(x1, 3)*
																		pow(x2, 3)*
																		pow(y2, 2) -
																		pow(r1, 4)*
																		pow(x2, 4)*
																		pow(y2, 2) +
																		pow(r1, 2)*
																		pow(x1, 2)*
																		pow(x2, 4)*
																		pow(y2, 2) +
																		pow(r1, 6)*
																		pow(y2, 4) +
																		2 * pow(r1, 5)*r2*
																		pow(y2, 4) +
																		2 * pow(r1, 4)*
																		pow(r2, 2)*
																		pow(y2, 4) -
																		2 * pow(r1, 4)*
																		pow(x1, 2)*
																		pow(y2, 4) -
																		2 * pow(r1, 3)*r2*
																		pow(x1, 2)*
																		pow(y2, 4) -
																		2 * pow(r1, 2)*
																		pow(r2, 2)*
																		pow(x1, 2)*
																		pow(y2, 4) +
																		pow(r1, 2)*
																		pow(x1, 4)*
																		pow(y2, 4) +
																		2 * pow(r1, 4)*x1*x2*
																		pow(y2, 4) -
																		2 * pow(r1, 2)*
																		pow(x1, 3)*x2*
																		pow(y2, 4) -
																		2 * pow(r1, 4)*
																		pow(x2, 2)*
																		pow(y2, 4) +
																		2 * pow(r1, 2)*
																		pow(x1, 2)*
																		pow(x2, 2)*
																		pow(y2, 4) -
																		pow(r1, 4)*
																		pow(y2, 6) +
																		pow(r1, 2)*
																		pow(x1, 2)*
																		pow(y2, 6))) /
																		(4 * pow(r2, 2)*
																			pow(x1, 2) +
																			8 * r1*r2*x1*x2 +
																			4 * pow(r1, 2)*
																			pow(x2, 2) +
																			4 * pow(r1, 2)*
																			pow(y2, 2) -
																			4 * pow(x1, 2)*
																			pow(y2, 2))) /
																			(2.*r1*y2);

	r = sqrt(x*x + y*y);

	return rg_Circle2D(x, y, r);
}


// second tangent circles of enclosing circle and two contained circles.
rg_Circle2D rg_Circle2D::solution2(const double& x1, const double& r1, const double& x2, const double& y2, const double& r2)
{
	double x, y, r;

	if (rg_EQ(r1, 0.))
	{
		x = x1 / 2.;
		y = (pow(r2, 2)*y2 + x1*x2*y2 -
			pow(x2, 2)*y2 -
			pow(y2, 3) +
			r2*sqrt(pow(r2, 2) -
				pow(x2, 2) - pow(y2, 2))
			*sqrt(pow(r2, 2) -
				pow(x1, 2) + 2 * x1*x2 -
				pow(x2, 2) - pow(y2, 2))
			) /
			(2.*(pow(r2, 2) -
				pow(y2, 2)));
		r = sqrt(x*x + y*y);

		return rg_Circle2D(x, y, r);
	}

	if (rg_EQ(y2, 0.))
	{
		x = (-(pow(r1, 2)*r2) -
			r1*pow(r2, 2) +
			r2*pow(x1, 2) +
			r1*pow(x2, 2)) /
			(2 * r2*x1 + 2 * r1*x2);


		y = sqrt(pow(r2, 2) -
			2 * pow(x2, 2) +
			pow(x2, 4) / pow(r2, 2) -
			(4 * pow(r1, 4)*pow(r2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (8 * pow(r1, 3)*
				pow(r2, 3)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (4 * pow(r1, 2)*
				pow(r2, 4)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (8 * pow(r1, 2)*
				pow(r2, 2)*pow(x1, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (8 * r1*pow(r2, 3)*
				pow(x1, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (4 * pow(r2, 2)*
				pow(x1, 4)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (4 * pow(r1, 4)*
				pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (16 * pow(r1, 3)*r2*
				pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (12 * pow(r1, 2)*
				pow(r2, 2)*pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (8 * pow(r1, 2)*
				pow(x1, 2)*pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (16 * r1*r2*pow(x1, 2)*
				pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			+ (4 * pow(x1, 4)*
				pow(x2, 2)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (12 * pow(r1, 2)*
				pow(x2, 4)) /
			pow(2 * r2*x1 + 2 * r1*x2, 2)\
			- (8 * pow(r1, 3)*
				pow(x2, 4)) /
				(r2*pow(2 * r2*x1 + 2 * r1*x2,
					2)) +
					(8 * r1*pow(x1, 2)*
						pow(x2, 4)) /
						(r2*pow(2 * r2*x1 + 2 * r1*x2,
							2)) +
							(4 * pow(r1, 2)*pow(x2, 6)) /
			(pow(r2, 2)*
				pow(2 * r2*x1 + 2 * r1*x2, 2)
				) -
				(4 * pow(r1, 2)*r2*x2) /
			(2 * r2*x1 + 2 * r1*x2) -
			(4 * r1*pow(r2, 2)*x2) /
			(2 * r2*x1 + 2 * r1*x2) +
			(4 * r2*pow(x1, 2)*x2) /
			(2 * r2*x1 + 2 * r1*x2) +
			(8 * r1*pow(x2, 3)) /
			(2 * r2*x1 + 2 * r1*x2) +
			(4 * pow(r1, 2)*pow(x2, 3)) /
			(r2*(2 * r2*x1 + 2 * r1*x2)) -
			(4 * pow(x1, 2)*pow(x2, 3)) /
			(r2*(2 * r2*x1 + 2 * r1*x2)) -
			(4 * r1*pow(x2, 5)) /
			(pow(r2, 2)*
			(2 * r2*x1 + 2 * r1*x2))) / 2.;

		r = sqrt(x*x + y*y);

		return rg_Circle2D(x, y, r);
	}


	x = (-4 * pow(r1, 2)*
		pow(r2, 2)*x1 -
		4 * r1*pow(r2, 3)*x1 +
		4 * pow(r2, 2)*
		pow(x1, 3) -
		4 * pow(r1, 3)*r2*x2 -
		4 * pow(r1, 2)*pow(r2, 2)*
		x2 + 4 * r1*r2*pow(x1, 2)*
		x2 + 4 * r1*r2*x1*
		pow(x2, 2) +
		4 * pow(r1, 2)*
		pow(x2, 3) +
		4 * pow(r1, 2)*x1*
		pow(y2, 2) +
		4 * r1*r2*x1*pow(y2, 2) -
		4 * pow(x1, 3)*
		pow(y2, 2) +
		4 * pow(r1, 2)*x2*
		pow(y2, 2) -
		4 * sqrt(-(pow(r1, 6)*
			pow(r2, 2)*
			pow(y2, 2)) -
			2 * pow(r1, 5)*
			pow(r2, 3)*pow(y2, 2)
			- pow(r1, 4)*
			pow(r2, 4)*pow(y2, 2)
			+ 2 * pow(r1, 4)*
			pow(r2, 2)*
			pow(x1, 2)*pow(y2, 2)
			+ 2 * pow(r1, 3)*
			pow(r2, 3)*
			pow(x1, 2)*pow(y2, 2)
			+ pow(r1, 2)*
			pow(r2, 4)*
			pow(x1, 2)*pow(y2, 2)
			- pow(r1, 2)*
			pow(r2, 2)*
			pow(x1, 4)*pow(y2, 2)
			- 2 * pow(r1, 4)*
			pow(r2, 2)*x1*x2*
			pow(y2, 2) +
			2 * pow(r1, 2)*
			pow(r2, 2)*
			pow(x1, 3)*x2*
			pow(y2, 2) +
			pow(r1, 6)*pow(x2, 2)*
			pow(y2, 2) +
			2 * pow(r1, 5)*r2*
			pow(x2, 2)*pow(y2, 2)
			+ 2 * pow(r1, 4)*
			pow(r2, 2)*
			pow(x2, 2)*pow(y2, 2)
			- 2 * pow(r1, 4)*
			pow(x1, 2)*
			pow(x2, 2)*pow(y2, 2)
			- 2 * pow(r1, 3)*r2*
			pow(x1, 2)*
			pow(x2, 2)*pow(y2, 2)
			- 2 * pow(r1, 2)*
			pow(r2, 2)*
			pow(x1, 2)*
			pow(x2, 2)*pow(y2, 2)
			+ pow(r1, 2)*
			pow(x1, 4)*
			pow(x2, 2)*pow(y2, 2)
			+ 2 * pow(r1, 4)*x1*
			pow(x2, 3)*pow(y2, 2)
			- 2 * pow(r1, 2)*
			pow(x1, 3)*
			pow(x2, 3)*pow(y2, 2)
			- pow(r1, 4)*
			pow(x2, 4)*pow(y2, 2)
			+ pow(r1, 2)*
			pow(x1, 2)*
			pow(x2, 4)*pow(y2, 2)
			+ pow(r1, 6)*
			pow(y2, 4) +
			2 * pow(r1, 5)*r2*
			pow(y2, 4) +
			2 * pow(r1, 4)*
			pow(r2, 2)*pow(y2, 4)
			- 2 * pow(r1, 4)*
			pow(x1, 2)*pow(y2, 4)
			- 2 * pow(r1, 3)*r2*
			pow(x1, 2)*pow(y2, 4)
			- 2 * pow(r1, 2)*
			pow(r2, 2)*
			pow(x1, 2)*pow(y2, 4)
			+ pow(r1, 2)*
			pow(x1, 4)*pow(y2, 4)
			+ 2 * pow(r1, 4)*x1*x2*
			pow(y2, 4) -
			2 * pow(r1, 2)*
			pow(x1, 3)*x2*
			pow(y2, 4) -
			2 * pow(r1, 4)*
			pow(x2, 2)*pow(y2, 4)
			+ 2 * pow(r1, 2)*
			pow(x1, 2)*
			pow(x2, 2)*pow(y2, 4)
			- pow(r1, 4)*
			pow(y2, 6) +
			pow(r1, 2)*pow(x1, 2)*
			pow(y2, 6))) /
			(2.*(4 * pow(r2, 2)*
				pow(x1, 2) +
				8 * r1*r2*x1*x2 +
				4 * pow(r1, 2)*
				pow(x2, 2) +
				4 * pow(r1, 2)*
				pow(y2, 2) -
				4 * pow(x1, 2)*pow(y2, 2)
				));
	y = (-(pow(r1, 2)*r2) -
		r1*pow(r2, 2) +
		r2*pow(x1, 2) +
		r1*pow(x2, 2) +
		r1*pow(y2, 2) +
		(4 * pow(r1, 2)*pow(r2, 3)*
			pow(x1, 2)) /
			(4 * pow(r2, 2)*
				pow(x1, 2) +
				8 * r1*r2*x1*x2 +
				4 * pow(r1, 2)*
				pow(x2, 2) +
				4 * pow(r1, 2)*
				pow(y2, 2) -
				4 * pow(x1, 2)*
				pow(y2, 2)) +
				(4 * r1*pow(r2, 4)*
					pow(x1, 2)) /
					(4 * pow(r2, 2)*
						pow(x1, 2) +
						8 * r1*r2*x1*x2 +
						4 * pow(r1, 2)*
						pow(x2, 2) +
						4 * pow(r1, 2)*
						pow(y2, 2) -
						4 * pow(x1, 2)*
						pow(y2, 2)) -
						(4 * pow(r2, 3)*
							pow(x1, 4)) /
							(4 * pow(r2, 2)*
								pow(x1, 2) +
								8 * r1*r2*x1*x2 +
								4 * pow(r1, 2)*
								pow(x2, 2) +
								4 * pow(r1, 2)*
								pow(y2, 2) -
								4 * pow(x1, 2)*
								pow(y2, 2)) +
								(8 * pow(r1, 3)*pow(r2, 2)*
									x1*x2) /
									(4 * pow(r2, 2)*
										pow(x1, 2) +
										8 * r1*r2*x1*x2 +
										4 * pow(r1, 2)*
										pow(x2, 2) +
										4 * pow(r1, 2)*
										pow(y2, 2) -
										4 * pow(x1, 2)*
										pow(y2, 2)) +
										(8 * pow(r1, 2)*pow(r2, 3)*
											x1*x2) /
											(4 * pow(r2, 2)*
												pow(x1, 2) +
												8 * r1*r2*x1*x2 +
												4 * pow(r1, 2)*
												pow(x2, 2) +
												4 * pow(r1, 2)*
												pow(y2, 2) -
												4 * pow(x1, 2)*
												pow(y2, 2)) -
												(8 * r1*pow(r2, 2)*
													pow(x1, 3)*x2) /
													(4 * pow(r2, 2)*
														pow(x1, 2) +
														8 * r1*r2*x1*x2 +
														4 * pow(r1, 2)*
														pow(x2, 2) +
														4 * pow(r1, 2)*
														pow(y2, 2) -
														4 * pow(x1, 2)*
														pow(y2, 2)) +
														(4 * pow(r1, 4)*r2*
															pow(x2, 2)) /
															(4 * pow(r2, 2)*
																pow(x1, 2) +
																8 * r1*r2*x1*x2 +
																4 * pow(r1, 2)*
																pow(x2, 2) +
																4 * pow(r1, 2)*
																pow(y2, 2) -
																4 * pow(x1, 2)*
																pow(y2, 2)) +
																(4 * pow(r1, 3)*pow(r2, 2)*
																	pow(x2, 2)) /
																	(4 * pow(r2, 2)*
																		pow(x1, 2) +
																		8 * r1*r2*x1*x2 +
																		4 * pow(r1, 2)*
																		pow(x2, 2) +
																		4 * pow(r1, 2)*
																		pow(y2, 2) -
																		4 * pow(x1, 2)*
																		pow(y2, 2)) -
																		(4 * pow(r1, 2)*r2*
																			pow(x1, 2)*pow(x2, 2))
		/ (4 * pow(r2, 2)*
			pow(x1, 2) +
			8 * r1*r2*x1*x2 +
			4 * pow(r1, 2)*
			pow(x2, 2) +
			4 * pow(r1, 2)*
			pow(y2, 2) -
			4 * pow(x1, 2)*
			pow(y2, 2)) -
			(4 * r1*pow(r2, 2)*
				pow(x1, 2)*pow(x2, 2))
		/ (4 * pow(r2, 2)*
			pow(x1, 2) +
			8 * r1*r2*x1*x2 +
			4 * pow(r1, 2)*
			pow(x2, 2) +
			4 * pow(r1, 2)*
			pow(y2, 2) -
			4 * pow(x1, 2)*
			pow(y2, 2)) -
			(8 * pow(r1, 2)*r2*x1*
				pow(x2, 3)) /
				(4 * pow(r2, 2)*
					pow(x1, 2) +
					8 * r1*r2*x1*x2 +
					4 * pow(r1, 2)*
					pow(x2, 2) +
					4 * pow(r1, 2)*
					pow(y2, 2) -
					4 * pow(x1, 2)*
					pow(y2, 2)) -
					(4 * pow(r1, 3)*
						pow(x2, 4)) /
						(4 * pow(r2, 2)*
							pow(x1, 2) +
							8 * r1*r2*x1*x2 +
							4 * pow(r1, 2)*
							pow(x2, 2) +
							4 * pow(r1, 2)*
							pow(y2, 2) -
							4 * pow(x1, 2)*
							pow(y2, 2)) -
							(4 * pow(r1, 2)*r2*
								pow(x1, 2)*pow(y2, 2))
		/ (4 * pow(r2, 2)*
			pow(x1, 2) +
			8 * r1*r2*x1*x2 +
			4 * pow(r1, 2)*
			pow(x2, 2) +
			4 * pow(r1, 2)*
			pow(y2, 2) -
			4 * pow(x1, 2)*
			pow(y2, 2)) -
			(4 * r1*pow(r2, 2)*
				pow(x1, 2)*pow(y2, 2))
		/ (4 * pow(r2, 2)*
			pow(x1, 2) +
			8 * r1*r2*x1*x2 +
			4 * pow(r1, 2)*
			pow(x2, 2) +
			4 * pow(r1, 2)*
			pow(y2, 2) -
			4 * pow(x1, 2)*
			pow(y2, 2)) +
			(4 * r2*pow(x1, 4)*
				pow(y2, 2)) /
				(4 * pow(r2, 2)*
					pow(x1, 2) +
					8 * r1*r2*x1*x2 +
					4 * pow(r1, 2)*
					pow(x2, 2) +
					4 * pow(r1, 2)*
					pow(y2, 2) -
					4 * pow(x1, 2)*
					pow(y2, 2)) -
					(4 * pow(r1, 3)*x1*x2*
						pow(y2, 2)) /
						(4 * pow(r2, 2)*
							pow(x1, 2) +
							8 * r1*r2*x1*x2 +
							4 * pow(r1, 2)*
							pow(x2, 2) +
							4 * pow(r1, 2)*
							pow(y2, 2) -
							4 * pow(x1, 2)*
							pow(y2, 2)) -
							(8 * pow(r1, 2)*r2*x1*x2*
								pow(y2, 2)) /
								(4 * pow(r2, 2)*
									pow(x1, 2) +
									8 * r1*r2*x1*x2 +
									4 * pow(r1, 2)*
									pow(x2, 2) +
									4 * pow(r1, 2)*
									pow(y2, 2) -
									4 * pow(x1, 2)*
									pow(y2, 2)) +
									(4 * r1*pow(x1, 3)*x2*
										pow(y2, 2)) /
										(4 * pow(r2, 2)*
											pow(x1, 2) +
											8 * r1*r2*x1*x2 +
											4 * pow(r1, 2)*
											pow(x2, 2) +
											4 * pow(r1, 2)*
											pow(y2, 2) -
											4 * pow(x1, 2)*
											pow(y2, 2)) -
											(4 * pow(r1, 3)*pow(x2, 2)*
												pow(y2, 2)) /
												(4 * pow(r2, 2)*
													pow(x1, 2) +
													8 * r1*r2*x1*x2 +
													4 * pow(r1, 2)*
													pow(x2, 2) +
													4 * pow(r1, 2)*
													pow(y2, 2) -
													4 * pow(x1, 2)*
													pow(y2, 2)) +
													(4 * r2*x1*
														sqrt(-(pow(r1, 6)*
															pow(r2, 2)*
															pow(y2, 2)) -
															2 * pow(r1, 5)*
															pow(r2, 3)*
															pow(y2, 2) -
															pow(r1, 4)*
															pow(r2, 4)*
															pow(y2, 2) +
															2 * pow(r1, 4)*
															pow(r2, 2)*
															pow(x1, 2)*
															pow(y2, 2) +
															2 * pow(r1, 3)*
															pow(r2, 3)*
															pow(x1, 2)*
															pow(y2, 2) +
															pow(r1, 2)*
															pow(r2, 4)*
															pow(x1, 2)*
															pow(y2, 2) -
															pow(r1, 2)*
															pow(r2, 2)*
															pow(x1, 4)*
															pow(y2, 2) -
															2 * pow(r1, 4)*
															pow(r2, 2)*x1*x2*
															pow(y2, 2) +
															2 * pow(r1, 2)*
															pow(r2, 2)*
															pow(x1, 3)*x2*
															pow(y2, 2) +
															pow(r1, 6)*
															pow(x2, 2)*
															pow(y2, 2) +
															2 * pow(r1, 5)*r2*
															pow(x2, 2)*
															pow(y2, 2) +
															2 * pow(r1, 4)*
															pow(r2, 2)*
															pow(x2, 2)*
															pow(y2, 2) -
															2 * pow(r1, 4)*
															pow(x1, 2)*
															pow(x2, 2)*
															pow(y2, 2) -
															2 * pow(r1, 3)*r2*
															pow(x1, 2)*
															pow(x2, 2)*
															pow(y2, 2) -
															2 * pow(r1, 2)*
															pow(r2, 2)*
															pow(x1, 2)*
															pow(x2, 2)*
															pow(y2, 2) +
															pow(r1, 2)*
															pow(x1, 4)*
															pow(x2, 2)*
															pow(y2, 2) +
															2 * pow(r1, 4)*x1*
															pow(x2, 3)*
															pow(y2, 2) -
															2 * pow(r1, 2)*
															pow(x1, 3)*
															pow(x2, 3)*
															pow(y2, 2) -
															pow(r1, 4)*
															pow(x2, 4)*
															pow(y2, 2) +
															pow(r1, 2)*
															pow(x1, 2)*
															pow(x2, 4)*
															pow(y2, 2) +
															pow(r1, 6)*
															pow(y2, 4) +
															2 * pow(r1, 5)*r2*
															pow(y2, 4) +
															2 * pow(r1, 4)*
															pow(r2, 2)*
															pow(y2, 4) -
															2 * pow(r1, 4)*
															pow(x1, 2)*
															pow(y2, 4) -
															2 * pow(r1, 3)*r2*
															pow(x1, 2)*
															pow(y2, 4) -
															2 * pow(r1, 2)*
															pow(r2, 2)*
															pow(x1, 2)*
															pow(y2, 4) +
															pow(r1, 2)*
															pow(x1, 4)*
															pow(y2, 4) +
															2 * pow(r1, 4)*x1*x2*
															pow(y2, 4) -
															2 * pow(r1, 2)*
															pow(x1, 3)*x2*
															pow(y2, 4) -
															2 * pow(r1, 4)*
															pow(x2, 2)*
															pow(y2, 4) +
															2 * pow(r1, 2)*
															pow(x1, 2)*
															pow(x2, 2)*
															pow(y2, 4) -
															pow(r1, 4)*
															pow(y2, 6) +
															pow(r1, 2)*
															pow(x1, 2)*
															pow(y2, 6))) /
															(4 * pow(r2, 2)*
																pow(x1, 2) +
																8 * r1*r2*x1*x2 +
																4 * pow(r1, 2)*
																pow(x2, 2) +
																4 * pow(r1, 2)*
																pow(y2, 2) -
																4 * pow(x1, 2)*
																pow(y2, 2)) +
																(4 * r1*x2*
																	sqrt(-(pow(r1, 6)*
																		pow(r2, 2)*
																		pow(y2, 2)) -
																		2 * pow(r1, 5)*
																		pow(r2, 3)*
																		pow(y2, 2) -
																		pow(r1, 4)*
																		pow(r2, 4)*
																		pow(y2, 2) +
																		2 * pow(r1, 4)*
																		pow(r2, 2)*
																		pow(x1, 2)*
																		pow(y2, 2) +
																		2 * pow(r1, 3)*
																		pow(r2, 3)*
																		pow(x1, 2)*
																		pow(y2, 2) +
																		pow(r1, 2)*
																		pow(r2, 4)*
																		pow(x1, 2)*
																		pow(y2, 2) -
																		pow(r1, 2)*
																		pow(r2, 2)*
																		pow(x1, 4)*
																		pow(y2, 2) -
																		2 * pow(r1, 4)*
																		pow(r2, 2)*x1*x2*
																		pow(y2, 2) +
																		2 * pow(r1, 2)*
																		pow(r2, 2)*
																		pow(x1, 3)*x2*
																		pow(y2, 2) +
																		pow(r1, 6)*
																		pow(x2, 2)*
																		pow(y2, 2) +
																		2 * pow(r1, 5)*r2*
																		pow(x2, 2)*
																		pow(y2, 2) +
																		2 * pow(r1, 4)*
																		pow(r2, 2)*
																		pow(x2, 2)*
																		pow(y2, 2) -
																		2 * pow(r1, 4)*
																		pow(x1, 2)*
																		pow(x2, 2)*
																		pow(y2, 2) -
																		2 * pow(r1, 3)*r2*
																		pow(x1, 2)*
																		pow(x2, 2)*
																		pow(y2, 2) -
																		2 * pow(r1, 2)*
																		pow(r2, 2)*
																		pow(x1, 2)*
																		pow(x2, 2)*
																		pow(y2, 2) +
																		pow(r1, 2)*
																		pow(x1, 4)*
																		pow(x2, 2)*
																		pow(y2, 2) +
																		2 * pow(r1, 4)*x1*
																		pow(x2, 3)*
																		pow(y2, 2) -
																		2 * pow(r1, 2)*
																		pow(x1, 3)*
																		pow(x2, 3)*
																		pow(y2, 2) -
																		pow(r1, 4)*
																		pow(x2, 4)*
																		pow(y2, 2) +
																		pow(r1, 2)*
																		pow(x1, 2)*
																		pow(x2, 4)*
																		pow(y2, 2) +
																		pow(r1, 6)*
																		pow(y2, 4) +
																		2 * pow(r1, 5)*r2*
																		pow(y2, 4) +
																		2 * pow(r1, 4)*
																		pow(r2, 2)*
																		pow(y2, 4) -
																		2 * pow(r1, 4)*
																		pow(x1, 2)*
																		pow(y2, 4) -
																		2 * pow(r1, 3)*r2*
																		pow(x1, 2)*
																		pow(y2, 4) -
																		2 * pow(r1, 2)*
																		pow(r2, 2)*
																		pow(x1, 2)*
																		pow(y2, 4) +
																		pow(r1, 2)*
																		pow(x1, 4)*
																		pow(y2, 4) +
																		2 * pow(r1, 4)*x1*x2*
																		pow(y2, 4) -
																		2 * pow(r1, 2)*
																		pow(x1, 3)*x2*
																		pow(y2, 4) -
																		2 * pow(r1, 4)*
																		pow(x2, 2)*
																		pow(y2, 4) +
																		2 * pow(r1, 2)*
																		pow(x1, 2)*
																		pow(x2, 2)*
																		pow(y2, 4) -
																		pow(r1, 4)*
																		pow(y2, 6) +
																		pow(r1, 2)*
																		pow(x1, 2)*
																		pow(y2, 6))) /
																		(4 * pow(r2, 2)*
																			pow(x1, 2) +
																			8 * r1*r2*x1*x2 +
																			4 * pow(r1, 2)*
																			pow(x2, 2) +
																			4 * pow(r1, 2)*
																			pow(y2, 2) -
																			4 * pow(x1, 2)*
																			pow(y2, 2))) /
																			(2.*r1*y2);
	r = sqrt(x*x + y*y);

	return rg_Circle2D(x, y, r);
}


















void rg_Circle2D::shrinkCircle( const         rg_Circle2D& c1,
                                const         rg_Circle2D& c2,
                                const         rg_Circle2D& c3,
                                rg_Circle2D&  c1Tilde,
                                rg_Circle2D&  c2Tilde,
                                rg_Point2D&   smallest)
{
	if( c1.getRadius() < c2.getRadius() )
	{
		if( c1.getRadius() < c3.getRadius() )
		{
			c1Tilde.setCircle(c2.getCenterPt(), c2.getRadius() - c1.getRadius());
			c2Tilde.setCircle(c3.getCenterPt(), c3.getRadius() - c1.getRadius());
			smallest = c1.getCenterPt();
		}
		else
		{
			c1Tilde.setCircle(c1.getCenterPt(), c1.getRadius() - c3.getRadius());
			c2Tilde.setCircle(c2.getCenterPt(), c2.getRadius() - c3.getRadius());
			smallest = c3.getCenterPt();
		}
	}
	else
	{
		if( c2.getRadius() < c3.getRadius() )
		{
			c1Tilde.setCircle(c1.getCenterPt(), c1.getRadius() - c2.getRadius());
			c2Tilde.setCircle(c3.getCenterPt(), c3.getRadius() - c2.getRadius());
			smallest = c2.getCenterPt();
		}
		else
		{
			c1Tilde.setCircle(c1.getCenterPt(), c1.getRadius() - c3.getRadius());
			c2Tilde.setCircle(c2.getCenterPt(), c2.getRadius() - c3.getRadius());
			smallest = c3.getCenterPt();
		}
	}
}

//this function returns tangent lines 
void rg_Circle2D::makeExteriorTangentLinesOf2Circles(const rg_Circle2D&	circle1,
											         const rg_Circle2D&	circle2,
											         rg_ImplicitEquation&	result1,
											         rg_ImplicitEquation&	result2)
{
	rg_Circle2D w1, w2;

	//the radius of w2 is less than that of w1
	if( circle1.getRadius() < circle2.getRadius() )
	{
		w1 = circle2;
		w2 = circle1;
	}
	else
	{
		w1 = circle1;
		w2 = circle2;
	}

	rg_Point2D c2cVector = w1.getCenterPt() - w2.getCenterPt();
	rg_REAL r = w1.getRadius() - w2.getRadius();
	rg_REAL length = c2cVector.magnitude();
	rg_REAL sine = r / length;
	rg_REAL cosine = sqrt( length*length - r*r ) / length;

	//rotate theta  /  -theta
	rg_Point2D normal1( c2cVector.getX() * cosine - c2cVector.getY() * sine ,
					 c2cVector.getX() * sine + c2cVector.getY() * cosine );
	rg_Point2D normal2( c2cVector.getX() * cosine + c2cVector.getY() * sine ,
					 c2cVector.getY() * cosine - c2cVector.getX() * sine );
	normal1 = normal1.getUnitVector();
	normal2 = normal2.getUnitVector();
	
	//rotate -PI/2  /  PI/2
	normal1.setPoint( normal1.getY() , -1 * normal1.getX() );
	normal2.setPoint( -1 * normal2.getY() , normal2.getX() );

	result1.setCoeff(1, 0, normal1.getX());
	result1.setCoeff(0, 1, normal1.getY());
	result1.setCoeff(0, 0, -1 * normal1.getX() * w2.getCenterPt().getX() 
						   - normal1.getY() * w2.getCenterPt().getY() 
						   + w2.getRadius()   );

	result2.setCoeff(1, 0, normal2.getX());
	result2.setCoeff(0, 1, normal2.getY());
	result2.setCoeff(0, 0, -1 * normal2.getX() * w2.getCenterPt().getX() 
						   - normal2.getY() * w2.getCenterPt().getY() 
						   + w2.getRadius()   );

}


void rg_Circle2D::makeInteriorTangentLinesOf2Circles(const rg_Circle2D& circle1,
    const rg_Circle2D& circle2,
    rg_ImplicitEquation& result1,
    rg_ImplicitEquation& result2)
{
    if (circle1.isIntersectWith(circle2))
    {
        result1.setDegree(0);
        result2.setDegree(0);

        return;
    }

    rg_Circle2D w1, w2;

    //the radius of w2 is less than that of w1
    if (circle1.getRadius() < circle2.getRadius())
    {
        w1 = circle2;
        w2 = circle1;
    }
    else
    {
        w1 = circle1;
        w2 = circle2;
    }

    rg_Point2D c2cVector = w1.getCenterPt() - w2.getCenterPt();
    rg_REAL r = w1.getRadius() + w2.getRadius();
    rg_REAL length = c2cVector.magnitude();
    rg_REAL sine = r / length;
    rg_REAL cosine = sqrt(length*length - r*r) / length;

    //rotate theta  /  -theta
    rg_Point2D normal1(c2cVector.getX() * cosine - c2cVector.getY() * sine,
                       c2cVector.getX() * sine + c2cVector.getY() * cosine);
    rg_Point2D normal2(c2cVector.getX() * cosine + c2cVector.getY() * sine,
                       c2cVector.getY() * cosine - c2cVector.getX() * sine);
    normal1 = normal1.getUnitVector();
    normal2 = normal2.getUnitVector();

    //rotate -PI/2  /  PI/2
    normal1.setPoint(normal1.getY(), -1 * normal1.getX());
    normal2.setPoint(-1 * normal2.getY(), normal2.getX());

    result1.setCoeff(1, 0, normal1.getX());
    result1.setCoeff(0, 1, normal1.getY());
    result1.setCoeff(0, 0, -1 * normal1.getX() * w2.getCenterPt().getX()
                           - normal1.getY() * w2.getCenterPt().getY()
                           - w2.getRadius());

    result2.setCoeff(1, 0, normal2.getX());
    result2.setCoeff(0, 1, normal2.getY());
    result2.setCoeff(0, 0, -1 * normal2.getX() * w2.getCenterPt().getX()
                           - normal2.getY() * w2.getCenterPt().getY()
                           - w2.getRadius());
}



rg_Point2D rg_Circle2D::makeExteriorTangentLinesInWPlane( const rg_Circle2D& c1, 
//rg_Point2D makeExteriorTangentLinesInWPlane( const rg_Circle2D& c1, 
										                  const rg_Circle2D& c2,
											              const rg_Circle2D& c3,
											              rg_ImplicitEquation& result1,
											              rg_ImplicitEquation& result2)
{
	rg_Circle2D c1Tilde, c2Tilde;
	rg_Point2D z3;

	//c3 has not the smallest radius between three circles
	//so we must make c3 to have smallest radius by swapping for convienience
	shrinkCircle(c1, c2, c3, c1Tilde, c2Tilde, z3);

	rg_Circle2D w1, w2;
	w1 = transformZ2W( c1Tilde, z3 );
	w2 = transformZ2W( c2Tilde, z3 );

	//need to be added generating two tangent lines
	makeExteriorTangentLinesOf2Circles(w1, w2, result1, result2);

	return z3;
}

//
// W = 1 / (z - z3)
//
rg_Circle2D rg_Circle2D::transformZ2W(const rg_Circle2D&	cTilde,
								   const rg_Point2D&	smallest)
{
	rg_REAL x1 = cTilde.getCenterPt().getX();
	rg_REAL y1 = cTilde.getCenterPt().getY();
	rg_REAL x3 = smallest.getX();
	rg_REAL y3 = smallest.getY();

	rg_REAL p1 = pow(x1 - x3, 2) + pow(y1 - y3, 2) - pow(cTilde.getRadius(), 2);

	return rg_Circle2D( (x1 - x3)/p1 , -(y1 - y3)/p1, cTilde.getRadius() / p1 );
}

//
// Z = 1 / w + z3
//
rg_Circle2D rg_Circle2D::transformW2Z(const rg_ImplicitEquation& line, const rg_Point2D& z3)
{
	//degenerate case must be considered
	//case c == 0;

	rg_REAL a = line.getCoeff(1, 0);
	rg_REAL b = line.getCoeff(0, 1);
	rg_REAL c = line.getCoeff(0, 0);

	rg_REAL tempX = -a/2./c + z3.getX();
	rg_REAL tempY =  b/2./c + z3.getY();
	rg_REAL rad = 1./2./c;

	return rg_Circle2D(tempX, tempY, rad);
}


rg_BOOL rg_Circle2D::isEqual( const rg_Circle2D& circle ) const
{
    if( rg_EQ( centerPt.getX(), circle.centerPt.getX() ) &&
        rg_EQ( centerPt.getY(), circle.centerPt.getY() ) &&
        rg_EQ( radius, circle.radius ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}

#ifdef PYVORONOI
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

using Circle2D = rg_Circle2D;
using Point2D = rg_Point2D;

void init_Circle2D(py::module& m) {
	py::class_<Circle2D, std::shared_ptr<Circle2D>>(m, "Circle2D")
		.def(py::init([]() { return new Circle2D(); }))
		.def(py::init([](const double& x, const double& y, const double& radius) { return new Circle2D(x, y, radius); }))
		.def(py::init([](const Point2D& center, const double& radius) { return new Circle2D(center, radius); }))
		.def(py::init([](const Circle2D& circle) { return new Circle2D(circle); }))
		.def_property("center", &Circle2D::getCenterPt, &Circle2D::setCenterPt)
		.def_property("x", &Circle2D::getX, &Circle2D::setX)
		.def_property("y", &Circle2D::getY, &Circle2D::setY)
		.def_property("radius", &Circle2D::getRadius, &Circle2D::setRadius)
		.def("get_center", &Circle2D::getCenterPt)
		.def("get_x", &Circle2D::getX)
		.def("get_y", &Circle2D::getY)
		.def("get_radius", &Circle2D::getRadius)
		.def("set_center", &Circle2D::setCenterPt)
		.def("set_x", &Circle2D::setX)
		.def("set_y", &Circle2D::setY)
		.def("set_radius", &Circle2D::setRadius)
		.def("set_circle", static_cast< void (Circle2D::*)(const Circle2D&) > 				(&Circle2D::setCircle))
		.def("set_circle", static_cast< void (Circle2D::*)(const Point2D&,  const rg_REAL&) > (&Circle2D::setCircle))
		.def("distance", static_cast< double (Circle2D::*)(const Circle2D&) const> 			(&Circle2D::distance))   
		.def("distance", static_cast< double (Circle2D::*)(const Point2D&)  const>  			(&Circle2D::distance))   
		.def(py::self == py::self)
		// .def(py::self != py::self)
		;
}
#endif