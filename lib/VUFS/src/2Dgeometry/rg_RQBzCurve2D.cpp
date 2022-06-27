/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_RQBzCurve2D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_RQBzCurve2D 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jun 1997    
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include "rg_RQBzCurve2D.h"
#include "rg_Point2D.h"
#include "rg_Line.h"
#include "rg_IntersectFunc.h"
#include "rg_MathFunc.h"
#include "rg_Matrix.h"

rg_RQBzCurve2D::rg_RQBzCurve2D() :rg_RBzCurve2D(rg_QUADRATIC)
{
}

rg_RQBzCurve2D::rg_RQBzCurve2D(const rg_RQBzCurve2D &curve) :rg_RBzCurve2D(curve)
{
}

rg_RQBzCurve2D::rg_RQBzCurve2D(const rg_Point2D *ctrlpt) :rg_RBzCurve2D(rg_QUADRATIC, ctrlpt)
{
}

rg_RQBzCurve2D::~rg_RQBzCurve2D()
{
}

/*
void rg_RQBzCurve2D::makeRQBezier(rg_CBzCurve2D &curve)
{
	rg_Point2D *pt = curve.getCtrlPoint();

	makeRQBezier(pt[0], 
		         rg_IntersectFunc::intersectLineVsLine(rg_Line<rg_Point2D>(pt[0], pt[1]), rg_Line<rg_Point2D>(pt[3], pt[2])), 
 		         pt[3], 
				 curve.evaluatePt(0.5));

	delete []pt;
}
*/
void rg_RQBzCurve2D::makeRQBezier(const rg_Point2D& b0, 
                             const rg_Point2D& t0, 
                             const rg_Point2D& b2, 
							 const rg_Point2D& t2,
                             const rg_Point2D& p) 
{
    rg_Point2D b1;
    rg_FLAG isLinear= rg_EQ( t0*t2, 0.0 );
    if ( isLinear ==rg_TRUE )
    {
        b1= (b0+b1)*0.5;
        rg_BzCurve2D::setCtrlPt(0,b0);
        rg_BzCurve2D::setCtrlPt(1,b1);
        rg_BzCurve2D::setCtrlPt(2,b2);
        
        rg_RBzCurve2D::setWeight(0,1.0);
        rg_RBzCurve2D::setWeight(1,1.0);
        rg_RBzCurve2D::setWeight(2,1.0);
    }
    else
    {
        b1=rg_IntersectFunc::intersectLineVsLine(rg_Line<rg_Point2D>(b0, b0+t0), rg_Line<rg_Point2D>(b2, b2+t2));
    	makeRQBezier(b0, 
                     b1, 
 		    		 b2, 
			    	 p);
    }
}

void rg_RQBzCurve2D::makeRQBezier(const rg_Point2D& b0, // first control point
                             const rg_Point2D& b1, // second control point
                             const rg_Point2D& b2, // third control point
                             const rg_Point2D& p)  // passing point on the quadratic curve
{
	rg_REAL b0x = b0.getX();
    rg_REAL b0y = b0.getY();
    rg_REAL b1x = b1.getX();
    rg_REAL b1y = b1.getY();
    rg_REAL b2x = b2.getX();
    rg_REAL b2y = b2.getY();
    rg_REAL px  =  p.getX();
    rg_REAL py  =  p.getY();

	rg_Matrix D(3, 3);
	D[0][0] = b0x;
	D[0][1] = b1x;
	D[0][2] = b2x;
	D[1][0] = b0y;
	D[1][1] = b1y;
	D[1][2] = b2y;
	D[2][0] = 1.0;
	D[2][1] = 1.0;
	D[2][2] = 1.0;

	rg_Matrix tau0 = D;
	tau0[0][0] = px;
	tau0[1][0] = py;

	rg_Matrix tau1 = D;
	tau1[0][1] = px;
	tau1[1][1] = py;
	
	rg_Matrix tau2 = D;
	tau2[0][2] = px;
	tau2[1][2] = py;
	
	rg_REAL  td  = D.determinant();
    rg_REAL  tw0 = tau0.determinant()/ td;
    rg_REAL  tw1 = tau1.determinant()/ td;
    rg_REAL  tw2 = tau2.determinant()/ td;
    rg_REAL  w1  = tw1/(2 * sqrt( abs(tw0 * tw2)) );

	rg_Point2D *ctrlpt = new rg_Point2D[getDegree()+1];
	ctrlpt[0] = rg_Point2D(b0x, b0y);
	ctrlpt[1] = rg_Point2D(b1x, b1y);
	ctrlpt[2] = rg_Point2D(b2x, b2y);
	setCtrlPts(ctrlpt);
    delete[] ctrlpt;

	weight[0] = 1.0;
	weight[1] =  w1;
	weight[2] = 1.0;
}


rg_REAL rg_RQBzCurve2D::getParameter4Point(const rg_Point2D &point) const
{
	// Supposing a given point, P(X, Y), 
	// and a rational quardratic Bezier curve, 
	// c(t) = (x(t)/w(t), y(t)/w(t)),
	// where x(t) = a2 t^2 + a1 t + a0,
	// y(t) = b2 t^2 + b1 t + b0,
	// and w(t) = d2 t^2 + d1 t + d0.

	// Since a given point lies on c(t)
	// the following equations are satisfied.
	// x(t) = X * w(t)			Eq(1)
	// and y(t) = Y * w(t).		Eq(2)

	// There, we find t to satisfy Eq(1) and Eq(2).

	rg_REAL x0 = ctrlPts[0].getX();
	rg_REAL x1 = ctrlPts[1].getX();
	rg_REAL x2 = ctrlPts[2].getX();
	rg_REAL y0 = ctrlPts[0].getY();
	rg_REAL y1 = ctrlPts[1].getY();
	rg_REAL y2 = ctrlPts[2].getY();
	rg_REAL w0 = weight[0];
	rg_REAL w1 = weight[1];
	rg_REAL w2 = weight[2];

    rg_REAL a2 =    w0*x0 - 2*w1*x1 + w2*x2;
	rg_REAL a1 = -2*w0*x0 + 2*w1*x1;
	rg_REAL a0 =    w0*x0;

	rg_REAL b2 =    w0*y0 - 2*w1*y1 + w2*y2;
	rg_REAL b1 = -2*w0*y0 + 2*w1*y1;
	rg_REAL b0 =    w0*y0;

	rg_REAL d2 =    w0 - 2*w1 + w2;
	rg_REAL d1 = -2*w0 + 2*w1;
	rg_REAL d0 =    w0;

	rg_REAL X = point.getX();
	rg_REAL Y = point.getY();

	rg_REAL *tx = rg_MathFunc::solveQuadraticEq(a2-X*d2, a1-X*d1, a0-X*d0);
	rg_REAL *ty = rg_MathFunc::solveQuadraticEq(b2-Y*d2, b1-Y*d1, b0-Y*d0);

/*
	if( rg_BTORexclusive(0.0, tx[0], 1.0) )
	{
		if( rg_EQ(tx[0], ty[0]) || rg_EQ(tx[0], ty[1]) )
		{
			return tx[0];
		}
		else 
		{
			return 0.0;
		}
	}
	
	else if( rg_BTORexclusive(0.0, tx[1], 1.0) )
	{
		if( rg_EQ(tx[1], ty[0]) || rg_EQ(tx[1], ty[1]) )
		{
			return tx[1];
		}
		else 
		{
			return 0.0;
		}
	}

	else 
	{
		return 0.0;
	}
*/
	if( rg_EQ(tx[0], ty[0]) 
	 && rg_BTORexclusive(0.0, tx[0], 1.0) ) return tx[0];

	else if ( rg_EQ(tx[0], ty[1])
	       && rg_BTORexclusive(0.0, tx[0], 1.0) ) return tx[0]; 

	else if ( rg_EQ(tx[1], ty[0])
		   && rg_BTORexclusive(0.0, tx[1], 1.0) ) return tx[1];

	else if ( rg_EQ(tx[1], ty[1])
		   && rg_BTORexclusive(0.0, tx[1], 1.0) ) return tx[1];

	else 
	{
		return 0.0;
	}
 }

rg_BoundingBox2D rg_RQBzCurve2D::makeBoundingBox() const
{
    rg_REAL   t  = weight[1]/(weight[1]+1);
    rg_Point2D  q0 = (1-t)*ctrlPts[0] + t*ctrlPts[1];
    rg_Point2D  q1 = (1-t)*ctrlPts[2] + t*ctrlPts[1];

    rg_BoundingBox2D box;
    box.contain(ctrlPts[0]);
    box.contain(q0);
    box.contain(q1);
    box.contain(ctrlPts[2]);

    return box;
}

rg_RQBzCurve2D& rg_RQBzCurve2D::operator=(const rg_RQBzCurve2D& curve)
{
	rg_BzCurve2D::setDegree(curve.degree);
	rg_BzCurve2D::setCtrlPts(curve.ctrlPts);
	rg_RBzCurve2D::setWeight(curve.weight);

	return *this;
}


