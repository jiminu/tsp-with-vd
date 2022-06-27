#include "FunctionsForVoronoiDiagram3D.h"
#include "rg_RelativeOp.h"

#include "rg_BallGenerator.h"
#include "VDVertex.h"
#include "VDCell.h"
using namespace V::GeometryTier;

#include "rg_TMatrix3D.h"



rg_INT  V::GeometryTier::classifyNormalSign(rg_REAL xNormal, rg_REAL yNormal, rg_REAL zNormal)
{
    rg_INT normalSign = 0;
    if ( rg_GE(xNormal, 0.0) )  {
        if ( rg_GE(yNormal, 0.0) ) {
            if ( rg_GE(zNormal, 0.0) )
                normalSign = XGE_YGE_ZGE;
            else
                normalSign = XGE_YGE_ZLT;
        }
        else  {
            if ( rg_GE(zNormal, 0.0) )
                normalSign = XGE_YLT_ZGE;
            else
                normalSign = XGE_YLT_ZLT;
        }
    }
    else  {
        if ( rg_GE(yNormal, 0.0) ) {
            if ( rg_GE(zNormal, 0.0) )
                normalSign = XLT_YGE_ZGE;
            else
                normalSign = XLT_YGE_ZLT;
        }
        else  {
            if ( rg_GE(zNormal, 0.0) )
                normalSign = XLT_YLT_ZGE;
            else
                normalSign = XLT_YLT_ZLT;
        }
    }

    return normalSign;
}



int V::GeometryTier::compareInternalVoids(const void* ts1, const void* ts2)
{
    VDVertex* internalVoid1 = *(VDVertex**) ts1;
    VDVertex* internalVoid2 = *(VDVertex**) ts2;

    // for decending order
    if ( rg_GT( internalVoid1->getRadiusOfTangentSphere(), internalVoid2->getRadiusOfTangentSphere() ) )
        return -1;
    else if ( rg_LT( internalVoid1->getRadiusOfTangentSphere(), internalVoid2->getRadiusOfTangentSphere() ) )
        return 1;
    else 
        return 0;
}



int V::GeometryTier::compareVDCellByID(const void* cell1, const void* cell2)
{
    VDVertex* vdCell1 = *(VDVertex**) cell1;
    VDVertex* vdCell2 = *(VDVertex**) cell2;

    //  for ascending order
    if ( vdCell1->getID() > vdCell2->getID() )
        return 1;
    else if ( vdCell1->getID() < vdCell2->getID() )
        return -1;
    else 
        return 0;
}



int V::GeometryTier::compareBallsByRadius(const void* ball1, const void* ball2)
{
    BallGenerator* ballGen1 = *(BallGenerator**) ball1;
    BallGenerator* ballGen2 = *(BallGenerator**) ball2;

    //  for ascending order
    if ( ballGen1->getRadius() > ballGen2->getRadius() )
        return 1;
    else if ( ballGen1->getRadius() < ballGen2->getRadius() )
        return -1;
    else 
        return 0;
}



int V::GeometryTier::compareBallsByRadiusInDescendingOrder(const void* ball1, const void* ball2)
{
    BallGenerator* ballGen1 = *(BallGenerator**) ball1;
    BallGenerator* ballGen2 = *(BallGenerator**) ball2;

    //  for descending order
    if ( ballGen1->getRadius() < ballGen2->getRadius() )
        return 1;
    else if ( ballGen1->getRadius() > ballGen2->getRadius() )
        return -1;
    else 
        return 0;
}



rg_INT V::GeometryTier::computeCircleTangentTo3CirclesOutside_GLOBAL( const rg_Circle2D& circle1, 
                                              const rg_Circle2D& circle2, 
                                              const rg_Circle2D& circle3, 
                                              rg_Circle2D* tangentCircle)
{
    rg_Circle2D circle[3];
    if (    rg_LT(circle1.getRadius(), circle2.getRadius()) 
         && rg_LT(circle1.getRadius(), circle3.getRadius()) )  {
        circle[0] = circle2;
        circle[1] = circle3;
        circle[2] = circle1;
    }
    else if (    rg_LT(circle2.getRadius(), circle1.getRadius()) 
              && rg_LT(circle2.getRadius(), circle3.getRadius()) )  {
        circle[0] = circle1;
        circle[1] = circle3;
        circle[2] = circle2;
    }
    else  {
        circle[0] = circle1;
        circle[1] = circle2;
        circle[2] = circle3;
    }

	double x1 = circle[0].getCenterPt().getX() - circle[2].getCenterPt().getX();
	double y1 = circle[0].getCenterPt().getY() - circle[2].getCenterPt().getY();
	double r1 = circle[0].getRadius()          - circle[2].getRadius();

	double x2 = circle[1].getCenterPt().getX() - circle[2].getCenterPt().getX();
	double y2 = circle[1].getCenterPt().getY() - circle[2].getCenterPt().getY();
	double r2 = circle[1].getRadius()          - circle[2].getRadius();
    

    double determinant1 = x1*y2 - x2*y1;
	double determinant2 = x1*r2 - x2*r1;
	double determinant3 = y1*r2 - y2*r1;

    //system of equations
	//eq1: x1*x + y1*y + r1*r = rhs1
	//eq2: x2*x + y2*y + r2*r = rhs2
	//eq4: x*x + y*y + z*z = r*r

	double rhs1 = (x1*x1+y1*y1-r1*r1)*0.5;
	double rhs2 = (x2*x2+y2*y2-r2*r2)*0.5;

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

		double a_x = (y2*rhs1 - y1*rhs2)/determinant1;
        double b_x = (y1*r2   - y2*r1)  /determinant1;
		double a_y = (x1*rhs2 - x2*rhs1)/determinant1;
        double b_y = (x2*r1   - x1*r2)  /determinant1;

		double a_r = b_x*b_x + b_y*b_y - 1.0;
        double b_r = 2.0*(a_x*b_x + a_y*b_y);
        double c_r = a_x*a_x + a_y*a_y;

		double det = b_r*b_r - 4.*a_r*c_r; //b^2 - 4ac
		if(det < 0) 
			return 0;

		double radius1 = (-b_r + sqrt(det))/(2.* a_r);
		double radius2 = (-b_r - sqrt(det))/(2.* a_r);

        double x;
        double y;
        if( radius1 > 0 && radius2 > 0 )  {

			x = a_x + b_x*radius1 + circle[2].getCenterPt().getX();
			y = a_y + b_y*radius1 + circle[2].getCenterPt().getY();
			tangentCircle[0].setCircle( rg_Point2D(x, y), radius1-circle[2].getRadius() );

			x = a_x + b_x*radius2 + circle[2].getCenterPt().getX();
			y = a_y + b_y*radius2 + circle[2].getCenterPt().getY();
			tangentCircle[1].setCircle( rg_Point2D(x, y), radius2-circle[2].getRadius() );

			return 2;
		}
		else if( radius1 > 0 )  {

			x = a_x + b_x*radius1 + circle[2].getCenterPt().getX();
			y = a_y + b_y*radius1 + circle[2].getCenterPt().getY();
			tangentCircle[0].setCircle( rg_Point2D(x, y), radius1-circle[2].getRadius() );

			return 1;
		}
		else if( radius2 > 0 )  {

			x = a_x + b_x*radius2 + circle[2].getCenterPt().getX();
			y = a_y + b_y*radius2 + circle[2].getCenterPt().getY();
			tangentCircle[0].setCircle( rg_Point2D(x, y), radius2-circle[2].getRadius() );

			return 1;
		}
		else  {
			return 0;
		}    
    }
    else if( !rg_ZERO(determinant2) )  {
        //system of equations
	    //eq1: x1*x + r1*r = rhs1 - y1*y 
	    //eq2: x2*x + r2*r = rhs2 - y2*y 
        //
        //  x = a_x + b_x*y
        //  r = a_r + b_r*y
        //
        //  eq3 <- x, r
        //  a_y*y^2 + b_y*y + c_y = 0

		double a_x = (r2*rhs1 - r1*rhs2)/determinant2;
        double b_x = (r1*y2   - r2*y1)  /determinant2;
		double a_r = (x1*rhs2 - x2*rhs1)/determinant2;
        double b_r = (x2*y1   - x1*y2)  /determinant2;

		double a_y = b_x*b_x + b_r*b_r - 1.0;
        double b_y = 2.0*(a_x*b_x + a_r*b_r);
        double c_y = a_x*a_x + a_r*a_r;

		double det = b_y*b_y - 4.*a_y*c_y; //b^2 - 4ac
		if(det < 0) 
			return 0;

		double y1 = (-b_y + sqrt(det))/(2.* a_y);
		double y2 = (-b_y - sqrt(det))/(2.* a_y);
        
		double radius1 = a_r + b_r*y1;
		double radius2 = a_r + b_r*y2;
        
        double x;
        double y;
        if( radius1 > 0 && radius2 > 0 )  {

			x = a_x + b_x*radius1 + circle[2].getCenterPt().getX();
			y = y1                + circle[2].getCenterPt().getY();
			tangentCircle[0].setCircle( rg_Point2D(x, y), radius1-circle[2].getRadius() );

			x = a_x + b_x*radius2 + circle[2].getCenterPt().getX();
			y = y2                + circle[2].getCenterPt().getY();
			tangentCircle[1].setCircle( rg_Point2D(x, y), radius2-circle[2].getRadius() );

			return 2;
		}
		else if( radius1 > 0 )  {

			x = a_x + b_x*radius1 + circle[2].getCenterPt().getX();
			y = y1                + circle[2].getCenterPt().getY();
			tangentCircle[0].setCircle( rg_Point2D(x, y), radius1-circle[2].getRadius() );

			return 1;
		}
		else if( radius2 > 0 )  {

			x = a_x + b_x*radius2 + circle[2].getCenterPt().getX();
			y = y2                + circle[2].getCenterPt().getY();
			tangentCircle[0].setCircle( rg_Point2D(x, y), radius2-circle[2].getRadius() );

			return 1;
		}
		else  {
			return 0;
		}        
    }
    else if( !rg_ZERO(determinant3) )  {
        //system of equations
	    //eq1: y1*y + r1*r = rhs1 - x1*x 
	    //eq2: y2*y + r2*r = rhs2 - x2*x 
        //
        //  y = a_y + b_y*x
        //  r = a_r + b_r*x
        //
        //  eq3 <- y, r
        //  a_x*x^2 + b_x*x + c_x = 0

		double a_y = (r2*rhs1 - r1*rhs2)/determinant2;
        double b_y = (r1*x2   - r2*x1)  /determinant2;
		double a_r = (y1*rhs2 - y2*rhs1)/determinant2;
        double b_r = (y2*x1   - y1*x2)  /determinant2;

		double a_x = b_y*b_y + b_r*b_r - 1.0;
        double b_x = 2.0*(a_y*b_y + a_r*b_r);
        double c_x = a_y*a_y + a_r*a_r;

		double det = b_x*b_x - 4.*a_x*c_x; //b^2 - 4ac
		if(det < 0) 
			return 0;

		double x1 = (-b_x + sqrt(det))/(2.* a_x);
		double x2 = (-b_x - sqrt(det))/(2.* a_x);
        
		double radius1 = a_r + b_r*x1;
		double radius2 = a_r + b_r*x2;
        
        double x;
        double y;
        if( radius1 > 0 && radius2 > 0 )  {

			x = x1                + circle[2].getCenterPt().getX();
			y = a_y + b_y*radius1 + circle[2].getCenterPt().getY();
			tangentCircle[0].setCircle( rg_Point2D(x, y), radius1-circle[2].getRadius() );

			x = x2                + circle[2].getCenterPt().getX();
			y = a_y + b_y*radius2 + circle[2].getCenterPt().getY();
			tangentCircle[1].setCircle( rg_Point2D(x, y), radius2-circle[2].getRadius() );

			return 2;
		}
		else if( radius1 > 0 )  {

			x = x1                + circle[2].getCenterPt().getX();
			y = a_y + b_y*radius1 + circle[2].getCenterPt().getY();
			tangentCircle[0].setCircle( rg_Point2D(x, y), radius1-circle[2].getRadius() );

			return 1;
		}
		else if( radius2 > 0 )  {

			x = x2                + circle[2].getCenterPt().getX();
			y = a_y + b_y*radius2 + circle[2].getCenterPt().getY();
			tangentCircle[0].setCircle( rg_Point2D(x, y), radius2-circle[2].getRadius() );

			return 1;
		}
		else  {
			return 0;
		}        
    }
    else  {
        return 0;
    }

}



rg_INT V::GeometryTier::computeMinTangentSpheresTo3SpheresOutside_GLOBAL( const Sphere& sphere1, 
                                                         const Sphere& sphere2, 
                                                         const Sphere& sphere3, 
                                                         Sphere* minTangentSphere)
{
    //Sphere gateBall[3] = { m_gateBallCCW1->getBall(), m_gateBallCCW2->getBall(), m_gateBallCCW3->getBall() }; 
    rg_Point3D centerOfGB[3] = { sphere1.getCenter(), sphere2.getCenter(), sphere3.getCenter() };
    rg_REAL    radiusOfGB[3] = { sphere1.getRadius(), sphere2.getRadius(), sphere3.getRadius() };



	rg_TMatrix3D trMatrix;
	rg_TMatrix3D invMatrix;

    rg_REAL length[3] = { centerOfGB[0].distance( centerOfGB[1] ),
                          centerOfGB[0].distance( centerOfGB[2] ),
                          centerOfGB[1].distance( centerOfGB[2] ) };

    //  normal of center plane for a linear, parabolic, hyperbolic, or elliptic edge.
    if ( rg_GT(length[0]+length[1], length[2]) && rg_GT(length[0]+length[2], length[1]) && rg_GT(length[1]+length[2], length[0]) )  {
	    rg_Point3D vec21 = centerOfGB[1] - centerOfGB[0];
	    rg_Point3D vec31 = centerOfGB[2] - centerOfGB[0];

        rg_Point3D normalOfCenterPlane;
	    normalOfCenterPlane = vec21.crossProduct(vec31);
        normalOfCenterPlane.normalize();


	    trMatrix.translate(-centerOfGB[0]);

        rg_FLAG bReverse = rg_EQ( normalOfCenterPlane.innerProduct( rg_Point3D(0.0, 0.0, 1.0) ), -1);
        if( !bReverse )
	        trMatrix.rotate(normalOfCenterPlane, rg_Point3D(0.0, 0.0, 1.0));
        else
            trMatrix.rotateY(rg_PI);

        if( !bReverse )
	        invMatrix.rotate(rg_Point3D(0.0, 0.0, 1.0), normalOfCenterPlane);
        else
            invMatrix.rotateY(-rg_PI);
	    invMatrix.translate( centerOfGB[0] );
    }
    //  normal of center plane for a circular edge.
    else  {
	    rg_Point3D vecBg1Bg2 = centerOfGB[1]         - centerOfGB[0];
        vecBg1Bg2.normalize();

	    trMatrix.translate(-centerOfGB[0]);

        rg_REAL innerByXAxis = vecBg1Bg2.innerProduct( rg_Point3D(1.0, 0.0, 0.0) );
        if ( rg_EQ(innerByXAxis, 1.0) || rg_EQ(innerByXAxis, -1.0) )  {
        }
        else {
	        trMatrix.rotate(vecBg1Bg2, rg_Point3D(0.0, 0.0, 1.0));
	        invMatrix.rotate(rg_Point3D(1.0, 0.0, 0.0), vecBg1Bg2);
        }
        
	    invMatrix.translate( centerOfGB[0] );
    }



    rg_Point3D  trGateCenter[3] = { trMatrix*centerOfGB[0], trMatrix*centerOfGB[1], trMatrix*centerOfGB[2] };
	rg_Circle2D transformedCircle1( trGateCenter[0].getX(), trGateCenter[0].getY(), radiusOfGB[0] );
	rg_Circle2D transformedCircle2( trGateCenter[1].getX(), trGateCenter[1].getY(), radiusOfGB[1] );
	rg_Circle2D transformedCircle3( trGateCenter[2].getX(), trGateCenter[2].getY(), radiusOfGB[2] );

	rg_INT numOfCC = 0;
	rg_Circle2D tangentCircle[2];
    numOfCC = computeCircleTangentTo3CirclesOutside_GLOBAL( transformedCircle1, transformedCircle2, 
                                                            transformedCircle3, tangentCircle);   

    if (numOfCC == 1)  {      
        minTangentSphere[0].setSphere( invMatrix*tangentCircle[0].getCenterPt(), tangentCircle[0].getRadius() );
    }
    else if (numOfCC == 2){
        if ( tangentCircle[0].getRadius() < tangentCircle[1].getRadius() )  {
            minTangentSphere[0].setSphere( invMatrix*tangentCircle[0].getCenterPt(), tangentCircle[0].getRadius() );
            minTangentSphere[1].setSphere( invMatrix*tangentCircle[1].getCenterPt(), tangentCircle[1].getRadius() );
        }
        else  {
            minTangentSphere[0].setSphere( invMatrix*tangentCircle[1].getCenterPt(), tangentCircle[1].getRadius() );
            minTangentSphere[1].setSphere( invMatrix*tangentCircle[0].getCenterPt(), tangentCircle[0].getRadius() );
        }
    }
    else  {}

    return numOfCC;        
}



Sphere V::GeometryTier::computeMinSphereTangentToTwoBalls(const Sphere& ball1, const Sphere& ball2)
{
    rg_REAL radius = ball1.distanceToSphere(ball2)*0.5;

    rg_Point3D vecV1V2 = ball2.getCenter() - ball1.getCenter();
    vecV1V2.normalize();

    rg_Point3D center = ball1.getCenter() + (ball1.getRadius()+radius)*vecV1V2;

    return Sphere(center, radius);
}


