#include "Plane.h"
#include "rg_RelativeOp.h"
#include "rg_TMatrix3D.h"


Plane::Plane()
: m_distanceFromOrigin(0.0)
{
}



Plane::Plane(const rg_Point3D& normal, const rg_REAL& distanceFromOrigin)
: m_normal(normal), m_distanceFromOrigin(distanceFromOrigin)
{
}



Plane::Plane(const rg_Point3D& normal, const rg_Point3D& passingPoint)
{
    definePlaneByNormalAndPassingPoint(normal, passingPoint);
}



Plane::Plane(const rg_Point3D& passingPoint1, const rg_Point3D& passingPoint2, const rg_Point3D& passingPoint3)
{
    definePlaneByThreePoints(passingPoint1, passingPoint2, passingPoint3);
}



Plane::Plane(const Plane& plane)
: m_normal(plane.m_normal), m_distanceFromOrigin(plane.m_distanceFromOrigin)
{
}



Plane::~Plane()
{
}




rg_Point3D Plane::getNormal() const
{
    return m_normal;
}



rg_REAL    Plane::getDistanceFromOrigin() const
{
    return m_distanceFromOrigin;
}

Plane Plane::getReversedPlane() const
{
	return Plane(-m_normal, -m_distanceFromOrigin);
}


void Plane::setNormal(const rg_Point3D& normal)
{
    m_normal = normal;
}



void Plane::setDistanceFromOrigin(const rg_REAL& distanceFromOrigin)
{
    m_distanceFromOrigin = distanceFromOrigin;
}



void Plane::setPlane(const rg_Point3D& normal, const rg_REAL& distanceFromOrigin)
{
    m_normal             = normal;
    m_distanceFromOrigin = distanceFromOrigin;
}



void Plane::reverseNormal()
{
    m_normal = -m_normal;
    m_distanceFromOrigin = -m_distanceFromOrigin;
}



void Plane::definePlaneByNormalAndPassingPoint(const rg_Point3D& normal, const rg_Point3D& passingPoint)
{
    m_normal             = normal;
    m_normal.normalize();
    m_distanceFromOrigin = m_normal.innerProduct( passingPoint );
}



void Plane::definePlaneByThreePoints(const rg_Point3D& passingPoint1, const rg_Point3D& passingPoint2, const rg_Point3D& passingPoint3)
{
    if ( !rg_Point3D::areThreePointsColinear( passingPoint1, passingPoint2, passingPoint3 ) )  {
	    rg_Point3D vec21 = passingPoint2 - passingPoint1;
	    rg_Point3D vec31 = passingPoint3 - passingPoint1;

	    m_normal             = (vec21.crossProduct(vec31)).getUnitVector();
        m_distanceFromOrigin = m_normal.innerProduct( passingPoint1 );
    }
    else {
        rg_Point3D vectorPt1ToPt2 = (passingPoint2 - passingPoint1).getUnitVector();

        rg_TMatrix3D transformMatrix;
        if( !rg_EQ( vectorPt1ToPt2.innerProduct( rg_Point3D(0.0, 0.0, 1.0) ), -1) ) {
	        transformMatrix.rotate(rg_Point3D(0.0, 0.0, 1.0), vectorPt1ToPt2);
        }
        else {
            transformMatrix.rotateY(-rg_PI);
        }

        m_normal             = transformMatrix*rg_Point3D(1.0, 0.0, 0.0);
        m_distanceFromOrigin = m_normal.innerProduct( passingPoint1 );
    }


    //  Modified by Youngsong Cho (2014. 11. 7)
    //  old version
    //
	//rg_Point3D vec21 = passingPoint2 - passingPoint1;
	//rg_Point3D vec31 = passingPoint3 - passingPoint1;

	//m_normal             = (vec21.crossProduct(vec31)).getUnitVector();
    //m_distanceFromOrigin = m_normal.innerProduct( passingPoint1 );

}




rg_REAL Plane::distanceFromPoint(const rg_Point3D& point) const
{
    rg_REAL distance = m_normal.innerProduct(point) - m_distanceFromOrigin;

    return distance;
}



rg_REAL Plane::computeRelativeAngleInCCW(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3)
{
    rg_Point3D vector1 = pt1-pt2;
    vector1.normalize();
    rg_Point3D vector2 = pt3-pt2;
    vector2.normalize();

	rg_REAL angle = vector1.innerProduct(vector2);

	if ( angle > 1.0 ) 
		angle = 0.0;
	else if ( angle < -1.0 )
		angle = 2.0;
    else 
		angle = 1.0 - angle;


    rg_Point3D cross = vector1.crossProduct( vector2 );
    cross.normalize();

    rg_REAL orientation = m_normal.innerProduct( cross );
	if( rg_POS( orientation ) )  //less than PI (cross product)
		return angle;
	else
		return (4.0 - angle);
}



rg_REAL Plane::computeAngleInCCW(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3) const
{
    rg_REAL dist = pt1.distance( pt3 );

    if ( rg_ZERO(dist, resNeg10) ) {
        return 0.0;
    }
//     if ( pt1 == pt3 ) {
//         return 0.0;
//     }

    rg_Point3D vector1 = pt1-pt2;
    vector1.normalize();
    rg_Point3D vector2 = pt3-pt2;
    vector2.normalize();

	rg_REAL angleBtwVecs = vector1.angle(vector2);

    rg_Point3D cross = vector1.crossProduct( vector2 );
    cross.normalize();

    rg_REAL orientation = m_normal.innerProduct( cross );

	if( rg_POS( orientation ) )  //less than PI (cross product)
		return angleBtwVecs;
	else
		return (2.0*rg_PI - angleBtwVecs);
}



rg_REAL Plane::computeProjectedAngleInCCW(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3) const
{
    if ( pt1 == pt3 ) {
        return 0.0;
    }

    rg_Point3D point[3];
    point[0] = projectPointOnPlane( pt1 );
    point[1] = projectPointOnPlane( pt2 );
    point[2] = projectPointOnPlane( pt3 );

    return computeAngleInCCW( point[0], point[1], point[2] );
}



rg_Point3D Plane::reflectPoint(const rg_Point3D& point) const
{
    rg_REAL signedDistFromPoint = m_normal.innerProduct(point) - m_distanceFromOrigin;
    rg_REAL distFromPoint       = rg_ABS( signedDistFromPoint );
    rg_Point3D reflection;
    if ( rg_ZERO( signedDistFromPoint ) ) {
        reflection = point;
    }
    else if ( signedDistFromPoint > 0 ) {
        reflection = point - 2.0*distFromPoint*m_normal;
    }
    else {
        reflection = point + 2.0*distFromPoint*m_normal;
    }

    return reflection;
}



rg_Point3D Plane::projectPointOnPlane(const rg_Point3D& point) const
{
    rg_REAL distFromPoint = m_normal.innerProduct(point) - m_distanceFromOrigin;

    rg_Point3D projectedPoint = point - (distFromPoint*m_normal);

    return projectedPoint;
}


rg_FLAG Plane::computeIntersectionWithPlane(const Plane& plane, Line3D& line) const
{
	// two planes are parallel
	rg_Point3D normalVec = plane.m_normal;
	if (rg_EQ(rg_ABS(m_normal.innerProduct(normalVec)), 1.0))
		return rg_FALSE;


	// otherwise they have an intersection as an infinite line
	rg_Point3D directionVec = m_normal.crossProduct(normalVec);
	// passing point "P" is represented by a linear combination of two normal vectors
	// let P = coeff1 * m_normal + coeff2 * normalVec

	rg_REAL innerProd = m_normal.innerProduct(normalVec);
	rg_REAL squardInnerProd = innerProd * innerProd;
	rg_REAL squardNorm1 = m_normal.squaredMagnitude();
	rg_REAL squardNorm2 = normalVec.squaredMagnitude();
	rg_REAL multSquardNorms = squardNorm1 * squardNorm2;
	rg_REAL distFromOrg1 = m_distanceFromOrigin;
	rg_REAL distFromOrg2 = plane.m_distanceFromOrigin;
	rg_REAL coeff1 = (distFromOrg2 * innerProd - distFromOrg1 * squardNorm2) 
		           / (squardInnerProd - multSquardNorms);
	rg_REAL coeff2 = (distFromOrg1 * innerProd - distFromOrg2 * squardNorm1) 
		           / (squardInnerProd - multSquardNorms);

	rg_Point3D passingPt = coeff1 * m_normal + coeff2 * normalVec;
	line.set(m_normal.crossProduct(normalVec), passingPt);


	return rg_TRUE;
}



Plane& Plane::operator =(const Plane& plane)
{
    if ( this == &plane )
        return *this;

    m_normal             = plane.m_normal;
    m_distanceFromOrigin = plane.m_distanceFromOrigin;

    return *this;
}



rg_BOOL Plane::operator==(const Plane& plane) const
{
    if ( m_normal == plane.m_normal ) {
        if ( rg_EQ(m_distanceFromOrigin, plane.m_distanceFromOrigin ) ){
            return rg_TRUE;
        }
    }
    else {
        if ( m_normal == -plane.m_normal ) {
            if ( rg_EQ(m_distanceFromOrigin, -plane.m_distanceFromOrigin ) ){
                return rg_TRUE;
            }
        }
    }

    return rg_FALSE;
}




rg_BOOL Plane::operator!=(const Plane& plane) const
{
    if ( operator==(plane) ) {
        return rg_FALSE;
    }
    else {
        return rg_TRUE;
    }
}



rg_FLAG Plane::isOnNormalSide(const rg_Point3D& point) const
{
	if(rg_POS(m_normal.innerProduct(point) - m_distanceFromOrigin))
		return rg_TRUE;
	else
		return rg_FALSE;
}



rg_FLAG Plane::isOnOppositeNormalSide(const rg_Point3D& point) const
{
	if(rg_NEG(m_normal.innerProduct(point) - m_distanceFromOrigin))
		return rg_TRUE;
	else
		return rg_FALSE;
}



rg_FLAG Plane::isOnThisPlane(const rg_Point3D& point) const
{
	if(rg_ZERO(m_normal.innerProduct(point) - m_distanceFromOrigin))
		return rg_TRUE;
	else
		return rg_FALSE;
}



rg_FLAG Plane::isOnNormalSideOrOnThisPlane(const rg_Point3D& point) const
{
	if(rg_NNEG(m_normal.innerProduct(point) - m_distanceFromOrigin))
		return rg_TRUE;
	else
		return rg_FALSE;
}



rg_FLAG Plane::isOnOppositeNormalSideOrOnThisPlane(const rg_Point3D& point) const
{
	if(rg_NPOS(m_normal.innerProduct(point) - m_distanceFromOrigin))
		return rg_TRUE;
	else
		return rg_FALSE;
}



rg_BOOL  Plane::isThereIntersectionWith(const Plane& plane) const
{
    if ( (*this) == plane ) {
        return rg_TRUE;
    }
    else {
        rg_REAL innerProduct = m_normal.innerProduct( plane.m_normal );
        if ( rg_EQ( innerProduct, 1.0, resNeg10) || rg_EQ( innerProduct, -1.0, resNeg10) ) {
		    return rg_FALSE;
        }
        else {
            return rg_TRUE;
        }
    }
}



rg_BOOL Plane::intersect(const Plane& plane, Line3D& line) const
{
    if ( operator==(plane) ) {
        return rg_TRUE;
    }
    else {
        if ( isParallelTo( plane ) ) {
            //  Two planes are parallel.
		    return rg_FALSE;
        }
        else {
	        // otherwise they have an intersection as an infinite line
	        rg_Point3D directionVec = m_normal.crossProduct(plane.m_normal);
	        // passing point "P" is represented by a linear combination of two normal vectors
	        // let P = coeff1 * m_normal + coeff2 * normalVec

	        rg_REAL innerProd       = m_normal.innerProduct( plane.m_normal );
	        rg_REAL squardInnerProd = innerProd * innerProd;
	        rg_REAL squardNorm1     = m_normal.squaredMagnitude();
	        rg_REAL squardNorm2     = plane.m_normal.squaredMagnitude();
	        rg_REAL multSquardNorms = squardNorm1 * squardNorm2;
	        rg_REAL distFromOrg1    = m_distanceFromOrigin;
	        rg_REAL distFromOrg2    = plane.m_distanceFromOrigin;
	        rg_REAL coeff1 = ((distFromOrg2 * innerProd) - (distFromOrg1 * squardNorm2))
		                      / (squardInnerProd - multSquardNorms);
	        rg_REAL coeff2 = ((distFromOrg1 * innerProd) - (distFromOrg2 * squardNorm1))
		                      / (squardInnerProd - multSquardNorms);

	        rg_Point3D passingPt = coeff1 * m_normal + coeff2 * plane.m_normal;

	        line.set(directionVec, passingPt);

	        return rg_TRUE;
        }
    }


    /*
    if ( operator==(plane) ) {
        return rg_TRUE;
    }
    else {
        rg_REAL innerProduct = m_normal.innerProduct( plane.m_normal );
        if ( isParallelTo( plane ) ) {
            //  Two planes are parallel.
		    return rg_FALSE;
        }
        else {
	        rg_Point3D directionVec = m_normal.crossProduct(plane.m_normal);
            directionVec.normalize();

            //                             /
            //                           /
            //         p1    p3        / 
            // --------*-----*-------*-----------------
            //               |     / psPt
            //               |   /
            //               | /
            //               * p2
            //              /
            //
            //   theta = angle of two intersecting plane
            //         = angle (p3, psPt, p2)
            //         = angle (p1, p2, p3)
            //   angle (p3, p1, p2) = rg_PI/2 - theta
            //
            //   Q) Is a vector from p1 to p3, v_13, perpendicular to directionVec? 
            //      If not, we can be find intersection by the following.
            //      A) Yes! 
            //          v_21 = plane.m_normal
            //          v_23 = m_normal
            //          therefore, the plane defined by p1, p2, p3 is perperdicular to directionVec.

            rg_Point3D p1 = projectPointOnPlane( rg_Point3D(0.0, 0.0, 0.0) );
            rg_Point3D p2 = plane.projectPointOnPlane( p1 );
            rg_Point3D p3 = projectPointOnPlane( p2 );

            rg_REAL    theta    = acos( innerProduct );
            rg_REAL    dist     = p1.distance( p2 );
            rg_REAL    movement = dist / cos( rg_PI*0.5 - theta );
            rg_Point3D vec13    = p3 - p1;
            vec13.normalize();

            rg_Point3D passingPoint = p1 + movement*vec13;

	        line.set(directionVec, passingPoint);




	        // otherwise they have an intersection as an infinite line
	        rg_Point3D directionVec1 = m_normal.crossProduct(plane.m_normal);
	        // passing point "P" is represented by a linear combination of two normal vectors
	        // let P = coeff1 * m_normal + coeff2 * normalVec

	        rg_REAL innerProd = m_normal.innerProduct( plane.m_normal );
	        rg_REAL squardInnerProd = innerProd * innerProd;
	        rg_REAL squardNorm1 = m_normal.squaredMagnitude();
	        rg_REAL squardNorm2 = plane.m_normal.squaredMagnitude();
	        rg_REAL multSquardNorms = squardNorm1 * squardNorm2;
	        rg_REAL distFromOrg1 = m_distanceFromOrigin;
	        rg_REAL distFromOrg2 = plane.m_distanceFromOrigin;
	        rg_REAL coeff1 = (distFromOrg2 * innerProd - distFromOrg1 * squardNorm2) 
		                   / (squardInnerProd - multSquardNorms);
	        rg_REAL coeff2 = (distFromOrg1 * innerProd - distFromOrg2 * squardNorm1) 
		                   / (squardInnerProd - multSquardNorms);

	        rg_Point3D passingPt = coeff1 * m_normal + coeff2 * plane.m_normal;

	        line.set(directionVec, passingPt);


            rg_REAL dista1[2] = { distanceFromPoint(passingPoint), plane.distanceFromPoint(passingPoint) };
            rg_REAL dista2[2] = { distanceFromPoint(passingPt), plane.distanceFromPoint(passingPt) };


	        return rg_TRUE;
        }
    }
    */
}



rg_INT Plane::intersect(const LineSegment3D& lineSegment, rg_Point3D& intersectionPoint) const
{
    return intersectLineSegment(lineSegment.getStartPt(), lineSegment.getEndPt(), intersectionPoint);
}



rg_INT Plane::intersectLineSegment(const rg_Point3D& pt1, const rg_Point3D& pt2, rg_Point3D& intersectionPoint) const
{
    rg_REAL distPt1 = distanceFromPoint(pt1);
    rg_REAL distPt2 = distanceFromPoint(pt2);

    if (distPt1*distPt2 > 0) {
        return 0;
    }

	// inserted by Joonghyun
	// has no intersection when the line segment or its vertex on the plane
    if (rg_ZERO(distPt1, resNeg10) || rg_ZERO(distPt2, resNeg10)) {
        return 0;
    }

    distPt1 = rg_ABS( distPt1 );
    distPt2 = rg_ABS( distPt2 );

    rg_REAL dist    = pt1.distance(pt2);

    rg_Point3D vecPt1Pt2 = pt2 - pt1;
    vecPt1Pt2.normalize();

    rg_REAL len = (distPt1*dist)/(distPt1+distPt2);

    intersectionPoint = pt1 + (len*vecPt1Pt2);

    return 1;	
}



rg_Point3D Plane::intersectLineSegment(const rg_Point3D& pt1, const rg_Point3D& pt2) const
{
    rg_REAL distPt1 = distanceFromPoint(pt1);
    rg_REAL distPt2 = distanceFromPoint(pt2);

    if (distPt1*distPt2 > 0) {
        return rg_Point3D();
    }

	// inserted by Joonghyun
	// has no intersection when the line segment or its vertex on the plane
    if (rg_ZERO(distPt1, resNeg10) || rg_ZERO(distPt2, resNeg10)) {
        return rg_Point3D();
    }

    distPt1 = rg_ABS( distPt1 );
    distPt2 = rg_ABS( distPt2 );

    rg_REAL dist    = pt1.distance(pt2);

    rg_Point3D vecPt1Pt2 = pt2 - pt1;
    vecPt1Pt2.normalize();

    rg_REAL len = (distPt1*dist)/(distPt1+distPt2);
    rg_Point3D intersect = pt1 + (len*vecPt1Pt2);

    return intersect;	
}



void Plane::projectLineSegmentOnPlane(const LineSegment3D& targetLineSegment, LineSegment3D& projection) const
{
	rg_Point3D sPt = projectPointOnPlane(targetLineSegment.getStartPt());
	rg_Point3D ePt = projectPointOnPlane(targetLineSegment.getEndPt());

	projection.setPoints(sPt, ePt);
}



rg_FLAG Plane::computeIntersectionWithLineSegment(const rg_Point3D& pt1, const rg_Point3D& pt2, rg_Point3D& intersecPt) const
{
    rg_REAL distPt1 = distanceFromPoint(pt1);
    rg_REAL distPt2 = distanceFromPoint(pt2);

    if (rg_POS(distPt1*distPt2)) {
        return rg_FALSE;
    }

	// inserted by Joonghyun
	// has no intersection when the line segment or its vertex on the plane
    if (rg_ZERO(distPt1) || rg_ZERO(distPt2)) {
        return rg_FALSE;
    }

    distPt1 = rg_ABS( distPt1 );
    distPt2 = rg_ABS( distPt2 );

    rg_REAL dist = pt1.distance(pt2);

    rg_Point3D vecPt1Pt2 = pt2 - pt1;
    vecPt1Pt2.normalize();

    rg_REAL len = (distPt1*dist)/(distPt1+distPt2);
    intersecPt = pt1 + (len*vecPt1Pt2);
	return rg_TRUE;
}

rg_FLAG Plane::computeIntersectionWithLineSegment(const LineSegment3D& linesegment, 
												  rg_Point3D& intersecPt, 
												  rg_REAL& intersecParam,
												  PosOnLineSegOrArc& pos) const
{
	rg_Point3D dirVec = linesegment.getDirVec();
	rg_REAL denom = m_normal.innerProduct(dirVec);
	if(rg_ZERO(denom))
		return rg_FALSE;

	rg_Point3D sPt = linesegment.getStartPt();
	intersecParam = (m_distanceFromOrigin - m_normal.innerProduct(sPt)) / denom;
	if(rg_LT(intersecParam, 0.0) || rg_GT(intersecParam, 1.0))
		return rg_FALSE;

	intersecPt = sPt + intersecParam * dirVec;
	if(rg_ZERO(intersecParam))
		pos = START_PT;
	else if(rg_EQ(intersecParam, 1.0))
		pos = END_PT;
	else
		pos = NONEXTREME_PT;

	return rg_TRUE;
}

rg_FLAG Plane::computeIntersectionWithLineSegment(const LineSegment3D& linesegment, 
												  rg_Point3D& intersecPt, 
												  rg_REAL& intersecParam) const
{
	rg_Point3D dirVec = linesegment.getDirVec();
	rg_REAL denom = m_normal.innerProduct(dirVec);
	if(rg_ZERO(denom))
		return rg_FALSE;

	rg_Point3D sPt = linesegment.getStartPt();
	intersecParam = (m_distanceFromOrigin - m_normal.innerProduct(sPt)) / denom;
	if(rg_LT(intersecParam, 0.0) || rg_GT(intersecParam, 1.0))
		return rg_FALSE;

	intersecPt = sPt + intersecParam * dirVec;
	return rg_TRUE;
}

rg_FLAG Plane::computeIntersectionWithLineSegment(const LineSegment3D& linesegment, 
												  rg_REAL& intersecParam) const
{
	rg_Point3D dirVec = linesegment.getDirVec();
	rg_REAL denom = m_normal.innerProduct(dirVec);
	if(rg_ZERO(denom))
		return rg_FALSE;

	rg_Point3D sPt = linesegment.getStartPt();
	intersecParam = (m_distanceFromOrigin - m_normal.innerProduct(sPt)) / denom;
	if(rg_LT(intersecParam, 0.0) || rg_GT(intersecParam, 1.0))
		return rg_FALSE;
	else
		return rg_TRUE;
}

rg_FLAG Plane::isParallelTo(const Plane& plane) const
{
    rg_REAL innerProduct = rg_ABS( m_normal.innerProduct(plane.m_normal) );
    if( rg_EQ(innerProduct, 1.0, resNeg10) ) {
		return rg_TRUE;
    }
    else {
		return rg_FALSE;
    }
}

rg_FLAG Plane::isCoincidentWith(const Plane& plane) const
{
	if(operator==(plane))
		return rg_TRUE;
	Plane reversedPlane = getReversedPlane();
	if(reversedPlane == plane)
		return rg_TRUE;
	else
		return rg_FALSE;
}

