#include "Arc3D.h"
#include "rg_RelativeOp.h"

#include "SquaredDistFuncLineSegArc.h"
#include <float.h>


Arc3D::Arc3D()
: Circle3D()
{
}



Arc3D::Arc3D(const Circle3D& circle, const rg_Point3D& startPt, const rg_Point3D& endPt)
: Circle3D(circle),
  m_startPoint(startPt), m_endPoint(endPt)
{
}



Arc3D::Arc3D(const Arc3D& arc)
: Circle3D(arc),
  m_startPoint(arc.m_startPoint), m_endPoint(arc.m_endPoint)
{
}



Arc3D::~Arc3D()
{
}


/*
rg_Point3D  Arc3D::getStartPoint() const
{
    return m_startPoint;
}



rg_Point3D  Arc3D::getEndPoint() const
{
    return m_endPoint;
}
*/



void        Arc3D::setStartPoint(const rg_Point3D& startPt)
{
    m_startPoint = startPt;
}



void        Arc3D::setEndPoint(const rg_Point3D& endPt)
{
    m_endPoint = endPt;
}



void        Arc3D::setPoints(const rg_Point3D& startPt, const rg_Point3D& endPt)
{
    m_startPoint = startPt;
    m_endPoint   = endPt;
}



void        Arc3D::setArc(const Circle3D& circle, const rg_Point3D& startPt, const rg_Point3D& endPt)
{
    Circle3D::operator =(circle);
    m_startPoint = startPt;
    m_endPoint   = endPt;
}



Arc3D&      Arc3D::operator =(const Arc3D& arc)
{
    if ( this == &arc ) {
        return *this;
    }

    Circle3D::operator =(arc);

    m_startPoint = arc.m_startPoint;
    m_endPoint   = arc.m_endPoint;

    return *this;
}



rg_BOOL Arc3D::operator==(const Arc3D& arc) const
{
    rg_BOOL isSameArc = rg_FALSE;
    if ( Circle3D::operator ==(arc) ) {
        if ( m_startPoint == m_endPoint && arc.m_startPoint == arc.m_endPoint ) {
            isSameArc = rg_TRUE;
        }
        else {
            rg_Point3D normal    = getNormal();
            rg_Point3D arcNormal = arc.getNormal();

            if ( normal == arcNormal )  {
                if (    m_startPoint == arc.m_startPoint 
                     && m_endPoint   == arc.m_endPoint ) {
                    isSameArc = rg_TRUE;
                }
            }
            else if ( normal == -arcNormal ) {
                if (    m_startPoint == arc.m_endPoint 
                     && m_endPoint   == arc.m_startPoint ) {
                    isSameArc = rg_TRUE;
                }
            }
            else {
                // do nothing
            }
        }
    }

    return isSameArc;
}



rg_BOOL Arc3D::isOnThisArc(const rg_Point3D& point) const
{
    if ( point == m_startPoint || point == m_endPoint ) {
        return rg_TRUE;
    }

    if ( isOnCircle(point) ) {
        if ( m_startPoint == m_endPoint ) {
            //  it means that this arc is a circle.
            return rg_TRUE;
        }

        Plane planeContainingArc = getPlaneContainingThisCircle();

        rg_Point3D center = getCenter();
        rg_REAL angleStartToEnd = planeContainingArc.computeAngleInCCW(m_startPoint, center, m_endPoint);
        rg_REAL angleStartToPt  = planeContainingArc.computeAngleInCCW(m_startPoint, center, point);

        if ( angleStartToEnd > angleStartToPt ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_BOOL Arc3D::isMinorArc() const
{
    rg_REAL angle = evaluateAngle();

    if ( angle < rg_PI ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_BOOL Arc3D::isContainedIn(const Arc3D& bigArc) const
{
    rg_BOOL isThisArcContainedInBigArc = rg_FALSE;

    if ( Circle3D::operator ==(bigArc) ) {
        Plane planeOfBigArc = bigArc.getPlaneContainingThisCircle();
        rg_REAL angleToEndOfBigArc    = planeOfBigArc.computeAngleInCCW(bigArc.m_startPoint, bigArc.m_center, bigArc.m_endPoint);
        rg_REAL angleToStartOfThisArc = planeOfBigArc.computeAngleInCCW(bigArc.m_startPoint, bigArc.m_center, m_startPoint);
        rg_REAL angleToEndOfThisArc   = planeOfBigArc.computeAngleInCCW(bigArc.m_startPoint, bigArc.m_center, m_endPoint);

        if ( rg_LE(angleToStartOfThisArc, angleToEndOfThisArc) ) {
            if ( rg_LE(angleToEndOfThisArc, angleToEndOfBigArc) ) {
                isThisArcContainedInBigArc = rg_TRUE;
            }
        }
    }

    return isThisArcContainedInBigArc;
}



ARCType Arc3D::getArcType() const
{
    rg_REAL angle = evaluateAngle();

    if ( angle < rg_PI ) {
        return MINOR_ARC;
    }
    else {
        return MAJOR_ARC;
    }
}


void Arc3D::reverse()
{
    rg_Point3D normal = getNormal();
    setNormal( -normal );

    rg_Point3D tempPoint = m_startPoint;
    m_startPoint = m_endPoint;
    m_endPoint   = tempPoint;
}



Arc3D Arc3D::getReversedArc() const
{
    Arc3D reversedArc(*this);
    reversedArc.reverse();

    return reversedArc;
}



Circle3D Arc3D::getCircle() const
{
	Circle3D circle(m_center, m_radius, m_normal);
	return circle;
}



rg_INT Arc3D::intersect(const Plane& plane, rg_Point3D* intersectPts) const
{
	rg_INT numIntersection = Circle3D::intersect(plane, intersectPts);

    if ( m_startPoint == m_endPoint ) {
        //  it means that this arc is a circle.
        return numIntersection;
    }

    if( numIntersection == 2 ) {
        rg_INT i = 0;
        for ( i=0; i<2; i++ ) {
            if ( intersectPts[i] == m_startPoint ) {
                intersectPts[i] = m_startPoint;
            }
            if ( intersectPts[i] == m_endPoint ) {
                intersectPts[i] = m_endPoint;
            }
        }

        rg_BOOL isPt1OnArc = isOnThisArc( intersectPts[0] );
        rg_BOOL isPt2OnArc = isOnThisArc( intersectPts[1] );

        if ( isPt1OnArc && isPt2OnArc ) {
            Plane planeContainingArc = getPlaneContainingThisCircle();

            rg_Point3D center = getCenter();
            rg_REAL angleStartToPt1 = planeContainingArc.computeAngleInCCW(m_startPoint, center, intersectPts[0]);
            rg_REAL angleStartToPt2 = planeContainingArc.computeAngleInCCW(m_startPoint, center, intersectPts[1]);
            if ( angleStartToPt1 > angleStartToPt2 ) {
                rg_Point3D temp = intersectPts[0];
                intersectPts[0] = intersectPts[1];
                intersectPts[1] = temp;
            }

            return 2;
        }
        else if ( isPt1OnArc ) {
            return 1;
        }
        else if ( isPt2OnArc ) {
			intersectPts[ 0 ] = intersectPts[ 1 ];
            return 1;
        }
        else {
            return 0;
        }
    }
    else if ( numIntersection == 1 ) {
        if ( intersectPts[0] == m_startPoint ) {
            intersectPts[0] = m_startPoint;
        }
        if ( intersectPts[0] == m_endPoint ) {
            intersectPts[0] = m_endPoint;
        }

        if ( isOnThisArc( intersectPts[0] ) ) {
            return 1; 
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }
}



rg_INT Arc3D::intersect(const Circle3D& circle, rg_Point3D* intersectPts) const
{
	rg_INT numIntersection = Circle3D::intersect(circle, intersectPts);

    if ( m_startPoint == m_endPoint ) {
        //  it means that this arc is a circle.
        return numIntersection;
    }

    if( numIntersection == 2 ) {
        rg_BOOL isPt1OnArc = isOnThisArc( intersectPts[0] );
        rg_BOOL isPt2OnArc = isOnThisArc( intersectPts[1] );

        if ( isPt1OnArc && isPt2OnArc ) {
            Plane planeContainingArc = getPlaneContainingThisCircle();

            rg_Point3D center = getCenter();
            rg_REAL angleStartToPt1 = planeContainingArc.computeAngleInCCW(m_startPoint, center, intersectPts[0]);
            rg_REAL angleStartToPt2 = planeContainingArc.computeAngleInCCW(m_startPoint, center, intersectPts[1]);
            if ( angleStartToPt1 > angleStartToPt2 ) {
                rg_Point3D temp = intersectPts[0];
                intersectPts[0] = intersectPts[1];
                intersectPts[1] = temp;
            }

            return 2;
        }
        else if ( isPt1OnArc ) {
            return 1;
        }
        else if ( isPt2OnArc ) {
			intersectPts[ 0 ] = intersectPts[ 1 ];
            return 1;
        }
        else {
            return 0;
        }
    }
    else if ( numIntersection == 1 ) {
        if ( isOnThisArc( intersectPts[0] ) ) {
            return 1; 
        }
        else {
            return 0;
        }
    }
    else {
        // In this case, numIntersection is 0 or INT_MAX.
        // INT_MAX means that this arc is completely contained in the circle.

        return numIntersection;
    }
}



rg_INT Arc3D::intersect(const Arc3D &arc, rg_Point3D* intersectPts, const rg_REAL tolerance) const
{    
	rg_INT numIntersections = 0;

    
    //  Special Case: Arc_A and arc_B are not coplanar.
    //                Arc_A : this arc
    //                Arc_B : input arc
    //  The start or end points of Arc_B is on Arc_A.
    //  The following codes may reduce a numerical error 
    //  when the intersections of two arcs are the start points, end points or both of them of an input arc.
    if ( Circle3D::operator!=(arc) ) {
        //  There may be 0, 1, or 2 intersections.
        if ( isOnThisArc( arc.m_startPoint ) ) {
            intersectPts[0] = arc.m_startPoint;
            numIntersections = 1;
        }
    
        if ( isOnThisArc( arc.m_endPoint ) ) {
            intersectPts[numIntersections] = arc.m_endPoint;
            numIntersections++;
        }
    
        Plane      section;
        rg_Point3D reflection;
        switch ( numIntersections ) {
            case 0:
                //  Two arcs do not intersect at extreme points of the input arc.
                //  However, we don't know yet whether two arcs intersect or not.
                numIntersections = 0;
                break;
    
            case 1:
                //  One of two extreme points of the input arc is on this arc.
                //  That is, two arcs intersect at one of two extreme points of the input arc.
                section.definePlaneByNormalAndPassingPoint( m_normal.crossProduct( arc.m_normal ), m_center );
                reflection = section.reflectPoint( intersectPts[0] );
    
                if ( arc.isOnThisArc(reflection) ) {
                    intersectPts[numIntersections] = reflection;
                    numIntersections++;
                }
                break;
    
            case 2:
                //  Two arcs intersect at both of two extreme points of the input arc.
                numIntersections = 2;
                break;
    
            default:
                break;
        }
    
        if ( numIntersections > 0 ) {
            return numIntersections;
        }
    }
    
    

    //  General Case:
    rg_Point3D intersectBtwCircles[2];
	rg_INT     numOfIntersectionsBtwCircles = Circle3D::intersect(arc, intersectBtwCircles, tolerance);

    Arc3D inputArc;
	switch( numOfIntersectionsBtwCircles ) {
	    case 0:
		    numIntersections = 0;
		    break;

	    case 1:
            if( isOnThisArc(intersectBtwCircles[0]) && arc.isOnThisArc(intersectBtwCircles[0]) ) {
                intersectPts[numIntersections] = intersectBtwCircles[0];
			    numIntersections = 1;
            }
		    break;

	    case 2:
		    if( isOnThisArc(intersectBtwCircles[0]) && arc.isOnThisArc(intersectBtwCircles[0]) ) {
                intersectPts[numIntersections] = intersectBtwCircles[0];
			    numIntersections = 1;
		    }

		    if( isOnThisArc(intersectBtwCircles[1]) && arc.isOnThisArc(intersectBtwCircles[1]) ) {
                intersectPts[numIntersections] = intersectBtwCircles[1];
			    numIntersections++;
		    }
	        break;

	    case rg_INT_INFINITY:
            //  Two arcs are based on the same circles. Therefore, they can be overlaped.
            inputArc = arc;
            if ( m_normal != arc.m_normal ) {
                inputArc.reverse();
            }


            if ( inputArc.isContainedIn( *this ) ) {
                //  The intersection of two arcs is identical to an input arc. 
                //  The input arc is a proper subset of this arc.
                intersectPts[0]  = inputArc.m_startPoint;
                intersectPts[1]  = inputArc.m_endPoint;

			    numIntersections = rg_INT_INFINITY;
            }
            else if ( isContainedIn( inputArc ) ) {
                //  The intersection of two arcs is identical to this arc. 
                //  This arc is a proper subset of the input arc.
                intersectPts[0]  = m_startPoint;
                intersectPts[1]  = m_endPoint;

			    numIntersections = rg_INT_INFINITY;
            }
            else if ( isOnThisArc(inputArc.m_startPoint) && isOnThisArc(inputArc.m_endPoint) ) {
                intersectPts[0]  = inputArc.m_startPoint;
                intersectPts[1]  = m_endPoint;
                intersectPts[2]  = m_startPoint;
                intersectPts[3]  = inputArc.m_endPoint;

			    numIntersections = rg_INT_INFINITY;
            }
            else if ( isOnThisArc(inputArc.m_startPoint) ) {
                //  Two arcs are partially overlapped from a start point of the input arc to an end point of this arc.
                intersectPts[0]  = inputArc.m_startPoint;
                intersectPts[1]  = m_endPoint;

			    numIntersections = rg_INT_INFINITY;
            }
            else if ( isOnThisArc(inputArc.m_endPoint) ) {
                //  Two arcs are partially overlapped from a start point of this arc to an end point of the input arc.
                intersectPts[0]  = m_startPoint;
                intersectPts[1]  = inputArc.m_endPoint;

			    numIntersections = rg_INT_INFINITY;
            }
		    else {
			    numIntersections = 0;
		    }

	        break;

	    default:
	        break;
	}

	return numIntersections;    
}




rg_REAL Arc3D::evaluateArcLength() const
{
    rg_REAL angle     = evaluateAngle();
    rg_REAL arcLength = angle*getRadius();

    return arcLength;
}



rg_REAL Arc3D::evaluateAngle() const
{
    rg_REAL dist = m_startPoint.distance( m_endPoint );

//     if ( m_startPoint == m_endPoint ) {
    if ( rg_ZERO(dist, resNeg10) ) {
        return 2.0*rg_PI;
    }
    
    Plane planeContainingArc = getPlaneContainingThisCircle();

    rg_Point3D center = getCenter();
    rg_REAL angleStartToEnd = planeContainingArc.computeAngleInCCW(m_startPoint, center, m_endPoint);

    return angleStartToEnd;
}



rg_Point3D  Arc3D::evaluateMiddlePoint() const
{
    rg_REAL angle = evaluateAngle();

    rg_Point3D center = getCenter();
    rg_REAL    radius = getRadius();
        
    rg_Point3D vecCenterToStart = m_startPoint - center;
    rg_Point3D vecCenterToEnd   = m_endPoint - center;
    rg_Point3D angularBisector = (vecCenterToStart + vecCenterToEnd)*0.5;
    angularBisector.normalize();

    rg_Point3D middlePoint;
    if ( angle <= rg_PI ) {
        middlePoint = center + (radius * angularBisector);
    }
    else {
        middlePoint = center - (radius * angularBisector);
    }

    return middlePoint;
}



void Arc3D::evaluatePointsOnArc(const rg_REAL& unitLength, rg_dList<rg_Point3D>& pointsOnArc)const
{
    if ( m_startPoint == m_endPoint ) {
        Circle3D::evaluatePointsOnCircle(unitLength, m_startPoint, pointsOnArc);
    }
    else {
        rg_REAL arcLength = evaluateArcLength();

        rg_INT  numLineSegment = (rg_INT) ceil( arcLength/unitLength );
        if ( numLineSegment == 1 ) {
            if ( isMinorArc() ) {
                pointsOnArc.add( m_startPoint );
                pointsOnArc.add( m_endPoint );
            }
            else {
                pointsOnArc.add( m_startPoint );
                pointsOnArc.add( evaluateMiddlePoint() );
                pointsOnArc.add( m_endPoint );
            }
        }
        else {
            rg_Point3D normal = getNormal();
            rg_Point3D center = getCenter();
            rg_REAL    radius = getRadius();

            rg_REAL stepAngle = evaluateAngle()/numLineSegment;
            rg_REAL magitudeOfTangent = radius*tan( stepAngle );

            pointsOnArc.add( m_startPoint );
            rg_Point3D currPoint = m_startPoint;
            for ( rg_INT i=0; i<(numLineSegment-1); i++ ) {
                rg_Point3D tangent = normal.crossProduct( currPoint-center );
                tangent.normalize();

                rg_Point3D sector = (currPoint + (magitudeOfTangent*tangent)) - center;
                sector.normalize();

                currPoint = center + (radius*sector);
                pointsOnArc.add( currPoint );
            }
            pointsOnArc.add( m_endPoint );
        }
    }
}



void Arc3D::evaluatePointsOnArc(const rg_INT&  depth, rg_Point3D*& pointsOnArc) const
{
    if ( depth == 0 ) {
        pointsOnArc = new rg_Point3D[2];
        pointsOnArc[0] = m_startPoint;
        pointsOnArc[1] = m_endPoint;
        return;
    }
    else if ( depth > 0 ) {
        rg_INT numLineSegment = (rg_INT) pow(2., depth);

        pointsOnArc = new rg_Point3D[numLineSegment+1];
        pointsOnArc[0]              = m_startPoint;
        pointsOnArc[numLineSegment] = m_endPoint;

        rg_Point3D normal = getNormal();
        rg_Point3D center = getCenter();
        rg_REAL    radius = getRadius();

        rg_REAL stepAngle = evaluateAngle()/numLineSegment;
        rg_REAL magitudeOfTangent = radius*tan( stepAngle );

        rg_Point3D currPoint = m_startPoint;
        for ( rg_INT i=1; i<numLineSegment; i++ ) {
            rg_Point3D tangent = normal.crossProduct( currPoint-center );
            tangent.normalize();

            rg_Point3D sector = (currPoint + (magitudeOfTangent*tangent)) - center;
            sector.normalize();

            currPoint = center + (radius*sector);
            pointsOnArc[i] = currPoint;
        }
    }
    else {
    }
}



void Arc3D::evaluatePointsOnArc(const rg_INT&  depth, rg_dList<rg_Point3D>& pointsOnArc) const
{
    if ( depth == 0 ) {
        pointsOnArc.add( m_startPoint );
        pointsOnArc.add( m_endPoint );
    }
    else {
        rg_INT numLineSegment = (rg_INT) pow(2., depth);

        pointsOnArc.add( m_startPoint );

        rg_Point3D normal = getNormal();
        rg_Point3D center = getCenter();
        rg_REAL    radius = getRadius();

        rg_REAL stepAngle = evaluateAngle()/numLineSegment;
        rg_REAL magitudeOfTangent = radius*tan( stepAngle );

        rg_Point3D currPoint = m_startPoint;
        for ( rg_INT i=1; i<numLineSegment; i++ ) {
            rg_Point3D tangent = normal.crossProduct( currPoint-center );
            tangent.normalize();

            rg_Point3D sector = (currPoint + (magitudeOfTangent*tangent)) - center;
            sector.normalize();

            currPoint = center + (radius*sector);
            pointsOnArc.add( currPoint );

        }

        pointsOnArc.add( m_endPoint );
    }
}



void Arc3D::evaluatePointsOnArcGivenNumPts(const rg_INT&  numPts, rg_Point3D*& pointsOnArc) const
{
    rg_INT numLineSegment = numPts-1;

    pointsOnArc = new rg_Point3D[numLineSegment+1];

    rg_Point3D normal = getNormal();
    rg_Point3D center = getCenter();
    rg_REAL    radius = getRadius();

	rg_Point3D uVector = (m_startPoint - center).getUnitVector();
	rg_Point3D vVector = normal.crossProduct(uVector);


    rg_REAL stepAngle;
    if ( m_startPoint == m_endPoint ) {
        stepAngle = 2.0*rg_PI/numLineSegment;
    }
    else {
        stepAngle= evaluateAngle()/numLineSegment;
    }


    pointsOnArc[0]       = m_startPoint;

    for ( rg_INT i=1; i<numLineSegment; i++ ) {
	    rg_REAL    theta = i*stepAngle;
        rg_Point3D point = center + radius * (cos(theta) * uVector + sin(theta) * vVector);
        pointsOnArc[i]   = point;
    }

    pointsOnArc[numLineSegment] = m_endPoint;

}



void Arc3D::evaluatePointsOnArcGivenNumPts(const rg_INT&  numPts, rg_dList<rg_Point3D>& pointsOnArc) const
{
    rg_INT numLineSegment = numPts-1;

    rg_Point3D normal = getNormal();
    rg_Point3D center = getCenter();
    rg_REAL    radius = getRadius();

	rg_Point3D uVector = (m_startPoint - center).getUnitVector();
	rg_Point3D vVector = normal.crossProduct(uVector);


    rg_REAL stepAngle;
    if ( m_startPoint == m_endPoint ) {
        stepAngle = 2.0*rg_PI/numLineSegment;
    }
    else {
        stepAngle= evaluateAngle()/numLineSegment;
    }

    pointsOnArc.add( m_startPoint );

    for ( rg_INT i=1; i<numLineSegment; i++ ) {
	    rg_REAL    theta = i*stepAngle;
        rg_Point3D point = center + radius * (cos(theta) * uVector + sin(theta) * vVector);
        pointsOnArc.add( point );
    }
    pointsOnArc.add( m_endPoint );
}



rg_REAL Arc3D::computeMinDistFromPoint(const rg_Point3D& targetPt, rg_Point3D& minPtOnArc) const
{
    rg_Point3D minPtOnCircle;
    rg_REAL minDistBasedOnCircle = Circle3D::computeMinDistFromPoint(targetPt, minPtOnCircle);

    rg_REAL minDistBasedOnArc = 0.0;


    if ( isOnThisArc(minPtOnCircle) ) {
        minDistBasedOnArc = minDistBasedOnCircle;
        minPtOnArc        = minPtOnCircle;
    }
    else {
        rg_REAL distFromStartToTarget = m_startPoint.distance( targetPt );
        rg_REAL distFromEndToTarget   = m_endPoint.distance( targetPt );

        if ( distFromStartToTarget <= distFromEndToTarget ) {
            minDistBasedOnArc = distFromStartToTarget;
            minPtOnArc        = m_startPoint;
        }
        else {
            minDistBasedOnArc = distFromEndToTarget;
            minPtOnArc        = m_endPoint;
        }
    }

    return minDistBasedOnArc;
}



rg_REAL Arc3D::computeMinDistanceFromLineSegment(const LineSegment3D& lineSegment, 
                                                 rg_Point3D& minPtOnArc, 
                                                 rg_Point3D& minPtOnLineSeg, 
                                                 PosOnLineSegOrArc& pos) const
{
	rg_dList<rg_Point3D> pointsOnArc;
	rg_dList<rg_Point3D> pointsOnLineSegment;
    rg_dList<PosOnLineSegOrArc> positionsOnLineSegment;

	locateLocalMinimaOfDistFromLineSegment(lineSegment, pointsOnArc, pointsOnLineSegment, positionsOnLineSegment);

	rg_REAL minDist = DBL_MAX;
	rg_REAL currDist = -1;
	pointsOnArc.reset4Loop();
	pointsOnLineSegment.reset4Loop();
	positionsOnLineSegment.reset4Loop();
	while(pointsOnArc.setNext4Loop() &&
		  pointsOnLineSegment.setNext4Loop() &&
		  positionsOnLineSegment.setNext4Loop())
	{
		rg_Point3D ptOnArc = pointsOnArc.getEntity();
		rg_Point3D ptOnLineSeg = pointsOnLineSegment.getEntity();
		currDist = ptOnArc.distance(ptOnLineSeg);
		if(rg_LT(currDist, minDist))
		{
			minDist = currDist;
			minPtOnArc = ptOnArc;
			minPtOnLineSeg = ptOnLineSeg;
			pos = positionsOnLineSegment.getEntity();
		}
	}

	return minDist;
}



rg_INT Arc3D::locateLocalMinimaOfDistFromLineSegment(const LineSegment3D& lineSegment, 
											         rg_Point3D*& pointsOnArc, 
											         rg_Point3D*& pointsOnLineSegment, 
											         PosOnLineSegOrArc*& positionsOnLineSegment) const
{
	rg_dList<rg_Point3D>   ptListOnArc,  ptListOnLineSegment;
	rg_dList<PosOnLineSegOrArc> posListOnLineSegment;
	rg_INT numOfPtsOnArcOrLineSeg =
	locateLocalMinimaOfDistFromLineSegment(lineSegment, 
		                                   ptListOnArc, 
										   ptListOnLineSegment, 
										   posListOnLineSegment);

	// Copy into array
	rg_INT num = ptListOnArc.getSize();
	pointsOnArc = new rg_Point3D[num];
	pointsOnLineSegment = new rg_Point3D[num];
	positionsOnLineSegment = new PosOnLineSegOrArc[num];
	
	rg_INDEX index = 0;
	ptListOnArc.reset4Loop();
	ptListOnLineSegment.reset4Loop();
	posListOnLineSegment.reset4Loop();

	while(ptListOnArc.setNext4Loop() &&
		  ptListOnLineSegment.setNext4Loop() &&
		  posListOnLineSegment.setNext4Loop())
	{
		pointsOnArc[index] = ptListOnArc.getEntity();
		pointsOnLineSegment[index] = ptListOnLineSegment.getEntity();
		positionsOnLineSegment[index++] = posListOnLineSegment.getEntity();
	}

	return numOfPtsOnArcOrLineSeg;
}

rg_INT Arc3D::locateLocalMinimaOfDistFromLineSegment(const LineSegment3D& lineSegment, 
											         rg_Point3D*& pointsOnArc, 
											         rg_Point3D*& pointsOnLineSegment, 
											         PosOnLineSegOrArc*& positionsOnArc,
											         PosOnLineSegOrArc*& positionsOnLineSegment) const
{
	rg_dList<rg_Point3D>   ptListOnArc,  ptListOnLineSegment;
	rg_dList<PosOnLineSegOrArc> posListOnArc, posListOnLineSegment;
	rg_INT numOfPtsOnArcOrLineSeg =
	locateLocalMinimaOfDistFromLineSegment(lineSegment, 
		                                   ptListOnArc, 
										   ptListOnLineSegment, 
										   posListOnArc,
										   posListOnLineSegment);

	// Copy into array
	rg_INT num = ptListOnArc.getSize();
	pointsOnArc = new rg_Point3D[num];
	pointsOnLineSegment = new rg_Point3D[num];
	positionsOnArc = new PosOnLineSegOrArc[num];
	positionsOnLineSegment = new PosOnLineSegOrArc[num];
	
	rg_INDEX index = 0;
	ptListOnArc.reset4Loop();
	ptListOnLineSegment.reset4Loop();
	posListOnArc.reset4Loop();
	posListOnLineSegment.reset4Loop();

	while(ptListOnArc.setNext4Loop() &&
		  ptListOnLineSegment.setNext4Loop() &&
		  posListOnArc.setNext4Loop() &&
		  posListOnLineSegment.setNext4Loop())
	{
		pointsOnArc[index] = ptListOnArc.getEntity();
		pointsOnLineSegment[index] = ptListOnLineSegment.getEntity();
		positionsOnArc[index] = posListOnArc.getEntity();
		positionsOnLineSegment[index++] = posListOnLineSegment.getEntity();
	}

	return numOfPtsOnArcOrLineSeg;
}

rg_INT Arc3D::locateLocalMinimaOfDistFromLineSegment(const LineSegment3D& lineSegment, 
													rg_dList<rg_Point3D>& pointsOnArc,
													rg_dList<rg_Point3D>& pointsOnLineSegment,
											        rg_dList<PosOnLineSegOrArc>& positionsOnLineSegment) const
{
	SquaredDistFuncLineSegArc sqdistFucn(lineSegment, *this);

	rg_dList<rg_REAL> localminima;
	rg_INT num = sqdistFucn.computeLocalMinima(localminima);

	localminima.reset4Loop();
	while(localminima.setNext4Loop())
	{
		rg_REAL param = localminima.getEntity();
		rg_Point3D* target = pointsOnLineSegment.add(lineSegment.evaluatePt(param));
		rg_Point3D minPtOnArc;
		computeMinDistFromPoint(*target, minPtOnArc);
		pointsOnArc.add(minPtOnArc);
		if(rg_LT(0.0, param) && rg_LT(param, 1.0))
			positionsOnLineSegment.add(NONEXTREME_PT);
		else if(rg_EQ(param, 0.0))
			positionsOnLineSegment.add(START_PT);
		else
			positionsOnLineSegment.add(END_PT);
	}
	return num;
}

rg_INT Arc3D::locateLocalMinimaOfDistFromLineSegment(const LineSegment3D& lineSegment, 
													 rg_dList<rg_Point3D>& pointsOnArc,
													 rg_dList<rg_Point3D>& pointsOnLineSegment,
													 rg_dList<PosOnLineSegOrArc>& positionsOnArc,
											         rg_dList<PosOnLineSegOrArc>& positionsOnLineSegment) const
{
	SquaredDistFuncLineSegArc sqdistFucn(lineSegment, *this);

	rg_dList<rg_REAL> localminima;
	rg_INT num = sqdistFucn.computeLocalMinima(localminima);

	localminima.reset4Loop();
	while(localminima.setNext4Loop())
	{
		rg_REAL param = localminima.getEntity();
		if(rg_LT(0.0, param) && rg_LT(param, 1.0))
			positionsOnLineSegment.add(NONEXTREME_PT);
		else if(rg_EQ(param, 0.0))
			positionsOnLineSegment.add(START_PT);
		else
			positionsOnLineSegment.add(END_PT);

		rg_Point3D* target = pointsOnLineSegment.add(lineSegment.evaluatePt(param));
		
		rg_Point3D minPtOnArc;
		computeMinDistFromPoint(*target, minPtOnArc);
		if(minPtOnArc == m_startPoint)
			positionsOnArc.add(START_PT);
		else if(minPtOnArc == m_endPoint)
			positionsOnArc.add(END_PT);
		else
			positionsOnArc.add(NONEXTREME_PT);

		pointsOnArc.add(minPtOnArc);
	}
	return num;
}

rg_INT Arc3D::locateInflectionPointsOfDistFromLineSegment(const LineSegment3D& lineSegment, 
		                                                  rg_Point3D*& pointsOnArc,
		                                                  rg_Point3D*& pointsOnLineSegment,
												          PosOnLineSegOrArc*& positionsOnArc,
														  PosOnLineSegOrArc*& positionsOnLineSegment) const
{
	rg_dList<rg_Point3D>   ptListOnArc,  ptListOnLineSegment;
	rg_dList<PosOnLineSegOrArc> posListOnArc, posListOnLineSegment;
	rg_INT numOfPtsOnArcOrLineSeg =
	locateInflectionPointsOfDistFromLineSegment(lineSegment, 
		                                        ptListOnArc, 
										        ptListOnLineSegment, 
										        posListOnArc,
										        posListOnLineSegment);

	// Copy into array
	rg_INT num = ptListOnArc.getSize();
	pointsOnArc = new rg_Point3D[num];
	pointsOnLineSegment = new rg_Point3D[num];
	positionsOnArc = new PosOnLineSegOrArc[num];
	positionsOnLineSegment = new PosOnLineSegOrArc[num];
	
	rg_INDEX index = 0;
	ptListOnArc.reset4Loop();
	ptListOnLineSegment.reset4Loop();
	posListOnArc.reset4Loop();
	posListOnLineSegment.reset4Loop();

	while(ptListOnArc.setNext4Loop() &&
		  ptListOnLineSegment.setNext4Loop() &&
		  posListOnArc.setNext4Loop() &&
		  posListOnLineSegment.setNext4Loop())
	{
		pointsOnArc[index] = ptListOnArc.getEntity();
		pointsOnLineSegment[index] = ptListOnLineSegment.getEntity();
		positionsOnArc[index] = posListOnArc.getEntity();
		positionsOnLineSegment[index++] = posListOnLineSegment.getEntity();
	}

	return numOfPtsOnArcOrLineSeg;
}

rg_INT Arc3D::locateInflectionPointsOfDistFromLineSegment(const LineSegment3D& lineSegment, 
		                                                  rg_dList<rg_Point3D>& pointsOnArc,
												          rg_dList<rg_Point3D>& pointsOnLineSegment,
												          rg_dList<PosOnLineSegOrArc>& positionsOnArc,
												          rg_dList<PosOnLineSegOrArc>& positionsOnLineSegment) const
{
	SquaredDistFuncLineSegArc sqdistFucn(lineSegment, *this);

	rg_dList<rg_REAL> inflectionPts;
	rg_INT num = sqdistFucn.computeInflectionPoints(inflectionPts);

	inflectionPts.reset4Loop();
	while(inflectionPts.setNext4Loop())
	{
		rg_REAL param = inflectionPts.getEntity();
		if(rg_LT(0.0, param) && rg_LT(param, 1.0))
			positionsOnLineSegment.add(NONEXTREME_PT);
		else if(rg_EQ(param, 0.0))
			positionsOnLineSegment.add(START_PT);
		else
			positionsOnLineSegment.add(END_PT);

		rg_Point3D* target = pointsOnLineSegment.add(lineSegment.evaluatePt(param));
		
		rg_Point3D minPtOnArc;
		computeMinDistFromPoint(*target, minPtOnArc);
		if(minPtOnArc == m_startPoint)
			positionsOnArc.add(START_PT);
		else if(minPtOnArc == m_endPoint)
			positionsOnArc.add(END_PT);
		else
			positionsOnArc.add(NONEXTREME_PT);

		pointsOnArc.add(minPtOnArc);
	}
	return num;
}



rg_INT Arc3D::locateLocalMinimaOfDistFromLineSegment(const LineSegment3D&                lineSegment,
                                                     rg_dList<DistEventOnLineSegAndArc>& distEventList )
{
	SquaredDistFuncLineSegArc sqdistFucn(lineSegment, *this);


    ///////////////////////////////////////////////////////////////////////////
    //
    //  locate local min.
    //
	rg_dList<rg_REAL> paramsOnLineSegForLocalMinima;
	sqdistFucn.computeLocalMinima(paramsOnLineSegForLocalMinima);

	paramsOnLineSegForLocalMinima.reset4Loop();
	while( paramsOnLineSegForLocalMinima.setNext4Loop() ) {
		rg_REAL     currParam  = paramsOnLineSegForLocalMinima.getEntity();

        rg_Point3D  pointOnLineSeg = lineSegment.evaluatePt(currParam);
        rg_Point3D  pointOnArc; 
		computeMinDistFromPoint(pointOnLineSeg, pointOnArc);

        PosOnLineSegOrArc positionOnLineSeg = UNKNOWN_POS;
        if ( rg_EQ(currParam, 0.0) || pointOnLineSeg == lineSegment.getStartPt() ) {
			positionOnLineSeg = START_PT;
        }
        else if ( rg_EQ(currParam, 1.0) || pointOnLineSeg == lineSegment.getEndPt() ){
			positionOnLineSeg = END_PT;
        }
        else {
			positionOnLineSeg = NONEXTREME_PT;
        }

        PosOnLineSegOrArc positionOnArc = UNKNOWN_POS;
        if ( pointOnArc == m_startPoint ) {
			positionOnArc = START_PT;
        }
        else if ( pointOnArc == m_endPoint ){
			positionOnArc = END_PT;
        }
        else {
			positionOnArc = NONEXTREME_PT;
        }


        distEventList.add( DistEventOnLineSegAndArc( DET_LOCAL_MIN, 
                                                     pointOnArc, 
                                                     pointOnLineSeg, 
                                                     positionOnArc, 
                                                     positionOnLineSeg ) );

	}
 
    return distEventList.getSize();
}



rg_INT Arc3D::locateMinEventsOnDistanceFromLineSegment(const LineSegment3D&                lineSegment,
                                                       rg_dList<DistEventOnLineSegAndArc>& distEventList )
{
	SquaredDistFuncLineSegArc sqdistFucn(lineSegment, *this);


    ///////////////////////////////////////////////////////////////////////////
    //
    //  locate local min.
    //
	rg_dList<rg_REAL> paramsOnLineSegForLocalMinima;
	sqdistFucn.computeLocalMinima(paramsOnLineSegForLocalMinima);

    rg_BOOL isLocalMinOnMidOfLineSegment = rg_FALSE;
	paramsOnLineSegForLocalMinima.reset4Loop();
	while( paramsOnLineSegForLocalMinima.setNext4Loop() ) {
		rg_REAL     currParam  = paramsOnLineSegForLocalMinima.getEntity();

        rg_Point3D  pointOnLineSeg = lineSegment.evaluatePt(currParam);
        rg_Point3D  pointOnArc; 
		computeMinDistFromPoint(pointOnLineSeg, pointOnArc);

        PosOnLineSegOrArc positionOnLineSeg = UNKNOWN_POS;
        if ( rg_EQ(currParam, 0.0) || pointOnLineSeg == lineSegment.getStartPt() ) {
			positionOnLineSeg = START_PT;
        }
        else if ( rg_EQ(currParam, 1.0) || pointOnLineSeg == lineSegment.getEndPt() ){
			positionOnLineSeg = END_PT;
        }
        else {
			positionOnLineSeg = NONEXTREME_PT;

            isLocalMinOnMidOfLineSegment = rg_TRUE;
        }

        PosOnLineSegOrArc positionOnArc = UNKNOWN_POS;
        if ( pointOnArc == m_startPoint ) {
			positionOnArc = START_PT;
        }
        else if ( pointOnArc == m_endPoint ){
			positionOnArc = END_PT;
        }
        else {
			positionOnArc = NONEXTREME_PT;
        }


        distEventList.add( DistEventOnLineSegAndArc( DET_LOCAL_MIN, 
                                                     pointOnArc, 
                                                     pointOnLineSeg, 
                                                     positionOnArc, 
                                                     positionOnLineSeg ) );

	}
    

    if ( isLocalMinOnMidOfLineSegment ) {
        return distEventList.getSize();
    }


    ///////////////////////////////////////////////////////////////////////////
    //
    //  locate inflection point.
    //
	rg_dList<rg_REAL> paramsOnLineSegForInflectionPts;
	rg_INT numInflectionPts = sqdistFucn.computeInflectionPoints(paramsOnLineSegForInflectionPts);

	paramsOnLineSegForInflectionPts.reset4Loop();
    while( paramsOnLineSegForInflectionPts.setNext4Loop() ) {
		rg_REAL currParam = paramsOnLineSegForInflectionPts.getEntity();

        rg_Point3D  pointOnLineSeg = lineSegment.evaluatePt(currParam);
        rg_Point3D  pointOnArc; 
		computeMinDistFromPoint(pointOnLineSeg, pointOnArc);

        PosOnLineSegOrArc positionOnLineSeg = UNKNOWN_POS;
        if ( rg_EQ(currParam, 0.0) ) {
			positionOnLineSeg = START_PT;
        }
        else if ( rg_EQ(currParam, 1.0) ){
			positionOnLineSeg = END_PT;
        }
        else {
			positionOnLineSeg = NONEXTREME_PT;

            isLocalMinOnMidOfLineSegment = rg_TRUE;
        }

        PosOnLineSegOrArc positionOnArc = UNKNOWN_POS;
        if ( pointOnArc == m_startPoint ) {
			positionOnArc = START_PT;
        }
        else if ( pointOnArc == m_endPoint ){
			positionOnArc = END_PT;
        }
        else {
			positionOnArc = NONEXTREME_PT;
        }


        distEventList.add( DistEventOnLineSegAndArc( DET_INFLECTION_POINT, 
                                                     pointOnArc, 
                                                     pointOnLineSeg, 
                                                     positionOnArc, 
                                                     positionOnLineSeg ) );
	}


    if ( numInflectionPts > 0 ) {
        return distEventList.getSize();
    }



    ///////////////////////////////////////////////////////////////////////////
    //
    //  locate root of 3rd derivative.
    //
	rg_dList<rg_REAL> paramsOnLineSegForRootOfThirdDerivative;
	sqdistFucn.computeRootOfThirdDerivative(paramsOnLineSegForRootOfThirdDerivative);

	paramsOnLineSegForRootOfThirdDerivative.reset4Loop();
    while( paramsOnLineSegForRootOfThirdDerivative.setNext4Loop() ) {
		rg_REAL currParam = paramsOnLineSegForRootOfThirdDerivative.getEntity();

        rg_Point3D  pointOnLineSeg = lineSegment.evaluatePt(currParam);
        rg_Point3D  pointOnArc; 
		computeMinDistFromPoint(pointOnLineSeg, pointOnArc);

        PosOnLineSegOrArc positionOnLineSeg = UNKNOWN_POS;
        if ( rg_EQ(currParam, 0.0) ) {
			positionOnLineSeg = START_PT;
        }
        else if ( rg_EQ(currParam, 1.0) ){
			positionOnLineSeg = END_PT;
        }
        else {
			positionOnLineSeg = NONEXTREME_PT;

            isLocalMinOnMidOfLineSegment = rg_TRUE;
        }

        PosOnLineSegOrArc positionOnArc = UNKNOWN_POS;
        if ( pointOnArc == m_startPoint ) {
			positionOnArc = START_PT;
        }
        else if ( pointOnArc == m_endPoint ){
			positionOnArc = END_PT;
        }
        else {
			positionOnArc = NONEXTREME_PT;
        }


        distEventList.add( DistEventOnLineSegAndArc( DET_ROOT_THIRD_DERIVATIVE, 
                                                     pointOnArc, 
                                                     pointOnLineSeg, 
                                                     positionOnArc, 
                                                     positionOnLineSeg ) );
	}


    return distEventList.getSize();

    /*
    rg_INT numMinDistEvents = 0;

    rg_Point3D*         pointsOnArc            = rg_NULL;
    rg_Point3D*         pointsOnLineSegment    = rg_NULL;
    PosOnLineSegOrArc*  positionsOnArc         = rg_NULL;
    PosOnLineSegOrArc*  positionsOnLineSegment = rg_NULL;

    rg_INT numLocalMin = locateLocalMinimaOfDistFromLineSegment( lineSegment, 
                                                                 pointsOnArc, 
                                                                 pointsOnLineSegment, 
                                                                 positionsOnArc,
                                                                 positionsOnLineSegment);
    numMinDistEvents += numLocalMin;

    rg_BOOL isLocalMinOnMidLineSegment = rg_FALSE;
    for ( rg_INT i=0; i<numLocalMin; i++ ) {
        distEventList.add( DistEventOnLineSegAndArc( DET_LOCAL_MIN, 
                                                     pointsOnArc[i], 
                                                     pointsOnLineSegment[i], 
                                                     positionsOnArc[i], 
                                                     positionsOnLineSegment[i] ) );

        if ( positionsOnLineSegment[i] == NONEXTREME_PT ) {
            isLocalMinOnMidLineSegment = rg_TRUE;
        }
    }

    delete [] positionsOnLineSegment;
    delete [] positionsOnArc;
    delete [] pointsOnLineSegment;
    delete [] pointsOnArc;


    if ( !isLocalMinOnMidLineSegment ) {
        pointsOnArc            = rg_NULL;
        pointsOnLineSegment    = rg_NULL;
        positionsOnArc         = rg_NULL;
        positionsOnLineSegment = rg_NULL;

        rg_INT numInflectionPoints = locateInflectionPointsOfDistFromLineSegment(lineSegment, 
                                                                                 pointsOnArc, 
                                                                                 pointsOnLineSegment, 
                                                                                 positionsOnArc,
                                                                                 positionsOnLineSegment);

        numMinDistEvents += numInflectionPoints;
        for ( i=0; i<numInflectionPoints; i++ ) {
            distEventList.add( DistEventOnLineSegAndArc( DET_INFLECTION_POINT, 
                                                         pointsOnArc[i], 
                                                         pointsOnLineSegment[i], 
                                                         positionsOnArc[i], 
                                                         positionsOnLineSegment[i] ) );
        }

        delete [] positionsOnLineSegment;
        delete [] positionsOnArc;
        delete [] pointsOnLineSegment;
        delete [] pointsOnArc;
    }


    return numMinDistEvents;
    */
}

