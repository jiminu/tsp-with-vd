#include "Circle3D.h"
#include "rg_RelativeOp.h"
#include "rg_TMatrix3D.h"
#include "SquaredDistFuncLineSegCircle.h"

#include <float.h>

#include <utility>
using namespace std;



Circle3D::Circle3D()
: m_radius(0.0)
{
}



Circle3D::Circle3D(const rg_Point3D& center, const rg_REAL& radius, const rg_Point3D& normal)
: m_center(center), m_radius(radius), m_normal(normal)
{
}



Circle3D::Circle3D(const Circle3D& circle)
: m_center(circle.m_center), m_radius(circle.m_radius), m_normal(circle.m_normal)
{
}



Circle3D::~Circle3D()
{
}




//rg_Point3D Circle3D::getCenter() const
//{
//    return m_center;
//}
//
//
//
//rg_REAL    Circle3D::getRadius() const
//{
//    return m_radius;
//}
//
//
//
//rg_Point3D Circle3D::getNormal() const
//{
//    return m_normal;
//}




void       Circle3D::setCenter(const rg_Point3D& center)
{
    m_center = center;
}



void       Circle3D::setRadius(const rg_REAL& radius)
{
    m_radius = radius;
}



void       Circle3D::setNormal(const rg_Point3D& normal)
{
    m_normal = normal;
}



void       Circle3D::setCircle(const rg_Point3D& center, const rg_REAL& radius, const rg_Point3D& normal)
{
    m_center = center;
    m_radius = radius;
    m_normal = normal;
}




Circle3D&  Circle3D::operator =(const Circle3D& circle)
{
    if ( this == &circle ) {
        return *this;
    }

    m_center = circle.m_center;
    m_radius = circle.m_radius;
    m_normal = circle.m_normal;

    return *this;
}



rg_BOOL    Circle3D::operator==(const Circle3D& circle) const
{
    if ( m_center == circle.m_center && rg_EQ(m_radius, circle.m_radius ) ) {
        if ( m_normal == circle.m_normal ) {
            return rg_TRUE;
        }
        else if ( m_normal == -circle.m_normal ) {
            return rg_TRUE;
        }
        else {
            return rg_FALSE;
        }            
    }
    else {
        return rg_FALSE;
    }
}




rg_BOOL    Circle3D::operator!=(const Circle3D& circle) const
{
    if ( operator==(circle) ) {
        return rg_FALSE;
    }
    else {
        return rg_TRUE;
    }
}



rg_BOOL    Circle3D::isOnCircle(const rg_Point3D& point, const rg_REAL tolerance) const
{
    Plane planeContainingThisCircle(m_normal, m_center);

    rg_REAL distFromPlane  = planeContainingThisCircle.distanceFromPoint(point);
    rg_REAL distFromCenter = m_center.distance( point );
        
    if ( rg_ZERO(distFromPlane, tolerance) && rg_EQ(distFromCenter, m_radius, tolerance) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_FLAG Circle3D::isInsideCircle(const rg_Point3D& targetPt, const rg_REAL tolerance) const
{
    Plane planeContainingThisCircle(m_normal, m_center);

    rg_REAL distFromPlane  = planeContainingThisCircle.distanceFromPoint(targetPt);
    rg_REAL distFromCenter = m_center.distance( targetPt );
        
    if ( rg_ZERO(distFromPlane, tolerance) && rg_LE(distFromCenter, m_radius, tolerance) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



void Circle3D::reverse()
{
    m_normal = -m_normal;
}



//Plane      Circle3D::getPlaneContainingThisCircle() const
//{
//    Plane planeContainingThisCircle(m_normal, m_center);
//
//    return planeContainingThisCircle;
//}



rg_INT Circle3D::intersect(const LineSegment3D& lineSegment, rg_Point3D* intersectPts) const
{
    rg_INT     numIntersection = 0;

    rg_REAL    param = 0.0;
    rg_Point3D projPtOnLine = lineSegment.computeProjectedPointForPoint3D(m_center, param);

    if ( projPtOnLine != m_center ) {
        rg_Point3D normal = projPtOnLine - m_center;
        normal.normalize();

        Plane planeContainingLine;
        planeContainingLine.definePlaneByNormalAndPassingPoint(normal, projPtOnLine);

        rg_Point3D intersectionBtwPlaneAndCirle[2];
        rg_INT     numIntersectionBtwPlaneAndCirle = intersect(planeContainingLine, intersectionBtwPlaneAndCirle);

        for ( rg_INT i=0; i<numIntersectionBtwPlaneAndCirle; i++ ) {
            if ( lineSegment.isOnLineSegment(intersectionBtwPlaneAndCirle[i]) ) {
                intersectPts[numIntersection] = intersectionBtwPlaneAndCirle[i];
                numIntersection++;
            }
        }
    }
    else {
        // line passes through the center of circle.
        Plane planeContainingThisCircle(m_normal, m_center);

        rg_Point3D startPt    = lineSegment.getStartPt();
        rg_Point3D endPt      = lineSegment.getEndPt();

        rg_REAL distToStartPt = planeContainingThisCircle.distanceFromPoint( startPt );
        rg_REAL distToEndPt   = planeContainingThisCircle.distanceFromPoint( endPt );

        if ( rg_ZERO(distToStartPt) && rg_ZERO(distToEndPt) ) {
            rg_FLAG isStartPtInCircle = isInsideCircle( startPt );
            rg_FLAG isEndPtInCircle   = isInsideCircle( endPt );

            if ( isStartPtInCircle == rg_FALSE ) {
                rg_Point3D vectorToStart = startPt - m_center;
                vectorToStart.normalize();

                rg_Point3D intersection = m_center + (m_radius*vectorToStart);

                intersectPts[numIntersection] = intersection;
                numIntersection++;
            }

            if ( isEndPtInCircle == rg_FALSE ) {
                rg_Point3D vectorToEnd = endPt - m_center;
                vectorToEnd.normalize();

                rg_Point3D intersection = m_center + (m_radius*vectorToEnd);

                intersectPts[numIntersection] = intersection;
                numIntersection++;
            }
        }
    }

    return numIntersection;
}



//  return: the number of intersections between a circle in 3D and plane 
//          NOTE: If the intersection between a circle in 3D and plane is circle itself,
//                the function return INT_MAX.
rg_INT Circle3D::intersect(const Plane& plane, rg_Point3D* intersectPts, const rg_REAL tolerance) const
{
    rg_INT numIntersections = 0;

	Plane planeOfCircle(m_normal, m_center);

    if ( planeOfCircle == plane ) {
        numIntersections = rg_INT_INFINITY; 
    }
    else if ( planeOfCircle.isParallelTo( plane ) ) {
        numIntersections = 0;
    }
    else {
        // Two planes intersect.

        Line3D line;
        planeOfCircle.intersect( plane, line );
    
	    rg_REAL distToLine = line.computeDistanceFromPoint3D(m_center);

        if ( rg_EQ(distToLine, m_radius, tolerance) ) {
	        // There is an intersection point between the circle and a plane.
	        rg_Point3D ptOnLine = line.computeProjectedPointForPoint3D(m_center);

		    intersectPts[ 0 ] = ptOnLine;

		    numIntersections = 1;
        } 
        else if ( rg_LT(distToLine, m_radius, tolerance) ) {
	        // There are two intersection points between the circle and a plane.

	        rg_Point3D ptOnLine = line.computeProjectedPointForPoint3D(m_center);
	        rg_Point3D transVec = line.getDirVec().getUnitVector();

	        rg_REAL    distToIntersecPts = sqrt( (m_radius*m_radius) - (distToLine*distToLine) );

		    intersectPts[ 0 ] = ptOnLine + (transVec * distToIntersecPts);
		    intersectPts[ 1 ] = ptOnLine - (transVec * distToIntersecPts);

		    numIntersections = 2;
        }
        else {
	        // There is no intersection between the circle and a plane.
        }
    }

    return numIntersections;
}



rg_INT Circle3D::intersect(const Circle3D& circle, rg_Point3D* intersectPts, const rg_REAL tolerance) const
{
    if ( (*this) == circle ) {
        return rg_INT_INFINITY;
    }
    

    rg_INT numIntersection = 0;
    Plane  planeOfThisCircle( m_normal,        m_center);
    Plane  planeOfInputCircle(circle.m_normal, circle.m_center);


    if ( planeOfThisCircle != planeOfInputCircle ) {
        rg_Point3D intersection[2];
        rg_INT     numIntersectionWithPlane = intersect(planeOfInputCircle, intersection, tolerance);

        for ( rg_INT i=0; i<numIntersectionWithPlane; i++ ) {
            if ( circle.isOnCircle( intersection[i] ) ) {
                intersectPts[numIntersection] = intersection[i];
                numIntersection++;
            }
        }
    }
    else {
        //  It means that two circles are on the same plane.
        rg_REAL distBtwCenters = m_center.distance(circle.m_center);
        rg_REAL sumOfTwoRadii  = m_radius + circle.m_radius;

        //  disc_A : disc of this circle
        //  disc_B : disc of input circle

        //  The closed disc_A and disc_B do not intersect.
        if ( rg_GT(distBtwCenters, sumOfTwoRadii, tolerance) ) {
            numIntersection = 0;
        }
        //  The closed disc_A and disc_B intersect at a single point in their outside.
        else if ( rg_EQ(distBtwCenters, sumOfTwoRadii, tolerance) ) {
            rg_Point3D vectorToCenter = circle.m_center - m_center;
            vectorToCenter.normalize();

            intersectPts[0] = m_center + (m_radius*vectorToCenter);
            numIntersection = 1;
        }
        //  The closed disc_A and disc_B intersect.
        else {
            //  The disc_B is completely contained in the disc_A.
            if (      rg_GT( m_radius, (circle.m_radius+distBtwCenters), tolerance) ) {
                numIntersection = 0;
            }
            //  The disc_A is completely contained in the disc_B.
            else if ( rg_LT( (m_radius+distBtwCenters), circle.m_radius, tolerance) ) {
                numIntersection = 0;
            }
            //  The disc_B is contained in the disc_A.
            //  And the closed disc_A and disc_B intersect at a single point in the inside of the disc_A.
            else if ( rg_EQ( m_radius, (circle.m_radius+distBtwCenters), tolerance) ) {
                rg_Point3D vectorToCenter = circle.m_center - m_center;
                vectorToCenter.normalize();

                intersectPts[0] = m_center + (m_radius*vectorToCenter);
                numIntersection = 1;
            }
            //  The disc_A is contained in the disc_B.
            //  And the closed disc_A and disc_B intersect at a single point in the inside of the disc_B.
            else if ( rg_EQ( (m_radius+distBtwCenters), circle.m_radius, tolerance) ) {
                rg_Point3D vectorToCenter = m_center - circle.m_center;
                vectorToCenter.normalize();

                intersectPts[0] = m_center + (m_radius*vectorToCenter);
                numIntersection = 1;
            }
            else {
                if ( rg_EQ(m_radius, circle.m_radius, tolerance) ) {
                    rg_Point3D vectorToCenter = circle.m_center - m_center;
                    rg_Point3D midPoint       = (m_center + circle.m_center)*0.5;

                    rg_REAL    distToMidPoint = m_center.distance(midPoint);

                    rg_Point3D vectorFromMidPtToIntersection = m_normal.crossProduct( vectorToCenter );
                    vectorFromMidPtToIntersection.normalize();

                    rg_REAL    distFromMidPtToIntersection   = sqrt( (m_radius*m_radius) - (distToMidPoint*distToMidPoint) );

                    intersectPts[0] = midPoint + (distFromMidPtToIntersection*vectorFromMidPtToIntersection);
                    intersectPts[1] = midPoint - (distFromMidPtToIntersection*vectorFromMidPtToIntersection);

                    numIntersection = 2;
                }
                else {
                    rg_Point3D vectorToCenter = circle.m_center - m_center;
                    vectorToCenter.normalize();

                    rg_REAL a = m_radius;
                    rg_REAL b = circle.m_radius;
                    rg_REAL c = distBtwCenters;

                    rg_REAL cosBeta = ( a*a + c*c - b*b )/(2.0*a*c);

                    rg_REAL    distToMidPoint = a*cosBeta;
                    rg_Point3D midPoint       = m_center + (distToMidPoint*vectorToCenter);

                    rg_Point3D vectorFromMidPtToIntersection = m_normal.crossProduct( vectorToCenter );
                    vectorFromMidPtToIntersection.normalize();

                    rg_REAL    distFromMidPtToIntersection   = sqrt( (m_radius*m_radius) - (distToMidPoint*distToMidPoint) );

                    intersectPts[0] = midPoint + (distFromMidPtToIntersection*vectorFromMidPtToIntersection);
                    intersectPts[1] = midPoint - (distFromMidPtToIntersection*vectorFromMidPtToIntersection);

                    numIntersection = 2;
                }
            }
        }
    }

    return numIntersection;
}



rg_REAL Circle3D::evaluateCircumference() const
{
    return (2.0*m_radius*rg_PI);
}



rg_REAL Circle3D::evaluateAngle(const rg_Point3D& pt1OnCircle, const rg_Point3D& pt2OnCircle) const
{
    Plane planeContainingThisCircle(m_normal, m_center);

    rg_REAL angle = planeContainingThisCircle.computeAngleInCCW( pt1OnCircle, m_center, pt2OnCircle );

    return angle;
}


rg_Point3D Circle3D::evaluateMiddlePoint(const rg_Point3D& point1, const rg_Point3D& point2) const
{
    rg_REAL angle = evaluateAngle(point1, point2);
    rg_Point3D angularBisector = point1 + point2;
    angularBisector.normalize();

    rg_Point3D middlePoint;
    if ( angle <= rg_PI ) {
        middlePoint = m_center + (m_radius * angularBisector);
    }
    else {
        middlePoint = m_center - (m_radius * angularBisector);
    }

    return middlePoint;
}



void Circle3D::evaluatePointsOnCircle(const rg_REAL& unitLength, const rg_Point3D& startPoint,
                                   rg_dList<rg_Point3D>& pointsOnCircle) const
{
    rg_INT  numLineSegment = (rg_INT) ceil( evaluateCircumference()/unitLength );
    if ( numLineSegment < 4 ) {
        numLineSegment = 4;
    }

    rg_REAL stepAngle = 2.0*rg_PI/numLineSegment;
    rg_REAL magitudeOfTangent = m_radius*tan( stepAngle );

    pointsOnCircle.add( startPoint );
    rg_Point3D currPoint = startPoint;
    for ( rg_INT i=0; i<(numLineSegment-1); i++ ) {
        rg_Point3D tangent = m_normal.crossProduct( currPoint-m_center );
        tangent.normalize();

        rg_Point3D sector = (currPoint + (magitudeOfTangent*tangent)) - m_center;
        sector.normalize();

        currPoint = m_center + (m_radius*sector);
        pointsOnCircle.add( currPoint );
    }
    pointsOnCircle.add( startPoint );
}



void    Circle3D::evaluatePointsOnCircle(list<rg_Point3D>& pointsOnCircle, const int& numPoints) const
{
    rg_REAL stepAngle = 2.0*rg_PI / numPoints;
    rg_REAL magitudeOfTangent = m_radius*tan(stepAngle);

    rg_Point3D startPoint;
    rg_Point3D normal = m_normal.getUnitVector();

    if (normal == rg_Point3D(0.0, 0.0, 1.0)) {
        startPoint.setPoint(m_radius, 0.0, 0.0);
    }
    else if (normal == rg_Point3D(0.0, 0.0, -1.0)) {
        startPoint.setPoint(-m_radius, 0.0, 0.0);
    }
    else {
        rg_TMatrix3D rotateMat;
        rotateMat.rotate(rg_Point3D(0.0, 0.0, 1.0), normal );
        startPoint = rotateMat*rg_Point3D(m_radius, 0.0, 0.0);
    }
    startPoint = m_center + startPoint;


    pointsOnCircle.push_back(startPoint);
    rg_Point3D currPoint = startPoint;
    for (rg_INT i = 0; i<(numPoints - 1); i++) {
        rg_Point3D tangent = m_normal.crossProduct(currPoint - m_center);
        tangent.normalize();

        rg_Point3D sector = (currPoint + (magitudeOfTangent*tangent)) - m_center;
        sector.normalize();

        currPoint = m_center + (m_radius*sector);
        pointsOnCircle.push_back(currPoint);
    }
}



rg_REAL Circle3D::computeMinDistanceFromLineSegment(const LineSegment3D& lineSegment, rg_Point3D& minPtOnCircle, rg_Point3D& minPtOnLineSeg, PosOnLineSegOrArc& posOnLineSeg) const
{
	// The Plane of circle
	rg_Point3D    unitNormalVec = m_normal.getUnitVector();
	Plane         planeOfCircle(unitNormalVec, m_center);

	// Project line segment onto the plane of circle
	LineSegment3D projectionOfLineSeg;
	planeOfCircle.projectLineSegmentOnPlane(lineSegment, projectionOfLineSeg);
	rg_Point3D sPtOnProjectionOfLS = projectionOfLineSeg.getStartPt();
	rg_Point3D ePtOnProjectionOfLS = projectionOfLineSeg.getEndPt();

	// Check whether line segment intersects the plane
	rg_Point3D   intersectionWithPlane;
	rg_REAL      paramOfPlaneIntersection;
	PosOnLineSegOrArc posOnLSForPlaneIntersection;
	rg_FLAG      bIntersectPlane = planeOfCircle.computeIntersectionWithLineSegment(lineSegment, 
		                                                                            intersectionWithPlane, 
																			        paramOfPlaneIntersection, 
																			        posOnLSForPlaneIntersection);

	rg_REAL minDist = DBL_MAX;
	rg_dList< pair<rg_Point3D, rg_REAL> > localMinimaMaxima;

	// 1. Line segment is perpendicular to the plane of the circle
	//    and intersects the plane
	if(sPtOnProjectionOfLS == ePtOnProjectionOfLS && bIntersectPlane)
	{
		// Plane intersection
		minPtOnLineSeg = intersectionWithPlane;
		minDist = computeMinDistFromPoint(minPtOnLineSeg, minPtOnCircle);
		if(rg_ZERO(paramOfPlaneIntersection))
			posOnLineSeg = START_PT;
		else if(rg_EQ(paramOfPlaneIntersection, 1.0))
			posOnLineSeg = END_PT;
		else
			posOnLineSeg = NONEXTREME_PT;
		return minDist;
	}
	// 2. Line segment is perpendicular to the plane of the circle
	//    but does not intersect the plane
	else if(sPtOnProjectionOfLS == ePtOnProjectionOfLS && !bIntersectPlane)
	{
		// Min. dist. point from circle center
		rg_REAL paramOfMinDistPt;
		lineSegment.computeMinDistFromPoint(m_center, minPtOnLineSeg, paramOfMinDistPt, posOnLineSeg);
		minDist = computeMinDistFromPoint(minPtOnLineSeg, minPtOnCircle);
		return minDist;
	}

	// 3. Line segment is not perpendicular to the plane of the circle
	rg_Point3D sPtOnLS = lineSegment.getStartPt();
	rg_Point3D ePtOnLS = lineSegment.getEndPt();	
	// Two end points of line segment
	pair<rg_Point3D, rg_REAL> slocalMinMax(sPtOnLS, 0.0);
	pair<rg_Point3D, rg_REAL> elocalMinMax(ePtOnLS, 1.0);
	localMinimaMaxima.add(slocalMinMax);
	localMinimaMaxima.add(elocalMinMax);
	// Min. dist. point from circle center
	rg_Point3D   minPt;
	rg_REAL      paramOfMinDistPt;
	PosOnLineSegOrArc posOnLSForMinDistPt;
	lineSegment.computeMinDistFromPoint(m_center, minPt, paramOfMinDistPt, posOnLSForMinDistPt);
	if(posOnLSForMinDistPt == NONEXTREME_PT && !rg_EQ(paramOfMinDistPt, paramOfPlaneIntersection))
	{
		pair<rg_Point3D, rg_REAL> localMinMax(minPt, paramOfMinDistPt);
		localMinimaMaxima.add(localMinMax);
	}
	// Plane intersection
	if(bIntersectPlane && posOnLSForPlaneIntersection == NONEXTREME_PT)
	{
		pair<rg_Point3D, rg_REAL> localMinMax(intersectionWithPlane, paramOfPlaneIntersection);
		localMinimaMaxima.add(localMinMax);
	}
	// Check whether the projection of line segment intersects circle
	rg_Point3D  intersectionWithCircle[ 2 ];
	rg_REAL     paramOfCircleIntersection[ 2 ];
	rg_INT      numOfIntersections = -1;
	numOfIntersections = computeIntersectionWith(projectionOfLineSeg, 
		                                         intersectionWithCircle, 
												 paramOfCircleIntersection);
	// Projection of line segment has intersection with circle
	if(numOfIntersections > 0)
	{
		// circle intersection
		rg_Point3D dirVec = lineSegment.getDirVec();		
		for(int i = 0;i < numOfIntersections;i++)
		{
			if(rg_EQ(paramOfCircleIntersection[ i ], 0.0) ||
			   rg_EQ(paramOfCircleIntersection[ i ], 1.0) ||
			   rg_EQ(paramOfCircleIntersection[ i ], paramOfMinDistPt) ||
			   rg_EQ(paramOfCircleIntersection[ i ], paramOfPlaneIntersection) )
			   continue;

			rg_Point3D liftUpOfIntersection = sPtOnLS + dirVec * paramOfCircleIntersection[ i ];
			pair<rg_Point3D, rg_REAL> localMinMax(liftUpOfIntersection, paramOfCircleIntersection[ i ]);
			localMinimaMaxima.add(localMinMax);
		}
	}

	// Sort the local minima and maxima
	// Note that the number of local minima and maxima is less than or equal to '6'.
	// (1) start and end points of line segment                          : 2
	// (2) intersection between line segment and plane                   : 1
	// (3) intersection between circle and projection of line segment    : 2
	// (4) minimum distance point between line segment and circle center : 1
	pair<rg_Point3D, rg_REAL>* currLocalMinMax = rg_NULL;
	rg_REAL currDist = -1.0;
	localMinimaMaxima.reset4Loop();
	while(localMinimaMaxima.setNext4Loop())
	{
		currLocalMinMax = localMinimaMaxima.getpEntity();
		rg_Point3D lcoalMinMaxPtOnLineSeg = currLocalMinMax->first;
		rg_Point3D localMinMaxPtOnCircle;
		currDist = computeMinDistFromPoint(lcoalMinMaxPtOnLineSeg, localMinMaxPtOnCircle);
		if(currDist < minDist)
		{
			minDist = currDist;
			minPtOnLineSeg = lcoalMinMaxPtOnLineSeg;
			minPtOnCircle  = localMinMaxPtOnCircle;
			if(rg_ZERO(currLocalMinMax->second))
				posOnLineSeg = START_PT;
			else if(rg_EQ(currLocalMinMax->second, 1.0))
				posOnLineSeg = END_PT;
			else
				posOnLineSeg = NONEXTREME_PT;
		}
	}
	return minDist;
}



rg_REAL Circle3D::computeMinDistFromPoint(const rg_Point3D& targetPt, rg_Point3D& minPtOnCircle) const
{
	// Project the target point onto the plane of circle
	Plane      planeOfCircle(m_normal, m_center);
	rg_Point3D projectionOfTarget = planeOfCircle.projectPointOnPlane(targetPt);

	if(projectionOfTarget != m_center) {
	    rg_Point3D dirVec = (projectionOfTarget - m_center).getUnitVector();
		minPtOnCircle     = m_center + m_radius * dirVec;

		return minPtOnCircle.distance(targetPt);
	}
	else {
		minPtOnCircle = m_center;
		rg_REAL hDist = m_radius;
		rg_REAL vDist = projectionOfTarget.distance(targetPt);
		return sqrt(hDist * hDist + vDist * vDist);
	}	
}



rg_INT Circle3D::computeIntersectionWith(const LineSegment3D& linesegment, rg_Point3D intersectPts[], rg_REAL param[]) const
{
	rg_REAL coeff[ 3 ] = {0., 0., 0.};
	rg_Point3D dirVec = linesegment.getDirVec();
	rg_Point3D sPt = linesegment.getStartPt();

	coeff[ 2 ] = dirVec.squaredMagnitude();
	coeff[ 1 ] = 2*sPt.innerProduct(dirVec) -2*m_center.innerProduct(dirVec);
	coeff[ 0 ] = sPt.squaredMagnitude() + m_center.squaredMagnitude() - 2*sPt.innerProduct(m_center) - m_radius * m_radius; 

	rg_REAL determinant = coeff[ 1 ] * coeff[ 1 ] - 4 * coeff[ 2 ] * coeff[ 0 ];
	if(rg_LT(determinant, 0.0))
		return 0;

	param[ 0 ] = (-coeff[ 1 ] - sqrt(determinant)) / (2 * coeff[ 2 ]);
	param[ 1 ] = (-coeff[ 1 ] + sqrt(determinant)) / (2 * coeff[ 2 ]);
	rg_INDEX numOfIntersections = 0;
	rg_FLAG  bSolution[ 2 ] = {rg_FALSE, rg_FALSE};
	if(rg_GE(param[ 0 ], 0.0) && rg_LE(param[ 0 ], 1.0))
	{
		intersectPts[ 0 ] = sPt + param[ 0 ] * dirVec;
		numOfIntersections++;
		bSolution[ 0 ] = rg_TRUE;
	}
	if(rg_GE(param[ 1 ], 0.0) && rg_LE(param[ 1 ], 1.0))
	{
		intersectPts[ 1 ] = sPt + param[ 1 ] * dirVec;
		numOfIntersections++;
		bSolution[ 1 ] = rg_TRUE;
	}

	if(numOfIntersections == 1)
	{
		if(bSolution[ 1 ])
		{
			param[ 0 ] = param[ 1 ];
			intersectPts[ 0 ] = intersectPts[ 1 ];
		}
	}


	return numOfIntersections;
}




rg_INT Circle3D::locatePointsWithLocalMinimaFromLineSegment(const LineSegment3D& lineSegment,
															rg_dList<rg_Point3D>& pointsOnCircle,
															rg_dList<rg_Point3D>& pointsOnLineSegment,
													        rg_dList<rg_REAL>&    paramsOnLineSegment) const
{
	SquaredDistFuncLineSegCircle sqdDistFunc(lineSegment, *this);
	sqdDistFunc.computeLocalMinima(paramsOnLineSegment);
	rg_INT numOfLocalMinima = paramsOnLineSegment.getSize();
	paramsOnLineSegment.reset4Loop();
	while(paramsOnLineSegment.setNext4Loop())
	{
		rg_Point3D* ptOnLineSeg
			= pointsOnLineSegment.add(lineSegment.evaluatePt(paramsOnLineSegment.getEntity()));
		rg_Point3D  ptOnCircle;
		computeMinDistFromPoint(*ptOnLineSeg, ptOnCircle);
		pointsOnCircle.add(ptOnCircle);
	}
	return numOfLocalMinima;
}

rg_INT Circle3D::locatePointsWithLocalMinimaFromLineSegment(const LineSegment3D& lineSegment,
															rg_dList<rg_Point3D>& pointsOnCircle,
															rg_dList<rg_Point3D>& pointsOnLineSegment,
															rg_dList<rg_REAL>&    paramsOnLineSegment,
															SquaredDistFuncLineSegCircle* sqdDistFunc) const
{
	sqdDistFunc->computeLocalMinima(paramsOnLineSegment);
	rg_INT numOfLocalMinima = paramsOnLineSegment.getSize();
	paramsOnLineSegment.reset4Loop();
	while(paramsOnLineSegment.setNext4Loop())
	{
		rg_Point3D* ptOnLineSeg 
			= pointsOnLineSegment.add(lineSegment.evaluatePt(paramsOnLineSegment.getEntity()));
		rg_Point3D  ptOnCircle;
		computeMinDistFromPoint(*ptOnLineSeg, ptOnCircle);
		pointsOnCircle.add(ptOnCircle);
	}
	return numOfLocalMinima;
}

rg_INT Circle3D::locatePointsWithLocalMinimaFromLineSegment(const LineSegment3D& lineSegment, 
															rg_dList<rg_Point3D>& pointsOnCircle, 
															rg_dList<rg_Point3D>& pointsOnLineSegment, 
															rg_dList<PosOnLineSegOrArc>& positionsOnCircle, 
															rg_dList<PosOnLineSegOrArc>& positionsOnLineSegment, 
															rg_dList<rg_REAL>& paramsOnLineSegment, 
															SquaredDistFuncLineSegCircle* sqdDistFunc) const
{
	sqdDistFunc->computeLocalMinima(paramsOnLineSegment);
	rg_INT numOfLocalMinima = paramsOnLineSegment.getSize();
	paramsOnLineSegment.reset4Loop();
	while(paramsOnLineSegment.setNext4Loop())
	{
		rg_REAL param = paramsOnLineSegment.getEntity();
		if(rg_EQ(param, 0.0))
		{
			positionsOnLineSegment.add(START_PT);
			rg_Point3D ptOnLineSeg = lineSegment.getStartPt();
			pointsOnLineSegment.add(ptOnLineSeg);
			positionsOnCircle.add(NONEXTREME_PT);
			rg_Point3D ptOnCircle;
			computeMinDistFromPoint(ptOnLineSeg, ptOnCircle);			
			pointsOnCircle.add(ptOnCircle);
		}
		else if(rg_EQ(param, 1.0))
		{
			positionsOnLineSegment.add(END_PT);
			rg_Point3D ptOnLineSeg = lineSegment.getEndPt();
			pointsOnLineSegment.add(ptOnLineSeg);
			positionsOnCircle.add(NONEXTREME_PT);
			rg_Point3D ptOnCircle;
			computeMinDistFromPoint(ptOnLineSeg, ptOnCircle);			
			pointsOnCircle.add(ptOnCircle);
		}
		else
		{
			positionsOnLineSegment.add(NONEXTREME_PT);
			rg_Point3D ptOnLineSeg = lineSegment.evaluatePt(paramsOnLineSegment.getEntity());
			pointsOnLineSegment.add(ptOnLineSeg);
			positionsOnCircle.add(NONEXTREME_PT);
			rg_Point3D ptOnCircle;
			computeMinDistFromPoint(ptOnLineSeg, ptOnCircle);
			pointsOnCircle.add(ptOnCircle);
		}
	}
	return numOfLocalMinima;
}

rg_INT Circle3D::locatePointsWithLocalMinimaFromLineSegment(const LineSegment3D& lineSegment, 
												  rg_Point3D*& pointsOnCircle, 
												  rg_Point3D*& pointsOnLineSegment, 
												  rg_REAL*&    paramsOnLineSegment,
												  PosOnLineSegOrArc*& positionsOnLineSegment,
												  rg_REAL*&    minDistances) const
{
	rg_INT numOfLocalMinimaMaxima = locatePointsWithLocalMinimaMaximaFromLineSegment(lineSegment,
		                                                                   pointsOnCircle,
																		   pointsOnLineSegment,
																		   paramsOnLineSegment,
																		   minDistances);

	// sort local minima and maxima according to their parameter on linesegment
	// (selection sort)
	rg_INT i = 0;
	for(i = 0;i < numOfLocalMinimaMaxima;i++)
	{
		rg_INT minIndex = i;
		rg_INT indexBnd = i + 1;
		for(rg_INT j = indexBnd;j < numOfLocalMinimaMaxima;j++)
		{
			if(rg_LT(paramsOnLineSegment[ j ], paramsOnLineSegment[minIndex]))
				minIndex = j;
		}
		// swap
		if(minIndex != i)
		{
			rg_Point3D tPt = pointsOnCircle[ i ];
			pointsOnCircle[ i ] = pointsOnCircle[minIndex];
			pointsOnCircle[minIndex] = tPt;
			tPt = pointsOnLineSegment[ i ];
			pointsOnLineSegment[ i ] = pointsOnLineSegment[minIndex];
			pointsOnLineSegment[minIndex] = tPt;
			rg_REAL tParam = paramsOnLineSegment[ i ];
			paramsOnLineSegment[ i ] = paramsOnLineSegment[minIndex];
			paramsOnLineSegment[minIndex] = tParam;
			rg_REAL tDist = minDistances[ i ];
			minDistances[ i ] = minDistances[minIndex];
			minDistances[minIndex] = tDist;
		}
	}

	rg_dList<rg_Point3D> ptsOnCircle, ptsOnLineSeg;
	rg_dList<rg_REAL>    paramsOnLineSeg, distances;

	// get only local minima
	if(rg_GT(minDistances[ 1 ], minDistances[ 0 ]))
	{
		ptsOnCircle.add(pointsOnCircle[ 0 ]);
		ptsOnLineSeg.add(pointsOnLineSegment[ 0 ]);
		paramsOnLineSeg.add(paramsOnLineSegment[ 0 ]);
		distances.add(minDistances[ 0 ]);
	}
	if(rg_GT(minDistances[numOfLocalMinimaMaxima - 2], minDistances[numOfLocalMinimaMaxima - 1]))
	{
		ptsOnCircle.add(pointsOnCircle[numOfLocalMinimaMaxima - 1]);
		ptsOnLineSeg.add(pointsOnLineSegment[numOfLocalMinimaMaxima - 1]);
		paramsOnLineSeg.add(paramsOnLineSegment[numOfLocalMinimaMaxima - 1]);
		distances.add(minDistances[numOfLocalMinimaMaxima - 1]);
	}

	rg_INT indexBnd = numOfLocalMinimaMaxima - 1;
	for(i = 1;i < indexBnd;i++)
	{
		if(rg_LT(minDistances[ i ], minDistances[i - 1]) && 
		   rg_LT(minDistances[ i ], minDistances[i + 1]))
		{
			ptsOnCircle.add(pointsOnCircle[ i ]);
			ptsOnLineSeg.add(pointsOnLineSegment[ i ]);
			paramsOnLineSeg.add(paramsOnLineSegment[ i ]);
			distances.add(minDistances[ i ]);
		}
	}
	delete [] pointsOnCircle;
	delete [] pointsOnLineSegment;
	delete [] paramsOnLineSegment;
	
	rg_INT numOfLocalMinima = paramsOnLineSeg.getSize();
	pointsOnCircle = new rg_Point3D[numOfLocalMinima];
	pointsOnLineSegment = new rg_Point3D[numOfLocalMinima];
	paramsOnLineSegment = new rg_REAL[numOfLocalMinima];
	positionsOnLineSegment = new PosOnLineSegOrArc[numOfLocalMinima];

	i = 0;
	ptsOnCircle.reset4Loop();
	while(ptsOnCircle.setNext4Loop())
		pointsOnCircle[ i++ ] = ptsOnCircle.getEntity();
	i = 0;
	ptsOnLineSeg.reset4Loop();
	while(ptsOnLineSeg.setNext4Loop())
		pointsOnLineSegment[ i++ ] = ptsOnLineSeg.getEntity();
	i = 0;
	paramsOnLineSeg.reset4Loop();
	while(paramsOnLineSeg.setNext4Loop())
	{
		paramsOnLineSegment[ i++ ] = paramsOnLineSeg.getEntity();
		if(rg_EQ(paramsOnLineSegment[ i ], 0.0))
			positionsOnLineSegment[ i ] = START_PT;
		else if(rg_EQ(paramsOnLineSegment[ i ], 1.0))
			positionsOnLineSegment[ i ] = END_PT;
		else
			positionsOnLineSegment[ i ] = NONEXTREME_PT;
	}
	return numOfLocalMinima;
}

// void Circle3D::makeQuarticPolyForLocalMiniMaxDistBTWLinesegNCircle(const LineSegment3D& linesegment, rg_QuarticPolynomial& QTPolynomial) const
// {
// 	// make distance function
// 	rg_REAL coeffOfDistFunc[ 6 ] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
// 	rg_Point3D sPt = linesegment.getStartPt();
// 	rg_Point3D ePt = linesegment.getEndPt();
// 	rg_REAL    dist   = - m_normalVec.innerProduct(m_center);
// 	
// 	rg_Point3D sPtMinusEPt = (sPt - ePt);
// 	rg_REAL    sqdMgntudeOfSPtMinusEPt = sPtMinusEPt.squaredMagnitude();
// 	rg_Point3D centerMinusSPt = (m_center - sPt);
// 	rg_REAL    sqdMgntudeOfCenterMinusSPt = centerMinusSPt.squaredMagnitude();
// 	rg_REAL    normalInnerProdSPt = m_normalVec.innerProduct(sPt);
// 	rg_REAL    normalInnerProdSPtPlusDist = normalInnerProdSPt + dist;
// 	rg_REAL    normalInnerProdSPtMinusEPt = m_normalVec.innerProduct(sPtMinusEPt);
// 	rg_REAL    sqdRadius = m_radius * m_radius;
// 	
// 	coeffOfDistFunc[ 5 ] = sqdMgntudeOfSPtMinusEPt;
// 	coeffOfDistFunc[ 4 ] = 2.0 * centerMinusSPt.innerProduct(sPtMinusEPt);
// 	coeffOfDistFunc[ 3 ] = sqdMgntudeOfCenterMinusSPt + sqdRadius;
// 	coeffOfDistFunc[ 2 ] = sqdMgntudeOfSPtMinusEPt - normalInnerProdSPtMinusEPt * normalInnerProdSPtMinusEPt;
// 	coeffOfDistFunc[ 1 ] = 2.0 * centerMinusSPt.innerProduct(sPtMinusEPt) + 2.0 * normalInnerProdSPtMinusEPt * normalInnerProdSPtPlusDist;
// 	coeffOfDistFunc[ 0 ] = sqdMgntudeOfCenterMinusSPt - normalInnerProdSPtPlusDist * normalInnerProdSPtPlusDist;
// 	
// 	// make quartic polynomial by equating the derivative of distance function to zero
// 	rg_REAL sqdCoeffOfDistFunc5 = coeffOfDistFunc[ 5 ] * coeffOfDistFunc[ 5 ];
// 	rg_REAL sqdCoeffOfDistFunc4 = coeffOfDistFunc[ 4 ] * coeffOfDistFunc[ 4 ];
// 	
// 	rg_REAL coeffOfQP[ 5 ] = {0.0, 0.0, 0.0, 0.0, 0.0};
// 	coeffOfQP[ 4 ] = 4.0 * sqdCoeffOfDistFunc5 * coeffOfDistFunc[ 2 ];
// 	coeffOfQP[ 3 ] = 4.0 * coeffOfDistFunc[ 5 ] * (coeffOfDistFunc[ 5 ] * coeffOfDistFunc[ 1 ] + coeffOfDistFunc[ 4 ] * coeffOfDistFunc[ 2 ]);
// 	coeffOfQP[ 2 ] = 4.0 * coeffOfDistFunc[ 5 ] * (coeffOfDistFunc[ 5 ] * coeffOfDistFunc[ 0 ] + coeffOfDistFunc[ 4 ] * coeffOfDistFunc[ 1 ]) +
// 		coeffOfDistFunc[ 2 ] * (sqdCoeffOfDistFunc4 - 4.0 * sqdRadius * coeffOfDistFunc[ 2 ]);
// 	coeffOfQP[ 1 ] = 4.0 * coeffOfDistFunc[ 5 ] * coeffOfDistFunc[ 4 ] * coeffOfDistFunc[ 0 ] +
// 		coeffOfDistFunc[ 1 ] * (sqdCoeffOfDistFunc4 - 4.0 * sqdRadius * coeffOfDistFunc[ 2 ]);
// 	coeffOfQP[ 0 ] = sqdCoeffOfDistFunc4 * coeffOfDistFunc[ 0 ] - sqdRadius * coeffOfDistFunc[ 1 ] * coeffOfDistFunc[ 1 ];
// 	QTPolynomial.setCoefficient(coeffOfQP);
// }

rg_INT Circle3D::locatePointsWithLocalMinimaMaximaFromLineSegment(const LineSegment3D& lineSegment, 
												        rg_Point3D*& pointsOnCircle, 
												        rg_Point3D*& pointsOnLineSegment, 
												        rg_REAL*&    paramsOnLineSegment,
														rg_REAL*&    minDistances) const
{
	rg_INT maxNumOfLocalMinimaMaxima = 6;
	pointsOnCircle = new rg_Point3D[maxNumOfLocalMinimaMaxima];
	pointsOnLineSegment = new rg_Point3D[maxNumOfLocalMinimaMaxima];
	paramsOnLineSegment = new rg_REAL[maxNumOfLocalMinimaMaxima];
	minDistances = new rg_REAL[maxNumOfLocalMinimaMaxima];
	rg_INT indexOfPtArray = 0;

	// The Plane of circle
	rg_Point3D    unitNormalVec = m_normal.getUnitVector();
	Plane         planeOfCircle(unitNormalVec, m_center);
	
	// Project line segment onto the plane of circle
	LineSegment3D projectionOfLineSeg;
	planeOfCircle.projectLineSegmentOnPlane(lineSegment, projectionOfLineSeg);
	rg_Point3D sPtOnProjectionOfLS = projectionOfLineSeg.getStartPt();
	rg_Point3D ePtOnProjectionOfLS = projectionOfLineSeg.getEndPt();
	
	// Check whether line segment intersects the plane
	rg_Point3D   intersectionWithPlane;
	rg_REAL      paramOfPlaneIntersection;
	PosOnLineSegOrArc posOnLSForPlaneIntersection;
	rg_FLAG      bIntersectPlane 
		= planeOfCircle.computeIntersectionWithLineSegment(lineSegment, 
		                                                   intersectionWithPlane, 
		                                                   paramOfPlaneIntersection, 
		                                                   posOnLSForPlaneIntersection);
	
	// 1. Line segment is perpendicular to the plane of the circle
	//    and intersects the plane
	if(sPtOnProjectionOfLS == ePtOnProjectionOfLS && bIntersectPlane)
	{
		// Plane intersection
		pointsOnLineSegment[indexOfPtArray] = intersectionWithPlane;
		paramsOnLineSegment[indexOfPtArray] = paramOfPlaneIntersection;
		minDistances[indexOfPtArray] = computeMinDistFromPoint(pointsOnLineSegment[indexOfPtArray], pointsOnCircle[indexOfPtArray]);		
		return 1;
	}
	// 2. Line segment is perpendicular to the plane of the circle
	//    but does not intersect the plane
	else if(sPtOnProjectionOfLS == ePtOnProjectionOfLS && !bIntersectPlane)
	{
		// Min. dist. point from circle center
		lineSegment.computeMinDistFromPoint(m_center, pointsOnLineSegment[indexOfPtArray], paramsOnLineSegment[indexOfPtArray]);
		minDistances[indexOfPtArray] = computeMinDistFromPoint(pointsOnLineSegment[indexOfPtArray], pointsOnCircle[indexOfPtArray]);
		return 1;
	}
	
	// 3. Line segment is not perpendicular to the plane of the circle
	rg_Point3D sPtOnLS = lineSegment.getStartPt();
	rg_Point3D ePtOnLS = lineSegment.getEndPt();	
	// Two end points of line segment
	pointsOnLineSegment[indexOfPtArray] = sPtOnLS;
	paramsOnLineSegment[indexOfPtArray++] = 0.0;
	pointsOnLineSegment[indexOfPtArray] = ePtOnLS;
	paramsOnLineSegment[indexOfPtArray++] = 1.0;

	// Min. dist. point from circle center
	rg_Point3D   minPt;
	rg_REAL      paramOfMinDistPt;
	PosOnLineSegOrArc posOnLSForMinDistPt;
	lineSegment.computeMinDistFromPoint(m_center, minPt, paramOfMinDistPt, posOnLSForMinDistPt);
	if(posOnLSForMinDistPt == NONEXTREME_PT && !rg_EQ(paramOfMinDistPt, paramOfPlaneIntersection))
	{
		pointsOnLineSegment[indexOfPtArray] = minPt;
		paramsOnLineSegment[indexOfPtArray++] = paramOfMinDistPt;
	}
	// Plane intersection
	if(bIntersectPlane && posOnLSForPlaneIntersection == NONEXTREME_PT)
	{
		pointsOnLineSegment[indexOfPtArray] = intersectionWithPlane;
		paramsOnLineSegment[indexOfPtArray++] = paramOfPlaneIntersection;
	}
	// Check whether the projection of line segment intersects circle
	rg_Point3D  intersectionWithCircle[ 2 ];
	rg_REAL     paramOfCircleIntersection[ 2 ];
	rg_INT      numOfIntersections = -1;
	numOfIntersections = computeIntersectionWith(projectionOfLineSeg, 
		                                         intersectionWithCircle, 
		                                         paramOfCircleIntersection);
	// Projection of line segment has intersection with circle
	if(numOfIntersections > 0)
	{
		// circle intersection
		rg_Point3D dirVec = lineSegment.getDirVec();		
		for(int i = 0;i < numOfIntersections;i++)
		{
			if(rg_EQ(paramOfCircleIntersection[ i ], 0.0) ||
			   rg_EQ(paramOfCircleIntersection[ i ], 1.0) ||
			   rg_EQ(paramOfCircleIntersection[ i ], paramOfMinDistPt) ||
			   rg_EQ(paramOfCircleIntersection[ i ], paramOfPlaneIntersection) )
				continue;
			
			pointsOnLineSegment[indexOfPtArray] = sPtOnLS + dirVec * paramOfCircleIntersection[ i ];
			paramsOnLineSegment[indexOfPtArray++] = paramOfCircleIntersection[ i ];
		}
	}
	// Compute points on circle corresponding to local minima or maxima for distance with linesegment
	for(rg_INT i = 0;i < indexOfPtArray;i++)
	{
		minDistances[ i ] = computeMinDistFromPoint(pointsOnLineSegment[ i ], pointsOnCircle[ i ]);
		//// test
		//rg_REAL param;
		//rg_Point3D projection = lineSegment.computeProjectedPointForPoint3D(pointsOnCircle[ i ], param);
	}
	//// test
	//rg_Point3D ptOnLS = lineSegment.evaluatePt(0.5271448580);
	//rg_Point3D ptOnC;
	//computeMinDistFromPoint(ptOnLS, ptOnC);
	//rg_Point3D projection = planeOfCircle.projectPointOnPlane(ptOnLS);
	//rg_REAL distance = sqrt((ptOnLS-projection).squaredMagnitude() + (projection.distance(m_center) - m_radius) * (projection.distance(m_center) - m_radius));

	return indexOfPtArray;
}

rg_INT Circle3D::locatePointsWithLocalMinimaMaximaFromLineSegment(const LineSegment3D& lineSegment, 
																  rg_dList<rg_Point3D>& pointsOnCircle, 
																  rg_dList<rg_Point3D>& pointsOnLineSegment, 
																  rg_dList<rg_REAL>& paramsOnLineSegment, 
																  rg_dList<rg_REAL>& minDistances) const
{
	// The Plane of circle
	rg_Point3D    unitNormalVec = m_normal.getUnitVector();
	Plane         planeOfCircle(unitNormalVec, m_center);
	
	// Project line segment onto the plane of circle
	LineSegment3D projectionOfLineSeg;
	planeOfCircle.projectLineSegmentOnPlane(lineSegment, projectionOfLineSeg);
	rg_Point3D sPtOnProjectionOfLS = projectionOfLineSeg.getStartPt();
	rg_Point3D ePtOnProjectionOfLS = projectionOfLineSeg.getEndPt();
	
	// Check whether line segment intersects the plane
	rg_Point3D   intersectionWithPlane;
	rg_REAL      paramOfPlaneIntersection;
	PosOnLineSegOrArc posOnLSForPlaneIntersection;
	rg_FLAG      bIntersectPlane 
		= planeOfCircle.computeIntersectionWithLineSegment(lineSegment, 
		intersectionWithPlane, 
		paramOfPlaneIntersection, 
		posOnLSForPlaneIntersection);
	
	// 1. Line segment is perpendicular to the plane of the circle
	//    and intersects the plane
	if(sPtOnProjectionOfLS == ePtOnProjectionOfLS && bIntersectPlane)
	{
		// Plane intersection
		rg_Point3D pointOnCircle;		
		computeMinDistFromPoint(intersectionWithPlane, pointOnCircle);
		pointsOnLineSegment.add(intersectionWithPlane);
		pointsOnCircle.add(pointOnCircle);
		paramsOnLineSegment.add(paramOfPlaneIntersection);
		return 1;
	}
	// 2. Line segment is perpendicular to the plane of the circle
	//    but does not intersect the plane
	else if(sPtOnProjectionOfLS == ePtOnProjectionOfLS && !bIntersectPlane)
	{
		// Min. dist. point from circle center
		rg_Point3D pointOnLineSegment, pointOnCircle;
		rg_REAL    paramOnLineSegment;
		lineSegment.computeMinDistFromPoint(m_center, pointOnLineSegment, paramOnLineSegment);
		computeMinDistFromPoint(pointOnLineSegment, pointOnCircle);
		pointsOnLineSegment.add(pointOnLineSegment);
		pointsOnCircle.add(pointOnCircle);
		paramsOnLineSegment.add(paramOnLineSegment);
		return 1;
	}
	
	// 3. Line segment is not perpendicular to the plane of the circle
	rg_Point3D sPtOnLS = lineSegment.getStartPt();
	rg_Point3D ePtOnLS = lineSegment.getEndPt();	
	// Two end points of line segment
	pointsOnLineSegment.add(sPtOnLS);
	pointsOnLineSegment.add(ePtOnLS);
	paramsOnLineSegment.add(0.0);
	paramsOnLineSegment.add(1.0);

	// Min. dist. point from circle center
	rg_Point3D   minPt;
	rg_REAL      paramOfMinDistPt;
	PosOnLineSegOrArc posOnLSForMinDistPt;
	lineSegment.computeMinDistFromPoint(m_center, minPt, paramOfMinDistPt, posOnLSForMinDistPt);
	if(posOnLSForMinDistPt == NONEXTREME_PT && 
	   !rg_EQ(paramOfMinDistPt, paramOfPlaneIntersection))
	{
		pointsOnLineSegment.add(minPt);
		paramsOnLineSegment.add(paramOfMinDistPt);
	}
	// Plane intersection
	if(bIntersectPlane && posOnLSForPlaneIntersection == NONEXTREME_PT)
	{
		pointsOnLineSegment.add(intersectionWithPlane);
		paramsOnLineSegment.add(paramOfPlaneIntersection);
	}
	// Check whether the projection of line segment intersects circle
	rg_Point3D  intersectionWithCircle[ 2 ];
	rg_REAL     paramOfCircleIntersection[ 2 ];
	rg_INT      numOfIntersections = -1;
	numOfIntersections = computeIntersectionWith(projectionOfLineSeg, 
		                                         intersectionWithCircle, 
		                                         paramOfCircleIntersection);
	// Projection of line segment has intersection with circle
	if(numOfIntersections > 0)
	{
		// circle intersection
		rg_Point3D dirVec = lineSegment.getDirVec();		
		for(int i = 0;i < numOfIntersections;i++)
		{
			if(rg_EQ(paramOfCircleIntersection[ i ], 0.0) ||
			   rg_EQ(paramOfCircleIntersection[ i ], 1.0) ||
			   rg_EQ(paramOfCircleIntersection[ i ], paramOfMinDistPt) ||
			   rg_EQ(paramOfCircleIntersection[ i ], paramOfPlaneIntersection) )
				continue;
			
			pointsOnLineSegment.add(sPtOnLS + dirVec * paramOfCircleIntersection[ i ]);
			paramsOnLineSegment.add(paramOfCircleIntersection[ i ]);
		}
	}
	// Compute points on circle corresponding to local minima or maxima for distance with linesegment
	pointsOnLineSegment.reset4Loop();
	while(pointsOnLineSegment.setNext4Loop())
	{
		rg_Point3D ptOnCircle;
		minDistances.add(computeMinDistFromPoint(pointsOnLineSegment.getEntity(), ptOnCircle));
		pointsOnCircle.add(ptOnCircle);
	}
	rg_INT numOfMinmaMaxima = paramsOnLineSegment.getSize();
	return numOfMinmaMaxima;
}

