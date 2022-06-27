#include "Arc2D.h"
#include <math.h>

Arc2D::Arc2D()
{
}



Arc2D::Arc2D(const rg_Circle2D& circle, const rg_Point2D& startPoint, const rg_Point2D& endPoint)
: rg_Circle2D( circle )
{
    m_startPoint = startPoint;
    m_endPoint   = endPoint;
}



Arc2D::Arc2D(const Arc2D& arc)
: rg_Circle2D( arc )
{
    m_startPoint = arc.m_startPoint;
    m_endPoint   = arc.m_endPoint;
}



Arc2D::~Arc2D()
{
}




Arc2D& Arc2D::operator =(const Arc2D& arc)
{
    if ( this == &arc ) {
        return *this;
    }

    rg_Circle2D::operator =(arc);
    m_startPoint = arc.m_startPoint;
    m_endPoint   = arc.m_endPoint;
        
    return *this;
}


    
rg_BOOL Arc2D::operator==(const Arc2D& arc) const
{
    if ( rg_Circle2D::operator==(arc) && m_startPoint == arc.m_startPoint && m_endPoint == arc.m_endPoint ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }        
}

    
rg_REAL Arc2D::angle() const
{
    if ( m_startPoint == m_endPoint ) {
        return 2*rg_PI;
    }
    else {
        rg_Point2D center = getCenterPt();
        rg_Point2D vecCenterToStart = m_startPoint - center;
        rg_Point2D vecCenterToEnd   = m_endPoint - center;

        rg_REAL angleStartToEnd = angleFromVec1toVec2(vecCenterToStart, vecCenterToEnd);

        return angleStartToEnd;
    }
}

rg_BOOL Arc2D::isOnThisArc(const rg_Point2D& point) const
{
    if ( point == m_startPoint || point == m_endPoint ) {
        return rg_TRUE;
    }

    rg_Point2D center = getCenterPt();

    if ( rg_EQ( center.distance(point), getRadius() ) ) {
        if ( m_startPoint == m_endPoint ) {
            //  it means that this arc is a circle.
            return rg_TRUE;
        }

        rg_Point2D vecCenterToStart = m_startPoint - center;
        rg_Point2D vecCenterToEnd   = m_endPoint - center;
        rg_Point2D vecCenterToPoint = point - center;

        rg_REAL angleStartToEnd = angleFromVec1toVec2(vecCenterToStart, vecCenterToEnd);
        rg_REAL angleStartToPt  = angleFromVec1toVec2(vecCenterToStart, vecCenterToPoint);

        if ( angleStartToEnd > angleStartToPt ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



rg_BOOL Arc2D::isContainedIn(const Arc2D& bigArc) const
{
    rg_BOOL isThisArcContainedInBigArc = rg_FALSE;

    if ( rg_Circle2D::operator ==(bigArc) ) {
        rg_Point2D center = getCenterPt();

        rg_Point2D vecCenterToBigStart = bigArc.m_startPoint - center;
        rg_Point2D vecCenterToBigEnd   = bigArc.m_endPoint - center;
        rg_Point2D vecCenterToStart    = m_startPoint - center;
        rg_Point2D vecCenterToEnd      = m_endPoint - center;

        rg_REAL angleStartToBigEnd = angleFromVec1toVec2(vecCenterToBigStart, vecCenterToBigEnd);
        rg_REAL angleStartToStart  = angleFromVec1toVec2(vecCenterToBigStart, vecCenterToStart);
        rg_REAL angleStartToEnd    = angleFromVec1toVec2(vecCenterToBigStart, vecCenterToEnd);

        if ( rg_LE(angleStartToStart, angleStartToEnd) ) {
            if ( rg_LE(angleStartToEnd, angleStartToBigEnd) ) {
                isThisArcContainedInBigArc = rg_TRUE;
            }
        }
    }

    return isThisArcContainedInBigArc;
}



rg_INT Arc2D::intersect(const Arc2D& arc, rg_Point2D intersection[]) const
{
	rg_INT numIntersections = 0;

    if ( !rg_Circle2D::isIntersectWith(arc) ) {
        return numIntersections;
    }

    if ( rg_Circle2D::operator ==(arc) ) {
        Arc2D inputArc = arc;

        if ( inputArc.isContainedIn( *this ) ) {
            intersection[0]  = inputArc.m_startPoint;
            intersection[1]  = inputArc.m_endPoint;

	        numIntersections = rg_INT_INFINITY;
        }
        else if ( isContainedIn( inputArc ) ) {
            //  The intersection of two arcs is identical to this arc. 
            //  This arc is a proper subset of the input arc.
            intersection[0]  = m_startPoint;
            intersection[1]  = m_endPoint;

	        numIntersections = rg_INT_INFINITY;
        }
        else if ( isOnThisArc(inputArc.m_startPoint) && isOnThisArc(inputArc.m_endPoint) ) {
            intersection[0]  = inputArc.m_startPoint;
            intersection[1]  = m_endPoint;
            intersection[2]  = m_startPoint;
            intersection[3]  = inputArc.m_endPoint;

	        numIntersections = rg_INT_INFINITY;
        }
        else if ( isOnThisArc(inputArc.m_startPoint) ) {
            //  Two arcs are partially overlapped from a start point of the input arc to an end point of this arc.
            intersection[0]  = inputArc.m_startPoint;
            intersection[1]  = m_endPoint;

	        numIntersections = rg_INT_INFINITY;
        }
        else if ( isOnThisArc(inputArc.m_endPoint) ) {
            //  Two arcs are partially overlapped from a start point of this arc to an end point of the input arc.
            intersection[0]  = m_startPoint;
            intersection[1]  = inputArc.m_endPoint;

	        numIntersections = rg_INT_INFINITY;
        }
        else {
	        numIntersections = 0;
        }
    }
    else {
	    rg_Point2D* intersectInBaseCircle = rg_Circle2D::getIntersectPt(arc);	

        if ( rg_Circle2D::isTangentTo(arc) ) {
            if ( this->isOnThisArc( intersectInBaseCircle[0] ) ) {
                intersection[0]  = intersectInBaseCircle[0];
                numIntersections = 1;
            }
        }
        else {
            rg_INT pos = 0;
            if ( this->isOnThisArc( intersectInBaseCircle[0] ) ) {
                intersection[numIntersections]  = intersectInBaseCircle[0];
                numIntersections++;
            }

            if ( this->isOnThisArc( intersectInBaseCircle[1] ) ) {
                intersection[numIntersections]  = intersectInBaseCircle[1];
                numIntersections++;
            }
        }

        delete [] intersectInBaseCircle;
    }    

    return numIntersections;
}




void Arc2D::evaluatePointsOnArcGivenResolution(const rg_INT&  circleResolution, list<rg_Point2D>& pointsOnArc) const
{
    rg_Point2D center = getCenterPt();
    rg_REAL    radius = getRadius();

    if ( m_startPoint == m_endPoint ) {

		for(int i = 0; i < circleResolution; i++) {
			double angle = 2*rg_PI*(double)i/(double)circleResolution;
			rg_Point2D point(center.getX() + radius*cos(angle), center.getY() + radius*sin(angle));  

            pointsOnArc.push_back( point );
		}		
        pointsOnArc.push_back( *(pointsOnArc.begin()) );
    }
    else {
        rg_Point2D startVec = m_startPoint - center;
        rg_Point2D endVec   = m_endPoint - center;
	    rg_Point2D uVector = startVec.getUnitVector();
        rg_Point2D vVector(-uVector.getY(), uVector.getX());

        rg_REAL    arcAngle = angleFromVec1toVec2(startVec, endVec);
        rg_REAL    angle    = 0.0;


        rg_INT i = 0;
        do {            
            rg_Point2D point = center + radius*(cos(angle)*uVector + sin(angle)*vVector);
            pointsOnArc.push_back( point );

            i++;
			angle = 2*rg_PI*(double)i/(double)circleResolution;
        } while ( angle < arcAngle );

        pointsOnArc.push_back( m_endPoint );
    }
}
