#include "rg_Triangle.h"
#include "rg_RelativeOp.h"
#include "Plane.h"


rg_Triangle::rg_Triangle()
{
}



rg_Triangle::rg_Triangle(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3)
{
    m_point[0] = pt1;
    m_point[1] = pt2;
    m_point[2] = pt3;
}



rg_Triangle::rg_Triangle(const rg_Triangle& triangle)
{
    m_point[0] = triangle.m_point[0];
    m_point[1] = triangle.m_point[1];
    m_point[2] = triangle.m_point[2];
}



rg_Triangle::~rg_Triangle()
{
}



rg_Triangle& rg_Triangle::operator =(const rg_Triangle& triangle)
{
    if ( this == &triangle ) {
        return *this;
    }

    m_point[0] = triangle.m_point[0];
    m_point[1] = triangle.m_point[1];
    m_point[2] = triangle.m_point[2];

    return *this;
}



rg_BOOL rg_Triangle::isPointOnTriangle(const rg_Point3D& point) const
{
    Plane basePlane(m_point[0], m_point[1], m_point[2]);
    rg_REAL dist = basePlane.distanceFromPoint( point );
    if ( !rg_ZERO( dist ) ) {
        return rg_FALSE;
    }

    rg_BOOL pointOnTriangle = rg_FALSE;
    rg_Point3D v1 = m_point[1] - m_point[0];
    rg_Point3D v2 = m_point[2] - m_point[0];
    rg_Point3D v  = point - m_point[0];


    double a[2][2] = { {v1.getX(), v2.getX()}, {v1.getY(), v2.getY()} };

    double det = a[0][0]*a[1][1] - a[0][1]*a[1][0];
    if ( rg_NZERO( det ) ) {
        double c[2] = { v.getX(), v.getY() };

        double x[2] = {0.0, 0.0};
        x[0] = (c[0]*a[1][1] - a[0][1]*c[1])/det;
        x[1] = (c[1]*a[0][0] - a[1][0]*c[0])/det;

        if ( rg_NNEG( x[0] ) && rg_NNEG( x[1] ) && rg_LE( x[0]+x[1], 1.0) ) {
            pointOnTriangle = rg_TRUE; 
        }

        return pointOnTriangle;
    }

    a[0][0] = v1.getY(); a[0][1] = v2.getY();
    a[1][0] = v1.getZ(); a[1][1] = v2.getZ();
    det     = a[0][0]*a[1][1] - a[0][1]*a[1][0];
    if ( rg_NZERO( det ) ) {
        double c[2] = { v.getY(), v.getZ() };

        double x[2] = {0.0, 0.0};
        x[0] = (c[0]*a[1][1] - a[0][1]*c[1])/det;
        x[1] = (c[1]*a[0][0] - a[1][0]*c[0])/det;

        if ( rg_NNEG( x[0] ) && rg_NNEG( x[1] ) && rg_LE( x[0]+x[1], 1.0) ) {
            pointOnTriangle = rg_TRUE; 
        }

        return pointOnTriangle;
    }

    a[0][0] = v1.getZ(); a[0][1] = v2.getZ();
    a[1][0] = v1.getX(); a[1][1] = v2.getX();
    det     = a[0][0]*a[1][1] - a[0][1]*a[1][0];
    if ( rg_NZERO( det ) ) {
        double c[2] = { v.getZ(), v.getX() };

        double x[2] = {0.0, 0.0};
        x[0] = (c[0]*a[1][1] - a[0][1]*c[1])/det;
        x[1] = (c[1]*a[0][0] - a[1][0]*c[0])/det;

        if ( rg_NNEG( x[0] ) && rg_NNEG( x[1] ) && rg_LE( x[0]+x[1], 1.0) ) {
            pointOnTriangle = rg_TRUE; 
        }

        return pointOnTriangle;
    }

    return pointOnTriangle;


    /*
    rg_Point3D v1 = m_point[1] - m_point[0];
    rg_Point3D v2 = m_point[2] - m_point[0];
    rg_Point3D v  = point - m_point[0];

    double a[3][2] = { {v1.getX(), v2.getX() },  {v1.getY(), v2.getY() }, {v1.getZ(), v2.getZ()} };
    double c[4] = {v.getX(), v.getY(), v.getZ() };

    int i=0;
    double det = a[i][0]*a[i+1][1] - a[i][1]*a[i+1][0];
    if ( rg_ZERO( det ) ) {
        return rg_FALSE;
    }

    double x1 = (c[i]*a[i+1][1] - a[i][1]*c[i+1])/det;
    double x2 = (c[i+1]*a[i][0] - a[i+1][0]*c[i])/det;

    rg_BOOL pointOnTriangle = rg_FALSE;
    if ( rg_NNEG( x1 ) && rg_NNEG( x2 ) && rg_LE( x1+x2, 1.0) ) {
        pointOnTriangle = rg_TRUE; 
    }

    return pointOnTriangle;
    */
}


    
rg_BOOL rg_Triangle::isCircleInTriangle(const Circle3D& circle) const
{
    rg_BOOL isCircleCompletelyIncludedInTriangle = rg_TRUE;

    if ( !isPointOnTriangle( circle.getCenter() ) ) {
        isCircleCompletelyIncludedInTriangle = rg_FALSE;
    }
    else {
        LineSegment3D edges[3] = { LineSegment3D(m_point[0], m_point[1]), 
                                   LineSegment3D(m_point[1], m_point[2]), 
                                   LineSegment3D(m_point[2], m_point[0]) };

        for ( rg_INT i=0; i<3; ++i ) {
            rg_REAL dist = edges[i].computeMinDistFromPoint( circle.getCenter() );
            if ( rg_LT( dist, circle.getRadius() ) ) {
                isCircleCompletelyIncludedInTriangle = rg_FALSE;
                break;
            }
        }
    }

    return isCircleCompletelyIncludedInTriangle;
}



bool    rg_Triangle::doesIntersectWith(const Sphere& sphere) const
{
    rg_REAL minimumDistanceToSphereCenter = computeMinDistFromPoint(sphere.getCenter());
    if (minimumDistanceToSphereCenter < sphere.getRadius()) {
        return true;
    }
    else {
        return false;
    }
}



double  rg_Triangle::area() const
{
    double length[3] = { m_point[0].distance(m_point[1]), m_point[1].distance(m_point[2]), m_point[2].distance(m_point[0]) };

    double k = (length[0] + length[1] + length[2])/2.0;
    return sqrt( k*(k-length[0])*(k-length[1])*(k-length[1]) );
}



rg_Point3D  rg_Triangle::computeNormalVector() const
{
    rg_Point3D vector[2] = { m_point[1] - m_point[0], m_point[2] - m_point[0] };
    rg_Point3D normal = vector[0].crossProduct(vector[1]);
    normal.normalize();

    return normal;
}



rg_INT  rg_Triangle::intersect(const LineSegment3D& lineSegment, rg_Point3D& intersectionPoint) const
{
    Plane basePlane(m_point[0], m_point[1], m_point[2]);

    rg_Point3D intersection;
    rg_INT numIntersect = basePlane.intersect( lineSegment, intersection );

    if ( numIntersect == 0 ) {
        return numIntersect;
    }

    if ( isPointOnTriangle(intersection) ) {
        intersectionPoint = intersection;
        return 1;
    }
    else {
        return 0;
    }

}


    
double  rg_Triangle::area(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3)
{
    double length[3] = { pt1.distance(pt2), pt2.distance(pt3), pt3.distance(pt1) };

    double k = (length[0] + length[1] + length[2])/2.0;
    return sqrt( k*(k-length[0])*(k-length[1])*(k-length[1]) );
}



rg_Triangle::Position    rg_Triangle::positionOnHyperplane( const Plane& hyperplane, const rg_REAL& res )
{
    int     numPtsInOpenLowerHalfSpace = 0;
    int     numPtsInOpenUpperHalfSpace = 0;
    int     numPtsOnHyperPlane         = 0;

    for ( int i=0; i<3; ++i ) {
        double distance = hyperplane.distanceFromPoint( m_point[i] );
        if ( rg_ZERO( distance, res ) ) {
            ++numPtsOnHyperPlane;
        }
        else if ( rg_POS( distance, res ) ) {
            ++numPtsInOpenUpperHalfSpace;
        }
        else {
            ++numPtsInOpenLowerHalfSpace;
        }
    }


    if ( numPtsInOpenLowerHalfSpace == 3 ) {
        return IN_LOWER_HALF_SPACE;
    }
    else if ( numPtsInOpenUpperHalfSpace == 3 ) {
        return IN_UPPER_HALF_SPACE;
    }
    else {
        return ON_HYPERPLANE;
    }
}



bool        rg_Triangle::intersectWithLowerHalfSpace( const Plane& hyperplane, list<rg_Point3D>& polygon, const rg_REAL& res )
{
    int     numPtsInOpenLowerHalfSpace = 0;
    int     numPtsInOpenUpperHalfSpace = 0;
    int     numPtsOnHyperPlane         = 0;

    bool    onHyperplane[3]     = { false, false, false };
    bool    onUpperHalfSpace[3] = { false, false, false };
    bool    onLowerHalfSpace[3] = { false, false, false };

    int i=0;
    for ( i=0; i<3; ++i ) {
        double distance = hyperplane.distanceFromPoint( m_point[i] );
        if ( rg_ZERO( distance, res ) ) {
            ++numPtsOnHyperPlane;
            onHyperplane[i] = true;
        }
        else if ( rg_POS( distance, res ) ) {
            ++numPtsInOpenUpperHalfSpace;
            onUpperHalfSpace[i] = true;
        }
        else {
            ++numPtsInOpenLowerHalfSpace;
            onLowerHalfSpace[i] = true;
        }
    }


    bool    hasIntersectionWithLowerHalfSpace = false;
    if ( numPtsInOpenUpperHalfSpace + numPtsOnHyperPlane == 3 ) {
        hasIntersectionWithLowerHalfSpace = false;
    }
    else if ( numPtsInOpenLowerHalfSpace + numPtsOnHyperPlane == 3 ) {
        hasIntersectionWithLowerHalfSpace = true;

        polygon.push_back( m_point[0] );
        polygon.push_back( m_point[1] );
        polygon.push_back( m_point[2] );
    }
    else if ( numPtsInOpenUpperHalfSpace == 1 && numPtsInOpenLowerHalfSpace == 2 ) {
        hasIntersectionWithLowerHalfSpace = true;

        if ( onUpperHalfSpace[0] ) {
            rg_Point3D intersection[2];
            hyperplane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[1] );

            polygon.push_back( intersection[0] );
            polygon.push_back( m_point[1]      );
            polygon.push_back( m_point[2]      );
            polygon.push_back( intersection[1] );
        }
        else if ( onUpperHalfSpace[1] ) {
            rg_Point3D intersection[2];
            hyperplane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[1] );

            polygon.push_back( intersection[0] );
            polygon.push_back( m_point[2]      );
            polygon.push_back( m_point[0]      );
            polygon.push_back( intersection[1] );
        }
        else { //if ( onUpperHalfSpace[2] ) {
            rg_Point3D intersection[2];
            hyperplane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[1] );

            polygon.push_back( intersection[0] );
            polygon.push_back( m_point[0]      );
            polygon.push_back( m_point[1]      );
            polygon.push_back( intersection[1] );
        }
    }
    else if ( numPtsInOpenUpperHalfSpace == 2 && numPtsInOpenLowerHalfSpace == 1 ) {
        hasIntersectionWithLowerHalfSpace = true;

        if ( onLowerHalfSpace[0] ) {
            rg_Point3D intersection[2];
            hyperplane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[1] );

            polygon.push_back( m_point[0]      );
            polygon.push_back( intersection[0] );
            polygon.push_back( intersection[1] );
        }
        else if ( onLowerHalfSpace[1] ) {
            rg_Point3D intersection[2];
            hyperplane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[1] );

            polygon.push_back( m_point[1]      );
            polygon.push_back( intersection[0] );
            polygon.push_back( intersection[1] );
        }
        else { //if ( onLowerHalfSpace[2] ) {
            rg_Point3D intersection[2];
            hyperplane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[1] );

            polygon.push_back( m_point[2]      );
            polygon.push_back( intersection[0] );
            polygon.push_back( intersection[1] );
        }
    }
    else {
        hasIntersectionWithLowerHalfSpace = true;
    }

    return hasIntersectionWithLowerHalfSpace;
}



bool        rg_Triangle::intersectWithLowerHalfSpace( const Plane& hyperplane, list<rg_Point3D>& polygon, LineSegment3D& lineSegToIntersect, const rg_REAL& res )
{
    int     numPtsInOpenLowerHalfSpace = 0;
    int     numPtsInOpenUpperHalfSpace = 0;
    int     numPtsOnHyperPlane         = 0;

    bool    onHyperplane[3]     = { false, false, false };
    bool    onUpperHalfSpace[3] = { false, false, false };
    bool    onLowerHalfSpace[3] = { false, false, false };

    int i=0;
    for ( i=0; i<3; ++i ) {
        double distance = hyperplane.distanceFromPoint( m_point[i] );
        if ( rg_ZERO( distance, res ) ) {
            ++numPtsOnHyperPlane;
            onHyperplane[i] = true;
        }
        else if ( rg_POS( distance, res ) ) {
            ++numPtsInOpenUpperHalfSpace;
            onUpperHalfSpace[i] = true;
        }
        else {
            ++numPtsInOpenLowerHalfSpace;
            onLowerHalfSpace[i] = true;
        }
    }


    bool    hasIntersectionWithLowerHalfSpace = false;
    if ( numPtsInOpenUpperHalfSpace + numPtsOnHyperPlane == 3 ) {
        hasIntersectionWithLowerHalfSpace = false;
    }
    else if ( numPtsInOpenLowerHalfSpace + numPtsOnHyperPlane == 3 ) {
        hasIntersectionWithLowerHalfSpace = true;

        polygon.push_back( m_point[0] );
        polygon.push_back( m_point[1] );
        polygon.push_back( m_point[2] );

        if ( numPtsOnHyperPlane == 2 ) {
            if ( onHyperplane[0] && onHyperplane[1] ) {
                lineSegToIntersect.setPoints( m_point[0], m_point[1] );
            }
            else if ( onHyperplane[1] && onHyperplane[2] ) {
                lineSegToIntersect.setPoints( m_point[1], m_point[2] );
            }
            else if ( onHyperplane[2] && onHyperplane[0] ) {
                lineSegToIntersect.setPoints( m_point[2], m_point[0] );
            }
        }
    }
    else if ( numPtsInOpenUpperHalfSpace == 1 && numPtsInOpenLowerHalfSpace == 2 ) {
        hasIntersectionWithLowerHalfSpace = true;

        rg_Point3D intersection[2];
        if ( onUpperHalfSpace[0] ) {
            hyperplane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[1] );

            polygon.push_back( intersection[0] );
            polygon.push_back( m_point[1]      );
            polygon.push_back( m_point[2]      );
            polygon.push_back( intersection[1] );
        }
        else if ( onUpperHalfSpace[1] ) {
            hyperplane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[1] );

            polygon.push_back( intersection[0] );
            polygon.push_back( m_point[2]      );
            polygon.push_back( m_point[0]      );
            polygon.push_back( intersection[1] );
        }
        else { //if ( onUpperHalfSpace[2] ) {
            hyperplane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[1] );

            polygon.push_back( intersection[0] );
            polygon.push_back( m_point[0]      );
            polygon.push_back( m_point[1]      );
            polygon.push_back( intersection[1] );
        }

        lineSegToIntersect.setPoints( intersection[0], intersection[1] );
    }
    else if ( numPtsInOpenUpperHalfSpace == 2 && numPtsInOpenLowerHalfSpace == 1 ) {
        hasIntersectionWithLowerHalfSpace = true;

        rg_Point3D intersection[2];
        if ( onLowerHalfSpace[0] ) {
            hyperplane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[1] );

            polygon.push_back( m_point[0]      );
            polygon.push_back( intersection[0] );
            polygon.push_back( intersection[1] );
        }
        else if ( onLowerHalfSpace[1] ) {
            hyperplane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[1] );

            polygon.push_back( m_point[1]      );
            polygon.push_back( intersection[0] );
            polygon.push_back( intersection[1] );
        }
        else { //if ( onLowerHalfSpace[2] ) {
            hyperplane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[0] );
            hyperplane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[1] );

            polygon.push_back( m_point[2]      );
            polygon.push_back( intersection[0] );
            polygon.push_back( intersection[1] );
        }

        lineSegToIntersect.setPoints( intersection[0], intersection[1] );
    }
    else {
        hasIntersectionWithLowerHalfSpace = true;
    }

    return hasIntersectionWithLowerHalfSpace;
}



bool        rg_Triangle::intersectWithPlane( const Plane& plane, LineSegment3D& lineSegToIntersect, const rg_REAL& res)
{
    int     numPtsInOpenLowerHalfSpace = 0;
    int     numPtsInOpenUpperHalfSpace = 0;
    int     numPtsOnHyperPlane         = 0;

    bool    onHyperplane[3]     = { false, false, false };
    bool    onUpperHalfSpace[3] = { false, false, false };
    bool    onLowerHalfSpace[3] = { false, false, false };

    int i=0;
    for ( i=0; i<3; ++i ) {
        double distance = plane.distanceFromPoint( m_point[i] );
        if ( rg_ZERO( distance, res ) ) {
            ++numPtsOnHyperPlane;
            onHyperplane[i] = true;
        }
        else if ( rg_POS( distance, res ) ) {
            ++numPtsInOpenUpperHalfSpace;
            onUpperHalfSpace[i] = true;
        }
        else {
            ++numPtsInOpenLowerHalfSpace;
            onLowerHalfSpace[i] = true;
        }
    }


    bool    hasIntersectionWithPlane = false;
    if ( numPtsInOpenUpperHalfSpace + numPtsOnHyperPlane == 3 ) {
        hasIntersectionWithPlane = false;
    }
    else if ( numPtsInOpenLowerHalfSpace == 3 ) {
        hasIntersectionWithPlane = false;
    }
    else if ( numPtsOnHyperPlane == 2 ) {
        hasIntersectionWithPlane = true;

        if ( onHyperplane[0] && onHyperplane[1] ) {
            lineSegToIntersect.setPoints( m_point[0], m_point[1] );
        }
        else if ( onHyperplane[1] && onHyperplane[2] ) {
            lineSegToIntersect.setPoints( m_point[1], m_point[2] );
        }
        else if ( onHyperplane[2] && onHyperplane[0] ) {
            lineSegToIntersect.setPoints( m_point[2], m_point[0] );
        }
    }
    else if ( numPtsInOpenUpperHalfSpace == 1 && numPtsInOpenLowerHalfSpace == 2 ) {
        hasIntersectionWithPlane = true;

        rg_Point3D intersection[2];
        if ( onUpperHalfSpace[0] ) {
            plane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[0] );
            plane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[1] );
        }
        else if ( onUpperHalfSpace[1] ) {
            plane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[0] );
            plane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[1] );
        }
        else { //if ( onUpperHalfSpace[2] ) {
            plane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[0] );
            plane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[1] );
        }

        lineSegToIntersect.setPoints( intersection[0], intersection[1] );
    }
    else if ( numPtsInOpenUpperHalfSpace == 2 && numPtsInOpenLowerHalfSpace == 1 ) {
        hasIntersectionWithPlane = true;

        rg_Point3D intersection[2];
        if ( onLowerHalfSpace[0] ) {
            plane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[0] );
            plane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[1] );
        }
        else if ( onLowerHalfSpace[1] ) {
            plane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[0] );
            plane.intersect( LineSegment3D( m_point[0], m_point[1]), intersection[1] );
        }
        else { //if ( onLowerHalfSpace[2] ) {
            plane.intersect( LineSegment3D( m_point[2], m_point[0]), intersection[0] );
            plane.intersect( LineSegment3D( m_point[1], m_point[2]), intersection[1] );
        }

        lineSegToIntersect.setPoints( intersection[0], intersection[1] );
    }
    else {
        hasIntersectionWithPlane = true;
    }

    return hasIntersectionWithPlane;
}



rg_REAL     rg_Triangle::computeMinDistFromPoint(const rg_Point3D& targetPt) const
{
    rg_REAL minimumDistance = DBL_MAX;

    Plane basePlane(m_point[0], m_point[1], m_point[2]);

    rg_Point3D ptProjectedToBase = basePlane.projectPointOnPlane(targetPt);
    if (isPointOnTriangle(ptProjectedToBase)) {
        minimumDistance = targetPt.distance(ptProjectedToBase);
    }
    else {
        LineSegment3D edge[3] = { LineSegment3D(m_point[0], m_point[1]), 
                                  LineSegment3D(m_point[1], m_point[2]),
                                  LineSegment3D(m_point[2], m_point[0]) };
        
        rg_REAL minDistToEdge = DBL_MAX;
        for (int i = 0; i < 3; ++i) {
            rg_REAL distToEdge = edge[i].computeMinDistFromPoint(targetPt);
            if (distToEdge < minDistToEdge) {
                minDistToEdge = distToEdge;
            }
        }

        minimumDistance = minDistToEdge;
    }

    return minimumDistance;
}


