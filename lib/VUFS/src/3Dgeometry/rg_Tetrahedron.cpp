#include "rg_Tetrahedron.h"
#include "Plane.h"
#include "rg_RelativeOp.h"



#include <vector>
using namespace std;


rg_Tetrahedron::rg_Tetrahedron()
{
}



rg_Tetrahedron::rg_Tetrahedron(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4)
{
    m_point[0] = pt1;
    m_point[1] = pt2; 
    m_point[2] = pt3;
    m_point[3] = pt4; 
}



rg_Tetrahedron::rg_Tetrahedron(const rg_Tetrahedron& tetrahedron)
{
    m_point[0] = tetrahedron.m_point[0];
    m_point[1] = tetrahedron.m_point[1]; 
    m_point[2] = tetrahedron.m_point[2];
    m_point[3] = tetrahedron.m_point[3]; 
}



rg_Tetrahedron::~rg_Tetrahedron()
{
}




rg_Point3D rg_Tetrahedron::getPoint(const rg_INT& i) const
{
    if ( i>=0 && i<4) {
        return m_point[i];
    }
    else {
        return rg_Point3D();
    }
}



void       rg_Tetrahedron::setPoint(const rg_INT& i, const rg_Point3D& pt)
{
    if ( i>=0 && i<4) {
        m_point[i] = pt;
    }
}




rg_Tetrahedron& rg_Tetrahedron::operator =(const rg_Tetrahedron& tetrahedron)
{
    if ( this == &tetrahedron ) {
        return *this;
    }

    m_point[0] = tetrahedron.m_point[0];
    m_point[1] = tetrahedron.m_point[1]; 
    m_point[2] = tetrahedron.m_point[2];
    m_point[3] = tetrahedron.m_point[3]; 

    return *this;
}




rg_BOOL rg_Tetrahedron::isPointInTetrahedron(const rg_Point3D& point) const
{
    rg_Point3D pt[7] = { m_point[0], m_point[1], m_point[2], m_point[3], m_point[0], m_point[1], m_point[2] };

    rg_BOOL isPointInThis = rg_TRUE;
    for ( rg_INT i=0; i<4; i++ ) {
        Plane plane( pt[i+1], pt[i+2], pt[i+3] );
        if ( plane.distanceFromPoint( pt[i] ) < 0.0 ) {
            plane.reverseNormal();
        }

        rg_REAL dist = plane.distanceFromPoint( point );
        if ( rg_NEG( dist ) ) {
            isPointInThis = rg_FALSE;
            break;
        }
    }

    return isPointInThis;
}



    
bool        rg_Tetrahedron::isLineSegmentInTetrahedron(const LineSegment3D& lineSegment) const
{
    if (    isPointInTetrahedron( lineSegment.getStartPt() ) 
         && isPointInTetrahedron( lineSegment.getEndPt()   ) ) {
        return true;
    }
    else {
        return false;
    }
}

      

bool        rg_Tetrahedron::isTriangleInTetrahedron(const rg_Triangle& triangle) const
{
    if (    isPointInTetrahedron( triangle.getPoint(0) ) 
         && isPointInTetrahedron( triangle.getPoint(1) )
         && isPointInTetrahedron( triangle.getPoint(2) ) ) {
        return true;
    }    
    else {
        return false;
    }
}

  

bool        rg_Tetrahedron::doesIntersectWith(const LineSegment3D& lineSegment) const
{
    if ( isLineSegmentInTetrahedron( lineSegment ) ) {
        return true;
    }
    else {
        for ( int i=0; i<4; ++i ) {
            rg_Triangle triangle = this->triangle(i);

            rg_Point3D  intersect;
            int numIntersection = triangle.intersect( lineSegment, intersect );
            if ( numIntersection != 0 ) {
                return true;
            }
        }

        return false;
    }
}

  

bool        rg_Tetrahedron::doesIntersectWith(const rg_Triangle& triangle) const
{
    if ( isTriangleInTetrahedron( triangle ) ) {
        return true;
    }
    else {
        LineSegment3D lineSegment[6] = { LineSegment3D( m_point[0], m_point[1] ), 
                                         LineSegment3D( m_point[0], m_point[2] ),
                                         LineSegment3D( m_point[0], m_point[3] ), 
                                         LineSegment3D( m_point[1], m_point[2] ),
                                         LineSegment3D( m_point[2], m_point[3] ), 
                                         LineSegment3D( m_point[3], m_point[1] ) };

        for ( int i=0; i<6; ++i ) {
            rg_Point3D  intersect;
            int numIntersection = triangle.intersect( lineSegment[i], intersect );
            if ( numIntersection != 0 ) {
                return true;
            }
        }

        return false;
    }
}



bool        rg_Tetrahedron::doesIntersectWith(const Sphere& sphere) const
{
    bool doesTetrahedronIntersectWithSphere = false;
    if (isPointInTetrahedron(sphere.getCenter())) {
        doesTetrahedronIntersectWithSphere = true;
    }
    else {
        for (int i = 0; i < 4; ++i) {
            rg_Triangle triangle = this->triangle(i);
            if (triangle.doesIntersectWith(sphere)) {
                doesTetrahedronIntersectWithSphere = true;
                break;
            }
        }
    }

    return doesTetrahedronIntersectWithSphere;
}



rg_Triangle rg_Tetrahedron::triangle(const rg_INT& i) const
{
    if ( i==0 ) {
        return rg_Triangle(m_point[1], m_point[2], m_point[3]);
    }
    else if ( i==1 ) {
        return rg_Triangle(m_point[0], m_point[3], m_point[2]);
    }
    else if ( i==2 ) {
        return rg_Triangle(m_point[0], m_point[1], m_point[3]);
    }
    else if ( i==3 ) {
        return rg_Triangle(m_point[0], m_point[2], m_point[1]);
    }
    else {
        return rg_Triangle();
    }
}



rg_REAL rg_Tetrahedron::computeSignedVolume() const
{
	// Get the coordinates of vertices
	rg_REAL x0 = m_point[0].getX();
	rg_REAL y0 = m_point[0].getY();
	rg_REAL z0 = m_point[0].getZ();

	rg_REAL x1 = m_point[1].getX();
	rg_REAL y1 = m_point[1].getY();
	rg_REAL z1 = m_point[1].getZ();

	rg_REAL x2 = m_point[2].getX();
	rg_REAL y2 = m_point[2].getY();
	rg_REAL z2 = m_point[2].getZ();

	rg_REAL x3 = m_point[3].getX();
	rg_REAL y3 = m_point[3].getY();
	rg_REAL z3 = m_point[3].getZ();

	rg_REAL signedVolume = 0.0;
	
	// Optimized by Maple 12
	// Before optimization: 23 additions 48 multiplications
	// Afeter optimization: 17 additions 16 multiplications 7 assignments

	rg_REAL t1, t2, t3, t4, t5, t6;
	t6 =  z0-z2;
	t5 =  z1-z0; 
	t4 =  z1-z3; 
	t3 =  z2-z1; 
	t2 =  z2-z3; 
	t1 = -z3+z0; 
	signedVolume 
		= ((-t5*y2-t6*y1-t3*y0)*x3+(t5*y3+t1*y1-t4*y0)*x2+(t6*y3-t1*y2+t2*y0)*x1+(t3*y3+t4*y2-t2*y1)*x0)/6.0;

	return signedVolume;
}


    
rg_Tetrahedron::Position    rg_Tetrahedron::positionOnHyperplane( const Plane& hyperplane, const rg_REAL& res )
{
    int     numPtsInOpenLowerHalfSpace = 0;
    int     numPtsInOpenUpperHalfSpace = 0;
    int     numPtsOnHyperPlane         = 0;

    for ( int i=0; i<4; ++i ) {
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


    if ( numPtsInOpenLowerHalfSpace == 4 ) {
        return IN_LOWER_HALF_SPACE;
    }
    else if ( numPtsInOpenUpperHalfSpace == 4 ) {
        return IN_UPPER_HALF_SPACE;
    }
    else {
        return ON_HYPERPLANE;
    }
}



bool        rg_Tetrahedron::hasIntersectionWithHyperplane(   const Plane& hyperplane, const rg_REAL& res )
{
    int     numPtsInOpenLowerHalfSpace = 0;
    int     numPtsInOpenUpperHalfSpace = 0;
    int     numPtsOnHyperPlane         = 0;

    for ( int i=0; i<4; ++i ) {
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


    if ( numPtsInOpenLowerHalfSpace==4 || numPtsInOpenUpperHalfSpace == 4 ) {
        return false;
    }
    else {
        return true;
    }
}



bool        rg_Tetrahedron::isInOpenLowerHalfSpace( const Plane& hyperplane, const rg_REAL& res )
{
    bool    isThisInLowerHalfSpace = true;
    for ( int i=0; i<4; ++i ) {
        if ( rg_NNEG( hyperplane.distanceFromPoint( m_point[i] ), res ) ) {
            isThisInLowerHalfSpace = false;
            break;
        }
    }

    return isThisInLowerHalfSpace;
}



bool        rg_Tetrahedron::isInOpenUpperHalfSpace( const Plane& hyperplane, const rg_REAL& res )
{
    bool    isThisInUpperHalfSpace = true;
    for ( int i=0; i<4; ++i ) {
        if ( rg_NPOS( hyperplane.distanceFromPoint( m_point[i] ), res ) ) {
            isThisInUpperHalfSpace = false;
            break;
        }
    }

    return isThisInUpperHalfSpace;
}



bool        rg_Tetrahedron::isInClosedLowerHalfSpace(      const Plane& hyperplane, const rg_REAL& res )
{
    bool    isThisInLowerHalfSpace = true;
    for ( int i=0; i<4; ++i ) {
        if ( rg_POS( hyperplane.distanceFromPoint( m_point[i] ), res ) ) {
            isThisInLowerHalfSpace = false;
            break;
        }
    }

    return isThisInLowerHalfSpace;
}



bool        rg_Tetrahedron::isInClosedUpperHalfSpace(      const Plane& hyperplane, const rg_REAL& res )
{
    bool    isThisInUpperHalfSpace = true;
    for ( int i=0; i<4; ++i ) {
        if ( rg_NEG( hyperplane.distanceFromPoint( m_point[i] ), res ) ) {
            isThisInUpperHalfSpace = false;
            break;
        }
    }

    return isThisInUpperHalfSpace;
}

    

bool        rg_Tetrahedron::intersectWithPlane( const Plane& hyperplane, list<rg_Point3D>& polygon, const rg_REAL& res )
{
    int     numPtsInOpenLowerHalfSpace = 0;
    int     numPtsInOpenUpperHalfSpace = 0;
    int     numPtsOnHyperPlane         = 0;

    bool    onHyperplane[4]     = { false, false, false, false };
    bool    onUpperHalfSpace[4] = { false, false, false, false };
    bool    onLowerHalfSpace[4] = { false, false, false, false };

    int i=0;
    for ( i=0; i<4; ++i ) {
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


    bool    hasThisIntersectionWithPlane = true;
    if ( numPtsInOpenLowerHalfSpace==4 || numPtsInOpenUpperHalfSpace == 4 ) {
        hasThisIntersectionWithPlane = false;
    }
    else if ( numPtsOnHyperPlane != 0 ) {
        //  special case!!! : the vertices of tetrahedron lie on the plane.
        for ( i=0; i<4; ++i ) {
            if ( onHyperplane[i] ) {
                polygon.push_back( m_point[i] );
            }
        }        
    }
    else {
        int numLineSeg = 0;
        rg_Point3D  lineSeg[4][2];
        if ( numPtsInOpenUpperHalfSpace == 1 ) {
            numLineSeg = 3;
            if ( onUpperHalfSpace[0] ) {
                lineSeg[0][0] = m_point[0]; lineSeg[0][1] = m_point[1]; 
                lineSeg[1][0] = m_point[0]; lineSeg[1][1] = m_point[2]; 
                lineSeg[2][0] = m_point[0]; lineSeg[2][1] = m_point[3]; 
            }
            else if ( onUpperHalfSpace[1] ) {
                lineSeg[0][0] = m_point[1]; lineSeg[0][1] = m_point[0]; 
                lineSeg[1][0] = m_point[1]; lineSeg[1][1] = m_point[3]; 
                lineSeg[2][0] = m_point[1]; lineSeg[2][1] = m_point[2]; 
            }
            else if ( onUpperHalfSpace[2] ) {
                lineSeg[0][0] = m_point[2]; lineSeg[0][1] = m_point[0]; 
                lineSeg[1][0] = m_point[2]; lineSeg[1][1] = m_point[1]; 
                lineSeg[2][0] = m_point[2]; lineSeg[2][1] = m_point[3]; 
            }
            else { // ( onUpperHalfSpace[0] ) {
                lineSeg[0][0] = m_point[3]; lineSeg[0][1] = m_point[0]; 
                lineSeg[1][0] = m_point[3]; lineSeg[1][1] = m_point[2]; 
                lineSeg[2][0] = m_point[3]; lineSeg[2][1] = m_point[1]; 
            }
        }
        else if ( numPtsInOpenLowerHalfSpace == 1 ) {
            numLineSeg = 3;
            if ( onLowerHalfSpace[0] ) {
                lineSeg[0][0] = m_point[0]; lineSeg[0][1] = m_point[1]; 
                lineSeg[1][0] = m_point[0]; lineSeg[1][1] = m_point[3]; 
                lineSeg[2][0] = m_point[0]; lineSeg[2][1] = m_point[2]; 
            }
            else if ( onLowerHalfSpace[1] ) {
                lineSeg[0][0] = m_point[1]; lineSeg[0][1] = m_point[0]; 
                lineSeg[1][0] = m_point[1]; lineSeg[1][1] = m_point[2]; 
                lineSeg[2][0] = m_point[1]; lineSeg[2][1] = m_point[3]; 
            }
            else if ( onLowerHalfSpace[2] ) {
                lineSeg[0][0] = m_point[2]; lineSeg[0][1] = m_point[0]; 
                lineSeg[1][0] = m_point[2]; lineSeg[1][1] = m_point[3]; 
                lineSeg[2][0] = m_point[2]; lineSeg[2][1] = m_point[1]; 
            }
            else { // ( onLowerHalfSpace[0] ) {
                lineSeg[0][0] = m_point[3]; lineSeg[0][1] = m_point[0]; 
                lineSeg[1][0] = m_point[3]; lineSeg[1][1] = m_point[1]; 
                lineSeg[2][0] = m_point[3]; lineSeg[2][1] = m_point[2]; 
            }
        }
        else { // ( numPtsInOpenUpperHalfSpace == 2 || numPtsInOpenLowerHalfSpace == 2 ) {
            numLineSeg = 4;
            if ( onUpperHalfSpace[0] && onUpperHalfSpace[1] ) {
                lineSeg[0][0] = m_point[0]; lineSeg[0][1] = m_point[2]; 
                lineSeg[1][0] = m_point[0]; lineSeg[1][1] = m_point[3]; 
                lineSeg[2][0] = m_point[1]; lineSeg[2][1] = m_point[3]; 
                lineSeg[3][0] = m_point[1]; lineSeg[3][1] = m_point[2]; 
            }
            else if ( onUpperHalfSpace[0] && onUpperHalfSpace[2] ) {
                lineSeg[0][0] = m_point[0]; lineSeg[0][1] = m_point[3]; 
                lineSeg[1][0] = m_point[0]; lineSeg[1][1] = m_point[1]; 
                lineSeg[2][0] = m_point[2]; lineSeg[2][1] = m_point[1]; 
                lineSeg[3][0] = m_point[2]; lineSeg[3][1] = m_point[3]; 
            }
            else if ( onUpperHalfSpace[0] && onUpperHalfSpace[3] ) {
                lineSeg[0][0] = m_point[0]; lineSeg[0][1] = m_point[1]; 
                lineSeg[1][0] = m_point[0]; lineSeg[1][1] = m_point[2]; 
                lineSeg[2][0] = m_point[3]; lineSeg[2][1] = m_point[2]; 
                lineSeg[3][0] = m_point[3]; lineSeg[3][1] = m_point[1]; 
            }
            else if ( onUpperHalfSpace[1] && onUpperHalfSpace[2] ) {
                lineSeg[0][0] = m_point[1]; lineSeg[0][1] = m_point[0]; 
                lineSeg[1][0] = m_point[1]; lineSeg[1][1] = m_point[3]; 
                lineSeg[2][0] = m_point[2]; lineSeg[2][1] = m_point[3]; 
                lineSeg[3][0] = m_point[2]; lineSeg[3][1] = m_point[0]; 
            }
            else if ( onUpperHalfSpace[1] && onUpperHalfSpace[3] ) {
                lineSeg[0][0] = m_point[1]; lineSeg[0][1] = m_point[2]; 
                lineSeg[1][0] = m_point[1]; lineSeg[1][1] = m_point[0]; 
                lineSeg[2][0] = m_point[3]; lineSeg[2][1] = m_point[0]; 
                lineSeg[3][0] = m_point[3]; lineSeg[3][1] = m_point[2]; 
            }
            else { // ( onUpperHalfSpace[2] && onUpperHalfSpace[3] ) {
                lineSeg[0][0] = m_point[2]; lineSeg[0][1] = m_point[0]; 
                lineSeg[1][0] = m_point[2]; lineSeg[1][1] = m_point[1]; 
                lineSeg[2][0] = m_point[3]; lineSeg[2][1] = m_point[1]; 
                lineSeg[3][0] = m_point[3]; lineSeg[3][1] = m_point[0]; 
            }
        }


        for ( i=0; i<numLineSeg; ++i ) {
            rg_Point3D intersection;
            hyperplane.intersect( LineSegment3D( lineSeg[i][0], lineSeg[i][1]), intersection );
            polygon.push_back( intersection );
        }
    }

    return hasThisIntersectionWithPlane;
}



rg_REAL rg_Tetrahedron::computeSignedVolume(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4)
{
	// Get the coordinates of vertices
	rg_REAL x0 = pt1.getX();
	rg_REAL y0 = pt1.getY();
	rg_REAL z0 = pt1.getZ();

	rg_REAL x1 = pt2.getX();
	rg_REAL y1 = pt2.getY();
	rg_REAL z1 = pt2.getZ();

	rg_REAL x2 = pt3.getX();
	rg_REAL y2 = pt3.getY();
	rg_REAL z2 = pt3.getZ();

	rg_REAL x3 = pt4.getX();
	rg_REAL y3 = pt4.getY();
	rg_REAL z3 = pt4.getZ();

	rg_REAL signedVolume = 0.0;
	
	// Optimized by Maple 12
	// Before optimization: 23 additions 48 multiplications
	// Afeter optimization: 17 additions 16 multiplications 7 assignments

	rg_REAL t1, t2, t3, t4, t5, t6;
	t6 =  z0-z2;
	t5 =  z1-z0; 
	t4 =  z1-z3; 
	t3 =  z2-z1; 
	t2 =  z2-z3; 
	t1 = -z3+z0; 
	signedVolume 
		= ((-t5*y2-t6*y1-t3*y0)*x3+(t5*y3+t1*y1-t4*y0)*x2+(t6*y3-t1*y2+t2*y0)*x1+(t3*y3+t4*y2-t2*y1)*x0)/6.0;

	return signedVolume;
}



rg_BOOL rg_Tetrahedron::checkOrientation() const
{
    Plane plane(m_point[1], m_point[2], m_point[3]);

    rg_REAL dist = plane.distanceFromPoint( m_point[0] );

    if ( dist >= 0.0 ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



