#include "Sphere.h"
#include "rg_RelativeOp.h"
//#include "FunctionsForVoronoiDiagram3D.h"
#include "Plane.h"
#include "rg_GeoFunc.h"
#include "rg_TMatrix3D.h"
#include "rg_Circle2D.h"



///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
Sphere::Sphere()
: m_center(), m_radius(0.0)
{
}

Sphere::Sphere(const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, const rg_REAL& radius)
: m_center(x, y, z), m_radius(radius)
{
}

Sphere::Sphere(const rg_Point3D& center, const rg_REAL& radius)
: m_center(center), m_radius(radius)
{
}

Sphere::Sphere(const Sphere& aSphere)
: m_center(aSphere.m_center), m_radius(aSphere.m_radius)
{
}

Sphere::~Sphere()
{
}






///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
//  convert to inline functions (by Youngsong Cho 2007.05.25)
//rg_Point3D Sphere::getCenter() const
//{
//    return m_center;
//}
//
//rg_REAL    Sphere::getRadius() const
//{
//    return m_radius;
//}


///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void Sphere::setCenter(const rg_REAL& x, const rg_REAL& y, const rg_REAL& z)
{
    m_center.setPoint(x, y, z);
}

void Sphere::setCenter(const rg_Point3D& center)
{
    m_center.setPoint(center);
}

void Sphere::setRadius(const rg_REAL& radius)
{
    m_radius = radius;
}

void Sphere::setSphere(const rg_Point3D& center, const rg_REAL& radius)
{
    m_center.setPoint(center);
    m_radius = radius;
}

void Sphere::setSphere(const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, const rg_REAL& radius)
{
    m_center.setPoint(x, y, z);
    m_radius = radius;
}


///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
Sphere& Sphere::operator =(const Sphere& aSphere)
{
    if ( this == &aSphere )
        return *this;

    m_center = aSphere.m_center;
    m_radius = aSphere.m_radius;

    return *this;
}

bool Sphere::operator==(const Sphere& aSphere) const
{
    if ( m_center == aSphere.m_center )
        if ( rg_EQ( m_radius, aSphere.m_radius) )
            return rg_TRUE;
        else
            return rg_FALSE;
    else
        return rg_FALSE;

}

bool Sphere::operator <(const Sphere& aSphere) const
{
    //  4개의 sphere에 대한 함수 rg_SphereSetVoronoiDiagram::evaluate_Voronoi_Vertex()의 결과에 대한 
    //  numerical error를 줄이기 위해 4개의 sphere를 정렬하는 데에만  사용함.
    if ( rg_LT(m_radius, aSphere.m_radius, rg_SYSTEM_RES) )
    {
        return rg_TRUE;
    }
    else if ( rg_EQ(m_radius, aSphere.m_radius, rg_SYSTEM_RES) )
    {
        if ( rg_LT(m_center.getX(), aSphere.m_center.getX(), rg_SYSTEM_RES) )
            return rg_TRUE;
        else if ( rg_EQ(m_center.getX(), aSphere.m_center.getX(), rg_SYSTEM_RES) )
        {
            if ( rg_LT(m_center.getY(), aSphere.m_center.getY(), rg_SYSTEM_RES) )
                return rg_TRUE;
            else if ( rg_EQ(m_center.getY(), aSphere.m_center.getY(), rg_SYSTEM_RES) )
            {
                if ( rg_LT(m_center.getZ(), aSphere.m_center.getZ(), rg_SYSTEM_RES) )
                    return rg_TRUE;
                else
                    return rg_FALSE;
            }
            else 
                return rg_FALSE;
        }
        else 
            return rg_FALSE;
    }
    else 
        return rg_FALSE;
}




bool Sphere::operator >(const Sphere& aSphere) const
{
    //  4개의 sphere에 대한 함수 rg_SphereSetVoronoiDiagram::evaluate_Voronoi_Vertex()의 결과에 대한 
    //  numerical error를 줄이기 위해 4개의 sphere를 정렬하는 데에만  사용함.
    if ( rg_GT(m_radius, aSphere.m_radius, rg_MATH_RES) )
    {
        return rg_TRUE;
    }
    else if ( rg_EQ(m_radius, aSphere.m_radius, rg_MATH_RES) )
    {
        if ( rg_GT(m_center.getX(), aSphere.m_center.getX(), rg_MATH_RES) )
            return rg_TRUE;
        else if ( rg_EQ(m_center.getX(), aSphere.m_center.getX(), rg_MATH_RES) )
        {
            if ( rg_GT(m_center.getY(), aSphere.m_center.getY(), rg_MATH_RES) )
                return rg_TRUE;
            else if ( rg_EQ(m_center.getY(), aSphere.m_center.getY(), rg_MATH_RES) )
            {
                if ( rg_GT(m_center.getZ(), aSphere.m_center.getZ(), rg_MATH_RES) )
                    return rg_TRUE;
                else
                    return rg_FALSE;
            }
            else 
                return rg_FALSE;
        }
        else 
            return rg_FALSE;
    }
    else 
        return rg_FALSE;
}

///////////////////////////////////////////////////////////////////////////////
//
//  geometric operators..
rg_REAL Sphere::distance(const Sphere& aSphere) const
{
    rg_REAL distanceOfTwoCenters = m_center.distance(aSphere.m_center);
    rg_REAL sumOfTwoRadii        = m_radius + aSphere.m_radius;

    if ( rg_GT( distanceOfTwoCenters, sumOfTwoRadii )  )
        return (distanceOfTwoCenters - sumOfTwoRadii);
    else
        return 0.0;
}



rg_REAL Sphere::distance(const rg_Point3D& point) const
{
    rg_REAL distBetCenterAndPt = m_center.distance(point);

    return (distBetCenterAndPt - m_radius);
}



rg_REAL Sphere::distance(rg_REAL* plane) const
{
    rg_REAL distBetCenterAndPlane 
                = m_center.getX()*plane[0] + m_center.getY()*plane[1] + m_center.getZ()*plane[2] + plane[3];


    if ( distBetCenterAndPlane > 0.0  )
        return (distBetCenterAndPlane - m_radius);
    else if  ( distBetCenterAndPlane < 0.0  )
        return (distBetCenterAndPlane + m_radius);
    else
        return m_radius;
}



rg_REAL Sphere::distanceToSphere(const Sphere& aSphere) const
{
    rg_REAL distanceOfTwoCenters = m_center.distance(aSphere.m_center);
    rg_REAL sumOfTwoRadii        = m_radius + aSphere.m_radius;

    return (distanceOfTwoCenters - sumOfTwoRadii);
}



rg_REAL Sphere::distanceToPlane(rg_REAL* plane) const
{
    rg_REAL distBetCenterAndPlane 
                = m_center.getX()*plane[0] + m_center.getY()*plane[1] + m_center.getZ()*plane[2] + plane[3];

    return distBetCenterAndPlane-m_radius;
}



rg_REAL Sphere::distanceBetweenCenterAndPlane(rg_REAL* plane) const
{
    rg_REAL distBetCenterAndPlane 
                = m_center.getX()*plane[0] + m_center.getY()*plane[1] 
                  + m_center.getZ()*plane[2] + plane[3];

    return distBetCenterAndPlane;
}



rg_Point3D Sphere::project(const rg_Point3D& point) const
{
    rg_Point3D vectorToPoint = point - m_center;
    vectorToPoint.normalize();

    rg_Point3D projectedPoint = m_center + (m_radius*vectorToPoint);

    return projectedPoint;
}


    
rg_BOOL Sphere::isEqual(const Sphere& aSphere, const rg_REAL& tolerance) const
{
    if ( m_center.isEqual( aSphere.m_center, tolerance ) ) {
        if ( rg_EQ( m_radius, aSphere.m_radius, tolerance) ) {
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

//  area and volume..
// rg_REAL Sphere::computeArea() const
// {
// 	rg_REAL area = 4.0*rg_PI*m_radius*m_radius;
// 	return area;
// }



rg_REAL Sphere::computeAreaOfSphericalTriangle(const rg_Point3D& pt1OnSphere, const rg_Point3D& pt2OnSphere, const rg_Point3D& pt3OnSphere) const
{
    /*
	Plane plane[3];
	plane[0].definePlaneByThreePoints(pt1OnSphere, pt2OnSphere, m_center);
	plane[1].definePlaneByThreePoints(pt2OnSphere, pt3OnSphere, m_center);
	plane[2].definePlaneByThreePoints(pt3OnSphere, pt1OnSphere, m_center);

    if ( plane[0].distanceFromPoint(pt3OnSphere) < 0.0 )
        plane[0].reverseNormal();
    if ( plane[1].distanceFromPoint(pt1OnSphere) < 0.0 )
        plane[1].reverseNormal();
    if ( plane[2].distanceFromPoint(pt2OnSphere) < 0.0 )
        plane[2].reverseNormal();


	rg_Point3D normal[3] = { plane[0].getNormal(), plane[1].getNormal(), plane[2].getNormal() };	
	rg_REAL    angle[3]  = { rg_PI-normal[0].angle(normal[1]), rg_PI-normal[1].angle(normal[2]), rg_PI-normal[2].angle(normal[0]) };

	rg_REAL sumOfAngles = angle[0] + angle[1] + angle[2];
	rg_REAL sphericalExcess = sumOfAngles - rg_PI;

	rg_REAL area = sphericalExcess*m_radius*m_radius;
*/
    rg_Point3D normal[3] = {m_center-pt1OnSphere, m_center-pt2OnSphere, m_center-pt3OnSphere};
	Plane plane[3];
	plane[0].definePlaneByNormalAndPassingPoint(normal[0], pt1OnSphere);
	plane[1].definePlaneByNormalAndPassingPoint(normal[1], pt2OnSphere);
	plane[2].definePlaneByNormalAndPassingPoint(normal[2], pt3OnSphere);

    rg_REAL    angle[3];
    rg_Point3D pt[5] = {pt1OnSphere, pt2OnSphere, pt3OnSphere, pt1OnSphere, pt2OnSphere};
    for (rg_INT i=0; i<3; i++) {
        rg_Point3D ptOnPlane[2] = { plane[i].projectPointOnPlane(pt[i+1]), plane[i].projectPointOnPlane(pt[i+2])};

        rg_Point3D vec[2] = { ptOnPlane[0]-pt[i], ptOnPlane[1]-pt[i]};
        angle[i] = vec[0].angle(vec[1]);
    }

	rg_REAL sumOfAngles = angle[0] + angle[1] + angle[2];
	rg_REAL sphericalExcess = sumOfAngles - rg_PI;

	rg_REAL area = sphericalExcess*m_radius*m_radius;

	return area;
}


// DSKIM: inline function
rg_REAL Sphere::computeAreaOfSphericalCap(const rg_REAL& h) const
{
	rg_REAL area = 2.0*rg_PI*m_radius*h;

	return area;
}



rg_REAL Sphere::computeAreaOfSphericalCap(const Plane& plane) const
{
	rg_REAL area = 0.0;
	rg_REAL dist = plane.distanceFromPoint( m_center );
	
	if ( rg_ABS(dist) >= m_radius )  {
		area = 0.0;
	}
	else {
		if ( dist < 0.0 ) {
			rg_REAL h = m_radius + dist;
			area = computeAreaOfSphericalCap( h );
		}
		else {
			rg_REAL h = m_radius - dist;
			area = computeArea() - computeAreaOfSphericalCap( h );
		}
	}

	return area;
}



rg_REAL Sphere::computeAreaOfSphericalCap(const Sphere& sphere) const
{
	Circle3D cirlce;
	SphereIntersectionType type = intersect(sphere, cirlce);
	Plane plane = cirlce.getPlaneContainingThisCircle();

	rg_REAL area = 0.;
	if(type == SI_INTERSECT)
		area = computeAreaOfSphericalCap(plane);

	return area;
}



// DSKIM: computeBoundaryAreaOfIntersectingGostoneVolume()
rg_REAL Sphere::computeIntersectionArea(const Sphere& sphere) const
{
	Circle3D cirlce;
	SphereIntersectionType type = intersect(sphere, cirlce);
	Plane plane = cirlce.getPlaneContainingThisCircle();

	rg_REAL area = 0.;
	if(type == SI_INTERSECT)
		area = (computeAreaOfSphericalCap(plane) + sphere.computeAreaOfSphericalCap(plane.getReversedPlane()));
	
	return area;
}



rg_REAL Sphere::ratioOfOverlappedSphere( const Sphere& sphere ) const
{
    rg_REAL distBetweenCenters = m_center.distance(sphere.m_center);

    if (rg_POS(m_radius + sphere.m_radius - distBetweenCenters))
    {
        rg_REAL smallRadius = -1.0;
        if (m_radius > sphere.m_radius)
            smallRadius = sphere.m_radius;
        else
            smallRadius = m_radius;

        return (m_radius + sphere.m_radius - distBetweenCenters) / (2.0 *smallRadius);
    }
    else
        return 0.0;
}


// rg_REAL Sphere::computeVolume() const
// {
// 	rg_REAL volume = 4.0*rg_PI*m_radius*m_radius*m_radius/3.0;
// 	return volume;
// }



// spherical triangular cone
rg_REAL Sphere::computeVolumeOfSphericalTriangle(const rg_Point3D& pt1OnSphere, const rg_Point3D& pt2OnSphere, const rg_Point3D& pt3OnSphere) const
{
	rg_REAL area = computeAreaOfSphericalTriangle(pt1OnSphere, pt2OnSphere, pt3OnSphere);

	//rg_REAL volume = computeVolume()*area/computeArea();
	rg_REAL volume = area*m_radius/3.0;

	return volume;
}


// spherical circular cone
// DSKIM: inline function
rg_REAL Sphere::computeVolumeOfSphericalCap(const rg_REAL& h) const
{
	rg_REAL volume = rg_PI*h*h*(3.0*m_radius - h)/3.0;

	return volume;
}



rg_REAL Sphere::computeVolumeOfSphericalCap(const Plane& plane) const
{
	rg_REAL volume = 0.0;
	rg_REAL dist = plane.distanceFromPoint( m_center );
	
	if ( rg_ABS(dist) >= m_radius )  {
		volume = 0.0;
	}
	else {
		if ( dist < 0.0 ) {
			rg_REAL h = m_radius + dist;
			volume = computeVolumeOfSphericalCap( h );
		}
		else {
			rg_REAL h = m_radius - dist;
			volume = computeVolume() - computeVolumeOfSphericalCap( h );
		}
	}

	return volume;
}



rg_REAL Sphere::computeVolumeBySphericalCapAndPlanePassedCenter(const Plane& planeForCap, const Plane& planeOnCenter) const
{
    rg_REAL chi = rg_PI - planeForCap.getNormal().angle( planeOnCenter.getNormal() );
    rg_REAL z   = rg_ABS( planeForCap.distanceFromPoint( m_center ) );
    rg_REAL r   = m_radius;

    rg_REAL radiusOfCap   = sqrt( r*r - z*z );

    rg_REAL volume = ( 2.0*r*r*r*acos( r*cos(chi)/radiusOfCap )
                       - z*(3.0*r*r - z*z)*acos( z/(radiusOfCap*tan(chi)) )
                       + z*z*sqrt(r*r - z*z/(sin(chi)*sin(chi)))/tan(chi) )/3.0;

    return volume;
}


rg_REAL Sphere::computeIntersectingVolumeA(const Sphere& aSphere) const
{
	rg_REAL distBTWCenters  = m_center.distance(aSphere.m_center);
	rg_REAL intersectingVolume = 0.0;
	intersectingVolume = rg_GeoFunc::computeXVolumeOfTwoSpheres(m_radius, aSphere.m_radius, distBTWCenters);

	return intersectingVolume;
}



rg_REAL Sphere::computeIntersectingVolume(const Sphere& aSphere) const
{
	// DSKIM: remove initialization
    rg_REAL intersectingVolume = 0.0;

	// DSKIM: remove this statement: this test should be done to verify the input data is valid or not.
	// Therefore, if this should be done, this should be done at the very beginning when loading the data
    rg_FLAG isContained = isContainedIn(aSphere);

    if ( isContained == rg_TRUE ) {
        intersectingVolume = computeVolume();
    }
    else if ( isContained == (rg_TRUE+rg_TRUE) ) {
        intersectingVolume = aSphere.computeVolume();    
    }
    else {
        if ( isThereIntersectionWith(aSphere) == rg_TRUE )  {
            rg_REAL distBtwCenters = m_center.distance(aSphere.m_center);

	        rg_REAL lenSide[3]     = {m_radius, aSphere.m_radius, distBtwCenters };	
	        rg_REAL halfPerimeter = (lenSide[0] + lenSide[1] + lenSide[2])*0.5;
	        rg_REAL squaredArea   = halfPerimeter*(halfPerimeter-lenSide[0])*(halfPerimeter-lenSide[1])*(halfPerimeter-lenSide[2]);
	        rg_REAL areaOfTriangle = sqrt( squaredArea );
	        
            rg_REAL radiusIntersectCircle = 2.0*areaOfTriangle/lenSide[2];

            rg_REAL distToIntersectPlane[2];
            //distToIntersectPlane[0] = sqrt(m_radius*m_radius - radiusIntersectCircle*radiusIntersectCircle);
			rg_REAL sqdDist = m_radius*m_radius - radiusIntersectCircle*radiusIntersectCircle;
			if(sqdDist > 0)
				distToIntersectPlane[0] = sqrt(sqdDist);
			else
				distToIntersectPlane[0] = 0.;

            distToIntersectPlane[1] = distBtwCenters-distToIntersectPlane[0];
            
            intersectingVolume += computeVolumeOfSphericalCap( m_radius-distToIntersectPlane[0] );
            intersectingVolume += aSphere.computeVolumeOfSphericalCap( aSphere.m_radius-distToIntersectPlane[1] );
        }
    }

    return intersectingVolume;

}

rg_REAL Sphere::computeIntersectingVolumeB(const Sphere& aSphere) const
{
	// DSKIM: remove initialization
//     rg_REAL intersectingVolume = 0.0;

	// DSKIM: remove this statement: this test should be done to verify the input data is valid or not.
	// Therefore, if this should be done, this should be done at the very beginning when loading the data
//     rg_FLAG isContained = isContainedIn(aSphere);

//     if ( isContained == rg_TRUE ) {
//         intersectingVolume = computeVolume();
//     }
//     else if ( isContained == (rg_TRUE+rg_TRUE) ) {
//         intersectingVolume = aSphere.computeVolume();    
//     }
//     else {
//         if ( isThereIntersectionWith(aSphere) == rg_TRUE )  {



            rg_REAL distBtwCenters = m_center.distance(aSphere.m_center);

	        rg_REAL lenSide[3]     = {m_radius, aSphere.m_radius, distBtwCenters };	
	        rg_REAL halfPerimeter = (lenSide[0] + lenSide[1] + lenSide[2])*0.5;
	        rg_REAL squaredArea   = halfPerimeter*(halfPerimeter-lenSide[0])*(halfPerimeter-lenSide[1])*(halfPerimeter-lenSide[2]);
	        rg_REAL areaOfTriangle = sqrt( squaredArea );
	        
            rg_REAL radiusIntersectCircle = 2.0*areaOfTriangle/lenSide[2];

            rg_REAL distToIntersectPlane[2];
            //distToIntersectPlane[0] = sqrt(m_radius*m_radius - radiusIntersectCircle*radiusIntersectCircle);
			rg_REAL sqdDist = m_radius*m_radius - radiusIntersectCircle*radiusIntersectCircle;
			if(sqdDist > 0)
				distToIntersectPlane[0] = sqrt(sqdDist);
			else
				distToIntersectPlane[0] = 0.;

            distToIntersectPlane[1] = distBtwCenters-distToIntersectPlane[0];
            
			rg_REAL intersectingVolume = computeVolumeOfSphericalCap( m_radius-distToIntersectPlane[0] );
            intersectingVolume += aSphere.computeVolumeOfSphericalCap( aSphere.m_radius-distToIntersectPlane[1] );



//         }
//     }

    return intersectingVolume;

}


rg_REAL Sphere::computeIntersectingVolumeWithSameRadiusSpheres(const Sphere& sphere1, const Sphere& sphere2) const
{
    //  A : this sphere
    //  B : sphere1
    //  C : sphere2
    rg_REAL a = sphere1.getCenter().distance( sphere2.getCenter() );  //  distance Btw centers of B & C
    rg_REAL b = sphere2.getCenter().distance( m_center );             //  distance Btw centers of C & A
    rg_REAL c = m_center.distance( sphere1.getCenter() );             //  distance Btw centers of A & B
    rg_REAL r = m_radius;
    rg_REAL q1 = b*b + c*c - a*a;
    rg_REAL q2 = c*c + a*a - b*b;
    rg_REAL q3 = a*a + b*b - c*c;
    rg_REAL sqrw_ABC = 2.0*r*r*(c*c*a*a* + c*c*b*b + a*a*b*b) - r*r*(a*a*a*a + b*b*b*b + c*c*c*c) - c*c*a*a*b*b;
    rg_REAL w_ABC = sqrt(sqrw_ABC);
    //rg_REAL w_ABC = sqrt(2.0*r*r*(c*c*a*a* + c*c*b*b + a*a*b*b) - r*r*(a*a*a*a + b*b*b*b + c*c*c*c) - c*c*a*a*b*b);
    
    rg_REAL volume = w_ABC/6.0
                     + (1.0/12.0)*(   (a*a*a - 12.0*a*r*r)*atan( 2.0*w_ABC/(a*q1) ) 
                                    + (b*b*b - 12.0*b*r*r)*atan( 2.0*w_ABC/(b*q2) ) 
                                    + (c*c*c - 12.0*c*r*r)*atan( 2.0*w_ABC/(c*q3) ) ) 
                     + (2.0/3.0)*r*r*r*( 2.0*atan(w_ABC/(r*q2)) + 2.0*atan(w_ABC/(r*q3)) + 2.0*atan(w_ABC/(r*q3)) );

    return volume;
}


rg_REAL Sphere::computeIntersectionVolumeWith(const Sphere& sphere1, const Sphere& sphere2) const
{
    rg_REAL volumeIntersectedByThreeSpheres  = 0.0;

	// distance between the 2nd sphere and 3rd sphere
    rg_REAL a    = sphere2.m_center.distance(m_center);  
    // distance between the 3rd sphere and 1st sphere
    rg_REAL b    = m_center.distance(sphere1.m_center);   
    // distance between the 1st sphere and 2nd sphere
    rg_REAL c    = sphere1.m_center.distance(sphere2.m_center);   

    rg_REAL r1   = sphere1.m_radius;  // r1   Radius of the 1st sphere
    rg_REAL r2   = sphere2.m_radius;  // r2   Radius of the 2nd sphere
    rg_REAL r3   = m_radius;  // r3   Radius of the 3rd sphere

    
    rg_REAL a_squared   = a*a;  // square of a
    rg_REAL b_squared   = b*b;  // square of b
    rg_REAL c_squared   = c*c;  // square of c

    rg_REAL r1_squared  = r1*r1;    // rr1  Square of r1, or r1*r1
    rg_REAL r2_squared  = r2*r2;    // rr2  Square of r2, or r2*r2
    rg_REAL r3_squared  = r3*r3;    // rr3  Square of r3, or r3*r3

   
    // in order to distinguish the case that whether relationship is intersected or not.

    //////////////////////////////////////////////////////
    //  w^2=                                            //
    //                                                  //
    //  0.5* [ [0.0,    c^2,    b^2,    r1^2    1.0 ]   //
    //         [c^2,    0.0,    a^2,    r2^2    1.0 ]   //
    //         [b^2,    a^2,    0.0,    r3^2    1.0 ]   //
    //         [r1^2,   r2^2,   r3^2,   0.0,    1.0 ]   //
    //         [1.0,    1.0,    1.0,    1.0,    0.0 ] ] //
    //                                                  //
    //////////////////////////////////////////////////////

    rg_REAL w_squared =  
    /*    
            c_squared   *a_squared  *r3_squared   
          + c_squared   *r2_squared *r3_squared  
          + a_squared   *r2_squared *b_squared   
          + a_squared   *r3_squared *r1_squared  
          - a_squared   *r2_squared *r3_squared  
          + r2_squared  *b_squared  *r3_squared  
          + r2_squared  *a_squared  *r1_squared  
          - c_squared   *r3_squared *r3_squared  
          - a_squared   *a_squared  *r1_squared   
          - r2_squared  *r2_squared *b_squared
          - c_squared   *c_squared  *r3_squared   
          - r2_squared  *b_squared  *b_squared   
          - a_squared   *r1_squared *r1_squared  
          - c_squared   *a_squared  *b_squared    
          + c_squared   *a_squared  *r1_squared   
          + c_squared   *r2_squared *b_squared   
          - c_squared   *r2_squared *r1_squared  
          + c_squared   *b_squared  *r3_squared   
          + c_squared   *r3_squared *r1_squared  
          + b_squared   *r2_squared *r1_squared 
          + b_squared   *a_squared  *r1_squared   
          - b_squared   *r3_squared *r1_squared;
    */
    /*
        //////////////////////////////////////////////////////////////////////////
        //opt.
            ((r1_squared-b_squared)*b_squared-r2_squared*b_squared)*r2_squared
          + ( r2_squared*b_squared-b_squared*r1_squared)*r3_squared
          + (-r1_squared*r1_squared+b_squared*r1_squared+(r1_squared+b_squared)
            *r2_squared+(r1_squared-r2_squared)*r3_squared-a_squared*r1_squared)*a_squared
          + ((-r1_squared+b_squared)*r2_squared+(r2_squared+r1_squared+b_squared-r3_squared)
            *r3_squared+(r3_squared-b_squared+r1_squared)*a_squared-r3_squared*c_squared)*c_squared;
        */
        computeW_SquaredForTripleIntersection( a_squared, b_squared, c_squared, r1_squared, r2_squared, r3_squared);

    if( w_squared > 0.0 )  // sphere 1, 2 and 3 are intersected simultaneously
    {        
        rg_REAL w  = sqrt(w_squared);

        rg_REAL e1 = (r2_squared-r3_squared)/a_squared;
        rg_REAL e2 = (r3_squared-r1_squared)/b_squared;
        rg_REAL e3 = (r1_squared-r2_squared)/c_squared;

        rg_REAL q1 = a*( b_squared+c_squared-a_squared + r2_squared+r3_squared-(2.0*r1_squared) + e1*(b_squared-c_squared) );
        rg_REAL q2 = b*( c_squared+a_squared-b_squared + r3_squared+r1_squared-(2.0*r2_squared) + e2*(c_squared-a_squared) );
        rg_REAL q3 = c*( a_squared+b_squared-c_squared + r1_squared+r2_squared-(2.0*r3_squared) + e3*(a_squared-b_squared) );

        volumeIntersectedByThreeSpheres = 
            + w/6.0
            - a/2.0*( r2_squared + r3_squared - a_squared*(1.0/6.0-e1*e1/2.0))
			         *atan2( 1.0, q1/(2.0*w) )
            - b/2.0*( r3_squared + r1_squared - b_squared*(1.0/6.0-e2*e2/2.0))
			         *atan2( 1.0, q2/(2.0*w) )
            - c/2.0*( r1_squared + r2_squared - c_squared*(1.0/6.0-e3*e3/2.0))
			         *atan2( 1.0, q3/(2.0*w) )
            + 2.0/3.0*r1_squared*r1*( atan2( 1.0, (r1*q2) / (b*w*(1.0-e2)) ) 
			+ atan2( 1.0, (r1*q3) / (c*w*(1.0+e3)) ) )
            + 2.0/3.0*r2_squared*r2*( atan2( 1.0, (r2*q3) / (c*w*(1.0-e3)) ) 
			+ atan2( 1.0, (r2*q1) / (a*w*(1.0+e1)) ) )
            + 2.0/3.0*r3_squared*r3*( atan2( 1.0, (r3*q1) / (a*w*(1.0-e1)) ) 
			+ atan2( 1.0, (r3*q2) / (b*w*(1.0+e2)) ) );   
		//opt.
    }
    else if( w_squared < 0.0 ) // three exception cases
    {
		rg_REAL t_squared = (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c);
        rg_REAL t  = 0.0;
        
 		if(t_squared  >= 0.0 )        
        {
            t  = sqrt( t_squared );
        }
        else
        {
            return volumeIntersectedByThreeSpheres;
        }
               
        rg_REAL t1_squared = (a+r2+r3)*(-a+r2+r3)*(a-r2+r3)*(a+r2-r3);
        rg_REAL t2_squared = (b+r3+r1)*(-b+r3+r1)*(b-r3+r1)*(b+r3-r1);
        rg_REAL t3_squared = (c+r1+r2)*(-c+r1+r2)*(c-r1+r2)*(c+r1-r2);
        
        if( t1_squared >= 0.0 && t2_squared >= 0.0 && t3_squared >= 0.0 )
        {
            rg_REAL t1 = sqrt( t1_squared );
            rg_REAL t2 = sqrt( t2_squared );
            rg_REAL t3 = sqrt( t3_squared );
        
            rg_REAL p1_plus  = ( ( b_squared - c_squared + r2_squared - r3_squared )
					           *( b_squared - c_squared + r2_squared - r3_squared ) 
                               + ( t + t1 )*( t + t1 ) )
				              /(4.0*a_squared) - r1_squared;
            rg_REAL p1_minus = ( ( b_squared - c_squared + r2_squared - r3_squared )
					           *( b_squared - c_squared + r2_squared - r3_squared ) 
                               + ( t - t1 )*( t - t1 ) )
				              /(4.0*a_squared) - r1_squared;
            rg_REAL p2_plus  = ( ( c_squared - a_squared + r3_squared - r1_squared )
					           *( c_squared - a_squared + r3_squared - r1_squared ) 
                               + ( t + t2 )*( t + t2 ) )
				              /(4.0*b_squared) - r2_squared;
            rg_REAL p2_minus = ( ( c_squared - a_squared + r3_squared - r1_squared )
					           *( c_squared - a_squared + r3_squared - r1_squared ) 
                               + ( t - t2 )*( t - t2 ) )
				              /(4.0*b_squared) - r2_squared;
            rg_REAL p3_plus  = ( ( a_squared - b_squared + r1_squared - r2_squared )
					           *( a_squared - b_squared + r1_squared - r2_squared ) 
                               + ( t + t3 )*( t + t3 ) )
				              /(4.0*c_squared) - r3_squared;
            rg_REAL p3_minus = ( ( a_squared - b_squared + r1_squared - r2_squared )
					           *( a_squared - b_squared + r1_squared - r2_squared ) 
                               + ( t - t3 )*( t - t3 ) )
				              /(4.0*c_squared) - r3_squared;

            // the 1st excetpoin case
            if(       ( p1_plus>0.0 && p1_minus>0.0 )  &&  
				      ( p2_plus>0.0 && p2_minus>0.0 )  &&  
				      ( p3_plus>0.0 && p3_minus>0.0 )  )        
            {
                volumeIntersectedByThreeSpheres = 0.0;
            }
            // the 2nd exception case
            else if(  ( p1_plus<0.0 && p1_minus<0.0 )  &&  
					  ( p2_plus>0.0 && p2_minus>0.0 )  &&  
					  ( p3_plus>0.0 && p3_minus>0.0 )  )   
            {
                volumeIntersectedByThreeSpheres = computeIntersectingVolume( sphere2 );
            }
            else if(  ( p1_plus>0.0 && p1_minus>0.0 )  &&  
					  ( p2_plus<0.0 && p2_minus<0.0 )  &&  
					  ( p3_plus>0.0 && p3_minus>0.0 )  )
            {
                volumeIntersectedByThreeSpheres = computeIntersectingVolume( sphere1 );
            }
            else if(  ( p1_plus>0.0 && p1_minus>0.0 )  &&  
					  ( p2_plus>0.0 && p2_minus>0.0 )  &&  
					  ( p3_plus<0.0 && p3_minus<0.0 )  )
            {
                volumeIntersectedByThreeSpheres = sphere1.computeIntersectingVolume( sphere2 );
            }
            // the 3rd exception case
            else if(  ( p1_plus>0.0 && p1_minus>0.0 )  &&  
					  ( p2_plus<0.0 && p2_minus<0.0 )  &&  
					  ( p3_plus<0.0 && p3_minus<0.0 )  )   
            {
                volumeIntersectedByThreeSpheres = sphere1.computeIntersectingVolume( sphere2 )
                                                + computeIntersectingVolume( sphere1 )
                                                - sphere1.computeVolume();
            }
            else if(  ( p1_plus<0.0 && p1_minus<0.0 )  &&  
					  ( p2_plus>0.0 && p2_minus>0.0 )  &&  
					  ( p3_plus<0.0 && p3_minus<0.0 )  )
            {
                volumeIntersectedByThreeSpheres = sphere1.computeIntersectingVolume( sphere2 )
                                                + computeIntersectingVolume( sphere2 )
                                                - sphere2.computeVolume();
            }
            else if(  ( p1_plus<0.0 && p1_minus<0.0 )  &&  
					  ( p2_plus<0.0 && p2_minus<0.0 )  &&  
					  ( p3_plus>0.0 && p3_minus>0.0 )  )
            {
                volumeIntersectedByThreeSpheres = computeIntersectingVolume( sphere1 )
                                                + computeIntersectingVolume( sphere2 )
                                                - computeVolume();
            }
        }
        else    //!(t1_squared >= 0.0 && t2_squared >= 0.0 && t3_squared >= 0.0)
        {
            // 구 1은 구 2와 3의 rg_REAL intersection안에 포함된다.
            if( t1_squared >= 0.0 && t2_squared < 0.0 && t3_squared < 0.0 )
            {
                volumeIntersectedByThreeSpheres = sphere1.computeVolume();
            }
            // 구 2은 구 1과 3의 rg_REAL intersection안에 포함된다.
            else if( t1_squared < 0.0 && t2_squared >= 0.0 && t3_squared < 0.0 )
            {
                volumeIntersectedByThreeSpheres = sphere2.computeVolume();
            }
            // 구 3은 구 1과 2의 rg_REAL intersection안에 포함된다.
            else if( t1_squared < 0.0 && t2_squared < 0.0 && t3_squared >= 0.0)
            {
                volumeIntersectedByThreeSpheres = computeVolume();
            }
            // 구 1는 구 2와 3 사이에서 서로 rg_REAL intersection하고 있다. 구 2와 3은 겹치지 않는다.
            else if( t1_squared < 0.0 && t2_squared >= 0.0 && t3_squared >= 0.0 )
            {
                volumeIntersectedByThreeSpheres = 0.0;
            }
            // 구 2는 구 1과 3 사이에서 서로 rg_REAL intersection하고 있다. 구 1과 3은 겹치지 않는다.
            else if( t1_squared >= 0.0 && t2_squared < 0.0 && t3_squared >= 0.0 )
            {
                volumeIntersectedByThreeSpheres = 0.0;
            }
            // 구 3는 구 1과 2 사이에서 서로 rg_REAL intersection하고 있다. 구 1과 2은 겹치지 않는다.  
            else if( t1_squared >= 0.0 && t2_squared >= 0.0 && t3_squared < 0.0 )
            {
                volumeIntersectedByThreeSpheres = 0.0;
            }
            // 구 세개가 다 떨어져 있다.
            else    //( t1_squared < 0.0 && t2_squared < 0.0 && t3_squared < 0.0 )  
            {
                volumeIntersectedByThreeSpheres = 0.0;
            }            
        }
    }
    else    //( w_squared == 0.0 ) sphere 1, 2 and 3 are intersected at SINGLE POINT simultaneously
    {
        volumeIntersectedByThreeSpheres = 0.0;
    }

    return volumeIntersectedByThreeSpheres;	
}

rg_REAL Sphere::computeSumPairwiseXVolOfBallWithBallSet(rg_dList<Sphere>& ballSet) const
{
    rg_REAL xVolume = 0.0;	
    ballSet.reset4Loop();
    while (ballSet.setNext4Loop())
    {
        Sphere currBall = ballSet.getEntity();
        //xVolume += ball1.computeIntersectingVolume( currBall );
        xVolume += rg_GeoFunc::computeXVolumeOfTwoSpheres(*this, currBall);
    }
    return xVolume;
}

rg_REAL Sphere::computeIntersectionVolumeWithFourSpheresOnlyIfThereIsCommonIntersection(const Sphere& sphere1, const Sphere& sphere2, const Sphere& sphere3) const
{
    rg_REAL volumeIntersectedByFourSpheres  = 0.0;

    rg_REAL r1 = sphere1.m_radius; 
    rg_REAL r2 = sphere2.m_radius;
    rg_REAL r3 = sphere3.m_radius;
    rg_REAL r4 = m_radius;

    rg_REAL a  = sphere2.m_center.distance(sphere3.m_center);
    rg_REAL b  = sphere3.m_center.distance(sphere1.m_center);
    rg_REAL c  = sphere1.m_center.distance(sphere2.m_center);
    rg_REAL f  = sphere1.m_center.distance(m_center);
    rg_REAL g  = sphere2.m_center.distance(m_center);
    rg_REAL h  = sphere3.m_center.distance(m_center);


    rg_REAL r1_squared   = r1*r1;
    rg_REAL r2_squared   = r2*r2;
    rg_REAL r3_squared   = r3*r3;
    rg_REAL r4_squared   = r4*r4;
    
    rg_REAL a_squared    = a*a;
    rg_REAL b_squared    = b*b;
    rg_REAL c_squared    = c*c;
    rg_REAL f_squared    = f*f;
    rg_REAL g_squared    = g*g;
    rg_REAL h_squared    = h*h;

    rg_REAL w_squared = computeW_SquaredForQuadrupleIntersection(a_squared, b_squared, c_squared,
		                                                         f_squared, g_squared, h_squared,
																 r1_squared, r2_squared, r3_squared, r4_squared);

   rg_REAL w = 0.0;

	if( w_squared >= 0 )
    {
        w = sqrt(w_squared);
    }
    else
    {
        return volumeIntersectedByFourSpheres;
    }

    rg_REAL w1_squared   = computeW_SquaredForTripleIntersection( a_squared, b_squared, c_squared, r1_squared, r2_squared, r3_squared );//123(23,31,12)
    rg_REAL w2_squared   = computeW_SquaredForTripleIntersection( g_squared, f_squared, c_squared, r1_squared, r2_squared, r4_squared );//124(24,41,12)
    rg_REAL w3_squared   = computeW_SquaredForTripleIntersection( h_squared, f_squared, b_squared, r1_squared, r3_squared, r4_squared );//134(34,41,13)
    rg_REAL w4_squared   = computeW_SquaredForTripleIntersection( h_squared, g_squared, a_squared, r2_squared, r3_squared, r4_squared );//234(34,42,23)
    
    // sphere 1, 2, 3 and 4 are intersected simultaneously
    if( w1_squared > 0.0 && w2_squared > 0.0 && w3_squared > 0.0 && w4_squared > 0.0 )
    {
        
    
        rg_REAL s1   = a*( b_squared + c_squared - a_squared + g_squared + h_squared - 2*f_squared + 
			             (g_squared-h_squared)*(b_squared-c_squared)/a_squared );
        rg_REAL s2   = b*( c_squared + a_squared - b_squared + h_squared + f_squared - 2*g_squared +
		                 (h_squared-f_squared)*(c_squared-a_squared)/b_squared );
        rg_REAL s3   = c*( a_squared + b_squared - c_squared + f_squared + g_squared - 2*h_squared + 
			             (f_squared-g_squared)*(a_squared-b_squared)/c_squared );
        rg_REAL s4   = h*( f_squared + b_squared - h_squared + a_squared + g_squared - 2*c_squared + 
			             (f_squared-b_squared)*(a_squared-g_squared)/h_squared ); 
        rg_REAL s5   = g*( h_squared + a_squared - g_squared + c_squared + f_squared - 2*b_squared + 
			             (c_squared-f_squared)*(h_squared-a_squared)/g_squared );
        rg_REAL s6   = f*( g_squared + c_squared - f_squared + b_squared + h_squared - 2*a_squared + 
			             (b_squared-h_squared)*(g_squared-c_squared)/f_squared );


        volumeIntersectedByFourSpheres = 
            - w/12.0
            - rg_PI*( r1_squared*r1 + r2_squared*r2 + r3_squared*r3 + r4_squared*r4 )/3.0
            + ( (r2_squared+r3_squared)*a + (r2_squared-r3_squared)*(r2_squared-r3_squared)/(2.0*a)
			   	    - a_squared*a/6.0 )*atan2( 1.0, s1 / (2.0*w) )/4.0
            + ( (r1_squared+r3_squared)*b + (r1_squared-r3_squared)*(r1_squared-r3_squared)/(2.0*b)
			   	    - b_squared*b/6.0 )*atan2( 1.0, s2 / (2.0*w) )/4.0
            + ( (r1_squared+r2_squared)*c + (r1_squared-r2_squared)*(r1_squared-r2_squared)/(2.0*c)
			   	    - c_squared*c/6.0 )*atan2( 1.0, s3 / (2.0*w) )/4.0
            + ( (r3_squared+r4_squared)*h + (r3_squared-r4_squared)*(r3_squared-r4_squared)/(2.0*h)
			   	    - h_squared*h/6.0 )*atan2( 1.0, s4 / (2.0*w) )/4.0
            + ( (r2_squared+r4_squared)*g + (r2_squared-r4_squared)*(r2_squared-r4_squared)/(2.0*g)
			   	    - g_squared*g/6.0 )*atan2( 1.0, s5 / (2.0*w) )/4.0
            + ( (r1_squared+r4_squared)*f + (r1_squared-r4_squared)*(r1_squared-r4_squared)/(2.0*f)
			   	    - f_squared*f/6.0 )*atan2( 1.0, s6 / (2.0*w) )/4.0
            + sphere3.computeIntersectionVolumeWith(sphere1,sphere2)/2.0
            + computeIntersectionVolumeWith(sphere1,sphere2)/2.0
            + computeIntersectionVolumeWith(sphere1,sphere3)/2.0
            + computeIntersectionVolumeWith(sphere2,sphere3)/2.0;
    }
    // sphere 1, 2, 3 and 4 are NOT intersected simultaneously ( no quadruple intersection )
    else if( w1_squared < 0.0 && w2_squared < 0.0 && w3_squared < 0.0 && w4_squared < 0.0 )
    {
        volumeIntersectedByFourSpheres = 0.0;
    }
    else    // ( exception cases )    
    {
        volumeIntersectedByFourSpheres = 0.0;
    }

    return volumeIntersectedByFourSpheres;
}

rg_REAL Sphere::computeW_SquaredForTripleIntersection( const rg_REAL& a_squared,  const rg_REAL& b_squared,  const rg_REAL& c_squared,
                                                       const rg_REAL& r1_squared, const rg_REAL& r2_squared, const rg_REAL& r3_squared ) const
{
    rg_REAL w_squared = 

        ((r1_squared-b_squared)*b_squared-r2_squared*b_squared)*r2_squared
      + ( r2_squared*b_squared-b_squared*r1_squared)*r3_squared
      + (-r1_squared*r1_squared+b_squared*r1_squared+(r1_squared+b_squared)
        *r2_squared+(r1_squared-r2_squared)*r3_squared-a_squared*r1_squared)*a_squared
      + ((-r1_squared+b_squared)*r2_squared+(r2_squared+r1_squared+b_squared-r3_squared)
        *r3_squared+(r3_squared-b_squared+r1_squared)*a_squared-r3_squared*c_squared)*c_squared;

    return w_squared;
}

rg_REAL Sphere::computeW_SquaredForQuadrupleIntersection( const rg_REAL& a_squared,  const rg_REAL& b_squared,  const rg_REAL& c_squared,
			                                              const rg_REAL& f_squared,  const rg_REAL& g_squared,  const rg_REAL& h_squared,
														  const rg_REAL& r1_squared, const rg_REAL& r2_squared, const rg_REAL& r3_squared, const rg_REAL& r4_squared) const
{
    rg_REAL w_squared =   
        - a_squared*g_squared*h_squared + b_squared*a_squared*f_squared 
        + a_squared*g_squared*b_squared + b_squared*g_squared*f_squared 
        + c_squared*b_squared*h_squared - b_squared*h_squared*f_squared 
        + c_squared*h_squared*f_squared + g_squared*a_squared*f_squared 
        + c_squared*g_squared*h_squared + c_squared*a_squared*f_squared 
        + c_squared*g_squared*b_squared - c_squared*g_squared*f_squared 
        + a_squared*h_squared*f_squared + c_squared*a_squared*h_squared 
        + g_squared*b_squared*h_squared - c_squared*a_squared*b_squared
        - c_squared*h_squared*h_squared - a_squared*a_squared*f_squared 
        - g_squared*g_squared*b_squared - c_squared*c_squared*h_squared 
        - g_squared*b_squared*b_squared - a_squared*f_squared*f_squared;

	return w_squared;
}

FourSpheresIntersectionType parseTypeOfIntersectionAmongFourSpheres(const Sphere& sphere1, 
                                                                    const Sphere& sphere2, 
                                                                    const Sphere& sphere3, 
                                                                    const Sphere& sphere4 )
{	
	rg_INT     intersectionPointInclusion[ 4 ][ 2 ] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};
	rg_Point3D intersectionPoint[ 4 ][ 2 ];
	rg_INT     numIntersectionPoint[ 4 ];
	doContainmentTestOfTripletSpheresIntersectionForFourthSphere(sphere1, 
		                                                         sphere2, 
																 sphere3, 
																 sphere4, 
																 intersectionPointInclusion,
																 intersectionPoint,
																 numIntersectionPoint);

	rg_INT intersectionCount[ 4 ] = {0, 0, 0, 0};
	rg_INT i;
	for(i = 0;i < 4;i++)
		intersectionCount[i] = intersectionPointInclusion[i][0] + intersectionPointInclusion[i][1];

	FourSpheresIntersectionType type;

	if(intersectionCount[0] == 0 && intersectionCount[1] == 0 && intersectionCount[2] == 0 && intersectionCount[3] == 0)
	{
		type = SI_QUADRUPLE_NO_INTERSECTION;
	}
	else if(intersectionCount[0] == 1 && intersectionCount[1] == 1 && intersectionCount[2] == 1 && intersectionCount[3] == 1) 
	{
		type = SI_QUADRUPLE_COMMON_INTERSECTION;
	}
	else
	{
		type = SI_QUADRUPLE_COMMON_INTERSECTION_YET_IGNORABLE;
	}

	return type;
}

void doContainmentTestOfTripletSpheresIntersectionForFourthSphere(const Sphere& sphere1, 
                                                                  const Sphere& sphere2, 
                                                                  const Sphere& sphere3, 
                                                                  const Sphere& sphere4,
													              rg_INT intersectionPointInclusion[][2],
																  rg_Point3D intersectionPoint[][2],
																  rg_INT numIntersectionPoint[])
{
	Sphere sphere[ 4 ] = {sphere1, sphere2, sphere3, sphere4};
	rg_INDEX indexForTripletSpheres[ 4 ][ 3 ] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
	rg_INDEX indexForFourthSphere[ 4 ]        = {3, 2, 1, 0};	

	rg_Point3D* intersection = new rg_Point3D[ 2 ];
	rg_INT      numIntersections = 0;
	rg_INT i;
	for (i = 0;i < 4;i++)
	{
		numIntersections = computeIntersectionPointsAmong(sphere[indexForTripletSpheres[i][0]],
			                                              sphere[indexForTripletSpheres[i][1]],
									                      sphere[indexForTripletSpheres[i][2]],
														  intersection);
		numIntersectionPoint[i] = numIntersections;
		intersectionPoint[i][0] = intersection[0];
		intersectionPoint[i][1] = intersection[1];

		if(sphere[indexForFourthSphere[i]].doesContain(intersection[0]))
			intersectionPointInclusion[i][0]++;

		if(sphere[indexForFourthSphere[i]].doesContain(intersection[1]))
			intersectionPointInclusion[i][1]++;
	}
	delete [] intersection;
}


bool    Sphere::hasIntersectionWith(const LineSegment3D& lineSegment) const
{
    rg_Point3D          footprint;
    rg_REAL             param;
    PosOnLineSegOrArc   pos;

    rg_REAL distance = lineSegment.computeMinDistFromPoint(m_center, footprint, param, pos);

    if (pos == NONEXTREME_PT && distance < m_radius) {
        return true;
    }
    else {
        return false;
    }
}



rg_FLAG Sphere::isThereIntersectionWith(const Sphere& anotherSphere) const
{
    rg_REAL distanceOfTwoCenters = m_center.distance(anotherSphere.m_center);
    rg_REAL sumOfTwoRadii        = m_radius + anotherSphere.m_radius;
    if ( rg_LT(distanceOfTwoCenters, sumOfTwoRadii ) ) 
        return rg_TRUE;
    else 
        return rg_FALSE;
}

rg_BOOL Sphere::doesContainIntersectionBTW(const Sphere& sphere1, const Sphere& sphere2) const
{
	// distance between the 2nd sphere and 3rd sphere
    rg_REAL a    = sphere2.m_center.distance(m_center);  
    // distance between the 3rd sphere and 1st sphere
    rg_REAL b    = m_center.distance(sphere1.m_center);   
    // distance between the 1st sphere and 2nd sphere
    rg_REAL c    = sphere1.m_center.distance(sphere2.m_center);   

    rg_REAL r1   = sphere1.m_radius;  // r1   Radius of the 1st sphere
    rg_REAL r2   = sphere2.m_radius;  // r2   Radius of the 2nd sphere
    rg_REAL r3   = m_radius;  // r3   Radius of the 3rd sphere
    
    rg_REAL a_squared   = a*a;  // square of a
    rg_REAL b_squared   = b*b;  // square of b
    rg_REAL c_squared   = c*c;  // square of c

    rg_REAL r1_squared  = r1*r1;    // rr1  Square of r1, or r1*r1
    rg_REAL r2_squared  = r2*r2;    // rr2  Square of r2, or r2*r2
    rg_REAL r3_squared  = r3*r3;    // rr3  Square of r3, or r3*r3

	rg_REAL t_squared = (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c);
	rg_REAL t1_squared = (a+r2+r3)*(-a+r2+r3)*(a-r2+r3)*(a+r2-r3);
	rg_REAL t2_squared = (b+r3+r1)*(-b+r3+r1)*(b-r3+r1)*(b+r3-r1);
	rg_REAL t3_squared = (c+r1+r2)*(-c+r1+r2)*(c-r1+r2)*(c+r1-r2);

	if(t_squared < 0 || t1_squared < 0 || t2_squared < 0 || t3_squared <0)
		return rg_FALSE;

	rg_REAL t  = sqrt( t_squared );
    rg_REAL t1 = sqrt( t1_squared );
    rg_REAL t2 = sqrt( t2_squared );
    rg_REAL t3 = sqrt( t3_squared );

    rg_REAL p1_plus  = ( ( b_squared - c_squared + r2_squared - r3_squared )
					   *( b_squared - c_squared + r2_squared - r3_squared ) 
                       + ( t + t1 )*( t + t1 ) )
				      /(4.0*a_squared) - r1_squared;
    rg_REAL p1_minus = ( ( b_squared - c_squared + r2_squared - r3_squared )
					   *( b_squared - c_squared + r2_squared - r3_squared ) 
                       + ( t - t1 )*( t - t1 ) )
				      /(4.0*a_squared) - r1_squared;
    rg_REAL p2_plus  = ( ( c_squared - a_squared + r3_squared - r1_squared )
					   *( c_squared - a_squared + r3_squared - r1_squared ) 
                       + ( t + t2 )*( t + t2 ) )
				      /(4.0*b_squared) - r2_squared;
    rg_REAL p2_minus = ( ( c_squared - a_squared + r3_squared - r1_squared )
					   *( c_squared - a_squared + r3_squared - r1_squared ) 
                       + ( t - t2 )*( t - t2 ) )
				      /(4.0*b_squared) - r2_squared;
    rg_REAL p3_plus  = ( ( a_squared - b_squared + r1_squared - r2_squared )
					   *( a_squared - b_squared + r1_squared - r2_squared ) 
                       + ( t + t3 )*( t + t3 ) )
				      /(4.0*c_squared) - r3_squared;
    rg_REAL p3_minus = ( ( a_squared - b_squared + r1_squared - r2_squared )
					   *( a_squared - b_squared + r1_squared - r2_squared ) 
                       + ( t - t3 )*( t - t3 ) )
				      /(4.0*c_squared) - r3_squared;

	if(  ( p1_plus>0.0 && p1_minus>0.0 )  &&  
	     ( p2_plus>0.0 && p2_minus>0.0 )  &&  
		 ( p3_plus<0.0 && p3_minus<0.0 )     )
	{
		return rg_TRUE;	
	}
	else
		return rg_FALSE;
}

rg_BOOL Sphere::doesContainIntersectionAmong(const Sphere& sphere1, const Sphere& sphere2, const Sphere& sphere3) const
{
	Sphere   sphere[ 4 ]     = {sphere1, sphere2, sphere3, *this};
	rg_INDEX indexForTripletSpheres[ 3 ][ 3 ] = {{0,1,3}, {0,2,3}, {1,2,3}};
	rg_INDEX indexOfFourthSphere[ 3 ]         = {2,       1,       0      };
	rg_Point3D* intersections = new rg_Point3D[ 2 ];
	rg_INT numOfIntersections =
		computeIntersectionPointsAmong(sphere[0], 
			                           sphere[1], 
									   sphere[2], intersections);

	rg_BOOL bOnePairIncluded = rg_TRUE;
	rg_BOOL bIncluded        = rg_TRUE;

	if(numOfIntersections != 2 ||
	   doesContain(intersections[ 0 ]) == rg_FALSE ||
	   doesContain(intersections[ 1 ]) == rg_FALSE   )
	   bOnePairIncluded = rg_FALSE;

	if(!bOnePairIncluded)
		bIncluded = rg_FALSE;
	else
	{
		rg_INT i;
		for(i = 0;i < 3;i++)
		{
			numOfIntersections =
				computeIntersectionPointsAmong(sphere[indexForTripletSpheres[i][0]], 
											   sphere[indexForTripletSpheres[i][1]], 
											   sphere[indexForTripletSpheres[i][2]], intersections);

			if(numOfIntersections != 2 || 
			   sphere[indexOfFourthSphere[i]].doesContain(intersections[ 0 ]) == rg_TRUE ||
			   sphere[indexOfFourthSphere[i]].doesContain(intersections[ 1 ]) == rg_TRUE   )
			{
				bIncluded = rg_FALSE;
				break;
			}
		}
	}

	delete [] intersections;

	return bIncluded;
}

rg_BOOL Sphere::doesContain(const rg_Point3D& point) const
{
    rg_REAL distance = m_center.distance(point);

    if ( rg_LE( distance, m_radius ) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}

rg_BOOL Sphere::isOnSphere(const rg_Point3D& point) const
{
    rg_REAL distance = m_center.distance(point);
	if(rg_EQ(distance, m_radius))
		return rg_TRUE;
	else
		return rg_FALSE;
}


rg_BOOL Sphere::isIncludedIn(const Sphere& anotherSphere) const
{
    rg_REAL distance = m_center.distance(anotherSphere.m_center);

    if ( (distance+m_radius) <= anotherSphere.m_radius )
        return rg_TRUE;
    else
        return rg_FALSE;
}


rg_FLAG Sphere::isContainedIn(const Sphere& anotherSphere) const
{
    rg_REAL distance = m_center.distance(anotherSphere.m_center);

    if ( (distance+m_radius) <= anotherSphere.m_radius )
        return rg_TRUE;
    else if ( (distance+anotherSphere.m_radius) <= m_radius )
        return rg_TRUE+rg_TRUE;
    else
        return rg_FALSE;
}


rg_BOOL Sphere::isIncludedInUnionOfSpheres(const Sphere& sphere1, const Sphere& sphere2) const
{
	// distance between the 2nd sphere and 3rd sphere
    rg_REAL a    = sphere2.m_center.distance(m_center);  
    // distance between the 3rd sphere and 1st sphere
    rg_REAL b    = m_center.distance(sphere1.m_center);   
    // distance between the 1st sphere and 2nd sphere
    rg_REAL c    = sphere1.m_center.distance(sphere2.m_center);   

    rg_REAL r1   = sphere1.m_radius;  // r1   Radius of the 1st sphere
    rg_REAL r2   = sphere2.m_radius;  // r2   Radius of the 2nd sphere
    rg_REAL r3   = m_radius;  // r3   Radius of the 3rd sphere
    
    rg_REAL a_squared   = a*a;  // square of a
    rg_REAL b_squared   = b*b;  // square of b
    rg_REAL c_squared   = c*c;  // square of c

    rg_REAL r1_squared  = r1*r1;    // rr1  Square of r1, or r1*r1
    rg_REAL r2_squared  = r2*r2;    // rr2  Square of r2, or r2*r2
    rg_REAL r3_squared  = r3*r3;    // rr3  Square of r3, or r3*r3

	rg_REAL t_squared = (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c);
	rg_REAL t1_squared = (a+r2+r3)*(-a+r2+r3)*(a-r2+r3)*(a+r2-r3);
	rg_REAL t2_squared = (b+r3+r1)*(-b+r3+r1)*(b-r3+r1)*(b+r3-r1);
	rg_REAL t3_squared = (c+r1+r2)*(-c+r1+r2)*(c-r1+r2)*(c+r1-r2);

	if(t_squared < 0 || t1_squared < 0 || t2_squared < 0 || t3_squared <0)
		return rg_FALSE;

	rg_REAL t  = sqrt( t_squared );
    rg_REAL t1 = sqrt( t1_squared );
    rg_REAL t2 = sqrt( t2_squared );
    rg_REAL t3 = sqrt( t3_squared );

    rg_REAL p1_plus  = ( ( b_squared - c_squared + r2_squared - r3_squared )
					   *( b_squared - c_squared + r2_squared - r3_squared ) 
                       + ( t + t1 )*( t + t1 ) )
				      /(4.0*a_squared) - r1_squared;
    rg_REAL p1_minus = ( ( b_squared - c_squared + r2_squared - r3_squared )
					   *( b_squared - c_squared + r2_squared - r3_squared ) 
                       + ( t - t1 )*( t - t1 ) )
				      /(4.0*a_squared) - r1_squared;
    rg_REAL p2_plus  = ( ( c_squared - a_squared + r3_squared - r1_squared )
					   *( c_squared - a_squared + r3_squared - r1_squared ) 
                       + ( t + t2 )*( t + t2 ) )
				      /(4.0*b_squared) - r2_squared;
    rg_REAL p2_minus = ( ( c_squared - a_squared + r3_squared - r1_squared )
					   *( c_squared - a_squared + r3_squared - r1_squared ) 
                       + ( t - t2 )*( t - t2 ) )
				      /(4.0*b_squared) - r2_squared;
    rg_REAL p3_plus  = ( ( a_squared - b_squared + r1_squared - r2_squared )
					   *( a_squared - b_squared + r1_squared - r2_squared ) 
                       + ( t + t3 )*( t + t3 ) )
				      /(4.0*c_squared) - r3_squared;
    rg_REAL p3_minus = ( ( a_squared - b_squared + r1_squared - r2_squared )
					   *( a_squared - b_squared + r1_squared - r2_squared ) 
                       + ( t - t3 )*( t - t3 ) )
				      /(4.0*c_squared) - r3_squared;

	if(  ( p1_plus<0.0 && p1_minus<0.0 )  &&  
		 ( p2_plus<0.0 && p2_minus<0.0 )  &&  
		 ( p3_plus>0.0 && p3_minus>0.0 )  )
		 return rg_TRUE;
	else
		return rg_FALSE;
}

rg_BOOL Sphere::isIncludedInUnionOfSpheres(const Sphere& sphere1, 
										   const Sphere& sphere2, 										   
										   const Sphere& sphere3) const
{
	Sphere   sphere[ 4 ]     = {sphere1, sphere2, sphere3, *this};

	rg_INT     intersectionPointInclusion[ 4 ][ 2 ] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};
	rg_Point3D intersectionPoint[ 4 ][ 2 ];
	rg_INT     numIntersectionPoint[ 4 ];
	doContainmentTestOfTripletSpheresIntersectionForFourthSphere(sphere[0], 
		                                                         sphere[1], 
																 sphere[2], 
																 sphere[3], 
																 intersectionPointInclusion,
																 intersectionPoint,
																 numIntersectionPoint);

	rg_INT intersectionCount[ 4 ] = {0, 0, 0, 0};
	rg_INT numPairIncluded = 0;
	rg_INT i;
	for(i = 0;i < 4;i++)
	{
		intersectionCount[i] = intersectionPointInclusion[i][0] + intersectionPointInclusion[i][1];
		if(intersectionCount[i] == 2)
			numPairIncluded++;
	}	
	
	rg_BOOL bIncluded = rg_TRUE;

	if(numPairIncluded == 3)
	{
		rg_INDEX indexForTripletSpheres[ 3 ][ 3 ] = {{0,1,3}, {0,2,3}, {1,2,3}};
		rg_INDEX indexForFourthSphere[ 3 ]        = {2,        1,       0     };
		rg_Point3D* intersections = new rg_Point3D[ 2 ];
		rg_INT numOfIntersections = 0;
		rg_BOOL bThreePairsIncluded = rg_TRUE;
		
		rg_INT i;
		for(i = 0;i < 3;i++)
		{
			numOfIntersections =
				computeIntersectionPointsAmong(sphere[indexForTripletSpheres[i][0]], 
											   sphere[indexForTripletSpheres[i][1]], 
											   sphere[indexForTripletSpheres[i][2]], intersections);

			if(numOfIntersections != 2 || 
			   sphere[indexForFourthSphere[i]].doesContain(intersections[ 0 ]) == rg_FALSE ||
			   sphere[indexForFourthSphere[i]].doesContain(intersections[ 1 ]) == rg_FALSE   )
			{
				bThreePairsIncluded = rg_FALSE;
				break;
			}
		}

		if(bThreePairsIncluded)
		{
			numOfIntersections =
				computeIntersectionPointsAmong(sphere[0], sphere[1], sphere[2], intersections);
			if(numOfIntersections != 2 ||
			   sphere[3].doesContain(intersections[0]) == rg_TRUE ||
			   sphere[3].doesContain(intersections[1]) == rg_TRUE   )
			{
			   bIncluded = rg_FALSE;
			}
			else
			{
				bIncluded = rg_TRUE;
			}
		}
		else
			bIncluded = rg_FALSE;

		delete [] intersections;

	}
	else if(numPairIncluded == 4)
	{
		rg_INDEX indexForTripletSpheres[ 4 ][ 3 ] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
		rg_INDEX indexForFourthSphere[ 4 ]        = {3, 2, 1, 0};	

		rg_Point3D intersection[2] = {intersectionPoint[0][0], intersectionPoint[0][1]};
		for(i = 0;i < 4;i++)
		{
			if(!sphere[indexForFourthSphere[i]].isOnSphere(intersectionPoint[i][0]) ||
			   !sphere[indexForFourthSphere[i]].isOnSphere(intersectionPoint[i][1]) ||
			   intersection[0] != intersectionPoint[i][0]                           ||
			   intersection[1] != intersectionPoint[i][1]                           )
			{
				bIncluded = rg_FALSE;
				break;
			}
		}
		
		Plane plane;
		plane.definePlaneByThreePoints(sphere[0].getCenter(), sphere[1].getCenter(), sphere[2].getCenter());
		if(plane.isOnThisPlane(sphere[3].getCenter()))
		{
			rg_REAL areaOfTriangle[ 4 ] = {0., 0., 0., 0.};
			for(i = 0;i < 4;i++)
			{
				areaOfTriangle[ i ] = 
                    rg_GeoFunc::computeTriangleArea(sphere[indexForTripletSpheres[i][0]].getCenter(), 
										 sphere[indexForTripletSpheres[i][1]].getCenter(), 
										 sphere[indexForTripletSpheres[i][2]].getCenter());
			}

			if(rg_EQ(areaOfTriangle[0], areaOfTriangle[1]+areaOfTriangle[2]+areaOfTriangle[3]))
			{
				bIncluded = rg_TRUE;
			}
			else
			{
				bIncluded = rg_FALSE;
			}
		}
		else
			bIncluded = rg_FALSE;
	}
	else
	{
		bIncluded = rg_FALSE;
	}

	return bIncluded;
}

void Sphere::makeContainingSphereWith(const Sphere& anotherSphere) 
{
    rg_FLAG isContained = isContainedIn(anotherSphere);
    if ( isContained == rg_TRUE )  {
    }
    else if ( isContained == rg_TRUE+rg_TRUE)  {
        m_center = anotherSphere.m_center;
        m_radius = anotherSphere.m_radius;
    }
    else  {
        rg_Point3D vecThisToAnother = anotherSphere.m_center - m_center;
        vecThisToAnother.normalize();


        rg_Point3D farPointOfThis    = m_center + (-m_radius*vecThisToAnother);
        rg_Point3D farPointOfAnother = anotherSphere.m_center + (anotherSphere.m_radius*vecThisToAnother);

        m_center = ( farPointOfThis + farPointOfAnother )/2.0;

        m_radius = farPointOfThis.distance(farPointOfAnother)/2.0;
    }
}



Sphere Sphere::computeContainingSphereOfThisAnd(const Sphere& anotherSphere) const
{
    rg_REAL distance = m_center.distance(anotherSphere.m_center);

    if ( rg_LE( distance+anotherSphere.m_radius, m_radius) )
        return *this;
    else if ( rg_LE( distance+m_radius, anotherSphere.m_radius) )
        return anotherSphere;
    else
    {
        rg_Point3D leftMostPoint;
        rg_Point3D rightMostPoint;

        rg_Point3D vectorFromLMPoRMP = (m_center - anotherSphere.m_center)/distance;

        leftMostPoint = m_center + (m_radius*vectorFromLMPoRMP);
        leftMostPoint = anotherSphere.m_center - (anotherSphere.m_radius*vectorFromLMPoRMP);

        rg_Point3D center = ( leftMostPoint + rightMostPoint )/2.0;

        rg_REAL radius = distance + m_radius + anotherSphere.m_radius;

        return Sphere(center, radius);
    }
}



rg_INT Sphere::intersect(const LineSegment3D& lineSegment, rg_Point3D* intersectionPoint) const
{
    rg_INT numIntersection = 0;

    rg_REAL    param = 0.0;
    rg_Point3D projPtOnLine = lineSegment.computeProjectedPointForPoint3D(m_center, param);

    if ( m_center != projPtOnLine ) {
        rg_Point3D normal = projPtOnLine - m_center;
        normal.normalize();

        Plane planeContainingLine(normal, projPtOnLine);

        Circle3D intersectingCircle;
        rg_BOOL isIntersect = intersect(planeContainingLine, intersectingCircle);

        if ( isIntersect ) {
            numIntersection = intersectingCircle.intersect(lineSegment, intersectionPoint);
        }
    }
    else {
        rg_Point3D startPt    = lineSegment.getStartPt();
        rg_Point3D endPt      = lineSegment.getEndPt();

        rg_FLAG isStartPtInsideSphere = doesContain( startPt );
        rg_FLAG isEndPtInsideSphere   = doesContain( endPt );

        if ( !isStartPtInsideSphere ) {
            rg_Point3D vectorToStart = startPt - m_center;
            vectorToStart.normalize();

            rg_Point3D intersection = m_center + (m_radius*vectorToStart);

            intersectionPoint[numIntersection] = intersection;
            numIntersection++;
        }

        if ( !isEndPtInsideSphere ) {
            rg_Point3D vectorToEnd = endPt - m_center;
            vectorToEnd.normalize();

            rg_Point3D intersection = m_center + (m_radius*vectorToEnd);

            intersectionPoint[numIntersection] = intersection;
            numIntersection++;
        }
    }

    return numIntersection;
}


rg_INT Sphere::intersect(const Circle3D& circle, rg_Point3D* intersectionPoint) const
{
    Plane planeContainingCircle = circle.getPlaneContainingThisCircle();

    Circle3D intersectCircle;
    rg_BOOL bIntersected = intersect(planeContainingCircle, intersectCircle);

	rg_INT numIntersections = 0;
	if(bIntersected)
	{
		numIntersections = intersectCircle.intersect(circle, intersectionPoint);
	}
	return numIntersections;
}


rg_INT Sphere::intersect(const Arc3D& arc, rg_Point3D* intersectionPoint) const
{
    Plane planeContainingArc = arc.getPlaneContainingThisCircle();

    Circle3D intersectingCircle;
    rg_BOOL isIntersect = intersect(planeContainingArc, intersectingCircle);

    rg_INT numIntersection = 0;
    if ( isIntersect ) {
        numIntersection = arc.intersect(intersectingCircle, intersectionPoint);
    }

    return numIntersection;
}



rg_BOOL Sphere::intersect(const Plane& plane, Circle3D& intersectingCircle) const
{
    rg_REAL dist = plane.distanceFromPoint(m_center);

    if ( rg_ABS(dist) > m_radius ) {
        return rg_FALSE;
    }

    rg_Point3D circleCenter = plane.projectPointOnPlane(m_center);

    rg_REAL    squaredRadius = (m_radius*m_radius) - (dist*dist);
    rg_REAL    circleRadius  = sqrt( squaredRadius );

    intersectingCircle.setCircle(circleCenter, circleRadius, plane.getNormal() );

    return rg_TRUE;
}



SphereIntersectionType Sphere::intersect(const Sphere& sphere, Circle3D& intersectingCircle) const
{
    rg_REAL distBtwCenters = m_center.distance(sphere.m_center);
    rg_REAL sumOfTwoRadii  = m_radius + sphere.m_radius;

    SphereIntersectionType intersectionType = SI_NOT_INTERSECT;
    if ( distBtwCenters > sumOfTwoRadii ) {
        intersectionType = SI_NOT_INTERSECT;
    }
    else if ( m_radius >= (distBtwCenters+sphere.m_radius) ) {
        intersectionType = SI_CONTAIN;
    }
    else if ( (distBtwCenters+m_radius) <= sphere.m_radius ) {
        intersectionType = SI_CONTAINED;
    }
    else {
        if ( m_radius == sphere.m_radius ) {
            //  compute a normal of the plane containing an intersecting circle.
            rg_Point3D circleNormal     = sphere.m_center - m_center;

            rg_Point3D circleCenter     = (m_center + sphere.m_center)*0.5;

            rg_REAL    halfMagitudeOfNormal = circleNormal.magnitude()*0.5;
            rg_REAL    squaredRadius        = (m_radius*m_radius) - (halfMagitudeOfNormal*halfMagitudeOfNormal);
            rg_REAL    circleRadius         = sqrt( squaredRadius );

            circleNormal.normalize();

            intersectingCircle.setCircle(circleCenter, circleRadius, circleNormal);
        }
        else {
            //  compute a normal of the plane containing an intersecting circle.
            rg_Point3D circleNormal = sphere.m_center - m_center;
            circleNormal.normalize();

            rg_REAL a = m_radius;
            rg_REAL b = sphere.m_radius;
            rg_REAL c = distBtwCenters;

            rg_REAL cosBeta = ( a*a + c*c - b*b )/(2.0*a*c);
            rg_REAL distFromCenterToCircleCenter = a*cosBeta;

            rg_Point3D circleCenter = m_center + (distFromCenterToCircleCenter*circleNormal);
            rg_REAL    circleRadius = sqrt( a*a - distFromCenterToCircleCenter*distFromCenterToCircleCenter);

            intersectingCircle.setCircle(circleCenter, circleRadius, circleNormal);

//             //  compute the radius of an intersecting circle.
// 	        rg_REAL lenSide[3]    = { distBtwCenters, m_radius, sphere.m_radius };	
// 	        rg_REAL halfPerimeter = (lenSide[0] + lenSide[1] + lenSide[2])*0.5;
// 	        rg_REAL squaredArea   = halfPerimeter*(halfPerimeter-lenSide[0])*(halfPerimeter-lenSide[1])*(halfPerimeter-lenSide[2]);
// 	        rg_REAL triangleArea  = sqrt( squaredArea );
// 
//             rg_REAL circleRadius  = triangleArea*2.0/distBtwCenters;
// 
//             //  compute the center of an intersecting circle.
//             rg_REAL squaredDistFromCircleCenter = (m_radius*m_radius) - (circleRadius*circleRadius);
//             rg_REAL distFromCircleCenter        = sqrt( squaredDistFromCircleCenter );
// 
//             rg_Point3D circleCenter;
//             if ( m_radius >= sphere.m_radius ) {
//                 circleCenter = m_center + (circleNormal*distFromCircleCenter);
//             }
//             else {
//                 rg_REAL distWithGreatCircleAsIntersection = sqrt( (sphere.m_radius*sphere.m_radius) - (m_radius*m_radius) );
//                 if ( distBtwCenters >= distWithGreatCircleAsIntersection ) {
//                     circleCenter = m_center + (circleNormal*distFromCircleCenter);
//                 }
//                 else {
//                     circleCenter = m_center - (circleNormal*distFromCircleCenter);
//                 }     
//             }
// 
//             intersectingCircle.setCircle(circleCenter, circleRadius, circleNormal);
        }
        intersectionType = SI_INTERSECT;
    }

    return intersectionType;
}


SphereIntersectionType Sphere::intersect(const Sphere& sphere) const
{
	rg_REAL distBtwCenters = m_center.distance(sphere.m_center);
	rg_REAL sumOfTwoRadii  = m_radius + sphere.m_radius;

	SphereIntersectionType intersectionType = SI_NOT_INTERSECT;
	if ( distBtwCenters > sumOfTwoRadii ) 
	{
		intersectionType = SI_NOT_INTERSECT;
	}
	else if ( m_radius >= (distBtwCenters+sphere.m_radius) ) 
	{
		intersectionType = SI_CONTAIN;
	}
	else if ( (distBtwCenters+m_radius) <= sphere.m_radius ) 
	{
		intersectionType = SI_CONTAINED;
	}
	else 
	{
		intersectionType = SI_INTERSECT;
	}

	return intersectionType;
}

rg_INT computeIntersectionPointsAmong(const Sphere& sphere1, 
									  const Sphere& sphere2, 
									  const Sphere& sphere3,
									  rg_Point3D* intersectionPoint)
{
	Circle3D intersectCircle12;
	SphereIntersectionType type = sphere1.intersect(sphere2, intersectCircle12);

	rg_INT numIntersections = 0;
	if(type == SI_INTERSECT)
	{
		numIntersections = sphere3.intersect(intersectCircle12, intersectionPoint);	
	}
	
	return numIntersections;	
}

rg_BOOL isIntersectionOfFirstTwoIncludedInUnionOfLastTwo(const Sphere& sphere1, 
                                                         const Sphere& sphere2, 
                                                         const Sphere& sphere3, 
                                                         const Sphere& sphere4)
{
	rg_INT     intersectionPointInclusion[ 4 ][ 2 ] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};
	rg_Point3D intersectionPoint[ 4 ][ 2 ];
	rg_INT     numIntersectionPoint[ 4 ];
	doContainmentTestOfTripletSpheresIntersectionForFourthSphere(sphere1, 
		                                                         sphere2, 
																 sphere3, 
																 sphere4, 
																 intersectionPointInclusion,
																 intersectionPoint,
																 numIntersectionPoint);

	rg_INT intersectionCount[ 4 ] = {0, 0, 0, 0};
	rg_INT numPairIncluded = 0;
	rg_INT i;
	for(i = 0;i < 4;i++)
	{
		intersectionCount[i] = intersectionPointInclusion[i][0] + intersectionPointInclusion[i][1];
		if(intersectionCount[i] == 2)
			numPairIncluded++;
	}
	
	Sphere   sphere[ 4 ]     = {sphere1, sphere2, sphere3, sphere4};
	rg_BOOL bIncluded = rg_TRUE;

	if(numPairIncluded == 2)
	{
		rg_INDEX indexForTripletSpheres[ 4 ][ 3 ] = {{0,1,2}, {0,1,3}, {1,2,3}, {0,2,3}};
		rg_INDEX indexForFourthSphere[ 4 ] = {3, 2, 0, 1};
		rg_Point3D* intersections = new rg_Point3D[ 2 ];
		rg_INT numOfIntersections = 0;

		rg_BOOL bTwoPairsIncluded = rg_TRUE;		
		
		for(i = 0;i < 3;i++)
		{
			numOfIntersections =
				computeIntersectionPointsAmong(sphere[indexForTripletSpheres[i][0]],
											   sphere[indexForTripletSpheres[i][1]], 
											   sphere[indexForTripletSpheres[i][2]], intersections);

			if(numOfIntersections != 2 || 
			   sphere[indexForFourthSphere[i]].doesContain(intersections[ 0 ]) == rg_FALSE ||
			   sphere[indexForFourthSphere[i]].doesContain(intersections[ 1 ]) == rg_FALSE   )
			{
				bTwoPairsIncluded = rg_FALSE;
				break;
			}
		}
		if(bTwoPairsIncluded)
		{
			numOfIntersections =
				computeIntersectionPointsAmong(sphere[indexForTripletSpheres[3][0]], 
											   sphere[indexForTripletSpheres[3][1]], 
											   sphere[indexForTripletSpheres[3][2]], intersections);

			if(numOfIntersections != 2 || 
			   sphere[indexForFourthSphere[3]].doesContain(intersections[ 0 ]) == rg_TRUE ||
			   sphere[indexForFourthSphere[3]].doesContain(intersections[ 1 ]) == rg_TRUE   )
			{
				bIncluded = rg_FALSE;
			}
		}
		else
			bIncluded = rg_FALSE;

		delete [] intersections;
	}
	else if(numPairIncluded == 4)
	{
		rg_INDEX indexForTripletSpheres[ 4 ][ 3 ] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
		rg_INDEX indexForFourthSphere[ 4 ]        = {3, 2, 1, 0};	

		rg_Point3D intersection[2] = {intersectionPoint[0][0], intersectionPoint[0][1]};
		for(i = 0;i < 4;i++)
		{
			if(!sphere[indexForFourthSphere[i]].isOnSphere(intersectionPoint[i][0]) ||
			   !sphere[indexForFourthSphere[i]].isOnSphere(intersectionPoint[i][1]) ||
			   intersection[0] != intersectionPoint[i][0]                           ||
			   intersection[1] != intersectionPoint[i][1]                           )
			{
				bIncluded = rg_FALSE;
				break;
			}
		}
		
		Plane plane;
		plane.definePlaneByThreePoints(sphere[0].getCenter(), sphere[1].getCenter(), sphere[2].getCenter());
		if(plane.isOnThisPlane(sphere[3].getCenter()))
		{
			rg_REAL areaOfTriangle[ 4 ] = {0., 0., 0., 0.};
			for(i = 0;i < 4;i++)
			{
				areaOfTriangle[ i ] = 
					rg_GeoFunc::computeTriangleArea(sphere[indexForTripletSpheres[i][0]].getCenter(), 
										 sphere[indexForTripletSpheres[i][1]].getCenter(), 
										 sphere[indexForTripletSpheres[i][2]].getCenter());
			}

			//rg_INT indexForDiagonal[ 4 ] = {-1, -1, -1, -1};
			if(rg_EQ(areaOfTriangle[0]+areaOfTriangle[1], areaOfTriangle[2]+areaOfTriangle[3]))
			{
	// 			indexForDiagonal[0] = 0;
	// 			indexForDiagonal[1] = 1;
	// 			indexForDiagonal[2] = 2;
	// 			indexForDiagonal[3] = 3;
				bIncluded = rg_TRUE;
			}
			else
			{
				bIncluded = rg_FALSE;
			}
		}
		else
			bIncluded = rg_FALSE;

	}
	else
		bIncluded = rg_FALSE;

	return bIncluded;
}

rg_INT computeSphereWithGivenRadiusAndTangentTo3Spheres(
                                                 const rg_REAL& radius,
                                                       Sphere*  sphere, 
                                                       Sphere*  tangentSphere )
{
    //  
    //  (x - s1.center.x)^2 + (y - s1.center.y)^2 + (z - s1.center.z)^2 = s1.radius + radius 
    //  (x - s2.center.x)^2 + (y - s2.center.y)^2 + (z - s2.center.z)^2 = s2.radius + radius 
    //  (x - s3.center.x)^2 + (y - s3.center.y)^2 + (z - s3.center.z)^2 = s3.radius + radius 
    //
    //  To make the system of equations simple, 
    //  translate s3.center into the origin.
    //
    //  system of equations after translating the center of s3 to the origin
    //  eq1: (x-x1)^2 + (y-y1)^2 + (z-z1)^2 = (r1)^2
    //  eq2: (x-x2)^2 + (y-y2)^2 + (z-z2)^2 = (r2)^2
    //  eq3: x^2      + y^2      + z^2      = (r3)^2

    rg_REAL x1 = sphere[0].m_center.getX() - sphere[2].m_center.getX();
	rg_REAL y1 = sphere[0].m_center.getY() - sphere[2].m_center.getY();
	rg_REAL z1 = sphere[0].m_center.getZ() - sphere[2].m_center.getZ();
	rg_REAL r1 = sphere[0].m_radius          + radius;

	rg_REAL x2 = sphere[1].m_center.getX() - sphere[2].m_center.getX();
	rg_REAL y2 = sphere[1].m_center.getY() - sphere[2].m_center.getY();
	rg_REAL z2 = sphere[1].m_center.getZ() - sphere[2].m_center.getZ();
	rg_REAL r2 = sphere[1].m_radius          + radius;

	rg_REAL r3 = sphere[2].m_radius          + radius;
    
    //system of equations
	//eq1-eq3: x1*x + y1*y + z1*z = rhs1((x1^2 + y1^2 + z1^2 + r3^2 - r1^2)/2)
	//eq2-eq3: x2*x + y2*y + z2*z = rhs2((x2^2 + y2^2 + z2^2 + r3^2 - r2^2)/2)

	rg_REAL rhs1 = (x1*x1 + y1*y1 + z1*z1 + r3*r3 - r1*r1)*0.5;
	rg_REAL rhs2 = (x2*x2 + y2*y2 + z2*z2 + r3*r3 - r2*r2)*0.5;

    rg_REAL determinant1 = x1*y2 - x2*y1;
	if( !rg_ZERO(determinant1) )  {
        //system of equations
	    //eq1: x1*x + y1*y = rhs1 - z1*z 
	    //eq2: x2*x + y2*y = rhs2 - z2*z 
        //

        if ( rg_ZERO(z1) && rg_ZERO(z2) ) {
		    //  x = (y2*rhs1 - y1*rhs2)/determinant1;
		    //  y = (x1*rhs2 - x2*rhs1)/determinant1;
            //
            //  eq3 <- x, y            
	        //  z^2      = (r3)^2 - x^2 - y^2       
		    rg_REAL x = (y2*rhs1 - y1*rhs2)/determinant1;
		    rg_REAL y = (x1*rhs2 - x2*rhs1)/determinant1;
            rg_REAL zz = (r3*r3) - (x*x) - (y*y);

            if ( rg_POS( zz ) ) {
                rg_REAL z = sqrt( zz );

                x += sphere[2].m_center.getX();
                y += sphere[2].m_center.getY();
                tangentSphere[0].m_center.setPoint(x, y, -z+sphere[2].m_center.getZ());
                tangentSphere[0].m_radius = radius;
			    tangentSphere[1].m_center.setPoint(x, y, z + sphere[2].m_center.getZ());
                tangentSphere[1].m_radius = radius;
                    
                return 2;
            }
        }
        else {
            //  x = a_x + b_x*z
            //  y = a_y + b_y*z
            //
            //  eq3 <- x, y
            //  a_z*z^2 + b_z*z + c_z = 0

		    rg_REAL a_x = (y2*rhs1 - y1*rhs2)/determinant1;
            rg_REAL b_x = (y1*z2   - y2*z1)  /determinant1;
		    rg_REAL a_y = (x1*rhs2 - x2*rhs1)/determinant1;
            rg_REAL b_y = (x2*z1   - x1*z2)  /determinant1;

		    rg_REAL a_z = b_x*b_x + b_y*b_y + 1.0;
            rg_REAL b_z = 2.0*(a_x*b_x + a_y*b_y);
            rg_REAL c_z = a_x*a_x + a_y*a_y - (r3*r3);

		    rg_REAL det = b_z*b_z - 4.*a_z*c_z; //b^2 - 4ac
		    if(det >= 0) {
                rg_REAL sqrtDet = sqrt(det); 
                if ( rg_ZERO( sqrtDet ) )  {
                    rg_REAL z = -b_z/(2.* a_z);
            
			        rg_REAL x = a_x + b_x*z + sphere[2].m_center.getX();
			        rg_REAL y = a_y + b_y*z + sphere[2].m_center.getY();
			        z += sphere[2].m_center.getZ();

			        tangentSphere[0].m_center.setPoint(x, y, z);
                    tangentSphere[0].m_radius = radius;

                    return 1;
                }
                else {
		            rg_REAL z1 = (-b_z + sqrtDet)/(2.* a_z);
		            rg_REAL z2 = (-b_z - sqrtDet)/(2.* a_z);

			        rg_REAL x = a_x + b_x*z1 + sphere[2].m_center.getX();
			        rg_REAL y = a_y + b_y*z1 + sphere[2].m_center.getY();
			        z1 += sphere[2].m_center.getZ();

			        tangentSphere[0].m_center.setPoint(x, y, z1);
                    tangentSphere[0].m_radius = radius;


			        x = a_x + b_x*z2 + sphere[2].m_center.getX();
			        y = a_y + b_y*z2 + sphere[2].m_center.getY();
			        z2 += sphere[2].m_center.getZ();

			        tangentSphere[1].m_center.setPoint(x, y, z2);
                    tangentSphere[1].m_radius = radius;

                    return 2;
                }
            }
        }
    }


	rg_REAL determinant2 = x1*z2 - x2*z1;
    if( !rg_ZERO(determinant2) )  {
        //system of equations
	    //eq1: x1*x + z1*z = rhs1 - y1*y 
	    //eq2: x2*x + z2*z = rhs2 - y2*y 
        //

        if ( rg_ZERO(y1) && rg_ZERO(y2) ) {
		    //  x = (z2*rhs1 - z1*rhs2)/determinant2;
		    //  z = (x1*rhs2 - x2*rhs1)/determinant2;
            //
            //  eq3 <- x, z            
	        //  y^2      = (r3)^2 - x^2 - z^2       
		    rg_REAL x = (z2*rhs1 - z1*rhs2)/determinant2;
		    rg_REAL z = (x1*rhs2 - x2*rhs1)/determinant2;
            rg_REAL yy = (r3*r3) - (x*x) - (z*z);

            if ( rg_POS( yy ) ) {
                rg_REAL y = sqrt( yy );

                x += sphere[2].m_center.getX();
                z += sphere[2].m_center.getZ();
			    tangentSphere[0].m_center.setPoint(x, -y+sphere[2].m_center.getY(), z);
                tangentSphere[0].m_radius = radius;
			    tangentSphere[1].m_center.setPoint(x, y+sphere[2].m_center.getY(), z);
                tangentSphere[1].m_radius = radius;
                    
                return 2;
            }
        }
        else {
            //  x = a_x + b_x*y
            //  z = a_z + b_z*y
            //
            //  eq3 <- x, z
            //  a_y*y^2 + b_y*y + c_y = 0
		    rg_REAL a_x = (z2*rhs1 - z1*rhs2)/determinant2;
            rg_REAL b_x = (z1*y2   - z2*y1)  /determinant2;
		    rg_REAL a_z = (x1*rhs2 - x2*rhs1)/determinant2;
            rg_REAL b_z = (x2*y1   - x1*y2)  /determinant2;

		    rg_REAL a_y = b_x*b_x + b_z*b_z + 1.0;
            rg_REAL b_y = 2.0*(a_x*b_x + a_z*b_z);
            rg_REAL c_y = a_x*a_x + a_z*a_z - (r3*r3);

		    rg_REAL det = b_y*b_y - 4.*a_y*c_y; //b^2 - 4ac
		    if(det >= 0) {
                rg_REAL sqrtDet = sqrt(det); 
                if ( rg_ZERO( sqrtDet ) )  {
                    rg_REAL y = -b_z/(2.* a_z);
            
			        rg_REAL x = a_x + b_x*y + sphere[2].m_center.getX();
			        rg_REAL z = a_z + b_z*y + sphere[2].m_center.getZ();
			        y += sphere[2].m_center.getY();

			        tangentSphere[0].m_center.setPoint(x, y, z);
                    tangentSphere[0].m_radius = radius;

                    return 1;
                }
                else {
		            rg_REAL y1 = (-b_y + sqrtDet)/(2.* a_y);
		            rg_REAL y2 = (-b_y - sqrtDet)/(2.* a_y);

			        rg_REAL x = a_x + b_x*y1 + sphere[2].m_center.getX();
			        rg_REAL z = a_z + b_z*y1 + sphere[2].m_center.getZ();
			        y1 += sphere[2].m_center.getY();

			        tangentSphere[0].m_center.setPoint(x, y1, z);
                    tangentSphere[0].m_radius = radius;


			        x = a_x + b_x*y2 + sphere[2].m_center.getX();
			        z = a_z + b_z*y2 + sphere[2].m_center.getZ();
			        y2 += sphere[2].m_center.getY();

			        tangentSphere[1].m_center.setPoint(x, y2, z);
                    tangentSphere[1].m_radius = radius;

                    return 2;
                }
            }
        }
    }

    
    rg_REAL determinant3 = y1*z2 - y2*z1;
    if( !rg_ZERO(determinant3) )  {
        //system of equations
	    //eq1: y1*y + z1*r = rhs1 - x1*x 
	    //eq2: y2*y + z2*r = rhs2 - x2*x 
        //
        if ( rg_ZERO(x1) && rg_ZERO(x2) ) {
		    //  y = (z2*rhs1 - z1*rhs2)/determinant3;
		    //  z = (y1*rhs2 - y2*rhs1)/determinant3;
            //
            //  eq3 <- y, z            
	        //  x^2      = (r3)^2 - y^2 - x^2       
		    rg_REAL y = (z2*rhs1 - z1*rhs2)/determinant3;
		    rg_REAL z = (y1*rhs2 - y2*rhs1)/determinant3;
            rg_REAL xx = (r3*r3) - (y*y) - (z*z);

            if ( rg_POS( xx ) ) {
                rg_REAL x = sqrt(xx);

                y += sphere[2].m_center.getY();
                z += sphere[2].m_center.getZ();
                tangentSphere[0].m_center.setPoint(-x+sphere[2].m_center.getX(), y, z);
                tangentSphere[0].m_radius = radius;
			    tangentSphere[1].m_center.setPoint(x+sphere[2].m_center.getX(), y, z);
                tangentSphere[1].m_radius = radius;
                    
                return 2;
            }
        }
        else {
            //  y = a_y + b_y*x
            //  z = a_z + b_z*x
            //
            //  eq3 <- y, z
            //  a_x*x^2 + b_x*x + c_x = 0

		    rg_REAL a_y = (z2*rhs1 - z1*rhs2)/determinant3;
            rg_REAL b_y = (z1*x2   - z2*x1)  /determinant3;
		    rg_REAL a_z = (y1*rhs2 - y2*rhs1)/determinant3;
            rg_REAL b_z = (y2*x1   - y1*x2)  /determinant3;

		    rg_REAL a_x = b_y*b_y + b_z*b_z + 1.0;
            rg_REAL b_x = 2.0*(a_y*b_y + a_z*b_z);
            rg_REAL c_x = a_y*a_y + a_z*a_z - (r3*r3);

		    rg_REAL det = b_x*b_x - 4.*a_x*c_x; //b^2 - 4ac
		    if(det >= 0) {
                rg_REAL sqrtDet = sqrt(det); 
                if ( rg_ZERO( sqrtDet ) )  {
                    rg_REAL x = -b_x/(2.* a_x);
            
			        rg_REAL y = a_y + b_y*x + sphere[2].m_center.getY();
			        rg_REAL z = a_z + b_z*x + sphere[2].m_center.getZ();
			        x += sphere[2].m_center.getX();

			        tangentSphere[0].m_center.setPoint(x, y, z);
                    tangentSphere[0].m_radius = radius;

                    return 1;
                }
                else {
		            rg_REAL x1 = (-b_x + sqrtDet)/(2.* a_x);
		            rg_REAL x2 = (-b_x - sqrtDet)/(2.* a_x);

			        rg_REAL y = a_y + b_y*x1 + sphere[2].m_center.getY();
			        rg_REAL z = a_z + b_z*x1 + sphere[2].m_center.getZ();
			        x1 += sphere[2].m_center.getX();

			        tangentSphere[0].m_center.setPoint(x1, y, z);
                    tangentSphere[0].m_radius = radius;


			        y = a_y + b_y*x2 + sphere[2].m_center.getY();
			        z = a_z + b_z*x2 + sphere[2].m_center.getZ();
			        x2 += sphere[2].m_center.getX();

			        tangentSphere[1].m_center.setPoint(x2, y, z);
                    tangentSphere[1].m_radius = radius;

                    return 2;
                }
            }
        }
    }

    return 0;
}



rg_INT computeSphereTangentTo4SpheresOutside(const Sphere& s1, 
                                             const Sphere& s2, 
                                             const Sphere& s3, 
                                             const Sphere& s4,
                                             Sphere* tangentSpheres)
{	
    //donguk's version (2004. 7. 5)
    //Modified by Youngsong Cho (2005. 8. 12)
    Sphere balls[4] = {s1, s2, s3, s4};
    sort4SpheresByRadiusInDescendingPowers( balls );

	if(    (balls[0].m_radius == balls[1].m_radius) 
        && (balls[1].m_radius == balls[2].m_radius) 
        && (balls[2].m_radius == balls[3].m_radius) )
	{
        return computeSphereTangentTo4SpheresWithSameRadiusOutside( balls[0], balls[1], balls[2], balls[3], tangentSpheres);
	}
	else
	{
        /*
        Sphere balls[4];
        if ( s1.m_radius < s2.m_radius )  {
            if ( s1.m_radius < s3.m_radius )  {
                if ( s1.m_radius < s4.m_radius )  {
                    balls[0] = s2;
                    balls[1] = s3;
                    balls[2] = s4;
                    balls[3] = s1;
                }
                else  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s3;
                    balls[3] = s4;
                }
            }
            else  {
                if ( s3.m_radius < s4.m_radius )  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s4;
                    balls[3] = s3;
                }
                else  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s3;
                    balls[3] = s4;
                }
            }
        }
        else  {
            if ( s2.m_radius < s3.m_radius )  {
                if ( s2.m_radius < s4.m_radius )  {
                    balls[0] = s1;
                    balls[1] = s3;
                    balls[2] = s4;
                    balls[3] = s2;
                }
                else  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s3;
                    balls[3] = s4;
                }
            }
            else  {
                if ( s3.m_radius < s4.m_radius )  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s4;
                    balls[3] = s3;
                }
                else  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s3;
                    balls[3] = s4;
                }
            }
        }
        */

		rg_REAL x1 = balls[0].m_center.getX() - balls[3].m_center.getX();
		rg_REAL y1 = balls[0].m_center.getY() - balls[3].m_center.getY();
		rg_REAL z1 = balls[0].m_center.getZ() - balls[3].m_center.getZ();
		rg_REAL r1 = balls[0].m_radius        - balls[3].m_radius;

		rg_REAL x2 = balls[1].m_center.getX() - balls[3].m_center.getX();
		rg_REAL y2 = balls[1].m_center.getY() - balls[3].m_center.getY();
		rg_REAL z2 = balls[1].m_center.getZ() - balls[3].m_center.getZ();
		rg_REAL r2 = balls[1].m_radius        - balls[3].m_radius;

		rg_REAL x3 = balls[2].m_center.getX() - balls[3].m_center.getX();
		rg_REAL y3 = balls[2].m_center.getY() - balls[3].m_center.getY();
		rg_REAL z3 = balls[2].m_center.getZ() - balls[3].m_center.getZ();
		rg_REAL r3 = balls[2].m_radius        - balls[3].m_radius;


        //rg_REAL determinant1 = getDeterminantOf33Matrix(x1,y1,z1,x2,y2,z2,x3,y3,z3);
		//rg_REAL determinant2 = getDeterminantOf33Matrix(x1,y1,r1,x2,y2,r2,x3,y3,r3);
		//rg_REAL determinant3 = getDeterminantOf33Matrix(x1,r1,z1,x2,r2,z2,x3,r3,z3);
		//rg_REAL determinant4 = getDeterminantOf33Matrix(r1,y1,z1,r2,y2,z2,r3,y3,z3);

        //rg_REAL determinant1 = x1*(y2*z3 - z2*y3) - y1*(x2*z3 - z2*x3) + z1*(x2*y3 - y2*x3);
        //rg_REAL determinant2 = x1*(y2*r3 - r2*y3) - y1*(x2*r3 - r2*x3) + r1*(x2*y3 - y2*x3);
        //rg_REAL determinant3 = x1*(r2*z3 - z2*r3) - r1*(x2*z3 - z2*x3) + z1*(x2*r3 - r2*x3);
        //rg_REAL determinant4 = r1*(y2*z3 - z2*y3) - y1*(r2*z3 - z2*r3) + z1*(r2*y3 - y2*r3);


        //system of equations
		//eq1: x1*x + y1*y + z1*z = rhs1
		//eq2: x2*x + y2*y + z2*z = rhs2
		//eq3: x3*x + y3*y + z3*z = rhs3
		//eq4: x*x + y*y + z*z = r*r

		rg_REAL rhs1 = (x1*x1+y1*y1+z1*z1-r1*r1)/2.0;
		rg_REAL rhs2 = (x2*x2+y2*y2+z2*z2-r2*r2)/2.0;
		rg_REAL rhs3 = (x3*x3+y3*y3+z3*z3-r3*r3)/2.0;


		//finding non-null minor
        rg_REAL determinant1 = x1*(y2*z3 - z2*y3) - y1*(x2*z3 - z2*x3) + z1*(x2*y3 - y2*x3);
		if( !rg_ZERO(determinant1) ) 
		{
			//a_r * r^2 + b_r * r + c_r = 0

			//x = a_x * r + b_x;
			//y = a_y * r + b_y;
			//z = a_z * r + b_z;

			rg_REAL a_x = ((z2*y3-y2*z3)*r1  +(z3*y1-z1*y3)*r2  +(z1*y2-z2*y1)*r3)  /determinant1;
			rg_REAL b_x = ((y2*z3-z2*y3)*rhs1+(z1*y3-z3*y1)*rhs2+(z2*y1-z1*y2)*rhs3)/determinant1;
			rg_REAL a_y = ((z3*x2-z2*x3)*r1  +(x3*z1-x1*z3)*r2  +(z2*x1-z1*x2)*r3)  /determinant1;
			rg_REAL b_y = ((z2*x3-z3*x2)*rhs1+(x1*z3-x3*z1)*rhs2+(z1*x2-z2*x1)*rhs3)/determinant1;
			rg_REAL a_z = ((y2*x3-x2*y3)*r1  +(x1*y3-x3*y1)*r2  +(x2*y1-y2*x1)*r3)  /determinant1;
			rg_REAL b_z = ((x2*y3-y2*x3)*rhs1+(x3*y1-x1*y3)*rhs2+(y2*x1-x2*y1)*rhs3)/determinant1;


			// r = (-b_r + sqrt(b_r^2 - 4. * a_r * c_r))/(2.* a_r);
			// r = (-b_r - sqrt(b_r^2 - 4. * a_r * c_r))/(2.* a_r);
			rg_REAL a_r = a_y*a_y+a_x*a_x+a_z*a_z-1.0;
			rg_REAL b_r = 2.0*b_x*a_x+2.0*b_y*a_y+2.0*b_z*a_z;
			rg_REAL c_r = b_x*b_x+b_z*b_z+b_y*b_y;

			rg_REAL det = b_r*b_r - 4.*a_r*c_r; //b^2 - 4ac
			if(det < 0) 
				return 0;

            rg_REAL sqrtDet = sqrt(det);
			rg_REAL radius1 = (-b_r + sqrtDet)/(2.* a_r);
			rg_REAL radius2 = (-b_r - sqrtDet)/(2.* a_r);


            if( radius1 > 0 && radius2 > 0 )  {

				rg_REAL x = a_x * radius1 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y * radius1 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z * radius1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z);
                tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

				x = a_x * radius2 + b_x + balls[3].m_center.getX();
				y = a_y * radius2 + b_y + balls[3].m_center.getY();
				z = a_z * radius2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z);
                tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;

				return 2;
			}
			else if( radius1 > 0 )  {

				rg_REAL x = a_x * radius1 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y * radius1 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z * radius1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z);
                tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

				return 1;
			}
			else if( radius2 > 0 )  {
				rg_REAL x = a_x * radius2 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y * radius2 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z * radius2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z);
                tangentSpheres[0].m_radius = radius2 - balls[3].m_radius;

				return 1;
			}
			else  {
				return 0;
			}

		}

        rg_REAL determinant2 = x1*(y2*r3 - r2*y3) - y1*(x2*r3 - r2*x3) + r1*(x2*y3 - y2*x3);
		if( !rg_ZERO(determinant2) ) 
		{
			//a_z * z^2 + b_z * z + c_z = 0

			//x = a_x * z + b_x;
			//y = a_y * z + b_y;
			//r = a_r * z + b_r;

			rg_REAL a_x = ((r2*y3-y2*r3)*z1  +(r3*y1-r1*y3)*z2  +(r1*y2-r2*y1)*z3)  /determinant2;
			rg_REAL b_x = ((y2*r3-r2*y3)*rhs1+(r1*y3-y1*r3)*rhs2+(r2*y1-y2*r1)*rhs3)/determinant2;
			rg_REAL a_y = ((r3*x2-r2*x3)*z1  +(x3*r1-x1*r3)*z2  +(r2*x1-r1*x2)*z3)  /determinant2;
			rg_REAL b_y = ((r2*x3-r3*x2)*rhs1+(x1*r3-x3*r1)*rhs2+(r1*x2-r2*x1)*rhs3)/determinant2;
			rg_REAL a_r = ((y2*x3-x2*y3)*z1  +(x1*y3-x3*y1)*z2  +(x2*y1-y2*x1)*z3)  /determinant2;
			rg_REAL b_r = ((x2*y3-y2*x3)*rhs1+(x3*y1-x1*y3)*rhs2+(y2*x1-x2*y1)*rhs3)/determinant2;

			rg_REAL a_z = a_x*a_x+a_y*a_y+1.0-a_r*a_r;
			rg_REAL b_z = 2.0*b_x*a_x+2.0*b_y*a_y-2.0*b_r*a_r;
			rg_REAL c_z = b_x*b_x+b_y*b_y-b_r*b_r;

			rg_REAL det = b_z*b_z - 4.0*a_z*c_z;
			if(det < 0)  
				return 0;


            rg_REAL sqrtDet = sqrt(det);
			rg_REAL zValue1 = (-b_z + sqrtDet)/(2.0*a_z);
			rg_REAL zValue2 = (-b_z - sqrtDet)/(2.0*a_z);

			rg_REAL rad1 = a_r*zValue1 + b_r;
			rg_REAL rad2 = a_r*zValue2 + b_r;

			if( rad1 > 0 && rad2 > 0)  {

				rg_REAL x = a_x*zValue1 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y*zValue1 + b_y + balls[3].m_center.getY();
                rg_REAL z = zValue1           + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				x = a_x*zValue2 + b_x + balls[3].m_center.getX();
				y = a_y*zValue2 + b_y + balls[3].m_center.getY();
                z = zValue2           + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z );
                tangentSpheres[1].m_radius = rad2 - balls[3].m_radius;

				return 2;
			}
			else if( rad1 > 0 )  {

				rg_REAL x = a_x*zValue1 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y*zValue1 + b_y + balls[3].m_center.getY();
                rg_REAL z = zValue1           + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				return 1;
			}
			else if( rad2 > 0 )  {

				rg_REAL x = a_x*zValue2 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y*zValue2 + b_y + balls[3].m_center.getY();
                rg_REAL z = zValue2           + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad2 - balls[3].m_radius;

				return 1;
			}
			else  {
				return 0;
			}
		}


        rg_REAL determinant3 = x1*(r2*z3 - z2*r3) - r1*(x2*z3 - z2*x3) + z1*(x2*r3 - r2*x3);
		if( !rg_ZERO(determinant3) ) 
		{
			//a_y * y^2 + b_y * y + c_y = 0

			//x = a_x * y + b_x;
			//r = a_r * y + b_r;
			//z = a_z * y + b_z;
			rg_REAL a_x = ((z2*r3-r2*z3)*y1  +(z3*r1-z1*r3)*y2  +(z1*r2-z2*r1)*y3)  /determinant3;
			rg_REAL b_x = ((r2*z3-z2*r3)*rhs1+(z1*r3-z3*r1)*rhs2+(z2*r1-z1*r2)*rhs3)/determinant3;
			rg_REAL a_r = ((z3*x2-z2*x3)*y1  +(x3*z1-x1*z3)*y2  +(z2*x1-z1*x2)*y3)  /determinant3;
			rg_REAL b_r = ((z2*x3-z3*x2)*rhs1+(x1*z3-x3*z1)*rhs2+(z1*x2-z2*x1)*rhs3)/determinant3;
			rg_REAL a_z = ((r2*x3-x2*r3)*y1  +(x1*r3-x3*r1)*y2  +(x2*r1-r2*x1)*y3)  /determinant3;
			rg_REAL b_z = ((x2*r3-r2*x3)*rhs1+(x3*r1-x1*r3)*rhs2+(r2*x1-x2*r1)*rhs3)/determinant3;

			rg_REAL a_y = a_x*a_x+a_z*a_z+1.0-a_r*a_r;
			rg_REAL b_y = 2.0*b_x*a_x+2.0*b_z*a_z-2.0*b_r*a_r;
			rg_REAL c_y = b_x*b_x-b_r*b_r+b_z*b_z;

			rg_REAL det = b_y*b_y - 4.0*a_y*c_y;
			if(det < 0)
				return 0;


            rg_REAL sqrtDet = sqrt(det);
			rg_REAL yValue1 = (-b_y + sqrtDet)/(2.0*a_y);
			rg_REAL yValue2 = (-b_y - sqrtDet)/(2.0*a_y);

			rg_REAL rad1 = a_r*yValue1 + b_r;
			rg_REAL rad2 = a_r*yValue2 + b_r;
			if( rad1 > 0 && rad2 > 0)  {

                rg_REAL x = a_x*yValue1 + b_x + balls[3].m_center.getX();
                rg_REAL y = yValue1           + balls[3].m_center.getY();
				rg_REAL z = a_z*yValue1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;


				x = a_x*yValue2 + b_x + balls[3].m_center.getX();
                y = yValue2           + balls[3].m_center.getY();
				z = a_z*yValue2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z );
                tangentSpheres[1].m_radius = rad2 - balls[3].m_radius;

				return 2;
			}
			else if( rad1 > 0 )  {

                rg_REAL x = a_x*yValue1 + b_x + balls[3].m_center.getX();
                rg_REAL y = yValue1           + balls[3].m_center.getY();
				rg_REAL z = a_z*yValue1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				return 1;
			}
			else if( rad2 > 0 )  {

				rg_REAL x = a_x*yValue2 + b_x + balls[3].m_center.getX();
                rg_REAL y = yValue2           + balls[3].m_center.getY();
				rg_REAL z = a_z*yValue2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z );
                tangentSpheres[1].m_radius = rad2 - balls[3].m_radius;

				return 1;
			}
            else  {
				return 0;
			}
		
		}


        rg_REAL determinant4 = r1*(y2*z3 - z2*y3) - y1*(r2*z3 - z2*r3) + z1*(r2*y3 - y2*r3);
		if( !rg_ZERO(determinant4) ) 
		{
			//a_x * x^2 + b_x * x + c_x = 0

			//r = a_r * x + b_r;
			//y = a_y * x + b_y;
			//z = a_z * x + b_z;
			rg_REAL a_r = ((z2*y3-y2*z3)*x1  +(z3*y1-z1*y3)*x2  +(z1*y2-z2*y1)*x3)  /determinant4;
			rg_REAL b_r = ((y2*z3-z2*y3)*rhs1+(z1*y3-z3*y1)*rhs2+(z2*y1-z1*y2)*rhs3)/determinant4;
			rg_REAL a_y = ((z3*r2-z2*r3)*x1  +(r3*z1-r1*z3)*x2  +(z2*r1-z1*r2)*x3)  /determinant4;
			rg_REAL b_y = ((z2*r3-z3*r2)*rhs1+(r1*z3-r3*z1)*rhs2+(z1*r2-z2*r1)*rhs3)/determinant4;
			rg_REAL a_z = ((y2*r3-r2*y3)*x1  +(r1*y3-r3*y1)*x2  +(r2*y1-y2*r1)*x3)  /determinant4;
			rg_REAL b_z = ((r2*y3-y2*r3)*rhs1+(r3*y1-r1*y3)*rhs2+(y2*r1-r2*y1)*rhs3)/determinant4;

			rg_REAL a_x = a_z*a_z+a_y*a_y+1.0-a_r*a_r;
			rg_REAL b_x = 2.0*b_y*a_y-2.0*b_r*a_r+2.0*b_z*a_z;
			rg_REAL c_x = b_z*b_z+b_y*b_y-b_r*b_r;

			rg_REAL det = b_x*b_x - 4.0*a_x*c_x;
			if(det < 0)
			{
				tangentSpheres = rg_NULL;
				return 0;
			}

            rg_REAL sqrtDet = sqrt(det);
			rg_REAL xValue1 = (-b_x + sqrtDet)/(2.0*a_x);
			rg_REAL xValue2 = (-b_x - sqrtDet)/(2.0*a_x);

			rg_REAL rad1 = a_r*xValue1 + b_r;
			rg_REAL rad2 = a_r*xValue2 + b_r;
			if( rad1 > 0 && rad2 > 0)  {

				rg_REAL x = xValue1           + balls[3].m_center.getX();
                rg_REAL y = a_y*xValue1 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z*xValue1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				x = xValue2           + balls[3].m_center.getX();
                y = a_y*xValue2 + b_y + balls[3].m_center.getY();
				z = a_z*xValue2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z );
                tangentSpheres[1].m_radius = rad2 - balls[3].m_radius;

                return 2;
			}
			else if( rad1 > 0 )  {

				rg_REAL x = xValue1           + balls[3].m_center.getX();
                rg_REAL y = a_y*xValue1 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z*xValue1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				return 1;
			}
			else if( rad2 > 0 )  {

				rg_REAL x = xValue2           + balls[3].m_center.getX();
                rg_REAL y = a_y*xValue2 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z*xValue2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad2 - balls[3].m_radius;

                return 1;
			}
			else  {
				return 0;
			}
		}
		
        
        return 0;
	}
	
    return 0;
}




rg_INT compute_tangent_sphere_of_three_ball_generators_outside_and_spherical_container_inside(const Sphere & b1, const Sphere & b2, const Sphere & b3, const Sphere & container, Sphere * tangentSpheres)
{
	// The spherical container corresponds to the index "2 (third element)".
	Sphere balls[4] = { b1, b2, container, b3 };

	rg_REAL x1 = balls[0].m_center.getX() - balls[3].m_center.getX();
	rg_REAL y1 = balls[0].m_center.getY() - balls[3].m_center.getY();
	rg_REAL z1 = balls[0].m_center.getZ() - balls[3].m_center.getZ();
	rg_REAL r1 = balls[0].m_radius - balls[3].m_radius;

	rg_REAL x2 = balls[1].m_center.getX() - balls[3].m_center.getX();
	rg_REAL y2 = balls[1].m_center.getY() - balls[3].m_center.getY();
	rg_REAL z2 = balls[1].m_center.getZ() - balls[3].m_center.getZ();
	rg_REAL r2 = balls[1].m_radius - balls[3].m_radius;

	rg_REAL x3 = balls[2].m_center.getX() - balls[3].m_center.getX();
	rg_REAL y3 = balls[2].m_center.getY() - balls[3].m_center.getY();
	rg_REAL z3 = balls[2].m_center.getZ() - balls[3].m_center.getZ();
	rg_REAL r3 = balls[2].m_radius + balls[3].m_radius;

	//system of equations
	//eq1: x1*x + y1*y + z1*z + r1*r = rhs1
	//eq2: x2*x + y2*y + z2*z + r2*r = rhs2
	//eq3: x3*x + y3*y + z3*z - r3*r = rhs3
	//eq4: x*x + y*y + z*z = r*r

	rg_REAL rhs1 = (x1*x1 + y1 * y1 + z1 * z1 - r1 * r1) / 2.0;
	rg_REAL rhs2 = (x2*x2 + y2 * y2 + z2 * z2 - r2 * r2) / 2.0;
	rg_REAL rhs3 = (x3*x3 + y3 * y3 + z3 * z3 - r3 * r3) / 2.0;

	//finding non-null minor
	rg_REAL determinant1 = (y1*z2 - y2 * z1)*x3 + (-y1 * z3 + y3 * z1)*x2 + (y2*z3 - y3 * z2)*x1;
	if (!rg_ZERO(determinant1))
	{
		//a_r * r^2 + b_r * r + c_r = 0

		//x = a_x * r + b_x;
		//y = a_y * r + b_y;
		//z = a_z * r + b_z;		

		rg_REAL a_x = ((y1*z2 - y2 * z1)*r3 + (y1*z3 - y3 * z1)*r2 + (-y2 * z3 + y3 * z2)*r1) / determinant1;
		rg_REAL b_x = ((y1*z2 - y2 * z1)*rhs3 + (-y1 * z3 + y3 * z1)*rhs2 + (y2*z3 - y3 * z2)*rhs1) / determinant1;
		rg_REAL a_y = ((-x1 * z2 + x2 * z1)*r3 + (-x1 * z3 + x3 * z1)*r2 + (x2*z3 - x3 * z2)*r1) / determinant1;
		rg_REAL b_y = ((-x1 * z2 + x2 * z1)*rhs3 + (x1*z3 - x3 * z1)*rhs2 + (-x2 * z3 + x3 * z2)*rhs1) / determinant1;
		rg_REAL a_z = ((y2*x1 - x2 * y1)*r3 + (y3*x1 - x3 * y1)*r2 + (-y3 * x2 + x3 * y2)*r1) / determinant1;
		rg_REAL b_z = ((y2*x1 - x2 * y1)*rhs3 + (-y3 * x1 + x3 * y1)*rhs2 + (y3*x2 - x3 * y2)*rhs1) / determinant1;

		//plug (x, y, z) in eq4 and
		//solve the quadratic equation in r
		// r = (-b_r + sqrt(b_r^2 - 4. * a_r * c_r))/(2.* a_r);
		// r = (-b_r - sqrt(b_r^2 - 4. * a_r * c_r))/(2.* a_r);

		rg_REAL a_r = a_x * a_x + a_y * a_y + a_z * a_z - 1.0;
		rg_REAL b_r = 2.0*a_x*b_x + 2.0*a_y*b_y + 2.0*a_z*b_z;
		rg_REAL c_r = b_x * b_x + b_y * b_y + b_z * b_z;

		rg_REAL det = b_r * b_r - 4.*a_r*c_r; //b^2 - 4ac
		if (det < 0)
			return 0;

		rg_REAL sqrtDet = sqrt(det);
		rg_REAL radius1 = (-b_r + sqrtDet) / (2.* a_r);
		rg_REAL radius2 = (-b_r - sqrtDet) / (2.* a_r);


		if (radius1 > 0 && radius2 > 0) 
		{
			// backtransform with radius1
			rg_REAL x = a_x * radius1 + b_x + balls[3].m_center.getX();
			rg_REAL y = a_y * radius1 + b_y + balls[3].m_center.getY();
			rg_REAL z = a_z * radius1 + b_z + balls[3].m_center.getZ();
			tangentSpheres[0].m_center.setPoint(x, y, z);
			tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

			// backtransform with radius2
			x = a_x * radius2 + b_x + balls[3].m_center.getX();
			y = a_y * radius2 + b_y + balls[3].m_center.getY();
			z = a_z * radius2 + b_z + balls[3].m_center.getZ();
			tangentSpheres[1].m_center.setPoint(x, y, z);
			tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;

			return 2;
		}
		else if (radius1 > 0) 
		{
			// backtransform with radius1
			rg_REAL x = a_x * radius1 + b_x + balls[3].m_center.getX();
			rg_REAL y = a_y * radius1 + b_y + balls[3].m_center.getY();
			rg_REAL z = a_z * radius1 + b_z + balls[3].m_center.getZ();
			tangentSpheres[0].m_center.setPoint(x, y, z);
			tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

			return 1;
		}
		else if (radius2 > 0) 
		{
			// backtransform with radius2
			rg_REAL x = a_x * radius2 + b_x + balls[3].m_center.getX();
			rg_REAL y = a_y * radius2 + b_y + balls[3].m_center.getY();
			rg_REAL z = a_z * radius2 + b_z + balls[3].m_center.getZ();
			tangentSpheres[1].m_center.setPoint(x, y, z);
			tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;

			return 1;
		}
		else 
		{
			return 0;
		}
	}

	rg_REAL determinant2 = (-x1 * y2 + x2 * y1)*r3 + (-x1 * y3 + x3 * y1)*r2 + (x2*y3 - x3 * y2)*r1;
	if (!rg_ZERO(determinant2))
	{
		//a_z * z^2 + b_z * z + c_z = 0

		//x = a_x * z + b_x;
		//y = a_y * z + b_y;
		//r = a_r * z + b_r;

		rg_REAL a_x = ((-y1 * z2 + y2 * z1)*r3 + (-y1 * z3 + y3 * z1)*r2 + (y2*z3 - y3 * z2)*r1) / determinant2;
		rg_REAL b_x = ((-y2 * rhs1 + rhs2 * y1)*r3 + (-y3 * rhs1 + rhs3 * y1)*r2 + (y3*rhs2 - rhs3 * y2)*r1) / determinant2;
		rg_REAL a_y = ((z2*x1 - x2 * z1)*r3 + (x1*z3 - x3 * z1)*r2 + (-x2 * z3 + x3 * z2)*r1) / determinant2;
		rg_REAL b_y = ((x2*rhs1 - rhs2 * x1)*r3 + (x3*rhs1 - rhs3 * x1)*r2 + (-x3 * rhs2 + rhs3 * x2)*r1) / determinant2;
		rg_REAL a_r = ((-y1 * z2 + y2 * z1)*x3 + (y1*z3 - y3 * z1)*x2 + (-y2 * z3 + y3 * z2)*x1) / determinant2;
		rg_REAL b_r = ((x1*y2 - x2 * y1)*rhs3 + (-x1 * y3 + x3 * y1)*rhs2 + (x2*y3 - x3 * y2)*rhs1) / determinant2;

		//plug (x, y, r) in eq4 and
		//solve the quadratic equation in z
		rg_REAL a_z = a_x * a_x + a_y * a_y + 1.0 - a_r * a_r;
		rg_REAL b_z = 2.0*a_x*b_x + 2.0*a_y*b_y - 2.0*a_r*b_r;
		rg_REAL c_z = b_x * b_x + b_y * b_y - b_r * b_r;

		rg_REAL det = b_z * b_z - 4.0*a_z*c_z;
		if (det < 0)
			return 0;
		
		rg_REAL sqrtDet = sqrt(det);
		rg_REAL zValue1 = (-b_z + sqrtDet) / (2.0*a_z);
		rg_REAL zValue2 = (-b_z - sqrtDet) / (2.0*a_z);

		rg_REAL radius1 = a_r * zValue1 + b_r;
		rg_REAL radius2 = a_r * zValue2 + b_r;

		if (radius1 > 0 && radius2 > 0) 
		{
			rg_REAL x = a_x * zValue1 + b_x + balls[3].m_center.getX();
			rg_REAL y = a_y * zValue1 + b_y + balls[3].m_center.getY();
			rg_REAL z = zValue1 + balls[3].m_center.getZ();
			tangentSpheres[0].m_center.setPoint(x, y, z);
			tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

			x = a_x * zValue2 + b_x + balls[3].m_center.getX();
			y = a_y * zValue2 + b_y + balls[3].m_center.getY();
			z = zValue2 + balls[3].m_center.getZ();
			tangentSpheres[1].m_center.setPoint(x, y, z);
			tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;

			return 2;
		}
		else if (radius1 > 0) 
		{
			rg_REAL x = a_x * zValue1 + b_x + balls[3].m_center.getX();
			rg_REAL y = a_y * zValue1 + b_y + balls[3].m_center.getY();
			rg_REAL z = zValue1 + balls[3].m_center.getZ();
			tangentSpheres[0].m_center.setPoint(x, y, z);
			tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

			return 1;
		}
		else if (radius2 > 0) 
		{
			rg_REAL x = a_x * zValue2 + b_x + balls[3].m_center.getX();
			rg_REAL y = a_y * zValue2 + b_y + balls[3].m_center.getY();
			rg_REAL z = zValue2 + balls[3].m_center.getZ();
			tangentSpheres[1].m_center.setPoint(x, y, z);
			tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;

			return 1;
		}
		else 
		{
			return 0;
		}
	}
	
	rg_REAL determinant3 = (-z2 * x1 + x2 * z1)*r3 + (-x1 * z3 + x3 * z1)*r2 + (x2*z3 - x3 * z2)*r1;
	if (!rg_ZERO(determinant3))
	{
		//a_y * y^2 + b_y * y + c_y = 0

		//x = a_x * y + b_x;		
		//z = a_z * y + b_z;
		//r = a_r * y + b_r;
		rg_REAL a_x = ((y1*z2 - y2 * z1)*r3 + (y1*z3 - y3 * z1)*r2 + (-y2 * z3 + y3 * z2)*r1) / determinant3;
		rg_REAL b_x = ((-z2 * rhs1 + rhs2 * z1)*r3 + (-z3 * rhs1 + rhs3 * z1)*r2 + (z3*rhs2 - rhs3 * z2)*r1) / determinant3;
		rg_REAL a_z = ((x1*y2 - x2 * y1)*r3 + (x1*y3 - x3 * y1)*r2 + (-x2 * y3 + x3 * y2)*r1) / determinant3;
		rg_REAL b_z = ((rhs1*x2 - rhs2 * x1)*r3 + (rhs1*x3 - rhs3 * x1)*r2 + (-rhs2 * x3 + rhs3 * x2)*r1) / determinant3;
		rg_REAL a_r = ((y1*z2 - y2 * z1)*x3 + (-y1 * z3 + y3 * z1)*x2 + (y2*z3 - y3 * z2)*x1) / determinant3;
		rg_REAL b_r = ((x1*z2 - x2 * z1)*rhs3 + (-x1 * z3 + x3 * z1)*rhs2 + (x2*z3 - x3 * z2)*rhs1) / determinant3;

		//plug (x, z, r) in eq4 and
		//solve the quadratic equation in y
		rg_REAL a_y = a_x * a_x + 1.0 + a_z * a_z - a_r * a_r;
		rg_REAL b_y = 2.0*a_x*b_x + 2.0*a_z*b_z - 2.0*a_r*b_r;
		rg_REAL c_y = b_x * b_x + b_z * b_z - b_r * b_r;

		rg_REAL det = b_y * b_y - 4.0*a_y*c_y;
		if (det < 0)
			return 0;

		rg_REAL sqrtDet = sqrt(det);
		rg_REAL yValue1 = (-b_y + sqrtDet) / (2.0*a_y);
		rg_REAL yValue2 = (-b_y - sqrtDet) / (2.0*a_y);

		rg_REAL radius1 = a_r * yValue1 + b_r;
		rg_REAL radius2 = a_r * yValue2 + b_r;
		if (radius1 > 0 && radius2 > 0) 
		{
			rg_REAL x = a_x * yValue1 + b_x + balls[3].m_center.getX();
			rg_REAL y = yValue1 + balls[3].m_center.getY();
			rg_REAL z = a_z * yValue1 + b_z + balls[3].m_center.getZ();
			tangentSpheres[0].m_center.setPoint(x, y, z);
			tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

			x = a_x * yValue2 + b_x + balls[3].m_center.getX();
			y = yValue2 + balls[3].m_center.getY();
			z = a_z * yValue2 + b_z + balls[3].m_center.getZ();
			tangentSpheres[1].m_center.setPoint(x, y, z);
			tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;

			return 2;
		}
		else if (radius1 > 0) 
		{
			rg_REAL x = a_x * yValue1 + b_x + balls[3].m_center.getX();
			rg_REAL y = yValue1 + balls[3].m_center.getY();
			rg_REAL z = a_z * yValue1 + b_z + balls[3].m_center.getZ();
			tangentSpheres[0].m_center.setPoint(x, y, z);
			tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

			return 1;
		}
		else if (radius2 > 0) 
		{
			rg_REAL x = a_x * yValue2 + b_x + balls[3].m_center.getX();
			rg_REAL y = yValue2 + balls[3].m_center.getY();
			rg_REAL z = a_z * yValue2 + b_z + balls[3].m_center.getZ();
			tangentSpheres[1].m_center.setPoint(x, y, z);
			tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;
			return 1;
		}
		else 
		{
			return 0;
		}
	}
	
	rg_REAL determinant4 = (-y1 * z2 + y2 * z1)*r3 + (-y1 * z3 + y3 * z1)*r2 + (y2*z3 - y3 * z2)*r1;
	if (!rg_ZERO(determinant4))
	{
		//a_x * x^2 + b_x * x + c_x = 0
				
		//y = a_y * x + b_y;
		//z = a_z * x + b_z;
		//r = a_r * x + b_r;
		rg_REAL a_y = ((x1*z2 - x2 * z1)*r3 + (x1*z3 - x3 * z1)*r2 + (-x2 * z3 + x3 * z2)*r1) / determinant4;
		rg_REAL b_y = ((-z2 * rhs1 + rhs2 * z1)*r3 + (-z3 * rhs1 + rhs3 * z1)*r2 + (z3*rhs2 - rhs3 * z2)*r1) / determinant4;
		rg_REAL a_z = ((-x1 * y2 + x2 * y1)*r3 + (-x1 * y3 + x3 * y1)*r2 + (x2*y3 - x3 * y2)*r1) / determinant4;
		rg_REAL b_z = ((y2*rhs1 - rhs2 * y1)*r3 + (y3*rhs1 - rhs3 * y1)*r2 + (-y3 * rhs2 + rhs3 * y2)*r1) / determinant4;
		rg_REAL a_r = ((-y1 * z2 + y2 * z1)*x3 + (y1*z3 - y3 * z1)*x2 + (-y2 * z3 + y3 * z2)*x1) / determinant4;
		rg_REAL b_r = ((y1*z2 - y2 * z1)*rhs3 + (-y1 * z3 + y3 * z1)*rhs2 + (y2*z3 - y3 * z2)*rhs1) / determinant4;

		//plug (x, z, r) in eq4 and
		//solve the quadratic equation in x
		rg_REAL a_x = 1.0 + a_y * a_y + a_z * a_z - a_r * a_r;
		rg_REAL b_x = 2.0*a_y*b_y + 2.0*a_z*b_z - 2.0*a_r*b_r;
		rg_REAL c_x = b_y * b_y + b_z * b_z - b_r * b_r;

		rg_REAL det = b_x * b_x - 4.0*a_x*c_x;
		if (det < 0)
			return 0;

		rg_REAL sqrtDet = sqrt(det);
		rg_REAL xValue1 = (-b_x + sqrtDet) / (2.0*a_x);
		rg_REAL xValue2 = (-b_x - sqrtDet) / (2.0*a_x);

		rg_REAL radius1 = a_r * xValue1 + b_r;
		rg_REAL radius2 = a_r * xValue2 + b_r;
		if (radius1 > 0 && radius2 > 0) 
		{
			rg_REAL x = xValue1 + balls[3].m_center.getX();
			rg_REAL y = a_y * xValue1 + b_y + balls[3].m_center.getY();
			rg_REAL z = a_z * xValue1 + b_z + balls[3].m_center.getZ();
			tangentSpheres[0].m_center.setPoint(x, y, z);
			tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

			x = xValue2 + balls[3].m_center.getX();
			y = a_y * xValue2 + b_y + balls[3].m_center.getY();
			z = a_z * xValue2 + b_z + balls[3].m_center.getZ();
			tangentSpheres[1].m_center.setPoint(x, y, z);
			tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;

			return 2;
		}
		else if (radius1 > 0) 
		{
			rg_REAL x = xValue1 + balls[3].m_center.getX();
			rg_REAL y = a_y * xValue1 + b_y + balls[3].m_center.getY();
			rg_REAL z = a_z * xValue1 + b_z + balls[3].m_center.getZ();
			tangentSpheres[0].m_center.setPoint(x, y, z);
			tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

			return 1;
		}
		else if (radius2 > 0) 
		{
			rg_REAL x = xValue2 + balls[3].m_center.getX();
			rg_REAL y = a_y * xValue2 + b_y + balls[3].m_center.getY();
			rg_REAL z = a_z * xValue2 + b_z + balls[3].m_center.getZ();
			tangentSpheres[1].m_center.setPoint(x, y, z);
			tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;

			return 1;
		}
		else 
		{
			return 0;
		}
	}
	
	return 0;
}

rg_INT computeSphereTangentTo4SpheresWithDiffRadiusOutside( const Sphere& s1,
                                                            const Sphere& s2, 
                                                            const Sphere& s3, 
                                                            const Sphere& s4,
                                                            Sphere* tangentSpheres)
{	
    //donguk's version (2004. 7. 5)
    //Modified by Youngsong Cho (2005. 8. 12)
    Sphere balls[4] = {s1, s2, s3, s4};
    //sort4SpheresByRadiusInDescendingPowers( balls );

    /*
	if(    (balls[0].m_radius == balls[1].m_radius) 
        && (balls[1].m_radius == balls[2].m_radius) 
        && (balls[2].m_radius == balls[3].m_radius) )
	{
        return computeSphereTangentTo4SpheresWithSameRadiusOutside( balls[0], balls[1], balls[2], balls[3], tangentSpheres);
	}
	else
	{
        if ( s1.m_radius < s2.m_radius )  {
            if ( s1.m_radius < s3.m_radius )  {
                if ( s1.m_radius < s4.m_radius )  {
                    balls[0] = s2;
                    balls[1] = s3;
                    balls[2] = s4;
                    balls[3] = s1;
                }
                else  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s3;
                    balls[3] = s4;
                }
            }
            else  {
                if ( s3.m_radius < s4.m_radius )  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s4;
                    balls[3] = s3;
                }
                else  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s3;
                    balls[3] = s4;
                }
            }
        }
        else  {
            if ( s2.m_radius < s3.m_radius )  {
                if ( s2.m_radius < s4.m_radius )  {
                    balls[0] = s1;
                    balls[1] = s3;
                    balls[2] = s4;
                    balls[3] = s2;
                }
                else  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s3;
                    balls[3] = s4;
                }
            }
            else  {
                if ( s3.m_radius < s4.m_radius )  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s4;
                    balls[3] = s3;
                }
                else  {
                    balls[0] = s1;
                    balls[1] = s2;
                    balls[2] = s3;
                    balls[3] = s4;
                }
            }
        }
    */
    {    
		rg_REAL x1 = balls[0].m_center.getX() - balls[3].m_center.getX();
		rg_REAL y1 = balls[0].m_center.getY() - balls[3].m_center.getY();
		rg_REAL z1 = balls[0].m_center.getZ() - balls[3].m_center.getZ();
		rg_REAL r1 = balls[0].m_radius        - balls[3].m_radius;

		rg_REAL x2 = balls[1].m_center.getX() - balls[3].m_center.getX();
		rg_REAL y2 = balls[1].m_center.getY() - balls[3].m_center.getY();
		rg_REAL z2 = balls[1].m_center.getZ() - balls[3].m_center.getZ();
		rg_REAL r2 = balls[1].m_radius        - balls[3].m_radius;

		rg_REAL x3 = balls[2].m_center.getX() - balls[3].m_center.getX();
		rg_REAL y3 = balls[2].m_center.getY() - balls[3].m_center.getY();
		rg_REAL z3 = balls[2].m_center.getZ() - balls[3].m_center.getZ();
		rg_REAL r3 = balls[2].m_radius        - balls[3].m_radius;


        //rg_REAL determinant1 = getDeterminantOf33Matrix(x1,y1,z1,x2,y2,z2,x3,y3,z3);
		//rg_REAL determinant2 = getDeterminantOf33Matrix(x1,y1,r1,x2,y2,r2,x3,y3,r3);
		//rg_REAL determinant3 = getDeterminantOf33Matrix(x1,r1,z1,x2,r2,z2,x3,r3,z3);
		//rg_REAL determinant4 = getDeterminantOf33Matrix(r1,y1,z1,r2,y2,z2,r3,y3,z3);

        //rg_REAL determinant1 = x1*(y2*z3 - z2*y3) - y1*(x2*z3 - z2*x3) + z1*(x2*y3 - y2*x3);
        //rg_REAL determinant2 = x1*(y2*r3 - r2*y3) - y1*(x2*r3 - r2*x3) + r1*(x2*y3 - y2*x3);
        //rg_REAL determinant3 = x1*(r2*z3 - z2*r3) - r1*(x2*z3 - z2*x3) + z1*(x2*r3 - r2*x3);
        //rg_REAL determinant4 = r1*(y2*z3 - z2*y3) - y1*(r2*z3 - z2*r3) + z1*(r2*y3 - y2*r3);


        //system of equations
		//eq1: x1*x + y1*y + z1*z = rhs1
		//eq2: x2*x + y2*y + z2*z = rhs2
		//eq3: x3*x + y3*y + z3*z = rhs3
		//eq4: x*x + y*y + z*z = r*r

		rg_REAL rhs1 = (x1*x1+y1*y1+z1*z1-r1*r1)/2.0;
		rg_REAL rhs2 = (x2*x2+y2*y2+z2*z2-r2*r2)/2.0;
		rg_REAL rhs3 = (x3*x3+y3*y3+z3*z3-r3*r3)/2.0;


		//finding non-null minor
        rg_REAL determinant1 = x1*(y2*z3 - z2*y3) - y1*(x2*z3 - z2*x3) + z1*(x2*y3 - y2*x3);
		if( !rg_ZERO(determinant1) ) 
		{
			//a_r * r^2 + b_r * r + c_r = 0

			//x = a_x * r + b_x;
			//y = a_y * r + b_y;
			//z = a_z * r + b_z;

			rg_REAL a_x = ((z2*y3-y2*z3)*r1  +(z3*y1-z1*y3)*r2  +(z1*y2-z2*y1)*r3)  /determinant1;
			rg_REAL b_x = ((y2*z3-z2*y3)*rhs1+(z1*y3-z3*y1)*rhs2+(z2*y1-z1*y2)*rhs3)/determinant1;
			rg_REAL a_y = ((z3*x2-z2*x3)*r1  +(x3*z1-x1*z3)*r2  +(z2*x1-z1*x2)*r3)  /determinant1;
			rg_REAL b_y = ((z2*x3-z3*x2)*rhs1+(x1*z3-x3*z1)*rhs2+(z1*x2-z2*x1)*rhs3)/determinant1;
			rg_REAL a_z = ((y2*x3-x2*y3)*r1  +(x1*y3-x3*y1)*r2  +(x2*y1-y2*x1)*r3)  /determinant1;
			rg_REAL b_z = ((x2*y3-y2*x3)*rhs1+(x3*y1-x1*y3)*rhs2+(y2*x1-x2*y1)*rhs3)/determinant1;


			// r = (-b_r + sqrt(b_r^2 - 4. * a_r * c_r))/(2.* a_r);
			// r = (-b_r - sqrt(b_r^2 - 4. * a_r * c_r))/(2.* a_r);
			rg_REAL a_r = a_y*a_y+a_x*a_x+a_z*a_z-1.0;
			rg_REAL b_r = 2.0*b_x*a_x+2.0*b_y*a_y+2.0*b_z*a_z;
			rg_REAL c_r = b_x*b_x+b_z*b_z+b_y*b_y;

			rg_REAL det = b_r*b_r - 4.*a_r*c_r; //b^2 - 4ac
			if(det < 0) 
				return 0;

            rg_REAL sqrtDet = sqrt(det);
			rg_REAL radius1 = (-b_r + sqrtDet)/(2.* a_r);
			rg_REAL radius2 = (-b_r - sqrtDet)/(2.* a_r);


            if( radius1 > 0 && radius2 > 0 )  {

				rg_REAL x = a_x * radius1 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y * radius1 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z * radius1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z);
                tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

				x = a_x * radius2 + b_x + balls[3].m_center.getX();
				y = a_y * radius2 + b_y + balls[3].m_center.getY();
				z = a_z * radius2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z);
                tangentSpheres[1].m_radius = radius2 - balls[3].m_radius;

				return 2;
			}
			else if( radius1 > 0 )  {

				rg_REAL x = a_x * radius1 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y * radius1 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z * radius1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z);
                tangentSpheres[0].m_radius = radius1 - balls[3].m_radius;

				return 1;
			}
			else if( radius2 > 0 )  {
				rg_REAL x = a_x * radius2 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y * radius2 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z * radius2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z);
                tangentSpheres[0].m_radius = radius2 - balls[3].m_radius;

				return 1;
			}
			else  {
				return 0;
			}

		}

        rg_REAL determinant2 = x1*(y2*r3 - r2*y3) - y1*(x2*r3 - r2*x3) + r1*(x2*y3 - y2*x3);
		if( !rg_ZERO(determinant2) ) 
		{
			//a_z * z^2 + b_z * z + c_z = 0

			//x = a_x * z + b_x;
			//y = a_y * z + b_y;
			//r = a_r * z + b_r;

			rg_REAL a_x = ((r2*y3-y2*r3)*z1  +(r3*y1-r1*y3)*z2  +(r1*y2-r2*y1)*z3)  /determinant2;
			rg_REAL b_x = ((y2*r3-r2*y3)*rhs1+(r1*y3-y1*r3)*rhs2+(r2*y1-y2*r1)*rhs3)/determinant2;
			rg_REAL a_y = ((r3*x2-r2*x3)*z1  +(x3*r1-x1*r3)*z2  +(r2*x1-r1*x2)*z3)  /determinant2;
			rg_REAL b_y = ((r2*x3-r3*x2)*rhs1+(x1*r3-x3*r1)*rhs2+(r1*x2-r2*x1)*rhs3)/determinant2;
			rg_REAL a_r = ((y2*x3-x2*y3)*z1  +(x1*y3-x3*y1)*z2  +(x2*y1-y2*x1)*z3)  /determinant2;
			rg_REAL b_r = ((x2*y3-y2*x3)*rhs1+(x3*y1-x1*y3)*rhs2+(y2*x1-x2*y1)*rhs3)/determinant2;

			rg_REAL a_z = a_x*a_x+a_y*a_y+1.0-a_r*a_r;
			rg_REAL b_z = 2.0*b_x*a_x+2.0*b_y*a_y-2.0*b_r*a_r;
			rg_REAL c_z = b_x*b_x+b_y*b_y-b_r*b_r;

			rg_REAL det = b_z*b_z - 4.0*a_z*c_z;
			if(det < 0)  
				return 0;


            rg_REAL sqrtDet = sqrt(det);
			rg_REAL zValue1 = (-b_z + sqrtDet)/(2.0*a_z);
			rg_REAL zValue2 = (-b_z - sqrtDet)/(2.0*a_z);

			rg_REAL rad1 = a_r*zValue1 + b_r;
			rg_REAL rad2 = a_r*zValue2 + b_r;

			if( rad1 > 0 && rad2 > 0)  {

				rg_REAL x = a_x*zValue1 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y*zValue1 + b_y + balls[3].m_center.getY();
                rg_REAL z = zValue1           + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				x = a_x*zValue2 + b_x + balls[3].m_center.getX();
				y = a_y*zValue2 + b_y + balls[3].m_center.getY();
                z = zValue2           + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z );
                tangentSpheres[1].m_radius = rad2 - balls[3].m_radius;

				return 2;
			}
			else if( rad1 > 0 )  {

				rg_REAL x = a_x*zValue1 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y*zValue1 + b_y + balls[3].m_center.getY();
                rg_REAL z = zValue1           + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				return 1;
			}
			else if( rad2 > 0 )  {

				rg_REAL x = a_x*zValue2 + b_x + balls[3].m_center.getX();
				rg_REAL y = a_y*zValue2 + b_y + balls[3].m_center.getY();
                rg_REAL z = zValue2           + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad2 - balls[3].m_radius;

				return 1;
			}
			else  {
				return 0;
			}
		}


        rg_REAL determinant3 = x1*(r2*z3 - z2*r3) - r1*(x2*z3 - z2*x3) + z1*(x2*r3 - r2*x3);
		if( !rg_ZERO(determinant3) ) 
		{
			//a_y * y^2 + b_y * y + c_y = 0

			//x = a_x * y + b_x;
			//r = a_r * y + b_r;
			//z = a_z * y + b_z;
			rg_REAL a_x = ((z2*r3-r2*z3)*y1  +(z3*r1-z1*r3)*y2  +(z1*r2-z2*r1)*y3)  /determinant3;
			rg_REAL b_x = ((r2*z3-z2*r3)*rhs1+(z1*r3-z3*r1)*rhs2+(z2*r1-z1*r2)*rhs3)/determinant3;
			rg_REAL a_r = ((z3*x2-z2*x3)*y1  +(x3*z1-x1*z3)*y2  +(z2*x1-z1*x2)*y3)  /determinant3;
			rg_REAL b_r = ((z2*x3-z3*x2)*rhs1+(x1*z3-x3*z1)*rhs2+(z1*x2-z2*x1)*rhs3)/determinant3;
			rg_REAL a_z = ((r2*x3-x2*r3)*y1  +(x1*r3-x3*r1)*y2  +(x2*r1-r2*x1)*y3)  /determinant3;
			rg_REAL b_z = ((x2*r3-r2*x3)*rhs1+(x3*r1-x1*r3)*rhs2+(r2*x1-x2*r1)*rhs3)/determinant3;

			rg_REAL a_y = a_x*a_x+a_z*a_z+1.0-a_r*a_r;
			rg_REAL b_y = 2.0*b_x*a_x+2.0*b_z*a_z-2.0*b_r*a_r;
			rg_REAL c_y = b_x*b_x-b_r*b_r+b_z*b_z;

			rg_REAL det = b_y*b_y - 4.0*a_y*c_y;
			if(det < 0)
				return 0;


            rg_REAL sqrtDet = sqrt(det);
			rg_REAL yValue1 = (-b_y + sqrtDet)/(2.0*a_y);
			rg_REAL yValue2 = (-b_y - sqrtDet)/(2.0*a_y);

			rg_REAL rad1 = a_r*yValue1 + b_r;
			rg_REAL rad2 = a_r*yValue2 + b_r;
			if( rad1 > 0 && rad2 > 0)  {

                rg_REAL x = a_x*yValue1 + b_x + balls[3].m_center.getX();
                rg_REAL y = yValue1           + balls[3].m_center.getY();
				rg_REAL z = a_z*yValue1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;


				x = a_x*yValue2 + b_x + balls[3].m_center.getX();
                y = yValue2           + balls[3].m_center.getY();
				z = a_z*yValue2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z );
                tangentSpheres[1].m_radius = rad2 - balls[3].m_radius;

				return 2;
			}
			else if( rad1 > 0 )  {

                rg_REAL x = a_x*yValue1 + b_x + balls[3].m_center.getX();
                rg_REAL y = yValue1           + balls[3].m_center.getY();
				rg_REAL z = a_z*yValue1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				return 1;
			}
			else if( rad2 > 0 )  {

				rg_REAL x = a_x*yValue2 + b_x + balls[3].m_center.getX();
                rg_REAL y = yValue2           + balls[3].m_center.getY();
				rg_REAL z = a_z*yValue2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z );
                tangentSpheres[1].m_radius = rad2 - balls[3].m_radius;

				return 1;
			}
            else  {
				return 0;
			}
		
		}


        rg_REAL determinant4 = r1*(y2*z3 - z2*y3) - y1*(r2*z3 - z2*r3) + z1*(r2*y3 - y2*r3);
		if( !rg_ZERO(determinant4) ) 
		{
			//a_x * x^2 + b_x * x + c_x = 0

			//r = a_r * x + b_r;
			//y = a_y * x + b_y;
			//z = a_z * x + b_z;
			rg_REAL a_r = ((z2*y3-y2*z3)*x1  +(z3*y1-z1*y3)*x2  +(z1*y2-z2*y1)*x3)  /determinant4;
			rg_REAL b_r = ((y2*z3-z2*y3)*rhs1+(z1*y3-z3*y1)*rhs2+(z2*y1-z1*y2)*rhs3)/determinant4;
			rg_REAL a_y = ((z3*r2-z2*r3)*x1  +(r3*z1-r1*z3)*x2  +(z2*r1-z1*r2)*x3)  /determinant4;
			rg_REAL b_y = ((z2*r3-z3*r2)*rhs1+(r1*z3-r3*z1)*rhs2+(z1*r2-z2*r1)*rhs3)/determinant4;
			rg_REAL a_z = ((y2*r3-r2*y3)*x1  +(r1*y3-r3*y1)*x2  +(r2*y1-y2*r1)*x3)  /determinant4;
			rg_REAL b_z = ((r2*y3-y2*r3)*rhs1+(r3*y1-r1*y3)*rhs2+(y2*r1-r2*y1)*rhs3)/determinant4;

			rg_REAL a_x = a_z*a_z+a_y*a_y+1.0-a_r*a_r;
			rg_REAL b_x = 2.0*b_y*a_y-2.0*b_r*a_r+2.0*b_z*a_z;
			rg_REAL c_x = b_z*b_z+b_y*b_y-b_r*b_r;

			rg_REAL det = b_x*b_x - 4.0*a_x*c_x;
			if(det < 0)
			{
				tangentSpheres = rg_NULL;
				return 0;
			}

            rg_REAL sqrtDet = sqrt(det);
			rg_REAL xValue1 = (-b_x + sqrtDet)/(2.0*a_x);
			rg_REAL xValue2 = (-b_x - sqrtDet)/(2.0*a_x);

			rg_REAL rad1 = a_r*xValue1 + b_r;
			rg_REAL rad2 = a_r*xValue2 + b_r;
			if( rad1 > 0 && rad2 > 0)  {

				rg_REAL x = xValue1           + balls[3].m_center.getX();
                rg_REAL y = a_y*xValue1 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z*xValue1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				x = xValue2           + balls[3].m_center.getX();
                y = a_y*xValue2 + b_y + balls[3].m_center.getY();
				z = a_z*xValue2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[1].m_center.setPoint( x, y, z );
                tangentSpheres[1].m_radius = rad2 - balls[3].m_radius;

                return 2;
			}
			else if( rad1 > 0 )  {

				rg_REAL x = xValue1           + balls[3].m_center.getX();
                rg_REAL y = a_y*xValue1 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z*xValue1 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad1 - balls[3].m_radius;

				return 1;
			}
			else if( rad2 > 0 )  {

				rg_REAL x = xValue2           + balls[3].m_center.getX();
                rg_REAL y = a_y*xValue2 + b_y + balls[3].m_center.getY();
				rg_REAL z = a_z*xValue2 + b_z + balls[3].m_center.getZ();
				tangentSpheres[0].m_center.setPoint( x, y, z );
                tangentSpheres[0].m_radius = rad2 - balls[3].m_radius;

                return 1;
			}
			else  {
				return 0;
			}
		}
		
        
        return 0;
	}
	
    return 0;
}




rg_INT computeSphereTangentTo4SpheresWithSameRadiusOutside( const Sphere& s1, 
                                                            const Sphere& s2, 
                                                            const Sphere& s3, 
                                                            const Sphere& s4,
                                                            Sphere* tangentSpheres)
{
    /*
	Sphere* arrSphere[4] = {s1, s2, s3, s4};
	
	// initilizing variables for calculating Vertex
	// 3 spheres are shrinked and translated by minimum sphere
    rg_Point3D centerOfSphere[4];
    for (int i=0; i<4; i++)  {
        centerOfSphere[i] = arrSphere[i]->getCenter();
    }
    */

	double x1 = s1.m_center.getX() - s4.m_center.getX();
	double y1 = s1.m_center.getY() - s4.m_center.getY();
	double z1 = s1.m_center.getZ() - s4.m_center.getZ();
	//double r1 = s1.m_radius        - s4.m_radius;

	double x2 = s2.m_center.getX() - s4.m_center.getX();
	double y2 = s2.m_center.getY() - s4.m_center.getY();
	double z2 = s2.m_center.getZ() - s4.m_center.getZ();
	//double r2 = s2.m_radius        - s4.m_radius;

	double x3 = s3.m_center.getX() - s4.m_center.getX();
	double y3 = s3.m_center.getY() - s4.m_center.getY();
	double z3 = s3.m_center.getZ() - s4.m_center.getZ();
	//double r3 = s3.m_radius        - s4.m_radius;

    /*
	double det_A =   ( x1 * ( ( y2 * z3 ) - ( z2 * y3 ) ) )
	 			   - ( y1 * ( ( x2 * z3 ) - ( z2 * x3 ) ) )
				   + ( z1 * ( ( x2 * y3 ) - ( y2 * x3 ) ) );
    */
    double det_A = rg_GeoFunc::getDeterminantOf33Matrix(x1, y1, z1, x2, y2, z2, x3, y3, z3);

		
	if ( det_A == 0 )
		return 0;

	// evalueate Final D-Sphere's X, Y, Z coordinates and radius
	double l1 = (x1*x1) + (y1*y1) + (z1*z1);
	double l2 = (x2*x2) + (y2*y2) + (z2*z2);
	double l3 = (x3*x3) + (y3*y3) + (z3*z3);

    double finalX = rg_GeoFunc::getDeterminantOf33Matrix(l1, y1, z1, l2, y2, z2, l3, y3, z3)/(2*det_A) + s4.m_center.getX();
	double finalY = rg_GeoFunc::getDeterminantOf33Matrix(x1, l1, z1, x2, l2, z2, x3, l3, z3)/(2*det_A) + s4.m_center.getY();
	double finalZ = rg_GeoFunc::getDeterminantOf33Matrix(x1, y1, l1, x2, y2, l2, x3, y3, l3)/(2*det_A) + s4.m_center.getZ();

    rg_Point3D center(finalX, finalY, finalZ);

	double radius = center.distance( s1.m_center ) - s4.m_radius;

    tangentSpheres[0].m_center = center;
    tangentSpheres[0].m_radius = radius;

	return 1;
}




rg_INT computeSphereWithItsCenterOnPlaneAndTangentTo3SpheresOutside( 
                                                         rg_REAL*      plane,
                                                         const Sphere& s1, 
                                                         const Sphere& s2, 
                                                         const Sphere& s3, 
                                                         Sphere* tangentSpheres)
{
    Sphere balls[3];
    if ( rg_LT( s1.m_radius, s2.m_radius ) )  {
        if ( rg_LT( s1.m_radius, s3.m_radius ) )  {
            balls[0] = s2;
            balls[1] = s3;
            balls[2] = s1;
        }
        else  {
            balls[0] = s1;
            balls[1] = s2;
            balls[2] = s3;
        }
    }
    else  {
        if ( rg_LT( s2.m_radius, s3.m_radius ) )  {
            balls[0] = s1;
            balls[1] = s3;
            balls[2] = s2;
        }
        else  {
            balls[0] = s1;
            balls[1] = s2;
            balls[2] = s3;
        }
    }

	rg_REAL x1 = balls[0].m_center.getX() - balls[2].m_center.getX();
	rg_REAL y1 = balls[0].m_center.getY() - balls[2].m_center.getY();
	rg_REAL z1 = balls[0].m_center.getZ() - balls[2].m_center.getZ();
	rg_REAL r1 = balls[0].m_radius        - balls[2].m_radius;

	rg_REAL x2 = balls[1].m_center.getX() - balls[2].m_center.getX();
	rg_REAL y2 = balls[1].m_center.getY() - balls[2].m_center.getY();
	rg_REAL z2 = balls[1].m_center.getZ() - balls[2].m_center.getZ();
	rg_REAL r2 = balls[1].m_radius        - balls[2].m_radius;

    //system of equations after translating the center of s3 to the origin
	//eq1: (x-x1)^2 + (y-y1)^2 + (z-z1)^2 = (r+r1)^2
	//eq2: (x-x2)^2 + (y-y2)^2 + (z-z2)^2 = (r+r2)^2
	//eq3: x^2      + y^2      + z^2      = r^2

    //system of equations
	//eq1-eq3: x1*x + y1*y = rhs1 - z1*z - r1*r (rhs1 = (x1^2 + y1^2 + z1^2 + r3^2 - r1^2)/2)
	//eq2-eq3: x2*x + y2*y = rhs2 - z2*z - r2*r (rhs2 = (x2^2 + y2^2 + z2^2 + r3^2 - r2^2)/2)

	rg_REAL rhs1 = (x1*x1 + y1*y1 + z1*z1 - r1*r1)*0.5;
	rg_REAL rhs2 = (x2*x2 + y2*y2 + z2*z2 - r2*r2)*0.5;


    rg_INT numTS = 0;
    rg_REAL determinant1 = x1*y2 - x2*y1;
    if ( rg_NZERO(determinant1) )  {
        //  x = a_x + b_x*z + c_x*r
        //  y = a_y + b_y*z + c_y*r
        rg_REAL a_x = (y2*rhs1-y1*rhs2);
        rg_REAL b_x = (-y2*z1+y1*z2);
        rg_REAL c_x = (-y2*r1+y1*r2);
        rg_REAL a_y = (x1*rhs2-x2*rhs1);
        rg_REAL b_y = (-x1*z2+x2*z1);
        rg_REAL c_y = (-x1*r2+x2*r1);

        //  a_z*z^2 + 2*(b_z*r + c_z)*z + d_z*r*r + e_z*r + f_z
        rg_REAL a_z = (b_x*b_x + b_y*b_y + determinant1*determinant1);
        rg_REAL b_z = (b_x*c_x + b_y*c_y);
        rg_REAL c_z = (a_y*b_y + a_x*b_x);
        rg_REAL d_z = (c_x*c_x + c_y*c_y - determinant1*determinant1);
        rg_REAL e_z = 2.0*(a_y*c_y + a_x*c_x);
        rg_REAL f_z = a_x*a_x + a_y*a_y;

        // a_r*r + b_r = c_r*sqrt      
        rg_REAL a_r = (-plane[1]*a_y*b_z/a_z - plane[2]*b_z/a_z + plane[0]*b_x + plane[1]*b_y - plane[0]*a_x*b_z/a_z);
        rg_REAL b_r = plane[3] - plane[0]*a_x*c_z/a_z - plane[1]*a_y*c_z/a_z - plane[2]*c_z/a_z + plane[0]*c_x + plane[1]*c_y;
        rg_REAL c_r = plane[0]*a_x + plane[1]*a_y + plane[2];

        // (a_r*r + b_r)^2 = (c_r*sqrt)^2
        // -> a_R*r^2 + b_R*r + c_R = 0
        rg_REAL a_R = (c_r*c_r*a_z*d_z-c_r*c_r*b_z*b_z+a_r*a_r);
        rg_REAL b_R = 2.0*(c_r*c_r*a_z*e_z-c_r*c_r*b_z*c_z+a_r*b_r);
        rg_REAL c_R = c_r*c_r*a_z*f_z-c_r*c_r*c_z*c_z+b_r*b_r;

        rg_REAL detR = b_R*b_R - 4.0*a_R*c_R;

        if ( detR < 0.0 )
            return 0;

        rg_REAL sqrtDetR = sqrt(detR);
        rg_REAL radius[2];
        radius[0] = (-b_R + sqrtDetR)/(2.0*a_R);
        radius[1] = (-b_R - sqrtDetR)/(2.0*a_R);
        
        if ( radius[0] >= 0.0 )  {
            rg_REAL detZ = (b_z*radius[0] + c_z)*(b_z*radius[0] + c_z)
                           - a_z*(d_z*radius[0]*radius[0] + e_z*radius[0] + f_z);

            if ( rg_ZERO( detZ ) ) {
                rg_REAL z = (-(b_z*radius[0] + c_z) )/a_z;

                rg_REAL x = a_x + b_x*z + c_x*radius[0];
                rg_REAL y = a_y + b_y*z + c_y*radius[0];

            }
            else if ( rg_POS( detZ ) ) {
                rg_REAL sqrtDetZ = sqrt(detZ);

                rg_REAL x[2];
                rg_REAL y[2];
                rg_REAL z[2];

                z[0] = (-(b_z*radius[0] + c_z) + sqrtDetZ )/a_z;
                z[1] = (-(b_z*radius[0] + c_z) - sqrtDetZ )/a_z;

                x[0] = a_x + b_x*z[0] + c_x*radius[0];
                x[1] = a_x + b_x*z[1] + c_x*radius[0];

                y[0] = a_y + b_y*z[0] + c_y*radius[0];
                y[1] = a_y + b_y*z[1] + c_y*radius[0];

            }
            else {}
        }

    }

    return 0;
}





void computePlaneTangentTo3SphereFromCCWOutside(Sphere* generator, rg_REAL* tangentPlane)
{
//  int i=0;
//
//	Sphere generator[3];
//	for(i=0; i<3; i++)
//	{
//		generator[i] = ((BallGeneratorCore*)cells[i]->getGenerator())->getBall();
//	}

	//make the index of the smallest generator  0.
	rg_REAL radius0 = generator[0].getRadius();
	rg_REAL radius1 = generator[1].getRadius();
	rg_REAL radius2 = generator[2].getRadius();

	if( radius1 <= radius2 && radius1 <= radius0 )
	{
		Sphere temp = generator[0];
		generator[0] = generator[1];
		generator[1] = generator[2];
		generator[2] = temp;
	}
	else if( radius2 <= radius1 && radius2 <= radius0 )
	{
		Sphere temp = generator[0];
		generator[0] = generator[2];
		generator[2] = generator[1];
		generator[1] = temp;
	}
	

	//normal vector of the plane passing through three centers of generators
	//toward infinity vertex
	rg_Point3D vec1 = generator[1].getCenter() - generator[0].getCenter();
	rg_Point3D vec2 = generator[2].getCenter() - generator[1].getCenter();
	rg_Point3D normalVec = vec1.crossProduct(vec2).getUnitVector();

	rg_REAL x1 = generator[0].getCenter().getX();
	rg_REAL y1 = generator[0].getCenter().getY();
	rg_REAL z1 = generator[0].getCenter().getZ();
	rg_REAL r1 = generator[0].getRadius();

	rg_REAL x2 = generator[1].getCenter().getX();
	rg_REAL y2 = generator[1].getCenter().getY();
	rg_REAL z2 = generator[1].getCenter().getZ();
	rg_REAL r2 = generator[1].getRadius();

	rg_REAL x3 = generator[2].getCenter().getX();
	rg_REAL y3 = generator[2].getCenter().getY();
	rg_REAL z3 = generator[2].getCenter().getZ();
	rg_REAL r3 = generator[2].getRadius();

	bool bPerturb_x = false;
	bool bPerturb_y = false;
	bool bPerturb_z = false;
	rg_REAL den = (y1*x3*z2+x1*y2*z3+x2*y3*z1-y3*x1*z2-x3*y2*z1-x2*y1*z3);
	if(rg_ZERO(den))
	{
		bPerturb_x = true;

		x1 = x1 + 1;
		x2 = x2 + 1;
		x3 = x3 + 1;
		den = (y1*x3*z2+x1*y2*z3+x2*y3*z1-y3*x1*z2-x3*y2*z1-x2*y1*z3);
	}
	if(rg_ZERO(den))
	{
		bPerturb_x = false;
		bPerturb_y = true;
		
		x1 = x1 - 1;
		x2 = x2 - 1;
		x3 = x3 - 1;

		y1 = y1 + 1;
		y2 = y2 + 1;
		y3 = y3 + 1;
		den = (y1*x3*z2+x1*y2*z3+x2*y3*z1-y3*x1*z2-x3*y2*z1-x2*y1*z3);
	}
	if(rg_ZERO(den))
	{
		bPerturb_y = false;
		bPerturb_z = true;

		y1 = y1 - 1;
		y2 = y2 - 1;
		y3 = y3 - 1;

		z1 = z1 + 1;
		z2 = z2 + 1;
		z3 = z3 + 1;
		den = (y1*x3*z2+x1*y2*z3+x2*y3*z1-y3*x1*z2-x3*y2*z1-x2*y1*z3);
	}

	rg_REAL a1 = -(-y3*z2+y2*z3-y2*z1+y1*z2+y3*z1-y1*z3);
	rg_REAL a2 = -(r2*y3*z1-r2*y1*z3+y2*z3*r1-y2*r3*z1-z2*y3*r1+z2*y1*r3);
	rg_REAL b1 = x2*z3-x2*z1-x3*z2+x1*z2-x1*z3+x3*z1;
	rg_REAL b2 = x2*z3*r1-x2*r3*z1+x3*r2*z1-x3*z2*r1-z3*x1*r2+r3*x1*z2;
	rg_REAL c1 = -(-x3*y2-x2*y1+x2*y3+y1*x3+x1*y2-y3*x1);
	rg_REAL c2 = -(x2*y3*r1-x2*y1*r3-y3*x1*r2-x3*y2*r1+x1*y2*r3+y1*x3*r2);

	rg_REAL d = 1/(c1*c1+b1*b1+a1*a1)*
			(-2.0*c1*c2-2.0*b1*b2-2.0*a1*a2+2.0*
			 sqrt(2.0*c1*c2*b1*b2+2.0*c1*c2*a1*a2+
				 2.0*b1*b2*a1*a2-c1*c1*b2*b2-c1*c1*a2*a2+
				 c1*c1*den*den-b1*b1*c2*c2-b1*b1*a2*a2+
				 b1*b1*den*den-a1*a1*c2*c2-a1*a1*b2*b2+
				 a1*a1*den*den))/2.0;

	rg_REAL a = (a1*d+a2)/den;
	rg_REAL b = (b1*d+b2)/den;
	rg_REAL c = (c1*d+c2)/den;

    // vector (a,b,c)가 세 구의 중심으로 정의되는 vector와 방향이 일치하는
    // 지를 검사하여야 한다.
	if( rg_Point3D(a,b,c).innerProduct(normalVec) < 0 )
	{

		d = 1/(c1*c1+b1*b1+a1*a1)*(-2.0*c1*c2-2.0*b1*b2-2.0*a1*a2-2.0*
			sqrt(2.0*c1*c2*b1*b2+2.0*c1*c2*a1*a2+2.0*b1*b2*a1*a2-
				c1*c1*b2*b2-c1*c1*a2*a2+c1*c1*den*den-
				b1*b1*c2*c2-b1*b1*a2*a2+b1*b1*den*den-
				a1*a1*c2*c2-a1*a1*b2*b2+a1*a1*den*den))/2.0;

		a = (a1*d+a2)/den;
		b = (b1*d+b2)/den;
		c = (c1*d+c2)/den;
	}

	if(bPerturb_x)
	{
		d = d + a;
	}
	else if(bPerturb_y)
	{
		d = d + b;
	}
	else if(bPerturb_z)
	{
		d = d + c;
	}
    tangentPlane[0] = a;
    tangentPlane[1] = b;
    tangentPlane[2] = c;
    tangentPlane[3] = d;
}



void sort3SpheresByRadiusInDescendingPowers( Sphere* spheres )
{
    if ( spheres[0] < spheres[1] )
        rg_SWAP(spheres[0], spheres[1]);

    if ( spheres[1] < spheres[2] )
        rg_SWAP(spheres[1], spheres[2]);

    if ( spheres[0] < spheres[1] )
        rg_SWAP(spheres[0], spheres[1]);
}

void sort4SpheresByRadiusInDescendingPowers( Sphere* spheres )
{
    if ( spheres[0] < spheres[1] )
        rg_SWAP(spheres[0], spheres[1]);


    if ( spheres[2] < spheres[3] )
        rg_SWAP(spheres[2], spheres[3]);


    if ( spheres[0] > spheres[2] )
    {
        if ( spheres[1] < spheres[2] )
        {
            rg_SWAP(spheres[1], spheres[2]);                

            if ( spheres[2] < spheres[3] )
                rg_SWAP(spheres[2], spheres[3]);                
        }
    }
    else 
    {
        rg_SWAP( spheres[0], spheres[2]);                
        rg_SWAP( spheres[1], spheres[3]);                

        if ( spheres[1] < spheres[2] )
        {
            rg_SWAP(spheres[1], spheres[2]);                

            if ( spheres[2] < spheres[3] )
                rg_SWAP(spheres[2], spheres[3]);                
        }
    }
}


     
vector<Sphere> Sphere::computeMinimumSphereTangentTo3Balls(const Sphere& ball1, const Sphere& ball2, const Sphere& ball3)
{
    rg_Point3D center[3] = {ball1.getCenter(), ball2.getCenter(), ball3.getCenter()};
    rg_REAL    radius[3] = {ball1.getRadius(), ball2.getRadius(), ball3.getRadius()};

    Plane centerPlane( center[0], center[1], center[2] );
    rg_Point3D normalOfCenterPlane = centerPlane.getNormal();

	rg_TMatrix3D trMatrix;
	rg_TMatrix3D invMatrix;

	trMatrix.translate(-center[0]);

    rg_FLAG bReverse = rg_EQ( normalOfCenterPlane.innerProduct( rg_Point3D(0.0, 0.0, 1.0) ), -1);
    if (!bReverse) {
        trMatrix.rotate(normalOfCenterPlane, rg_Point3D(0.0, 0.0, 1.0));
    }
    else {
        trMatrix.rotateY(rg_PI);
    }

    if (!bReverse) {
        invMatrix.rotate(rg_Point3D(0.0, 0.0, 1.0), normalOfCenterPlane);
    }
    else {
        invMatrix.rotateY(-rg_PI);
    }

	invMatrix.translate( center[0] );


    rg_Point3D  trGateCenter[3] = { trMatrix*center[0], trMatrix*center[1], trMatrix*center[2] };
	rg_Circle2D transformedCircle1( trGateCenter[0].getX(), trGateCenter[0].getY(), radius[0] );
	rg_Circle2D transformedCircle2( trGateCenter[1].getX(), trGateCenter[1].getY(), radius[1] );
	rg_Circle2D transformedCircle3( trGateCenter[2].getX(), trGateCenter[2].getY(), radius[2] );

	rg_INT numOfCC = 0;
	rg_Circle2D tangentCircle[2];
    numOfCC = computeCircleTangentTo3CirclesOutside( transformedCircle1, transformedCircle2, 
                                                            transformedCircle3, tangentCircle);   




    vector<Sphere> minimumTangentSpheres;
    if (numOfCC == 1)  {      
        minimumTangentSpheres.push_back( Sphere( invMatrix*tangentCircle[0].getCenterPt(), tangentCircle[0].getRadius() ) );
    }
    else if (numOfCC == 2){
        if ( tangentCircle[0].getRadius() < tangentCircle[1].getRadius() )  {
            minimumTangentSpheres.push_back( Sphere( invMatrix*tangentCircle[0].getCenterPt(), tangentCircle[0].getRadius() ) );
            minimumTangentSpheres.push_back( Sphere( invMatrix*tangentCircle[1].getCenterPt(), tangentCircle[1].getRadius() ) );
        }
        else  {
            minimumTangentSpheres.push_back( Sphere( invMatrix*tangentCircle[1].getCenterPt(), tangentCircle[1].getRadius() ) );
            minimumTangentSpheres.push_back( Sphere( invMatrix*tangentCircle[0].getCenterPt(), tangentCircle[0].getRadius() ) );
        }
    }
    else  {}

    return minimumTangentSpheres;
}


