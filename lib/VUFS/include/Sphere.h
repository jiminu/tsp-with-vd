#ifndef _SPHERE_H
#define _SPHERE_H

#include "rg_Const.h"
#include "rg_Point3D.h"
#include "Plane.h"
#include "Circle3D.h"
#include "LineSegment3D.h"
#include "Arc3D.h"

#include <vector>
using namespace std;



enum SphereIntersectionType {SI_NOT_INTERSECT, SI_INTERSECT, SI_CONTAIN, SI_CONTAINED};

enum FourSpheresIntersectionType {SI_QUADRUPLE_NO_INTERSECTION, 
                                  SI_QUADRUPLE_COMMON_INTERSECTION,
                                  SI_QUADRUPLE_COMMON_INTERSECTION_YET_IGNORABLE};

class Sphere
{
protected:
    rg_Point3D m_center;
    rg_REAL    m_radius;

public:
    //  constructor & deconstructor..
    Sphere();
    Sphere(const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, const rg_REAL& radius);
    Sphere(const rg_Point3D& center, const rg_REAL& radius);
    Sphere(const Sphere& aSphere);
    ~Sphere();

    //  get functions.. 
    inline rg_Point3D getCenter() const {return m_center;};
    inline rg_REAL    getRadius() const {return m_radius;};

    //  set functions..
    void setCenter(const rg_REAL& x, const rg_REAL& y, const rg_REAL& z);
    void setCenter(const rg_Point3D& center);
    void setRadius(const rg_REAL& radius);
	inline void resizeRadiusBy(const rg_REAL& delta) {m_radius += delta;}
    void setSphere(const rg_Point3D& center, const rg_REAL& radius);
    void setSphere(const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, const rg_REAL& radius);

    //  operator overloading..
    Sphere& operator =(const Sphere& aSphere);
    bool operator==(const Sphere& aSphere) const;
    bool operator <(const Sphere& aSphere) const;
    bool operator >(const Sphere& aSphere) const;

    //  geometric operators..
    rg_REAL distance(const Sphere& aSphere) const;
    rg_REAL distance(const rg_Point3D& point) const;
    rg_REAL distance(rg_REAL* plane) const;
    rg_REAL distanceToSphere(const Sphere& aSphere) const;
    rg_REAL distanceToPlane(rg_REAL* plane) const;
    rg_REAL distanceBetweenCenterAndPlane(rg_REAL* plane) const;

	rg_Point3D project(const rg_Point3D& point) const;
		
    rg_BOOL isEqual(const Sphere& aSphere, const rg_REAL& tolerance = resNeg6) const;

    //  area and volume..
	// jhryu
	inline rg_REAL computeArea() const {return 4.0*rg_PI*m_radius*m_radius;}
	//rg_REAL computeArea() const;
	rg_REAL computeAreaOfSphericalTriangle(const rg_Point3D& pt1OnSphere, const rg_Point3D& pt2OnSphere, const rg_Point3D& pt3OnSphere) const;
	rg_REAL computeAreaOfSphericalCap(const rg_REAL& h) const;
	rg_REAL computeAreaOfSphericalCap(const Plane& plane) const;
	// jhryu
	rg_REAL computeAreaOfSphericalCap(const Sphere& sphere) const;
	rg_REAL computeIntersectionArea(const Sphere& sphere) const;
    rg_REAL ratioOfOverlappedSphere(const Sphere& sphere) const;
	// jhryu
	inline rg_REAL computeVolume() const {return 4.0*rg_PI*m_radius*m_radius*m_radius/3.0;}
	//rg_REAL computeVolume() const;
	rg_REAL computeVolumeOfSphericalTriangle(const rg_Point3D& pt1OnSphere, const rg_Point3D& pt2OnSphere, const rg_Point3D& pt3OnSphere) const;
	rg_REAL computeVolumeOfSphericalCap(const rg_REAL& h) const;
	rg_REAL computeVolumeOfSphericalCap(const Plane& plane) const;
	rg_REAL computeVolumeBySphericalCapAndPlanePassedCenter(const Plane& planeForCap, const Plane& planeOnCenter) const;

	rg_REAL computeIntersectingVolumeA(const Sphere& aSphere) const;
    rg_REAL computeIntersectingVolume(const Sphere& aSphere) const;
	rg_REAL computeIntersectingVolumeB(const Sphere& aSphere) const;
    rg_REAL computeIntersectingVolumeWithSameRadiusSpheres(const Sphere& sphere1, const Sphere& sphere2) const;
	rg_REAL computeIntersectionVolumeWith(const Sphere& sphere1, const Sphere& sphere2) const;
    rg_REAL computeSumPairwiseXVolOfBallWithBallSet(rg_dList<Sphere>& ballSet) const;
	rg_REAL computeIntersectionVolumeWithFourSpheresOnlyIfThereIsCommonIntersection(const Sphere& sphere1, const Sphere& sphere2, const Sphere& sphere3) const;
		rg_REAL    computeW_SquaredForTripleIntersection( const rg_REAL& a_squared,  const rg_REAL& b_squared,  const rg_REAL& c_squared,
                                                          const rg_REAL& r1_squared, const rg_REAL& r2_squared, const rg_REAL& r3_squared ) const;

		rg_REAL    computeW_SquaredForQuadrupleIntersection( const rg_REAL& a_squared,  const rg_REAL& b_squared,  const rg_REAL& c_squared,
		                                                     const rg_REAL& f_squared,  const rg_REAL& g_squared,  const rg_REAL& h_squared,
		  										             const rg_REAL& r1_squared, const rg_REAL& r2_squared, const rg_REAL& r3_squared, const rg_REAL& r4_squared ) const;

		friend FourSpheresIntersectionType    parseTypeOfIntersectionAmongFourSpheres(const Sphere& sphere1, 
                                                                                      const Sphere& sphere2, 
                                                                                      const Sphere& sphere3, 
                                                                                      const Sphere& sphere4 );

		friend void                           doContainmentTestOfTripletSpheresIntersectionForFourthSphere(const Sphere& sphere1, 
                                                                                                           const Sphere& sphere2, 
                                                                                                           const Sphere& sphere3, 
                                                                                                           const Sphere& sphere4,
																										   rg_INT intersectionPointInclusion[][2],
																										   rg_Point3D intersectionPoint[][2],
																										   rg_INT numIntersectionPoint[]);



    bool    hasIntersectionWith(const LineSegment3D& lineSegment) const;

    rg_FLAG isThereIntersectionWith(const Sphere& anotherSphere) const;
    rg_FLAG isContainedIn(const Sphere& anotherSphere) const;
	rg_BOOL doesContainIntersectionBTW(const Sphere& sphere1, const Sphere& sphere2) const;
	rg_BOOL doesContainIntersectionAmong(const Sphere& sphere1, const Sphere& sphere2, const Sphere& sphere3) const;
    rg_BOOL doesContain(const rg_Point3D& point) const;
	rg_BOOL isOnSphere(const rg_Point3D& point) const;
	rg_BOOL isIncludedIn(const Sphere& anotherSphere) const;
	rg_BOOL isIncludedInUnionOfSpheres(const Sphere& sphere1, const Sphere& sphere2) const;
	rg_BOOL isIncludedInUnionOfSpheres(const Sphere& sphere1, const Sphere& sphere2, const Sphere& sphere3) const;

    void   makeContainingSphereWith(const Sphere& anotherSphere);
    Sphere computeContainingSphereOfThisAnd(const Sphere& anotherSphere) const;

    rg_INT                 intersect(const LineSegment3D& lineSegment, rg_Point3D* intersectionPoint) const;
	rg_INT                 intersect(const Circle3D& cirlce, rg_Point3D* intersectionPoint) const;
    rg_INT                 intersect(const Arc3D& arc, rg_Point3D* intersectionPoint) const;
    rg_BOOL                intersect(const Plane& plane, Circle3D& intersectingCircle) const;
    SphereIntersectionType intersect(const Sphere& sphere, Circle3D& intersectingCircle) const;
	SphereIntersectionType intersect(const Sphere& sphere) const;

	friend rg_INT computeIntersectionPointsAmong(const Sphere& sphere1, 
		                                         const Sphere& sphere2, 
												 const Sphere& sphere3,
												 rg_Point3D* intersectionPoint);

	friend rg_BOOL isIntersectionOfFirstTwoIncludedInUnionOfLastTwo(const Sphere& sphere1, 
                                                                    const Sphere& sphere2, 
                                                                    const Sphere& sphere3, 
                                                                    const Sphere& sphere4);

	friend rg_INT computeSphereWithGivenRadiusAndTangentTo3Spheres(
                                                 const rg_REAL& radius,
                                                       Sphere*  sphere, 
                                                       Sphere*  tangentSphere );
    friend void   computePlaneTangentTo3SphereFromCCWOutside(Sphere*  generator, 
                                                           rg_REAL* tangentPlane);

    friend rg_INT computeSphereTangentTo4SpheresOutside( const Sphere& s1, 
                                                         const Sphere& s2, 
                                                         const Sphere& s3, 
                                                         const Sphere& s4,
                                                         Sphere* tangentSpheres);

    //  for test of numerical error.



	friend rg_INT compute_tangent_sphere_of_three_ball_generators_outside_and_spherical_container_inside(const Sphere& b1,
		                                                                                                 const Sphere& b2,
		                                                                                                 const Sphere& b3,
		                                                                                                 const Sphere& container,
		                                                                                                 Sphere* tangentSpheres);


    friend rg_INT computeSphereTangentTo4SpheresWithDiffRadiusOutside( const Sphere& s1, 
                                                                       const Sphere& s2, 
                                                                       const Sphere& s3, 
                                                                       const Sphere& s4,
                                                                       Sphere* tangentSpheres);

    friend rg_INT computeSphereTangentTo4SpheresWithSameRadiusOutside( 
                                                         const Sphere& s1, 
                                                         const Sphere& s2, 
                                                         const Sphere& s3, 
                                                         const Sphere& s4,
                                                         Sphere* tangentSpheres);
    
    friend rg_INT computeSphereWithItsCenterOnPlaneAndTangentTo3SpheresOutside( 
                                                         rg_REAL*      plane,
                                                         const Sphere& s1, 
                                                         const Sphere& s2, 
                                                         const Sphere& s3, 
                                                         Sphere* tangentSpheres);

    friend void sort3SpheresByRadiusInDescendingPowers( Sphere* spheres );
    friend void sort4SpheresByRadiusInDescendingPowers( Sphere* spheres );


    static vector<Sphere> computeMinimumSphereTangentTo3Balls(const Sphere& ball1, const Sphere& ball2, const Sphere& ball3);

};


#endif
