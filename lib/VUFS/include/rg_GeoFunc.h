#ifndef _GEOFUNC_IS_INCLUDED
#define _GEOFUNC_IS_INCLUDED

#include "rg_Const.h"
#include "rg_Point2D.h"
#include "rg_Point3D.h"
#include "rg_Line.h"
#include "Sphere.h"
#include "Plane.h"

#include "rg_Line2D.h"
#include "Ellipse2D.h"
#include "rg_Circle2D.h"


template<class T> void rg_SWAP( T& item1, T& item2  );

template<class T> 
void rg_SWAP( T& item1, T& item2  )
{
    T tempItem( item1 );

    item1 = item2;
    item2 = tempItem;
}




class rg_GeoFunc
{
public:
    static rg_Line<rg_Point2D>  lineOffset(const rg_Line<rg_Point2D> &line, const rg_REAL &offset);
    static rg_Point2D           getNormal(const rg_Point2D &point);
    static rg_REAL              distanceLineVsPointl(const rg_Line<rg_Point2D> &line, rg_Point2D &point);

    static rg_Point3D*        samplePtsOnArc( const rg_Point3D& center,
                                  const rg_Point3D& localXAxis,
							      const rg_Point3D& localYAxis,
							      const rg_REAL&  radius,
							      const rg_REAL&  startDegree,
							      const rg_REAL&  endDegree  ,
							      const rg_INT&   size );
    static rg_Point3D         getUnitNormal( const rg_Point3D& pt1,
                                   const rg_Point3D& pt2,
							       const rg_Point3D& pt3 );


	static inline bool left_turn(const rg_Point2D& firstPt, const rg_Point2D& secondPt, const rg_Point2D& thirdPt)
	{
		rg_Point2D first2Second = secondPt - firstPt;
		rg_Point2D first2Thrid = thirdPt - firstPt;
		
		if (first2Second * first2Thrid > 0)
			return true;
		else
			return false;
	}


	static inline rg_REAL compute_orientation(const rg_Point2D & vtx1, const rg_Point2D & vtx2, const rg_Point2D & vtx3)
	{
		rg_REAL x1 = vtx1.getX();
		rg_REAL y1 = vtx1.getY();
		rg_REAL x2 = vtx2.getX();
		rg_REAL y2 = vtx2.getY();
		rg_REAL x3 = vtx3.getX();
		rg_REAL y3 = vtx3.getY();

		rg_REAL determinant = (x2 * y3 - x3 * y2)
			- (x1 * y3 - x3 * y1)
			+ (x1 * y2 - x2 * y1);

		return determinant;
	}

    //                   _* pt3
    //                   /|
    //                  / 
    //                 / 
    //    *---------->*  
    //   pt1           pt2
    static rg_REAL calculateAngle(const rg_Point3D& pt1, 
                                  const rg_Point3D& pt2,
				                  const rg_Point3D& pt3 );

    //       _* vect2
    //       /|
    //      / 
    //     / 
    //    *---------->*  
    //  origin        vect1
    static rg_REAL calculateAngle(const rg_Point3D& vect1, 
                                  const rg_Point3D& vect2 );

    static rg_REAL calculateAngle(const rg_Point2D& vect1,
                                  const rg_Point2D& vect2);

    static inline rg_REAL angleFromVec1toVec2(const rg_Point2D& vector1, const rg_Point2D& vector2, const rg_REAL& tolerance)
    {
        rg_Point2D vec1 = vector1.getUnitVector();
        rg_Point2D vec2 = vector2.getUnitVector();
        rg_Point2D pt2pt(vec1 - vec2);
        rg_REAL length = pt2pt.magnitude();

        rg_REAL cosine = (2. - length * length) / (2.);

        //  Lee, mokwon modify the following process at 2017.06.24.
        //  The input parameter of acos function should be in the interval [-1, +1].
        if (cosine < -1.0) {
            cosine = -1.0;
        }

        if (cosine > 1.0) {
            cosine = 1.0;
        }
        /////////////////////////////////////////////////////////////////////////////

        //    if( vec1 * vec2 > 0.0 )  //less than PI (cross product)
        if (rg_NNEG(vec1 * vec2, tolerance))  //less than PI (cross product)
        {
            return acos(cosine);
        }
        else
        {
            return 2.*rg_PI - acos(cosine);
        }
    }

    static rg_Point2D angle_bisecting_vector(const rg_Point2D& vector1, const rg_Point2D& vector2, const rg_REAL& tolerance);

    static rg_REAL getNormalizedAngle(const rg_REAL &angle);
    // connvert angles 
    //   such that 0<= start < end < 4*PI
    static void convertNormalizedAngle(rg_REAL& start,
                                       rg_REAL& end );
    static rg_REAL* chordLength(const rg_INT &n, const rg_Point3D* const ptsPassedThrough);
    static rg_REAL* centripetal(const rg_INT &n, const rg_Point3D* const ptsPassedThrough);

    static rg_REAL computeTriangleArea(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3);
    static rg_REAL computeAreaOfTriangle(rg_Point3D pt1, rg_Point3D pt2 ,rg_Point3D pt3);
    static rg_REAL computeTetrahedronVolume(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4);

    static rg_REAL computeInternalAngle(const rg_Point2D& pt1, const rg_Point2D& pt2, const rg_Point2D& pt3, const rg_INDEX& pntIndex = 1);
    static rg_REAL computeInternalAngle(const rg_REAL& edgeLen1, const rg_REAL& edgeLen2, const rg_REAL& edgeLen3, const rg_INDEX& edgeIndex = 1);
    static rg_REAL computeSignedDihedralAngle(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4);
    static rg_REAL computeAbsoluteDihedralAngle(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4);

    static rg_REAL computeSignedAreaOfTriangle(const rg_Point2D& vtx0, const rg_Point2D& vtx1, const rg_Point2D& vtx2);
    static rg_REAL computeAreaOfTriangle(const rg_REAL& edgeLen1, const rg_REAL& edgeLen2, const rg_REAL& edgeLen3);


    static rg_REAL getDeterminantOf33Matrix( const rg_REAL& mat11, const rg_REAL& mat12, const rg_REAL& mat13,
                                             const rg_REAL& mat21, const rg_REAL& mat22, const rg_REAL& mat23,
                                             const rg_REAL& mat31, const rg_REAL& mat32, const rg_REAL& mat33 );

    static rg_INT  compute_intersection_between_two_circles(const rg_Circle2D& circle1,
                                                            const rg_Circle2D& circle2,
                                                            rg_Point2D* intersectionPt);

    static rg_INT compute_tangent_circles_of_two_circles_with_given_radius( const rg_Circle2D& circle1,
                                                                            const rg_Circle2D& circle2,
                                                                            const rg_INT&      containerCircleIndex,
                                                                            const rg_REAL&     givenRadius,
                                                                            rg_Circle2D* tangentCircle);

    static rg_REAL  computeLengthOfTriangleGivenTwoLengthesNAngle(const rg_REAL& length1, 
															      const rg_REAL& length2, 
															      const rg_REAL& angle);

    // Functions for computing volume, area and intersection among sphere
    static rg_INT computeIntersectionPointAmongThreeSpheres(const Sphere& sphere1,
			                                                const Sphere& sphere2,
												            const Sphere& sphere3,
												            rg_Point3D*&  intersectionPoint);

    static inline rg_REAL computeSignedVolumeOfTetrahedron(const rg_Point3D& point0, 
				                                           const rg_Point3D& point1, 
											               const rg_Point3D& point2, 
											               const rg_Point3D& point3)
    {
	    // Get the coordinates of vertices
	    rg_REAL x0 = point0.getX();
	    rg_REAL y0 = point0.getY();
	    rg_REAL z0 = point0.getZ();

	    rg_REAL x1 = point1.getX();
	    rg_REAL y1 = point1.getY();
	    rg_REAL z1 = point1.getZ();

	    rg_REAL x2 = point2.getX();
	    rg_REAL y2 = point2.getY();
	    rg_REAL z2 = point2.getZ();

	    rg_REAL x3 = point3.getX();
	    rg_REAL y3 = point3.getY();
	    rg_REAL z3 = point3.getZ();

	    rg_REAL signedVolume = 0.0;

    // 	signedVolume = (x1 * y2 * z3 - x1 * y3 * z2 
    // 		          + y1 * z2 * x3 - y1 * x2 * z3 
    // 		          + z1 * x2 * y3 - z1 * y2 * x3 
    // 		          - x0 * y2 * z3 + x0 * y3 * z2 
    // 		          - x0 * y1 * z2 + x0 * y1 * z3 
    // 		          - x0 * z1 * y3 + x0 * z1 * y2 
    // 		          + y0 * x2 * z3 - y0 * z2 * x3 
    // 		          + y0 * x1 * z2 - y0 * x1 * z3 
    // 		          + y0 * z1 * x3 - y0 * z1 * x2 
    // 		          - z0 * x2 * y3 + z0 * y2 * x3 
    // 		          - z0 * x1 * y2 + z0 * x1 * y3 
    // 		          - z0 * y1 * x3 + z0 * y1 * x2) / 6.0;

	
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

    static rg_REAL computeXVolumeOfTwoSpheres(const Sphere& sphere1, const Sphere& sphere2);
    static rg_REAL computeXVolumeOfTwoSpheres(const rg_REAL& radius1,
	                                          const rg_REAL& radius2,
	                                          const rg_REAL& distBTWCenters);
    static rg_REAL computeVolumeOfGoStone(const Sphere& sphere1, const Sphere& sphere2);
    
    static rg_REAL computeXAreaOfTwoSphere(const Sphere& sphere1, const Sphere& sphere2);
    static rg_REAL computeXAreaOfTwoSphere(const rg_REAL& radius1,
                                           const rg_REAL& radius2,
                                           const rg_REAL& distBTWCenters);
    static rg_REAL computeAreaOfGoStone(const Sphere& sphere1, const Sphere& sphere2);

    static rg_REAL computeXCircleAreaOfTwoSphere(const Sphere& sphere1, const Sphere& sphere2);
    static rg_REAL computeXCircleAreaOfTwoSphere(const rg_REAL& radius1,
                                                 const rg_REAL& radius2,
                                                 const rg_REAL& distBTWCenters);

    static rg_REAL computeXCircleRadiusOfTwoSphere(const Sphere& sphere1, const Sphere& sphere2);
    static rg_REAL computeXCircleRadiusOfTwoSphere(const rg_REAL& radius1,
                                                   const rg_REAL& radius2,
                                                   const rg_REAL& distBTWCenters);

    static rg_REAL computeSignedVolumeOfCutAwayTriangularPrism(const Plane& referencePlane, const rg_Point3D& vrtx1, const rg_Point3D& vrtx2, const rg_Point3D& vrtx3);

    
    //static rg_REAL computeSignedVolumeOfCutAwayTriangularPrism(const Plane& referencePlane, const rg_Point3D& vrtx1, const rg_Point3D& vrtx2, const rg_Point3D& vrtx3);


    static rg_Point2D   compute_intersection_between_two_lines(const rg_Point2D& startPointOfLineSegment1, const rg_Point2D& endPointOfLineSegment1, const rg_Point2D& startPointOfLineSegment2, const rg_Point2D& endPointOfLineSegment2, rg_REAL& parameterOfIntersectionForLineSeg1, rg_REAL& parameterOfIntersectionForLineSeg2);
    static rg_Point2D   compute_intersection_between_two_lines(const rg_Point2D& startPointOfLineSegment1, const rg_Point2D& endPointOfLineSegment1, const rg_Point2D& startPointOfLineSegment2, const rg_Point2D& endPointOfLineSegment2, rg_REAL& parameterOfIntersectionForLineSeg1, rg_REAL& parameterOfIntersectionForLineSeg2, bool& bTwoLinesAreParallel);
    static rg_Point2D   compute_intersection_between_two_lines(const rg_Point2D& startPointOfLineSegment1, const rg_Point2D& endPointOfLineSegment1, const rg_Point2D& startPointOfLineSegment2, const rg_Point2D& endPointOfLineSegment2, bool& bTwoLinesAreParallel);
    static bool         compute_intersection_between_two_line_segments(const rg_Line2D& lineSegment1, const rg_Line2D& lineSegment2, rg_Point2D& intersectionPt);
    static bool         compute_intersection_between_two_line_segments(const rg_Point2D& startPointOfLineSegment1, const rg_Point2D& endPointOfLineSegment1, const rg_Point2D& startPointOfLineSegment2, const rg_Point2D& endPointOfLineSegment2, rg_Point2D& intersectionPt);
    static rg_Point2D   compute_intersection_between_two_rays(const rg_Point2D& startPointOfRay1, const rg_Point2D& directionVecOfRay1, const rg_Point2D& startPointOfRay2, const rg_Point2D& directionVecOfRay2, rg_REAL& parameterOfIntersectionForRay1, rg_REAL& parameterOfIntersectionForRay2);
    static rg_Point2D   compute_intersection_between_two_rays(const rg_Point2D& startPointOfRay1, const rg_Point2D& directionVecOfRay1, const rg_Point2D& startPointOfRay2, const rg_Point2D& directionVecOfRay2);

    static rg_INT     compute_intersection_between_circle_and_line_segment(const rg_Circle2D& circle, const rg_Line2D& lineSegment, rg_Point2D intersection[]);

    static bool      compute_intersection_between_circle_and_ellipse(const rg_Circle2D& circle, Ellipse2D& ellipse, rg_Point2D intersection[]);

    static rg_INT     compute_intersection_between_circle_and_point(const rg_Circle2D& circle, const rg_Point2D point);

    static rg_Line2D compute_bisector_line_between_two_line_segments(const rg_Line2D& lineSegment1, const rg_Line2D& lineSegment2);

    static rg_Line2D compute_bisector_line_between_two_points(const rg_Point2D& point1, const rg_Point2D& point2);

static rg_INT compute_intersection_between_parabola_and_line(const rg_Point2D& focusOfParabola, const rg_Line2D&  directrixOfParabola, const rg_Line2D& line, list<rg_Point2D>& intersectionPts );

    static rg_INT compute_intersection_between_ellipse_and_line(const Ellipse2D& ellipse, const rg_Line2D& lineSegment, rg_Point2D intersection[]);

    static rg_Point2D compute_a_point_on_line_whose_distance_with_anchor_point_on_line_is_equal_to_distance_between_the_point_and_its_footprint_of_ellipse(const rg_Point2D seedPoint, const rg_Line2D& line, const rg_Point2D& anchorPointOnLine, const Ellipse2D& ellipse, const rg_REAL& tolerance, const rg_INT& maximumNumberOfIterations);

    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const Ellipse2D& ellipse2, const Ellipse2D& ellipse3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const Ellipse2D& ellipse2, const rg_Line2D& lineSeg3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const Ellipse2D& ellipse2, const rg_Point2D& point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const rg_Line2D& lineSeg2, const rg_Line2D& lineSeg3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const rg_Line2D& lineSeg2, const rg_Point2D& point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const rg_Point2D&  point2, const rg_Point2D& point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const rg_Line2D& lineSeg1, const rg_Line2D& lineSeg2, const rg_Line2D& lineSeg3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const rg_Line2D& lineSeg1, const rg_Line2D& lineSeg2, const rg_Point2D&  point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const rg_Line2D& lineSeg1, const rg_Point2D&  point2, const rg_Point2D&  point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const rg_Point2D&  point1, const rg_Point2D&  point2, const rg_Point2D&  point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);

    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D&  lineSeg3, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const rg_Circle2D& disk1, const rg_Line2D& lineSeg2, const rg_Line2D&  lineSeg3, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);
    static rg_Circle2D compute_empty_tangent_cirlce_of_three_generators(const rg_Circle2D& disk1, const rg_Line2D& lineSeg2, const rg_Point2D& point3, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION = 10e-3, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION = 10);

	static rg_INT compute_tangent_circles_of_point_and_line_with_given_radius(const rg_Point2D & point, const rg_Line2D & line, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle);

	static rg_INT compute_tangent_circles_of_point_and_linesegment_with_given_radius(const rg_Point2D & point, const rg_Line2D & lineSegment, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle);

	// This function locates tangent circles on the same side as a circle center w.r.t. line
	static rg_INT compute_tangent_circles_of_circle_and_line_with_given_radius(const rg_Circle2D& circle, const rg_Line2D& line, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle);

	static rg_INT compute_tangent_circles_of_circle_and_linesegment_with_given_radius(const rg_Circle2D& circle, const rg_Line2D& lineSegment, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle);

	static rg_INT compute_tangent_circles_of_two_lines_not_parallel_to_each_other_with_given_radius(const rg_Line2D& line1, const rg_Line2D& line2, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle);

	static rg_INT compute_tangent_circles_of_two_linesegments_not_parallel_to_each_other_with_given_radius(const rg_Line2D& lineSegment1, const rg_Line2D& lineSegment2, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle);

    // NOTE: The care will be taken for query point on the polygon boundary later on.
    static bool point_in_polygon(const rg_Point2D& queryPoint, const list<rg_Point2D>& verticesOnShell);
};

inline rg_Point2D rg_GeoFunc::angle_bisecting_vector(const rg_Point2D& vector1, const rg_Point2D& vector2, const rg_REAL& tolerance)
{
    rg_Point2D angleBisector = (vector1 + vector2) / 2.0;
    rg_REAL crossProduct = vector1 * vector2;

    // Angle between two vectors is greater than PI
    if (rg_NEG(crossProduct, tolerance))
    {
        angleBisector = -angleBisector;
    }
    // equal to PI
    else if (rg_ZERO(crossProduct, tolerance))
    {
        angleBisector.setPoint(-vector1.getY(), vector1.getX());
    }
    else {}
    
    return angleBisector.getUnitVector();
}

#endif


