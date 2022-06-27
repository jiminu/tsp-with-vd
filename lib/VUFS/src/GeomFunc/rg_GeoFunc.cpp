//#include "defconst.h"
#include "rg_GeoFunc.h"
//#include "rg_Point2D.h"
#include "rg_RelativeOp.h"
#include "rg_MathFunc.h"
//#include <math.h>
#include "Plane.h"
#include "rg_QuadraticPolynomial.h"
#include "rg_ComplexNumber.h"

#include <cmath>
using namespace std;


rg_Line<rg_Point2D> rg_GeoFunc::lineOffset(const rg_Line<rg_Point2D> &line, const rg_REAL &offset)
//rg_Line<rg_Point2D> rg_GeoFunc::lineOffset(rg_Line<rg_Point2D> &line, const rg_REAL &offset)
{   
	rg_REAL px = line.getEP().getX()-line.getSP().getX();
	rg_REAL py = line.getEP().getY()-line.getSP().getY();
	rg_Point2D startPendP(px, py);
	rg_Point2D normal( rg_GeoFunc::getNormal(startPendP) );

	// make orthonormal vector
	rg_REAL root = sqrt(pow(normal.getX(), 2) + pow(normal.getY(), 2));
	rg_Point2D orthonormal( normal.getX()/root, normal.getY()/root );
	                
	px = line.getSP().getX() + offset*orthonormal.getX();
	py = line.getSP().getY() + offset*orthonormal.getY();
	rg_REAL nx = line.getEP().getX() + offset*orthonormal.getX();
	rg_REAL ny = line.getEP().getY() + offset*orthonormal.getY();
		                                            
	return rg_Line<rg_Point2D>(rg_Point2D(px, py), rg_Point2D(nx, ny));	
}

rg_Point2D rg_GeoFunc::getNormal(const rg_Point2D &point)
{               
	rg_REAL px = point.getY();
	rg_REAL py = -point.getX();
                      
	return rg_Point2D(px, py);
}

//		    * pt2
//         /| 
//        / |
//       /  |
//      /   |
// pt3 *----* pt1
// following function returns 
//       unit normal vector for plane contains three points(pt1,pt2,pt3)

rg_Point3D  rg_GeoFunc::getUnitNormal( const rg_Point3D& pt1,
                                       const rg_Point3D& pt2,
							           const rg_Point3D& pt3 )
{
    rg_Point3D vector1=pt2-pt1;
    rg_Point3D vector2=pt3-pt2;

    rg_Point3D normal=vector1*vector2;
    rg_Point3D output=normal.getUnitVector();

    return output;
}


rg_REAL rg_GeoFunc::distanceLineVsPointl(const rg_Line<rg_Point2D> &line, rg_Point2D &point)
{
	// pLine equation : ax + by + c = 0
	rg_REAL a;
	rg_REAL b;
	rg_REAL c;
	
	// define the line equation of one line 
	if( rg_ZERO(line.getEP().getX() - line.getSP().getX()) ) 
	{
		a = 1;
		b = 0;
		c = -line.getSP().getX();
	}
	
	else if( rg_ZERO(line.getEP().getY() - line.getSP().getY()) )
	{
		a = 0;
		b = -1;
		c = line.getSP().getY();
	}
	
	else
	{
		a = (line.getEP().getY()-line.getSP().getY())
		   /(line.getEP().getX()-line.getSP().getX());
		b = -1;                  
		c = -a*line.getSP().getX() + line.getSP().getY();
	}

	rg_REAL distance = fabs(a * point.getX() + b * point.getY() + c)/sqrt(a*a + b*b);

	rg_Point2D pt1( line.getEP() - line.getSP() ); 
	rg_Point2D pt2( point        - line.getSP() );

//	if( crossProduct(pt1, pt2) < 0 ) return distance;
	if( pt1*pt2 < 0 ) return distance;
	else return -distance;

//	return fabs(a * point.getX() + b * point.getY() + c)/sqrt(a*a + b*b);
}

rg_Point3D*        rg_GeoFunc::samplePtsOnArc( const rg_Point3D& center,
                                               const rg_Point3D& localXAxis,
							                   const rg_Point3D& localYAxis,
							                   const rg_REAL&  radius,
							                   const rg_REAL&  startDegree,
							                   const rg_REAL&  endDegree   ,
                                               const rg_INT&   size )
{
    const rg_Point3D& unitLocalXAxis=localXAxis.getUnitVector();
    const rg_Point3D& unitLocalYAxis=localYAxis.getUnitVector();

    const rg_REAL stepRadian=(endDegree-startDegree)/(size-1);
    rg_REAL currentRadian=startDegree;

    rg_Point3D* output= new rg_Point3D[size];

    for( rg_INT i=0; i < size; i++ )
    {
        output[i]=center + radius*cos(currentRadian)*unitLocalXAxis
                         + radius*sin(currentRadian)*unitLocalYAxis;
        currentRadian+=stepRadian;
    }

    return output;
}
        
//                   _* pt3
//                   /|
//                  / 
//                 / 
//    *---------->*  
//   pt1           pt2
rg_REAL rg_GeoFunc::calculateAngle(const rg_Point3D& pt1, 
                                   const rg_Point3D& pt2,
				                   const rg_Point3D& pt3 )
{
	rg_Point3D vet1 = getUnitVector(pt2-pt1);
	rg_Point3D vet2 = getUnitVector(pt3-pt2);

	rg_REAL angle = 0.0;
    rg_REAL  pi=4.0*atan(1.0);
	angle = vet1.innerProduct(vet2);
	if ( rg_EQ(angle, 1.0) )
	{
		angle = 0.0;
	}
	else if ( rg_EQ(angle, -1.0) )
	{
		angle = -pi;
	}
	else
	{
		angle = acos(angle);
	}

	return angle;
}

//       _* vect2
//       /|
//      / 
//     / 
//    *---------->*  
//               vect1
rg_REAL rg_GeoFunc::calculateAngle(const rg_Point3D& vect1, 
                                   const rg_Point3D& vect2 )
{

	rg_REAL angle = 0.0;
	//rg_REAL	RADIAN = 57.2957795130823;
    rg_REAL  pi=4.0*atan(1.0);
	rg_Point3D vector1=getUnitVector(vect1);
	rg_Point3D vector2=getUnitVector(vect2);

	angle = vector1.innerProduct(vector2);
	if ( rg_EQ(angle, 1.0) )
	{
		angle = 0.0;
	}
	else if ( rg_EQ(angle, -1.0) )
	{
		angle = -pi;
	}
	else
	{
		angle = acos(angle);
	}

	return angle;
}


rg_REAL rg_GeoFunc::calculateAngle(const rg_Point2D& vect1,
    const rg_Point2D& vect2)
{
    rg_REAL angle = 0.0;
    //rg_REAL	RADIAN = 57.2957795130823;
    rg_REAL  pi = 4.0*atan(1.0);
    rg_Point2D vector1 = vect1.getUnitVector();
    rg_Point2D vector2 = vect2.getUnitVector();

    angle = vector1.operator%(vector2);
    if (rg_EQ(angle, 1.0))
    {
        angle = 0.0;
    }
    else if (rg_EQ(angle, -1.0))
    {
        angle = -pi;
    }
    else
    {
        angle = acos(angle);
    }

    return angle;
}


rg_REAL rg_GeoFunc::getNormalizedAngle(const rg_REAL &angle)
{
    const rg_REAL doublePi=2*4.0*atan(1.0);
    rg_REAL step=floor(angle/(doublePi))*doublePi;
    rg_REAL output=angle-step;

    return output;
}

void rg_GeoFunc::convertNormalizedAngle(rg_REAL& start,
                                        rg_REAL& end )
{
    const rg_REAL doublePi=2*4.0*atan(1.0);
    start=rg_GeoFunc::getNormalizedAngle(start);
    end=rg_GeoFunc::getNormalizedAngle(end);

    if ( rg_GE(start, end) )
    {
        end=end+doublePi;
    }
}

rg_REAL* rg_GeoFunc::chordLength(const rg_INT &n, const rg_Point3D* const ptsPassedThrough)
{
	rg_REAL d = 0.0;
	rg_INT i = 0;
	for(i=1; i<n; i++) 
        d += ptsPassedThrough[i].distance(ptsPassedThrough[i-1]);

	///////////////////////////
	rg_REAL *u = new rg_REAL[n];

	u[0] = 0;
	u[n-1] = 1;

	for(i=1; i<n-1; i++)
        u[i] = u[i-1] 
               + ptsPassedThrough[i].distance(ptsPassedThrough[i-1])/d;

	return u;
}

rg_REAL* rg_GeoFunc::centripetal(const rg_INT &n, const rg_Point3D* const ptsPassedThrough)
{
	rg_REAL d = 0;
	rg_INT i = 0;
	for (i=1; i<n; i++) 
        d += sqrt(ptsPassedThrough[i].distance(ptsPassedThrough[i-1]));

	///////////////////////////
	rg_REAL *u = new rg_REAL[n];

	u[0] = 0;
	u[n-1] = 1;

	for (i=1; i<n-1; i++)
        u[i] = u[i-1] 
               + sqrt(ptsPassedThrough[i].distance(ptsPassedThrough[i-1]))/d;

	return u;
}


rg_REAL rg_GeoFunc::computeTriangleArea(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3)
{
    rg_REAL edgeLength[3];
    edgeLength[0] = pt1.distance( pt2 );
    edgeLength[1] = pt2.distance( pt3 );
    edgeLength[2] = pt3.distance( pt1 );

    rg_REAL halfPerimeter = (edgeLength[0]+edgeLength[1]+edgeLength[2])*0.5;

    rg_REAL faceArea = sqrt( halfPerimeter*(halfPerimeter-edgeLength[0])*(halfPerimeter-edgeLength[1])*(halfPerimeter-edgeLength[2]) );

    return faceArea;
}

rg_REAL rg_GeoFunc::computeAreaOfTriangle(rg_Point3D pt1, rg_Point3D pt2 ,rg_Point3D pt3)
{
    rg_Point3D vec1, vec2;
    rg_REAL triArea = 0.0 ;

    vec1 = pt2 - pt1;
    vec2 = pt3 - pt2;
    triArea =  vec1.crossProduct(vec2).magnitude() / 2.0; 

    return triArea;
}

rg_REAL rg_GeoFunc::computeTetrahedronVolume(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4)
{
    rg_REAL area = computeTriangleArea(pt1, pt2, pt3);

    Plane plane;
    plane.definePlaneByThreePoints(pt1, pt2, pt3);
    rg_REAL height = plane.distanceFromPoint( pt4 );
    height = rg_ABS( height );

    rg_REAL tetrahedronVolume = area*height/3.0;

    return tetrahedronVolume;
}


rg_REAL rg_GeoFunc::computeInternalAngle(const rg_Point2D& pt1, const rg_Point2D& pt2, const rg_Point2D& pt3, const rg_INDEX& pntIndex)
{
	rg_REAL len12 = pt1.distance(pt2);
	rg_REAL len23 = pt2.distance(pt3);
	rg_REAL len31 = pt3.distance(pt1);

	rg_REAL angle = -1.;

	switch(pntIndex)
	{
	case 1:
		angle = acos((len31 * len31 + len12 * len12 - len23 * len23) / (2. * len31 * len12));
		break;
	case 2:
		angle = acos((len12 * len12 + len23 * len23 - len31 * len31) / (2. * len12 * len23));
		break;
	case 3:
		angle = acos((len23 * len23 + len31 * len31 - len12 * len12) / (2. * len23 * len31));
		break;
	default:
		break;
	}

	return angle;
}


rg_REAL rg_GeoFunc::computeInternalAngle(const rg_REAL& edgeLen1, const rg_REAL& edgeLen2, const rg_REAL& edgeLen3, const rg_INDEX& edgeIndex)
{
	rg_REAL angle = -1.;

	switch(edgeIndex)
	{
	case 1:
		angle = acos((edgeLen2 * edgeLen2 + edgeLen3 * edgeLen3 - edgeLen1 * edgeLen1) / (2. * edgeLen2 * edgeLen3));		
		break;
	case 2:
		angle = acos((edgeLen3 * edgeLen3 + edgeLen1 * edgeLen1 - edgeLen2 * edgeLen2) / (2. * edgeLen3 * edgeLen1));		
		break;
	case 3:
		angle = acos((edgeLen1 * edgeLen1 + edgeLen2 * edgeLen2 - edgeLen3 * edgeLen3) / (2. * edgeLen1 * edgeLen2));		
		break;
	default:
		break;
	}

	return angle;
}

//Refer to the following article for the definition of dihedral angle
//IUPAC-IUB Commission on Biochemical Nomenclature. 
//Abbreviations and Symbols for the Description of the Conformation of Polypeptide Chains. 
//Tentative Rules(1969), Biochemistry Journal, 1971, 121, 577-585.
//range of angle (-180, +180]

rg_REAL rg_GeoFunc::computeSignedDihedralAngle(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4)
{
	rg_Point3D normal = pt2 - pt3;
	normal.normalize();
	rg_Point3D passPt = (pt2 + pt3) / 2.0;
	Plane plane(normal, passPt);

	rg_Point3D projectionOfPt1 = plane.projectPointOnPlane(pt1);
	rg_Point3D projectionOfPt4 = plane.projectPointOnPlane(pt4);
	rg_Point3D projectionOfPt2 = passPt;

	rg_REAL angle = plane.computeAngleInCCW(projectionOfPt4, projectionOfPt2, projectionOfPt1);

	if(rg_GT(angle, rg_PI))
	{
		angle = - (2.0 * rg_PI - angle);
	}
	return angle;
}

rg_REAL rg_GeoFunc::computeAbsoluteDihedralAngle(const rg_Point3D& pt1, const rg_Point3D& pt2, const rg_Point3D& pt3, const rg_Point3D& pt4)
{
    rg_Point3D normal = pt2 - pt3;
    normal.normalize();
    rg_Point3D passPt = (pt2 + pt3) / 2.0;
    Plane plane(normal, passPt);

    rg_Point3D projectionOfPt1 = plane.projectPointOnPlane(pt1);
    rg_Point3D projectionOfPt4 = plane.projectPointOnPlane(pt4);
    rg_Point3D projectionOfPt2 = passPt;

    rg_REAL angle = plane.computeAngleInCCW(projectionOfPt4, projectionOfPt2, projectionOfPt1);

    return angle;
}

rg_REAL rg_GeoFunc::computeSignedAreaOfTriangle(const rg_Point2D& vtx0, 
                                                const rg_Point2D& vtx1, 
						                        const rg_Point2D& vtx2)
{
	rg_REAL x0 = vtx0.getX();
	rg_REAL y0 = vtx0.getY();
	rg_REAL x1 = vtx1.getX();
	rg_REAL y1 = vtx1.getY();
	rg_REAL x2 = vtx2.getX();
	rg_REAL y2 = vtx2.getY();		

	rg_REAL determinant = (x1 * y2 - x2 * y1) 
		                - (x0 * y2 - x2 * y0) 
						+ (x0 * y1 - x1 * y0);

	rg_REAL signedArea = determinant * 0.5;
	
	return signedArea;
}

rg_REAL rg_GeoFunc::computeAreaOfTriangle(const rg_REAL& edgeLen1, const rg_REAL& edgeLen2, const rg_REAL& edgeLen3)
{
	rg_REAL s = 0.5 * (edgeLen1 + edgeLen2 + edgeLen3);
	rg_REAL area = 0.0;
	rg_REAL sqdArea =  s * (s - edgeLen1) * (s - edgeLen2) * (s - edgeLen3);
	if(sqdArea > 0.)
		area = sqrt(sqdArea);
	return area;
}


rg_REAL rg_GeoFunc::getDeterminantOf33Matrix( const rg_REAL& mat11, const rg_REAL& mat12, const rg_REAL& mat13,
                                              const rg_REAL& mat21, const rg_REAL& mat22, const rg_REAL& mat23,
                                              const rg_REAL& mat31, const rg_REAL& mat32, const rg_REAL& mat33 )
{
	rg_REAL det_Value =   mat11 * ( ( mat22 * mat33 ) - ( mat23 * mat32 ) )
	 				    - mat12 * ( ( mat21 * mat33 ) - ( mat23 * mat31 ) ) 
				        + mat13 * ( ( mat21 * mat32 ) - ( mat22 * mat31 ) );
	return det_Value;
}


rg_REAL rg_GeoFunc::computeLengthOfTriangleGivenTwoLengthesNAngle(const rg_REAL& length1, const rg_REAL& length2, const rg_REAL& angle)
{
	rg_REAL len = sqrt(length1 * length1 + length2 * length2 - 2. * length1 * length2 * cos(angle));
	return len;
}

// Please consult a technical report ... for the explanation of this function

rg_INT rg_GeoFunc::computeIntersectionPointAmongThreeSpheres(const Sphere& sphere1, 
															 const Sphere& sphere2, 
															 const Sphere& sphere3, 
															 rg_Point3D*& intersectionPoint)
{
	// Get the coordinates of center points and the radii
	rg_Point3D center1 = sphere1.getCenter();
	rg_Point3D center2 = sphere2.getCenter();
	rg_Point3D center3 = sphere3.getCenter();

	rg_REAL    x1      = center1.getX();
	rg_REAL    y1      = center1.getY();
	rg_REAL    z1      = center1.getZ();

	rg_REAL    x2      = center2.getX();
	rg_REAL    y2      = center2.getY();
	rg_REAL    z2      = center2.getZ();

	rg_REAL    x3      = center3.getX();
	rg_REAL    y3      = center3.getY();
	rg_REAL    z3      = center3.getZ();

	// For handling the offset volume and area, the radii are inflated by "m_betaValue".
	rg_REAL    r1      = sphere1.getRadius();
	rg_REAL    r2      = sphere2.getRadius();
	rg_REAL    r3      = sphere3.getRadius();

	// 1. compute normal vectors of planes P1 and P2 which contain intersection circles	
	// between sphere1 and sphere3
	rg_Point3D normalVec13(-2*x1+2*x3, -2*y1+2*y3, -2*z1+2*z3);
	// between sphere1 and sphere2
	rg_Point3D normalVec12(-2*x1+2*x2, -2*y1+2*y2, -2*z1+2*z2);

	// 2. compute constant terms of planes
	rg_REAL    constant13 = -x1*x1-y1*y1-z1*z1+r1*r1 + x3*x3+y3*y3+z3*z3-r3*r3;
	rg_REAL    constant12 = -x1*x1-y1*y1-z1*z1+r1*r1 + x2*x2+y2*y2+z2*z2-r2*r2;

	// 3. compute direction vector of the intersection line L12 between P1 and P2
	rg_Point3D dirVec = normalVec13.crossProduct(normalVec12);

	// 4. compute passing point of the intersection line L12
	rg_REAL innerProd    = normalVec13.innerProduct(normalVec12);
	rg_REAL sqdNorm13    = normalVec13.innerProduct(normalVec13);
	rg_REAL sqdNorm12    = normalVec12.innerProduct(normalVec12);
	rg_REAL denom = innerProd * innerProd - sqdNorm13 * sqdNorm12;
	rg_REAL coeff13 = (constant12 * innerProd - constant13 * sqdNorm12) / denom;
	rg_REAL coeff12 = (constant13 * innerProd - constant12 * sqdNorm13) / denom;

	rg_Point3D passingPt = coeff13 * normalVec13 + coeff12 * normalVec12;

	// 5. compute intersection among three sphere by solving the quadratic equation
	// which is made of L12 and sphere1
	
	// passing point
	rg_REAL px = passingPt.getX();
	rg_REAL py = passingPt.getY();
	rg_REAL pz = passingPt.getZ();

	// direction vector
	rg_REAL vx = dirVec.getX();
	rg_REAL vy = dirVec.getY();
	rg_REAL vz = dirVec.getZ();

	// solve quadratic equation
	rg_REAL a2 = vx * vx + vz * vz + vy * vy;
	rg_REAL a1 = 2 * (py - y1) * vy + 2 * (px - x1) * vx + 2 * (pz - z1) * vz;
	rg_REAL a0 = (px-x1)*(px-x1)-r1*r1+(pz-z1)*(pz-z1)+(py-y1)*(py-y1);

	rg_REAL param0 = 0.;
	rg_REAL param1 = 0.;

	if(a2 <= 1.0e-6)
	{
		param0 = param1 = - a0 / a1;
	}
	else
	{
		rg_REAL determinant = a1 * a1 - 4 * a2 * a0;
		rg_REAL sRootDeterminant = 0.0;
		if(determinant > 0.0)
			sRootDeterminant = sqrt(determinant);
		param0 = (-a1 + sRootDeterminant) / (2 * a2);
		param1 = (-a1 - sRootDeterminant) / (2 * a2);
	}

	intersectionPoint[0] = passingPt + param0 * dirVec;
	intersectionPoint[1] = passingPt + param1 * dirVec;

	rg_INT numIntersections = 2;

	return numIntersections;
}

// inline rg_REAL rg_GeoFunc::computeSignedVolumeOfTetrahedron(const rg_Point3D& point0, 
// 														    const rg_Point3D& point1, 
// 															const rg_Point3D& point2, 
// 															const rg_Point3D& point3)
// {
// 	// Get the coordinates of vertices
// 	rg_REAL x0 = point0.getX();
// 	rg_REAL y0 = point0.getY();
// 	rg_REAL z0 = point0.getZ();
// 
// 	rg_REAL x1 = point1.getX();
// 	rg_REAL y1 = point1.getY();
// 	rg_REAL z1 = point1.getZ();
// 
// 	rg_REAL x2 = point2.getX();
// 	rg_REAL y2 = point2.getY();
// 	rg_REAL z2 = point2.getZ();
// 
// 	rg_REAL x3 = point3.getX();
// 	rg_REAL y3 = point3.getY();
// 	rg_REAL z3 = point3.getZ();
// 
// 	rg_REAL signedVolume = 0.0;
// 
// // 	signedVolume = (x1 * y2 * z3 - x1 * y3 * z2 
// // 		          + y1 * z2 * x3 - y1 * x2 * z3 
// // 		          + z1 * x2 * y3 - z1 * y2 * x3 
// // 		          - x0 * y2 * z3 + x0 * y3 * z2 
// // 		          - x0 * y1 * z2 + x0 * y1 * z3 
// // 		          - x0 * z1 * y3 + x0 * z1 * y2 
// // 		          + y0 * x2 * z3 - y0 * z2 * x3 
// // 		          + y0 * x1 * z2 - y0 * x1 * z3 
// // 		          + y0 * z1 * x3 - y0 * z1 * x2 
// // 		          - z0 * x2 * y3 + z0 * y2 * x3 
// // 		          - z0 * x1 * y2 + z0 * x1 * y3 
// // 		          - z0 * y1 * x3 + z0 * y1 * x2) / 6.0;
// 
// 	
// 	// Optimized by Maple 12
// 	// Before optimization: 23 additions 48 multiplications
// 	// Afeter optimization: 17 additions 16 multiplications 7 assignments
// 
// 	rg_REAL t1, t2, t3, t4, t5, t6;
// 	t6 =  z0-z2;
// 	t5 =  z1-z0; 
// 	t4 =  z1-z3; 
// 	t3 =  z2-z1; 
// 	t2 =  z2-z3; 
// 	t1 = -z3+z0; 
// 	signedVolume 
// 		= ((-t5*y2-t6*y1-t3*y0)*x3+(t5*y3+t1*y1-t4*y0)*x2+(t6*y3-t1*y2+t2*y0)*x1+(t3*y3+t4*y2-t2*y1)*x0)/6.0;
// 
// 	return signedVolume;
// }


rg_REAL rg_GeoFunc::computeXVolumeOfTwoSpheres(const Sphere& sphere1, const Sphere& sphere2)
{
    rg_Point3D center1 = sphere1.getCenter();
    rg_Point3D center2 = sphere2.getCenter();

    rg_REAL    distBTWCenters  = center1.distance(center2);
    rg_REAL    radius1 = sphere1.getRadius();
    rg_REAL    radius2 = sphere2.getRadius();

    rg_REAL    XVolume = 0.0;
    XVolume = computeXVolumeOfTwoSpheres(radius1, radius2, distBTWCenters);
    return XVolume;
}

rg_REAL rg_GeoFunc::computeXVolumeOfTwoSpheres(const rg_REAL& radius1,
	                                           const rg_REAL& radius2,
	                                           const rg_REAL& distBTWCenters)
{
	rg_REAL distDiff = fabs(radius1 - radius2);
	rg_REAL sumOfTwoRadii  = radius1 + radius2;
	rg_REAL XVolume = 0.0;

	if(distBTWCenters <= distDiff)
	{
		rg_REAL minRadius = (radius1 < radius2) ? radius1 : radius2;
		XVolume = 4.0 * rg_PI * minRadius * minRadius * minRadius / 3.0;
	}
	else if( (distDiff < distBTWCenters) && (distBTWCenters < sumOfTwoRadii) )
	{
		rg_REAL    radius1_squared = radius1 * radius1;
		rg_REAL    radius2_squared = radius2 * radius2;

		// Formula by Scheraga
		XVolume
			= rg_PI* ( 2.0/3.0*( (radius1_squared*radius1)+(radius2_squared*radius2)+(distBTWCenters*distBTWCenters*distBTWCenters*0.125) ) 
			- 0.25*( 2.0*distBTWCenters*(radius1_squared + radius2_squared) + (radius1_squared - radius2_squared)*(radius1_squared - radius2_squared)/distBTWCenters ) );

		//rg_REAL radius1_squared = radius1 * radius1;
		//rg_REAL radius2_squared = radius2 * radius2;
		//rg_REAL radius1_cubic   = radius1_squared * radius1;
		//rg_REAL radius2_cubic   = radius2_squared * radius2;

		//XVolume += ( 2.0 * rg_PI * (radius1_cubic + radius2_cubic + distBTWCenters * distBTWCenters * distBTWCenters / 8.0 ) / 3.0 );
		//XVolume -= ( rg_PI * (radius1_squared + radius2_squared) * distBTWCenters / 2.0 );
		//XVolume -= ( rg_PI * (radius1_squared - radius2_squared) * (radius1_squared - radius2_squared) / (4.0 * distBTWCenters) );
	}
	else
	{
		XVolume = 0.0;
	}

	return XVolume;
}

rg_REAL rg_GeoFunc::computeVolumeOfGoStone(const Sphere& sphere1, const Sphere& sphere2)
{
	// Get the coordinates of center points and the radii
	rg_Point3D center1 = sphere1.getCenter();
	rg_Point3D center2 = sphere2.getCenter();
	
	rg_REAL    distBTWCenters  = center1.distance(center2);
	rg_REAL    r1      = sphere1.getRadius();
	rg_REAL    r2      = sphere2.getRadius();
	rg_REAL    r1_squared = r1 * r1;
	rg_REAL    r2_squared = r2 * r2;

	// Formula by Scheraga
	rg_REAL volumeOfGoStone
	= rg_PI* ( 2.0/3.0*( (r1_squared*r1)+(r2_squared*r2)+(distBTWCenters*distBTWCenters*distBTWCenters*0.125) ) 
             - 0.25*( 2.0*distBTWCenters*(r1_squared + r2_squared) + (r1_squared - r2_squared)*(r1_squared - r2_squared)/distBTWCenters ) );

	return volumeOfGoStone;
}

rg_REAL rg_GeoFunc::computeXAreaOfTwoSphere(const Sphere& sphere1, const Sphere& sphere2)
{
    rg_Point3D center1 = sphere1.getCenter();
    rg_Point3D center2 = sphere2.getCenter();

    rg_REAL    distBTWCenters  = center1.distance(center2);
    rg_REAL    radius1 = sphere1.getRadius();
    rg_REAL    radius2 = sphere2.getRadius();

    rg_REAL    bndryAreaOfXVol = 0.0;
    bndryAreaOfXVol = computeXAreaOfTwoSphere(radius1, radius2, distBTWCenters);
    return bndryAreaOfXVol;
}

rg_REAL rg_GeoFunc::computeXAreaOfTwoSphere(const rg_REAL& radius1,
                                            const rg_REAL& radius2,
                                            const rg_REAL& distBTWCenters)
{
    rg_REAL distDiff = fabs(radius1 - radius2);
    rg_REAL sumOfTwoRadii  = radius1 + radius2;
    rg_REAL bndryAreaOfXVol = 0.0;

    if(distBTWCenters <= distDiff)
    {
        rg_REAL minRadius = (radius1 < radius2) ? radius1 : radius2;
        bndryAreaOfXVol = 4.0 * rg_PI * minRadius * minRadius;
    }
    else if( (distDiff < distBTWCenters) && (distBTWCenters < sumOfTwoRadii) )
    {
        rg_REAL    radius1_squared = radius1 * radius1;
        rg_REAL    radius2_squared = radius2 * radius2;

        // Formula by Scheraga
        bndryAreaOfXVol
            =  rg_PI * ( 2.0*( radius1_squared + radius2_squared ) - distBTWCenters*(radius1+radius2) - (radius1-radius2)*(radius1_squared-radius2_squared)/distBTWCenters );
    }
    else
    {
        bndryAreaOfXVol = 0.0;
    }

    return bndryAreaOfXVol;
}

rg_REAL rg_GeoFunc::computeAreaOfGoStone(const Sphere& sphere1, const Sphere& sphere2)
{
	// Get the coordinates of center points and the radii
	rg_Point3D center1 = sphere1.getCenter();
	rg_Point3D center2 = sphere2.getCenter();
	rg_REAL distBTWCenters = center1.distance(center2);

// 	rg_REAL    x1      = center1.getX();
// 	rg_REAL    y1      = center1.getY();
// 	rg_REAL    z1      = center1.getZ();
// 
// 	rg_REAL    x2      = center2.getX();
// 	rg_REAL    y2      = center2.getY();
// 	rg_REAL    z2      = center2.getZ();

	rg_REAL    r1      = sphere1.getRadius();
	rg_REAL    r2      = sphere2.getRadius();
	rg_REAL    r1_squared = r1 * r1;
	rg_REAL    r2_squared = r2 * r2;

	// Get the plane which contains intersection circle
	// between sphere1 and sphere2
// 	rg_REAL    constant12 = x1*x1+y1*y1+z1*z1-r1*r1 -x2*x2-y2*y2-z2*z2+r2*r2;
// 	rg_REAL    n12x = -2*x1+2*x2;
// 	rg_REAL    n12y = -2*y1+2*y2;
// 	rg_REAL    n12z = -2*z1+2*z2;
// 	rg_REAL    denom = sqrt(n12x * n12x + n12y * n12y + n12z * n12z);
// 
// 	rg_REAL    constant21 = -x1*x1-y1*y1-z1*z1+r1*r1 +x2*x2+y2*y2+z2*z2-r2*r2;
// 	rg_REAL    n21x = 2*x1-2*x2;
// 	rg_REAL    n21y = 2*y1-2*y2;
// 	rg_REAL    n21z = 2*z1-2*z2;
// 
// 	// Compute the heights of two spherical caps
// 	rg_REAL    signedDist12 = (n12x * x1 + n12y * y1 + n12z * z1 + constant12) / denom;
// 	rg_REAL    signedDist21 = (n21x * x2 + n21y * y2 + n21z * z2 + constant21) / denom;
// 
// 	// Compute the areas of two spherical caps
// 	rg_REAL areaOfSphericalCap12 = 0.;
// 	rg_REAL areaOfSphericalCap21 = 0.;
// 
// 	if(signedDist12 < 0)
// 	{
// 		areaOfSphericalCap12 
// 			= TWO_PI * r1 * (r1 + signedDist12);
// 	}
// 	else
// 	{
// 		areaOfSphericalCap12 
// 			= FOUR_PI * r1 * r1
// 			- TWO_PI * r1 * (r1 - signedDist12);
// 	}
// 
// 	if(signedDist21 < 0)
// 	{
// 		areaOfSphericalCap21 
// 			= TWO_PI * r2 * (r2 + signedDist21);
// 	}
// 	else
// 	{
// 		areaOfSphericalCap21 
// 			= FOUR_PI * r2 * r2
// 			- TWO_PI * r2 * (r2 - signedDist21);
// 	}
// 
// 	rg_REAL areaOfGoStone 
// 		= areaOfSphericalCap12 + areaOfSphericalCap21;
// 
// 	return areaOfGoStone;		

	// Formula by Scheraga	
	rg_REAL areaOfGoStone 
		=  rg_PI * ( 2.0*( r1_squared + r2_squared ) - distBTWCenters*(r1+r2) - (r1-r2)*(r1_squared-r2_squared)/distBTWCenters ); 
	return areaOfGoStone;
}


rg_REAL rg_GeoFunc::computeXCircleAreaOfTwoSphere(const Sphere& sphere1, const Sphere& sphere2)
{
    rg_Point3D center1 = sphere1.getCenter();
    rg_Point3D center2 = sphere2.getCenter();

    rg_REAL    distBTWCenters  = center1.distance(center2);
    rg_REAL    radius1 = sphere1.getRadius();
    rg_REAL    radius2 = sphere2.getRadius();

    rg_REAL    xCircleArea = 0.0;
    xCircleArea = computeXCircleAreaOfTwoSphere(radius1, radius2, distBTWCenters);
    return xCircleArea;
}

rg_REAL rg_GeoFunc::computeXCircleAreaOfTwoSphere(const rg_REAL& radius1,
                                                  const rg_REAL& radius2,
                                                  const rg_REAL& distBTWCenters)
{
    rg_REAL xCircleRadius = computeXCircleRadiusOfTwoSphere(radius1, radius2, distBTWCenters);
    
    rg_REAL xCircleArea = 0.0;
    xCircleArea = xCircleRadius * xCircleRadius * rg_PI;    

    return xCircleArea;
}


rg_REAL rg_GeoFunc::computeXCircleRadiusOfTwoSphere(const Sphere& sphere1, const Sphere& sphere2)
{
    rg_Point3D center1 = sphere1.getCenter();
    rg_Point3D center2 = sphere2.getCenter();

    rg_REAL    distBTWCenters  = center1.distance(center2);
    rg_REAL    radius1 = sphere1.getRadius();
    rg_REAL    radius2 = sphere2.getRadius();

    rg_REAL xCircleRadius = 0.0;
    xCircleRadius = computeXCircleRadiusOfTwoSphere(radius1, radius2, distBTWCenters);

    return xCircleRadius;
}

rg_REAL rg_GeoFunc::computeXCircleRadiusOfTwoSphere(const rg_REAL& radius1,
                                                    const rg_REAL& radius2,
                                                    const rg_REAL& distBTWCenters)
{
    rg_REAL distDiff = fabs(radius1 - radius2);
    rg_REAL sumOfTwoRadii  = radius1 + radius2;
    rg_REAL xCircleRadius = 0.0;

    if(distBTWCenters <= distDiff)
    {
        xCircleRadius = 0.0;
    }
    else if( (distDiff < distBTWCenters) && (distBTWCenters < sumOfTwoRadii) )
    {
        rg_REAL    radius1_squared = radius1 * radius1;
        rg_REAL    radius2_squared = radius2 * radius2;
        rg_REAL    distBTWCenters_squared = distBTWCenters * distBTWCenters;
        rg_REAL    temp = (radius1_squared - radius2_squared + distBTWCenters_squared);
        rg_REAL    numerator = temp * temp;

        rg_REAL    tSquaredRadius = radius1_squared - numerator / ( 4.0 * distBTWCenters_squared );

        if(rg_POS(tSquaredRadius))
            xCircleRadius = sqrt( tSquaredRadius );
        else
            xCircleRadius = 0.0;
    }
    else
    {
        xCircleRadius = 0.0;
    }

    return xCircleRadius;
}


rg_REAL rg_GeoFunc::computeSignedVolumeOfCutAwayTriangularPrism(const Plane& referencePlane, const rg_Point3D& vrtx1, const rg_Point3D& vrtx2, const rg_Point3D& vrtx3)
{
	rg_REAL volume = 0.;
	// check the sign of volume
	rg_Point3D vec1 = vrtx3 - vrtx2;
	rg_Point3D vec2 = vrtx1 - vrtx2;
	rg_Point3D normalVecOfTri = vec1.crossProduct(vec2);
	normalVecOfTri.normalize();
	rg_REAL sign = 0.;
	rg_Point3D normalVecOfPlane = referencePlane.getNormal();
	if(normalVecOfPlane.innerProduct(normalVecOfTri) >= 0)
		sign = 1.;
	else
		sign = -1.;

	// compute the volume of triangular prism
	rg_REAL dist[ 3 ];
	rg_Point3D vrtx[ 3 ], prjVrtx[ 3 ];
	vrtx[ 0 ] = vrtx1;
	vrtx[ 1 ] = vrtx2;
	vrtx[ 2 ] = vrtx3;

    rg_INT i = 0;
	for(i = 0;i < 3;i++)
	{
		dist[ i ] = rg_ABS(referencePlane.distanceFromPoint(vrtx[ i ]));
	}
	rg_REAL maxDist = dist[ 0 ];
	rg_INT  maxIndex = 0;
	for(i = 1;i < 3;i++)
	{
		if(maxDist < dist[ i ])
		{
			maxDist = dist[ i ];
			maxIndex = i;
		}
	}
	rg_INT index1, index2;
	if(maxIndex == 0)
	{
		index1 = 1;
		index2 = 2;
	}
	else if(maxIndex == 1)
	{
		index1 = 0;
		index2 = 2;
	}
	else
	{
		index1 = 0;
		index2 = 1;
	}
	for(i = 0;i < 3;i++) {
		//prjVrtx[ i ] = vrtx[ i ] + normalVecOfPlane * (dist[maxIndex] - dist[ i ]);
        prjVrtx[ i ] = referencePlane.projectPointOnPlane( vrtx[i] );
    }

	rg_REAL heightOfPrism  = dist[maxIndex];
	rg_REAL areaOfTriangle = rg_GeoFunc::computeTriangleArea(prjVrtx[ 0 ], prjVrtx[ 1 ], prjVrtx[ 2 ]);
	rg_REAL volumeOfPrism  = heightOfPrism * areaOfTriangle;
	if( rg_NPOS( volumeOfPrism ) )
		return 0.;

	// compute the volume of pyramid
	if( rg_EQ(dist[maxIndex], dist[index1]) && rg_EQ(dist[maxIndex], dist[index2]) ) {
		volume = volumeOfPrism;
    }
    else {
	    rg_REAL heightOfTrapezoid = prjVrtx[index1].distance(prjVrtx[index2]);
	    if( rg_ZERO(heightOfTrapezoid) ) {
		    volume = volumeOfPrism;
        }
        else {
	        rg_REAL areaOfTrapezoid = heightOfTrapezoid * (2. * dist[maxIndex] - dist[index1] - dist[index2]) / 2.0;
	        rg_REAL heightOfPyramid = 2. * areaOfTriangle / heightOfTrapezoid;
	        rg_REAL volumeOfPyramid = heightOfPyramid * areaOfTrapezoid / 3.0;

	        volume = volumeOfPrism - volumeOfPyramid;
        }
    }
	    
    return (sign * volume);
} 

rg_Point2D rg_GeoFunc::compute_intersection_between_two_lines(const rg_Point2D& startPointOfLineSegment1, const rg_Point2D& endPointOfLineSegment1,
    const rg_Point2D& startPointOfLineSegment2, const rg_Point2D& endPointOfLineSegment2,
    rg_REAL& parameterOfIntersectionForLineSeg1, rg_REAL& parameterOfIntersectionForLineSeg2)
{
    rg_REAL startX1 = startPointOfLineSegment1.getX();
    rg_REAL startY1 = startPointOfLineSegment1.getY();
    rg_REAL endX1 = endPointOfLineSegment1.getX();
    rg_REAL endY1 = endPointOfLineSegment1.getY();

    rg_REAL startX2 = startPointOfLineSegment2.getX();
    rg_REAL startY2 = startPointOfLineSegment2.getY();
    rg_REAL endX2 = endPointOfLineSegment2.getX();
    rg_REAL endY2 = endPointOfLineSegment2.getY();

    parameterOfIntersectionForLineSeg1 = (endX2 * startY1 - endX2 * startY2 - endY2 * startX1 + endY2 * startX2 + startX1 * startY2 - startX2 * startY1) / (endX1 * endY2 - endX1 * startY2 - endX2 * endY1 + endX2 * startY1 + endY1 * startX2 - endY2 * startX1 + startX1 * startY2 - startX2 * startY1);
    parameterOfIntersectionForLineSeg2 = (endX1 * startY1 - endX1 * startY2 - endY1 * startX1 + endY1 * startX2 + startX1 * startY2 - startX2 * startY1) / (endX1 * endY2 - endX1 * startY2 - endX2 * endY1 + endX2 * startY1 + endY1 * startX2 - endY2 * startX1 + startX1 * startY2 - startX2 * startY1);

    return ((1. - parameterOfIntersectionForLineSeg1)*startPointOfLineSegment1 + parameterOfIntersectionForLineSeg1*endPointOfLineSegment1);
}



rg_Point2D rg_GeoFunc::compute_intersection_between_two_lines(const rg_Point2D& startPointOfLineSegment1, const rg_Point2D& endPointOfLineSegment1,
    const rg_Point2D& startPointOfLineSegment2, const rg_Point2D& endPointOfLineSegment2,
    rg_REAL& parameterOfIntersectionForLineSeg1, rg_REAL& parameterOfIntersectionForLineSeg2,
    bool& bTwoLinesAreParallel)
{
    rg_Point2D lineSeg1Vec = endPointOfLineSegment1 - startPointOfLineSegment1;
    rg_Point2D lineSeg2Vec = endPointOfLineSegment2 - startPointOfLineSegment2;

    if (rg_ZERO(lineSeg1Vec.operator*(lineSeg2Vec)))
    {
        bTwoLinesAreParallel = true;
        return rg_Point2D(DBL_MAX, DBL_MAX);
    }
    else
    {
        bTwoLinesAreParallel = false;
        return compute_intersection_between_two_lines(startPointOfLineSegment1, endPointOfLineSegment1, startPointOfLineSegment2, endPointOfLineSegment2, parameterOfIntersectionForLineSeg1, parameterOfIntersectionForLineSeg2);
    }
}






rg_Point2D rg_GeoFunc::compute_intersection_between_two_lines(const rg_Point2D& startPointOfLineSegment1, const rg_Point2D& endPointOfLineSegment1, const rg_Point2D& startPointOfLineSegment2, const rg_Point2D& endPointOfLineSegment2, bool& bTwoLinesAreParallel)
{
    rg_REAL parameterOfIntersectionForLineSeg1 = DBL_MAX;
    rg_REAL parameterOfIntersectionForLineSeg2 = DBL_MAX;

    return compute_intersection_between_two_lines(startPointOfLineSegment1, endPointOfLineSegment1, startPointOfLineSegment2, endPointOfLineSegment2, parameterOfIntersectionForLineSeg1, parameterOfIntersectionForLineSeg2, bTwoLinesAreParallel);
}


rg_Point2D rg_GeoFunc::compute_intersection_between_two_rays(const rg_Point2D& startPointOfRay1, const rg_Point2D& directionVecOfRay1,
    const rg_Point2D& startPointOfRay2, const rg_Point2D& directionVecOfRay2,
    rg_REAL& parameterOfIntersectionForRay1, rg_REAL& parameterOfIntersectionForRay2)
{
    rg_REAL d1x = directionVecOfRay1.getX();
    rg_REAL d1y = directionVecOfRay1.getY();
    rg_REAL startX1 = startPointOfRay1.getX();
    rg_REAL startY1 = startPointOfRay1.getY();

    rg_REAL d2x = directionVecOfRay2.getX();
    rg_REAL d2y = directionVecOfRay2.getY();
    rg_REAL startX2 = startPointOfRay2.getX();
    rg_REAL startY2 = startPointOfRay2.getY();

    parameterOfIntersectionForRay1 = (d2x * startY1 - d2x * startY2 - d2y * startX1 + d2y * startX2) / (d1x * d2y - d1y * d2x);
    parameterOfIntersectionForRay2 = (d1x * startY1 - d1x * startY2 - d1y * startX1 + d1y * startX2) / (d1x * d2y - d1y * d2x);

    return (startPointOfRay1 + parameterOfIntersectionForRay1*directionVecOfRay1);
}

rg_Point2D rg_GeoFunc::compute_intersection_between_two_rays(const rg_Point2D& startPointOfRay1, const rg_Point2D& directionVecOfRay1,
    const rg_Point2D& startPointOfRay2, const rg_Point2D& directionVecOfRay2)
{
    rg_REAL parameterOfIntersectionForRay1 = -1.;
    rg_REAL parameterOfIntersectionForRay2 = -1.;
    return compute_intersection_between_two_rays(startPointOfRay1, directionVecOfRay1, startPointOfRay2, directionVecOfRay2, parameterOfIntersectionForRay1, parameterOfIntersectionForRay2);
}


rg_INT rg_GeoFunc::compute_intersection_between_circle_and_line_segment(const rg_Circle2D& circle, const rg_Line2D& lineSegment, rg_Point2D intersection[])
{
    rg_Point2D startPoint = lineSegment.getSP();
    rg_Point2D endPoint = lineSegment.getEP();
    rg_Point2D center = circle.getCenterPt();
    rg_REAL radius = circle.getRadius();

    rg_REAL startX = startPoint.getX();
    rg_REAL startY = startPoint.getY();
    rg_REAL endX = endPoint.getX();
    rg_REAL endY = endPoint.getY();
    rg_REAL centerX = center.getX();
    rg_REAL centerY = center.getY();

    rg_REAL coeff[3] = { 0.0, 0.0, 0.0 };
    coeff[0] = pow(startX - centerX, 2.0) + pow(startY - centerY, 2.0) - radius * radius;;
    coeff[1] = 2.0 * (startX - centerX) * (-startX + endX) + 2.0 * (startY - centerY) * (-startY + endY);
    coeff[2] = pow(-startX + endX, 2.0) + pow(-startY + endY, 2.0);

    rg_QuadraticPolynomial qp(coeff);
    rg_ComplexNumber* root = qp.solve();
    rg_INT numberOfIntersections = 0;

    for (rg_INT i = 0; i < 2; i++)
    {
        if (root[i].isPureRealNumber())
        {
            rg_REAL r = root[i].getRealNumber();
            if (rg_GE(r, 0.0) && rg_LE(r, 1.0))
            {
                intersection[numberOfIntersections] = (1. - r)*startPoint + r*endPoint;
                numberOfIntersections++;
            }
        }
    }

    if (root != rg_NULL)
        delete[] root;

    return numberOfIntersections;
}


bool rg_GeoFunc::compute_intersection_between_circle_and_ellipse(const rg_Circle2D& circle, Ellipse2D& ellipse, rg_Point2D intersection[])
{
    rg_Point2D center = circle.getCenterPt();
    rg_REAL    radius = circle.getRadius();
    rg_Point2D footPrintOfCircleCenterOntoEllipse;
    ellipse.compute_perpendicular_footprint_of_point_onto_ellipse(center, footPrintOfCircleCenterOntoEllipse);

    if (rg_GT(center.distance(footPrintOfCircleCenterOntoEllipse), radius))
        return false;
    else
        return true;
}

rg_INT rg_GeoFunc::compute_intersection_between_circle_and_point(const rg_Circle2D& circle, const rg_Point2D point)
{
    rg_Point2D centerPoint = circle.getCenterPt();
    rg_REAL radius = circle.getRadius();
    if (rg_LT(centerPoint.distance(point), radius))
        return true;
    else
        return false;
}


rg_Line2D rg_GeoFunc::compute_bisector_line_between_two_line_segments(const rg_Line2D& lineSegment1, const rg_Line2D& lineSegment2)
{
    bool bTwoLinesAreParallel = false;
    rg_Point2D intersecetion = lineSegment1.compute_intersection_with_line(lineSegment2, bTwoLinesAreParallel);
	rg_Point2D dirVec = ( lineSegment1.evaluateVector().getUnitVector() + lineSegment2.evaluateVector().getUnitVector() ).getUnitVector();

    if (bTwoLinesAreParallel)
    {
        intersecetion = (lineSegment1.getSP() + lineSegment2.getSP()) / 2.0;
		dirVec = lineSegment1.evaluateVector().getUnitVector();
    }    

    return rg_Line2D(intersecetion, intersecetion + dirVec);
}


rg_Line2D rg_GeoFunc::compute_bisector_line_between_two_points(const rg_Point2D& point1, const rg_Point2D& point2)
{
    rg_Line2D  line(point1, point2);
    rg_Point2D midPoint = (point1 + point2) / 2.0;
    rg_Point2D pointOnBisecetor = midPoint + line.getNormalVector().getUnitVector();

    return rg_Line2D(midPoint, pointOnBisecetor);
}

rg_INT rg_GeoFunc::compute_intersection_between_ellipse_and_line(const Ellipse2D& ellipse, const rg_Line2D& lineSegment, rg_Point2D intersection[])
{
    const rg_INT NUM_COEFF_ELLIPSE_EQ = 6;
    rg_REAL ellipseCoefficient[NUM_COEFF_ELLIPSE_EQ];
    const_cast<Ellipse2D&>(ellipse).get_coefficients_of_ellipse_equation(ellipseCoefficient);
    rg_Point2D startPoint = lineSegment.getSP();
    rg_Point2D endPoint = lineSegment.getEP();

    rg_REAL startX = startPoint.getX();
    rg_REAL startY = startPoint.getY();
    rg_REAL endX = endPoint.getX();
    rg_REAL endY = endPoint.getY();
    rg_REAL A = ellipseCoefficient[5]; rg_REAL B = ellipseCoefficient[4]; rg_REAL C = ellipseCoefficient[3];
    rg_REAL D = ellipseCoefficient[2]; rg_REAL E = ellipseCoefficient[1]; rg_REAL F = ellipseCoefficient[0];

    rg_REAL polynomialCoeff[3];
    polynomialCoeff[2] = A * pow(endX - startX, 2.0) + B * (endX - startX) * (endY - startY) + C * pow(endY - startY, 2.0);
    polynomialCoeff[1] = 2 * A * startX * (endX - startX) + B * startX * (endY - startY) + B * (endX - startX) * startY + 2 * C * startY * (endY - startY) + D * (endX - startX) + E * (endY - startY);
    polynomialCoeff[0] = A * startX * startX + B * startX * startY + C * startY * startY + D * startX + E * startY + F;

    rg_QuadraticPolynomial quadraticPolynomial(polynomialCoeff);
    rg_ComplexNumber* solution = rg_NULL;
    solution = quadraticPolynomial.solve();
    rg_INT numberOfRealRoots = 0;
    if (solution[0].isPureRealNumber())
    {
        intersection[numberOfRealRoots] = (1. - solution[0].getRealNumber())*startPoint + solution[0].getRealNumber()*endPoint;
        numberOfRealRoots++;
    }
    if (solution[1].isPureRealNumber())
    {
        intersection[numberOfRealRoots] = (1. - solution[1].getRealNumber())*startPoint + solution[1].getRealNumber()*endPoint;
        numberOfRealRoots++;
    }

    if (solution != rg_NULL)
        delete[] solution;

    return numberOfRealRoots;
}

rg_Point2D rg_GeoFunc::compute_a_point_on_line_whose_distance_with_anchor_point_on_line_is_equal_to_distance_between_the_point_and_its_footprint_of_ellipse(const rg_Point2D seedPoint, const rg_Line2D& line, const rg_Point2D& anchorPointOnLine, const Ellipse2D& ellipse, const rg_REAL& tolerance, const rg_INT& maximumNumberOfIterations)
{
    // test: VVetex is little bit away from perpendicular bisector between reflex vertex and edge (line segment) because this VVertex is not yet refined
    //rg_Point2D testFootPrint;
    //line.compute_perpendicular_footprint_of_point_onto_entire_line(seedPoint, testFootPrint);
    //rg_REAL dist = testFootPrint.distance(seedPoint);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rg_Point2D currentPoint;
    rg_Point2D targetPoint = seedPoint;
    rg_INT numberOfIterations = 0;
    do
    {
        currentPoint = targetPoint;
        rg_Point2D footPrintOfCurrentPointOntoEllipse;
        ellipse.compute_perpendicular_footprint_of_point_onto_ellipse(currentPoint, footPrintOfCurrentPointOntoEllipse);
        rg_Line2D bisectorLineBetweenFootPrintAndAnchorPoint = rg_GeoFunc::compute_bisector_line_between_two_points(footPrintOfCurrentPointOntoEllipse, anchorPointOnLine);
        bool bTwoLinesAreParallel = false;
        targetPoint = bisectorLineBetweenFootPrintAndAnchorPoint.compute_intersection_with_line(line, bTwoLinesAreParallel);
        ++numberOfIterations;
    } while (targetPoint.distance(currentPoint) >= tolerance && numberOfIterations <= maximumNumberOfIterations);

    return (currentPoint + targetPoint) / 2.0;


    /*
    We cannot gurantee that seedPoint is closer to anchorPoint than footPrintOfLine...
    Therefore, we cannot find the point by bisection method...????

    rg_Point2D footPrintOfEllipse;
    ellipse.compute_perpendicular_footprint_of_point_onto_ellipse(seedPoint, footPrintOfEllipse);
    rg_Point2D footPrintOfLine;
    line.compute_perpendicular_footprint_of_point_onto_entire_line(footPrintOfEllipse, footPrintOfLine);
    rg_REAL distance2AnchorPointOnLine  = seedPoint.distance(anchorPointOnLine);
    rg_REAL distance2FootPrintOfEllipse = seedPoint.distance(footPrintOfEllipse);

    if (rg_EQ(distance2AnchorPointOnLine, distance2FootPrintOfEllipse))
    return seedPoint;

    rg_Point2D oldMidPoint = (anchorPointOnLine + footPrintOfLine) / 2.0;
    rg_Point2D newMidPoint;

    do
    {
    newMidPoint = (oldMidPoint + footPrintOfLine) / 2.0;
    rg_Point2D temporaryOldMidPoint = newMidPoint;
    ellipse.compute_perpendicular_footprint_of_point_onto_ellipse(newMidPoint, footPrintOfEllipse);

    line.compute_perpendicular_footprint_of_point_onto_entire_line(footPrintOfEllipse, footPrintOfLine);

    distance2AnchorPointOnLine = newMidPoint.distance(anchorPointOnLine);
    distance2FootPrintOfEllipse = newMidPoint.distance(footPrintOfEllipse);

    if (rg_GT(distance2AnchorPointOnLine, distance2FootPrintOfEllipse))
    {
    newMidPoint = (oldMidPoint + newMidPoint) / 2.0;
    }
    else
    {
    newMidPoint = (footPrintOfLine + newMidPoint) / 2.0;
    }
    oldMidPoint = temporaryOldMidPoint;
    } while (rg_NZERO(oldMidPoint.distance(newMidPoint)));

    return (oldMidPoint + newMidPoint) / 2.0;
    */
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const Ellipse2D& ellipse2, const Ellipse2D& ellipse3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Point2D targetCenter = seedCircle.getCenterPt();
    rg_Point2D currentCenter;
    rg_Circle2D refinedTangentCircle;
    rg_INT numIterations = 0;
    do
    {
        currentCenter = targetCenter;
        rg_Point2D footprint[3];
        ellipse1.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[0]);
        ellipse2.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[1]);
        ellipse3.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[2]);
        refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
        targetCenter = refinedTangentCircle.getCenterPt();
        ++numIterations;
    } while (targetCenter.distance(currentCenter) >= TOLERANCE_FOR_VERTEX_CORRECTION && numIterations <= MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION);

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const Ellipse2D& ellipse2, const rg_Line2D& lineSeg3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Point2D targetCenter = seedCircle.getCenterPt();
    rg_Point2D currentCenter;
    rg_Circle2D refinedTangentCircle;
    rg_INT numIterations = 0;
    do
    {
        currentCenter = targetCenter;
        rg_Point2D footprint[3];
        ellipse1.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[0]);
        ellipse2.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[1]);
        lineSeg3.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[2]);
        refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
        targetCenter = refinedTangentCircle.getCenterPt();
        ++numIterations;
    } while (targetCenter.distance(currentCenter) >= TOLERANCE_FOR_VERTEX_CORRECTION && numIterations <= MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION);

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const Ellipse2D& ellipse2, const rg_Point2D& point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Point2D targetCenter = seedCircle.getCenterPt();
    rg_Point2D currentCenter;
    rg_Circle2D refinedTangentCircle;
    rg_INT numIterations = 0;
    do
    {
        currentCenter = targetCenter;
        rg_Point2D footprint[3];
        ellipse1.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[0]);
        ellipse2.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[1]);
        point3.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[2]);
        refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
        targetCenter = refinedTangentCircle.getCenterPt();
        ++numIterations;
    } while (targetCenter.distance(currentCenter) >= TOLERANCE_FOR_VERTEX_CORRECTION && numIterations <= MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION);

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const rg_Line2D& lineSeg2, const rg_Line2D& lineSeg3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Point2D targetCenter = seedCircle.getCenterPt();
    rg_Point2D currentCenter;
    rg_Circle2D refinedTangentCircle;
    rg_INT numIterations = 0;
    do
    {
        currentCenter = targetCenter;
        rg_Point2D footprint[3];
        ellipse1.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[0]);
        lineSeg2.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[1]);
        lineSeg3.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[2]);
        refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
        targetCenter = refinedTangentCircle.getCenterPt();
        ++numIterations;
    } while (targetCenter.distance(currentCenter) >= TOLERANCE_FOR_VERTEX_CORRECTION && numIterations <= MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION);

    //rg_REAL dist1, dist2, dist3;
    //rg_Point2D footprint[3];
    //ellipse1.compute_perpendicular_footprint_of_point_onto_ellipse(targetCenter, footprint[0]);
    //lineSeg2.compute_footprint_of_point_onto_line_segment(targetCenter, footprint[1]);
    //lineSeg3.compute_footprint_of_point_onto_line_segment(targetCenter, footprint[2]);
    //dist1 = targetCenter.distance(footprint[0]);
    //dist2 = targetCenter.distance(footprint[1]);
    //dist3 = targetCenter.distance(footprint[2]);

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const rg_Line2D& lineSeg2, const rg_Point2D& point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Point2D targetCenter = seedCircle.getCenterPt();
    rg_Point2D currentCenter;
    rg_Circle2D refinedTangentCircle;
    rg_INT numIterations = 0;
    do
    {
        currentCenter = targetCenter;
        rg_Point2D footprint[3];
        ellipse1.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[0]);
        lineSeg2.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[1]);
        point3.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[2]);
        refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
        targetCenter = refinedTangentCircle.getCenterPt();
        ++numIterations;
    } while (targetCenter.distance(currentCenter) >= TOLERANCE_FOR_VERTEX_CORRECTION && numIterations <= MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION);

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const Ellipse2D& ellipse1, const rg_Point2D&  point2, const rg_Point2D& point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Point2D targetCenter = seedCircle.getCenterPt();
    rg_Point2D currentCenter;
    rg_Circle2D refinedTangentCircle;
    rg_INT numIterations = 0;
    do
    {
        currentCenter = targetCenter;
        rg_Point2D footprint[3];
        ellipse1.compute_perpendicular_footprint_of_point_onto_ellipse(currentCenter, footprint[0]);
        point2.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[1]);
        point3.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[2]);
        refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
        targetCenter = refinedTangentCircle.getCenterPt();
        ++numIterations;
    } while (targetCenter.distance(currentCenter) >= TOLERANCE_FOR_VERTEX_CORRECTION && numIterations <= MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION);

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const rg_Line2D& lineSeg1, const rg_Line2D& lineSeg2, const rg_Line2D& lineSeg3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    //rg_Point2D targetCenter = seedCircle.getCenterPt();
    //rg_Point2D currentCenter;
    //rg_Circle2D refinedTangentCircle;
    //rg_INT numIterations = 0;
    //do
    //{
    //	currentCenter = targetCenter;
    //	rg_Point2D footprint[3];
    //	lineSeg1.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[0]);
    //	lineSeg2.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[1]);
    //	lineSeg3.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[2]);
    //	refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
    //	targetCenter = refinedTangentCircle.getCenterPt();
    //	++numIterations;
    //} while (targetCenter.distance(currentCenter) >= TOLERANCE_FOR_VERTEX_CORRECTION && numIterations <= MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION);
    //return refinedTangentCircle;

    rg_Line2D bisector12 = rg_GeoFunc::compute_bisector_line_between_two_line_segments(lineSeg1.get_reversed_line2D(), lineSeg2);
    rg_Line2D bisector23 = rg_GeoFunc::compute_bisector_line_between_two_line_segments(lineSeg2.get_reversed_line2D(), lineSeg3);
    bool bTwoLinesAreParllel = false;
    rg_Point2D refinedCenter = bisector12.compute_intersection_with_line(bisector23, bTwoLinesAreParllel);
    rg_REAL    refinedRadius = lineSeg1.getDistance(refinedCenter);
    rg_Circle2D refinedTangentCircle(refinedCenter, refinedRadius);

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const rg_Line2D& lineSeg1, const rg_Line2D& lineSeg2, const rg_Point2D&  point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Point2D targetCenter = seedCircle.getCenterPt();
    rg_Point2D currentCenter;
    rg_Circle2D refinedTangentCircle;
    rg_INT numIterations = 0;
    do
    {
        currentCenter = targetCenter;
        rg_Point2D footprint[3];
        lineSeg1.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[0]);
        lineSeg2.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[1]);
        point3.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[2]);
        refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
        targetCenter = refinedTangentCircle.getCenterPt();
        ++numIterations;
    } while (targetCenter.distance(currentCenter) >= TOLERANCE_FOR_VERTEX_CORRECTION && numIterations <= MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION);

    double radius = refinedTangentCircle.getCenterPt().distance( point3 );
    refinedTangentCircle.setRadius( radius );

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const rg_Line2D& lineSeg1, const rg_Point2D&  point2, const rg_Point2D&  point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Point2D targetCenter = seedCircle.getCenterPt();
    rg_Point2D currentCenter;
    rg_Circle2D refinedTangentCircle;
    rg_INT numIterations = 0;
    do
    {
        currentCenter = targetCenter;
        rg_Point2D footprint[3];
        lineSeg1.compute_footprint_of_point_onto_line_segment(currentCenter, footprint[0]);
        point2.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[1]);
        point3.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[2]);
        refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);
        targetCenter = refinedTangentCircle.getCenterPt();
        ++numIterations;
    } while (targetCenter.distance(currentCenter) >= TOLERANCE_FOR_VERTEX_CORRECTION && numIterations <= MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION);

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const rg_Point2D&  point1, const rg_Point2D&  point2, const rg_Point2D&  point3, const rg_Circle2D& seedCircle, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Point2D currentCenter = seedCircle.getCenterPt();

    rg_Circle2D refinedTangentCircle;
    rg_Point2D footprint[3];
    point1.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[0]);
    point2.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[1]);
    point3.compute_perpendicular_footprint_of_point_onto_this_point(currentCenter, footprint[2]);
    refinedTangentCircle.setCircleWithThreePassingPoints(footprint[0], footprint[1], footprint[2]);

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const rg_Circle2D& disk1, const rg_Circle2D& disk2, const rg_Line2D&  lineSeg3, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Circle2D refinedTangentCircle;

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const rg_Circle2D& disk1, const rg_Line2D& lineSeg2, const rg_Line2D&  lineSeg3, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Circle2D refinedTangentCircle;

    return refinedTangentCircle;
}

rg_Circle2D rg_GeoFunc::compute_empty_tangent_cirlce_of_three_generators(const rg_Circle2D& disk1, const rg_Line2D& lineSeg2, const rg_Point2D& point3, const rg_REAL TOLERANCE_FOR_VERTEX_CORRECTION, const rg_INT  MAXIMUM_ITERATIONS_FOR_VERTEX_CORRECTION)
{
    rg_Circle2D refinedTangentCircle;

    return refinedTangentCircle;
}

rg_INT rg_GeoFunc::compute_intersection_between_two_circles( const rg_Circle2D & circle1, const rg_Circle2D & circle2, rg_Point2D * intersectionPt )
{
    rg_REAL radius1 = circle1.getRadius();
    rg_REAL radius2 = circle2.getRadius();
    rg_Point2D center1 = circle1.getCenterPt();
    rg_Point2D center2 = circle2.getCenterPt();

    rg_Point2D vecCenter12 = center2 - center1;
    rg_Point2D vecOrthogonalToVecCenter12(vecCenter12.getY(), -vecCenter12.getX());
    rg_REAL    radiu1_plus_radius2 = radius1 + radius2;
    rg_REAL    radiu1_minus_radius2 = radius1 - radius2;
    rg_REAL    normOfVecCenter12_squared = vecCenter12.magnitudeSquare();
    rg_REAL    parameterOfVecCenter12 = (radiu1_plus_radius2 * radiu1_minus_radius2 / normOfVecCenter12_squared + 1.0) * 0.5;
    rg_REAL    parameterOfVecOrthogonalToVecCenter12_squared = (-normOfVecCenter12_squared + radiu1_plus_radius2*radiu1_plus_radius2) * (normOfVecCenter12_squared - radiu1_minus_radius2*radiu1_minus_radius2) / (4.0 * normOfVecCenter12_squared * normOfVecCenter12_squared);

    if (parameterOfVecOrthogonalToVecCenter12_squared < 0.0)
    {
        return 0;
    }
	else if (parameterOfVecOrthogonalToVecCenter12_squared == 0.0)
	{
		intersectionPt[0] = center1 + parameterOfVecCenter12 * vecCenter12;
		intersectionPt[1] = center1 + parameterOfVecCenter12 * vecCenter12;

		return 1;
	}
    else
    {
        rg_REAL  parameterOfVecOrthogonalToVecCenter12[2] = { sqrt(parameterOfVecOrthogonalToVecCenter12_squared), -sqrt(parameterOfVecOrthogonalToVecCenter12_squared) };

        intersectionPt[0] = center1 + parameterOfVecCenter12 * vecCenter12 + parameterOfVecOrthogonalToVecCenter12[0] * vecOrthogonalToVecCenter12;
        intersectionPt[1] = center1 + parameterOfVecCenter12 * vecCenter12 + parameterOfVecOrthogonalToVecCenter12[1] * vecOrthogonalToVecCenter12;

        return 2;
    }
}

rg_INT rg_GeoFunc::compute_tangent_circles_of_two_circles_with_given_radius( const rg_Circle2D & circle1, const rg_Circle2D & circle2, const rg_INT & containerCircleIndex, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle )
{
    rg_Circle2D circle1_offSetByGivenRadius = circle1;
    rg_Circle2D circle2_offSetByGivenRadius = circle2;

    switch (containerCircleIndex)
    {
    case 1:
        circle1_offSetByGivenRadius.setRadius(circle1.getRadius() - givenRadius);
        circle2_offSetByGivenRadius.setRadius(circle2.getRadius() + givenRadius);
        break;
    case 2:
        circle1_offSetByGivenRadius.setRadius(circle1.getRadius() + givenRadius);
        circle2_offSetByGivenRadius.setRadius(circle2.getRadius() - givenRadius);
        break;
    default:
        circle1_offSetByGivenRadius.setRadius(circle1.getRadius() + givenRadius);
        circle2_offSetByGivenRadius.setRadius(circle2.getRadius() + givenRadius);
        break;
    }

    rg_Point2D* intersectionPt = new rg_Point2D[2];
    rg_INT numIntersectionPts = compute_intersection_between_two_circles(circle1_offSetByGivenRadius, circle2_offSetByGivenRadius, intersectionPt);

    rg_INT numTangentCircles = 0;
    switch (numIntersectionPts)
    {
    case 1:
        tangentCircle[0] = rg_Circle2D(intersectionPt[0], givenRadius);
        numTangentCircles = 1;
        break;
    case 2:
        tangentCircle[0] = rg_Circle2D(intersectionPt[0], givenRadius);
        tangentCircle[1] = rg_Circle2D(intersectionPt[1], givenRadius);
        numTangentCircles = 2;
        break;
    default:
        break;
    }

    if (intersectionPt != rg_NULL)
        delete[] intersectionPt;

    return numTangentCircles;
}

rg_INT rg_GeoFunc::compute_tangent_circles_of_point_and_line_with_given_radius(const rg_Point2D & point, const rg_Line2D & line, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle)
{
	rg_REAL coeffX = 0.0, coeffY = 0.0, coeffConst = 0.0;
	line.get_coefficients_of_implicit_form_of_line_equation_in_normalized(coeffX, coeffY, coeffConst);

	rg_REAL pointX = point.getX();
	rg_REAL pointY = point.getY();

	rg_REAL constForPointTrans = coeffConst + coeffX * pointX + coeffY * pointY;

	// The point is on the line
	if (rg_ZERO(fabs(constForPointTrans), rg_MATH_RES))
	{
		rg_REAL x1 = pointX + coeffX * givenRadius;
		rg_REAL y1 = pointY + coeffY * givenRadius;
		rg_REAL x2 = pointX - coeffX * givenRadius;
		rg_REAL y2 = pointY - coeffY * givenRadius;

		tangentCircle[0].setCircle(rg_Point2D(x1, y1), givenRadius);
		tangentCircle[1].setCircle(rg_Point2D(x2, y2), givenRadius);

		return 2;
	}

	if (constForPointTrans < 0.0)
	{
		coeffX = -coeffX;
		coeffY = -coeffY;
		constForPointTrans = -constForPointTrans;
	}

	rg_REAL cTemp = constForPointTrans - givenRadius;
	rg_REAL root  = givenRadius * givenRadius - cTemp * cTemp;
	
	// no tangent circle exists
	if (root < 0.0)
	{
		return 0;
	}

	// two tangent circles
	if (root > 0.0)
	{
		root = sqrt(root);
		rg_REAL xConst = pointX - coeffX * cTemp;
		rg_REAL yConst = pointY - coeffY * cTemp;
		rg_REAL xVar = coeffY * root;
		rg_REAL yVar = coeffX * root;

		rg_REAL x1 = xConst + xVar;
		rg_REAL y1 = yConst - yVar;
		rg_REAL x2 = xConst - xVar;
		rg_REAL y2 = yConst + yVar;

		tangentCircle[0].setCircle(rg_Point2D(x1, y1), givenRadius);
		tangentCircle[1].setCircle(rg_Point2D(x2, y2), givenRadius);

		return 2;
	}
	// one tangent circle
	else
	{
		rg_REAL x = pointX - coeffX * cTemp;
		rg_REAL y = pointY - coeffY * cTemp;

		tangentCircle[0].setCircle(rg_Point2D(x, y), givenRadius);

		return 1;
	}
}

rg_INT rg_GeoFunc::compute_tangent_circles_of_point_and_linesegment_with_given_radius(const rg_Point2D & point, const rg_Line2D & lineSegment, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle)
{
	rg_INT numTangentCircles = 0;
	numTangentCircles = rg_GeoFunc::compute_tangent_circles_of_point_and_line_with_given_radius(point, lineSegment, givenRadius, tangentCircle);

	switch (numTangentCircles)
	{
	case 1:
	{
		rg_REAL parameter = DBL_MAX;
		lineSegment.compute_perpendicular_footprint_of_point_onto_entire_line(tangentCircle[0].getCenterPt(), parameter);
		if (0.0 <= parameter && parameter <= 1.0)
			return 1;
		else
			return 0;
	}
	break;
	case 2:
	{
		bool tangentToLineSeg[2] = { false, false };
		rg_INT numTangentCirclesOfPointAndLineSeg = 0;
		rg_REAL parameter = DBL_MAX;
		lineSegment.compute_perpendicular_footprint_of_point_onto_entire_line(tangentCircle[0].getCenterPt(), parameter);
		if (0.0 <= parameter && parameter <= 1.0)
		{
			tangentToLineSeg[0] = true;
			++numTangentCirclesOfPointAndLineSeg;
		}

		parameter = DBL_MAX;
		lineSegment.compute_perpendicular_footprint_of_point_onto_entire_line(tangentCircle[1].getCenterPt(), parameter);
		if (0.0 <= parameter && parameter <= 1.0)
		{
			tangentToLineSeg[1] = true;
			++numTangentCirclesOfPointAndLineSeg;
		}

		if (numTangentCirclesOfPointAndLineSeg == 1 && tangentToLineSeg[1] == true)
		{
			tangentCircle[0] = tangentCircle[1];
		}

		return numTangentCirclesOfPointAndLineSeg;
	}
	break;
	default:
		return 0;
	}
}

rg_INT rg_GeoFunc::compute_tangent_circles_of_circle_and_line_with_given_radius(const rg_Circle2D & circle, const rg_Line2D & line, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle)
{
	rg_REAL coeffX = 0.0, coeffY = 0.0, coeffConst = 0.0;
	line.get_coefficients_of_implicit_form_of_line_equation_in_normalized(coeffX, coeffY, coeffConst);

	rg_Point2D circleCenter = circle.getCenterPt();
	rg_REAL centerX = circleCenter.getX();
	rg_REAL centerY = circleCenter.getY();

	rg_REAL constForPointTrans = coeffConst + coeffX * centerX + coeffY * centerY;

	if (constForPointTrans < 0.0)
	{
		coeffX = -coeffX;
		coeffY = -coeffY;
		constForPointTrans = -constForPointTrans;
	}

	rg_REAL cTemp = constForPointTrans - givenRadius;
	rg_REAL rTemp = givenRadius + circle.getRadius();

	rg_REAL root = rTemp * rTemp - cTemp * cTemp;

	// no tangent circle exists
	if (root < 0.0)
	{
		return 0;
	}

	// two tangent circles
	if (root > 0.0)
	{
		root = sqrt(root);
		rg_REAL xConst = centerX - coeffX * cTemp;
		rg_REAL yConst = centerY - coeffY * cTemp;
		rg_REAL xVar = coeffY * root;
		rg_REAL yVar = coeffX * root;

		rg_REAL x1 = xConst - xVar;
		rg_REAL y1 = yConst + yVar;
		rg_REAL x2 = xConst + xVar;
		rg_REAL y2 = yConst - yVar;

		tangentCircle[0].setCircle(rg_Point2D(x1, y1), givenRadius);
		tangentCircle[1].setCircle(rg_Point2D(x2, y2), givenRadius);

		return 2;
	}
	// one tangent circle
	else
	{
		rg_REAL x = centerX + coeffX * cTemp;
		rg_REAL y = centerY + coeffY * cTemp;

		tangentCircle[0].setCircle(rg_Point2D(x, y), givenRadius);

		return 1;
	}
}

rg_INT rg_GeoFunc::compute_tangent_circles_of_circle_and_linesegment_with_given_radius(const rg_Circle2D & circle, const rg_Line2D & lineSegment, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle)
{
	rg_INT numTangentCircles = 0;
	numTangentCircles = rg_GeoFunc::compute_tangent_circles_of_circle_and_line_with_given_radius(circle, lineSegment, givenRadius, tangentCircle);

	switch (numTangentCircles)
	{
	case 1:
	{
		rg_REAL parameter = DBL_MAX;
		lineSegment.compute_perpendicular_footprint_of_point_onto_entire_line(tangentCircle[0].getCenterPt(), parameter);
		if (0.0 <= parameter && parameter <= 1.0)
			return 1;
		else
			return 0;
	}
	break;
	case 2:
	{
		bool tangentToLineSeg[2] = { false, false };
		rg_INT numTangentCirclesOfCircleAndLineSeg = 0;
		rg_REAL parameter = DBL_MAX;
		lineSegment.compute_perpendicular_footprint_of_point_onto_entire_line(tangentCircle[0].getCenterPt(), parameter);
		if (0.0 <= parameter && parameter <= 1.0)
		{
			tangentToLineSeg[0] = true;
			++numTangentCirclesOfCircleAndLineSeg;
		}

		parameter = DBL_MAX;
		lineSegment.compute_perpendicular_footprint_of_point_onto_entire_line(tangentCircle[1].getCenterPt(), parameter);
		if (0.0 <= parameter && parameter <= 1.0)
		{
			tangentToLineSeg[1] = true;
			++numTangentCirclesOfCircleAndLineSeg;
		}

		if (numTangentCirclesOfCircleAndLineSeg == 1 && tangentToLineSeg[1] == true)
		{
			tangentCircle[0] = tangentCircle[1];
		}

		return numTangentCirclesOfCircleAndLineSeg;
	}
	break;
	default:
		return 0;
	}
}

rg_INT rg_GeoFunc::compute_tangent_circles_of_two_lines_not_parallel_to_each_other_with_given_radius(const rg_Line2D & line1, const rg_Line2D & line2, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle)
{
	rg_REAL coeffX1 = 0.0, coeffY1 = 0.0, coeffConst1 = 0.0;
	line1.get_coefficients_of_implicit_form_of_line_equation_in_normalized(coeffX1, coeffY1, coeffConst1);
	rg_REAL coeffX2 = 0.0, coeffY2 = 0.0, coeffConst2 = 0.0;
	line2.get_coefficients_of_implicit_form_of_line_equation_in_normalized(coeffX2, coeffY2, coeffConst2);
	
	rg_REAL determinant = coeffX2 * coeffY1 - coeffX1 * coeffY2;

	// Two lines are parallel.
	if (rg_ZERO(fabs(determinant), rg_MATH_RES))
	{
		return 0;
	}

	rg_INT numTangentCircles = 4;

	rg_REAL line1CoeffSumOfSquared_Sqrt = sqrt(coeffX1*coeffX1 + coeffY1*coeffY1);
	rg_REAL line2CoeffSumOfSquared_Sqrt = sqrt(coeffX2*coeffX2 + coeffY2*coeffY2);

	// sol. 1
	rg_REAL line1_const_radius_line1CoeffSumOfSquared_Sqrt = coeffConst1 + givenRadius * line1CoeffSumOfSquared_Sqrt;
	rg_REAL line2_const_radius_line2CoeffSumOfSquared_Sqrt = coeffConst2 + givenRadius * line2CoeffSumOfSquared_Sqrt;
	rg_REAL line1_const_radius_line2CoeffSumOfSquared_Sqrt = coeffConst1 + givenRadius * line2CoeffSumOfSquared_Sqrt;
	rg_REAL line2_const_radius_line1CoeffSumOfSquared_Sqrt = coeffConst2 + givenRadius * line1CoeffSumOfSquared_Sqrt;
	rg_REAL x1 = (coeffY2*line1_const_radius_line2CoeffSumOfSquared_Sqrt - coeffY1 * line2_const_radius_line1CoeffSumOfSquared_Sqrt) / determinant;
	rg_REAL y1 = (coeffX1*line2_const_radius_line2CoeffSumOfSquared_Sqrt - coeffX2 * line1_const_radius_line1CoeffSumOfSquared_Sqrt) / determinant;
	tangentCircle[0].setCircle(rg_Point2D(x1, y1), givenRadius);

	// sol. 2
	        line1_const_radius_line1CoeffSumOfSquared_Sqrt = coeffConst1 - givenRadius * line1CoeffSumOfSquared_Sqrt;
			line2_const_radius_line2CoeffSumOfSquared_Sqrt = coeffConst2 - givenRadius * line2CoeffSumOfSquared_Sqrt;
			line1_const_radius_line2CoeffSumOfSquared_Sqrt = coeffConst1 + givenRadius * line2CoeffSumOfSquared_Sqrt;
			line2_const_radius_line1CoeffSumOfSquared_Sqrt = coeffConst2 + givenRadius * line1CoeffSumOfSquared_Sqrt;
	rg_REAL x2 = (coeffY2*line1_const_radius_line2CoeffSumOfSquared_Sqrt - coeffY1 * line2_const_radius_line1CoeffSumOfSquared_Sqrt) / determinant;
	rg_REAL y2 = (coeffX1*line2_const_radius_line2CoeffSumOfSquared_Sqrt - coeffX2 * line1_const_radius_line1CoeffSumOfSquared_Sqrt) / determinant;
	tangentCircle[1].setCircle(rg_Point2D(x2, y2), givenRadius);

	// sol. 3
	        line1_const_radius_line1CoeffSumOfSquared_Sqrt = coeffConst1 + givenRadius * line1CoeffSumOfSquared_Sqrt;
			line2_const_radius_line2CoeffSumOfSquared_Sqrt = coeffConst2 + givenRadius * line2CoeffSumOfSquared_Sqrt;
			line1_const_radius_line2CoeffSumOfSquared_Sqrt = coeffConst1 - givenRadius * line2CoeffSumOfSquared_Sqrt;
			line2_const_radius_line1CoeffSumOfSquared_Sqrt = coeffConst2 - givenRadius * line1CoeffSumOfSquared_Sqrt;
	rg_REAL x3 = (coeffY2*line1_const_radius_line2CoeffSumOfSquared_Sqrt - coeffY1 * line2_const_radius_line1CoeffSumOfSquared_Sqrt) / determinant;
	rg_REAL y3 = (coeffX1*line2_const_radius_line2CoeffSumOfSquared_Sqrt - coeffX2 * line1_const_radius_line1CoeffSumOfSquared_Sqrt) / determinant;
	tangentCircle[2].setCircle(rg_Point2D(x3, y3), givenRadius);

	// sol. 4
	        line1_const_radius_line1CoeffSumOfSquared_Sqrt = coeffConst1 - givenRadius * line1CoeffSumOfSquared_Sqrt;
			line2_const_radius_line2CoeffSumOfSquared_Sqrt = coeffConst2 - givenRadius * line2CoeffSumOfSquared_Sqrt;
			line1_const_radius_line2CoeffSumOfSquared_Sqrt = coeffConst1 - givenRadius * line2CoeffSumOfSquared_Sqrt;
			line2_const_radius_line1CoeffSumOfSquared_Sqrt = coeffConst2 - givenRadius * line1CoeffSumOfSquared_Sqrt;
	rg_REAL x4 = (coeffY2*line1_const_radius_line2CoeffSumOfSquared_Sqrt - coeffY1 * line2_const_radius_line1CoeffSumOfSquared_Sqrt) / determinant;
	rg_REAL y4 = (coeffX1*line2_const_radius_line2CoeffSumOfSquared_Sqrt - coeffX2 * line1_const_radius_line1CoeffSumOfSquared_Sqrt) / determinant;
	tangentCircle[3].setCircle(rg_Point2D(x4, y4), givenRadius);

	return numTangentCircles;
}

rg_INT rg_GeoFunc::compute_tangent_circles_of_two_linesegments_not_parallel_to_each_other_with_given_radius(const rg_Line2D& lineSegment1, const rg_Line2D& lineSegment2, const rg_REAL & givenRadius, rg_Circle2D * tangentCircle)
{
	rg_INT numTangentCircles = 0;
	tangentCircle = new rg_Circle2D[4];
	numTangentCircles = rg_GeoFunc::compute_tangent_circles_of_two_lines_not_parallel_to_each_other_with_given_radius(lineSegment1, lineSegment2, givenRadius, tangentCircle);

	bool tangentToLineSeg[4] = { false, false , false, false };
	if (numTangentCircles > 0)
	{
		numTangentCircles = 0;

		rg_REAL parameter1 = -1.0;
		lineSegment1.project(tangentCircle[0].getCenterPt(), parameter1);
		rg_REAL parameter2 = -1.0;
		lineSegment2.project(tangentCircle[0].getCenterPt(), parameter2);
		if ( (0.0 <= parameter1 && parameter1 <= 1.0) && (0.0 <= parameter2 && parameter2 <= 1.0) )
		{
			tangentToLineSeg[0] = true;
			++numTangentCircles;
		}
		parameter1 = -1.0;
		lineSegment1.project(tangentCircle[1].getCenterPt(), parameter1);
		parameter2 = -1.0;
		lineSegment2.project(tangentCircle[1].getCenterPt(), parameter2);
		if ((0.0 <= parameter1 && parameter1 <= 1.0) && (0.0 <= parameter2 && parameter2 <= 1.0))
		{
			tangentToLineSeg[1] = true;
			++numTangentCircles;
		}
		parameter1 = -1.0;
		lineSegment1.project(tangentCircle[2].getCenterPt(), parameter1);
		parameter2 = -1.0;
		lineSegment2.project(tangentCircle[2].getCenterPt(), parameter2);
		if ((0.0 <= parameter1 && parameter1 <= 1.0) && (0.0 <= parameter2 && parameter2 <= 1.0))
		{
			tangentToLineSeg[2] = true;
			++numTangentCircles;
		}
		parameter1 = -1.0;
		lineSegment1.project(tangentCircle[3].getCenterPt(), parameter1);
		parameter2 = -1.0;
		lineSegment2.project(tangentCircle[3].getCenterPt(), parameter2);
		if ((0.0 <= parameter1 && parameter1 <= 1.0) && (0.0 <= parameter2 && parameter2 <= 1.0))
		{
			tangentToLineSeg[3] = true;
			++numTangentCircles;
		}
	}

	switch (numTangentCircles)
	{
	case 1:
	case 2:
	case 3:
		{
			rg_INT j = 0;
			for (rg_INT i = 0; i <= 3; ++i)
			{
				if (tangentToLineSeg[i])
				{
					tangentCircle[j] = tangentCircle[i];
					++j;
				}
			}
		}
		return numTangentCircles;
	default:
		return 0;
	}
}

bool rg_GeoFunc::point_in_polygon(const rg_Point2D& queryPoint, const list<rg_Point2D>& verticesOnShell)
{
	rg_Point2D vtx0 = verticesOnShell.back();
	// Get test bit for above/below X-axis (horizontal line containing the query point).
	bool yflag0 = (vtx0.getY() >= queryPoint.getY());

	bool inside_flag = false;
	for(const auto& vertex : verticesOnShell)
	{
		rg_Point2D vtx1 = vertex;
		bool yflag1 = (vtx1.getY() >= queryPoint.getY());
		// Check if endpoints straddle (are on opposite sides) of X-axis.
		// If so, +X ray could intersect this edge.
		if (yflag0 != yflag1) 
		{
			bool xflag0 = (vtx0.getX() >= queryPoint.getX());
			// Check if endpoints are on same side of the Y axis (vertical line containing the query point).
			// If so, it's easy to test if edge hits or misses.
			if (xflag0 == (vtx1.getX() >= queryPoint.getX())) 
			{
				// If edge's X values both right of the point, must hit.
				if (xflag0) 
					inside_flag = !inside_flag;
			}
			else 
			{
				// Compute intersection of polygon edge with +X ray.
				// If >= query point's X, the ray hits the edge.
				rg_REAL xCoodOfIntersection = (vtx1.getX() - (vtx1.getY() - queryPoint.getY()) * (vtx0.getX() - vtx1.getX()) / (vtx0.getY() - vtx1.getY()));
				if (xCoodOfIntersection >= queryPoint.getX())
				{
					inside_flag = !inside_flag;
				}
			}
		}
		// Move to next pair of vertices, retaining info.
		yflag0 = yflag1;
		vtx0 = vtx1;
	}
	return inside_flag;
}
