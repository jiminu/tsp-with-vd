//********************************************************************
//
//	  FILENAME    : rg_Line3D.cpp
//	  
//    DESCRIPTION : 
//           This consists of the implementation of class rg_Line3D.      
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 11 Jul. 1997    
//
//    HISTROY     :
//     1. By Young-Song Cho.   19 Jul. 1997
//         make : 
//            rg_FLAG isParallelToLine(const rg_Line3D &line) const;
//            rg_FLAG isPerpendicularToLine(const rg_Line3D &line) const;
//            rg_FLAG isThisPointOnLine(const rg_Point3D &pt) const;
//     2. By Tae-bum Jang   1998. 7.13
//         make : 
//            rg_Point3D getUnitDirection() const;
//            rg_REAL getDistance(const rg_Point3D& pt) const;
//

//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#include "rg_Line3D.h"
#include "rg_RelativeOp.h"

// Construction/Destruction
rg_Line3D::rg_Line3D()
{
}

rg_Line3D::rg_Line3D(const rg_Point3D &sPt, const rg_Point3D &ePt)
{
 	sp = sPt;
    ep = ePt;
}

rg_Line3D::rg_Line3D(const rg_Line3D &line)
{
    sp = line.sp;
    ep = line.ep;
}

rg_Line3D::~rg_Line3D()
{
}
            
	// Access elements	                       
rg_Point3D rg_Line3D::getSP() const
{
    return sp;
}

rg_Point3D rg_Line3D::getEP() const  
{
    return ep;
}

//	In 2D the following function is meaningful but in 3D is not.
//	rg_REAL  getTangent() const;
	
void rg_Line3D::setSP(const rg_Point3D &pt)
{
    sp = pt;
}

void rg_Line3D::setEP(const rg_Point3D &pt)
{
    ep = pt;
}

void rg_Line3D::setLine(const rg_Point3D &sPt, const rg_Point3D &ePt)
{
    sp = sPt;
    ep = ePt;
}

void rg_Line3D::setLine(const rg_Line3D &line)
{
    sp = line.sp;
    ep = line.ep;
}

// Operations
rg_REAL rg_Line3D::getLength() const
{
    return ep.distance(sp);
}

rg_REAL rg_Line3D::getDistance(const rg_Point3D& pt) const
{
    rg_REAL t=0.0;
    if ( rg_Line3D::isThisPointOnLine(pt,t) )
    {
        return 0.0;
    }

    // calculated cosine of the angle of start point
    rg_Point3D   roofVectorFromSP( pt-sp );
    rg_Point3D   lineVectorFromSP( ep-sp );
    rg_Point3D   unitRoofVectorFromSP = roofVectorFromSP.getUnitVector();
    rg_Point3D   unitLineVectorFromSP = lineVectorFromSP.getUnitVector();
    rg_REAL roofLengthFromSP            = roofVectorFromSP.magnitude();
    rg_REAL cosineFromSP = unitRoofVectorFromSP % unitLineVectorFromSP;

    // calculated cosine of the angle of end point
    rg_Point3D   roofVectorFromEP( pt-ep );
    rg_Point3D   lineVectorFromEP( sp-ep );
    rg_Point3D   unitRoofVectorFromEP = roofVectorFromEP.getUnitVector();
    rg_Point3D   unitLineVectorFromEP = lineVectorFromEP.getUnitVector();
    rg_REAL roofLengthFromEP            = roofVectorFromEP.magnitude();
    rg_REAL cosineFromEP = unitRoofVectorFromEP % unitLineVectorFromEP;

    rg_REAL output=0.0;

    if ( rg_LT(cosineFromSP,0.0) ) // The projection of point on the line is out of start point
    {
        output=pt.distance(sp);
    }
    else if ( rg_LT(cosineFromEP,0.0) ) // The projection of point on line is out of end point
    {
        output=pt.distance(ep);
    }
    else    // projection of point on line is between start point and end point
    {
        rg_REAL sineFromSP   = sqrt(1-cosineFromSP*cosineFromSP);
        output = roofLengthFromSP*sineFromSP;
    }

    return output;
}



rg_Point3D rg_Line3D::getUnitDirection() const
{
    rg_Point3D output=(ep-sp).getUnitVector();

    return output;
}

rg_FLAG rg_Line3D::isParallelToLine(const rg_Line3D &line) const
{
    if ( (ep - sp).getUnitVector() == (line.ep - line.sp).getUnitVector() )
    {
        return rg_TRUE;
    }
    else
        return rg_FALSE;
}

rg_FLAG rg_Line3D::isColinearWith(const rg_Line3D &line) const
{
    if (   ( isParallelToLine(line) == rg_TRUE )  
        && ( (ep - sp).getUnitVector() == (ep - line.sp).getUnitVector() ) )
    {
        return rg_TRUE;
    }
    else
        return rg_FALSE;
}


rg_FLAG rg_Line3D::isPerpendicularToLine(const rg_Line3D &line) const
{
    if ( !( (ep - sp)%(line.ep - line.sp) ) )
        return rg_TRUE;
    else
        return rg_FALSE;
}

rg_FLAG rg_Line3D::isThisPointOnLine(const rg_Point3D &pt, rg_REAL &u) const
{
    //  Parametric Form of rg_Line
    //      l(u) = (1-u)*sp + u*ep
    // 
    //  Is This point, pt, on rg_Line?
    //
    //    (1-u)*sp + u*ep = pt
    //    u = (pt - sp)/(ep - sp)
    //    0. <= u <= 1. ?

    rg_REAL u1 = (pt.getX() - sp.getX())/(ep.getX() - sp.getX());
    rg_REAL u2 = (pt.getY() - sp.getY())/(ep.getY() - sp.getY());
    rg_REAL u3 = (pt.getZ() - sp.getZ())/(ep.getZ() - sp.getZ());

    if ( rg_EQ(u1, u2) && rg_EQ(u2, u3) && rg_EQ(u3, u1)
         && rg_BTOR( 0.0, u1, 1.0) )
    {
        u = u1;
        return rg_TRUE;
    }
    else
    {
        u = -1.;    
        return rg_FALSE;
    }
}

// Operator Overloading
rg_FLAG rg_Line3D::operator==(const rg_Line3D &line) const
{
    if ( ((sp == line.sp) && (ep == line.ep))
         && ((sp == line.ep) && (ep == line.sp)) )
    {
        return rg_TRUE;
    }
    else
        return rg_FALSE;
}

// Intersection
rg_Point3D* rg_Line3D::intersectLine(const rg_Line3D &line, INTERSECT_TYPE &type) const
{
    //  Parametric Form of rg_Line
    //      l(u) = (1-u)*sp + u*ep
    //
    //  Intersection of line and line
    //      (1-u)*sp + u*ep = (1-v)*line.sp + v*line.ep
    //      u*(ep - sp) + v*(line.sp - line.ep) = line.sp - sp

    rg_REAL* param = solveSimultaneousEq(
                      (ep - sp),
                      (line.sp - line.ep),
                      (line.sp - sp) );

    rg_Point3D* intersectPt = rg_NULL;
    if ( param != rg_NULL )
    {
        if ( !rg_BTOR(0., param[0], 1.) || !rg_BTOR(0., param[1], 1.) )
        {
            type = NON_INTERSECT;

            return rg_NULL;
        }

        intersectPt = new rg_Point3D(( (1-param[0]) * sp ) + (param[0] * ep) );
        type = POINT_INTERSECT;

        delete [] param;

        return intersectPt;
    }
    else
    {
        if ( isParallelToLine(line) )
        {
            rg_REAL u = 0.;
            rg_REAL v = 0.;
            if ( line.isThisPointOnLine(sp, u) && line.isThisPointOnLine(ep, v) )
            {
                intersectPt    = new rg_Point3D[2];
                intersectPt[0] = sp;
                intersectPt[1] = ep;

                type = LINE_INTERSECT;
            }
            else if ( line.isThisPointOnLine(sp, u) )
            {
                if ( rg_EQ( u, 0.) || rg_EQ(u, 1.) )
                {
                    intersectPt    = new rg_Point3D(sp);
                    type = POINT_INTERSECT;

                }
                else
                {
                    intersectPt    = new rg_Point3D[2];
                    intersectPt[0] = sp;
                    if ( isThisPointOnLine(line.sp, v) )
                        intersectPt[1] = line.sp;
                    else 
                        intersectPt[1] = line.ep;
                    type = LINE_INTERSECT;
                }
            }
            else if ( line.isThisPointOnLine(ep, v) )
            {
                if ( rg_EQ( v, 0.) || rg_EQ(v, 1.) )
                {
                    intersectPt    = new rg_Point3D(ep);
                    type = POINT_INTERSECT;

                }
                else
                {
                    intersectPt    = new rg_Point3D[2];
                    intersectPt[0] = ep;
                    if ( isThisPointOnLine(line.sp, u) )
                        intersectPt[1] = line.sp;
                    else 
                        intersectPt[1] = line.ep;
                    type = LINE_INTERSECT;
                }
            }
            else
                type = NON_INTERSECT;
        }
        else
            type = NON_INTERSECT;

        return intersectPt;
    }
}

//  uA + vB = C
//  find u, v 
//  The first returned values is u, the second v.
//  After you use this function, you must check the validity of u and v
//      ? 0.0 <= u, v <= 1.0
rg_REAL* solveSimultaneousEq(const rg_Point3D &A, const rg_Point3D &B, const rg_Point3D &C) 
{
    //  1. u*A.x + v*B.x = C.x 
    //  2. u*A.y + v*B.y = C.y
    //  3. u*A.z + v*B.z = C.z

    //  Using 1, 2, find u and v
    //  | A.x  B.x | |u|   | C.x |
    //  | A.y  B.y | |v| = | C.y |


    
    rg_REAL* param = new rg_REAL[2];

    rg_REAL denominator1 = ( A.getX()*B.getY() ) - ( B.getX()*A.getY() );
    rg_REAL denominator2 = ( A.getX()*B.getZ() ) - ( B.getX()*A.getZ() );
    rg_REAL denominator3 = ( A.getY()*B.getZ() ) - ( B.getY()*A.getZ() );
    rg_REAL u = 0.0;
    rg_REAL v = 0.0;

    //  if ( rg_NZERO( denominator1 )  )
    if ( rg_NZERO( denominator1, rg_SYSTEM_RES )  )
    {
        u = ( (C.getX()*B.getY()) - (B.getX()*C.getY()) )
            / denominator1;
        v = ( (A.getX()*C.getY()) - (C.getX()*A.getY()) )
            / denominator1;

        if ( rg_EQ( (u*A.getZ())+(v*B.getZ()), C.getZ() ) )
        {
            param[0] = u;
            param[1] = v;

            return param;
        }
    }
//    else if ( rg_NZERO( denominator2 )  )
    else if ( rg_NZERO( denominator2, rg_SYSTEM_RES )  )
    {
        u = ( (C.getX()*B.getZ()) - (B.getX()*C.getZ()) )
            / denominator2;
        v = ( (A.getX()*C.getZ()) - (C.getX()*A.getZ()) )
            / denominator2;

        if ( rg_EQ( (u*A.getY())+(v*B.getY()), C.getY() ) )
        {
            param[0] = u;
            param[1] = v;

            return param;
        }
    }
//    else if ( rg_NZERO( denominator3 )  )
    else if ( rg_NZERO( denominator3, rg_SYSTEM_RES )  )
    {
        u = ( (C.getY()*B.getZ()) - (B.getY()*C.getZ()) )
            / denominator3;
        v = ( (A.getY()*C.getZ()) - (C.getY()*A.getZ()) )
            / denominator3;

        if ( rg_EQ( (u*A.getX())+(v*B.getX()), C.getX() ) )
        {
            param[0] = u;
            param[1] = v;

            return param;
        }
    }
    else ;

    delete [] param;
    
    return rg_NULL;    
}


 



