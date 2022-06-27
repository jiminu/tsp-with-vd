//********************************************************************
//
//	  FILENAME    : rg_Plane3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of class rg_Point3D.      
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 8 Jul. 1997    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#include "rg_Plane3D.h"
#include "rg_RelativeOp.h"

//  Constructor & Destructor.--------------------------------------------------
rg_Plane3D::rg_Plane3D()
{
}

rg_Plane3D::rg_Plane3D(const rg_Point3D& pVector, const rg_Point3D &uVector, const rg_Point3D &vVector)
{
    posVec = pVector;
    uVec   = uVector;
    vVec   = vVector;
}

rg_Plane3D::rg_Plane3D(const rg_Plane3D &plane)
{
    posVec = plane.posVec;
    uVec   = plane.uVec;
    vVec   = plane.vVec;
}

rg_Plane3D::~rg_Plane3D()
{
}

//  Access Functions.----------------------------------------------------------
rg_Point3D rg_Plane3D::getPosVector() const
{
    return posVec;
}

rg_Point3D rg_Plane3D::getUVector() const
{
    return uVec;
}

rg_Point3D rg_Plane3D::getVVector() const
{
    return vVec;
}

void rg_Plane3D::setPosVector(const rg_Point3D &pVector)
{
    posVec = pVector;
}

void rg_Plane3D::setUVector(const rg_Point3D &uVector)
{
    uVec = uVector;
}

void rg_Plane3D::setVVector(const rg_Point3D &vVector)
{
    vVec = vVector;
}

void rg_Plane3D::setPlane(const rg_Point3D &pVector, const rg_Point3D &uVector, const rg_Point3D &vVector)
{
    posVec = pVector;
    uVec  = uVector;
    vVec  = vVector;
}

void rg_Plane3D::setPlane(const rg_Point3D* const cornerPt)
{
    rg_INDEX nearOrigin  = -1;
    rg_INDEX farOrigin   = -1;
    rg_REAL  minDistance = 1000000;
    rg_REAL  maxDistance = -1;	

    rg_INDEX i = 0; 
    for (i=0; i<4; i++)
    {
        rg_REAL dist = cornerPt[i].distance(rg_Point3D(0., 0., 0.));
        if ( rg_LE(dist, minDistance ) )
        {
            minDistance = dist;
            nearOrigin  = i;
        }
        else if ( rg_GE(dist, maxDistance) )
        {
            maxDistance = dist;
            farOrigin   = i;
        }
    }

    rg_Point3D pt[2];
    rg_INDEX k =0;
    for (i=0; i<4; i++)
    {
        if (i != nearOrigin && i != farOrigin)
        {
            pt[k] = cornerPt[i];
            k++;
        }
    }

    if ( (pt[0] + pt[1] - cornerPt[nearOrigin]) != cornerPt[farOrigin])  
        return;
    
    posVec = cornerPt[nearOrigin];
    uVec   = pt[0] - posVec;
    vVec   = pt[1] - posVec;
}

void rg_Plane3D::setPlane(const rg_Plane3D &plane)
{
    posVec = plane.posVec;
    uVec  = plane.uVec;
    vVec  = plane.vVec;
}

//  Evaluate rg_Point3D.------------------------------------------------------------
rg_Point3D rg_Plane3D::evaluatePt(const rg_PARAMETER &u, const rg_PARAMETER &v) const
{
    return posVec + (u*uVec) + (v*vVec);
}

//  Operations.----------------------------------------------------------------
rg_FLAG rg_Plane3D::isParallelTo(const rg_Line3D &line) const
{
    rg_Point3D vector     = line.getEP() - line.getSP();
    rg_Point3D crsProduct = uVec*vVec;

    if ( crsProduct%vector )
        return rg_FALSE;
    else
        return rg_TRUE;
    
}

rg_FLAG rg_Plane3D::isParallelTo(const rg_Plane3D &plane) const
{
    rg_Point3D crsProduct = uVec*vVec;

    if ( (crsProduct%plane.uVec) && (crsProduct%plane.vVec) )
        return rg_FALSE;
    else
        return rg_TRUE;
}

rg_FLAG rg_Plane3D::isPerpendicularTo(const rg_Line3D &line) const
{
    rg_Point3D vector     = line.getEP() - line.getSP();

    if ( (uVec%vector) && (vVec%vector) )
        return rg_FALSE;
    else
        return rg_TRUE;
}

rg_FLAG rg_Plane3D::isPerpendicularTo(const rg_Plane3D &plane) const
{
    if ( (uVec * vVec)%(plane.uVec * plane.vVec) )
        return rg_FALSE;
    else 
        return rg_TRUE;
}

rg_FLAG rg_Plane3D::isThisPointOnPlane(const rg_Point3D &pt) const
{
    //  A + uB + vC = pt
    //  -> uB + vC = pt - A
    rg_REAL* param; // param[0] = u, param[1] = v

    param = solveSimultaneousEq(uVec, vVec, pt-posVec);

    if ( rg_GE(param[0], 0.0) && rg_LE(param[0], 1.0)
         && rg_GE(param[1], 0.0) && rg_LE(param[1], 1.0) )
    {
        delete [] param;

        return rg_TRUE;
    }
    else
    {
        if (param != rg_NULL)
            delete [] param;

        return rg_FALSE;
    }
}

//  Intersect. ----------------------------------------------------------------

rg_sListByPtr* rg_Plane3D::intersectLine(const rg_Line3D &line, rg_INT &numOfintersectPt) const
{
    //  rg_Line3D has two points which represent start point and end point of line.
    //  Parametric form of rg_Line3D form
    //      line(t) = (1-t)*sp + t*ep
    //              = sp + t*(ep - sp)
    //
    //  Intersection of plane and line
    //      plane(u, v) = line(t)
    //      A + uB + vC = sp + t*(ep - sp)
    //

    rg_sListByPtr* intersectPt = new rg_sListByPtr;
    numOfintersectPt = 0;

    rg_Point3D crsPrct = uVec*vVec;
    rg_Point3D sp      = line.getSP();
    rg_Point3D ep      = line.getEP();

    if ( rg_ZERO( crsPrct%(ep - sp) ) )
    {
        rg_FLAG isSpOnPlane = isThisPointOnPlane(sp);
        rg_FLAG isEpOnPlane = isThisPointOnPlane(ep);

        if ( isSpOnPlane && isEpOnPlane )
        {
            //  This line is on plane.
            rg_Point3D* intersect1 = new rg_Point3D( sp );
            intersectPt->append(intersect1);

            rg_Point3D* intersect2 = new rg_Point3D( ep );
            intersectPt->append(intersect2);

            numOfintersectPt = 2;

            return intersectPt;
        }
        else 
        {
            rg_Line3D comparedLine[4];
            comparedLine[0].setLine(posVec, posVec + uVec);
            comparedLine[1].setLine(posVec, posVec + vVec);
            comparedLine[2].setLine(posVec + uVec, posVec + uVec + vVec);
            comparedLine[3].setLine(posVec + vVec, posVec + uVec + vVec);

            for (rg_INDEX i=0; i<4; i++)
            {
                INTERSECT_TYPE type;
                rg_Point3D* intersect1 = line.intersectLine(comparedLine[i], type);
                if ( intersect1 != rg_NULL )
                {
                    rg_Point3D* intersect2;
                    if ( isSpOnPlane ) 
                    {
                        intersect2  = new rg_Point3D( sp );
                        if ( (*intersect1) == (*intersect2) )
                        {
                            intersectPt->append(intersect1);

                            numOfintersectPt = 1;
                        }
                        else
                        {
                            intersectPt->append(intersect1);
                            intersectPt->append(intersect2);

                            numOfintersectPt = 2;
                        }
    
                        return intersectPt;
                    }
                    else if ( isEpOnPlane ) 
                    {
                        intersect2  = new rg_Point3D( ep );
                        if ( (*intersect1) == (*intersect2) )
                        {
                            intersectPt->append(intersect1);

                            numOfintersectPt = 1;
                        }
                        else
                        {
                            intersectPt->append(intersect1);
                            intersectPt->append(intersect2);

                            numOfintersectPt = 2;
                        }

                        return intersectPt;
                    }
                    else
                    {
                        intersectPt->append(intersect1);

                        numOfintersectPt++;
                    }
                }   // end of if( intersect1 != rg_NULL ) 
            }   // end of for(    
            return intersectPt;
        }
    }   // if ( rg_ZERO( crsPrct%(ep - sp) ) )
    else 
    {
        rg_REAL param = ( (sp%crsPrct) - (posVec%crsPrct) )
                     / ( (ep - sp)%crsPrct );

        if ( rg_GE(param, 0.0) && rg_LE(param, 1.) )
        {
            rg_Point3D* intersect = new rg_Point3D( sp + (param*(ep - sp)) );

            intersectPt->append(intersect);

            numOfintersectPt = 1;

            return intersectPt;
        }
        else
        {
            numOfintersectPt = 0;

            return intersectPt = rg_NULL;
        }
    }
}

rg_sListByPtr* rg_Plane3D::intersectPlane(const rg_Plane3D &plane, rg_INT &numOfintersectPt) const
{   
    rg_sListByPtr* intersectPt = new rg_sListByPtr;
    numOfintersectPt = 0;

    rg_Point3D  crsPrct = plane.uVec*plane.vVec;

    if ( rg_ZERO(uVec%crsPrct) && rg_ZERO(vVec%crsPrct) )
    {
    //  Two planes are parallel.
    //  This part is not made.
        return rg_NULL;
    }
    else if ( rg_ZERO(uVec%crsPrct) && rg_NZERO(vVec%crsPrct) )
    {
    //  A vector of u parametric direction is parallel to either 
    //  one of s parametric direction or the other of t parametric 
    //  direction
        rg_REAL vParam = ( (plane.posVec%crsPrct) - (posVec%crsPrct) )
                      / (vVec%crsPrct);

        rg_Line3D line( evaluatePt(0, vParam), evaluatePt(1, vParam) );

        intersectPt = plane.intersectLine(line, numOfintersectPt);

        return intersectPt;
    }               
    else if ( rg_NZERO(uVec%crsPrct) && rg_ZERO(vVec%crsPrct) )
    {     
    //  A vector of v parametric direction is parallel to either 
    //  one of s parametric direction or the other of t parametric 
    //  direction
        rg_REAL uParam = ( (plane.posVec%crsPrct) - (posVec%crsPrct) )
                      / (uVec%crsPrct);

        rg_Line3D line( evaluatePt(uParam, 0), evaluatePt(uParam, 1) );

        intersectPt = plane.intersectLine(line, numOfintersectPt);

        return intersectPt;
    }
    else
    {
        //  A + uB + vC = D + sE + tF
        rg_REAL param[8][2];

        rg_Point3D crsPrct1 = uVec*vVec;
        rg_Point3D crsPrct2 = plane.uVec*plane.vVec;

        //  v = 0, u?
        param[0][0] 
            = (plane.posVec%crsPrct2 - posVec%crsPrct2) / (uVec%crsPrct2);
        param[0][1] = 0.0;

        //  v = 1, u?
        param[1][0] 
            = (plane.posVec%crsPrct2 - posVec%crsPrct2 - vVec%crsPrct2) 
              / (uVec%crsPrct2);
        param[1][1] = 1.0;

        //  u = 0, v?
        param[2][0] = 0.0;
        param[2][1]
            = (plane.posVec%crsPrct2 - posVec%crsPrct2) / (vVec%crsPrct2);

        //  u = 1, v?
        param[3][0] = 1.0;
        param[3][1] 
            = (plane.posVec%crsPrct2 - posVec%crsPrct2 - uVec%crsPrct2)
              / (vVec%crsPrct2);

        //  t = 0, s?
        param[4][0] 
            = (posVec%crsPrct1 - plane.posVec%crsPrct1) / (plane.uVec%crsPrct1);
        param[4][1] = 0.0;

        //  t = 1, s?
        param[5][0] 
            = (posVec%crsPrct1 - plane.posVec%crsPrct1 - plane.vVec%crsPrct1)
              / (plane.uVec%crsPrct1);
        param[5][1] = 1.0;

        //  s = 0, t?
        param[6][0] = 0.0;
        param[6][1]
            = (posVec%crsPrct1 - plane.posVec%crsPrct1) / (plane.vVec%crsPrct1);
        
        //  s = 0, t?
        param[7][0] = 1.0;
        param[7][1]
            = (posVec%crsPrct1 - plane.posVec%crsPrct1 - plane.uVec%crsPrct1)
              / (plane.vVec%crsPrct1);

        //  Find the parameters,s and t, of the compared plane
        //  , using the obtained u and v. 
        rg_INDEX i = 0;
	for (i=0; i<4; i++)
        {
            //  check whether both u and v are valid or not.
            if ( rg_GE(param[i][0], 0.0) && rg_LE(param[i][0], 1.0) 
                 && rg_GE(param[i][1], 0.0) && rg_LE(param[i][1], 1.0) )               
            {
                rg_Point3D leftPart = posVec + (param[i][0]*uVec)
                                 + (param[i][1]*vVec) - plane.posVec;
            
                //  stValue[0] = s, stValue[1] = t.
                rg_REAL* stValue = solveSimultaneousEq(plane.uVec, plane.vVec, leftPart);
                
                //  check whether both s and t are valid or not.
                if ( (stValue == rg_NULL) ||
                     !( rg_GE(stValue[0], 0.0) && rg_LE(stValue[0], 1.0) 
                        && rg_GE(stValue[1], 0.0) && rg_LE(stValue[1], 1.0) ) )             
                {
                    continue;
                }
                
                rg_Point3D intersect( plane.posVec + (stValue[0]*plane.uVec) 
                                 + (stValue[1]*plane.vVec) );
                
                rg_Point3D* prevPt = (rg_Point3D*)intersectPt->getFirstNode();
                if ( (prevPt == rg_NULL) || 
                     ( (*prevPt != intersect) && ( numOfintersectPt<2 ) ) )
                {
                    rg_Point3D* point = new rg_Point3D( intersect );
                    intersectPt->append(point);

                    numOfintersectPt++;
                }

                delete [] stValue;            
            }
        }

        //  Find the parameters,u and v, of the plane
        //  , using the obtained s and t. 
        for (i=4; i<8; i++)
        {
            //  check whether both s and t are valid or not.
            if ( rg_GE(param[i][0], 0.0) && rg_LE(param[i][0], 1.0) 
                 && rg_GE(param[i][1], 0.0) && rg_LE(param[i][1], 1.0) )
            {
                rg_Point3D rightPart = posVec + (param[i][0]*uVec)
                                 + (param[i][1]*vVec) - plane.posVec;
            
                //  stValue[0] = u, stValue[1] = v.
                rg_REAL* stValue = solveSimultaneousEq(uVec, vVec, rightPart);
                
                //  check whether both u and v are valid or not.
                if ( (stValue == rg_NULL) ||
                     !( rg_GE(stValue[0], 0.0) && rg_LE(stValue[0], 1.0) 
                       && rg_GE(stValue[1], 0.0) && rg_LE(stValue[1], 1.0) ) )             
                {
                    continue;
                }

                rg_Point3D intersect( posVec + (stValue[0]*uVec) + (stValue[1]*vVec) );
                
                rg_Point3D* prevPt = (rg_Point3D*)intersectPt->getFirstNode();
                if ( (prevPt == rg_NULL) || 
                     ( (*prevPt != intersect) && ( numOfintersectPt<2 ) ) )
                {
                    rg_Point3D* point = new rg_Point3D( intersect );
                    intersectPt->append(point);

                    numOfintersectPt++;
                }

                delete [] stValue;            
            }
        }

        return intersectPt;
    }
}

//  Operator Overloading.------------------------------------------------------
rg_Plane3D& rg_Plane3D::operator =(const rg_Plane3D &plane)
{
    posVec = plane.posVec;
    uVec  = plane.uVec;
    vVec  = plane.vVec;

    return *this;
}

rg_FLAG  rg_Plane3D::operator==(const rg_Plane3D &plane) const
{
    rg_Point3D ptsOnThisPlane[4]     
        = { posVec, (posVec + uVec),
            (posVec + uVec + vVec), (posVec + vVec) };
    rg_Point3D ptsOnComparedPlane[4] 
        = { plane.posVec, (plane.posVec + plane.uVec),
            (plane.posVec + plane.uVec + plane.vVec), 
            (plane.posVec + plane.vVec) };

    for (rg_INDEX i = 0; i<4; i++)
    {
        if (ptsOnThisPlane[0] == ptsOnComparedPlane[i])
        {
            if (ptsOnThisPlane[2] == ptsOnComparedPlane[(i+2)%4] )
            {
                if (   (ptsOnThisPlane[1] == ptsOnComparedPlane[ (i+1)%4 ]
                       && ptsOnThisPlane[3] == ptsOnComparedPlane[ (i+3)%4 ])
                    || (ptsOnThisPlane[1] == ptsOnComparedPlane[ (i+3)%4 ]
                       && ptsOnThisPlane[3] == ptsOnComparedPlane[ (i+1)%4 ]) )
                {
                    return rg_TRUE;
                }
            }
        }
    }

    return rg_FALSE;
}



