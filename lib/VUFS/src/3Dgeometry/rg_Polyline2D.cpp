#include "rg_Const.h"
#include "rg_Polyline2D.h"
#include "rg_Point2D.h"
#include "rg_Line2D.h"
#include "rg_RelativeOp.h"

rg_Polyline2D::rg_Polyline2D()
: numOfPts(0), pts(rg_NULL)
{
}

rg_Polyline2D::rg_Polyline2D(const rg_INT            &tNumOfPts,
                             const rg_Point2D* const  tPts)
:numOfPts(tNumOfPts)
{
    if ( numOfPts > 0 )
    {
        pts=new rg_Point2D[numOfPts];
        for(rg_INT i=0; i < numOfPts; i++ )
        {
            pts[i]=tPts[i];
        }
    }
    else
    {
        pts=rg_NULL;
    }
}

rg_Polyline2D::rg_Polyline2D(const rg_Polyline2D&     temp)
:numOfPts(temp.numOfPts)
{
    if ( numOfPts > 0 )
    {
        pts=new rg_Point2D[numOfPts];
        for(rg_INT i=0; i < numOfPts; i++ )
        {
            pts[i]=temp.pts[i];
        }
    }
    else
    {
        pts=rg_NULL;
    }
}

rg_Polyline2D::~rg_Polyline2D()
{
    removeAll();
}

void rg_Polyline2D::removeAll()
{
    if ( numOfPts > 0 )
    {
        delete[] pts;
    }
    numOfPts=0;
    pts=rg_NULL;
}

rg_INT        rg_Polyline2D::getNumOfPts() const
{
    return numOfPts;
}

rg_Point2D*   rg_Polyline2D::getPts()  const
{
    if ( numOfPts < 1 )
    {
        return rg_NULL;
    }

    rg_Point2D* output=new rg_Point2D[numOfPts];
    for( rg_INT i=0; i < numOfPts; i++)
    {
        output[i]=pts[i];
    }

    return output;
}
rg_Point2D    rg_Polyline2D::getPt(const rg_INT &index)  const
{
    return pts[index];
}
rg_Point2D*   rg_Polyline2D::accessPts() const
{
    return pts;
}

void          rg_Polyline2D::setNumOfPts( const rg_INT& tNumOfPts)
{
    removeAll();
    numOfPts=tNumOfPts;

    if ( numOfPts < 1 )
    {
        pts=rg_NULL;
        return;
    }

    pts=new rg_Point2D[numOfPts];

}
        
void          rg_Polyline2D::setPts(const rg_Point2D* const tPts)
{
    if ( numOfPts < 1 || tPts == rg_NULL)
    {
        return;
    }

    for( rg_INT i=0; i < numOfPts; i++ )
    {
        pts[i]=tPts[i];
    }
}


void          rg_Polyline2D::setPt(const rg_INT     &index,
                                   const rg_Point2D &tPt)
{
    if ( numOfPts < 1 )
    {
        return;
    }
    pts[index]=tPt;

}

void          rg_Polyline2D::setAll(const rg_Polyline2D& temp)
{
    removeAll();
    
    numOfPts=temp.numOfPts;

    if ( numOfPts < 1 )
    {
        pts=rg_NULL;
        return;
    }
    pts=new rg_Point2D[numOfPts];
    for( rg_INT i=0; i < numOfPts; i++ )
    {
        pts[i]=temp.pts[i];
    }
}


rg_Polyline2D& rg_Polyline2D::operator=(const rg_Polyline2D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }
    setAll(temp);

    return *this;
}
    
rg_FLAG   rg_Polyline2D::isNull() const
{
    if ( numOfPts < 1 )
    {
        return rg_TRUE;
    }
    else
    {
        return rg_FALSE;
    }
}

rg_REAL   rg_Polyline2D::evaluateDistance(const rg_Point2D& pt) const
{
    if ( isNull() )
    {
        return rg_MAX_REAL;
    }
    if ( numOfPts == 1 )
    {
        rg_REAL  distance=pts[0].distance(pt);
        return distance;
    }

    rg_REAL distance=rg_MAX_REAL;

    for( rg_INT i=1; i < numOfPts; i++ )
    {
        rg_Point2D prevPt=pts[i-1];
        rg_Point2D nextPt=pts[i];
        rg_Line2D currentLine(prevPt,nextPt);

        rg_REAL currentDistance=currentLine.getDistance(pt);

        if ( rg_LT( currentDistance, distance ) )
        {
            distance=currentDistance;
        }
    }
    return distance;
}

rg_FLAG   rg_Polyline2D::isOverlapped(const rg_Point2D& pt, 
                                      const rg_REAL& tolerance ) const
{
    if ( isNull() )
    {
        return rg_FALSE;
    }
/*
    if ( numOfPts == 1 )
    {
        rg_REAL  distance=pts[0].distance(pt);
        return ( rg_LE ( distance, tolerance ) );
    }

    rg_REAL distance=rg_MAX_REAL;

    for( rg_INT i=1; i < numOfPts; i++ )
    {
        rg_Point2D prevPt=pts[i-1];
        rg_Point2D nextPt=pts[i];
        rg_Line2D currentLine(prevPt,nextPt);

        rg_REAL currentDistance=currentLine.getDistance(pt);

        if ( rg_LT( currentDistance, distance ) )
        {
            distance=currentDistance;
        }
    }
*/
    rg_REAL distance= evaluateDistance(pt);
    return rg_LE( distance, tolerance) ;
}

        
