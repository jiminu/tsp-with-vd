#include "rg_Const.h"
#include "rg_Polyline3D.h"
#include "rg_Point3D.h"

rg_Polyline3D::rg_Polyline3D()
: numOfPts(0), pts(rg_NULL)
{
}

rg_Polyline3D::rg_Polyline3D(const rg_INT            &tNumOfPts,
                             const rg_Point3D* const  tPts)
:numOfPts(tNumOfPts)
{
    if ( numOfPts > 0 )
    {
        pts=new rg_Point3D[numOfPts];
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

rg_Polyline3D::rg_Polyline3D(const rg_Polyline3D&     temp)
:numOfPts(temp.numOfPts)
{
    if ( numOfPts > 0 )
    {
        pts=new rg_Point3D[numOfPts];
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

rg_Polyline3D::~rg_Polyline3D()
{
    removeAll();
}

void rg_Polyline3D::removeAll()
{
    if ( numOfPts > 0 )
    {
        delete[] pts;
    }
    numOfPts=0;
    pts=rg_NULL;
}

rg_INT        rg_Polyline3D::getNumOfPts() const
{
    return numOfPts;
}

rg_Point3D*   rg_Polyline3D::getPts()  const
{
    if ( numOfPts < 1 )
    {
        return rg_NULL;
    }

    rg_Point3D* output=new rg_Point3D[numOfPts];
    for( rg_INT i=0; i < numOfPts; i++)
    {
        output[i]=pts[i];
    }

    return output;
}

rg_Point3D*   rg_Polyline3D::accessPts() const
{
    return pts;
}

void          rg_Polyline3D::setNumOfPts( const rg_INT& tNumOfPts)
{
    removeAll();
    numOfPts=tNumOfPts;

    if ( numOfPts < 1 )
    {
        pts=rg_NULL;
        return;
    }

    pts=new rg_Point3D[numOfPts];

}
        
void          rg_Polyline3D::setPts(const rg_Point3D* const tPts)
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


void          rg_Polyline3D::setPt(const rg_INT     &index,
                                   const rg_Point3D &tPt)
{
    if ( numOfPts < 1 )
    {
        return;
    }
    pts[index]=tPt;

}

void          rg_Polyline3D::setAll(const rg_Polyline3D& temp)
{
    removeAll();
    
    numOfPts=temp.numOfPts;

    if ( numOfPts < 1 )
    {
        pts=rg_NULL;
        return;
    }
    pts=new rg_Point3D[numOfPts];
    for( rg_INT i=0; i < numOfPts; i++ )
    {
        pts[i]=temp.pts[i];
    }
}


rg_Polyline3D& rg_Polyline3D::operator=(const rg_Polyline3D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }
    setAll(temp);

    return *this;
}
    

