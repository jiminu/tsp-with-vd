#include "rg_Outline3D.h"

rg_Outline3D::rg_Outline3D()
: rg_NURBSplineCurve3D()
, passingPts(rg_NULL), numOfPassingPts(0)
{
}

rg_Outline3D::rg_Outline3D(const rg_Outline3D& temp)
: rg_NURBSplineCurve3D(temp)
, numOfPassingPts(temp.numOfPassingPts)
{
    passingPts=new rg_Point3D[numOfPassingPts];

    for( rg_INT i=0; i < numOfPassingPts; i++ )
    {
        passingPts[i]=temp.passingPts[i];
    }
}


rg_Outline3D::~rg_Outline3D()
{
    if ( passingPts != rg_NULL )
    {
        delete[] passingPts;
    }
}

// get function
rg_Point3D* rg_Outline3D::getPassingPts() const
{
    if ( passingPts == rg_NULL )
    {
        return rg_NULL;
    }

    rg_Point3D* output=new rg_Point3D[numOfPassingPts];
    for( rg_INT i=0; i < numOfPassingPts; i++ )
    {
        output[i]=passingPts[numOfPassingPts];
    }
    return output;
}

rg_Point3D  rg_Outline3D::getPassingPt(const rg_INT& index) const
{
    return passingPts[index];
}

rg_INT    rg_Outline3D::getNumOfPassingPts() const
{
    return numOfPassingPts;
}


// set function
void   rg_Outline3D::setPassingPts( rg_Point3D* tPassingPts,
                               const rg_INT&  tNumOfPassingPts)
{
    if ( passingPts != rg_NULL )
    {
        delete[] passingPts;
    }

    numOfPassingPts=tNumOfPassingPts;
    passingPts=new rg_Point3D[numOfPassingPts];

    for( rg_INT i=0; i < numOfPassingPts; i++ )
    {
        passingPts[i]=tPassingPts[i];
    }
}

void   rg_Outline3D::setPassingPt( const rg_INT& index,
                              const rg_Point3D& pt )
{
    passingPts[index]=pt;
}

// make curve
void rg_Outline3D::fitCurve( const rg_INT& tOrder,
                        const rg_Point3D& startTangent,
                        const rg_Point3D& endTangent)
{
    rg_NUBSplineCurve3D::interpolateWithEndCondition( numOfPassingPts,
                                                      passingPts,
                                                      startTangent,
                                                      endTangent);
    rg_INT numOfCtrlPts=rg_BSplineCurve3D::getNumOfCtrlPts();
    
    rg_REAL* newWeights=new rg_REAL[numOfCtrlPts];

    for( rg_INT i=0 ; i < numOfCtrlPts; i++ )
    {
        newWeights[i]=1.0;
    }
    rg_NURBSplineCurve3D::setWeightVector(newWeights);
    delete[] newWeights;
}
             

// operator overloading
rg_Outline3D& rg_Outline3D::operator=(const rg_Outline3D& temp)
{
    if ( &temp == this )
    {
        return *this;
    }

    rg_NURBSplineCurve3D::operator=(temp);
    if ( passingPts != rg_NULL )
    {
        delete[] passingPts;
    }
    passingPts=new rg_Point3D[numOfPassingPts];

    for( rg_INT i=0; i < numOfPassingPts; i++ )
    {
        passingPts[i]=temp.passingPts[i];
    }
    return *this;
}    



