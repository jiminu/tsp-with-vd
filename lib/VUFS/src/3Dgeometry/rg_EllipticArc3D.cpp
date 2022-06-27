/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_EllipticArc3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_EllipticArc3D 
//           which define ellipse and its property. 
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Taeboom Jang
//    START DATE  : 2001.  1.  29. 
//
//    HISTORY     :
//
//    
//            Copyright (c) CAD/CAM Lab.    	  	
//
/////////////////////////////////////////////////////////////////////
#include "rg_EllipticArc3D.h"
#include <math.h>
#include "rg_GeoFunc.h"

rg_EllipticArc3D::rg_EllipticArc3D()
: centerPt(0.0,0.0,0.0), localXAxis(1.0,0.0,0.0), localYAxis(0.0,1.0,0.0), startAngle(0.0), endAngle(0.0)
{
}

rg_EllipticArc3D::rg_EllipticArc3D(const rg_EllipticArc3D &temp)
: centerPt(temp.centerPt), localXAxis(temp.localXAxis), localYAxis(temp.localYAxis), startAngle(temp.startAngle), endAngle(temp.endAngle)
{
}

rg_EllipticArc3D::~rg_EllipticArc3D()
{
}

	
//get functions
rg_Point3D  rg_EllipticArc3D::getCenterPt()   const
{
    return centerPt;
}

rg_Point3D  rg_EllipticArc3D::getLocalXAxis() const
{
    return localXAxis;
}

rg_Point3D  rg_EllipticArc3D::getLocalYAxis() const
{
    return localYAxis;
}

rg_REAL     rg_EllipticArc3D::getStartAngle() const
{
    return startAngle;
}

rg_REAL     rg_EllipticArc3D::getEndAngle()   const
{
    return endAngle;
}
	
//set functions
void rg_EllipticArc3D::setCenterPt(const rg_Point3D &tCenterPt)
{
    centerPt=tCenterPt;
}

void rg_EllipticArc3D::setLocalXAxis(const rg_Point3D &tLocalXAxis)
{
    localXAxis=tLocalXAxis;
}

void rg_EllipticArc3D::setLocalYAxis(const rg_Point3D &tLocalYAxis)
{
    localYAxis=tLocalYAxis;
}

void rg_EllipticArc3D::setStartAngle(const rg_REAL &tStartAngle)
{
    startAngle=tStartAngle;
}

void rg_EllipticArc3D::setEndAngle(const rg_REAL &tEndAngle)
{
    endAngle=tEndAngle;
}

void rg_EllipticArc3D::setEllipticArc3D(const rg_EllipticArc3D &temp)
{
    centerPt=temp.centerPt;
    localXAxis=temp.localXAxis;
    localYAxis=temp.localYAxis;
    startAngle=temp.startAngle;
    endAngle=temp.endAngle;
}

//evaluate
rg_Point3D rg_EllipticArc3D::evaluatePt(const rg_REAL& angle) const 
{
    const rg_REAL normAngle=rg_GeoFunc::getNormalizedAngle(angle);

    rg_Point3D output(centerPt);

    const rg_REAL cosAngle=cos(normAngle);
    const rg_REAL sinAngle=sin(normAngle);

    output=output+localXAxis*cosAngle;
    output=output+localYAxis*sinAngle;

    return output;
}

rg_Point3D rg_EllipticArc3D::evaluateDerivative(const rg_REAL& angle) const
{
    const rg_REAL normAngle=rg_GeoFunc::getNormalizedAngle(angle);

    rg_Point3D output(0.0, 0.0, 0.0);
    const rg_REAL cosAngle=cos(normAngle);
    const rg_REAL sinAngle=sin(normAngle);

    output+=-sinAngle*localXAxis;
    output+= cosAngle*localYAxis;

    return output;
}
// ETC
rg_FLAG    rg_EllipticArc3D::isOn(const rg_Point3D& pt) const 
{
    rg_Point3D* foci=getFoci();
    rg_Point3D  majorAxis=getMajorAxis();
    rg_REAL     lengthOfMajorAxis=2.0*majorAxis.magnitude();
    rg_REAL     distance=   foci[0].distance(pt) 
                       + foci[1].distance(pt);
    rg_FLAG     output=rg_EQ(distance, lengthOfMajorAxis);
    delete[] foci;

    return output;
}


rg_Point3D* rg_EllipticArc3D::getFoci() const
{
    rg_Point3D  majorAxis=this->getMajorAxis();
    rg_Point3D  minorAxis=this->getMinorAxis();
    
    rg_REAL     halfLengthOfMajorAxis=majorAxis.magnitude();
    rg_REAL     halfLengthOfMinorAxis=minorAxis.magnitude();
    rg_Point3D  unitMajorAxis=majorAxis.getUnitVector();

    rg_REAL     halfDistanceBtFoci= sqrt (  pow(halfLengthOfMajorAxis,2.0) 
                                       - pow(halfLengthOfMinorAxis,2.0));
    
    rg_Point3D* output=new rg_Point3D[2];

    output[0]=centerPt+  halfDistanceBtFoci* unitMajorAxis;
    output[1]=centerPt-  halfDistanceBtFoci* unitMajorAxis;

    return output;
}

rg_Point3D  rg_EllipticArc3D::getMajorAxis() const
{
    rg_REAL lengthOfLocalX=localXAxis.magnitude();
    rg_REAL lengthOfLocalY=localYAxis.magnitude();

    rg_Point3D output;
    output = ( rg_GE(lengthOfLocalX, lengthOfLocalY) ) ? localXAxis:localYAxis;

    return output;
}


rg_Point3D  rg_EllipticArc3D::getMinorAxis() const
{
    rg_REAL lengthOfLocalX=localXAxis.magnitude();
    rg_REAL lengthOfLocalY=localYAxis.magnitude();

    rg_Point3D output;
    output = ( rg_LE(lengthOfLocalX, lengthOfLocalY) ) ? localXAxis:localYAxis;

    return output;
}

rg_REAL rg_EllipticArc3D::getAngle(const rg_Point3D& pt) const
{
    rg_REAL localXScale=localXAxis.magnitude();
    rg_REAL localYScale=localYAxis.magnitude();

    rg_Point3D translatedPoint=pt-centerPt;

    rg_REAL localXCoord= translatedPoint%(localXAxis.getUnitVector());
    rg_REAL localYCoord= translatedPoint%(localYAxis.getUnitVector());

    rg_REAL noramlizedX=localXCoord/localXScale;
    rg_REAL noramlizedY=localYCoord/localYScale;

    rg_REAL output=atan2(noramlizedY,noramlizedX);

    return output;
}





//operators	
rg_EllipticArc3D& rg_EllipticArc3D::operator=(const rg_EllipticArc3D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }

    setEllipticArc3D(temp);
    return *this;
}


