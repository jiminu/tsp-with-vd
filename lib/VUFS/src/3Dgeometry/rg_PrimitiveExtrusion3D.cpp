//******************************************************************************
//
//	  FILENAME    : rg_PrimitiveExtrusion3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_PrimitiveExtrusion3D 
//           which will be the base of all kinds of rg_Extrusion3D
//                          
//	  CLASS NAME  : rg_PrimitiveExtrusion3D
//
//    BASE CLASS  : rg_Geometry
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 8.22
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea    	  	
//
//******************************************************************************
#include "rg_PrimitiveExtrusion3D.h"

rg_PrimitiveExtrusion3D::rg_PrimitiveExtrusion3D()
:rg_Geometry()
{
    height=0.0;
}

rg_PrimitiveExtrusion3D::rg_PrimitiveExtrusion3D(const GeometryType& tType,const unsigned rg_INT& tID)
:rg_Geometry(tType,tID)
{
    height=0.0;
}

rg_PrimitiveExtrusion3D::rg_PrimitiveExtrusion3D(const rg_PrimitiveExtrusion3D& temp)
:rg_Geometry(temp)
{
    height=temp.height;
    outerLoop=temp.outerLoop;
    listOfInnerLoop=temp.listOfInnerLoop;
}

rg_PrimitiveExtrusion3D::~rg_PrimitiveExtrusion3D()
{
}

// get function
rg_REAL  rg_PrimitiveExtrusion3D::getHeight() const
{
    return height;
}

rg_INT   rg_PrimitiveExtrusion3D::getNumOfPointInOuterLoop() const
{
    return outerLoop.getSize();
}

rg_INT   rg_PrimitiveExtrusion3D::getNumOfPointInInnerLoop(const rg_INT& index) const
{
    return listOfInnerLoop[index].getSize();
}

rg_INT   rg_PrimitiveExtrusion3D::getNumOfInnerLoop() const
{
    return listOfInnerLoop.getSize();
}

rg_Array< rg_Point3D > rg_PrimitiveExtrusion3D::getOuterLoop() const
{
    return outerLoop;
}

rg_Array< rg_Point3D > rg_PrimitiveExtrusion3D::getInnerLoop( const rg_INT& index) const
{
    return listOfInnerLoop[index];
}

rg_Point3D rg_PrimitiveExtrusion3D::getPointInOuterLoop(const rg_INT& index) const
{
    return outerLoop[index];
}

rg_Array<rg_Point3D> rg_PrimitiveExtrusion3D::getOuterLoopInUpperPolygon() const
{
    rg_Array<rg_Point3D> output;
    rg_INT numOfPoints=outerLoop.getSize();
    output.setSize(numOfPoints);
    rg_REAL  halfHeight=height/2.0;

    for( rg_INT i=0; i < numOfPoints; i++ )
    {
        output[i]=outerLoop[i]+halfHeight*rg_Point3D(0.0,0.0,1.0);
    }
    return output;
}

rg_Array<rg_Point3D> rg_PrimitiveExtrusion3D::getInnerLoopInUpperPolygon(const rg_INT& index) const
{
    rg_Array<rg_Point3D> output;
    rg_INT numOfPoints=listOfInnerLoop[index].getSize();
    output.setSize(numOfPoints);
    rg_REAL  halfHeight=height/2.0;

    for( rg_INT i=0; i < numOfPoints; i++ )
    {
        output[i]=listOfInnerLoop[index][i]+halfHeight*rg_Point3D(0.0,0.0,1.0);
    }

    return output;
}

rg_Array<rg_Point3D> rg_PrimitiveExtrusion3D::getOuterLoopInLowerPolygon() const
{
    rg_Array<rg_Point3D> output;
    rg_INT numOfPoints=outerLoop.getSize();
    output.setSize(numOfPoints);
    rg_REAL  halfHeight=height/2.0;

    for( rg_INT i=0; i < numOfPoints; i++ )
    {
        output[i]=outerLoop[numOfPoints-1-i]-halfHeight*rg_Point3D(0.0,0.0,1.0);
    }
    return output;
}

rg_Array<rg_Point3D> rg_PrimitiveExtrusion3D::getInnerLoopInLowerPolygon(const rg_INT& index) const
{
    rg_Array<rg_Point3D> output;
    rg_INT numOfPoints=listOfInnerLoop[index].getSize();
    output.setSize(numOfPoints);
    rg_REAL  halfHeight=height/2.0;

    for( rg_INT i=0; i < numOfPoints; i++ )
    {
        output[i]=listOfInnerLoop[index][numOfPoints-1-i]-halfHeight*rg_Point3D(0.0,0.0,1.0);
    }

    return output;
}

rg_Point3D rg_PrimitiveExtrusion3D::getPointInInnerLoop(const rg_INT& innerLoopIndex,
                                              const rg_INT& ptIndex) const
{
    return listOfInnerLoop[innerLoopIndex][ptIndex];
}

rg_Array< rg_Array<rg_Point3D> > rg_PrimitiveExtrusion3D::getListOfInnerLoopInUpperPolygon() const
{
    rg_INT numOfInnerLoop=rg_PrimitiveExtrusion3D::getNumOfInnerLoop();
    rg_Array< rg_Array<rg_Point3D> > output(numOfInnerLoop);
    for( rg_INT i=0; i < numOfInnerLoop; i++ )
    {
        rg_Array<rg_Point3D> innerLoopInUpperPolygon
                    = rg_PrimitiveExtrusion3D::getInnerLoopInUpperPolygon(i);
        output.setAt(i,innerLoopInUpperPolygon);
    }
    return output;;
}

rg_Array< rg_Array<rg_Point3D> > rg_PrimitiveExtrusion3D::getListOfInnerLoopInLowerPolygon() const
{
    rg_INT numOfInnerLoop=rg_PrimitiveExtrusion3D::getNumOfInnerLoop();
    rg_Array< rg_Array<rg_Point3D> > output(numOfInnerLoop);
    for( rg_INT i=0; i < numOfInnerLoop; i++ )
    {
        rg_Array<rg_Point3D> innerLoopInLowerPolygon
                    = rg_PrimitiveExtrusion3D::getInnerLoopInLowerPolygon(i);
        output.setAt(i,innerLoopInLowerPolygon);
    }
    return output;
}

// set function
void   rg_PrimitiveExtrusion3D::setHeight(const rg_REAL& tHeight)
{
    height=tHeight;
}

void   rg_PrimitiveExtrusion3D::setOuterLoop(const rg_Array<rg_Point3D>& tOuterLoop)
{
    outerLoop=tOuterLoop;
}

void   rg_PrimitiveExtrusion3D::setInnerLoop(const rg_INT& index,
                                        const rg_Array<rg_Point3D>& tInnerLoop)
{
    listOfInnerLoop[index]=tInnerLoop;
}

void   rg_PrimitiveExtrusion3D::setListOfInnerLoop(const rg_Array< rg_Array<rg_Point3D> >& tListInnerLoops)
{
    listOfInnerLoop=tListInnerLoops;
}

void   rg_PrimitiveExtrusion3D::setPrimitiveExtrusion(const rg_PrimitiveExtrusion3D& temp)
{
    rg_Geometry::setID(temp.getID());
    rg_Geometry::setType(temp.getType());
    height=temp.height;
    listOfInnerLoop=temp.listOfInnerLoop;
    outerLoop=temp.outerLoop;
}


rg_PrimitiveExtrusion3D& rg_PrimitiveExtrusion3D::operator=(const rg_PrimitiveExtrusion3D& temp)
{
    if ( &temp == this )
    {
        return *this;
    }

    rg_PrimitiveExtrusion3D::setPrimitiveExtrusion(temp);
    return *this;
}


