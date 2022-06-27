//******************************************************************************
//
//	  FILENAME    : rg_PrimitiveExtrusion3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_PrimitiveExtrusion3D 
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

#ifndef _RG_PRIMITIVEEXTRUSION3D_H
#define _RG_PRIMITIVEEXTRUSION3D_H

#include "rg_Array.h"
#include "rg_Point3D.h"
#include "rg_Geometry.h"
#include "rg_Const.h"

class rg_PrimitiveExtrusion3D: virtual public rg_Geometry
{
protected:
     rg_Array< rg_Point3D > outerLoop;
     rg_Array< rg_Array<rg_Point3D> > listOfInnerLoop;
     rg_REAL  height;
public:
    rg_PrimitiveExtrusion3D();
    rg_PrimitiveExtrusion3D(const GeometryType& tType,const unsigned rg_INT& tID=0);
    rg_PrimitiveExtrusion3D(const rg_PrimitiveExtrusion3D& temp);
    ~rg_PrimitiveExtrusion3D();

    // get function
    rg_INT   getNumOfPointInOuterLoop() const;
    rg_INT   getNumOfPointInInnerLoop(const rg_INT& index) const;
    rg_INT   getNumOfInnerLoop() const;
    rg_Array<rg_Point3D> getOuterLoop() const;
    rg_Array<rg_Point3D> getInnerLoop( const rg_INT& index) const;
    rg_Array<rg_Point3D> getOuterLoopInUpperPolygon() const;
    rg_Array<rg_Point3D> getInnerLoopInUpperPolygon(const rg_INT& index) const;
    rg_Array<rg_Point3D> getOuterLoopInLowerPolygon() const;
    rg_Array<rg_Point3D> getInnerLoopInLowerPolygon(const rg_INT& index) const;

    rg_Array< rg_Array<rg_Point3D> > getListOfInnerLoopInUpperPolygon() const;
    rg_Array< rg_Array<rg_Point3D> > getListOfInnerLoopInLowerPolygon() const;


    rg_Point3D getPointInOuterLoop(const rg_INT& index) const;
    rg_Point3D getPointInInnerLoop(const rg_INT& innerLoopIndex,
                              const rg_INT& ptIndex) const;
    rg_REAL  getHeight() const;

    // set function
    void   setOuterLoop(const rg_Array<rg_Point3D>& tOuterLoop);
    void   setInnerLoop(const rg_INT& index,
                        const rg_Array<rg_Point3D>& tInnerLoop);
    void   setListOfInnerLoop(const rg_Array< rg_Array<rg_Point3D> >& tInnererLoops);
    void   setPrimitiveExtrusion(const rg_PrimitiveExtrusion3D& temp);
    void   setHeight(const rg_REAL& tHeight);

    rg_PrimitiveExtrusion3D& operator=(const rg_PrimitiveExtrusion3D& temp);
};

#endif


