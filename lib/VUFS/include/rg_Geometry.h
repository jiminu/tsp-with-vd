//******************************************************************************
//
//    FILENAME    : rg_Geometry.h
//    
//    DESCRIPTION : 
//           This is the interface of the class rg_Geometry 
//           which will be the base class of all kinds of Geometries. 
//                          
//    CLASS NAME  : rg_Geometry
//
//    BASE CLASS  : NONE
//      
//
//    AUTHOR      : Deok-Soo Kim, Tae-Bum Jang
//    START DATE  : 1998. 4.2
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering
//                          Hanyang University, Seoul Korea         
//
//******************************************************************************

#ifndef _GEOMETRY_INCLUDED
#define _GEOMETRY_INCLUDED

#include "rg_Const.h"
//#include "Defconst.h"

enum GeometryType
{
    GEOMETRY,
    PRIMITIVE_CYLINDER,
    PRIMITIVE_FRUSTUM,
    PRIMITIVE_BOX,
    CYLINDER,
    FRUSTUM,
    BOX,
    PARTIAL_CYLINDER,
    LINK,
    PRIMITIVE_EXTRUSION,
    EXTRUSION,
    REVOLVED_SURFACE,
    OUT_LINE
};

class rg_Geometry
{
protected:
    unsigned rg_INT id;
    GeometryType type;
public:
////    Constructor & Destructor.
    rg_Geometry( const unsigned rg_INT& tID=0);
    rg_Geometry( const GeometryType& tType,const unsigned rg_INT &tID=0);
    rg_Geometry( const rg_Geometry &temp);
    virtual ~rg_Geometry();

////    Get Functions.
    unsigned rg_INT getID() const;
    GeometryType getType() const;


////    Set Functions.
    void setID(const unsigned rg_INT &newID);
    void setType(const GeometryType& tType);

    void setGeometry(const rg_Geometry& temp);
    rg_Geometry& operator=(const rg_Geometry& temp);
////    Abstract Virtual functions
};

#endif


