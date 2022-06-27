//******************************************************************************
//
//    FILENAME    : rg_Geometry.cpp
//    
//    DESCRIPTION : 
//           This is the implementation of the class rg_Geometry 
//           which will be the base class of all kinds of Geometrys. 
//                          
//	  CLASS NAME  : rg_Geometry
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

#include "rg_Geometry.h"



rg_Geometry::rg_Geometry( const unsigned rg_INT& tID)
:id(tID),type(GEOMETRY)
{
}

rg_Geometry::rg_Geometry( const GeometryType& tType,const unsigned rg_INT &tID)
:id(tID), type(tType)
{
}

rg_Geometry::rg_Geometry( const rg_Geometry &temp )
:id(temp.id), type(temp.type)
{
}

rg_Geometry::~rg_Geometry()
{
}


////    Get Functions.
unsigned rg_INT rg_Geometry::getID() const
{
    return id;
}

GeometryType rg_Geometry::getType() const
{
    return type;
}


////    Set Functions.
void rg_Geometry::setID(const unsigned rg_INT &newID)
{
    id=newID;
}

void rg_Geometry::setType(const GeometryType& tType)
{
    type=tType;
}

void rg_Geometry::setGeometry(const rg_Geometry& temp)
{
    id=temp.id;
    type=temp.type;
}

rg_Geometry& rg_Geometry::operator=(const rg_Geometry& temp)
{
    if ( this == &temp )
    {
        return*this;
    }

    id=temp.id;
    type=temp.type;

    return *this;
}

