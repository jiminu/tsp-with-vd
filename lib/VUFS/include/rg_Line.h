//********************************************************************
//
//	  FILENAME    : rg_Line.h
//	  
//    DESCRIPTION : 
//           This consists of the interface and implementation
//           of class rg_Line.      
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Deok-Soo Kim, Soon-Woong Lee
//    START DATE  : 8 Mar. 1997    
//
//    HISTORY:
//   1. Add : By Young-Song Cho,  11 Jul. 1997
//              rg_FLAG operator==(const rg_Line<PT> &line) const;
//   2. Update: By Tae-Bum Jang, 1997.7.21
//            	rg_REAL getLength() ;
//           -> rg_REAL   getLength() const;
//   3. Add : By Tae-Bum Jang, 1997. 7.21
//              rg_REAL   getDistance(const PT& point);
//   4. Add : By Tae-Bum Jang, 1997. 7.21
//              #include <math.h>
//                [Why? "rg_REAL   getDistance(const PT& point)"  use "sqrt()" in math.h
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_LINE_H
#define _RG_LINE_H

#include "rg_Const.h"
#include "rg_RelativeOp.h"
#include <math.h>

template <class PT>
class rg_Line
{
protected:
	PT sp;
	PT ep;
public: 

	// Construction/Destruction
	rg_Line();
	rg_Line(const PT &sPt, const PT &ePt);
	rg_Line(const rg_Line<PT> &line);
	~rg_Line();            
              
	// Access elements	                       
	PT      getSP() const;
	PT      getEP() const;  

//	In 2D the following function is meaningful but in 3D is not.
//	rg_REAL  getTangent() const;
	
	void   setSP(const PT &pt);
	void   setEP(const PT &pt);
	void   setLine(const PT &sPt, const PT &ePt);
	void   setLine(const rg_Line<PT> &line);

	// Operations
//	rg_REAL getLength() ;
    rg_REAL   getLength() const;
    rg_REAL   getDistance(const PT& point) const;

    // Operator Overloading
    rg_FLAG operator==(const rg_Line<PT> &line) const;
//	In 2D the following function is meaningful but in 3D is not.
//	rg_INT    isLeftSide(const rg_Line<PT> &line) const;
};

////////////////////////////////////////////////////////////////////////
//  Implementation of template class


////////////////////////////////////////////////////////////////////////
// Construction and Destruction

template <class PT>
rg_Line<PT>::rg_Line()
{
}

template <class PT>
rg_Line<PT>::rg_Line(const PT &sPt, const PT &ePt)
{
	sp = sPt;
	ep = ePt;
}

template <class PT>
rg_Line<PT>::rg_Line(const rg_Line<PT> &line)
{
	sp = line.sp;
	ep = line.ep;
}

template <class PT>
rg_Line<PT>::~rg_Line()
{
}

////////////////////////////////////////////////////////////////////////
// Access elements
template <class PT>
PT rg_Line<PT>::getSP() const
{
	return sp;
}

template <class PT>
PT rg_Line<PT>::getEP() const
{
	return ep;
}

/*
template <class PT>
rg_REAL rg_Line::getTangent() const
{
	return (endP.getY()-startP.getY())/(endP.getX()-startP.getX());
}
*/

template <class PT>
void rg_Line<PT>::setSP(const PT &pt)
{
	sp = pt;
}

template <class PT>
void rg_Line<PT>::setEP(const PT &pt)
{
	ep = pt;
}

template <class PT>
void rg_Line<PT>::setLine(const PT &sPt, const PT &ePt)
{
	sp.setPoint(sPt);
	ep.setPoint(ePt);
}

template <class PT>
void rg_Line<PT>::setLine(const rg_Line<PT> &line)
{
	sp = line.sp;
	ep = line.ep;
}
                  

////////////////////////////////////////////////////////////////////////
// Operations
template <class PT>
rg_REAL rg_Line<PT>::getLength() const
{
	PT     length(ep - sp);
	return length.magnitude();
}

// Following is the illustration of 
//                    " rg_REAL rg_Line<PT>::getDistance(const PT& point) const "
//
//                  * point
//                 /|
//                / |
//               /  |     
//         roof /   | Distance between point and line (sp,ep)
//             /    |
//            /     |  
//           *______+ ___________* 
//           sp                 ep
//    Distance = (length of roof) * ( sine of angle between roof and this line )
template <class PT>
rg_REAL rg_Line<PT>::getDistance(const PT& point) const
{
    PT   roofVector( point-sp );
    PT   lineVector( ep   -sp );
    PT   unitRoofVector        = roofVector.getUnitVector();
    PT   unitLineVector        = lineVector.getUnitVector();
    rg_REAL roofLength            = roofVector.magnitude();

    rg_REAL cosine = unitRoofVector % unitLineVector;
    rg_REAL sine   = sqrt(1-cosine*cosine);

    rg_REAL distance = roofLength*sine;

    return distance;
}


////////////////////////////////////////////////////////////////////////
// Operator Overloading
template <class PT>
rg_FLAG rg_Line<PT>::operator==(const rg_Line<PT> &line) const
{
    if ( (sp == line.sp && ep == line.ep)
         || (sp == line.ep && ep == line.sp) )
    {
        return rg_TRUE;
    }
    else 
        return rg_FALSE;
}

/* 
template <class PT>
rg_INT rg_Line::isLeftSide(const rg_Line<PT> &line) const
{
	PT vector1(endP.getX()-startP.getX(), endP.getY()-startP.getY());
	PT vector2(line.getEndP().getX()-line.getStartP().getX(), 
		            line.getEndP().getY()-line.getStartP().getY());

	rg_REAL det = vector1.getX()*vector2.getY() - vector1.getY()*vector2.getX();
	if(rg_POS(det, rg_SYSTEM_RES)) return rg_TRUE;
	else return rg_FALSE;
}
*/

#endif


