//********************************************************************
//
//	  FILENAME    : rg_RBzSurface3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_BzSurface3D
//           which define a rational Bezier rg_Surface and its property. 
//                          
//	  CLASS NAME  : rg_RBzSurface3D
//
//    BASE CLASS  : rg_BzSurface3D
//      
//    AUTHOR      : Deok-Soo Kim, Taeboom Jang
//
//    HISTORY     : 	
//                   1. insert following funtions by Taeboom, Jang at 2000.2.21
//                      rg_BzSurface3D getUHodograph() const;
//                      rg_BzSurface3D getVHodograph() const;
//
//
//    START DATE  : 1998. 8.13
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_RBZSURFACE3D_H
#define _RG_RBZSURFACE3D_H

#include "rg_BzSurface3D.h"
#include "rg_Const.h"
#include "rg_Point3D.h"
#include "rg_Matrix.h"

class rg_RBzSurface3D : public rg_BzSurface3D
{
protected:
	rg_REAL** weight_vector;

public:
    ////  Constructor & Destructor
    rg_RBzSurface3D();
    rg_RBzSurface3D( const rg_DEGREE& uDegree, const rg_DEGREE& vDegree );
    rg_RBzSurface3D( const rg_DEGREE& uDegree,
	            const rg_DEGREE& vDegree,
			    const rg_Point3D** ctrlPts );
    rg_RBzSurface3D( const rg_DEGREE& uDegree,
	            const rg_DEGREE& vDegree,
			    const rg_Point3D** ctrlPts ,
                const rg_REAL**  weightVector);

    rg_RBzSurface3D( const rg_RBzSurface3D& surface );
    virtual ~rg_RBzSurface3D();

    ////  Access elements
    rg_REAL   getWeight( const rg_INDEX& i, const rg_INDEX& j ) const;

    //return pointer to the copy of control points.
    rg_REAL**  getWeightVector() const;


    void    setOneWeight( const rg_INDEX& i,
		                  const rg_INDEX& j,
					      const rg_REAL& weight );
    void    setWeightVector(const rg_REAL** weightVector);
    void    setDegree( const rg_DEGREE& newUDegree,
                       const rg_DEGREE& newVDegree );
    void    setSurface( const rg_RBzSurface3D& temp);
    // 
    rg_RBzSurface3D getUHodograph() const;
    rg_RBzSurface3D getVHodograph() const;
    //  Operations
    virtual rg_Point3D evaluatePt( const rg_PARAMETER& u,
		                      const rg_PARAMETER& v );

    // operator overloading
    rg_RBzSurface3D& operator=(const rg_RBzSurface3D& temp);

};

#endif


