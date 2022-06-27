//********************************************************************
//
//	  FILENAME    : rg_BzSurface3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_BzSurface3D
//           which define a Bezier rg_Surface and its property. 
//                          
//	  CLASS NAME  : rg_BzSurface3D
//
//    BASE CLASS  : rg_Surface
//      
//    AUTHOR      : Deok-Soo Kim, Dong-Gyou Lee
//
//    HISTORY     : 
//                   1. insert following functions by Taeboom, Jang 
//                      void       setSurface(const rg_BzSurface3D& temp);
//                      rg_BzSurface3D& operator=(const rg_BzSurface3D& temp);
//                   
//                   2. insert following functions by Youngsong Cho
//                      rg_Point3D** evaluatePt(const rg_INT& numIsoparaCrvOfU,
//		                                     const rg_INT& numIsoparaCrvOfV,
//						                     const rg_INT& resolution);  
//
//                   3. insert following functions by Taeboom Jang
//                      rg_BzSurface3D getUHodograph() const;
//                      rg_BzSurface3D getVHodograph() const;
//		                          
//						          
//
//    START DATE  : 23 Mar. 1998    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#ifndef _RG_BZSURFACE3D_INCLUDED
#define _RG_BZSURFACE3D_INCLUDED

#include "rg_Surface.h"
#include "rg_Const.h"
#include "rg_Point3D.h"
#include "rg_Matrix.h"

class rg_BzSurface3D : public rg_Surface
{
protected:
    rg_DEGREE degreeOfU;   
	rg_DEGREE degreeOfV;   
    
	rg_Point3D** ctrlNet;

public:
    ////  Constructor & Destructor
    rg_BzSurface3D();
    rg_BzSurface3D( const rg_DEGREE& uDegree, const rg_DEGREE& vDegree );
    rg_BzSurface3D( const rg_DEGREE& uDegree,
		       const rg_DEGREE& vDegree,
			   const rg_Point3D** ctrlPts );
    rg_BzSurface3D( const rg_BzSurface3D& surface );
    virtual ~rg_BzSurface3D();

    ////  Access elements
	rg_FLAG isValidSurface() const;
    rg_DEGREE  getDegreeOfU() const;
	rg_DEGREE  getDegreeOfV() const;
    rg_Point3D getCtrlPt( const rg_INDEX& i, const rg_INDEX& j ) const;

    //return pointer to the copy of control points.
    rg_Point3D** getControlNet() const;

    void    setDegree( const rg_DEGREE& uDegree, const rg_DEGREE& vDegree );
    void    setCtrlPt( const rg_INDEX& i,
		                        const rg_INDEX& j,
								const rg_Point3D& pt );
    void    setControlNet( const rg_INT& numOfRows,
		                   const rg_INT& numOfCols, 
		                   rg_Point3D**  ctrlPts);
    void    setSurface(const rg_BzSurface3D& temp);

    //  Operations
    virtual rg_Point3D evaluatePt( const rg_PARAMETER& u,
		                      const rg_PARAMETER& v );
	rg_Point3D** evaluatePt(const rg_INT& numIsoparaCrvOfU,
		                 const rg_INT& numIsoparaCrvOfV,
						 const rg_INT& resolution);  
    // Derivative 
    rg_BzSurface3D getUHodograph() const;
    rg_BzSurface3D getVHodograph() const;

//  rg_Point3D**     deCasteljau( const rg_PARAMETER& u );
    rg_REAL        bernstein( const rg_INDEX&     i,
		                   const rg_DEGREE&    dgr,
						   const rg_PARAMETER& t);
    rg_REAL       factorial( const rg_REAL& n );

	//	Conversion between power basis & bezier form.
   	void powerToBezierSurface( const rg_DEGREE& uDegree,
		                       const rg_DEGREE& vDegree, 
			                   const rg_REAL paramValuesOfU[],
							   const rg_REAL paramValuesOfV[],
							   const rg_Matrix& powerCoeff );

    // operator overloading
    rg_BzSurface3D& operator=(const rg_BzSurface3D& temp);

};

#endif


