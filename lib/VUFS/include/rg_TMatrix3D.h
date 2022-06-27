//********************************************************************
//
//	  FILENAME    : rg_TMatrix3D.h
//	  
//    DESCRIPTION : 
//           This consists of the interface of class rg_TMatrix3D.
//           rg_TMatrix3D is for transformaion in 3 dimension    
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Tae-Bum Jang
//    START DATE  : 1997.  7. 10    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//    History     : 1. Insert 
//                            void rotateX(const rg_REAL& cosine,
//                                         const rg_REAL& sine   )
//                            void rotateY(const rg_REAL& cosine,
//                                         const rg_REAL& sine   )
//                            void rotateZ(const rg_REAL& cosine,
//                                         const rg_REAL& sine   )
//					2. Updated by JungHyun Ryu
//							  void rg_TMatrix3D::rotateX( const rg_REAL& theta )
//							  void rg_TMatrix3D::rotateY( const rg_REAL& theta )
//							  void rg_TMatrix3D::rotateZ( const rg_REAL& theta )
//							  void rg_TMatrix3D::rotateX(const rg_REAL& cosine,
//							                        const rg_REAL& sine   )
//							  void rg_TMatrix3D::rotateY(const rg_REAL& cosine,
//							                        const rg_REAL& sine   )
//					          void rg_TMatrix3D::rotateZ(const rg_REAL& cosine,
//							                        const rg_REAL& sine   )
//							  void rg_TMatrix3D::projectX( )
//							  void rg_TMatrix3D::projectY( )
//							  void rg_TMatrix3D::projectZ( )
//							  void rg_TMatrix3D::translate( const rg_Point3D& movement)
//							  void rg_TMatrix3D::reflectToPoint(const rg_Point3D& center)
//							  void rg_TMatrix3D::scale( const rg_REAL& xScale,
//							                       const rg_REAL& yScale,
//												   const rg_REAL& zScale )
//       
//                   3. Insert by Tae-boom Jang
//                           void convertAxis( const rg_Point3D& fromOrigin, const rg_Point3D& fromXAxis, const rg_Point3D fromZAxis,
//                                             const rg_Point3D& fromOrigin, const rg_Point3D& fromXAxis, const rg_Point3D fromZAxis );
//
//
//*********************************************************************



#ifndef _RG_TMATRIX3D_H
#define _RG_TMATRIX3D_H

#include <math.h>
#include "rg_Const.h"
#include "rg_Point3D.h"
#include "rg_Matrix.h"
#include "rg_TMatrix2D.h"
// rg_TMatrix3D is 
//   the abbreviation of "Transformating rg_Matrix"
// Moreover, rg_TMatrix3D is "Homgeneous".

class rg_TMatrix3D: public rg_Matrix
{
public:
//// constructors & destructor    /////////
    rg_TMatrix3D();
    rg_TMatrix3D(const rg_TMatrix3D    &tMat);
    rg_TMatrix3D(const rg_Matrix       &tMat);
    rg_TMatrix3D(const rg_TMatrix2D    &tMat);
    ~rg_TMatrix3D();

//// set function                 ////////
    void reset();

//// Transfromation               ////////
    void rotateX(const rg_REAL& theta);
    void rotateY(const rg_REAL& theta);
    void rotateZ(const rg_REAL& theta);
    void rotateX(const rg_REAL& cosine,
                 const rg_REAL& sine   );
    void rotateY(const rg_REAL& cosine,
                 const rg_REAL& sine   );
    void rotateZ(const rg_REAL& cosine,
                 const rg_REAL& sine   );
    void rotate( const rg_Point3D& fromVector,
                 const rg_Point3D&   toVector );

    void rotateArbitraryAxis(const rg_Point3D& axis,
                             const rg_REAL&  theta);

    
    void projectX();
    void projectY();
    void projectZ();
    void project( const rg_REAL& theta,
                  const rg_REAL&   phi );
    void translate(const rg_Point3D& movement);

    void reflectToPoint(const rg_Point3D& center=rg_Point3D(0.,0.,0.));

	void scale( const rg_REAL& globalScale );
    void scale( const rg_REAL& xScale,
                const rg_REAL& yScale,
                const rg_REAL& zScale );    

    void convertAxis( const rg_Point3D& oldOrigin, const rg_Point3D& oldXAxis, const rg_Point3D& oldZAxis,
                      const rg_Point3D& newOrigin, const rg_Point3D& newXAxis, const rg_Point3D& newZAxis );
    //// operator overloading         ////////
    rg_Point3D   operator*(const rg_Point3D&  tPoint) const;
    rg_Matrix  operator*(const rg_Matrix&  tMat) const;

    static rg_REAL arcTan(const rg_REAL& y,const rg_REAL& x);
};

// arcTangent



#endif


