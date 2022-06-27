#include "rg_TMatrix3D.h"
#include "rg_Const.h"
#include "rg_RelativeOp.h"

#include <iostream>//test

//// constructors & destructor    /////////
rg_TMatrix3D::rg_TMatrix3D()
 :rg_Matrix( 4 , 4)
{
    rg_Matrix::makeIdentity();
}

rg_TMatrix3D::rg_TMatrix3D(const rg_TMatrix2D    &tMat)
:rg_Matrix( 4 , 4)
{
    rg_Matrix::makeIdentity();
    
    //row 0-1
    for(rg_INT i=0; i < 2; i++)
    {
        element[i][0]=tMat.getElement(i,0);
        element[i][1]=tMat.getElement(i,1);
        //col 2 has no change
        element[i][3]=tMat.getElement(i,2);
    }
    
    //row 2 has no change

    //row 3
    element[3][0]=tMat.getElement(2,0);
    element[3][1]=tMat.getElement(2,1);
    //col 2 has no change
    element[3][3]=tMat.getElement(2,2);
}


rg_TMatrix3D::rg_TMatrix3D(const rg_TMatrix3D      &tMat)
 :rg_Matrix(tMat)
{
}

rg_TMatrix3D::rg_TMatrix3D(const rg_Matrix& tMat)
 :rg_Matrix(tMat)
{
}

rg_TMatrix3D::~rg_TMatrix3D()
{
}

//// set function                 ////////

void rg_TMatrix3D::reset()
{
    if( element != rg_NULL )
    {
       makeIdentity();
    }
}

//// Transfromation               ////////
void rg_TMatrix3D::rotateX( const rg_REAL& theta )
{
    const rg_REAL thetaRad=theta;

	rg_REAL mapping[ 2 ][ 2 ];
    mapping[0][0] =  cos(thetaRad);
    mapping[0][1] = -sin(thetaRad);
    mapping[1][0] =  sin(thetaRad);
    mapping[1][1] =  cos(thetaRad);
	
	rg_REAL willBeChanged[ 2 ][ 4 ];
	rg_INT i = 0;
	for(i = 0;i < 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			willBeChanged[ i ][ j ] = (*this)[i + 1][ j ];
		}
	}

	for(i = 1;i <= 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			(*this)[ i ][ j ] =  mapping[i - 1][ 0 ] * willBeChanged[ 0 ][ j ]
                               + mapping[i - 1][ 1 ] * willBeChanged[ 1 ][ j ];
		}
	}
}

void rg_TMatrix3D::rotateY( const rg_REAL& theta )
{
    const rg_REAL thetaRad=theta;


    rg_REAL mapping[ 2 ][ 2 ];
	mapping[0][0] =  cos(thetaRad);
	mapping[0][1] =  sin(thetaRad);
	mapping[1][0] = -sin(thetaRad);
  	mapping[1][1] =  cos(thetaRad);
	
	
	rg_REAL willBeChanged[ 2 ][ 4 ];
	rg_INT i = 0;
	for(i = 0;i < 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			willBeChanged[ i ][ j ] = (*this)[i * 2][ j ];
		}
	}

	for(i = 0;i < 3;i+=2)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			(*this)[ i ][ j ] =  mapping[i / 2][ 0 ] * willBeChanged[ 0 ][ j ]
                               + mapping[i / 2][ 1 ] * willBeChanged[ 1 ][ j ];

		}
	}

}

void rg_TMatrix3D::rotateZ( const rg_REAL& theta )
{
    const rg_REAL thetaRad=theta;

    rg_REAL mapping[ 2 ][ 2 ];
    mapping[0][0] =  cos(thetaRad);
    mapping[0][1] = -sin(thetaRad);
    mapping[1][0] =  sin(thetaRad);
    mapping[1][1] =  cos(thetaRad);
	
	rg_REAL willBeChanged[ 2 ][ 4 ];
	rg_INT i = 0;
	for(i = 0;i < 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			willBeChanged[ i ][ j ] = (*this)[ i ][ j ];
		}
	}

	for(i = 0;i < 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			(*this)[ i ][ j ] =   mapping[ i ][ 0 ] * willBeChanged[ 0 ][ j ]
                                + mapping[ i ][ 1 ] * willBeChanged[ 1 ][ j ];
		}
	}
}

void rg_TMatrix3D::rotateX(const rg_REAL& cosine,
                      const rg_REAL& sine   )
{
    if (  rg_NE( cosine*cosine+sine*sine , 1.0 ) )
    {
        return;
    }

	//    
	rg_REAL mapping[ 2 ][ 2 ];
    mapping[0][0] =  cosine;
    mapping[0][1] = -sine;
    mapping[1][0] =  sine;
    mapping[1][1] =  cosine;
	
	rg_REAL willBeChanged[ 2 ][ 4 ];
	rg_INT i = 0;
	for(i = 0;i < 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			willBeChanged[ i ][ j ] = (*this)[i + 1][ j ];
		}
	}

	for(i = 1;i <= 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			(*this)[ i ][ j ] =  mapping[i - 1][ 0 ] * willBeChanged[ 0 ][ j ] 
                               + mapping[i - 1][ 1 ] * willBeChanged[ 1 ][ j ];
		}
	}
}

void rg_TMatrix3D::rotateY(const rg_REAL& cosine,
                      const rg_REAL& sine   )
{
    if ( rg_NE( cosine*cosine+sine*sine , 1.0 ) )
    {
        return;
    }

	//
    rg_REAL mapping[ 2 ][ 2 ];
	mapping[0][0] =  cosine;
	mapping[0][1] =  sine;
	mapping[1][0] = -sine;
  	mapping[1][1] =  cosine;
	
	
	rg_REAL willBeChanged[ 2 ][ 4 ];
	rg_INT i = 0;
	for(i = 0;i < 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			willBeChanged[ i ][ j ] = (*this)[i * 2][ j ];
		}
	}

	for(i = 0;i < 3;i+=2)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			(*this)[ i ][ j ] =  mapping[i / 2][ 0 ] * willBeChanged[ 0 ][ j ]
                               + mapping[i / 2][ 1 ] * willBeChanged[ 1 ][ j ];
		}
	}

}

void rg_TMatrix3D::rotateZ(const rg_REAL& cosine,
                      const rg_REAL& sine   )
{
    if ( rg_NE( cosine*cosine+sine*sine , 1.0 ) )
    {
        return;
    }

    rg_REAL mapping[ 2 ][ 2 ];
    mapping[0][0] =  cosine;
    mapping[0][1] = -sine;
    mapping[1][0] =  sine;
    mapping[1][1] =  cosine;
	
	rg_REAL willBeChanged[ 2 ][ 4 ];
	rg_INT i = 0;
	for(i = 0;i < 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			willBeChanged[ i ][ j ] = (*this)[ i ][ j ];
		}
	}

	for(i = 0;i < 2;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			(*this)[ i ][ j ] =  mapping[ i ][ 0 ] * willBeChanged[ 0 ][ j ]
                               + mapping[ i ][ 1 ] * willBeChanged[ 1 ][ j ];
		}
	}
}


void rg_TMatrix3D::rotate( const rg_Point3D& fromVector,
                      const rg_Point3D&   toVector )
{

    if (    fromVector == rg_Point3D(0.,0.,0.)
         ||   toVector == rg_Point3D(0.,0.,0.)  
         || fromVector == toVector        )
    {
        return;
    }

    rg_Point3D unitFrom=fromVector/ (fromVector.magnitude());
    rg_Point3D unitTo=toVector/(toVector.magnitude());

    const rg_REAL theta = acos(unitFrom%unitTo);
    rg_Point3D axis=unitFrom*unitTo;
    const rg_REAL PHI=atan(1.0)*4.0;
    if ( rg_EQ( theta, 0.0 ) )
    {
        return;
    } 
    else if( rg_EQ( theta, PHI) )
    {
        return;
    }
    else
    {
        rotateArbitraryAxis(axis,theta);
    }
}

void rg_TMatrix3D::rotateArbitraryAxis(const rg_Point3D& axis,
                                  const rg_REAL&   theta)
{
    rg_REAL axisX=axis.getX();
    rg_REAL axisY=axis.getY();
    rg_REAL axisZ=axis.getZ();

    rg_REAL cosOfAxisX=0.0;
    rg_REAL sinOfAxisX=0.0;
    rg_REAL cosOfAxisY=0.0;
    rg_REAL sinOfAxisY=0.0;
    rg_REAL degreeZ=theta;

    if ( rg_LT( axisZ, 0.0 ) )
    {
        axisX=-axisX;
        axisY=-axisY;
        axisZ=-axisZ;
        degreeZ=-degreeZ;
    }

    if (   rg_NE( axisX, 0.0 )
        || rg_NE( axisY, 0.0 ) )
    {
        const rg_REAL distanceOnlyYZ=sqrt( axisY*axisY + axisZ*axisZ );
        const rg_REAL distance=axis.magnitude();
        
        if ( rg_EQ( distanceOnlyYZ , 0.0 ) )
        {
            cosOfAxisX =1.0;
            sinOfAxisY =0.0;
        }
        else
        {
            cosOfAxisX=axisZ/distanceOnlyYZ;
            sinOfAxisX=axisY/distanceOnlyYZ;
        }
        cosOfAxisY=distanceOnlyYZ/distance;
        sinOfAxisY=-axisX/distance;
    }
    else
    {
        rotateZ( degreeZ );
        return;
    }
    rotateX( cosOfAxisX,  sinOfAxisX );
    rotateY( cosOfAxisY,  sinOfAxisY );
    rotateZ( degreeZ );
    rotateY( cosOfAxisY, -sinOfAxisY );
    rotateX( cosOfAxisX, -sinOfAxisX);
}

void rg_TMatrix3D::projectX( )
{
    rg_TMatrix3D mapping;

	for(rg_INDEX i = 0;i < 4;i++)
		(*this)[ 0 ][ i ] = 0.0;
}

void rg_TMatrix3D::projectY()
{
    rg_TMatrix3D mapping;

	for(rg_INDEX i = 0;i < 4;i++)
		(*this)[ 1 ][ i ] = 0.0;

}

void rg_TMatrix3D::projectZ()
{
    rg_TMatrix3D mapping;

	for(rg_INDEX i = 0;i < 4;i++)
		(*this)[ 2 ][ i ] = 0.0;
}

void rg_TMatrix3D::project( const rg_REAL& theta,
                       const rg_REAL&   phi )
{
    rotateX(theta);
    rotateY(  phi);
    projectZ();
}

void rg_TMatrix3D::translate( const rg_Point3D& movement )
{
	rg_REAL mapping[ 3 ];

    mapping[0] = movement.getX();
    mapping[1] = movement.getY();
    mapping[2] = movement.getZ();

	for(rg_INDEX i = 0;i < 3;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			(*this)[ i ][ j ] =  (*this)[ i ][ j ] 
                               + mapping[ i ] * (*this)[ 3 ][ j ];
		}
	}
}

void rg_TMatrix3D::reflectToPoint(const rg_Point3D& center)
{
    translate(-center);// move center to origin

	for(rg_INDEX i = 0;i < 3;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			(*this)[ i ][ j ] = - (*this)[ i ][ j ];
		}
	}
	//

    translate(center);//move origin to cener
}

void rg_TMatrix3D::scale( const rg_REAL& globalScale )
{
	(*this)[ 3 ][ 3 ] = 1/globalScale;
}

void rg_TMatrix3D::scale( const rg_REAL& xScale,
                     const rg_REAL& yScale,
                     const rg_REAL& zScale )
{

	rg_REAL scale[ 3 ];
	scale[ 0 ] = xScale;
	scale[ 1 ] = yScale;
	scale[ 2 ] = zScale;

	for(rg_INDEX i = 0;i < 3;i++)
	{
		for(rg_INDEX j = 0;j < 4;j++)
		{
			(*this)[ i ][ j ] = scale[ i ] * (*this)[ i ][ j ];
		}
	}

}
void rg_TMatrix3D::convertAxis( const rg_Point3D& oldOrigin, const rg_Point3D& oldXAxis, const rg_Point3D& oldZAxis,
                           const rg_Point3D& newOrigin, const rg_Point3D& newXAxis, const rg_Point3D& newZAxis )
{

    rg_TMatrix3D rotation;
    rotation.rotate(oldZAxis,newZAxis);

    rg_Point3D transientXAxis=rotation*oldXAxis;
    rotation.rotate(transientXAxis,newXAxis);

/*
    oldYAxis=oldZAxis*oldXAxis;
    newYAxis=newZAxis*newXAxis;
    
    rg_Matrix A(3,3);
    A(0,0)=oldXAxis.getX();
    A(0,1)=oldXAxis.getY();
    A(0,2)=oldXAxis.getZ();
    A(1,0)=oldYAxis.getX();
    A(1,1)=oldYAxis.getY();    
    A(1,1)=oldYAxis.getZ();    
    A(2,0)=oldZAxis.getX();
    A(2,1)=oldZAxis.getY();
    A(2,2)=oldZAxis.getZ();
    rg_Matrix B(3,3);
    
    B(0,0)=newXAxis.getX();
    B(1,0)=newXAxis.getY();
    B(2,0)=newXAxis.getZ();
    B(0,1)=newYAxis.getX();
    B(1,1)=newYAxis.getY();    
    B(2,1)=newYAxis.getZ();    
    B(0,2)=newZAxis.getX();
    B(1,2)=newZAxis.getY();
    B(2,2)=newZAxis.getZ();    
    
    rg_Matrix C= B*A;
    
    rg_TMatrix3D rotation;
    for(rg_INT i=0;i < 3; i++ )
    {
        for(rg_INT j=0; j < 3; j++ )
        {
            rotation[i][j]=C[i][j];
        }
    }
*/
    rg_TMatrix3D translation;
    translation.translate(newOrigin-oldOrigin);
    (*this)=translation*rotation*(*this);
}


//// operator overloading         ////////
rg_Point3D rg_TMatrix3D::operator*(const rg_Point3D& tPoint) const
{
    rg_Matrix hPt(4,1);// abbrieviation of " Homogenous point"
    hPt[0][0]=tPoint.getX();
    hPt[1][0]=tPoint.getY();
    hPt[2][0]=tPoint.getZ();
    hPt[3][0]=1.;

    rg_Matrix result= (*(this)) * hPt;
    if ( result[3][0] != 0. )
    {
        for( rg_INT i=0; i < 4; i++)
        {
            result[i][0]=result[i][0]/result[3][0];
        }
    }

    rg_Point3D output;
    output.setX(result[0][0]);
    output.setY(result[1][0]);
    output.setZ(result[2][0]);

    return output;
}

rg_Matrix rg_TMatrix3D::operator*(const rg_Matrix&  tMat) const
{
    rg_Matrix lftOperand( (*this) );
    
    return lftOperand*tMat;
}

rg_REAL rg_TMatrix3D::arcTan(const rg_REAL& y,const rg_REAL& x)
{
    if (    rg_EQ(y, 0.0)
         && rg_EQ(x, 0.0) )
    {
         return 0.0;
    }else
    {
        return atan2(y,x);
    }
}

    

    
