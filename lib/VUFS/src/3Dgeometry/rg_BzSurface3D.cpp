//********************************************************************
//
//	  FILENAME    : rg_BzSurface3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_BzSurface3D
//           which defines a Bezier rg_Surface and represents 
//           its properties. 
//                          
//	  CLASS NAME  : rg_BzSurface3D
//
//    BASE CLASS  : rg_Surface
//
//    AUTHOR      : Deok-Soo Kim, Dong-Gyou Lee
//
//    HISTORY     : 	
//
//    START DATE  : 23 Mar. 1998    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#include "rg_BzSurface3D.h"
#include "rg_CurveSurfaceFunc.h"
#include "rg_RelativeOp.h"

//  Constructor & Destructor
rg_BzSurface3D::rg_BzSurface3D()
{
    degreeOfU=-1;
    degreeOfV=-1;
    ctrlNet = rg_NULL;
}

rg_BzSurface3D::rg_BzSurface3D( const rg_DEGREE& uDegree, const rg_DEGREE& vDegree )
    : degreeOfU( uDegree ), degreeOfV( vDegree )
{
    ctrlNet = new rg_Point3D*[degreeOfU + 1];

	for( rg_INDEX i=0; i <= degreeOfU; i++ )
		ctrlNet[i] = new rg_Point3D[degreeOfV + 1];
}

rg_BzSurface3D::rg_BzSurface3D( const rg_DEGREE& uDegree,
					  const rg_DEGREE& vDegree,
					  const rg_Point3D** ctrlPts )
    : degreeOfU( uDegree ), degreeOfV( vDegree )
{
	ctrlNet = new rg_Point3D*[degreeOfU + 1];

	rg_INDEX i = 0;
	for( i=0; i <= degreeOfU; i++ )
		ctrlNet[i] = new rg_Point3D[degreeOfV + 1];

    for ( i=0; i <= degreeOfU; i++ )
	{
		for( rg_INDEX j=0; j <= degreeOfV; j++ )
			ctrlNet[i][j] = ctrlPts[i][j];
	}
}

rg_BzSurface3D::rg_BzSurface3D( const rg_BzSurface3D& surface )
    : degreeOfU( surface.degreeOfU ), degreeOfV( surface.degreeOfV )
{
	ctrlNet = new rg_Point3D*[degreeOfU + 1];
	rg_INDEX i = 0; 
	for( i=0; i <= degreeOfU; i++ )
		ctrlNet[i] = new rg_Point3D[degreeOfV + 1];

    for ( i=0; i <= degreeOfU; i++ )
	{
		for( rg_INDEX j=0; j <= degreeOfV; j++ )
			ctrlNet[i][j] = surface.ctrlNet[i][j];
	}
}

rg_BzSurface3D::~rg_BzSurface3D()
{
    if( ctrlNet )
	{	
		for( rg_INDEX i=0; i<=degreeOfU; i++ )
			delete[] ctrlNet[i];

        delete [] ctrlNet;
	}
}

//  Access elements
rg_FLAG rg_BzSurface3D::isValidSurface() const
{
	if ( ctrlNet )  {
		return rg_TRUE;
	}  
	else  {
		return rg_FALSE;
	}
}

rg_DEGREE rg_BzSurface3D::getDegreeOfU() const
{
    return degreeOfU;
}

rg_DEGREE rg_BzSurface3D::getDegreeOfV() const
{
    return degreeOfV;
}


rg_Point3D rg_BzSurface3D::getCtrlPt( const rg_INDEX& i, const rg_INDEX& j ) const
{
    return ctrlNet[i][j];
}

rg_Point3D** rg_BzSurface3D::getControlNet() const
{
    if ( ctrlNet == rg_NULL )
    {
        return rg_NULL;
    }
	rg_Point3D** returnNet = new rg_Point3D*[degreeOfU + 1];
	rg_INDEX i=0;
	for( i=0; i <= degreeOfU; i++ )
		returnNet[i] = new rg_Point3D[degreeOfV + 1];

    for ( i=0; i <= degreeOfU; i++ )
	{
		for( rg_INDEX j=0; j <= degreeOfV; j++ )
			returnNet[i][j] = ctrlNet[i][j];
	}
    return returnNet;
}

void rg_BzSurface3D::setDegree( const rg_DEGREE& uDegree, const rg_DEGREE& vDegree )
{
    degreeOfU = uDegree;
	degreeOfV = vDegree;

    if( ctrlNet )
	{	
		for( rg_INDEX j=0; j <= degreeOfV; j++ )
			delete[] ctrlNet[j];

        delete [] ctrlNet;
	}

	ctrlNet = new rg_Point3D*[degreeOfU + 1];

	for( rg_INDEX i=0; i <= degreeOfU; i++ )
		ctrlNet[i] = new rg_Point3D[degreeOfV + 1];
}

void rg_BzSurface3D::setCtrlPt( const rg_INDEX& i,
								    const rg_INDEX& j,
									const rg_Point3D& pt )
{
    ctrlNet[i][j] = pt;
}

void rg_BzSurface3D::setControlNet( const rg_INT& numOfRows,
		                       const rg_INT& numOfCols, 
		                       rg_Point3D**  ctrlPts)
{
    degreeOfU = numOfRows - 1;
	degreeOfV = numOfCols - 1;

    if( ctrlNet )
	{	
		for( rg_INDEX j=0; j <= degreeOfV; j++ )
			delete[] ctrlNet[j];

        delete [] ctrlNet;
	}

	ctrlNet = new rg_Point3D*[degreeOfU + 1];

	for( rg_INDEX i=0; i <= degreeOfU; i++ )  {
		ctrlNet[i] = new rg_Point3D[degreeOfV + 1];
		for( rg_INDEX j=0; j<= degreeOfV; j++)  {
			ctrlNet[i][j] = ctrlPts[i][j];
		}
	}
}

void    rg_BzSurface3D::setSurface(const rg_BzSurface3D& temp)
{
    rg_INT newUDegree=temp.degreeOfU;
    rg_INT newVDegree=temp.degreeOfV;

    rg_BzSurface3D::setDegree(newUDegree,newVDegree);

    for ( rg_INT i=0; i < newUDegree+1; i++ )
    {
        for( rg_INT j=0; j < newVDegree+1; j++ )
        {
            rg_Point3D ctrlPt=temp.getCtrlPt(i,j);
            rg_BzSurface3D::setCtrlPt(i,j,ctrlPt);
        }
    }
}

//  Operations
rg_Point3D rg_BzSurface3D::evaluatePt( const rg_PARAMETER& u,
		                     const rg_PARAMETER& v )

{
    rg_Point3D ptOnSurface;

    for ( rg_INDEX i=0; i <= degreeOfU; i++ )
	{
		for ( rg_INDEX j=0; j <= degreeOfV; j++ )
			ptOnSurface += ctrlNet[i][j]*bernstein(i, degreeOfU, u)*bernstein(j, degreeOfV, v);
	}

    return ptOnSurface;
}

rg_Point3D** rg_BzSurface3D::evaluatePt(const rg_INT& numIsoparaCrvOfU,
		                        const rg_INT& numIsoparaCrvOfV,
								const rg_INT& resolution)
{
	rg_INDEX i=0;
	rg_INDEX j=0;
	
	rg_Point3D** ptOnSurface;
	rg_INT numIsoParamCrv = numIsoparaCrvOfU+numIsoparaCrvOfV;
	ptOnSurface = new rg_Point3D*[numIsoParamCrv];
	for(i=0; i<numIsoParamCrv; i++)  {
		ptOnSurface[i] = new rg_Point3D[resolution];
	}
	
	rg_REAL incrementOfRes = 1.0/(resolution-1);
	rg_REAL incrementOfU   = 1.0/(numIsoparaCrvOfU-1);
	rg_REAL incrementOfV   = 1.0/(numIsoparaCrvOfV-1);

	rg_PARAMETER u = 0.0;
	rg_PARAMETER v = 0.0;

	i = 0;
	for ( u = 0.0; rg_LE(u, 1.0); u+=incrementOfU)  {
		j=0;
		for ( v = 0.0; rg_LE(v, 1.0); v+=incrementOfRes)  {
			ptOnSurface[i][j] = evaluatePt(u, v);
			j++;
		}
		i++;
	}

	u = 0.0;
	v = 0.0;
//	i = 0;
	for ( v = 0.0; rg_LE(v, 1.0); v+=incrementOfV)  {
		j=0;
		for ( u = 0.0; rg_LE(u, 1.0); u+=incrementOfRes)  {
			ptOnSurface[i][j] = evaluatePt(u, v);
			j++;
		}
		i++;
	}
	
	return ptOnSurface;
}
// Hodograph

rg_BzSurface3D rg_BzSurface3D::getUHodograph() const
{
    const rg_INT uDegree=degreeOfU-1;
    const rg_INT vDegree=degreeOfV;
    rg_BzSurface3D output(uDegree,vDegree);

    for( rg_INT j=0; j < vDegree+1; j++ )
    {
        for( rg_INT i=0; i < uDegree+1; i++ )
        {
            rg_Point3D pt = ctrlNet[i+1][j]-ctrlNet[i][j];
            output.setCtrlPt( i, j, pt );
        }
    }
    return output;
}

rg_BzSurface3D rg_BzSurface3D::getVHodograph() const
{
    const rg_INT uDegree=degreeOfU;
    const rg_INT vDegree=degreeOfV-1;
    rg_BzSurface3D output(uDegree,vDegree);

    for( rg_INT i=0; i < uDegree+1; i++ )
    {
        for( rg_INT j=0; j < vDegree+1; j++ )
        {
            output.setCtrlPt( i,j, ctrlNet[i][j+1]-ctrlNet[i][j]);
        }
    }
    return output;
}

/*
rg_Point3D** rg_BzSurface3D::deCasteljau( const rg_PARAMETER& u )
{
    rg_Point3D** b = new rg_Point3D*[degree+1];
    for(rg_INDEX i=0; i<=degree; i++)
	    b[i] = new rg_Point3D[degree+1];

    // range of parameter value : [0, 1] 
    // if not satisfied, exit.
    for(rg_INDEX j=0; j<=degree; j++)
	    b[0][j] = ctrlNet[j];

    for(rg_INDEX r=1; r<=degree; r++)
	    for(rg_INDEX j=0; j<=degree-r; j++)
		    b[r][j] = (1-u)*b[r-1][j] + u*b[r-1][j+1];

    return b;
}
*/
rg_REAL rg_BzSurface3D::bernstein( const rg_INDEX&     i,
						   const rg_DEGREE&    dgr,
						   const rg_PARAMETER& t )
{
    rg_REAL B;

    if( i < 0 || dgr-i < 0 ) 
        return 0;	
    else
    {
	    B = factorial((rg_REAL)dgr) / factorial((rg_REAL)i) / factorial((rg_REAL)(dgr-i));
	    return B * pow(t, i) * pow(1-t, dgr-i);
    }

}

rg_REAL rg_BzSurface3D::factorial( const rg_REAL& n )
{
    if( n ) 
        return n * factorial(n-1);
    else 
        return 1;
}


//*****************************************************************************
//
//    FUNCTION    : powerToBezierSurface
//    DESCRIPTION : 
//                  This function is for conversion between bezier and 
//                  power basis form.
//                                    
//    AUTHOR      : Dong-Gyou Lee
//    START DATE  : 24 Mar. 1998   
//    REFERENCE   : The NURBS Book, Les Piegl & Wayne Tiller, p.277 
//
//*****************************************************************************

void rg_BzSurface3D::powerToBezierSurface( const rg_DEGREE& uDegree,
									  const rg_DEGREE& vDegree, 
									  const rg_REAL paramValuesOfU[],
									  const rg_REAL paramValuesOfV[],
									  const rg_Matrix& powerCoeff )
{
	degreeOfU = uDegree;
	degreeOfV = vDegree;

    if( ctrlNet )
	{	
		for( rg_INDEX j=0; j <= degreeOfV; j++ )
			delete[] ctrlNet[j];

        delete [] ctrlNet;
	}

	ctrlNet = new rg_Point3D*[degreeOfU + 1];

	rg_INT i = 0;
	for( i=0; i <= degreeOfU; i++ )
		ctrlNet[i] = new rg_Point3D[degreeOfV + 1];

	rg_REAL upperParamOfU = paramValuesOfU[1];
	rg_REAL lowerParamOfU = paramValuesOfU[0];
	
	rg_REAL upperParamOfV = paramValuesOfV[1];
	rg_REAL lowerParamOfV = paramValuesOfV[0];

	rg_Matrix Mpi( rg_CurveSurfaceFunc::powerToBezierMatrix( degreeOfU + 1 ) );
	rg_Matrix Rpi( rg_CurveSurfaceFunc::reparameterMatrix( (degreeOfU+1),
		        upperParamOfU-lowerParamOfU,
				lowerParamOfU ) );

	rg_Matrix Mqi( rg_CurveSurfaceFunc::powerToBezierMatrix( degreeOfV + 1 ) );
	rg_Matrix Rqi( rg_CurveSurfaceFunc::reparameterMatrix( (degreeOfV+1),
		        upperParamOfV-lowerParamOfV,
				lowerParamOfV ) );

	rg_Matrix Mqt = Mqi.transpose();
	rg_Matrix Rqt = Rqi.transpose();

	rg_Matrix xMatrix( degreeOfU+1, degreeOfV+1 );
	rg_Matrix yMatrix( degreeOfU+1, degreeOfV+1 );
	rg_Matrix zMatrix( degreeOfU+1, degreeOfV+1 );

	for( i=0; i <= degreeOfU; i++ )
	{
		for( rg_INDEX j=0; j <= degreeOfV; j++)
		{
			xMatrix[i][j] = powerCoeff.getElement( i*(degreeOfV+1)+j, 0 );
			yMatrix[i][j] = powerCoeff.getElement( i*(degreeOfV+1)+j, 1 );
			zMatrix[i][j] = powerCoeff.getElement( i*(degreeOfV+1)+j, 2 );
		}
	}

	rg_Matrix resultOfX = (Mpi * Rpi) * xMatrix * ( Rqt * Mqt );
	rg_Matrix resultOfY = (Mpi * Rpi) * yMatrix * ( Rqt * Mqt );
	rg_Matrix resultOfZ = (Mpi * Rpi) * zMatrix * ( Rqt * Mqt );
    
	for( i=0; i <= degreeOfU; i++ )
	{
		for( rg_INDEX j=0; j <= degreeOfV; j++)
		{
			ctrlNet[i][j].setX( resultOfX[i][j] );
			ctrlNet[i][j].setY( resultOfY[i][j] );
			ctrlNet[i][j].setZ( resultOfZ[i][j] );
		}
	}
}


rg_BzSurface3D& rg_BzSurface3D::operator=(const rg_BzSurface3D& temp)
{
    if ( this == &temp )
    {
        return *this;
    }

    rg_BzSurface3D::setSurface(temp);
    return *this;
}
			




