//********************************************************************
//
//	  FILENAME    : rg_BzCurve3D.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_BzCurve3D
//           which defines a Bezier rg_Curve in 3-D and represents 
//           its properties. 
//                          
//	  CLASS NAME  : rg_BzCurve3D
//
//    BASE CLASS  : rg_Curve
//      
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//
//    HISTORY     : 	
//	        By Dong-Gyou Lee 18 Mar. 1998
//                  void powerToBezierCurve( const rg_DEGREE& dgr,
//                                           const rg_REAL[] paramValues,
//                                           const rg_Matrix& powerCoeff )
//
//    START DATE  : 9 Jul. 1997    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************

#include "rg_ListByPtr.h"
#include "rg_BzCurve3D.h"
#include "rg_CubicPolynomial.h"

#include "rg_CurveSurfaceFunc.h"
#include "rg_MathFunc.h"

//  Constructor & Destructor
rg_BzCurve3D::rg_BzCurve3D()
{
    ctrlPts = rg_NULL;
}

rg_BzCurve3D::rg_BzCurve3D(const rg_DEGREE &dgr)
    : degree(dgr)
{
    ctrlPts = new rg_Point3D[degree + 1];
}

rg_BzCurve3D::rg_BzCurve3D(const rg_DEGREE &dgr, const rg_Point3D* tCtrlPts)
    : degree(dgr)
{
    ctrlPts = new rg_Point3D [degree + 1];
    for (rg_INDEX i = 0; i<=degree; i++)
        ctrlPts[i] = tCtrlPts[i];
}

rg_BzCurve3D::rg_BzCurve3D(const rg_BzCurve3D &curve)
    : degree( curve.degree )
{
    ctrlPts = new rg_Point3D [degree + 1];
    for (rg_INDEX i = 0; i<=degree; i++)
        ctrlPts[i] = curve.ctrlPts[i];
}

rg_BzCurve3D::~rg_BzCurve3D()
{
    if (ctrlPts != rg_NULL)
        delete [] ctrlPts;
}

//  Access elements
rg_DEGREE rg_BzCurve3D::getDegree() const
{
    return degree;
}

rg_DEGREE rg_BzCurve3D::getOrder() const
{
    return degree+1;
}


rg_Point3D rg_BzCurve3D::getCtrlPt(const rg_INDEX &i) const
{
    return ctrlPts[i];
}

rg_Point3D* rg_BzCurve3D::getCtrlPts() const
{
    rg_Point3D* ctrlPolygon = new rg_Point3D[degree + 1];

    for (rg_INDEX i=0; i<=degree; i++)
        ctrlPolygon[i] = ctrlPts[i];

    return ctrlPolygon;
}
rg_BzCurve2D rg_BzCurve3D::evaluateBzCurve2D() const
{
	rg_BzCurve2D output;

	output.setDegree(degree);

	for( rg_INT i=0; i < degree+1; i++ )
	{
		output.setCtrlPt(i, ctrlPts[i].evaluatePt2D() );
	}

	return output;
}

void rg_BzCurve3D::setDegree(const rg_DEGREE &dgr)
{
    degree = dgr;

	if( ctrlPts )
		delete[] ctrlPts;

    ctrlPts = new rg_Point3D[degree + 1];
}

void rg_BzCurve3D::setOrder(const rg_DEGREE &tOrder)
{
    degree = tOrder-1;

	if( ctrlPts )
		delete[] ctrlPts;

    ctrlPts = new rg_Point3D[degree + 1];
}


void rg_BzCurve3D::setCtrlPt(const rg_INDEX &i, const rg_Point3D &pt)
{
    ctrlPts[i] = pt;
}

void rg_BzCurve3D::setCtrlPts(const rg_INT &numOfCtrlPt, const rg_Point3D* tCtrlPts)
{
//    if ( (degree + 1) != numOfCtrlPt )
//        return;

    degree=numOfCtrlPt-1;
    if (ctrlPts != rg_NULL)
        delete [] ctrlPts;

    ctrlPts = new rg_Point3D[degree + 1];
    for (rg_INDEX i = 0; i<=degree; i++)
        ctrlPts[i] = tCtrlPts[i];
}

void rg_BzCurve3D::setCurve(const rg_BzCurve3D &curve)
{
    degree = curve.degree;

    if (ctrlPts != rg_NULL)
        delete [] ctrlPts;

    ctrlPts = new rg_Point3D[degree + 1];
    for (rg_INDEX i = 0; i<=degree; i++)
        ctrlPts[i] = curve.ctrlPts[i];
}

//  Operations
rg_Point3D rg_BzCurve3D::evaluatePt(const rg_PARAMETER &u) const
{
    rg_Point3D ptOnCurve;

    for (rg_INDEX i=0; i<=degree; i++)
        ptOnCurve += ctrlPts[i]*bernstein(i, degree, u);

    return ptOnCurve;
}

rg_BzCurve3D rg_BzCurve3D::makeDerivative()
{
    rg_Point3D* ctrlPts = new rg_Point3D[degree];
    for(rg_INDEX i=0; i<=degree-1; i++)
        ctrlPts[i] = degree*(ctrlPts[i+1]-ctrlPts[i]);
	
    return rg_BzCurve3D(degree-1, ctrlPts);
}

void rg_BzCurve3D::raiseDegree(const rg_DEGREE& raisingTimes)
{
    rg_DEGREE   oldDegree=degree;
    degree=oldDegree+raisingTimes;
    rg_Point3D* oldCtrlPts=ctrlPts;
    ctrlPts= new rg_Point3D[degree+1];

    for( rg_INT i=0; i < degree+1; i++ )
    {
        ctrlPts[i]=rg_Point3D( 0.0, 0.0, 0.0);

        rg_INT minJ=( i-raisingTimes < 0 ) ? 0:i-raisingTimes;
        rg_INT maxJ=( i > oldDegree ) ? oldDegree:i;
        for( rg_INT j= minJ ; j <= maxJ; j++ )
        {
            ctrlPts[i]+= oldCtrlPts[j]
                          *rg_MathFunc::combination(oldDegree,j)
                          *rg_MathFunc::combination(raisingTimes,i-j)
                          /rg_MathFunc::combination(degree,i);
        }
    }
    delete []oldCtrlPts;
}

rg_Point3D** rg_BzCurve3D::deCasteljau(const rg_PARAMETER &u)
{
    rg_Point3D** b = new rg_Point3D*[degree+1];
    for(rg_INDEX i=0; i<=degree; i++)
	    b[i] = new rg_Point3D[degree+1];

    // range of parameter value : [0, 1] 
    // if not satisfied, exit.
    for(rg_INDEX j=0; j<=degree; j++)
	    b[0][j] = ctrlPts[j];

    for(rg_INDEX r=1; r<=degree; r++)
	    for(rg_INDEX j=0; j<=degree-r; j++)
		    b[r][j] = (1-u)*b[r-1][j] + u*b[r-1][j+1];

    return b;
}

rg_Polynomial rg_BzCurve3D::bernstein(const rg_DEGREE &n, const rg_INDEX &i) const
{
	// define t
	rg_REAL    *t    = new rg_REAL[2];
	t[0] = 0.0;
	t[1] = 1.0;

	// define (1-t)
	rg_REAL    *t_1  = new rg_REAL[2];
	t_1[0] =  1.0;
    t_1[1] = -1.0;

	rg_Polynomial t_i(1, t);
	rg_Polynomial t_1_n_i(1, t_1);
    delete[] t;
    delete[] t_1;

	// define the Bernstein polynomial : (n, i) t^i (1-t)^(n-i)
	rg_Polynomial result = rg_MathFunc::combination(n, i)*t_i.power(i)*t_1_n_i.power(n-i);
    
    return result;
}

rg_Polynomial* rg_BzCurve3D::convertBzCurve2Polynomial() const
{
    rg_REAL  *x = new rg_REAL[degree+1];
    rg_REAL  *y = new rg_REAL[degree+1];
	rg_REAL  *z = new rg_REAL[degree+1];
    
	rg_INT i = 0;
	for(i = 0; i < degree+1; i++)
	{
		x[i] = ctrlPts[i].getX();
		y[i] = ctrlPts[i].getY();
		z[i] = ctrlPts[i].getZ();
	}

    rg_Polynomial *ct = new rg_Polynomial[3];	
    ct[0] = x[0]*bernstein(degree, 0);     // x(t)
    ct[1] = y[0]*bernstein(degree, 0);     // y(t)
	ct[2] = y[0]*bernstein(degree, 0);     // z(t)

    for(i = 1; i < degree + 1; i++)
    {
        ct[0] = ct[0] + x[i]*bernstein(degree, i);
        ct[1] = ct[1] + y[i]*bernstein(degree, i);
		ct[2] = ct[2] + z[i]*bernstein(degree, i);
    }

    delete[] x;
    delete[] y;
	delete[] z;

    return ct;
}

rg_REAL rg_BzCurve3D::bernstein(const rg_INDEX &i, const rg_DEGREE &dgr, const rg_PARAMETER &t) const
{
    rg_REAL B;

    if( i < 0 || dgr-i < 0 ) 
        return 0;	
    else
    {
	    B = factorial((rg_REAL)dgr) / factorial((rg_REAL)i) / factorial((rg_REAL)(dgr-i));	    return B * pow(t, i) * pow(1-t, dgr-i);
    }

}

rg_REAL rg_BzCurve3D::factorial(const rg_REAL &n) const
{
    if( n ) 
        return n * factorial(n-1);
    else 
        return 1;
}

//  Intersection 
rg_sListByPtr* rg_BzCurve3D::intersectOfCubicBezierAndPlane(const rg_Plane3D &plane) 
{
    if (degree != rg_CUBIC)
        return rg_NULL;

    //--------------------------------------------------------------
    //  Cubic Bezier curve
    //          - 3    3
    //   B(t) = >     B (t) P 
    //          - i=0  i     i
    //
    //  Implicit form of cubic Bezier curve
    //                                                      
    //   B(t) = (P3 - 3*P2 + 3*P1 - P0)*t^3  + (3*P2 - 6*P1 + 3*P0)*t^2
    //               + (3*P1 - 3*P0)*t + P0
    //
    //  rg_Point3D : Q(u, v) = A + u*B + v*C
    //
    //  Intersection of cubic Bezier curve and plane
    //      B(t) = Q(u, v)
    //
    //    inner product (B * C)
    //
    //      C3*t^3  + C2*t^2  + C1*t + C0 = 0
    //--------------------------------------------------------------

    rg_Point3D coeffOfImplicitForm[4] 
            = { (ctrlPts[0] - plane.getPosVector()), 
                (3*ctrlPts[1] - 3*ctrlPts[0]), 
                (3*ctrlPts[2] - 6*ctrlPts[1] + 3*ctrlPts[0]),
                (ctrlPts[3] - 3*ctrlPts[2] + 3*ctrlPts[1] - ctrlPts[0]) };

    rg_Point3D crsProduct = plane.getUVector()*plane.getVVector();

    rg_REAL* coeff = new rg_REAL[4];
    
    for(rg_INDEX i=0; i<4; i++)
    {
        coeff[i] = coeffOfImplicitForm[i]%crsProduct;
    }

    rg_sListByPtr* intersectPoint = new rg_sListByPtr;
    if ( rg_NZERO( coeff[3] ) )
    {
        rg_Polynomial intersectEqOfCubicBezierPlane(rg_CUBIC);
        intersectEqOfCubicBezierPlane.setCoefficient(coeff);

        rg_ComplexNumber* root = intersectEqOfCubicBezierPlane.solve();
        for (rg_INDEX j=0; j<rg_CUBIC; j++)
        {
            if ( rg_ZERO( root[j].getImaginaryNumber() ) ) 
            {
                rg_REAL realRoot = root[j].getRealNumber();
                if ( rg_BTOR(0.0, realRoot, 1.0) )
                {
					rg_Point3D tPt=evaluatePt(realRoot);
                    rg_Point3D* pt = new rg_Point3D(tPt  );
                    intersectPoint->append(pt);
                }
            }
        }
        delete [] root;
    }
    else if ( rg_NZERO( coeff[2] ) )
    {
        rg_Polynomial intersectEqOfCubicBezierPlane(rg_QUADRATIC);
        for (rg_INT i=0; i<=rg_QUADRATIC; i++)
            intersectEqOfCubicBezierPlane.setCoefficient(i, coeff[i]);

        rg_ComplexNumber* root = intersectEqOfCubicBezierPlane.solve();
        for (rg_INDEX j=0; j<rg_QUADRATIC; j++)
        {
            if ( rg_ZERO( root[j].getImaginaryNumber() ) ) 
            {
                rg_REAL realRoot = root[j].getRealNumber();
                if ( rg_BTOR(0.0, realRoot, 1.0) )
                {
					rg_Point3D tPt=evaluatePt(realRoot);
                    rg_Point3D* pt = new rg_Point3D(tPt);
                    intersectPoint->append(pt);
                }
            }
        }
        delete [] root;
    }
    else if ( rg_NZERO( coeff[1] ) )
    {
        rg_REAL realRoot = -coeff[0]/coeff[1];
        if ( rg_BTOR(0.0, realRoot, 1.0) )
        {
			rg_Point3D tPt=evaluatePt(realRoot);
            rg_Point3D* pt = new rg_Point3D(tPt);
            intersectPoint->append(pt);
/*
            rg_Point3D* pt = new rg_Point3D( evaluatePt(realRoot) );
            intersectPoint->append(pt);
*/
        }
    }
    else ;
//    rg_CubicPolynomial intersectEqOfCubicBezierPlane(coeff);

    delete [] coeff;

    return intersectPoint;
}

rg_dList<rg_Point3D> rg_BzCurve3D::intersectWithPlaneForCubic(const rg_Plane3D& plane) const
{
    rg_dList<rg_Point3D> output;

    if (degree != rg_CUBIC)
    {
        return output;
    }

    //--------------------------------------------------------------
    //  Cubic Bezier curve
    //          - 3    3
    //   B(t) = >     B (t) P 
    //          - i=0  i     i
    //
    //  Implicit form of cubic Bezier curve
    //                                                      
    //   B(t) = (P3 - 3*P2 + 3*P1 - P0)*t^3  + (3*P2 - 6*P1 + 3*P0)*t^2
    //               + (3*P1 - 3*P0)*t + P0
    //
    //  rg_Point3D : Q(u, v) = A + u*B + v*C
    //
    //  Intersection of cubic Bezier curve and plane
    //      B(t) = Q(u, v)
    //
    //    inner product (B * C)
    //
    //      C3*t^3  + C2*t^2  + C1*t + C0 = 0
    //--------------------------------------------------------------

    rg_Point3D coeffOfImplicitForm[4] 
            = { (ctrlPts[0] - plane.getPosVector()), 
                (3*ctrlPts[1] - 3*ctrlPts[0]), 
                (3*ctrlPts[2] - 6*ctrlPts[1] + 3*ctrlPts[0]),
                (ctrlPts[3] - 3*ctrlPts[2] + 3*ctrlPts[1] - ctrlPts[0]) };

    rg_Point3D crsProduct = plane.getUVector()*plane.getVVector();

    rg_REAL* coeff = new rg_REAL[4];
    
    for(rg_INDEX i=0; i<4; i++)
    {
        coeff[i] = coeffOfImplicitForm[i]%crsProduct;
    }

    if ( rg_NZERO( coeff[3] ) )
    {
        rg_Polynomial intersectEqOfCubicBezierPlane(rg_CUBIC);
        intersectEqOfCubicBezierPlane.setCoefficient(coeff);

        rg_ComplexNumber* root = intersectEqOfCubicBezierPlane.solve();
        for (rg_INDEX j=0; j<rg_CUBIC; j++)
        {
            if ( rg_ZERO( root[j].getImaginaryNumber() ) ) 
            {
                rg_REAL realRoot = root[j].getRealNumber();
                if ( rg_BTOR(0.0, realRoot, 1.0) )
                {
                    rg_Point3D pt = evaluatePt(realRoot) ;
                    output.add(pt);
                }
            }
        }
        delete [] root;
    }
    else if ( rg_NZERO( coeff[2] ) )
    {
        rg_Polynomial intersectEqOfCubicBezierPlane(rg_QUADRATIC);
        for (rg_INT i=0; i<=rg_QUADRATIC; i++)
            intersectEqOfCubicBezierPlane.setCoefficient(i, coeff[i]);

        rg_ComplexNumber* root = intersectEqOfCubicBezierPlane.solve();
        for (rg_INDEX j=0; j<rg_QUADRATIC; j++)
        {
            if ( rg_ZERO( root[j].getImaginaryNumber() ) ) 
            {
                rg_REAL realRoot = root[j].getRealNumber();
                if ( rg_BTOR(0.0, realRoot, 1.0) )
                {
                    rg_Point3D pt =  evaluatePt(realRoot);
                    output.add(pt);
                }
            }
        }
        delete [] root;
    }
    else if ( rg_NZERO( coeff[1] ) )
    {
        rg_REAL realRoot = -coeff[0]/coeff[1];
        if ( rg_BTOR(0.0, realRoot, 1.0) )
        {
            rg_Point3D pt =  evaluatePt(realRoot);
            output.add(pt);
        }
    }
    else ;
//    rg_CubicPolynomial intersectEqOfCubicBezierPlane(coeff);

    delete [] coeff;

    return output;

}

//*****************************************************************************
//
//    FUNCTION    : powerToBezierCurve
//    DESCRIPTION : 
//                  This function is for conversion between bezier and 
//                  power basis form.
//                                    
//                  form : Cj(s) = [Bi,p(s)][Pi] = [si]Mp[Pi] 
//                               = [ui]RpMp[Pi]  = [ui][ai]    = Cj(u)
//
//                         So, [Pi] = MiRi[ai] 
//
//    AUTHOR      : Dong-Gyou Lee
//    START DATE  : 18 Mar. 1998   
//    REFERENCE   : The NURBS Book, Les Piegl & Wayne Tiller, p.271 
//
//*****************************************************************************


void rg_BzCurve3D::powerToBezierCurve( const rg_DEGREE& dgr,
		                            const rg_REAL paramValues[],
						            const rg_Matrix& powerCoeff )
{
	degree = dgr;
	
	if( ctrlPts )
		delete[] ctrlPts;

    ctrlPts = new rg_Point3D[degree + 1];

	rg_REAL upperParam = paramValues[1];
	rg_REAL lowerParam = paramValues[0];

	rg_Matrix Mi( rg_CurveSurfaceFunc::powerToBezierMatrix( degree + 1 ) );
	rg_Matrix Ri( rg_CurveSurfaceFunc::reparameterMatrix( (degree+1), (upperParam-lowerParam), lowerParam ) );

	rg_Matrix resultMatrix = (Mi * Ri) * powerCoeff;
    
	for ( rg_INDEX i=0; i <= degree; i++ )
	{
		ctrlPts[i].setX( resultMatrix[i][0] );
		ctrlPts[i].setY( resultMatrix[i][1] );
		ctrlPts[i].setZ( resultMatrix[i][2] );
	}
}
			


