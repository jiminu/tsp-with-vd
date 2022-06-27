//********************************************************************
//
//	  FILENAME    : BSplineCurve.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_BSplineCurve3D 
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun 1996    
//
//    HISTORY     :
//          BY Young-Song Cho.  13 Jul. 1997
//              rg_INT rg_BSplineCurve3D::getNumOfKnotSpan() const
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//*********************************************************************


#include <math.h>

#include "rg_RelativeOp.h"

#include "rg_BSplineCurve3D.h"

////	Constructor & Destructor
rg_BSplineCurve3D::rg_BSplineCurve3D()
: order(0), numOfCtrlPts(0)
{
	ctrlPts = rg_NULL;
}

rg_BSplineCurve3D::rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                                const rg_Planarity    &newPlanarity )
: rg_Curve(newID, newPlanarity), order(0), numOfCtrlPts(0)
{
	ctrlPts = rg_NULL;
}

rg_BSplineCurve3D::rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                                const rg_Planarity    &newPlanarity,
                                const rg_ORDER        &newOrder )
:rg_Curve(newID, newPlanarity), order(newOrder), numOfCtrlPts(0),ctrlPts(rg_NULL) 
{
}

rg_BSplineCurve3D::rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                                const rg_Planarity    &newPlanarity,
                                const rg_INT          &num )
:rg_Curve(newID, newPlanarity), order(4), numOfCtrlPts(num)
{
	ctrlPts = new rg_Point3D [num];
}

rg_BSplineCurve3D::rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                                const rg_Planarity    &newPlanarity,
                                const rg_INT          &num, 
                                rg_Point3D*              newControlP)
:rg_Curve(newID, newPlanarity), order(4), numOfCtrlPts(num)
{
	ctrlPts = new rg_Point3D [num];

	for (rg_INT i=0; i<num; i++)
		ctrlPts[i] = newControlP[i];
}

rg_BSplineCurve3D::rg_BSplineCurve3D( const rg_ORDER &newOrder )
: order(newOrder), numOfCtrlPts(0)
{
	ctrlPts = rg_NULL;
}

rg_BSplineCurve3D::rg_BSplineCurve3D( const rg_INT &num )
: order(4), numOfCtrlPts(num)
{
    if ( numOfCtrlPts > 0 )
    {
    	ctrlPts = new rg_Point3D [num];
    }
    else
    {
        ctrlPts=rg_NULL;
    }
}

rg_BSplineCurve3D::rg_BSplineCurve3D( const rg_INT &num, 
                                rg_Point3D*     newControlP )
: order(4), numOfCtrlPts(num)
{
	ctrlPts = new rg_Point3D [num];

	for (rg_INT i=0; i<num; i++)
		ctrlPts[i] = newControlP[i];
}

rg_BSplineCurve3D::rg_BSplineCurve3D( const rg_ORDER &newOrder, 
                                const rg_INT   &num, 
                                rg_Point3D*       newControlP )
: order(newOrder), numOfCtrlPts(num)
{
	ctrlPts = new rg_Point3D [num];

	for (rg_INT i=0; i<num; i++)
		ctrlPts[i] = newControlP[i];
}

////  Constructor      : March 13 1997
rg_BSplineCurve3D::rg_BSplineCurve3D( const unsigned rg_INT &newID, 
                                const rg_Planarity    &newPlanarity,
                                const rg_ORDER        &newOrder,
                                const rg_INT          &num, 
                                rg_Point3D*             newControlP )
: rg_Curve(newID, newPlanarity), order(newOrder), numOfCtrlPts(num)
{
    ctrlPts = new rg_Point3D[numOfCtrlPts];

    for(rg_INT i=0; i<numOfCtrlPts; i++)
        ctrlPts[i] = newControlP[i];
}

////  Copy Constructor March 13 1997
rg_BSplineCurve3D::rg_BSplineCurve3D( const rg_BSplineCurve3D &curve)
: rg_Curve(curve.getID(), curve.getPlanarity()), 
  order(curve.order)
{
    numOfCtrlPts = curve.numOfCtrlPts;

	if ( numOfCtrlPts > 0 ) 
	{
		ctrlPts = new rg_Point3D[numOfCtrlPts];
	    for(rg_INT i=0; i<numOfCtrlPts; i++)
		    ctrlPts[i] = curve.ctrlPts[i];
	}
	else
	{
		ctrlPts=rg_NULL;
	}
	
}

rg_BSplineCurve3D::~rg_BSplineCurve3D()
{
	if (ctrlPts != rg_NULL)
		delete [] ctrlPts;
}

////	Get Functions.
rg_ORDER rg_BSplineCurve3D::getOrder() const
{
	return order;
}

rg_INT rg_BSplineCurve3D::getNumOfCtrlPts() const
{
	return numOfCtrlPts;
}

rg_Point3D rg_BSplineCurve3D::getCtrlPt( const rg_INT &i ) const
{
	return ctrlPts[i];
}

rg_Point3D* rg_BSplineCurve3D::getCtrlPts() const
{
	if ( isNull() )
	{
		return rg_NULL;
	}

    rg_Point3D* ctrlplygn = new rg_Point3D[numOfCtrlPts];

    for(rg_INT i=0; i<numOfCtrlPts; i++)
        ctrlplygn[i] = ctrlPts[i];

    return ctrlplygn;
}

rg_Point3D*      rg_BSplineCurve3D::accessCtrlPts() const
{
    return ctrlPts;
}

rg_REAL rg_BSplineCurve3D::getKnotValue( const rg_INT &i ) const
{
	rg_INT n      = numOfCtrlPts;
	rg_INT iOrder = (rg_INT) order;

	if ( i < (n + iOrder) )
	{
		if ( i < iOrder )
			return  0.;
		else if ( i < n )
			return (i - iOrder + 1.) / (n - iOrder + 1.);
		else
			return 1.;
	}
	else
		return -1.;
}

rg_INT rg_BSplineCurve3D::getNumOfKnotSpan() const
{
    return numOfCtrlPts + order - 2 * (order-1) - 1;
}

rg_Point3D  rg_BSplineCurve3D::getStartPoint() const
{
    return ctrlPts[0];
}

rg_Point3D  rg_BSplineCurve3D::getEndPoint() const
{
    const rg_INT lastIndex=getNumOfCtrlPts() -1;
    return ctrlPts[lastIndex];
}

////	Set Functions.
void rg_BSplineCurve3D::setOrder( const rg_ORDER &newOrder )
{
	order = newOrder;
}

void rg_BSplineCurve3D::setNumOfCtrlPts( const rg_INT &numOfCtrlPt )
{
	if (ctrlPts != rg_NULL)
		delete [] ctrlPts;

	numOfCtrlPts = numOfCtrlPt;

	ctrlPts = new rg_Point3D [numOfCtrlPt];
}

void rg_BSplineCurve3D::setCtrlPt( const rg_INT   &i, 
										 const rg_Point3D &newPt )
{
	ctrlPts[i] = newPt;
}

void rg_BSplineCurve3D::setCtrlPts( const rg_INT &numOfCtrlPt, 
								    rg_Point3D*     newCtrlPlygn )
{
	if (ctrlPts != rg_NULL)
		delete [] ctrlPts;

	numOfCtrlPts = numOfCtrlPt;

	if ( numOfCtrlPts > 0 )
	{
		ctrlPts = new rg_Point3D [numOfCtrlPt];

		for (rg_INT i=0; i<numOfCtrlPt; i++)
			ctrlPts[i] = newCtrlPlygn[i];
	}
	else
	{
		ctrlPts=rg_NULL;
	}
}

////	Evaluating Function
rg_REAL rg_BSplineCurve3D::evaluateBasisFunc( const rg_INT  &index, 
                                        const rg_REAL &param,
                                        const rg_INT  &Order ) const
{
	rg_INT n      = numOfCtrlPts-1;
//    if (Order == -1)
//        Order = (rg_INT) order;
	
	if (    ( index == 0 && rg_EQ(param, getKnotValue(0)) )
		 || ( index == n && rg_EQ(param, getKnotValue(n+Order)) ) )
	{
		return 1.0;
	}
	else if (    rg_LT(param, getKnotValue(index))
		      || rg_GE(param, getKnotValue(index + Order)) )
	{
		return 0.0;
	}
	else
	{}

	rg_REAL Uleft  = 0.;
	rg_REAL Uright = 0.;
	rg_REAL temp   = 0.;
	rg_REAL saved  = 0.;

	rg_REAL* triN = new rg_REAL [Order];
	rg_INT i = 0;
	for (i=0; i<Order; i++)
	{
		if (    rg_GE(param, getKnotValue(index + i))
			 && rg_LT(param, getKnotValue(index + i + 1)) )
		{
			triN[i] = 1.0;
		}
		else 
			triN[i] = 0.0;
	}

	for (i=1; i<Order; i++)
	{
		if (triN[0] == 0.0) 
			saved = 0.0;
		else 
		{
			saved = ( (param - getKnotValue(index)) * triN[0] ) 
				    / ( getKnotValue(index+i) - getKnotValue(index) );
		}
		
		for (rg_INT j=0; j<(Order-i); j++)
		{
			Uleft  = getKnotValue(index + j + 1);
			Uright = getKnotValue(index + j + i + 1);

			if (triN[j+1] == 0.0)
			{
				triN[j] = saved;
				saved   = 0.0;
			}
			else
			{
				temp    = triN[j+1] / (Uright - Uleft);
				triN[j] = saved + (Uright - param)*temp;
				saved   = (param - Uleft) * temp;
			}
		}
	}
	
	rg_REAL returnValue = triN[0];

	delete [] triN;

	return returnValue;
}

//  April  3 1997 : Modified
////////////////////////////////////////////////////////////////// 
rg_Point3D rg_BSplineCurve3D::evaluatePt( const rg_REAL &u ) const
{

    rg_Point3D ptOnCurve;

    rg_REAL basisValue = 0.;
	for (rg_INT i=0; i<numOfCtrlPts; i++)
	{
		basisValue = evaluateBasisFunc(i, u, order);

        if ( rg_NZERO(basisValue) )
            ptOnCurve += ctrlPts[i]*basisValue;
	}

	return ptOnCurve;
}

//  April  3 1997 : Modified
////////////////////////////////////////////////////////////////// 
rg_Point3D* rg_BSplineCurve3D::evaluatePtsInEvenParameter( const rg_INT &numfPtOnCurve ) const
{
	rg_Point3D* evaluatedPoint = new rg_Point3D [numfPtOnCurve];

	rg_REAL increment = 1. / (numfPtOnCurve-1);
	rg_REAL u         = 0.;

	for (rg_INT i=0; i<(numfPtOnCurve); i++)
	{
		evaluatedPoint[i] = evaluatePt(u);
		u += increment;
	}

	return evaluatedPoint;
}

rg_Polyline2D rg_BSplineCurve3D::makePolyline2DInEvenParameter(const rg_INT &numOfPts) const
{
    if ( numOfPts < 1 )
    {
        return rg_Polyline2D();
    }
    rg_Point3D* tPts=evaluatePtsInEvenParameter(numOfPts);
    rg_Point2D*  pts=new rg_Point2D[numOfPts];
    for( rg_INT i=0; i < numOfPts; i++ )
    {
        pts[i]=tPts[i].evaluatePt2D();
    }
    
    rg_Polyline2D output(numOfPts, pts);
    delete[] pts;
    delete[] tPts;

    return output;
}

rg_Polyline3D rg_BSplineCurve3D::makePolyline3DInEvenParameter(const rg_INT &numOfPts) const
{
    if ( numOfPts < 1 )
    {
        return rg_Polyline3D();
    }
    rg_Point3D* pts=evaluatePtsInEvenParameter(numOfPts);
    
    rg_Polyline3D output(numOfPts, pts);
    delete[] pts;

    return output;
}

bool    rg_BSplineCurve3D::isPlanarCurve(rg_Point3D& normal) const
{
	if ( numOfCtrlPts<3 )
		return rg_TRUE;

	rg_Point3D vec1;
	rg_Point3D vec2;
	rg_Point3D nVec;

	vec1 = ctrlPts[1] - ctrlPts[0];
	vec2 = ctrlPts[2] - ctrlPts[0];
	nVec = vec1.crossProduct(vec2);
	rg_REAL d = nVec.getX()*ctrlPts[0].getX()
		       + nVec.getY()*ctrlPts[0].getY()
			   + nVec.getZ()*ctrlPts[0].getZ();

	for (rg_INT i=3; i<numOfCtrlPts; i++)
	{
		rg_REAL exist;
		exist = nVec.getX()*ctrlPts[i].getX()
		        + nVec.getY()*ctrlPts[i].getY()
			    + nVec.getZ()*ctrlPts[i].getZ()
		        - d;
		if ( rg_NZERO(exist, resNeg5) )
		{
			return rg_FALSE;	
		}
 	}

	normal = nVec;

	return rg_TRUE;
}

void    rg_BSplineCurve3D::transform(const rg_TMatrix3D& transform)
{
    if ( isNull() == rg_TRUE )
    {
        return;
    }

    for( rg_INT i=0; i < numOfCtrlPts; i++ )
    {
        ctrlPts[i]=transform*ctrlPts[i];
    }
}

rg_Point3D  rg_BSplineCurve3D::evaluateCenterOfCtrlPts() const
{
    rg_Point3D output(0.0,0.0,0.0);

    if ( isNull() != rg_TRUE )
    {
        for( rg_INT i=0; i < numOfCtrlPts; i++ )
        {
            output+=ctrlPts[i];
        }
        output=output/(rg_REAL)numOfCtrlPts;
    }

    return output;
}

void   rg_BSplineCurve3D::reverseTrace()
{
    rg_INT halfNumOfCtrlPts=numOfCtrlPts/2;
    for( rg_INT i=0; i < halfNumOfCtrlPts; i++ )
    {
        const rg_INT before=i;
        const rg_INT after=numOfCtrlPts-1-i;
        const rg_Point3D temp=ctrlPts[before];
        ctrlPts[before]=ctrlPts[after];
        ctrlPts[after]=temp;
    }
}



rg_REAL     rg_BSplineCurve3D::evaluateCurveLengthByCtrlPts() const
{
    rg_REAL curveLengthByCtrlPts = 0.0;
    for (rg_INT i = 0; i < numOfCtrlPts - 1; ++i) {
        curveLengthByCtrlPts += ctrlPts[i].distance(ctrlPts[i + 1]);
    }

    return curveLengthByCtrlPts;
}

////	Derivative of B-Spline rg_Curve

//  April 7 1997 : Modified.
//      Because of copy constructor and = operator overloading
rg_BSplineCurve3D rg_BSplineCurve3D::makeDerivative()
{
	rg_REAL degreeDividedByKnotDelta = 0.;

	rg_BSplineCurve3D derivative(numOfCtrlPts-1);

    //  Set the order Of derivative curve.
    derivative.order = order -1;

    //  Set control polygon of derivative curve.
    rg_Point3D newControlPoint;
	for(rg_INT i=0; i<numOfCtrlPts-1; i++)
	{
		degreeDividedByKnotDelta = (order - 1)
			                       / (getKnotValue(i+order) - getKnotValue(i+1));

		newControlPoint = degreeDividedByKnotDelta
			              * (ctrlPts[i+1] - ctrlPts[i]);

		derivative.setCtrlPt(i, newControlPoint);
	}

	return derivative;
}

//--------------------------------------------------------------------
//	Given parameter,u, solve the curvature in XY-plane.
rg_REAL rg_BSplineCurve3D::getCurvatureInXY( const rg_REAL &u )
{
	rg_BSplineCurve3D firstDeriv  = makeDerivative();
	rg_BSplineCurve3D secondDeriv = firstDeriv.makeDerivative();

	rg_Point3D ptOnFirstDeri  = firstDeriv.evaluatePt(u);
	rg_Point3D ptOnSecondDeri = secondDeriv.evaluatePt(u);

	rg_REAL signedCurvature = 0.;
	rg_REAL xPrime          = ptOnFirstDeri.getX();
	rg_REAL yPrime          = ptOnFirstDeri.getY();
	rg_REAL xTwoPrime       = ptOnSecondDeri.getX();  
	rg_REAL yTwoPrime       = ptOnSecondDeri.getY();  

	signedCurvature = ( xTwoPrime * yPrime - yTwoPrime * xPrime )
		              / pow( (xPrime*xPrime + yPrime*yPrime), 1.5 );
	
	return signedCurvature;
}

//--------------------------------------------------------------------
//	Given parameter,u, solve the curvature in XZ-plane.
rg_REAL rg_BSplineCurve3D::getCurvatureInXZ( const rg_REAL &u )
{
	rg_BSplineCurve3D firstDeriv  = makeDerivative();
	rg_BSplineCurve3D secondDeriv = firstDeriv.makeDerivative();

	rg_Point3D ptOnFirstDeri  = firstDeriv.evaluatePt(u);
	rg_Point3D ptOnSecondDeri = secondDeriv.evaluatePt(u);

	rg_REAL signedCurvature = 0.;
	rg_REAL xPrime          = ptOnFirstDeri.getX();
	rg_REAL zPrime          = ptOnFirstDeri.getZ();
	rg_REAL xTwoPrime       = ptOnSecondDeri.getX();  
	rg_REAL zTwoPrime       = ptOnSecondDeri.getZ();  

	signedCurvature = ( xTwoPrime * zPrime - zTwoPrime * xPrime )
		              / pow( (xPrime*xPrime + zPrime*zPrime), 1.5 );
	
	return signedCurvature;
}

//--------------------------------------------------------------------
//	Given parameter,u, solve the curvature in YZ-plane.
rg_REAL rg_BSplineCurve3D::getCurvatureInYZ( const rg_REAL &u )
{
	rg_BSplineCurve3D firstDeriv  = makeDerivative();
	rg_BSplineCurve3D secondDeriv = firstDeriv.makeDerivative();

	rg_Point3D ptOnFirstDeri  = firstDeriv.evaluatePt(u);
	rg_Point3D ptOnSecondDeri = secondDeriv.evaluatePt(u);

	rg_REAL signedCurvature = 0.;
	rg_REAL yPrime          = ptOnFirstDeri.getY();
	rg_REAL zPrime          = ptOnFirstDeri.getZ();
	rg_REAL yTwoPrime       = ptOnSecondDeri.getY();  
	rg_REAL zTwoPrime       = ptOnSecondDeri.getZ();  

	signedCurvature = ( yTwoPrime * zPrime - zTwoPrime * yPrime )
		              / pow( (yPrime*yPrime + zPrime*zPrime), 1.5 );
	
	return signedCurvature;
}



////    Operator Overloading

//  April 7 1997 : made.
rg_BSplineCurve3D& rg_BSplineCurve3D::operator =(const rg_BSplineCurve3D &curve)
{
    if ( this == &curve )
        return *this;

    rg_Curve::setID(curve.getID());
    rg_Curve::setPlanarity(curve.getPlanarity());

    order = curve.order;

    numOfCtrlPts = curve.numOfCtrlPts;

    if (ctrlPts != rg_NULL)
        delete [] ctrlPts;

    ctrlPts = new rg_Point3D[numOfCtrlPts];

    for(rg_INT i=0; i<numOfCtrlPts; i++)
        ctrlPts[i] = curve.ctrlPts[i];

    return *this;
}

rg_FLAG rg_BSplineCurve3D::isNull() const
{
    if ( ctrlPts == NULL )
    {
        return rg_TRUE;
    }
    else
    {
        return rg_FALSE;
    }
}

rg_INT   rg_BSplineCurve3D::findIndexOfNearestCtrlPt(const rg_Point3D& pt) const
{

    if ( numOfCtrlPts < 1 )
    {
        return -1;
    }
    

    rg_REAL  min=(pt-(this->getCtrlPt(0))).magnitude();
    rg_INT   indexOfMin=0;
    for( rg_INT i=1; i < numOfCtrlPts; i++)
    {
        rg_REAL candidate=pt.distance(this->getCtrlPt(i));
        if ( rg_GT(min, candidate) )
        {
            min=candidate;
            indexOfMin=i;
        }
    }
    return indexOfMin;
}

rg_BoundingBox3D rg_BSplineCurve3D::makeBoundingBox3D() const
{
	rg_BoundingBox3D box;
	box.contain(numOfCtrlPts,ctrlPts);
	return box;
}

rg_BoundingBox2D rg_BSplineCurve3D::makeBoundingBox2D() const
{
	rg_BoundingBox2D box;
    for( rg_INT i=0; i < numOfCtrlPts; i++)
    {
	    box.contain(ctrlPts[i].evaluatePt2D());
    }
	return box;
}


void   rg_BSplineCurve3D::removeAll()
{
    order=0;
    if ( numOfCtrlPts > 0 )
    {
        delete[] ctrlPts;
        ctrlPts=rg_NULL;
    }
    numOfCtrlPts=0;
    ctrlPts=rg_NULL;
    rg_Curve::removeAll();
}

