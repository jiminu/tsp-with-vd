//********************************************************************
//
//	  FILENAME    : NURBSplineSurface.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_NURBSplineSurface3D 
//
//    AUTHOR      : Young-Song Cho
//    START DATE  : 21 Jun 1996    
//
//            Copyright (c) CAD/CAM Lab.    	  	
//
//********************************************************************


#include <math.h>
#include "rg_RelativeOp.h"
#include "rg_CurveSurfaceFunc.h"

#include "rg_NUBSplineCurve3D.h"

#include "rg_NURBSplineSurface3D.h"
#include "rg_GeoFunc.h"
#include "rg_MathFunc.h"

////    Constructor & Destructor.
rg_NURBSplineSurface3D::rg_NURBSplineSurface3D()
                         : rg_NUBSplineSurface3D()
{
    weight_vectors = rg_NULL;
}

rg_NURBSplineSurface3D::rg_NURBSplineSurface3D( const rg_INT   &row, 
                                          const rg_INT   &col, 
                                          const rg_ORDER &uOrder, 
                                          const rg_ORDER &vOrder )
: rg_NUBSplineSurface3D( row, col, uOrder, vOrder )
{
    weight_vectors = rg_NULL;
}

rg_NURBSplineSurface3D::rg_NURBSplineSurface3D(const unsigned rg_INT &newID, 
                                         const rg_Planarity    &newPlanarity, 
                                         const rg_INT          &row, 
                                         const rg_INT          &col, 
                                         const rg_ORDER        &uOrder, 
                                         const rg_ORDER        &vOrder,
                                         rg_Point3D**             ctrlNet,
                                         rg_REAL*               uKnotVector,
                                         rg_REAL*               vKnotVector,
                                         rg_REAL**              weightVector)
:  rg_NUBSplineSurface3D( newID,
                       newPlanarity,
                       row,
                       col,
                       uOrder,
                       vOrder,
                       ctrlNet,
                       uKnotVector,
                       vKnotVector )
{
    weight_vectors = new rg_REAL*[row];
    for (rg_INT i=0; i<row; i++)
    {
        weight_vectors[i] = new rg_REAL[col];
        for (rg_INT j=0; j<col; j++)
            weight_vectors[i][j] = weightVector[i][j];
    }
}

rg_NURBSplineSurface3D::rg_NURBSplineSurface3D(const rg_NURBSplineSurface3D &surface )
:  rg_NUBSplineSurface3D( surface.getID(),
                       surface.getPlanarity(),
                       surface.getRowOfControlNet(),
                       surface.getColumnOfControlNet(),
                       surface.getOrderOfU(),
                       surface.getOrderOfV(),
                       surface.getControlNet(),
                       surface.getKnotVectorOfU(),
                       surface.getKnotVectorOfV() )
{
    rg_INT row = surface.getRowOfControlNet();
    rg_INT col = surface.getColumnOfControlNet();

    weight_vectors = new rg_REAL*[row];
    for (rg_INT i=0; i<row; i++)
    {
        weight_vectors[i] = new rg_REAL[col];
        for (rg_INT j=0; j<col; j++)
            weight_vectors[i][j] = surface.getWeight(i,j);
    }
}    

rg_NURBSplineSurface3D::~rg_NURBSplineSurface3D()
{
    rg_INT row = getRowOfControlNet();

    if (weight_vectors != rg_NULL)
    {
        for (rg_INT i=0; i<row; i++)
            delete [] weight_vectors[i];
        delete [] weight_vectors;
    }
}

////    Get Functions.
rg_REAL rg_NURBSplineSurface3D::getWeight(const rg_INDEX &i, const rg_INDEX &j) const
{
    return weight_vectors[i][j];
}

rg_REAL** rg_NURBSplineSurface3D::getWeightVector() const
{
    return weight_vectors;
}

rg_NURBSplineCurve3D rg_NURBSplineSurface3D::getUIsoparamCurve(const rg_PARAMETER &u)
{
    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

    rg_ORDER uOrder = rg_BSplineSurface3D::getOrderOfU();
    rg_ORDER vOrder = rg_BSplineSurface3D::getOrderOfV();

    rg_Point3D* ctrlPtOfCurve = new rg_Point3D[col];
    rg_REAL*  weightOfCurve = new rg_REAL[col];

    for (rg_INDEX j=0; j<col; j++)
    {
        ctrlPtOfCurve[j] = 0.0;
        weightOfCurve[j] = 0.0;
        for (rg_INDEX i=0; i<row; i++)
        {
            ctrlPtOfCurve[j] += rg_NUBSplineSurface3D::evaluateBasisFuncU(i, u, uOrder)
                                * weight_vectors[i][j] 
                                * getPointOnControlNet(i, j);
            weightOfCurve[j] += rg_NUBSplineSurface3D::evaluateBasisFuncU(i, u, uOrder)
                                * weight_vectors[i][j];    
        }
        ctrlPtOfCurve[j] = ctrlPtOfCurve[j]/weightOfCurve[j];
    }

    rg_NURBSplineCurve3D isoparametricCurve(vOrder);
    isoparametricCurve.setCtrlPts(col, ctrlPtOfCurve);
    isoparametricCurve.setKnotVector(
                         vOrder + col,
                         rg_NUBSplineSurface3D::getKnotVectorOfV() );
    isoparametricCurve.setWeightVector( weightOfCurve );

    delete [] ctrlPtOfCurve;
    delete [] weightOfCurve;

    return isoparametricCurve;
}

rg_NURBSplineCurve3D rg_NURBSplineSurface3D::getVIsoparamCurve(const rg_PARAMETER &v)
{
    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

    rg_ORDER uOrder = rg_BSplineSurface3D::getOrderOfU();
    rg_ORDER vOrder = rg_BSplineSurface3D::getOrderOfV();

    rg_Point3D* ctrlPtOfCurve = new rg_Point3D[row];
    rg_REAL*  weightOfCurve = new rg_REAL[row];

    for (rg_INDEX i=0; i<row; i++)
    {
        ctrlPtOfCurve[i] = 0.0;
        weightOfCurve[i] = 0.0;
        for (rg_INDEX j=0; j<col; j++)
        {
            ctrlPtOfCurve[i] += rg_NUBSplineSurface3D::evaluateBasisFuncV(j, v, vOrder)
                                * weight_vectors[i][j] 
                                * getPointOnControlNet(i, j);
            weightOfCurve[i] += rg_NUBSplineSurface3D::evaluateBasisFuncV(j, v, vOrder)
                                * weight_vectors[i][j];    
        }
        ctrlPtOfCurve[i] = ctrlPtOfCurve[i]/weightOfCurve[i];
    }

    rg_NURBSplineCurve3D isoparametricCurve(uOrder);
    isoparametricCurve.setCtrlPts(row, ctrlPtOfCurve);
    isoparametricCurve.setKnotVector(
                         uOrder + row, 
                         rg_NUBSplineSurface3D::getKnotVectorOfU() );
    isoparametricCurve.setWeightVector( weightOfCurve );

    delete [] ctrlPtOfCurve;
    delete [] weightOfCurve;

    return isoparametricCurve;
}

void rg_NURBSplineSurface3D::getPiecewiseSurfaceInPowerForm(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_REAL**** & polyCoeffOfWeight) const
{
	rg_INT p = getOrderOfU() - 1;                     // degree in U-dir
	rg_INT q = getOrderOfV() - 1;                     // degree in V-dir
	rg_INT n = getRowOfControlNet() - 1;              // No. of rows of control net - 1
	rg_INT m = getColumnOfControlNet() - 1;           // No. of column of control net - 1

	rg_INT r = n + p + 1;                             // no. of knotvector in U-dir - 1
	rg_INT s = m + q + 1;                             // no. of knotvector in V-dir - 1

	rg_INT*** allPossiblePathInU = new rg_INT** [p + 1]; // all possible paths from (p + 1) graph primitives in U-dir
	rg_INT*** allPossiblePathInV = new rg_INT** [q + 1]; // all possible paths from (q + 1) graph primitives in V-dir
	rg_INT NoOfNonZeroLengthKnotSpanInU = r - 2 * p;
	rg_INT NoOfNonZeroLengthKnotSpanInV = s - 2 * q;

	// find all the possible paths from (p + 1) and (q + 1) graph primitives

	rg_INT* noOfAllPossiblePathsInEachGraphInU = new rg_INT[p + 1];
	rg_INT* noOfAllPossiblePathsInEachGraphInV = new rg_INT[q + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		allPossiblePathInU[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, p - i);
		noOfAllPossiblePathsInEachGraphInU[ i ] = rg_MathFunc::combination(p, i);
	}
	for(i = 0;i < q + 1;i++)
	{
		allPossiblePathInV[ i ] = rg_MathFunc::enumerateZeroOneSequenceRevised(i, q - i);
		noOfAllPossiblePathsInEachGraphInV[ i ] = rg_MathFunc::combination(q, i);
	}

	//------------------------------------------

	// initialization of coefficient of piecewise polyomial surface

	//rg_INT noOfPiecewiseSurface = NoOfNonZeroLengthKnotSpanInU * NoOfNonZeroLengthKnotSpanInV;

	polyCoeffOfX = new rg_REAL*** [NoOfNonZeroLengthKnotSpanInU];
	polyCoeffOfY = new rg_REAL*** [NoOfNonZeroLengthKnotSpanInU];
	polyCoeffOfZ = new rg_REAL*** [NoOfNonZeroLengthKnotSpanInU];
	polyCoeffOfWeight = new rg_REAL*** [NoOfNonZeroLengthKnotSpanInU];

//	rg_REAL** tempPolyCoeff;

	for(i = 0;i < NoOfNonZeroLengthKnotSpanInU;i++)
	{
		polyCoeffOfX[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];
		polyCoeffOfY[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];
		polyCoeffOfZ[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];
		polyCoeffOfWeight[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];

		for(rg_INT j = 0;j < NoOfNonZeroLengthKnotSpanInV;j++)
		{
			polyCoeffOfX[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfY[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfZ[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfWeight[ i ][ j ] = new rg_REAL* [p + 1];

			for(rg_INT k = 0;k < p + 1;k++)
			{
				polyCoeffOfX[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfY[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfZ[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfWeight[ i ][ j ][ k ] = new rg_REAL [q + 1];

				for(rg_INT l = 0;l < q + 1;l++)
				{
					polyCoeffOfX[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfY[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfZ[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfWeight[ i ][ j ][ k ][ l ] = 0.0;
				}

			}
		}
	}

	
	// generating the piecewise surface in a power form

	rg_Point3D** ctrlNet = getControlNet();
	rg_INT indexOfPolySurfaceInU = 0;  // index for tracing surface subpatch in U-dir
	rg_INT indexOfPolySurfaceInV = 0;  // index for tracing surface subpatch in V-dir
	//rg_INT indexOFBasisFuncInU = 0;    // index for tracing B-spline basis function in U-dir
	//rg_INT indexOFBasisFuncInV = 0;    // index for tracing B-spline basis function in V-dir

	for(i = p;i <= r - p - 1;i++)
	{
		for(rg_INT j = q;j <= s - q - 1;j++)
		{
			for(rg_INT k = i - p;k <= i;k++)
			{
				for(rg_INT l = j - q;l <= j;l++)
				{
					rg_REAL** tempPolyCoeff = getDistributionPolynomialInOneGraphPrimitive(i ,j ,allPossiblePathInU[k - i + p] ,allPossiblePathInV[l - j + q] ,noOfAllPossiblePathsInEachGraphInU[k - i + p] ,noOfAllPossiblePathsInEachGraphInV[l - j + q]);

					// nested loop for computing multiplication tensor product basis function with control points
					for(rg_INT indexOFBasisFuncInU = 0;indexOFBasisFuncInU < p + 1;indexOFBasisFuncInU++)
					{
						for(rg_INT indexOFBasisFuncInV = 0;indexOFBasisFuncInV < q + 1;indexOFBasisFuncInV++)
						{
							polyCoeffOfX[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (weight_vectors[ k ][ l ] * ctrlNet[ k ][ l ].getX() * tempPolyCoeff[indexOFBasisFuncInU][indexOFBasisFuncInV]);
							polyCoeffOfY[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (weight_vectors[ k ][ l ] * ctrlNet[ k ][ l ].getY() * tempPolyCoeff[indexOFBasisFuncInU][indexOFBasisFuncInV]);
							polyCoeffOfZ[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (weight_vectors[ k ][ l ] * ctrlNet[ k ][ l ].getZ() * tempPolyCoeff[indexOFBasisFuncInU][indexOFBasisFuncInV]);
							polyCoeffOfWeight[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (weight_vectors[ k ][ l ] * tempPolyCoeff[indexOFBasisFuncInU][indexOFBasisFuncInV]);
						}
					}
					//------------------------------------------------------------------------------------------

					for(rg_INT indexU = 0;indexU < p;indexU++)
					{
						delete[] tempPolyCoeff[indexU];
					}
					delete[] tempPolyCoeff;
				}
			}
			indexOfPolySurfaceInV++;
		}
		indexOfPolySurfaceInU++;
	}

}


////    Set Functions.
void rg_NURBSplineSurface3D::setControlNetAndWeight(const rg_INDEX& tRow,
                                                 const rg_INDEX& tCol )
{
    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet(); 

    if (weight_vectors != rg_NULL)
    {
        for (rg_INT i=0; i<row; i++)
            delete [] weight_vectors[i];
        delete [] weight_vectors;
    }

    rg_BSplineSurface3D::setControlNet(tRow,tCol);

    weight_vectors = new rg_REAL*[tRow];
    for (rg_INT i=0; i<tRow; i++)
    {
        weight_vectors[i] = new rg_REAL[tCol];
    }
}


void rg_NURBSplineSurface3D::setWeight(const rg_INDEX &i, const rg_INDEX &j, const rg_REAL &weight)
{
    weight_vectors[i][j] = weight;
}

void rg_NURBSplineSurface3D::setWeightVector(rg_REAL** weightVector)
{

    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet(); 
    if (weight_vectors != rg_NULL)
    {
        for (rg_INT i=0; i<row; i++)
            delete [] weight_vectors[i];
        delete [] weight_vectors;
    }

    weight_vectors = new rg_REAL*[row];
    for (rg_INT i=0; i<row; i++)
    {
        weight_vectors[i] = new rg_REAL[col];
        for (rg_INT j=0; j<col; j++)
            weight_vectors[i][j] = weightVector[i][j];
    }
}



void rg_NURBSplineSurface3D::setSurface(const rg_RBzSurface3D& surface)
{
    const rg_INT uOrder=surface.getDegreeOfU()+1;
    const rg_INT vOrder=surface.getDegreeOfV()+1;
    rg_BSplineSurface3D::setOrderOfU(uOrder);
    rg_BSplineSurface3D::setOrderOfV(vOrder);

    rg_Point3D** ctrlPts=surface.getControlNet();
    rg_REAL**  weights=surface.getWeightVector();

    rg_BSplineSurface3D::setControlNet(uOrder,vOrder,ctrlPts);
    rg_NURBSplineSurface3D::setWeightVector(weights);
    rg_INT i = 0;
    for( i=0; i < uOrder; i++ )
    {
        delete[] ctrlPts[i];
        delete[] weights[i];
    }

    delete[] ctrlPts;
    delete[] weights;
    
    const rg_INT numOfUKnots=2*uOrder;
    rg_REAL* uKnots=new rg_REAL[numOfUKnots];
    for ( i=0; i <uOrder; i++ )
    {
        uKnots[i]=0.0;
        uKnots[uOrder+i]=1.0;
    }

    const rg_INT numOfVKnots=2*vOrder;
    rg_REAL* vKnots=new rg_REAL[numOfVKnots];
    for ( i=0; i <vOrder; i++ )
    {
        vKnots[i]=0.0;
        vKnots[uOrder+i]=1.0;
    }

    rg_NUBSplineSurface3D::setKnotVectorOfU(numOfUKnots,uKnots);
    rg_NUBSplineSurface3D::setKnotVectorOfV(numOfVKnots,vKnots);
    
    delete[] uKnots;
    delete[] vKnots;
    
}

////    Operating & Calculating.
//
//  April  3 1997 : Modified
////////////////////////////////////////////////////////////////// 
rg_Point3D rg_NURBSplineSurface3D::evaluatePt(const rg_PARAMETER &u, const rg_PARAMETER &v)
{
    if (rg_NUBSplineSurface3D::getKnotVectorOfU() == rg_NULL)
	{
        rg_NUBSplineSurface3D::setInitialKnotVectorOfU();
	}
	if (rg_NUBSplineSurface3D::getKnotVectorOfV() == rg_NULL)
	{
		rg_NUBSplineSurface3D::setInitialKnotVectorOfV();
	}
	
	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();

    rg_INDEX validRow = 0;
    rg_INDEX validCol = 0;

    rg_INT first = uOrder-1;
    rg_INT last  = rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT middle;

    while ( last >= first )
    {
        middle = (first + last)/2;

        if ( rg_LT(u, rg_NUBSplineSurface3D::getKnotValueOfU(middle)) )
            last = middle;
        else if ( rg_GT(u, rg_NUBSplineSurface3D::getKnotValueOfU(middle + 1)) )
            first = middle + 1;
        else
        {
			if ( middle != first )
			{
	            while ( rg_EQ(getKnotValueOfU(middle), getKnotValueOfU(middle+1)) )
		            middle++;
	            validRow = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
			else
			{
	            validRow = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
			break;
        }
    }

    first = vOrder-1;
    last  = rg_BSplineSurface3D::getColumnOfControlNet();
    while ( last >= first )
    {
        middle = (first + last)/2;

        if ( rg_LT(v, rg_NUBSplineSurface3D::getKnotValueOfV(middle)) )
            last = middle;
        else if ( rg_GT(v, rg_NUBSplineSurface3D::getKnotValueOfV(middle + 1)) )
            first = middle + 1;
        else
        {
			if ( middle != first )
			{
	            while ( rg_EQ(getKnotValueOfV(middle), getKnotValueOfV(middle+1)) )
		            middle++;
	            validCol = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
			else
			{
	            validCol = middle; // modify : By Young-Song Cho  14 Aug. 1997
			}
            break;
        }
    }

    rg_REAL* nonZeroBasisU = evaluateMultiBasisFuncU( validRow, u, uOrder );
    rg_REAL* nonZeroBasisV = evaluateMultiBasisFuncV( validCol, v, vOrder );

	rg_Point3D ptOnSurface;
    rg_REAL  rationalPart = 0.0;
    for (rg_INDEX i=validRow-uOrder+1; i<=validRow; i++)
    {
        for (rg_INDEX j=validCol-vOrder+1; j<=validCol; j++)
        {
            rationalPart += nonZeroBasisU[i - validRow + uOrder - 1]
                            * nonZeroBasisV[j - validCol + vOrder - 1]
                            * weight_vectors[i][j];

            ptOnSurface += nonZeroBasisU[i - validRow + uOrder - 1]
                           * nonZeroBasisV[j - validCol + vOrder - 1]
                           * weight_vectors[i][j]
                           * rg_BSplineSurface3D::getPointOnControlNet(i,j);               
        }
    }

    delete [] nonZeroBasisU;
    delete [] nonZeroBasisV;

    return ptOnSurface/rationalPart;
/*
    for	(rg_INT i=0; i<row; i++)
	{
		uBasisValue = evaluateBasisFuncU(i, u, uOrder);

        if ( rg_NZERO(uBasisValue) ) 
        {
    		for (rg_INT j=0; j<col; j++)
	    	{
                vBasisValue = evaluateBasisFuncV(j, v, vOrder);
                if ( rg_NZERO(vBasisValue) )
                {
                    rationalPart += uBasisValue * vBasisValue * weight_vectors[i][j]; 

                    ptOnSurface += rg_BSplineSurface3D::getPointOnControlNet(i,j)
                               * uBasisValue * vBasisValue * weight_vectors[i][j];
                }
            }
		}
	}

	return ptOnSurface/rationalPart;
*/
}

////    Derivative.
rg_Point3D rg_NURBSplineSurface3D::derivativeSurfaceOfU(const rg_PARAMETER &u, const rg_PARAMETER &v) 
{
//
//           A(u,v)                Au(u,v) - wu(u,v)S(u,v)
// S(u,v) = -------- ,  Su(u,v) = -------------------------
//           w(u,v)                        w(u,v)
//
//____________________________________________________________

    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();

    //  To obtain the derivatie basis function in the u direction.
    rg_NUBSplineCurve3D forUBasisFunc(row-1);
    forUBasisFunc.rg_BSplineCurve3D::setOrder(uOrder-1);

    rg_REAL* tempKnot = new rg_REAL [row + uOrder - 2];
    rg_INT i = 0;
    for (i=0; i<(row+uOrder-2); i++)
        tempKnot[i] = rg_NUBSplineSurface3D::getKnotValueOfU(i+1);

    forUBasisFunc.rg_NUBSplineCurve3D::setKnotVector(row+uOrder-2, tempKnot);
    delete [] tempKnot;

    rg_Point3D AprimeU;        //  Au(u, v)
    rg_Point3D SU;             //  S(u, v)
    rg_REAL  wprimeU = 0.0;  //  wu(u, v)
    rg_REAL  wU      = 0.0;  //  w(u, v)

    rg_REAL  uBasisFunc = 0.0;
    rg_REAL  vBasisFunc = 0.0;
    rg_REAL  uKnotDelta = 0.0;

    SU = evaluatePt(u, v);

    for (i=0; i<row-1; i++)
    {
        uBasisFunc = forUBasisFunc.rg_NUBSplineCurve3D::evaluateBasisFunc(i, u, uOrder-1);
        uKnotDelta = rg_NUBSplineSurface3D::getKnotValueOfU(i+uOrder)
                     - rg_NUBSplineSurface3D::getKnotValueOfU(i+1);

        if ( rg_NZERO(uBasisFunc) )
        {
            for (rg_INT j=0; j<col; j++)
            {
                vBasisFunc = evaluateBasisFuncV(j, v, vOrder);
                if ( rg_NZERO(vBasisFunc) )
                {
                    AprimeU += uBasisFunc*vBasisFunc
                               * (weight_vectors[i+1][j]*getPointOnControlNet(i+1, j)
                                 - weight_vectors[i][j]*getPointOnControlNet(i, j) )
                               / uKnotDelta;
                    wprimeU += uBasisFunc*vBasisFunc
                               * (weight_vectors[i+1][j] - weight_vectors[i][j]);
                }
            }
        }
    }
    AprimeU = (uOrder-1)*AprimeU;
    wprimeU = (uOrder-1)*wprimeU;

    for (i=0; i<row; i++)
    {
        uBasisFunc = evaluateBasisFuncU(i, u, uOrder);
        if ( rg_NZERO(uBasisFunc) )
        {
            for (rg_INT j=0; j<col; j++)
                wU += uBasisFunc
                      * evaluateBasisFuncV(j, v, vOrder)
                      * weight_vectors[i][j];
        }
    }

    rg_Point3D derivative = (AprimeU - wprimeU*SU) / wU;

    return derivative;
}

rg_Point3D rg_NURBSplineSurface3D::derivativeSurfaceOfV(const rg_PARAMETER &u, const rg_PARAMETER &v) 
{
    //
    //           A(u,v)                Av(u,v) - wv(u,v)S(u,v)
    // S(u,v) = -------- ,  Sv(u,v) = -------------------------
    //           w(u,v)                        w(u,v)
    //
    //____________________________________________________________
    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();

    //  To obtain the derivatie basis function in the v direction.
    rg_NUBSplineCurve3D forVBasisFunc(col-1);
    forVBasisFunc.rg_BSplineCurve3D::setOrder(vOrder-1);

    rg_REAL* tempKnot = new rg_REAL [col + vOrder - 2];
    rg_INT i = 0;
    for (i=0; i<(col+vOrder-2); i++)
        tempKnot[i] = rg_NUBSplineSurface3D::getKnotValueOfV(i+1);

    forVBasisFunc.rg_NUBSplineCurve3D::setKnotVector(col+vOrder-2, tempKnot);
    delete [] tempKnot;

    rg_Point3D AprimeU;        //  Av(u, v)
    rg_Point3D SU;             //  S(u, v)
    rg_REAL  wprimeU = 0.0;  //  wv(u, v)
    rg_REAL  wU      = 0.0;  //  w(u, v)

    rg_REAL  uBasisFunc = 0.0;
    rg_REAL  vBasisFunc = 0.0;
    rg_REAL  vKnotDelta = 0.0;

    SU = evaluatePt(u, v);

//    for (i=0; i<row-1; i++)
    for (rg_INT j=0; j<col-1; j++)
    {
        vBasisFunc = forVBasisFunc.rg_NUBSplineCurve3D::evaluateBasisFunc(j, v, vOrder-1);
        vKnotDelta = rg_NUBSplineSurface3D::getKnotValueOfV(i+vOrder)
                     - rg_NUBSplineSurface3D::getKnotValueOfV(i+1);

        if ( rg_NZERO(vBasisFunc) )
        {
            for (i=0; i<row; i++)       
//          for (rg_INT j=0; j<col; j++)
            {
                uBasisFunc = evaluateBasisFuncU(i, u, uOrder);
                if ( rg_NZERO(uBasisFunc) )
                {
                    AprimeU += uBasisFunc*vBasisFunc
                               * (weight_vectors[i][j+1]*getPointOnControlNet(i, j+1)
                                 - weight_vectors[i][j]*getPointOnControlNet(i, j) )
                               / vKnotDelta;
                    wprimeU += uBasisFunc*vBasisFunc
                               * (weight_vectors[i][j+1] - weight_vectors[i][j]);
                }
            }
        }
    }
    AprimeU = (vOrder-1)*AprimeU;
    wprimeU = (vOrder-1)*wprimeU;

    for (i=0; i<row; i++)
    {
        uBasisFunc = evaluateBasisFuncU(i, u, uOrder);
        if ( rg_NZERO(uBasisFunc) )
        {
            for (rg_INT j=0; j<col; j++)
                wU += uBasisFunc
                      * evaluateBasisFuncV(j, v, vOrder)
                      * weight_vectors[i][j];
        }
    }

    rg_Point3D derivative = (AprimeU - wprimeU*SU) / wU;

    return derivative;
    
}
/*
rg_Point3D rg_NURBSplineSurface3D::derivativeSurfaceOfUV(const rg_REAL &u, const rg_REAL &v) 
{
    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();
    
}
*/
////    rg_Surface Construction Technique.
void rg_NURBSplineSurface3D::skinnedSurface(
                              const rg_INT                     &n,
                              const rg_NURBSplineCurve3D* const sectionCurves,
                              rg_REAL*                          parameterization)
{
    //  1. Section curves must have the same order, knotVector.
    ////////////////////////////////////////////////////////////

    //  2. v-directional curves are interpolated.
    /////////////////////////////////////////////
    rg_NURBSplineCurve3D* vInterpolatingCurves;

    //  vInterpolatingCurves is generated with new operator.
    vInterpolatingCurves = vDirectionalCurveInterpolation4Skinning(
                                                       n, 
                                                       sectionCurves,
                                                       parameterization);
    
    //  3. set the skinned surface.
    ///////////////////////////////
    rg_INT nU = sectionCurves[0].rg_BSplineCurve3D::getNumOfCtrlPts(); 
    rg_INT nV = vInterpolatingCurves[0].rg_BSplineCurve3D::getNumOfCtrlPts(); 
    
    //      3.1 set the orders of u and v direction.
    ////////////////////////////////////////////////
    rg_INT uOrder = (rg_INT)(sectionCurves[0].rg_BSplineCurve3D::getOrder());
    rg_INT vOrder = (rg_INT)(vInterpolatingCurves[0].rg_BSplineCurve3D::getOrder());
    
    rg_BSplineSurface3D::setOrderOfU( uOrder );  //  sectionCurves[0].rg_BSplineCurve3D::getOrder() ); 
    rg_BSplineSurface3D::setOrderOfV( vOrder );  //  vInterpolatingCurves[0].rg_BSplineCurve3D::getOrder() ); 

    //      3.2 set the control net.
    ////////////////////////////////
    rg_Point3D** newCtrlNt = new rg_Point3D*[nU];
    rg_INT i = 0;
    for (i=0; i<nU; i++)
    {
        newCtrlNt[i] = new rg_Point3D[nV];
        for (rg_INT j=0; j<nV; j++)
        {
            newCtrlNt[i][j] = vInterpolatingCurves[i].rg_BSplineCurve3D::getCtrlPt(j); 
        }
    }
    rg_BSplineSurface3D::setControlNet(nU, nV, newCtrlNt);

    //      3.3 set the knot vector of u and v direction.
    /////////////////////////////////////////////////////
    rg_NUBSplineSurface3D::setKnotVectorOfU(
                           nU+uOrder, 
                           sectionCurves[0].rg_NUBSplineCurve3D::getKnotVector() );

    rg_NUBSplineSurface3D::setKnotVectorOfV(
                           nV+vOrder, 
                           vInterpolatingCurves[0].rg_NUBSplineCurve3D::getKnotVector() );
    
    //      3.4 set the weight.
    ///////////////////////////
    weight_vectors = new rg_REAL* [nU];
    for (i=0; i<nU; i++)
    {
        weight_vectors[i] = new rg_REAL [nV];
        weight_vectors[i][0] = sectionCurves[0].getWeight(i)
                               * vInterpolatingCurves[i].getWeight(0);
        for (rg_INT j=1; j<nV-1; j++)
        {
            weight_vectors[i][j] = sectionCurves[j-1].getWeight(i)
                                   * vInterpolatingCurves[i].getWeight(j);
        }
        weight_vectors[i][nV-1] = sectionCurves[nV-3].getWeight(i)
                               * vInterpolatingCurves[i].getWeight(nV-1);
    }         
}
     
rg_NURBSplineCurve3D* rg_NURBSplineSurface3D::vDirectionalCurveInterpolation4Skinning(
                        const rg_INT                     &n,
                        const rg_NURBSplineCurve3D* const sectionCurves,
                        rg_REAL*                          parameterization )
{
    rg_INT numOfCtrlPtOfVCurve = sectionCurves[0].rg_BSplineCurve3D::getNumOfCtrlPts(); 

    rg_NURBSplineCurve3D* vDirectionalCurve 
                            = new rg_NURBSplineCurve3D[numOfCtrlPtOfVCurve];

    //  construct the point sets to pass through. 
    rg_Point3D** ptsPassed = new rg_Point3D*[numOfCtrlPtOfVCurve];
    rg_INT i = 0;
    for(i=0; i<numOfCtrlPtOfVCurve; i++)
    {
        ptsPassed[i] = new rg_Point3D[n];
        for (rg_INT j=0; j<n; j++)
            ptsPassed[i][j] = sectionCurves[j].getCtrlPt(i);
    }

    if ( parameterization == rg_NULL )
    {
        //  parameterize the point sets.
        rg_REAL** paramOfVCurve = new rg_REAL*[numOfCtrlPtOfVCurve];
        for (i=0; i<numOfCtrlPtOfVCurve; i++)
            paramOfVCurve[i] = rg_GeoFunc::chordLength(n, ptsPassed[i]);

        //  obtain the parameters to use curve interpolation. 
        parameterization = new rg_REAL [n];
        for (rg_INT j=0; j<n; j++)
        {
            rg_REAL totalParam = 0.0;
            for (i=0; i<numOfCtrlPtOfVCurve; i++)
                totalParam += paramOfVCurve[i][j];
            parameterization[j] = totalParam/numOfCtrlPtOfVCurve;
        }
        for (i=0; i<numOfCtrlPtOfVCurve; i++)
            delete [] paramOfVCurve[i];
        delete [] paramOfVCurve;
    }

    //  v-directional curves interpolation for skinning.
    for (i=0; i<numOfCtrlPtOfVCurve; i++)
        vDirectionalCurve[i].curveInterpolation(n, ptsPassed[i], 4, parameterization);

    for (i=0; i<numOfCtrlPtOfVCurve; i++)
        delete [] ptsPassed[i];
    delete [] ptsPassed;
    
    return vDirectionalCurve;
}

////    Operator Overloading.
//  April 7 1997 : made.
rg_NURBSplineSurface3D& rg_NURBSplineSurface3D::operator =(const rg_NURBSplineSurface3D &surface)
{
    if (this == &surface)
        return *this;

    rg_Surface::setID( surface.rg_Surface::getID() );
    rg_Surface::setPlanarity( surface.rg_Surface::getPlanarity() );

    if (weight_vectors != rg_NULL)
    {
        rg_INT n = rg_BSplineSurface3D::getRowOfControlNet();
        for (rg_INT i=0; i<n; i++)
            delete [] weight_vectors[i];
        delete [] weight_vectors;
    }

    rg_INT nU = surface.rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT nV = surface.rg_BSplineSurface3D::getColumnOfControlNet();

    rg_INT uOrder = (rg_INT) surface.rg_BSplineSurface3D::getOrderOfU();
    rg_INT vOrder = (rg_INT) surface.rg_BSplineSurface3D::getOrderOfV();

    rg_BSplineSurface3D::setControlNet(nU, nV,
                                    surface.rg_BSplineSurface3D::getControlNet());

    rg_BSplineSurface3D::setOrderOfSurface(uOrder, vOrder);

    rg_NUBSplineSurface3D::setKnotVectorOfU( nU+uOrder,
                        surface.rg_NUBSplineSurface3D::getKnotVectorOfU());
    rg_NUBSplineSurface3D::setKnotVectorOfV( nV+vOrder,
                        surface.rg_NUBSplineSurface3D::getKnotVectorOfV());

    weight_vectors = new rg_REAL* [nU];
    for (rg_INT i=0; i<nU; i++)
    {
        weight_vectors[i] = new rg_REAL[nV];
        for (rg_INT j=0; j<nV; j++)
            weight_vectors[i][j] = surface.weight_vectors[i][j];
    }

    return *this;
}

//*****************************************************************************
//
//    FUNCTION    : makeTabulatedCylinder
//    DESCRIPTION : 
//                  Representation of special surface using NURBS surface form.	 
//                  A tabulated cylinder is a surface formed by moving a line 
//                  segment the generatrix parallel to itself along a curve 
//                  called the directrix.
//    FORM        :
//                   X(u,v) = CX(u) + v * (LX - CX(0)) 
//                   Y(u,v) = CY(u) + v * (LY - CY(0)) 
//                   Z(u,v) = CZ(u) + v * (LZ - CZ(0)) 
//                     0 <= u <= 1, 0 <= v <= 1 
//
//    AUTHOR      : Dong-Gyou Lee
//    START DATE  : 8 Apr. 1998   
//    REFERENCE   : IGES ver. 4.0, p.107
//
//*****************************************************************************
rg_FLAG rg_NURBSplineSurface3D::makeTabulatedCylinder( const rg_NURBSplineCurve3D& directionalCurve,
												 const rg_Point3D& endPoint )
{
	rg_INT newUOrder = directionalCurve.getOrder();
	rg_INT newRows   = directionalCurve.getNumOfCtrlPts();
	rg_INT newUKnots = newUOrder + newRows;

//	v-direction is a line!!!
	rg_INT newVOrder = 2;
	rg_INT newCols   = 2;
	rg_INT newVKnots = newVOrder + newCols;

	//	set new control net.
	rg_Point3D** newControlNet = new rg_Point3D*[newRows];
	
	rg_INDEX i = 0;
	for( i=0; i < newRows; i++ )
		newControlNet[i] = new rg_Point3D[newCols];

	for( i=0; i < newRows; i++ )
	{
		newControlNet[i][0] = directionalCurve.getCtrlPt(i);
		newControlNet[i][1] = directionalCurve.getCtrlPt(i) 
			                + ( endPoint - directionalCurve.getCtrlPt(0) );
	}

	//	set net knot vectors.
	rg_REAL* newUKnotVector = new rg_REAL[newUKnots];

	for( i=0; i < newUKnots; i++ )
		newUKnotVector[i] = directionalCurve.getKnotValue(i);

	rg_REAL* newVKnotVector = new rg_REAL[newVKnots];
	
	newVKnotVector[0] = 0.;
	newVKnotVector[1] = 0.;
	newVKnotVector[2] = 1.;
	newVKnotVector[3] = 1.;

	//	set new weights.
	rg_REAL** newWeight = new rg_REAL*[newRows];

	for( i=0; i < newRows; i++ )
		newWeight[i] = new rg_REAL[newCols];

	for( i=0; i < newRows; i++ )
	{
		for( rg_INDEX j=0; j < newCols; j++ )
			newWeight[i][j] = 1.0;
	}

	//	set surface!!!!
	rg_BSplineSurface3D::setOrderOfSurface( newUOrder, newVOrder );
	rg_BSplineSurface3D::setControlNet( newRows, newCols, newControlNet );
	rg_NUBSplineSurface3D::setKnotVectorOfU( newUKnots, newUKnotVector );
	rg_NUBSplineSurface3D::setKnotVectorOfV( newVKnots, newVKnotVector );
	rg_NURBSplineSurface3D::setWeightVector( newWeight );

	//	delete data.
	for( i=0; i < newRows; i++ )
	{
		delete[] newWeight[i];
		delete[] newControlNet[i];
	}

	delete[] newWeight;
	delete[] newVKnotVector;
	delete[] newUKnotVector;
	delete[] newControlNet;

	return rg_TRUE;
}



