//********************************************************************
//
//	  FILENAME    : NUBSplineSurface.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_NUBSplineSurface3D 
//
//    AUTHOR      : Deok-Soo Kim, Young-Song Cho
//    START DATE  : 21 Jun 1996    
//
//    HISTORY     :
//          BY Young-Song Cho.  14 Aug. 1997
//            make : rg_REAL* evaluateMultiBasisFuncU( const rg_INDEX     &uKnotIndex,
//                                                  const rg_PARAMETER &u,
//                                                  const rg_ORDER     &uOrder)
//            make : rg_REAL* evaluateMultiBasisFuncV( const rg_INDEX     &vKnotIndex,
//                                                  const rg_PARAMETER &v,
//                                                  const rg_ORDER     &vOrder)
//            modify : rg_Point3D evaluatePt( const rg_PARAMETER &u, 
//                                       const rg_PARAMETER &v )
//          BY Young-Song Cho.   14 Oct. 1997
//            make : rg_REAL  getGaussianCurvature( const rg_PARAMETER &u, 
//                                               const rg_PARAMETER &v)
//
//          By Dong-Gyou Lee 24 Mar. 1998
//                	void powerSplineToNUBSplineSurface( const rg_DEGREE& uDegree,
//													    const rg_DEGREE& vDegree,
//														const rg_INT& numOfUPatch,
//														const rg_INT& numOfVPatch,
//														rg_REAL** paramValuesOfU,
//														rg_REAL** paramValuesOfV,
//														rg_Matrix** powerCoeff )
//
//                  void bzSurfacesToNUBSplineSurface( const rg_INT& numOfUPatch,
//													   const rg_INT& numOfVPatch,
//													   rg_BzSurface3D** bzSurfaces )
//           By Dong-Gyou Lee 26 Mar. 1998
//                  void reparameterizationKnotVector() 
//
//            Copyright (c) CAD/CAM Lab. Department of Industrial Engineering   
//                          Hanyang University, Seoul Korea	  	
//
//*********************************************************************

#include <math.h>
#include <stdio.h>
#include "rg_RelativeOp.h"
#include "rg_CurveSurfaceFunc.h"

#include "rg_BzSurface3D.h"
#include "rg_NUBSplineSurface3D.h"
#include "rg_GeoFunc.h"
#include "rg_MathFunc.h"

////	Constructor & Destructor.----------------------------------------------
rg_NUBSplineSurface3D::rg_NUBSplineSurface3D()
: rg_BSplineSurface3D()
{
	u_knot_vector = rg_NULL;
	v_knot_vector = rg_NULL;
}

rg_NUBSplineSurface3D::rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
                                        const rg_Planarity    &newPlanarity )
: rg_BSplineSurface3D(newID, newPlanarity)
{
	u_knot_vector = rg_NULL;
	v_knot_vector = rg_NULL;
}

rg_NUBSplineSurface3D::rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
                                        const rg_Planarity    &newPlanarity, 
                                        const rg_INT          &row, 
                                        const rg_INT          &col )
: rg_BSplineSurface3D(newID, newPlanarity, row, col)
{
	u_knot_vector = rg_NULL;
	v_knot_vector = rg_NULL;
}

rg_NUBSplineSurface3D::rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
                                        const rg_Planarity    &newPlanarity, 
                                        const rg_ORDER        &uOrder, 
                                        const rg_ORDER        &vOrder )
: rg_BSplineSurface3D(newID, newPlanarity, uOrder, vOrder)
{
	u_knot_vector = rg_NULL;
	v_knot_vector = rg_NULL;
}

rg_NUBSplineSurface3D::rg_NUBSplineSurface3D( const rg_INT &row, 
                                        const rg_INT &col )
: rg_BSplineSurface3D(row, col)
{
	u_knot_vector = rg_NULL;
	v_knot_vector = rg_NULL;
}

rg_NUBSplineSurface3D::rg_NUBSplineSurface3D( const rg_ORDER &uOrder, 
                                        const rg_ORDER &vOrder )
: rg_BSplineSurface3D(uOrder, vOrder)
{
	u_knot_vector = rg_NULL;
	v_knot_vector = rg_NULL;
}

rg_NUBSplineSurface3D::rg_NUBSplineSurface3D( const rg_INT          &row, 
                                        const rg_INT          &col,
                                        const rg_ORDER &uOrder, 
                                        const rg_ORDER &vOrder )
: rg_BSplineSurface3D(row, col, uOrder, vOrder)
{
	u_knot_vector = rg_NULL;
	v_knot_vector = rg_NULL;
}

rg_NUBSplineSurface3D::rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
                                        const rg_Planarity    &newPlanarity, 
                                        const rg_INT          &row, 
                                        const rg_INT          &col, 
                                        const rg_ORDER        &uOrder, 
                                        const rg_ORDER        &vOrder )
: rg_BSplineSurface3D(newID, newPlanarity, row, col, uOrder, vOrder)
{
	u_knot_vector = rg_NULL;
	v_knot_vector = rg_NULL;
}

//// Constructor     : March 13 1997
rg_NUBSplineSurface3D::rg_NUBSplineSurface3D( const unsigned rg_INT &newID, 
                                        const rg_Planarity    &newPlanarity, 
                                        const rg_INT          &row, 
                                        const rg_INT          &col, 
                                        const rg_ORDER        &uOrder, 
                                        const rg_ORDER        &vOrder,
                                        rg_Point3D**             ctrlNet,
                                        rg_REAL*               uKnotVector,
                                        rg_REAL*               vKnotVector )
: rg_BSplineSurface3D( newID,
                    newPlanarity,
                    row,
                    col,
                    uOrder,
                    vOrder,
                    ctrlNet )
{
    rg_INT nU = (rg_INT) (row + uOrder);
    rg_INT nV = (rg_INT) (col + vOrder);

    u_knot_vector = new rg_REAL [nU];
    v_knot_vector = new rg_REAL [nV];

    rg_INT i = 0;
	for(i=0; i<nU; i++)
        u_knot_vector[i] = uKnotVector[i];

    for( i=0; i<nV; i++)
        v_knot_vector[i] = vKnotVector[i];
}

//// Copy Constructor     : March 13 1997
rg_NUBSplineSurface3D::rg_NUBSplineSurface3D(const rg_NUBSplineSurface3D &surface)
:  rg_BSplineSurface3D( surface.getID(),
                     surface.getPlanarity(),
                     surface.getRowOfControlNet(),
                     surface.getColumnOfControlNet(),
                     surface.getOrderOfU(),
                     surface.getOrderOfV(),
                     surface.getControlNet() )
{
    rg_INT nU = (rg_INT) (surface.getRowOfControlNet() + surface.getOrderOfU());
    rg_INT nV = (rg_INT) (surface.getColumnOfControlNet() + surface.getOrderOfV());

    u_knot_vector = new rg_REAL [nU];
    v_knot_vector = new rg_REAL [nV];

    rg_INT i = 0;
	for(i=0; i<nU; i++)
        u_knot_vector[i] = surface.u_knot_vector[i];

    for( i=0; i<nV; i++)
        v_knot_vector[i] = surface.v_knot_vector[i];
}

rg_NUBSplineSurface3D::~rg_NUBSplineSurface3D()
{
	if (u_knot_vector != rg_NULL)
		delete [] u_knot_vector;

	if (v_knot_vector != rg_NULL)
		delete [] v_knot_vector;
}

////	Get Functions.---------------------------------------------------------
rg_REAL rg_NUBSplineSurface3D::getKnotValueOfU( const rg_INDEX &kIndex ) const
{
	return u_knot_vector[kIndex];
}

rg_REAL* rg_NUBSplineSurface3D::getKnotVectorOfU() const
{
    return u_knot_vector;
}

rg_REAL rg_NUBSplineSurface3D::getKnotValueOfV( const rg_INDEX &kIndex ) const
{
	return v_knot_vector[kIndex];
}

rg_REAL* rg_NUBSplineSurface3D::getKnotVectorOfV() const
{
    return v_knot_vector;
}

rg_INT rg_NUBSplineSurface3D::getNumOfNonZeroLengthKnotSpanOfKnotVectorU() const
{
	rg_INT NumOfNonZeroLengthKnotSpan = 0;
	rg_INT n = getRowOfControlNet() - 1;
	rg_INT p = getOrderOfU() - 1;
	rg_INT i = p;
	rg_INT r = n + p + 1;

	while(i <= r - p - 1)
	{
		if( rg_NE(u_knot_vector[i], u_knot_vector[i+1]) )
			NumOfNonZeroLengthKnotSpan++;
		i++;
	}
	return NumOfNonZeroLengthKnotSpan;
}

rg_INT rg_NUBSplineSurface3D::getNumOfNonZeroLengthKnotSpanOfKnotVectorV() const
{
	rg_INT NumOfNonZeroLengthKnotSpan = 0;
	rg_INT m = getColumnOfControlNet() - 1;
	rg_INT q = getOrderOfV() - 1;
	rg_INT i = q;
	rg_INT s = m + q + 1;

	while(i <= s - q - 1)
	{
		if( rg_NE(v_knot_vector[i], v_knot_vector[i+1]) )
			NumOfNonZeroLengthKnotSpan++;
		i++;
	}
	return NumOfNonZeroLengthKnotSpan;
}

rg_REAL* rg_NUBSplineSurface3D::getDistinctKnotValuesU() const
{
	rg_INT n = getRowOfControlNet() - 1;
	rg_INT p = getOrderOfU() - 1;
	rg_INT i = p;
	rg_INT r = n + p + 1;

	rg_INT NumOfNonZeroLengthKnotSpan = getNumOfNonZeroLengthKnotSpanOfKnotVectorU();
	rg_REAL* distinctKnotValue = new rg_REAL[NumOfNonZeroLengthKnotSpan + 1];

	rg_INT index = 0;
	while((index < NumOfNonZeroLengthKnotSpan + 1) && (i <= r))
	{
		if(rg_NE(u_knot_vector[ i ], u_knot_vector[i + 1]))
		{
			distinctKnotValue[index] = u_knot_vector[ i ];
			index++;
		}
		i++;
	}
	return distinctKnotValue;
}

rg_REAL* rg_NUBSplineSurface3D::getDistinctKnotValuesV() const
{
	rg_INT m = getColumnOfControlNet() - 1;
	rg_INT q = getOrderOfV() - 1;
	rg_INT i = q;
	rg_INT s = m + q + 1;

	rg_INT NumOfNonZeroLengthKnotSpan = getNumOfNonZeroLengthKnotSpanOfKnotVectorV();
	rg_REAL* distinctKnotValue = new rg_REAL[NumOfNonZeroLengthKnotSpan + 1];

	rg_INT index = 0;
	while((index < NumOfNonZeroLengthKnotSpan + 1) && (i <= s))
	{
		if(rg_NE(v_knot_vector[ i ], v_knot_vector[i + 1]))
		{
			distinctKnotValue[index] = v_knot_vector[ i ];
			index++;
		}
		i++;
	}
	return distinctKnotValue;
}

rg_NUBSplineCurve3D rg_NUBSplineSurface3D::getUIsoparamCurve(const rg_PARAMETER &u) 
{
    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

    rg_ORDER uOrder = rg_BSplineSurface3D::getOrderOfU();
    rg_ORDER vOrder = rg_BSplineSurface3D::getOrderOfV();

    rg_Point3D* ctrlPtOfCurve = new rg_Point3D[col];

    rg_INDEX validRow = 0;
    rg_INT   first    = uOrder-1;
    rg_INT   last     = row;
    rg_INT   middle;

    while ( last >= first )
    {
        middle = (first + last)/2;

        if ( rg_LT(u, u_knot_vector[middle]) )
            last = middle;
        else if ( rg_GT(u, u_knot_vector[middle + 1]) )
            first = middle;
        else
        {
//            validRow = middle - uOrder + 1;
            validRow = middle; // modify : By Young-Song Cho  14 Aug. 1997
            break;
        }
    }

    for (rg_INDEX j=0; j<col; j++)
    {
        ctrlPtOfCurve[j] = 0.;

        rg_REAL* nonZeroBasisU = evaluateMultiBasisFuncU( validRow, u, uOrder );

        for (rg_INDEX i=validRow-uOrder+1; i<=validRow; i++)
        {
            ctrlPtOfCurve[j] += rg_BSplineSurface3D::getPointOnControlNet(i,j) 
                                * nonZeroBasisU[i - validRow + uOrder - 1];
        }

        delete [] nonZeroBasisU;
    }

    rg_NUBSplineCurve3D isoparametricCurve(vOrder, col, ctrlPtOfCurve, v_knot_vector);

    delete [] ctrlPtOfCurve;
    return isoparametricCurve;
}

rg_NUBSplineCurve3D rg_NUBSplineSurface3D::getVIsoparamCurve(const rg_PARAMETER &v) 
{
    rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

    rg_ORDER uOrder = rg_BSplineSurface3D::getOrderOfU();
    rg_ORDER vOrder = rg_BSplineSurface3D::getOrderOfV();

    rg_Point3D* ctrlPtOfCurve = new rg_Point3D[row];

    rg_INDEX validCol = 0;
    rg_INT   first    = vOrder-1;
    rg_INT   last     = col;
    rg_INT   middle;

    while ( last >= first )
    {
        middle = (first + last)/2;

        if ( rg_LT(v, v_knot_vector[middle]) )
            last = middle;
        else if ( rg_GT(v, v_knot_vector[middle + 1]) )
            first = middle;
        else
        {
            validCol = middle;
            break;
        }
    }

    for (rg_INDEX i=0; i<row; i++)
    {
        ctrlPtOfCurve[i] = 0.;

        rg_REAL* nonZeroBasisV = evaluateMultiBasisFuncV( validCol, v, vOrder );

        for (rg_INDEX j=validCol-vOrder+1; j<=validCol; j++)
        {
            ctrlPtOfCurve[i] += rg_BSplineSurface3D::getPointOnControlNet(i,j) 
                                * nonZeroBasisV[j - validCol + vOrder - 1];
        }

        delete [] nonZeroBasisV;
    }

/*
    for (rg_INDEX i=0; i<row; i++)
    {
        ctrlPtOfCurve[i] = 0.0;
        for (rg_INDEX j=0; j<col; j++)
        {
            ctrlPtOfCurve[i] += evaluateBasisFuncV(j, v, vOrder)
                                * getPointOnControlNet(i, j);
        }
    }
*/
    rg_NUBSplineCurve3D isoparametricCurve(uOrder, row, ctrlPtOfCurve, u_knot_vector);

    delete [] ctrlPtOfCurve;
    return isoparametricCurve;
}

rg_Point3D rg_NUBSplineSurface3D::getUnitNormalVector(const rg_PARAMETER &u, const rg_PARAMETER &v)
{
    rg_NUBSplineSurface3D S_u = derivativeSurfaceOfU();
    rg_NUBSplineSurface3D S_v = derivativeSurfaceOfV();

    rg_Point3D S_u_pt = S_u.evaluatePt(u, v);
    rg_Point3D S_v_pt = S_v.evaluatePt(u, v);

    rg_Point3D unitNormal = S_u_pt*S_v_pt;

    return unitNormal/unitNormal.magnitude();
}

rg_REAL rg_NUBSplineSurface3D::getGaussianCurvature(const rg_PARAMETER &u, 
                                              const rg_PARAMETER &v)
{
    rg_NUBSplineSurface3D S_u = derivativeSurfaceOfU();
    rg_NUBSplineSurface3D S_v = derivativeSurfaceOfV();

    rg_NUBSplineSurface3D S_uu = S_u.derivativeSurfaceOfU();
    rg_NUBSplineSurface3D S_uv = S_u.derivativeSurfaceOfV();
    rg_NUBSplineSurface3D S_vv = S_v.derivativeSurfaceOfV();

    rg_Point3D S_u_pt  = S_u.evaluatePt(u, v);
    rg_Point3D S_v_pt  = S_v.evaluatePt(u, v);
    rg_Point3D S_uu_pt = S_uu.evaluatePt(u, v);
    rg_Point3D S_uv_pt = S_uv.evaluatePt(u, v);
    rg_Point3D S_vv_pt = S_vv.evaluatePt(u, v);

    rg_Point3D unitNormal = S_u_pt*S_v_pt;
    unitNormal = unitNormal/unitNormal.magnitude();

    rg_REAL GaussianCurvature =   ( (unitNormal%S_uu_pt)*(unitNormal%S_vv_pt) 
                                 - (unitNormal%S_uv_pt)*(unitNormal%S_uv_pt) )
                             / ( (S_u_pt%S_u_pt)*(S_v_pt%S_v_pt)
                                 - (S_u_pt%S_v_pt)*(S_u_pt%S_v_pt) );

    return GaussianCurvature;
}

rg_REAL** rg_NUBSplineSurface3D::getDistributionPolynomialInOneGraphPrimitive(rg_INT contributionKnotSpanInU, rg_INT contributionKnotSpanInV, rg_INT** graphInU, rg_INT** graphInV, rg_INT noOfAllPossiblePathsInU, rg_INT noOfAllPossiblePathsInV) const
{
	rg_INT p = getOrderOfU() - 1;
	rg_INT q = getOrderOfV() - 1;

	rg_REAL** distributionPolyCoeff = new rg_REAL* [p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		distributionPolyCoeff[ i ] = new rg_REAL [q + 1];
	}

	// Expanding sum of multiplications of linear terms according to graph primitive in U-dir

	rg_REAL* distributionPolyCoeffInU = new rg_REAL[p + 1];
	rg_REAL* tempPolyCoeffInU = new rg_REAL[p + 1];

	for(i = 0;i < p + 1;i++)
		distributionPolyCoeffInU[ i ] = 0.0;

	for(i = 0;i < noOfAllPossiblePathsInU;i++)
	{
		rg_INT degree = 0;
		rg_INT index = contributionKnotSpanInU;

		degree++;

		if(graphInU[ i ][ 0 ] == 0)
		{
			tempPolyCoeffInU[ 0 ] = - u_knot_vector[index] / (u_knot_vector[index + degree] - u_knot_vector[index]);
			tempPolyCoeffInU[ 1 ] = 1.0 / (u_knot_vector[index + degree] - u_knot_vector[index]);
		}
		else
		{
			index--;
			tempPolyCoeffInU[ 0 ] = u_knot_vector[index + degree + 1] / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]);
			tempPolyCoeffInU[ 1 ] = - 1.0 / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]);
		}
		//tempPolyCoeffInU[ 0 ] = coeffOfLinearTerms[ i ][ 0 ][ 0 ];
		//tempPolyCoeffInU[ 1 ] = coeffOfLinearTerms[ i ][ 0 ][ 1 ];
		for(rg_INT j = 1;j <= p - 1;j++)
		{
			degree++;

			if(graphInU[ i ][ j ] == 0)
			{
				tempPolyCoeffInU[j + 1] = tempPolyCoeffInU[ j ] * (1.0 / (u_knot_vector[index + degree] - u_knot_vector[index]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempPolyCoeffInU[ k ] = tempPolyCoeffInU[k - 1] * (1.0 / (u_knot_vector[index + degree] - u_knot_vector[index])) + tempPolyCoeffInU[ k ] * (- u_knot_vector[index] / (u_knot_vector[index + degree] - u_knot_vector[index]));
				}
				tempPolyCoeffInU[ 0 ] *= (- u_knot_vector[index] / (u_knot_vector[index + degree] - u_knot_vector[index]));
			}
			else
			{
				index--;
				tempPolyCoeffInU[j + 1] = tempPolyCoeffInU[ j ] * (- 1.0 / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempPolyCoeffInU[ k ] = tempPolyCoeffInU[k - 1] * (- 1.0 / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1])) + tempPolyCoeffInU[ k ] * (u_knot_vector[index + degree + 1] / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]));
				}
				tempPolyCoeffInU[ 0 ] *= (u_knot_vector[index + degree + 1] / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]));
			}
		}
		for(rg_INT l = 0;l < p + 1;l++)
			distributionPolyCoeffInU[ l ] += tempPolyCoeffInU[ l ]; 
	}

	delete[] tempPolyCoeffInU;

	//---------------------------------------------------------------------------------------


	
	
	// Expanding sum of multiplications of linear terms according to graph primitive in V-dir

	rg_REAL* distributionPolyCoeffInV = new rg_REAL[q + 1];
	rg_REAL* tempPolyCoeffInV = new rg_REAL[q + 1];

	for(i = 0;i < q + 1;i++)
		distributionPolyCoeffInV[ i ] = 0.0;

	for(i = 0;i < noOfAllPossiblePathsInV;i++)
	{
		rg_INT degree = 0;
		rg_INT index = contributionKnotSpanInV;

		degree++;

		if(graphInV[ i ][ 0 ] == 0)
		{
			tempPolyCoeffInV[ 0 ] = - v_knot_vector[index] / (v_knot_vector[index + degree] - v_knot_vector[index]);
			tempPolyCoeffInV[ 1 ] = 1.0 / (v_knot_vector[index + degree] - v_knot_vector[index]);
		}
		else
		{
			index--;
			tempPolyCoeffInV[ 0 ] = v_knot_vector[index + degree + 1] / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]);
			tempPolyCoeffInV[ 1 ] = - 1.0 / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]);
		}
		//tempPolyCoeffInV[ 0 ] = coeffOfLinearTerms[ i ][ 0 ][ 0 ];
		//tempPolyCoeffInV[ 1 ] = coeffOfLinearTerms[ i ][ 0 ][ 1 ];
		for(rg_INT j = 1;j <= q - 1;j++)
		{
			degree++;

			if(graphInV[ i ][ j ] == 0)
			{
				tempPolyCoeffInV[j + 1] = tempPolyCoeffInV[ j ] * (1.0 / (v_knot_vector[index + degree] - v_knot_vector[index]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempPolyCoeffInV[ k ] = tempPolyCoeffInV[k - 1] * (1.0 / (v_knot_vector[index + degree] - v_knot_vector[index])) + tempPolyCoeffInV[ k ] * (- v_knot_vector[index] / (v_knot_vector[index + degree] - v_knot_vector[index]));
				}
				tempPolyCoeffInV[ 0 ] *= (- v_knot_vector[index] / (v_knot_vector[index + degree] - v_knot_vector[index]));
			}
			else
			{
				index--;
				tempPolyCoeffInV[j + 1] = tempPolyCoeffInV[ j ] * (- 1.0 / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempPolyCoeffInV[ k ] = tempPolyCoeffInV[k - 1] * (- 1.0 / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1])) + tempPolyCoeffInV[ k ] * (v_knot_vector[index + degree + 1] / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]));
				}
				tempPolyCoeffInV[ 0 ] *= (v_knot_vector[index + degree + 1] / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]));
			}
		}
		for(rg_INT l = 0;l < q + 1;l++)
			distributionPolyCoeffInV[ l ] += tempPolyCoeffInV[ l ]; 
	}

	delete[] tempPolyCoeffInV;

	//---------------------------------------------------------------------------------------


	// Computing the tensor product of each disribution polynomial

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < q + 1;j++)
		{
			distributionPolyCoeff[ i ][ j ] = distributionPolyCoeffInU[ i ] * distributionPolyCoeffInV[ j ];
		}
	}

	delete[] distributionPolyCoeffInU;
	delete[] distributionPolyCoeffInV;

	//---------------------------------------------------------------------------------------

	return distributionPolyCoeff;
}

void rg_NUBSplineSurface3D::getPiecewiseSurfaceInPowerForm(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ) const
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

	//rg_REAL** tempPolyCoeff;

	for(i = 0;i < NoOfNonZeroLengthKnotSpanInU;i++)
	{
		polyCoeffOfX[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];
		polyCoeffOfY[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];
		polyCoeffOfZ[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];

		for(rg_INT j = 0;j < NoOfNonZeroLengthKnotSpanInV;j++)
		{
			polyCoeffOfX[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfY[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfZ[ i ][ j ] = new rg_REAL* [p + 1];

			for(rg_INT k = 0;k < p + 1;k++)
			{
				polyCoeffOfX[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfY[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfZ[ i ][ j ][ k ] = new rg_REAL [q + 1];

				for(rg_INT l = 0;l < q + 1;l++)
				{
					polyCoeffOfX[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfY[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfZ[ i ][ j ][ k ][ l ] = 0.0;
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

	for(i = p;i <= r - p - 1 && indexOfPolySurfaceInU < NoOfNonZeroLengthKnotSpanInU;i++)
	{
		indexOfPolySurfaceInV = 0;
		for(rg_INT j = q;j <= s - q - 1 && indexOfPolySurfaceInV < NoOfNonZeroLengthKnotSpanInV;j++)
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
							polyCoeffOfX[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (ctrlNet[ k ][ l ].getX() * tempPolyCoeff[indexOFBasisFuncInU][indexOFBasisFuncInV]);
							polyCoeffOfY[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (ctrlNet[ k ][ l ].getY() * tempPolyCoeff[indexOFBasisFuncInU][indexOFBasisFuncInV]);
							polyCoeffOfZ[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (ctrlNet[ k ][ l ].getZ() * tempPolyCoeff[indexOFBasisFuncInU][indexOFBasisFuncInV]);
						}
					}
					//------------------------------------------------------------------------------------------

					for(rg_INT indexU = 0;indexU < p + 1;indexU++)
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

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraphInU[ i ];j++)
		{
			delete[] allPossiblePathInU[ i ][ j ];
		}
		delete[] allPossiblePathInU[ i ];
	}
	delete[] allPossiblePathInU;
	for(i = 0;i < q + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraphInV[ i ];j++)
		{
			delete[] allPossiblePathInV[ i ][ j ];
		}
		delete[] allPossiblePathInV[ i ];
	}
	delete[] allPossiblePathInV;
	delete[] noOfAllPossiblePathsInEachGraphInU;
	delete[] noOfAllPossiblePathsInEachGraphInV;
}

void rg_NUBSplineSurface3D::makePiecewisePowerFormPolyUsingDE(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_REAL****** & truncatedBasisPolys) const
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
	truncatedBasisPolys = new rg_REAL***** [NoOfNonZeroLengthKnotSpanInU];

	//rg_REAL** tempPolyCoeff;

	for(i = 0;i < NoOfNonZeroLengthKnotSpanInU;i++)
	{
		polyCoeffOfX[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];
		polyCoeffOfY[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];
		polyCoeffOfZ[ i ] = new rg_REAL** [NoOfNonZeroLengthKnotSpanInV];
		truncatedBasisPolys[ i ] = new rg_REAL**** [NoOfNonZeroLengthKnotSpanInV];

		for(rg_INT j = 0;j < NoOfNonZeroLengthKnotSpanInV;j++)
		{
			polyCoeffOfX[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfY[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfZ[ i ][ j ] = new rg_REAL* [p + 1];
			truncatedBasisPolys[ i ][ j ] = new rg_REAL*** [p + 1];

			for(rg_INT k = 0;k < p + 1;k++)
			{
				polyCoeffOfX[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfY[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfZ[ i ][ j ][ k ] = new rg_REAL [q + 1];
				truncatedBasisPolys[ i ][ j ][ k ] = new rg_REAL** [q + 1];

				for(rg_INT l = 0;l < q + 1;l++)
				{
					polyCoeffOfX[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfY[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfZ[ i ][ j ][ k ][ l ] = 0.0;
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

	for(i = p;i <= r - p - 1 && indexOfPolySurfaceInU < NoOfNonZeroLengthKnotSpanInU;i++)
	{
		indexOfPolySurfaceInV = 0;
		for(rg_INT j = q;j <= s - q - 1 && indexOfPolySurfaceInV < NoOfNonZeroLengthKnotSpanInV;j++)
		{
			for(rg_INT k = i - p;k <= i;k++)
			{
				for(rg_INT l = j - q;l <= j;l++)
				{
					//rg_REAL** tempPolyCoeff = getDistributionPolynomialInOneGraphPrimitive(i ,j ,allPossiblePathInU[k - i + p] ,allPossiblePathInV[l - j + q] ,noOfAllPossiblePathsInEachGraphInU[k - i + p] ,noOfAllPossiblePathsInEachGraphInV[l - j + q]);
					makeTruncatedBasisPolynomial(truncatedBasisPolys[i - p][j - q][k - (i - p)][l - (j - q)], i ,j ,allPossiblePathInU[k - i + p] ,allPossiblePathInV[l - j + q] ,noOfAllPossiblePathsInEachGraphInU[k - i + p] ,noOfAllPossiblePathsInEachGraphInV[l - j + q]);

					// nested loop for computing multiplication tensor product basis function with control points
					for(rg_INT indexOFBasisFuncInU = 0;indexOFBasisFuncInU < p + 1;indexOFBasisFuncInU++)
					{
						for(rg_INT indexOFBasisFuncInV = 0;indexOFBasisFuncInV < q + 1;indexOFBasisFuncInV++)
						{
							polyCoeffOfX[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (ctrlNet[ k ][ l ].getX() * truncatedBasisPolys[i - p][j - q][k - (i - p)][l - (j - q)][indexOFBasisFuncInU][indexOFBasisFuncInV]);
							polyCoeffOfY[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (ctrlNet[ k ][ l ].getY() * truncatedBasisPolys[i - p][j - q][k - (i - p)][l - (j - q)][indexOFBasisFuncInU][indexOFBasisFuncInV]);
							polyCoeffOfZ[indexOfPolySurfaceInU][indexOfPolySurfaceInV][indexOFBasisFuncInU][indexOFBasisFuncInV] += (ctrlNet[ k ][ l ].getZ() * truncatedBasisPolys[i - p][j - q][k - (i - p)][l - (j - q)][indexOFBasisFuncInU][indexOFBasisFuncInV]);
						}
					}
					//------------------------------------------------------------------------------------------
				}
			}
			indexOfPolySurfaceInV++;
		}
		indexOfPolySurfaceInU++;
	}

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraphInU[ i ];j++)
		{
			delete[] allPossiblePathInU[ i ][ j ];
		}
		delete[] allPossiblePathInU[ i ];
	}
	delete[] allPossiblePathInU;
	for(i = 0;i < q + 1;i++)
	{
		for(rg_INT j = 0;j < noOfAllPossiblePathsInEachGraphInV[ i ];j++)
		{
			delete[] allPossiblePathInV[ i ][ j ];
		}
		delete[] allPossiblePathInV[ i ];
	}
	delete[] allPossiblePathInV;
	delete[] noOfAllPossiblePathsInEachGraphInU;
	delete[] noOfAllPossiblePathsInEachGraphInV;
}

void rg_NUBSplineSurface3D::makeTruncatedBasisPolynomial(rg_REAL** & truncatedBasisPolynomial, rg_INT contributionKnotSpanInU, rg_INT contributionKnotSpanInV, rg_INT** graphInU, rg_INT** graphInV, rg_INT noOfAllPossiblePathsInU, rg_INT noOfAllPossiblePathsInV) const
{
	rg_INT p = getOrderOfU() - 1;
	rg_INT q = getOrderOfV() - 1;

	truncatedBasisPolynomial = new rg_REAL* [p + 1];

	rg_INT i = 0;
	for(i = 0;i < p + 1;i++)
	{
		truncatedBasisPolynomial[ i ] = new rg_REAL [q + 1];
	}

	// Expanding sum of multiplications of linear terms according to graph primitive in U-dir

	rg_REAL* distributionPolyCoeffInU = new rg_REAL[p + 1];
	rg_REAL* tempPolyCoeffInU = new rg_REAL[p + 1];

	for(i = 0;i < p + 1;i++)
		distributionPolyCoeffInU[ i ] = 0.0;

	for(i = 0;i < noOfAllPossiblePathsInU;i++)
	{
		rg_INT degree = 0;
		rg_INT index = contributionKnotSpanInU;

		degree++;

		if(graphInU[ i ][ 0 ] == 0)
		{
			tempPolyCoeffInU[ 0 ] = - u_knot_vector[index] / (u_knot_vector[index + degree] - u_knot_vector[index]);
			tempPolyCoeffInU[ 1 ] = 1.0 / (u_knot_vector[index + degree] - u_knot_vector[index]);
		}
		else
		{
			index--;
			tempPolyCoeffInU[ 0 ] = u_knot_vector[index + degree + 1] / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]);
			tempPolyCoeffInU[ 1 ] = - 1.0 / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]);
		}
		//tempPolyCoeffInU[ 0 ] = coeffOfLinearTerms[ i ][ 0 ][ 0 ];
		//tempPolyCoeffInU[ 1 ] = coeffOfLinearTerms[ i ][ 0 ][ 1 ];
		for(rg_INT j = 1;j <= p - 1;j++)
		{
			degree++;

			if(graphInU[ i ][ j ] == 0)
			{
				tempPolyCoeffInU[j + 1] = tempPolyCoeffInU[ j ] * (1.0 / (u_knot_vector[index + degree] - u_knot_vector[index]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempPolyCoeffInU[ k ] = tempPolyCoeffInU[k - 1] * (1.0 / (u_knot_vector[index + degree] - u_knot_vector[index])) + tempPolyCoeffInU[ k ] * (- u_knot_vector[index] / (u_knot_vector[index + degree] - u_knot_vector[index]));
				}
				tempPolyCoeffInU[ 0 ] *= (- u_knot_vector[index] / (u_knot_vector[index + degree] - u_knot_vector[index]));
			}
			else
			{
				index--;
				tempPolyCoeffInU[j + 1] = tempPolyCoeffInU[ j ] * (- 1.0 / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempPolyCoeffInU[ k ] = tempPolyCoeffInU[k - 1] * (- 1.0 / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1])) + tempPolyCoeffInU[ k ] * (u_knot_vector[index + degree + 1] / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]));
				}
				tempPolyCoeffInU[ 0 ] *= (u_knot_vector[index + degree + 1] / (u_knot_vector[index + degree + 1] - u_knot_vector[index + 1]));
			}
		}
		for(rg_INT l = 0;l < p + 1;l++)
			distributionPolyCoeffInU[ l ] += tempPolyCoeffInU[ l ]; 
	}

	delete[] tempPolyCoeffInU;

	//---------------------------------------------------------------------------------------


	
	
	// Expanding sum of multiplications of linear terms according to graph primitive in V-dir

	rg_REAL* distributionPolyCoeffInV = new rg_REAL[q + 1];
	rg_REAL* tempPolyCoeffInV = new rg_REAL[q + 1];

	for(i = 0;i < q + 1;i++)
		distributionPolyCoeffInV[ i ] = 0.0;

	for(i = 0;i < noOfAllPossiblePathsInV;i++)
	{
		rg_INT degree = 0;
		rg_INT index = contributionKnotSpanInV;

		degree++;

		if(graphInV[ i ][ 0 ] == 0)
		{
			tempPolyCoeffInV[ 0 ] = - v_knot_vector[index] / (v_knot_vector[index + degree] - v_knot_vector[index]);
			tempPolyCoeffInV[ 1 ] = 1.0 / (v_knot_vector[index + degree] - v_knot_vector[index]);
		}
		else
		{
			index--;
			tempPolyCoeffInV[ 0 ] = v_knot_vector[index + degree + 1] / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]);
			tempPolyCoeffInV[ 1 ] = - 1.0 / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]);
		}
		//tempPolyCoeffInV[ 0 ] = coeffOfLinearTerms[ i ][ 0 ][ 0 ];
		//tempPolyCoeffInV[ 1 ] = coeffOfLinearTerms[ i ][ 0 ][ 1 ];
		for(rg_INT j = 1;j <= q - 1;j++)
		{
			degree++;

			if(graphInV[ i ][ j ] == 0)
			{
				tempPolyCoeffInV[j + 1] = tempPolyCoeffInV[ j ] * (1.0 / (v_knot_vector[index + degree] - v_knot_vector[index]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempPolyCoeffInV[ k ] = tempPolyCoeffInV[k - 1] * (1.0 / (v_knot_vector[index + degree] - v_knot_vector[index])) + tempPolyCoeffInV[ k ] * (- v_knot_vector[index] / (v_knot_vector[index + degree] - v_knot_vector[index]));
				}
				tempPolyCoeffInV[ 0 ] *= (- v_knot_vector[index] / (v_knot_vector[index + degree] - v_knot_vector[index]));
			}
			else
			{
				index--;
				tempPolyCoeffInV[j + 1] = tempPolyCoeffInV[ j ] * (- 1.0 / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]));
				for(rg_INT k = j;k >= 1;k--)
				{
					tempPolyCoeffInV[ k ] = tempPolyCoeffInV[k - 1] * (- 1.0 / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1])) + tempPolyCoeffInV[ k ] * (v_knot_vector[index + degree + 1] / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]));
				}
				tempPolyCoeffInV[ 0 ] *= (v_knot_vector[index + degree + 1] / (v_knot_vector[index + degree + 1] - v_knot_vector[index + 1]));
			}
		}
		for(rg_INT l = 0;l < q + 1;l++)
			distributionPolyCoeffInV[ l ] += tempPolyCoeffInV[ l ]; 
	}

	delete[] tempPolyCoeffInV;

	//---------------------------------------------------------------------------------------


	// Computing the tensor product of each disribution polynomial

	for(i = 0;i < p + 1;i++)
	{
		for(rg_INT j = 0;j < q + 1;j++)
		{
			truncatedBasisPolynomial[ i ][ j ] = distributionPolyCoeffInU[ i ] * distributionPolyCoeffInV[ j ];
		}
	}

	delete[] distributionPolyCoeffInU;
	delete[] distributionPolyCoeffInV;

	//---------------------------------------------------------------------------------------
}

void rg_NUBSplineSurface3D::makePiecewiseSurfaceInPowerFormUsingKR(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ) const
{
	rg_DEGREE p = getOrderOfU() - 1;
	rg_DEGREE q = getOrderOfV() - 1;
	rg_DEGREE n = getRowOfControlNet() - 1;
	rg_DEGREE m = getColumnOfControlNet() - 1;
	rg_DEGREE r = n + p + 1;
	rg_DEGREE s = m + q + 1;

	rg_INT numOfKnotSpanInU = getNumOfNonZeroLengthKnotSpanOfKnotVectorU(); // the number of nonzero length knot span in direction U
    rg_INT numOfInteriorKnotInU = numOfKnotSpanInU-1;

	rg_INT numOfKnotSpanInV = getNumOfNonZeroLengthKnotSpanOfKnotVectorV(); // the number of nonzero length knot span in direction V
    rg_INT numOfInteriorKnotInV = numOfKnotSpanInV-1;

	rg_BzSurface3D** bzSurfacePatches = decomposeSurfaceIntoBezierPatchesUsingKnotRefinement();
	//rg_BzSurface3D** bzSurfacePatches = decomposeSurfaceIntoBezierPatches();

	//polyCoeffOfX;
	//polyCoeffOfY;
	//polyCoeffOfZ;

	polyCoeffOfX = new rg_REAL*** [numOfKnotSpanInU];
	polyCoeffOfY = new rg_REAL*** [numOfKnotSpanInU];
	polyCoeffOfZ = new rg_REAL*** [numOfKnotSpanInU];

	rg_INT i = 0;
	for(i = 0;i < numOfKnotSpanInU;i++)
	{
		polyCoeffOfX[ i ] = new rg_REAL** [numOfKnotSpanInV];
		polyCoeffOfY[ i ] = new rg_REAL** [numOfKnotSpanInV];
		polyCoeffOfZ[ i ] = new rg_REAL** [numOfKnotSpanInV];

		for(rg_INT j = 0;j < numOfKnotSpanInV;j++)
		{
			polyCoeffOfX[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfY[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfZ[ i ][ j ] = new rg_REAL* [p + 1];

			for(rg_INT k = 0;k < p + 1;k++)
			{
				polyCoeffOfX[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfY[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfZ[ i ][ j ][ k ] = new rg_REAL [q + 1];

				for(rg_INT l = 0;l < q + 1;l++)
				{
					polyCoeffOfX[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfY[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfZ[ i ][ j ][ k ][ l ] = 0.0;
				}

			}
		}
	}	

	rg_REAL* distinctKnotValuesU = getDistinctKnotValuesU();
	rg_REAL* distinctKnotValuesV = getDistinctKnotValuesV();
	rg_Matrix P_x(p + 1, q + 1);
	rg_Matrix P_y(p + 1, q + 1);
	rg_Matrix P_z(p + 1, q + 1);
	rg_Matrix M_p = rg_CurveSurfaceFunc::bezierToPowerMatrix(p + 1);
	rg_Matrix M_q = rg_CurveSurfaceFunc::bezierToPowerMatrix(q + 1);

	for(i = 0;i < numOfKnotSpanInU;i++)
	{
		for(rg_INT j = 0;j < numOfKnotSpanInV;j++)
		{

			rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
										   1/(distinctKnotValuesU[i + 1] - distinctKnotValuesU[ i ]),
										   - distinctKnotValuesU[ i ] /(distinctKnotValuesU[i + 1] - distinctKnotValuesU[ i ]));

			rg_Matrix R_q = rg_CurveSurfaceFunc::reparameterMatrix(q + 1, 
										   1/(distinctKnotValuesV[j + 1] - distinctKnotValuesV[ j ]),
										   - distinctKnotValuesV[ j ] /(distinctKnotValuesV[j + 1] - distinctKnotValuesV[ j ]));
			rg_INT k = 0;
			for(k = 0;k < p + 1;k++)
			{
				for(rg_INT l = 0;l < q + 1;l++)
				{
					P_x[ k ][ l ] = bzSurfacePatches[ i ][ j ].getCtrlPt(k, l).getX();
					P_y[ k ][ l ] = bzSurfacePatches[ i ][ j ].getCtrlPt(k, l).getY();
					P_z[ k ][ l ] = bzSurfacePatches[ i ][ j ].getCtrlPt(k, l).getZ();
				}
			}

			rg_Matrix coefficientOfX = R_p * M_p * P_x * M_q.transpose() * R_q.transpose();
			rg_Matrix coefficientOfY = R_p * M_p * P_y * M_q.transpose() * R_q.transpose();
			rg_Matrix coefficientOfZ = R_p * M_p * P_z * M_q.transpose() * R_q.transpose();


			for(k = 0;k < p + 1;k++)
			{
				for(rg_INT l = 0;l < q + 1;l++)
				{
					polyCoeffOfX[ i ][ j ][ k ][ l ] = coefficientOfX[ k ][ l ];
					polyCoeffOfY[ i ][ j ][ k ][ l ] = coefficientOfY[ k ][ l ];
					polyCoeffOfZ[ i ][ j ][ k ][ l ] = coefficientOfZ[ k ][ l ];
				}
			}
		}
	}

	delete[] distinctKnotValuesU;
	delete[] distinctKnotValuesV;

	for(i = 0;i < numOfKnotSpanInU;i++)
		delete[] bzSurfacePatches[ i ];
	delete[] bzSurfacePatches;
}

#include <time.h>

void rg_NUBSplineSurface3D::makePiecewiseSurfaceInPowerFormUsingTE(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ)
{
	rg_DEGREE p = getOrderOfU() - 1;
	rg_DEGREE q = getOrderOfV() - 1;
	rg_DEGREE n = getRowOfControlNet() - 1;
	rg_DEGREE m = getColumnOfControlNet() - 1;
	rg_DEGREE r = n + p + 1;
	rg_DEGREE s = m + q + 1;

	rg_INT numOfKnotSpanInU = r - 2 * p;//getNumOfNonZeroLengthKnotSpanOfKnotVectorU(); // the number of nonzero length knot span in direction U
    rg_INT numOfInteriorKnotInU = numOfKnotSpanInU-1;

	rg_INT numOfKnotSpanInV = s - 2 * q;//getNumOfNonZeroLengthKnotSpanOfKnotVectorV(); // the number of nonzero length knot span in direction V
    rg_INT numOfInteriorKnotInV = numOfKnotSpanInV-1;

	// Calculate partial derivatives of surface to degree p and degree q

	rg_NUBSplineSurface3D** derivatives = new rg_NUBSplineSurface3D* [p + 1];

	derivatives[ 0 ] = new rg_NUBSplineSurface3D[q + 1];
	derivatives[ 0 ][ 0 ] = (*this);
	rg_INT i = 0;
	for(i = 1;i <= q;i++)
		derivatives[ 0 ][ i ] = derivatives[ 0 ][i - 1].derivativeSurfaceOfV();
	//for(i = 1;i < p;i++)
	//	derivatives[ i ][ 0 ] = derivatives[i - 1][ 0 ].derivativeSurfaceOfU();


	for(i = 1;i <= p;i++)
	{
		derivatives[ i ] = new rg_NUBSplineSurface3D[q + 1];
		derivatives[ i ][ 0 ] = derivatives[i - 1][ 0 ].derivativeSurfaceOfU();
		for(rg_INT j = 1;j <= q;j++)
			derivatives[ i ][ j ] = derivatives[i - 1][j - 1].derivativeSurfaceOfUV();
	}
	// ------------------------------------------------------------------

	polyCoeffOfX = new rg_REAL*** [numOfKnotSpanInU];
	polyCoeffOfY = new rg_REAL*** [numOfKnotSpanInU];
	polyCoeffOfZ = new rg_REAL*** [numOfKnotSpanInU];

	for(i = 0;i < numOfKnotSpanInU;i++)
	{
		polyCoeffOfX[ i ] = new rg_REAL** [numOfKnotSpanInV];
		polyCoeffOfY[ i ] = new rg_REAL** [numOfKnotSpanInV];
		polyCoeffOfZ[ i ] = new rg_REAL** [numOfKnotSpanInV];

		for(rg_INT j = 0;j < numOfKnotSpanInV;j++)
		{
			polyCoeffOfX[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfY[ i ][ j ] = new rg_REAL* [p + 1];
			polyCoeffOfZ[ i ][ j ] = new rg_REAL* [p + 1];

			for(rg_INT k = 0;k < p + 1;k++)
			{
				polyCoeffOfX[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfY[ i ][ j ][ k ] = new rg_REAL [q + 1];
				polyCoeffOfZ[ i ][ j ][ k ] = new rg_REAL [q + 1];

				for(rg_INT l = 0;l < q + 1;l++)
				{
					polyCoeffOfX[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfY[ i ][ j ][ k ][ l ] = 0.0;
					polyCoeffOfZ[ i ][ j ][ k ][ l ] = 0.0;
				}

			}
		}
	}

	rg_REAL*** TaylorExpansion = new rg_REAL** [ 3 ];
	for(i = 0;i < 3;i++)
	{
		TaylorExpansion[ i ] = new rg_REAL* [p + 1];
		rg_INT j = 0;
		for(j = 0;j <= p;j++)
			TaylorExpansion[ i ][ j ] = new rg_REAL[q + 1];
		for(j = 0;j <= p;j++)
			for(rg_INT k = 0;k <= q;k++)
				TaylorExpansion[ i ][ j ][ k ] = 0.0;
	}

	for(i = p;i <= r - p - 1;i++)
	{
		for(rg_INT j = q;j <= s - q - 1;j++)
		{
			/*
			rg_Point3D evaluatedPt = evaluatePt(u_knot_vector[ i ], v_knot_vector[ j ]);
			TaylorExpansion[ 0 ][ 0 ][ 0 ] = evaluatedPt.getX();
			TaylorExpansion[ 1 ][ 0 ][ 0 ] = evaluatedPt.getY();
			TaylorExpansion[ 2 ][ 0 ][ 0 ] = evaluatedPt.getZ();
			*/

			rg_INT k = 0;
			for(k = 0;k <= p + q;k++)
			{
				for(rg_INT l = 0;l <= k;l++)
				{
					if( ( (0 <= k - l) && (k - l <= p) ) && ( (0 <= l) && (l <= q) ) )
					{
						rg_REAL middleKnotU = (u_knot_vector[ i ] + u_knot_vector[i + 1]) / 2;
						rg_REAL middleKnotV = (v_knot_vector[ j ] + v_knot_vector[j + 1]) / 2;
						
						//rg_Point3D derivative = derivatives[k - l][ l ].evaluatePt(u_knot_vector[ i ], v_knot_vector[ j ]);
						rg_Point3D derivative = derivatives[k - l][ l ].evaluatePt(middleKnotU, middleKnotV);
						
						// Simple symbolic computation is required!!-------------------------------------------------------
						rg_REAL* multiplicationOfLinearPoly1 = new rg_REAL[k - l + 1];
						rg_INT m = 0;
						for(m = 0;m <= k - l;m++) {
							//multiplicationOfLinearPoly1[ m ] = rg_MathFunc::combination(k - l, m) * pow(- u_knot_vector[ i ], k - l - m);
							multiplicationOfLinearPoly1[ m ] = rg_MathFunc::combination(k - l, m) * pow(- middleKnotU, k - l - m);
						}

						rg_REAL* multiplicationOfLinearPoly2 = new rg_REAL[l + 1];
						for(m = 0;m <= l;m++) {
							//multiplicationOfLinearPoly2[ m ] = rg_MathFunc::combination(l, m) * pow(- v_knot_vector[ j ], l - m);
							multiplicationOfLinearPoly2[ m ] = rg_MathFunc::combination(l, m) * pow(- middleKnotV, l - m);
						}

						for(m = 0;m <= k - l;m++)
						{
							for(rg_INT n = 0;n <= l;n++)
							{
								TaylorExpansion[ 0 ][ m ][ n ] += (1 / rg_MathFunc::factorial( k )) * (rg_MathFunc::combination(k, l) * multiplicationOfLinearPoly1[ m ] * multiplicationOfLinearPoly2[ n ] * derivative.getX());
								TaylorExpansion[ 1 ][ m ][ n ] += (1 / rg_MathFunc::factorial( k )) * (rg_MathFunc::combination(k, l) * multiplicationOfLinearPoly1[ m ] * multiplicationOfLinearPoly2[ n ] * derivative.getY());
								TaylorExpansion[ 2 ][ m ][ n ] += (1 / rg_MathFunc::factorial( k )) * (rg_MathFunc::combination(k, l) * multiplicationOfLinearPoly1[ m ] * multiplicationOfLinearPoly2[ n ] * derivative.getZ());
							}
						}
						// ------------------------------------------------------------------------------------------------
					}
				}
			}

			// Complete polynomai surface in power form of one knot region
			for(k = 0;k <= p;k++)
			{
				for(rg_INT l = 0;l <= q;l++)
				{
					polyCoeffOfX[i - p][j - q][ k ][ l ] = TaylorExpansion[ 0 ][ k ][ l ];
					polyCoeffOfY[i - p][j - q][ k ][ l ] = TaylorExpansion[ 1 ][ k ][ l ];
					polyCoeffOfZ[i - p][j - q][ k ][ l ] = TaylorExpansion[ 2 ][ k ][ l ];
				}
			}

			for(k = 0;k <= p;k++)
			{
				for(rg_INT l = 0;l <= q;l++)
				{
					TaylorExpansion[ 0 ][ k ][ l ] = 0.0;
					TaylorExpansion[ 1 ][ k ][ l ] = 0.0;
					TaylorExpansion[ 2 ][ k ][ l ] = 0.0;
				}
			}
		}
	}

	for(i = 0;i < 3;i++)
	{
		for(rg_INT j = 0;j <= p;j++)
			delete[] TaylorExpansion[ i ][ j ];
		delete[] TaylorExpansion[ i ];
	}
	delete[] TaylorExpansion;
}



void rg_NUBSplineSurface3D::makeUpdatedSurfaceInPowerBasisUsingDE(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_REAL****** & truncatedBasisPolys, rg_INT** indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts)
{
	rg_DEGREE p = getOrderOfU() - 1;
	rg_DEGREE q = getOrderOfV() - 1;
	rg_DEGREE n = getRowOfControlNet() - 1;
	rg_DEGREE m = getColumnOfControlNet() - 1;
	rg_DEGREE r = n + p + 1;
	rg_DEGREE s = m + q + 1;

	for(rg_INT i = 0;i < numOfChangedCtrlPts;i++)
	{
		for(rg_INT j = indexOfChangedCtrlPts[ i ][ 0 ];j < indexOfChangedCtrlPts[ i ][ 0 ] + p + 1;j++)
		{
			for(rg_INT k = indexOfChangedCtrlPts[ i ][ 1 ];k < indexOfChangedCtrlPts[ i ][ 1 ] + q + 1;k++)
			{
				if( ((p <= j) && (j <= r - p - 1)) && ((q <= k) && (k <= s - q - 1)) )
				{
					rg_Point3D ctrlPt = getPointOnControlNet(indexOfChangedCtrlPts[ i ][ 0 ], indexOfChangedCtrlPts[ i ][ 1 ]);
					for(rg_INT l = 0;l < p + 1;l++)
					{
						for(rg_INT m = 0;m < q + 1;m++)
						{
							polyCoeffOfX[j - p][k - q][ l ][ m ] += (ctrlPt.getX() * truncatedBasisPolys[j - p][k - q][indexOfChangedCtrlPts[ i ][ 0 ] - (j - p)][indexOfChangedCtrlPts[ i ][ 1 ] - (k - q)][ l ][ m ]);
							polyCoeffOfY[j - p][k - q][ l ][ m ] += (ctrlPt.getY() * truncatedBasisPolys[j - p][k - q][indexOfChangedCtrlPts[ i ][ 0 ] - (j - p)][indexOfChangedCtrlPts[ i ][ 1 ] - (k - q)][ l ][ m ]);
							polyCoeffOfZ[j - p][k - q][ l ][ m ] += (ctrlPt.getZ() * truncatedBasisPolys[j - p][k - q][indexOfChangedCtrlPts[ i ][ 0 ] - (j - p)][indexOfChangedCtrlPts[ i ][ 1 ] - (k - q)][ l ][ m ]);
						}
					}
				}
			}
		}
	}
}

void rg_NUBSplineSurface3D::makeUpdatedSurfaceInPowerBasisUsingKR(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_INT** indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts)
{
	rg_INT i = 0;
	for(i = 0;i < numOfChangedCtrlPts;i++)
	{
		rg_DEGREE p = getOrderOfU() - 1;
		rg_DEGREE q = getOrderOfV() - 1;
		rg_DEGREE n = getRowOfControlNet() - 1;
		rg_DEGREE m = getColumnOfControlNet() - 1;
		rg_DEGREE r = n + p + 1;
		rg_DEGREE s = m + q + 1;

		rg_REAL* distinctKnotValuesInU = getDistinctKnotValuesU();
		rg_REAL* distinctKnotValuesInV = getDistinctKnotValuesV();
		rg_Matrix P_x(p + 1, q + 1);
		rg_Matrix P_y(p + 1, q + 1);
		rg_Matrix P_z(p + 1, q + 1);
		rg_Matrix M_p = rg_CurveSurfaceFunc::bezierToPowerMatrix(p + 1);
		rg_Matrix M_q = rg_CurveSurfaceFunc::bezierToPowerMatrix(q + 1);

		for(rg_INT j = indexOfChangedCtrlPts[ i ][ 0 ];j < indexOfChangedCtrlPts[ i ][ 0 ] + p + 1;j++)
		{
			for(rg_INT k = indexOfChangedCtrlPts[ i ][ 1 ];k < indexOfChangedCtrlPts[ i ][ 1 ] + q + 1;k++)
			{
				if((rg_NZERO(u_knot_vector[j + 1] - u_knot_vector[ j ])) && (rg_NZERO(v_knot_vector[k + 1] - v_knot_vector[ k ])))
				{
					// knot refinement for conversion to Bernstein basis---------------------------------------------------
					rg_NUBSplineSurface3D duplicatedSurface = (*this);

					// In U-direction
					rg_INT multiplicityPreInU = 1;
					rg_INT l = j;
					while(l >= 1)
					{
						if(rg_ZERO(u_knot_vector[ j ] - u_knot_vector[l - 1]))
							multiplicityPreInU++;
						else
							break;
						l--;
					}
					rg_INT multiplicityNextInU = 1;
					l = j + 1;
					while(l <= r - p - 2)
					{
						if(rg_ZERO(u_knot_vector[j + 1] - u_knot_vector[l + 1]))
							multiplicityNextInU++;
						else
							break;
						l++;
					}
					rg_INT numOfInsertingKnotValuesPrevInU = rg_Max(0, p - multiplicityPreInU);
					rg_INT numOfInsertingKnotValuesNextInU = rg_Max(0, p - multiplicityNextInU);
					rg_REAL* insertingKnotValuesPrevInU = new rg_REAL [numOfInsertingKnotValuesPrevInU];
					rg_REAL* insertingKnotValuesNextInU = new rg_REAL [numOfInsertingKnotValuesNextInU];
					for(l = 0;l < numOfInsertingKnotValuesPrevInU;l++)
						insertingKnotValuesPrevInU[ l ] = u_knot_vector[ j ];
					for(l = 0;l < numOfInsertingKnotValuesNextInU;l++)
						insertingKnotValuesNextInU[ l ] = u_knot_vector[j + 1];					
					// In V-direction
					rg_INT multiplicityPreInV = 1;
					rg_INT m = k;
					while(m >= 1)
					{
						if(rg_ZERO(v_knot_vector[ k ] - v_knot_vector[m - 1]))
							multiplicityPreInV++;
						else
							break;
						m--;
					}
					rg_INT multiplicityNextInV = 1;
					m = k + 1;
					while(m <= s - q - 2)
					{
						if(rg_ZERO(v_knot_vector[k + 1] - v_knot_vector[m + 1]))
							multiplicityNextInV++;
						else
							break;
						m++;
					}
					rg_INT numOfInsertingKnotValuesPrevInV = rg_Max(0, q - multiplicityPreInV);
					rg_INT numOfInsertingKnotValuesNextInV = rg_Max(0, q - multiplicityNextInV);
					rg_REAL* insertingKnotValuesPrevInV = new rg_REAL [numOfInsertingKnotValuesPrevInV];
					rg_REAL* insertingKnotValuesNextInV = new rg_REAL [numOfInsertingKnotValuesNextInV];
					for(m = 0;m < numOfInsertingKnotValuesPrevInV;m++)
						insertingKnotValuesPrevInV[ m ] = v_knot_vector[ k ];
					for(m = 0;m < numOfInsertingKnotValuesNextInV;m++)
						insertingKnotValuesNextInV[ m ] = v_knot_vector[k + 1];


					//if((numOfInsertingKnotValuesPrevInU > 0) && (numOfInsertingKnotValuesPrevInV > 0))
					//	duplicatedSurface.knotRefinement(insertingKnotValuesPrevInU, numOfInsertingKnotValuesPrevInU, insertingKnotValuesPrevInV, numOfInsertingKnotValuesPrevInV);
					//if((numOfInsertingKnotValuesNextInU > 0) && (numOfInsertingKnotValuesNextInV > 0))
					//	duplicatedSurface.knotRefinement(insertingKnotValuesNextInU, numOfInsertingKnotValuesNextInU, insertingKnotValuesNextInV, numOfInsertingKnotValuesNextInV);

					duplicatedSurface.knotRefinementOfU(insertingKnotValuesPrevInU, numOfInsertingKnotValuesPrevInU);
					duplicatedSurface.knotRefinementOfU(insertingKnotValuesNextInU, numOfInsertingKnotValuesNextInU);
					duplicatedSurface.knotRefinementOfV(insertingKnotValuesPrevInV, numOfInsertingKnotValuesPrevInV);
					duplicatedSurface.knotRefinementOfV(insertingKnotValuesNextInV, numOfInsertingKnotValuesNextInV);

					//----------------------------------------------------------------------------------------------------

					// Conversion of a Bezier patch to a Power basis patch------------------------------------------------
/*
					rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
												   1/(distinctKnotValuesInU[j + 1] - distinctKnotValuesInU[ j ]),
												   - distinctKnotValuesInU[ j ] /(distinctKnotValuesInU[j + 1] - distinctKnotValuesInU[ j ]));

					rg_Matrix R_q = rg_CurveSurfaceFunc::reparameterMatrix(q + 1,
												   1/(distinctKnotValuesInV[k + 1] - distinctKnotValuesInV[ k ]),
												   - distinctKnotValuesInV[ k ] /(distinctKnotValuesInV[k + 1] - distinctKnotValuesInV[ k ]));
*/
					rg_Matrix R_p = rg_CurveSurfaceFunc::reparameterMatrix(p + 1, 
												   1/(distinctKnotValuesInU[j - p + 1] - distinctKnotValuesInU[j - p]),
												   - distinctKnotValuesInU[j - p] /(distinctKnotValuesInU[j - p + 1] - distinctKnotValuesInU[j - p]));

					rg_Matrix R_q = rg_CurveSurfaceFunc::reparameterMatrix(q + 1,
												   1/(distinctKnotValuesInV[k - q + 1] - distinctKnotValuesInV[k - q]),
												   - distinctKnotValuesInV[k - q] /(distinctKnotValuesInV[k - q + 1] - distinctKnotValuesInV[k - q]));

					// determine the index of starting Bezier points
					rg_INT indexOfStartCtrlPtOfRow, indexOfStartCtrlPtOfCol;
					
					if(j <= p + 1)
						indexOfStartCtrlPtOfRow = (j - p) * p;
					else
						indexOfStartCtrlPtOfRow = (j - p) + (p - 1); // fucking index!!

					if(k <= q + 1)
						indexOfStartCtrlPtOfCol = (k - q) * q;
					else
						indexOfStartCtrlPtOfCol = (k - q) + (q - 1); // fucking index!!

					for(l = 0;l < p + 1;l++)
					{
						for(m = 0;m < q + 1;m++)
						{
							P_x[ l ][ m ] = duplicatedSurface.getPointOnControlNet(indexOfStartCtrlPtOfRow + l, indexOfStartCtrlPtOfCol + m).getX();
							P_y[ l ][ m ] = duplicatedSurface.getPointOnControlNet(indexOfStartCtrlPtOfRow + l, indexOfStartCtrlPtOfCol + m).getY();
							P_z[ l ][ m ] = duplicatedSurface.getPointOnControlNet(indexOfStartCtrlPtOfRow + l, indexOfStartCtrlPtOfCol + m).getZ();
						}
					}

					rg_Matrix coefficientOfX = R_p * M_p * P_x * M_q.transpose() * R_q.transpose();
					rg_Matrix coefficientOfY = R_p * M_p * P_y * M_q.transpose() * R_q.transpose();
					rg_Matrix coefficientOfZ = R_p * M_p * P_z * M_q.transpose() * R_q.transpose();

					for(l = 0;l < p + 1;l++)
					{
						for(m = 0;m < q + 1;m++)
						{
							polyCoeffOfX[j - p][k - q][ l ][ m ] = coefficientOfX[ l ][ m ];
							polyCoeffOfY[j - p][k - q][ l ][ m ] = coefficientOfY[ l ][ m ];
							polyCoeffOfZ[j - p][k - q][ l ][ m ] = coefficientOfZ[ l ][ m ];
						}
					}
					//----------------------------------------------------------------------------------------------------
					delete[] insertingKnotValuesPrevInU;
					delete[] insertingKnotValuesNextInU;
					delete[] insertingKnotValuesPrevInV;
					delete[] insertingKnotValuesNextInV;
				}
			}
		}
		delete[] distinctKnotValuesInU;
		delete[] distinctKnotValuesInV;
	}

}

void rg_NUBSplineSurface3D::makeUpdatedSurfaceInPowerBasisUsingTE(rg_REAL**** & polyCoeffOfX, rg_REAL**** & polyCoeffOfY, rg_REAL**** & polyCoeffOfZ, rg_INT** indexOfChangedCtrlPts, rg_INT numOfChangedCtrlPts)
{
	rg_INT i = 0;
	for(i = 0;i < numOfChangedCtrlPts;i++)
	{
		rg_DEGREE p = getOrderOfU() - 1;
		rg_DEGREE q = getOrderOfV() - 1;
		rg_DEGREE n = getRowOfControlNet() - 1;
		rg_DEGREE m = getColumnOfControlNet() - 1;
		rg_DEGREE r = n + p + 1;
		rg_DEGREE s = m + q + 1;

		// Calculate partial derivatives of surface to degree p and degree q

		rg_NUBSplineSurface3D** derivatives = new rg_NUBSplineSurface3D* [p + 1];

		derivatives[ 0 ] = new rg_NUBSplineSurface3D[q + 1];
		derivatives[ 0 ][ 0 ] = (*this);
		rg_INT j = 0;
		for(j = 1;j <= q;j++)
			derivatives[ 0 ][ j ] = derivatives[ 0 ][j - 1].derivativeSurfaceOfV();
		//for(i = 1;i < p;i++)
		//	derivatives[ i ][ 0 ] = derivatives[i - 1][ 0 ].derivativeSurfaceOfU();


		for(j = 1;j <= p;j++)
		{
			derivatives[ j ] = new rg_NUBSplineSurface3D[q + 1];
			derivatives[ j ][ 0 ] = derivatives[j - 1][ 0 ].derivativeSurfaceOfU();
			for(rg_INT k = 1;k <= q;k++)
				derivatives[ j ][ k ] = derivatives[j - 1][k - 1].derivativeSurfaceOfUV();
		}
		// ------------------------------------------------------------------

		rg_REAL*** TaylorExpansion = new rg_REAL** [ 3 ];
		for(j = 0;j < 3;j++)
		{
			TaylorExpansion[ j ] = new rg_REAL* [p + 1];
			rg_INT k = 0;
			for(k = 0;k <= p;k++)
				TaylorExpansion[ j ][ k ] = new rg_REAL[q + 1];
			for(k = 0;k <= p;k++)
				for(rg_INT l = 0;l <= q;l++)
					TaylorExpansion[ j ][ k ][ l ] = 0.0;
		}

		for(j = indexOfChangedCtrlPts[ i ][ 0 ];j < indexOfChangedCtrlPts[ i ][ 0 ] + p + 1;j++)
		{
			for(rg_INT k = indexOfChangedCtrlPts[ i ][ 1 ];k < indexOfChangedCtrlPts[ i ][ 1 ] + q + 1;k++)
			{
				if((rg_NZERO(u_knot_vector[j + 1] - u_knot_vector[ j ])) && (rg_NZERO(v_knot_vector[k + 1] - v_knot_vector[ k ])))
				{
					/*
					rg_Point3D evaluatedPt = evaluatePt(u_knot_vector[ j ], v_knot_vector[ k ]);
					TaylorExpansion[ 0 ][ 0 ][ 0 ] = evaluatedPt.getX();
					TaylorExpansion[ 1 ][ 0 ][ 0 ] = evaluatedPt.getY();
					TaylorExpansion[ 2 ][ 0 ][ 0 ] = evaluatedPt.getZ();
					*/

					rg_INT l = 0;
					for(l = 0;l <= p + q;l++)
					{
						for(rg_INT m = 0;m <= l;m++)
						{
							if( ( (0 <= l - m) && (l - m <= p) ) && ( (0 <= m) && (m <= q) ) )
							{
								rg_REAL middleKnotU = (u_knot_vector[ j ] + u_knot_vector[j + 1]) / 2;
								rg_REAL middleKnotV = (v_knot_vector[ k ] + v_knot_vector[k + 1]) / 2;

								//rg_Point3D derivative = derivatives[l - m][ m ].evaluatePt(u_knot_vector[ j ], v_knot_vector[ k ]);
								rg_Point3D derivative = derivatives[l - m][ m ].evaluatePt(middleKnotU, middleKnotV);
								
								// Simple symbolic computation is required!!-------------------------------------------------------
								rg_REAL* multiplicationOfLinearPoly1 = new rg_REAL[l - m + 1];
								rg_INT n = 0;
								for(n = 0;n <= l - m;n++) {
									//multiplicationOfLinearPoly1[ n ] = rg_MathFunc::combination(l - m, n) * pow(- u_knot_vector[ j ], l - m - n);
									multiplicationOfLinearPoly1[ n ] = rg_MathFunc::combination(l - m, n) * pow(- middleKnotU, l - m - n);
								}

								rg_REAL* multiplicationOfLinearPoly2 = new rg_REAL[m + 1];
								for(n = 0;n <= m;n++) {
									//multiplicationOfLinearPoly2[ n ] = rg_MathFunc::combination(m, n) * pow(- v_knot_vector[ k ], m - n);
									multiplicationOfLinearPoly2[ n ] = rg_MathFunc::combination(m, n) * pow(- middleKnotV, m - n);
								}

								for(n = 0;n <= l - m;n++)
								{
									for(rg_INT o = 0;o <= m;o++)
									{
										TaylorExpansion[ 0 ][ n ][ o ] += (1 / rg_MathFunc::factorial( l )) * (rg_MathFunc::combination(l, m) * multiplicationOfLinearPoly1[ n ] * multiplicationOfLinearPoly2[ o ] * derivative.getX());
										TaylorExpansion[ 1 ][ n ][ o ] += (1 / rg_MathFunc::factorial( l )) * (rg_MathFunc::combination(l, m) * multiplicationOfLinearPoly1[ n ] * multiplicationOfLinearPoly2[ o ] * derivative.getY());
										TaylorExpansion[ 2 ][ n ][ o ] += (1 / rg_MathFunc::factorial( l )) * (rg_MathFunc::combination(l, m) * multiplicationOfLinearPoly1[ n ] * multiplicationOfLinearPoly2[ o ] * derivative.getZ());
									}
								}
								// ------------------------------------------------------------------------------------------------
							}
						}
					}

					// Complete polynomai surface in power form of one knot region
					for(l = 0;l <= p;l++)
					{
						for(rg_INT m = 0;m <= q;m++)
						{
							polyCoeffOfX[j - p][k - q][ l ][ m ] = TaylorExpansion[ 0 ][ l ][ m ];
							polyCoeffOfY[j - p][k - q][ l ][ m ] = TaylorExpansion[ 1 ][ l ][ m ];
							polyCoeffOfZ[j - p][k - q][ l ][ m ] = TaylorExpansion[ 2 ][ l ][ m ];
						}
					}

					for(l = 0;l <= p;l++)
					{
						for(rg_INT m = 0;m <= q;m++)
						{
							TaylorExpansion[ 0 ][ l ][ m ] = 0.0;
							TaylorExpansion[ 1 ][ l ][ m ] = 0.0;
							TaylorExpansion[ 2 ][ l ][ m ] = 0.0;
						}
					}
				}
			}
		}
		for(j = 0;j < 3;j++)
		{
			for(rg_INT k = 0;k <= p;k++)
				delete[] TaylorExpansion[ j ][ k ];
			delete[] TaylorExpansion[ j ];
		}
		delete[] TaylorExpansion;		
	}
}

////	Set Functions.---------------------------------------------------------
void rg_NUBSplineSurface3D::setInitialKnotVectorOfU()
{
	rg_INT row	   = rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();

	u_knot_vector = new rg_REAL[row + uOrder];

	for (rg_INT i=0; i<(row+uOrder); i++)
	{
		if ( i < uOrder )
			u_knot_vector[i] = 0.;
		else if ( i < row )
			u_knot_vector[i] = (i - uOrder + 1.) / (row - uOrder + 1.);
		else
			u_knot_vector[i] = 1.;
	}
}

void rg_NUBSplineSurface3D::setInitialKnotVectorOfV()
{
	rg_INT col    = rg_BSplineSurface3D::getColumnOfControlNet();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();

	v_knot_vector = new rg_REAL [col + vOrder];

	for (rg_INT i=0; i<(col+vOrder); i++)
	{
		if ( i < vOrder )
			v_knot_vector[i] = 0.;
		else if ( i < col )
			v_knot_vector[i] = (i - vOrder + 1.) / (col - vOrder + 1.);
		else
			v_knot_vector[i] = 1.;
	}
}

void rg_NUBSplineSurface3D::setKnotValueOfU( const rg_INDEX &kIndex, 
                                          const rg_REAL  &knotValue )
{
	u_knot_vector[kIndex] = knotValue;
}

void rg_NUBSplineSurface3D::setKnotValueOfV( const rg_INDEX &kIndex, 
                                          const rg_REAL  &knotValue )
{
	v_knot_vector[kIndex] = knotValue;
}

void rg_NUBSplineSurface3D::setKnotVectorOfU( const rg_INT &numOfKnot,
                                           const rg_REAL* const newKnotVector )
{
    rg_INT row	   = rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();

    if ( numOfKnot != (row + uOrder) )
        return;

    if ( u_knot_vector != rg_NULL )
        delete [] u_knot_vector;

    u_knot_vector = new rg_REAL[row + uOrder];
    for (rg_INT i=0; i<row+uOrder; i++)
	    u_knot_vector[i] = newKnotVector[i];
}

void rg_NUBSplineSurface3D::setKnotVectorOfV( const rg_INT &numOfKnot,
                                           const rg_REAL* const newKnotVector )
{
	rg_INT col    = rg_BSplineSurface3D::getColumnOfControlNet();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();

    if ( numOfKnot != (col + vOrder) )
        return;

    if ( v_knot_vector != rg_NULL )
        delete [] v_knot_vector;

    v_knot_vector = new rg_REAL[numOfKnot];
    for (rg_INT j=0; j<numOfKnot; j++)
	    v_knot_vector[j] = newKnotVector[j];
}

////	Operating & Calculating.-----------------------------------------------
rg_REAL rg_NUBSplineSurface3D::evaluateBasisFuncU( const rg_INDEX     &index, 
                                             const rg_PARAMETER &u,
                                             const rg_ORDER     &uOrder)
{
	rg_INT row    = rg_BSplineSurface3D::getRowOfControlNet() - 1;
//    if (uOrder == -1)
//	    uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();

    if (    (index == 0 && rg_EQ(u, u_knot_vector[0]) )
	     || (index == row && rg_EQ(u, u_knot_vector[row + uOrder])) )
    {
	    return 1.0;
    }
    else if (    rg_LT(u, u_knot_vector[index])
	          || rg_GE(u, u_knot_vector[index + uOrder]) )
//    else if ( rg_EQ(u, u_knot_vector[index + uOrder]) )
    {
	    return 0.0;
    }
    else 
    {}

    rg_REAL Uleft  = 0.;
    rg_REAL Uright = 0.;
    rg_REAL temp   = 0.;
    rg_REAL saved  = 0.;

    rg_REAL* triN = new rg_REAL [uOrder];

    rg_INT i = 0;
	for (i=0; i<uOrder; i++)
    {
	    if (    rg_LT(u, u_knot_vector[index + i + 1]) 
	         && rg_GE(u, u_knot_vector[index + i]) ) 
	    {
		      triN[i] = 1.0;
	    }
	    else 
		    triN[i] = 0.0;
    }

    for (i=1; i<uOrder; i++)
    {
	    if ( rg_EQ(triN[0], 0.0) )
		    saved = 0.0;
	    else 
	    {
		    saved = ( (u - u_knot_vector[index]) * triN[0] )
			        / ( u_knot_vector[index+i] - u_knot_vector[index] );
	    }

	    for (rg_INT j=0; j<(uOrder-i); j++)
	    {
		    Uleft  = u_knot_vector[index + j + 1];
		    Uright = u_knot_vector[index + j + i + 1];

		    if ( rg_EQ(triN[j+1], 0.0) )
		    {
			    triN[j] = saved;
			    saved   = 0.0;
		    }
		    else
		    {
			    temp    = triN[j+1] / (Uright - Uleft);
			    triN[j] = saved + (Uright - u) * temp;
			    saved   = (u - Uleft) * temp;
		    }
	    }
    }

    rg_REAL returnValue = triN[0];
    delete [] triN;

    return returnValue;
}

rg_REAL rg_NUBSplineSurface3D::evaluateBasisFuncV( const rg_INDEX     &index, 
                                             const rg_PARAMETER &v,
                                             const rg_ORDER     &vOrder)
{
	rg_INT col    = rg_BSplineSurface3D::getColumnOfControlNet() - 1;
//    if (vOrder == -1)
//	    vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();

    if (    ( index == 0 && rg_EQ(v, v_knot_vector[0]) )
	     || ( index == col && rg_EQ(v, v_knot_vector[col + vOrder]) ) )
    {
	    return 1.0;
    }
    else if (    rg_LT( v, v_knot_vector[index] )
	          || rg_GE( v, v_knot_vector[index + vOrder] ) )
    {
	    return 0.0;
    }
    else 
    {}

    rg_REAL Uleft  = 0.;
    rg_REAL Uright = 0.;
    rg_REAL temp   = 0.;
    rg_REAL saved  = 0.;

    rg_REAL* triN = new rg_REAL [vOrder];

    rg_INT i = 0;
	for (i=0; i<vOrder; i++)
    {
	    if (    rg_LT( v, v_knot_vector[index + i + 1] ) 
	         && rg_GE( v, v_knot_vector[index + i] ) ) 
	    {
		      triN[i] = 1.0;
	    }
	    else 
		    triN[i] = 0.0;
    }

    for (i=1; i<vOrder; i++)
    {
	    if ( rg_EQ(triN[0], 0.0) )
		    saved = 0.0;
	    else 
	    {
		    saved = ( (v - v_knot_vector[index]) * triN[0] )
			        / ( v_knot_vector[index+i] - v_knot_vector[index] );
	    }

	    for (rg_INT j=0; j<(vOrder-i); j++)
	    {
		    Uleft  = v_knot_vector[index + j + 1];
		    Uright = v_knot_vector[index + j + i + 1];

		    if ( rg_EQ(triN[j+1], 0.0) )
		    {
			    triN[j] = saved;
			    saved   = 0.0;
		    }
		    else
		    {
			    temp    = triN[j+1] / (Uright - Uleft);
			    triN[j] = saved + (Uright - v) * temp;
			    saved   = (v - Uleft) * temp;
		    }
	    }
    }

    rg_REAL returnValue = triN[0];
    delete [] triN;

    return returnValue;
}

rg_REAL* rg_NUBSplineSurface3D::evaluateMultiBasisFuncU( const rg_INDEX     &uKnotIndex,
                                                   const rg_PARAMETER &u,
                                                   const rg_ORDER     &uOrder)
{
    rg_REAL* nonZeroBasis = new rg_REAL[uOrder];

    rg_REAL  tempBasis  = 0.0;
    rg_INDEX basisIndex = 0; 

    nonZeroBasis[0] = 1.0;
    //  N            (u)
    //   knotIndex, 1   
    for (rg_INDEX basisOrder=2; basisOrder<=uOrder; basisOrder++)
    {
        for (rg_INDEX j=basisOrder-1; j>=0; j--)
        {
            basisIndex = uKnotIndex - (basisOrder - j) + 1;
            if ( j == (basisOrder-1) )
            {
                nonZeroBasis[j] 
//                tempBasis
                    = (u - u_knot_vector[basisIndex]) * nonZeroBasis[j-1] 
                      / (u_knot_vector[basisIndex + basisOrder - 1] - u_knot_vector[basisIndex]);
            }
            else if ( j != 0 ) // else if ( j != 0) )
            {
                nonZeroBasis[j] 
//                tempBasis
                    = (u - u_knot_vector[basisIndex]) * nonZeroBasis[j-1] 
                      / (u_knot_vector[basisIndex +  basisOrder - 1] - u_knot_vector[basisIndex])
                      + 
                      (u_knot_vector[basisIndex + basisOrder] - u) * nonZeroBasis[j]
                      / (u_knot_vector[basisIndex + basisOrder] - u_knot_vector[basisIndex + 1]);       
            }
            else // if (j == 0)
            {
                nonZeroBasis[j]
//                tempBasis
                    = (u_knot_vector[basisIndex + basisOrder] - u) * nonZeroBasis[0]
                      / (u_knot_vector[basisIndex + basisOrder] - u_knot_vector[basisIndex + 1]);
            }

//            nonZeroBasis[j] = tempBasis;
        }
    }

    return nonZeroBasis;
}

rg_REAL* rg_NUBSplineSurface3D::evaluateMultiBasisFuncV( const rg_INDEX     &vKnotIndex,
                                                   const rg_PARAMETER &v,
                                                   const rg_ORDER     &vOrder)
{
    rg_REAL* nonZeroBasis = new rg_REAL[vOrder];

    rg_REAL  tempBasis  = 0.0;
    rg_INDEX basisIndex = 0;

    nonZeroBasis[0] = 1.0;
    //  N            (u)
    //   knotIndex, 1   
    for (rg_INDEX basisOrder=2; basisOrder<=vOrder; basisOrder++)
    {
        for (rg_INDEX j=basisOrder-1; j>=0; j--)
        {
            basisIndex = vKnotIndex-(basisOrder-j)+1;
            if ( j == (basisOrder-1) )
            {
                nonZeroBasis[j] 
//                tempBasis
                    = (v - v_knot_vector[basisIndex]) * nonZeroBasis[j-1] 
                      / (v_knot_vector[basisIndex + basisOrder - 1] - v_knot_vector[basisIndex]);
            }
            else if ( j != 0 ) // else if ( j != 0) )
            {
                nonZeroBasis[j] 
//                tempBasis
                    = (v - v_knot_vector[basisIndex]) * nonZeroBasis[j-1] 
                      / (v_knot_vector[basisIndex +  basisOrder - 1] - v_knot_vector[basisIndex])
                      + 
                      (v_knot_vector[basisIndex + basisOrder] - v) * nonZeroBasis[j]
                      / (v_knot_vector[basisIndex + basisOrder] - v_knot_vector[basisIndex + 1]);       
            }
            else // if (j == 0)
            {
                nonZeroBasis[j]
//                tempBasis
                    = (v_knot_vector[basisIndex + basisOrder] - v) * nonZeroBasis[0]
                      / (v_knot_vector[basisIndex + basisOrder] - v_knot_vector[basisIndex + 1]);
            }

//            nonZeroBasis[j] = tempBasis;
        }
    }

    return nonZeroBasis;
}

//  April  3 1997 : Modified
////////////////////////////////////////////////////////////////// 
rg_Point3D rg_NUBSplineSurface3D::evaluatePt( const rg_PARAMETER &u, 
                                      const rg_PARAMETER &v )
{
    if (u_knot_vector == rg_NULL)
        setInitialKnotVectorOfU();
	
	if (v_knot_vector == rg_NULL)
        setInitialKnotVectorOfV();
  
    rg_Point3D ptOnSurface;

	rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();

    rg_REAL uBasisValue = 0.;
	rg_REAL vBasisValue = 0.;

    rg_INDEX validRow = 0;
    rg_INDEX validCol = 0;

    rg_INT first = uOrder-1;
    rg_INT last  = row;
    rg_INT middle;

    while ( last >= first )
    {
        middle = (first + last)/2;

        if ( rg_LT(u, u_knot_vector[middle]) )
            last = middle;
        else if ( rg_GT(u, u_knot_vector[middle + 1]) )
            first = middle + 1;
        else
        {
			if ( middle != first )
			{
	            while ( rg_EQ(u_knot_vector[middle], u_knot_vector[middle+1]) )
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
    last  = col;
    while ( last >= first )
    {
        middle = (first + last)/2;

        if ( rg_LT(v, v_knot_vector[middle]) )
            last = middle;
        else if ( rg_GT(v, v_knot_vector[middle + 1]) )
            first = middle + 1;
        else
        {
			if ( middle != first )
			{
	            while ( rg_EQ(v_knot_vector[middle], v_knot_vector[middle+1]) )
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
/*
    for (rg_INDEX i=validRow; i<(validRow+uOrder); i++)
    {
        uBasisValue = evaluateBasisFuncU(i, u, uOrder);

        for (rg_INDEX j=validCol; j<(validCol+vOrder); j++)
        {
            vBasisValue = evaluateBasisFuncV(j, v, vOrder);

            ptOnSurface += rg_BSplineSurface3D::getPointOnControlNet(i,j) 
                           * uBasisValue * vBasisValue;              
        }
    }

    return ptOnSurface;
*/

    rg_REAL* nonZeroBasisU = evaluateMultiBasisFuncU( validRow, u, uOrder );
    rg_REAL* nonZeroBasisV = evaluateMultiBasisFuncV( validCol, v, vOrder );

    for (rg_INDEX i=validRow-uOrder+1; i<=validRow; i++)
    {
        for (rg_INDEX j=validCol-vOrder+1; j<=validCol; j++)
        {
            ptOnSurface += rg_BSplineSurface3D::getPointOnControlNet(i,j) 
                           * nonZeroBasisU[i - validRow + uOrder - 1]
                           * nonZeroBasisV[j - validCol + vOrder - 1];              
        }
    }

    delete [] nonZeroBasisU;
    delete [] nonZeroBasisV;

    return ptOnSurface;
}

////    Derivative.------------------------------------------------------------
rg_NUBSplineSurface3D rg_NUBSplineSurface3D::derivativeSurfaceOfU() const
{
	rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();
    
    rg_NUBSplineSurface3D derivativeSurface( row-1, col, uOrder-1, vOrder );

    rg_Point3D newCtrlPt;
    rg_INT i = 0;
	for (i=0; i<row-1; i++)
    {
        for (rg_INT j=0; j<col; j++)
        {
            newCtrlPt 
                = ( uOrder-1 )
                  * (getPointOnControlNet(i+1,j) - getPointOnControlNet(i,j))
                  / (u_knot_vector[i+uOrder] - u_knot_vector[i+1]);

            derivativeSurface.setPointOnControlNet(i, j, newCtrlPt);
        }
    }     

    derivativeSurface.u_knot_vector = new rg_REAL [row + uOrder - 2];
    for (i=0; i<(row+uOrder-2); i++)
        derivativeSurface.u_knot_vector[i] = u_knot_vector[i+1];

    derivativeSurface.v_knot_vector = new rg_REAL [col + vOrder]; 
    for (rg_INT j=0; j<(col+vOrder); j++)
        derivativeSurface.v_knot_vector[j] = v_knot_vector[j];

    return derivativeSurface;
}

rg_NUBSplineSurface3D rg_NUBSplineSurface3D::derivativeSurfaceOfV() const
{
	rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();
    
    rg_NUBSplineSurface3D derivativeSurface( row, col-1, uOrder, vOrder-1 );

    rg_Point3D newCtrlPt;
    rg_INT i = 0;
	for (i=0; i<row; i++)
    {
        for (rg_INT j=0; j<col-1; j++)
        {
            newCtrlPt 
                = ( vOrder-1 )
                  * (getPointOnControlNet(i,j+1) - getPointOnControlNet(i,j))
                  / ( v_knot_vector[j+vOrder] - v_knot_vector[j+1] );

            derivativeSurface.setPointOnControlNet(i, j, newCtrlPt);
        }
    }     

    derivativeSurface.u_knot_vector = new rg_REAL [row + uOrder];
    for (i=0; i<(row+uOrder); i++)
        derivativeSurface.u_knot_vector[i] = u_knot_vector[i];

    derivativeSurface.v_knot_vector = new rg_REAL [col + vOrder - 2]; 
    for (rg_INT j=0; j<(col+vOrder-2); j++)
        derivativeSurface.v_knot_vector[j] = v_knot_vector[j+1];

    return derivativeSurface;
}

rg_NUBSplineSurface3D rg_NUBSplineSurface3D::derivativeSurfaceOfUV() const
{
	rg_INT row = rg_BSplineSurface3D::getRowOfControlNet();
	rg_INT col = rg_BSplineSurface3D::getColumnOfControlNet();

	rg_INT uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();
	rg_INT vOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfV();
    
    rg_NUBSplineSurface3D derivativeSurface( row-1, col-1, uOrder-1, vOrder-1 );

    rg_Point3D newCtrlPt;
    rg_INT i = 0;
	for (i=0; i<row-1; i++)
    {
        for (rg_INT j=0; j<col-1; j++)
        {
            newCtrlPt 
                = ( uOrder-1 )*( vOrder - 1) 
                  * (getPointOnControlNet(i+1,j+1) - getPointOnControlNet(i+1,j)
                    - getPointOnControlNet(i, j+1) + getPointOnControlNet(i,j))
                  / ( (u_knot_vector[i+uOrder] - u_knot_vector[i+1])
                      * (v_knot_vector[j+vOrder] - v_knot_vector[j+1]) );

            derivativeSurface.setPointOnControlNet(i, j, newCtrlPt);
        }
    }     

    derivativeSurface.u_knot_vector = new rg_REAL [row + uOrder - 2];
    for (i=0; i<(row+uOrder-2); i++)
        derivativeSurface.u_knot_vector[i] = u_knot_vector[i+1];

    derivativeSurface.v_knot_vector = new rg_REAL [col + vOrder - 2]; 
    for (rg_INT j=0; j<(col+vOrder-2); j++)
        derivativeSurface.v_knot_vector[j] = v_knot_vector[j+1];

    return derivativeSurface;
}

rg_INT rg_NUBSplineSurface3D::findTheCorrespondingKnotSpanInU(const rg_REAL& insertingKnot)
{
    rg_INT n = getRowOfControlNet() - 1;
	rg_INT p = getOrderOfU() - 1;

	if(rg_EQ(insertingKnot, u_knot_vector[n + 1]))
		return n;

	rg_INT low = p;
	rg_INT high = n + 1;
	rg_INT mid = (low + high) / 2;

	while(rg_LT(insertingKnot, u_knot_vector[mid]) || rg_GE(insertingKnot, u_knot_vector[mid + 1]))
	{
		if(rg_LT(insertingKnot, u_knot_vector[mid]))
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

rg_INT rg_NUBSplineSurface3D::findTheCorrespondingKnotSpanInV(const rg_REAL& insertingKnot)
{
    rg_INT m = getColumnOfControlNet() - 1;
	rg_INT q = getOrderOfV() - 1;

	if(rg_EQ(insertingKnot, v_knot_vector[m + 1]))
		return m;

	rg_INT low = q;
	rg_INT high = m + 1;
	rg_INT mid = (low + high) / 2;

	while(rg_LT(insertingKnot, v_knot_vector[mid]) || rg_GE(insertingKnot, v_knot_vector[mid + 1]))
	{
		if(rg_LT(insertingKnot, v_knot_vector[mid]))
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

void rg_NUBSplineSurface3D::knotRefinement(rg_REAL* &insertingKnotValuesInU, const rg_INT& numOfinsertingKnotValuesInU,
	                rg_REAL* &insertingKnotValuesInV, const rg_INT& numOfinsertingKnotValuesInV)
{
	if(numOfinsertingKnotValuesInU ==0 || numOfinsertingKnotValuesInV == 0)
		return;

	rg_Point3D** ctrlNet = getControlNet();

    rg_INT n = getRowOfControlNet() - 1;
	rg_INT p = getOrderOfU() - 1;
	rg_INT r = n + p + 1;

	rg_INT m = getColumnOfControlNet() - 1;
	rg_INT q = getOrderOfV() - 1;
	rg_INT s = m + q + 1;

	// Knot refinement process for U-direction

	rg_INT a = findTheCorrespondingKnotSpanInU(insertingKnotValuesInU[ 0 ]);
	rg_INT b = findTheCorrespondingKnotSpanInU(insertingKnotValuesInU[numOfinsertingKnotValuesInU - 1]) + 1;

	// construct new knot vector in U-direction
	// and copy unchanged knot values

	rg_REAL* newKnotVectorInU = new rg_REAL[r + 1 + numOfinsertingKnotValuesInU];
	rg_INT i = 0;
	for(i = 0;i <= a;i++)
		newKnotVectorInU[ i ] = u_knot_vector[ i ];
	for(i = b + p;i <= r;i++)
		newKnotVectorInU[i + numOfinsertingKnotValuesInU] = u_knot_vector[ i ];

	// construct new control points
	// and copy unchanged control points

	rg_INT numOfRow = n + 1 + numOfinsertingKnotValuesInU;	

	rg_Point3D** newCtrlNetU = new rg_Point3D* [numOfRow];
	for(i = 0;i < numOfRow;i++)
		newCtrlNetU[ i ] = new rg_Point3D [m + 1];

	rg_INT row = 0;
	for(row = 0;row <= m;row++)
	{
		rg_INT j = 0;
		for(j = 0;j <= a - p;j++)
			newCtrlNetU[ j ][row] = ctrlNet[ j ][row];
		for(j = b - 1;j <= n;j++)
			newCtrlNetU[j + numOfinsertingKnotValuesInU][row] = ctrlNet[ j ][row];
	}

	// computing the changed control points
	// and inserting new knot values

	i = b + p - 1;
	rg_INT k = b + p + numOfinsertingKnotValuesInU - 1;

	rg_INT j = 0;
	for(j = numOfinsertingKnotValuesInU - 1;j >= 0;j--)
	{
		while(rg_LE(insertingKnotValuesInU[ j ], u_knot_vector[ i ]) && i > a)
		{
			for(row = 0;row <= m;row++)
				newCtrlNetU[k - p - 1][row] = ctrlNet[i - p - 1][row];
			newKnotVectorInU[ k ] = u_knot_vector[ i ];
			k--;
			i--;
		}
		for(row = 0;row <= m;row++)
			newCtrlNetU[k - p - 1][row] = newCtrlNetU[k - p][row];
		for(rg_INT l = 1;l <= p;l++)
		{
			rg_REAL alpha = newKnotVectorInU[k + l] - insertingKnotValuesInU[ j ];

			if(rg_ZERO(alpha))
			{
				for(row = 0;row <= m;row++)
					newCtrlNetU[k - p + l - 1][row] = newCtrlNetU[k - p + l][row];
			}
			else
			{
				alpha = alpha / (newKnotVectorInU[k + l] - u_knot_vector[i - p + l]);
				for(row = 0;row <= m;row++)
				{					
				    newCtrlNetU[k - p + l - 1][row] = alpha * newCtrlNetU[k - p + l - 1][row] 
     			                            + (1.0 - alpha) * newCtrlNetU[k - p + l][row];
				}
			}
		}
		newKnotVectorInU[ k ] = insertingKnotValuesInU[ j ];
		k--;
	}

	n = numOfRow - 1;
	r = numOfRow - 1 + p + 1;

	// Knot refinement process for V-direction

	a = findTheCorrespondingKnotSpanInV(insertingKnotValuesInV[ 0 ]);
	b = findTheCorrespondingKnotSpanInV(insertingKnotValuesInV[numOfinsertingKnotValuesInV - 1]) + 1;

	// construct new knot vector in V-direction
	// and copy unchanged knot values

	rg_REAL* newKnotVectorInV = new rg_REAL[s + 1 + numOfinsertingKnotValuesInV];
	for(i = 0;i <= a;i++)
		newKnotVectorInV[ i ] = v_knot_vector[ i ];
	for(i = b + q;i <= s;i++)
		newKnotVectorInV[i + numOfinsertingKnotValuesInV] = v_knot_vector[ i ];

	// construct new control points
	// and copy unchanged control points

	rg_INT numOfCol = m + 1 + numOfinsertingKnotValuesInV;

	rg_Point3D** newCtrlNet = new rg_Point3D* [numOfRow];
	for(i = 0;i < numOfRow;i++)
		newCtrlNet[ i ] = new rg_Point3D [numOfCol];

	rg_INT col = 0;
	for(col = 0;col <= n;col++)
	{
		rg_INT j = 0;
		for(j = 0;j <= a - q;j++)
			newCtrlNet[col][ j ] = newCtrlNetU[col][ j ];
		for(j = b - 1;j <= m;j++)
			newCtrlNet[col][j + numOfinsertingKnotValuesInV] = newCtrlNetU[col][ j ];
	}

	// computing the changed control points
	// and inserting new knot values

	i = b + q - 1;
	k = b + q + numOfinsertingKnotValuesInV - 1;

	for(j = numOfinsertingKnotValuesInV - 1;j >= 0;j--)
	{
		while(rg_LE(insertingKnotValuesInV[ j ], v_knot_vector[ i ]) && i > a)
		{
			for(col = 0;col <= n;col++)
				newCtrlNet[col][k - q - 1] = newCtrlNetU[col][i - q - 1];
			newKnotVectorInV[ k ] = v_knot_vector[ i ];
			k--;
			i--;
		}
		for(col = 0;col <= n;col++)
			newCtrlNet[col][k - q - 1] = newCtrlNet[col][k - q];
		for(rg_INT l = 1;l <= q;l++)
		{
			rg_REAL alpha = newKnotVectorInV[k + l] - insertingKnotValuesInV[ j ];

			if(rg_ZERO(alpha))
			{
				for(col = 0;col <= n;col++)
					newCtrlNet[col][k - q + l - 1] = newCtrlNet[col][k - q + l];
			}
			else
			{
				alpha = alpha / (newKnotVectorInV[k + l] - v_knot_vector[i - q + l]);
				for(col = 0;col <= n;col++)
				{					
				    newCtrlNet[col][k - q + l - 1] = alpha * newCtrlNet[col][k - q + l - 1] 
					                       + (1.0 - alpha) * newCtrlNet[col][k - q + l];
				}
			}
		}
		newKnotVectorInV[ k ] = insertingKnotValuesInV[ j ];
		k--;
	}

	// set new knot vector and control points

	delete[] u_knot_vector;
	delete[] v_knot_vector;
	for(i = 0;i < numOfRow;i++)
		delete[] newCtrlNetU[ i ];

	u_knot_vector = newKnotVectorInU;
	v_knot_vector = newKnotVectorInV;

	setControlNet(numOfRow, numOfCol, newCtrlNet);
}

void rg_NUBSplineSurface3D::knotRefinementOfU(rg_REAL* &insertingKnotValuesInU, const rg_INT& numOfinsertingKnotValuesInU)
{
	if(numOfinsertingKnotValuesInU == 0)
		return;

	rg_Point3D** ctrlNet = getControlNet();

    rg_INT n = getRowOfControlNet() - 1;
	rg_INT p = getOrderOfU() - 1;
	rg_INT r = n + p + 1;

	rg_INT m = getColumnOfControlNet() - 1;

	// Knot refinement process for U-direction

	rg_INT a = findTheCorrespondingKnotSpanInU(insertingKnotValuesInU[ 0 ]);
	rg_INT b = findTheCorrespondingKnotSpanInU(insertingKnotValuesInU[numOfinsertingKnotValuesInU - 1]) + 1;

	// construct new knot vector in U-direction
	// and copy unchanged knot values

	rg_REAL* newKnotVectorInU = new rg_REAL[r + 1 + numOfinsertingKnotValuesInU];
	rg_INT i = 0;
	for(i = 0;i <= a;i++)
		newKnotVectorInU[ i ] = u_knot_vector[ i ];
	for(i = b + p;i <= r;i++)
		newKnotVectorInU[i + numOfinsertingKnotValuesInU] = u_knot_vector[ i ];

	// construct new control points
	// and copy unchanged control points

	rg_INT numOfRow = n + 1 + numOfinsertingKnotValuesInU;	

	rg_Point3D** newCtrlNetU = new rg_Point3D* [numOfRow];
	for(i = 0;i < numOfRow;i++)
		newCtrlNetU[ i ] = new rg_Point3D [m + 1];

	rg_INT row = 0;
	for(row = 0;row <= m;row++)
	{
		rg_INT j = 0;
		for(j = 0;j <= a - p;j++)
			newCtrlNetU[ j ][row] = ctrlNet[ j ][row];
		for(j = b - 1;j <= n;j++)
			newCtrlNetU[j + numOfinsertingKnotValuesInU][row] = ctrlNet[ j ][row];
	}

	// computing the changed control points
	// and inserting new knot values

	i = b + p - 1;
	rg_INT k = b + p + numOfinsertingKnotValuesInU - 1;

	for(rg_INT j = numOfinsertingKnotValuesInU - 1;j >= 0;j--)
	{
		while(rg_LE(insertingKnotValuesInU[ j ], u_knot_vector[ i ]) && i > a)
		{
			for(row = 0;row <= m;row++)
				newCtrlNetU[k - p - 1][row] = ctrlNet[i - p - 1][row];
			newKnotVectorInU[ k ] = u_knot_vector[ i ];
			k--;
			i--;
		}
		for(row = 0;row <= m;row++)
			newCtrlNetU[k - p - 1][row] = newCtrlNetU[k - p][row];
		for(rg_INT l = 1;l <= p;l++)
		{
			rg_REAL alpha = newKnotVectorInU[k + l] - insertingKnotValuesInU[ j ];

			if(rg_ZERO(alpha))
			{
				for(row = 0;row <= m;row++)
					newCtrlNetU[k - p + l - 1][row] = newCtrlNetU[k - p + l][row];
			}
			else
			{
				alpha = alpha / (newKnotVectorInU[k + l] - u_knot_vector[i - p + l]);
				for(row = 0;row <= m;row++)
				{					
				    newCtrlNetU[k - p + l - 1][row] = alpha * newCtrlNetU[k - p + l - 1][row] 
     			                            + (1.0 - alpha) * newCtrlNetU[k - p + l][row];
				}
			}
		}
		newKnotVectorInU[ k ] = insertingKnotValuesInU[ j ];
		k--;
	}
	// set new knot vector and control points

	delete[] u_knot_vector;

	u_knot_vector = newKnotVectorInU;

	setControlNet(numOfRow, m + 1, newCtrlNetU);
}

void rg_NUBSplineSurface3D::knotRefinementOfV(rg_REAL* &insertingKnotValuesInV, const rg_INT& numOfinsertingKnotValuesInV)
{
	if(numOfinsertingKnotValuesInV == 0)
		return;

	rg_Point3D** ctrlNet = getControlNet();

    rg_INT n = getRowOfControlNet() - 1;
	
	rg_INT m = getColumnOfControlNet() - 1;
	rg_INT q = getOrderOfV() - 1;
	rg_INT s = m + q + 1;

	// Knot refinement process for V-direction

	rg_INT a = findTheCorrespondingKnotSpanInV(insertingKnotValuesInV[ 0 ]);
	rg_INT b = findTheCorrespondingKnotSpanInV(insertingKnotValuesInV[numOfinsertingKnotValuesInV - 1]) + 1;

	// construct new knot vector in V-direction
	// and copy unchanged knot values

	rg_REAL* newKnotVectorInV = new rg_REAL[s + 1 + numOfinsertingKnotValuesInV];
	rg_INT i = 0;
	for(i = 0;i <= a;i++)
		newKnotVectorInV[ i ] = v_knot_vector[ i ];
	for(i = b + q;i <= s;i++)
		newKnotVectorInV[i + numOfinsertingKnotValuesInV] = v_knot_vector[ i ];

	// construct new control points
	// and copy unchanged control points

	rg_INT numOfCol = m + 1 + numOfinsertingKnotValuesInV;

	rg_Point3D** newCtrlNetV = new rg_Point3D* [n + 1];
	for(i = 0;i < n + 1;i++)
		newCtrlNetV[ i ] = new rg_Point3D [numOfCol];

	rg_INT col = 0;
	for(col = 0;col <= n;col++)
	{
		rg_INT j = 0;
		for(j = 0;j <= a - q;j++)
			newCtrlNetV[col][ j ] = ctrlNet[col][ j ];
		for(j = b - 1;j <= m;j++)
			newCtrlNetV[col][j + numOfinsertingKnotValuesInV] = ctrlNet[col][ j ];
	}

	// computing the changed control points
	// and inserting new knot values

	i = b + q - 1;
	rg_INT k = b + q + numOfinsertingKnotValuesInV - 1;

	for(rg_INT j = numOfinsertingKnotValuesInV - 1;j >= 0;j--)
	{
		while(rg_LE(insertingKnotValuesInV[ j ], v_knot_vector[ i ]) && i > a)
		{
			for(col = 0;col <= n;col++)
				newCtrlNetV[col][k - q - 1] = ctrlNet[col][i - q - 1];
			newKnotVectorInV[ k ] = v_knot_vector[ i ];
			k--;
			i--;
		}
		for(col = 0;col <= n;col++)
			newCtrlNetV[col][k - q - 1] = newCtrlNetV[col][k - q];
		for(rg_INT l = 1;l <= q;l++)
		{
			rg_REAL alpha = newKnotVectorInV[k + l] - insertingKnotValuesInV[ j ];

			if(rg_ZERO(alpha))
			{
				for(col = 0;col <= n;col++)
					newCtrlNetV[col][k - q + l - 1] = newCtrlNetV[col][k - q + l];
			}
			else
			{
				alpha = alpha / (newKnotVectorInV[k + l] - v_knot_vector[i - q + l]);
				for(col = 0;col <= n;col++)
				{					
				    newCtrlNetV[col][k - q + l - 1] = alpha * newCtrlNetV[col][k - q + l - 1] 
					                       + (1.0 - alpha) * newCtrlNetV[col][k - q + l];
				}
			}
		}
		newKnotVectorInV[ k ] = insertingKnotValuesInV[ j ];
		k--;
	}

	// set new knot vector and control points

	delete[] v_knot_vector;

	v_knot_vector = newKnotVectorInV;

	setControlNet(n + 1, numOfCol, newCtrlNetV);
}


// Assume that the multiplicity of interior knot = 1.

rg_BzSurface3D** rg_NUBSplineSurface3D::decomposeSurfaceIntoBezierPatches() const
{
	rg_DEGREE p = getOrderOfU() - 1;
	rg_DEGREE q = getOrderOfV() - 1;
	rg_DEGREE n = getRowOfControlNet() - 1;
	rg_DEGREE m = getColumnOfControlNet() - 1;
	//rg_DEGREE r = n + p + 1;
	//rg_DEGREE s = m + q + 1;

	// compute the number of Bezier surface patches to be generated
	
	rg_INT numOfKnotSpanInU = getNumOfNonZeroLengthKnotSpanOfKnotVectorU(); // the number of nonzero length knot span in direction U
    rg_INT numOfInteriorKnotInU = numOfKnotSpanInU-1;

	rg_INT numOfKnotSpanInV = getNumOfNonZeroLengthKnotSpanOfKnotVectorV(); // the number of nonzero length knot span in direction V
    rg_INT numOfInteriorKnotInV = numOfKnotSpanInV-1;

	//rg_INT numOfBezierPathes = numOfKnotSpanInU * numOfKnotSpanInV;

	//rg_Point3D*** ctrlNetOfBezierPatches = new rg_Point3D**[numOfBezierPathes];
	/*
	for(rg_INT i = 0;i < numOfBezierPathes;i++)
	{
		ctrlNetOfBezierPatches[ i ] = new rg_Point3D*[p + 1];
		for(rg_INT j = 0;j <= p;j++)
		{
			ctrlNetOfBezierPatches[ i ][ j ] = new rg_Point3D[q + 1];
		}
	}
	*/
	
	rg_INT numOfBezierStrips = numOfKnotSpanInV;
	rg_Point3D*** ctrlNetOfBezierStrips = new rg_Point3D** [numOfBezierStrips];

	rg_INT i = 0;
	for(i = 0;i < numOfBezierStrips;i++)
	{
		ctrlNetOfBezierStrips[ i ] = new rg_Point3D*[q + 1];
		for(rg_INT j = 0;j <= p;j++)
		{
			ctrlNetOfBezierStrips[ i ][ j ] = new rg_Point3D[n + 1];
		}
	}

	// U-direction
	rg_INT a = p;
	rg_INT b = p + 1;
	rg_INT indexOfBezierStrips = 0;

	rg_Point3D** control_net = getControlNet();
	for(i = 0;i <= p;i++)
		for(rg_INT row = 0;row <= m;row++)
			ctrlNetOfBezierStrips[indexOfBezierStrips][ i ][row] = control_net[ i ][row];

	while(b < m)
	{
		i = b;
		while(b < m && rg_EQ(u_knot_vector[b + 1], u_knot_vector[ b ])) b++;
		rg_INT mult = b - i + 1;
		if(mult < p)
		{
			// compute and store alphas
			rg_REAL numer = u_knot_vector[ b ] - u_knot_vector[ a ];
			rg_REAL* alphas = new rg_REAL[p - mult];
			rg_INT j = 0;
			for(j = p;j > mult;j--)
				alphas[j - mult - 1] = numer / (u_knot_vector[a + j] - u_knot_vector[ a ]);
			rg_DEGREE numOfKnotInsertion = p - mult;
			for(j = 1;j <= numOfKnotInsertion;j++)
			{
				rg_INT save = numOfKnotInsertion - j;
				rg_INT s = mult + j;
				for(rg_INT k = p;k >= s;k--)
				{
					rg_REAL alpha = alphas[k - s];
					for(rg_INT row = 0;row <= m;row++)
						ctrlNetOfBezierStrips[indexOfBezierStrips][ k ][row]
							= alpha * ctrlNetOfBezierStrips[indexOfBezierStrips][ k ][row] + (1.0 - alpha) * ctrlNetOfBezierStrips[indexOfBezierStrips][k - 1][row];
				}
				if(b < m)
					for(rg_INT row = 0;row <= m;row++)
						ctrlNetOfBezierStrips[indexOfBezierStrips + 1][save][row]
							= ctrlNetOfBezierStrips[indexOfBezierStrips][ p ][row];
			}
			delete[] alphas;
		}
		indexOfBezierStrips++;
		if(b < m)
		{
			for(i = p - mult;i <= p;i++)
				for(rg_INT row = 0;row <= m;row++)
					ctrlNetOfBezierStrips[indexOfBezierStrips][ i ][row] = control_net[b - p + i][row];
			a = b;
			b++;
		}
	}

	// V-direction
	
	rg_INT numOfBezierPathes = numOfKnotSpanInU * numOfKnotSpanInV;

	rg_Point3D*** ctrlNetOfBezierPatches = new rg_Point3D**[numOfBezierPathes];

	for(i = 0;i < numOfBezierPathes;i++)
	{
		ctrlNetOfBezierPatches[ i ] = new rg_Point3D*[p + 1];
		for(rg_INT j = 0;j <= p;j++)
		{
			ctrlNetOfBezierPatches[ i ][ j ] = new rg_Point3D[q + 1];
		}
	}
	
	a = q;
	b = q + 1;
	rg_INT indexOfBezierPatches = 0;

	n = p * (numOfBezierStrips - 1) + (p + 1) * 1 - 1;
	
	indexOfBezierStrips = 0;
	for(i = 0;i <= q;i++)
		for(rg_INT col = 0;col <= n;col++)
		{
			ctrlNetOfBezierPatches[indexOfBezierPatches][col][ i ] = ctrlNetOfBezierStrips[indexOfBezierStrips][col][i % (p + 1)];
			if((i % (p + 1)) == p)
				indexOfBezierStrips++;
		}

	while(b < n)
	{
		i = b;
		while(b < n && rg_EQ(v_knot_vector[b + 1], v_knot_vector[ b ])) b++;
		rg_INT mult = b - i + 1;
		if(mult < q)
		{
			// compute and store alphas
			rg_REAL numer = v_knot_vector[ b ] - v_knot_vector[ a ];
			rg_REAL* alphas = new rg_REAL[q - mult];
			rg_INT j = 0;
			for(j = q;j > mult;j--)
				alphas[j - mult - 1] = numer / (v_knot_vector[a + j] - v_knot_vector[ a ]);
			rg_DEGREE numOfKnotInsertion = q - mult;
			for(j = 1;j <= numOfKnotInsertion;j++)
			{
				rg_INT save = numOfKnotInsertion - j;
				rg_INT s = mult + j;
				for(rg_INT k = q;k >= s;k--)
				{
					rg_REAL alpha = alphas[k - s];
					for(rg_INT col = 0;col <= n;col++)
						ctrlNetOfBezierPatches[indexOfBezierPatches][col][ k ]
							= alpha * ctrlNetOfBezierPatches[indexOfBezierPatches][col][ k ] + (1.0 - alpha) * ctrlNetOfBezierPatches[indexOfBezierPatches][col][k - 1];
				}
				if(b < n)
					for(rg_INT col = 0;col <= n;col++)
						ctrlNetOfBezierPatches[indexOfBezierPatches + 1][col][save]
							= ctrlNetOfBezierPatches[indexOfBezierPatches][col][ q ];
			}
			delete[] alphas;
		}
		indexOfBezierPatches++;
		if(b < n)
		{
			indexOfBezierStrips = 0;
			for(i = q - mult;i <= q;i++)
				for(rg_INT col = 0;col <= n;col++)
				{
					ctrlNetOfBezierPatches[indexOfBezierPatches][col][ i ] = ctrlNetOfBezierStrips[indexOfBezierStrips][col][b - q + i];
					if((col % (p + 1)) == p)
						indexOfBezierStrips++;
				}
			a = b;
			b++;
		}
	}


	// construct Bezier surface patches
	rg_BzSurface3D** BezierSufaceList = new rg_BzSurface3D* [numOfKnotSpanInU];
	for(i = 0;i < numOfKnotSpanInU;i++)
		BezierSufaceList[ i ] = new rg_BzSurface3D [numOfKnotSpanInV];

	indexOfBezierPatches = 0;
	for(i = 0;i < numOfKnotSpanInU;i++)
	{
		for(rg_INT j = 0;j < numOfKnotSpanInV;j++)
		{
			BezierSufaceList[ i ][ j ].setDegree(p, q);
			for(rg_INT k = 0;k <= p;k++)
				for(rg_INT l = 0;l <= q;l++)
					BezierSufaceList[ i ][ j ].setCtrlPt(k, l, ctrlNetOfBezierPatches[indexOfBezierPatches][ k ][ l ]);
			indexOfBezierPatches++;
		}
	}

	// free dynamic memorry allocations
	
	for(i = 0;i < m;i++)
	{		
		for(rg_INT j = 0;j <= n;j++)
		{			
			delete[] ctrlNetOfBezierPatches[ i ][ j ];
		}
		delete[] ctrlNetOfBezierPatches[ i ];
	}
	delete[] ctrlNetOfBezierPatches;

	return BezierSufaceList;
}

rg_BzSurface3D** rg_NUBSplineSurface3D::decomposeSurfaceIntoBezierPatchesUsingKnotRefinement() const
{
	rg_DEGREE p = getOrderOfU() - 1;
	rg_DEGREE q = getOrderOfV() - 1;
	rg_DEGREE n = getRowOfControlNet() - 1;
	rg_DEGREE m = getColumnOfControlNet() - 1;
	//rg_DEGREE r = n + p + 1;
	rg_DEGREE sizeOfKnotVectorU = n + p + 1 + 1;
	//rg_DEGREE s = m + q + 1;
	rg_DEGREE sizeOfKnotVectorV = m + q + 1 + 1;

	// compute the number of Bezier surface patches to be generated
	
	rg_INT numOfKnotSpanInU = getNumOfNonZeroLengthKnotSpanOfKnotVectorU(); // the number of nonzero length knot span in direction U
    rg_INT numOfInteriorKnotInU = numOfKnotSpanInU-1;

	rg_INT numOfKnotSpanInV = getNumOfNonZeroLengthKnotSpanOfKnotVectorV(); // the number of nonzero length knot span in direction V
    rg_INT numOfInteriorKnotInV = numOfKnotSpanInV-1;

	rg_INT numOfBezierPathes = numOfKnotSpanInU * numOfKnotSpanInV;

	rg_NUBSplineSurface3D duplicatedSurface = (*this);

	rg_REAL** multiplicityU = new rg_REAL* [numOfInteriorKnotInU];
	rg_REAL** multiplicityV = new rg_REAL* [numOfInteriorKnotInV];

	rg_INT i = 0;
	for(i = 0;i < numOfInteriorKnotInU;i++)
		multiplicityU[ i ] = new rg_REAL[ 2 ];
	for(i = 0;i < numOfInteriorKnotInV;i++)
		multiplicityV[ i ] = new rg_REAL[ 2 ];

    i = 0;               //  i traces the array multiplicity.
    rg_INT j = p + 1;    //  j traces knotVector.
	rg_INT numOfInsertingKnotValuesU = 0;
	while(i < numOfInteriorKnotInU)
	{
		rg_INT s = 1;
		while( (j < sizeOfKnotVectorU - 1 - p) && rg_EQ(u_knot_vector[j], u_knot_vector[j+1]))
		{
			s++;
			j++;
		}
		multiplicityU[i][ 0 ] = u_knot_vector[j];   // replicated knot value in the interior knot span
		multiplicityU[i][ 1 ] = s;               // multiplicity of replicated knot value in the interior knot span
		numOfInsertingKnotValuesU += (p - s);
		i++;
        j++;
	}

	i = 0;               //  i traces the array multiplicity.
    j = q + 1;    //  j traces knotVector.
	rg_INT numOfInsertingKnotValuesV = 0;
	while(i < numOfInteriorKnotInV)
	{
		rg_INT s = 1;
		while( (j < sizeOfKnotVectorV - 1 - q) && rg_EQ(v_knot_vector[j], v_knot_vector[j+1]))
		{
			s++;
			j++;
		}
		multiplicityV[i][ 0 ] = v_knot_vector[j];   // replicated knot value in the interior knot span
		multiplicityV[i][ 1 ] = s;                  // multiplicity of replicated knot value in the interior knot span
		numOfInsertingKnotValuesV += (q - s);
		i++;
        j++;
	}

    rg_REAL* insertingKnotValuesU = new rg_REAL [numOfInsertingKnotValuesU];
	i = 0;
	rg_INT indexOfinsertingKnotValuesU = 0;
	while((i < numOfInteriorKnotInU) && (indexOfinsertingKnotValuesU < numOfInsertingKnotValuesU))
	{
		rg_INT s = 1;
		while(s <= p - multiplicityU[ i ][ 1 ])
		{
			insertingKnotValuesU[indexOfinsertingKnotValuesU] = multiplicityU[ i ][ 0 ];
			indexOfinsertingKnotValuesU++;
			s++;
		}
		i++;
	}

    rg_REAL* insertingKnotValuesV = new rg_REAL [numOfInsertingKnotValuesV];
	i = 0;
	rg_INT indexOfinsertingKnotValuesV = 0;
	while((i < numOfInteriorKnotInV) && (indexOfinsertingKnotValuesV < numOfInsertingKnotValuesV))
	{
		rg_INT s = 1;
		while(s <= q - multiplicityV[ i ][ 1 ])
		{
			insertingKnotValuesV[indexOfinsertingKnotValuesV] = multiplicityV[ i ][ 0 ];
			indexOfinsertingKnotValuesV++;
			s++;
		}
		i++;
	}

	//if(numOfInsertingKnotValuesU > 0 || numOfInsertingKnotValuesV >0)
	if(numOfInsertingKnotValuesU > 0 || numOfInsertingKnotValuesV > 0)
		duplicatedSurface.knotRefinement(insertingKnotValuesU, numOfInsertingKnotValuesU, insertingKnotValuesV, numOfInsertingKnotValuesV);
		

	// construct Bezier surface patches
	rg_BzSurface3D** BezierSufaceList = new rg_BzSurface3D* [numOfKnotSpanInU];
	for(i = 0;i < numOfKnotSpanInU;i++)
		BezierSufaceList[ i ] = new rg_BzSurface3D [numOfKnotSpanInV];

	for(i = 0;i < numOfKnotSpanInU;i++)
		for(rg_INT j = 0;j < numOfKnotSpanInV;j++)
		{
			BezierSufaceList[ i ][ j ].setDegree(p, q);
			for(rg_INT k = 0;k <= p;k++)
				for(rg_INT l = 0;l <= q;l++)
					BezierSufaceList[ i ][ j ].setCtrlPt(k, l, duplicatedSurface.getPointOnControlNet(i * p + k, j * q + l));
		}

	for(i = 0;i < numOfInteriorKnotInU;i++)
		delete [] multiplicityU[ i ];
	delete[] multiplicityU;
	
	for(i = 0;i < numOfInteriorKnotInV;i++)
		delete [] multiplicityV[ i ];
	delete[] multiplicityV;

	delete[] insertingKnotValuesU;
	delete[] insertingKnotValuesV;

	return BezierSufaceList;
}

////////////////////////////////////////////////////
////	Skinned rg_Surface.						////---------------------------
////////////////////////////////////////////////////

void rg_NUBSplineSurface3D::makeSkinnedSurface( const rg_INT         &noOfSectionC, 
                                             rg_NUBSplineCurve3D* sectionCurves ) 
{
	rg_sListByPtr* controlPolygons = new rg_sListByPtr [noOfSectionC];
	
	setUOrdersOfSkinnedSurface(noOfSectionC, sectionCurves, controlPolygons);

	rg_INT numberOfNewKnot = 0;

	setNewKnotVectorInSkinnedSurface(noOfSectionC, sectionCurves, numberOfNewKnot);
	
	modifyControlPolygonsWithKnotInsertion(numberOfNewKnot, noOfSectionC, 
		                                    sectionCurves, controlPolygons);

	makeControlNetForSkinnedSurface(numberOfNewKnot, noOfSectionC, controlPolygons);

	delete [] controlPolygons;
}

void rg_NUBSplineSurface3D::setUOrdersOfSkinnedSurface( const rg_INT         &noOfSectionC,
                                                     rg_NUBSplineCurve3D* sectionCurves, 
                                                     rg_sListByPtr*            modifiedControlPolygons )
{
	rg_INT	noOfControlPointSectionC = 0;

	rg_INT i = 0;
	for (i=0; i<noOfSectionC; i++)
	{
		noOfControlPointSectionC = sectionCurves[i].getNumOfCtrlPts();
		for (rg_INT j=0; j<noOfControlPointSectionC; j++)
		{
			rg_Point3D* controlP = new rg_Point3D( sectionCurves[i].getCtrlPt(j) );
			modifiedControlPolygons[i].append(controlP);
		}
	}
	
	rg_BSplineSurface3D::setOrderOfU(sectionCurves[0].getOrder());
	rg_BSplineSurface3D::setOrderOfV(4);
}

void rg_NUBSplineSurface3D::setNewKnotVectorInSkinnedSurface( const rg_INT         &noOfSectionC, 
                                                           rg_NUBSplineCurve3D* sectionCurves, 
                                                           rg_INT               &noOfNewKnot )
{
	rg_dListByPtr newKnotVector;
	rg_REAL* tempKnot;

	rg_REAL* noOfKnotVector = new rg_REAL [noOfSectionC];
	rg_INT i = 0;
	for (i=0; i<noOfSectionC; i++)
	{
		noOfKnotVector[i] = sectionCurves[i].getNumOfCtrlPts() 
							+ sectionCurves[i].getOrder();
	}

	for (i=0; i<noOfKnotVector[0]; i++)
	{
		tempKnot  = new rg_REAL;
		*tempKnot = sectionCurves[0].getKnotValue(i);
		newKnotVector.append(tempKnot);
	}
	
	rg_INT    uOrder = (rg_INT) rg_BSplineSurface3D::getOrderOfU();
	rg_dNodeByPtr* curNode;
	rg_dNodeByPtr* nextNode;
	rg_REAL*  insertingValue;
	for (i=1; i<noOfSectionC; i++)
	{
		for (rg_INT j=uOrder; j<(noOfKnotVector[i]-uOrder); j++)
		{
			
			newKnotVector.reset();
			curNode = newKnotVector.getFirstNode();

			rg_INT  noOfnewKnotVector = newKnotVector.getNumberOfEntity();
			rg_REAL comparedKnot      = sectionCurves[i].getKnotValue(j);

			for (rg_INT k=0; k<(noOfnewKnotVector-uOrder); k++)
			{
				nextNode = curNode->getNext();

				rg_REAL* curEntity  = (rg_REAL*) curNode->getEntity();
				rg_REAL* nextEntity = (rg_REAL*) nextNode->getEntity();

				if ( rg_GT(comparedKnot, *curEntity) && rg_LT(comparedKnot, *nextEntity) )
				{
					insertingValue  = new rg_REAL;
					*insertingValue = comparedKnot;

					newKnotVector.insert(insertingValue, curNode);
					
					break;
				}	
/*				else if ( rg_EQ(*curEntity, comparedKnot) )
				{
		
					rg_INT noOfSameKnotInNewList = 1;
					rg_INT noOfSameKnotInOldList = 1;
					
					//		Find the number of the same knot Value in the new knot vector
					//		which will be the u_knot_vector of the skinned surface.			
					while ( curEntity = (rg_REAL*) nextNode->getEntity() )
					{
					if ( rg_EQ(*curEntity, comparedKnot) )
							noOfSameKnotInNewList++;
						else 
							break;
						curNode = nextNode;
					}

					//		Find the number of the same knot Value in the old knot vector
					//		which is the knot vector the ith section curve.
					for(rg_INT newj=j+1; newj<(noOfKnotVector[i]-uOrder); newj++)
					{
						if ( rg_EQ(*curEntity, sectionCurves[i].getKnotValue(newj)) )	
							noOfSameKnotInOldList++;
						else 
							break;
					}

					if (noOfSameKnotInNewList >= noOfSameKnotInOldList)
						i += noOfSameKnotInOldList;
					else
					{
						rg_REAL** insertingValues 
						             = new rg_REAL* [noOfSameKnotInOldList-noOfSameKnotInNewList];

						for(rg_INT l=0; l<(noOfSameKnotInOldList-noOfSameKnotInNewList); l++)
						{
							*insertingValues[l] = sectionCurves[i].getKnotValue(j); 

							newKnotVector.insert(insertingValues[l], curNode);
							curNode = curNode->getNext();
						}
					}

					break;
				}  //	else if ( rg_EQ(*curKnotOfNewList, ...
				else 
*/					curNode = nextNode;

			}  //	while ( curKnotOfNewList = ...
		}  //	for (rg_INT j=rg_BSplineSurface3D::getOrderOfU() ...
	}  //	for (i=1; i<noOfSectionC; ... 

	noOfNewKnot = newKnotVector.getNumberOfEntity();

	if ( u_knot_vector != rg_NULL )
		delete [] u_knot_vector;

	u_knot_vector = new rg_REAL [noOfNewKnot];

	rg_REAL* newKnot = rg_NULL;
	
	newKnotVector.reset();

	i = 0;
	while ( newKnot = (rg_REAL*) newKnotVector.getNextEntity() )
	{
		u_knot_vector[i] = *newKnot;
		i++;
	}

	delete [] noOfKnotVector;
}

void rg_NUBSplineSurface3D::knotInsertionInSkinnedSurface ( 
                                        rg_sListByPtr* curKnotVector, 
                                        rg_sListByPtr &curControlPolygon, 
                                        const rg_REAL &insertingKnot )
{
	curKnotVector->reset();
	curControlPolygon.reset();

	rg_INT  uOrder = (rg_INT)(rg_BSplineSurface3D::getOrderOfU());
	rg_REAL alpha  = 0.;

	rg_REAL* inserting     = rg_NULL;
	rg_REAL* curNodeValue  = rg_NULL;
	rg_REAL* nextNodeValue = rg_NULL;

	rg_sNodeByPtr* curNodeOfKnotV  = rg_NULL;
	rg_sNodeByPtr* nextNodeOfKnotV = rg_NULL;

	rg_sNodeByPtr* curNode		  = rg_NULL;
	rg_sNodeByPtr* influencedNode = rg_NULL;
	rg_sNodeByPtr* killNode		  = rg_NULL;

	rg_Point3D* curPoint = rg_NULL;
	rg_Point3D* oldPoint = rg_NULL;

	curNodeOfKnotV = curKnotVector->getFirstNode();

	rg_INT noOfCurKnotVector = curKnotVector->getNumberOfEntity();
	
	rg_INT influencedKnot = 0;
	for(influencedKnot=0; influencedKnot<noOfCurKnotVector; influencedKnot++)
	{
		curNodeValue  = (rg_REAL*) curNodeOfKnotV->getEntity();

		nextNodeOfKnotV = curNodeOfKnotV->getNext();
		nextNodeValue   = (rg_REAL*) nextNodeOfKnotV->getEntity();

		if ( rg_EQ(insertingKnot, *curNodeValue) )
			return;
		else if ( rg_GT(insertingKnot, *curNodeValue) && rg_LT(insertingKnot, *nextNodeValue) )
			break;
		else;
		
		curNodeOfKnotV = curNodeOfKnotV->getNext();
	}

	oldPoint = new rg_Point3D[uOrder];
	curNode  = curControlPolygon.getFirstNode();

	//	set the first location of the infleneced Points by knot insertion.
	rg_INT k = 0;
	for (k=0; k<(influencedKnot-(uOrder-1)); k++)
		curNode = curNode->getNext();
		
	//	get the original control point to find the new control point.
	influencedNode = curNode;
	for ( k=0; k<uOrder; k++)
	{
		curPoint       = (rg_Point3D*) influencedNode->getEntity();
		oldPoint[k]    = *curPoint;
		influencedNode = influencedNode->getNext();
	}

	//	kill original control point to be modified. 
	for (k=0; k<(uOrder-2); k++)
	{
		killNode = curNode->getNext();
		curControlPolygon.killNodeOnly(killNode);
	}

	//	insert the new control point.

	rg_INT noOfNewKnot = curKnotVector->getNumberOfEntity();

	rg_REAL* tempKnotVector = new rg_REAL [noOfNewKnot];
	rg_REAL* newKnot        = rg_NULL;

	rg_INT count = 0;
	curKnotVector->reset();

	while ( newKnot = (rg_REAL*) curKnotVector->getNextEntity() )
	{
		tempKnotVector[count] = *newKnot;
		count++;
	}

	rg_REAL x = 0.;
	rg_REAL y = 0.;
	rg_REAL z = 0.;
	for ( k=0; k<(uOrder-1); k++)
	{
		alpha = ( insertingKnot - tempKnotVector[influencedKnot - uOrder + 2 + k] )
				/ ( tempKnotVector[influencedKnot + 1 + k] 
				    - tempKnotVector[influencedKnot - uOrder + 2 + k] );

		x = alpha * oldPoint[k+1].getX() + (1-alpha) * oldPoint[k].getX();
		y = alpha * oldPoint[k+1].getY() + (1-alpha) * oldPoint[k].getY();
		z = alpha * oldPoint[k+1].getZ() + (1-alpha) * oldPoint[k].getZ();

		rg_Point3D* newPoint = new rg_Point3D(x, y, z);
		curControlPolygon.insert(newPoint, curNode);
		curNode = curNode->getNext();
	}

	delete [] tempKnotVector;
	delete [] oldPoint;

	inserting  = new rg_REAL;
	*inserting = insertingKnot;
	curKnotVector->insert(inserting, curNodeOfKnotV);
}

void rg_NUBSplineSurface3D::modifyControlPolygonsWithKnotInsertion( 
                                             const rg_INT        &noOfNewKnot, 
                                             const rg_INT        &noOfSectionC, 
                                             rg_NUBSplineCurve3D* sectionCurves, 
                                             rg_sListByPtr*            controlPolygons )
{
	rg_REAL* noOfKnotVector = new rg_REAL [noOfSectionC];
	rg_INT i = 0;
	for (i=0; i<noOfSectionC; i++)
	{
		noOfKnotVector[i] = sectionCurves[i].getNumOfCtrlPts() 
							+ sectionCurves[i].getOrder();
	}

	rg_INT	uOrder       = (rg_INT)(rg_BSplineSurface3D::getOrderOfU());
	rg_REAL* tempKnot = rg_NULL;

	rg_sListByPtr* tempOldKnotVector = new rg_sListByPtr;
	
	for (i=0; i<noOfSectionC; i++)
	{
		rg_INT j = 0;
		for (j=0; j<noOfKnotVector[i]; j++)
		{
			tempKnot  = new rg_REAL;
			*tempKnot = sectionCurves[i].getKnotValue(j);
			tempOldKnotVector->append(tempKnot);
		}

		for (j= uOrder; j<(noOfNewKnot-uOrder); j++)
		{
			knotInsertionInSkinnedSurface (tempOldKnotVector, controlPolygons[i], u_knot_vector[j]);			
		}

		tempOldKnotVector->killList();
	}

	delete tempOldKnotVector;
	delete [] noOfKnotVector;	
}

void rg_NUBSplineSurface3D::setKnotVectorInVDirection(	const rg_INT &noOfSectionC, 
                                                    rg_sListByPtr*    controlPolygons)
{
	rg_INT n = controlPolygons[0].getNumberOfEntity();
    rg_Point3D** ptSet4KnotVec = new rg_Point3D*[n];
    rg_INT i = 0;
	for(i=0; i<n; i++)
        ptSet4KnotVec[i] = new rg_Point3D[noOfSectionC];

    rg_INT j = 0;
	for(j=0; j<noOfSectionC; j++)
    {
		controlPolygons[j].reset();
        i = 0;
        rg_Point3D* pt;
        while( pt = (rg_Point3D*)controlPolygons[j].getNextEntity() )
        {
            ptSet4KnotVec[i][j] = (*pt);
            i++;
        }
    }

    rg_REAL** parameterization = new rg_REAL*[n];
    for(i=0; i<n; i++)
        parameterization[i] = rg_GeoFunc::chordLength(noOfSectionC, ptSet4KnotVec[i]);

    rg_REAL* averageParam = new rg_REAL[noOfSectionC];
    for (j=0; j<noOfSectionC; j++)
    {
        rg_REAL knotValue = 0.0;
        for (i=0; i<n; i++)
            knotValue += parameterization[i][j];

        averageParam[j] = knotValue/n;
    }

	rg_ORDER vOrder = rg_BSplineSurface3D::getOrderOfV();
    rg_REAL* knotVecV = new rg_REAL[noOfSectionC + 2*(vOrder-1)];
    for(j=0; j<(vOrder-1); j++)
    {
        knotVecV[j] = 0.0;
        knotVecV[noOfSectionC + 2*(vOrder - 1) - 1 - j] = 1.0;
    }


    for(j = vOrder-1; j<(noOfSectionC + vOrder - 1); j++)
        knotVecV[j] = averageParam[j-vOrder+1];

	setKnotVectorOfV(noOfSectionC + 2*(vOrder - 1), knotVecV);

    delete [] knotVecV;
    delete [] averageParam;

    for (i=0; i<n; i++)
        delete [] parameterization[i];
    delete [] parameterization;

    for (i=0; i<n; i++)
        delete [] ptSet4KnotVec[i];
    delete [] ptSet4KnotVec;    
}

//              |
//              |
//              |     0....0....0....0....0
//            ^ |
//          j | |  
//          v   |
//				|
//               ---------------------------------
//							--->
//                          i u
/////////////////////////////////////////////////////////////////
//          section curve,       u, i, row
//          v-directional curve, v, j, column
//
/////////////////////////////////////////////////////////////////
	
void rg_NUBSplineSurface3D::makeControlNetForSkinnedSurface( 
                                        const rg_INT &noOfKnotVector,
                                        const rg_INT &noOfSectionC, 
                                        rg_sListByPtr*     controlPolygons)
{
    //  row is the number of control points of one section curve.
	rg_INT row = noOfKnotVector - rg_BSplineSurface3D::getOrderOfU();

	//linked list reset.
	rg_INT i = 0;
	for (i=0; i<noOfSectionC; i++)
		controlPolygons[i].reset();
	
	rg_Point3D* tempPoint          = rg_NULL;
	rg_Point3D* tempControlPolygon = rg_NULL;

	rg_NUBSplineCurve3D* vSectionCurves = new rg_NUBSplineCurve3D [row];

	for (i=0; i<row; i++)
	{
		tempControlPolygon = new rg_Point3D[noOfSectionC];
	
		for (rg_INT j=0; j<noOfSectionC; j++)
		{
			tempPoint = (rg_Point3D*) controlPolygons[j].getNextEntity();
			tempControlPolygon[j] = *tempPoint;
		}
        rg_Point3D zeroTangent;
		vSectionCurves[i].interpolateWithEndCondition( noOfSectionC,
                                                       tempControlPolygon,
                                                       zeroTangent,
                                                       zeroTangent);

        delete [] tempControlPolygon;
	}

    //  col is the number of control points of one v-directional curve.
	rg_INT col = vSectionCurves[0].getNumOfCtrlPts();

	//	set control_net.
    rg_Point3D** controlNet = new rg_Point3D* [row];

    for ( i=0; i<row; i++)
    {
        controlNet[i] = new rg_Point3D[col];
        for (rg_INT j=0; j<col; j++)
            controlNet[i][j] = vSectionCurves[i].getCtrlPt(j);
    }

	rg_BSplineSurface3D::setControlNet(row, col, controlNet);

	//	set v_knot_vector.
	setKnotVectorInVDirection(noOfSectionC, controlPolygons);

    delete [] vSectionCurves;

}

////    Operator Overloading.--------------------------------------------------
//  April 7 1997 : made.
rg_NUBSplineSurface3D& rg_NUBSplineSurface3D::operator =(const rg_NUBSplineSurface3D &surface)
{
    if ( this == &surface )
        return *this;

    rg_Surface::setID( surface.rg_Surface::getID() );
    rg_Surface::setPlanarity( surface.rg_Surface::getPlanarity() );

    rg_INT nU = surface.rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT nV = surface.rg_BSplineSurface3D::getColumnOfControlNet();

    rg_INT uOrder = (rg_INT) surface.rg_BSplineSurface3D::getOrderOfU();
    rg_INT vOrder = (rg_INT) surface.rg_BSplineSurface3D::getOrderOfV();

    rg_BSplineSurface3D::setControlNet(nU, nV,
                                    surface.rg_BSplineSurface3D::getControlNet());

    rg_BSplineSurface3D::setOrderOfSurface(uOrder, vOrder);

    if (u_knot_vector != rg_NULL)
        delete [] u_knot_vector;
    if (v_knot_vector != rg_NULL)
        delete [] v_knot_vector;

    u_knot_vector = new rg_REAL [nU+uOrder];
    for (rg_INT i=0; i<nU+uOrder; i++)
        u_knot_vector[i] = surface.u_knot_vector[i];
    v_knot_vector = new rg_REAL [nV+vOrder];
    for (rg_INT j=0; j<nV+vOrder; j++)
        v_knot_vector[j] = surface.v_knot_vector[j];

    return *this;
}

///	File-Out Function.

//void rg_NUBSplineSurface3D::fileOut(const char* fileName)
//{
//	FILE* fout;
//
//	if ((fout = fopen(fileName, "a+")) == rg_NULL)
//		return;
//	
//	rg_INT u_closeState = 0;
//	rg_INT v_closeState = 0;
//	
//	rg_INT noOfUKnots = (rg_BSplineSurface3D::getRowOfControlNet()) 
//		             + (rg_BSplineSurface3D::getOrderOfU());
//	rg_INT noOfVKnots = (rg_BSplineSurface3D::getColumnOfControlNet()) 
//		             + (rg_BSplineSurface3D::getOrderOfV());
//	
////	rg_BSplineSurface3D::fileOut(fout, fileName);	
//	rg_BSplineSurface3D::fileOut(fileName);	
//
//	fprintf( fout, "u_closed\n" );
//	fprintf( fout, "%d\n", u_closeState );
//	fprintf( fout, "v_closed\n" );
//	fprintf( fout, "%d\n", v_closeState );
//
//	fprintf( fout, "u_knots\n" );
//	rg_INT i = 0;
//	for (i=0; i<noOfUKnots; i++)
//		fprintf( fout, "%9.8f ", u_knot_vector[i] );
//	
//	fprintf( fout, "\n" );
//
//	fprintf( fout, "v_knots\n" );
//	for (i=0; i<noOfVKnots; i++)
//		fprintf( fout, "%9.8f ", v_knot_vector[i] );
//	
//	fprintf( fout, "\n" );
//	
//	fclose( fout );
//}

//	Conversion between power spline & b-spline form.
void rg_NUBSplineSurface3D::powerSplineToNUBSplineSurface( const rg_DEGREE& uDegree,
													    const rg_DEGREE& vDegree,
														const rg_INT& numOfUPatch,
														const rg_INT& numOfVPatch,
														rg_REAL** paramValuesOfU,
														rg_REAL** paramValuesOfV,
														rg_Matrix** powerCoeff )
{
	rg_BzSurface3D** surfacePatch = new rg_BzSurface3D* [numOfUPatch];

	rg_INT i = 0;
	for( i=0; i < numOfUPatch; i++ )
		surfacePatch[i] = new rg_BzSurface3D[numOfVPatch];

	for( i=0; i < numOfUPatch; i++ )
	{
		rg_REAL intervalOfU[2];

		intervalOfU[0] = paramValuesOfU[i][0];
		intervalOfU[1] = paramValuesOfU[i][1];

		for( rg_INDEX j=0; j < numOfVPatch; j++ )
		{
			rg_REAL intervalOfV[2];

			intervalOfV[0] = paramValuesOfV[j][0];
			intervalOfV[1] = paramValuesOfV[j][1];

			surfacePatch[i][j].powerToBezierSurface( uDegree, vDegree,
													 intervalOfU, intervalOfV,
													 powerCoeff[i][j] );
		}
	}

	bzSurfacesToNUBSplineSurface( numOfUPatch, numOfVPatch, surfacePatch );

	for( i=0; i < numOfUPatch; i++ )
		delete[] surfacePatch[i];

	delete[] surfacePatch;
}

void rg_NUBSplineSurface3D::bzSurfacesToNUBSplineSurface( const rg_INT& numOfUPatch,
													   const rg_INT& numOfVPatch,
													   rg_BzSurface3D** bzSurfaces )
{
	if( numOfUPatch <= 0 || numOfVPatch <= 0 )
		return;

	if( numOfUPatch > 1 && numOfVPatch > 1 )
	{
		for( rg_INDEX i=0; i < (numOfUPatch - 1); i++ )
		{
			for( rg_INDEX j=0; j < (numOfVPatch - 1); j++ )
			{
				if( bzSurfaces[i][j].getDegreeOfU() 
					!= bzSurfaces[i][j+1].getDegreeOfU()
					|| bzSurfaces[i][j].getDegreeOfV() 
					!= bzSurfaces[i][j+1].getDegreeOfV())
					return;
			}
		}
	}
	
	// set order & control points.
	rg_ORDER uOrder = bzSurfaces[0][0].getDegreeOfU() + 1;	
	rg_ORDER vOrder = bzSurfaces[0][0].getDegreeOfV() + 1;	

	rg_INT numOfRowCtrlPt = uOrder * numOfUPatch;
	rg_INT numOfColCtrlPt = vOrder * numOfVPatch;

	rg_BSplineSurface3D::setOrderOfU( uOrder );
	rg_BSplineSurface3D::setOrderOfV( vOrder );

	rg_BSplineSurface3D::setControlNet( numOfRowCtrlPt, numOfColCtrlPt );
	
	for( rg_INDEX i=0; i < numOfUPatch; i++ )
	{
		for( rg_INDEX u=0; u < uOrder; u++ )
		{
			for( rg_INDEX j=0; j < numOfVPatch; j++ )
			{
				for( rg_INDEX v=0; v < vOrder; v++ )
					rg_BSplineSurface3D::setPointOnControlNet( (i*uOrder)+u,
					                                        (j*vOrder)+v,
														    bzSurfaces[i][j].getCtrlPt( u, v ) );
			}
		}
	}

	// initialize knot.
	setInitialKnotVectorOfU();
	setInitialKnotVectorOfV();

	// set remaining knots.
	if( numOfUPatch > 1 )
	{
		for( rg_INDEX i=0; i < (numOfUPatch-1); i++ )
		{
			rg_REAL newKnot = (rg_REAL)(i+1)/(rg_REAL)numOfUPatch;
		
			for( rg_INDEX u=0; u < uOrder; u++ )
				setKnotValueOfU( (i+1)*uOrder+u, newKnot );
		}
	}

	if( numOfVPatch > 1 )
	{
		for( rg_INDEX j=0; j < (numOfVPatch-1); j++ )
		{
			rg_REAL newKnot = (rg_REAL)(j+1)/(rg_REAL)numOfVPatch;
		
			for( rg_INDEX v=0; v < vOrder; v++ )
				setKnotValueOfV( (j+1)*vOrder+v, newKnot );
		}
	}
}

////	Reparameterization
////	This function reparameterize knots from 0 to 1.
void rg_NUBSplineSurface3D::reparameterizationKnotVector()
{
    rg_INT nU     = rg_BSplineSurface3D::getRowOfControlNet();
    rg_INT uOrder = rg_BSplineSurface3D::getOrderOfU();
	rg_INT nV     = rg_BSplineSurface3D::getColumnOfControlNet();
    rg_INT vOrder = rg_BSplineSurface3D::getOrderOfV();
	
	rg_INT numOfUKnots = nU + uOrder;
	rg_INT numOfVKnots = nV + vOrder;

	if( rg_NE( u_knot_vector[0], 0. ) || rg_NE( u_knot_vector[numOfUKnots-1], 1. ) )
	{
		rg_REAL firstKnot = u_knot_vector[0];

		rg_INT i = 0;
		for( i=0; i < numOfUKnots; i++ )
			u_knot_vector[i] = u_knot_vector[i] - firstKnot;
	
		rg_REAL interval = u_knot_vector[numOfUKnots - 1] - u_knot_vector[0];

		for( i=0; i < numOfUKnots; i++ )
			u_knot_vector[i] = u_knot_vector[i] / interval;
	}

	if( rg_NE( v_knot_vector[0], 0. ) || rg_NE( v_knot_vector[numOfVKnots-1], 1. ) )
	{
		rg_REAL firstKnot = v_knot_vector[0];

		rg_INT i = 0;
		for( i=0; i < numOfVKnots; i++ )
			v_knot_vector[i] = v_knot_vector[i] - firstKnot;
	
		rg_REAL interval = v_knot_vector[numOfVKnots - 1] - v_knot_vector[0];

		for( i=0; i < numOfVKnots; i++ )
			v_knot_vector[i] = v_knot_vector[i] / interval;
	}
}




