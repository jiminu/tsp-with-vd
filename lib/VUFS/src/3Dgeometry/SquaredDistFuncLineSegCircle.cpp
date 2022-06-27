#include "SquaredDistFuncLineSegCircle.h"
#include "rg_MathFunc.h"

SquaredDistFuncLineSegCircle::SquaredDistFuncLineSegCircle()
: DistFunc()
{
}

SquaredDistFuncLineSegCircle::SquaredDistFuncLineSegCircle(const SquaredDistFuncLineSegCircle& sDistFunc)
: DistFunc(sDistFunc)
{
	for(rg_INT i = 0;i < NUM_COEFF_SDIST_FUNC_BTW_LINESEG_CIRCLE;i++)
		m_coeff[ i ] = sDistFunc.m_coeff[ i ];

	m_radius = sDistFunc.m_radius;
}

SquaredDistFuncLineSegCircle::SquaredDistFuncLineSegCircle(const LineSegment3D& lineSegment, 
														   const Circle3D& circle, 
														   rg_REAL lower /* = 0.0 */, 
														   rg_REAL upper /* = 1.0 */)
														   : DistFunc(LINESEGMENT_CIRCLE, lower, upper)
{
 	rg_Point3D sPt = lineSegment.getStartPt();
 	rg_Point3D ePt = lineSegment.getEndPt();
 	rg_Point3D center = circle.getCenter();
 	rg_Point3D normalVec = circle.getNormal().getUnitVector();
 	rg_REAL    dist   = - normalVec.innerProduct(center);
 
 	rg_Point3D sPtMinusEPt = (sPt - ePt);
 	rg_REAL    sqdMgntudeOfSPtMinusEPt = sPtMinusEPt.squaredMagnitude();
 	rg_Point3D centerMinusSPt = (center - sPt);
 	rg_REAL    sqdMgntudeOfCenterMinusSPt = centerMinusSPt.squaredMagnitude();
 	rg_REAL    normalInnerProdSPt = normalVec.innerProduct(sPt);
 	rg_REAL    normalInnerProdSPtPlusDist = normalInnerProdSPt + dist;
 	rg_REAL    normalInnerProdSPtMinusEPt = normalVec.innerProduct(sPtMinusEPt);
 	           m_radius = circle.getRadius();
 	rg_REAL    sqdRadius = m_radius * m_radius;
 
 	m_coeff[ 6 ] = sqdMgntudeOfSPtMinusEPt;
 	m_coeff[ 5 ] = 2.0 * centerMinusSPt.innerProduct(sPtMinusEPt);
 	m_coeff[ 4 ] = sqdMgntudeOfCenterMinusSPt + sqdRadius;
	m_coeff[ 3 ] = -2.0 * m_radius;
 	m_coeff[ 2 ] = sqdMgntudeOfSPtMinusEPt - normalInnerProdSPtMinusEPt * normalInnerProdSPtMinusEPt;
 	m_coeff[ 1 ] = 2.0 * centerMinusSPt.innerProduct(sPtMinusEPt) + 2.0 * normalInnerProdSPtMinusEPt * normalInnerProdSPtPlusDist;
 	m_coeff[ 0 ] = sqdMgntudeOfCenterMinusSPt - normalInnerProdSPtPlusDist * normalInnerProdSPtPlusDist;

	for(rg_INT i = 0;i < NUM_COEFF_SDIST_FUNC_BTW_LINESEG_CIRCLE;i++)
	{
		if(rg_ZERO(m_coeff[ i ], resNeg15))
		{
			// Force the coefficient to be equal to zero
			m_coeff[ i ] = 0.0;
		}
	}
}

SquaredDistFuncLineSegCircle::~SquaredDistFuncLineSegCircle()
{
}

rg_INT SquaredDistFuncLineSegCircle::computeLocalMinimaInsideInterval(rg_dList<rg_REAL>& localMinima)
{
	// compute local minima and maxima
	rg_dList<rg_REAL> localExtrema;
	rg_dList<LocalPtStatus> localStatuses;
	computeLocalExrema(localExtrema, localStatuses);
	
	// sort local minima and maxima accroding to their parameter value
	sortLocalExtremaAccordingToParam(localExtrema, localStatuses);

	rg_INDEX indexBnd = m_numOfLocalExtrema - 1;
	// check out whether two interval points are local extremum or not.
	m_statusOfSortedLocalExtrema[ 0 ] = isThisPtLocalMinOrMax(m_sortedLocalExtrema[ 0 ], m_statusOfSortedLocalExtrema[ 0 ]);
	m_statusOfSortedLocalExtrema[indexBnd] = isThisPtLocalMinOrMax(m_sortedLocalExtrema[indexBnd], m_statusOfSortedLocalExtrema[indexBnd]);
	
	for(rg_INT i = 1;i < indexBnd;i++)
	{
		m_statusOfSortedLocalExtrema[ i ] = isThisPtLocalMinOrMax(m_sortedLocalExtrema[ i ], m_statusOfSortedLocalExtrema[ i ]);
		if(m_statusOfSortedLocalExtrema[ i ] == LOCAL_MIN)
		{
			if(!isThisPtIn(m_sortedLocalExtrema[ i ], localMinima))
				localMinima.add(m_sortedLocalExtrema[ i ]);
		}
	}
	return localMinima.getSize();
}

LocalPtStatus SquaredDistFuncLineSegCircle::getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point) const
{
	if(m_numOfLocalExtrema == 0)
		return LOCAL_UNKNOWN;
	rg_REAL minDist = DBL_MAX;
	rg_INDEX minIndex = -1;
	for(rg_INT i = 0;i < m_numOfLocalExtrema;i++)
	{
		rg_REAL currDist = rg_ABS(m_sortedLocalExtrema[ i ] - point);
		if(!rg_ZERO(currDist) && rg_LT(currDist, minDist))
		{
			minDist = currDist;
			minIndex = i;
		}
	}
	if((0 <= minIndex) && (minIndex <= m_numOfLocalExtrema - 1))
		return m_statusOfSortedLocalExtrema[minIndex];
	else
		return LOCAL_UNKNOWN;
}

// LocalPtStatus SquaredDistFuncLineSegCircle::getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point, rg_REAL& fVal) const
// {
// 	if(m_numOfLocalExtrema == 0)
// 	{
// 		fVal = DBL_MAX;
// 		return LOCAL_UNKNOWN;
// 	}
// 	rg_REAL minDist = DBL_MAX;
// 	rg_INDEX minIndex = -1;
// 	for(rg_INT i = 0;i < m_numOfLocalExtrema;i++)
// 	{
// 		rg_REAL currDist = rg_ABS(m_sortedLocalExtrema[ i ] - point);
// 		if(rg_LT(currDist, minDist))
// 		{
// 			minDist = currDist;
// 			minIndex = i;
// 		}
// 	}
// 	fVal = SquaredDistFuncLineSegCircle::evaluate(m_sortedLocalExtrema[minIndex]);
// 	return m_statusOfSortedLocalExtrema[minIndex];
// }

LocalPtStatus SquaredDistFuncLineSegCircle::getStatusOfNearestExtremum(const rg_REAL& point) const
{
	if(m_numOfLocalExtrema == 0)
		return LOCAL_UNKNOWN;
	rg_REAL minDist = DBL_MAX;
	rg_INDEX minIndex = -1;
	for(rg_INT i = 0;i < m_numOfLocalExtrema;i++)
	{
		rg_REAL currDist = rg_ABS(m_sortedLocalExtrema[ i ] - point);
		if(!rg_ZERO(currDist) && rg_LT(currDist, minDist))
		{
			minDist = currDist;
			minIndex = i;
		}
	}
	if((0 <= minIndex) && (minIndex <= m_numOfLocalExtrema - 1))
		return m_statusOfSortedLocalExtrema[minIndex];
	else
		return LOCAL_UNKNOWN;
}

rg_INT SquaredDistFuncLineSegCircle::computeLocalExremaInsideInterval(rg_dList<rg_REAL>& localExtrema, rg_dList<LocalPtStatus>& localStatuses) const
{
	// compute the roots for zeros of derivative of squared distance function
	rg_dList<rg_REAL> roots;
	computeCriticalPoints(roots);

	roots.reset4Loop();
	while(roots.setNext4Loop())
	{
		rg_REAL currRoot = roots.getEntity();
		if(!isThisPtIn(currRoot, localExtrema))
		{
			localExtrema.add(currRoot);
			localStatuses.add(LOCAL_EXTREMUM);
		}
	}
	return localExtrema.getSize();
}

rg_INT SquaredDistFuncLineSegCircle::computeLocalMinima(rg_dList<rg_REAL>& localMinima)
{
	// compute local minima and maxima and
	// zero and one parameters are also included
	rg_dList<rg_REAL> localExtrema;
	rg_dList<LocalPtStatus> localStatuses;
	computeLocalExrema(localExtrema, localStatuses);
	
	// sort local minima and maxima accroding to their parameter value
	sortLocalExtremaAccordingToParam(localExtrema, localStatuses);

	for(rg_INT i = 0;i < m_numOfLocalExtrema;i++)
	{
		m_statusOfSortedLocalExtrema[ i ] = isThisPtLocalMinOrMax(m_sortedLocalExtrema[ i ], m_statusOfSortedLocalExtrema[ i ]);
		if(m_statusOfSortedLocalExtrema[ i ] == LOCAL_MIN)
		{
			if(!isThisPtIn(m_sortedLocalExtrema[ i ], localMinima))
				localMinima.add(m_sortedLocalExtrema[ i ]);
		}
	}
	return localMinima.getSize();
}

rg_INT SquaredDistFuncLineSegCircle::computeLocalExrema(rg_dList<rg_REAL>& localExtrema, rg_dList<LocalPtStatus>& localStatuses) const
{
	// compute the roots for zeros of derivative of squared distance function
	rg_dList<rg_REAL> roots;
	computeCriticalPoints(roots);

	rg_FLAG isMinParamIn = rg_FALSE;
	rg_FLAG isMaxParamIn  = rg_FALSE;
	roots.reset4Loop();
	while(roots.setNext4Loop())
	{
		rg_REAL currRoot = roots.getEntity();
		if(!isThisPtIn(currRoot, localExtrema))
		{
			if(rg_EQ(currRoot, m_interval[ 0 ]))
			{
				currRoot = m_interval[ 0 ];
				isMinParamIn = rg_TRUE;
			}
			else if(rg_EQ(currRoot, m_interval[ 1 ]))
			{
				currRoot = m_interval[ 1 ];
				isMaxParamIn = rg_TRUE;
			}
			else{}// do nothing
			localExtrema.add(currRoot);
			localStatuses.add(LOCAL_EXTREMUM);
		}
	}

	if(!isMinParamIn)
	{
		localExtrema.add(m_interval[ 0 ]);
		localStatuses.add(LOCAL_UNKNOWN);
	}
	if(!isMaxParamIn)
	{
		localExtrema.add(m_interval[ 1 ]);
		localStatuses.add(LOCAL_UNKNOWN);
	}

	return localExtrema.getSize();
}

void SquaredDistFuncLineSegCircle::computeCriticalPoints(rg_dList<rg_REAL>& roots) const
{
	rg_REAL coeffOfQP[ 5 ] = {0.0, 0.0, 0.0, 0.0, 0.0};

	// Since squared distance function has no variable(parameter) inside square root function,
	// the roots of derivative may be obtained by differentiating the quadratic equation and solving
	// the corresponding linear equation.
	if(rg_ZERO(m_coeff[ 2 ], resNeg15) && rg_ZERO(m_coeff[ 1 ], resNeg15))
	{
		coeffOfQP[ 0 ] = m_coeff[ 5 ];
		coeffOfQP[ 1 ] = 2.0 * m_coeff[ 6 ];
		if(!rg_ZERO(coeffOfQP[ 1 ], resNeg15))
		{
			rg_REAL root = -coeffOfQP[ 0 ] / coeffOfQP[ 1 ];
			if(rg_BTOR(m_interval[ 0 ], root, m_interval[ 1 ]))
			{
				if(rg_EQ(root, m_interval[ 0 ]))
					root = m_interval[ 0 ];
				if(rg_EQ(root, m_interval[ 1 ]))
					root = m_interval[ 1 ];
				roots.add(root);
			}
		}
	}
	// We square the derivative and solve the resulting quartic polynomial equation
	// NOTE: The degree of the squared derivative could be reduced depending on the situations.
	else
	{
		rg_REAL sqdCoeffOfDistFunc5 = m_coeff[ 6 ] * m_coeff[ 6 ];
		rg_REAL sqdCoeffOfDistFunc4 = m_coeff[ 5 ] * m_coeff[ 5 ];
		rg_REAL sqdRadius = m_radius * m_radius;
		
		coeffOfQP[ 4 ] = 4.0 * sqdCoeffOfDistFunc5 * m_coeff[ 2 ];
		coeffOfQP[ 3 ] = 4.0 * m_coeff[ 6 ] * (m_coeff[ 6 ] * m_coeff[ 1 ] + m_coeff[ 5 ] * m_coeff[ 2 ]);
		coeffOfQP[ 2 ] = 4.0 * m_coeff[ 6 ] * (m_coeff[ 6 ] * m_coeff[ 0 ] + m_coeff[ 5 ] * m_coeff[ 1 ]) +
			m_coeff[ 2 ] * (sqdCoeffOfDistFunc4 - 4.0 * sqdRadius * m_coeff[ 2 ]);
		coeffOfQP[ 1 ] = 4.0 * m_coeff[ 6 ] * m_coeff[ 5 ] * m_coeff[ 0 ] +
			m_coeff[ 1 ] * (sqdCoeffOfDistFunc4 - 4.0 * sqdRadius * m_coeff[ 2 ]);
		coeffOfQP[ 0 ] = sqdCoeffOfDistFunc4 * m_coeff[ 0 ] - sqdRadius * m_coeff[ 1 ] * m_coeff[ 1 ];
		
		rg_Polynomial derivative( 4 );
		derivative.setCoefficient(coeffOfQP);
		derivative.reconstruct();

		rg_DEGREE dg = derivative.getDegree();
		rg_ComplexNumber *root_param = rg_NULL;
		switch (dg)
		{
 		case 2:
			{
				rg_REAL *coeff2 = derivative.getCoefficient();
				rg_QuadraticPolynomial qDerivative(coeff2); 
				root_param = qDerivative.solve();
			}
			break;
		case 3:
			{
				rg_REAL *coeff3 = derivative.getCoefficient();
				rg_CubicPolynomial cDerivative(coeff3);
				root_param = cDerivative.solve();
			}
			break;
		case 4:
			{
				rg_REAL *coeff4 = derivative.getCoefficient();
				rg_QuarticPolynomial quarticDerivative(coeff4);
				root_param = derivative.solve();
			}
			break;
		}
		// Then, we filter out the roots that do not make the derivative to be zero.
		if(root_param)
		{
			for(rg_INT i = 0; i < dg; i++)
			{
				// 0 <= root_param <=1 && root_param is a real number
				if(root_param[i].isPureRealNumber() && rg_BTOR(m_interval[ 0 ], root_param[i].getRealNumber(), m_interval[ 1 ]))
				{
					rg_REAL r = root_param[i].getRealNumber();

					/*
					rg_FLAG isPositive1 = rg_POS( r * (m_coeff[ 2 ] * r + m_coeff[ 1 ]) + m_coeff[ 0 ] ) ;
					//rg_FLAG isPositive1 = rg_NNEG( r * (m_coeff[ 2 ] * r + m_coeff[ 1 ]) + m_coeff[ 0 ] ) ;
					rg_FLAG isPositiveOrZero2 = rg_FALSE;
					rg_REAL determinant1 = ( 2.0 * m_coeff[ 2 ] * r + m_coeff[ 1 ] );
					rg_REAL determinant2 = ( 2.0 * m_coeff[ 6 ] * r + m_coeff[ 5 ] );
					if( rg_POS(determinant1 * determinant2) ||
					   (rg_ZERO(determinant1) && rg_ZERO(determinant2)) )
					   isPositiveOrZero2 = rg_TRUE;

					if(isPositive1 && isPositiveOrZero2)
						roots.add( r );
					*/

					rg_REAL denom = r * (m_coeff[ 2 ] * r + m_coeff[ 1 ]) + m_coeff[ 0 ];
					// squared distance function is not differentiable at this point
					if(rg_ZERO(denom))
						roots.add( r );
					// differentiable
					else if(rg_POS(denom))
					{
						rg_REAL fVal = 2.0 * m_coeff[ 6 ] * r + m_coeff[ 5 ] 
							         - m_radius * (2.0 * m_coeff[ 2 ] * r + m_coeff[ 1 ]) / sqrt(denom);
						// derivative == 0
						if(rg_ZERO(fVal))
						{
							roots.add( r );
						}
					}
					// squared distance function is not defined
					else
					{// do nothing
					}
				}
			}
		}
	}
}

void SquaredDistFuncLineSegCircle::sortLocalExtremaAccordingToParam(rg_dList<rg_REAL>& unsortedLocalExtrema, rg_dList<LocalPtStatus>& unsortedLocalPtStatuses)
{
	// Copy elements into array
	m_numOfLocalExtrema = unsortedLocalExtrema.getSize();
	m_sortedLocalExtrema = new rg_REAL[m_numOfLocalExtrema];
	m_statusOfSortedLocalExtrema = new LocalPtStatus[m_numOfLocalExtrema];
	rg_INT index = 0;
	unsortedLocalExtrema.reset4Loop();
	unsortedLocalPtStatuses.reset4Loop();
	while (unsortedLocalExtrema.setNext4Loop() && unsortedLocalPtStatuses.setNext4Loop())
	{
		m_sortedLocalExtrema[index] = unsortedLocalExtrema.getEntity();
		m_statusOfSortedLocalExtrema[index++] = unsortedLocalPtStatuses.getEntity();
	}
	// Sort elements by selection sort
	for(rg_INT i = 0;i < m_numOfLocalExtrema;i++)
	{
		rg_INT minIndex = i;
		rg_INT indexBnd = i + 1;
		for(rg_INT j = indexBnd;j < m_numOfLocalExtrema;j++)
		{
			if(rg_LT(m_sortedLocalExtrema[ j ], m_sortedLocalExtrema[minIndex]))
				minIndex = j;
		}
		if(minIndex != i)
		{
			rg_REAL tLocalMinMax = m_sortedLocalExtrema[ i ];
			m_sortedLocalExtrema[ i ] = m_sortedLocalExtrema[minIndex];
			m_sortedLocalExtrema[minIndex] = tLocalMinMax;

			LocalPtStatus tLocalStatus = m_statusOfSortedLocalExtrema[ i ];
			m_statusOfSortedLocalExtrema[ i ] = m_statusOfSortedLocalExtrema[minIndex];
			m_statusOfSortedLocalExtrema[minIndex] = tLocalStatus;
		}
	}
}

void SquaredDistFuncLineSegCircle::computeSqdDerivativeForSqdDistFunc(rg_QuarticPolynomial& derivative) const
{
	rg_REAL coeffOfQP[ 5 ] = {0.0, 0.0, 0.0, 0.0, 0.0};

	if(rg_ZERO(m_coeff[ 2 ], resNeg15) && rg_ZERO(m_coeff[ 1 ], resNeg15))
	{		
		coeffOfQP[ 0 ] = m_coeff[ 5 ] * m_coeff[ 5 ];
		coeffOfQP[ 1 ] = 4.0 * m_coeff[ 5 ] * m_coeff[ 6 ];
		coeffOfQP[ 2 ] = 4.0 * m_coeff[ 6 ] * m_coeff[ 6 ];;
	}
	else
	{
		rg_REAL sqdCoeffOfDistFunc5 = m_coeff[ 6 ] * m_coeff[ 6 ];
		rg_REAL sqdCoeffOfDistFunc4 = m_coeff[ 5 ] * m_coeff[ 5 ];
		rg_REAL sqdRadius = m_radius * m_radius;
		
		coeffOfQP[ 4 ] = 4.0 * sqdCoeffOfDistFunc5 * m_coeff[ 2 ];
		coeffOfQP[ 3 ] = 4.0 * m_coeff[ 6 ] * (m_coeff[ 6 ] * m_coeff[ 1 ] + m_coeff[ 5 ] * m_coeff[ 2 ]);
		coeffOfQP[ 2 ] = 4.0 * m_coeff[ 6 ] * (m_coeff[ 6 ] * m_coeff[ 0 ] + m_coeff[ 5 ] * m_coeff[ 1 ]) +
			m_coeff[ 2 ] * (sqdCoeffOfDistFunc4 - 4.0 * sqdRadius * m_coeff[ 2 ]);
		coeffOfQP[ 1 ] = 4.0 * m_coeff[ 6 ] * m_coeff[ 5 ] * m_coeff[ 0 ] +
			m_coeff[ 1 ] * (sqdCoeffOfDistFunc4 - 4.0 * sqdRadius * m_coeff[ 2 ]);
		coeffOfQP[ 0 ] = sqdCoeffOfDistFunc4 * m_coeff[ 0 ] - sqdRadius * m_coeff[ 1 ] * m_coeff[ 1 ];		
	}
	derivative.setCoefficient(coeffOfQP);
}

void SquaredDistFuncLineSegCircle::computeDerivativeForSqdDistFunc(rg_Polynomial& derivative, rg_DEGREE& degree) const
{
	rg_REAL sqdCoeffOfDistFunc5 = m_coeff[ 6 ] * m_coeff[ 6 ];
	rg_REAL sqdCoeffOfDistFunc4 = m_coeff[ 5 ] * m_coeff[ 5 ];
	rg_REAL sqdRadius = m_radius * m_radius;
	
	rg_REAL coeffOfQP[ 5 ] = {0.0, 0.0, 0.0, 0.0, 0.0};
	
	coeffOfQP[ 4 ] = 4.0 * sqdCoeffOfDistFunc5 * m_coeff[ 2 ];
	coeffOfQP[ 3 ] = 4.0 * m_coeff[ 6 ] * (m_coeff[ 6 ] * m_coeff[ 1 ] + m_coeff[ 5 ] * m_coeff[ 2 ]);
	coeffOfQP[ 2 ] = 4.0 * m_coeff[ 6 ] * (m_coeff[ 6 ] * m_coeff[ 0 ] + m_coeff[ 5 ] * m_coeff[ 1 ]) +
		m_coeff[ 2 ] * (sqdCoeffOfDistFunc4 - 4.0 * sqdRadius * m_coeff[ 2 ]);
	coeffOfQP[ 1 ] = 4.0 * m_coeff[ 6 ] * m_coeff[ 5 ] * m_coeff[ 0 ] +
		m_coeff[ 1 ] * (sqdCoeffOfDistFunc4 - 4.0 * sqdRadius * m_coeff[ 2 ]);
 	coeffOfQP[ 0 ] = sqdCoeffOfDistFunc4 * m_coeff[ 0 ] - sqdRadius * m_coeff[ 1 ] * m_coeff[ 1 ];

	degree = 0;
	for(rg_INT i = 1;i <= 4;i++)
	{
		if(coeffOfQP[ i ] > 0)
			degree++;
		else
			break;
	}		
}

LocalPtStatus SquaredDistFuncLineSegCircle::isThisPtLocalMinOrMax(const rg_REAL& point, const LocalPtStatus& status) const
{	
	rg_REAL denom = point * (m_coeff[ 2 ] * point + m_coeff[ 1 ]) + m_coeff[ 0 ];

	// derivative == 0 && second derivative exists
	if(status == LOCAL_EXTREMUM && rg_POS(denom) )
	{
		rg_REAL numer = 2.0 * m_coeff[ 2 ] * point + m_coeff[ 1 ];
		rg_REAL eval = 2.0 * m_coeff[ 6 ] + 
					   0.5 * m_radius * numer * numer / (pow(denom, 1.5)) 
					 - 2.0 * m_radius * m_coeff[ 2 ] / sqrt(denom);
		if(rg_POS(eval))
			return LOCAL_MIN;
		else if (rg_NEG(eval))
			return LOCAL_MAX;
		else
			;// do nothing
	}

	// SquaredDistFuncLineSegArc::isThisIntervalPtLocalMinOrMax()와 동일한 방법으로 수정하는 것을 고려 할 것.
	// (가장 가까이 있는 extreumum의 status를 이용하는 방법)
	// 가장 가까이 있는 extreumum의 status가 결정되기 전 일 수도 있기 때문에, 현재는 주변(delta neighborhood)의 함수값을 이용하고 있음.

	// Otherwise, we need to look into the distance function value in delta neighborhood of a given point
	// Search interval which contains current point
	rg_INT indexOfCurrPt = -1;
	rg_INT indexBnd = m_numOfLocalExtrema - 1;
	for(rg_INT i = 0;i < indexBnd;i++)
	{
		if(rg_LE(m_sortedLocalExtrema[ i ], point) && rg_LT(point, m_sortedLocalExtrema[i + 1]))
		{
			indexOfCurrPt = i;
			break;
		}
	}
	// m_sortedLocalExtrema[ 0 ] == 0.0
	if(rg_EQ(point, m_sortedLocalExtrema[ 0 ]))
		indexOfCurrPt = 0;
	// m_sortedLocalExtrema[indexBnd] == 1.0
	else if(rg_EQ(point, m_sortedLocalExtrema[indexBnd]))
		indexOfCurrPt = indexBnd;
	else
		;// do nothing
	
	rg_REAL distToNearestLocalExtremum = -1.0;
	rg_INT  indexOfNearestLocalExtremum = -2;
	rg_REAL Delta_Neighborhood = -1.0;
	if(1 <= indexOfCurrPt && indexOfCurrPt <= indexBnd - 1)
	{
		// Find the distance to nearest local extremum and its index	
		rg_REAL distToLeftExtremum  = rg_ABS(point - m_sortedLocalExtrema[indexOfCurrPt]);
		rg_REAL distToRightExtremum = rg_ABS(m_sortedLocalExtrema[indexOfCurrPt + 1] - point);

		if(rg_LT(distToLeftExtremum, distToRightExtremum))
		{
			distToNearestLocalExtremum  = distToLeftExtremum;
			indexOfNearestLocalExtremum = indexOfCurrPt;
		}
		else
		{
			distToNearestLocalExtremum = distToRightExtremum;
			indexOfNearestLocalExtremum = indexOfCurrPt + 1;
		}
		// Find the delta neighborhood of a given point
		Delta_Neighborhood = distToNearestLocalExtremum / 2.0;
		if(rg_ZERO(Delta_Neighborhood))
		{
			rg_REAL leftDelta   = rg_ABS((point - m_sortedLocalExtrema[indexOfCurrPt - 1]) / 2.0);
			rg_REAL rightDelta  = rg_ABS((m_sortedLocalExtrema[indexOfCurrPt + 1] - point) / 2.0);
			Delta_Neighborhood = (leftDelta < rightDelta) ? leftDelta : rightDelta;		
		}
		else
			;// do nothing

		// Comparison of distance function value
		if(rg_LE(evaluate(point), evaluate(point + Delta_Neighborhood)) &&
		   rg_LE(evaluate(point), evaluate(point - Delta_Neighborhood))    )
		{
			return LOCAL_MIN;
		}
		else if(rg_GE(evaluate(point), evaluate(point + Delta_Neighborhood)) &&
			    rg_GE(evaluate(point), evaluate(point - Delta_Neighborhood))    )
		{
			return LOCAL_MAX;
		}
		else
			return LOCAL_NONE;
	}
	// indexOfCurrPt == 0 or indexOfCurrPt == indexBnd
	else
	{
		if(indexOfCurrPt == 0 )
		{
			// Find the distance to nearest local extremum and its index	
			rg_REAL distToLeftExtremum  = rg_ABS(point - m_sortedLocalExtrema[indexOfCurrPt]);
			rg_REAL distToRightExtremum = rg_ABS(m_sortedLocalExtrema[indexOfCurrPt + 1] - point);
			
			if(rg_LT(distToLeftExtremum, distToRightExtremum))
			{
				distToNearestLocalExtremum = distToLeftExtremum;
				indexOfNearestLocalExtremum = indexOfCurrPt;
			}
			else
			{
				distToNearestLocalExtremum = distToRightExtremum;
				indexOfNearestLocalExtremum = indexOfCurrPt + 1;
			}
			// Find the delta neighborhood of a given point
			Delta_Neighborhood = distToNearestLocalExtremum / 2.0;
			if(rg_ZERO(Delta_Neighborhood))
			{
				Delta_Neighborhood = rg_ABS((m_sortedLocalExtrema[ 1 ] - point) / 2.0);
				// Comparison of distance function value
				if(rg_LE(evaluate(point), evaluate(point + Delta_Neighborhood)))
					return LOCAL_MIN;
				else
					return LOCAL_MAX;
			}
			else
			{
				// Comparison of distance function value
				if(rg_LE(evaluate(point), evaluate(point + Delta_Neighborhood)) &&
				   rg_LE(evaluate(point), evaluate(point - Delta_Neighborhood))    )
				{
					return LOCAL_MIN;
				}
				else if(rg_GE(evaluate(point), evaluate(point + Delta_Neighborhood)) &&
					rg_GE(evaluate(point), evaluate(point - Delta_Neighborhood))    )
				{
					return LOCAL_MAX;
				}
				else
					return LOCAL_NONE;
			}
		}
		else if(indexOfCurrPt == indexBnd)
		{
			// Find the distance to nearest local extremum and its index	
			distToNearestLocalExtremum = 0.0;
			indexOfNearestLocalExtremum = indexOfCurrPt;

			// Find the delta neighborhood of a given point
			Delta_Neighborhood = rg_ABS((point - m_sortedLocalExtrema[indexOfCurrPt - 1]) / 2.0);
			// Comparison of distance function value
			if(rg_LE(evaluate(point), evaluate(point - Delta_Neighborhood)))
				return LOCAL_MIN;
			else
				return LOCAL_MAX;
		}
		else
			;// do nothing
	}

    return LOCAL_UNKNOWN;
}

LocalPtStatus SquaredDistFuncLineSegCircle::isThisPtLocalMinOrMax(const rg_REAL& point) const
{	
	rg_INT i = 0;
	for(i = 0;i < m_numOfLocalExtrema;i++)
	{
		if(rg_EQ(point, m_sortedLocalExtrema[ i ]))
		{
			return m_statusOfSortedLocalExtrema[ i ];
		}
	}
	LocalPtStatus status = LOCAL_UNKNOWN;

	rg_REAL denom = point * (m_coeff[ 2 ] * point + m_coeff[ 1 ]) + m_coeff[ 0 ];

	// derivative == 0 && second derivative exists
	if(status == LOCAL_EXTREMUM && rg_POS(denom) )
	{
		rg_REAL numer = 2.0 * m_coeff[ 2 ] * point + m_coeff[ 1 ];
		rg_REAL eval = 2.0 * m_coeff[ 6 ] + 
					   0.5 * m_radius * numer * numer / (pow(denom, 1.5)) 
					 - 2.0 * m_radius * m_coeff[ 2 ] / sqrt(denom);
		if(rg_POS(eval))
			return LOCAL_MIN;
		else if (rg_NEG(eval))
			return LOCAL_MAX;
	}

	// Otherwise, we need to look into the distance function value in delta neighborhood of a given point
	// Search interval which contains current point
	rg_INT indexOfCurrPt = -1;
	rg_INT indexBnd = m_numOfLocalExtrema - 1;
	for(i = 0;i < indexBnd;i++)
	{
		if(rg_LE(m_sortedLocalExtrema[ i ], point) && rg_LT(point, m_sortedLocalExtrema[i + 1]))
		{
			indexOfCurrPt = i;
			break;
		}
	}
	// m_sortedLocalExtrema[ 0 ] == 0.0
	if(rg_EQ(point, m_sortedLocalExtrema[ 0 ]))
		indexOfCurrPt = 0;
	// m_sortedLocalExtrema[indexBnd] == 1.0
	else if(rg_EQ(point, m_sortedLocalExtrema[indexBnd]))
		indexOfCurrPt = indexBnd;
	else
		; // do nothing
	
	rg_REAL distToNearestLocalExtremum = -1.0;
	rg_INT  indexOfNearestLocalExtremum = -2;
	rg_REAL Delta_Neighborhood = -1.0;
	if(1 <= indexOfCurrPt && indexOfCurrPt <= indexBnd - 1)
	{
		// Find the distance to nearest local extremum and its index	
		rg_REAL distToLeftExtremum  = rg_ABS(point - m_sortedLocalExtrema[indexOfCurrPt]);
		rg_REAL distToRightExtremum = rg_ABS(m_sortedLocalExtrema[indexOfCurrPt + 1] - point);

		if(rg_LT(distToLeftExtremum, distToRightExtremum))
		{
			distToNearestLocalExtremum  = distToLeftExtremum;
			indexOfNearestLocalExtremum = indexOfCurrPt;
		}
		else
		{
			distToNearestLocalExtremum = distToRightExtremum;
			indexOfNearestLocalExtremum = indexOfCurrPt + 1;
		}
		// Find the delta neighborhood of a given point
		Delta_Neighborhood = distToNearestLocalExtremum / 2.0;
		if(rg_ZERO(Delta_Neighborhood))
		{
			rg_REAL leftDelta   = rg_ABS((point - m_sortedLocalExtrema[indexOfCurrPt - 1]) / 2.0);
			rg_REAL rightDelta  = rg_ABS((m_sortedLocalExtrema[indexOfCurrPt + 1] - point) / 2.0);
			Delta_Neighborhood = (leftDelta < rightDelta) ? leftDelta : rightDelta;		
		}

		// Comparison of distance function value
		if(rg_LE(evaluate(point), evaluate(point + Delta_Neighborhood)) &&
		   rg_LE(evaluate(point), evaluate(point - Delta_Neighborhood))    )
		{
			return LOCAL_MIN;
		}
		else if(rg_GE(evaluate(point), evaluate(point + Delta_Neighborhood)) &&
			    rg_GE(evaluate(point), evaluate(point - Delta_Neighborhood))    )
		{
			return LOCAL_MAX;
		}
		else
			return LOCAL_NONE;
	}
	// indexOfCurrPt == 0 or indexOfCurrPt == indexBnd
	else
	{
		if(indexOfCurrPt == 0 )
		{
			// Find the distance to nearest local extremum and its index	
			rg_REAL distToLeftExtremum  = rg_ABS(point - m_sortedLocalExtrema[indexOfCurrPt]);
			rg_REAL distToRightExtremum = rg_ABS(m_sortedLocalExtrema[indexOfCurrPt + 1] - point);
			
			if(rg_LT(distToLeftExtremum, distToRightExtremum))
			{
				distToNearestLocalExtremum = distToLeftExtremum;
				indexOfNearestLocalExtremum = indexOfCurrPt;
			}
			else
			{
				distToNearestLocalExtremum = distToRightExtremum;
				indexOfNearestLocalExtremum = indexOfCurrPt + 1;
			}
			// Find the delta neighborhood of a given point
			Delta_Neighborhood = distToNearestLocalExtremum / 2.0;
			if(rg_ZERO(Delta_Neighborhood))
			{
				Delta_Neighborhood = rg_ABS((m_sortedLocalExtrema[ 1 ] - point) / 2.0);
				// Comparison of distance function value
				if(rg_LE(evaluate(point), evaluate(point + Delta_Neighborhood)))
					return LOCAL_MIN;
				else
					return LOCAL_MAX;
			}
			else
			{
				// Comparison of distance function value
				if(rg_LE(evaluate(point), evaluate(point + Delta_Neighborhood)) &&
				   rg_LE(evaluate(point), evaluate(point - Delta_Neighborhood))    )
				{
					return LOCAL_MIN;
				}
				else if(rg_GE(evaluate(point), evaluate(point + Delta_Neighborhood)) &&
					rg_GE(evaluate(point), evaluate(point - Delta_Neighborhood))    )
				{
					return LOCAL_MAX;
				}
				else
					return LOCAL_NONE;
			}
		}
		else if(indexOfCurrPt == indexBnd)
		{
			// Find the distance to nearest local extremum and its index	
			distToNearestLocalExtremum = 0.0;
			indexOfNearestLocalExtremum = indexOfCurrPt;

			// Find the delta neighborhood of a given point
			Delta_Neighborhood = rg_ABS((point - m_sortedLocalExtrema[indexOfCurrPt - 1]) / 2.0);
			// Comparison of distance function value
			if(rg_LE(evaluate(point), evaluate(point - Delta_Neighborhood)))
				return LOCAL_MIN;
			else
				return LOCAL_MAX;
		}
		else
			; // do nothing
	}
}

rg_REAL SquaredDistFuncLineSegCircle::evaluate(const rg_REAL& point) const
{
	if(rg_BTOR(m_interval[ 0 ], point ,m_interval[ 1 ]))
	{
		rg_REAL argOfSqrt = point * (m_coeff[ 2 ] * point + m_coeff[ 1 ]) + m_coeff[ 0 ];
		rg_REAL sqrtVal = 0.0;
		if(rg_POS(argOfSqrt))
			sqrtVal = sqrt(argOfSqrt);

		rg_REAL dist = point*(m_coeff[ 6 ] * point + m_coeff[ 5 ]) + m_coeff[ 4 ] 
					 + m_coeff[ 3 ] * sqrtVal;
		return dist;
	}
	else
		return DBL_MAX;
}

rg_INT SquaredDistFuncLineSegCircle::computeInflectionPoints(rg_dList<rg_REAL>& inflectionPts) const
{
	if(!rg_ZERO(m_coeff[ 6 ]))
	{
		// construct a cubic eq.
		rg_REAL constTerm =	m_radius * (4. * m_coeff[ 0 ] * m_coeff[ 2 ] - m_coeff[ 1 ] * m_coeff[ 1 ]);
		constTerm = constTerm * constTerm;
		rg_ComplexNumber cubicEq(constTerm / (16. * m_coeff[ 6 ] * m_coeff[ 6 ]), 0.0);
		rg_ComplexNumber* cRoot = rg_MathFunc::root(cubicEq, 3);
		for(rg_INT i = 0;i < 3;i++)
		{
			// construct a quadratic eq.
			rg_ComplexNumber* coeff2 = new rg_ComplexNumber[ 3 ];
			coeff2[ 0 ] = rg_ComplexNumber(m_coeff[ 0 ], 0.) - cRoot[ i ];
			coeff2[ 1 ] = rg_ComplexNumber(m_coeff[ 1 ], 0.);
			coeff2[ 2 ] = rg_ComplexNumber(m_coeff[ 2 ], 0.);
			rg_INT solNum = -1;
			rg_ComplexNumber* qRoot = rg_MathFunc::solveQuadraticEq(coeff2[ 2 ], coeff2[ 1 ], coeff2[ 0 ], solNum);
			for(rg_INT j = 0;j < solNum;j++)
			{
				if(qRoot[ j ].isPureRealNumber())
				{
					rg_REAL r = qRoot[ j ].getRealNumber();
					if(rg_BTOR(m_interval[ 0 ], r, m_interval[ 1 ]))
					{
						// check whether current root makes a second derivative of squared dist. func. to be zero
						rg_REAL denom = r * (r * m_coeff[ 2 ] + m_coeff[ 1 ]) + m_coeff[ 0 ];
						if(!rg_ZERO(denom))
						{
							rg_REAL numer = (2. * m_coeff[ 2 ] * r + m_coeff[ 1 ]);
							numer = numer * numer;
							rg_REAL val = 
								   2. * m_coeff[ 6 ] 
								- (2. * m_coeff[ 2 ] * m_radius) / sqrt(denom) 
								+ m_radius * numer / (2. * pow(denom, 1.5));
							if(rg_ZERO(val))
							{
								if(rg_EQ(r, m_interval[ 0 ]))
									r = m_interval[ 0 ];
								if(rg_EQ(r, m_interval[ 1 ]))
									r = m_interval[ 1 ];
								inflectionPts.add(r);
							}
						}
					}
				}
			}
			delete [] qRoot;
		}
		delete [] cRoot;
	}
	return inflectionPts.getSize();
}

rg_INT SquaredDistFuncLineSegCircle::computeRootOfThirdDerivative(rg_dList<rg_REAL>& root) const
{
	rg_REAL denom = 6.0 * m_coeff[ 2 ] * m_coeff[ 2 ] * m_coeff[ 0 ] * m_radius
		          - (3.0 / 2.0) * m_coeff[ 2 ] * m_coeff[ 1 ] * m_coeff[ 1 ] * m_radius;

	if(!rg_ZERO(denom))
	{
		rg_REAL numerator = (3.0 / 4.0) * m_coeff[ 1 ] * m_coeff[ 1 ] * m_coeff[ 1 ] * m_radius
			              - 3.0 * m_coeff[ 2 ] * m_coeff[ 1 ] * m_coeff[ 0 ] * m_radius;
		rg_REAL r = numerator / denom;
		if(rg_BTOR(m_interval[ 0 ], r, m_interval[ 1 ]))
			root.add( r );
	}
	return root.getSize();
}

SquaredDistFuncLineSegCircle& SquaredDistFuncLineSegCircle::operator =(const SquaredDistFuncLineSegCircle& sDistFunc)
{
	if(this == &sDistFunc)
		return *this;

	DistFunc::operator =(sDistFunc);

	for(rg_INT i = 0;i < NUM_COEFF_SDIST_FUNC_BTW_LINESEG_CIRCLE;i++)
		m_coeff[ i ] = sDistFunc.m_coeff[ i ];

	m_radius = sDistFunc.m_radius;
	
	return *this;
}

