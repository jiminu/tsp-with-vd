#include "SquaredDistFuncLineSegPoint.h"
#include "float.h"

SquaredDistFuncLineSegPoint::SquaredDistFuncLineSegPoint()
: DistFunc()
{
}

SquaredDistFuncLineSegPoint::SquaredDistFuncLineSegPoint(const SquaredDistFuncLineSegPoint& sDistFunc)
: DistFunc(sDistFunc)
{
	m_distFunc = sDistFunc.m_distFunc;
}

SquaredDistFuncLineSegPoint::SquaredDistFuncLineSegPoint(const LineSegment3D& linesegment, 
														 const rg_Point3D& point, 
														 rg_REAL lower /* = 0.0 */, 
														 rg_REAL upper /* = 1.0 */)
														 : DistFunc(LINESEGMENT_POINT, lower, upper)
{
	rg_Point3D sPt = linesegment.getStartPt();
	rg_Point3D ePt = linesegment.getEndPt();

	rg_REAL* coeff = new rg_REAL[ 3 ];
	coeff[ 0 ] = (sPt - point).squaredMagnitude();
	coeff[ 1 ] = 2.0 * (sPt - point).innerProduct(ePt - sPt);
	coeff[ 2 ] = (ePt - sPt).squaredMagnitude();

	for(rg_INT i = 0;i < NUM_COEFF_SDIST_FUNC_BTW_LINESEG_POINT;i++)
	{
		if(rg_ZERO(coeff[ i ], resNeg15))
		{
			// Force the coefficient to be equal to zero
			coeff[ i ] = 0.0;
		}
	}

	m_distFunc.setCoefficient(coeff);

	delete[] coeff;	
}

SquaredDistFuncLineSegPoint::~SquaredDistFuncLineSegPoint()
{
}

rg_INT SquaredDistFuncLineSegPoint::computeLocalMinimaInsideInterval(rg_dList<rg_REAL>& localMinima)
{
// 	rg_INT numOfLocalMinima = computeLocalExremaWithoutConsideringBoundary();
// 	for(rg_INT i = 0;i < numOfLocalMinima;i++)
// 	{
// 		if(m_statusOfSortedLocalExtrema[ i ] == LOCAL_MIN)
// 			localMinima.add(m_sortedLocalExtrema[ i ]);
// 	}
// 	return localMinima.getSize();
	rg_INT numOfLocalMinima = computeLocalExrema();
	if(numOfLocalMinima == 3)
	{
		if(m_statusOfSortedLocalExtrema[ 1 ] == LOCAL_MIN)
			localMinima.add(m_sortedLocalExtrema[ 1 ]);
	}
	return localMinima.getSize();
}

rg_INT SquaredDistFuncLineSegPoint::computeLocalExremaWithoutConsideringBoundary()
{
	rg_REAL* coeff = m_distFunc.getCoefficient();
	if(rg_ZERO(coeff[ 2 ], resNeg15))
	{
		m_numOfLocalExtrema = 0;
		return m_numOfLocalExtrema;
	}

	rg_REAL localExtremum = - coeff[ 1 ] / (2.0 * coeff[ 2 ]);
	if(rg_BTOR(m_interval[ 0 ], localExtremum, m_interval[ 1 ]))
	{
		m_numOfLocalExtrema = 1;
		m_sortedLocalExtrema = new rg_REAL[m_numOfLocalExtrema];
		m_statusOfSortedLocalExtrema = new LocalPtStatus[m_numOfLocalExtrema];
		if(rg_EQ(localExtremum, m_interval[ 0 ]))
			localExtremum = m_interval[ 0 ];
		else if(rg_EQ(localExtremum, m_interval[ 1 ]))
			localExtremum = m_interval[ 1 ];
		else
			{} // do nothing
		m_sortedLocalExtrema[ 0 ] = localExtremum;
		if(rg_POS(coeff[ 2 ]))
		{
			m_statusOfSortedLocalExtrema[ 0 ] = LOCAL_MIN;
		}
		if(rg_NEG(coeff[ 2 ]))
		{
			m_statusOfSortedLocalExtrema[ 0 ] = LOCAL_MAX;
		}
	}
	else
	{
		m_numOfLocalExtrema = 0;
	}
	return m_numOfLocalExtrema;
}

rg_INT SquaredDistFuncLineSegPoint::computeLocalMinima(rg_dList<rg_REAL>& localMinima)
{
	rg_INT numOfLocalMinima = computeLocalExrema();
	for(rg_INT i = 0;i < numOfLocalMinima;i++)
	{
		if(m_statusOfSortedLocalExtrema[ i ] == LOCAL_MIN)
			localMinima.add(m_sortedLocalExtrema[ i ]);
	}
	return localMinima.getSize();
}

rg_INT SquaredDistFuncLineSegPoint::computeLocalExrema()
{
	rg_REAL* coeff = m_distFunc.getCoefficient();
	if(rg_ZERO(coeff[ 2 ], resNeg15))
	{
		m_numOfLocalExtrema = 0;
		return m_numOfLocalExtrema;
	}

	rg_REAL localExtremum = - coeff[ 1 ] / (2.0 * coeff[ 2 ]);
	if(rg_BTOR(m_interval[ 0 ], localExtremum ,m_interval[ 1 ]))
	{
		if(rg_EQ(localExtremum, m_interval[ 0 ]) || rg_EQ(localExtremum, m_interval[ 1 ]))
		{
			m_numOfLocalExtrema = 2;
			m_sortedLocalExtrema = new rg_REAL[m_numOfLocalExtrema];
			m_statusOfSortedLocalExtrema = new LocalPtStatus[m_numOfLocalExtrema];
			m_sortedLocalExtrema[ 0 ] = m_interval[ 0 ];
			m_sortedLocalExtrema[ 1 ] = m_interval[ 1 ];
			rg_REAL fVal0 = evaluate(m_sortedLocalExtrema[ 0 ]);
			rg_REAL fVal1 = evaluate(m_sortedLocalExtrema[ 1 ]);
			if(rg_GE(fVal0, fVal1))
			{
				m_statusOfSortedLocalExtrema[ 0 ] = LOCAL_MAX;
				m_statusOfSortedLocalExtrema[ 1 ] = LOCAL_MIN;
			}
			else
			{
				m_statusOfSortedLocalExtrema[ 0 ] = LOCAL_MIN;
				m_statusOfSortedLocalExtrema[ 1 ] = LOCAL_MAX;
			}
		}
		else
		{
			m_numOfLocalExtrema = 3;
			m_sortedLocalExtrema = new rg_REAL[m_numOfLocalExtrema];
			m_statusOfSortedLocalExtrema = new LocalPtStatus[m_numOfLocalExtrema];

			m_sortedLocalExtrema[ 0 ] = m_interval[ 0 ];
			m_sortedLocalExtrema[ 1 ] = localExtremum;
			m_sortedLocalExtrema[ 2 ] = m_interval[ 1 ];
			if(rg_POS(coeff[ 2 ]))
			{
				m_statusOfSortedLocalExtrema[ 0 ] = LOCAL_MAX;
				m_statusOfSortedLocalExtrema[ 1 ] = LOCAL_MIN;
				m_statusOfSortedLocalExtrema[ 2 ] = LOCAL_MAX;
			}
			if(rg_NEG(coeff[ 2 ]))
			{
				m_statusOfSortedLocalExtrema[ 0 ] = LOCAL_MIN;
				m_statusOfSortedLocalExtrema[ 1 ] = LOCAL_MAX;
				m_statusOfSortedLocalExtrema[ 2 ] = LOCAL_MIN;
			}
		}
	}
	else
	{
		m_numOfLocalExtrema = 2;
		m_sortedLocalExtrema = new rg_REAL[m_numOfLocalExtrema];
		m_statusOfSortedLocalExtrema = new LocalPtStatus[m_numOfLocalExtrema];

		m_sortedLocalExtrema[ 0 ] = m_interval[ 0 ];
		m_sortedLocalExtrema[ 1 ] = m_interval[ 1 ];
		rg_REAL fVal0 = evaluate(m_sortedLocalExtrema[ 0 ]);
		rg_REAL fVal1 = evaluate(m_sortedLocalExtrema[ 1 ]);
		if(rg_GE(fVal0, fVal1))
		{
			m_statusOfSortedLocalExtrema[ 0 ] = LOCAL_MAX;
			m_statusOfSortedLocalExtrema[ 1 ] = LOCAL_MIN;
		}
		else
		{
			m_statusOfSortedLocalExtrema[ 0 ] = LOCAL_MIN;
			m_statusOfSortedLocalExtrema[ 1 ] = LOCAL_MAX;
		}
	}
	return m_numOfLocalExtrema;
}


LocalPtStatus SquaredDistFuncLineSegPoint::getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point) const
{
// 	if(m_numOfLocalExtrema == 3)
// 	{
// 		return m_statusOfSortedLocalExtrema[ 1 ];
// 	}
// 	else
// 		return LOCAL_UNKNOWN;
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

// LocalPtStatus SquaredDistFuncLineSegPoint::getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point, rg_REAL& fVal) const
// {
// 	if(m_numOfLocalExtrema == 3)
// 	{
// 		fVal = SquaredDistFuncLineSegPoint::evaluate(m_sortedLocalExtrema[ 1 ]);
// 		return m_statusOfSortedLocalExtrema[ 1 ];
// 	}
// 	else
// 	{
// 		fVal = DBL_MAX;
// 		return LOCAL_UNKNOWN;
// 	}
// }

LocalPtStatus SquaredDistFuncLineSegPoint::getStatusOfNearestExtremum(const rg_REAL& point) const
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

rg_REAL SquaredDistFuncLineSegPoint::evaluate(const rg_REAL& point) const
{
	if(rg_BTOR(m_interval[ 0 ], point ,m_interval[ 1 ]))
		return m_distFunc.evaluatePolynomial(point);
	else
		return DBL_MAX;
}

// rg_INT SquaredDistFuncLineSegPoint::computeInflectionPoints(rg_dList<rg_REAL>& inflectionPts) const
// {
// 	return 0;
// }

SquaredDistFuncLineSegPoint& SquaredDistFuncLineSegPoint::operator =(const SquaredDistFuncLineSegPoint& distFunc)
{
	if(this == &distFunc)
		return *this;

	DistFunc::operator =(distFunc);
	m_distFunc = distFunc.m_distFunc;	
	return *this;
}


