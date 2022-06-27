#include "DistFunc.h"

rg_FLAG DistFunc::isThisPtIn(const rg_REAL& point, rg_dList<rg_REAL>& ptList) const
{
	ptList.reset4Loop();
	while(ptList.setNext4Loop())
	{
		if(rg_EQ(ptList.getEntity(), point))
			return rg_TRUE;
	}
	return rg_FALSE;
}

DistFunc::DistFunc()
{
	m_type = UNDEFINED;
	m_interval[ 0 ] = -DBL_MAX;
	m_interval[ 1 ] =  DBL_MAX;
	m_sortedLocalExtrema = rg_NULL;
	m_statusOfSortedLocalExtrema = rg_NULL;
	m_numOfLocalExtrema = 0;
}

DistFunc::DistFunc(const DistFunc& distFunc)
{
	m_type = distFunc.m_type;
	m_interval[ 0 ] = distFunc.m_interval[ 0 ];
	m_interval[ 1 ] = distFunc.m_interval[ 1 ];

	if(m_sortedLocalExtrema)
		delete [] m_sortedLocalExtrema;

	if(m_statusOfSortedLocalExtrema)
		delete [] m_statusOfSortedLocalExtrema;
	
	rg_INT i = 0;
	m_numOfLocalExtrema = distFunc.m_numOfLocalExtrema;
	m_sortedLocalExtrema = new rg_REAL[m_numOfLocalExtrema];
	for(i = 0;i < m_numOfLocalExtrema;i++)
		m_sortedLocalExtrema[ i ] = distFunc.m_sortedLocalExtrema[ i ];
	
	m_statusOfSortedLocalExtrema = new LocalPtStatus[m_numOfLocalExtrema];
	for(i = 0;i < m_numOfLocalExtrema;i++)
		m_statusOfSortedLocalExtrema[ i ] = distFunc.m_statusOfSortedLocalExtrema[ i ];
}

DistFunc::DistFunc(const InvolvedEntityType& type, 
				   const rg_REAL& lower, 
				   const rg_REAL& upper)
{
	m_type = type;
	m_interval[ 0 ] = lower;
	m_interval[ 1 ] = upper;

	m_sortedLocalExtrema = rg_NULL;
	m_statusOfSortedLocalExtrema = rg_NULL;
	m_numOfLocalExtrema = 0;
}

DistFunc::~DistFunc()
{
	m_type = UNDEFINED;
	m_interval[ 0 ] = -DBL_MAX;
	m_interval[ 1 ] =  DBL_MAX;

	if(m_sortedLocalExtrema)
	{
		delete [] m_sortedLocalExtrema;
		m_numOfLocalExtrema = 0;
	}
	if(m_statusOfSortedLocalExtrema)
		delete [] m_statusOfSortedLocalExtrema;
}

InvolvedEntityType DistFunc::getType() const
{
	return m_type;
}

void DistFunc::getInterval(rg_REAL interval[])
{
	interval[ 0 ] = m_interval[ 0 ];
	interval[ 1 ] = m_interval[ 1 ];
}

// rg_REAL DistFunc::evaluate(const rg_REAL& point) const
// {
// 	if(rg_BTOR(m_interval[ 0 ], point ,m_interval[ 1 ]))
// 		return -DBL_MAX;
// 	else
// 		return DBL_MAX;
// }
// 
// rg_REAL DistFunc::evaluateDerivative(const rg_REAL& point) const
// {
// 	if(rg_BTOR(m_interval[ 0 ], point ,m_interval[ 1 ]))
// 		return -DBL_MAX;
// 	else
// 		return DBL_MAX;
// }
// 
// rg_REAL DistFunc::evaluateDerivativeOfDerivative(const rg_REAL& point) const
// {
// 	if(rg_BTOR(m_interval[ 0 ], point ,m_interval[ 1 ]))
// 		return -DBL_MAX;
// 	else
// 		return DBL_MAX;
// }
// 
// rg_INT DistFunc::computeLocalMinimaWithoutConsideringBoundary(rg_dList<rg_REAL>& localMinima)
// {
// 	return localMinima.getSize();
// }
// 
// rg_INT DistFunc::computeLocalMinima(rg_dList<rg_REAL>& localMinima)
// {
// 	return localMinima.getSize();
// }
//
// LocalPtStatus DistFunc::getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point) const
// {
// 	return LOCAL_UNKNOWN;
// }
// 
// LocalPtStatus DistFunc::getStatusOfNearestExtremumWithoutConsideringBoundary(const rg_REAL& point, rg_REAL& fVal) const
// {
// 	fVal = DBL_MAX;
// 	return LOCAL_UNKNOWN;
// }
// LocalPtStatus DistFunc::getStatusOfNearestExtremum(const rg_REAL& point) const
// {
// 	return LOCAL_UNKNOWN;
// }

void DistFunc::setType(const InvolvedEntityType& type)
{
	m_type = type;
}

void DistFunc::setInterval(const rg_REAL& lower, const rg_REAL& upper)
{
	m_interval[ 0 ] = lower;
	m_interval[ 1 ] = upper;
}

DistFunc& DistFunc::operator =(const DistFunc& distFunc)
{
	if(this == &distFunc)
		return *this;

	m_type = distFunc.m_type;
	m_interval[ 0 ] = distFunc.m_interval[ 0 ];
	m_interval[ 1 ] = distFunc.m_interval[ 1 ];

	if(m_sortedLocalExtrema)
		delete [] m_sortedLocalExtrema;

	if(m_statusOfSortedLocalExtrema)
		delete [] m_statusOfSortedLocalExtrema;
	
	rg_INT i = 0;
	m_numOfLocalExtrema = distFunc.m_numOfLocalExtrema;
	m_sortedLocalExtrema = new rg_REAL[m_numOfLocalExtrema];
	for(i = 0;i < m_numOfLocalExtrema;i++)
		m_sortedLocalExtrema[ i ] = distFunc.m_sortedLocalExtrema[ i ];
	
	m_statusOfSortedLocalExtrema = new LocalPtStatus[m_numOfLocalExtrema];
	for(i = 0;i < m_numOfLocalExtrema;i++)
		m_statusOfSortedLocalExtrema[ i ] = distFunc.m_statusOfSortedLocalExtrema[ i ];	

	return *this;
}


