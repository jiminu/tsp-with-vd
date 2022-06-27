#include "SquaredDistFuncLineSegArc.h"

SquaredDistFuncLineSegArc::SquaredDistFuncLineSegArc()
{
	m_numOfIntervals = -1;
	m_intervalPoint = rg_NULL;
	m_sDistFunc = rg_NULL;
}

SquaredDistFuncLineSegArc::SquaredDistFuncLineSegArc(const SquaredDistFuncLineSegArc& sDistFunc)
{
	setIntervalPts_sDistFuncs(sDistFunc.m_numOfIntervals, 
		                      sDistFunc.m_intervalPoint, 
							  sDistFunc.m_sDistFunc);
}

SquaredDistFuncLineSegArc::SquaredDistFuncLineSegArc(const LineSegment3D& lineSegment, 
													 const Arc3D& arc, 
													 const rg_REAL& lower /* = 0.0 */, 
													 const rg_REAL& upper /* = 1.0 */)
: DistFunc(LINESEGMENT_ARC, lower, upper)
{
	rg_Point3D startPtOfArc = arc.getStartPoint();
	rg_Point3D endPtOfArc   = arc.getEndPoint();

	// Associate each point of a line segment with 
	// (1) a point on an arc except two end points: (the region R_M)
	// (2) start point of an arc                  : (the region R_S)
	// (3) end   point of an arc                  : (the region R_E)

	TrisectorOfSpaceByArc trisector(arc);
	rg_dList<rg_REAL> params;
	rg_dList<EffectiveRegion> regions;
	m_numOfIntervals = trisector.identifyEffectiveRegionsOfLineSegment(lineSegment, params, regions);
	
	m_intervalPoint = new rg_REAL[m_numOfIntervals + 1];
	rg_INT index = 0;
	params.reset4Loop();
	while(params.setNext4Loop())
	{
		m_intervalPoint[index++] = params.getEntity();
	}
	m_sDistFunc = new DistFunc* [m_numOfIntervals];	
	
	// Define distance functions for each part of a line segment
	index = 0;
	regions.reset4Loop();
	while(regions.setNext4Loop())
	{
		EffectiveRegion currRegion = regions.getEntity();
		switch(currRegion)
		{
		case Region_M:
			{
				SquaredDistFuncLineSegCircle* sDistFuncLineSegCircle = 
					new SquaredDistFuncLineSegCircle(lineSegment, 
													 arc,
													 m_intervalPoint[index], 
													 m_intervalPoint[index + 1]);
				m_sDistFunc[index++] = sDistFuncLineSegCircle;
				break;
			}			
		case Region_S:
			{
				SquaredDistFuncLineSegPoint* qDistFunc = 
					new SquaredDistFuncLineSegPoint(lineSegment, 
													startPtOfArc,
				 									m_intervalPoint[index], 
													m_intervalPoint[index + 1]);
				m_sDistFunc[index++] = qDistFunc;
				break;
			}			
		case Region_E:
			{
				SquaredDistFuncLineSegPoint* qDistFunc 
					= new SquaredDistFuncLineSegPoint(lineSegment, 
													  endPtOfArc, 
				 									  m_intervalPoint[index], 
													  m_intervalPoint[index + 1]);
				m_sDistFunc[index++] = qDistFunc;
				break;
			}
        default:
            break;
		}
	}
}

SquaredDistFuncLineSegArc::~SquaredDistFuncLineSegArc()
{
	if(m_intervalPoint)
	{
		delete [] m_intervalPoint;		
	}
	if(m_sDistFunc)
	{
		for(rg_INT i = 0;i < m_numOfIntervals;i++)
			delete [] m_sDistFunc[ i ];
		delete [] m_sDistFunc;
	}
	m_numOfIntervals = 0;
}

void SquaredDistFuncLineSegArc::setIntervalPts_sDistFuncs(const rg_INT& numOfIntervals, 
														  const rg_REAL* intervalPoint, 
														  DistFunc** sDistFunc)
{
	if(m_intervalPoint)
		delete [] m_intervalPoint;

	if(m_sDistFunc)
	{
		for(rg_INT i = 0;i < m_numOfIntervals;i++)
			delete [] m_sDistFunc[ i ];
		delete [] m_sDistFunc;
	}

	rg_INT i = 0;
	m_numOfIntervals = numOfIntervals;
	m_intervalPoint = new rg_REAL[m_numOfIntervals + 1];
	for(i = 0;i <= m_numOfIntervals;i++)
		m_intervalPoint[ i ] = intervalPoint[ i ];

	m_sDistFunc = new DistFunc* [m_numOfIntervals];
	for(i = 0;i < m_numOfIntervals;i++)
	{
		InvolvedEntityType type = sDistFunc[ i ]->getType();
		if(type == LINESEGMENT_POINT)
			m_sDistFunc[ i ] = new SquaredDistFuncLineSegPoint((*(SquaredDistFuncLineSegPoint*)sDistFunc[ i ]));
		// (type == LINESEGMENT_CIRCLE)
		else
			m_sDistFunc[ i ] = new SquaredDistFuncLineSegCircle((*(SquaredDistFuncLineSegCircle*)sDistFunc[ i ]));
	}
}

rg_INT  SquaredDistFuncLineSegArc::computeLocalMinima(rg_dList<rg_REAL>& localMinima)
{
	// Check inside interval of each distance function
	rg_INT i = 0;
	for(i = 0;i < m_numOfIntervals;i++)
	{
		rg_dList<rg_REAL> tLocalMinima;
		InvolvedEntityType type = m_sDistFunc[ i ]->getType();
		if(type == LINESEGMENT_POINT)
			((SquaredDistFuncLineSegPoint*)m_sDistFunc[ i ])->computeLocalMinimaInsideInterval(tLocalMinima);
		else if(type == LINESEGMENT_CIRCLE)
			((SquaredDistFuncLineSegCircle*)m_sDistFunc[ i ])->computeLocalMinimaInsideInterval(tLocalMinima);

		localMinima.appendTail(tLocalMinima);
	}

	//rg_REAL delta_neighborhood = computeDeltaNeighborhood(localMinima);

	// Check the breakpoints including parameter = 0.0 and 1.0
	for(i = 0;i <= m_numOfIntervals;i++)
	{
		if(!isThisPtIn(m_intervalPoint[ i ], localMinima) &&
			//isThisPtLocalMinOrMax(m_intervalPoint[ i ], delta_neighborhood) == LOCAL_MIN)
			isThisIntervalPtLocalMinOrMax(m_intervalPoint[ i ]) == LOCAL_MIN)
		{
			localMinima.add(m_intervalPoint[ i ]);
		}
	}
	return localMinima.getSize();
}

//LocalPtStatus SquaredDistFuncLineSegArc::isThisPtLocalMinOrMax(const rg_REAL& point, const rg_REAL& delta_neighborhood) const
LocalPtStatus SquaredDistFuncLineSegArc::isThisIntervalPtLocalMinOrMax(const rg_REAL& point) const
{
	if(rg_LE(m_intervalPoint[ 1 ], point) && rg_LT(point, m_intervalPoint[m_numOfIntervals]))
	{
		return isThisIntervalPtLocalMinOrMax(point, rg_FALSE);
	}
	// point == m_intervalPoint[ 0 ] or m_intervalPoint[m_numOfIntervals + 1]
	else
	{
		return isThisIntervalPtLocalMinOrMax(point, rg_TRUE);
	}
	/*
	if(rg_LT(m_intervalPoint[ 0 ], point) && rg_LT(point, m_intervalPoint[m_numOfIntervals]))
	{
		rg_REAL fVal       = evaluate(point);
		rg_REAL fVal_minus = evaluate(point - delta_neighborhood);
		rg_REAL fVal_plus  = evaluate(point + delta_neighborhood);

		if(rg_LE(fVal, fVal_minus) && rg_LE(fVal, fVal_plus))
			return LOCAL_MIN;
		else if(rg_GE(fVal, fVal_minus) && rg_GE(fVal, fVal_plus))
			return LOCAL_MAX;
		else
			return LOCAL_UNKNOWN;
	}
	*/
	// point == m_intervalPoint[ 0 ] or m_intervalPoint[m_numOfIntervals]
// 	else
// 	{
// 		return isThisPtLocalMinOrMax(point);
// 	}
	/*
	else if(rg_EQ(point, m_intervalPoint[ 0 ]))
	{
		rg_REAL fVal       = evaluate(point);
		rg_REAL fVal_plus  = evaluate(point + delta_neighborhood);

		if(rg_LE(fVal, fVal_plus))
			return LOCAL_MIN;
		else if(rg_GE(fVal, fVal_plus))
			return LOCAL_MAX;
		else
			return LOCAL_UNKNOWN;
	}
	else if(rg_EQ(point, m_intervalPoint[m_numOfIntervals]))
	{
		rg_REAL fVal       = evaluate(point);
		rg_REAL fVal_minus = evaluate(point - delta_neighborhood);

		if(rg_LE(fVal, fVal_minus))
			return LOCAL_MIN;
		else if(rg_GE(fVal, fVal_minus))
			return LOCAL_MAX;
		else
			return LOCAL_UNKNOWN;
	}
	else
	{
		// do nothing
	}
	*/
}

rg_REAL SquaredDistFuncLineSegArc::computeDeltaNeighborhood(rg_dList<rg_REAL>& localMinima) const
{
	// merge local minima with break points
	rg_INT i = 0;
	rg_dList<rg_REAL> ptList;
	for(i = 0;i < m_numOfIntervals;i++)
	{
		ptList.add(m_intervalPoint[ i ]);
		localMinima.reset4Loop();
		while(localMinima.setNext4Loop())
		{
			rg_REAL localMinimum = localMinima.getEntity();
			if(rg_LT(m_intervalPoint[ i ], localMinimum) && rg_LT(localMinimum, m_intervalPoint[i + 1]))
				ptList.add(localMinimum);
			else
				break;
		}
	}
	ptList.add(m_intervalPoint[m_numOfIntervals]);
	
	// compute minimum delta
	rg_REAL delta_neighborhood = DBL_MAX;
	rg_INT size = ptList.getSize();
	rg_dNode<rg_REAL>* currNode = ptList.getFirstpNode();
	for(i = 0;i < size - 1;currNode = currNode->getNext(),i++)
	{
		rg_REAL currDelta = rg_ABS(currNode->getEntity() - currNode->getNext()->getEntity());
		if(rg_LT(currDelta, delta_neighborhood))
			delta_neighborhood = currDelta;
	}

	if(delta_neighborhood != DBL_MAX)
		return delta_neighborhood / 2.0;
	else
		return 0.5;
}

LocalPtStatus SquaredDistFuncLineSegArc::isThisPtLocalMinOrMax(const rg_REAL& point) const
{
	// Get the status of nearest local extremum for a given point
	rg_INT indexBnd = m_numOfIntervals - 1;
	rg_INT i = 0;
	LocalPtStatus status;
	for(i = 0;i < indexBnd;i++)
	{
		if(rg_LE(m_intervalPoint[ i ], point) && rg_LT(point, m_intervalPoint[i + 1]))
		{
			InvolvedEntityType type = m_sDistFunc[ i ]->getType();
			if(type == LINESEGMENT_POINT)
				status = ((SquaredDistFuncLineSegPoint*)m_sDistFunc[ i ])->getStatusOfNearestExtremumWithoutConsideringBoundary(point);
			else if(type == LINESEGMENT_CIRCLE)
				status = ((SquaredDistFuncLineSegCircle*)m_sDistFunc[ i ])->getStatusOfNearestExtremumWithoutConsideringBoundary(point);
		}
	}

	if(rg_LE(m_intervalPoint[m_numOfIntervals - 1], point) && rg_LE(point, m_intervalPoint[m_numOfIntervals]))
	{
		InvolvedEntityType type = m_sDistFunc[ i ]->getType();
		if(type == LINESEGMENT_POINT)
			status = ((SquaredDistFuncLineSegPoint*)m_sDistFunc[ i ])->getStatusOfNearestExtremumWithoutConsideringBoundary(point);
		else if(type == LINESEGMENT_CIRCLE)
			status = ((SquaredDistFuncLineSegCircle*)m_sDistFunc[ i ])->getStatusOfNearestExtremumWithoutConsideringBoundary(point);
	}

	if(status == LOCAL_MIN)
		return LOCAL_MAX;
	else if(status == LOCAL_MAX)
		return LOCAL_MIN;
	//status == LOCAL_UNKNOWN
	else
	{
		if(rg_EQ(point, m_intervalPoint[ 0 ]) || rg_EQ(point, m_intervalPoint[m_numOfIntervals]))
		{
			rg_REAL fVal       = evaluate(point);
			rg_REAL fVal_neighborhood = evaluate((m_intervalPoint[ 0 ] + m_intervalPoint[m_numOfIntervals]) / 2.0);

			if(rg_LT(fVal, fVal_neighborhood))
				return LOCAL_MIN;
			else if(rg_GT(fVal, fVal_neighborhood))
				return LOCAL_MAX;
			else
				return LOCAL_UNKNOWN;
			}
		else
			return LOCAL_UNKNOWN;
	}
}

LocalPtStatus SquaredDistFuncLineSegArc::isThisIntervalPtLocalMinOrMax(const rg_REAL& point, const rg_FLAG& isStartOrEndIntervalPt) const
{
	if(isStartOrEndIntervalPt)
	{
		// Get the status of nearest local extremum for a given point
		rg_INDEX distFuncIndex = -1;
		rg_INDEX intervalIndex = -1;
		if(rg_EQ(point, m_intervalPoint[ 0 ]))
		{
			distFuncIndex = 0;
			intervalIndex = 0;
		}
		else
		{
			distFuncIndex = m_numOfIntervals - 1;
			intervalIndex = m_numOfIntervals;
		}

		LocalPtStatus status;
		InvolvedEntityType type = m_sDistFunc[distFuncIndex]->getType();
		if(type == LINESEGMENT_POINT)
			status = ((SquaredDistFuncLineSegPoint*)m_sDistFunc[distFuncIndex])->getStatusOfNearestExtremum(point);
		else if(type == LINESEGMENT_CIRCLE)
			status = ((SquaredDistFuncLineSegCircle*)m_sDistFunc[distFuncIndex])->getStatusOfNearestExtremum(point);

		if(status == LOCAL_MIN)
			return LOCAL_MAX;
		else if(status == LOCAL_MAX)
			return LOCAL_MIN;
		//status == LOCAL_UNKNOWN
		else
		{
			rg_REAL fVal       = evaluate(point);
			rg_REAL fVal_neighborhood;
			if(intervalIndex == 0)
				fVal_neighborhood = evaluate(m_intervalPoint[intervalIndex + 1]);
			// index == m_numOfIntervals
			else
				fVal_neighborhood = evaluate(m_intervalPoint[intervalIndex - 1]);

			if(rg_LT(fVal, fVal_neighborhood))
				return LOCAL_MIN;
			else if(rg_GT(fVal, fVal_neighborhood))
				return LOCAL_MAX;
			else
				return LOCAL_UNKNOWN;
		}
	}
	else
	{
		rg_INDEX intervalIndex = -1;		
		rg_INT indexBnd = m_numOfIntervals - 1;
		for(rg_INT i = 1;i <= indexBnd;i++)
		{
			if(rg_EQ(m_intervalPoint[ i ], point))
			{
				intervalIndex = i;
				break;
			}
		}
		rg_INDEX distFuncIndex = intervalIndex - 1;
		LocalPtStatus lStatus, rStatus;

		InvolvedEntityType type = m_sDistFunc[distFuncIndex]->getType();
		if(type == LINESEGMENT_POINT)
			lStatus = ((SquaredDistFuncLineSegPoint*)m_sDistFunc[distFuncIndex])->getStatusOfNearestExtremum(point);
		else if(type == LINESEGMENT_CIRCLE)
			lStatus = ((SquaredDistFuncLineSegCircle*)m_sDistFunc[distFuncIndex])->getStatusOfNearestExtremum(point);

		distFuncIndex = intervalIndex;
		type = m_sDistFunc[distFuncIndex]->getType();
		if(type == LINESEGMENT_POINT)
			rStatus = ((SquaredDistFuncLineSegPoint*)m_sDistFunc[distFuncIndex])->getStatusOfNearestExtremum(point);
		else if(type == LINESEGMENT_CIRCLE)
			rStatus = ((SquaredDistFuncLineSegCircle*)m_sDistFunc[distFuncIndex])->getStatusOfNearestExtremum(point);

		if(lStatus == LOCAL_MAX && rStatus == LOCAL_MAX)
			return LOCAL_MIN;
		else if(lStatus == LOCAL_MIN && rStatus == LOCAL_MIN)
			return LOCAL_MAX;
		else
			return LOCAL_UNKNOWN;
	}
}

rg_REAL SquaredDistFuncLineSegArc::evaluate(const rg_REAL& point) const
{
	rg_INT indexBnd = m_numOfIntervals - 1;
	rg_INT i = 0;
	for(i = 0;i < indexBnd;i++)
	{
		if(rg_LE(m_intervalPoint[ i ], point) && rg_LT(point, m_intervalPoint[i + 1]))
		{
			InvolvedEntityType type = m_sDistFunc[ i ]->getType();
			if(type == LINESEGMENT_POINT)
				return ((SquaredDistFuncLineSegPoint*)m_sDistFunc[ i ])->evaluate(point);
			else if(type == LINESEGMENT_CIRCLE)
				return ((SquaredDistFuncLineSegCircle*)m_sDistFunc[ i ])->evaluate(point);
		}
	}
	if(rg_LE(m_intervalPoint[m_numOfIntervals - 1], point) && rg_LE(point, m_intervalPoint[m_numOfIntervals]))
	{
		InvolvedEntityType type = m_sDistFunc[ i ]->getType();
		if(type == LINESEGMENT_POINT)
			return ((SquaredDistFuncLineSegPoint*)m_sDistFunc[ i ])->evaluate(point);
		else if(type == LINESEGMENT_CIRCLE)
			return ((SquaredDistFuncLineSegCircle*)m_sDistFunc[ i ])->evaluate(point);
	}
	return DBL_MAX;
}

rg_INT SquaredDistFuncLineSegArc::computeInflectionPoints(rg_dList<rg_REAL>& inflectionPts) const
{
	// Check inside interval of each distance function
	for(rg_INT i = 0;i < m_numOfIntervals;i++)
	{
		rg_dList<rg_REAL> tInflectionPts;
		InvolvedEntityType type = m_sDistFunc[ i ]->getType();
		if(type == LINESEGMENT_CIRCLE)
			((SquaredDistFuncLineSegCircle*)m_sDistFunc[ i ])->computeInflectionPoints(tInflectionPts);

		inflectionPts.appendTail(tInflectionPts);
	}
	return inflectionPts.getSize();	
}

rg_INT SquaredDistFuncLineSegArc::computeRootOfThirdDerivative(rg_dList<rg_REAL>& root) const
{
	// Check inside interval of each distance function
	for(rg_INT i = 0;i < m_numOfIntervals;i++)
	{
		rg_dList<rg_REAL> tRoot;
		InvolvedEntityType type = m_sDistFunc[ i ]->getType();
		if(type == LINESEGMENT_CIRCLE)
			((SquaredDistFuncLineSegCircle*)m_sDistFunc[ i ])->computeRootOfThirdDerivative(tRoot);

		root.appendTail(tRoot);
	}
	return root.getSize();
}

SquaredDistFuncLineSegArc& SquaredDistFuncLineSegArc::operator=(const SquaredDistFuncLineSegArc& sDistFunc)
{
	if(this == &sDistFunc)
		return *this;
	
	DistFunc::operator =(sDistFunc);

	setIntervalPts_sDistFuncs(sDistFunc.m_numOfIntervals, 
		                      sDistFunc.m_intervalPoint, 
							  sDistFunc.m_sDistFunc);
	
	return *this;
}


