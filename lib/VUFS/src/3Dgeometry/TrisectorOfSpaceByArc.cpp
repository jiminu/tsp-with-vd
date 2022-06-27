#include "TrisectorOfSpaceByArc.h"
#include "rg_TMatrix3D.h"


TrisectorOfSpaceByArc::TrisectorOfSpaceByArc()
{
}

TrisectorOfSpaceByArc::TrisectorOfSpaceByArc(const Arc3D& arc)
{
	m_arcType = arc.getArcType();
	rg_Point3D center    = arc.getCenter();
	rg_Point3D stPt      = arc.getStartPoint();
	rg_Point3D edPt      = arc.getEndPoint();
	rg_Point3D normalVec = arc.getNormal();
	rg_Point3D vecCenterToSt = stPt - center;
	rg_Point3D vecCenterToEd = edPt - center;

	rg_Point3D ptOnAxisOfArc = center + normalVec * arc.getRadius();
	rg_Point3D angularBisector = computeAngularBisector(vecCenterToSt, vecCenterToEd, normalVec);
	rg_Point3D ptOnBisector = center - arc.getRadius() * angularBisector;
	
	m_splittor[ 0 ].definePlaneByThreePoints(ptOnAxisOfArc, stPt, center);
	m_splittor[ 1 ].definePlaneByThreePoints(ptOnAxisOfArc, edPt, center);
	m_splittor[ 2 ].definePlaneByThreePoints(ptOnAxisOfArc, ptOnBisector, center);
}

rg_INT TrisectorOfSpaceByArc::identifyEffectiveRegionsOfLineSegment(const LineSegment3D& lineSeg,
													   		        rg_dList<rg_REAL>& params, 
															        rg_dList<EffectiveRegion>& regions)
{
	// check the effective regions of start and end points
	EffectiveRegion regionOfStPt, regionOfEdPt;
	NumberOfOccupiedRegion numOfOccupiedRegions = 
		getEffectiveRegion(lineSeg, regionOfStPt, regionOfEdPt);

	EffectiveRegion effectiveRegion[ 4 ];
	rg_INT          numOfEffectiveRegions = 0;
	rg_dList<rg_REAL> sortedParams;

	// Line segment occupies 1 effective region
	if(numOfOccupiedRegions == ONE)
	{
		effectiveRegion[ 0 ] = regionOfStPt;
		numOfEffectiveRegions = 1;
	}
	// Line segment occupies 1, 2 or 3 regions
	else
	{
		rg_INDEX sumOfSplittorIndices;
		rg_dList<rg_REAL> unsortedParams;
		EffectiveRegion firstIntersectedSplittor;
		// compute intersections with splittors
		computeIntersectionWithSplittors(lineSeg, 
			                             unsortedParams, 
								 		 sumOfSplittorIndices, 
										 firstIntersectedSplittor);
		// Sort parameters in ascending order
		sortIntersectionParams(unsortedParams, sortedParams);				
		rg_INT numOfIntersections = sortedParams.getSize();

		// 0 intersection
		// 1 effective region
		// This case can occur only when the arc is a major arc
		// OnSplittor_SM/Region_M or OnSplittor_ME/Region_M
		if(numOfIntersections == 0)
		{
			effectiveRegion[ 0 ] = Region_M;
			numOfEffectiveRegions = 1;
		}
		else
		{			
			locateEffectiveRegions(effectiveRegion, 
				                   numOfEffectiveRegions, 
								   regionOfStPt, 
								   regionOfEdPt, 
								   numOfIntersections, 
								   sumOfSplittorIndices, 
								   firstIntersectedSplittor);
		}
	}
	// Assign effective regions to appropriate intervals
	rg_INDEX regionIndex = 0;
	rg_REAL  tParam;
	params.add(0.0);
	sortedParams.reset4Loop();
	while(sortedParams.setNext4Loop())
	{
		tParam = sortedParams.getEntity();
		if(rg_LT(0.0, tParam) && rg_LT(tParam, 1.0))
		{
			params.add(tParam);
			regions.add(effectiveRegion[regionIndex]);
		}
		regionIndex++;
	}
	params.add(1.0);
	regions.add(effectiveRegion[regionIndex]);

	return regions.getSize();
}

NumberOfOccupiedRegion TrisectorOfSpaceByArc::getEffectiveRegion(const LineSegment3D& lineSeg, 
												                 EffectiveRegion& regionOfStPt, 
												                 EffectiveRegion& regionOfEdPt)
{
	rg_Point3D stPt = lineSeg.getStartPt();
	rg_Point3D edPt = lineSeg.getEndPt();
	regionOfStPt = getEffectiveRegion(stPt);
	regionOfEdPt = getEffectiveRegion(edPt);

	NumberOfOccupiedRegion numOfRegions = TWO_THREE_OR_TO_BE_VERIFIED;

	// Start and end points belong to the same region/splittor/axis
	if(regionOfStPt == regionOfEdPt)
	{
		if(m_arcType == MINOR_ARC)
			numOfRegions = ONE;
		// In this case, maximum three intersections with splittors are possible.
		// Hence, we need to check that possibility.
		else if(m_arcType == MAJOR_ARC && regionOfStPt == Region_M)
			numOfRegions = TWO_THREE_OR_TO_BE_VERIFIED;
		else
			numOfRegions = ONE;
	}
	// Start and end points belong to different region/splittor/axis
	else
	{
		if(regionOfStPt == OnAxisOfArc || regionOfEdPt == OnAxisOfArc)
		{
			numOfRegions = ONE;
		}		
		//else if(regionOfStPt != OnAxisOfArc && regionOfEdPt != OnAxisOfArc)
		else
		{
			if(m_arcType == MINOR_ARC)
			{
				// Two different splittors
				if( (regionOfStPt == OnSplittor_SM || regionOfStPt == OnSplittor_ME || regionOfStPt == OnSplittor_ES) &&
					(regionOfEdPt == OnSplittor_SM || regionOfEdPt == OnSplittor_ME || regionOfEdPt == OnSplittor_ES)    )
					numOfRegions = ONE;
				
				// One splittor and one region
				if(  regionOfStPt == OnSplittor_SM && (regionOfEdPt == Region_S || regionOfEdPt == Region_M) ||
					 regionOfStPt == OnSplittor_ME && (regionOfEdPt == Region_M || regionOfEdPt == Region_E) ||
					 regionOfStPt == OnSplittor_ES && (regionOfEdPt == Region_E || regionOfEdPt == Region_S)   )
					numOfRegions = ONE;

				// One splittor and one region
				if(  regionOfEdPt == OnSplittor_SM && (regionOfStPt == Region_S || regionOfStPt == Region_M) ||
					 regionOfEdPt == OnSplittor_ME && (regionOfStPt == Region_M || regionOfStPt == Region_E) ||
					 regionOfEdPt == OnSplittor_ES && (regionOfStPt == Region_E || regionOfStPt == Region_S)   )
					numOfRegions = ONE;
			}
			// m_arcType == MAJOR_ARC
			else
			{
				// Two different splittors
				if( ((regionOfStPt == OnSplittor_SM || regionOfStPt == OnSplittor_ES) &&
					 (regionOfEdPt == OnSplittor_SM || regionOfEdPt == OnSplittor_ES)   ) ||
					((regionOfStPt == OnSplittor_ME || regionOfStPt == OnSplittor_ES) &&
					 (regionOfEdPt == OnSplittor_ME || regionOfEdPt == OnSplittor_ES)   )   )
					numOfRegions = ONE;

				// One splittor and one region
				if(  regionOfStPt == OnSplittor_SM &&  regionOfEdPt == Region_S                              ||
					 regionOfStPt == OnSplittor_ME &&  regionOfEdPt == Region_E                              ||
					 regionOfStPt == OnSplittor_ES && (regionOfEdPt == Region_E || regionOfEdPt == Region_S)   )
					numOfRegions = ONE;

				// One splittor and one region
				if(  regionOfEdPt == OnSplittor_SM &&  regionOfStPt == Region_S                              ||
					 regionOfEdPt == OnSplittor_ME &&  regionOfStPt == Region_E                              ||
					 regionOfEdPt == OnSplittor_ES && (regionOfStPt == Region_E || regionOfStPt == Region_S)   )
					numOfRegions = ONE;
			}
		}
		// Otherwise
		if(numOfRegions != ONE)
			numOfRegions = TWO_THREE_OR_TO_BE_VERIFIED;
	}

	// If either of two extreme points belongs to either splittors or axis
	// locate their corresponding effective regions.
	if(numOfRegions == ONE)
		lcoateEffectiveRegion(regionOfStPt, regionOfEdPt);

	return numOfRegions;
}

void TrisectorOfSpaceByArc::lcoateEffectiveRegion(EffectiveRegion& regionOfStPt, 
										          EffectiveRegion& regionOfEdPt)
{
	if((regionOfStPt == Region_S || regionOfStPt == Region_M || regionOfStPt == Region_E) &&
	   (regionOfEdPt == Region_S || regionOfEdPt == Region_M || regionOfEdPt == Region_E)   )
	   return;
	else
	{
		if(regionOfStPt == OnSplittor_SM)
		{
			if(regionOfEdPt == OnSplittor_ES || regionOfEdPt == OnSplittor_SM || regionOfEdPt == OnAxisOfArc)
				regionOfStPt = regionOfEdPt = Region_S;
			else if(regionOfEdPt == OnSplittor_ME)
				regionOfStPt = regionOfEdPt = Region_M;
			// regionOfEdPt == Region_S or Region_M
			else
				regionOfStPt = regionOfEdPt;
		}
		else if(regionOfStPt == OnSplittor_ME)
		{
			if(regionOfEdPt == OnSplittor_ES || regionOfEdPt == OnSplittor_ME || regionOfEdPt == OnAxisOfArc)
				regionOfStPt = regionOfEdPt = Region_E;
			else if(regionOfEdPt == OnSplittor_SM)
				regionOfStPt = regionOfEdPt = Region_M;
			// regionOfEdPt == Region_M or Region_E
			else
				regionOfStPt = regionOfEdPt;
		}
		else if(regionOfStPt == OnSplittor_ES)
		{
			if(regionOfEdPt == OnSplittor_ME || regionOfEdPt == OnSplittor_ES || regionOfEdPt == OnAxisOfArc)
				regionOfStPt = regionOfEdPt = Region_E;
			else if(regionOfEdPt == OnSplittor_SM)
				regionOfStPt = regionOfEdPt = Region_S;
			// regionOfEdPt == Region_E or Region_S
			else
				regionOfStPt = regionOfEdPt;
		}
		// regionOfStPt == OnAxis
		else
		{
			if(regionOfEdPt == OnSplittor_ME || regionOfEdPt == OnSplittor_ES)
				regionOfStPt = regionOfEdPt = Region_E;
			else if(regionOfEdPt == OnSplittor_SM || regionOfEdPt == OnAxisOfArc)
				regionOfStPt = regionOfEdPt = Region_S;
			// regionOfEdPt == Region_S, Region_M or Region_E
			else
				regionOfStPt = regionOfEdPt;								
		}
	}
}

EffectiveRegion TrisectorOfSpaceByArc::getEffectiveRegion(const rg_Point3D& pt)
{
	rg_INT index[ 4 ] = {0, 1, 2, 0};
	rg_INT i = 0;
	for(i = 0;i < 3;i++)
	{
		if(m_arcType == MAJOR_ARC && i == 0)
		{
			if(!(m_splittor[index[ i ]].isOnOppositeNormalSideOrOnThisPlane(pt) && 
				 m_splittor[index[i+1]].isOnNormalSideOrOnThisPlane(pt)))
				 return REGION[ i ];
		}
		else if(m_splittor[index[ i ]].isOnNormalSide(pt) && 
		        m_splittor[index[i+1]].isOnOppositeNormalSide(pt))
		{
			return REGION[ i ];
		}
	}
	for(i = 0;i < 3;i++)
	{
		if(m_splittor[index[ i ]].isOnThisPlane(pt) && 
		   m_splittor[index[i+1]].isOnThisPlane(pt))
		{
			return OnAxisOfArc;
		}
	}
	for(i = 0;i < 3;i++)
	{
		if(m_splittor[ i ].isOnThisPlane(pt))
		{
			return REGION[3 + i];
		}
	}

    return RegionUnknown;
}

rg_INT TrisectorOfSpaceByArc::computeIntersectionWithSplittors(const LineSegment3D& lineSeg,
															   rg_dList<rg_REAL>& unsortedParams, 
															   rg_INDEX& sumOfSplittorIndices, 
															   EffectiveRegion& firstIntersectedSplittor)
{
	sumOfSplittorIndices = 0;
	rg_FLAG    isFirstIntersection = rg_TRUE;
	rg_Point3D tPt;
	rg_REAL    tParam;
	rg_INDEX index[ 5 ] = {0, 1, 2, 0, 1};

	for(rg_INT i = 0;i < 3;i++)
	{
		if(!m_splittor[index[ i ]].computeIntersectionWithLineSegment(lineSeg, tPt, tParam))
			continue;
			
		rg_FLAG isValidIntersection  =  rg_FALSE;
		if(m_arcType == MINOR_ARC)
		{
			if(m_splittor[index[i + 1]].isOnOppositeNormalSideOrOnThisPlane(tPt) && 
			   m_splittor[index[i + 2]].isOnNormalSideOrOnThisPlane(tPt))
				isValidIntersection = rg_TRUE;
		}
		//(m_arcType == MAJOR_ARC)
		else
		{
			if(i == 0)
			{
				if(m_splittor[index[i + 1]].isOnNormalSideOrOnThisPlane(tPt) && 
				   m_splittor[index[i + 2]].isOnNormalSideOrOnThisPlane(tPt))
					isValidIntersection = rg_TRUE;
			}
			else if(i == 1)
			{
				if(m_splittor[index[i + 1]].isOnOppositeNormalSideOrOnThisPlane(tPt) && 
				   m_splittor[index[i + 2]].isOnOppositeNormalSideOrOnThisPlane(tPt))
					isValidIntersection = rg_TRUE;
			}
			else
			{
				if(m_splittor[index[i + 1]].isOnOppositeNormalSideOrOnThisPlane(tPt) && 
				   m_splittor[index[i + 2]].isOnNormalSideOrOnThisPlane(tPt))
					isValidIntersection = rg_TRUE;
			}
		}

		if(isValidIntersection)
		{
			// Check if the intersection parameter is zero or one
			if(rg_ZERO(tParam) || rg_EQ(tParam, 1.0))
				continue;
			// ////////////////

			// Check if the intersection parameter is already found
			rg_FLAG isAlreadyFound = rg_FALSE;
			unsortedParams.reset4Loop();
			while(unsortedParams.setNext4Loop())
			{
				if(rg_EQ(unsortedParams.getEntity(), tParam))
				{
					isAlreadyFound = rg_TRUE;
					break;
				}
			}
			if(!isAlreadyFound)
			{
				unsortedParams.add(tParam);
				sumOfSplittorIndices += i;
			}
			if(isFirstIntersection)
			{
				firstIntersectedSplittor = REGION[3 + index[ i ]];
				isFirstIntersection = rg_FALSE;
			}
		}
	}
	return unsortedParams.getSize();
}

void TrisectorOfSpaceByArc::locateEffectiveRegions(EffectiveRegion effectiveRegion[], 
												   rg_INT& numOfEffectiveRegions, 
                                                   const EffectiveRegion& regionOfStPt,
									               const EffectiveRegion& regionOfEdPt,
												   const rg_INT& numOfIntersections, 
												   const rg_INDEX sumOfSplittorIndices, 
												   const EffectiveRegion& firstIntersectedSplittor)
{
	// 1 intersection
	// 0 < intersectionParam < 1
	// 2 effective regions
	if(numOfIntersections == 1)
	{
		if(firstIntersectedSplittor == OnSplittor_SM)
		{
			if(regionOfStPt == OnSplittor_ES || regionOfStPt == Region_S)
			{
				effectiveRegion[ 0 ] = Region_S;
				effectiveRegion[ 1 ] = Region_M;
			}
			// (regionOfStPt == OnSplittor_ME || regionOfStPt == Region_M)
			else
			{
				effectiveRegion[ 0 ] = Region_M;
				effectiveRegion[ 1 ] = Region_S;
			}
		}
		else if(firstIntersectedSplittor == OnSplittor_ME)
		{
			if(regionOfStPt == OnSplittor_SM || regionOfStPt == Region_M)
			{
				effectiveRegion[ 0 ] = Region_M;
				effectiveRegion[ 1 ] = Region_E;
			}
			// (regionOfStPt == OnSplittor_ES || regionOfStPt == Region_E)
			else
			{
				effectiveRegion[ 0 ] = Region_E;
				effectiveRegion[ 1 ] = Region_M;
			}
		}
		// firstIntersectedSplittor == OnSplittor_ES
		else
		{
			if(regionOfStPt == OnSplittor_ME || regionOfStPt == Region_E)
			{
				effectiveRegion[ 0 ] = Region_E;
				effectiveRegion[ 1 ] = Region_S;
			}
			// (regionOfStPt == OnSplittor_SM || regionOfStPt == Region_S)
			else
			{
				effectiveRegion[ 0 ] = Region_S;
				effectiveRegion[ 1 ] = Region_E;
			}
		}
		numOfEffectiveRegions = 2;
	}
	// 2 or 3 intersections
	else
	{
		// 2 intersections
		// 2 or 3 effective regions
		if(numOfIntersections == 2)
		{
			if(sumOfSplittorIndices == 1)
			{
				if(firstIntersectedSplittor == OnSplittor_SM)
				{					
					if(regionOfStPt == OnSplittor_SM)
					{
						effectiveRegion[ 0 ] = Region_M;
						effectiveRegion[ 1 ] = Region_E;
						numOfEffectiveRegions = 2;
					}
					// regionOfStPt == Region_S
					else
					{
						effectiveRegion[ 0 ] = Region_S;
						effectiveRegion[ 1 ] = Region_M;
						effectiveRegion[ 2 ] = Region_E;
						numOfEffectiveRegions = 3;
					}
				}
				//(firstIntersectedSplittor == OnSplittor_ME)
				else
				{
					if(regionOfStPt == OnSplittor_ME)
					{
						effectiveRegion[ 0 ] = Region_M;
						effectiveRegion[ 1 ] = Region_S;
						numOfEffectiveRegions = 2;
					}
					// regionOfStPt == Region_E
					else
					{
						effectiveRegion[ 0 ] = Region_E;
						effectiveRegion[ 1 ] = Region_M;
						effectiveRegion[ 2 ] = Region_S;
						numOfEffectiveRegions = 3;
					}
				}
			}
			else if(sumOfSplittorIndices == 2)
			{
				if(firstIntersectedSplittor == OnSplittor_ES)
				{
					if(regionOfStPt == OnSplittor_SM)
					{
						effectiveRegion[ 0 ] = Region_S;
						effectiveRegion[ 1 ] = Region_E;
						numOfEffectiveRegions = 2;
					}
					// regionOfStPt == Region_M
					else
					{
						effectiveRegion[ 0 ] = Region_M;
						effectiveRegion[ 1 ] = Region_S;
						effectiveRegion[ 2 ] = Region_E;
						numOfEffectiveRegions = 3;
					}
				}
				//(firstIntersectedSplittor == OnSplittor_SM)
				else
				{
					if(regionOfStPt == OnSplittor_ES)
					{
						effectiveRegion[ 0 ] = Region_S;
						effectiveRegion[ 1 ] = Region_M;
						numOfEffectiveRegions = 2;
					}
					// regionOfStPt == Region_E
					else
					{
						effectiveRegion[ 0 ] = Region_E;
						effectiveRegion[ 1 ] = Region_S;
						effectiveRegion[ 2 ] = Region_M;
						numOfEffectiveRegions = 3;
					}
				}
			}
			// sumOfSplittorIndices == 3
			else
			{
				if(firstIntersectedSplittor == OnSplittor_ME)
				{
					if(regionOfStPt == OnSplittor_ME)
					{
						effectiveRegion[ 0 ] = Region_E;
						effectiveRegion[ 1 ] = Region_S;
						numOfEffectiveRegions = 2;
					}
					// regionOfStPt == Region_M
					else
					{
						effectiveRegion[ 0 ] = Region_M;
						effectiveRegion[ 1 ] = Region_E;
						effectiveRegion[ 2 ] = Region_S;
						numOfEffectiveRegions = 3;
					}
				}
				//(firstIntersectedSplittor == OnSplittor_ES)
				else
				{
					if(regionOfStPt == OnSplittor_ES)
					{
						effectiveRegion[ 0 ] = Region_E;
						effectiveRegion[ 1 ] = Region_M;
						numOfEffectiveRegions = 2;
					}
					// regionOfStPt == Region_S
					else
					{
						effectiveRegion[ 0 ] = Region_S;
						effectiveRegion[ 1 ] = Region_E;
						effectiveRegion[ 2 ] = Region_M;
						numOfEffectiveRegions = 3;
					}
				}
			}
		}
		// 3 intersections
		// 2, 3 or 4 effective regions
		// (For a major arc only)
		else
		{
			if(firstIntersectedSplittor == OnSplittor_SM)
			{
				if(regionOfStPt == OnSplittor_SM && regionOfEdPt == OnSplittor_ME)
				{
					effectiveRegion[ 0 ] = Region_S;
					effectiveRegion[ 1 ] = Region_E;
					numOfEffectiveRegions = 2;
				}
				else if(regionOfStPt == OnSplittor_SM && regionOfEdPt == Region_M)
				{
					effectiveRegion[ 0 ] = Region_S;
					effectiveRegion[ 1 ] = Region_E;
					effectiveRegion[ 2 ] = Region_M;
					numOfEffectiveRegions = 3;
				}
				else if(regionOfStPt == Region_M && regionOfEdPt == OnSplittor_ME)
				{
					effectiveRegion[ 0 ] = Region_M;
					effectiveRegion[ 1 ] = Region_S;
					effectiveRegion[ 2 ] = Region_E;
					numOfEffectiveRegions = 3;
				}
				else
				{
					effectiveRegion[ 0 ] = Region_M;
					effectiveRegion[ 1 ] = Region_S;
					effectiveRegion[ 2 ] = Region_E;
					effectiveRegion[ 3 ] = Region_M;
					numOfEffectiveRegions = 4;
				}
			}
			else
			{
				if(regionOfStPt == OnSplittor_ME && regionOfEdPt == OnSplittor_SM)
				{
					effectiveRegion[ 0 ] = Region_E;
					effectiveRegion[ 1 ] = Region_S;
					numOfEffectiveRegions = 2;
				}
				else if(regionOfStPt == OnSplittor_ME && regionOfEdPt == Region_M)
				{
					effectiveRegion[ 0 ] = Region_E;
					effectiveRegion[ 1 ] = Region_S;
					effectiveRegion[ 2 ] = Region_M;
					numOfEffectiveRegions = 3;
				}
				else if(regionOfStPt == Region_M && regionOfEdPt == OnSplittor_SM)
				{
					effectiveRegion[ 0 ] = Region_M;
					effectiveRegion[ 1 ] = Region_E;
					effectiveRegion[ 2 ] = Region_S;
					numOfEffectiveRegions = 3;
				}
				else
				{
					effectiveRegion[ 0 ] = Region_M;
					effectiveRegion[ 1 ] = Region_E;
					effectiveRegion[ 2 ] = Region_S;
					effectiveRegion[ 3 ] = Region_M;
					numOfEffectiveRegions = 4;
				}
			}
		}			
	}
}

void TrisectorOfSpaceByArc::sortIntersectionParams(rg_dList<rg_REAL>& unsortedParams, rg_dList<rg_REAL>& sortedParams)
{
	rg_INT numIntersections = unsortedParams.getSize();
	rg_REAL* tParams = new rg_REAL[numIntersections];

	rg_INT i = 0;
	unsortedParams.reset4Loop();
	while(unsortedParams.setNext4Loop())
	{
		tParams[ i++ ] = unsortedParams.getEntity();
	}
	for(i = 0;i < numIntersections;i++)
	{
		rg_INT minIndex = i;
		for(rg_INT j = i+1;j < numIntersections;j++)
		{				
			if(rg_LT(tParams[ j ], tParams[minIndex]))
				minIndex = j;
		}
		if(minIndex != i)
		{
			rg_REAL temp = tParams[ i ];
			tParams[ i ] = tParams[minIndex];
			tParams[minIndex] = temp;
		}
	}
	
	for(i = 0;i < numIntersections;i++)
		sortedParams.add(tParams[ i ]);
	delete [] tParams;
}

rg_INT TrisectorOfSpaceByArc::identifyRelativePosOfLineSegment_old(const LineSegment3D& lineSeg,
															       rg_dList<rg_REAL>& params, 
															       rg_dList<EffectiveRegion>& regions)
{
	EffectiveRegion regionOfStPt, regionOfEdPt;
	NumberOfOccupiedRegion numOfOccupiedRegions = 
		getEffectiveRegion(lineSeg, regionOfStPt, regionOfEdPt);

	// Occupy 1 effective region
	if(numOfOccupiedRegions == ONE)
	{
		params.add(0.0);
		params.add(1.0);

		if(regionOfStPt == Region_M)
		{
			regions.add(Region_M);
		}
		else if(regionOfStPt == Region_E      || 
			    regionOfStPt == OnSplittor_ME || 
				regionOfStPt == OnSplittor_ES)
		{
			regions.add(Region_E);
		}
		else
		{
			regions.add(Region_S);
		}
	}
	// Occupy 1, 2 or 3 regions
	else
	{
		rg_INDEX sumOfSplittorIndices;
		rg_dList<rg_REAL> unsortedParams;
		EffectiveRegion firstIntersectedSplittor;
		// compute intersections with splittors
		computeIntersectionWithSplittors(lineSeg, 
			                             unsortedParams, 
										 sumOfSplittorIndices, 
										 firstIntersectedSplittor);
		// sort parameters in ascending order
		rg_dList<rg_REAL> sortedParams;
		sortIntersectionParams(unsortedParams, sortedParams);
				
		EffectiveRegion REGION[ 4 ];

		rg_INT numOfIntersections = sortedParams.getSize();
		// 0 intersection
		// This case occurs only when the arc is a major arc
		if(numOfIntersections == 0)
		{
			params.add(0.0);
			params.add(1.0);
			regions.add(regionOfStPt);
		}
		// 1 intersection
		else if(numOfIntersections == 1)
		{
			rg_REAL intersectionParam = sortedParams.getFirstEntity();
			// 2 regions
			if(rg_LT(0.0, intersectionParam) && rg_LT(intersectionParam, 1.0))
			{
				params.add(0.0);
				params.add(intersectionParam);
				params.add(1.0);
				regions.add(regionOfStPt);
				regions.add(regionOfEdPt);
			}
			// 1 region
			else if(rg_ZERO(intersectionParam))
			{
				params.add(0.0);
				params.add(1.0);
				regions.add(regionOfEdPt);
			}
			// 1 region
			else
			{
				params.add(0.0);
				params.add(1.0);
				regions.add(regionOfStPt);
			}
		}
		else
		{
			// 2 intersections
			if(numOfIntersections == 2)
			{
				if(sumOfSplittorIndices == 1)
				{
					if(regionOfStPt == Region_S)
					{					
						REGION[ 0 ] = Region_S;
						REGION[ 1 ] = Region_M;
						REGION[ 2 ] = Region_E;
					}
					else
					{
						REGION[ 0 ] = Region_E;
						REGION[ 1 ] = Region_M;
						REGION[ 2 ] = Region_S;
					}
				}
				else if(sumOfSplittorIndices == 2)
				{
					if(regionOfStPt == Region_E)
					{
						REGION[ 0 ] = Region_E;
						REGION[ 1 ] = Region_S;
						REGION[ 2 ] = Region_M;
					}
					else
					{
						REGION[ 0 ] = Region_M;
						REGION[ 1 ] = Region_S;
						REGION[ 2 ] = Region_E;
					}
				}
				// sumOfSplittorIndices == 3
				else
				{
					if(regionOfStPt == Region_M)
					{
						REGION[ 0 ] = Region_M;
						REGION[ 1 ] = Region_E;
						REGION[ 2 ] = Region_S;
					}
					else
					{
						REGION[ 0 ] = Region_S;
						REGION[ 1 ] = Region_E;
						REGION[ 2 ] = Region_M;
					}
				}
			}
			// 3 intersections
			// (For a major arc only)
			else
			{
				if(firstIntersectedSplittor == OnSplittor_SM)
				{
					REGION[ 0 ] = Region_M;
					REGION[ 1 ] = Region_S;
					REGION[ 2 ] = Region_E;
					REGION[ 3 ] = Region_M;
				}
				else
				{
					REGION[ 0 ] = Region_M;
					REGION[ 1 ] = Region_E;
					REGION[ 2 ] = Region_S;
					REGION[ 3 ] = Region_M;
				}
			}

			rg_INDEX regionIndex = 0;
			if(regionOfStPt == OnSplittor_SM ||
			   regionOfStPt == OnSplittor_ME ||	
			   regionOfStPt == OnSplittor_ES   )
			   regionIndex++;
			rg_REAL  tParam;
			params.add(0.0);
			sortedParams.reset4Loop();
			while(sortedParams.setNext4Loop())
			{
				tParam = sortedParams.getEntity();
				if(rg_LT(0.0, tParam) && rg_LT(tParam, 1.0))
				{
					params.add(tParam);
					regions.add(REGION[regionIndex++]);
				}
			}
			params.add(1.0);
			regions.add(REGION[regionIndex]);			
		}
	}
	return regions.getSize();
}


rg_Point3D TrisectorOfSpaceByArc::computeAngularBisector(const rg_Point3D& vec1, const rg_Point3D& vec2, const rg_Point3D& normalVec)
{
	rg_Point3D angleBisector;
	rg_Point3D vec;
	rg_REAL innerProd;

	angleBisector = vec1 + vec2;
	if(angleBisector != rg_Point3D(0., 0., 0.))
	{		
		angleBisector.normalize();
		vec = vec1.crossProduct(vec2);
		innerProd = normalVec.innerProduct(vec);
	}
	else
	{
		rg_TMatrix3D rot;
		rot.rotateArbitraryAxis(normalVec, rg_PI/2);
		angleBisector = rot * vec1;
		angleBisector.normalize();
		innerProd = 1.0;
	}
	if(rg_NEG(innerProd))
		return (-angleBisector);
	else
		return angleBisector;
}


