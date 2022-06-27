////////////////////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : TrisectorOfSpaceByArc.h
//	  
//    DESCRIPTION : 
//           This is the implementation of the class TrisectorOfSpaceByArc
//           and is used for representing the trisectioning of the space by an arc in 3d.
//
//    Reference: CAD ±³¼ö´Ô ³í¹®
//               Trisector, Splittor, effective region °³³ä comment
//
//    AUTHOR      : Ryu, Jooghyun
//    START DATE  : Jan. 14, 2010
//
//              
//           Copyright ¨Ï 2010 by Voronoi Diagram Research Center, Hanyang University
//
/////////////////////////////////////////////////////////////////////////////////////


#ifndef TRISECTOR_OF_SPCACE_BY_ARC
#define TRISECTOR_OF_SPCACE_BY_ARC

#include "Arc3D.h"
//#include "UtilityFuncsForMoleSurf.h"


enum EffectiveRegion {RegionUnknown, Region_M = 0, Region_E, Region_S, 
                      OnSplittor_SM, OnSplittor_ME, OnSplittor_ES, 
					  OnAxisOfArc};

const EffectiveRegion REGION[ 7 ] = {Region_M, Region_E, Region_S, 
		                             OnSplittor_SM, OnSplittor_ME, OnSplittor_ES, 
						             OnAxisOfArc};

enum NumberOfOccupiedRegion {ONE, TWO_THREE_OR_TO_BE_VERIFIED};

class TrisectorOfSpaceByArc
{
private:
	Plane   m_splittor[ 3 ]; // [ 0 ]: Splittor_SM, [ 1 ]: Splittor_ME, [ 2 ]; Splittor_ES
	ARCType m_arcType;       // minor or major arc

public:
	TrisectorOfSpaceByArc();
	TrisectorOfSpaceByArc(const Arc3D& arc);
	
	rg_INT identifyEffectiveRegionsOfLineSegment(const LineSegment3D& lineSeg,
		                                         rg_dList<rg_REAL>& params, 
		                                         rg_dList<EffectiveRegion>& regions);
	
		NumberOfOccupiedRegion getEffectiveRegion(const LineSegment3D& lineSeg, 
												  EffectiveRegion& regionOfStPt, 
												  EffectiveRegion& regionOfEdPt);	
			EffectiveRegion getEffectiveRegion(const rg_Point3D& pt);
		void             lcoateEffectiveRegion(EffectiveRegion& regionOfStPt, 
			                                   EffectiveRegion& regionOfEdPt);

		rg_INT computeIntersectionWithSplittors(const LineSegment3D& lineSeg,
											    rg_dList<rg_REAL>& unsortedParams, 
												rg_INDEX& sumOfSplittorIndices,
												EffectiveRegion& firstIntersectedSplittor);

		void   locateEffectiveRegions(EffectiveRegion effectiveRegion[], 
			                          rg_INT& numOfEffectiveRegions,
									  const EffectiveRegion& regionOfStPt,
									  const EffectiveRegion& regionOfEdPt,
			                          const rg_INT& numOfIntersections,  
									  const rg_INDEX sumOfSplittorIndices,
									  const EffectiveRegion& firstIntersectedSplittor);

		void sortIntersectionParams(rg_dList<rg_REAL>& unsortedParams, 
			                        rg_dList<rg_REAL>& sortedParams);

	rg_INT identifyRelativePosOfLineSegment_old(const LineSegment3D& lineSeg,
		                                        rg_dList<rg_REAL>& params, 
		                                        rg_dList<EffectiveRegion>& regions);

    rg_Point3D computeAngularBisector(const rg_Point3D& vec1, const rg_Point3D& vec2, const rg_Point3D& normalVec);

};

#endif
