#ifndef _ARC3D_H
#define _ARC3D_H

#include "rg_Point3D.h"
#include "Circle3D.h"
#include "Plane.h"
#include "rg_dList.h"
#include "DistEventOnLineSegAndArc.h"


enum ARCType {MINOR_ARC, MAJOR_ARC};

class Arc3D : public Circle3D
{
private:
    rg_Point3D  m_startPoint;
    rg_Point3D  m_endPoint;

public:
    Arc3D();
    Arc3D(const Circle3D& circle, const rg_Point3D& startPt, const rg_Point3D& endPt);
    Arc3D(const Arc3D& arc);
    ~Arc3D();

    inline rg_Point3D  getStartPoint() const { return m_startPoint; }
    inline rg_Point3D  getEndPoint() const   { return m_endPoint; }

    void        setStartPoint(const rg_Point3D& startPt);
    void        setEndPoint(const rg_Point3D& endPt);
    void        setPoints(const rg_Point3D& startPt, const rg_Point3D& endPt);
    void        setArc(const Circle3D& circle, const rg_Point3D& startPt, const rg_Point3D& endPt);

    Arc3D&      operator =(const Arc3D& arc);
    rg_BOOL     operator==(const Arc3D& arc) const;

    rg_BOOL     isOnThisArc(const rg_Point3D& point) const;
    rg_BOOL     isMinorArc() const;
    rg_BOOL     isContainedIn(const Arc3D& bigArc) const;
	ARCType     getArcType() const;

    void        reverse();
    Arc3D       getReversedArc() const;
	Circle3D    getCircle() const;

    rg_INT      intersect(const Plane& plane, rg_Point3D* intersectPts) const;
    rg_INT      intersect(const Circle3D& circle, rg_Point3D* intersectPts) const;
	rg_INT      intersect(const Arc3D& arc, rg_Point3D* intersectPts, const rg_REAL tolerance = resNeg6) const;
    
    rg_REAL     evaluateAngle() const;
    rg_REAL     evaluateArcLength() const;
    rg_Point3D  evaluateMiddlePoint() const;

    void        evaluatePointsOnArc(const rg_REAL& unitLength, rg_dList<rg_Point3D>& pointsOnArc) const;
    void        evaluatePointsOnArc(const rg_INT&  depth, rg_Point3D*& pointsOnArc) const;
    void        evaluatePointsOnArc(const rg_INT&  depth, rg_dList<rg_Point3D>& pointsOnArc) const;
    void        evaluatePointsOnArcGivenNumPts(const rg_INT&  numPts, rg_Point3D*& pointsOnArc) const;
    void        evaluatePointsOnArcGivenNumPts(const rg_INT&  numPts, rg_dList<rg_Point3D>& pointsOnArc) const;


	rg_REAL computeMinDistFromPoint(const rg_Point3D& targetPt, rg_Point3D& minPtOnArc) const;
	rg_REAL computeMinDistanceFromLineSegment(const LineSegment3D& lineSegment, rg_Point3D& minPtOnArc, rg_Point3D& minPtOnLineSeg, PosOnLineSegOrArc& pos) const;


	rg_INT locateLocalMinimaOfDistFromLineSegment(const LineSegment3D& lineSegment, 
		                                          rg_Point3D*& pointsOnArc,
											      rg_Point3D*& pointsOnLineSegment,
											      PosOnLineSegOrArc*& positionsOnLineSegment) const;
	
	rg_INT locateLocalMinimaOfDistFromLineSegment(const LineSegment3D& lineSegment, 
		                                          rg_Point3D*& pointsOnArc,
		                                          rg_Point3D*& pointsOnLineSegment,
												  PosOnLineSegOrArc*& positionsOnArc,
											      PosOnLineSegOrArc*& positionsOnLineSegment) const;

	rg_INT locateLocalMinimaOfDistFromLineSegment(const LineSegment3D& lineSegment, 
		                                          rg_dList<rg_Point3D>& pointsOnArc,
		                                          rg_dList<rg_Point3D>& pointsOnLineSegment,
											      rg_dList<PosOnLineSegOrArc>& positionsOnLineSegment) const;

	rg_INT locateLocalMinimaOfDistFromLineSegment(const LineSegment3D& lineSegment, 
		                                          rg_dList<rg_Point3D>& pointsOnArc,
		                                          rg_dList<rg_Point3D>& pointsOnLineSegment,
												  rg_dList<PosOnLineSegOrArc>& positionsOnArc,
											      rg_dList<PosOnLineSegOrArc>& positionsOnLineSegment) const;

	rg_INT locateInflectionPointsOfDistFromLineSegment(const LineSegment3D& lineSegment, 
		                                               rg_Point3D*& pointsOnArc,
		                                               rg_Point3D*& pointsOnLineSegment,
												       PosOnLineSegOrArc*& positionsOnArc,
											           PosOnLineSegOrArc*& positionsOnLineSegment) const;

	rg_INT locateInflectionPointsOfDistFromLineSegment(const LineSegment3D& lineSegment, 
		                                               rg_dList<rg_Point3D>& pointsOnArc,
													   rg_dList<rg_Point3D>& pointsOnLineSegment,
													   rg_dList<PosOnLineSegOrArc>& positionsOnArc,
													   rg_dList<PosOnLineSegOrArc>& positionsOnLineSegment) const;


    rg_INT locateLocalMinimaOfDistFromLineSegment(  const LineSegment3D&                lineSegment,
                                                    rg_dList<DistEventOnLineSegAndArc>& distEventList );

    rg_INT locateMinEventsOnDistanceFromLineSegment(const LineSegment3D&                lineSegment,
                                                    rg_dList<DistEventOnLineSegAndArc>& distEventList );
};

#endif


