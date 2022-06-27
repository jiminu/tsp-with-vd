#include "LineSegment3D.h"

LineSegment3D::LineSegment3D()
: Line3D()
{
}



LineSegment3D::LineSegment3D(const rg_Point3D& sPt, const rg_Point3D& ePt)
: Line3D(ePt - sPt, sPt)
{
	m_sPoint = sPt;
	m_ePoint = ePt;
}



LineSegment3D::LineSegment3D(const LineSegment3D& lineSegment)
: Line3D( lineSegment )
{
	m_sPoint = lineSegment.m_sPoint;
	m_ePoint = lineSegment.m_ePoint;
}



LineSegment3D::~LineSegment3D()
{
}


	

rg_Point3D LineSegment3D::getStartPt() const
{
	return m_sPoint;
}



rg_Point3D LineSegment3D::getEndPt() const
{
	return m_ePoint;
}



rg_REAL LineSegment3D::getLength() const
{
	return (m_ePoint - m_sPoint).distance();
}



void LineSegment3D::setPoints(const rg_Point3D& sPt, const rg_Point3D& ePt)
{
	Line3D::set(ePt - sPt, sPt);

	m_sPoint = sPt;
	m_ePoint = ePt;
}



LineSegment3D& LineSegment3D::operator =(const LineSegment3D& lineSegment)
{
    if ( this == &lineSegment ) {
        return *this;
    }

    Line3D::operator=(lineSegment);

    m_sPoint = lineSegment.m_sPoint;
    m_ePoint = lineSegment.m_ePoint;

    return *this;
}



rg_BOOL LineSegment3D::isOnLineSegment(const rg_Point3D& point, const rg_REAL& res) const
{
    if ( Line3D::isOnLine( point, res ) ) {
        rg_REAL param = (point.getX() - m_ePoint.getX())/(m_sPoint.getX() - m_ePoint.getX());

        if ( rg_GE(param, 0.0, res) && rg_LE(param, 1.0, res) ) {
            return rg_TRUE;
        }
        else {
            return rg_FALSE;
        }
    }
    else {
        return rg_FALSE;
    }
}



rg_Point3D LineSegment3D::evaluatePt(const rg_REAL& param) const
{
	return m_sPoint + param * (m_ePoint - m_sPoint);
}



rg_Point3D LineSegment3D::computeProjectedPointForPoint3D(const rg_Point3D& target, rg_REAL& param) const
{
	rg_Point3D projection = Line3D::computeProjectedPointForPoint3D(target, param);

	return projection;
}



rg_REAL LineSegment3D::computeMinDistFromPoint(const rg_Point3D& targetPt, 
											   rg_Point3D& minPt, 
											   rg_REAL& minPtParam, 
											   PosOnLineSegOrArc& pos) const
{
	rg_REAL minDist = Line3D::computeMinDistFromPoint(targetPt, minPt, minPtParam);
	if(rg_LE(minPtParam, 0.0))
	{
		minPtParam = 0.0;
		minPt = m_sPoint;
		pos = START_PT;
		return m_sPoint.distance(targetPt);
	}
	else if(rg_GE(minPtParam, 1.0))
	{
		minPtParam = 1.0;
		minPt = m_ePoint;
		pos = END_PT;
		return m_ePoint.distance(targetPt);
	}
	pos = NONEXTREME_PT;
	return minDist;
}



rg_REAL LineSegment3D::computeMinDistFromPoint(const rg_Point3D& targetPt) const
{
    rg_Point3D ptWithMinDist;
    rg_REAL    paramOnPtWithMinDist;

	rg_REAL minDist = Line3D::computeMinDistFromPoint(targetPt, ptWithMinDist, paramOnPtWithMinDist);
	
    if( rg_LE(paramOnPtWithMinDist, 0.0) ){
		return m_sPoint.distance(targetPt);
	}
	else if( rg_GE(paramOnPtWithMinDist, 1.0) )	{
		return m_ePoint.distance(targetPt);
	}
    else {
	    return minDist;
    }
}



rg_REAL LineSegment3D::computeMinDistFromPoint(const rg_Point3D& targetPt, 
											   rg_Point3D& minPt, 
											   rg_REAL& minPtParam) const
{
	rg_REAL minDist = Line3D::computeMinDistFromPoint(targetPt, minPt, minPtParam);
	if(rg_LE(minPtParam, 0.0))
	{
		minPtParam = 0.0;
		minPt = m_sPoint;
		return m_sPoint.distance(targetPt);
	}
	else if(rg_GE(minPtParam, 1.0))
	{
		minPtParam = 1.0;
		minPt = m_ePoint;
		return m_ePoint.distance(targetPt);
	}
	return minDist;
}

rg_REAL LineSegment3D::computeMaxDistFromPoint(const rg_Point3D& targetPt, 
											   rg_Point3D& maxPt, 
											   PosOnLineSegOrArc& pos) const
{
	rg_REAL sDist = m_sPoint.distance(targetPt);
	rg_REAL eDist = m_ePoint.distance(targetPt);
	if(rg_GE(sDist, eDist))
	{
		maxPt = m_sPoint;
		pos = START_PT;
		return sDist;
	}
	else
	{
		maxPt = m_ePoint;
		pos = END_PT;
		return eDist;
	}
}

