#include "Line3D.h"
#include "rg_RelativeOp.h"


Line3D::Line3D()
{
}



Line3D::Line3D(const rg_Point3D& directionVec, const rg_Point3D& passingPt)
{
	m_directionVec  = directionVec;
	m_passingPt     = passingPt;
}



Line3D::Line3D(const Line3D& line)
{
	m_directionVec  = line.m_directionVec;
	m_passingPt     = line.m_passingPt;
}



Line3D::~Line3D()
{
}



rg_Point3D Line3D::getDirVec() const
{
	return m_directionVec;
}



rg_Point3D Line3D::getPassPt() const
{
	return m_passingPt;
}



void Line3D::set(const rg_Point3D& directionVec, const rg_Point3D& passingPt)
{
	m_directionVec  = directionVec;
	m_passingPt     = passingPt;
}



void Line3D::setDirVec(const rg_Point3D& directionVec)
{
	m_directionVec = directionVec;
}



void Line3D::setPassPt(const rg_Point3D& passingPt)
{
	m_passingPt = passingPt;
}



Line3D& Line3D::operator =(const Line3D &line)
{
    if(this == &line) {
		return *this;
    }

	m_directionVec  = line.m_directionVec;
	m_passingPt     = line.m_passingPt;

	return *this;
}



rg_BOOL Line3D::isOnLine(const rg_Point3D& point, const rg_REAL& res) const
{
    rg_Point3D vectorFromPassingToPoint = point - m_passingPt;
    vectorFromPassingToPoint.normalize();

    rg_Point3D unitDirectionVector = m_directionVec.getUnitVector();

    rg_REAL innerProduct = unitDirectionVector.innerProduct(vectorFromPassingToPoint);

    if ( rg_EQ(innerProduct, 1.0, res) || rg_EQ(innerProduct, -1.0, res) ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



rg_Point3D Line3D::computeProjectedPointForPoint3D(const rg_Point3D& point) const
{
	rg_Point3D vec      = point - m_passingPt;
	rg_Point3D unitDir  = m_directionVec.getUnitVector();

	rg_REAL    distance       = vec.innerProduct(unitDir);
    rg_Point3D projectedPoint = m_passingPt + (distance * unitDir);

	return projectedPoint;	
}



rg_Point3D Line3D::computeProjectedPointForPoint3D(const rg_Point3D& point, rg_REAL& dist) const
{
	rg_Point3D vec     = point - m_passingPt;
	rg_Point3D unitDir = m_directionVec.getUnitVector();

	dist = vec.innerProduct(unitDir);
    rg_Point3D projectedPoint = m_passingPt + (dist * unitDir);

	return projectedPoint;	
}



rg_REAL Line3D::computeDistanceFromPoint3D(const rg_Point3D& point) const
{
	rg_Point3D projectedPt = computeProjectedPointForPoint3D(point);
	rg_REAL    distance    = point.distance(projectedPt);

	return distance;
}



rg_REAL Line3D::computeMinDistFromPoint(const rg_Point3D& targetPt, 
										rg_Point3D& minPt, 
										rg_REAL& minPtParam) const
{
	minPtParam =   (targetPt - m_passingPt).innerProduct(m_directionVec) 
                 / m_directionVec.innerProduct(m_directionVec);
	minPt      = m_passingPt + minPtParam * m_directionVec;

	rg_REAL minDist = targetPt.distance(minPt);

	return minDist;
}



