#ifndef ARC2D_H
#define ARC2D_H

#include "rg_Circle2D.h"
#include <list>
using namespace std;


class Arc2D : public rg_Circle2D
{
private:
    rg_Point2D  m_startPoint;
    rg_Point2D  m_endPoint;

public:
    Arc2D();
    Arc2D(const rg_Circle2D& circle, const rg_Point2D& startPoint, const rg_Point2D& endPoint);
    Arc2D(const Arc2D& arc);
    ~Arc2D();

    inline rg_Point2D getStartPoint() const { return m_startPoint; }
    inline rg_Point2D getEndPoint() const { return m_endPoint; }

    inline void setStartPoint(const rg_Point2D& startPoint) { m_startPoint = startPoint; }
    inline void setEndPoint(const rg_Point2D& endPoint) { m_endPoint = endPoint; }

    Arc2D& operator =(const Arc2D& arc);
    rg_BOOL operator==(const Arc2D& arc) const;


    rg_REAL angle() const;

    rg_BOOL isOnThisArc(const rg_Point2D& point) const;
    rg_BOOL isContainedIn(const Arc2D& bigArc) const;
    rg_INT  intersect(const Arc2D& arc, rg_Point2D intersection[]) const;

    void evaluatePointsOnArcGivenResolution(const rg_INT&  resolution, list<rg_Point2D>& pointsOnArc) const;

};

#endif

