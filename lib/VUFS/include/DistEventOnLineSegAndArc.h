#ifndef _DISTEVENTONLINESEGANDARC_H
#define _DISTEVENTONLINESEGANDARC_H

#include "rg_Point3D.h"
#include "LineSegment3D.h"

enum DistEventType { DET_UNKNOWN, DET_LOCAL_MIN, DET_LOCAL_MAX, 
                     DET_INFLECTION_POINT,
                     DET_ROOT_THIRD_DERIVATIVE };

class DistEventOnLineSegAndArc
{
private:
    DistEventType       m_eventType;
    rg_Point3D          m_pointOnArc;
    rg_Point3D          m_pointOnLineSegment;
    PosOnLineSegOrArc   m_positionOnArc;
    PosOnLineSegOrArc   m_positionOnLineSegment;

public:
    DistEventOnLineSegAndArc();
    DistEventOnLineSegAndArc(const DistEventType&     eventType, 
                             const rg_Point3D&        pointOnArc,
                             const rg_Point3D&        pointOnLineSegment,
                             const PosOnLineSegOrArc& positionOnArc,
                             const PosOnLineSegOrArc& positionOnLineSegment);
    DistEventOnLineSegAndArc(const DistEventOnLineSegAndArc& distEvent);
    ~DistEventOnLineSegAndArc();

    DistEventType       getEventType() const;
    rg_Point3D          getPointOnArc() const;
    rg_Point3D          getPointOnLineSegment() const;
    PosOnLineSegOrArc   getPositionOnArc() const;
    PosOnLineSegOrArc   getPositionOnLineSegment() const;

    void    setEventType(            const DistEventType&     eventType);
    void    setPointOnArc(           const rg_Point3D&        pointOnArc);
    void    setPointOnLineSegment(   const rg_Point3D&        pointOnLineSegment);
    void    setPositionOnArc(        const PosOnLineSegOrArc& positionOnArc);
    void    setPositionOnLineSegment(const PosOnLineSegOrArc& positionOnLineSegment);
    void    setDistEvent(const DistEventType&     eventType, 
                         const rg_Point3D&        pointOnArc,
                         const rg_Point3D&        pointOnLineSegment,
                         const PosOnLineSegOrArc& positionOnArc,
                         const PosOnLineSegOrArc& positionOnLineSegment);

    DistEventOnLineSegAndArc& operator =(const DistEventOnLineSegAndArc& distEvent);
};

#endif

