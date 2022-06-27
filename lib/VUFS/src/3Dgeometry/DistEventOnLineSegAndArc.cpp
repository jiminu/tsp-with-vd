#include "DistEventOnLineSegAndArc.h"

DistEventOnLineSegAndArc::DistEventOnLineSegAndArc()
: m_eventType(DET_UNKNOWN),
  m_positionOnArc(UNKNOWN_POS),
  m_positionOnLineSegment(UNKNOWN_POS)
{
}



DistEventOnLineSegAndArc::DistEventOnLineSegAndArc(const DistEventType&     eventType, 
                                                   const rg_Point3D&        pointOnArc,
                                                   const rg_Point3D&        pointOnLineSegment,
                                                   const PosOnLineSegOrArc& positionOnArc,
                                                   const PosOnLineSegOrArc& positionOnLineSegment)
: m_eventType(eventType),
  m_pointOnArc(pointOnArc),
  m_pointOnLineSegment(pointOnLineSegment),
  m_positionOnArc(positionOnArc),
  m_positionOnLineSegment(positionOnLineSegment)
{
}



DistEventOnLineSegAndArc::DistEventOnLineSegAndArc(const DistEventOnLineSegAndArc& distEvent)
: m_eventType(distEvent.m_eventType),
  m_pointOnArc(distEvent.m_pointOnArc),
  m_pointOnLineSegment(distEvent.m_pointOnLineSegment),
  m_positionOnArc(distEvent.m_positionOnArc),
  m_positionOnLineSegment(distEvent.m_positionOnLineSegment)
{
}



DistEventOnLineSegAndArc::~DistEventOnLineSegAndArc()
{
}




DistEventType       DistEventOnLineSegAndArc::getEventType() const
{
    return m_eventType;
}



rg_Point3D          DistEventOnLineSegAndArc::getPointOnArc() const
{
    return m_pointOnArc;
}



rg_Point3D          DistEventOnLineSegAndArc::getPointOnLineSegment() const
{
    return m_pointOnLineSegment;
}



PosOnLineSegOrArc   DistEventOnLineSegAndArc::getPositionOnArc() const
{
    return m_positionOnArc;
}



PosOnLineSegOrArc   DistEventOnLineSegAndArc::getPositionOnLineSegment() const
{
    return m_positionOnLineSegment;
}




void DistEventOnLineSegAndArc::setEventType(const DistEventType&     eventType)
{
    m_eventType = eventType;
}



void DistEventOnLineSegAndArc::setPointOnArc(const rg_Point3D&        pointOnArc)
{
    m_pointOnArc = pointOnArc;
}



void DistEventOnLineSegAndArc::setPointOnLineSegment(const rg_Point3D&        pointOnLineSegment)
{
    m_pointOnLineSegment = pointOnLineSegment;
}



void DistEventOnLineSegAndArc::setPositionOnArc(const PosOnLineSegOrArc& positionOnArc)
{
    m_positionOnArc = positionOnArc;
}



void DistEventOnLineSegAndArc::setPositionOnLineSegment(const PosOnLineSegOrArc& positionOnLineSegment)
{
    m_positionOnLineSegment = positionOnLineSegment;
}



void DistEventOnLineSegAndArc::setDistEvent(const DistEventType&     eventType, 
                                            const rg_Point3D&        pointOnArc,
                                            const rg_Point3D&        pointOnLineSegment,
                                            const PosOnLineSegOrArc& positionOnArc,
                                            const PosOnLineSegOrArc& positionOnLineSegment)
{
    m_eventType             = eventType;
    m_pointOnArc            = pointOnArc;
    m_pointOnLineSegment    = pointOnLineSegment;
    m_positionOnArc         = positionOnArc;
    m_positionOnLineSegment = positionOnLineSegment;
}




DistEventOnLineSegAndArc& DistEventOnLineSegAndArc::operator =(const DistEventOnLineSegAndArc& distEvent)
{
    if ( this == &distEvent ) {
        return *this;
    }

    m_eventType             = distEvent.m_eventType;
    m_pointOnArc            = distEvent.m_pointOnArc;
    m_pointOnLineSegment    = distEvent.m_pointOnLineSegment;
    m_positionOnArc         = distEvent.m_positionOnArc;
    m_positionOnLineSegment = distEvent.m_positionOnLineSegment;

    return *this;
}



