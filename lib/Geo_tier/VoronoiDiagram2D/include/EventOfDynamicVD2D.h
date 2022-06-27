#ifndef _EVENT_OF_DYNAMIC_VD2D_H_
#define _EVENT_OF_DYNAMIC_VD2D_H_

namespace V {
namespace GeometryTier {


class EventOfDynamicVD2D
{
public:
    enum EventType { EDGE_FLIP, DISK_COLLISION, DISK_VELOCITY_CHANGE, DISK_HOPPING };

protected:
	double      m_occurringTime;
	EventType   m_eventType;

public:
    inline virtual ~EventOfDynamicVD2D() {};

	//setter
    inline void set_occurring_time(const double& occuringTime) { m_occurringTime = occuringTime; };

	//getter
    inline double     get_occurring_clock() const { return m_occurringTime; };
    inline EventType  get_event_type()     const   { return m_eventType; };

protected:
    void copy_from(const EventOfDynamicVD2D& eventOfDVD);
};


inline void EventOfDynamicVD2D::copy_from(const EventOfDynamicVD2D& eventOfDVD)
{
    m_eventType     = eventOfDVD.m_eventType;
    m_occurringTime = eventOfDVD.m_occurringTime;
}



}
}

#endif