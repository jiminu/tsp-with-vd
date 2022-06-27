#ifndef _DISK_VELOCITY_CHANGE_EVENT_
#define _DISK_VELOCITY_CHANGE_EVENT_

#include "EventOfDynamicVD2D.h"
#include "Generator2D.h"

namespace V {
namespace GeometryTier {


class DiskVelocityChangeEvent : public EventOfDynamicVD2D
{
private:
	Generator2D* m_DiskGenerator;
	rg_Point2D   m_VelocityVectorAfterVelocityChange;
	rg_Point2D   m_VelocityVectorBeforeVelocityChange;
	rg_Point2D   m_DiskCenterAtVelocityChangeTime;
    bool         m_VelocityVectorSetting;
	
public:
	//constructor
	DiskVelocityChangeEvent();
    DiskVelocityChangeEvent(Generator2D* diskGenerator);
    DiskVelocityChangeEvent(const double& occuringTime, Generator2D* diskGenerator);
    DiskVelocityChangeEvent(const double& occuringTime, Generator2D* diskGenerator, const rg_Point2D& newVelocityVector);
	DiskVelocityChangeEvent(const DiskVelocityChangeEvent& diskVelocityChangeEvent);

	//deconstructor
	~DiskVelocityChangeEvent();

	//operator
	DiskVelocityChangeEvent& operator=(const DiskVelocityChangeEvent& diskVelocityChangeEvent);

    virtual void copy_from(const DiskVelocityChangeEvent& diskVelocityChangeEvent);

	//setter
    inline void set_velocity_changed_generator(Generator2D* diskGenerator) { m_DiskGenerator = diskGenerator;};
	
    inline void set_velocity_vector_after_change(const rg_Point2D& velocityVector) { m_VelocityVectorAfterVelocityChange = velocityVector; };

    inline void set_velocity_vector_before_change(const rg_Point2D& velocityVector){ m_VelocityVectorBeforeVelocityChange = velocityVector;};
   
    inline void set_disk_center_at_velocity_change_time(const rg_Point2D& center) {m_DiskCenterAtVelocityChangeTime = center; };

    inline void set_velocity_vector_setting(const bool& velocityVectorSetting) { m_VelocityVectorSetting = velocityVectorSetting; };
	
    //getter
    inline Generator2D* get_disk_generator() const { return m_DiskGenerator; };

    inline const rg_Point2D& get_future_velocity_vector_of_disk() const { return m_VelocityVectorAfterVelocityChange; };
  
    inline const rg_Point2D& get_past_velocity_vector_of_disk() const { return m_VelocityVectorBeforeVelocityChange; };

    inline const rg_Point2D& get_center_of_disk() const { return m_DiskCenterAtVelocityChangeTime; };
  
    inline bool  get_velocity_vector_setting() const { return m_VelocityVectorSetting; };

};

}
}

#endif //_DISK_VELOCITY_CHANGE_EVENT_
