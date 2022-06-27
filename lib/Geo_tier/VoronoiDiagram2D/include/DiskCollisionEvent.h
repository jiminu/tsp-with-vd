#ifndef _DISK_COLLISION_EVENT_H_
#define _DISK_COLLISION_EVENT_H_

#include "EventOfDynamicVD2D.h"
#include "Generator2D.h"
#include "DynamicDisk.h"
#include "rg_Point2D.h"

namespace V {
namespace GeometryTier {


class DiskCollisionEvent : public EventOfDynamicVD2D
{
private:
	Generator2D* m_CollidedDiskGenerator[2];
	rg_Point2D   m_MotionVectorsAfterCollision[2];
	rg_Point2D   m_MotionVectorsBeforeCollision[2];
	rg_Point2D   m_DiskCentersAtCollisionTime[2];
    bool         m_VelocityVectorSetting;
	
public:
	//constructor
	DiskCollisionEvent();
    DiskCollisionEvent(Generator2D* diskGenerator1, Generator2D* diskGenerator2);
    DiskCollisionEvent(const double& occuringTime, Generator2D* diskGenerator1, Generator2D* diskGenerator2);
	DiskCollisionEvent(const DiskCollisionEvent& diskCollisionEvent);

	//deconstructor
	~DiskCollisionEvent();

	//operator
	DiskCollisionEvent& operator=(const DiskCollisionEvent& diskCollisionEvent);

	//setter
    inline void set_two_collided_generators(Generator2D* diskGenerator1, Generator2D* diskGenerator2) { m_CollidedDiskGenerator[0] = diskGenerator1; 
                                                                                                        m_CollidedDiskGenerator[1] = diskGenerator2;};
	
    inline void set_velocity_vectors_after_collision(const rg_Point2D& motionVector1, const rg_Point2D& motionVector2) { m_MotionVectorsAfterCollision[0] = motionVector1;
                                                                                                                       m_MotionVectorsAfterCollision[1] = motionVector2;};


    inline void set_velocity_vectors_before_collision(const rg_Point2D& motionVector1, const rg_Point2D& motionVector2){ m_MotionVectorsBeforeCollision[0] = motionVector1;
                                                                                                                       m_MotionVectorsBeforeCollision[1] = motionVector2;};
   
    inline void set_disk_centers_at_collision_time(const rg_Point2D& center1, const rg_Point2D& center2){ m_DiskCentersAtCollisionTime[0] = center1;
                                                                                                          m_DiskCentersAtCollisionTime[1] = center2;};

    inline void set_velocity_vector_setting(const bool& velocityVectorSetting) { m_VelocityVectorSetting = velocityVectorSetting; };
	
    //getter
    inline Generator2D* get_disk_generator1() const { return m_CollidedDiskGenerator[0]; };
    inline Generator2D* get_disk_generator2() const { return m_CollidedDiskGenerator[1]; };

    inline DynamicDisk* get_dynamic_disk1() const { return static_cast<DynamicDisk*>(get_disk_generator1()->getUserData()); };
    inline DynamicDisk* get_dynamic_disk2() const { return static_cast<DynamicDisk*>(get_disk_generator2()->getUserData()); };

    inline const rg_Point2D& get_future_velocity_vector_of_disk1() const { return m_MotionVectorsAfterCollision[0]; };
    inline const rg_Point2D& get_future_velocity_vector_of_disk2() const { return m_MotionVectorsAfterCollision[1]; };

    inline const rg_Point2D& get_past_velocity_vector_of_disk1() const { return m_MotionVectorsBeforeCollision[0]; };
    inline const rg_Point2D& get_past_velocity_vector_of_disk2() const { return m_MotionVectorsBeforeCollision[1]; };

    inline const rg_Point2D& get_center_of_disk1() const { return m_DiskCentersAtCollisionTime[0]; };
    inline const rg_Point2D& get_center_of_disk2() const { return m_DiskCentersAtCollisionTime[1]; };

    inline bool  get_velocity_vector_setting() const { return m_VelocityVectorSetting; };

private:
    void copy_from(const DiskCollisionEvent& diskCollisionEvent);
};


}
}

#endif