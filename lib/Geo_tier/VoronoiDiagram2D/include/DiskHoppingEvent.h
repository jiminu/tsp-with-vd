#ifndef _DELETE_N_INSERT_DISK_EVENT_H_
#define _DELETE_N_INSERT_DISK_EVENT_H_

#include "EventOfDynamicVD2D.h"
#include "Generator2D.h"
#include "rg_Point2D.h"

namespace V {
namespace GeometryTier {


class DiskHoppingEvent : public EventOfDynamicVD2D
{
private:
    Generator2D* m_TargetGenerator;
    Generator2D* m_InWhatContainerThisGeneratorIsDeletedNInserted[2];
    rg_Point2D   m_DeletedNInsertedPosition[2];
    bool		 m_doesThisEventHappen;

    rg_Point2D   m_VelocityVectorBeforeCollisionOfTargetGenerator;
    rg_Point2D   m_VelocityVectorAfterCollisionOfTargetGenerator;
    bool         m_VelocityVectorSetting;

public:
    //constructor
    DiskHoppingEvent();
    DiskHoppingEvent(const double& occuringTime, 
                           Generator2D* targetGenerator,
                           Generator2D* containerInWhichGeneratorIsDeleted, 
                           Generator2D* containerInWhichGeneratorIsInserted);
    DiskHoppingEvent(const DiskHoppingEvent& diskHoppingEvent);

    //deconstructor
    ~DiskHoppingEvent();

    //operator
    DiskHoppingEvent& operator=(const DiskHoppingEvent& diskHoppingEvent);

    //setter
    inline void set_target_generator(Generator2D* targetGenerator) {
        m_TargetGenerator = targetGenerator;
    };


    inline void set_in_what_container_this_generator_is_deleted_n_inserted(
        Generator2D* containerInWhichGeneratorIsDeleted, 
        Generator2D* containerInWhichGeneratorIsInserted) {
        m_InWhatContainerThisGeneratorIsDeletedNInserted[0] = containerInWhichGeneratorIsDeleted;
        m_InWhatContainerThisGeneratorIsDeletedNInserted[1] = containerInWhichGeneratorIsInserted;
    };


    inline void set_deleted_n_inserted_position(
        const rg_Point2D& deletedPosition, const rg_Point2D& insertedPosition) {
        m_DeletedNInsertedPosition[0] = deletedPosition;
        m_DeletedNInsertedPosition[1] = insertedPosition;
    };


    inline void set_velocity_vector_change_position(const rg_Point2D& position) {
        m_DeletedNInsertedPosition[0] = position;
        m_DeletedNInsertedPosition[1] = position;
    };


    inline void set_inserted_position(const rg_Point2D& insertedPosition) {
        m_DeletedNInsertedPosition[1] = insertedPosition;
    };

    inline void set_does_this_event_happen(const bool& doesThisEventHappen) { m_doesThisEventHappen = doesThisEventHappen; };

    inline void set_velocity_vector_after_collision(const rg_Point2D& velocityVector) { m_VelocityVectorAfterCollisionOfTargetGenerator = velocityVector; };

    inline void set_velocity_vector_before_collision(const rg_Point2D& velocityVector) { m_VelocityVectorBeforeCollisionOfTargetGenerator = velocityVector; };

    inline void set_velocity_vector_setting(const bool& velocityVectorSetting) { m_VelocityVectorSetting = velocityVectorSetting; };


    //getter
    inline Generator2D* get_target_generator() const { return m_TargetGenerator; };
    inline Generator2D* get_in_what_container_this_generator_is_deleted() const { return m_InWhatContainerThisGeneratorIsDeletedNInserted[0]; };
    inline Generator2D* get_in_what_container_this_generator_is_inserted() const { return m_InWhatContainerThisGeneratorIsDeletedNInserted[1]; };
    inline const rg_Point2D& get_deleted_position() const { return m_DeletedNInsertedPosition[0]; };
    inline const rg_Point2D& get_inserted_position() const { return m_DeletedNInsertedPosition[1]; };
    inline bool get_does_this_event_happen() { return m_doesThisEventHappen; };
    inline const rg_Point2D& get_future_velocity_vector_of_disk() const { return m_VelocityVectorAfterCollisionOfTargetGenerator; };
    inline const rg_Point2D& get_past_velocity_vector_of_disk()   const { return m_VelocityVectorBeforeCollisionOfTargetGenerator; };
    inline bool  get_velocity_vector_setting() const { return m_VelocityVectorSetting; };

    void copy_from(const DiskHoppingEvent& diskHoppingEvent);
};

}
}


#endif