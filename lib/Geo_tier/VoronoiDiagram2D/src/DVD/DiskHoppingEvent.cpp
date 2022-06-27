#include "DiskHoppingEvent.h"
using namespace V::GeometryTier;

DiskHoppingEvent::DiskHoppingEvent()
{
    m_eventType       = DISK_HOPPING;
    m_occurringTime   = -1.0;
    m_TargetGenerator = NULL;
    m_InWhatContainerThisGeneratorIsDeletedNInserted[0] = NULL;
    m_InWhatContainerThisGeneratorIsDeletedNInserted[1] = NULL;
    m_VelocityVectorSetting = false;
}


DiskHoppingEvent::DiskHoppingEvent(const double & occuringTime, Generator2D * targetGenerator, Generator2D * containerInWhichGeneratorIsDeleted, Generator2D * containerInWhichGeneratorIsInserted)
{
    m_eventType       = DISK_HOPPING;
    m_occurringTime   = occuringTime;
    m_TargetGenerator = targetGenerator;
    m_InWhatContainerThisGeneratorIsDeletedNInserted[0] = containerInWhichGeneratorIsDeleted;
    m_InWhatContainerThisGeneratorIsDeletedNInserted[1] = containerInWhichGeneratorIsInserted;
    m_VelocityVectorSetting = false;
}


DiskHoppingEvent::DiskHoppingEvent(const DiskHoppingEvent & diskHoppingEvent)
{
    copy_from(diskHoppingEvent);
}


DiskHoppingEvent::~DiskHoppingEvent()
{
}


DiskHoppingEvent & DiskHoppingEvent::operator=(const DiskHoppingEvent & diskHoppingEvent)
{
    if (this != &diskHoppingEvent)
    {
        copy_from(diskHoppingEvent);
    }

    return *this;
}


void DiskHoppingEvent::copy_from(const DiskHoppingEvent & diskHoppingEvent)
{
    EventOfDynamicVD2D::copy_from(diskHoppingEvent);

    m_TargetGenerator = diskHoppingEvent.m_TargetGenerator;

    m_InWhatContainerThisGeneratorIsDeletedNInserted[0]
        = diskHoppingEvent.m_InWhatContainerThisGeneratorIsDeletedNInserted[0];

    m_InWhatContainerThisGeneratorIsDeletedNInserted[1] 
        = diskHoppingEvent.m_InWhatContainerThisGeneratorIsDeletedNInserted[1];

    m_DeletedNInsertedPosition[0] 
        = diskHoppingEvent.m_DeletedNInsertedPosition[0];

    m_DeletedNInsertedPosition[1]
        = diskHoppingEvent.m_DeletedNInsertedPosition[1];

    m_doesThisEventHappen = diskHoppingEvent.m_doesThisEventHappen;

    m_VelocityVectorAfterCollisionOfTargetGenerator 
        = diskHoppingEvent.m_VelocityVectorAfterCollisionOfTargetGenerator;
  
    m_VelocityVectorBeforeCollisionOfTargetGenerator
        = diskHoppingEvent.m_VelocityVectorBeforeCollisionOfTargetGenerator;

    m_VelocityVectorSetting = diskHoppingEvent.m_VelocityVectorSetting;
}
