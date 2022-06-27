#include "DiskVelocityChangeEvent.h"
using namespace V::GeometryTier;

DiskVelocityChangeEvent::DiskVelocityChangeEvent()
{
    m_eventType                = DISK_VELOCITY_CHANGE;
    m_occurringTime            = -1.0;
    m_DiskGenerator            = NULL;
    m_VelocityVectorSetting    = false;
}


DiskVelocityChangeEvent::DiskVelocityChangeEvent(Generator2D* diskGenerator)
{
    m_eventType                = DISK_VELOCITY_CHANGE;
    m_occurringTime            = -1.0;
    m_DiskGenerator            = diskGenerator;
    m_VelocityVectorSetting    = false;
}


DiskVelocityChangeEvent::DiskVelocityChangeEvent(const double& occuringTime, Generator2D* diskGenerator)
{
    m_eventType                = DISK_VELOCITY_CHANGE;
    m_occurringTime            = occuringTime;
    m_DiskGenerator            = diskGenerator;
    m_VelocityVectorSetting    = false;
}


DiskVelocityChangeEvent::DiskVelocityChangeEvent(const DiskVelocityChangeEvent& diskVelocityChangeEvent)
{
    copy_from(diskVelocityChangeEvent);
}


DiskVelocityChangeEvent::DiskVelocityChangeEvent(const double& occuringTime, Generator2D* diskGenerator, const rg_Point2D& newVelocityVector)
    :m_VelocityVectorAfterVelocityChange(newVelocityVector)
{
    m_eventType                = DISK_VELOCITY_CHANGE;
    m_occurringTime            = occuringTime;
    m_DiskGenerator            = diskGenerator;
    m_VelocityVectorSetting    = false;
}


DiskVelocityChangeEvent::~DiskVelocityChangeEvent()
{
}


void DiskVelocityChangeEvent::copy_from(const DiskVelocityChangeEvent& diskVelocityChangeEvent)
{
    EventOfDynamicVD2D::copy_from(diskVelocityChangeEvent);

    m_DiskGenerator                       = diskVelocityChangeEvent.m_DiskGenerator;
    m_VelocityVectorAfterVelocityChange   = diskVelocityChangeEvent.m_VelocityVectorAfterVelocityChange;
    m_VelocityVectorBeforeVelocityChange  = diskVelocityChangeEvent.m_VelocityVectorBeforeVelocityChange;
    m_DiskCenterAtVelocityChangeTime      = diskVelocityChangeEvent.m_DiskCenterAtVelocityChangeTime;
    m_VelocityVectorSetting               = diskVelocityChangeEvent.m_VelocityVectorSetting;
}

DiskVelocityChangeEvent& DiskVelocityChangeEvent::operator=(const DiskVelocityChangeEvent& diskVelocityChangeEvent)
{
    if (this != &diskVelocityChangeEvent)
    {
        copy_from(diskVelocityChangeEvent);
    }
    
    return *this;
}

