#include "DiskCollisionEvent.h"
using namespace V::GeometryTier;

DiskCollisionEvent::DiskCollisionEvent()
{
	m_eventType                    = DISK_COLLISION;
	m_occurringTime                = -1.0;
	m_CollidedDiskGenerator[0]     = NULL;
	m_CollidedDiskGenerator[1]     = NULL;
    m_VelocityVectorSetting        = false;
}


DiskCollisionEvent::DiskCollisionEvent(Generator2D* diskGenerator1, Generator2D* diskGenerator2)
{
    m_eventType                = DISK_COLLISION;
    m_occurringTime            = -1.0;
    m_CollidedDiskGenerator[0] = diskGenerator1;
    m_CollidedDiskGenerator[1] = diskGenerator2;
    m_VelocityVectorSetting    = false;
}


DiskCollisionEvent::DiskCollisionEvent(const double& occuringTime, Generator2D* diskGenerator1, Generator2D* diskGenerator2)
{
	m_eventType                    = DISK_COLLISION;
	m_occurringTime                = occuringTime;
    m_CollidedDiskGenerator[0]     = diskGenerator1;
    m_CollidedDiskGenerator[1]     = diskGenerator2;
    m_VelocityVectorSetting        = false;
}



DiskCollisionEvent::DiskCollisionEvent(const DiskCollisionEvent& diskCollisionEvent)
{
	copy_from(diskCollisionEvent);
}



DiskCollisionEvent::~DiskCollisionEvent()
{
}



void DiskCollisionEvent::copy_from(const DiskCollisionEvent& diskCollisionEvent)
{
	EventOfDynamicVD2D::copy_from(diskCollisionEvent);

	m_CollidedDiskGenerator[0]          = diskCollisionEvent.m_CollidedDiskGenerator[0];
	m_CollidedDiskGenerator[1]          = diskCollisionEvent.m_CollidedDiskGenerator[1];
	m_MotionVectorsAfterCollision[0]    = diskCollisionEvent.m_MotionVectorsAfterCollision[0];
	m_MotionVectorsAfterCollision[1]    = diskCollisionEvent.m_MotionVectorsAfterCollision[1];
	m_MotionVectorsBeforeCollision[0]   = diskCollisionEvent.m_MotionVectorsBeforeCollision[0];
	m_MotionVectorsBeforeCollision[1]   = diskCollisionEvent.m_MotionVectorsBeforeCollision[1];
	m_DiskCentersAtCollisionTime[0]     = diskCollisionEvent.m_DiskCentersAtCollisionTime[0];
	m_DiskCentersAtCollisionTime[1]     = diskCollisionEvent.m_DiskCentersAtCollisionTime[1];
	m_VelocityVectorSetting               = diskCollisionEvent.m_VelocityVectorSetting;

}


DiskCollisionEvent& DiskCollisionEvent::operator=(const DiskCollisionEvent& diskCollisionEvent)
{
	if (this != &diskCollisionEvent)
	{
		copy_from(diskCollisionEvent);
	}

   
	return *this;
}


