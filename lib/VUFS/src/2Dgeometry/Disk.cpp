#include "Disk.h"

Disk::Disk()
{
	m_id = -1;
}

Disk::Disk(const rg_INT& id)
{ 
	m_id = id;
}

Disk::Disk(const rg_Point2D& center, const rg_REAL& radius) : rg_Circle2D(center, radius)
{
	m_id = -1;
}

Disk::Disk(const rg_REAL& xCoord, const rg_REAL& yCoord, const rg_REAL& radius) : rg_Circle2D(xCoord, yCoord, radius)
{
	m_id = -1;
}

Disk::Disk(const rg_Circle2D& circle) : rg_Circle2D(circle)
{
	m_id = -1;
}

Disk::Disk(const rg_INT& id, const rg_Point2D& center, const rg_REAL& radius) : rg_Circle2D(center, radius)
{
	m_id = id;
}

Disk::Disk(const rg_INT& id, const rg_REAL& xCoord, const rg_REAL& yCoord, const rg_REAL& radius) : rg_Circle2D(xCoord, yCoord, radius)
{ 
	m_id = id;
}

Disk::Disk(const rg_INT& id, const rg_Circle2D& circle) : rg_Circle2D(circle)
{
	m_id = id;
}

Disk::Disk(const Disk& disk) : rg_Circle2D(disk)
{
	m_id = disk.m_id;
}

Disk::~Disk()
{
}

Disk& Disk::operator=(const Disk& disk)
{
	if (this == &disk) {
		return *this;
	}

	m_id = disk.m_id;
	setCenterPt(disk.getCenterPt());
	setRadius(disk.getRadius());

	return *this;
}
