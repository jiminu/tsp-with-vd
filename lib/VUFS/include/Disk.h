#ifndef _DISK_H
#define _DISK_H

#include "rg_Circle2D.h"

class Disk : public rg_Circle2D
{
private:
	rg_INT m_id;

public:
	Disk();
	Disk(const rg_INT& id);
	Disk(const rg_Point2D& center, const rg_REAL& radius);
	Disk(const rg_REAL& xCoord, const rg_REAL& yCoord, const rg_REAL& radius);
	Disk(const rg_Circle2D& circle);
	Disk(const rg_INT& id, const rg_Point2D& center, const rg_REAL& radius);
	Disk(const rg_INT& id, const rg_REAL& xCoord, const rg_REAL& yCoord, const rg_REAL& radius);
	Disk(const rg_INT& id, const rg_Circle2D& circle);
	Disk(const Disk& disk);
	virtual ~Disk();

	inline rg_INT getID() const
	{
		return m_id;
	}

	inline void setID(const rg_INT& id)
	{
		m_id = id;
	}

	inline void setDisk(const rg_INT& id, const rg_Point2D& center, const rg_REAL& radius)
	{
		m_id = id;
		setCircle(center, radius);
	}

	inline void setDisk(rg_INT& id, rg_Circle2D& disk)
	{
		m_id = id;
		setCircle(disk);
	}

	Disk& operator=(const Disk& disk);
};

#endif