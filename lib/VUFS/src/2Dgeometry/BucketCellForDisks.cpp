#include "BucketCellForDisks.h"

BucketCellForDisks::BucketCellForDisks()
{
	m_Xindex = -1;
	m_Yindex = -1;
}


BucketCellForDisks::BucketCellForDisks(const BucketCellForDisks& bucketCellForDisk)
{
	m_Xindex		  = bucketCellForDisk.m_Xindex;
	m_Yindex		  = bucketCellForDisk.m_Yindex;
	m_ContainedDisks  = bucketCellForDisk.m_ContainedDisks;

}


BucketCellForDisks::~BucketCellForDisks()
{
}


BucketCellForDisks& BucketCellForDisks::operator=(const BucketCellForDisks& bucketCellForDisk)
{
	if (this != &bucketCellForDisk)
	{
		m_Xindex		    = bucketCellForDisk.m_Xindex;
		m_Yindex		    = bucketCellForDisk.m_Yindex;
		m_ContainedDisks	= bucketCellForDisk.m_ContainedDisks;
	}

	return *this;
}

