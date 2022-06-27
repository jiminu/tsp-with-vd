#include "BucketCellForDisk.h"
using namespace V::GeometryTier;

BucketCellForDisk::BucketCellForDisk(void)
{
    m_containingGenerator   = NULL;
    //m_isChecked = false;
}



BucketCellForDisk::BucketCellForDisk( Generator2D* const containingGenerator )
{
    m_containingGenerator   = containingGenerator;
    //m_isChecked         = false;
}



BucketCellForDisk::BucketCellForDisk( const BucketCellForDisk& bucketCell )
{
    m_containingGenerator   = bucketCell.m_containingGenerator;
    //m_isChecked         = bucketCell.m_isChecked;
}



BucketCellForDisk::~BucketCellForDisk(void)
{
}


