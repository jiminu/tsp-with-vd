#include "BucketCellIndex.h"
using namespace V::GeometryTier;

BucketCellIndex::BucketCellIndex()
: m_i(-1), m_j(-1), m_k(-1)
{}



BucketCellIndex::BucketCellIndex(const rg_INT& i, const rg_INT& j, const rg_INT& k)
: m_i(i), m_j(j), m_k(k)
{}



BucketCellIndex::BucketCellIndex(const BucketCellIndex& cellIndex)
: m_i(cellIndex.m_i), m_j(cellIndex.m_j), m_k(cellIndex.m_k)
{}



BucketCellIndex::~BucketCellIndex()
{}



void   BucketCellIndex::setI(const rg_INT& i)  
{ 
    m_i = i; 
}



void   BucketCellIndex::setJ(const rg_INT& j)  
{ 
    m_j = j; 
}



void   BucketCellIndex::setK(const rg_INT& k)  
{ 
    m_k = k; 
}



void   BucketCellIndex::setCellIndex(const rg_INT& i, const rg_INT& j, const rg_INT& k)
{
    m_i = i; 
    m_j = j; 
    m_k = k; 
}



BucketCellIndex& BucketCellIndex::operator =(const BucketCellIndex& cellIndex)
{
    if ( this == &cellIndex )
        return *this;

    m_i = cellIndex.m_i;
    m_j = cellIndex.m_j;
    m_k = cellIndex.m_k;

    return *this;
}

