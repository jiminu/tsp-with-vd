//#include "StdAfx.h"
#include "NeighborInfo.h"
using namespace V::GeometryTier;


NeighborInfo::NeighborInfo(void)
{
    m_ID = -1;
    m_totalDist = 0.0;
    m_minDist   = DBL_MAX;
    m_maxDist   = -1.0;
}


NeighborInfo::~NeighborInfo(void)
{
}
