#include "NavigatingBucketCell.h"
using namespace V::GeometryTier;


NavigatingBucketCell::NavigatingBucketCell()
: m_index(), m_status(rg_UNKNOWN), m_distBetPandVi(0.), m_prevCell(rg_UNKNOWN)  
{
    for (rg_INT i=0; i<4; i++)
        m_knownFaceMask[i] = 0;
}



NavigatingBucketCell::NavigatingBucketCell(const BucketCellIndex& index, 
                                           const rg_FLAG&         status)
: m_index(index), m_status(status)
{}




NavigatingBucketCell::NavigatingBucketCell(
                     const BucketCellIndex& index, 
                     const rg_FLAG&         status, 
                     const rg_REAL&         dist)
: m_index(index), m_status(status), m_distBetPandVi(dist)  
{
}



NavigatingBucketCell::NavigatingBucketCell(
                     const BucketCellIndex& index, 
                     const rg_FLAG&         status, 
                     const rg_REAL&         dist,
                     const rg_FLAG&         prevCell,
                           rg_INT*          knownFaceMask)
: m_index(index), m_status(status), m_distBetPandVi(dist), m_prevCell(prevCell)  
{
    for (rg_INT i=0; i<4; i++)
        m_knownFaceMask[i] = knownFaceMask[i];
}



NavigatingBucketCell::NavigatingBucketCell(const NavigatingBucketCell& naviCell)
: m_index(naviCell.m_index), m_status(naviCell.m_status), 
  m_distBetPandVi(naviCell.m_distBetPandVi), m_prevCell(naviCell.m_prevCell)  
{
    for (rg_INT i=0; i<4; i++)
        m_knownFaceMask[i] = naviCell.m_knownFaceMask[i];
}



NavigatingBucketCell::~NavigatingBucketCell()
{
}




BucketCellIndex NavigatingBucketCell::getIndex() const
{
    return m_index;
}




rg_FLAG NavigatingBucketCell::getStatus() const
{
    return m_status;
}



rg_REAL NavigatingBucketCell::getDistanceVi() const
{
    return m_distBetPandVi;
}




rg_FLAG NavigatingBucketCell::getPrevCell() const
{
    return m_prevCell;
}



rg_INT  NavigatingBucketCell::getKnownFaceMask(const rg_INT& i) const
{
    return m_knownFaceMask[i];
}





void    NavigatingBucketCell::setIndex(const BucketCellIndex& index)
{
    m_index = index;
}



void    NavigatingBucketCell::setStatus(const rg_FLAG& status)
{
    m_status = status;
}



void    NavigatingBucketCell::setDistanceVi(const rg_REAL& dist)
{
    m_distBetPandVi = dist;
}




void    NavigatingBucketCell::setPrevCell(const rg_FLAG& prevCell)
{
    m_prevCell = prevCell;
}



void    NavigatingBucketCell::setKnownFaceMask(rg_INT* knownFaceMask)
{
    for (rg_INT i=0; i<4; i++)
        m_knownFaceMask[i] = knownFaceMask[i];
}




void    NavigatingBucketCell::setCell(const BucketCellIndex& index, 
                                      const rg_FLAG&         status, 
                                      const rg_REAL&         dist,
                                      const rg_FLAG&         prevCell,
                                            rg_INT*          knownFaceMask)
{
    m_index         = index;
    m_status        = status;
    m_distBetPandVi = dist;
    m_prevCell      = prevCell;
    
    for (rg_INT i=0; i<4; i++)
        m_knownFaceMask[i] = knownFaceMask[i];
}




NavigatingBucketCell& NavigatingBucketCell::operator =(const NavigatingBucketCell& naviCell)
{
    if ( this == &naviCell )
        return *this;

    m_index         = naviCell.m_index;
    m_status        = naviCell.m_status;
    m_distBetPandVi = naviCell.m_distBetPandVi;
    m_prevCell      = naviCell.m_prevCell;
    
    for (rg_INT i=0; i<4; i++)
        m_knownFaceMask[i] = naviCell.m_knownFaceMask[i];

    return *this;
}



void NavigatingBucketCell::computeVisFromPrevCell(rg_REAL* distVis, const rg_Point3D& unitLength) const
{

    switch ( m_prevCell )  {
        case I_MINUS:
            distVis[0]= m_knownFaceMask[0];
            distVis[1]= m_distBetPandVi + unitLength.getX();
            distVis[2]= m_distBetPandVi + unitLength.getX() + unitLength.getY();
            distVis[3]= m_knownFaceMask[1];
            distVis[4]= m_knownFaceMask[2];
            distVis[5]= distVis[1] + unitLength.getZ();
            distVis[6]= distVis[2] + unitLength.getZ();
            distVis[7]= m_knownFaceMask[3];
            break;

        case J_MINUS:
            distVis[0]= m_knownFaceMask[1];
            distVis[1]= m_knownFaceMask[0];
            distVis[2]= m_distBetPandVi + unitLength.getX() + unitLength.getY();
            distVis[3]= m_distBetPandVi                     + unitLength.getY();
            distVis[4]= m_knownFaceMask[3];
            distVis[5]= m_knownFaceMask[2];
            distVis[6]= distVis[2] + unitLength.getZ();
            distVis[7]= distVis[3] + unitLength.getZ();
            break;

        case K_MINUS:
            distVis[0]= m_knownFaceMask[0];
            distVis[1]= m_knownFaceMask[1];
            distVis[2]= m_knownFaceMask[2];
            distVis[3]= m_knownFaceMask[3];
            distVis[4]= distVis[0] + unitLength.getZ();
            distVis[5]= distVis[1] + unitLength.getZ();
            distVis[6]= distVis[2] + unitLength.getZ();
            distVis[7]= distVis[3] + unitLength.getZ();
            break;

        case I_PLUS:
            distVis[0]= m_distBetPandVi;
            distVis[1]= m_knownFaceMask[0];
            distVis[2]= m_knownFaceMask[1];
            distVis[3]= m_distBetPandVi                     + unitLength.getY();
            distVis[4]= distVis[0] + unitLength.getZ();
            distVis[5]= m_knownFaceMask[2];
            distVis[6]= m_knownFaceMask[3];
            distVis[7]= distVis[3] + unitLength.getZ();
            break;

        case J_PLUS:
            distVis[0]= m_distBetPandVi;
            distVis[1]= m_distBetPandVi + unitLength.getX();
            distVis[2]= m_knownFaceMask[1];
            distVis[3]= m_knownFaceMask[0];
            distVis[4]= distVis[0] + unitLength.getZ();
            distVis[5]= distVis[1] + unitLength.getZ();
            distVis[6]= m_knownFaceMask[3];
            distVis[7]= m_knownFaceMask[2];
            break;

        case K_PLUS:
            distVis[0]= m_distBetPandVi;
            distVis[1]= m_distBetPandVi + unitLength.getX();
            distVis[2]= m_distBetPandVi + unitLength.getX() + unitLength.getY();
            distVis[3]= m_distBetPandVi                     + unitLength.getY();
            distVis[4]= m_knownFaceMask[0];
            distVis[5]= m_knownFaceMask[1];
            distVis[6]= m_knownFaceMask[2];
            distVis[7]= m_knownFaceMask[3];
            break;

        default:
            break;
    }
}



