// lusterFromQT.cpp: implementation of the ClusterFromQT class.
//
//////////////////////////////////////////////////////////////////////
#include "ClusterFromQT.h"
using namespace V::GeometryTier;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ClusterFromQT::ClusterFromQT()
{
    m_type = UNK;
}

ClusterFromQT::ClusterFromQT( const ClusterType& typeOfCluster )
{
    m_type = typeOfCluster;
}

ClusterFromQT::~ClusterFromQT()
{

}

rg_BOOL ClusterFromQT::isBetaCellInCluster( BetaCell* aBetaCell )
{
    rg_BOOL isCellInCluster = rg_TRUE;

    Hash4BetaCell::iterator betaCell_i = m_betaCells.find( aBetaCell->getID() );
    
    if ( betaCell_i == m_betaCells.end() ) {
        isCellInCluster = rg_FALSE;
    }

    return isCellInCluster;
}