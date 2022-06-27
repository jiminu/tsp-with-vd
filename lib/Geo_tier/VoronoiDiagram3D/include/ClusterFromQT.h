#ifndef _CLUSTERFROMQT_H
#define _CLUSTERFROMQT_H

#include "BetaCell.h"

#include <map>
using namespace std;

namespace V {

namespace GeometryTier {


typedef map<int, BetaCell*> Hash4BetaCell;

enum ClusterType {UNK, HELIX, BETASTRAND};

class ClusterFromQT  
{
private:
    Hash4BetaCell   m_betaCells;
    ClusterType     m_type;
    
public:
    ClusterFromQT();
    ClusterFromQT(const ClusterType& typeOfCluster);
    ~ClusterFromQT();
    
    inline  Hash4BetaCell*  getBetaCells()      {return &m_betaCells;}
    inline  ClusterType     getTypeOfCluster()  {return m_type;}


    inline  void            setTyperOfCluster(const ClusterType& typeOfCluster) {m_type = typeOfCluster;}
    inline  void            addCellIntoHash(BetaCell* aBetaCell)                {m_betaCells.insert( Hash4BetaCell::value_type( aBetaCell->getID(), aBetaCell ) );}


    rg_BOOL isBetaCellInCluster(BetaCell* aBetaCell);
};

} // namespace GeometryTier

} // namespace V


#endif
