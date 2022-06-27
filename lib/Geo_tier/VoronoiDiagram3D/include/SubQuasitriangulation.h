#ifndef _SUBQUASITRIANGULATION_H
#define _SUBQUASITRIANGULATION_H

#include "BetaUniverse.h"
#include "ClusterFromQT.h"
#include "rg_dList.h"
#include "rg_Molecule.h"

namespace V {

namespace GeometryTier {


class SubQuasitriangulation  
{
private:
    BetaUniverse*           m_betaUniverse;
    rg_dList<ClusterFromQT> m_clusters;    
    
public:
    SubQuasitriangulation();
    SubQuasitriangulation(BetaUniverse* betaUniverse);
    ~SubQuasitriangulation();
    
    inline  rg_dList<ClusterFromQT>*  getClusters()      {return &m_clusters;}
    inline  BetaUniverse*             getBetaUniverse()  {return m_betaUniverse;}
    ClusterType                       getClusterType(const rg_INT& index) const;

    inline  void setBetaUniverse(BetaUniverse* betaUniverse) {m_betaUniverse = betaUniverse;}


    void  makeClustersForSecondaryStructures(Molecule* aMolecule);

    void  huntCellsInCluster4GivenBetaValue(    const rg_REAL& beta, rg_dList< rg_dList< BetaCell*> >&   betaCellList ) const;
    void  huntFacesInCluster4GivenBetaValue(    const rg_REAL& beta, rg_dList< rg_dList< BetaFace*> >&   betaFaceList ) const;
    void  huntEdgesInCluster4GivenBetaValue(    const rg_REAL& beta, rg_dList< rg_dList< BetaEdge*> >&   betaEdgeList ) const;
    void  huntVerticesInCluster4GivenBetaValue( const rg_REAL& beta, rg_dList< rg_dList< BetaVertex*> >& betaVertexList ) const;
        
    void  huntFacesOnBoundaryOfCluster4GivenBetaValue(    const rg_REAL& beta, rg_dList< rg_dList< BetaFace*> >&   betaFaceList ) const;
    void  huntEdgesOnBoundaryOfCluster4GivenBetaValue(    const rg_REAL& beta, rg_dList< rg_dList< BetaEdge*> >&   betaEdgeList ) const;
    void  huntVerticesOnBoundaryOfCluster4GivenBetaValue( const rg_REAL& beta, rg_dList< rg_dList< BetaVertex*> >& betaVertexList ) const;


private:
    void  huntCellsInCluster4GivenBetaValue(    const rg_REAL& beta, rg_dList<ClusterFromQT>&   clusterList ) const;
};

} // namespace GeometryTier

} // namespace V


#endif
