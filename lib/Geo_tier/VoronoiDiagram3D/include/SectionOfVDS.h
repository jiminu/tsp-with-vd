#ifndef _SECTIONOFVDS_H
#define _SECTIONOFVDS_H



#include "rg_SphereSetVoronoiDiagram.h"

#include "VDCell.h"
#include "VDFace.h"
#include "VDEdge.h"

#include "rg_Point3D.h"
#include "rg_dList.h"

//#include "Pair.h"
#include "rg_WingedEdgeDataStructure.h"


#include "GateForSectionOfVDS.h"



namespace V {

namespace GeometryTier {


class SectionOfVDS
{
private:
    rg_dList<VertexOnSectionOfVDS> m_vertexList;
    rg_dList<EdgeOnSectionOfVDS>   m_edgeList;
    rg_dList<FaceOnSectionOfVDS>   m_faceList;

    rg_SphereSetVoronoiDiagram*    m_VDS;
    rg_REAL m_plane[4];


public:
    SectionOfVDS();
    SectionOfVDS( const SectionOfVDS& section);
    ~SectionOfVDS();

    rg_dList<VertexOnSectionOfVDS>* getVertexList();
    rg_dList<EdgeOnSectionOfVDS>*   getEdgeList();
    rg_dList<FaceOnSectionOfVDS>*   getFaceList();

    void connectVDS( rg_SphereSetVoronoiDiagram* aVDS );

    void setPlane(rg_REAL* plane);

    void addVertex( const VertexOnSectionOfVDS& vertex );
    void addEdge( const EdgeOnSectionOfVDS& edge );
    void addFace( const FaceOnSectionOfVDS& face );

    void computeSectionOfVDS(rg_SphereSetVoronoiDiagram* aVDS, rg_REAL* plane);
        void    initializeComputingSectionOfVDS(rg_dList< GateForSectionOfVDS >& gateList);
        VDEdge* findNewVertex(GateForSectionOfVDS* gate, rg_FLAG& vertexDeterminant);
        rg_FLAG exploreNewVertex( VDEdge*                          geoForNewVertex,
                                  const rg_FLAG&                   vertexDeterminant,
                                  GateForSectionOfVDS*             gate,
                                  rg_dList< GateForSectionOfVDS >& gateList);
    void computeVertexCoordnates();
    void adjustFacesIncidentToInfiniteVertex();
    void checkTopology();

};

} // namespace GeometryTier

} // namespace V


#endif

