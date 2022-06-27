#ifndef _GATEFORSECTIONOFVDS_H
#define _GATEFORSECTIONOFVDS_H

#include "rg_Const.h"
#include "VDCell.h"
#include "VDLoop.h"
#include "VDPartialEdge.h"
#include <utility>
using namespace std;
//#include "Pair.h"
#include "rg_WingedEdgeDataStructure.h"

namespace V {

namespace GeometryTier {


typedef WVertex< pair<VDEdge*, rg_Point3D>, VDLoop*, VDCell*> VertexOnSectionOfVDS;
typedef WEdge< pair<VDEdge*, rg_Point3D>, VDLoop*, VDCell*>   EdgeOnSectionOfVDS;
typedef WFace< pair<VDEdge*, rg_Point3D>, VDLoop*, VDCell*>   FaceOnSectionOfVDS;



class GateForSectionOfVDS
{
private:
    VDCell* m_geoForLeftFace;
    VDCell* m_geoForRightFace;

    GateForSectionOfVDS* m_gateForLeftHand;
    GateForSectionOfVDS* m_gateForRightHand;
    GateForSectionOfVDS* m_gateForLeftLeg;
    GateForSectionOfVDS* m_gateForRightLeg;

    VertexOnSectionOfVDS* m_startVertex;

    VDPartialEdge*      m_geoEdge;
    EdgeOnSectionOfVDS* m_edge;  

public:
    rg_INT m_ID;

    GateForSectionOfVDS();
    GateForSectionOfVDS( VDCell* geoForLeftFace,
                         VDCell* geoForRightFace,
                         VertexOnSectionOfVDS* startVertex,
                         VDPartialEdge* geoEdge);
    GateForSectionOfVDS(const GateForSectionOfVDS& gate);
    ~GateForSectionOfVDS();

    VDCell*               getGeometryForLeftFace() const;
    VDCell*               getGeometryForRightFace() const;
    GateForSectionOfVDS*  getGateForLeftHand() const;
    GateForSectionOfVDS*  getGateForRightHand() const;
    GateForSectionOfVDS*  getGateForLeftLeg() const;
    GateForSectionOfVDS*  getGateForRightLeg() const;
    VertexOnSectionOfVDS* getStartVertex() const;
    VDPartialEdge*        getGeometryForEdge() const;
    EdgeOnSectionOfVDS*   getEdge() const;

    void setGeometryForLeftFace(VDCell* geoForLeftFace);
    void setGeometryForRightFace(VDCell* geoForRightFace);
    void setGateForLeftHand(GateForSectionOfVDS* gateForLeftHand);
    void setGateForRightHand(GateForSectionOfVDS* gateForRightHand);
    void setGateForLeftLeg(GateForSectionOfVDS* gateForLeftLeg);
    void setGateForRightLeg(GateForSectionOfVDS* gateForRightLeg);
    void setStartVertex(VertexOnSectionOfVDS* startVertex);
    void setGeometryForEdge(VDPartialEdge* geoEdge);
    void setEdge(EdgeOnSectionOfVDS* edge);

    GateForSectionOfVDS& operator =(const GateForSectionOfVDS& gate);
};

} // namespace GeometryTier

} // namespace V


#endif

