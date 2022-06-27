#include "GateForSectionOfVDS.h"
using namespace V::GeometryTier;

GateForSectionOfVDS::GateForSectionOfVDS()
: m_geoForLeftFace(  rg_NULL),
  m_geoForRightFace( rg_NULL),
  m_gateForLeftHand( rg_NULL),
  m_gateForRightHand(rg_NULL),
  m_gateForLeftLeg(  rg_NULL),
  m_gateForRightLeg( rg_NULL),
  m_startVertex(     rg_NULL),
  m_geoEdge(         rg_NULL),
  m_edge(            rg_NULL)
{
}



GateForSectionOfVDS::GateForSectionOfVDS( VDCell* geoForLeftFace,
                                          VDCell* geoForRightFace,
                                          VertexOnSectionOfVDS* startVertex,
                                          VDPartialEdge* geoEdge)
: m_geoForLeftFace(  geoForLeftFace),
  m_geoForRightFace( geoForRightFace),
  m_gateForLeftHand( rg_NULL),
  m_gateForRightHand(rg_NULL),
  m_gateForLeftLeg(  rg_NULL),
  m_gateForRightLeg( rg_NULL),
  m_startVertex(     startVertex),
  m_geoEdge(         geoEdge),
  m_edge(            rg_NULL)
{
}



GateForSectionOfVDS::GateForSectionOfVDS(const GateForSectionOfVDS& gate)
: m_geoForLeftFace(  gate.m_geoForLeftFace),
  m_geoForRightFace( gate.m_geoForRightFace),
  m_gateForLeftHand( gate.m_gateForLeftHand),
  m_gateForRightHand(gate.m_gateForRightHand),
  m_gateForLeftLeg(  gate.m_gateForLeftLeg),
  m_gateForRightLeg( gate.m_gateForRightLeg),
  m_startVertex(     gate.m_startVertex),
  m_geoEdge(         gate.m_geoEdge),
  m_edge(            gate.m_edge)
{
}



GateForSectionOfVDS::~GateForSectionOfVDS()
{
    m_geoForLeftFace   = rg_NULL;
    m_geoForRightFace  = rg_NULL;
    m_gateForLeftHand  = rg_NULL;
    m_gateForRightHand = rg_NULL;
    m_gateForLeftLeg   = rg_NULL;
    m_gateForRightLeg  = rg_NULL;
    m_startVertex      = rg_NULL;
    m_geoEdge          = rg_NULL;
    m_edge             = rg_NULL;
}




VDCell* GateForSectionOfVDS::getGeometryForLeftFace() const
{
    return m_geoForLeftFace;
}



VDCell* GateForSectionOfVDS::getGeometryForRightFace() const
{
    return m_geoForRightFace;
}



GateForSectionOfVDS*  GateForSectionOfVDS::getGateForLeftHand() const
{
    return m_gateForLeftHand;
}



GateForSectionOfVDS*  GateForSectionOfVDS::getGateForRightHand() const
{
    return m_gateForRightHand;
}



GateForSectionOfVDS*  GateForSectionOfVDS::getGateForLeftLeg() const
{
    return m_gateForLeftLeg;
}



GateForSectionOfVDS*  GateForSectionOfVDS::getGateForRightLeg() const
{
    return m_gateForRightLeg;
}



VertexOnSectionOfVDS* GateForSectionOfVDS::getStartVertex() const
{
    return m_startVertex;
}



VDPartialEdge* GateForSectionOfVDS::getGeometryForEdge() const
{
    return m_geoEdge;
}



EdgeOnSectionOfVDS*   GateForSectionOfVDS::getEdge() const
{
    return m_edge;
}




void GateForSectionOfVDS::setGeometryForLeftFace(VDCell* geoForLeftFace)
{
    m_geoForLeftFace = geoForLeftFace;
}



void GateForSectionOfVDS::setGeometryForRightFace(VDCell* geoForRightFace)
{
    m_geoForRightFace = geoForRightFace;
}



void GateForSectionOfVDS::setGateForLeftHand(GateForSectionOfVDS* gateForLeftHand)
{
    m_gateForLeftHand = gateForLeftHand;
}



void GateForSectionOfVDS::setGateForRightHand(GateForSectionOfVDS* gateForRightHand)
{
    m_gateForRightHand = gateForRightHand;
}



void GateForSectionOfVDS::setGateForLeftLeg(GateForSectionOfVDS* gateForLeftLeg)
{
    m_gateForLeftLeg = gateForLeftLeg;
}



void GateForSectionOfVDS::setGateForRightLeg(GateForSectionOfVDS* gateForRightLeg)
{
    m_gateForRightLeg = gateForRightLeg;
}



void GateForSectionOfVDS::setStartVertex(VertexOnSectionOfVDS* startVertex)
{
    m_startVertex = startVertex;
}



void GateForSectionOfVDS::setGeometryForEdge(VDPartialEdge* geoEdge)
{
    m_geoEdge = geoEdge;
}



void GateForSectionOfVDS::setEdge(EdgeOnSectionOfVDS* edge)
{
    m_edge = edge;
}




GateForSectionOfVDS& GateForSectionOfVDS::operator =(const GateForSectionOfVDS& gate)
{
    if ( this == &gate )
        return *this;

    m_geoForLeftFace   = gate.m_geoForLeftFace;
    m_geoForRightFace  = gate.m_geoForRightFace;
    m_gateForLeftHand  = gate.m_gateForLeftHand;
    m_gateForRightHand = gate.m_gateForRightHand;
    m_gateForLeftLeg   = gate.m_gateForLeftLeg;
    m_gateForRightLeg  = gate.m_gateForRightLeg;
    m_startVertex      = gate.m_startVertex;
    m_geoEdge          = gate.m_geoEdge;
    m_edge             = gate.m_edge;

    return *this;
}



