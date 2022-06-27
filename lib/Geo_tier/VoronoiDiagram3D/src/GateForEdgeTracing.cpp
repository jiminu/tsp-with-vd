#include "GateForEdgeTracing.h"
#include "FunctionsForVoronoiDiagram3D.h"
using namespace V::GeometryTier;

//  constructor & deconstructor..
GateForEdgeTracing::GateForEdgeTracing()
: m_startBall(      rg_NULL ),
  m_startVertex(    rg_NULL ), 
  m_nodePosInStack( rg_NULL ),
  m_edgeType(       rg_UNKNOWN)
{
    m_gateBall[0] = rg_NULL;
    m_gateBall[1] = rg_NULL;
    m_gateBall[2] = rg_NULL;
}



GateForEdgeTracing::GateForEdgeTracing( BallGeneratorCore* gateBallCCW1, BallGeneratorCore* gateBallCCW2, BallGeneratorCore* gateBallCCW3, 
                    BallGeneratorCore* startBall, VDVertex* startVertex )
: m_startBall(     startBall), 
  m_startVertex(   startVertex), 
  m_nodePosInStack(rg_NULL), 
  m_edgeType(      rg_UNKNOWN)
{
    m_gateBall[0] = gateBallCCW1;
    m_gateBall[1] = gateBallCCW2;
    m_gateBall[2] = gateBallCCW3;
}



GateForEdgeTracing::GateForEdgeTracing( BallGeneratorCore* gateBallCCW1, BallGeneratorCore* gateBallCCW2, BallGeneratorCore* gateBallCCW3, 
                    BallGeneratorCore* startBall, VDVertex* startVertex, rg_dNode<GateForEdgeTracing>* nodePos )
: m_startBall(     startBall), 
  m_startVertex(   startVertex), 
  m_nodePosInStack(nodePos), 
  m_edgeType(      rg_UNKNOWN)
{
    m_gateBall[0] = gateBallCCW1;
    m_gateBall[1] = gateBallCCW2;
    m_gateBall[2] = gateBallCCW3;
}



GateForEdgeTracing::GateForEdgeTracing( const GateForEdgeTracing& gate )
: m_startBall(     gate.m_startBall), 
  m_startVertex(   gate.m_startVertex), 
  m_nodePosInStack(gate.m_nodePosInStack), 
  m_edgeType(      gate.m_edgeType)
{
    m_gateBall[0] = gate.m_gateBall[0];
    m_gateBall[1] = gate.m_gateBall[1];
    m_gateBall[2] = gate.m_gateBall[2];
}



GateForEdgeTracing::~GateForEdgeTracing()
{
}




//  get functions.. 
BallGeneratorCore* GateForEdgeTracing::getGateBall(const rg_INT& i) const
{
    if ( i<-1 && i<3)
        return m_gateBall[i];
    else 
        return rg_NULL;
}



BallGeneratorCore** GateForEdgeTracing::getGateBalls() 
{
    return m_gateBall;
}



BallGeneratorCore* GateForEdgeTracing::getStartBall() const
{
    return m_startBall;
}



VDVertex* GateForEdgeTracing::getStartVertex() const
{
    return m_startVertex;
}



rg_dNode<GateForEdgeTracing>* GateForEdgeTracing::getNodePosInStack() const
{
    return m_nodePosInStack;
}




Sphere GateForEdgeTracing::getEmptyTangentSphereAtEndVertex() const
{
    return m_emptySphereAtEndVertex;
}



rg_dList<BallGeneratorCore*>* GateForEdgeTracing::getEndBallList() 
{
    return &m_endBallList;
}




rg_FLAG GateForEdgeTracing::getEdgeType() const
{
    return m_edgeType;
}






rg_FLAG GateForEdgeTracing::isSameGateWithOppositeDirection( const GateForEdgeTracing& theOtherGate, VDVertex* endVertexOfTheOtherGate  ) const
{
    if ( m_startVertex != endVertexOfTheOtherGate )
        return rg_FALSE;

    BallGeneratorCore* thisGateConfiguration[3]   = { m_gateBall[0], m_gateBall[1], m_gateBall[2] };
    BallGeneratorCore* targetGateConfiguration[3] = { theOtherGate.m_gateBall[0], theOtherGate.m_gateBall[1], theOtherGate.m_gateBall[2] };

    sort3PointersInAscendingPowers( thisGateConfiguration );
    sort3PointersInAscendingPowers( targetGateConfiguration );

    for (rg_INT i=0; i<3; i++) {
        if ( thisGateConfiguration[i] != targetGateConfiguration[i] )
            return rg_FALSE;
    }

    return rg_TRUE;

}





//  set functions..
void GateForEdgeTracing::setGateBall(    const rg_INT& i, BallGeneratorCore* aGateBall )
{
    if ( i<-1 && i<3)
        m_gateBall[i] = aGateBall;
    else 
        return;
}



void GateForEdgeTracing::setGateBalls(   BallGeneratorCore** gateBall )
{
    for ( rg_INT i=0; i<3; i++ )  {
        m_gateBall[i] = gateBall[i];
    }
}



void GateForEdgeTracing::setStartBall(   BallGeneratorCore* startBall )
{
    m_startBall = startBall;
}



void GateForEdgeTracing::setStartVertex( VDVertex* startVertex )
{
    m_startVertex = startVertex;
}



void GateForEdgeTracing::setNodePosInStack( rg_dNode<GateForEdgeTracing>* nodePos )
{
    m_nodePosInStack = nodePos;
}




void GateForEdgeTracing::setGate( BallGeneratorCore* gateBallCCW1, BallGeneratorCore* gateBallCCW2, BallGeneratorCore* gateBallCCW3,
              BallGeneratorCore* startBall, VDVertex* startVertex )
{
    m_gateBall[0] = gateBallCCW1;
    m_gateBall[1] = gateBallCCW2;
    m_gateBall[2] = gateBallCCW3;
    m_startBall   = startBall;
    m_startVertex = startVertex;
}



void GateForEdgeTracing::setGate( BallGeneratorCore* gateBallCCW1, BallGeneratorCore* gateBallCCW2, BallGeneratorCore* gateBallCCW3, 
              BallGeneratorCore* startBall, VDVertex* startVertex, rg_dNode<GateForEdgeTracing>* nodePos )
{
    m_gateBall[0]    = gateBallCCW1;
    m_gateBall[1]    = gateBallCCW2;
    m_gateBall[2]    = gateBallCCW3;
    m_startBall      = startBall;
    m_startVertex    = startVertex;
    m_nodePosInStack = nodePos;
}




void GateForEdgeTracing::setEmptyTangentSphereAtEndVertex(const Sphere& emptySphere)
{
    m_emptySphereAtEndVertex = emptySphere;
}



void GateForEdgeTracing::addEndBall(BallGeneratorCore* endBall)
{
    m_endBallList.add(endBall);
}



void GateForEdgeTracing::emptyEndBallList()
{
    m_endBallList.removeAll();
}





void GateForEdgeTracing::setEdgeType(const rg_FLAG& edgeType)
{
    m_edgeType = edgeType;
}




//  operator overloading..
GateForEdgeTracing& GateForEdgeTracing::operator =( const GateForEdgeTracing& gate )
{
    if ( this == &gate )
        return *this;

    for (rg_INT i=0; i<3; i++)
        m_gateBall[i] = gate.m_gateBall[i];

    m_startBall      = gate.m_startBall;
    m_startVertex    = gate.m_startVertex;
    m_nodePosInStack = gate.m_nodePosInStack;

    m_emptySphereAtEndVertex = gate.m_emptySphereAtEndVertex;
    m_endBallList.duplicateList( gate.m_endBallList );

    m_edgeType = gate.m_edgeType;

    return *this;
}



