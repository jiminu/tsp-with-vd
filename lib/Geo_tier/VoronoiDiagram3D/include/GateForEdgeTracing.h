#ifndef _GATEFOREDGETRACING_H
#define _GATEFOREDGETRACING_H

#include "ConstForVoronoiDiagram3D.h"
#include "rg_dList.h"
#include "Sphere.h"

namespace V {

namespace GeometryTier {


class BallGeneratorCore;
class VDVertex;



class GateForEdgeTracing
{
private:
    BallGeneratorCore*   m_gateBall[3];
    BallGeneratorCore*   m_startBall;
    VDVertex*        m_startVertex;
    rg_dNode<GateForEdgeTracing>*  m_nodePosInStack;

	Sphere           m_emptySphereAtEndVertex;
    rg_dList<BallGeneratorCore*> m_endBallList;

	rg_FLAG          m_edgeType;

public:
    //  constructor & deconstructor..
    GateForEdgeTracing();
    GateForEdgeTracing( BallGeneratorCore* gateBallCCW1, BallGeneratorCore* gateBallCCW2, BallGeneratorCore* gateBallCCW3, 
                        BallGeneratorCore* startBall, VDVertex* startVertex );
    GateForEdgeTracing( BallGeneratorCore* gateBallCCW1, BallGeneratorCore* gateBallCCW2, BallGeneratorCore* gateBallCCW3, 
                        BallGeneratorCore* startBall, VDVertex* startVertex, rg_dNode<GateForEdgeTracing>* nodePos );
    GateForEdgeTracing( const GateForEdgeTracing& gate );
    ~GateForEdgeTracing();

    //  get functions.. 
    BallGeneratorCore*   getGateBall(const rg_INT& i) const;
    BallGeneratorCore**  getGateBalls();
    BallGeneratorCore*   getStartBall() const;
    VDVertex*        getStartVertex() const;  
    rg_dNode<GateForEdgeTracing>* getNodePosInStack() const;

	Sphere           getEmptyTangentSphereAtEndVertex() const;
	rg_dList<BallGeneratorCore*>* getEndBallList();

    rg_FLAG          getEdgeType() const;


    
    rg_FLAG isSameGateWithOppositeDirection( const GateForEdgeTracing& theOtherGate, VDVertex* endVertexOfTheOtherGate  ) const;


    //  set functions..
    void setGateBall(    const rg_INT& i, BallGeneratorCore* aGateBall );
    void setGateBalls(   BallGeneratorCore** gateBall );
    void setStartBall(   BallGeneratorCore* startBall );
    void setStartVertex( VDVertex* startVertex );
    void setNodePosInStack( rg_dNode<GateForEdgeTracing>* nodePos );

    void setGate( BallGeneratorCore* gateBallCCW1, BallGeneratorCore* gateBallCCW2, BallGeneratorCore* gateBallCCW3, 
                  BallGeneratorCore* startBall, VDVertex* startVertex );
    void setGate( BallGeneratorCore* gateBallCCW1, BallGeneratorCore* gateBallCCW2, BallGeneratorCore* gateBallCCW3, 
                  BallGeneratorCore* startBall, VDVertex* startVertex, rg_dNode<GateForEdgeTracing>* nodePos );

	void setEmptyTangentSphereAtEndVertex(const Sphere& emptySphere);
	void addEndBall(BallGeneratorCore* endBall);
    void emptyEndBallList();

    
	void setEdgeType(const rg_FLAG& edgeType);

    //  operator overloading..
    GateForEdgeTracing& operator =( const GateForEdgeTracing& gate );

};

} // namespace GeometryTier

} // namespace V


#endif 


