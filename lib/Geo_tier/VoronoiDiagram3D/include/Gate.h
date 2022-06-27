#ifndef _GATE_H
#define _GATE_H

#include "ConstForVoronoiDiagram3D.h"
#include "rg_dList.h"
#include "Sphere.h"
//#include "Pair.h"
#include "Plane.h"

#include <utility>
using namespace std;


namespace V {

namespace GeometryTier {


class BallGenerator;
class VDVertex;

class Gate
{
private:
    BallGenerator*   m_gateBallCCW1;
    BallGenerator*   m_gateBallCCW2;
    BallGenerator*   m_gateBallCCW3;

    BallGenerator*   m_startBall;
    VDVertex*		 m_startVertex;
    rg_dNode<Gate>*  m_nodePosInStack;

	BallGenerator*   m_endBall;

    //rg_dList<BallGenerator*> m_endBallList; // for robustness
    rg_dList< pair<BallGenerator*, Sphere> > m_endBallList; // for robustness
	rg_FLAG                  m_edgeType;
	Sphere                   m_emptySphereAtEndVertex;
    rg_FLAG                  m_numVertexBy4Balls;

    rg_INT          m_numExcludedBalls;
    BallGenerator** m_excludedBall;
    
	//  2006. 8. 17  By Youngsong Cho    
	rg_Point3D  m_normalOfCenterPlane;
	rg_Point3D  m_normalOfEdgePlane;
	rg_Point3D  m_axisPoint;
	rg_Vector3D m_vecAxisToStart;
	Sphere      m_minTangentSphere[2]; 
	rg_REAL     m_minAngularDistance;
	rg_REAL     m_maxValidAngle;

    rg_REAL     m_tolerance;

    Plane       m_centerPlane;
    rg_REAL     m_distanceBetVsAndPc;

public:
    //  constructor & deconstructor..
    Gate();
    Gate( BallGenerator* gateBallCCW1, BallGenerator* gateBallCCW2, BallGenerator* gateBallCCW3, BallGenerator* startBall, VDVertex* startVertex );
    Gate( BallGenerator* gateBallCCW1, BallGenerator* gateBallCCW2, BallGenerator* gateBallCCW3, BallGenerator* startBall, VDVertex* startVertex, rg_dNode<Gate>* nodePos );
    Gate( const Gate& aGate );
    ~Gate();

    //  get functions.. 
    BallGenerator*   getGateBallCCW1() const;
    BallGenerator*   getGateBallCCW2() const;
    BallGenerator*   getGateBallCCW3() const;

    BallGenerator*   getStartBall() const;
    VDVertex* getStartVertex() const;
    
    rg_dNode<Gate>* getNodePosInStack() const;

	rg_FLAG getEdgeType() const;
	BallGenerator* getEndBall() const;
	Sphere  getEmptyTangentSphereAtEndVertex() const;

    rg_FLAG isTheSameGate( const Gate& theOtherGate, VDVertex* endVertexOfTheOtherGate  ) const;
    rg_FLAG isTheSameGate( BallGenerator* gateBallCCW1, BallGenerator* gateBallCCW2, BallGenerator* gateBallCCW3, VDVertex* startVertex  ) const;

    
    //rg_dList<BallGenerator*>* getEndBallList();
    rg_dList< pair<BallGenerator*, Sphere> >* getEndBallList();
    BallGenerator*            getSingleEndBall() const;
    void addEndBall(BallGenerator* endBall);
    void emptyEndBallList();
    rg_INT          getNumBallsTangentToSphere() const;
    BallGenerator** getBallsTangentToSphere();
    rg_Point3D*     getNormalizedTangentPointsBetweenBallsAndSphere();
    void   getBallsTangentToSphere(rg_dList<BallGenerator*>& ballsTangentToSphere);

    rg_FLAG         getNumVerticesBy4Balls() const;

    rg_INT          getNumExcludedBalls() const;
    BallGenerator** getExcludedBalls() const;
    void            setExcludedBalls(const rg_INT& numExcludedBalls );
    void            setExcludedBalls( const rg_INT& numExcludedBalls, BallGenerator** excludedBalls);
    void            setOneExcludedBall(const rg_INT& i, BallGenerator* oneExcludedBall);
    void            setCheckOfAllExcludedBalls(const rg_FLAG& check);

    //  set functions..
    void setGateBallCCW1( BallGenerator* gateBallCCW1 );
    void setGateBallCCW2( BallGenerator* gateBallCCW2 );
    void setGateBallCCW3( BallGenerator* gateBallCCW3 );

    void setStartBall( BallGenerator* startBall );
    void setStartVertex( VDVertex* startVertex );

    void setNodePosInStack( rg_dNode<Gate>* nodePos );

    void setGate( BallGenerator* gateBallCCW1, BallGenerator* gateBallCCW2, BallGenerator* gateBallCCW3, BallGenerator* startBall, VDVertex* startVertex );
    void setGate( BallGenerator* gateBallCCW1, BallGenerator* gateBallCCW2, BallGenerator* gateBallCCW3, BallGenerator* startBall, VDVertex* startVertex, rg_dNode<Gate>* nodePos );


	void setEdgeType(const rg_FLAG& edgeType);
	void setEndBall(BallGenerator* endBall);
	void setEmptyTangentSphereAtEndVertex(const Sphere& emptySphere);

    //  operator overloading..
    Gate& operator =( const Gate& aGate );


	//  2006. 8. 17  By Youngsong Cho
    rg_INT     getNumOfEndBall() const;
    rg_Point3D getNormalOfCenterPlane() const;
    inline rg_Point3D getNormalOfEdgePlane() const { return m_normalOfEdgePlane; }
    rg_REAL    getCoeffOfCenterPlane() const;
    Sphere     getMinTangentSphere() const;
    Sphere     getMinTangentSphere(const rg_INT& i) const;

    rg_Point3D getAxisPoint() const;
    Plane      getCenterPlane() const;
    Plane      getEdgePlane() const;


	void    obtainGeometricProperties();
    void    determineEdgeType();

	rg_REAL angleCCWOfTwoUnitVectors(const rg_Point3D& normal, 
                                     const rg_Point3D& vector1, 
                                     const rg_Point3D& vector2);
	rg_REAL angleCCW(const rg_Point2D& vector1, const rg_Point2D& vector2);
    void    findEmptyTangentSphereWithMinAngleForParabolicEdge(rg_dList<BallGenerator*>& candidateBalls);
    void    findEmptyTangentSphereWithMinAngleForLinearEdge(rg_dList<BallGenerator*>& candidateBalls);
    void    findEmptyTangentSphereWithMinAngleForEllipticEdge(rg_dList<BallGenerator*>& candidateBalls);

    static  bool doesDefineEllipticCurve(const Sphere& ball1, const Sphere& ball2, const Sphere& ball3);
};

} // namespace GeometryTier

} // namespace V


#endif

