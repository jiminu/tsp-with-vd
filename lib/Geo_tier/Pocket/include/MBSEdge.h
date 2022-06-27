#ifndef _MBSEDGE_H
#define _MBSEDGE_H

#include "rg_Const.h"
#include "TopologicalEntity.h"
#include "Ball.h"

class MBSShell;
class MBSFace;
class MBSVertex;

namespace V{
    namespace GeometryTier{
class BetaEdge;
    }
}

using namespace V::GeometryTier;

class MBSEdge : public TopologicalEntity
{
private:
    MBSFace*     m_leftFace;
    MBSFace*     m_rightFace;

    MBSEdge*     m_leftHand;
    MBSEdge*	 m_rightHand;
    MBSEdge*     m_leftLeg;
    MBSEdge*	 m_rightLeg;

    MBSVertex*   m_startVertex;
    MBSVertex*   m_endVertex;

    BetaEdge*   m_originalBetaEdge;
    rg_BOOL     m_isArtificial;
    
    rg_BOOL     m_visited;

public:
    MBSEdge();
    MBSEdge(const rg_INT& ID);
    MBSEdge(const rg_INT& ID, MBSVertex* startVertex, MBSVertex* endVertex);
    MBSEdge(const MBSEdge& edge);
    ~MBSEdge();

    MBSFace*     getLeftFace() const;
    MBSFace*     getRightFace() const;
    MBSEdge*     getLeftHand() const;
    MBSEdge*     getRightHand() const;
    MBSEdge*     getLeftLeg() const;
    MBSEdge*     getRightLeg() const;
    MBSVertex*   getStartVertex() const;
    MBSVertex*   getEndVertex() const;

    BetaEdge*    getOriginalBetaEdge() const;
    rg_BOOL      isArtificial() const;
    
    MBSShell*    getShell() const;

    inline rg_BOOL  isVisited() const { return m_visited; }
    inline void     isVisited(const rg_BOOL& visited) { m_visited = visited; }

    void        setLeftFace(MBSFace* leftFace);
    void        setRightFace(MBSFace* rightFace);
    void        setLeftHand(MBSEdge* leftHand);
    void        setRightHand(MBSEdge* rightHand);
    void        setLeftLeg(MBSEdge* leftLeg);
    void        setRightLeg(MBSEdge* rightLeg);
    void        setStartVertex(MBSVertex* startVertex);
    void        setEndVertex(MBSVertex* endVertex);

    void        setOriginalBetaEdge( BetaEdge* originalBetaEdge );
    void        isArtificial( const rg_BOOL& isArtificial );
    
    MBSEdge&     operator =(const MBSEdge& edge);

    
    void        searchBoundingVertices( rg_dList<MBSVertex*>& vertexList );
    rg_FLAG     searchAdjacentEdges( rg_dList<MBSEdge*>& edgeList );
    void        searchAdjacentFaces( rg_dList<MBSFace*>& faceList );
};

#endif
