#ifndef _BETAVERTEX_H
#define _BETAVERTEX_H

#include "rg_Const.h"
#include "Sphere.h"
#include "Ball.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"
#include "rg_BetaSpan.h"



namespace V {

namespace GeometryTier {

class BetaEdge;
class BetaFace;
class BetaCell;
class VDCell;

class BetaVertex : public TopologicalEntity
{
private:
    BetaEdge* m_firstEdge;      // one of edges in the biggest world of vertex
    BetaCell* m_firstCell;      // one of cells in the biggest world of vertex
    
    rg_BetaSpan  m_betaSpan;
    rg_BOOL   m_visited;

    rg_INT    m_depth;

    Ball*     m_ball;

    VDCell*   m_vCell;

public:
    BetaVertex();
    BetaVertex(const rg_INT& ID);
    BetaVertex(const rg_INT& ID, Ball* ball);
    BetaVertex(const BetaVertex& betaVtx);
    ~BetaVertex();

    inline rg_BOOL isVisited() const { return m_visited; }
    inline void    isVisited(const rg_BOOL& visited) { m_visited = visited; }
    rg_BOOL        isVirtual() const;

    inline  rg_INT getDepth() const { return m_depth; }
    inline  void   setDepth(const rg_INT& depth) { m_depth = depth; }

    inline VDCell* getVCell() const { return m_vCell; }

    BetaEdge* getFirstEdge() const;
    BetaCell* getFirstCell() const;
    rg_BetaSpan  getBetaSpan() const;
    inline const rg_BetaSpan*  getpBetaSpan() const { return &m_betaSpan; }
    rg_REAL   getMinValueForValidSimplex() const;

    Ball*     getBallProperty() const;
    Sphere    getBall() const;
    void*     getProperty();

    void      setFirstEdge(BetaEdge* firstEdge);
    void      setFirstCell(BetaCell* firstCell);
    void      setBetaSpan(const rg_BetaSpan& betaSpan);
    void      setBall(const Sphere& ball);
    void      setProperty(void* property);
    void      setBallProperty(Ball* ball);

	// JHRYU
	void         shiftBetaSpan(const rg_REAL& delta);

    BetaVertex& operator =(const BetaVertex& betaVtx);

    inline rg_INT  getBoundingState(const rg_REAL& beta) const { return m_betaSpan.findBoundingState(beta); }

    rg_BOOL isOnConvexHull() const;

    void    computeBetaSpan();

    ///////////////////////////////////////////////////////////////////////////
    //  Quasi-operator
    //  
    //  1. Q(v, Y | W): Query primitives for intra-world.
    void searchVerticesInIntraWorld( rg_dList<BetaVertex*>& vertexList, BetaFace* world = rg_NULL) const;
    void searchEdgesInIntraWorld(    rg_dList<BetaEdge*>&   edgeList,   BetaFace* world = rg_NULL) const;
    void searchFacesInIntraWorld(    rg_dList<BetaFace*>&   faceList,   BetaFace* world = rg_NULL) const;
    void searchCellsInIntraWorld(    rg_dList<BetaCell*>&   cellList,   BetaFace* world = rg_NULL) const;
    void searchFiniteVerticesInIntraWorld( rg_dList<BetaVertex*>& finiteVertexList, BetaFace* world = rg_NULL) const;
    void searchFiniteEdgesInIntraWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world = rg_NULL) const;
    void searchFiniteFacesInIntraWorld(    rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world = rg_NULL) const;
    void searchFiniteCellsInIntraWorld(    rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world = rg_NULL) const;

    //  2. Query primitives for inter-world.
    void searchSmallWorlds(      rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL) const;
    void searchWholeSmallWorlds( rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL) const;
//    void searchBigWorlds(        rg_dList<BetaFace*>& bigWorldList, BetaFace* world = rg_NULL) const;
//    void inquireWholeBigWorldsAtVertex( rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL);


    void searchVerticesInSmallWorld( rg_dList<BetaVertex*>& vertexList, BetaFace* world = rg_NULL) const;
    void searchEdgesInSmallWorld(    rg_dList<BetaEdge*>&   edgeList,   BetaFace* world = rg_NULL) const;
    void searchFacesInSmallWorld(    rg_dList<BetaFace*>&   faceList,   BetaFace* world = rg_NULL) const;
    void searchCellsInSmallWorld(    rg_dList<BetaCell*>&   cellList,   BetaFace* world = rg_NULL) const;
    void searchFiniteVerticesInSmallWorld( rg_dList<BetaVertex*>& finiteVertexList, BetaFace* world = rg_NULL) const;
    void searchFiniteEdgesInSmallWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world = rg_NULL) const;
    void searchFiniteFacesInSmallWorld(    rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world = rg_NULL) const;
    void searchFiniteCellsInSmallWorld(    rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world = rg_NULL) const;


    //  3. Q(f, Y): Query primitives for whole-world.
    void searchVerticesInWholeWorld( rg_dList<BetaVertex*>& vertexList ) const;
    void searchEdgesInWholeWorld(    rg_dList<BetaEdge*>&   edgeList   ) const;
    void searchFacesInWholeWorld(    rg_dList<BetaFace*>&   faceList   ) const;
    void searchCellsInWholeWorld(    rg_dList<BetaCell*>&   cellList   ) const;
    void searchFiniteVerticesInWholeWorld( rg_dList<BetaVertex*>& finiteVertexList ) const;
    void searchFiniteEdgesInWholeWorld(    rg_dList<BetaEdge*>&   finiteEdgeList   ) const;
    void searchFiniteFacesInWholeWorld(    rg_dList<BetaFace*>&   finiteFaceList   ) const;
    void searchFiniteCellsInWholeWorld(    rg_dList<BetaCell*>&   finiteCellList   ) const;


    void searchFiniteVerticesInSmallerWorldViaGateEdge(BetaEdge* gateEdge, rg_dList<BetaVertex*>& finiteVertexList   ) const;
    void searchFiniteEdgesInSmallerWorldViaGateEdge(   BetaEdge* gateEdge, rg_dList<BetaEdge*>&   finiteEdgeList   ) const;
    void searchFiniteFacesInSmallerWorldViaGateEdge(   BetaEdge* gateEdge, rg_dList<BetaFace*>&   finiteFaceList   ) const;
    void searchFiniteCellsInSmallerWorldViaGateEdge(   BetaEdge* gateEdge, rg_dList<BetaCell*>&   finiteCellList   ) const;

    //
    //  END of Quasi-operator
    ///////////////////////////////////////////////////////////////////////////};



    ///////////////////////////////////////////////////////////////////////////
    //  Beta-operator
    //  
    void searchVerticesInBetaComplex(const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const;
    void searchEdgesInBetaComplex(   const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList   ) const;
    void searchFacesInBetaComplex(   const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList   ) const;
    void searchCellsInBetaComplex(   const rg_REAL& beta, rg_dList<BetaCell*>&   betaCellList   ) const;

    void searchExtraneousVerticesInBetaComplex(const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const;
    void searchExtraneousEdgesInBetaComplex(   const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList   ) const;
    void searchExtraneousFacesInBetaComplex(   const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList   ) const;
    void searchExtraneousCellsInBetaComplex(   const rg_REAL& beta, rg_dList<BetaCell*>&   betaCellList   ) const;

    void searchVerticesInBetaComplex(const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaVertex*>& betaVertexList ) const;
    void searchEdgesInBetaComplex(   const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaEdge*>& betaEdgeList ) const;
    void searchFacesInBetaComplex(   const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaFace*>& betaFaceList ) const;
    void searchCellsInBetaComplex(   const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaCell*>& betaCellList ) const;
    
    void searchVerticesInBetaShape( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const;
    void searchEdgesInBetaShape(    const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList   ) const;
    void searchFacesInBetaShape(    const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList   ) const;


    rg_INT countNumGroupsOfFaceConnectedExtraneousCellsInBetaComplex( const rg_REAL& beta ) const;
    rg_INT countNumGroupsOfFaceConnectedInteriorCellsInBetaComplex(   const rg_REAL& beta ) const;
    void   searchGroupsOfFaceConnectedExtraneousCellsInBetaComplex(   const rg_REAL& beta, rg_dList< rg_dList<BetaCell*> >& groupsOfExtraneousCell ) const;
    void   searchGroupsOfFaceConnectedInteriorCellsInBetaComplex(     const rg_REAL& beta, rg_dList< rg_dList<BetaCell*> >& groupsOfInteriorCell ) const;

    void   searchGroupsOfFaceConnectedExtraneousCellsInGivenWorld(   const rg_REAL& beta, rg_dList< rg_dList<BetaCell*> >& groupsOfExtraneousCell, BetaFace* world = rg_NULL ) const;

//     void searchFacesAtStarInBetaComplex( const rg_REAL& beta, rg_dList<BetaFace*>&   faceList   ) const; 
//     void searchCellsAtStarInBetaComplex( const rg_REAL& beta, rg_dList<BetaFace*>&   faceList   ) const; 
//     void searchFacesAtStarInBetaShape( const rg_REAL& beta, rg_dList<BetaFace*>&   faceList   ) const; 
    //
    //  END of Quasi-operator
    ///////////////////////////////////////////////////////////////////////////};
    

    void connectVCell(   VDCell* v_cell);
    void disconnectVCell(VDCell* v_cell);

};

} // namespace GeometryTier

} // namespace V


#endif


