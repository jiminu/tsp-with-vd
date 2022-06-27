#ifndef _BETAFACE_H
#define _BETAFACE_H

#include "rg_Const.h"
#include "Sphere.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"
#include "rg_BetaSpan.h"


namespace V {

namespace GeometryTier {


class BetaVertex;
class BetaEdge;
class BetaCell;

class VDEdge;

const rg_INT EIWDS_NUM_EDGE_ON_FACE = 3;


class BetaFace : public TopologicalEntity
{
private:
    BetaCell* m_leftCell;
    BetaCell* m_rightCell;
    BetaEdge* m_edge[EIWDS_NUM_EDGE_ON_FACE];
    rg_BOOL   m_orientation[EIWDS_NUM_EDGE_ON_FACE];

    rg_BetaSpan  m_betaSpan;

    rg_BOOL   m_visited;

    VDEdge*   m_vEdge;

public:
    BetaFace();
    BetaFace(const rg_INT& ID);
	BetaFace(BetaEdge** edge);
    BetaFace(BetaCell* leftCell, BetaCell* rightCell);
    BetaFace(const rg_INT& ID, BetaCell* leftCell, BetaCell* rightCell);
    BetaFace(const BetaFace& betaFace);
    ~BetaFace();

    inline rg_BOOL isVisited() const { return m_visited; }
    inline void    isVisited(const rg_BOOL& visited) { m_visited = visited; }
    rg_BOOL        isVirtual() const;

    inline VDEdge* getVEdge() const { return m_vEdge; }

    BetaCell*  getLeftCell() const;
    BetaCell*  getRightCell() const;
    BetaCell*  getOppositeCell(BetaCell* incidentCell) const;  

	BetaEdge** getEdges();
    BetaEdge*  getEdge(const rg_INT& i) const;
    rg_INT     getPosOfEdge(BetaEdge* edge) const;
    rg_BOOL*   getEdgeOrientations();
    rg_BOOL    getEdgeOrientation(const rg_INT& i) const;
    rg_BOOL    getEdgeOrientation(const BetaEdge* const edge) const;
    rg_BetaSpan   getBetaSpan() const;
    inline const rg_BetaSpan*  getpBetaSpan() const { return &m_betaSpan; }
    rg_REAL    getMinValueForValidSimplex() const;

    rg_INT      getDepth() const;


	void setLeftCell(BetaCell*  leftCell);
    void setRightCell(BetaCell* rightCell);
	void setEdges(BetaEdge** edges);
	void setEdge(const rg_INT& i, BetaEdge* edge);
    void setEdge(const rg_INT& i, const rg_BOOL& edgeOrientation, BetaEdge* edge);
    void setEdgeOrientations(rg_BOOL* edgeOrientation);
    void setEdgeOrientation(const rg_INT& i, const rg_BOOL& edgeOrientation);
	void addEdge(BetaEdge* edge);
    void setBetaSpan(const rg_BetaSpan& betaSpan);
    
	// JHRYU
	void         shiftBetaSpan(const rg_REAL& delta);

    BetaFace& operator=(const BetaFace& betaFace);

    inline rg_INT  getBoundingState(const rg_REAL& beta) const { return m_betaSpan.findBoundingState(beta); }

    rg_BOOL   isIsolatedFace() const;
    rg_BOOL   isOnConvexHull() const;
    rg_BOOL   isSingularEllipticFace(const rg_REAL& beta) const;

    rg_BOOL   isIn2AdjacencyAnomaly() const;
    rg_BOOL   isIn3AdjacencyAnomaly() const;
    rg_BOOL   isIn4AdjacencyAnomaly() const;

    rg_BOOL   isThere(   const BetaVertex* const vertex) const;
    rg_BOOL   isThere(   const BetaEdge* const   edge)   const;
    rg_BOOL   isAdjacent(const BetaFace* const   face)   const;

    BetaEdge* findSharingEdgeWith(const BetaFace* const adjacentFace) const;
    BetaEdge* findEdge(const BetaVertex* const vertex1, const BetaVertex* const vertex2) const; 
    void      findEdgeIncidentToVertex(const BetaVertex* const vertex, BetaEdge** edge) const; 
    BetaEdge*   findMateEdgeOfVertex(const BetaVertex* const vertex) const;
    BetaVertex* findMateVertexOfEdge(const BetaEdge* const edge) const;
    BetaFace* findFaceIncidentToEdgeInCell(const BetaEdge* const edge, const BetaCell* const cell) const;

    void      computeBetaSpan();
private:
    void      computeBetaSpanOfNonEllipticFace(const Sphere& minTangentSphere);
    void      computeBetaSpanOfEllipticFace(const Sphere& minTangentSphere, const Sphere& maxTangentSphere);
        rg_BOOL isAttached(const rg_FLAG& isEllipticFaceAttached, const Sphere& minTangentSphere) const;
        rg_BOOL isNonEllipticFaceAttached() const;
public:
    ///////////////////////////////////////////////////////////////////////////
    //  Quasi-operator
    //  
    //  1. Q(f, Y | W): Query primitives for intra-world.
    void searchVerticesInIntraWorld( BetaVertex** vertexArray,          BetaFace* world = rg_NULL) const;
    void searchVerticesInIntraWorld( rg_dList<BetaVertex*>& vertexList, BetaFace* world = rg_NULL) const;
    void searchEdgesInIntraWorld(    rg_dList<BetaEdge*>&   edgeList,   BetaFace* world = rg_NULL) const;
    void searchFacesInIntraWorld(    rg_dList<BetaFace*>&   faceList,   BetaFace* world = rg_NULL) const;
    void searchCellsInIntraWorld(    rg_dList<BetaCell*>&   cellList,   BetaFace* world = rg_NULL) const;

    void searchFiniteVerticesInIntraWorld( rg_dList<BetaVertex*>& finiteVertexList, BetaFace* world = rg_NULL) const;
    void searchFiniteEdgesInIntraWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world = rg_NULL) const;
    void searchFiniteFacesInIntraWorld(    rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world = rg_NULL) const;
    void searchFiniteCellsInIntraWorld(    rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world = rg_NULL) const;

//     //  2. Query primitives for inter-world.
//     void inquireSmallWorldsAtFace( rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL);

    //  3. Q(f, Y): Query primitives for whole-world.
    void searchVerticesInWholeWorld( rg_dList<BetaVertex*>& vertexList ) const;
    void searchEdgesInWholeWorld(    rg_dList<BetaEdge*>&   edgeList   ) const;
    void searchFacesInWholeWorld(    rg_dList<BetaFace*>&   faceList   ) const;
    void searchCellsInWholeWorld(    rg_dList<BetaCell*>&   cellList   ) const;

    void searchFiniteVerticesInWholeWorld( rg_dList<BetaVertex*>& finiteVertexList ) const;
    void searchFiniteEdgesInWholeWorld(    rg_dList<BetaEdge*>&   finiteEdgeList   ) const;
    void searchFiniteFacesInWholeWorld(    rg_dList<BetaFace*>&   finiteFaceList   ) const;
    void searchFiniteCellsInWholeWorld(    rg_dList<BetaCell*>&   finiteCellList   ) const;
    //
    //  END of Quasi-operator
    ///////////////////////////////////////////////////////////////////////////};


    ///////////////////////////////////////////////////////////////////////////
    //  Beta-operator
    //  
    void searchFacesInBetaComplex( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const;
    void searchCellsInBetaComplex( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList ) const;
    
    void searchFacesInBetaShape(   const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const;
    //
    //  END of Quasi-operator
    ///////////////////////////////////////////////////////////////////////////};


    void connectVEdge(   VDEdge* v_edge);
    void disconnectVEdge(VDEdge* v_edge);
};

} // namespace GeometryTier

} // namespace V


#endif


