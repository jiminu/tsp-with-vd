#ifndef _BETAEDGE_H
#define _BETAEDGE_H

#include "rg_Const.h"
#include "Sphere.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"
#include "rg_BetaSpan.h"


namespace V {

namespace GeometryTier {


class BetaVertex;
class BetaFace;
class BetaCell;
class VDFace;

class BetaEdge : public TopologicalEntity
{
private:
    BetaVertex*         m_startVertex;
    BetaVertex*         m_endVertex;
    BetaFace*           m_firstFace;
    rg_dList<BetaFace*> m_smallWorlds;

    rg_BetaSpan            m_betaSpan;

    rg_BOOL             m_visited;

    VDFace*             m_vFace;

public:
    BetaEdge();
    BetaEdge(BetaVertex* startVertex, BetaVertex* endVertex);
    BetaEdge(const rg_INT& ID, BetaVertex* startVertex, BetaVertex* endVertex);
    BetaEdge(const rg_INT& ID, BetaVertex* startVertex, BetaVertex* endVertex, BetaFace* firstFace);
    BetaEdge(const BetaEdge& betaEdge);
    ~BetaEdge();

    inline rg_BOOL isVisited() const { return m_visited; }
    inline void    isVisited(const rg_BOOL& visited) { m_visited = visited; }
    rg_BOOL        isVirtual() const;

    inline VDFace* getVFace() const { return m_vFace; }

    BetaVertex* getStartVertex() const;
    BetaVertex* getEndVertex() const;
    BetaFace*   getFirstFace() const;

    rg_BOOL     isGateToSmallWorlds() const;
    rg_INT      getNumOfSmallWorlds() const;
    rg_dList<BetaFace*>* getSmallWorlds();

    rg_INT      getDepth() const;

    rg_BetaSpan    getBetaSpan() const;
    inline const rg_BetaSpan*  getpBetaSpan() const { return &m_betaSpan; }
    rg_REAL     getMinValueForValidSimplex() const;

    void        setStartVertex(BetaVertex* qtVertex);
    void        setEndVertex(BetaVertex* qtVertex);
    void        setFirstFace(BetaFace* firstFace);
    void        addSmallWorld(BetaFace* smallWorld);
    void        setBetaSpan(const rg_BetaSpan& betaSpan);

	// JHRYU
	void         shiftBetaSpan(const rg_REAL& delta);

    BetaEdge& operator=(const BetaEdge& betaEdge);
    
    inline rg_INT  getBoundingState(const rg_REAL& beta) const { return m_betaSpan.findBoundingState(beta); }

    rg_BOOL     isOnConvexHull() const;
    rg_BOOL     isThere(const BetaVertex* const vertex) const;

    BetaCell*   getCCWNextCell(const BetaFace* const currFace) const;
    BetaFace*   getCCWNextFace(const BetaFace* const currFace) const;
    BetaCell*   getCWNextCell(const BetaFace* const currFace) const;
    BetaFace*   getCWNextFace(const BetaFace* const currFace) const;

    void        computeBetaSpan();
private:
        rg_BetaSpan computeBetaSpan(const Sphere& minTangentSphere, BetaFace* world = rg_NULL) const;
        rg_BOOL  isAttached(const Sphere& minTangentSphere, BetaFace* world = rg_NULL) const;
public:

    ///////////////////////////////////////////////////////////////////////////
    //  Quasi-operator
    //  
    //  1. Q(e, Y | W): Query primitives for intra-world.
    void searchVerticesInIntraWorld( rg_dList<BetaVertex*>& vertexList, BetaFace* world = rg_NULL) const;
    void searchEdgesInIntraWorld(    rg_dList<BetaEdge*>&   edgeList,   BetaFace* world = rg_NULL) const;
    void searchFacesInIntraWorld(    rg_dList<BetaFace*>&   faceList,   BetaFace* world = rg_NULL) const;
    void searchCellsInIntraWorld(    rg_dList<BetaCell*>&   cellList,   BetaFace* world = rg_NULL) const;
    void searchFiniteEdgesInIntraWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world = rg_NULL) const;
    void searchFiniteFacesInIntraWorld(    rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world = rg_NULL) const;
    void searchFiniteCellsInIntraWorld(    rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world = rg_NULL) const;

    //  2. Query primitives for inter-world.
    void      searchSmallWorlds( rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL) const;
    BetaFace* searchBigWorld( BetaFace* world ) const;

    //  3. Q(e, Y): Query primitives for whole-world.
    void searchVerticesInWholeWorld( rg_dList<BetaVertex*>& vertexList ) const;
    void searchEdgesInWholeWorld(    rg_dList<BetaEdge*>&   edgeList   ) const;
    void searchFacesInWholeWorld(    rg_dList<BetaFace*>&   faceList   ) const;
    void searchCellsInWholeWorld(    rg_dList<BetaCell*>&   cellList   ) const;
    void searchFiniteEdgesInWholeWorld(    rg_dList<BetaEdge*>&   finiteEdgeList   ) const;
    void searchFiniteFacesInWholeWorld(    rg_dList<BetaFace*>&   finiteFaceList   ) const;
    void searchFiniteCellsInWholeWorld(    rg_dList<BetaCell*>&   finiteCellList   ) const;

    //
    //  END of Quasi-operator
    ///////////////////////////////////////////////////////////////////////////};



    ///////////////////////////////////////////////////////////////////////////
    //  Beta-operator
    //  
    void searchEdgesInBetaComplex( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList ) const;
    void searchFacesInBetaComplex( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const;
    void searchCellsInBetaComplex( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList ) const;
    
    void searchExtraneousEdgesInBetaComplex( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList ) const;
    void searchExtraneousFacesInBetaComplex( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const;
    void searchExtraneousCellsInBetaComplex( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList ) const;

    void searchEdgesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaEdge*>& betaEdgeList ) const;
    void searchFacesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaFace*>& betaFaceList ) const;
    void searchCellsInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaCell*>& betaCellList ) const;

    void searchEdgesInBetaShape(   const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList ) const;
    void searchFacesInBetaShape(   const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) const;

    //
    //  END of Quasi-operator
    ///////////////////////////////////////////////////////////////////////////};


    void connectVFace(   VDFace* v_face);
    void disconnectVFace(VDFace* v_face);
};

} // namespace GeometryTier

} // namespace V


#endif


