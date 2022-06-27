#ifndef _BETACELL_H
#define _BETACELL_H

#include "rg_Const.h"
#include "Sphere.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"
#include "rg_BetaSpan.h"

const rg_INT EIWDS_NUM_FACE_ON_CELL = 4;
const rg_INT EIWDS_NUM_VERTEX_ON_CELL = 4;


namespace V {

namespace GeometryTier {


class BetaVertex;
class BetaEdge;
class BetaFace;

class VDVertex;

class BetaCell : public TopologicalEntity
{
private:
    BetaFace*   m_face[EIWDS_NUM_FACE_ON_CELL];
    BetaVertex* m_vertex[EIWDS_NUM_VERTEX_ON_CELL];
    rg_BetaSpan    m_betaSpan;

    Sphere      m_minTangentSphere;
    rg_BOOL     m_visited;

    VDVertex*   m_vVertex;

public:
    BetaCell();
    BetaCell(const rg_INT& ID);
    BetaCell(const BetaCell& betaCell);
    ~BetaCell();
    

    inline rg_BOOL isVisited() const { return m_visited; }
    inline void    isVisited(const rg_BOOL& visited) { m_visited = visited; }
    rg_BOOL        isVirtual() const;

    inline VDVertex* getVVertex() const { return m_vVertex; }

	BetaFace**   getFaces();
    BetaFace*    getFace(const rg_INT& i) const;	
    rg_INT       getFacePos(const BetaFace* const face) const;
    BetaCell*    getNeighborCell(const rg_INT& i) const;	
	BetaVertex** getVertices();
    BetaVertex*  getVertex(const rg_INT& i) const;	
    rg_INT       getVertexPos(const BetaVertex* const vertex) const;
    inline Sphere       getMinTangentSphere() const { return m_minTangentSphere; }
    rg_BetaSpan     getBetaSpan() const;
    inline const rg_BetaSpan*  getpBetaSpan() const { return &m_betaSpan; }
    rg_REAL      getMinValueForValidSimplex() const;

    rg_INT       getDepth() const;

	void         setFaces(BetaFace** faces);
	void         setFace(const rg_INT& i, BetaFace* face);
	void         setVertices(BetaVertex** vertices);
	void         setVertex(const rg_INT& i, BetaVertex* vertex);
    void         setFaceAndItsMateVertex(const rg_INT& i, BetaFace* face, BetaVertex* vertex);
    void         addFaceAndItsMateVertex(BetaFace* face, BetaVertex* vertex);
    void         setBetaSpan(const rg_BetaSpan& betaSpan);

	// JHRYU
	void         shiftBetaSpan(const rg_REAL& delta);

	void setMinTangentSphere(const Sphere& sphere);

    BetaCell& operator=(const BetaCell& betaCell);

    inline rg_INT  getBoundingState(const rg_REAL& beta) const { return m_betaSpan.findBoundingState(beta); }

    rg_BOOL     isThere(BetaVertex* vertex) const;
    rg_BOOL     isThere(BetaEdge* edge) const;
    rg_BOOL     isThere(BetaFace* face) const;
    rg_BOOL     isAdjacent(const BetaCell* const  cell) const;
    rg_BOOL     isIntersectedAtEdgeWith(BetaFace* face) const;
    rg_BOOL     isIntersectedAtEdgeWith(BetaCell* cell) const;

    //rg_BOOL     isOnConvexHull() const;
    rg_INT      isMultiplicity() const;
    BetaCell*   getNeighborCellWithMultiplicity() const;

    BetaFace*   findSharingFaceWith(const BetaCell* const neighborCell) const;
    BetaEdge*   findEdge(const BetaVertex* const vertex1, const BetaVertex* const vertex2) const; 
    void        findFacesIncidentToVertex(const BetaVertex* const vertex, BetaFace** face) const; 
    void        findEdgesIncidentToVertex(const BetaVertex* const vertex, BetaEdge** edge) const; 
    void        findFacesIncidentToEdge(  const BetaEdge*   const edge,   BetaFace** face) const; 
    rg_REAL     computeAngleBtwFacesIncidnetToEdge( const BetaEdge*   const edge ) const;

	void        getNeighborCells(BetaCell** neighborCells) const;  
	BetaCell*   getNeighborCell(const BetaFace* const face) const;
    BetaVertex* getMateVertex(const BetaFace* const face) const;
    BetaFace*   getMateFace(const BetaVertex* const vertex) const;
    BetaCell*   getMateCell(const rg_INT& i_vtx) const;
    BetaCell*   getMateCell(const BetaVertex* const vertex) const;

	// JHRYU
    void        updateRadiusOfMinimumTangentSphere(const rg_REAL& deltaOfBallRadii);
    void        computeBetaSpan();
	rg_REAL     computeSignedVolume_old() const;
	// JHRYU
    rg_REAL     computeSignedVolume() const;

    ///////////////////////////////////////////////////////////////////////////
    //  Quasi-operator
    //  
    //  1. Q(c, Y | W): Query primitives for intra-world.
    void searchVerticesInIntraWorld(       rg_dList<BetaVertex*>& vertexList, BetaFace* world = rg_NULL) const;
    void searchEdgesInIntraWorld(          rg_dList<BetaEdge*>&   edgeList,   BetaFace* world = rg_NULL) const;
    void searchFacesInIntraWorld(          rg_dList<BetaFace*>&   faceList,   BetaFace* world = rg_NULL) const;
    void searchCellsInIntraWorld(          rg_dList<BetaCell*>&   cellList,   BetaFace* world = rg_NULL) const;

    void searchFiniteVerticesInIntraWorld( rg_dList<BetaVertex*>& finiteVertexList, BetaFace* world = rg_NULL) const;
    void searchFiniteEdgesInIntraWorld(    rg_dList<BetaEdge*>&   finiteEdgeList,   BetaFace* world = rg_NULL) const;
    void searchFiniteFacesInIntraWorld(    rg_dList<BetaFace*>&   finiteFaceList,   BetaFace* world = rg_NULL) const;
    void searchFiniteCellsInIntraWorld(    rg_dList<BetaCell*>&   finiteCellList,   BetaFace* world = rg_NULL) const;

//     //  2. Query primitives for inter-world.
//     void inquireSmallWorldsAtCell(  rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL);

    //  3. Q(c, Y): Query primitives for whole-world.
    void searchVerticesInWholeWorld(       rg_dList<BetaVertex*>& vertexList ) const;
    void searchEdgesInWholeWorld(          rg_dList<BetaEdge*>&   edgeList   ) const;
    void searchFacesInWholeWorld(          rg_dList<BetaFace*>&   faceList   ) const;
    void searchCellsInWholeWorld(          rg_dList<BetaCell*>&   cellList   ) const;

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
    void searchCellsInBetaComplex( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList ) const;
    
    void searchEdgesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaEdge*>& betaEdgeList ) const;
    void searchFacesInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaFace*>& betaFaceList ) const;
    void searchCellsInBetaComplex( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaCell*>& betaCellList ) const;

    //
    //  END of Quasi-operator
    ///////////////////////////////////////////////////////////////////////////};
    
    
    // function to connect with vertex in VD
    void connectVVertex(   VDVertex* v_vtx);
    void disconnectVVertex(VDVertex* v_vtx);

};

} // namespace GeometryTier

} // namespace V


#endif


