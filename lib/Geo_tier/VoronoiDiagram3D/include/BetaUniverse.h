#ifndef _BETAUNIVERSE_H
#define _BETAUNIVERSE_H

#include "rg_dList.h"

#include "BetaVertex.h"
#include "BetaEdge.h"
#include "BetaFace.h"
#include "BetaCell.h"

#include "Ball.h"

#include "rg_QuasiTriangulation.h"

#include <list>
#include <set>
using namespace std;

// JKKIM ///////////////////////////////////////////
//#include "SimplexQueryingMachineInQTByKDTree.h"
////////////////////////////////////////////////////

// #include <fstream>
// using namespace std;


namespace V {

namespace GeometryTier {


class BetaUniverse
{
private:
    rg_dList<Ball>       m_balls;

    rg_dList<BetaCell>   m_cells;
    rg_dList<BetaFace>   m_faces;
    rg_dList<BetaEdge>   m_edges;
    rg_dList<BetaVertex> m_vertices;

    rg_INT               m_numFiniteCells;
    rg_INT               m_numFiniteFaces;
    rg_INT               m_numFiniteEdges;
    rg_INT               m_numFiniteVertices;

    BetaCell**           m_sortedFiniteCells;
    BetaFace**           m_sortedFiniteFaces;
    BetaEdge**           m_sortedFiniteEdges;
    BetaVertex**         m_sortedFiniteVertices;

	// JKKIM ///////////////////////////////////////////
	//SimplexQueryingMachineInQTByKDTree* m_simplexQueryingMachine;
	////////////////////////////////////////////////////

public:
    ////  for time analysis of BC
    rg_REAL m_time4BC[4];
    rg_REAL m_time4BetaSpan[4];


    BetaUniverse();
    ~BetaUniverse();

	// JHRYU
    inline rg_dList<Ball>*         getBallList()   {return &m_balls;}
    inline rg_dList<BetaCell>*     getCellList()   {return &m_cells;}
    inline rg_dList<BetaFace>*     getFaceList()   {return &m_faces;}
    inline rg_dList<BetaEdge>*     getEdgeList()   {return &m_edges;}
    inline rg_dList<BetaVertex>*   getVertexList() {return &m_vertices;}

    inline const rg_dList<Ball>&         getBallList() const   {return m_balls;}
    inline const rg_dList<BetaCell>&     getCellList() const   {return m_cells;}
    inline const rg_dList<BetaFace>&     getFaceList() const   {return m_faces;}
    inline const rg_dList<BetaEdge>&     getEdgeList() const   {return m_edges;}
    inline const rg_dList<BetaVertex>&   getVertexList() const {return m_vertices;}


    inline rg_INT getNumBalls() const              {return m_balls.getSize();}
    inline rg_INT getNumCells() const              {return m_numFiniteCells;}
    inline rg_INT getNumFaces() const              {return m_numFiniteFaces;}
    inline rg_INT getNumEdges() const              {return m_numFiniteEdges;}
    inline rg_INT getNumVertices() const           {return m_numFiniteVertices;}

    inline BetaCell**   getSortedCells()           {return m_sortedFiniteCells;}
    inline BetaFace**   getSortedFaces()           {return m_sortedFiniteFaces;}
    inline BetaEdge**   getSortedEdges()           {return m_sortedFiniteEdges;}
    inline BetaVertex** getSortedVertices()        {return m_sortedFiniteVertices;}

	
	// JHRYU NOTE: VDW Volume 
    void rearrangeSimplexIDs();
        void rearrangeCellIDs();
        void rearrangeFaceIDs();
        void rearrangeEdgeIDs();
        void rearrangeVertexIDs();

	void enlargeRadiiOfBallsNUpdateBetaSpan(const rg_REAL& deltaOfBallRadii);
	void shrinkRadiiOfBallsNUpdateBetaSpan(const rg_REAL& deltaOfBallRadii);
	void resizeRadiiOfBallsNUpdateBetaSpan(const rg_REAL& deltaOfBallRadii);
		void updateBetaSpan(const rg_REAL& deltaOfBallRadii);
			void updateRadiiOfMinimumTangentSpheresForBetaCells(const rg_REAL& deltaOfBallRadii);
		
		void updateBetaSpan2(const rg_REAL& deltaOfBallRadii);
			void shiftBetaSpan(const rg_REAL& deltaOfBallRadii);

	void enlargeRadiiOfBalls(const rg_REAL& deltaOfBallRadii);
	void shrinkRadiiOfBalls(const rg_REAL& deltaOfBallRadii);
	


    void construct(QuasiTriangulation& qtInIWDS);
    void clean();

    BetaVertex* findVirtualVertex() const;


	void  huntCellsInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const;
	void  huntFacesInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const;
	void  huntEdgesInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const;
	void  huntVerticesInBetaComplex( const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const;

	// JHRYU
	void  huntSignificantCellsInBetaComplex(    const rg_REAL& beta, rg_dList<BetaCell*>&   betaCellList, const rg_REAL& tolerance = rg_MATH_RES ) const;

	void  huntCellsInBetaComplex(    const rg_REAL& beta, rg_dList<BetaCell*>&   betaCellList ) const;
	void  huntFacesInBetaComplex(    const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList ) const;
	void  huntEdgesInBetaComplex(    const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList ) const;
	void  huntVerticesInBetaComplex( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const;

	void  huntFacesInBetaShape(    const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList ) const;
	void  huntEdgesInBetaShape(    const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList ) const;
	void  huntVerticesInBetaShape( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) const;


    rg_INT countNumGroupsOfFaceConnectedExtraneousCellsInBetaComplex( const rg_REAL& beta );
    rg_INT countNumGroupsOfFaceConnectedInteriorCellsInBetaComplex(   const rg_REAL& beta ) const;
    void   searchGroupOfOutmostExtraneousCellsInBetaComplex(   const rg_REAL& beta, rg_dList<BetaCell*>& outmostExtraneousCells );
    void   searchGroupsOfFaceConnectedExtraneousCellsInBetaComplex(   const rg_REAL& beta, rg_dList< rg_dList<BetaCell*> >& groupsOfExtraneousCell );
    void   searchGroupsOfFaceConnectedInteriorCellsInBetaComplex(     const rg_REAL& beta, rg_dList< rg_dList<BetaCell*> >& groupsOfInteriorCell ) const;


	//////////////////////////////////////////////////////////////////////////////
    //  for Quasi-triangulation by Y. Cho 20130319
	rg_INT  findAll2AdjacencyAnomaly( list< pair<BetaCell*, BetaCell*> >& anomaliesWith2Adjacency ) const; 		
	rg_INT  findAll3AdjacencyAnomaly( list< pair<BetaCell*, BetaCell*> >& anomaliesWith3Adjacency ) const; 
	rg_INT  findAll4AdjacencyAnomaly( list< pair<BetaCell*, BetaCell*> >& anomaliesWith4Adjacency ) const; 
	rg_INT  findAllSmallWorlds(       list< BetaFace* >& smallWorlds ) const; 
	rg_INT  findAllIsolatedTriangles( list< BetaFace* >& isolatedTriangles ) const; 

    rg_INT  findFaceConnectedCellsInWorld( const BetaFace* const world, list< BetaCell* >& cellsInWorld ) const;

	//// JKKIM //////////////////////////////////////////////////////////////////////////
	//// Sixteen primitive query operators //////////////////////////////////////////////		
//	inline void huntCellsForGivenBoundingState(	const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaCell*>& betaCellList) {m_simplexQueryingMachine->huntCellsForGivenBoundingState( beta, boundingState, betaCellList);}
//	inline void huntFacesForGivenBoundingState(	const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaFace*>& betaFaceList) {m_simplexQueryingMachine->huntFacesForGivenBoundingState( beta, boundingState, betaFaceList);}
//	inline void huntEdgesForGivenBoundingState(	const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaEdge*>& betaEdgeList) {m_simplexQueryingMachine->huntEdgesForGivenBoundingState( beta, boundingState, betaEdgeList);}
//	inline void huntVerticesForGivenBoundingState( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaVertex*>& betaVertexList) {m_simplexQueryingMachine->huntVerticesForGivenBoundingState( beta, boundingState, betaVertexList);}
//	
//	//
//	inline void huntCellsForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaCell*>& betaCellList) {m_simplexQueryingMachine->huntCellsForGivenChangedState(fromBeta, toBeta, fromBoundingState, toBoundingState, betaCellList);}
//	inline void huntFacesForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaFace*>& betaFaceList) {m_simplexQueryingMachine->huntFacesForGivenChangedState(fromBeta, toBeta, fromBoundingState, toBoundingState, betaFaceList);}
//	inline void huntEdgesForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaEdge*>& betaEdgeList) {m_simplexQueryingMachine->huntEdgesForGivenChangedState(fromBeta, toBeta, fromBoundingState, toBoundingState, betaEdgeList);}
//	inline void huntVerticesForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaVertex*>& betaVertexList) {m_simplexQueryingMachine->huntVerticesForGivenChangedState(fromBeta, toBeta, fromBoundingState, toBoundingState, betaVertexList);}
//
//
//	// 10 high-level query(HQ)
//	// HQ1
//	inline void huntFacesInSqueezedHull( rg_dList<BetaFace*>& betaFaceList) {m_simplexQueryingMachine->huntFacesInSqueezedHull(betaFaceList);}
//	inline void huntEdgesInSqueezedHull( rg_dList<BetaEdge*>& betaEdgeList) {m_simplexQueryingMachine->huntEdgesInSqueezedHull(betaEdgeList);}
//	inline void huntVerticesInSqueezedHull(	rg_dList<BetaVertex*>& betaVertexList) {m_simplexQueryingMachine->huntVerticesInSqueezedHull(betaVertexList);}
//
//	//// HQ2
//	inline void huntCellsInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList ) {m_simplexQueryingMachine->huntCellsInBetaComplexInKDTree(beta, betaCellList);}
//	inline void huntFacesInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList ) {m_simplexQueryingMachine->huntFacesInBetaComplexInKDTree(beta, betaFaceList);}
//	inline void huntEdgesInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList ) {m_simplexQueryingMachine->huntEdgesInBetaComplexInKDTree(beta, betaEdgeList);}
//	inline void huntVerticesInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) {m_simplexQueryingMachine->huntVerticesInBetaComplexInKDTree(beta, betaVertexList);}
//
//	//// HQ3
//	inline void huntFacesInBetaShapeInKDTree( const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList ) {m_simplexQueryingMachine->huntFacesInBetaShapeInKDTree(beta, betaFaceList);}
//	inline void huntEdgesInBetaShapeInKDTree( const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList ) {m_simplexQueryingMachine->huntEdgesInBetaShapeInKDTree(beta, betaEdgeList);}
//	inline void huntVerticesInBetaShapeInKDTree( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) {m_simplexQueryingMachine->huntVerticesInBetaShapeInKDTree(beta, betaVertexList);}
//	//
//	// HQ4
//	inline void huntFacesInRegularShapeInKDTree( const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList ) {m_simplexQueryingMachine->huntFacesInRegularShapeInKDTree(beta, betaFaceList);}
//	inline void huntEdgesInRegularShapeInKDTree( const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList ) {m_simplexQueryingMachine->huntEdgesInRegularShapeInKDTree(beta, betaEdgeList);}
//	inline void huntVerticesInRegularShapeInKDTree(	const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList ) {m_simplexQueryingMachine->huntVerticesInRegularShapeInKDTree(beta, betaVertexList);}
//
//	// HQ5
//	inline void huntCellsInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaCell*>& betaCellListWithCellsForFromBeta ) {m_simplexQueryingMachine->huntCellsInBetaComplexForIncresingBeta(fromBeta, toBeta, betaCellListWithCellsForFromBeta);}
//	inline void huntFacesInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaFace*>& betaFaceListWithFacesForFromBeta ) {m_simplexQueryingMachine->huntFacesInBetaComplexForIncresingBeta(fromBeta, toBeta, betaFaceListWithFacesForFromBeta);}
//	inline void huntEdgesInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaEdge*>& betaEdgeListWithEdgesForFromBeta ) {m_simplexQueryingMachine->huntEdgesInBetaComplexForIncresingBeta(fromBeta, toBeta, betaEdgeListWithEdgesForFromBeta);}
//	inline void huntVerticesInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaVertex*>& betaVertexListWithVerticesForFromBeta ) {m_simplexQueryingMachine->huntVerticesInBetaComplexForIncresingBeta(fromBeta, toBeta, betaVertexListWithVerticesForFromBeta);}

	// HQ6 to HQ10 have to be implemented.

	///////////////////////////////////////////////////////////////////////////////////////////////////////


private:
	void computeBetaSpan();
        void computerBetaSpanOfCells();
        void computerBetaSpanOfFaces();
        void computerBetaSpanOfEdges();
        void computerBetaSpanOfVertices();

    void setDepthOfWorlds();
        void setDepthOfVerticesInWorlds();
            void setDepthOfVerticesInRootWorld();
            void setDepthOfVerticesInSmallWorlds();

    void generateSortedFiniteSimplexes();
        void generateSortedFiniteCells();
        void generateSortedFiniteFaces();
        void generateSortedFiniteEdges();
        void generateSortedFiniteVertices();

        
        friend int compareBetaSpanOfCells(const void* ts1, const void* ts2);
        friend int compareBetaSpanOfFaces(const void* ts1, const void* ts2);
        friend int compareBetaSpanOfEdges(const void* ts1, const void* ts2);
        friend int compareBetaSpanOfVertices(const void* ts1, const void* ts2);

//     void rearrangeSimplexIDs();
//         void rearrangeCellIDs();
//         void rearrangeFaceIDs();
//         void rearrangeEdgeIDs();
//         void rearrangeVertexIDs();
    ///////////////////////////////////////////////////////////////////////////
    //
    //  Converter
    //
    //  1. IWDS -> eIWDS
    void convertFromIWDSToEIWDS(QuasiTriangulation& qtInIWDS);
        BetaCell**   createCellsAndMakeRefTableForCells(QuasiTriangulation& qtInIWDS);
        BetaVertex** createVerticesAndMakeRefTableForVertices(QuasiTriangulation& qtInIWDS);
        void              setTopologyBetweenCellsAndVertices( QuasiTriangulation& qtInIWDS, 
                                                              BetaCell**     refTableForCells, 
                                                              BetaVertex**   refTableForVertices );
        void createFacesAndSetTopologyBetweenCellsAndFaces( QuasiTriangulation& qtInIWDS, 
                                                            BetaCell**     refTableForCells );
        void createEdgesAndSetTopologyAmongFacesEdgesAndVertices( QuasiTriangulation& qtInIWDS, 
                                                                  BetaCell**     refTableForCells,
                                                                  BetaVertex**   refTableForVertices );
            void setTopologyFromFaceToEdge( QuasiTriangulation& qtInIWDS, 
                                            BetaCell**          refTableForCells,
                                            BetaVertex**        refTableForVertices,
                                            const QTEdge&       edgeInIWDS,
                                            BetaEdge*           edgeInEIWDS );
            void setTopologyFromVerticesToEdges();

        BetaFace** createIsolatedFacesAndMakeRefTableForFaces( QuasiTriangulation& qtInIWDS );

        void setTopologyBetweenGateEdgesAndSmallWorlds( QuasiTriangulation& qtInIWDS, 
                                                        BetaCell**     refTableForCells,
                                                        BetaFace**     refTableForIsolatedFaces,
                                                        BetaVertex**   refTableForVertices );

            BetaEdge* setTopologyBetweenGateEdgeAndSmallWorldAsSetOfCells(
                                                                      BetaVertex**   vertexOfGate,
                                                                      BetaEdge*      gateEdge,
                                                                      QTTetrahedron*      cellInSmallWorldOfIWDS,
                                                                      QuasiTriangulation& qtInIWDS, 
                                                                      BetaCell**     refTableForCells );
            void setTopologyBetweenGateEdgeAndSmallWorldAsIsolatedFace( BetaEdge*      gateEdge,
                                                                        QTTetrahedron*      cellInSmallWorldOfIWDS,
                                                                        QuasiTriangulation& qtInIWDS, 
                                                                        BetaFace**     refTableForIsolatedFaces,
                                                                        BetaVertex**   refTableForVertices );

        void linkVD( QuasiTriangulation& qtInIWDS, BetaCell** refTableForCells, BetaVertex** refTableForVertices);
            void linkVVertexAndBetaCell( QuasiTriangulation& qtInIWDS, BetaCell** refTableForCells );
            void linkVCellAndBetaVertex( QuasiTriangulation& qtInIWDS, BetaVertex** refTableForVertices );
            void linkVEdgeAndBetaFace();
            void linkVFaceAndBetaEdge();

        void checkLinkVD2BC();

public:

    /*
    ///////////////////////////////////////////////////////////////////////////
    //
    //  Quasi-operator
    //
    //  1. Query primitives for intra-world.
    //
    //  1.1. Q(v, Y | W): Query primitives for intra-world.
    void inquireVerticesOfVertexInIntraWorld(BetaVertex* givenVertex, rg_dList<BetaVertex*>& vertexList, BetaFace* world = rg_NULL);
    void inquireEdgesOfVertexInIntraWorld(   BetaVertex* givenVertex, rg_dList<BetaEdge*>&   edgeList,   BetaFace* world = rg_NULL);
    void inquireFacesOfVertexInIntraWorld(   BetaVertex* givenVertex, rg_dList<BetaFace*>&   faceList,   BetaFace* world = rg_NULL);
    void inquireCellsOfVertexInIntraWorld(   BetaVertex* givenVertex, rg_dList<BetaCell*>&   cellList,   BetaFace* world = rg_NULL);
    
    //  1.2. Q(e, Y | W): Query primitives for intra-world.
    void inquireVerticesOfEdgeInIntraWorld( BetaEdge* givenEdge, rg_dList<BetaVertex*>& vertexList, BetaFace* world = rg_NULL);
    void inquireEdgesOfEdgeInIntraWorld(    BetaEdge* givenEdge, rg_dList<BetaEdge*>&   edgeList,   BetaFace* world = rg_NULL);
    void inquireFacesOfEdgeInIntraWorld(    BetaEdge* givenEdge, rg_dList<BetaFace*>&   faceList,   BetaFace* world = rg_NULL);
    void inquireCellsOfEdgeInIntraWorld(    BetaEdge* givenEdge, rg_dList<BetaCell*>&   cellList,   BetaFace* world = rg_NULL);
    
    //  1.3. Q(f, Y | W): Query primitives for intra-world.
    void inquireVerticesOfFaceInIntraWorld( BetaFace* givenFace, rg_dList<BetaVertex*>& vertexList, BetaFace* world = rg_NULL);
    void inquireEdgesOfFaceInIntraWorld(    BetaFace* givenFace, rg_dList<BetaEdge*>&   edgeList,   BetaFace* world = rg_NULL);
    void inquireFacesOfFaceInIntraWorld(    BetaFace* givenFace, rg_dList<BetaFace*>&   faceList,   BetaFace* world = rg_NULL);
    void inquireCellsOfFaceInIntraWorld(    BetaFace* givenFace, rg_dList<BetaCell*>&   cellList,   BetaFace* world = rg_NULL);
    
    //  1.4. Q(c, Y | W): Query primitives for intra-world.
    void inquireVerticesOfCellInIntraWorld(BetaCell* givenCell, rg_dList<BetaVertex*>& vertexList, BetaFace* world = rg_NULL);
    void inquireEdgesOfCellInIntraWorld(   BetaCell* givenCell, rg_dList<BetaEdge*>&   edgeList,   BetaFace* world = rg_NULL);
    void inquireFacesOfCellInIntraWorld(   BetaCell* givenCell, rg_dList<BetaFace*>&   faceList,   BetaFace* world = rg_NULL);
    void inquireCellsOfCellInIntraWorld(   BetaCell* givenCell, rg_dList<BetaCell*>&   cellList,   BetaFace* world = rg_NULL);
    
    //  2. Query primitives for inter-world.
    //
    void inquireSmallWorldsAtVertex(BetaVertex* givenVertex, rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL);
    void inquireSmallWorldsAtEdge(  BetaEdge*   givenEdge,   rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL);
    void inquireSmallWorldsAtFace(  BetaFace*   givenFace,   rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL);
    void inquireSmallWorldsAtCell(  BetaCell*   givenCell,   rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL);

    void inquireWholeSmallWorldsAtVertex(BetaVertex* givenVertex, rg_dList<BetaFace*>& smallWorldList, BetaFace* world = rg_NULL);

    //  3. Query primitives for whole-world.
    //
    //  3.1. Q(v, Y): Query primitives for whole-world.
    void inquireVerticesOfVertex(BetaVertex* givenVertex, rg_dList<BetaVertex*>& vertexList );
    void inquireEdgesOfVertex(   BetaVertex* givenVertex, rg_dList<BetaEdge*>&   edgeList   );
    void inquireFacesOfVertex(   BetaVertex* givenVertex, rg_dList<BetaFace*>&   faceList   );
    void inquireCellsOfVertex(   BetaVertex* givenVertex, rg_dList<BetaCell*>&   cellList   );
    
    //  3.2. Q(e, Y): Query primitives for whole-world.
    void inquireVerticesOfEdge( BetaEdge* givenEdge, rg_dList<BetaVertex*>& vertexList );
    void inquireEdgesOfEdge(    BetaEdge* givenEdge, rg_dList<BetaEdge*>&   edgeList   );
    void inquireFacesOfEdge(    BetaEdge* givenEdge, rg_dList<BetaFace*>&   faceList   );
    void inquireCellsOfEdge(    BetaEdge* givenEdge, rg_dList<BetaCell*>&   cellList   );
    
    //  3.3. Q(f, Y): Query primitives for whole-world.
    void inquireVerticesOfFace( BetaFace* givenFace, rg_dList<BetaVertex*>& vertexList );
    void inquireEdgesOfFace(    BetaFace* givenFace, rg_dList<BetaEdge*>&   edgeList   );
    void inquireFacesOfFace(    BetaFace* givenFace, rg_dList<BetaFace*>&   faceList   );
    void inquireCellsOfFace(    BetaFace* givenFace, rg_dList<BetaCell*>&   cellList   );
    
    //  3.4. Q(c, Y): Query primitives for whole-world.
    void inquireVerticesOfCell(BetaCell* givenCell, rg_dList<BetaVertex*>& vertexList );
    void inquireEdgesOfCell(   BetaCell* givenCell, rg_dList<BetaEdge*>&   edgeList   );
    void inquireFacesOfCell(   BetaCell* givenCell, rg_dList<BetaFace*>&   faceList   );
    void inquireCellsOfCell(   BetaCell* givenCell, rg_dList<BetaCell*>&   cellList   );

    //
    //  END of Quasi-operator
    //  
    ///////////////////////////////////////////////////////////////////////////
    */
    

//     void reportQT(ofstream& fout);
//     void reportBetaSpan();
// 
//     ///////////////////////////////////////////////////////////////////////////
//     //
//     //  Modules for testing quasi-operators.
//     //
//     void testQuery_v_C(ofstream& fout);
//     void testQuery_v_F(ofstream& fout);
//     void testQuery_v_E(ofstream& fout);
//     void testQuery_v_V(ofstream& fout);
// 
//     void testQuery_e_C(ofstream& fout);
//     void testQuery_e_F(ofstream& fout);
//     void testQuery_e_E(ofstream& fout);
// 
//     //
//     ///////////////////////////////////////////////////////////////////////////
};


int compareBetaSpanOfCells(const void* ts1, const void* ts2);
int compareBetaSpanOfFaces(const void* ts1, const void* ts2);
int compareBetaSpanOfEdges(const void* ts1, const void* ts2);
int compareBetaSpanOfVertices(const void* ts1, const void* ts2);


} // namespace GeometryTier

} // namespace V


#endif

