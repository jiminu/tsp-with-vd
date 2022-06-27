#ifndef _SIMPLEXQUERYINGMACHINEINQTBYKDTREE_H
#define _SIMPLEXQUERYINGMACHINEINQTBYKDTREE_H

#include "rg_dList.h"

#include "BetaVertex.h"
#include "BetaEdge.h"
#include "BetaFace.h"
#include "BetaCell.h"

#include "KDTree.h"


namespace V {

namespace GeometryTier {


class BetaUniverse;

class SimplexQueryingMachineInQTByKDTree
{
private:
	KDTree<BetaCell*>	m_cellsInKDTree;
	KDTree<BetaFace*>	m_facesInKDTree;
	KDTree<BetaEdge*>	m_edgesInKDTree;
	KDTree<BetaVertex*>	m_verticesInKDTree;

	rg_dList<BetaFace*>	m_ellipticFacesOfAnomalyCase;
    rg_dList<BetaEdge*>	m_gateEedges;

public:
	SimplexQueryingMachineInQTByKDTree();	
	SimplexQueryingMachineInQTByKDTree(BetaUniverse* betaUniverse);
//	SimplexQueryingMachineInQTByKDTree(rg_dList<BetaCell>&  cellsInBetaUniverse, rg_dList<BetaFace>& facesInBetaUniverse,
//									   rg_dList<BetaEdge>& edgesInBetaUniverse, rg_dList<BetaVertex>& verticesInBetaUniverse);	
	~SimplexQueryingMachineInQTByKDTree();


	void makeKDTreeForSearchingSimplex(BetaUniverse* betaUniverse);
//	void makeKDTreeForSearchingSimplex(rg_dList<BetaCell>&  cellsInBetaUniverse, rg_dList<BetaFace>& facesInBetaUniverse,
//									   rg_dList<BetaEdge>& edgesInBetaUniverse, rg_dList<BetaVertex>& verticesInBetaUniverse);

	
	// Sixteen primitive query operators //////////////////////////////////////////////		
	void huntCellsForGivenBoundingState(	const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaCell*>& betaCellList);
	void huntFacesForGivenBoundingState(	const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaFace*>& betaFaceList);
	void huntEdgesForGivenBoundingState(	const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaEdge*>& betaEdgeList);
	void huntVerticesForGivenBoundingState( const rg_REAL& beta, const rg_INT& boundingState, rg_dList<BetaVertex*>& betaVertexList);
	// argument ordering
	
	void huntCellsForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaCell*>& betaCellList);
	void huntFacesForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaFace*>& betaFaceList);
	void huntEdgesForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaEdge*>& betaEdgeList);
	void huntVerticesForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, rg_dList<BetaVertex*>& betaVertexList);


	// 10 high-level query(HQ)
	// HQ1
	void huntFacesInSqueezedHull( rg_dList<BetaFace*>& betaFaceList);
	void huntEdgesInSqueezedHull( rg_dList<BetaEdge*>& betaEdgeList);
	void huntVerticesInSqueezedHull(	rg_dList<BetaVertex*>& betaVertexList);

	// HQ2
	void huntCellsInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaCell*>& betaCellList );
	void huntFacesInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaFace*>& betaFaceList );
	void huntEdgesInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaEdge*>& betaEdgeList );
	void huntVerticesInBetaComplexInKDTree( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList );

	// HQ3
	void huntFacesInBetaShapeInKDTree( const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList );
	void huntEdgesInBetaShapeInKDTree( const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList );
	void huntVerticesInBetaShapeInKDTree( const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList );
	
	// HQ4
	void huntFacesInRegularShapeInKDTree( const rg_REAL& beta, rg_dList<BetaFace*>&   betaFaceList );
	void huntEdgesInRegularShapeInKDTree( const rg_REAL& beta, rg_dList<BetaEdge*>&   betaEdgeList );
	void huntVerticesInRegularShapeInKDTree(	const rg_REAL& beta, rg_dList<BetaVertex*>& betaVertexList );

	// HQ5
	void huntCellsInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaCell*>& betaCellListWithCellsForFromBeta );
	void huntFacesInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaFace*>& betaFaceListWithFacesForFromBeta );
	void huntEdgesInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaEdge*>& betaEdgeListWithEdgesForFromBeta );
	void huntVerticesInBetaComplexForIncresingBeta( const rg_REAL& fromBeta, const rg_REAL& toBeta, rg_dList<BetaVertex*>& betaVertexListWithVerticesForFromBeta );

	// HQ6 to HQ10 have to be implemented.

	/////////////////////////////////////////////////////////////////////////////////////////////////////

private:	
		void makeKDTreeOfCells(rg_dList<BetaCell>*  cellsInBetaUniverse);			
			void getBetaCellsForKDTree(rg_dList<BetaCell>*  cellsInBetaUniverse, rg_dList<BetaCell*>& cellList, rg_dList<double*>& valueList);
		void makeKDTreeOfFaces(rg_dList<BetaFace>* facesInBetaUniverse);		
			void getBetaFacesForKDTreeAndSetAnomalyCase(rg_dList<BetaFace>* facesInBetaUniverse, rg_dList<BetaFace*>& faceList, rg_dList<double*>& valueList);
		void makeKDTreeOfEdges(rg_dList<BetaEdge>* edgesInBetaUniverse);		
			void getBetaEdgesForKDTreeAndAddGateEdges(rg_dList<BetaEdge>* edgesInBetaUniverse, rg_dList<BetaEdge*>& edgeList, rg_dList<double*>& valueList);
		void makeKDTreeOfVertices(rg_dList<BetaVertex>* verticesInBetaUniverse);
			void getBetaVerticesForKDTree(rg_dList<BetaVertex>* verticesInBetaUniverse, rg_dList<BetaVertex*>& vertexList, rg_dList<double*>& valueList);

	void makeRangeForGivenBoundingState(const rg_REAL& beta, const rg_INT& boundingState, KDTreeInterval*& range);
	void makeRangeForGivenChangedState( const rg_REAL& fromBeta, const rg_REAL& toBeta, const rg_INT& fromBoundingState, const rg_INT& toBoundingState, KDTreeInterval*& range);
};

} // namespace GeometryTier

} // namespace V


#endif
