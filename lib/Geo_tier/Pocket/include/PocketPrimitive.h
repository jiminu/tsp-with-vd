#ifndef _POCKETPRIMITIVE_H__
#define _POCKETPRIMITIVE_H__

#include "rg_Const.h"
#include "rg_dList.h"

#include "TopologicalEntity.h"
#include "BetaFace.h"
#include "BetaEdge.h"
#include "BetaVertex.h"

#include "ManifoldizedBetaShape.h"


#include "ChemicalPropertiesOfPocket.h"
#include "rg_Atom.h"


//for Beta-Mol
//#include "Particle.h"


class PocketPrimitive : public TopologicalEntity
{
private:
	//boundary of pocket primitive
	//BetaFace*			        m_faceOnOuterBetaShape;
	
    //rg_dList< pair<MBSEdge*, rg_FLAG> > m_boundaryOfPocketPrimitive;

	//components of pocket primitive
	rg_dList<MBSFace*>			m_MBSFaces;


	rg_dList<MBSVertex*>		m_MBSVertices;
    rg_dList< void* >           m_ballProperties;
    
	
    rg_REAL                     m_maxDepth;
    rg_REAL                     m_aveDepth;
	rg_REAL                     m_sumOfAreaOfMBSFaces; // CMKIM ADDED ON 090509

	ChemicalPropertiesOfPocket	m_chemicalPropertiesOfPocket;

public:
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Constructor
	PocketPrimitive();
	//PocketPrimitive(const rg_INT& ID, BetaFace* faceOnOuterBShape);
	//PocketPrimitive(BetaFace* faceOnOuterBShape, Ridge* ridge1, Ridge* ridge2, Ridge* ridge3 );
	PocketPrimitive(const PocketPrimitive& tempPocketPrimitive);
	~PocketPrimitive();

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Get Functions
	//BetaFace*         			getFacesOnOuterBetaShape();
	rg_dList<MBSFace*>*			getMBSFaces();
	rg_dList<MBSVertex*>*		getMBSVertices();
		void fillMBSVertexList();

    
    //rg_dList< pair<MBSEdge*, rg_FLAG> >* getBoundaryOfPocketPrimitive();
    
    rg_INT  getNumMBSFaces() const;

    rg_REAL getMaximumDepth() const;
    rg_REAL getAverageDepth() const;

    rg_dList< void* >*  getBallPropertyList();
		void fillBallPropertyList();

    ////////////////// CMKIM ADDED ON 090509 //////////////////
    rg_REAL            getAreaOfMBSFaces();
        rg_REAL computeAreaOfMBSFace( MBSFace* givenMBSFace );
        rg_FLAG isAreaOfMBSFaceZero( MBSFace* givenMBSFace );
    ///////////////////////////////////////////////////////////

	ChemicalPropertiesOfPocket*	getpChemicalPropertiesOfPocket();
		void computeChemicalProperties();

	

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Set Functions
	//void setFaceOnOuterBetaShape( BetaFace* tempFaceOfOuterBetaShape);
	
    void addMBSFaceToThisPrimitive(MBSFace* newFace);
    void addMBSVertexToThisPrimitive(MBSVertex* newVertex);
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Operator Overloadings
	PocketPrimitive& operator =(const PocketPrimitive& tempPocketPrimitive);

	rg_FLAG isItLengthZeroEdge( MBSEdge* edge );

//     void addAndRefineBoundary( Ridge* ridges, BetaEdge** boundingBCEdges, rg_BOOL* orientOfBoundingBCEdges );
//         void refineBoundary( rg_dList< pair< MBSEdge*, rg_FLAG > >* innerRidges  );
// 
// 
//     void addBoundaryByRidgeAndOrientation(Ridge* givenRidge, const rg_FLAG& orientOfGivenRidge);
//     void reverseBoundary();
//         void reverseOrientationOfBoundary();
//     void refineBoundaryByRemovingRedundentEdges();
//     void refineBoundaryByRemovingLoop();
//         void removeLoop( rg_dList< MBSVertex* >& redundantVertices );
//         rg_FLAG isItLengthZeroEdge( MBSEdge* edge );
//         void checkAndRemoveZeroLengthLoop( MBSVertex* nextSideVertex, rg_dNode< pair<MBSEdge*, rg_FLAG> >*& currNode );
//         rg_dNode< pair<MBSEdge*, rg_FLAG> >* getPrevRealEdgeNode( rg_dNode< pair<MBSEdge*, rg_FLAG> >* currNode );
//         
//     void addSimplicesOfBoundary();
// 
//     rg_FLAG isThisVertexInBoundary(MBSVertex* givenVertex);
/*
    void makeBoundaryByAligningEdgesInInnerRidge();
    
    void initializeBoundaryByValidFirstEdge(BetaComplex* betaComplex, const rg_REAL& innerBetaVal);


    void makePreliminaryPrimitiveByAddingSimplexesInBoundary();

    void becomeMatePrimitive(rg_dList<PocketPrimitive*>& listOfMatePrimitives);

    void rearrangeReversedRidge();

 */
    void removeAll();
	void sortPDBAtomsWRTAtomSerialNum();
		void swapNode( rg_dNode< void* >* node1, rg_dNode< void* >* node2 );
		
};

#endif