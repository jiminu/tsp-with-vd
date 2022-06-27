#ifndef _POCKET_H__
#define _POCKET_H__

#include "PocketPrimitive.h"
#include "ConstForPocket.h"


//namespace V{
//    namespace GeometryTier{
//        class BetaUniverse;
//    }
//}

class V::GeometryTier::BetaUniverse;

class ManifoldizedBetaShape;

const rg_REAL DEFAULT_INNER_BETA_VALUE = 1.4;
const rg_REAL DEFAULT_OUTER_BETA_VALUE = 3.0;

class Pocket
{
private:
	BetaUniverse*	m_qusiTriangulationForBeta;
    rg_REAL 					m_outerBetaVal;
    
    rg_REAL                     m_innerBetaVal;
	ManifoldizedBetaShape		m_innerBetaShape;

    
	PocketPrimitive*        	m_pocketPrimitives;     //ordered by ID of outerBetaShape Faces
    PocketPrimitive**			m_sortedPocketPrimitives;

	rg_INT						m_numOfPocketPrimitives;
    
public:
    
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Constructor
	Pocket();
	Pocket( const rg_REAL& outerBetaVal, ManifoldizedBetaShape* innerBetaShape );
	
	~Pocket();
	
    void clear();
	/////////////////////////////////////////////////////////////////////////////////
	//
	// Get Functions
	BetaUniverse*	getQusiTriangulationForBeta();
	rg_REAL						getOuterBetaVal() const;
    rg_REAL                     getInnerBetaVal() const;
    ManifoldizedBetaShape*      getpInnerBetaShape();
	
	PocketPrimitive*	        getPocketPrimitives();
    PocketPrimitive**			getSortedPocketPrimitives();

	rg_INT						getNumOfPocketPrimitives() const;
    

	/////////////////////////////////////////////////////////////////////////////////
	//
	// Set Functions
	void setQuasiTriangulationForBeta(BetaUniverse* tempQusiTriangulationForBeta);
	void setOuterBetaVal(const rg_REAL& outerBetaVal);
    void setInnerBetaVal(const rg_REAL& innerBetaVal);
    void setInnerBetaShape(ManifoldizedBetaShape* tempInnerBetaShape);


    
    
    void constructPocket();
	void constructPocket( BetaUniverse* innerBetaShape, const rg_REAL& outerBetaVal, const rg_REAL& innerBetaVal );
		void markVerticesOnPocket();
            void markVerticesOfVertexList( rg_dList< MBSVertex* >* vertices );
                rg_REAL findUpperBoundOfRegularInterval( BetaVertex* givenVertex );
		void markFacesOnPocket();
            void markFacesOfFaceList( rg_dList< MBSFace* >* faces );
        void generatePocketPrimitiveByClusteringFaces();
            void clusterFaces( rg_dList< MBSFace* >* faces, rg_dList< PocketPrimitive >& pocketPrimitiveList );
        void markVerticesToUnvisited();
        void markFacesOnCurrentPockets();

    
    void sortPocketPrimitives( const rg_INT& compareWRT = COMPARE_WRT_AREA_OF_FACES );

    void removeAll();

	
	void shavePocketPrimitives();
        void shrinkPocketPrimitives();
            void shrinkPocketPrimitive( PocketPrimitive* currPocketPrimitive );
                rg_BOOL isThisFaceOnBoundaryOfPocket( MBSFace* currFace );
                    void findRegularFacesBoundedByVertexOfGivenFace( MBSFace* givenFace, rg_dList< MBSFace* >& resultIncidentFaces );
                        void findAllVerticesHavingIdenticalGeometry( MBSVertex* givenVertex, rg_dList< MBSVertex* >& resultIdenticalGeometryVertices );

        void removeMinorEdgeConnectedFaceGroup();
            void removeMinorEdgeConnectedFaceGroupForEachPocketPrimitive( PocketPrimitive* currPocketPrimitive );
                rg_REAL computeSumOfAreaOnGivenFaceList( rg_dList< MBSFace* >* currFaceList );
                    rg_REAL computeAreaOfMBSFace( MBSFace* givenBetaFace );
        void enlargePocketPrimitives();
            void enlargePocketPrimitive( PocketPrimitive* currPocketPrimitive );

            
            
    void reportAtomsOfPockets( const char* filePath, const int& numOfPockets );

    
    
    
    
    rg_INT computeNumOfVoids();
    
        
       
    
 
};

#endif