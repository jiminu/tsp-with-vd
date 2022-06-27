#ifndef _AUGMENTEDBETACOMPLEX_H
#define _AUGMENTEDBETACOMPLEX_H

#include "rg_dList.h"
#include "BetaUniverse.h"

#include "AugmentedBCVertex.h"
#include "AugmentedBCEdge.h"
#include "AugmentedBCFace.h"
#include "AugmentedBCCell.h"
#include "Ball.h"
#include "PointerToSimplex.h"
#include "NonManifoldCorner.h"

#include "ManifoldizedBetaShape.h"

#include <map>
#include <fstream>
using namespace std;



const rg_INT FALSE_EXTERIOR_SIMPLEX = -2;


class AugmentedBetaComplex
{
private:
    BetaUniverse*                   m_betaUniverse;
	rg_REAL                         m_betaValue;

	rg_dList< AugmentedBCCell >		m_augmentedBCCells;
	rg_dList< AugmentedBCFace >		m_augmentedBCFaces;
	rg_dList< AugmentedBCEdge >		m_augmentedBCEdges;
	rg_dList< AugmentedBCVertex >	m_augmentedBCVertices;

   	rg_dList< AugmentedBCCell >		m_artificialAugmentedBCCells;
	rg_dList< AugmentedBCFace >		m_artificialAugmentedBCFaces;
	rg_dList< AugmentedBCEdge >		m_artificialAugmentedBCEdges;
	rg_dList< AugmentedBCVertex >	m_artificialAugmentedBCVertices;

    rg_dList<Ball>					m_balls;


	rg_INT*		r_numExtraneousCellGroupAtVertex;
	//ofstream	r_fout;
	rg_INT      r_numVerticesDifferentNumConesInputAndTT_D_V;
	
public:
	AugmentedBetaComplex();
	~AugmentedBetaComplex();

	///////////////////////////////////////////////////////////////////////////
	//   
    //  GET functions
    BetaUniverse* getBetaUniverse() const;
	rg_REAL	getBetaValue() const;

    rg_INT  getNumAugmentedBCCells() const;
	rg_INT  getNumAugmentedBCFaces() const;
	rg_INT  getNumAugmentedBCEdges() const;
	rg_INT  getNumAugmentedBCVertices() const;
    
    rg_INT  getNumArtificialAugmentedBCCells() const;
	rg_INT  getNumArtificialAugmentedBCFaces() const;
	rg_INT  getNumArtificialAugmentedBCEdges() const;
	rg_INT  getNumArtificialAugmentedBCVertices() const;

    rg_INT  getNumAllAugmentedBCCells() const;
	rg_INT  getNumAllAugmentedBCFaces() const;
	rg_INT  getNumAllAugmentedBCEdges() const;
	rg_INT  getNumAllAugmentedBCVertices() const;    
    
    rg_dList<AugmentedBCCell>*		getAugmentedBCCells();
	rg_dList<AugmentedBCFace>*      getAugmentedBCFaces();
	rg_dList<AugmentedBCEdge>*      getAugmentedBCEdges();
	rg_dList<AugmentedBCVertex>*    getAugmentedBCVertices();

    rg_dList<AugmentedBCCell>*		getArtificialAugmentedBCCells();
	rg_dList<AugmentedBCFace>*      getArtificialAugmentedBCFaces();
	rg_dList<AugmentedBCEdge>*      getArtificialAugmentedBCEdges();
	rg_dList<AugmentedBCVertex>*    getArtificialAugmentedBCVertices();	

    rg_dList<Ball>* getBalls();

    rg_INT getMaxIndexOfAugmentedBCVertices();
    rg_INT getMaxIndexOfAugmentedBCEdges();
    rg_INT getMaxIndexOfAugmentedBCFaces();
    rg_INT getMaxIndexOfAugmentedBCCells();



	///////////////////////////////////////////////////////////////////////////
	//   
    //  SET functions    
    void setBetaUniverse( BetaUniverse* betaUniverse );
	void setBetaValue(const rg_REAL& value);


    void clearMyself();


	///////////////////////////////////////////////////////////////////////////
	//   
    //  MAKE augmented beta-complex
	//void construct( const rg_REAL& betaValue, BetaUniverse* betaUniverse ); 

    //  The following member function is fixed by Youngsong Cho 2011-04-23
    //void construct( const rg_REAL& betaValue, BetaUniverse* betaUniverse, const rg_INT& deleteOption = CONSTRUCT_WITH_WHOLE_BETA_SHAPE );
    void construct( const rg_REAL& betaValue, BetaUniverse* betaUniverse, const rg_INT& constructOption = CONSTRUCT_WITH_WHOLE_BETA_SHAPE );

	void extractBoundaryOfAugmentedBetaComplex( ManifoldizedBetaShape& manifoldizedBS);

    void adjustIDOfAugmentedBetaComplex();
    void adjustIDOfAugmentedBCVertices();
    void adjustIDOfAugmentedBCEdges();
    void adjustIDOfAugmentedBCFaces();
    void adjustIDOfAugmentedBCCells();
	void adjustIDOfBalls();

    


    void reportResultOfManifoldization(const double& totalTime);
    void checkCopy( BetaUniverse& BC );
    void checkAllEdgesIncidentTo2Faces(ofstream& fout);

    void testTheIDOfAugmentedBetaComplex( BetaUniverse& BC );    


private:
    void duplicateOriginalQuasiTriangulationForBetaToAugmentedBetaComplex( BetaUniverse* betaUniverse );
        void duplicateOriginalBetaVerticesToAugmentedBCVertices(rg_dList<BetaVertex>* betaVertices,
																AugmentedBCVertex** ptOfAugVertices,
																AugmentedBCEdge**   ptOfAugEdges,
																AugmentedBCCell**   ptOfAugCells);
        void duplicateOrigianlBetaEdgesToAugmentedBCEdges( rg_dList<BetaEdge>* betaEdges, 
														   AugmentedBCVertex** ptOfAugVertices, 
														   AugmentedBCEdge**   ptOfAugEdges, 
														   AugmentedBCFace**   ptOfAugFaces);
        void duplicateOriginalBetaFacesToAugmentedBCFaces( rg_dList<BetaFace>* betaFaces, 
                                                 AugmentedBCEdge**   ptOfAugEdgesToBetaEdges, 
                                                 AugmentedBCFace**   ptOfAugFacesToBetaFaces, 
                                                 AugmentedBCCell**   ptOfAugCellsToBetaCells);
        void duplicateOriginalBetaCellsToAugmentedBCCells( rg_dList<BetaCell>* betaCells,
												 AugmentedBCVertex** ptOfAugVertices,
                                                 AugmentedBCFace**   ptOfAugFacesToBetaFaces, 
                                                 AugmentedBCCell**   ptOfAugCellsToBetaCells);


    //////////////////////////////////////////////////////////////////////////////////////
    //  Youngsong Cho 2011-04-22 
    void duplicateBetaComplexToAugmentedBetaComplex();
        void createAugmentedBetaSimplexesAndRegisterToMap(  map<BetaVertex*, AugmentedBCVertex*>& mapBetaVertexToAugmentedBetaVertex,
                                                            map<BetaEdge*,   AugmentedBCEdge*>&   mapBetaEdgeToAugmentedBetaEdge,
                                                            map<BetaFace*,   AugmentedBCFace*>&   mapBetaFaceToAugmentedBetaFace,
                                                            map<BetaCell*,   AugmentedBCCell*>&   mapBetaCellToAugmentedBetaCell );

        void duplicateBetaVerticesToAugmentedBCVertices(    const map<BetaVertex*, AugmentedBCVertex*>& mapBetaVertexToAugmentedBetaVertex,
                                                            const map<BetaEdge*,   AugmentedBCEdge*>&   mapBetaEdgeToAugmentedBetaEdge,
                                                            const map<BetaCell*,   AugmentedBCCell*>&   mapBetaCellToAugmentedBetaCell );

        void duplicateBetaEdgesToAugmentedBCEdges(          const map<BetaVertex*, AugmentedBCVertex*>& mapBetaVertexToAugmentedBetaVertex,
                                                            const map<BetaEdge*,   AugmentedBCEdge*>&   mapBetaEdgeToAugmentedBetaEdge,
                                                            const map<BetaFace*,   AugmentedBCFace*>&   mapBetaFaceToAugmentedBetaFace);
    
        void duplicateBetaFacesToAugmentedBCFaces(          const map<BetaEdge*,   AugmentedBCEdge*>&   mapBetaEdgeToAugmentedBetaEdge,
                                                            const map<BetaFace*,   AugmentedBCFace*>&   mapBetaFaceToAugmentedBetaFace,
                                                            const map<BetaCell*,   AugmentedBCCell*>&   mapBetaCellToAugmentedBetaCell );

        void duplicateBetaCellsToAugmentedBCCells(          const map<BetaVertex*, AugmentedBCVertex*>& mapBetaVertexToAugmentedBetaVertex,
                                                            const map<BetaFace*,   AugmentedBCFace*>&   mapBetaFaceToAugmentedBetaFace,
                                                            const map<BetaCell*,   AugmentedBCCell*>&   mapBetaCellToAugmentedBetaCell );

        void writeDuplicate();
        void writeNonManifoldCorner( rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges );
    //
    //////////////////////////////////////////////////////////////////////////////////////

    void setEdgeAsSingularIfEdgeIncidentToExtraneousFacesOnly();
		

    void indexExteriorRegionOfEdgesAndFaces();



    void searchNonManifoldConfigurationsForAllVertices( rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMVertices, 
                                                        rg_dList< NonManifoldCorner >& NonManifoldCornerList) const;

        void generateNonManifoldVertexConfiguration(AugmentedBCVertex*             currVertex, 
                                                    rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMVertices,
                                                    rg_dList< NonManifoldCorner >& NonManifoldCornerList)  const;
            void findSingularEdgesIncidentToVertex( AugmentedBCVertex*  givenVertex, rg_dList< BetaEdge* >& singularEdgeList ) const;
            void findSingularFacesIncidentToVertex( AugmentedBCVertex*  givenVertex, rg_dList< BetaFace* >& singularFaceList ) const;
            void searchGroupsOfFaceConnectedExtraneousCellsInBetaComplex( AugmentedBCVertex* givenVertex, rg_dList< rg_dList<BetaCell*> >& groupsOfExtraneousCell ) const;
            void searchGroupsOfFaceConnectedInteriorCellsInBetaComplex( AugmentedBCVertex* givenVertex, rg_dList< rg_dList<BetaCell*> >& groupsOfInteriorCell ) const;
	
			void generateNMVertexCornersInEachExtraneousCellGroup( AugmentedBCVertex*               givenVertex, 
                                                                   rg_dList< rg_dList<BetaCell*> >& extraneousCellGroupList, 
																   rg_dList<NonManifoldCorner*>*    ptToNMCornerOfNMVertices,
																   rg_dList< NonManifoldCorner >&   nonManifoldCornerList ) const;
				void findSingularEdgesIncidentToVertexInCellGroup( AugmentedBCVertex* givenVertex, rg_dList< BetaCell* >* givenCellGroup, rg_dList< BetaEdge* >& singularEdgesIncidentToGivenVertex) const;
				void findBetaShapeFacesIncidentToVertexInCellGroup(AugmentedBCVertex* givenVertex, rg_dList< BetaCell* >* givenCellGroup, rg_dList< BetaFace* >& singularFacesIncidentToGivenVertex) const;				
				void separateFacesIntoEdgeConnectedComponent( AugmentedBCVertex* givenVertex, rg_dList < BetaFace* >& facesToBeSeparatedByEdgeConnection, rg_dList< rg_dList < BetaFace* > >& edgeConnectedFaceGroupList ) const;
					rg_BOOL isTheseFacesConnectedWithEdgesWhichBoundedByGivenVertex(AugmentedBCVertex* givenVertex, BetaFace* sourceFace, BetaFace* targetFace) const;
				void generateNMVertexCornersWithSingularEdges( AugmentedBCVertex*             givenVertex, 
                                                               rg_dList<BetaEdge*>&           singularEdgeList,             
                                                               rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMVertices,
															   rg_dList< NonManifoldCorner >& nonManifoldCornerList ) const;
				void generateNMVertexCornersWithEdgeConnectedFaceGroup( AugmentedBCVertex*               givenVertex, 
                                                                        rg_dList< rg_dList<BetaFace*> >& edgeConnectedFaceList, 
                                                                        rg_dList<NonManifoldCorner*>*    ptToNMCornerOfNMVertices,
                                                                        rg_dList< NonManifoldCorner >&   nonManifoldCornerList ) const;
					BetaFace* findRegularFaceInGroupOfFaceList( rg_dList< rg_dList < BetaFace* > >& edgeConnectedFaceGroupList ) const;
					BetaFace* findRegularFaceInFaceList(rg_dList < BetaFace* >* faceList) const;
					void createPtOfSingularFaceOrCellBoundedByRegularFace( AugmentedBCFace* givenFace, PointerToSimplex& ptToSingularFaceOrCell ) const;
	

	void searchNonManifoldConfigurationsForAllEdges( rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMEdges, 
                                                     rg_dList< NonManifoldCorner >& NonManifoldCornerList);
        void generateNonManifoldEdgeConfiguration( AugmentedBCEdge*               curEdge, 
                                                   rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMEdges,
                                                   rg_dList< NonManifoldCorner >& NonManifoldCornerList);
            void generateNonManifoldEdgeConfigurationByOrderOfFacesWRTGivenEdge(AugmentedBCEdge*               givenEdge, 
                                                                                rg_dList< BetaFace* >&         incidentFaceList,
																				rg_dList<NonManifoldCorner*>*  ptToNMCornerOfNMEdges,
																				rg_dList< NonManifoldCorner >& NonManifoldCornerList) ;
                


                
    void detachAndRemoveExtraneousSimplexes();
        void setEdgeAsExtraneousIfEdgeIncidentToExtraneousFacesOnly();

        void detachExtraneousEdgesFromAllVertices();
		void detachExtraneousCellsFromAllVertices();
        void detachExtraneousFacesFromAllEdges();
        void detachExtraneousCellsFromAllFaces();
        
        void removeExtraneousCells();
        void removeExtraneousFaces();
        void removeExtraneousEdges();
        void removeExtraneousVertices();
    
    
    void augment(rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges,
                 rg_dList< NonManifoldCorner >& NonManifoldCornerList);
        void deletePTToNMCornerOfGivenNMCorner( NonManifoldCorner* givenNMCorner, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges,
                                                rg_dList< NonManifoldCorner >& NonManifoldCornerList  );
        NonManifoldCorner resolve_EE_V_Configuration( const NonManifoldCorner& givenNMConf, 
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices );
            void O_EE_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCEdge* edgeSimplex1, AugmentedBCEdge* edgeSimplex2, AugmentedBCEdge*& raisedNMSimplex, AugmentedBCFace*& raisedSimplex1, AugmentedBCFace*& raisedSimplex2);
        
        NonManifoldCorner resolve_EF_V_Configuration( const NonManifoldCorner& givenNMConf, 
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges );
            void O_EF_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCEdge* edgeSimplex1, AugmentedBCFace* edgeSimplex2, AugmentedBCEdge*& raisedNMSimplex, AugmentedBCFace*& raisedSimplex1, AugmentedBCCell*& raisedSimplex2);
        
        NonManifoldCorner resolve_ET_V_Configuration( const NonManifoldCorner& givenNMConf, 
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices );
            void O_ET_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCEdge* givenEdge, AugmentedBCFace* givenRegularFace, AugmentedBCEdge*& raisedNMEdge, AugmentedBCFace*& raisedFace, AugmentedBCCell*& raisedTetra);
        
        NonManifoldCorner resolve_FF_V_Configuration( const NonManifoldCorner& givenNMConf, 
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges );
            void O_FF_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCFace* givenFace1, AugmentedBCFace* givenFace2, AugmentedBCEdge*& raisedNMEdge, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2);
        
        NonManifoldCorner resolve_FT_V_Configuration( const NonManifoldCorner& givenNMConf, 
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges );
            void O_FT_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCFace* givenFace, AugmentedBCFace* givenRegularFace, AugmentedBCEdge*& raisedNMEdge, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2);
        
        NonManifoldCorner resolve_TT_V_Configuration( const NonManifoldCorner& givenNMConf );
            void O_TT_V_operator(AugmentedBCVertex* givenVertex, AugmentedBCFace* givenFirstRegularFace, AugmentedBCFace* givenSecondRegularFace, AugmentedBCEdge*& raisedNMEdge, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2);

        
        NonManifoldCorner resolve_FF_E_Configuration( const NonManifoldCorner& givenNMConf, 
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges );
            void O_FF_E_operator(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenFace1, AugmentedBCFace* givenFace2, AugmentedBCFace*& raisedNMEdge, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2);
        
        NonManifoldCorner resolve_FT_E_Configuration( const NonManifoldCorner& givenNMConf, 
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges);
            void O_FT_E_operator(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenFace, AugmentedBCFace* givenRegularFace, AugmentedBCFace*& raisedNMFace, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2);
        
        NonManifoldCorner resolve_TF_E_Configuration( const NonManifoldCorner& givenNMConf, 
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices,
                                                      rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges);
            void O_TF_E_operator(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenRegularFace, AugmentedBCFace* givenFace, AugmentedBCFace*& raisedNMFace, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2);
        
        NonManifoldCorner resolve_TT_E_Configuration( const NonManifoldCorner& givenNMConf );
            void O_TT_E_operator(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenRegularFace1, AugmentedBCFace* givenRegularFace2, AugmentedBCFace*& raisedNMFace, AugmentedBCCell*& raisedTetra1, AugmentedBCCell*& raisedTetra2);

        
        NonManifoldCorner resolve_TT_D_V_Configuration( const NonManifoldCorner& givenNMConf );
			void O_TT_D_V_operator(AugmentedBCVertex* givenVertex);


        rg_INT findIndexOfCommonExteriorRegion( AugmentedBCFace* givenFace1, AugmentedBCFace* givenFace2 );

        void setIndexOfExteriorRegionOfSingularFace( const rg_INT& indexOfExteriorRegion, AugmentedBCFace* singularFace );
        void setIndexOfExteriorRegionOfCellExceptGivenFace( const rg_INT& indexOfExteriorRegion, AugmentedBCFace* exceptThisFace, AugmentedBCCell* givenCell );


        AugmentedBCFace* findRegularFaceTouchedExactExterior( AugmentedBCVertex* givenVertex,     AugmentedBCFace* givenRegularFace,             const rg_INT& indexOfExteriorRegion );
		AugmentedBCFace* findRegularFaceBoundedByGivenEdge( AugmentedBCEdge* givenEdge, AugmentedBCFace* givenRegularFace );
        
        
        void updateNMCornersForLowerDimensionSimplexOfGivenSimplex(AugmentedBCEdge* givenEdge, AugmentedBCFace* raisedFace, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices);
        void updateNMCornersForLowerDimensionSimplexOfGivenSimplex(AugmentedBCFace* givenFace, AugmentedBCCell* raisedTetra, 
                                                                   rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMVertices, rg_dList<NonManifoldCorner*>* ptToNMCornerOfNMEdges);
            void replaceGivenSimplexToRaisedSimplex( void* original, void* substitute, rg_dList<NonManifoldCorner*>* currPTToNMConfOfNMVertices );
        
        AugmentedBCVertex* replicateGivenVertex(AugmentedBCVertex* originalVertex);
        AugmentedBCEdge*   addArtificialEdge(AugmentedBCVertex* givenVertex, AugmentedBCVertex* artificialVertex);
        AugmentedBCFace*   addArtificialFaceBy2Edges(AugmentedBCEdge* e1, AugmentedBCEdge* e2, AugmentedBCVertex* commonVertex);
        AugmentedBCFace*   addArtificialFaceBy3Edges(AugmentedBCEdge* e1, AugmentedBCEdge* e2, AugmentedBCEdge* e3);
        AugmentedBCCell*   addArtificialCellBy1EdgeAnd1Face(AugmentedBCEdge* givenEdge, AugmentedBCFace* givenFace, AugmentedBCVertex* commonVertex);
        AugmentedBCCell*   addArtificialCellByRegularAndSingularFaces(AugmentedBCFace* givenFirstFace, AugmentedBCFace* givenSecondFace, AugmentedBCEdge* commonEdge);

        AugmentedBCFace* findIncidentRegularFaceOfGivenCellAndVertex(AugmentedBCCell* givenCell, AugmentedBCVertex* givenVertex);
        AugmentedBCFace* findRegularFaceOfGivenCellByBackwardOrientation(AugmentedBCFace* givenRegularFace, AugmentedBCEdge* givenEdge);
        AugmentedBCFace* findRegularFaceOfGivenCellByForwardOrientation(AugmentedBCFace* givenRegularFace, AugmentedBCEdge* givenEdge);

        void addArtificialCellsByGivenVertexAndFaceGroup(AugmentedBCVertex* givenVertex, rg_dList<AugmentedBCFace*>* faceGroup);
            void findBridgeFaceByIncidentFaceAndCellNotThisBoundingFace(AugmentedBCVertex* givenVertex, AugmentedBCFace* incidentFace, AugmentedBCCell* givenCell, AugmentedBCFace* noThisBoundingFace,            
                                                                        AugmentedBCEdge*& outputCommonEdge, AugmentedBCFace*& outputBoundingFace);
            void setClassOfAugmentedBCEdgesToInteriorIfEdgesBoundRegularFaces(AugmentedBCVertex* givenVertex, rg_dList<AugmentedBCFace*>* faceGroup);

        void fixClassOfEdgesAndVertices();

		void groupFacesIncidentToVertexByEdge(AugmentedBCVertex* givenVertex, 
                                              rg_dList< rg_dList<AugmentedBCFace*> >& faceGroupList);



    
	//extract the boundary of augmented beta-complex.
    void initializePointers( MBSVertex**& ptMBSVertices,       MBSEdge**& ptMBSEdges, MBSFace**& ptMBSFaces );
    void duplicateBetaShapeSimplexesOfAugBCtoMBS( MBSVertex** ptMBSVertices, MBSEdge** ptMBSEdges, MBSFace** ptMBSFaces,
                                                    rg_dList< MBSVertex* >& mbsVertices, rg_dList< MBSEdge* >& mbsEdges, rg_dList< MBSFace* >& mbsFaces );
        void duplicateVerticesOfBetaShapeInAugBCtoMBS( MBSVertex** ptMBSVertices, rg_dList< MBSVertex* >& mbsVertices );
        void duplicateEdgesOfBetaShapeInAugBCtoMBS( MBSVertex** ptMBSVertices, MBSEdge** ptMBSEdges,            rg_dList< MBSEdge* >& mbsEdges );
        void duplicateFacesOfBetaShapeInAugBCtoMBS( MBSEdge** ptMBSEdges, MBSFace** ptMBSFaces,            rg_dList< MBSFace* >& mbsFaces );
    void buildTopologyOfMBSEdges( MBSVertex** ptMBSVertices, MBSEdge** ptMBSEdges, MBSFace** ptMBSFaces );
        void buildTopologyOfMBSEdge( AugmentedBCEdge* currAugmentedBCEdge, MBSVertex** ptMBSVertices, MBSEdge** ptMBSEdges, MBSFace** ptMBSFaces );
    void divideMBSSimplexesIntoShells( rg_dList< MBSVertex* >& mbsVertices, rg_dList< MBSEdge* >& mbsEdges, rg_dList< MBSFace* >& mbsFaces,
                                       rg_dList< MBSShell* >& mbsShells );
    void divideMBSShellsIntoBodies( rg_dList< MBSShell* > mbsShells, ManifoldizedBetaShape& manifoldizedBS );
        rg_INT markToVerticesWhichBodyTheyBelongAndReturnNumberOfBodies( rg_INT* mbsVerticesBelongToThisBody );
    
    void divideExteriorShellFromInteriorShells( ManifoldizedBetaShape& manifoldizedBS );
        void searchVerticesWhichBoundExteriorShell( rg_dList< BetaVertex* >& verticesBoundExteriorShell );
		

    void findAugRegularFacesIncidentToAugEdge( AugmentedBCEdge* givenAugEdge, rg_dList< AugmentedBCFace* >& incidentFaces );
	
		void findAugInteriorCellsIncidentToAugVertex( AugmentedBCVertex* givenAugVertex, rg_dList< AugmentedBCCell* >& incidentCells );
        void findAugRegularFacesIncidentToAugVertex( AugmentedBCVertex* givenAugVertex, rg_dList< AugmentedBCFace* >& incidentFaces );
public:
	void checkComparingNumExtraneousCellGroupAndNumRegularFaceGroup();
	void reportNumNMCorners(rg_dList< NonManifoldCorner >& nmCorners);
	rg_INT getNumVerticesDifferentNumConesInputAndTT_D_V();
	rg_INT getNumVerticesOutputConeIsNotOne();
	
	
	rg_REAL findUpperValueOfBoundingState( const rg_BetaSpan& betaSpan, rg_INT givenBoundingState);
	
    
    
    
    
    
    
    AugmentedBCFace* createSingularEllipticFace( AugmentedBCEdge* givenEdge, AugmentedBCFace* singularEllipticFace );
        void fixBoundingStateOfBoundingEdgesOfFace( AugmentedBCFace* givenFace );
        void fixBoundingStateOfBoundingVerticesOfFace( AugmentedBCFace* givenFace );
        
        
        
        
        
    void searchFacesIncidentToEdgeInBetaComplex( BetaEdge*        givenEdge, rg_dList< BetaFace* >& incidentFacesInBetaComplex )     const;
    void searchFacesIncidentToEdgeInBetaShape(   AugmentedBCEdge* givenEdge, rg_dList< BetaFace* >& facesIncidentToEdgeInBetaShape ) const;

    
    
    void printBetaCell(ofstream& fout, BetaCell* betaCell, const rg_REAL& beta);
    void printBetaFace(ofstream& fout, BetaFace* betaFace, const rg_REAL& beta);
    char getBoundingState(const rg_INT& boundingState);

    
    
    void sortFacesCCWOfEdgeByGeometricComparison( AugmentedBCEdge*  givenEdge, rg_dList< BetaFace* >& facesIncidentToEdge ) const;
        
        
    void setSingularEdgesToFalseExterior();
        void setSingularEdgeToFalseExterior( AugmentedBCEdge* currEdge );
    void setSingularFacesToFalseExterior();
    void setSingularVerticesToFalseExterior();
    
    
    
    
    void adjustIDOfQTForBeta( BetaUniverse* betaUniverse );
        void adjustIDOfQTForBetaVertices( BetaUniverse* betaUniverse );
        void adjustIDOfQTForBetaEdges( BetaUniverse* betaUniverse );
        void adjustIDOfQTForBetaFaces( BetaUniverse* betaUniverse );
        void adjustIDOfQTForBetaCells( BetaUniverse* betaUniverse );
        void adjustIDOfBallsInQTForBeta( BetaUniverse* betaUniverse );
        void printBetaEdge(ofstream& fout, BetaEdge* betaEdge, const rg_REAL& beta);
        
        
        
        
    
	
};

#endif

