#ifndef VORONOIDIAGRAMCIC_H
#define VORONOIDIAGRAMCIC_H

#include "VoronoiDiagram2DC.h"
#include "rg_Triplet.h"

#include <list>
#include <fstream>
using namespace std;


namespace V {
namespace GeometryTier {

enum CONTAINER_REDUCTION_CONSTRAINT
{
    STRONG_CONSTRAINT_INCLUDED,
    WEEK_CONSTRAINT_AT_LEAST_INTERSECTED
};

class VoronoiDiagramCIC : public VoronoiDiagram2DC
{
protected: 
    rg_Circle2D         m_container;
    list<rg_Circle2D>   m_disks;
    Generator2D*        m_containerGenerator;
	
public:
    //BEGIN - add constructor, operator= : [2016. 11.10. cysong]
    VoronoiDiagramCIC();
    VoronoiDiagramCIC(const VoronoiDiagramCIC& CIC);
    virtual ~VoronoiDiagramCIC();
    
    VoronoiDiagramCIC& operator=(const VoronoiDiagramCIC& CIC);

    void        clear();

    Generator2D*    insertDiskIntoVD( const rg_Circle2D& disk );
    Generator2D*    insertDiskIntoVD( pair< rg_Circle2D, void* >& diskNUserDataPair);
    Generator2D*    insertDiskIntoVD( const rg_Triplet< rg_Circle2D, int, void* > & diskAndIDAndUserData);
    
    bool            isThisDiskInContainer( const rg_Circle2D& disk );
    void            expandRadiusOfContainer( const rg_Circle2D& disk );
    void        clear_except_current_container();
    //END - add constructor, operator= : [2016. 11.10. cysong]

    void        constructVoronoiDiagramCIC( list<rg_Circle2D>& circleset );
    void        constructVoronoiDiagramCIC( list< pair< rg_Circle2D, void* > >& circleNUserDataPairSet );
    void        constructVoronoiDiagramCIC( list< rg_Triplet<rg_Circle2D, int, void*> >& circleAndIDAndUserDataSet );


    void        constructVoronoiDiagramCIC_noContainerInInput( list<rg_Circle2D>& circleset_without_container );
    //////////////////////////////////////////////////////////////////////////////////////////
    // For big VD
    
    void        constructVoronoiDiagramCIC( list<Generator2D*>& generators ); //for big VD
    void        constructVoronoiDiagramCIC( Generator2D* universeGenerator ); // for big VD
    void        constructVoronoiDiagramCIC_hierarchy( list< rg_Circle2D >& hierarchicalDisks ); // for big VD
    void        constructVoronoiDiagramCIC_hierarchy( list< pair< rg_Circle2D, void* > >& hierarchicalDisks ); // for big VD
    
    //
    //////////////////////////////////////////////////////////////////////////////////////////

    void        updateVoronoiDiagram_insertion( Generator2D* const generator );
    void        updateVoronoiDiagram_removal( Generator2D* const generator );
    void        deleteContainer_1f();

    //////////////////////////////////////////////////////////////////////////////////////////
    // For big VD

    void        updateVoronoiDiagram_insertion_addToGeneratorList( Generator2D* const generator, list<VEdge2D*>& newEdges, list<VEdge2D*>& removedEdges, list<VEdge2D*>& cuttedEdges );
    void        updateVoronoiDiagram_removal_deleteFromGeneratorList( Generator2D* const generator, list<VEdge2D*>& newEdges, list<VEdge2D*>& removedEdges, list<VEdge2D*>& influencedEdges );

    void        getAllVertices_hierarchy( list<VVertex2D*>& vertices ) const;
    void        getAllEdges_hierarchy( list<VEdge2D*>& edges ) const;
    void        getAllFaces_hierarchy( list<VFace2D*>& faces ) const;
    void        getAllGenerators_hierarchy( list<Generator2D*>& generators ) const;
    void        getAllVoronoiDiagrams_hierarchy( list<VoronoiDiagramCIC*>& VDs, const bool& includingThisVDItself = true) const;

    void        updateAllVEdges_hierarchy();

    //
    //////////////////////////////////////////////////////////////////////////////////////////
    // For VVertex computation
    void        updateCircumcirclesOfNewVVertices(list<VVertex2D*>& newVVertices);
	
    
    void        nutcracking( Generator2D* generator );

    bool        thisEdgeIsIntersectedByOffsetDisk(VEdge2D* edge, const double betaValue);

    inline const Generator2D* getContainerGenerator() const { return m_generators.front(); } ;
    inline void               setContainerGenerator(Generator2D* containerGenerator) { m_containerGenerator = containerGenerator; };
	bool        changeRadiusOfContainer( const double& targetRadius );
    void        printGenerators( );
	void        printGenerators(const string& fname);
    void        writeTopology(ofstream& fout);

    VFace2D*    findVFaceContainingQueryPoint( const rg_Point2D& pt ) const;
	
	inline void get_generators_without_container(list<Generator2D*>& diskGenerators) const // by Joonghyun on November 07, 2017
	{
		list<VFace2D*>::const_iterator i_face = m_VFaces.begin();
		for (; i_face != m_VFaces.end(); ++i_face)
		{
			if ((*i_face)->isInfinite())
				continue;

			Generator2D* currGenerator = (*i_face)->getGenerator();

			if(currGenerator != m_containerGenerator)
				diskGenerators.push_back(currGenerator);
		}
	}

	inline void get_generators(list<Generator2D*>& diskGenerators) const // by Joonghyun on October 28, 2017
	{
		list<VFace2D*>::const_iterator i_face = m_VFaces.begin();
		for (; i_face != m_VFaces.end(); ++i_face) 
		{
			if ((*i_face)->isInfinite())
				continue;

			Generator2D* currGenerator = (*i_face)->getGenerator();
			diskGenerators.push_back(currGenerator);
		}
	}

	inline rg_Circle2D getContainer() const // by Joonghyun on October 28, 2017
	{
		return m_generators.front()->getDisk();
	}

	inline Generator2D* getGeneratorOfContainer() // by Joonghyun on October 28, 2017
	{
		return m_generators.front();
	}

	inline void getGeneratorsWhichShareVEdgesWithContainer(list<Generator2D*>& generatorsSharingVEdgesWithContainer) // by Joonghyun on April 9, 2018
	{
		//VFace2D* VFaceOfContainer = m_generators.front()->getOuterFace();
        VFace2D* VFaceOfContainer = m_generators.front()->getInnerFace();
		list<VEdge2D*> boundingVEdgesForVFaceOfContinaer;
		VFaceOfContainer->getBoundaryVEdges(boundingVEdgesForVFaceOfContinaer);

		list<VEdge2D*>::iterator i_VEdge = boundingVEdgesForVFaceOfContinaer.begin();
		for (; i_VEdge != boundingVEdgesForVFaceOfContinaer.end(); ++i_VEdge)
		{
			VEdge2D* currVEdge = (*i_VEdge);
			Generator2D* currDiskGen = getFaceOfOppositeSide(VFaceOfContainer, currVEdge)->getGenerator();
			generatorsSharingVEdgesWithContainer.push_back(currDiskGen);
		}
	}

protected: //[2016.09.15 - cysong] private -> protected
    
    //BEGIN - copy, clear function : [2016. 11.10. cysong]

    void copyFrom(const VoronoiDiagramCIC& VD_CIC);

    void        getAllVEntities(const VoronoiDiagramCIC& VD_CIC,
		                        list<VoronoiDiagramCIC*>& allVDs_Input,
	                            list<Generator2D*>&       allGeneratorsInHierarcicalVD_Input,
	                            list<VFace2D*>&           allFacesInHierarcicalVD_Input,
	                            list<VEdge2D*>&           allEdgesInHierarcicalVD_Input,
	                            list<VVertex2D*>&         allVerticesInHierarcicalVD_Input) const;
  
    void        createNLinkVDsAndThoseEntities(const list<VoronoiDiagramCIC*>& allVDs_Input, VoronoiDiagramCIC* outerMostVD,
                                               const list<Generator2D*>& allGeneratorsInHierarcicalVD_Input,
                                               const list<VFace2D*>& allFacesInHierarcicalVD_Input,
                                               const list<VEdge2D*>& allEdgesInHierarcicalVD_Input,
                                               const list<VVertex2D*>& allVerticesInHierarcicalVD_Input,
                                               unordered_map<VoronoiDiagramCIC*, VoronoiDiagramCIC*>& VDMapFromOldToNew,
                                               unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew,
                                               unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew,
                                               unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew,
                                               unordered_map<VVertex2D*, VVertex2D*>& vertexMapFromOldToNew) const;


    void        createNLinkVoronoiDiagrams(const list<VoronoiDiagramCIC*>& allVDs, VoronoiDiagramCIC* outerMostVD, unordered_map<VoronoiDiagramCIC*, VoronoiDiagramCIC*>& VDMapFromOldToNew) const;
   
    void        copyTopologyUnrelatedObjects(const VoronoiDiagramCIC& VD_CIC);
   
    void        insertCreatedEntitiesToVDs(
                        const list<VoronoiDiagramCIC*>& allVDs_Input,
                        const unordered_map<VoronoiDiagramCIC*, VoronoiDiagramCIC*>& VDMapFromOldToNew,
                        const unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew,
                        const unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew,
                        const unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew,
                        const unordered_map<VVertex2D*, VVertex2D*>& vertexMapFromOldToNew);
    
    void        connectCreatedEntities(const list<Generator2D*>&       allGeneratorsInHierarcicalVD_Input,
                                       const list<VFace2D*>&           allFacesInHierarcicalVD_Input,
                                       const list<VEdge2D*>&           allEdgesInHierarcicalVD_Input,
                                       const list<VVertex2D*>&         allVerticesInHierarcicalVD_Input,
                                       const unordered_map<VoronoiDiagramCIC*, VoronoiDiagramCIC*>& VDMapFromOldToNew,
                                       const unordered_map<Generator2D*, Generator2D*>&             generatorMapFromOldToNew,
                                       const unordered_map<VFace2D*, VFace2D*>&                     faceMapFromOldToNew,
                                       const unordered_map<VEdge2D*, VEdge2D*>&                     edgeMapFromOldToNew,
                                       const unordered_map<VVertex2D*, VVertex2D*>&                 vertexMapFromOldToNew) const;

    void        connectCreatedEntitiesInGenerators(
                        const list<Generator2D*>& allGeneratorsInHierarcicalVD_Input,
		                const unordered_map<VoronoiDiagramCIC*, VoronoiDiagramCIC*>& VDMapFromOldToNew,
		                const unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew,
		                const unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew) const;
   

    //END - copy, clear function : [2016. 11.10. cysong]



    void        constructVoronoiDiagramCIC_TOI_Method();
    void        createGenerators( list<rg_Circle2D>& circleset );
    void        createGenerators( list< pair< rg_Circle2D, void* > >& circleNUserDataPairSet );
    void        createGenerators( list< rg_Triplet<rg_Circle2D, int, void*> >& circleAndIDAndUserDataSet );
    void        createGenerators_DoublePrecision(list<rg_Circle2D>& circleset);
    void        createGenerators_DoublePrecision(list< pair< rg_Circle2D, void* > >& circleNUserDataPairSet);
    void        createGenerators_DoublePrecision(list< rg_Triplet<rg_Circle2D, int, void*> >& circleAndIDAndUserDataSet);

    //////////////////////////////////////////////////////////////////////////////////////////
    // For big VD

    void        createGenerators_hierarchy( list< pair< rg_Circle2D, void* > >& hierarchicalDisks );
    void        createGenerators_hierarchy( list< rg_Circle2D >& hierarchicalDisks );
    void        collectContainerGenerators( list<Generator2D*>& containers ) const;

    //
    //////////////////////////////////////////////////////////////////////////////////////////

    void        constructSeedVoronoiDiagram();
    void        updateVoronoiDiagramCIC_TOI_Method( Generator2D* const generator );
    
    bool        isAnomalizingEdge( VEdge2D* const incidentEdge, Generator2D* const newGenerator );
    void        calculateTangentCircles( const rg_Circle2D& enclosingCircle, const rg_Circle2D& circle1, const rg_Circle2D& circle2, rg_Circle2D& result1, rg_Circle2D& result2 );
    rg_Circle2D solution1( const double& x1, const double& r1,
                           const double& x2, const double& y2, const double& r2 );
    rg_Circle2D solution2( const double& x1, const double& r1,
                           const double& x2, const double& y2, const double& r2 );

    void        computeCoordOfNewVVertices( list<VVertex2D*>& newVVertices );
    virtual void updateCircumcircle( VVertex2D* const vertex );
    void        getThreeGeneratorsDefiningVertex( VVertex2D*    tVertex,
                                                  Generator2D*& generator1,
                                                  Generator2D*& generator2,
                                                  Generator2D*& generator3 );
    rg_Point2D  getPassingPtOnEdge( VFace2D* lFace, VFace2D* rFace, VVertex2D* sVertex, VVertex2D* eVertex );
	int			computeCircumcirclesForDeletion(VFace2D* prevFace, VFace2D* currFace, VFace2D* nextFace, rg_Circle2D& circumCircle1, rg_Circle2D& circumCircle2);
	bool        isCorrectCircumcircle(rg_Circle2D& circumcircle, VFace2D* prevFace, VFace2D* currFace, VFace2D* nextFace);





    bool        increaseRadiusOfContainer( Generator2D* increasedContainer );
    void        collectPossiblyFlippingEdgesOnBoundaryEdgesOfContainer_increasing( Generator2D* container, rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges );
    bool        isThisEdgeFlippingEdge_increasingContainer( VEdge2D* edge, VFace2D* faceOfContainer, rg_Circle2D& circumcircle );
    bool        isThisCircumcircleInContainer( const rg_Circle2D& circumcircle, Generator2D* container );
    void        flipThisEdgeNFindTwoIncidentTriplets_increasing( VFace2D* faceOfContainer, pair<VEdge2D*, rg_Circle2D>& edgeNCircle, VEdge2D*& incidentTriplet1, VEdge2D*& incidentTriplet2 );
    void        reflectTheFlipOnTheIncidentTriplet_increasingContainer( rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges, VEdge2D*& edgeOfIncidentTriplet, VFace2D*& faceOfContainer );


    bool        reduceRadiusOfContainer( Generator2D* reducedContainer );
    bool        ThisIsPossibleContainerReduction( Generator2D* reducedContainer, const CONTAINER_REDUCTION_CONSTRAINT& constraint );
    void        collectPossiblyFlippingEdgesOnQuillEdgesOfContainer_reducing( Generator2D* container, rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges );
    bool        isThisEdgeFlippingEdge_reducingContainer( VEdge2D* edge, VFace2D* faceOfContainer, rg_Circle2D& circumcircle );
    void        flipThisEdgeNFindTwoIncidentQuillEdges( VFace2D* faceOfContainer, VEdge2D* edgeToBeFlipped, VEdge2D*& incidentQuillEdge1, VEdge2D*& incidentQuillEdge2 );
    void        reflectTheFlipOnTheIncidentQuillEdge_reducingContainer( rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges, VEdge2D*& edgeOfIncidentTriplet, VFace2D*& faceOfContainer );



    bool        flippingTest_outside( VEdge2D* const outerBoundaryEdge, VFace2D* const faceToBeRemoved, list<Generator2D*>& outerNeighbor, rg_Circle2D& circumcircle );
    bool        flippingTest_inside(  VEdge2D* const innerBoundaryEdge, VFace2D* const faceToBeRemoved, list<Generator2D*>& innerNeighbor, rg_Circle2D& circumcircle );
    void        reflectTheFlipOnTheIncidentTriplet_nutBreaking( rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges, VEdge2D*& edgeOfIncidentTriplet, list<Generator2D*>& outerNeighbor, VFace2D* faceToBeRemoved );
    Generator2D* findGeneratorWhichBlockThisEdge( VEdge2D* const edge, VFace2D* const nutFace, list<Generator2D*>& boundaryGen, rg_Circle2D& circumcircle );

};

} // GeometryTier
} // V

#endif


