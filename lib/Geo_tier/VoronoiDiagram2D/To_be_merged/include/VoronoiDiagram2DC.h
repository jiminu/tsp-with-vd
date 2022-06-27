#ifndef VORONOIDIAGRAM2DC_H
#define VORONOIDIAGRAM2DC_H


#include "VEdge2D.h"
#include "VFace2D.h"
#include "VVertex2D.h"
#include "Generator2D.h"

#include "rg_Point2D.h"
#include "rg_Circle2D.h"
#include "rg_BoundingBox2D.h"
#include "BucketForDisks.h"
//#include "Priority_Q_VEdge.h"

#include <set>
#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {

class Priority_Q_VEdge;

const bool IS_BUCKET_USED = true;
const bool IS_INLINE_USED = true;
const double RATIO_OF_CHANGING_TO_BUCKET = -1.0; // If you do not use bucket and Voronoi diagram structure simultaneously, set the value negative.

class VoronoiDiagram2DC
{
private:
    list<VEdge2D>       m_VEdges;
    list<VFace2D>       m_VFaces;
    list<VVertex2D>     m_VVertices;
    list<Generator2D>   m_generators;
    list<Generator2D>   m_phantomGenerators;
    list<rg_Circle2D>   m_phantomCircles;

    static const int    INFINITE_CIRCLE_RADIUS = 1e6;

    BucketForDisks      m_bucket;

    double time1, time2, time3, time4;
    double time2_1, time2_2, time2_3, time2_4, time2_5, time2_6;
    double time2_2_1, time2_2_2, time2_2_3, time2_2_4;
    double phantom1, phantom2, phantom3, phantom4, phantom5, phantom6;



public:
    VoronoiDiagram2DC();
    VoronoiDiagram2DC( list<rg_Circle2D>& circleset );
    //VoronoiDiagram2DC( const VoronoiDiagram2DC& VD2DC );
    ~VoronoiDiagram2DC();

    void         getVoronoiEdges(    list<const VEdge2D*>& VEdgesList ) const;
    void         getVoronoiFaces(    list<const VFace2D*>& VFacesList ) const;
    void         getVoronoiVertices( list<const VVertex2D*>& VVerticesList ) const;
    void         getGenerators(      list<const Generator2D*>& generatorsList ) const;
    void         getPhantomGenerator(list<const Generator2D*>& phantomGeneratorsList ) const;


    /////////////// for experiment ///////////////
    void         getFaces( list<VFace2D*>& faces );
    void         getVertices( list<VVertex2D*>& vertices );
    //////////////////////////////////////////////

    // For Dynamic Voronoi Diagram
    void         getBoundedVoronoiEdges( list<VEdge2D*>& boundedEdgeList );

    int          getNumberOfVVertices();
    int          getNumberOfVEdges();
    int          getNumberOfVFaces();

    void         createSortedGeneratorSet(      list<rg_Circle2D>& circleset );
    //void         setGenerators(      list<Generator2D>& generators );

    VoronoiDiagram2DC& operator=(    const VoronoiDiagram2DC& VD2DP );

    void         constructVoronoiDiagram( list<rg_Circle2D>& circleset );
    void         constructVoronoiDiagram();

    // For Dynamic Voronoi Diagram //
    void         constructVoronoiDiagramWithoutPhantomDiskRemoval( list<rg_Circle2D>& circleset );
    void         constructVoronoiDiagramWithoutPhantomDiskRemoval();
    /////////////////////////////////

    VFace2D*     findVFaceContainingQueryPoint( const rg_Point2D& pt ) ;
    void         clear();
    void         updateGeometry();
    void         printComputationTime( ofstream& fout );

private:
    VVertex2D*   createVVertex( const VVertex2D& vertex );
    VEdge2D*     createVEdge(   const VEdge2D& edge );
    VFace2D*     createVFace(   const VFace2D& face );

    void         insertPhantomGenerators();
    void         getBoundingBoxNMaxRadius( rg_BoundingBox2D& boundingBox, double& maxRadius );
    void         computeSugiharaCoordOfPhantomGenerators( const rg_BoundingBox2D& boundingBox,
                                                          rg_Point2D& firstCoord,
                                                          rg_Point2D& secondCoord,
                                                          rg_Point2D& thirdCoord );
    void         computeCoordOfPhantomGenerators_beforeOffset( const rg_BoundingBox2D& boundingBox,
                                                               rg_Point2D& firstCoord,
                                                               rg_Point2D& secondCoord,
                                                               rg_Point2D& thirdCoord );
    void         transformCoordinatesToMakeOffsetTriangleByMaxRadius( const double& maxRadius,
                                                                      rg_Point2D& firstCoord,
                                                                      rg_Point2D& secondCoord,
                                                                      rg_Point2D& thirdCoord );
    void         setPhantomGenerators( const double& maxRadius,
                                       rg_Point2D& firstCoord,
                                       rg_Point2D& secondCoord,
                                       rg_Point2D& thirdCoord );
    void         constructVoronoiDiagramForPhantomGenerators();
    void         updateVoronoiDiagramWithNewGenerator( Generator2D* const newGenerator );
    Generator2D* findClosestGeneratorToNewGenerator(   Generator2D* const newGenerator );
    void         removeAllExtraneousVVerticesAndVEdges();


    void         findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop_EdgeSplitPossible( Generator2D* const newGenerator,
                                                                                      Generator2D* const closestGenerator,
                                                                                      list<VEdge2D*>& intersectingVEdges,
                                                                                      list<VVertex2D*>& fictitiousVVertices);
    void         resetAnomalyTestDoneToggle( list<VEdge2D*>& anomalyTestedVEdges );
    void         findAnomalizingEdgeAmongIncidentVEdgesAndSplit( VVertex2D* const currRedVVertex,
                                                                 Generator2D* const newGenerator,
                                                                 list<VVertex2D*>& fictitiousVVertices,
                                                                 list<VEdge2D*>& anomalyTestedVEdges );
    void         findRedVVertexAmongNeighbor( VVertex2D* const currRedVVertex,
                                              Generator2D* const newGenerator,
                                              list<VVertex2D*>& blueVVertices,
                                              list<VVertex2D*>& redVVertices,
                                              list<VVertex2D*>& fictitiousVVertices,
                                              list<VEdge2D*>& anomalyTestedVEdges );
    void         findCrossingEdges( list<VVertex2D*>& redVVertices, list<VEdge2D*>& intersectingVEdges );
    void         reconfigureByConvertingBlueVVerticesToWhite( list<VVertex2D*>& blueVVertices );
    void         makeVEdgeLoopForNewGeneratorAndConnectToCurrVD(    Generator2D* const newGenerator,
                                                                    list<VEdge2D*>& intersectingVEdges,
                                                                    list<VVertex2D*>& newVVertices,
                                                                    list<VEdge2D*>& newVEdges );
    void         connectCurrVDToNewVEdgeLoop( list<VEdge2D*>& newVEdges );
    void         computeCoordOfNewVVertices(  list<VVertex2D*>& newVVertices );

    VVertex2D*   findFirstRedVertex( Generator2D* const newGenerator, Generator2D* const closestGenerator );
    bool         hasCycleOccurred(   VVertex2D* const candidateForRedVVertex );
    bool         hasOldRegionSplit(  VVertex2D* const candidateForRedVVertex );
    bool         isAnomalizingEdge(  VEdge2D* const incidentEdge, Generator2D* const newGenerator );
    void         mergeSplitVEdgesByFictitiousVVertex( list<VVertex2D*>& fictitousVVertices );
    rg_BOOL      areTwoVerticesDefinedBySameGeneratorTriplet( VVertex2D* const firstVertex, VVertex2D* const secondVertex );
    VVertex2D*   splitVEdgeAtFictitiousVVertex( VEdge2D* const currVEdge );
    void         updateEdge( VEdge2D* target );
    rg_Point2D   getTangentVector( VVertex2D* tVertex, VFace2D* lFace, VFace2D* rFace );
    rg_Point2D   getPassingPtOnEdge( VFace2D* lFace, VFace2D* rFace, VVertex2D* sVertex, VVertex2D* eVertex );
    void         wavePropagation_ver1( Generator2D* newGenerator, 
                                       list<VVertex2D*>& redVVertices,
                                       list<VVertex2D*>& blueVVertices,
                                       list<VVertex2D*>& fictitiousVVertices );
    void         wavePropagation_ver2( Generator2D* newGenerator, 
                                       list<VVertex2D*>& redVVertices,
                                       list<VVertex2D*>& blueVVertices,
                                       list<VVertex2D*>& fictitiousVVertices );


    static bool  compareCircleDescendingorder( const rg_Circle2D& circle1, const rg_Circle2D& circle2 );

private:
    void        removePhantomGenerators();
        void        collectEntitiesInfluencedByPhantomRemoval( set<VFace2D*>& phantomFaces, list<VEdge2D*>& virtualEdges, list<VEdge2D*>& unboundedEdges_before_phantom_removal, list<VEdge2D*>& edgesDefinedByPhantomNInputDisk );
        void        putRedTagsOnEntitiesToBeRemoved( list<VEdge2D*>& virtualEdges, list<VEdge2D*>& unboundedEdges_before_phantom_removal );
        void        adjustTopologyOfInterfaceBetweenInputDisksNPhantomDisks( VFace2D* virtualFace, list<VEdge2D*>& unboundedEdges_before_phantom_removal, set<VFace2D*>& phantomFaces, list<VEdge2D*>& edgesDefinedByPhantomNInputDisk );
        void        removePhantomFaces( set<VFace2D*>& phantomFaces );

        void        correctTopologyInUnboundedRegion();
        void        removeEdgesDefinedByOnlyPhantomGenerators();
        void        replaceRegionsOfPhantomGeneratorsWithVirtualRegion();
        void        removeVoronoiRegionsOfPhantomGenerators();
        void        makeMoreThanTwoEdgesBetweenRealAndVirtualRegionIntoOne();

        void        makeCorrectTopologyOfTheInterfaceByFlippingEdges( VFace2D* virtualFace );
        void        assignCoordinateOfInfiniteVertices();

        void        collectPossiblyFlippingEdgesOnTheInterface_StoredInPriorityQueue( VFace2D* virtualFace, Priority_Q_VEdge& possiblyFlippingEdges );
        void        removeRootNodeOfPriorityQueue( rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges, rg_dNode< pair<VEdge2D*, rg_Circle2D> >*& nodeWithSmallestCircle );
        //void        flipThisEdgeNFindTwoIncidentTriplets( VFace2D* virtualFace, pair<VEdge2D*, rg_Circle2D>& edgeNCircle, VEdge2D*& incidentTriplet1, VEdge2D*& incidentTriplet2, rg_dList<VEdge2D*>& edgesToUpdateGeometry );
        void        flipThisEdgeNFindTwoIncidentTriplets( VFace2D* virtualFace, pair<VEdge2D*, rg_Circle2D>& edgeNCircle, VEdge2D*& incidentTriplet1, VEdge2D*& incidentTriplet2 );
        void        reflectTheFlipOnTheIncidentTriplet( Priority_Q_VEdge& possiblyFlippingEdges, VEdge2D*& edgeOfIncidentTriplet );


        VFace2D*    getMatingFace( VVertex2D* tVertex, VEdge2D* tEdge );
        bool        isFlippingEdgeOnVirtualRegion( VEdge2D* edge, rg_Circle2D& circumcircle );
        double      angleFromVec1ToVec2( const rg_Point2D& vector1, const rg_Point2D& vector2 );

        void        collectUnboundedEdges( list<VEdge2D*>& unboundedEdges );
        bool        isCorrectCircumcircle( rg_Circle2D& circumcircle, VFace2D* prevFace, VFace2D* currFace, VFace2D* nextFace );
        pair<VEdge2D*, rg_Circle2D> dequeueEdgeWithMinTangentCircleFromPossiblyFlippingEdges( rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges );


        // topology handling
        //     Vertex

        inline VEdge2D*     getCCWNextEdgeOnVertex( VEdge2D* tEdge, VVertex2D* tVertex ) {
            VEdge2D* CCWNextEdge = NULL;
            if( tVertex == tEdge->getStartVertex() ) {
                CCWNextEdge = tEdge->getLeftLeg();
            }
            else if( tVertex == tEdge->getEndVertex() ) {
                CCWNextEdge = tEdge->getRightHand();
            }

            return CCWNextEdge;
        }



        inline VEdge2D*     getCWNexteEdgeOnVertex( VEdge2D* tEdge, VVertex2D* tVertex ) {
            VEdge2D* CWNextEdge = NULL;
            if( tVertex == tEdge->getStartVertex() ) {
                CWNextEdge = tEdge->getRightLeg();
            }
            else if( tVertex == tEdge->getEndVertex() ) {
                CWNextEdge = tEdge->getLeftHand();
            }

            return CWNextEdge;
        }



        inline VEdge2D*     getCCWNextEdgeOnFace( VEdge2D* tEdge, VFace2D* tFace ) {
            VEdge2D* CCWNextEdge = NULL;
            if( tFace == tEdge->getLeftFace() ) {
                CCWNextEdge = tEdge->getLeftHand();
            }
            else if( tFace == tEdge->getRightFace() ) {
                CCWNextEdge = tEdge->getRightLeg();
            }

            return CCWNextEdge;
        }



        inline VEdge2D*     getCWNextEdgeOnFace( VEdge2D* tEdge, VFace2D* tFace ) {
            VEdge2D* CWNextEdge = NULL;
            if( tFace == tEdge->getLeftFace() ) {
                CWNextEdge = tEdge->getLeftLeg();
            }
            else if( tFace == tEdge->getRightFace() ) {
                CWNextEdge = tEdge->getRightHand();
            }

            return CWNextEdge;
        }
            


        inline VFace2D*     getFaceOfOppositeSide( VFace2D* tFace, VEdge2D* tEdge )
        {
            VFace2D* oppositeSideFace = NULL;
            if( tFace == tEdge->getRightFace() )
            {
                oppositeSideFace = tEdge->getLeftFace();
            }
            else if( tFace == tEdge->getLeftFace() )
            {
                oppositeSideFace = tEdge->getRightFace();
            }

            return oppositeSideFace;
        }



        void    makeWEDS( const char* fname, list<rg_Circle2D>& circleList );

        
};

}
}

#endif


