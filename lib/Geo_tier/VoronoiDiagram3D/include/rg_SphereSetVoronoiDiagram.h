#ifndef _SPHERESETVORONOIDIAGRAM_H
#define _SPHERESETVORONOIDIAGRAM_H

#include "ConstForVoronoiDiagram3D.h"

#include "rg_dList.h"

#include "VoronoiDiagram3D.h"
#include "rg_BallGenerator.h"
#include "VIDIC.h"
//#include "Pair.h"
#include "rg_Circle2D.h"
#include "rg_ImplicitEquation.h"
#include "rg_NURBSplineSurface3D.h"
#include "rg_RevolvedSurface.h"

#include "rg_Triple.h"
#include "NavigatingBucketCell.h"

#include "rg_QuasiTriangulation.h"

#include <time.h>

#include <utility>
#include <fstream>
#include <string>
using namespace std;


const rg_INT NUM_TIMES_FOR_ANALYSIS = 17;

namespace V {

namespace GeometryTier {


class rg_SphereSetVoronoiDiagram : public VoronoiDiagram3D
{
private:
    rg_dList< BallGenerator > m_GlobalGeneratorList;

    rg_Point3D m_minPointOfBoundingBoxForBalls;
    rg_Point3D m_maxPointOfBoundingBoxForBalls;

    rg_Point3D m_minPointOfBoundingBoxForVD;
    rg_Point3D m_maxPointOfBoundingBoxForVD;

    rg_REAL    m_AvgRadiusOfBalls;

    rg_FLAG    m_hasProbeTangibility;
    rg_REAL    m_probeRadius;
    rg_REAL    m_sizeOfBucketElement;

    //  for debugging.
    rg_INT ON_DEBUG;
    ofstream foutEVDS;
//    rg_Point3D** excludingPlangOfInitialGates;


public:

    rg_REAL    m_minRadiusOfBalls;
    rg_REAL    m_maxRadiusOfBalls;

    VDVertex*    startVertexOnSolventTrajectory;
    VDVertex*    endVertexOnSolventTrajectory;


    clock_t    start, finish;   
    clock_t    start_in, start_in2, finish_in;
    clock_t    computingTime[NUM_TIMES_FOR_ANALYSIS];
                //  computingTime[0]:  initializing EVD(S)
                //  computingTime[1]:  finding end vertex
                //  computingTime[2]:    preparing end vertex finding
                //  computingTime[3]:    finding end vertex of parabolic or hyperbolic edge.
                //  computingTime[4]:      finding end vertex of parabolic or hyperbolic edge in front feasible region
                //  computingTime[5]:      finding end vertex of parabolic or hyperbolic edge from front feasible region
                //  computingTime[6]:      finding end vertex of parabolic or hyperbolic edge in rear feasible region
                //  computingTime[7]:    finding end vertex of linear edge.
                //  computingTime[8]:      finding end vertex of linear edge in front feasible region
                //  computingTime[9]:      finding end vertex of linear edge from front feasible region
                //  computingTime[10]:     finding end vertex of linear edge in rear feasible region
                //  computingTime[11]:    finding end vertex of circular or elliptic edge
                //  computingTime[12]: exploring edge
                //  computingTime[13]: setting local topology
                //  computingTime[14]: setting global topology
                //  computingTime[15]: construct small worlds.
                //  computingTime[16]: construct VD(B).
    clock_t    timeToBeExtracted[2];
                //  timeToBeExtracted: time to count candidate balls in non-elliptic edge.

    rg_FLAG  m_bStartCheck;
    rg_FLAG  m_isFrontFeasibleRegion;

    /*
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //
    rg_INT   m_numCandidateBallsOfNonElliptic[2][2][3];
                //  m_numCandidateBallsOfNonElliptic[0]:  balls intersecting feasible region
                //  m_numCandidateBallsOfNonElliptic[1]:  balls actually checked 
                //  m_numCandidateBallsOfNonElliptic[i][0]:  for parabolic or hyperbolic edge
                //  m_numCandidateBallsOfNonElliptic[i][1]:  for linear edge
                //  m_numCandidateBallsOfNonElliptic[i][j][0]:  balls in only front feasible region
                //  m_numCandidateBallsOfNonElliptic[i][j][1]:  balls from front feasible region
                //  m_numCandidateBallsOfNonElliptic[i][j][2]:  balls in only rear feasible region
    rg_dList<BallGenerator*> m_feasibleBalls;

    rg_REAL  m_avgRadiusOfSphereForRf[2][3];
                //  m_avgRadiusOfSphereForRf[0]:  for parabolic or hyperbolic edge
                //  m_avgRadiusOfSphereForRf[1]:  for linear edge
                //  m_avgRadiusOfSphereForRf[j][0]:  balls in only front feasible region
                //  m_avgRadiusOfSphereForRf[j][1]:  balls from front feasible region
                //  m_avgRadiusOfSphereForRf[j][2]:  balls in only rear feasible region
    rg_INT   m_numEdges[2][3];
    rg_INT   m_numEndBallsInPosition[2][3];
                //  m_numEdges[0]:  parabolic or hyperbolic edge
                //  m_numEdges[1]:  linear edge
                //  m_numEdges[j][0]:  with end ball in only front feasible region
                //  m_numEdges[j][1]:  with end ball from front feasible region
                //  m_numEdges[j][2]:  with end ball in only rear feasible region
    rg_INT   m_numIntersectingBallsWithFrontFR[2];
    rg_INT*  m_frequencyOfCandidateBall;


    rg_INT   m_numCandidateBallsOfElliptic[2];
                //  m_numCandidateBallsOfElliptic[0]:  balls intersecting feasible region
                //  m_numCandidateBallsOfElliptic[1]:  balls actually checked 

    rg_INT   m_numBallsForUnboundedEdge[2];
    rg_INT   m_numSearchedBucketElements[2][4];
    //
    //////////  FOR ANALYSIS OF EDGE-TRACING  //////////
    //*/    
    
    //  constructor & deconstructor..
    rg_SphereSetVoronoiDiagram();
    rg_SphereSetVoronoiDiagram(const rg_SphereSetVoronoiDiagram& VD);
    ~rg_SphereSetVoronoiDiagram();

    //  get functions.. 
    rg_dList< BallGenerator >* getGlobalGeneratorList();
    inline const rg_dList< BallGenerator >& getGlobalGeneratorList() const {return m_GlobalGeneratorList;}

    rg_INT     getNumOfBalls() const;
    rg_Point3D getMinPointOfBoundingBoxForBalls() const;
    rg_Point3D getMaxPointOfBoundingBoxForBalls() const;
    rg_Point3D getMinPointOfBoundingBoxForVoronoiDiagram() const;
    rg_Point3D getMaxPointOfBoundingBoxForVoronoiDiagram() const;
    rg_REAL    getAverageRadiusOfBalls() const;

    rg_FLAG    isProbeTangibilityCalculated() const;
    rg_REAL    getProbeRadius() const;

    inline     rg_REAL getSizeOfBucketElement() const { return m_sizeOfBucketElement; }


    //  set functions..
    void setAvgRadiusOfBalls(const rg_REAL& avgRadius);
    void setMinPointOfBoundingBoxForBalls(const rg_Point3D& minPt);
    void setMaxPointOfBoundingBoxForBalls(const rg_Point3D& maxPt);
    void setBucket(const rg_REAL& loadFactor);

    BallGenerator*  addBallGenerator(const BallGenerator& aBall);
    rg_BOOL         addSphere(const Sphere& sphere, void* property = rg_NULL, const rg_INT& IDFromInput = -1);
    rg_BOOL         addSphereWithoutItsValidityCheck(const Sphere& sphere, void* property = rg_NULL, const rg_INT& IDFromInput = -1);



    rg_INT getNumOfVoronoiCells();
    rg_INT getNumOfUnboundedVoronoiCells();
    rg_INT getNumOfVoronoiFaces();
    rg_INT getNumOfVoronoiEdges();
    rg_INT getNumOfVoronoiVertices();


    void duplicate(const rg_SphereSetVoronoiDiagram& origin);
    void duplicateSubsetByComplement(const rg_SphereSetVoronoiDiagram& origin, const rg_dList<VDEdge*>& edgesInSubset);
    void clean();


    void makeEntityMapsForSubset( //const rg_REAL& betaValue, 
                                  const rg_dList<VDEdge*>&             edgesOnComponent, 
                                  map<VDCell*, VDCell*>&               cellMap,
                                  map<VDFace*, VDFace*>&               faceMap,
                                  map<VDLoop*, VDLoop*>&               loopMap,
                                  map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                  map<VDEdge*, VDEdge*>&               edgeMap,
                                  map<VDVertex*, VDVertex*>&           vertexMap);
    void setTopologyOfVerticesInSubset( const map<VDVertex*, VDVertex*>&           vertexMap,
                                        const map<VDEdge*, VDEdge*>&               edgeMap );
    void setTopologyOfEdgesInSubset(    const map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                        const map<VDEdge*, VDEdge*>&               edgeMap,
                                        const map<VDVertex*, VDVertex*>&           vertexMap);
    void setTopologyOfPartialEdgesInSubset( const map<VDLoop*, VDLoop*>&           loopMap,
                                            const map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                            const map<VDEdge*, VDEdge*>&               edgeMap);
    void setTopologyOfLoopsInSubset(    const map<VDFace*, VDFace*>&               faceMap,
                                        const map<VDLoop*, VDLoop*>&               loopMap,
                                        const map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap);
    void setTopologyOfFacesInSubset(    const map<VDCell*, VDCell*>&               cellMap,
                                        const map<VDFace*, VDFace*>&               faceMap,
                                        const map<VDLoop*, VDLoop*>&               loopMap);
    void setTopologyOfCellsInSubset(    const map<VDCell*, VDCell*>&               cellMap,
                                        const map<VDFace*, VDFace*>&               faceMap);

private:
    void duplicateGenerators(const rg_SphereSetVoronoiDiagram& origin);
    void perturbGenerators(map<BallGenerator*, rg_Point3D>& geneCenterMap);


public:

    //  operator overloading..


    ///////////////////////////////////////////////////////////////////////////
    //
    //  constructing functions of sphere set voronoi diagram
    rg_INT constructVoronoiDiagramOfSpheresByEdgeTracing();





    ///////////////////////////////////////////////
    //  by Youngsong Cho (2006. 7. 20)
    rg_INT constructSphereVoronoiDiagramRobustly();
        void constructBigWorldRobustly();
        rg_FLAG locateInitialVoronoiVertexInBigWorldRobustly(VIDIC& theVIDIC, rg_dList<Gate>& theGateStack ); 
            void makeBallConfigurationForVirtualVertex(VIDIC& theVIDIC, rg_dList<Gate>& theGateStack ); 
        rg_FLAG constructThisWorldRobustly(VDLoop* thisWorld, VIDIC& theVIDIC, rg_dList<Gate>& theGateStack); 
            void findEndVertexRobustly(Gate& currGate);
                void findEndVertexForParabolicEdgeRobustly(Gate& currGate);
                void findEndVertexForLinearEdgeRobustly(Gate&    currGate);
                void findEndVertexForEllipticEdgeRobustly(Gate&  currGate);
            void defineUnboundedEdge(VDLoop* thisWorld, VIDIC& theVIDIC, rg_dList<Gate>& theGateStack, const Gate& currGate);
            void defineBoundedEdge(VDLoop* thisWorld, VIDIC& theVIDIC, rg_dList<Gate>& theGateStack, const Gate& currGate);
            void solveDegeneracyAtEndVertexOfCurrentEdge(VDLoop* thisWorld, VIDIC& theVIDIC, rg_dList<Gate>& theGateStack, Gate& currGate);
                void updateLocalTopologyByEdge(VDLoop* thisWorld, const rg_FLAG& edgeType, 
                                               VDVertex* startVertex, VDVertex* endVertex, BallGenerator** ballToDefineEdge ); 
                void updateLocalTopologyByEdge(VDLoop* thisWorld, const Gate& currGate, VDVertex* endVertex); 
        void constructSmallWorldRobustly();
        void updateGlobalTopology();
            void defineOrderOfPartialEdgesInAllLoops();
                void makeAdditiveOuterLoopOfFace(VDFace* ptrFace, rg_dList<VDPartialEdge*>& multipleOuterLoop);
                void makeVDFaceOfAdditiveOuterLoop(VDFace* ptrFace, rg_dList<VDPartialEdge*>& multipleOuterLoop);
            void setTopologyOnInfinity();


    //
    ///////////////////////////////////////////////


    rg_INT constructSphereSetVoronoiDiagramByEdgeTracing();
        void    findInitialVertexForEdgeTracing( VIDIC& theVIDIC, rg_dList<Gate>& theGateStack );        
        rg_FLAG findEndVertexBasedOnAngularDistance( Gate& gate );


    rg_INT constructSphereVoronoiDiagramByAcceleratedEdgeTracing();
        //  0. find initial voronoi vertex.
        //       0.0  construct initial voronoi vertex by the selected 4 balls 
        //       0.1  do intersection test of the found vertex.
        //       0.2  create and push 4 gates into the GATE-STACK.
        //       0.3  insert the found vertex into the VIDIC.
    
        void shrinkBallsByMinimumRaidus();
        void recoverBallsByMinimumRaidus();
        
        rg_FLAG constructBigWorld(); 
            rg_FLAG   findInitialVertexInBigWorld( VIDIC& theVIDIC, rg_dList<Gate>& theGateStack );

            void defineInitialVoronoiVertex(       BallGenerator** surroundingBalls,
                                             const Sphere&         tangentSphere,
                                                   VIDIC&          theVIDIC,
                                                   rg_dList<Gate>& theGateStack );     

        rg_FLAG findEndVertexBasedOnAngularDistanceInAcceleration( Gate& gate );

            rg_FLAG findEndVertexOfParabolicAndHyperbolicEdge(       Gate&       gate, 
                                                               const rg_Point3D& normalOfCenterPlane, 
                                                               const rg_Point3D& normalOfEdgePlane,
                                                               const rg_Point3D& axisPoint,
                                                               const Sphere&     minSphere);
                void     findEmptyTangentSphereWithMinAngleInParabolicAndHyperbolicEdge(
                                                                  pair<BallGenerator*, Sphere>& endBall,
                                                                  rg_REAL&                      minAngularDistance,
                                                                  rg_dList<BallGenerator*>&     candidateBalls,
                                                                  Sphere*                       gateBall,
                                                                  BallGenerator*                startBall,   
                                                                  const rg_Point3D&             coordStartVertex,
					                                              const rg_Point3D&             normalOfEdgePlane,
					                                              const rg_Vector3D&            vecPaVs,
					                                              const rg_Point3D&             axisPoint,
                                                                  const rg_REAL&                maxValidAngle);
                void initializePropagationInRearFeasibleRegionOfNonEllipticEdge(
                                                               rg_dList< NavigatingBucketCell >& queueOfCells,
                                                               rg_Point3D&                   unitLength,
                                                               Sphere*                       gateBall,
                                                               rg_REAL*                      gatePlane );
                rg_Point3D computeCenterOfApproximatingSphereForR1( rg_Point3D* points,
                                                                    const rg_Point3D& normalOfCenterPlane,
                                                                    const rg_REAL&    coeffCenterPlane );
            rg_FLAG findEndVertexOfLinearEdgeInAcceleration(       Gate&       gate, 
                                                             const rg_Point3D& normalOfCenterPlane, 
                                                             const Sphere&     minSphere);
                void     findEmptyTangentSphereWithMinAngleInLinearEdge(
                                                                  pair<BallGenerator*, Sphere>& endBall,
                                                                  rg_REAL&                      minAngularDistance,
                                                                  rg_dList<BallGenerator*>&     candidateBalls,
                                                                  Sphere*                       gateBall,
                                                                  BallGenerator*                startBall,   
                                                                  const rg_Point3D&             coordStartVertex,
					                                              const rg_Point3D&             normalOfCenterPlane);
            rg_FLAG findEndVertexOfCircularAndEllipticEdgeInAcceleration(       Gate&       gate, 
                                                                          const rg_Point3D& normalOfCenterPlane, 
                                                                          const rg_Point3D& normalOfEdgePlane,
                                                                          const rg_Point3D& axisPoint,
                                                                          Sphere*           tangentSphereInCP);
                void     findEmptyTangentSphereWithMinAngleInEllipticAndCircularEdge(
                                                                  pair<BallGenerator*, Sphere>& endBall,
                                                                  rg_REAL&                      minAngularDistance,
                                                                  rg_dList<BallGenerator*>&     candidateBalls,
                                                                  Sphere*                       gateBall,
                                                                  BallGenerator*                startBall,   
                                                                  const rg_Point3D&             coordStartVertex,
					                                              const rg_Point3D&             normalOfEdgePlane,
					                                              const rg_Vector3D&            vecPaVs,
					                                              const rg_Point3D&             axisPoint);

            void    propagateBucketCells(const BucketCellIndex& cell, 
                                         const rg_REAL&         distance,
                                         const rg_Point3D&      unitLength, 
                                         rg_dList< rg_Triple< rg_INT, BucketCellIndex, rg_REAL > >* cellList);
            void    propagateBucketCellsOnPlane( const NavigatingBucketCell&       cell, 
                                                 const rg_Point3D&                 unitLength, 
                                                 rg_dList< NavigatingBucketCell >& queueOfCells);
            void    propagateBucketCellsInPlusPlane( const NavigatingBucketCell&       cell, 
                                                     const rg_INT&                     normalSign,
                                                     rg_dList< NavigatingBucketCell >& queueOfCells);

            void    computeTangentPlaneOfGate(const Gate& gate, rg_REAL* tangentPlane);
            rg_FLAG computeEdgePlaneForAngleDistanceAccordingToEdgeType(
                                           const Gate&       gate, 
                                                 rg_Point3D& normalOfCenterPlane,
                                                 rg_Point3D& normalOfEdgePlane,
                                                 rg_Point3D& axisPoint,
                                                 Sphere*     tangentSphereOnCP);
                rg_REAL angleCCW(const rg_Point2D& vector1, const rg_Point2D& vector2);
                rg_REAL angleCCW(const rg_Point3D& normal, const rg_Point3D& vector1, const rg_Point3D& vector2);
                rg_REAL angleCCWOfTwoUnitVectors(const rg_Point3D& normal, const rg_Point3D& vector1, const rg_Point3D& vector2);          
        VDVertex* exploreEdgeToBeMadeByCurrentGate( const Gate           &gate,
                                                          VIDIC          &theVIDIC,
                                                          rg_dList<Gate> &theGateStack);



        void makeVirtualTwoBigBrothers();

        void constructSmallWorlds();
            void      collectIsolatedBalls(rg_dList<BallGenerator*>& isolatedBalls);
            VDFace*   findBigBrothers( BallGenerator* smallBrother, BallGenerator** bigBrother );
            //VDVertex* findInitialVertexInSmallWorld(BallGenerator*  smallBrother, BallGenerator** bigBrother,
            //                                        VIDIC&          theVIDIC,     rg_dList<Gate>& theGateStack );
            rg_FLAG   findInitialVertexInSmallWorld(BallGenerator*  smallBrother, BallGenerator** bigBrother,
                                                    VDVertex*&      ptrVertex,    VIDIC&          theVIDIC, 
                                                    rg_dList<Gate>& theGateStack );
            void      constructSmallWorldWithoutVertices(BallGenerator*  smallBrother, 
                                                         BallGenerator** bigBrother,
                                                         VDFace*         faceByBigBrothers);
                void setTopologyByTwoBigBallsInSmallWorldWithoutVertices(
                                                              BallGenerator** bigBrother,
                                                              VDFace*         faceByBigBrothers,
                                                              VDFace**        facesBySmallAndTwoBigBrothers);

            void      setTopologyBySingleSmallBall(BallGenerator*  smallBrother, 
                                                   BallGenerator** bigBrother,
                                                   VDFace*         faceByBigBrothers);
            void      traceEdgesInSmallWorld(VDFace*         faceByBigBrothers,
                                             VIDIC&          theVIDIC,     
                                             rg_dList<Gate>& theGateStack );
                void  setLocalTopologyByEdgeInSmallWorld( VDFace*         faceByBigBrothers,
                                                          const Gate&     gate, 
                                                          VDVertex*       endVertex,
                                                          rg_dList<VDLoop*>& loopsInSmallWorld);
                void  arrangePartialEdgesInLoops( rg_dList<VDLoop*>& loopsInSmallWorld );



    VDEdge*   setLocalTopologyByEndVertex( const Gate &gate, VDVertex* endVertex ); 
    void setTopology();

	void setGeometry();
        void calculateBoundingBoxForVoronoiVertices();
		void setPropertyOfVoronoiEdges();
		void setGeometryOfVoronoiEdges();
        void setEndVertexOfInfiniteEdge();
            void evaluateEndVertexOfInfiniteEdge(const rg_Point3D& minPoint,
                                                 const rg_Point3D& maxPoint, 
                                                 Sphere*     balls, 
                                                 rg_dList<rg_Point3D>& candidatePointList);
            rg_INT evaluateInfiniteVertexOnOneBoundingFace(      rg_REAL* roots,
                                                           const rg_REAL& a_ij, 
                                                           const rg_REAL& b_ir, const rg_REAL& b_jr, 
                                                           const rg_REAL& c_ie, const rg_REAL& c_je,
                                                           const rg_REAL& boundingValue);
        void defineEndVertexOfInfiniteEdge(VDEdge*           infiniteEdge, 
                                           Sphere*           ballsInCCW,
                                           const rg_Point3D& minPtOfBox,
                                           const rg_Point3D& maxPtOfBox);
        void defindGeometryOfVoronoiEdge(VDEdge* currEdge, Sphere* ballsInCCW);

    void setAttribute();
        void setBoundednessOfFacesAndCells();
    

    VDVertex*       makeVoronoiVertex( const rg_FLAG& onInfinity, const Sphere& tangentSphere);
    VDEdge*         makeVoronoiEdge( VDVertex* startVertex, VDVertex* endVertex ); 
    VDPartialEdge** make3VoronoiPartialEdges( VDEdge* theEdge );        
    void            makeVoronoiFaceAndItsOuterLoop( const Gate &gate, VDPartialEdge** partialEdges ); 
    void            makeVoronoiFaceAndItsOuterLoop( BallGenerator** ballToDefineEdge , VDPartialEdge** partialEdges ); 







            
    ///////////////////////////////////////////////////////////////////////////
    //
    //  Computation of Geometric Entities for Voronoi Diagram
    //
    //    1. voronoi vertex computation 
    rg_INT evaluate_Voronoi_Vertex(Sphere*  s1, Sphere* s2, Sphere* s3, Sphere* s4,
                                   Sphere*& tangentSpheres);
    rg_INT evaluate_Voronoi_Vertex_SameRadius(Sphere* s1, Sphere* s2, Sphere* s3, Sphere* s4, 
                                              Sphere*& tangentSpheres);
        rg_REAL getDeterminant( const rg_REAL& mat11, const rg_REAL& mat12, const rg_REAL& mat13,
        					    const rg_REAL& mat21, const rg_REAL& mat22, const rg_REAL& mat23,
		        			    const rg_REAL& mat31, const rg_REAL& mat32, const rg_REAL& mat33 );
        


    ///////////////////////////////////////////////////////////////////////////
    //
    //    2. voronoi edge computation 	by Donguk Kim
    rg_FLAG evaluateVoronoiEdge(       rg_RBzCurve3D*& curve,
                                 const rg_Point3D&     startPt,
                                 const rg_Point3D&     endPt,
                                 const Sphere&         s1,
                                 const Sphere&         s2,
                                 const Sphere&         s3 );
	    rg_Point3D getTangentVectorOfEdgeAtVertex(const rg_Point3D& vertexPt,
											      const Sphere& s1,
											      const Sphere& s2,
											      const Sphere& s3 );

	    rg_INT makeCircumcircle(const rg_Circle2D& circle1,
							    const rg_Circle2D& circle2,
							    const rg_Circle2D& circle3,
							    rg_Circle2D&	result1, 
							    rg_Circle2D&	result2);
	    void shrinkCircle(const		rg_Circle2D& c1,
					      const		rg_Circle2D& c2,
					      const		rg_Circle2D& c3,
					      rg_Circle2D&	c1Tilde,
					      rg_Circle2D&	c2Tilde,
					      rg_Point2D&		c3Tilde);
	    rg_Point2D makeTangentLinesInWPlane(const rg_Circle2D&		c1, 
										    const rg_Circle2D&		c2,
										    const rg_Circle2D&		c3,
										    rg_ImplicitEquation&		result1,
										    rg_ImplicitEquation&		result2);

	    void makeTangentLinesOf2Circles(const rg_Circle2D&	circle1,
									    const rg_Circle2D&	circle2,
									    rg_ImplicitEquation&	result1,
									    rg_ImplicitEquation&	result2);

	    rg_Circle2D transformW2Z(const rg_ImplicitEquation& line,
							     const rg_Point2D& z3);

	    rg_Circle2D transformZ2W(const rg_Circle2D&	cTilde,
							     const rg_Point2D&	c3Tilde);


        rg_INT makeCircumcircleOnCenterPlane(const Sphere& s1,
											 const Sphere& s2,
											 const Sphere& s3,
											 Sphere& result1,
											 Sphere& result2);

    //  topological operators..

    //  geometric operators..

    Sphere evaluateSmallestTangentSphereFor3Spheres(const rg_Point3D& normalOfPlaneBy3SphereCenters,
                                                    const Sphere&     sphereInCCW1,
                                                    const Sphere&     sphereInCCW2,
                                                    const Sphere&     sphereInCCW3);

    void fineInternalVoids(const rg_Point3D& minPtOfBoundingBox,
                           const rg_Point3D& maxPtOfBoundingBox,
                                 rg_INT&     numOfInternalVoids,
                                 VDVertex**& internalVoids);
    rg_dList<VDFace*>* findInteractionInterface();
    //  rg_REAL distancesBetweenCenters[4] 
    //  rg_REAL distancesBetweenBalls[4]  
    //      distances..[0] : sum of distance
    //      distances..[1] : average of distance
    //      distances..[2] : maxinum distance
    //      distances..[3] : mininum distance
    void analyzeDistanceBetweenGroups( rg_dList<VDFace*>* facesWithTwoGeneratorInDifferentGroup,
                                       rg_REAL*           distancesBetweenCenters,
                                       rg_REAL*           distancesBetweenBalls );

	rg_dList<rg_NURBSplineSurface3D>* computeLinkBlends(const rg_REAL& radius);
	void computeLinkBlendToEdge(VDEdge* theEdge, 
								const rg_REAL& radius, 
								const PositionOnEdge& position, 
								rg_RevolvedSurface& computedSurf);

	rg_dList<rg_NURBSplineSurface3D>* computeRollingBlends(const rg_REAL& radius);
	void computeRollingBlendToFace(VDFace* theFace, 
								   const rg_REAL& radius, 
								   rg_dList<rg_NURBSplineSurface3D>* rollingBlends);
	rg_Point3D computeTangentProbe(VDEdge* theEdge,
								   const rg_REAL& radius,
								   const PositionOnEdge& position);

	//curve를 xy-plane에 projection시킨 후  curve를 생성한다. (수치오차때문에)
	void makeConicByProjection(const rg_Point3D& sPt, const rg_Point3D& tangentStart,
							   const rg_Point3D& ePt, const rg_Point3D& tangentEnd,
							   const rg_Point3D& passPt,
							   rg_RBzCurve3D& curve);

    void defineProbeTangibility(const rg_REAL& radiusOfProbe);
        //  link blend 존재유무 검사.
        void defineProbeTangibleEdges(const rg_REAL& radiusOfProbe);
        //  rolling blend 존재유무 검사.
        void defineProbeTangibleFaces(const rg_REAL& radiusOfProbe);
//    void findProbeTangibleBalls(const rg_REAL& radiusOfProbe);
//        void findProbeTangibleEdges(const rg_REAL& radiusOfProbe);
//        void findFacesWithProbeTangibleBalls(const rg_REAL& radiusOfProbe);
        void resetTangibility();

    rg_dList< pair<VDEdge*, rg_FLAG> >* findSolventTrajectory();

    

    //  check Voronoi diagram
    rg_FLAG isEmptySphere(const Sphere& aSphere, const rg_REAL& res=rg_MATH_RES);
    rg_FLAG isValid();
    
    void    perturbBallGenerators();

    //  file I/O
    void   readGenerators(const char* filename);
    rg_INT readGeneratorsInPDB(const char* filename);


    //  Quasi-Triangulation
    void checkLinkageVD2QT();

    void convertToQuasiTriangulation(QuasiTriangulation& quasiTriangulation);
        QTVertex**      initializeVerticesInQT(QuasiTriangulation& quasiTriangulation);
        QTTetrahedron** initializeTetrahedraInQT(QuasiTriangulation& quasiTriangulation);
        void            setTopologyBetweenVerticesAndTetrahedraInQT(QTVertex**      verticesInQT, 
                                                                    QTTetrahedron** tetrahedraInQT);
        void            setGatesInQT(QTVertex**          verticesInQT, 
                                     QTTetrahedron**     tetrahedraInQT,
                                     QuasiTriangulation& quasiTriangulation);
        void            linkVDandQT( QTVertex**      verticesInQT, QTTetrahedron** tetrahedraInQT);
            void            linkVCellAndQTVertex( QTVertex**      verticesInQT);
            void            linkVVertexAndQTCell( QTTetrahedron** tetrahedraInQT);

    void fileOutDelaunayEdgeLength( ofstream& fout );
    rg_FLAG reportConstructionEVDS(ofstream& fout);
    void reportVDSforQT(ofstream& fout);

    void reportNumericalErrorInTangentSphereComputationOf4Balls(ofstream& fout);


    void analyzeSphereVoronoiDiagram(ofstream& fout);
    void analyzeNumNeighborAtoms(ofstream& fout);
    void analyzeFrequencyOfNumEdgesPerFace(ofstream& fout);
    void analyzeDistanceNeighborBalls(ofstream& fout);
    void analyzeDistanceBtwBoundedNeighborBalls(ofstream& fout);
    

};

} // namespace GeometryTier

} // namespace V



#endif
