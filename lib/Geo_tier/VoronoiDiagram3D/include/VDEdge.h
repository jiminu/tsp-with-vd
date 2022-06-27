#ifndef _VDEDGE_H
#define _VDEDGE_H

#include <fstream>
using namespace std;

#include "ConstForVoronoiDiagram3D.h"
#include "TopologicalEntity.h"
#include "rg_dList.h"

#include "rg_Curve.h"
#include "rg_Point3D.h"
#include "Sphere.h"
#include "Plane.h"

namespace V {

namespace GeometryTier {


class VDCell;
class VDFace;
class VDPartialEdge;
class VDVertex;

class BetaFace;

class VDEdge : public TopologicalEntity
{
private:
    VDVertex*      m_startVertex;
    VDVertex*      m_endVertex;

    VDPartialEdge* m_partEdge;


    rg_Curve*      m_edgeEquation;

    rg_INT         m_numOfPointsOnEdge;
    rg_Point3D*    m_pointsOnEdge;

	//donguk: minimum distance between a point on the curve and a generator.
	rg_REAL		   m_minRadius; 
	PositionOnEdge m_minPosition; //START/END/MID
	PositionOnEdge m_maxPosition; //START/END

    rg_FLAG        m_isOnInfinity;
    rg_FLAG        m_ProbeTangibility; 


    //  youngsong: geometric property of edge
	rg_FLAG        m_edgeType;
    Sphere         m_minMaxTangentSphereByBalls[2];
    Plane          m_edgePlane;
    rg_Point3D     m_pointForAngleDist;


    //  
    rg_BOOL        m_visited;
    BetaFace*      m_betaFace;


public:
    //  constructor & deconstructor..
    VDEdge();

    VDEdge( const rg_INT&        ID, 
                  VDVertex*      startVertex, 
                  VDVertex*      endVertex );
    
    VDEdge(       VDVertex*      startVertex, 
                  VDVertex*      endVertex, 
                  VDPartialEdge* partialEdge, 
                  rg_Curve*      edgeEquation, 
            const rg_INT&        numPtsOnEdge, 
                  rg_Point3D*    ptsOnEdge );
    
    VDEdge( const rg_INT&        ID,      
                  VDVertex*      startVertex, 
                  VDVertex*      endVertex, 
                  VDPartialEdge* partialEdge, 
                  rg_Curve*      edgeEquation, 
            const rg_INT&        numPtsOnEdge, 
                  rg_Point3D*    ptsOnEdge );
    
    VDEdge( const TopologicalEntity& aTopoEntity,    
                  VDVertex*          startVertex, 
                  VDVertex*          endVertex, 
                  VDPartialEdge*     partialEdge, 
                  rg_Curve*          edgeEquation, 
            const rg_INT&            numPtsOnEdge, 
                  rg_Point3D*        ptsOnEdge );

    VDEdge( const VDEdge& anEdge );
    
    ~VDEdge();

    //  get functions.. 
    inline BetaFace* getBetaFace() const { return m_betaFace; }
    inline rg_BOOL   isVisited()   const { return m_visited; }
    inline void      isVisited(const rg_BOOL& visited) { m_visited = visited; }

    rg_BOOL        isVirtual() const;

    VDVertex*      getStartVertex() const;
    VDVertex*      getEndVertex() const;
    VDVertex*      getOppositeVertex(const VDVertex* const vertex) const;

    VDCell*        getStartVCell() const;
    VDCell*        getEndVCell() const;


    VDPartialEdge* getPartialEdge() const;
    
    rg_Curve*      getEdgeEquation() const;
    
    rg_INT         getNumberOfPointsOnEdge() const;
    rg_Point3D*    getPointsOnEdge() const;
    rg_Point3D     getPointOnEdge(const rg_INT& i) const;
    rg_Point3D     getPointOnEdgeByEquation(const rg_REAL& t) const;
	rg_REAL		   getMinRadius() const;
	PositionOnEdge getMinPosition() const;
	PositionOnEdge getMaxPosition() const;
	rg_REAL		   getMaxRadius() const;
    
	rg_FLAG        getEdgeType() const;
    Sphere         getMinMaxTangentSphereByBalls(const rg_INT& i) const;
    Plane          getEdgePlane() const;
    rg_Point3D     getPointForAngleDistance() const;

    rg_FLAG        isOnInfinity() const;
    rg_FLAG        isBoundedEdge() const;
    rg_FLAG        isBounded() const;

    rg_FLAG        isTangible() const;
    void           isTangible(const rg_FLAG& tangibility);

	//donguk: whether or not this curve is monotone with respect to the distance to the spheres.
	//다시말해 curve 가 세 개의 중심으로 이루어지는 plane을 통과하느냐의 여부.
	rg_FLAG		   isMonotone() const;
	

    //rg_BOOL        isIncidentTo(VDFace* v_face) const;

    //  set functions..
    void setStartVertex( VDVertex* startVertex );
    void setEndVertex( VDVertex* endVertex );

    void setPartialEdge( VDPartialEdge* partialEdge );

	void setMinRadius(const rg_REAL& radius);
	void setMinPosition(const PositionOnEdge& position);
	void setMaxPosition(const PositionOnEdge& position);
    void setEdgeEquation( rg_Curve* edgeEquation );
    void setPointsOnEdge( const rg_INT&     numPtsOnEdge,
                                rg_Point3D* ptsOnEdge );

	void setInfinity();


    void setEdge(       VDVertex*      startVertex, 
                        VDVertex*      endVertex, 
                        VDPartialEdge* partialEdge, 
                        rg_Curve*      edgeEquation, 
                  const rg_INT&        numPtsOnEdge, 
                        rg_Point3D*    ptsOnEdge );

    void setEdge( const rg_INT&        ID,      
                        VDVertex*      startVertex, 
                        VDVertex*      endVertex, 
                        VDPartialEdge* partialEdge, 
                        rg_Curve*      edgeEquation, 
                  const rg_INT&        numPtsOnEdge, 
                        rg_Point3D*    ptsOnEdge ); 


	void setEdgeType(const rg_FLAG& edgeType);
    void setMinMaxTangentSphereByBalls(Sphere* minMaxSphere);
    void setMinMaxTangentSphereByBalls(const Sphere& sphere1, const Sphere& sphere2);
    void setEdgePlane(const Plane& edgePlane);
    void setPointForAngleDistance(const rg_Point3D& point);

    //  operator overloading..
    VDEdge& operator =( const VDEdge& anEdge );

    void fileOutVDEdge( ofstream& fout );

    //  topological operators..
    
    //  topological quires..
    //  output : ARRAY -> VDCell* onlyEdgeSharingCellsInCCW[3]; 
    void               inquireIntoOnlyEdgeSharingCellsInCCW( VDCell** onlyEdgeSharingCellsInCCW ) const;


    void               inquireIntoOnlyEdgeSharingCellsInCCW( rg_dList<VDCell*>& onlyEdgeSharingCellsInCCW ) const;

    rg_dList<VDCell*>* inquireIntoOnlyEdgeSharingCellsInCCW();


    void inquireIncidentEdgesAtStartVertex( rg_dList<VDEdge*>& incidentEdgesAtStartVertex );
    void inquireIncidentEdgesAtEndVertex( rg_dList<VDEdge*>& incidentEdgesAtEndVertex );

    void  inquireBoundingVertices(rg_dList<VDVertex*>& vertexList);
    void  inquireIncidentEdges( rg_dList<VDEdge*>& edgeList );
    void  inquireIncidentFaces(rg_dList<VDFace*>& faceList) const;
    void  inquireIncidentCells(rg_dList<VDCell*>& cellList) const;

    rg_BOOL isBoundingCell( VDCell* currCell ) const;
    rg_BOOL isIncidentTo(   VDFace* currFace ) const;

    VDCell* getMateCellOfStartVertex() const;
    VDCell* getMateCellOfEndVertex() const;


    //rg_INT searchIncidentCells(rg_dList<VDCell*>& cellList) const;
    rg_INT getNumCellsToDefineThisEdge() const;

    rg_BOOL findPositionOfBlockingSphericalProbeFrom(VDVertex* startingVtx, 
                                                     const rg_REAL& radiusOfBlockingSphericalProbe,
                                                     Sphere& position) const;
     

    void  connectBetaFace(   BetaFace* b_face);
    void  disconnectBetaFace(BetaFace* b_face);
};

} // namespace GeometryTier

} // namespace V


#endif

