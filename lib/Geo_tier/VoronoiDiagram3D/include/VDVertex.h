#ifndef _VDVERTEX_H
#define _VDVERTEX_H

#include "ConstForVoronoiDiagram3D.h"
#include "TopologicalEntity.h"

#include "rg_Point3D.h"
#include "rg_dList.h"

namespace V {

namespace GeometryTier {


class VDCell;
class VDFace;
class VDEdge;

class Gate;

class QTTetrahedron;
class BetaCell;

class VDVertex : public TopologicalEntity
{
private:
    VDEdge*    m_incidentEdges[ NUM_INCIDENT_EDGES_OF_ONE_VERTEX ];
    Gate**     m_gates;

    rg_FLAG    m_isOnInfinity;

    rg_Point3D m_point;
    rg_REAL    m_radiusOfTangentSphere;
    rg_FLAG    m_ProbeTangibility; 

    rg_BOOL        m_visited;
    QTTetrahedron* m_qCell;
    BetaCell*      m_betaCell;


public:
    //  constructor & deconstructor..
    VDVertex();
    VDVertex( const rg_INT&     ID);
    VDVertex( const rg_FLAG&    onInfinity, 
              const rg_Point3D& point, 
              const rg_REAL&    radiusOfTS );

    VDVertex(       Gate**      gates, 
              const rg_FLAG&    onInfinity, 
              const rg_Point3D& point, 
              const rg_REAL&    radiusOfTS );

    VDVertex( const rg_INT&     ID, 
              const rg_FLAG&    onInfinity, 
              const rg_Point3D& point, 
              const rg_REAL&    radiusOfTS );

    VDVertex( const rg_INT&     ID, 
                    Gate**      gates, 
              const rg_FLAG&    onInfinity, 
              const rg_Point3D& point, 
              const rg_REAL&    radiusOfTS );

    VDVertex( const TopologicalEntity& aTopoEntity, 
              const rg_FLAG&           onInfinity, 
              const rg_Point3D&        point, 
              const rg_REAL&           radiusOfTS );

    VDVertex( const TopologicalEntity& aTopoEntity, 
                    Gate**             gates, 
              const rg_FLAG&           onInfinity, 
              const rg_Point3D&        point, 
              const rg_REAL&           radiusOfTS );

    VDVertex( const VDVertex& aVertex );
    
    ~VDVertex();


    //  get functions.. 
    inline QTTetrahedron* getQCell() const { return m_qCell; }
    inline BetaCell*      getBetaCell() const { return m_betaCell; }
    inline rg_BOOL        isVisited()   const { return m_visited; }
    inline void           isVisited(const rg_BOOL& visited) { m_visited = visited; }


    VDEdge*    getIncidentEdge( const rg_INT& i ) const;
    VDEdge**   getAllIncidentEdges();

    Gate*      getGate( const rg_INT& i ) const;
    Gate**     getAllGates() const;

    rg_FLAG    isOnInfinity() const;

    rg_Point3D getPoint() const;
    rg_REAL    getRadiusOfTangentSphere() const;

    rg_FLAG    isConnectedWithUnboundedEdge() const;

    rg_FLAG    isTangible() const;
    void       isTangible(const rg_FLAG& tangibility);

    //  set functions..
    void    setIncidentEdge( const rg_INT& i, VDEdge* anEdge );
    void    setAllIncidentEdge( VDEdge** incidentEdges );
    rg_FLAG setIncidentEdge( VDEdge* anEdge );

    void    setGate( const rg_INT& i, Gate* aGate );
    void    setAllGates( Gate** gates );
    rg_FLAG setGate( Gate* aGate );

    void isOnInfinity( const rg_FLAG& onInfinity );

    void setPoint( const rg_Point3D& point );
    void setRadiusOfTangentSphere( const rg_REAL& radiusOfTS );

    void removeAllGates();

    //  operator overloading..
    VDVertex& operator =(const VDVertex& aVertex );
    

    //  topological operators..
    void inquireAllCellsToDefineThisVertex(VDCell** cellToDefineVertex) const;

    //  topological quires..
    void  inquireNeighborVertices(rg_dList<VDVertex*>& vertexList) const;
    void  inquireIncidentEdges(rg_dList<VDEdge*>& edgeList) const;
    void  inquireIncidentFaces(rg_dList<VDFace*>& faceList) const;
    void  inquireIncidentCells(rg_dList<VDCell*>& cellList) const;

    rg_BOOL isBoundingCell( VDCell* currCell ) const;
    rg_BOOL isIncidentTo(   VDEdge* currEdge ) const;

    rg_BOOL hasTwinVertex() const;
    VDVertex* getTwinVertex() const;
    void    getMateCellForIncidentEdges(rg_dList<VDCell*>& cellList) const;

    //rg_INT searchIncidentCells(rg_dList<VDCell*>& cellList) const;
    rg_INT getNumCellsToDefineThisVertex() const;
    rg_INT getNumIncidentEdges() const;

    // function to connect with cell in IWDS or eIWDS
    void connectQCell(   QTTetrahedron* q_cell);
    void disconnectQCell(QTTetrahedron* q_cell);

    void connectBetaCell(   BetaCell* b_cell);
    void disconnectBetaCell(BetaCell* b_cell);
};

} // namespace GeometryTier

} // namespace V


#endif

