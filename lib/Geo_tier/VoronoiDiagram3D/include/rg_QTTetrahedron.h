#ifndef _QTTETRAHEDRON_H
#define _QTTETRAHEDRON_H

#include "Sphere.h"
#include "TopologicalEntity.h"
#include "rg_QTEdge.h"
#include "rg_QTFace.h"

namespace V {

namespace GeometryTier {


//const rg_INT MULTIPLICITY_ANOMALY = 2;
const rg_INT DEGENERACY_ANOMALY   = 1;
const rg_INT MULTIPLICITY_ANOMALY_BY_TWO_COMMON_FACE   = 2;
const rg_INT MULTIPLICITY_ANOMALY_BY_THREE_COMMON_FACE = 3;
const rg_INT MULTIPLICITY_ANOMALY_BY_FOUR_COMMON_FACE  = 4;

class QTVertex;
class VDVertex;

class QTTetrahedron : public TopologicalEntity
{
private:
    QTVertex*      m_vertex[4];
    QTTetrahedron* m_neighbor[4];

    Sphere         m_emptyTangentSphere;
    rg_BOOL        m_visited;

    VDVertex*      m_vVertex;

public:

    QTTetrahedron();
    QTTetrahedron(const rg_INT& id);
    QTTetrahedron(const Sphere& emptyTSphere);
    QTTetrahedron(const rg_INT& id, const Sphere& emptyTSphere);
    QTTetrahedron(const QTTetrahedron& tetrahedron);
    ~QTTetrahedron();

    inline VDVertex* getVVertex() const { return m_vVertex; }

    QTVertex*       getVertex(const rg_INDEX& i) const;
    QTVertex**      getVertices();
    QTTetrahedron*  getNeighbor(const rg_INDEX& i) const;
    QTTetrahedron** getNeighbors();
    Sphere          getEmptyTangentSphere() const;

    QTTetrahedron*  getMateTetrahedron(const rg_INT& pos) const;
    QTTetrahedron*  getMateTetrahedron(QTVertex* vertex) const;

    rg_INT          getPosOfVertex(QTVertex* vertex) const;
    rg_INT          getMateVertexPosOfIncidentFace(QTTetrahedron* neighbor, const rg_INT& mateVertexPos);
    rg_INT          getMatePosInNeighbor(QTTetrahedron* neighbor, const rg_INT& neighborPos) const;
    rg_INT          getPosOfNeighbor(QTTetrahedron* neighbor) const;
    QTTetrahedron*  getNeighborInMultiplicityAnomaly() const;
    rg_INT          getEdgesNotSharedWithNeighbor(QTTetrahedron* neighbor, QTEdge* edges);

    inline rg_BOOL isVisited() const { return m_visited; }
    inline void    isVisited(const rg_BOOL& visited) { m_visited = visited; }
    rg_BOOL        isVirtual() const;

    rg_FLAG isThere(QTVertex* vertex) const;
    rg_FLAG isThere(QTVertex* vertex1, QTVertex* vertex2) const;

    rg_FLAG isThereThisNeighbor(QTTetrahedron* neighbor) const;

    rg_FLAG isInfiniteTetrahedron() const;
    rg_BOOL isIsolatedFace() const;

    rg_FLAG isAnomaly() const;
    rg_INT  isMultiplicity(QTTetrahedron*& neighbor) const;


    void setVertex(const rg_INDEX& i, QTVertex* vertex);
    void setVertices(QTVertex* vertex1, QTVertex* vertex2, 
                     QTVertex* vertex3, QTVertex* vertex4);
    void setVertices(QTVertex** vertices);
    void setNeighbor(const rg_INDEX& i, QTTetrahedron* neighbor);
    void setNeighbors(QTTetrahedron* neighbor1, QTTetrahedron* neighbor2, 
                      QTTetrahedron* neighbor3, QTTetrahedron* neighbor4);
    void setNeighbors(QTTetrahedron** neighbors);
    void setEmptyTangentSphere(const Sphere& emptyTSphere);

    QTTetrahedron& operator =(const QTTetrahedron& tetrahedron);


    /////  For Queres  /////
    void findCommonEdgesWithNeighbor(QTTetrahedron* neighbor, rg_dList<QTEdge>& edgeList);
    void findCommonFacesWithNeighbor(QTTetrahedron* neighbor, rg_dList<QTFace>& faceList);


    void inquireBoundingVerticesInThisWorld(rg_dList<QTVertex*>& vertexList);
    void inquireBoundingEdgesInThisWorld(rg_dList<QTEdge>& edgeList);
    void inquireBoundingFacesInThisWorld(rg_dList<QTFace>& faceList);
    void inquireIncidentTetrahedraInThisWorld(rg_dList<QTTetrahedron*>& tetrahedronList);


    // function to connect with vertex in VD
    void connectVVertex(   VDVertex* v_vtx);
    void disconnectVVertex(VDVertex* v_vtx);

};

} // namespace GeometryTier

} // namespace V


#endif
