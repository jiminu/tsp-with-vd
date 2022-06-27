#ifndef _QTFACE_H
#define _QTFACE_H

#include "TopologicalEntity.h"
#include "rg_dList.h"

namespace V {

namespace GeometryTier {


class QTVertex;
class QTEdge;
class QTTetrahedron;

class QTFace : public TopologicalEntity
{
private:
    QTTetrahedron* m_tetrahedron;
    QTVertex*      m_mateVertex;

public:
    QTFace();
    QTFace(QTTetrahedron* tetrahedron, QTVertex* vertex);
    QTFace(const QTFace& face);
    ~QTFace();

    QTTetrahedron* getTetrahedron() const;
    QTVertex*      getMateVertex() const;
    QTFace         getTwinFace() const;
    QTTetrahedron* getIncidentTetrahedron() const;

    void setTetrahedron(QTTetrahedron* tetrahedron);
    void setMateVertex(QTVertex* vertex);
    void setQTFace(QTTetrahedron* tetrahedron, QTVertex* vertex);

    QTFace& operator =(const QTFace& face);

    rg_FLAG operator==(const QTFace& face);


    rg_BOOL isThere(const QTEdge& edge) const;
    rg_FLAG getEdgeOrientation(const QTEdge& edge) const;

    ///////////////////////////////////////////////////////////////////////////////
    //    For Queres  
    void inquireBoundingVerticesInThisWorld(QTVertex** vertexArray) const;
    void inquireBoundingVerticesInThisWorld(rg_dList<QTVertex*>& vertexList) const;
    void inquireBoundingEdgesInThisWorld(rg_dList<QTEdge>& edgeList) const;
    //void inquireIncidentFacesInThisWorld(rg_dList<QTFace>& faceList);
    void inquireIncidentTetrahedraInThisWorld(rg_dList<QTTetrahedron*>& tetrahedronList) const;
    
};

} // namespace GeometryTier

} // namespace V


#endif
