#ifndef _QTEDGE_H
#define _QTEDGE_H

#include "TopologicalEntity.h"
#include "rg_dList.h"
#include "rg_QTFace.h"

namespace V {

namespace GeometryTier {


class QTVertex;
class QTTetrahedron;

class QTEdge : public TopologicalEntity
{
private:
    QTTetrahedron* m_tetrahedron;
    QTVertex*      m_startVertex;
    QTVertex*      m_endVertex;
    
public:
    QTEdge();
    QTEdge(QTTetrahedron* tetrahedron, QTVertex* stVertex, QTVertex* edVertex);
    QTEdge(const QTEdge& edge);
    ~QTEdge();

    QTTetrahedron* getTetrahedron() const;
    QTVertex*      getStartVertex() const;
    QTVertex*      getEndVertex() const;

	void setTetrahedron(QTTetrahedron* tetrahedron);
	void setStartVertex(QTVertex* stVertex);
	void setEndVertex(QTVertex* edVertex);
	void setEdge(QTTetrahedron* tetrahedron, QTVertex* stVertex, QTVertex* edVertex);


	QTEdge& operator =(const QTEdge& edge);
    rg_FLAG operator==(const QTEdge& edge);

    QTFace getCCWNextFace(const QTFace& currFace) const;
    QTFace getCWNextFace( const QTFace& currFace) const;
    QTTetrahedron* getCCWNextCell(QTTetrahedron* currCell = rg_NULL) const;
    QTTetrahedron* getCWNextCell( QTTetrahedron* currCell = rg_NULL) const;


    /////  For Queres  /////
    void inquireBoundingVerticesInThisWorld(rg_dList<QTVertex*>& vertexList);
    void inquireIncidentEdgesInThisWorld(rg_dList<QTEdge>& edgeList);
    void inquireIncidentFacesInThisWorld(rg_dList<QTFace>& faceList);
    void inquireIncidentTetrahedraInThisWorld(rg_dList<QTTetrahedron*>& tetrahedronList);
    

};

} // namespace GeometryTier

} // namespace V


#endif
