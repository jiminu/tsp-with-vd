#ifndef _QTGATE_H
#define _QTGATE_H


#include "rg_Const.h"
#include "rg_dList.h"

namespace V {

namespace GeometryTier {


class QTVertex;
class QTTetrahedron;

class QTGate
{
private:
    QTVertex*                m_vertex[2];
    QTTetrahedron*           m_bigTetrahedron;
    rg_dList<QTTetrahedron*> m_smallTetrahedra;

public:
    QTGate();
    QTGate(QTVertex* vertex1, QTVertex* vertex2);
    QTGate(const QTGate& gate);
    ~QTGate();

    QTVertex*                 getVertex(const rg_INDEX& i);
    QTVertex**                getVertices();
    QTTetrahedron*            getBigTetrahedron();
    rg_dList<QTTetrahedron*>* getSmallTetrahedra();

    void setVertex(const rg_INDEX& i, QTVertex* vertex);
    void setVertices(QTVertex* vertex1, QTVertex* vertex2);
    void setVertices(QTVertex** vertices);
    void setBigTetrahedron(QTTetrahedron* bigTetrahedron);
    void addSmallTetrahedron(QTTetrahedron* smallTetrahedron);

    QTGate& operator =(const QTGate& gate);
    rg_FLAG operator==(const QTGate& gate);
};

} // namespace GeometryTier

} // namespace V


#endif
