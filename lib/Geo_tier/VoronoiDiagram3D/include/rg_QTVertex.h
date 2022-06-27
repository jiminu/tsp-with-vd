#ifndef _QTVERTEX_H
#define _QTVERTEX_H

#include "Sphere.h"
#include "TopologicalEntity.h"

#include "rg_dList.h"
#include "Ball.h"


namespace V {

namespace GeometryTier {


class QTTetrahedron;
class VDCell;

class QTVertex : public TopologicalEntity
{
private:
    QTTetrahedron*  m_firstTetrahedron;

    Ball*           m_ball;

    rg_BOOL         m_visited;

    VDCell* m_vCell;

public:
    QTVertex();
    QTVertex(const rg_INT& id);
    QTVertex(const rg_INT& id, Ball* ball);
    QTVertex(const QTVertex& vertex);
    ~QTVertex();

    inline VDCell* getVCell() const { return m_vCell; }

    QTTetrahedron* getFirstTetrahedron() const;
    Sphere         getBall() const;
    void*          getProperty() const;
    rg_INT         getIDFromInput() const;
    Ball*          getBallProperty();
    
    inline rg_BOOL isVisited() const { return m_visited; }
    inline void    isVisited(const rg_BOOL& visited) { m_visited = visited; }
    rg_BOOL        isVirtual() const;

    void   setFirstTetrahedron(QTTetrahedron* tetrahedron);
    void   setBall(Ball* ball);

    QTVertex& operator =(const QTVertex& vertex);

    void inquireIncidentTetrahedra(rg_dList<QTTetrahedron*>& incidentTetrahedra);

    void connectVCell(   VDCell* v_cell);
    void disconnectVCell(VDCell* v_cell);
};

} // namespace GeometryTier

} // namespace V


#endif
