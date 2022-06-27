#ifndef QTVERTEX_H_
#define QTVERTEX_H_

#include "rg_Circle2D.h"


#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {

class QTEdge;
class QTFace;

class QTVertex
{
private:
    int             m_ID;
    QTEdge*         m_firstEdge;
    rg_Circle2D*    m_circle;

public:
    QTVertex(void);
    QTVertex( const int& ID );
    QTVertex( const int& ID, QTEdge* const firstEdge );
    QTVertex( const int& ID, rg_Circle2D* const circle );
    QTVertex( const int& ID, rg_Circle2D* const circle, QTEdge* const firstEdge );
    QTVertex( const QTVertex& vertex );
    ~QTVertex(void);

    inline int          getID() const           { return m_ID; }
    inline QTEdge*      getFirstEdge() const    { return m_firstEdge; }
    inline rg_Circle2D* getCircle() const       { return m_circle; }

    inline void         setID( const int& ID )                  { m_ID = ID; }
    inline void         setFirstEdge( QTEdge* const edge )      { m_firstEdge = edge; }
    inline void         setCircle( rg_Circle2D* const circle )  { m_circle = circle; }

    QTVertex            operator=( const QTVertex& vertex );
    bool                operator==( const QTVertex& vertex ) const;

    void                getIncidentEdges( list<QTEdge*>& incidentEdges ) const;
    void                getIncidentFaces( list<QTFace*>& incidentFaces ) const;
    void                getAdjacentVertices( list<QTVertex*>& adjacentVertices ) const;

    bool                isInfinite() const;
};

}
}

#endif


