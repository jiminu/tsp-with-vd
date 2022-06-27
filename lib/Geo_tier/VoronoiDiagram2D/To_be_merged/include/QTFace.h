#ifndef QTFACE_H_
#define QTFACE_H_



#include <list>
using namespace std;

namespace BULL2D{

namespace GeometryTier {

class QTEdge;
class QTVertex;

class QTFace
{
private:
    int         m_ID;
    QTEdge*     m_firstEdge;

public:
    QTFace(void);
    QTFace( const int& ID );
    QTFace( const int& ID, QTEdge* const firstEdge );
    QTFace( const QTFace& face );
    ~QTFace(void);

    inline int          getID() const { return m_ID; }
    inline QTEdge*      getFirstEdge() const { return m_firstEdge; }

    inline void         setID( const int& ID ) { m_ID = ID; }
    inline void         setFirstEdge( QTEdge* const edge ) { m_firstEdge = edge; }

    QTFace&             operator=( const QTFace& face );
    bool                operator==( const QTFace& face ) const;

    void                getBoundaryEdges( list<QTEdge*>& boundaryEdges ) const;
    void                getBoundaryVertices( list<QTVertex*>& boundaryVertices ) const;
    void                getAdjacentFaces( list<QTFace*>& adjacentFaces ) const;

    bool                isInfinite() const;
};

}
}

#endif



