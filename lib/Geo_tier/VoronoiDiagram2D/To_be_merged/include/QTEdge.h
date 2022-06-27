#ifndef QTEDGE_H_
#define QTEDGE_H_


#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {

class QTVertex;
class QTFace;

class QTEdge
{
private:
    int         m_ID;

    QTVertex*   m_startVertex;
    QTVertex*   m_endVertex;

    QTFace*     m_leftFace;
    QTFace*     m_rightFace;

    QTEdge*     m_leftHand;
    QTEdge*     m_rightHand;
    QTEdge*     m_leftLeg;
    QTEdge*     m_rightLeg;

public:
    QTEdge(void);
    QTEdge( const int& ID );
    QTEdge( const int& ID, QTVertex* const startVertex, QTVertex* const endVertex );
    QTEdge( const QTEdge& edge );
    ~QTEdge(void);

    inline int          getID() const           { return m_ID; }
    inline QTVertex*    getStartVertex() const  { return m_startVertex; }
    inline QTVertex*    getEndVertex() const    { return m_endVertex; }
    inline QTFace*      getLeftFace() const     { return m_leftFace; }
    inline QTFace*      getRightFace() const    { return m_rightFace; }
    inline QTEdge*      getLeftHand() const     { return m_leftHand; }
    inline QTEdge*      getRightHand() const    { return m_rightHand; }
    inline QTEdge*      getLeftLeg() const      { return m_leftLeg; }
    inline QTEdge*      getRightLeg() const     { return m_rightLeg; }

    inline void         setID( const int& ID ) { m_ID = ID; }
    inline void         setStartVertex( QTVertex* const startVertex ) { m_startVertex = startVertex; }
    inline void         setEndVertex( QTVertex* const endVertex ) { m_endVertex = endVertex; }
    inline void         setLeftFace( QTFace* const leftFace ) { m_leftFace = leftFace; }
    inline void         setRightFace( QTFace* const rightFace ) { m_rightFace = rightFace; }
    inline void         setLeftHand( QTEdge* const leftHand ) { m_leftHand = leftHand; }
    inline void         setRightHand( QTEdge* const rightHand ) { m_rightHand = rightHand; }
    inline void         setLeftLeg( QTEdge* const leftLeg ) { m_leftLeg = leftLeg; }
    inline void         setRightLeg( QTEdge* const rightLeg ) { m_rightLeg = rightLeg; }

    QTEdge&             operator=( const QTEdge& edge );
    bool                operator==( const QTEdge& edge ) const;

    void                getIncidentVertices( list<QTVertex*> incidentVertices );
    void                getIncidentFaces( list<QTFace*> incidentFaces );
    void                getAdjacentEdges( list<QTEdge*> adjacentEdges );
    
    bool                isGoingToInfinite() const;
};

}
}

#endif

