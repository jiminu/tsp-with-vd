#ifndef VEDGE2D_H
#define VEDGE2D_H


#include "rg_RQBzCurve2D.h"

#include <list>
using namespace std;

namespace BULL2D {

namespace GeometryTier {


class VVertex2D;
class VFace2D;
class Generator2D;

enum STATUS_OF_VEDGE
{
    RED_E,        // edge in the red ocean
    PURPLE_E,     // intersecting edge
    WHITE_E       // default color
};

class VEdge2D
{
private:
    int                 m_ID;

    VVertex2D*          m_startVertex;
    VVertex2D*          m_endVertex;

    VFace2D*            m_leftFace;
    VFace2D*            m_rightFace;

    VEdge2D*            m_leftHand;
    VEdge2D*            m_rightHand;
    VEdge2D*            m_leftLeg;
    VEdge2D*            m_rightLeg;

    STATUS_OF_VEDGE     m_status;
    rg_RQBzCurve2D      m_geometry;

    bool                m_isAnomalyTestDone;
    bool                m_isAlreadyCandidateForFlippingInPhantomRemoval;

public:
    VEdge2D();
    VEdge2D( const int& ID );
    VEdge2D( const int& ID, VVertex2D* const startVertex, VVertex2D* const endVertex );
    VEdge2D( const VEdge2D& edge );
    ~VEdge2D();

    inline int              getID() const { return m_ID; }
    inline VVertex2D*       getStartVertex() const { return m_startVertex; }
    inline VVertex2D*       getEndVertex() const { return m_endVertex; }
    inline VFace2D*         getLeftFace() const { return m_leftFace; }
    inline VFace2D*         getRightFace() const { return m_rightFace; }
    inline VEdge2D*         getLeftHand() const { return m_leftHand; }
    inline VEdge2D*         getRightHand() const { return m_rightHand; }
    inline VEdge2D*         getLeftLeg() const { return m_leftLeg; }
    inline VEdge2D*         getRightLeg() const { return m_rightLeg; }
    inline STATUS_OF_VEDGE  getStatus() const { return m_status; }
    inline rg_RQBzCurve2D   getGeometry() const { return m_geometry; }
    //JKKIM_ADDED
    VFace2D*                getMateFace(VVertex2D* vertex) const;


    inline void             setID( const int& ID ) { m_ID = ID; }
    inline void             setStartVertex( VVertex2D* const startVertex ) { m_startVertex = startVertex; }
    inline void             setEndVertex( VVertex2D* const endVertex ) { m_endVertex = endVertex; }
    inline void             setLeftFace( VFace2D* const leftFace ) { m_leftFace = leftFace; }
    inline void             setRightFace( VFace2D* const rightFace ) { m_rightFace = rightFace; }
    inline void             setLeftHand( VEdge2D* const leftHand ) { m_leftHand = leftHand; }
    inline void             setRightHand( VEdge2D* const rightHand ) { m_rightHand = rightHand; }
    inline void             setLeftLeg( VEdge2D* const leftLeg ) { m_leftLeg = leftLeg; }
    inline void             setRightLeg( VEdge2D* const rightLeg ) { m_rightLeg = rightLeg; }
    inline void             setStatus( const STATUS_OF_VEDGE& status ) { m_status = status; }
    void                    setGeometry(const rg_Point2D& sp, 
                                        const rg_Point2D& tvs, 
                                        const rg_Point2D& ep, 
                                        const rg_Point2D& tve, 
                                        const rg_Point2D& passPt);
    void                    setGeometryLine();

    VEdge2D&                operator=( const VEdge2D& edge );
    bool                    operator==( const VEdge2D& edge ) const;
    
    void                    getIncidentVVertices( list<VVertex2D*>& incidentVVerticesList ) const;
    void                    getIncidentVFaces( list<VFace2D*>& incidentVFacesList ) const;
    void                    getAdjacentVEdges( list<VEdge2D*>& adjacentVEdgesList ) const;
    VVertex2D*              getOppositVVertex( VVertex2D* vertex ) const;
    VFace2D*                getOppositeVFace( VFace2D* face ) const;

    void                    setTopology( VVertex2D* const startVertex, VVertex2D* const endVertex,
                                         VFace2D* const leftFace, VFace2D* const rightFace,
                                         VEdge2D* const leftHand, VEdge2D* const rightHand,
                                         VEdge2D* const leftLeg, VEdge2D* const rightLeg );
    void                    setTopology( VFace2D* const leftFace, VFace2D* const rightFace,
                                         VEdge2D* const leftHand, VEdge2D* const rightHand,
                                         VEdge2D* const leftLeg, VEdge2D* const rightLeg );
    void                    setTopology( VEdge2D* const leftHand, VEdge2D* const rightHand,
                                         VEdge2D* const leftLeg, VEdge2D* const rightLeg );

    bool    isInfinite() const;
    bool    isBounded() const;
    bool    isUnBounded() const;
    void    flip();

    inline bool     isAnomalyTestDone() const { return m_isAnomalyTestDone; }
    inline void     setTrueAnomalyTestDone()  { m_isAnomalyTestDone = true; }
    inline void     setFalseAnomalyTestDone() { m_isAnomalyTestDone = false; }
    inline bool     isAlreadyCandidateForFlippingInPhantomRemoval() const { return m_isAlreadyCandidateForFlippingInPhantomRemoval; }
    inline void     setTrueCandidateForFlippingInPhantomRemoval()  { m_isAlreadyCandidateForFlippingInPhantomRemoval = true; }
    inline void     setFalseCandidateForFlippingInPhantomRemoval() { m_isAlreadyCandidateForFlippingInPhantomRemoval = false; }




    bool            isMemberOf( const VFace2D* const face ) const;
    bool            isAdjacentTo( const VEdge2D* const edge ) const;
    bool            isIncidentTo( const VVertex2D* const vertex ) const;

    void            get4GeneratorsDefineEdgeEquationAndStartAndEndVertices( list<Generator2D*>& generatorList );
};

}
}

#endif


