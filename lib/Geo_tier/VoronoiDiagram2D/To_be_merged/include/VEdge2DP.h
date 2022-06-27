#ifndef VEDGE2DP_H
#define VEDGE2DP_H

#include "rg_Point2D.h"



#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {

class VVertex2DP;
class VFace2DP;

enum STATUS_OF_VEDGEP
{
    RED_E_P,        // vertex in the 
    WHITE_E_P
};

class VEdge2DP
{
private:
    int             m_ID;
    VVertex2DP*      m_startVertex;
    VVertex2DP*      m_endVertex;
    VFace2DP*        m_leftFace;
    VFace2DP*        m_rightFace;
    VEdge2DP*        m_leftHand;
    VEdge2DP*        m_rightHand;
    VEdge2DP*        m_leftLeg;
    VEdge2DP*        m_rightLeg;
    STATUS_OF_VEDGEP m_status;

public:
    VEdge2DP();
    VEdge2DP( const int& ID );
    VEdge2DP( const int& ID, VVertex2DP* const startVertex, VVertex2DP* const endVertex );
    VEdge2DP( const VEdge2DP& edge );
    ~VEdge2DP();

    inline int              getID() const { return m_ID; }
    inline VVertex2DP*      getStartVertex() const { return m_startVertex; }
    inline VVertex2DP*      getEndVertex() const { return m_endVertex; }
    inline VFace2DP*        getLeftFace() const { return m_leftFace; }
    inline VFace2DP*        getRightFace() const { return m_rightFace; }
    inline VEdge2DP*        getLeftHand() const { return m_leftHand; }
    inline VEdge2DP*        getRightHand() const { return m_rightHand; }
    inline VEdge2DP*        getLeftLeg() const { return m_leftLeg; }
    inline VEdge2DP*        getRightLeg() const { return m_rightLeg; }
    inline STATUS_OF_VEDGEP getStatus() const { return m_status; }

    VFace2DP*               getMateFace(VVertex2DP* vertex) const;

    inline void             setID( const int& ID ) { m_ID = ID; }
    inline void             setStartVertex( VVertex2DP* const startVertex ) { m_startVertex = startVertex; }
    inline void             setEndVertex( VVertex2DP* const endVertex ) { m_endVertex = endVertex; }
    inline void             setLeftFace( VFace2DP* const leftFace ) { m_leftFace = leftFace; }
    inline void             setRightFace( VFace2DP* const rightFace ) { m_rightFace = rightFace; }
    inline void             setLeftHand( VEdge2DP* const leftHand ) { m_leftHand = leftHand; }
    inline void             setRightHand( VEdge2DP* const rightHand ) { m_rightHand = rightHand; }
    inline void             setLeftLeg( VEdge2DP* const leftLeg ) { m_leftLeg = leftLeg; }
    inline void             setRightLeg( VEdge2DP* const rightLeg ) { m_rightLeg = rightLeg; }
    inline void             setStatus( const STATUS_OF_VEDGEP& status ) { m_status = status; }

    VEdge2DP&                operator=( const VEdge2DP& edge );
    
    void                    getIncidentVVertices( list<VVertex2DP*>& incidentVVerticesList ) const;
    void                    getIncidentVFaces( list<VFace2DP*>& incidentVFacesList ) const;
    void                    getAdjacentVEdges( list<VEdge2DP*>& adjacentVEdgesList ) const;
    VVertex2DP*             getOppositeVVertex( VVertex2DP* vertex ) const;
    VFace2DP*               getOppositeVFace( VFace2DP* face ) const;


    void                    setTopology( VVertex2DP* const startVertex, VVertex2DP* const endVertex,
                                         VFace2DP* const leftFace, VFace2DP* const rightFace,
                                         VEdge2DP* const leftHand, VEdge2DP* const rightHand,
                                         VEdge2DP* const leftLeg, VEdge2DP* const rightLeg );
    void                    setTopology( VFace2DP* const leftFace, VFace2DP* const rightFace,
                                         VEdge2DP* const leftHand, VEdge2DP* const rightHand,
                                         VEdge2DP* const leftLeg, VEdge2DP* const rightLeg );
    void                    setTopology( VEdge2DP* const leftHand, VEdge2DP* const rightHand,
                                         VEdge2DP* const leftLeg, VEdge2DP* const rightLeg );

    bool    isInfinite() const;
    bool    isBounded() const;
    bool    isUnbounded() const;

    bool    isMemberOf( const VFace2DP* const face ) const;
    bool    isAdjacentTo( const VEdge2DP* const edge ) const;
    bool    isIncidentTo( const VVertex2DP* const vertex ) const;

};

} // namespace GeometryTier
} // namespace BULL2D



#endif


