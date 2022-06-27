#ifndef _VDPARTIALEDGE_H
#define _VDPARTIALEDGE_H

#include "ConstForVoronoiDiagram3D.h"
#include "TopologicalEntity.h"

namespace V {

namespace GeometryTier {


class VDLoop;
class VDEdge;

class VDPartialEdge : public TopologicalEntity
{
private:
    VDLoop*        m_loop;
    VDEdge*        m_originalEdge;

    VDPartialEdge* m_nextPartEdgeInRadialCycle;
    VDPartialEdge* m_nextPartEdgeInLoop;
    VDPartialEdge* m_prevPartEdgeInLoop;

    rg_FLAG        m_orientationInLoop;

    rg_BOOL        m_visited;

public:
    //  constructor & deconstructor..
    VDPartialEdge();
    VDPartialEdge( const rg_INT&        ID, 
                         VDEdge*        oriEdge );
    VDPartialEdge( const rg_INT&        ID, 
                         VDLoop*        loop, 
                         VDEdge*        oriEdge,
                         VDPartialEdge* nextPartEdgeInRadialCycle,
                         VDPartialEdge* nextPartEdgeInLoop, 
                         VDPartialEdge* prevPartEdgeInLoop,
                   const rg_FLAG&       orientation );
    VDPartialEdge( const TopologicalEntity& aTopoEntity, 
                         VDLoop*            loop, 
                         VDEdge*            oriEdge,
                         VDPartialEdge*     nextPartEdgeInRadialCycle,
                         VDPartialEdge*     nextPartEdgeInLoop, 
                         VDPartialEdge*     prevPartEdgeInLoop,
                   const rg_FLAG&           orientation );
    VDPartialEdge( const VDPartialEdge& aPartialEdge );
    ~VDPartialEdge();

    //  get functions.. 
    inline rg_BOOL   isVisited()   const { return m_visited; }
    inline void      isVisited(const rg_BOOL& visited) { m_visited = visited; }


    VDLoop*        getLoop() const;
    VDEdge*        getOriginalEdge() const;
    
    VDPartialEdge* getNextPartialEdgeInRadialCycle() const;
    VDPartialEdge* getNextPartialEdgeInLoop() const;
    VDPartialEdge* getPreviousPartialEdgeInLoop() const;

    rg_FLAG        isRightOrientationInLoop() const;

    //  set functions..
    void setLoop( VDLoop* loop );
    void setOriginalEdge( VDEdge* oriEdge );

    void setNextPartialEdgeInRadialCycle( VDPartialEdge* nextPartEdge );
    void setNextPartialEdgeInLoop( VDPartialEdge* nextPartEdge );
    void setPreviousPartialEdgeInLoop( VDPartialEdge* prevPartEdge );

    void addNextPartialEdgeInLoop( const rg_FLAG& orientation, VDLoop* loop, VDPartialEdge* nextPartEdge );
    void addPreviousPartialEdgeInLoop( const rg_FLAG& orientation, VDLoop* loop, VDPartialEdge* prevPartEdge );

    void moveThisToNextPartialEdgeInLoopOf(const rg_FLAG& orientation, VDPartialEdge* fiducialPartEdge);
    void moveThisToPreviousPartialEdgeInLoopOf(const rg_FLAG& orientation, VDPartialEdge* fiducialPartEdge);

    void insertAsNextPartialEdgeInLoop( VDPartialEdge* nextPartEdge );
    void insertAsPreviousPartialEdgeInLoop( VDPartialEdge* prevPartEdge );
    void moveToNextOf( VDPartialEdge* prEdge );
    void moveToPreviousOf( VDPartialEdge* prEdge );

    void isRightOrientationInLoop( const rg_FLAG& orientation );

    void setPartialEdge(       VDLoop*        loop, 
                               VDEdge*        oriEdge,
                               VDPartialEdge* nextPartEdgeInRadialCycle,
                               VDPartialEdge* nextPartEdgeInLoop, 
                               VDPartialEdge* prevPartEdgeInLoop,
                         const rg_FLAG&       orientation );
    void setPartialEdge( const rg_INT&        ID,       
                               VDLoop*        loop, 
                               VDEdge*        oriEdge,
                               VDPartialEdge* nextPartEdgeInRadialCycle,
                               VDPartialEdge* nextPartEdgeInLoop, 
                               VDPartialEdge* prevPartEdgeInLoop,
                         const rg_FLAG&       orientation );


    //  operator overloading..
    VDPartialEdge& operator =( const VDPartialEdge& aPartialEdge );

    //  topological operators..

};

} // namespace GeometryTier

} // namespace V


#endif
