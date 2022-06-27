#include "VDPartialEdge.h"
using namespace V::GeometryTier;


///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
VDPartialEdge::VDPartialEdge()
: m_loop(                      rg_NULL ), 
  m_originalEdge(              rg_NULL ),
  m_nextPartEdgeInRadialCycle( rg_NULL ),
  m_orientationInLoop(         rg_UNKNOWN )
{
    m_nextPartEdgeInLoop = this;
    m_prevPartEdgeInLoop = this;

    m_visited = rg_FALSE;
}

VDPartialEdge::VDPartialEdge( const rg_INT&        ID, 
                                    VDEdge*        oriEdge )
: TopologicalEntity(           ID ),
  m_loop(                      rg_NULL ), 
  m_originalEdge(              oriEdge ),
  m_nextPartEdgeInRadialCycle( rg_NULL ),
  m_nextPartEdgeInLoop(        rg_NULL ),
  m_prevPartEdgeInLoop(        rg_NULL ),
  m_orientationInLoop(         rg_UNKNOWN )
{
    m_visited = rg_FALSE;
}

VDPartialEdge::VDPartialEdge( const rg_INT&        ID, 
                                    VDLoop*        loop, 
                                    VDEdge*        oriEdge,
                                    VDPartialEdge* nextPartEdgeInRadialCycle,
                                    VDPartialEdge* nextPartEdgeInLoop, 
                                    VDPartialEdge* prevPartEdgeInLoop,
                              const rg_FLAG&       orientation )
: TopologicalEntity(           ID ),
  m_loop(                      loop ), 
  m_originalEdge(              oriEdge ),
  m_nextPartEdgeInRadialCycle( nextPartEdgeInRadialCycle ),
  m_nextPartEdgeInLoop(        nextPartEdgeInLoop ),
  m_prevPartEdgeInLoop(        prevPartEdgeInLoop ),
  m_orientationInLoop (        orientation )
{
    m_visited = rg_FALSE;
}

VDPartialEdge::VDPartialEdge( const TopologicalEntity& aTopoEntity, 
                                    VDLoop*            loop, 
                                    VDEdge*            oriEdge,
                                    VDPartialEdge*     nextPartEdgeInRadialCycle,
                                    VDPartialEdge*     nextPartEdgeInLoop, 
                                    VDPartialEdge*     prevPartEdgeInLoop,
                              const rg_FLAG&           orientation )
: TopologicalEntity(           aTopoEntity ),
  m_loop(                      loop ), 
  m_originalEdge(              oriEdge ),
  m_nextPartEdgeInRadialCycle( nextPartEdgeInRadialCycle ),
  m_nextPartEdgeInLoop(        nextPartEdgeInLoop ),
  m_prevPartEdgeInLoop(        prevPartEdgeInLoop ),
  m_orientationInLoop(        orientation )
{
    m_visited = rg_FALSE;
}

VDPartialEdge::VDPartialEdge( const VDPartialEdge& aPartialEdge )
: TopologicalEntity(           aPartialEdge ),
  m_loop(                      aPartialEdge.m_loop ), 
  m_originalEdge(              aPartialEdge.m_originalEdge ),
  m_nextPartEdgeInRadialCycle( aPartialEdge.m_nextPartEdgeInRadialCycle ),
  m_nextPartEdgeInLoop(        aPartialEdge.m_nextPartEdgeInLoop ),
  m_prevPartEdgeInLoop(        aPartialEdge.m_prevPartEdgeInLoop ),
  m_orientationInLoop(         aPartialEdge.m_orientationInLoop )
{
    m_visited = aPartialEdge.m_visited;
}

VDPartialEdge::~VDPartialEdge()
{
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
VDLoop* VDPartialEdge::getLoop() const
{
    return m_loop;
}

VDEdge* VDPartialEdge::getOriginalEdge() const
{
    return m_originalEdge;
}


VDPartialEdge* VDPartialEdge::getNextPartialEdgeInRadialCycle() const
{
    return m_nextPartEdgeInRadialCycle;
}

VDPartialEdge* VDPartialEdge::getNextPartialEdgeInLoop() const
{
    return m_nextPartEdgeInLoop;
}

VDPartialEdge* VDPartialEdge::getPreviousPartialEdgeInLoop() const
{
    return m_prevPartEdgeInLoop;
}


rg_FLAG VDPartialEdge::isRightOrientationInLoop() const
{
    return m_orientationInLoop;
}


///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void VDPartialEdge::setLoop( VDLoop* loop )
{
    m_loop = loop;
}

void VDPartialEdge::setOriginalEdge( VDEdge* oriEdge )
{
    m_originalEdge = oriEdge;
}


void VDPartialEdge::setNextPartialEdgeInRadialCycle( VDPartialEdge* nextPartEdge )
{
    m_nextPartEdgeInRadialCycle = nextPartEdge;
}

void VDPartialEdge::setNextPartialEdgeInLoop( VDPartialEdge* nextPartEdge )
{
    m_nextPartEdgeInLoop = nextPartEdge;
}

void VDPartialEdge::setPreviousPartialEdgeInLoop( VDPartialEdge* prevPartEdge )
{
    m_prevPartEdgeInLoop = prevPartEdge;
}

void VDPartialEdge::addNextPartialEdgeInLoop( const rg_FLAG& orientation, VDLoop* loop, VDPartialEdge* nextPartEdge )
{
    nextPartEdge->m_orientationInLoop = orientation;
    nextPartEdge->m_loop              = loop;

    nextPartEdge->m_nextPartEdgeInLoop = m_nextPartEdgeInLoop;
    nextPartEdge->m_prevPartEdgeInLoop = this;

    m_nextPartEdgeInLoop->m_prevPartEdgeInLoop = nextPartEdge;
    this->m_nextPartEdgeInLoop                 = nextPartEdge;
}

void VDPartialEdge::addPreviousPartialEdgeInLoop( const rg_FLAG& orientation, VDLoop* loop, VDPartialEdge* prevPartEdge )
{
    prevPartEdge->m_orientationInLoop = orientation;
    prevPartEdge->m_loop              = loop;

    prevPartEdge->m_nextPartEdgeInLoop = this;
    prevPartEdge->m_prevPartEdgeInLoop = m_prevPartEdgeInLoop;

    m_prevPartEdgeInLoop->m_nextPartEdgeInLoop = prevPartEdge;
    this->m_prevPartEdgeInLoop                 = prevPartEdge;
}

void VDPartialEdge::insertAsNextPartialEdgeInLoop( VDPartialEdge* nextPartEdge )
{
    nextPartEdge->m_nextPartEdgeInLoop = m_nextPartEdgeInLoop;
    nextPartEdge->m_prevPartEdgeInLoop = this;

    m_nextPartEdgeInLoop->m_prevPartEdgeInLoop = nextPartEdge;
    this->m_nextPartEdgeInLoop                 = nextPartEdge;
}

void VDPartialEdge::insertAsPreviousPartialEdgeInLoop( VDPartialEdge* prevPartEdge )
{
    prevPartEdge->m_nextPartEdgeInLoop = this;
    prevPartEdge->m_prevPartEdgeInLoop = m_prevPartEdgeInLoop;

    m_prevPartEdgeInLoop->m_nextPartEdgeInLoop = prevPartEdge;
    this->m_prevPartEdgeInLoop                 = prevPartEdge;
}

void VDPartialEdge::moveToNextOf( VDPartialEdge* prEdge )
{
    if ( this == prEdge->m_nextPartEdgeInLoop )
        return;

    m_prevPartEdgeInLoop->m_nextPartEdgeInLoop = m_nextPartEdgeInLoop;
    m_nextPartEdgeInLoop->m_prevPartEdgeInLoop = m_prevPartEdgeInLoop;

    
    m_nextPartEdgeInLoop = prEdge->m_nextPartEdgeInLoop;
    m_prevPartEdgeInLoop = prEdge;

    m_nextPartEdgeInLoop->m_prevPartEdgeInLoop = this;
    prEdge->m_nextPartEdgeInLoop     = this;

}

void VDPartialEdge::moveToPreviousOf( VDPartialEdge* prEdge )
{
    if ( this == prEdge->m_prevPartEdgeInLoop )
        return;

    m_prevPartEdgeInLoop->m_nextPartEdgeInLoop = m_nextPartEdgeInLoop;
    m_nextPartEdgeInLoop->m_prevPartEdgeInLoop = m_prevPartEdgeInLoop;

    
    m_nextPartEdgeInLoop = prEdge;
    m_prevPartEdgeInLoop = prEdge->m_prevPartEdgeInLoop;

    m_prevPartEdgeInLoop->m_nextPartEdgeInLoop = this;
    prEdge->m_prevPartEdgeInLoop               = this;
}

void VDPartialEdge::moveThisToNextPartialEdgeInLoopOf(const rg_FLAG& orientation, VDPartialEdge* fiducialPartEdge)
{
    m_orientationInLoop = orientation;

    if ( this == fiducialPartEdge->m_nextPartEdgeInLoop )
        return;

    m_prevPartEdgeInLoop->m_nextPartEdgeInLoop = m_nextPartEdgeInLoop;
    m_nextPartEdgeInLoop->m_prevPartEdgeInLoop = m_prevPartEdgeInLoop;

    
    m_nextPartEdgeInLoop = fiducialPartEdge->m_nextPartEdgeInLoop;
    m_prevPartEdgeInLoop = fiducialPartEdge;

    m_nextPartEdgeInLoop->m_prevPartEdgeInLoop = this;
    fiducialPartEdge->m_nextPartEdgeInLoop     = this;
}

void VDPartialEdge::moveThisToPreviousPartialEdgeInLoopOf(const rg_FLAG& orientation, VDPartialEdge* fiducialPartEdge)
{
    m_orientationInLoop = orientation;

    if ( this == fiducialPartEdge->m_prevPartEdgeInLoop )
        return;

    m_prevPartEdgeInLoop->m_nextPartEdgeInLoop = m_nextPartEdgeInLoop;
    m_nextPartEdgeInLoop->m_prevPartEdgeInLoop = m_prevPartEdgeInLoop;

    
    m_nextPartEdgeInLoop = fiducialPartEdge;
    m_prevPartEdgeInLoop = fiducialPartEdge->m_prevPartEdgeInLoop;

    m_prevPartEdgeInLoop->m_nextPartEdgeInLoop = this;
    fiducialPartEdge->m_prevPartEdgeInLoop     = this;
}


void VDPartialEdge::isRightOrientationInLoop( const rg_FLAG& orientation )
{
    m_orientationInLoop = orientation;
}


void VDPartialEdge::setPartialEdge(       VDLoop*        loop, 
                                          VDEdge*        oriEdge,
                                          VDPartialEdge* nextPartEdgeInRadialCycle,
                                          VDPartialEdge* nextPartEdgeInLoop, 
                                          VDPartialEdge* prevPartEdgeInLoop,
                                    const rg_FLAG&       orientation )
{
    m_loop                      = loop;
    m_originalEdge              = oriEdge;

    m_nextPartEdgeInRadialCycle = nextPartEdgeInRadialCycle;
    m_nextPartEdgeInLoop        = nextPartEdgeInLoop;
    m_prevPartEdgeInLoop        = prevPartEdgeInLoop;

    m_orientationInLoop         = orientation;
}

void VDPartialEdge::setPartialEdge( const rg_INT&        ID,       
                                          VDLoop*        loop, 
                                          VDEdge*        oriEdge,
                                          VDPartialEdge* nextPartEdgeInRadialCycle,
                                          VDPartialEdge* nextPartEdgeInLoop, 
                                          VDPartialEdge* prevPartEdgeInLoop,
                                    const rg_FLAG&       orientation )
{
    m_ID                        = ID;

    m_loop                      = loop;
    m_originalEdge              = oriEdge;

    m_nextPartEdgeInRadialCycle = nextPartEdgeInRadialCycle;
    m_nextPartEdgeInLoop        = nextPartEdgeInLoop;
    m_prevPartEdgeInLoop        = prevPartEdgeInLoop;

    m_orientationInLoop         = orientation;
}


///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
VDPartialEdge& VDPartialEdge::operator =( const VDPartialEdge& aPartialEdge )
{
    if ( this == &aPartialEdge )
        return *this;

    m_ID                        = aPartialEdge.m_ID;

    m_loop                      = aPartialEdge.m_loop;
    m_originalEdge              = aPartialEdge.m_originalEdge;

    m_nextPartEdgeInRadialCycle = aPartialEdge.m_nextPartEdgeInRadialCycle;
    m_nextPartEdgeInLoop        = aPartialEdge.m_nextPartEdgeInLoop;
    m_prevPartEdgeInLoop        = aPartialEdge.m_prevPartEdgeInLoop;

    m_orientationInLoop         = aPartialEdge.m_orientationInLoop;


    m_visited = aPartialEdge.m_visited;

    return *this;
}


///////////////////////////////////////////////////////////////////////////////
//
//  topological operators..

