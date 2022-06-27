#include "VDLoop.h"

#include "VDFace.h"
#include "VDPartialEdge.h"
#include "VDEdge.h"
#include "VDVertex.h"
using namespace V::GeometryTier;


///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
VDLoop::VDLoop()
: m_face(rg_NULL), m_partEdge(rg_NULL), m_isOuterLoop(rg_UNKNOWN), m_numOfBoundingEdges(0)
{
}
    
VDLoop::VDLoop(const rg_INT& ID, VDFace* face)
: TopologicalEntity(ID), m_face(face), m_partEdge(rg_NULL), m_isOuterLoop(rg_UNKNOWN), m_numOfBoundingEdges(0)
{
}

VDLoop::VDLoop(const rg_INT& ID, VDFace* face, const rg_FLAG& bIsOuterLoop )
: TopologicalEntity(ID), m_face(face), m_partEdge(rg_NULL), m_isOuterLoop(bIsOuterLoop), m_numOfBoundingEdges(0)
{
}

VDLoop::VDLoop(const rg_INT& ID, 
               VDFace* face, VDPartialEdge* partEdge, const rg_FLAG& bIsOuterLoop)
: TopologicalEntity(ID), m_face(face), m_partEdge(partEdge), m_isOuterLoop(bIsOuterLoop), m_numOfBoundingEdges(0)
{
}

VDLoop::VDLoop(const TopologicalEntity& aTopoEntity, 
               VDFace* face, VDPartialEdge* partEdge, const rg_FLAG& bIsOuterLoop)
: TopologicalEntity(aTopoEntity), m_face(face), m_partEdge(partEdge), m_isOuterLoop(bIsOuterLoop), m_numOfBoundingEdges(0)
{
}

VDLoop::VDLoop(const VDLoop& aLoop)
: TopologicalEntity(aLoop), 
  m_face(aLoop.m_face), m_partEdge(aLoop.m_partEdge), m_isOuterLoop(aLoop.m_isOuterLoop), m_numOfBoundingEdges(aLoop.m_numOfBoundingEdges)
{
}

VDLoop::~VDLoop()
{
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
VDFace* VDLoop::getFace() const
{
    return m_face;
}

VDPartialEdge* VDLoop::getPartialEdge() const
{
    return m_partEdge;
}

rg_INT VDLoop::getNumOfBoundingEdges() const
{
    return m_numOfBoundingEdges;
}

rg_FLAG VDLoop::isOuterLoop() const
{
    return m_isOuterLoop;
}

rg_FLAG VDLoop::isClosedLoop() const
{
    VDPartialEdge* currPartEdge     = m_partEdge;

    do
    {
        if ( currPartEdge->getOriginalEdge()->isBoundedEdge() == rg_FALSE )
            return rg_FALSE;

        currPartEdge = currPartEdge->getNextPartialEdgeInLoop();
    } while ( currPartEdge != m_partEdge );

    return rg_TRUE;
}



rg_BOOL VDLoop::isSingleEdge() const
{
    if ( m_partEdge->getOriginalEdge()->getStartVertex() == rg_NULL ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}

///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void VDLoop::setFace(VDFace* face)
{
    m_face = face;
}

void VDLoop::setPartialEdge(VDPartialEdge* partEdge)
{
    m_partEdge = partEdge;
}

void VDLoop::isOuterLoop(const rg_FLAG& bIsOuterLoop)
{
    m_isOuterLoop = bIsOuterLoop;
}


void VDLoop::setLoop(VDFace* face, VDPartialEdge* partEdge, const rg_FLAG& bIsOuterLoop)
{
    m_face        = face;
    m_partEdge    = partEdge;
    m_isOuterLoop = bIsOuterLoop;
}

void VDLoop::setLoop(const rg_INT& ID, VDFace* face, VDPartialEdge* partEdge, const rg_FLAG& bIsOuterLoop)
{
    m_ID          = ID;

    m_face        = face;
    m_partEdge    = partEdge;
    m_isOuterLoop = bIsOuterLoop;
}

rg_FLAG VDLoop::addPartialEdgeByEdgeTracing(VDPartialEdge *const partEdgeToAdd)
{
    ////////////////////////////////////////////////////////////
    //
    //  face와 face의 outer loop는 gate의 edge가 정의될 때 같이 정의된다.
    //    m_partEdge는 loop를 구성하는 partial edge 중 
    //    orientation이 참인 방향으로 제일 끝에 있는 partial edge이다.

    m_numOfBoundingEdges++;

    if ( m_partEdge == rg_NULL )  
    {
        //  loop에 처음 partial edge가 할당될 때.
        partEdgeToAdd->isRightOrientationInLoop( rg_TRUE );
        partEdgeToAdd->setLoop( this );

        m_partEdge = partEdgeToAdd;

        return rg_TRUE;
    }
    else
    {
        VDEdge* edgeToAdd    = partEdgeToAdd->getOriginalEdge();

        VDPartialEdge* currentPartEdge = m_partEdge;

        //  추가할 partial edge를 loop를 구성하는 partial edge들과 
        //  vertex를 검사하여 적합한 위치에 넣는다. 
        do 
        {
            VDEdge* currentEdge = currentPartEdge->getOriginalEdge();

            //  currentEdge : =====>
            //  edgeToAdd   : ----->

            //   <-----o=====>
            if ( currentEdge->getStartVertex() == edgeToAdd->getStartVertex() )
            {
                //   -----loop--->
                //   <-----o=====>
                if ( currentPartEdge->isRightOrientationInLoop() == rg_TRUE )
                {
                    currentPartEdge->addPreviousPartialEdgeInLoop( rg_FALSE, this, partEdgeToAdd );
                    return rg_TRUE;               
                }
                //   <----loop----
                //   <-----o=====>
                else if ( currentPartEdge->isRightOrientationInLoop() == rg_FALSE )
                {
                    currentPartEdge->addNextPartialEdgeInLoop( rg_TRUE, this, partEdgeToAdd );

                    if ( currentPartEdge == m_partEdge )
                        m_partEdge = partEdgeToAdd;

                    return rg_TRUE;               
                }
                //   <-----o=====>
                else // if ( currentPartEdge->isRightOrientationInLoop() == rg_UNKNOWN )
                {
                    currentPartEdge->addPreviousPartialEdgeInLoop( rg_UNKNOWN, this, partEdgeToAdd );
                    return rg_FALSE;               
                }
            }
            //   ----->o=====>
            else if ( currentEdge->getStartVertex() == edgeToAdd->getEndVertex() )
            {
                //   -----loop--->
                //   ----->o=====>
                if ( currentPartEdge->isRightOrientationInLoop() == rg_TRUE )
                {
                    currentPartEdge->addPreviousPartialEdgeInLoop( rg_TRUE, this, partEdgeToAdd );
                    return rg_TRUE;               
                }
                //   <----loop----
                //   ----->o=====>
                else if ( currentPartEdge->isRightOrientationInLoop() == rg_FALSE )
                {
                    currentPartEdge->addNextPartialEdgeInLoop( rg_FALSE, this, partEdgeToAdd );

                    if ( currentPartEdge == m_partEdge )
                        m_partEdge = partEdgeToAdd;

                    return rg_TRUE;               
                }
                //   ----->o=====>
                else // if ( currentPartEdge->isRightOrientationInLoop() == rg_UNKNOWN )
                {
                    currentPartEdge->addPreviousPartialEdgeInLoop( rg_UNKNOWN, this, partEdgeToAdd );
                    return rg_FALSE;               
                }
            }
            //   =====>o----->
            else if ( currentEdge->getEndVertex() == edgeToAdd->getStartVertex() )
            {
                //   -----loop--->
                //   =====>o----->
                if ( currentPartEdge->isRightOrientationInLoop() == rg_TRUE )
                {
                    currentPartEdge->addNextPartialEdgeInLoop( rg_TRUE, this, partEdgeToAdd );

                    //  처음 만들어진 partial edge의 orientation이 정방향이고,
                    //  이후, 이 partial edge에 정방향으로 연결된 partial edge만이 
                    //  m_partEdge가 된다.
                    m_partEdge = partEdgeToAdd;
                    return rg_TRUE;               
                }
                //   <----loop----
                //   =====>o----->
                else if ( currentPartEdge->isRightOrientationInLoop() == rg_FALSE )
                {
                    currentPartEdge->addPreviousPartialEdgeInLoop( rg_FALSE, this, partEdgeToAdd );
                    return rg_TRUE;               
                }
                //   =====>o----->
                else // if ( currentPartEdge->isRightOrientationInLoop() == rg_UNKNOWN )
                {
                    currentPartEdge->addNextPartialEdgeInLoop( rg_UNKNOWN, this, partEdgeToAdd );
                    return rg_FALSE;               
                }
            }
            //   =====>o<-----
            else if ( currentEdge->getEndVertex() == edgeToAdd->getEndVertex() )
            {
                //   -----loop--->
                //   =====>o<-----
                if ( currentPartEdge->isRightOrientationInLoop() == rg_TRUE )
                {
                    currentPartEdge->addNextPartialEdgeInLoop( rg_FALSE, this, partEdgeToAdd );

                    if ( currentPartEdge == m_partEdge )
                        m_partEdge = partEdgeToAdd;

                    return rg_TRUE;               
                }
                //   <----loop----
                //   =====>o<-----
                else if ( currentPartEdge->isRightOrientationInLoop() == rg_FALSE )
                {
                    currentPartEdge->addPreviousPartialEdgeInLoop( rg_TRUE, this, partEdgeToAdd );
                    return rg_TRUE;               
                }
                //   =====>o<-----
                else // if ( currentPartEdge->isRightOrientationInLoop() == rg_UNKNOWN )
                {
                    currentPartEdge->addNextPartialEdgeInLoop( rg_UNKNOWN, this, partEdgeToAdd );
                    return rg_FALSE;               
                }
            }
            else
            {
                currentPartEdge = currentPartEdge->getNextPartialEdgeInLoop();
            }                
        } while ( currentPartEdge != m_partEdge );


        //  기존의 loop를 구성하는 partial edge들과 연결되지 않을 때, 
        //  추가할 partial edge의 orientation을 당장은 알 수 없고,
        //  m_partEdge의 다음에 놓는다.
        //  즉, orientation이 결정되지 않은 partial edge들은 모두 
        //  m_partEdge의 다음에 놓인다.
        m_partEdge->addNextPartialEdgeInLoop( rg_UNKNOWN, this, partEdgeToAdd );

        return rg_FALSE;
    }
}

void VDLoop::addPartialEdgeToLoop(VDPartialEdge *const partEdgeToAdd)
{
    m_numOfBoundingEdges++;
    partEdgeToAdd->setLoop( this );

    if ( m_partEdge == rg_NULL )  
    {
        m_partEdge = partEdgeToAdd;
    }
    else
    {
        m_partEdge->insertAsNextPartialEdgeInLoop( partEdgeToAdd );

        /*
        VDEdge* edgeToAdd    = partEdgeToAdd->getOriginalEdge();

        VDPartialEdge* startPrEdge = m_partEdge;
        VDPartialEdge* currPrEdge  = startPrEdge;

        rg_FLAG isAdded = rg_FALSE;
        do 
        {
            VDEdge* currEdge = currPrEdge->getOriginalEdge();

            //  currentEdge : =====
            //  edgeToAdd   : -----

            //   -----o=====>
            if (    currEdge->getStartVertex() == edgeToAdd->getStartVertex() 
                 || currEdge->getStartVertex() == edgeToAdd->getEndVertex()  )
            {
                //   -----loop--->
                //   ------o=====>
                if ( currPrEdge->isRightOrientationInLoop() == rg_TRUE )  
                {
                    currPrEdge->insertAsPreviousPartialEdgeInLoop( partEdgeToAdd );
                    isAdded = rg_TRUE;
                    break;
                }
                //   <----loop----
                //   ------o=====>
                else  
                {
                    currPrEdge->insertAsNextPartialEdgeInLoop( partEdgeToAdd );
                    isAdded = rg_TRUE;
                    break;
                }
            }
            //   =====>o------
            else if (    currEdge->getEndVertex() == edgeToAdd->getStartVertex() 
                      || currEdge->getEndVertex() == edgeToAdd->getEndVertex()  )
            {
                //   -----loop--->
                //   =====>o------
                if ( currPrEdge->isRightOrientationInLoop() == rg_TRUE )  
                {
                    currPrEdge->insertAsNextPartialEdgeInLoop( partEdgeToAdd );
                    isAdded = rg_TRUE;
                    break;
                }
                //   <----loop----
                //   =====>o------
                else  
                {
                    currPrEdge->insertAsPreviousPartialEdgeInLoop( partEdgeToAdd );
                    isAdded = rg_TRUE;
                    break;
                }
            }
            else
            {
                currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
            }                

        } while ( startPrEdge != currPrEdge );
        
        if ( isAdded )
        {
            if ( m_partEdge->isRightOrientationInLoop() == rg_TRUE )
            {
                if (    m_partEdge->getOriginalEdge()->getStartVertex() == partEdgeToAdd->getOriginalEdge()->getStartVertex() 
                     || m_partEdge->getOriginalEdge()->getStartVertex() == partEdgeToAdd->getOriginalEdge()->getEndVertex() )
                    m_partEdge = partEdgeToAdd;
            }
            else
            {
                if (    m_partEdge->getOriginalEdge()->getEndVertex() == partEdgeToAdd->getOriginalEdge()->getStartVertex() 
                     || m_partEdge->getOriginalEdge()->getEndVertex() == partEdgeToAdd->getOriginalEdge()->getEndVertex() )
                    m_partEdge = partEdgeToAdd;
            }
        }
        else
        {
            m_partEdge->insertAsPreviousPartialEdgeInLoop( partEdgeToAdd );
        }
        */
    }
}



void VDLoop::addPartialEdgeOnInfinityToLoop(VDPartialEdge *const partEdgeToAdd)
{
    m_numOfBoundingEdges++;
    partEdgeToAdd->setLoop( this );

    if ( m_partEdge == rg_NULL )  
    {
        m_partEdge = partEdgeToAdd;
    }
    else
    {
        VDPartialEdge* currPartEdge     = m_partEdge;

        do
        {
            if ( currPartEdge->getOriginalEdge()->isBoundedEdge() == rg_FALSE )
                break;

            currPartEdge = currPartEdge->getNextPartialEdgeInLoop();
        } while ( currPartEdge != m_partEdge );

        currPartEdge->insertAsNextPartialEdgeInLoop( partEdgeToAdd );
    }
}



rg_FLAG VDLoop::defineOrderOfPartialEdgesInLoop( rg_dList<VDPartialEdge*>& multipleOuterLoop )
{
    VDPartialEdge* startPartEdge    = m_partEdge;
    VDPartialEdge* currPartEdge     = startPartEdge;

    //  loop내의 partial edge들이 모두 방향성을 가지고 연결되었는가?
    rg_FLAG bIsAllEdgeOrientable = rg_TRUE;
    do 
    {
        if ( currPartEdge->isRightOrientationInLoop() == rg_UNKNOWN )
        {
            bIsAllEdgeOrientable = rg_FALSE;
            break;
        }
        currPartEdge = currPartEdge->getNextPartialEdgeInLoop();
    } while ( currPartEdge != startPartEdge );

    if ( bIsAllEdgeOrientable == rg_TRUE )
        return rg_TRUE;


    //  unbounded edge의 수가 3 이상일 경우, 다른 방법을 통해 방향성을 결정해야함.
    //  현재 고려하지 못했음. by Youngsong Cho, 2004-4-16.
    rg_INT numOfInfiniteEdges = 0;
    currPartEdge     = m_partEdge;
    do
    {
        if ( currPartEdge->getOriginalEdge()->isBoundedEdge() == rg_FALSE )
            numOfInfiniteEdges++;

        currPartEdge = currPartEdge->getNextPartialEdgeInLoop();
    } while ( currPartEdge != m_partEdge );

    if (   numOfInfiniteEdges > 2
        || numOfInfiniteEdges == m_numOfBoundingEdges )
        return rg_FALSE;
    


    //  lastCounterCWPartEdge: m_partEdge를 기준으로 CCW 방향(=next)으로 방향을 가진 마지막 partial edge
    //  lastCWPartEdge: m_partEdge를 기준으로 CW 방향(=previous)으로 방향을 가진 마지막 partial edge
    VDPartialEdge* lastCounterCWPartEdge = m_partEdge;
    VDPartialEdge* lastCWPartEdge        = m_partEdge;

    currPartEdge = m_partEdge->getPreviousPartialEdgeInLoop();
    while ( currPartEdge->isRightOrientationInLoop() != rg_UNKNOWN )
    {
        lastCWPartEdge = currPartEdge;

        currPartEdge = currPartEdge->getPreviousPartialEdgeInLoop();
    }

    
    //  방향을 가지지 않는 partial edge는 lastCounterCWPartEdge와 lastCWPartEdge사이에 놓인다.
    //  아래의 조건은 경우에 따라 무한 루프를 발생시킬 수 있음.
    while ( lastCounterCWPartEdge->getNextPartialEdgeInLoop() != lastCWPartEdge )
    {
        VDEdge* lastCounterCWEdge = lastCounterCWPartEdge->getOriginalEdge();
        VDEdge* lastCWEdge        = lastCWPartEdge->getOriginalEdge();

        rg_FLAG bIsMultiOuterLoop = rg_FALSE;

        if (    lastCounterCWEdge->getEndVertex()   == lastCWEdge->getStartVertex() 
             || lastCounterCWEdge->getEndVertex()   == lastCWEdge->getEndVertex() 
             || lastCounterCWEdge->getStartVertex() == lastCWEdge->getStartVertex() 
             || lastCounterCWEdge->getStartVertex() == lastCWEdge->getEndVertex()  )
        {
            currPartEdge = lastCounterCWPartEdge->getNextPartialEdgeInLoop();
            do
            {
                multipleOuterLoop.add( currPartEdge );

                currPartEdge = currPartEdge->getNextPartialEdgeInLoop();
            } while ( currPartEdge != lastCWPartEdge );

            lastCounterCWPartEdge->setNextPartialEdgeInLoop( lastCWPartEdge );
            lastCWPartEdge->setPreviousPartialEdgeInLoop( lastCounterCWPartEdge );

            return rg_FALSE;
        }

        currPartEdge = lastCounterCWPartEdge->getNextPartialEdgeInLoop();        
        while ( currPartEdge != lastCWPartEdge )
        {
            VDEdge* currEdge = currPartEdge->getOriginalEdge();

            //  =====>  : lastCounterCWPartEdge or lastCWPartEdge
            //  ----->  : currPartEdge

            if ( lastCWEdge->getStartVertex() == currEdge->getStartVertex() )
            {
                //   -----loop--->
                //   <-----o=====>
                if ( lastCWPartEdge->isRightOrientationInLoop() )
                    currPartEdge->moveThisToPreviousPartialEdgeInLoopOf( rg_FALSE, lastCWPartEdge );
                //   <----loop----
                //   <-----o=====>
                else
                    currPartEdge->moveThisToPreviousPartialEdgeInLoopOf( rg_TRUE, lastCWPartEdge );

                lastCWPartEdge = currPartEdge;
                break;
            }
            else if ( lastCWEdge->getStartVertex() == currEdge->getEndVertex() )
            {
                //   -----loop--->
                //   ----->o=====>
                if ( lastCWPartEdge->isRightOrientationInLoop() )
                    currPartEdge->moveThisToPreviousPartialEdgeInLoopOf( rg_TRUE, lastCWPartEdge );
                //   <----loop----
                //   ----->o=====>
                else
                    currPartEdge->moveThisToPreviousPartialEdgeInLoopOf( rg_FALSE, lastCWPartEdge );

                lastCWPartEdge = currPartEdge;
                break;
            }
            else if ( lastCWEdge->getEndVertex() == currEdge->getStartVertex() )
            {
                //   <----loop----
                //   <-----o<=====
                if ( lastCWPartEdge->isRightOrientationInLoop() )
                    currPartEdge->moveThisToPreviousPartialEdgeInLoopOf( rg_TRUE, lastCWPartEdge );
                //   ----loop---->
                //   <-----o<=====
                else
                    currPartEdge->moveThisToPreviousPartialEdgeInLoopOf( rg_FALSE, lastCWPartEdge );

                lastCWPartEdge = currPartEdge;
                break;
            }
            else if ( lastCWEdge->getEndVertex() == currEdge->getEndVertex() )
            {
                //   <----loop----
                //   ----->o<=====
                if ( lastCWPartEdge->isRightOrientationInLoop() )
                    currPartEdge->moveThisToPreviousPartialEdgeInLoopOf( rg_FALSE, lastCWPartEdge );
                //   ----loop---->
                //   ----->o<=====
                else
                    currPartEdge->moveThisToPreviousPartialEdgeInLoopOf( rg_TRUE, lastCWPartEdge );

                lastCWPartEdge = currPartEdge;
                break;
            }
            else if ( lastCounterCWEdge->getEndVertex() == currEdge->getStartVertex() )
            {
                //   -----loop--->
                //   =====>o----->
                if ( lastCounterCWPartEdge->isRightOrientationInLoop() )
                    currPartEdge->moveThisToNextPartialEdgeInLoopOf( rg_TRUE, lastCounterCWPartEdge );
                //   <----loop----
                //   =====>o----->
                else
                    currPartEdge->moveThisToNextPartialEdgeInLoopOf( rg_FALSE, lastCounterCWPartEdge );

                lastCounterCWPartEdge = currPartEdge;
                break;
            }
            else if ( lastCounterCWEdge->getEndVertex() == currEdge->getEndVertex() )
            {
                //   -----loop--->
                //   =====>o<-----
                if ( lastCounterCWPartEdge->isRightOrientationInLoop() )
                    currPartEdge->moveThisToNextPartialEdgeInLoopOf( rg_FALSE, lastCounterCWPartEdge );
                //   <----loop----
                //   =====>o<-----
                else
                    currPartEdge->moveThisToNextPartialEdgeInLoopOf( rg_TRUE, lastCounterCWPartEdge );

                lastCounterCWPartEdge = currPartEdge;
                break;
            }
            else if ( lastCounterCWEdge->getStartVertex() == currEdge->getStartVertex() )
            {
                //   <---loop----
                //   <=====o----->
                if ( lastCounterCWPartEdge->isRightOrientationInLoop() )
                    currPartEdge->moveThisToNextPartialEdgeInLoopOf( rg_FALSE, lastCounterCWPartEdge );
                //   ----loop---->
                //   <=====o----->
                else
                    currPartEdge->moveThisToNextPartialEdgeInLoopOf( rg_TRUE, lastCounterCWPartEdge );

                lastCounterCWPartEdge = currPartEdge;
                break;
            }
            else if ( lastCounterCWEdge->getStartVertex() == currEdge->getEndVertex() )
            {
                //   <---loop----
                //   <=====o<----
                if ( lastCounterCWPartEdge->isRightOrientationInLoop() )
                    currPartEdge->moveThisToNextPartialEdgeInLoopOf( rg_TRUE, lastCounterCWPartEdge );
                //   ----loop---->
                //   <=====o<-----
                else
                    currPartEdge->moveThisToNextPartialEdgeInLoopOf( rg_FALSE, lastCounterCWPartEdge );

                lastCounterCWPartEdge = currPartEdge;
                break;
            }
            else
            {
                currPartEdge = currPartEdge->getNextPartialEdgeInLoop();
            }
        }
    }

    return rg_TRUE;
/*
        VDPartialEdge* prevPartEdge = currPartEdge->getPreviousPartialEdgeInLoop();

        VDEdge* currEdge = currPartEdge->getOriginalEdge();
        VDEdge* prevEdge = currPartEdge->getPreviousPartialEdgeInLoop()->getOriginalEdge();

        if ( currPartEdge->isRightOrientationInLoop() == rg_TRUE )
        {
            if ( currEdge->getStartVertex() == prevEdge->getStartVertex() )
            {
                prevPartEdge->isRightOrientationInLoop( rg_FALSE );
            }
            else if ( currEdge->getStartVertex() == prevEdge->getEndVertex() )
            {
                prevPartEdge->isRightOrientationInLoop( rg_TRUE );
            }
            else
            {
                isLiningUpCompleted = rg_FALSE;
                break;
            }
        }
        else if ( currPartEdge->isRightOrientationInLoop() == rg_FALSE )
        {
            if ( currEdge->getEndVertex() == prevEdge->getStartVertex() )
            {
                prevPartEdge->isRightOrientationInLoop( rg_FALSE );
            }
            else if ( currEdge->getEndVertex() == prevEdge->getEndVertex() )
            {
                prevPartEdge->isRightOrientationInLoop( rg_TRUE );
            }
            else
            {
                isLiningUpCompleted = rg_FALSE;
                break;
            }
        }
        else {}

        currPartEdge = prevPartEdge;
    }
    
    if ( isLiningUpCompleted )
        return rg_TRUE;
    else
        return rg_FALSE;
*/

}

rg_FLAG VDLoop::reorderOfPartialEdgesInLoop( rg_dList<VDPartialEdge*>& multipleOuterLoop )
{
    VDEdge* currEdge = m_partEdge->getOriginalEdge();

    VDVertex* vertexOfCurrEdge = rg_NULL;
    VDVertex* vertexOfNextEdge = rg_NULL;
    VDVertex* oppositeVertexOfNextEdge = rg_NULL;

    VDPartialEdge* startPrEdge        = m_partEdge;
    VDPartialEdge* endConnectedPrEdge = rg_NULL;



    rg_FLAG isStartEdgeConnected = rg_FALSE;

    VDPartialEdge* currPrEdge  = startPrEdge;
    VDPartialEdge* nextPrEdge  = rg_NULL;
    VDPartialEdge* prevPrEdge  = rg_NULL;


    ////////////////////////////////////////////////////////////////
    //
    //  find the END EDGE to be connected.
    //
    if ( startPrEdge->isRightOrientationInLoop() == rg_TRUE )
        vertexOfCurrEdge = currPrEdge->getOriginalEdge()->getEndVertex();
    else
        vertexOfCurrEdge = currPrEdge->getOriginalEdge()->getStartVertex();

    do
    {
        nextPrEdge  = currPrEdge->getNextPartialEdgeInLoop();
    

        if ( nextPrEdge->isRightOrientationInLoop() == rg_TRUE ) {
            vertexOfNextEdge         = nextPrEdge->getOriginalEdge()->getStartVertex();
            oppositeVertexOfNextEdge = nextPrEdge->getOriginalEdge()->getEndVertex();
        }
        else  {
            vertexOfNextEdge         = nextPrEdge->getOriginalEdge()->getEndVertex();
            oppositeVertexOfNextEdge = nextPrEdge->getOriginalEdge()->getStartVertex();
        }


        if ( vertexOfCurrEdge == vertexOfNextEdge )
        {
            currPrEdge         = nextPrEdge;
            endConnectedPrEdge = nextPrEdge;

            vertexOfCurrEdge   = oppositeVertexOfNextEdge;

            isStartEdgeConnected = rg_TRUE;
        }
        else
        {
            endConnectedPrEdge = currPrEdge;
            break;
        }

    } while ( startPrEdge != currPrEdge );




    ////////////////////////////////////////////////////////////////
    //
    //  find the START EDGE to be connected.
    //
    if (    endConnectedPrEdge != startPrEdge 
         && endConnectedPrEdge != startPrEdge->getPreviousPartialEdgeInLoop() )
    {
        if ( startPrEdge->isRightOrientationInLoop() == rg_TRUE )
            vertexOfCurrEdge = currPrEdge->getOriginalEdge()->getStartVertex();
        else
            vertexOfCurrEdge = currPrEdge->getOriginalEdge()->getEndVertex();

        VDVertex* vertexOfPrevEdge = rg_NULL;
        VDVertex* oppositeVertexOfPrevEdge = rg_NULL;
        currPrEdge  = startPrEdge;

        do
        {
            prevPrEdge = currPrEdge->getPreviousPartialEdgeInLoop();
    

            if ( prevPrEdge->isRightOrientationInLoop() == rg_TRUE )  {
                vertexOfPrevEdge         = prevPrEdge->getOriginalEdge()->getEndVertex();
                oppositeVertexOfPrevEdge = prevPrEdge->getOriginalEdge()->getStartVertex();
            }
            else  {
                vertexOfPrevEdge         = prevPrEdge->getOriginalEdge()->getStartVertex();
                oppositeVertexOfPrevEdge = prevPrEdge->getOriginalEdge()->getEndVertex();
            }

            if ( vertexOfCurrEdge == vertexOfPrevEdge )
            {
                currPrEdge = prevPrEdge;
                m_partEdge = prevPrEdge;

                vertexOfCurrEdge = oppositeVertexOfPrevEdge;
            }
            else
            {
                break;
            }

        } while ( endConnectedPrEdge != currPrEdge );
    }
    else
    {
    
    }

    //  ==>  The partial edges from startPrEdge to endConnectedPrEdge are connected
    //
    /////////////////////////////////////////////////////////////////////////////////

    startPrEdge = m_partEdge;
    VDVertex* startVertexOfLoop = rg_NULL;
    if ( startPrEdge->isRightOrientationInLoop() == rg_TRUE )
        startVertexOfLoop = startPrEdge->getOriginalEdge()->getStartVertex();
    else
        startVertexOfLoop = startPrEdge->getOriginalEdge()->getEndVertex();

    while ( rg_TRUE )
    {
        //////////////////////////////////////////////////////////
        //  There are two types of faces: bounded face and unbounded face.
        //  The bounded face consists of one outer loop 
        //  with the connected partial edges of bounded edges.
        //  The unbounded face consists of only one unbounded edges with two extream vertices on infinity.
        //  It may define the degenerate case.
        if ( startPrEdge == endConnectedPrEdge )
        {
            if ( isStartEdgeConnected )
                return rg_TRUE;
            else
            {
                currPrEdge = startPrEdge->getNextPartialEdgeInLoop();
                do  
                {
                    multipleOuterLoop.add( currPrEdge );
                    m_numOfBoundingEdges--;
                    currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
                } while ( startPrEdge != currPrEdge );        
        
                startPrEdge->setPreviousPartialEdgeInLoop( endConnectedPrEdge );
                endConnectedPrEdge->setNextPartialEdgeInLoop( startPrEdge );

                return rg_FALSE;
            }
        }

        //////////////////////////////////////////////////////////
        //  The face is unbounded and consists of one outer loop.
        if ( endConnectedPrEdge == startPrEdge->getPreviousPartialEdgeInLoop() )
            return rg_TRUE;


        //////////////////////////////////////////////////////////
        //  The face consists of more than one outer loops.
        VDVertex* endVertexOfLoop = rg_NULL;
        if ( endConnectedPrEdge->isRightOrientationInLoop() == rg_TRUE )
            endVertexOfLoop = endConnectedPrEdge->getOriginalEdge()->getEndVertex();
        else
            endVertexOfLoop = endConnectedPrEdge->getOriginalEdge()->getStartVertex();

        if (    startVertexOfLoop == endVertexOfLoop 
             || (    startPrEdge->getOriginalEdge()->isBoundedEdge() == rg_FALSE 
                  && endConnectedPrEdge->getOriginalEdge()->isBoundedEdge() == rg_FALSE ) 
           )
        {
            currPrEdge = endConnectedPrEdge->getNextPartialEdgeInLoop();
            do  
            {
                multipleOuterLoop.add( currPrEdge );
                m_numOfBoundingEdges--;

                currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
            } while ( startPrEdge != currPrEdge );        
        
            startPrEdge->setPreviousPartialEdgeInLoop( endConnectedPrEdge );
            endConnectedPrEdge->setNextPartialEdgeInLoop( startPrEdge );

            return rg_FALSE;
        }



        currPrEdge = endConnectedPrEdge;
        nextPrEdge = endConnectedPrEdge->getNextPartialEdgeInLoop();

        if ( currPrEdge->isRightOrientationInLoop() == rg_TRUE )
            vertexOfCurrEdge = currPrEdge->getOriginalEdge()->getEndVertex();
        else
            vertexOfCurrEdge = currPrEdge->getOriginalEdge()->getStartVertex();

        if ( nextPrEdge->isRightOrientationInLoop() == rg_TRUE )
            vertexOfNextEdge = nextPrEdge->getOriginalEdge()->getStartVertex();
        else
            vertexOfNextEdge = nextPrEdge->getOriginalEdge()->getEndVertex();


        VDPartialEdge* prEdgeDisconnected = startPrEdge->getPreviousPartialEdgeInLoop();

        rg_FLAG connectionStatus = 0;
        //  connectionStatus = 0 : disconnected
        //  connectionStatus = 1 : connect with endConnectedPrEdge
        //  connectionStatus = 2 : connect with next prEdge of endConnectedPrEdge in loop
        //  connectionStatus = 3 : connect with endConnectedPrEdge and next prEdge of endConnectedPrEdge in loop

        do
        {
            VDVertex* startVertex = prEdgeDisconnected->getOriginalEdge()->getStartVertex();
            VDVertex* endVertex   = prEdgeDisconnected->getOriginalEdge()->getEndVertex();

            //  _____ : endConnectedPrEdge
            //  ===== : endConnectedPrEdge->getNextPartialEdgeInLoop()
            //  ----- : prEdgeDisconnected
            //  *     : vertexOfCurrEdge(vertex of endConnectedPrEdge)
            //  o     : vertexOfNextEdge

            //   -------loop------>
            //   _____*----->o=====
            if ( prEdgeDisconnected->isRightOrientationInLoop() == rg_TRUE )
            {
                if ( vertexOfCurrEdge == startVertex && vertexOfNextEdge == endVertex )
                    connectionStatus = 3;
                else if ( vertexOfCurrEdge == startVertex )
                    connectionStatus = 1;
                else if ( vertexOfNextEdge == endVertex )
                    connectionStatus = 2;
                else 
                    connectionStatus = 0;
            }
            //   -------loop------>
            //   _____*<-----o=====
            else
            {
                if ( vertexOfCurrEdge == endVertex && vertexOfNextEdge == startVertex )
                    connectionStatus = 3;
                else if ( vertexOfCurrEdge == endVertex )
                    connectionStatus = 1;
                else if ( vertexOfNextEdge == startVertex )
                    connectionStatus = 2;
                else 
                    connectionStatus = 0;
            }

            if ( connectionStatus == 0 )
                prEdgeDisconnected = prEdgeDisconnected->getPreviousPartialEdgeInLoop();
            else if ( connectionStatus == 1 )
            {
                prEdgeDisconnected->moveToNextOf(endConnectedPrEdge);
                endConnectedPrEdge = prEdgeDisconnected;
                break;
            }
            else if ( connectionStatus == 2 )
            {
                prEdgeDisconnected->moveToNextOf(endConnectedPrEdge);
                prEdgeDisconnected = startPrEdge->getPreviousPartialEdgeInLoop();
            }
            else // if ( connectionStatus == 3 )
            {
                prEdgeDisconnected->moveToNextOf(endConnectedPrEdge);
                endConnectedPrEdge = prEdgeDisconnected->getNextPartialEdgeInLoop();

                currPrEdge  = endConnectedPrEdge;
                do
                {
                    nextPrEdge  = currPrEdge->getNextPartialEdgeInLoop();
    
                    if ( currPrEdge->isRightOrientationInLoop() == rg_TRUE )
                        vertexOfCurrEdge = currPrEdge->getOriginalEdge()->getEndVertex();
                    else
                        vertexOfCurrEdge = currPrEdge->getOriginalEdge()->getStartVertex();

                    if ( nextPrEdge->isRightOrientationInLoop() == rg_TRUE )
                        vertexOfNextEdge = nextPrEdge->getOriginalEdge()->getStartVertex();
                    else
                        vertexOfNextEdge = nextPrEdge->getOriginalEdge()->getEndVertex();


                    if ( vertexOfCurrEdge == vertexOfNextEdge )
                    {
                        currPrEdge         = nextPrEdge;
                        endConnectedPrEdge = nextPrEdge;
                    }
                    else
                    {
                        endConnectedPrEdge = currPrEdge;
                        break;
                    }
                } while ( startPrEdge != currPrEdge );
                break;
            }
        } while ( prEdgeDisconnected != endConnectedPrEdge );
    }
    return rg_TRUE;
}


rg_FLAG VDLoop::constructLoopCycle( rg_dList<VDPartialEdge*>& multipleOuterLoop )
{
    VDPartialEdge* startPrEdgeConnected = m_partEdge;
    VDPartialEdge* endPrEdgeConnected   = m_partEdge;

    VDPartialEdge* currPrEdge = rg_NULL;
    VDEdge*        currEdge   = rg_NULL;

    VDVertex* connectingVertex = rg_NULL;

    /////////////////////////////////////////////////////////////////////////
    //
    //  find END partial edge to be connected
    //
    if ( endPrEdgeConnected->isRightOrientationInLoop() )
        connectingVertex = endPrEdgeConnected->getOriginalEdge()->getEndVertex();
    else
        connectingVertex = endPrEdgeConnected->getOriginalEdge()->getStartVertex();

    do  {   
        if ( connectingVertex->isOnInfinity() )
            break;

        rg_FLAG isFoundConnectedPrEdge = rg_FALSE;

        currPrEdge = endPrEdgeConnected->getNextPartialEdgeInLoop();
        do  {
            currEdge = currPrEdge->getOriginalEdge();

            if ( connectingVertex == currEdge->getStartVertex() )  {
                currPrEdge->moveToNextOf(endPrEdgeConnected);
                
                connectingVertex   = currEdge->getEndVertex();
                
                isFoundConnectedPrEdge = rg_TRUE;
                break;
            }
            else if ( connectingVertex == currEdge->getEndVertex() )  {
                currPrEdge->moveToNextOf(endPrEdgeConnected);

                connectingVertex   = currEdge->getStartVertex();

                isFoundConnectedPrEdge = rg_TRUE;
                break;
            }
            else  {
                currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
            }

        } while ( currPrEdge != m_partEdge );

        if ( isFoundConnectedPrEdge )  
            endPrEdgeConnected = endPrEdgeConnected->getNextPartialEdgeInLoop();
        else
            break;

    } while ( endPrEdgeConnected != m_partEdge );


    /////////////////////////////////////////////////////////////////////
    //
    //  1. the BOUNDED face with a single outer loop
    //
    if ( m_partEdge == endPrEdgeConnected && !connectingVertex->isOnInfinity() )
        return rg_TRUE;



    /////////////////////////////////////////////////////////////////////////
    //
    //  find START partial edge to be connected
    //
    if ( startPrEdgeConnected->isRightOrientationInLoop() )
        connectingVertex = startPrEdgeConnected->getOriginalEdge()->getStartVertex();
    else
        connectingVertex = startPrEdgeConnected->getOriginalEdge()->getEndVertex();

    do  {   
        if ( connectingVertex->isOnInfinity() )
            break;

        rg_FLAG isFoundConnectedPrEdge = rg_FALSE;

        currPrEdge = startPrEdgeConnected->getPreviousPartialEdgeInLoop();
        while ( currPrEdge != endPrEdgeConnected )  {
            currEdge = currPrEdge->getOriginalEdge();

            if ( connectingVertex == currEdge->getStartVertex() )  {
                currPrEdge->moveToPreviousOf(startPrEdgeConnected);
                
                connectingVertex   = currEdge->getEndVertex();
                
                isFoundConnectedPrEdge = rg_TRUE;
                break;
            }
            else if ( connectingVertex == currEdge->getEndVertex() )  {
                currPrEdge->moveToPreviousOf(startPrEdgeConnected);

                connectingVertex   = currEdge->getStartVertex();

                isFoundConnectedPrEdge = rg_TRUE;
                break;
            }
            else  {
                currPrEdge = currPrEdge->getPreviousPartialEdgeInLoop();
            }
        } 

        if ( isFoundConnectedPrEdge )  
            startPrEdgeConnected = startPrEdgeConnected->getPreviousPartialEdgeInLoop();
        else
            break;

    } while ( startPrEdgeConnected != endPrEdgeConnected );

    /////////////////////////////////////////////////////////////////////
    //
    //  2. the UNBOUNDED face with a single outer loop
    //
    if ( endPrEdgeConnected->getNextPartialEdgeInLoop() == startPrEdgeConnected )
        return rg_TRUE;


    /////////////////////////////////////////////////////////////////////
    //
    //  3. the face with multiple outer loops
    //
    currPrEdge = endPrEdgeConnected->getNextPartialEdgeInLoop();
    do  
    {
        multipleOuterLoop.add( currPrEdge );
        m_numOfBoundingEdges--;
        currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
    } while ( currPrEdge != startPrEdgeConnected );        

    startPrEdgeConnected->setPreviousPartialEdgeInLoop( endPrEdgeConnected );
    endPrEdgeConnected->setNextPartialEdgeInLoop( startPrEdgeConnected );

    return rg_FALSE;
}

///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
VDLoop& VDLoop::operator =(const VDLoop& aLoop)
{
    if ( this == &aLoop )
        return *this;

    m_ID          = aLoop.m_ID;

    m_face        = aLoop.m_face;
    m_partEdge    = aLoop.m_partEdge;
    m_isOuterLoop = aLoop.m_isOuterLoop;

    m_numOfBoundingEdges = aLoop.m_numOfBoundingEdges;

    return *this;
}


///////////////////////////////////////////////////////////////////////////////
//
//  topological operators..
void VDLoop::collectBoundingPartialEdges( rg_dList<VDPartialEdge*>& boundingPartialEdges)
{
    VDPartialEdge* startPrEdge = m_partEdge;
    VDPartialEdge* currPrEdge  = startPrEdge;

    do {
        boundingPartialEdges.add( currPrEdge );

        currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
    } while ( currPrEdge != startPrEdge );
}



void VDLoop::collectBoundingPartialEdgesInCCW( rg_dList<VDPartialEdge*>& boundingPartialEdges)
{
    VDPartialEdge* startPrEdge = m_partEdge;
    VDPartialEdge* currPrEdge  = startPrEdge;

    do {
        boundingPartialEdges.add( currPrEdge );

        currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
    } while ( currPrEdge != startPrEdge );
}



void VDLoop::collectBoundingPartialEdgesInCW( rg_dList<VDPartialEdge*>& boundingPartialEdges)
{
    VDPartialEdge* startPrEdge = m_partEdge;
    VDPartialEdge* currPrEdge  = startPrEdge;

    do {
        boundingPartialEdges.add( currPrEdge );

        currPrEdge = currPrEdge->getPreviousPartialEdgeInLoop();
    } while ( currPrEdge != startPrEdge );
}



rg_BOOL VDLoop::searchBoundingEdges( rg_dList<VDEdge*>& boundingEdges ) const
{
    rg_BOOL isSearchCompleted = rg_TRUE;

    VDPartialEdge* startPrEdge = m_partEdge;
    VDPartialEdge* currPrEdge  = startPrEdge;

    do {
        boundingEdges.add( currPrEdge->getOriginalEdge() );

        currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
        if ( currPrEdge == rg_NULL ) {
            isSearchCompleted = rg_FALSE;
            break;
        }
    } while ( currPrEdge != startPrEdge );

    return isSearchCompleted;
}

