//#include "StdAfx.h"
#include "ComponentBSFace.h"
using namespace V::GeometryTier;


ComponentBSFace::ComponentBSFace(void)
    :m_ID(-1)
    , m_betaValue(0.0)
{
    m_isConsiderable = true;
}



ComponentBSFace::ComponentBSFace( const int& ID )
    :m_ID(ID)
    , m_betaValue( 0.0 )
{
    m_isConsiderable = true;
}



ComponentBSFace::ComponentBSFace( const ComponentBSFace& component )
{
    m_ID        = component.m_ID;
    m_betaFaces = component.m_betaFaces;
    m_isConsiderable = component.m_isConsiderable;
    m_betaValue = component.m_betaValue;
}



ComponentBSFace::~ComponentBSFace(void)
{
}



void ComponentBSFace::addBetaFace( FaceBU2D* betaFace )
{
    //betaFace->setID(m_ID);
    m_betaFaces.add(betaFace);
}



LoopOfComponent* ComponentBSFace::setOuterLoopAsLongestLoop()
{
    // set outer loop via length of loop

    if(m_loops.getSize() == 1)
    {
        m_outerLoop = m_loops.getFirstpEntity();
        m_outerLoop->isOuterLoop(true);

        return m_outerLoop;
    }


    LoopOfComponent* outerLoop = NULL;
    double maxLength = -1.0;
    double currLength = -1.0;

    m_loops.reset4Loop();
    while(m_loops.setNext4Loop())
    {
        LoopOfComponent* currLoop = m_loops.getpEntity();

        currLoop->isOuterLoop(false);
        currLength = currLoop->computeLength();

        if( currLength > maxLength )
        {
            m_outerLoop = currLoop;
            maxLength = currLength;
        }
    }

    m_outerLoop->isOuterLoop(true);
    return m_outerLoop;
}



void ComponentBSFace::computeArea()
{
    m_area = 0.0;
    double currArea = 0.0;

    int numBetaFaces = m_betaFaces.getSize();
    m_area = (double)numBetaFaces / 2.0;

    //m_betaFaces.reset4Loop();
    //while(m_betaFaces.setNext4Loop()) {
    //    FaceBU2D* currFace = m_betaFaces.getEntity();
    //    currArea = currFace->computeSignedArea();

    //    if(currArea < 0.0)
    //    {
    //        currArea = currArea * -1.0;
    //    }

    //    m_area += currArea;
    //}
}



void ComponentBSFace::computePerimeter( const double& betaValue )
{
    double lengthOfDiagonal = sqrt( 2.0 );
    m_perimeter = 0.0;

    rg_dList<EdgeBU2D*> boundaryEdges;

    m_betaFaces.reset4Loop();
    while ( m_betaFaces.setNext4Loop() ) {
        FaceBU2D* currFace = m_betaFaces.getEntity();
        rg_dList<EdgeBU2D*> currBoundingEdges;
        currFace->getBoundingEdges( currBoundingEdges );

        currBoundingEdges.reset4Loop();
        while ( currBoundingEdges.setNext4Loop() ) {
            EdgeBU2D* currEdge = currBoundingEdges.getEntity();
            if ( currEdge->getBoundingState( betaValue ) == REGULAR_SIMPLEX && !currEdge->isVisited() )
            {
                currEdge->isVisited( true );

                if ( currEdge->getStartVertex()->getCoord().getX() == currEdge->getEndVertex()->getCoord().getX()
                    || currEdge->getStartVertex()->getCoord().getY() == currEdge->getEndVertex()->getCoord().getY() ) {
                    m_perimeter += 1.0;
                }
                else {
                    m_perimeter += lengthOfDiagonal;
                }
                boundaryEdges.add( currEdge );
            }
        }
    }

    boundaryEdges.reset4Loop();
    while ( boundaryEdges.setNext4Loop() ) {
        EdgeBU2D* currEdge = boundaryEdges.getEntity();
        currEdge->isVisited( false );
    }
}



void ComponentBSFace::computeCentroid()
{
    m_betaFaces.reset4Loop();
    while(m_betaFaces.setNext4Loop()) {
        FaceBU2D* currFace = m_betaFaces.getEntity();

        rg_dList<VertexBU2D*> boundingVertices;
        currFace->getBoundingVertices(boundingVertices);

        boundingVertices.reset4Loop();
        while(boundingVertices.setNext4Loop()) {
            VertexBU2D* currVertex = boundingVertices.getEntity();

            currVertex->isVisited(false);
        }
    }

    double totalSumXCoord = 0.0;
    double totalSumYCoord = 0.0;
    int    NumVerticesOnComponent = 0;
    rg_dList<VertexBU2D*> visitedVertices;

    m_betaFaces.reset4Loop();
    while(m_betaFaces.setNext4Loop()) {
        FaceBU2D* currFace = m_betaFaces.getEntity();

        rg_dList<VertexBU2D*> boundingVertices;
        currFace->getBoundingVertices(boundingVertices);

        boundingVertices.reset4Loop();
        while(boundingVertices.setNext4Loop()) {
            VertexBU2D* currVertex = boundingVertices.getEntity();

            if(!currVertex->isVisited())
            {
                totalSumXCoord += currVertex->getCoord().getX();
                totalSumYCoord += currVertex->getCoord().getY();
                NumVerticesOnComponent++;
                currVertex->isVisited(true);
                visitedVertices.add(currVertex);
            }
        }
    }

    m_centroid.setX(totalSumXCoord/(double)NumVerticesOnComponent);
    m_centroid.setY(totalSumYCoord/(double)NumVerticesOnComponent);


    visitedVertices.reset4Loop();
    while(visitedVertices.setNext4Loop()) {
        VertexBU2D* currVertex = visitedVertices.getEntity();
        currVertex->isVisited(false);
    }
}



void ComponentBSFace::computeCircularity_with_area_N_Perimeter()
{
    double radiusOfApproximatedCircle = sqrt( m_area / rg_PI );
    double boundaryLengthOfApproximatedCircle = 2.0 * rg_PI * radiusOfApproximatedCircle;

    m_circularity = boundaryLengthOfApproximatedCircle / m_perimeter;
}



void ComponentBSFace::computeCircularity(const double& betaValue)
{
    m_betaFaces.reset4Loop();
    while(m_betaFaces.setNext4Loop()) {
        FaceBU2D* currFace = m_betaFaces.getEntity();
        rg_dList<EdgeBU2D*> currBoundingEdges;
        currFace->getBoundingEdges(currBoundingEdges);

        currBoundingEdges.reset4Loop();
        while(currBoundingEdges.setNext4Loop()) {
            EdgeBU2D* currEdge = currBoundingEdges.getEntity();
            currEdge->isVisited(false);
        }
    }


    rg_dList<EdgeBU2D*> boundaryEdges;

    m_betaFaces.reset4Loop();
    while(m_betaFaces.setNext4Loop()) {
        FaceBU2D* currFace = m_betaFaces.getEntity();
        rg_dList<EdgeBU2D*> currBoundingEdges;
        currFace->getBoundingEdges(currBoundingEdges);

        currBoundingEdges.reset4Loop();
        while(currBoundingEdges.setNext4Loop()) {
            EdgeBU2D* currEdge = currBoundingEdges.getEntity();
            if( currEdge->getBoundingState(betaValue) == REGULAR_SIMPLEX && !currEdge->isVisited() )
            {
                currEdge->isVisited(true);
                boundaryEdges.add(currEdge);
            }
        }
    }


    double lengthOfComponentBoundary = 0.0;

    // reset the visited record
    boundaryEdges.reset4Loop();
    while(boundaryEdges.setNext4Loop())
    {
        EdgeBU2D* currEdge = boundaryEdges.getEntity();
        double currLength = currEdge->getStartVertex()->getCoord().distance( currEdge->getEndVertex()->getCoord() );
        lengthOfComponentBoundary += currLength;

        boundaryEdges.getEntity()->isVisited(false);
    }


    double radiusOfApproximatedCircle = sqrt(m_area / rg_PI);
    double boundaryLengthOfApproximatedCircle = 2.0 * rg_PI * radiusOfApproximatedCircle;



    m_circularity = boundaryLengthOfApproximatedCircle / lengthOfComponentBoundary;
}



void ComponentBSFace::computeLoops()
{
    rg_dList<EdgeBU2D*> boundaryEdges;
    m_betaFaces.reset4Loop();
    while ( m_betaFaces.setNext4Loop() )
    {
        FaceBU2D* currFace = m_betaFaces.getEntity();
        rg_dList<EdgeBU2D*> currBoundingEdges;
        currFace->getBoundingEdges( currBoundingEdges );

        currBoundingEdges.reset4Loop();
        while ( currBoundingEdges.setNext4Loop() )
        {
            EdgeBU2D* currEdge = currBoundingEdges.getEntity();
            if ( currEdge->getBoundingState( m_betaValue ) == REGULAR_SIMPLEX && !currEdge->isVisited() )
            {
                currEdge->isVisited( true );
                boundaryEdges.add( currEdge );
            }
        }
    }



    // reset the visited record
    boundaryEdges.reset4Loop();
    while ( boundaryEdges.setNext4Loop() )
    {
        boundaryEdges.getEntity()->isVisited( false );
    }



    // separate edges into connected edge group
    bool        isThereAnotherLoop = false;
    EdgeBU2D*   seedEdgeOfLoop = boundaryEdges.getHead()->getEntity();

    do
    {
        LoopOfComponent* currLoop = m_loops.add( LoopOfComponent() );
        currLoop->setComponent( this );

        isThereAnotherLoop = false;
        propagateLoopWithOrientation( seedEdgeOfLoop, currLoop );

        boundaryEdges.reset4Loop();
        while ( boundaryEdges.setNext4Loop() )
        {
            EdgeBU2D* currBoundaryEdge = boundaryEdges.getEntity();
            if ( !currBoundaryEdge->isVisited() )
            {
                isThereAnotherLoop = true;
                seedEdgeOfLoop = currBoundaryEdge;
                break;
            }
        }
    } while ( isThereAnotherLoop );



    // reset the visited record
    boundaryEdges.reset4Loop();
    while ( boundaryEdges.setNext4Loop() )
    {
        boundaryEdges.getEntity()->isVisited( false );
    }

    classifyLoopsAsAnOuter_N_innerLoops();
}



LoopOfComponent* ComponentBSFace::findLoopIncludingThisEdge( EdgeBU2D* edge )
{
    map<EdgeBU2D*, LoopOfComponent*>::iterator i_edgeToLoop;
    // Added by Joonghyun on November 18, 2020
    //unordered_map<EdgeBU2D*, LoopOfComponent*>::iterator i_edgeToLoop;
    i_edgeToLoop = m_mapEdgeToLoop.find(edge);
    LoopOfComponent* loopIncludingThisEdge;

    if( i_edgeToLoop == m_mapEdgeToLoop.end() )
    {
        loopIncludingThisEdge = NULL;
    }
    else
    {
        loopIncludingThisEdge = i_edgeToLoop->second;
    }

    return loopIncludingThisEdge;

}



void ComponentBSFace::insertEdgeLoopPair( EdgeBU2D* edge, LoopOfComponent* loop )
{
    m_mapEdgeToLoop[edge] = loop;
}



bool ComponentBSFace::isThereNeighborList( const int& ID, NeighborInfo*& neighbor )
{
    if(m_neighbors.getSize() == 0)
    {
        return false;
    }

    m_neighbors.reset4Loop();
    while(m_neighbors.setNext4Loop()) {
        NeighborInfo* currNeighbor = m_neighbors.getpEntity();

        if( currNeighbor->m_ID == ID )
        {
            neighbor = currNeighbor;
            return true;
        }
    }

    neighbor = NULL;
    return false;
}


void ComponentBSFace::classifyLoopsAsAnOuter_N_innerLoops()
{
    m_outerLoop = NULL;
    if ( m_loops.getSize() == 1 )
    {
        m_outerLoop = m_loops.getFirstpEntity();
        m_outerLoop->isOuterLoop( true );
    }
    else {
        m_loops.reset4Loop();
        while ( m_loops.setNext4Loop() )
        {
            LoopOfComponent* currLoop = m_loops.getpEntity();

            if ( m_outerLoop != NULL ) {
                currLoop->isOuterLoop( false );
                m_innerLoops.add( currLoop );
            }
            else {
                double signedArea = currLoop->computeSignedArea();
                if ( signedArea > 0.0 ) {
                    m_outerLoop = currLoop;
                    m_outerLoop->isOuterLoop( true );
                }
                else {
                    currLoop->isOuterLoop( false );
                    m_innerLoops.add( currLoop );
                }
            }
        }
    }
}



void ComponentBSFace::propagateLoopWithOrientation( const EdgeBU2D* const seedEdge, LoopOfComponent* const loop )
{

    // 이 방법은 loop의 바깥쪽에서 작대기를 반시계방향으로 돌리면서 찾는 방법
    // 시계방향으로 뒤지면 component 정보가 필요없고 
    // 함수가 단순해질 수 있다.
    ComponentBSFace* currComponent = loop->getComponent();


    VertexBU2D* anchorVertex = NULL;
    //if ( seedEdge->getLeftFace()->getBoundingState( m_betaValue ) == INTERIOR_SIMPLEX ) {
    //    anchorVertex = seedEdge->getStartVertex();
    //}
    //else {
    //    anchorVertex = seedEdge->getEndVertex();
    //}


    EdgeBU2D* currEdge = const_cast<EdgeBU2D*>( seedEdge );
    do
    {
        if ( currEdge->getBoundingState( m_betaValue ) == REGULAR_SIMPLEX ) {
            if ( currEdge->getLeftFace()->getBoundingState( m_betaValue ) == INTERIOR_SIMPLEX ) {
                anchorVertex = currEdge->getEndVertex();
                loop->addEdge( currEdge );
                currComponent->insertEdgeLoopPair( currEdge, loop );
                currEdge->isVisited( true );

                currEdge = currEdge->getLeftHand();
            }
            else {
                anchorVertex = currEdge->getStartVertex();
                loop->addEdge( currEdge );
                currComponent->insertEdgeLoopPair( currEdge, loop );
                currEdge->isVisited( true );

                currEdge = currEdge->getRightLeg();
            }
        }
        else { // currEdge is INTERIOR_SIMPLEX
            if ( anchorVertex == currEdge->getStartVertex() )
            {
                currEdge = currEdge->getRightLeg();
            }
            else
            {
                currEdge = currEdge->getLeftHand();
            }
        }
    } while ( currEdge != seedEdge );
}



// void ComponentBSFace::fillUp( const double& betaValue )
// {
//     // bring boundary edges
//     rg_dList<EdgeBU2D*> boundaryEdges;
//     m_betaFaces.reset4Loop();
//     while(m_betaFaces.setNext4Loop())
//     {
//         FaceBU2D* currFace = m_betaFaces.getEntity();
//         rg_dList<EdgeBU2D*> currBoundingEdges;
//         currFace->getBoundingEdges(currBoundingEdges);
//         
//         currBoundingEdges.reset4Loop();
//         while(currBoundingEdges.setNext4Loop())
//         {
//             EdgeBU2D* currEdge = currBoundingEdges.getEntity();
//             if(currEdge->getBoundingState(betaValue) == REGULAR_SIMPLEX && !currEdge->isVisited())
//             {
//                 currEdge->isVisited(TRUE);
//                 boundaryEdges.add(currEdge);
//             }
//         }
//     }
// 
//     // separate edges into connected edge group
// 
//     // 
// 
// 
//     // reset the visited record
//     boundaryEdges.reset4Loop();
//     while(boundaryEdges.setNext4Loop());
//     {
//          boundaryEdges.getEntity()->isVisited(FALSE);
//     }
// }


