#include "VoronoiDiagram2DP.h"

#include <vector>
using namespace std;

using namespace BULL2D::GeometryTier;

VoronoiDiagram2DP::VoronoiDiagram2DP(void)
{
}



VoronoiDiagram2DP::VoronoiDiagram2DP( list<rg_Point2D>& pointset )
{
    int currID = 1;
    for( list<rg_Point2D>::iterator it_point = pointset.begin() ; it_point != pointset.end() ; it_point++, currID++ ) {
        m_generators.push_back( Generator2DP( currID, &(*it_point) ) );
    }
}



VoronoiDiagram2DP::VoronoiDiagram2DP( const VoronoiDiagram2DP& VD2DP )
{

}



VoronoiDiagram2DP::~VoronoiDiagram2DP(void)
{
}



void VoronoiDiagram2DP::getVoronoiEdges( list<const VEdge2DP*>& VEdgesList ) const
{
    for( list<VEdge2DP>::const_iterator it_edge = m_VEdges.begin() ; it_edge != m_VEdges.end() ; it_edge++ ) {
        const VEdge2DP* currEdge = &(*it_edge);
        VEdgesList.push_back( currEdge );
    }
}



void VoronoiDiagram2DP::getVoronoiFaces( list<const VFace2DP*>& VFacesList ) const
{
    for( list<VFace2DP>::const_iterator it_face = m_VFaces.begin() ; it_face != m_VFaces.end() ; it_face++ ) {
        const VFace2DP* currFace = &(*it_face);
        VFacesList.push_back( currFace );
    }
}



void VoronoiDiagram2DP::getVoronoiVertices( list<const VVertex2DP*>& VVerticesList ) const
{
    for( list<VVertex2DP>::const_iterator it_vertex = m_VVertices.begin() ; it_vertex != m_VVertices.end() ; it_vertex++ ) {
        const VVertex2DP* currVertex = &(*it_vertex);
        VVerticesList.push_back( currVertex );
    }
}



void VoronoiDiagram2DP::getGenerators( list<const Generator2DP*>& generatorsList ) const
{
    for( list<Generator2DP>::const_iterator it_generator = m_generators.begin() ; it_generator != m_generators.end() ; it_generator++ ) {
        const Generator2DP* currGenerator = &(*it_generator);
        generatorsList.push_back( currGenerator );
    }
}



void VoronoiDiagram2DP::setGenerators( list<rg_Point2D>& pointset )
{
    m_generators.clear();

    int currID = 1;
    for( list<rg_Point2D>::iterator it_point = pointset.begin() ; it_point != pointset.end() ; it_point++, currID++ ) {
        m_generators.push_back( Generator2DP( currID, &(*it_point) ) );
    }
}



void VoronoiDiagram2DP::setGenerators( list<Generator2DP>& generators )
{
    m_generators.clear();
    m_generators = generators;
}



void VoronoiDiagram2DP::constructVoronoiDiagram( list<rg_Point2D>& pointset )
{
    setGenerators( pointset );
    constructVoronoiDiagram();
}



void VoronoiDiagram2DP::constructVoronoiDiagram() // input generator list as argument
{
    makePhantomGenerators();
    constructVoronoiDiagramForPhantomGenerators();

    for( list<Generator2DP>::iterator it_generator = m_generators.begin() ; it_generator != m_generators.end() ; it_generator++ ) {
        Generator2DP* currGenerator = &(*it_generator);
        updateVoronoiDiagramWithNewGenerator( currGenerator );
    }

    // 원래는 루프 안에 있어야 하는데 반복하는 시간이 걸려서 댕글링 되있는 것을 이피션시를 위해서 한번에 제거 
    removeAllExtraneousVVerticesAndVEdges();
}



VFace2DP* VoronoiDiagram2DP::findVFaceContainingQueryPoint( const rg_Point2D& pt ) // 첫번째 가까운게 나오면 그리로 넘어가도록
{
    VFace2DP*    VFaceContainingPoint = &( *m_VFaces.rbegin() );
    double      minDistancePt2GeneratorOnVFace = pt.distance( *VFaceContainingPoint->getGenerator()->getLocation() );

    bool        isThereCloserGeneratorToInputPoint = false;

    do {
        list<Generator2DP*>  neighbor;
        double              distancePt2CurrNeighbor = 0.0;
        Generator2DP*        currNeighbor = NULL;
        isThereCloserGeneratorToInputPoint = false;

        VFaceContainingPoint->getGenerator()->getNeighborGenerators( neighbor );

        for( list<Generator2DP*>::iterator it_neighbor = neighbor.begin() ; it_neighbor != neighbor.end() ; it_neighbor++ ) {
            currNeighbor            = *it_neighbor;
            distancePt2CurrNeighbor = pt.distance( *currNeighbor->getLocation() );

            if( distancePt2CurrNeighbor < minDistancePt2GeneratorOnVFace ) {
                VFaceContainingPoint                = currNeighbor->getAssignedVFace();
                minDistancePt2GeneratorOnVFace      = distancePt2CurrNeighbor;
                isThereCloserGeneratorToInputPoint  = true;
            }
        }
    } while ( isThereCloserGeneratorToInputPoint );

    return VFaceContainingPoint;
}



void VoronoiDiagram2DP::clear()
{
    m_VEdges.clear();
    m_VFaces.clear();
    m_VVertices.clear();

    m_phantomGenerators.clear();
}



VVertex2DP* VoronoiDiagram2DP::createVVertex( const VVertex2DP& vertex )
{
    m_VVertices.push_back( vertex );

    return &(*m_VVertices.rbegin());
}



VEdge2DP* VoronoiDiagram2DP::createVEdge( const VEdge2DP& edge )
{
    m_VEdges.push_back( edge );

    return &(*m_VEdges.rbegin());
}



VFace2DP* VoronoiDiagram2DP::createVFace( const VFace2DP& face )
{
    m_VFaces.push_back( face );

    return &(*m_VFaces.rbegin());
}



void VoronoiDiagram2DP::makePhantomGenerators()
{
    // The triangle of phantom generator is a equilateral triangle
    rg_Point2D minPt;
    rg_Point2D maxPt;

    list<Generator2DP>::iterator it_generator = m_generators.begin();
    Generator2DP currGenerator = *it_generator++;

    minPt = *currGenerator.getLocation();
    maxPt = minPt;

    for( ; it_generator != m_generators.end() ; it_generator++ ) {
        currGenerator = *it_generator;
        if( currGenerator.getLocation()->getX() < minPt.getX() ) {
            minPt.setX( currGenerator.getLocation()->getX() );
        }
        if( currGenerator.getLocation()->getY() < minPt.getY() ) {
            minPt.setY( currGenerator.getLocation()->getY() );
        }
        if( currGenerator.getLocation()->getX() > maxPt.getX() ) {
            maxPt.setX( currGenerator.getLocation()->getX() );
        }
        if( currGenerator.getLocation()->getY() > maxPt.getY() ) {
            maxPt.setY( currGenerator.getLocation()->getY() );
        }
    }

    double lengthOfX = maxPt.getX() - minPt.getX();
    double lengthOfY = maxPt.getY() - minPt.getY();
    double length = 0.0;

    if( lengthOfX > lengthOfY ) {
        length = lengthOfX;
    }
    else {
        length = lengthOfY;
    }

    rg_Point2D* coordinate = new rg_Point2D[3];

    ///////////////////////////////////////////////////////////
    //                                                       //
    // Coordinates of phantom generators                     //
    // Concepts and Applications of Voronoi Diagrams - p.244 //
    //                                                       //
    ///////////////////////////////////////////////////////////

    coordinate[0].setX(minPt.getX() + 0.5 * length);
    coordinate[0].setY(minPt.getY() + (3.0 * sqrt(2.0)/2.0 + 0.5) * length);

    coordinate[1].setX(minPt.getX() + (-3.0 * sqrt(6.0)/4.0 + 0.5) * length);
    coordinate[1].setY(minPt.getY() + (-3.0 * sqrt(2.0)/4.0 + 0.5) * length);

    coordinate[2].setX(minPt.getX() + (3.0 * sqrt(6.0)/4.0 + 0.5) * length);
    coordinate[2].setY(minPt.getY() + (-3.0 * sqrt(2.0)/4.0 + 0.5) * length);


    m_phantomGenerators.push_back(Generator2DP(0, &coordinate[0]));
    m_phantomGenerators.push_back(Generator2DP(0, &coordinate[1]));
    m_phantomGenerators.push_back(Generator2DP(0, &coordinate[2]));
}



void VoronoiDiagram2DP::constructVoronoiDiagramForPhantomGenerators()
{
    //////////////////////////////////////////////////////////////////////////////////
    //                                                                              //
    //                                   process                                    //
    //                                                                              //
    //  1. make elements of voronoi diagram of phantom generators                   //
    //  2. set topology between elements of voronoi diagram of phantom generators   //
    //  3. define coordinate of vertices on voronoi diagram of phantom generators   //
    //      3.1 define coordinate of vertex defined by three phantom generators     //
    //      3.2 define coordinate of infinite vertices                              //
    //                                                                              //
    //////////////////////////////////////////////////////////////////////////////////

    list<Generator2DP>::iterator it_phantoms = m_phantomGenerators.begin();
    Generator2DP* phantoms[3] = {&(*it_phantoms++), &(*it_phantoms++), &(*it_phantoms++)};


    //  1. make elements of voronoi diagram of phantom generators

    VVertex2DP* vertices[4]   = { createVVertex( VVertex2DP(0) ),        // vertex defined by three phantom generators
                                  createVVertex( VVertex2DP(1) ),        // infinite vertex on edge0
                                  createVVertex( VVertex2DP(2) ),        // infinite vertex on edge1
                                  createVVertex( VVertex2DP(3) ) };      // infinite vertex on edge2

    VEdge2DP*   edges[6]      = { createVEdge( VEdge2DP(0) ),            // edge made by phantom0 and phantom1     
                                  createVEdge( VEdge2DP(1) ),            // edge made by phantom1 and phantom2
                                  createVEdge( VEdge2DP(2) ),            // edge made by phantom2 and phantom0
                                  createVEdge( VEdge2DP(3) ),            // artificial boundary for unbounded face1
                                  createVEdge( VEdge2DP(4) ),            // artificial boundary for unbounded face2 
                                  createVEdge( VEdge2DP(5) ) };          // artificial boundary for unbounded face3 

    VFace2DP* faces[4]        = { createVFace( VFace2DP(0) ),            // infinite space defined by artificial boundary for voronoi diagram
                                  createVFace( VFace2DP(1) ),            // face containing phantom0 
                                  createVFace( VFace2DP(2) ),            // face containing phantom1
                                  createVFace( VFace2DP(3) ) };          // face containing phantom2



    //  2. set topology between elements of voronoi diagram of phantom generators

    phantoms[0]->setAssignedVFace(faces[1]);
    phantoms[1]->setAssignedVFace(faces[2]);
    phantoms[2]->setAssignedVFace(faces[3]);

    vertices[0]->setFirstVEdge(edges[0]);
    vertices[1]->setFirstVEdge(edges[0]);
    vertices[2]->setFirstVEdge(edges[1]);
    vertices[3]->setFirstVEdge(edges[2]);

    edges[0]->setTopology(vertices[0],vertices[1],faces[2],faces[1],edges[4],edges[3],edges[1],edges[2]);
    edges[1]->setTopology(vertices[0],vertices[2],faces[3],faces[2],edges[5],edges[4],edges[2],edges[0]);
    edges[2]->setTopology(vertices[0],vertices[3],faces[1],faces[3],edges[3],edges[5],edges[0],edges[1]);
    edges[3]->setTopology(vertices[3],vertices[1],faces[1],faces[0],edges[0],edges[4],edges[2],edges[5]);
    edges[4]->setTopology(vertices[1],vertices[2],faces[2],faces[0],edges[1],edges[5],edges[0],edges[3]);
    edges[5]->setTopology(vertices[2],vertices[3],faces[3],faces[0],edges[2],edges[3],edges[1],edges[4]);

    faces[0]->setFirstVEdge(edges[3]);
    faces[1]->setGenerator(phantoms[0]);
    faces[1]->setFirstVEdge(edges[0]);
    faces[2]->setGenerator(phantoms[1]);
    faces[2]->setFirstVEdge(edges[1]);
    faces[3]->setGenerator(phantoms[2]);
    faces[3]->setFirstVEdge(edges[2]);



    //  3. define coordinate of vertices on voronoi diagram of phantom generators
    //      3.1. define coordinate in voronoi diagram

    rg_Point2D locationOfVertex0 = computeCircumcenter(*(phantoms[0]->getLocation()), *(phantoms[1]->getLocation()), *(phantoms[2]->getLocation()));
    vertices[0]->setLocation(locationOfVertex0);


    //      3.2. define coordinate of infinite vertices

    rg_Point2D centerOfPoint0And1 = (*(phantoms[0]->getLocation()) + *(phantoms[1]->getLocation())) / 2.0;
    rg_Point2D centerOfPoint1And2 = (*(phantoms[1]->getLocation()) + *(phantoms[2]->getLocation())) / 2.0;
    rg_Point2D centerOfPoint2And0 = (*(phantoms[2]->getLocation()) + *(phantoms[0]->getLocation())) / 2.0;

    rg_Point2D vectorOfVertex1 = centerOfPoint0And1 - vertices[0]->getLocation();
    rg_Point2D vectorOfVertex2 = centerOfPoint1And2 - vertices[0]->getLocation();
    rg_Point2D vectorOfVertex3 = centerOfPoint2And0 - vertices[0]->getLocation();

    rg_Point2D locationOfVertex1 = centerOfPoint0And1 + 5.0 * vectorOfVertex1;
    rg_Point2D locationOfVertex2 = centerOfPoint1And2 + 5.0 * vectorOfVertex2;
    rg_Point2D locationOfVertex3 = centerOfPoint2And0 + 5.0 * vectorOfVertex3;

    vertices[1]->setLocation(locationOfVertex1);
    vertices[2]->setLocation(locationOfVertex2);
    vertices[3]->setLocation(locationOfVertex3); 
}



void VoronoiDiagram2DP::updateVoronoiDiagramWithNewGenerator( Generator2DP* const newGenerator )
{
    Generator2DP* closestGenerator = findClosestGeneratorToNewGenerator( newGenerator );

    // closetGenerator 와 newGenerator 의 coordinate이 같은지 확인
    // 같은 점이면 아예 입력을 받지 않도록
    if( (closestGenerator->getLocation()->getX() == newGenerator->getLocation()->getX()) && 
        (closestGenerator->getLocation()->getY() == newGenerator->getLocation()->getY()) ) {
            return;
    }

    list<VEdge2DP*>   intersectingVEdges;
    findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop( newGenerator, closestGenerator, intersectingVEdges );

    list<VVertex2DP*> newVVertices;
    list<VEdge2DP*>   newVEdges;
    makeVEdgeLoopForNewGeneratorAndConnectToCurrVD( newGenerator, intersectingVEdges, newVVertices, newVEdges );

    connectCurrVDToNewVEdgeLoop( newVEdges );

    computeCoordOfNewVVertices( newVVertices );
}



Generator2DP* VoronoiDiagram2DP::findClosestGeneratorToNewGenerator( Generator2DP* const newGenerator )
{
    VFace2DP* closestVFace = findVFaceContainingQueryPoint( *newGenerator->getLocation() );
    return closestVFace->getGenerator();
}



void VoronoiDiagram2DP::removeAllExtraneousVVerticesAndVEdges()
{
    // while 문으로...
    for( list<VVertex2DP>::iterator it_vertex = m_VVertices.begin() ; it_vertex != m_VVertices.end() ; ) { 
        VVertex2DP currVertex = *it_vertex;
        if( currVertex.getStatus() == RED_V_P ) {
            it_vertex = m_VVertices.erase( it_vertex );
        }
        else {
            it_vertex++;
        }
    }

    for( list<VEdge2DP>::iterator it_edge = m_VEdges.begin() ; it_edge != m_VEdges.end() ;) {
        VEdge2DP currEdge = *it_edge;
        if( currEdge.getStatus() == RED_E_P ) {
            it_edge = m_VEdges.erase( it_edge );
        }
        else {
            it_edge++;
        }
    }
}



void VoronoiDiagram2DP::findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop( Generator2DP* const newGenerator, Generator2DP* const closestGenerator, list<VEdge2DP*>& intersectingVEdges )
{
    list<VVertex2DP*> redVVertices;
    VVertex2DP*  firstRedVertex = findFirstRedVertex( newGenerator, closestGenerator );

    firstRedVertex->setStatus( RED_V_P );
    redVVertices.push_back( firstRedVertex );

    list<VVertex2DP*> blueVVertices;
    list<VVertex2DP*>::iterator it_redVVertex = redVVertices.begin();
    for( ; it_redVVertex != redVVertices.end() ; it_redVVertex++ ) {
        VVertex2DP* currRedVVertex = *it_redVVertex;
        list<VVertex2DP*> neighborsOfCurrRedVVertex;
        currRedVVertex->getAdjacent3VVertices( neighborsOfCurrRedVVertex );

        for( list<VVertex2DP*>::iterator it_neighbor = neighborsOfCurrRedVVertex.begin() ; it_neighbor != neighborsOfCurrRedVVertex.end() ; it_neighbor++ ) {
            VVertex2DP* currNeighbor = *it_neighbor;

            if( currNeighbor->getStatus() != WHITE_V_P ) {
                continue;
            }

            double HValue = currNeighbor->computeHValue( newGenerator );

            if( HValue >= 0) {
                currNeighbor->setStatus( BLUE_V_P );
                blueVVertices.push_back( currNeighbor );
                continue;
            }

            if( hasCycleOccurred( currNeighbor ) ) {
                currNeighbor->setStatus( BLUE_V_P );
                blueVVertices.push_back( currNeighbor );
                continue;
            }

            if( hasOldRegionSplit( currNeighbor ) ) { // split 일어나는 횟수 count   
                currNeighbor->setStatus( BLUE_V_P );
                blueVVertices.push_back( currNeighbor );
                continue;
            }

            currNeighbor->setStatus( RED_V_P );
            redVVertices.push_back( currNeighbor );
        }
    }

    it_redVVertex = redVVertices.begin();
    for( ; it_redVVertex != redVVertices.end() ; it_redVVertex++ ) {
        VVertex2DP* currRedVVertex = *it_redVVertex;
        list<VEdge2DP*> incidentVEdges;
        currRedVVertex->getIncident3VEdges( incidentVEdges );

        for( list<VEdge2DP*>::iterator it_incidentVEdge = incidentVEdges.begin() ; it_incidentVEdge != incidentVEdges.end() ; it_incidentVEdge++ ) {
            VEdge2DP* currIncidentVEdge = *it_incidentVEdge;

            if( currIncidentVEdge->getStatus() != WHITE_E_P ) {
                continue;
            }

            if( currIncidentVEdge->getStartVertex()->getStatus() != currIncidentVEdge->getEndVertex()->getStatus() ) {
                intersectingVEdges.push_back( currIncidentVEdge );
            }
            else {
                currIncidentVEdge->setStatus( RED_E_P );
            }
        }
    }

    for( list<VVertex2DP*>::iterator it_blueVVertex = blueVVertices.begin() ; it_blueVVertex != blueVVertices.end() ; it_blueVVertex++ ) {
        VVertex2DP* currBlueVVertex = *it_blueVVertex;
        currBlueVVertex->setStatus( WHITE_V_P );
    }


}



void VoronoiDiagram2DP::makeVEdgeLoopForNewGeneratorAndConnectToCurrVD( Generator2DP* const newGenerator, list<VEdge2DP*>& intersectingVEdges, list<VVertex2DP*>& newVVertices, list<VEdge2DP*>& newVEdges )
{
    int currPos = 0;
    vector<VEdge2DP*> orderedIntersectingVEdges;
    vector<VFace2DP*> CCWVFacesOfOrderedVEdges;
    vector<VVertex2DP*> newVVerticesOnOrderedVEdges;

    // make intersecting VEdge list in CCW order and create newVVertices on the ordered VEdges
    orderedIntersectingVEdges.push_back( *intersectingVEdges.begin() );

    while( currPos < intersectingVEdges.size()) { 
        VEdge2DP* currEdge = orderedIntersectingVEdges[currPos];

        VFace2DP* CCWFace;

        if( currEdge->getStartVertex()->getStatus() == RED_V_P ) {
            CCWFace = currEdge->getLeftFace();
        }
        else {
            CCWFace = currEdge->getRightFace();
        }

        CCWVFacesOfOrderedVEdges.push_back( CCWFace );
        newVVerticesOnOrderedVEdges.push_back( createVVertex( VVertex2DP(m_VVertices.rbegin()->getID()+1) ) );

        newVVertices.push_back( newVVerticesOnOrderedVEdges[currPos] );


        list<VEdge2DP*>::iterator i_edge = intersectingVEdges.begin();
        for(; i_edge != intersectingVEdges.end(); i_edge++) {
            VEdge2DP* currIntersectingEdge = *i_edge;

            if(currIntersectingEdge == currEdge) {
                continue;
            }

            if(currIntersectingEdge->getRightFace() == CCWFace || currIntersectingEdge->getLeftFace() == CCWFace) {
                orderedIntersectingVEdges.push_back( currIntersectingEdge );
                currPos++;
                break;
            }            
        }        
    }
    newVVerticesOnOrderedVEdges.push_back(  newVVerticesOnOrderedVEdges[0] );


    // create new VEdges and connect new VEdges to new VVertices
    int i = 0;
    for ( i = 0 ; i < newVVerticesOnOrderedVEdges.size()-1 ; i++ ) {
        newVVerticesOnOrderedVEdges[i]->setFirstVEdge( orderedIntersectingVEdges[i] );

        VEdge2DP* newEdge = createVEdge( VEdge2DP(m_VEdges.rbegin()->getID() + 1) );
        newEdge->setStartVertex( newVVerticesOnOrderedVEdges[i] );
        newEdge->setEndVertex(   newVVerticesOnOrderedVEdges[i+1] );

        newVEdges.push_back( newEdge );    
    }

    vector<VEdge2DP*> orderedNewVEdges;
    i = 0;
    orderedNewVEdges.push_back( *newVEdges.rbegin() );
    list<VEdge2DP*>::iterator i_edge;
    for ( i_edge= newVEdges.begin(); i_edge != newVEdges.end(); i_edge++, i++ ) {
        orderedNewVEdges.push_back( *i_edge );
    }
    orderedNewVEdges.push_back( *newVEdges.begin() );


    // BOOKKEEPING
    VFace2DP* newVFace = createVFace( VFace2DP(m_VFaces.rbegin()->getID() + 1) );
    newGenerator->setAssignedVFace( newVFace );
    newVFace->setGenerator( newGenerator );
    newVFace->setFirstVEdge( *(newVEdges.begin()) );   

    for ( i = 1 ; i <= newVEdges.size() ; i++ ) {
        orderedNewVEdges[i]->setRightHand( orderedIntersectingVEdges[i] );
        orderedNewVEdges[i]->setRightLeg(  orderedIntersectingVEdges[i-1] );

        orderedNewVEdges[i]->setLeftHand(  orderedNewVEdges[i+1] );
        orderedNewVEdges[i]->setLeftLeg(   orderedNewVEdges[i-1] );

        orderedNewVEdges[i]->setRightFace( CCWVFacesOfOrderedVEdges[i-1] );
        orderedNewVEdges[i]->setLeftFace(  newVFace );
    }
}



void VoronoiDiagram2DP::connectCurrVDToNewVEdgeLoop( list<VEdge2DP*>& newVEdges )
{
    for( list<VEdge2DP*>::iterator it_newVEdges = newVEdges.begin() ; it_newVEdges != newVEdges.end() ; it_newVEdges++ ) {
        VEdge2DP* currEdge = *it_newVEdges;

        // set topology of right leg of current VEdge
        if( currEdge->getRightFace() == currEdge->getRightLeg()->getRightFace() ) {
            currEdge->getRightLeg()->setEndVertex( currEdge->getStartVertex() );
            currEdge->getRightLeg()->setRightHand( currEdge );
        }
        else {
            currEdge->getRightLeg()->setStartVertex( currEdge->getStartVertex() );
            currEdge->getRightLeg()->setLeftLeg(     currEdge );
        }


        // set topology of right hand of current VEdge
        if( currEdge->getRightFace() == currEdge->getRightHand()->getRightFace() ) {
            currEdge->getRightHand()->setRightLeg( currEdge );
        }
        else {
            currEdge->getRightHand()->setLeftHand( currEdge );
        }

        // reset first VEdge pointer for all VFaces incident to new VFace
        currEdge->getRightFace()->setFirstVEdge( currEdge );
    }
}



void VoronoiDiagram2DP::computeCoordOfNewVVertices( list<VVertex2DP*>& newVVertices )
{
    for( list<VVertex2DP*>::iterator it_vertex = newVVertices.begin() ; it_vertex != newVVertices.end() ; it_vertex++ ) {
        VVertex2DP* currVVertex = *it_vertex;

        list<Generator2DP*> definingGenerators;
        currVVertex->getDefining3Generators( definingGenerators );

        list<Generator2DP*>::iterator it_generator = definingGenerators.begin();
        Generator2DP* definingGenerator[3] = { *it_generator++, *it_generator++, *it_generator };

        rg_Point2D coordOfVVertex = computeCircumcenter( *definingGenerator[0]->getLocation(),
                                                         *definingGenerator[1]->getLocation(),
                                                         *definingGenerator[2]->getLocation() );
        currVVertex->setLocation( coordOfVVertex );
    }
}



VVertex2DP* VoronoiDiagram2DP::findFirstRedVertex( Generator2DP* const newGenerator, Generator2DP* const closestGenerator )
{
    list<VVertex2DP*> boundaryVVerticesOfFaceOfClosestGenerator;
    closestGenerator->getAssignedVFace()->getBoundaryVVertices( boundaryVVerticesOfFaceOfClosestGenerator );

    double      minHValue = DBL_MAX;
    VVertex2DP* firstRedVertex = NULL;

    for( list<VVertex2DP*>::iterator it_boundaryVertex = boundaryVVerticesOfFaceOfClosestGenerator.begin() ; it_boundaryVertex != boundaryVVerticesOfFaceOfClosestGenerator.end() ; it_boundaryVertex++ ) {
        VVertex2DP* currVVertex = *it_boundaryVertex;

        double currHValue = currVVertex->computeHValue( newGenerator );
        if( currHValue < minHValue ) {
            minHValue       = currHValue;
            firstRedVertex  = currVVertex;
        }
    }

    return firstRedVertex;
}



bool VoronoiDiagram2DP::hasCycleOccurred( VVertex2DP* const candidateForRedVVertex )
{
    int numRedVVertexAmongNeighbor = 0;
    list<VVertex2DP*> neighborVertices;
    candidateForRedVVertex->getAdjacent3VVertices( neighborVertices );

    for( list<VVertex2DP*>::iterator it_neighbor = neighborVertices.begin() ; it_neighbor != neighborVertices.end() ; it_neighbor++ ) {
        VVertex2DP* currNeighbor = *it_neighbor;

        if( currNeighbor->getStatus() == RED_V_P ) {
            numRedVVertexAmongNeighbor++;
        }
    }

    if( numRedVVertexAmongNeighbor > 1) {
        return true;
    }
    else {
        return false;
    }
}



bool VoronoiDiagram2DP::hasOldRegionSplit( VVertex2DP* const candidateForRedVVertex )
{
    list<VFace2DP*> incidentVFaces;
    candidateForRedVVertex->getIncident3VFaces( incidentVFaces );

    for( list<VFace2DP*>::iterator it_face = incidentVFaces.begin() ; it_face != incidentVFaces.end() ; it_face++ ) {
        VFace2DP* currVFace = *it_face;

        list<VVertex2DP*> boundaryVVerticesOfCurrVFace;
        currVFace->getBoundaryVVertices( boundaryVVerticesOfCurrVFace );



        // Count the total number of red VVertices on boundary of the current VFace
        list<VVertex2DP*>::iterator it_boundaryVVertex = boundaryVVerticesOfCurrVFace.begin();
        int numRedVVerticesOnBoundaryVVertices = 0;

        for( ; it_boundaryVVertex != boundaryVVerticesOfCurrVFace.end() ; it_boundaryVVertex++ ) {
            VVertex2DP* currBoundaryVVertex = *it_boundaryVVertex;
            if( currBoundaryVVertex->getStatus() == RED_V_P ) {
                numRedVVerticesOnBoundaryVVertices++;
            }
        }
        numRedVVerticesOnBoundaryVVertices++;



        // There is no split if the number of red VVertices is one.
        if( numRedVVerticesOnBoundaryVVertices == 1 ) {
            continue;
        }



        // Make the boundary vertices list start at the currBndVtx. 
        // 이 블럭을 함수로 따로
        int size = boundaryVVerticesOfCurrVFace.size();
        int i = 0;
        while ( i < size ) {
            VVertex2DP* currBndVtx = *boundaryVVerticesOfCurrVFace.begin();
            if( candidateForRedVVertex != currBndVtx ) {
                boundaryVVerticesOfCurrVFace.pop_front();
                boundaryVVerticesOfCurrVFace.push_back( currBndVtx );
            }
            else {
                break;
            }
        }

        it_boundaryVVertex = boundaryVVerticesOfCurrVFace.begin();
        it_boundaryVVertex++;
        int numConnectedRedVVertices = 1;



        // Count the connected red vertices CCW
        for( ; it_boundaryVVertex != boundaryVVerticesOfCurrVFace.end() ; it_boundaryVVertex++ ) {
            VVertex2DP* currBndVtx = *it_boundaryVVertex;
            if( currBndVtx->getStatus() == RED_V_P ) {
                numConnectedRedVVertices++;
            }
            else {
                break;
            }
        }

        it_boundaryVVertex = boundaryVVerticesOfCurrVFace.end();
        it_boundaryVVertex--;



        // Count the connected red vertices CW
        for( ; it_boundaryVVertex != boundaryVVerticesOfCurrVFace.begin() ; it_boundaryVVertex-- ) {
            VVertex2DP* currBndVtx = *it_boundaryVVertex;
            if( currBndVtx->getStatus() == RED_V_P ) {
                numConnectedRedVVertices++;
            }
            else {
                break;
            }
        }

        if( numRedVVerticesOnBoundaryVVertices != numConnectedRedVVertices ) {
            return true;
        }
    }

    return false;
}





rg_Point2D VoronoiDiagram2DP::computeCircumcenter( const rg_Point2D& pt1, const rg_Point2D& pt2, const rg_Point2D& pt3 )
{
    double x1   = pt1.getX();
    double y1   = pt1.getY();
    double x2   = pt2.getX();
    double y2   = pt2.getY();
    double x3   = pt3.getX();
    double y3   = pt3.getY();

    double a    = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
    double b_x  = (x1*x1 + y1*y1) * (y3 - y2) + (x2*x2 + y2*y2) * (y1 - y3) + (x3*x3 + y3*y3) * (y2 - y1);
    double b_y  = (x1*x1 + y1*y1) * (x2 - x3) + (x2*x2 + y2*y2) * (x3 - x1) + (x3*x3 + y3*y3) * (x1 - x2);

    double  XCoord = - b_x / (2*a);
    double  YCoord = - b_y / (2*a);

    return rg_Point2D(XCoord, YCoord);
}
