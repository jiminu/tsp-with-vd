#include "VoronoiDiagram2DC.h"
#include "rg_Point3D.h"
#include "Priority_Q_VEdge.h"
using namespace BULL2D::GeometryTier;



#include <vector>
#include <fstream>
#include <time.h>
using namespace std;





VoronoiDiagram2DC::VoronoiDiagram2DC(void)
{
}



VoronoiDiagram2DC::VoronoiDiagram2DC( list<rg_Circle2D>& circleset )
{
    int currID = 1;
    for( list<rg_Circle2D>::iterator it_point = circleset.begin() ; it_point != circleset.end() ; it_point++, currID++ ) {
        m_generators.push_back( Generator2D( currID, &(*it_point) ) );
    }
}



// VoronoiDiagram2DC::VoronoiDiagram2DC( const VoronoiDiagram2DC& VD2DC )
// {
// 
// }



VoronoiDiagram2DC::~VoronoiDiagram2DC(void)
{

}



void VoronoiDiagram2DC::getVoronoiEdges( list<const VEdge2D*>& VEdgesList ) const
{
    for( list<VEdge2D>::const_iterator it_edge = m_VEdges.begin() ; it_edge != m_VEdges.end() ; it_edge++ ) {
        const VEdge2D* currEdge = &(*it_edge);
        VEdgesList.push_back( currEdge );
    }
}



void VoronoiDiagram2DC::getVoronoiFaces( list<const VFace2D*>& VFacesList ) const
{
    for( list<VFace2D>::const_iterator it_face = m_VFaces.begin() ; it_face != m_VFaces.end() ; it_face++ ) {
        const VFace2D* currFace = &(*it_face);
        VFacesList.push_back( currFace );
    }
}



void VoronoiDiagram2DC::getVoronoiVertices( list<const VVertex2D*>& VVerticesList ) const
{
    for( list<VVertex2D>::const_iterator it_vertex = m_VVertices.begin() ; it_vertex != m_VVertices.end() ; it_vertex++ ) {
        const VVertex2D* currVertex = &(*it_vertex);
        VVerticesList.push_back( currVertex );
    }
}



void VoronoiDiagram2DC::getGenerators( list<const Generator2D*>& generatorsList ) const
{
    for( list<Generator2D>::const_iterator it_generator = m_generators.begin() ; it_generator != m_generators.end() ; it_generator++ ) {
        const Generator2D* currGenerator = &(*it_generator);
        generatorsList.push_back( currGenerator );
    }
}



void VoronoiDiagram2DC::getBoundedVoronoiEdges( list<VEdge2D*>& boundedEdgeList )
{
    for( list<VEdge2D>::iterator it_edge = m_VEdges.begin() ; it_edge != m_VEdges.end() ; it_edge++ ) {
        VEdge2D* currEdge = &(*it_edge);

        if( currEdge->isBounded() ) {
            boundedEdgeList.push_back( currEdge );
        }        
    }
}



void VoronoiDiagram2DC::createSortedGeneratorSet( list<rg_Circle2D>& circleset )
{

    circleset.sort(compareCircleDescendingorder);
    m_generators.clear();

    int currID = 1;
    for( list<rg_Circle2D>::iterator it_point = circleset.begin() ; it_point != circleset.end() ; it_point++, currID++ ) {
        m_generators.push_back( Generator2D( currID, &(*it_point) ) );
    }


}



// void VoronoiDiagram2DC::setGenerators( list<Generator2D>& generators )
// {
//     m_generators.clear();
//     m_generators = generators;
// }



void VoronoiDiagram2DC::constructVoronoiDiagram( list<rg_Circle2D>& circleset )
{
//     ofstream fout2("computing_with_sorting.txt");
//     clock_t startTime, endTime;
//     startTime = clock();

    createSortedGeneratorSet( circleset );
    
    if( IS_BUCKET_USED ) {
        m_bucket.setBucketCell( m_generators );
    }

    constructVoronoiDiagram();

//     endTime = clock();
//     double computingTime = endTime - startTime;
//     fout2 << computingTime <<endl;
//     fout2.close();
}



void VoronoiDiagram2DC::constructVoronoiDiagram()
{
    time2_1 = 0.0;
    time2_2 = 0.0;
    time2_3 = 0.0;
    time2_4 = 0.0;
    time2_5 = 0.0;
    time2_6 = 0.0;

    time2_2_1 = 0.0;
    time2_2_2 = 0.0;
    time2_2_3 = 0.0;
    time2_2_4 = 0.0;

    clock_t startTime, endTime;
    startTime = clock();

    insertPhantomGenerators();
    constructVoronoiDiagramForPhantomGenerators();

    endTime = clock();
    time1 = endTime - startTime;

    
    for( list<Generator2D>::iterator it_generator = m_generators.begin() ; it_generator != m_generators.end() ; it_generator++ ) {
        Generator2D* currGenerator = &(*it_generator);
        updateVoronoiDiagramWithNewGenerator( currGenerator );

        if( IS_BUCKET_USED ) {
            m_bucket.addGenerator( currGenerator );
        }
    }



    startTime = clock();

    removePhantomGenerators();

    endTime = clock();
    time3 = endTime - startTime;



    startTime = clock();

    removeAllExtraneousVVerticesAndVEdges();
        
    endTime = clock();
    time4 = endTime - startTime;
}



VFace2D* VoronoiDiagram2DC::findVFaceContainingQueryPoint( const rg_Point2D& pt )
{
    VFace2D*        VFaceContainingPoint = NULL;

//     int numOfTotalDisks = m_generators.size();
//     int numOfCurrentAddedDisks = m_VFaces.size() - 4;    
//     double currRationOfDiskAddition = (double)numOfCurrentAddedDisks / (double)numOfTotalDisks;
//     if( currRationOfDiskAddition > RATIO_OF_CHANGING_TO_BUCKET && IS_BUCKET_USED ) {

    if( IS_BUCKET_USED )
    {
        Generator2D* closeGenerator = m_bucket.findGoodCloseGenerator( pt );
        if( closeGenerator != NULL )
        {
            VFaceContainingPoint = closeGenerator->getAssignedVFace();
        }
        else
        {
            VFaceContainingPoint = &( *m_VFaces.rbegin() );
        }
    }
    else
    {
        VFaceContainingPoint = &( *m_VFaces.rbegin() );
    }



    rg_Circle2D*    currDisk = VFaceContainingPoint->getGenerator()->getDisk();
    double          minDistancePt2GeneratorOnVFace = pt.distance( currDisk->getCenterPt() ) - currDisk->getRadius();



    if( IS_INLINE_USED ) {
        VEdge2D*        beginningEdgeOfFace = VFaceContainingPoint->getFirstVEdge();
        VEdge2D*        currBoundaryEdge = beginningEdgeOfFace;
        VFace2D*        currOppositeFace = getFaceOfOppositeSide( VFaceContainingPoint, currBoundaryEdge );

        Generator2D*    neighborGenerator = NULL;
        double          distanceFromNeighbor2InputPt = DBL_MAX;


        do {
            currOppositeFace    = getFaceOfOppositeSide( VFaceContainingPoint, currBoundaryEdge );
            neighborGenerator   = currOppositeFace->getGenerator();

            if( neighborGenerator != NULL ) {
                currDisk = neighborGenerator->getDisk();
                distanceFromNeighbor2InputPt = pt.distance( currDisk->getCenterPt() ) - currDisk->getRadius();
            }

            if( distanceFromNeighbor2InputPt < minDistancePt2GeneratorOnVFace ) {
                VFaceContainingPoint        = currOppositeFace;
                minDistancePt2GeneratorOnVFace = distanceFromNeighbor2InputPt;
                beginningEdgeOfFace = currBoundaryEdge;
            }

            currBoundaryEdge = getCCWNextEdgeOnFace( currBoundaryEdge, VFaceContainingPoint );

        } while ( currBoundaryEdge != beginningEdgeOfFace );
    }


    else {
        bool            isThereCloserGeneratorToInputPoint = false;
        do {
            list<Generator2D*>  neighbor;
            double              distancePt2CurrNeighbor = 0.0;
            Generator2D*        currNeighbor = NULL;
            isThereCloserGeneratorToInputPoint = false;

            VFaceContainingPoint->getGenerator()->getNeighborGenerators( neighbor );

            for( list<Generator2D*>::iterator it_neighbor = neighbor.begin() ; it_neighbor != neighbor.end() ; it_neighbor++ ) {
                currNeighbor            = *it_neighbor;
                currDisk                = currNeighbor->getDisk();
                distancePt2CurrNeighbor = pt.distance( currDisk->getCenterPt() ) - currDisk->getRadius();

                if( distancePt2CurrNeighbor < minDistancePt2GeneratorOnVFace ) {
                    VFaceContainingPoint                = currNeighbor->getAssignedVFace();
                    minDistancePt2GeneratorOnVFace      = distancePt2CurrNeighbor;
                    isThereCloserGeneratorToInputPoint  = true;
                }
            }
        } while ( isThereCloserGeneratorToInputPoint );
    }

    return VFaceContainingPoint;
}



void VoronoiDiagram2DC::clear()
{
    m_VEdges.clear();
    m_VFaces.clear();
    m_VVertices.clear();
    m_generators.clear();
    m_phantomCircles.clear();
    m_phantomGenerators.clear();
}



VVertex2D* VoronoiDiagram2DC::createVVertex( const VVertex2D& vertex )
{
    m_VVertices.push_back( vertex );

    return &(*m_VVertices.rbegin());
}



VEdge2D* VoronoiDiagram2DC::createVEdge( const VEdge2D& edge )
{
    m_VEdges.push_back( edge );

    return &(*m_VEdges.rbegin());
}



VFace2D* VoronoiDiagram2DC::createVFace( const VFace2D& face )
{
    m_VFaces.push_back( face );

    return &(*m_VFaces.rbegin());
}



void VoronoiDiagram2DC::insertPhantomGenerators()
{
    double              maxRadius = 0.0;
    rg_BoundingBox2D    boundingBox;

    getBoundingBoxNMaxRadius( boundingBox, maxRadius);


    rg_Point2D coordinate[3];
    //computeSugiharaCoordOfPhantomGenerators( boundingBox, coordinate[0], coordinate[1], coordinate[2] );
    computeCoordOfPhantomGenerators_beforeOffset( boundingBox, coordinate[0], coordinate[1], coordinate[2] );
    
    
    transformCoordinatesToMakeOffsetTriangleByMaxRadius( maxRadius, coordinate[0], coordinate[1], coordinate[2] );


    setPhantomGenerators( maxRadius, coordinate[0], coordinate[1], coordinate[2] );
}



void VoronoiDiagram2DC::getBoundingBoxNMaxRadius( rg_BoundingBox2D& boundingBox, double& maxRadius )
{
    if( IS_BUCKET_USED ) {
        boundingBox.setMinPt( m_bucket.getMinPt() );
        boundingBox.setMaxPt( m_bucket.getMaxPt() );
    }
    else {
        for( list<Generator2D>::iterator it_generator = m_generators.begin() ; it_generator != m_generators.end() ; it_generator++ ) {
            rg_Circle2D currDisk = *it_generator->getDisk();
            boundingBox.updateBoxByAddingCircle( currDisk );
        }
    }
    
    maxRadius = (*m_generators.begin()).getDisk()->getRadius();
}



void VoronoiDiagram2DC::computeSugiharaCoordOfPhantomGenerators( const rg_BoundingBox2D& boundingBox, rg_Point2D& firstCoord, rg_Point2D& secondCoord, rg_Point2D& thirdCoord )
{
    double lengthOfX = boundingBox.evaluateXLength();
    double lengthOfY = boundingBox.evaluateYLength();
    double length = 0.0;

    if( lengthOfX > lengthOfY ) {
        length = lengthOfX;
    }
    else {
        length = lengthOfY;
    }

    rg_Point2D minPt = boundingBox.getMinPt();
    rg_Point2D maxPt = boundingBox.getMaxPt();

    firstCoord.setX( minPt.getX() + 0.5 * length );
    firstCoord.setY( minPt.getY() + (3.0 * sqrt(2.0)/2.0 + 0.5) * length );

    secondCoord.setX( minPt.getX() + (-3.0 * sqrt(6.0)/4.0 + 0.5) * length );
    secondCoord.setY( minPt.getY() + (-3.0 * sqrt(2.0)/4.0 + 0.5) * length );

    thirdCoord.setX( minPt.getX() + (3.0 * sqrt(6.0)/4.0 + 0.5) * length );
    thirdCoord.setY( minPt.getY() + (-3.0 * sqrt(2.0)/4.0 + 0.5) * length );

}



void VoronoiDiagram2DC::computeCoordOfPhantomGenerators_beforeOffset( const rg_BoundingBox2D& boundingBox, rg_Point2D& firstCoord, rg_Point2D& secondCoord, rg_Point2D& thirdCoord )
{
    double      expansionRatio                      = 10.1;
    double      radiusOfCircumcircleOfBoundingBox   = boundingBox.evaluateLongestLength() / 2.0;
    double      expandedRadius                      = radiusOfCircumcircleOfBoundingBox * expansionRatio;
    rg_Point2D  centerPtOfBoundingBox               = boundingBox.getCenterPt();

    firstCoord.setX( centerPtOfBoundingBox.getX() );
    firstCoord.setY( centerPtOfBoundingBox.getY() + expandedRadius );

    secondCoord.setX( centerPtOfBoundingBox.getX() - ( expandedRadius * sqrt(3.0) / 2.0 ) );
    secondCoord.setY( centerPtOfBoundingBox.getY() - ( expandedRadius / 2.0 ) );

    thirdCoord.setX( centerPtOfBoundingBox.getX() + ( expandedRadius * sqrt(3.0) / 2.0 ) );
    thirdCoord.setY( centerPtOfBoundingBox.getY() - ( expandedRadius / 2.0 ) );
}


void VoronoiDiagram2DC::transformCoordinatesToMakeOffsetTriangleByMaxRadius( const double& maxRadius, rg_Point2D& firstCoord, rg_Point2D& secondCoord, rg_Point2D& thirdCoord )
{
    firstCoord  = firstCoord + rg_Point2D( 0.0, 2.0 * maxRadius );
    secondCoord = secondCoord + rg_Point2D( -sqrt(3.0) * maxRadius, -maxRadius );
    thirdCoord  = thirdCoord + rg_Point2D( sqrt(3.0) * maxRadius, -maxRadius );
}



void VoronoiDiagram2DC::setPhantomGenerators( const double& maxRadius, rg_Point2D& firstCoord, rg_Point2D& secondCoord, rg_Point2D& thirdCoord )
{
    m_phantomCircles.push_back( rg_Circle2D( firstCoord, maxRadius ) );
    m_phantomCircles.push_back( rg_Circle2D( secondCoord, maxRadius ) );
    m_phantomCircles.push_back( rg_Circle2D( thirdCoord, maxRadius ) );

    list<rg_Circle2D>::iterator it_phantom = m_phantomCircles.begin();
    rg_Circle2D* phantomCircle[3] = { &(*it_phantom++), &(*it_phantom++), &(*it_phantom) };

    m_phantomGenerators.push_back( Generator2D( -2, phantomCircle[0] ) );
    m_phantomGenerators.push_back( Generator2D( -1, phantomCircle[1] ) );
    m_phantomGenerators.push_back( Generator2D(  0, phantomCircle[2] ) );
}



void VoronoiDiagram2DC::constructVoronoiDiagramForPhantomGenerators()
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


    list<Generator2D>::iterator it_phantoms = m_phantomGenerators.begin();
    Generator2D* phantoms[3] = {&(*it_phantoms++), &(*it_phantoms++), &(*it_phantoms)};


    //  1. make elements of voronoi diagram of phantom generators

    VVertex2D* vertices[4]   = { createVVertex( VVertex2D(0) ),        // vertex defined by three phantom generators
        createVVertex( VVertex2D(1) ),        // infinite vertex on edge0
        createVVertex( VVertex2D(2) ),        // infinite vertex on edge1
        createVVertex( VVertex2D(3) ) };      // infinite vertex on edge2

    VEdge2D*   edges[6]      = { createVEdge( VEdge2D(0) ),            // edge made by phantom0 and phantom1     
        createVEdge( VEdge2D(1) ),            // edge made by phantom1 and phantom2
        createVEdge( VEdge2D(2) ),            // edge made by phantom2 and phantom0
        createVEdge( VEdge2D(3) ),            // artificial boundary for unbounded face1
        createVEdge( VEdge2D(4) ),            // artificial boundary for unbounded face2 
        createVEdge( VEdge2D(5) ) };          // artificial boundary for unbounded face3 

    VFace2D* faces[4]        = { createVFace( VFace2D(0) ),            // infinite space defined by artificial boundary for voronoi diagram
        createVFace( VFace2D(1) ),            // face containing phantom0 
        createVFace( VFace2D(2) ),            // face containing phantom1
        createVFace( VFace2D(3) ) };          // face containing phantom2



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

    rg_Circle2D tangentCircle1, tangentCircle2;
    rg_Circle2D::makeCircumcircle( *phantoms[0]->getDisk(), *phantoms[1]->getDisk(), *phantoms[2]->getDisk(), tangentCircle1, tangentCircle2);
    rg_Point2D locationOfVertex0 = tangentCircle1.getCenterPt();
    vertices[0]->setLocation(locationOfVertex0);


    //      3.2. define coordinate of infinite vertices

    rg_Point2D centerOfPoint0And1 = ( phantoms[0]->getDisk()->getCenterPt() + phantoms[1]->getDisk()->getCenterPt() ) / 2.0;
    rg_Point2D centerOfPoint1And2 = ( phantoms[1]->getDisk()->getCenterPt() + phantoms[2]->getDisk()->getCenterPt() ) / 2.0;
    rg_Point2D centerOfPoint2And0 = ( phantoms[2]->getDisk()->getCenterPt() + phantoms[0]->getDisk()->getCenterPt() ) / 2.0;

    rg_Point2D vectorOfVertex1 = centerOfPoint0And1 - locationOfVertex0;
    rg_Point2D vectorOfVertex2 = centerOfPoint1And2 - locationOfVertex0;
    rg_Point2D vectorOfVertex3 = centerOfPoint2And0 - locationOfVertex0;

    rg_Point2D locationOfVertex1 = centerOfPoint0And1 + 5.0 * vectorOfVertex1;
    rg_Point2D locationOfVertex2 = centerOfPoint1And2 + 5.0 * vectorOfVertex2;
    rg_Point2D locationOfVertex3 = centerOfPoint2And0 + 5.0 * vectorOfVertex3;

    vertices[1]->setLocation(locationOfVertex1);
    vertices[2]->setLocation(locationOfVertex2);
    vertices[3]->setLocation(locationOfVertex3); 
}



void VoronoiDiagram2DC::updateVoronoiDiagramWithNewGenerator( Generator2D* const newGenerator )
{

    
    clock_t startTime = clock();
    Generator2D* closestGenerator = findClosestGeneratorToNewGenerator( newGenerator );
    
    clock_t endTime = clock();
    time2_1 += endTime - startTime;


    startTime = clock();
    list<VEdge2D*>   intersectingVEdges;
    list<VVertex2D*> fictitiousVVertices;
    findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop_EdgeSplitPossible( newGenerator, closestGenerator, intersectingVEdges, fictitiousVVertices );

    endTime = clock();
    time2_2 += endTime - startTime;



    startTime = clock();
    list<VVertex2D*> newVVertices;
    list<VEdge2D*>   newVEdges;
    makeVEdgeLoopForNewGeneratorAndConnectToCurrVD( newGenerator, intersectingVEdges, newVVertices, newVEdges );
    
    endTime = clock();
    time2_3 += endTime - startTime;



    startTime = clock();
    connectCurrVDToNewVEdgeLoop( newVEdges );
    
    endTime = clock();
    time2_4 += endTime - startTime;



    startTime = clock();
    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );
    
    endTime = clock();
    time2_5 += endTime - startTime;



    startTime = clock();
    computeCoordOfNewVVertices( newVVertices );
    
    endTime = clock();
    time2_6 += endTime - startTime;
}



Generator2D* VoronoiDiagram2DC::findClosestGeneratorToNewGenerator( Generator2D* const newGenerator )
{
    VFace2D* closestVFace = findVFaceContainingQueryPoint( newGenerator->getDisk()->getCenterPt() );
    return closestVFace->getGenerator();
}



void VoronoiDiagram2DC::removeAllExtraneousVVerticesAndVEdges()
{
    // while 문으로...
    for( list<VVertex2D>::iterator it_vertex = m_VVertices.begin() ; it_vertex != m_VVertices.end() ; ) { 
        VVertex2D currVertex = *it_vertex;
        if( currVertex.getStatus() == RED_V || currVertex.getStatus() == YELLOW_V ) {
            it_vertex = m_VVertices.erase( it_vertex );
        }
        else {
            it_vertex++;
        }
    }

    for( list<VEdge2D>::iterator it_edge = m_VEdges.begin() ; it_edge != m_VEdges.end() ;) {
        VEdge2D currEdge = *it_edge;
        if( currEdge.getStatus() == RED_E ) {
            it_edge = m_VEdges.erase( it_edge );
        }
        else {
            it_edge++;
        }
    }
}



void VoronoiDiagram2DC::findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop_EdgeSplitPossible( Generator2D* const newGenerator, Generator2D* const closestGenerator, list<VEdge2D*>& intersectingVEdges, list<VVertex2D*>& fictitiousVVertices )
{


    clock_t startTime = clock();
    list<VVertex2D*> redVVertices;
    VVertex2D*  firstRedVertex = findFirstRedVertex( newGenerator, closestGenerator );
       

    firstRedVertex->setStatus( RED_V );
    redVVertices.push_back( firstRedVertex );

    clock_t endTime = clock();
    time2_2_1 += endTime - startTime;


    startTime = clock();
    list<VVertex2D*> blueVVertices;
    wavePropagation_ver1( newGenerator, redVVertices, blueVVertices, fictitiousVVertices );

    endTime = clock();
    time2_2_2 += endTime - startTime;


    startTime = clock();
    findCrossingEdges( redVVertices, intersectingVEdges );
    
    endTime = clock();
    time2_2_3 += endTime - startTime;


    startTime = clock();
    reconfigureByConvertingBlueVVerticesToWhite( blueVVertices );
    
    endTime = clock();
    time2_2_4 += endTime - startTime;

}



void VoronoiDiagram2DC::findAnomalizingEdgeAmongIncidentVEdgesAndSplit( VVertex2D* const currRedVVertex, Generator2D* const newGenerator, list<VVertex2D*>& fictitiousVVertices, list<VEdge2D*>& anomalyTestedVEdges )
{
    if( IS_INLINE_USED ) {
        VEdge2D* beginningEdgeOfVertex = currRedVVertex->getFirstVEdge();
        VEdge2D* nextIncidentEdge = beginningEdgeOfVertex;
        VEdge2D* currIncidentEdge = NULL;

        do 
        {
            currIncidentEdge = nextIncidentEdge;
            nextIncidentEdge = getCCWNextEdgeOnVertex( currIncidentEdge, currRedVVertex );

            if( currIncidentEdge->isAnomalyTestDone() == false ) {
                currIncidentEdge->setTrueAnomalyTestDone();
                anomalyTestedVEdges.push_back( currIncidentEdge );

                if( isAnomalizingEdge( currIncidentEdge, newGenerator ) ) {
                    VVertex2D* fictitiousVVertex = splitVEdgeAtFictitiousVVertex( currIncidentEdge );
                    beginningEdgeOfVertex = currRedVVertex->getFirstVEdge();
                    fictitiousVVertices.push_back( fictitiousVVertex );
                }
            }

        } while ( nextIncidentEdge != beginningEdgeOfVertex );
    }

    else {
        list<VEdge2D*> incidentVEdges;
        currRedVVertex->getIncident3VEdges( incidentVEdges );

        for( list<VEdge2D*>::iterator it_incidentVEdge = incidentVEdges.begin() ; it_incidentVEdge != incidentVEdges.end() ; it_incidentVEdge++ ) {
            VEdge2D* currVEdge = *it_incidentVEdge;
            if( currVEdge->isAnomalyTestDone() == false ) {
                currVEdge->setTrueAnomalyTestDone();
                anomalyTestedVEdges.push_back( currVEdge );

                if( isAnomalizingEdge( currVEdge, newGenerator ) ) {
                    VVertex2D* fictitiousVVertex = splitVEdgeAtFictitiousVVertex( currVEdge );
                    fictitiousVVertices.push_back( fictitiousVVertex );
                }
            }
        }
    }
    
}



void VoronoiDiagram2DC::findRedVVertexAmongNeighbor( VVertex2D* const currRedVVertex, Generator2D* const newGenerator, list<VVertex2D*>& blueVVertices, list<VVertex2D*>& redVVertices, list<VVertex2D*>& fictitiousVVertices, list<VEdge2D*>& anomalyTestedVEdges )
{
    if( IS_INLINE_USED ) {
        VVertex2D* currNeighbor = NULL;
        VEdge2D*   beginningEdgeOfVertex = currRedVVertex->getFirstVEdge();
        VEdge2D*   nextIncidentEdge = beginningEdgeOfVertex;
        VEdge2D*   currIncidentEdge = NULL;

        do 
        {
            currIncidentEdge = nextIncidentEdge;
            nextIncidentEdge = getCCWNextEdgeOnVertex( currIncidentEdge, currRedVVertex );

            currNeighbor = currIncidentEdge->getOppositVVertex( currRedVVertex );

            if( currNeighbor->getStatus() != WHITE_V ) {
                continue;
            }

            double MUValue = currNeighbor->computeMUValue( newGenerator );

            if( MUValue >= 0) {
                currNeighbor->setStatus( BLUE_V );
                blueVVertices.push_back( currNeighbor );
                continue;
            }

            if( hasCycleOccurred( currNeighbor ) ) {
                findAnomalizingEdgeAmongIncidentVEdgesAndSplit( currNeighbor, newGenerator, fictitiousVVertices, anomalyTestedVEdges );
                if( hasCycleOccurred( currNeighbor ) ) {
                    currNeighbor->setStatus( BLUE_V );
                    blueVVertices.push_back( currNeighbor );
                    continue;
                }

            }

            if( hasOldRegionSplit( currNeighbor ) ) {  
                currNeighbor->setStatus( BLUE_V );
                blueVVertices.push_back( currNeighbor );
                continue;
            }

            currNeighbor->setStatus( RED_V );
            redVVertices.push_back( currNeighbor );

        } while ( nextIncidentEdge != beginningEdgeOfVertex );
    }

    else {
        list<VVertex2D*> neighborsOfCurrRedVVertex;
        currRedVVertex->getAdjacent3VVertices( neighborsOfCurrRedVVertex );

        for( list<VVertex2D*>::iterator it_neighbor = neighborsOfCurrRedVVertex.begin() ; it_neighbor != neighborsOfCurrRedVVertex.end() ; it_neighbor++ ) {
            VVertex2D* currNeighbor = *it_neighbor;

            if( currNeighbor->getStatus() != WHITE_V ) {
                continue;
            }

            double MUValue = currNeighbor->computeMUValue( newGenerator );

            if( MUValue >= 0) {
                currNeighbor->setStatus( BLUE_V );
                blueVVertices.push_back( currNeighbor );
                continue;
            }

            if( hasCycleOccurred( currNeighbor ) ) {
                findAnomalizingEdgeAmongIncidentVEdgesAndSplit( currNeighbor, newGenerator, fictitiousVVertices, anomalyTestedVEdges );
                if( hasCycleOccurred( currNeighbor ) ) {
                    currNeighbor->setStatus( BLUE_V );
                    blueVVertices.push_back( currNeighbor );
                    continue;
                }

            }

            if( hasOldRegionSplit( currNeighbor ) ) {   
                currNeighbor->setStatus( BLUE_V );
                blueVVertices.push_back( currNeighbor );
                continue;
            }

            currNeighbor->setStatus( RED_V );
            redVVertices.push_back( currNeighbor );
        }
    }

}



void VoronoiDiagram2DC::findCrossingEdges( list<VVertex2D*>& redVVertices, list<VEdge2D*>& intersectingVEdges )
{
    list<VVertex2D*>::iterator it_redVVertex = redVVertices.begin();
    for( ; it_redVVertex != redVVertices.end() ; it_redVVertex++ ) {
        VVertex2D* currRedVVertex = *it_redVVertex;
        list<VEdge2D*> incidentVEdges;
        currRedVVertex->getIncident3VEdges( incidentVEdges );

        for( list<VEdge2D*>::iterator it_incidentVEdge = incidentVEdges.begin() ; it_incidentVEdge != incidentVEdges.end() ; it_incidentVEdge++ ) {
            VEdge2D* currIncidentVEdge = *it_incidentVEdge;

            if( currIncidentVEdge->getStatus() != WHITE_E ) {
                continue;
            }

            if( currIncidentVEdge->getStartVertex()->getStatus() != currIncidentVEdge->getEndVertex()->getStatus() ) {
                currIncidentVEdge->setStatus( PURPLE_E ); // PURPLE_E
                intersectingVEdges.push_back( currIncidentVEdge );
            }
            else {
                currIncidentVEdge->setStatus( RED_E );
            }
        }
    }
}



void VoronoiDiagram2DC::reconfigureByConvertingBlueVVerticesToWhite( list<VVertex2D*>& blueVVertices )
{
    for( list<VVertex2D*>::iterator it_blueVVertex = blueVVertices.begin() ; it_blueVVertex != blueVVertices.end() ; it_blueVVertex++ ) {
        VVertex2D* currBlueVVertex = *it_blueVVertex;
        currBlueVVertex->setStatus( WHITE_V );
    }
}



void VoronoiDiagram2DC::makeVEdgeLoopForNewGeneratorAndConnectToCurrVD( Generator2D* const newGenerator, list<VEdge2D*>& intersectingVEdges, list<VVertex2D*>& newVVertices, list<VEdge2D*>& newVEdges )
{
    int currPos = 0;
    vector<VEdge2D*> orderedIntersectingVEdges;
    vector<VFace2D*> CCWVFacesOfOrderedVEdges;
    vector<VVertex2D*> newVVerticesOnOrderedVEdges;

    // make intersecting VEdge list in CCW order and create newVVertices on the ordered VEdges
    orderedIntersectingVEdges.push_back( *intersectingVEdges.begin() );

    while( currPos < intersectingVEdges.size()) { 
        VEdge2D* currEdge = orderedIntersectingVEdges[currPos];

        VFace2D* CCWFace;

        if( currEdge->getStartVertex()->getStatus() == RED_V ) {
//        if( currEdge->getRightFace() == CCWFace ) {
            CCWFace = currEdge->getLeftFace();
        }
        else {
            CCWFace = currEdge->getRightFace();
        }

        CCWVFacesOfOrderedVEdges.push_back( CCWFace );
        VVertex2D* currNewVVertex = createVVertex( VVertex2D(m_VVertices.rbegin()->getID()+1) );
        currNewVVertex->setStatus( PURPLE_V );                                                          // make new VVertex purple to know which vertex is new VVertex without coordinate
        newVVerticesOnOrderedVEdges.push_back( currNewVVertex );

        newVVertices.push_back( newVVerticesOnOrderedVEdges[currPos] );


        do 
        {
            if( CCWFace == currEdge->getLeftFace() ) {
                currEdge = currEdge->getLeftLeg();
            }
            else {
                currEdge = currEdge->getRightHand();
            }
        } while ( currEdge->getStatus() != PURPLE_E );

        orderedIntersectingVEdges.push_back( currEdge );
        currPos++;
    }
    newVVerticesOnOrderedVEdges.push_back(  newVVerticesOnOrderedVEdges[0] );


    // reconfigure intersecting VEdges WHITE_E
    for( list<VEdge2D*>::iterator it_intersectingEdge = intersectingVEdges.begin() ; it_intersectingEdge != intersectingVEdges.end() ; it_intersectingEdge++ ) {
        VEdge2D* currIntersectingEdge = *it_intersectingEdge;
        currIntersectingEdge->setStatus( WHITE_E );
    }


    // create new VEdges and connect new VEdges to new VVertices
    int i = 0;
    for ( i = 0 ; i < newVVerticesOnOrderedVEdges.size()-1 ; i++ ) {
        if( orderedIntersectingVEdges[i]->getStartVertex()->getStatus() == YELLOW_V ) {
            newVVerticesOnOrderedVEdges[i]->setFirstVEdge( orderedIntersectingVEdges[i]->getLeftLeg() );
        }
        else {
            newVVerticesOnOrderedVEdges[i]->setFirstVEdge( orderedIntersectingVEdges[i] );
        }
        

        VEdge2D* newEdge = createVEdge( VEdge2D(m_VEdges.rbegin()->getID() + 1) );
        newEdge->setStartVertex( newVVerticesOnOrderedVEdges[i] );
        newEdge->setEndVertex(   newVVerticesOnOrderedVEdges[i+1] );
        newVEdges.push_back( newEdge );    
    }

    vector<VEdge2D*> orderedNewVEdges;
    i = 0;
    orderedNewVEdges.push_back( *newVEdges.rbegin() );
    list<VEdge2D*>::iterator i_edge;
    for ( i_edge= newVEdges.begin(); i_edge != newVEdges.end(); i_edge++, i++ ) {
        orderedNewVEdges.push_back( *i_edge );
    }
    orderedNewVEdges.push_back( *newVEdges.begin() );


    // BOOKKEEPING
    VFace2D* newVFace = createVFace( VFace2D(m_VFaces.rbegin()->getID() + 1) );
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



void VoronoiDiagram2DC::connectCurrVDToNewVEdgeLoop( list<VEdge2D*>& newVEdges )
{
    for( list<VEdge2D*>::iterator it_newVEdges = newVEdges.begin() ; it_newVEdges != newVEdges.end() ; it_newVEdges++ ) {
        VEdge2D* currEdge = *it_newVEdges;

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



void VoronoiDiagram2DC::computeCoordOfNewVVertices( list<VVertex2D*>& newVVertices )
{
    /*
    while( newVVertices.size() != 0) {
        VVertex2D* currVVertex = *newVVertices.begin();
        newVVertices.pop_front();

        list<Generator2D*> definingGenerators;
        currVVertex->getDefining3Generators( definingGenerators );

        list<Generator2D*>::iterator it_generator = definingGenerators.begin();
        Generator2D definingGenerator[3] = { **it_generator++, **it_generator++, **it_generator };

        rg_Circle2D firstCircle;
        rg_Circle2D secondCircle;

        int numCircumcircles = 0;
        numCircumcircles = rg_Circle2D::makeCircumcircle( *definingGenerator[0].getDisk(),
                                                          *definingGenerator[1].getDisk(),
                                                          *definingGenerator[2].getDisk(),
                                                          firstCircle,
                                                          secondCircle );


        switch( numCircumcircles )  {
        case 0:
            cout << "Cannot compute coordinate of newVVertex!" << endl;
            exit(1);
            break;


        case 1:
            currVVertex->setLocation( firstCircle.getCenterPt() );
            currVVertex->setStatus( WHITE_V );
            break;


        case 2:
            // Check current VVertex has white neighbor VVertex to compute coordinate via coordinate of white neighbor
            // except anomaly vertex
            list<VEdge2D*> incidentVEdges;
            currVVertex->getIncident3VEdges( incidentVEdges );

            rg_BOOL  isThereWhiteNeighborVertex                    = false;
            rg_BOOL  isWhiteNeighborVertexDefinedBySame3Generators = false;
            VEdge2D* referenceEdgeForComputingCoordinate = rg_NULL;

            list<VEdge2D*>::iterator it_incidentVEdge = incidentVEdges.begin();
            for(; it_incidentVEdge != incidentVEdges.end() ; it_incidentVEdge++ ) {
                VEdge2D* currEdge = *it_incidentVEdge;

                if( currEdge->getStartVertex()->getStatus() != currEdge->getEndVertex()->getStatus() ) {               
                    isThereWhiteNeighborVertex          = true;
                    referenceEdgeForComputingCoordinate = currEdge;

                    if( areTwoVerticesDefinedBySameGeneratorTriplet( currEdge->getStartVertex(), currEdge->getEndVertex() ) ) {
                        isWhiteNeighborVertexDefinedBySame3Generators = true;
                    }
                    break;
                }
            }            


            if( isThereWhiteNeighborVertex ) {
               VVertex2D* referenceVertex = NULL;
                if( referenceEdgeForComputingCoordinate->getStartVertex() == currVVertex ) {
                    referenceVertex = referenceEdgeForComputingCoordinate->getEndVertex();
                }
                else {
                    referenceVertex = referenceEdgeForComputingCoordinate->getStartVertex();
                }

                if( !isWhiteNeighborVertexDefinedBySame3Generators ) {

                    rg_Point2D  centerOfReferenceGenerator;
                    rg_Point2D  baseVector;
                    if( referenceEdgeForComputingCoordinate->getStartVertex() == currVVertex ) {
                        centerOfReferenceGenerator = referenceEdgeForComputingCoordinate->getRightFace()->getGenerator()->getDisk()->getCenterPt();
                        baseVector = referenceEdgeForComputingCoordinate->getEndVertex()->getLocation() - centerOfReferenceGenerator;
                    }
                    else {
                        centerOfReferenceGenerator = referenceEdgeForComputingCoordinate->getLeftFace()->getGenerator()->getDisk()->getCenterPt();
                        baseVector = referenceEdgeForComputingCoordinate->getStartVertex()->getLocation() - centerOfReferenceGenerator;
                    }
                    rg_Point2D vecToFirstCircle  = firstCircle.getCenterPt()  - centerOfReferenceGenerator;
                    rg_Point2D vecToSecondCircle = secondCircle.getCenterPt() - centerOfReferenceGenerator;


                    rg_REAL angleToFirstCircle  = angleFromVec1toVec2( baseVector, vecToFirstCircle );
                    rg_REAL angleToSecondCircle = angleFromVec1toVec2( baseVector, vecToSecondCircle );

                    if( angleToFirstCircle < angleToSecondCircle ) {
                        currVVertex->setLocation( firstCircle.getCenterPt() );
                    }
                    else {
                        currVVertex->setLocation( secondCircle.getCenterPt() );
                    }
                    currVVertex->setStatus( WHITE_V );
                    break;
                }
                else {
                    double distanceFromRefVtxToFirstCircle  = referenceVertex->getLocation().distance( firstCircle.getCenterPt() );
                    double distanceFromRefVtxToSecondCircle = referenceVertex->getLocation().distance( secondCircle.getCenterPt() );

                    if( distanceFromRefVtxToFirstCircle > distanceFromRefVtxToSecondCircle ) {
                        currVVertex->setLocation( firstCircle.getCenterPt() );
                    }
                    else {
                        currVVertex->setLocation( secondCircle.getCenterPt() );
                    }
                    currVVertex->setStatus( WHITE_V );

                }
            }     
            else {
                newVVertices.push_back(currVVertex);
            }

            break;
        }
    }
    */

    while( newVVertices.size() != 0) {
        VVertex2D* currVVertex = *newVVertices.begin();
        newVVertices.pop_front();

        list<Generator2D*> definingGenerators;
        currVVertex->getDefining3Generators( definingGenerators );

        list<Generator2D*>::iterator it_generator = definingGenerators.begin();
        Generator2D definingGenerator[3] = { **it_generator++, **it_generator++, **it_generator };

        rg_Circle2D firstCircle;
        rg_Circle2D secondCircle;

        int numCircumcircles = 0;
        numCircumcircles = rg_Circle2D::makeCircumcircle( *definingGenerator[0].getDisk(),
            *definingGenerator[1].getDisk(),
            *definingGenerator[2].getDisk(),
            firstCircle,
            secondCircle );


        switch( numCircumcircles )
        {

        case 0:
            cout << "Cannot compute coordinate of newVVertex!" << endl;
            exit(1);
            break;


        case 1:
            currVVertex->setLocation( firstCircle.getCenterPt() );
            currVVertex->setStatus( WHITE_V );
            break;


        case 2:
            // Check current VVertex has white neighbor VVertex to compute coordinate via coordinate of white neighbor
            // except anomaly vertex
            list<VEdge2D*> incidentVEdges;
            currVVertex->getIncident3VEdges( incidentVEdges );

            rg_BOOL hasWhiteNeighborVVertex = false;
            VEdge2D* currIncidentEdge = rg_NULL;

            list<VEdge2D*>::iterator it_incidentVEdge = incidentVEdges.begin();
            for(; it_incidentVEdge != incidentVEdges.end() ; it_incidentVEdge++ ) {

                currIncidentEdge = *it_incidentVEdge;

                // 양 끝의 색이 같다면 
                if( currIncidentEdge->getStartVertex()->getStatus() == currIncidentEdge->getEndVertex()->getStatus() ) {
                    continue;
                }
                // 다르다면 
                else {
                    // 같은 
                    if( areTwoVerticesDefinedBySameGeneratorTriplet( currIncidentEdge->getStartVertex(), currIncidentEdge->getEndVertex() ) ) {
                        // definingGeneratorCircles
                        continue;
                    }
                    else {
                        hasWhiteNeighborVVertex = true;
                        break;
                    }

                }

            }            


            if( !hasWhiteNeighborVVertex ) {
                newVVertices.push_back(currVVertex);
                continue;
            }

            rg_Point2D  centerOfReferenceGenerator;
            rg_Point2D  baseLine;
            rg_Point2D  firstLine;
            rg_Point2D  secondLine;

            if( currIncidentEdge->getStartVertex() == currVVertex ) {
                centerOfReferenceGenerator = currIncidentEdge->getRightFace()->getGenerator()->getDisk()->getCenterPt();
                baseLine = currIncidentEdge->getEndVertex()->getLocation() - centerOfReferenceGenerator;
            }
            else {
                centerOfReferenceGenerator = currIncidentEdge->getLeftFace()->getGenerator()->getDisk()->getCenterPt();
                baseLine = currIncidentEdge->getStartVertex()->getLocation() - centerOfReferenceGenerator;
            }

            firstLine = firstCircle.getCenterPt() - centerOfReferenceGenerator;
            secondLine = secondCircle.getCenterPt() - centerOfReferenceGenerator;

            rg_REAL firstAngle = angleFromVec1toVec2( baseLine, firstLine );
            rg_REAL secondAngle = angleFromVec1toVec2( baseLine, secondLine );

            if( firstAngle < secondAngle ) {
                currVVertex->setLocation( firstCircle.getCenterPt() );
            }
            else {
                currVVertex->setLocation( secondCircle.getCenterPt() );
            }

            currVVertex->setStatus( WHITE_V );
            break;

        }

    }
}



VVertex2D* VoronoiDiagram2DC::findFirstRedVertex( Generator2D* const newGenerator, Generator2D* const closestGenerator )
{
    list<VVertex2D*> boundaryVVerticesOfFaceOfClosestGenerator;
    closestGenerator->getAssignedVFace()->getBoundaryVVertices( boundaryVVerticesOfFaceOfClosestGenerator );

    double      minMUValue = DBL_MAX;
    VVertex2D*  firstRedVertex = NULL;

    for( list<VVertex2D*>::iterator it_boundaryVertex = boundaryVVerticesOfFaceOfClosestGenerator.begin() ; it_boundaryVertex != boundaryVVerticesOfFaceOfClosestGenerator.end() ; it_boundaryVertex++ ) {
        VVertex2D* currVVertex = *it_boundaryVertex;

        double currHValue = currVVertex->computeMUValue( newGenerator );
        if( currHValue < minMUValue ) {
            minMUValue       = currHValue;
            firstRedVertex  = currVVertex;
        }
    }

    if( minMUValue >= 0) {
        list<VEdge2D*> boundaryVEdges;
        closestGenerator->getAssignedVFace()->getBoundaryVEdges( boundaryVEdges );

        for( list<VEdge2D*>::iterator it_boundaryEdge = boundaryVEdges.begin() ; it_boundaryEdge != boundaryVEdges.end() ; it_boundaryEdge++ ) {
            VEdge2D* currVEdge = *it_boundaryEdge;
            if( isAnomalizingEdge( currVEdge, newGenerator ) ) {
                VVertex2D* fictitiousVVertex = splitVEdgeAtFictitiousVVertex( currVEdge );
                firstRedVertex = fictitiousVVertex;
                break;
            }
        }
    }

    return firstRedVertex;
}



bool VoronoiDiagram2DC::hasCycleOccurred( VVertex2D* const candidateForRedVVertex )
{
    int numRedVVertexAmongNeighbor = 0;

    if( IS_INLINE_USED ) {
        VVertex2D* currNeighbor = NULL;
        VEdge2D*   beginningEdgeOfVertex = candidateForRedVVertex->getFirstVEdge();
        VEdge2D*   nextIncidentEdge = beginningEdgeOfVertex;
        VEdge2D*   currIncidentEdge = NULL;

        do 
        {
            currIncidentEdge = nextIncidentEdge;
            nextIncidentEdge = getCCWNextEdgeOnVertex( currIncidentEdge, candidateForRedVVertex );

            currNeighbor = currIncidentEdge->getOppositVVertex( candidateForRedVVertex );

            if( currNeighbor->getStatus() == RED_V ) {
                numRedVVertexAmongNeighbor++;
            }
        }while( nextIncidentEdge != beginningEdgeOfVertex );
    }

    else {
        list<VVertex2D*> neighborVertices;
        candidateForRedVVertex->getAdjacent3VVertices( neighborVertices );

        for( list<VVertex2D*>::iterator it_neighbor = neighborVertices.begin() ; it_neighbor != neighborVertices.end() ; it_neighbor++ ) {
            VVertex2D* currNeighbor = *it_neighbor;

            if( currNeighbor->getStatus() == RED_V ) {
                numRedVVertexAmongNeighbor++;
            }
        }
    }
    

    if( numRedVVertexAmongNeighbor > 1) {
        return true;
    }
    else {
        return false;
    }
}



bool VoronoiDiagram2DC::hasOldRegionSplit( VVertex2D* const candidateForRedVVertex )
{


    if( IS_INLINE_USED ) {
        VEdge2D* beginningEdgeOfVertex = candidateForRedVVertex->getFirstVEdge();
        VFace2D* currIncidentFace = NULL;
        VEdge2D* currIncidentEdge = NULL;
        VEdge2D* nextIncidentEdge = beginningEdgeOfVertex;

        do 
        {
            currIncidentEdge = nextIncidentEdge;
            nextIncidentEdge = getCCWNextEdgeOnVertex( currIncidentEdge, candidateForRedVVertex );

            if( candidateForRedVVertex == currIncidentEdge->getStartVertex() ) {
                currIncidentFace = currIncidentEdge->getLeftFace();
            }
            else {
                currIncidentFace = currIncidentEdge->getRightFace();
            }



            candidateForRedVVertex->setStatus( RED_V );


            int     numComponentsOfBoundaryVVertices = 0;
            bool    isCurrVVertexRed = false;
            bool    isPrevVVertexRed = false;

            VEdge2D*   beginningEdgeOfIncidentFace = currIncidentFace->getFirstVEdge();
            VVertex2D* firstBoundaryVertex = NULL;

            VVertex2D* currBoundaryVertex = NULL;
            VEdge2D*   currBoundaryEdge = NULL;
            VEdge2D*   nextBoundaryEdge = beginningEdgeOfIncidentFace;

            if( currIncidentFace == beginningEdgeOfIncidentFace->getLeftFace() ) {
                firstBoundaryVertex = beginningEdgeOfIncidentFace->getEndVertex();
            }
            else {
                firstBoundaryVertex = beginningEdgeOfIncidentFace->getStartVertex();
            }


            if( firstBoundaryVertex->getStatus() == RED_V ) {
                isCurrVVertexRed = true;
            }
            isPrevVVertexRed = isCurrVVertexRed;
            nextBoundaryEdge = getCCWNextEdgeOnFace( beginningEdgeOfIncidentFace, currIncidentFace);

            do 
            {
                currBoundaryEdge = nextBoundaryEdge;
                nextBoundaryEdge = getCCWNextEdgeOnFace( currBoundaryEdge, currIncidentFace );

                if( currIncidentFace == currBoundaryEdge->getLeftFace() ) {
                    currBoundaryVertex = currBoundaryEdge->getEndVertex();
                }
                else {
                    currBoundaryVertex = currBoundaryEdge->getStartVertex();
                }


                if( currBoundaryVertex->getStatus() == RED_V ) {
                    isCurrVVertexRed = true;
                }
                else {
                    isCurrVVertexRed = false;
                }


                if( isCurrVVertexRed != isPrevVVertexRed ) {
                    numComponentsOfBoundaryVVertices++;
                }

                isPrevVVertexRed = isCurrVVertexRed;


            } while ( nextBoundaryEdge != beginningEdgeOfIncidentFace );


            if( firstBoundaryVertex->getStatus() == RED_V ) {
                isCurrVVertexRed = true;
            }
            else {
                isCurrVVertexRed = false;
            }


            if( isCurrVVertexRed != isPrevVVertexRed ) {
                numComponentsOfBoundaryVVertices++;
            }
            

            candidateForRedVVertex->setStatus( WHITE_V );

            if( numComponentsOfBoundaryVVertices == 2 || numComponentsOfBoundaryVVertices == 4 || numComponentsOfBoundaryVVertices == 6 ) {
                continue;
            }
            else {
                return true;
            }

        } while ( nextIncidentEdge != beginningEdgeOfVertex );
    }

    else {
        // would-be red VVertex의 incident한 cell 3개를 가지고 온다
        list<VFace2D*> incidentVFaces;
        candidateForRedVVertex->getIncident3VFaces( incidentVFaces );


        // 각 cell마다 split이 일어났는지 확인한다.
        for( list<VFace2D*>::iterator it_face = incidentVFaces.begin() ; it_face != incidentVFaces.end() ; it_face++ ) {
            VFace2D* currVFace = *it_face;

            // would-be red VVertex를 임시로 red로 만든다
            candidateForRedVVertex->setStatus( RED_V );

            list<VVertex2D*> boundaryVVerticesOfCurrVFace;
            currVFace->getBoundaryVVertices( boundaryVVerticesOfCurrVFace );

            list<VVertex2D*>::iterator it_boundaryVVertex = boundaryVVerticesOfCurrVFace.begin();
            int     numComponentsOfBoundaryVVertices = 0;
            bool    isCurrVVertexRed = false;

            //if( (*it_boundaryVVertex)->getStatus() == RED_V || (*it_boundaryVVertex)->getStatus() == YELLOW_V ) {
            if( (*it_boundaryVVertex)->getStatus() == RED_V ) {
                isCurrVVertexRed = true;
            }
            bool    isPrevVVertexRed = isCurrVVertexRed;
            it_boundaryVVertex++;


            for( ; it_boundaryVVertex != boundaryVVerticesOfCurrVFace.end() ; it_boundaryVVertex++ ) {
                VVertex2D* currBoundaryVVertex = *it_boundaryVVertex;
                if( currBoundaryVVertex->getStatus() == RED_V ) {
                    isCurrVVertexRed = true;
                }
                else {
                    isCurrVVertexRed = false;
                }

                if( isCurrVVertexRed != isPrevVVertexRed ) {
                    numComponentsOfBoundaryVVertices++;
                }

                isPrevVVertexRed = isCurrVVertexRed;
            }

            VVertex2D* firstVVertex =   *boundaryVVerticesOfCurrVFace.begin();
            VVertex2D* lastVVertex =    *boundaryVVerticesOfCurrVFace.rbegin();

            if( firstVVertex->getStatus() == RED_V ) {
                isCurrVVertexRed = true;
            }
            else {
                isCurrVVertexRed = false;
            }

            if( isCurrVVertexRed != isPrevVVertexRed ) {
                numComponentsOfBoundaryVVertices++;
            }

            candidateForRedVVertex->setStatus( WHITE_V );

            //if( numComponentsOfBoundaryVVertices == 2 || numComponentsOfBoundaryVVertices == 4 ) {
            if( numComponentsOfBoundaryVVertices == 2 || numComponentsOfBoundaryVVertices == 4 || numComponentsOfBoundaryVVertices == 6 ) {
                continue;
            }
            else {
                return true;
            }
        }
    }

    

    return false;

}



bool VoronoiDiagram2DC::isAnomalizingEdge( VEdge2D* const incidentEdge, Generator2D* const newGenerator )
{
    if( incidentEdge->getStartVertex()->isFictitious() || incidentEdge->getEndVertex()->isFictitious() ) {
        return false;
    }

    rg_Circle2D disk[3];
    
    disk[0] = *newGenerator->getDisk();
    disk[1] = *incidentEdge->getLeftFace()->getGenerator()->getDisk();
    disk[2] = *incidentEdge->getRightFace()->getGenerator()->getDisk();

    rg_Circle2D firstCircle;
    rg_Circle2D secondCircle;

    int numCircumcircles = 0;
    numCircumcircles = rg_Circle2D::makeCircumcircle( disk[0], disk[1], disk[2], firstCircle, secondCircle );


    rg_Point2D center0;
    rg_Point2D CWBoundary;
    rg_Point2D CCWBoundary;

    rg_Point2D firstLocation;
    rg_Point2D secondLocation;

    rg_Point3D CWBoundary3D;
    rg_Point3D CCWBoundary3D;
    rg_Point3D firstLocation3D;
    rg_Point3D secondLocation3D;

    rg_Point3D CW_First_Test;
    rg_Point3D CCW_First_Test;
    rg_Point3D CW_Second_Test;
    rg_Point3D CCW_Second_Test;

    switch (numCircumcircles) {

    case 0:

        return false;



    case 1:

        return false;



    case 2:



        if( incidentEdge->getLeftFace()->getGenerator()->getDisk()->getRadius() < incidentEdge->getRightFace()->getGenerator()->getDisk()->getRadius() ) {
            center0     = incidentEdge->getRightFace()->getGenerator()->getDisk()->getCenterPt();
            CWBoundary  = incidentEdge->getEndVertex()->getLocation() - center0;
            CCWBoundary = incidentEdge->getStartVertex()->getLocation() - center0;
        }
        else {
            center0     = incidentEdge->getLeftFace()->getGenerator()->getDisk()->getCenterPt() ;
            CWBoundary  = incidentEdge->getStartVertex()->getLocation() - center0;
            CCWBoundary = incidentEdge->getEndVertex()->getLocation() - center0;
        }

        firstLocation = firstCircle.getCenterPt();
        secondLocation = secondCircle.getCenterPt();
        firstLocation   = firstLocation - center0;
        secondLocation  = secondLocation - center0;

        CWBoundary3D     = rg_Point3D( CWBoundary );
        CCWBoundary3D    = rg_Point3D( CCWBoundary );
        firstLocation3D  = rg_Point3D( firstLocation );
        secondLocation3D = rg_Point3D( secondLocation );
        
        CW_First_Test = CWBoundary3D.crossProduct( firstLocation3D );
        CCW_First_Test = CCWBoundary3D.crossProduct( firstLocation3D );
        CW_Second_Test = CWBoundary3D.crossProduct( secondLocation3D );
        CCW_Second_Test = CCWBoundary3D.crossProduct( secondLocation3D );

        if( CW_First_Test.getZ() < 0 ) {
            return false;
        }
        if( CW_Second_Test.getZ() < 0 ) {
            return false;
        }
        if( CCW_First_Test.getZ() > 0 ) {
            return false;
        }
        if( CCW_Second_Test.getZ() > 0 ) {
            return false;
        }

        return true;


    default:
        return false;
    }


}



void VoronoiDiagram2DC::mergeSplitVEdgesByFictitiousVVertex( list<VVertex2D*>& fictitousVVertices )
{
    for( list<VVertex2D*>::iterator it_fictitousVVertex = fictitousVVertices.begin() ; it_fictitousVVertex != fictitousVVertices.end() ; it_fictitousVVertex++ ) {
        VVertex2D* currFictitiousVVertex = *it_fictitousVVertex;

        VEdge2D* edgeBeforeSplitVertex = currFictitiousVVertex->getFirstVEdge();
        VEdge2D* edgeAfterSplitVertex = edgeBeforeSplitVertex->getRightHand();

        if( edgeAfterSplitVertex->getLeftHand()->getLeftLeg() == edgeAfterSplitVertex ) {
            edgeAfterSplitVertex->getLeftHand()->setLeftLeg( edgeBeforeSplitVertex );
        }
        else {
            edgeAfterSplitVertex->getLeftHand()->setRightHand( edgeBeforeSplitVertex );
        }

        if( edgeAfterSplitVertex->getRightHand()->getRightLeg() == edgeAfterSplitVertex ) {
            edgeAfterSplitVertex->getRightHand()->setRightLeg( edgeBeforeSplitVertex );
        }
        else {
            edgeAfterSplitVertex->getRightHand()->setLeftHand( edgeBeforeSplitVertex );
        }

        edgeBeforeSplitVertex->setEndVertex( edgeAfterSplitVertex->getEndVertex() );
        edgeBeforeSplitVertex->setLeftHand( edgeAfterSplitVertex->getLeftHand() );
        edgeBeforeSplitVertex->setRightHand( edgeAfterSplitVertex->getRightHand() );

        edgeAfterSplitVertex->getEndVertex()->setFirstVEdge( edgeBeforeSplitVertex );

        edgeAfterSplitVertex->setStatus( RED_E );        
    }
}



rg_BOOL VoronoiDiagram2DC::areTwoVerticesDefinedBySameGeneratorTriplet( VVertex2D* const firstVertex, VVertex2D* const secondVertex )
{
    list<Generator2D*> generatorsOfFirstVertex;
    firstVertex->getDefining3Generators( generatorsOfFirstVertex );

    set<Generator2D*> setOfGeneratorsOfFirstVertex;
    setOfGeneratorsOfFirstVertex.insert( generatorsOfFirstVertex.begin(), generatorsOfFirstVertex.end() );


    list<Generator2D*> generatorsOfSecondVertex;
    secondVertex->getDefining3Generators( generatorsOfSecondVertex );

    rg_BOOL areTwoVerticesDefinedBySameGenerators = rg_TRUE;
    for ( list<Generator2D*>::iterator i_gen=generatorsOfSecondVertex.begin(); i_gen!=generatorsOfSecondVertex.end(); i_gen++ ) {
        Generator2D* currGenerator = *i_gen;

        if ( setOfGeneratorsOfFirstVertex.find(currGenerator) == setOfGeneratorsOfFirstVertex.end() ) {
            areTwoVerticesDefinedBySameGenerators = rg_FALSE;
            break;
        }
    }

    return areTwoVerticesDefinedBySameGenerators;

//     list<Generator2D*> definingCircles1;
//     list<Generator2D*> definingCircles2;
// 
//     firstVertex->getDefining3Generators( definingCircles1 );
//     secondVertex->getDefining3Generators( definingCircles2 );
// 
//     list<Generator2D*>::iterator it_defining1;
//     list<Generator2D*>::iterator it_defining2;
//     rg_BOOL globalRet = rg_FALSE;
// 
//     for( it_defining1 = definingCircles1.begin() ; it_defining1 != definingCircles1.end() ; it_defining1++ ) {
//         Generator2D* currGenerator1 = *it_defining1;
//         rg_BOOL localRet = rg_FALSE;
// 
//         for( it_defining2 = definingCircles2.begin() ; it_defining2 != definingCircles2.end() ; it_defining2++ ) {
//             Generator2D* currGenerator2 = *it_defining2;
// 
//             if( currGenerator1 == currGenerator2 ) {
//                 definingCircles2.erase( it_defining2 );
//                 localRet = rg_TRUE;
//                 break;
//             }
//         }
// 
//         if( localRet == rg_FALSE ) {
//             break;
//         }
//     }
// 
//     if( definingCircles2.size() == 0 ) {
//         globalRet = rg_TRUE;
//     }
// 
//     return globalRet;
}



bool VoronoiDiagram2DC::compareCircleDescendingorder( const rg_Circle2D& circle1, const rg_Circle2D& circle2 )
{
    return circle1.getRadius() > circle2.getRadius();
}



void VoronoiDiagram2DC::getPhantomGenerator( list<const Generator2D*>& phantomGeneratorsList ) const
{
    for( list<Generator2D>::const_iterator it_generator = m_phantomGenerators.begin() ; it_generator != m_phantomGenerators.end() ; it_generator++ ) {
        const Generator2D* currGenerator = &(*it_generator);
        phantomGeneratorsList.push_back( currGenerator );
    }
}

VVertex2D* VoronoiDiagram2DC::splitVEdgeAtFictitiousVVertex( VEdge2D* const currVEdge )
{
    VVertex2D* fictitiousVVertex = createVVertex( VVertex2D(m_VVertices.rbegin()->getID()+1) );
    fictitiousVVertex->setStatus( YELLOW_V );
    fictitiousVVertex->setFirstVEdge( currVEdge );
    fictitiousVVertex->setFictitious();

    VEdge2D* dividedVEdge = createVEdge( VEdge2D(m_VEdges.rbegin()->getID() + 1) );

    if( currVEdge->getLeftHand()->getLeftLeg() == currVEdge ) {
        currVEdge->getLeftHand()->setLeftLeg( dividedVEdge );
    }
    else {
        currVEdge->getLeftHand()->setRightHand( dividedVEdge );
    }

    if( currVEdge->getRightHand()->getRightLeg() == currVEdge ) {
        currVEdge->getRightHand()->setRightLeg( dividedVEdge );
    }
    else {
        currVEdge->getRightHand()->setLeftHand( dividedVEdge );
    }

    dividedVEdge->setLeftHand(      currVEdge->getLeftHand() );
    dividedVEdge->setRightHand(     currVEdge->getRightHand() );
    dividedVEdge->setLeftLeg(       currVEdge );
    dividedVEdge->setRightLeg(      currVEdge );
    dividedVEdge->setStartVertex(   fictitiousVVertex );
    dividedVEdge->setEndVertex(     currVEdge->getEndVertex() );
    dividedVEdge->setLeftFace(      currVEdge->getLeftFace() );
    dividedVEdge->setRightFace(     currVEdge->getRightFace() );

    currVEdge->setEndVertex(        fictitiousVVertex );
    currVEdge->setLeftHand(         dividedVEdge );
    currVEdge->setRightHand(        dividedVEdge );

    dividedVEdge->getEndVertex()->setFirstVEdge( dividedVEdge );

    return fictitiousVVertex;
}

void VoronoiDiagram2DC::updateEdge( VEdge2D* target )
{
    rg_Point2D sp, ep, tvs, tve, passPoint;

    if( rg_EQ( target->getLeftFace()->getGenerator()->getDisk()->getRadius(),
        target->getRightFace()->getGenerator()->getDisk()->getRadius(), resNeg6 ) )
    {
        target->setGeometryLine();	//updated to a line segment
    }
    else
    {
        sp = target->getStartVertex()->getLocation();
        ep = target->getEndVertex()->getLocation();
        tvs = getTangentVector( target->getStartVertex(), 
            target->getLeftFace(),
            target->getRightFace() );
        tve = getTangentVector( target->getEndVertex(), 
            target->getLeftFace(),
            target->getRightFace() );
        passPoint = getPassingPtOnEdge( target->getLeftFace(), 
            target->getRightFace(),
            target->getStartVertex(),
            target->getEndVertex());

        target->setGeometry(sp, tvs, ep, tve, passPoint); 	
        if( rg_EQ(target->getGeometry().getWeight(1) , 1.) )
        {
            target->setGeometryLine();
        }
    }
}

rg_Point2D VoronoiDiagram2DC::getTangentVector( VVertex2D* tVertex, VFace2D* lFace, VFace2D* rFace )
{
    rg_Point2D sp, vec1, vec2;
    sp = tVertex->getLocation();
    vec1 = (lFace->getGenerator()->getDisk()->getCenterPt() - sp).getUnitVector();
    vec2 = (rFace->getGenerator()->getDisk()->getCenterPt() - sp).getUnitVector();

    if( rg_ZERO( (vec1+vec2).magnitude() ) )
    {
        return rg_Point2D( -vec1.getY(), vec1.getX() );
    }
    else
    {
        return vec1 + vec2;
    }
}

rg_Point2D VoronoiDiagram2DC::getPassingPtOnEdge( VFace2D* lFace, VFace2D* rFace, VVertex2D* sVertex, VVertex2D* eVertex )
{
    rg_Circle2D c1, c2;	//c1.r > c2.r
    if( lFace->getGenerator()->getDisk()->getRadius() > rFace->getGenerator()->getDisk()->getRadius() )
    {
        c1 = *lFace->getGenerator()->getDisk();
        c2 = *rFace->getGenerator()->getDisk();
    }
    else
    {
        c2 = *lFace->getGenerator()->getDisk();
        c1 = *rFace->getGenerator()->getDisk();
    }

    rg_Point2D cp = ( c1.getCenterPt() + c2.getCenterPt() ) / 2.;
    rg_Point2D vector = (cp - c2.getCenterPt()).getUnitVector();

    rg_REAL A = c2.getCenterPt().getX() - c1.getCenterPt().getX();
    rg_REAL B = c2.getCenterPt().getY() - c1.getCenterPt().getY();
    rg_REAL r = c1.getRadius() - c2.getRadius();

    rg_REAL t = (r*r - A*A - B*B)/(A*vector.getX()+B*vector.getY() - r)/2.;

    //Is t always positive?
    // --> t should be nonnegative.
    rg_Point2D passingPt = c2.getCenterPt() + t * vector;

    return passingPt;
}



void VoronoiDiagram2DC::updateGeometry()
{
    list<VEdge2D>::iterator it_edges;
    for( it_edges = m_VEdges.begin() ; it_edges != m_VEdges.end() ; it_edges++ ) {
        VEdge2D* currEdge = &(*it_edges);

        if( currEdge->isInfinite() ) {
            continue;
        }

        updateEdge(currEdge);
    }
}



void VoronoiDiagram2DC::wavePropagation_ver1( Generator2D* newGenerator, list<VVertex2D*>& redVVertices, list<VVertex2D*>& blueVVertices, list<VVertex2D*>& fictitiousVVertices )
{
    list<VEdge2D*>   anomalyTestedVEdges;
    list<VVertex2D*>::iterator it_redVVertex = redVVertices.begin();

    for( ; it_redVVertex != redVVertices.end() ; it_redVVertex++ ) {
        VVertex2D* currRedVVertex = *it_redVVertex;

        findAnomalizingEdgeAmongIncidentVEdgesAndSplit( currRedVVertex, newGenerator, fictitiousVVertices, anomalyTestedVEdges );

        findRedVVertexAmongNeighbor( currRedVVertex, newGenerator, blueVVertices, redVVertices, fictitiousVVertices, anomalyTestedVEdges );
    }

    resetAnomalyTestDoneToggle( anomalyTestedVEdges );
}



void VoronoiDiagram2DC::wavePropagation_ver2( Generator2D* newGenerator, list<VVertex2D*>& redVVertices, list<VVertex2D*>& blueVVertices, list<VVertex2D*>& fictitiousVVertices )
{
    list<VEdge2D*>   anomalyTestedVEdges;
    list<VVertex2D*>::iterator it_redVVertex = redVVertices.begin();

    for( ; it_redVVertex != redVVertices.end() ; it_redVVertex++ ) {
        VVertex2D* currRedVVertex = *it_redVVertex;

        list<VEdge2D*> incidentVEdges;
        currRedVVertex->getIncident3VEdges( incidentVEdges );

        for( list<VEdge2D*>::iterator it_incidentVEdge = incidentVEdges.begin() ; it_incidentVEdge != incidentVEdges.end() ; it_incidentVEdge++ ) {
            VEdge2D* currVEdge = *it_incidentVEdge;

            if( isAnomalizingEdge( currVEdge, newGenerator ) ) {
                VVertex2D* fictitiousVVertex = splitVEdgeAtFictitiousVVertex( currVEdge );
                fictitiousVVertices.push_back( fictitiousVVertex );
                continue;
            }
            else {
                VVertex2D* couldBeRedVertex = currVEdge->getOppositVVertex( currRedVVertex );
                if( couldBeRedVertex->getStatus() != WHITE_V ) {
                    continue;
                }

                double MUValue = couldBeRedVertex->computeMUValue( newGenerator );

                if( MUValue >= 0 ) {
                    couldBeRedVertex->setStatus( BLUE_V );
                    blueVVertices.push_back( couldBeRedVertex );
                    continue;
                }

                if( hasCycleOccurred( couldBeRedVertex ) ) {
                    findAnomalizingEdgeAmongIncidentVEdgesAndSplit( couldBeRedVertex, newGenerator, fictitiousVVertices, anomalyTestedVEdges );
                    if( hasCycleOccurred( couldBeRedVertex ) ) {
                        couldBeRedVertex->setStatus( BLUE_V );
                        blueVVertices.push_back( couldBeRedVertex );
                        continue;
                    }
                }

                if( hasOldRegionSplit( couldBeRedVertex ) ) { 
                    couldBeRedVertex->setStatus( BLUE_V );
                    blueVVertices.push_back( couldBeRedVertex );
                    continue;
                }

                couldBeRedVertex->setStatus( RED_V );
                redVVertices.push_back( couldBeRedVertex );
            }
        }
    }

    resetAnomalyTestDoneToggle( anomalyTestedVEdges );
}

void VoronoiDiagram2DC::removePhantomGenerators()
{
    VFace2D*            virtualFace = &(*m_VFaces.begin());
    set<VFace2D*>       phantomFaces; 

    list<VEdge2D*> virtualEdges;
    list<VEdge2D*> unboundedEdges_before_phantom_removal;
    list<VEdge2D*> edgesDefinedByPhantomNInputDisk;

    collectEntitiesInfluencedByPhantomRemoval( phantomFaces, virtualEdges, unboundedEdges_before_phantom_removal, edgesDefinedByPhantomNInputDisk );

    putRedTagsOnEntitiesToBeRemoved( virtualEdges, unboundedEdges_before_phantom_removal );
    
    adjustTopologyOfInterfaceBetweenInputDisksNPhantomDisks( virtualFace, unboundedEdges_before_phantom_removal, phantomFaces, edgesDefinedByPhantomNInputDisk );
    
    removePhantomFaces( phantomFaces );
    
    makeCorrectTopologyOfTheInterfaceByFlippingEdges( virtualFace );
    
    assignCoordinateOfInfiniteVertices();
}



void VoronoiDiagram2DC::collectEntitiesInfluencedByPhantomRemoval( set<VFace2D*>& phantomFaces, list<VEdge2D*>& virtualEdges, list<VEdge2D*>& unboundedEdges_before_phantom_removal, list<VEdge2D*>& edgesDefinedByPhantomNInputDisk )
{
    const int numDummyCircles = 3;
    list<VFace2D>::iterator it_face = m_VFaces.begin();
    VFace2D* virtualFace = &(*it_face);
    it_face++;

    for( int i = 0 ; i < numDummyCircles ; i++, it_face++ ) {
        VFace2D* currPhantomRegion = &(*it_face);
        phantomFaces.insert( currPhantomRegion );
    }



    virtualFace->getBoundaryVEdges( virtualEdges );
    list<VEdge2D*>::iterator i_virtualEdge = virtualEdges.begin();
    for( ; i_virtualEdge != virtualEdges.end() ; i_virtualEdge++ ) {
        VEdge2D* currVirtualEdge = *i_virtualEdge;
        if( currVirtualEdge->getRightFace() == virtualFace )
        {
            unboundedEdges_before_phantom_removal.push_back( currVirtualEdge->getLeftHand() );
        }
        else
        {
            unboundedEdges_before_phantom_removal.push_back( currVirtualEdge->getRightLeg() );
        }
    }



    set<VFace2D*>::iterator NOT_PHANTOM = phantomFaces.end();
    set<VFace2D*>::iterator i_phantom = phantomFaces.begin();
    for( ; i_phantom != phantomFaces.end() ; i_phantom++ ) {
        list<VEdge2D*> boundaryEdgeOfPhantomRegion;
        (*i_phantom)->getBoundaryVEdges( boundaryEdgeOfPhantomRegion );

        list<VEdge2D*>::iterator i_bndEdge = boundaryEdgeOfPhantomRegion.begin();
        for( ; i_bndEdge != boundaryEdgeOfPhantomRegion.end() ; i_bndEdge++ ) {
            VEdge2D* currEdge = *i_bndEdge;
            if( currEdge->getLeftFace()->getGenerator() == NULL || currEdge->getRightFace()->getGenerator() == NULL )
            {
                continue;
            }
            else if( phantomFaces.find(currEdge->getLeftFace()) != NOT_PHANTOM && phantomFaces.find(currEdge->getRightFace()) != NOT_PHANTOM )
            {
                continue;
            }

            edgesDefinedByPhantomNInputDisk.push_back( currEdge );
        }
    }



//     list<VEdge2D>::iterator it_edge = m_VEdges.begin();
//     for( ; it_edge != m_VEdges.end() ; it_edge++ ) {
//         VEdge2D* currEdge = &(*it_edge);
//         if( currEdge->getStatus() == RED_E ) {
//             continue;
//         }
// 
//         if( currEdge->isInfinite() ) {
//             virtualEdges.push_back( currEdge );
//         }
//         else if( currEdge->isUnBounded() ) {
//             unboundedEdges_before_phantom_removal.push_back( currEdge );
//         }
//         else if( (currEdge->getLeftFace()->isUnBounded() && currEdge->getRightFace()->isBounded()) ||
//             (currEdge->getLeftFace()->isBounded() && currEdge->getRightFace()->isUnBounded()) ) {
//                 edgesDefinedByPhantomNInputDisk.push_back( currEdge );
//         }
//     }
}



void VoronoiDiagram2DC::putRedTagsOnEntitiesToBeRemoved( list<VEdge2D*>& virtualEdges, list<VEdge2D*>& unboundedEdges_before_phantom_removal )
{
    list<VEdge2D*>::iterator it_infiniteEdge = virtualEdges.begin();
    for( ; it_infiniteEdge != virtualEdges.end() ; it_infiniteEdge++ ) {
        VEdge2D* currInfiniteEdge = *it_infiniteEdge;
        currInfiniteEdge->setStatus( RED_E );
    }

    list<VEdge2D*>::iterator it_unboundedEdge = unboundedEdges_before_phantom_removal.begin();
    for( ; it_unboundedEdge != unboundedEdges_before_phantom_removal.end() ; it_unboundedEdge++ ) {
        VEdge2D* currUnboundedEdge = *it_unboundedEdge;
        currUnboundedEdge->setStatus( RED_E );
        if( currUnboundedEdge->getStartVertex()->isInfinite() ) {
            currUnboundedEdge->getStartVertex()->setStatus( RED_V );
        }
        else {
            currUnboundedEdge->getEndVertex()->setStatus( RED_V );
        }
    }
}



void VoronoiDiagram2DC::adjustTopologyOfInterfaceBetweenInputDisksNPhantomDisks( VFace2D* virtualFace, list<VEdge2D*>& unboundedEdges_before_phantom_removal, set<VFace2D*>& phantomFaces, list<VEdge2D*>& edgesDefinedByPhantomNInputDisk )
{
    list<VEdge2D*>::iterator it_unboundedEdge = unboundedEdges_before_phantom_removal.begin();
    for( ; it_unboundedEdge != unboundedEdges_before_phantom_removal.end() ; it_unboundedEdge++ ) {
        VEdge2D* currUnboundedEdge = *it_unboundedEdge;
        VEdge2D* mergedEdge = createVEdge( VEdge2D(m_VEdges.rbegin()->getID() + 1) );

        if( currUnboundedEdge->getStartVertex()->getStatus() == RED_V ) {
            currUnboundedEdge->getEndVertex()->setStatus( RED_V );
            VEdge2D* leftHand = currUnboundedEdge->getLeftHand();
            VEdge2D* rightHand = currUnboundedEdge->getRightHand();
            leftHand->setStatus( RED_E );
            rightHand->setStatus( RED_E );


            if( currUnboundedEdge == leftHand->getRightHand() ) {
                mergedEdge->setStartVertex( leftHand->getStartVertex() );
                mergedEdge->setRightLeg( leftHand->getRightLeg() );
                mergedEdge->setLeftLeg( leftHand->getLeftLeg() );
                mergedEdge->setRightFace( virtualFace );
                mergedEdge->setLeftFace( leftHand->getLeftFace() );

                mergedEdge->getStartVertex()->setFirstVEdge( mergedEdge );
                mergedEdge->getLeftFace()->setFirstVEdge( mergedEdge );

                if( leftHand == leftHand->getRightLeg()->getRightHand() ) {
                    leftHand->getRightLeg()->setRightHand( mergedEdge );
                }
                else {
                    leftHand->getRightLeg()->setLeftLeg( mergedEdge );
                }

                if( leftHand == leftHand->getLeftLeg()->getLeftHand() ) {
                    leftHand->getLeftLeg()->setLeftHand( mergedEdge );
                }
                else {
                    leftHand->getLeftLeg()->setRightLeg( mergedEdge );
                }
            }
            else {
                mergedEdge->setStartVertex( leftHand->getEndVertex() );
                mergedEdge->setRightLeg( leftHand->getLeftHand() );
                mergedEdge->setLeftLeg( leftHand->getRightHand() );
                mergedEdge->setRightFace( virtualFace );
                mergedEdge->setLeftFace( leftHand->getRightFace() );

                mergedEdge->getStartVertex()->setFirstVEdge( mergedEdge );
                mergedEdge->getLeftFace()->setFirstVEdge( mergedEdge );

                if( leftHand == leftHand->getLeftHand()->getRightHand() ) {
                    leftHand->getLeftHand()->setRightHand( mergedEdge );
                }
                else {
                    leftHand->getLeftHand()->setLeftLeg( mergedEdge );
                }

                if( leftHand == leftHand->getRightHand()->getLeftHand() ) {
                    leftHand->getRightHand()->setLeftHand( mergedEdge );
                }
                else {
                    leftHand->getRightHand()->setRightLeg( mergedEdge );
                }
            }


            if( currUnboundedEdge == rightHand->getLeftHand() ) {
                mergedEdge->setEndVertex( rightHand->getStartVertex() );
                mergedEdge->setRightHand( rightHand->getLeftLeg() );
                mergedEdge->setLeftHand( rightHand->getRightLeg() );

                mergedEdge->getEndVertex()->setFirstVEdge( mergedEdge );

                if( rightHand == rightHand->getLeftLeg()->getLeftHand() ) {
                    rightHand->getLeftLeg()->setLeftHand( mergedEdge );
                }
                else {
                    rightHand->getLeftLeg()->setRightLeg( mergedEdge );
                }

                if( rightHand == rightHand->getRightLeg()->getRightHand() ) {
                    rightHand->getRightLeg()->setRightHand( mergedEdge );
                }
                else {
                    rightHand->getRightLeg()->setLeftLeg( mergedEdge );
                }
            }
            else {
                mergedEdge->setEndVertex( rightHand->getEndVertex() );
                mergedEdge->setRightHand( rightHand->getRightHand() );
                mergedEdge->setLeftHand( rightHand->getLeftHand() );

                mergedEdge->getEndVertex()->setFirstVEdge( mergedEdge );

                if( rightHand == rightHand->getRightHand()->getLeftHand() ) {
                    rightHand->getRightHand()->setLeftHand( mergedEdge );
                }
                else {
                    rightHand->getRightHand()->setRightLeg( mergedEdge );
                }

                if( rightHand == rightHand->getLeftHand()->getRightHand() ) {
                    rightHand->getLeftHand()->setRightHand( mergedEdge );
                }
                else {
                    rightHand->getLeftHand()->setLeftLeg( mergedEdge );
                }
            }
        }
        else {
            currUnboundedEdge->getStartVertex()->setStatus( RED_V );
            VEdge2D* rightLeg = currUnboundedEdge->getRightLeg();
            VEdge2D* leftLeg = currUnboundedEdge->getLeftLeg();
            rightLeg->setStatus( RED_E );
            leftLeg->setStatus( RED_E );


            if( currUnboundedEdge == rightLeg->getRightHand() ) {
                mergedEdge->setStartVertex( rightLeg->getStartVertex() );
                mergedEdge->setRightLeg( rightLeg->getRightLeg() );
                mergedEdge->setLeftLeg( rightLeg->getLeftLeg() );
                mergedEdge->setRightFace( virtualFace );
                mergedEdge->setLeftFace( rightLeg->getLeftFace() );

                mergedEdge->getStartVertex()->setFirstVEdge( mergedEdge );
                mergedEdge->getLeftFace()->setFirstVEdge( mergedEdge );

                if( rightLeg == rightLeg->getRightLeg()->getRightHand() ) {
                    rightLeg->getRightLeg()->setRightHand( mergedEdge );
                }
                else {
                    rightLeg->getRightLeg()->setLeftLeg( mergedEdge );
                }

                if( rightLeg == rightLeg->getLeftLeg()->getLeftHand() ) {
                    rightLeg->getLeftLeg()->setLeftHand( mergedEdge );
                }
                else {
                    rightLeg->getLeftLeg()->setRightLeg( mergedEdge );
                }
            }
            else {
                mergedEdge->setStartVertex( rightLeg->getEndVertex() );
                mergedEdge->setRightLeg( rightLeg->getLeftHand() );
                mergedEdge->setLeftLeg( rightLeg->getRightHand() );
                mergedEdge->setRightFace( virtualFace );
                mergedEdge->setLeftFace( rightLeg->getRightFace() );

                mergedEdge->getStartVertex()->setFirstVEdge( mergedEdge );
                mergedEdge->getLeftFace()->setFirstVEdge( mergedEdge );

                if( rightLeg == rightLeg->getLeftHand()->getRightHand() ) {
                    rightLeg->getLeftHand()->setRightHand( mergedEdge );
                }
                else {
                    rightLeg->getLeftHand()->setLeftLeg( mergedEdge );
                }

                if( rightLeg == rightLeg->getRightHand()->getLeftHand() ) {
                    rightLeg->getRightHand()->setLeftHand( mergedEdge );
                }
                else {
                    rightLeg->getRightHand()->setRightLeg( mergedEdge );
                }
            }


            if( currUnboundedEdge == leftLeg->getLeftHand() ) {
                mergedEdge->setEndVertex( leftLeg->getStartVertex() );
                mergedEdge->setRightHand( leftLeg->getLeftLeg() );
                mergedEdge->setLeftHand( leftLeg->getRightLeg() );

                mergedEdge->getEndVertex()->setFirstVEdge( mergedEdge );

                if( leftLeg == leftLeg->getLeftLeg()->getLeftHand() ) {
                    leftLeg->getLeftLeg()->setLeftHand( mergedEdge );
                }
                else {
                    leftLeg->getLeftLeg()->setRightLeg( mergedEdge );
                }

                if( leftLeg == leftLeg->getRightLeg()->getRightHand() ) {
                    leftLeg->getRightLeg()->setRightHand( mergedEdge );
                }
                else {
                    leftLeg->getRightLeg()->setLeftLeg( mergedEdge );
                }
            }
            else {
                mergedEdge->setEndVertex( leftLeg->getEndVertex() );
                mergedEdge->setRightHand( leftLeg->getRightHand() );
                mergedEdge->setLeftHand( leftLeg->getLeftHand() );

                mergedEdge->getEndVertex()->setFirstVEdge( mergedEdge );

                if( leftLeg == leftLeg->getRightHand()->getLeftHand() ) {
                    leftLeg->getRightHand()->setLeftHand( mergedEdge );
                }
                else {
                    leftLeg->getRightHand()->setRightLeg( mergedEdge );
                }

                if( leftLeg == leftLeg->getLeftHand()->getRightHand() ) {
                    leftLeg->getLeftHand()->setRightHand( mergedEdge );
                }
                else {
                    leftLeg->getLeftHand()->setLeftLeg( mergedEdge );
                }
            }
        }

        virtualFace->setFirstVEdge( mergedEdge );
        edgesDefinedByPhantomNInputDisk.push_back( mergedEdge );
    }

    set<VFace2D*>::iterator NOT_PHANTOM = phantomFaces.end();

    list<VEdge2D*>::iterator it_outsideBoundingEdge = edgesDefinedByPhantomNInputDisk.begin();
    for( ; it_outsideBoundingEdge != edgesDefinedByPhantomNInputDisk.end() ; it_outsideBoundingEdge++ ) {
        VEdge2D* currOutsideBoundingEdge = *it_outsideBoundingEdge;

        if( currOutsideBoundingEdge->getStatus() == RED_E ) {
            continue;
        }
        
        if( phantomFaces.find(currOutsideBoundingEdge->getLeftFace()) != NOT_PHANTOM ) {
            currOutsideBoundingEdge->setLeftFace( virtualFace );
        }
        else {
            currOutsideBoundingEdge->setRightFace( virtualFace );
        }

        virtualFace->setFirstVEdge( currOutsideBoundingEdge );
    }
}



void VoronoiDiagram2DC::removePhantomFaces( set<VFace2D*>& phantomFaces )
{
    set<VFace2D*>::iterator it_phantomFace = phantomFaces.begin();
    for( ; it_phantomFace != phantomFaces.end() ; it_phantomFace++ ) {
        VFace2D phantomFace = *(*it_phantomFace);
        m_VFaces.remove( phantomFace );
    }
}



void VoronoiDiagram2DC::correctTopologyInUnboundedRegion()
{
    VEdge2D* currEdge = NULL;

    // fix faces without the first edge.
    list<VEdge2D>::iterator it_edge = m_VEdges.begin();
    for( ; it_edge != m_VEdges.end() ; it_edge++ ) {
        currEdge = &(*it_edge);

        VVertex2D* startVertex = currEdge->getStartVertex();
        if( startVertex != NULL && startVertex->getFirstVEdge() == NULL ) {
            startVertex->setFirstVEdge( currEdge );
        }

        VVertex2D* endVertex = currEdge->getEndVertex();
        if( endVertex != NULL && endVertex->getFirstVEdge() == NULL ) {
            endVertex->setFirstVEdge( currEdge );
        }

        VFace2D* leftFace = currEdge->getLeftFace();
        if( leftFace != NULL && leftFace->getFirstVEdge() == NULL ) {
            leftFace->setFirstVEdge( currEdge );
        }

        VFace2D* rightFace = currEdge->getRightFace();
        if( rightFace != NULL && rightFace->getFirstVEdge() == NULL ) {
            rightFace->setFirstVEdge( currEdge );
        }
    }
}

void VoronoiDiagram2DC::removeEdgesDefinedByOnlyPhantomGenerators()
{
    set<VFace2D*> phantomRegions;

    // find dummy disks
    const int numDummyCircles = 3;
    int dummy_i = 0;
    list<VFace2D>::iterator it_face = m_VFaces.begin();
    while( it_face != m_VFaces.end() && dummy_i < numDummyCircles ) {
        VFace2D* currFace = &(*it_face);

        if( currFace->isUnBounded() ) {
            phantomRegions.insert( currFace );
            dummy_i++;
        }

        it_face++;
    }

    // collect edges defined by phantom generators
    list<VEdge2D*> edgesDefinedByOnlyPhantom;

    set<VFace2D*>::iterator NOT_PHANTOM = phantomRegions.end();

    list<VEdge2D>::iterator it_edge = m_VEdges.begin();
    for( ; it_edge != m_VEdges.end() ; it_edge++ ) {
        VEdge2D* currEdge = &(*it_edge);

        VFace2D* leftFace = currEdge->getLeftFace();
        VFace2D* rightFace = currEdge->getRightFace();

        if(     phantomRegions.find( leftFace )  != NOT_PHANTOM
            &&  phantomRegions.find( rightFace ) != NOT_PHANTOM ) {
                edgesDefinedByOnlyPhantom.push_back( currEdge );
        }
    }

    // remove edges defined by phantom generators
    list<VVertex2D*> infiniteVertices;

    list<VEdge2D*>::iterator it_unboundedEdge = edgesDefinedByOnlyPhantom.begin();
    for( ; it_unboundedEdge != edgesDefinedByOnlyPhantom.end() ; it_unboundedEdge++ ) {
        VEdge2D*  unboundedEdge = *it_unboundedEdge;
        VVertex2D* infiniteVertex = NULL;

        if( unboundedEdge->getEndVertex()->isInfinite() ) {
            infiniteVertex = unboundedEdge->getEndVertex();

            VEdge2D* leftLeg = unboundedEdge->getLeftLeg();
            VEdge2D* rightLeg = unboundedEdge->getRightLeg();


            if( unboundedEdge == rightLeg->getRightHand() ) {
                rightLeg->setRightHand( leftLeg );
            }
            else {
                rightLeg->setLeftLeg( leftLeg );
            }


            if( unboundedEdge == leftLeg->getLeftHand() ) {
                leftLeg->setLeftHand( rightLeg );
            }
            else {
                leftLeg->setRightLeg( rightLeg );
            }


            if( unboundedEdge == unboundedEdge->getStartVertex()->getFirstVEdge() ) {
                unboundedEdge->getStartVertex()->setFirstVEdge( leftLeg );
            }
        }
        else {
            infiniteVertex = unboundedEdge->getStartVertex();

            VEdge2D* leftHand = unboundedEdge->getLeftHand();
            VEdge2D* rightHand = unboundedEdge->getRightHand();


            if( unboundedEdge == rightHand->getRightLeg() ) {
                rightHand->setRightLeg( leftHand );
            }
            else {
                rightHand->setLeftHand( leftHand );
            }


            if( unboundedEdge ==  leftHand->getRightHand() ) {
                leftHand->setRightHand( rightHand );
            }
            else {
                leftHand->setLeftLeg( rightHand );
            }


            if( unboundedEdge == unboundedEdge->getEndVertex()->getFirstVEdge() ) {
                unboundedEdge->getEndVertex()->setFirstVEdge( leftHand );
            }
        }

        m_VEdges.remove(*unboundedEdge);
        m_VVertices.remove(*infiniteVertex);
    }
}

void VoronoiDiagram2DC::replaceRegionsOfPhantomGeneratorsWithVirtualRegion()
{
    set<VFace2D*> phantomRegions;
    // find V-region of phantom generator
    const int numDummyCircles = 3;
    list<VFace2D>::iterator it_face = m_VFaces.begin();
    VFace2D* virtualFace = &(*it_face);
    it_face++;

    for( int i = 0 ; i < numDummyCircles ; i++, it_face++ ) {
        VFace2D* currPhantomRegion = &(*it_face);
        phantomRegions.insert( currPhantomRegion );
    }

    // collect edges defined by phantom generators
    set<VFace2D*>::iterator NOT_PHANTOM = phantomRegions.end();

    list<VEdge2D>::iterator it_edge = m_VEdges.begin();
    for( ; it_edge != m_VEdges.end() ; it_edge++ ) {
        VEdge2D* currEdge = &(*it_edge);

        if( currEdge->getStatus() == RED_E ) {
            continue;
        }


    }


}

void VoronoiDiagram2DC::removeVoronoiRegionsOfPhantomGenerators()
{

}



void VoronoiDiagram2DC::makeMoreThanTwoEdgesBetweenRealAndVirtualRegionIntoOne()
{

}



void VoronoiDiagram2DC::makeCorrectTopologyOfTheInterfaceByFlippingEdges( VFace2D* virtualFace )
{
    Priority_Q_VEdge possiblyFlippingEdges;
    collectPossiblyFlippingEdgesOnTheInterface_StoredInPriorityQueue( virtualFace, possiblyFlippingEdges );

    rg_dList<VEdge2D*> edgesToUpdateGeometry;

    while( possiblyFlippingEdges.size() > 0 ) {
        pair<VEdge2D*, rg_Circle2D> edgeToDefineMinTangentCircle = possiblyFlippingEdges.pop();
        edgeToDefineMinTangentCircle.first->setFalseCandidateForFlippingInPhantomRemoval();

        VEdge2D* edgeOfIncidentTriplet[2] = {NULL, NULL};
        flipThisEdgeNFindTwoIncidentTriplets( virtualFace, edgeToDefineMinTangentCircle, edgeOfIncidentTriplet[0], edgeOfIncidentTriplet[1] );
        
        reflectTheFlipOnTheIncidentTriplet( possiblyFlippingEdges, edgeOfIncidentTriplet[0]);
        reflectTheFlipOnTheIncidentTriplet( possiblyFlippingEdges, edgeOfIncidentTriplet[1]);
    }
}


void VoronoiDiagram2DC::assignCoordinateOfInfiniteVertices()
{
    rg_BoundingBox2D boundingBox( m_bucket.getMinPt(), m_bucket.getMaxPt() );

    list<VVertex2D>::iterator it_vertex = m_VVertices.begin();
    for( ; it_vertex != m_VVertices.end() ; it_vertex++ ) {
        VVertex2D* currVtx = &(*it_vertex);

        if( currVtx->getStatus() == RED_V ) {
            continue;
        }

        boundingBox.contain( currVtx->getLocation() );
    }

    rg_Point2D center = boundingBox.getCenterPt();
    double     radius = boundingBox.evaluateLongestLength();


//     rg_BoundingBox2D boundingBox;
//     VFace2D* virtualFace = &(*m_VFaces.begin());
// 
//     list<VFace2D>::iterator it_face = m_VFaces.begin();
//     for( ; it_face != m_VFaces.end() ; it_face++ ) {
//         VFace2D* currFace = &(*it_face);
// 
//         if( currFace == virtualFace ) {
//             continue;
//         }
// 
//         rg_Circle2D* circle = currFace->getGenerator()->getDisk();
//         double x = circle->getCenterPt().getX();
//         double y = circle->getCenterPt().getY();
//         double r = circle->getRadius();
// 
//         rg_BoundingBox2D box( rg_Point2D(x-r, y-r), rg_Point2D(x+r, y+r) );
//         boundingBox.contain( box );
//     }
// 
//     list<VVertex2D>::iterator it_vertex = m_VVertices.begin();
//     for( ; it_vertex != m_VVertices.end() ; it_vertex++ ) {
//         VVertex2D* currVtx = &(*it_vertex);
// 
//         if( currVtx->getStatus() == RED_V ) {
//             continue;
//         }
// 
//         rg_Point2D point = currVtx->getLocation();
//         boundingBox.contain( point );
//     }
// 
//     rg_Point2D center = boundingBox.getCenterPt();
//     double     radius = 1.5*boundingBox.evaluateLongestLength();


    // update geometry of edge going on infinite


    list<VEdge2D*> unboundedEdges;

    //collect unbounded edges
    list<VEdge2D>::iterator it_edge = m_VEdges.begin();
    for( ; it_edge != m_VEdges.end() ; it_edge++ ) {
        VEdge2D* currEdge = &(*it_edge);

        if( currEdge->getStatus() == RED_E ) {
            continue;
        }

        if( !currEdge->isUnBounded() ) {
            continue;
        }

        unboundedEdges.push_back( currEdge );
    }


    list<VEdge2D*>::iterator it_unboundedEdge = unboundedEdges.begin();
    for( ; it_unboundedEdge != unboundedEdges.end() ; it_unboundedEdge++ ) {
        VEdge2D* currEdge = *it_unboundedEdge;


        VVertex2D* vtxOnInfinite = NULL;
        VVertex2D* vtxOnFinite   = NULL;
        rg_Point2D ptOnFinite;

        // check which vertex is on infinite and finite
        if( currEdge->getStartVertex()->isInfinite() ) {
            vtxOnInfinite = currEdge->getStartVertex();
            vtxOnFinite   = currEdge->getEndVertex();
        }
        else {
            vtxOnInfinite = currEdge->getEndVertex();
            vtxOnFinite   = currEdge->getStartVertex();
        }

        ptOnFinite = vtxOnFinite->getLocation();

        rg_Circle2D leftCircle  = *currEdge->getLeftFace()->getGenerator()->getDisk();
        rg_Circle2D rightCircle = *currEdge->getRightFace()->getGenerator()->getDisk();

//         double givenRadius = leftCircle.getCenterPt().distance( rightCircle.getCenterPt() );
//         double distanceToCenterOfTangentCircle = ptOnFinite.distance( leftCircle.getCenterPt() ) - leftCircle.getRadius();
// 
//         if( distanceToCenterOfTangentCircle > givenRadius ) {
//             //givenRadius = distanceToCenterOfTangentCircle * 1.2;
//             givenRadius = distanceToCenterOfTangentCircle + 5.0;
//         }

        double givenRadius = radius - center.distance( ptOnFinite );

        rg_Circle2D tangentCircle[2];
        int numTangentCircle = computeCircleTangentTo2CirclesGivenRadius( leftCircle, rightCircle, givenRadius, tangentCircle );


        rg_Point2D vectorToRight = rightCircle.getCenterPt() - tangentCircle[0].getCenterPt();
        rg_Point2D vectorToLeft  = leftCircle.getCenterPt() - tangentCircle[0].getCenterPt();

        if( vtxOnInfinite == currEdge->getEndVertex() ) {
            if( vectorToRight * vectorToLeft < 0 ) {
                vtxOnInfinite->setLocation( tangentCircle[0].getCenterPt() );
            }
            else {
                vtxOnInfinite->setLocation( tangentCircle[1].getCenterPt() );
            }
        }
        else {
            if( vectorToRight * vectorToLeft < 0 ) {
                vtxOnInfinite->setLocation( tangentCircle[1].getCenterPt() );
            }
            else {
                vtxOnInfinite->setLocation( tangentCircle[0].getCenterPt() );
            }
        }
    }
}




void VoronoiDiagram2DC::collectPossiblyFlippingEdgesOnTheInterface_StoredInPriorityQueue( VFace2D* virtualFace, Priority_Q_VEdge& possiblyFlippingEdges )
{
    list<VEdge2D*> edgesToBoundVirtualFace;
    virtualFace->getBoundaryVEdges( edgesToBoundVirtualFace );

    list<VEdge2D*>::iterator it_boundingEdge = edgesToBoundVirtualFace.begin();
    for( ; it_boundingEdge != edgesToBoundVirtualFace.end() ; it_boundingEdge++ ) {
        VEdge2D* currEdge = *it_boundingEdge;

        rg_Circle2D circumcircle;
        bool isFlippingEdge = isFlippingEdgeOnVirtualRegion( currEdge, circumcircle );

        if( isFlippingEdge ) {
            possiblyFlippingEdges.push( pair<VEdge2D*, rg_Circle2D>(currEdge, circumcircle) );
            currEdge->setTrueCandidateForFlippingInPhantomRemoval();
        }
    }

}

VFace2D* VoronoiDiagram2DC::getMatingFace( VVertex2D* tVertex, VEdge2D* tEdge )
{
    if( tVertex == tEdge->getStartVertex() )
    {
        return getFaceOfOppositeSide( tEdge->getLeftFace(), 
            tEdge->getLeftHand() );

    }
    else
    {
        return getFaceOfOppositeSide( tEdge->getRightFace(), 
            tEdge->getRightLeg() );

    }
}



bool VoronoiDiagram2DC::isFlippingEdgeOnVirtualRegion( VEdge2D* edge, rg_Circle2D& circumcircle )
{
    VFace2D* virtualFace = &(*m_VFaces.begin());

    VFace2D* currFace = rg_NULL;
    VFace2D* nextFace = rg_NULL;
    VFace2D* prevFace = rg_NULL;
//     VEdge2D* prevEdge = rg_NULL;
//     VEdge2D* nextEdge = rg_NULL;
//     VVertex2D* prevStartVtx = rg_NULL;
//     VVertex2D* prevEndVtx   = rg_NULL;
//     VVertex2D* nextStartVtx = rg_NULL;
//     VVertex2D* nextEndVtx   = rg_NULL;

    if ( edge->getLeftFace() == virtualFace ) {
        currFace = edge->getRightFace();
        prevFace = getMatingFace( edge->getStartVertex(), edge );
        nextFace = getMatingFace( edge->getEndVertex(),   edge );
    }
    else {
        currFace = edge->getLeftFace();
        prevFace = getMatingFace( edge->getEndVertex(),   edge );
        nextFace = getMatingFace( edge->getStartVertex(), edge );
    }

//     if ( edge->getLeftFace() == virtualFace ) {
//         currFace = edge->getRightFace();
//         nextFace = getMatingFace( edge->getStartVertex(), edge );
//         prevFace = getMatingFace( edge->getEndVertex(),   edge );
//         prevEdge = edge->getRightLeg();
//         nextEdge = edge->getRightHand();
// 
//         prevEndVtx = edge->getStartVertex();
//         nextEndVtx = edge->getEndVertex();
// 
//         prevStartVtx = prevEdge->getOppositVVertex( prevEndVtx );
//         nextStartVtx = nextEdge->getOppositVVertex( nextEndVtx );
//     }
//     else {
//         currFace = edge->getLeftFace();
//         nextFace = getMatingFace( edge->getEndVertex(),   edge );
//         prevFace = getMatingFace( edge->getStartVertex(), edge );
//         prevEdge = edge->getLeftHand();
//         nextEdge = edge->getLeftLeg();
// 
//         prevEndVtx = edge->getEndVertex();
//         nextEndVtx = edge->getStartVertex();
// 
//         prevStartVtx = prevEdge->getOppositVVertex( prevEndVtx );
//         nextStartVtx = nextEdge->getOppositVVertex( nextEndVtx );
//     }

    rg_BOOL isFlippingEdge = rg_FALSE;

    if ( nextFace == virtualFace || prevFace == virtualFace ) {
        return isFlippingEdge;
    }
    if ( nextFace == prevFace || nextFace == currFace || prevFace == currFace ) {
        return isFlippingEdge;
    }


    rg_Circle2D generator[3] = { *prevFace->getGenerator()->getDisk(), *currFace->getGenerator()->getDisk(), *nextFace->getGenerator()->getDisk() };
    rg_Point2D  center[3]    = {generator[0].getCenterPt(), generator[1].getCenterPt(), generator[2].getCenterPt()};
    rg_Circle2D tangentCircle[2];
    rg_INT numCC = rg_Circle2D::makeCircumcircle( generator[0], generator[1], generator[2], tangentCircle[0], tangentCircle[1]);

    vector< rg_Circle2D > validCircumcircle;
    for ( rg_INT i=0; i<numCC; ++i ) {
        if( tangentCircle[i].getRadius() < INFINITE_CIRCLE_RADIUS ) {
            validCircumcircle.push_back( tangentCircle[i] );
        }
    }

    switch ( validCircumcircle.size() ) {
    case 0:
        // do nothing
        break;


    case 1:
        if( isCorrectCircumcircle( validCircumcircle[0], prevFace, currFace, nextFace ) ) {
            circumcircle = validCircumcircle[0];
            isFlippingEdge = true;
        }
        break;


    case 2:
        if( isCorrectCircumcircle( validCircumcircle[0], prevFace, currFace, nextFace ) ) {
            circumcircle = validCircumcircle[0];
            isFlippingEdge = true;
        }
        else if( isCorrectCircumcircle( validCircumcircle[1], prevFace, currFace, nextFace ) ) {
            circumcircle = validCircumcircle[1];
            isFlippingEdge = true;
        }
        else {
            int stop3 = 3;
        }
        break;

    default:
        break;
    }
//     bool isFirstTangetCircleTooBig = false;
//     bool isSecondTangetCircleTooBig = false;
// 
//     switch(numCC) {
//     case 1:
//         if(tangentCircle[0].getRadius() > 1e6) {
//             numCC = 0;
//         }
//         break;
//     case 2:
//         if(tangentCircle[0].getRadius() > 1e6) {
//             numCC = numCC - 1;
//             isFirstTangetCircleTooBig = true;
//         }
//         if(tangentCircle[1].getRadius() >1e6) {
//             numCC = numCC - 1;
//             isSecondTangetCircleTooBig = true;
//         }
//         if(numCC == 1){
//             if(isFirstTangetCircleTooBig) {
//                 tangentCircle[0] = tangentCircle[1];
//             }
//         }
//     }
// 
//     int stop0 = 3;
//     
//     for ( rg_INT i=0; i<numCC; i++ ) {
//         rg_Point2D centerOfTangentCircle = tangentCircle[i].getCenterPt();
// 
//         // vertex좌표 말고 vertex를 정의하는 triplet들이 tangent circle을 만든 combination과 같은지 판단
//         if ( centerOfTangentCircle == prevStartVtx->getLocation() || centerOfTangentCircle == nextStartVtx->getLocation() ) {
//             continue;
//         }
// 
//         //  prev vs curr
//         rg_Point2D vecCurrToPrev = (center[0] - center[1]).getUnitVector();
//         rg_Point2D vecPrevInf( vecCurrToPrev.getY(), -vecCurrToPrev.getX());
// 
//         if ( generator[1].getRadius() > generator[0].getRadius() ) {
//             rg_Point2D focus = center[1];
// 
//             rg_Point2D vecFocusToTan     = centerOfTangentCircle - focus;
//             rg_Point2D vecFocusToPrevEnd = prevEndVtx->getLocation() - focus;
// 
//             rg_REAL angleFromInfToTan = angleFromVec1toVec2(vecPrevInf, vecFocusToTan);
//             rg_REAL angleFromInfToPrevEnd = angleFromVec1toVec2(vecPrevInf, vecFocusToPrevEnd);
// 
//             if ( rg_GE(angleFromInfToTan, angleFromInfToPrevEnd ) ){
//                 continue;
//             }
//         }
//         else {
//             rg_Point2D focus = center[0];
// 
//             rg_Point2D vecFocusToTan     = centerOfTangentCircle - focus;
//             rg_Point2D vecFocusToPrevEnd = prevEndVtx->getLocation() - focus;
// 
//             rg_REAL angleFromPrevEndToTan = angleFromVec1toVec2(vecFocusToPrevEnd, vecFocusToTan);
//             rg_REAL angleFromPrevEndToInf = angleFromVec1toVec2(vecFocusToPrevEnd, vecPrevInf);
// 
//             if ( rg_GE(angleFromPrevEndToTan, angleFromPrevEndToInf ) ) {
//                 continue;
//             }
//         }
// 
// 
//         //  curr vs next
//         rg_Point2D vecCurrToNext = (center[2] - center[1]).getUnitVector();
//         rg_Point2D vecNextInf( -vecCurrToNext.getY(), vecCurrToNext.getX());
// 
//         if ( generator[1].getRadius() > generator[2].getRadius() ) {
//             rg_Point2D focus = center[1];
// 
//             rg_Point2D vecFocusToTan     = centerOfTangentCircle - focus;
//             rg_Point2D vecFocusToNextEnd = nextEndVtx->getLocation() - focus;
// 
//             rg_REAL angleFromNextEndToTan = angleFromVec1toVec2(vecFocusToNextEnd, vecFocusToTan);
//             rg_REAL angleFromNextEndToInf = angleFromVec1toVec2(vecFocusToNextEnd, vecNextInf);
// 
//             if ( rg_GE(angleFromNextEndToTan, angleFromNextEndToInf ) ) {
//                 continue;
//             }
//         }
//         else {
//             rg_Point2D focus = center[2];
// 
//             rg_Point2D vecFocusToTan     = centerOfTangentCircle - focus;
//             rg_Point2D vecFocusToNextEnd = nextEndVtx->getLocation() - focus;
// 
//             rg_REAL angleFromNextInfToTan = angleFromVec1toVec2(vecNextInf, vecFocusToTan);
//             rg_REAL angleFromNextInfToNextEnd = angleFromVec1toVec2(vecNextInf, vecFocusToNextEnd);
// 
//             if ( rg_GE(angleFromNextInfToTan, angleFromNextInfToNextEnd ) ) {
//                 continue;
//             }
//         }
// 
//         ////////////////////////
// //         int LIMIT_RADIUS = 1e6;
// //         if( tangentCircle[i].getRadius() > LIMIT_RADIUS ) {
// //             break;
// //         }
//         ////////////////////////
// 
// 
//         isFlippingEdge = rg_TRUE;
//         circumcircle   = tangentCircle[i];
//         break;
//     }

    return isFlippingEdge;
}


double VoronoiDiagram2DC::angleFromVec1ToVec2( const rg_Point2D& vector1, const rg_Point2D& vector2 )
{
    rg_Point2D vec1 = vector1.getUnitVector();
    rg_Point2D vec2 = vector2.getUnitVector();
    rg_Point2D pt2pt( vec1 - vec2 );



    double length = pt2pt.magnitude();

    double cosine = (2. - length*length)/(2.);

    if( vec1 * vec2 > 0.0 )  //less than PI (cross product)
    {
        return acos( cosine );
    }
    else
    {
        return 2.*rg_PI - acos( cosine );
    }
}



void VoronoiDiagram2DC::removeRootNodeOfPriorityQueue( rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges, rg_dNode< pair<VEdge2D*, rg_Circle2D> >*& nodeWithSmallestCircle )
{


    double minRadius = DBL_MAX;
    rg_dNode< pair<VEdge2D*, rg_Circle2D> >* currNode = possiblyFlippingEdges.getFirstpNode();

    for( int i = 0 ; i < possiblyFlippingEdges.getSize() ; i++, currNode = currNode->getNext() ) {
        pair<VEdge2D*, rg_Circle2D>* currEdgeAndCircle = currNode->getpEntity();

        double radius = currEdgeAndCircle->second.getRadius();
        if( radius < minRadius ) {
            minRadius = radius;
            nodeWithSmallestCircle = currNode;
        }
    }

}



//void VoronoiDiagram2DC::flipThisEdgeNFindTwoIncidentTriplets( VFace2D* virtualFace, pair<VEdge2D*, rg_Circle2D>& edgeNCircle, VEdge2D*& incidentTriplet1, VEdge2D*& incidentTriplet2, rg_dList<VEdge2D*>& edgesToUpdateGeometry )
void VoronoiDiagram2DC::flipThisEdgeNFindTwoIncidentTriplets( VFace2D* virtualFace, pair<VEdge2D*, rg_Circle2D>& edgeNCircle, VEdge2D*& incidentTriplet1, VEdge2D*& incidentTriplet2 )
{
    VEdge2D*        flippingEdge = edgeNCircle.first;
    rg_Circle2D     circumcircle = edgeNCircle.second;

    flippingEdge->flip();

    VVertex2D*      startVtx = flippingEdge->getStartVertex();
    VVertex2D*      endVtx   = flippingEdge->getEndVertex();
    VFace2D*        matingFaceOfStartVtx = getMatingFace( startVtx, flippingEdge );

    if( matingFaceOfStartVtx == virtualFace ) {
        startVtx->setLocation( circumcircle.getCenterPt() );
        incidentTriplet1 = flippingEdge->getLeftHand();
        incidentTriplet2 = flippingEdge->getRightHand();
    }
    else {
        endVtx->setLocation( circumcircle.getCenterPt() );
        incidentTriplet1 = flippingEdge->getLeftLeg();
        incidentTriplet2 = flippingEdge->getRightLeg();
    }
}



void VoronoiDiagram2DC::reflectTheFlipOnTheIncidentTriplet( Priority_Q_VEdge& possiblyFlippingEdges, VEdge2D*& edgeOfIncidentTriplet )
{
    rg_Circle2D circumcircle;
    bool isFlippingEdge = isFlippingEdgeOnVirtualRegion( edgeOfIncidentTriplet, circumcircle );

    if( isFlippingEdge )
    {
        if( !edgeOfIncidentTriplet->isAlreadyCandidateForFlippingInPhantomRemoval() )
        {
            possiblyFlippingEdges.push( pair<VEdge2D*, rg_Circle2D>( edgeOfIncidentTriplet, circumcircle ) );
            edgeOfIncidentTriplet->setTrueCandidateForFlippingInPhantomRemoval();
        }
        else
        {
            possiblyFlippingEdges.changeCircumCircle( edgeOfIncidentTriplet, circumcircle );
        }
    }
    else
    {
        if( edgeOfIncidentTriplet->isAlreadyCandidateForFlippingInPhantomRemoval() )
        {
            possiblyFlippingEdges.kill( edgeOfIncidentTriplet );
            edgeOfIncidentTriplet->setFalseCandidateForFlippingInPhantomRemoval();
        }
    }
}



bool VoronoiDiagram2DC::isCorrectCircumcircle( rg_Circle2D& circumcircle, VFace2D* prevFace, VFace2D* currFace, VFace2D* nextFace )
{
    rg_Point2D vecCircumcircleToPrevFace = prevFace->getGenerator()->getDisk()->getCenterPt() - circumcircle.getCenterPt();
    rg_Point2D vecCircumcircleToCurrFace = currFace->getGenerator()->getDisk()->getCenterPt() - circumcircle.getCenterPt();
    rg_Point2D vecCircumcircleToNextFace = nextFace->getGenerator()->getDisk()->getCenterPt() - circumcircle.getCenterPt();
    

    if ( rg_POS(vecCircumcircleToNextFace*vecCircumcircleToCurrFace) && rg_POS(vecCircumcircleToCurrFace*vecCircumcircleToPrevFace) ) {
        return true;
    }
    else {
        return false;
    }
}

pair<VEdge2D*, rg_Circle2D> VoronoiDiagram2DC::dequeueEdgeWithMinTangentCircleFromPossiblyFlippingEdges( rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges )
{
    rg_dNode< pair<VEdge2D*, rg_Circle2D> >* nodeWithSmallestCircle = NULL;

    double minRadius = DBL_MAX;
    rg_dNode< pair<VEdge2D*, rg_Circle2D> >* currNode = possiblyFlippingEdges.getFirstpNode();

    for( int i = 0 ; i < possiblyFlippingEdges.getSize() ; i++, currNode = currNode->getNext() ) {
        pair<VEdge2D*, rg_Circle2D>* currEdgeAndCircle = currNode->getpEntity();

        double radius = currEdgeAndCircle->second.getRadius();
        if( radius < minRadius ) {
            minRadius = radius;
            nodeWithSmallestCircle = currNode;
        }
    }

    possiblyFlippingEdges.setHead( nodeWithSmallestCircle );

    return possiblyFlippingEdges.popFront();
}



int VoronoiDiagram2DC::getNumberOfVVertices()
{
    VFace2D* virtualFace = &(*m_VFaces.begin());
    list<VVertex2D*> infiniteVertices;
    virtualFace->getBoundaryVVertices( infiniteVertices );

    int     totalNumOfVertices      = m_VVertices.size();
    int     numOfInfiniteVertices   = infiniteVertices.size();
    int     numOfVVerticesWithoutInfiniteVertices = totalNumOfVertices - numOfInfiniteVertices;

    return numOfVVerticesWithoutInfiniteVertices;
}



int VoronoiDiagram2DC::getNumberOfVEdges()
{
    VFace2D* virtualFace = &(*m_VFaces.begin());
    list<VEdge2D*> VirtualEdges;
    virtualFace->getBoundaryVEdges( VirtualEdges );

    int     totalNumOfEdges         = m_VEdges.size();
    int     numOfVirtualEdges       = VirtualEdges.size();
    int     numOfVEdgesWithoutVirtualEdges = totalNumOfEdges - numOfVirtualEdges;

    return numOfVEdgesWithoutVirtualEdges;
}



int VoronoiDiagram2DC::getNumberOfVFaces()
{
    int totalNumOfFaces               = m_VFaces.size();
    int numOfVFacesWithoutVirtualFace = totalNumOfFaces - 1;

    return numOfVFacesWithoutVirtualFace;
}



void VoronoiDiagram2DC::resetAnomalyTestDoneToggle( list<VEdge2D*>& anomalyTestedVEdges )
{
    list<VEdge2D*>::iterator it_testedEdge = anomalyTestedVEdges.begin();
    for( ; it_testedEdge != anomalyTestedVEdges.end(); it_testedEdge++ ) {
        VEdge2D* currTestedEdge = *it_testedEdge;
        currTestedEdge->setFalseAnomalyTestDone();
    }
}




void VoronoiDiagram2DC::constructVoronoiDiagramWithoutPhantomDiskRemoval( list<rg_Circle2D>& circleset )
{
    createSortedGeneratorSet( circleset );

    if( IS_BUCKET_USED ) {
        m_bucket.setBucketCell( m_generators );
    }

    constructVoronoiDiagramWithoutPhantomDiskRemoval();
}




void VoronoiDiagram2DC::constructVoronoiDiagramWithoutPhantomDiskRemoval()
{
    insertPhantomGenerators();

    constructVoronoiDiagramForPhantomGenerators();

    for( list<Generator2D>::iterator it_generator = m_generators.begin() ; it_generator != m_generators.end() ; it_generator++ ) {
        Generator2D* currGenerator = &(*it_generator);
        updateVoronoiDiagramWithNewGenerator( currGenerator );

        if( IS_BUCKET_USED ) {
            m_bucket.addGenerator( currGenerator );
        }
    }

    removeAllExtraneousVVerticesAndVEdges();
}

void VoronoiDiagram2DC::printComputationTime( ofstream& fout )
{
    time2_2     = time2_2_1 + time2_2_2 + time2_2_3 + time2_2_4;
    time2       = time2_1 + time2_2 + time2_3 + time2_4 + time2_5 + time2_6; 

    fout << time1 << "\t" << time2 << "\t" << time2_1 << "\t" << time2_2 << "\t" << time2_2_1 << "\t" << time2_2_2 << "\t" << time2_2_3 << "\t" << time2_2_4 << "\t" << time2_3 << "\t" << time2_4 << "\t" << time2_5 << "\t" << time2_6 << "\t" << time3 << "\t" << time4 << "\t" << phantom1 << "\t" << phantom2 << "\t" << phantom3 << "\t" << phantom4 << "\t" << phantom5 << "\t" << phantom6 << endl;
}

void VoronoiDiagram2DC::makeWEDS( const char* fname, list<rg_Circle2D>& circleList )
{
    string      fileNameWithPath = (string)fname;
    ifstream    fin;
    fin.open(fileNameWithPath);

    char        buffer[100];
    const char* seps = " \t\n[]";
    char*       token = NULL;

    fin.getline(buffer, 100);
    token = strtok(buffer, seps);
    
    int numOfFaces = atoi(token);
    token = strtok(NULL, seps);

    int numOfVertices = atoi(token);
    fin.getline(buffer, 100);
    fin.getline(buffer, 100);

    

    // make face and generator
    // connect face and generator
    // make face list to array
    VFace2D** faceArray = new VFace2D*[numOfFaces];

    m_VFaces.push_back(VFace2D(0));

    for( int i = 1 ; i < numOfFaces ; i++ ) {
        m_VFaces.push_back(VFace2D(i));
        VFace2D* currFace = &(*m_VFaces.rbegin());

        double x    = atof(strtok(buffer, seps));
        double y    = atof(strtok(NULL, seps));
        double r    = atof(strtok(NULL, seps));

        circleList.push_back( rg_Circle2D(x, y, r) );
        rg_Circle2D* currCircle = &(*circleList.rbegin());

        m_generators.push_back( Generator2D(i, currCircle, currFace) );
        Generator2D* currGene = &(*m_generators.rbegin());

        currFace->setGenerator(currGene);
    }
    



}

void VoronoiDiagram2DC::getFaces( list<VFace2D*>& faces )
{
    list<VFace2D>::iterator i_face = m_VFaces.begin();
    for( ; i_face != m_VFaces.end() ; i_face++ ) {
        VFace2D* currFace = &(*i_face);
        faces.push_back( currFace );
    }
}

void VoronoiDiagram2DC::getVertices( list<VVertex2D*>& vertices )
{
    list<VVertex2D>::iterator i_vertex = m_VVertices.begin();
    for( ; i_vertex != m_VVertices.end() ; i_vertex++ ) {
        VVertex2D* currVertex = &(*i_vertex);
        vertices.push_back( currVertex );
    }
}






