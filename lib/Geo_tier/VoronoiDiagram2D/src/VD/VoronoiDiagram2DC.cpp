#include "VoronoiDiagram2DC.h"
#include "rg_Point3D.h"
#include "rg_Point2D.h"
#include "rg_GeoFunc.h"
#include "Priority_Q_VEdge.h"
#include "VertexGenerator2D.h"
using namespace V::GeometryTier;

#include <algorithm>
#include <vector>
#include <fstream>
#include <time.h>
#include <map>
using namespace std;


VoronoiDiagram2DC::VoronoiDiagram2DC(void)
{
	m_directAddressTableFromDiskIDToDiskArrangement = rg_NULL; // disk arragement in double precision by Joonghyun March, 17
}



VoronoiDiagram2DC::VoronoiDiagram2DC( list<rg_Circle2D>& circleset )
{
    int currID = 1;
    for( list<rg_Circle2D>::iterator i_disk = circleset.begin() ; i_disk != circleset.end() ; ++i_disk, ++currID ) {
        //m_generators.push_back( Generator2D( currID, *it_point ) );
        Generator2D* newGen = createGenerator();
        newGen->setID( currID );
        newGen->setDisk( *i_disk );
    }

	m_directAddressTableFromDiskIDToDiskArrangement = rg_NULL; // disk arragement in double precision by Joonghyun March, 17
}



VoronoiDiagram2DC::VoronoiDiagram2DC( const VoronoiDiagram2DC& VD2DC )
{
	copyVoronoiDiagramFrom(VD2DC);
}



VoronoiDiagram2DC::~VoronoiDiagram2DC(void)
{
    if (!m_generators.empty())
    {
        clear();
    }
    
    if (m_directAddressTableFromDiskIDToDiskArrangement != rg_NULL)
    {
        delete[] m_directAddressTableFromDiskIDToDiskArrangement;
        m_directAddressTableFromDiskIDToDiskArrangement = rg_NULL;
    }
}



void VoronoiDiagram2DC::getVoronoiEdges( list<const VEdge2D*>& VEdgesList ) const
{
    for( list<VEdge2D*>::const_iterator it_edge = m_VEdges.begin() ; it_edge != m_VEdges.end() ; ++it_edge ) {
        const VEdge2D* currEdge = *it_edge;
        VEdgesList.push_back( currEdge );
    }
}



void VoronoiDiagram2DC::getVoronoiFaces( list<const VFace2D*>& VFacesList ) const
{
    for( list<VFace2D*>::const_iterator it_face = m_VFaces.begin() ; it_face != m_VFaces.end() ; ++it_face ) {
        const VFace2D* currFace = *it_face;
        VFacesList.push_back( currFace );
    }
}



void VoronoiDiagram2DC::getVoronoiVertices( list<const VVertex2D*>& VVerticesList ) const
{
    for( list<VVertex2D*>::const_iterator it_vertex = m_VVertices.begin() ; it_vertex != m_VVertices.end() ; ++it_vertex ) {
        const VVertex2D* currVertex = *it_vertex;
        VVerticesList.push_back( currVertex );
    }
}



void VoronoiDiagram2DC::getGenerators( list<const Generator2D*>& generatorsList ) const
{
    for( list<Generator2D*>::const_iterator it_generator = m_generators.begin() ; it_generator != m_generators.end() ; ++it_generator ) {
        const Generator2D* currGenerator = *it_generator;
        generatorsList.push_back( currGenerator );
    }
}




void VoronoiDiagram2DC::getVoronoiEdges(list<VEdge2D*>& VEdgesList) const
{
    for (list<VEdge2D*>::const_iterator it_edge = m_VEdges.begin(); it_edge != m_VEdges.end(); ++it_edge) {
        const VEdge2D* currEdge = *it_edge;
        VEdgesList.push_back(const_cast<VEdge2D*>(currEdge));
    }
}



void VoronoiDiagram2DC::getVoronoiFaces(list<VFace2D*>& VFacesList) const
{
    for (list<VFace2D*>::const_iterator it_face = m_VFaces.begin(); it_face != m_VFaces.end(); ++it_face) {
        const VFace2D* currFace = *it_face;
        VFacesList.push_back(const_cast<VFace2D*>(currFace));
    }
}



void VoronoiDiagram2DC::getVoronoiVertices(list<VVertex2D*>& VVerticesList) const
{
    for (list<VVertex2D*>::const_iterator it_vertex = m_VVertices.begin(); it_vertex != m_VVertices.end(); ++it_vertex) {
        const VVertex2D* currVertex = *it_vertex;
        VVerticesList.push_back(const_cast<VVertex2D*>(currVertex));
    }
}


void VoronoiDiagram2DC::getGenerators(list<Generator2D*>& generatorsList) const
{
    for (list<Generator2D*>::const_iterator it_generator = m_generators.begin(); it_generator != m_generators.end(); ++it_generator) {
        const Generator2D* currGenerator = *it_generator;
        generatorsList.push_back(const_cast<Generator2D*>(currGenerator));
    }
}




void VoronoiDiagram2DC::setGenerators( list<rg_Circle2D>& circleset ) //sorting ?ÔøΩÍ∏∞Ôø? Îπ†Ïßê
    // createSortedGeneratorSet
{
    circleset.sort(compareCircleDescendingorder);
    m_generators.clear();

    int currID = 1;
    for( list<rg_Circle2D>::iterator i_disk = circleset.begin() ; i_disk != circleset.end() ; ++i_disk, ++currID ) {
        rg_Circle2D currDisk = *i_disk;
        float xCoord = currDisk.getCenterPt().getX();
        float yCoord = currDisk.getCenterPt().getY();
        float radius = currDisk.getRadius();

        rg_Circle2D currDisk_single_precision( xCoord, yCoord, radius );
        Generator2D* newGenerator = createGenerator();
        newGenerator->setID(currID);
        newGenerator->setDisk(currDisk_single_precision);
    }
}



void VoronoiDiagram2DC::setGenerators( list< pair< rg_Circle2D, void* > >& circleNUserDataPairSet )
{
    circleNUserDataPairSet.sort(compareCirclePairDescendingorder);
    m_generators.clear();

    int currID = 1;
    for( list< pair< rg_Circle2D, void* > >::iterator i_disk = circleNUserDataPairSet.begin() ; i_disk != circleNUserDataPairSet.end() ; ++i_disk, ++currID ) {
        pair< rg_Circle2D, void* > currDiskPair = *i_disk;
        float xCoord = currDiskPair.first.getCenterPt().getX();
        float yCoord = currDiskPair.first.getCenterPt().getY();
        float radius = currDiskPair.first.getRadius();

        rg_Circle2D currDisk_single_precision( xCoord, yCoord, radius );
        Generator2D* newGenerator = createGenerator();
        newGenerator->setID(currID);
        newGenerator->setDisk(currDisk_single_precision);
        newGenerator->setUserData(currDiskPair.second);
    }
}



void VoronoiDiagram2DC::setGenerators( list< rg_Triplet< rg_Circle2D, int, void* > >& circleAndIDAndUserDataSet )
{
    circleAndIDAndUserDataSet.sort(compareCircleTripletDescendingorder);
    m_generators.clear();

    for( list< rg_Triplet< rg_Circle2D, int, void* > >::iterator i_disk = circleAndIDAndUserDataSet.begin() ; i_disk != circleAndIDAndUserDataSet.end() ; ++i_disk ) {
        rg_Triplet< rg_Circle2D, int, void* > currDiskTriplet = *i_disk;
        float xCoord = currDiskTriplet.m_first.getCenterPt().getX();
        float yCoord = currDiskTriplet.m_first.getCenterPt().getY();
        float radius = currDiskTriplet.m_first.getRadius();

        rg_Circle2D currDisk_single_precision( xCoord, yCoord, radius );
        Generator2D* newGenerator = createGenerator();
        newGenerator->setID(currDiskTriplet.m_second);
        newGenerator->setDisk(currDisk_single_precision);
        newGenerator->setUserData(currDiskTriplet.m_third);
    }
}



void VoronoiDiagram2DC::setGenerators( list<Generator2D*>& generators )
{
    m_generators.clear();
    m_generators = generators;
}



void VoronoiDiagram2DC::constructVoronoiDiagram( list<rg_Circle2D>& circleset )
{
    setGenerators( circleset );
    
    if( IS_BUCKET_USED ) {
        m_bucket.setBucketCell( m_generators );
    }

    constructVoronoiDiagram();
}



void VoronoiDiagram2DC::constructVoronoiDiagram( list< pair< rg_Circle2D, void* > >& circleNUserDataPairSet )
{
    setGenerators( circleNUserDataPairSet );

    if( IS_BUCKET_USED ) {
        m_bucket.setBucketCell( m_generators );
    }

    constructVoronoiDiagram();
}



void VoronoiDiagram2DC::constructVoronoiDiagram( list< rg_Triplet< rg_Circle2D, int, void* > >& circleAndIDAndUserDataSet )
{
    setGenerators( circleAndIDAndUserDataSet );

    if( IS_BUCKET_USED ) {
        m_bucket.setBucketCell( m_generators );
    }

    constructVoronoiDiagram();
}



void VoronoiDiagram2DC::constructVoronoiDiagramWithoutSorting( list<rg_Circle2D>& circleset )
{
    setGeneratorsWithoutSorting( circleset );

    if( IS_BUCKET_USED ) {
        m_bucket.setBucketCell( m_generators );
    }

    constructVoronoiDiagram();
}



void VoronoiDiagram2DC::constructVoronoiDiagram()
{
    insertPhantomGenerators();

    constructVoronoiDiagramForPhantomGenerators();

    for( list<Generator2D*>::iterator i_gen = m_generators.begin() ; i_gen != m_generators.end() ; ++i_gen ) {
        Generator2D* currGenerator = *i_gen;
        updateVoronoiDiagramWithNewGenerator( currGenerator );

        if( IS_BUCKET_USED ) {
            m_bucket.addGenerator( currGenerator );
        }
    }

    removePhantomGenerators();

    removeAllExtraneousVVerticesAndVEdges();
}



void VoronoiDiagram2DC::constructVoronoiDiagramWithoutPhantomDiskRemoval( list<rg_Circle2D>& circleset )
{
    setGenerators(circleset);

    if( IS_BUCKET_USED ) {
        m_bucket.setBucketCell( m_generators );
    }

    constructVoronoiDiagramWithoutPhantomDiskRemoval();
}



void VoronoiDiagram2DC::constructVoronoiDiagramWithoutPhantomDiskRemoval( list< pair< rg_Circle2D, void* > >& circleNUserDataPairSet )
{
    setGenerators(circleNUserDataPairSet);

    if( IS_BUCKET_USED ) {
        m_bucket.setBucketCell( m_generators );
    }

    constructVoronoiDiagramWithoutPhantomDiskRemoval();
}



void VoronoiDiagram2DC::constructVoronoiDiagramWithoutPhantomDiskRemoval( list< rg_Triplet< rg_Circle2D, int, void* > >& circleAndIDAndUserDataSet )
{
    setGenerators(circleAndIDAndUserDataSet);

    if( IS_BUCKET_USED ) {
        m_bucket.setBucketCell( m_generators );
    }

    constructVoronoiDiagramWithoutPhantomDiskRemoval();
}



void VoronoiDiagram2DC::constructVoronoiDiagramWithoutPhantomDiskRemoval()
{
    insertPhantomGenerators();

    constructVoronoiDiagramForPhantomGenerators();

    int numInsertedDisks = 0;
    for( list<Generator2D*>::iterator i_gen = m_generators.begin() ; i_gen != m_generators.end() ; ++i_gen ) {
        Generator2D* currGenerator = *i_gen;
        updateVoronoiDiagramWithNewGenerator( currGenerator );

        if( IS_BUCKET_USED ) {
            m_bucket.addGenerator( currGenerator );
        }

        ++numInsertedDisks;
        if ( numInsertedDisks % 5000 == 0 ) {
            //AfxMessageBox(_T("10000 complete"));
            removeAllExtraneousVVerticesAndVEdges();
        }
    }

    removeAllExtraneousVVerticesAndVEdges();
}



VFace2D* VoronoiDiagram2DC::findVFaceContainingQueryPoint( const rg_Point2D& pt ) const
{
    VFace2D* VFaceContainingPoint = NULL;
    if (IS_BUCKET_USED)
    {
        Generator2D* closeGenerator = m_bucket.findGoodCloseGenerator(pt);
        if (closeGenerator != NULL)
        {
            VFaceContainingPoint = closeGenerator->getOuterFace();
        }
        else
        {
            VFaceContainingPoint = m_VFaces.back();
        }
    }
    else
    {
        VFaceContainingPoint = m_VFaces.back();
    }

    rg_Circle2D     disk = VFaceContainingPoint->getGenerator()->getDisk();
	rg_Circle2D*    currDisk = &disk;
		
    double          minDistancePt2GeneratorOnVFace;
    if ( (VoronoiDiagram2DC*)VFaceContainingPoint->getGenerator()->getInnerVD() == this ) {
        minDistancePt2GeneratorOnVFace = currDisk->getRadius() - pt.distance( currDisk->getCenterPt() );
    }
    else {
        minDistancePt2GeneratorOnVFace = pt.distance( currDisk->getCenterPt() ) - currDisk->getRadius();
    }

    VEdge2D*        beginningEdgeOfFace = VFaceContainingPoint->getFirstVEdge();
    VEdge2D*        currBoundaryEdge = beginningEdgeOfFace;
    VFace2D*        currOppositeFace = getFaceOfOppositeSide( VFaceContainingPoint, currBoundaryEdge );

    Generator2D*    neighborGenerator = NULL;
    double          distanceFromNeighbor2InputPt = DBL_MAX;


    do {
        currOppositeFace    = getFaceOfOppositeSide( VFaceContainingPoint, currBoundaryEdge );
        neighborGenerator   = currOppositeFace->getGenerator();

        if( neighborGenerator != NULL ) {
            //currDisk = &neighborGenerator->getDisk();
            rg_Circle2D     disk = neighborGenerator->getDisk();
            currDisk = &disk;

            if ( (VoronoiDiagram2DC*)neighborGenerator->getInnerVD() == this ) {
                distanceFromNeighbor2InputPt = currDisk->getRadius() - pt.distance( currDisk->getCenterPt() );
            }
            else {
                distanceFromNeighbor2InputPt = pt.distance( currDisk->getCenterPt() ) - currDisk->getRadius();
            }            
        }

        if( distanceFromNeighbor2InputPt < minDistancePt2GeneratorOnVFace ) {
            VFaceContainingPoint        = currOppositeFace;
            minDistancePt2GeneratorOnVFace = distanceFromNeighbor2InputPt;
            beginningEdgeOfFace = currBoundaryEdge;
        }

        currBoundaryEdge = getCCWNextEdgeOnFace( currBoundaryEdge, VFaceContainingPoint );

    } while ( currBoundaryEdge != beginningEdgeOfFace );

    return VFaceContainingPoint;
}


VVertex2D* VoronoiDiagram2DC::findClosestVVertexToPoint(const rg_Point2D& pt) const
{
    VFace2D*  faceContainingPoint = findVFaceContainingQueryPoint(pt);

    list<VVertex2D*> boundaryVertices;
    faceContainingPoint->getBoundaryVVertices(boundaryVertices);

    VVertex2D* closestVertex = rg_NULL;
    double     minDistance = DBL_MAX;
    for (list<VVertex2D*>::iterator i_vtx = boundaryVertices.begin(); i_vtx != boundaryVertices.end(); ++i_vtx) {
        VVertex2D* currVtx = (*i_vtx);
        if (currVtx->isInfinite()) continue;

        double dist = currVtx->getCircumcircle().getCenterPt().distance(pt);
        if (dist < minDistance) {
            closestVertex = currVtx;
            minDistance = dist;
        }
    }

    return closestVertex;
}


VVertex2D* VoronoiDiagram2DC::findVisibleClosestVVertexToPoint(const rg_Point2D& pt) const
{
    //1. find intersecting disks
    vector<Generator2D*> intersectingNeighbors;

    VFace2D*   faceContainingPoint = findVFaceContainingQueryPoint(pt);
 
    list<VEdge2D*> boundaryEdges;
    faceContainingPoint->getBoundaryVEdges(boundaryEdges);

    for (list<VEdge2D*>::const_iterator it_Edge = boundaryEdges.begin(); it_Edge != boundaryEdges.end(); ++it_Edge)
    {
        VEdge2D* currVEdge = *it_Edge;

        if (!currVEdge->is_passable(0.0))
        {
            VFace2D* firstNeighborFace = NULL;

            if (currVEdge->getLeftFace() == faceContainingPoint)
            {
                firstNeighborFace = currVEdge->getRightFace();
            }
            else //(currVEdge->getRightFace() == faceContainingPoint)
            {
                firstNeighborFace = currVEdge->getLeftFace();
            }

            intersectingNeighbors.push_back(firstNeighborFace->getGenerator());
        }
    }

    VVertex2D* closestVertex = rg_NULL;

    if (intersectingNeighbors.size() < 2)
    {
        closestVertex = findClosestVVertexToPoint(pt);
    }
    else
    {
        //2. find interval including query point
        rg_Point2D centerPtOfFace = faceContainingPoint->getGenerator()->getDisk().getCenterPt();

        vector<double> angleOfIntersectingNeighbors;

        const int numOfIntersectingNeighbors = intersectingNeighbors.size();

        for (int i = 0; i < numOfIntersectingNeighbors; ++i)
        {
            rg_Point2D centerPtOfCurrFirstNeighbor = intersectingNeighbors[i]->getDisk().getCenterPt();
            rg_Point2D vecToFirstNeigbor           = centerPtOfCurrFirstNeighbor - centerPtOfFace;
            double     angle                       = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vecToFirstNeigbor);

            angleOfIntersectingNeighbors.push_back(angle);
        }

        sort(angleOfIntersectingNeighbors.begin(), angleOfIntersectingNeighbors.end());
        angleOfIntersectingNeighbors.push_back(angleOfIntersectingNeighbors[0] + 2 * rg_PI);

        rg_Point2D vecToQuenryPt  = pt - centerPtOfFace;
        double     angleOfQueryPt = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vecToQuenryPt);
        if (angleOfQueryPt < angleOfIntersectingNeighbors[0]) {
            angleOfQueryPt += 2 * rg_PI;
        }

        pair<double, double> intervalOfAngleIncludingQueryPt(DBL_MAX, DBL_MAX);

        for (int i = 0; i < numOfIntersectingNeighbors; ++i)
        {
            if ((angleOfQueryPt >= angleOfIntersectingNeighbors[i]) && (angleOfQueryPt < angleOfIntersectingNeighbors[i + 1]))
            {
                intervalOfAngleIncludingQueryPt.first  = angleOfIntersectingNeighbors[i];
                intervalOfAngleIncludingQueryPt.second = angleOfIntersectingNeighbors[i + 1];
                break;
            }
        }

        //3. find vertices in the angle interval
        list<VVertex2D*> verticesInAngleInterval;

        list<VVertex2D*> boundaryVertices;
        faceContainingPoint->getBoundaryVVertices(boundaryVertices);

        for (list<VVertex2D*>::const_iterator it_Vertex = boundaryVertices.begin(); it_Vertex != boundaryVertices.end(); ++it_Vertex)
        {
            VVertex2D* currBoundaryVertex = *it_Vertex;
            if (currBoundaryVertex->isInfinite())
            {
                continue;
            }
            
            rg_Point2D vecToCurrVertexPt   = currBoundaryVertex->getLocation() - centerPtOfFace;
            double     angleOfCurrVertexPt = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vecToCurrVertexPt);

            if (angleOfCurrVertexPt < intervalOfAngleIncludingQueryPt.first) {
                angleOfCurrVertexPt += 2 * rg_PI;
            }

            if ((angleOfCurrVertexPt >= intervalOfAngleIncludingQueryPt.first) && (angleOfCurrVertexPt < intervalOfAngleIncludingQueryPt.second))
            {
                verticesInAngleInterval.push_back(currBoundaryVertex);
            }
        }

        //4. find vertices containing the input pt
        list<VVertex2D*> verticesOfCircucirclesContainingPt;
        for (list<VVertex2D*>::iterator i_vtx = verticesInAngleInterval.begin(); i_vtx != verticesInAngleInterval.end(); ++i_vtx)
        {
            VVertex2D* currVtx = (*i_vtx);
            if (currVtx->isInfinite()) continue;

            if (currVtx->getCircumcircle().distance(pt) < 0) {
                verticesOfCircucirclesContainingPt.push_back(currVtx);
            }
        }


        if (verticesOfCircucirclesContainingPt.size() == 0)
        {
            //5. find min distance vertices in the angle interval
            double     minDistance = DBL_MAX;
            for (list<VVertex2D*>::iterator i_vtx = verticesInAngleInterval.begin(); i_vtx != verticesInAngleInterval.end(); ++i_vtx)
            {
                VVertex2D* currVtx = (*i_vtx);
                if (currVtx->isInfinite()) continue;

                double dist = currVtx->getCircumcircle().distance(pt);
                if (dist < minDistance) {
                    closestVertex = currVtx;
                    minDistance = dist;
                }
            }
        }
        else
        {
            //5. find min distance vertices in the angle interval
            double     minDistance = DBL_MAX;
            for (list<VVertex2D*>::iterator i_vtx = verticesOfCircucirclesContainingPt.begin(); i_vtx != verticesOfCircucirclesContainingPt.end(); ++i_vtx)
            {
                VVertex2D* currVtx = (*i_vtx);
                if (currVtx->isInfinite()) continue;

                double dist = currVtx->getCircumcircle().getCenterPt().distance(pt);
                if (dist < minDistance) {
                    closestVertex = currVtx;
                    minDistance = dist;
                }
            }
        }
    }

    return closestVertex;
}


VVertex2D* VoronoiDiagram2DC::findVisibleNearestVVertexToPoint(const rg_Point2D& queryPoint, const double& probeRadius /*= 0.0*/) const
{
    //1. find intersecting disks
    list<double> anglesOfIntersectingNeighbors;

    VFace2D* faceContainingPoint = findVFaceContainingQueryPoint(queryPoint);
    rg_Point2D centerOfGeneratorWhoseVFaceContainingQueryPoint = faceContainingPoint->getGenerator()->getDisk().getCenterPt();
 
    list<VEdge2D*> boundaryEdges;
    faceContainingPoint->getBoundaryVEdges(boundaryEdges);

    for (list<VEdge2D*>::const_iterator it_Edge = boundaryEdges.begin(); 
         it_Edge != boundaryEdges.end(); 
         ++it_Edge)
    {
        VEdge2D* currVEdge = *it_Edge;

        if (!currVEdge->radius_interval_is_calculated())
        {
            currVEdge->compute_radius_interval_of_tangent_circles();
        }

        if (!currVEdge->is_passable(probeRadius))
        {
            Generator2D* firstNeighborGenerator = NULL;

            if (currVEdge->getLeftFace() == faceContainingPoint)
            {
                firstNeighborGenerator = currVEdge->getRightFace()->getGenerator();
            }
            else //(currVEdge->getRightFace() == faceContainingPoint)
            {
                firstNeighborGenerator = currVEdge->getLeftFace()->getGenerator();
            }

            rg_Point2D vectorFromCenterToNeighborCenter = firstNeighborGenerator->getDisk().getCenterPt() 
                                                        - centerOfGeneratorWhoseVFaceContainingQueryPoint;

            double angleOfFirstNeighborGenerator = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vectorFromCenterToNeighborCenter);
            anglesOfIntersectingNeighbors.push_back(angleOfFirstNeighborGenerator);
        }
    }
    

    if (anglesOfIntersectingNeighbors.empty())
    {
        list<VVertex2D*> boundingVVerticesOfFaceIncludingQueryPt;
        faceContainingPoint->getBoundaryVVertices(boundingVVerticesOfFaceIncludingQueryPt);

        VVertex2D* visibleNearestVVertex = boundingVVerticesOfFaceIncludingQueryPt.front();
        double     minDistance           = visibleNearestVVertex->getCircumcircle().distance(queryPoint);

        for (list<VVertex2D*>::const_iterator it_VVertex = next(boundingVVerticesOfFaceIncludingQueryPt.begin());
            it_VVertex != boundingVVerticesOfFaceIncludingQueryPt.end();
            ++it_VVertex)
        {
            VVertex2D* currVertex   = *it_VVertex;
            double     currDistance = currVertex->getCircumcircle().distance(queryPoint);

            if (currDistance < minDistance)
            {
                visibleNearestVVertex = currVertex;
                minDistance           = currDistance;
            }
        }

        return visibleNearestVVertex;
    }
    else //i.e. intersection with neighbor disks exists.
    {    
        //2. adjust intersection angle interval and angle of query point 

        if (anglesOfIntersectingNeighbors.size() == 1)
        {
            double onlyOneIntersectingAngle = anglesOfIntersectingNeighbors.front();
            if (onlyOneIntersectingAngle >= rg_PI)
            {
                anglesOfIntersectingNeighbors.push_front(onlyOneIntersectingAngle - rg_PI);
                anglesOfIntersectingNeighbors.push_back(onlyOneIntersectingAngle + rg_PI);
            }
            else
            {
                anglesOfIntersectingNeighbors.push_back(onlyOneIntersectingAngle + rg_PI);
                anglesOfIntersectingNeighbors.push_back(onlyOneIntersectingAngle + 2.0 * rg_PI);
            }
        }
        else if (anglesOfIntersectingNeighbors.size() > 1)
        {
            anglesOfIntersectingNeighbors.sort();
            anglesOfIntersectingNeighbors.push_back(anglesOfIntersectingNeighbors.front() + 2.0 * rg_PI);
        }

        rg_Point2D vectorFromCenterToQueryPoint = queryPoint - centerOfGeneratorWhoseVFaceContainingQueryPoint;
        double angleOfQueryPoint = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vectorFromCenterToQueryPoint);

        if (angleOfQueryPoint < anglesOfIntersectingNeighbors.front())
        {
            angleOfQueryPoint += 2.0 * rg_PI;
        }


        //3. find where the query point is included 
        double angleIntervalIncludingQueryPoint[2] = { DBL_MAX, DBL_MAX };
        for (list<double>::const_iterator it_Angle = anglesOfIntersectingNeighbors.begin();
            it_Angle != prev(anglesOfIntersectingNeighbors.end());
            ++it_Angle)
        {
            double currAngle = *it_Angle;
            double nextAngle = *next(it_Angle);

            if (angleOfQueryPoint >= currAngle && angleOfQueryPoint < nextAngle)
            {
                angleIntervalIncludingQueryPoint[0] = currAngle;
                angleIntervalIncludingQueryPoint[1] = nextAngle;
            }
        }

        //4. find all VVertices in the interval which contains the query point
        list<VVertex2D*> VVerticesInSameAngleInterval;

        list<VVertex2D*> boundaryVVertices;
        faceContainingPoint->getBoundaryVVertices(boundaryVVertices);

        for (list<VVertex2D*>::const_iterator it_VVertex = boundaryVVertices.begin();
            it_VVertex != boundaryVVertices.end();
            ++it_VVertex)
        {
            VVertex2D* currVertex = *it_VVertex;

            rg_Point2D vectorFromCenterToCurrVtx = currVertex->getLocation() - centerOfGeneratorWhoseVFaceContainingQueryPoint;
            double angleOfCurrVtx = angleFromVec1toVec2(rg_Point2D(1.0, 0.0), vectorFromCenterToCurrVtx);

            if (angleOfCurrVtx >= angleIntervalIncludingQueryPoint[0]
                && angleOfCurrVtx < angleIntervalIncludingQueryPoint[1])
            {
                VVerticesInSameAngleInterval.push_back(currVertex);
            }
            else if (angleOfCurrVtx + 2.0*rg_PI >= angleIntervalIncludingQueryPoint[0]
                && angleOfCurrVtx + 2.0*rg_PI < angleIntervalIncludingQueryPoint[1])
            {
                VVerticesInSameAngleInterval.push_back(currVertex);
            }
        }


        //5. find nearest VVertex
        VVertex2D* visibleNearestVVertex = VVerticesInSameAngleInterval.front();
        double     minDistance           = visibleNearestVVertex->getCircumcircle().distance(queryPoint);

        for (list<VVertex2D*>::const_iterator it_VVertex = next(VVerticesInSameAngleInterval.begin());
            it_VVertex != VVerticesInSameAngleInterval.end();
            ++it_VVertex)
        {
            VVertex2D* currVertex = *it_VVertex;
            double     currDistance = currVertex->getCircumcircle().distance(queryPoint);

            if (currDistance < minDistance)
            {
                visibleNearestVVertex = currVertex;
                minDistance = currDistance;
            }
        }

        return visibleNearestVVertex;
    }
}



void VoronoiDiagram2DC::clear()
{

    /*
	for (list<VVertex2D*>::iterator i_vtx = m_VVertices.begin(); i_vtx != m_VVertices.end(); ++i_vtx) {
		VVertex2D* vertex = *i_vtx;
		delete vertex;
	}

	for (list<VEdge2D*>::iterator i_edge = m_VEdges.begin(); i_edge != m_VEdges.end(); ++i_edge) {
		VEdge2D* edge = *i_edge;
		delete edge;
	}

	for (list<VFace2D*>::iterator i_face = m_VFaces.begin(); i_face != m_VFaces.end(); ++i_face) {
		VFace2D* face = *i_face;
		delete face;
	}

	for (list<Generator2D*>::iterator i_gen = m_generators.begin(); i_gen != m_generators.end(); ++i_gen) {
		Generator2D* generator = *i_gen;

        //«ˆ¿Á VD¿« container∞° æ∆¥— ∞ÊøÏ∏∏ delete∏¶ Ω√≈¥
        //«ˆ¿Á VD¿« container (NULL¿Œ ªÛ≈¬)¥¬ ¿ÃπÃ ªÛ¿ßø°º≠ delete∞° µ 
		
		VoronoiDiagram2DC* currentInnerVD = (VoronoiDiagram2DC*)generator->getInnerVD();
		if (currentInnerVD != this)
		{
			//if(currentInnerVD != NULL) // Added by Joonghyun on August 14, 2020
			delete generator;
		}
	}

	for (list<Generator2D*>::iterator i_phantom = m_phantomGenerators.begin(); i_phantom != m_phantomGenerators.end(); ++i_phantom) {
		Generator2D* phantomGen = *i_phantom;
		delete phantomGen;
	}


	m_VVertices.clear();
	m_VEdges.clear();
	m_VFaces.clear();
	m_generators.clear();
	m_phantomCircles.clear();
	m_phantomGenerators.clear();

	m_bucket.destroyBucket();
    */
    
	clearTopologyUnrelatedObjects();
	clearVEntities();
	clearGenerators();
    

}



void VoronoiDiagram2DC::clearVEntities()
{
	for (list<VVertex2D*>::iterator i_vtx = m_VVertices.begin(); i_vtx != m_VVertices.end(); ++i_vtx) {
		VVertex2D* vertex = *i_vtx;
		delete vertex;
	}

	for (list<VEdge2D*>::iterator i_edge = m_VEdges.begin(); i_edge != m_VEdges.end(); ++i_edge) {
		VEdge2D* edge = *i_edge;
		delete edge;
	}

	for (list<VFace2D*>::iterator i_face = m_VFaces.begin(); i_face != m_VFaces.end(); ++i_face) {
		VFace2D* face = *i_face;
		delete face;
	}

	m_VVertices.clear();
	m_VEdges.clear();
	m_VFaces.clear();
}


void VoronoiDiagram2DC::clearGenerators()
{
	for (list<Generator2D*>::iterator i_gen = m_generators.begin(); i_gen != m_generators.end(); ++i_gen) {
		Generator2D* generator = *i_gen;

		//ÔøΩÔøΩÔøΩÔøΩ VDÔøΩÔøΩ containerÔøΩÔøΩ ÔøΩ∆¥ÔøΩ ÔøΩÔøΩÔø??deleteÔøΩÔøΩ ÔøΩÔøΩ≈¥
		//ÔøΩÔøΩÔøΩÔøΩ VDÔøΩÔøΩ container (NULLÔøΩÔøΩ ÔøΩÔøΩÔøΩÔøΩ)ÔøΩÔøΩ ÔøΩÃπÔøΩ ÔøΩÔøΩÔøΩÔøΩÔøΩÔøΩÔøΩÔøΩ deleteÔøΩÔøΩ ÔøΩÔøΩ
		VoronoiDiagram2DC* currentInnerVD = (VoronoiDiagram2DC*)generator->getInnerVD();
		if (currentInnerVD != this)
		{
			delete generator;
		}
	}

	for (list<Generator2D*>::iterator i_phantom = m_phantomGenerators.begin(); i_phantom != m_phantomGenerators.end(); ++i_phantom) {
		Generator2D* phantomGen = *i_phantom;
		delete phantomGen;
	}

	m_generators.clear();
	m_phantomGenerators.clear();
}


void VoronoiDiagram2DC::clearGeneratorsWithoutDeletingDynamicMemoryOfTheGenerators()
{
    m_generators.clear();
}

void VoronoiDiagram2DC::clearTopologyUnrelatedObjects()
{
	m_phantomCircles.clear();
	m_bucket.destroyBucket();
}


VVertex2D* VoronoiDiagram2DC::createVVertex( VVertex2D* const vertex )
{
    m_VVertices.push_back( vertex );

    return m_VVertices.back();
}



VEdge2D* VoronoiDiagram2DC::createVEdge( VEdge2D* const edge )
{
    m_VEdges.push_back( edge );

    return m_VEdges.back();
}



VFace2D* VoronoiDiagram2DC::createVFace( VFace2D* const face )
{
    m_VFaces.push_back( face );

    return m_VFaces.back();
}



VVertex2D* VoronoiDiagram2DC::createVertex( const int& id )
{
    VVertex2D* newVertex = new VVertex2D(id);
    m_VVertices.push_back(newVertex);
    return newVertex;
}



VEdge2D* VoronoiDiagram2DC::createEdge( const int& id )
{
    VEdge2D* newEdge = new VEdge2D(id);
    m_VEdges.push_back(newEdge);
    return newEdge;
}



VFace2D* VoronoiDiagram2DC::createFace( const int& id )
{
    VFace2D* newFace = new VFace2D(id);
    m_VFaces.push_back(newFace);
    return newFace;
}



Generator2D * VoronoiDiagram2DC::createGenerator()
{
    Generator2D* generator = new Generator2D();
    m_generators.push_back(generator);
    return generator;
}



bool VoronoiDiagram2DC::isThisPhantomGenerator( const Generator2D * generator ) const
{
    bool b_isThisPhantomGenerator = false;

    list<const Generator2D*> phantomGenerators;
    getPhantomGenerator(phantomGenerators);

    list<const Generator2D*>::iterator i_phantom;
    for ( i_phantom = phantomGenerators.begin(); i_phantom != phantomGenerators.end(); ++i_phantom) {

        const Generator2D* aPhantom = *i_phantom;

        if (aPhantom == generator){
            b_isThisPhantomGenerator = true;
            break;
        }
    }

    return b_isThisPhantomGenerator;
}



VFace2D* VoronoiDiagram2DC::getVirtualFace() const
{
    return m_VFaces.front();
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
    for( list<Generator2D*>::iterator i_gen = m_generators.begin() ; i_gen != m_generators.end() ; ++i_gen ) {
        Generator2D* currGenerator = *i_gen;
        rg_Circle2D currDisk = currGenerator->getDisk();
        boundingBox.updateBoxByAddingCircle( currDisk );

        if( currDisk.getRadius() > maxRadius ) {
            maxRadius = currDisk.getRadius();
        }
    }
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
    double      expansionRatio                      = 3.1;
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
    rg_Circle2D phantomCircle[3] = { *it_phantom++, *it_phantom++, *it_phantom };

    Generator2D* phantom1 = new Generator2D();
    Generator2D* phantom2 = new Generator2D();
    Generator2D* phantom3 = new Generator2D();

    phantom1->setID(-3);
    phantom2->setID(-2);
    phantom3->setID(-1);

    phantom1->setDisk( phantomCircle[0] );
    phantom2->setDisk( phantomCircle[1] );
    phantom3->setDisk( phantomCircle[2] );

    m_phantomGenerators.push_back(phantom1);
    m_phantomGenerators.push_back(phantom2);
    m_phantomGenerators.push_back(phantom3);
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


    list<Generator2D*>::iterator it_phantoms = m_phantomGenerators.begin();
    Generator2D* phantoms[3] = {*it_phantoms++, *it_phantoms++, *it_phantoms};


    //  1. make elements of voronoi diagram of phantom generators

    VVertex2D* vertices[4]   = { createVertex(0),        // vertex defined by three phantom generators
        createVertex(1),        // infinite vertex on edge0
        createVertex(2),        // infinite vertex on edge1
        createVertex(3) };      // infinite vertex on edge2

    VEdge2D*   edges[6]      = { createEdge(0),            // edge made by phantom0 and phantom1     
        createEdge(1),            // edge made by phantom1 and phantom2
        createEdge(2),            // edge made by phantom2 and phantom0
        createEdge(3),            // artificial boundary for unbounded face1
        createEdge(4),            // artificial boundary for unbounded face2 
        createEdge(5) };          // artificial boundary for unbounded face3 

    VFace2D* faces[4]        = { createFace(0),            // infinite space defined by artificial boundary for voronoi diagram
        createFace(1),            // face containing phantom0 
        createFace(2),            // face containing phantom1
        createFace(3) };          // face containing phantom2



    //  2. set topology between elements of voronoi diagram of phantom generators

    phantoms[0]->setOuterFace(faces[1]);
    phantoms[1]->setOuterFace(faces[2]);
    phantoms[2]->setOuterFace(faces[3]);

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
    rg_Circle2D::makeCircumcircle( phantoms[0]->getDisk(), phantoms[1]->getDisk(), phantoms[2]->getDisk(), tangentCircle1, tangentCircle2);
    //rg_Point2D locationOfVertex0 = tangentCircle1.getCenterPt();
    //vertices[0]->setLocation(locationOfVertex0);
    vertices[0]->setCircumcircle(tangentCircle1);


    //      3.2. define coordinate of infinite vertices

    rg_Point2D centerOfPoint0And1 = ( phantoms[0]->getDisk().getCenterPt() + phantoms[1]->getDisk().getCenterPt() ) / 2.0;
    rg_Point2D centerOfPoint1And2 = ( phantoms[1]->getDisk().getCenterPt() + phantoms[2]->getDisk().getCenterPt() ) / 2.0;
    rg_Point2D centerOfPoint2And0 = ( phantoms[2]->getDisk().getCenterPt() + phantoms[0]->getDisk().getCenterPt() ) / 2.0;

    rg_Point2D locationOfVertex0 = tangentCircle1.getCenterPt();
    rg_Point2D vectorOfVertex1   = centerOfPoint0And1 - locationOfVertex0;
    rg_Point2D vectorOfVertex2   = centerOfPoint1And2 - locationOfVertex0;
    rg_Point2D vectorOfVertex3   = centerOfPoint2And0 - locationOfVertex0;

    rg_Point2D locationOfVertex1 = centerOfPoint0And1 + 5.0 * vectorOfVertex1;
    rg_Point2D locationOfVertex2 = centerOfPoint1And2 + 5.0 * vectorOfVertex2;
    rg_Point2D locationOfVertex3 = centerOfPoint2And0 + 5.0 * vectorOfVertex3;

    vertices[1]->setLocation(locationOfVertex1);
    vertices[2]->setLocation(locationOfVertex2);
    vertices[3]->setLocation(locationOfVertex3); 
}



void VoronoiDiagram2DC::updateVoronoiDiagramWithNewGenerator( Generator2D* const newGenerator )
{
    Generator2D* closestGenerator = findClosestGeneratorToNewGenerator( newGenerator );
    
    
    list<VEdge2D*>   intersectingVEdges;
    list<VVertex2D*> fictitiousVVertices;
    findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop_EdgeSplitPossible( newGenerator, closestGenerator, intersectingVEdges, fictitiousVVertices );
    

    list<VVertex2D*> newVVertices;
    list<VEdge2D*>   newVEdges;
    makeVEdgeLoopForNewGeneratorAndConnectToCurrVD( newGenerator, intersectingVEdges, newVVertices, newVEdges );
    
    
    
    connectCurrVDToNewVEdgeLoop( newVEdges );
    
    
    
    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );
    
    
    
    computeCoordOfNewVVertices( newVVertices );

	// compute tangent circle for contact map by Joonghyun March, 17
	updateTangentCircleForContactMap(newVEdges);
}



void VoronoiDiagram2DC::updateVoronoiDiagramWithNewGenerator_new_influenced_removed_edges( Generator2D * const newGenerator, list<VEdge2D*>& newEdges, list<VEdge2D*>& influencedEdges, list<VEdge2D*>& removedEdges )
{
    Generator2D* closestGenerator = findClosestGeneratorToNewGenerator( newGenerator );

    clock_t endTime = clock();

    list<VEdge2D*>   intersectingVEdges;
    list<VVertex2D*> fictitiousVVertices;
    findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop_EdgeSplitPossible_influenced_removed_edge( newGenerator, closestGenerator, intersectingVEdges, fictitiousVVertices, influencedEdges, removedEdges );
    
    list<VVertex2D*> newVVertices;
    list<VEdge2D*>   newVEdges;
    makeVEdgeLoopForNewGeneratorAndConnectToCurrVD( newGenerator, intersectingVEdges, newVVertices, newVEdges );

    connectCurrVDToNewVEdgeLoop( newVEdges );

    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );

    computeCoordOfNewVVertices( newVVertices );

    // compute tangent circle for contact map by Joonghyun March, 17
    updateTangentCircleForContactMap(newVEdges);

    VFace2D* newFace = newGenerator->getOuterFace();
    for ( list<VEdge2D*>::iterator i_edge = newVEdges.begin(); i_edge != newVEdges.end(); ++i_edge ) {
        newEdges.push_back(*i_edge);

        if ( (*i_edge)->getLeftFace() == newFace ) {
            influencedEdges.push_back((*i_edge)->getRightHand());
        }
        else {
            influencedEdges.push_back((*i_edge)->getLeftLeg());
        }
    }
}



void VoronoiDiagram2DC::updateVoronoiDiagramWithNewGenerator_for_initial_arrangement_of_diskpacking( Generator2D * const newGenerator, list<VVertex2D*>& redVertices, list<VEdge2D*>& redEdges, list<VVertex2D*>& newVertices, list<VEdge2D*>& newEdges, list<VEdge2D*>& influencedEdges )
{
    Generator2D* closestGenerator = findClosestGeneratorToNewGenerator( newGenerator );

    clock_t endTime = clock();

    list<VEdge2D*>   intersectingVEdges;
    list<VVertex2D*> fictitiousVVertices;
    findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop_red_entities( newGenerator, closestGenerator, intersectingVEdges, fictitiousVVertices, redVertices, redEdges );

    list<VVertex2D*> newVVertices;
    list<VEdge2D*>   newVEdges;
    makeVEdgeLoopForNewGeneratorAndConnectToCurrVD( newGenerator, intersectingVEdges, newVVertices, newVEdges );

    connectCurrVDToNewVEdgeLoop( newVEdges );

    mergeSplitVEdgesByFictitiousVVertex( fictitiousVVertices );

    computeCoordOfNewVVertices( newVVertices );

    // compute tangent circle for contact map by Joonghyun March, 17
    updateTangentCircleForContactMap(newVEdges);

    VFace2D* newFace = newGenerator->getOuterFace();
    for ( list<VEdge2D*>::iterator i_edge = newVEdges.begin(); i_edge != newVEdges.end(); ++i_edge ) {
        newEdges.push_back(*i_edge);

        if ( (*i_edge)->getLeftFace() == newFace ) {
            influencedEdges.push_back((*i_edge)->getRightHand());
        }
        else {
            influencedEdges.push_back((*i_edge)->getLeftLeg());
        }
    }

    for ( list<VVertex2D*>::iterator i_vtx = newVVertices.begin(); i_vtx != newVVertices.end(); ++i_vtx ) {
        VVertex2D* vertex = *i_vtx;
        newVertices.push_back(vertex);
    }
}



Generator2D* VoronoiDiagram2DC::findClosestGeneratorToNewGenerator( Generator2D* const newGenerator )
{
    VFace2D* closestVFace = findVFaceContainingQueryPoint( newGenerator->getDisk().getCenterPt() );
    return closestVFace->getGenerator();
}



void VoronoiDiagram2DC::removeAllExtraneousVVerticesAndVEdges()
{
    list<VVertex2D*>::iterator i_vtx = m_VVertices.begin();
    while ( i_vtx != m_VVertices.end() ) {
        VVertex2D* vertex = *i_vtx;
        if ( vertex->getStatus() == RED_V || vertex->getStatus() == YELLOW_V ) {
            i_vtx = m_VVertices.erase(i_vtx);
            delete vertex;
        }
        else {
            ++i_vtx;
        }
    }

    list<VEdge2D*>::iterator i_edge = m_VEdges.begin();
    while ( i_edge != m_VEdges.end() ) {
        VEdge2D* edge = *i_edge;
        if ( edge->getStatus() == RED_E ) {
            i_edge = m_VEdges.erase(i_edge);
            delete edge;
        }
        else {
            ++i_edge;
        }
    }
}



void VoronoiDiagram2DC::findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop_EdgeSplitPossible( Generator2D* const newGenerator, Generator2D* const closestGenerator, list<VEdge2D*>& intersectingVEdges, list<VVertex2D*>& fictitiousVVertices )
{
    list<VVertex2D*> redVVertices;
    VVertex2D*  firstRedVertex = findFirstRedVertex( newGenerator, closestGenerator );
       

    firstRedVertex->setStatus( RED_V );
    redVVertices.push_back( firstRedVertex );

    
    
    list<VVertex2D*> blueVVertices;
    wavePropagation_ver1( newGenerator, redVVertices, blueVVertices, fictitiousVVertices );
    
    
    
    findCrossingEdges( redVVertices, intersectingVEdges );
    
    
    
    reconfigureByConvertingBlueVVerticesToWhite( blueVVertices );
}



void VoronoiDiagram2DC::findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop_EdgeSplitPossible_influenced_removed_edge( Generator2D * const newGenerator, Generator2D * const closestGenerator, list<VEdge2D*>& intersectingVEdges, list<VVertex2D*>& fictitiousVVertices, list<VEdge2D*>& influencedEdges, list<VEdge2D*>& removedEdges )
{
    list<VVertex2D*> redVVertices;
    VVertex2D*  firstRedVertex = findFirstRedVertex( newGenerator, closestGenerator );

    firstRedVertex->setStatus( RED_V );
    redVVertices.push_back( firstRedVertex );

    list<VVertex2D*> blueVVertices;
    wavePropagation_ver1( newGenerator, redVVertices, blueVVertices, fictitiousVVertices );

    //findCrossingEdges( redVVertices, intersectingVEdges );
    findCrossingEdges_red_edge( redVVertices, intersectingVEdges, removedEdges );

    reconfigureByConvertingBlueVVerticesToWhite( blueVVertices );
}



void VoronoiDiagram2DC::findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop_red_entities( Generator2D * const newGenerator, Generator2D * const closestGenerator, list<VEdge2D*>& intersectingVEdges, list<VVertex2D*>& fictitiousVVertices, list<VVertex2D*>& redVertices, list<VEdge2D*>& redEdges )
{
    list<VVertex2D*> redVVertices;
    VVertex2D*  firstRedVertex = findFirstRedVertex( newGenerator, closestGenerator );

    firstRedVertex->setStatus( RED_V );
    redVVertices.push_back( firstRedVertex );

    list<VVertex2D*> blueVVertices;
    wavePropagation_ver1( newGenerator, redVVertices, blueVVertices, fictitiousVVertices );

    findCrossingEdges_red_edge( redVVertices, intersectingVEdges, redEdges );

    reconfigureByConvertingBlueVVerticesToWhite( blueVVertices );

    for ( list<VVertex2D*>::iterator i_redVtx = redVVertices.begin(); i_redVtx != redVVertices.end(); ++i_redVtx ) {
        VVertex2D* vertex = *i_redVtx;
        redVertices.push_back(vertex);
    }
}



void VoronoiDiagram2DC::findAnomalizingEdgeAmongIncidentVEdgesAndSplit( VVertex2D* const currRedVVertex, Generator2D* const newGenerator, list<VVertex2D*>& fictitiousVVertices, list<VEdge2D*>& anomalyTestedVEdges )
{
    VEdge2D* beginningEdgeOfVertex  = currRedVVertex->getFirstVEdge();
    VEdge2D* nextIncidentEdge       = beginningEdgeOfVertex;
    VEdge2D* currIncidentEdge       = NULL;

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



void VoronoiDiagram2DC::findRedVVertexAmongNeighbor( VVertex2D* const currRedVVertex, Generator2D* const newGenerator, list<VVertex2D*>& blueVVertices, list<VVertex2D*>& redVVertices, list<VVertex2D*>& fictitiousVVertices, list<VEdge2D*>& anomalyTestedVEdges )
{
    VVertex2D* currNeighbor             = NULL;
    VEdge2D*   beginningEdgeOfVertex    = currRedVVertex->getFirstVEdge();
    VEdge2D*   nextIncidentEdge         = beginningEdgeOfVertex;
    VEdge2D*   currIncidentEdge         = NULL;

    do 
    {
        currIncidentEdge = nextIncidentEdge;
        nextIncidentEdge = getCCWNextEdgeOnVertex( currIncidentEdge, currRedVVertex );

        currNeighbor = currIncidentEdge->getOppositVVertex( currRedVVertex );

        if( currNeighbor->getStatus() != WHITE_V ) {
            continue;
        }

        double MUValue = computeMUValue( currNeighbor, newGenerator );

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



void VoronoiDiagram2DC::findCrossingEdges_red_edge( list<VVertex2D*>& redVVertices, list<VEdge2D*>& intersectingVEdges, list<VEdge2D*>& redEdges )
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
                redEdges.push_back( currIncidentVEdge );
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
        VVertex2D* currNewVVertex = createVertex( m_VVertices.back()->getID()+1 );
        //VVertex2D* currNewVVertex = createVVertex( VVertex2D(m_VVertices.rbegin()->getID()+1) );
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
        

        VEdge2D* newEdge = createEdge( m_VEdges.back()->getID() + 1 );
        //VEdge2D* newEdge = createVEdge( VEdge2D(m_VEdges.rbegin()->getID() + 1) );
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
    VFace2D* newVFace = createFace( m_VFaces.back()->getID() + 1 );
    //VFace2D* newVFace = createVFace( VFace2D(m_VFaces.rbegin()->getID() + 1) );
    newGenerator->setOuterFace( newVFace );
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
        numCircumcircles = rg_Circle2D::makeCircumcircle( definingGenerator[0].getDisk(),
            definingGenerator[1].getDisk(),
            definingGenerator[2].getDisk(),
            firstCircle,
            secondCircle );


        switch( numCircumcircles )
        {

        case 0:
            cout << "Cannot compute coordinate of newVVertex!" << endl;
            exit(1);
            break;


        case 1:
            //currVVertex->setLocation( firstCircle.getCenterPt() );
            currVVertex->setCircumcircle( firstCircle );
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

                // ???ÔøΩÏùò ?ÔøΩÏù¥ Í∞ôÎã§Ôø??
                if( currIncidentEdge->getStartVertex()->getStatus() == currIncidentEdge->getEndVertex()->getStatus() ) {
                    continue;
                }
                // ?ÔøΩÎ•¥?ÔøΩÎ©¥ 
                else {
                    // Í∞ôÔøΩ? 
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
                centerOfReferenceGenerator = currIncidentEdge->getRightFace()->getGenerator()->getDisk().getCenterPt();
                baseLine = currIncidentEdge->getEndVertex()->getLocation() - centerOfReferenceGenerator;
            }
            else {
                centerOfReferenceGenerator = currIncidentEdge->getLeftFace()->getGenerator()->getDisk().getCenterPt();
                baseLine = currIncidentEdge->getStartVertex()->getLocation() - centerOfReferenceGenerator;
            }

            firstLine = firstCircle.getCenterPt() - centerOfReferenceGenerator;
            secondLine = secondCircle.getCenterPt() - centerOfReferenceGenerator;

            rg_REAL firstAngle = angleFromVec1toVec2( baseLine, firstLine );
            rg_REAL secondAngle = angleFromVec1toVec2( baseLine, secondLine );

            if( firstAngle < secondAngle ) {
                //currVVertex->setLocation( firstCircle.getCenterPt() );
                currVVertex->setCircumcircle( firstCircle );
            }
            else {
                //currVVertex->setLocation( secondCircle.getCenterPt() );
                currVVertex->setCircumcircle( secondCircle );
            }

            currVVertex->setStatus( WHITE_V );
            break;

        }

    }
}



double VoronoiDiagram2DC::computeMUValue( VVertex2D * const vertex, Generator2D * const newGenerator )
{
    list<Generator2D*> definingGenerators;
    vertex->getDefining3Generators( definingGenerators );

    if ( definingGenerators.size() < 3 ) {
        return DBL_MAX;
    }

    Generator2D* definingGenerator = definingGenerators.front();
    double       rho = vertex->getCircumcircle().getRadius();
    double       delta = 0.0;
    if ( newGenerator->getType() == Generator2D::DISK_G ) {
        delta = newGenerator->getDisk().getCenterPt().distance( vertex->getCircumcircle().getCenterPt() ) - newGenerator->getDisk().getRadius();
    }
    else if ( newGenerator->getType() == Generator2D::VERTEX_G ) {
        delta = ( (VertexGenerator2D*)newGenerator )->get_point().distance( vertex->getCircumcircle().getCenterPt() );
    }
    double       rhoSquare = rho * rho;

    double mu = ( delta - rho ) / rhoSquare;

    return mu;
}



VVertex2D* VoronoiDiagram2DC::findFirstRedVertex( Generator2D* const newGenerator, Generator2D* const closestGenerator )
{
    list<VVertex2D*> boundaryVVerticesOfFaceOfClosestGenerator;
    if ( (VoronoiDiagram2DC*)(closestGenerator->getInnerVD()) == this ) 
        closestGenerator->getInnerFace()->getBoundaryVVertices( boundaryVVerticesOfFaceOfClosestGenerator );
    else
        closestGenerator->getOuterFace()->getBoundaryVVertices( boundaryVVerticesOfFaceOfClosestGenerator );

    double      minMUValue = DBL_MAX;
    VVertex2D*  firstRedVertex = NULL;

    for( list<VVertex2D*>::iterator it_boundaryVertex = boundaryVVerticesOfFaceOfClosestGenerator.begin() ; it_boundaryVertex != boundaryVVerticesOfFaceOfClosestGenerator.end() ; it_boundaryVertex++ ) {
        VVertex2D* currVVertex = *it_boundaryVertex;

        double currMUValue = computeMUValue( currVVertex, newGenerator );
        if( currMUValue < minMUValue ) {
            minMUValue      = currMUValue;
            firstRedVertex  = currVVertex;
        }
    }

    if( minMUValue >= 0) {
        list<VEdge2D*> boundaryVEdges;

        if ( (VoronoiDiagram2DC*)(closestGenerator->getInnerVD()) == this ) 
            closestGenerator->getInnerFace()->getBoundaryVEdges( boundaryVEdges );
        else
            closestGenerator->getOuterFace()->getBoundaryVEdges( boundaryVEdges );

        for( list<VEdge2D*>::iterator it_boundaryEdge = boundaryVEdges.begin() ; it_boundaryEdge != boundaryVEdges.end() ; it_boundaryEdge++ ) {
            VEdge2D* currVEdge = *it_boundaryEdge;
            if( this->isAnomalizingEdge( currVEdge, newGenerator ) ) {
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

    VVertex2D* currNeighbor             = NULL;
    VEdge2D*   beginningEdgeOfVertex    = candidateForRedVVertex->getFirstVEdge();
    VEdge2D*   nextIncidentEdge         = beginningEdgeOfVertex;
    VEdge2D*   currIncidentEdge         = NULL;

    do 
    {
        currIncidentEdge = nextIncidentEdge;
        nextIncidentEdge = getCCWNextEdgeOnVertex( currIncidentEdge, candidateForRedVVertex );

        currNeighbor = currIncidentEdge->getOppositVVertex( candidateForRedVVertex );

        if( currNeighbor->getStatus() == RED_V ) {
            numRedVVertexAmongNeighbor++;
        }
    }while( nextIncidentEdge != beginningEdgeOfVertex );

    if( numRedVVertexAmongNeighbor > 1) {
        return true;
    }
    else {
        return false;
    }
}



bool VoronoiDiagram2DC::hasOldRegionSplit( VVertex2D* const candidateForRedVVertex )
{
    VEdge2D* beginningEdgeOfVertex  = candidateForRedVVertex->getFirstVEdge();
    VFace2D* currIncidentFace       = NULL;
    VEdge2D* currIncidentEdge       = NULL;
    VEdge2D* nextIncidentEdge       = beginningEdgeOfVertex;

    do 
    {
        currIncidentEdge            = nextIncidentEdge;
        nextIncidentEdge            = getCCWNextEdgeOnVertex( currIncidentEdge, candidateForRedVVertex );

        if( candidateForRedVVertex == currIncidentEdge->getStartVertex() ) {
            currIncidentFace        = currIncidentEdge->getLeftFace();
        }
        else {
            currIncidentFace        = currIncidentEdge->getRightFace();
        }

        // would-be red VVertexÔø???ÔøΩÏãúÔø??redÔø??ÎßåÎì†??
        candidateForRedVVertex->setStatus( RED_V );

        int     numComponentsOfBoundaryVVertices = 0;
        bool    isCurrVVertexRed    = false;
        bool    isPrevVVertexRed    = false;

        VEdge2D*   beginningEdgeOfIncidentFace = currIncidentFace->getFirstVEdge();
        VVertex2D* firstBoundaryVertex = NULL;

        VVertex2D* currBoundaryVertex   = NULL;
        VEdge2D*   currBoundaryEdge     = NULL;
        VEdge2D*   nextBoundaryEdge     = beginningEdgeOfIncidentFace;

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

        if( numComponentsOfBoundaryVVertices % 2 == 0 ) {
            continue;
        }
        else {
            return true;
        }

    } while ( nextIncidentEdge != beginningEdgeOfVertex );

    return false;
}



bool VoronoiDiagram2DC::isAnomalizingEdge( VEdge2D* const incidentEdge, Generator2D* const newGenerator )
{
    if( incidentEdge->getStartVertex()->isFictitious() || incidentEdge->getEndVertex()->isFictitious() ) {
        return false;
    }

    rg_Circle2D disk[3];
    
    disk[0] = newGenerator->getDisk();
    disk[1] = incidentEdge->getLeftFace()->getGenerator()->getDisk();
    disk[2] = incidentEdge->getRightFace()->getGenerator()->getDisk();

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

    bool thisEdgeIsAnAnomalyzingEdge = false;

    switch (numCircumcircles) {
    case 2:
        if( incidentEdge->getLeftFace()->getGenerator()->getDisk().getRadius() < incidentEdge->getRightFace()->getGenerator()->getDisk().getRadius() ) {
            center0         = incidentEdge->getRightFace()->getGenerator()->getDisk().getCenterPt();
            CWBoundary      = incidentEdge->getEndVertex()->getLocation() - center0;
            CCWBoundary     = incidentEdge->getStartVertex()->getLocation() - center0;
        }
        else {
            center0         = incidentEdge->getLeftFace()->getGenerator()->getDisk().getCenterPt() ;
            CWBoundary      = incidentEdge->getStartVertex()->getLocation() - center0;
            CCWBoundary     = incidentEdge->getEndVertex()->getLocation() - center0;
        }

        firstLocation       = firstCircle.getCenterPt();
        secondLocation      = secondCircle.getCenterPt();
        firstLocation       = firstLocation - center0;
        secondLocation      = secondLocation - center0;

        CWBoundary3D        = rg_Point3D( CWBoundary );
        CCWBoundary3D       = rg_Point3D( CCWBoundary );
        firstLocation3D     = rg_Point3D( firstLocation );
        secondLocation3D    = rg_Point3D( secondLocation );
        
        CW_First_Test       = CWBoundary3D.crossProduct( firstLocation3D );
        CCW_First_Test      = CCWBoundary3D.crossProduct( firstLocation3D );
        CW_Second_Test      = CWBoundary3D.crossProduct( secondLocation3D );
        CCW_Second_Test     = CCWBoundary3D.crossProduct( secondLocation3D );

        if( CW_First_Test.getZ() < 0 ) {
            break;
        }
        if( CW_Second_Test.getZ() < 0 ) {
            break;
        }
        if( CCW_First_Test.getZ() > 0 ) {
            break;
        }
        if( CCW_Second_Test.getZ() > 0 ) {
            break;
        }

        thisEdgeIsAnAnomalyzingEdge = true;

    default:
        break;
    }

    return thisEdgeIsAnAnomalyzingEdge;
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

        currFictitiousVVertex->setStatus( RED_V );
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
	return rg_GT(circle1.getRadius(), circle2.getRadius());    
}



bool VoronoiDiagram2DC::compareCircleAscendingOrder(const rg_Circle2D& circle1, const rg_Circle2D& circle2)
{
    return rg_LT(circle1.getRadius(), circle2.getRadius());
}



bool VoronoiDiagram2DC::compareCirclePairDescendingorder( const pair< rg_Circle2D, void* >& circle1, const pair< rg_Circle2D, void* >& circle2 )
{
    return rg_GT(circle1.first.getRadius(), circle2.first.getRadius());    
}



bool VoronoiDiagram2DC::compareCirclePairAscendingOrder( const pair<rg_Circle2D, void*>& circle1, const pair<rg_Circle2D, void*>& circle2 )
{
    return rg_LT( circle1.first.getRadius(), circle2.first.getRadius() );
}



bool VoronoiDiagram2DC::compareCircleTripletDescendingorder( const rg_Triplet< rg_Circle2D, int, void* >& circle1, const rg_Triplet< rg_Circle2D, int, void* >& circle2 )
{
    return rg_GT(circle1.m_first.getRadius(), circle2.m_first.getRadius());    
}



bool VoronoiDiagram2DC::compareGeneratorNonAscendingOrder(const Generator2D* const generator1, const Generator2D* const generator2)
{
	return rg_GT(generator1->getDisk().getRadius(), generator2->getDisk().getRadius());
}



void VoronoiDiagram2DC::getPhantomGenerator( list<const Generator2D*>& phantomGeneratorsList ) const
{
    for( list<Generator2D*>::const_iterator it_generator = m_phantomGenerators.begin() ; it_generator != m_phantomGenerators.end() ; ++it_generator ) {
        const Generator2D* currGenerator = *it_generator;
        phantomGeneratorsList.push_back( currGenerator );
    }
}


void VoronoiDiagram2DC::getPhantomGenerator(list<Generator2D*>& phantomGeneratorsList) const
{
    for (list<Generator2D*>::const_iterator it_generator = m_phantomGenerators.begin(); it_generator != m_phantomGenerators.end(); ++it_generator) {
        const Generator2D* currGenerator = *it_generator;
        phantomGeneratorsList.push_back(const_cast<Generator2D*>(currGenerator));
    }
}



VVertex2D* VoronoiDiagram2DC::splitVEdgeAtFictitiousVVertex( VEdge2D* const currVEdge )
{
    VVertex2D* fictitiousVVertex = createVertex( m_VVertices.back()->getID()+1 );
    fictitiousVVertex->setStatus( YELLOW_V );
    fictitiousVVertex->setFirstVEdge( currVEdge );
    fictitiousVVertex->setFictitious();

    VEdge2D* dividedVEdge = createEdge( m_VEdges.back()->getID() + 1 );

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

    if( rg_EQ( target->getLeftFace()->getGenerator()->getDisk().getRadius(),
        target->getRightFace()->getGenerator()->getDisk().getRadius(), resNeg6 ) )
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
        passPoint = this->getPassingPtOnEdge( target->getLeftFace(), 
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
//     rg_Point2D sp, vec1, vec2;
//     sp = tVertex->getLocation();
//     vec1 = (lFace->getGenerator()->getDisk()->getCenterPt() - sp).getUnitVector();
//     vec2 = (rFace->getGenerator()->getDisk()->getCenterPt() - sp).getUnitVector();
// 
//     if( rg_ZERO( (vec1+vec2).magnitude() ) )
//     {
//         return rg_Point2D( -vec1.getY(), vec1.getX() );
//     }
//     else
//     {
//         return vec1 + vec2;
//     }

    rg_Point2D sp, vec1, vec2;
    sp = tVertex->getLocation();
    vec1 = (lFace->getGenerator()->getDisk().getCenterPt() - sp).getUnitVector();
    vec2 = (rFace->getGenerator()->getDisk().getCenterPt() - sp).getUnitVector();

    rg_Point2D tangentVector;

    if ( (VoronoiDiagram2DC*)lFace->getGenerator()->getInnerVD() == this ) {
        tangentVector = vec2 - vec1;
    }
    else if ( (VoronoiDiagram2DC*)rFace->getGenerator()->getInnerVD() == this ) {
        tangentVector = vec1 - vec2;
    }
    else {
        tangentVector = vec1 + vec2;
    }

    if( rg_ZERO( tangentVector.magnitude() ) ) {
        // Control pointÔø??Ï∞æÍ∏∞ ?ÔøΩÌï¥ tangent vectorÔø??Í≥ÑÏÇ∞???ÔøΩÎäî (+) (-)Ôø???ÔøΩÌôï??Í≥ÑÏÇ∞???ÔøΩÏöî???ÔøΩÎã§.
        tangentVector = rg_Point2D( -vec1.getY(), vec1.getX() );
    }

    return tangentVector;
}



rg_Point2D VoronoiDiagram2DC::getPassingPtOnEdge( VFace2D* lFace, VFace2D* rFace, VVertex2D* sVertex, VVertex2D* eVertex )
{
    rg_Circle2D c1, c2;	//c1.r > c2.r

    if( lFace->getGenerator()->getDisk().getRadius() > rFace->getGenerator()->getDisk().getRadius() )
    {
        c1 = lFace->getGenerator()->getDisk();
        c2 = rFace->getGenerator()->getDisk();
    }
    else
    {
        c2 = lFace->getGenerator()->getDisk();
        c1 = rFace->getGenerator()->getDisk();
    }

    rg_Point2D cp = ( c1.getCenterPt() + c2.getCenterPt() ) / 2.;

    return cp + ( c2.getCenterPt() - c1.getCenterPt() ).getUnitVector() * ( c1.getRadius() / 2.0 - c2.getRadius() / 2.0 );
}



void VoronoiDiagram2DC::removeREDVVerticesAndVEdges()
{
    list<VVertex2D*>::iterator i_vtx = m_VVertices.begin();
    while (i_vtx != m_VVertices.end()) {
        VVertex2D* vertex = *i_vtx;
        if (vertex->getStatus() == RED_V || vertex->getStatus() == YELLOW_V) {
            i_vtx = m_VVertices.erase(i_vtx);
            delete vertex;
        }
        else {
            ++i_vtx;
        }
    }

    list<VEdge2D*>::iterator i_edge = m_VEdges.begin();
    while (i_edge != m_VEdges.end()) {
        VEdge2D* edge = *i_edge;
        if (edge->getStatus() == RED_E) {
            i_edge = m_VEdges.erase(i_edge);
            delete edge;
        }
        else {
            ++i_edge;
        }
    }
}


void VoronoiDiagram2DC::updateGeometry()
{
    list<VEdge2D*>::iterator it_edges;
    for( it_edges = m_VEdges.begin() ; it_edges != m_VEdges.end() ; ++it_edges ) {
        VEdge2D* currEdge = *it_edges;

        if( currEdge->isInfinite() ) {
            continue;
        }

        updateEdge(currEdge);
    }
}



void VoronoiDiagram2DC::updateGeometryOfInsertedGenerator(Generator2D* insertedGenerator) 
{
    list<VEdge2D*> boundaryEdges;
    insertedGenerator->getOuterFace()->getBoundaryVEdges(boundaryEdges);
    list<VEdge2D*> quillEdges;
    insertedGenerator->getOuterFace()->getQuillVEdges(quillEdges);
    for ( list<VEdge2D*>::iterator i_edge = boundaryEdges.begin(); i_edge != boundaryEdges.end(); ++i_edge ) {
        VEdge2D* currEdge = *i_edge;
        updateEdge(currEdge);
    }
    for ( list<VEdge2D*>::iterator i_edge = quillEdges.begin(); i_edge != quillEdges.end(); ++i_edge ) {
        VEdge2D* currEdge = *i_edge;
        updateEdge(currEdge);
    }
}

void VoronoiDiagram2DC::updateGeometryOfDeletedGenerator(Generator2D* removedGenerator, list<VEdge2D*>& influencedEdges) 
{
    for (VEdge2D* currEdge : influencedEdges){
        updateEdge(currEdge);
    }
    // list<VEdge2D*> boundaryEdges;
    // removedGenerator->getOuterFace()->getBoundaryVEdges(boundaryEdges);
    // list<VEdge2D*> quillEdges;
    // removedGenerator->getOuterFace()->getQuillVEdges(quillEdges);

    // for ( list<VEdge2D*>::iterator i_edge = boundaryEdges.begin(); i_edge != boundaryEdges.end(); ++i_edge ) {
    //     VEdge2D* currEdge = *i_edge;
    //     updateEdge(currEdge);
    // }
    // for ( list<VEdge2D*>::iterator i_edge = quillEdges.begin(); i_edge != quillEdges.end(); ++i_edge ) {
    //     VEdge2D* currEdge = *i_edge;
    //     updateEdge(currEdge);
    // }
}



void VoronoiDiagram2DC::wavePropagation_ver1( Generator2D* newGenerator, list<VVertex2D*>& redVVertices, list<VVertex2D*>& blueVVertices, list<VVertex2D*>& fictitiousVVertices )
{
    list<VEdge2D*> anomalyTestedVEdges;
    list<VVertex2D*>::iterator it_redVVertex = redVVertices.begin();

    for( ; it_redVVertex != redVVertices.end() ; it_redVVertex++ ) {
        VVertex2D* currRedVVertex = *it_redVVertex;

        findAnomalizingEdgeAmongIncidentVEdgesAndSplit( currRedVVertex, newGenerator, fictitiousVVertices, anomalyTestedVEdges );

        findRedVVertexAmongNeighbor( currRedVVertex, newGenerator, blueVVertices, redVVertices, fictitiousVVertices, anomalyTestedVEdges );
    }

    resetAnomalyTestDoneToggle( anomalyTestedVEdges );
}



void VoronoiDiagram2DC::wavePropagation_ver2( Generator2D* newGenerator, list<VVertex2D*>& redVVertices, list<VVertex2D*>& blueVVertices, list<VVertex2D*>& fictitiousVVertices, list<VEdge2D*>& anomalyTestedVEdges )
{
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

                double MUValue = computeMUValue( couldBeRedVertex, newGenerator );

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

                if( hasOldRegionSplit( couldBeRedVertex ) ) { // split ?ÔøΩÏñ¥?ÔøΩÎäî ?ÔøΩÏàò count   
                    couldBeRedVertex->setStatus( BLUE_V );
                    blueVVertices.push_back( couldBeRedVertex );
                    continue;
                }

                couldBeRedVertex->setStatus( RED_V );
                redVVertices.push_back( couldBeRedVertex );
            }
        }
    }
}

void VoronoiDiagram2DC::removePhantomGenerators()
{
    VFace2D* virtualFace = m_VFaces.front();
    set<VFace2D*> phantomFaces; 

    list<VEdge2D*> virtualEdges;
    list<VEdge2D*> unboundedEdges_before_phantom_removal;
    // unboundedEdges_before_phantom_removal
    // AREdges, BridgeEdges, 
    list<VEdge2D*> edgesDefinedByPhantomNInputDisk;


    collectEntitiesInfluencedByPhantomRemoval( phantomFaces, virtualEdges, unboundedEdges_before_phantom_removal, edgesDefinedByPhantomNInputDisk );

    putRedTagsOnEntitiesTobeRemoved( virtualEdges, unboundedEdges_before_phantom_removal );

    adjustTopologyOfInterfaceBetweenInputDisksNPhantomDisks( virtualFace, unboundedEdges_before_phantom_removal, phantomFaces, edgesDefinedByPhantomNInputDisk );

    removePhantomFaces( phantomFaces );

    makeCorrectTopologyOfTheInterfaceByFlippingEdges( virtualFace );

    assignCoordinateOfInfiniteVertices();
}



void VoronoiDiagram2DC::collectEntitiesInfluencedByPhantomRemoval( set<VFace2D*>& phantomFaces, list<VEdge2D*>& virtualEdges, list<VEdge2D*>& unboundedEdges_before_phantom_removal, list<VEdge2D*>& edgesDefinedByPhantomNInputDisk )
{
    const int numDummyCircles = 3;
    list<VFace2D*>::iterator i_face = m_VFaces.begin();
    ++i_face;

    for( int i = 0 ; i < numDummyCircles ; ++i, ++i_face ) {
        VFace2D* currPhantomRegion = *i_face;
        phantomFaces.insert( currPhantomRegion );
    }

    list<VEdge2D*>::iterator i_edge = m_VEdges.begin();
    for( ; i_edge != m_VEdges.end() ; ++i_edge ) {
        VEdge2D* currEdge = *i_edge;
        if( currEdge->getStatus() == RED_E ) {
            continue;
        }

        if( currEdge->isInfinite() ) {
            virtualEdges.push_back( currEdge );
        }
        else if( currEdge->isUnBounded() ) {
            unboundedEdges_before_phantom_removal.push_back( currEdge );
        }
        else if( (currEdge->getLeftFace()->isUnBounded() && currEdge->getRightFace()->isBounded()) ||
            (currEdge->getLeftFace()->isBounded() && currEdge->getRightFace()->isUnBounded()) ) {
                edgesDefinedByPhantomNInputDisk.push_back( currEdge );
        }
    }
}



void VoronoiDiagram2DC::putRedTagsOnEntitiesTobeRemoved( list<VEdge2D*>& virtualEdges, list<VEdge2D*>& unboundedEdges_before_phantom_removal )
{
    list<VEdge2D*>::iterator it_infiniteEdge = virtualEdges.begin();
    for( ; it_infiniteEdge != virtualEdges.end() ; ++it_infiniteEdge ) {
        VEdge2D* currInfiniteEdge = *it_infiniteEdge;
        currInfiniteEdge->setStatus( RED_E );
    }

    list<VEdge2D*>::iterator it_unboundedEdge = unboundedEdges_before_phantom_removal.begin();
    for( ; it_unboundedEdge != unboundedEdges_before_phantom_removal.end() ; ++it_unboundedEdge ) {
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
    for( ; it_unboundedEdge != unboundedEdges_before_phantom_removal.end() ; ++it_unboundedEdge ) {
        VEdge2D* currUnboundedEdge = *it_unboundedEdge;
        VEdge2D* mergedEdge = createEdge( m_VEdges.back()->getID() + 1 );

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
    for( ; it_outsideBoundingEdge != edgesDefinedByPhantomNInputDisk.end() ; ++it_outsideBoundingEdge ) {
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
    for( ; it_phantomFace != phantomFaces.end() ; ++it_phantomFace ) {
        VFace2D* phantomFace = *it_phantomFace;
        m_VFaces.remove( phantomFace );
        delete phantomFace;
    }
}



void VoronoiDiagram2DC::correctTopologyInUnboundedRegion()
{
    // fix faces without the first edge.
    for( list<VEdge2D*>::iterator i_edge = m_VEdges.begin(); i_edge != m_VEdges.end() ; ++i_edge ) {
        VEdge2D* currEdge = *i_edge;

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
    list<VFace2D*>::iterator i_face = m_VFaces.begin();
    while( i_face != m_VFaces.end() && dummy_i < numDummyCircles ) {
        VFace2D* currFace = *i_face;

        if( currFace->isUnBounded() ) {
            phantomRegions.insert( currFace );
            ++dummy_i;
        }

        ++i_face;
    }

    // collect edges defined by phantom generators
    list<VEdge2D*> edgesDefinedByOnlyPhantom;

    set<VFace2D*>::iterator NOT_PHANTOM = phantomRegions.end();

    list<VEdge2D*>::iterator i_edge = m_VEdges.begin();
    for( ; i_edge != m_VEdges.end() ; ++i_edge ) {
        VEdge2D* currEdge = *i_edge;

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
    for( ; it_unboundedEdge != edgesDefinedByOnlyPhantom.end() ; ++it_unboundedEdge ) {
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

        m_VEdges.remove(unboundedEdge);
        m_VVertices.remove(infiniteVertex);

        delete unboundedEdge;
        delete infiniteVertex;
    }
}

void VoronoiDiagram2DC::replaceRegionsOfPhantomGeneratorsWithVirtualRegion()
{
    set<VFace2D*> phantomRegions;
    // find V-region of phantom generator
    const int numDummyCircles = 3;
    list<VFace2D*>::iterator it_face = m_VFaces.begin();
    VFace2D* virtualFace = *it_face;
    it_face++;

    for( int i = 0 ; i < numDummyCircles ; ++i, ++it_face ) {
        VFace2D* currPhantomRegion = *it_face;
        phantomRegions.insert( currPhantomRegion );
    }

    // collect edges defined by phantom generators
    set<VFace2D*>::iterator NOT_PHANTOM = phantomRegions.end();

    list<VEdge2D*>::iterator it_edge = m_VEdges.begin();
    for( ; it_edge != m_VEdges.end() ; ++it_edge ) {
        VEdge2D* currEdge = *it_edge;

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
        edgeToDefineMinTangentCircle.first->isCandidateToBeFlippedInPhantomRemoval(false);

        VEdge2D* edgeOfIncidentTriplet[2] = {NULL, NULL};
        flipThisEdgeNFindTwoIncidentTriplets( virtualFace, edgeToDefineMinTangentCircle, edgeOfIncidentTriplet[0], edgeOfIncidentTriplet[1] );
        reflectTheFlipOnTheIncidentTriplet( possiblyFlippingEdges, edgeOfIncidentTriplet[0], virtualFace );
        reflectTheFlipOnTheIncidentTriplet( possiblyFlippingEdges, edgeOfIncidentTriplet[1], virtualFace );
    }
}


void VoronoiDiagram2DC::assignCoordinateOfInfiniteVertices()
{
    rg_BoundingBox2D boundingBox( m_bucket.getMinPt(), m_bucket.getMaxPt() );

    VFace2D* virtualFace = getVirtualFace();
    list<VVertex2D*> infiniteVertices;
    virtualFace->getBoundaryVVertices( infiniteVertices );

    for( list<VVertex2D*>::iterator i_vtx = infiniteVertices.begin(); i_vtx != infiniteVertices.end(); ++i_vtx ) {
        VVertex2D* currVtx = *i_vtx;
        boundingBox.contain( currVtx->getLocation() );
    }

    rg_Point2D center = boundingBox.getCenterPt();
    double     radius = 1.5*boundingBox.evaluateLongestLength();

    m_circleEnclosingVertices.setCenterPt(center);
    m_circleEnclosingVertices.setRadius(radius);





    // update geometry of unbounded edge
    list<VEdge2D*> unboundedEdges;
    virtualFace->getQuillVEdges( unboundedEdges );

    list<VEdge2D*>::iterator it_unboundedEdge = unboundedEdges.begin();
    for( ; it_unboundedEdge != unboundedEdges.end() ; ++it_unboundedEdge ) {
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

        rg_Circle2D leftCircle  = currEdge->getLeftFace()->getGenerator()->getDisk();
        rg_Circle2D rightCircle = currEdge->getRightFace()->getGenerator()->getDisk();

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



// Previous version which uses linked list as a priority queue
/*
void VoronoiDiagram2DC::collectPossiblyFlippingEdgesOnTheInterface_StoredInPriorityQueue( VFace2D* virtualFace, rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges )
{
    list<VEdge2D*> edgesToBoundVirtualFace;
    virtualFace->getBoundaryVEdges( edgesToBoundVirtualFace );

    list<VEdge2D*>::iterator it_boundingEdge = edgesToBoundVirtualFace.begin();
    for( ; it_boundingEdge != edgesToBoundVirtualFace.end() ; ++it_boundingEdge ) {
        VEdge2D* currEdge = *it_boundingEdge;

        rg_Circle2D circumcircle;
        bool isFlippingEdge = isFlippingEdgeOnVirtualRegion( currEdge, circumcircle, virtualFace );

        if( isFlippingEdge ) {
            possiblyFlippingEdges.add( pair<VEdge2D*, rg_Circle2D>(currEdge, circumcircle) );
        }
    }

}
*/


void VoronoiDiagram2DC::collectPossiblyFlippingEdgesOnTheInterface_StoredInPriorityQueue( VFace2D * virtualFace, Priority_Q_VEdge & possiblyFlippingEdges )
{
    list<VEdge2D*> edgesToBoundVirtualFace;
    virtualFace->getBoundaryVEdges( edgesToBoundVirtualFace );

    
    for( list<VEdge2D*>::iterator i_bndEdge = edgesToBoundVirtualFace.begin(); i_bndEdge != edgesToBoundVirtualFace.end() ; ++i_bndEdge ) {
        VEdge2D* currEdge = *i_bndEdge;

        rg_Circle2D circumcircle;
        bool isFlippingEdge = isFlippingEdgeOnVirtualRegion( currEdge, circumcircle, virtualFace );

        if( isFlippingEdge ) {
            possiblyFlippingEdges.push( pair<VEdge2D*, rg_Circle2D>(currEdge, circumcircle), circumcircle.getRadius() );
            currEdge->isCandidateToBeFlippedInPhantomRemoval(true);
        }
    }
}



bool VoronoiDiagram2DC::isFlippingEdgeOnVirtualRegion( VEdge2D* edge, rg_Circle2D& circumcircle, VFace2D* virtualFace )
{
    //VFace2D* virtualFace = &(*m_VFaces.begin());

    VFace2D* currFace = rg_NULL;
    VFace2D* nextFace = rg_NULL;
    VFace2D* prevFace = rg_NULL;

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

    if( currFace->getGenerator() == rg_NULL || prevFace->getGenerator() == rg_NULL || nextFace->getGenerator() == rg_NULL ) {
        return isFlippingEdge;
    }
    if ( nextFace == virtualFace || prevFace == virtualFace ) {
        return isFlippingEdge;
    }
    if ( nextFace == prevFace || nextFace == currFace || prevFace == currFace ) {
        return isFlippingEdge;
    }

    rg_Circle2D tangentCircle[2];
    rg_INT numCC = computeCircumcirclesForDeletion( prevFace, currFace, nextFace, tangentCircle[0], tangentCircle[1] );
    

    // else,
    // ...

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
    }

    return isFlippingEdge;
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



void VoronoiDiagram2DC::flipThisEdgeNFindTwoIncidentTriplets( VFace2D* virtualFace, pair<VEdge2D*, rg_Circle2D>& edgeNCircle, VEdge2D*& incidentTriplet1, VEdge2D*& incidentTriplet2 )
{
    VEdge2D*        flippingEdge = edgeNCircle.first;
    rg_Circle2D     circumcircle = edgeNCircle.second;

    flippingEdge->flip();

    VVertex2D*      startVtx = flippingEdge->getStartVertex();
    VVertex2D*      endVtx   = flippingEdge->getEndVertex();
    VFace2D*        matingFaceOfStartVtx = getMatingFace( startVtx, flippingEdge );

    if( matingFaceOfStartVtx == virtualFace ) {
        startVtx->setCircumcircle( circumcircle );
        incidentTriplet1 = flippingEdge->getLeftHand();
        incidentTriplet2 = flippingEdge->getRightHand();

		// compute tangent circle for contact map by Joonghyun March, 17
		updateTangentCircleForContactMap(startVtx);
    }
    else {
        endVtx->setCircumcircle( circumcircle );
        incidentTriplet1 = flippingEdge->getLeftLeg();
        incidentTriplet2 = flippingEdge->getRightLeg();

		// compute tangent circle for contact map by Joonghyun March, 17
		updateTangentCircleForContactMap(endVtx);
    }

	// compute tangent circle for contact map by Joonghyun March, 17
	updateTangentCircleForContactMap(flippingEdge);
}


// Previous version which uses linked list as a priority queue

void VoronoiDiagram2DC::reflectTheFlipOnTheIncidentTriplet_deletion( rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges, VEdge2D*& edgeOfIncidentTriplet, VFace2D*& virtualFace )
{
    rg_dNode< pair<VEdge2D*, rg_Circle2D> >* updateNode = NULL;
    rg_dNode< pair<VEdge2D*, rg_Circle2D> >* currNode   = possiblyFlippingEdges.getFirstpNode();
    for( int i = 0 ; i < possiblyFlippingEdges.getSize() ; i++, currNode = currNode->getNext() ) {
        VEdge2D* currEdge = currNode->getpEntity()->first;

        if( currEdge == edgeOfIncidentTriplet ) {
            updateNode = currNode;
            break;
        }
    }

    rg_Circle2D circumcircle;
    bool isFlippingEdge = isFlippingEdgeOnVirtualRegion( edgeOfIncidentTriplet, circumcircle, virtualFace );

    if( isFlippingEdge ) {
        if( updateNode == NULL ) {
            possiblyFlippingEdges.add( pair<VEdge2D*, rg_Circle2D>( edgeOfIncidentTriplet, circumcircle ) );
        }
        else {
            currNode->getpEntity()->second = circumcircle;
        }
    }
    else {
        if( updateNode != NULL ) {
            possiblyFlippingEdges.kill( updateNode );
        }
    }
}



void VoronoiDiagram2DC::reflectTheFlipOnTheIncidentTriplet( Priority_Q_VEdge & possiblyFlippingEdges, VEdge2D *& edgeOfIncidentTriplet, VFace2D* const virtualFace )
{
    rg_Circle2D circumcircle;
    bool isFlippingEdge = isFlippingEdgeOnVirtualRegion( edgeOfIncidentTriplet, circumcircle, virtualFace );

    if( isFlippingEdge )
    {
        if( !edgeOfIncidentTriplet->isCandidateToBeFlippedInPhantomRemoval() )
        {
            possiblyFlippingEdges.push( pair<VEdge2D*, rg_Circle2D>( edgeOfIncidentTriplet, circumcircle ), circumcircle.getRadius() );
            edgeOfIncidentTriplet->isCandidateToBeFlippedInPhantomRemoval(true);
        }
        else
        {
            possiblyFlippingEdges.changeCircumCircle( edgeOfIncidentTriplet, circumcircle );
        }
    }
    else
    {
        if( edgeOfIncidentTriplet->isCandidateToBeFlippedInPhantomRemoval() )
        {
            possiblyFlippingEdges.kill( edgeOfIncidentTriplet );
            edgeOfIncidentTriplet->isCandidateToBeFlippedInPhantomRemoval(false);
        }
    }
}



bool VoronoiDiagram2DC::isCorrectCircumcircle( rg_Circle2D& circumcircle, VFace2D* prevFace, VFace2D* currFace, VFace2D* nextFace )
{
    rg_Point2D vecCircumcircleToPrevFace = prevFace->getGenerator()->getDisk().getCenterPt() - circumcircle.getCenterPt();
    rg_Point2D vecCircumcircleToCurrFace = currFace->getGenerator()->getDisk().getCenterPt() - circumcircle.getCenterPt();
    rg_Point2D vecCircumcircleToNextFace = nextFace->getGenerator()->getDisk().getCenterPt() - circumcircle.getCenterPt();
    
	double angleNextToPrev = angleFromVec1toVec2(vecCircumcircleToNextFace, vecCircumcircleToPrevFace);
	double angleNextToCurr = angleFromVec1toVec2(vecCircumcircleToNextFace, vecCircumcircleToCurrFace);

	if (angleNextToPrev > angleNextToCurr) {
		return true;
	}
	else {
		return false;
	}
}



pair<VEdge2D*, rg_Circle2D> VoronoiDiagram2DC::dequeueRootNode_DiskTripletWithMinTangentCircle( rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges )
{
    rg_dNode< pair<VEdge2D*, rg_Circle2D> >* nodeWithSmallestCircle = NULL;

    double minRadius = DBL_MAX;
    rg_dNode< pair<VEdge2D*, rg_Circle2D> >* currNode = possiblyFlippingEdges.getFirstpNode();

    for( int i = 0 ; i < possiblyFlippingEdges.getSize() ; i++, currNode = currNode->getNext() ) {
        pair<VEdge2D*, rg_Circle2D>* currEdgeAndCircle = currNode->getpEntity();

        if ( currEdgeAndCircle->first->getStatus() == RED_E ) {
            continue;
        }

        double radius = currEdgeAndCircle->second.getRadius();
        if( radius < minRadius ) {
            minRadius = radius;
            nodeWithSmallestCircle = currNode;
        }
    }

    possiblyFlippingEdges.setHead( nodeWithSmallestCircle );

    return possiblyFlippingEdges.popFront();
}



pair<VEdge2D*, rg_Circle2D> VoronoiDiagram2DC::dequeueRootNode_VEdgeWithCircumcircleFarFromGeneratorToBeRemoved( rg_dList<pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges, const Generator2D* const generatorToBeRemoved )
{
    rg_dNode< pair<VEdge2D*, rg_Circle2D> >* nodeWithCircumcircleFarFromGeneratorToBeRemoved = NULL;

    double maxDistance = -DBL_MAX;
    rg_dNode< pair<VEdge2D*, rg_Circle2D> >* currNode = possiblyFlippingEdges.getFirstpNode();

    for( int i = 0 ; i < possiblyFlippingEdges.getSize() ; i++, currNode = currNode->getNext() ) {
        pair<VEdge2D*, rg_Circle2D>* currEdgeAndCircle = currNode->getpEntity();

        if ( currEdgeAndCircle->first->getStatus() == RED_E ) {
            continue;
        }

        rg_Circle2D circumcircle = currEdgeAndCircle->second;
        double currDistance = circumcircle.getCenterPt().distance(generatorToBeRemoved->getDisk().getCenterPt()) - circumcircle.getRadius();
        if( currDistance > maxDistance ) {
            maxDistance = currDistance;
            nodeWithCircumcircleFarFromGeneratorToBeRemoved = currNode;
        }
    }

    possiblyFlippingEdges.setHead( nodeWithCircumcircleFarFromGeneratorToBeRemoved );

    return possiblyFlippingEdges.popFront();
}



int VoronoiDiagram2DC::getNumberOfVVertices() const
{
    VFace2D* virtualFace = m_VFaces.front();
    list<VVertex2D*> infiniteVertices;
    virtualFace->getBoundaryVVertices( infiniteVertices );

    int     totalNumOfVertices      = m_VVertices.size();
    int     numOfInfiniteVertices   = infiniteVertices.size();
    int     numOfVVerticesWithoutInfiniteVertices = totalNumOfVertices - numOfInfiniteVertices;

    return numOfVVerticesWithoutInfiniteVertices;
}



int VoronoiDiagram2DC::getNumberOfVEdges() const
{
    VFace2D* virtualFace = m_VFaces.front();
    list<VEdge2D*> VirtualEdges;
    virtualFace->getBoundaryVEdges( VirtualEdges );

    int     totalNumOfEdges         = m_VEdges.size();
    int     numOfVirtualEdges       = VirtualEdges.size();
    int     numOfVEdgesWithoutVirtualEdges = totalNumOfEdges - numOfVirtualEdges;

    return numOfVEdgesWithoutVirtualEdges;
}



int VoronoiDiagram2DC::getNumberOfVFaces() const
{
    int totalNumOfFaces               = m_VFaces.size();
    int numOfVFacesWithoutVirtualFace = totalNumOfFaces - 1;

    return numOfVFacesWithoutVirtualFace;
}


void VoronoiDiagram2DC::updateTangentCircleForContactMap(VVertex2D* const vertex)
{
	if (m_directAddressTableFromDiskIDToDiskArrangement == rg_NULL)
		return;

	rg_INT      numberOfCircumcircles = 0;
	rg_Circle2D result1, result2;

	Generator2D* definingGenerator[3];
	//getThreeGeneratorsDefiningVertex(VD_CIC, VVertex, definingGenerator[0], definingGenerator[1], definingGenerator[2]);
	getThreeGeneratorsDefiningVertex(vertex, definingGenerator[0], definingGenerator[1], definingGenerator[2]);

	rg_INT index[3] = { definingGenerator[0]->getID(), definingGenerator[1]->getID(), definingGenerator[2]->getID() };
	for (rg_INT i = 0; i < 3; i++)
	{
		if (index[i] < 0)
			index[i] = index[i] + 1;
	}

#ifdef _DEBUG
	if (index[0] == -1 && index[1] == 2 && index[2] == 5)
		int here = 1;
#endif


	rg_Circle2D circle[3] = { m_directAddressTableFromDiskIDToDiskArrangement[index[0]]->getCircle(),
		                      m_directAddressTableFromDiskIDToDiskArrangement[index[1]]->getCircle(),
		                      m_directAddressTableFromDiskIDToDiskArrangement[index[2]]->getCircle() };

	if ( (VoronoiDiagram2DC*)definingGenerator[0]->getInnerVD() == this ) {
		rg_Circle2D::calculateTangentCircles(circle[0], circle[1], circle[2], result1, result2);
		numberOfCircumcircles = -1;
	}
	else if ( (VoronoiDiagram2DC*)definingGenerator[1]->getInnerVD() == this ) {
		rg_Circle2D::calculateTangentCircles(circle[1], circle[0], circle[2], result1, result2);
		numberOfCircumcircles = -1;
	}
	else if ( (VoronoiDiagram2DC*)definingGenerator[2]->getInnerVD() == this ) {
		rg_Circle2D::calculateTangentCircles(circle[2], circle[1], circle[0], result1, result2);
		numberOfCircumcircles = -1;
	}
	else {
		//numberOfCircumcircles = rg_Circle2D::makeCircumcircle(circle[0], circle[1], circle[2], result1, result2);
		//DiskPackingSolver?ÔøΩÏÑú??tolerance?Ôø? ?ÔøΩÏπò ?ÔøΩÌÇ§Ôø???ÔøΩÌï¥ Ôø?Í≤ΩÌï®.
		rg_REAL tolerance = 0.0;
		numberOfCircumcircles = rg_Circle2D::makeCircumcircle(circle[0], circle[1], circle[2], result1, result2, tolerance);
	}

	if (numberOfCircumcircles == 0) {
		// error
		//printGenerators();
		exit(1);
	}

	rg_Point2D vector1 = circle[0].getCenterPt() - result1.getCenterPt();
	rg_Point2D vector2 = circle[1].getCenterPt() - result1.getCenterPt();
	rg_Point2D vector3 = circle[2].getCenterPt() - result1.getCenterPt();

	if ( (VoronoiDiagram2DC*)definingGenerator[0]->getInnerVD() == this ) {
		vector1 = -1.0 * vector1;
	}
	if ( (VoronoiDiagram2DC*)definingGenerator[1]->getInnerVD() == this ) {
		vector2 = -1.0 * vector2;
	}
	if ( (VoronoiDiagram2DC*)definingGenerator[2]->getInnerVD() == this ) {
		vector3 = -1.0 * vector3;
	}

	//double angle12 = angleFromVec1toVec2(vector1, vector2);
	//double angle13 = angleFromVec1toVec2(vector1, vector3);
	//DiskPackingSolver?ÔøΩÏÑú??tolerance?Ôø? ?ÔøΩÏπò ?ÔøΩÌÇ§Ôø???ÔøΩÌï¥ Ôø?Í≤ΩÌï®.
	rg_REAL tolerance = 0.0;
	double angle12 = rg_GeoFunc::angleFromVec1toVec2(vector1, vector2, tolerance);
	double angle13 = rg_GeoFunc::angleFromVec1toVec2(vector1, vector3, tolerance);

	rg_Circle2D circumCircle;

	switch (numberOfCircumcircles)
	{
	case 1:
		circumCircle = result1;
		break;
	case 2:
		if (angle12 < angle13) {
			circumCircle = result1;
		}
		else {
			circumCircle = result2;
		}
		break;
	default: // numberOfCircumcircles == -1, tangent outside to two circles and inside to the infinite circle.
		if (angle12 < angle13) {
			circumCircle = result1;
		}
		else {
			circumCircle = result2;
		}
		break;
	}

	vertex->setTangentCircle( circumCircle );
}


void VoronoiDiagram2DC::resetAnomalyTestDoneToggle( list<VEdge2D*>& anomalyTestedVEdges )
{
    list<VEdge2D*>::iterator it_testedEdge = anomalyTestedVEdges.begin();
    for( ; it_testedEdge != anomalyTestedVEdges.end(); it_testedEdge++ ) {
        VEdge2D* currTestedEdge = *it_testedEdge;
        currTestedEdge->setFalseAnomalyTestDone();
    }
}

void VoronoiDiagram2DC::setGeneratorsWithoutSorting( list<rg_Circle2D>& circleset )
{
    m_generators.clear();

    int currID = 1;
    for( list<rg_Circle2D>::iterator i_disk = circleset.begin() ; i_disk != circleset.end() ; ++i_disk, ++currID ) {
        rg_Circle2D currDisk = *i_disk;
        float xCoord = currDisk.getCenterPt().getX();
        float yCoord = currDisk.getCenterPt().getY();
        float radius = currDisk.getRadius();

        rg_Circle2D currDisk_single_precision( xCoord, yCoord, radius );
        Generator2D* newGen = createGenerator();
        newGen->setID( currID );
        newGen->setDisk( currDisk_single_precision );
    }
}

//VoronoiDiagram2DC & VoronoiDiagram2DC::operator=(const VoronoiDiagram2DC & VD2DC)
//{
//	if(this != &VD2DC)
//		copyVoronoiDiagramFrom(VD2DC);
//
//	return *this;
//}

Generator2D* VoronoiDiagram2DC::deleteOneGenerator( Generator2D* generator ) //generator -> thisGenerator
{
    VFace2D* VFaceToBeRemoved = generator->getOuterFace();
    list<VEdge2D*> boundingVEdges;
    VFaceToBeRemoved->getBoundaryVEdges(boundingVEdges);

    int numBoundingVEdges = boundingVEdges.size();
    

	// VEdges to compute tangent circle for contact map by Joonghyun March, 17
	list<VEdge2D*> VEdgesToRecomputeTangentCircle;

    rg_dList< pair<VEdge2D*, rg_Circle2D> > possiblyFlippingEdges;

    switch (numBoundingVEdges)
    {
    case 0: // error
    case 1: // error
        exit(1); // LOG COMMENT
        break;


    case 2: // anomaly
		//AfxMessageBox(_T("Anomaly face removal attempted."));
        removeVFaceBoundedByTwoVEdges(VFaceToBeRemoved);

        break;


    case 3:
        removeVFaceBoundedByThreeVEdges(VFaceToBeRemoved);
        
        break;


    default:
        // collect possibly flipped edges
        collectPossiblyFlippingEdgesOnTheInterfaceForBoolean_StoredInPriorityQueue(VFaceToBeRemoved, possiblyFlippingEdges); //possiblyFlippingEdgesQ
        
		bool		isThereVEdgesWithIdenticalDiskTriplet = false;
		VEdge2D*	firstVEdgeWithIdenticalDiskTriplet    = rg_NULL;
		VEdge2D*	secondVEdgeWithIdenticalDiskTriplet   = rg_NULL;

        //pair<VEdge2D*, rg_Circle2D> edgeToDefineMinTangentCircle = dequeueRootNode_DiskTripletWithMinTangentCircle(possiblyFlippingEdges);
        set<Generator2D*> allTrappedGenerators;

        while ( numBoundingVEdges > 3 ) {
            //pair<VEdge2D*, rg_Circle2D> edgeToDefineMinTangentCircle = dequeueRootNode_DiskTripletWithMinTangentCircle(possiblyFlippingEdges);
            pair<VEdge2D*, rg_Circle2D> edgeToDefineMinTangentCircle = dequeueRootNode_VEdgeWithCircumcircleFarFromGeneratorToBeRemoved(possiblyFlippingEdges, generator);
            

            VEdge2D* prevQuillEdge_CCW = rg_NULL;
            VEdge2D* nextQuillEdge_CCW = rg_NULL;
            if ( thereIsVEdgeWhichHasIdenticalDiskTripletBetweenAdjacentVEdgesOnBoundaryOfInputVFace( edgeToDefineMinTangentCircle.first, VFaceToBeRemoved, prevQuillEdge_CCW, nextQuillEdge_CCW ) ) {
                //list<VEdge2D*>  trappingEdges;
                //VFace2D*        trappingFace;
                //findTrappingVEdges( VFaceToBeRemoved, trappingEdges, trappingFace );
                //
                //VEdge2D* trappingEdge1  = trappingEdges.front();
                //VEdge2D* trappingEdge2  = trappingEdges.back();
                //
                //VEdge2D* prevQuillEdge_CCW = trappingEdge1;
                //VEdge2D* nextQuillEdge_CCW = trappingEdge2;
                //
                //if ( (VoronoiDiagram2DC*)trappingFace->getGenerator()->getInnerVD() == this ) {
                //    int numVEdgesPrevQuillEdgeToNextQuillEdge = 0;
                //    int numVEdgesNextQuillEdgeToPrevQuillEdge = 0;
                //
                //    VEdge2D* nextEdge_CCW = getCCWNextEdgeOnFace(prevQuillEdge_CCW, VFaceToBeRemoved);
                //    while ( nextEdge_CCW != nextQuillEdge_CCW ) {
                //        ++numVEdgesPrevQuillEdgeToNextQuillEdge;
                //        nextEdge_CCW = getCCWNextEdgeOnFace(nextEdge_CCW, VFaceToBeRemoved);
                //    }
                //
                //    numVEdgesNextQuillEdgeToPrevQuillEdge = numBoundingVEdges - numVEdgesPrevQuillEdgeToNextQuillEdge - 2;
                //
                //    if ( numVEdgesNextQuillEdgeToPrevQuillEdge < numVEdgesPrevQuillEdgeToNextQuillEdge ) {
                //        prevQuillEdge_CCW = trappingEdge2;
                //        nextQuillEdge_CCW = trappingEdge1;
                //    }
                //}

                set<Generator2D*> currTrappedGenerators;
                VFace2D* trappingFace = getFaceOfOppositeSide(VFaceToBeRemoved, prevQuillEdge_CCW);
                bool thisIsProperOrientationOfTrappingEdges = findTrappedGeneratorsAndEntities(VFaceToBeRemoved, trappingFace, prevQuillEdge_CCW, nextQuillEdge_CCW, currTrappedGenerators);

                if ( !thisIsProperOrientationOfTrappingEdges && trappingFace != m_VFaces.front() ) {
                    currTrappedGenerators.clear();
                    VEdge2D* tempEdge = prevQuillEdge_CCW;
                    prevQuillEdge_CCW = nextQuillEdge_CCW;
                    nextQuillEdge_CCW = tempEdge;
                    findTrappedGeneratorsAndEntities(VFaceToBeRemoved, trappingFace, prevQuillEdge_CCW, nextQuillEdge_CCW, currTrappedGenerators);
                }

                int numTrappedEdgeAmongBoundaryEdge = 0;
                VEdge2D* nextEdge_CCW = getCCWNextEdgeOnFace(prevQuillEdge_CCW, VFaceToBeRemoved);
                while ( nextEdge_CCW != nextQuillEdge_CCW ) {
                    ++numTrappedEdgeAmongBoundaryEdge;
                    nextEdge_CCW = getCCWNextEdgeOnFace(nextEdge_CCW, VFaceToBeRemoved);
                }

                numBoundingVEdges -= numTrappedEdgeAmongBoundaryEdge;

                temporarilyRemoveTrappedGenerators(currTrappedGenerators);

                allTrappedGenerators.insert(currTrappedGenerators.begin(), currTrappedGenerators.end());

                VEdge2D* mergedVEdge = rg_NULL;
                mergedVEdge = mergeTwoTrappingVEdges(prevQuillEdge_CCW, nextQuillEdge_CCW, VFaceToBeRemoved);

                --numBoundingVEdges;

                reflectMergedEdge( mergedVEdge, VFaceToBeRemoved, possiblyFlippingEdges );

				// VEdges to compute tangent circle for contact map by Joonghyun March, 17
				VEdgesToRecomputeTangentCircle.push_back(mergedVEdge);

                continue;
            }
            else {
                VEdge2D* edgeOfIncidentTriplet[2] = { NULL, NULL };
                flipThisEdgeNFindTwoIncidentTriplets(VFaceToBeRemoved, edgeToDefineMinTangentCircle, edgeOfIncidentTriplet[0], edgeOfIncidentTriplet[1]);

                // decrement the counter by one
                --numBoundingVEdges;

                for (int j = 0; j < 2; j++) {
                    if (edgeOfIncidentTriplet[j]->isBounded()) {
                        reflectTheFlipOnTheIncidentTriplet_deletion(possiblyFlippingEdges, edgeOfIncidentTriplet[j], VFaceToBeRemoved);
                    }
                }

				// VEdges to compute tangent circle for contact map by Joonghyun March, 17
				VEdgesToRecomputeTangentCircle.push_back(edgeToDefineMinTangentCircle.first);
            }
        }

        if ( numBoundingVEdges == 2 ) {
            removeVFaceBoundedByTwoVEdges(VFaceToBeRemoved);
        }
        else {
            removeVFaceBoundedByThreeVEdges(VFaceToBeRemoved);
        }        

        // insert temporarily removed generators into VD
        set<Generator2D*>::iterator i_gen;
        for ( i_gen = allTrappedGenerators.begin(); i_gen != allTrappedGenerators.end(); ++i_gen ) {
            Generator2D* currGen = *i_gen;
            updateVoronoiDiagramWithNewGenerator(currGen);
        }

        break;
    }

	// compute tangent circle for contact map by Joonghyun March, 17
	updateTangentCircleForContactMap(VEdgesToRecomputeTangentCircle);

    removeAllExtraneousVVerticesAndVEdges();


    return generator;

    /*



		//isThereVEdgesWithIdenticalDiskTriplet = findTwoVEdgesWhichHaveIdenticalDiskTriplet( VFaceToBeRemoved, firstVEdgeWithIdenticalDiskTriplet, secondVEdgeWithIdenticalDiskTriplet );
		// thereAreTwoVEdgesWithIdenticalDiskQuadrupletIncludingThisGenerator
		// firstVEdgeWithIdenticalDiskTriplet -> firstVEdgeWithIdenticalDiskTripletInCCW

		/*if (isThereVEdgesWithIdenticalDiskTriplet) {
			removeVFaceBoundedByTwoVEdgeWithIdenticalDIskTripletAndTrappedRegion(VFaceToBeRemoved, firstVEdgeWithIdenticalDiskTriplet, secondVEdgeWithIdenticalDiskTriplet);
			
		}
		else {
			while (numBoundingVEdges > 3) {
				//isThereVEdgesWithIdenticalDiskTriplet = findTwoVEdgesWhichHaveIdenticalDiskTriplet(VFaceToBeRemoved, firstVEdgeWithIdenticalDiskTriplet, secondVEdgeWithIdenticalDiskTriplet);
				//if (isThereVEdgesWithIdenticalDiskTriplet ) {
				//	removeVFaceBoundedByTwoVEdgeWithIdenticalDIskTripletAndTrappedRegion(VFaceToBeRemoved, firstVEdgeWithIdenticalDiskTriplet, secondVEdgeWithIdenticalDiskTriplet);
				//	//break;
				//}

				pair<VEdge2D*, rg_Circle2D> edgeToDefineMinTangentCircle = dequeueRootNode_DiskTripletWithMinTangentCircle(possiblyFlippingEdges);
                
                if (thereIsVEdgeWhichHasIdenticalDiskTripletBetweenAdjacentVEdgesOnBoundaryOfInputVFace(edgeToDefineMinTangentCircle.first, VFaceToBeRemoved)) {
                    list<VEdge2D*> trappingEdges;
                    VFace2D*       trappingFace;
                    findTrappingVEdges(VFaceToBeRemoved, trappingEdges, trappingFace);

                    //findTrappedGenerators
                    //removeTrappedGenerators
                    //change the varialble numBoundingVEdges
                    //refletThemergedEdgeNadjacentVEdges
                    continue;
                }


				//if ( isThereVEdgesWithIdenticalDiskTriplet &&
				//	(edgeToDefineMinTangentCircle.first == firstVEdgeWithIdenticalDiskTriplet ||
				//	 edgeToDefineMinTangentCircle.first == secondVEdgeWithIdenticalDiskTriplet) )
				//{
				//	dequeueRootNode_DiskTripletWithMinTangentCircle(possiblyFlippingEdges);
				//	isThereVEdgesWithIdenticalDiskTriplet = false;
				//	
				//	VEdge2D* EdgeToBeFlipped = findEdgeToBeFlippedBetweenTwoEdgesWithIdenticalDiskTriplet(VFaceToBeRemoved, firstVEdgeWithIdenticalDiskTriplet, secondVEdgeWithIdenticalDiskTriplet);
				//	// ?ÔøΩÎñ§ edgeÔø? flip?ÔøΩÔøΩ? ?ÔøΩÎã®?ÔøΩÏ£º??ÏΩîÎìú
				//	edgeToDefineMinTangentCircle.first = EdgeToBeFlipped;
				//}

				VEdge2D* edgeOfIncidentTriplet[2] = { NULL, NULL };
				flipThisEdgeNFindTwoIncidentTriplets(VFaceToBeRemoved, edgeToDefineMinTangentCircle, edgeOfIncidentTriplet[0], edgeOfIncidentTriplet[1]);

				

				// decrement the counter by one
				--numBoundingVEdges;

				/*if (!isThereVEdgesWithIdenticalDiskTriplet) {
					isThereVEdgesWithIdenticalDiskTriplet = findTwoVEdgesWhichHaveIdenticalDiskTriplet(VFaceToBeRemoved, firstVEdgeWithIdenticalDiskTriplet, secondVEdgeWithIdenticalDiskTriplet);
				}

				if (isThereVEdgesWithIdenticalDiskTriplet && numBoundingVEdges > 3) {
					removeVFaceBoundedByTwoVEdgeWithIdenticalDIskTripletAndTrappedRegion(VFaceToBeRemoved, firstVEdgeWithIdenticalDiskTriplet, secondVEdgeWithIdenticalDiskTriplet);
				}

				for (int j = 0; j < 2; j++) {
					if (edgeOfIncidentTriplet[j]->isBounded()) {
						reflectTheFlipOnTheIncidentTriplet(possiblyFlippingEdges, edgeOfIncidentTriplet[j], VFaceToBeRemoved);
					}
				}
			}

			if (!isThereVEdgesWithIdenticalDiskTriplet) {
				removeVFaceBoundedByThreeVEdges(VFaceToBeRemoved);
			}
			
		//}        
    }


    removeAllExtraneousVVerticesAndVEdges();

    return generator;

    */
}


Generator2D* V::GeometryTier::VoronoiDiagram2DC::deleteOneGenerator( Generator2D* generator, list<VEdge2D*>& influencedEdges )
{
    VFace2D* VFaceToBeRemoved = generator->getOuterFace();
    list<VEdge2D*>  boundingVEdges;
    VFaceToBeRemoved->getBoundaryVEdges( boundingVEdges );
    int numBoundingVEdges = boundingVEdges.size();

    list<VEdge2D*> quillEdges;
    VFaceToBeRemoved->getQuillVEdges( quillEdges );
    // for ( list<VEdge2D*>::iterator i_edge = boundingVEdges.begin(); i_edge != boundingVEdges.end(); ++i_edge ) {
    //     VEdge2D* currEdge = *i_edge;
    //     VEdge2D* currQuillEdge = currEdge->getLeftFace() == VFaceToBeRemoved ? currEdge->getRightHand() : currEdge->getLeftLeg();
    //     quillEdges.push_back( currQuillEdge );
    // }

    rg_dList< pair<VEdge2D*, rg_Circle2D> > possiblyFlippingEdges;
    switch ( numBoundingVEdges )
    {
    case 0: // error
    case 1: // error
        exit( 1 ); // LOG COMMENT
        break;

    case 2: // anomaly
        //AfxMessageBox(_T("Anomaly face removal attempted."));
        removeVFaceBoundedByTwoVEdges( VFaceToBeRemoved );
        break;

    case 3:
        removeVFaceBoundedByThreeVEdges( VFaceToBeRemoved );
        break;

    default:
        // collect possibly flipped edges
        collectPossiblyFlippingEdgesOnTheInterfaceForBoolean_StoredInPriorityQueue( VFaceToBeRemoved, possiblyFlippingEdges ); //possiblyFlippingEdgesQ

        while ( numBoundingVEdges > 3 ) {
            pair<VEdge2D*, rg_Circle2D> edgeToDefineMinTangentCircle = dequeueRootNode_VEdgeWithCircumcircleFarFromGeneratorToBeRemoved( possiblyFlippingEdges, generator );

            VEdge2D* edgeOfIncidentTriplet[2] = { NULL, NULL };
            flipThisEdgeNFindTwoIncidentTriplets( VFaceToBeRemoved, edgeToDefineMinTangentCircle, edgeOfIncidentTriplet[0], edgeOfIncidentTriplet[1] );
            // flipThisEdgeNFindTwoIncidentBoundaryEdges
            for ( int j = 0; j < 2; j++ ) {
                if ( edgeOfIncidentTriplet[j]->isBounded() ) {
                    reflectTheFlipOnTheIncidentTriplet_deletion( possiblyFlippingEdges, edgeOfIncidentTriplet[j], VFaceToBeRemoved );
                    // checkFlippabilityOfThisEdge_calcFlipTime_doBookkeepingWithPriorityQ
                    // Ω«¡¶∑Œ¥¬ circumcircle¿Ã input generator∏¶ ¿‚æ∆∏‘¥¬ ¡§µµ∞° key∞° µ«¡ˆ∏∏
                    // input generator¿« π›¡ˆ∏ß¿Ã Ω√∞£ø° linear«œ∞‘ ¡ŸæÓµÁ¥Ÿ∞Ìƒ°∞Ì flipTime¿ª priorityQ¿« key∑Œ æ¥¥Ÿ«œ¿⁄.
                    // distance¥¬ Ω√∞£∞˙ monotonic function¿Ã¥Ÿ.
                    // Bookkeeping: 1) enqueue, 2) dequeue, 3) bubble up/down
                }
            }
            // decrement the counter by one
            --numBoundingVEdges;
        }

        if ( numBoundingVEdges == 2 ) {
            removeVFaceBoundedByTwoVEdges( VFaceToBeRemoved );
        }
        else {
            removeVFaceBoundedByThreeVEdges( VFaceToBeRemoved );
        }
        break;
    }

    for ( list<VEdge2D*>::iterator i_edge = boundingVEdges.begin(); i_edge != boundingVEdges.end(); ++i_edge ) {
        VEdge2D* currBoundingEdge = *i_edge;
        if ( currBoundingEdge->getStatus() == RED_E ) {
            continue;
        }

        influencedEdges.push_back( currBoundingEdge );
    }

    for ( list<VEdge2D*>::iterator i_edge = quillEdges.begin(); i_edge != quillEdges.end(); ++i_edge ) {
        VEdge2D* currQuillEdge = *i_edge;
        if ( currQuillEdge->getStatus() == RED_E ) {
            continue;
        }

        influencedEdges.push_back( currQuillEdge );
    }

    removeAllExtraneousVVerticesAndVEdges();
    return generator;
}


Generator2D* VoronoiDiagram2DC::deleteOneGenerator_deleteFromGeneratorList(Generator2D* generator)
{
    m_generators.remove(generator);
    return deleteOneGenerator(generator);
}


Generator2D* VoronoiDiagram2DC::deleteOneGenerator_deleteFromGeneratorList(Generator2D* generator, list<VEdge2D*>& influencedEdges)
{
    m_generators.remove(generator);
    return deleteOneGenerator(generator, influencedEdges);
}


void VoronoiDiagram2DC::insertOneGenerator(Generator2D * generator)
{
	updateVoronoiDiagramWithNewGenerator(generator);
	removeAllExtraneousVVerticesAndVEdges();
}



void VoronoiDiagram2DC::insertOneGenerator_addToGeneratorList(Generator2D* generator)
{
    m_generators.push_back(generator);
    insertOneGenerator(generator);
}



Generator2D * VoronoiDiagram2DC::insertDiskIntoVD( const rg_Circle2D & disk )
{

    int lastID = m_generators.back()->getID();
    Generator2D* newGen  = createGenerator();
    newGen->setID( lastID+1 );
    newGen->setDisk( disk );
    //m_generators.push_back(Generator2D(lastID+1, disk));

    Generator2D* currGenerator = m_generators.back();
    updateVoronoiDiagramWithNewGenerator( currGenerator );

    removeAllExtraneousVVerticesAndVEdges();

    //updateDiskID2VFaceMap_AfterDiskInsertion( currGenerator );

    return currGenerator;
}



Generator2D * VoronoiDiagram2DC::insertDiskIntoVD(pair<rg_Circle2D, void*>& diskNUserDataPair)
{
    int lastID = m_generators.back()->getID();
    Generator2D* newGen = createGenerator();
    newGen->setID( lastID + 1 );
    newGen->setDisk( diskNUserDataPair.first );
    newGen->setUserData( diskNUserDataPair.second );

    Generator2D* currGenerator = m_generators.back();
    updateVoronoiDiagramWithNewGenerator(currGenerator);

    removeAllExtraneousVVerticesAndVEdges();

    //updateDiskID2VFaceMap_AfterDiskInsertion( currGenerator );

    return currGenerator;
}



Generator2D * VoronoiDiagram2DC::insertDiskIntoVD(const rg_Triplet<rg_Circle2D, int, void*>& diskAndIDAndUserData)
{
    Generator2D* newGen = createGenerator();
    newGen->setID( diskAndIDAndUserData.m_second );
    newGen->setDisk( diskAndIDAndUserData.m_first );
    newGen->setUserData( diskAndIDAndUserData.m_third );

    updateVoronoiDiagramWithNewGenerator(newGen);

    removeAllExtraneousVVerticesAndVEdges();

    //updateDiskID2VFaceMap_AfterDiskInsertion( currGenerator );

    return newGen;
}



Generator2D * VoronoiDiagram2DC::insertDiskIntoVD_new_influenced_removed_edge( const rg_Triplet<rg_Circle2D, int, void*>& diskAndIDAndUserData, list<VEdge2D*>& newEdges, list<VEdge2D*>& influencedEdges, list<VEdge2D*>& removedEdges )
{
    Generator2D* newGen = createGenerator();
    newGen->setID( diskAndIDAndUserData.m_second );
    newGen->setDisk( diskAndIDAndUserData.m_first );
    newGen->setUserData( diskAndIDAndUserData.m_third );

    //updateVoronoiDiagramWithNewGenerator(currGenerator);
    updateVoronoiDiagramWithNewGenerator_new_influenced_removed_edges(newGen, newEdges, influencedEdges, removedEdges);

    //removeAllExtraneousVVerticesAndVEdges();

    return newGen;
}



Generator2D * VoronoiDiagram2DC::insertDiskIntoVD_for_initial_arrangement_of_diskpacking( const rg_Triplet<rg_Circle2D, int, void*>& diskAndIDAndUserData, list<VVertex2D*>& redVertices, list<VEdge2D*>& redEdges, list<VVertex2D*>& newVertices, list<VEdge2D*>& newEdges, list<VEdge2D*>& influencedEdges )
{
    Generator2D* newGen = createGenerator();
    newGen->setID( diskAndIDAndUserData.m_second );
    newGen->setDisk( diskAndIDAndUserData.m_first );
    newGen->setUserData( diskAndIDAndUserData.m_third );
    
    updateVoronoiDiagramWithNewGenerator_for_initial_arrangement_of_diskpacking(newGen, redVertices, redEdges, newVertices, newEdges, influencedEdges);

    return newGen;
}



void VoronoiDiagram2DC::collectPossiblyFlippingEdgesOnTheInterfaceForBoolean_StoredInPriorityQueue( VFace2D* toBeRemovedFace, rg_dList< pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges )
{
    list<VEdge2D*> edgesToBoundVirtualFace;
    toBeRemovedFace->getBoundaryVEdges( edgesToBoundVirtualFace );

    list<VEdge2D*>::iterator it_boundingEdge = edgesToBoundVirtualFace.begin();
    for( ; it_boundingEdge != edgesToBoundVirtualFace.end() ; ++it_boundingEdge ) {
        VEdge2D* currEdge = *it_boundingEdge;
        // Unbounded or virtual edges can not be flipped.
        if( !currEdge->isBounded() ) {
            continue;
        }

        rg_Circle2D circumcircle;
        bool isFlippingEdge = isFlippingEdgeOnVirtualRegion( currEdge, circumcircle, toBeRemovedFace );

        if( isFlippingEdge ) {
            possiblyFlippingEdges.add( pair<VEdge2D*, rg_Circle2D>(currEdge, circumcircle) );
        }
    }
}



bool VoronoiDiagram2DC::findTwoVEdgesWhichHaveIdenticalDiskTriplet(VFace2D* const faceToBeRemoved, VEdge2D *& edge1, VEdge2D *& edge2)
{
	bool isThereVEdgesWhichHaveIdenticalDiskTriplet = false;

	VFace2D* prevOppositeFace = rg_NULL;
	VFace2D* currOppositeFace = rg_NULL;
	VFace2D* nextOppositeFace = rg_NULL;

	VEdge2D* firstEdgeOfFaceToBeRemoved = faceToBeRemoved->getFirstVEdge();


	// Find first triplet
	VEdge2D* currBoundaryEdge = firstEdgeOfFaceToBeRemoved;
	VEdge2D* prevBoundaryEdge = rg_NULL;
	VEdge2D* nextBoundaryEdge = rg_NULL;

	if (currBoundaryEdge->getLeftFace() == faceToBeRemoved) {
		prevBoundaryEdge = currBoundaryEdge->getLeftLeg();
		nextBoundaryEdge = currBoundaryEdge->getLeftHand();

		currOppositeFace = currBoundaryEdge->getRightFace();
	}
	else {
		prevBoundaryEdge = currBoundaryEdge->getRightHand();
		nextBoundaryEdge = currBoundaryEdge->getRightLeg();

		currOppositeFace = currBoundaryEdge->getLeftFace();
	}

	if (prevBoundaryEdge->getLeftFace() == faceToBeRemoved) {
		prevOppositeFace = prevBoundaryEdge->getRightFace();
	}
	else {
		prevOppositeFace = prevBoundaryEdge->getLeftFace();
	}

	if (nextBoundaryEdge->getLeftFace() == faceToBeRemoved) {
		nextOppositeFace = nextBoundaryEdge->getRightFace();
	}
	else {
		nextOppositeFace = nextBoundaryEdge->getLeftFace();
	}




	// Find two edges with identical disk triplet
	do {
		prevBoundaryEdge = currBoundaryEdge;
		currBoundaryEdge = nextBoundaryEdge;
		
		if (nextBoundaryEdge->getLeftFace() == faceToBeRemoved) {
			nextBoundaryEdge = nextBoundaryEdge->getLeftHand();
		}
		else {
			nextBoundaryEdge = nextBoundaryEdge->getRightLeg();
		}



		VFace2D* faceToBeLeftOutFromTriplet = prevOppositeFace;
		VFace2D* faceToBeCameInTheTriplet;

		if (nextBoundaryEdge->getLeftFace() == faceToBeRemoved) {
			faceToBeCameInTheTriplet = nextBoundaryEdge->getRightFace();
		}
		else {
			faceToBeCameInTheTriplet = nextBoundaryEdge->getLeftFace();
		}


		if (faceToBeLeftOutFromTriplet == faceToBeCameInTheTriplet) {
			isThereVEdgesWhichHaveIdenticalDiskTriplet = true;
			edge1 = prevBoundaryEdge;
			edge2 = currBoundaryEdge;
			break;
		}

	} while (currBoundaryEdge == firstEdgeOfFaceToBeRemoved);

	return isThereVEdgesWhichHaveIdenticalDiskTriplet;
}



void VoronoiDiagram2DC::removeVFaceBoundedByThreeVEdges( VFace2D*& toBeRemovedVFace )
{
    VEdge2D*    toBeRemovedEdge1    = rg_NULL;
    VEdge2D*    toBeRemovedEdge2    = rg_NULL;
    VEdge2D*    toBeRemovedEdge3    = rg_NULL;

    VEdge2D*    starEdge1           = rg_NULL;
    VEdge2D*    starEdge2           = rg_NULL;
    VEdge2D*    starEdge3           = rg_NULL;

    VFace2D*    oppositeFace1       = rg_NULL;
    VFace2D*    oppositeFace2       = rg_NULL;
    VFace2D*    oppositeFace3       = rg_NULL;


    // find toBeRemovedEdges and starEdges of toBeRemovedVFace

    toBeRemovedEdge1 = toBeRemovedVFace->getFirstVEdge();
    if( toBeRemovedEdge1->getLeftFace() == toBeRemovedVFace ) {
        toBeRemovedEdge2    = toBeRemovedEdge1->getLeftHand();
        starEdge1           = toBeRemovedEdge1->getRightHand();
        oppositeFace1       = toBeRemovedEdge1->getRightFace();
        toBeRemovedEdge1->getEndVertex()->setStatus(RED_V);
    }
    else {
        toBeRemovedEdge2    = toBeRemovedEdge1->getRightLeg();
        starEdge1           = toBeRemovedEdge1->getLeftLeg();
        oppositeFace1       = toBeRemovedEdge1->getLeftFace();
        toBeRemovedEdge1->getStartVertex()->setStatus(RED_V);
    }

    if( toBeRemovedEdge2->getLeftFace() == toBeRemovedVFace ) {
        toBeRemovedEdge3    = toBeRemovedEdge2->getLeftHand();
        starEdge2           = toBeRemovedEdge2->getRightHand();
        oppositeFace2       = toBeRemovedEdge2->getRightFace();
        toBeRemovedEdge2->getEndVertex()->setStatus(RED_V);
    }
    else {
        toBeRemovedEdge3    = toBeRemovedEdge2->getRightLeg();
        starEdge2           = toBeRemovedEdge2->getLeftLeg();
        oppositeFace2       = toBeRemovedEdge2->getLeftFace();
        toBeRemovedEdge2->getStartVertex()->setStatus(RED_V);
    }

    if( toBeRemovedEdge3->getLeftFace() == toBeRemovedVFace ) {
        starEdge3           = toBeRemovedEdge3->getRightHand();
        oppositeFace3       = toBeRemovedEdge3->getRightFace();
        toBeRemovedEdge3->getEndVertex()->setStatus(RED_V);
    }
    else {
        starEdge3           = toBeRemovedEdge3->getLeftLeg();
        oppositeFace3       = toBeRemovedEdge3->getLeftFace();
        toBeRemovedEdge3->getStartVertex()->setStatus(RED_V);
    }

    

    // bookkeeping
    VVertex2D* newVVertex = createVertex( m_VVertices.back()->getID()+1 );
    newVVertex->setFirstVEdge(starEdge1);

    oppositeFace1->setFirstVEdge(starEdge1);
    oppositeFace2->setFirstVEdge(starEdge2);
    oppositeFace3->setFirstVEdge(starEdge3);

    if( starEdge1->getLeftHand() == toBeRemovedEdge1 ) {
        starEdge1->setLeftHand(starEdge3);
        starEdge1->setRightHand(starEdge2);
        starEdge1->setEndVertex(newVVertex);
    }
    else {
        starEdge1->setRightLeg(starEdge3);
        starEdge1->setLeftLeg(starEdge2);
        starEdge1->setStartVertex(newVVertex);
    }

    if( starEdge2->getLeftHand() == toBeRemovedEdge2 ) {
        starEdge2->setLeftHand(starEdge1);
        starEdge2->setRightHand(starEdge3);
        starEdge2->setEndVertex(newVVertex);
    }
    else {
        starEdge2->setRightLeg(starEdge1);
        starEdge2->setLeftLeg(starEdge3);
        starEdge2->setStartVertex(newVVertex);
    }

    if( starEdge3->getLeftHand() == toBeRemovedEdge3 ) {
        starEdge3->setLeftHand(starEdge2);
        starEdge3->setRightHand(starEdge1);
        starEdge3->setEndVertex(newVVertex);
    }
    else {
        starEdge3->setRightLeg(starEdge2);
        starEdge3->setLeftLeg(starEdge1);
        starEdge3->setStartVertex(newVVertex);
    }


    // RED MARK
    toBeRemovedEdge1->setStatus(RED_E);
    toBeRemovedEdge2->setStatus(RED_E);
    toBeRemovedEdge3->setStatus(RED_E);

    // REMOVE FACE
    //m_generators.remove(*toBeRemovedVFace->getGenerator());
    m_VFaces.remove(toBeRemovedVFace);
    delete toBeRemovedVFace;

    // COORDINATE OF THE NEW VVERTEX
    updateCircumcircle(newVVertex);

	// compute tangent circle for contact map by Joonghyun March, 17
	// Because the defining generators of VEdge do not change, we need not update tangent circle for contact map

	//updateTangentCircleForContactMap(starEdge1);
	//updateTangentCircleForContactMap(starEdge2);
	//updateTangentCircleForContactMap(starEdge3);
}

void VoronoiDiagram2DC::updateCircumcircle( VVertex2D* const vertex )
{
    // non-infinite vertex
    if( !vertex->isInfinite() ) {
        rg_INT      numberOfCircumcircles = 0;
        rg_Circle2D result1, result2;

        Generator2D* definingGenerator[3];
        getThreeGeneratorsDefiningVertex( vertex, definingGenerator[0], definingGenerator[1], definingGenerator[2] );

        rg_Circle2D circle[3] = {definingGenerator[0]->getDisk(), definingGenerator[1]->getDisk(), definingGenerator[2]->getDisk() };

        numberOfCircumcircles = rg_Circle2D::makeCircumcircle( circle[0], circle[1], circle[2], result1, result2 );

        if ( numberOfCircumcircles == 0 ) {
            // error
            exit(1);
        }

        rg_Point2D vector1 = circle[0].getCenterPt() - result1.getCenterPt();
        rg_Point2D vector2 = circle[1].getCenterPt() - result1.getCenterPt();
        rg_Point2D vector3 = circle[2].getCenterPt() - result1.getCenterPt();

        double angle12 = angleFromVec1toVec2(vector1, vector2);
        double angle13 = angleFromVec1toVec2(vector1, vector3);

        switch ( numberOfCircumcircles )
        {
        case 1:
            // Although there is only one circumcircle, we should ignore the result with wrong orientation.
            //if( angle12 < angle13 ) {
                vertex->setCircumcircle( result1 );
                vertex->setStatus(WHITE_V);
            //}
            /*else {
                            AfxMessageBox( "Vertex Coordinate: a circumcircle has wrong orientation.");
                exit(1);
            }*/

            break;



        case 2:
            if ( angle12 < angle13 ) {
                vertex->setCircumcircle( result1 );
                vertex->setStatus(WHITE_V);
            }
            else {
                vertex->setCircumcircle( result2 );
                vertex->setStatus(WHITE_V);
            }

            break;



        default:
            exit(1); // error

            break;
        }
    }


    // infinite vertex
    else {
        VEdge2D* unboundedEdge = rg_NULL;

        list<VEdge2D*> incidentEdges;
        vertex->getIncident3VEdges(incidentEdges);

        list<VEdge2D*>::iterator i_edge = incidentEdges.begin();
        for( ; i_edge != incidentEdges.end() ; ++i_edge ) {
            VEdge2D* currEdge = *i_edge;

            if( currEdge->isUnBounded() ) {
                unboundedEdge = currEdge;
                break;
            }
        }
        
        rg_Point2D center = m_circleEnclosingVertices.getCenterPt();
        double     radius = m_circleEnclosingVertices.getRadius();

        VVertex2D* vtxOnInfinite = NULL;
        VVertex2D* vtxOnFinite   = NULL;
        rg_Point2D ptOnFinite;

        // check which vertex is on infinite and finite
        if( unboundedEdge->getStartVertex()->isInfinite() ) {
            vtxOnInfinite = unboundedEdge->getStartVertex();
            vtxOnFinite   = unboundedEdge->getEndVertex();
        }
        else {
            vtxOnInfinite = unboundedEdge->getEndVertex();
            vtxOnFinite   = unboundedEdge->getStartVertex();
        }

        ptOnFinite = vtxOnFinite->getLocation();

        rg_Circle2D leftCircle  = unboundedEdge->getLeftFace()->getGenerator()->getDisk();
        rg_Circle2D rightCircle = unboundedEdge->getRightFace()->getGenerator()->getDisk();

        double givenRadius = radius - center.distance( ptOnFinite );

        rg_Circle2D tangentCircle[2];
        int numTangentCircle = computeCircleTangentTo2CirclesGivenRadius( leftCircle, rightCircle, givenRadius, tangentCircle );


        rg_Point2D vectorToRight = rightCircle.getCenterPt() - tangentCircle[0].getCenterPt();
        rg_Point2D vectorToLeft  = leftCircle.getCenterPt() - tangentCircle[0].getCenterPt();

        if( vtxOnInfinite == unboundedEdge->getEndVertex() ) {
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



void VoronoiDiagram2DC::getThreeGeneratorsDefiningVertex( VVertex2D* tVertex, Generator2D*& generator1, Generator2D*& generator2, Generator2D*& generator3 )
{
    VEdge2D* firstEdge = tVertex->getFirstVEdge();
    VEdge2D* tempEdge  = rg_NULL;

    if ( firstEdge->getStartVertex() == tVertex ) {
        tempEdge = firstEdge->getRightLeg();

        generator1 = firstEdge->getRightFace()->getGenerator();
        generator2 = firstEdge->getLeftFace()->getGenerator();
        VFace2D* faceOfGenerator3 = getFaceOfOppositeSide( firstEdge->getRightFace(), tempEdge );
        generator3 = faceOfGenerator3->getGenerator();
    }
    else {
        tempEdge = firstEdge->getLeftHand();

        generator1 = firstEdge->getLeftFace()->getGenerator();
        generator2 = firstEdge->getRightFace()->getGenerator();
        VFace2D* faceOfGenerator3 = getFaceOfOppositeSide( firstEdge->getLeftFace(), tempEdge );
        generator3 = faceOfGenerator3->getGenerator();
    }
}



void VoronoiDiagram2DC::removeVFaceBoundedByTwoVEdges( VFace2D*& toBeRemovedVFace )
{
    VEdge2D* firstVEdge     = rg_NULL;
    VEdge2D* frontQuillEdge = rg_NULL;
    VEdge2D* rearQuillEdge  = rg_NULL;

    bool isFrontQuillEdgeSameDirectionWithFirstVEdge = false;
    bool isRearQuillEdgeSameDirectionWithFirstVEdge = false;

    VEdge2D* mergedEdge = rg_NULL;

    list<VEdge2D*>::iterator i_edge;

    // get quill edges of the VFace
    // headEdge --> frontQuillEdge
    // tailEdge --> rearQuillEdge
    firstVEdge = toBeRemovedVFace->getFirstVEdge();

    if( firstVEdge->getLeftFace() == toBeRemovedVFace ) {
        rearQuillEdge = firstVEdge->getRightLeg();
        frontQuillEdge = firstVEdge->getRightHand();
    }
    else {
        rearQuillEdge = firstVEdge->getLeftLeg();
        frontQuillEdge = firstVEdge->getLeftHand();
    }

    // merge two star edges
    if( firstVEdge->getStartVertex() == rearQuillEdge->getEndVertex() ) {
        isRearQuillEdgeSameDirectionWithFirstVEdge = true;
    }
    if( firstVEdge->getEndVertex() == frontQuillEdge->getStartVertex() ) {
        isFrontQuillEdgeSameDirectionWithFirstVEdge = true;
    }

    mergedEdge = createEdge( m_VEdges.back()->getID() + 1 );
    if( isRearQuillEdgeSameDirectionWithFirstVEdge ) {
        mergedEdge->setLeftFace(    rearQuillEdge->getLeftFace());
        mergedEdge->setRightFace(   rearQuillEdge->getRightFace());
        mergedEdge->setStartVertex( rearQuillEdge->getStartVertex());
        mergedEdge->setRightLeg(    rearQuillEdge->getRightLeg());
        mergedEdge->setLeftLeg(     rearQuillEdge->getLeftLeg());
    }
    else {
        mergedEdge->setLeftFace(    rearQuillEdge->getRightFace());
        mergedEdge->setRightFace(   rearQuillEdge->getLeftFace());
        mergedEdge->setStartVertex( rearQuillEdge->getEndVertex());
        mergedEdge->setRightLeg(    rearQuillEdge->getLeftHand());
        mergedEdge->setLeftLeg(     rearQuillEdge->getRightHand());
    }

    if( isFrontQuillEdgeSameDirectionWithFirstVEdge ) {
        mergedEdge->setEndVertex(   frontQuillEdge->getEndVertex());
        mergedEdge->setRightHand(   frontQuillEdge->getRightHand());
        mergedEdge->setLeftHand(    frontQuillEdge->getLeftHand());
    }
    else {
        mergedEdge->setEndVertex(   frontQuillEdge->getStartVertex());
        mergedEdge->setRightHand(   frontQuillEdge->getLeftLeg());
        mergedEdge->setLeftHand(    frontQuillEdge->getRightLeg());
    }

    // bookkeeping
    mergedEdge->getLeftFace()->setFirstVEdge(   mergedEdge);
    mergedEdge->getRightFace()->setFirstVEdge(  mergedEdge);
    mergedEdge->getStartVertex()->setFirstVEdge(mergedEdge);
    mergedEdge->getEndVertex()->setFirstVEdge(  mergedEdge);


    if(mergedEdge->getLeftLeg()->getLeftHand() == rearQuillEdge) {
        mergedEdge->getLeftLeg()->setLeftHand(  mergedEdge);
    }
    else {
        mergedEdge->getLeftLeg()->setRightLeg(  mergedEdge);
    }


    if(mergedEdge->getRightLeg()->getRightHand() == rearQuillEdge) {
        mergedEdge->getRightLeg()->setRightHand(mergedEdge);
    }
    else {
        mergedEdge->getRightLeg()->setLeftLeg(  mergedEdge);
    }


    if(mergedEdge->getRightHand()->getRightLeg() == frontQuillEdge) {
        mergedEdge->getRightHand()->setRightLeg(mergedEdge);
    }
    else {
        mergedEdge->getRightHand()->setLeftHand(mergedEdge);
    }


    if(mergedEdge->getLeftHand()->getLeftLeg() == frontQuillEdge) {
        mergedEdge->getLeftHand()->setLeftLeg(  mergedEdge);
    }
    else {
        mergedEdge->getLeftHand()->setRightHand(mergedEdge);
    }


    // delete the VFace and generator
    rearQuillEdge->setStatus(RED_E);
    frontQuillEdge->setStatus(RED_E);
    firstVEdge->getStartVertex()->setStatus(RED_V);
    firstVEdge->getEndVertex()->setStatus(RED_V);

    list<VEdge2D*> boundingVEdges;
    toBeRemovedVFace->getBoundaryVEdges(boundingVEdges);

    i_edge = boundingVEdges.begin();
    for( ; i_edge != boundingVEdges.end() ; ++i_edge ) {
        VEdge2D* currEdge = *i_edge;
        currEdge->setStatus(RED_E);
    }

    //m_generators.remove(*toBeRemovedVFace->getGenerator());
    m_VFaces.remove(toBeRemovedVFace);
    delete toBeRemovedVFace;

	// compute tangent circle for contact map by Joonghyun March, 17
	updateTangentCircleForContactMap(mergedEdge);
}



void VoronoiDiagram2DC::removeVFaceBoundedByTwoVEdges(VFace2D *& toBeRemovedVFace, VEdge2D *& mergedVEdge)
{
    VEdge2D* firstVEdge = rg_NULL;
    VEdge2D* frontQuillEdge = rg_NULL;
    VEdge2D* rearQuillEdge = rg_NULL;

    bool isFrontQuillEdgeSameDirectionWithFirstVEdge = false;
    bool isRearQuillEdgeSameDirectionWithFirstVEdge = false;

    VEdge2D* mergedEdge = rg_NULL;

    list<VEdge2D*>::iterator i_edge;

    // get quill edges of the VFace
    // headEdge --> frontQuillEdge
    // tailEdge --> rearQuillEdge
    firstVEdge = toBeRemovedVFace->getFirstVEdge();

    if (firstVEdge->getLeftFace() == toBeRemovedVFace) {
        rearQuillEdge = firstVEdge->getRightLeg();
        frontQuillEdge = firstVEdge->getRightHand();
    }
    else {
        rearQuillEdge = firstVEdge->getLeftLeg();
        frontQuillEdge = firstVEdge->getLeftHand();
    }

    // merge two star edges
    if (firstVEdge->getStartVertex() == rearQuillEdge->getEndVertex()) {
        isRearQuillEdgeSameDirectionWithFirstVEdge = true;
    }
    if (firstVEdge->getEndVertex() == frontQuillEdge->getStartVertex()) {
        isFrontQuillEdgeSameDirectionWithFirstVEdge = true;
    }

    mergedEdge = createEdge(m_VEdges.back()->getID() + 1);
    if (isRearQuillEdgeSameDirectionWithFirstVEdge) {
        mergedEdge->setLeftFace(rearQuillEdge->getLeftFace());
        mergedEdge->setRightFace(rearQuillEdge->getRightFace());
        mergedEdge->setStartVertex(rearQuillEdge->getStartVertex());
        mergedEdge->setRightLeg(rearQuillEdge->getRightLeg());
        mergedEdge->setLeftLeg(rearQuillEdge->getLeftLeg());
    }
    else {
        mergedEdge->setLeftFace(rearQuillEdge->getRightFace());
        mergedEdge->setRightFace(rearQuillEdge->getLeftFace());
        mergedEdge->setStartVertex(rearQuillEdge->getEndVertex());
        mergedEdge->setRightLeg(rearQuillEdge->getLeftHand());
        mergedEdge->setLeftLeg(rearQuillEdge->getRightHand());
    }

    if (isFrontQuillEdgeSameDirectionWithFirstVEdge) {
        mergedEdge->setEndVertex(frontQuillEdge->getEndVertex());
        mergedEdge->setRightHand(frontQuillEdge->getRightHand());
        mergedEdge->setLeftHand(frontQuillEdge->getLeftHand());
    }
    else {
        mergedEdge->setEndVertex(frontQuillEdge->getStartVertex());
        mergedEdge->setRightHand(frontQuillEdge->getLeftLeg());
        mergedEdge->setLeftHand(frontQuillEdge->getRightLeg());
    }

    // bookkeeping
    mergedEdge->getLeftFace()->setFirstVEdge(mergedEdge);
    mergedEdge->getRightFace()->setFirstVEdge(mergedEdge);
    mergedEdge->getStartVertex()->setFirstVEdge(mergedEdge);
    mergedEdge->getEndVertex()->setFirstVEdge(mergedEdge);


    if (mergedEdge->getLeftLeg()->getLeftHand() == rearQuillEdge) {
        mergedEdge->getLeftLeg()->setLeftHand(mergedEdge);
    }
    else {
        mergedEdge->getLeftLeg()->setRightLeg(mergedEdge);
    }


    if (mergedEdge->getRightLeg()->getRightHand() == rearQuillEdge) {
        mergedEdge->getRightLeg()->setRightHand(mergedEdge);
    }
    else {
        mergedEdge->getRightLeg()->setLeftLeg(mergedEdge);
    }


    if (mergedEdge->getRightHand()->getRightLeg() == frontQuillEdge) {
        mergedEdge->getRightHand()->setRightLeg(mergedEdge);
    }
    else {
        mergedEdge->getRightHand()->setLeftHand(mergedEdge);
    }


    if (mergedEdge->getLeftHand()->getLeftLeg() == frontQuillEdge) {
        mergedEdge->getLeftHand()->setLeftLeg(mergedEdge);
    }
    else {
        mergedEdge->getLeftHand()->setRightHand(mergedEdge);
    }


    // delete the VFace and generator
    rearQuillEdge->setStatus(RED_E);
    frontQuillEdge->setStatus(RED_E);
    firstVEdge->getStartVertex()->setStatus(RED_V);
    firstVEdge->getEndVertex()->setStatus(RED_V);

    list<VEdge2D*> boundingVEdges;
    toBeRemovedVFace->getBoundaryVEdges(boundingVEdges);

    i_edge = boundingVEdges.begin();
    for (; i_edge != boundingVEdges.end(); ++i_edge) {
        VEdge2D* currEdge = *i_edge;
        currEdge->setStatus(RED_E);
    }

    //m_generators.remove(*toBeRemovedVFace->getGenerator());
    m_VFaces.remove(toBeRemovedVFace);
    delete toBeRemovedVFace;

    // compute tangent circle for contact map by Joonghyun March, 17
    updateTangentCircleForContactMap(mergedEdge);

    mergedVEdge = mergedEdge;
}



void VoronoiDiagram2DC::removeVFaceBoundedByTwoVEdgeWithIdenticalDIskTripletAndTrappedRegion(VFace2D*& toBeRemovedVFace, VEdge2D* edge1, VEdge2D* edge2)
{
	set<Generator2D*>  trappedGenerators;
	list<VEdge2D*>		edgesOnTrappedRegion;
	list<VVertex2D*>	verticesOnTrappedRegion;

	findTrappedGeneratorsNMarkTheirBoundaryEntities(toBeRemovedVFace, edge1, edge2, trappedGenerators, edgesOnTrappedRegion, verticesOnTrappedRegion);
	// find and mark the entitiees

	removeVFacesDefinedByTemporarilyDeletedGenerators(trappedGenerators);
	// removeVFacesDefinedByTemporarilyDeletedGenerators(trappedGenerators);

	mergeTwoQuillVEdgesOfTrappedRegion(toBeRemovedVFace, edge1, edge2);


	removeVFaceBoundedByThreeVEdges(toBeRemovedVFace);

	
	updateVoronoiDiagramWithTemporarilyDeletedGenerators(trappedGenerators);
}



void VoronoiDiagram2DC::findTrappedGeneratorsNMarkTheirBoundaryEntities(VFace2D* toBeRemovedVFace, VEdge2D* edge1WithIdenticalDiskTriplet, VEdge2D* edge2WithIdenticalDiskTriplet, set<Generator2D*>& trappedGenerators, list<VEdge2D*>& edgesOnTrappedRegion, list<VVertex2D*>& verticesOnTrappedRegion)
{
	VEdge2D* edge1_CCW;
	VEdge2D* edge2_CCW;

	if (edge1WithIdenticalDiskTriplet->getLeftFace() == toBeRemovedVFace)
	{
		if (edge1WithIdenticalDiskTriplet->getLeftHand() == edge2WithIdenticalDiskTriplet) {
			edge1_CCW = edge1WithIdenticalDiskTriplet;
			edge2_CCW = edge2WithIdenticalDiskTriplet;
		}
		else {
			edge1_CCW = edge2WithIdenticalDiskTriplet;
			edge2_CCW = edge1WithIdenticalDiskTriplet;
		}
	}
	else {
		if (edge1WithIdenticalDiskTriplet->getRightLeg() == edge2WithIdenticalDiskTriplet) {
			edge1_CCW = edge1WithIdenticalDiskTriplet;
			edge2_CCW = edge2WithIdenticalDiskTriplet;
		}
		else {
			edge1_CCW = edge2WithIdenticalDiskTriplet;
			edge2_CCW = edge1WithIdenticalDiskTriplet;
		}
	}

	VEdge2D* prevEdge_CCW;
	VEdge2D* nextEdge_CCW;

	if (edge1_CCW->getLeftFace() == toBeRemovedVFace) {
		prevEdge_CCW = edge1_CCW->getLeftLeg();
	}
	else {
		prevEdge_CCW = edge1_CCW->getRightHand();
	}

	if (edge2_CCW->getLeftFace() == toBeRemovedVFace) {
		nextEdge_CCW = edge2_CCW->getLeftHand();
	}
	else {
		nextEdge_CCW = edge2_CCW->getRightLeg();
	}

	VFace2D* trappingFace1 = prevEdge_CCW->getLeftFace();
	VFace2D* trappingFace2 = prevEdge_CCW->getRightFace();

	set<VEdge2D*> edgesOnTrappedRegion2;
	set<VVertex2D*> verticesOnTrappedRegion2;
	set<Generator2D*> trappedGenerators2;

	if (prevEdge_CCW->getLeftFace() == toBeRemovedVFace) {
		prevEdge_CCW->getStartVertex()->setStatus(RED_V);
		verticesOnTrappedRegion2.insert(prevEdge_CCW->getStartVertex());
	}
	else
	{
		prevEdge_CCW->getEndVertex()->setStatus(RED_V);
		verticesOnTrappedRegion2.insert(prevEdge_CCW->getEndVertex());
	}

	if (nextEdge_CCW->getLeftFace() == toBeRemovedVFace) {
		nextEdge_CCW->getEndVertex()->setStatus(RED_V);
	}
	else
	{
		nextEdge_CCW->getStartVertex()->setStatus(RED_V);
	}

	set<VVertex2D*>::iterator i_vtx;
	for (i_vtx = verticesOnTrappedRegion2.begin(); i_vtx != verticesOnTrappedRegion2.end(); ++i_vtx) {
		VVertex2D* currTrappedVertex = *i_vtx;

		list<VEdge2D*> incidentEdge;
		currTrappedVertex->getIncident3VEdges(incidentEdge);

		list<VEdge2D*>::iterator i_incidentEdge;
		for (i_incidentEdge = incidentEdge.begin(); i_incidentEdge != incidentEdge.end(); ++i_incidentEdge) {
			VEdge2D* currIncidentEdge = *i_incidentEdge;

			if (currIncidentEdge == prevEdge_CCW || currIncidentEdge == nextEdge_CCW) {
				continue;
			}
			if (currIncidentEdge->getStatus() == RED_E) {
				continue;
			}

			VVertex2D* currOppositeVertex = currIncidentEdge->getOppositVVertex(currTrappedVertex);
			if (currOppositeVertex->getStatus() != RED_V) {
				currOppositeVertex->setStatus(RED_V);
				verticesOnTrappedRegion2.insert(currOppositeVertex);
			}

			VFace2D* currLeftFace  = currIncidentEdge->getLeftFace();
			VFace2D* currRightFace = currIncidentEdge->getRightFace();

			if (currLeftFace != trappingFace1 && currLeftFace != trappingFace2) {
				trappedGenerators.insert(currLeftFace->getGenerator());
			}

			if (currRightFace != trappingFace1 && currRightFace != trappingFace2) {
				trappedGenerators.insert(currRightFace->getGenerator());
			}

			currIncidentEdge->setStatus(RED_E);
		}
	}
}



void VoronoiDiagram2DC::removeVFacesDefinedByTemporarilyDeletedGenerators(set<Generator2D*>& trappedGenerators)
{
	set<Generator2D*>::iterator i_gen;
	for (i_gen = trappedGenerators.begin(); i_gen != trappedGenerators.end(); ++i_gen) {
		VFace2D* currFace = (*i_gen)->getOuterFace();
		m_VFaces.remove(currFace);
        delete currFace;
	}
}



void VoronoiDiagram2DC::mergeTwoQuillVEdgesOfTrappedRegion(VFace2D*& toBeRemovedVFace, VEdge2D* edge1, VEdge2D* edge2)
{
	VEdge2D* edge1_CCW;
	VEdge2D* edge2_CCW;

	if (edge1->getLeftFace() == toBeRemovedVFace)
	{
		if (edge1->getLeftHand() == edge2) {
			edge1_CCW = edge1;
			edge2_CCW = edge2;
		}
		else {
			edge1_CCW = edge2;
			edge2_CCW = edge1;
		}
	}
	else {
		if (edge1->getRightLeg() == edge2) {
			edge1_CCW = edge1;
			edge2_CCW = edge2;
		}
		else {
			edge1_CCW = edge2;
			edge2_CCW = edge1;
		}
	}

	VEdge2D* prevEdge_CCW;
	VEdge2D* nextEdge_CCW;

	if (edge1_CCW->getLeftFace() == toBeRemovedVFace) {
		prevEdge_CCW = edge1_CCW->getLeftLeg();
	}
	else {
		prevEdge_CCW = edge1_CCW->getRightHand();
	}

	if (edge2_CCW->getLeftFace() == toBeRemovedVFace) {
		nextEdge_CCW = edge2_CCW->getLeftHand();
	}
	else {
		nextEdge_CCW = edge2_CCW->getRightLeg();
	}


	VEdge2D* mergedEdge = createEdge(m_VEdges.back()->getID() + 1);
	if (nextEdge_CCW->getLeftFace() == toBeRemovedVFace) {
		mergedEdge->setStartVertex(nextEdge_CCW->getStartVertex());
		mergedEdge->setLeftLeg(nextEdge_CCW->getLeftLeg());
		mergedEdge->setRightLeg(nextEdge_CCW->getRightLeg());
		mergedEdge->setLeftFace(nextEdge_CCW->getLeftFace());
		mergedEdge->setRightFace(nextEdge_CCW->getRightFace());
	}
	else {
		mergedEdge->setStartVertex(nextEdge_CCW->getEndVertex());
		mergedEdge->setLeftLeg(nextEdge_CCW->getRightHand());
		mergedEdge->setRightLeg(nextEdge_CCW->getLeftHand());
		mergedEdge->setLeftFace(nextEdge_CCW->getRightFace());
		mergedEdge->setRightFace(nextEdge_CCW->getLeftFace());
	}


	if (prevEdge_CCW->getLeftFace() == toBeRemovedVFace) {
		mergedEdge->setEndVertex(prevEdge_CCW->getEndVertex());
		mergedEdge->setLeftHand(prevEdge_CCW->getLeftHand());
		mergedEdge->setRightHand(prevEdge_CCW->getRightHand());
	}
	else {
		mergedEdge->setEndVertex(prevEdge_CCW->getEndVertex());
		mergedEdge->setLeftHand(prevEdge_CCW->getRightLeg());
		mergedEdge->setRightHand(prevEdge_CCW->getLeftLeg());
	}


	mergedEdge->getStartVertex()->setFirstVEdge(mergedEdge);
	mergedEdge->getEndVertex()->setFirstVEdge(mergedEdge);
	mergedEdge->getLeftFace()->setFirstVEdge(mergedEdge);
	mergedEdge->getRightFace()->setFirstVEdge(mergedEdge);


	if (mergedEdge->getLeftLeg()->getLeftHand() == nextEdge_CCW) {
		mergedEdge->getLeftLeg()->setLeftHand(mergedEdge);
	}
	else {
		mergedEdge->getLeftLeg()->setRightLeg(mergedEdge);
	}

	if (mergedEdge->getRightLeg()->getRightHand() == nextEdge_CCW) {
		mergedEdge->getRightLeg()->setRightHand(mergedEdge);
	}
	else {
		mergedEdge->getRightLeg()->setLeftLeg(mergedEdge);
	}

	if (mergedEdge->getLeftHand()->getRightHand() == prevEdge_CCW) {
		mergedEdge->getLeftHand()->setRightHand(mergedEdge);
	}
	else {
		mergedEdge->getLeftHand()->setLeftLeg(mergedEdge);
	}

	if (mergedEdge->getRightHand()->getLeftHand() == prevEdge_CCW) {
		mergedEdge->getRightHand()->setLeftHand(mergedEdge);
	}
	else {
		mergedEdge->getRightHand()->setRightLeg(mergedEdge);
	}


	prevEdge_CCW->setStatus(RED_E);
	nextEdge_CCW->setStatus(RED_E);
}



void VoronoiDiagram2DC::updateVoronoiDiagramWithTemporarilyDeletedGenerators(set<Generator2D*>& trappedGenerators)
{
	//trappedGenerators.sort(compareGeneratorNonAscendingOrder);

	set<Generator2D*>::iterator i_gen;
	for (i_gen = trappedGenerators.begin(); i_gen != trappedGenerators.end(); ++i_gen) {
		Generator2D* currGenerator = *i_gen;
		updateVoronoiDiagramWithNewGenerator(currGenerator);
	}
}



int VoronoiDiagram2DC::computeCircumcirclesForDeletion( VFace2D* prevFace, VFace2D* currFace, VFace2D* nextFace, rg_Circle2D& circumCircle1, rg_Circle2D& circumCircle2 )
{
    rg_Circle2D generator[3] = { prevFace->getGenerator()->getDisk(), currFace->getGenerator()->getDisk(), nextFace->getGenerator()->getDisk() };
    return rg_Circle2D::makeCircumcircle( generator[0], generator[1], generator[2], circumCircle1, circumCircle2);
}



bool VoronoiDiagram2DC::findTrappedGeneratorsAndEntities(VFace2D* VFaceToBeRemoved, VFace2D* trappingFace, VEdge2D* prevEdge_CCW, VEdge2D* nextEdge_CCW, set<Generator2D*>& trappedGenerators)
{
    //prevEdge_CCW->setStatus(RED_E);
    //nextEdge_CCW->setStatus(RED_E);

    bool trappedGeneratorsIncludeInfinite = false;
    list<VVertex2D*>    trappedVertices;

    if ( prevEdge_CCW->getLeftFace() == VFaceToBeRemoved ) {
        trappedVertices.push_back(prevEdge_CCW->getEndVertex());
    }
    else {
        trappedVertices.push_back(prevEdge_CCW->getStartVertex());
    }

    if ( nextEdge_CCW->getLeftFace() == VFaceToBeRemoved ) {
        trappedVertices.push_back(nextEdge_CCW->getStartVertex());
    }
    else {
        trappedVertices.push_back(nextEdge_CCW->getEndVertex());
    }

    list<VEdge2D*>      trappedEdges;
    set<VFace2D*>       trappedFaces;
    VFace2D*            infiniteFace = m_VFaces.front();

    if ( trappingFace == infiniteFace ) {
        infiniteFace = rg_NULL;
    }

    list<VVertex2D*>::iterator i_vtx;
    for ( i_vtx = trappedVertices.begin(); i_vtx != trappedVertices.end(); ++i_vtx ) {
        VVertex2D* currTrappedVertex = *i_vtx;

        if ( currTrappedVertex->getStatus() == RED_V ) {
            continue;
        }
        currTrappedVertex->setStatus(RED_V);

        list<VEdge2D*> incidentEdges;
        currTrappedVertex->getIncident3VEdges(incidentEdges);

        list<VEdge2D*>::iterator i_incidentEdge;
        for ( i_incidentEdge = incidentEdges.begin(); i_incidentEdge != incidentEdges.end(); ++i_incidentEdge ) {
            VEdge2D* currIncidentEdge = *i_incidentEdge;

            if ( currIncidentEdge == prevEdge_CCW || currIncidentEdge == nextEdge_CCW ) {
                continue;
            }
            if ( currIncidentEdge->getStatus() == RED_E ) {
                continue;
            }

            trappedEdges.push_back(currIncidentEdge);

            VVertex2D* currOppositeVertex = currIncidentEdge->getOppositVVertex(currTrappedVertex);
            if ( currOppositeVertex->getStatus() != RED_V ) {
                trappedVertices.push_back(currOppositeVertex);
            }

            VFace2D* currLeftFace  = currIncidentEdge->getLeftFace();
            VFace2D* currRightFace = currIncidentEdge->getRightFace();
            //VFace2D* virtualFace = &(m_VFaces.front());

            if ( currLeftFace == infiniteFace || currRightFace == infiniteFace ) {
                trappedGeneratorsIncludeInfinite = true;
                break;
            }

            if (currLeftFace != VFaceToBeRemoved && currLeftFace != trappingFace) {
                trappedGenerators.insert(currLeftFace->getGenerator());
            }

            if (currRightFace != VFaceToBeRemoved && currRightFace != trappingFace) {
                trappedGenerators.insert(currRightFace->getGenerator());
            }

            currIncidentEdge->setStatus(RED_E);
        }

        if ( trappedGeneratorsIncludeInfinite ) {
            break;
        }
    }



    if ( trappedGeneratorsIncludeInfinite ) {
        list<VVertex2D*>::iterator i_trappedVtx;
        for ( i_trappedVtx = trappedVertices.begin(); i_trappedVtx != trappedVertices.end(); ++i_trappedVtx ) {
            VVertex2D* currTrappedVertex = *i_trappedVtx;
            currTrappedVertex->setStatus(WHITE_V);
        }

        list<VEdge2D*>::iterator i_trappedEdge;
        for ( i_trappedEdge = trappedEdges.begin(); i_trappedEdge != trappedEdges.end(); ++i_trappedEdge ) {
            VEdge2D* currTrappedEdge = *i_trappedEdge;
            currTrappedEdge->setStatus(WHITE_E);
        }

        trappedVertices.clear();
        trappedEdges.clear();
        trappedGenerators.clear();
    }

    return !trappedGeneratorsIncludeInfinite;
}



bool VoronoiDiagram2DC::thereIsVEdgeWhichHasIdenticalDiskTripletBetweenAdjacentVEdgesOnBoundaryOfInputVFace( const VEdge2D * edge, const VFace2D * face, VEdge2D*& prevTrappingEdge, VEdge2D*& nextTrappingEdge )
{
    VEdge2D* prevVEdge      = rg_NULL;
    VEdge2D* prevprevVEdge  = rg_NULL;
    VEdge2D* nextVEdge      = rg_NULL;
    VEdge2D* nextnextVEdge  = rg_NULL;

    prevVEdge       = getCCWPrevEdgeOnFace(const_cast<VEdge2D*>(edge), const_cast<VFace2D*>(face));
    prevprevVEdge   = getCCWPrevEdgeOnFace(prevVEdge, const_cast<VFace2D*>(face));
    nextVEdge       = getCCWNextEdgeOnFace(const_cast<VEdge2D*>(edge), const_cast<VFace2D*>(face));
    nextnextVEdge   = getCCWNextEdgeOnFace(nextVEdge, const_cast<VFace2D*>(face));

    VFace2D* prevprevOppositeVFace  = getFaceOfOppositeSide(const_cast<VFace2D*>(face), prevprevVEdge);
    VFace2D* prevOppositeVFace      = getFaceOfOppositeSide(const_cast<VFace2D*>(face), prevVEdge);
    VFace2D* currOppositeVFace      = getFaceOfOppositeSide(const_cast<VFace2D*>(face), const_cast<VEdge2D*>(edge));
    VFace2D* nextOppositeVFace      = getFaceOfOppositeSide(const_cast<VFace2D*>(face), nextVEdge);
    VFace2D* nextnextOppositeVFace  = getFaceOfOppositeSide(const_cast<VFace2D*>(face), nextnextVEdge);

    if( prevprevOppositeVFace == nextOppositeVFace ) { // if prevVEdge has identical triplet
        prevTrappingEdge = nextVEdge;
        nextTrappingEdge = prevprevVEdge;
        return true;
    }
    else if( nextnextOppositeVFace == prevOppositeVFace ) { // if nextVEdge has identical triplet
        prevTrappingEdge = nextnextVEdge;
        nextTrappingEdge = prevVEdge;
        return true;
    }
    else {
        return false;
    }
}

void VoronoiDiagram2DC::findTrappingVEdges( VFace2D * faceToBeRemoved, list<VEdge2D*>& trappingEdges, VFace2D *& trappingFace )
{
    trappingFace = rg_NULL;

    list<VEdge2D*> boundingVEdges;
    faceToBeRemoved->getBoundaryVEdges(boundingVEdges);

    map<VFace2D*, VEdge2D*> oppositeFace_N_bndEdge;

    list<VEdge2D*>::iterator i_bndEdge = boundingVEdges.begin();
    for ( ; i_bndEdge != boundingVEdges.end(); ++i_bndEdge ) {
        VEdge2D* currBndEdge        = *i_bndEdge;
        VFace2D* currOppositeFace   = getFaceOfOppositeSide(faceToBeRemoved, currBndEdge);

        map<VFace2D*, VEdge2D*>::iterator i_map_oppositeFace_N_bndEdge = oppositeFace_N_bndEdge.find(currOppositeFace);
        if ( i_map_oppositeFace_N_bndEdge != oppositeFace_N_bndEdge.end() ) {
            trappingFace = i_map_oppositeFace_N_bndEdge->first;
            trappingEdges.push_back(i_map_oppositeFace_N_bndEdge->second);
            trappingEdges.push_back(currBndEdge);
            break;
        }
        else {
            oppositeFace_N_bndEdge.insert(make_pair(currOppositeFace, currBndEdge));
        }

        /*if ( trappingFace == rg_NULL ) {
            map<VFace2D*, VEdge2D*>::iterator i_map_oppositeFace_N_bndEdge = oppositeFace_N_bndEdge.find(currOppositeFace);
            if ( i_map_oppositeFace_N_bndEdge != oppositeFace_N_bndEdge.end() ) {
                trappingFace = i_map_oppositeFace_N_bndEdge->first;
                trappingEdges.push_back(i_map_oppositeFace_N_bndEdge->second);
                trappingEdges.push_back(currBndEdge);
            }
            else {
                oppositeFace_N_bndEdge.insert(make_pair(currOppositeFace, currBndEdge));
            }
        }
        else {
            if ( currOppositeFace == trappingFace ) {
                trappingEdges.push_back(currBndEdge); 
            }
        }*/
    }
}



void VoronoiDiagram2DC::temporarilyRemoveTrappedGenerators( set<Generator2D*>& trappedGenerators )
{
    set<Generator2D*>::iterator i_gen;
    for (i_gen = trappedGenerators.begin(); i_gen != trappedGenerators.end(); ++i_gen) {
        VFace2D* currFace = (*i_gen)->getOuterFace();
        m_VFaces.remove(currFace);
        delete currFace;
    }

}



VEdge2D* VoronoiDiagram2DC::mergeTwoTrappingVEdges( VEdge2D * prevQuillEdge_CCW, VEdge2D* nextQuillEdge_CCW, VFace2D* toBeRemovedVFace )
{
    /*VEdge2D* prevEdge_CCW;
    VEdge2D* nextEdge_CCW;

    if (edge1_CCW->getLeftFace() == toBeRemovedVFace) {
        prevEdge_CCW = edge1_CCW->getLeftLeg();
    }
    else {
        prevEdge_CCW = edge1_CCW->getRightHand();
    }

    if (edge2_CCW->getLeftFace() == toBeRemovedVFace) {
        nextEdge_CCW = edge2_CCW->getLeftHand();
    }
    else {
        nextEdge_CCW = edge2_CCW->getRightLeg();
    }*/


    VEdge2D* mergedEdge = createEdge(m_VEdges.back()->getID() + 1);
    if (prevQuillEdge_CCW->getLeftFace() == toBeRemovedVFace) {
        mergedEdge->setStartVertex(prevQuillEdge_CCW->getStartVertex());
        mergedEdge->setLeftLeg(prevQuillEdge_CCW->getLeftLeg());
        mergedEdge->setRightLeg(prevQuillEdge_CCW->getRightLeg());
        mergedEdge->setLeftFace(prevQuillEdge_CCW->getLeftFace());
        mergedEdge->setRightFace(prevQuillEdge_CCW->getRightFace());
    }
    else {
        mergedEdge->setStartVertex(prevQuillEdge_CCW->getEndVertex());
        mergedEdge->setLeftLeg(prevQuillEdge_CCW->getRightHand());
        mergedEdge->setRightLeg(prevQuillEdge_CCW->getLeftHand());
        mergedEdge->setLeftFace(prevQuillEdge_CCW->getRightFace());
        mergedEdge->setRightFace(prevQuillEdge_CCW->getLeftFace());
    }


    if (nextQuillEdge_CCW->getLeftFace() == toBeRemovedVFace) {
        mergedEdge->setEndVertex(nextQuillEdge_CCW->getEndVertex());
        mergedEdge->setLeftHand(nextQuillEdge_CCW->getLeftHand());
        mergedEdge->setRightHand(nextQuillEdge_CCW->getRightHand());
    }
    else {
        mergedEdge->setEndVertex(nextQuillEdge_CCW->getStartVertex());
        mergedEdge->setLeftHand(nextQuillEdge_CCW->getRightLeg());
        mergedEdge->setRightHand(nextQuillEdge_CCW->getLeftLeg());
    }


    mergedEdge->getStartVertex()->setFirstVEdge(mergedEdge);
    mergedEdge->getEndVertex()->setFirstVEdge(mergedEdge);
    mergedEdge->getLeftFace()->setFirstVEdge(mergedEdge);
    mergedEdge->getRightFace()->setFirstVEdge(mergedEdge);


    if (mergedEdge->getLeftLeg()->getLeftHand() == prevQuillEdge_CCW) {
        mergedEdge->getLeftLeg()->setLeftHand(mergedEdge);
    }
    else {
        mergedEdge->getLeftLeg()->setRightLeg(mergedEdge);
    }

    if (mergedEdge->getRightLeg()->getRightHand() == prevQuillEdge_CCW) {
        mergedEdge->getRightLeg()->setRightHand(mergedEdge);
    }
    else {
        mergedEdge->getRightLeg()->setLeftLeg(mergedEdge);
    }

    if (mergedEdge->getLeftHand()->getRightHand() == nextQuillEdge_CCW) {
        mergedEdge->getLeftHand()->setRightHand(mergedEdge);
    }
    else {
        mergedEdge->getLeftHand()->setLeftLeg(mergedEdge);
    }

    if (mergedEdge->getRightHand()->getLeftHand() == nextQuillEdge_CCW) {
        mergedEdge->getRightHand()->setLeftHand(mergedEdge);
    }
    else {
        mergedEdge->getRightHand()->setRightLeg(mergedEdge);
    }


    prevQuillEdge_CCW->setStatus(RED_E);
    nextQuillEdge_CCW->setStatus(RED_E);

    return mergedEdge;
}

void VoronoiDiagram2DC::reflectMergedEdge( VEdge2D * mergedVEdge, VFace2D* toBeRemovedVFace, rg_dList<pair<VEdge2D*, rg_Circle2D> >& possiblyFlippingEdges )
{
    rg_Circle2D circumcircle;
    bool isFlippingEdge = isFlippingEdgeOnVirtualRegion( mergedVEdge, circumcircle, toBeRemovedVFace );

    if( isFlippingEdge ) {
        possiblyFlippingEdges.add( pair<VEdge2D*, rg_Circle2D>(mergedVEdge, circumcircle) );
    }
}



//rg_dList< pair<VEdge2D*, rg_Circle2D> > possiblyFlippingEdges;
//collectPossiblyFlippingEdgesOnTheInterface_StoredInPriorityQueue( virtualFace, possiblyFlippingEdges );
//
//rg_dList<VEdge2D*> edgesToUpdateGeometry;
//
//while( possiblyFlippingEdges.getSize() > 0 ) {
//    pair<VEdge2D*, rg_Circle2D> edgeToDefineMinTangentCircle = dequeueRootNode_DiskTripletWithMinTangentCircle( possiblyFlippingEdges );
//
//    VEdge2D* edgeOfIncidentTriplet[2] = {NULL, NULL};
//    flipThisEdgeNFindTwoIncidentTriplets( virtualFace, edgeToDefineMinTangentCircle, edgeOfIncidentTriplet[0], edgeOfIncidentTriplet[1] );
//
//    for( int j = 0 ; j < 2 ; j++ ) {
//        reflectTheFlipOnTheIncidentTriplet( possiblyFlippingEdges, edgeOfIncidentTriplet[j]);
//    }
//}


//VEdge2D* VoronoiDiagram2DC::findEdgeToBeFlippedBetweenTwoEdgesWithIdenticalDiskTriplet(VFace2D* faceToBeRemoved, VEdge2D* edge1, VEdge2D* edge2)
//{
//	VEdge2D* quillEdge1 = rg_NULL;
//	VEdge2D* quillEdge2 = rg_NULL;
//	VEdge2D* quillEdge3 = rg_NULL;
//
//	VFace2D* oppositeFaceOverEdge1 = rg_NULL;
//	VFace2D* oppositeFaceOverEdge2 = rg_NULL;
//	VFace2D* oppositeFaceDefinedByOtherDiskTriplet = rg_NULL;
//
//	VEdge2D* prevEdgeOfEdge1_CCW = rg_NULL;
//	VEdge2D* nextEdgeOfEdge2_CCW = rg_NULL;
//
//	rg_Point2D endPointOfQuillEdge1;
//	rg_Point2D endPointOfQuillEdge2;
//	rg_Point2D endPointOfQuillEdge3;
//
//	rg_Point2D endPointOfNextEdge;
//	rg_Point2D endPointOfPrevEdge
//}



void VoronoiDiagram2DC::copyVoronoiDiagramFrom(const VoronoiDiagram2DC& VD2DC)
{
	clear();

	list<Generator2D*>       allGenerators_Input;
    list<Generator2D*>       allPhantomGenerators_Input;
	list<VFace2D*>           allFaces_Input;
	list<VEdge2D*>           allEdges_Input;
	list<VVertex2D*>         allVertices_Input;

    getAllVEntities(VD2DC, allGenerators_Input, allPhantomGenerators_Input, allFaces_Input, allEdges_Input, allVertices_Input);

	unordered_map<Generator2D*, Generator2D*>             generatorMapFromOldToNew;
	unordered_map<VFace2D*, VFace2D*>                     faceMapFromOldToNew;
	unordered_map<VEdge2D*, VEdge2D*>                     edgeMapFromOldToNew;
	unordered_map<VVertex2D*, VVertex2D*>                 vertexMapFromOldToNew;

    createNLinkEntities(allGenerators_Input, allPhantomGenerators_Input, allFaces_Input, allEdges_Input,
        allVertices_Input, generatorMapFromOldToNew, faceMapFromOldToNew, edgeMapFromOldToNew, vertexMapFromOldToNew);

    insertCreatedEntitiesToMemberData(allGenerators_Input, allPhantomGenerators_Input, allFaces_Input, allEdges_Input,
                                     allVertices_Input, generatorMapFromOldToNew, faceMapFromOldToNew, edgeMapFromOldToNew, vertexMapFromOldToNew);

    connectCreatedEntities(allGenerators_Input,  allPhantomGenerators_Input,  allFaces_Input,  allEdges_Input,
                           allVertices_Input, generatorMapFromOldToNew, faceMapFromOldToNew, edgeMapFromOldToNew, vertexMapFromOldToNew);

    copyTopologyUnrelatedObjects(VD2DC);
}



void VoronoiDiagram2DC::copyTopologyUnrelatedObjects(const VoronoiDiagram2DC& VD2DC)
{
    m_bucket                                        = VD2DC.m_bucket;
    m_circleEnclosingVertices                       = VD2DC.m_circleEnclosingVertices;
    m_phantomCircles                                = VD2DC.m_phantomCircles;

    m_directAddressTableFromDiskIDToDiskArrangement = VD2DC.m_directAddressTableFromDiskIDToDiskArrangement; // disk arragement in double precision by Joonghyun March, 17
}


void VoronoiDiagram2DC::getAllVEntities(const VoronoiDiagram2DC& VD2DC, list<Generator2D*>& generators, list<Generator2D*>& phantomGenerators, list<VFace2D*>& faces, list<VEdge2D*>& edges, list<VVertex2D*>& vertices)
{
    VD2DC.getGenerators(generators);
    VD2DC.getPhantomGenerator(generators);
    VD2DC.getVoronoiFaces(faces);
    VD2DC.getVoronoiEdges(edges);
    VD2DC.getVoronoiVertices(vertices);
}


void VoronoiDiagram2DC::createNLinkEntities(
    const list<Generator2D*>& oldGenerators, 
    const list<Generator2D*>& oldPhantomGenerators, 
    const list<VFace2D*>& oldFaces, 
    const list<VEdge2D*>& oldEdges, 
    const list<VVertex2D*>& oldVertices, 
    unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew, 
    unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew, 
    unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew, 
    unordered_map<VVertex2D*, VVertex2D*>& vertexMapFromOldToNew) const
{
    createNLinkGenerators(oldGenerators,        generatorMapFromOldToNew);
    createNLinkGenerators(oldPhantomGenerators, generatorMapFromOldToNew);
    createNLinkFaces(oldFaces, faceMapFromOldToNew);
    createNLinkEdges(oldEdges, edgeMapFromOldToNew);
    createNLinkVertices(oldVertices, vertexMapFromOldToNew);
}


void VoronoiDiagram2DC::createNLinkGenerators(const list<Generator2D*>& generators, unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew) const
{
	for (list<Generator2D*>::const_iterator i_gen = generators.begin(); i_gen != generators.end(); ++i_gen)
	{
		Generator2D* currGen = *i_gen;
		Generator2D* newGen  = new Generator2D(*currGen);
        generatorMapFromOldToNew.insert(pair<Generator2D*, Generator2D*>(currGen, newGen));
	}
    generatorMapFromOldToNew.insert(pair<Generator2D*, Generator2D*>(NULL, NULL));
}


void VoronoiDiagram2DC::createNLinkFaces(const list<VFace2D*>& faces, unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew) const
{
	for (list<VFace2D*>::const_iterator i_face = faces.begin(); i_face != faces.end(); ++i_face)
	{
		VFace2D* currFace = *i_face;
		VFace2D* newFace  = new VFace2D(*currFace);
        faceMapFromOldToNew.insert(pair<VFace2D*, VFace2D*>(currFace, newFace));
	}
    faceMapFromOldToNew.insert(pair<VFace2D*, VFace2D*>(NULL, NULL));
}


void VoronoiDiagram2DC::createNLinkEdges(const list<VEdge2D*>& edges, unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew) const
{
	for (list<VEdge2D*>::const_iterator i_edge = edges.begin(); i_edge != edges.end(); ++i_edge)
	{
		VEdge2D* currEdge = *i_edge;
		VEdge2D* newEdge = new VEdge2D(*currEdge);
        edgeMapFromOldToNew.insert(pair<VEdge2D*, VEdge2D*>(currEdge, newEdge));
	}
    edgeMapFromOldToNew.insert(pair<VEdge2D*, VEdge2D*>(NULL, NULL));
}


void VoronoiDiagram2DC::createNLinkVertices(const list<VVertex2D*>& vertices, unordered_map<VVertex2D*, VVertex2D*>& vertexMapFromOldToNew) const
{
	for (list<VVertex2D*>::const_iterator i_vtx = vertices.begin(); i_vtx != vertices.end(); ++i_vtx)
	{
		VVertex2D* currVtx = *i_vtx;
		VVertex2D* newVtx = new VVertex2D(*currVtx);
        vertexMapFromOldToNew.insert(pair<VVertex2D*, VVertex2D*>(currVtx, newVtx));
	}
    vertexMapFromOldToNew.insert(pair<VVertex2D*, VVertex2D*>(NULL, NULL));
}


void VoronoiDiagram2DC::insertCreatedEntitiesToMemberData(
    const list<Generator2D*>& oldGenerators, 
    const list<Generator2D*>& oldPhantomGenerators, 
    const list<VFace2D*>& oldFaces, 
    const list<VEdge2D*>& oldEdges, 
    const list<VVertex2D*>& oldVertices, 
    const unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew, 
    const unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew, 
    const unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew, 
    const unordered_map<VVertex2D*, VVertex2D*>& vertexMapFromOldToNew)
{
    insertCreatedGenerators(       oldGenerators,         generatorMapFromOldToNew);
    insertCreatedPhantomGenerators(oldPhantomGenerators,  generatorMapFromOldToNew);
    insertCreatedFaces(            oldFaces,              faceMapFromOldToNew);
    insertCreatedEdges(            oldEdges,              edgeMapFromOldToNew);
    insertCreatedVertices(         oldVertices,           vertexMapFromOldToNew);
}


void VoronoiDiagram2DC::insertCreatedPhantomGenerators(const list<Generator2D*>& oldPhantomGenerators, const unordered_map<Generator2D*, Generator2D*>& phantomGeneratorMapFromOldToNew)
{
	for (list<Generator2D*>::const_iterator i_gen = oldPhantomGenerators.begin(); i_gen != oldPhantomGenerators.end(); ++i_gen)
	{
		Generator2D* currGen = *i_gen;
        m_phantomGenerators.push_back(phantomGeneratorMapFromOldToNew.at(currGen));
	}
}


void VoronoiDiagram2DC::insertCreatedGenerators(const list<Generator2D*>& oldGenerators, const unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew)
{
	for (list<Generator2D*>::const_iterator i_gen = oldGenerators.begin(); i_gen != oldGenerators.end(); ++i_gen)
	{
		Generator2D* currGen = *i_gen;
		m_generators.push_back(generatorMapFromOldToNew.at(currGen));
	}
}


void VoronoiDiagram2DC::insertCreatedFaces(const list<VFace2D*>& oldFaces, const unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew)
{
	for (list<VFace2D*>::const_iterator i_face = oldFaces.begin(); i_face != oldFaces.end(); ++i_face)
	{
		VFace2D* currFace = *i_face;
		m_VFaces.push_back(faceMapFromOldToNew.at(currFace));
	}
}


void VoronoiDiagram2DC::insertCreatedEdges(const list<VEdge2D*>& oldEdges, const unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew)
{
	for (list<VEdge2D*>::const_iterator i_edge = oldEdges.begin(); i_edge != oldEdges.end(); ++i_edge)
	{
        VEdge2D* currEdge = *i_edge;
		m_VEdges.push_back(edgeMapFromOldToNew.at(currEdge));
	}
}


void VoronoiDiagram2DC::insertCreatedVertices(const list<VVertex2D*>& oldVertices, const unordered_map<VVertex2D*, VVertex2D*>& vertexMapFromOldToNew)
{
	for (list<VVertex2D*>::const_iterator i_vtx = oldVertices.begin(); i_vtx != oldVertices.end(); ++i_vtx)
	{
		VVertex2D* currVtx = *i_vtx;
		m_VVertices.push_back(vertexMapFromOldToNew.at(currVtx));
	}
}


void VoronoiDiagram2DC::connectCreatedEntities(const list<Generator2D*>& oldGenerators, 
                                               const list<Generator2D*>& oldPhantomGenerators, 
                                               const list<VFace2D*>& oldFaces, 
                                               const list<VEdge2D*>& oldEdges, 
                                               const list<VVertex2D*>& oldVertices, 
                                               const unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew, 
                                               const unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew, 
                                               const unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew, 
                                               const unordered_map<VVertex2D*, VVertex2D*>& vertexMapFromOldToNew) const
{
	connectCreatedEntitiesInGenerators(oldGenerators,        generatorMapFromOldToNew, faceMapFromOldToNew);
    connectCreatedEntitiesInGenerators(oldPhantomGenerators, generatorMapFromOldToNew, faceMapFromOldToNew);
	connectCreatedEntitiesInFaces(oldFaces,                  generatorMapFromOldToNew, faceMapFromOldToNew, edgeMapFromOldToNew);
	connectCreatedEntitiesInEdges(oldEdges,                  faceMapFromOldToNew,      edgeMapFromOldToNew, vertexMapFromOldToNew);
	connectCreatedEntitiesInVertices(oldVertices,            edgeMapFromOldToNew,      vertexMapFromOldToNew);
}


void VoronoiDiagram2DC::connectCreatedEntitiesInGenerators(const list<Generator2D*>& oldGenerators, 
                                                           const unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew, 
                                                           const unordered_map<VFace2D*, VFace2D*>&         faceMapFromOldToNew) const
{
	for (list<Generator2D*>::const_iterator i_gen = oldGenerators.begin();
		 i_gen != oldGenerators.end();
		 ++i_gen)
	{
		Generator2D* currOldGen = *i_gen;
		Generator2D* currNEWGen = generatorMapFromOldToNew.at(currOldGen);
        currNEWGen->setOuterFace(faceMapFromOldToNew.at(currOldGen->getOuterFace()));
	}
}


void VoronoiDiagram2DC::connectCreatedEntitiesInFaces(
    const list<VFace2D*>& oldFaces, 
    const unordered_map<Generator2D*, Generator2D*>& generatorMapFromOldToNew,
    const unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew,
    const unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew) const
{
	for (list<VFace2D*>::const_iterator i_face = oldFaces.begin(); i_face != oldFaces.end(); ++i_face)
	{
		VFace2D* currOldFace = *i_face;
        VFace2D* currNEWFace = faceMapFromOldToNew.at(currOldFace);

        currNEWFace->setFirstVEdge(edgeMapFromOldToNew.at(currOldFace->getFirstVEdge()));
        currNEWFace->setGenerator(generatorMapFromOldToNew.at(currOldFace->getGenerator()));
	}
}


void VoronoiDiagram2DC::connectCreatedEntitiesInEdges(
    const list<VEdge2D*>& oldEdges, 
    const unordered_map<VFace2D*, VFace2D*>& faceMapFromOldToNew, 
    const unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew, 
    const unordered_map<VVertex2D*, VVertex2D*>& vertexMapFromOldToNew) const
{
     for (list<VEdge2D*>::const_iterator i_edge = oldEdges.begin(); i_edge != oldEdges.end(); ++i_edge)
    {
        VEdge2D* currOldEdge    = *i_edge;
        VEdge2D* currNEWEdge    = edgeMapFromOldToNew.at(currOldEdge);

        currNEWEdge->setStartVertex(vertexMapFromOldToNew.at( currOldEdge->getStartVertex()));
        currNEWEdge->setEndVertex(vertexMapFromOldToNew.at(currOldEdge->getEndVertex()));

        currNEWEdge->setRightHand(edgeMapFromOldToNew.at(currOldEdge->getRightHand()));
        currNEWEdge->setLeftHand(edgeMapFromOldToNew.at(currOldEdge->getLeftHand()));
        currNEWEdge->setRightLeg(edgeMapFromOldToNew.at(currOldEdge->getRightLeg()));
        currNEWEdge->setLeftLeg(edgeMapFromOldToNew.at(currOldEdge->getLeftLeg()));

        currNEWEdge->setRightFace(faceMapFromOldToNew.at(currOldEdge->getRightFace()));
        currNEWEdge->setLeftFace(faceMapFromOldToNew.at(currOldEdge->getLeftFace()));
    }
}


void VoronoiDiagram2DC::connectCreatedEntitiesInVertices(
    const list<VVertex2D*>& oldVertices, 
    const unordered_map<VEdge2D*, VEdge2D*>& edgeMapFromOldToNew, 
    const unordered_map<VVertex2D*, VVertex2D*>& vertexMapFromOldToNew) const
{
	for (list<VVertex2D*>::const_iterator i_vtx = oldVertices.begin(); i_vtx != oldVertices.end(); ++i_vtx)
	{
		VVertex2D* currOldVertex = *i_vtx;
		VVertex2D* currNEWVertex = vertexMapFromOldToNew.at(currOldVertex);
		currNEWVertex->setFirstVEdge(edgeMapFromOldToNew.at(currOldVertex->getFirstVEdge()));
	}
}


VoronoiDiagram2DC& VoronoiDiagram2DC::operator=(const VoronoiDiagram2DC& VD2DC)
{
    if (this == &VD2DC){
        return *this;
    }

    clear();
    copyVoronoiDiagramFrom(VD2DC);

    return *this;
}



bool VoronoiDiagram2DC::operator==(const VoronoiDiagram2DC& VD2DC) const
{
    // 0_1. Check whether faces and generators have right connections
    if (!does_vfaces_N_generators_have_right_connection() || !VD2DC.does_vfaces_N_generators_have_right_connection())
    {
        return false;
    }


    // 0_2. Check number of generators // VVertex // VEdges // VFaces 
    if (getNumberOfVVertices() != VD2DC.getNumberOfVVertices() ||
        getNumberOfVEdges() != VD2DC.getNumberOfVEdges() ||
        getNumberOfVFaces() != VD2DC.getNumberOfVFaces())
    {
        return false;
    }


    // 1. Generator Mapping : by comparing center and radius
    //if fail, no same
    map < Generator2D*, Generator2D* > linker_gen_to_gen;

    list<Generator2D*> generators_VD1;
    list<Generator2D*> generators_VD2;

    getGenerators(generators_VD1);
    getPhantomGenerator(generators_VD1);

    VD2DC.getGenerators(generators_VD2);
    VD2DC.getPhantomGenerator(generators_VD2);


    if (!are_same_generators(generators_VD1, generators_VD2, linker_gen_to_gen))
    {
        return false;
    }

    // 2. Compare face topology :
    //get a face from each VD2DC, check the neighbor generators with order
    list<VFace2D*> VFaces1;
    list<VFace2D*> VFaces2;

    getVoronoiFaces(VFaces1);
    VD2DC.getVoronoiFaces(VFaces2);

    if (!are_same_vface_topologies(VFaces1, VFaces2, linker_gen_to_gen))
    {
        return false;
    }


    // 3. Compare edge topology :
    // get a edge from each VD2DC, check the both side's generators
    //1) edge Î¶¨Ïä§?ÔøΩÔøΩ? Ôø??VD2DC?ÔøΩÏÑú Ôø??ÔøΩÏò®??
    //2) Í∞ôÔøΩ? ÏßùÏù¥ ?ÔøΩÎäîÔø?Ôø??Ï∞æÎäî??: edge?Ôø? Ôø???????ÔøΩÏùÑ ?ÔøΩÏùò?ÔøΩÎäî generatorÔø??ÎπÑÍµê

    list<VEdge2D*> VEdges1;
    list<VEdge2D*> VEdges2;

    getVoronoiEdges(VEdges1);
    VD2DC.getVoronoiEdges(VEdges2);

    if (!are_same_vedge_topologies(VEdges1, VEdges2, linker_gen_to_gen))
    {
        return false;
    }


    // 4. Compare vertex topology :
    // get a vertex from each VD2DC, check the neigibor generators with order
    //1) vertex Î¶¨Ïä§?ÔøΩÔøΩ? Ôø??VD2DC?ÔøΩÏÑú Ôø??ÔøΩÏò®??
    //2) Í∞ôÔøΩ? ÏßùÏù¥ ?ÔøΩÎäîÔø?Ôø??Ï∞æÎäî?? : vertex??neighbor???ÔøΩÎäî face?ÔøΩÏùò generatorÔø??ÎπÑÍµê

    list<VVertex2D*> VVertices1;
    list<VVertex2D*> VVertices2;

    getVoronoiVertices(VVertices1);
    VD2DC.getVoronoiVertices(VVertices2);

    if (!are_same_vvertex_topologies_new(VVertices1, VVertices2, linker_gen_to_gen))
    {
        return false;
    }


    return true;
}



bool VoronoiDiagram2DC::operator!=(const VoronoiDiagram2DC& VD2DC) const
{
    return !operator==(VD2DC);
}





bool VoronoiDiagram2DC::does_vfaces_N_generators_have_right_connection() const
{
    list<Generator2D*> generators;
    getGenerators(generators);

    //check if vface has a right generator	
    for (list<Generator2D*>::const_iterator iter_gen = generators.begin();
        iter_gen != generators.end();
        iter_gen++)
    {
        Generator2D* gen_a = *iter_gen;
        Generator2D* gen_b = gen_a->getOuterFace()->getGenerator();

        if (gen_a != gen_b)
        {
            return false;
        }
    }

    list<Generator2D*> phantomGenerators;
    getPhantomGenerator(phantomGenerators);

    //check if vface has a right generator	
    for (list<Generator2D*>::const_iterator it_PhantomGenerator = phantomGenerators.begin();
        it_PhantomGenerator != phantomGenerators.end();
        it_PhantomGenerator++)
    {
        Generator2D* gen_a = *it_PhantomGenerator;
        Generator2D* gen_b = gen_a->getOuterFace()->getGenerator();

        if (gen_a != gen_b)
        {
            return false;
        }
    }

    return true;
}



bool VoronoiDiagram2DC::are_same_generators(const list<Generator2D*>& generators1, const list<Generator2D*>& generators2, map < Generator2D*, Generator2D* >& linker_gen_to_gen) const
{
    for (list<Generator2D*>::const_iterator iter_gen1 = generators1.begin();
        iter_gen1 != generators1.end();
        iter_gen1++)
    {
        list<Generator2D*>::const_iterator iter_gen2 = generators2.begin();

        while (iter_gen2 != generators2.end())
        {
            Generator2D* gen1 = *iter_gen1;
            Generator2D* gen2 = *iter_gen2;

            if ((gen1->getDisk()) == (gen2->getDisk()))
            {
                linker_gen_to_gen.insert(pair<Generator2D*, Generator2D*>(gen1, gen2));
                break;
            }
            else
            {
                iter_gen2++;
            }
        }

        if (iter_gen2 == generators2.end())
        {
            return false;
        }
    }

    //for a virtual vface
    linker_gen_to_gen.insert(pair<Generator2D*, Generator2D*>(nullptr, nullptr));

    return true;
}



bool VoronoiDiagram2DC::are_same_vface_topologies(const list<VFace2D*>& VFaces1, const list<VFace2D*>& VFaces2, const map < Generator2D*, Generator2D* >& linker_gen_to_gen) const
{
    //1)Ôø??VD??faceÔø??Î∞õÏïÑ?ÔøΩÎã§.	
    //2)face??generatorÔø??ÎπÑÍµê?ÔøΩÏó¨ Í∞ôÔøΩ? faceÔø??Ï∞æÎäî?? 
    // - face??generatorÔø? Í∞ôÏúºÔø?? ?ÔøΩÏùåÔø??Í∞ôÔøΩ? faceÔø??Ï∞æÎäî?? (Î™®Îì† faceÔø? Í∞ôÏïÑ???? 
    //		- face Ï£ºÔøΩ???face??generator?ÔøΩÏù¥ Í∞ôÔøΩ?Ôø?Ôø???ÔøΩÎã®?ÔøΩÎã§.
    //			- Í∞ôÔøΩ? faceÔø??Ï∞æÎäî??
    //			- Í∞ôÔøΩ? faceÔø????ÔøΩÏûë?ÔøΩÏÑú, Ï£ºÔøΩ? face??Î™®Îì† faceÔø? Í∞ôÔøΩ?Ôø? ?ÔøΩÎã®?ÔøΩÎã§
    //			- Î™®Îëê Í∞ôÏúºÔø???ÔøΩÏùå face 2Í∞úÔøΩ? Ï∞æÏïÑ Î∞òÎ≥µ, ?ÔøΩÎ•¥Ôø???ÔøΩÏÉÅ???ÔøΩÎ•¥?ÔøΩÍ≥† ?ÔøΩÎã®?ÔøΩÎã§. 
    // - ?ÔøΩÏùåÔø??faceÔø? Í∞ôÔøΩ? ?ÔøΩÏúºÔø?? topologyÔø? ?ÔøΩÎ•¥?ÔøΩÍ≥† ?ÔøΩÎã®?ÔøΩÎã§.


    for (list<VFace2D*>::const_iterator iter_face1 = VFaces1.begin();
        iter_face1 != VFaces1.end();
        iter_face1++)
    {
        list<VFace2D*>::const_iterator iter_face2 = VFaces2.begin();

        while (iter_face2 != VFaces2.end())
        {
            VFace2D* face1 = *iter_face1;
            VFace2D* face2 = *iter_face2;

            //check whether two are same generators
            if (linker_gen_to_gen.at(face1->getGenerator()) != face2->getGenerator())
            {
                iter_face2++;
                continue;
            }

            //check whether container or not
            if ( ( (VoronoiDiagram2DC*)face1->getGenerator()->getInnerVD() == this ) || ( (VoronoiDiagram2DC*)face2->getGenerator()->getInnerVD() == this ) )
            {
                break;
            }
            else
            {
                list<VFace2D*> adjacentVFaces1;
                face1->getAdjacentVFaces(adjacentVFaces1); //clock-wise

                list<VFace2D*> adjacentVFaces2;
                face2->getAdjacentVFaces(adjacentVFaces2); //clock-wise

                list<VFace2D*>::const_iterator iter_firstSameVFace = adjacentVFaces2.end();

                //find first same generator
                for (list<VFace2D*>::const_iterator iter_adjFace2 = adjacentVFaces2.begin();
                    iter_adjFace2 != adjacentVFaces2.end();
                    iter_adjFace2++)
                {
                    VFace2D* adjVFace1 = adjacentVFaces1.front();
                    VFace2D* adjVFace2 = *iter_adjFace2;

                    if (linker_gen_to_gen.at(adjVFace1->getGenerator()) == adjVFace2->getGenerator())
                    {
                        iter_firstSameVFace = iter_adjFace2;
                        break;
                    }
                }

                //if there is no same VFace
                if (iter_firstSameVFace == adjacentVFaces2.end())
                {
                    return false;
                }


                //from the first same generator, check whether the remaining generator is same
                list<VFace2D*>::const_iterator iter_adjFace1 = ++adjacentVFaces1.begin();
                list<VFace2D*>::const_iterator iter_adjFace2 = iter_firstSameVFace;
                iter_adjFace2++;

                while (iter_adjFace2 != adjacentVFaces2.end())
                {
                    VFace2D* adjVFace1 = *(iter_adjFace1);
                    VFace2D* adjVFace2 = *(iter_adjFace2);

                    if (linker_gen_to_gen.at(adjVFace1->getGenerator()) != adjVFace2->getGenerator())
                    {
                        return false;
                    }
                    else
                    {
                        iter_adjFace1++;
                        iter_adjFace2++;
                    }
                }


                if (adjacentVFaces2.begin() != iter_firstSameVFace)
                {
                    iter_adjFace2 = adjacentVFaces2.begin();

                    while (iter_adjFace2 != iter_firstSameVFace)
                    {
                        VFace2D* adjVFace1 = *(iter_adjFace1);
                        VFace2D* adjVFace2 = *(iter_adjFace2);

                        if (linker_gen_to_gen.at(adjVFace1->getGenerator()) != adjVFace2->getGenerator())
                        {
                            return false;
                        }
                        else
                        {
                            iter_adjFace1++;
                            iter_adjFace2++;
                        }
                    }
                }


                break;
            }
        }

        if (iter_face2 == VFaces2.end())
        {
            return false;
        }
    }

    return true;
}



bool VoronoiDiagram2DC::are_same_vedge_topologies(const list<VEdge2D*>& VEdges1, const list<VEdge2D*>& VEdges2, const map < Generator2D*, Generator2D* >& linker_gen_to_gen) const
{
    for (list<VEdge2D*>::const_iterator iter_edge1 = VEdges1.begin();
        iter_edge1 != VEdges1.end();
        iter_edge1++)
    {
        list<VEdge2D*>::const_iterator iter_edge2 = VEdges2.begin();

        while (iter_edge2 != VEdges2.end())
        {
            VEdge2D* edge1 = *iter_edge1;
            VEdge2D* edge2 = *iter_edge2;

            if (linker_gen_to_gen.at(edge1->getLeftFace()->getGenerator()) != edge2->getLeftFace()->getGenerator() &&
                linker_gen_to_gen.at(edge1->getRightFace()->getGenerator()) != edge2->getLeftFace()->getGenerator())
            {
                iter_edge2++;
                continue;
            }

            if (linker_gen_to_gen.at(edge1->getLeftFace()->getGenerator()) != edge2->getRightFace()->getGenerator() &&
                linker_gen_to_gen.at(edge1->getRightFace()->getGenerator()) != edge2->getRightFace()->getGenerator())
            {
                iter_edge2++;
                continue;
            }

            if (linker_gen_to_gen.at(edge1->getMateFace(edge1->getStartVertex())->getGenerator()) != edge2->getMateFace(edge2->getStartVertex())->getGenerator() &&
                linker_gen_to_gen.at(edge1->getMateFace(edge1->getEndVertex())->getGenerator()) != edge2->getMateFace(edge2->getStartVertex())->getGenerator())
            {
                iter_edge2++;
                continue;
            }

            if (linker_gen_to_gen.at(edge1->getMateFace(edge1->getStartVertex())->getGenerator()) != edge2->getMateFace(edge2->getEndVertex())->getGenerator() &&
                linker_gen_to_gen.at(edge1->getMateFace(edge1->getEndVertex())->getGenerator()) != edge2->getMateFace(edge2->getEndVertex())->getGenerator())
            {
                iter_edge2++;
                continue;
            }

            break;
        }

        if (iter_edge2 == VEdges2.end())
        {
            return false;
        }
    }

    return true;
}



bool VoronoiDiagram2DC::are_same_vvertex_topologies(const list<VVertex2D*>& VVertices1, const list<VVertex2D*>& VVertices2, const map < Generator2D*, Generator2D* >& linker_gen_to_gen) const
{
    for (list<VVertex2D*>::const_iterator iter_vertex1 = VVertices1.begin();
        iter_vertex1 != VVertices1.end();
        iter_vertex1++)
    {
        list<VVertex2D*>::const_iterator iter_vertex2 = VVertices2.begin();

        while (iter_vertex2 != VVertices2.end())
        {
            VVertex2D* vertex1 = *iter_vertex1;
            VVertex2D* vertex2 = *iter_vertex2;

            list<VFace2D*> incidentVFaces1;
            vertex1->getIncident3VFaces(incidentVFaces1); //clock-wise

            list<VFace2D*> incidentVFaces2;
            vertex2->getIncident3VFaces(incidentVFaces2); //clock-wise

            list<VFace2D*>::const_iterator iter_firstSameVFace = incidentVFaces2.end();

            //find first same generator
            for (list<VFace2D*>::const_iterator iter_incFace2 = incidentVFaces2.begin();
                iter_incFace2 != incidentVFaces2.end();
                iter_incFace2++)
            {
                VFace2D* incVFace1 = incidentVFaces1.front();
                VFace2D* incVFace2 = *iter_incFace2;

                if (linker_gen_to_gen.at(incVFace1->getGenerator()) == incVFace2->getGenerator())
                {
                    iter_firstSameVFace = iter_incFace2;
                    break;
                }
            }

            //if there is no same VFace
            if (iter_firstSameVFace == incidentVFaces2.end())
            {
                iter_vertex2++;
                continue;
            }


            //from the first same generator, check whether the remaining generator is same
            // generator 3Í∞úÔøΩ? ?ÔøΩÍ∞ô?Ôø? Í±∞Î©¥, ?ÔøΩÏÑúÍπåÔøΩ? ?ÔøΩÏù∏?ÔøΩÏïº ??!
            // 
            bool areSameIncidentVFacesUntilNow = true;

            list<VFace2D*>::const_iterator iter_incFace1 = ++incidentVFaces1.begin();
            list<VFace2D*>::const_iterator iter_incFace2 = iter_firstSameVFace;
            iter_incFace2++;

            while (iter_incFace2 != incidentVFaces2.end())
            {
                VFace2D* incVFace1 = *(iter_incFace1);
                VFace2D* incVFace2 = *(iter_incFace2);

                if (linker_gen_to_gen.at(incVFace1->getGenerator()) != incVFace2->getGenerator())
                {
                    areSameIncidentVFacesUntilNow = false;
                    break;
                }
                else
                {
                    iter_incFace1++;
                    iter_incFace2++;
                }
            }


            //if there is no same VFace
            if (!areSameIncidentVFacesUntilNow)
            {
                iter_vertex2++;
                continue;
            }


            if (incidentVFaces2.begin() != iter_firstSameVFace)
            {
                iter_incFace2 = incidentVFaces2.begin();

                while (iter_incFace2 != iter_firstSameVFace)
                {
                    VFace2D* incVFace1 = *(iter_incFace1);
                    VFace2D* incVFace2 = *(iter_incFace2);

                    if (linker_gen_to_gen.at(incVFace1->getGenerator()) != incVFace2->getGenerator())
                    {
                        areSameIncidentVFacesUntilNow = false;
                        break;
                    }
                    else
                    {
                        iter_incFace1++;
                        iter_incFace2++;
                    }
                }

                //if there is no same VFace
                if (!areSameIncidentVFacesUntilNow)
                {
                    iter_vertex2++;
                    continue;
                }
            }

            break;
        }

        if (iter_vertex2 == VVertices2.end())
        {
            return false;
        }
    }

    return true;
}




bool VoronoiDiagram2DC::are_same_vvertex_topologies_new(const list<VVertex2D*>& VVertices1, const list<VVertex2D*>& VVertices2, const map < Generator2D*, Generator2D* >& linker_gen_to_gen) const
{

    //vface 3Ôø??generatorÔø? Î™®Îëê Í∞ôÔøΩ?Ôø?Ôø???ÔøΩÏù∏
    // Î™®Îëê Í∞ôÏúºÔø???ÔøΩÏÑúÔø???ÔøΩÏù∏ 
    //  - ?ÔøΩÏÑúÔø? Í∞ôÏúºÔø??VVertices ?ÔøΩÍ∞ú???ÔøΩÏùåÍ±∏Î°ú ?ÔøΩÏñ¥Ôø??
    //  - ?ÔøΩÏÑúÔø? ?ÔøΩÎ•¥Ôø??2Î≤àÏß∏ VVertices ?ÔøΩÏùåÍ≤ÉÏúºÔø???ÔøΩÏñ¥Ôø??
    // ?ÔøΩÎÇò?ÔøΩÎèÑ ?ÔøΩÎ•¥Ôø??2Î≤àÏß∏ VVertices ?ÔøΩÏùåÍ≤ÉÏúºÔø???ÔøΩÏñ¥Ôø??

    for (list<VVertex2D*>::const_iterator iter_vertex1 = VVertices1.begin();
        iter_vertex1 != VVertices1.end();
        iter_vertex1++)
    {
        list<VVertex2D*>::const_iterator iter_vertex2 = VVertices2.begin();

        while (iter_vertex2 != VVertices2.end())
        {
            bool doesFindSameVtx = false;

            VVertex2D* vertex1 = *iter_vertex1;
            VVertex2D* vertex2 = *iter_vertex2;

            list<VFace2D*> incidentVFaces1;
            vertex1->getIncident3VFaces(incidentVFaces1); //clock-wise

            list<VFace2D*> incidentVFaces2;
            vertex2->getIncident3VFaces(incidentVFaces2); //clock-wise

            list<VFace2D*>::iterator iter_face1 = incidentVFaces1.begin();
            list<VFace2D*>::iterator iter_face2 = incidentVFaces2.end();

            int countOfSameVFace = 0;

            while (iter_face1 != incidentVFaces1.end())
            {
                VFace2D* incVFace1 = *iter_face1;

                iter_face2 = incidentVFaces2.begin();

                while (iter_face2 != incidentVFaces2.end())
                {
                    VFace2D* incVFace2 = *iter_face2;

                    if (linker_gen_to_gen.at(incVFace1->getGenerator()) == incVFace2->getGenerator())
                    {
                        countOfSameVFace++;
                        break;
                    }

                    iter_face2++;
                }

                iter_face1++;
            }

            switch (countOfSameVFace)
            {
            case 3:
            {
                iter_face1 = incidentVFaces1.begin();
                iter_face2 = incidentVFaces2.begin();

                VFace2D* incVFace1 = *iter_face1;
                VFace2D* incVFace2 = NULL;

                //find first same generator
                while (iter_face2 != incidentVFaces2.end())
                {
                    incVFace2 = *iter_face2;

                    if (linker_gen_to_gen.at(incVFace1->getGenerator()) == incVFace2->getGenerator())
                    {
                        break;
                    }

                    iter_face2++;
                }

                // ?ÔøΩÎÇò ?ÔøΩÏñ¥Í∞îÏùÑ ??end?ÔøΩÎ©¥ beginÔø??ÎπÑÍµê
                iter_face1++;
                iter_face2++;

                if (iter_face2 == incidentVFaces2.end())
                {
                    iter_face2 = incidentVFaces2.begin();
                }

                // ÎπÑÍµê
                incVFace1 = *iter_face1;
                incVFace2 = *iter_face2;

                if (linker_gen_to_gen.at(incVFace1->getGenerator()) == incVFace2->getGenerator())
                {
                    doesFindSameVtx = true;
                }

                break;
            }

            default:
            {
                doesFindSameVtx = false;
            }
                break;

            }

            if (doesFindSameVtx)
            {
                break;
            }
            else
            {
                iter_vertex2++;
            }
        }

        if (iter_vertex2 == VVertices2.end())
        {
            return false;
        }
    }

    return true;
}



bool VoronoiDiagram2DC::mergeSecondVFaceToTheFirstIfTheyAreAdjacent( const VFace2D * const VCellToRemain, const VFace2D * const VCellToBeMerged )
{
    list<VEdge2D*> boundaryVEdgesOfVCellToBeMerged;
    VCellToBeMerged->getBoundaryVEdges( boundaryVEdgesOfVCellToBeMerged );

    VEdge2D* VEdgeDefinedByTwoInputVFaces = rg_NULL;
    list<VEdge2D*>::iterator i_bndEdge;
    for ( list<VEdge2D*>::iterator i_bndEdge = boundaryVEdgesOfVCellToBeMerged.begin() ; i_bndEdge != boundaryVEdgesOfVCellToBeMerged.end() ; ++i_bndEdge) {
        VEdge2D* currBoundaryEdge = *i_bndEdge;

        // There is no anomaly case when we merge VFaces of children generators on EdgeGenerator.
        // That means two VFaces can share at most one VEdges.
        // This function is for PolygonVD.
        // If you want to use this function for general cases, you should consider anomaly cases.
        if( currBoundaryEdge->getLeftFace() == VCellToRemain || currBoundaryEdge->getRightFace() == VCellToRemain ) {
            VEdgeDefinedByTwoInputVFaces = currBoundaryEdge;
            break;
        }
    }

    if( VEdgeDefinedByTwoInputVFaces == rg_NULL ) {
        return false;
    }

    VEdge2D* rightHand = VEdgeDefinedByTwoInputVFaces->getRightHand();
    VEdge2D* leftHand  = VEdgeDefinedByTwoInputVFaces->getLeftHand();

    // Connect topology from merged edge (mergedHand) to old Voronoi diagram
    VEdge2D* mergedHand = createEdge(m_VEdges.back()->getID() + 1);
    mergedHand->setLeftFace(const_cast<VFace2D*>(VCellToRemain));

    if( rightHand->getLeftHand() == VEdgeDefinedByTwoInputVFaces ) {
        mergedHand->setStartVertex( rightHand->getStartVertex() );
        mergedHand->setLeftLeg(     rightHand->getLeftLeg() );
        mergedHand->setRightLeg(    rightHand->getRightLeg() );
        mergedHand->setRightFace(   rightHand->getRightFace() );
    }
    else {
        mergedHand->setStartVertex( rightHand->getEndVertex() );
        mergedHand->setLeftLeg(     rightHand->getRightHand() );
        mergedHand->setRightLeg(    rightHand->getLeftHand() );
        mergedHand->setRightFace(   rightHand->getLeftFace() );
    }

    if( leftHand->getLeftLeg() == VEdgeDefinedByTwoInputVFaces ) {
        mergedHand->setEndVertex(   leftHand->getEndVertex() );
        mergedHand->setLeftHand(    leftHand->getLeftHand() );
        mergedHand->setRightHand(   leftHand->getRightHand() );
    }
    else {
        mergedHand->setEndVertex(   leftHand->getStartVertex() );
        mergedHand->setLeftHand(    leftHand->getRightLeg() );
        mergedHand->setRightHand(   leftHand->getLeftLeg() );
    }

    // Connect topology old Voronoi diagram to merged edges (mergedHand)
    mergedHand->getStartVertex()->setFirstVEdge( mergedHand );
    mergedHand->getEndVertex()->setFirstVEdge( mergedHand );
    mergedHand->getLeftFace()->setFirstVEdge( mergedHand );
    mergedHand->getRightFace()->setFirstVEdge( mergedHand );

    if( mergedHand->getLeftLeg()->getLeftHand() == rightHand ) {
        mergedHand->getLeftLeg()->setLeftHand( mergedHand );
    }
    else {
        mergedHand->getLeftLeg()->setRightLeg( mergedHand );
    }

    if( mergedHand->getRightLeg()->getRightHand() == rightHand ) {
        mergedHand->getRightLeg()->setRightHand( mergedHand );
    }
    else {
        mergedHand->getRightLeg()->setLeftLeg( mergedHand );
    }

    if( mergedHand->getLeftHand()->getRightHand() == leftHand ) {
        mergedHand->getLeftHand()->setRightHand( mergedHand );
    }
    else {
        mergedHand->getLeftHand()->setLeftLeg( mergedHand );
    }

    if( mergedHand->getRightHand()->getLeftHand() == leftHand ) {
        mergedHand->getRightHand()->setLeftHand( mergedHand );
    }
    else {
        mergedHand->getRightHand()->setRightLeg( mergedHand );
    }




    VEdge2D* leftLeg   = VEdgeDefinedByTwoInputVFaces->getLeftLeg();
    VEdge2D* rightLeg  = VEdgeDefinedByTwoInputVFaces->getRightLeg();

    // Connect topology from merged edge (mergedLeg) to old Voronoi diagram
    VEdge2D* mergedLeg = createEdge(m_VEdges.back()->getID() + 1);
    mergedLeg->setLeftFace(const_cast<VFace2D*>(VCellToRemain));

    if( leftLeg->getLeftHand() == VEdgeDefinedByTwoInputVFaces ) {
        mergedLeg->setStartVertex( leftLeg->getStartVertex() );
        mergedLeg->setLeftLeg(     leftLeg->getLeftLeg() );
        mergedLeg->setRightLeg(    leftLeg->getRightLeg() );
        mergedLeg->setRightFace(   leftLeg->getRightFace() );
    }
    else {
        mergedLeg->setStartVertex( leftLeg->getEndVertex() );
        mergedLeg->setLeftLeg(     leftLeg->getRightHand() );
        mergedLeg->setRightLeg(    leftLeg->getLeftHand() );
        mergedLeg->setRightFace(   leftLeg->getLeftFace() );
    }

    if( rightLeg->getLeftLeg() == VEdgeDefinedByTwoInputVFaces ) {
        mergedLeg->setEndVertex(   rightLeg->getEndVertex() );
        mergedLeg->setLeftHand(    rightLeg->getLeftHand() );
        mergedLeg->setRightHand(   rightLeg->getRightHand() );
    }
    else {
        mergedLeg->setEndVertex(   rightLeg->getStartVertex() );
        mergedLeg->setLeftHand(    rightLeg->getRightLeg() );
        mergedLeg->setRightHand(   rightLeg->getLeftLeg() );
    }



    // Connect topology old Voronoi diagram to merged edges (mergedLeg)
    mergedLeg->getStartVertex()->setFirstVEdge( mergedLeg );
    mergedLeg->getEndVertex()->setFirstVEdge( mergedLeg );
    mergedLeg->getLeftFace()->setFirstVEdge( mergedLeg );
    mergedLeg->getRightFace()->setFirstVEdge( mergedLeg );

    if( mergedLeg->getLeftLeg()->getLeftHand() == leftLeg ) {
        mergedLeg->getLeftLeg()->setLeftHand( mergedLeg );
    }
    else {
        mergedLeg->getLeftLeg()->setRightLeg( mergedLeg );
    }

    if( mergedLeg->getRightLeg()->getRightHand() == leftLeg ) {
        mergedLeg->getRightLeg()->setRightHand( mergedLeg );
    }
    else {
        mergedLeg->getRightLeg()->setLeftLeg( mergedLeg );
    }

    if( mergedLeg->getLeftHand()->getRightHand() == rightLeg ) {
        mergedLeg->getLeftHand()->setRightHand( mergedLeg );
    }
    else {
        mergedLeg->getLeftHand()->setLeftLeg( mergedLeg );
    }

    if( mergedLeg->getRightHand()->getLeftHand() == rightLeg ) {
        mergedLeg->getRightHand()->setLeftHand( mergedLeg );
    }
    else {
        mergedLeg->getRightHand()->setRightLeg( mergedLeg );
    }



    // color the removed entities to RED
    VEdgeDefinedByTwoInputVFaces->getStartVertex()->setStatus( RED_V );
    VEdgeDefinedByTwoInputVFaces->getEndVertex()->setStatus( RED_V );

    VEdgeDefinedByTwoInputVFaces->setStatus( RED_E );
    rightHand->setStatus( RED_E );
    leftHand->setStatus(  RED_E );
    leftLeg->setStatus(   RED_E );
    rightLeg->setStatus(  RED_E );

    // change pointer of incident VFace of boundary VEdges on VCellToBeMerged
    for( i_bndEdge = boundaryVEdgesOfVCellToBeMerged.begin() ; i_bndEdge != boundaryVEdgesOfVCellToBeMerged.end() ; ++i_bndEdge ) {
        VEdge2D* currBoundaryEdge = *i_bndEdge;

        if( currBoundaryEdge->getStatus() == RED_E ) {
            continue;
        }

        if( currBoundaryEdge->getLeftFace() == VCellToBeMerged ) {
            currBoundaryEdge->setLeftFace( const_cast<VFace2D*>(VCellToRemain) );
        }
        else {
            currBoundaryEdge->setRightFace( const_cast<VFace2D*>(VCellToRemain) );
        }
    }

    Generator2D* currGenerator = VCellToBeMerged->getGenerator();
    m_generators.remove(currGenerator);
    delete currGenerator;

    m_VFaces.remove(const_cast<VFace2D*>(VCellToBeMerged));
    delete VCellToBeMerged;

    // move below function (removeAllExtraneousVVerticesAndVEdges) to outside
    // It is better to call this function after we finish merging process.
    //removeAllExtraneousVVerticesAndVEdges();

    return true;
}


