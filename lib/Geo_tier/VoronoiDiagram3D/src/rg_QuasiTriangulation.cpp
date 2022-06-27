#include "rg_QuasiTriangulation.h"
#include "rg_Point3D.h"
#include "rg_TMatrix3D.h"
#include "rg_QTEdge.h"
#include "rg_QTFace.h"
#include "ConstForVoronoiDiagram3D.h"
using namespace V::GeometryTier;


#include <string>
#include <ctime>
#include <iomanip>
using namespace std;


QuasiTriangulation::QuasiTriangulation()
{
    m_isLinkedVD = rg_FALSE;
}



QuasiTriangulation::~QuasiTriangulation()
{
}



rg_INT QuasiTriangulation::getNumBalls() const
{
    return m_balls.getSize();
}



rg_INT QuasiTriangulation::getNumQTVertices() const
{
    return m_vertices.getSize();
}



rg_INT QuasiTriangulation::getNumQTTetrahedron() const
{
    return m_tetrahedra.getSize();
}



rg_INT QuasiTriangulation::getNumQTGates() const
{
    return m_gates.getSize();
}



rg_INT QuasiTriangulation::getNumFiniteCells() 
{
    QTVertex* virtualVertex = rg_NULL;

    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() ) {
        QTVertex* currVtx = m_vertices.getpEntity();
        if ( currVtx->isVirtual() ) {
            virtualVertex = currVtx;
            break;
        }
    }


    rg_dList<QTTetrahedron*> cellsIncidentToVirtualVertex;
    inquireCellsOfVertex( virtualVertex, cellsIncidentToVirtualVertex );

    rg_INT numFiniteCells = m_tetrahedra.getSize() - cellsIncidentToVirtualVertex.getSize();

    return numFiniteCells;
}



rg_dList<QTVertex>* QuasiTriangulation::getVertices() 
{
    return &m_vertices;
}



rg_dList<QTTetrahedron>* QuasiTriangulation::getTetrahedra() 
{
    return &m_tetrahedra;
}



QTGateList* QuasiTriangulation::getGates() 
{
    return &m_gates;
}



rg_INT QuasiTriangulation::getAbsoluteIDOfVertex(QTVertex* vertex)
{
    rg_INT numVertices = m_vertices.getSize();
    rg_INT vertexID    = vertex->getID();

    while ( vertexID < 0 ) {
        vertexID += numVertices;
    }

    return vertexID;
}



Ball* QuasiTriangulation::addBall(const Ball& ball)
{
    return m_balls.add( ball );
}



QTVertex* QuasiTriangulation::addVertex(const QTVertex& vertex)
{
    return m_vertices.add( vertex );
}



QTTetrahedron* QuasiTriangulation::addTetrahedron(const QTTetrahedron& tetrahedron)
{
    return m_tetrahedra.add( tetrahedron );
}



QTGate* QuasiTriangulation::addGate(const QTGate& gate)
{
    return m_gates.addGate( gate );
}



rg_FLAG QuasiTriangulation::isOnConvexHull(QTVertex* vertex)
{
    QTVertex* infiniteVertex = m_vertices.getFirstpEntity();
    if ( vertex == infiniteVertex )
        return rg_FALSE;

    rg_dList<QTTetrahedron*> incidentTetrahedra;
    inquireIncidentTetrahedraOfVertexInThisWorld(vertex, incidentTetrahedra);

    QTTetrahedron* currTetrahedron = rg_NULL;
    incidentTetrahedra.reset4Loop();
    while ( incidentTetrahedra.setNext4Loop() )  {
        currTetrahedron = incidentTetrahedra.getEntity();
            
        if ( currTetrahedron->isThere( infiniteVertex ) )
            return rg_TRUE;
    }

    return rg_FALSE;
}



rg_FLAG QuasiTriangulation::isOnConvexHull(const QTEdge& edge)
{
    QTVertex* infiniteVertex = m_vertices.getFirstpEntity();

    if ( edge.getStartVertex() == infiniteVertex || edge.getEndVertex() == infiniteVertex )
        return rg_FALSE;

    rg_dList<QTTetrahedron*> incidentTetrahedra;
    inquireIncidentTetrahedraOfEdgeInThisWorld(edge, incidentTetrahedra);

    QTTetrahedron* currTetrahedron = rg_NULL;
    incidentTetrahedra.reset4Loop();
    while ( incidentTetrahedra.setNext4Loop() )  {
        currTetrahedron = incidentTetrahedra.getEntity();
            
        if ( currTetrahedron->isThere( infiniteVertex ) )
            return rg_TRUE;
    }

    return rg_FALSE;
}



rg_FLAG QuasiTriangulation::isOnConvexHull(const QTFace& face)
{
    QTVertex* infiniteVertex = m_vertices.getFirstpEntity();

    if ( face.getMateVertex() == infiniteVertex )
        return rg_TRUE;

    if ( face.getTetrahedron()->isThere(infiniteVertex) )
        return rg_FALSE;

    if ( face.getTwinFace().getMateVertex() == infiniteVertex)
        return rg_TRUE;
    else 
        return rg_FALSE;
}



rg_FLAG QuasiTriangulation::isOnConvexHull(QTTetrahedron* tetrahedron)
{
    QTVertex* infiniteVertex = m_vertices.getFirstpEntity();

    if ( tetrahedron->isThere(infiniteVertex) )
        return rg_TRUE;
    else
        return rg_FALSE;
}



void QuasiTriangulation::resetIDs()
{
    rg_INT vertexID = 1;

    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() ) {
        QTVertex* currVtx = m_vertices.getpEntity();

        if ( !currVtx->isVirtual() ) {
            currVtx->setID( vertexID );
            vertexID++;
        }
        else {
            currVtx->setID( 0 );
        }
    }


    rg_INT cellID = 0;
    m_tetrahedra.reset4Loop();
    while ( m_tetrahedra.setNext4Loop() ) {
        QTTetrahedron* currCell = m_tetrahedra.getpEntity();

        if ( !currCell->isVirtual() ) {
            currCell->setID( cellID );
            cellID++;
        }
    }

    m_tetrahedra.reset4Loop();
    while ( m_tetrahedra.setNext4Loop() ) {
        QTTetrahedron* currCell = m_tetrahedra.getpEntity();

        if ( currCell->isVirtual() ) {
            currCell->setID( cellID );
            cellID++;
        }
    }
}



///////////////////////////////////////////////////////////////////////////////
//   For Queres  
//      T -> V, E, F, T
//
void QuasiTriangulation::inquireBoundingVerticesOfTetrahedronInThisWorld(QTTetrahedron*       givenTetrahedron, 
                                                                         rg_dList<QTVertex*>& vertexList)
{
    if ( givenTetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        for ( rg_INT i=0; i<3; i++ )  
            vertexList.add( givenTetrahedron->getVertex(i) );
    }
    else  {
        for ( rg_INT i=0; i<4; i++ )  
            vertexList.add( givenTetrahedron->getVertex(i) );
    }
}



void QuasiTriangulation::inquireBoundingEdgesOfTetrahedronInThisWorld(QTTetrahedron*    givenTetrahedron, 
                                                                      rg_dList<QTEdge>& edgeList)
{
    QTVertex** vertices = givenTetrahedron->getVertices();
    if ( givenTetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        edgeList.add( QTEdge(givenTetrahedron, vertices[0], vertices[1]) );
        edgeList.add( QTEdge(givenTetrahedron, vertices[1], vertices[2]) );
        edgeList.add( QTEdge(givenTetrahedron, vertices[2], vertices[0]) );
    }
    else  {
        edgeList.add( QTEdge(givenTetrahedron, vertices[0], vertices[1]) );
        edgeList.add( QTEdge(givenTetrahedron, vertices[0], vertices[2]) );
        edgeList.add( QTEdge(givenTetrahedron, vertices[0], vertices[3]) );
        edgeList.add( QTEdge(givenTetrahedron, vertices[1], vertices[2]) );
        edgeList.add( QTEdge(givenTetrahedron, vertices[2], vertices[3]) );
        edgeList.add( QTEdge(givenTetrahedron, vertices[3], vertices[1]) );
    }
}



void QuasiTriangulation::inquireBoundingFacesOfTetrahedronInThisWorld(QTTetrahedron*    givenTetrahedron, 
                                                                      rg_dList<QTFace>& faceList)
{
    QTVertex** vertices = givenTetrahedron->getVertices();
    if ( givenTetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        faceList.add( QTFace(givenTetrahedron, vertices[3]) );        
    }
    else  {
        faceList.add( QTFace(givenTetrahedron, vertices[0]) );
        faceList.add( QTFace(givenTetrahedron, vertices[1]) );
        faceList.add( QTFace(givenTetrahedron, vertices[2]) );
        faceList.add( QTFace(givenTetrahedron, vertices[3]) );
    }
}



void QuasiTriangulation::inquireIncidentTetrahedraOfTetrahedronInThisWorld(QTTetrahedron*            givenTetrahedron, 
                                                                           rg_dList<QTTetrahedron*>& tetrahedronList)
{
    if ( givenTetrahedron->isAnomaly() != DEGENERACY_ANOMALY )  {
        for ( rg_INT i=0; i<4; i++ )  {
            tetrahedronList.addWithoutSame( givenTetrahedron->getNeighbor(i) );
        }
    }
}



///////////////////////////////////////////////////////////////////////////////
//   For Queres  
//      F -> V, E, F, T
//
void QuasiTriangulation::inquireBoundingVerticesOfFaceInThisWorld(const QTFace& givenFace, rg_dList<QTVertex*>& vertexList)
{
    QTTetrahedron* tetrahedron = givenFace.getTetrahedron();
    
    if ( tetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        for ( rg_INT i=0; i<3; i++ )  
            vertexList.add( tetrahedron->getVertex(i) );
    }
    else  {
        QTVertex** vertices   = tetrahedron->getVertices();
        QTVertex*  mateVertex = givenFace.getMateVertex();
        if ( mateVertex == vertices[0] )  {
            vertexList.add(vertices[1]);
            vertexList.add(vertices[2]);
            vertexList.add(vertices[3]);
        }
        else if  ( mateVertex == vertices[1] )  {
            vertexList.add(vertices[0]);
            vertexList.add(vertices[3]);
            vertexList.add(vertices[2]);
        }
        else if  ( mateVertex == vertices[2] )  {
            vertexList.add(vertices[0]);
            vertexList.add(vertices[1]);
            vertexList.add(vertices[3]);
        }
        else  {
            vertexList.add(vertices[0]);
            vertexList.add(vertices[2]);
            vertexList.add(vertices[1]);
        }
    }
}



void QuasiTriangulation::inquireBoundingEdgesOfFaceInThisWorld(const QTFace& givenFace, rg_dList<QTEdge>& edgeList)
{
    rg_dList<QTVertex*> vertexList;
    inquireBoundingVerticesOfFaceInThisWorld(givenFace, vertexList);
    QTVertex**     vertices    = vertexList.getArray();
    QTTetrahedron* tetrahedron = givenFace.getTetrahedron();

    edgeList.add( QTEdge(tetrahedron, vertices[0], vertices[1]) );
    edgeList.add( QTEdge(tetrahedron, vertices[1], vertices[2]) );
    edgeList.add( QTEdge(tetrahedron, vertices[2], vertices[0]) );

    delete [] vertices;
}



void QuasiTriangulation::inquireIncidentFacesOfFaceInThisWorld(const QTFace& givenFace, rg_dList<QTFace>& faceList)
{
    rg_dList<QTEdge> edgeList;
    inquireBoundingEdgesOfFaceInThisWorld(givenFace, edgeList);
    QTEdge currEdge;
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() )  {
        currEdge = edgeList.getEntity();

        rg_dList<QTFace> faceListIncidentTeEdge;
        inquireIncidentFacesOfEdgeInThisWorld(currEdge, faceListIncidentTeEdge);
        faceListIncidentTeEdge.reset4Loop();
        while ( faceListIncidentTeEdge.setNext4Loop() )  {
            QTFace currFace = faceListIncidentTeEdge.getEntity();

            if ( currFace == givenFace )
                continue;
            else
                faceList.add(currFace);            
        }
    }
}



void QuasiTriangulation::inquireIncidentTetrahedraOfFaceInThisWorld(const QTFace& givenFace, rg_dList<QTTetrahedron*>& tetrahedronList)
{
    QTTetrahedron* tetrahedron = givenFace.getTetrahedron();
    tetrahedronList.add( tetrahedron );
    tetrahedronList.add( tetrahedron->getMateTetrahedron(givenFace.getMateVertex()) );
}



///////////////////////////////////////////////////////////////////////////////
//   For Queres  
//       E -> V, E, F, T
//
void QuasiTriangulation::inquireBoundingVerticesOfEdgeInThisWorld(const QTEdge& givenEdge, rg_dList<QTVertex*>& vertexList)
{
    vertexList.add( givenEdge.getStartVertex() );
    vertexList.add( givenEdge.getEndVertex() );
}



void QuasiTriangulation::inquireIncidentEdgesOfEdgeInThisWorld(const QTEdge& givenEdge, rg_dList<QTEdge>& edgeList)
{
    QTTetrahedron* tetrahedron = givenEdge.getTetrahedron();
    QTVertex*      startVertex = givenEdge.getStartVertex();
    QTVertex*      endVertex   = givenEdge.getEndVertex();
    if ( tetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        QTVertex* thirdVertex = rg_NULL;
        for ( rg_INT i=0; i<3; i++ )  {
            QTVertex* currVertex = tetrahedron->getVertex(i);
            if ( currVertex != startVertex  && currVertex != endVertex )  {
                thirdVertex = currVertex;
                break;
            }
        }
        
        edgeList.add( QTEdge(tetrahedron, endVertex, thirdVertex) );
        edgeList.add( QTEdge(tetrahedron, thirdVertex, startVertex) );

        return;
    }

    
    QTTetrahedron* currTetrahedron = tetrahedron;
    do  {
        QTFace currFace(currTetrahedron, startVertex);
        rg_dList<QTVertex*> vertexList;
        inquireBoundingVerticesOfFaceInThisWorld(currFace, vertexList);
        QTVertex** vertices = vertexList.getArray();

        if ( endVertex == vertices[0] )  {
            edgeList.add( QTEdge(currTetrahedron, endVertex, vertices[1]) );
            edgeList.add( QTEdge(currTetrahedron, vertices[1], startVertex) );

            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[2]);
        }
        else if ( endVertex == vertices[1] )  {
            edgeList.add( QTEdge(currTetrahedron, endVertex, vertices[2]) );
            edgeList.add( QTEdge(currTetrahedron, vertices[2], startVertex) );

            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[0]);
        }
        else  {
            edgeList.add( QTEdge(currTetrahedron, endVertex, vertices[0]) );
            edgeList.add( QTEdge(currTetrahedron, vertices[0], startVertex) );
            
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[1]);
        }

        delete [] vertices;
    } while ( tetrahedron != currTetrahedron );
}



void QuasiTriangulation::inquireIncidentFacesOfEdgeInThisWorld(const QTEdge& givenEdge, rg_dList<QTFace>& faceList)
{
    QTTetrahedron* tetrahedron = givenEdge.getTetrahedron();
    if ( tetrahedron->isAnomaly() == DEGENERACY_ANOMALY )  {
        faceList.add( QTFace(tetrahedron, rg_NULL) );        
        return;
    }

    
    QTVertex*  startVertex = givenEdge.getStartVertex();
    QTVertex*  endVertex   = givenEdge.getEndVertex();

    QTTetrahedron* currTetrahedron = tetrahedron;
    do  {
        QTFace currFace(currTetrahedron, startVertex);
        rg_dList<QTVertex*> vertexList;
        inquireBoundingVerticesOfFaceInThisWorld(currFace, vertexList);
        QTVertex** vertices = vertexList.getArray();

        if ( endVertex == vertices[0] )  {
            faceList.add( QTFace(currTetrahedron, vertices[2]) );        
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[2]);
        }
        else if ( endVertex == vertices[1] )  {
            faceList.add( QTFace(currTetrahedron, vertices[0]) );        
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[0]);
        }
        else  {
            faceList.add( QTFace(currTetrahedron, vertices[1]) );        
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[1]);
        }

        delete [] vertices;
    } while ( tetrahedron != currTetrahedron );
}



void QuasiTriangulation::inquireIncidentTetrahedraOfEdgeInThisWorld(const QTEdge& givenEdge, rg_dList<QTTetrahedron*>& tetrahedronList)
{
    QTTetrahedron* tetrahedron = givenEdge.getTetrahedron();
    if ( tetrahedron->isAnomaly() == DEGENERACY_ANOMALY )
        return;

    
    QTTetrahedron* currTetrahedron = tetrahedron;
    QTVertex*  startVertex = givenEdge.getStartVertex();
    QTVertex*  endVertex   = givenEdge.getEndVertex();
    do  {
        tetrahedronList.add(currTetrahedron);

        QTFace currFace(currTetrahedron, startVertex);

        rg_dList<QTVertex*> vertexList;
        inquireBoundingVerticesOfFaceInThisWorld(currFace, vertexList);
        QTVertex** vertices = vertexList.getArray();

        if ( endVertex == vertices[0] )  
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[2]);
        else if ( endVertex == vertices[1] )  
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[0]);
        else  
            currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[1]);

        delete [] vertices;
    } while ( tetrahedron != currTetrahedron );
}



///////////////////////////////////////////////////////////////////////////////
//   For Queres  
//       V -> V, E, F, T
//
void QuasiTriangulation::inquireNeighborVerticesOfVertexInThisWorld(QTVertex* givenVertex, rg_dList<QTVertex*>& vertexList)
{
}



void QuasiTriangulation::inquireIncidentEdgesOfVertexInThisWorld(QTVertex* givenVertex, rg_dList<QTEdge>& edgeList)
{
}



void QuasiTriangulation::inquireIncidentFacesOfVertexInThisWorld(QTVertex* givenVertex, rg_dList<QTFace>& faceList)
{
    rg_dList<QTTetrahedron*> tetrahedronList;
    inquireIncidentTetrahedraOfVertexInThisWorld(givenVertex, tetrahedronList);

    QTTetrahedron* currTetrahedron = rg_NULL;
    tetrahedronList.reset4Loop();
    while ( tetrahedronList.setNext4Loop() )  {
        currTetrahedron = tetrahedronList.getEntity();
    
        QTVertex** vertices = currTetrahedron->getVertices();
        for ( rg_INT i=0; i<4; i++ )  {
            if ( vertices[i] != givenVertex )
                faceList.addWithoutSame( QTFace( currTetrahedron, vertices[i] ) );
        }
    }
}



void QuasiTriangulation::inquireIncidentTetrahedraOfVertexInThisWorld(QTVertex* givenVertex, rg_dList<QTTetrahedron*>& tetrahedronList)
{
    tetrahedronList.add( givenVertex->getFirstTetrahedron() );

    rg_INT i=0;
    rg_dNode<QTTetrahedron*>* startTetrahedronNode = tetrahedronList.getFirstpNode();    
    rg_dNode<QTTetrahedron*>* currTetrahedronNode  = startTetrahedronNode;    
    do  {
        QTTetrahedron** neighbor = currTetrahedronNode->getEntity()->getNeighbors();
        for ( i=0; i<4; i++ ) {
            if ( neighbor[i]->isThere(givenVertex) )
                tetrahedronList.addWithoutSame( neighbor[i] );
        }

        currTetrahedronNode = currTetrahedronNode->getNext();
    } while ( startTetrahedronNode != currTetrahedronNode );
}








Sphere QuasiTriangulation::computeMaxTangentSphereOnQTTetraheron(QTTetrahedron const* ptrTetrahedron)
{
    return ptrTetrahedron->getEmptyTangentSphere();        
}



Sphere QuasiTriangulation::computeMinSphereTangentToTwoBalls(const Sphere& ball1, const Sphere& ball2)
{
    rg_REAL radius = ball1.distanceToSphere(ball2)*0.5;

    rg_Point3D vecV1V2 = ball2.getCenter() - ball1.getCenter();
    vecV1V2.normalize();

    rg_Point3D center = ball1.getCenter() + (ball1.getRadius()+radius)*vecV1V2;

    return Sphere(center, radius);

    /*
    rg_Point3D vecV1V2 = ball2.getCenter() - ball1.getCenter();
    vecV1V2.normalize();

    rg_Point3D ptsOnBall[2];
    ptsOnBall[0] = ball1.getCenter() + (ball1.getRadius()*vecV1V2);
    ptsOnBall[1] = ball2.getCenter() - (ball2.getRadius()*vecV1V2);

    rg_Point3D center = (ptsOnBall[0] + ptsOnBall[1])/2.0;
    rg_REAL    radius;

    rg_REAL distBetCenters    = ball1.getCenter().distance(ball2.getCenter());
    rg_REAL sumRadiiOfTwoBall = ball1.getRadius() + ball2.getRadius();
    if ( distBetCenters > sumRadiiOfTwoBall )
        radius = ptsOnBall[0].distance(ptsOnBall[1])/2.0;
    else
        radius = -ptsOnBall[0].distance(ptsOnBall[1])/2.0;

    return Sphere(center, radius);
    */
}




///////////////////////////////////////////////////////////////////////////
//
//  Quasi-operator
//
//  1. Query primitives for intra-world.
//
//  1.1. Q(v, Y): Query primitives for intra-world.
//
void QuasiTriangulation::inquireVerticesOfVertexInIntraWorld(QTVertex* givenVertex, rg_dList<QTVertex*>& vertexList, QTTetrahedron* world)
{
    QTTetrahedron* tetrahedronOnWorld = rg_NULL;
    if ( world == rg_NULL ) {
        tetrahedronOnWorld = givenVertex->getFirstTetrahedron();
    }
    else {
        tetrahedronOnWorld = world;
    }


    rg_dList<QTTetrahedron*> visitedCellList;
    visitedCellList.add( tetrahedronOnWorld );
    
    int numVisitedCell = 0;
    rg_dNode<QTTetrahedron*>* currCellNode = visitedCellList.getFirstpNode();
    while ( numVisitedCell < visitedCellList.getSize() ) {
        QTTetrahedron* currCell = currCellNode->getEntity();

        rg_dList<QTVertex*> boundingVertices;
        inquireVerticesOfCellInIntraWorld(currCell, boundingVertices);
        
        QTVertex* currVertex = rg_NULL;
        boundingVertices.reset4Loop();
        while ( boundingVertices.setNext4Loop() ) {
            currVertex = boundingVertices.getEntity();

            if ( currVertex == givenVertex ) {
                continue;
            }

            if ( !currVertex->isVisited() ) {
                currVertex->isVisited(rg_TRUE);
                vertexList.add( currVertex );
            }

            QTTetrahedron* mateCell = currCell->getMateTetrahedron( currVertex );
            if ( !mateCell->isVisited() ) {
                visitedCellList.add( mateCell );
            }
        }

        currCell->isVisited( rg_TRUE );
        numVisitedCell++;

        currCellNode = currCellNode->getNext();
    }


    visitedCellList.reset4Loop();
    while ( visitedCellList.setNext4Loop() ) {
        visitedCellList.getEntity()->isVisited( rg_FALSE );
    }

    vertexList.reset4Loop();
    while ( vertexList.setNext4Loop() ) {
        vertexList.getEntity()->isVisited( rg_FALSE );
    }
}



void QuasiTriangulation::inquireEdgesOfVertexInIntraWorld(   QTVertex* givenVertex, rg_dList<QTEdge>& edgeList, QTTetrahedron* world)
{
    QTTetrahedron* tetrahedronOnWorld = rg_NULL;
    if ( world == rg_NULL ) {
        tetrahedronOnWorld = givenVertex->getFirstTetrahedron();
    }
    else {
        tetrahedronOnWorld = world;
    }


    rg_dList<QTVertex*>      visitedVertexList;
    rg_dList<QTTetrahedron*> visitedCellList;
    visitedCellList.add( tetrahedronOnWorld );


    int numVisitedCell = 0;
    rg_dNode<QTTetrahedron*>* currCellNode = visitedCellList.getFirstpNode();
    while ( numVisitedCell < visitedCellList.getSize() ) {
        QTTetrahedron* currCell = currCellNode->getEntity();

        rg_dList<QTVertex*> boundingVertices;
        inquireVerticesOfCellInIntraWorld(currCell, boundingVertices);
        
        QTVertex* currVertex = rg_NULL;
        boundingVertices.reset4Loop();
        while ( boundingVertices.setNext4Loop() ) {
            currVertex = boundingVertices.getEntity();

            if ( currVertex == givenVertex ) {
                continue;
            }

            if ( !currVertex->isVisited() ) {
                currVertex->isVisited(rg_TRUE);
                visitedVertexList.add( currVertex );

                edgeList.add( QTEdge(currCell, givenVertex, currVertex) );
            }

            QTTetrahedron* mateCell = currCell->getMateTetrahedron( currVertex );
            if ( !mateCell->isVisited() ) {
                visitedCellList.add( mateCell );
            }
        }

        currCell->isVisited( rg_TRUE );
        numVisitedCell++;

        currCellNode = currCellNode->getNext();
    }


    visitedCellList.reset4Loop();
    while ( visitedCellList.setNext4Loop() ) {
        visitedCellList.getEntity()->isVisited( rg_FALSE );
    }

    visitedVertexList.reset4Loop();
    while ( visitedVertexList.setNext4Loop() ) {
        visitedVertexList.getEntity()->isVisited( rg_FALSE );
    }
}



void QuasiTriangulation::inquireFacesOfVertexInIntraWorld(   QTVertex* givenVertex, rg_dList<QTFace>& faceList, QTTetrahedron* world)
{
    QTTetrahedron* tetrahedronOnWorld = rg_NULL;
    if ( world == rg_NULL ) {
        tetrahedronOnWorld = givenVertex->getFirstTetrahedron();
    }
    else {
        tetrahedronOnWorld = world;
    }


    rg_dList<QTTetrahedron*> stackOfCellsToBeVisted;
    stackOfCellsToBeVisted.add( tetrahedronOnWorld );

    rg_dList<QTTetrahedron*> visitedCellList;


    while ( stackOfCellsToBeVisted.getSize() > 0 ) {
        QTTetrahedron* currCell = stackOfCellsToBeVisted.popFront();

        if ( currCell->isVisited() ) {
            continue;
        }


        rg_dList<QTVertex*> boundingVertices;
        inquireVerticesOfCellInIntraWorld(currCell, boundingVertices);
        
        QTVertex* currVertex = rg_NULL;
        boundingVertices.reset4Loop();
        while ( boundingVertices.setNext4Loop() ) {
            currVertex = boundingVertices.getEntity();

            if ( currVertex == givenVertex ) {
                continue;
            }

            QTTetrahedron* mateCell = currCell->getMateTetrahedron( currVertex );
            if ( !mateCell->isVisited() ) {
                faceList.add( QTFace(currCell, currVertex) );
                stackOfCellsToBeVisted.add( mateCell );
            }
        }

        currCell->isVisited( rg_TRUE );
        visitedCellList.add( currCell );
    }


    visitedCellList.reset4Loop();
    while ( visitedCellList.setNext4Loop() ) {
        visitedCellList.getEntity()->isVisited( rg_FALSE );
    }
}



void QuasiTriangulation::inquireCellsOfVertexInIntraWorld(   QTVertex* givenVertex, rg_dList<QTTetrahedron*>& tetrahedronList, QTTetrahedron* world)
{
    QTTetrahedron* tetrahedronOnWorld = rg_NULL;
    if ( world == rg_NULL ) {
        tetrahedronOnWorld = givenVertex->getFirstTetrahedron();
    }
    else {
        tetrahedronOnWorld = world;
    }


    rg_dList<QTTetrahedron*> stackOfCellsToBeVisted;
    tetrahedronOnWorld->isVisited( rg_TRUE );
    stackOfCellsToBeVisted.add( tetrahedronOnWorld );

    while ( stackOfCellsToBeVisted.getSize() > 0 ) {
        QTTetrahedron* currCell = stackOfCellsToBeVisted.popFront();
        tetrahedronList.add( currCell );

        rg_dList<QTVertex*> boundingVertices;
        inquireVerticesOfCellInIntraWorld(currCell, boundingVertices);
        
        QTVertex* currVertex = rg_NULL;
        boundingVertices.reset4Loop();
        while ( boundingVertices.setNext4Loop() ) {
            currVertex = boundingVertices.getEntity();

            if ( currVertex == givenVertex ) {
                continue;
            }

            QTTetrahedron* mateCell = currCell->getMateTetrahedron( currVertex );
            if ( !mateCell->isVisited() ) {
                mateCell->isVisited( rg_TRUE );
                stackOfCellsToBeVisted.add( mateCell );
            }
        }
    }


    tetrahedronList.reset4Loop();
    while ( tetrahedronList.setNext4Loop() ) {
        tetrahedronList.getEntity()->isVisited( rg_FALSE );
    }
}




//  1.2. Q(e, Y): Query primitives for intra-world.
//
void QuasiTriangulation::inquireVerticesOfEdgeInIntraWorld(const QTEdge& givenEdge, rg_dList<QTVertex*>& vertexList, QTTetrahedron* world)
{
    vertexList.add( givenEdge.getStartVertex() );
    vertexList.add( givenEdge.getEndVertex() );
}



void QuasiTriangulation::inquireEdgesOfEdgeInIntraWorld(   const QTEdge& givenEdge, rg_dList<QTEdge>& edgeList, QTTetrahedron* world)
{
    QTTetrahedron* tetrahedronOnWorld = rg_NULL;
    if ( world == rg_NULL ) {
        tetrahedronOnWorld = givenEdge.getTetrahedron();
    }
    else {
        tetrahedronOnWorld = world;
    }


    QTVertex*  startVertex = givenEdge.getStartVertex();
    QTVertex*  endVertex   = givenEdge.getEndVertex();

    if ( tetrahedronOnWorld->isAnomaly() == DEGENERACY_ANOMALY )  {
        QTFace currFace(tetrahedronOnWorld, startVertex);
        inquireEdgesOfFaceInIntraWorld(currFace, edgeList);
    }
    else {   
        QTTetrahedron* currTetrahedron = tetrahedronOnWorld;
        do  {
            QTFace currFace(currTetrahedron, startVertex);
            rg_dList<QTVertex*> vertexList;
            inquireVerticesOfFaceInIntraWorld(currFace, vertexList);
            QTVertex** vertices = vertexList.getArray();

            if ( endVertex == vertices[0] )  {
                edgeList.add( QTEdge(currTetrahedron, endVertex, vertices[1]) );        
                edgeList.add( QTEdge(currTetrahedron, vertices[1], startVertex) );        

                currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[2]);
            }
            else if ( endVertex == vertices[1] )  {
                edgeList.add( QTEdge(currTetrahedron, endVertex, vertices[2]) );        
                edgeList.add( QTEdge(currTetrahedron, vertices[2], startVertex) );      
                
                currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[0]);
            }
            else  if ( endVertex == vertices[2] ) {
                edgeList.add( QTEdge(currTetrahedron, endVertex, vertices[0]) );        
                edgeList.add( QTEdge(currTetrahedron, vertices[0], startVertex) );       
                
                currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[1]);
            }

            delete [] vertices;
        } while ( tetrahedronOnWorld != currTetrahedron );
    }
}



void QuasiTriangulation::inquireFacesOfEdgeInIntraWorld(   const QTEdge& givenEdge, rg_dList<QTFace>& faceList, QTTetrahedron* world)
{
    QTTetrahedron* tetrahedronOnWorld = rg_NULL;
    if ( world == rg_NULL ) {
        tetrahedronOnWorld = givenEdge.getTetrahedron();
    }
    else {
        tetrahedronOnWorld = world;
    }


    if ( tetrahedronOnWorld->isAnomaly() == DEGENERACY_ANOMALY )  {
        faceList.add( QTFace(tetrahedronOnWorld, rg_NULL) );        
    }
    else {
    
        QTVertex*  startVertex = givenEdge.getStartVertex();
        QTVertex*  endVertex   = givenEdge.getEndVertex();

        QTTetrahedron* currTetrahedron = tetrahedronOnWorld;
        do  {
            QTFace currFace(currTetrahedron, startVertex);
            rg_dList<QTVertex*> vertexList;
            inquireVerticesOfFaceInIntraWorld(currFace, vertexList);
            QTVertex** vertices = vertexList.getArray();

            if ( endVertex == vertices[0] )  {
                faceList.add( QTFace(currTetrahedron, vertices[2]) );        
                currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[2]);
            }
            else if ( endVertex == vertices[1] )  {
                faceList.add( QTFace(currTetrahedron, vertices[0]) );        
                currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[0]);
            }
            else  if ( endVertex == vertices[2] ) {
                faceList.add( QTFace(currTetrahedron, vertices[1]) );        
                currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[1]);
            }

            delete [] vertices;
        } while ( tetrahedronOnWorld != currTetrahedron );
    }
}



void QuasiTriangulation::inquireCellsOfEdgeInIntraWorld(   const QTEdge& givenEdge, rg_dList<QTTetrahedron*>& tetrahedronList, QTTetrahedron* world)
{
    QTTetrahedron* tetrahedronOnWorld = rg_NULL;
    if ( world == rg_NULL ) {
        tetrahedronOnWorld = givenEdge.getTetrahedron();
    }
    else {
        tetrahedronOnWorld = world;
    }


    if ( tetrahedronOnWorld->isAnomaly() != DEGENERACY_ANOMALY )  {
        QTVertex*  startVertex = givenEdge.getStartVertex();
        QTVertex*  endVertex   = givenEdge.getEndVertex();

        QTTetrahedron* currTetrahedron = tetrahedronOnWorld;
        do  {
            tetrahedronList.add( currTetrahedron );

            QTFace currFace(currTetrahedron, startVertex);
            rg_dList<QTVertex*> vertexList;
            inquireVerticesOfFaceInIntraWorld(currFace, vertexList);
            QTVertex** vertices = vertexList.getArray();

            if ( endVertex == vertices[0] )  {
                currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[2]);
            }
            else if ( endVertex == vertices[1] )  {
                currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[0]);
            }
            else if ( endVertex == vertices[2] ) {
                currTetrahedron = currTetrahedron->getMateTetrahedron(vertices[1]);
            }

            delete [] vertices;
        } while ( tetrahedronOnWorld != currTetrahedron );
    }
}




//  1.3. Q(f, Y): Query primitives for intra-world.
//
void QuasiTriangulation::inquireVerticesOfFaceInIntraWorld(const QTFace& givenFace, rg_dList<QTVertex*>& vertexList, QTTetrahedron* world)
{
    QTTetrahedron* tetrahedron = givenFace.getTetrahedron();
    
    if ( tetrahedron->isAnomaly() != DEGENERACY_ANOMALY )  {
        QTVertex** vertices   = tetrahedron->getVertices();
        QTVertex*  mateVertex = givenFace.getMateVertex();
        if ( mateVertex == vertices[0] )  {
            vertexList.add(vertices[1]);
            vertexList.add(vertices[2]);
            vertexList.add(vertices[3]);
        }
        else if  ( mateVertex == vertices[1] )  {
            vertexList.add(vertices[0]);
            vertexList.add(vertices[3]);
            vertexList.add(vertices[2]);
        }
        else if  ( mateVertex == vertices[2] )  {
            vertexList.add(vertices[0]);
            vertexList.add(vertices[1]);
            vertexList.add(vertices[3]);
        }
        else  {
            vertexList.add(vertices[0]);
            vertexList.add(vertices[2]);
            vertexList.add(vertices[1]);
        }
    }
    else  {
        for ( rg_INT i=0; i<3; i++ )  {
            vertexList.add( tetrahedron->getVertex(i) );
        }
    }
}



void QuasiTriangulation::inquireEdgesOfFaceInIntraWorld(   const QTFace& givenFace, rg_dList<QTEdge>& edgeList, QTTetrahedron* world)
{
    rg_dList<QTVertex*> vertexList;
    inquireVerticesOfFaceInIntraWorld(givenFace, vertexList);
    QTVertex**     vertices    = vertexList.getArray();
    QTTetrahedron* tetrahedron = givenFace.getTetrahedron();

    edgeList.add( QTEdge(tetrahedron, vertices[0], vertices[1]) );
    edgeList.add( QTEdge(tetrahedron, vertices[1], vertices[2]) );
    edgeList.add( QTEdge(tetrahedron, vertices[2], vertices[0]) );

    delete [] vertices;
}



void QuasiTriangulation::inquireFacesOfFaceInIntraWorld(   const QTFace& givenFace, rg_dList<QTFace>& faceList, QTTetrahedron* world)
{
    const int numCellsInSolution = 2;
    QTTetrahedron* solutionSpace[numCellsInSolution] = { givenFace.getTetrahedron(), 
                                                         givenFace.getIncidentTetrahedron() };

    QTVertex*    vertexOnFace = givenFace.getMateVertex();

    for ( int i_cell=0; i_cell<numCellsInSolution; i_cell++ ) {
        for ( int i_vertex=0; i_vertex<4; i_vertex++ ) {
            QTVertex* vertexOnCell = solutionSpace[i_cell]->getVertex(i_vertex);
            if ( vertexOnFace != vertexOnCell ) {
                faceList.add( QTFace(solutionSpace[i_cell], vertexOnCell) );
            }
        }
    }
}



void QuasiTriangulation::inquireCellsOfFaceInIntraWorld(   const QTFace& givenFace, rg_dList<QTTetrahedron*>& tetrahedronList, QTTetrahedron* world)
{
    tetrahedronList.add( givenFace.getTetrahedron() );
    tetrahedronList.add( givenFace.getIncidentTetrahedron() );
}




//  1.4. Q(c, Y): Query primitives for intra-world.
//
void QuasiTriangulation::inquireVerticesOfCellInIntraWorld(QTTetrahedron* givenCell, rg_dList<QTVertex*>& vertexList, QTTetrahedron* world)
{
    QTVertex** vertices = givenCell->getVertices();

    if ( givenCell->isAnomaly() != DEGENERACY_ANOMALY )  {
        vertexList.add( vertices[0] );
        vertexList.add( vertices[1] );
        vertexList.add( vertices[2] );
        vertexList.add( vertices[3] );
    }
    else  {
        vertexList.add( vertices[0] );
        vertexList.add( vertices[1] );
        vertexList.add( vertices[2] );
    }
}



void QuasiTriangulation::inquireEdgesOfCellInIntraWorld(   QTTetrahedron* givenCell, rg_dList<QTEdge>& edgeList, QTTetrahedron* world)
{
    QTVertex** vertices = givenCell->getVertices();

    if ( givenCell->isAnomaly() != DEGENERACY_ANOMALY )  {
        edgeList.add( QTEdge(givenCell, vertices[0], vertices[1]) );
        edgeList.add( QTEdge(givenCell, vertices[0], vertices[2]) );
        edgeList.add( QTEdge(givenCell, vertices[0], vertices[3]) );
        edgeList.add( QTEdge(givenCell, vertices[1], vertices[2]) );
        edgeList.add( QTEdge(givenCell, vertices[2], vertices[3]) );
        edgeList.add( QTEdge(givenCell, vertices[3], vertices[1]) );
    }
    else  {
        edgeList.add( QTEdge(givenCell, vertices[0], vertices[1]) );
        edgeList.add( QTEdge(givenCell, vertices[1], vertices[2]) );
        edgeList.add( QTEdge(givenCell, vertices[2], vertices[0]) );
    }
}



void QuasiTriangulation::inquireFacesOfCellInIntraWorld(   QTTetrahedron* givenCell, rg_dList<QTFace>& faceList, QTTetrahedron* world)
{
    QTVertex** vertices = givenCell->getVertices();

    if ( givenCell->isAnomaly() != DEGENERACY_ANOMALY )  {
        faceList.add( QTFace(givenCell, vertices[0]) );
        faceList.add( QTFace(givenCell, vertices[1]) );
        faceList.add( QTFace(givenCell, vertices[2]) );
        faceList.add( QTFace(givenCell, vertices[3]) );
    }
    else  {
        faceList.add( QTFace(givenCell, vertices[3]) );        
    }
}



void QuasiTriangulation::inquireCellsOfCellInIntraWorld(   QTTetrahedron* givenCell, rg_dList<QTTetrahedron*>& tetrahedronList, QTTetrahedron* world)
{
    if ( givenCell->isAnomaly() != DEGENERACY_ANOMALY )  {
        for ( rg_INT i=0; i<4; i++ )  {
            tetrahedronList.add( givenCell->getNeighbor(i) );
        }
    }
}




//  2. Query primitives for inter-world.
//
void QuasiTriangulation::inquireSmallWorldsAtVertex(QTVertex* givenVertex, rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world)
{
    rg_dList<QTTetrahedron*> cellList;
    inquireCellsOfVertexInIntraWorld(givenVertex, cellList, world);
    cellList.reset4Loop();
    while ( cellList.setNext4Loop() ) {
        cellList.getEntity()->isVisited(rg_TRUE);
    }


    rg_dList<QTEdge> edgeList;
    inquireEdgesOfVertexInIntraWorld(givenVertex, edgeList, world);

    QTEdge* currEdge = rg_NULL;
    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        currEdge = edgeList.getpEntity();
    
        QTGate* currGate = m_gates.find( QTGate(currEdge->getStartVertex(), currEdge->getEndVertex()) );
        if ( currGate != rg_NULL ) {
            if ( currGate->getBigTetrahedron()->isVisited() ) {
                rg_dList<QTTetrahedron*>* smallWorldInEdge = currGate->getSmallTetrahedra();
                smallWorldList.append( *smallWorldInEdge );
            }
        }
    }


    cellList.reset4Loop();
    while ( cellList.setNext4Loop() ) {
        cellList.getEntity()->isVisited(rg_FALSE);
    }
}



void QuasiTriangulation::inquireSmallWorldsAtEdge(  const QTEdge&  givenEdge,   rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world)
{
    QTGate* gate = m_gates.find( QTGate(givenEdge.getStartVertex(), givenEdge.getEndVertex()) );
    if ( gate != rg_NULL ) {
        rg_dList<QTTetrahedron*> cellList;
        inquireCellsOfEdgeInIntraWorld( givenEdge, cellList, world);

        cellList.reset4Loop();
        while ( cellList.setNext4Loop() ) {
            cellList.getEntity()->isVisited(rg_TRUE);
        }

        if ( gate->getBigTetrahedron()->isVisited() ) {
            rg_dList<QTTetrahedron*>* smallWorldInEdge = gate->getSmallTetrahedra();
            smallWorldList.append( *smallWorldInEdge );
        }

        cellList.reset4Loop();
        while ( cellList.setNext4Loop() ) {
            cellList.getEntity()->isVisited(rg_FALSE);
        }
    }
}



void QuasiTriangulation::inquireSmallWorldsAtFace(  const QTFace&  givenFace,   rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world)
{
    rg_dList<QTEdge> edgeList;
    inquireEdgesOfFaceInIntraWorld( givenFace, edgeList, world);

    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        QTEdge currEdge = edgeList.getEntity();

        inquireSmallWorldsAtEdge( currEdge, smallWorldList, world);
    }
}



void QuasiTriangulation::inquireSmallWorldsAtCell(  QTTetrahedron* givenCell,   rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world)
{
    rg_dList<QTEdge> edgeList;
    inquireEdgesOfCellInIntraWorld( givenCell, edgeList, world);

    edgeList.reset4Loop();
    while ( edgeList.setNext4Loop() ) {
        QTEdge currEdge = edgeList.getEntity();

        inquireSmallWorldsAtEdge( currEdge, smallWorldList, world);
    }
}




void QuasiTriangulation::inquireWholeSmallWorldsAtVertex(QTVertex* givenVertex, rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world)
{
    QTTetrahedron* givenWorld = rg_NULL;
    if ( world == rg_NULL ) {
        givenWorld = givenVertex->getFirstTetrahedron();
    }
    else {
        givenWorld = world;
    }

    rg_dList<QTTetrahedron*> worldList;
    worldList.add( givenWorld );

    while ( worldList.getSize() > 0 ) {
        QTTetrahedron* currWorld = worldList.popFront();

        rg_dList<QTTetrahedron*> cellList;
        inquireCellsOfVertexInIntraWorld(givenVertex, cellList, currWorld);

        cellList.reset4Loop();
        while ( cellList.setNext4Loop() ) {
            cellList.getEntity()->isVisited(rg_TRUE);
        }


        rg_dList<QTEdge> edgeList;
        inquireEdgesOfVertexInIntraWorld(givenVertex, edgeList, currWorld);

        QTEdge* currEdge = rg_NULL;
        edgeList.reset4Loop();
        while ( edgeList.setNext4Loop() ) {
            currEdge = edgeList.getpEntity();
    
            QTGate* currGate = m_gates.find( QTGate(currEdge->getStartVertex(), currEdge->getEndVertex()) );
            if ( currGate != rg_NULL ) {
                if ( currGate->getBigTetrahedron()->isVisited() ) {
                    rg_dList<QTTetrahedron*>* smallWorldInEdge = currGate->getSmallTetrahedra();
                    smallWorldList.append( *smallWorldInEdge );
                    worldList.append( *smallWorldInEdge );
                }
            }
        }


        cellList.reset4Loop();
        while ( cellList.setNext4Loop() ) {
            cellList.getEntity()->isVisited(rg_FALSE);
        }
    }
}




//  3. Query primitives for whole-world.
//
void QuasiTriangulation::inquireVerticesOfVertex(QTVertex* givenVertex, rg_dList<QTVertex*>& vertexList           )
{
    rg_dList<QTTetrahedron*> contributingWorlds;
    contributingWorlds.add( givenVertex->getFirstTetrahedron() );

    inquireWholeSmallWorldsAtVertex(givenVertex, contributingWorlds, givenVertex->getFirstTetrahedron() );

    while ( contributingWorlds.getSize() > 0 ) {
        QTTetrahedron* currWorld = contributingWorlds.popFront();

        inquireVerticesOfVertexInIntraWorld(givenVertex, vertexList, currWorld);
    }
}



void QuasiTriangulation::inquireEdgesOfVertex(   QTVertex* givenVertex, rg_dList<QTEdge>& edgeList               )
{
    rg_dList<QTTetrahedron*> contributingWorlds;
    contributingWorlds.add( givenVertex->getFirstTetrahedron() );

    inquireWholeSmallWorldsAtVertex(givenVertex, contributingWorlds, givenVertex->getFirstTetrahedron() );

    while ( contributingWorlds.getSize() > 0 ) {
        QTTetrahedron* currWorld = contributingWorlds.popFront();

        inquireEdgesOfVertexInIntraWorld(givenVertex, edgeList, currWorld);
    }
}



void QuasiTriangulation::inquireFacesOfVertex(   QTVertex* givenVertex, rg_dList<QTFace>& faceList               )
{
    rg_dList<QTTetrahedron*> contributingWorlds;
    contributingWorlds.add( givenVertex->getFirstTetrahedron() );

    inquireWholeSmallWorldsAtVertex(givenVertex, contributingWorlds, givenVertex->getFirstTetrahedron() );

    while ( contributingWorlds.getSize() > 0 ) {
        QTTetrahedron* currWorld = contributingWorlds.popFront();

        inquireFacesOfVertexInIntraWorld(givenVertex, faceList, currWorld);
    }
}



void QuasiTriangulation::inquireCellsOfVertex(   QTVertex* givenVertex, rg_dList<QTTetrahedron*>& tetrahedronList ) 
{
    rg_dList<QTTetrahedron*> contributingWorlds;
    contributingWorlds.add( givenVertex->getFirstTetrahedron() );

    inquireWholeSmallWorldsAtVertex(givenVertex, contributingWorlds, givenVertex->getFirstTetrahedron() );

    while ( contributingWorlds.getSize() > 0 ) {
        QTTetrahedron* currWorld = contributingWorlds.popFront();

        inquireCellsOfVertexInIntraWorld(givenVertex, tetrahedronList, currWorld);
    }
}




void QuasiTriangulation::inquireVerticesOfEdge(const QTEdge& givenEdge, rg_dList<QTVertex*>& vertexList           )
{
    rg_dList<QTTetrahedron*> contributingWorlds;
    contributingWorlds.add( givenEdge.getTetrahedron() );

    inquireSmallWorldsAtEdge(givenEdge, contributingWorlds, givenEdge.getTetrahedron() );

    while ( contributingWorlds.getSize() > 0 ) {
        QTTetrahedron* currWorld = contributingWorlds.popFront();

        inquireVerticesOfEdgeInIntraWorld(givenEdge, vertexList, currWorld);
    }
}



void QuasiTriangulation::inquireEdgesOfEdge(   const QTEdge& givenEdge, rg_dList<QTEdge>& edgeList               )
{
    rg_dList<QTTetrahedron*> contributingWorlds;
    contributingWorlds.add( givenEdge.getTetrahedron() );

    inquireSmallWorldsAtEdge(givenEdge, contributingWorlds, givenEdge.getTetrahedron() );

    while ( contributingWorlds.getSize() > 0 ) {
        QTTetrahedron* currWorld = contributingWorlds.popFront();

        inquireEdgesOfEdgeInIntraWorld(givenEdge, edgeList, currWorld);
    }
}



void QuasiTriangulation::inquireFacesOfEdge(   const QTEdge& givenEdge, rg_dList<QTFace>& faceList               )
{
    rg_dList<QTTetrahedron*> contributingWorlds;
    contributingWorlds.add( givenEdge.getTetrahedron() );

    inquireSmallWorldsAtEdge(givenEdge, contributingWorlds, givenEdge.getTetrahedron() );

    while ( contributingWorlds.getSize() > 0 ) {
        QTTetrahedron* currWorld = contributingWorlds.popFront();

        inquireFacesOfEdgeInIntraWorld(givenEdge, faceList, currWorld);
    }
}



void QuasiTriangulation::inquireCellsOfEdge(   const QTEdge& givenEdge, rg_dList<QTTetrahedron*>& tetrahedronList )
{
    rg_dList<QTTetrahedron*> contributingWorlds;
    contributingWorlds.add( givenEdge.getTetrahedron() );

    inquireSmallWorldsAtEdge(givenEdge, contributingWorlds, givenEdge.getTetrahedron() );

    while ( contributingWorlds.getSize() > 0 ) {
        QTTetrahedron* currWorld = contributingWorlds.popFront();

        inquireCellsOfEdgeInIntraWorld(givenEdge, tetrahedronList, currWorld);
    }
}




void QuasiTriangulation::inquireVerticesOfFace(const QTFace& givenFace, rg_dList<QTVertex*>& vertexList           )
{
    inquireVerticesOfFaceInIntraWorld(givenFace, vertexList);
}



void QuasiTriangulation::inquireEdgesOfFace(   const QTFace& givenFace, rg_dList<QTEdge>& edgeList               )
{
    inquireEdgesOfFaceInIntraWorld(givenFace, edgeList);
}



void QuasiTriangulation::inquireFacesOfFace(   const QTFace& givenFace, rg_dList<QTFace>& faceList               )
{
    inquireFacesOfFaceInIntraWorld(givenFace, faceList);
}



void QuasiTriangulation::inquireCellsOfFace(   const QTFace& givenFace, rg_dList<QTTetrahedron*>& tetrahedronList )
{
    inquireCellsOfFaceInIntraWorld(givenFace, tetrahedronList);
}




void QuasiTriangulation::inquireVerticesOfCell(QTTetrahedron* givenCell, rg_dList<QTVertex*>& vertexList           )
{
    inquireVerticesOfCellInIntraWorld(givenCell, vertexList);
}



void QuasiTriangulation::inquireEdgesOfCell(   QTTetrahedron* givenCell, rg_dList<QTEdge>& edgeList               )
{
    inquireEdgesOfCellInIntraWorld(givenCell, edgeList);
}



void QuasiTriangulation::inquireFacesOfCell(   QTTetrahedron* givenCell, rg_dList<QTFace>& faceList               )
{
    inquireFacesOfCellInIntraWorld(givenCell, faceList);
}



void QuasiTriangulation::inquireCellsOfCell(   QTTetrahedron* givenCell, rg_dList<QTTetrahedron*>& tetrahedronList )
{
    inquireCellsOfCellInIntraWorld(givenCell, tetrahedronList);
}




//
//  END of Quasi-operator
//  
///////////////////////////////////////////////////////////////////////////




rg_INT QuasiTriangulation::computeMinSphereTangentToThreeBalls(const Sphere& ball1, 
                                                               const Sphere& ball2, 
                                                               const Sphere& ball3, 
                                                               Sphere* tangentSphere)
{
    return makeCircumcircleOnCenterPlane(ball1, ball2, ball3, tangentSphere[0], tangentSphere[1]);
}



rg_INT QuasiTriangulation::makeCircumcircleOnCenterPlane(const Sphere& s1,
														 const Sphere& s2,
														 const Sphere& s3,
														 Sphere& result1,
														 Sphere& result2)
{
	rg_Point3D vec21 = s2.getCenter()-s1.getCenter();
	rg_Point3D vec31 = s3.getCenter()-s1.getCenter();
	rg_Point3D planeNormal = vec21.crossProduct(vec31);


    rg_Vector3D zeroVec(0.0, 0.0, 0.0);
    if ( planeNormal == zeroVec )  {
        planeNormal = vec21;
        planeNormal.normalize();

	    rg_TMatrix3D transMat;   //translate and then rotate
        if( !rg_EQ( planeNormal.innerProduct( rg_Point3D(0.0, 0.0, 1.0) ), -1) )
	        transMat.rotate(planeNormal, rg_Point3D(0,0,1)); //rotate from planeNormal to z-axis
        else
            transMat.rotateY(rg_PI);

        transMat.rotateY(rg_PI*0.5);
        planeNormal = transMat*planeNormal;
    }


	rg_TMatrix3D transformMat;   //translate and then rotate
	transformMat.translate( -(s1.getCenter()) ); //translate 
	transformMat.rotate(planeNormal, rg_Point3D(0,0,1)); //rotate from planeNormal to z-axis
	
	rg_TMatrix3D invTransformMat;
	invTransformMat.rotate(rg_Point3D(0,0,1), planeNormal);
	invTransformMat.translate( s1.getCenter() );

    
	rg_Circle2D transformedCircle1((transformMat*s1.getCenter()).getX(),
								   (transformMat*s1.getCenter()).getY(),
								   s1.getRadius() );
	rg_Circle2D transformedCircle2((transformMat*s2.getCenter()).getX(),
								   (transformMat*s2.getCenter()).getY(),
								   s2.getRadius() );
	rg_Circle2D transformedCircle3((transformMat*s3.getCenter()).getX(),
								   (transformMat*s3.getCenter()).getY(),
								   s3.getRadius() );
	rg_INT numOfCC = 0;
	rg_Circle2D tangentCircle1, tangentCircle2;
	numOfCC = makeCircumcircle(transformedCircle1,transformedCircle2,transformedCircle3,
							   tangentCircle1, tangentCircle2);

	result1.setCenter(invTransformMat*rg_Point3D(tangentCircle1.getCenterPt().getX(),tangentCircle1.getCenterPt().getY(),0));
	result2.setCenter(invTransformMat*rg_Point3D(tangentCircle2.getCenterPt().getX(),tangentCircle2.getCenterPt().getY(),0));
	result1.setRadius(tangentCircle1.getRadius());
	result2.setRadius(tangentCircle2.getRadius());

	return numOfCC;
}



rg_INT QuasiTriangulation::makeCircumcircle(const rg_Circle2D& circle1,
									        const rg_Circle2D& circle2,
									        const rg_Circle2D& circle3,
								            rg_Circle2D&	result1, 
								            rg_Circle2D&	result2)
{
	rg_INT numOfCC = 0;

	rg_ImplicitEquation tangentLine1(1);
	rg_ImplicitEquation tangentLine2(1);

	rg_Point2D z3;
	z3 = makeTangentLinesInWPlane(circle1, circle2, circle3, tangentLine1, tangentLine2);

	rg_REAL smallest = 0;
	if( z3 == circle1.getCenterPt() )
	{
		smallest = circle1.getRadius();
	}
	else if ( z3 == circle2.getCenterPt() )
	{
		smallest = circle2.getRadius();
	}
	else
	{
		smallest = circle3.getRadius();
	}

	//하나인지 두개인지 판단 (point location problem)
	rg_REAL sign1 = tangentLine1.evaluateImpEquation(0, 0);
	rg_REAL sign2 = tangentLine2.evaluateImpEquation(0, 0);

	if( sign1 * sign2 < 0 ) //alpha region
	{
		numOfCC = 1;
		if( sign1 > 0 )
		{
			result1 = transformW2Z(tangentLine1, z3);		//+ -
			result1.setRadius( result1.getRadius() - smallest );
		}
		else
		{
			result1 = transformW2Z(tangentLine2, z3);		//- +
			result1.setRadius( result1.getRadius() - smallest );
		}
	}
	else if( sign1 >0 && sign2 >0 ) //gamma region
	{
		numOfCC = 2;
		result1 = transformW2Z(tangentLine1, z3);			//+ +
		result1.setRadius( result1.getRadius() - smallest );

		result2 = transformW2Z(tangentLine2, z3);
		result2.setRadius( result2.getRadius() - smallest );
	}
	else if( sign1 < 0 && sign2 < 0 )						//- -
	{
		numOfCC = 0;
	}
	else if( rg_EQ(sign1, 0., resNeg15) && sign2 > 0 )						//0 +
	{
		numOfCC = 1;
		result1 = transformW2Z(tangentLine2, z3);
		result1.setRadius( result1.getRadius() - smallest );
	}
	else if( sign1 > 0 && rg_EQ(sign2, 0., resNeg15) )						//+ 0
	{
		numOfCC = 1;
		result1 = transformW2Z(tangentLine1, z3);
		result1.setRadius( result1.getRadius() - smallest );
	}
	else												//0 0
	{													//0 -
		numOfCC = 0;									//- 0
	}

	return numOfCC;
}

//this member function shrink radii of circles
//and c3Tilde becomes a point
void QuasiTriangulation::shrinkCircle(const		rg_Circle2D& c1,
                                      const		rg_Circle2D& c2,
							          const		rg_Circle2D& c3,
							          rg_Circle2D&	c1Tilde,
							          rg_Circle2D&	c2Tilde,
							          rg_Point2D&		c3Tilde)
{
	if( c1.getRadius() < c2.getRadius() )
	{
		if( c1.getRadius() < c3.getRadius() )
		{
			c1Tilde.setCircle(c2.getCenterPt(), c2.getRadius() - c1.getRadius());
			c2Tilde.setCircle(c3.getCenterPt(), c3.getRadius() - c1.getRadius());
			c3Tilde = c1.getCenterPt();
		}
		else
		{
			c1Tilde.setCircle(c1.getCenterPt(), c1.getRadius() - c3.getRadius());
			c2Tilde.setCircle(c2.getCenterPt(), c2.getRadius() - c3.getRadius());
			c3Tilde = c3.getCenterPt();
		}
	}
	else
	{
		if( c2.getRadius() < c3.getRadius() )
		{
			c1Tilde.setCircle(c1.getCenterPt(), c1.getRadius() - c2.getRadius());
			c2Tilde.setCircle(c3.getCenterPt(), c3.getRadius() - c2.getRadius());
			c3Tilde = c2.getCenterPt();
		}
		else
		{
			c1Tilde.setCircle(c1.getCenterPt(), c1.getRadius() - c3.getRadius());
			c2Tilde.setCircle(c2.getCenterPt(), c2.getRadius() - c3.getRadius());
			c3Tilde = c3.getCenterPt();
		}
	}
}

rg_Point2D QuasiTriangulation::makeTangentLinesInWPlane(const rg_Circle2D&		c1, 
										                const rg_Circle2D&		c2,
											            const rg_Circle2D&		c3,
											            rg_ImplicitEquation&	result1,
											            rg_ImplicitEquation&	result2)
{
	rg_Circle2D c1Tilde, c2Tilde;
	rg_Point2D z3;

	//c3 has not the smallest radius between three circles
	//so we must make c3 to have smallest radius by swapping for convienience
	shrinkCircle(c1, c2, c3, c1Tilde, c2Tilde, z3);

	rg_Circle2D w1, w2;
	w1 = transformZ2W( c1Tilde, z3 );
	w2 = transformZ2W( c2Tilde, z3 );

	//need to be added generating two tangent lines
	makeTangentLinesOf2Circles(w1, w2, result1, result2);

	return z3;
}

//this function returns tangent lines 
void QuasiTriangulation::makeTangentLinesOf2Circles(const rg_Circle2D&   circle1,
									 		        const rg_Circle2D&   circle2,
											        rg_ImplicitEquation& result1,
											        rg_ImplicitEquation& result2)
{
	rg_Circle2D w1, w2;

	//the radius of w2 is less than that of w1
	if( circle1.getRadius() < circle2.getRadius() )
	{
		w1 = circle2;
		w2 = circle1;
	}
	else
	{
		w1 = circle1;
		w2 = circle2;
	}

	rg_Point2D c2cVector = w1.getCenterPt() - w2.getCenterPt();
	rg_REAL r = w1.getRadius() - w2.getRadius();
	rg_REAL length = c2cVector.magnitude();
	rg_REAL sine = r / length;
	rg_REAL cosine = sqrt( length*length - r*r ) / length;

	//rotate theta  /  -theta
	rg_Point2D normal1( c2cVector.getX() * cosine - c2cVector.getY() * sine ,
					 c2cVector.getX() * sine + c2cVector.getY() * cosine );
	rg_Point2D normal2( c2cVector.getX() * cosine + c2cVector.getY() * sine ,
					 c2cVector.getY() * cosine - c2cVector.getX() * sine );
	normal1 = normal1.getUnitVector();
	normal2 = normal2.getUnitVector();
	
	//rotate -PI/2  /  PI/2
	normal1.setPoint( normal1.getY() , -1 * normal1.getX() );
	normal2.setPoint( -1 * normal2.getY() , normal2.getX() );

	result1.setCoeff(1, 0, normal1.getX());
	result1.setCoeff(0, 1, normal1.getY());
	result1.setCoeff(0, 0, -1 * normal1.getX() * w2.getCenterPt().getX() 
						   - normal1.getY() * w2.getCenterPt().getY() 
						   + w2.getRadius()   );

	result2.setCoeff(1, 0, normal2.getX());
	result2.setCoeff(0, 1, normal2.getY());
	result2.setCoeff(0, 0, -1 * normal2.getX() * w2.getCenterPt().getX() 
						   - normal2.getY() * w2.getCenterPt().getY() 
						   + w2.getRadius()   );

}

//
// Z = 1 / w + z3
//
rg_Circle2D QuasiTriangulation::transformW2Z(const rg_ImplicitEquation& line,
								             const rg_Point2D&          z3)
{
	//degenerate case must be considered
	//case c == 0;

	rg_REAL a = line.getCoeff(1, 0);
	rg_REAL b = line.getCoeff(0, 1);
	rg_REAL c = line.getCoeff(0, 0);

	rg_REAL tempX = -a/2./c + z3.getX();
	rg_REAL tempY =  b/2./c + z3.getY();
	rg_REAL rad = 1./2./c;

	return rg_Circle2D(tempX, tempY, rad);
}

//
// W = 1 / (z - z3)
//
rg_Circle2D QuasiTriangulation::transformZ2W(const rg_Circle2D& cTilde,
								             const rg_Point2D&  c3Tilde)
{
	rg_REAL x1 = cTilde.getCenterPt().getX();
	rg_REAL y1 = cTilde.getCenterPt().getY();
	rg_REAL x3 = c3Tilde.getX();
	rg_REAL y3 = c3Tilde.getY();

	rg_REAL p1 = pow(x1 - x3, 2.) + pow(y1 - y3, 2.) - pow(cTilde.getRadius(), 2.);

	return rg_Circle2D( (x1 - x3)/p1 , -(y1 - y3)/p1, cTilde.getRadius() / p1 );
}






void QuasiTriangulation::reportQuasiTriangulation( ofstream& fout )
{
    rg_INT numVertex = m_vertices.getSize();
    fout << "QT-VERTEX:\t" << numVertex << endl;

    QTVertex* currVertex = rg_NULL;
    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() )  {
        currVertex = m_vertices.getpEntity();
        
        rg_Point3D pt = currVertex->getBall().getCenter();
        fout << currVertex->getID() << "\t"; 
        fout << currVertex->getFirstTetrahedron()->getID() << "\t";
        fout << pt.getX() << "\t" << pt.getY() << "\t" << pt.getZ() << "\t"; 
        fout << currVertex->getBall().getRadius() << endl;
    }
    fout << endl;

    fout << "QT-TETRAHEDRON:\t" << m_tetrahedra.getSize() << endl;
    QTTetrahedron* currTetra = rg_NULL;
    m_tetrahedra.reset4Loop();
    while ( m_tetrahedra.setNext4Loop() )  {
        currTetra = m_tetrahedra.getpEntity();

        QTVertex**      vertices = currTetra->getVertices();
        QTTetrahedron** neighbors = currTetra->getNeighbors();

		rg_INT i=0;
        fout << currTetra->getID() << " :\t";
        for (i=0; i<4; i++) {
            if ( vertices[i] != rg_NULL )
                fout << vertices[i]->getID() << "\t";
            else
                fout << "0" << "\t";
        }
        fout << "|\t";
        for (i=0; i<4; i++) {
            if ( neighbors[i] != rg_NULL )
                fout << neighbors[i]->getID() << "\t";
            else
                fout << "-1" << "\t";
        }

        if ( !currTetra->isInfiniteTetrahedron() && currTetra->getVertex(3) != rg_NULL )  {        
            rg_Point3D pt = currTetra->getEmptyTangentSphere().getCenter();
            fout << pt.getX() << "\t" << pt.getY() << "\t" << pt.getZ() << "\t"; 
            fout << currTetra->getEmptyTangentSphere().getRadius();
        }
        fout << endl;
    }
    fout << endl;

    fout << "QT-GATE:\t" << m_gates.getSize() << endl;
    rg_dList<QTGate>* gates = m_gates.getGates();
    QTGate* currGate = rg_NULL;
    gates->reset4Loop();
    while ( gates->setNext4Loop() )  {
        currGate = gates->getpEntity();

        fout << currGate->getVertex(0)->getID() << "\t" << currGate->getVertex(1)->getID() << "\t";


        rg_INT depth = 0;
        rg_INT ID[2] = { currGate->getVertex(0)->getID(), currGate->getVertex(1)->getID() };
        while ( rg_TRUE ){
            if ( ID[0]<0 && ID[1]<0 )  {
                ID[0] += numVertex;
                ID[1] += numVertex;
                depth++;
            }
            else  {
                break;
            }
        }

        fout << "D-" << depth << "\t";


        fout << currGate->getBigTetrahedron()->getID() << "\t";

        rg_dList<QTTetrahedron*>* smallTetra = currGate->getSmallTetrahedra();
        smallTetra->reset4Loop();
        while ( smallTetra->setNext4Loop() )  {
            currTetra = smallTetra->getEntity();

            fout << currTetra->getID() << "\t";
        }
        fout << endl;

            
        
    }

    fout << endl;

}



void QuasiTriangulation::analyzeAnomaly(ofstream& fout)
{
    rg_INT numSingularity = m_gates.getSize();
    rg_INT numDegeneracy  = 0;
    rg_INT numMultiplicityByTwoCommonFace = 0;
    rg_INT numMultiplicityByThreeCommonFace = 0;
    rg_INT numMultiplicityByFourCommonFace = 0;

    QTTetrahedron* currTetra = rg_NULL;
    m_tetrahedra.reset4Loop();
    while ( m_tetrahedra.setNext4Loop() )  {
        currTetra = m_tetrahedra.getpEntity();

        rg_INT anomaly = currTetra->isAnomaly();
        switch( anomaly ) {
        case DEGENERACY_ANOMALY:
            numDegeneracy++;
        	break;
        case MULTIPLICITY_ANOMALY_BY_TWO_COMMON_FACE:
            numMultiplicityByTwoCommonFace++;
        	break;
        case MULTIPLICITY_ANOMALY_BY_THREE_COMMON_FACE:
        	numMultiplicityByThreeCommonFace++;
            break;
        case MULTIPLICITY_ANOMALY_BY_FOUR_COMMON_FACE:
        	numMultiplicityByFourCommonFace++;
        	break;
        default:
        	break;
        }
    }


    fout << "###  Num. Anomaly  ###" << endl;
    fout << "#_multi_2" << "\t"
         << "#_multi_3" << "\t"
         << "#_multi_4" << "\t"
         << "#_singularity" << "\t"
         << "#_degeneracy" << endl;
    fout << numMultiplicityByTwoCommonFace << "\t"
         << numMultiplicityByThreeCommonFace << "\t"
         << numMultiplicityByFourCommonFace << "\t"
         << numSingularity << "\t"
         << numDegeneracy << endl;
}




void QuasiTriangulation::readQuasiTriangulationFile( const string& qtfFilename )
{
	ifstream fin( qtfFilename.c_str() );
    fin.seekg(0, ios::end);

    streampos sp = fin.tellg();
    char* theQTFFile = new char[sp];

    fin.seekg(0, ios::beg);
    fin.read(theQTFFile, sp);
    
    streamsize rp = fin.gcount();

    unsigned int byteNumOfPDBFile = (unsigned int) sp;
    unsigned int byteNumToBeRead  = (unsigned int) rp;

    unsigned int numLineFeed = byteNumOfPDBFile - byteNumToBeRead;

    if ( numLineFeed == 0 )
        numLineFeed = INT_MAX;

    unsigned int numEmptyLine = 0;
    for ( unsigned int posInF=0; posInF<byteNumOfPDBFile; posInF++ )  {
        if ( theQTFFile[posInF]=='\n' && theQTFFile[posInF+1]=='\n' )
            numEmptyLine++;
    }
    numLineFeed -= numEmptyLine;

	const char*	seps = "\n";
	const char*	sepsInToken = " \t";
    char*	    token = NULL;
    char*	    word = NULL;

    unsigned int posInFile = 0;
    unsigned int numLine   = 0;


    QTVertex**      qtVertex      = rg_NULL;
    QTTetrahedron** qtTetrahedron = rg_NULL;

    int i=0;
    int numQTEntity = 0;
    int loadStateInQT = DEALWITH_NO_QTENTITY_IN_QTF;
    numLine++;
    token = strtok(theQTFFile, seps);   

    while ( token != NULL && numLine <= numLineFeed )
    {		
        switch (loadStateInQT) {
        case DEALWITH_QTVERTEX_IN_QTF:
            break;

        case DEALWITH_QTTETRAHEDRON_IN_QTF:
            break;

        case DEALWITH_QTGATE_IN_QTF:
            break;

        case DEALWITH_NO_QTENTITY_IN_QTF:
            numQTEntity = changeLoadStateInQT(loadStateInQT, token);

            switch (loadStateInQT) {
            case DEALWITH_QTVERTEX_IN_QTF:
                qtVertex = new QTVertex*[numQTEntity+1];
                for(i=0; i<numQTEntity+1; i++) 
                    qtVertex[i] = m_vertices.add( QTVertex(i) );
                break;

            case DEALWITH_QTTETRAHEDRON_IN_QTF:
                qtTetrahedron = new QTTetrahedron*[numQTEntity+1];
                for(i=0; i<numQTEntity+1; i++) 
                    qtTetrahedron[i] = m_tetrahedra.add( QTTetrahedron(i) );
                break;

            default:
                break;
            }

            break;

        default:
            break;
        }


        token = strtok(NULL, seps);   
    }
    
    delete [] theQTFFile;
}




int QuasiTriangulation::changeLoadStateInQT(int& loadStateInQT, char* oneLineOfQTF) 
{

    int   pos = 0;
    char substring[256];

    getNextSubstring(oneLineOfQTF, pos, substring);

    
    if ( strcmp(substring, "QTVERTEX") == 0 )  {
        loadStateInQT = DEALWITH_QTVERTEX_IN_QTF;

        getNextSubstring(oneLineOfQTF, pos, substring);
        return atoi(substring);
    }
    else if ( strcmp(substring, "END_QTVERTEX") == 0 )  {
        loadStateInQT = DEALWITH_NO_QTENTITY_IN_QTF;
    }
    else if ( strcmp(substring, "QTTETRAHEDRON") == 0 )  {
        loadStateInQT = DEALWITH_QTTETRAHEDRON_IN_QTF;

        getNextSubstring(oneLineOfQTF, pos, substring);
        return atoi(substring);
    }
    else if ( strcmp(substring, "END_QTTETRAHEDRON") == 0 )  {
        loadStateInQT = DEALWITH_NO_QTENTITY_IN_QTF;
    }
    else if ( strcmp(substring, "QTGATE") == 0 )  {
        loadStateInQT = DEALWITH_QTGATE_IN_QTF;

        getNextSubstring(oneLineOfQTF, pos, substring);
        return atoi(substring);
    }
    else if ( strcmp(substring, "END_QTGATE") == 0 )  {
        loadStateInQT = DEALWITH_NO_QTENTITY_IN_QTF;
    }
    else {
        loadStateInQT = DEALWITH_NO_QTENTITY_IN_QTF;
    }

    return 0;
}



void QuasiTriangulation::passWhiteSpace(char* source, int& pos)
{
    int sourceLength    = strlen(source);
    for ( ; pos<sourceLength; pos++ )  {
        if ( source[pos] != ' ' && source[pos] != '\t' )
            break;
    }
}



void QuasiTriangulation::getNextSubstring(char* source, int& pos, char* destination)
{
    passWhiteSpace(source, pos);
    
    int sourceLength    = strlen(source);
    int isInValidChar   = 0;
    int i = 0;

    for ( ; pos<sourceLength; pos++ )  {
        if ( source[pos] != ' ' && source[pos] != '\t' )
            destination[i++] = source[pos];
        else
            break;
    }
    destination[i] = 0;

}



void QuasiTriangulation::loadQTVertex(QTVertex** qtVertex, char* oneLineOfQTF)
{
    int   pos = 0;
    char substring[256];
    getNextSubstring(oneLineOfQTF, pos, substring);
    
    QTVertex* currVertex = qtVertex[ atoi(substring) ];

    getNextSubstring(oneLineOfQTF, pos, substring);
    rg_REAL x = atof(substring);

    getNextSubstring(oneLineOfQTF, pos, substring);
    rg_REAL y = atof(substring);

    getNextSubstring(oneLineOfQTF, pos, substring);
    rg_REAL z = atof(substring);

    getNextSubstring(oneLineOfQTF, pos, substring);
    rg_REAL radius = atof(substring);

    Ball* ball = m_balls.add(Ball(Sphere(x, y, z, radius), rg_NULL));
    currVertex->setBall( ball );
}



void QuasiTriangulation::loadQTTetrahedron(QTVertex** qtVertex, QTTetrahedron** qtTetrahedron, char* oneLineOfQTF)
{
    int   pos = 0;
    char substring[256];
    getNextSubstring(oneLineOfQTF, pos, substring);
    
    QTTetrahedron* currTetrahedron = qtTetrahedron[ atoi(substring) ];

    rg_INT vertex[4]   = {-1, -1, -1, -1};
    rg_INT neighbor[4] = {-1, -1, -1, -1};

    getNextSubstring(oneLineOfQTF, pos, substring);
    vertex[0] = atoi(substring);
    getNextSubstring(oneLineOfQTF, pos, substring);
    vertex[1] = atoi(substring);
    getNextSubstring(oneLineOfQTF, pos, substring);
    vertex[2] = atoi(substring);
    getNextSubstring(oneLineOfQTF, pos, substring);
    vertex[3] = atoi(substring);
    
    getNextSubstring(oneLineOfQTF, pos, substring);
    neighbor[0] = atoi(substring);
    getNextSubstring(oneLineOfQTF, pos, substring);
    neighbor[1] = atoi(substring);
    getNextSubstring(oneLineOfQTF, pos, substring);
    neighbor[2] = atoi(substring);
    getNextSubstring(oneLineOfQTF, pos, substring);
    neighbor[3] = atoi(substring);

    int i=0;
    for ( i=0; i<4; i++ )  {
        if ( vertex[i] != -1 ) 
            currTetrahedron->setVertex(i, qtVertex[ vertex[i] ]);
        else
            currTetrahedron->setVertex(i, rg_NULL );

        if ( neighbor[i] != -1 ) 
            currTetrahedron->setNeighbor(i, qtTetrahedron[ neighbor[i] ]);
        else
            currTetrahedron->setNeighbor(i, rg_NULL );
    }


}



void QuasiTriangulation::loadQTGate(QTVertex** qtVertex, QTTetrahedron** qtTetrahedron, char* oneLineOfQTF)
{
    int   pos = 0;
    char substring[256];
    getNextSubstring(oneLineOfQTF, pos, substring);
}



/*
//void QuasiTriangulation::reportQuasiTriangulation( ofstream& fout, const string& source )
void QuasiTriangulation::writeQuasiTriangulationFile( ofstream& fout, const string& source )
{
    time_t ltime;
    tzset();
    time( &ltime );
    char today[128];
    strcpy(today, ctime( &ltime ));
    today[ strlen(today)-1 ] = '\0';

    int    length = source.length();
    string sourceType = source.substr(length-3, 3);
    
    fout << "HEADER Quasi-Triangulation File Format (QTF)  v1.0"        << endl;
    fout << "HEADER Copyright by Voronoi Diagram Research Center"       << endl;
    //fout << "HEADER contact: http://voronoi.hanyang.ac.kr"              << endl;
    fout << "HEADER Computed by Edge-tracing algorithm v1.0"            << endl;
    fout << "HEADER " << today                                          << endl;
    fout << "SOURCE ";
    if ( sourceType=="ent" || sourceType=="pdb")
        fout << "PDB\t" << source << endl;
    else
        fout << "PLAIN\t" << source << endl;


    rg_INT numVertex = m_vertices.getSize();
    fout << "QTVERTEX\t" << (numVertex-1) << endl;

    QTVertex* currVertex = rg_NULL;
    m_vertices.reset4Loop();
    m_vertices.setNext4Loop();
    while ( m_vertices.setNext4Loop() )  {
        currVertex = m_vertices.getpEntity();
        
        rg_Point3D pt = currVertex->getBall().getCenter();
        fout << currVertex->getID() << "\t"; 
        //fout << currVertex->getFirstTetrahedron()->getID() << "\t";
        fout << pt.getX() << "\t" << pt.getY() << "\t" << pt.getZ() << "\t"; 
        fout << currVertex->getBall().getRadius() << "\t";

        if ( sourceType=="ent" || sourceType=="pdb")  {
            PDBAtom* pdbAtom = (PDBAtom*) currVertex->getProperty();
            fout << "ATOM  " << setw(5) << pdbAtom->getAtomSerial()
                 << " "     << pdbAtom->getAtomName()
                 << " "     << pdbAtom->getResidue()->getResidueName()
                 << " "     << pdbAtom->getChainID() << endl;
        }
        else  {
            fout << endl;
        }

    }
    fout << "END_QTVERTEX\t" << endl;
    //fout << endl;

    fout << "QTTETRAHEDRON\t" << m_tetrahedra.getSize() << endl;
    QTTetrahedron* currTetra = rg_NULL;
    m_tetrahedra.reset4Loop();
    while ( m_tetrahedra.setNext4Loop() )  {
        currTetra = m_tetrahedra.getpEntity();

        QTVertex**      vertices = currTetra->getVertices();
        QTTetrahedron** neighbors = currTetra->getNeighbors();

        fout << (currTetra->getID()+1) << "\t";
        rg_INT i;
        for (i=0; i<4; i++) {
            if ( vertices[i] != rg_NULL )
                fout << vertices[i]->getID() << "\t";
            else
                fout << "-1" << "\t";
        }
        for (i=0; i<4; i++) {
            if ( neighbors[i] != rg_NULL )
                fout << (neighbors[i]->getID()+1) << "\t";
            else
                fout << "-1" << "\t";
        }

        if ( !currTetra->isInfiniteTetrahedron() && currTetra->getVertex(3) != rg_NULL )  {        
            rg_Point3D pt = currTetra->getEmptyTangentSphere().getCenter();
            fout << pt.getX() << "\t" << pt.getY() << "\t" << pt.getZ() << "\t"; 
            fout << currTetra->getEmptyTangentSphere().getRadius();
        }
        fout << endl;
    }
    fout << "END_QTTETRAHEDRON\t" << endl;
    //fout << endl;


    if ( m_gates.getSize() > 0 )  {

        int i_gate = 1;
        fout << "QTGATE\t" << m_gates.getSize() << endl;
        rg_dList<QTGate>* gates = m_gates.getGates();
        QTGate* currGate = rg_NULL;
        gates->reset4Loop();
        while ( gates->setNext4Loop() )  {
            currGate = gates->getpEntity();

            fout << i_gate << "\t";
            fout << (currGate->getBigTetrahedron()->getID()+1) << "\t";
            fout << currGate->getVertex(0)->getID()            << "\t"; 
            fout << currGate->getVertex(1)->getID()            << "\t";

            rg_dList<QTTetrahedron*>* smallTetra = currGate->getSmallTetrahedra();
            smallTetra->reset4Loop();
            while ( smallTetra->setNext4Loop() )  {
                currTetra = smallTetra->getEntity();

                fout << (currTetra->getID()+1) << "\t";
            }
            fout << endl;

            i_gate++;
        }
        fout << "END_QTGATE\t" << endl;
    }
    fout << "END" << endl;

    //fout << endl;

}
*/

