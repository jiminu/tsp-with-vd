#include "SolidObject.h"
using namespace V::GeometryTier;


SolidObject::SolidObject()
{}



SolidObject::~SolidObject()
{}




rg_dList< VertexOfSolid >* SolidObject::getVertices() 
{
    return &m_vertices;
}



rg_dList< EdgeOfSolid >*   SolidObject::getEdges() 
{
    return &m_edges;
}



rg_dList< FaceOfSolid >*   SolidObject::getFaces() 
{
    return &m_faces;
}




VertexOfSolid* SolidObject::splitEdge(EdgeOfSolid* edge, const rg_Point3D& splitPt)
{
    VertexOfSolid* startVertex = edge->getStartVertex();
    VertexOfSolid* endVertex   = edge->getEndVertex();
    VertexOfSolid* insertedVertex = m_vertices.add( VertexOfSolid(m_vertices.getSize(), splitPt) );

    EdgeOfSolid*   insertedEdge = m_edges.add( EdgeOfSolid(m_edges.getSize(), startVertex, insertedVertex) );
 
    insertedEdge->setEdge( edge->getRightFace(), edge->getLeftFace(), edge, edge, edge->getRightLeg(), edge->getLeftLeg() );

    if ( edge->getLeftLeg()->getStartVertex() == startVertex )
        edge->getLeftLeg()->setRightLeg(insertedEdge);
    else
        edge->getLeftLeg()->setLeftHand(insertedEdge);

    if ( edge->getRightLeg()->getStartVertex() == startVertex )
        edge->getRightLeg()->setLeftLeg(insertedEdge);
    else
        edge->getRightLeg()->setRightHand(insertedEdge);


    startVertex->setFirstEdge( insertedEdge );
    insertedVertex->setFirstEdge( insertedEdge );

    edge->setStartVertex( insertedVertex );
    edge->setRightLeg( insertedEdge );
    edge->setLeftLeg( insertedEdge );

    return insertedVertex;
}



void SolidObject::makeEdge(FaceOfSolid* face, VertexOfSolid* startVertex, VertexOfSolid* endVertex)
{}



void SolidObject::makeUnitOctahedron()
{
    //  make vertices
    VertexOfSolid* ptrVertex[6];
    ptrVertex[0] = m_vertices.add( VertexOfSolid(0, rg_Point3D(0.0, 0.0, 1.0)) );
    ptrVertex[1] = m_vertices.add( VertexOfSolid(1, rg_Point3D(1.0, 0.0, 0.0)) );
    ptrVertex[2] = m_vertices.add( VertexOfSolid(2, rg_Point3D(0.0, 1.0, 0.0)) );
    ptrVertex[3] = m_vertices.add( VertexOfSolid(3, rg_Point3D(-1.0, 0.0, 0.0)) );
    ptrVertex[4] = m_vertices.add( VertexOfSolid(4, rg_Point3D(0.0, -1.0, 0.0)) );
    ptrVertex[5] = m_vertices.add( VertexOfSolid(5, rg_Point3D(0.0, 0.0, -1.0)) );

    //  make edges
    EdgeOfSolid*  ptrEdges[12];
    ptrEdges[0]  = m_edges.add( EdgeOfSolid( 0, ptrVertex[0], ptrVertex[1] ) );
    ptrEdges[1]  = m_edges.add( EdgeOfSolid( 1, ptrVertex[0], ptrVertex[2] ) );
    ptrEdges[2]  = m_edges.add( EdgeOfSolid( 2, ptrVertex[0], ptrVertex[3] ) );
    ptrEdges[3]  = m_edges.add( EdgeOfSolid( 3, ptrVertex[0], ptrVertex[4] ) );
    ptrEdges[4]  = m_edges.add( EdgeOfSolid( 4, ptrVertex[1], ptrVertex[2] ) );
    ptrEdges[5]  = m_edges.add( EdgeOfSolid( 5, ptrVertex[2], ptrVertex[3] ) );
    ptrEdges[6]  = m_edges.add( EdgeOfSolid( 6, ptrVertex[3], ptrVertex[4] ) );
    ptrEdges[7]  = m_edges.add( EdgeOfSolid( 7, ptrVertex[4], ptrVertex[1] ) );
    ptrEdges[8]  = m_edges.add( EdgeOfSolid( 8, ptrVertex[5], ptrVertex[1] ) );
    ptrEdges[9]  = m_edges.add( EdgeOfSolid( 9, ptrVertex[5], ptrVertex[2] ) );
    ptrEdges[10] = m_edges.add( EdgeOfSolid( 10, ptrVertex[5], ptrVertex[3] ) );
    ptrEdges[11] = m_edges.add( EdgeOfSolid( 11, ptrVertex[5], ptrVertex[4] ) );

    //  make faces
    FaceOfSolid*  ptrFaces[8];
    ptrFaces[0]  = m_faces.add( FaceOfSolid( 0, ptrEdges[0] ) );
    ptrFaces[1]  = m_faces.add( FaceOfSolid( 1, ptrEdges[1] ) );
    ptrFaces[2]  = m_faces.add( FaceOfSolid( 2, ptrEdges[2] ) );
    ptrFaces[3]  = m_faces.add( FaceOfSolid( 3, ptrEdges[3] ) );
    ptrFaces[4]  = m_faces.add( FaceOfSolid( 4, ptrEdges[9] ) );
    ptrFaces[5]  = m_faces.add( FaceOfSolid( 5, ptrEdges[10] ) );
    ptrFaces[6]  = m_faces.add( FaceOfSolid( 6, ptrEdges[11] ) );
    ptrFaces[7]  = m_faces.add( FaceOfSolid( 7, ptrEdges[8] ) );


    //  construct topology
    ptrVertex[0]->setFirstEdge(ptrEdges[0]);
    ptrVertex[1]->setFirstEdge(ptrEdges[4]);
    ptrVertex[2]->setFirstEdge(ptrEdges[5]);
    ptrVertex[3]->setFirstEdge(ptrEdges[6]);
    ptrVertex[4]->setFirstEdge(ptrEdges[7]);
    ptrVertex[5]->setFirstEdge(ptrEdges[8]);

    
    ptrEdges[0]->setEdge( ptrFaces[3], ptrFaces[0], ptrEdges[7], ptrEdges[4], ptrEdges[3], ptrEdges[1] );
    ptrEdges[1]->setEdge( ptrFaces[0], ptrFaces[1], ptrEdges[4], ptrEdges[5], ptrEdges[0], ptrEdges[2] );
    ptrEdges[2]->setEdge( ptrFaces[1], ptrFaces[2], ptrEdges[5], ptrEdges[6], ptrEdges[1], ptrEdges[3] );
    ptrEdges[3]->setEdge( ptrFaces[2], ptrFaces[3], ptrEdges[6], ptrEdges[7], ptrEdges[2], ptrEdges[0] );
    
    ptrEdges[4]->setEdge( ptrFaces[4], ptrFaces[0], ptrEdges[9], ptrEdges[1], ptrEdges[8], ptrEdges[0] );
    ptrEdges[5]->setEdge( ptrFaces[5], ptrFaces[1], ptrEdges[10], ptrEdges[2], ptrEdges[9], ptrEdges[1] );
    ptrEdges[6]->setEdge( ptrFaces[6], ptrFaces[2], ptrEdges[11], ptrEdges[3], ptrEdges[10], ptrEdges[2] );
    ptrEdges[7]->setEdge( ptrFaces[7], ptrFaces[3], ptrEdges[8], ptrEdges[0], ptrEdges[11], ptrEdges[3] );

    ptrEdges[8]->setEdge( ptrFaces[4], ptrFaces[7], ptrEdges[4], ptrEdges[7], ptrEdges[9], ptrEdges[11] );
    ptrEdges[9]->setEdge( ptrFaces[5], ptrFaces[4], ptrEdges[5], ptrEdges[4], ptrEdges[10], ptrEdges[8] );
    ptrEdges[10]->setEdge( ptrFaces[6], ptrFaces[5], ptrEdges[6], ptrEdges[5], ptrEdges[11], ptrEdges[9] );
    ptrEdges[11]->setEdge( ptrFaces[7], ptrFaces[6], ptrEdges[7], ptrEdges[6], ptrEdges[8], ptrEdges[10] );

}




void SolidObject::makeSphericalTriangularMesh(const rg_INT&     levelOfResolution, 
                                              const rg_Point3D& center, 
                                              const rg_REAL&    radius)
{
    if ( levelOfResolution < 1 )
        return;

    //  1st level of resolution : octahedron
    makeUnitOctahedron();

    VertexOfSolid* currVertex = rg_NULL;
    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() )  {
        currVertex = m_vertices.getpEntity();

        rg_Point3D pos = currVertex->getGeometry();
        pos = radius*pos + center;

        currVertex->setGeometry(pos);
    }

    if ( levelOfResolution == 1 )
        return;


    //  raise up the level of resolution
    for ( rg_INT level = 2; level <= levelOfResolution; level++ )  {

        //  split an edge in two.
        //      after splitting edge, the face consists of 6 edges.
        rg_INT numVertexOfCurrLevel = m_vertices.getSize();

        EdgeOfSolid* currEdge = rg_NULL;
        rg_INT numEdgeOfCurrLevel = m_edges.getSize();
        rg_dNode<EdgeOfSolid>* currEdgeNode  = m_edges.getFirstpNode();

        rg_INT i=0;
		for ( i=0; i<numEdgeOfCurrLevel; i++, currEdgeNode=currEdgeNode->getNext())  {

            currEdge = currEdgeNode->getpEntity();
            
            rg_Point3D startPt = currEdge->getStartVertex()->getGeometry();
            rg_Point3D endPt   = currEdge->getEndVertex()->getGeometry();

            rg_Point3D insertedPt = (startPt + endPt)/2.0;
            insertedPt.normalize();
            insertedPt = radius*insertedPt;

            VertexOfSolid* insertedVertex = splitEdge(currEdge, insertedPt);
        }

        //  split a face in four.
        FaceOfSolid* currFace = rg_NULL;

        rg_INT numFaceOfCurrLevel = m_faces.getSize();
        rg_dNode<FaceOfSolid>* currFaceNode = m_faces.getFirstpNode();

        for ( i=0; i<numFaceOfCurrLevel; i++, currFaceNode=currFaceNode->getNext())  {

            currFace = currFaceNode->getpEntity();

            ///////////////////////////////////////////////////////////////////


            EdgeOfSolid* boundingEdges[6];
            rg_INT       currEdgeIndex = 0;
            
            EdgeOfSolid* startEdge = currFace->getFirstEdge();
            if ( currFace == startEdge->getLeftFace() )  {
                //if ( startEdge->getStartVertex()->getID() == level )
                if ( startEdge->getStartVertex()->getID() < numVertexOfCurrLevel )                    
                    startEdge = startEdge->getLeftHand();
            }
            else  {
                //if ( startEdge->getEndVertex()->getID() == level )
                if ( startEdge->getEndVertex()->getID() < numVertexOfCurrLevel )
                    startEdge = startEdge->getRightLeg();
            }
            currEdge  = startEdge;

            do  {
                boundingEdges[currEdgeIndex++] = currEdge;

                if ( currFace == currEdge->getLeftFace() )
                    currEdge = currEdge->getLeftHand();
                else
                    currEdge = currEdge->getRightLeg();

            } while ( startEdge != currEdge );
            

            
            ///////////////////////////////////////////////////////////////////
            EdgeOfSolid* ptrNewEdge[3] = {rg_NULL, rg_NULL, rg_NULL};
            FaceOfSolid* ptrNewFace[3] = {rg_NULL, rg_NULL, rg_NULL};
            currEdgeIndex = 0;

			rg_INT j=0;
            for ( j=0; j<6; j+=2)  {
                VertexOfSolid* startVertex = rg_NULL;
                VertexOfSolid* endVertex   = rg_NULL;

                //if ( boundingEdges[i]->getStartVertex()->getID() == level )
                if ( boundingEdges[j]->getStartVertex()->getID() >= numVertexOfCurrLevel )
                    startVertex = boundingEdges[j]->getStartVertex();
                else
                    startVertex = boundingEdges[j]->getEndVertex();

                //if ( boundingEdges[i+1]->getStartVertex()->getID() == level )
                if ( boundingEdges[j+1]->getStartVertex()->getID() >= numVertexOfCurrLevel )
                    endVertex = boundingEdges[j+1]->getStartVertex();
                else
                    endVertex = boundingEdges[j+1]->getEndVertex();

                ptrNewEdge[currEdgeIndex] = m_edges.add( EdgeOfSolid( m_edges.getSize(), startVertex, endVertex ) );
                ptrNewFace[currEdgeIndex] = m_faces.add( FaceOfSolid( m_faces.getSize(), ptrNewEdge[currEdgeIndex] ) );

                if ( boundingEdges[j]->getLeftFace() == currFace )  {
                    boundingEdges[j]->setLeftLeg(ptrNewEdge[currEdgeIndex]);
                    boundingEdges[j]->setLeftFace(ptrNewFace[currEdgeIndex]);
                }
                else  {
                    boundingEdges[j]->setRightHand(ptrNewEdge[currEdgeIndex]);
                    boundingEdges[j]->setRightFace(ptrNewFace[currEdgeIndex]);
                }

                if ( boundingEdges[j+1]->getLeftFace() == currFace )  {
                    boundingEdges[j+1]->setLeftHand(ptrNewEdge[currEdgeIndex]);
                    boundingEdges[j+1]->setLeftFace(ptrNewFace[currEdgeIndex]);
                }
                else  {
                    boundingEdges[j+1]->setRightLeg(ptrNewEdge[currEdgeIndex]);
                    boundingEdges[j+1]->setRightFace(ptrNewFace[currEdgeIndex]);
                }

                ptrNewEdge[currEdgeIndex]->setEdge(ptrNewFace[currEdgeIndex], currFace, boundingEdges[j+1], rg_NULL, boundingEdges[j], rg_NULL);
                currEdgeIndex++;
            }

            currFace->setFirstEdge(ptrNewEdge[0]);

            ptrNewEdge[0]->setLeftHand( ptrNewEdge[1] );
            ptrNewEdge[0]->setLeftLeg( ptrNewEdge[2] );

            ptrNewEdge[1]->setLeftHand( ptrNewEdge[2] );
            ptrNewEdge[1]->setLeftLeg( ptrNewEdge[0] );

            ptrNewEdge[2]->setLeftHand( ptrNewEdge[0] );
            ptrNewEdge[2]->setLeftLeg( ptrNewEdge[1] );
        }
    }
}



void SolidObject::scaleSphericalTriangularMesh(const rg_REAL& radius)
{
    VertexOfSolid* currVertex = rg_NULL;
    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() )  {
        currVertex = m_vertices.getpEntity();

        rg_Point3D pos = currVertex->getGeometry();
        pos.normalize();
        pos = radius*pos;

        currVertex->setGeometry(pos);
    }
}



void SolidObject::translateSphericalTriangularMesh(const rg_Point3D& delta)
{
    VertexOfSolid* currVertex = rg_NULL;
    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() )  {
        currVertex = m_vertices.getpEntity();

        rg_Point3D pos = currVertex->getGeometry();
        pos += delta;

        currVertex->setGeometry(pos);
    }
}



void SolidObject::writeIndexedFacedSetOfVRML(char* firstIndention, char* indentionUnit, ofstream& fout)
{
    fout << firstIndention << "geometry IndexedFaceSet {" << endl;
	fout << firstIndention << indentionUnit << "coord Coordinate {" << endl;
	fout << firstIndention << indentionUnit << indentionUnit << "point [" << endl;

    VertexOfSolid* currVertex = rg_NULL;
    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() )  {
        currVertex = m_vertices.getpEntity();

        rg_Point3D pos = currVertex->getGeometry();
	    fout << firstIndention << indentionUnit << indentionUnit << indentionUnit;
        fout << pos.getX() << "  " << pos.getY() << "  " << pos.getZ() << "," << endl;
    }

	fout << firstIndention << indentionUnit << indentionUnit << "]" << endl;
	fout << firstIndention << indentionUnit << "}" << endl;
	fout << firstIndention << indentionUnit << "coordIndex [" << endl;

    FaceOfSolid* currFace = rg_NULL;
    m_faces.reset4Loop();
    while ( m_faces.setNext4Loop() )  {
        currFace = m_faces.getpEntity();

        rg_dList<VertexOfSolid*> boundingVertices;
        currFace->getIncidentVertexList(boundingVertices);

	    fout << firstIndention << indentionUnit << indentionUnit;
        boundingVertices.reset4Loop();
        while ( boundingVertices.setNext4Loop() )  {
            currVertex = boundingVertices.getEntity();

            fout << currVertex->getID() << ",  ";            
        }
	    fout << "-1," << endl;
    }

	fout << firstIndention << indentionUnit << "]" << endl;
	fout << firstIndention << indentionUnit << "normal Normal {" << endl;
	fout << firstIndention << indentionUnit << indentionUnit << "vector [" << endl;

    
    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() )  {
        currVertex = m_vertices.getpEntity();

        rg_Point3D pos = currVertex->getGeometry();
        pos.normalize();
	    fout << firstIndention << indentionUnit << indentionUnit << indentionUnit;
        fout << pos.getX() << "  " << pos.getY() << "  " << pos.getZ() << "," << endl;
    }


	fout << firstIndention << indentionUnit << indentionUnit<< "]" << endl;
	fout << firstIndention << indentionUnit << "}" << endl;


    fout << firstIndention << "}" << endl;

}


