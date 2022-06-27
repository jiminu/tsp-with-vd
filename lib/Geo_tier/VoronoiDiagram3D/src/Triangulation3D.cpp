#include "Triangulation3D.h"
using namespace V::GeometryTier;

#include <float.h>

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "rg_RelativeOp.h"


Triangulation3D::Triangulation3D()
{
//    FOUT.open("test.txt");
}



Triangulation3D::~Triangulation3D()
{
//    FOUT.close();
} 



void Triangulation3D::reset()
{
    m_vertices.removeAll();
    m_tetrahedra.removeAll();
}



rg_INT Triangulation3D::getNumOfVertices() const
{
    return m_vertices.getSize();
}



rg_INT Triangulation3D::getNumOfTetrahedra() const
{
    return m_tetrahedra.getSize();
}




rg_dList<T3DVertex>* Triangulation3D::getVertices()
{
    return &m_vertices;
}



rg_dList<T3DTetrahedron>* Triangulation3D::getTetrahedra()
{
    return &m_tetrahedra;
}




T3DVertex* Triangulation3D::addVertex(const T3DVertex& vertex)
{
    return m_vertices.add( vertex );
}



T3DTetrahedron* Triangulation3D::addTetrahedron(const T3DTetrahedron& tetrahedron)
{
    return m_tetrahedra.add( tetrahedron );
}



rg_FLAG Triangulation3D::locateFaceByVertexID(T3DFace& faceToLocate, const rg_INT& vertexID1, const rg_INT& vertexID2, const rg_INT& vertexID3)
{
    T3DTetrahedron* currTetra = rg_NULL;
    m_tetrahedra.reset4Loop();
    while ( m_tetrahedra.setNext4Loop() )  {
        currTetra = m_tetrahedra.getpEntity();

        if ( currTetra->locateFace(faceToLocate, vertexID1, vertexID2, vertexID3) )  
            return rg_TRUE;
    }

    return rg_FALSE;
}


void Triangulation3D::tetrahedralizePointsOnSphere(const rg_FLAG& initCondition, const rg_INT& numPts, rg_Point3D* points)
{
    //loadRandomData(numPts);
    //loadTestData();
    loadPoints(numPts, points);

    
	//  make FaceQueue;
    T3DVertex* convergingVertex = rg_NULL;
    rg_dList< rg_dList<GateForT3D>* > listOfMakingLine;
    initializeTetrahedralizationOfPointsOnSphere(initCondition, convergingVertex, listOfMakingLine);

    rg_INT countLine = 0;
    rg_INT count = 0;
    rg_dList<GateForT3D>* currLine = rg_NULL;
    while ( listOfMakingLine.getSize() > 0 )  {
        currLine = listOfMakingLine.popFront();

//        FOUT << endl;
//        reportMakingLine(currLine);
//        FOUT << endl;

        rg_dNode<GateForT3D>* currNode = currLine->getFirstpNode();
        while ( currLine->getSize() > 0 ) {
            rg_dNode<GateForT3D>* nextNode = currNode->getNext();

            if ( count == 89 )
                int aaa = 1;

            GateForT3D* currGate = currNode->getpEntity();

            T3DVertex* vertex[2]  = {currGate->getVertex(), nextNode->getpEntity()->getVertex()};
            T3DVertex* mateVertex = currGate->getMateVertex();
            
//            FOUT << countLine << "  " << count << "  v" << vertex[0]->getID() << "  v" << vertex[1]->getID();

            T3DVertex* vertexToDefineAdjacentFace = findVertexToDefineAdjacentFaceOnConvexHull(vertex[1], vertex[0], mateVertex);
            
//            FOUT  << "  (v" << vertexToDefineAdjacentFace->getID() << ")  \t";

            if ( vertexToDefineAdjacentFace->isChecked() == NOT_VISITED )  {
                updateMakingLineByNewVertex(convergingVertex, vertexToDefineAdjacentFace, currNode, currLine);
                reportMakingLine(currLine);

            }
            else if ( vertexToDefineAdjacentFace->isChecked() == ON_PROCESS ) {
                updateMakingLineByVertexOnLine(convergingVertex, vertexToDefineAdjacentFace, 
                                               currNode, currLine, listOfMakingLine);
                
                reportMakingLine(currLine);
            }

//            FOUT << endl;
            count++;
        }

        delete currLine;     

        countLine++;
//        FOUT << endl << endl;
    }

    T3DTetrahedron* currTetra = m_tetrahedra.getFirstpEntity();
    convergingVertex->setFirstTetrahedron( currTetra );
}



void Triangulation3D::loadPoints(const rg_INT& numPts, rg_Point3D* points)
{
    rg_REAL maxABSValue = 0.0;

    for ( rg_INT i=0; i<numPts; i++ )  {
        m_vertices.add( T3DVertex(i, points[i]) );
        
        rg_REAL x = points[i].getX();
        rg_REAL y = points[i].getY();
        rg_REAL z = points[i].getZ();

        if ( rg_ABS(x) > maxABSValue )
            maxABSValue = x;
        if ( rg_ABS(y) > maxABSValue )
            maxABSValue = y;
        if ( rg_ABS(z) > maxABSValue )
            maxABSValue = z;
    }

    m_decimalPoint = (rg_INT) log10(maxABSValue);


}



void Triangulation3D::initializeTetrahedralizationOfPointsOnSphere(const rg_FLAG& initCondition, 
                                                                   T3DVertex*& convergingVertex,
                                                          rg_dList< rg_dList<GateForT3D>* >& listOfMakingLine)
{
    if ( initCondition == GIVEN_FACE_ON_CONVEXHULL )  {

        T3DVertex* vertexForInitMakingLine[3] = {rg_NULL, rg_NULL, rg_NULL};
        rg_INT i = 0;
        m_vertices.reset4Loop();
        while( m_vertices.setNext4Loop() && i<3 ) {
            vertexForInitMakingLine[i++] = m_vertices.getpEntity();
        }

        rg_dList<GateForT3D>* initMakingLine = makeInitialMakingLine(convergingVertex, vertexForInitMakingLine);

        listOfMakingLine.add( initMakingLine );
    }
    else if ( initCondition == NO_INIT_CONDITION )  {
        T3DVertex* vertexForInitMakingLine[3] = {rg_NULL, rg_NULL, rg_NULL};

        /////////////////////////////////////////////////////////////////////////
        //
        //  for test
        rg_INT numVertices = m_vertices.getSize();
        T3DVertex** vertexArray = new T3DVertex*[numVertices];
        rg_INT i=0;
        m_vertices.reset4Loop();
        while ( m_vertices.setNext4Loop() )
            vertexArray[i++] = m_vertices.getpEntity();

        rg_FLAG isFound = rg_FALSE;
        for ( i=0; i<(numVertices-2); i++ )  {
            for (rg_INT j=1; j<(numVertices-1); j++ )  {
                for (rg_INT k=2; k<numVertices; k++ )  {
                    if ( k == i || k == j )
                        continue;
                
                    rg_FLAG isValid = do3VerticesDefineConvexHullFace(vertexArray[i], vertexArray[j], vertexArray[k]);

                    if ( do3VerticesDefineConvexHullFace(vertexArray[i], vertexArray[j], vertexArray[k]) )  {
                        vertexForInitMakingLine[0] = vertexArray[i];
                        vertexForInitMakingLine[1] = vertexArray[j];
                        vertexForInitMakingLine[2] = vertexArray[k];
                        isFound = rg_TRUE;
                        break;
                    }
                }
                if ( isFound ) break;
            }
            if ( isFound ) break;
        }


        if ( isFound )  {   

            rg_dList<GateForT3D>* initMakingLine = makeInitialMakingLine(convergingVertex, vertexForInitMakingLine);
            listOfMakingLine.add( initMakingLine );

//            FOUT << "Initialization: ";
//            FOUT << convergingVertex->getID() << "  (";
//            for ( int i=0; i<3; i++ )
//                FOUT << vertexForInitMakingLine[i]->getID() << "  ";
//            FOUT << ")" << endl;
//            FOUT << "\t";
//            GateForT3D* currGate = rg_NULL;
//            initMakingLine->reset4Loop();
//            while ( initMakingLine->setNext4Loop() )  {
//                currGate = initMakingLine->getpEntity();
//                FOUT << "v" << currGate->getVertex()->getID() << " - ";
//            }
//            FOUT << endl;

        }
        
        delete [] vertexArray;
    }
    else  {
    }
}



rg_dList<GateForT3D>* Triangulation3D::makeInitialMakingLine(T3DVertex*& convergingVertex, 
                                                             T3DVertex** vertexForInitMakingLine)
{
    convergingVertex = vertexForInitMakingLine[0];

    for ( rg_INT i=0; i<3; i++ )
        vertexForInitMakingLine[i]->setCheck(ON_PROCESS);
    

    rg_dList<GateForT3D>* initMakingLine = new rg_dList<GateForT3D>;
    initMakingLine->add( GateForT3D(vertexForInitMakingLine[0], vertexForInitMakingLine[2], rg_NULL) );
    initMakingLine->add( GateForT3D(vertexForInitMakingLine[1], vertexForInitMakingLine[0], rg_NULL) );
    initMakingLine->add( GateForT3D(vertexForInitMakingLine[2], vertexForInitMakingLine[1], rg_NULL) );

    
    T3DVertex* terminatingVertex     = vertexForInitMakingLine[1];
    T3DVertex* vertexForAdjacentFace = rg_NULL;

    while ( rg_TRUE ) {
        rg_dNode<GateForT3D>* currNode = initMakingLine->getLastpNode();
        GateForT3D*           currGate = currNode->getpEntity();

        T3DVertex* vertex[2]  = {currGate->getVertex(), convergingVertex};

        vertexForAdjacentFace = findVertexToDefineAdjacentFaceOnConvexHull(vertex[1], vertex[0], 
                                                                           currGate->getMateVertex());
        
        if ( vertexForAdjacentFace == terminatingVertex )  {              
            currGate->setMateVertex( vertex[1] );            
            convergingVertex->setCheck( COMPLETED );
            initMakingLine->popFront();
            break;
        }
        else {       
            vertexForAdjacentFace->setCheck( ON_PROCESS );
            initMakingLine->insertAfter( GateForT3D(vertexForAdjacentFace, vertex[0], rg_NULL), currNode );
            
            currGate->setMateVertex( vertex[1] );
        }        
    }

    return initMakingLine;
}



T3DVertex* Triangulation3D::findVertexToDefineAdjacentFaceOnConvexHull(T3DVertex* vertex1, T3DVertex* vertex2, T3DVertex* mateVertex)
{
    T3DVertex* vertexToDefineAdjacentFace = rg_NULL;

	rg_dNode<T3DVertex>* startNode = m_vertices.getFirstpNode();
	rg_dNode<T3DVertex>* currNode  = startNode;
	do  {
		T3DVertex* currVertex = currNode->getpEntity();

        if ( currVertex->isChecked() == COMPLETED ||
             currVertex == vertex1 || currVertex == vertex2 || currVertex == mateVertex )  
        {
		    currNode = currNode->getNext();
			continue;
        }
                
        if ( do3VerticesDefineConvexHullFace(vertex1, vertex2, currVertex) )  {
            vertexToDefineAdjacentFace = currVertex;
            break;
        }

		currNode = currNode->getNext();
	} while ( currNode != startNode);

    return vertexToDefineAdjacentFace;
}



rg_FLAG Triangulation3D::do3VerticesDefineConvexHullFace(T3DVertex* vertex1, T3DVertex* vertex2, T3DVertex* vertex3) 
{
    rg_Point3D point[3] = { vertex1->getPoint(), vertex2->getPoint(), vertex3->getPoint()};

	rg_dNode<T3DVertex>* startNode = m_vertices.getFirstpNode();
	rg_dNode<T3DVertex>* currNode  = startNode;
	do  {
		T3DVertex* currVertex = currNode->getpEntity();
		
        if ( currVertex == vertex1 || currVertex == vertex2 || currVertex == vertex3 )  {
		    currNode = currNode->getNext();
			continue;
        }

        rg_Point3D pointToCheck = currVertex->getPoint();

        //FOUT << "\t\t";
        //FOUT << vertex1->getID() << " " << vertex2->getID() << " " 
        //     << vertex3->getID() << " : " << currVertex->getID() << "  >>  ";
        if ( isLastPointInPlusHalfSpaceByFirstThreePoints(point, pointToCheck) )  {
            return rg_FALSE;
        }

		//  test whether startPt, endPt, and thirdPt is on boundary of convex hull.
		//  if yes, return currVertex.

		currNode = currNode->getNext();
	} while ( currNode != startNode);

    return rg_TRUE;
}



rg_FLAG Triangulation3D::isLastPointInPlusHalfSpaceByFirstThreePoints(rg_Point3D* point, const rg_Point3D& lastPoint)
{
    //  NOTE : -1 <= x, y, z <= 1
    //rg_REAL scaleFactor = 8 - m_decimalPoint;
    //rg_REAL scalingForValidSignificantFigure = pow(10, scaleFactor);
    rg_REAL scalingForValidSignificantFigure = 1.0e8;
    
    rg_INT intPt1[3] = {0, 0, 0};
    rg_INT intPt2[3] = {0, 0, 0};
    rg_INT intPt3[3] = {0, 0, 0};
    rg_INT intPt4[3] = {0, 0, 0};
    
    intPt1[0] = (rg_INT) (point[0].getX()*scalingForValidSignificantFigure + 0.5);
    intPt1[1] = (rg_INT) (point[0].getY()*scalingForValidSignificantFigure + 0.5);
    intPt1[2] = (rg_INT) (point[0].getZ()*scalingForValidSignificantFigure + 0.5);

    intPt2[0] = (rg_INT) (point[1].getX()*scalingForValidSignificantFigure + 0.5);
    intPt2[1] = (rg_INT) (point[1].getY()*scalingForValidSignificantFigure + 0.5);
    intPt2[2] = (rg_INT) (point[1].getZ()*scalingForValidSignificantFigure + 0.5);

    intPt3[0] = (rg_INT) (point[2].getX()*scalingForValidSignificantFigure + 0.5);
    intPt3[1] = (rg_INT) (point[2].getY()*scalingForValidSignificantFigure + 0.5);
    intPt3[2] = (rg_INT) (point[2].getZ()*scalingForValidSignificantFigure + 0.5);
    
    intPt4[0] = (rg_INT) (lastPoint.getX()*scalingForValidSignificantFigure + 0.5);
    intPt4[1] = (rg_INT) (lastPoint.getY()*scalingForValidSignificantFigure + 0.5);
    intPt4[2] = (rg_INT) (lastPoint.getZ()*scalingForValidSignificantFigure + 0.5);


    
    ////////////////////////////////////////////////////////////////////
    //
    //  Floating-point Arithmetic
    //
    rg_REAL realPt1[3];
    rg_REAL realPt2[3];
    rg_REAL realPt3[3];
    rg_INT i=0;
	for ( i=0; i<3; i++ )  {
        realPt1[i] = (rg_REAL) (intPt2[i] - intPt1[i]);
        realPt2[i] = (rg_REAL) (intPt3[i] - intPt1[i]);
        realPt3[i] = (rg_REAL) (intPt4[i] - intPt1[i]);
    }

    //               | a11  a12  a13 |  a11  a12
    // determinant = | a21  a22  a23 |  a21  a22
    //               | a31  a32  a33 |  a31  a32
    rg_REAL termForDet[6];
    termForDet[0] = realPt1[0]*realPt2[1]*realPt3[2];
    termForDet[1] = realPt1[1]*realPt2[2]*realPt3[0];
    termForDet[2] = realPt1[2]*realPt2[0]*realPt3[1];
    termForDet[3] = -1.0*realPt1[2]*realPt2[1]*realPt3[0];
    termForDet[4] = -1.0*realPt1[0]*realPt2[2]*realPt3[1];
    termForDet[5] = -1.0*realPt1[1]*realPt2[0]*realPt3[2];

    rg_REAL detOfFPA = 0.0;
    rg_REAL maxABS   = 0.0;
    for ( i=0; i<6; i++ )  {
        detOfFPA += termForDet[i];

        if ( rg_ABS(termForDet[i]) > maxABS )
            maxABS = rg_ABS(termForDet[i]);
    }


    rg_REAL errorUpperBound = 18*maxABS*10e-8;
    //rg_REAL errorUpperBound = 18*maxABS;

    //FOUT << "FTA  " << detOfFPA << " ( " << errorUpperBound << " )\t";

    if ( rg_ABS( detOfFPA ) > errorUpperBound )  {

        //FOUT << endl;
        
        if ( detOfFPA > 0.0 )  {
            return rg_TRUE;
        }
        else{
            return rg_FALSE;
        }
    }


    ////////////////////////////////////////////////////////////////////
    //
    //  Exact Arithmetic based on Symbolic Perturbation
    //
    IntegerBy96BIT exactPt1[3];
    IntegerBy96BIT exactPt2[3];
    IntegerBy96BIT exactPt3[3];
    IntegerBy96BIT exactPt4[3];
    for ( i=0; i<3; i++)  {
        exactPt1[i].setIntegerBy96Bit( intPt1[i] );
        exactPt2[i].setIntegerBy96Bit( intPt2[i] );
        exactPt3[i].setIntegerBy96Bit( intPt3[i] );
        exactPt4[i].setIntegerBy96Bit( intPt4[i] );
    }

    IntegerBy96BIT detOfEA;
    detOfEA = computeDeterminantForConvexHullTestBySymbolicPerturbation(
                        exactPt1, exactPt2, exactPt3, exactPt4);

    //FOUT << "EA  " << detOfEA.convertToREAL() << endl;

    if ( detOfEA.isGTZERO() == rg_TRUE )  {
        return rg_TRUE;
    }
    else  {
        return rg_FALSE;
    }
}



void Triangulation3D::updateMakingLineByNewVertex(T3DVertex* convergingVertex, 
                                                  T3DVertex* vertexForAdjacentFace, 
                                                  rg_dNode<GateForT3D>*& currNode,
                                                  rg_dList<GateForT3D>* makingLine)
{
    rg_dNode<GateForT3D>* nextNode     = currNode->getNext();
    GateForT3D*           currGate     = currNode->getpEntity();
    T3DVertex*            ptrVertex[2] = {currGate->getVertex(),  nextNode->getpEntity()->getVertex() };


    T3DTetrahedron* newTetrahedron = rg_NULL;
    newTetrahedron = m_tetrahedra.add( T3DTetrahedron(m_tetrahedra.getSize(), convergingVertex, 
                                               vertexForAdjacentFace, ptrVertex[0], ptrVertex[1]) );

    T3DTetrahedron* currTetrahedron = currGate->getTetrahedron();
    if ( currTetrahedron != rg_NULL )  {
        currTetrahedron->setNeighbor( currGate->getMateVertex(), newTetrahedron );
        newTetrahedron->setNeighbor( vertexForAdjacentFace, currTetrahedron );
    }

    makingLine->insertAfter( GateForT3D(vertexForAdjacentFace, ptrVertex[0], newTetrahedron), currNode );
    vertexForAdjacentFace->setCheck(ON_PROCESS);

    currGate->setMateVertex( ptrVertex[1] );
    currGate->setTetrahedron( newTetrahedron );

    currNode = nextNode;


//    FOUT << "NEW VERTEX  ";
}




void Triangulation3D::updateMakingLineByVertexOnLine(T3DVertex* convergingVertex, 
                                                     T3DVertex* vertexForAdjacentFace, 
                                                     rg_dNode<GateForT3D>*& currNode,
                                                     rg_dList<GateForT3D>* makingLine,
                                                     rg_dList< rg_dList<GateForT3D>* >& listOfMakingLine)
{
    rg_dNode<GateForT3D>* ptrNode[4]   = {currNode->getPrev(),      currNode,                 currNode->getNext(),      currNode->getNext()->getNext() };
    GateForT3D*           ptrGate[4]   = {ptrNode[0]->getpEntity(), ptrNode[1]->getpEntity(), ptrNode[2]->getpEntity(), ptrNode[3]->getpEntity() };
    T3DVertex*            ptrVertex[4] = {ptrGate[0]->getVertex(),  ptrGate[1]->getVertex(),  ptrGate[2]->getVertex(),  ptrGate[3]->getVertex() };
        
    GateForT3D* currGate   = ptrGate[1];

    if ( vertexForAdjacentFace == ptrVertex[3] && vertexForAdjacentFace == ptrVertex[0] )  {
        T3DTetrahedron* newTetrahedron = rg_NULL;
        newTetrahedron = m_tetrahedra.add( T3DTetrahedron(m_tetrahedra.getSize(), convergingVertex, 
                                                 vertexForAdjacentFace, ptrVertex[1], ptrVertex[2]) );

        T3DTetrahedron* currTetrahedron = ptrGate[1]->getTetrahedron();
        if ( currTetrahedron != rg_NULL )  {            
            newTetrahedron->setNeighbor( ptrVertex[3], currTetrahedron );
            currTetrahedron->setNeighbor( ptrGate[1]->getMateVertex(), newTetrahedron );
        }
        currTetrahedron = ptrGate[2]->getTetrahedron();
        if ( currTetrahedron != rg_NULL )  {            
            newTetrahedron->setNeighbor( ptrVertex[1], currTetrahedron );
            currTetrahedron->setNeighbor( ptrGate[2]->getMateVertex(), newTetrahedron );
        }
        currTetrahedron = ptrGate[3]->getTetrahedron();
        if ( currTetrahedron != rg_NULL )  {            
            newTetrahedron->setNeighbor( ptrVertex[2], currTetrahedron );
            currTetrahedron->setNeighbor( ptrGate[3]->getMateVertex(), newTetrahedron );
        }
        
        for ( rg_INT i=1; i<4; i++ )  {
            ptrVertex[i]->setCheck( COMPLETED );
            makingLine->kill( ptrNode[i] );
        }

//        FOUT << "v* = v_i-1 = v_i+2  \t";
    }
    else if ( vertexForAdjacentFace == ptrVertex[3] )  {
        T3DTetrahedron* newTetrahedron = rg_NULL;
        newTetrahedron = m_tetrahedra.add( T3DTetrahedron(m_tetrahedra.getSize(), convergingVertex, 
                                                 vertexForAdjacentFace, ptrVertex[1], ptrVertex[2]) );

        T3DTetrahedron* currTetrahedron = ptrGate[1]->getTetrahedron();
        if ( currTetrahedron != rg_NULL )  {            
            newTetrahedron->setNeighbor( vertexForAdjacentFace, currTetrahedron );
            currTetrahedron->setNeighbor( ptrGate[1]->getMateVertex(), newTetrahedron );
        }
        currTetrahedron = ptrGate[2]->getTetrahedron();
        if ( currTetrahedron != rg_NULL )  {            
            newTetrahedron->setNeighbor( ptrVertex[1], currTetrahedron );
            currTetrahedron->setNeighbor( ptrGate[2]->getMateVertex(), newTetrahedron );
        }
        
        ptrGate[1]->setTetrahedron( newTetrahedron );
        
        ptrGate[1]->setMateVertex( ptrVertex[2] );

        ptrVertex[2]->setCheck(COMPLETED);
        makingLine->kill( ptrNode[2] );

        currNode = ptrNode[3];


//        FOUT << "v* = v_i+2  \t";

    }
    else if ( vertexForAdjacentFace == ptrVertex[0] )  {
        T3DTetrahedron* newTetrahedron = rg_NULL;
        newTetrahedron = m_tetrahedra.add( T3DTetrahedron(m_tetrahedra.getSize(), convergingVertex, 
                                                 vertexForAdjacentFace, ptrVertex[1], ptrVertex[2]) );

        T3DTetrahedron* currTetrahedron = currGate->getTetrahedron();
        if ( currTetrahedron != rg_NULL )  {            
            newTetrahedron->setNeighbor( vertexForAdjacentFace, currTetrahedron );
            currTetrahedron->setNeighbor( currGate->getMateVertex(), newTetrahedron );
        }
        currTetrahedron = ptrGate[0]->getTetrahedron();
        if ( currTetrahedron != rg_NULL )  {            
            newTetrahedron->setNeighbor( ptrVertex[2], currTetrahedron );
            currTetrahedron->setNeighbor( ptrGate[0]->getMateVertex(), newTetrahedron );
        }
        ptrGate[0]->setTetrahedron( newTetrahedron );

        ptrGate[0]->setMateVertex( ptrVertex[1] );

        ptrVertex[1]->setCheck(COMPLETED);
        makingLine->kill( ptrNode[1] );

        currNode = ptrNode[2];


//        FOUT << "v* = v_i-1  \t";

    }
    else  {
        T3DTetrahedron* newTetrahedron = rg_NULL;
        newTetrahedron = m_tetrahedra.add( T3DTetrahedron(m_tetrahedra.getSize(), convergingVertex, 
                                                 vertexForAdjacentFace, ptrVertex[1], ptrVertex[2]) );

        T3DTetrahedron* currTetrahedron = currGate->getTetrahedron();
        if ( currTetrahedron != rg_NULL )  {            
            newTetrahedron->setNeighbor( vertexForAdjacentFace, currTetrahedron );
            currTetrahedron->setNeighbor( currGate->getMateVertex(), newTetrahedron );
        }
        ptrGate[1]->setTetrahedron( newTetrahedron );
        ptrGate[1]->setMateVertex( ptrVertex[2] );
        
        rg_dNode<GateForT3D>* nodeForVertexToDefineAdjacentFace = rg_NULL;
        rg_dNode<GateForT3D>* iNode = ptrNode[3]->getNext();
        while ( rg_TRUE )  {
            if ( iNode->getpEntity()->getVertex() == vertexForAdjacentFace )  {
                nodeForVertexToDefineAdjacentFace = iNode;
                break;
            }
            iNode = iNode->getNext();
        }

        rg_dList<GateForT3D>* newList = new rg_dList<GateForT3D>;
        iNode = nodeForVertexToDefineAdjacentFace;
        do  {
            newList->add( iNode->getEntity() );
            iNode = iNode->getNext();
        } while ( iNode != ptrNode[2] );

        iNode = nodeForVertexToDefineAdjacentFace;
        do  {
            makingLine->kill( iNode->getNext() );
        } while ( iNode->getNext() != ptrNode[2] );
        
        vertexForAdjacentFace->setCheck( ON_PROCESS );
        nodeForVertexToDefineAdjacentFace->getpEntity()->setTetrahedron( newTetrahedron );
        nodeForVertexToDefineAdjacentFace->getpEntity()->setMateVertex( ptrVertex[1] );

        listOfMakingLine.add( newList );

        currNode = ptrNode[2];

        
//        FOUT << "SPLIT  \t" << endl;
//        reportMakingLine( newList );

    }

}



void Triangulation3D::reportMakingLine(rg_dList<GateForT3D>* makingLine)
{
//    FOUT << endl;
//    FOUT << "\t";
//    GateForT3D* currGate = rg_NULL;
//    makingLine->reset4Loop();
//    while ( makingLine->setNext4Loop() )  {
//        currGate = makingLine->getpEntity();
//        FOUT << "v" << currGate->getVertex()->getID() << "(" << currGate->getVertex()->isChecked() << ")" << " - ";
//    }
}



void Triangulation3D::triangulatePointsGivenFaceOnConvexHull(const rg_INT& numPts, rg_Point3D* points)
{
    //loadPoints(numPts, points);

    // for test
    loadRandomData(numPts);
    //loadTestData();

    
    
	//  make FaceQueue;
	rg_dList<T3DGate> gateQueue;
	initializeTriangulationGivenFaceOnConvexHull(gateQueue);
	
//    FOUT << endl;
//    FOUT << "End of initialization" << endl;
//    FOUT << endl << endl << endl;
//    rg_INT i = 0;
//    m_vertices.reset4Loop();
//    while (m_vertices.setNext4Loop())  {
//        rg_Point3D pt = m_vertices.getpEntity()->getPoint();
//        FOUT << "v" << i << " :\t" << pt.getX() << "  " << pt.getY() << "  " << pt.getZ() << endl;
//        i++;
//    }
//    FOUT << endl << endl;



    rg_INT countGate = 0;
	while ( gateQueue.getSize() > 0 )  {

        T3DGate currGate = gateQueue.popFront();
      
//        FOUT << "gate " << countGate << " (" << gateQueue.getSize() << " ) " << " : ";
//	    T3DVertex** vertex = currGate.getVertices();
//        FOUT << vertex[0]->getID() << " " << vertex[1]->getID() << " " << vertex[2]->getID();
//        T3DVertex** vertexForCHFace = currGate.getVerticesForNextTetrahedron();
//        FOUT << "\t|\t ";
//        for ( i=0; i<3; i++ )  {
//            if ( vertexForCHFace[i] != rg_NULL )
//                FOUT << "v" << vertexForCHFace[i]->getID() << "  ";
//            else
//                FOUT << "-1  ";
//        }   
//        FOUT << endl;


        if ( currGate.getIncidentNeighbor() != rg_NULL )  {
//            FOUT << "\tdetermined before" << endl;
            continue;
        }

        
        findVertexToDefineConvexHullFaceWithEdgeOfCurrFace(currGate);
        
        selectVertexToMakeConvexHullFace(currGate, gateQueue);

        countGate++;
	}

    T3DTetrahedron* firstTetrahedron = m_tetrahedra.getFirstpEntity();
    T3DTetrahedron* mateTetrahedron = firstTetrahedron->getNeighbor(0);
    rg_INT          neighborPos = mateTetrahedron->locateNeighborPos( firstTetrahedron );
    mateTetrahedron->setNeighbor(neighborPos, rg_NULL );
    
    m_tetrahedra.popFront();
       
//    FOUT << "V :\t" << m_vertices.getSize() << endl;
//    FOUT << "T :\t" << m_tetrahedra.getSize() << endl;
}




void Triangulation3D::initializeTriangulationGivenFaceOnConvexHull(rg_dList<T3DGate>& gateQueue)
{
    T3DVertex* vertexOnInitialTetrahedron[3] = {rg_NULL, rg_NULL, rg_NULL};

    /*
    rg_INT i = 0;
    m_vertices.reset4Loop();
    while( m_vertices.setNext4Loop() && i<3 ) {
        vertexOnInitialTetrahedron[i++] = m_vertices.getpEntity();
    }
    */

    /////////////////////////////////////////////////////////////////////////
    //
    //  for test
    rg_INT numVertices = m_vertices.getSize();
    T3DVertex** vertexArray = new T3DVertex*[numVertices];
    rg_INT i=0;
    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() )
        vertexArray[i++] = m_vertices.getpEntity();

    rg_FLAG isFound = rg_FALSE;
    for ( i=0; i<(numVertices-2); i++ )  {
        for (rg_INT j=1; j<(numVertices-1); j++ )  {
            for (rg_INT k=2; k<numVertices; k++ )  {
                if ( k == i || k == j )
                    continue;
                
                rg_FLAG isValid = do3VerticesDefineConvexHullFace(vertexArray[i], vertexArray[j], vertexArray[k]);

                if ( do3VerticesDefineConvexHullFace(vertexArray[i], vertexArray[j], vertexArray[k]) )  {
                    vertexOnInitialTetrahedron[0] = vertexArray[i];
                    vertexOnInitialTetrahedron[1] = vertexArray[j];
                    vertexOnInitialTetrahedron[2] = vertexArray[k];
                    isFound = rg_TRUE;
                    break;
                }
            }

            if ( isFound )
                break;
        }

        if ( isFound )
            break;
    }

    
	T3DTetrahedron* initialTetrahedron 
        = m_tetrahedra.add( T3DTetrahedron( m_tetrahedra.getSize(), 
                                            rg_NULL, vertexOnInitialTetrahedron[0], 
                                            vertexOnInitialTetrahedron[1], vertexOnInitialTetrahedron[2]) );
    //
    //  for test
    /////////////////////////////////////////////////////////////////////////
    /*
    //  make initial tetrahedron
	T3DTetrahedron* initialTetrahedron 
        = m_tetrahedra.add( T3DTetrahedron( m_tetrahedra.getSize(), 
                                            rg_NULL, vertexOnInitialTetrahedron[0], 
                                            vertexOnInitialTetrahedron[2], vertexOnInitialTetrahedron[1]) );
    */

    for ( i=0; i<4; i++ )
        vertexOnInitialTetrahedron[i]->setFirstTetrahedron(initialTetrahedron);

	//  make gateQueue;
	gateQueue.add( T3DGate(initialTetrahedron, 0) );
}




void Triangulation3D::findVertexToDefineConvexHullFaceWithEdgeOfCurrFace(T3DGate& currGate)
{
	T3DVertex** vertex = currGate.getVertices();

	T3DVertex*  vertexOnCurrFace[4];
	vertexOnCurrFace[0] = vertex[0];
	vertexOnCurrFace[1] = vertex[1];
	vertexOnCurrFace[2] = vertex[2];
	vertexOnCurrFace[3] = vertex[0];

    //for ( rg_INT i=0; i<3; i++ )
    //    vertex[i]->setCheck(rg_TRUE);

    //  edge 1 : vertexOnCurrFace[0] -> vertexOnCurrFace[1]
    //  edge 2 : vertexOnCurrFace[1] -> vertexOnCurrFace[2]
    //  edge 3 : vertexOnCurrFace[2] -> vertexOnCurrFace[3]

    rg_INT i = 0;
	for ( i=0; i<3; i++ )  {

        //FOUT << "\tedge " << vertexOnCurrFace[i]->getID() << " -> " << vertexOnCurrFace[i+1]->getID() << endl;
        
        if ( currGate.isValidVertexForNextTetrahedron(i) == rg_TRUE )  {
            //FOUT << "\t\t";
            //FOUT << vertexOnCurrFace[i]->getID() << " " << vertexOnCurrFace[i+1]->getID() << " ";
            //if ( currGate.getVertexForNextTetrahedron(i) != rg_NULL )
            //    FOUT << currGate.getVertexForNextTetrahedron(i)->getID() << " :  determined before" << endl;
            //else
            //    FOUT << "-1" << " :  determined before" << endl;

            continue;
        }
        
	    rg_dNode<T3DVertex>* startNode = m_vertices.getFirstpNode();
	    rg_dNode<T3DVertex>* currNode  = startNode;
	    do  {
		    T3DVertex* currVertex = currNode->getpEntity();

            if ( currVertex == vertex[0] || currVertex == vertex[1] || currVertex == vertex[2] 
                 || currVertex->isChecked() == rg_TRUE )  {
		        currNode = currNode->getNext();
			    continue;
            }
            

            
            if ( do3VerticesDefineConvexHullFace(vertexOnCurrFace[i], vertexOnCurrFace[i+1], currVertex) )  {
                currGate.setVertexForNextTetrahedron(i, currVertex);
                break;
            }

            //FOUT << endl;

		    //  test whether startPt, endPt, and thirdPt is on boundary of convex hull.
		    //  if yes, return currVertex.

		    currNode = currNode->getNext();
	    } while ( currNode != startNode);

        //FOUT << endl;
    }
}




void Triangulation3D::selectVertexToMakeConvexHullFace(T3DGate&     currGate, 
                                                       rg_dList<T3DGate>& gateQueue)
{

    T3DVertex** vertex = currGate.getVertices();
    T3DVertex** vertexForCHFace = currGate.getVerticesForNextTetrahedron();
    // vertexForCHFace[0] for edge (v[0] - v[1])
    // vertexForCHFace[1] for edge (v[1] - v[2])
    // vertexForCHFace[2] for edge (v[2] - v[0])

    ///////////////////////////////////////////////////////////////////////////
    //
    /*
    FOUT << "\tGate:\t";
    FOUT << "T" << currGate.getTetrahedron()->getID() 
                << "(" << currGate.getMateVertexPos() << " ) " << "\tV: ";
    for (int i=0; i<3; i++ )  {
        if ( vertex[i] != rg_NULL )
            FOUT << "v" << vertex[i]->getID() << "  ";
        else
            FOUT << "-1  ";
    }   
    FOUT << "\t|\tV4nextT: ";
    for ( i=0; i<3; i++ )  {
        if ( vertexForCHFace[i] != rg_NULL )
            FOUT << "v" << vertexForCHFace[i]->getID() << "  ";
        else
            FOUT << "-1  ";
    }   
    FOUT << endl;
    */
    //
    ///////////////////////////////////////////////////////////////////////////




	if ( vertexForCHFace[0] == vertexForCHFace[1] && vertexForCHFace[0] == vertexForCHFace[2] )  {

        makeTetrahedron(currGate, 0);
        
        vertexForCHFace[0]->setCheck(rg_TRUE);
	}
	else if ( vertexForCHFace[0] == vertexForCHFace[1] && vertexForCHFace[0] != rg_NULL )  {
        // vertexForCHFace[0] for edge (v[0] - v[1])
        // vertexForCHFace[1] for edge (v[1] - v[2])
        
        //vertex[1]->setCheck(rg_TRUE);

        T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 0);
        T3DGate*        ptrGate        = rg_NULL;

        if ( currGate.getIndexOfPriorEdge() == 2 )          
            ptrGate = makeGate(gateQueue, newTetrahedron, 2, rg_TRUE);
        else  
            ptrGate = makeGate(gateQueue, newTetrahedron, 2, rg_FALSE);

        if ( ptrGate != rg_NULL )  {
            ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[2] );
            ptrGate->setIndexOfPriorEdge(1);
        }
        else
            vertexForCHFace[1]->setCheck(rg_TRUE);

	}
	else if ( vertexForCHFace[1] == vertexForCHFace[2] && vertexForCHFace[1] != rg_NULL )  {
        // vertexForCHFace[1] for edge (v[1] - v[2])
        // vertexForCHFace[2] for edge (v[2] - v[0])

        //vertex[2]->setCheck(rg_TRUE);

        T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 1);
        T3DGate*        ptrGate        = rg_NULL;

        if ( currGate.getIndexOfPriorEdge() == 0 )          
            ptrGate = makeGate(gateQueue, newTetrahedron, 3, rg_TRUE);
        else  
            ptrGate = makeGate(gateQueue, newTetrahedron, 3, rg_FALSE);

        if ( ptrGate != rg_NULL )  {
            ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[0] );
            ptrGate->setIndexOfPriorEdge(1);
        }
        else
            vertexForCHFace[2]->setCheck(rg_TRUE);

	}
	else if ( vertexForCHFace[0] == vertexForCHFace[2] && vertexForCHFace[0] != rg_NULL )  {
        // vertexForCHFace[0] for edge (v[0] - v[1])
        // vertexForCHFace[2] for edge (v[2] - v[0])
        //vertex[0]->setCheck(rg_TRUE);
 
        T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 2);
        T3DGate*        ptrGate        = rg_NULL;

        if ( currGate.getIndexOfPriorEdge() == 1 )          
            ptrGate = makeGate(gateQueue, newTetrahedron, 1, rg_TRUE);
        else  
            ptrGate = makeGate(gateQueue, newTetrahedron, 1, rg_FALSE);

        if ( ptrGate != rg_NULL )  {
            ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[1] );
            ptrGate->setIndexOfPriorEdge(1);
        }
        else
            vertexForCHFace[0]->setCheck(rg_TRUE);
	}
	else  {
        if ( currGate.getIndexOfPriorEdge() == 0 )  { 
		    if ( vertexForCHFace[1] != rg_NULL )  {
                // vertexForCHFace[1] for edge (v[1] - v[2])

                T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 1);
                T3DGate*        ptrGate        = rg_NULL;
                
                ptrGate = makeGate(gateQueue, newTetrahedron, 3, rg_TRUE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[0] );
                    ptrGate->setIndexOfPriorEdge(1);
                }

                ptrGate = makeGate(gateQueue, newTetrahedron, 2, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[2] );
                    ptrGate->setIndexOfPriorEdge(2);
                }
		    }
		    else if ( vertexForCHFace[2] != rg_NULL )  {
                // vertexForCHFace[2] for edge (v[2] - v[0])

                T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 2);
                T3DGate*        ptrGate        = rg_NULL;
                
                ptrGate = makeGate(gateQueue, newTetrahedron, 3, rg_TRUE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[0] );
                    ptrGate->setIndexOfPriorEdge(1);
                }

                ptrGate = makeGate(gateQueue, newTetrahedron, 1, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[1] );
                    ptrGate->setIndexOfPriorEdge(0);
                }
		    }            
        }
        else if ( currGate.getIndexOfPriorEdge() == 1 )  { 
		    if ( vertexForCHFace[0] != rg_NULL )  {
                // vertexForCHFace[0] for edge (v[0] - v[1])

                T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 0);

                T3DGate* ptrGate = makeGate(gateQueue, newTetrahedron, 1, rg_TRUE);

                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[1] );
                    ptrGate->setIndexOfPriorEdge(1);
                }

                ptrGate = makeGate(gateQueue, newTetrahedron, 2, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[2] );
                    ptrGate->setIndexOfPriorEdge(0);
                }
		    }
		    else if ( vertexForCHFace[2] != rg_NULL )  {
                // vertexForCHFace[2] for edge (v[2] - v[0])

                T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 2);
                T3DGate*        ptrGate        = rg_NULL;
                
                ptrGate = makeGate(gateQueue, newTetrahedron, 1, rg_TRUE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[1] );
                    ptrGate->setIndexOfPriorEdge(1);
                }

                ptrGate = makeGate(gateQueue, newTetrahedron, 3, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[0] );
                    ptrGate->setIndexOfPriorEdge(2);
                }
		    }            
        }
        else if ( currGate.getIndexOfPriorEdge() == 2 )  { 
		    if ( vertexForCHFace[0] != rg_NULL )  {
                // vertexForCHFace[0] for edge (v[0] - v[1])

                T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 0);

                T3DGate* ptrGate = makeGate(gateQueue, newTetrahedron, 2, rg_TRUE);

                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[2] );
                    ptrGate->setIndexOfPriorEdge(1);
                }

                ptrGate = makeGate(gateQueue, newTetrahedron, 1, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[1] );
                    ptrGate->setIndexOfPriorEdge(2);
                }
		    }
		    else if ( vertexForCHFace[1] != rg_NULL )  {
                // vertexForCHFace[1] for edge (v[1] - v[2])

                T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 1);

                T3DGate* ptrGate = makeGate(gateQueue, newTetrahedron, 2, rg_TRUE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[2] );
                    ptrGate->setIndexOfPriorEdge(1);
                }

                ptrGate = makeGate(gateQueue, newTetrahedron, 3, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[0] );
                    ptrGate->setIndexOfPriorEdge(0);
                }
		    }
        }
        else  {
		    if ( vertexForCHFace[0] != rg_NULL )  {
                // vertexForCHFace[0] for edge (v[0] - v[1])

                T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 0);

                T3DGate* ptrGate = makeGate(gateQueue, newTetrahedron, 1, rg_FALSE);

                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[1] );
                    ptrGate->setIndexOfPriorEdge(2);
                }

                ptrGate = makeGate(gateQueue, newTetrahedron, 2, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[2] );
                    ptrGate->setIndexOfPriorEdge(0);
                }
		    }
		    else if ( vertexForCHFace[1] != rg_NULL )  {
                // vertexForCHFace[1] for edge (v[1] - v[2])

                T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 1);

                T3DGate* ptrGate = makeGate(gateQueue, newTetrahedron, 2, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[2] );
                    ptrGate->setIndexOfPriorEdge(2);
                }

                ptrGate = makeGate(gateQueue, newTetrahedron, 3, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[0] );
                    ptrGate->setIndexOfPriorEdge(0);
                }
		    }
		    else if ( vertexForCHFace[2] != rg_NULL )  {
                // vertexForCHFace[2] for edge (v[2] - v[0])

                T3DTetrahedron* newTetrahedron = makeTetrahedron(currGate, 2);

                T3DGate* ptrGate = makeGate(gateQueue, newTetrahedron, 3, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[0] );
                    ptrGate->setIndexOfPriorEdge(2);
                }

                ptrGate = makeGate(gateQueue, newTetrahedron, 1, rg_FALSE);
                if ( ptrGate != rg_NULL )  {
                    ptrGate->setVertexForNextTetrahedron(1, vertexForCHFace[1] );
                    ptrGate->setIndexOfPriorEdge(0);
                }
		    }
        }
	}

    /*

    T3DTetrahedron* lastTetrahedron = m_tetrahedra.getLastpEntity();
    FOUT << "\tmake Tetrahedron:\tT" << lastTetrahedron->getID() << "\t";
    for ( i=0; i<4; i++ )  {
        T3DVertex* vertexLT = lastTetrahedron->getVertex(i);
        if ( vertexLT != rg_NULL )
            FOUT << vertexLT->getID() << "  ";
        else
            FOUT << "-1  ";
    }
    FOUT << endl << endl;
    */
}



T3DTetrahedron* Triangulation3D::makeTetrahedron(T3DGate& currGate, const rg_INT& edgePos)
{
	T3DVertex** vertexOnCurrFace           = currGate.getVertices();
	T3DVertex*  vertexToMakeConvexHullFace = currGate.getVertexForNextTetrahedron(edgePos);

	T3DTetrahedron* newTetrahedron = m_tetrahedra.add( T3DTetrahedron(m_tetrahedra.getSize()) );
	newTetrahedron->setVertices( vertexToMakeConvexHullFace, vertexOnCurrFace[0],
								 vertexOnCurrFace[1],        vertexOnCurrFace[2] );

    T3DTetrahedron* mateTetrahedron = currGate.getTetrahedron();
	mateTetrahedron->setNeighbor( currGate.getMateVertexPos(), newTetrahedron );
	newTetrahedron->setNeighbor( 0, mateTetrahedron );

    vertexToMakeConvexHullFace->setFirstTetrahedron(mateTetrahedron);
    
    return newTetrahedron;
}



T3DGate* Triangulation3D::makeGate(rg_dList<T3DGate>& gateQueue, 
                                   T3DTetrahedron* tetrahedron, 
                                   const rg_INT& mateVertexPos, const rg_FLAG& priority)
{
    T3DGate* ptrGate = rg_NULL;
    T3DGate  newGate(tetrahedron, mateVertexPos);

//    FOUT << "\t\tnew gate\t";
//	T3DVertex** vertexOnNewGate = newGate.getVertices();
//    FOUT << vertexOnNewGate[0]->getID() << " " << vertexOnNewGate[1]->getID() << " " << vertexOnNewGate[2]->getID();

    
    T3DGate* sameGate = isThereSameGate(gateQueue, newGate);
    if ( sameGate != rg_NULL )  {
        T3DTetrahedron* mateTetrahedron = sameGate->getTetrahedron();
	    mateTetrahedron->setNeighbor( sameGate->getMateVertexPos(), tetrahedron );
	    tetrahedron->setNeighbor( mateVertexPos, mateTetrahedron );

//        FOUT << "\tdon't insert" << endl;    
    }
    else  {
        if ( priority )
            ptrGate = gateQueue.addHead(newGate);
        else
            ptrGate = gateQueue.add(newGate);

//        FOUT << endl;
    }

    return ptrGate;
}


T3DGate* Triangulation3D::isThereSameGate(rg_dList<T3DGate>& gateQueue, const T3DGate& newGate)
{
    T3DGate* currGate = rg_NULL;

    gateQueue.reset4Loop();
    while ( gateQueue.setNext4Loop() )  {
        currGate = gateQueue.getpEntity();

        if ( currGate->compareGateByVertices(newGate) )
            return currGate;
    }

    return rg_NULL;
}



T3DGate* Triangulation3D::isThereSameGate(rg_dList<T3DGate>& gateQueue, T3DVertex** vertexOfGate)
{
    T3DGate* currGate = rg_NULL;

    gateQueue.reset4Loop();
    while ( gateQueue.setNext4Loop() )  {
        currGate = gateQueue.getpEntity();

        if ( currGate->compareGateByVertices(vertexOfGate) )
            return currGate;
    }

    return rg_NULL;
}



rg_INT Triangulation3D::exploreConnectionBetweenGatesAndVertices(rg_dList<T3DGate>& gateQueue, 
                                                                 T3DGate&     currGate)
{
    rg_INT  highScoreOfConnection    = -3;
    rg_INT  posOfVertexWithHighScore = 0;

    T3DVertex** vertex          = currGate.getVertices();
    T3DVertex** vertexForCHFace = currGate.getVerticesForNextTetrahedron();

    T3DVertex* nextGate[6][3] = { {vertexForCHFace[0], vertex[1], vertex[2]},
                                  {vertexForCHFace[0], vertex[2], vertex[0]},
                                  {vertexForCHFace[1], vertex[0], vertex[1]},
                                  {vertexForCHFace[1], vertex[2], vertex[0]},
                                  {vertexForCHFace[2], vertex[1], vertex[2]},
                                  {vertexForCHFace[2], vertex[0], vertex[1]} };
    T3DGate* oldGate[6];
    rg_INT   score[6] = {-1, -1, -1, -1, -1, -1};    
    rg_INT i=0;
	for ( i=0; i<6; i++ )  {
        if ( nextGate[i][0] == rg_NULL )
            continue;

        oldGate[i] = isThereSameGate(gateQueue, nextGate[i]);
        if ( oldGate[i] == rg_NULL )
            score[i] = 0;
        else
            score[i] = 1;
    }

    rg_INT scoreByEdge[3] = { score[0]+score[1], score[2]+score[3], score[4]+score[5] };

    for ( i=0; i<3; i++ )  {
        if ( highScoreOfConnection < scoreByEdge[i] )  {
            highScoreOfConnection    = scoreByEdge[i];
            posOfVertexWithHighScore = i;
        }
    }

    return posOfVertexWithHighScore;
}


void Triangulation3D::loadRandomData(const rg_INT& numVertices)
{
    ofstream fout;
    fout.open("T3Ddata.txt");
    fout << numVertices << endl;

    rg_INT  upperBound = 200;
    rg_INT  translate  = -100;
    rg_REAL radius     = 100.0;
    srand( (unsigned)time(NULL));
    for ( rg_INT i=0; i<numVertices; i++)  {  
	    rg_INT x = rand()%upperBound + translate;
		rg_INT y = rand()%upperBound + translate;
		rg_INT z = rand()%upperBound + translate;

        rg_Point3D point(x, y, z);
        point.normalize();
        //point = point*radius;

        m_vertices.add( T3DVertex(i, point) );

        fout << i << "\t" << point.getX() << "\t" << point.getY() << "\t" << point.getZ() << endl;
    }

    fout.close();
}




void Triangulation3D::loadTestData()
{
    
    const int numPts = 6;
    rg_Point3D point[numPts];
    point[0].setPoint( 0.51519000503999712,  0.24497744320502832,  0.82131925036955922);
    point[1].setPoint( -0.14663161591296409,  0.92829314969163801,  0.34171771603229956);
    point[2].setPoint( 0.92829314969163779,  0.34171771603229983, -0.14663161591296409);
    point[3].setPoint( 0.24497744320502809,  0.82131925036955944,  0.51519000503999712);
    point[4].setPoint( 0.82131925036955922,  0.51519000503999746,  0.24497744320502815);
    point[5].setPoint( 0.34171771603229956, -0.14663161591296406,  0.92829314969163790);
    /*
    const int numPts = 4;
    rg_Point3D point[numPts];
    point[0].setPoint( 0.51519000503999712,  0.24497744320502832,  0.82131925036955922);
    point[1].setPoint( -0.14663161591296409,  0.92829314969163801,  0.34171771603229956);
    point[2].setPoint( 0.92829314969163779,  0.34171771603229983, -0.14663161591296409);
    point[3].setPoint( 0.34171771603229956, -0.14663161591296406,  0.92829314969163790);
*/

    
    for ( rg_INT i=0; i<numPts; i++)  {  
        m_vertices.add( T3DVertex(i, point[i]) );
    }

    testDiscrimentFunction();
//    FOUT << endl << endl;
}



void Triangulation3D::testDiscrimentFunction()
{
    rg_INT numPts = m_vertices.getSize();
    rg_Point3D* pt = new rg_Point3D[numPts];
    rg_INT i = 0;
    m_vertices.reset4Loop();
    while( m_vertices.setNext4Loop() )  {
        pt[i++] = m_vertices.getpEntity()->getPoint();
    }

    for ( i=0; i<(numPts-2); i++)  {
        for (rg_INT j=i+1; j<(numPts-1); j++)  {
            for (rg_INT k=j+1; k<numPts; k++)  {
                for (rg_INT l=0; l<numPts; l++ )  {
                    if ( l == i || l == j || l == k )
                        continue;
                    
                    rg_REAL detOriFPA  = computeDetermineOriginalFloatingPointArithmetic(pt[i], pt[j], pt[k], pt[l]);
                    rg_REAL detFullFPA = computeDetermineFullFloatingPointArithmetic(pt[i], pt[j], pt[k], pt[l]);
                    rg_REAL errorUB;
                    rg_REAL detFiniteFPA = computeDetermineFiniteFloatingPointArithmetic(errorUB, pt[i], pt[j], pt[k], pt[l]);

//                    FOUT << i << " " << j << " " << k << " : " << l << "  >>  ";
//                    FOUT << detOriFPA << "\t";
//                    FOUT << detFullFPA << "\t";
//                    FOUT << detFiniteFPA << "  ( " << errorUB << " )" << endl;
                }

//                FOUT << endl;
            }
        }
    }
}



rg_REAL Triangulation3D::computeDetermineOriginalFloatingPointArithmetic(const rg_Point3D& pt1, const rg_Point3D& pt2, 
                                                            const rg_Point3D& pt3, const rg_Point3D& pt4)
{
    //  NOTE : -1 <= x, y, z <= 1
    rg_REAL scalingForValidSignificantFigure = 1.0e8;
    
    rg_REAL scaledPt1[3] = {0, 0, 0};
    rg_REAL scaledPt2[3] = {0, 0, 0};
    rg_REAL scaledPt3[3] = {0, 0, 0};
    rg_REAL scaledPt4[3] = {0, 0, 0};
    
    scaledPt1[0] = pt1.getX();
    scaledPt1[1] = pt1.getY();
    scaledPt1[2] = pt1.getZ();

    scaledPt2[0] = pt2.getX();
    scaledPt2[1] = pt2.getY();
    scaledPt2[2] = pt2.getZ();

    scaledPt3[0] = pt3.getX();
    scaledPt3[1] = pt3.getY();
    scaledPt3[2] = pt3.getZ();
    
    scaledPt4[0] = pt4.getX();
    scaledPt4[1] = pt4.getY();
    scaledPt4[2] = pt4.getZ();


    
    ////////////////////////////////////////////////////////////////////
    //
    //  Floating-point Arithmetic
    //
    rg_REAL realPt1[3];
    rg_REAL realPt2[3];
    rg_REAL realPt3[3];
    rg_INT i=0;
	for ( i=0; i<3; i++ )  {
        realPt1[i] = scaledPt2[i] - scaledPt1[i];
        realPt2[i] = scaledPt3[i] - scaledPt1[i];
        realPt3[i] = scaledPt4[i] - scaledPt1[i];
    }

    //               | a11  a12  a13 |  a11  a12
    // determinant = | a21  a22  a23 |  a21  a22
    //               | a31  a32  a33 |  a31  a32
    rg_REAL termForDet[6];
    termForDet[0] = realPt1[0]*realPt2[1]*realPt3[2];
    termForDet[1] = realPt1[1]*realPt2[2]*realPt3[0];
    termForDet[2] = realPt1[2]*realPt2[0]*realPt3[1];
    termForDet[3] = -1.0*realPt1[2]*realPt2[1]*realPt3[0];
    termForDet[4] = -1.0*realPt1[0]*realPt2[2]*realPt3[1];
    termForDet[5] = -1.0*realPt1[1]*realPt2[0]*realPt3[2];

    rg_REAL detOfFullFPA = 0.0;
    //rg_REAL maxABS   = 0.0;
    for ( i=0; i<6; i++ )  {
        detOfFullFPA += termForDet[i];

        //if ( rg_ABS(termForDet[i]) > maxABS )
        //    maxABS = rg_ABS(termForDet[i]);
    }

    return detOfFullFPA;
}


rg_REAL Triangulation3D::computeDetermineFullFloatingPointArithmetic(const rg_Point3D& pt1, const rg_Point3D& pt2, 
                                                                     const rg_Point3D& pt3, const rg_Point3D& pt4)
{
    //  NOTE : -1 <= x, y, z <= 1
    rg_REAL scalingForValidSignificantFigure = 1.0e8;
    
    rg_REAL scaledPt1[3] = {0, 0, 0};
    rg_REAL scaledPt2[3] = {0, 0, 0};
    rg_REAL scaledPt3[3] = {0, 0, 0};
    rg_REAL scaledPt4[3] = {0, 0, 0};
    
    scaledPt1[0] = pt1.getX()*scalingForValidSignificantFigure;
    scaledPt1[1] = pt1.getY()*scalingForValidSignificantFigure;
    scaledPt1[2] = pt1.getZ()*scalingForValidSignificantFigure;

    scaledPt2[0] = pt2.getX()*scalingForValidSignificantFigure;
    scaledPt2[1] = pt2.getY()*scalingForValidSignificantFigure;
    scaledPt2[2] = pt2.getZ()*scalingForValidSignificantFigure;

    scaledPt3[0] = pt3.getX()*scalingForValidSignificantFigure;
    scaledPt3[1] = pt3.getY()*scalingForValidSignificantFigure;
    scaledPt3[2] = pt3.getZ()*scalingForValidSignificantFigure;
    
    scaledPt4[0] = pt4.getX()*scalingForValidSignificantFigure;
    scaledPt4[1] = pt4.getY()*scalingForValidSignificantFigure;
    scaledPt4[2] = pt4.getZ()*scalingForValidSignificantFigure;


    
    ////////////////////////////////////////////////////////////////////
    //
    //  Floating-point Arithmetic
    //
    rg_REAL realPt1[3];
    rg_REAL realPt2[3];
    rg_REAL realPt3[3];
    rg_INT i=0;
	for ( i=0; i<3; i++ )  {
        realPt1[i] = scaledPt2[i] - scaledPt1[i];
        realPt2[i] = scaledPt3[i] - scaledPt1[i];
        realPt3[i] = scaledPt4[i] - scaledPt1[i];
    }

    //               | a11  a12  a13 |  a11  a12
    // determinant = | a21  a22  a23 |  a21  a22
    //               | a31  a32  a33 |  a31  a32
    rg_REAL termForDet[6];
    termForDet[0] = realPt1[0]*realPt2[1]*realPt3[2];
    termForDet[1] = realPt1[1]*realPt2[2]*realPt3[0];
    termForDet[2] = realPt1[2]*realPt2[0]*realPt3[1];
    termForDet[3] = -1.0*realPt1[2]*realPt2[1]*realPt3[0];
    termForDet[4] = -1.0*realPt1[0]*realPt2[2]*realPt3[1];
    termForDet[5] = -1.0*realPt1[1]*realPt2[0]*realPt3[2];

    rg_REAL detOfFullFPA = 0.0;
    //rg_REAL maxABS   = 0.0;
    for ( i=0; i<6; i++ )  {
        detOfFullFPA += termForDet[i];

        //if ( rg_ABS(termForDet[i]) > maxABS )
        //    maxABS = rg_ABS(termForDet[i]);
    }

    return detOfFullFPA;
}



rg_REAL Triangulation3D::computeDetermineFiniteFloatingPointArithmetic(rg_REAL& errorUpperBound,
                                                                       const rg_Point3D& pt1, const rg_Point3D& pt2, 
                                                                       const rg_Point3D& pt3, const rg_Point3D& pt4)
{
    //  NOTE : -1 <= x, y, z <= 1
    rg_REAL scalingForValidSignificantFigure = 1.0e8;
    
    rg_INT intPt1[3] = {0, 0, 0};
    rg_INT intPt2[3] = {0, 0, 0};
    rg_INT intPt3[3] = {0, 0, 0};
    rg_INT intPt4[3] = {0, 0, 0};
    
    intPt1[0] = (rg_INT) (pt1.getX()*scalingForValidSignificantFigure);
    intPt1[1] = (rg_INT) (pt1.getY()*scalingForValidSignificantFigure);
    intPt1[2] = (rg_INT) (pt1.getZ()*scalingForValidSignificantFigure);

    intPt2[0] = (rg_INT) (pt2.getX()*scalingForValidSignificantFigure);
    intPt2[1] = (rg_INT) (pt2.getY()*scalingForValidSignificantFigure);
    intPt2[2] = (rg_INT) (pt2.getZ()*scalingForValidSignificantFigure);

    intPt3[0] = (rg_INT) (pt3.getX()*scalingForValidSignificantFigure);
    intPt3[1] = (rg_INT) (pt3.getY()*scalingForValidSignificantFigure);
    intPt3[2] = (rg_INT) (pt3.getZ()*scalingForValidSignificantFigure);
    
    intPt4[0] = (rg_INT) (pt4.getX()*scalingForValidSignificantFigure);
    intPt4[1] = (rg_INT) (pt4.getY()*scalingForValidSignificantFigure);
    intPt4[2] = (rg_INT) (pt4.getZ()*scalingForValidSignificantFigure);


    
    ////////////////////////////////////////////////////////////////////
    //
    //  Floating-point Arithmetic
    //
    rg_REAL realPt1[3];
    rg_REAL realPt2[3];
    rg_REAL realPt3[3];
    rg_INT i=0;
	for ( i=0; i<3; i++ )  {
        realPt1[i] = (rg_REAL) (intPt2[i] - intPt1[i]);
        realPt2[i] = (rg_REAL) (intPt3[i] - intPt1[i]);
        realPt3[i] = (rg_REAL) (intPt4[i] - intPt1[i]);
    }

    //               | a11  a12  a13 |  a11  a12
    // determinant = | a21  a22  a23 |  a21  a22
    //               | a31  a32  a33 |  a31  a32
    rg_REAL termForDet[6];
    termForDet[0] = realPt1[0]*realPt2[1]*realPt3[2];
    termForDet[1] = realPt1[1]*realPt2[2]*realPt3[0];
    termForDet[2] = realPt1[2]*realPt2[0]*realPt3[1];
    termForDet[3] = -1.0*realPt1[2]*realPt2[1]*realPt3[0];
    termForDet[4] = -1.0*realPt1[0]*realPt2[2]*realPt3[1];
    termForDet[5] = -1.0*realPt1[1]*realPt2[0]*realPt3[2];

    rg_REAL detOfFPA = 0.0;
    rg_REAL maxABS   = 0.0;
    for ( i=0; i<6; i++ )  {
        detOfFPA += termForDet[i];

        if ( rg_ABS(termForDet[i]) > maxABS )
            maxABS = rg_ABS(termForDet[i]);
    }

    errorUpperBound = 18*maxABS;//*1.0e-8;

    return detOfFPA;
}

