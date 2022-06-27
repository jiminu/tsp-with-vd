#include "GeometricConverter.h"
using namespace V::GeometryTier;

GeometricConverter::GeometricConverter()
{
}



GeometricConverter::~GeometricConverter()
{
}



void GeometricConverter::convertQuasiTriangulationIntoVoronoiDiagram(const BetaUniverse& QTIneIWDS, rg_SphereSetVoronoiDiagram& VD)
{
    map<BetaCell*, VDVertex*>  mapCellInBUToVertexInVD;
    map<BetaFace*, VDEdge*>    mapFaceInBUToEdgeInVD;
    map<BetaEdge*, VDFace*>    mapEdgeInBUToFaceInVD;    
    map<BetaVertex*, VDCell*>  mapVertexInBUToCellInVD;


    //  make V-entities corresponding to q-entities
    makeVVerticesCorrespondingToQCells( QTIneIWDS, mapCellInBUToVertexInVD,  VD );
    makeVEdgesCorrespondingToQFaces(    QTIneIWDS, mapFaceInBUToEdgeInVD,    VD );
    makeVFacesCorrespondingToQEdges(    QTIneIWDS, mapEdgeInBUToFaceInVD,   VD );
    makeVCellsCorrespondingToQVertices( QTIneIWDS, mapVertexInBUToCellInVD, VD );

    setTopologyOfVVerticesViaQCells( mapCellInBUToVertexInVD, mapFaceInBUToEdgeInVD);
    setTopologyOfVEdgesViaQFaces(    mapFaceInBUToEdgeInVD,   mapCellInBUToVertexInVD, mapEdgeInBUToFaceInVD );
    setTopologyOfVFacesViaQEdges(    mapEdgeInBUToFaceInVD,   mapFaceInBUToEdgeInVD,   mapVertexInBUToCellInVD );
    setTopologyOfVCellsViaQVertices( mapVertexInBUToCellInVD, mapEdgeInBUToFaceInVD,   VD );

    VD.calculateBoundingBoxForVoronoiVertices();
    VD.setEndVertexOfInfiniteEdge();

}



    //  Quasi-triangulation in eIWDS --> Voronoi diagram of spheres
void GeometricConverter::makeVVerticesCorrespondingToQCells( 
                                     const BetaUniverse&         QTIneIWDS,
                                     map<BetaCell*, VDVertex*>&  mapCellInBUToVertexInVD,
                                     rg_SphereSetVoronoiDiagram& VD )
{
    const rg_dList<BetaCell>& q_cells = QTIneIWDS.getCellList();

    q_cells.reset4Loop();
    while ( q_cells.setNext4Loop() ) {
        BetaCell* currQCell = q_cells.getpEntity();

        VDVertex* currVVertex = VD.createVertex( VDVertex(currQCell->getID()) );

        currVVertex->isOnInfinity(             currQCell->isVirtual() );
        currVVertex->setPoint(                 currQCell->getMinTangentSphere().getCenter() );
        currVVertex->setRadiusOfTangentSphere( currQCell->getMinTangentSphere().getRadius() );
        currVVertex->connectBetaCell(          currQCell );

        mapCellInBUToVertexInVD.insert( make_pair( currQCell, currVVertex ) );
    }
}



void GeometricConverter::makeVEdgesCorrespondingToQFaces( 
                                     const BetaUniverse&         QTIneIWDS,
                                     map<BetaFace*, VDEdge*>&    mapFaceInBUToEdgeInVD,
                                     rg_SphereSetVoronoiDiagram& VD )
{
    const rg_dList<BetaFace>& q_faces = QTIneIWDS.getFaceList();

    q_faces.reset4Loop();
    while ( q_faces.setNext4Loop() ) {
        BetaFace* currQFace = q_faces.getpEntity();

        VDEdge* currVEdge  = VD.createEdge( VDEdge( currQFace->getID(), rg_NULL, rg_NULL ) );
        if ( currQFace->isVirtual() ) {
            currVEdge->setInfinity();
        }
        currVEdge->connectBetaFace( currQFace );

        VD.createPrEdge( currVEdge );

        mapFaceInBUToEdgeInVD.insert( make_pair(currQFace, currVEdge) );
    }
}




void GeometricConverter::makeVFacesCorrespondingToQEdges( 
                                     const BetaUniverse&         QTIneIWDS,
                                     map<BetaEdge*, VDFace*>&    mapEdgeInBUToFaceInVD,    
                                     rg_SphereSetVoronoiDiagram& VD )
{
    const rg_dList<BetaEdge>& q_edges = QTIneIWDS.getEdgeList();

    q_edges.reset4Loop();
    while ( q_edges.setNext4Loop() ) {
        BetaEdge* currQEdge = q_edges.getpEntity();

        VDFace* currVFace = VD.createFace( VDFace(currQEdge->getID()) );
        if ( currQEdge->isVirtual() ) {
            currVFace->setInfinity();
        }
        currVFace->connectBetaEdge( currQEdge );


        VDLoop* outerLoopOfCurrVFace = VD.createLoop( VDLoop(currVFace->getID(), currVFace, rg_TRUE) );
        currVFace->addLoop( outerLoopOfCurrVFace );


        rg_INT numSmallWorldInCurrQEdge = currQEdge->getNumOfSmallWorlds();
        for ( rg_INT i=0; i<numSmallWorldInCurrQEdge; i++ ) {
            VDLoop* innerLoopOfCurrVFace = VD.createLoop( VDLoop(currVFace->getID()+i+1, currVFace, rg_FALSE) );
            currVFace->addLoop( innerLoopOfCurrVFace );
        }


        mapEdgeInBUToFaceInVD.insert( make_pair( currQEdge, currVFace ) );        
    }
}




void GeometricConverter::makeVCellsCorrespondingToQVertices( 
                                     const BetaUniverse&         QTIneIWDS,
                                     map<BetaVertex*, VDCell*>&  mapVertexInBUToCellInVD,
                                     rg_SphereSetVoronoiDiagram& VD )
{
    const rg_dList<BetaVertex>& q_vertices = QTIneIWDS.getVertexList();

    q_vertices.reset4Loop();
    while ( q_vertices.setNext4Loop() ) {
        BetaVertex* currQVtx = q_vertices.getpEntity();

        VDCell* currVCell = rg_NULL;

        if ( !currQVtx->isVirtual() ) {
            Ball* currBallInBU = currQVtx->getBallProperty();

            BallGenerator* currBallInVD = VD.addBallGenerator( BallGenerator( currBallInBU->getID(), 
                                                                              currBallInBU->getGeometry(), 
                                                                              currBallInBU->getProperty(), 
                                                                              currBallInBU->getIDFromInput()) );
            currVCell = currBallInVD->getCell();
        }
        else {
            currVCell = VD.createCell( VDCell(currQVtx->getID()) );
        }

        mapVertexInBUToCellInVD.insert( make_pair( currQVtx, currVCell ) );
    }
}




void GeometricConverter::setTopologyOfVVerticesViaQCells( const map<BetaCell*, VDVertex*>&  mapCellInBUToVertexInVD,
                                                          const map<BetaFace*, VDEdge*>&    mapFaceInBUToEdgeInVD)
{
    map<BetaCell*, VDVertex*>::const_iterator i_qCell;
    for ( i_qCell = mapCellInBUToVertexInVD.begin(); i_qCell != mapCellInBUToVertexInVD.end(); i_qCell++ ) {
        BetaCell* currQCell   = i_qCell->first;
        VDVertex* currVVertex = i_qCell->second;

        for ( rg_INT i=0; i<4; i++ ) {
            BetaFace* currQFace  = currQCell->getFace( i );
            VDEdge*   currVEdge = mapFaceInBUToEdgeInVD.find( currQFace )->second;

            currVVertex->setIncidentEdge( i, currVEdge );
        }
    }
}



void GeometricConverter::setTopologyOfVEdgesViaQFaces( 
                                 const map<BetaFace*, VDEdge*>&    mapFaceInBUToEdgeInVD,
                                 const map<BetaCell*, VDVertex*>&  mapCellInBUToVertexInVD,
                                 const map<BetaEdge*, VDFace*>&    mapEdgeInBUToFaceInVD )
{
    map<BetaCell*, VDVertex*>::const_iterator i_qCell;

    map<BetaFace*, VDEdge*>::const_iterator i_qFace;
    for ( i_qFace = mapFaceInBUToEdgeInVD.begin(); i_qFace != mapFaceInBUToEdgeInVD.end(); i_qFace++ ) {
        BetaFace* currQFace = i_qFace->first;
        VDEdge*   currVEdge = i_qFace->second;

        //  VDEdge <--> VDVertex
        VDVertex* startVVtx = mapCellInBUToVertexInVD.find( currQFace->getLeftCell() )->second;
        VDVertex* endVVtx   = mapCellInBUToVertexInVD.find( currQFace->getRightCell() )->second;
        currVEdge->setStartVertex( startVVtx );
        currVEdge->setEndVertex(   endVVtx   );


        if ( currVEdge->getID() == 2951 || currVEdge->getID() == 2950 || currVEdge->getID() == 6653 || currVEdge->getID() == 7717  ){
            int stop = 1;
        }


        ////  VDEdge <--> VDFace
        //VDPartialEdge* currPrEdge[3] = { currVEdge->getPartialEdge(), 
        //                                 currPrEdge[0]->getNextPartialEdgeInRadialCycle(), 
        //                                 currPrEdge[1]->getNextPartialEdgeInRadialCycle() };
        //for ( rg_INT i=0; i<3; i++ ) {
        //    BetaEdge* currQEdge = currQFace->getEdge(i);
        //    VDFace*   currVFace = mapEdgeInBUToFaceInVD.find( currQEdge )->second;

        //    if ( currQEdge->isGateToSmallWorlds() ) continue;

        //    VDLoop*   outerLoopOfCurrVFace = currVFace->getOuterLoop();
        //    currPrEdge[i]->setLoop( outerLoopOfCurrVFace );
        //}
    }
}



void GeometricConverter::setTopologyOfVFacesViaQEdges( 
                                 const map<BetaEdge*, VDFace*>&    mapEdgeInBUToFaceInVD,
                                 const map<BetaFace*, VDEdge*>&    mapFaceInBUToEdgeInVD,
                                 const map<BetaVertex*, VDCell*>&  mapVertexInBUToCellInVD )
{
    map<BetaEdge*, VDFace*>::const_iterator i_qEdge;
    for ( i_qEdge = mapEdgeInBUToFaceInVD.begin(); i_qEdge != mapEdgeInBUToFaceInVD.end(); i_qEdge++ ) {
        BetaEdge* currQEdge = i_qEdge->first;
        VDFace*   currVFace = i_qEdge->second;


        //  VDFace <--> VDCell
        VDCell*   leftVCell  = mapVertexInBUToCellInVD.find( currQEdge->getStartVertex() )->second;
        VDCell*   rightVCell = mapVertexInBUToCellInVD.find( currQEdge->getEndVertex() )->second;
        currVFace->setLeftCell(  leftVCell );
        currVFace->setRightCell( rightVCell );

        //  VDFace <--> VDEdge
        //  set topology of outer loop of currVFace from the topology of  big world of currQEdge
        VDLoop*   outerLoopOfCurrVFace = currVFace->getOuterLoop();

        rg_dList<BetaFace*> faceList;
        currQEdge->searchFacesInIntraWorld( faceList );

        faceList.reset4Loop();
        while ( faceList.setNext4Loop() ) {
            BetaFace* currQFace = faceList.getEntity();
            VDEdge*   currVEdge = mapFaceInBUToEdgeInVD.find( currQFace )->second;

            VDPartialEdge* currPrEdge[3] = { currVEdge->getPartialEdge(), 
                                             currPrEdge[0]->getNextPartialEdgeInRadialCycle(), 
                                             currPrEdge[1]->getNextPartialEdgeInRadialCycle() };

            rg_INT         posCurrQEdge  = currQFace->getPosOfEdge(currQEdge);

            currPrEdge[ posCurrQEdge ]->setLoop( outerLoopOfCurrVFace );
            currPrEdge[ posCurrQEdge ]->isRightOrientationInLoop( currQFace->getEdgeOrientation(posCurrQEdge) );

            outerLoopOfCurrVFace->setPartialEdge( currPrEdge[ posCurrQEdge ] );
        }

        rg_dNode<BetaFace*>* currNode = faceList.getFirstpNode();
        for (rg_INT i=0; i<faceList.getSize(); i++, currNode=currNode->getNext() ) {
            BetaFace* qFace[3]        = {currNode->getPrev()->getEntity(),  currNode->getEntity(),             currNode->getNext()->getEntity()};
            rg_INT    posCurrQEdge[3] = {qFace[0]->getPosOfEdge(currQEdge), qFace[1]->getPosOfEdge(currQEdge), qFace[2]->getPosOfEdge(currQEdge)};

            VDPartialEdge* prEdge[3] = {rg_NULL, rg_NULL, rg_NULL};

            for ( rg_INT j=0; j<3; j++ ) {
                VDEdge*   currVEdge = mapFaceInBUToEdgeInVD.find( qFace[j] )->second;
                prEdge[j] = currVEdge->getPartialEdge();
                if ( posCurrQEdge[j] == 1 ) {
                    prEdge[j] = currVEdge->getPartialEdge()->getNextPartialEdgeInRadialCycle();
                }
                else if ( posCurrQEdge[j] == 2 ) {
                    prEdge[j] = currVEdge->getPartialEdge()->getNextPartialEdgeInRadialCycle()->getNextPartialEdgeInRadialCycle();
                }
            }

            prEdge[1]->setPreviousPartialEdgeInLoop( prEdge[0] );
            prEdge[1]->setNextPartialEdgeInLoop(     prEdge[2] );
        }


        //  set topology of inner loop of currVFace from the topology of  small worlds of currQEdge
        rg_INT IDOfInnerLoop = 1;
        rg_dList<BetaFace*>* smallWorlds = currQEdge->getSmallWorlds();
        smallWorlds->reset4Loop();
        while ( smallWorlds->setNext4Loop() ) {
            BetaFace* currQFaceInSmallWorld = smallWorlds->getEntity();

            VDLoop*   innerLoopOfCurrVFace = currVFace->getLoop( IDOfInnerLoop );

            rg_dList<BetaFace*> faceListInSmallWorld;
            currQEdge->searchFacesInIntraWorld( faceListInSmallWorld, currQFaceInSmallWorld );

            faceListInSmallWorld.reset4Loop();
            while ( faceListInSmallWorld.setNext4Loop() ) {
                BetaFace* currQFace = faceListInSmallWorld.getEntity();
                VDEdge*   currVEdge = mapFaceInBUToEdgeInVD.find( currQFace )->second;

                VDPartialEdge* currPrEdge[3] = { currVEdge->getPartialEdge(), 
                                                 currPrEdge[0]->getNextPartialEdgeInRadialCycle(), 
                                                 currPrEdge[1]->getNextPartialEdgeInRadialCycle() };

                rg_INT         posCurrQEdge  = currQFace->getPosOfEdge(currQEdge);

                currPrEdge[ posCurrQEdge ]->setLoop( innerLoopOfCurrVFace );
                currPrEdge[ posCurrQEdge ]->isRightOrientationInLoop( currQFace->getEdgeOrientation(posCurrQEdge) );

                innerLoopOfCurrVFace->setPartialEdge( currPrEdge[ posCurrQEdge ] );
            }

            rg_dNode<BetaFace*>* currNode = faceListInSmallWorld.getFirstpNode();
            for (rg_INT i=0; i<faceListInSmallWorld.getSize(); i++, currNode=currNode->getNext() ) {
                BetaFace* qFace[3]        = {currNode->getPrev()->getEntity(),  currNode->getEntity(),             currNode->getNext()->getEntity()};
                rg_INT    posCurrQEdge[3] = {qFace[0]->getPosOfEdge(currQEdge), qFace[1]->getPosOfEdge(currQEdge), qFace[2]->getPosOfEdge(currQEdge)};

                VDPartialEdge* prEdge[3] = {rg_NULL, rg_NULL, rg_NULL};

                for ( rg_INT j=0; j<3; j++ ) {
                    VDEdge*   currVEdge = mapFaceInBUToEdgeInVD.find( qFace[j] )->second;
                    prEdge[j] = currVEdge->getPartialEdge();
                    if ( posCurrQEdge[j] == 1 ) {
                        prEdge[j] = currVEdge->getPartialEdge()->getNextPartialEdgeInRadialCycle();
                    }
                    else if ( posCurrQEdge[j] == 2 ) {
                        prEdge[j] = currVEdge->getPartialEdge()->getNextPartialEdgeInRadialCycle()->getNextPartialEdgeInRadialCycle();
                    }
                }

                prEdge[1]->setPreviousPartialEdgeInLoop( prEdge[0] );
                prEdge[1]->setNextPartialEdgeInLoop(     prEdge[2] );
            }

        }
    }
}



void GeometricConverter::setTopologyOfVCellsViaQVertices( 
                                            const map<BetaVertex*, VDCell*>&  mapVertexInBUToCellInVD,
                                            const map<BetaEdge*, VDFace*>&    mapEdgeInBUToFaceInVD,
                                            rg_SphereSetVoronoiDiagram& VD )
{
    rg_dList< VDFace >* facesInVD = VD.getGlobalFaceList();

    facesInVD->reset4Loop();
    while ( facesInVD->setNext4Loop() ) {
        VDFace* currVFace = facesInVD->getpEntity();

        VDCell* leftVCell  = currVFace->getLeftCell();
        VDCell* rightVCell = currVFace->getRightCell();

        leftVCell->addBoundingFace( currVFace );
        rightVCell->addBoundingFace( currVFace );
    }




    VDCell* virtualCell = rg_NULL;
    rg_dList< VDCell >* cellsInVD = VD.getGlobalCellList();
    cellsInVD->reset4Loop();
    while ( cellsInVD->setNext4Loop() ) {
        VDCell* currVCell = cellsInVD->getpEntity();
        if ( currVCell->getGenerator() == rg_NULL ) {
            virtualCell = currVCell;
            break;
        }
    }


    rg_dList<VDFace*>* facesOnInfinity = virtualCell->getBoungindFaces();
    facesOnInfinity->reset4Loop();
    while ( facesOnInfinity->setNext4Loop() ) {
        VDFace* currVFace = facesOnInfinity->getEntity();

        if ( currVFace->getLeftCell() == virtualCell ) {
            currVFace->getRightCell()->isBounded( rg_FALSE );
        }
        else {
            currVFace->getLeftCell()->isBounded( rg_FALSE );
        }
    }

}



