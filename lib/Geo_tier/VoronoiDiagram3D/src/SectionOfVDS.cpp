#include "SectionOfVDS.h"
using namespace V::GeometryTier;

#include <fstream>
using namespace std;


SectionOfVDS::SectionOfVDS()
: m_VDS(rg_NULL)
{
    for(rg_INT i=0; i<4; i++)
        m_plane[i] = 0.0;
}



SectionOfVDS::~SectionOfVDS()
{
    m_VDS = rg_NULL;
}




rg_dList<VertexOnSectionOfVDS>* SectionOfVDS::getVertexList()
{
    return &m_vertexList;
}



rg_dList<EdgeOnSectionOfVDS>*   SectionOfVDS::getEdgeList()
{
    return &m_edgeList;
}



rg_dList<FaceOnSectionOfVDS>*   SectionOfVDS::getFaceList()
{
    return &m_faceList;
}




void SectionOfVDS::connectVDS( rg_SphereSetVoronoiDiagram* aVDS )
{
    m_VDS = aVDS;
}




void SectionOfVDS::setPlane(rg_REAL* plane)
{
    for(rg_INT i=0; i<4; i++)
        m_plane[i] = plane[i];    
}




void SectionOfVDS::addVertex( const VertexOnSectionOfVDS& vertex )
{
    m_vertexList.add( vertex );
}



void SectionOfVDS::addEdge( const EdgeOnSectionOfVDS& edge )
{
    m_edgeList.add( edge );
}



void SectionOfVDS::addFace( const FaceOnSectionOfVDS& face )
{
    m_faceList.add( face );
}




void SectionOfVDS::computeSectionOfVDS(rg_SphereSetVoronoiDiagram* aVDS, rg_REAL* plane)
{
    m_VDS = aVDS;
    for(rg_INT i=0; i<4; i++)
        m_plane[i] = plane[i];    

    rg_dList< GateForSectionOfVDS > gateList;
    
    initializeComputingSectionOfVDS( gateList );

    rg_dNode< GateForSectionOfVDS >* gateNode = rg_NULL;
    GateForSectionOfVDS*             currGate = rg_NULL;
    for (gateNode=gateList.getHead(); gateNode!=gateList.getTail(); gateNode=gateNode->getNext() )  {
        currGate = gateNode->getpEntity();
        if ( currGate->getStartVertex() == rg_NULL )
            continue;

        if ( currGate->m_ID == 123 )
            int aaa = 1;

        rg_FLAG vertexDeterminant = rg_TRUE;
        VDEdge* geoForNewVertex = findNewVertex( currGate, vertexDeterminant );

        exploreNewVertex( geoForNewVertex, vertexDeterminant, currGate, gateList);  
        
        if ( currGate->getEdge()->getID() == 192 )  // gate 123
            int aaa = 1;

        if ( m_faceList.getLastpEntity()->getID() == 42 )
            int bbb = 2;
    }
    adjustFacesIncidentToInfiniteVertex();
    computeVertexCoordnates();
//    checkTopology();

}



void SectionOfVDS::initializeComputingSectionOfVDS(rg_dList< GateForSectionOfVDS >& gateList)
{
    rg_Vector3D normal(m_plane[0], m_plane[1], m_plane[2]);
    rg_dList< VDEdge >* globalEdges = m_VDS->getGlobalEdgeList();

    VertexOnSectionOfVDS* ptrVertex = m_vertexList.add( VertexOnSectionOfVDS() );
    pair<VDEdge*, rg_Point3D> geoVertex( (VDEdge*)rg_NULL, rg_Point3D() );
    ptrVertex->setGeometry( geoVertex );
    ptrVertex->setID( 0 );

    
    VDEdge* initialEdge = rg_NULL;
    VDEdge* currEdge    = rg_NULL;
   	globalEdges->reset4Loop();
    while ( globalEdges->setNext4Loop() )  {
        currEdge = globalEdges->getpEntity();

        if ( currEdge->isOnInfinity() )
            continue;

        double distPAndStart = normal.innerProduct(currEdge->getStartVertex()->getPoint()) + m_plane[3];
        double distPAndEnd   = normal.innerProduct(currEdge->getEndVertex()->getPoint()) + m_plane[3];
        double determinant = distPAndStart*distPAndEnd;

        if ( rg_LE( determinant, 0.0 ) )  {
            initialEdge = currEdge;

            VertexOnSectionOfVDS* ptrVertex = m_vertexList.add( VertexOnSectionOfVDS() );
            pair<VDEdge*, rg_Point3D> geoVertex( initialEdge, rg_Point3D() );
            ptrVertex->setGeometry( geoVertex );
            ptrVertex->setID( m_vertexList.getSize()-1 );

            
            GateForSectionOfVDS* ptrGates[3];
            rg_INT               i_gate = 0;

            if ( distPAndEnd >= distPAndStart )  {
                VDPartialEdge* startPrEdge = initialEdge->getPartialEdge();
                VDPartialEdge* currPrEdge  = startPrEdge;

                do  {
                    VDFace* currFace = currPrEdge->getLoop()->getFace();
                    if ( currPrEdge->isRightOrientationInLoop() )  {
                        ptrGates[i_gate] = gateList.add( GateForSectionOfVDS( currFace->getLeftCell(), 
                                                                                currFace->getRightCell(), 
                                                                                ptrVertex, currPrEdge) );
                    }
                    else  {
                        ptrGates[i_gate] = gateList.add( GateForSectionOfVDS( currFace->getRightCell(), 
                                                                                currFace->getLeftCell(), 
                                                                                ptrVertex, currPrEdge) );
                    }
                    i_gate++;
                    currPrEdge = currPrEdge->getNextPartialEdgeInRadialCycle();
                } while ( currPrEdge != startPrEdge );

            }
            else  {
                i_gate = 2;
                VDPartialEdge* startPrEdge = initialEdge->getPartialEdge();
                VDPartialEdge* currPrEdge  = startPrEdge;

                do  {
                    VDFace* currFace = currPrEdge->getLoop()->getFace();
                    if ( currPrEdge->isRightOrientationInLoop() )  {
                        ptrGates[i_gate] = gateList.add( GateForSectionOfVDS( currFace->getRightCell(), 
                                                                                currFace->getLeftCell(), 
                                                                                ptrVertex, currPrEdge) );
                    }
                    else  {
                        ptrGates[i_gate] = gateList.add( GateForSectionOfVDS( currFace->getLeftCell(), 
                                                                                currFace->getRightCell(), 
                                                                                ptrVertex, currPrEdge) );
                    }
                    i_gate--;
                    currPrEdge = currPrEdge->getNextPartialEdgeInRadialCycle();
                } while ( currPrEdge != startPrEdge );
            }
            ptrGates[0]->setGateForLeftLeg( ptrGates[2] );
            ptrGates[0]->setGateForRightLeg( ptrGates[1] );

            ptrGates[1]->setGateForLeftLeg( ptrGates[0] );
            ptrGates[1]->setGateForRightLeg( ptrGates[2] );

            ptrGates[2]->setGateForLeftLeg( ptrGates[1] );
            ptrGates[2]->setGateForRightLeg( ptrGates[0] );
            /*
            ptrGates[0]->setGateForLeftLeg( ptrGates[1] );
            ptrGates[0]->setGateForRightLeg( ptrGates[2] );

            ptrGates[1]->setGateForLeftLeg( ptrGates[2] );
            ptrGates[1]->setGateForRightLeg( ptrGates[0] );

            ptrGates[2]->setGateForLeftLeg( ptrGates[0] );
            ptrGates[2]->setGateForRightLeg( ptrGates[1] );
            */
            for ( rg_INT i=0; i<3; i++ )
                ptrGates[i]->m_ID = i;

            break;
        }
    }   
}




VDEdge* SectionOfVDS::findNewVertex(GateForSectionOfVDS* gate, rg_FLAG& vertexDeterminant)
{
    VDEdge*        usedEdge    = gate->getGeometryForEdge()->getOriginalEdge();
    VDPartialEdge* startPrEdge = gate->getGeometryForEdge();
    VDPartialEdge* currPrEdge  = startPrEdge;

    do  {
        
        VDEdge* currEdge = currPrEdge->getOriginalEdge();

        if ( currEdge->isOnInfinity() || usedEdge == currEdge )  {
            currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
            continue;
        }

        rg_Point3D start = currEdge->getStartVertex()->getPoint();
        rg_Point3D end   = currEdge->getEndVertex()->getPoint();

        double distPAndStart = (m_plane[0]*start.getX()) + (m_plane[1]*start.getY()) + (m_plane[2]*start.getZ()) + m_plane[3];
        double distPAndEnd   = (m_plane[0]*end.getX()) + (m_plane[1]*end.getY()) + (m_plane[2]*end.getZ()) + m_plane[3];
        double determinant = distPAndStart*distPAndEnd;

        if ( rg_LE( determinant, 0.0 ) )  {
            if ( distPAndEnd >= distPAndStart )  
                vertexDeterminant = rg_TRUE;
            else 
                vertexDeterminant = rg_FALSE;

            return currEdge;
        }

        currPrEdge = currPrEdge->getNextPartialEdgeInLoop();
    } while ( startPrEdge != currPrEdge );

    return rg_NULL;
}



rg_FLAG SectionOfVDS::exploreNewVertex( VDEdge*                          geoForNewVertex,
                                        const rg_FLAG&                   vertexDeterminant,
                                        GateForSectionOfVDS*             gate,
                                        rg_dList< GateForSectionOfVDS >& gateList)
{
    //  unbounded edge.
    if ( geoForNewVertex == rg_NULL )  {
        VertexOnSectionOfVDS* ptrVertex = m_vertexList.getFirstpEntity();
        
        //  2. define current edge
        EdgeOnSectionOfVDS* ptrEdge = m_edgeList.add( EdgeOnSectionOfVDS() );  
        ptrEdge->setGeometry( gate->getGeometryForEdge()->getLoop() );
        ptrEdge->setID( m_edgeList.getSize()-1 );
        gate->setEdge( ptrEdge );

            //  2.1 set start and end vertices of current edge.
            ptrEdge->setStartVertex( gate->getStartVertex() );
            ptrEdge->setEndVertex( ptrVertex );

            gate->getStartVertex()->setFirstEdge( ptrEdge );
            ptrVertex->setFirstEdge( ptrEdge );

            //  2.2 set left leg of current edge.
            EdgeOnSectionOfVDS* leftLeg = gate->getGateForLeftLeg()->getEdge();
            if ( leftLeg != rg_NULL )  {
                ptrEdge->setLeftLeg( leftLeg );

                if ( ptrEdge->getStartVertex() == leftLeg->getEndVertex() )  {
                    leftLeg->setLeftHand( ptrEdge );
                    ptrEdge->setLeftFace( leftLeg->getLeftFace() );
                }
                else  {
                    leftLeg->setRightLeg( ptrEdge );
                    ptrEdge->setLeftFace( leftLeg->getRightFace() );
                }
            }

            //  2.3 set right leg of current edge.
            EdgeOnSectionOfVDS* rightLeg = gate->getGateForRightLeg()->getEdge();
            if ( rightLeg != rg_NULL )  {
                ptrEdge->setRightLeg( rightLeg );

                if ( ptrEdge->getStartVertex() == rightLeg->getEndVertex() )  {
                    rightLeg->setRightHand( ptrEdge );
                    ptrEdge->setRightFace( rightLeg->getRightFace() );
                }
                else  {
                    rightLeg->setLeftLeg( ptrEdge );
                    ptrEdge->setRightFace( rightLeg->getLeftFace() );
                }
            }
            //  2.4 set left hand leg of current edge.
            ptrEdge->setLeftHand( rg_NULL );
            //  2.5 set right hand leg of current edge.
            ptrEdge->setRightHand( rg_NULL );


            //  2.6 set left face of current edge
            if ( ptrEdge->getLeftFace() == rg_NULL )  {
                FaceOnSectionOfVDS* ptrLeftFace = m_faceList.add( FaceOnSectionOfVDS() );  
                ptrLeftFace->setGeometry( gate->getGeometryForLeftFace() );
                ptrLeftFace->setFirstEdge( ptrEdge );
                ptrLeftFace->setID( m_faceList.getSize()-1 );

                ptrEdge->setLeftFace( ptrLeftFace );
            }

            //  2.7 set right face of current edge
            if ( ptrEdge->getRightFace() == rg_NULL )  {
                FaceOnSectionOfVDS* ptrRightFace = m_faceList.add( FaceOnSectionOfVDS() );  
                ptrRightFace->setGeometry( gate->getGeometryForRightFace() );
                ptrRightFace->setFirstEdge( ptrEdge );
                ptrRightFace->setID( m_faceList.getSize()-1 );

                ptrEdge->setRightFace( ptrRightFace );
            }

            
        return rg_TRUE;
    }



    //  bounded edge.
    //  chech whether new vertex is predefined or not.
    VertexOnSectionOfVDS* predefinedVertex = rg_NULL;
    VertexOnSectionOfVDS* currVertex       = rg_NULL;
    m_vertexList.reset4Loop();
    while ( m_vertexList.setNext4Loop() )  {
        currVertex = m_vertexList.getpEntity();

        if ( currVertex->getpGeometry()->first == geoForNewVertex )  {
            predefinedVertex = currVertex;
            break;
        }
    }



    if ( predefinedVertex == rg_NULL )  {
    //  new vertex is really NEW!.
    //  make 2 gates and add to gateList

        //  1. define end vertex.
        VertexOnSectionOfVDS* ptrVertex = m_vertexList.add( VertexOnSectionOfVDS() );       
        ptrVertex->setGeometry( pair<VDEdge*, rg_Point3D>( geoForNewVertex, rg_Point3D() ) );
        ptrVertex->setID( m_vertexList.getSize()-1 );


        //  2. define current edge
        EdgeOnSectionOfVDS* ptrEdge = m_edgeList.add( EdgeOnSectionOfVDS() );  
        ptrEdge->setGeometry( gate->getGeometryForEdge()->getLoop() );
        ptrEdge->setID( m_edgeList.getSize()-1 );
        gate->setEdge( ptrEdge );

            //  2.1 set start and end vertices of current edge.
            ptrEdge->setStartVertex( gate->getStartVertex() );
            ptrEdge->setEndVertex( ptrVertex );

            gate->getStartVertex()->setFirstEdge( ptrEdge );
            ptrVertex->setFirstEdge( ptrEdge );

            //  2.2 set left leg of current edge.
            EdgeOnSectionOfVDS* leftLeg = gate->getGateForLeftLeg()->getEdge();
            if ( leftLeg != rg_NULL )  {
                ptrEdge->setLeftLeg( leftLeg );

                if ( ptrEdge->getStartVertex() == leftLeg->getEndVertex() )  {
                    leftLeg->setLeftHand( ptrEdge );
                    ptrEdge->setLeftFace( leftLeg->getLeftFace() );
                }
                else  {
                    leftLeg->setRightLeg( ptrEdge );
                    ptrEdge->setLeftFace( leftLeg->getRightFace() );
                }
            }

            //  2.3 set right leg of current edge.
            EdgeOnSectionOfVDS* rightLeg = gate->getGateForRightLeg()->getEdge();
            if ( rightLeg != rg_NULL )  {
                ptrEdge->setRightLeg( rightLeg );

                if ( ptrEdge->getStartVertex() == rightLeg->getEndVertex() )  {
                    rightLeg->setRightHand( ptrEdge );
                    ptrEdge->setRightFace( rightLeg->getRightFace() );
                }
                else  {
                    rightLeg->setLeftLeg( ptrEdge );
                    ptrEdge->setRightFace( rightLeg->getLeftFace() );
                }
            }

            //  2.4 set left face of current edge
            if ( ptrEdge->getLeftFace() == rg_NULL )  {
                FaceOnSectionOfVDS* ptrLeftFace = m_faceList.add( FaceOnSectionOfVDS() );  
                ptrLeftFace->setGeometry( gate->getGeometryForLeftFace() );
                ptrLeftFace->setFirstEdge( ptrEdge );
                ptrLeftFace->setID( m_faceList.getSize()-1 );

                ptrEdge->setLeftFace( ptrLeftFace );
            }

            //  2.5 set right face of current edge
            if ( ptrEdge->getRightFace() == rg_NULL )  {
                FaceOnSectionOfVDS* ptrRightFace = m_faceList.add( FaceOnSectionOfVDS() );  
                ptrRightFace->setGeometry( gate->getGeometryForRightFace() );
                ptrRightFace->setFirstEdge( ptrEdge );
                ptrRightFace->setID( m_faceList.getSize()-1 );

                ptrEdge->setRightFace( ptrRightFace );
            }


        // 3. make 2 gates
        rg_INT i_gate = 0;
        GateForSectionOfVDS  newGates[3];
        GateForSectionOfVDS* ptrNewGates[3] = {rg_NULL, rg_NULL, rg_NULL};
        rg_INT               validityNewGates[3] = {rg_TRUE, rg_TRUE, rg_TRUE};
        if ( vertexDeterminant )  {
            VDPartialEdge* startPrEdge = geoForNewVertex->getPartialEdge();
            VDPartialEdge* currPrEdge  = startPrEdge;

            do  {
                if ( gate->getGeometryForEdge()->getLoop() == currPrEdge->getLoop() )  {
                    validityNewGates[i_gate] = rg_FALSE;
                }
                else  {                
                    VDFace* currFace = currPrEdge->getLoop()->getFace();
                    if ( currPrEdge->isRightOrientationInLoop() )  {
                        newGates[i_gate].setGeometryForLeftFace( currFace->getLeftCell() ); 
                        newGates[i_gate].setGeometryForRightFace( currFace->getRightCell() ); 
                        newGates[i_gate].setStartVertex( ptrVertex );
                        newGates[i_gate].setGeometryForEdge( currPrEdge );
                    }
                    else  {
                        newGates[i_gate].setGeometryForLeftFace( currFace->getRightCell() ); 
                        newGates[i_gate].setGeometryForRightFace( currFace->getLeftCell() ); 
                        newGates[i_gate].setStartVertex( ptrVertex );
                        newGates[i_gate].setGeometryForEdge( currPrEdge );
                    }
                }
                i_gate++;

                currPrEdge = currPrEdge->getNextPartialEdgeInRadialCycle();
            } while ( currPrEdge != startPrEdge );
            
        }
        else  {
            i_gate = 2;
            VDPartialEdge* startPrEdge = geoForNewVertex->getPartialEdge();
            VDPartialEdge* currPrEdge  = startPrEdge;

            do  {
                if ( gate->getGeometryForEdge()->getLoop() == currPrEdge->getLoop() )  {
                    validityNewGates[i_gate] = rg_FALSE;
                }
                else  {                
                    VDFace* currFace = currPrEdge->getLoop()->getFace();
                    if ( currPrEdge->isRightOrientationInLoop() )  {
                        newGates[i_gate].setGeometryForLeftFace( currFace->getRightCell() ); 
                        newGates[i_gate].setGeometryForRightFace( currFace->getLeftCell() ); 
                        newGates[i_gate].setStartVertex( ptrVertex );
                        newGates[i_gate].setGeometryForEdge( currPrEdge );
                    }
                    else  {
                        newGates[i_gate].setGeometryForLeftFace( currFace->getLeftCell() ); 
                        newGates[i_gate].setGeometryForRightFace( currFace->getRightCell() ); 
                        newGates[i_gate].setStartVertex( ptrVertex );
                        newGates[i_gate].setGeometryForEdge( currPrEdge );
                    }
                }
                i_gate--;

                currPrEdge = currPrEdge->getNextPartialEdgeInRadialCycle();
            } while ( currPrEdge != startPrEdge );            
        }

        for ( rg_INT i=0; i<3; i++ )  {
            if ( validityNewGates[i] )  {
                ptrNewGates[i] = gateList.add( newGates[i] );
            }
        }

        if ( validityNewGates[0] == rg_FALSE )  {
            gate->setGateForLeftHand( ptrNewGates[1] );
            gate->setGateForRightHand( ptrNewGates[2] );
            ptrNewGates[1]->setGateForLeftLeg(gate);
            ptrNewGates[1]->setGateForRightLeg(ptrNewGates[2]);
            ptrNewGates[2]->setGateForLeftLeg(ptrNewGates[1]);
            ptrNewGates[2]->setGateForRightLeg(gate);

            ptrNewGates[1]->m_ID = gateList.getSize()-2;
            ptrNewGates[2]->m_ID = gateList.getSize()-1;
        }
        else if ( validityNewGates[1] == rg_FALSE )  {
            gate->setGateForLeftHand( ptrNewGates[2] );
            gate->setGateForRightHand( ptrNewGates[0] );
            ptrNewGates[0]->setGateForLeftLeg(ptrNewGates[2]);
            ptrNewGates[0]->setGateForRightLeg(gate);
            ptrNewGates[2]->setGateForLeftLeg(gate);
            ptrNewGates[2]->setGateForRightLeg(ptrNewGates[0]);

            ptrNewGates[0]->m_ID = gateList.getSize()-2;
            ptrNewGates[2]->m_ID = gateList.getSize()-1;
        }
        else  {
            gate->setGateForLeftHand( ptrNewGates[0] );
            gate->setGateForRightHand( ptrNewGates[1] );
            ptrNewGates[0]->setGateForLeftLeg(gate);
            ptrNewGates[0]->setGateForRightLeg(ptrNewGates[1]);
            ptrNewGates[1]->setGateForLeftLeg(ptrNewGates[0]);
            ptrNewGates[1]->setGateForRightLeg(gate);

            ptrNewGates[0]->m_ID = gateList.getSize()-2;
            ptrNewGates[1]->m_ID = gateList.getSize()-1;
        }

        return rg_TRUE;
    }
    //  there is the predefined vertex.
    else  {
        VDLoop* currLoop = gate->getGeometryForEdge()->getLoop();

        EdgeOnSectionOfVDS* predefinedEdge = rg_NULL;
        EdgeOnSectionOfVDS* currEdge       = rg_NULL;
        m_edgeList.reset4Loop();
        while ( m_edgeList.setNext4Loop() )  {
            currEdge = m_edgeList.getpEntity();

            if ( *(currEdge->getpGeometry()) == currLoop )  {
                predefinedEdge = currEdge;
                break;
            }
        }
        
        //  Even though current vertex is predefined, 
        //  edge is not defined.
        if ( predefinedEdge == rg_NULL )  {

            GateForSectionOfVDS* predefinedGate = rg_NULL;
            GateForSectionOfVDS* currGate = rg_NULL;
            gateList.reset4Loop();
            while ( gateList.setNext4Loop() )  {
                currGate = gateList.getpEntity();

                if ( currGate == gate )
                    continue;

                if ( currLoop == currGate->getGeometryForEdge()->getLoop() )  {
                    predefinedGate = currGate;
                    break;
                }
            }

            VertexOnSectionOfVDS* ptrVertex = predefinedGate->getStartVertex();  
            
            //  1. define current edge
            EdgeOnSectionOfVDS* ptrEdge = m_edgeList.add( EdgeOnSectionOfVDS() );  
            ptrEdge->setGeometry( currLoop );
            ptrEdge->setID( m_edgeList.getSize()-1 );
            gate->setEdge( ptrEdge );

            //  1.1 set start and end vertices of current edge.
            ptrEdge->setStartVertex( gate->getStartVertex() );
            ptrEdge->setEndVertex( ptrVertex );

            gate->getStartVertex()->setFirstEdge( ptrEdge );
            ptrVertex->setFirstEdge( ptrEdge );
            predefinedGate->setStartVertex( rg_NULL );

            gate->setGateForLeftHand( predefinedGate->getGateForRightLeg() );
            gate->setGateForRightHand( predefinedGate->getGateForLeftLeg() );

            if ( ptrVertex == predefinedGate->getGateForLeftLeg()->getStartVertex() )  {            
                predefinedGate->getGateForRightLeg()->setGateForRightHand( gate );
                predefinedGate->getGateForLeftLeg()->setGateForRightLeg( gate );
            }
            else if ( ptrVertex == predefinedGate->getGateForRightLeg()->getStartVertex() )  {
                predefinedGate->getGateForRightLeg()->setGateForLeftLeg( gate );
                predefinedGate->getGateForLeftLeg()->setGateForLeftHand( gate );
            }
            else;

            //  1.2 set left leg of current edge.
            EdgeOnSectionOfVDS* leftLeg = gate->getGateForLeftLeg()->getEdge();
            if ( leftLeg != rg_NULL )  {
                ptrEdge->setLeftLeg( leftLeg );

                if ( ptrEdge->getStartVertex() == leftLeg->getEndVertex() )  {
                    leftLeg->setLeftHand( ptrEdge );
                    ptrEdge->setLeftFace( leftLeg->getLeftFace() );
                }
                else  {
                    leftLeg->setRightLeg( ptrEdge );
                    ptrEdge->setLeftFace( leftLeg->getRightFace() );
                }
            }

            //  1.3 set right leg of current edge.
            EdgeOnSectionOfVDS* rightLeg = gate->getGateForRightLeg()->getEdge();
            if ( rightLeg != rg_NULL )  {
                ptrEdge->setRightLeg( rightLeg );

                if ( ptrEdge->getStartVertex() == rightLeg->getEndVertex() )  {
                    rightLeg->setRightHand( ptrEdge );
                    ptrEdge->setRightFace( rightLeg->getRightFace() );
                }
                else  {
                    rightLeg->setLeftLeg( ptrEdge );
                    ptrEdge->setRightFace( rightLeg->getLeftFace() );
                }
            }

            //  1.4 set left hand leg of current edge.
            EdgeOnSectionOfVDS* leftHand = gate->getGateForLeftHand()->getEdge();
            if ( leftHand != rg_NULL )  {
                ptrEdge->setLeftHand( leftHand );

                if ( ptrEdge->getEndVertex() == leftHand->getStartVertex() )  {
                    leftHand->setLeftLeg( ptrEdge );
                    ptrEdge->setLeftFace( leftHand->getLeftFace() );
                }
                else  {
                    leftHand->setRightHand( ptrEdge );
                    ptrEdge->setLeftFace( leftHand->getRightFace() );
                }
            }

            //  1.5 set right hand leg of current edge.
            EdgeOnSectionOfVDS* rightHand = gate->getGateForRightHand()->getEdge();
            if ( rightHand != rg_NULL )  {
                ptrEdge->setRightHand( rightHand );

                if ( ptrEdge->getEndVertex() == rightHand->getStartVertex() )  {
                    rightHand->setRightLeg( ptrEdge );
                    ptrEdge->setRightFace( rightHand->getRightFace() );
                }
                else  {
                    rightHand->setLeftHand( ptrEdge );
                    ptrEdge->setRightFace( rightHand->getLeftFace() );
                }
            }

            //  1.6 set left face of current edge
            if ( ptrEdge->getLeftFace() == rg_NULL )  {
                FaceOnSectionOfVDS* ptrLeftFace = m_faceList.add( FaceOnSectionOfVDS() );  
                ptrLeftFace->setGeometry( gate->getGeometryForLeftFace() );
                ptrLeftFace->setFirstEdge( ptrEdge );
                ptrLeftFace->setID( m_faceList.getSize()-1 );

                ptrEdge->setLeftFace( ptrLeftFace );
            }

            //  1.7 set right face of current edge
            if ( ptrEdge->getRightFace() == rg_NULL )  {
                FaceOnSectionOfVDS* ptrRightFace = m_faceList.add( FaceOnSectionOfVDS() );  
                ptrRightFace->setGeometry( gate->getGeometryForRightFace() );
                ptrRightFace->setFirstEdge( ptrEdge );
                ptrRightFace->setID( m_faceList.getSize()-1 );

                ptrEdge->setRightFace( ptrRightFace );
            }

        }

        return rg_TRUE;
    }

    return rg_TRUE;
}




void SectionOfVDS::computeVertexCoordnates()
{
    int countInfiniteVertex = 0;
    VertexOnSectionOfVDS* ptrCurrVertex = rg_NULL;
    m_vertexList.reset4Loop();
    while ( m_vertexList.setNext4Loop() )  {
        ptrCurrVertex = m_vertexList.getpEntity();

        VDEdge* ptrEdgeVD = ptrCurrVertex->getpGeometry()->first;

        if ( ptrEdgeVD == rg_NULL )  {
            countInfiniteVertex++;
            continue;
        }
        rg_Point3D start = ptrEdgeVD->getStartVertex()->getPoint();
        rg_Point3D end   = ptrEdgeVD->getEndVertex()->getPoint();

        rg_REAL distPAndStart = (m_plane[0]*start.getX()) + (m_plane[1]*start.getY()) + (m_plane[2]*start.getZ()) + m_plane[3];
        rg_REAL distPAndEnd   = (m_plane[0]*end.getX()) + (m_plane[1]*end.getY()) + (m_plane[2]*end.getZ()) + m_plane[3];

        rg_REAL ratio = rg_ABS(distPAndStart)/rg_ABS(distPAndEnd-distPAndStart);
        rg_Point3D pt = (end*ratio) + (start*(1-ratio));


        ptrCurrVertex->setGeometry( pair<VDEdge*, rg_Point3D>(ptrEdgeVD, pt) );
    }
}



void SectionOfVDS::checkTopology()
{
    ofstream fout("test-1yat.txt");
    EdgeOnSectionOfVDS* currEdge = rg_NULL;
    m_edgeList.reset4Loop();
    while ( m_edgeList.setNext4Loop() )  {
        currEdge = m_edgeList.getpEntity();
        
        if (    currEdge->getStartVertex()->getID() == 0 
             || currEdge->getEndVertex()->getID() == 0 )
             continue;

            fout << currEdge->getID() << "\t"
                 << currEdge->getStartVertex()->getID() << "\t"
                 << currEdge->getEndVertex()->getID()   << "\t"
                 << currEdge->getLeftFace()->getID()    << "\t"
                 << currEdge->getRightFace()->getID()   << "\t"
                 << currEdge->getLeftHand()->getID()    << "\t"
                 << currEdge->getRightHand()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getID()     << "\t"
                 << currEdge->getRightLeg()->getID()    << "\t"
                 << currEdge->getGeometry()->getID() << "\t"
                 << currEdge->getStartVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getEndVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getLeftFace()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightFace()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftHand()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightHand()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getGeometry()->getID()     << "\t"
                 << currEdge->getRightLeg()->getGeometry()->getID()    << endl;
    }    
    
    /*
    ofstream fout("test-1yat.txt");

    rg_INT inValidFaceID[7] = { 42, 43, 44, 49, 67, 72, 80 };

    EdgeOnSectionOfVDS* currEdge = rg_NULL;
    m_edgeList.reset4Loop();
    while ( m_edgeList.setNext4Loop() )  {
        currEdge = m_edgeList.getpEntity();

        if (    currEdge->getLeftFace()->getID() == inValidFaceID[0]
            || currEdge->getRightFace()->getID() == inValidFaceID[0] )  
        {
            fout << currEdge->getID() << "\t"
                 << currEdge->getStartVertex()->getID() << "\t"
                 << currEdge->getEndVertex()->getID()   << "\t"
                 << currEdge->getLeftFace()->getID()    << "\t"
                 << currEdge->getRightFace()->getID()   << "\t"
                 << currEdge->getLeftHand()->getID()    << "\t"
                 << currEdge->getRightHand()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getID()     << "\t"
                 << currEdge->getRightLeg()->getID()    << "\t"
                 << currEdge->getGeometry()->getID() << "\t"
                 << currEdge->getStartVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getEndVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getLeftFace()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightFace()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftHand()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightHand()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getGeometry()->getID()     << "\t"
                 << currEdge->getRightLeg()->getGeometry()->getID()    << endl;
        }
    }    

    fout << endl;
    m_edgeList.reset4Loop();
    while ( m_edgeList.setNext4Loop() )  {
        currEdge = m_edgeList.getpEntity();

        if (    currEdge->getLeftFace()->getID() == inValidFaceID[1]
            || currEdge->getRightFace()->getID() == inValidFaceID[1] )  
        {
            fout << currEdge->getID() << "\t"
                 << currEdge->getStartVertex()->getID() << "\t"
                 << currEdge->getEndVertex()->getID()   << "\t"
                 << currEdge->getLeftFace()->getID()    << "\t"
                 << currEdge->getRightFace()->getID()   << "\t"
                 << currEdge->getLeftHand()->getID()    << "\t"
                 << currEdge->getRightHand()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getID()     << "\t"
                 << currEdge->getRightLeg()->getID()    << "\t"
                 << currEdge->getGeometry()->getID() << "\t"
                 << currEdge->getStartVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getEndVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getLeftFace()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightFace()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftHand()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightHand()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getGeometry()->getID()     << "\t"
                 << currEdge->getRightLeg()->getGeometry()->getID()    << endl;
        }
    }    
    fout << endl;
    m_edgeList.reset4Loop();
    while ( m_edgeList.setNext4Loop() )  {
        currEdge = m_edgeList.getpEntity();

        if (    currEdge->getLeftFace()->getID() == inValidFaceID[2]
            || currEdge->getRightFace()->getID() == inValidFaceID[2] )  
        {
            fout << currEdge->getID() << "\t"
                 << currEdge->getStartVertex()->getID() << "\t"
                 << currEdge->getEndVertex()->getID()   << "\t"
                 << currEdge->getLeftFace()->getID()    << "\t"
                 << currEdge->getRightFace()->getID()   << "\t"
                 << currEdge->getLeftHand()->getID()    << "\t"
                 << currEdge->getRightHand()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getID()     << "\t"
                 << currEdge->getRightLeg()->getID()    << "\t"
                 << currEdge->getGeometry()->getID() << "\t"
                 << currEdge->getStartVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getEndVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getLeftFace()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightFace()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftHand()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightHand()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getGeometry()->getID()     << "\t"
                 << currEdge->getRightLeg()->getGeometry()->getID()    << endl;
        }
    }    
    fout << endl;
    m_edgeList.reset4Loop();
    while ( m_edgeList.setNext4Loop() )  {
        currEdge = m_edgeList.getpEntity();

        if (    currEdge->getLeftFace()->getID() == inValidFaceID[3]
            || currEdge->getRightFace()->getID() == inValidFaceID[3] )  
        {
            fout << currEdge->getID() << "\t"
                 << currEdge->getStartVertex()->getID() << "\t"
                 << currEdge->getEndVertex()->getID()   << "\t"
                 << currEdge->getLeftFace()->getID()    << "\t"
                 << currEdge->getRightFace()->getID()   << "\t"
                 << currEdge->getLeftHand()->getID()    << "\t"
                 << currEdge->getRightHand()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getID()     << "\t"
                 << currEdge->getRightLeg()->getID()    << "\t"
                 << currEdge->getGeometry()->getID() << "\t"
                 << currEdge->getStartVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getEndVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getLeftFace()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightFace()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftHand()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightHand()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getGeometry()->getID()     << "\t"
                 << currEdge->getRightLeg()->getGeometry()->getID()    << endl;
        }
    }    
    fout << endl;
    m_edgeList.reset4Loop();
    while ( m_edgeList.setNext4Loop() )  {
        currEdge = m_edgeList.getpEntity();

        if (    currEdge->getLeftFace()->getID() == inValidFaceID[4]
            || currEdge->getRightFace()->getID() == inValidFaceID[4] )  
        {
            fout << currEdge->getID() << "\t"
                 << currEdge->getStartVertex()->getID() << "\t"
                 << currEdge->getEndVertex()->getID()   << "\t"
                 << currEdge->getLeftFace()->getID()    << "\t"
                 << currEdge->getRightFace()->getID()   << "\t"
                 << currEdge->getLeftHand()->getID()    << "\t"
                 << currEdge->getRightHand()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getID()     << "\t"
                 << currEdge->getRightLeg()->getID()    << "\t"
                 << currEdge->getGeometry()->getID() << "\t"
                 << currEdge->getStartVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getEndVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getLeftFace()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightFace()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftHand()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightHand()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getGeometry()->getID()     << "\t"
                 << currEdge->getRightLeg()->getGeometry()->getID()    << endl;
        }
    }    
    fout << endl;
    m_edgeList.reset4Loop();
    while ( m_edgeList.setNext4Loop() )  {
        currEdge = m_edgeList.getpEntity();

        if (    currEdge->getLeftFace()->getID() == inValidFaceID[5]
            || currEdge->getRightFace()->getID() == inValidFaceID[5] )  
        {
            fout << currEdge->getID() << "\t"
                 << currEdge->getStartVertex()->getID() << "\t"
                 << currEdge->getEndVertex()->getID()   << "\t"
                 << currEdge->getLeftFace()->getID()    << "\t"
                 << currEdge->getRightFace()->getID()   << "\t"
                 << currEdge->getLeftHand()->getID()    << "\t"
                 << currEdge->getRightHand()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getID()     << "\t"
                 << currEdge->getRightLeg()->getID()    << "\t"
                 << currEdge->getGeometry()->getID() << "\t"
                 << currEdge->getStartVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getEndVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getLeftFace()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightFace()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftHand()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightHand()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getGeometry()->getID()     << "\t"
                 << currEdge->getRightLeg()->getGeometry()->getID()    << endl;
        }
    }    
    fout << endl;
    m_edgeList.reset4Loop();
    while ( m_edgeList.setNext4Loop() )  {
        currEdge = m_edgeList.getpEntity();

        if (    currEdge->getLeftFace()->getID() == inValidFaceID[6]
            || currEdge->getRightFace()->getID() == inValidFaceID[6] )  
        {
            fout << currEdge->getID() << "\t"
                 << currEdge->getStartVertex()->getID() << "\t"
                 << currEdge->getEndVertex()->getID()   << "\t"
                 << currEdge->getLeftFace()->getID()    << "\t"
                 << currEdge->getRightFace()->getID()   << "\t"
                 << currEdge->getLeftHand()->getID()    << "\t"
                 << currEdge->getRightHand()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getID()     << "\t"
                 << currEdge->getRightLeg()->getID()    << "\t"
                 << currEdge->getGeometry()->getID() << "\t"
                 << currEdge->getStartVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getEndVertex()->getGeometry().first->getID() << "\t" 
                 << currEdge->getLeftFace()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightFace()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftHand()->getGeometry()->getID()    << "\t"
                 << currEdge->getRightHand()->getGeometry()->getID()   << "\t"
                 << currEdge->getLeftLeg()->getGeometry()->getID()     << "\t"
                 << currEdge->getRightLeg()->getGeometry()->getID()    << endl;
        }
    }    


    fout << endl << endl;
    VertexOnSectionOfVDS* currVertex = rg_NULL;
    m_vertexList.reset4Loop();
    while ( m_vertexList.setNext4Loop() )  {
        currVertex = m_vertexList.getpEntity();

        if (    currVertex->getID() == 59
             || currVertex->getID() == 62
             || currVertex->getID() == 63
             || currVertex->getID() == 71
             || currVertex->getID() == 74
           )
        {
            rg_dList< EdgeOnSectionOfVDS* > edgeList;
            currVertex->getIncidentWEdgeList( edgeList );
            edgeList.reset4Loop();
            while ( edgeList.setNext4Loop() )  {
                currEdge = edgeList.getEntity();

                fout << currEdge->getID() << "\t"
                     << currEdge->getStartVertex()->getID() << "\t"
                     << currEdge->getEndVertex()->getID()   << "\t"
                     << currEdge->getLeftFace()->getID()    << "\t"
                     << currEdge->getRightFace()->getID()   << "\t"
                     << currEdge->getLeftHand()->getID()    << "\t"
                     << currEdge->getRightHand()->getID()   << "\t"
                     << currEdge->getLeftLeg()->getID()     << "\t"
                     << currEdge->getRightLeg()->getID()    << "\t"
                     << currEdge->getGeometry()->getID() << "\t"
                     << currEdge->getStartVertex()->getGeometry().first->getID() << "\t" 
                     << currEdge->getEndVertex()->getGeometry().first->getID() << "\t" 
                     << currEdge->getLeftFace()->getGeometry()->getID()    << "\t"
                     << currEdge->getRightFace()->getGeometry()->getID()   << "\t"
                     << currEdge->getLeftHand()->getGeometry()->getID()    << "\t"
                     << currEdge->getRightHand()->getGeometry()->getID()   << "\t"
                     << currEdge->getLeftLeg()->getGeometry()->getID()     << "\t"
                     << currEdge->getRightLeg()->getGeometry()->getID()    << endl;
            }
        }
    }
    */
    
    rg_INT inValidFace = 0;
    FaceOnSectionOfVDS* ptrCurrFace = rg_NULL;
    m_faceList.reset4Loop();
    while ( m_faceList.setNext4Loop() )  {
        ptrCurrFace = m_faceList.getpEntity();

        if ( ptrCurrFace->getFirstEdge() == rg_NULL )
            continue;

        rg_dList<VertexOnSectionOfVDS*> vertexList;

        /*
        if (    ptrCurrFace->getID() == 42 
             || ptrCurrFace->getID() == 43
             || ptrCurrFace->getID() == 44
             || ptrCurrFace->getID() == 49
             || ptrCurrFace->getID() == 67
             || ptrCurrFace->getID() == 72
             || ptrCurrFace->getID() == 80
           )  {           
            continue;
        }
        */

        ptrCurrFace->getIncidentVertexList(vertexList);

        if ( vertexList.getSize() > 0 )  {
            int aaa = 0;
        }
        else  {
            int bbb = 0;
        }

    }

    int bbb = 1;
    
}



void SectionOfVDS::adjustFacesIncidentToInfiniteVertex()
{
    /*
    VertexOnSectionOfVDS* infiniteVertex = m_vertexList.getFirstpEntity();

    EdgeOnSectionOfVDS*   firstEdge = infiniteVertex->getFirstEdge();
    EdgeOnSectionOfVDS*   currEdge  = firstEdge;    

    rg_INT i_removingFace = 0;
    FaceOnSectionOfVDS*   currFace = rg_NULL;
    do  {
        currFace = currEdge->getRightFace();

        EdgeOnSectionOfVDS* movingEdge = currEdge;    
        do {
            if ( currFace == movingEdge->getLeftFace() )
                movingEdge = movingEdge->getLeftHand();
            else
                movingEdge = movingEdge->getRightLeg();
        } while ( movingEdge->getEndVertex() != infiniteVertex);

        currFace->setFirstEdge( rg_NULL );
        currEdge = movingEdge;
        i_removingFace++;

    } while ( firstEdge != currEdge );
    */

    rg_INT numUnboundFace = 0;

    EdgeOnSectionOfVDS* currEdge = rg_NULL;
    m_edgeList.reset4Loop();
    while ( m_edgeList.setNext4Loop() )  {
        currEdge = m_edgeList.getpEntity();
        
        if (    currEdge->getStartVertex()->getID() == 0 
             || currEdge->getEndVertex()->getID()   == 0 )  {
            currEdge->getLeftFace()->setFirstEdge(rg_NULL);
            currEdge->getRightFace()->setFirstEdge(rg_NULL);
        }

    }

    rg_INT numInValidFace = 0;
    EdgeOnSectionOfVDS*   firstEdge = rg_NULL;
    FaceOnSectionOfVDS*   currFace  = rg_NULL;
    m_faceList.reset4Loop();
    while ( m_faceList.setNext4Loop() )  {
        currFace = m_faceList.getpEntity();

        if ( currFace->getFirstEdge() == rg_NULL )
            continue;

        firstEdge = currFace->getFirstEdge();
        currEdge  = firstEdge;

        rg_FLAG isValid = rg_TRUE;
        do  {
            if ( currFace == currEdge->getLeftFace() )
                currEdge = currEdge->getLeftHand();
            else if ( currFace == currEdge->getRightFace() )  
                currEdge = currEdge->getRightLeg();
            else  {
                isValid = rg_FALSE;
                break;
            }

        } while ( firstEdge != currEdge ); 
        
        if ( isValid == rg_FALSE )  {
            FaceOnSectionOfVDS* inValidFace = rg_NULL;
            if ( currFace->getGeometry() == currEdge->getLeftFace()->getGeometry() )  {
                inValidFace = currEdge->getLeftFace();
            }
            else if ( currFace->getGeometry() == currEdge->getRightFace()->getGeometry() )  { 
                inValidFace = currEdge->getRightFace();
            }
            else
                continue;
            
            EdgeOnSectionOfVDS* tempEdge = rg_NULL;
            m_edgeList.reset4Loop();
            while ( m_edgeList.setNext4Loop() )  {
                tempEdge = m_edgeList.getpEntity();
                if ( tempEdge->getLeftFace() == inValidFace )
                    tempEdge->setLeftFace( currFace );
                else if ( tempEdge->getRightFace() == inValidFace )
                    tempEdge->setRightFace( currFace );
                else;
            }

            inValidFace->setFirstEdge( rg_NULL );
            numInValidFace++;
        }
    }

    int i = 0;
}

