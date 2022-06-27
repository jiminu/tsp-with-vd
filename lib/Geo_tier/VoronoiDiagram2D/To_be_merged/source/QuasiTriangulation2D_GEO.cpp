#include "QuasiTriangulation2D_GEO.h"
using namespace BULL2D::GeometryTier;


QuasiTriangulation2D_GEO::QuasiTriangulation2D_GEO(void)
{
}



QuasiTriangulation2D_GEO::~QuasiTriangulation2D_GEO(void)
{
}

    

// void QuasiTriangulation2D_GEO::construct(CircleVoronoiDiagram& circleSetVD)
// {
//     VertexQT2D** mappingTableFromFaceInVD2VertexInQT = createAndInitializeVerticesInQT2D( circleSetVD);
//     FaceQT2D**   mappingTableFromVertexInVD2FaceInQT = createAndinitializeFacesInQT2D(    circleSetVD);
// 
// 
//     setTopologyBetweenVerticesAndTetrahedraInQT( circleSetVD, 
//                                                  mappingTableFromFaceInVD2VertexInQT, 
//                                                  mappingTableFromVertexInVD2FaceInQT);
// 
// 
//     delete [] mappingTableFromVertexInVD2FaceInQT;
//     delete [] mappingTableFromFaceInVD2VertexInQT;
// }
// 
// 
//    
// void QuasiTriangulation2D_GEO::construct(CircleSetVoronoiDiagram& circleSetVD)
// {
//     VertexQT2D** mappingTableFromFaceInVD2VertexInQT = createAndInitializeVerticesInQT2D( circleSetVD);
//     FaceQT2D**   mappingTableFromVertexInVD2FaceInQT = createAndinitializeFacesInQT2D(    circleSetVD);
// 
// 
//     setTopologyBetweenVerticesAndTetrahedraInQT( circleSetVD, 
//                                                  mappingTableFromFaceInVD2VertexInQT, 
//                                                  mappingTableFromVertexInVD2FaceInQT);
// 
// 
//     delete [] mappingTableFromVertexInVD2FaceInQT;
//     delete [] mappingTableFromFaceInVD2VertexInQT;
// }

void QuasiTriangulation2D_GEO::construct( const VoronoiDiagram2DC& circleSetVD )
{
	map<const VFace2D*, VertexQT2D*> mappingTableFromFaceInVD2VertexInQT;
	map<const VVertex2D*, FaceQT2D*>   mappingTableFromVertexInVD2FaceInQT;

	createAndInitializeVerticesInQT2D( circleSetVD, mappingTableFromFaceInVD2VertexInQT);
	createAndinitializeFacesInQT2D(    circleSetVD, mappingTableFromVertexInVD2FaceInQT);


	setTopologyBetweenVerticesAndTetrahedraInQT( circleSetVD, 
		mappingTableFromFaceInVD2VertexInQT, 
		mappingTableFromVertexInVD2FaceInQT);

}



void QuasiTriangulation2D_GEO::getFacesIncidentToVertex(VertexQT2D* vertex, rg_dList<FaceQT2D*>& faceList) const
{
    rg_dList<FaceSCDS*> incidentFaces;
    vertex->getIncidentFace( incidentFaces );

    incidentFaces.reset4Loop();
    while ( incidentFaces.setNext4Loop() ) {
        faceList.add( (FaceQT2D*) incidentFaces.getEntity() );
    }
}



void QuasiTriangulation2D_GEO::getBoundingVerticesOfFace(FaceQT2D* face, rg_dList<VertexQT2D*>&  vertexList) const
{
    VertexSCDS** vertices = face->getVertices();

    vertexList.add( (VertexQT2D*)vertices[0] );
    vertexList.add( (VertexQT2D*)vertices[1] );
    vertexList.add( (VertexQT2D*)vertices[2] );
}



void QuasiTriangulation2D_GEO::getNeighboringFacesOfFace(FaceQT2D* face, rg_dList<FaceQT2D*>&    faceList) const
{
    FaceSCDS** neighbors = face->getNeighbors();

    faceList.add( (FaceQT2D*)neighbors[0] );
    faceList.add( (FaceQT2D*)neighbors[1] );
    faceList.add( (FaceQT2D*)neighbors[2] );
}



void QuasiTriangulation2D_GEO::getBoundingVerticesOfFace(FaceQT2D* face, VertexQT2D** vertexArray) const
{
    VertexSCDS** vertices = face->getVertices();

    for ( rg_INT i=0; i<3; i++ ) {
        vertexArray[i] = (VertexQT2D*)vertices[i];
    }
}


void QuasiTriangulation2D_GEO::getNeighboringFacesOfFace(FaceQT2D* face, FaceQT2D**   faceArray) const
{
    FaceSCDS** neighbors = face->getNeighbors();

    for ( rg_INT i=0; i<3; i++ ) {
        faceArray[i] = (FaceQT2D*)neighbors[i];
    }
}



// VertexQT2D** QuasiTriangulation2D_GEO::createAndInitializeVerticesInQT2D(CircleVoronoiDiagram& circleSetVD)
// {
// 	rg_INT          numFacesInVD = circleSetVD.getNumFaces();
// 	rg_dList<Face>& facesInVD    = circleSetVD.accessFaceList();
// 
//     
// 	VertexQT2D**    mappingTableFromFaceInVD2VertexInQT = new VertexQT2D*[numFacesInVD];
// 	
//     for(rg_INT i = 0;i < numFacesInVD; i++) {
// 		mappingTableFromFaceInVD2VertexInQT[ i ] = rg_NULL;
//     }
// 
// 
// 
//     facesInVD.reset4Loop();
//     while ( facesInVD.setNext4Loop() ) {
//         Face* currFaceInVD = facesInVD.getpEntity();
// 
//         rg_INT      ID          = currFaceInVD->getID();
//         VertexQT2D* currVtxInQT = rg_NULL;
//         if ( circleSetVD.isVirtualRegion(currFaceInVD) ) {
//             currVtxInQT = m_vertices.addHead( VertexQT2D(ID, rg_NULL) );  
//         }
//         else {
//             rg_Circle2D generator = currFaceInVD->getGenerator();
//             Disc*       ptrDisc   = m_discs.add( Disc(ID, generator) );
//             
//             currVtxInQT = m_vertices.add( VertexQT2D(ID, ptrDisc) );  
//         }
//         mappingTableFromFaceInVD2VertexInQT[ID] = currVtxInQT;
//     }
// 
//     return mappingTableFromFaceInVD2VertexInQT;
// }
// 
// 
// 
// VertexQT2D** QuasiTriangulation2D_GEO::createAndInitializeVerticesInQT2D(CircleSetVoronoiDiagram& circleSetVD)
// {
// 	rg_INT             numFacesInVD = circleSetVD.getNumFaces();
// 	rg_dList<FaceVDC>& facesInVD    = circleSetVD.accessFaceList();
// 
//     
// 	VertexQT2D**    mappingTableFromFaceInVD2VertexInQT = new VertexQT2D*[numFacesInVD];
// 	
//     for(rg_INT i = 0;i < numFacesInVD; i++) {
// 		mappingTableFromFaceInVD2VertexInQT[ i ] = rg_NULL;
//     }
// 
// 
// 
//     facesInVD.reset4Loop();
//     while ( facesInVD.setNext4Loop() ) {
//         FaceVDC* currFaceInVD = facesInVD.getpEntity();
// 
//         rg_INT      ID          = currFaceInVD->getID();
//         VertexQT2D* currVtxInQT = rg_NULL;
//         if ( currFaceInVD->isOnInfinity() ) {
//             currVtxInQT = m_vertices.addHead( VertexQT2D(ID, rg_NULL) );  
//         }
//         else {
//             rg_Circle2D generator = currFaceInVD->getGeneratorGeometry();
//             Disc*       ptrDisc   = m_discs.add( Disc(ID, generator) );
//             
//             currVtxInQT = m_vertices.add( VertexQT2D(ID, ptrDisc) );  
//         }
//         mappingTableFromFaceInVD2VertexInQT[ID] = currVtxInQT;
//     }
// 
//     return mappingTableFromFaceInVD2VertexInQT;
// }


// 2014-04-15 JKKIM ADDED ///////////////////////////////////////////////////////////////////
void QuasiTriangulation2D_GEO::createAndInitializeVerticesInQT2D(const VoronoiDiagram2DC& circleSetVD, map<const VFace2D*, VertexQT2D*>& mappingTableFromFaceInVD2VertexInQT)
{
	list<const VFace2D*> VFaces;
	circleSetVD.getVoronoiFaces(VFaces);

	list<const VFace2D*>::iterator it_vFaces = VFaces.begin();
	for(; it_vFaces != VFaces.end(); it_vFaces++) {
		const VFace2D* currVFace = *it_vFaces;

		rg_INT      ID          = currVFace->getID();
		VertexQT2D* currVtxInQT = rg_NULL;
		if ( currVFace->isInfinite()) {
			currVtxInQT = m_vertices.addHead( VertexQT2D(ID, rg_NULL) );  
		}
		else {
			rg_Circle2D generator = *(currVFace->getGenerator()->getDisk());
			Disc*       ptrDisc   = m_discs.add( Disc(ID, generator) );

			currVtxInQT = m_vertices.add( VertexQT2D(ID, ptrDisc) );  
		}
		mappingTableFromFaceInVD2VertexInQT.insert(map<const VFace2D*, VertexQT2D*>::value_type(currVFace, currVtxInQT));
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////

// FaceQT2D** QuasiTriangulation2D_GEO::createAndinitializeFacesInQT2D(CircleVoronoiDiagram& circleSetVD)
// {
// 	rg_INT            numVerticesInVD = circleSetVD.getNumVertices();
// 	rg_dList<Vertex>& verticesInVD    = circleSetVD.accessVertexList();
// 
//     
// 	FaceQT2D**        mappingTableFromVertexInVD2FaceInQT = new FaceQT2D*[numVerticesInVD];
// 	
//     for(rg_INT i = 0;i < numVerticesInVD; i++) {
// 		mappingTableFromVertexInVD2FaceInQT[ i ] = rg_NULL;
//     }
// 
// 
// 
//     verticesInVD.reset4Loop();
//     while ( verticesInVD.setNext4Loop() ) {
//         Vertex*     currVtxInVD = verticesInVD.getpEntity();
// 
//         rg_INT      ID                 = currVtxInVD->getID();
//         rg_Circle2D emptyTangentCircle = currVtxInVD->getCircumcircle();
// 
//         FaceQT2D*   currFaceInQT       = m_faces.add( FaceQT2D(ID, emptyTangentCircle) );  
// 
//         mappingTableFromVertexInVD2FaceInQT[ID] = currFaceInQT;
//     }
// 
//     return mappingTableFromVertexInVD2FaceInQT;
// }
// 
// 
// 
// FaceQT2D** QuasiTriangulation2D_GEO::createAndinitializeFacesInQT2D(CircleSetVoronoiDiagram& circleSetVD)
// {
// 	rg_INT               numVerticesInVD = circleSetVD.getNumVertices();
// 	rg_dList<VertexVDC>& verticesInVD    = circleSetVD.accessVertexList();
// 
//     
// 	FaceQT2D**        mappingTableFromVertexInVD2FaceInQT = new FaceQT2D*[numVerticesInVD];
// 	
//     for(rg_INT i = 0;i < numVerticesInVD; i++) {
// 		mappingTableFromVertexInVD2FaceInQT[ i ] = rg_NULL;
//     }
// 
// 
// 
//     verticesInVD.reset4Loop();
//     while ( verticesInVD.setNext4Loop() ) {
//         VertexVDC*  currVtxInVD = verticesInVD.getpEntity();
// 
//         rg_INT      ID                 = currVtxInVD->getID();
//         rg_Circle2D emptyTangentCircle = currVtxInVD->getCircumcircle();
// 
//         FaceQT2D*   currFaceInQT       = m_faces.add( FaceQT2D(ID, emptyTangentCircle) );  
// 
//         mappingTableFromVertexInVD2FaceInQT[ID] = currFaceInQT;
//     }
// 
//     return mappingTableFromVertexInVD2FaceInQT;
// }

// 2014-04-15 JKKIM ADDED ///////////////////////////////////////////////////////////////////
void   QuasiTriangulation2D_GEO::createAndinitializeFacesInQT2D( const VoronoiDiagram2DC& circleSetVD, map<const VVertex2D*, FaceQT2D*>& mappingTableFromVertexInVD2FaceInQT)
{
	list<const VVertex2D*> VVertices;
	circleSetVD.getVoronoiVertices(VVertices);
	
	list<const VVertex2D*>::iterator it_VVertices = VVertices.begin();
	for(; it_VVertices != VVertices.end(); it_VVertices++) {
		const VVertex2D* currVVertex = *it_VVertices;		

		rg_INT      ID                 = currVVertex->getID();
		rg_Circle2D emptyTangentCircle(currVVertex->getLocation(), currVVertex->computeRadiusOfTangentCircle());

		FaceQT2D*   currFaceInQT       = m_faces.add( FaceQT2D(ID, emptyTangentCircle) );  

		mappingTableFromVertexInVD2FaceInQT.insert(map<const VVertex2D*, FaceQT2D*>::value_type(currVVertex, currFaceInQT));		
	}
	
}
/////////////////////////////////////////////////////////////////////////////////////////////



// void QuasiTriangulation2D_GEO::setTopologyBetweenVerticesAndTetrahedraInQT(CircleVoronoiDiagram& circleSetVD,
//                                                                        VertexQT2D** mappingTableFromFaceInVD2VertexInQT,
//                                                                        FaceQT2D**   mappingTableFromVertexInVD2FaceInQT)
// {
//     rg_INT i=0;
// 
// 	rg_dList<Vertex>& verticesInVD = circleSetVD.accessVertexList();
// 
//     verticesInVD.reset4Loop();
//     while ( verticesInVD.setNext4Loop() ) {
//         Vertex* currVtxInVD = verticesInVD.getpEntity();
//         rg_INT  vtxID       = currVtxInVD->getID();
// 
//         //if ( currVtxInVD->isExist() == rg_FALSE ) {
//         //    continue;
//         //}
// 
//         WingedEdge* incidentEdgeOfCurrVtx[3] = { currVtxInVD->getFirstEdge(),
//                                                  incidentEdgeOfCurrVtx[0]->getCCWNextEdgeOnVertex( currVtxInVD ), 
//                                                  incidentEdgeOfCurrVtx[1]->getCCWNextEdgeOnVertex( currVtxInVD )  }; 
//         Vertex*     adjacentVtx[3]           = {rg_NULL, rg_NULL, rg_NULL};
//         Face*       mateFaceOfAdjacentVtx[3] = {rg_NULL, rg_NULL, rg_NULL};
//         for ( i=0; i<3; i++ ) {
//             if ( incidentEdgeOfCurrVtx[i]->isStartVertex( currVtxInVD ) ) {
//                 adjacentVtx[i] = incidentEdgeOfCurrVtx[i]->getEndVertex();
//             }
//             else {
//                 adjacentVtx[i] = incidentEdgeOfCurrVtx[i]->getStartVertex();
//             }
// 
//             mateFaceOfAdjacentVtx[i] = incidentEdgeOfCurrVtx[i]->getMateFace(adjacentVtx[i]);
//         }
// 
// 
//         //  set topology of face in QuasiTriangulation2D.
//         FaceQT2D*   currFaceInQT   = mappingTableFromVertexInVD2FaceInQT[vtxID];
//         VertexQT2D* boundingVtx[3] = { mappingTableFromFaceInVD2VertexInQT[ mateFaceOfAdjacentVtx[0]->getID() ], 
//                                        mappingTableFromFaceInVD2VertexInQT[ mateFaceOfAdjacentVtx[1]->getID() ], 
//                                        mappingTableFromFaceInVD2VertexInQT[ mateFaceOfAdjacentVtx[2]->getID() ] };
//         currFaceInQT->setVertices(  boundingVtx[0], boundingVtx[1], boundingVtx[2] ); 
//         currFaceInQT->setNeighbors( mappingTableFromVertexInVD2FaceInQT[ adjacentVtx[0]->getID() ],
//                                     mappingTableFromVertexInVD2FaceInQT[ adjacentVtx[1]->getID() ],
//                                     mappingTableFromVertexInVD2FaceInQT[ adjacentVtx[2]->getID() ] );
// 
// 
//         //  set topology of face in QuasiTriangulation2D.
//         for ( i=0; i<3; i++ ) {
//             boundingVtx[i]->setFirstFace( currFaceInQT );
//         }
//     }
// }
// 
// 
// 
// void QuasiTriangulation2D_GEO::setTopologyBetweenVerticesAndTetrahedraInQT(CircleSetVoronoiDiagram& circleSetVD,
//                                                                        VertexQT2D** mappingTableFromFaceInVD2VertexInQT,
//                                                                        FaceQT2D**   mappingTableFromVertexInVD2FaceInQT)
// {
//     rg_INT i=0;
// 
// 	rg_dList<VertexVDC>& verticesInVD = circleSetVD.accessVertexList();
// 
//     verticesInVD.reset4Loop();
//     while ( verticesInVD.setNext4Loop() ) {
//         VertexVDC* currVtxInVD = verticesInVD.getpEntity();
//         rg_INT     vtxID       = currVtxInVD->getID();
// 
//         EdgeVDC*   incidentEdgeOfCurrVtx[3] = { currVtxInVD->getFirstEdge(),
//                                                 currVtxInVD->getCCWNextEdge(incidentEdgeOfCurrVtx[0]),
//                                                 currVtxInVD->getCCWNextEdge(incidentEdgeOfCurrVtx[1]) };
//         VertexVDC* adjacentVtx[3]           = {rg_NULL, rg_NULL, rg_NULL};
//         FaceVDC*   mateFaceOfAdjacentVtx[3] = {rg_NULL, rg_NULL, rg_NULL};
//         for ( i=0; i<3; i++ ) {
//             if ( incidentEdgeOfCurrVtx[i]->isStartVertex( currVtxInVD ) ) {
//                 adjacentVtx[i] = incidentEdgeOfCurrVtx[i]->getEndVertex();
//             }
//             else {
//                 adjacentVtx[i] = incidentEdgeOfCurrVtx[i]->getStartVertex();
//             }
// 
//             mateFaceOfAdjacentVtx[i] = incidentEdgeOfCurrVtx[i]->getMateFace(adjacentVtx[i]);
//         }
// 
// 
//         //  set topology of face in QuasiTriangulation2D.
//         FaceQT2D*   currFaceInQT   = mappingTableFromVertexInVD2FaceInQT[vtxID];
//         VertexQT2D* boundingVtx[3] = { mappingTableFromFaceInVD2VertexInQT[ mateFaceOfAdjacentVtx[0]->getID() ], 
//                                        mappingTableFromFaceInVD2VertexInQT[ mateFaceOfAdjacentVtx[1]->getID() ], 
//                                        mappingTableFromFaceInVD2VertexInQT[ mateFaceOfAdjacentVtx[2]->getID() ] };
//         currFaceInQT->setVertices(  boundingVtx[0], boundingVtx[1], boundingVtx[2] ); 
//         currFaceInQT->setNeighbors( mappingTableFromVertexInVD2FaceInQT[ adjacentVtx[0]->getID() ],
//                                     mappingTableFromVertexInVD2FaceInQT[ adjacentVtx[1]->getID() ],
//                                     mappingTableFromVertexInVD2FaceInQT[ adjacentVtx[2]->getID() ] );
// 
// 
//         //  set topology of face in QuasiTriangulation2D.
//         for ( i=0; i<3; i++ ) {
//             boundingVtx[i]->setFirstFace( currFaceInQT );
//         }
//     }
// }


// 2014-04-15 JKKIM ADDED ///////////////////////////////////////////////////////////////////
void QuasiTriangulation2D_GEO::setTopologyBetweenVerticesAndTetrahedraInQT( const VoronoiDiagram2DC& circleSetVD, 
																		    const map<const VFace2D*, VertexQT2D*>& mappingTableFromFaceInVD2VertexInQT,
																		    const map<const VVertex2D*, FaceQT2D*>&   mappingTableFromVertexInVD2FaceInQT )
{
	rg_INT i=0;

	list<const VVertex2D*> VVertices;
	circleSetVD.getVoronoiVertices(VVertices);

	list<const VVertex2D*>::iterator it_VVertices = VVertices.begin();
	for(; it_VVertices != VVertices.end(); it_VVertices++) {
		const VVertex2D* currVVertex = *it_VVertices;		
		rg_INT     vtxID       = currVVertex->getID();

		list<VEdge2D*> VEdges;
		currVVertex->getIncident3VEdges(VEdges);
				
		VVertex2D* adjacentVtx[3]           = {rg_NULL, rg_NULL, rg_NULL};
		VFace2D*   mateFaceOfAdjacentVtx[3] = {rg_NULL, rg_NULL, rg_NULL};

		list<VEdge2D*>::iterator it_VEdge = VEdges.begin();
		for(i = 0;it_VEdge != VEdges.end(); it_VEdge++, i++	) {		
			VEdge2D* currVEdge = *it_VEdge;

			if ( currVEdge->getStartVertex() == currVVertex ) {
				adjacentVtx[i] = currVEdge->getEndVertex();
			}
			else {
				adjacentVtx[i] = currVEdge->getStartVertex();
			}

			mateFaceOfAdjacentVtx[i] = currVEdge->getMateFace(adjacentVtx[i]);
		}


		//  set topology of face in QuasiTriangulation2D.		
		FaceQT2D*   currFaceInQT   = (*mappingTableFromVertexInVD2FaceInQT.find( currVVertex )).second;
		VertexQT2D* boundingVtx[3] = { (*mappingTableFromFaceInVD2VertexInQT.find( mateFaceOfAdjacentVtx[0] )).second, 
			(*mappingTableFromFaceInVD2VertexInQT.find( mateFaceOfAdjacentVtx[1] )).second, 
			(*mappingTableFromFaceInVD2VertexInQT.find( mateFaceOfAdjacentVtx[2] )).second };
		currFaceInQT->setVertices(  boundingVtx[0], boundingVtx[1], boundingVtx[2] ); 
		currFaceInQT->setNeighbors( (*mappingTableFromVertexInVD2FaceInQT.find( adjacentVtx[0] )).second,
			(*mappingTableFromVertexInVD2FaceInQT.find( adjacentVtx[1] )).second,
			(*mappingTableFromVertexInVD2FaceInQT.find( adjacentVtx[2] )).second );


		//  set topology of face in QuasiTriangulation2D.
		for ( i=0; i<3; i++ ) {
			boundingVtx[i]->setFirstFace( currFaceInQT );
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////