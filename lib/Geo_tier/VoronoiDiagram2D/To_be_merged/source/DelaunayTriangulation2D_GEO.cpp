#include "DelaunayTriangulation2D_GEO.h"
using namespace BULL2D::GeometryTier;


DelaunayTriangulation2D_GEO::DelaunayTriangulation2D_GEO(void)
{
}



DelaunayTriangulation2D_GEO::~DelaunayTriangulation2D_GEO(void)
{
}

    



void DelaunayTriangulation2D_GEO::construct( VoronoiDiagram2DP& pointSetVD )
{
	map<const VFace2DP*, VertexQT2D*> mappingTableFromFaceInVD2VertexInQT;
	map<const VVertex2DP*, FaceQT2D*>   mappingTableFromVertexInVD2FaceInQT;

	createAndInitializeVerticesInQT2D( pointSetVD, mappingTableFromFaceInVD2VertexInQT);
	createAndinitializeFacesInQT2D(    pointSetVD, mappingTableFromVertexInVD2FaceInQT);


	setTopologyBetweenVerticesAndTetrahedraInQT( pointSetVD, 
		mappingTableFromFaceInVD2VertexInQT, 
		mappingTableFromVertexInVD2FaceInQT);

}



void DelaunayTriangulation2D_GEO::getFacesIncidentToVertex(VertexQT2D* vertex, rg_dList<FaceQT2D*>& faceList) const
{
    rg_dList<FaceSCDS*> incidentFaces;
    vertex->getIncidentFace( incidentFaces );

    incidentFaces.reset4Loop();
    while ( incidentFaces.setNext4Loop() ) {
        faceList.add( (FaceQT2D*) incidentFaces.getEntity() );
    }
}



void DelaunayTriangulation2D_GEO::getBoundingVerticesOfFace(FaceQT2D* face, rg_dList<VertexQT2D*>&  vertexList) const
{
    VertexSCDS** vertices = face->getVertices();

    vertexList.add( (VertexQT2D*)vertices[0] );
    vertexList.add( (VertexQT2D*)vertices[1] );
    vertexList.add( (VertexQT2D*)vertices[2] );
}



void DelaunayTriangulation2D_GEO::getNeighboringFacesOfFace(FaceQT2D* face, rg_dList<FaceQT2D*>&    faceList) const
{
    FaceSCDS** neighbors = face->getNeighbors();

    faceList.add( (FaceQT2D*)neighbors[0] );
    faceList.add( (FaceQT2D*)neighbors[1] );
    faceList.add( (FaceQT2D*)neighbors[2] );
}



void DelaunayTriangulation2D_GEO::getBoundingVerticesOfFace(FaceQT2D* face, VertexQT2D** vertexArray) const
{
    VertexSCDS** vertices = face->getVertices();

    for ( rg_INT i=0; i<3; i++ ) {
        vertexArray[i] = (VertexQT2D*)vertices[i];
    }
}


void DelaunayTriangulation2D_GEO::getNeighboringFacesOfFace(FaceQT2D* face, FaceQT2D**   faceArray) const
{
    FaceSCDS** neighbors = face->getNeighbors();

    for ( rg_INT i=0; i<3; i++ ) {
        faceArray[i] = (FaceQT2D*)neighbors[i];
    }
}



// 2014-04-15 JKKIM ADDED ///////////////////////////////////////////////////////////////////
void DelaunayTriangulation2D_GEO::createAndInitializeVerticesInQT2D(VoronoiDiagram2DP& pointSetVD, map<const VFace2DP*, VertexQT2D*>& mappingTableFromFaceInVD2VertexInQT)
{
	list<const VFace2DP*> VFaces;
	pointSetVD.getVoronoiFaces(VFaces);

	list<const VFace2DP*>::iterator it_vFaces = VFaces.begin();
	for(; it_vFaces != VFaces.end(); it_vFaces++) {
		const VFace2DP* currVFace = *it_vFaces;

		rg_INT      ID          = currVFace->getID();
		VertexQT2D* currVtxInQT = rg_NULL;
		if ( currVFace->isInfinite()) {
			currVtxInQT = m_vertices.addHead( VertexQT2D(ID, rg_NULL) );  
		}
		else {
			rg_Circle2D generator(*(currVFace->getGenerator()->getLocation()), 0.0);
			Disc*       ptrDisc   = m_discs.add( Disc(ID, generator) );

			currVtxInQT = m_vertices.add( VertexQT2D(ID, ptrDisc) );  
		}
		mappingTableFromFaceInVD2VertexInQT.insert(map<const VFace2DP*, VertexQT2D*>::value_type(currVFace, currVtxInQT));
	}
}


// 2014-04-15 JKKIM ADDED ///////////////////////////////////////////////////////////////////
void   DelaunayTriangulation2D_GEO::createAndinitializeFacesInQT2D(   VoronoiDiagram2DP& pointSetVD, map<const VVertex2DP*, FaceQT2D*>& mappingTableFromVertexInVD2FaceInQT)
{
	list<const VVertex2DP*> VVertices;
	pointSetVD.getVoronoiVertices(VVertices);
	
	list<const VVertex2DP*>::iterator it_VVertices = VVertices.begin();
	for(; it_VVertices != VVertices.end(); it_VVertices++) {
		const VVertex2DP* currVVertex = *it_VVertices;		

		rg_INT      ID                 = currVVertex->getID();
		rg_Circle2D emptyTangentCircle(currVVertex->getLocation(), currVVertex->computeRadiusOfTangentCircle());

		FaceQT2D*   currFaceInQT       = m_faces.add( FaceQT2D(ID, emptyTangentCircle) );  

		mappingTableFromVertexInVD2FaceInQT.insert(map<const VVertex2DP*, FaceQT2D*>::value_type(currVVertex, currFaceInQT));		
	}
	
}


// 2014-04-15 JKKIM ADDED ///////////////////////////////////////////////////////////////////
void DelaunayTriangulation2D_GEO::setTopologyBetweenVerticesAndTetrahedraInQT( VoronoiDiagram2DP& pointSetVD, 
																		const map<const VFace2DP*, VertexQT2D*>& mappingTableFromFaceInVD2VertexInQT,
																		const map<const VVertex2DP*, FaceQT2D*>&   mappingTableFromVertexInVD2FaceInQT )
{
	rg_INT i=0;

	list<const VVertex2DP*> VVertices;
	pointSetVD.getVoronoiVertices(VVertices);

	list<const VVertex2DP*>::iterator it_VVertices = VVertices.begin();
	for(; it_VVertices != VVertices.end(); it_VVertices++) {
		const VVertex2DP* currVVertex = *it_VVertices;		
		rg_INT     vtxID       = currVVertex->getID();

		list<VEdge2DP*> VEdges;
		currVVertex->getIncident3VEdges(VEdges);
				
		VVertex2DP* adjacentVtx[3]           = {rg_NULL, rg_NULL, rg_NULL};
		VFace2DP*   mateFaceOfAdjacentVtx[3] = {rg_NULL, rg_NULL, rg_NULL};

		list<VEdge2DP*>::iterator it_VEdge = VEdges.begin();
		for(i = 0;it_VEdge != VEdges.end(); it_VEdge++, i++	) {		
			VEdge2DP* currVEdge = *it_VEdge;

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