#include "BetaUniverse2D.h"
using namespace BULL2D::GeometryTier;




BetaUniverse2D::BetaUniverse2D()
{
	m_sortedFiniteVertices = rg_NULL;
	m_sortedFiniteEdges    = rg_NULL;
	m_sortedFiniteFaces    = rg_NULL;
}




BetaUniverse2D::~BetaUniverse2D()
{
    if ( m_sortedFiniteFaces != rg_NULL ) {
        delete [] m_sortedFiniteFaces;
    }
    if ( m_sortedFiniteEdges != rg_NULL ) {
        delete [] m_sortedFiniteEdges;
    }
    if ( m_sortedFiniteVertices != rg_NULL ) {
        delete [] m_sortedFiniteVertices;
    }
}



void BetaUniverse2D::construct(QuasiTriangulation2D_GEO& QT2D)
{
    convert( QT2D );

    computeBetaSpan();

    computeBetaInterval();

    generateSortedFiniteSimplexes();

    rearrangeSimplexIDs();
}



void BetaUniverse2D::clean()
{
    if ( m_sortedFiniteFaces != rg_NULL ) {
        delete [] m_sortedFiniteFaces;
    }
    if ( m_sortedFiniteEdges != rg_NULL ) {
        delete [] m_sortedFiniteEdges;
    }
    if ( m_sortedFiniteVertices != rg_NULL ) {
        delete [] m_sortedFiniteVertices;
    }

    m_numFiniteFaces    = 0;
    m_numFiniteEdges    = 0;
    m_numFiniteVertices = 0;

    m_faces.removeAll();
    m_edges.removeAll();
    m_vertices.removeAll();
    m_discs.removeAll();

    m_betaInterval = Interval_REAL();
}



// Query functions
void  BetaUniverse2D::huntFacesInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const
{
    startIndex = 0;
    endIndex   = 0;

    if ( beta < m_sortedFiniteFaces[0]->getLowerBoundForValidSimplex() ) {
        startIndex = 0;
        endIndex   = -1;
    }
    else {
        for ( rg_INT i = 1; i < m_numFiniteFaces; i++) {
            if ( m_sortedFiniteFaces[i]->getLowerBoundForValidSimplex() > beta ) {
                break;
            }
            else {
                endIndex = i;
            }
        }
    }
}



void  BetaUniverse2D::huntEdgesInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const
{
    startIndex = 0;
    endIndex   = 0;

    if ( beta < m_sortedFiniteEdges[0]->getLowerBoundForValidSimplex() ) {
        startIndex = 0;
        endIndex   = -1;
    }
    else {
        for ( rg_INT i = 1; i < m_numFiniteFaces; i++) {
            if ( m_sortedFiniteEdges[i]->getLowerBoundForValidSimplex() > beta ) {
                break;
            }
            else {
                endIndex = i;
            }
        }
    }
}



void  BetaUniverse2D::huntVerticesInBetaComplex( const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const
{
    startIndex = 0;
    endIndex   = m_numFiniteVertices-1;
}




void  BetaUniverse2D::huntFacesInBetaComplex(    const rg_REAL& beta, rg_dList<FaceBU2D*>&   betaFaces ) const
{
    for ( rg_INT i=0; i < m_numFiniteFaces; i++ ) {
        if ( beta >= m_sortedFiniteFaces[i]->getLowerBoundForValidSimplex() ) {
            betaFaces.add( m_sortedFiniteFaces[i] );
		}
        else {
            break;
        }
    }
}



void  BetaUniverse2D::huntEdgesInBetaComplex(    const rg_REAL& beta, rg_dList<EdgeBU2D*>&   betaEdges ) const
{
    for ( rg_INT i=0; i < m_numFiniteEdges; i++ ) {
        if ( beta >= m_sortedFiniteEdges[i]->getLowerBoundForValidSimplex() ) {
            betaEdges.add( m_sortedFiniteEdges[i] );
        }
        else {
            break;
        }
    }
}



void  BetaUniverse2D::huntVerticesInBetaComplex( const rg_REAL& beta, rg_dList<VertexBU2D*>& betaVertices ) const
{
    for ( rg_INT i=0; i < m_numFiniteVertices; i++ ) { 
		betaVertices.add( m_sortedFiniteVertices[i] );
    }
}




void  BetaUniverse2D::huntEdgesInBetaShape(    const rg_REAL& beta, rg_dList<EdgeBU2D*>&   betaEdges ) const
{
    for ( rg_INT i=0; i < m_numFiniteEdges; i++ ) {
		rg_INT state = m_sortedFiniteEdges[i]->getBoundingState(beta);

        if( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
			betaEdges.add(m_sortedFiniteEdges[i]);
        }
    }
}



void  BetaUniverse2D::huntVerticesInBetaShape( const rg_REAL& beta, rg_dList<VertexBU2D*>& betaVertices ) const
{
    for ( rg_INT i=0; i < m_numFiniteVertices; i++ ) {
		rg_INT state = m_sortedFiniteVertices[i]->getBoundingState(beta);

        if( state == SINGULAR_SIMPLEX || state == REGULAR_SIMPLEX ) {
			betaVertices.add(m_sortedFiniteVertices[i]);
        }
    }
}



void    BetaUniverse2D::reportBetaUniverse2D(const string& reportFilename)
{
    ofstream fout;
    fout.open( reportFilename.c_str() );

    // disc
    fout << "### Discs : " << m_discs.getSize() << endl;
    m_discs.reset4Loop();
    while ( m_discs.setNext4Loop() ) {
        Disc* currDisc = m_discs.getpEntity();
        rg_Circle2D circle = currDisc->getGeometry();

        fout << "d" << currDisc->getID() << "\t";
        fout << circle.getCenterPt().getX() << "\t";
        fout << circle.getCenterPt().getY() << "\t";
        fout << circle.getRadius() << endl;

    }
    fout << endl << endl;

    // vertex
    fout << "### Vertex : " << m_vertices.getSize() << endl;
    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() ) {
        VertexBU2D* currVtx = m_vertices.getpEntity();

        if ( currVtx->isVirtual() ) {
            continue;
        }

    	BetaSpan betaspan = currVtx->getBetaSpan();


        fout << "v" << currVtx->getID() << "\t";
        fout << "d" << currVtx->getDisc()->getID() << "\t";

        rg_INT   numBI = betaspan.getNumBetaInterval();
        rg_INT*  boundingState = betaspan.getBoundingStates();
        rg_REAL* betaInterval  = betaspan.getBetaIntervals();
        rg_INT   i = 0;
        for ( i=0; i<numBI; i++ ) {
            fout << betaInterval[i] << "\t";
            switch ( boundingState[i] ) {
                case EXTRANEOUS_SIMPLEX:
                    fout << "E" << "\t";
                    break;
                case SINGULAR_SIMPLEX:
                    fout << "S" << "\t";
                    break;
                case REGULAR_SIMPLEX:
                    fout << "R" << "\t";
                    break;
                case INTERIOR_SIMPLEX:
                    fout << "I" << "\t";
                    break;
            }
        }
        fout << betaInterval[numBI] << endl;
    }
    fout << endl << endl;

    // edge
    fout << "### Edge : " << m_edges.getSize() << endl;
    m_edges.reset4Loop();
    while ( m_edges.setNext4Loop() ) {
        EdgeBU2D* currEdge = m_edges.getpEntity();

        if ( currEdge->isVirtual() ) {
            continue;
        }

    	BetaSpan betaspan = currEdge->getBetaSpan();

        fout << "e" << currEdge->getID() << "\t";
        fout << "v" << currEdge->getStartVertex()->getID() << "->";
        fout << "v" << currEdge->getEndVertex()->getID() << "\t";
        fout << "lf" << currEdge->getLeftFace()->getID() << "->";
        fout << "rf" << currEdge->getRightFace()->getID() << "\t";

        rg_INT   numBI = betaspan.getNumBetaInterval();
        rg_INT*  boundingState = betaspan.getBoundingStates();
        rg_REAL* betaInterval  = betaspan.getBetaIntervals();
        rg_INT   i = 0;
        for ( i=0; i<numBI; i++ ) {
            fout << betaInterval[i] << "\t";
            switch ( boundingState[i] ) {
                case EXTRANEOUS_SIMPLEX:
                    fout << "E" << "\t";
                    break;
                case SINGULAR_SIMPLEX:
                    fout << "S" << "\t";
                    break;
                case REGULAR_SIMPLEX:
                    fout << "R" << "\t";
                    break;
                case INTERIOR_SIMPLEX:
                    fout << "I" << "\t";
                    break;
            }
        }
        fout << betaInterval[numBI] << endl;
    }
    fout << endl << endl;

    // face
    fout << "### Face : " << m_faces.getSize() << endl;
    m_faces.reset4Loop();
    while ( m_faces.setNext4Loop() ) {
        FaceBU2D* currFace = m_faces.getpEntity();

        if ( currFace->isVirtual() ) {
            continue;
        }

    	BetaSpan betaspan = currFace->getBetaSpan();

        fout << "f" << currFace->getID() << "\t";

        rg_dList<VertexBU2D*> boundingVertices;
        currFace->getBoundingVertices(boundingVertices);
        fout << "v" << boundingVertices.getFirstEntity()->getID() << "->";
        fout << "v" << boundingVertices.getSecondEntity()->getID() << "->";
        fout << "v" << boundingVertices.getLastEntity()->getID() << "\t";


        rg_INT   numBI = betaspan.getNumBetaInterval();
        rg_INT*  boundingState = betaspan.getBoundingStates();
        rg_REAL* betaInterval  = betaspan.getBetaIntervals();
        rg_INT   i = 0;
        for ( i=0; i<numBI; i++ ) {
            fout << betaInterval[i] << "\t";
            switch ( boundingState[i] ) {
                case EXTRANEOUS_SIMPLEX:
                    fout << "E" << "\t";
                    break;
                case SINGULAR_SIMPLEX:
                    fout << "S" << "\t";
                    break;
                case REGULAR_SIMPLEX:
                    fout << "R" << "\t";
                    break;
                case INTERIOR_SIMPLEX:
                    fout << "I" << "\t";
                    break;
            }
        }
        fout << betaInterval[numBI] << endl;
    }
    fout << endl << endl;

    fout.close();
}



void BetaUniverse2D::convert(QuasiTriangulation2D_GEO& QT2D)
{
	map<FaceSCDS*, FaceBU2D*>		  mappingTableFromFaceInQTToFaceInBU;
	map<VertexSCDS*, VertexBU2D*>   mappingTableFromVtxInQTToVtxInBU;
	


    createAndInitializeFacesInBU2D( QT2D, mappingTableFromFaceInQTToFaceInBU );
    createAndInitializeVerticesInBU2D( QT2D, mappingTableFromVtxInQTToVtxInBU );

    rg_INT i=0;
	rg_INT              numFacesInQT2D = QT2D.getNumFaces();
	rg_dList<FaceQT2D>& facesInQT2D    = QT2D.getFaces();
    
	map<FaceSCDS*, EdgeBU2D**>   edgesInBU;
		
	facesInQT2D.reset4Loop();
	while ( facesInQT2D.setNext4Loop() ) {
		FaceQT2D* currFaceInQT2D = facesInQT2D.getpEntity();

		EdgeBU2D** currEdgeInBU = new EdgeBU2D*[3];
		for(int i = 0; i <3; i++) {
			currEdgeInBU[i] = rg_NULL;
		}

		edgesInBU.insert(map<FaceSCDS*, EdgeBU2D**>::value_type(currFaceInQT2D, currEdgeInBU));		
	}


    facesInQT2D.reset4Loop();
    while ( facesInQT2D.setNext4Loop() ) {
        FaceQT2D* currFaceInQT2D = facesInQT2D.getpEntity();
		
        rg_INT    faceID         = currFaceInQT2D->getID();
        FaceBU2D* currFaceInBU   = (*mappingTableFromFaceInQTToFaceInBU.find(currFaceInQT2D)).second;

        VertexSCDS** vertexInQT    = currFaceInQT2D->getVertices();
        VertexBU2D*  vertexInBU[4] = { (*mappingTableFromVtxInQTToVtxInBU.find( vertexInQT[1] )).second,
                                       (*mappingTableFromVtxInQTToVtxInBU.find( vertexInQT[2] )).second,
                                       (*mappingTableFromVtxInQTToVtxInBU.find( vertexInQT[0] )).second,
                                       (*mappingTableFromVtxInQTToVtxInBU.find( vertexInQT[1] )).second };

        FaceSCDS**   neighborInQT    = currFaceInQT2D->getNeighbors();
        FaceBU2D*    neighborInBU[3] = { (*mappingTableFromFaceInQTToFaceInBU.find( neighborInQT[0])).second,
                                         (*mappingTableFromFaceInQTToFaceInBU.find( neighborInQT[1])).second,
                                         (*mappingTableFromFaceInQTToFaceInBU.find( neighborInQT[2])).second };

		
		EdgeBU2D** currEdgeInBU = (*edgesInBU.find(currFaceInQT2D)).second;
        for ( i=0; i<3; i++ ) {
			if(currEdgeInBU[i] == rg_NULL) {
				currEdgeInBU[i] = m_edges.add( EdgeBU2D( m_edges.getSize() ) );

				currEdgeInBU[i]->setStartVertex(   vertexInBU[i] );
				currEdgeInBU[i]->setEndVertex(     vertexInBU[i+1] );
				vertexInBU[i]->setFirstEdge(            currEdgeInBU[i] );
				vertexInBU[i+1]->setFirstEdge(          currEdgeInBU[i] );

				currEdgeInBU[i]->connectLeftFace(  currFaceInBU    );
				currEdgeInBU[i]->connectRightFace( neighborInBU[i] );


				rg_INT posInNeighbor = currFaceInQT2D->getThisPosInNeighbor(i, neighborInQT[i]);

				EdgeBU2D** NeighborEdgeInBU = (*edgesInBU.find(neighborInQT[i])).second;
				NeighborEdgeInBU[posInNeighbor] = currEdgeInBU[i];            
			}
        }

		

		


        EdgeBU2D*    edgeInBU[5] = { currEdgeInBU[0], currEdgeInBU[1], currEdgeInBU[2], 
                                     currEdgeInBU[0], currEdgeInBU[1] };
        for ( i=0; i<3; i++ ) {
            if ( currFaceInBU == edgeInBU[i]->getLeftFace() ) {
                edgeInBU[i]->connectLeftHand( edgeInBU[i+1] );
                edgeInBU[i]->connectLeftLeg(  edgeInBU[i+2] );
            }
            else {
                edgeInBU[i]->connectRightHand( edgeInBU[i+2] );
                edgeInBU[i]->connectRightLeg(  edgeInBU[i+1] );
            }
        }
    }


	facesInQT2D.reset4Loop();
	while ( facesInQT2D.setNext4Loop() ) {
		FaceQT2D* currFaceInQT2D = facesInQT2D.getpEntity();

		EdgeBU2D** currEdgeInBU = (*edgesInBU.find(currFaceInQT2D)).second;

		delete[] currEdgeInBU;
	}

	edgesInBU.clear();
}



void BetaUniverse2D::computeBetaSpan()
{
    //  compute beta-spans of faces
	m_faces.reset4Loop();
	while ( m_faces.setNext4Loop() ) {
		FaceBU2D* currFace = m_faces.getpEntity();
		currFace->computeBetaSpan();
	}

    //  compute beta-spans of edges
	m_edges.reset4Loop();
	while ( m_edges.setNext4Loop() ) {
		EdgeBU2D* currEdge = m_edges.getpEntity();

#if RG_DEBUG
        if ( currEdge->getID() == 8 || currEdge->getID() == 8 ) {
            int stop = 1;
        }
#endif

		currEdge->computeBetaSpan();
	}

    //  compute beta-spans of vertices
	m_vertices.reset4Loop();
	while ( m_vertices.setNext4Loop() ) {
		VertexBU2D* currVtx = m_vertices.getpEntity();
		currVtx->computeBetaSpan();
	}
}



void BetaUniverse2D::generateSortedFiniteSimplexes()
{
    generateSortedFiniteFaces();
    generateSortedFiniteEdges();
    generateSortedFiniteVertices();
}



void BetaUniverse2D::generateSortedFiniteFaces()
{
    m_numFiniteFaces = 0;

    FaceBU2D* currFace = rg_NULL;
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		currFace = m_faces.getpEntity();

        if ( !currFace->isVirtual() ) {
            m_numFiniteFaces++;
        }
	}

    m_sortedFiniteFaces = new FaceBU2D*[m_numFiniteFaces];
    rg_INT i = 0;
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		currFace = m_faces.getpEntity();

        if ( !currFace->isVirtual() ) {
            m_sortedFiniteFaces[i] = currFace;
            i++;
        }
	}

    
    qsort( (void *) m_sortedFiniteFaces, m_numFiniteFaces, sizeof(FaceBU2D*), compareBetaSpanOfFaces);
}



void BetaUniverse2D::generateSortedFiniteEdges()
{
    m_numFiniteEdges = 0;

    EdgeBU2D* currEdge = rg_NULL;
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		currEdge = m_edges.getpEntity();

        if ( !currEdge->isVirtual() ) {
            m_numFiniteEdges++;
        }
	}

    m_sortedFiniteEdges = new EdgeBU2D*[m_numFiniteEdges];
    rg_INT i = 0;
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		currEdge = m_edges.getpEntity();

        if ( !currEdge->isVirtual() ) {
            m_sortedFiniteEdges[i] = currEdge;
            i++;
        }
	}


    qsort( (void *) m_sortedFiniteEdges, m_numFiniteEdges, sizeof(EdgeBU2D*), compareBetaSpanOfEdges);
}



void BetaUniverse2D::generateSortedFiniteVertices()
{
    m_numFiniteVertices = 0;

    VertexBU2D* currVertex = rg_NULL;
	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		currVertex = m_vertices.getpEntity();

        if ( !currVertex->isVirtual() ) {
            m_numFiniteVertices++;
        }
	}

    m_sortedFiniteVertices = new VertexBU2D*[m_numFiniteVertices];
    rg_INT i = 0;
	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		currVertex = m_vertices.getpEntity();

        if ( !currVertex->isVirtual() ) {
            m_sortedFiniteVertices[i] = currVertex;
            i++;
        }
	}


    qsort( (void *) m_sortedFiniteVertices, m_numFiniteVertices, sizeof(VertexBU2D*), compareBetaSpanOfVertices);
}





void BetaUniverse2D::createAndInitializeFacesInBU2D( QuasiTriangulation2D_GEO& QT2D,    map<FaceSCDS*, FaceBU2D*>&     mappingTableFromFaceInQTToFaceInBU )
{	
	rg_dList<FaceQT2D>& facesInQT2D    = QT2D.getFaces();
 	
    facesInQT2D.reset4Loop();
    while ( facesInQT2D.setNext4Loop() ) {
        FaceQT2D*   currFaceInQT2D = facesInQT2D.getpEntity();

        rg_INT      ID                 = currFaceInQT2D->getID();
        rg_Circle2D emptyTangentCircle = currFaceInQT2D->getEmptyTangentCircle();

        FaceBU2D*   currFaceInBU       = m_faces.add( FaceBU2D(ID, emptyTangentCircle) );  

        mappingTableFromFaceInQTToFaceInBU.insert(map<FaceSCDS*, FaceBU2D*>::value_type(currFaceInQT2D, currFaceInBU));
    }


}


void BetaUniverse2D::createAndInitializeVerticesInBU2D( QuasiTriangulation2D_GEO& QT2D, map<VertexSCDS*, VertexBU2D*>& mappingTableFromVtxInQTToVtxInBU )
{	
	rg_dList<VertexQT2D>& verticesInQT2D = QT2D.getVertices();
    
	
    verticesInQT2D.reset4Loop();
    while ( verticesInQT2D.setNext4Loop() ) {
        VertexQT2D*   currVtxInQT2D = verticesInQT2D.getpEntity();

        rg_INT      ID          = currVtxInQT2D->getID();
        VertexBU2D* currVtxInBU = rg_NULL;
        if ( currVtxInQT2D->isVirtual() ) {
            currVtxInBU = m_vertices.addHead( VertexBU2D(ID, rg_NULL) );  
        }
        else {
            rg_Circle2D disc    = currVtxInQT2D->getDiscGeometry();
            Disc*       ptrDisc = m_discs.add( Disc(ID, disc) );

            currVtxInBU = m_vertices.add( VertexBU2D(ID, ptrDisc) );  
        }

        mappingTableFromVtxInQTToVtxInBU.insert(map<VertexSCDS*, VertexBU2D*>::value_type(currVtxInQT2D, currVtxInBU));
    }	   
}



int BULL2D::GeometryTier::compareBetaSpanOfFaces(const void* ts1, const void* ts2)
{
    FaceBU2D* face1 = *(FaceBU2D**) ts1;
    FaceBU2D* face2 = *(FaceBU2D**) ts2;


    rg_REAL lowerBoundOfFace1 = face1->getLowerBoundForValidSimplex();
    rg_REAL lowerBoundOfFace2 = face2->getLowerBoundForValidSimplex();
    if ( lowerBoundOfFace1 > lowerBoundOfFace2 ) {
        return 1;
    }
    else if ( lowerBoundOfFace1 < lowerBoundOfFace2 ) {
        return -1;
    }
    else {
        return 0;
    }
}



int BULL2D::GeometryTier::compareBetaSpanOfEdges(const void* ts1, const void* ts2)
{
    EdgeBU2D* edge1 = *(EdgeBU2D**) ts1;
    EdgeBU2D* edge2 = *(EdgeBU2D**) ts2;


    rg_REAL lowerBoundOfEdge1 = edge1->getLowerBoundForValidSimplex();
    rg_REAL lowerBoundOfEdge2 = edge2->getLowerBoundForValidSimplex();
    if ( lowerBoundOfEdge1 > lowerBoundOfEdge2 ) {
        return 1;
    }
    else if ( lowerBoundOfEdge1 < lowerBoundOfEdge2 ) {
        return -1;
    }
    else {
        return 0;
    }
}



int BULL2D::GeometryTier::compareBetaSpanOfVertices(const void* ts1, const void* ts2)
{
    VertexBU2D* vertex1 = *(VertexBU2D**) ts1;
    VertexBU2D* vertex2 = *(VertexBU2D**) ts2;

    rg_REAL lowerBoundOfVtx1 = vertex1->getLowerBoundForValidSimplex();
    rg_REAL lowerBoundOfVtx2 = vertex2->getLowerBoundForValidSimplex();
    if ( lowerBoundOfVtx1 > lowerBoundOfVtx2 ) {
        return 1;
    }
    else if ( lowerBoundOfVtx1 < lowerBoundOfVtx2 ) {
        return -1;
    }
    else {
        return 0;
    }
}


    
void BetaUniverse2D::computeBetaInterval()
{
    rg_REAL lowerBeta = rg_REAL_INFINITY;
	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		VertexBU2D* currVertex = m_vertices.getpEntity();

        if ( currVertex->isVirtual() ) {
            continue;
        }

        rg_REAL lower = currVertex->getLowerBoundForValidSimplex();
        if ( lower < lowerBeta ) {
            lowerBeta = lower;
        }
    }


    rg_REAL upperBeta = -rg_REAL_INFINITY;
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		FaceBU2D* currFace = m_faces.getpEntity();

        if ( currFace->isVirtual() ) {
            continue;
        }

        rg_REAL lower = currFace->getLowerBoundForValidSimplex();
        if ( lower > upperBeta ) {
            upperBeta = lower;
        }
    }

    m_betaInterval.setLowerAndUpperValue(lowerBeta, upperBeta);
}



void BetaUniverse2D::rearrangeSimplexIDs()
{
    rearrangeFaceIDs();
    rearrangeEdgeIDs();
    rearrangeVertexIDs();
}



void BetaUniverse2D::rearrangeFaceIDs()
{
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		m_faces.getpEntity()->setID( -1 );
	}

    rg_INT i=0; 
    for ( i=0; i<m_numFiniteFaces; i++ ) {
        m_sortedFiniteFaces[i]->setID( i );
    }

    rg_INT ID = m_numFiniteFaces;
    FaceBU2D* currFace = rg_NULL;
	m_faces.reset4Loop();
	while( m_faces.setNext4Loop() )  {
		currFace = m_faces.getpEntity();

        if ( currFace->getID() == -1 ) {
            currFace->setID( ID++ );
        }
	}
}



void BetaUniverse2D::rearrangeEdgeIDs()
{
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		m_edges.getpEntity()->setID( -1 );
	}

    rg_INT i=0; 
    for ( i=0; i<m_numFiniteEdges; i++ ) {
        m_sortedFiniteEdges[i]->setID( i );
    }

    rg_INT ID = m_numFiniteEdges;
    EdgeBU2D* currEdge = rg_NULL;
	m_edges.reset4Loop();
	while( m_edges.setNext4Loop() )  {
		currEdge = m_edges.getpEntity();

        if ( currEdge->getID() == -1 ) {
            currEdge->setID( ID++ );
        }
	}
}



void BetaUniverse2D::rearrangeVertexIDs()
{
	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		m_vertices.getpEntity()->setID( -1 );
	}

    rg_INT i=0; 
    for ( i=0; i<m_numFiniteVertices; i++ ) {
        m_sortedFiniteVertices[i]->setID( i );
    }

    rg_INT ID = m_numFiniteVertices;
    VertexBU2D* currVertex = rg_NULL;
	m_vertices.reset4Loop();
	while( m_vertices.setNext4Loop() )  {
		currVertex = m_vertices.getpEntity();

        if ( currVertex->getID() == -1 ) {
            currVertex->setID( ID++ );
        }
	}
}