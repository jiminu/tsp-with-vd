#include "BetaUniverse2D.h"
#ifdef _COMP_TIME_BC_IMAGE
#include <time.h>
#endif
#include <fstream>
using namespace std;
using namespace V::GeometryTier;


/*******************************************************************/
/*   Moved above all the other functions by CYSONG at July16, 2020.*/
/*******************************************************************/
int compareBetaSpanOfFaces(const void* ts1, const void* ts2)
{
    FaceBU2D* face1 = *(FaceBU2D**)ts1;
    FaceBU2D* face2 = *(FaceBU2D**)ts2;


    rg_REAL lowerBoundOfFace1 = face1->getLowerBoundForValidSimplex();
    rg_REAL lowerBoundOfFace2 = face2->getLowerBoundForValidSimplex();
    if (lowerBoundOfFace1 > lowerBoundOfFace2) {
        return 1;
    }
    else if (lowerBoundOfFace1 < lowerBoundOfFace2) {
        return -1;
    }
    else {
        return 0;
    }
}



int compareBetaSpanOfEdges(const void* ts1, const void* ts2)
{
    EdgeBU2D* edge1 = *(EdgeBU2D**)ts1;
    EdgeBU2D* edge2 = *(EdgeBU2D**)ts2;


    rg_REAL lowerBoundOfEdge1 = edge1->getLowerBoundForValidSimplex();
    rg_REAL lowerBoundOfEdge2 = edge2->getLowerBoundForValidSimplex();
    if (lowerBoundOfEdge1 > lowerBoundOfEdge2) {
        return 1;
    }
    else if (lowerBoundOfEdge1 < lowerBoundOfEdge2) {
        return -1;
    }
    else {
        return 0;
    }
}



int compareBetaSpanOfVertices(const void* ts1, const void* ts2)
{
    VertexBU2D* vertex1 = *(VertexBU2D**)ts1;
    VertexBU2D* vertex2 = *(VertexBU2D**)ts2;

    rg_REAL lowerBoundOfVtx1 = vertex1->getLowerBoundForValidSimplex();
    rg_REAL lowerBoundOfVtx2 = vertex2->getLowerBoundForValidSimplex();
    if (lowerBoundOfVtx1 > lowerBoundOfVtx2) {
        return 1;
    }
    else if (lowerBoundOfVtx1 < lowerBoundOfVtx2) {
        return -1;
    }
    else {
        return 0;
    }
}

/**************************************************/









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



void BetaUniverse2D::construct(QuasiTriangulation2D& QT2D)
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

    	rg_BetaSpan betaspan = currVtx->getBetaSpan();


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

    	rg_BetaSpan betaspan = currEdge->getBetaSpan();

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

    	rg_BetaSpan betaspan = currFace->getBetaSpan();

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



void BetaUniverse2D::convert(QuasiTriangulation2D& QT2D)
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





void BetaUniverse2D::createAndInitializeFacesInBU2D( QuasiTriangulation2D& QT2D,    map<FaceSCDS*, FaceBU2D*>&     mappingTableFromFaceInQTToFaceInBU )
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


void BetaUniverse2D::createAndInitializeVerticesInBU2D( QuasiTriangulation2D& QT2D, map<VertexSCDS*, VertexBU2D*>& mappingTableFromVtxInQTToVtxInBU )
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



void BetaUniverse2D::constructForImageProcessing( const vector< vector< vector< int > > >& colorDataSet, const int & lengthOfXDir, const int & lengthOfYDir, const bool& darkerPixelBecomesBiggerDisk )
{
#ifdef _COMP_TIME_BC_IMAGE
    int compTime_creatEntities  = 0;
    int compTime_checkSlope     = 0;
    int compTime_collect_vtx_outside    = 0;
    int compTime_collect_edge_outside   = 0;
    int compTime_createInnerFace_tangentCircle  = 0;
    int compTime_bookkeeping                    = 0;
    int compTime_betaSpan                       = 0;
    int compTime_sortEntities                   = 0;
    int compTime_betaInterval                   = 0;

    clock_t startTime, endTime;
#endif


    rg_REAL  left_line_coordinate = - ((rg_REAL)lengthOfXDir / 2.0) + 0.5;
    rg_REAL  bottom_line_coordinate = -((rg_REAL)lengthOfYDir / 2.0) + 0.5;
    
    vector< vector< Disc* > >       discs;
    vector< vector< VertexBU2D* > > vertices;
    
    discs.resize(    lengthOfXDir );
    vertices.resize( lengthOfXDir );

#ifdef _COMP_TIME_BC_IMAGE
    startTime = clock();
#endif
    
    int ID = 0;
    for ( int i = 0; i < lengthOfXDir; ++i ) {
        discs[i].resize(    lengthOfYDir );
        vertices[i].resize( lengthOfYDir );
    }

    for ( int j = 0; j < lengthOfYDir; ++j ) {
        for ( int i = 0; i < lengthOfXDir; ++i, ++ID ) {
            rg_FLOAT x = left_line_coordinate + (rg_FLOAT)i;
            rg_FLOAT y = bottom_line_coordinate + (rg_FLOAT)j;
            rg_FLOAT r = calculateRadiusOfPixel( colorDataSet[i][j][0], colorDataSet[i][j][1], colorDataSet[i][j][2], darkerPixelBecomesBiggerDisk );

            discs[i][j]     = m_discs.add( Disc(ID, rg_Circle2D(x, y, r) ) );
            vertices[i][j]  = m_vertices.add( VertexBU2D(ID, discs[i][j]) );
        }
    }

    VertexBU2D* infiniteVertex = m_vertices.addHead( VertexBU2D( -1, rg_NULL ) ); // infinite vertex in BetaUniverse2D

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_creatEntities = endTime - startTime;

    startTime = clock();
#endif

    vector< vector< bool > > b_slopeOfEdgeIsTiltedToTheRightUp;
    b_slopeOfEdgeIsTiltedToTheRightUp.resize( lengthOfXDir - 1 );
    for ( int i = 0; i < lengthOfXDir - 1; ++i ) {
        b_slopeOfEdgeIsTiltedToTheRightUp[i].resize( lengthOfYDir - 1 );
        for ( int j = 0; j < lengthOfYDir - 1; ++j ) {
            Disc* disk_leftBottom   = discs[i][j];
            Disc* disk_rightBottom  = discs[i+1][j];
            Disc* disk_leftTop      = discs[i][j+1];
            Disc* disk_rightTop     = discs[i+1][j+1];

            b_slopeOfEdgeIsTiltedToTheRightUp[i][j] = diagonal_is_positive( disk_leftBottom, disk_rightBottom, disk_leftTop, disk_rightTop );
        }
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_checkSlope = endTime - startTime;
    
    startTime = clock();
#endif

    // vertices on outside
    int numVerticesOnOutside = 2 * (lengthOfXDir + lengthOfYDir) - 4;
    vector< VertexBU2D* > verticesOnOutside;
    verticesOnOutside.resize( numVerticesOnOutside );
    
    int i_vtx_outside = 0;
    for ( int i = 0; i < lengthOfXDir; ++i, ++i_vtx_outside ) {
        verticesOnOutside[i_vtx_outside] = vertices[i][0];
    }
    
    for ( int j = 1; j < lengthOfYDir - 1; ++j, ++i_vtx_outside ) {
        verticesOnOutside[i_vtx_outside] = vertices[lengthOfXDir-1][j];
    }

    for ( int i = lengthOfXDir - 1; i >= 0; --i, ++i_vtx_outside ) {
        verticesOnOutside[i_vtx_outside] = vertices[i][lengthOfYDir-1];
    }

    for ( int j = lengthOfYDir - 2; j > 0; --j, ++i_vtx_outside ) {
        verticesOnOutside[i_vtx_outside] = vertices[0][j];
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_collect_vtx_outside = endTime - startTime;

    startTime = clock();
#endif

    // horizontal edges
    int ID_edge = 0;

    vector< vector< EdgeBU2D* > > horizontalEdges;
    horizontalEdges.resize( lengthOfXDir - 1 );
    for ( int i = 0; i < lengthOfXDir - 1; ++i ) {
        horizontalEdges[i].resize( lengthOfYDir );
        for ( int j = 0; j < lengthOfYDir; ++j, ++ID_edge ) {
            horizontalEdges[i][j] = m_edges.add( EdgeBU2D(ID_edge) );
        }
    }

    // vertical edges
    vector< vector< EdgeBU2D* > > verticalEdges;
    verticalEdges.resize( lengthOfXDir );
    for ( int i = 0; i < lengthOfXDir; ++i ) {
        verticalEdges[i].resize( lengthOfYDir - 1 );
        for ( int j = 0; j < lengthOfYDir - 1; ++j, ++ID_edge ) {
            verticalEdges[i][j] = m_edges.add( EdgeBU2D(ID_edge) );
        }
    }

    // slope edges
    vector< vector< EdgeBU2D* > > slopeEdges;
    slopeEdges.resize( lengthOfXDir - 1 );
    for ( int i = 0; i < lengthOfXDir - 1; ++i ) {
        slopeEdges[i].resize( lengthOfYDir - 1 );
        for ( int j = 0; j < lengthOfYDir - 1; ++j, ++ID_edge ) {
            slopeEdges[i][j] = m_edges.add( EdgeBU2D(ID_edge) );
        }
    }

    // unbounded edges
    int numUnboundedEdges = 2*(lengthOfXDir + lengthOfYDir) - 4;
    vector< EdgeBU2D* > unboundedEdges;
    unboundedEdges.resize( numUnboundedEdges + 1 );
    for ( int i = 0; i < numUnboundedEdges; ++i, ++ID_edge ) {
        unboundedEdges[i] = m_edges.add( EdgeBU2D(ID_edge) );
    }
    unboundedEdges[numUnboundedEdges] = unboundedEdges[0];

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_creatEntities = compTime_creatEntities + endTime - startTime;

    startTime = clock();
#endif


    // edges on outside
    int numEdgesOnOutside = numVerticesOnOutside;
    vector< EdgeBU2D* > edgesOnOutside;
    edgesOnOutside.resize( numEdgesOnOutside );

    int i_edge_outside = 0;
    for ( int i = 0; i < lengthOfXDir - 1; ++i, ++i_edge_outside ) {
        edgesOnOutside[i_edge_outside] = horizontalEdges[i][0];
    }

    for ( int j = 0; j < lengthOfYDir - 1; ++j, ++i_edge_outside ) {
        edgesOnOutside[i_edge_outside] = verticalEdges[lengthOfXDir-1][j];
    }

    for ( int i = lengthOfXDir - 2; i >= 0; --i, ++i_edge_outside ) {
        edgesOnOutside[i_edge_outside] = horizontalEdges[i][lengthOfYDir-1];
    }

    for ( int j = lengthOfYDir - 2; j >= 0; --j, ++i_edge_outside ) {
        edgesOnOutside[i_edge_outside] = verticalEdges[0][j];
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_collect_edge_outside = endTime - startTime;

    startTime = clock();
#endif


    // inner faces
    int ID_face = 0;
    vector< vector< vector< FaceBU2D* > > > innerFaces;
    innerFaces.resize( lengthOfXDir - 1 );
    for ( int i = 0; i < lengthOfXDir - 1; ++i ) {
        innerFaces[i].resize( lengthOfYDir - 1 );
        for ( int j = 0; j < lengthOfYDir - 1; ++j, ++ID_face ) {
            innerFaces[i][j].resize( 2 );

            //rg_Circle2D tangentCircle_lowerFace;
            //rg_Circle2D tangentCircle_upperFace;
            //if ( b_slopeOfEdgeIsTiltedToTheRightUp[i][j] ) {
            //    rg_Circle2D::makeCircumcircle( discs[i][j]->getGeometry(), discs[i+1][j+1]->getGeometry(), discs[i+1][j]->getGeometry(), tangentCircle_lowerFace, rg_Circle2D() );
            //    rg_Circle2D::makeCircumcircle( discs[i][j]->getGeometry(), discs[i][j+1]->getGeometry(), discs[i+1][j+1]->getGeometry(), tangentCircle_upperFace, rg_Circle2D() );
            //}
            //else {
            //    rg_Circle2D::makeCircumcircle( discs[i][j]->getGeometry(), discs[i][j+1]->getGeometry(), discs[i+1][j]->getGeometry(), tangentCircle_lowerFace, rg_Circle2D() );
            //    rg_Circle2D::makeCircumcircle( discs[i][j+1]->getGeometry(), discs[i+1][j+1]->getGeometry(), discs[i+1][j]->getGeometry(), tangentCircle_upperFace, rg_Circle2D() );
            //}
            //innerFaces[i][j][0] = m_faces.add( FaceBU2D(ID_face, tangentCircle_lowerFace) );
            //++ID_face;
            //innerFaces[i][j][1] = m_faces.add( FaceBU2D(ID_face, tangentCircle_upperFace) );

            double r_leftBottom     = discs[i][j]->getGeometry().getRadius();
            double r_rightBottom    = discs[i+1][j]->getGeometry().getRadius();
            double r_leftTop        = discs[i][j+1]->getGeometry().getRadius();
            double r_rightTop       = discs[i+1][j+1]->getGeometry().getRadius();

            rg_Circle2D tangentCircle_lowerFace;
            rg_Circle2D tangentCircle_upperFace;
            if ( b_slopeOfEdgeIsTiltedToTheRightUp[i][j] ) {
                tangentCircle_lowerFace.setRadius( calculateRadiusOfTangentCircle_quadraticFormula(r_rightBottom, r_leftBottom, r_rightTop) );
                tangentCircle_upperFace.setRadius( calculateRadiusOfTangentCircle_quadraticFormula(r_leftTop, r_rightTop, r_leftBottom) );
            }
            else {
                tangentCircle_lowerFace.setRadius( calculateRadiusOfTangentCircle_quadraticFormula(r_leftBottom, r_leftTop, r_rightBottom) );
                tangentCircle_upperFace.setRadius( calculateRadiusOfTangentCircle_quadraticFormula(r_rightTop, r_rightBottom, r_leftTop) );
            }
            innerFaces[i][j][0] = m_faces.add( FaceBU2D(ID_face, tangentCircle_lowerFace) );
            ++ID_face;
            innerFaces[i][j][1] = m_faces.add( FaceBU2D(ID_face, tangentCircle_upperFace) );
        }
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_createInnerFace_tangentCircle = endTime - startTime;

    startTime = clock();
#endif

    // infinite faces
    int numInfiniteFaces = numUnboundedEdges;
    vector< FaceBU2D* > infiniteFaces;
    infiniteFaces.resize( numInfiniteFaces );
    for ( int i = 0; i < numInfiniteFaces; ++i, ++ID_face ) {
        infiniteFaces[i] = m_faces.add( FaceBU2D(ID_face, rg_Circle2D()) );
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_creatEntities = compTime_creatEntities + endTime - startTime;

    startTime = clock();
#endif

    // BOOKKEEPING
    // vertices to edges
    for ( int i = 0; i < lengthOfXDir - 1; ++i ) {
        for ( int j = 0; j < lengthOfYDir; ++j ) {
            vertices[i][j]->setFirstEdge( horizontalEdges[i][j] );
        }
    }

    for ( int j = 0; j < lengthOfYDir; ++j ) {
        vertices[lengthOfXDir-1][j]->setFirstEdge( horizontalEdges[lengthOfXDir-2][j] );
    }

    infiniteVertex->setFirstEdge( unboundedEdges[0] );

    // edge to vertices
    // horizontal edges
    for ( int i = 0; i < lengthOfXDir - 1; ++i ) {
        for ( int j = 0; j < lengthOfYDir; ++j ) {
            horizontalEdges[i][j]->setStartVertex( vertices[i][j] );
            horizontalEdges[i][j]->setEndVertex(   vertices[i+1][j] );
        }
    }

    // vertical edges
    for ( int i = 0; i < lengthOfXDir; ++i ) {
        for ( int j = 0; j < lengthOfYDir - 1; ++j ) {
            verticalEdges[i][j]->setStartVertex( vertices[i][j] );
            verticalEdges[i][j]->setEndVertex(   vertices[i][j+1] );
        }
    }

    // slope edges
    for ( int i = 0; i < lengthOfXDir - 1; ++i ) {
        for ( int j = 0; j < lengthOfYDir - 1; ++j ) {
            if ( b_slopeOfEdgeIsTiltedToTheRightUp[i][j] ) {
                slopeEdges[i][j]->setStartVertex( vertices[i][j] );
                slopeEdges[i][j]->setEndVertex(   vertices[i+1][j+1] );
            }
            else {
                slopeEdges[i][j]->setStartVertex( vertices[i][j+1] );
                slopeEdges[i][j]->setEndVertex(   vertices[i+1][j] );
            }
        }
    }

    // unbounded edges
    for ( int i = 0; i < numUnboundedEdges; ++i ) {
        unboundedEdges[i]->setStartVertex( verticesOnOutside[i] );
        unboundedEdges[i]->setEndVertex(   infiniteVertex );
    }

    // inner grid cell
    for ( int i = 0; i < lengthOfXDir - 1; ++i ) {
        for ( int j = 0; j < lengthOfYDir - 1; ++j ) {
            EdgeBU2D* topEdge       = horizontalEdges[i][j+1];
            EdgeBU2D* bottomEdge    = horizontalEdges[i][j];
            EdgeBU2D* leftEdge      = verticalEdges[i][j];
            EdgeBU2D* rightEdge     = verticalEdges[i+1][j];
            EdgeBU2D* slopeEdge     = slopeEdges[i][j];

            FaceBU2D* lowerFace     = innerFaces[i][j][0];
            FaceBU2D* upperFace     = innerFaces[i][j][1];
            
            upperFace->setFirstEdge( topEdge );
            lowerFace->setFirstEdge( bottomEdge );

            if ( b_slopeOfEdgeIsTiltedToTheRightUp[i][j] ) {
                topEdge->setRightHand(      slopeEdge );
                topEdge->setRightLeg(       leftEdge );
                topEdge->setRightFace(      upperFace );

                bottomEdge->setLeftHand(    rightEdge );
                bottomEdge->setLeftLeg(     slopeEdge );
                bottomEdge->setLeftFace(    lowerFace );

                leftEdge->setRightHand(     topEdge );
                leftEdge->setRightLeg(      slopeEdge );
                leftEdge->setRightFace(     upperFace );

                rightEdge->setLeftHand(     slopeEdge );
                rightEdge->setLeftLeg(      bottomEdge );
                rightEdge->setLeftFace(     lowerFace );

                slopeEdge->setRightHand(    rightEdge );
                slopeEdge->setRightLeg(     bottomEdge );
                slopeEdge->setLeftHand(     topEdge );
                slopeEdge->setLeftLeg(      leftEdge );
                slopeEdge->setRightFace(    lowerFace );
                slopeEdge->setLeftFace(     upperFace );
            }
            else {
                topEdge->setRightHand(      rightEdge );
                topEdge->setRightLeg(       slopeEdge );
                topEdge->setRightFace(      upperFace );

                bottomEdge->setLeftHand(    slopeEdge );
                bottomEdge->setLeftLeg(     leftEdge );
                bottomEdge->setLeftFace(    lowerFace );

                leftEdge->setRightHand(     slopeEdge );
                leftEdge->setRightLeg(      bottomEdge );
                leftEdge->setRightFace(     lowerFace );

                rightEdge->setLeftHand(     topEdge );
                rightEdge->setLeftLeg(      slopeEdge );
                rightEdge->setLeftFace(     upperFace );

                slopeEdge->setRightHand(    bottomEdge );
                slopeEdge->setRightLeg(     leftEdge );
                slopeEdge->setLeftHand(     rightEdge );
                slopeEdge->setLeftLeg(      topEdge );
                slopeEdge->setRightFace(    lowerFace );
                slopeEdge->setLeftFace(     upperFace );
            }
        }
    }

    // infinite faces
    for ( int i = 0; i < numInfiniteFaces; ++i ) {
        EdgeBU2D* prevEdge_CCW  = unboundedEdges[i];
        EdgeBU2D* nextEdge_CCW  = unboundedEdges[i+1];
        EdgeBU2D* edge_finite   = edgesOnOutside[i];

        FaceBU2D* infiniteFace  = infiniteFaces[i];

        infiniteFace->setFirstEdge( edge_finite );

        prevEdge_CCW->setLeftHand( nextEdge_CCW );
        prevEdge_CCW->setLeftLeg(  edge_finite );
        prevEdge_CCW->setLeftFace( infiniteFace );

        nextEdge_CCW->setRightHand( prevEdge_CCW );
        nextEdge_CCW->setRightLeg(  edge_finite );
        nextEdge_CCW->setRightFace( infiniteFace );

        if ( edge_finite->getLeftFace() == rg_NULL ) {
            edge_finite->setLeftHand( prevEdge_CCW );
            edge_finite->setLeftLeg(  nextEdge_CCW );
            edge_finite->setLeftFace( infiniteFace );
        }
        else {
            edge_finite->setRightLeg(  prevEdge_CCW );
            edge_finite->setRightHand( nextEdge_CCW );
            edge_finite->setRightFace( infiniteFace );
        }
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_bookkeeping = endTime - startTime;

    startTime = clock();
#endif


    //computeBetaSpan();
    computeBetaSpanForImageProcessing(lengthOfXDir, lengthOfYDir);



#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_betaSpan = endTime - startTime;

    startTime = clock();
#endif


#ifdef _WRITE_IMAGE_TOPOLOGY
    ofstream fout("imageTopology.txt");
    fout << lengthOfXDir << " X " << lengthOfYDir << endl;
    for ( int i = 0; i < lengthOfXDir; ++i ) {
        for ( int j = 0; j < lengthOfYDir; ++j ) {
            fout << colorDataSet[i][j][0] << "\t" << colorDataSet[i][j][1] << "\t" << colorDataSet[i][j][2] << endl;
        }
    }
    fout << endl;

    fout << "The number of vertices: " << m_vertices.getSize() << endl;
    m_vertices.reset4Loop();
    while ( m_vertices.setNext4Loop() ) {
        VertexBU2D* currVertex = m_vertices.getpEntity();

        fout << currVertex->getID() << "\t" << currVertex->getFirstEdge()->getID() << "\t";
        if ( currVertex->isVirtual() ) {
            fout << "INFINITE VERTEX" << endl;
            continue;
        }

        fout << currVertex->getCircle().getCenterPt().getX() << "\t" << 
            currVertex->getCircle().getCenterPt().getY() << "\t" << 
            currVertex->getCircle().getRadius() << endl;
    }
    fout << endl;

    m_edges.reset4Loop();
    while ( m_edges.setNext4Loop() ) {
        EdgeBU2D* currEdge = m_edges.getpEntity();
        fout << currEdge->getID() << "\t" << 
            currEdge->getStartVertex()->getID() << "\t" << 
            currEdge->getEndVertex()->getID() << "\t" << 
            currEdge->getRightHand()->getID() << "\t" << 
            currEdge->getRightLeg()->getID() << "\t" << 
            currEdge->getLeftHand()->getID() << "\t" << 
            currEdge->getLeftLeg()->getID() << "\t" << 
            currEdge->getRightFace()->getID() << "\t" << 
            currEdge->getLeftFace()->getID() << endl;
    }
    fout << endl;

    m_faces.reset4Loop();
    while ( m_faces.setNext4Loop() ) {
        FaceBU2D* currFace = m_faces.getpEntity();
        fout << currFace->getID() << "\t" << 
            currFace->getFirstEdge()->getID() << endl; //<< currFace->getBetaSpan().getBetaInterval(1) << endl;
    }

    fout.close();     
#endif

    generateSortedFiniteSimplexes();

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_sortEntities = endTime - startTime;

    startTime = clock();
#endif

    computeBetaInterval();

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_betaInterval = endTime - startTime;

    ofstream fout_compTime("compTimeForBC_construction.txt");
    fout_compTime << "Entity creation: " << compTime_creatEntities << endl;
    fout_compTime << "Inner face creation with computing tangent circle: " << compTime_createInnerFace_tangentCircle << endl;
    fout_compTime << "Checking slope with 4 disks: " << compTime_checkSlope << endl;
    fout_compTime << "Collecting outside vertices: " << compTime_collect_vtx_outside << endl;
    fout_compTime << "Collecting outside edges: " << compTime_collect_edge_outside << endl;
    fout_compTime << "Bookkeeping: " << compTime_bookkeeping << endl;
    fout_compTime << "Computing betaspan: " << compTime_betaSpan << endl;
    fout_compTime << "Sorting entities with betaspan: " << compTime_sortEntities << endl;
    fout_compTime << "Computing betainterval: " << compTime_betaInterval << endl;

    fout_compTime.close();
#endif

}


rg_FLOAT BetaUniverse2D::calculateRadiusOfPixel( const int & r, const int & g, const int & b, const bool& darkerPixelBecomesBiggerDisk )
{
    rg_FLOAT    red     = (rg_FLOAT)r / 255.0f;
    rg_FLOAT    green   = (rg_FLOAT)g / 255.0f;
    rg_FLOAT    blue    = (rg_FLOAT)b / 255.0f;
        
    rg_FLOAT    gray    = 0.2989 * red + 0.5870 * green + 0.1140 * blue;

    return darkerPixelBecomesBiggerDisk ? gray / 2.0 : ( 1.0 - gray ) / 2.0;
}


void BetaUniverse2D::constructForImageProcessing(const vector< vector< unsigned char > >& colorDataSet, const int& lengthOfXDir, const int& lengthOfYDir, const bool& darkerPixelBecomesBiggerDisk)
{
#ifdef _COMP_TIME_BC_IMAGE
    int compTime_creatEntities = 0;
    int compTime_checkSlope = 0;
    int compTime_collect_vtx_outside = 0;
    int compTime_collect_edge_outside = 0;
    int compTime_createInnerFace_tangentCircle = 0;
    int compTime_bookkeeping = 0;
    int compTime_betaSpan = 0;
    int compTime_sortEntities = 0;
    int compTime_betaInterval = 0;

    clock_t startTime, endTime;
#endif


    rg_REAL  left_line_coordinate = -((rg_REAL)lengthOfXDir / 2.0) + 0.5;
    rg_REAL  bottom_line_coordinate = -((rg_REAL)lengthOfYDir / 2.0) + 0.5;

    vector< vector< Disc* > >       discs;
    vector< vector< VertexBU2D* > > vertices;

    discs.resize(lengthOfXDir);
    vertices.resize(lengthOfXDir);

#ifdef _COMP_TIME_BC_IMAGE
    startTime = clock();
#endif

    int ID = 0;
    for (int i = 0; i < lengthOfXDir; ++i) {
        discs[i].resize(lengthOfYDir);
        vertices[i].resize(lengthOfYDir);
    }

    for (int j = 0; j < lengthOfYDir; ++j) {
        for (int i = 0; i < lengthOfXDir; ++i, ++ID) {
            rg_FLOAT x = left_line_coordinate + (rg_FLOAT)i;
            rg_FLOAT y = bottom_line_coordinate + (rg_FLOAT)j;
            rg_FLOAT r = calculateRadiusOfPixel(colorDataSet[i][j], darkerPixelBecomesBiggerDisk);

            discs[i][j] = m_discs.add(Disc(ID, rg_Circle2D(x, y, r)));
            vertices[i][j] = m_vertices.add(VertexBU2D(ID, discs[i][j]));
        }
    }

    VertexBU2D* infiniteVertex = m_vertices.addHead(VertexBU2D(-1, rg_NULL)); // infinite vertex in BetaUniverse2D

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_creatEntities = endTime - startTime;

    startTime = clock();
#endif

    vector< vector< bool > > b_slopeOfEdgeIsTiltedToTheRightUp;
    b_slopeOfEdgeIsTiltedToTheRightUp.resize(lengthOfXDir - 1);
    for (int i = 0; i < lengthOfXDir - 1; ++i) {
        b_slopeOfEdgeIsTiltedToTheRightUp[i].resize(lengthOfYDir - 1);
        for (int j = 0; j < lengthOfYDir - 1; ++j) {
            Disc* disk_leftBottom = discs[i][j];
            Disc* disk_rightBottom = discs[i + 1][j];
            Disc* disk_leftTop = discs[i][j + 1];
            Disc* disk_rightTop = discs[i + 1][j + 1];

            b_slopeOfEdgeIsTiltedToTheRightUp[i][j] = diagonal_is_positive(disk_leftBottom, disk_rightBottom, disk_leftTop, disk_rightTop);
        }
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_checkSlope = endTime - startTime;

    startTime = clock();
#endif

    // vertices on outside
    int numVerticesOnOutside = 2 * (lengthOfXDir + lengthOfYDir) - 4;
    vector< VertexBU2D* > verticesOnOutside;
    verticesOnOutside.resize(numVerticesOnOutside);

    int i_vtx_outside = 0;
    for (int i = 0; i < lengthOfXDir; ++i, ++i_vtx_outside) {
        verticesOnOutside[i_vtx_outside] = vertices[i][0];
    }

    for (int j = 1; j < lengthOfYDir - 1; ++j, ++i_vtx_outside) {
        verticesOnOutside[i_vtx_outside] = vertices[lengthOfXDir - 1][j];
    }

    for (int i = lengthOfXDir - 1; i >= 0; --i, ++i_vtx_outside) {
        verticesOnOutside[i_vtx_outside] = vertices[i][lengthOfYDir - 1];
    }

    for (int j = lengthOfYDir - 2; j > 0; --j, ++i_vtx_outside) {
        verticesOnOutside[i_vtx_outside] = vertices[0][j];
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_collect_vtx_outside = endTime - startTime;

    startTime = clock();
#endif

    // horizontal edges
    int ID_edge = 0;

    vector< vector< EdgeBU2D* > > horizontalEdges;
    horizontalEdges.resize(lengthOfXDir - 1);
    for (int i = 0; i < lengthOfXDir - 1; ++i) {
        horizontalEdges[i].resize(lengthOfYDir);
        for (int j = 0; j < lengthOfYDir; ++j, ++ID_edge) {
            horizontalEdges[i][j] = m_edges.add(EdgeBU2D(ID_edge));
        }
    }

    // vertical edges
    vector< vector< EdgeBU2D* > > verticalEdges;
    verticalEdges.resize(lengthOfXDir);
    for (int i = 0; i < lengthOfXDir; ++i) {
        verticalEdges[i].resize(lengthOfYDir - 1);
        for (int j = 0; j < lengthOfYDir - 1; ++j, ++ID_edge) {
            verticalEdges[i][j] = m_edges.add(EdgeBU2D(ID_edge));
        }
    }

    // slope edges
    vector< vector< EdgeBU2D* > > slopeEdges;
    slopeEdges.resize(lengthOfXDir - 1);
    for (int i = 0; i < lengthOfXDir - 1; ++i) {
        slopeEdges[i].resize(lengthOfYDir - 1);
        for (int j = 0; j < lengthOfYDir - 1; ++j, ++ID_edge) {
            slopeEdges[i][j] = m_edges.add(EdgeBU2D(ID_edge));
        }
    }

    // unbounded edges
    int numUnboundedEdges = 2 * (lengthOfXDir + lengthOfYDir) - 4;
    vector< EdgeBU2D* > unboundedEdges;
    unboundedEdges.resize(numUnboundedEdges + 1);
    for (int i = 0; i < numUnboundedEdges; ++i, ++ID_edge) {
        unboundedEdges[i] = m_edges.add(EdgeBU2D(ID_edge));
    }
    unboundedEdges[numUnboundedEdges] = unboundedEdges[0];

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_creatEntities = compTime_creatEntities + endTime - startTime;

    startTime = clock();
#endif


    // edges on outside
    int numEdgesOnOutside = numVerticesOnOutside;
    vector< EdgeBU2D* > edgesOnOutside;
    edgesOnOutside.resize(numEdgesOnOutside);

    int i_edge_outside = 0;
    for (int i = 0; i < lengthOfXDir - 1; ++i, ++i_edge_outside) {
        edgesOnOutside[i_edge_outside] = horizontalEdges[i][0];
    }

    for (int j = 0; j < lengthOfYDir - 1; ++j, ++i_edge_outside) {
        edgesOnOutside[i_edge_outside] = verticalEdges[lengthOfXDir - 1][j];
    }

    for (int i = lengthOfXDir - 2; i >= 0; --i, ++i_edge_outside) {
        edgesOnOutside[i_edge_outside] = horizontalEdges[i][lengthOfYDir - 1];
    }

    for (int j = lengthOfYDir - 2; j >= 0; --j, ++i_edge_outside) {
        edgesOnOutside[i_edge_outside] = verticalEdges[0][j];
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_collect_edge_outside = endTime - startTime;

    startTime = clock();
#endif


    // inner faces
    int ID_face = 0;
    vector< vector< vector< FaceBU2D* > > > innerFaces;
    innerFaces.resize(lengthOfXDir - 1);
    for (int i = 0; i < lengthOfXDir - 1; ++i) {
        innerFaces[i].resize(lengthOfYDir - 1);
        for (int j = 0; j < lengthOfYDir - 1; ++j, ++ID_face) {
            innerFaces[i][j].resize(2);

            //rg_Circle2D tangentCircle_lowerFace;
            //rg_Circle2D tangentCircle_upperFace;
            //if ( b_slopeOfEdgeIsTiltedToTheRightUp[i][j] ) {
            //    rg_Circle2D::makeCircumcircle( discs[i][j]->getGeometry(), discs[i+1][j+1]->getGeometry(), discs[i+1][j]->getGeometry(), tangentCircle_lowerFace, rg_Circle2D() );
            //    rg_Circle2D::makeCircumcircle( discs[i][j]->getGeometry(), discs[i][j+1]->getGeometry(), discs[i+1][j+1]->getGeometry(), tangentCircle_upperFace, rg_Circle2D() );
            //}
            //else {
            //    rg_Circle2D::makeCircumcircle( discs[i][j]->getGeometry(), discs[i][j+1]->getGeometry(), discs[i+1][j]->getGeometry(), tangentCircle_lowerFace, rg_Circle2D() );
            //    rg_Circle2D::makeCircumcircle( discs[i][j+1]->getGeometry(), discs[i+1][j+1]->getGeometry(), discs[i+1][j]->getGeometry(), tangentCircle_upperFace, rg_Circle2D() );
            //}
            //innerFaces[i][j][0] = m_faces.add( FaceBU2D(ID_face, tangentCircle_lowerFace) );
            //++ID_face;
            //innerFaces[i][j][1] = m_faces.add( FaceBU2D(ID_face, tangentCircle_upperFace) );

            double r_leftBottom = discs[i][j]->getGeometry().getRadius();
            double r_rightBottom = discs[i + 1][j]->getGeometry().getRadius();
            double r_leftTop = discs[i][j + 1]->getGeometry().getRadius();
            double r_rightTop = discs[i + 1][j + 1]->getGeometry().getRadius();

            rg_Circle2D tangentCircle_lowerFace;
            rg_Circle2D tangentCircle_upperFace;
            if (b_slopeOfEdgeIsTiltedToTheRightUp[i][j]) {
                tangentCircle_lowerFace.setRadius(calculateRadiusOfTangentCircle_quadraticFormula(r_rightBottom, r_leftBottom, r_rightTop));
                tangentCircle_upperFace.setRadius(calculateRadiusOfTangentCircle_quadraticFormula(r_leftTop, r_rightTop, r_leftBottom));
            }
            else {
                tangentCircle_lowerFace.setRadius(calculateRadiusOfTangentCircle_quadraticFormula(r_leftBottom, r_leftTop, r_rightBottom));
                tangentCircle_upperFace.setRadius(calculateRadiusOfTangentCircle_quadraticFormula(r_rightTop, r_rightBottom, r_leftTop));
            }
            innerFaces[i][j][0] = m_faces.add(FaceBU2D(ID_face, tangentCircle_lowerFace));
            ++ID_face;
            innerFaces[i][j][1] = m_faces.add(FaceBU2D(ID_face, tangentCircle_upperFace));
        }
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_createInnerFace_tangentCircle = endTime - startTime;

    startTime = clock();
#endif

    // infinite faces
    int numInfiniteFaces = numUnboundedEdges;
    vector< FaceBU2D* > infiniteFaces;
    infiniteFaces.resize(numInfiniteFaces);
    for (int i = 0; i < numInfiniteFaces; ++i, ++ID_face) {
        infiniteFaces[i] = m_faces.add(FaceBU2D(ID_face, rg_Circle2D()));
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_creatEntities = compTime_creatEntities + endTime - startTime;

    startTime = clock();
#endif

    // BOOKKEEPING
    // vertices to edges
    for (int i = 0; i < lengthOfXDir - 1; ++i) {
        for (int j = 0; j < lengthOfYDir; ++j) {
            vertices[i][j]->setFirstEdge(horizontalEdges[i][j]);
        }
    }

    for (int j = 0; j < lengthOfYDir; ++j) {
        vertices[lengthOfXDir - 1][j]->setFirstEdge(horizontalEdges[lengthOfXDir - 2][j]);
    }

    infiniteVertex->setFirstEdge(unboundedEdges[0]);

    // edge to vertices
    // horizontal edges
    for (int i = 0; i < lengthOfXDir - 1; ++i) {
        for (int j = 0; j < lengthOfYDir; ++j) {
            horizontalEdges[i][j]->setStartVertex(vertices[i][j]);
            horizontalEdges[i][j]->setEndVertex(vertices[i + 1][j]);
        }
    }

    // vertical edges
    for (int i = 0; i < lengthOfXDir; ++i) {
        for (int j = 0; j < lengthOfYDir - 1; ++j) {
            verticalEdges[i][j]->setStartVertex(vertices[i][j]);
            verticalEdges[i][j]->setEndVertex(vertices[i][j + 1]);
        }
    }

    // slope edges
    for (int i = 0; i < lengthOfXDir - 1; ++i) {
        for (int j = 0; j < lengthOfYDir - 1; ++j) {
            if (b_slopeOfEdgeIsTiltedToTheRightUp[i][j]) {
                slopeEdges[i][j]->setStartVertex(vertices[i][j]);
                slopeEdges[i][j]->setEndVertex(vertices[i + 1][j + 1]);
            }
            else {
                slopeEdges[i][j]->setStartVertex(vertices[i][j + 1]);
                slopeEdges[i][j]->setEndVertex(vertices[i + 1][j]);
            }
        }
    }

    // unbounded edges
    for (int i = 0; i < numUnboundedEdges; ++i) {
        unboundedEdges[i]->setStartVertex(verticesOnOutside[i]);
        unboundedEdges[i]->setEndVertex(infiniteVertex);
    }

    // inner grid cell
    for (int i = 0; i < lengthOfXDir - 1; ++i) {
        for (int j = 0; j < lengthOfYDir - 1; ++j) {
            EdgeBU2D* topEdge = horizontalEdges[i][j + 1];
            EdgeBU2D* bottomEdge = horizontalEdges[i][j];
            EdgeBU2D* leftEdge = verticalEdges[i][j];
            EdgeBU2D* rightEdge = verticalEdges[i + 1][j];
            EdgeBU2D* slopeEdge = slopeEdges[i][j];

            FaceBU2D* lowerFace = innerFaces[i][j][0];
            FaceBU2D* upperFace = innerFaces[i][j][1];

            upperFace->setFirstEdge(topEdge);
            lowerFace->setFirstEdge(bottomEdge);

            if (b_slopeOfEdgeIsTiltedToTheRightUp[i][j]) {
                topEdge->setRightHand(slopeEdge);
                topEdge->setRightLeg(leftEdge);
                topEdge->setRightFace(upperFace);

                bottomEdge->setLeftHand(rightEdge);
                bottomEdge->setLeftLeg(slopeEdge);
                bottomEdge->setLeftFace(lowerFace);

                leftEdge->setRightHand(topEdge);
                leftEdge->setRightLeg(slopeEdge);
                leftEdge->setRightFace(upperFace);

                rightEdge->setLeftHand(slopeEdge);
                rightEdge->setLeftLeg(bottomEdge);
                rightEdge->setLeftFace(lowerFace);

                slopeEdge->setRightHand(rightEdge);
                slopeEdge->setRightLeg(bottomEdge);
                slopeEdge->setLeftHand(topEdge);
                slopeEdge->setLeftLeg(leftEdge);
                slopeEdge->setRightFace(lowerFace);
                slopeEdge->setLeftFace(upperFace);
            }
            else {
                topEdge->setRightHand(rightEdge);
                topEdge->setRightLeg(slopeEdge);
                topEdge->setRightFace(upperFace);

                bottomEdge->setLeftHand(slopeEdge);
                bottomEdge->setLeftLeg(leftEdge);
                bottomEdge->setLeftFace(lowerFace);

                leftEdge->setRightHand(slopeEdge);
                leftEdge->setRightLeg(bottomEdge);
                leftEdge->setRightFace(lowerFace);

                rightEdge->setLeftHand(topEdge);
                rightEdge->setLeftLeg(slopeEdge);
                rightEdge->setLeftFace(upperFace);

                slopeEdge->setRightHand(bottomEdge);
                slopeEdge->setRightLeg(leftEdge);
                slopeEdge->setLeftHand(rightEdge);
                slopeEdge->setLeftLeg(topEdge);
                slopeEdge->setRightFace(lowerFace);
                slopeEdge->setLeftFace(upperFace);
            }
        }
    }

    // infinite faces
    for (int i = 0; i < numInfiniteFaces; ++i) {
        EdgeBU2D* prevEdge_CCW = unboundedEdges[i];
        EdgeBU2D* nextEdge_CCW = unboundedEdges[i + 1];
        EdgeBU2D* edge_finite = edgesOnOutside[i];

        FaceBU2D* infiniteFace = infiniteFaces[i];

        infiniteFace->setFirstEdge(edge_finite);

        prevEdge_CCW->setLeftHand(nextEdge_CCW);
        prevEdge_CCW->setLeftLeg(edge_finite);
        prevEdge_CCW->setLeftFace(infiniteFace);

        nextEdge_CCW->setRightHand(prevEdge_CCW);
        nextEdge_CCW->setRightLeg(edge_finite);
        nextEdge_CCW->setRightFace(infiniteFace);

        if (edge_finite->getLeftFace() == rg_NULL) {
            edge_finite->setLeftHand(prevEdge_CCW);
            edge_finite->setLeftLeg(nextEdge_CCW);
            edge_finite->setLeftFace(infiniteFace);
        }
        else {
            edge_finite->setRightLeg(prevEdge_CCW);
            edge_finite->setRightHand(nextEdge_CCW);
            edge_finite->setRightFace(infiniteFace);
        }
    }

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_bookkeeping = endTime - startTime;

    startTime = clock();
#endif


    //computeBetaSpan();
    computeBetaSpanForImageProcessing(lengthOfXDir, lengthOfYDir);



#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_betaSpan = endTime - startTime;

    startTime = clock();
#endif


#ifdef _WRITE_IMAGE_TOPOLOGY
    ofstream fout("imageTopology.txt");
    fout << lengthOfXDir << " X " << lengthOfYDir << endl;
    for (int i = 0; i < lengthOfXDir; ++i) {
        for (int j = 0; j < lengthOfYDir; ++j) {
            fout << colorDataSet[i][j][0] << "\t" << colorDataSet[i][j][1] << "\t" << colorDataSet[i][j][2] << endl;
        }
    }
    fout << endl;

    fout << "The number of vertices: " << m_vertices.getSize() << endl;
    m_vertices.reset4Loop();
    while (m_vertices.setNext4Loop()) {
        VertexBU2D* currVertex = m_vertices.getpEntity();

        fout << currVertex->getID() << "\t" << currVertex->getFirstEdge()->getID() << "\t";
        if (currVertex->isVirtual()) {
            fout << "INFINITE VERTEX" << endl;
            continue;
        }

        fout << currVertex->getCircle().getCenterPt().getX() << "\t" <<
            currVertex->getCircle().getCenterPt().getY() << "\t" <<
            currVertex->getCircle().getRadius() << endl;
    }
    fout << endl;

    m_edges.reset4Loop();
    while (m_edges.setNext4Loop()) {
        EdgeBU2D* currEdge = m_edges.getpEntity();
        fout << currEdge->getID() << "\t" <<
            currEdge->getStartVertex()->getID() << "\t" <<
            currEdge->getEndVertex()->getID() << "\t" <<
            currEdge->getRightHand()->getID() << "\t" <<
            currEdge->getRightLeg()->getID() << "\t" <<
            currEdge->getLeftHand()->getID() << "\t" <<
            currEdge->getLeftLeg()->getID() << "\t" <<
            currEdge->getRightFace()->getID() << "\t" <<
            currEdge->getLeftFace()->getID() << endl;
    }
    fout << endl;

    m_faces.reset4Loop();
    while (m_faces.setNext4Loop()) {
        FaceBU2D* currFace = m_faces.getpEntity();
        fout << currFace->getID() << "\t" <<
            currFace->getFirstEdge()->getID() << endl; //<< currFace->getBetaSpan().getBetaInterval(1) << endl;
    }

    fout.close();
#endif

    generateSortedFiniteSimplexes();

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_sortEntities = endTime - startTime;

    startTime = clock();
#endif

    computeBetaInterval();

#ifdef _COMP_TIME_BC_IMAGE
    endTime = clock();
    compTime_betaInterval = endTime - startTime;

    ofstream fout_compTime("compTimeForBC_construction.txt");
    fout_compTime << "Entity creation: " << compTime_creatEntities << endl;
    fout_compTime << "Inner face creation with computing tangent circle: " << compTime_createInnerFace_tangentCircle << endl;
    fout_compTime << "Checking slope with 4 disks: " << compTime_checkSlope << endl;
    fout_compTime << "Collecting outside vertices: " << compTime_collect_vtx_outside << endl;
    fout_compTime << "Collecting outside edges: " << compTime_collect_edge_outside << endl;
    fout_compTime << "Bookkeeping: " << compTime_bookkeeping << endl;
    fout_compTime << "Computing betaspan: " << compTime_betaSpan << endl;
    fout_compTime << "Sorting entities with betaspan: " << compTime_sortEntities << endl;
    fout_compTime << "Computing betainterval: " << compTime_betaInterval << endl;

    fout_compTime.close();
#endif

}



rg_FLOAT BetaUniverse2D::calculateRadiusOfPixel(const unsigned char& intensity, const bool& darkerPixelBecomesBiggerDisk)
{
    rg_FLOAT graylevelIntensity = (rg_FLOAT)intensity / 255.0f;
    return darkerPixelBecomesBiggerDisk ? graylevelIntensity / 2.0 : (1.0 - graylevelIntensity) / 2.0;
}

bool BetaUniverse2D::diagonal_is_positive( const Disc* const disk_leftBottom, const Disc* const disk_rightBottom, const Disc* const disk_leftTop, const Disc* const disk_rightTop )
{
    //rg_Circle2D tempCircle[2];
    //rg_Circle2D::makeCircumcircle(disk_leftBottom->getGeometry(), disk_leftTop->getGeometry(), disk_rightBottom->getGeometry(), tempCircle[0], tempCircle[1]);
    
    double baseRadius   = disk_leftBottom->getGeometry().getRadius();
    double r_n          = disk_rightBottom->getGeometry().getRadius() - baseRadius;
    double r_p          = disk_leftTop->getGeometry().getRadius() - baseRadius;

    double r_n_pow_2 = r_n * r_n;
    double r_p_pow_2 = r_p * r_p;
    double r_n_pow_3 = r_n_pow_2 * r_n;
    double r_p_pow_3 = r_p_pow_2 * r_p;
    double r_n_pow_4 = r_n_pow_2 * r_n_pow_2;
    double r_p_pow_4 = r_p_pow_2 * r_p_pow_2;

    double a = 4.0 * ( 1.0 - r_n_pow_2 - r_p_pow_2 );
    double b = 2.0 * ( r_n + r_p - r_n_pow_3 - r_p_pow_3 );
    double c = 2.0 * ( r_n_pow_2 + r_p_pow_2 ) - r_n_pow_4 - r_p_pow_4 - 2.0;

    double enlarged_radius = ( sqrt( b*b - a*c ) - b ) / a;
    double r_tangentCircle = enlarged_radius - baseRadius;

    double x_tangentCircle = - ( r_n_pow_2 + 2.0 * r_n * enlarged_radius - 1.0 ) / 2.0;
    double y_tangentCircle = - ( r_p_pow_2 + 2.0 * r_p * enlarged_radius - 1.0 ) / 2.0;

    double distance_square  = ( 1.0 - x_tangentCircle ) * ( 1.0 - x_tangentCircle ) + ( 1.0 - y_tangentCircle ) * ( 1.0 - y_tangentCircle );
    double sum_radii        = r_tangentCircle + disk_rightTop->getGeometry().getRadius();
    double sum_radii_square = sum_radii * sum_radii;

    return (distance_square < sum_radii_square) ? true : false;
}



void BetaUniverse2D::updateBetaSpan( const int& lengthOfX, const int& lengthOfY )
{
    computeBetaSpanForImageProcessing(lengthOfX, lengthOfY);

    if ( m_sortedFiniteFaces != rg_NULL ) {
        delete [] m_sortedFiniteFaces;
    }
    if ( m_sortedFiniteEdges != rg_NULL ) {
        delete [] m_sortedFiniteEdges;
    }
    if ( m_sortedFiniteVertices != rg_NULL ) {
        delete [] m_sortedFiniteVertices;
    }

    generateSortedFiniteSimplexes();
    computeBetaInterval();
}



rg_FLOAT BetaUniverse2D::calculateRadiusOfTangentCircle_quadraticFormula( const double & baseRadius, const double & prevRadius_CCW, const double & nextRadius_CCW )
{
    double r_n = nextRadius_CCW - baseRadius;
    double r_p = prevRadius_CCW - baseRadius;

    double r_n_pow_2 = r_n * r_n;
    double r_p_pow_2 = r_p * r_p;
    double r_n_pow_3 = r_n_pow_2 * r_n;
    double r_p_pow_3 = r_p_pow_2 * r_p;
    double r_n_pow_4 = r_n_pow_2 * r_n_pow_2;
    double r_p_pow_4 = r_p_pow_2 * r_p_pow_2;

    double a = 4.0 * ( 1.0 - r_n_pow_2 - r_p_pow_2 );
    double b = 2.0 * ( r_n + r_p - r_n_pow_3 - r_p_pow_3 );
    double c = 2.0 * ( r_n_pow_2 + r_p_pow_2 ) - r_n_pow_4 - r_p_pow_4 - 2.0;

    double enlarged_radius = ( sqrt( b*b - a*c ) - b ) / a;

    return enlarged_radius - baseRadius;
}



void BetaUniverse2D::computeBetaSpanForImageProcessing( const int& lengthOfX, const int& lengthOfY )
{
    m_faces.reset4Loop();
    while ( m_faces.setNext4Loop() ) {
        FaceBU2D* currFace = m_faces.getpEntity();
        currFace->computeBetaSpan();
    }

    int numHorizontalEdges  = ( lengthOfX - 1 ) * lengthOfY;
    int numVerticalEdges    = lengthOfX * ( lengthOfY - 1 );
    int numDiagonalEdges    = ( lengthOfX - 1 ) * ( lengthOfY - 1 );

    m_edges.reset4Loop();
    m_edges.setNext4Loop();
    for ( int i = 0; i < numHorizontalEdges; ++i, m_edges.setNext4Loop() ) {
        EdgeBU2D* currHorizontalEdge = m_edges.getpEntity();
        currHorizontalEdge->computeBetaSpan_horizontal_or_vertical_edge();
    }

    for ( int i = 0; i < numVerticalEdges; ++i, m_edges.setNext4Loop() ) {
        EdgeBU2D* currVerticalEdge = m_edges.getpEntity();
        currVerticalEdge->computeBetaSpan_horizontal_or_vertical_edge();
    }

    for ( int i = 0; i < numDiagonalEdges; ++i, m_edges.setNext4Loop() ) {
        EdgeBU2D* currDiagonalEdge = m_edges.getpEntity();
        currDiagonalEdge->computeBetaSpan_diagonal_edge();
    }
}

void BetaUniverse2D::constructForImageProcessing( const vector<vector<int>>& grayDataSet, const bool& darkerPixelBecomesBiggerDisk )
{
    int lengthOfX = grayDataSet.size();
    int lengthOfY = grayDataSet[0].size();
    vector< vector< vector< int > > > dataset;
    dataset.resize( lengthOfX );
    for ( int i = 0; i < lengthOfX; ++i ) {
        dataset[i].resize( lengthOfY );
        for ( int j = 0; j < lengthOfY; ++j ) {
            dataset[i][j].resize( 3 );
            dataset[i][j][0] = grayDataSet[i][j];
            dataset[i][j][1] = grayDataSet[i][j];
            dataset[i][j][2] = grayDataSet[i][j];
        }
    }

    constructForImageProcessing( dataset, lengthOfX, lengthOfY, darkerPixelBecomesBiggerDisk );
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