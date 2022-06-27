#include "EdgeWEDS.h"
#include "VertexWEDS.h"
#include "FaceWEDS.h"



EdgeWEDS::EdgeWEDS()
: m_startVertex(rg_NULL), m_endVertex(rg_NULL), 
  m_leftFace(rg_NULL),    m_rightFace(rg_NULL), 
  m_leftHand(rg_NULL),    m_rightHand(rg_NULL), 
  m_leftLeg(rg_NULL),     m_rightLeg(rg_NULL), 
  m_visited(rg_FALSE)
{
}



EdgeWEDS::EdgeWEDS(const rg_INT& ID)
: TopologicalEntity(ID), 
  m_startVertex(rg_NULL), m_endVertex(rg_NULL), 
  m_leftFace(rg_NULL),    m_rightFace(rg_NULL), 
  m_leftHand(rg_NULL),    m_rightHand(rg_NULL), 
  m_leftLeg(rg_NULL),     m_rightLeg(rg_NULL), 
  m_visited(rg_FALSE)
{
}



EdgeWEDS::EdgeWEDS(const rg_INT& ID, VertexWEDS* startVertex, VertexWEDS* endVertex)
: TopologicalEntity(ID), 
  m_startVertex(startVertex), m_endVertex(endVertex), 
  m_leftFace(rg_NULL),        m_rightFace(rg_NULL), 
  m_leftHand(rg_NULL),        m_rightHand(rg_NULL), 
  m_leftLeg(rg_NULL),         m_rightLeg(rg_NULL), 
  m_visited(rg_FALSE)
{
}



EdgeWEDS::EdgeWEDS(const rg_INT& ID, VertexWEDS* startVertex, VertexWEDS* endVertex,
                                     FaceWEDS*   leftFace,    FaceWEDS*   rightFace,                              
                                     EdgeWEDS*   leftHand,    EdgeWEDS*   rightHand,               
                                     EdgeWEDS*   leftLeg,     EdgeWEDS*   rightLeg)
: TopologicalEntity(ID), 
  m_startVertex(startVertex), m_endVertex(endVertex), 
  m_leftFace(leftFace),       m_rightFace(rightFace), 
  m_leftHand(leftHand),       m_rightHand(rightHand), 
  m_leftLeg(leftLeg),         m_rightLeg(rightLeg), 
  m_visited(rg_FALSE)
{
}



EdgeWEDS::EdgeWEDS(const EdgeWEDS& edge)
: TopologicalEntity(edge), 
  m_startVertex(edge.m_startVertex), m_endVertex(edge.m_endVertex), 
  m_leftFace(edge.m_leftFace),       m_rightFace(edge.m_rightFace), 
  m_leftHand(edge.m_leftHand),       m_rightHand(edge.m_rightHand), 
  m_leftLeg(edge.m_leftLeg),         m_rightLeg(edge.m_rightLeg), 
  m_visited(edge.m_visited)
{
}



EdgeWEDS::~EdgeWEDS()
{
}



void EdgeWEDS::setEdge( VertexWEDS* startVertex, VertexWEDS* endVertex,
                        FaceWEDS*   leftFace,    FaceWEDS*   rightFace,                              
                        EdgeWEDS*   leftHand,    EdgeWEDS*   rightHand,               
                        EdgeWEDS*   leftLeg,     EdgeWEDS*   rightLeg)
{
    m_startVertex = startVertex;
    m_endVertex   = endVertex;
    m_leftFace    = leftFace;
    m_rightFace   = rightFace;
    m_leftHand    = leftHand;
    m_rightHand   = rightHand;
    m_leftLeg     = leftLeg;
    m_rightLeg    = rightLeg;
}



EdgeWEDS& EdgeWEDS::operator=(const EdgeWEDS& edge)
{
    if ( this == &edge ) {
        return *this;
    }

    TopologicalEntity::operator =(edge);
    m_startVertex = edge.m_startVertex;
    m_endVertex   = edge.m_endVertex;
    m_leftFace    = edge.m_leftFace;
    m_rightFace   = edge.m_rightFace;
    m_leftHand    = edge.m_leftHand;
    m_rightHand   = edge.m_rightHand;
    m_leftLeg     = edge.m_leftLeg;
    m_rightLeg    = edge.m_rightLeg;
    m_visited     = edge.m_visited;

    return *this;
}




//  Topological operators
void EdgeWEDS::connectRightHand( EdgeWEDS* rightHand)
{
    m_rightHand = rightHand;

    if ( m_endVertex == rightHand->m_startVertex ) {
        rightHand->m_rightLeg = this;
    }
    else {
        rightHand->m_leftHand = this;
    }
}



void EdgeWEDS::connectLeftHand(  EdgeWEDS* leftHand)
{
    m_leftHand = leftHand;

    if ( m_endVertex == leftHand->m_startVertex ) {
        leftHand->m_leftLeg   = this;
    }
    else {
        leftHand->m_rightHand = this;
    }
}



void EdgeWEDS::connectRightLeg(  EdgeWEDS* rightLeg)
{
    m_rightLeg = rightLeg;

    if ( m_startVertex == rightLeg->m_startVertex ) {
        rightLeg->m_leftLeg = this;
    }
    else {
        rightLeg->m_rightHand = this;
    }
}



void EdgeWEDS::connectLeftLeg(   EdgeWEDS* leftLeg)
{
    m_leftLeg = leftLeg;

    if ( m_startVertex == leftLeg->m_startVertex ) {
        leftLeg->m_rightLeg = this;
    }
    else {
        leftLeg->m_leftHand = this;
    }
}




void EdgeWEDS::connectRightFace( FaceWEDS* rightFace)
{
    m_rightFace = rightFace;

    if ( rightFace != rg_NULL ) {
        rightFace->setFirstEdge( this );
    }
}



void EdgeWEDS::connectLeftFace(  FaceWEDS* leftFace)
{
    m_leftFace = leftFace;

    if ( leftFace != rg_NULL ) {
        leftFace->setFirstEdge( this );
    }
}



void EdgeWEDS::connectStartVertex( VertexWEDS* vertex )
{
    m_startVertex = vertex;

    if ( vertex != rg_NULL ) {
        vertex->setFirstEdge( this );
    }
}



void EdgeWEDS::connectEndVertex( VertexWEDS* vertex )
{
    m_endVertex = vertex;

    if ( vertex != rg_NULL ) {
        vertex->setFirstEdge( this );
    }
}



void EdgeWEDS::reverse()
{
    VertexWEDS* startVertex = m_startVertex;
    VertexWEDS* endVertex   = m_endVertex;

    EdgeWEDS*   leftHand    = m_leftHand;
    EdgeWEDS*   rightHand   = m_rightHand;
    EdgeWEDS*   leftLeg     = m_leftLeg;
    EdgeWEDS*   rightLeg    = m_rightLeg;

    FaceWEDS*   leftFace    = m_leftFace;
    FaceWEDS*   rightFace   = m_rightFace;


    m_startVertex = endVertex;
    m_endVertex   = startVertex;

    m_leftHand    = rightLeg;
    m_rightHand   = leftLeg;
    m_leftLeg     = rightHand;
    m_rightLeg    = leftHand;

    m_leftFace    = rightFace;
    m_rightFace   = leftFace;
}


    
void EdgeWEDS::flip()
{
	EdgeWEDS*   RH = m_rightHand;
	EdgeWEDS*   LH = m_leftHand;
	EdgeWEDS*   RL = m_rightLeg;
	EdgeWEDS*   LL = m_leftLeg;
	
	VertexWEDS* SV = m_startVertex;
	VertexWEDS* EV = m_endVertex;

	FaceWEDS*   RF = m_rightFace;
	FaceWEDS*   LF = m_leftFace;
	
	//move first edge on two faces
	if( this == m_leftFace->getFirstEdge() )	{
		if( LH != NULL) {
			m_leftFace->setFirstEdge( LH );
		}
		else if( LL != NULL) {
			m_leftFace->setFirstEdge( LL );
		}
		else {
			return;
		}
	}

	if( this == m_rightFace->getFirstEdge() ) {
		if(RH != NULL) {
			m_rightFace->setFirstEdge( RH );
		}
		else if(RL != NULL) {
			m_rightFace->setFirstEdge( RL );
		}
		else {
			return; //error
		}
	}
	
	//move first edge on two vertices (start/end)
	SV->setFirstEdge( this );
	EV->setFirstEdge( this );
	
	//change vertex(start or end) of incident edges
	//left hand and right leg does not change start or end vertex
	if( RH->getStartVertex() == EV) {
        RH->m_startVertex = SV;
        RH->connectRightLeg( RL );
	}
    else { //end	
        RH->m_endVertex = SV;
		RH->connectLeftHand( RL );
	}
	
	if(LL->getStartVertex() == SV) {
        LL->m_startVertex = EV;
		LL->connectRightLeg( LH );
	}
	else {
        LL->m_endVertex = EV;
		LL->connectLeftHand( LH );
	}

	//change right and left face
	if( LF == LH->getLeftFace() ) {
		m_rightFace = LH->getRightFace();
	}
	else {
		m_rightFace = LH->getLeftFace();
	}
	
	if( RF == RL->getRightFace() ) {
		m_leftFace = RL->getLeftFace();
	}
	else {
		m_leftFace = RL->getRightFace();
	}
			
	//update topology of two legs and two hands
    connectRightHand( LH );
    connectLeftHand(  LL );
    connectRightLeg(  RH );
    connectLeftLeg(   RL );
}



rg_BOOL EdgeWEDS::isConnected(EdgeWEDS* otherEdge) const
{
    if ( otherEdge->isIncidentTo( m_startVertex ) || otherEdge->isIncidentTo( m_endVertex ) ){
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}




VertexWEDS* EdgeWEDS::getVertexOfLeftHand() const
{
    if ( m_endVertex == m_leftHand->m_startVertex ) {
        return m_leftHand->m_endVertex;
    }
    else {
        return m_leftHand->m_startVertex;
    }
}



VertexWEDS* EdgeWEDS::getVertexOfRightHand() const
{
    if ( m_endVertex == m_rightHand->m_startVertex ) {
        return m_rightHand->m_endVertex;
    }
    else {
        return m_rightHand->m_startVertex;
    }
}


    
rg_INT EdgeWEDS::getNeighborVertices(rg_dList<VertexWEDS*>& vertexList) const
{
    rg_dList<VertexWEDS*> neighborOfVtxs;
    m_startVertex->getNeighborVertices( neighborOfVtxs );
    m_endVertex->getNeighborVertices(   neighborOfVtxs );

    neighborOfVtxs.reset4Loop();
    while ( neighborOfVtxs.setNext4Loop() ) {
        VertexWEDS* currVtx = neighborOfVtxs.getEntity();
        if ( currVtx != m_startVertex && currVtx != m_endVertex ) {
            vertexList.add( currVtx );
        }
    }

    return vertexList.getSize();
}



rg_INT EdgeWEDS::getIncidentEdges(rg_dList<EdgeWEDS*>& edgeList) const
{
    rg_dList<EdgeWEDS*> edgesIncidentToVtxs;
    m_startVertex->getIncidentEdges( edgesIncidentToVtxs );
    m_endVertex->getIncidentEdges(   edgesIncidentToVtxs );

    edgesIncidentToVtxs.reset4Loop();
    while ( edgesIncidentToVtxs.setNext4Loop() ) {
        EdgeWEDS* currEdge = edgesIncidentToVtxs.getEntity();
        if ( currEdge != this ) {
            edgeList.add( currEdge );
        }
    }

    return edgeList.getSize();
}



rg_INT EdgeWEDS::getIncidentFaces(rg_dList<FaceWEDS*>& faceList) const
{
    m_startVertex->getIncidentFaces( faceList );

    rg_dList<FaceWEDS*> facesIncidentToEndVertex;
    m_endVertex->getIncidentFaces( facesIncidentToEndVertex );

    facesIncidentToEndVertex.reset4Loop();
    while ( facesIncidentToEndVertex.setNext4Loop() ) {
        FaceWEDS* currFace = facesIncidentToEndVertex.getEntity();

        if ( currFace != m_leftFace && currFace != m_rightFace ) {
            faceList.add( currFace );
        }
    }

    return faceList.getSize();
}


rg_INT EdgeWEDS::getEdgesInStar(rg_dList<EdgeWEDS*>& edgeList) const
{
    rg_dList<EdgeWEDS*> edgesInStarOfVtxs;
    m_startVertex->getEdgesInStar( edgesInStarOfVtxs );
    m_endVertex->getEdgesInStar( edgesInStarOfVtxs );

    edgesInStarOfVtxs.reset4Loop();
    while ( edgesInStarOfVtxs.setNext4Loop() ) {
        EdgeWEDS* currEdge = edgesInStarOfVtxs.getEntity();

        if (    currEdge != m_leftHand 
             && currEdge != m_rightHand 
             && currEdge != m_leftLeg 
             && currEdge != m_rightLeg
             && currEdge != this) {
            edgeList.add( currEdge );
        }
    }

    edgeList.add( m_leftHand );
    edgeList.add( m_rightHand );
    edgeList.add( m_leftLeg );
    edgeList.add( m_rightLeg );

    return edgeList.getSize();
}



rg_INT EdgeWEDS::getFacesInStar(rg_dList<FaceWEDS*>& faceList) const
{
    return getIncidentFaces(faceList);
}



