#include "WingedEdgeDataStructure.h"


#include <set>
using namespace std;


WingedEdgeDataStructure::WingedEdgeDataStructure()
: m_number_of_faces(0),
  m_number_of_edges(0),
  m_number_of_vertices(0)
{
}



WingedEdgeDataStructure::WingedEdgeDataStructure(const WingedEdgeDataStructure& weds)
: m_number_of_faces(0),
  m_number_of_edges(0),
  m_number_of_vertices(0)
{
    duplicate(weds);
}



WingedEdgeDataStructure::~WingedEdgeDataStructure()
{
    clear();
}



void    WingedEdgeDataStructure::clear()
{
    for (list<Face*>::iterator i_face=m_faces.begin(); i_face!=m_faces.end(); ++i_face ) {
        delete (*i_face);
    }

    for (list<Edge*>::iterator i_edge=m_edges.begin(); i_edge!=m_edges.end(); ++i_edge ) {
        delete (*i_edge);
    }

    for (list<Vertex*>::iterator i_vx=m_vertices.begin(); i_vx!=m_vertices.end(); ++i_vx ) {
        delete (*i_vx);
    }


    m_faces.clear();
    m_edges.clear();
    m_vertices.clear();

    m_number_of_faces = 0;
    m_number_of_edges = 0;
    m_number_of_vertices = 0;
}

    

WingedEdgeDataStructure& WingedEdgeDataStructure::operator =(const WingedEdgeDataStructure& weds)
{
    if ( this != &weds ) {
        duplicate(weds);
    }

    return *this;
}



WingedEdgeDataStructure::Face*           WingedEdgeDataStructure::create_face(        const WingedEdgeDataStructure::Face&        face )
{
    m_faces.push_back( new Face(face) );
    ++m_number_of_faces;

    return m_faces.back();
}



WingedEdgeDataStructure::Edge*           WingedEdgeDataStructure::create_edge(        const WingedEdgeDataStructure::Edge&        edge )
{
    m_edges.push_back( new Edge(edge) );
    ++m_number_of_edges;

    return m_edges.back();
}



WingedEdgeDataStructure::Vertex*         WingedEdgeDataStructure::create_vertex(      const WingedEdgeDataStructure::Vertex&      vertex )
{
    m_vertices.push_back( new Vertex(vertex) );
    ++m_number_of_vertices;

    return m_vertices.back();
}



WingedEdgeDataStructure::Face*   WingedEdgeDataStructure::create_face()
{
    int ID = 0;
    if (!m_faces.empty()) {
        ID = m_faces.back()->getID() + 1;
    }
    m_faces.push_back( new Face(ID) );
    ++m_number_of_faces;

    return m_faces.back();
}



WingedEdgeDataStructure::Edge*   WingedEdgeDataStructure::create_edge()
{
    int ID = 0;
    if (!m_edges.empty()) {
        ID = m_edges.back()->getID() + 1;
    }
    m_edges.push_back( new Edge(ID) );
    ++m_number_of_edges;

    return m_edges.back();
}



WingedEdgeDataStructure::Vertex* WingedEdgeDataStructure::create_vertex()
{
    int ID = 0;
    if (!m_vertices.empty()) {
        ID = m_vertices.back()->getID() + 1;
    }

    m_vertices.push_back( new Vertex(ID) );
    ++m_number_of_vertices;

    return m_vertices.back();
}



void            WingedEdgeDataStructure::concatenate_face(   const WingedEdgeDataStructure::Face*        const face )
{
    m_faces.push_back( const_cast<Face*>(face) ); 
    ++m_number_of_faces;
}



void            WingedEdgeDataStructure::concatenate_edge(   const WingedEdgeDataStructure::Edge*        const edge )
{
    m_edges.push_back( const_cast<Edge*>(edge) ); 
    ++m_number_of_edges;
}



void            WingedEdgeDataStructure::concatenate_vertex( const WingedEdgeDataStructure::Vertex*      const vertex )
{
    m_vertices.push_back( const_cast<Vertex*>(vertex) ); 
    ++m_number_of_vertices;
}



void            WingedEdgeDataStructure::remove_face(        const WingedEdgeDataStructure::Face*        const face )
{
    if ( face != rg_NULL ) {
        m_faces.remove( const_cast<Face*>(face) );
        --m_number_of_faces;
        delete face;
    }
}



void            WingedEdgeDataStructure::remove_edge(        const WingedEdgeDataStructure::Edge*        const edge )
{
    if ( edge != rg_NULL ) {
        m_edges.remove( const_cast<Edge*>(edge) );
        --m_number_of_edges;
        delete edge;
    }
}



void            WingedEdgeDataStructure::remove_vertex(      const WingedEdgeDataStructure::Vertex*      const vertex )
{
    if ( vertex != rg_NULL ) {
        m_vertices.remove( const_cast<Vertex*>(vertex) );
        --m_number_of_vertices;
        delete vertex;
    }
}



void    WingedEdgeDataStructure::duplicate(const WingedEdgeDataStructure& weds)
{
    clear();

    map<Face*,        Face*>        faceMap;
    map<Edge*,        Edge*>        edgeMap;
    map<Vertex*,      Vertex*>      vertexMap;
    make_entity_map_for_duplication(    weds, faceMap, edgeMap, vertexMap );

    duplicate_topology(                 weds, faceMap, edgeMap, vertexMap );
}



void            WingedEdgeDataStructure::make_entity_map_for_duplication( 
                                                 const WingedEdgeDataStructure&     weds, 
                                                 map<Face*,        Face*>&              faceMap,
                                                 map<Edge*,        Edge*>&              edgeMap,
                                                 map<Vertex*,      Vertex*>&            vertexMap)
{
    //  make map for V-faces
    list<Face*>::const_iterator i_face;
    for ( i_face=weds.m_faces.begin(); i_face!=weds.m_faces.end(); ++i_face ) {
        Face* originFace = const_cast<Face*>( (*i_face) );
        Face* currFace   = this->create_face();
        currFace->setID( originFace->getID() );

        faceMap.insert( make_pair( originFace, currFace ) );
    }
    faceMap.insert( make_pair( (Face*)rg_NULL, (Face*)rg_NULL ) );


    //  make map for V-edges
    list<Edge*>::const_iterator i_edge;
    for ( i_edge=weds.m_edges.begin(); i_edge!=weds.m_edges.end(); ++i_edge ) {
        Edge* originEdge = const_cast<Edge*>( (*i_edge) );
        Edge* currEdge   = this->create_edge();
        currEdge->setID( originEdge->getID() );

        edgeMap.insert( make_pair( originEdge, currEdge ) );
    }
    edgeMap.insert( make_pair( (Edge*)rg_NULL, (Edge*)rg_NULL ) );


    //  make map for V-vertices
    list<Vertex*>::const_iterator i_vtx;
    for ( i_vtx=weds.m_vertices.begin(); i_vtx!=weds.m_vertices.end(); ++i_vtx ) {
        Vertex* originVtx = const_cast<Vertex*>( (*i_vtx) );
        Vertex* currVtx   = this->create_vertex();
        currVtx->setID( originVtx->getID() );

        vertexMap.insert( make_pair( originVtx, currVtx ) );
    }
    vertexMap.insert( make_pair( (Vertex*)rg_NULL, (Vertex*)rg_NULL ) );
}



void            WingedEdgeDataStructure::duplicate_topology(     
                                           const WingedEdgeDataStructure&     weds, 
                                           const map<Face*,        Face*>&        faceMap,
                                           const map<Edge*,        Edge*>&        edgeMap,
                                           const map<Vertex*,      Vertex*>&      vertexMap)
{
    list<Face*>::const_iterator i_face;
    for ( i_face=weds.m_faces.begin(); i_face!=weds.m_faces.end(); i_face++ ) {
        Face* originFace = const_cast<Face*>( (*i_face) );
        Face* currFace   = faceMap.find( originFace )->second;

        currFace->set_first_edge( edgeMap.find( originFace->first_edge() )->second );
    }


    list<Edge*>::const_iterator i_edge;
    for ( i_edge=weds.m_edges.begin(); i_edge!=weds.m_edges.end(); i_edge++ ) {
        Edge* originEdge = const_cast<Edge*>( (*i_edge) );
        Edge* currEdge   = edgeMap.find( originEdge )->second;

        currEdge->set_edge( vertexMap.find( originEdge->start_vertex()    )->second,
                            vertexMap.find( originEdge->end_vertex()      )->second,
                            faceMap.find(   originEdge->right_face()      )->second,
                            faceMap.find(   originEdge->left_face()       )->second,
                            edgeMap.find(   originEdge->right_hand_edge() )->second,
                            edgeMap.find(   originEdge->left_hand_edge()  )->second,
                            edgeMap.find(   originEdge->right_leg_edge()  )->second,
                            edgeMap.find(   originEdge->left_leg_edge()   )->second );

    }


    list<Vertex*>::const_iterator i_vtx;
    for ( i_vtx=weds.m_vertices.begin(); i_vtx!=weds.m_vertices.end(); i_vtx++ ) {
        Vertex* originVtx = const_cast<Vertex*>( (*i_vtx) );
        Vertex* currVtx   = vertexMap.find( originVtx )->second;

        currVtx->set_first_edge( edgeMap.find( originVtx->first_edge() )->second );
    }
}





///////////////////////////////////////////////////////////////////////////////
//
// class WingedEdgeDataStructure::Face : public TopologicalEntity
//
WingedEdgeDataStructure::Face::Face() 
: TopologicalEntity(0),
  m_firstEdge(rg_NULL)
{
}



WingedEdgeDataStructure::Face::Face(const rg_INT& ID)
: TopologicalEntity(ID),
  m_firstEdge(rg_NULL)
{
}



WingedEdgeDataStructure::Face::Face(const WingedEdgeDataStructure::Face& face)
: TopologicalEntity( face ),
  m_firstEdge(       face.m_firstEdge )
{
}



WingedEdgeDataStructure::Face::~Face()
{
}



WingedEdgeDataStructure::Face&   WingedEdgeDataStructure::Face::operator =(const WingedEdgeDataStructure::Face& face)
{
    if ( this != &face ) {
        TopologicalEntity::operator=( face );
        m_firstEdge = face.m_firstEdge;
    }

    return *this;
}



bool    WingedEdgeDataStructure::Face::is_incident_to(const WingedEdgeDataStructure::Vertex* const vertex) const
{
    WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
    WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    do  {
        if ( this == currEdge->left_face() )  {
            if ( vertex == currEdge->end_vertex() ) {
                return true;
            }
            else {
                currEdge = currEdge->left_hand_edge();
            }
        }
        else if ( this == currEdge->right_face() )  {
            if ( vertex == currEdge->start_vertex() ) {
                return true;
            }
            else {            
                currEdge = currEdge->right_leg_edge();
            }
        }
        else {
            break;
        }

    } while ( startEdge != currEdge );   

    return false;
}



bool    WingedEdgeDataStructure::Face::is_incident_to(const WingedEdgeDataStructure::Edge* const edge) const
{
    return edge->is_member_of(this);
}



bool    WingedEdgeDataStructure::Face::is_adjacent_to(const WingedEdgeDataStructure::Face* const face) const
{
    WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
    WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    do  {
        if ( this == currEdge->left_face() )  {
            if ( face == currEdge->right_face() ) {
                return true;
            }
            else {
                currEdge = currEdge->left_hand_edge();
            }
        }
        else if ( this == currEdge->right_face() )  {
            if ( face == currEdge->left_face() ) {
                return true;
            }
            else {            
                currEdge = currEdge->right_leg_edge();
            }
        }
        else {
            break;
        }

    } while ( startEdge != currEdge );   

    return false;
}



WingedEdgeDataStructure::Edge*   WingedEdgeDataStructure::Face::find_edge(
                                            const WingedEdgeDataStructure::Vertex* const vertex1, 
                                            const WingedEdgeDataStructure::Vertex* const vertex2 ) const
{
    /*
    WingedEdgeDataStructure::Edge* edgeWanted = rg_NULL;

    list<WingedEdgeDataStructure::Edge*> edgesIncidentToV1;
    vertex1->find_incident_edges( edgesIncidentToV1 );

    for ( list<Edge*>::iterator i_edge=edgesIncidentToV1.begin(); i_edge!=edgesIncidentToV1.end(); ++i_edge ) {
        Edge* currEdge = *i_edge;
        if ( this->is_incident_to( currEdge ) && currEdge->is_incident_to( vertex2 ) ) {
            edgeWanted = currEdge;
            break;
        }
    }

    return edgeWanted;
    */


    
    list<WingedEdgeDataStructure::Edge*> boundingEdges;
    find_bounding_edges( boundingEdges);

    for ( list<Edge*>::iterator i_edge=boundingEdges.begin(); i_edge!=boundingEdges.end(); ++i_edge ) {
        Edge* currEdge = *i_edge;
        if ( currEdge->is_incident_to( vertex1 ) && currEdge->is_incident_to( vertex2 ) ) {
            return currEdge;
        }
    }

    return rg_NULL;
}



rg_INT  WingedEdgeDataStructure::Face::number_of_bounding_vertices() const
{
    list<WingedEdgeDataStructure::Vertex*> boundingVertices;
    return find_bounding_vertices( boundingVertices );
}



rg_INT  WingedEdgeDataStructure::Face::number_of_bounding_edges() const
{
    list<WingedEdgeDataStructure::Edge*> boundingEdges;
    return find_bounding_edges( boundingEdges );
}



rg_INT  WingedEdgeDataStructure::Face::number_of_adjacent_faces() const
{
    list<WingedEdgeDataStructure::Face*> adjacentFaces;
    return find_adjacent_faces( adjacentFaces );
}




rg_INT  WingedEdgeDataStructure::Face::find_bounding_vertices( list<WingedEdgeDataStructure::Vertex*>& boundingVertices) const
{
    WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
    WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    unsigned int numBoundingVertices = 0;

    do  {
        if ( this == currEdge->left_face() )  {
            boundingVertices.push_back( currEdge->end_vertex() );
            ++numBoundingVertices;

            currEdge = currEdge->left_hand_edge();
        }
        else if ( this == currEdge->right_face() )  {
            boundingVertices.push_back( currEdge->start_vertex() );
            ++numBoundingVertices;

            currEdge = currEdge->right_leg_edge();
        }
        else  {
            break;
        }
    } while ( currEdge != rg_NULL && currEdge != startEdge );   

    return numBoundingVertices;
}



rg_INT  WingedEdgeDataStructure::Face::find_bounding_edges(    list<WingedEdgeDataStructure::Edge*>&   boundingEdges) const
{
    WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
    WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    unsigned int numBoundingEdges = 0;
    do  {
        boundingEdges.push_back( currEdge );
        ++numBoundingEdges;

        if ( this == currEdge->left_face() ) {
            currEdge = currEdge->left_hand_edge();
        }
        else if ( this == currEdge->right_face() )  {
            currEdge = currEdge->right_leg_edge();
        }
        else {
            break;
        }
    } while ( currEdge != rg_NULL && currEdge != startEdge );   

    return numBoundingEdges;
}



rg_INT  WingedEdgeDataStructure::Face::find_adjacent_faces(    list<WingedEdgeDataStructure::Face*>&   adjacentFaces) const
{
    WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
    WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    unsigned int numAdjacentFaces = 0;
    do  {
        WingedEdgeDataStructure::Face* rightFace = currEdge->left_face();
        WingedEdgeDataStructure::Face* leftFace  = currEdge->right_face();

        if ( this == leftFace )  {
            if ( rightFace != rg_NULL ) {
                adjacentFaces.push_back( rightFace );
                ++numAdjacentFaces;
            }

            currEdge = currEdge->left_hand_edge();
        }
        else if ( this == rightFace )  {
            if ( leftFace != rg_NULL ) {
                adjacentFaces.push_back( leftFace );
                ++numAdjacentFaces;
            }

            currEdge = currEdge->right_leg_edge();
        }
        else  {
            break;
        }
    } while ( currEdge != rg_NULL && currEdge != startEdge );   

    return numAdjacentFaces;
}





///////////////////////////////////////////////////////////////////////////////
//
// class WingedEdgeDataStructure::Edge : public TopologicalEntity
//
WingedEdgeDataStructure::Edge::Edge()
: TopologicalEntity(0),
  m_startVertex(rg_NULL ),
  m_endVertex(  rg_NULL ), 
  m_rightFace(  rg_NULL ),
  m_leftFace(   rg_NULL ),
  m_rightHand(  rg_NULL ),
  m_leftHand(   rg_NULL ),
  m_rightLeg(   rg_NULL ),
  m_leftLeg(    rg_NULL )
{
}



WingedEdgeDataStructure::Edge::Edge(const rg_INT& ID)
: TopologicalEntity(ID),
  m_startVertex(rg_NULL ),
  m_endVertex(  rg_NULL ), 
  m_rightFace(  rg_NULL ),
  m_leftFace(   rg_NULL ),
  m_rightHand(  rg_NULL ),
  m_leftHand(   rg_NULL ),
  m_rightLeg(   rg_NULL ),
  m_leftLeg(    rg_NULL )
{
}



WingedEdgeDataStructure::Edge::Edge(const WingedEdgeDataStructure::Edge& edge)
: TopologicalEntity( edge ),
  m_startVertex(     edge.m_startVertex ),
  m_endVertex(       edge.m_endVertex ), 
  m_rightFace(       edge.m_rightFace ),
  m_leftFace(        edge.m_leftFace ),
  m_rightHand(       edge.m_rightHand ),
  m_leftHand(        edge.m_leftHand ),
  m_rightLeg(        edge.m_rightLeg ),
  m_leftLeg(         edge.m_leftLeg )
{
}



WingedEdgeDataStructure::Edge::~Edge()
{
}



WingedEdgeDataStructure::Edge& WingedEdgeDataStructure::Edge::operator =(const WingedEdgeDataStructure::Edge& edge)
{
    if ( this != &edge ) {
        TopologicalEntity::operator=( edge );

        m_startVertex = edge.m_startVertex;
        m_endVertex   = edge.m_endVertex; 
        m_rightFace   = edge.m_rightFace;
        m_leftFace    = edge.m_leftFace;
        m_rightHand   = edge.m_rightHand;
        m_leftHand    = edge.m_leftHand;
        m_rightLeg    = edge.m_rightLeg;
        m_leftLeg     = edge.m_leftLeg;
    }

    return *this;
}



WingedEdgeDataStructure::Vertex* WingedEdgeDataStructure::Edge::opposite_vertex(const WingedEdgeDataStructure::Vertex* const vertex) const
{
    if ( vertex == m_startVertex ) {
        return m_endVertex;
    }
    else if ( vertex == m_endVertex ) {
        return m_startVertex;
    }
    else {
        return rg_NULL;
    }
}

    

WingedEdgeDataStructure::Face*   WingedEdgeDataStructure::Edge::opposite_face(const WingedEdgeDataStructure::Face* const face) const
{
    if ( face == m_rightFace ) {
        return m_leftFace;
    }
    else if ( face == m_leftFace ) {
        return m_rightFace;
    }
    else {
        return rg_NULL;
    }
}



void    WingedEdgeDataStructure::Edge::set_edge(  
                    const WingedEdgeDataStructure::Vertex* const startVertex ,  const WingedEdgeDataStructure::Vertex* const endVertex,
                    const WingedEdgeDataStructure::Face*   const rightFace,     const WingedEdgeDataStructure::Face*   const leftFace,
                    const WingedEdgeDataStructure::Edge*   const rightHandEdge, const WingedEdgeDataStructure::Edge*   const leftHandEdge,
                    const WingedEdgeDataStructure::Edge*   const rightLegEdge,  const WingedEdgeDataStructure::Edge*   const leftLegEdge)
{
    m_startVertex   = const_cast<WingedEdgeDataStructure::Vertex*>(startVertex);
    m_endVertex     = const_cast<WingedEdgeDataStructure::Vertex*>(endVertex);

    m_rightFace     = const_cast<WingedEdgeDataStructure::Face*>(rightFace);
    m_leftFace      = const_cast<WingedEdgeDataStructure::Face*>(leftFace);

    m_rightHand     = const_cast<WingedEdgeDataStructure::Edge*>(rightHandEdge);
    m_leftHand      = const_cast<WingedEdgeDataStructure::Edge*>(leftHandEdge);
    m_rightLeg      = const_cast<WingedEdgeDataStructure::Edge*>(rightLegEdge);
    m_leftLeg       = const_cast<WingedEdgeDataStructure::Edge*>(leftLegEdge);
}



rg_INT  WingedEdgeDataStructure::Edge::find_bounding_vertices( list<WingedEdgeDataStructure::Vertex*>& boundingVertices) const
{
    if ( m_startVertex != rg_NULL )
        boundingVertices.push_back( m_startVertex );
    if (m_endVertex != rg_NULL)
        boundingVertices.push_back( m_endVertex );

    return 2;
}



rg_INT  WingedEdgeDataStructure::Edge::find_adjacent_edges(    list<WingedEdgeDataStructure::Edge*>&   adjacentEdges) const
{
    list<WingedEdgeDataStructure::Edge*> edgesIncidentToTwoVertices;
    m_startVertex->find_incident_edges( edgesIncidentToTwoVertices );
    m_endVertex->find_incident_edges(   edgesIncidentToTwoVertices );

    unsigned int numAdjacentEdges = 0;
    list<WingedEdgeDataStructure::Edge*>::iterator i_edge;
    for (i_edge=edgesIncidentToTwoVertices.begin(); i_edge!=edgesIncidentToTwoVertices.end(); ++i_edge ) {
        if ( *i_edge != this ) {
            adjacentEdges.push_back( *i_edge );
            ++numAdjacentEdges;
        }
    }

    return numAdjacentEdges;
}



rg_INT  WingedEdgeDataStructure::Edge::find_incident_faces(    list<WingedEdgeDataStructure::Face*>&   incidentFaces) const
{
    if ( m_rightFace != rg_NULL ) {
        incidentFaces.push_back( m_rightFace );
    }

    if ( m_leftFace != rg_NULL ) {
        incidentFaces.push_back( m_leftFace );
    }

    return 2;
}




void    WingedEdgeDataStructure::Edge::reverse()
{
    WingedEdgeDataStructure::Vertex* startVertex = m_startVertex;
    WingedEdgeDataStructure::Vertex* endVertex   = m_endVertex;

    WingedEdgeDataStructure::Edge*   leftHand    = m_leftHand;
    WingedEdgeDataStructure::Edge*   rightHand   = m_rightHand;
    WingedEdgeDataStructure::Edge*   leftLeg     = m_leftLeg;
    WingedEdgeDataStructure::Edge*   rightLeg    = m_rightLeg;

    WingedEdgeDataStructure::Face*   leftFace    = m_leftFace;
    WingedEdgeDataStructure::Face*   rightFace   = m_rightFace;


    m_startVertex = endVertex;
    m_endVertex   = startVertex;

    m_leftHand    = rightLeg;
    m_rightHand   = leftLeg;
    m_leftLeg     = rightHand;
    m_rightLeg    = leftHand;

    m_leftFace    = rightFace;
    m_rightFace   = leftFace;
}

    

void    WingedEdgeDataStructure::Edge::flip()
{
	WingedEdgeDataStructure::Edge*   RH = m_rightHand;
	WingedEdgeDataStructure::Edge*   LH = m_leftHand;
	WingedEdgeDataStructure::Edge*   RL = m_rightLeg;
	WingedEdgeDataStructure::Edge*   LL = m_leftLeg;
	
	WingedEdgeDataStructure::Vertex* SV = m_startVertex;
	WingedEdgeDataStructure::Vertex* EV = m_endVertex;

	WingedEdgeDataStructure::Face*   RF = m_rightFace;
	WingedEdgeDataStructure::Face*   LF = m_leftFace;
	
	//move first edge on two faces
    if( this == m_leftFace->first_edge() )	{
		if( LH != rg_NULL) {
			m_leftFace->set_first_edge( LH );
		}
		else if( LL != rg_NULL) {
			m_leftFace->set_first_edge( LL );
		}
		else {
			return;
		}
	}

	if( this == m_rightFace->first_edge() ) {
		if(RH != rg_NULL) {
			m_rightFace->set_first_edge( RH );
		}
		else if(RL != rg_NULL) {
			m_rightFace->set_first_edge( RL );
		}
		else {
			return; //error
		}
	}
	
	//move first edge on two vertices (start/end)
	SV->set_first_edge( this );
	EV->set_first_edge( this );
	
	//change vertex(start or end) of incident edges
	//left hand and right leg does not change start or end vertex
    if( RH->m_startVertex == EV) {
        RH->m_startVertex = SV;
        RH->connect_right_leg_edge( RL );
	}
    else { //end	
        RH->m_endVertex = SV;
		RH->connect_left_hand_edge( RL );
	}
	
	if(LL->m_startVertex == SV) {
        LL->m_startVertex = EV;
		LL->connect_right_leg_edge( LH );
	}
	else {
        LL->m_endVertex = EV;
		LL->connect_left_hand_edge( LH );
	}

	//change right and left face
	if( LF == LH->m_leftFace ) {
		m_rightFace = LH->m_rightFace;
	}
	else {
		m_rightFace = LH->m_leftFace;
	}
	
	if( RF == RL->m_rightFace ) {
		m_leftFace = RL->m_leftFace;
	}
	else {
		m_leftFace = RL->m_rightFace;
	}
			
	//update topology of two legs and two hands
    connect_right_hand_edge( LH );
    connect_left_hand_edge(  LL );
    connect_right_leg_edge(  RH );
    connect_left_leg_edge(   RL );
}


void    WingedEdgeDataStructure::Edge::connect_right_hand_edge( WingedEdgeDataStructure::Edge* const rightHand)
{
    m_rightHand = rightHand;

    if ( m_rightHand != rg_NULL ) {
        if ( m_endVertex == rightHand->m_startVertex ) {
            rightHand->m_rightLeg = this;
        }
        else {
            rightHand->m_leftHand = this;
        }
    }
}



void    WingedEdgeDataStructure::Edge::connect_left_hand_edge(  WingedEdgeDataStructure::Edge* const leftHand)
{
    m_leftHand = leftHand;

    if ( m_leftHand != rg_NULL ) {
        if ( m_endVertex == leftHand->m_startVertex ) {
            leftHand->m_leftLeg   = this;
        }
        else {
            leftHand->m_rightHand = this;
        }
    }
}



void    WingedEdgeDataStructure::Edge::connect_right_leg_edge(  WingedEdgeDataStructure::Edge* const rightLeg)
{
    m_rightLeg = rightLeg;

    if ( m_rightLeg != rg_NULL ) {
        if ( m_startVertex == rightLeg->m_startVertex ) {
            rightLeg->m_leftLeg = this;
        }
        else {
            rightLeg->m_rightHand = this;
        }
    }
}



void    WingedEdgeDataStructure::Edge::connect_left_leg_edge(   WingedEdgeDataStructure::Edge* const leftLeg)
{
    m_leftLeg = leftLeg;

    if ( m_leftLeg != rg_NULL ) {
        if ( m_startVertex == leftLeg->m_startVertex ) {
            leftLeg->m_rightLeg = this;
        }
        else {
            leftLeg->m_leftHand = this;
        }
    }
}


    
rg_INT  WingedEdgeDataStructure::Edge::find_edges_in_star(     list<WingedEdgeDataStructure::Edge*>&   edgesInStar) const
{
    list<WingedEdgeDataStructure::Edge*> edgesInStarOfTwoVertices;
    m_startVertex->find_edges_in_star( edgesInStarOfTwoVertices );
    m_endVertex->find_edges_in_star( edgesInStarOfTwoVertices );

    set<WingedEdgeDataStructure::Edge*> setOfEdgesInStar;
    setOfEdgesInStar.insert( edgesInStarOfTwoVertices.begin(), edgesInStarOfTwoVertices.end() );

    setOfEdgesInStar.erase( const_cast<Edge*>(this) );

    edgesInStar.insert( edgesInStar.begin(), setOfEdgesInStar.begin(), setOfEdgesInStar.end() );

    return edgesInStar.size();
}


///////////////////////////////////////////////////////////////////////////////
//
// class WingedEdgeDataStructure::Vertex : public TopologicalEntity
//
WingedEdgeDataStructure::Vertex::Vertex()
: TopologicalEntity(0),
  m_firstEdge(rg_NULL)
{
}



WingedEdgeDataStructure::Vertex::Vertex(const rg_INT& ID)
: TopologicalEntity(ID),
  m_firstEdge(rg_NULL)
{
}



WingedEdgeDataStructure::Vertex::Vertex(const WingedEdgeDataStructure::Vertex& vertex)
: TopologicalEntity( vertex ),
  m_firstEdge( vertex.m_firstEdge )
{
}



WingedEdgeDataStructure::Vertex::~Vertex()
{
}



WingedEdgeDataStructure::Vertex& WingedEdgeDataStructure::Vertex::operator =(const WingedEdgeDataStructure::Vertex& vertex)
{
    if ( this != &vertex ) {
        TopologicalEntity::operator=( vertex );
        m_firstEdge = vertex.m_firstEdge;
    }

    return *this;
}



bool    WingedEdgeDataStructure::Vertex::is_adjacent_to(const WingedEdgeDataStructure::Vertex* const vertex) const
{   
	WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
	WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    do {
        if ( this == currEdge->start_vertex() ) {
            if ( vertex == currEdge->end_vertex() ) {
                return true;
            }
            currEdge = currEdge->left_leg_edge();
        }
        else if ( this == currEdge->end_vertex() ) {
            if ( vertex == currEdge->start_vertex() ) {
                return true;
            }
            currEdge = currEdge->right_hand_edge();
        }
        else {
            break;
        }
    } while (startEdge != currEdge);

    return false;
}



bool    WingedEdgeDataStructure::Vertex::is_on_boundary() const
{
    bool    isVertexOnBoundary = false;
    list<WingedEdgeDataStructure::Edge*>   incidentEdges;
    find_incident_edges( incidentEdges );

    list<WingedEdgeDataStructure::Edge*>::iterator i_edge;
    for ( i_edge=incidentEdges.begin(); i_edge!=incidentEdges.end(); ++i_edge ) {
        if ( (*i_edge)->is_on_boundary() ) {
            isVertexOnBoundary = true;
            break;
        }
    }

    return isVertexOnBoundary;
}



rg_INT  WingedEdgeDataStructure::Vertex::number_of_adjacent_vertices() const
{
    list<WingedEdgeDataStructure::Vertex*> adjacentVertices;
    return find_adjacent_vertices(adjacentVertices);
}



rg_INT  WingedEdgeDataStructure::Vertex::number_of_incident_edges() const
{
    list<WingedEdgeDataStructure::Edge*> incidentEdges;
    return find_incident_edges(incidentEdges);
}



rg_INT  WingedEdgeDataStructure::Vertex::number_of_incident_faces() const
{
    list<WingedEdgeDataStructure::Face*> incidentFaces;
    return find_incident_faces(incidentFaces);
}



rg_INT  WingedEdgeDataStructure::Vertex::find_adjacent_vertices( list<WingedEdgeDataStructure::Vertex*>& adjacentVertices) const
{
	WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
	WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    unsigned int numAdjacentVertices = 0;
    do {
        if ( this == currEdge->start_vertex() ) {
            adjacentVertices.push_back( currEdge->end_vertex() );
            ++numAdjacentVertices;

            currEdge = currEdge->left_leg_edge();
        }
        else if ( this == currEdge->end_vertex() ) {
            adjacentVertices.push_back( currEdge->start_vertex() );
            ++numAdjacentVertices;

            currEdge = currEdge->right_hand_edge();
        }
        else {
            break;
        }
    } while (startEdge != currEdge);

    return numAdjacentVertices;
}



rg_INT  WingedEdgeDataStructure::Vertex::find_incident_edges(    list<WingedEdgeDataStructure::Edge*>&   incidentEdges) const
{
	WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
	WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    unsigned int numIncidentEdges = 0;
    do {
        incidentEdges.push_back( currEdge );
        ++numIncidentEdges;

        if ( this == currEdge->start_vertex() ) {
            currEdge = currEdge->left_leg_edge();
        }
        else if ( this == currEdge->end_vertex() ) {
            currEdge = currEdge->right_hand_edge();
        }
        else {
            break;
        }
    } while (startEdge != currEdge);

    return numIncidentEdges;
}



rg_INT  WingedEdgeDataStructure::Vertex::find_incident_faces(    list<WingedEdgeDataStructure::Face*>&   incidentFaces) const
{
	WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
	WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    unsigned int numIncidentFaces = 0;
    do {
        if ( this == currEdge->start_vertex() ) {
            if ( currEdge->left_face() != rg_NULL ) {
                incidentFaces.push_back( currEdge->left_face() );
                ++numIncidentFaces;
            }
            currEdge = currEdge->left_leg_edge();
        }
        else if ( this == currEdge->end_vertex() ) {
            if ( currEdge->right_face() != rg_NULL ) {
                incidentFaces.push_back( currEdge->right_face() );
                ++numIncidentFaces;
            }
            currEdge = currEdge->right_hand_edge();
        }
    } while (currEdge != rg_NULL && startEdge != currEdge);

    return numIncidentFaces;
}



rg_INT  WingedEdgeDataStructure::Vertex::find_edges_in_star(     list<WingedEdgeDataStructure::Edge*>&   edgesInStar) const
{
    list<WingedEdgeDataStructure::Face*> facesIncidentToVtx;
    find_incident_faces( facesIncidentToVtx );

    set<WingedEdgeDataStructure::Edge*> setOfEdgesInStar;
    list<WingedEdgeDataStructure::Face*>::iterator i_face;
    for ( i_face=facesIncidentToVtx.begin(); i_face!=facesIncidentToVtx.end(); ++i_face ) {

        list<WingedEdgeDataStructure::Edge*> boundingEdges;
        (*i_face)->find_bounding_edges( boundingEdges );

        setOfEdgesInStar.insert( boundingEdges.begin(), boundingEdges.end() );
    }
    edgesInStar.insert( edgesInStar.begin(), setOfEdgesInStar.begin(), setOfEdgesInStar.end() );

    return edgesInStar.size();
}



rg_INT  WingedEdgeDataStructure::Vertex::find_edges_in_shell(    list<WingedEdgeDataStructure::Edge*>&   edgesInShell) const
{
	WingedEdgeDataStructure::Edge* startEdge = m_firstEdge;
	WingedEdgeDataStructure::Edge* currEdge  = startEdge;

    unsigned int numEdgesInShellOfVtx = 0;
    do {
        WingedEdgeDataStructure::Face* currFace = rg_NULL;
        WingedEdgeDataStructure::Edge* nextEdge = rg_NULL;

        if ( this->is_start_vertex_of( currEdge ) ) {
            currFace = currEdge->left_face();
            nextEdge = currEdge->left_leg_edge();
        }
        else {
            currFace = currEdge->right_face();
            nextEdge = currEdge->right_hand_edge();
        }


        if ( currFace != rg_NULL ) {
            WingedEdgeDataStructure::Edge* currEdgeToTraceFace  = currEdge;
            WingedEdgeDataStructure::Edge* endEdgeToTraceFace   = nextEdge;
            do  {
                if ( !currEdgeToTraceFace->is_incident_to( this ) ) {
                    edgesInShell.push_back( currEdgeToTraceFace );
                    ++numEdgesInShellOfVtx;
                }

                if ( currFace == currEdgeToTraceFace->left_face() ) {
                    currEdgeToTraceFace = currEdgeToTraceFace->left_hand_edge();
                }
                else if ( currFace == currEdgeToTraceFace->right_face() )  {
                    currEdgeToTraceFace = currEdgeToTraceFace->right_leg_edge();
                }
                else {
                    break;
                }
            } while ( currEdgeToTraceFace != rg_NULL && currEdgeToTraceFace != endEdgeToTraceFace );   
        }

        currEdge = nextEdge;

    } while (currEdge != rg_NULL && startEdge != currEdge);


    return numEdgesInShellOfVtx;
}




