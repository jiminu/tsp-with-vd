#include "SimplificationByGarlan.h"

#include <fstream>
#include <algorithm>
#include <set>
using namespace std;


SimplificationByGarlan::SimplificationByGarlan(void)
{
}


SimplificationByGarlan::~SimplificationByGarlan(void)
{
}



void    SimplificationByGarlan::run_simplification( PolygonMeshModel*   model, const double& limitOfQEM )
{
    m_model = model;
    m_limitOfQEM = limitOfQEM;



    list<PolygonMeshModel::Shell*> all_shells;
    m_model->all_shells( all_shells );

    //ofstream fout("simplify_edge.txt");
    //fout << "BEFORE " << endl;
    //fout << "# shells   : " << model->number_of_shells() << endl;
    //fout << "# faces    : " << model->number_of_faces() << endl;
    //fout << "# edges    : " << model->number_of_edges() << endl;
    //fout << "# vertices : " << model->number_of_vertices() << endl << endl;


    list<PolygonMeshModel::Shell*>::iterator i_shell;
    for ( i_shell=all_shells.begin(); i_shell!=all_shells.end(); ++i_shell ) {
        //fout << "shell " << (*i_shell)->getID() << "\t" << (*i_shell)->number_of_faces() << "\t->\t";
        //if ( (*i_shell)->getID() != 347 ) {
        //    continue;
        //}

        simplify_shell( *i_shell );

        //fout << (*i_shell)->number_of_faces() << endl;
    }


    //fout << "AFTER  " << endl;
    //fout << "# shells   : " << model->number_of_shells() << endl;
    //fout << "# faces    : " << model->number_of_faces() << endl;
    //fout << "# edges    : " << model->number_of_edges() << endl;
    //fout << "# vertices : " << model->number_of_vertices() << endl << endl;
    //fout.close();
}



void    SimplificationByGarlan::simplify_shell( PolygonMeshModel::Shell* shell)
{
    //ofstream fout;
    //fout.open("simplify_detail.txt", ios::app);
    //fout << endl;
    //fout << "Shell " << shell->getID() << endl;
    //fout << "\t# faces    : " << shell->number_of_faces() << endl;
    //fout << "\t# edges    : " << shell->number_of_edges() << endl;
    //fout << "\t# vertices : " << shell->number_of_vertices() << endl;



    map<PolygonMeshModel::Edge*, SimplificationByGarlan::Edge>      edgeMap;
    map<PolygonMeshModel::Vertex*, SimplificationByGarlan::Vertex>  vertexMap;

    list<SimplificationByGarlan::Edge*> edgesSortedByCost;
    initialize( shell, edgeMap, vertexMap, edgesSortedByCost );


    //fout << "initial edges sorted by cost " << endl;
    //int count = 0;
    //list<SimplificationByGarlan::Edge*>::iterator it_sime;
    //for ( it_sime=edgesSortedByCost.begin(); it_sime!=edgesSortedByCost.end(); ++it_sime ) {
    //    SimplificationByGarlan::Edge* edge = *it_sime;
    //    PolygonMeshModel::Edge* model_edge  = edge->pedge();

    //    fout << ++count << "\tedge " << model_edge->getID() << "\tin shell " << shell->getID();
    //    fout << "\tcost = " << edge->cost();
    //    fout << endl;
    //}
    //fout << endl << endl;


    
    while ( !edgesSortedByCost.empty() && shell->number_of_faces() > 4 ) {
        SimplificationByGarlan::Edge* edgeToCollapse = edgesSortedByCost.front();
        edgesSortedByCost.pop_front();


        PolygonMeshModel::Edge* model_edge  = edgeToCollapse->pedge();
        rg_Point3D              positionOfNewVertex = edgeToCollapse->position_for_new_vertex();

        //fout << "\t" << "collapse e_" << model_edge->getID() << "\tin shell " << shell->getID();
        //fout << "\t" << "cost : " << edgeToCollapse->cost();        
        //fout << "\t" << edgesSortedByCost.size();

        if ( edgeToCollapse->cost() > m_limitOfQEM ) {
            //fout << "\t break !!! " << endl;
            break;
        }


        if ( model_edge->is_on_boundary() ) {
            //fout << "\t on boundary" << endl;
            continue;
        }

        if ( m_model->does_nonmanifold_happen_by_edge_collapse( model_edge ) ) {
            //fout << "\t nonmanifold " << endl;
            continue;
        }

        if ( m_model->does_folding_face_happen_by_edge_collapse( model_edge, positionOfNewVertex ) ) 
        {
            //fout << "\t folding" << endl;
            continue;
        }

        //fout << "\t collapse " << endl;
        //fout << "  sv v_" << model_edge->start_vertex()->getID();
        //fout << "  ev v_" << model_edge->end_vertex()->getID();
        //fout << "  rf f_" << model_edge->right_face()->getID();
        //fout << "  lf f_" << model_edge->left_face()->getID();
        //fout << "  rh e_" << model_edge->right_hand_edge()->getID();
        //fout << "  lh e_" << model_edge->left_hand_edge()->getID();
        //fout << "  rl e_" << model_edge->right_leg_edge()->getID();
        //fout << "  ll e_" << model_edge->left_leg_edge()->getID();

        list<SimplificationByGarlan::Edge*> edgesInfluencedByEdgeCollapse;
        collect_edges_influenced_by_edge_collapse( 
                            edgeToCollapse, edgeMap, edgesInfluencedByEdgeCollapse);

        remove_influenced_edges_from_sorted_list(
                            edgeToCollapse, edgeMap, edgesSortedByCost, edgesInfluencedByEdgeCollapse);
        
        SimplificationByGarlan::Vertex* newVertex   = edgeToCollapse->end_vertex();


        m_model->collapse_edge( model_edge, positionOfNewVertex );


        newVertex->set_Q( edgeToCollapse );
        update_influenced_edges( vertexMap, edgesInfluencedByEdgeCollapse);
        insert_influenced_edges_into_sorted_list( edgesSortedByCost, edgesInfluencedByEdgeCollapse);


        //fout << "influenced edges" << endl;
        //list<SimplificationByGarlan::Edge*>::iterator i_edge;
        //for ( i_edge=edgesInfluencedByEdgeCollapse.begin(); i_edge!=edgesInfluencedByEdgeCollapse.end(); ++i_edge ) {
        //    fout << "\t" << "collapse e_" << (*i_edge)->pedge()->getID() << "\tin shell " << shell->getID();
        //    fout << "\t" << "cost : " << (*i_edge)->cost();        
        //    fout << "\t" << endl;
        //}
        //fout << endl;
        //fout << endl;

        //fout << "edges sorted by cost" << endl;
        //for ( i_edge=edgesSortedByCost.begin(); i_edge!=edgesSortedByCost.end(); ++i_edge ) {
        //    fout << "\t" << "collapse e_" << (*i_edge)->pedge()->getID() << "\tin shell " << shell->getID();
        //    fout << "\t" << "cost : " << (*i_edge)->cost();        
        //    fout << "\t" << endl;
        //}
        //fout << endl;
        //fout << endl;

    }
    //fout << endl;

    //fout << "Result of Collapsing Edges" << endl;
    //fout << "\t# faces    : " << shell->number_of_faces() << endl;
    //fout << "\t# edges    : " << shell->number_of_edges() << endl;
    //fout << "\t# vertices : " << shell->number_of_vertices() << endl;
    //fout << endl;

    //fout << "Remaining edges : " << edgesSortedByCost.size() << endl;
    ////list<SimplificationByGarlan::Edge*>::iterator it_sime;
    //for ( it_sime=edgesSortedByCost.begin(); it_sime!=edgesSortedByCost.end(); ++it_sime ) {
    //    SimplificationByGarlan::Edge* edge = *it_sime;
    //    PolygonMeshModel::Edge* model_edge  = edge->pedge();

    //    fout << "\tedge " << model_edge->getID() << "\tin shell " << shell->getID();
    //    fout << "\tcost = " << edge->cost();
    //    fout << endl;
    //}
    //fout << endl << endl;
    //fout.close();


    m_model->release_deleted_entities();

    
}


        

void    SimplificationByGarlan::initialize( PolygonMeshModel::Shell* shell,
                            map<PolygonMeshModel::Edge*, SimplificationByGarlan::Edge>& edgeMap,
                            map<PolygonMeshModel::Vertex*, SimplificationByGarlan::Vertex>& vertexMap,
                            list<SimplificationByGarlan::Edge*>& edgesSortedByCost)
{
    list<PolygonMeshModel::Vertex*> verticesInShell;
    int numVertices = shell->find_bounding_vertices( verticesInShell );

    list<PolygonMeshModel::Vertex*>::iterator i_vtx;
    for ( i_vtx=verticesInShell.begin(); i_vtx!=verticesInShell.end(); ++i_vtx ) {
        PolygonMeshModel::Vertex* currVtx = *i_vtx;

        SimplificationByGarlan::Vertex simVtx( currVtx );
        vertexMap.insert( make_pair(currVtx, simVtx) );
    }



    list<PolygonMeshModel::Edge*> edgesInShell;
    int numEdges = shell->find_bounding_edges( edgesInShell );

    list<PolygonMeshModel::Edge*>::iterator i_edge;
    for ( i_edge=edgesInShell.begin(); i_edge!=edgesInShell.end(); ++i_edge ) {
        PolygonMeshModel::Edge* currEdge = *i_edge;

        SimplificationByGarlan::Vertex* startSimVtx = &(vertexMap.find( currEdge->start_vertex() )->second);
        SimplificationByGarlan::Vertex* endSimVtx   = &(vertexMap.find( currEdge->end_vertex() )->second);

        SimplificationByGarlan::Edge simEdge( startSimVtx, endSimVtx, currEdge );

        edgeMap.insert( make_pair(currEdge, simEdge) );
    }

        
    vector<SimplificationByGarlan::Edge*> edges4Sorting(numEdges);
    int i=0;
    map<PolygonMeshModel::Edge*, SimplificationByGarlan::Edge>::iterator it_edge;
    for ( it_edge=edgeMap.begin(); it_edge!=edgeMap.end(); ++it_edge, ++i ) {
        edges4Sorting[i] = &(it_edge->second);
    }
    sort( edges4Sorting.begin(), edges4Sorting.end(), Edge::LessThan );


    edgesSortedByCost.insert( edgesSortedByCost.begin(), edges4Sorting.begin(), edges4Sorting.end() );

    
    
}



void    SimplificationByGarlan::collect_edges_influenced_by_edge_collapse(
                            SimplificationByGarlan::Edge* edgeToCollapse,
                            const map<PolygonMeshModel::Edge*, SimplificationByGarlan::Edge>& edgeMap,
                            list<SimplificationByGarlan::Edge*>& edgesInfluencedByEdgeCollapse)
{
    PolygonMeshModel::Edge* modelEdge = edgeToCollapse->pedge();

    list<PolygonMeshModel::Edge*> edgesInStar;
    modelEdge->find_edges_in_star( edgesInStar );

    list<PolygonMeshModel::Edge*>::iterator i_edge;
    for ( i_edge=edgesInStar.begin(); i_edge!=edgesInStar.end(); ++i_edge ) {
        edgesInfluencedByEdgeCollapse.push_back( const_cast<Edge*>( &(edgeMap.find( *i_edge )->second) ) );
    }
}

        
        
void    SimplificationByGarlan::remove_influenced_edges_from_sorted_list(
                            SimplificationByGarlan::Edge*               edgeToCollapse,
                            const map<PolygonMeshModel::Edge*, SimplificationByGarlan::Edge>& edgeMap,
                            list<SimplificationByGarlan::Edge*>&        edgesSortedByCost, 
                            list<SimplificationByGarlan::Edge*>&        edgesInfluencedByEdgeCollapse)
{
    list<SimplificationByGarlan::Edge*>::const_iterator i_edge;
    for ( i_edge=edgesInfluencedByEdgeCollapse.begin(); i_edge!=edgesInfluencedByEdgeCollapse.end(); ++i_edge ) {
        edgesSortedByCost.remove( *i_edge );
    }

    PolygonMeshModel::Edge* modelEdge = edgeToCollapse->pedge();
    PolygonMeshModel::Edge* edgeToBeRemoved[2] = { modelEdge->right_leg_edge(), modelEdge->left_leg_edge() };
    for ( int i=0; i<2; ++i ) {
        edgesInfluencedByEdgeCollapse.remove( const_cast<Edge*>( &(edgeMap.find( edgeToBeRemoved[i] )->second) ) );
    }
}


void    SimplificationByGarlan::update_influenced_edges(
                            const map<PolygonMeshModel::Vertex*, SimplificationByGarlan::Vertex>& vertexMap,
                            list<SimplificationByGarlan::Edge*>& edgesInfluencedByEdgeCollapse)
{
    list<SimplificationByGarlan::Edge*>::iterator i_edge;
    for ( i_edge=edgesInfluencedByEdgeCollapse.begin(); i_edge!=edgesInfluencedByEdgeCollapse.end(); ++i_edge ) {
        SimplificationByGarlan::Edge* currEdge = *i_edge;

        PolygonMeshModel::Edge* model_edge = currEdge->pedge();
        SimplificationByGarlan::Vertex* startVtx = const_cast<Vertex*>( &(vertexMap.find( model_edge->start_vertex() )->second) );
        SimplificationByGarlan::Vertex* endVtx = const_cast<Vertex*>( &(vertexMap.find( model_edge->end_vertex() )->second) );

        currEdge->set_start_vertex( startVtx );
        currEdge->set_end_vertex(   endVtx   );
        currEdge->initialize();
    }
}



void    SimplificationByGarlan::insert_influenced_edges_into_sorted_list(
                            list<SimplificationByGarlan::Edge*>& edgesSortedByCost, 
                            list<SimplificationByGarlan::Edge*>& edgesInfluencedByEdgeCollapse)
{
    vector<SimplificationByGarlan::Edge*> influencedEdgesToSort;
    influencedEdgesToSort.insert( influencedEdgesToSort.begin(), edgesInfluencedByEdgeCollapse.begin(), edgesInfluencedByEdgeCollapse.end() );
    sort( influencedEdgesToSort.begin(), influencedEdgesToSort.end(), Edge::LessThan );

    int numEdgeToInsert = influencedEdgesToSort.size();
    int pos = 0;

    list<SimplificationByGarlan::Edge*>::iterator i_edge=edgesSortedByCost.begin();
    while ( pos < numEdgeToInsert && i_edge!=edgesSortedByCost.end() ) {
        if ( (*i_edge)->cost() > influencedEdgesToSort[pos]->cost() ) {
            edgesSortedByCost.insert( i_edge, influencedEdgesToSort[pos] );
            ++pos;
        }
        else {
            ++i_edge;
        }
    }
    for ( int i=pos; i<numEdgeToInsert; ++i ) {
        edgesSortedByCost.push_back( influencedEdgesToSort[i] );
    }
    
}

///////////////////////////////////////////////////////////////////////////////
//
// class SimplificationByGarlan::Edge 
//
SimplificationByGarlan::Edge::Edge()
: m_pedge( rg_NULL )
{
    m_vertex[0] = rg_NULL;
    m_vertex[1] = rg_NULL;
}



SimplificationByGarlan::Edge::Edge(Vertex* start, Vertex* end, PolygonMeshModel::Edge* pedge)
: m_pedge(pedge)
{
    m_vertex[0] = start;
    m_vertex[1] = end;

    initialize();
}



SimplificationByGarlan::Edge::Edge(const SimplificationByGarlan::Edge& edge)
: m_pedge( edge.m_pedge )
, m_cost(  edge.m_cost )
, m_position4newVertex( edge.m_position4newVertex )
{
    m_vertex[0] = edge.m_vertex[0];
    m_vertex[1] = edge.m_vertex[1];
}



SimplificationByGarlan::Edge::~Edge()
{
}



SimplificationByGarlan::Edge& SimplificationByGarlan::Edge::operator =(const SimplificationByGarlan::Edge& edge)
{
    if ( this != &edge ) {
        m_vertex[0] = edge.m_vertex[0];
        m_vertex[1] = edge.m_vertex[1];
        m_pedge     = edge.m_pedge;
        m_cost      = edge.m_cost;
        m_position4newVertex = edge.m_position4newVertex;
    }

    return *this;
}



void    SimplificationByGarlan::Edge::initialize()
{
    if ( m_vertex[0] == rg_NULL || m_vertex[1] == rg_NULL || m_pedge == rg_NULL ) {
        return;
    }



    //////////////////////////////////////////////////////////////////
    //  calculate Q for edge
    double Q[10];
    for ( int i=0; i<10; ++i ) {
        Q[i] = m_vertex[0]->Q(i) + m_vertex[1]->Q(i);
    }


    //////////////////////////////////////////////////////////////////
    //  calculate position of new vertex
	//A = vQv
	//	dA/dx =      2q00x + (q01+q10)y + (q02+20)z + (q03+q30) 
	//	      =      2q00x +      2q01y +     2q02z +      2q03
	//	       
	//	dA/dy = (q10+q01)x +      2q11y + (q12+21)z + (q13+q31)  
	//	      =      2q01x +      2q11y +     2q12z +      2q13
	//	       
	//	dA/dz = (q20+q02)x + (q21+q12)y +     2q22z + (q23+q32)  
	//	      =      2q02x +      2q12y +     2q22z +      2q23		 

    double dQ[9] = { 2*Q[0], 2*Q[1], 2*Q[2], 2*Q[1], 2*Q[4], 2*Q[5], 2*Q[2], 2*Q[5], 2*Q[7] };
    double c[3]  = { 2*Q[3], 2*Q[6], 2*Q[8] };

	double determinent =   dQ[2]*dQ[4]*dQ[6] - dQ[1]*dQ[5]*dQ[6] - dQ[2]*dQ[3]*dQ[7] 
                         + dQ[0]*dQ[5]*dQ[7] + dQ[1]*dQ[3]*dQ[8] - dQ[0]*dQ[4]*dQ[8];
		
	if ( rg_ZERO(determinent) ) {
        rg_Point3D startPoint = m_pedge->start_vertex()->coordinate();
		rg_Point3D endPoint   = m_pedge->end_vertex()->coordinate();
			
		m_position4newVertex = (startPoint + endPoint)/2.0;
	}
	else {
		double x = -(  dQ[5]*dQ[7]*c[0] - dQ[4]*dQ[8]*c[0] - dQ[2]*dQ[7]*c[1]
			         + dQ[1]*dQ[8]*c[1] + dQ[2]*dQ[4]*c[2] - dQ[1]*dQ[5]*c[2]);
			           
		double y =  (  dQ[5]*dQ[6]*c[0] - dQ[3]*dQ[8]*c[0] - dQ[2]*dQ[6]*c[1]
			         + dQ[0]*dQ[8]*c[1] + dQ[2]*dQ[3]*c[2] - dQ[0]*dQ[5]*c[2]);
			           
		double z = -(  dQ[4]*dQ[6]*c[0] - dQ[3]*dQ[7]*c[0] - dQ[1]*dQ[6]*c[1]
			         + dQ[0]*dQ[7]*c[1] + dQ[1]*dQ[3]*c[2] - dQ[0]*dQ[4]*c[2]); 
            
        m_position4newVertex.setPoint(x/determinent, y/determinent, z/determinent);
    }   


    //////////////////////////////////////////////////////////////////
    //  calculate cost of edge
	double x = m_position4newVertex.getX();
	double y = m_position4newVertex.getY();
	double z = m_position4newVertex.getZ();
		
	m_cost   =   Q[0]*x*x   + Q[4]*y*y   + Q[7]*z*z   
               + 2*Q[1]*x*y + 2*Q[5]*y*z + 2*Q[2]*x*z 
               + 2*Q[3]*x   + 2*Q[6]*y   + 2*Q[8]*z
               + Q[9];

}



    
bool SimplificationByGarlan::Edge::LessThan( SimplificationByGarlan::Edge* e1, SimplificationByGarlan::Edge* e2 ) 
{
    return (e1->cost() < e2->cost()) ? true : false;
}


///////////////////////////////////////////////////////////////////////////////
//
// class SimplificationByGarlan::Vertex 
//
SimplificationByGarlan::Vertex::Vertex()
: m_pvertex(rg_NULL)
{
	for(int i=0; i<10; ++i)  {
		m_Q[i] = 0.0;
	}
}



SimplificationByGarlan::Vertex::Vertex(PolygonMeshModel::Vertex* pvertex)
: m_pvertex(pvertex)
{
	for(int i=0; i<10; ++i)  {
		m_Q[i] = 0.0;
	}

    initialize();
}



SimplificationByGarlan::Vertex::Vertex(const SimplificationByGarlan::Vertex& vertex)
: m_pvertex(vertex.m_pvertex)
{
	for(int i=0; i<10; ++i)  {
        m_Q[i] = vertex.m_Q[i];
	}
}



SimplificationByGarlan::Vertex::~Vertex()
{
}



SimplificationByGarlan::Vertex& SimplificationByGarlan::Vertex::operator =(const SimplificationByGarlan::Vertex& vertex)
{
    if ( this != &vertex ) {
        m_pvertex = vertex.m_pvertex;
	    for(int i=0; i<10; ++i)  {
            m_Q[i] = vertex.m_Q[i];
	    }
    }

    return *this;
}



void    SimplificationByGarlan::Vertex::set_Q(SimplificationByGarlan::Edge* edge)
{
    for ( int i=0; i<10; ++i ) {
        m_Q[i] = edge->start_vertex()->m_Q[i] + edge->end_vertex()->m_Q[i];
    }
}



void    SimplificationByGarlan::Vertex::initialize()
{
	for(int i=0; i<10; ++i)  {
		m_Q[i] = 0.0;
	}

    list<PolygonMeshModel::Face*> facesIncidentToVtx;
    m_pvertex->find_incident_faces( facesIncidentToVtx );

    list<PolygonMeshModel::Face*>::iterator i_face;
    for ( i_face=facesIncidentToVtx.begin(); i_face!=facesIncidentToVtx.end(); ++i_face ) {
        PolygonMeshModel::Face* currFace = *i_face;
        Plane                   facePlane = currFace->plane();

		double a = facePlane.getNormal().getX();
		double b = facePlane.getNormal().getY();
		double c = facePlane.getNormal().getZ();
		double d = -facePlane.getDistanceFromOrigin();
			
		m_Q[0] += a*a;
		m_Q[1] += a*b;
		m_Q[2] += a*c;
		m_Q[3] += a*d;
		m_Q[4] += b*b;
		m_Q[5] += b*c;
		m_Q[6] += b*d;
		m_Q[7] += c*c;
		m_Q[8] += c*d;
		m_Q[9] += d*d;
    }
}



