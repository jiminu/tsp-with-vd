#include "Graph.h"

#include <list>
#include <set>
#include <map>
using namespace std;



Graph::Graph()
: m_number_of_nodes(0)
, m_number_of_arcs(0)
{
}



Graph::Graph(const Graph& graph)
: m_number_of_nodes(graph.m_number_of_nodes)
, m_number_of_arcs(graph.m_number_of_arcs)
{
}



Graph::~Graph()
{
    clear();
}



void    Graph::clear()
{
    for (list<Node*>::iterator i_node = m_nodes.begin(); i_node != m_nodes.end(); ++i_node) {
        delete (*i_node);
    }

    for (list<Arc*>::iterator i_arc = m_arcs.begin(); i_arc != m_arcs.end(); ++i_arc) {
        delete (*i_arc);
    }
}



Graph& Graph::operator =(const Graph& graph)
{
    if (this != &graph) {
    }

    return *this;
}



const list< Graph::Node* >&    Graph::get_all_nodes() const
{
    return m_nodes;
}



const list< Graph::Arc* >&     Graph::get_all_arcs() const
{
    return m_arcs;
}



unsigned int            Graph::number_of_nodes() const
{
    return m_number_of_nodes;
}



unsigned int            Graph::number_of_arcs() const
{
    return m_number_of_arcs;
}



Graph::Node*    Graph::create_node()
{
    int ID = 0;
    if (!m_nodes.empty()) {
        ID = m_nodes.back()->getID() + 1;
    }
    m_nodes.push_back(new Node(ID));
    ++m_number_of_nodes;

    return m_nodes.back();
}



Graph::Arc*     Graph::create_arc()
{
    int ID = 0;
    if (!m_arcs.empty()) {
        ID = m_arcs.back()->getID() + 1;
    }
    m_arcs.push_back(new Arc(ID));
    ++m_number_of_arcs;

    return m_arcs.back();
}



void    Graph::remove_node(const Graph::Node* const node)
{
    if (node != rg_NULL) {
        m_nodes.remove(const_cast<Node*>(node));
        --m_number_of_nodes;
        delete node;
    }
}



void    Graph::remove_arc(const Graph::Arc* const arc)
{
    if (arc != rg_NULL) {
        m_arcs.remove(const_cast<Arc*>(arc));
        --m_number_of_arcs;
        delete arc;
    }
}



void    Graph::divide_into_subgraphs_connected_by_arcs(list<Graph>& subgraphs)
{
    list< set<Arc*> > componentsOfConnectedArcs;

    set<Node*> visitedNodes;
    for (list<Node*>::iterator i_node = m_nodes.begin(); i_node != m_nodes.end(); ++i_node) {
        Node* firstNode = *i_node;

        if (visitedNodes.find(firstNode) != visitedNodes.end()) {
            continue;
        }

        componentsOfConnectedArcs.push_back(set<Arc*>());
        set<Arc*>* currComponent = &componentsOfConnectedArcs.back();

        list<Node*> nodeQueue;
        nodeQueue.push_back(firstNode);
        while (!nodeQueue.empty()) {
            Node* currNode = nodeQueue.front();
            nodeQueue.pop_front();

            visitedNodes.insert(currNode);

            const list< Arc* >& incidentArcs = currNode->incident_arcs();
            for (list<Arc*>::const_iterator i_arc = incidentArcs.begin(); i_arc != incidentArcs.end(); ++i_arc) {
                Arc* currArc = *i_arc;

                currComponent->insert(currArc);

                Node* neighbor = currArc->opposite_node(currNode);
                if (visitedNodes.find(neighbor) == visitedNodes.end()) {
                    nodeQueue.push_back(neighbor);
                }
            }
        }
    }


    for (list< set<Arc*> >::iterator i_com = componentsOfConnectedArcs.begin(); i_com != componentsOfConnectedArcs.end(); ++i_com) {
        list<Arc*> arcsInSubgraph;
        arcsInSubgraph.insert(arcsInSubgraph.begin(), i_com->begin(), i_com->end());

        subgraphs.push_back(Graph());
        Graph* subgraph = &subgraphs.back();

        make_subgraph(arcsInSubgraph, *subgraph);
    }

}



void    Graph::make_subgraph(const list<Arc*>& arcs, Graph& subgraph)
{
    //  make nodes of subgraph
    map<Node*, Node*> nodeMap;

    list<Arc*>::const_iterator i_arc;
    for (i_arc = arcs.begin(); i_arc != arcs.end(); ++i_arc) {
        Node* node[2] = { (*i_arc)->start_node(), (*i_arc)->end_node() };
        for (int i = 0; i < 2; ++i) {
            if (nodeMap.find(node[i]) == nodeMap.end()) {
                Node* nodeOfSubgraph = subgraph.create_node();
                nodeOfSubgraph->set_node_data(node[i]->node_data());

                nodeMap.insert(make_pair(node[i], nodeOfSubgraph));                
            }
        }
    }


    for (i_arc = arcs.begin(); i_arc != arcs.end(); ++i_arc) {
        Node* node[2]           = { (*i_arc)->start_node(),        (*i_arc)->end_node() };
        Node* nodeOfSubgraph[2] = { nodeMap.find(node[0])->second, nodeMap.find(node[1])->second };

        Arc* arcOfSubgraph = subgraph.create_arc();
        arcOfSubgraph->set_start_node(nodeOfSubgraph[0]);
        arcOfSubgraph->set_end_node(nodeOfSubgraph[1]);

        nodeOfSubgraph[0]->add_arc(arcOfSubgraph);
        nodeOfSubgraph[1]->add_arc(arcOfSubgraph);
    }
}



///////////////////////////////////////////////////////////////////////////////
//
// class Graph::Node 
//

Graph::Node::Node()
: TopologicalEntity()
, m_node_data(rg_NULL)
{
}



Graph::Node::Node(const int& ID)
: TopologicalEntity(ID)
, m_node_data(rg_NULL)
{
}



Graph::Node::Node(const Graph::Node& node)
: TopologicalEntity(node)
, m_node_data(node.m_node_data)
{
    m_incident_arcs = node.m_incident_arcs;
}



Graph::Node::~Node()
{
}



Graph::Node& Graph::Node::operator =(const Graph::Node& node)
{
    if (this != &node) {
        TopologicalEntity::operator=(node);
        m_incident_arcs = node.m_incident_arcs;
        m_node_data = node.m_node_data;
    }

    return *this;
}



unsigned int            Graph::Node::number_of_incident_arcs() const
{
    return m_incident_arcs.size();
}



const list< Graph::Arc* >&     Graph::Node::incident_arcs() const
{
    return m_incident_arcs;
}



void                    Graph::Node::add_arc(Graph::Arc* arc)
{
    m_incident_arcs.push_back(arc);
}



void                    Graph::Node::remove_arc(Graph::Arc* arc)
{
    m_incident_arcs.remove(arc);
}



void*   Graph::Node::node_data() const
{
    return m_node_data;
}


void    Graph::Node::set_node_data(void* nodeData)
{
    m_node_data = nodeData;
}


///////////////////////////////////////////////////////////////////////////////
//
// class Graph::Arc 
//
Graph::Arc::Arc()
: TopologicalEntity()
, m_start_node(rg_NULL)
, m_end_node(rg_NULL)
{
}



Graph::Arc::Arc(const int& ID)
: TopologicalEntity(ID)
, m_start_node(rg_NULL)
, m_end_node(rg_NULL)
{
}



Graph::Arc::Arc(const Graph::Arc& arc)
: TopologicalEntity(arc)
, m_start_node(arc.m_start_node)
, m_end_node(arc.m_end_node)
{
}



Graph::Arc::~Arc()
{
}



Graph::Arc&    Graph::Arc::operator =(const Graph::Arc& arc)
{
    if (this != &arc) {
        TopologicalEntity::operator=(arc);
        m_start_node = arc.m_start_node;
        m_end_node   = arc.m_end_node;
    }

    return *this;
}



Graph::Node*   Graph::Arc::start_node() const
{
    return m_start_node;
}



Graph::Node*   Graph::Arc::end_node()  const
{
    return m_end_node;
}



void    Graph::Arc::set_start_node(Graph::Node* node)
{
    m_start_node = node;
}



void    Graph::Arc::set_end_node(Graph::Node* node)
{
    m_end_node = node;
}



Graph::Node*   Graph::Arc::opposite_node(Graph::Node* node)
{
    if (node == m_start_node) {
        return m_end_node;
    }
    else if (node == m_end_node) {
        return m_start_node;
    }
    else {
        return rg_NULL;
    }
}

