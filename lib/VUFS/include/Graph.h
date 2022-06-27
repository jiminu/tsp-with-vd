#ifndef _GRAPH_H
#define _GRAPH_H


#include "TopologicalEntity.h"

#include <list>
using namespace std;



class Graph
{
public:
    class Node;
    class Arc;

protected:
    list< Node* >   m_nodes;
    list< Arc* >    m_arcs;

    unsigned int    m_number_of_nodes;
    unsigned int    m_number_of_arcs;

public:
    Graph();
    Graph(const Graph& graph);
    ~Graph();

    void    clear();
    Graph& operator =(const Graph& graph);

    const list< Node* >&    get_all_nodes() const;
    const list< Arc* >&     get_all_arcs() const;

    unsigned int            number_of_nodes() const;
    unsigned int            number_of_arcs() const;

    virtual Graph::Node*    create_node();
    virtual Graph::Arc*     create_arc();

    void                    remove_node(const Graph::Node* const node);
    void                    remove_arc(const Graph::Arc* const arc);

    void                    divide_into_subgraphs_connected_by_arcs(list<Graph>& subgraphs);

    void                    make_subgraph(const list<Arc*>& arcs, Graph& subgraph);
};


class Graph::Node : public TopologicalEntity
{
private:
    list<Arc*>  m_incident_arcs;
    void*       m_node_data;

public:
    Node();
    Node(const int& ID);
    Node(const Node& node);
    ~Node();
    
    Node& operator =(const Node& node);

    unsigned int            number_of_incident_arcs() const;
    const list< Arc* >&     incident_arcs() const;

    void                    add_arc(Arc* arc);
    void                    remove_arc(Arc* arc);

    void*                   node_data() const;
    void                    set_node_data(void* nodeData);
};


class Graph::Arc : public TopologicalEntity
{
private:
    Node*   m_start_node;
    Node*   m_end_node;

public:
    Arc();
    Arc(const int& ID);
    Arc(const Arc& arc);
    ~Arc();

    Arc& operator =(const Arc& arc);

    Node*   start_node() const;
    Node*   end_node()  const;

    void    set_start_node(Node* node);
    void    set_end_node(Node* node);

    Node*   opposite_node(Node* node);
};


#endif


