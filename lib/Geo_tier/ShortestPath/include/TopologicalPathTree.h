#ifndef _TOPOLOGICAL_PATH_TREE_
#define _TOPOLOGICAL_PATH_TREE_

#include "TPTNode.h"
#include "EdgeBU2D.h"
#include <unordered_map>
#include <vector>
using namespace std;

class TopologicalPathTree
{
private:
    list<TPTNode>  m_Nodes;
    list<TPTNode*> m_LeafNodes;
    list<TPTNode*> m_SortedLeafNodes;

    TPTNode*       m_RootNode;

public:
    //constructor
    TopologicalPathTree();
    TopologicalPathTree(const TopologicalPathTree& TPT);

    //destructor;
    ~TopologicalPathTree();

    //equal operation
    TopologicalPathTree& operator=(const TopologicalPathTree& TPT);

    //setter
    inline void set_nodes(const list<TPTNode>& nodes)   { m_Nodes.clear(); m_Nodes = nodes; };
    inline void set_leaf_nodes(const list<TPTNode*>& leafNodes) { m_LeafNodes.clear(); m_LeafNodes = leafNodes; };
    inline void set_root(TPTNode* rootNode) { m_RootNode = rootNode; };

    //getter
    void get_nodes(list<TPTNode*>& nodes) const;
    inline void get_leaf_nodes(list<TPTNode*>& leafNodes) const { leafNodes = m_LeafNodes; };

    inline TPTNode* get_root_node() const { return m_RootNode; };

    //function
    TPTNode* create_node(const TPTNode& node);
    void     add_leaf_node(TPTNode* node);
    bool     does_this_node_have_ancestor_having_this_face(TPTNode* node, FaceBU2D* face) const;
    void     get_topological_path_from_this_leaf_node(TPTNode* leafNode, list<FaceBU2D*>& topologicalQTFacePath) const;
    void     get_topological_paths_in_non_decreasing_order_by(const int& numOfPaths, vector<list<FaceBU2D*>>& topologicalQTFacePaths);
 
    static bool compare_two_leaf_nodes_in_non_decreasing_order(TPTNode* node1, TPTNode* node2);

private:
    void copy_from(const TopologicalPathTree& TPT);
    void create_nodes_from_prev_nodes(const list<TPTNode*>& prevNodes, unordered_map<TPTNode*, TPTNode*>& nodeLinkerFromPrevToCurr);
    void modify_leaf_node(const list<TPTNode*>& prevNodes, list<TPTNode*>& currNodes, const unordered_map<TPTNode*, TPTNode*>& nodeLinkerFromPrevToCurr);
    void modify_node_topology(const unordered_map<TPTNode*, TPTNode*>& nodeLinkerFromPrevToCurr);
    void modify_root_node(TPTNode* oldRootNode, const unordered_map<TPTNode*, TPTNode*>& nodeLinkerFromPrevToCurr);
    void sort_leaf_node_by_its_depth();

    void clear();
};


#endif