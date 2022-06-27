#include "TopologicalPathTree.h"


TopologicalPathTree::TopologicalPathTree()
{
    m_RootNode = NULL;
}

TopologicalPathTree::TopologicalPathTree(const TopologicalPathTree& TPT)
{
    copy_from(TPT);
}

TopologicalPathTree::~TopologicalPathTree()
{
}

TopologicalPathTree& TopologicalPathTree::operator=(const TopologicalPathTree& TPT)
{
    if (this == &TPT)
    {
        return *this;
    }

    clear();
    copy_from(TPT);

    return *this;
}

void TopologicalPathTree::get_nodes(list<TPTNode*>& nodes) const
{
    for (list<TPTNode>::const_iterator it_Node = m_Nodes.begin(); it_Node != m_Nodes.end(); it_Node++)
    {
        TPTNode* currNode = const_cast<TPTNode*>(&(*it_Node));
        nodes.push_back(currNode);
    }
}

TPTNode* TopologicalPathTree::create_node(const TPTNode& node)
{
    m_Nodes.push_back(TPTNode(node, m_Nodes.size()));

    TPTNode* newNode = (&(*m_Nodes.rbegin()));

    if (m_Nodes.size() == 1)
    {
        set_root(newNode);
    }

    return newNode;
}


void TopologicalPathTree::add_leaf_node(TPTNode* node)
{
    m_LeafNodes.push_back(node);
}

void TopologicalPathTree::clear()
{
    m_Nodes.clear();
    m_LeafNodes.clear();
    m_SortedLeafNodes.clear();

    m_RootNode = NULL;
}

void TopologicalPathTree::copy_from(const TopologicalPathTree& TPT)
{
    clear();

    unordered_map<TPTNode*, TPTNode*> NodeLinkerFromPrevToCurr;
    
    list<TPTNode*> prevNodes;
    TPT.get_nodes(prevNodes);

    //1. make link for node
    create_nodes_from_prev_nodes(prevNodes, NodeLinkerFromPrevToCurr);

    //2. change node topology from prev to curr
    modify_node_topology(NodeLinkerFromPrevToCurr);

    //3. change node topology from prev to curr
    modify_root_node(TPT.get_root_node(), NodeLinkerFromPrevToCurr);

    //4.  change node topology from prev to curr
    modify_leaf_node(TPT.m_LeafNodes, m_LeafNodes, NodeLinkerFromPrevToCurr);
    modify_leaf_node(TPT.m_SortedLeafNodes, m_SortedLeafNodes, NodeLinkerFromPrevToCurr);
}

void TopologicalPathTree::create_nodes_from_prev_nodes(const list<TPTNode*>& prevNodes, unordered_map<TPTNode*, TPTNode*>& nodeLinkerFromPrevToCurr)
{

    for (list<TPTNode*>::const_iterator it_Node = prevNodes.begin(); it_Node != prevNodes.end(); it_Node++)
    {
        TPTNode* prevNode = *it_Node;
        TPTNode* newNode  = create_node(*prevNode);
        nodeLinkerFromPrevToCurr.insert(make_pair(prevNode, newNode));

    }
}


void TopologicalPathTree::modify_leaf_node(const list<TPTNode*>& prevNodes, list<TPTNode*>& currNodes, const unordered_map<TPTNode*, TPTNode*>& nodeLinkerFromPrevToCurr)
{
    for (list<TPTNode*>::const_iterator it_Node = prevNodes.begin(); it_Node != prevNodes.end(); it_Node++)
    {
        TPTNode* prevNode = *it_Node;
        currNodes.push_back(nodeLinkerFromPrevToCurr.at(prevNode));
    }
}


void TopologicalPathTree::modify_node_topology(const unordered_map<TPTNode*, TPTNode*>& nodeLinkerFromPrevToCurr)
{
    //1. Chamge in m_SortedLeafNodes
    for (list<TPTNode>::iterator it_Node = m_Nodes.begin(); it_Node != m_Nodes.end(); it_Node++)
    {
        TPTNode* newNode = &(*it_Node);

        //1. children node
        list <TPTNode*> oldChildrenNode;
        newNode->get_children_nodes(oldChildrenNode);

        list <TPTNode*> newChildrenNode;
        
        for (list<TPTNode*>::const_iterator it_OldNode = oldChildrenNode.begin(); it_OldNode != oldChildrenNode.end(); it_OldNode++)
        {
            TPTNode* oldNode = *it_OldNode;
            newChildrenNode.push_back(nodeLinkerFromPrevToCurr.at(oldNode));
        }

        //2. parent node
        TPTNode* newParentNode = NULL;

        if (newNode->get_parent_node() != NULL)
        {
            newParentNode = nodeLinkerFromPrevToCurr.at(newNode->get_parent_node());
        }

        //3. set nodes
        newNode->set_children_nodes(newChildrenNode);
        newNode->set_parent_node(newParentNode);
    }


}


void TopologicalPathTree::modify_root_node(TPTNode* oldRootNode, const unordered_map<TPTNode*, TPTNode*>& nodeLinkerFromPrevToCurr)
{
    m_RootNode = nodeLinkerFromPrevToCurr.at(oldRootNode);
}


void TopologicalPathTree::sort_leaf_node_by_its_depth()
{
    m_SortedLeafNodes = m_LeafNodes;

    m_SortedLeafNodes.sort(compare_two_leaf_nodes_in_non_decreasing_order);
}



bool TopologicalPathTree::does_this_node_have_ancestor_having_this_face(TPTNode* node, FaceBU2D* face) const
{
    bool doesThisNodeHaveAncestorHavingThisFace = false;

    TPTNode* currNode = node->get_parent_node();

    while (currNode != NULL)
    {
        if (currNode->get_QT_face() == face)
        {
            doesThisNodeHaveAncestorHavingThisFace = true;
            break;
        }

        currNode = currNode->get_parent_node();
    }

    return doesThisNodeHaveAncestorHavingThisFace;
}


void TopologicalPathTree::get_topological_path_from_this_leaf_node(TPTNode* leafNode, list<FaceBU2D*>& topologicalQTFacePath) const
{
    TPTNode* currNode = leafNode;

    while (currNode != NULL)
    {
        FaceBU2D* currFace = currNode->get_QT_face();
        topologicalQTFacePath.push_front(currFace);

        currNode = currNode->get_parent_node();
    }
}


void TopologicalPathTree::get_topological_paths_in_non_decreasing_order_by(const int & numOfPaths, vector<list<FaceBU2D*>>& topologicalQTFacePaths)
{
    if (m_SortedLeafNodes.empty())
    {
        sort_leaf_node_by_its_depth();
    }

    int count = 0;

    for (list<TPTNode*>::const_iterator it_Node = m_SortedLeafNodes.begin();
         it_Node != m_SortedLeafNodes.end();
         it_Node++)
    {
        if (count < numOfPaths)
        {
            TPTNode* currLeafNode = *it_Node;

            list<FaceBU2D*> topologicalQTFacePath;
            get_topological_path_from_this_leaf_node(currLeafNode, topologicalQTFacePath);
            topologicalQTFacePaths.push_back(topologicalQTFacePath);

            count++;
        }
    }
}


bool TopologicalPathTree::compare_two_leaf_nodes_in_non_decreasing_order(TPTNode * node1, TPTNode * node2)
{
    if (node1->get_depth() < node2->get_depth())
    {
        return true;
    }
    else
    {
        return false;
    }
}
