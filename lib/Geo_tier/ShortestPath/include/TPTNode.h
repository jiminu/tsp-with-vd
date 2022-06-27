#ifndef _TPT_NODE_
#define _TPT_NODE_

#include "FaceBU2D.h"
#include <list>
using namespace std;

class TPTNode
{
private:
    TPTNode*         m_ParentNode;
    list<TPTNode*>   m_ChildrenNodes;

    FaceBU2D*        m_QTFace;

    int              m_Idx;
    int              m_Depth;

public:
    //constructor
    TPTNode();
    TPTNode(FaceBU2D* face);
    TPTNode(FaceBU2D* face, TPTNode* parentNode);
    TPTNode(FaceBU2D* face, TPTNode* parentNode, const list<TPTNode*>& childrenNodes);
    TPTNode(const TPTNode& node, const int& idx);
    TPTNode(const TPTNode& node);

    //destructor;
    ~TPTNode();

    //equal operation
    TPTNode& operator=(const TPTNode& node);

    //setter
    inline void set_parent_node(TPTNode* parentNode) { m_ParentNode = parentNode; };
    inline void set_QT_face(FaceBU2D* qtFace)        { m_QTFace = qtFace; };
    inline void set_children_nodes(const list<TPTNode*> childrenNode) { m_ChildrenNodes.clear();  m_ChildrenNodes = childrenNode; };
    inline void set_depth(const int& depth) { m_Depth = depth; };
    //getter
    inline TPTNode*  get_parent_node() const { return m_ParentNode; };
    inline FaceBU2D* get_QT_face() const { return m_QTFace; };
    inline void      get_children_nodes(list<TPTNode*>& childrenNodes) const { childrenNodes = m_ChildrenNodes; };
    inline int       get_depth() const { return m_Depth; };

    inline void      add_a_child(TPTNode* childNode) { m_ChildrenNodes.push_back(childNode); };
    void clear();

private:
    void copy_from(const TPTNode& node);

};


#endif