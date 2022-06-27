#include "TPTNode.h"

TPTNode::TPTNode()
{
    m_ParentNode = NULL;
    m_QTFace     = NULL;
    m_Idx        = -1;
    m_Depth      = 0;
}

TPTNode::TPTNode(FaceBU2D* face)
{
    m_ParentNode = NULL;
    m_QTFace     = face;
    m_Idx        = -1;
    m_Depth = 0;
}

TPTNode::TPTNode(FaceBU2D* face, TPTNode* parentNode)
{
    m_ParentNode = parentNode;
    m_QTFace     = face;
    m_Idx        = -1;

    if (parentNode != NULL)
    {
        m_Depth = parentNode->m_Depth + 1;
    }
    else
    {
        m_Depth = 0;
    }
}

TPTNode::TPTNode(FaceBU2D* face, TPTNode* parentNode, const list<TPTNode*>& childrenNodes)
{
    m_ParentNode = parentNode;
    m_QTFace     = face;
    m_ChildrenNodes = childrenNodes;
    m_Idx        = -1;

    if (parentNode != NULL)
    {
        m_Depth = parentNode->m_Depth + 1;
    }
    else
    {
        m_Depth = 0;
    }
}


TPTNode::TPTNode(const TPTNode& node, const int& idx)
{
    copy_from(node);
    m_Idx = idx;
}

TPTNode::TPTNode(const TPTNode& node)
{
    copy_from(node);
}


TPTNode::~TPTNode()
{
}

TPTNode& TPTNode::operator=(const TPTNode& node)
{
    if (this == &node)
    {
        return *this;
    }

    clear();
    copy_from(node);

    return *this;
}

void TPTNode::clear()
{
    m_ParentNode = NULL;
    m_QTFace     = NULL;

    m_ChildrenNodes.clear();
}

void TPTNode::copy_from(const TPTNode& node)
{
    m_ChildrenNodes.clear();

    m_ParentNode    = node.m_ParentNode;
    m_QTFace        = node.m_QTFace;
    m_ChildrenNodes = node.m_ChildrenNodes;
    m_Idx           = node.m_Idx;
    m_Depth         = node.m_Depth;
}

