#ifndef _PRIORITY_Q_VEDGE_H
#define _PRIORITY_Q_VEDGE_H

#include "Priority_Q.h"
#include "VEdge2D.h"
#include "rg_Circle2D.h"

#include <map>
using namespace std;

namespace BULL2D {
namespace GeometryTier {


typedef pair<VEdge2D*, rg_Circle2D> edgeNCirclePair;
typedef PQ_Node<edgeNCirclePair>    PQ_Node_VEdge;

class Priority_Q_VEdge
{
private:
    Priority_Q< edgeNCirclePair >   m_priorityQ;
    map<VEdge2D*, PQ_Node_VEdge*>   m_mapBetweenEdgeAndNode;

public:
    Priority_Q_VEdge();
    Priority_Q_VEdge( const edgeNCirclePair& aPair );
    Priority_Q_VEdge( const Priority_Q_VEdge& pq_vedge );
    ~Priority_Q_VEdge();

    edgeNCirclePair*    top();
    edgeNCirclePair     pop();
    edgeNCirclePair*    push( const edgeNCirclePair& aPair );
    PQ_Node_VEdge*      getNode( VEdge2D* edge );

    int                 size() const;
    bool                empty() const;

    void                changeCircumCircle( VEdge2D* edge, const rg_Circle2D& circle );
    void                kill( VEdge2D* edge );
    void                clear();

private:
    PQ_Node_VEdge*      topNode();
    int                 getIndex( VEdge2D* edge );
    edgeNCirclePair     killNodeAt( const int& idx );
};



Priority_Q_VEdge::Priority_Q_VEdge()
{

}



Priority_Q_VEdge::Priority_Q_VEdge( const edgeNCirclePair& aPair )
{
    PQ_Node_VEdge* node = m_priorityQ.push( aPair, aPair.second.getRadius() );
    m_mapBetweenEdgeAndNode.insert( pair<VEdge2D*, PQ_Node_VEdge*>(aPair.first, node) );
}



Priority_Q_VEdge::Priority_Q_VEdge( const Priority_Q_VEdge& pq_vedge )
{
    m_priorityQ = pq_vedge.m_priorityQ;

    m_mapBetweenEdgeAndNode.clear();
    map<VEdge2D*, PQ_Node_VEdge*>::const_iterator i_map;
    for( i_map = pq_vedge.m_mapBetweenEdgeAndNode.begin() ; i_map != pq_vedge.m_mapBetweenEdgeAndNode.end() ; i_map++ ) {
        m_mapBetweenEdgeAndNode[(*i_map).first] = (*i_map).second;
    }

}



Priority_Q_VEdge::~Priority_Q_VEdge()
{
    m_mapBetweenEdgeAndNode.clear();
}



inline edgeNCirclePair* Priority_Q_VEdge::top()
{
    return m_priorityQ.top();
}



inline edgeNCirclePair Priority_Q_VEdge::pop()
{
    edgeNCirclePair aPair = m_priorityQ.pop();
    m_mapBetweenEdgeAndNode.erase(aPair.first);

    return aPair;
}



inline edgeNCirclePair* Priority_Q_VEdge::push( const edgeNCirclePair& aPair )
{
    PQ_Node_VEdge* node = m_priorityQ.push( aPair, aPair.second.getRadius() );
    m_mapBetweenEdgeAndNode[aPair.first] = node;

    return node->getpEntity();
}



inline PQ_Node_VEdge* Priority_Q_VEdge::topNode()
{
    return m_priorityQ.topNode();
}



inline PQ_Node_VEdge* Priority_Q_VEdge::getNode( VEdge2D* edge )
{
    if( m_mapBetweenEdgeAndNode.find(edge) != m_mapBetweenEdgeAndNode.end() )
    {
        return m_mapBetweenEdgeAndNode[edge];
    }
    else
    {
        return NULL;
    }
}



inline int Priority_Q_VEdge::size() const
{
    return m_priorityQ.size();
}



inline bool Priority_Q_VEdge::empty() const
{
    return m_priorityQ.empty();
}



inline void Priority_Q_VEdge::changeCircumCircle( VEdge2D* edge, const rg_Circle2D& circle )
{
    PQ_Node_VEdge* node = getNode( edge );

    if( node != NULL )
    {
        node->getpEntity()->second = circle;
        m_priorityQ.changeKeyValue( node->getKey(), circle.getRadius() );
    }
}



inline edgeNCirclePair Priority_Q_VEdge::killNodeAt( const int& idx )
{
    edgeNCirclePair aPair = m_priorityQ.killNodeAt( idx );
    m_mapBetweenEdgeAndNode.erase( aPair.first );

    return aPair;
}



inline void Priority_Q_VEdge::clear()
{
    m_priorityQ.clear();

    map<VEdge2D*, PQ_Node_VEdge*>::iterator i_map = m_mapBetweenEdgeAndNode.begin();
    for( ; i_map != m_mapBetweenEdgeAndNode.end() ; i_map++ ) {
        (*i_map).first->setFalseCandidateForFlippingInPhantomRemoval();
    }

    m_mapBetweenEdgeAndNode.clear();

}



inline int Priority_Q_VEdge::getIndex( VEdge2D* edge )
{
    if( m_mapBetweenEdgeAndNode.find( edge ) != m_mapBetweenEdgeAndNode.end() )
    {
        return m_mapBetweenEdgeAndNode[edge]->getIndex();
    }
    else
    {
        return -1;
    }
}



inline void Priority_Q_VEdge::kill( VEdge2D* edge )
{
    int idx = getIndex( edge );
    killNodeAt( idx );
}

} // namespace GeometryTier
} // namespace BULL2D

#endif


