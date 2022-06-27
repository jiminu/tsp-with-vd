#ifndef _PRIORITY_Q_VEDGE_H
#define _PRIORITY_Q_VEDGE_H

#include "PriorityQueue.h"
#include "VEdge2D.h"
#include "rg_Circle2D.h"

#include <map>
using namespace std;

namespace V {
namespace GeometryTier {


typedef pair<VEdge2D*, rg_Circle2D> edgeNCirclePair;
typedef PQ_Node<edgeNCirclePair>    PQ_Node_VEdge;

class Priority_Q_VEdge : public PriorityQueue<edgeNCirclePair>
{
private:
    map<VEdge2D*, PQ_Node_VEdge*>       m_mapBetweenEdgeAndNode;

public:
    Priority_Q_VEdge();
    Priority_Q_VEdge( const Priority_Q_VEdge& pq_vedge );
    ~Priority_Q_VEdge();

    PQ_Node_VEdge*      getNode( VEdge2D* edge );
    void                changeCircumCircle( VEdge2D* edge, const rg_Circle2D& circle );
    void                kill( VEdge2D* edge );

private:
    int                 getIndex( VEdge2D* edge );
};



Priority_Q_VEdge::Priority_Q_VEdge()
{
}



Priority_Q_VEdge::Priority_Q_VEdge( const Priority_Q_VEdge& pq_vedge )
{
    copyFrom(pq_vedge);

    m_mapBetweenEdgeAndNode.clear();
    map<VEdge2D*, PQ_Node_VEdge*>::const_iterator i_map;
    for( i_map = pq_vedge.m_mapBetweenEdgeAndNode.begin() ; i_map != pq_vedge.m_mapBetweenEdgeAndNode.end() ; ++i_map ) {
        m_mapBetweenEdgeAndNode[(*i_map).first] = (*i_map).second;
    }
}



Priority_Q_VEdge::~Priority_Q_VEdge()
{
    clear();
    m_mapBetweenEdgeAndNode.clear();
}



inline PQ_Node_VEdge* Priority_Q_VEdge::getNode( VEdge2D* edge )
{
    if( m_mapBetweenEdgeAndNode.find(edge) != m_mapBetweenEdgeAndNode.end() )
        return m_mapBetweenEdgeAndNode[edge];
    else
        return NULL;
}



inline void Priority_Q_VEdge::changeCircumCircle( VEdge2D* edge, const rg_Circle2D& circle )
{
    PQ_Node_VEdge* node = getNode( edge );

    if( node != NULL )
    {
        int idx = node->getIndex();
        changeKeyValue( idx, circle.getRadius() );
    }
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
    m_mapBetweenEdgeAndNode.erase( edge );
}


} // GeometryTier
} // V


#endif


