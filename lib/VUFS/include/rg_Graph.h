#ifndef _RG_GRAPH_H_
#define _RG_GRAPH_H_

#include "rg_Const.h"
#include "rg_dList.h"
#include "TopologicalEntity.h"

template<class GNodeElement, class GArcElement>
class rg_GArc;

template<class GNodeElement, class GArcElement>
class rg_GNode : public TopologicalEntity
{
private:
    rg_dList< rg_GArc<GNodeElement, GArcElement>* > m_incidentArc;
    GNodeElement                      m_element;

public:
    rg_GNode();
    rg_GNode(const rg_INT& ID);
    rg_GNode(const rg_INT& ID, const GNodeElement& element);
    rg_GNode(const rg_GNode<GNodeElement, GArcElement>& node);
    ~rg_GNode();

    inline rg_INT       getDegree() const { return m_incidentArc.getSize(); }
    inline GNodeElement getElement() const { return m_element; }
    inline const rg_dList< rg_GArc<GNodeElement, GArcElement>* >& getIncidentArcs() const { return m_incidentArc; }

    inline void   setElement(const GNodeElement& element) { m_element = element; }
    inline void   addIncidentArc(rg_GArc<GNodeElement, GArcElement>* arc) { m_incidentArc.add( arc ); }

    rg_GNode<GNodeElement, GArcElement>& operator =(const rg_GNode<GNodeElement, GArcElement>& node);

    rg_BOOL isAdjacent( const rg_GNode<GNodeElement, GArcElement>* node ) const;
    rg_INT  searchAdjacentNodes( rg_dList< rg_GNode<GNodeElement, GArcElement>* >& adjacentNodes ) const;
};


template<class GNodeElement, class GArcElement>
class rg_GArc : public TopologicalEntity
{
private:
    rg_GNode<GNodeElement, GArcElement>* m_origin;
    rg_GNode<GNodeElement, GArcElement>* m_destination;
    GArcElement             m_element;

public:
    rg_GArc();
    rg_GArc(const rg_INT& ID);
    rg_GArc(const rg_INT& ID, const GArcElement& element);
    rg_GArc(const rg_INT& ID, rg_GNode<GNodeElement, GArcElement>* origin, rg_GNode<GNodeElement, GArcElement>* destination, const GArcElement& element);
    rg_GArc(const rg_GArc<GNodeElement, GArcElement>& arc);
    ~rg_GArc();

    inline rg_GNode<GNodeElement, GArcElement>* getOrigin() const { return m_origin; }
    inline rg_GNode<GNodeElement, GArcElement>* getDestination() const { return m_destination; }
    inline GNodeElement            getElement() const { return m_element; }

    inline void setOrigin(rg_GNode<GNodeElement, GArcElement>* origin) { m_origin = origin; }
    inline void setDestination(rg_GNode<GNodeElement, GArcElement>* destination) { m_destination = destination; }
    inline void setElement(const GArcElement& element) { m_element = element; }

    rg_GArc<GNodeElement, GArcElement>& operator =(const rg_GArc<GNodeElement, GArcElement>& arc);

    rg_BOOL isOrigin( rg_GNode<GNodeElement, GArcElement>* node ) const;
    rg_BOOL isDestination( rg_GNode<GNodeElement, GArcElement>* node ) const;
};


template<class GNodeElement, class GArcElement>
class rg_Graph
{
private:
    rg_dList< rg_GNode<GNodeElement, GArcElement> > m_nodes;
    rg_dList< rg_GArc<GNodeElement, GArcElement> >  m_arcs;

public:
    rg_Graph();
    ~rg_Graph();

    inline rg_INT getNumNodes() const { return m_nodes.getSize(); }
    inline rg_INT getNumArcs() const { return m_arcs.getSize(); }

    inline const rg_dList< rg_GNode<GNodeElement, GArcElement> >& getNodes() const { return m_nodes; }
    inline const rg_dList< rg_GArc<GNodeElement, GArcElement> >&   getArcs() const  { return m_arcs; }

    void insertNode( const rg_GNode<GNodeElement, GArcElement>& node );
    void insertArc(  const rg_GArc<GNodeElement, GArcElement>& arc );

};



////////////////////////////////////////////////////////////////////////////////
// 
//  member functions of rg_GNode
//
template<class GNodeElement, class GArcElement>
rg_GNode<GNodeElement, GArcElement>::rg_GNode()
{
}



template<class GNodeElement, class GArcElement>
rg_GNode<GNodeElement, GArcElement>::rg_GNode(const rg_INT& ID)
: TopologicalEntity(ID)
{
}



template<class GNodeElement, class GArcElement>
rg_GNode<GNodeElement, GArcElement>::rg_GNode(const rg_INT& ID, const GNodeElement& element)
: TopologicalEntity(ID),
  m_element(element)
{
}



template<class GNodeElement, class GArcElement>
rg_GNode<GNodeElement, GArcElement>::rg_GNode(const rg_GNode<GNodeElement, GArcElement>& node)
: TopologicalEntity( node.m_ID ),
  m_element(node.m_element)
{
}



template<class GNodeElement, class GArcElement>
rg_GNode<GNodeElement, GArcElement>::~rg_GNode()
{
}



template<class GNodeElement, class GArcElement>
rg_GNode<GNodeElement, GArcElement>& rg_GNode<GNodeElement, GArcElement>::operator =(const rg_GNode<GNodeElement, GArcElement>& node)
{
    if ( this == &node ) {
        return *this;
    }

    TopologicalEntity::operator=(node);
    m_element = node.m_element;
    
    return *this;
}



template<class GNodeElement, class GArcElement>
rg_BOOL rg_GNode<GNodeElement, GArcElement>::isAdjacent( const rg_GNode<GNodeElement, GArcElement>* node ) const
{
    rg_BOOL isAdjacentNode = rg_FALSE;

    m_incidentArc.reset4Loop();
    while ( m_incidentArc.setNext4Loop() ) {
        rg_GArc<GNodeElement, GArcElement>* currArc = m_incidentArc.getEntity();

        if ( node == currArc->getOrigin() || node == currArc->getDestination() ) {
            isAdjacentNode = rg_TRUE;
            break;
        }
    }
    return isAdjacentNode;
}


template<class GNodeElement, class GArcElement>
rg_INT  rg_GNode<GNodeElement, GArcElement>::searchAdjacentNodes( rg_dList< rg_GNode<GNodeElement, GArcElement>* >& adjacentNodes ) const
{
    m_incidentArc.reset4Loop();
    while ( m_incidentArc.setNext4Loop() ) {
        rg_GArc<GNodeElement, GArcElement>* currArc = m_incidentArc.getEntity();

        if ( this == currArc->getOrigin() ) {
            adjacentNodes.add( currArc->getDestination() );
        }
        else {
            adjacentNodes.add( currArc->getOrigin() );
        }
    }

    return adjacentNodes.getSize();
}




////////////////////////////////////////////////////////////////////////////////
// 
//  member functions of rg_GArc
//

template<class GNodeElement, class GArcElement>
rg_GArc<GNodeElement, GArcElement>::rg_GArc()
{
    m_origin      = rg_NULL;
    m_destination = rg_NULL;
}



template<class GNodeElement, class GArcElement>
rg_GArc<GNodeElement, GArcElement>::rg_GArc(const rg_INT& ID)
: TopologicalEntity(ID)
{
    m_origin      = rg_NULL;
    m_destination = rg_NULL;
}



template<class GNodeElement, class GArcElement>
rg_GArc<GNodeElement, GArcElement>::rg_GArc(const rg_INT& ID, const GArcElement& element)
: TopologicalEntity(ID),
  m_element( element )
{
    m_origin      = rg_NULL;
    m_destination = rg_NULL;
}



template<class GNodeElement, class GArcElement>
rg_GArc<GNodeElement, GArcElement>::rg_GArc(const rg_INT& ID, rg_GNode<GNodeElement, GArcElement>* origin, rg_GNode<GNodeElement, GArcElement>* destination, const GArcElement& element)
: TopologicalEntity(ID),
  m_element( element )
{
    if ( origin == rg_NULL ) {
        m_origin = rg_NULL;
    }
    else {
        m_origin = origin;
        m_origin->addIncidentArc( this );
    }

    if ( destination == rg_NULL ) {
        m_destination = rg_NULL;
    }
    else {
        m_destination = destination;
        m_origin->addIncidentArc( this );
    }
}



template<class GNodeElement, class GArcElement>
rg_GArc<GNodeElement, GArcElement>::rg_GArc(const rg_GArc<GNodeElement, GArcElement>& arc)
: TopologicalEntity(arc.m_ID),
  m_element( arc.m_element )
{
    m_origin = arc.m_origin;
    m_destination = arc.m_destination;
}



template<class GNodeElement, class GArcElement>
rg_GArc<GNodeElement, GArcElement>::~rg_GArc()
{
}



template<class GNodeElement, class GArcElement>
rg_GArc<GNodeElement, GArcElement>& rg_GArc<GNodeElement, GArcElement>::operator =(const rg_GArc<GNodeElement, GArcElement>& arc)
{
    if ( this == &node ) {
        return *this;
    }

    TopologicalEntity::operator=(node);
    m_element = node.m_element;

    m_origin = arc.m_origin;
    m_destination = arc.m_destination;

    return *this;
}



template<class GNodeElement, class GArcElement>
rg_BOOL rg_GArc<GNodeElement, GArcElement>::isOrigin( rg_GNode<GNodeElement, GArcElement>* node ) const
{
    if ( m_origin == node ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



template<class GNodeElement, class GArcElement>
rg_BOOL rg_GArc<GNodeElement, GArcElement>::isDestination( rg_GNode<GNodeElement, GArcElement>* node ) const
{
    if ( m_destination == node ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}




////////////////////////////////////////////////////////////////////////////////
// 
//  member functions of rg_Graph
//

template<class GNodeElement, class GArcElement>
rg_Graph<GNodeElement, GArcElement>::rg_Graph()
{
}



template<class GNodeElement, class GArcElement>
rg_Graph<GNodeElement, GArcElement>::~rg_Graph()
{
}





template<class GNodeElement, class GArcElement>
void rg_Graph<GNodeElement, GArcElement>::insertNode( const rg_GNode<GNodeElement, GArcElement>& node )
{
    m_nodes.add( node );
}



template<class GNodeElement, class GArcElement>
void rg_Graph<GNodeElement, GArcElement>::insertArc(  const rg_GArc<GNodeElement, GArcElement>& arc )
{
    m_arcs.add( arc );
}





#endif


