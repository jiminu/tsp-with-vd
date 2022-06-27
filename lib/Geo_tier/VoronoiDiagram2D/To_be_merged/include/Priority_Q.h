#ifndef _PRIORITY_Q_H
#define _PRIORITY_Q_H

#include <vector>
using namespace std;

namespace BULL2D{
namespace GeometryTier {


template <class T>
class Priority_Q;

template <class T>
class PQ_Node
{
    friend class Priority_Q<T>;

private:
    T           m_entity;
    double      m_key;
    int         m_index;

public:
    PQ_Node();
    PQ_Node( const T& entity, const double& key, const int& idx );
    PQ_Node( const PQ_Node<T>& node );
    ~PQ_Node();

    T           getEntity() const;
    T*          getpEntity();
    double      getKey() const;
    int         getIndex() const;

    void        setEntity( const T& entity);
    void        setKey( const double& key);
    void        setIndex( const int& idx);

    PQ_Node<T>& operator=( const PQ_Node<T>& pqNode );
    bool        operator==( const PQ_Node<T>& pqNode ) const;
    bool        operator!=( const PQ_Node<T>& pqNode ) const;
    bool        operator<( const PQ_Node<T>& pqNode ) const;
    bool        operator>( const PQ_Node<T>& pqNode ) const;
};




template<class T>
class Priority_Q
{
private:
    vector< PQ_Node<T>* >   m_container;

public:
    Priority_Q();
    Priority_Q( const T& entity, const double& key );
    Priority_Q( const Priority_Q& pq );
    ~Priority_Q();

    T*                      top();
    T                       pop();
    PQ_Node<T>*             push( const T& entity, const double& key );
    PQ_Node<T>*             push( const pair<T, double>& aPair );
    PQ_Node<T>*             topNode();
    PQ_Node<T>*             getNodeAt( const int& idx );

    int                     size() const;
    bool                    empty() const;

    void                    changeKeyValue( const int& idx, const double& key );
    T                       killNodeAt( const int& idx );
    void                    clear();

    Priority_Q<T>&          operator=( const Priority_Q<T>& p_q );

private:
    int                     computeIndexOfParent( const int& idx );
    int                     computeIndexOfLeftChild( const int& idx );
    int                     computeIndexOfRightChild( const int& idx );

    void                    bubbleUp( const int& idx );
    void                    bubbleDown( const int& idx );
    void                    swapNode( const int& idx1, const int& idx2 );
    int                     getNumChildNode( const int& idx );
};



////////////////////////////////////////////////
// Priority queue class using template define //
////////////////////////////////////////////////

// Node for Priority queue class using template

template<class T>
PQ_Node<T>::PQ_Node()
:m_entity(), m_key(-1.0), m_index(-1)
{
    
}



template <class T>
PQ_Node<T>::PQ_Node( const T& entity, const double& key, const int& idx )
:m_entity(entity), m_key(key), m_index(idx)
{

}



template <class T>
PQ_Node<T>::PQ_Node( const PQ_Node<T>& node )
:m_entity(node.m_entity), m_key(node.m_key), m_index(node.m_index)
{
    
}



template <class T>
PQ_Node<T>::~PQ_Node()
{

}



template <class T>
inline T PQ_Node<T>::getEntity() const
{
    return m_entity;
}



template <class T>
inline T* PQ_Node<T>::getpEntity()
{
    return &m_entity;
}



template <class T>
inline double PQ_Node<T>::getKey() const
{
    return m_key;
}



template <class T>
inline int PQ_Node<T>::getIndex() const
{
    return m_index;
}



template <class T>
inline void PQ_Node<T>::setEntity( const T& entity )
{
    m_entity = entity;
}



template <class T>
inline void PQ_Node<T>::setKey( const double& key )
{
    m_key = key;
}



template <class T>
inline void PQ_Node<T>::setIndex( const int& idx )
{
    m_index = idx;
}



template <class T>
inline PQ_Node<T>& PQ_Node<T>::operator=( const PQ_Node<T>& pqNode )
{
    if( this == *pqNode )
    {
        return *this;
    }

    m_entity    = pqNode.m_entity;
    m_key       = pqNode.m_key;
    m_index     = pqNode.m_index;

    return *this;
}



template <class T>
inline bool PQ_Node<T>::operator==( const PQ_Node<T>& pqNode ) const
{
    if( this == &pqNode )
    {
        return true;
    }

    if( m_entity == pqNode.m_entity &&
        m_key == pqNode.m_key &&
        m_index == pqNode.m_index )
    {
        return true;
    }
    else
    {
        return false;
    }
}



template <class T>
inline bool PQ_Node<T>::operator!=( const PQ_Node<T>& pqNode ) const
{
    if( *this == pqNode )
    {
        return false;
    }
    else
    {
        return true;
    }
}



template <class T>
inline bool PQ_Node<T>::operator<( const PQ_Node<T>& pqNode ) const
{
    return ( m_key < pqNode.m_key );
//     if( m_key < pqNode.m_key )
//     {
//         return true;
//     }
//     else
//     {
//         return false;
//     }
}



template <class T>
inline bool PQ_Node<T>::operator>( const PQ_Node<T>& pqNode ) const
{
    return ( m_key > pqNode.m_key );
}



///////////////////////////////////////////////////////////
// Implementation of Priority queue class using template //
///////////////////////////////////////////////////////////



template<class T>
Priority_Q<T>::Priority_Q()
{

}



template<class T>
Priority_Q<T>::Priority_Q( const T& entity, const double& key )
{
    PQ_Node<T>* pqNode = new PQ_Node<T>(entity, key, 0);
    m_container.push_back( pqNode );
}



template<class T>
Priority_Q<T>::Priority_Q( const Priority_Q& pq )
{
    clear();
    int Q_size = pq.size();

    for( int i = 0 ; i < Q_size ; i++ ) {
        T&      entity  = pq.m_container[i]->getEntity();
        double  key     = pq.m_container[i]->getKey();
        int     idx     = pq.m_container[i]->getIndex();

        PQ_Node<T>* node = new PQ_Node<T>( entity, key, idx );

        m_container.pusdh_back( node );
    }
}



template<class T>
Priority_Q<T>::~Priority_Q()
{
    clear();
}



template<class T>
inline T* Priority_Q<T>::top()
{
    ////////
    if( m_container.size() == 0 )
    {
        return NULL;
    }
    ////////


    return m_container.front()->getpEntity();
}



template<class T>
inline T Priority_Q<T>::pop()
{
    int bottomIndex = m_container.size() - 1;

    T topEntity = m_container[0]->m_entity;
    swapNode( 0, bottomIndex );

    delete m_container[bottomIndex];
    m_container.pop_back();

    bubbleDown( 0 );

    return topEntity;
}



template<class T>
inline PQ_Node<T>* Priority_Q<T>::push( const T& entity, const double& key )
{
    int idx = m_container.size();
    PQ_Node<T>* node = new PQ_Node<T>( entity, key, idx );

    //////
    m_container.push_back( node );
    //////

    bubbleUp( idx );

    return node;
}



template<class T>
inline PQ_Node<T>* Priority_Q<T>::push( const pair<T, double>& aPair )
{
    int idx = m_container.size();
    PQ_Node<T>* node = new PQ_Node<T>( aPair.first, aPair.second, idx );

    //////
    m_container.push_back( node );
    //////

    bubbleUp( idx );

    return node;
}



template<class T>
inline PQ_Node<T>* Priority_Q<T>::topNode()
{
    //////
    if( m_container.size() == 0 )
    {
        return NULL;
    }
    //////


    return m_container[0];
}



template<class T>
inline PQ_Node<T>* Priority_Q<T>::getNodeAt( const int& idx )
{
    return m_container[idx];
}



template<class T>
inline int Priority_Q<T>::size() const
{
    return m_container.size();
}



template<class T>
inline bool Priority_Q<T>::empty() const
{
    return m_container.empty();
}



template<class T>
inline void Priority_Q<T>::changeKeyValue( const int& idx, const double& key )
{
    double currKey = m_container[idx]->m_key;
    m_container[idx]->m_key = key;

    if( key < currKey )
    {
        bubbleUp( idx );
    }
    else
    {
        bubbleDown( idx );
    }
}



template<class T>
inline T Priority_Q<T>::killNodeAt( const int& idx )
{
    if( idx < 0 || idx >= m_container.size() )
    {
        exit(1);
    }

    int bottomIndex = m_container.size() - 1;

    T entity = m_container[idx]->getEntity();
    swapNode( idx, bottomIndex );

    delete m_container[bottomIndex];
    m_container.pop_back();

    bubbleDown( idx );

    return entity;
}



template<class T>
inline void Priority_Q<T>::clear()
{
    int Q_size = m_container.size();

    for( int i = 0 ; i < Q_size ; i++ ) {
        delete m_container[i];
    }

    m_container.clear();
}



template<class T>
Priority_Q<T>& Priority_Q<T>::operator=( const Priority_Q<T>& p_q )
{
    if( this == &p_q )
    {
        return *this;
    }

    clear();
    for( int i = 0 ; i < p_q.size() ; i++ ) {
        T       entity = p_q.m_container[i]->m_entity;
        double  key    = p_q.m_container[i]->m_key;
        int     idx    = p_q.m_container[i]->m_index;

        PQ_Node<T>* node = new PQ_Node<T>(entity, key, idx);
        m_container.push_back( node );
    }

    return *this;
}



template<class T>
inline int Priority_Q<T>::computeIndexOfParent( const int& idx )
{
    return (int)( (idx-1) / 2 );
}



template<class T>
inline int Priority_Q<T>::computeIndexOfLeftChild( const int& idx )
{
    return idx * 2 + 1;
}



template<class T>
inline int Priority_Q<T>::computeIndexOfRightChild( const int& idx )
{
    return idx * 2 + 2;
}



template<class T>
void Priority_Q<T>::bubbleUp( const int& idx )
{
    if( idx < 0 || idx > (m_container.size() - 1) )
    {
        exit(1);
    }

    int currIndex       = idx;
    int indexOfParent   = 0;
    while( currIndex > 0 ) // Index of top is 0
    {
        indexOfParent = computeIndexOfParent( currIndex );

        if( m_container[currIndex]->m_key < m_container[indexOfParent]->m_key )
        {
            swapNode( currIndex, indexOfParent );
            currIndex = indexOfParent;
        }
        else
        {
            break;
        }
    }
}



template<class T>
void Priority_Q<T>::bubbleDown( const int& idx )
{
    if( idx < 0 || idx > (m_container.size() - 1) )
    {
        exit(1);
    }

    int     currIndex         = idx;
    int     indexOfLeftChild  = 0;
    int     indexOfRightChild = 0;

    int     queueSize = m_container.size();
    bool    canBeChanged = true;
    while( canBeChanged )
    {
        int     numChild = getNumChildNode( currIndex );
        int     indexOfSmallerChild = 0;

        switch ( numChild )
        {

        case 0:
            canBeChanged = false;
        	break;


        case 1:
            indexOfLeftChild = computeIndexOfLeftChild( currIndex );
            if( m_container[indexOfLeftChild]->m_key < m_container[currIndex]->m_key )
            {
                swapNode( currIndex, indexOfLeftChild );
                currIndex = indexOfLeftChild;
            }
            canBeChanged = false;
            break;


        case 2:
            indexOfLeftChild  = computeIndexOfLeftChild( currIndex );
            indexOfRightChild = indexOfLeftChild + 1;
            
            if( m_container[indexOfLeftChild]->m_key < m_container[indexOfRightChild]->m_key )
            {
                indexOfSmallerChild = indexOfLeftChild;
            }
            else
            {
                indexOfSmallerChild = indexOfRightChild;
            }

            if( m_container[indexOfSmallerChild]->m_key < m_container[currIndex]->m_key )
            {
                swapNode( currIndex, indexOfSmallerChild );
                currIndex = indexOfSmallerChild;
            }
            break;

        default:
            exit(1);
            break;
        }
    }
}



template<class T>
inline void Priority_Q<T>::swapNode( const int& idx1, const int& idx2 )
{
    PQ_Node<T>* tempNode = m_container[idx1];
    m_container[idx1] = m_container[idx2];
    m_container[idx2] = tempNode;

    m_container[idx1]->m_index = idx1;
    m_container[idx2]->m_index = idx2;
}



template<class T>
inline int Priority_Q<T>::getNumChildNode( const int& idx )
{
    if( idx < 0 || idx > (m_container.size() - 1) )
    {
        return 0;
    }

    int numChild = 0;
    int indexOfRightChild = computeIndexOfRightChild( idx );

    if( indexOfRightChild < m_container.size() )
    {
        numChild = 2;
    }
    else if( indexOfRightChild == m_container.size() )
    {
        numChild = 1;
    }
    else
    {
        numChild = 0;
    }

    return numChild;
}


} // namespace GeometryTier
} // namespace BULL2D

#endif

