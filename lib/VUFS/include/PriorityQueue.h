#ifndef _PRIORITY_QUEUE_H
#define _PRIORITY_QUEUE_H
#include <cstddef>
#include <functional>
#include <vector>
#include <list>
using namespace std;

#include "rg_RelativeOp.h"

template <typename T>
class PQ_Node
{
	template<typename _T, typename _COMPARATOR> 
    friend class PriorityQueue;

	template<typename _T, typename _COMPARATOR, typename _Hash, typename _EqualPred> 
    friend class EntityAccessiblePriorityQ;

private:
    T           m_entity;
    double      m_key;
    int         m_index;

public:
    PQ_Node();
    PQ_Node( const T& entity, const double& key, const int& idx );
    PQ_Node( const PQ_Node<T>& node );
    ~PQ_Node();

    const T&    getEntity()  const;           //cy song :  T getEntity() -> const T&    getEntity() const; 
    const T*    getpEntity() const;           //cy song :  Don't need it...?
    double      getKey() const;
    int         getIndex() const;

//     void        setEntity( const T& entity);   //cy song :  Don't need it...?
//     void        setKey( const double& key);
//     void        setIndex( const int& idx);

    PQ_Node<T>& operator=( const PQ_Node<T>& pqNode );
    bool        operator==( const PQ_Node<T>& pqNode ) const;
    bool        operator!=( const PQ_Node<T>& pqNode ) const;
    bool        operator<( const PQ_Node<T>& pqNode ) const;
    bool        operator>( const PQ_Node<T>& pqNode ) const;
};



template<typename T, typename COMPARATOR = less<PQ_Node<T>>>
class PriorityQueue
{
protected:
    vector< PQ_Node<T>* >   m_container;
	COMPARATOR              m_comp; //comparator

public:
	typedef typename vector< PQ_Node<T>* >::iterator       iterator;
	typedef typename vector< PQ_Node<T>* >::const_iterator const_iterator;

    PriorityQueue();
    PriorityQueue( const T& entity, const double& key );
    PriorityQueue( const PriorityQueue& pq );
    ~PriorityQueue();

    //const T&                top() const;                                     //cysong T* -> const T& const
    T                       pop();
    const PQ_Node<T>*       push( const T& entity, const double& key );  //cysong PQ_Node<T>*  -> const  PQ_Node<T>* 
    const PQ_Node<T>*       push( const pair<T, double>& aPair );        //cysong PQ_Node<T>*  -> const  PQ_Node<T>* 
    const PQ_Node<T>*       topNode() const;                             //cysong PQ_Node<T>* topNode() -> const PQ_Node<T>* topNode() const
    int                     size() const;
    bool                    empty() const;
    void                    findEntitiesHavingSamekey(const double& key, list<T>& outputEntities, const double& res = resNeg4);
    void                    clear();

	iterator                begin();
	iterator                end();

	const_iterator          begin() const;
	const_iterator          end()   const;

    PriorityQueue&          operator=( const PriorityQueue& pq );


protected:
	const PQ_Node<T>*       getNodeAt(const int& idx) const;           //cysong PQ_Node<T>* getNodeAt( const int& idx ) -> const PQ_Node<T>* getNodeAt( const int& idx ) const
	void                    changeKeyValue(const int& idx, const double& key);
	void                    changeEntity(const int& idx, const T& entity);
	T                       killNodeAt(const int& idx);
	void                    copyFrom(const PriorityQueue& pq);

private:
    int                     computeIndexOfParent( const int& idx );
    int                     computeIndexOfLeftChild( const int& idx );
    int                     computeIndexOfRightChild( const int& idx );

    void                    bubbleUpOrDownOrNot(const int& idx);
	void                    bubbleUp( const int& idx );
    void                    bubbleDown( const int& idx );
    void                    swapNode( const int& idx1, const int& idx2 );
    int                     getNumChildNode( const int& idx );

};




///////////////////////////////////////////////////////////
// Implementation of Priority queue class using template //
///////////////////////////////////////////////////////////



////////////////////////////////////////////////
// Priority queue class using template define //
////////////////////////////////////////////////

// Node for Priority queue class using template

template<typename T>
PQ_Node<T>::PQ_Node()
	:m_entity(), m_key(-1.0), m_index(-1)
{

}



template <typename T>
PQ_Node<T>::PQ_Node(const T& entity, const double& key, const int& idx)
	:m_entity(entity), m_key(key), m_index(idx)
{

}



template <typename T>
PQ_Node<T>::PQ_Node(const PQ_Node<T>& node)
	:m_entity(node.m_entity), m_key(node.m_key), m_index(node.m_index)
{

}



template <typename T>
PQ_Node<T>::~PQ_Node()
{

}



template <typename T>
const T& PQ_Node<T>::getEntity() const
{
	return m_entity;
}



template <typename T>
const T* PQ_Node<T>::getpEntity() const
{
	return &m_entity;
}



template <typename T>
double PQ_Node<T>::getKey() const
{
	return m_key;
}



template <typename T>
int PQ_Node<T>::getIndex() const
{
	return m_index;
}



// template <typename T>
// inline void PQ_Node<T>::setEntity( const T& entity )
// {
//     m_entity = entity;
// }
// 
// 
// 
// template <typename T>
// inline void PQ_Node<T>::setKey( const double& key )
// {
//     m_key = key;
// }
// 
// 
// 
// template <typename T>
// inline void PQ_Node<T>::setIndex( const int& idx )
// {
//     m_index = idx;
// }



template <typename T>
PQ_Node<T>& PQ_Node<T>::operator=(const PQ_Node<T>& pqNode)
{
	if (this == &pqNode)   //cy song : if (this == *pqNode) -> (this == &pqNode)
	{
		return *this;
	}

	m_entity = pqNode.m_entity;
	m_key = pqNode.m_key;
	m_index = pqNode.m_index;

	return *this;
}



template <typename T>
bool PQ_Node<T>::operator==(const PQ_Node<T>& pqNode) const
{
	if (this == &pqNode)
	{
		return true;
	}

	if (m_entity == pqNode.m_entity &&
		m_key == pqNode.m_key &&
		m_index == pqNode.m_index)
	{
		return true;
	}
	else
	{
		return false;
	}
}



template <typename T>
bool PQ_Node<T>::operator!=(const PQ_Node<T>& pqNode) const
{
	if (*this == pqNode)
	{
		return false;
	}
	else
	{
		return true;
	}
}



template <typename T>
bool PQ_Node<T>::operator<(const PQ_Node<T>& pqNode) const
{
	return (m_key < pqNode.m_key);
	//     if( m_key < pqNode.m_key )
	//     {
	//         return true;
	//     }
	//     else
	//     {
	//         return false;
	//     }
}



template <typename T>
bool PQ_Node<T>::operator>(const PQ_Node<T>& pqNode) const
{
	return (m_key > pqNode.m_key);
	//     if( m_key > pqNode.m_key )
	//     {
	//         return true;
	//     }
	//     else
	//     {
	//         return false;
	//     }
}







template<typename T, typename COMPARATOR>
PriorityQueue<T, COMPARATOR>::PriorityQueue()
{
}



template<typename T, typename COMPARATOR>
PriorityQueue<T, COMPARATOR>::PriorityQueue( const T& entity, const double& key )
{
    PQ_Node<T>* pqNode = new PQ_Node<T>(entity, key, 0);
    m_container.push_back( pqNode );
}



template<typename T, typename COMPARATOR>
PriorityQueue<T, COMPARATOR>::PriorityQueue( const PriorityQueue& pq )
{
	copyFrom(pq);
}



template<typename T, typename COMPARATOR>
PriorityQueue<T, COMPARATOR>::~PriorityQueue()
{
    clear();
}



template<typename T, typename COMPARATOR>
T PriorityQueue<T, COMPARATOR>::pop()
{
	if (m_container.size() == 0){
		return T();
	}

    int bottomIndex = m_container.size() - 1;

    T topEntity = m_container[0]->m_entity;
    swapNode( 0, bottomIndex );

    delete m_container[bottomIndex];
    m_container.pop_back();

    bubbleDown( 0 );

    return topEntity;
}



template<typename T, typename COMPARATOR>
const PQ_Node<T>* PriorityQueue<T, COMPARATOR>::push( const T& entity, const double& key )
{
    int idx = m_container.size();
    PQ_Node<T>* node = new PQ_Node<T>( entity, key, idx );
	m_container.push_back(node); //cy song

    bubbleUp( idx );

    return node;
}



template<typename T, typename COMPARATOR>
const PQ_Node<T>* PriorityQueue<T, COMPARATOR>::push( const pair<T, double>& aPair )
{
    int idx = m_container.size();
    PQ_Node<T>* node = new PQ_Node<T>( aPair.first, aPair.second, idx );
	m_container.push_back(node); //cy song

    bubbleUp( idx );

    return node;
}



template<typename T, typename COMPARATOR>
const PQ_Node<T>* PriorityQueue<T, COMPARATOR>::topNode() const 
{
	if (m_container.size() == 0)
    {
		return NULL;
	}

    return m_container[0];
}



template<typename T, typename COMPARATOR>
const PQ_Node<T>* PriorityQueue<T, COMPARATOR>::getNodeAt( const int& idx ) const 
{
	if (idx < 0 || idx >= m_container.size())
	{
		//exit(1);
        return NULL; //20170413. chanyoung song
	}

    return m_container[idx];
}



template<typename T, typename COMPARATOR>
int PriorityQueue<T, COMPARATOR>::size() const
{
    return m_container.size();
}



template<typename T, typename COMPARATOR>
bool PriorityQueue<T, COMPARATOR>::empty() const
{
    return m_container.empty();
}



template<typename T, typename COMPARATOR>
void PriorityQueue<T, COMPARATOR>::changeKeyValue( const int& idx, const double& key )
{
	if (idx < 0 || idx >= m_container.size()){
		//exit(1);
        return; //20170413. chanyoung song
	}

	double currKey = m_container[idx]->m_key;
    m_container[idx]->m_key = key;

	//cy song : do bubble up or down 
    if (key != currKey)
	{
        bubbleUpOrDownOrNot(idx);
	}
}


template<typename T, typename COMPARATOR>
void PriorityQueue<T, COMPARATOR>::changeEntity(const int& idx, const T& entity)
{
	if (idx < 0 || idx >= m_container.size()){
		//exit(1);
        return; //20170413. chanyoung song
	}

	m_container[idx]->m_entity = entity;
}



template<typename T, typename COMPARATOR>
T PriorityQueue<T, COMPARATOR>::killNodeAt( const int& idx )
{
    if( idx < 0 || idx >= m_container.size() )
    {
        //exit(1);
        return T(); //20170413. chanyoung song
    }

    T entity = m_container[idx]->getEntity();
    
    int bottomIndex = m_container.back()->getIndex();

    if (idx != bottomIndex)
    {
        swapNode(idx, bottomIndex);

        delete m_container[bottomIndex];
        m_container.pop_back();

        bubbleUpOrDownOrNot(idx);
         // 올라갈수 있으면 올라감
         // 내려갈수 있으면 내려감
    }
    else
    {
        delete m_container[bottomIndex];
        m_container.pop_back();
    }



    // 선택해야하는 인덱스를 찾기 : output_idx
    // 현재 노드의 자식 개수 여부 확인
    //  1) 자식이 2개 있으면 그 중 m_comp로 왼쪽 선택해서 output_idx로
    //  2) 자식이 1개 있으면 그것을 output_idx로
    //  3) 자식이 없으면 output_idx = -1

    // output_idx가 -1이면 break;
    // output_idx가 -1이 아니면, 현재 노드와 비교해서 m_comp에 따라 내려갈지 결정


    // 그대로라면 현재 노드의 부모의 존재 여부 확인
    // 1) 없으면 그대로 ouput_idx = -1로 하고 break; 
    // 2) 있으면 현재 노드의 값과 비교서 올라갈지를 결정 
   
    return entity;
}



template<typename T, typename COMPARATOR>
void PriorityQueue<T, COMPARATOR>::clear()
{
    int Q_size = m_container.size();

    for( int i = 0 ; i < Q_size ; i++ ) {
        delete m_container[i];
    }

    m_container.clear();
}



template<typename T, typename COMPARATOR>
typename PriorityQueue<T, COMPARATOR>::const_iterator PriorityQueue<T, COMPARATOR>::begin() const
{
	return m_container.begin();
}


template<typename T, typename COMPARATOR>
typename PriorityQueue<T, COMPARATOR>::const_iterator PriorityQueue<T, COMPARATOR>::end() const
{
	return m_container.end();
}


template<typename T, typename COMPARATOR>
typename PriorityQueue<T, COMPARATOR>::iterator PriorityQueue<T, COMPARATOR>::begin()
{
	return m_container.begin();
}


template<typename T, typename COMPARATOR>
typename PriorityQueue<T, COMPARATOR>::iterator PriorityQueue<T, COMPARATOR>::end()
{
	return m_container.end();
}


template<typename T, typename COMPARATOR>
PriorityQueue<T, COMPARATOR>& PriorityQueue<T, COMPARATOR>::operator=(const PriorityQueue<T, COMPARATOR>& pq)
{
    if( this == &pq ){
        return *this;
    }

    clear();
	copyFrom(pq);

    return *this;
}



template<typename T, typename COMPARATOR>
int PriorityQueue<T, COMPARATOR>::computeIndexOfParent( const int& idx )
{
    return (int)( (idx-1) / 2 );
}



template<typename T, typename COMPARATOR>
int PriorityQueue<T, COMPARATOR>::computeIndexOfLeftChild( const int& idx )
{
    return idx * 2 + 1;
}



template<typename T, typename COMPARATOR>
int PriorityQueue<T, COMPARATOR>::computeIndexOfRightChild( const int& idx )
{
    return idx * 2 + 2;
}


template<typename T, typename COMPARATOR /*= less<PQ_Node<T>>*/>
void PriorityQueue<T, COMPARATOR>::bubbleUpOrDownOrNot(const int& idx)
{
    bubbleUp(idx);
    bubbleDown(idx);
}


template<typename T, typename COMPARATOR>
void PriorityQueue<T, COMPARATOR>::bubbleUp( const int& idx )
{
    if( idx < 0 || idx >= m_container.size() )
    {
        //exit(1);
        return; //20170413. chanyoung song
    }

    int currIndex       = idx;
    int indexOfParent   = 0;
    while( currIndex > 0 ) // Index of top is 0
    {
        indexOfParent = computeIndexOfParent( currIndex );

		//cy song : m_container[currIndex]->m_key < m_container[indexOfParent]->m_key 
		// -->  
        if(m_comp(*m_container[currIndex], *m_container[indexOfParent]) )
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



template<typename T, typename COMPARATOR>
void PriorityQueue<T, COMPARATOR>::bubbleDown( const int& idx )
{
    if( idx < 0 || idx >= m_container.size() )
    {
        //exit(1);
        return; //20170413. chanyoung song
    }

    int     currIndex         = idx;
    int     indexOfLeftChild  = 0;
    int     indexOfRightChild = 0;

    int     queueSize = m_container.size();
    bool    canBeChanged = true;
    while( canBeChanged )
    {
        int     numChild = getNumChildNode( currIndex );
        int     indexOfTargetChild = 0;

		switch (numChild)
		{

		case 0:
		{
			canBeChanged = false;
			break;
		}

		case 1:
		{
			indexOfLeftChild = computeIndexOfLeftChild(currIndex);

			//cy song : m_container[currIndex]->m_key > m_container[indexOfParent]->m_key 
			// -->
			if (m_comp(*m_container[indexOfLeftChild], *m_container[currIndex]))
			{
				swapNode(currIndex, indexOfLeftChild);
				currIndex = indexOfLeftChild;
			}
			canBeChanged = false;
			break;
		}

		case 2:
		{
			indexOfLeftChild = computeIndexOfLeftChild(currIndex);
			indexOfRightChild = indexOfLeftChild + 1;


			//cy song : m_container[indexOfLeftChild]->m_key < m_container[indexOfRightChild]->m_key
			// -->
			if (m_comp(*m_container[indexOfLeftChild], *m_container[indexOfRightChild]))
			{
				indexOfTargetChild = indexOfLeftChild;
			}
			else
			{
				indexOfTargetChild = indexOfRightChild;
			}

			//cy song : m_container[indexOfSmallerChild]->m_key < m_container[currIndex]->m_key
			// -->
			if (m_comp(*m_container[indexOfTargetChild], *m_container[currIndex]))
			{
				swapNode(currIndex, indexOfTargetChild);
				currIndex = indexOfTargetChild;
			}
			else{
				canBeChanged = false;
			}

			break;
		}

		default:
		{
            return;
			//exit(1);
			break;
		}

        }
    }
}



template<typename T, typename COMPARATOR>
void PriorityQueue<T, COMPARATOR>::swapNode( const int& idx1, const int& idx2 )
{
    PQ_Node<T>* tempNode = m_container[idx1];
    m_container[idx1] = m_container[idx2];
    m_container[idx2] = tempNode;

    m_container[idx1]->m_index = idx1;
    m_container[idx2]->m_index = idx2;
}



template<typename T, typename COMPARATOR>
int PriorityQueue<T, COMPARATOR>::getNumChildNode( const int& idx )
{
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




template<typename T, typename COMPARATOR>
void PriorityQueue<T, COMPARATOR>::findEntitiesHavingSamekey(const double& key, list<T>& outputEntities, const double& res /*=resNeg4*/)
{
	enum tag { WHITE, GRAY, BLACK };

	
	vector<tag> tags;

	for (int i = 0; i < m_container.size(); i++)
	{
		tags.push_back(WHITE);
	}



	int currIdx = 0;

	while (!((currIdx == 0) && (tags[0] == BLACK)))
	{
        if (rg_NE(m_container[currIdx]->m_key, key, res))
		//if (m_container[currIdx]->m_key != key)
		{
			tags[currIdx] = BLACK;
		}

		switch (tags[currIdx])
		{
		case WHITE:
		{
			outputEntities.push_back(m_container[currIdx]->m_entity);
			
			if (currIdx * 2 + 1 < m_container.size())
			{
				tags[currIdx] = GRAY;
				currIdx = currIdx * 2 + 1;
			}
			else
			{
				tags[currIdx] = BLACK;
				currIdx = (currIdx - 1) / 2;
			}

			break;
		}

		case GRAY:
		{
			if (currIdx * 2 + 2 < m_container.size())
			{
				tags[currIdx] = BLACK;
				currIdx = currIdx * 2 + 2;
			}
			else
			{
				tags[currIdx] = BLACK;
				currIdx = (currIdx - 1) / 2;
			}
			
			break;
		}

		case BLACK:
		{
			if (currIdx == 0)
			{
				currIdx = 0;
			}
			else
			{
				currIdx = (currIdx - 1) / 2;
			}
			
			break;
		}

		default:
			break;
		} 
	}

}




template<typename T, typename COMPARATOR /*= less<PQ_Node<T>>*/>
void PriorityQueue<T, COMPARATOR>::copyFrom(const PriorityQueue& pq)
{
	int Q_size = pq.size();

	for (int i = 0; i < Q_size; i++) {
		T&      entity = pq.m_container[i]->m_entity;   //cysong T -> T&
		double  key = pq.m_container[i]->m_key;
		int     idx = pq.m_container[i]->m_index;

		PQ_Node<T>* node = new PQ_Node<T>(entity, key, idx);
		m_container.push_back(node);
	}

	m_comp = pq.m_comp;
}


#endif

