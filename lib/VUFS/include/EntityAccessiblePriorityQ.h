#ifndef _ENTITY_ACCESSIBLE_PRIORITY_Q
#define _ENTITY_ACCESSIBLE_PRIORITY_Q

//////////////////////////////////////////////////////////////////////////////
/*     2015.11.26. CY Song                                                  */
/*     Condition : Entity should be unique.                                 */ 
/*                 If not, it does not work according to what you expect    */
//////////////////////////////////////////////////////////////////////////////  


#include "PriorityQueue.h"

#if defined(__GNUC__) || defined (__GNUG__)
    #define GCC_VERSION (  __GNUC__ * 10000 \
                         + __GNUC_MINOR__ * 100 \
                         + __GNUC_PATCHLEVEL__)    
    #if( GCC_VERSION >= 40801)
            #define USING_UNORDERED_MAP
    #endif
#elif defined(_MSC_VER)
    #if( _MSC_VER >= 1800)
        #define USING_UNORDERED_MAP
    #endif
#endif

#ifdef USING_UNORDERED_MAP
#include <unordered_map>
#include <functional>
#else
#include <map>
#endif


#include <iostream>
using namespace std;


#ifdef USING_UNORDERED_MAP
template <typename T, typename COMPARATOR = less<PQ_Node<T>>, typename UMAP_Hash = hash<T>, typename UMAP_Pred = equal_to<T>> 
class EntityAccessiblePriorityQ: public PriorityQueue<T, COMPARATOR>
#else
template <typename T, typename COMPARATOR = less<PQ_Node<T>>, typename MAP_Comparator = less<T>>
class EntityAccessiblePriorityQ : public PriorityQueue<T, COMPARATOR>
#endif 
{
private:
	typedef PQ_Node<T> PQNode;


#ifdef USING_UNORDERED_MAP
	unordered_map<T, const PQNode*, UMAP_Hash, UMAP_Pred>  m_bridgeBetweenEntityAndNodeInPriorityQ;
#else
	map<T, const PQNode*, MAP_Comparator>  m_bridgeBetweenEntityAndNodeInPriorityQ;
#endif 


public:
	//constructor
	EntityAccessiblePriorityQ();
	EntityAccessiblePriorityQ(const T& entity, const double& key);
	EntityAccessiblePriorityQ(const EntityAccessiblePriorityQ& priorityQ);


	//deconstructor
	~EntityAccessiblePriorityQ();


	//operator
	EntityAccessiblePriorityQ& operator=(const EntityAccessiblePriorityQ& priorityQ);


	//function
	const PQNode*       push(const T& entity, const double& key);  //cysong PQ_Node<T>*  -> const  PQ_Node<T>* 
	const PQNode*       push(const pair<T, double>& aPair);        //cysong PQ_Node<T>*  -> const  PQ_Node<T>* 
	T                   pop();

	void		killNodeRelatedWith(const T& entity);
	void        replaceKey(const T& entity, const double& key);
	void        replaceEntity(const T& originalEntity, const T& newEntity);
	bool        doesHave(const T& entity) const;
    double      findKey(const T& entity)  const;

	void		clear();
	void        copyFrom(const EntityAccessiblePriorityQ& priorityQ);
};


#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::EntityAccessiblePriorityQ()
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::EntityAccessiblePriorityQ()
#endif 
    : PriorityQueue<T, COMPARATOR>()
{
}


 
#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::EntityAccessiblePriorityQ(const T& entity, const double& key)
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::EntityAccessiblePriorityQ(const T& entity, const double& key)
#endif 	
	: PriorityQueue<T, COMPARATOR>(entity, key)
{
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::EntityAccessiblePriorityQ(const EntityAccessiblePriorityQ& priorityQ)
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::EntityAccessiblePriorityQ(const EntityAccessiblePriorityQ& priorityQ)
#endif 	
{
	copyFrom(priorityQ);
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::~EntityAccessiblePriorityQ()
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::~EntityAccessiblePriorityQ()
#endif 	
{
	clear();
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>& EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::operator=(const EntityAccessiblePriorityQ& priorityQ)
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>& EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::operator=(const EntityAccessiblePriorityQ& priorityQ)
#endif 	
{
	if (this != &priorityQ)
	{
        clear();
        copyFrom(priorityQ);
	}

	return *this;
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	const PQ_Node<T>* EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::push(const pair<T, double>& aPair)
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	const PQ_Node<T>* EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::push(const pair<T, double>& aPair)
#endif 	
{
	const PQNode* node = PriorityQueue<T,COMPARATOR>::push(aPair);
	m_bridgeBetweenEntityAndNodeInPriorityQ.insert(pair<T, const PQNode*>(aPair.first, node));

	return node;
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	const PQ_Node<T>* EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::push(const T& entity, const double& key)
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	const PQ_Node<T>* EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::push(const T& entity, const double& key)
#endif 	
{
	const PQNode* node = PriorityQueue<T, COMPARATOR>::push(entity, key);
	m_bridgeBetweenEntityAndNodeInPriorityQ.insert(pair<T, const PQNode*>(entity, node));

	return node;
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	T EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::pop()
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	T EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::pop()
#endif 
{
	//T& entity = PriorityQueue::pop();
    T entity = PriorityQueue<T, COMPARATOR>::pop(); //20170413. chanyoung song
	m_bridgeBetweenEntityAndNodeInPriorityQ.erase(entity);

	return entity;
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	void EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::killNodeRelatedWith(const T& entity)
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	void EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::killNodeRelatedWith(const T& entity)
#endif 
{
    try {
        const PQNode* node = m_bridgeBetweenEntityAndNodeInPriorityQ.at(entity);

		PriorityQueue<T, COMPARATOR>::killNodeAt(node->getIndex());
        m_bridgeBetweenEntityAndNodeInPriorityQ.erase(entity);
    }
    catch (const std::out_of_range& oor) {
        std::cerr << "Out of Range error: " << oor.what() << '\n';
    }
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	void EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::replaceKey(const T& entity, const double& key)
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	void EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::replaceKey(const T& entity, const double& key)
#endif 
{
    try {
        const PQNode* node = m_bridgeBetweenEntityAndNodeInPriorityQ.at(entity);

		PriorityQueue<T, COMPARATOR>::changeKeyValue(node->getIndex(), key);
    }
    catch (const std::out_of_range& oor) {
        std::cerr << "Out of Range error: " << oor.what() << '\n';
    }
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	void EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::replaceEntity(const T& originalEntity, const T& newEntity)
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	void EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::replaceEntity(const T& originalEntity, const T& newEntity)
#endif 
{
    try {
        const PQNode* node = m_bridgeBetweenEntityAndNodeInPriorityQ.at(originalEntity);
		m_bridgeBetweenEntityAndNodeInPriorityQ.erase(originalEntity);
		m_bridgeBetweenEntityAndNodeInPriorityQ.insert(make_pair(newEntity, node));

		PriorityQueue<T, COMPARATOR>::changeEntity(node->getIndex(), newEntity);
    }
    catch (const std::out_of_range& oor) {
        std::cerr << "Out of Range error: " << oor.what() << '\n';
    }
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	bool EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::doesHave(const T& entity) const
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	bool EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::doesHave(const T& entity) const
#endif
{
	if (m_bridgeBetweenEntityAndNodeInPriorityQ.find(entity) != m_bridgeBetweenEntityAndNodeInPriorityQ.end())
	{
		return true;
	}
	else
	{
		return false;
	}

	/*  
	   [CYSONG, April27, 20] This is LINEAR SCAN algorithm so that this function is commented out.
		

#ifdef USING_UNORDERED_MAP

		for (typename unordered_map<T, const PQNode*, Hash>::const_iterator iter_pos = m_bridgeBetweenEntityAndNodeInPriorityQ.begin();
		iter_pos != m_bridgeBetweenEntityAndNodeInPriorityQ.end();
		iter_pos++)
		{
			if ((*iter_pos).first == entity){
				b_doesHave = true;
				break;
			}
		}

		

#else 

		for (typename map<T, const PQNode*>::const_iterator iter_pos = m_bridgeBetweenEntityAndNodeInPriorityQ.begin();
			iter_pos != m_bridgeBetweenEntityAndNodeInPriorityQ.end();
			iter_pos++)
		{
			if ((*iter_pos).first == entity)
			{
				b_doesHave = true;
				break;
			}
		}
#endif 

	return b_doesHave;
	*/
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	double EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::findKey(const T& entity) const
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	double EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::findKey(const T& entity) const
#endif
{
    try {
        return (m_bridgeBetweenEntityAndNodeInPriorityQ.at(entity))->getKey();      // vector::at throws an out-of-range
    }
    catch (const std::out_of_range& oor) {
        std::cerr << "Out of Range error: " << oor.what() << '\n';
    }
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	void EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::copyFrom(const EntityAccessiblePriorityQ& priorityQ)
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	void EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::copyFrom(const EntityAccessiblePriorityQ& priorityQ)
#endif
{
	PriorityQueue<T, COMPARATOR>::copyFrom(priorityQ);

	for (typename EntityAccessiblePriorityQ::iterator iter_node = PriorityQueue<T, COMPARATOR>::begin();
		iter_node != PriorityQueue<T, COMPARATOR>::end();
		iter_node++)
	{
		T&       entity  = (*iter_node)->m_entity;   //cysong T -> T&
		PQNode*  nodePtr = (*iter_node);

		m_bridgeBetweenEntityAndNodeInPriorityQ.insert(pair<T, PQNode*>(entity, nodePtr));
	}
}



#ifdef USING_UNORDERED_MAP
	template <typename T, typename COMPARATOR, typename UMAP_Hash, typename UMAP_Pred>
	void EntityAccessiblePriorityQ<T, COMPARATOR, UMAP_Hash, UMAP_Pred>::clear()
#else
	template <typename T, typename COMPARATOR, typename MAP_Comparator>
	void EntityAccessiblePriorityQ<T, COMPARATOR, MAP_Comparator>::clear()
#endif
{
    m_bridgeBetweenEntityAndNodeInPriorityQ.clear(); 
	PriorityQueue<T, COMPARATOR>::clear();
}


#endif

