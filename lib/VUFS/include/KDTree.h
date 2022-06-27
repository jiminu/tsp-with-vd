///////////////////////////////////////////////////////////////////////
//                                                                        
//    FILENAME    : KDTree.h
//    
//    DESCRIPTION : 
//           This consists of the definition and interface
//                                 of class KDTree using template.  
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Jae-Kwan Kim
//    START DATE  : 2010. 02. 26    
//
//            Copyright (c) Voronoi Diagram Research Center
//////////////////////////////////////////////////////////////////////

#ifndef _KDTREE_
#define _KDTREE_

#include "rg_dList.h"
#include "KDTreeInterval.h"

#include <algorithm>
#include <fstream>
using namespace std;

template<class T>
class node_less 
{
public:
    static int axis;
    bool operator()(const pair<T, double*>& x, const pair<T, double*>& y) {
        return x.second[axis] < y.second[axis];
    }
};


template <class T>
class KDNode
{
private:
    int             m_dimension;
    int             m_axis;
    KDTreeInterval* m_range;
    double*         m_values;
    
    KDNode*     m_parent;
    KDNode*     m_leftChild;
    KDNode*     m_rightChild;
    
    T           m_entity;

public:
    //////// constructor ////////////////////////
    KDNode();
    KDNode(const int& dimension, const int& axis, double* values, const T& entity);
    KDNode(const KDNode<T>& kdNode);
    ~KDNode();


    //////// get functions //////////////////////
    int             getDimension() const;
    int             getAxis() const;
    KDTreeInterval* getRange();
    double*         getValues();

    KDNode<T>*  getParent();
    KDNode<T>*  getLeftChild();
    KDNode<T>*  getRightChild();

    T           getEntity();
    T*          getpEntity();


    //////// set functions /////////////////////
    void setDimension(const int& dimension);
    void setAxis(const int& axis);
    void setRange(KDTreeInterval* range);
    void setValues(double* valueInNode);

    void setParent(KDNode<T>* parent);
    void setLeftChild(KDNode<T>* leftChild);
    void setRightChild(KDNode<T>* rightChild);
    
    void setEntity(const T& entity);	
};


template <class T>
class KDTree
{
private:
    int         m_dimension;
    KDNode<T>*  m_rootNode;

    int         m_numberOfNodes;
    
    
public:
    //////// constructor ////////////////////////
    KDTree();    
//	KDTree(const KDTree<T>& kdTree);
    KDTree(rg_dList<T>& entities, rg_dList<double*>& values, const int& dimension);
    ~KDTree();
    
    
	void setAllEntitiesFromList(rg_dList<T>& entities, rg_dList<double*>& values, const int& dimension);

    //////// search finctions ///////////////////
    void        getEntitiesInGivenRange(KDTreeInterval* range, KDNode<T>* rootNode, rg_dList<T>& entities);
	//void        getEntitiesInGivenRange(KDTreeInterval* range, rg_dList<T>& entities);
    KDNode<T>*  getRootNode();
    int         getSize() const;
    
    
private:
    KDNode<T>* makeTree(pair<T, double*>* inputData, const int& numOfData, const int& axis);
    void getSortedDataByGivenAxis(pair<T, double*>* inputData, 
                                  const int& startIndex, const int& endIndex, 
                                  const int& axis,
                                  pair<T, double*>*& data, int& numOfData);
    
    void quickSort(pair<T, double*>* data, const int& startOfIndex, const int& endOfIndex, const int& axis);
    
    void computeRangeOfNode(KDNode<T>* node, KDTreeInterval*& range);

    void addAllEntitiesGivenTree(KDNode<T>* rootNode, rg_dList<T>& entities);

    void killNode(KDNode<T>* node);

//	void duplicateTree(const KDTree<T>& kdTree);
    
};




template <class T>
KDNode<T>::KDNode()
{
    m_dimension     = 0;
    m_axis          = 0;
    m_range         = NULL;
    m_values        = NULL;
    
    m_parent        = NULL;
    m_leftChild     = NULL;
    m_rightChild    = NULL;
}



template <class T>
KDNode<T>::KDNode(const int& dimension, const int& axis, double* values, const T& entity)
{
    m_dimension     = dimension;
    m_axis          = axis;
    
    m_values = new double[m_dimension];

    int i = 0;
    for(i = 0; i < m_dimension; i++) {
        m_values[i] = values[i];
    }

    m_range         = NULL;

    m_entity        = entity;

    m_parent        = NULL;
    m_leftChild     = NULL;
    m_rightChild    = NULL;
}



template <class T>
KDNode<T>::KDNode(const KDNode<T>& kdNode)
{
    m_dimension      = kdNode.m_dimension;
    m_axis           = kdNode.m_axis;


    m_range = new KDTreeInterval[m_dimension];

    int i = 0;
    for(i = 0; i < m_dimension; i++) {        
        m_range[i].setLowerAndUpperValue(m_range[i].getLowerValue(), kdNode.m_range[i].getUpperValue());        
    }

    
    m_values = new double[m_dimension];
    
    for(i = 0; i < m_dimension; i++) {
        m_values[i] = kdNode.m_values[i];
    }    
    
    m_parent         = kdNode.m_parent;
    m_leftChild      = kdNode.m_leftChild;
    m_rightChild     = kdNode.m_rightChild;
    
    m_entity         = kdNode.m_entity;
}


template <class T>
KDNode<T>::~KDNode()
{
    if(m_range != NULL) {
        delete[] m_range;
    }

    if(m_values != NULL) {
        delete[] m_values;
    }
}




template <class T>
int KDNode<T>::getDimension() const
{
    return m_dimension;
}


template <class T>
int KDNode<T>::getAxis() const
{
    return m_axis;
}


template <class T>
KDTreeInterval* KDNode<T>::getRange()
{
    return m_range;
}


template <class T>
double* KDNode<T>::getValues()
{
    return m_values;
}



template <class T>
KDNode<T>* KDNode<T>::getParent()
{
    return m_parent;
}


template <class T>
KDNode<T>* KDNode<T>::getLeftChild()
{
    return m_leftChild;
}


template <class T>
KDNode<T>* KDNode<T>::getRightChild()
{
    return m_rightChild;
}



template <class T>
T KDNode<T>::getEntity()
{
    return m_entity;
}


template <class T>
T* KDNode<T>::getpEntity()
{
    return &m_entity;
}




template <class T>
void KDNode<T>::setDimension(const int& dimension)
{
    m_dimension = dimension;
}



template <class T>
void KDNode<T>::setAxis(const int& axis)
{
    m_axis = axis;
}



template <class T>
void KDNode<T>::setRange(KDTreeInterval* range)
{
    if(m_range != NULL) {
        delete[] m_range;
    }

    m_range = new KDTreeInterval[m_dimension];

    int i = 0;
    for(i = 0; i < m_dimension; i++) {
        m_range[i].setLowerAndUpperValue(range[i].getLowerValue(), range[i].getUpperValue());        
    }
}



template <class T>
void KDNode<T>::setValues(double* values)
{
    if(m_values != NULL) {
        delete[] m_values;
    }


    m_values = new double[m_dimension];
    
    int i = 0;
    for(i = 0; i < m_dimension; i++) {
        m_values[i] = values[i];
    }
}



template <class T>
void KDNode<T>::setParent(KDNode<T>* parent)
{
    m_parent = parent;
}


template <class T>
void KDNode<T>::setLeftChild(KDNode<T>* leftChild)
{
    m_leftChild = leftChild;
}

template <class T>
void KDNode<T>::setRightChild(KDNode<T>* rightChild)
{
    m_rightChild = rightChild;
}


template <class T>
void KDNode<T>::setEntity(const T& entity)
{
    m_entity = entity;
}









template <class T>
KDTree<T>::KDTree()
{
    m_dimension = 0;
    m_rootNode  = NULL;
}    


//template <class T>
//KDTree<T>::KDTree( const KDTree<T>& kdTree )
//{
//	m_dimension = 0;
//    m_rootNode  = NULL;
//    duplicateTree(kdTree);
//}



template <class T>
KDTree<T>::KDTree(rg_dList<T>& entities, rg_dList<double*>& values, const int& dimension)
{
	setAllEntitiesFromList(entities, values, dimension);
}






template <class T>
KDTree<T>::~KDTree()
{
    if(m_rootNode != NULL) {
        killNode(m_rootNode);
    }
}


template <class T>
void KDTree<T>::setAllEntitiesFromList( rg_dList<T>& entities, rg_dList<double*>& values, const int& dimension )
{
    m_numberOfNodes = 0;
    m_dimension = dimension;

    
    pair<T, double*>* inputData;
    inputData = new pair<T, double*>[entities.getSize()];
    
    int index = 0;
    
    entities.reset4Loop();
    values.reset4Loop();
    
    while(entities.setNext4Loop() && values.setNext4Loop()) {
        T currEntity     = entities.getEntity();
        double* valueArr = values.getEntity();
        
        inputData[index].first  = currEntity;
        inputData[index].second = valueArr;
        
        index++;
    }

    node_less<T>::axis = 0;
    sort(&inputData[0], &inputData[entities.getSize()], node_less<T>() );

    m_rootNode = makeTree(inputData, entities.getSize(), 0);
 
    delete[] inputData;




// 	m_numberOfNodes = 0;
//     m_dimension = dimension;
// 
//     
//     pair<T, double*>* inputData;
//     inputData = new pair<T, double*>[entities.getSize()];
//     
//     int index = 0;
//     
//     entities.reset4Loop();
//     values.reset4Loop();
//     
//     while(entities.setNext4Loop() && values.setNext4Loop()) {
//         T currEntity     = entities.getEntity();
//         double* valueArr = values.getEntity();
//         
//         inputData[index].first  = currEntity;
//         inputData[index].second = valueArr;
//         
//         index++;
//     }
//     
//     int numOfData;
//     pair<T, double*>* sortedData;
//     
//     getSortedDataByGivenAxis(inputData, 0, entities.getSize() - 1, 0, sortedData, numOfData);
//     m_rootNode = makeTree(sortedData, numOfData, 0);
// 
//     delete[] inputData;
}



template <class T>
int KDTree<T>::getSize() const
{
    return m_numberOfNodes;    
}



template <class T>
void KDTree<T>::getEntitiesInGivenRange(KDTreeInterval* range, KDNode<T>* rootNode, rg_dList<T>& entities)
{
    KDTreeInterval* rangeOfRootNode = rootNode->getRange();
    
    int intersectionStatus = IN_RANGE;

    int i = 0;
    
    for(i = 0; i < m_dimension; i++) {        
        switch(range[i].getIntersectionStatus(rangeOfRootNode[i])) {
            case INTERSECTION:
                intersectionStatus = INTERSECTION;
                break;
            case IN_RANGE:
                break;
            case OUT_OF_RANGE:
                intersectionStatus = OUT_OF_RANGE;
                break;
            default:
                break;       
        }
        
        if(intersectionStatus == OUT_OF_RANGE) {
            break;
        }
    }

    bool isThisNodeInRang = true;
        
    double* currNodeValue = rootNode->getValues();

    switch(intersectionStatus) {
        case INTERSECTION:
            for(i = 0; i < m_dimension; i++) {
                if(range[i].isInInterval(currNodeValue[i]) == rg_FALSE) {
                    isThisNodeInRang = false;
                    break;
                }
            }
            
            if(isThisNodeInRang == true) {
                entities.addTail(rootNode->getEntity());
            }
            
            
            if(rootNode->getLeftChild() != NULL) {
                getEntitiesInGivenRange(range, rootNode->getLeftChild(), entities);
            }
            
            if(rootNode->getRightChild() != NULL) {
                getEntitiesInGivenRange(range, rootNode->getRightChild(), entities);
            }
            break;
        case IN_RANGE:
            addAllEntitiesGivenTree(rootNode, entities);
            break;
        case OUT_OF_RANGE:            
            break;
        default:
            break;       
    }
}



template <class T>
void KDTree<T>::addAllEntitiesGivenTree( KDNode<T>* rootNode, rg_dList<T>& entities )
{
     KDNode<T>* currNode = rootNode;

     entities.addTail(currNode->getEntity());

     if(currNode->getLeftChild() != NULL) {
        addAllEntitiesGivenTree(currNode->getLeftChild(), entities);
     }

     if(currNode->getRightChild() != NULL) {
         addAllEntitiesGivenTree(currNode->getRightChild(), entities);
     }
}



template <class T>
KDNode<T>* KDTree<T>::makeTree(pair<T, double*>* inputData, const int& numOfData, const int& axis)
{
    KDNode<T>* currNode;

    m_numberOfNodes++;

    int midIndex = numOfData/2;        
    currNode = new KDNode<T>(m_dimension, axis, inputData[midIndex].second, inputData[midIndex].first);
    
    int axisForChild = (axis + 1) % m_dimension;        

    pair<T, double*>* leftData  = NULL;
    pair<T, double*>* rightData = NULL;
    int numOfLeftData  = 0;
    int numOfRightData = 0;

    if(midIndex - 1 >= 0) {
        getSortedDataByGivenAxis(inputData, 0, midIndex - 1, axisForChild, leftData, numOfLeftData);
        KDNode<T>* leftChild = makeTree(leftData, numOfLeftData, axisForChild);
        leftChild->setParent(currNode);
        currNode->setLeftChild(leftChild);
    }
    
    if(midIndex + 1 < numOfData) {
        getSortedDataByGivenAxis(inputData, midIndex + 1, numOfData - 1, axisForChild, rightData, numOfRightData);
        KDNode<T>* rightChild = makeTree(rightData, numOfRightData, axisForChild);
        rightChild->setParent(currNode);
        currNode->setRightChild(rightChild);
    }


    KDTreeInterval* range;    
    computeRangeOfNode(currNode, range);

    currNode->setRange(range);


    delete[] leftData;
    delete[] rightData;


    return currNode;    
}   



template <class T>
void KDTree<T>::getSortedDataByGivenAxis(pair<T, double*>* inputData, 
                                         const int& startIndex, const int& endIndex, 
                                         const int& axis,
                                         pair<T, double*>*& data, int& numOfData)
{
    int i = 0;
    
    numOfData = endIndex - startIndex + 1;
    data = new pair<T, double*>[numOfData];

    for(i = 0; i < numOfData; i++) {
        data[i] = inputData[i + startIndex];
    }


    node_less<T>::axis = axis;
    sort(&data[0], &data[numOfData], node_less<T>() );

    //quickSort(data, 0, numOfData - 1, axis);
}



template <class T>
void KDTree<T>::computeRangeOfNode(KDNode<T>* node, KDTreeInterval*& range)
{
    double* currNodeValues = node->getValues();

    range = new KDTreeInterval[m_dimension];

    int i = 0;

    for(i = 0; i < m_dimension; i++) {        
        range[i].setLowerAndUpperValue(currNodeValues[i], currNodeValues[i]);    
        range[i].setStatusOfOpenClose(CLOSE_CLOSE);
    }
    
        
    if(node->getLeftChild() != NULL) {
        KDTreeInterval* rangeForLeft  = node->getLeftChild()->getRange();
        
        double min, max;
 
        for(i = 0; i < m_dimension; i++) {
            min = range[i].getLowerValue() < rangeForLeft[i].getLowerValue() ? range[i].getLowerValue() : rangeForLeft[i].getLowerValue();
            max = range[i].getUpperValue() > rangeForLeft[i].getUpperValue() ? range[i].getUpperValue() : rangeForLeft[i].getUpperValue();
            
            range[i].setLowerAndUpperValue(min, max);            
        }
    }
    
    if(node->getRightChild() != NULL) {
        KDTreeInterval* rangeForRight  = node->getRightChild()->getRange();
        
        double min, max;
        
        for(i = 0; i < m_dimension; i++) {
            min = range[i].getLowerValue() < rangeForRight[i].getLowerValue() ? range[i].getLowerValue() : rangeForRight[i].getLowerValue();
            max = range[i].getUpperValue() > rangeForRight[i].getUpperValue() ? range[i].getUpperValue() : rangeForRight[i].getUpperValue();
            
            range[i].setLowerAndUpperValue(min, max);         
        }
    }
}


template <class T>
void KDTree<T>::quickSort(pair<T, double*>* data, const int& startOfIndex, const int& endOfIndex, const int& axis)
{
    if (endOfIndex -  startOfIndex < 1)
    {
        return;
    }
     

    pair<T, double*> tempData;
    
    if (endOfIndex - startOfIndex == 1 )
    {
        if (data[startOfIndex].second[axis] > data[endOfIndex].second[axis])
        {
            tempData = data[startOfIndex];
            data[startOfIndex] = data[endOfIndex];
            data[endOfIndex] = tempData;
        }
        else
        {
        }		
    }
    else if (endOfIndex - startOfIndex > 1)
    {
        int pivot = (rand() % (endOfIndex - startOfIndex + 1)) + startOfIndex;
        tempData = data[startOfIndex];
        data[startOfIndex] = data[pivot];
        data[pivot] = tempData;
        
        int i = startOfIndex + 1;
        
        for (int j = startOfIndex + 1 ; j < endOfIndex + 1 ; j++)
        {
            if (data[startOfIndex].second[axis] > data[j].second[axis])
            {
                tempData = data[j];
                data[j] = data[i];
                data[i] = tempData;
                i++;
            }
            else
            {				
            }			
        }
        
        tempData = data[i-1];
        data[i-1] = data[startOfIndex];
        data[startOfIndex] = tempData; 
        
        
        quickSort(data, startOfIndex, i - 2, axis);        
        quickSort(data, i , endOfIndex, axis);
        
    }
    else
    {
    }
}



template <class T>
void KDTree<T>::killNode(KDNode<T>* node) 
{
    if(node->getLeftChild() != NULL) {
        killNode(node->getLeftChild());
    }

    if(node->getRightChild() != NULL) {
        killNode(node->getRightChild());
    }

    delete node;
}


template <class T>
KDNode<T>* KDTree<T>::getRootNode()
{
    return m_rootNode;    
}



//template <class T>
//void KDTree<T>::duplicateTree( const KDTree<T>& kdTree )
//{
// 	if ( size != 0 )
//         killNode(m_rootNode);
// 
//     rg_dNode<T> *ptr=thisList.head;
// 
//     for( rg_INDEX i=0; i < thisList.size; i++)
//     {
//         addTail(ptr->entity);
//         ptr=ptr->next;
//     }
//}

#endif
