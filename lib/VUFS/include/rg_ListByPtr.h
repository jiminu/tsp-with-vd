   /*********************************************************/
   /** This is HEADER FILE for Single & Double Linked List **/
   /*********************************************************/

#ifndef  _RG_LISTBYPTR_H
#define  _RG_LISTBYPTR_H
#include <stdlib.h>
//#include "..\header\rg_Const.h"
#include "rg_Const.h"

// NODE class for a single-linked list class
class rg_sNodeByPtr
{
   friend class rg_sListByPtr;

   private:
      rg_sNodeByPtr *next;
      void  *e;
   public :
   // CONSTRUCTORS AND DESTRUCTOR
	  rg_sNodeByPtr();
      //rg_sNodeByPtr(void* a, rg_sNodeByPtr* n = 0);
      rg_sNodeByPtr(void *a); //, rg_sNodeByPtr* n = 0);
	  ~rg_sNodeByPtr();
      void*    getEntity();
      rg_sNodeByPtr*   getNext();
      void     setNext(rg_sNodeByPtr *thisNode);
};

class rg_sListByPtr
{
   private:
      rg_sNodeByPtr *head;
      rg_sNodeByPtr *current;
   public :
   // CONSTRUCTORS AND DESTRUCTOR
      rg_sListByPtr();
      rg_sListByPtr(void *a);
      ~rg_sListByPtr();
   
   ////////// get functions
      void*    getNextEntity();
      void*    getNextEntity(void *currEntity);
      void*    getNextEntity(rg_sNodeByPtr *currNode);

      rg_sNodeByPtr*   getFirstNode();
      rg_sNodeByPtr*   getSecondNode(); 
      rg_sNodeByPtr*   getLastNode();

      void*    getFirstEntity();
      void*    getSecondEntity();
      void*    getLastEntity();
      
      rg_INT      getNumberOfEntity();
	  rg_sListByPtr*   duplicateListWithSameElements();
   ////////// set functions
      void     setHead(rg_sNodeByPtr *thisNode);
      void     setCurrent(rg_sNodeByPtr *thisNode);
      void     reset();
      void     resetAll();

   ////////// etc functions
      rg_INT      append(void *entity);
      rg_INT      appendIfNotExist(void *entity);
      rg_INT      insert(void *entity, rg_sNodeByPtr *afterThisNode);
      rg_INT      inSequence(void *a, void *b);//tests if a comes before b
      rg_INT      areNeighbors(void *e1, void *e2);
      rg_INT      exists(void *a);

      void     concatenateElementOnlyNotNodes(rg_sListByPtr *list2);//the nodes of list2 unchanged
      void     concatenate(rg_sListByPtr *list2);				 //list2 is completely deleted
	  void     deConcatenateNodeOnly(rg_sListByPtr *list2);//the nodes of list2 unchanged
	  void     deConcatenate(rg_sListByPtr *list2);              //list2 is completely deleted

      void     killList();     //delete all nodes and entities
      void     killNodes();    //delete all nodes only
      void     killEntities(); //delete all entities only
      void     killOneNodeNentity(rg_sNodeByPtr *thisNode); //delete a node and an entity
      void     killOneNodeNentity(void *a);		 //delete a node and an entity
      void     killNodeOnly(rg_sNodeByPtr *thisNode);
      void     killNodeOnly(void *a);

};

/*----- Double Linked List Class Define . ----------------------------*/
class rg_dNodeByPtr
{
   friend class rg_dListByPtr;
   private: 
      rg_dNodeByPtr *prev;
      rg_dNodeByPtr *next;
      void  *e;
   public :
      rg_dNodeByPtr(void *a, rg_dNodeByPtr *p=rg_NULL, rg_dNodeByPtr *n=rg_NULL);
	  ~rg_dNodeByPtr();
      rg_dNodeByPtr*    getNext();
      rg_dNodeByPtr*    getPrev();
      void*     getEntity();
      void      setNext(rg_dNodeByPtr *thisNode);
      void      setPrev(rg_dNodeByPtr *thisNode);
};

class rg_dListByPtr
{
   private:
      rg_dNodeByPtr *head;
      rg_dNodeByPtr *current;
   public :
   //////////////// constructor & destructor
      rg_dListByPtr();
      rg_dListByPtr(void *a);
      ~rg_dListByPtr();

   ////////// get functions
      void*    getNextEntity();
      void*    getNextEntity(void *currEntity);
      void*    getNextEntity(rg_dNodeByPtr *currNode);

      void*    getPrevEntity();
      void*    getPrevEntity(void *now);
      void*    getPrevEntity(rg_dNodeByPtr *now);

      void*    getFirstEntityNdeleteNode(); //get data element and delete the node

      void*    getFirstEntity();
      void*    getSecondEntity();
      void*    getLastEntity();

      rg_dNodeByPtr*   getFirstNode();
      rg_dNodeByPtr*   getSecondNode(); 
      rg_dNodeByPtr*   getLastNode();

      //body is not defined yet.
      rg_INT      getNumberOfEntity();


	  rg_dListByPtr*   duplicateListWithSameElements();

   /////////// set functions
      void     setHead(rg_dNodeByPtr *thisNode);
      void     setCurrent(rg_dNodeByPtr *thisNode);

      void     reset();
      void     resetAll();

      rg_INT      insert(void *a, rg_dNodeByPtr *now);
      rg_INT      append(void *a);

      //list1 + list2
      void     concatenateElementOnlyNotNodes(rg_dListByPtr *list2);//the nodes of list2 unchanged
      void     concatenate(rg_dListByPtr *list2);				 //list2 is completely deleted
	  void     deConcatenateNodeOnly(rg_dListByPtr *list2);//the nodes of list2 unchanged
	  void     deConcatenate(rg_dListByPtr *list2);              //list2 is completely deleted
	   
      rg_INT      deleteEntity(rg_dNodeByPtr *node);

      rg_INT      inSequence(void *a, void *b);//tests if a comes before b
      rg_INT      areNeighbors(void *e1, void *e2);
      rg_INT      exists(void *a);

      void     killList();     //delete all nodes and entities
      void     killNodes();    //delete all nodes only
      void     killEntities(); //delete all entities only
      void     killOneNodeNentity(rg_dNodeByPtr *dk); //delete a node and an entity
      void     killOneNodeNentity(void *a);   //delete a node and an entity
      void     killNodeOnly(rg_dNodeByPtr *dk);
      void     killNodeOnly(void *a);

      // considering this list as if a QUEUE
      void     enQueue(void*);
      void*    deQueue();

      // considering this list as if a STACK
      void     push(void*);
      void*    pop();
      void     remove(void*);  //killOneNodeNentity

      // common operations for both STACK & QUEUE
      void*    getTop();
      void*    getBottom();
};


#endif
