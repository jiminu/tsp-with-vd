   /**********************************************************/
   /** This Program is Single &Double Linked List Module    **/
   /**********************************************************/

#include <iostream>
#include <stdlib.h>

//#include "..\header\rg_Const.h"
//#include "..\header\list.h"
#include "rg_ListByPtr.h"

rg_sNodeByPtr::rg_sNodeByPtr()
{
   e = rg_NULL;
   next = rg_NULL;	  //or this
}
rg_sNodeByPtr::rg_sNodeByPtr(void* a) //, rg_sNodeByPtr* n)
{
   e = a;
   //next = n;
   next = this;
}
rg_sNodeByPtr::~rg_sNodeByPtr()
{
}

//get functions
void*    rg_sNodeByPtr::getEntity()
{ 
   return e;
}
rg_sNodeByPtr* rg_sNodeByPtr::getNext()
{ 
   return next;
}
// set functions
void   rg_sNodeByPtr::setNext(rg_sNodeByPtr* thisNode) 
{ 
   next= thisNode;
}
////////////////////////////////////////////////
// member functions for a sinle-linked list
rg_sListByPtr::rg_sListByPtr()
{
   head = rg_NULL;
   current = rg_NULL;
}
rg_sListByPtr::rg_sListByPtr(void* a)
{ 
   head = new rg_sNodeByPtr(a);
   head->next = head;
}
rg_sListByPtr::~rg_sListByPtr() 
{ 
   killList(); 
}
rg_sNodeByPtr* rg_sListByPtr::getFirstNode() 
{ 
   return head;
}
rg_sNodeByPtr* rg_sListByPtr::getSecondNode() 
{
   if (head==rg_NULL) 
      return rg_NULL;
   else 
      return head->next;
}
rg_sNodeByPtr*   rg_sListByPtr::getLastNode()
{
   if (head == rg_NULL)
      return rg_NULL;

   rg_sNodeByPtr *curr = head;
   rg_sNodeByPtr* prev;
   do
   {
      prev = curr;
      curr = curr->next;
   } while(curr != head);
   
   return prev;
}

   ////////// set functions
void   rg_sListByPtr::setHead(rg_sNodeByPtr *thisNode) 
{ 
   head = thisNode; 
}
void   rg_sListByPtr::setCurrent(rg_sNodeByPtr *thisNode) 
{ 
   current = thisNode; 
}
void   rg_sListByPtr::reset() 
{ 
   current = rg_NULL; 
}
void   rg_sListByPtr::resetAll() 
{ 
   head = rg_NULL;  
   current = rg_NULL; 
}
rg_INT rg_sListByPtr::append(void* entity)
{
   //reset();
   
   if (head != rg_NULL)
   {
       rg_sNodeByPtr* lastNode = getLastNode();

       rg_sNodeByPtr* tempNode = new rg_sNodeByPtr(entity);
	   tempNode->next = lastNode->next;
       lastNode->next = tempNode;
   }
   else
   {
       head = new rg_sNodeByPtr(entity);
       //head->next = head;  done in rg_sNodeByPtr constructor
   }
   return rg_NULL;
}
//concatenate the data of list2 only into list1.
//the nodes of list2 and the list2 itself is not changed.
void rg_sListByPtr::concatenateElementOnlyNotNodes(rg_sListByPtr *list2)
{
   list2->reset();
   void *data = rg_NULL;
   while (data = (void *)list2->getNextEntity())
   {
	  append(data);
   }
   return;
}
//concatenate the nodes and elements of list2 into list1, and delete 
//the list2 itself only. NOT the nodes and elements of list2
void rg_sListByPtr::concatenate(rg_sListByPtr* list2)
{
   rg_sNodeByPtr* topNodeOfList1 = head;
   rg_sNodeByPtr* topNodeOfList2 = list2->getFirstNode();

   if ((topNodeOfList1) && (topNodeOfList2))
   {
      rg_sNodeByPtr* botNodeOfList1 = getLastNode();
      rg_sNodeByPtr* botNodeOfList2 = list2->getLastNode();
      botNodeOfList1->setNext(topNodeOfList2);
      botNodeOfList2->setNext(topNodeOfList1);
   }
   else if ((! topNodeOfList1) && (topNodeOfList2))
   {
      setHead(topNodeOfList2);
   }
   //else if ((topNodeOfList1) && (! topNodeOfList2))
   // or if ((! topNodeOfList1) && (! topNodeOfList2))
      //do nothing

   //if not resetAll, the nodes and entities of list2 will be deleted   
   list2->resetAll(); 
   delete list2;
   return;
}
//the nodes and elements of list2 unchanged
//only the nodes of list1 corresponding to elements of list2 are deleted
void rg_sListByPtr::deConcatenateNodeOnly(rg_sListByPtr *list2)
{			
   if (list2)
   {
	  rg_sNodeByPtr *headOfList2 = list2->getFirstNode();
	  rg_sNodeByPtr *currNodeOfList2 = headOfList2;
	  do
	  {
		 rg_sNodeByPtr *nextNodeOfList2 = currNodeOfList2->getNext();
		 killNodeOnly(currNodeOfList2);//kill from "this" list1
		 currNodeOfList2 = nextNodeOfList2;
	  }while (currNodeOfList2 != headOfList2);
   }
   return;
}
void rg_sListByPtr::deConcatenate(rg_sListByPtr *list2)//list2 is completely deleted
{
   deConcatenateNodeOnly(list2);
   list2->killList();
   delete list2;
   return;
}

rg_INT rg_sListByPtr::appendIfNotExist(void* entity)
{
   if (!exists(entity))
   {
      append(entity);
      return 1;
   }
   else 
      return rg_NULL;
}
rg_INT rg_sListByPtr::insert(void* entity, rg_sNodeByPtr* afterThisNode)
{
   rg_sNodeByPtr* tempNode = new rg_sNodeByPtr(entity);
   tempNode->next = afterThisNode->next;
   afterThisNode->next = tempNode;

   return rg_NULL;
}

void*     rg_sListByPtr::getFirstEntity()
{
   return head->e;
}

void*     rg_sListByPtr::getSecondEntity()
{
   return head->next->e;
}
void*    rg_sListByPtr::getLastEntity()
{ 
   return getLastNode()->getEntity(); 
}

void rg_sListByPtr::killNodes()
{
   rg_sNodeByPtr* curr = head;
   if (curr != rg_NULL)
   {
      do
      {
         rg_sNodeByPtr* prev = curr;
         curr = curr->next;
         delete prev;
      } while (curr != head);
      
      head = rg_NULL;
   }
   return;
}

void rg_sListByPtr::killEntities()
{
   rg_sNodeByPtr* curr = head;
   if (curr != rg_NULL)
   {
      do
      {
         delete curr->e;
         curr = curr->next;
      } while (curr != head);
   }
   return;
}

void rg_sListByPtr::killList()
{
   rg_sNodeByPtr* curr = head;
   if (curr != rg_NULL)
   {
      do
      {
         rg_sNodeByPtr* prev = curr;
         curr = curr->next;
         delete prev->e;
         delete prev;
      } while (curr != head);
      
      head = rg_NULL;
   }
   return;
}

void* rg_sListByPtr::getNextEntity()
{
   if (head == rg_NULL) 
      return rg_NULL;

   if (current == rg_NULL)
      current = head;
   else
   {
      current = current->next;
      current = (current == head) ? rg_NULL : current;
   }

   return current ? current->e : rg_NULL;
}

void*     rg_sListByPtr::getNextEntity(void* currEntity)
{
   rg_sNodeByPtr* currNode = head;

   do
   {
      if(currNode->e == currEntity)
         return currNode->next->e;
      
      currNode = currNode->next;
   } while (currNode != head);

   //temp 96.1.22
   return currNode->next->e;
}

void*     rg_sListByPtr::getNextEntity(rg_sNodeByPtr *currNode)
{
   return  currNode->next->e;
}

rg_INT rg_sListByPtr::getNumberOfEntity()
{
   rg_INT count = 0;

   if (head == rg_NULL)
        return count;
   else
   {
      rg_sNodeByPtr* tmp = head;
      //rg_sNodeByPtr* top = tmp;
      do
      {
         tmp = tmp->next;
         count++;
      } while(tmp != head);
      
      return count;
   }
}
//NOTE that the elements are not duplicated
//only the nodes are new nodes
rg_sListByPtr*   rg_sListByPtr::duplicateListWithSameElements()
{
   rg_sListByPtr *slist = new rg_sListByPtr;
   
   void *elem=rg_NULL;
   reset();
   while(elem = getNextEntity())
   {
   	  slist->append(elem);
   }  
   return slist;
}

void rg_sListByPtr::killOneNodeNentity(rg_sNodeByPtr *thisNode)
{
   rg_sNodeByPtr* now = head;
   if (now == rg_NULL)
      return;

   if(head== head->next)
   {
      delete   head->e;
      delete  head;
      head= rg_NULL;
      return;
   }

   rg_sNodeByPtr* prv = getLastNode();
   rg_INT Tag=0;

   do
   {
      if(now == thisNode)
      {
         if (thisNode == head) 
            head= thisNode->next;
         Tag = 1;
         prv->next = now->next;
         delete now->e;
         delete now;
      }
      else
      {
         prv = prv->next;
         now = now->next;
      }
   } while (Tag!=1);
}

void rg_sListByPtr::killOneNodeNentity(void* a)
{
   rg_sNodeByPtr* now = head;
   if (now == rg_NULL)
      return;
   if(head== head->next)
   {
      delete  head->e;
      delete  head;
      head= rg_NULL;
      return;
   }
   rg_sNodeByPtr* prv = getLastNode();
   rg_INT Tag=0;
   do
   {
      if(now->e == a)
      {
         if (now == head) 
            head= now->next;
         Tag = 1;
         prv->next = now->next;
         delete now->e;
         delete now;
      }
      else
      {
         prv = prv->next;
         now = now->next;
      }
   } while (Tag!=1);
}

void rg_sListByPtr::killNodeOnly(rg_sNodeByPtr *thisNode)
{
   rg_sNodeByPtr* now = head;
   if (now == rg_NULL)
      return;

   if(head == head->next)
   {
      delete  head;
      head= rg_NULL;
      return;
   }

   rg_sNodeByPtr* prv = getLastNode();
   rg_INT Tag=0;

   do
   {
      if(now == thisNode)
      {
         if (thisNode == head) 
            head = thisNode->next;
         Tag = 1;
         prv->next = now->next;
         delete now;
      }
      else
      {
         prv = prv->next;
         now = now->next;
      }
   } while (Tag!=1);
}

void rg_sListByPtr::killNodeOnly(void* a)
{
   if (head == rg_NULL)
      return;
   
   if ((head == head->next) && (head->getEntity() == a))
   {
      delete  head;
      head= rg_NULL;
      return;
   }

   rg_sNodeByPtr* now = head;
   rg_sNodeByPtr* prv = getLastNode();
   rg_INT Tag=0;
   do
   {
      if(now->e == a)
      {
         if (now == head)
            head= now->next;
         Tag = 1;
         prv->next = now->next;
         delete now;
      }
      else
      {
         prv = prv->next;
         now = now->next;
      }
   } while (Tag!=1);
}

rg_INT      rg_sListByPtr::inSequence(void* a, void* b)
{
   rg_INT count = 0;
   rg_INT aPosition = 0;
   rg_INT bPosition = 0;
   current = head;
   do
   { 
      count++;
      if(current->e == a)         
         aPosition = count;
      else if(current->e == b)    
         bPosition = count;
      
	  if (aPosition * bPosition) break;

      current = current->next;
   } while(current != head);

   if (aPosition < bPosition) 
      return 1;
   else                      
      return 0;
}
rg_INT      rg_sListByPtr::areNeighbors(void* e1, void* e2)
{
   if ((e1 == getNextEntity(e2)) || (e2 == getNextEntity(e1)))
       return 1;
   else 
       return 0;
}

/////////////////////////////////////////////////
// Double-linked lists

// rg_REAL-linked list NODE members
rg_dNodeByPtr::rg_dNodeByPtr(void* a, rg_dNodeByPtr* p, rg_dNodeByPtr* n)
{
   e = a;
   prev = p;
   next = n;
}
rg_dNodeByPtr::~rg_dNodeByPtr()
{
}
rg_dNodeByPtr* rg_dNodeByPtr::getNext()
{ 
   return next;
}
rg_dNodeByPtr* rg_dNodeByPtr::getPrev()
{ 
   return prev;
}
void* rg_dNodeByPtr::getEntity()
{ 
   return e;
}
void  rg_dNodeByPtr::setNext(rg_dNodeByPtr* thisNode) 
{ 
   next= thisNode;
}
void  rg_dNodeByPtr::setPrev(rg_dNodeByPtr* thisNode) 
{ 
   prev= thisNode;
}

   //////////////// constructor & destructor
rg_dListByPtr::rg_dListByPtr() { head = rg_NULL; }
rg_dListByPtr::rg_dListByPtr(void* a)
{
   head = new rg_dNodeByPtr(a);
   head->prev = head->next = head;
}
rg_dListByPtr::~rg_dListByPtr() { killList(); }

///////////////// get functions
rg_dNodeByPtr*   rg_dListByPtr::getLastNode()
{
   if (head == rg_NULL)
      return rg_NULL;
   else
      return head->prev;
}

/*----- Insert "a" behind now ----------------------------------------*/
rg_INT rg_dListByPtr::insert(void* a, rg_dNodeByPtr* now)
{
   if (head)
   {
                rg_dNodeByPtr* now_aft= now->next;
                rg_dNodeByPtr* tmp = new rg_dNodeByPtr(a, now, now->next);
                now->next = tmp;
                now_aft->prev = tmp;
   }
   else 
   {
                head = new rg_dNodeByPtr(a);
                head->prev = head->next = head;
   }
   return 0;
}

rg_INT rg_dListByPtr::append(void* a)
{
   if (head)
   {
      head->prev = new rg_dNodeByPtr(a, head->prev, head);
      head->prev->prev->next = head->prev;
   }
   else
   {
      head = new rg_dNodeByPtr(a);
      head->prev = head->next = head;
   }
   return 0;
}

rg_INT rg_dListByPtr::getNumberOfEntity()
{
    rg_INT count= 0;

    if (head == rg_NULL)
    {
       return count;
    }
    else
    {
       rg_dNodeByPtr* tmp = head;
       rg_dNodeByPtr* top = tmp;
       do
       {
          tmp = tmp->next;
          count++;
       } while(tmp!= top);
 
       return count;
    }
}
//NOTE that the elements are not duplicated
//only the nodes are new nodes
rg_dListByPtr*   rg_dListByPtr::duplicateListWithSameElements()
{
   rg_dListByPtr *dlist = new rg_dListByPtr;
   
   void *elem=rg_NULL;
   reset();
   while(elem = getNextEntity())
   {
   	  dlist->append(elem);
   }  
   return dlist;
}
void* rg_dListByPtr::getFirstEntityNdeleteNode()
{
   rg_dNodeByPtr* f = head;
   void* r = f->e;

   if (f == head->next)
      head = rg_NULL;
   else
   {
      head = f->prev->next = f->next;
      f->next->prev = f->prev ;
   }
   delete f;
   return r;
}

void*     rg_dListByPtr::getFirstEntity()
{
        return head->e;
}

void*     rg_dListByPtr::getSecondEntity()
{
        return head->next->e;
}

rg_INT      rg_dListByPtr::deleteEntity(rg_dNodeByPtr *node)
{
        if (node == head)
        {
                head = node->next;
        }
        node->prev->next = node->next;
        node->next->prev = node->prev;
        delete node->e;
        delete node;
        return 0;
}

void rg_dListByPtr::killNodes()
{
   rg_dNodeByPtr* curr = head;
   if (curr != rg_NULL)
   {
      do
      {
         rg_dNodeByPtr* prev = curr;
         curr = curr->next;
         delete prev;
      } while (curr != head);
      
      head = rg_NULL;
   }
   return;
}
void rg_dListByPtr::killEntities()
{
   rg_dNodeByPtr* curr = head;
   if (curr != rg_NULL)
   {
      do
      {
         delete curr->e;
         curr = curr->next;
      } while (curr != head);
   }
   return;
}

void rg_dListByPtr::killList()
{
   rg_dNodeByPtr* curr = head;
   if (curr != rg_NULL)
   {
      do
      {
         rg_dNodeByPtr* prev = curr;
         curr = curr->next;
         delete prev->e;
         delete prev;
      } while (curr != head);
      
      head = rg_NULL;
   }
   return;
}

void* rg_dListByPtr::getPrevEntity()
{
   if (current == rg_NULL)
   {
      current = head;
   }
   else
   {
      current = current->prev;
      current = (current == head) ? rg_NULL : current;
   }
   return current ? current->e : rg_NULL;
}


void* rg_dListByPtr::getNextEntity()
{
   if (head == rg_NULL) return rg_NULL;

   if (current == rg_NULL)
      current = head;
   else
   {
      current = current->next;
      current = (current == head) ? rg_NULL : current;
   }
   return current ? current->e : rg_NULL;
}

void*     rg_dListByPtr::getNextEntity(void *currEntity)
{
   rg_dNodeByPtr* currNode = head;

   do
   {
      if(currNode->e == currEntity)
         return currNode->next->e;
      
      currNode = currNode->next;
   } while (currNode != head);

   //temp 96.1.22
   return currNode->next->e;
}

void*     rg_dListByPtr::getNextEntity(rg_dNodeByPtr *currNode)
{
   return  currNode->next->e;
}

rg_dNodeByPtr* rg_dListByPtr::getFirstNode() 
{
   return head;
}
rg_dNodeByPtr* rg_dListByPtr::getSecondNode() 
{
   if (head==rg_NULL) 
      return rg_NULL;
   else 
      return head->next;
}
void*     rg_dListByPtr::getPrevEntity(void* cur)
{
   rg_dNodeByPtr* ll = head;

   do
   {
          if(ll->e == cur)
          {
                return ll->prev->e;
          }
      ll = ll->next;
   } while (ll != head);

//temp 96
   return ll->next->e;
}

void*     rg_dListByPtr::getPrevEntity(rg_dNodeByPtr *cur)
{
   return  cur->prev->e;
}

void rg_dListByPtr::killOneNodeNentity(rg_dNodeByPtr *node)
{
   rg_dNodeByPtr* ll = head;
   if (ll == rg_NULL) return;
   if(head== head->next)
   {
      delete head->e;
      delete head;
      head= rg_NULL;
      return;
   }

   if (node == head)
   {
      head = node->next;
   }
   node->prev->next = node->next;
   node->next->prev = node->prev;
   delete node->e;
   delete node;
}

void rg_dListByPtr::killOneNodeNentity(void* a)
{
   rg_dNodeByPtr* ll = head;
   if (ll == rg_NULL)
      return;
   if(head== head->next)
   {
                delete   head->e;
                delete  head;
                head= rg_NULL;
                return;
   }

   rg_INT   Tag=0;
   do
   {
          if(ll->e == a)
          {
                if(ll == head)
                {
                        head = ll->next;
                }
                ll->prev->next = ll->next;
                ll->next->prev = ll->prev;
                delete ll->e;
        delete ll;
                Tag=1;
          }
          else ll = ll->next;
   } while (Tag!=1);
}

void* rg_dListByPtr::getLastEntity()
{ 
   return getLastNode()->getEntity(); 
}
void rg_dListByPtr::killNodeOnly(rg_dNodeByPtr *node)
{
   if (head == rg_NULL) 
      return;
   
   if(head== head->next)
   {
      delete  head;
      head= rg_NULL;
      return;
   }

   if (node == head)
   {
      head = node->next;
      node->prev->next = node->next;
      node->next->prev = node->prev;
      delete node;
      return;
   }
}
   /////////// set functions
void   rg_dListByPtr::setHead(rg_dNodeByPtr *thisNode) 
{ 
   head = thisNode; 
}
void   rg_dListByPtr::setCurrent(rg_dNodeByPtr *thisNode) 
{ 
   current = thisNode; 
}
void rg_dListByPtr::reset() 
{ 
   current = rg_NULL; 
}
void rg_dListByPtr::resetAll() 
{ 
   head = rg_NULL; 
   current = rg_NULL; 
}
//concatenate the data of list2 only into list1.
//the nodes of list2 and the list2 itself is not changed.
void rg_dListByPtr::concatenateElementOnlyNotNodes(rg_dListByPtr *list2)
{
   list2->reset();
   void *data = rg_NULL;
   while (data = (void *)list2->getNextEntity())
   {
	  append(data);
   }
   return;
}
//concatenate the nodes and elements of list2 into list1, and delete 
//the list2 itself only. NOT the nodes and elements of list2
void rg_dListByPtr::concatenate(rg_dListByPtr *list2)
{
   rg_dNodeByPtr *topNodeOfList1 = head; 
   rg_dNodeByPtr *topNodeOfList2 = list2->getFirstNode(); 

   if ((topNodeOfList1) && (topNodeOfList2))
   {
      rg_dNodeByPtr *botNodeOfList1 = head->getPrev(); 
      rg_dNodeByPtr *botNodeOfList2 = list2->getLastNode(); 
   
      topNodeOfList1->setPrev(botNodeOfList2);
      botNodeOfList2->setNext(topNodeOfList1);
      botNodeOfList1->setNext(topNodeOfList2);
      topNodeOfList2->setPrev(botNodeOfList1);
   }
   else if ((! topNodeOfList1) && (topNodeOfList2))
   {
      setHead(topNodeOfList2);
   }
   //else if ((topNodeOfList1) && (! topNodeOfList2))
   // or if ((! topNodeOfList1) && (! topNodeOfList2))
      //do nothing

   //if not resetAll, the nodes and entities of list2 will be deleted   
   list2->resetAll(); 
   delete list2;
}
//
void rg_dListByPtr::deConcatenateNodeOnly(rg_dListByPtr *list2)//the nodes of list2 unchanged
{
   if (list2)
   {
	  rg_dNodeByPtr *headOfList2 = list2->getFirstNode();
	  rg_dNodeByPtr *currNodeOfList2 = headOfList2;
	  do
	  {
		 rg_dNodeByPtr *nextNodeOfList2 = currNodeOfList2->getNext();
		 killNodeOnly(currNodeOfList2->getEntity());//kill from "this" list1
		 currNodeOfList2 =	nextNodeOfList2;
	  }while (currNodeOfList2 != headOfList2);
   }
   return;
}
void rg_dListByPtr::deConcatenate(rg_dListByPtr *list2)//list2 is completely deleted
{
   deConcatenateNodeOnly(list2);
   list2->killList();
   delete list2;
   return;
}

void rg_dListByPtr::killNodeOnly(void* a)
{

   if (head == rg_NULL)
      return;

   if ((head == head->next) && (head->getEntity() == a))
   {
      delete  head;
      head= rg_NULL;
      return;
   }

   rg_dNodeByPtr* ll = head;
   rg_INT Tag=0;
   do
   {
          if(ll->e == a)
          {
                if(ll == head)
                {
                        head = ll->next;
                }
                ll->prev->next = ll->next;
                ll->next->prev = ll->prev;
                delete ll;
                Tag=1;
          }
          else ll = ll->next;
   } while (Tag!=1);
}

rg_INT      rg_dListByPtr::inSequence(void* a, void* b)
{
   rg_INT count = 0;
   rg_INT aPosition = 0;
   rg_INT bPosition = 0;
   current = head;
   do
   { 
      count++;
      if(current->e == a)         
         aPosition = count;
      else if(current->e == b)    
         bPosition = count;
      
	  if (aPosition * bPosition) break;

      current = current->next;
   } while(current != head);

   if (aPosition < bPosition) 
      return 1;
   else                      
      return 0;
}

// For STACK
void* rg_dListByPtr::pop()
// delete and return the top of rg_dListByPtr
{
        if (head == rg_NULL)
        {
                //cerr << "It's rg_NULL In POP\n";
                return rg_NULL;
        }
        else if (head->next == head)
        {
                void* elem = head->e;
                delete head;
                head = rg_NULL;
                return elem;
        }
        else
        {
                rg_dNodeByPtr*   tmp= head;
                head = tmp->next;
                head->prev = tmp->prev;
                tmp->prev->next = head;
                void*     elem= tmp->e;
                delete   tmp;
 
                return elem;
        }
}

void rg_dListByPtr::push(void* a)
//insert to the top of rg_dListByPtr
{
        append(a);
        head = head->prev;
}

void rg_dListByPtr::remove(void* a)
{
        killOneNodeNentity(a);
}

void rg_dListByPtr::enQueue(void* a)
//append to the list
{
        append(a); // Append entity to Last of List
}

void* rg_dListByPtr::deQueue()
// delete and return the top of rg_dListByPtr
{
        return pop(); // Delete Top of List and return that...
}

rg_INT      rg_sListByPtr::exists(void* a)
{
        if (head == rg_NULL) return 0;

        reset();
        void*      tmp;
        while(tmp = this->getNextEntity())
        {
                if (tmp == a) return 1;
        }
        return 0;
}
void*    rg_dListByPtr::getTop() { return getFirstEntity();}
void*    rg_dListByPtr::getBottom() { return getLastEntity();}
rg_INT      rg_dListByPtr::exists(void* a)
{
        if (head == rg_NULL) return 0;

        reset();
        void*      tmp;
        while(tmp = this->getNextEntity())
        {
                if (tmp == a) return 1;
        }
        return 0;
}

rg_INT      rg_dListByPtr::areNeighbors(void* e1, void* e2)
{
        if ((e1 == getNextEntity(e2)) || (e2 == getNextEntity(e1)))
           return 1;
        else 
           return 0;
}


