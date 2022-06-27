#ifndef _SORTFUNC_EXT
#define _SORTFUNC_EXT


#include "rg_Const.h"

template <class T>
void Swap(T& a, T& b);


// QuickSort accepts an array and two range parameters
template <class T>
void QuickSort(T A[], rg_INT low, rg_INT high)
{
   // local variables holding the mid index of the range,
   // its value A[mid] and the scanning indices
   T  pivot;
   rg_INT scanUp, scanDown;
   rg_INT mid;

   // if the range is not at least two elements, return  
   if (high - low <= 0)
	return;
   else 
   // if sublist has two elements, compare them and
   // exchange their values if necessary
   if (high - low == 1)
   {           
      if (A[high] < A[low]) 
         Swap(A[low], A[high]);
      return;
   }
   
   // get the mid index and assign its value to pivot
   mid = (low + high)/2;
   pivot = A[mid];

   // exchange the pivot and the low end of the range
   // and initialize the indices scanUp and scanDown.
   Swap(A[mid], A[low]);
   scanUp = low + 1;
   scanDown = high;
      
   // manage the indices to locate elements that are in
   // the wrong sublist; stop when scanDown < scanUp 
   do 
   {
      // move up lower sublist; stop when scanUp enters
   	  // upper sublist or identifies an element > pivot
   	  while (scanUp <= scanDown && A[scanUp] <= pivot)
         scanUp++;
 
      // scan down upper sublist; stop when scanDown locates 
      // an element <= pivot; we guarantee we stop at A[low]
      while (pivot < A[scanDown])
         scanDown--;
	         
      // if indices are still in their sublists, then they
      // identify two elements in wrong sublists. exchange
      if (scanUp < scanDown)
         Swap(A[scanUp], A[scanDown]);
   } 
   while (scanUp < scanDown);

   // copy pivot to index (scanDown) that partitions sublists 
   A[low] = A[scanDown];
   A[scanDown] = pivot;

   // if the lower sublist (low to scanDown-1) has 2 or more 
   // elements, make the recursive call
   if (low < scanDown-1)
      QuickSort(A, low, scanDown-1);

   // if higher sublist (scanDown+1 to high) has 2 or more 
   // elements, make the recursive call
   if (scanDown+1 < high)
      QuickSort(A, scanDown+1, high);
}

// the "QuickSort" function should have the "Swap" function
// this function exchange two type of value.
template <class T>
void Swap(T& a, T& b)
{
	T temp;
	temp = a;
	a    = b;
	b    = temp;
}




//#include "../Collections/sortFunc.cpp"
//
//template <class T>
//void QuickSort(T A[], rg_INT low, rg_INT high);
//
//template <class T>
//void Swap(T& a, T& b);


#endif // sort functions 

