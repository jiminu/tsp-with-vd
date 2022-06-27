#ifndef _SARRAY_H_TEMPLATE
#define _SARRAY_H_TEMPLATE

#include "rg_Const.h"

template <class T>
class rg_SignArray
{
private:
	T 	   *arr;
	rg_INT		sign;
	rg_INT		size;
public:
	rg_SignArray();
	rg_SignArray( const rg_INT &sg, 
		    const rg_INT &sn );
	~rg_SignArray();

	void setArray( const rg_INT &sg, const rg_INT &sn);
	T &operator []( rg_INT i );
	T operator []( rg_INT i ) const;
};

#endif

//////////////////////////////////////////////////////////
//  Implementation of template class

template <class T>
rg_SignArray<T>::rg_SignArray()
{
	arr  = rg_NULL;
	size = 0;
	sign = 0;
}

template <class T>
rg_SignArray<T>::rg_SignArray( const rg_INT &sg, 
			    const rg_INT &sn )
{
	sign = sg;
	size = sn;
	arr  = new T[size];
}


template <class T>
rg_SignArray<T>::~rg_SignArray()
{
	if(arr != rg_NULL) delete arr;
}

template <class T>
void rg_SignArray<T>::setArray( const rg_INT &sg, const rg_INT &sn)
{
	if(arr!=rg_NULL) return;

	sign = sg;
	size = sn;
	arr  = new T[size];
}

template <class T>
T &rg_SignArray<T>::operator []( rg_INT i )
{
	return arr[i-sign];
}

template <class T>
T rg_SignArray<T>::operator []( rg_INT i ) const
{
	return arr[i-sign];
}

