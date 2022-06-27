#ifndef _TARRAY_H_
#define _TARRAY_H_
#include "rg_Const.h"

template <class T>
class rg_Array
{
private:
	T* arr;
	rg_INT size;
	rg_FLAG checkBound(const rg_INT& index) const;

public:
	//constructors
	rg_Array();
	rg_Array(const rg_Array<T> &array);
	rg_Array(const rg_INT& tSize);
	rg_Array(const T* array, const rg_INT& tSize);

	//destructor
	~rg_Array();

	//gets
	rg_INT getSize() const;
	rg_INT getUpperBound() const;
	T getAt(const rg_INT& index) const;
    T*  getPointerOfArray() const;
    T&  elementAt(const rg_INT& index) const;

	//sets
	void setSize(const rg_INT& tSize);
	void setAt(const rg_INT& index, const T& value);
	void setAll(const T& value);

	//etc
	void removeAll();
	void add(const T& value);
	void append(const rg_Array<T> &array);
	void insertAt(const rg_INT &index, const T &value);
	void insertAt(const rg_INT &index, const rg_Array<T> &array);
	void insertAt(const rg_INT &index, const T* array, const rg_INT &tSize);
	void removeAt(const rg_INT &index);
	void copy(const rg_Array<T> &array);
	void copy(const T* const array, const rg_INT &tSize);
    void reverse();

	//overloaded operator
	rg_Array<T>& operator=(const rg_Array<T> &array);
	T& operator[](rg_INT index);
	operator T* (void) const;

};

template <class T>
rg_FLAG rg_Array<T>::checkBound(const rg_INT& index) const
{
	if( index > size-1 || index <0 )
	{
		return rg_FALSE;	//rg_FALSE
	}
	else
	{
		return rg_TRUE;   //rg_TRUE
	}
}

//constructor
template <class T>
rg_Array<T>::rg_Array(void)
{
	arr=rg_NULL;
	size=0;
}

template <class T>
rg_Array<T>::rg_Array(const rg_Array<T> &array)
{
	size=array.size;
	arr = new T[size];

	for(rg_INT i=0; i<size; i++)
	{
		arr[i]=array.arr[i];
	}
}

template <class T>
rg_Array<T>::rg_Array(const rg_INT& tSize)
{
	size=tSize;
	arr=new T[size];
}

template <class T>
rg_Array<T>::rg_Array(const T* array, const rg_INT& tSize)
{
	size=tSize;
	arr=new T [size];
	for(rg_INT i=0; i<size; i++)
	{
		arr[i]=array[i];
	}
}


//destructor
template <class T>
rg_Array<T>::~rg_Array(void)
{
	if(arr!=rg_NULL)
	{
		delete [] arr;
	}
}


//operator overloading
template <class T>
rg_Array<T>& rg_Array<T>::operator=(const rg_Array<T>& array)
{
	if(&array==this)  
	{
		return *this;
	}
	else if(array.size == 0)	//CASE: assigned by rg_NULL array
	{
		if(arr!=rg_NULL)
		{
			delete [] arr;
		}
		size=0;
		return *this;
	}

	if(arr!=rg_NULL)
	{
		delete [] arr;
	}

	size=array.size;
	arr=new T[size];
	for(rg_INT i=0; i<size; i++)
	{
		arr[i]=array.arr[i];
	}
	return *this;
}


template <class T>
T& rg_Array<T>::operator[](rg_INT index)
{
	return arr[index];
}

template <class T>
rg_Array<T>::operator T* (void) const
{
	return arr;
}


//gets
template <class T>
rg_INT rg_Array<T>::getSize() const
{
	return size;
}

template <class T>
rg_INT rg_Array<T>::getUpperBound() const
{
	return size-1;
}

template <class T>
T rg_Array<T>::getAt(const rg_INT& index) const
{
	return arr[index];
}
template <class T>
T&  rg_Array<T>::elementAt(const rg_INT& index) const
{
    return arr[index];
}

template <class T>
T* rg_Array<T>::getPointerOfArray() const
{
    return arr;
}

//sets
template <class T>
void rg_Array<T>::setSize(const rg_INT& tSize)
{
	rg_Array<T> temp(tSize);

	rg_INT low;
	low=(size< tSize ? size : tSize);

	rg_INT i = 0;
	for(i=0; i< low; i++)
	{
		temp[i]=arr[i];
	}
	delete [] arr;
	arr=new T[tSize];

	for(i=0; i<tSize; i++)
	{
		arr[i]=temp[i];
	}
	size=tSize;
}

template <class T>
void rg_Array<T>::setAt(const rg_INT& index, const T& value)
{
	if( checkBound(index) )
	{
		arr[index]=value;
	}
}

template <class T>
void rg_Array<T>::setAll(const T& value)
{
	for(rg_INT i=0; i<size; i++)
	{
		arr[i]=value;
	}
}


//etc
template <class T>
void rg_Array<T>::removeAll()
{
	if(arr!=rg_NULL)
	{
		delete [] arr;
		arr=rg_NULL;
		size=0;
	}
}
	

template <class T>
void rg_Array<T>::add(const T &value)
{
	if(size==0)
	{
		size=1;
		arr=new T [1];
		arr[0]=value;
		return;
	}
	
	T* newArray= new T [size+1];

	//assign existed value to new array
	for(rg_INT i=0; i<size; i++)
	{
		newArray[i]=arr[i];
	}
	newArray[size]=value;

	size++;
	delete [] arr;
	arr=newArray;
}

template <class T>
void rg_Array<T>::append(const rg_Array<T>& array)
{
	if(&array==this)
	{
		rg_Array<T> tempArray(array);
		append(tempArray);
	}
	
	T* newArray = new T [size+array.size];

	//assign existed value to new array
	rg_INT i = 0;
	for(i=0; i<size; i++)
	{
		newArray[i]= arr[i];
	}
	
	//assign new value
	for(i=size; i<size+array.size; i++)
	{
		newArray[i]=array[i-size];
	}

	size+=array.size;
	delete [] arr;
	arr=newArray;
}


template <class T>
void rg_Array<T>::insertAt(const rg_INT &index, const T &value)
{
	if(	!checkBound(index) )
	{
		return;
	}

	//CASE : insert at rg_NULL array
	if(arr==rg_NULL && index==0)
	{
		size=1;
		arr=new T [size];
		arr[0]=value;
		return;
	}

	T* newArray = new T [size+1];
	rg_INT i = 0;
	for(i=0; i<index; i++)
	{
		newArray[i]=arr[i];
	}
	newArray[index]=value;
	for(i=index+1; i<size+1; i++)
	{
		newArray[i]=arr[i-1];
	}

	size++;
	delete [] arr;
	arr=newArray;
}


template <class T>
void rg_Array<T>::insertAt(const rg_INT &index, const rg_Array<T> &array)
{
	if(	!checkBound(index) )
	{
		return;
	}

	if(arr==rg_NULL && index==0)
	{
		size=array.size;
		arr=new T [size];
		for(rg_INT i=0; i<size; i++)
		{
			arr[i]=array[i];
		}
	}

	T* newArray = new T [size+ array.size];
	rg_INT i = 0;
	for(i=0; i<index; i++)	
	{
		newArray[i]=arr[i];
	}
	for(i=index; i<index+array.size; i++)
	{
		newArray[i]=array[i-index];
	}
	for(i=index+array.size; i<size+array.size; i++)
	{
		newArray[i]=arr[i-array.size];
	}

	size+=array.size;
	delete [] arr;
	arr=newArray;
}


template <class T>
void rg_Array<T>::insertAt(const rg_INT &index, const T* const array, const rg_INT &tSize)
{
	if(	!checkBound(index) )
	{
		return;
	}

	if(arr==rg_NULL)
	{
		size=tSize;
		arr=new T [size];
		for(rg_INT i=0; i<size; i++)
		{
			arr[i]=array[i];
		}
	}

	T* newArray = new T [size+ tSize];
	rg_INT i = 0;
	for(i=0; i<index; i++)	
	{
		newArray[i]=arr[i];
	}
	for(i=index; i<index+tSize; i++)
	{
		newArray[i]=array[i-index];
	}
	for(i=index+tSize; i<size+tSize; i++)
	{
		newArray[i]=arr[i-tSize];
	}

	size+=tSize;
	delete [] arr;
	arr=newArray;
}


template <class T>
void rg_Array<T>::removeAt(const rg_INT &index)
{
	if(	!checkBound(index) || size==0 )
	{
		return;
	}

	if(size==1)
	{
		delete [] arr;
		size==0;
	}
	
	T* newArray= new T [size-1];

	rg_INT i = 0;
	for(i=0; i<index; i++)
	{
		newArray[i]=arr[i];
	}
	for(i=index; i<size-1; i++)
	{
		newArray[i]=arr[i+1];
	}

	size--;
	delete [] arr;
	arr=newArray;
}

template <class T>
void rg_Array<T>::copy(const rg_Array<T> &array)
{
	if(&array==this)
	{
		return;
	}

	if(arr!=rg_NULL)
	{
		delete [] arr;
	}
	size=array.size;

	if(size > 0)
	{
		arr=new T [size];
		for(rg_INT i=0; i<size; i++)
		{
			arr[i]=array[i];
		}
	}
	else
	{
		arr=rg_NULL;
	}
}


template <class T>
void rg_Array<T>::copy(const T* const array, const rg_INT &tSize)
{
	if(arr!=rg_NULL)
	{
		delete [] arr;
	}

	size=tSize;
	arr=new T [size];
	for(rg_INT i=0; i<size; i++)
	{
		arr[i]=array[i];
	}
}

template <class T>
void rg_Array<T>::reverse()
{
    const rg_INT middle= ( size%2 == 0 ) ? size/2:size/2+1;

    for( rg_INT i=0; i < middle; i++ )
    {
        T temp=arr[i];
        arr[i]=arr[size-1-i];
        arr[size-1-i]=temp;
    }
}


#endif

