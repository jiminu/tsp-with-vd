#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "rg_Const.h"

//  This is the interface of the class Vector 
//  which stores and manipulates data in one-dimensional array 

class Vector
{
private:
	rg_INT     m_size;
	rg_REAL*   m_element;

public:
	//////////////////////////////////////////////////////////////////////////
	// constructors and destructor
	//////////////////////////////////////////////////////////////////////////
	Vector();
	Vector(const rg_REAL* element, const rg_INT& size);
	Vector(const rg_REAL& oneDimPoint);
	Vector(const Vector& tVector);
	Vector(const rg_INT& tSize);
	~Vector();

	//////////////////////////////////////////////////////////////////////////
	// get and set functions
	//////////////////////////////////////////////////////////////////////////
	rg_INT   getSize() const;
	rg_REAL* getElement();
	rg_INT   getElement(rg_REAL*& element);
	rg_REAL  getElement(const rg_INT& i) const;
	void     setSize(const rg_INT& tSize);

	rg_REAL  innerProduct(const Vector& tVector) const;
	rg_REAL  getVectorNorm() const;
	rg_REAL  getSquaredNorm() const;
	Vector   getNormalizedForm() const;
	void     normalize();

	//////////////////////////////////////////////////////////////////////////
	// overloaded operators
	//////////////////////////////////////////////////////////////////////////
	rg_REAL& operator[](const rg_INT& Index);
	Vector&  operator=(const Vector& tVector);
	Vector   operator+(const Vector& tVector) const;
	Vector   operator-() const;
	Vector   operator-(const Vector& tVector) const;
	Vector   operator*(const Vector& tVector) const;

	friend Vector operator*(const rg_REAL& scalar, const Vector& tVector);
	friend Vector operator*(const Vector& tVector, const rg_REAL& scalar);
};

#endif