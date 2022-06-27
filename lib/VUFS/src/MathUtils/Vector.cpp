#include "Vector.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////////
// constructors and destructor
//////////////////////////////////////////////////////////////////////////
Vector::Vector()
{
	m_size = 0;
	m_element = rg_NULL;
}

Vector::Vector(const rg_REAL* element, const rg_INT& size)
{
	m_size = size;
	m_element = new rg_REAL[m_size];

	rg_INDEX i;
	for (i = 0;i < m_size;i++)
	{
		m_element[ i ] = element[ i ];
	}
}

Vector::Vector(const rg_REAL& oneDimPoint)
{
	m_size = 1;
	m_element = new rg_REAL[m_size];
	m_element[ 0 ] = oneDimPoint;
}

Vector::Vector(const Vector& tVector)
{
	m_size = tVector.m_size;
	m_element = new rg_REAL[m_size];

	for(rg_INT i = 0;i < m_size;i++)
		m_element[ i ] = tVector.m_element[ i ];
}

Vector::Vector(const rg_INT& tSize)
{
	m_size = tSize;
	m_element = new rg_REAL[m_size];
	for(rg_INT i = 0;i < m_size;i++)
		m_element[ i ] = 0.;
}

Vector::~Vector()
{
	if(m_element != rg_NULL)
	{
		delete[] m_element;
		m_element = rg_NULL;
	}
	m_size = 0;
}

//////////////////////////////////////////////////////////////////////////
// get and set functions
//////////////////////////////////////////////////////////////////////////
rg_INT Vector::getSize() const
{
	return m_size;
}

rg_REAL* Vector::getElement()
{
	return m_element;
}

rg_INT Vector::getElement(rg_REAL*& element)
{
	for(rg_INT i = 0;i < m_size;i++)
		element[ i ] = m_element[ i ];

	return m_size;
}

rg_REAL Vector::getElement(const rg_INT& i) const
{
	return m_element[ i ];
}

void Vector::setSize(const rg_INT& tSize)
{
	m_size = tSize;
	if(m_element != rg_NULL)
		delete [] m_element;

	m_element = new rg_REAL[m_size];
	for(rg_INT i = 0;i < m_size;i++)
		m_element[ i ] = 0.;
}

rg_REAL Vector::innerProduct(const Vector& tVector) const
{
	rg_REAL result = 0.;
	for(rg_INT i = 0;i < m_size;i++)
		result += (m_element[ i ] * tVector.m_element[ i ]);
	return result;
}

rg_REAL Vector::getVectorNorm() const
{
	rg_REAL norm = innerProduct(*this);

	return sqrt(norm);
}

rg_REAL Vector::getSquaredNorm() const
{
	rg_REAL norm = innerProduct(*this);

	return norm;
}

Vector Vector::getNormalizedForm() const
{
	Vector result(m_size);
	rg_REAL mag = getVectorNorm();

	for(rg_INT i = 0;i < m_size;i++)
		result.m_element[ i ] = m_element[ i ] / mag;

	return result;
}

void Vector::normalize()
{
	rg_REAL mag = getVectorNorm();

	for(rg_INT i = 0;i < m_size;i++)
		m_element[ i ] = m_element[ i ] / mag;
}

//////////////////////////////////////////////////////////////////////////
// overloaded operators
//////////////////////////////////////////////////////////////////////////
rg_REAL& Vector::operator[](const rg_INT& Index)
{
	return m_element[Index];
}

Vector& Vector::operator=(const Vector& tVector)
{
	if(this == &tVector)
	{
		return *this;
	}
	else
	{
		setSize(tVector.m_size);
		for(rg_INT i = 0;i < m_size;i++)
			m_element[ i ] = tVector.m_element[ i ];

		return *this;
	}
}

Vector Vector::operator+(const Vector& tVector) const
{
	Vector result(m_size);

	for(rg_INT i = 0;i < m_size;i++)
		result.m_element[ i ] = m_element[ i ] + tVector.m_element[ i ];

	return result;
}

Vector Vector::operator-() const
{
	Vector result(m_size);

	for(rg_INT i = 0;i < m_size;i++)
		result.m_element[ i ] = - m_element[ i ];

	return result;
}

Vector Vector::operator-(const Vector& tVector) const
{
	Vector result(m_size);

	for(rg_INT i = 0;i < m_size;i++)
		result.m_element[ i ] = m_element[ i ] - tVector.m_element[ i ];

	return result;
}

Vector Vector::operator*(const Vector& tVector) const
{
	Vector result(m_size);

	for(rg_INT i = 0;i < m_size;i++)
		result.m_element[ i ] = m_element[ i ] * tVector.m_element[ i ];

	return result;
}

Vector operator*(const rg_REAL& scalar, const Vector& tVector)
{
	rg_INT m_size = tVector.getSize();
	Vector result(m_size);

	for(rg_INT i = 0;i < m_size;i++)
		result.m_element[ i ] = tVector.m_element[ i ] * scalar;

	return result;
}

Vector operator*(const Vector& tVector, const rg_REAL& scalar)
{
	return scalar * tVector;
}
