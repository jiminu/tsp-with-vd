/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_ComplexNumber.cpp
//	  
//    DESCRIPTION : 
//           This is the implementation of the class rg_ComplexNumber 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 18 Dev 1996    
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include "rg_ComplexNumber.h"

rg_ComplexNumber::rg_ComplexNumber()
{
	realNumber      = 0;
	imaginaryNumber = 0;
}

rg_ComplexNumber::rg_ComplexNumber(const rg_REAL &r, const rg_REAL &i)
{
	realNumber      = r;
	imaginaryNumber = i;
}

rg_ComplexNumber::rg_ComplexNumber(const rg_ComplexNumber &other)
{
	realNumber      = other.realNumber;
	imaginaryNumber = other.imaginaryNumber;
}

rg_ComplexNumber::~rg_ComplexNumber()
{

}

void rg_ComplexNumber::setRealNumber(const rg_REAL &r)
{
	realNumber = r;
}

void rg_ComplexNumber::setImaginaryNumber(const rg_REAL &i)
{
	imaginaryNumber = i;
}

void rg_ComplexNumber::setComplexNumber(const rg_ComplexNumber &other)
{
	realNumber      = other.realNumber;
	imaginaryNumber = other.imaginaryNumber;
}

rg_REAL rg_ComplexNumber::getRealNumber() const
{
	return realNumber;
}

rg_REAL rg_ComplexNumber::getImaginaryNumber() const
{
	return imaginaryNumber;
}

rg_ComplexNumber rg_ComplexNumber::getComplexNumber() const
{
	return rg_ComplexNumber(realNumber, imaginaryNumber);
}

rg_INT rg_ComplexNumber::isZero() const
{
	if(   rg_ZERO(realNumber)
	   && rg_ZERO(imaginaryNumber) )
	{
		return 1;
	}
	else return 0;
}

rg_INT rg_ComplexNumber::isPureRealNumber() const
{
	if( rg_ZERO(imaginaryNumber))
		return 1;
	else return 0;
}

rg_ComplexNumber rg_ComplexNumber::operator+(const rg_ComplexNumber &other)
{
	return rg_ComplexNumber(realNumber + other.realNumber, imaginaryNumber + other.imaginaryNumber);	
}

rg_ComplexNumber rg_ComplexNumber::operator-(const rg_ComplexNumber &other)
{
	return rg_ComplexNumber(realNumber - other.realNumber, imaginaryNumber - other.imaginaryNumber);
}

rg_ComplexNumber rg_ComplexNumber::operator*(const rg_ComplexNumber &other)
{
	rg_REAL rNumber = realNumber*other.realNumber - imaginaryNumber*other.imaginaryNumber;
	rg_REAL iNumber = realNumber*other.imaginaryNumber + imaginaryNumber*other.realNumber;

	return rg_ComplexNumber(rNumber, iNumber);
}

rg_ComplexNumber rg_ComplexNumber::operator/(const rg_ComplexNumber &other)
{
	rg_REAL rNumber = (realNumber*other.realNumber + imaginaryNumber*other.imaginaryNumber)/(other.realNumber*other.realNumber + other.imaginaryNumber*other.imaginaryNumber);
	//rg_REAL iNumber = (realNumber*other.imaginaryNumber - imaginaryNumber*other.realNumber)/(other.realNumber*other.realNumber + other.imaginaryNumber*other.imaginaryNumber);
	rg_REAL iNumber = (imaginaryNumber*other.realNumber - realNumber*other.imaginaryNumber)/(other.realNumber*other.realNumber + other.imaginaryNumber*other.imaginaryNumber);
	// modified by Joonghyun Feb. 20, 2010

	return rg_ComplexNumber(rNumber, iNumber);
}

rg_ComplexNumber rg_ComplexNumber::operator+(const rg_REAL &r)
{
	return rg_ComplexNumber(realNumber + r, imaginaryNumber);
}

rg_ComplexNumber rg_ComplexNumber::operator-(const rg_REAL &r)
{
	return rg_ComplexNumber(realNumber - r, imaginaryNumber);
}

rg_ComplexNumber rg_ComplexNumber::operator*(const rg_REAL &r)
{
	return rg_ComplexNumber(realNumber * r, imaginaryNumber * r);
}

rg_ComplexNumber rg_ComplexNumber::operator/(const rg_REAL &r)
{
	return rg_ComplexNumber(realNumber / r, imaginaryNumber / r);
}

rg_ComplexNumber operator*(const rg_REAL &r, rg_ComplexNumber &other)
{ 
	return other * r;
}

rg_ComplexNumber operator/(const rg_REAL &r, rg_ComplexNumber &other)
{ 
	rg_ComplexNumber temp(r, 0);
	return temp / other;
}

//unary operator
rg_ComplexNumber operator-(const rg_ComplexNumber &other)
{
	return rg_ComplexNumber(-other.realNumber, -other.imaginaryNumber);
}

rg_INT rg_ComplexNumber::operator==(const rg_ComplexNumber &other)
{
	if(rg_EQ(realNumber, other.realNumber) && rg_EQ(imaginaryNumber, other.imaginaryNumber))
		return 1;
	else return 0;
}

/*
//      | a=0    a<0    a>0
// -----+--------------------
// b=0  |  0      a      a
// b<0  |  bi    a-bi   a-bi
// b>0  |  bi    a+bi   a+bi 
ostream& operator<<(ostream& os, const rg_ComplexNumber &other)
{
	if(rg_NZERO(other.realNumber))
	{
		if(rg_ZERO(other.imaginaryNumber))
			os << other.realNumber;
		else if(rg_NEG(other.imaginaryNumber))
			os << other.realNumber << other.imaginaryNumber << "i";
		else
			os << other.realNumber << "+" << other.imaginaryNumber << "i";
	}
	
	else
	{
		if(rg_ZERO(other.imaginaryNumber))
			os << 0;
		else
			os << other.imaginaryNumber << "i";
	}

	return os;
}
*/


