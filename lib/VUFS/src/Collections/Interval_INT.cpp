#include "Interval_INT.h"

Interval_INT::Interval_INT()
{	
    m_lowerValue = INT_MINUS_INFINITE;
    m_upperValue = INT_PLUS_INFINITE;
}



Interval_INT::Interval_INT(const rg_INT& lowerValue, const rg_INT& upperValue)
{
    m_lowerValue = lowerValue;
    m_upperValue = upperValue;
}



Interval_INT::Interval_INT(const Interval_INT& interval)
{
	m_upperValue = interval.m_upperValue;
	m_lowerValue = interval.m_lowerValue;
}



Interval_INT::~Interval_INT()
{
}



rg_INT Interval_INT::getUpperValue() const
{
	return m_upperValue;
}



rg_INT Interval_INT::getLowerValue() const
{
	return m_lowerValue;
}



void Interval_INT::setUpperValue(const rg_INT& value)
{
	m_upperValue = value;
}



void Interval_INT::setLowerValue(const rg_INT& value)
{
	m_lowerValue = value;
}


	
void Interval_INT::setLowerAndUpperValue(const rg_INT& lowerValue, const rg_INT& upperValue)
{
	m_upperValue = upperValue;
	m_lowerValue = lowerValue;
}



Interval_INT& Interval_INT::operator =(const Interval_INT& interval)
{
	if( this == &interval )
		return *this;
	
	m_upperValue = interval.m_upperValue;
	m_lowerValue = interval.m_lowerValue;
  
	return *this;
}


