#include "Interval_REAL.h"
using namespace BULL2D::GeometryTier;


Interval_REAL::Interval_REAL()
: m_upperValue(REAL_MINUS_INFINITE), m_lowerValue(REAL_MINUS_INFINITE)
{	
}



Interval_REAL::Interval_REAL(const rg_REAL& lowerValue, const rg_REAL& upperValue)
: m_lowerValue(lowerValue), m_upperValue(upperValue)
{
}



Interval_REAL::Interval_REAL(const Interval_REAL& interval)
{
	m_upperValue = interval.m_upperValue;
	m_lowerValue = interval.m_lowerValue;
}



Interval_REAL::~Interval_REAL()
{
}



rg_REAL Interval_REAL::getUpperValue() const
{
	return m_upperValue;
}



rg_REAL Interval_REAL::getLowerValue() const
{
	return m_lowerValue;
}



void Interval_REAL::setUpperValue(const rg_REAL& value)
{
	m_upperValue = value;
}



void Interval_REAL::setLowerValue(const rg_REAL& value)
{
	m_lowerValue = value;
}


	
void Interval_REAL::setLowerAndUpperValue(const rg_REAL& lowerValue, const rg_REAL& upperValue)
{
	m_lowerValue = lowerValue;
	m_upperValue = upperValue;
}



rg_FLAG Interval_REAL::isInIntervalWithOpenOpenEnd(const rg_REAL& value) const
{
    if ( (m_lowerValue < value) && (value < m_upperValue ) )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_FLAG Interval_REAL::isInIntervalWithOpenCloseEnd(const rg_REAL& value) const
{
    if ( m_lowerValue == m_upperValue )  {
        if ( m_lowerValue == value )
            return rg_TRUE;
        else
            return rg_FALSE;
    }
    else {    
        if ( (m_lowerValue < value) && (value <= m_upperValue ) )
            return rg_TRUE;
        else
            return rg_FALSE;
    }
}



rg_FLAG Interval_REAL::isInIntervalWithCloseOpenEnd(const rg_REAL& value) const
{
    if ( m_lowerValue == m_upperValue )  {
        if ( m_lowerValue == value )
            return rg_TRUE;
        else
            return rg_FALSE;
    }
    else {    
        if ( (m_lowerValue <= value) && (value < m_upperValue ) )
            return rg_TRUE;
        else
            return rg_FALSE;
    }
}



rg_FLAG Interval_REAL::isInIntervalWithCloseCloseEnd(const rg_REAL& value) const
{
    if ( (m_lowerValue <= value) && (value <= m_upperValue ) )
        return rg_TRUE;
    else
        return rg_FALSE;
}



Interval_REAL& Interval_REAL::operator =(const Interval_REAL& interval)
{
	if( this == &interval )
		return *this;
	
	m_upperValue = interval.m_upperValue;
	m_lowerValue = interval.m_lowerValue;
  
	return *this;
}


