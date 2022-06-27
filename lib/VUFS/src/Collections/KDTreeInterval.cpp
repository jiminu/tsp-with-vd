#include "KDTreeInterval.h"

KDTreeInterval::KDTreeInterval()
: m_upperValue(KDTREE_REAL_PLUS_INFINITE), m_lowerValue(KDTREE_REAL_MINUS_INFINITE), m_statusOfOpenClose(CLOSE_CLOSE)
{	
}

KDTreeInterval::KDTreeInterval(const rg_REAL& lowerValue, const rg_REAL& upperValue)
: m_lowerValue(lowerValue), m_upperValue(upperValue), m_statusOfOpenClose(CLOSE_CLOSE)
{
}

KDTreeInterval::KDTreeInterval(const rg_REAL& lowerValue, const rg_REAL& upperValue, const int& statusOfOpenClose)
: m_lowerValue(lowerValue), m_upperValue(upperValue), m_statusOfOpenClose(statusOfOpenClose)
{
}

KDTreeInterval::KDTreeInterval(const KDTreeInterval& interval)
{
	m_upperValue = interval.m_upperValue;
	m_lowerValue = interval.m_lowerValue;
	m_statusOfOpenClose = interval.m_statusOfOpenClose;
}



KDTreeInterval::~KDTreeInterval()
{
}



rg_REAL KDTreeInterval::getUpperValue() const
{
	return m_upperValue;
}



rg_REAL KDTreeInterval::getLowerValue() const
{
	return m_lowerValue;
}



int KDTreeInterval::getStatusOfOpenClose() const
{
    return m_statusOfOpenClose;
}



void KDTreeInterval::setUpperValue(const rg_REAL& value)
{
	m_upperValue = value;
}



void KDTreeInterval::setLowerValue(const rg_REAL& value)
{
	m_lowerValue = value;
}


	
void KDTreeInterval::setLowerAndUpperValue(const rg_REAL& lowerValue, const rg_REAL& upperValue)
{
	m_lowerValue = lowerValue;
	m_upperValue = upperValue;
}


void KDTreeInterval::setStatusOfOpenClose(const int& statusOfOpenClose)
{
    m_statusOfOpenClose = statusOfOpenClose;
}


rg_FLAG KDTreeInterval::isInInterval(const rg_REAL& value) const
{
    switch(m_statusOfOpenClose) {
        case OPEN_OPEN:
            return isInIntervalWithOpenOpenEnd(value);
            break;
        
        case OPEN_CLOSE:
            return isInIntervalWithOpenCloseEnd(value);
            break;
        
        case CLOSE_OPEN:
            return isInIntervalWithCloseOpenEnd(value);
            break;
        
        case CLOSE_CLOSE:
            return isInIntervalWithCloseCloseEnd(value);
            break;
        
        default:
            return isInIntervalWithCloseCloseEnd(value);
            break;
    }
}


rg_INT  KDTreeInterval::getIntersectionStatus(const KDTreeInterval& interval) const
{
    if(interval.getLowerValue() >= m_lowerValue && interval.getUpperValue() <= m_upperValue) {
        if(interval.getLowerValue() == m_lowerValue) {
            if((interval.getStatusOfOpenClose() == CLOSE_OPEN || interval.getStatusOfOpenClose() == CLOSE_CLOSE) &&
               (m_statusOfOpenClose             == OPEN_CLOSE || m_statusOfOpenClose             == OPEN_OPEN  )) {
                return INTERSECTION;
            }
        }
        
        if(interval.getUpperValue() == m_upperValue) {
            if((interval.getStatusOfOpenClose() == OPEN_CLOSE || interval.getStatusOfOpenClose() == CLOSE_CLOSE) &&
               (m_statusOfOpenClose             == CLOSE_OPEN || m_statusOfOpenClose             == OPEN_OPEN  )) {
                return INTERSECTION;
            }
        }
            
        return IN_RANGE;
    }
    else if(interval.getUpperValue() <= m_lowerValue || interval.getLowerValue() >= m_upperValue) {
        if(interval.getUpperValue() == m_lowerValue) {
            if((interval.getStatusOfOpenClose() == CLOSE_CLOSE || interval.getStatusOfOpenClose() == OPEN_CLOSE ) &&
               (m_statusOfOpenClose             == CLOSE_OPEN  || m_statusOfOpenClose             == CLOSE_CLOSE)) {
                return INTERSECTION;
            }
        }
        
        if(interval.getLowerValue() == m_upperValue) {
            if((interval.getStatusOfOpenClose() == CLOSE_OPEN || interval.getStatusOfOpenClose() == CLOSE_CLOSE) &&
               (m_statusOfOpenClose             == OPEN_CLOSE || m_statusOfOpenClose             == OPEN_CLOSE )) {
                return INTERSECTION;
            }
        }
                
        return OUT_OF_RANGE;
    }
    else {
        return INTERSECTION;
    }
    
    
  
}


rg_FLAG KDTreeInterval::isInInterval(const KDTreeInterval& interval) const
{   
    if(interval.getLowerValue() == m_lowerValue) {
        if((interval.getStatusOfOpenClose() == CLOSE_OPEN || interval.getStatusOfOpenClose() == CLOSE_CLOSE) &&
           (m_statusOfOpenClose             == OPEN_CLOSE || m_statusOfOpenClose             == OPEN_OPEN  )) {
            return rg_FALSE;
        }
    }
    
    if(interval.getUpperValue() == m_upperValue) {
        if((interval.getStatusOfOpenClose() == OPEN_CLOSE || interval.getStatusOfOpenClose() == CLOSE_CLOSE) &&
           (m_statusOfOpenClose             == CLOSE_OPEN || m_statusOfOpenClose             == OPEN_OPEN  )) {
            return rg_FALSE;
        }
    }
           
        
    if(interval.getLowerValue() >= m_lowerValue && interval.getUpperValue() <= m_upperValue) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }  
}


rg_FLAG KDTreeInterval::isIntersectInterval( const KDTreeInterval& interval ) const
{
    if(interval.getUpperValue() == m_lowerValue) {
        if((interval.getStatusOfOpenClose() == CLOSE_CLOSE || interval.getStatusOfOpenClose() == OPEN_CLOSE ) &&
           (m_statusOfOpenClose             == CLOSE_OPEN  || m_statusOfOpenClose             == CLOSE_CLOSE)) {
            return rg_TRUE;
        }
    }
    
    if(interval.getLowerValue() == m_upperValue) {
        if((interval.getStatusOfOpenClose() == CLOSE_OPEN || interval.getStatusOfOpenClose() == CLOSE_CLOSE) &&
           (m_statusOfOpenClose             == OPEN_CLOSE || m_statusOfOpenClose             == OPEN_CLOSE )) {
            return rg_TRUE;
        }
    }
    
    
    if(interval.getUpperValue() <= m_lowerValue || interval.getLowerValue() >= m_upperValue) {
        return rg_FALSE;
    }
    else {
        return rg_TRUE;
    } 
}


rg_FLAG KDTreeInterval::isInIntervalWithOpenOpenEnd(const rg_REAL& value) const
{
    if ( (m_lowerValue < value) && (value < m_upperValue ) )
        return rg_TRUE;
    else
        return rg_FALSE;
}



rg_FLAG KDTreeInterval::isInIntervalWithOpenCloseEnd(const rg_REAL& value) const
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



rg_FLAG KDTreeInterval::isInIntervalWithCloseOpenEnd(const rg_REAL& value) const
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



rg_FLAG KDTreeInterval::isInIntervalWithCloseCloseEnd(const rg_REAL& value) const
{
    if ( (m_lowerValue <= value) && (value <= m_upperValue ) )
        return rg_TRUE;
    else
        return rg_FALSE;
}



KDTreeInterval& KDTreeInterval::operator=( const KDTreeInterval& kdTreeInterval )
{
	if(this == &kdTreeInterval)
		return *this;

	m_upperValue = kdTreeInterval.m_upperValue;
	m_lowerValue = kdTreeInterval.m_lowerValue;
	m_statusOfOpenClose = kdTreeInterval.m_statusOfOpenClose;

	return *this;
}
