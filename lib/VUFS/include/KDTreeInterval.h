#ifndef _KDTREEINTERVAL_H__
#define _KDTREEINTERVAL_H__

#include "rg_Const.h"
#include <float.h>

const rg_REAL KDTREE_REAL_PLUS_INFINITE  = DBL_MAX;
const rg_REAL KDTREE_REAL_MINUS_INFINITE = -DBL_MAX;

const int OPEN_OPEN   = 0;
const int OPEN_CLOSE  = 1;
const int CLOSE_OPEN  = 2;
const int CLOSE_CLOSE = 3;

const int INTERSECTION = 0;
const int IN_RANGE     = 1;
const int OUT_OF_RANGE = 2;

class KDTreeInterval
{
private:
	rg_REAL m_upperValue;
	rg_REAL m_lowerValue;
	rg_INT  m_statusOfOpenClose;

public:

    KDTreeInterval();
    KDTreeInterval(const rg_REAL& lowerValue, const rg_REAL& upperValue);
	KDTreeInterval(const rg_REAL& lowerValue, const rg_REAL& upperValue, const int& statusOfOpenClose);
	KDTreeInterval(const KDTreeInterval& interval);
	~KDTreeInterval();
	
	rg_REAL getUpperValue() const;
	rg_REAL getLowerValue() const;
	int     getStatusOfOpenClose() const;
	
	void setUpperValue(const rg_REAL& value);
	void setLowerValue(const rg_REAL& value);
    void setLowerAndUpperValue(const rg_REAL& lowerValue, const rg_REAL& upperValue);
    
    void setStatusOfOpenClose(const int& statusOfOpenClose);

	
	//rg_BOOL
	rg_FLAG isInInterval(const rg_REAL& value) const;	
	//
	rg_FLAG isInInterval(const KDTreeInterval& interval) const;
    rg_FLAG isIntersectInterval(const KDTreeInterval& interval) const;
	// 함수명 변경
    rg_INT  getIntersectionStatus(const KDTreeInterval& interval) const;

	KDTreeInterval& operator =(const KDTreeInterval& kdTreeInterval);
	
private:
	rg_FLAG isInIntervalWithOpenOpenEnd(const rg_REAL& value) const;
    rg_FLAG isInIntervalWithOpenCloseEnd(const rg_REAL& value) const;
    rg_FLAG isInIntervalWithCloseOpenEnd(const rg_REAL& value) const;
    rg_FLAG isInIntervalWithCloseCloseEnd(const rg_REAL& value) const;
};
#endif
