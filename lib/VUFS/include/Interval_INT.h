#ifndef _INTERVAL_INT_H__
#define _INTERVAL_INT_H__

#include "rg_Const.h"
#include <limits.h>

const rg_INT INT_PLUS_INFINITE  = INT_MAX;
const rg_INT INT_MINUS_INFINITE = INT_MIN;


class Interval_INT
{
private:
	rg_INT m_upperValue;
	rg_INT m_lowerValue;

public:

    Interval_INT();
	Interval_INT(const rg_INT& lowerValue, const rg_INT& upperValue);
	Interval_INT(const Interval_INT& interval);
	~Interval_INT();
	
	rg_INT getUpperValue() const;
	rg_INT getLowerValue() const;
	
	void setUpperValue(const rg_INT& value);
	void setLowerValue(const rg_INT& value);
	void setLowerAndUpperValue(const rg_INT& lowerValue, const rg_INT& upperValue);
	
	Interval_INT& operator =(const Interval_INT& interval);
};


#endif

