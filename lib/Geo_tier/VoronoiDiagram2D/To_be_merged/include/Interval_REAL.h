#ifndef _INTERVAL_REAL_H__
#define _INTERVAL_REAL_H__

#include "rg_Const.h"
#include <float.h>

const rg_REAL REAL_PLUS_INFINITE  = DBL_MAX;
const rg_REAL REAL_MINUS_INFINITE = -DBL_MAX;

namespace BULL2D {
namespace GeometryTier {

class Interval_REAL
{
private:
	rg_REAL m_upperValue;
	rg_REAL m_lowerValue;

public:

    Interval_REAL();
	Interval_REAL(const rg_REAL& lowerValue, const rg_REAL& upperValue);
	Interval_REAL(const Interval_REAL& interval);
	~Interval_REAL();
	
	rg_REAL getUpperValue() const;
	rg_REAL getLowerValue() const;
	
	void setUpperValue(const rg_REAL& value);
	void setLowerValue(const rg_REAL& value);
    void setLowerAndUpperValue(const rg_REAL& lowerValue, const rg_REAL& upperValue);

 
    rg_FLAG isInIntervalWithOpenOpenEnd(const rg_REAL& value) const;
    rg_FLAG isInIntervalWithOpenCloseEnd(const rg_REAL& value) const;
    rg_FLAG isInIntervalWithCloseOpenEnd(const rg_REAL& value) const;
    rg_FLAG isInIntervalWithCloseCloseEnd(const rg_REAL& value) const;

	Interval_REAL& operator =(const Interval_REAL& interval);
};

} // namespace GeometryTier
} // namespace BULL2D

#endif

