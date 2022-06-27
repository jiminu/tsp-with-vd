#include "TimeStatistics.h"

TimeStatistics::TimeStatistics()
{
}



TimeStatistics::TimeStatistics(const int& numTimes)
{
    m_TotalTime.resize(numTimes, 0.0);
    m_StartClockTime.resize(numTimes, 0.0);
}



TimeStatistics::TimeStatistics(const TimeStatistics& statistics)
    : m_TotalTime(statistics.m_TotalTime), m_StartClockTime(statistics.m_StartClockTime)
{
}



TimeStatistics::~TimeStatistics()
{
}



TimeStatistics& TimeStatistics::operator =(const TimeStatistics& statistics)
{
    if (this != &statistics) 
    {
        m_TotalTime     = statistics.m_TotalTime;
        m_StartClockTime = statistics.m_StartClockTime;
    }

    return *this;
}

