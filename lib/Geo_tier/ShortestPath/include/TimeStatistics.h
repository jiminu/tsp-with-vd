#ifndef _TIME_STATISTICS_H_
#define _TIME_STATISTICS_H_

#include <float.h>
#include <time.h>
#include <vector>
using namespace std;


class TimeStatistics
{
private:
    vector<double> m_TotalTime;
    vector<double> m_StartClockTime;
public:
    TimeStatistics();
    TimeStatistics(const int& numTimes);
    TimeStatistics(const TimeStatistics& statistics);
    ~TimeStatistics();

    TimeStatistics& operator =(const TimeStatistics& statistics);

    int     size() const;
    void    reset(const int& numTimes = 0);
    double  time(const unsigned int& it) const;
    void    setTime(const unsigned int& it, const double& time);
    void    setTime(const double& time);

    void    start_clock(const unsigned int& it);
    void    end_clock(const unsigned int& it);

    void    start_clock();
    void    end_clock();
};

inline int TimeStatistics::size() const
{
    return m_TotalTime.size();
}


inline double TimeStatistics::time(const unsigned int& it) const 
{
    if (it < m_TotalTime.size())
    {
        return m_TotalTime[it];
    }
    else
    {
        return DBL_MAX;
    }
}

inline void TimeStatistics::setTime(const unsigned int& it, const double& inputTime)
{
    if (it < m_TotalTime.size()) 
    {
        m_TotalTime[it] = inputTime;
    }
}


inline void TimeStatistics::setTime(const double& time)
{
    m_TotalTime.push_back(time);
}


inline void TimeStatistics::reset(const int& numTimes) 
{
    m_TotalTime.clear();
    m_StartClockTime.clear();

    if (numTimes > 0) 
    {
        m_TotalTime.resize(numTimes, 0.0);
        m_StartClockTime.resize(numTimes, 0.0);
    }
}


inline void TimeStatistics::start_clock(const unsigned int& it)
{
    if (it < m_StartClockTime.size())
    {
        m_StartClockTime[it] = ((double)clock() / CLOCKS_PER_SEC);
    }
}


inline void TimeStatistics::end_clock(const unsigned int& it)
{
    if (it < m_StartClockTime.size())
    {
        double startClock = m_StartClockTime[it];
        double endClock   = ((double)clock() / CLOCKS_PER_SEC);

        double totalTime = time(it) + (endClock - startClock);
        
        setTime(it, totalTime);
    }
}


inline void TimeStatistics::start_clock()
{
    m_StartClockTime.push_back(((double)clock() / CLOCKS_PER_SEC));
}



inline void TimeStatistics::end_clock()
{
    double startClock = m_StartClockTime.back();
    double endClock = ((double)clock() / CLOCKS_PER_SEC);

    double totalTime = endClock - startClock;

    setTime(totalTime);
}




#endif

