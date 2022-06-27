#ifndef _TIME_STATISTICS_H_
#define _TIME_STATISTICS_H_

#include <float.h>
#include <string>
#include <vector>
#include <chrono>
using namespace std;
using namespace std::chrono;

template <class T>
class TimeStatistics
{
private:
    vector<int>                      m_ElapsedTimes;
    vector<steady_clock::time_point> m_StartClocks;
    vector<string>                   m_ProcessName;

public:
    TimeStatistics();
    TimeStatistics(const int& numTimes);
    TimeStatistics(const TimeStatistics& statistics);
    TimeStatistics(TimeStatistics&& statistics);
    ~TimeStatistics();

    TimeStatistics& operator =(const TimeStatistics& statistics);
    TimeStatistics& operator =(TimeStatistics&& statistics);

    inline int     size() const;
    inline void    reset(const int& numTimes = 0);
    inline int     time(const unsigned int& it) const;
    inline double  timeInSec(const unsigned int& it) const;
    inline void    setTime(const unsigned int& it, const int& time);
    inline void    setTime(const int& time);

    inline void    start_clock(const unsigned int& it);
    inline void    end_clock(const unsigned int& it);

    inline void    start_clock();
    inline void    end_clock();

    inline string  process_name(const unsigned int& it) const;
    inline void    set_process_name(const unsigned int& it, const string& processName);
    inline void    set_process_name(const string& processName);
    
private:
    void copy_from(const TimeStatistics& statistics);
    void move_from(TimeStatistics&& statistics);
};


template <class T>
TimeStatistics<T>::TimeStatistics()
{
}


template <class T>
TimeStatistics<T>::TimeStatistics(const int& numTimes)
{
    m_ElapsedTimes.resize(numTimes, 0);
    //m_StartClocks.resize(numTimes, steady_clock::time_point::min());
    m_StartClocks.resize(numTimes);
    m_ProcessName.resize(numTimes);
}


template <class T>
TimeStatistics<T>::TimeStatistics(const TimeStatistics& statistics)
{
    copy_from(statistics);
}


template <class T>
TimeStatistics<T>::TimeStatistics(TimeStatistics&& statistics)
{
    move_from(move(statistics));
}



template <class T>
TimeStatistics<T>::~TimeStatistics()
{
}


template <class T>
TimeStatistics<T>& TimeStatistics<T>::operator =(const TimeStatistics& statistics)
{
    if (this != &statistics)
    {
        copy_from(statistics);
    }

    return *this;
}


template <class T>
TimeStatistics<T>& TimeStatistics<T>::operator=(TimeStatistics&& statistics)
{
    if (this != &statistics)
    {
        move_from(move(statistics));
    }

    return *this;
}


template <class T>
void TimeStatistics<T>::copy_from(const TimeStatistics& statistics)
{
    m_ElapsedTimes = statistics.m_ElapsedTimes;
    m_StartClocks  = statistics.m_StartClocks;
    m_ProcessName  = statistics.m_ProcessName;
}


template <class T>
void TimeStatistics<T>::move_from(TimeStatistics&& statistics)
{
    m_ElapsedTimes = move(statistics.m_ElapsedTimes);
    m_StartClocks  = move(statistics.m_StartClocks);
    m_ProcessName  = move(statistics.m_ProcessName);
}


template <class T>
inline int TimeStatistics<T>::size() const
{
    return m_ElapsedTimes.size();
}


template <class T>
inline int TimeStatistics<T>::time(const unsigned int& it) const
{
    if (it < m_ElapsedTimes.size())
    {
        return m_ElapsedTimes[it];
    }
    else
    {
        return INT_MAX;
    }
}


template <class T>
inline double TimeStatistics<T>::timeInSec(const unsigned int& it) const
{
    if (it < m_ElapsedTimes.size())
    {
        return m_ElapsedTimes[it]/(double)T::den;
    }
    else
    {
        return DBL_MAX;
    }
}

template <class T>
inline void TimeStatistics<T>::setTime(const unsigned int& it, const int& inputTime)
{
    if (it < m_ElapsedTimes.size()) 
    {
        m_ElapsedTimes[it] = inputTime;
    }
}

template <class T>
inline void TimeStatistics<T>::setTime(const int& time)
{
    m_ElapsedTimes.push_back(time);
}


template <class T>
inline void TimeStatistics<T>::reset(const int& numTimes)
{
    m_ElapsedTimes.clear();
    m_StartClocks.clear();
    m_ProcessName.clear();

    if (numTimes > 0) 
    {
        m_ElapsedTimes.resize(numTimes, 0);
        //m_StartClocks.resize(numTimes, steady_clock::time_point::min());
        m_StartClocks.resize(numTimes);
        m_ProcessName.resize(numTimes);
    }
}


template <class T>
inline void TimeStatistics<T>::start_clock(const unsigned int& it)
{
    if (it < m_StartClocks.size())
    {
        m_StartClocks[it] = steady_clock::now();
    }
}


template <class T>
inline void TimeStatistics<T>::end_clock(const unsigned int& it)
{
    if (it < m_StartClocks.size())
    {
        steady_clock::time_point startClock = m_StartClocks[it];
        steady_clock::time_point endClock   = steady_clock::now();

        int totalTime = time(it) + (int)(duration_cast<duration<double, T>>(endClock - startClock)).count();
        
        setTime(it, totalTime);
    }
}


template <class T>
void TimeStatistics<T>::start_clock()
{
    m_StartClocks.push_back(steady_clock::now());
}


template <class T>
inline void TimeStatistics<T>::end_clock()
{
    steady_clock::time_point startClock = m_StartClocks.back();
    steady_clock::time_point endClock   = steady_clock::now();

    int totalTime = (int)(duration_cast<duration<double, T>>(endClock - startClock)).count();

    setTime(totalTime);
}


template <class T>
inline string  TimeStatistics<T>::process_name(const unsigned int& it) const
{
    if (it < m_ProcessName.size())
    {
        return m_ProcessName[it];
    }
    else
    {
        return string("");
    }
}


template <class T>
inline void TimeStatistics<T>::set_process_name(const unsigned int& it, const string& processName)
{
    if (it < m_ProcessName.size())
    {
         m_ProcessName[it] = processName;
    }
}


template <class T>
inline void TimeStatistics<T>::set_process_name(const string& processName)
{
    m_ProcessName.push_back(processName);
}

#endif

