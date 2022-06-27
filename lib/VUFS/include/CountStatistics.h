#ifndef _COUNT_STATISTICS_H_
#define _COUNT_STATISTICS_H_

#include <float.h>
#include <vector>
#include <climits>
using namespace std;


class CountStatistics
{
private:    
    vector<int>    m_Counts;

public:
    CountStatistics();
    CountStatistics(const CountStatistics& statistics);
    ~CountStatistics();

    CountStatistics& operator =(const CountStatistics& statistics);

    int     size() const;
    void    reset(const int& numCounts = 0);
    int     count(const unsigned int& it) const;
    void    setCount(const unsigned int& it, const int& counts);
    void    setCount(const int& counts);

    void    addCount(const unsigned int& it);
    void    addCount();
};


inline int CountStatistics::size() const
{
    return m_Counts.size();
}



inline int  CountStatistics::count(const unsigned int& it) const 
{
    if (it < m_Counts.size())
    {
        return m_Counts[it];
    }
    else
    {
        return INT_MAX;
    }
}


inline void CountStatistics::setCount(const int& counts)
{
    m_Counts.push_back(counts);
}


inline void    CountStatistics::setCount(const unsigned int& it, const int& counts)
{
    if (it < m_Counts.size()) 
    {
        m_Counts[it] = counts;
    }
}


inline void    CountStatistics::reset(const int& numCounts) 
{
    m_Counts.clear();

    if (numCounts > 0) 
    {
        m_Counts.resize(numCounts, 0);
    }
}


inline void CountStatistics::addCount()
{
    if (!m_Counts.empty())
    {
		m_Counts.back() += 1;
    }
}


inline void CountStatistics::addCount(const unsigned int& it)
{
    if (it < m_Counts.size())
    {
        m_Counts[it] = m_Counts[it] + 1;
    }
}



#endif

