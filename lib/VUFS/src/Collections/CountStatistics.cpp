#include "CountStatistics.h"

CountStatistics::CountStatistics()
{
}



CountStatistics::CountStatistics(const CountStatistics& statistics)
    : m_Counts(statistics.m_Counts)
{
}



CountStatistics::~CountStatistics()
{
}



CountStatistics& CountStatistics::operator =(const CountStatistics& statistics)
{
    if (this != &statistics) 
    {
        m_Counts                       = statistics.m_Counts;
    }

    return *this;
}

