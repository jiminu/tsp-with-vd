#ifndef _PAIR_H
#define _PAIR_H

namespace V {

namespace GeometryTier {


template<class T1, class T2>
class Pair
{
public:
    T1 m_first;
    T2 m_second;

    Pair();
    Pair(const T1& first, const T2& second);
    Pair(const Pair& aPair);
    ~Pair();

    Pair& operator =(const Pair& aPair);
};

template<class T1, class T2>
Pair<T1, T2>::Pair()
: m_first( T1() ), m_second( T2() )
{
}

template<class T1, class T2>
Pair<T1, T2>::Pair(const T1& first, const T2& second)
: m_first( first ), m_second( second )
{
}

template<class T1, class T2>
Pair<T1, T2>::Pair(const Pair<T1, T2>& aPair)
: m_first( aPair.m_first ), m_second( aPair.m_second )
{
}

template<class T1, class T2>
Pair<T1, T2>::~Pair()
{
}

template<class T1, class T2>
Pair<T1, T2>& Pair<T1, T2>::operator =(const Pair<T1, T2>& aPair)
{
    if ( this == &aPair )
        return *this;

    m_first  = aPair.m_first;
    m_second = aPair.m_second;

    return *this;
}


} // namespace GeometryTier

} // namespace V


#endif

