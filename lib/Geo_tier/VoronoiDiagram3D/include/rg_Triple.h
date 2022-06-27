#ifndef _TRIPLE_H
#define _TRIPLE_H

namespace V {

namespace GeometryTier {


template<class T1, class T2, class T3>
class rg_Triple
{
public:
    T1 m_first;
    T2 m_second;
	T3 m_third;

    rg_Triple();
    rg_Triple(const T1& first, const T2& second, const T3& third);
    rg_Triple(const rg_Triple& aTriple);
    ~rg_Triple();

    void setTriple(const T1& first, const T2& second, const T3& third);

    rg_Triple& operator =(const rg_Triple& aTriple);
};




template<class T1, class T2, class T3>
rg_Triple<T1, T2, T3>::rg_Triple()
: m_first( T1() ), m_second( T2() ), m_third( T3() )
{}



template<class T1, class T2, class T3>
rg_Triple<T1, T2, T3>::rg_Triple(const T1& first, const T2& second, const T3& third)
: m_first( first ), m_second( second ), m_third( third )
{}



template<class T1, class T2, class T3>
rg_Triple<T1, T2, T3>::rg_Triple(const rg_Triple<T1, T2, T3>& aTriple)
: m_first( aTriple.m_first ), m_second( aTriple.m_second ), m_third( aTriple.m_third )
{}



template<class T1, class T2, class T3>
rg_Triple<T1, T2, T3>::~rg_Triple()
{}



template<class T1, class T2, class T3>
void rg_Triple<T1, T2, T3>::setTriple(const T1& first, const T2& second, const T3& third)
{
    m_first  = first;
    m_second = second;
    m_third  = third;
}



template<class T1, class T2, class T3>
rg_Triple<T1, T2, T3>& rg_Triple<T1, T2, T3>::operator =(const rg_Triple<T1, T2, T3>& aTriple)
{
    if ( this == &aTriple )
        return *this;

    m_first  = aTriple.m_first;
    m_second = aTriple.m_second;
    m_third  = aTriple.m_third;

    return *this;
}


} // namespace GeometryTier

} // namespace V


#endif

