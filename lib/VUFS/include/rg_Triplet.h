#ifndef RG_TRIPLET_H
#define RG_TRIPLET_H


template<class T1, class T2, class T3>
class rg_Triplet
{
public:
    T1 m_first;
    T2 m_second;
    T3 m_third;

    rg_Triplet();
    rg_Triplet(const T1& first, const T2& second, const T3& third);
    rg_Triplet(const rg_Triplet& triplet);
    ~rg_Triplet();

    rg_Triplet& operator =(const rg_Triplet& triplet);
};

template<class T1, class T2, class T3>
rg_Triplet<T1, T2, T3>::rg_Triplet()
: m_first( T1() ), m_second( T2() ), m_third( T3() )
{
}

template<class T1, class T2, class T3>
rg_Triplet<T1, T2, T3>::rg_Triplet(const T1& first, const T2& second, const T3& third)
: m_first( first ), m_second( second ), m_third( third )
{
}

template<class T1, class T2, class T3>
rg_Triplet<T1, T2, T3>::rg_Triplet(const rg_Triplet<T1, T2, T3>& triplet)
: m_first( triplet.m_first ), m_second( triplet.m_second ), m_third( triplet.m_third )
{
}

template<class T1, class T2, class T3>
rg_Triplet<T1, T2, T3>::~rg_Triplet()
{
}

template<class T1, class T2, class T3>
rg_Triplet<T1, T2, T3>& rg_Triplet<T1, T2, T3>::operator =(const rg_Triplet<T1, T2, T3>& triplet)
{
    if ( this == &triplet )
        return *this;

    m_first  = triplet.m_first;
    m_second = triplet.m_second;
    m_third  = triplet.m_third;

    return *this;
}


#endif
