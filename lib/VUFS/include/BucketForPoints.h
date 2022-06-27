#ifndef BUCKETFORPOINTS_H
#define BUCKETFORPOINTS_H

#include "rg_Point3D.h"
#include "AxisAlignedBox.h"
#include "rg_Const.h"

#include <vector>
#include <list>
using namespace std;


template <class T>
class BucketForPoints_Element
{
private:
    rg_Point3D  m_point;
    T           m_userData;

public:
    BucketForPoints_Element();
    BucketForPoints_Element(const rg_Point3D& point, const T& userData);
    BucketForPoints_Element(const BucketForPoints_Element<T>& element);
    ~BucketForPoints_Element();

    BucketForPoints_Element<T>& operator =(const BucketForPoints_Element<T>& element);

    rg_Point3D  point() const;
    T           user_data() const;

    void        set(const rg_Point3D& point, const T& userData);
};


template <class T>
class BucketForPoints
{
//public:
//    class Element<T>;
//
private:
    vector< list< BucketForPoints_Element<T> > >    m_data;
    AxisAlignedBox                  m_range;
    unsigned int                    m_size;

public:
    BucketForPoints();
    BucketForPoints(const BucketForPoints<T>& bucket);
    ~BucketForPoints();

    BucketForPoints& operator =(const BucketForPoints<T>& bucket);

    void    clear();

    void    initialize(const unsigned int& size, const AxisAlignedBox& range);
    void    initialize(const list< pair<rg_Point3D, T> >& points);

    void    insert(const rg_Point3D& point, const T& userData);       

    bool    has_same_data(const rg_Point3D& point, T& userData, const rg_REAL& res=rg_MATH_RES) const;
    unsigned int    find_user_data(const rg_Point3D& point, list<T>& userData, const rg_REAL& res=rg_MATH_RES) const;

private:
    unsigned int    find_bucket_index(const rg_Point3D& point) const;
};





template <class T>
BucketForPoints<T>::BucketForPoints()
: m_size(0)
{
}



template <class T>
BucketForPoints<T>::BucketForPoints(const BucketForPoints<T>& bucket)
: m_range( bucket.m_range ),
  m_size( bucket.m_size )
{
    m_data.resize( m_size );

    for ( int i=0; i<m_data.size(); ++i ) {
        m_data[i] = bucket.m_data[i];
    }
}



template <class T>
BucketForPoints<T>::~BucketForPoints()
{
}



template <class T>
BucketForPoints<T>& BucketForPoints<T>::operator =(const BucketForPoints<T>& bucket)
{
    if ( this != &bucket ) {
        m_data.clear();

        m_range = bucket.m_range;
        m_size  = bucket.m_size;

        m_data.resize( m_size );
        for ( int i=0; i<m_data.size(); ++i ) {
            m_data[i] = bucket.m_data[i];
        }
    }

    return *this;
}


    
template <class T>
void    BucketForPoints<T>::clear()
{
    m_data.clear();
    m_range.setAxisAlignedBox(  rg_Point3D(rg_REAL_INFINITY, rg_REAL_INFINITY, rg_REAL_INFINITY), 
                                rg_Point3D(-rg_REAL_INFINITY, -rg_REAL_INFINITY, -rg_REAL_INFINITY) );
    m_size = 0;
}



template <class T>
void    BucketForPoints<T>::initialize(const unsigned int& size, const AxisAlignedBox& range)
{
    clear();

    m_range = range;
    m_size  = size;
    m_data.resize( m_size+1 );
}



template <class T>
void    BucketForPoints<T>::initialize(const list< pair<rg_Point3D, T> >& points)
{
    typename list< pair<rg_Point3D, T> >::const_iterator i_pt;
    for ( i_pt=points.begin(); i_pt!=points.end(); ++i_pt ) {
        m_range.update( i_pt->first );
    }

    m_size = points.size();
    m_data.resize( m_size+1 );

    for ( i_pt=points.begin(); i_pt!=points.end(); ++i_pt ) {
        insert( i_pt->first, i_pt->second );
    }
}


    
template <class T>
bool    BucketForPoints<T>::has_same_data(const rg_Point3D& point, T& userData, const rg_REAL& res) const
{
    bool bucketHasSamePoint = false;
    unsigned int index = find_bucket_index( point );

    typename list< BucketForPoints_Element<T> >::const_iterator i_data;
    for ( i_data=m_data.at(index).begin(); i_data!=m_data.at(index).end(); ++i_data ) {
        if ( point.isEqual( i_data->point(), res ) ) {
            bucketHasSamePoint = true;
            userData = i_data->user_data();
            break;
        }
    }

    return bucketHasSamePoint;
}



template <class T>
unsigned int    BucketForPoints<T>::find_user_data(const rg_Point3D& point, list<T>& userData, const rg_REAL& res) const
{
    unsigned int numSamePoint = 0;
    unsigned int index = find_bucket_index( point );

    typename list< BucketForPoints_Element<T> >::const_iterator i_data;
    for ( i_data=m_data.at(index).begin(); i_data!=m_data.at(index).end(); ++i_data ) {
        if ( point.isEqual( i_data->point(), res ) ) {
            userData.push_back( i_data->user_data() );
            ++numSamePoint;
        }
    }

    return numSamePoint;
}



template <class T>
void    BucketForPoints<T>::insert(const rg_Point3D& point, const T& userData)    
{
    unsigned int index = find_bucket_index( point );
    m_data.at(index).push_back( BucketForPoints_Element<T>(point, userData) );
}



template <class T>
unsigned int    BucketForPoints<T>::find_bucket_index(const rg_Point3D& point) const
{
    rg_Point3D minPt = m_range.getMinPoint();
    rg_Point3D maxPt = m_range.getMaxPoint();
    rg_Point3D differenceBtwMinAndMax = m_range.getSize();

    double sumOfCoordRatio   =   0.0;
    sumOfCoordRatio +=  (point.getX() - minPt.getX())/differenceBtwMinAndMax.getX();
    sumOfCoordRatio +=  (point.getY() - minPt.getY())/differenceBtwMinAndMax.getY();
    sumOfCoordRatio +=  (point.getZ() - minPt.getZ())/differenceBtwMinAndMax.getZ();
    sumOfCoordRatio /=  3.0f;

    sumOfCoordRatio *= m_size;

    unsigned int index = (unsigned int)sumOfCoordRatio;
    if ( index<0 || index>m_size ) {
        int stop = 1;
    }

    return index;

}









template <class T>
BucketForPoints_Element<T>::BucketForPoints_Element()
{
}



template <class T>
BucketForPoints_Element<T>::BucketForPoints_Element(const rg_Point3D& point, const T& userData)
: m_point( point ),
  m_userData( userData )
{
}



template <class T>
BucketForPoints_Element<T>::BucketForPoints_Element(const BucketForPoints_Element<T>& element)
: m_point( element.m_point ),
  m_userData( element.m_userData )
{
}



template <class T>
BucketForPoints_Element<T>::~BucketForPoints_Element()
{
}



template <class T>
BucketForPoints_Element<T>& BucketForPoints_Element<T>::operator =(const BucketForPoints_Element<T>& element)
{
    if ( this != &element ) {
        m_point    = element.m_point;
        m_userData = element.m_userData;
    }

    return *this;
}



template <class T>
rg_Point3D  BucketForPoints_Element<T>::point() const
{
    return m_point;
}



template <class T>
T           BucketForPoints_Element<T>::user_data() const
{
    return m_userData;
}



template <class T>
void        BucketForPoints_Element<T>::set(const rg_Point3D& point, const T& userData)
{
    m_point    = point;
    m_userData = userData;
}




#endif



